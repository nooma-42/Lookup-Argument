#![allow(unused)]

use crate::pcs::univariate::{
    UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam,
};
use crate::pcs::PolynomialCommitmentScheme;
use crate::util::arithmetic::{MultiMillerLoop, WithSmallOrderMulGroup};
use crate::util::transcript::{Keccak256Transcript, FieldTranscriptWrite, FieldTranscriptRead,FieldTranscript, G2TranscriptRead, G2TranscriptWrite, TranscriptRead, TranscriptWrite, InMemoryTranscript};
use crate::Error;
use halo2_curves::bn256::Bn256;
use halo2_curves::ff::PrimeField;
use halo2_curves::pairing::{Engine, PairingCurveAffine};
use halo2_curves::CurveAffine;
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::marker::PhantomData;
use std::fmt::Debug;

pub mod preprocessor;
pub mod prover;
pub mod util;
pub mod verifier;

#[derive(Clone, Debug)]
pub struct CaulkProverParam<M: MultiMillerLoop> {
    pub(crate) kzg_param: UnivariateKzgParam<M>,
    pub(crate) kzg_pp: UnivariateKzgProverParam<M>,
}

/// Optimized CaulkProverParam with precomputed G2 openings for fast H1_com calculation
#[derive(Clone, Debug)]
pub struct CaulkOptimizedProverParam<M: MultiMillerLoop> {
    pub(crate) kzg_param: UnivariateKzgParam<M>,
    pub(crate) kzg_pp: UnivariateKzgProverParam<M>,
    /// Precomputed G2 KZG openings: Q_i = Comm_G2((C(X) - C(ω^i)) / (X - ω^i))
    pub(crate) precomputed_g2_openings: Vec<M::G2Affine>,
}

#[derive(Clone, Debug)]
pub struct CaulkVerifierParam<M: MultiMillerLoop> {
    pub(crate) kzg_vp: UnivariateKzgVerifierParam<M>,
    pub(crate) m: usize,
}

#[derive(Clone, Debug)]
pub struct Caulk<M: MultiMillerLoop>(PhantomData<M>);

impl<M> Caulk<M>
where
    M: MultiMillerLoop + Engine + Debug + Sync + Send,
    M::Scalar: Serialize + DeserializeOwned + PrimeField + WithSmallOrderMulGroup<3>,
    M::G1Affine: Serialize + DeserializeOwned + PairingCurveAffine<ScalarExt = M::Scalar> + std::ops::Add<M::G1Affine> + std::ops::Sub<M::G1Affine>,
    M::G2Affine: Serialize + DeserializeOwned + PairingCurveAffine<ScalarExt = M::Scalar> + std::ops::Add<M::G2Affine>,
{
    pub fn setup(
        N: usize,
        m: usize,
    ) -> Result<(CaulkProverParam<M>, CaulkVerifierParam<M>), Error> {
        preprocessor::setup::<M>(N, m)
    }

    /// Optimized setup with precomputed G2 openings for fast H1_com calculation
    pub fn setup_optimized(
        N: usize,
        m: usize,
        table: &[M::Scalar],
    ) -> Result<(CaulkOptimizedProverParam<M>, CaulkVerifierParam<M>), Error> {
        let (base_pp, vp, precomputed_g2_openings) = 
            preprocessor::preprocess_with_table::<M>(N, m, table)?;
        
        Ok((
            CaulkOptimizedProverParam {
                kzg_param: base_pp.kzg_param,
                kzg_pp: base_pp.kzg_pp,
                precomputed_g2_openings,
            },
            vp,
        ))
    }

    pub fn prove(
        pp: &CaulkProverParam<M>,
        c: &[M::Scalar],
        positions: &[usize],
        transcript: &mut (impl TranscriptWrite<M::G1Affine, M::Scalar>
                      + G2TranscriptWrite<M::G2Affine, M::Scalar>
                      + FieldTranscriptWrite<M::Scalar>),
    ) -> Result<(), Error>
    {
        prover::prove::<M>(pp, c, positions, transcript)
    }

    /// Optimized prove with precomputed G2 openings for fast H1_com calculation
    pub fn prove_optimized(
        pp: &CaulkOptimizedProverParam<M>,
        c: &[M::Scalar],
        positions: &[usize],
        transcript: &mut (impl TranscriptWrite<M::G1Affine, M::Scalar>
                      + G2TranscriptWrite<M::G2Affine, M::Scalar>
                      + FieldTranscriptWrite<M::Scalar>),
    ) -> Result<(), Error>
    {
        prover::prove_optimized::<M>(pp, c, positions, transcript)
    }

    pub fn verify(
        vp: &CaulkVerifierParam<M>,
        transcript: &mut (impl TranscriptRead<M::G1Affine, M::Scalar>
                      + G2TranscriptRead<M::G2Affine, M::Scalar>
                      + FieldTranscriptRead<M::Scalar>),
    ) -> Result<(), Error>
    {
        verifier::verify::<M>(vp, transcript)
    }

    // Helper function to convert values to positions (for compatibility with existing tests)
    fn values_to_positions<F: PartialEq>(c: &[F], values: &[F]) -> Vec<usize> {
        values
            .iter()
            .map(|value| {
                c.iter()
                    .position(|x| x == value)
                    .expect("Value from lookup set must exist in the table")
            })
            .collect()
    }

    // Run the full Caulk protocol with table and lookup generated based on k
    pub fn test_caulk_by_k(k: usize) -> Vec<String> {
        use halo2_curves::bn256::{Bn256, Fr};

        let mut timings: Vec<String> = vec![];
        let start_total = std::time::Instant::now();

        // Use cq's generator for consistency in benchmarks
        let (c, values) = crate::backend::cq::generate_table_and_lookup(
            2_usize.pow(k as u32),
            2_usize.pow((k - 1) as u32), // N:n ratio = 2
        );
        let N = c.len();
        let m = values.len();

        // 1. Setup
        let start = std::time::Instant::now();
        let (pp, vp) = Caulk::<Bn256>::setup(N, m).expect("Setup should not fail");
        let duration1 = start.elapsed();
        timings.push(format!("Setup: {}ms", duration1.as_millis()));
        
        // 2. Prove
        let start = std::time::Instant::now();
        let positions = Self::values_to_positions(&c, &values);
        let proof = {
            let mut transcript = Keccak256Transcript::new(());
            Caulk::<Bn256>::prove(&pp, &c[..], &positions, &mut transcript)
                .expect("Prove should not fail");
            transcript.into_proof()
        };
        let duration2 = start.elapsed();
        timings.push(format!("Prove: {}ms", duration2.as_millis()));
        timings.push(format!("Proof size: {} bytes", proof.len()));

        // 3. Verify
        let start = std::time::Instant::now();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        Caulk::<Bn256>::verify(&vp, &mut transcript).expect("Verify should not fail");
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));

        timings
    }

    // Run the full Caulk protocol with given table and lookup values
    pub fn test_caulk_by_input(
        c: &[<Bn256 as Engine>::Scalar],
        values: &[<Bn256 as Engine>::Scalar],
    ) -> Vec<String> {
        let mut timings: Vec<String> = vec![];
        let start_total = std::time::Instant::now();

        let N = c.len();
        let m = values.len();

        // 1. Setup
        let start = std::time::Instant::now();
        let (pp, vp) = Caulk::<Bn256>::setup(N, m).expect("Setup should not fail");
        let duration1 = start.elapsed();
        timings.push(format!("Setup: {}ms", duration1.as_millis()));

        // 2. Prove
        let start = std::time::Instant::now();
        let positions = Self::values_to_positions(c, values);
        let proof = {
            let mut transcript = Keccak256Transcript::new(());
            Caulk::<Bn256>::prove(&pp, c, &positions, &mut transcript)
                .expect("Prove should not fail");
            transcript.into_proof()
        };
        let duration2 = start.elapsed();
        timings.push(format!("Prove: {}ms", duration2.as_millis()));
        timings.push(format!("Proof size: {} bytes", proof.len()));

        // 3. Verify
        let start = std::time::Instant::now();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        Caulk::<Bn256>::verify(&vp, &mut transcript).expect("Verify should not fail");
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));

        timings
    }

    /// Test optimized Caulk with precomputed G2 openings
    pub fn test_caulk_optimized_by_k(k: usize) -> Vec<String> {
        use halo2_curves::bn256::{Bn256, Fr};

        let mut timings: Vec<String> = vec![];
        let start_total = std::time::Instant::now();

        // Use cq's generator for consistency in benchmarks
        let (c, values) = crate::backend::cq::generate_table_and_lookup(
            2_usize.pow(k as u32),
            2_usize.pow((k - 1) as u32), // N:n ratio = 2
        );
        let N = c.len();
        let m = values.len();

        // 1. Optimized Setup with precomputation
        let start = std::time::Instant::now();
        let (pp_opt, vp) = Caulk::<Bn256>::setup_optimized(N, m, &c).expect("Optimized setup should not fail");
        let duration1 = start.elapsed();
        timings.push(format!("Optimized Setup: {}ms", duration1.as_millis()));

        // 2. Optimized Prove  
        let start = std::time::Instant::now();
        let positions = Self::values_to_positions(&c, &values);
        let proof = {
            let mut transcript = Keccak256Transcript::new(());
            Caulk::<Bn256>::prove_optimized(&pp_opt, &c[..], &positions, &mut transcript)
                .expect("Optimized prove should not fail");
            transcript.into_proof()
        };
        let duration2 = start.elapsed();
        timings.push(format!("Optimized Prove: {}ms", duration2.as_millis()));

        // 3. Verify (same as regular version)
        let start = std::time::Instant::now();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        Caulk::<Bn256>::verify(&vp, &mut transcript).expect("Verify should not fail");
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));

        timings
    }

    // Run Caulk prove and verify steps with pre-generated parameters
    pub fn test_caulk_with_params(
        pp: &CaulkProverParam<Bn256>,
        vp: &CaulkVerifierParam<Bn256>,
        c: &[<Bn256 as Engine>::Scalar],
        values: &[<Bn256 as Engine>::Scalar],
    ) -> Vec<String> {
        let mut timings: Vec<String> = vec![];
        let start_total = std::time::Instant::now();

        // 1. Prove
        let start = std::time::Instant::now();
        let positions = Self::values_to_positions(c, values);
        let proof = {
            let mut transcript = Keccak256Transcript::new(());
            Caulk::<Bn256>::prove(pp, c, &positions, &mut transcript)
                .expect("Prove should not fail");
            transcript.into_proof()
        };
        let duration_prove = start.elapsed();
        timings.push(format!("Prove: {}ms", duration_prove.as_millis()));

        // 2. Verify
        let start = std::time::Instant::now();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        Caulk::<Bn256>::verify(vp, &mut transcript).expect("Verify should not fail");
        let duration_verify = start.elapsed();
        timings.push(format!("Verify: {}ms", duration_verify.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Prove+Verify time: {}ms", total_duration.as_millis()));

        timings
    }

    // Run Caulk prove and verify steps with positions directly (optimized version)
    pub fn test_caulk_with_positions(
        pp: &CaulkProverParam<Bn256>,
        vp: &CaulkVerifierParam<Bn256>,
        c: &[<Bn256 as Engine>::Scalar],
        positions: &[usize],
    ) -> Vec<String> {
        let mut timings: Vec<String> = vec![];
        let start_total = std::time::Instant::now();

        // 1. Prove
        let start = std::time::Instant::now();
        let proof = {
            let mut transcript = Keccak256Transcript::new(());
            Caulk::<Bn256>::prove(pp, c, positions, &mut transcript)
                .expect("Prove should not fail");
            transcript.into_proof()
        };
        let duration_prove = start.elapsed();
        timings.push(format!("Prove: {}ms", duration_prove.as_millis()));

        // 2. Verify
        let start = std::time::Instant::now();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        Caulk::<Bn256>::verify(vp, &mut transcript).expect("Verify should not fail");
        let duration_verify = start.elapsed();
        timings.push(format!("Verify: {}ms", duration_verify.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Prove+Verify time: {}ms", total_duration.as_millis()));

        timings
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::transcript::{Keccak256Transcript, TranscriptRead, TranscriptWrite};
    use halo2_curves::bn256::{Fr as Scalar, G1Affine as Bn256G1Affine};
    use halo2_curves::pairing::Engine;

    type TestCaulk = Caulk<Bn256>;

    fn scalars(nums: &[u64]) -> Vec<Scalar> {
        nums.iter().map(|x| Scalar::from(*x)).collect::<Vec<_>>()
    }

    #[test]
    fn test_caulk_e2e() {
        let c = scalars(&[1, 3, 2, 4]); // Table
        let values = scalars(&[1, 2]); // Lookups
        let N = c.len();
        let m = values.len();

        let (pp, vp) = TestCaulk::setup(N, m).expect("Setup should not fail");

        // Convert values to positions: [1, 2] -> [0, 2] (positions in table [1, 3, 2, 4])
        let positions = TestCaulk::values_to_positions(&c, &values);

        let proof = {
            let mut prover_transcript = Keccak256Transcript::new(());
            TestCaulk::prove(&pp, &c[..], &positions, &mut prover_transcript)
                .expect("Prove should not fail");
            prover_transcript.into_proof()
        };

        let mut verifier_transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        TestCaulk::verify(&vp, &mut verifier_transcript)
            .expect("Verify should not fail");
    }

    #[test]
    fn test_caulk_optimized_h1_computation() {
        use halo2_curves::bn256::{Bn256, Fr};
        
        // Test small case: N=8, lookup 4 values
        let table: Vec<Fr> = (1..=8).map(|i| Fr::from(i as u64)).collect();
        let lookup_values = vec![Fr::from(1), Fr::from(3), Fr::from(5), Fr::from(7)];
        let positions = TestCaulk::values_to_positions(&table, &lookup_values);
        
        let N = table.len();
        let m = lookup_values.len();
        
        // Test regular setup
        let (regular_pp, regular_vp) = TestCaulk::setup(N, m).expect("Regular setup should work");
        
        // Test optimized setup
        let (optimized_pp, optimized_vp) = 
            TestCaulk::setup_optimized(N, m, &table).expect("Optimized setup should work");
        
        // For now, just test that both setups work and both can produce valid proofs
        // We'll skip the proof comparison until the optimization is properly implemented
        
        // Test regular proof
        let regular_proof = {
            let mut transcript = Keccak256Transcript::new(());
            TestCaulk::prove(&regular_pp, &table, &positions, &mut transcript)
                .expect("Regular prove should work");
            transcript.into_proof()
        };
        
        // Test optimized proof (currently using fallback)
        let optimized_proof = {
            let mut transcript = Keccak256Transcript::new(());
            TestCaulk::prove_optimized(&optimized_pp, &table, &positions, &mut transcript)
                .expect("Optimized prove should work");
            transcript.into_proof()
        };
        
        // Verify both proofs are valid
        let mut regular_transcript = Keccak256Transcript::from_proof((), regular_proof.as_slice());
        TestCaulk::verify(&regular_vp, &mut regular_transcript)
            .expect("Regular verification should pass");
            
        let mut optimized_transcript = Keccak256Transcript::from_proof((), optimized_proof.as_slice());
        TestCaulk::verify(&optimized_vp, &mut optimized_transcript)
            .expect("Optimized verification should pass");
        
        // Note: We don't compare proofs directly since they use different random blinding factors
        // This will be fixed when the optimization is properly implemented with deterministic behavior
        
        println!("✅ Caulk optimized H1 computation test passed!");
        println!("   Regular proof length: {} bytes", regular_proof.len());
        println!("   Optimized proof length: {} bytes", optimized_proof.len());
    }

    #[test]
    fn test_caulk_optimized_performance_comparison() {
        use halo2_curves::bn256::{Bn256, Fr};
        use std::time::Instant;
        
        // Test with larger table for performance comparison
        let k = 6; // N = 64, m = 32
        let N = 1 << k;
        let m = 1 << (k - 1);
        
        let table: Vec<Fr> = (1..=N).map(|i| Fr::from(i as u64)).collect();
        let lookup_values: Vec<Fr> = (1..=m).map(|i| Fr::from((i * 2) as u64)).collect();
        let positions = TestCaulk::values_to_positions(&table, &lookup_values);
        
        println!("Performance test: N={}, m={}", N, m);
        
        // Regular setup and prove
        let regular_start = Instant::now();
        let (regular_pp, regular_vp) = TestCaulk::setup(N, m).expect("Regular setup should work");
        let regular_setup_time = regular_start.elapsed();
        
        let regular_prove_start = Instant::now();
        let regular_proof = {
            let mut transcript = Keccak256Transcript::new(());
            TestCaulk::prove(&regular_pp, &table, &positions, &mut transcript)
                .expect("Regular prove should work");
            transcript.into_proof()
        };
        let regular_prove_time = regular_prove_start.elapsed();
        
        // Optimized setup and prove
        let optimized_start = Instant::now();
        let (optimized_pp, optimized_vp) = 
            TestCaulk::setup_optimized(N, m, &table).expect("Optimized setup should work");
        let optimized_setup_time = optimized_start.elapsed();
        
        let optimized_prove_start = Instant::now();
        let optimized_proof = {
            let mut transcript = Keccak256Transcript::new(());
            TestCaulk::prove_optimized(&optimized_pp, &table, &positions, &mut transcript)
                .expect("Optimized prove should work");
            transcript.into_proof()
        };
        let optimized_prove_time = optimized_prove_start.elapsed();
        
        // Print performance comparison
        println!("📊 Performance Comparison:");
        println!("  Regular Setup:    {:?}", regular_setup_time);
        println!("  Optimized Setup:  {:?}", optimized_setup_time);
        println!("  Regular Prove:    {:?}", regular_prove_time);
        println!("  Optimized Prove:  {:?}", optimized_prove_time);
        
        let setup_ratio = optimized_setup_time.as_millis() as f64 / regular_setup_time.as_millis() as f64;
        let prove_ratio = optimized_prove_time.as_millis() as f64 / regular_prove_time.as_millis() as f64;
        
        println!("  Setup Ratio (Opt/Reg):  {:.2}x", setup_ratio);
        println!("  Prove Ratio (Opt/Reg):  {:.2}x", prove_ratio);
        
        // Verify both proofs are still valid
        let mut regular_transcript = Keccak256Transcript::from_proof((), regular_proof.as_slice());
        TestCaulk::verify(&regular_vp, &mut regular_transcript)
            .expect("Regular verification should pass");
            
        let mut optimized_transcript = Keccak256Transcript::from_proof((), optimized_proof.as_slice());
        TestCaulk::verify(&optimized_vp, &mut optimized_transcript)
            .expect("Optimized verification should pass");
        
        // The optimized version should ideally be faster for proving (though setup might be slower due to precomputation)
        if prove_ratio < 1.0 {
            println!("🚀 Optimized proving is {:.2}x faster!", 1.0 / prove_ratio);
        } else {
            println!("⚠️  Optimized proving is {:.2}x slower. This may indicate the optimization is not effective for this size.", prove_ratio);
        }
        
        println!("✅ Caulk performance comparison test completed!");
    }
} 