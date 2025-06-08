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
    fn test_caulk_by_input() {
        let c = scalars(&[1, 3, 2, 4, 5, 8, 7, 6]); // Table size 8
        let values = scalars(&[1, 2, 3, 4]); // Lookup size 4
        let timings = TestCaulk::test_caulk_by_input(&c, &values);
        println!("Caulk by Input Timings:");
        for timing in timings {
            println!("  {}", timing);
        }
    }

    #[test]
    fn test_caulk_with_params() {
        let c = scalars(&[1, 3, 2, 4]); // Table N=4
        let values = scalars(&[1, 2]); // Lookups m=2
        let N = c.len();
        let m = values.len();

        // Pre-generate params
        let (pp, vp) = TestCaulk::setup(N, m).expect("Setup should not fail");

        let timings = TestCaulk::test_caulk_with_params(&pp, &vp, &c, &values);
        println!("Caulk with Params Timings (Prove/Verify only):");
        for timing in timings {
            println!("  {}", timing);
        }
    }

    #[test]
    fn test_caulk_by_k() {
        let timings = TestCaulk::test_caulk_by_k(4);
        println!("Caulk by K Timings:");
        for timing in timings {
            println!("  {}", timing);
        }
    }

    #[test]
    fn test_caulk_with_positions() {
        let c = scalars(&[1, 3, 2, 4]); // Table N=4: positions [0, 1, 2, 3]
        let positions = vec![0, 2]; // Direct positions for values [1, 2]
        let N = c.len();
        let m = positions.len();

        // Pre-generate params
        let (pp, vp) = TestCaulk::setup(N, m).expect("Setup should not fail");

        let timings = TestCaulk::test_caulk_with_positions(&pp, &vp, &c, &positions);
        println!("Caulk with Direct Positions Timings (Optimized):");
        for timing in timings {
            println!("  {}", timing);
        }
    }
} 