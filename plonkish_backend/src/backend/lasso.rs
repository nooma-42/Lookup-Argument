// plonkish_backend/src/backend/lasso.rs
use crate::pcs::multilinear::kzg::MultilinearKzg;
use crate::{
    pcs::{self, Evaluation, Point, PolynomialCommitmentScheme}, // Import PolynomialCommitmentScheme and Point
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::{Field, PrimeField}, // Removed: powers (unused)
        // expression::Query, // Removed: unused
        transcript::{
            FieldTranscriptRead, FieldTranscriptWrite, InMemoryTranscript, Keccak256Transcript,
        },
        // timer::Timer, // Removed: unused
    },
    Error,
};
use halo2_curves::ff::WithSmallOrderMulGroup;
use std::{fmt::Debug, hash::Hash, marker::PhantomData}; // Removed HashMap, Instant

pub mod preprocessor;
pub mod prover;
pub mod sostable;
pub mod util;
pub mod verifier;

// --- Structures moved to bottom for clarity, using definitions consistent with prover/verifier ---

#[derive(Clone, Debug)]
pub struct Lasso<F, Pcs>(PhantomData<(F, Pcs)>);

// --- Core Trait Implementation ---

impl<F, Pcs> Lasso<F, Pcs>
where
    F: PrimeField + WithSmallOrderMulGroup<3> + Hash + Eq + Send + Sync,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    // Bounds match PolynomialCommitmentScheme associated types
    Pcs::Commitment:
        Default + Clone + Debug + PartialEq + Send + Sync + AsRef<[Pcs::CommitmentChunk]>,
    Pcs::CommitmentChunk: Clone + Debug + Default + Send + Sync,
    // Pcs::Proof removed
    Pcs::ProverParam: Clone + Debug + Send + Sync,
    Pcs::VerifierParam: Clone + Debug + Send + Sync,
    Pcs::Param: Clone + Debug + Send + Sync,
{
    /// Preprocess the lookup information.
    pub fn preprocess(
        pcs_param_raw: &Pcs::Param,
        info: &LassoInfo<F, Pcs>,
    ) -> Result<
        (
            prover::LassoProverParam<F, Pcs>,
            verifier::LassoVerifierParam<F, Pcs>,
        ),
        Error,
    > {
        // Use verifier::LassoVerifierParam
        preprocessor::preprocess(pcs_param_raw, info)
    }

    /// Create a Lasso proof.
    pub fn prove(
        pp: prover::LassoProverParam<F, Pcs>,
        witness: Vec<F>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<LassoProof<F, Pcs>, Error> {
        prover::prove::<F, Pcs>(pp, witness, transcript)
    }

    /// Verify a Lasso proof.
    pub fn verify(
        vp: verifier::LassoVerifierParam<F, Pcs>, // Use verifier::LassoVerifierParam
        proof: &LassoProof<F, Pcs>,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(), Error> {
        verifier::verify::<F, Pcs>(vp, proof, transcript)
    }
}

// --- Benchmark Integration ---
// ... existing benchmark/test setup ...
// type BenchPcs = Zeromorph<UnivariateKzg<Bn256>>;
// ... existing benchmark/test setup ...

// --- Tests ---
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        pcs::multilinear::MultilinearKzg,
        // util::arithmetic::powers, // Removed: unused
        util::test::std_rng, // Use std_rng from test utils
        util::transcript::{InMemoryTranscript, Keccak256Transcript},
    };
    use halo2_curves::bn256::{Bn256, Fr}; // Removed G1Affine (unused)
    use rand::rngs::OsRng; // Use OsRng for setup, Rng for witness gen
                           // use std::iter; // Removed: unused
    use std::time::Instant; // Keep Instant for timing

    // Define type alias for KZG
    type Kzg = MultilinearKzg<Bn256>;
    type DefaultLasso = Lasso<Fr, Kzg>;

    #[test]
    fn test_lasso_simple_lookup() {
        let mut rng = std_rng(); // Use std_rng for witness generation
        let start_total = Instant::now();

        // Parameters
        let l = 3; // Subtable size 2^3 = 8
        let c = 1; // 1 chunk
        let k = 1; // 1 subtable type per chunk
                   // let alpha = c * k; // Derived inside preprocess
        let subtable_size = 1 << l;
        let m: usize = 8; // Witness size - specify type
        let logm = m.next_power_of_two().ilog2() as usize;

        // Simple subtable T = [0, 1, ..., 7]
        let subtable: Vec<Fr> = (0..subtable_size).map(|i| Fr::from(i as u64)).collect();
        let subtables = vec![subtable.clone()];

        // Simple witness a = [3, 1, 4, 1, 5, 0, 2, 7] (m=8)
        let witness: Vec<Fr> = [3u64, 1, 4, 1, 5, 0, 2, 7] // Use values within [0, 7]
            .iter()
            .map(|&x| Fr::from(x))
            .collect();
        assert_eq!(witness.len(), m);

        // Setup PCS
        let max_vars = logm.max(l).max(l + 1); // Max degree needed for all polynomials (incl. GP)
        let kzg_degree = 1 << max_vars;
        let kzg_params = Kzg::setup(kzg_degree, 0, &mut OsRng).expect("Kzg setup failed"); // Use OsRng for setup

        // Create LassoInfo
        let info = LassoInfo {
            l,
            c,
            k,
            subtables,
            witness: witness.clone(),
            pcs_param_raw: kzg_params.clone(), // Add raw params
        };

        // Preprocessing
        println!("Starting preprocessing...");
        let start = Instant::now();
        let (pp, vp) = DefaultLasso::preprocess(&kzg_params, &info).expect("Preprocessing failed");
        println!("Preprocessing: {}ms", start.elapsed().as_millis());

        // Proving
        println!("Starting proving...");
        let start = Instant::now();
        let mut transcript_p = Keccak256Transcript::<Fr>::new(b"lasso_test");
        let proof =
            DefaultLasso::prove(pp, witness.clone(), &mut transcript_p).expect("Proving failed");
        let proof_bytes = transcript_p.into_proof();
        println!("Prove: {}ms", start.elapsed().as_millis());
        println!("Proof size: {} bytes", proof_bytes.len());

        // Verification
        println!("Starting verification...");
        let start = Instant::now();
        let mut transcript_v = Keccak256Transcript::<Fr>::new(b"lasso_test");
        let result = DefaultLasso::verify(vp, &proof, &mut transcript_v);
        println!("Verify: {}ms", start.elapsed().as_millis());

        assert!(result.is_ok(), "Verification failed: {:?}", result.err());
        println!("Verification successful!");
        println!("Total time: {}ms", start_total.elapsed().as_millis());
    }

    #[test]
    fn test_lasso_small_m_pad() {
        let mut rng = std_rng(); // Use std_rng for witness generation
                                 // Test case where witness length m is not a power of 2
        let l = 3; // Subtable size 8
        let c = 1;
        let k = 1;
        // let alpha = c * k;
        let subtable_size = 1 << l;
        let subtable: Vec<Fr> = (0..subtable_size).map(|i| Fr::from(i as u64)).collect();
        let subtables = vec![subtable.clone()];

        // Witness a = [3, 1, 4] (m=3)
        let witness: Vec<Fr> = [3u64, 1, 4].iter().map(|&x| Fr::from(x)).collect();
        let m: usize = witness.len(); // 3 - specify type
        let logm = m.next_power_of_two().ilog2() as usize; // log2(4) = 2

        let max_vars = logm.max(l).max(l + 1);
        let kzg_degree = 1 << max_vars;
        let kzg_params = Kzg::setup(kzg_degree, 0, &mut OsRng).expect("Kzg setup failed (padding)"); // Use OsRng for setup

        let info = LassoInfo {
            l,
            c,
            k,
            subtables,
            witness: witness.clone(),
            pcs_param_raw: kzg_params.clone(), // Add raw params
        };

        println!("Starting preprocessing (padding)...");
        let start = Instant::now();
        let (pp, vp) =
            DefaultLasso::preprocess(&kzg_params, &info).expect("Preprocessing failed (padding)");
        println!("Preprocessing (padding): {}ms", start.elapsed().as_millis());

        println!("Starting proving (padding)...");
        let start = Instant::now();
        let mut transcript_p = Keccak256Transcript::<Fr>::new(b"lasso_pad_test");
        // Pass original witness to prove, padding is handled internally based on pp
        let proof = DefaultLasso::prove(pp, witness.clone(), &mut transcript_p)
            .expect("Proving failed (padding)");
        let proof_bytes = transcript_p.into_proof();
        println!("Prove (padding): {}ms", start.elapsed().as_millis());
        println!("Proof size (padding): {} bytes", proof_bytes.len());

        println!("Starting verification (padding)...");
        let start = Instant::now();
        let mut transcript_v = Keccak256Transcript::<Fr>::new(b"lasso_pad_test");
        let result = DefaultLasso::verify(vp, &proof, &mut transcript_v);
        println!("Verify (padding): {}ms", start.elapsed().as_millis());

        assert!(
            result.is_ok(),
            "Verification failed (padding): {:?}",
            result.err()
        );
        println!("Verification successful (padding)!");
    }

    // TODO: Add more tests (e.g., multiple chunks, multiple k, edge cases)
}

// --- Data Structures (Consistent Definitions) ---

/// Lasso public parameters info provided by the user.
#[derive(Clone, Debug)]
pub struct LassoInfo<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub l: usize, // Subtable size = 2^l
    pub c: usize, // Number of chunks
    pub k: usize, // Number of subtable types per chunk (usually k=1 for simple lookup)
    // Note: alpha = c * k is derived
    pub subtables: Vec<Vec<F>>,    // Subtables T_0, ..., T_{alpha-1}
    pub witness: Vec<F>,           // Witness assignments a = [a_0, ..., a_{m-1}]
    pub pcs_param_raw: Pcs::Param, // Raw PCS params before trimming
}

// --- Proof Structures ---

/// Lasso proof message round 1
#[derive(Clone, Debug)]
pub struct Message1<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) a_comm: Pcs::Commitment,
    pub(crate) logm: usize,
    pub(crate) dim_comm: Vec<Pcs::Commitment>,
    pub(crate) _marker: PhantomData<F>, // Keep phantom data if needed for generics not used in fields
}

/// Lasso proof message round 2
#[derive(Clone, Debug)]
pub struct Message2<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) a_eval: F,
    pub(crate) E_comm: Vec<Pcs::Commitment>,
    pub(crate) read_comm: Vec<Pcs::Commitment>,
    pub(crate) final_comm: Vec<Pcs::Commitment>,
    // Placeholder for PCS proofs related to `a` opening - handled implicitly by transcript
    pub(crate) placeholder_proofs: PhantomData<Pcs>,
}

/// Lasso proof message round 3
#[derive(Clone, Debug)]
pub struct Message3<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) rz: Point<F, Pcs::Polynomial>, // Point uses Polynomial assoc type
    pub(crate) E_eval: Vec<F>,
    // Placeholder for PCS proofs related to `E` openings at rz
    pub(crate) placeholder_proofs: PhantomData<Pcs>,
}

/// Lasso proof message round 4
#[derive(Clone, Debug)]
pub struct Message4<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) S0_comm: Vec<Pcs::Commitment>,
    pub(crate) S_comm: Vec<Pcs::Commitment>,
    pub(crate) RS_comm: Vec<Pcs::Commitment>,
    pub(crate) WS_comm: Vec<Pcs::Commitment>,
    pub(crate) _marker: PhantomData<F>,
}

/// Data associated with Grand Product checks (PCS openings)
#[derive(Clone, Debug)]
pub struct GrandProductData<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) f_0_r: F,   // Evaluation f(0, r')
    pub(crate) f_1_r: F,   // Evaluation f(1, r')
    pub(crate) f_r_0: F,   // Evaluation f(r', 0)
    pub(crate) f_r_1: F,   // Evaluation f(r', 1)
    pub(crate) product: F, // Evaluation f(1, 0..0) = Root_P
    // Placeholder for PCS proofs - handled implicitly by transcript
    pub(crate) placeholder_proofs: PhantomData<Pcs>,
}

/// Lasso proof message round 5
#[derive(Clone, Debug)]
pub struct Message5<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) r_prime: Vec<Point<F, Pcs::Polynomial>>, // alpha points [b, r'] for S0
    pub(crate) r_prime2: Vec<Point<F, Pcs::Polynomial>>, // alpha points [b, r'] for S
    pub(crate) r_prime3: Vec<Point<F, Pcs::Polynomial>>, // alpha points [b, r'] for RS
    pub(crate) r_prime4: Vec<Point<F, Pcs::Polynomial>>, // alpha points [b, r'] for WS
    pub(crate) S0_data: Vec<GrandProductData<F, Pcs>>,
    pub(crate) S_data: Vec<GrandProductData<F, Pcs>>,
    pub(crate) RS_data: Vec<GrandProductData<F, Pcs>>,
    pub(crate) WS_data: Vec<GrandProductData<F, Pcs>>,
    pub(crate) E_eval2: Vec<F>,    // E_i(r3')
    pub(crate) dim_eval: Vec<F>,   // dim_i(r3')
    pub(crate) read_eval: Vec<F>,  // read_i(r3')
    pub(crate) final_eval: Vec<F>, // final_i(r2')
    // Placeholder for PCS batch proofs - handled implicitly by transcript
    pub(crate) placeholder_proofs: PhantomData<Pcs>,
}

/// The final Lasso proof.
#[derive(Clone, Debug)]
pub struct LassoProof<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) msg1: Message1<F, Pcs>,
    pub(crate) msg2: Message2<F, Pcs>,
    pub(crate) msg3: Message3<F, Pcs>,
    pub(crate) msg4: Message4<F, Pcs>,
    pub(crate) msg5: Message5<F, Pcs>,
}
