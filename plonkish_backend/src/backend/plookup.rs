use crate::{
    pcs::PolynomialCommitmentScheme,
    poly::univariate::UnivariatePolynomial,
    util::{
        arithmetic::PrimeField,
        transcript::{InMemoryTranscript, Keccak256Transcript, TranscriptRead, TranscriptWrite},
        Deserialize, DeserializeOwned, Serialize,
    },
    Error,
};
use halo2_curves::ff::WithSmallOrderMulGroup;
use std::{fmt::Debug, hash::Hash, marker::PhantomData, time::Instant};
use crate::backend::cq::generate_table_and_lookup;

pub mod preprocessor;
pub mod prover;
pub mod util;
pub mod verifier;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlookupProverParam<F, Pcs>
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F>,
{
    pcs: Pcs::ProverParam,
    g: F,
    table: Vec<F>,
    lookup: Vec<F>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlookupVerifierParam<F, Pcs>
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F>,
{
    pcs: Pcs::VerifierParam,
    g: F,
    table: Vec<F>,
}

#[derive(Clone, Debug)]
pub struct PlookupInfo<F> {
    k: u32,
    table: Vec<F>,
    lookup: Vec<F>,
}

#[derive(Clone, Debug)]
pub struct Plookup<F, Pcs>(PhantomData<F>, PhantomData<Pcs>);

impl<F, Pcs> Plookup<F, Pcs>
where
    F: PrimeField + WithSmallOrderMulGroup<3> + Hash + Serialize + DeserializeOwned,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = UnivariatePolynomial<F>>,
{
    fn preprocess(
        param: &Pcs::Param,
        info: &PlookupInfo<F>,
    ) -> Result<(PlookupProverParam<F, Pcs>, PlookupVerifierParam<F, Pcs>), Error> {
        preprocessor::preprocess(param, info)
    }

    fn prove(
        pp: PlookupProverParam<F, Pcs>,
        transcript: &mut (impl TranscriptWrite<Pcs::CommitmentChunk, F> + InMemoryTranscript),
    ) -> Result<(), Error> {
        prover::prove(pp, transcript)
    }

    fn verify(
        vp: PlookupVerifierParam<F, Pcs>,
        transcript: &mut (impl TranscriptRead<Pcs::CommitmentChunk, F> + InMemoryTranscript),
    ) -> Result<(), Error> {
        verifier::verify(vp, transcript)
    }
}

// Concrete implementation for BN256/Fr
use crate::pcs::univariate::UnivariateKzg;
use halo2_curves::bn256::{Bn256, Fr};

impl Plookup<Fr, UnivariateKzg<Bn256>> {
    // Run the full Plookup protocol with given table and lookup
    pub fn test_plookup_by_input(table: Vec<Fr>, lookup: Vec<Fr>) -> Vec<String> {
        let mut timings: Vec<String> = vec![];

        let start_total = Instant::now();

        // Calculate k value based on max of table and lookup size
        let size = std::cmp::max(table.len(), lookup.len()).next_power_of_two();
        let k = (size as f64).log2() as u32;
        let n = 1 << k;

        // 1. Setup
        let info = PlookupInfo {
            k,
            table: table.clone(),
            lookup: lookup.clone(),
        };

        let start = Instant::now();
        let mut rng = crate::util::test::std_rng();
        let param = UnivariateKzg::<Bn256>::setup(n * 4, 1, &mut rng).unwrap();
        let (pp, vp) = Self::preprocess(&param, &info).unwrap();
        let duration1 = start.elapsed();
        timings.push(format!("Setup and preprocess: {}ms", duration1.as_millis()));

        // 2. Prove
        let start = Instant::now();
        let mut transcript = Keccak256Transcript::new(());
        Self::prove(pp, &mut transcript).unwrap();
        let proof = transcript.into_proof();
        let duration2 = start.elapsed();
        timings.push(format!("Prove: {}ms", duration2.as_millis()));

        // 3. Verify
        let start = Instant::now();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        Self::verify(vp, &mut transcript).unwrap();
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));

        timings
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        pcs::univariate::UnivariateKzg,
        util::{test::std_rng, transcript::Keccak256Transcript},
    };
    use halo2_curves::bn256::{Bn256, Fr};

    type Pcs = UnivariateKzg<Bn256>;
    type Pb = Plookup<Fr, Pcs>;

    #[test]
    fn test_e2e() {
        let k = 2;
        let n = 1 << k;
        let lookup = vec![Fr::one(), Fr::one(), Fr::from(2)];
        let table = vec![Fr::one(), Fr::from(2), Fr::from(3), Fr::from(4)];
        let info = PlookupInfo { k, table, lookup };
        let mut rng = std_rng();

        // Setup
        let (pp, vp) = {
            let param = Pcs::setup(n * 4, 1, &mut rng).unwrap();
            Pb::preprocess(&param, &info).unwrap()
        };

        let mut transcript = Keccak256Transcript::new(());
        // Prove
        Pb::prove(pp, &mut transcript).unwrap();
        let proof = transcript.into_proof();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        // Verify
        Pb::verify(vp, &mut transcript).unwrap();
    }

    #[test]
    fn test_plookup_by_input() {
        let table_size = 2_usize.pow(6);
        let lookup_size = 4;

        // Generate table and lookup
        let (table, lookup) = generate_table_and_lookup(table_size, lookup_size);

        // Convert to Fr type since we're using the Bn256 curve
        let table_fr = table.iter().map(|&v| v.clone()).collect();
        let lookup_fr = lookup.iter().map(|&v| v.clone()).collect();

        let timings = Pb::test_plookup_by_input(table_fr, lookup_fr);

        // Print all timing information
        for timing in &timings {
            println!("{}", timing);
        }
    }
}
