use halo2_proofs::transcript;
use rand::{rngs::OsRng, RngCore};
use std::{fmt::Debug, hash::Hash, marker::PhantomData, ops::Div};

use halo2_curves::{bn256::Bn256, ff::WithSmallOrderMulGroup};

use crate::{
    backend::{
        plookup::{
            preprocessor::preprocess,
            prover::prove,
        },
        PlonkishBackend, PlonkishCircuit, PlonkishCircuitInfo,
    },
    poly::Polynomial,
    poly::univariate::UnivariatePolynomial,
    pcs::{
        PolynomialCommitmentScheme,
        univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam},
    },
    util::{
        arithmetic::{Field, PrimeField, MultiMillerLoop, root_of_unity},
        expression::Expression,
        test::std_rng,
        Deserialize, DeserializeOwned, Itertools, Serialize,
        transcript::{InMemoryTranscript, TranscriptRead, TranscriptWrite, Keccak256Transcript},
    },
    Error,
};


pub mod preprocessor;
pub mod prover;
pub mod verifier;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlookupProverParam <F, Pcs>
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F>,
{
    pub(crate) pcs: Pcs::ProverParam,
    pub(crate) g: F,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlookupVerifierParam<F, Pcs>
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F>,
{
    pub(crate) pcs: Pcs::VerifierParam,
    pub(crate) g: F,
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
    ) -> Result<
        (
            PlookupProverParam<F, Pcs>,
            PlookupVerifierParam<F, Pcs>,
        ), 
        Error> {
        preprocess(param, info)
    }

    fn prove(
        pp: PlookupProverParam<F, Pcs>,
        info: &PlookupInfo<F>,
        transcript: &mut (impl TranscriptWrite<Pcs::CommitmentChunk, F> + InMemoryTranscript),
    ) -> Result<(), Error> {
        prove(pp, info, transcript)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::HashSet;
    use halo2_curves::bn256::Fr;
    use crate::util::transcript::{FieldTranscriptRead, FieldTranscriptWrite};
    use num_bigint::BigUint;

    type Pcs = UnivariateKzg<Bn256>;
    type Pb = Plookup<Fr, Pcs>;

    #[test]
    fn test_e2e() {
        let k = 2;
        let n = 1 << k;
        let lookup = vec![Fr::one(), Fr::one(), Fr::from(2)];
        let table = vec![Fr::one(), Fr::from(2), Fr::from(3), Fr::from(4)];
        let info = PlookupInfo{k, table, lookup};
        let mut rng = std_rng();

        // Setup
        let (pp, vp) = {
            let param = Pcs::setup(n, 1, &mut rng).unwrap();
            Pb::preprocess(&param, &info).unwrap()
        };

        let mut transcript = Keccak256Transcript::new(());
        // Prove
        Pb::prove(pp, &info, &mut transcript).unwrap();
    }
}