use halo2_proofs::transcript;
use rand::{rngs::OsRng, RngCore};
use std::{fmt::Debug, hash::Hash, marker::PhantomData, ops::Div};

use halo2_curves::{bn256::Bn256, ff::WithSmallOrderMulGroup};

use crate::{
    backend::{
        plookup::{
            preprocessor::preprocess,
            prover::{sorted_by_table, compute_inner_polynomial, compute_quotient_polynomial, aggregate_poly, aggregate_field},
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
    // pub(crate) lookup: Vec<(Expression<F>, Expression<F>)>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlookupVerifierParam<F, Pcs>
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F>,
{
    pub(crate) pcs: Pcs::VerifierParam,
    pub(crate) g: F,
    // pub(crate) preprocess_comms: Vec<Pcs::Commitment>,
}

#[derive(Clone, Debug)]
pub struct PlookupInfo<F> {
    k: u32,
    table: Vec<F>,
    witness: Vec<F>,
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
        let order = 1 << info.k;
        let table = info.table.clone();
        assert_eq!(table.len(), order);
        let mut witness = info.witness.clone();
        assert!(witness.len() < order);

        // round 1
        // pad witness to length order-1
        let last = witness[witness.len()-1];
        while witness.len() < order-1 {
            witness.push(last);
        }
        let s = sorted_by_table(&table, &witness);
        assert_eq!(s.len(), order*2-1);
        let h1 = s[..order].to_vec();
        let h2 = s[order-1..].to_vec();
        let t_poly = UnivariatePolynomial::lagrange(table.clone());
        let mut witness_clone = witness.clone();
        witness_clone.push(last); // pad to length order for polynomial
        let f_poly = UnivariatePolynomial::lagrange(witness_clone);
        let h1_poly = UnivariatePolynomial::lagrange(h1);
        let h2_poly = UnivariatePolynomial::lagrange(h2);
        // write f_comm, h1_comm, h2_comm
        Pcs::commit_and_write(&pp.pcs, &f_poly, transcript).unwrap();
        Pcs::commit_and_write(&pp.pcs, &h1_poly, transcript).unwrap();
        Pcs::commit_and_write(&pp.pcs, &h2_poly, transcript).unwrap();

        // round 2
        let mut challenges: Vec<F> = Vec::with_capacity(2);
        challenges.extend(transcript.squeeze_challenges(2));
        let beta = &challenges[0];
        let gamma = &challenges[1];
        let z_poly = compute_inner_polynomial(beta, gamma, &witness, &table, &s);
        // write z_comm
        Pcs::commit_and_write(&pp.pcs, &z_poly, transcript).unwrap();

        // round 3
        let delta = &transcript.squeeze_challenge();
        let q_poly = compute_quotient_polynomial(
            beta, gamma, delta, &t_poly, &f_poly, &h1_poly, &h2_poly, &z_poly,
        );
        // write q_comm
        Pcs::commit_and_write(&pp.pcs, &q_poly, transcript).unwrap();

        // round 4
        let zeta = transcript.squeeze_challenge();
        let f_eval = f_poly.evaluate(&zeta);
        let h1_eval = h1_poly.evaluate(&zeta);
        let h2_eval = h2_poly.evaluate(&zeta);
        let z_eval = z_poly.evaluate(&zeta);
        let q_eval = q_poly.evaluate(&zeta);
        transcript.write_field_elements(vec![
            &f_eval, &h1_eval, &h2_eval, &z_eval, &q_eval
        ]).unwrap();
        let g_zeta = pp.g * zeta;
        let h1_g_eval = h1_poly.evaluate(&g_zeta);
        let h2_g_eval = h2_poly.evaluate(&g_zeta);
        let z_g_eval = z_poly.evaluate(&g_zeta);
        transcript.write_field_elements(vec![
            &h1_g_eval, &h2_g_eval, &z_g_eval
        ]).unwrap();

        // round 5
        let eps = &transcript.squeeze_challenge();
        let agg_witness_comm = {
            let agg_poly = aggregate_poly(eps, vec![
                &f_poly, &h1_poly, &h2_poly, &z_poly, &q_poly
            ]);
            let agg_eval = aggregate_field(eps, vec![
                &f_eval, &h1_eval, &h2_eval, &z_eval, &q_eval
            ]);
            let numer = agg_poly + -agg_eval;
            let denom = UnivariatePolynomial::monomial(vec![-zeta, F::ONE]);
            let agg_witness = numer.div(&denom);
            Pcs::commit_and_write(&pp.pcs, &agg_witness, transcript).unwrap();
        };
        let agg_g_witness_comm = {
            let agg_poly = aggregate_poly(eps, vec![
                &h1_poly, &h2_poly, &z_poly
            ]);
            let agg_eval = aggregate_field(eps, vec![
                &h1_g_eval, &h2_g_eval, &z_g_eval
            ]);
            let numer = agg_poly + -agg_eval;
            let denom = UnivariatePolynomial::monomial(vec![-g_zeta, F::ONE]);
            let agg_witness = numer.div(&denom);
            Pcs::commit_and_write(&pp.pcs, &agg_witness, transcript).unwrap();
        };
        Ok(())
        // (*transcript.into_proof()).to_vec()
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

    #[test]
    fn test_e2e() {
        let lookup = vec![Fr::one(), Fr::one(), Fr::from(2)];
        let table = vec![Fr::one(), Fr::from(2), Fr::from(3), Fr::from(4)];
        let m = lookup.len();
        let mut rng = std_rng();

        // Setup
        let (pp, vp) = {
            let mut rng = OsRng;
            let poly_size = m;
            print!("poly_size: {:?}\n", poly_size);
            let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
            Pcs::trim(&param, poly_size, 1).unwrap()
        };
        print!("lookup: {:?}\n", lookup);

        // Commit and open
        // let proof = {
        //     let mut transcript = Keccak256Transcript::new(());
        //     let poly = <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::monomial(lookup.clone());
        //     print!("coeffs: {:?}\n", poly.coeffs());
        //     let comm = Pcs::commit_and_write(&pp, &poly, &mut transcript).unwrap();
        //     let point = <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::squeeze_point(m, &mut transcript);
        //     let eval = poly.evaluate(&point);

        //     // Use the correct method to write the field element
        //     transcript.write_field_element(&eval).unwrap();

        //     Pcs::open(&pp, &poly, &comm, &point, &eval, &mut transcript).unwrap();
        //     transcript.into_proof()
        // };

        // Verify
        // let result = {
        //     let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        //     Pcs::verify(
        //         &vp,
        //         &Pcs::read_commitment(&vp, &mut transcript).unwrap(),
        //         &<Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::squeeze_point(m, &mut transcript),

        //         // Use the correct method to read the field element
        //         &transcript.read_field_element().unwrap(),
        //         &mut transcript,
        //     )
        // };
        // assert_eq!(result, Ok(()));

    }
}