use rand::rngs::OsRng;
use std::{fmt::Debug, marker::PhantomData};

use halo2_curves::{bn256::{ pairing, Bn256, Fr, G1Affine, G2Affine, G2Prepared, Gt, G1, G2}, pairing::MillerLoopResult};

use crate::{
    poly::Polynomial,
    poly::univariate::UnivariatePolynomial,
    backend::baloo::preprocessor::preprocess,
    pcs::{
        PolynomialCommitmentScheme,
        Additive,
        univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam, UnivariateKzgCommitment},
    },
    util::{
        arithmetic::{Field, PrimeField, root_of_unity, variable_base_msm, barycentric_weights},
        test::std_rng,
        transcript::{InMemoryTranscript, TranscriptRead, TranscriptWrite, Keccak256Transcript},
    }
};

pub mod preprocessor;
pub mod prover;
pub mod verifier;
pub mod util;

use prover::Prover;
use verifier::Verifier;

#[cfg(test)]
mod tests {
    use super::*;
    use halo2_curves::bn256::Fr;
    use crate::util::transcript::{FieldTranscriptRead, FieldTranscriptWrite, G2TranscriptRead, G2TranscriptWrite};
    type Pcs = UnivariateKzg<Bn256>;
    use std::cmp::max;
    #[test]
    fn test_baloo() {
        let table = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4), Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8)];
        let lookup = vec![Fr::from(4), Fr::from(3), Fr::from(5), Fr::from(2)];

        let m = lookup.len();
        let t = table.len();
        let poly_size = max(t, m).next_power_of_two() * 2;
        let d = poly_size - 2;

        // 1. setup
        let (param, pp, vp) = preprocess(t, m).unwrap();
        assert_eq!(poly_size, 2_usize.pow(pp.k() as u32));

        // 2. generate proof
        let prover = Prover::new(&table, &param, &pp);
        let proof = prover.prove(&lookup);
        println!("proof: {:?}", proof);

        let scalar_0 = Fr::from(0 as u64);
        let scalar_1 = Fr::from(1 as u64);

        // 3.1 prepare for verifier
        // z_h(x) = X^t - 1, [-1, 0, ..., 0, 1], t-1 0s in between
        let z_h_poly_coeffs = vec![scalar_1.neg()].into_iter().chain(vec![scalar_0; t - 1]).chain(vec![scalar_1]).collect();
        let z_h_poly = UnivariatePolynomial::monomial(z_h_poly_coeffs);
        // [z_h(x)]1
        let z_h_comm_1 = Pcs::commit_monomial(&pp, &z_h_poly.coeffs());
        // t(x)
        let t_poly = UnivariatePolynomial::lagrange(table.clone()).ifft();
        // [t(x)]1
        let t_comm_1 = Pcs::commit_monomial(&pp, &t_poly.coeffs());

        // φ(x)
        let phi_poly = UnivariatePolynomial::lagrange(lookup.clone()).ifft();
        let phi_comm_1 = Pcs::commit_monomial(&pp, &phi_poly.coeffs());
        // todo: cached all [x^s]1, [x^s]2?
        // X^m
        let x_m_exponent_poly = UnivariatePolynomial::monomial(vec![scalar_0; m].into_iter().chain(vec![scalar_1]).collect());
        // [X^m]1
        let x_m_exponent_poly_comm_1 = Pcs::commit_monomial(&pp, &x_m_exponent_poly.clone().coeffs());

        // X^(d-m+1)
        let coeffs_x_exponent_poly = vec![scalar_0; d - m + 1].into_iter().chain(vec![scalar_1]).collect();
        let x_exponent_poly = UnivariatePolynomial::monomial(coeffs_x_exponent_poly);
        // [X^(d-m+1)]2
        let x_exponent_poly_comm_2 = Pcs::commit_monomial_g2(&param, &x_exponent_poly.coeffs());
        println!("x_exponent_poly_comm_2: {:?}", x_exponent_poly_comm_2);

        // X^(d-m+2)
        let coeffs_x_exponent_poly_2 = vec![scalar_0; d - m + 2].into_iter().chain(vec![scalar_1]).collect();
        let x_exponent_poly_2 = UnivariatePolynomial::monomial(coeffs_x_exponent_poly_2);
        // [X^(d-m+2)]1
        let x_exponent_poly_2_comm_1 = Pcs::commit_monomial(&pp, &x_exponent_poly_2.coeffs());
        // [X^(d-m+2)]2
        let x_exponent_poly_2_comm_2 = Pcs::commit_monomial_g2(&param, &x_exponent_poly_2.coeffs());

        // 3.2 verifier to verify
        let verifier = Verifier::new(&vp);
        verifier.verify(
            &proof,
            &t_comm_1,
            &z_h_comm_1,
            &phi_comm_1,
            &x_m_exponent_poly_comm_1,
            &x_exponent_poly_comm_2,
            &x_exponent_poly_2_comm_1,
            &x_exponent_poly_2_comm_2,
            m
        );

        println!("Finished to verify: baloo");
    }

}
