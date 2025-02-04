use rand::rngs::OsRng;
use std::{fmt::Debug, marker::PhantomData};

use halo2_curves::{bn256::{ pairing, Bn256, Fr, G1Affine, G2Affine, G2Prepared, Gt, G1, G2}, pairing::MillerLoopResult};

use crate::{
    poly::Polynomial,
    poly::univariate::UnivariatePolynomial,
    backend::cq::preprocessor::preprocess,
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
    use std::time::Instant;
    #[test]
    fn test_cq() {
        let mut table = vec![];
        for k in 1..=2_usize.pow(6) {
            table.push(Fr::from(k as u64));
        }
        // let table = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4), Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8)];
        let lookup = vec![Fr::from(4), Fr::from(3), Fr::from(5), Fr::from(2)];
    
        let m = lookup.len();
        let t = table.len();
        let poly_size = max(t, m).next_power_of_two() * 2;

        // 1. setup
        let start = Instant::now();
        let (param, pp, vp, q_t_comm_poly_coeffs) = preprocess(t, m, &table).unwrap();
        assert_eq!(poly_size, 2_usize.pow(pp.k() as u32));
        let duration1 = start.elapsed();
        // timings.push(format!("k={k}, setup and preprocess time: {}ms",duration1.as_millis()));
        println!("\n ------------Setup and preprocess: {}ms----------- \n",duration1.as_millis());

        // 2. generate proof
        let start = Instant::now();
        let prover = Prover::new(&table, &param, &pp);
        let proof = prover.prove(&lookup, &q_t_comm_poly_coeffs);
        println!("proof: {:?}", proof);
        let duration2 = start.elapsed();
        // timings.push(format!("k={k}, setup and preprocess time: {}ms",duration1.as_millis()));
        println!("\n ------------prove: {}ms----------- \n",duration2.as_millis());


        let scalar_0 = Fr::from(0 as u64);
        let scalar_1 = Fr::from(1 as u64);

        // 3 verifier to verify
        // 3.1 prepare for verifier
        // z_v(x) = X^N - 1, [-1, 0, ..., 0, 1], t-1 0s in between
        let z_v_poly_coeffs = vec![scalar_1.neg()].into_iter().chain(vec![scalar_0; t - 1]).chain(vec![scalar_1]).collect();
        let z_v_poly = UnivariatePolynomial::monomial(z_v_poly_coeffs);
        // [z_h(x)]2
        let z_v_comm_2 = Pcs::commit_monomial_g2(&param, &z_v_poly.coeffs());
        // t(x)
        let t_poly = UnivariatePolynomial::lagrange(table.clone()).ifft();
        // [t(x)]2
        let t_comm_2 = Pcs::commit_monomial_g2(&param, &t_poly.coeffs());

        // [X^{N - 1 - (n - 2)}]2
        // x_exponent_poly_comm_2
        let x_exponent_order = t - 1 - (m - 2);
        let x_exponent_values_in_coeff = vec![scalar_0; x_exponent_order].into_iter().chain(vec![scalar_1]).collect();
        let x_exponent_poly = UnivariatePolynomial::monomial(x_exponent_values_in_coeff);
        // commit x_exponent_poly
        let x_exponent_poly_comm_2 = Pcs::commit_monomial_g2(&param, &x_exponent_poly.coeffs());

        // 3.2 verifier to verify
        let start = Instant::now();
        let verifier = Verifier::new(&vp);
        verifier.verify(
            &proof,
            &t_comm_2,
            &z_v_comm_2,
            &x_exponent_poly_comm_2,
            m,
            t
        );
        let duration3 = start.elapsed();
        println!("\n ------------verify: {}ms----------- \n",duration3.as_millis());

        println!("Finished to verify: cq");
    }

}
