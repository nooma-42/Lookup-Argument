use rand::rngs::OsRng;
use std::{fmt::Debug, marker::PhantomData, collections::HashSet, ops::Mul};

use halo2_curves::{bn256::{multi_miller_loop, pairing, Bn256, Fr, G1Affine, G2Affine, G2Prepared, Gt, G1, G2}, pairing::MillerLoopResult};

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
        transcript::{InMemoryTranscript, TranscriptRead, TranscriptWrite, Keccak256Transcript, FieldTranscript, FieldTranscriptRead, FieldTranscriptWrite, G2TranscriptRead, G2TranscriptWrite},
    }
};


pub fn lagrange_interp(h_i_values: &[Fr], t_values_from_lookup: &[Fr]) -> UnivariatePolynomial<Fr> {
    assert!(h_i_values.len() == t_values_from_lookup.len());

    let vanishing_poly = UnivariatePolynomial::vanishing(h_i_values, Fr::one());
    let mut bary_centric_weights = vec![Fr::one(); h_i_values.len()];
    let mut sum = UnivariatePolynomial::monomial(vec![Fr::zero()]);
    for (idx, h_i) in h_i_values.iter().enumerate() {
        for (jdx, h_j) in h_i_values.iter().enumerate() {
            if jdx == idx {
                continue;
            }
            bary_centric_weights[idx] = bary_centric_weights[idx] * (h_i - h_j).invert().unwrap();
        }
        let y_i = t_values_from_lookup[idx];
        // x - x_i
        let v_poly = UnivariatePolynomial::monomial(vec![-h_i, Fr::one()]);
        let (v_poly_inv, _) = vanishing_poly.div_rem(&v_poly);
        let accu = &v_poly_inv * (y_i * bary_centric_weights[idx]);
        sum += accu;
    }
    sum
}

pub fn multi_pairing(g1: &[G1Affine], g2: &[G2Affine]) -> Gt {
    assert_eq!(g1.len(), g2.len(), "Input slices must have the same length");

    let g2_prepared: Vec<G2Prepared> = g2.iter().map(|&g| G2Prepared::from_affine(g)).collect();
    let terms: Vec<(&G1Affine, &G2Prepared)> = g1.iter().zip(g2_prepared.iter()).collect();

    let u = multi_miller_loop(&terms);
    u.final_exponentiation()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{collections::HashSet, ops::{Add, Mul}};
    use halo2_curves::bn256::Fr;
    use num_integer::Roots;
    use crate::util::transcript::{FieldTranscript, FieldTranscriptRead, FieldTranscriptWrite, G2TranscriptRead, G2TranscriptWrite};
    type Pcs = UnivariateKzg<Bn256>;

    #[test]
    fn test_interpolation_poly() {
        let h_i = vec![Fr::from(2), Fr::from(3), Fr::from(4)];
        let t_values_from_lookup = vec![Fr::from(4), Fr::from(9), Fr::from(16)];
        // f(x) = x^2
        let poly = lagrange_interp(&h_i, &t_values_from_lookup);
        assert_eq!(poly.coeffs(), vec![Fr::from(0), Fr::from(0), Fr::from(1)]);
        assert_eq!(poly.evaluate(&Fr::from(5)), Fr::from(25));
    }


    #[test]
    fn test_multi_pairing_only_g2() {
        let mut rng = OsRng;
        let g1 = G1::generator();
        let g2 = G2::generator();
        let g1_affine = G1Affine::from(g1);
        let g2_affine = G2Affine::from(g2);
        let srs_size = 1 << 8;
        print!("srs_size: {:?}\n", srs_size);
        let param = Pcs::setup(srs_size, 1, &mut rng).unwrap();
        let (pp, vp) = Pcs::trim(&param, srs_size, 1).unwrap();
        let mut transcript = Keccak256Transcript::new(());
        let alpha = Fr::from(2 as u64);
        let scalar_1 = Fr::from(1 as u64);
        let test_poly = UnivariatePolynomial::monomial(vec![Fr::from(1 as u64), Fr::from(2 as u64), Fr::from(3 as u64), Fr::from(4 as u64)]);
        let x_alpha_poly = UnivariatePolynomial::monomial(vec![-alpha, scalar_1]);
        let test_comm_1 = Pcs::commit_and_write(&pp, &test_poly, &mut transcript).unwrap();
        print!("test_poly_comm_1: {:?}\n", test_comm_1);
        let test_comm_1_affine: G1Affine = test_comm_1.to_affine();
        let test_at_alpha = test_poly.evaluate(&alpha);
        let quotient_poly = &(test_poly + test_at_alpha.neg()) / &x_alpha_poly;
        print!("quotient_poly: {:?}\n", quotient_poly);
        let quotient_comm_1 = Pcs::commit_and_write(&pp, &quotient_poly, &mut transcript).unwrap();
        print!("quotient_comm_1: {:?}\n", quotient_comm_1);
        let quotient_comm_1_affine: G1Affine = quotient_comm_1.to_affine();
        let lhs = pairing(&quotient_comm_1_affine, &vp.s_g2());
        let rhs_affine = test_comm_1_affine.add(&quotient_comm_1_affine.mul(alpha)).add(g1_affine.mul(test_at_alpha.neg())).into();
        let rhs = pairing(&rhs_affine, &g2_affine);
        assert_eq!(lhs, rhs);
        let rhs_terms_g1 = vec![test_comm_1_affine, quotient_comm_1_affine.mul(alpha).into(), g1_affine.mul(test_at_alpha.neg()).into()];
        let rhs_terms_g2 = vec![g2_affine, g2_affine, g2_affine];
        let rhs = multi_pairing(&rhs_terms_g1, &rhs_terms_g2);
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn test_multi_pairing_with_s_g2() {
        let mut rng = OsRng;
        let g1 = G1::generator();
        let g2 = G2::generator();
        let g1_affine = G1Affine::from(g1);
        let g2_affine = G2Affine::from(g2);
        let srs_size = 1 << 8;
        print!("srs_size: {:?}\n", srs_size);
        let param = Pcs::setup(srs_size, 1, &mut rng).unwrap();
        let (pp, vp) = Pcs::trim(&param, srs_size, 1).unwrap();
        let mut transcript = Keccak256Transcript::new(());
        let alpha = Fr::from(2 as u64);
        let scalar_0 = Fr::from(0 as u64);
        let scalar_1 = Fr::from(1 as u64);
        let test_poly = UnivariatePolynomial::monomial(vec![Fr::from(1 as u64), Fr::from(2 as u64), Fr::from(3 as u64), Fr::from(4 as u64)]);
        let x_poly = UnivariatePolynomial::monomial(vec![scalar_0, scalar_1]);
        let test_x_poly = test_poly.poly_mul(x_poly.clone());
        let x_alpha_poly = UnivariatePolynomial::monomial(vec![-alpha, scalar_1]);

        let test_comm_1 = Pcs::commit_and_write(&pp, &test_poly, &mut transcript).unwrap();
        let test_comm_1_affine: G1Affine = test_comm_1.to_affine();
        let test_x_comm_1 = Pcs::commit_and_write(&pp, &test_x_poly, &mut transcript).unwrap();
        let test_x_comm_1_affine: G1Affine = test_x_comm_1.to_affine();

        let test_x_at_alpha = test_x_poly.evaluate(&alpha);
        let x_quotient_poly = &(test_x_poly + test_x_at_alpha.neg()) / &x_alpha_poly;

        let x_quotient_comm_1 = Pcs::commit_and_write(&pp, &x_quotient_poly, &mut transcript).unwrap();
        let x_quotient_comm_1_affine: G1Affine = x_quotient_comm_1.to_affine();

        let lhs = pairing(&x_quotient_comm_1_affine, &vp.s_g2());
        let rhs_affine = test_x_comm_1_affine.add(&x_quotient_comm_1_affine.mul(alpha)).add(g1_affine.mul(test_x_at_alpha.neg())).into();
        let rhs = pairing(&rhs_affine, &g2_affine);
        assert_eq!(lhs, rhs);

        let rhs_terms_g1 = vec![test_comm_1_affine, x_quotient_comm_1_affine.mul(alpha).into(), g1_affine.mul(test_x_at_alpha.neg()).into()];
        let rhs_terms_g2 = vec![vp.s_g2(), g2_affine, g2_affine];
        let rhs = multi_pairing(&rhs_terms_g1, &rhs_terms_g2);
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn test_pairing() {
        let mut rng = OsRng;
        let a = Fr::random(&mut rng);
        let b = Fr::random(&mut rng);

        let mut g1 = G1::generator();
        g1 = g1.mul(a);

        let mut g2 = G2::generator();
        g1 = g1.mul(b);
        let pair_ab = pairing(&G1Affine::from(g1), &G2Affine::from(g2));

        g1 = G1::generator();
        g1 = g1.mul(b);

        g2 = G2::generator();
        g1 = g1.mul(a);

        let pair_ba = pairing(&G1Affine::from(g1), &G2Affine::from(g2));

        assert_eq!(pair_ab, pair_ba);
        println!("pairing: {:?}", pair_ab);
    }
}
