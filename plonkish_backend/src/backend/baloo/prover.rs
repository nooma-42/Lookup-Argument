use rand::rngs::OsRng;
use num_integer::Roots;
use std::{fmt::Debug, collections::HashSet, ops::Mul};
use halo2_curves::{bn256::{pairing, Bn256, Fr, G1Affine, G2Affine, G1, G2}, pairing::MillerLoopResult};
use crate::{
    poly::{Polynomial, univariate::UnivariatePolynomial},
    backend::baloo::{
      util::{lagrange_interp, multi_pairing},
      preprocessor::preprocess,
    },
    pcs::{
        PolynomialCommitmentScheme,
        Additive,
        univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam, UnivariateKzgCommitment},
    },
    util::{
        arithmetic::{Field, PrimeField, root_of_unity, variable_base_msm, barycentric_weights},
        transcript::{InMemoryTranscript, TranscriptWrite, Keccak256Transcript, FieldTranscript, FieldTranscriptRead, FieldTranscriptWrite, G2TranscriptRead, G2TranscriptWrite},
    }
};

type Pcs = UnivariateKzg<Bn256>;

pub struct Baloo<F> {
    table: Vec<F>,
}

impl Baloo<Fr>
{
    pub fn new(table: Vec<Fr>) -> Self {
        Baloo { table }
    }
    /*
    How to calculate ξ(x)(or v(x) in code)

    H = [1, ω, ω^2, ω^3, ω^4, ω^5, ω^6, ω^7]
    table = [1, 2, 3, 4, 5, 6, 7, 8]
    lookup = [3, 7, 3, 4]
    t_I:  [3, 4, 7] # choose non-duplicated elements from lookup
    I:  [2, 3, 6] # get indexes from table
    s ∈ I = [2, 3, 6]
    H_I = {ω^s} = [ω^2, ω^3, ω^6]
    k = len(I) = 3
    vanishing polynomial z_I(X) = (X - H_I_0)(X - H_I_1)(X - H_I_2)
                                = (X - ω^2)(X - ω^3)(X - ω^6)
    M * t = lookup
    M = [
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
    ]
    m = len(lookup) # 4
    col[0] = M[0].index(1)
    col[1] = M[1].index(1)
    col[2] = M[2].index(1)
    col[3] = M[3].index(1)
    col: [2, 6, 2, 3]
    xi = H_I[col_i] = [ω^2, ω^6, ω^2, ω^3]
    Interpolation with xi and get polynomial: ξ(x)
    */
    pub fn round1(&self, lookup: Vec<Fr>) -> Vec<Fr>
    {
        let table = self.table.clone();
        let m = lookup.len();
        let t = table.len();

        // Setup
        let mut rng = OsRng;
        let poly_size = m;
        print!("poly_size: {:?}\n", poly_size);
        // TODO: proper srs size
        let srs_size = 1 << 8;
        print!("srs_size: {:?}\n", srs_size);
        let param = Pcs::setup(srs_size, 1, &mut rng).unwrap();
        let (pp, vp) = Pcs::trim(&param, srs_size, 1).unwrap();
        print!("lookup: {:?}\n", lookup);

        // Commit and open
        let mut transcript = Keccak256Transcript::new(());
        // commit ξ(x) on G1
        let phi_poly = UnivariatePolynomial::lagrange(lookup.clone()).ifft();
        print!("coeffs: {:?}\n", phi_poly.coeffs());
        let phi_comm_1 = Pcs::commit_and_write(&pp, &phi_poly, &mut transcript).unwrap();

        // remove duplicated elements
        let t_values_from_lookup: HashSet<_> = lookup.clone().into_iter().collect();
        print!("t_values_from_lookup: {:?}\n", t_values_from_lookup);
        // I: the index of t_values_from_lookup elements in sub table t_I
        let i_values: Vec<_> = t_values_from_lookup.iter().map(|elem| table.iter().position(|&x| x == *elem).unwrap()).collect();
        print!("i_values: {:?}\n", i_values);
        let log_m = m.sqrt();
        let v_root_of_unity = root_of_unity::<Fr>(log_m);

        // H_I = {ξ_i} , i = [1...k], ξ(Xi)
        let h_i: Vec<_> = i_values.iter().map(|&i| {
            let i_as_u64 = i as u64;
            println!("i_as_u64: {:?}", i_as_u64);
            v_root_of_unity.pow([i_as_u64])
        }).collect();
        print!("h_i: {:?}\n", h_i);
        // TODO: optimize interpolation polynomial with https://github.com/gy001/hypercube/blob/main/univarization/src/unipoly.rs#L391
        // refer to barycentric_weights in arithmetic.rs
        let t_values_from_lookup_set : Vec<Fr>= t_values_from_lookup.clone().into_iter().collect();
        print!("t_values_from_lookup_set: {:?}\n", t_values_from_lookup_set);
        let t_i_poly = lagrange_interp(&h_i, &t_values_from_lookup_set);

        let z_i_poly = UnivariatePolynomial::vanishing(&h_i, Fr::one());
        assert_eq!(z_i_poly.evaluate(&h_i[0]), Fr::zero());
        assert_eq!(z_i_poly.evaluate(&h_i[1]), Fr::zero());
        assert_eq!(z_i_poly.evaluate(&h_i[2]), Fr::zero());
        let z_i_comm_2 = Pcs::commit_monomial_g2(&param, &z_i_poly.coeffs());
        transcript.write_commitment_g2(&z_i_comm_2.clone().to_affine()).unwrap();
        print!("z_i_comm_2: {:?}\n", z_i_comm_2);

        let mut col_values = Vec::new();
        let mut v_values = Vec::new();
        for i in 0..m {
            // find the index of 1 in jth row of M
            let col_i = t_values_from_lookup.iter().position(|&x| x == lookup[i]).unwrap();
            col_values.push(col_i);
            let col_i_root = h_i[col_i];
            // Note: v = 1 / col_i_root in paper
            // Here we use different construction that does not affect the verification
            let v = col_i_root;
            v_values.push(v);
        }
        // ξ(x) polynomial
        let v_poly = UnivariatePolynomial::lagrange(v_values.clone()).ifft();
        let v_comm_1 = Pcs::commit_and_write(&pp, &v_poly, &mut transcript).unwrap();
        print!("v_comm_1: {:?}\n", v_comm_1);

        // round 2
        let alpha = transcript.squeeze_challenge();
        let beta = transcript.squeeze_challenge();
        // let z_i_poly = self::z_i_poly.clone();
        // let v_poly = self::v_poly.clone();
        // let t_i_poly = self::t_i_poly.clone();
        // let col = self::col.clone();
        // let h_i = self::h_i.clone();
        // let phi_poly = self::phi_poly.clone();
        // let m = self::m;

        let scalar_0 = Fr::from(0);
        let scalar_1 = Fr::from(1);

        let zero_poly = UnivariatePolynomial::monomial(vec![scalar_0]);

        // [-1, 0, 0, ..., 1], m - 1 0s in between
        let z_v_values: Vec<Fr> = vec![scalar_1.neg()]
            .into_iter()
            .chain((0..m - 1).map(|_| scalar_0))
            .chain(vec![scalar_1])
            .collect();
        assert_eq!(z_v_values.len(), m + 1);
        // X^m - 1
        let z_v_poly = UnivariatePolynomial::monomial(z_v_values);
        assert_eq!(z_v_poly.evaluate(&v_root_of_unity), scalar_0);

        // z_I(0)
        let z_i_at_0 = z_i_poly.evaluate(&scalar_0);
        // calculate D(X) = Σ_{0, m-1} μ_i(α) * τ^_{col(i)}(X)
        let mut d_poly: UnivariatePolynomial<Fr> = zero_poly.clone();
        let mut e_poly: UnivariatePolynomial<Fr> = zero_poly.clone();
        for i in 0..m {
            // col(i)
            let col_i = col_values[i];
            // ω^i
            let v_root = v_root_of_unity.pow([i as u64]);
            print!("v_root: {:?}\n", v_root);
            // X - ω^i
            let v_root_poly = UnivariatePolynomial::monomial(vec![-v_root, scalar_1]);
            // ξ_i
            let col_i_root = h_i[col_i];
            // X - ξ_i
            let x_root_poly = UnivariatePolynomial::monomial(vec![-col_i_root, scalar_1]);
            // Lagrange polynomial on V: μ_i(X)
            // z_v_poly / v_root_poly * v_root / Fr::from(m as u64);
            let mu_poly = &(&z_v_poly / &v_root_poly) * (v_root * (Fr::from(m as u64).invert().unwrap()));
            // Normalized Lagrange Polynomial: τ_col(i)(X) / τ_col(i)(0)
            // z_i_poly / x_root_poly * (-col_i_root) / z_i_at_0;
            let normalized_lag_poly = &(&z_i_poly / &x_root_poly) * (col_i_root.neg() * (z_i_at_0.invert().unwrap()));
            // μ_i(α)
            let mu_poly_at_alpha = mu_poly.evaluate(&alpha);
            // Normalized Lagrange Polynomial at β: τ_col(i)(β) / τ_col(i)(0)
            let normalized_lag_poly_at_beta = normalized_lag_poly.evaluate(&beta);
            // D(X) = Σ_i(μ_i(α) * normalized_lag_poly)
            d_poly += &normalized_lag_poly * mu_poly_at_alpha;
            // E(X) = Σ_i(μ_i(X) * normalized_lag_poly(β))
            e_poly += &mu_poly * normalized_lag_poly_at_beta;
        }
        print!("d_poly: {:?}\n", d_poly);
        print!("e_poly: {:?}\n", e_poly);

        // D(X) * t_I(X)
        let d_t_poly = d_poly.poly_mul(t_i_poly.clone());
        // φ(α)
        let phi_poly_at_alpha = phi_poly.evaluate(&alpha);
        print!("d_t_poly: {:?}\n", d_t_poly);
        print!("phi_poly_at_alpha: {:?}\n", phi_poly_at_alpha);

        // Q_D(X), R(X) = (D(X) * t_I(X) - φ(α)) / z_I(X)
        let (q_d_poly, r_poly) = (d_t_poly + phi_poly_at_alpha.neg()).div_rem(&z_i_poly);
        print!("q_d_poly: {:?}\n", q_d_poly);
        print!("r_poly: {:?}\n", r_poly);
        assert_eq!(r_poly.evaluate(&scalar_0), scalar_0);

        let z_i_at_beta = z_i_poly.evaluate(&beta);
        // Q_E(X) = (E(X) * (β - v(X)) + v(X) * z_I(β) / z_I(0)) / z_V(X)
        // (e_poly * (beta - v_poly) + v_poly * z_i_at_beta / z_i_at_0) / z_v_poly;
        let q_e_poly = {
            // let beta_poly = UnivariatePolynomial::lagrange(vec![beta; m]);
            // beta - v_poly
            let aaa: UnivariatePolynomial<Fr> = &v_poly * scalar_1.neg() + beta;
            // e_poly * (beta - v_poly)
            let bbb = &e_poly.poly_mul(aaa);
            // z_i_at_beta / z_i_at_0
            let ccc = z_i_at_beta.mul(z_i_at_0.invert().unwrap());
            // v_poly * z_i_at_beta / z_i_at_0
            let ddd = &v_poly * ccc;
            // e_poly * (beta - v_poly) + v_poly * z_i_at_beta / z_i_at_0
            &(bbb + ddd) / &z_v_poly
        };
        // e_poly.mul(aaa.add(v_poly.mul(z_i_at_beta.div(z_i_at_0)))).div(z_v_poly);
        print!("q_e_poly: {:?}\n", q_e_poly);

        // π2 = ([D]1 = [D(x)]1, [R]1 = [R(x)]1, [Q2]1 = [Q_D(x)]1)
        let t_i_comm_1 = Pcs::commit_and_write(&pp, &t_i_poly, &mut transcript).unwrap();
        let d_comm_1 = Pcs::commit_and_write(&pp, &d_poly, &mut transcript).unwrap();
        let r_comm_1 = Pcs::commit_and_write(&pp, &r_poly, &mut transcript).unwrap();
        let q_d_comm_1 = Pcs::commit_and_write(&pp, &q_d_poly, &mut transcript).unwrap();
        // π3 = ([E]1 = [E(x)]1, [Q1]1 = [Q_E(x)]1)
        let e_comm_1 = Pcs::commit_and_write(&pp, &e_poly, &mut transcript).unwrap();
        let q_e_comm_1 = Pcs::commit_and_write(&pp, &q_e_poly, &mut transcript).unwrap();

        print!("d_comm_1: {:?}\n", d_comm_1);
        print!("r_comm_1: {:?}\n", r_comm_1);
        print!("q_d_comm_1: {:?}\n", q_d_comm_1);
        print!("e_comm_1: {:?}\n", e_comm_1);
        print!("q_e_comm_1: {:?}\n", q_e_comm_1);

        // self::z_V_poly = z_V_poly;
        // self::D_poly = D_poly;
        // self::E_poly = E_poly;
        // self::R_poly = R_poly;
        // self::Q_D_poly = Q_D_poly;
        // self::Q_E_poly = Q_E_poly;

        // TODO
        // return message2(d_comm_1, r_comm_1, q_d_comm_1, e_comm_1, q_e_comm_1);

        // round 3
        let gamma: Fr = transcript.squeeze_challenge();
        let zeta: Fr = transcript.squeeze_challenge();
        let gamma_2 = gamma.mul(gamma);
        let gamma_3 = gamma_2.mul(gamma);

        let d = srs_size - 2;
        // calculate v1, v2, v3, v4, v5
        // v1 = e(α)
        let v1 = e_poly.evaluate(&alpha);
        // v2 = φ(α)
        let v2 = phi_poly.evaluate(&alpha);
        // v3 = z_I(0)
        let v3 = z_i_at_0;
        // v4 = z_I(β)
        let v4 = z_i_poly.evaluate(&beta);
        // v5 = e(ζ)
        let v5 = e_poly.evaluate(&zeta);

        // calculate D(β)
        let d_beta = d_poly.evaluate(&beta);
        // z_V(ζ)
        let z_v_zeta = z_v_poly.evaluate(&zeta);
        // beta - v_poly
        let beta_sub_v_poly = &v_poly * scalar_1.neg() + beta;

        // P_D(X) = D(β) * t_I(X) - φ(α) - R(X) - z_I(β) * Q_D(X)
        let p_d_poly = &(&t_i_poly * d_beta + v2.neg()) - (&r_poly + &q_d_poly * v4);
        // P_E(X) = E(ζ) * (β - v(X)) + v(X) * z_I(β) / z_I(0) - z_V(ζ) * Q_E(X)
        let p_e_poly = &(&(&beta_sub_v_poly * v5) + &v_poly * (v4.mul(v3.invert().unwrap()))) - &q_e_poly * z_v_zeta;
        // X^(d-m+1)
        let coeffs = vec![scalar_0; d - m + 1].into_iter().chain(vec![scalar_1]).collect();
        let x_exponent_poly = UnivariatePolynomial::monomial(coeffs);
        // calculate [w1]1, [w2]1, [w2]1, [w4]1
        // X - α
        let x_alpha_poly = UnivariatePolynomial::monomial(vec![-alpha, scalar_1]);
        // calculate w1 = X^(d-m+1) * (E(X) - E(α) + (φ(X) - φ(α))γ) / X - α
        let mut w1 = &(&(e_poly.clone() + v1.neg()) + &(phi_poly.clone() + v2.neg()) * gamma) / &x_alpha_poly;
        w1 = w1.poly_mul(x_exponent_poly.clone());
        // calculate polynomial X
        let x_poly = UnivariatePolynomial::monomial(vec![scalar_0, scalar_1]);
        // X^m
        let x_m_exponent_poly = UnivariatePolynomial::monomial(vec![scalar_0; m].into_iter().chain(vec![scalar_1]).collect());
        // calculate w2 = (z_I(X) - v3 / X + γ * R(X) / X +  γ^2 * X^(d-m+1) * (z_I(X) - X^m)) + γ^3 * X^(d-m+1) * R(X)
        let w2 = &(
                &(&(z_i_poly.clone() + v3.neg()) / &x_poly)
                + &(&(&r_poly * gamma) / &x_poly)
            ) + &x_exponent_poly.poly_mul(
                &(&(&z_i_poly - x_m_exponent_poly.clone()) * gamma_2)
                + &r_poly * gamma_3
            );
        print!("w2: {:?}\n", w2.degree());
        // let w2 = (z_I_poly - v3) / x_poly + R_poly * gamma / x_poly + x_exponent_poly * ( (z_I_poly - x_m_exponent_poly) * gamma ** 2 + R_poly * gamma ** 3 )
        // calculate X - β
        let x_beta_poly = UnivariatePolynomial::monomial(vec![-beta, scalar_1]);
        // calculate w3 = (D(X) - E(α) + (z_I(X) - z_I(β))γ + P_D(X)γ^2) / X - β
        // v1 = E(α) == D(β)
        let w3 = &(&(&(d_poly + v1.neg()) + &(z_i_poly.clone() + v4.neg()) * gamma) + &p_d_poly * gamma_2) / &x_beta_poly;
        // calculate X - ζ
        let x_zeta_poly = UnivariatePolynomial::monomial(vec![-zeta, scalar_1]);
        // calculate w4 = (E(X) - E(ζ) + P_E(X)γ) / X - ζ
        // v5 = E(ζ)

        let w4 = &(&(e_poly + v5.neg()) + &p_e_poly * gamma) / &x_zeta_poly;
        // let w44 = &p_e_poly / &x_zeta_poly;
        // let w4 = &(e_poly + v5.neg()) / &x_zeta_poly;

        // calculate [w1]1, [w2]1, [w2]1, [w4]1
        let w1_comm_1 = Pcs::commit_and_write(&pp, &w1, &mut transcript).unwrap();
        let w2_comm_1 = Pcs::commit_and_write(&pp, &w2, &mut transcript).unwrap();
        let w3_comm_1 = Pcs::commit_and_write(&pp, &w3, &mut transcript).unwrap();
        let w4_comm_1 = Pcs::commit_and_write(&pp, &w4, &mut transcript).unwrap();

        print!("w1_comm_1: {:?}\n", w1_comm_1);
        print!("w2_comm_1: {:?}\n", w2_comm_1);
        print!("w3_comm_1: {:?}\n", w3_comm_1);
        print!("w4_comm_1: {:?}\n", w4_comm_1);

        // caulk+ calculate w5, w6
        let t_poly = UnivariatePolynomial::lagrange(table.clone()).ifft();
        assert_eq!(t_poly.evaluate(&scalar_1), Fr::from(1 as u64));
        assert_eq!(t_poly.evaluate(&v_root_of_unity), Fr::from(2 as u64));
        assert_eq!(t_poly.evaluate(&v_root_of_unity.pow([2 as u64])), Fr::from(3 as u64));
        assert_eq!(t_poly.evaluate(&v_root_of_unity.pow([3 as u64])), Fr::from(4 as u64));
        let t_comm_1 = Pcs::commit_and_write(&pp, &t_poly, &mut transcript).unwrap();
        // z_h_poly = X^t - 1, [-1, 0, ..., 0, 1], t-1 0s in between
        let z_h_poly_coeffs = vec![scalar_1.neg()].into_iter().chain(vec![scalar_0; t - 1]).chain(vec![scalar_1]).collect();
        let z_h_poly = UnivariatePolynomial::monomial(z_h_poly_coeffs);
        let z_h_comm_1 = Pcs::commit_and_write(&pp, &z_h_poly, &mut transcript).unwrap();

        // calculate barycentric_weights
        let bc_weights = barycentric_weights(&h_i);
        let log_t = t.sqrt();
        let t_root_of_unity = root_of_unity::<Fr>(log_t);

        // w5_poly = (t_poly - t_I_poly) / z_I_poly
        let w5_poly_direct = &(&t_poly - &t_i_poly) / &z_i_poly;
        // todo: move this calculation to common inputs
        // q_t_poly_i = (t_poly - table[i])/X-root_of_unity^i
        let q_t_polys: Vec<_> = table.clone().into_iter().enumerate().map(|(i, point)| &(t_poly.clone() + point.neg()) / &(x_poly.clone() + t_root_of_unity.pow([i as u64]).neg())).collect();
        // optimize w5_poly = bc_weights[0] * q_t_polys[i_values[0]] + bc_weights[1] * q_t_polys[i_values[1]] + ... + bc_weights[h_i.len()-1] * q_t_polys[i_values[h_i.len()-1]]
        let w5_poly = bc_weights.clone().into_iter().enumerate().map(|(i, weight)| &q_t_polys[i_values[i]] * weight).reduce(|acc, poly| &acc + poly).unwrap();
        assert_eq!(w5_poly, w5_poly_direct);
        // w6_poly = z_H_poly * (bc_weights[0] * 1/X-h_i[0] + bc_weights[1] * 1/X-h_i[1] + ... + bc_weights[h_i.len()-1] * 1/X-h_i[h_i.len()-1])
        let denom_polys: Vec<_> = h_i.into_iter().map(|root| UnivariatePolynomial::monomial(vec![root.neg(), scalar_1])).collect();
        let w6_poly = bc_weights.clone().into_iter().enumerate().map(|(i, weight)| &(&z_h_poly / &denom_polys[i]) * weight).reduce(|acc, poly| &acc + poly).unwrap();
        // w6_poly = z_H_poly / z_I_poly
        let w6_poly_direct = &z_h_poly / &z_i_poly;
        assert_eq!(w6_poly, w6_poly_direct);

        let w5_comm_1 = Pcs::commit_and_write(&pp, &w5_poly, &mut transcript).unwrap();
        let w6_comm_1 = Pcs::commit_and_write(&pp, &w6_poly, &mut transcript).unwrap();
        // print!("w5_comm_1: {:?}\n", w5_comm_1);
        // print!("w6_comm_1: {:?}\n", w6_comm_1);

        // verification
        let g1 = G1::generator();
        let g2 = G2::generator();
        let g1_affine = G1Affine::from(g1);
        let g2_affine = G2Affine::from(g2);
        let z_i_comm_2_affine = z_i_comm_2.to_affine();
        // 1. verify subtable
        // subtable_lhs = ec_lincomb([
        //     (t_comm_1, scalar_one),
        //     (t_I_comm_1, -scalar_one),
        //     (z_H_comm_1, gamma),
        // ])
        // subtable_rhs = ec_lincomb([
        //     (w5_comm_1, scalar_one),
        //     (w6_comm_1, gamma),
        // ])
        let subtable_msm_lhs = variable_base_msm(
            &[scalar_1, -scalar_1, gamma],
            &[t_comm_1.clone().to_affine(), t_i_comm_1.clone().to_affine(), z_h_comm_1.clone().to_affine()]
        ).into();
        let subtable_msm_rhs = variable_base_msm(
            &[scalar_1, gamma],
            &[w5_comm_1.clone().to_affine(), w6_comm_1.clone().to_affine()]
        ).into();
        let subtable_pairing_lhs = pairing(&subtable_msm_lhs, &g2_affine.clone());
        let subtable_pairing_rhs = pairing(&subtable_msm_rhs, &z_i_comm_2_affine.clone());
        assert_eq!(subtable_pairing_lhs, subtable_pairing_rhs);
        println!("Finished to verify: subtable");

        // 2. verify w1 for X = α
        // # w1 = X^(d-m+1) * (E(X) - e(α) + (φ(X) - a(α))γ) / X - α
        // let coeffs = vec![scalar_0; d - m + 1].into_iter().chain(vec![scalar_1]).collect();
        // let x_exponent_poly = UnivariatePolynomial::monomial(coeffs);
        let x_exponent_poly_comm_2 = Pcs::commit_monomial_g2(&param, &x_exponent_poly.coeffs());
        println!("x_exponent_poly_comm_2: {:?}", x_exponent_poly_comm_2);
        let w1_comm_1_affine: G1Affine = w1_comm_1.to_affine();
        // calculate left hand side pairing
        let w1_lhs = pairing(&w1_comm_1_affine, &vp.s_g2());
        // calculate right hand side pairing
        let w1_rhs1: G1Affine = variable_base_msm(
            &[scalar_1, -v1, gamma, -gamma.mul(v2)],
            &[e_comm_1.clone().to_affine(), g1_affine.clone(), phi_comm_1.clone().to_affine(), g1_affine.clone()]
        ).into();
        let w1_rhs2: G1Affine = variable_base_msm(
            &[alpha],
            &[w1_comm_1_affine.clone()]
        ).into();
        let x_exponent_poly_comm_2_affine = x_exponent_poly_comm_2.to_affine();
        print!("w1_rhs2: {:?}\n", w1_rhs2);
        assert_eq!(vp.g2(), g2_affine);
        let g1_terms = vec![w1_rhs1, w1_rhs2];
        let g2_terms = vec![x_exponent_poly_comm_2_affine, g2_affine];
        let w1_rhs = multi_pairing(&g1_terms, &g2_terms);

        assert_eq!(w1_lhs, w1_rhs);

        println!("Finished to verify: w1");

        // 3. verify w2 for X = 0
        // # X^(d-m+2)
        let coeffs_x_exponent_poly_2 = vec![scalar_0; d - m + 2].into_iter().chain(vec![scalar_1]).collect();
        let x_exponent_poly_2 = UnivariatePolynomial::monomial(coeffs_x_exponent_poly_2);
        let x_exponent_poly_2_comm_1 = Pcs::commit_and_write(&pp, &x_exponent_poly_2, &mut transcript).unwrap();
        let x_exponent_poly_2_comm_2 = Pcs::commit_monomial_g2(&param, &x_exponent_poly_2.coeffs());
        // X^m
        let x_m_exponent_poly_comm_1 = Pcs::commit_and_write(&pp, &x_m_exponent_poly.clone(), &mut transcript).unwrap();
        // to affine
        let x_exponent_poly_2_comm_1_affine = x_exponent_poly_2_comm_1.to_affine();
        let x_exponent_poly_2_comm_2_affine = x_exponent_poly_2_comm_2.to_affine();
        let x_m_exponent_poly_comm_1_affine = x_m_exponent_poly_comm_1.to_affine();
        let r_comm_1_affine = r_comm_1.to_affine();
        let w2_comm_1_affine: G1Affine = w2_comm_1.to_affine();

        // w2_rhs1 = ec_lincomb([
        //     (b.G1, scalar_one),
        //     (x_exp_poly_2_comm_1, gamma ** 2),
        // ])
        let w2_rhs1 = variable_base_msm(
            &[gamma_2, scalar_1],
            &[x_exponent_poly_2_comm_1_affine.clone(), g1_affine.clone()]
        ).into();

        // w2_rhs2 = ec_lincomb([
        //     (R_comm_1, gamma ** 3),
        //     (x_m_exp_poly_comm_1, -gamma ** 2),
        // ])
        let w2_rhs2 = variable_base_msm(
            &[gamma_3, -gamma_2],
            &[r_comm_1_affine.clone(), x_m_exponent_poly_comm_1_affine.clone()],
        ).into();

        // w2_rhs3 = ec_lincomb([
        //     (R_comm_1, gamma),
        //     (b.G1, -v3),
        // ])
        let w2_rhs3 = variable_base_msm(
            &[gamma, -v3],
            &[r_comm_1_affine.clone(), g1_affine.clone()]
        ).into();
        // assert b.pairing(x2, w2_comm_1) == b.pairing(z_I_comm_2, w2_rhs1) * b.pairing(
        //     x_exp_poly_2_comm_2, w2_rhs2) * b.pairing(b.G2, w2_rhs3), "w2 paring check failed"
        // print("Finished to verify: w2")

        let w2_lhs = pairing(&w2_comm_1_affine, &vp.s_g2());
        let w2_pairing_g1_terms = vec![w2_rhs1, w2_rhs2, w2_rhs3];
        let w2_pairing_g2_terms = vec![z_i_comm_2_affine.clone(), x_exponent_poly_2_comm_2_affine.clone(), g2_affine.clone()];
        let rhs = multi_pairing(&w2_pairing_g1_terms, &w2_pairing_g2_terms);
        assert_eq!(w2_lhs, rhs);
        println!("Finished to verify: w2");

        // 4. verify w3 for X = β
        let w3_comm_1_affine: G1Affine = w3_comm_1.to_affine();
        // # calculate commitment [P_D(X)]1
        let p_d_comm_1_affine: G1Affine = variable_base_msm(
            &[v1, -v2, -scalar_1, -v4],
            &[t_i_comm_1.clone().to_affine(), g1_affine.clone(), r_comm_1_affine.clone(), q_d_comm_1.clone().to_affine()]
        ).into();
        let w3_rhs1 = variable_base_msm(
            &[scalar_1, beta, -v1 - gamma.mul(v4), gamma_2],
            &[d_comm_1.clone().to_affine(), w3_comm_1_affine.clone(), g1_affine.clone(), p_d_comm_1_affine.clone()]
        ).into();
        let w3_rhs2 = variable_base_msm(
            &[gamma],
            &[g1_affine.clone()]
        ).into();

        let w3_lhs = pairing(&w3_comm_1_affine, &vp.s_g2());
        let w3_pairing_g1_terms = vec![w3_rhs1, w3_rhs2];
        let w3_pairing_g2_terms = vec![g2_affine.clone(), z_i_comm_2_affine.clone()];
        let w3_rhs = multi_pairing(&w3_pairing_g1_terms, &w3_pairing_g2_terms);
        assert_eq!(w3_lhs, w3_rhs);
        println!("Finished to verify: w3");

        // 5. verify w4 for X = ζ
        let p_e_comm_1_affine: G1Affine = variable_base_msm(
            &[v5.mul(beta), -v5 + v4 * v3.invert().unwrap(), -z_v_zeta],
            &[g1_affine.clone(), v_comm_1.clone().to_affine(), q_e_comm_1.clone().to_affine()]
        ).into();
        let w4_comm_1_affine: G1Affine = w4_comm_1.to_affine();
        let w4_msm_rhs = variable_base_msm(
            &[scalar_1, zeta, gamma, -v5],
            &[e_comm_1.clone().to_affine(), w4_comm_1_affine.clone(), p_e_comm_1_affine.clone(), g1_affine.clone()]
        ).into();
        let w4_pairing_lhs = pairing(&w4_comm_1_affine, &vp.s_g2());
        let w4_pairing_rhs = pairing(&w4_msm_rhs, &g2_affine.clone());
        assert_eq!(w4_pairing_lhs, w4_pairing_rhs);
        println!("Finished to verify: w4");
        vec![Fr::one(), Fr::one()]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2_curves::bn256::Fr;
    use crate::util::transcript::{FieldTranscriptRead, FieldTranscriptWrite, G2TranscriptRead, G2TranscriptWrite};
    type Pcs = UnivariateKzg<Bn256>;

    #[test]
    fn test_baloo_proof() {
        let lookup = vec![Fr::from(3), Fr::from(2), Fr::from(3), Fr::from(4)];
        let table = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let baloo = Baloo::new(table);
        let result = baloo.round1(lookup);
        println!("result: {:?}", result);
    }

    #[test]
    fn test_verify() {

        let lookup = vec![Fr::one(), Fr::one()];
        let table = vec![Fr::one(), Fr::one()];
        let m = lookup.len();

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
        let proof = {
            let mut transcript = Keccak256Transcript::new(());
            let poly = <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::monomial(lookup.clone());
            print!("coeffs: {:?}\n", poly.coeffs());
            let comm = Pcs::commit_and_write(&pp, &poly, &mut transcript).unwrap();
            let point = <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::squeeze_point(m, &mut transcript);
            let eval = poly.evaluate(&point);

            // Use the correct method to write the field element
            transcript.write_field_element(&eval).unwrap();

            Pcs::open(&pp, &poly, &comm, &point, &eval, &mut transcript).unwrap();
            transcript.into_proof()
        };

        // Verify
        let result = {
            let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
            Pcs::verify(
                &vp,
                &Pcs::read_commitment(&vp, &mut transcript).unwrap(),
                &<Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::squeeze_point(m, &mut transcript),

                // Use the correct method to read the field element
                &transcript.read_field_element().unwrap(),
                &mut transcript,
            )
        };
        assert_eq!(result, Ok(()));

    }

    #[test]
    fn test_transcript() {
        let lookup = vec![Fr::one(), Fr::one()];
        let table = vec![Fr::one(), Fr::from(2)];
        let m = lookup.len();
        let poly = UnivariatePolynomial::monomial(lookup.clone());
        let poly_table = UnivariatePolynomial::monomial(table.clone());
        // Setup
        let mut rng = OsRng;
        let poly_size = m;
        print!("poly_size: {:?}\n", poly_size);
        let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
        let (pp, vp) = Pcs::trim(&param, poly_size, 1).unwrap();

        let mut transcript = Keccak256Transcript::new(());
        transcript.write_field_element(&Fr::from(1 as u64)).unwrap();
        transcript.write_field_element(&Fr::from(2 as u64)).unwrap();
        transcript.write_field_element(&Fr::from(3 as u64)).unwrap();

        let comm = <Pcs as PolynomialCommitmentScheme<Fr>>::commit(&pp, &poly).unwrap();
        println!("comm: {:?}", comm);
        println!("comm affine: {:?}", comm.clone().to_affine());
        transcript.write_commitment(&comm.clone().to_affine()).unwrap();
        let comm_1 = Pcs::commit_monomial(&pp, &poly.coeffs());
        println!("comm_1: {:?}", comm_1);
        assert_eq!(comm_1, comm);

        // g2
        let comm_2 = Pcs::commit_monomial_g2(&param, &poly.coeffs());
        println!("comm_2: {:?}", comm_2);
        let comm_2_affine = comm_2.clone().to_affine();
        transcript.write_commitment_g2(&comm_2_affine).unwrap();

        let comm_table = Pcs::commit_and_write(&pp, &poly_table, &mut transcript).unwrap();
        println!("comm_table: {:?}", comm_table);

        let proof = transcript.into_proof();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());

        let a: Fr = transcript.read_field_element().unwrap();
        println!("a: {:?}", a);
        assert_eq!(a, Fr::from(1 as u64));
        let b: Fr = transcript.read_field_element().unwrap();
        println!("b: {:?}", b);
        assert_eq!(b, Fr::from(2 as u64));
        let c: Fr = transcript.read_field_element().unwrap();
        println!("c: {:?}", c);
        assert_eq!(c, Fr::from(3 as u64));

        let comm_back = Pcs::read_commitment(&vp, &mut transcript).unwrap();
        println!("comm_back: {:?}", comm_back);

        // g2
        let comm_table_g2: G2Affine = transcript.read_commitment_g2().unwrap();
        println!("comm_table_g2: {:?}", comm_table_g2);

        let comm_table_back = Pcs::read_commitment(&vp, &mut transcript).unwrap();
        println!("comm_table_back: {:?}", comm_table_back);
        assert_eq!(comm, comm_back);
        assert_eq!(comm_table, comm_table_back);
        assert_eq!(comm_2.to_affine(), comm_table_g2);
    }

}
