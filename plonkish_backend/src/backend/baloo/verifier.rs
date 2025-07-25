use crate::{
    backend::baloo::{
        key::BalooVerifierKey,
        util::{log_2, multi_pairing},
    },
    pcs::{
        univariate::{
            UnivariateKzg, UnivariateKzgCommitment,
            UnivariateKzgVerifierParam,
        }, PolynomialCommitmentScheme,
    },
    poly::{univariate::UnivariatePolynomial, Polynomial},
    util::{
        arithmetic::{root_of_unity, variable_base_msm, Field},
        transcript::{
            FieldTranscript, FieldTranscriptRead, G2TranscriptRead, InMemoryTranscript, Keccak256Transcript,
        },
    },
};
use halo2_curves::bn256::{pairing, Bn256, Fr, G1Affine, G2Affine, G1, G2};
use std::ops::Mul;

type Pcs = UnivariateKzg<Bn256>;

pub struct Verifier<'b> {
    vk: &'b BalooVerifierKey,
}

impl Verifier<'_> {
    pub fn new(vk: &BalooVerifierKey) -> Verifier<'_> {
        Verifier { vk }
    }

    pub fn verify(
        &self,
        proof: &Vec<u8>,
        phi_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_m_exponent_poly_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_comm_2: &UnivariateKzgCommitment<G2Affine>,
        x_exponent_poly_2_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_2_comm_2: &UnivariateKzgCommitment<G2Affine>,
        m: usize,
    ) -> bool {
        let scalar_0 = Fr::from(0_u64);
        let scalar_1 = Fr::from(1_u64);
        let vp = &self.vk.vp;
        
        // Extract preprocessed commitments from verifier key
        let t_comm_1 = &self.vk.table_comm;
        let z_h_comm_1 = &self.vk.z_h_comm;
        
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());

        // read pi_1 = (v_comm_1.clone(), z_i_comm_2.clone(), t_i_comm_1.clone());
        let v_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();
        // g2
        let z_i_comm_2: G2Affine = transcript.read_commitment_g2().unwrap();

        let t_i_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let alpha: Fr = transcript.squeeze_challenge();

        // read pi_2 = (d_comm_1.clone(), r_comm_1.clone(), q_d_comm_1.clone());
        let d_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let r_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let q_d_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let beta: Fr = transcript.squeeze_challenge();

        // read pi_3 = (e_comm_1.clone(), q_e_comm_1.clone());
        let e_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let q_e_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let gamma: Fr = transcript.squeeze_challenge();
        let zeta: Fr = transcript.squeeze_challenge();
        let gamma_2 = gamma.mul(gamma);
        let gamma_3 = gamma_2.mul(gamma);

        // read pi_4 = (v1, v2, v3, v4, v5, a_comm_1.clone(), w1_comm_1.clone(), w2_comm_1.clone(), w3_comm_1.clone(), w4_comm_1.clone());
        let v1: Fr = transcript.read_field_element().unwrap();

        let v2: Fr = transcript.read_field_element().unwrap();

        let v3: Fr = transcript.read_field_element().unwrap();

        let v4: Fr = transcript.read_field_element().unwrap();

        let v5: Fr = transcript.read_field_element().unwrap();

        let a_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let w1_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let w2_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let w3_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let w4_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        // Construct X^m - 1, [-1, 0, 0, ..., 1], m - 1 0s in between
        let z_v_values: Vec<Fr> = vec![scalar_1.neg()]
            .into_iter()
            .chain((0..m - 1).map(|_| scalar_0))
            .chain(vec![scalar_1])
            .collect();
        assert_eq!(z_v_values.len(), m + 1);

        let log_m = log_2(m);
        let v_root_of_unity = root_of_unity::<Fr>(log_m);
        // X^m - 1
        let z_v_poly = UnivariatePolynomial::monomial(z_v_values);
        assert_eq!(z_v_poly.evaluate(&v_root_of_unity), scalar_0);
        let z_v_zeta = z_v_poly.evaluate(&zeta);

        /************
        Verification
        ************/
        let g1 = G1::generator();
        let g2 = G2::generator();
        let g1_affine = G1Affine::from(g1);
        let g2_affine = G2Affine::from(g2);
        let z_i_comm_2_affine = z_i_comm_2;
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
            &[
                t_comm_1.clone().to_affine(),
                t_i_comm_1.clone().to_affine(),
                z_h_comm_1.clone().to_affine(),
            ],
        )
        .into();
        let subtable_msm_rhs = a_comm_1.clone().to_affine();
        let subtable_pairing_lhs = pairing(&subtable_msm_lhs, &g2_affine.clone());
        let subtable_pairing_rhs = pairing(&subtable_msm_rhs, &z_i_comm_2_affine.clone());
        assert_eq!(subtable_pairing_lhs, subtable_pairing_rhs);
        //println!("Finished to verify: subtable");

        // 2. verify w1 for X = α
        // # w1 = X^(d-m+1) * (E(X) - e(α) + (φ(X) - a(α))γ) / X - α
        let w1_comm_1_affine: G1Affine = w1_comm_1.to_affine();
        // calculate left hand side pairing
        let w1_lhs = pairing(&w1_comm_1_affine, &vp.s_g2());
        // calculate right hand side pairing
        let w1_rhs1: G1Affine = variable_base_msm(
            &[scalar_1, -v1, gamma, -gamma.mul(v2)],
            &[
                e_comm_1.clone().to_affine(),
                g1_affine,
                phi_comm_1.clone().to_affine(),
                g1_affine,
            ],
        )
        .into();
        let w1_rhs2: G1Affine = variable_base_msm(&[alpha], &[w1_comm_1_affine]).into();
        let x_exponent_poly_comm_2_affine = x_exponent_poly_comm_2.clone().to_affine();
        assert_eq!(vp.g2(), g2_affine);
        let g1_terms = vec![w1_rhs1, w1_rhs2];
        let g2_terms = vec![x_exponent_poly_comm_2_affine, g2_affine];
        let w1_rhs = multi_pairing(&g1_terms, &g2_terms);

        assert_eq!(w1_lhs, w1_rhs);

        //println!("Finished to verify: w1");

        // 3. verify w2 for X = 0
        // to affine
        let x_exponent_poly_2_comm_1_affine = x_exponent_poly_2_comm_1.clone().to_affine();
        let x_exponent_poly_2_comm_2_affine = x_exponent_poly_2_comm_2.clone().to_affine();
        let x_m_exponent_poly_comm_1_affine = x_m_exponent_poly_comm_1.clone().to_affine();
        let r_comm_1_affine = r_comm_1.to_affine();
        let w2_comm_1_affine: G1Affine = w2_comm_1.to_affine();

        // w2_rhs1 = ec_lincomb([
        //     (b.G1, scalar_one),
        //     (x_exp_poly_2_comm_1, gamma ** 2),
        // ])
        let w2_rhs1 = variable_base_msm(
            &[gamma_2, scalar_1],
            &[x_exponent_poly_2_comm_1_affine, g1_affine],
        )
        .into();

        // w2_rhs2 = ec_lincomb([
        //     (R_comm_1, gamma ** 3),
        //     (x_m_exp_poly_comm_1, -gamma ** 2),
        // ])
        let w2_rhs2 = variable_base_msm(
            &[gamma_3, -gamma_2],
            &[
                r_comm_1_affine,
                x_m_exponent_poly_comm_1_affine,
            ],
        )
        .into();

        // w2_rhs3 = ec_lincomb([
        //     (R_comm_1, gamma),
        //     (b.G1, -v3),
        // ])
        let w2_rhs3 =
            variable_base_msm(&[gamma, -v3], &[r_comm_1_affine, g1_affine]).into();
        // assert b.pairing(x2, w2_comm_1) == b.pairing(z_I_comm_2, w2_rhs1) * b.pairing(
        //     x_exp_poly_2_comm_2, w2_rhs2) * b.pairing(b.G2, w2_rhs3), "w2 paring check failed"
        // print("Finished to verify: w2")

        let w2_lhs = pairing(&w2_comm_1_affine, &vp.s_g2());
        let w2_pairing_g1_terms = vec![w2_rhs1, w2_rhs2, w2_rhs3];
        let w2_pairing_g2_terms = vec![
            z_i_comm_2_affine,
            x_exponent_poly_2_comm_2_affine,
            g2_affine,
        ];
        let rhs = multi_pairing(&w2_pairing_g1_terms, &w2_pairing_g2_terms);
        assert_eq!(w2_lhs, rhs);
        //println!("Finished to verify: w2");

        // 4. verify w3 for X = β
        let w3_comm_1_affine: G1Affine = w3_comm_1.to_affine();
        // # calculate commitment [P_D(X)]1
        let p_d_comm_1_affine: G1Affine = variable_base_msm(
            &[v1, -v2, -scalar_1, -v4],
            &[
                t_i_comm_1.clone().to_affine(),
                g1_affine,
                r_comm_1_affine,
                q_d_comm_1.clone().to_affine(),
            ],
        )
        .into();
        let w3_rhs1 = variable_base_msm(
            &[scalar_1, beta, -v1 - gamma.mul(v4), gamma_2],
            &[
                d_comm_1.clone().to_affine(),
                w3_comm_1_affine,
                g1_affine,
                p_d_comm_1_affine,
            ],
        )
        .into();
        let w3_rhs2 = variable_base_msm(&[gamma], &[g1_affine]).into();

        let w3_lhs = pairing(&w3_comm_1_affine, &vp.s_g2());
        let w3_pairing_g1_terms = vec![w3_rhs1, w3_rhs2];
        let w3_pairing_g2_terms = vec![g2_affine, z_i_comm_2_affine];
        let w3_rhs = multi_pairing(&w3_pairing_g1_terms, &w3_pairing_g2_terms);
        assert_eq!(w3_lhs, w3_rhs);
        //println!("Finished to verify: w3");

        // 5. verify w4 for X = ζ
        let p_e_comm_1_affine: G1Affine = variable_base_msm(
            &[v5.mul(beta), -v5 + v4 * v3.invert().unwrap(), -z_v_zeta],
            &[
                g1_affine,
                v_comm_1.clone().to_affine(),
                q_e_comm_1.clone().to_affine(),
            ],
        )
        .into();
        let w4_comm_1_affine: G1Affine = w4_comm_1.to_affine();
        let w4_msm_rhs = variable_base_msm(
            &[scalar_1, zeta, gamma, -v5],
            &[
                e_comm_1.clone().to_affine(),
                w4_comm_1_affine,
                p_e_comm_1_affine,
                g1_affine,
            ],
        )
        .into();
        let w4_pairing_lhs = pairing(&w4_comm_1_affine, &vp.s_g2());
        let w4_pairing_rhs = pairing(&w4_msm_rhs, &g2_affine.clone());
        assert_eq!(w4_pairing_lhs, w4_pairing_rhs);
        //println!("Finished to verify: w4");

        true
    }
}

// Legacy verifier for backward compatibility
pub struct LegacyVerifier<'b> {
    vp: &'b UnivariateKzgVerifierParam<Bn256>,
}

impl LegacyVerifier<'_> {
    pub fn new(vp: &UnivariateKzgVerifierParam<Bn256>) -> LegacyVerifier<'_> {
        LegacyVerifier { vp }
    }

    pub fn verify(
        &self,
        proof: &Vec<u8>,
        t_comm_1: &UnivariateKzgCommitment<G1Affine>,
        z_h_comm_1: &UnivariateKzgCommitment<G1Affine>,
        phi_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_m_exponent_poly_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_comm_2: &UnivariateKzgCommitment<G2Affine>,
        x_exponent_poly_2_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_2_comm_2: &UnivariateKzgCommitment<G2Affine>,
        m: usize,
    ) -> bool {
        let scalar_0 = Fr::from(0_u64);
        let scalar_1 = Fr::from(1_u64);
        let vp = self.vp;
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());

        // read pi_1 = (v_comm_1.clone(), z_i_comm_2.clone(), t_i_comm_1.clone());
        let v_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();
        // g2
        let z_i_comm_2: G2Affine = transcript.read_commitment_g2().unwrap();

        let t_i_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let alpha: Fr = transcript.squeeze_challenge();

        // read pi_2 = (d_comm_1.clone(), r_comm_1.clone(), q_d_comm_1.clone());
        let d_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let r_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let q_d_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let beta: Fr = transcript.squeeze_challenge();

        // read pi_3 = (e_comm_1.clone(), q_e_comm_1.clone());
        let e_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let q_e_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let gamma: Fr = transcript.squeeze_challenge();
        let zeta: Fr = transcript.squeeze_challenge();
        let gamma_2 = gamma.mul(gamma);
        let gamma_3 = gamma_2.mul(gamma);

        // read pi_4 = (v1, v2, v3, v4, v5, a_comm_1.clone(), w1_comm_1.clone(), w2_comm_1.clone(), w3_comm_1.clone(), w4_comm_1.clone());
        let v1: Fr = transcript.read_field_element().unwrap();

        let v2: Fr = transcript.read_field_element().unwrap();

        let v3: Fr = transcript.read_field_element().unwrap();

        let v4: Fr = transcript.read_field_element().unwrap();

        let v5: Fr = transcript.read_field_element().unwrap();

        let a_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let w1_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let w2_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let w3_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        let w4_comm_1 = Pcs::read_commitment(vp, &mut transcript).unwrap();

        // Construct X^m - 1, [-1, 0, 0, ..., 1], m - 1 0s in between
        let z_v_values: Vec<Fr> = vec![scalar_1.neg()]
            .into_iter()
            .chain((0..m - 1).map(|_| scalar_0))
            .chain(vec![scalar_1])
            .collect();
        assert_eq!(z_v_values.len(), m + 1);

        let log_m = log_2(m);
        let v_root_of_unity = root_of_unity::<Fr>(log_m);
        // X^m - 1
        let z_v_poly = UnivariatePolynomial::monomial(z_v_values);
        assert_eq!(z_v_poly.evaluate(&v_root_of_unity), scalar_0);
        let z_v_zeta = z_v_poly.evaluate(&zeta);

        /************
        Verification
        ************/
        let g1 = G1::generator();
        let g2 = G2::generator();
        let g1_affine = G1Affine::from(g1);
        let g2_affine = G2Affine::from(g2);
        let z_i_comm_2_affine = z_i_comm_2;
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
            &[
                t_comm_1.clone().to_affine(),
                t_i_comm_1.clone().to_affine(),
                z_h_comm_1.clone().to_affine(),
            ],
        )
        .into();
        let subtable_msm_rhs = a_comm_1.clone().to_affine();
        let subtable_pairing_lhs = pairing(&subtable_msm_lhs, &g2_affine.clone());
        let subtable_pairing_rhs = pairing(&subtable_msm_rhs, &z_i_comm_2_affine.clone());
        assert_eq!(subtable_pairing_lhs, subtable_pairing_rhs);
        //println!("Finished to verify: subtable");

        // 2. verify w1 for X = α
        // # w1 = X^(d-m+1) * (E(X) - e(α) + (φ(X) - a(α))γ) / X - α
        let w1_comm_1_affine: G1Affine = w1_comm_1.to_affine();
        // calculate left hand side pairing
        let w1_lhs = pairing(&w1_comm_1_affine, &vp.s_g2());
        // calculate right hand side pairing
        let w1_rhs1: G1Affine = variable_base_msm(
            &[scalar_1, -v1, gamma, -gamma.mul(v2)],
            &[
                e_comm_1.clone().to_affine(),
                g1_affine,
                phi_comm_1.clone().to_affine(),
                g1_affine,
            ],
        )
        .into();
        let w1_rhs2: G1Affine = variable_base_msm(&[alpha], &[w1_comm_1_affine]).into();
        let x_exponent_poly_comm_2_affine = x_exponent_poly_comm_2.clone().to_affine();
        assert_eq!(vp.g2(), g2_affine);
        let g1_terms = vec![w1_rhs1, w1_rhs2];
        let g2_terms = vec![x_exponent_poly_comm_2_affine, g2_affine];
        let w1_rhs = multi_pairing(&g1_terms, &g2_terms);

        assert_eq!(w1_lhs, w1_rhs);

        //println!("Finished to verify: w1");

        // 3. verify w2 for X = 0
        // to affine
        let x_exponent_poly_2_comm_1_affine = x_exponent_poly_2_comm_1.clone().to_affine();
        let x_exponent_poly_2_comm_2_affine = x_exponent_poly_2_comm_2.clone().to_affine();
        let x_m_exponent_poly_comm_1_affine = x_m_exponent_poly_comm_1.clone().to_affine();
        let r_comm_1_affine = r_comm_1.to_affine();
        let w2_comm_1_affine: G1Affine = w2_comm_1.to_affine();

        // w2_rhs1 = ec_lincomb([
        //     (b.G1, scalar_one),
        //     (x_exp_poly_2_comm_1, gamma ** 2),
        // ])
        let w2_rhs1 = variable_base_msm(
            &[gamma_2, scalar_1],
            &[x_exponent_poly_2_comm_1_affine, g1_affine],
        )
        .into();

        // w2_rhs2 = ec_lincomb([
        //     (R_comm_1, gamma ** 3),
        //     (x_m_exp_poly_comm_1, -gamma ** 2),
        // ])
        let w2_rhs2 = variable_base_msm(
            &[gamma_3, -gamma_2],
            &[
                r_comm_1_affine,
                x_m_exponent_poly_comm_1_affine,
            ],
        )
        .into();

        // w2_rhs3 = ec_lincomb([
        //     (R_comm_1, gamma),
        //     (b.G1, -v3),
        // ])
        let w2_rhs3 =
            variable_base_msm(&[gamma, -v3], &[r_comm_1_affine, g1_affine]).into();
        // assert b.pairing(x2, w2_comm_1) == b.pairing(z_I_comm_2, w2_rhs1) * b.pairing(
        //     x_exp_poly_2_comm_2, w2_rhs2) * b.pairing(b.G2, w2_rhs3), "w2 paring check failed"
        // print("Finished to verify: w2")

        let w2_lhs = pairing(&w2_comm_1_affine, &vp.s_g2());
        let w2_pairing_g1_terms = vec![w2_rhs1, w2_rhs2, w2_rhs3];
        let w2_pairing_g2_terms = vec![
            z_i_comm_2_affine,
            x_exponent_poly_2_comm_2_affine,
            g2_affine,
        ];
        let rhs = multi_pairing(&w2_pairing_g1_terms, &w2_pairing_g2_terms);
        assert_eq!(w2_lhs, rhs);
        //println!("Finished to verify: w2");

        // 4. verify w3 for X = β
        let w3_comm_1_affine: G1Affine = w3_comm_1.to_affine();
        // # calculate commitment [P_D(X)]1
        let p_d_comm_1_affine: G1Affine = variable_base_msm(
            &[v1, -v2, -scalar_1, -v4],
            &[
                t_i_comm_1.clone().to_affine(),
                g1_affine,
                r_comm_1_affine,
                q_d_comm_1.clone().to_affine(),
            ],
        )
        .into();
        let w3_rhs1 = variable_base_msm(
            &[scalar_1, beta, -v1 - gamma.mul(v4), gamma_2],
            &[
                d_comm_1.clone().to_affine(),
                w3_comm_1_affine,
                g1_affine,
                p_d_comm_1_affine,
            ],
        )
        .into();
        let w3_rhs2 = variable_base_msm(&[gamma], &[g1_affine]).into();

        let w3_lhs = pairing(&w3_comm_1_affine, &vp.s_g2());
        let w3_pairing_g1_terms = vec![w3_rhs1, w3_rhs2];
        let w3_pairing_g2_terms = vec![g2_affine, z_i_comm_2_affine];
        let w3_rhs = multi_pairing(&w3_pairing_g1_terms, &w3_pairing_g2_terms);
        assert_eq!(w3_lhs, w3_rhs);
        //println!("Finished to verify: w3");

        // 5. verify w4 for X = ζ
        let p_e_comm_1_affine: G1Affine = variable_base_msm(
            &[v5.mul(beta), -v5 + v4 * v3.invert().unwrap(), -z_v_zeta],
            &[
                g1_affine,
                v_comm_1.clone().to_affine(),
                q_e_comm_1.clone().to_affine(),
            ],
        )
        .into();
        let w4_comm_1_affine: G1Affine = w4_comm_1.to_affine();
        let w4_msm_rhs = variable_base_msm(
            &[scalar_1, zeta, gamma, -v5],
            &[
                e_comm_1.clone().to_affine(),
                w4_comm_1_affine,
                p_e_comm_1_affine,
                g1_affine,
            ],
        )
        .into();
        let w4_pairing_lhs = pairing(&w4_comm_1_affine, &vp.s_g2());
        let w4_pairing_rhs = pairing(&w4_msm_rhs, &g2_affine.clone());
        assert_eq!(w4_pairing_lhs, w4_pairing_rhs);
        //println!("Finished to verify: w4");

        true
    }
}
