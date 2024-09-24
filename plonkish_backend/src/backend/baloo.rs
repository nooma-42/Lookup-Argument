use rand::rngs::OsRng;
use std::{fmt::Debug, marker::PhantomData};

use halo2_curves::bn256::{Bn256, Fr};

use crate::{
    poly::Polynomial,
    poly::univariate::UnivariatePolynomial,
    backend::baloo::preprocessor::preprocess,
    pcs::{
        PolynomialCommitmentScheme,
        univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam},
    },
    util::{
        arithmetic::{Field, PrimeField, MultiMillerLoop, root_of_unity},
        test::std_rng,
        Deserialize, DeserializeOwned, Itertools, Serialize,
        transcript::{InMemoryTranscript, TranscriptRead, TranscriptWrite, Keccak256Transcript},
    }
};


pub mod preprocessor;
pub mod prover;
pub mod verifier;

#[derive(Clone, Debug)]
pub struct BalooProverParam //<F, Pcs>
// where
//     F: PrimeField,
//     Pcs: PolynomialCommitmentScheme<F>,
{
    pub(crate) num_vars: usize,
}

#[derive(Clone, Debug)]
pub struct BalooVerifierParam<F, Pcs>
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F>,
{
    // [z_H_comm_1, t_comm_1]
    pub(crate) preprocess_comms: Vec<Pcs::Commitment>,
}
pub struct Baloo<F> {
    table: Vec<F>,
    // round1: (&self, &UnivariateKzgProverParam<M>, Vec<F>, Vec<F>) -> Vec<F>
}

impl<F: Field> Baloo<F>
{
    pub fn new(table: Vec<F>) -> Self {
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
    fn round1<Pcs, T>(&self, lookup: Vec<F>) -> Vec<F>
    where
        Pcs: PolynomialCommitmentScheme<F, Polynomial = UnivariatePolynomial<F>>,
        T: TranscriptRead<Pcs::CommitmentChunk, F>
            + TranscriptWrite<Pcs::CommitmentChunk, F>
            + InMemoryTranscript<Param = ()>,
    {
        // type ProverParam = UnivariateKzgProverParam<M>;
        let m = lookup.len();
        let mut rng = std_rng();
        // let phi_values = lookup.iter().map(|&x| PrimeField::<F>::from(x)).collect::<Vec<F>>();
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
            let mut transcript = T::new(());
            let poly = Pcs::Polynomial::monomial(
                lookup.clone(),
            );
            print!("coeffs: {:?}\n", poly.coeffs());
            let comm = Pcs::commit_and_write(&pp, &poly, &mut transcript).unwrap();
            let point = Pcs::Polynomial::squeeze_point(m, &mut transcript);
            let eval = poly.evaluate(&point);
            transcript.write_field_element(&eval).unwrap();
            Pcs::open(&pp, &poly, &comm, &point, &eval, &mut transcript).unwrap();
            transcript.into_proof()
        };
        // Verify
        let result = {
            let mut transcript = T::from_proof((), proof.as_slice());
            Pcs::verify(
                &vp,
                &Pcs::read_commitment(&vp, &mut transcript).unwrap(),
                &<Pcs as PolynomialCommitmentScheme<F>>::Polynomial::squeeze_point(m, &mut transcript),
                &transcript.read_field_element().unwrap(),
                &mut transcript,
            )
        };
        assert_eq!(result, Ok(()));

        // let phi_poly = Pcs::new(lookup);
        // commit phi(X) on G1
        // let phi_comm_1 = Self::setup.commit_g1(&pp, &phi_poly);
        // let phi_comm_1 = Pcs::commit(&pp, &phi_poly);
        // remove duplicated elements
        // let t_values_from_lookup: HashSet<_> = lookup.into_iter().collect();
        // // I: the index of t_values_from_lookup elements in sub table t_I
        // let i_values: Vec<_> = t_values_from_lookup.iter().map(|elem| table.iter().position(|&x| x == *elem).unwrap()).collect();
        // let roots_of_unity = Vec::new(); // TODO: get roots of unity
        // // H_I = {ξ_i} , i = [1, k], ξ(Xi)
        // let h_i = i_values.iter().map(|&i| roots_of_unity[i]).collect();
        // let h_i_interp_poly = InterpolationPoly::new(&h_i, &t_values_from_lookup); // TODO: implement interpolation polynomial
        // // Multiple all root polynomial: (X - I_0)(X - I_1)(X - I_2)...
        // let mut col_values = Vec::new();
        // let mut v_values = Vec::new();
        // for i in 0..m {
        //     // find the index of 1 in jth row of M
        //     let col_i = t_values_from_lookup.iter().position(|&x| x == lookup[i]).unwrap();
        //     col_values.push(col_i);
        //     let col_i_root = h_i[col_i];
        //     // Note: v = 1 / col_i_root in paper
        //     // Here we use different construction that does not affect the verification
        //     let v = col_i_root;
        //     v_values.push(v);
        // }
        // // Vanishing polynomial in coefficient form
        // Self::z_i_poly = h_i_interp_poly.vanishing_poly();
        // assert_eq!(col_values.iter().map(|&elem| t_values_from_lookup[elem]).collect::<Vec<_>>(), lookup);
        // // ξ(x) polynomial
        // let v_poly = Polynomial::new(v_values, Basis::Lagrange);
        // // Refer to section 5. cannot use FFT due to it's not in multiplicative subgroup
        // let t_i_poly = InterpolationPoly::new(&h_i, &t_values_from_lookup);
        // Self::t_i_poly = t_i_poly.poly();

        // Self::col = col_values.clone();
        // Self::h_i = h_i;
        // Self::phi_poly = phi_poly.ifft();
        // Self::phi_comm_1 = Self::setup.commit_g1(&Self::phi_poly);
        // Self::v_poly = v_poly.ifft();

        // // Commit
        // // π1 = ([z_i]_2 = [z_i(x)]_2, [v]_1 = [v(x)]_1, t = [t(x)]_1)
        // Self::z_i_comm_2 = Self::setup.commit_g2(&Self::z_i_poly);
        // Self::v_comm_1 = Self::setup.commit_g1(&Self::v_poly);
        // Self::t_i_comm_1 = Self::setup.commit_g1(&Self::t_i_poly);

        // Message1::new(Self::z_i_comm_2.clone(), Self::v_comm_1.clone(), Self::t_i_comm_1.clone(), Self::phi_comm_1.clone())
        Vec::new()
    }

    // type Pcs = Pcs;
    // type ProverParam = BalooProverParam;
    // type VerifierParam = BalooVerifierParam<F, Pcs>;

}

fn lagrange_interp(h_i_values: &[Fr], t_values_from_lookup: &[Fr]) -> UnivariatePolynomial<Fr> {
    assert!(h_i_values.len() == t_values_from_lookup.len());
    println!("h_i_values: {:?}", h_i_values);

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::{collections::HashSet, ops::{Add, Mul}};
    use bitvec::vec;
    use halo2_curves::bn256::Fr;
    use num_integer::Roots;
    use crate::util::transcript::{FieldTranscriptRead, FieldTranscriptWrite};

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
    fn test_e2e() {
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
    fn test_prover() {
        let lookup = vec![Fr::from(3), Fr::from(2), Fr::from(3), Fr::from(4)];
        let table = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        // let m = 16;
        let m = lookup.len();

        // Setup
        let mut rng = OsRng;
        let poly_size = m;
        print!("poly_size: {:?}\n", poly_size);
        // TODO: large srs size
        let srs_size = 1 << 8;
        print!("srs_size: {:?}\n", srs_size);
        let param = Pcs::setup(srs_size, 1, &mut rng).unwrap();
        let (pp, vp) = Pcs::trim(&param, srs_size, 1).unwrap();
        print!("lookup: {:?}\n", lookup);

        // Commit and open
        let mut transcript = Keccak256Transcript::new(());
        // commit phi(X) on G1
        let phi_poly = UnivariatePolynomial::lagrange(lookup.clone()).ifft();
        print!("coeffs: {:?}\n", phi_poly.coeffs());
        let phi_comm_1 = Pcs::commit_and_write(&pp, &phi_poly, &mut transcript).unwrap();

        // remove duplicated elements
        let t_values_from_lookup: HashSet<_> = lookup.clone().into_iter().collect();
        print!("t_values_from_lookup: {:?}\n", t_values_from_lookup);
        // I: the index of t_values_from_lookup elements in sub table t_I
        let i_values: Vec<_> = t_values_from_lookup.iter().map(|elem| table.iter().position(|&x| x == *elem).unwrap()).collect();
        print!("i_values: {:?}\n", i_values);
        // todo
        let log_m = m.sqrt();
        let v_root_of_unity = root_of_unity::<Fr>(log_m);
        println!("v_root_of_unity: {:?}", v_root_of_unity);
        assert_eq!(v_root_of_unity.pow([m as u64]), Fr::one());
        for i in 1..m+1 {
            println!("i: {:?}, v_root_of_unity^i: {:?}", i, v_root_of_unity.pow([i as u64]));
        }

        // H_I = {ξ_i} , i = [1...k], ξ(Xi)
        let h_i: Vec<_> = i_values.iter().map(|&i| {
            let i_as_u64 = i as u64;
            println!("i_as_u64: {:?}", i_as_u64);
            v_root_of_unity.pow([i_as_u64])
        }).collect();
        print!("h_i: {:?}\n", h_i);
        // TODO: optimize interpolation polynomial with https://github.com/gy001/hypercube/blob/main/univarization/src/unipoly.rs#L391
        // refer to barycentric_weights in arithmetic.rs
        let t_values_vec: Vec<Fr> = t_values_from_lookup.iter().cloned().collect();
        let t_i_poly = lagrange_interp(&h_i, &t_values_vec);
        let t_commit_1 = Pcs::commit_and_write(&pp, &t_i_poly, &mut transcript);

        let z_i_poly = UnivariatePolynomial::vanishing(&h_i, Fr::one());
        let z_i_comm_2 = Pcs::commit_monomial_g2(&param, &z_i_poly.coeffs());
        // TODO
        // transcript.write_commitments(&z_i_comm_2.0.x.c0).unwrap();
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
        // TODO: get randomness from transcript
        let alpha = Fr::from(2 as u64);
        let beta = Fr::from(3 as u64);
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
        // let t_commit_1 = Pcs::commit_and_write(&pp, &t_i_poly, &mut transcript);
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
        let alpha = Fr::from(2 as u64);
        let beta = Fr::from(3 as u64);
        let gamma = Fr::from(4 as u64);
        let zeta = Fr::from(5 as u64);

        let d = srs_size - 2;
        // calculate v1, v2, v3, v4, v5
        // v1 = e(α)
        let v1 = e_poly.evaluate(&alpha);
        // v2 = a(α)
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
        // todo: use fft to do division?
        let x_alpha_poly = UnivariatePolynomial::monomial(vec![-alpha, scalar_1]);
        // calculate w1 = (E(X) - e(α) + (φ(X) - a(α))γ) / X - α
        let w1 = &(&(e_poly.clone() + v1.neg()) + &(phi_poly.clone() + v2.neg()) * gamma) / &x_alpha_poly;
        // calculate polynomial X
        let x_poly = UnivariatePolynomial::monomial(vec![scalar_0, scalar_1]);
        // X^m
        let x_m_exponent_poly = UnivariatePolynomial::monomial(vec![scalar_0; m].into_iter().chain(vec![scalar_1]).collect());
        // calculate w2 = (z_I(X) - z_I(0) / X + γ * R(X) / X +  γ^2 * X^(d-m+1) * (z_I(X) - X^m)) + γ^3 * X^(d-m+1) * R(X)
        // todo: use fft to do division?
        let w2 = (
                &(
                    &(&(z_i_poly.clone() + v3.neg()) / &x_poly)
                    + &(&(&r_poly * gamma) / &x_poly)
                )
                + &x_exponent_poly
            ).poly_mul(
                &(&(&z_i_poly - x_m_exponent_poly) * gamma.mul(gamma))
                + &r_poly * gamma.mul(gamma).mul(gamma)
            );
        print!("w2: {:?}\n", w2.degree());
        // let w2 = (z_I_poly - v3) / x_poly + R_poly * gamma / x_poly + x_exponent_poly * ( (z_I_poly - x_m_exponent_poly) * gamma ** 2 + R_poly * gamma ** 3 )
        // calculate X - β
        let x_beta_poly = UnivariatePolynomial::monomial(vec![-beta, scalar_1]);
        // calculate w3 = (D(X) - E(α) + (z_I(X) - z_I(β))γ + P_D(X)γ^2) / X - β
        // v1 = E(α) == D(β)
        let w3 = &(&(&(d_poly + v1.neg()) + &(z_i_poly + v4.neg()) * gamma) + &p_d_poly * gamma.mul(gamma)) / &x_beta_poly;
        // calculate X - ζ
        let x_zeta_poly = UnivariatePolynomial::monomial(vec![-zeta, scalar_1]);
        // calculate w4 = (E(X) - E(ζ) + P_E(X)γ) / X - ζ
        // v5 = E(ζ)
        // todo: use fft to do division?
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

        // todo: caulk+ calculate w5, w6
        // w5_poly = (t_poly - t_I_poly) / z_I_poly
        // w6_poly = z_H_poly / z_I_poly
        // let w5_poly = &(&t_poly - &t_i_poly) / &z_i_poly;
        // let w6_poly = &z_h_poly / &z_i_poly;

        // let w5_comm_1 = Pcs::commit_and_write(&pp, &w5_poly, &mut transcript).unwrap();
        // let w6_comm_1 = Pcs::commit_and_write(&pp, &w6_poly, &mut transcript).unwrap();

        // print!("w5_comm_1: {:?}\n", w5_comm_1);
        // print!("w6_comm_1: {:?}\n", w6_comm_1);
    }
}