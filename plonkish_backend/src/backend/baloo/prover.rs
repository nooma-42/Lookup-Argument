use crate::{
    backend::baloo::util::{log_2, pow_2},
    backend::baloo::preprocessor::preprocess,
    pcs::{
        univariate::{
            UnivariateKzg, UnivariateKzgCommitment, UnivariateKzgParam, UnivariateKzgProverParam,
        }, PolynomialCommitmentScheme,
    },
    poly::{univariate::UnivariatePolynomial, Polynomial},
    util::{
        arithmetic::{barycentric_weights, root_of_unity, Field},
        transcript::{
            FieldTranscript, FieldTranscriptRead, FieldTranscriptWrite,
            G2TranscriptRead, G2TranscriptWrite, InMemoryTranscript, Keccak256Transcript,
            TranscriptWrite,
        },
    },
};
use halo2_curves::bn256::{Bn256, Fr, G1Affine, G2Affine};
use std::{collections::HashSet, ops::Mul};
use rand::rngs::OsRng;

type Pcs = UnivariateKzg<Bn256>;
type Scalar = Fr;

pub struct Prover<'b> {
    table: &'b Vec<Fr>,
    param: &'b UnivariateKzgParam<Bn256>,
    pp: &'b UnivariateKzgProverParam<Bn256>,
    d: usize,
}

impl Prover<'_> {
    pub fn new<'a>(
        table: &'a Vec<Fr>,
        param: &'a UnivariateKzgParam<Bn256>,
        pp: &'a UnivariateKzgProverParam<Bn256>,
    ) -> Prover<'a> {
        let d = (1 << pp.k()) - 2;
        Prover {
            table,
            param,
            pp,
            d,
        }
    }

    /// Compute polynomial multilication naively in O(n^2)
    fn naive_multiplication(coeffs0: &[Scalar], coeffs1: &[Scalar]) -> Vec<Scalar> {
        let mut c = vec![Scalar::zero(); coeffs0.len() + coeffs1.len() - 1];

        for i in 0..coeffs0.len() {
            for j in 0..coeffs1.len() {
                c[i + j] += coeffs0[i] * coeffs1[j];
            }
        }
        c
    }

    /// Compute polynomial multilication wrapper
    pub fn multiplication(coeffs0: &[Scalar], coeffs1: &[Scalar]) -> Vec<Scalar> {
        // TODO: to select different multiplication algs based
        // on the degree of polynomials.
        Self::naive_multiplication(coeffs0, coeffs1)
    }

    fn eval_rec(
        tree: &Vec<Vec<Vec<Scalar>>>,
        k: usize,
        base: usize,
        f: &[Scalar],
        u: &[Scalar],
    ) -> Vec<Scalar> {
        let n = u.len();
        // println!("eval_rec> k={}, base={}, n={}", k, base, n);
        // println!("f={}", scalar_vector_to_string(&f.to_vec()));
        // println!("u={}", scalar_vector_to_string(&u.to_vec()));

        if k == 0 {
            return f.to_vec();
        }

        // println!("division> k-1={}, left={}, right={}, n={}", k-1, base+0, base+1, n);
        let (q0, r0) = Self::division(f, &tree[k - 1][2 * base]);
        let (q1, r1) = Self::division(f, &tree[k - 1][2 * base + 1]);

        let (u0, u1) = u.split_at(n / 2);
        // println!("u0={}", scalar_vector_to_string(&u0.to_vec()));
        // println!("u1={}", scalar_vector_to_string(&u1.to_vec()));
        let mut rs0: Vec<Scalar> = Self::eval_rec(tree, k - 1, base * 2, &r0, u0);
        let mut rs1: Vec<Scalar> = Self::eval_rec(tree, k - 1, base * 2 + 1, &r1, u1);
        rs0.append(&mut rs1);
        rs0
    }

    // TODO: https://en.wikipedia.org/wiki/Synthetic_division
    fn division(dividend: &[Scalar], divisor: &[Scalar]) -> (Vec<Scalar>, Vec<Scalar>) {
        // Important: if dividend.len() < divisor.len(), then the quotient is zero and
        // the remainder is the dividend.
        if dividend.len() < divisor.len() {
            return (vec![Scalar::zero()], dividend.to_vec());
        }

        let mut quotient = vec![Scalar::zero(); dividend.len() - divisor.len() + 1];
        let mut remainder = dividend.to_vec();

        for i in (0..quotient.len()).rev() {
            quotient[i] =
                remainder[i + divisor.len() - 1] * divisor[divisor.len() - 1].invert().unwrap();
            for j in 0..divisor.len() {
                remainder[i + j] -= quotient[i] * divisor[j];
            }
        }

        // Remove leading zeros
        while remainder.len() > 1 && remainder[remainder.len() - 1] == Scalar::zero() {
            remainder.pop();
        }

        (quotient, remainder)
    }

    /// Compute subproduct tree in O(M(n) * log(n)) time, where O(M(n)) is the
    /// asymptotic complexity of multiplication, and equal to O(nlog(n)) if using
    /// FFT-based fast multiplication.
    ///
    /// Return a vector of levels, each level is a vector of polynomials,
    /// and each polynomial is a vector of coefficients.
    fn contruct_subproduct_tree(domain: &[Scalar]) -> Vec<Vec<Vec<Scalar>>> {
        let n = domain.len();
        assert!(n.is_power_of_two());

        let mut tree = Vec::new();
        let mut level = Vec::new();
        for u in domain.iter() {
            level.push(vec![-*u, Scalar::one()]);
        }
        tree.push(level.clone());

        for k in 0..log_2(n) {
            let mut new_level = Vec::new();
            for i in 0..(n >> (k + 1)) {
                let left = &level[2 * i];
                let right = &level[2 * i + 1];
                let poly = Self::multiplication(left, right);
                new_level.push(poly);
            }
            tree.push(new_level.clone());
            level = new_level;
        }
        assert_eq!(tree.len(), log_2(n) + 1);

        // for i in 0..tree.len() {
        //     println!("tree[{}]=", i);
        //     let level = &tree[i];
        //     for j in 0..level.len() {
        //         println!("tree[{}][{}]={}", i, j, scalar_vector_to_string(&level[j]));
        //     }
        // }

        tree
    }

    /// Compute f(X) = ∑i c_i * z(X) / (X - u_i) in O(M(n) * log(n)) time
    ///
    /// The algorithm is the core of fast interpolation,
    /// from Modern Computer Algebra, Algorithm 10.9, "Linear Combination for linear moduli".
    ///
    /// # Arguments
    ///
    /// - k: `log(n)`
    /// - base: the base index of the subproduct tree
    /// - c: a vector of coefficients of size n
    /// - u: a vector of points of size n
    ///
    /// return: a polynomial in coefficients
    ///
    /// # Example
    ///
    /// TODO: I use this alg. to compute `z'(X)`. Is there
    /// any faster algorithm? It is mentioned in [TAB20] that `z'(X)`
    /// can be computed in O(n) time.
    ///
    fn linear_combination_linear_moduli_fix(
        tree: &Vec<Vec<Vec<Scalar>>>,
        k: usize,
        base: usize,
        c: &[Scalar],
        u: &[Scalar],
    ) -> Vec<Scalar> {
        let n = u.len();

        // println!("lc_fix> k={}, base={}, n={}", k, base, n);
        // println!("c={}", scalar_vector_to_string(&c.to_vec()));
        // println!("u={}", scalar_vector_to_string(&u.to_vec()));

        if k == 0 {
            assert_eq!(c.len(), 1);
            assert_eq!(u.len(), 1);
            return vec![c[0]];
        }

        assert!(n.is_power_of_two());
        assert_eq!(n, pow_2(k));

        let node0 = &tree[k - 1][2 * base];
        let node1 = &tree[k - 1][2 * base + 1];

        let (c0, c1) = c.split_at(n / 2);
        let (u0, u1) = u.split_at(n / 2);

        let r0 = Self::linear_combination_linear_moduli_fix(tree, k - 1, 2 * base, c0, u0);
        let r1 = Self::linear_combination_linear_moduli_fix(tree, k - 1, 2 * base + 1, c1, u1);

        let poly0 = UnivariatePolynomial::monomial(Self::multiplication(node1, &r0));
        let poly1 = UnivariatePolynomial::monomial(Self::multiplication(node0, &r1));
        (&poly0 + &poly1).coeffs().to_vec()
    }

    /// Polynomial interpolation in O(M(n) * log(n)) time.
    /// If O(M(n)) is O(nlog(n)) by using fast polynomial multiplication,
    /// then the asymtotics is O(nlog^2(n)).
    ///
    /// # Arguments
    ///
    /// - evals: a vector of evaluations of size n
    /// - domain: a domain of size n, n must be a power of 2
    ///
    /// # Return
    ///
    /// - a polynomial (coefficients) of degree (n-1)
    ///
    pub fn compute_coeffs_from_evals_fast_2(evals: &[Scalar], domain: &[Scalar]) -> Vec<Scalar> {
        let n = domain.len();
        assert!(n.is_power_of_two());
        assert_eq!(evals.len(), n);

        // 1. building up subproduct tree

        let tree = Self::contruct_subproduct_tree(domain);
        // 2. construct a polynomial with linear moduli

        let f_derivative = Self::linear_combination_linear_moduli_fix(
            &tree,
            log_2(n),
            0,
            &vec![Scalar::one(); n],
            domain,
        );

        let f_derivative_at_u = Self::eval_rec(&tree, log_2(n), 0, &f_derivative, domain);

        let mut bary_centric_weights: Vec<Scalar> = f_derivative_at_u
            .iter()
            .map(|e| e.invert().unwrap())
            .collect();

        let bary_centric_weights2: Vec<Scalar> = f_derivative_at_u
            .iter()
            .map(|e| e.invert().unwrap())
            .collect();
        assert_eq!(bary_centric_weights, bary_centric_weights2);
        // println!("bary_centric_weights={}", scalar_vector_to_string(&bary_centric_weights));
        // println!("bary_centric_weights2={}", scalar_vector_to_string(&UniPolynomial::barycentric_weights(domain)));

        bary_centric_weights
            .iter_mut()
            .enumerate()
            .for_each(|(i, w)| *w *= evals[i]);

        // println!("before linear_combination_linear_moduli_fix");
        

        // {
        //     let z_poly = UniPolynomial::from_coeffs(&tree[log_2(n)][0]);
        //     println!("z_poly={}", scalar_vector_to_string(&z_poly.coeffs));
        //     let mut g = Self::zero();
        //     for i in 0..n {
        //         let (gi, _) = z_poly.div_by_linear_divisor(&domain[i]);
        //         g = Self::add(&g, &gi, false);
        //     }
        //     println!("g={}", scalar_vector_to_string(&g.coeffs));
        // }

        Self::linear_combination_linear_moduli_fix(
            &tree,
            log_2(n),
            0,
            &bary_centric_weights,
            domain,
        )
    }

    pub fn prove(&self, lookup: &Vec<Fr>) -> Vec<u8> {
        let table = self.table.clone();
        let param = self.param.clone();
        let pp = self.pp.clone();
        let d = self.d;

        let m = lookup.len();
        let t = table.len();

        /************
          Round 1
        ************/
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
        ξ = H_I[col_i] = [ω^2, ω^6, ω^2, ω^3]
        Interpolation with ξ and get polynomial: ξ(x)
        */

        // initialize transcript
        let mut transcript = Keccak256Transcript::new(());
        // φ(x)
        let phi_poly = UnivariatePolynomial::lagrange(lookup.clone()).ifft();
        // t(x)
        let t_poly = UnivariatePolynomial::lagrange(table.clone()).ifft();

        // remove duplicated elements
        let t_values_from_lookup: HashSet<_> = lookup.clone().into_iter().collect();
        // I: the index of t_values_from_lookup elements in sub table t_I
        let i_values: Vec<_> = t_values_from_lookup
            .iter()
            .map(|elem| table.iter().position(|&x| x == *elem).unwrap())
            .collect();
        let log_m = log_2(m);
        let v_root_of_unity = root_of_unity::<Fr>(log_m);
        // cache all roots of unity
        let log_t = log_2(t);
        let t_root_of_unity = root_of_unity::<Fr>(log_t);
        let t_roots_of_unity = (0..t)
            .map(|i| t_root_of_unity.pow([i as u64]))
            .collect::<Vec<Fr>>();
        let v_roots_of_unity = (0..m)
            .map(|i| v_root_of_unity.pow([i as u64]))
            .collect::<Vec<Fr>>();
        // H_I = {ξ_i} , i = [1...k], ξ(Xi)
        let h_i: Vec<_> = i_values.iter().map(|&i| t_roots_of_unity[i]).collect();
        // TODO: optimize interpolation polynomial with https://github.com/gy001/hypercube/blob/main/univarization/src/unipoly.rs#L391
        // refer to barycentric_weights in arithmetic.rs
        let t_values_from_lookup_set: Vec<Fr> = t_values_from_lookup.clone().into_iter().collect();
        // let t_i_poly = lagrange_interp(&h_i, &t_values_from_lookup_set);
        let t_i_poly_coeffs =
            Self::compute_coeffs_from_evals_fast_2(&t_values_from_lookup_set, &h_i);
        let t_i_poly = UnivariatePolynomial::monomial(t_i_poly_coeffs);
        let z_i_poly = UnivariatePolynomial::vanishing(&h_i, Fr::one());
        // sanity check
        for (i, &root) in h_i.iter().enumerate() {
            assert_eq!(z_i_poly.evaluate(&root), Fr::zero());
            assert_eq!(t_i_poly.evaluate(&root), t_values_from_lookup_set[i]);
            assert_eq!(t_i_poly.evaluate(&root), t_poly.evaluate(&root));
        }

        let mut col_values = Vec::new();
        let mut v_values = Vec::new();
        for i in 0..m {
            // find the index of 1 in jth row of M
            let col_i = t_values_from_lookup
                .iter()
                .position(|&x| x == lookup[i])
                .unwrap();
            col_values.push(col_i);
            let col_i_root = h_i[col_i];
            // Note: v = 1 / col_i_root in paper
            // Here we use different construction that does not affect the verification
            let v = col_i_root;
            v_values.push(v);
        }
        // ξ(x) polynomial
        let v_poly = UnivariatePolynomial::lagrange(v_values.clone()).ifft();

        // [ξ(x)]1
        let v_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &v_poly, &mut transcript).unwrap();
        // [z_I(x)]2
        let z_i_comm_2: UnivariateKzgCommitment<G2Affine> =
            Pcs::commit_monomial_g2(&param, z_i_poly.coeffs());
        transcript
            .write_commitment_g2(&z_i_comm_2.clone().to_affine())
            .unwrap();
        // [t(x)]1
        let t_i_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &t_i_poly, &mut transcript).unwrap();

        let alpha = transcript.squeeze_challenge();

        // π1 = ([ξ(x)]1, [z_I(x)]2, [t(x)]1)
        let pi_1 = (v_comm_1.clone(), z_i_comm_2.clone(), t_i_comm_1.clone());

        /************
          Round 2
        ************/
        /*
        Calculate μ_i(X), i = [0, m - 1], by interpolating it in V(multiplicative subgroup, order is m)
        col_i = col[i]
        v_i = V[i] # ω^i
        v_col_i = V[col_i] # ω^col_i
        μ_i(X) = z_V(X) / (z_V'(v_i) * (X - v_i)) = v_i / m * (X^m - 1) / (X - v_i)

        Calculate Normalized Lagrange Polynomial: τ_col(i)(X) / τ_col(i)(0):
        ξ_{col(i)} = H_I[col_i] # h_i
        normalized_lag_poly = τ_col(i)(X) / τ_col(i)(0)
                            = z_I(X) / z_I(0) * (0 - ξ_{col(i)}) / (X - ξ_{col(i)})

        Calculate D(X), E(X)
        D(X) = Σ_i(μ_i(α) * normalized_lag_poly)
        E(X) = Σ_i(μ_i(X) * normalized_lag_poly(β))

        Calculate Q_D(X) and R(X): Theorem 5 (Inner Product Polynomial Relation) on Baloo paper
        D(X) * t_I(X) - φ(α) - R(X) = z_I(X) * Q_D(X)
        Q_D(X), R(X) = (D(X) * t_I(X) - φ(α)) / z_I(X) # polynomial division with remainder

        Calculate Q_E(X)
        1) Baloo paper uses this construction
        v_i = 1 / H_I[col_i]
        v(X): interpolate polynomial with v_i values
        E(X) * (βv(X) - 1) + z_I(β) / z_I(0) = z_V(X) * Q_E(X)
        Q_E(X) = (E(X) * (βv(X) - 1) + z_I(β) / z_I(0)) / z_V(X)

        2) Our code uses this optimized construction:
        v_i = H_I[col_i]
        v(X): interpolate polynomial with v_i values
        E(X) * (β - v(X)) + v(X) * z_I(β) / z_I(0) = z_V(X) * Q_E(X)
        Q_E(X) = (E(X) * (β - v(X)) + v(X) * z_I(β) / z_I(0)) / z_V(X)
        */

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
        for i in 0..m {
            // col(i)
            let col_i = col_values[i];
            // ω^i
            let v_root = v_roots_of_unity[i];
            // X - ω^i
            let v_root_poly = UnivariatePolynomial::monomial(vec![-v_root, scalar_1]);
            // ξ_i
            let col_i_root = h_i[col_i];
            // X - ξ_i
            let x_root_poly = UnivariatePolynomial::monomial(vec![-col_i_root, scalar_1]);
            // Lagrange polynomial on V: μ_i(X)
            // z_v_poly / v_root_poly * v_root / Fr::from(m as u64);
            let mu_poly =
                &(&z_v_poly / &v_root_poly) * (v_root * (Fr::from(m as u64).invert().unwrap()));
            // Normalized Lagrange Polynomial: τ_col(i)(X) / τ_col(i)(0)
            // z_i_poly / x_root_poly * (-col_i_root) / z_i_at_0;
            let normalized_lag_poly =
                &(&z_i_poly / &x_root_poly) * (col_i_root.neg() * (z_i_at_0.invert().unwrap()));
            // μ_i(α)
            let mu_poly_at_alpha = mu_poly.evaluate(&alpha);
            // D(X) = Σ_i(μ_i(α) * normalized_lag_poly)
            d_poly += &normalized_lag_poly * mu_poly_at_alpha;
        }

        // D(X) * t_I(X)
        let d_t_poly = d_poly.poly_mul(t_i_poly.clone());
        // φ(α)
        let phi_poly_at_alpha = phi_poly.evaluate(&alpha);

        // Q_D(X), R(X) = (D(X) * t_I(X) - φ(α)) / z_I(X)
        let (q_d_poly, r_poly) = (d_t_poly + phi_poly_at_alpha.neg()).div_rem(&z_i_poly);
        assert_eq!(r_poly.evaluate(&scalar_0), scalar_0);

        // π2 = ([D]1 = [D(x)]1, [R]1 = [R(x)]1, [Q2]1 = [Q_D(x)]1)
        let d_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &d_poly, &mut transcript).unwrap();
        let r_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &r_poly, &mut transcript).unwrap();
        let q_d_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &q_d_poly, &mut transcript).unwrap();

        let beta = transcript.squeeze_challenge();

        // calculate E(X) = Σ_i(μ_i(X) * normalized_lag_poly(β))
        let mut e_poly: UnivariatePolynomial<Fr> = zero_poly.clone();
        for i in 0..m {
            // col(i)
            let col_i = col_values[i];
            // ω^i
            let v_root = v_roots_of_unity[i];
            // X - ω^i
            let v_root_poly = UnivariatePolynomial::monomial(vec![-v_root, scalar_1]);
            // ξ_i
            let col_i_root = h_i[col_i];
            // X - ξ_i
            let x_root_poly = UnivariatePolynomial::monomial(vec![-col_i_root, scalar_1]);
            // Lagrange polynomial on V: μ_i(X)
            // z_v_poly / v_root_poly * v_root / Fr::from(m as u64);
            let mu_poly =
                &(&z_v_poly / &v_root_poly) * (v_root * (Fr::from(m as u64).invert().unwrap()));
            // Normalized Lagrange Polynomial: τ_col(i)(X) / τ_col(i)(0)
            // z_i_poly / x_root_poly * (-col_i_root) / z_i_at_0;
            let normalized_lag_poly =
                &(&z_i_poly / &x_root_poly) * (col_i_root.neg() * (z_i_at_0.invert().unwrap()));
            // Normalized Lagrange Polynomial at β: τ_col(i)(β) / τ_col(i)(0)
            let normalized_lag_poly_at_beta = normalized_lag_poly.evaluate(&beta);
            // E(X) = Σ_i(μ_i(X) * normalized_lag_poly(β))
            e_poly += &mu_poly * normalized_lag_poly_at_beta;
        }

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

        // π3 = ([E]1 = [E(x)]1, [Q1]1 = [Q_E(x)]1)
        let e_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &e_poly, &mut transcript).unwrap();
        let q_e_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &q_e_poly, &mut transcript).unwrap();

        let gamma: Fr = transcript.squeeze_challenge();
        let zeta: Fr = transcript.squeeze_challenge();
        let gamma_2 = gamma.mul(gamma);
        let gamma_3 = gamma_2.mul(gamma);

        // π2 = ([D]1, [R]1, [Q2]1)
        let pi_2 = (d_comm_1.clone(), r_comm_1.clone(), q_d_comm_1.clone());
        // π3 = ([E]1, [Q1]1)
        let pi_3 = (e_comm_1.clone(), q_e_comm_1.clone());

        /************
          Round 3: optimize with linear combination of polynomials
        ************/
        /*
        Calculate w1, w2, w3, w4
        */

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
        let p_e_poly = &(&(&beta_sub_v_poly * v5) + &v_poly * (v4.mul(v3.invert().unwrap())))
            - &q_e_poly * z_v_zeta;
        // X^(d-m+1)
        let coeffs = vec![scalar_0; d - m + 1]
            .into_iter()
            .chain(vec![scalar_1])
            .collect();
        let x_exponent_poly = UnivariatePolynomial::monomial(coeffs);
        // calculate [w1]1, [w2]1, [w2]1, [w4]1
        // X - α
        let x_alpha_poly = UnivariatePolynomial::monomial(vec![-alpha, scalar_1]);
        // calculate w1 = X^(d-m+1) * (E(X) - E(α) + (φ(X) - φ(α))γ) / X - α
        let mut w1 = &(&(e_poly.clone() + v1.neg()) + &(phi_poly.clone() + v2.neg()) * gamma)
            / &x_alpha_poly;
        w1 = w1.poly_mul(x_exponent_poly.clone());
        // calculate polynomial X
        let x_poly = UnivariatePolynomial::monomial(vec![scalar_0, scalar_1]);
        // X^m
        let x_m_exponent_poly = UnivariatePolynomial::monomial(
            vec![scalar_0; m]
                .into_iter()
                .chain(vec![scalar_1])
                .collect(),
        );
        // calculate w2 = (z_I(X) - v3 / X + γ * R(X) / X +  γ^2 * X^(d-m+1) * (z_I(X) - X^m)) + γ^3 * X^(d-m+1) * R(X)
        let w2 = &(&(&(z_i_poly.clone() + v3.neg()) / &x_poly) + &(&(&r_poly * gamma) / &x_poly))
            + &x_exponent_poly.poly_mul(
                &(&(&z_i_poly - x_m_exponent_poly.clone()) * gamma_2) + &r_poly * gamma_3,
            );
        // calculate X - β
        let x_beta_poly = UnivariatePolynomial::monomial(vec![-beta, scalar_1]);
        // calculate w3 = (D(X) - E(α) + (z_I(X) - z_I(β))γ + P_D(X)γ^2) / X - β
        // v1 = E(α) == D(β)
        let w3 = &(&(&(d_poly + v1.neg()) + &(z_i_poly.clone() + v4.neg()) * gamma)
            + &p_d_poly * gamma_2)
            / &x_beta_poly;
        // calculate X - ζ
        let x_zeta_poly = UnivariatePolynomial::monomial(vec![-zeta, scalar_1]);
        // calculate w4 = (E(X) - E(ζ) + P_E(X)γ) / X - ζ
        // v5 = E(ζ)
        let w4 = &(&(e_poly + v5.neg()) + &p_e_poly * gamma) / &x_zeta_poly;

        // caulk+ calculate w5, w6
        // z_h_poly = X^t - 1, [-1, 0, ..., 0, 1], t-1 0s in between
        let z_h_poly_coeffs = vec![scalar_1.neg()]
            .into_iter()
            .chain(vec![scalar_0; t - 1])
            .chain(vec![scalar_1])
            .collect();
        let z_h_poly = UnivariatePolynomial::monomial(z_h_poly_coeffs);

        // calculate barycentric_weights
        let bc_weights = barycentric_weights(&h_i);

        // w5_poly = (t_poly - t_I_poly) / z_I_poly
        let w5_poly_direct = &(&t_poly - &t_i_poly) / &z_i_poly;
        // q_t_poly_i = (t_poly - table[i])/X-root_of_unity^i
        let q_t_polys: Vec<_> = table
            .clone()
            .into_iter()
            .enumerate()
            .map(|(i, point)| {
                &(t_poly.clone() + point.neg())
                    / &(x_poly.clone() + t_root_of_unity.pow([i as u64]).neg())
            })
            .collect();
        // optimize w5_poly = bc_weights[0] * q_t_polys[i_values[0]] + bc_weights[1] * q_t_polys[i_values[1]] + ... + bc_weights[h_i.len()-1] * q_t_polys[i_values[h_i.len()-1]]
        let w5_poly = bc_weights
            .clone()
            .into_iter()
            .enumerate()
            .map(|(i, weight)| &q_t_polys[i_values[i]] * weight)
            .reduce(|acc, poly| &acc + poly)
            .unwrap();
        assert_eq!(w5_poly, w5_poly_direct);
        // w6_poly = z_H_poly * (bc_weights[0] * 1/X-h_i[0] + bc_weights[1] * 1/X-h_i[1] + ... + bc_weights[h_i.len()-1] * 1/X-h_i[h_i.len()-1])
        let denom_polys: Vec<_> = h_i
            .into_iter()
            .map(|root| UnivariatePolynomial::monomial(vec![root.neg(), scalar_1]))
            .collect();
        let w6_poly = bc_weights
            .clone()
            .into_iter()
            .enumerate()
            .map(|(i, weight)| &(&z_h_poly / &denom_polys[i]) * weight)
            .reduce(|acc, poly| &acc + poly)
            .unwrap();
        // w6_poly = z_H_poly / z_I_poly
        let w6_poly_direct = &z_h_poly / &z_i_poly;
        assert_eq!(w6_poly, w6_poly_direct);

        // Compress Caulk+ proof.
        let a_poly = &w5_poly.clone() + &w6_poly.clone() * gamma;

        // write v1, v2, v3, v4, v5 to transcript
        transcript.write_field_element(&v1).unwrap();
        transcript.write_field_element(&v2).unwrap();
        transcript.write_field_element(&v3).unwrap();
        transcript.write_field_element(&v4).unwrap();
        transcript.write_field_element(&v5).unwrap();

        // [a]1
        let a_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &a_poly, &mut transcript).unwrap();
        // calculate [w1]1, [w2]1, [w2]1, [w4]1 and write to transcript
        let w1_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &w1, &mut transcript).unwrap();
        let w2_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &w2, &mut transcript).unwrap();
        let w3_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &w3, &mut transcript).unwrap();
        let w4_comm_1: UnivariateKzgCommitment<G1Affine> =
            Pcs::commit_and_write(&pp, &w4, &mut transcript).unwrap();

        // π4 = (v1, v2, v3, v4, v5, [a]1, [w1]1, [w2]1, [w3]1, [w4]1)
        let pi_4 = (
            v1,
            v2,
            v3,
            v4,
            v5,
            a_comm_1.clone(),
            w1_comm_1.clone(),
            w2_comm_1.clone(),
            w3_comm_1.clone(),
            w4_comm_1.clone(),
        );

        // generate proof from transcript
        transcript.into_proof()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::transcript::{
        FieldTranscriptRead, FieldTranscriptWrite, G2TranscriptRead, G2TranscriptWrite,
    };
    use halo2_curves::bn256::Fr;
    type Pcs = UnivariateKzg<Bn256>;

    #[test]
    fn test_baloo() {
        let table = vec![
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
            Fr::from(6),
            Fr::from(7),
            Fr::from(8),
        ];
        let lookup = vec![Fr::from(4), Fr::from(3), Fr::from(5), Fr::from(2)];
        let m = lookup.len();
        let t = table.len();
        // 1. setup
        let (param, pp, _) = preprocess(t, m).unwrap();
        // 2. generate proof
        let prover = Prover::new(&table, &param, &pp);
        let proof = prover.prove(&lookup);
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
            let poly =
                <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::monomial(lookup.clone());
            print!("coeffs: {:?}\n", poly.coeffs());
            let comm = Pcs::commit_and_write(&pp, &poly, &mut transcript).unwrap();
            let point = <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::squeeze_point(
                m,
                &mut transcript,
            );
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
                &<Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::squeeze_point(
                    m,
                    &mut transcript,
                ),
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

        let alpha: Fr = transcript.squeeze_challenge();

        let comm = <Pcs as PolynomialCommitmentScheme<Fr>>::commit(&pp, &poly).unwrap();
        println!("comm: {:?}", comm);
        println!("comm affine: {:?}", comm.clone().to_affine());
        transcript
            .write_commitment(&comm.clone().to_affine())
            .unwrap();
        let comm_1 = Pcs::commit_monomial(&pp, &poly.coeffs());
        println!("comm_1: {:?}", comm_1);
        assert_eq!(comm_1, comm);

        let beta: Fr = transcript.squeeze_challenge();

        // g2
        let comm_2 = Pcs::commit_monomial_g2(&param, &poly.coeffs());
        println!("comm_2: {:?}", comm_2);
        let comm_2_affine = comm_2.clone().to_affine();
        transcript.write_commitment_g2(&comm_2_affine).unwrap();

        let comm_table = Pcs::commit_and_write(&pp, &poly_table, &mut transcript).unwrap();
        println!("comm_table: {:?}", comm_table);

        let zeta: Fr = transcript.squeeze_challenge();

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

        let alpha_verifier: Fr = transcript.squeeze_challenge();
        assert_eq!(alpha, alpha_verifier);

        let comm_verifier = Pcs::read_commitment(&vp, &mut transcript).unwrap();
        println!("comm_verifier: {:?}", comm_verifier);

        let beta_verifier: Fr = transcript.squeeze_challenge();
        assert_eq!(beta, beta_verifier);
        // g2
        let comm_table_g2: G2Affine = transcript.read_commitment_g2().unwrap();
        println!("comm_table_g2: {:?}", comm_table_g2);

        let comm_table_verifier = Pcs::read_commitment(&vp, &mut transcript).unwrap();
        println!("comm_table_verifier: {:?}", comm_table_verifier);
        assert_eq!(comm, comm_verifier);
        assert_eq!(comm_table, comm_table_verifier);
        assert_eq!(comm_2.to_affine(), comm_table_g2);

        let zeta_verifier: Fr = transcript.squeeze_challenge();
        assert_eq!(zeta, zeta_verifier);
    }
}
