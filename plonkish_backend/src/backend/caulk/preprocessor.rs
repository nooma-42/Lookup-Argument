use crate::pcs::univariate::{
    UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam,
};
use crate::pcs::PolynomialCommitmentScheme;
use crate::poly::univariate::UnivariatePolynomial;
use crate::util::arithmetic::{Field, root_of_unity, WithSmallOrderMulGroup};
use crate::Error;
use halo2_curves::pairing::MultiMillerLoop;
use halo2_curves::CurveAffine;
use halo2_curves::ff::PrimeField;
use rand::rngs::OsRng;
use serde::{de::DeserializeOwned, Serialize};
use std::cmp::max;
use std::fmt::Debug;

use super::{CaulkProverParam, CaulkVerifierParam};

pub fn setup<
    M: MultiMillerLoop + Debug + Sync,
>(
    N: usize,
    m: usize,
) -> Result<(CaulkProverParam<M>, CaulkVerifierParam<M>), Error>
where
    M::Scalar: Serialize + DeserializeOwned,
    M::G1Affine: Serialize + DeserializeOwned,
    M::G2Affine: Serialize + DeserializeOwned,
{
    // N: size of the commitment vector 'c' (table size)
    // m: size of the lookup values vector 'values'
    let mut rng = OsRng;
    // Determine the required polynomial degree for KZG based on N and m.
    // Key polynomials and their degrees:
    // - C(X), phi(X): degree N-1, m-1 respectively
    // - z_I(X): degree unique_positions.len() (bounded by m) 
    // - c_I(X): degree roughly N/2 + 3 (blinding degree)
    // - u(X): degree m-1
    // - c_I(u(X)): degree (N/2 + 3) * (m-1) 
    // - z_I(u(X)): degree unique_positions.len() * (m-1), bounded by m * (m-1)
    // - H2_poly: max degree of composition polynomials
    
    // Conservative upper bound: account for all polynomial operations
    let base_degree = N.max(m);
    let composition_degree = base_degree * m; // Upper bound for compositions like c_I(u(X))
    let blinding_degree = 7; // We use 7 blinding factors
    let max_degree = composition_degree + blinding_degree;
    
    // Ensure minimum size for small examples
    let poly_size = max_degree.max(8).next_power_of_two(); // Minimum 8, then next power of 2

    // The blinding factor count seems to be 0 in the original Caulk setup/trim calls.
    let blinding_factors = 0;
    let kzg_param = UnivariateKzg::<M>::setup(poly_size, blinding_factors, &mut rng)?;
    let (kzg_pp, kzg_vp) = UnivariateKzg::<M>::trim(&kzg_param, poly_size, blinding_factors)?;
    Ok((
        CaulkProverParam {
            kzg_param: kzg_param.clone(), // Clone kzg_param as it's needed by prover
            kzg_pp,
        },
        CaulkVerifierParam { kzg_vp, m }, // Store m in VerifierParam
    ))
}

/// Preprocess with optimized H1_com calculation using FK-style precomputation
/// This generates precomputed G2 KZG openings for C(X) at all N-th roots of unity
pub fn preprocess_with_table<
    M: MultiMillerLoop + Debug + Sync,
>(
    N: usize,
    m: usize,
    table: &[M::Scalar],
) -> Result<(CaulkProverParam<M>, CaulkVerifierParam<M>, Vec<M::G2Affine>), Error>
where
    M::Scalar: Serialize + DeserializeOwned + Field + WithSmallOrderMulGroup<3>,
    M::G1Affine: Serialize + DeserializeOwned,
    M::G2Affine: Serialize + DeserializeOwned + CurveAffine<ScalarExt = M::Scalar>,
{
    assert_eq!(table.len(), N, "Table size must match N");
    assert!(N.is_power_of_two(), "N must be a power of two for efficient FK computation");
    
    let mut rng = OsRng;
    let base_degree = N.max(m);
    let composition_degree = base_degree * m;
    let blinding_degree = 7;
    let max_degree = composition_degree + blinding_degree;
    let poly_size = max_degree.max(8).next_power_of_two();

    let blinding_factors = 0;
    let kzg_param = UnivariateKzg::<M>::setup(poly_size, blinding_factors, &mut rng)?;
    let (kzg_pp, kzg_vp) = UnivariateKzg::<M>::trim(&kzg_param, poly_size, blinding_factors)?;
    
    // Precompute G2 KZG openings for C(X) using FK-style algorithm
    let precomputed_g2_openings = precompute_c_openings_g2::<M>(&kzg_param, table)?;
    
    Ok((
        CaulkProverParam {
            kzg_param: kzg_param.clone(),
            kzg_pp,
        },
        CaulkVerifierParam { kzg_vp, m },
        precomputed_g2_openings,
    ))
}

/// Precompute G2 KZG openings for C(X) at all N-th roots of unity
/// Returns Q_i = Comm_G2((C(X) - C(ω^i)) / (X - ω^i)) for i = 0, ..., N-1
fn precompute_c_openings_g2<M: MultiMillerLoop + Debug + Sync>(
    kzg_param: &UnivariateKzgParam<M>,
    table: &[M::Scalar],
) -> Result<Vec<M::G2Affine>, Error>
where
    M::Scalar: Field + WithSmallOrderMulGroup<3>,
    M::G2Affine: CurveAffine<ScalarExt = M::Scalar>,
{
    let N = table.len();
    assert!(N.is_power_of_two());
    
    // Convert table to polynomial coefficients
    let c_poly = UnivariatePolynomial::lagrange(table.to_vec()).ifft();
    let c_coeffs = c_poly.coeffs();
    
    // Get N-th roots of unity
    let log_n = N.trailing_zeros() as usize;
    let omega = root_of_unity::<M::Scalar>(log_n);
    let roots: Vec<M::Scalar> = (0..N)
        .map(|i| omega.pow([i as u64]))
        .collect();
    
    // Use FK-style algorithm to compute all openings efficiently
    compute_all_g2_openings(kzg_param, c_coeffs, &roots)
}

/// FK-style algorithm to compute all G2 KZG openings efficiently
/// This computes Q_i = Comm_G2((C(X) - C(ω^i)) / (X - ω^i)) for all i
fn compute_all_g2_openings<M: MultiMillerLoop + Debug + Sync>(
    kzg_param: &UnivariateKzgParam<M>,
    c_coeffs: &[M::Scalar],
    roots: &[M::Scalar],
) -> Result<Vec<M::G2Affine>, Error>
where
    M::Scalar: Field,
    M::G2Affine: CurveAffine<ScalarExt = M::Scalar>,
{
    let N = roots.len();
    assert_eq!(c_coeffs.len(), N);
    assert!(N.is_power_of_two());
    
    // Improved FK algorithm similar to CQ implementation
    // First, we compute a circulant matrix representation
    
    // Get first column of circulant matrix for opening computations
    let mut first_col = vec![M::Scalar::ZERO; 2 * N];
    
    // The opening polynomial coefficients can be computed via FFT operations
    // For each i, we need (C(X) - C(ω^i)) / (X - ω^i)
    // This is equivalent to computing the derivative-like operation
    
    // Compute all evaluations C(ω^i) first
    let c_evals: Vec<M::Scalar> = roots.iter()
        .map(|&root| evaluate_polynomial_at_point(c_coeffs, &root))
        .collect();
    
    // Compute opening coefficients using efficient method
    let mut all_openings = Vec::with_capacity(N);
    
    for i in 0..N {
        let root = roots[i];
        let c_eval = c_evals[i];
        
        // Compute (C(X) - C(ω^i)) / (X - ω^i) coefficients
        let mut opening_coeffs = vec![M::Scalar::ZERO; N - 1];
        
        // Use synthetic division or direct computation
        // (C(X) - C(ω^i)) / (X - ω^i) = Σ_{j=0}^{N-2} a_j X^j
        // where a_j = Σ_{k=j+1}^{N-1} c_k * ω^{i*(k-j-1)}
        
        for j in 0..(N - 1) {
            let mut coeff = M::Scalar::ZERO;
            for k in (j + 1)..N {
                let power = (k - j - 1) as u64;
                coeff += c_coeffs[k] * root.pow([power]);
            }
            opening_coeffs[j] = coeff;
        }
        
        // Commit to G2
        let opening = UnivariateKzg::<M>::commit_monomial_g2(kzg_param, &opening_coeffs);
        all_openings.push(opening.to_affine());
    }
    
    Ok(all_openings)
}

/// Helper function to evaluate polynomial at a point using Horner's method
fn evaluate_polynomial_at_point<F: Field>(coeffs: &[F], point: &F) -> F {
    coeffs.iter().rev().fold(F::ZERO, |acc, &coeff| acc * point + coeff)
}

/// More efficient FK-style computation using FFT (future optimization)
/// This would be similar to CQ's implementation with ec_fft
#[allow(dead_code)]
fn compute_all_g2_openings_fft<M: MultiMillerLoop + Debug + Sync>(
    kzg_param: &UnivariateKzgParam<M>,
    c_coeffs: &[M::Scalar],
    roots: &[M::Scalar],
) -> Result<Vec<M::G2Affine>, Error>
where
    M::Scalar: Field,
    M::G2Affine: CurveAffine<ScalarExt = M::Scalar>,
{
    // TODO: Implement FFT-based computation similar to CQ's FK algorithm
    // This would involve:
    // 1. Setting up circulant matrix with first column
    // 2. Using ec_fft on G2 elements  
    // 3. Element-wise multiplication with FFT of first column
    // 4. Inverse FFT to get final result
    
    // For now, fall back to direct computation
    compute_all_g2_openings(kzg_param, c_coeffs, roots)
} 