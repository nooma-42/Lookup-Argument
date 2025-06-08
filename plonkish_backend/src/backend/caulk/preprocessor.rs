// src/preprocessor.rs

use crate::pcs::univariate::{
    UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam,
};
use crate::pcs::PolynomialCommitmentScheme;
use crate::poly::univariate::UnivariatePolynomial;
use crate::util::arithmetic::{radix2_fft, root_of_unity, root_of_unity_inv, Field, WithSmallOrderMulGroup};
use crate::Error;
use halo2_curves::ff::PrimeField;
use halo2_curves::group::prime::PrimeCurveAffine;
use halo2_curves::group::Group;
use halo2_curves::pairing::MultiMillerLoop;
use halo2_curves::CurveAffine;
use rand::rngs::OsRng;
use serde::{de::DeserializeOwned, Serialize};
use std::cmp::max;
use std::fmt::Debug;

use super::{CaulkProverParam, CaulkVerifierParam};

// A helper function to determine a safe degree bound for the SRS,
// closely following the logic of the original arkworks implementation.
fn get_degree_bound(n: usize, m: usize) -> usize {
    // The degree of H2(X) and related polynomials can be up to O(m^2).
    // The composition of blinded polynomials like c_I(u(X)) can lead to degrees
    // around m^2 + O(m). For example, (m+2)*(m+2) = m^2+4m+4.
    // We need a bound that covers this.
    // 2*m*m is a safe upper bound used by the original implementation.
    let m_squared_bound = 2 * m.saturating_mul(m);

    // Let's add a small constant padding to handle potential off-by-one degree issues
    // from combinations of blinding factors, especially in small cases.
    let padding = 8;
    
    // The degree bound must be at least N to commit to C(X).
    // We also enforce a minimum degree to handle very small test cases,
    // ensuring there's enough room for blinding polynomials.
    let min_degree = 64; // Increased minimum degree to be safer.

    max(n, m_squared_bound + padding).max(min_degree)
}


pub fn setup<M: MultiMillerLoop + Debug + Sync>(
    N: usize,
    m: usize,
) -> Result<(CaulkProverParam<M>, CaulkVerifierParam<M>), Error>
where
    M::Scalar: Serialize + DeserializeOwned + WithSmallOrderMulGroup<3>,
    M::G1Affine: Serialize + DeserializeOwned,
    M::G2Affine: Serialize + DeserializeOwned,
{
    let mut rng = OsRng;

    // Use a robust method to determine the SRS size.
    let degree_bound = get_degree_bound(N, m);
    let poly_size = degree_bound.next_power_of_two();
    
    let blinding_factors = 0;
    let kzg_param = UnivariateKzg::<M>::setup(poly_size, blinding_factors, &mut rng)?;
    let (kzg_pp, kzg_vp) = UnivariateKzg::<M>::trim(&kzg_param, poly_size, blinding_factors)?;
    
    Ok((
        CaulkProverParam {
            kzg_param: kzg_param.clone(),
            kzg_pp,
        },
        CaulkVerifierParam { kzg_vp, m },
    ))
}

pub fn preprocess_with_table<M: MultiMillerLoop + Debug + Sync>(
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
    assert!(
        N.is_power_of_two(),
        "N must be a power of two for efficient FK computation"
    );

    let mut rng = OsRng;
    
    // Use the same robust method to determine the SRS size.
    let degree_bound = get_degree_bound(N, m);
    let poly_size = degree_bound.next_power_of_two();
    
    let blinding_factors = 0;
    let kzg_param = UnivariateKzg::<M>::setup(poly_size, blinding_factors, &mut rng)?;
    let (kzg_pp, kzg_vp) = UnivariateKzg::<M>::trim(&kzg_param, poly_size, blinding_factors)?;

    let precomputed_g2_openings = precompute_c_openings_g2_optimized::<M>(&kzg_param, table)?;

    Ok((
        CaulkProverParam {
            kzg_param: kzg_param.clone(),
            kzg_pp,
        },
        CaulkVerifierParam { kzg_vp, m },
        precomputed_g2_openings,
    ))
}

// ... The rest of the file remains the same ...

/// Optimized FK-style computation using existing FFT library
/// This implements the circulant matrix approach with G2 FFT operations
fn precompute_c_openings_g2_optimized<M: MultiMillerLoop + Debug + Sync>(
    kzg_param: &UnivariateKzgParam<M>,
    table: &[M::Scalar],
) -> Result<Vec<M::G2Affine>, Error>
where
    M::Scalar: Field + WithSmallOrderMulGroup<3>,
    M::G2Affine: CurveAffine<ScalarExt = M::Scalar>,
{
    let N = table.len();
    assert!(
        N.is_power_of_two(),
        "Table size must be power of two for FK algorithm"
    );

    // Convert table to polynomial coefficients using existing FFT
    let c_poly = UnivariatePolynomial::lagrange(table.to_vec()).ifft();
    let mut c_coeffs = c_poly.coeffs().to_vec();

    // Pad coefficients to match domain size N
    c_coeffs.resize(N, M::Scalar::ZERO);
    
    // Try to use FK algorithm with existing FFT library
    if let Ok(result) = fk_g2::<M>(&mut c_coeffs, kzg_param) {
        return Ok(result);
    }
    
    // Fallback to direct computation if FK fails
    let log_n = N.trailing_zeros() as usize;
    let root_of_unity_n = root_of_unity::<M::Scalar>(log_n);
    let mut results = Vec::with_capacity(N);
    
    // For each i, compute KZG opening at ω^i
    for i in 0..N {
        let omega_i = root_of_unity_n.pow([i as u64]);
        
        // Compute (C(X) - C(ω^i)) / (X - ω^i)
        let c_eval = c_poly.evaluate(&omega_i);
        
        // Create polynomial C(X) - C(ω^i)
        let mut shifted_coeffs = c_poly.coeffs().to_vec();
        if shifted_coeffs.is_empty() {
            shifted_coeffs.push(-c_eval);
        } else {
            shifted_coeffs[0] -= c_eval;
        }
        
        // Divide by (X - ω^i) using synthetic division
        let quotient_coeffs = synthetic_division(&shifted_coeffs, omega_i);
        
        // Commit to the quotient polynomial using G2
        let opening = UnivariateKzg::<M>::commit_monomial_g2(kzg_param, &quotient_coeffs);
        results.push(opening.to_affine());
    }
    
    Ok(results)
}

/// Synthetic division: divide polynomial by (X - root)
/// This is more efficient than using polynomial division for linear divisors
fn synthetic_division<F: Field>(coeffs: &[F], root: F) -> Vec<F> {
    if coeffs.is_empty() {
        return vec![];
    }
    
    if coeffs.len() == 1 {
        // Constant polynomial divided by linear gives zero polynomial
        return vec![];
    }
    
    let mut result = vec![F::ZERO; coeffs.len() - 1];
    
    // Start from the highest degree coefficient
    let last_coeff = *coeffs.last().unwrap();
    let result_len = result.len();
    result[result_len - 1] = last_coeff;
    
    // Work backwards through the coefficients using Horner's method
    for i in (1..coeffs.len() - 1).rev() {
        result[i - 1] = coeffs[i] + result[i] * root;
    }
    
    result
}

/// FK algorithm for G2 commitments using existing FFT library
/// Computes all KZG opening proofs Q_i = Comm_G2((C(X) - C(ω^i)) / (X - ω^i))
fn fk_g2<M: MultiMillerLoop + Debug + Sync>(
    coeffs: &mut Vec<M::Scalar>,
    kzg_param: &UnivariateKzgParam<M>,
) -> Result<Vec<M::G2Affine>, Error>
where
    M::Scalar: Field + WithSmallOrderMulGroup<3>,
    M::G2Affine: CurveAffine<ScalarExt = M::Scalar>,
{
    let N = coeffs.len();
    assert!(N.is_power_of_two(), "Coefficients length must be power of two");

    // Step 1: Get first column of circulant matrix in length 2*N
    let mut first_col = coeffs.to_vec();
    first_col[0] = M::Scalar::ZERO; // Zero out the first element
    
    // Pad with zeros to double size
    first_col.splice(0..0, std::iter::repeat(M::Scalar::ZERO).take(N));
    // Rotate right by 1
    first_col.rotate_right(1);

    // Step 2: Get G2 powers with inverse ordering
    let g2_powers = kzg_param.powers_of_s_g2();
    let mut inv_g2_powers = g2_powers[..N-1].to_vec();
    inv_g2_powers.reverse();
    inv_g2_powers.push(M::G2Affine::identity()); // Add identity element
    
    // Pad with identity elements to double size
    let g2_identity_vals = vec![M::G2Affine::identity(); N];
    inv_g2_powers.extend(g2_identity_vals);

    // Step 3: Apply FK transform using existing FFT library
    let log_size = (2 * N).trailing_zeros() as usize;
    let omega = root_of_unity::<M::Scalar>(log_size);
    let omega_inv = root_of_unity_inv::<M::Scalar>(log_size);
    
    // Right hand side: G2_FFT(inv_g2_powers) using radix2_fft
    let mut rhs = inv_g2_powers.iter().map(|x| (*x).into()).collect::<Vec<M::G2>>();
    radix2_fft(&mut rhs, omega, log_size);
    
    // Middle hand side: FFT(first_col) using radix2_fft
    let mut mhs = first_col;
    radix2_fft(&mut mhs, omega, log_size);
    
    // Element-wise multiplication in frequency domain
    let m_r_hs: Vec<M::G2> = rhs
        .iter()
        .zip(mhs.iter())
        .map(|(&r, &m)| {
            let mut result = r;
            result *= m;
            result
        })
        .collect();
    
    // Step 4: Inverse G2 FFT using radix2_fft
    let mut result = m_r_hs;
    radix2_fft(&mut result, omega_inv, log_size);
    
    // Apply IFFT scaling
    let size_inv = M::Scalar::from((2 * N) as u64).invert().unwrap();
    result.iter_mut().for_each(|x| *x *= size_inv);
    
    // Convert back to affine and extract first N values
    let mut result_affine: Vec<M::G2Affine> = result.into_iter().map(|x| x.into()).collect();
    result_affine.truncate(N);
    
    // Scale by roots of unity
    let log_n = N.trailing_zeros() as usize;
    let root_of_unity_n = root_of_unity::<M::Scalar>(log_n);
    let roots: Vec<M::Scalar> = (0..N)
        .map(|i| root_of_unity_n.pow([i as u64]))
        .collect();
    
    for i in 0..N {
        let scale = roots[i] * (M::Scalar::from(N as u64).invert().unwrap());
        let mut scaled: M::G2 = result_affine[i].into();
        scaled *= scale;
        result_affine[i] = scaled.into();
    }
    
    Ok(result_affine)
}


/// Helper function to check if a number is power of two
fn is_power_of_two(n: usize) -> bool {
    (n & (n - 1)) == 0 && n != 0
}