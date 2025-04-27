// plonkish_backend/src/backend/lasso/util.rs
use crate::{
    poly::multilinear::MultilinearPolynomial,
    util::arithmetic::{Field, PrimeField},
    Error,
};
// Existing functions: calculate_multiplicities, vec_to_mle, hash_tuple

/// Converts a value (assumed to be from the simple range table [0..N-1])
/// back into its decomposed subtable indices. Inverse of the simple g.
/// value = idx^(1) + idx^(2)*2^l + ... + idx^(c)*2^(l*(c-1))
pub(crate) fn get_index(value: usize, l: usize, c: usize) -> Result<Vec<usize>, Error> {
    if l == 0 || c == 0 {
        return Err(Error::InvalidSnark("l and c must be > 0".to_string()));
    }
    let subtable_size = 1 << l;
    let mut val = value;
    let mut index = Vec::with_capacity(c);
    for _ in 0..c {
        index.push(val % subtable_size);
        val /= subtable_size;
    }
     if val != 0 {
         // This means the original value was too large for the N=2^(l*c) range
         return Err(Error::InvalidSnark(format!("Value {} too large for range 2^({})", value, l*c)));
     }
    Ok(index)
}

/// Helper to convert Vec<F> to MLE and return both.
pub(crate) fn get_poly_and_vec<F: PrimeField>(
    vec: Vec<F>,
    num_vars: usize,
) -> Result<(MultilinearPolynomial<F>, Vec<F>), Error> {
    if vec.len() != (1 << num_vars) {
         return Err(Error::InvalidSnark(format!(
             "get_poly_and_vec: Input vector length {} mismatch with num_vars {}", vec.len(), num_vars
         )));
    }
    let poly = MultilinearPolynomial::new(vec.clone()); // Clone needed as vec is consumed by new
    Ok((poly, vec))
}


/// The 'g' function for the simple range table [0..N-1].
/// g(E_1, ..., E_alpha) = Sum_{j=0}^{c-1} Sum_{i=0}^{k-1} E_{j*k+i} * 2^(l*j)
/// For k=1: g(E_0, ..., E_{c-1}) = E_0 + E_1*2^l + ... + E_{c-1}*2^(l*(c-1))
/// Note: Python example had a +1, assuming table [1..N]. This assumes [0..N-1]. Adjust if needed.
/// NOTE: This function is specific to the simple range table structure.
/// For general SOS tables, this needs to be replaced with the appropriate g function.
pub(crate) fn g_func_simple_range<F: PrimeField>(
    E_evals: &[F], // E_0(r'), ..., E_{alpha-1}(r')
    l: usize,
    c: usize,
    k: usize,
) -> F {
    let alpha = c * k;
    if E_evals.len() != alpha {
        // Consider returning a Result<F, Error> instead of panicking
        panic!("g_func_simple_range: Incorrect number of evaluations. Expected {}, got {}.", alpha, E_evals.len());
    }
    let mut result = F::ZERO;
    let mut power_of_2_l = F::ONE;
    let base = F::from(2u64).pow([l as u64, 0, 0, 0]); // 2^l

    for j in 0..c {
        let mut chunk_sum = F::ZERO;
        for i in 0..k {
            let index = j * k + i;
            if index < E_evals.len() {
                 chunk_sum += E_evals[index];
            }
            // else: Error? Or assume missing evals are 0? Lasso paper implies alpha evals.
        }
        result += chunk_sum * power_of_2_l;
        power_of_2_l *= base;
    }
    result
}

/// The 'g' function *polynomial* for the simple range table.
/// Input: E_poly[0], ..., E_poly[alpha-1]
/// Output: g(E_0, ..., E_{alpha-1}) as a polynomial
/// NOTE: This function is specific to the simple range table structure.
/// For general SOS tables, this needs to be replaced with the appropriate g polynomial constructor.
pub(crate) fn g_poly_simple_range<F: PrimeField>(
    E_polys: &[MultilinearPolynomial<F>],
    l: usize,
    c: usize,
    k: usize,
) -> Result<MultilinearPolynomial<F>, Error> {
    /* // Comment out body to avoid Add/Mul implementation for MultilinearPolynomial
    if E_polys.len() != c * k {
        return Err(Error::InvalidSnark("Incorrect number of E polys for g_poly".to_string()));
    }

    let mut result_poly = MultilinearPolynomial::zero(); // Assume all E have same num_vars (logm)
    let mut power_of_2_l = F::ONE;
    let base = F::from(2u64).pow([l as u64, 0, 0, 0]); // 2^l

    for j in 0..c {
        let mut chunk_poly = MultilinearPolynomial::zero();
        for i in 0..k {
            let index = j * k + i;
            // Ensure the index is within bounds
            if index < E_polys.len() {
                // If chunk_poly is the zero polynomial, initialize it with the first clone
                // Otherwise, add the cloned polynomial
                if chunk_poly.evals().is_empty() || chunk_poly.evals().iter().all(|v: &F| v.is_zero().into()) { // More robust check for zero poly
                     chunk_poly = E_polys[index].clone();
                } else {
                    chunk_poly = chunk_poly + E_polys[index].clone(); // Need to clone polys for addition
                }
            } // else: what to do if index is out of bounds? The initial check should prevent this.
        }
        result_poly = result_poly + chunk_poly * power_of_2_l;
        power_of_2_l *= base;
    }
    Ok(result_poly)
    */
    Err(Error::InvalidSnark("g_poly_simple_range is currently disabled".to_string())) // Return error instead
}

/// Placeholder hash function H described in Lasso paper section 3.1
/// H(x, y, z) = Hash(gamma, tau, x, y, z)
/// Needs a proper cryptographic hash implementation (e.g., using Poseidon or Keccak)
pub(crate) fn hash_tuple<F: PrimeField>(
    tuple: (F, F, F),
    gamma: &F, // Challenge, currently unused
    tau: &F,   // Challenge, currently unused
) -> F {
    // TODO: Replace with actual hash implementation (e.g. Keccak/Poseidon)
    // Simple non-cryptographic combination for testing purposes.
    tuple.0 + (*gamma * tuple.1) + (*tau * *gamma * tuple.2) // Example: Combine elements linearly
}