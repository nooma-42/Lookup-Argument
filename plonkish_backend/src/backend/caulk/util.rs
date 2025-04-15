use crate::poly::univariate::UnivariatePolynomial;
use crate::util::arithmetic::{root_of_unity, Field, PrimeField};
use std::collections::HashSet;

pub fn is_power_of_2(n: usize) -> bool {
    n != 0 && (n & (n - 1)) == 0
}

pub fn get_roots<F: Field + PrimeField>(N: usize) -> Vec<F> {
    assert!(is_power_of_2(N));
    let root = root_of_unity::<F>(N.ilog2() as usize);
    (0..N).map(|i| root.pow(&[i as u64])).collect()
}

pub fn get_vanishing_poly<F: PrimeField>(n: usize) -> UnivariatePolynomial<F> {
    let mut coeffs = vec![F::ZERO; n + 1];
    coeffs[0] = -F::ONE;
    coeffs[n] = F::ONE;
    UnivariatePolynomial::monomial(coeffs)
}

pub fn get_unique_positions<F: PrimeField>(c: &[F], values: &[F]) -> Vec<usize> {
    let mut result = HashSet::new();
    for value in values {
        // This will panic if a value in `values` is not present in `c`.
        // Consider returning Result or handling the error appropriately if this is possible.
        result.insert(c.iter().position(|x| x == value).unwrap());
    }
    result.into_iter().collect()
} 