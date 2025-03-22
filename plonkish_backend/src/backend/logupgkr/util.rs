use std::collections::HashMap;
use crate::poly::multilinear::MultilinearPolynomial;
use crate::util::arithmetic::PrimeField;

// Helper function to generate all binary combinations of length n
pub fn generate_binary_combinations(length: u32) -> Vec<Vec<bool>> {
    let total_combinations = 1 << length;
    let mut combinations = Vec::with_capacity(total_combinations);

    for i in 0..total_combinations {
        let mut combination = Vec::with_capacity(length as usize);
        for j in 0..length {
            combination.push((i & (1 << j)) != 0);
        }
        combinations.push(combination);
    }

    combinations
}

pub fn binary_to_usize(binary: &[bool]) -> usize {
    binary.iter().fold(0, |acc, &b| (acc << 1) | (b as usize))
}

pub fn p<F: PrimeField>(x: &[bool], y: &[bool], m: &MultilinearPolynomial<F>) -> F {
    if y.iter().all(|&value| value) {
        let x_field: Vec<F> = x.iter().map(|&b| F::from(b as u64)).collect();
        m.evaluate(&x_field)
    } else {
        -F::from(1u64)
    }
}

pub fn q<F: PrimeField>(x: &[bool], y: &[bool], t: &MultilinearPolynomial<F>, w: &[MultilinearPolynomial<F>], a: F) -> F {
    if y.iter().all(|&value| value) {
        let x_field: Vec<F> = x.iter().map(|&b| F::from(b as u64)).collect();
        a - t.evaluate(&x_field)
    } else {
        let y_index = binary_to_usize(y);
        let x_field: Vec<F> = x.iter().map(|&b| F::from(b as u64)).collect();
        a - w[y_index].evaluate(&x_field)
    }
}

pub fn create_multilinear_poly<F: PrimeField>(map: HashMap<Vec<bool>, F>) -> MultilinearPolynomial<F> {
    let num_vars = map.keys().next().unwrap().len();
    let mut evals = vec![F::from(0u64); 1 << num_vars];
    for (input, value) in map {
        let index = binary_to_usize(&input);
        evals[index] = value;
    }
    MultilinearPolynomial::new(evals)
}