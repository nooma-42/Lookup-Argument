use std::collections::HashMap;
use crate::poly::multilinear::MultilinearPolynomial;
use crate::util::arithmetic::PrimeField;
use halo2_curves::bn256::Fr;

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

pub fn convert_to_logupgkr_format(
    lookup: Vec<Fr>,
    table: Vec<Fr>
) -> (HashMap<Vec<bool>, Fr>, HashMap<Vec<bool>, Fr>, Vec<HashMap<Vec<bool>, Fr>>) {
    if table.len() & (table.len() - 1) != 0 {
        panic!("table.len() must be a power of 2, you must pad it");
    }
    let num_bits = (table.len() as f64).log2().ceil() as usize;
    
    let mut t_values = HashMap::new();
    for (i, &value) in table.iter().enumerate() {
        let binary = index_to_binary(i, num_bits);
        t_values.insert(binary, value);
    }
    
    // Assume all lookup values are in the same column, and only one column
    let mut w_values_map = HashMap::new();
    for (i, &value) in lookup.iter().enumerate() {
        let binary = index_to_binary(i, num_bits);
        w_values_map.insert(binary, value);
    }
    let w_values = vec![w_values_map];
    
    let mut value_counts = HashMap::new();
    for &value in &lookup {
        *value_counts.entry(value).or_insert(0u64) += 1;
    }
    let mut m_values = HashMap::new();
    for (i, &value) in table.iter().enumerate() {
        let count = value_counts.get(&value).cloned().unwrap_or(0);
        let binary = index_to_binary(i, num_bits);
        m_values.insert(binary, Fr::from(count));
    }
    
    (m_values, t_values, w_values)
}

fn index_to_binary(index: usize, num_bits: usize) -> Vec<bool> {
    let mut binary = Vec::with_capacity(num_bits);
    let mut index = index;
    
    for _ in 0..num_bits {
        binary.push(index & 1 == 1);
        index >>= 1;
    }
    binary.reverse();
    
    binary
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2_curves::bn256::Fr;

    #[test]
    fn test_basic_conversion() {
        // Test case 1: Simple lookup with power-of-2 sized table
        let lookup = vec![Fr::one(), Fr::one(), Fr::from(2u64)];
        let table = vec![Fr::one(), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)];
        
        let (m_values, t_values, w_values) = convert_to_logupgkr_format(lookup, table);
        
        // Verify t_values (table)
        assert_eq!(t_values.get(&vec![false, false]), Some(&Fr::one()));
        assert_eq!(t_values.get(&vec![false, true]), Some(&Fr::from(2u64)));
        assert_eq!(t_values.get(&vec![true, false]), Some(&Fr::from(3u64)));
        assert_eq!(t_values.get(&vec![true, true]), Some(&Fr::from(4u64)));
        
        // Verify m_values (multiplicities)
        assert_eq!(m_values.get(&vec![false, false]), Some(&Fr::from(2u64))); // 1 appears twice
        assert_eq!(m_values.get(&vec![false, true]), Some(&Fr::from(1u64))); // 2 appears once
        assert_eq!(m_values.get(&vec![true, false]), Some(&Fr::from(0u64))); // 3 appears zero times
        assert_eq!(m_values.get(&vec![true, true]), Some(&Fr::from(0u64))); // 4 appears zero times
        
        // Verify w_values (witness)
        assert_eq!(w_values.len(), 1); // Single column
        let w_column = &w_values[0];
        assert_eq!(w_column.get(&vec![false, false]), Some(&Fr::one()));
        assert_eq!(w_column.get(&vec![false, true]), Some(&Fr::one()));
        assert_eq!(w_column.get(&vec![true, false]), Some(&Fr::from(2u64)));
    }

    #[test]
    #[should_panic(expected = "table.len() must be a power of 2")]
    fn test_non_power_of_two_table() {
        let lookup = vec![Fr::one(), Fr::one(), Fr::from(2u64)];
        let table = vec![Fr::one(), Fr::from(2u64), Fr::from(3u64)]; // Length 3 is not a power of 2
        
        convert_to_logupgkr_format(lookup, table);
    }

    #[test]
    fn test_empty_lookup() {
        let lookup: Vec<Fr> = vec![];
        let table = vec![Fr::one(), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)];
        
        let (m_values, t_values, w_values) = convert_to_logupgkr_format(lookup, table);
        
        // All multiplicities should be 0
        for value in m_values.values() {
            assert_eq!(*value, Fr::zero());
        }
        
        // w_values should have one empty column
        assert_eq!(w_values.len(), 1);
        assert!(w_values[0].is_empty());
    }

    #[test]
    fn test_lookup_with_repeated_values() {
        let lookup = vec![Fr::from(2u64), Fr::from(2u64), Fr::from(2u64), Fr::from(1u64)];
        let table = vec![Fr::one(), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)];
        
        let (m_values, t_values, w_values) = convert_to_logupgkr_format(lookup, table);
        
        // Check multiplicities
        assert_eq!(m_values.get(&vec![false, false]), Some(&Fr::from(1u64))); // 1 appears once
        assert_eq!(m_values.get(&vec![false, true]), Some(&Fr::from(3u64))); // 2 appears three times
        assert_eq!(m_values.get(&vec![true, false]), Some(&Fr::from(0u64))); // 3 appears zero times
        assert_eq!(m_values.get(&vec![true, true]), Some(&Fr::from(0u64))); // 4 appears zero times
    }
}