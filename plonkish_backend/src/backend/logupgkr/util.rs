use std::collections::HashMap;
use crate::poly::multilinear::MultilinearPolynomial;
use crate::util::arithmetic::PrimeField;
use halo2_curves::bn256::Fr;
use itertools::izip;

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
        let x_field: Vec<F> = x.iter().map(|&b| if b { F::ONE } else { F::ZERO }).collect();
        m.evaluate(&x_field)
    } else {
        -F::ONE
    }
}

pub fn q<F: PrimeField>(x: &[bool], y: &[bool], t: &MultilinearPolynomial<F>, w: &[MultilinearPolynomial<F>], a: F) -> F {
    let y_index = binary_to_usize(y);
    let total_cols = w.len() + 1; // +1 for the t column
    let max_index = if total_cols.is_power_of_two() {
        total_cols - 1
    } else {
        (1 << ((total_cols.ilog2() + 1) as usize)) - 1
    };
    
    let x_field: Vec<F> = x.iter().map(|&b| if b { F::ONE } else { F::ZERO }).collect();
    
    if y_index == max_index {
        // The highest index corresponds to the table polynomial
        a - t.evaluate(&x_field)
    } else if y_index < w.len() {
        // Valid witness polynomial index
        a - w[y_index].evaluate(&x_field)
    } else {
        // Padding case - treat as zero contribution
        a
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

/// Verifier calculates the expected value of p(x), using the multilinear extension formula.
///
/// # Arguments
/// * `x_full` - the full challenge point `(x_row, x_col)` returned by GKR protocol.
/// * `m_at_x_row` - the trusted `m(x_row)` evaluation.
/// * `num_vars_row` - the number of variables `n`.
///
/// # Returns
/// The expected evaluation of p(x_full).
pub fn evaluate_p_at_x<F: PrimeField>(
    x_full: &[F],
    m_at_x_row: F,
    num_vars_row: usize,
) -> F {
    let x_col = &x_full[num_vars_row..];

    // L_k(x_col, 1) multilinear extension ∏ x_col_i
    let one = F::ONE;
    let l_k_at_1 = x_col.iter().copied().product::<F>();

    // p(x) = L_k(x_col, 1) * m(x_row) + (1 - L_k(x_col, 1)) * (-1)
    let one_minus_lk = one - l_k_at_1;
    l_k_at_1 * m_at_x_row - one_minus_lk
}

/// Verifier calculates the expected value of q(x), using the multilinear extension formula.
///
/// # Arguments
/// * `x_full` - the full challenge point `(x_row, x_col)` returned by GKR protocol.
/// * `t_at_x_row` - the trusted `t(x_row)` evaluation.
/// * `w_at_x_rows` - the trusted `w_i(x_row)` evaluation list.
/// * `a` - the random challenge obtained from transcript.
/// * `num_vars_row` - the number of variables `n`.
///
/// # Returns
/// The expected evaluation of q(x_full).
pub fn evaluate_q_at_x<F: PrimeField>(
    x_full: &[F],
    t_at_x_row: F,
    w_at_x_rows: &[F],
    a: F,
    num_vars_row: usize,
) -> F {
    let x_col = &x_full[num_vars_row..];
    let num_vars_col = x_col.len();
    
    // pre-calculate all y_j ∈ H_k\{1} corresponding to L_k(x_col, y_j)
    // H_k\{1} corresponds to index 0 to 2^k - 2
    // assume w_0, w_1, ... corresponds to y_0, y_1, ...
    // the boolean representation of H_k is {0,1}^k
    let h_k_bool = generate_binary_combinations(num_vars_col as u32);

    let one = F::ONE;
    let two_inv = F::from(2).invert().unwrap();

    // calculate the general function of L_k(x_col, y)
    let evaluate_lagrange_basis = |y_bool: &[bool]| -> F {
        assert_eq!(y_bool.len(), num_vars_col);
        izip!(x_col.iter(), y_bool.iter())
            .map(|(challenge_i, y_i_bool)| {
                // convert boolean {0, 1} to field element {0, 1}
                let y_i_field = if *y_i_bool { one } else { F::ZERO };
                (one - *challenge_i) * (one - y_i_field) + *challenge_i * y_i_field
            })
            .product::<F>()
    };

    // calculate the second part of q(x): ∑ L_k(x_col, y_j) * (a - w_j(x_row))
    // assume the order of `w_at_x_rows` matches the order of non-all-true vectors in `h_k_bool`
    // the order of `h_k_bool` is 00..0, 10..0, 01..0, ...
    // we need a deterministic mapping i(y) -> witness_index
    // for simplicity, we assume the binary value of y is the index of witness
    // M = 2^k - 1, the length of w_polys is M
    // we need to be careful about the index mapping
    
    // Updated logic to match the q function
    let total_cols = w_at_x_rows.len() + 1; // +1 for the t column
    let max_index = if total_cols.is_power_of_two() {
        total_cols - 1
    } else {
        (1 << num_vars_col) - 1
    };

    // Use the binary indexing to correctly map y vectors to indices
    let mut q_eval = F::ZERO;
    
    for y_bool in h_k_bool.iter() {
        let y_index = binary_to_usize(y_bool);
        let l_k_at_y = evaluate_lagrange_basis(y_bool);
        
        if y_index == max_index {
            // The highest index corresponds to the table polynomial
            q_eval += l_k_at_y * (a - t_at_x_row);
        } else if y_index < w_at_x_rows.len() {
            // Valid witness polynomial index
            q_eval += l_k_at_y * (a - w_at_x_rows[y_index]);
        } else {
            // Padding case - treat as zero contribution (a - a = 0), so add a
            q_eval += l_k_at_y * a;
        }
    }

    q_eval
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