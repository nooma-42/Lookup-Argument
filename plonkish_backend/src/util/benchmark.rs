use halo2_curves::bn256::Fr;
use halo2_curves::ff::PrimeField;

/// Generate values for range check benchmarking
/// 
/// This function creates test data for range check operations in lookup arguments.
/// - `k`: Log of the range size (range will be 2^k values: 0 to 2^k-1)  
/// - `n_to_n_ratio`: The ratio N:n where N is table size and n is lookup size
/// 
/// Returns a vector of Fr elements representing values to be looked up in the range [0, 2^k-1]
pub fn generate_range_check_values(k: usize, n_to_n_ratio: usize) -> Vec<Fr> {
    let range_size = 1 << k; // 2^k - full range without artificial restrictions
    let lookup_size = range_size / n_to_n_ratio;
    
    // Generate values to check (within the range)
    (0..lookup_size)
        .map(|i| Fr::from((i % range_size) as u64))
        .collect()
}

/// Generate invalid values for soundness testing
/// 
/// This function creates test data that should be rejected by range check operations.
/// - `k`: Log of the range size (valid range is [0, 2^k-1])
/// - `n_to_n_ratio`: The ratio N:n where N is table size and n is lookup size
/// 
/// Returns a vector of Fr elements representing values outside the valid range
pub fn generate_invalid_values(k: usize, n_to_n_ratio: usize) -> Vec<Fr> {
    let range_size = 1 << k;
    let lookup_size = range_size / n_to_n_ratio;
    
    // Generate values outside the valid range [0, 2^k-1]
    (0..lookup_size)
        .map(|i| Fr::from((range_size + i + 1) as u64)) // All values >= 2^k
        .collect()
}

/// Generate test data for addition operations
/// 
/// This function creates test data for addition lookup tables used in Lasso.
/// - `k`: Log of the range size for operands
/// - `n_to_n_ratio`: The ratio N:n where N is table size and n is lookup size
/// 
/// Returns a vector of (a, b, sum) tuples where all values are in range [0, 2^k-1]
pub fn generate_add_operation_data(k: usize, n_to_n_ratio: usize) -> Vec<(Fr, Fr, Fr)> {
    let range_size = 1 << k;
    let lookup_size = range_size / n_to_n_ratio;
    
    let mut add_cases = Vec::new();
    
    for i in 0..lookup_size {
        // Generate operands that won't overflow the range
        let a = (i / 2) % (range_size / 2);
        let b = i % (range_size / 2);
        let sum = a + b;  // For completeness test, sum should be the actual sum
        
        // Skip if sum exceeds the range (would be invalid for the table)
        if sum >= range_size {
            continue;
        }
        
        add_cases.push((
            Fr::from(a as u64),
            Fr::from(b as u64),
            Fr::from(sum as u64)
        ));
    }
    
    add_cases
}

/// Generate invalid test data for addition operations (soundness testing)
/// 
/// This function creates test data with incorrect addition results that should be rejected.
/// - `k`: Log of the range size for operands
/// - `n_to_n_ratio`: The ratio N:n where N is table size and n is lookup size
/// 
/// Returns a vector of (a, b, wrong_sum) tuples where wrong_sum != a + b
pub fn generate_invalid_add_operation_data(k: usize, n_to_n_ratio: usize) -> Vec<(Fr, Fr, Fr)> {
    let range_size = 1 << k;
    let lookup_size = range_size / n_to_n_ratio;
    
    let mut invalid_cases = Vec::new();
    
    for i in 0..lookup_size {
        let a = i % (range_size / 2);
        let b = (i + 1) % (range_size / 2);
        let correct_sum = (a + b) % range_size;
        let wrong_sum = (correct_sum + 1 + i) % range_size; // Intentionally wrong
        
        invalid_cases.push((
            Fr::from(a as u64),
            Fr::from(b as u64),
            Fr::from(wrong_sum as u64)
        ));
    }
    
    invalid_cases
}

/// Generate a table for range check benchmarking
/// 
/// This function creates the lookup table containing all values in the range [0, 2^k-1]
/// - `k`: Log of the range size
/// 
/// Returns a vector of Fr elements representing the complete range table
pub fn generate_range_check_table(k: usize) -> Vec<Fr> {
    let range_size = 1 << k; // 2^k
    
    (0..range_size)
        .map(|i| Fr::from(i as u64))
        .collect()
}

/// Generate both table and lookup values for range check benchmarking
/// 
/// This is a convenience function that generates both the complete range table
/// and the lookup values for benchmarking.
/// - `k`: Log of the range size
/// - `n_to_n_ratio`: The ratio N:n where N is table size and n is lookup size
/// 
/// Returns a tuple (table, lookup_values) where:
/// - table: Complete range [0, 2^k-1] 
/// - lookup_values: Values to be looked up, sized according to ratio
pub fn generate_range_check_data(k: usize, n_to_n_ratio: usize) -> (Vec<Fr>, Vec<Fr>) {
    let table = generate_range_check_table(k);
    let lookup_values = generate_range_check_values(k, n_to_n_ratio);
    (table, lookup_values)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test] 
    fn test_generate_range_check_values() {
        let k = 4; // 16 values in range [0, 15]
        let n_to_n_ratio = 2;
        
        let values = generate_range_check_values(k, n_to_n_ratio);
        
        // Should have 16/2 = 8 lookup values
        assert_eq!(values.len(), 8);
        
        // All values should be in range [0, 15]
        for value in values {
            let val_u64 = u64::from(value.to_repr().as_ref()[0]);
            assert!(val_u64 < 16);
        }
    }

    #[test]
    fn test_generate_invalid_values() {
        let k = 3; // Valid range [0, 7]
        let n_to_n_ratio = 4;
        
        let values = generate_invalid_values(k, n_to_n_ratio);
        
        // Should have 8/4 = 2 invalid values
        assert_eq!(values.len(), 2);
        
        // All values should be outside valid range [0, 7]
        for value in values {
            let val_u64 = u64::from(value.to_repr().as_ref()[0]);
            assert!(val_u64 >= 8); // Should be >= 2^k
        }
    }

    #[test]
    fn test_generate_add_operation_data() {
        let k = 4; // Range [0, 15]
        let n_to_n_ratio = 8;
        
        let cases = generate_add_operation_data(k, n_to_n_ratio);
        
        // Should have 16/8 = 2 cases
        assert_eq!(cases.len(), 2);
        
        // Verify all cases have valid additions
        for (a, b, sum) in cases {
            let a_u64 = u64::from(a.to_repr().as_ref()[0]);
            let b_u64 = u64::from(b.to_repr().as_ref()[0]);
            let sum_u64 = u64::from(sum.to_repr().as_ref()[0]);
            
            // All values should be in range [0, 15]
            assert!(a_u64 < 16);
            assert!(b_u64 < 16);
            assert!(sum_u64 < 16);
            
            // Sum should be correct (modulo 16)
            assert_eq!((a_u64 + b_u64) % 16, sum_u64);
        }
    }

    #[test]
    fn test_generate_invalid_add_operation_data() {
        let k = 3; // Range [0, 7]
        let n_to_n_ratio = 4;
        
        let cases = generate_invalid_add_operation_data(k, n_to_n_ratio);
        
        // Should have 8/4 = 2 invalid cases
        assert_eq!(cases.len(), 2);
        
        // Verify all cases have incorrect additions
        for (a, b, wrong_sum) in cases {
            let a_u64 = u64::from(a.to_repr().as_ref()[0]);
            let b_u64 = u64::from(b.to_repr().as_ref()[0]);
            let wrong_sum_u64 = u64::from(wrong_sum.to_repr().as_ref()[0]);
            
            // All values should be in range [0, 7]
            assert!(a_u64 < 8);
            assert!(b_u64 < 8);
            assert!(wrong_sum_u64 < 8);
            
            // Sum should be incorrect (not equal to a + b mod 8)
            assert_ne!((a_u64 + b_u64) % 8, wrong_sum_u64);
        }
    }
    
    #[test]
    fn test_generate_range_check_table() {
        let k = 3; // 8 values in range [0, 7]
        
        let table = generate_range_check_table(k);
        
        // Should have exactly 2^3 = 8 values
        assert_eq!(table.len(), 8);
        
        // Should contain values 0 through 7
        for (i, value) in table.iter().enumerate() {
            let val_u64 = u64::from(value.to_repr().as_ref()[0]);
            assert_eq!(val_u64, i as u64);
        }
    }
    
    #[test]
    fn test_generate_range_check_data() {
        let k = 3; // 8 values in range [0, 7] 
        let n_to_n_ratio = 4;
        
        let (table, lookup_values) = generate_range_check_data(k, n_to_n_ratio);
        
        // Table should have 2^3 = 8 values
        assert_eq!(table.len(), 8);
        
        // Lookup should have 8/4 = 2 values
        assert_eq!(lookup_values.len(), 2);
    }
}