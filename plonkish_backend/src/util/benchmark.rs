use halo2_curves::bn256::Fr;

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
            let val_u64: u64 = value.into();
            assert!(val_u64 < 16);
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
            let val_u64: u64 = (*value).into();
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