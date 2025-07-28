use std::time::{Duration, SystemTime};
use halo2_curves::ff::PrimeField;

/// Types of tests that can be performed on lookup argument systems
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TestType {
    /// Performance benchmarking (existing functionality)
    Performance,
    /// Soundness testing - invalid lookups should fail verification
    Soundness,
    /// Completeness testing - valid lookups should always pass verification
    Completeness,
}

impl std::fmt::Display for TestType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TestType::Performance => write!(f, "Performance"),
            TestType::Soundness => write!(f, "Soundness"),
            TestType::Completeness => write!(f, "Completeness"),
        }
    }
}

/// Result of a test execution
#[derive(Debug, Clone)]
pub struct TestResult {
    pub test_type: TestType,
    pub system_name: String,
    pub operation: String,
    pub k_value: usize,
    pub n_to_n_ratio: usize,
    pub passed: bool,
    pub setup_time: Option<Duration>,
    pub prove_time: Option<Duration>, 
    pub verify_time: Option<Duration>,
    pub proof_size: Option<u64>,
    pub error_message: Option<String>,
    pub timestamp: SystemTime,
}

impl TestResult {
    pub fn new(
        test_type: TestType,
        system_name: String,
        operation: String,
        k_value: usize,
        n_to_n_ratio: usize,
    ) -> Self {
        Self {
            test_type,
            system_name,
            operation,
            k_value,
            n_to_n_ratio,
            passed: false,
            setup_time: None,
            prove_time: None,
            verify_time: None,
            proof_size: None,
            error_message: None,
            timestamp: SystemTime::now(),
        }
    }

    pub fn success(mut self, setup_time: Duration, prove_time: Duration, verify_time: Duration, proof_size: u64) -> Self {
        self.passed = true;
        self.setup_time = Some(setup_time);
        self.prove_time = Some(prove_time);
        self.verify_time = Some(verify_time);
        self.proof_size = Some(proof_size);
        self
    }

    pub fn failure(mut self, error_message: String) -> Self {
        self.passed = false;
        self.error_message = Some(error_message);
        self
    }

    pub fn total_time(&self) -> Option<Duration> {
        match (self.setup_time, self.prove_time, self.verify_time) {
            (Some(setup), Some(prove), Some(verify)) => Some(setup + prove + verify),
            _ => None,
        }
    }
}

/// Test case generator for soundness and completeness testing
pub struct TestCaseGenerator;

impl TestCaseGenerator {
    /// Generate valid values for completeness testing
    pub fn generate_valid_range_values<F: PrimeField>(k: usize, n_to_n_ratio: usize) -> Vec<F> {
        let range_size = 1 << k;
        let lookup_size = range_size / n_to_n_ratio;
        
        // Generate all valid values within the range
        (0..lookup_size)
            .map(|i| F::from((i % range_size) as u64))
            .collect()
    }

    /// Generate invalid values for soundness testing
    pub fn generate_invalid_range_values<F: PrimeField>(k: usize, n_to_n_ratio: usize) -> Vec<F> {
        let range_size = 1 << k;
        let lookup_size = range_size / n_to_n_ratio;
        
        // Generate values outside the valid range [0, 2^k-1]
        let mut invalid_values = Vec::new();
        
        // Add values that are too large
        for i in 0..lookup_size {
            let invalid_value = range_size + i + 1; // Values >= 2^k
            invalid_values.push(F::from(invalid_value as u64));
        }
        
        invalid_values
    }

    /// Generate valid addition test cases for completeness testing
    pub fn generate_valid_add_values<F: PrimeField>(k: usize, n_to_n_ratio: usize) -> Vec<(F, F, F)> {
        let range_size = 1 << k;
        let lookup_size = range_size / n_to_n_ratio;
        
        let mut valid_cases = Vec::new();
        
        for i in 0..lookup_size {
            let a = (i / 2) % (range_size / 2); // Ensure a + b doesn't overflow
            let b = i % (range_size / 2);
            let sum = (a + b) % range_size; // Stay within range
            
            valid_cases.push((
                F::from(a as u64),
                F::from(b as u64), 
                F::from(sum as u64)
            ));
        }
        
        valid_cases
    }

    /// Generate invalid addition test cases for soundness testing
    pub fn generate_invalid_add_values<F: PrimeField>(k: usize, n_to_n_ratio: usize) -> Vec<(F, F, F)> {
        let range_size = 1 << k;
        let lookup_size = range_size / n_to_n_ratio;
        
        let mut invalid_cases = Vec::new();
        
        for i in 0..lookup_size {
            let a = i % (range_size / 2);
            let b = (i + 1) % (range_size / 2);
            let correct_sum = (a + b) % range_size;
            let wrong_sum = (correct_sum + 1 + i) % range_size; // Intentionally wrong result
            
            invalid_cases.push((
                F::from(a as u64),
                F::from(b as u64),
                F::from(wrong_sum as u64)
            ));
        }
        
        invalid_cases
    }
}

/// Test statistics collector
#[derive(Debug, Default)]
pub struct TestStatistics {
    pub total_tests: usize,
    pub passed_tests: usize,
    pub failed_tests: usize,
    pub performance_tests: usize,
    pub soundness_tests: usize,
    pub completeness_tests: usize,
}

impl TestStatistics {
    pub fn add_result(&mut self, result: &TestResult) {
        self.total_tests += 1;
        
        if result.passed {
            self.passed_tests += 1;
        } else {
            self.failed_tests += 1;
        }
        
        match result.test_type {
            TestType::Performance => self.performance_tests += 1,
            TestType::Soundness => self.soundness_tests += 1,
            TestType::Completeness => self.completeness_tests += 1,
        }
    }

    pub fn success_rate(&self) -> f64 {
        if self.total_tests == 0 {
            0.0
        } else {
            self.passed_tests as f64 / self.total_tests as f64
        }
    }

    pub fn soundness_pass_rate(&self) -> f64 {
        if self.soundness_tests == 0 {
            0.0
        } else {
            // For soundness tests, we expect them to FAIL (detect invalid inputs)
            // So a "passed" soundness test means it correctly rejected invalid input
            self.passed_tests as f64 / self.soundness_tests as f64
        }
    }

    pub fn completeness_pass_rate(&self) -> f64 {
        if self.completeness_tests == 0 {
            0.0
        } else {
            // For completeness tests, we expect them to PASS (accept valid inputs)
            self.passed_tests as f64 / self.completeness_tests as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2_curves::bn256::Fr;

    #[test]
    fn test_valid_range_values_generation() {
        let k = 4; // Range [0, 15]
        let n_to_n_ratio = 2;
        
        let values = TestCaseGenerator::generate_valid_range_values::<Fr>(k, n_to_n_ratio);
        
        assert_eq!(values.len(), 8); // 16/2 = 8
        
        // All values should be in valid range [0, 15]
        for value in values {
            let val_u64 = u64::from(value.to_repr().as_ref()[0]);
            assert!(val_u64 < 16);
        }
    }

    #[test]
    fn test_invalid_range_values_generation() {
        let k = 3; // Range [0, 7]  
        let n_to_n_ratio = 4;
        
        let values = TestCaseGenerator::generate_invalid_range_values::<Fr>(k, n_to_n_ratio);
        
        assert_eq!(values.len(), 2); // 8/4 = 2
        
        // All values should be outside valid range [0, 7]
        for value in values {
            let val_u64 = u64::from(value.to_repr().as_ref()[0]);
            assert!(val_u64 >= 8);
        }
    }

    #[test]
    fn test_valid_add_values_generation() {
        let k = 4; // Values in range [0, 15]
        let n_to_n_ratio = 8;
        
        let cases = TestCaseGenerator::generate_valid_add_values::<Fr>(k, n_to_n_ratio);
        
        assert_eq!(cases.len(), 2); // 16/8 = 2
        
        // Verify all cases have valid additions
        for (a, b, sum) in cases {
            let a_u64 = u64::from(a.to_repr().as_ref()[0]);
            let b_u64 = u64::from(b.to_repr().as_ref()[0]);
            let sum_u64 = u64::from(sum.to_repr().as_ref()[0]);
            
            assert_eq!((a_u64 + b_u64) % 16, sum_u64);
        }
    }

    #[test]
    fn test_statistics() {
        let mut stats = TestStatistics::default();
        
        let result1 = TestResult::new(TestType::Performance, "Lasso".to_string(), "range".to_string(), 8, 2)
            .success(Duration::from_millis(10), Duration::from_millis(20), Duration::from_millis(5), 1024);
        
        let result2 = TestResult::new(TestType::Soundness, "Lasso".to_string(), "range".to_string(), 8, 2)
            .failure("Correctly rejected invalid input".to_string());
        
        stats.add_result(&result1);
        stats.add_result(&result2);
        
        assert_eq!(stats.total_tests, 2);
        assert_eq!(stats.passed_tests, 1);
        assert_eq!(stats.failed_tests, 1);
        assert_eq!(stats.performance_tests, 1);
        assert_eq!(stats.soundness_tests, 1);
    }
}