use std::fmt::Debug;

use halo2_curves::bn256::{Bn256, Fr, G1Affine, G2Affine, G1};
use halo2_curves::group::{Group, Curve};

use crate::{
    backend::cq::generate_table_and_lookup,
    pcs::{
        univariate::{
            UnivariateKzg, UnivariateKzgCommitment, UnivariateKzgParam, UnivariateKzgProverParam,
            UnivariateKzgVerifierParam,
        }, PolynomialCommitmentScheme,
    },
    poly::univariate::UnivariatePolynomial,
    poly::Polynomial,
    util::{
        arithmetic::Field,
        transcript::InMemoryTranscript,
    },
    Error,
};

pub mod key;
pub mod preprocessor;
pub mod prover;
pub mod util;
pub mod verifier;

// Specific implementation for Bn256 curves
#[derive(Clone, Debug)]
pub struct BalooProverParam {
    param: UnivariateKzgParam<Bn256>,
    pp: UnivariateKzgProverParam<Bn256>,
    table: Vec<Fr>,
    d: usize,
}

#[derive(Clone, Debug)]
pub struct BalooVerifierParam {
    vp: UnivariateKzgVerifierParam<Bn256>,
}

#[derive(Clone, Debug)]
pub struct BalooInfo {
    table: Vec<Fr>,
    lookup: Vec<Fr>,
}

#[derive(Clone, Debug)]
pub struct Baloo;

impl Baloo {
    // Main preprocess function that matches the original API
    pub fn preprocess(
        t: usize,
        m: usize,
    ) -> Result<
        (
            UnivariateKzgParam<Bn256>,
            UnivariateKzgProverParam<Bn256>,
            UnivariateKzgVerifierParam<Bn256>,
        ),
        Error,
    > {
        preprocessor::preprocess(t, m)
    }

    // Optimized table-specific preprocessing that achieves O(m) proving complexity
    pub fn preprocess_for_table(
        table: &[Fr],
        m_max: usize,
    ) -> Result<(key::BalooProverKey, key::BalooVerifierKey), Error> {
        preprocessor::preprocess_for_table(table, m_max)
    }

    // Alternative preprocess function that takes an info struct
    pub fn preprocess_with_info(
        info: &BalooInfo,
    ) -> Result<(BalooProverParam, BalooVerifierParam), Error> {
        let m = info.lookup.len();
        let t = info.table.len();

        let (param, pp, vp) = preprocessor::preprocess(t, m)?;

        let poly_size = std::cmp::max(t, m).next_power_of_two() * 2;
        let d = poly_size - 2;

        Ok((
            BalooProverParam {
                param: param.clone(),
                pp,
                table: info.table.clone(),
                d,
            },
            BalooVerifierParam { vp },
        ))
    }

    pub fn prove(
        table: &Vec<Fr>,
        param: &UnivariateKzgParam<Bn256>,
        pp: &UnivariateKzgProverParam<Bn256>,
        lookup: &Vec<Fr>,
    ) -> Vec<u8> {
        // Legacy prove function - uses the original Prover interface for backward compatibility
        // For optimal O(m) performance, use prove_with_key instead
        let prover = prover::Prover::new(table, param, pp);
        prover.prove(lookup)
    }

    // Optimized prove function using preprocessed BalooProverKey (O(m) complexity)
    pub fn prove_with_key(pk: &key::BalooProverKey, lookup: &[Fr]) -> Vec<u8> {
        let prover = prover::OptimizedProver::new(pk);
        prover.prove(lookup)
    }

    pub fn verify_with_key(
        vk: &key::BalooVerifierKey,
        proof: &[u8],
        phi_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_m_exponent_poly_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_comm_2: &UnivariateKzgCommitment<G2Affine>,
        x_exponent_poly_2_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_2_comm_2: &UnivariateKzgCommitment<G2Affine>,
        m: usize,
    ) -> bool {
        let verifier = verifier::Verifier::new(vk);
        verifier.verify(
            &proof.to_vec(),
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
            m,
        )
    }

    pub fn prove_with_param(pp: &BalooProverParam, lookup: &Vec<Fr>) -> Vec<u8> {
        // Convert BalooProverParam to legacy API call
        Baloo::prove(&pp.table, &pp.param, &pp.pp, lookup)
    }

    pub fn verify(
        vp: &UnivariateKzgVerifierParam<Bn256>,
        proof: &Vec<u8>,
        t_comm_1: &UnivariateKzgCommitment<G1Affine>,
        z_h_comm_1: &UnivariateKzgCommitment<G1Affine>,
        phi_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_m_exponent_poly_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_comm_2: &UnivariateKzgCommitment<G2Affine>,
        x_exponent_poly_2_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_2_comm_2: &UnivariateKzgCommitment<G2Affine>,
        m: usize,
    ) -> bool {
        let verifier = verifier::LegacyVerifier::new(vp);
        verifier.verify(
            proof,
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
            m,
        )
    }

    pub fn verify_with_param(
        vp: &BalooVerifierParam,
        proof: &Vec<u8>,
        t_comm_1: &UnivariateKzgCommitment<G1Affine>,
        z_h_comm_1: &UnivariateKzgCommitment<G1Affine>,
        phi_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_m_exponent_poly_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_comm_2: &UnivariateKzgCommitment<G2Affine>,
        x_exponent_poly_2_comm_1: &UnivariateKzgCommitment<G1Affine>,
        x_exponent_poly_2_comm_2: &UnivariateKzgCommitment<G2Affine>,
        m: usize,
    ) -> bool {
        let verifier = verifier::LegacyVerifier::new(&vp.vp);
        verifier.verify(
            proof,
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
            m,
        )
    }

    // Helper method to prepare verification data
    pub fn prepare_verification_data(
        param: &UnivariateKzgParam<Bn256>,
        pp: &UnivariateKzgProverParam<Bn256>,
        table: &[Fr],
        lookup: &[Fr],
        m: usize,
        t: usize,
        d: usize,
    ) -> (
        UnivariateKzgCommitment<G1Affine>, // t_comm_1
        UnivariateKzgCommitment<G1Affine>, // z_h_comm_1
        UnivariateKzgCommitment<G1Affine>, // phi_comm_1
        UnivariateKzgCommitment<G1Affine>, // x_m_exponent_poly_comm_1
        UnivariateKzgCommitment<G2Affine>, // x_exponent_poly_comm_2
        UnivariateKzgCommitment<G1Affine>, // x_exponent_poly_2_comm_1
        UnivariateKzgCommitment<G2Affine>, // x_exponent_poly_2_comm_2
    ) {
        let scalar_0 = Fr::zero();
        let scalar_1 = Fr::one();

        let z_h_poly_coeffs = vec![scalar_1.neg()]
            .into_iter()
            .chain(vec![scalar_0; t - 1])
            .chain(vec![scalar_1])
            .collect();
        let z_h_poly = UnivariatePolynomial::monomial(z_h_poly_coeffs);
        let z_h_comm_1 = UnivariateKzg::<Bn256>::commit_monomial(pp, z_h_poly.coeffs());

        let t_poly = UnivariatePolynomial::lagrange(table.to_vec()).ifft();
        let t_comm_1 = UnivariateKzg::<Bn256>::commit_monomial(pp, t_poly.coeffs());

        let phi_poly = UnivariatePolynomial::lagrange(lookup.to_vec()).ifft();
        let phi_comm_1 = UnivariateKzg::<Bn256>::commit_monomial(pp, phi_poly.coeffs());

        // X^m
        let x_m_exponent_poly = UnivariatePolynomial::monomial(
            vec![scalar_0; m]
                .into_iter()
                .chain(vec![scalar_1])
                .collect(),
        );
        let x_m_exponent_poly_comm_1 =
            UnivariateKzg::<Bn256>::commit_monomial(pp, x_m_exponent_poly.clone().coeffs());

        // X^(d-m+1)
        let coeffs_x_exponent_poly = vec![scalar_0; d - m + 1]
            .into_iter()
            .chain(vec![scalar_1])
            .collect();
        let x_exponent_poly = UnivariatePolynomial::monomial(coeffs_x_exponent_poly);
        // [X^(d-m+1)]2
        let x_exponent_poly_comm_2 =
            UnivariateKzg::<Bn256>::commit_monomial_g2(param, x_exponent_poly.coeffs());

        // X^(d-m+2)
        let coeffs_x_exponent_poly_2 = vec![scalar_0; d - m + 2]
            .into_iter()
            .chain(vec![scalar_1])
            .collect();
        let x_exponent_poly_2 = UnivariatePolynomial::monomial(coeffs_x_exponent_poly_2);
        let x_exponent_poly_2_comm_1 =
            UnivariateKzg::<Bn256>::commit_monomial(pp, x_exponent_poly_2.coeffs());
        let x_exponent_poly_2_comm_2 =
            UnivariateKzg::<Bn256>::commit_monomial_g2(param, x_exponent_poly_2.coeffs());

        (
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
        )
    }

    // Run the full Baloo protocol with given table and lookup
    pub fn test_baloo_by_input(table: Vec<Fr>, lookup: Vec<Fr>) -> Vec<String> {
        let mut timings: Vec<String> = vec![];

        let start_total = std::time::Instant::now();

        let m = lookup.len();
        let t = table.len();
        let poly_size = std::cmp::max(t, m).next_power_of_two() * 2;
        let d = poly_size - 2;

        // 1. Setup
        let start = std::time::Instant::now();
        let (param, pp, vp) = Baloo::preprocess(t, m).unwrap();
        let duration1 = start.elapsed();
        timings.push(format!("Setup and preprocess: {}ms", duration1.as_millis()));

        // 2. Generate proof
        let start = std::time::Instant::now();
        let proof = Baloo::prove(&table, &param, &pp, &lookup);
        let duration2 = start.elapsed();
        timings.push(format!("Prove: {}ms", duration2.as_millis()));
        timings.push(format!("Proof size: {} bytes", proof.len()));

        // 3. Prepare verification data
        let start = std::time::Instant::now();
        let (
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
        ) = Baloo::prepare_verification_data(&param, &pp, &table, &lookup, m, t, d);

        // 4. Verify
        let result = Baloo::verify(
            &vp,
            &proof,
            &t_comm_1,
            &z_h_comm_1,
            &phi_comm_1,
            &x_m_exponent_poly_comm_1,
            &x_exponent_poly_comm_2,
            &x_exponent_poly_2_comm_1,
            &x_exponent_poly_2_comm_2,
            m,
        );

        assert!(result);
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));

        timings
    }

    // Run the optimized Baloo protocol with O(m) proving complexity
    pub fn test_baloo_optimized_by_input(table: Vec<Fr>, lookup: Vec<Fr>) -> Vec<String> {
        let mut timings: Vec<String> = vec![];
        
        let start_total = std::time::Instant::now();
        
        let m = lookup.len();
        let t = table.len();
        
        // 1. Optimized preprocessing with table-specific precomputation
        let start = std::time::Instant::now();
        let (pk, vk) = Baloo::preprocess_for_table(&table, m).unwrap();
        let duration1 = start.elapsed();
        timings.push(format!("Optimized Setup+Preprocess: {}ms", duration1.as_millis()));
        
        // 2. Generate proof using optimized prover (O(m) complexity)
        let start = std::time::Instant::now();
        let proof = Baloo::prove_with_key(&pk, &lookup);
        let duration2 = start.elapsed();
        timings.push(format!("Optimized Prove: {}ms", duration2.as_millis()));
        timings.push(format!("Proof size: {} bytes", proof.len()));
        
        // 3. Prepare verification data
        let start = std::time::Instant::now();
        let poly_size = std::cmp::max(t, m).next_power_of_two() * 2;
        let d = poly_size - 2;
        let (
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
        ) = Baloo::prepare_verification_data(&pk.param, &pk.pp, &table, &lookup, m, t, d);
        
        // 4. Verify using optimized verifier
        let result = Baloo::verify_with_key(
            &vk,
            &proof,
            &phi_comm_1,
            &x_m_exponent_poly_comm_1,
            &x_exponent_poly_comm_2,
            &x_exponent_poly_2_comm_1,
            &x_exponent_poly_2_comm_2,
            m,
        );
        
        assert!(result);
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));
        
        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));
        
        timings
    }

    // Run the full Baloo protocol with table and lookup generated based on k (legacy version)
    pub fn test_baloo_by_k(k: usize) -> Vec<String> {
        let (table, lookup) =
            generate_table_and_lookup(2_usize.pow(k as u32), 2_usize.pow((k - 1) as u32));
        Self::test_baloo_by_input(table, lookup)
    }
    
    // Run the optimized Baloo protocol with table and lookup generated based on k
    pub fn test_baloo_optimized_by_k(k: usize) -> Vec<String> {
        let (table, lookup) =
            generate_table_and_lookup(2_usize.pow(k as u32), 2_usize.pow((k - 1) as u32));
        Self::test_baloo_optimized_by_input(table, lookup)
    }

    /// Test Baloo with k parameter and N:n ratio using unified range check data
    pub fn test_baloo_by_k_with_ratio(k: usize, n_to_n_ratio: usize) -> Vec<String> {
        let (table, lookup) = crate::util::benchmark::generate_range_check_data(k, n_to_n_ratio);
        Self::test_baloo_by_input(table, lookup)
    }
    
    /// Test optimized Baloo with k parameter and N:n ratio using unified range check data
    pub fn test_baloo_optimized_by_k_with_ratio(k: usize, n_to_n_ratio: usize) -> Vec<String> {
        let (table, lookup) = crate::util::benchmark::generate_range_check_data(k, n_to_n_ratio);
        Self::test_baloo_optimized_by_input(table, lookup)
    }

    fn extract_prove_time(timings: &[String]) -> Option<u64> {
        for timing in timings {
            if timing.contains("Prove:") || timing.contains("Optimized Prove:") {
                // Extract number from string like "Prove: 123ms" or "Optimized Prove: 123ms"
                let parts: Vec<&str> = timing.split(':').collect();
                if parts.len() >= 2 {
                    let time_part = parts[1].trim().replace("ms", "");
                    return time_part.parse().ok();
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    type Pcs = UnivariateKzg<Bn256>;
    use std::cmp::max;
    use std::time::Instant;

    #[test]
    fn test_baloo() {
        // Generate table and lookup
        let (table, lookup) = generate_table_and_lookup(8, 4);

        let m = lookup.len();
        let t = table.len();
        let poly_size = max(t, m).next_power_of_two() * 2;
        let d = poly_size - 2;

        // 1. Setup
        let start = Instant::now();
        let (param, pp, vp) = Baloo::preprocess(t, m).unwrap();
        let duration1 = start.elapsed();
        println!(
            "\n ------------Setup and preprocess: {}ms----------- \n",
            duration1.as_millis()
        );

        // 2. Generate proof
        let start = Instant::now();
        let proof = Baloo::prove(&table, &param, &pp, &lookup);
        let duration2 = start.elapsed();

        // 3. Prepare verification data
        let start = Instant::now();
        let (
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
        ) = Baloo::prepare_verification_data(&param, &pp, &table, &lookup, m, t, d);

        // 4. Verify the proof
        let result = Baloo::verify(
            &vp,
            &proof,
            &t_comm_1,
            &z_h_comm_1,
            &phi_comm_1,
            &x_m_exponent_poly_comm_1,
            &x_exponent_poly_comm_2,
            &x_exponent_poly_2_comm_1,
            &x_exponent_poly_2_comm_2,
            m,
        );

        let duration3 = start.elapsed();
        println!(
            "\n ------------verify: {}ms----------- \n",
            duration3.as_millis()
        );

        assert!(result);
        println!("Finished to verify: baloo");
    }

    #[test]
    fn test_baloo_performance_comparison() {
        // Test performance comparison between legacy and optimized implementations
        println!("\n🚀 Testing Baloo Performance Comparison");
        println!("=======================================");
        
        for k in [8usize, 10usize] {
            println!("\n📊 Testing with K = {} (table size = {}, lookup size = {})", k, 2_usize.pow(k as u32), 2_usize.pow((k-1) as u32));
            
            // Test legacy implementation
            println!("\n🔄 Legacy Implementation:");
            let legacy_timings = Baloo::test_baloo_by_k(k);
            for timing in &legacy_timings {
                println!("  {}", timing);
            }
            
            // Test optimized implementation  
            println!("\n⚡ Optimized Implementation:");
            let optimized_timings = Baloo::test_baloo_optimized_by_k(k);
            for timing in &optimized_timings {
                println!("  {}", timing);
            }
            
            // Extract and compare proving times
            let legacy_prove_time = Baloo::extract_prove_time(&legacy_timings);
            let optimized_prove_time = Baloo::extract_prove_time(&optimized_timings);
            
            if let (Some(legacy), Some(optimized)) = (legacy_prove_time, optimized_prove_time) {
                let speedup = legacy as f64 / optimized as f64;
                println!("\n💨 Performance Improvement:");
                println!("  Legacy prove time: {}ms", legacy);
                println!("  Optimized prove time: {}ms", optimized);
                println!("  Speedup: {:.2}x", speedup);
                
                if speedup > 1.0 {
                    println!("  ✅ Optimization successful!");
                } else {
                    println!("  ⚠️  No significant improvement");
                }
            }
            
            println!("\n{}", "=".repeat(50));
        }
    }


    #[test]
    fn test_baloo_with_info() {
        // Generate table and lookup
        let (table, lookup) = generate_table_and_lookup(8, 4);

        let m = lookup.len();
        let t = table.len();
        let poly_size = max(t, m).next_power_of_two() * 2;
        let d = poly_size - 2;

        // 1. Setup using the info struct API
        let start = Instant::now();

        let info = BalooInfo {
            table: table.clone(),
            lookup: lookup.clone(),
        };

        let (pp, vp) = Baloo::preprocess_with_info(&info).unwrap();
        let duration1 = start.elapsed();
        println!(
            "\n ------------Setup and preprocess with info: {}ms----------- \n",
            duration1.as_millis()
        );

        // 2. Generate proof using the new API
        let start = Instant::now();
        let proof = Baloo::prove_with_param(&pp, &lookup);
        let duration2 = start.elapsed();

        // 3. Prepare verification data and verify
        let start = Instant::now();

        let (
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
        ) = Baloo::prepare_verification_data(&pp.param, &pp.pp, &table, &lookup, m, t, d);

        // Verify the proof
        let result = Baloo::verify_with_param(
            &vp,
            &proof,
            &t_comm_1,
            &z_h_comm_1,
            &phi_comm_1,
            &x_m_exponent_poly_comm_1,
            &x_exponent_poly_comm_2,
            &x_exponent_poly_2_comm_1,
            &x_exponent_poly_2_comm_2,
            m,
        );

        let duration3 = start.elapsed();
        println!(
            "\n ------------verify with param: {}ms----------- \n",
            duration3.as_millis()
        );

        assert!(result);
        println!("Finished to verify: baloo with info");
    }

    #[test]
    fn test_baloo_by_input() {
        let table_size = 2_usize.pow(6);
        let lookup_size = 4;

        // Generate table and lookup
        let (table, lookup) = generate_table_and_lookup(table_size, lookup_size);

        let timings = Baloo::test_baloo_by_input(table, lookup);

        // Print all timing information
        for timing in &timings {
            println!("{}", timing);
        }
    }
}
