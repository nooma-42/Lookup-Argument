use rand::rngs::OsRng;
use std::{fmt::Debug, marker::PhantomData};
use serde::{Deserialize, Serialize};

use halo2_curves::{bn256::{ pairing, Bn256, Fr, G1Affine, G2Affine, G2Prepared, Gt, G1, G2}, pairing::MillerLoopResult};

use crate::{
    poly::Polynomial,
    poly::univariate::UnivariatePolynomial,
    backend::cq::generate_table_and_lookup,
    pcs::{
        PolynomialCommitmentScheme,
        Additive,
        univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam, UnivariateKzgCommitment},
    },
    util::{
        arithmetic::{Field, PrimeField, root_of_unity, variable_base_msm, barycentric_weights},
        test::std_rng,
        transcript::{InMemoryTranscript, TranscriptRead, TranscriptWrite, Keccak256Transcript},
    },
    Error,
};

pub mod preprocessor;
pub mod prover;
pub mod verifier;
pub mod util;

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
        m: usize
    ) -> Result<
        (
            UnivariateKzgParam<Bn256>,
            UnivariateKzgProverParam<Bn256>,
            UnivariateKzgVerifierParam<Bn256>,
        ),
        Error
    > {
        preprocessor::preprocess(t, m)
    }

    // Alternative preprocess function that takes an info struct
    pub fn preprocess_with_info(
        info: &BalooInfo,
    ) -> Result<
        (
            BalooProverParam,
            BalooVerifierParam,
        ),
        Error
    > {
        let m = info.lookup.len();
        let t = info.table.len();
        
        let (param, pp, vp) = 
            preprocessor::preprocess(t, m)?;
        
        let poly_size = std::cmp::max(t, m).next_power_of_two() * 2;
        let d = poly_size - 2;
        
        Ok((
            BalooProverParam {
                param: param.clone(),
                pp,
                table: info.table.clone(),
                d,
            },
            BalooVerifierParam {
                vp,
            },
        ))
    }

    pub fn prove(
        table: &Vec<Fr>,
        param: &UnivariateKzgParam<Bn256>,
        pp: &UnivariateKzgProverParam<Bn256>,
        lookup: &Vec<Fr>,
    ) -> Vec<u8> {
        let prover = prover::Prover::new(table, param, pp);
        prover.prove(lookup)
    }

    pub fn prove_with_param(
        pp: &BalooProverParam,
        lookup: &Vec<Fr>,
    ) -> Vec<u8> {
        let table_vec = pp.table.clone(); // Clone to avoid lifetime issues
        let lookup_vec = lookup.clone();
        let prover = prover::Prover::new(&table_vec, &pp.param, &pp.pp);
        prover.prove(&lookup_vec)
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
        let verifier = verifier::Verifier::new(vp);
        verifier.verify(
            proof,
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
            m
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
        let verifier = verifier::Verifier::new(&vp.vp);
        verifier.verify(
            proof,
            t_comm_1,
            z_h_comm_1,
            phi_comm_1,
            x_m_exponent_poly_comm_1,
            x_exponent_poly_comm_2,
            x_exponent_poly_2_comm_1,
            x_exponent_poly_2_comm_2,
            m
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
        
        let z_h_poly_coeffs = vec![scalar_1.neg()].into_iter().chain(vec![scalar_0; t - 1]).chain(vec![scalar_1]).collect();
        let z_h_poly = UnivariatePolynomial::monomial(z_h_poly_coeffs);
        let z_h_comm_1 = UnivariateKzg::<Bn256>::commit_monomial(pp, &z_h_poly.coeffs());
        
        let t_poly = UnivariatePolynomial::lagrange(table.to_vec()).ifft();
        let t_comm_1 = UnivariateKzg::<Bn256>::commit_monomial(pp, &t_poly.coeffs());

        let phi_poly = UnivariatePolynomial::lagrange(lookup.to_vec()).ifft();
        let phi_comm_1 = UnivariateKzg::<Bn256>::commit_monomial(pp, &phi_poly.coeffs());
        
        // X^m
        let x_m_exponent_poly = UnivariatePolynomial::monomial(vec![scalar_0; m].into_iter().chain(vec![scalar_1]).collect());
        let x_m_exponent_poly_comm_1 = UnivariateKzg::<Bn256>::commit_monomial(pp, &x_m_exponent_poly.clone().coeffs());

        // X^(d-m+1)
        let coeffs_x_exponent_poly = vec![scalar_0; d - m + 1].into_iter().chain(vec![scalar_1]).collect();
        let x_exponent_poly = UnivariatePolynomial::monomial(coeffs_x_exponent_poly);
        // [X^(d-m+1)]2
        let x_exponent_poly_comm_2 = UnivariateKzg::<Bn256>::commit_monomial_g2(param, &x_exponent_poly.coeffs());

        // X^(d-m+2)
        let coeffs_x_exponent_poly_2 = vec![scalar_0; d - m + 2].into_iter().chain(vec![scalar_1]).collect();
        let x_exponent_poly_2 = UnivariatePolynomial::monomial(coeffs_x_exponent_poly_2);
        let x_exponent_poly_2_comm_1 = UnivariateKzg::<Bn256>::commit_monomial(pp, &x_exponent_poly_2.coeffs());
        let x_exponent_poly_2_comm_2 = UnivariateKzg::<Bn256>::commit_monomial_g2(param, &x_exponent_poly_2.coeffs());
        
        (
            t_comm_1, 
            z_h_comm_1, 
            phi_comm_1, 
            x_m_exponent_poly_comm_1, 
            x_exponent_poly_comm_2, 
            x_exponent_poly_2_comm_1, 
            x_exponent_poly_2_comm_2
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
        
        // 3. Prepare verification data
        let start = std::time::Instant::now();
        let (
            t_comm_1, 
            z_h_comm_1, 
            phi_comm_1, 
            x_m_exponent_poly_comm_1, 
            x_exponent_poly_comm_2, 
            x_exponent_poly_2_comm_1, 
            x_exponent_poly_2_comm_2
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
            m
        );
        
        assert!(result);
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));
        
        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));
        
        timings
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
        println!("\n ------------Setup and preprocess: {}ms----------- \n", duration1.as_millis());

        // 2. Generate proof
        let start = Instant::now();
        let proof = Baloo::prove(&table, &param, &pp, &lookup);
        println!("proof: {:?}", proof);
        let duration2 = start.elapsed();
        println!("\n ------------prove: {}ms----------- \n", duration2.as_millis());

        // 3. Prepare verification data
        let start = Instant::now();
        let (
            t_comm_1, 
            z_h_comm_1, 
            phi_comm_1, 
            x_m_exponent_poly_comm_1, 
            x_exponent_poly_comm_2, 
            x_exponent_poly_2_comm_1, 
            x_exponent_poly_2_comm_2
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
            m
        );
        
        let duration3 = start.elapsed();
        println!("\n ------------verify: {}ms----------- \n", duration3.as_millis());

        assert!(result);
        println!("Finished to verify: baloo");
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
        println!("\n ------------Setup and preprocess with info: {}ms----------- \n", duration1.as_millis());

        // 2. Generate proof using the new API
        let start = Instant::now();
        let proof = Baloo::prove_with_param(&pp, &lookup);
        println!("proof: {:?}", proof);
        let duration2 = start.elapsed();
        println!("\n ------------prove with param: {}ms----------- \n", duration2.as_millis());

        // 3. Prepare verification data and verify
        let start = Instant::now();
        
        let (
            t_comm_1, 
            z_h_comm_1, 
            phi_comm_1, 
            x_m_exponent_poly_comm_1, 
            x_exponent_poly_comm_2, 
            x_exponent_poly_2_comm_1, 
            x_exponent_poly_2_comm_2
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
            m
        );
        
        let duration3 = start.elapsed();
        println!("\n ------------verify with param: {}ms----------- \n", duration3.as_millis());

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
