use rand::rngs::OsRng;
use std::{fmt::Debug, marker::PhantomData};

use halo2_curves::{bn256::{ pairing, Bn256, Fr, G1Affine, G2Affine, G2Prepared, Gt, G1, G2}, pairing::MillerLoopResult};

use crate::{
    poly::Polynomial,
    poly::univariate::UnivariatePolynomial,
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

use prover::Prover;
use verifier::Verifier;

// Specific implementation for Bn256 curves
#[derive(Clone, Debug)]
pub struct CqProverParam {
    param: UnivariateKzgParam<Bn256>,
    pp: UnivariateKzgProverParam<Bn256>,
    table: Vec<Fr>,
    q_t_comm_poly_coeffs: Vec<G1>,
}

#[derive(Clone, Debug)]
pub struct CqVerifierParam {
    vp: UnivariateKzgVerifierParam<Bn256>,
}

#[derive(Clone, Debug)]
pub struct CqInfo {
    table: Vec<Fr>,
    lookup: Vec<Fr>,
}

// Function to generate table and lookup for benchmarking
pub fn generate_table_and_lookup(table_size: usize, lookup_size: usize) -> (Vec<Fr>, Vec<Fr>) {
    let mut table = Vec::with_capacity(table_size);
    for i in 1..=table_size {
        table.push(Fr::from(i as u64));
    }
    
    // Generate lookup (random values from table)
    let mut lookup = Vec::with_capacity(lookup_size);
    for i in 0..lookup_size {
        // Use simple pattern for predictable lookups: select values from table modulo table_size
        let idx = (i % table_size) + 1;
        lookup.push(Fr::from(idx as u64));
    }
    
    (table, lookup)
}

#[derive(Clone, Debug)]
pub struct Cq;

impl Cq {
    // Main preprocess function that matches the original API
    pub fn preprocess(
        t: usize, 
        m: usize, 
        table: &Vec<Fr>
    ) -> Result<
        (
            UnivariateKzgParam<Bn256>,
            UnivariateKzgProverParam<Bn256>,
            UnivariateKzgVerifierParam<Bn256>,
            Vec<G1>, // q_t_comm_poly_coeffs
        ),
        Error
    > {
        preprocessor::preprocess(t, m, table)
    }

    // Alternative preprocess function that takes an info struct
    pub fn preprocess_with_info(
        info: &CqInfo,
    ) -> Result<
        (
            CqProverParam,
            CqVerifierParam,
            Vec<G1>, // q_t_comm_poly_coeffs
        ),
        Error
    > {
        let m = info.lookup.len();
        let t = info.table.len();
        
        let (param, pp, vp, q_t_comm_poly_coeffs) = 
            preprocessor::preprocess(t, m, &info.table)?;
        
        Ok((
            CqProverParam {
                param: param.clone(),
                pp,
                table: info.table.clone(),
                q_t_comm_poly_coeffs: q_t_comm_poly_coeffs.clone(),
            },
            CqVerifierParam {
                vp,
            },
            q_t_comm_poly_coeffs,
        ))
    }

    pub fn prove(
        table: &Vec<Fr>,
        param: &UnivariateKzgParam<Bn256>,
        pp: &UnivariateKzgProverParam<Bn256>,
        lookup: &Vec<Fr>,
        q_t_comm_poly_coeffs: &Vec<G1>,
    ) -> Vec<u8> {
        let prover = Prover::new(table, param, pp);
        prover.prove(lookup, q_t_comm_poly_coeffs)
    }

    pub fn prove_with_param(
        pp: &CqProverParam,
        lookup: &Vec<Fr>,
    ) -> Vec<u8> {
        let prover = Prover::new(&pp.table, &pp.param, &pp.pp);
        prover.prove(lookup, &pp.q_t_comm_poly_coeffs)
    }

    pub fn verify(
        vp: &UnivariateKzgVerifierParam<Bn256>,
        proof: &Vec<u8>,
        t_comm_2: &UnivariateKzgCommitment<G2Affine>,
        z_v_comm_2: &UnivariateKzgCommitment<G2Affine>,
        x_exponent_poly_comm_2: &UnivariateKzgCommitment<G2Affine>,
        m: usize,
        t: usize,
    ) -> bool {
        let verifier = Verifier::new(vp);
        verifier.verify(
            proof,
            t_comm_2,
            z_v_comm_2,
            x_exponent_poly_comm_2,
            m,
            t
        )
    }
    
    pub fn verify_with_param(
        vp: &CqVerifierParam,
        proof: &Vec<u8>,
        t_comm_2: &UnivariateKzgCommitment<G2Affine>,
        z_v_comm_2: &UnivariateKzgCommitment<G2Affine>,
        x_exponent_poly_comm_2: &UnivariateKzgCommitment<G2Affine>,
        m: usize,
        t: usize,
    ) -> bool {
        let verifier = Verifier::new(&vp.vp);
        verifier.verify(
            proof,
            t_comm_2,
            z_v_comm_2,
            x_exponent_poly_comm_2,
            m,
            t
        )
    }
    
    // Helper method to prepare verification data
    pub fn prepare_verification_data(
        param: &UnivariateKzgParam<Bn256>,
        table: &Vec<Fr>,
        m: usize,
        t: usize,
    ) -> (
        UnivariateKzgCommitment<G2Affine>, // t_comm_2
        UnivariateKzgCommitment<G2Affine>, // z_v_comm_2
        UnivariateKzgCommitment<G2Affine>, // x_exponent_poly_comm_2
    ) {
        let scalar_0 = Fr::zero();
        let scalar_1 = Fr::one();
        
        // z_v(x) = X^N - 1, [-1, 0, ..., 0, 1], t-1 0s in between
        let z_v_poly_coeffs = vec![scalar_1.neg()].into_iter().chain(vec![scalar_0; t - 1]).chain(vec![scalar_1]).collect();
        let z_v_poly = UnivariatePolynomial::monomial(z_v_poly_coeffs);
        // [z_h(x)]2
        let z_v_comm_2 = UnivariateKzg::<Bn256>::commit_monomial_g2(param, &z_v_poly.coeffs());
        
        // t(x)
        let t_poly = UnivariatePolynomial::lagrange(table.clone()).ifft();
        // [t(x)]2
        let t_comm_2 = UnivariateKzg::<Bn256>::commit_monomial_g2(param, &t_poly.coeffs());

        // [X^{N - 1 - (n - 2)}]2
        // x_exponent_poly_comm_2
        let x_exponent_order = t - 1 - (m - 2);
        let x_exponent_values_in_coeff = vec![scalar_0; x_exponent_order].into_iter().chain(vec![scalar_1]).collect();
        let x_exponent_poly = UnivariatePolynomial::monomial(x_exponent_values_in_coeff);
        // commit x_exponent_poly
        let x_exponent_poly_comm_2 = UnivariateKzg::<Bn256>::commit_monomial_g2(param, &x_exponent_poly.coeffs());
        
        (t_comm_2, z_v_comm_2, x_exponent_poly_comm_2)
    }
    
    // Run the full CQ protocol with given table and lookup
    pub fn test_cq_by_input(table: Vec<Fr>, lookup: Vec<Fr>) -> Vec<String> {
        let mut timings: Vec<String> = vec![];
        
        let start_total = std::time::Instant::now();
        
        let m = lookup.len();
        let t = table.len();
        
        // 1. Setup using the original API
        let start = std::time::Instant::now();
        let (param, pp, vp, q_t_comm_poly_coeffs) = Cq::preprocess(t, m, &table).unwrap();
        let duration1 = start.elapsed();
        timings.push(format!("Setup and preprocess: {}ms", duration1.as_millis()));
        
        // 2. Generate proof
        let start = std::time::Instant::now();
        let proof = Cq::prove(&table, &param, &pp, &lookup, &q_t_comm_poly_coeffs);
        let duration2 = start.elapsed();
        timings.push(format!("Prove: {}ms", duration2.as_millis()));
        
        // 3. Prepare verification data
        let start = std::time::Instant::now();
        let (t_comm_2, z_v_comm_2, x_exponent_poly_comm_2) = 
            Cq::prepare_verification_data(&param, &table, m, t);
        
        // 4. Verify
        let result = Cq::verify(
            &vp,
            &proof,
            &t_comm_2,
            &z_v_comm_2,
            &x_exponent_poly_comm_2,
            m,
            t
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
    use crate::util::transcript::{FieldTranscriptRead, FieldTranscriptWrite, G2TranscriptRead, G2TranscriptWrite};
    type Pcs = UnivariateKzg<Bn256>;
    use std::cmp::max;
    use std::time::Instant;
    
    #[test]
    fn test_cq() {
        let mut table = vec![];
        for k in 1..=2_usize.pow(6) {
            table.push(Fr::from(k as u64));
        }
        let lookup = vec![Fr::from(4), Fr::from(3), Fr::from(5), Fr::from(2)];
    
        let m = lookup.len();
        let t = table.len();
        let poly_size = max(t, m).next_power_of_two() * 2;
        
        // 1. Setup using the original API
        let start = Instant::now();
        let (param, pp, vp, q_t_comm_poly_coeffs) = Cq::preprocess(t, m, &table).unwrap();
        assert_eq!(poly_size, 2_usize.pow(pp.k() as u32));
        let duration1 = start.elapsed();
        println!("\n ------------Setup and preprocess: {}ms----------- \n", duration1.as_millis());

        // 2. Generate proof using the original API
        let start = Instant::now();
        let proof = Cq::prove(&table, &param, &pp, &lookup, &q_t_comm_poly_coeffs);
        println!("proof: {:?}", proof);
        let duration2 = start.elapsed();
        println!("\n ------------prove: {}ms----------- \n", duration2.as_millis());

        // 3. Verify using the original API
        let start = Instant::now();
        
        // Prepare verification data
        let (t_comm_2, z_v_comm_2, x_exponent_poly_comm_2) = 
            Cq::prepare_verification_data(&param, &table, m, t);
        
        // Verify the proof
        let result = Cq::verify(
            &vp,
            &proof,
            &t_comm_2,
            &z_v_comm_2,
            &x_exponent_poly_comm_2,
            m,
            t
        );
        
        let duration3 = start.elapsed();
        println!("\n ------------verify: {}ms----------- \n", duration3.as_millis());

        assert!(result);
        println!("Finished to verify: cq");
    }
    
    #[test]
    fn test_cq_with_info() {
        // Generate table and lookup
        let (table, lookup) = generate_table_and_lookup(2_usize.pow(6), 4);
    
        let m = lookup.len();
        let t = table.len();
        let poly_size = max(t, m).next_power_of_two() * 2;
        
        // 1. Setup using the info struct API
        let start = Instant::now();
        
        let info = CqInfo {
            table: table.clone(),
            lookup: lookup.clone(),
        };
        
        let (pp, vp, q_t_comm_poly_coeffs) = Cq::preprocess_with_info(&info).unwrap();
        let duration1 = start.elapsed();
        println!("\n ------------Setup and preprocess with info: {}ms----------- \n", duration1.as_millis());

        // 2. Generate proof using the new API
        let start = Instant::now();
        let proof = Cq::prove_with_param(&pp, &lookup);
        println!("proof: {:?}", proof);
        let duration2 = start.elapsed();
        println!("\n ------------prove with param: {}ms----------- \n", duration2.as_millis());

        // 3. Verify using the new API
        let start = Instant::now();
        
        // Prepare verification data
        let (t_comm_2, z_v_comm_2, x_exponent_poly_comm_2) = 
            Cq::prepare_verification_data(&pp.param, &table, m, t);
        
        // Verify the proof
        let result = Cq::verify_with_param(
            &vp,
            &proof,
            &t_comm_2,
            &z_v_comm_2,
            &x_exponent_poly_comm_2,
            m,
            t
        );
        
        let duration3 = start.elapsed();
        println!("\n ------------verify with param: {}ms----------- \n", duration3.as_millis());

        assert!(result);
        println!("Finished to verify: cq with info");
    }
    
    #[test]
    fn test_cq_by_input() {
        let table_size = 2_usize.pow(6);
        let lookup_size = 4;
        
        // Generate table and lookup
        let (table, lookup) = generate_table_and_lookup(table_size, lookup_size);
        
        let timings = Cq::test_cq_by_input(table, lookup);
        
        // Print all timing information
        for timing in &timings {
            println!("{}", timing);
        }
    }
}
