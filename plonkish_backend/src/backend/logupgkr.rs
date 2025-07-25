// src/logupgkr.rs

pub mod prover;
pub mod verifier;
pub mod util;

use crate::{
    pcs::{multilinear::MultilinearKzg, PolynomialCommitmentScheme, Evaluation},
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::PrimeField,
        transcript::{FiatShamirTranscript, InMemoryTranscript, TranscriptRead, TranscriptWrite, Keccak256Transcript},
    },
    Error,
};
use halo2_curves::bn256::{Bn256, Fr, G1Affine, G2Affine}; // 直接導入具體類型
use std::{marker::PhantomData, time::Instant, collections::HashMap, io::Cursor};
use rand::rngs::OsRng;

// 類型別名現在是具體的，不再有泛型
type Pcs = MultilinearKzg<Bn256>;
type OwnedTranscript = Keccak256Transcript<Cursor<Vec<u8>>>;

// --- 公共結構體定義 (無泛型) ---

#[derive(Clone, Debug)]
pub struct VerifierParam {
    pub pcs_vp: <Pcs as PolynomialCommitmentScheme<Fr>>::VerifierParam,
    pub num_vars: usize,
    pub num_vars_row: usize,
    pub m_comm: <Pcs as PolynomialCommitmentScheme<Fr>>::Commitment,
    pub t_comm: <Pcs as PolynomialCommitmentScheme<Fr>>::Commitment,
    pub w_comms: Vec<<Pcs as PolynomialCommitmentScheme<Fr>>::Commitment>,
}

#[derive(Clone)]
pub struct ProverParam {
    pub pcs_pp: <Pcs as PolynomialCommitmentScheme<Fr>>::ProverParam,
    pub m_poly: MultilinearPolynomial<Fr>,
    pub t_poly: MultilinearPolynomial<Fr>,
    pub w_polys: Vec<MultilinearPolynomial<Fr>>,
}

// --- 主結構體和核心 API (無泛型) ---
pub struct LogupGkr {
    _marker: PhantomData<()>,
}

impl LogupGkr {
    pub fn setup(
        m_poly: &MultilinearPolynomial<Fr>,
        t_poly: &MultilinearPolynomial<Fr>,
        w_polys: &[MultilinearPolynomial<Fr>],
    ) -> Result<(ProverParam, VerifierParam), Error> {
        let num_vars_row = t_poly.num_vars();
        let num_vars_col = if w_polys.is_empty() { 0 } else { 
            let total_cols = w_polys.len() + 1; // +1 for the t column
            // We need ceil(log2(total_cols)) bits to index all columns
            // For powers of 2, this is exactly log2, but for non-powers we need to round up
            if total_cols.is_power_of_two() {
                total_cols.ilog2() as usize
            } else {
                (total_cols.ilog2() + 1) as usize
            }
        };
        let num_vars = num_vars_row + num_vars_col;
        assert_eq!(m_poly.num_vars(), num_vars_row, "Mismatch in m_poly num_vars");
        for w_poly in w_polys {
            assert_eq!(w_poly.num_vars(), num_vars_row, "Mismatch in w_poly num_vars");
        }
        let num_polys_to_commit = 2 + w_polys.len();
        let max_poly_size = 1 << num_vars_row;
        let pcs_param = Pcs::setup(max_poly_size, num_polys_to_commit, OsRng)?;
        let (pcs_pp, pcs_vp) = Pcs::trim(&pcs_param, max_poly_size, num_polys_to_commit)?;
        let m_comm = Pcs::commit(&pcs_pp, m_poly)?;
        let t_comm = Pcs::commit(&pcs_pp, t_poly)?;
        let w_comms = Pcs::batch_commit(&pcs_pp, w_polys.iter())?;
        let pp = ProverParam { pcs_pp, m_poly: m_poly.clone(), t_poly: t_poly.clone(), w_polys: w_polys.to_vec() };
        let vp = VerifierParam { pcs_vp, num_vars, num_vars_row, m_comm, t_comm, w_comms };
        Ok((pp, vp))
    }

    pub fn prove(pp: &ProverParam) -> Result<Vec<u8>, Error> {
        let mut transcript = OwnedTranscript::default();
        prover::prove(pp, &mut transcript)?;
        Ok(transcript.into_proof())
    }
    
    pub fn verify(vp: &VerifierParam, proof: &[u8]) -> Result<(), Error> {
        let mut transcript = OwnedTranscript::from_proof((), proof);
        verifier::verify(vp, &mut transcript)
    }

    /// 一個端到端的測試函數
    pub fn test_logupgkr(
        m_values: HashMap<Vec<bool>, Fr>,
        t_values: HashMap<Vec<bool>, Fr>,
        w_values: Vec<HashMap<Vec<bool>, Fr>>,
    ) -> Vec<String> {
        let mut timings: Vec<String> = vec![];
        let start_total = Instant::now();

        let m_poly = util::create_multilinear_poly(m_values);
        let t_poly = util::create_multilinear_poly(t_values);
        let w_polys: Vec<_> = w_values.into_iter().map(util::create_multilinear_poly).collect();
        
        println!("Running setup...");
        let start = Instant::now();
        let (pp, vp) = Self::setup(&m_poly, &t_poly, &w_polys).unwrap();
        let duration1 = start.elapsed();
        timings.push(format!("Setup: {}ms", duration1.as_millis()));

        println!("Generating proof...");
        let start = Instant::now();
        let proof = Self::prove(&pp).unwrap();
        let duration2 = start.elapsed();
        timings.push(format!("Prove: {}ms", duration2.as_millis()));
        timings.push(format!("Proof size: {} bytes", proof.len()));

        println!("Verifying proof...");
        let start = Instant::now();
        Self::verify(&vp, &proof).unwrap();
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));

        timings
    }

    /// Test LogupGKR with k parameter and N:n ratio using unified range check data
    pub fn test_logupgkr_by_k_with_ratio(k: usize, n_to_n_ratio: usize) -> Vec<String> {
        // Generate unified range check data
        let (table, lookup) = crate::util::benchmark::generate_range_check_data(k, n_to_n_ratio);
        
        // Convert to logupgkr format using the utility function
        let (m_values, t_values, w_values) = util::convert_to_logupgkr_format(lookup, table);
        
        // Run the test with converted data
        Self::test_logupgkr(m_values, t_values, w_values)
    }
}




#[cfg(test)]
mod test {
    use super::*;
    use std::collections::HashMap;
    use halo2_curves::bn256::Fr;

    #[test]
    fn logupgkr_e2e_test_kzg() {
        let m_values = HashMap::from([
            (vec![false, false], Fr::from(3u64)),
            (vec![false, true], Fr::from(3u64)),
            (vec![true, false], Fr::from(1u64)),
            (vec![true, true], Fr::from(1u64)),
        ]);

        let t_values = HashMap::from([
            (vec![false, false], Fr::from(1u64)),
            (vec![false, true], Fr::from(2u64)),
            (vec![true, false], Fr::from(3u64)),
            (vec![true, true], Fr::from(4u64)),
        ]);

        let w_values = vec![
            HashMap::from([
                (vec![false, false], Fr::from(1u64)),
                (vec![false, true], Fr::from(2u64)),
                (vec![true, false], Fr::from(3u64)),
                (vec![true, true], Fr::from(1u64)),
            ]),
            HashMap::from([
                (vec![false, false], Fr::from(2u64)),
                (vec![false, true], Fr::from(1u64)),
                (vec![true, false], Fr::from(4u64)),
                (vec![true, true], Fr::from(2u64)),
            ]),
        ];

        let timings = LogupGkr::test_logupgkr(m_values, t_values, w_values);

        println!("\n--- LogupGkr with KZG E2E Test ---");
        for timing in &timings {
            println!("{}", timing);
        }
    }
}