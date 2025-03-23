use crate::{
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::PrimeField,
        transcript::{InMemoryTranscript, Keccak256Transcript, FieldTranscriptWrite, FieldTranscriptRead},
        Deserialize, DeserializeOwned, Serialize,
    },
    Error,
};
use halo2_curves::ff::WithSmallOrderMulGroup;
use std::{fmt::Debug, hash::Hash, marker::PhantomData, time::Instant, collections::HashMap};

pub mod preprocessor;
pub mod prover;
pub mod verifier;
pub mod util;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LogupGkrProverParam<F>
where
    F: PrimeField,
{
    m_poly: MultilinearPolynomial<F>,
    t_poly: MultilinearPolynomial<F>,
    w_polys: Vec<MultilinearPolynomial<F>>,
    a: F,
    ps: MultilinearPolynomial<F>,
    qs: MultilinearPolynomial<F>,
    p_0s: Vec<Option<F>>,
    q_0s: Vec<Option<F>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LogupGkrVerifierParam<F>
where
    F: PrimeField,
{
    num_vars: usize,
    p_0s: Vec<Option<F>>,
    q_0s: Vec<Option<F>>,
}

#[derive(Clone, Debug)]
pub struct LogupGkrInfo<F> {
    m_values: HashMap<Vec<bool>, F>,
    t_values: HashMap<Vec<bool>, F>,
    w_values: Vec<HashMap<Vec<bool>, F>>,
    a: F,
}

#[derive(Clone, Debug)]
pub struct LogupGkr<F>(PhantomData<F>);

impl<F> LogupGkr<F>
where
    F: PrimeField + WithSmallOrderMulGroup<3> + Hash + Serialize + DeserializeOwned,
{
    pub fn preprocess(
        info: &LogupGkrInfo<F>,
    ) -> Result<(LogupGkrProverParam<F>, LogupGkrVerifierParam<F>), Error> {
        preprocessor::preprocess(info)
    }

    pub fn prove(
        pp: LogupGkrProverParam<F>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(), Error> {
        prover::prove(pp, transcript)
    }

    pub fn verify(
        vp: LogupGkrVerifierParam<F>,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(), Error> {
        verifier::verify(vp, transcript)
    }
}

// Concrete implementation for BN256/Fr
use halo2_curves::bn256::{Bn256, Fr};

impl LogupGkr<Fr> {
    // Run the full LogupGkr protocol with given parameters
    pub fn test_logupgkr(
        m_values: HashMap<Vec<bool>, Fr>,
        t_values: HashMap<Vec<bool>, Fr>,
        w_values: Vec<HashMap<Vec<bool>, Fr>>,
        a: Fr,
    ) -> Vec<String> {
        let mut timings: Vec<String> = vec![];
        let start_total = Instant::now();

        // 1. Setup
        let info = LogupGkrInfo {
            m_values,
            t_values,
            w_values,
            a,
        };

        let start = Instant::now();
        let (pp, vp) = Self::preprocess(&info).unwrap();
        let duration1 = start.elapsed();
        timings.push(format!("Preprocess: {}ms", duration1.as_millis()));

        // 2. Prove
        let start = Instant::now();
        let mut transcript = Keccak256Transcript::new(());
        Self::prove(pp, &mut transcript).unwrap();
        let proof = transcript.into_proof();
        let duration2 = start.elapsed();
        timings.push(format!("Prove: {}ms", duration2.as_millis()));

        // 3. Verify
        let start = Instant::now();
        let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
        Self::verify(vp, &mut transcript).unwrap();
        let duration3 = start.elapsed();
        timings.push(format!("Verify: {}ms", duration3.as_millis()));

        let total_duration = start_total.elapsed();
        timings.push(format!("Total time: {}ms", total_duration.as_millis()));

        timings
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::backend::logupgkr::util::{generate_binary_combinations, create_multilinear_poly};
    use std::collections::HashMap;
    use halo2_curves::bn256::Fr;

    #[test]
    fn logupgkr_e2e_test1() {
        /* 
        witness (w_values): (w1, w2) 2 cols
        1, 2
        2, 1
        3, 4
        1, 2

        multiplicity (m_values): 
        3, # 3 1s in witness
        3, # 3 2s in witness
        1, # 1 3s in witness 
        1, # 1 4s in witness

        table (t_values):
        1, # 1 is first element in table being looked up
        2, # 2 is second element in table being looked up
        3, # 3 is third element in table being looked up
        4, # 4 is fourth element in table being looked up
         */
        // Create example parameter values
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

        let a = Fr::from(5u64);

        let timings = LogupGkr::test_logupgkr(m_values, t_values, w_values, a);

        // Print all timing information
        for timing in &timings {
            println!("{}", timing);
        }
    }

    #[test]
    fn logupgkr_e2e_test2() {
        /* 
        witness (w_values): (w1) 1 cols
        1
        2
        3
        1

        multiplicity (m_values): 
        2, # 2 1s in witness
        1, # 1 2s in witness
        1, # 1 3s in witness 
        1, # 1 4s in witness

        table (t_values):
        1, # 1 is first element in table being looked up
        2, # 2 is second element in table being looked up
        3, # 3 is third element in table being looked up
        4, # 4 is fourth element in table being looked up
        */

        let m_values = HashMap::from([
            (vec![false, false], Fr::from(2u64)),
            (vec![false, true], Fr::from(1u64)),
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
        ];

        let a = Fr::from(5u64);

        let timings = LogupGkr::test_logupgkr(m_values, t_values, w_values, a);
        
        
    }
}
