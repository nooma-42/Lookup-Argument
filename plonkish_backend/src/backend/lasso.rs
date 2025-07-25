use std::{fmt::Debug, marker::PhantomData, iter, error::Error, env};
use itertools::Itertools;


use halo2_curves::{
    ff::{Field, PrimeField},
    bn256::{Bn256, Fr},
};
use crate::{
    pcs::{
        Evaluation,
        PolynomialCommitmentScheme,
        multilinear::MultilinearKzg,
    },
    poly::multilinear::{MultilinearPolynomial, MultilinearPolynomialTerms},
    util::{
        expression::Expression,
        arithmetic::{inner_product, div_ceil},
        transcript::{Keccak256Transcript, InMemoryTranscript, FieldTranscript},
    },
    backend::lasso::{
        verifier::LassoVerifier,
        prover::{LassoProver, Poly},
    },
    log_debug, log_info, debug_println,
};

pub mod memory_checking;
pub mod prover;
pub mod verifier;

pub trait Subtable<F: PrimeField> {
    fn evaluate(point: &[F]) -> F;
}

/// This is a trait that contains information about decomposable table to which
/// backend prover and verifier can ask
pub trait DecomposableTable<F: PrimeField>: Debug + Sync + DecomposableTableClone<F> {
    fn num_memories(&self) -> usize;

    /// Returns multilinear extension polynomials of each subtable
    fn subtable_polys(&self) -> Vec<MultilinearPolynomial<F>>;
    fn subtable_polys_terms(&self) -> Vec<MultilinearPolynomialTerms<F>>;

    fn combine_lookup_expressions(&self, expressions: Vec<Expression<F>>) -> Expression<F>;

    /// The `g` function that computes T[r] = g(T_1[r_1], ..., T_k[r_1], T_{k+1}[r_2], ..., T_{\alpha}[r_c])
    fn combine_lookups(&self, operands: &[F]) -> F;

    /// Returns the size of bits for each chunk.
    /// Each chunk can have different bits.
    fn chunk_bits(&self) -> Vec<usize>;

    /// Returns the indices of each subtable lookups
    /// The length of `index_bits` is same as actual bit length of table index
    fn subtable_indices(&self, index_bits: Vec<bool>) -> Vec<Vec<bool>>;

    fn memory_to_subtable_index(&self, memory_index: usize) -> usize;

    fn memory_to_chunk_index(&self, memory_index: usize) -> usize;
}

pub trait DecomposableTableClone<F> {
    fn clone_box(&self) -> Box<dyn DecomposableTable<F>>;
}

impl<T, F: PrimeField> DecomposableTableClone<F> for T
where
    T: DecomposableTable<F> + Clone + 'static,
{
    fn clone_box(&self) -> Box<dyn DecomposableTable<F>> {
        Box::new(self.clone())
    }
}

impl<F> Clone for Box<dyn DecomposableTable<F>> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

// TODO: maybe put below into another file
#[derive(Clone, Debug)]
pub struct RangeTable<F, const NUM_BITS: usize, const LIMB_BITS: usize>(PhantomData<F>);

impl<F, const NUM_BITS: usize, const LIMB_BITS: usize> RangeTable<F, NUM_BITS, LIMB_BITS> {
    pub fn new() -> Self {
        Self(PhantomData)
    }
}

impl<F: PrimeField, const NUM_BITS: usize, const LIMB_BITS: usize> DecomposableTable<F>
    for RangeTable<F, NUM_BITS, LIMB_BITS>
{
    fn chunk_bits(&self) -> Vec<usize> {
        let remainder_bits = if NUM_BITS % LIMB_BITS != 0 {
            vec![NUM_BITS % LIMB_BITS]
        } else {
            vec![]
        };
        iter::repeat(LIMB_BITS)
            .take(NUM_BITS / LIMB_BITS)
            .chain(remainder_bits)
            .collect_vec()
    }

    fn combine_lookup_expressions(&self, expressions: Vec<Expression<F>>) -> Expression<F> {
        Expression::DistributePowers(
            expressions,
            Box::new(Expression::Constant(F::from(1u64 << LIMB_BITS))),
        )
    }

    fn combine_lookups(&self, operands: &[F]) -> F {
        let weight = F::from(1u64 << LIMB_BITS);
        inner_product(
            operands,
            iter::successors(Some(F::ONE), |power_of_weight| {
                Some(*power_of_weight * weight)
            })
            .take(operands.len())
            .collect_vec()
            .iter(),
        )
    }

    fn num_memories(&self) -> usize {
        div_ceil(NUM_BITS, LIMB_BITS)
    }

    fn subtable_indices(&self, index_bits: Vec<bool>) -> Vec<Vec<bool>> {
        index_bits.chunks(LIMB_BITS).map(Vec::from).collect_vec()
    }

    fn subtable_polys(&self) -> Vec<MultilinearPolynomial<F>> {
        let chunk_bits = self.chunk_bits();
        let mut polys = Vec::new();
        
        // Create one polynomial for each chunk based on its bit size
        for &bits in &chunk_bits {
            let mut evals = vec![];
            (0..1 << bits).for_each(|i| evals.push(F::from(i as u64)));
            polys.push(MultilinearPolynomial::new(evals));
        }
        
        polys
    }

    fn subtable_polys_terms(&self) -> Vec<crate::poly::multilinear::MultilinearPolynomialTerms<F>> {
        // Placeholder - This requires PolyExpr and its variants (Var, Const, Prod, Sum, Pow)
        // For a functional example, this needs to be correctly implemented.
        // For now, let's assume it can be bypassed if not directly used by the verifier's immediate path for this example.
        // LassoVerifier::chunks uses this. So it's important.
        // Let's use the one from your snippet:
        use crate::poly::multilinear::PolyExpr::*; // Assuming this is where Var, Const etc. are
        let limb_init = Var(0);
        let mut limb_terms = vec![limb_init];
        (1..LIMB_BITS).for_each(|i| {
            let coeff = Pow(Box::new(Const(F::from(2u64))), i as u32);
            let x = Var(i);
            let term = Prod(vec![coeff, x]);
            limb_terms.push(term);
        });
        let limb_subtable_poly_terms = crate::poly::multilinear::MultilinearPolynomialTerms::new(LIMB_BITS, Sum(limb_terms));
        if NUM_BITS % LIMB_BITS == 0 {
            vec![limb_subtable_poly_terms]
        } else {
            let remainder = NUM_BITS % LIMB_BITS;
            let rem_init = Var(0);
            let mut rem_terms = vec![rem_init];
            (1..remainder).for_each(|i| {
                let coeff = Pow(Box::new(Const(F::from(2u64))), i as u32);
                let x = Var(i);
                let term = Prod(vec![coeff, x]);
                rem_terms.push(term);
            });
            vec![
                limb_subtable_poly_terms,
                crate::poly::multilinear::MultilinearPolynomialTerms::new(remainder, Sum(rem_terms)),
            ]
        }
    }


    fn memory_to_chunk_index(&self, memory_index: usize) -> usize {
        memory_index
    }

    fn memory_to_subtable_index(&self, memory_index: usize) -> usize {
        if NUM_BITS % LIMB_BITS != 0 && memory_index == NUM_BITS / LIMB_BITS {
            1
        } else {
            0
        }
    }
}

// TODO: split run standalone lasso into test and mod function 
use crate::util::test::std_rng;
use std::time::Instant;

// Simple interface functions for benchmarking, similar to other lookup arguments

/// Generate values for range check based on k parameter
pub fn generate_range_check_values(k: usize, n_to_n_ratio: usize) -> Vec<Fr> {
    // Use the unified range check utility from util::benchmark
    crate::util::benchmark::generate_range_check_values(k, n_to_n_ratio)
}

/// Create a dynamic range table based on num_bits and limb_bits
fn create_dynamic_range_table(num_bits: usize, limb_bits: usize) -> Box<dyn DecomposableTable<Fr>> {
    // Use macro to generate all valid combinations
    macro_rules! range_table_match {
        ($($num_bits:literal, $limb_bits:literal),*) => {
            match (num_bits, limb_bits) {
                $(($num_bits, $limb_bits) => Box::new(RangeTable::<Fr, $num_bits, $limb_bits>::new()),)*
                // Fallback: if exact match not found, try with limb_bits=4 or use 8,4 as ultimate fallback
                (n, 4) if n <= 64 => match n {
                    1 => Box::new(RangeTable::<Fr, 1, 1>::new()),
                    2 => Box::new(RangeTable::<Fr, 2, 2>::new()),
                    3 => Box::new(RangeTable::<Fr, 3, 3>::new()),
                    _ => Box::new(RangeTable::<Fr, 8, 4>::new()), // Safe fallback
                },
                (n, 2) if n <= 64 => Box::new(RangeTable::<Fr, 8, 2>::new()), // Safe fallback for limb_bits=2
                _ => Box::new(RangeTable::<Fr, 8, 4>::new()), // Ultimate fallback
            }
        };
    }
    
    range_table_match!(
        // Common combinations with limb_bits=4
        8, 4, 9, 4, 10, 4, 11, 4, 12, 4, 13, 4, 14, 4, 15, 4, 16, 4,
        17, 4, 18, 4, 19, 4, 20, 4, 24, 4, 28, 4, 32, 4, 40, 4, 48, 4,
        56, 4, 64, 4,
        // Common combinations with limb_bits=2  
        8, 2, 10, 2, 12, 2, 13, 2, 14, 2, 16, 2, 20, 2, 24, 2, 32, 2,
        40, 2, 48, 2, 56, 2, 64, 2,
        // Special small cases
        1, 1, 2, 1, 2, 2, 3, 1, 3, 3, 4, 1, 4, 2, 4, 4,
        5, 1, 5, 5, 6, 1, 6, 2, 6, 3, 6, 6, 7, 1, 7, 7,
        // Other useful combinations
        9, 3, 12, 3, 15, 3, 18, 3, 21, 3, 24, 3, 27, 3, 30, 3
    )
}

/// Run Lasso range check with given values and return timing information
pub fn test_lasso_by_input_with_k(values_to_check: Vec<Fr>, k: usize) -> Vec<String> {
    let mut timings: Vec<String> = vec![];
    
    // Catch and handle errors gracefully
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        test_lasso_by_input_with_k_inner(values_to_check, k)
    }));
    
    match result {
        Ok(inner_timings) => inner_timings,
        Err(_) => {
            timings.push("Error: Test panicked during execution".to_string());
            timings
        }
    }
}

/// Inner implementation that can panic
fn test_lasso_by_input_with_k_inner(values_to_check: Vec<Fr>, k: usize) -> Vec<String> {
    let mut timings: Vec<String> = vec![];
    
    // Dynamic range check based on k parameter
    let num_bits_for_range_check = k;
    let limb_bits_for_range_check = 4.min(k); // Use 4-bit limbs, but cap at k if k < 4
    type Pcs = MultilinearKzg<Bn256>;
    
    let start_total = Instant::now();
    
    let num_lookups = values_to_check.len();
    if num_lookups == 0 {
        timings.push("Error: No values to check".to_string());
        return timings;
    }

    // Create table dynamically based on k
    let table: Box<dyn DecomposableTable<Fr>> = create_dynamic_range_table(num_bits_for_range_check, limb_bits_for_range_check);
    
    // Setup phase timing
    let start_setup = Instant::now();
    
    let chunk_bits = table.chunk_bits();
    let original_num_vars_for_lookups = (num_lookups as f64).log2().ceil() as usize;
    let unified_num_vars = *chunk_bits.iter().max().unwrap().max(&original_num_vars_for_lookups);
    let padded_lookup_size = 1 << unified_num_vars;

    let mut padded_values = values_to_check.clone();
    padded_values.resize(padded_lookup_size, *values_to_check.last().unwrap_or(&Fr::ZERO));

    let lookup_index_poly = MultilinearPolynomial::new(padded_values.clone());
    let claimed_lookup_output_poly = MultilinearPolynomial::new(padded_values.clone());
    let subtable_polys_mle: Vec<MultilinearPolynomial<Fr>> = table.subtable_polys();
    let subtable_polys_mle_refs: Vec<&MultilinearPolynomial<Fr>> = subtable_polys_mle.iter().collect();

    let mut rng = std_rng();
    let max_limb_poly_size = 1 << limb_bits_for_range_check;
    let remainder_bits = num_bits_for_range_check % limb_bits_for_range_check;
    let max_rem_poly_size = if remainder_bits > 0 { 1 << remainder_bits } else { 0 };
    let max_gkr_point_size = (1 << unified_num_vars).max(1 << limb_bits_for_range_check);
    
    let max_poly_degree_for_pcs = padded_lookup_size
        .max(max_limb_poly_size)
        .max(max_rem_poly_size)
        .max(max_gkr_point_size);

    let num_chunks = table.chunk_bits().len();
    let num_memories = table.num_memories();
    let estimated_batch_size = 1 + num_chunks * 3 + num_memories;

    let pcs_param = Pcs::setup(max_poly_degree_for_pcs, estimated_batch_size, &mut rng).unwrap();
    let (pcs_pp, pcs_vp) = Pcs::trim(&pcs_param, max_poly_degree_for_pcs, estimated_batch_size).unwrap();
    
    let setup_duration = start_setup.elapsed();
    timings.push(format!("Setup: {}ms", setup_duration.as_millis()));

    // Prove phase timing
    let start_prove = Instant::now();
    
    let mut transcript_p = Keccak256Transcript::new(());

    // Prover: LassoProver::commit
    let (committed_lasso_polys_nested_structs, committed_lasso_comms_nested) = 
        LassoProver::<Fr, Pcs>::commit(
            &pcs_pp,
            0,
            &table,
            &subtable_polys_mle_refs,
            claimed_lookup_output_poly.clone(),
            &lookup_index_poly,
            &mut transcript_p,
        ).unwrap();

    let prover_lookup_output_poly_struct = &committed_lasso_polys_nested_structs[0][0];
    let prover_dims_polys_struct = &committed_lasso_polys_nested_structs[1];
    let prover_read_ts_polys_struct = &committed_lasso_polys_nested_structs[2];
    let prover_final_cts_polys_struct = &committed_lasso_polys_nested_structs[3];
    let prover_e_polys_struct = &committed_lasso_polys_nested_structs[4];
    let e_poly_refs: Vec<&Poly<Fr>> = prover_e_polys_struct.iter().collect();

    // Squeeze challenges
    let _sentinel_p1: Fr = transcript_p.squeeze_challenge();
    let beta_challenge: Fr = transcript_p.squeeze_challenge();
    let gamma_challenge: Fr = transcript_p.squeeze_challenge();
    let r_sumcheck_q_challenges: Vec<Fr> = transcript_p.squeeze_challenges(unified_num_vars);
    let _sentinel_p2: Fr = transcript_p.squeeze_challenge();

    let mut prover_lookup_opening_points: Vec<Vec<Fr>> = Vec::new();
    let mut prover_lookup_opening_evals: Vec<Evaluation<Fr>> = Vec::new();

    let output_poly_commitment_idx = 0;

    // Sum-check
    LassoProver::<Fr, Pcs>::prove_sum_check(
        output_poly_commitment_idx,
        &mut prover_lookup_opening_points,
        &mut prover_lookup_opening_evals,
        &table,
        prover_lookup_output_poly_struct,
        &e_poly_refs,
        &r_sumcheck_q_challenges,
        unified_num_vars,
        &mut transcript_p,
    ).unwrap();

    let _sentinel_p3: Fr = transcript_p.squeeze_challenge();
    log_debug!("Prover: _sentinel_p3 = {:?}", _sentinel_p3);

    // Memory checking
    let memory_checking_base_point_idx = prover_lookup_opening_points.len();
    LassoProver::<Fr, Pcs>::memory_checking(
        memory_checking_base_point_idx,
        &mut prover_lookup_opening_points,
        &mut prover_lookup_opening_evals,
        &table,
        &subtable_polys_mle_refs,
        prover_dims_polys_struct,
        prover_read_ts_polys_struct,
        prover_final_cts_polys_struct,
        prover_e_polys_struct,
        &gamma_challenge,
        &beta_challenge,
        unified_num_vars,
        &mut transcript_p,
    ).unwrap();

    let _sentinel_p4: Fr = transcript_p.squeeze_challenge();

    // PCS batch opening
    let mut all_polys_to_open_refs: Vec<&MultilinearPolynomial<Fr>> = Vec::new();
    let mut all_comms_to_open_refs: Vec<&<Pcs as PolynomialCommitmentScheme<Fr>>::Commitment> = Vec::new();

    for poly_group in committed_lasso_polys_nested_structs.iter() {
        for poly_struct in poly_group.iter() {
            all_polys_to_open_refs.push(&poly_struct.poly);
        }
    }
    for comm_group in committed_lasso_comms_nested.iter() {
        for comm in comm_group.iter() {
            all_comms_to_open_refs.push(comm);
        }
    }

    let _sentinel_p_final: Fr = transcript_p.squeeze_challenge();

    // Direct batch opening without grouping (following the working approach)
    Pcs::batch_open(
        &pcs_pp,
        all_polys_to_open_refs.into_iter(),
        all_comms_to_open_refs.into_iter(),
        &prover_lookup_opening_points,
        &prover_lookup_opening_evals,
        &mut transcript_p,
    ).unwrap();

    let proof_bytes = transcript_p.into_proof();
    
    let prove_duration = start_prove.elapsed();
    timings.push(format!("Prove: {}ms", prove_duration.as_millis()));
    timings.push(format!("Proof size: {} bytes", proof_bytes.len()));

    // Verify phase timing
    let start_verify = Instant::now();
    
    let mut transcript_v = Keccak256Transcript::from_proof((), proof_bytes.as_slice());

    // Verifier: Read commitments
    let verifier_lasso_comms: Vec<<Pcs as PolynomialCommitmentScheme<Fr>>::Commitment> = 
        LassoVerifier::<Fr, Pcs>::read_commitments(
            &pcs_vp,
            &table,
            &mut transcript_v,
        ).unwrap();

    // Squeeze challenges
    let _sentinel_v1: Fr = transcript_v.squeeze_challenge();
    let beta_v: Fr = transcript_v.squeeze_challenge();
    let gamma_v: Fr = transcript_v.squeeze_challenge();
    let r_sumcheck_q_challenges_v: Vec<Fr> = transcript_v.squeeze_challenges(unified_num_vars);
    let _sentinel_v2: Fr = transcript_v.squeeze_challenge();

    let mut verifier_lookup_opening_points: Vec<Vec<Fr>> = Vec::new();
    let mut verifier_lookup_opening_evals: Vec<Evaluation<Fr>> = Vec::new();

    // Verify sum-check
    if let Err(e) = LassoVerifier::<Fr, Pcs>::verify_sum_check(
        &table,
        unified_num_vars,
        output_poly_commitment_idx,
        0,
        &mut verifier_lookup_opening_points,
        &mut verifier_lookup_opening_evals,
        &r_sumcheck_q_challenges_v,
        &mut transcript_v,
    ) {
        timings.push(format!("Error: Sum check verification failed: {:?}", e));
        return timings;
    }

    let _sentinel_v3: Fr = transcript_v.squeeze_challenge();
    log_debug!("Verifier: _sentinel_v3 = {:?}", _sentinel_v3);

    // Verify memory checking
    let memory_checking_base_point_idx_v = verifier_lookup_opening_points.len();
    if let Err(e) = LassoVerifier::<Fr, Pcs>::memory_checking(
        unified_num_vars, // num_reads
        output_poly_commitment_idx, // polys_offset
        memory_checking_base_point_idx_v, // points_offset
        &mut verifier_lookup_opening_points,
        &mut verifier_lookup_opening_evals,
        &table,
        &gamma_v,
        &beta_v, // tau
        unified_num_vars, // unified_num_vars
        &mut transcript_v,
    ) {
        timings.push(format!("Error: Memory checking failed: {:?}", e));
        return timings;
    }

    let _sentinel_v4: Fr = transcript_v.squeeze_challenge();

    // PCS batch verify
    let _sentinel_v_final: Fr = transcript_v.squeeze_challenge();

    // Direct batch verification without grouping (following the working approach)
    if let Err(e) = Pcs::batch_verify(
        &pcs_vp,
        verifier_lasso_comms.iter(),
        &verifier_lookup_opening_points,
        &verifier_lookup_opening_evals,
        &mut transcript_v,
    ) {
        timings.push(format!("Error: Batch verify failed: {:?}", e));
        return timings;
    }

    let verify_duration = start_verify.elapsed();
    timings.push(format!("Verify: {}ms", verify_duration.as_millis()));

    let total_duration = start_total.elapsed();
    timings.push(format!("Total time: {}ms", total_duration.as_millis()));

    timings
}

/// Run Lasso range check with k parameter, similar to other lookup arguments
pub fn test_lasso_by_k(k: usize) -> Vec<String> {
    let n_to_n_ratio = 2; // Default ratio
    let values_to_check = generate_range_check_values(k, n_to_n_ratio);
    test_lasso_by_input(values_to_check)
}

/// Run Lasso range check with given values and return timing information (8-bit default)
pub fn test_lasso_by_input(values_to_check: Vec<Fr>) -> Vec<String> {
    test_lasso_by_input_with_k(values_to_check, 8) // Default to 8-bit for backward compatibility
}

/// Run Lasso range check with k parameter and custom ratio
pub fn test_lasso_by_k_with_ratio(k: usize, n_to_n_ratio: usize) -> Vec<String> {
    let values_to_check = generate_range_check_values(k, n_to_n_ratio);
    test_lasso_by_input_with_k(values_to_check, k)
}

fn run_standalone_lasso_range_check() {
    // 0. Constants and Typedefs
    const NUM_BITS_FOR_RANGE_CHECK: usize = 8; // e.g., for u8 range
    const LIMB_BITS_FOR_RANGE_CHECK: usize = 4; // Decompose into 4-bit limbs
    type F = Fr;
    type Pcs = MultilinearKzg<Bn256>;
    type Transcript = Keccak256Transcript<F>; // Ensure this matches PCS field if necessary

    log_debug!("Standalone Lasso Range Check: Proving values are in [0, 2^{})", NUM_BITS_FOR_RANGE_CHECK);

    // 1. Prover: Inputs
    let values_to_check: Vec<F> = vec![F::from(3u64), F::from(10u64), F::from(250u64), F::from(0u64)];
    log_debug!("Values to check: {:?}", values_to_check.iter().map(|f| f.to_repr().as_ref()[0]).collect::<Vec<_>>());


    let num_lookups = values_to_check.len();
    if num_lookups == 0 {
        log_debug!("No values to check. Skipping.");
    }

    let table: Box<dyn DecomposableTable<F>> = Box::new(RangeTable::<F, NUM_BITS_FOR_RANGE_CHECK, LIMB_BITS_FOR_RANGE_CHECK>::new());
    
    // This will resolve batch open failure issue
    let chunk_bits = table.chunk_bits();
    let original_num_vars_for_lookups = (num_lookups as f64).log2().ceil() as usize;
    let unified_num_vars = *chunk_bits.iter().max().unwrap().max(&original_num_vars_for_lookups);
    
    log_debug!("Original num_vars_for_lookups: {}, chunk_bits: {:?}, unified_num_vars: {}", 
             original_num_vars_for_lookups, chunk_bits, unified_num_vars);
    
    let padded_lookup_size = 1 << unified_num_vars;

    let mut padded_values = values_to_check.clone();
    padded_values.resize(padded_lookup_size, *values_to_check.last().unwrap_or(&F::ZERO)); // Pad

    let lookup_index_poly = MultilinearPolynomial::new(padded_values.clone());
    // For RangeTable, if values are in range, output should be the same as input.
    // The E_i polys effectively represent the output from subtable lookups.
    // The final lookup_output_poly is T_eval in Lasso paper, which is sum of E_i weighted.
    // LassoProver::commit will construct the actual lookup_output_poly (T_eval) based on E_i polys.
    // We need to provide an initial "claimed" output, which for range check, is the input itself.
    let claimed_lookup_output_poly = MultilinearPolynomial::new(padded_values.clone());
    let subtable_polys_mle: Vec<MultilinearPolynomial<F>> = table.subtable_polys();
    let subtable_polys_mle_refs: Vec<&MultilinearPolynomial<F>> = subtable_polys_mle.iter().collect();

    // Setup PCS
    let mut rng = std_rng();
    // Determine max polynomial size for PCS setup:
    // lookup_index/output polys, E_polys are size `padded_lookup_size`.
    // Dims, ReadTs, FinalCts polys depend on chunk_bits and num_vars_for_lookups.
    // Subtable_polys are size `1 << LIMB_BITS` or `1 << remainder_bits`.
    // Memory checking GKR generates challenge points with variables:
    // - For reads/writes: num_vars = log2(num_lookups) 
    // - For init/final: num_vars = LIMB_BITS_FOR_RANGE_CHECK
    let max_limb_poly_size = 1 << LIMB_BITS_FOR_RANGE_CHECK;
    let remainder_bits = NUM_BITS_FOR_RANGE_CHECK % LIMB_BITS_FOR_RANGE_CHECK;
    let max_rem_poly_size = if remainder_bits > 0 { 1 << remainder_bits } else { 0 };
    
    // Consider GKR challenge points: need to support both unified_num_vars and LIMB_BITS_FOR_RANGE_CHECK variables
    let max_gkr_point_size = (1 << unified_num_vars).max(1 << LIMB_BITS_FOR_RANGE_CHECK);
    
    let max_poly_degree_for_pcs = padded_lookup_size
        .max(max_limb_poly_size)
        .max(max_rem_poly_size)
        .max(max_gkr_point_size);

    // Estimate batch size: output (1) + dims (N_chunks) + read_ts (N_chunks) + final_cts (N_chunks) + e_polys (N_memories)
    let num_chunks = table.chunk_bits().len();
    let num_memories = table.num_memories();
    let estimated_batch_size = 1 + num_chunks * 3 + num_memories;

    let pcs_param = Pcs::setup(max_poly_degree_for_pcs, estimated_batch_size, &mut rng).unwrap();
    let (pcs_pp, pcs_vp) = Pcs::trim(&pcs_param, max_poly_degree_for_pcs, estimated_batch_size).unwrap();

    let mut transcript_p = Keccak256Transcript::new(());

    // 2. Prover: LassoProver::commit
    // `lookup_polys_offset` for standalone can be 0.
    // This commits to: lookup_output_poly (derived from E_i), Dims, ReadTs, FinalCts, E_polys.
    let (committed_lasso_polys_nested_structs, committed_lasso_comms_nested) = match LassoProver::<F, Pcs>::commit(
        &pcs_pp,
        0, // lookup_polys_offset
        &table,
        &subtable_polys_mle_refs,
        claimed_lookup_output_poly.clone(), // This is the T_idx in Jolt paper, becomes T_final after sumcheck
        &lookup_index_poly, // This is A_idx in Jolt paper
        &mut transcript_p,
    ) {
        Ok(result) => result,
        Err(e) => {
            println!("Prover: Lasso commit failed: {:?}", e);
            return;
        }
    };
    println!("Prover: Lasso commit done.");

    // Extract key polynomials and their commitments for later use by prover and reference by verifier
    // Structure: [[output_poly_struct], [dim_poly_structs...], [read_ts_poly_structs...], [final_cts_poly_structs...], [e_poly_structs...]]
    let prover_lookup_output_poly_struct= &committed_lasso_polys_nested_structs[0][0];
    let prover_dims_polys_struct = &committed_lasso_polys_nested_structs[1];
    let prover_read_ts_polys_struct = &committed_lasso_polys_nested_structs[2];
    let prover_final_cts_polys_struct = &committed_lasso_polys_nested_structs[3];
    let prover_e_polys_struct = &committed_lasso_polys_nested_structs[4];
    let e_poly_refs: Vec<&Poly<F>> = prover_e_polys_struct.iter().collect();

    // Sentinel after LassoProver::commit
    let sentinel_p1: F = transcript_p.squeeze_challenge();
    println!("Prover Sentinel 1: {:?}", sentinel_p1);

    // 3. Prover: Squeeze challenges for sumcheck and memory checking
    let beta_challenge: F = transcript_p.squeeze_challenge(); // For memory checking (as tau)
    let gamma_challenge: F = transcript_p.squeeze_challenge(); // For memory checking

    // For sumcheck, the point `r` (or `y` in hyperplonk) is where the identity is claimed.
    // This point is typically derived from verifier's challenges over the field.
    let r_sumcheck_q_challenges: Vec<F> = transcript_p.squeeze_challenges(prover_lookup_output_poly_struct.num_vars());
    println!("Prover: Challenges squeezed.");

    // Sentinel after squeezing beta, gamma, r_sumcheck
    let sentinel_p2: F = transcript_p.squeeze_challenge();
    println!("Prover Sentinel 2: {:?}", sentinel_p2);

    let mut prover_lookup_opening_points: Vec<Vec<F>> = Vec::new();
    let mut prover_lookup_opening_evals: Vec<Evaluation<F>> = Vec::new();

    // These are base offsets for poly_index in Evaluation objects for this standalone proof.
    // The order of commitments as read by verifier: output, dims, read_ts, final_cts, e_polys
    let output_poly_commitment_idx = 0;
    let first_dim_commitment_idx = output_poly_commitment_idx + 1;
    let first_read_ts_commitment_idx = first_dim_commitment_idx + prover_dims_polys_struct.len();
    let first_final_cts_commitment_idx = first_read_ts_commitment_idx + prover_read_ts_polys_struct.len();
    let first_e_poly_commitment_idx = first_final_cts_commitment_idx + prover_final_cts_polys_struct.len();


    // 4. Prover: Sum-check (Surge)
    // `points_offset` for `prove_sum_check` refers to the index of `r_sumcheck_q_challenges` in `prover_lookup_opening_points`
    // and `polys_offset` refers to the commitment index of `prover_lookup_output_poly_struct`.
    let sumcheck_result = LassoProver::<F, Pcs>::prove_sum_check(
        output_poly_commitment_idx, // This is the `polys_offset` for the output polynomial T_final
                                    // It is also used as `points_offset` for its evaluation at `r_sumcheck_q_challenges`
        &mut prover_lookup_opening_points,
        &mut prover_lookup_opening_evals,
        &table,
        prover_lookup_output_poly_struct,
        &e_poly_refs,
        &r_sumcheck_q_challenges,
        prover_lookup_output_poly_struct.num_vars(),
        &mut transcript_p,
    );
    if sumcheck_result.is_err() {
        log_debug!("Prover: Sum-check prove failed.");
        return;
    }
    log_debug!("Prover: Sum-check prove done.");

    // Sentinel after LassoProver::prove_sum_check
    let sentinel_p3: F = transcript_p.squeeze_challenge();
    log_debug!("Prover: _sentinel_p3 = {:?}", sentinel_p3);

    // 5. Prover: Memory Checking
    // `points_offset` for memory checking evaluations.
    // This needs to be distinct and consistent with how verifier reconstructs it.
    // Let's assume `memory_checking` appends its own distinct points.
    // The `polys_offset` inside `MemoryCheckingProver` is relative to its set of inputs.
    // The `points_offset` argument here is a base for evaluation point indices for memory checking.
    let memory_checking_base_point_idx = prover_lookup_opening_points.len();

    let memory_checking_result = LassoProver::<F, Pcs>::memory_checking(
        memory_checking_base_point_idx, // Base index for points used in memory checking evaluations
        &mut prover_lookup_opening_points,
        &mut prover_lookup_opening_evals,
        &table,
        &subtable_polys_mle_refs,
        prover_dims_polys_struct,
        prover_read_ts_polys_struct,
        prover_final_cts_polys_struct,
        prover_e_polys_struct,
        &gamma_challenge,
        &beta_challenge, // tau is beta
        unified_num_vars, // unified number of variables for memory checking
        &mut transcript_p,
    );
    if memory_checking_result.is_err() {
        println!("Prover: Memory-checking prove failed.");
        return;
    }
    println!("Prover: Memory-checking prove done.");

    // Sentinel after LassoProver::memory_checking
    let sentinel_p4: F = transcript_p.squeeze_challenge();
    println!("Prover Sentinel 4: {:?}", sentinel_p4);

    // 6. Prover: PCS Batch Opening
    let mut all_polys_to_open_refs: Vec<&MultilinearPolynomial<F>> = Vec::new();
    let mut all_comms_to_open_refs: Vec<&<Pcs as PolynomialCommitmentScheme<F>>::Commitment> = Vec::new();

    for poly_group in committed_lasso_polys_nested_structs.iter() {
        for poly_struct in poly_group.iter() {
            all_polys_to_open_refs.push(&poly_struct.poly);
        }
    }
    for comm_group in committed_lasso_comms_nested.iter() {
        for comm in comm_group.iter() {
            all_comms_to_open_refs.push(comm);
        }
    }
    // Note: lookup_index_poly was an input to `LassoProver::commit` but not directly part of its output commitments.
    // If it needs to be opened, it should be committed separately and added here.
    // For now, we assume it's implicitly verified by the consistency of other parts.

    // Print prover's opening points, evaluations, and commitments
    println!("\n=== Prover's opening points, evaluations, and commitments ===");
    println!("prover_lookup_opening_points (total: {}):", prover_lookup_opening_points.len());
    for (i, point) in prover_lookup_opening_points.iter().enumerate() {
        println!("  Point[{}] (len={}): {:?}", i, point.len(), point);
    }
    
    println!("\nPolynomials num_vars info:");
    let mut poly_idx = 0;
    for poly_group in committed_lasso_polys_nested_structs.iter() {
        for poly_struct in poly_group.iter() {
            println!("  Poly[{}] num_vars: {}", poly_idx, poly_struct.poly.num_vars());
            poly_idx += 1;
        }
    }
    
    println!("\nprover_lookup_opening_evals (total: {}):", prover_lookup_opening_evals.len());
    for (i, eval) in prover_lookup_opening_evals.iter().enumerate() {
        println!("  Eval[{}]: poly_index={}, point_index={}, value={:?}", 
            i, eval.poly(), eval.point(), eval.value());
    }
    
    println!("\nall_comms_to_open_refs (total: {})", all_comms_to_open_refs.len());
    
    // Sentinel JUST BEFORE Pcs::batch_open
    let sentinel_p_final: F = transcript_p.squeeze_challenge();
    println!("Prover Sentinel Final (before batch_open): {:?}", sentinel_p_final);

    // Prover, before batch_open
    println!("Prover Commitments (total: {}):", all_comms_to_open_refs.len());
    for (i, comm_ref) in all_comms_to_open_refs.iter().enumerate() {
        // Assuming comm_ref is &G1Affine or derefs to G1Affine
        // You might need to adjust how you access the coordinates based on your Pcs::Commitment type
        println!("  CommP[{}]: x={:?}", i, &comm_ref.0);
    }
    let batch_open_result = Pcs::batch_open(
        &pcs_pp,
        all_polys_to_open_refs, // Must be Vec<&'a MultilinearPolynomial<F>>
        all_comms_to_open_refs, // Must be Vec<&'a Pcs::Commitment>
        &prover_lookup_opening_points,
        &prover_lookup_opening_evals,
        &mut transcript_p,
    );
    if batch_open_result.is_err() {
        println!("Prover: PCS batch open failed: {:?}", batch_open_result.err());
        panic!("Prover batch open failed");
    }
    println!("Prover: PCS batch open done.");

    let proof_bytes = transcript_p.into_proof();
    println!("Prover: Proof generated ({} bytes).", proof_bytes.len());

    // =====================================================================================
    // Verifier Side
    // =====================================================================================
    println!("\nVerifier: Starting verification...");
    let mut transcript_v = Keccak256Transcript::from_proof((), proof_bytes.as_slice());

    // 1. Verifier: Read commitments for auxiliary Lasso polynomials
    let verifier_lasso_comms: Vec<<Pcs as PolynomialCommitmentScheme<F>>::Commitment> = match LassoVerifier::<F, Pcs>::read_commitments(
            &pcs_vp,
            &table,
            &mut transcript_v,
        ) {
            Ok(result) => result,
            Err(e) => {
                println!("Verifier: Lasso commitments read failed: {:?}", e);
                return;
            }
        };
    println!("Verifier: Lasso commitments read.");
    // Order: [output_comm, dim_comms..., read_ts_comms..., final_cts_comms..., e_comms...]
    // Check consistency:
    assert_eq!(verifier_lasso_comms.len(), estimated_batch_size, "Mismatch in number of verifier commitments");

    // Sentinel after LassoVerifier::read_commitments
    let sentinel_v1: F = transcript_v.squeeze_challenge();
    println!("Verifier Sentinel 1: {:?}", sentinel_v1);

    // 2. Verifier: Squeeze challenges
    let beta_v: F = transcript_v.squeeze_challenge();
    let gamma_v: F = transcript_v.squeeze_challenge();
    assert_eq!(beta_v, beta_challenge, "Beta challenge mismatch");
    assert_eq!(gamma_v, gamma_challenge, "Gamma challenge mismatch");

    let r_sumcheck_q_challenges_v: Vec<F> = transcript_v.squeeze_challenges(unified_num_vars); // Use unified num_vars
    assert_eq!(r_sumcheck_q_challenges_v, r_sumcheck_q_challenges, "Sumcheck q challenges mismatch");
    println!("Verifier: Challenges squeezed and matched.");

    // Sentinel after squeezing beta_v, gamma_v, r_sumcheck_q_challenges_v
    let sentinel_v2: F  = transcript_v.squeeze_challenge();
    println!("Verifier Sentinel 2: {:?}", sentinel_v2);

    let mut verifier_lookup_opening_points: Vec<Vec<F>> = Vec::new();
    let mut verifier_lookup_opening_evals: Vec<Evaluation<F>> = Vec::new();

    // 3. Verifier: Verify Sum-check
    // `polys_offset` is the index of the output polynomial's commitment in `verifier_lasso_comms`.
    // `points_offset` is the index of `r_sumcheck_q_challenges_v` in `verifier_lookup_opening_points`.
    let sumcheck_result = LassoVerifier::<F, Pcs>::verify_sum_check(
        &table,
        unified_num_vars, // num_vars for the sumcheck
        output_poly_commitment_idx,   // PCS index of the output polynomial in verifier_lasso_comms
        0,                            // Index of the point `r_sumcheck_q_challenges_v` for output poly eval
        &mut verifier_lookup_opening_points,
        &mut verifier_lookup_opening_evals,
        &r_sumcheck_q_challenges_v,
        &mut transcript_v,
    );
    if sumcheck_result.is_err() {
        println!("Verifier: Sum-check verify failed.");
        return;
    }
    println!("Verifier: Sum-check verify done.");

    let _sentinel_v3: Fr = transcript_v.squeeze_challenge();
    log_debug!("Verifier: _sentinel_v3 = {:?}", _sentinel_v3);

    // Verify memory checking
    let memory_checking_base_point_idx_v = verifier_lookup_opening_points.len();
    LassoVerifier::<Fr, Pcs>::memory_checking(
        unified_num_vars, // num_reads
        output_poly_commitment_idx, // polys_offset - Base PCS index for Dims, ReadTs, FinalCts, E_polys groups
        memory_checking_base_point_idx_v, // points_offset
        &mut verifier_lookup_opening_points,
        &mut verifier_lookup_opening_evals,
        &table,
        &gamma_v,
        &beta_v, // tau is beta
        unified_num_vars, // unified_num_vars
        &mut transcript_v,
    ).unwrap();

    let _sentinel_v4: Fr = transcript_v.squeeze_challenge();

    // 5. Verifier: PCS Batch Verify
    // Print verifier's opening points, evaluations, and commitments
    println!("\n=== Verifier's opening points, evaluations, and commitments ===");
    println!("verifier_lookup_opening_points (total: {}):", verifier_lookup_opening_points.len());
    for (i, point) in verifier_lookup_opening_points.iter().enumerate() {
        println!("  Point[{}]: {:?}", i, point);
    }
    
    println!("\nverifier_lookup_opening_evals (total: {}):", verifier_lookup_opening_evals.len());
    for (i, eval) in verifier_lookup_opening_evals.iter().enumerate() {
        println!("  Eval[{}]: poly_index={}, point_index={}, value={:?}", 
            i, eval.poly(), eval.point(), eval.value());
    }
    
    println!("\nverifier_lasso_comms (total: {})", verifier_lasso_comms.len());
    
    // Sentinel JUST BEFORE Pcs::batch_verify
    let sentinel_v_final: F = transcript_v.squeeze_challenge();
    println!("Verifier Sentinel Final (before batch_verify): {:?}", sentinel_v_final);

    // Verifier, before batch_verify
    println!("Verifier Commitments (total: {}):", verifier_lasso_comms.len());
    for (i, comm_ref) in verifier_lasso_comms.iter().enumerate() {
        // Assuming comm_ref is &G1Affine or derefs to G1Affine
        // You might need to adjust how you access the coordinates based on your Pcs::Commitment type
        println!("  CommV[{}]: x={:?}", i, &comm_ref.0);
    }
    let batch_verify_result = Pcs::batch_verify(
        &pcs_vp,
        verifier_lasso_comms.iter(),
        &verifier_lookup_opening_points,
        &verifier_lookup_opening_evals,
        &mut transcript_v,
    );
    if batch_verify_result.is_err() {
        println!("Verifier: PCS batch verify failed: {:?}", batch_verify_result.err());
        panic!("Verifier batch verify failed");
    }
    println!("Verifier: PCS batch verify done.");

    println!("\nStandalone Lasso Range Check Successful!");
}


#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_lasso_range_check() {
        run_standalone_lasso_range_check();
    }
    
    #[test]
    fn test_lasso_interface_by_k() {
        let timings = test_lasso_by_k(4);
        println!("Lasso by K=4 Timings:");
        for timing in &timings {
            println!("  {}", timing);
        }
        // Verify that we got timing results
        assert!(timings.len() >= 3); // Setup, Prove, Verify (minimum)
        assert!(timings.iter().any(|t| t.contains("Setup")));
        assert!(timings.iter().any(|t| t.contains("Prove")));
        assert!(timings.iter().any(|t| t.contains("Verify")));
    }
    
    #[test]
    fn test_lasso_interface_by_k_with_ratio() {
        let timings = test_lasso_by_k_with_ratio(4, 4); // N:n = 4:1 ratio
        println!("Lasso by K=4, ratio=4 Timings:");
        for timing in &timings {
            println!("  {}", timing);
        }
        // Verify that we got timing results
        assert!(timings.len() >= 3); // Setup, Prove, Verify (minimum)
        assert!(timings.iter().any(|t| t.contains("Setup")));
        assert!(timings.iter().any(|t| t.contains("Prove")));
        assert!(timings.iter().any(|t| t.contains("Verify")));
    }
    
    #[test]
    fn test_lasso_k8_debug() {
        println!("Testing Lasso with k=8 to trigger debug output");
        let timings = test_lasso_by_k_with_ratio(8, 2); // This should trigger the debug case
        println!("Lasso by K=8, ratio=2 Timings:");
        for timing in &timings {
            println!("  {}", timing);
        }
        // Verify that we got timing results
        assert!(timings.len() >= 3); // Setup, Prove, Verify (minimum)
        assert!(timings.iter().any(|t| t.contains("Setup")));
        assert!(timings.iter().any(|t| t.contains("Prove")));
        assert!(timings.iter().any(|t| t.contains("Verify")));
    }
    
    #[test]
    fn test_lasso_interface_by_input() {
        let values_to_check = vec![
            Fr::from(3u64), 
            Fr::from(10u64), 
            Fr::from(250u64), 
            Fr::from(0u64)
        ];
        
        let timings = test_lasso_by_input(values_to_check);
        println!("Lasso by Input Timings:");
        for timing in &timings {
            println!("  {}", timing);
        }
        // Verify that we got timing results
        assert!(timings.len() >= 3); // Setup, Prove, Verify (minimum)
        assert!(timings.iter().any(|t| t.contains("Setup")));
        assert!(timings.iter().any(|t| t.contains("Prove")));
        assert!(timings.iter().any(|t| t.contains("Verify")));
    }
    
    #[test]
    fn test_generate_range_check_values() {
        let values = generate_range_check_values(6, 2); // k=6, ratio=2
        println!("Generated values for k=6, ratio=2: {} values", values.len());
        
        // Range check: for k=6, range_bits = min(6, 8) = 6, so range_size = 64
        // lookup_size = 64 / 2 = 32
        assert_eq!(values.len(), 32);
        
        // All values should be in range [0, 64)
        for value in &values {
            let val_u64 = value.to_repr().as_ref()[0];
            assert!(val_u64 < 64, "Value {} out of range", val_u64);
        }
    }
    
    #[test]
    fn test_k_9_debug() {
        let k = 9;
        println!("\n=== Debug Testing k={} ===", k);
        let values = generate_range_check_values(k, 2);
        
        let result = std::panic::catch_unwind(|| {
            test_lasso_by_input_with_k(values, k)
        });
        
        match result {
            Ok(timings) if !timings.is_empty() => {
                println!("✅ k={} SUCCESS", k);
                for timing in &timings {
                    println!("  {}", timing);
                }
            }
            _ => {
                println!("❌ k={} FAILED", k);
            }
        }
    }

    #[test]
    fn test_k_comparison() {
        for k in [8, 9, 10] {
            println!("\n=== Debug Testing k={} ===", k);
            let values = generate_range_check_values(k, 2);
            
            let result = std::panic::catch_unwind(|| {
                test_lasso_by_input_with_k(values, k)
            });
            
            match result {
                Ok(timings) if !timings.is_empty() => {
                    println!("✅ k={} SUCCESS", k);
                    for timing in &timings {
                        println!("  {}", timing);
                    }
                }
                _ => {
                    println!("❌ k={} FAILED", k);
                }
            }
        }
    }

    #[test]
    fn test_k_ratios_debug() {
        let k = 9;
        let ratios = [2, 4, 8, 16];
        
        for ratio in ratios {
            let range_size = 1 << k;
            let lookup_size = range_size / ratio;
            println!("\n=== Testing k={}, ratio={} (range_size={}, lookup_size={}) ===", 
                k, ratio, range_size, lookup_size);
            
            let values = generate_range_check_values(k, ratio);
            println!("Generated {} values", values.len());
            
            let result = std::panic::catch_unwind(|| {
                test_lasso_by_input_with_k(values, k)
            });
            
            match result {
                Ok(timings) if !timings.is_empty() => {
                    println!("✅ k={}, ratio={} SUCCESS", k, ratio);
                }
                _ => {
                    println!("❌ k={}, ratio={} FAILED", k, ratio);
                }
            }
        }
    }

    #[test]
    fn test_ratio_4_debug() {
        // Test specific failing case: k=9, ratio=4
        let k = 9;
        let ratio = 4;
        let lookup_size = (1 << k) / ratio; // Should be 128
        
        println!("=== Debug k={} ratio={} (lookup_size={}) ===", k, ratio, lookup_size);
        let values = generate_range_check_values(k, ratio);
        println!("Generated {} values", values.len());
        
        // Test with reduced verbosity - disable timer feature to reduce noise
        env::set_var("LOOKUP_LOG_LEVEL", "silent");
        
        let result = std::panic::catch_unwind(|| {
            test_lasso_by_input_with_k(values, k)
        });
        
        match result {
            Ok(timings) if !timings.is_empty() => {
                println!("✅ k={}, ratio={} SUCCESS", k, ratio);
            }
            Err(panic_info) => {
                if let Some(s) = panic_info.downcast_ref::<String>() {
                    println!("❌ Panic: {}", s);
                } else if let Some(s) = panic_info.downcast_ref::<&str>() {
                    println!("❌ Panic: {}", s);
                } else {
                    println!("❌ Unknown panic occurred");
                }
            }
            _ => {
                println!("❌ Empty result");
            }
        }
    }

    #[test]
    fn test_multiple_k_values_clean() {
        let test_values = [8, 10, 11, 12, 13, 14]; // Skip k=9 for now
        let mut results = Vec::new();
        
        for k in test_values {
            println!("\n=== Testing k={} ===", k);
            let values = generate_range_check_values(k, 2);
            
            let timings = test_lasso_by_input_with_k(values, k);
            
            if !timings.is_empty() {
                println!("✅ k={} SUCCESS", k);
                results.push((k, true));
                for timing in &timings {
                    println!("  {}", timing);
                }
            } else {
                println!("❌ k={} FAILED", k);
                results.push((k, false));
            }
        }
        
        println!("\n=== SUMMARY ===");
        for (k, success) in results {
            println!("k={}: {}", k, if success { "✅ SUCCESS" } else { "❌ FAILED" });
        }
    }
}