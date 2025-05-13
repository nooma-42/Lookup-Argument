use std::{fmt::Debug, marker::PhantomData, iter, error::Error};
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
        let mut evals = vec![];
        (0..1 << LIMB_BITS).for_each(|i| evals.push(F::from(i as u64)));
        let limb_subtable_poly = MultilinearPolynomial::new(evals);
        if NUM_BITS % LIMB_BITS != 0 {
            let remainder = NUM_BITS % LIMB_BITS;
            let mut evals_rem = vec![]; // Renamed to avoid conflict
            (0..1 << remainder).for_each(|i| {
                evals_rem.push(F::from(i as u64));
            });
            let rem_subtable_poly = MultilinearPolynomial::new(evals_rem);
            vec![limb_subtable_poly, rem_subtable_poly]
        } else {
            vec![limb_subtable_poly]
        }
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


fn run_standalone_lasso_range_check() {
    // 0. Constants and Typedefs
    const NUM_BITS_FOR_RANGE_CHECK: usize = 8; // e.g., for u8 range
    const LIMB_BITS_FOR_RANGE_CHECK: usize = 4; // Decompose into 4-bit limbs
    type F = Fr;
    type Pcs = MultilinearKzg<Bn256>;
    type Transcript = Keccak256Transcript<F>; // Ensure this matches PCS field if necessary

    println!("Standalone Lasso Range Check: Proving values are in [0, 2^{})", NUM_BITS_FOR_RANGE_CHECK);

    // 1. Prover: Inputs
    let values_to_check: Vec<F> = vec![F::from(3u64), F::from(10u64), F::from(250u64), F::from(0u64)];
    println!("Values to check: {:?}", values_to_check.iter().map(|f| f.to_repr().as_ref()[0]).collect::<Vec<_>>());


    let num_lookups = values_to_check.len();
    if num_lookups == 0 {
        println!("No values to check. Skipping.");
    }
    let num_vars_for_lookups = (num_lookups as f64).log2().ceil() as usize; // num_vars for MLEs
    let padded_lookup_size = 1 << num_vars_for_lookups;

    let mut padded_values = values_to_check.clone();
    padded_values.resize(padded_lookup_size, *values_to_check.last().unwrap_or(&F::ZERO)); // Pad

    let lookup_index_poly = MultilinearPolynomial::new(padded_values.clone());
    // For RangeTable, if values are in range, output should be the same as input.
    // The E_i polys effectively represent the output from subtable lookups.
    // The final lookup_output_poly is T_eval in Lasso paper, which is sum of E_i weighted.
    // LassoProver::commit will construct the actual lookup_output_poly (T_eval) based on E_i polys.
    // We need to provide an initial "claimed" output, which for range check, is the input itself.
    let claimed_lookup_output_poly = MultilinearPolynomial::new(padded_values.clone());


    let table: Box<dyn DecomposableTable<F>> = Box::new(RangeTable::<F, NUM_BITS_FOR_RANGE_CHECK, LIMB_BITS_FOR_RANGE_CHECK>::new());
    let subtable_polys_mle: Vec<MultilinearPolynomial<F>> = table.subtable_polys();
    let subtable_polys_mle_refs: Vec<&MultilinearPolynomial<F>> = subtable_polys_mle.iter().collect();

    // Setup PCS
    let mut rng = std_rng();
    // Determine max polynomial size for PCS setup:
    // lookup_index/output polys, E_polys are size `padded_lookup_size`.
    // Dims, ReadTs, FinalCts polys depend on chunk_bits and num_vars_for_lookups.
    // Subtable_polys are size `1 << LIMB_BITS` or `1 << remainder_bits`.
    let max_limb_poly_size = 1 << LIMB_BITS_FOR_RANGE_CHECK;
    let remainder_bits = NUM_BITS_FOR_RANGE_CHECK % LIMB_BITS_FOR_RANGE_CHECK;
    let max_rem_poly_size = if remainder_bits > 0 { 1 << remainder_bits } else { 0 };
    let max_poly_degree_for_pcs = padded_lookup_size
        .max(max_limb_poly_size)
        .max(max_rem_poly_size);

    // Estimate batch size: output (1) + dims (N_chunks) + read_ts (N_chunks) + final_cts (N_chunks) + e_polys (N_memories)
    let num_chunks = table.chunk_bits().len();
    let num_memories = table.num_memories();
    let estimated_batch_size = 1 + num_chunks * 3 + num_memories + 1; // +1 for lookup_index_poly if committed separately

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


    // 3. Prover: Squeeze challenges for sumcheck and memory checking
    let beta_challenge: F = transcript_p.squeeze_challenge(); // For memory checking (as tau)
    let gamma_challenge: F = transcript_p.squeeze_challenge(); // For memory checking

    // For sumcheck, the point `r` (or `y` in hyperplonk) is where the identity is claimed.
    // This point is typically derived from verifier's challenges over the field.
    let r_sumcheck_q_challenges: Vec<F> = transcript_p.squeeze_challenges(prover_lookup_output_poly_struct.num_vars());
    println!("Prover: Challenges squeezed.");

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
    LassoProver::<F, Pcs>::prove_sum_check(
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
    println!("Prover: Sum-check prove done.");

    // 5. Prover: Memory Checking
    // `points_offset` for memory checking evaluations.
    // This needs to be distinct and consistent with how verifier reconstructs it.
    // Let's assume `memory_checking` appends its own distinct points.
    // The `polys_offset` inside `MemoryCheckingProver` is relative to its set of inputs.
    // The `points_offset` argument here is a base for evaluation point indices for memory checking.
    let memory_checking_base_point_idx = prover_lookup_opening_points.len();

    LassoProver::<F, Pcs>::memory_checking(
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
        &mut transcript_p,
    );
    println!("Prover: Memory-checking prove done.");

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

    Pcs::batch_open(
        &pcs_pp,
        all_polys_to_open_refs.into_iter(), // Must be Iterator<&'a MultilinearPolynomial<F>>
        all_comms_to_open_refs.into_iter(), // Must be Iterator<&'a Pcs::Commitment>
        &prover_lookup_opening_points,
        &prover_lookup_opening_evals,
        &mut transcript_p,
    );
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
    assert_eq!(verifier_lasso_comms.len(), estimated_batch_size -1 , "Mismatch in number of verifier commitments");


    // 2. Verifier: Squeeze challenges
    let beta_v: F = transcript_v.squeeze_challenge();
    let gamma_v: F = transcript_v.squeeze_challenge();
    assert_eq!(beta_v, beta_challenge, "Beta challenge mismatch");
    assert_eq!(gamma_v, gamma_challenge, "Gamma challenge mismatch");

    let r_sumcheck_q_challenges_v: Vec<F> = transcript_v.squeeze_challenges(lookup_index_poly.num_vars()); // Assuming output_poly has same num_vars
    assert_eq!(r_sumcheck_q_challenges_v, r_sumcheck_q_challenges, "Sumcheck q challenges mismatch");
    println!("Verifier: Challenges squeezed and matched.");

    let mut verifier_lookup_opening_points: Vec<Vec<F>> = Vec::new();
    let mut verifier_lookup_opening_evals: Vec<Evaluation<F>> = Vec::new();

    // 3. Verifier: Verify Sum-check
    // `polys_offset` is the index of the output polynomial's commitment in `verifier_lasso_comms`.
    // `points_offset` is the index of `r_sumcheck_q_challenges_v` in `verifier_lookup_opening_points`.
    LassoVerifier::<F, Pcs>::verify_sum_check(
        &table,
        lookup_index_poly.num_vars(), // num_vars for the sumcheck
        output_poly_commitment_idx,   // PCS index of the output polynomial in verifier_lasso_comms
        0,                            // Index of the point `r_sumcheck_q_challenges_v` for output poly eval
        &mut verifier_lookup_opening_points,
        &mut verifier_lookup_opening_evals,
        &r_sumcheck_q_challenges_v,
        &mut transcript_v,
    );
    println!("Verifier: Sum-check verify done.");

    // 4. Verifier: Verify Memory Checking
    // `polys_offset` needs to be the base index for dim/read_ts/final_cts commitments in `verifier_lasso_comms`.
    // `points_offset` is the base index for points related to memory checking in `verifier_lookup_opening_points`.
    let memory_checking_base_point_idx_v = verifier_lookup_opening_points.len();
    LassoVerifier::<F, Pcs>::memory_checking(
        padded_lookup_size, // num_reads
        first_dim_commitment_idx, // Base PCS index for Dims, ReadTs, FinalCts, E_polys groups
        memory_checking_base_point_idx_v,
        &mut verifier_lookup_opening_points,
        &mut verifier_lookup_opening_evals,
        &table,
        &gamma_v,
        &beta_v, // tau is beta
        &mut transcript_v,
    );
    println!("Verifier: Memory-checking verify done.");

    // 5. Verifier: PCS Batch Verify
    Pcs::batch_verify(
        &pcs_vp,
        verifier_lasso_comms.iter(),
        &verifier_lookup_opening_points,
        &verifier_lookup_opening_evals,
        &mut transcript_v,
    );
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
}