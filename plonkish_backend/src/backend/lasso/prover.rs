// plonkish_backend/src/backend/lasso/prover.rs
use super::{
    LassoProof, Message1, Message2, Message3, Message4, Message5,
    GrandProductData, util::{g_func_simple_range, hash_tuple}
};
use crate::{
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::{Field, PrimeField},
        transcript::{FieldTranscriptWrite, TranscriptWrite, TranscriptRead}, // Use TranscriptWrite for squeeze
        expression::{CommonPolynomial, Expression, Query, Rotation},
        // timer::Timer,
    },
    Error,
    piop::sum_check::{self, classic::{ClassicSumCheck, CoefficientsProver}, VirtualPolynomial},
    pcs::PolynomialCommitmentScheme,
    pcs::{Evaluation, Point}, // Import Point type
};
use std::{marker::PhantomData, fmt::Debug, hash::Hash};
use halo2_curves::ff::WithSmallOrderMulGroup;

// --- Virtual Polynomial for h = eq(r, x) * g(E_1(x), ..., E_alpha(x)) ---
// Remove HSumcheckPoly struct and impl


// --- Virtual Polynomial for Grand Product Check ---
// Remove GrandProductSumcheckPoly struct and impl


// --- Main Prove Function ---

pub fn prove<F, Pcs>(
    mut pp: LassoProverParam<F, Pcs>,
    witness: Vec<F>, // Witness is actually not used directly here, it was used in preprocessing
    transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
) -> Result<LassoProof<F, Pcs>, Error>
where
    F: PrimeField + WithSmallOrderMulGroup<3> + Hash,
    F: PrimeField + WithSmallOrderMulGroup<3> + Hash + Eq + Send + Sync,
    // Explicit bound for transcript needed due to Pcs::open constraint
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    // Bounds match PolynomialCommitmentScheme associated types
    Pcs::Commitment: Default + Clone + Debug + PartialEq + Send + Sync + AsRef<[Pcs::CommitmentChunk]>,
    Pcs::CommitmentChunk: Clone + Debug + Default + Send + Sync,
    Pcs::ProverParam: Clone + Debug + Send + Sync,
    Pcs::VerifierParam: Clone + Debug + Send + Sync,
    Pcs::Param: Clone + Debug + Send + Sync,
{
    // Cast transcript to FieldTranscriptWrite where needed for field operations
    let transcript_f = transcript.borrow_mut();

    // --- Round 1: Commitments and Challenges ---
    let msg1 = Message1 {
        a_comm: pp.a_comm.clone(),
        logm: pp.logm,
        dim_comm: pp.dim_comm.clone(),
        _marker: PhantomData, // Add marker
    };
    transcript_f.write_serializable(b"msg1", &msg1)?; // Prover writes msg1
    let r: Point<F, Pcs::Polynomial> = transcript_f.squeeze_challenges(pp.logm);

    // --- Round 2: Evaluate 'a', Commit E, read, final ---
    // Evaluate a(r) and open the commitment.
    let a_eval = pp.a_poly.evaluate(&r);
    // Use non-mutable pcs_param
    Pcs::open(&pp.pcs_param, &pp.a_poly, &pp.a_comm, &r, &a_eval, transcript)?;

    let msg2 = Message2 {
        a_eval,
        E_comm: pp.E_comm.clone(),
        read_comm: pp.read_comm.clone(),
        final_comm: pp.final_comm.clone(),
        placeholder_proofs: PhantomData, // Add placeholder
    };
    transcript_f.write_serializable(b"msg2", &msg2)?; // Prover writes msg2

    // --- Round 3: Sumcheck for h = eq(r, x) * g(E_1(x), ..., E_alpha(x)) ---
    // Define the expression for H
    let mut g_expr = Expression::Constant(F::ZERO);
    let mut power_of_2_l = F::ONE;
    let base = F::from(2u64).pow([pp.l as u64, 0, 0, 0]); // 2^l
    let e_poly_offset = 0; // Assuming E polys start at index 0 in the polys vector passed to VirtualPolynomial
    for j in 0..pp.c {
        let mut chunk_expr = Expression::Constant(F::ZERO);
        for i in 0..pp.k {
            let idx = j * pp.k + i;
            // Use PolynomialIndex instead of integer offset
            chunk_expr = chunk_expr + Expression::Polynomial(Query::new(idx, Rotation::cur()));
        }
        // Multiply by scalar constant correctly
        g_expr = g_expr + (chunk_expr * Expression::Constant(power_of_2_l));
        power_of_2_l *= base;
    }
    let h_expr = Expression::CommonPolynomial(CommonPolynomial::EqXY(0)) * g_expr;

    // Define the virtual polynomial for sumcheck
    // The `polys` vector should contain references to all polynomials used in the expression.
    // In this case, just the E_poly vector.
    let e_poly_refs: Vec<&MultilinearPolynomial<F>> = pp.E_poly.iter().collect();
    let virtual_poly_h = VirtualPolynomial::new(
        &h_expr,
        e_poly_refs, // Pass refs to E polynomials
        &[], // No external challenges used in H expression itself
        &[r.clone()], // Pass challenge point `r` as the first `y` vector for EqXY(0)
    );

    // Run sumcheck prover
    let (_h_sumcheck_proof, rz, h_final_eval) = ClassicSumCheck::<CoefficientsProver<F>, F>::prove(
        &(), // Prover param for sumcheck (currently empty tuple)
        pp.logm, // num_vars
        virtual_poly_h, // Pass the constructed VirtualPolynomial
        a_eval, // claimed sum H(r) = a(r)
        transcript_f, // Use the FieldTranscriptWrite capable transcript
    )?;
    // The third element returned by prove is the evaluation of the virtual polynomial at the random point rz
    let E_eval_at_rz = h_final_eval; // Renamed for clarity

    // Open E commitments at rz
    // We need the actual evaluations E_i(rz) first.
    let mut E_eval = Vec::with_capacity(pp.alpha);
    for i in 0..pp.alpha {
        let eval = pp.E_poly[i].evaluate(&rz);
        E_eval.push(eval);
        // Use non-mutable pcs_param
        Pcs::open(&pp.pcs_param, &pp.E_poly[i], &pp.E_comm[i], &rz, &eval, transcript)?;
    }

    // Check consistency: Sumcheck provides final evaluation at rz
    // Need to recalculate H(rz) using the opened E_eval values
    let eq_r_rz = MultilinearPolynomial::eq_xy_evaluate(&r, &rz);
    let g_at_rz = g_func_simple_range(&E_eval, pp.l, pp.c, pp.k);
    let h_eval_check = eq_r_rz * g_at_rz;
    if h_eval_check != E_eval_at_rz { // Compare calculated H(rz) with sumcheck result
         return Err(Error::InvalidSnark(format!(
             "H sumcheck prover final eval mismatch: expected {:?}, got {:?}", h_eval_check, E_eval_at_rz
         )));
    }

    let msg3 = Message3 {
        // h_sumcheck_proof, // Removed: proof is implicitly in transcript
        rz: rz.clone(), // Sumcheck returns the point rz
        E_eval, // Use the evaluations opened above
        placeholder_proofs: PhantomData,
    };
    transcript_f.write_serializable(b"msg3", &msg3)?; // Prover writes msg3
    let (tau, gamma) = {
        let challenges = transcript_f.squeeze_challenges(2);
        (challenges[0], challenges[1])
    };

    // --- Round 4: Grand Product Polynomials and Commitments ---
    // Build and Commit S0, S, RS, WS polynomials
    // Removed parallelization
    let mut S0_poly = Vec::with_capacity(pp.alpha);
    let mut S0_comm = Vec::with_capacity(pp.alpha);
    let mut S_poly = Vec::with_capacity(pp.alpha);
    let mut S_comm = Vec::with_capacity(pp.alpha);
    let mut RS_poly = Vec::with_capacity(pp.alpha);
    let mut RS_comm = Vec::with_capacity(pp.alpha);
    let mut WS_poly = Vec::with_capacity(pp.alpha);
    let mut WS_comm = Vec::with_capacity(pp.alpha);

    for i in 0..pp.alpha {
        let chunk_idx = i / pp.k;
        // Ensure final_poly exists and has correct evals
        let final_vals = pp.final_poly.get(i).ok_or(Error::InvalidSnark("Missing final_poly".to_string()))?.evals();
        let E_vals = pp.E_poly.get(i).ok_or(Error::InvalidSnark("Missing E_poly".to_string()))?.evals();
        let read_vals = pp.read_poly.get(i).ok_or(Error::InvalidSnark("Missing read_poly".to_string()))?.evals();
        let write_vals = pp.write_poly.get(i).ok_or(Error::InvalidSnark("Missing write_poly".to_string()))?.evals();
        let dim_values: Vec<F> = pp.dim_poly.get(chunk_idx).ok_or(Error::InvalidSnark("Missing dim_poly".to_string()))?.evals().to_vec();

        let subtable_domain_size = 1 << pp.l; // Domain size for S0/S is 2^l
        let access_domain_size = 1 << pp.logm; // Domain size for RS/WS is 2^logm

        let s0_tuples: Vec<(F, F, F)> = (0..subtable_domain_size)
            .map(|idx| (F::from(idx as u64), pp.subtables[i][idx], F::ZERO)) // Use F::ZERO for timestamp
            .collect();
        let poly_s0 = build_grand_product_poly(&s0_tuples, pp.l, &tau, &gamma)?;
        // Use non-mutable pcs_param
        let comm_s0 = Pcs::commit(&pp.pcs_param, &poly_s0)?;

        let s_tuples: Vec<(F, F, F)> = (0..subtable_domain_size)
             .map(|idx| (F::from(idx as u64), pp.subtables[i][idx], final_vals[idx])) // Uses final_poly evals
             .collect();
        let poly_s = build_grand_product_poly(&s_tuples, pp.l, &tau, &gamma)?;
        // Use non-mutable pcs_param
        let comm_s = Pcs::commit(&pp.pcs_param, &poly_s)?;

        let rs_tuples: Vec<(F, F, F)> = (0..access_domain_size)
             .map(|access| {
                // Handle padding for dim_values, E_vals, read_vals if access >= witness.len()
                let idx = access.min(witness.len().saturating_sub(1)); // Use last valid index if padded
                (dim_values[idx], E_vals[idx], read_vals[idx])
            })
             .collect();
        let poly_rs = build_grand_product_poly(&rs_tuples, pp.logm, &tau, &gamma)?; // Use logm for num_vars_base
        // Use non-mutable pcs_param
        let comm_rs = Pcs::commit(&pp.pcs_param, &poly_rs)?;

        let ws_tuples: Vec<(F, F, F)> = (0..access_domain_size)
             .map(|access| {
                 // Handle padding
                 let idx = access.min(witness.len().saturating_sub(1));
                 (dim_values[idx], E_vals[idx], write_vals[idx])
             })
             .collect();
        let poly_ws = build_grand_product_poly(&ws_tuples, pp.logm, &tau, &gamma)?; // Use logm for num_vars_base
        // Use non-mutable pcs_param
        let comm_ws = Pcs::commit(&pp.pcs_param, &poly_ws)?;

        S0_poly.push(poly_s0);
        S0_comm.push(comm_s0);
        S_poly.push(poly_s);
        S_comm.push(comm_s);
        RS_poly.push(poly_rs);
        RS_comm.push(comm_rs);
        WS_poly.push(poly_ws);
        WS_comm.push(comm_ws);
    }

    let msg4 = Message4 { S0_comm: S0_comm.clone(), S_comm: S_comm.clone(), RS_comm: RS_comm.clone(), WS_comm: WS_comm.clone(), _marker: PhantomData }; // Add marker
    transcript_f.write_serializable(b"msg4", &msg4)?; // Prover writes msg4

    // --- Round 5: Grand Product Sumchecks and Final Evaluations ---
    let mut r_prime = Vec::with_capacity(pp.alpha);
    let mut r_prime2 = Vec::with_capacity(pp.alpha);
    let mut r_prime3 = Vec::with_capacity(pp.alpha);
    let mut r_prime4 = Vec::with_capacity(pp.alpha);
    let mut S0_data = Vec::with_capacity(pp.alpha);
    let mut S_data = Vec::with_capacity(pp.alpha);
    let mut RS_data = Vec::with_capacity(pp.alpha);
    let mut WS_data = Vec::with_capacity(pp.alpha);
    let mut E_eval2 = Vec::with_capacity(pp.alpha);
    let mut dim_eval = Vec::with_capacity(pp.alpha);
    let mut read_eval = Vec::with_capacity(pp.alpha);
    let mut final_eval = Vec::with_capacity(pp.alpha);

    // Store grand product polys temporarily
    pp.S0_poly = Some(S0_poly);
    pp.S_poly = Some(S_poly);
    pp.RS_poly = Some(RS_poly);
    pp.WS_poly = Some(WS_poly);

    for i in 0..pp.alpha {
        // S0 Sumcheck
        let (_proof0, rp0, _eval0) = run_grand_product_sumcheck::<_, Pcs>(&pp.S0_poly.as_ref().unwrap()[i], pp.l, transcript)?;
        let data0 = generate_grand_product_data::<_, Pcs>(&pp.pcs_param, &pp.S0_poly.as_ref().unwrap()[i], &S0_comm[i], &rp0, &tau, &gamma, transcript)?;
        r_prime.push(rp0); // rp0 is the point [b_rand, r_prime_rand...] from sumcheck
        S0_data.push(data0);

        // S Sumcheck
        let (_proof2, rp2, _eval2) = run_grand_product_sumcheck::<_, Pcs>(&pp.S_poly.as_ref().unwrap()[i], pp.l, transcript)?;
        let data2 = generate_grand_product_data::<_, Pcs>(&pp.pcs_param, &pp.S_poly.as_ref().unwrap()[i], &S_comm[i], &rp2, &tau, &gamma, transcript)?;
        r_prime2.push(rp2);
        S_data.push(data2);

        // RS Sumcheck
        let (_proof3, rp3, _eval3) = run_grand_product_sumcheck::<_, Pcs>(&pp.RS_poly.as_ref().unwrap()[i], pp.logm, transcript)?;
        let data3 = generate_grand_product_data::<_, Pcs>(&pp.pcs_param, &pp.RS_poly.as_ref().unwrap()[i], &RS_comm[i], &rp3, &tau, &gamma, transcript)?;
        r_prime3.push(rp3);
        RS_data.push(data3);

        // WS Sumcheck
        let (_proof4, rp4, _eval4) = run_grand_product_sumcheck::<_, Pcs>(&pp.WS_poly.as_ref().unwrap()[i], pp.logm, transcript)?;
        let data4 = generate_grand_product_data::<_, Pcs>(&pp.pcs_param, &pp.WS_poly.as_ref().unwrap()[i], &WS_comm[i], &rp4, &tau, &gamma, transcript)?;
        r_prime4.push(rp4);
        WS_data.push(data4);
    }

    // --- Final Openings ---
    // These points r_prime2, r_prime3 etc. are the full points [b, r'_...] from sumcheck.
    // The openings need to be at the 'x' part of these points.
    // Let's adjust the points used for opening.

    // Create batch opening structures
    let mut batch_open_polys = Vec::new();
    let mut batch_open_commits = Vec::new();
    let mut batch_open_points = Vec::new();
    let mut batch_open_evals = Vec::new();

    for i in 0..pp.alpha {
        let chunk_idx = i / pp.k;
        // Points for evaluation should exclude the first element 'b' from sumcheck point
        // r_primeX has size num_vars_base + 1. We need the last num_vars_base elements.
        let rp2_x = &r_prime2[i][1..]; // Point for final_poly (size l)
        let rp3_x = &r_prime3[i][1..]; // Point for E, dim, read (size logm)

        // Evaluate E, dim, read at rp3_x and add to batch
        let eval_e = pp.E_poly[i].evaluate(rp3_x);
        E_eval2.push(eval_e);
        batch_open_polys.push(&pp.E_poly[i]);
        batch_open_commits.push(&pp.E_comm[i]);
        batch_open_points.push(rp3_x.to_vec());
        batch_open_evals.push(eval_e);

        let eval_dim = pp.dim_poly[chunk_idx].evaluate(rp3_x);
        dim_eval.push(eval_dim);
        batch_open_polys.push(&pp.dim_poly[chunk_idx]);
        batch_open_commits.push(&pp.dim_comm[chunk_idx]);
        batch_open_points.push(rp3_x.to_vec());
        batch_open_evals.push(eval_dim);

        let eval_read = pp.read_poly[i].evaluate(rp3_x);
        read_eval.push(eval_read);
        batch_open_polys.push(&pp.read_poly[i]);
        batch_open_commits.push(&pp.read_comm[i]);
        batch_open_points.push(rp3_x.to_vec());
        batch_open_evals.push(eval_read);

        // Evaluate final at rp2_x and add to batch
        let eval_final = pp.final_poly[i].evaluate(rp2_x);
        final_eval.push(eval_final);
        batch_open_polys.push(&pp.final_poly[i]);
        batch_open_commits.push(&pp.final_comm[i]);
        batch_open_points.push(rp2_x.to_vec());
        batch_open_evals.push(eval_final);
    }

    // Perform batch opening
    Pcs::batch_open(
        &pp.pcs_param,
        &batch_open_polys,
        &batch_open_commits,
        &batch_open_points,
        &batch_open_evals,
        transcript,
    )?;


    let msg5 = Message5 {
        // gp_sumcheck_proofs, // Removed: proofs implicitly in transcript
        r_prime, r_prime2, r_prime3, r_prime4, // These are the full points [b, r'_...]
        S0_data, S_data, RS_data, WS_data,
        E_eval2, dim_eval, read_eval, final_eval,
        placeholder_proofs: PhantomData,
    };
    transcript_f.write_serializable(b"msg5", &msg5)?; // Prover writes msg5

    Ok(LassoProof { msg1, msg2, msg3, msg4, msg5 })
}


// Helper for grand product polynomial f(b, x)
// f(0, x) = H(tuple(x))
// f(1, x) = f(x, 0) * f(x, 1) -> This is the recursive relation, not the direct evaluation.
// Let's follow Spartan approach (Figure 5):
// f_evals = [H(t_0), ..., H(t_{N-1}), P_0, P_1, ..., P_{N/2-1}, ..., Root, 0] where N = 2^num_vars_base
// Total size is 2N = 2^(num_vars_base + 1)
// f(0, x) corresponds to indices 0 to N-1. H(t_x) is stored at index x.
// f(1, x) corresponds to indices N to 2N-1. Product P_y is stored at N + y. Root is at 2N-2. 2N-1 is 0.
fn build_grand_product_poly<F: PrimeField + Hash>( // Added Hash bound for tuple hashing
    tuples: &[(F, F, F)],
    num_vars_base: usize,
    tau: &F,
    gamma: &F,
) -> Result<MultilinearPolynomial<F>, Error> {
    let n_pow = 1 << num_vars_base; // N = 2^num_vars_base
    let total_size = n_pow * 2; // 2N = 2^(num_vars_base + 1)
    let num_provided_tuples = tuples.len();

    // If num_vars_base is 0, n_pow is 1, total_size is 2.
    if num_vars_base == 0 {
        if num_provided_tuples > 1 {
            return Err(Error::InvalidSnark(format!(
                "build_grand_product_poly (0-var): Expected 0 or 1 tuple, got {}", num_provided_tuples
            )));
        }
        let mut f_evals = Vec::with_capacity(total_size);
        // f(0, empty) = H(tuple_0) or H(default) if no tuple provided
        let h_val = if num_provided_tuples == 1 {
            hash_tuple(tuples[0], gamma, tau)
        } else {
            // What's the default hash for padding/empty? Let's hash (0,0,0)
            hash_tuple((F::ZERO, F::ZERO, F::ZERO), gamma, tau)
        };
        f_evals.push(h_val);
        // f(1, empty) = 0 according to paper
        f_evals.push(F::ZERO);
        return Ok(MultilinearPolynomial::new(f_evals));
    }

    // Must have n_pow tuples provided for non-zero vars case.
    if num_provided_tuples != n_pow {
        return Err(Error::InvalidSnark(format!(
            "build_grand_product_poly: Expected {} tuples for {} vars, got {}",
            n_pow, num_vars_base, num_provided_tuples
        )));
    }

    let mut f_evals = vec![F::ZERO; total_size];

    // Compute f(0, x) = H(tuple(x)) for x in hypercube {0,1}^num_vars_base
    // Store H(t_x) at index x
    for i in 0..n_pow {
        f_evals[i] = hash_tuple(tuples[i], gamma, tau);
    }

    // Compute f(1, x) using recursive product based on paper's Figure 5 / Spartan approach
    // f_evals[N..2N] holds the product tree values.
    // Level d products (size N/2^d) start at index N + N/2 + ... + N/2^(d-1).
    // Level 1 (leaves): f(0, x) values are already in f_evals[0..N]
    // Level 2 (products of pairs): Compute P_i = H(t_{2i}) * H(t_{2i+1})
    // Store P_i at index N + i (for i = 0 to N/2 - 1)
    let mut current_layer_start_idx = 0;
    let mut next_layer_write_idx = n_pow; // Start writing products at index N
    let mut layer_size = n_pow;

    while layer_size > 1 {
        let next_layer_size = layer_size / 2;
        for i in 0..next_layer_size {
            let left_child_idx = current_layer_start_idx + 2 * i;
            let right_child_idx = left_child_idx + 1;
            f_evals[next_layer_write_idx + i] = f_evals[left_child_idx] * f_evals[right_child_idx];
        }
        current_layer_start_idx = next_layer_write_idx; // Next layer reads from where we just wrote
        next_layer_write_idx += next_layer_size; // Update write index for the *next* next layer
        layer_size = next_layer_size;
    }
    // The root is now at index 2N - 2.
    // The value at 2N - 1 should remain 0 (padding / convention).

     // Final check for size
     if f_evals.len() != total_size {
          return Err(Error::InvalidSnark(format!("Grand product poly final size mismatch: {} vs {}", f_evals.len(), total_size)));
     }

    Ok(MultilinearPolynomial::new(f_evals))
}


// Helper to run sumcheck for grand product
// Runs sumcheck on F(b, x) = (1-b)f(0,x) + b*f(1,x) which is equivalent to f(b,x) directly.
// Claimed sum is 0.
fn run_grand_product_sumcheck<F, Pcs>(
    f_poly: &MultilinearPolynomial<F>, // This is f(b, x)
    num_vars_base: usize,
    transcript: &mut impl FieldTranscriptWrite<F>,
) -> Result<(Point<F, Pcs::Polynomial>, F), Error> // Return (point, final_eval)
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
{
    let num_vars = num_vars_base + 1; // Sumcheck over b and x
    // Define the expression for F(b, x) - it's just the polynomial itself
    let expr = Expression::Polynomial(Query::new(0, Rotation::cur()));
    let polys = [f_poly]; // Pass the f(b,x) polynomial

    // Define the virtual polynomial for sumcheck
    let virtual_poly_gp = VirtualPolynomial::new(
        &expr,
        &polys,
        &[], // No external challenges
        &[], // No eq_xy terms
    );

    // Run sumcheck
    let (gp_proof, gp_point, gp_final_eval) = ClassicSumCheck::<CoefficientsProver<F>, F>::prove(
        &(), // prover param for sumcheck
        num_vars, // Num vars = base + 1 (for b)
        virtual_poly_gp,
        F::ZERO, // claimed sum is zero
        transcript, // Pass FieldTranscriptWrite
    )?;

    Ok((gp_point, gp_final_eval))
}


// Helper to generate GrandProductData and run Pcs::open
// The point r_prime_full is the full point [b', r'_1, ..., r'_n] from sumcheck.
// We need to evaluate f at points related to r_prime_full:
// f(0, r') = f(0, r'_1, ..., r'_n)
// f(1, r') = f(1, r'_1, ..., r'_n)
// f(r', 0) = f(b', r'_1, ..., r'_{n-1}, 0)
// f(r', 1) = f(b', r'_1, ..., r'_{n-1}, 1)
// product = f(1, 1, ..., 1, 0) <-- Check this definition. Spartan Fig 5 uses Root_P = f(1,0,...,0). Let's use Root_P.
// The "product" required by the verifier check H0*Hw = Hf*Hr refers to the overall product
// of the leaves, which is Root_P = f(1, 0, ..., 0).
fn generate_grand_product_data<F, Pcs>(
    pcs_param: &Pcs::ProverParam,
    f_poly: &Pcs::Polynomial, // This is f(b, x_1, ..., x_n)
    f_comm: &Pcs::Commitment, // Commitment to f(b, x)
    r_prime_full: &[F], // Point [b', r'_1, ..., r'_n] from sumcheck (size n+1)
    tau: &F, // Needed for potential hash checks? Not directly used here.
    gamma: &F, // Needed for potential hash checks? Not directly used here.
    transcript: &mut impl FieldTranscriptWrite<F>,
) -> Result<GrandProductData<F, Pcs>, Error>
where
    F: PrimeField + Hash,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    Pcs::Commitment: Default + Clone + Debug + PartialEq + Send + Sync + AsRef<[Pcs::CommitmentChunk]>,
    Pcs::CommitmentChunk: Clone + Debug + Default + Send + Sync,
    Pcs::ProverParam: Clone + Debug + Send + Sync,
{
    let num_vars = r_prime_full.len(); // n + 1
    if num_vars == 0 {
         return Err(Error::InvalidSnark("generate_grand_product_data called with 0 vars".to_string()));
    }
    let num_vars_base = num_vars - 1; // n
    let r_prime_x = &r_prime_full[1..]; // r' = [r'_1, ..., r'_n]

    // Points for evaluation
    let point_0r: Point<F, Pcs::Polynomial> = [[F::ZERO].as_slice(), r_prime_x].concat();
    let point_1r: Point<F, Pcs::Polynomial> = [[F::ONE].as_slice(), r_prime_x].concat();
    let point_r0: Point<F, Pcs::Polynomial> = if num_vars_base > 0 {
        [r_prime_full[0..num_vars_base].as_ref(), &[F::ZERO]].concat() // [b', r'_1, ..., r'_{n-1}, 0]
    } else {
        vec![F::ZERO] // Special case n=0: point is [0]
    };
    let point_r1: Point<F, Pcs::Polynomial> = if num_vars_base > 0 {
        [r_prime_full[0..num_vars_base].as_ref(), &[F::ONE]].concat() // [b', r'_1, ..., r'_{n-1}, 1]
    } else {
        vec![F::ONE] // Special case n=0: point is [1]
    };

    // Point for Root_P: f(1, 0, ..., 0)
    let mut point_root_p = vec![F::ZERO; num_vars];
    point_root_p[0] = F::ONE; // b = 1

    // Evaluate f at these points
    let f_0_r = f_poly.evaluate(&point_0r);
    let f_1_r = f_poly.evaluate(&point_1r);
    let f_r_0 = f_poly.evaluate(&point_r0);
    let f_r_1 = f_poly.evaluate(&point_r1);
    let product = f_poly.evaluate(&point_root_p); // Root_P = f(1, 0..0)

    // Open commitment at these points (use batch open?)
    // For now, open individually as verifier also batch verifies everything together later.
    // This might add overhead. Consider collecting all openings for a single batch call.
    Pcs::open(pcs_param, f_poly, f_comm, &point_0r, &f_0_r, transcript)?;
    Pcs::open(pcs_param, f_poly, f_comm, &point_1r, &f_1_r, transcript)?;
    Pcs::open(pcs_param, f_poly, f_comm, &point_r0, &f_r_0, transcript)?;
    Pcs::open(pcs_param, f_poly, f_comm, &point_r1, &f_r_1, transcript)?;
    Pcs::open(pcs_param, f_poly, f_comm, &point_root_p, &product, transcript)?;


    Ok(GrandProductData {
        f_0_r,
        f_1_r,
        f_r_0,
        f_r_1,
        product, // Store Root_P = f(1, 0..0) here
        placeholder_proofs: PhantomData, // PCS proofs are written to transcript
    })
}


// Add temporary storage for grand product polys in ProverParam
#[derive(Clone, Debug)]
pub struct LassoProverParam<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) l: usize,
    pub(crate) c: usize,
    pub(crate) k: usize,
    pub(crate) alpha: usize,
    pub(crate) logm: usize,
    pub(crate) subtables: Vec<Vec<F>>,
    pub(crate) a_poly: MultilinearPolynomial<F>,
    pub(crate) dim_poly: Vec<MultilinearPolynomial<F>>,
    pub(crate) E_poly: Vec<MultilinearPolynomial<F>>,
    pub(crate) read_poly: Vec<MultilinearPolynomial<F>>,
    pub(crate) write_poly: Vec<MultilinearPolynomial<F>>, // Needed for WS polynomial
    pub(crate) final_poly: Vec<MultilinearPolynomial<F>>,
    pub(crate) a_comm: Pcs::Commitment,
    pub(crate) dim_comm: Vec<Pcs::Commitment>,
    pub(crate) E_comm: Vec<Pcs::Commitment>,
    pub(crate) read_comm: Vec<Pcs::Commitment>,
    pub(crate) final_comm: Vec<Pcs::Commitment>,
    pub(crate) pcs_param: Pcs::ProverParam,
    // Temporary storage for grand product polys needed across rounds
    pub(crate) S0_poly: Option<Vec<MultilinearPolynomial<F>>>,
    pub(crate) S_poly: Option<Vec<MultilinearPolynomial<F>>>,
    pub(crate) RS_poly: Option<Vec<MultilinearPolynomial<F>>>,
    pub(crate) WS_poly: Option<Vec<MultilinearPolynomial<F>>>,
}

