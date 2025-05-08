// plonkish_backend/src/backend/lasso/verifier.rs
use super::{
    LassoProof, Message1, Message2, Message3, Message4, Message5,
    GrandProductData, util::{g_func_simple_range, hash_tuple},
};
use crate::{
    poly::multilinear::MultilinearPolynomial,
    util::{
        arithmetic::{Field, PrimeField},
        transcript::FieldTranscriptRead,
        expression::{CommonPolynomial, Expression, Query, Rotation},
    },
    Error,
    piop::sum_check::{self, classic::{ClassicSumCheck, CoefficientsVerifier}, VirtualPolynomial},
    pcs::PolynomialCommitmentScheme,
    pcs::{Evaluation, Point},
};
use halo2_curves::ff::WithSmallOrderMulGroup;
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;
use std::marker::PhantomData;

#[derive(Clone, Debug)]
pub struct LassoVerifierParam<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>> {
    pub(crate) l: usize,
    pub(crate) c: usize,
    pub(crate) k: usize,
    pub(crate) alpha: usize,
    pub(crate) logm: usize,
    pub(crate) subtables: Vec<Vec<F>>,
    pub(crate) a_comm: Pcs::Commitment,
    pub(crate) dim_comm: Vec<Pcs::Commitment>,
    pub(crate) E_comm: Vec<Pcs::Commitment>,
    pub(crate) read_comm: Vec<Pcs::Commitment>,
    pub(crate) final_comm: Vec<Pcs::Commitment>,
    pub(crate) pcs_param: Pcs::VerifierParam,
    pub(crate) _marker: PhantomData<Pcs::CommitmentChunk>,
}

pub fn verify<F, Pcs>(
    vp: LassoVerifierParam<F, Pcs>,
    proof: &LassoProof<F, Pcs>,
    transcript: &mut impl FieldTranscriptRead<F>,
) -> Result<(), Error>
where
    F: PrimeField + WithSmallOrderMulGroup<3> + Hash + Eq + Send + Sync,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    Pcs::Commitment: Clone + Default + Debug + PartialEq + Send + Sync + AsRef<[Pcs::CommitmentChunk]>,
    Pcs::CommitmentChunk: Clone + Debug + Default + Send + Sync,
    Pcs::VerifierParam: Clone + Debug + Send + Sync,
{
    // Collect openings for batch verification
    let mut batch_commitments = Vec::new();
    let mut batch_points_map = HashMap::new(); // Reuse points by value
    let mut batch_points = Vec::new();
    let mut batch_evaluations = Vec::new();

    let mut add_opening = |comm: &Pcs::Commitment, point: Point<F, Pcs::Polynomial>, eval: F| {
        let point_idx = *batch_points_map.entry(point.clone()).or_insert_with(|| {
            batch_points.push(point);
            batch_points.len() - 1
        });
        batch_commitments.push(comm);
        batch_evaluations.push(Evaluation::new(batch_commitments.len() - 1, point_idx, eval));
    };

    let r: Point<F, Pcs::Polynomial> = transcript.squeeze_challenges(vp.logm);

    // Add opening for 'a'
    add_opening(&proof.msg1.a_comm, r.clone(), proof.msg2.a_eval);

    let (tau, gamma) = {
        let challenges = transcript.squeeze_challenges(2);
        (challenges[0], challenges[1])
    };

    let mut g_expr = Expression::Constant(F::ZERO);
    let mut power_of_2_l = F::ONE;
    let base = F::from(2u64).pow([vp.l as u64, 0, 0, 0]);
    for j in 0..vp.c {
        let mut chunk_expr = Expression::Constant(F::ZERO);
        for i in 0..vp.k {
            let idx = j * vp.k + i;
            chunk_expr = chunk_expr + Expression::Polynomial(Query::new(idx, Rotation::cur()));
        }
        g_expr = g_expr + (chunk_expr * Expression::Constant(power_of_2_l));
        power_of_2_l *= base;
    }
    let h_expr = Expression::CommonPolynomial(CommonPolynomial::EqXY(0)) * g_expr;

    let (h_final_eval_calc, h_challenges) = ClassicSumCheck::<CoefficientsVerifier<F>, F>::verify(
        &(),
        vp.logm,
        h_expr.degree(),
        proof.msg2.a_eval,
        &[r.clone()],
        transcript,
    )?;

    if h_challenges != proof.msg3.rz {
        return Err(Error::InvalidSnark("H sumcheck challenge mismatch".to_string()));
    }

    let eq_r_rz = MultilinearPolynomial::eq_xy_evaluate(&r, &proof.msg3.rz);
    let g_at_rz = g_func_simple_range(&proof.msg3.E_eval, vp.l, vp.c, vp.k);
    let h_eval_expected = eq_r_rz * g_at_rz;
    if h_final_eval_calc != h_eval_expected {
        return Err(Error::InvalidSnark(format!(
            "H sumcheck final evaluation mismatch: expected {:?}, got {:?}",
            h_eval_expected, h_final_eval_calc
        )));
    }

    // Add E openings at rz to batch collection
    for i in 0..vp.alpha {
        add_opening(&proof.msg2.E_comm[i], proof.msg3.rz.clone(), proof.msg3.E_eval[i]);
    }

    let mut table_polys = Vec::with_capacity(vp.alpha);
    for i in 0..vp.alpha {
        table_polys.push(MultilinearPolynomial::new(vp.subtables[i].clone()));
    }

    for i in 0..vp.alpha {
        let chunk_idx = i / vp.k;
        let rp2_x = &proof.msg5.r_prime2[i][1..];
        let rp3_x = &proof.msg5.r_prime3[i][1..];
        let rp4_x = &proof.msg5.r_prime4[i][1..];
        let rp_x = &proof.msg5.r_prime[i][1..];

        let identity_eval_rp = evaluate_identity::<F>(rp_x);
        let table_eval_rp = table_polys[i].evaluate(rp_x);

        let expected_s0_f0r = hash_tuple((identity_eval_rp, table_eval_rp, F::ZERO), &gamma, &tau);
        if expected_s0_f0r != proof.msg5.S0_data[i].f_0_r {
             return Err(Error::InvalidSnark(format!("S0 hash check failed (f_0_r) for subtable {}", i)));
        }

        let identity_eval_rp2 = evaluate_identity::<F>(rp2_x);
        let table_eval_rp2 = table_polys[i].evaluate(rp2_x);
        let expected_s_f0r = hash_tuple((identity_eval_rp2, table_eval_rp2, proof.msg5.final_eval[i]), &gamma, &tau);
        if expected_s_f0r != proof.msg5.S_data[i].f_0_r {
            return Err(Error::InvalidSnark(format!("S hash check failed (f_0_r) for subtable {}", i)));
        }

        let expected_rs_f0r = hash_tuple((proof.msg5.dim_eval[i], proof.msg5.E_eval2[i], proof.msg5.read_eval[i]), &gamma, &tau);
        if expected_rs_f0r != proof.msg5.RS_data[i].f_0_r {
            return Err(Error::InvalidSnark(format!("RS hash check failed (f_0_r) for subtable {}", i)));
        }

        // Add openings for E, dim, read at rp3_x
        let rp3_x_vec = proof.msg5.r_prime3[i][1..].to_vec();
        add_opening(&vp.E_comm[i], rp3_x_vec.clone(), proof.msg5.E_eval2[i]);
        add_opening(&vp.dim_comm[chunk_idx], rp3_x_vec.clone(), proof.msg5.dim_eval[i]);
        add_opening(&vp.read_comm[i], rp3_x_vec.clone(), proof.msg5.read_eval[i]);

        // Add opening for final at rp2_x
        let rp2_x_vec = proof.msg5.r_prime2[i][1..].to_vec();
        add_opening(&vp.final_comm[i], rp2_x_vec, proof.msg5.final_eval[i]);
    }

    // Perform the batch verification
    Pcs::batch_verify(&vp.pcs_param, &batch_commitments, &batch_points, &batch_evaluations, transcript)?;

    Ok(())
}

fn evaluate_identity<F: PrimeField>(point: &[F]) -> F {
    let mut res = F::ZERO;
    let mut power_of_2 = F::ONE;
    for x_i in point {
        res += *x_i * power_of_2;
        power_of_2 = power_of_2.double();
    }
    res
}

fn verify_grand_product_sumcheck<F, Pcs>(
    num_vars_base: usize,
    r_prime_full: &Point<F, Pcs::Polynomial>,
    data: &GrandProductData<F, Pcs>,
    transcript: &mut impl FieldTranscriptRead<F>,
) -> Result<(F, Vec<F>), Error>
where
    F: PrimeField + WithSmallOrderMulGroup<3>,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    Pcs::Commitment: Clone + Default + Debug + PartialEq + Send + Sync + AsRef<[Pcs::CommitmentChunk]>,
    Pcs::CommitmentChunk: Clone + Debug + Default + Send + Sync,
    Pcs::VerifierParam: Clone + Debug + Send + Sync,
{
    let num_vars = num_vars_base + 1;
    let gp_degree = 1;

    let (gp_final_eval_calc, gp_challenges) = ClassicSumCheck::<CoefficientsVerifier<F>, F>::verify(
         &(),
         num_vars,
         gp_degree,
         F::ZERO,
         &[],
         transcript,
     )?;

    if gp_challenges != *r_prime_full {
         return Err(Error::InvalidSnark(format!("Grand product sumcheck verification failed (challenge mismatch)")));
    }

    let b_prime = r_prime_full[0];
    let f_at_r_prime_expected = data.f_0_r * (F::ONE - b_prime) + data.f_1_r * b_prime;

    if gp_final_eval_calc != f_at_r_prime_expected {
         return Err(Error::InvalidSnark(format!("Grand product sumcheck verification failed (final eval mismatch: expected {:?}, got {:?})", f_at_r_prime_expected, gp_final_eval_calc)));
    }

    Ok((gp_final_eval_calc, gp_challenges))
}
