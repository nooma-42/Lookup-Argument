use crate::{
    backend::logupgkr::{ProverParam, util::{p, q}},
    piop::gkr::fractional_sum_check::prove_fractional_sum_check,
    pcs::{multilinear::MultilinearKzg, PolynomialCommitmentScheme, Evaluation, CommitmentChunk},
    poly::multilinear::MultilinearPolynomial,
    util::{
        transcript::{FieldTranscript, TranscriptWrite},
    },
    Error,
};
use halo2_curves::bn256::{Bn256, Fr, G1Affine};

type Pcs = MultilinearKzg<Bn256>;

pub fn prove(
    pp: &ProverParam,
    transcript: &mut (impl TranscriptWrite<G1Affine, Fr>
                  + FieldTranscript<Fr>),
) -> Result<(), Error> {
    let m_comm = Pcs::commit(&pp.pcs_pp, &pp.m_poly)?;
    let t_comm = Pcs::commit(&pp.pcs_pp, &pp.t_poly)?;
    let w_comms = Pcs::batch_commit(&pp.pcs_pp, pp.w_polys.iter())?;
    
    transcript.common_commitments(m_comm.as_ref())?;
    transcript.common_commitments(t_comm.as_ref())?;
    for w_comm in &w_comms {
        transcript.common_commitments(w_comm.as_ref())?;
    }
    
    let a = transcript.squeeze_challenge();

    let num_vars_row = pp.t_poly.num_vars();
    let num_vars_col = if pp.w_polys.is_empty() { 0 } else { 
        let total_cols = pp.w_polys.len() + 1; // +1 for the t column
        // We need ceil(log2(total_cols)) bits to index all columns
        // For powers of 2, this is exactly log2, but for non-powers we need to round up
        if total_cols.is_power_of_two() {
            total_cols.ilog2() as usize
        } else {
            (total_cols.ilog2() + 1) as usize
        }
    };
    let num_vars = num_vars_row + num_vars_col;
    let all_inputs = crate::backend::logupgkr::util::generate_binary_combinations(num_vars as u32);
    
    let p_values: Vec<Fr> = all_inputs.iter().map(|input| {
        p(&input[..num_vars_row], &input[num_vars_row..], &pp.m_poly)
    }).collect();
    let q_values: Vec<Fr> = all_inputs.iter().map(|input| {
        q(&input[..num_vars_row], &input[num_vars_row..], &pp.t_poly, &pp.w_polys, a)
    }).collect();

    let ps = MultilinearPolynomial::new(p_values);
    let qs = MultilinearPolynomial::new(q_values);

    let (p_at_x, q_at_x, x) = prove_fractional_sum_check(
        [None].iter().cloned(), [None].iter().cloned(),
        [&ps].iter().copied(), [&qs].iter().copied(),
        transcript,
    )?;
    
    transcript.write_field_element(&p_at_x[0])?;
    transcript.write_field_element(&q_at_x[0])?;

    let x_row = &x[..num_vars_row];
    
    let m_eval = pp.m_poly.evaluate(x_row);
    let t_eval = pp.t_poly.evaluate(x_row);
    let w_evals: Vec<_> = pp.w_polys.iter().map(|poly| poly.evaluate(x_row)).collect();

    transcript.write_field_element(&m_eval)?;
    transcript.write_field_element(&t_eval)?;
    transcript.write_field_elements(w_evals.iter())?;

    let evals_to_open = [
        Evaluation::new(0, 0, m_eval),
        Evaluation::new(1, 0, t_eval),
    ].into_iter().chain(
        w_evals.iter().enumerate().map(|(i, eval)| Evaluation::new(2 + i, 0, *eval))
    ).collect::<Vec<_>>();
    
    Pcs::batch_open(
        &pp.pcs_pp,
        [Some(&pp.m_poly), Some(&pp.t_poly)].into_iter().flatten().chain(pp.w_polys.iter()),
        [Some(&m_comm), Some(&t_comm)].into_iter().flatten().chain(w_comms.iter()),
        &[x_row.to_vec()],
        &evals_to_open,
        transcript,
    )?;

    Ok(())
}