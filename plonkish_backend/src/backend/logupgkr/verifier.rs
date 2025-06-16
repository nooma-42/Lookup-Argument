// src/logupgkr/verifier.rs

use crate::{
    backend::logupgkr::{VerifierParam, util},
    piop::gkr::fractional_sum_check::verify_fractional_sum_check,
    pcs::{multilinear::MultilinearKzg, PolynomialCommitmentScheme, Evaluation, CommitmentChunk},
    util::{
        transcript::{FieldTranscript, TranscriptRead},
    },
    Error,
};
use halo2_curves::bn256::{Bn256, Fr, G1Affine};

type Pcs = MultilinearKzg<Bn256>;

pub fn verify(
    vp: &VerifierParam,
    transcript: &mut (impl TranscriptRead<G1Affine, Fr>
                  + FieldTranscript<Fr>),
) -> Result<(), Error> {
    transcript.common_commitments(vp.m_comm.as_ref())?;
    transcript.common_commitments(vp.t_comm.as_ref())?;
    for w_comm in &vp.w_comms {
        transcript.common_commitments(w_comm.as_ref())?;
    }

    let a = transcript.squeeze_challenge();

    let (p_at_x_claimed, q_at_x_claimed, x) = verify_fractional_sum_check(
        vp.num_vars, [None].iter().cloned(), [None].iter().cloned(),
        transcript,
    )?;
    
    let p_at_x_from_proof = transcript.read_field_element()?;
    let q_at_x_from_proof = transcript.read_field_element()?;

    if p_at_x_from_proof != p_at_x_claimed[0] || q_at_x_from_proof != q_at_x_claimed[0] {
         return Err(Error::InvalidSnark("IOP claimed values mismatch with proof".to_string()));
    }

    let x_row = &x[..vp.num_vars_row];
    
    let m_at_x_row = transcript.read_field_element()?;
    let t_at_x_row = transcript.read_field_element()?;
    let w_at_x_rows = transcript.read_field_elements(vp.w_comms.len())?;

    let evals_to_verify = [
        Evaluation::new(0, 0, m_at_x_row),
        Evaluation::new(1, 0, t_at_x_row),
    ].into_iter().chain(
        w_at_x_rows.iter().enumerate().map(|(i, eval)| Evaluation::new(2 + i, 0, *eval))
    ).collect::<Vec<_>>();
    
    Pcs::batch_verify(
        &vp.pcs_vp,
        [Some(&vp.m_comm), Some(&vp.t_comm)].into_iter().flatten().chain(vp.w_comms.iter()),
        &[x_row.to_vec()],
        &evals_to_verify,
        transcript
    )?;

    let p_at_x_expected = util::evaluate_p_at_x(&x, m_at_x_row, vp.num_vars_row);
    let q_at_x_expected = util::evaluate_q_at_x(&x, t_at_x_row, &w_at_x_rows, a, vp.num_vars_row);
    
    if p_at_x_from_proof != p_at_x_expected || q_at_x_from_proof != q_at_x_expected {
        return Err(Error::InvalidSnark("Final evaluation check failed".to_string()));
    }

    Ok(())
}