use crate::{
    backend::logupgkr::LogupGkrProverParam,
    piop::gkr::fractional_sum_check::prove_fractional_sum_check,
    util::{
        arithmetic::PrimeField,
        transcript::FieldTranscriptWrite,
    },
    Error,
};

pub fn prove<F>(
    pp: LogupGkrProverParam<F>,
    transcript: &mut impl FieldTranscriptWrite<F>,
) -> Result<(), Error>
where
    F: PrimeField,
{
    // Extract parameters from prover param
    let p_0s = pp.p_0s;
    let q_0s = pp.q_0s;
    let ps = pp.ps;
    let qs = pp.qs;
    
    // Run the fractional sum check protocol
    prove_fractional_sum_check::<F>(
        p_0s,
        q_0s,
        [&ps].iter().copied(),
        [&qs].iter().copied(),
        transcript,
    )?;
    
    Ok(())
}