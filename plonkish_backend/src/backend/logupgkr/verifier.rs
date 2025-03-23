use crate::{
    backend::logupgkr::LogupGkrVerifierParam,
    piop::gkr::fractional_sum_check::verify_fractional_sum_check,
    util::{
        arithmetic::PrimeField,
        transcript::FieldTranscriptRead,
    },
    Error,
};

pub fn verify<F>(
    vp: LogupGkrVerifierParam<F>,
    transcript: &mut impl FieldTranscriptRead<F>,
) -> Result<(), Error>
where
    F: PrimeField,
{
    // Run the fractional sum check verification
    verify_fractional_sum_check::<F>(
        vp.num_vars,
        vp.p_0s,
        vp.q_0s,
        transcript,
    )?;
    
    Ok(())
}