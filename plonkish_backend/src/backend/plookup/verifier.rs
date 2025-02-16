use std::ops::Div;

use crate::{
    backend::plookup::{PlookupInfo, PlookupVerifierParam},
    pcs::PolynomialCommitmentScheme,
    poly::univariate::UnivariatePolynomial,
    util::{arithmetic::PrimeField, transcript::TranscriptWrite},
    halo2_curves::ff::WithSmallOrderMulGroup,
    Error,
};

pub(super) fn prove<
    F: PrimeField + WithSmallOrderMulGroup<3>, 
    Pcs: PolynomialCommitmentScheme<F, Polynomial = UnivariatePolynomial<F>>
>(
    vp: PlookupVerifierParam<F, Pcs>,
    transcript: &impl TranscriptWrite<Pcs::CommitmentChunk, F>,
) -> bool {
    true
}