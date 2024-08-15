use crate::{
    util::arithmetic::PrimeField,
    pcs::PolynomialCommitmentScheme,
    backend::baloo::{
        preprocessor::preprocess,
        // prover::prove,
        // verifier::verify,
    }
};

pub mod preprocessor;
pub mod prover;
pub mod verifier;

#[derive(Clone, Debug)]
pub struct BalooProverParam //<F, PCS>
// where
//     F: PrimeField,
//     PCS: PolynomialCommitmentScheme<F>,
{
    pub(crate) num_vars: usize,
}

#[derive(Clone, Debug)]
pub struct BalooVerifierParam<F, PCS>
where
    F: PrimeField,
    PCS: PolynomialCommitmentScheme<F>,
{
    // [z_H_comm_1, t_comm_1]
    pub(crate) preprocess_comms: Vec<PCS::Commitment>,
}
