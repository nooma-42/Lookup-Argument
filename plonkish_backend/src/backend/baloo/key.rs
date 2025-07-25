use halo2_curves::bn256::{Fr, G1Affine};
use crate::pcs::univariate::{UnivariateKzgProverParam, UnivariateKzgVerifierParam, UnivariateKzgCommitment, UnivariateKzgParam};

/// Contains preprocessed data for the Baloo prover, specific to a given lookup table.
/// This structure moves all O(t log t) computations from the online proving phase
/// to the offline preprocessing phase.
#[derive(Debug, Clone)]
pub struct BalooProverKey {
    /// Universal KZG prover parameters
    pub pp: UnivariateKzgProverParam<halo2_curves::bn256::Bn256>,
    /// Universal KZG parameters (needed for G2 commitments)
    pub param: UnivariateKzgParam<halo2_curves::bn256::Bn256>,
    /// Commitment to the complete lookup table polynomial T(X)
    pub table_comm: UnivariateKzgCommitment<G1Affine>,
    /// Commitment to the vanishing polynomial of the full group Z_H(X) = X^t - 1
    pub z_h_comm: UnivariateKzgCommitment<G1Affine>,
    /// Precomputed witness commitments for each table element
    /// table_element_proofs[i] = [(T(X) - table[i]) / (X - ω^i)]_1
    pub table_element_proofs: Vec<G1Affine>,
    /// Precomputed witness commitments for vanishing polynomial evaluation
    /// subgroup_element_proofs[i] = [Z_H(X) / (X - ω^i)]_1
    pub subgroup_element_proofs: Vec<G1Affine>,
    /// The original lookup table (needed for indexing during proving)
    pub table: Vec<Fr>,
    /// Table size
    pub t: usize,
    /// Maximum degree bound for polynomials
    pub d: usize,
}

/// Contains preprocessed data for the Baloo verifier, specific to a given lookup table.
#[derive(Debug, Clone)]
pub struct BalooVerifierKey {
    /// Universal KZG verifier parameters
    pub vp: UnivariateKzgVerifierParam<halo2_curves::bn256::Bn256>,
    /// Universal KZG parameters (needed for G2 operations)
    pub param: UnivariateKzgParam<halo2_curves::bn256::Bn256>,
    /// Commitment to the complete lookup table polynomial T(X)
    pub table_comm: UnivariateKzgCommitment<G1Affine>,
    /// Commitment to the vanishing polynomial of the full group Z_H(X) = X^t - 1
    pub z_h_comm: UnivariateKzgCommitment<G1Affine>,
    /// Table size
    pub t: usize,
    /// Maximum degree bound for polynomials
    pub d: usize,
}