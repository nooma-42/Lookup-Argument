use crate::pcs::univariate::{
    UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam,
};
use crate::pcs::PolynomialCommitmentScheme;
use crate::Error;
use halo2_curves::pairing::MultiMillerLoop;
use rand::rngs::OsRng;
use serde::{de::DeserializeOwned, Serialize};
use std::cmp::max;
use std::fmt::Debug;

use super::{CaulkProverParam, CaulkVerifierParam};

pub fn setup<
    M: MultiMillerLoop + Debug + Sync,
>(
    N: usize,
    m: usize,
) -> Result<(CaulkProverParam<M>, CaulkVerifierParam<M>), Error>
where
    M::Scalar: Serialize + DeserializeOwned,
    M::G1Affine: Serialize + DeserializeOwned,
    M::G2Affine: Serialize + DeserializeOwned,
{
    // N: size of the commitment vector 'c' (table size)
    // m: size of the lookup values vector 'values'
    let mut rng = OsRng;
    // Determine the required polynomial degree for KZG based on N and m.
    // Key polynomials and their degrees:
    // - C(X), phi(X): degree N-1, m-1 respectively
    // - z_I(X): degree unique_positions.len() (bounded by m) 
    // - c_I(X): degree roughly N/2 + 3 (blinding degree)
    // - u(X): degree m-1
    // - c_I(u(X)): degree (N/2 + 3) * (m-1) 
    // - z_I(u(X)): degree unique_positions.len() * (m-1), bounded by m * (m-1)
    // - H2_poly: max degree of composition polynomials
    
    // Conservative upper bound: account for all polynomial operations
    let base_degree = N.max(m);
    let composition_degree = base_degree * m; // Upper bound for compositions like c_I(u(X))
    let blinding_degree = 7; // We use 7 blinding factors
    let max_degree = composition_degree + blinding_degree;
    
    // Ensure minimum size for small examples
    let poly_size = max_degree.max(8).next_power_of_two(); // Minimum 8, then next power of 2

    // The blinding factor count seems to be 0 in the original Caulk setup/trim calls.
    let blinding_factors = 0;
    let kzg_param = UnivariateKzg::<M>::setup(poly_size, blinding_factors, &mut rng)?;
    let (kzg_pp, kzg_vp) = UnivariateKzg::<M>::trim(&kzg_param, poly_size, blinding_factors)?;
    Ok((
        CaulkProverParam {
            kzg_param: kzg_param.clone(), // Clone kzg_param as it's needed by prover
            kzg_pp,
        },
        CaulkVerifierParam { kzg_vp, m }, // Store m in VerifierParam
    ))
} 