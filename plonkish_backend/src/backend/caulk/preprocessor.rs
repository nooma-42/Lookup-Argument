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
    // Need to commit polynomials up to degree N-1 (for C(X)) and potentially others related to m.
    // Crucially, the prover computes compositions like c_I(u(X)), which drastically increase degree.
    // deg(c_I) is roughly N/2 + blinding_degree (e.g., 3)
    // deg(u) is m-1
    // deg(c_I(u(X))) is roughly (N/2 + 3) * (m-1)
    // H2_poly = (z_I(u(X)) + (c_I(u(X)) - phi(X))*chi) / Z_v_m(X)
    // Degree of H2 is roughly deg(c_I(u(X))) - deg(Z_v_m) = (N/2+3)*(m-1) - m
    let max_degree_composition = ((N / 2).max(1) + 3) * (m.max(1) - 1); // Approximate max degree from composition
    let max_degree_h2 = max_degree_composition.saturating_sub(m);
    let max_degree_commit = max(N.max(1) - 1, max_degree_h2); // Consider C(X), phi(X) (m-1), H1 (~N/2), H2, p1, p2, u
    let poly_size = (max_degree_commit + 1).next_power_of_two(); // Degree bound for KZG setup

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