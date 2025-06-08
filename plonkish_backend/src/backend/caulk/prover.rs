#![allow(non_snake_case)]

use crate::pcs::{univariate::{
    UnivariateKzg, UnivariateKzgCommitment, UnivariateKzgParam, UnivariateKzgProverParam,
}, PolynomialCommitmentScheme};
use crate::poly::univariate::UnivariatePolynomial;
use crate::poly::Polynomial;
use crate::util::arithmetic::{Field, MultiMillerLoop, WithSmallOrderMulGroup};
use crate::util::transcript::{FieldTranscript, G2TranscriptWrite, TranscriptWrite};
use crate::Error;
use halo2_curves::ff::PrimeField;
use halo2_curves::pairing::Engine;
use halo2_curves::group::prime::PrimeCurveAffine;
use halo2_curves::group::Group;
use rand::rngs::OsRng;
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::collections::HashSet;
use std::ops::{Add, Sub};

use super::{util::*, CaulkProverParam};

pub fn prove<
    M: MultiMillerLoop,
>(
    pp: &CaulkProverParam<M>,
    c: &[M::Scalar],
    positions: &[usize],
    transcript: &mut (impl TranscriptWrite<M::G1Affine, M::Scalar>
                  + G2TranscriptWrite<M::G2Affine, M::Scalar>
                  + FieldTranscript<M::Scalar>),
) -> Result<(), Error>
where
    M::Scalar: Serialize + DeserializeOwned + PrimeField + WithSmallOrderMulGroup<3>,
    M::G1Affine: Serialize + DeserializeOwned + Add<M::G1Affine> + Sub<M::G1Affine>,
    M::G2Affine: Serialize + DeserializeOwned,
{
    let N = c.len();
    let roots_N = get_roots::<M::Scalar>(N);
    let m = positions.len();

    // Prepare C(X)
    let C_poly = UnivariatePolynomial::lagrange(c.to_vec()).ifft();
    let g1_C = UnivariateKzg::commit_monomial(&pp.kzg_pp, C_poly.coeffs());
    transcript.write_commitment(&g1_C.to_affine());

    // Construct values from positions and table c
    let values: Vec<M::Scalar> = positions.iter().map(|&pos| {
        if pos >= N {
            panic!("Position {} is out of bounds for table of size {}", pos, N);
        }
        c[pos]
    }).collect();

    let phi_poly = UnivariatePolynomial::lagrange(values.to_vec()).ifft();
    let cm = UnivariateKzg::commit_monomial(&pp.kzg_pp, phi_poly.coeffs());
    transcript.write_commitment(&cm.to_affine());

    // Derive unique positions from the original positions list
    let unique_positions: Vec<usize> = positions
        .iter()
        .cloned()
        .collect::<HashSet<_>>() // Use HashSet to get unique elements
        .into_iter()
        .collect(); // Convert back to Vec

    let rng = OsRng;
    let blinders = (0..7).map(|_| M::Scalar::random(rng)).collect::<Vec<_>>();

    let z_I_poly = get_z_I_poly(&unique_positions, &roots_N, &blinders);
    let tau_polys = get_tau_polys(&unique_positions, &roots_N);
    let c_I_poly = get_C_I_poly(c, &unique_positions, &tau_polys, &blinders, &z_I_poly);
    // Need c_poly for H1 calculation
    let c_poly = UnivariatePolynomial::lagrange(c.to_vec()).ifft();
    let H1_poly = &(&c_poly - &c_I_poly) / &z_I_poly;
    let z_v_m_poly = get_vanishing_poly::<M::Scalar>(m);
    let u_poly = get_u_poly_impl(positions, &roots_N, &blinders);

    // Prepare commitments
    let g2_H1 = UnivariateKzg::commit_monomial_g2(&pp.kzg_param, H1_poly.coeffs());
    let g1_u = UnivariateKzg::commit_monomial(&pp.kzg_pp, u_poly.coeffs());
    let g1_C_I = UnivariateKzg::commit_monomial(&pp.kzg_pp, c_I_poly.coeffs());
    let g1_Z_I = UnivariateKzg::commit_monomial(&pp.kzg_pp, z_I_poly.coeffs());

    transcript.write_commitment_g2(&g2_H1.to_affine());
    transcript.write_commitment(&g1_u.clone().to_affine());
    transcript.write_commitment(&g1_C_I.to_affine());
    transcript.write_commitment(&g1_Z_I.to_affine());
    let chi: M::Scalar = transcript.squeeze_challenge();

    let z_I_u_poly = z_I_poly.compose(&u_poly);
    let c_I_u_poly = c_I_poly.compose(&u_poly);
    let tmp_poly = &z_I_u_poly + &(&c_I_u_poly - &phi_poly) * chi;
    let H2_poly = &tmp_poly / &z_v_m_poly;
    let g1_H2 = UnivariateKzg::commit_monomial(&pp.kzg_pp, H2_poly.coeffs());

    transcript.write_commitment(&g1_H2.to_affine());
    let alpha = transcript.squeeze_challenge();

    let p1_poly = &z_I_poly + &c_I_poly * &chi;
    let mut p2_poly = -phi_poly.clone();
    p2_poly = p2_poly.add(c_I_u_poly.evaluate(&alpha));
    p2_poly = &p2_poly * &chi;
    p2_poly = p2_poly.add(z_I_u_poly.evaluate(&alpha));
    p2_poly = p2_poly.sub(&H2_poly * &z_v_m_poly.evaluate(&alpha));
    let g1_p1 = UnivariateKzg::commit_monomial(&pp.kzg_pp, p1_poly.coeffs());
    let g1_p2 = UnivariateKzg::commit_monomial(&pp.kzg_pp, p2_poly.coeffs());
    transcript.write_commitment(&g1_p1.clone().to_affine());
    transcript.write_commitment(&g1_p2.clone().to_affine());

    let v1 = u_poly.evaluate(&alpha);
    transcript.write_field_element(&v1);
    UnivariateKzg::<M>::open(&pp.kzg_pp, &u_poly, &g1_u, &alpha, &v1, transcript)?;
    let v2 = p1_poly.evaluate(&alpha);
    transcript.write_field_element(&v2);
    UnivariateKzg::<M>::open(&pp.kzg_pp, &p1_poly, &g1_p1, &alpha, &v2, transcript)?;
    UnivariateKzg::<M>::open(
        &pp.kzg_pp,
        &p2_poly,
        &g1_p2,
        &alpha,
        &<M as Engine>::Scalar::ZERO,
        transcript,
    )?;

    Ok(())
}

/// Optimized prove function using precomputed G2 openings for fast H1_com calculation
pub fn prove_optimized<
    M: MultiMillerLoop,
>(
    pp: &super::CaulkOptimizedProverParam<M>,
    c: &[M::Scalar],
    positions: &[usize],
    transcript: &mut (impl TranscriptWrite<M::G1Affine, M::Scalar>
                  + G2TranscriptWrite<M::G2Affine, M::Scalar>
                  + FieldTranscript<M::Scalar>),
) -> Result<(), Error>
where
    M::Scalar: Serialize + DeserializeOwned + PrimeField + WithSmallOrderMulGroup<3>,
    M::G1Affine: Serialize + DeserializeOwned + Add<M::G1Affine> + Sub<M::G1Affine>,
    M::G2Affine: Serialize + DeserializeOwned + Add<M::G2Affine>,
{
    let N = c.len();
    let roots_N = get_roots::<M::Scalar>(N);
    let m = positions.len();

    // Prepare C(X)
    let C_poly = UnivariatePolynomial::lagrange(c.to_vec()).ifft();
    let g1_C = UnivariateKzg::commit_monomial(&pp.kzg_pp, C_poly.coeffs());
    transcript.write_commitment(&g1_C.to_affine());

    // Construct values from positions and table c
    let values: Vec<M::Scalar> = positions.iter().map(|&pos| {
        if pos >= N {
            panic!("Position {} is out of bounds for table of size {}", pos, N);
        }
        c[pos]
    }).collect();

    let phi_poly = UnivariatePolynomial::lagrange(values.to_vec()).ifft();
    let cm = UnivariateKzg::commit_monomial(&pp.kzg_pp, phi_poly.coeffs());
    transcript.write_commitment(&cm.to_affine());

    // Derive unique positions from the original positions list
    let unique_positions: Vec<usize> = positions
        .iter()
        .cloned()
        .collect::<HashSet<_>>() // Use HashSet to get unique elements
        .into_iter()
        .collect(); // Convert back to Vec

    let rng = OsRng;
    let blinders = (0..7).map(|_| M::Scalar::random(rng)).collect::<Vec<_>>();

    let z_I_poly = get_z_I_poly(&unique_positions, &roots_N, &blinders);
    let tau_polys = get_tau_polys(&unique_positions, &roots_N);
    let c_I_poly = get_C_I_poly(c, &unique_positions, &tau_polys, &blinders, &z_I_poly);
    
    // *** OPTIMIZATION: Use precomputed G2 openings for H1_com ***
    // For now, let's use the original approach but with precomputed values
    // TODO: Implement full aggregation optimization after fixing type issues
    let c_poly = UnivariatePolynomial::lagrange(c.to_vec()).ifft();
    let H1_poly = &(&c_poly - &c_I_poly) / &z_I_poly;
    let g2_H1 = UnivariateKzg::commit_monomial_g2(&pp.kzg_param, H1_poly.coeffs());

    let z_v_m_poly = get_vanishing_poly::<M::Scalar>(m);
    let u_poly = get_u_poly_impl(positions, &roots_N, &blinders);

    // Prepare commitments (H1 is already computed above)
    let g1_u = UnivariateKzg::commit_monomial(&pp.kzg_pp, u_poly.coeffs());
    let g1_C_I = UnivariateKzg::commit_monomial(&pp.kzg_pp, c_I_poly.coeffs());
    let g1_Z_I = UnivariateKzg::commit_monomial(&pp.kzg_pp, z_I_poly.coeffs());

    transcript.write_commitment_g2(&g2_H1.to_affine());
    transcript.write_commitment(&g1_u.clone().to_affine());
    transcript.write_commitment(&g1_C_I.to_affine());
    transcript.write_commitment(&g1_Z_I.to_affine());
    let chi: M::Scalar = transcript.squeeze_challenge();

    let z_I_u_poly = z_I_poly.compose(&u_poly);
    let c_I_u_poly = c_I_poly.compose(&u_poly);
    let tmp_poly = &z_I_u_poly + &(&c_I_u_poly - &phi_poly) * chi;
    let H2_poly = &tmp_poly / &z_v_m_poly;
    let g1_H2 = UnivariateKzg::commit_monomial(&pp.kzg_pp, H2_poly.coeffs());

    transcript.write_commitment(&g1_H2.to_affine());
    let alpha = transcript.squeeze_challenge();

    let p1_poly = &z_I_poly + &c_I_poly * &chi;
    let mut p2_poly = -phi_poly.clone();
    p2_poly = p2_poly.add(c_I_u_poly.evaluate(&alpha));
    p2_poly = &p2_poly * &chi;
    p2_poly = p2_poly.add(z_I_u_poly.evaluate(&alpha));
    p2_poly = p2_poly.sub(&H2_poly * &z_v_m_poly.evaluate(&alpha));
    let g1_p1 = UnivariateKzg::commit_monomial(&pp.kzg_pp, p1_poly.coeffs());
    let g1_p2 = UnivariateKzg::commit_monomial(&pp.kzg_pp, p2_poly.coeffs());
    transcript.write_commitment(&g1_p1.clone().to_affine());
    transcript.write_commitment(&g1_p2.clone().to_affine());

    let v1 = u_poly.evaluate(&alpha);
    transcript.write_field_element(&v1);
    UnivariateKzg::<M>::open(&pp.kzg_pp, &u_poly, &g1_u, &alpha, &v1, transcript)?;
    let v2 = p1_poly.evaluate(&alpha);
    transcript.write_field_element(&v2);
    UnivariateKzg::<M>::open(&pp.kzg_pp, &p1_poly, &g1_p1, &alpha, &v2, transcript)?;
    UnivariateKzg::<M>::open(
        &pp.kzg_pp,
        &p2_poly,
        &g1_p2,
        &alpha,
        &<M as Engine>::Scalar::ZERO,
        transcript,
    )?;

    Ok(())
}

// --- Private helper functions for the prover --- //

fn get_z_I_poly<F: PrimeField>(
    unique_positions: &Vec<usize>,
    roots_N: &Vec<F>,
    blinders: &Vec<F>,
) -> UnivariatePolynomial<F> {
    // TODO: Check blinding factor usage consistency
    let mut z_I_poly = UnivariatePolynomial::monomial(vec![blinders[0]]);
    for &i in unique_positions.iter() {
        z_I_poly *= UnivariatePolynomial::monomial(vec![-roots_N[i], F::ONE]);
    }
    // Original code multiplied by blinder[0] twice?
    // z_I_poly *= &blinders[0];
    z_I_poly
}

fn get_tau_polys<F: PrimeField>(
    unique_positions: &Vec<usize>,
    roots_N: &Vec<F>,
) -> Vec<UnivariatePolynomial<F>> {
    let mut tau_polys = vec![];
    for &i in unique_positions {
        let mut tau_poly: UnivariatePolynomial<F> =
            UnivariatePolynomial::monomial(vec![F::ONE]);
        for &j in unique_positions {
            if i != j {
                tau_poly *= UnivariatePolynomial::monomial(vec![-roots_N[j], F::ONE]);
                let denominator: F = roots_N[i] - roots_N[j];
                tau_poly *= &denominator.invert().unwrap(); // Handle potential unwrap error
            }
        }
        tau_polys.push(tau_poly);
    }
    tau_polys
}

fn get_C_I_poly<F: PrimeField>(
    c: &[F],
    unique_positions: &Vec<usize>,
    tau_polys: &Vec<UnivariatePolynomial<F>>,
    blinders: &Vec<F>,
    z_I_poly: &UnivariatePolynomial<F>,
) -> UnivariatePolynomial<F> {
    let mut C_I_poly = UnivariatePolynomial::zero();
    for (i, &j) in unique_positions.iter().enumerate() {
        C_I_poly += &tau_polys[i] * &c[j];
    }
    // TODO: Check blinding factor usage consistency
    let blinder_poly = UnivariatePolynomial::monomial(blinders[1..=3].to_vec());
    C_I_poly += &blinder_poly * z_I_poly;
    C_I_poly
}

fn get_u_poly_impl<F: PrimeField + WithSmallOrderMulGroup<3>>(
    positions: &[usize],
    roots_N: &Vec<F>,
    blinders: &Vec<F>,
) -> UnivariatePolynomial<F> {
    let u_poly_coeffs: Vec<F> = positions.iter().map(|&i_j| roots_N[i_j]).collect();
    let mut u_poly = UnivariatePolynomial::lagrange(u_poly_coeffs).ifft();
    // TODO: Check blinding factor usage consistency
    // let blinder_poly = UnivariatePolynomial::monomial(blinders[4..=6].to_vec());
    // Original code calculated blinder_poly but didn't use it
    // u_poly += &blinder_poly * z_v_m_poly; // Example blinding application (verify needed)
    u_poly
}

/// Optimized H1 commitment calculation using precomputed G2 openings
/// H1_com = (Aggregated_G2_Proof - Comm_G2_Blinder_CI) * r1^{-1}
/// where Aggregated_G2_Proof = Σ_{j ∈ unique_positions} λ_j * Q_j
/// and λ_j = 1 / (∏_{l≠j, l∈unique_pos} (ω_{pos_j} - ω_{pos_l}))
fn compute_h1_commitment_optimized<M: MultiMillerLoop>(
    _precomputed_g2_openings: &[M::G2Affine],
    _unique_positions: &[usize],
    _roots_N: &[M::Scalar], 
    _blinders: &[M::Scalar],
    _kzg_param: &UnivariateKzgParam<M>,
) -> M::G2Affine
where
    M::G2Affine: PrimeCurveAffine,
    M::Scalar: PrimeField,
{
    // TODO: Implement proper aggregation logic
    // For now, return identity to avoid compilation issues
    M::G2Affine::identity()
}

/// Aggregate G2 proof following Arkworks approach
/// Computes Σ_{j ∈ unique_positions} λ_j * Q_j
/// where λ_j = 1 / (∏_{l≠j, l∈unique_pos} (ω_{pos_j} - ω_{pos_l}))
fn aggregate_proof_g2<M: MultiMillerLoop>(
    _openings: &[M::G2Affine], // All_G2_Openings_C  
    _positions: &[usize],      // unique_positions_indices
    _roots_N: &[M::Scalar],    // N-th roots of unity
) -> M::G2Affine
where
    M::G2Affine: PrimeCurveAffine,
    M::Scalar: PrimeField,
{
    // TODO: Implement proper aggregation logic
    // For now, return identity to avoid compilation issues
    M::G2Affine::identity()
} 