#![allow(non_snake_case)]

use crate::pcs::{univariate::{UnivariateKzg, UnivariateKzgCommitment, UnivariateKzgVerifierParam}, PolynomialCommitmentScheme};
use crate::util::arithmetic::{Field, MultiMillerLoop, WithSmallOrderMulGroup};
use crate::util::transcript::{FieldTranscript, G2TranscriptRead, TranscriptRead};
use crate::Error;
use halo2_curves::ff::PrimeField;
use halo2_curves::group::prime::PrimeCurveAffine;
use halo2_curves::pairing::{Engine, PairingCurveAffine};
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::ops::{Add, Sub};

use super::{util::*, CaulkVerifierParam};

pub fn verify<
    M: MultiMillerLoop + Engine + std::marker::Sync + std::marker::Send + std::fmt::Debug,
>(
    vp: &CaulkVerifierParam<M>,
    transcript: &mut (impl TranscriptRead<M::G1Affine, M::Scalar>
                  + G2TranscriptRead<M::G2Affine, M::Scalar>
                  + FieldTranscript<M::Scalar>),
) -> Result<(), Error>
where
    M::Scalar: Serialize + DeserializeOwned + PrimeField + WithSmallOrderMulGroup<3>,
    M::G1Affine: Serialize + DeserializeOwned + Clone,
    M::G2Affine: Serialize + DeserializeOwned + Clone,
{
    let g1_C = transcript.read_commitment()?;
    let cm = transcript.read_commitment()?;
    let g2_H1 = transcript.read_commitment_g2()?;
    let g1_u = transcript.read_commitment()?;
    let g1_C_I = transcript.read_commitment()?;
    let g1_Z_I = transcript.read_commitment()?;
    let chi = transcript.squeeze_challenge();
    let g1_H2 = transcript.read_commitment()?;
    let alpha = transcript.squeeze_challenge();
    let g1_p1 = transcript.read_commitment()?;
    let g1_p2 = transcript.read_commitment()?;

    let v1 = transcript.read_field_element()?;
    UnivariateKzg::<M>::verify(
        &vp.kzg_vp,
        &UnivariateKzgCommitment(g1_u),
        &alpha,
        &v1,
        transcript,
    )?;

    let v2 = transcript.read_field_element()?;
    UnivariateKzg::<M>::verify(
        &vp.kzg_vp,
        &UnivariateKzgCommitment(g1_p1),
        &alpha,
        &v2,
        transcript,
    )?;

    UnivariateKzg::<M>::verify(
        &vp.kzg_vp,
        &UnivariateKzgCommitment(g1_p2),
        &alpha,
        &<M as Engine>::Scalar::ZERO,
        transcript,
    )?;

    let z_v_m_poly = get_vanishing_poly::<M::Scalar>(vp.m);
    let _z_v_m_at_alpha = z_v_m_poly.evaluate(&alpha);

    let g1_c_proj = M::G1::from(g1_C);
    let g1_c_i_proj = M::G1::from(g1_C_I);
    let diff_proj = g1_c_proj - g1_c_i_proj;
    let C_minus_C_I: M::G1Affine = diff_proj.into();
    let g2_gen = M::G2Affine::generator();
    let lhs1 = M::pairing(&C_minus_C_I, &g2_gen);
    let rhs1 = M::pairing(&g1_Z_I.into(), &g2_H1.into());
    if lhs1 != rhs1 {
        return Err(crate::Error::InvalidSnark("Pairing check 1 failed".to_string()));
    }

    Ok(())
} 