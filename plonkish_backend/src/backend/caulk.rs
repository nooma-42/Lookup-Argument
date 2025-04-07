#![allow(non_snake_case)]
#![allow(unused)]

use crate::pcs::univariate::{
    UnivariateKzg, UnivariateKzgCommitment, UnivariateKzgParam, UnivariateKzgProverParam,
    UnivariateKzgVerifierParam,
};
use crate::pcs::PolynomialCommitmentScheme;
use crate::poly::univariate::UnivariatePolynomial;
use crate::poly::Polynomial;
use crate::util::arithmetic::{root_of_unity, Field, MultiMillerLoop, WithSmallOrderMulGroup};
use crate::util::transcript::{
    FieldTranscript, FieldTranscriptRead, FieldTranscriptWrite, G2TranscriptRead,
    G2TranscriptWrite, InMemoryTranscript, TranscriptRead, TranscriptWrite,
};
use crate::Error;
use halo2_curves::ff::PrimeField;
use halo2_curves::group::prime::PrimeCurveAffine;
use halo2_curves::pairing::Engine;
use rand::rngs::OsRng;
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::collections::HashSet;
use std::marker::PhantomData;
use std::ops::{Add, Sub};

#[derive(Clone, Debug)]
pub struct Caulk<M: MultiMillerLoop>(PhantomData<M>);

#[derive(Clone, Debug)]
pub struct CaulkParam<M: MultiMillerLoop> {
    kzg_param: UnivariateKzgParam<M>,
    kzg_pp: UnivariateKzgProverParam<M>,
    kzg_vp: UnivariateKzgVerifierParam<M>,
}

impl<M> Caulk<M>
where
    M: MultiMillerLoop,
    M::Scalar: Serialize + DeserializeOwned,
    M::G1Affine: Serialize + DeserializeOwned,
    M::G2Affine: Serialize + DeserializeOwned,
{
    pub fn setup(k: i32) -> CaulkParam<M> {
        let mut rng = OsRng;
        let poly_size = 1 << k;
        let kzg_param = UnivariateKzg::<M>::setup(poly_size, 0, &mut rng).unwrap();
        let (kzg_pp, kzg_vp) = UnivariateKzg::<M>::trim(&kzg_param, poly_size, 0).unwrap();
        CaulkParam {
            kzg_param,
            kzg_pp,
            kzg_vp,
        }
    }
}

fn is_power_of_2(n: usize) -> bool {
    n != 0 && (n & (n - 1)) == 0
}

impl<M> Caulk<M>
where
    M: MultiMillerLoop,
    M::Scalar: Serialize + DeserializeOwned + PrimeField + WithSmallOrderMulGroup<3>,
    M::G1Affine: Serialize + DeserializeOwned + Add<M::G1Affine> + Sub<M::G1Affine>,
    M::G2Affine: Serialize + DeserializeOwned,
{
    pub fn prove(
        param: &CaulkParam<M>,
        c: &[M::Scalar],
        values: &[M::Scalar],
        transcript: &mut (impl TranscriptWrite<M::G1Affine, M::Scalar>
        + G2TranscriptWrite<M::G2Affine, M::Scalar>),
    ) -> Result<(), Error> {
        let N = c.len();
        let roots_N = Self::get_roots(N);
        let m = values.len();
        let roots_m = Self::get_roots(m);

        // prepare C(X)
        let C_poly = UnivariatePolynomial::lagrange(c.to_vec()).ifft();
        let g1_C = UnivariateKzg::commit_monomial(&param.kzg_pp, C_poly.coeffs());
        transcript.write_commitment(&g1_C.to_affine());

        let phi_poly = UnivariatePolynomial::lagrange(values.to_vec()).ifft();
        let cm = UnivariateKzg::commit_monomial(&param.kzg_pp, phi_poly.coeffs());
        transcript.write_commitment(&cm.to_affine());
        let positions: Vec<usize> = values
            .iter()
            .map(|value| c.iter().position(|x| x == value).unwrap())
            .collect();
        let unique_positions: Vec<usize> = positions
            .iter()
            .cloned()
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        let rng = OsRng;
        let blinders = (0..7).map(|_| M::Scalar::random(rng)).collect::<Vec<_>>();

        let c_poly = UnivariatePolynomial::lagrange(c.to_vec()).ifft();
        let z_I_poly = Self::get_z_I_poly(&unique_positions, &roots_N, &blinders);
        let tau_polys = Self::get_tau_polys(&unique_positions, &roots_N);
        let c_I_poly = Self::get_C_I_poly(c, &unique_positions, &tau_polys, &blinders, &z_I_poly);
        let H1_poly = &(&c_poly - &c_I_poly) / &z_I_poly;
        let z_v_m_poly = Self::get_vanishing_poly(m);
        let u_poly = Self::get_u_poly(&positions, &roots_N, &blinders, &z_v_m_poly);

        // prepare commitments
        let g2_H1 = UnivariateKzg::commit_monomial_g2(&param.kzg_param, H1_poly.coeffs());
        let g1_u = UnivariateKzg::commit_monomial(&param.kzg_pp, u_poly.coeffs());
        let g1_C_I = UnivariateKzg::commit_monomial(&param.kzg_pp, c_I_poly.coeffs());
        let g1_Z_I = UnivariateKzg::commit_monomial(&param.kzg_pp, z_I_poly.coeffs());

        transcript.write_commitment_g2(&g2_H1.clone().to_affine());
        transcript.write_commitment(&g1_u.clone().to_affine());
        transcript.write_commitment(&g1_C_I.to_affine());
        transcript.write_commitment(&g1_Z_I.to_affine());
        let chi: M::Scalar = transcript.squeeze_challenge();

        let z_I_u_poly = z_I_poly.compose(&u_poly);
        let c_I_u_poly = c_I_poly.compose(&u_poly);
        let tmp_poly = &z_I_u_poly + &(&c_I_u_poly - &phi_poly) * chi;
        let H2_poly = &tmp_poly / &z_v_m_poly;
        let g1_H2 = UnivariateKzg::commit_monomial(&param.kzg_pp, H2_poly.coeffs());

        transcript.write_commitment(&g1_H2.to_affine());
        let alpha = transcript.squeeze_challenge();

        let p1_poly = &z_I_poly + &c_I_poly * &chi;
        let mut p2_poly = -phi_poly.clone();
        p2_poly = p2_poly.add(c_I_u_poly.evaluate(&alpha));
        p2_poly = &p2_poly * &chi;
        p2_poly = p2_poly.add(z_I_u_poly.evaluate(&alpha));
        p2_poly = p2_poly.sub(&H2_poly * &z_v_m_poly.evaluate(&alpha));
        let g1_p1 = UnivariateKzg::commit_monomial(&param.kzg_pp, p1_poly.coeffs());
        let g1_p2 = UnivariateKzg::commit_monomial(&param.kzg_pp, p2_poly.coeffs());
        transcript.write_commitment(&g1_p1.clone().to_affine());
        transcript.write_commitment(&g1_p2.clone().to_affine());

        let v1 = u_poly.evaluate(&alpha);
        transcript.write_field_element(&v1);
        UnivariateKzg::open(&param.kzg_pp, &u_poly, &g1_u, &alpha, &v1, transcript)?;
        let v2 = p1_poly.evaluate(&v1);
        transcript.write_field_element(&v2);
        UnivariateKzg::open(&param.kzg_pp, &p1_poly, &g1_p1, &alpha, &v2, transcript)?;
        UnivariateKzg::open(
            &param.kzg_pp,
            &p2_poly,
            &g1_p2,
            &alpha,
            &<M as Engine>::Scalar::ZERO,
            transcript,
        )?;

        Ok(())
    }

    fn get_vanishing_poly(n: usize) -> UnivariatePolynomial<M::Scalar> {
        let mut coeffs = vec![M::Scalar::ZERO; n + 1];
        coeffs[0] = -M::Scalar::ONE;
        coeffs[n] = M::Scalar::ONE;
        UnivariatePolynomial::monomial(coeffs)
    }

    fn get_roots(N: usize) -> Vec<M::Scalar> {
        assert!(is_power_of_2(N));
        let root = root_of_unity::<M::Scalar>(N.ilog2() as usize);
        (0..N).map(|i| root.pow(&[i as u64])).collect()
    }

    fn get_unique_positions(c: &[M::Scalar], values: &[M::Scalar]) -> Vec<usize> {
        let mut result = HashSet::new();
        for value in values {
            result.insert(c.iter().position(|x| x == value).unwrap());
        }
        result.into_iter().collect()
    }

    fn get_z_I_poly(
        unique_positions: &Vec<usize>,
        roots_N: &Vec<M::Scalar>,
        blinders: &Vec<M::Scalar>,
    ) -> UnivariatePolynomial<M::Scalar> {
        let mut z_I_poly = UnivariatePolynomial::monomial(vec![blinders[0]]);
        for &i in unique_positions.iter() {
            z_I_poly *= UnivariatePolynomial::monomial(vec![-roots_N[i], M::Scalar::ONE]);
        }
        z_I_poly *= &blinders[0];
        z_I_poly
    }

    fn get_tau_polys(
        unique_positions: &Vec<usize>,
        roots_N: &Vec<M::Scalar>,
    ) -> Vec<UnivariatePolynomial<M::Scalar>> {
        let mut tau_polys = vec![];
        for &i in unique_positions {
            let mut tau_poly: UnivariatePolynomial<M::Scalar> =
                UnivariatePolynomial::monomial(vec![M::Scalar::ONE]);
            for &j in unique_positions {
                if i != j {
                    tau_poly *= UnivariatePolynomial::monomial(vec![-roots_N[j], M::Scalar::ONE]);
                    let denominator: M::Scalar = roots_N[i] - roots_N[j];
                    tau_poly *= &denominator.invert().unwrap();
                }
            }
            tau_polys.push(tau_poly);
        }

        tau_polys
    }

    fn get_C_I_poly(
        c: &[M::Scalar],
        unique_positions: &Vec<usize>,
        tau_polys: &Vec<UnivariatePolynomial<M::Scalar>>,
        blinders: &Vec<M::Scalar>,
        z_I_poly: &UnivariatePolynomial<M::Scalar>,
    ) -> UnivariatePolynomial<M::Scalar> {
        let mut C_I_poly = UnivariatePolynomial::zero();
        for (i, &j) in unique_positions.iter().enumerate() {
            C_I_poly += &tau_polys[i] * &c[j];
        }
        let blinder_poly = UnivariatePolynomial::monomial(blinders[1..=3].to_vec());
        C_I_poly += &blinder_poly * z_I_poly;
        C_I_poly
    }

    fn get_u_poly(
        positions: &Vec<usize>,
        roots_N: &Vec<M::Scalar>,
        blinders: &Vec<M::Scalar>,
        z_v_m_poly: &UnivariatePolynomial<M::Scalar>,
    ) -> UnivariatePolynomial<M::Scalar> {
        let u_poly_coeffs: Vec<M::Scalar> = positions.iter().map(|&i_j| roots_N[i_j]).collect();
        let mut u_poly = UnivariatePolynomial::lagrange(u_poly_coeffs).ifft();
        let blinder_poly = UnivariatePolynomial::monomial(blinders[4..=6].to_vec());

        u_poly
    }

    pub fn verify(
        param: &CaulkParam<M>,
        m: usize,
        transcript: &mut (impl TranscriptRead<M::G1Affine, M::Scalar>
        + G2TranscriptRead<M::G2Affine, M::Scalar>),
    ) -> Result<(), Error> {
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
        UnivariateKzg::verify(
            &param.kzg_vp,
            &UnivariateKzgCommitment(g1_u),
            &alpha,
            &v1,
            transcript,
        )?;

        let v2 = transcript.read_field_element()?;
        UnivariateKzg::verify(
            &param.kzg_vp,
            &UnivariateKzgCommitment(g1_p1),
            &alpha,
            &v2,
            transcript,
        )?;

        UnivariateKzg::verify(
            &param.kzg_vp,
            &UnivariateKzgCommitment(g1_p2),
            &alpha,
            &<M as Engine>::Scalar::ZERO,
            transcript,
        )?;

        let g1_P1 = g1_Z_I + (g1_C_I * chi).into();

        let z_v_m_poly = Self::get_vanishing_poly(m);
        let mut g1_P2 = g1_H2 * -z_v_m_poly.evaluate(&alpha);
        g1_P2 += M::G1Affine::generator() * v2;
        g1_P2 -= cm * chi;

        let a = (g1_C - g1_C_I).into();
        let lhs = M::pairing(&a, &M::G2Affine::generator());
        let rhs = M::pairing(&g1_Z_I, &g2_H1);
        assert!(lhs == rhs);

        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::transcript::{InMemoryTranscript, Keccak256Transcript};
    use halo2_curves::bn256::Bn256;
    use halo2_curves::pairing::Engine;

    type Scalar = <Bn256 as Engine>::Scalar;
    type G1Affine = <Bn256 as Engine>::G1Affine;

    #[test]
    fn test_caulk() {
        let param = Caulk::<Bn256>::setup(10);
        let mut transcipt = Keccak256Transcript::new(());
        Caulk::<Bn256>::prove(
            &param,
            &scalars(&[1, 3, 2, 4]),
            &scalars(&[1, 2]),
            &mut transcipt,
        ).expect("Should not fail");
        let proof = transcipt.into_proof();
        let mut transcipt = Keccak256Transcript::from_proof((), &proof);
        Caulk::verify(&param, 2, &mut transcipt)
            .expect("Should not fail");
    }

    fn scalars(nums: &[u64]) -> Vec<Scalar> {
        nums.iter().map(|x| Scalar::from(*x)).collect::<Vec<_>>()
    }
}