use rand::{rngs::OsRng, Rng};
use rand::RngCore;
use std::{collections::HashSet, fmt::Debug, hash::Hash, iter, marker::PhantomData};

use bincode::Error;
use halo2_proofs::transcript;
use halo2_curves::bn256::Bn256;

use crate::pcs::univariate;
use crate::{
    poly::Polynomial,
    poly::univariate::{UnivariateBasis::*, UnivariatePolynomial},
    backend::baloo::preprocessor::preprocess,
    pcs::{
        PolynomialCommitmentScheme,
        univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam},
    },
    util::{
        arithmetic::{Field, PrimeField, MultiMillerLoop},
        test::std_rng,
        Deserialize, DeserializeOwned, Itertools, Serialize,
        transcript::{InMemoryTranscript, TranscriptRead, TranscriptWrite, Keccak256Transcript},
    }
};

use super::PlonkishBackend;

pub mod preprocessor;
pub mod prover;
pub mod verifier;

#[derive(Clone, Debug)]
pub struct BalooProverParam //<F, Pcs>
// where
//     F: PrimeField,
//     Pcs: PolynomialCommitmentScheme<F>,
{
    pub(crate) num_vars: usize,
}

#[derive(Clone, Debug)]
pub struct BalooVerifierParam<F, Pcs>
where
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F>,
{
    // [z_H_comm_1, t_comm_1]
    pub(crate) preprocess_comms: Vec<Pcs::Commitment>,
}
pub struct Baloo<F> {
    table: Vec<F>,
    // round1: (&self, &UnivariateKzgProverParam<M>, Vec<F>, Vec<F>) -> Vec<F>
}

impl<F: Field> Baloo<F>
{
    pub fn new(table: Vec<F>) -> Self {
        Baloo { table }
    }
    /*
    How to calculate ξ(x)(or v(x) in code)

    H = [1, ω, ω^2, ω^3, ω^4, ω^5, ω^6, ω^7]
    table = [1, 2, 3, 4, 5, 6, 7, 8]
    lookup = [3, 7, 3, 4]
    t_I:  [3, 4, 7] # choose non-duplicated elements from lookup
    I:  [2, 3, 6] # get indexes from table
    s ∈ I = [2, 3, 6]
    H_I = {ω^s} = [ω^2, ω^3, ω^6]
    k = len(I) = 3
    vanishing polynomial z_I(X) = (X - H_I_0)(X - H_I_1)(X - H_I_2)
                                = (X - ω^2)(X - ω^3)(X - ω^6)
    M * t = lookup
    M = [
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
    ]
    m = len(lookup) # 4
    col[0] = M[0].index(1)
    col[1] = M[1].index(1)
    col[2] = M[2].index(1)
    col[3] = M[3].index(1)
    col: [2, 6, 2, 3]
    xi = H_I[col_i] = [ω^2, ω^6, ω^2, ω^3]
    Interpolation with xi and get polynomial: ξ(x)
    */
    fn round1<Pcs, T>(&self, lookup: Vec<F>) -> Vec<F>
    where
        Pcs: PolynomialCommitmentScheme<F>,
        T: TranscriptRead<Pcs::CommitmentChunk, F>
            + TranscriptWrite<Pcs::CommitmentChunk, F>
            + InMemoryTranscript<Param = ()>,
    {
        // type ProverParam = UnivariateKzgProverParam<M>;
        let m = lookup.len();
        let mut rng = std_rng();
        // let phi_values = lookup.iter().map(|&x| PrimeField::<F>::from(x)).collect::<Vec<F>>();
        // Setup
        let (pp, vp) = {
            let mut rng = OsRng;
            let poly_size = 1 << m;
            let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
            Pcs::trim(&param, poly_size, 1).unwrap()
        };
        // Commit and open
        let proof = {
            let mut transcript = T::new(());
            let poly = <UnivariatePolynomial<F>>::lagrange(
                lookup.clone(),
            );
        };
        // let phi_poly = Pcs::new(lookup);
        // commit phi(X) on G1
        // let phi_comm_1 = Self::setup.commit_g1(&pp, &phi_poly);
        // let phi_comm_1 = Pcs::commit(&pp, &phi_poly);
        // remove duplicated elements
        // let t_values_from_lookup: HashSet<_> = lookup.into_iter().collect();
        // // I: the index of t_values_from_lookup elements in sub table t_I
        // let i_values: Vec<_> = t_values_from_lookup.iter().map(|elem| table.iter().position(|&x| x == *elem).unwrap()).collect();
        // let roots_of_unity = Vec::new(); // TODO: get roots of unity
        // // H_I = {ξ_i} , i = [1, k], ξ(Xi)
        // let h_i = i_values.iter().map(|&i| roots_of_unity[i]).collect();
        // let h_i_interp_poly = InterpolationPoly::new(&h_i, &t_values_from_lookup); // TODO: implement interpolation polynomial
        // // Multiple all root polynomial: (X - I_0)(X - I_1)(X - I_2)...
        // let mut col_values = Vec::new();
        // let mut v_values = Vec::new();
        // for i in 0..m {
        //     // find the index of 1 in jth row of M
        //     let col_i = t_values_from_lookup.iter().position(|&x| x == lookup[i]).unwrap();
        //     col_values.push(col_i);
        //     let col_i_root = h_i[col_i];
        //     // Note: v = 1 / col_i_root in paper
        //     // Here we use different construction that does not affect the verification
        //     let v = col_i_root;
        //     v_values.push(v);
        // }
        // // Vanishing polynomial in coefficient form
        // Self::z_i_poly = h_i_interp_poly.vanishing_poly();
        // assert_eq!(col_values.iter().map(|&elem| t_values_from_lookup[elem]).collect::<Vec<_>>(), lookup);
        // // ξ(x) polynomial
        // let v_poly = Polynomial::new(v_values, Basis::Lagrange);
        // // Refer to section 5. cannot use FFT due to it's not in multiplicative subgroup
        // let t_interp_poly = InterpolationPoly::new(&h_i, &t_values_from_lookup);
        // Self::t_i_poly = t_interp_poly.poly();

        // Self::col = col_values.clone();
        // Self::h_i = h_i;
        // Self::phi_poly = phi_poly.ifft();
        // Self::phi_comm_1 = Self::setup.commit_g1(&Self::phi_poly);
        // Self::v_poly = v_poly.ifft();

        // // Commit
        // // π1 = ([z_i]_2 = [z_i(x)]_2, [v]_1 = [v(x)]_1, t = [t(x)]_1)
        // Self::z_i_comm_2 = Self::setup.commit_g2(&Self::z_i_poly);
        // Self::v_comm_1 = Self::setup.commit_g1(&Self::v_poly);
        // Self::t_i_comm_1 = Self::setup.commit_g1(&Self::t_i_poly);

        // Message1::new(Self::z_i_comm_2.clone(), Self::v_comm_1.clone(), Self::t_i_comm_1.clone(), Self::phi_comm_1.clone())
        Vec::new()
    }

    // type Pcs = Pcs;
    // type ProverParam = BalooProverParam;
    // type VerifierParam = BalooVerifierParam<F, Pcs>;

}

#[cfg(test)]
mod tests {
    use super::*;
    type Pcs = UnivariateKzg<Bn256>;

    #[test]
    fn test_prover() {
        let lookup = vec![Field::ONE];
        let table = vec![Field::ONE];
        let m = lookup.len();
        let mut rng = std_rng();
        // let phi_values = lookup.iter().map(|&x| PrimeField::<F>::from(x)).collect::<Vec<F>>();
        let param = Pcs::setup(m, 1, &mut rng).unwrap();
        let baloo = Baloo{ table };

        let pp = param;
        let msg_1 = baloo.round1::<Pcs, Keccak256Transcript<_>>(lookup);
    }
}