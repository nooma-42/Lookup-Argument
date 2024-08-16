use rand::RngCore;
use std::{collections::HashSet, fmt::Debug, hash::Hash, iter, marker::PhantomData};

use bincode::Error;
use halo2_proofs::transcript;

use crate::{
    pcs::PolynomialCommitmentScheme,
    backend::baloo::{
        preprocessor::preprocess,
        // prover::prove,
        // verifier::verify,
    },
    util::{
        arithmetic::PrimeField,
        test::std_rng,
        Deserialize, DeserializeOwned, Itertools, Serialize,
    },
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

#[derive(Clone, Debug)]
pub struct Baloo<Pcs>
{
    setup: Setup,
    m: usize,
    phi_poly: Polynomial,
    phi_comm_1: G1Commitment,
    h_i: Vec<RootOfUnity>,
    z_i_poly: Polynomial,
    col: Vec<usize>,
    v_poly: Polynomial,
    t_i_poly: Polynomial,
    roots_of_unity_n: Vec<RootOfUnity>,
    table: Vec<Scalar>,
}


impl<F, Pcs> PlonkishBackend<F> for Baloo<Pcs>
where
    F: PrimeField + Hash + Serialize + DeserializeOwned,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = UnivariatePolynomial<F>>,
{
    type Pcs = Pcs;
    type ProverParam = BalooProverParam;
    type VerifierParam = BalooVerifierParam<F, Pcs>;

    fn setup(
        n: usize,
        rng: impl RngCore
    ) -> Result<Pcs::Param, Error> {
        Pcs::setup(n, 1, &mut rng).unwrap();
    }

    fn preprocess(
        param: &Pcs::Param,
        n: usize
    ) -> Result<(Self::ProverParam, Self::VerifierParam), Error> {
        preprocess(param, n);
    }

    fn prove(
        pp: &Self::ProverParam,
        table: Vec<F>,
        lookup: Vec<F>,
        transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>
    ) -> Result<(), Error> {
        // t_poly in in lagrange basis
        let msg1 = round1(pp, lookup, table);
        let (alpha, beta) = transcript.round_1(msg_1);

        let msg2 = round2(pp);
        let (gamma, zeta) = transcript.round_2(msg_2);

        let msg3 = round3(pp);

        Ok(())
    }

    fn verify(
            vp: &Self::VerifierParam,
            instances: &[Vec<F>],
            transcript: &mut impl crate::util::transcript::TranscriptRead<crate::pcs::CommitmentChunk<F, Self::Pcs>, F>,
            rng: impl RngCore,
        ) -> Result<(), crate::Error> {
            Ok(())


    }
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

fn round1(pp: &Self::ProverParam, lookup: Vec<F>, table: Vec<F>) -> Vec<F> {
    let m = len(lookup);
    // let phi_values = lookup.iter().map(|&x| PrimeField::<F>::from(x)).collect::<Vec<F>>();
    let phi_poly = UnivariatePolynomial::new(lookup);
    // commit phi(X) on G1
    let phi_comm_1 = Pcs::commit(&pp, &phi_poly);
    // remove duplicated elements
    let t_values_from_lookup: HashSet<_> = lookup.into_iter().collect();
    // I: the index of t_values_from_lookup elements in sub table t_I
    let i_values: Vec<_> = t_values_from_lookup.iter().map(|elem| table.iter().position(|&x| x == *elem).unwrap()).collect();
    let roots_of_unity = Vec::new(); // TODO: get roots of unity
    // H_I = {ξ_i} , i = [1, k], ξ(Xi)
    let h_i = i_values.iter().map(|&i| roots_of_unity[i]).collect();
    let h_i_interp_poly = InterpolationPoly::new(&h_i, &t_values_from_lookup); // TODO: implement interpolation polynomial
    // Multiple all root polynomial: (X - I_0)(X - I_1)(X - I_2)...
    let mut col_values = Vec::new();
    let mut v_values = Vec::new();
    for i in 0..m {
        // find the index of 1 in jth row of M
        let col_i = t_values_from_lookup.iter().position(|&x| x == lookup[i]).unwrap();
        col_values.push(col_i);
        let col_i_root = h_i[col_i];
        // Note: v = 1 / col_i_root in paper
        // Here we use different construction that does not affect the verification
        let v = col_i_root;
        v_values.push(v);
    }
    // Vanishing polynomial in coefficient form
    self.z_i_poly = h_i_interp_poly.vanishing_poly();
    assert_eq!(col_values.iter().map(|&elem| t_values_from_lookup[elem]).collect::<Vec<_>>(), lookup);
    // ξ(x) polynomial
    let v_poly = Polynomial::new(v_values, Basis::Lagrange);
    // Refer to section 5. cannot use FFT due to it's not in multiplicative subgroup
    let t_interp_poly = InterpolationPoly::new(&h_i, &t_values_from_lookup);
    self.t_i_poly = t_interp_poly.poly();

    self.col = col_values.clone();
    self.h_i = h_i;
    self.phi_poly = phi_poly.ifft();
    self.phi_comm_1 = setup.commit_g1(&self.phi_poly);
    self.v_poly = v_poly.ifft();

    // Commit
    // π1 = ([z_i]_2 = [z_i(x)]_2, [v]_1 = [v(x)]_1, t = [t(x)]_1)
    self.z_i_comm_2 = setup.commit_g2(&self.z_i_poly);
    self.v_comm_1 = setup.commit_g1(&self.v_poly);
    self.t_i_comm_1 = setup.commit_g1(&self.t_i_poly);

    Message1::new(self.z_i_comm_2.clone(), self.v_comm_1.clone(), self.t_i_comm_1.clone(), self.phi_comm_1.clone())

}

// Placeholder types and methods for the sake of example
struct Setup;
struct Scalar(usize);
struct Polynomial;
struct G1Commitment;
struct RootOfUnity;
struct Message1;

impl Polynomial {
    fn new(values: Vec<Scalar>, basis: Basis) -> Self {
        // Implementation here
        Polynomial
    }

    fn ifft(&self) -> Self {
        // Implementation here
        Polynomial
    }
}


impl InterpolationPoly {
    fn new(h_i: &[RootOfUnity], t_values: &[Scalar]) -> Self {
        // Implementation here
        InterpolationPoly
    }

    fn vanishing_poly(&self) -> Polynomial {
        // Implementation here
        Polynomial
    }

    fn poly(&self) -> Polynomial {
        // Implementation here
        Polynomial
    }
}
