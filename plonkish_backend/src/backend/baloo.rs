use rand::rngs::OsRng;
use std::{fmt::Debug, marker::PhantomData};

use halo2_curves::bn256::Bn256;

use crate::{
    poly::Polynomial,
    poly::univariate::UnivariatePolynomial,
    backend::baloo::preprocessor::preprocess,
    pcs::{
        PolynomialCommitmentScheme,
        univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam},
    },
    util::{
        arithmetic::{Field, PrimeField, MultiMillerLoop, root_of_unity},
        test::std_rng,
        Deserialize, DeserializeOwned, Itertools, Serialize,
        transcript::{InMemoryTranscript, TranscriptRead, TranscriptWrite, Keccak256Transcript},
    }
};


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
        Pcs: PolynomialCommitmentScheme<F, Polynomial = UnivariatePolynomial<F>>,
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
            let poly_size = m;
            print!("poly_size: {:?}\n", poly_size);
            let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
            Pcs::trim(&param, poly_size, 1).unwrap()
        };
        print!("lookup: {:?}\n", lookup);
        // Commit and open
        let proof = {
            let mut transcript = T::new(());
            let poly = Pcs::Polynomial::monomial(
                lookup.clone(),
            );
            print!("coeffs: {:?}\n", poly.coeffs());
            let comm = Pcs::commit_and_write(&pp, &poly, &mut transcript).unwrap();
            let point = Pcs::Polynomial::squeeze_point(m, &mut transcript);
            let eval = poly.evaluate(&point);
            transcript.write_field_element(&eval).unwrap();
            Pcs::open(&pp, &poly, &comm, &point, &eval, &mut transcript).unwrap();
            transcript.into_proof()
        };
        // Verify
        let result = {
            let mut transcript = T::from_proof((), proof.as_slice());
            Pcs::verify(
                &vp,
                &Pcs::read_commitment(&vp, &mut transcript).unwrap(),
                &<Pcs::Polynomial as Polynomial<F>>::squeeze_point(m, &mut transcript),
                &transcript.read_field_element().unwrap(),
                &mut transcript,
            )
        };
        assert_eq!(result, Ok(()));

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
    use std::collections::HashSet;
    use halo2_curves::bn256::Fr;
    use crate::util::transcript::{FieldTranscriptRead, FieldTranscriptWrite};
    use num_bigint::BigUint;

    type Pcs = UnivariateKzg<Bn256>;

    #[test]
    fn test_e2e() {
        let lookup = vec![Fr::one(), Fr::one()];
        let table = vec![Fr::one(), Fr::one()];
        let m = lookup.len();
        let mut rng = std_rng();

        // Setup
        let (pp, vp) = {
            let mut rng = OsRng;
            let poly_size = m;
            print!("poly_size: {:?}\n", poly_size);
            let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
            Pcs::trim(&param, poly_size, 1).unwrap()
        };
        print!("lookup: {:?}\n", lookup);

        // Commit and open
        let proof = {
            let mut transcript = Keccak256Transcript::new(());
            let poly = <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::monomial(lookup.clone());
            print!("coeffs: {:?}\n", poly.coeffs());
            let comm = Pcs::commit_and_write(&pp, &poly, &mut transcript).unwrap();
            let point = <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::squeeze_point(m, &mut transcript);
            let eval = poly.evaluate(&point);

            // Use the correct method to write the field element
            transcript.write_field_element(&eval).unwrap();

            Pcs::open(&pp, &poly, &comm, &point, &eval, &mut transcript).unwrap();
            transcript.into_proof()
        };

        // Verify
        let result = {
            let mut transcript = Keccak256Transcript::from_proof((), proof.as_slice());
            Pcs::verify(
                &vp,
                &Pcs::read_commitment(&vp, &mut transcript).unwrap(),
                &<Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial::squeeze_point(m, &mut transcript),

                // Use the correct method to read the field element
                &transcript.read_field_element().unwrap(),
                &mut transcript,
            )
        };
        assert_eq!(result, Ok(()));

    }

    #[test]
    fn test_round1() {
        let lookup = vec![Fr::one(), Fr::one()];
        let table = vec![Fr::one(), Fr::one()];
        let m = lookup.len();
        let mut rng = std_rng();

        // Setup
        let mut rng = OsRng;
        let poly_size = m;
        print!("poly_size: {:?}\n", poly_size);
        let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
        let (pp, vp) = Pcs::trim(&param, poly_size, 1).unwrap();
        print!("lookup: {:?}\n", lookup);

        // Commit and open
        let mut transcript = Keccak256Transcript::new(());
        // commit phi(X) on G1
        let phi_poly = <Pcs as PolynomialCommitmentScheme<Fr>>::Polynomial ::lagrange(lookup.clone());
        print!("coeffs: {:?}\n", phi_poly.coeffs());
        let phi_comm_1 = Pcs::commit_and_write(&pp, &phi_poly, &mut transcript).unwrap();
        // remove duplicated elements
        let t_values_from_lookup: HashSet<_> = lookup.into_iter().collect();
        // I: the index of t_values_from_lookup elements in sub table t_I
        let i_values: Vec<_> = t_values_from_lookup.iter().map(|elem| table.iter().position(|&x| x == *elem).unwrap()).collect();
        let root_of_unity = root_of_unity::<Fr>(m);
        // H_I = {ξ_i} , i = [1...k], ξ(Xi)
        let h_i: Vec<_> = i_values.iter().map(|&i| {
            let i_as_u64 = i as u64;
            root_of_unity.pow(&[i_as_u64][..])
        }).collect();
        // let h_i_interp_poly = InterpolationPoly::new(&h_i, &t_values_from_lookup); // TODO: implement interpolation polynomial
        let z_i_poly = UnivariatePolynomial::vanishing(&h_i, Fr::one());
        let z_i_comm_2 = Pcs::commit_monomial_g2(&param, &z_i_poly.coeffs());
        print!("z_i_comm_2: {:?}\n", z_i_comm_2);


    }
}