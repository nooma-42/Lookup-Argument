
use crate::{
    backend::baloo::{BalooProverParam, BalooVerifierParam},
    pcs::PolynomialCommitmentScheme,
    poly::Polynomial,
    util::{
        arithmetic::PrimeField,
        test::std_rng,
    },
    poly::univariate::{UnivariateBasis, UnivariatePolynomial},
    Error,
};

pub fn preprocess<F: PrimeField, Pcs: PolynomialCommitmentScheme<F>>(
    param: &Pcs::Param,
    n: usize,
) -> Result<
    (
        Pcs::ProverParam,
        BalooVerifierParam<F, Pcs>,
    ),
    Error,
    > {
        // let n = 1 << circuit_info.k;
        let mut rng = std_rng();
        // let param = Pcs::setup(n, 1, &mut rng).unwrap();
        let (pp, _) = Pcs::trim(&param, n, 1).unwrap();

        // todo
        let z_h_poly = Pcs::Polynomial::rand(n, &mut rng);
        // todo
        let t_poly = Pcs::Polynomial::rand(n, &mut rng);
        let z_h_comm_1 = Pcs::commit(&pp, &z_h_poly).unwrap();
        let t_comm_1 = Pcs::commit(&pp, &t_poly).unwrap();
        // Alternatively, if you know the number of elements you will store, you can use `with_capacity`
        let mut commitments: Vec<Pcs::Commitment> = Vec::with_capacity(2); // Adjust the capacity as needed
        commitments.push(z_h_comm_1);
        commitments.push(t_comm_1);
        let vp = BalooVerifierParam {
            preprocess_comms: commitments,
        };
    Ok((pp, vp))
}

#[cfg(test)]
mod tests {
    use crate::{
        backend::baloo::{preprocessor::preprocess, BalooVerifierParam},
        pcs::PolynomialCommitmentScheme,
        util::{
            arithmetic::PrimeField,
            test::std_rng,
        },
    };
    use crate::pcs::univariate::UnivariateKzg;
    use halo2_curves::bn256::{Bn256, Fr}; // Use Fr for the field type

    type Pcs = UnivariateKzg<Bn256>;

    #[test]
    fn test_preprocess() {
        let mut rng = std_rng();
        let k = 10;
        let n = 1 << k;
        let param = Pcs::setup(n, 1, &mut rng).unwrap();
        let (pp, vp) = preprocess::<Fr, Pcs>(&param, n).unwrap(); // Explicitly specify Fr for F

        assert!(true);
    }
}
