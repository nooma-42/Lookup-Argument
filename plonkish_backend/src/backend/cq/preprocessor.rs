use rand::rngs::OsRng;
use std::cmp::max;
use crate::{
    pcs::PolynomialCommitmentScheme,
    Error,
};
use halo2_curves::bn256::Bn256;
use crate::pcs::univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam};

type Pcs = UnivariateKzg<Bn256>;

pub fn preprocess(
    t: usize,
    m: usize,
) -> Result<
    (
        UnivariateKzgParam<Bn256>,
        UnivariateKzgProverParam<Bn256>,
        UnivariateKzgVerifierParam<Bn256>,
    ), Error>
{
    let mut rng = OsRng;
    let poly_size = max(t.next_power_of_two() * 2, m.next_power_of_two() * 2);
    let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
    let (pp, vp) = Pcs::trim(&param, poly_size, 1).unwrap();

    Ok((param, pp, vp))
}

#[cfg(test)]
mod tests {
    use crate::backend::cq::preprocessor::preprocess;

    #[test]
    fn test_preprocess() {
        let (param, pp, vp) = preprocess(10, 10).unwrap();
        println!("param: {:?}", param);
        println!("pp: {:?}", pp);
        println!("vp: {:?}", vp);
    }
}
