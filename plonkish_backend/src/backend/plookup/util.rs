use crate::{
    pcs::PolynomialCommitmentScheme,
    poly::univariate::UnivariatePolynomial,
    util::arithmetic::PrimeField,
};

pub(super) fn aggregate_poly<F: PrimeField>(
    scalar: &F,
    polys: Vec<&UnivariatePolynomial<F>>,
) -> UnivariatePolynomial<F> {
    assert!(polys.len()>0);

    let mut result = polys[0].clone();
    let mut power = scalar.clone();
    for poly in polys[1..].iter() {
        result += *poly * power;
        power *= scalar;
    }
    result
}

pub(super) fn aggregate_field<F: PrimeField>(
    scalar: &F,
    nums: Vec<&F>,
) -> F {
    assert!(nums.len()>0);

    let mut result = nums[0].clone();
    let mut power = scalar.clone();
    for num in nums[1..].iter() {
        result += **num * power;
        power *= scalar;
    }
    result
}

pub(super) fn aggregate_comm<
    F: PrimeField,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = UnivariatePolynomial<F>>
>(
    scalar: &F,
    comms: Vec<&Pcs::Commitment>,
) -> Pcs::Commitment {
    assert!(comms.len()>0);

    let mut result = comms[0].clone();
    let mut power = scalar.clone();
    for comm in comms[1..].iter() {
        // TODO: result = ec_lincomb([(result, 1), (comm, power)]);
        power *= scalar;
    }
    result
}