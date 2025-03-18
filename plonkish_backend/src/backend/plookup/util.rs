use crate::{poly::univariate::UnivariatePolynomial, util::arithmetic::PrimeField};

pub(super) fn aggregate_poly<F: PrimeField>(
    scalar: &F,
    polys: Vec<&UnivariatePolynomial<F>>,
) -> UnivariatePolynomial<F> {
    assert!(!polys.is_empty());

    let mut result = polys[0].clone();
    let mut power = *scalar;
    for poly in polys[1..].iter() {
        result += *poly * power;
        power *= scalar;
    }
    result
}

pub(super) fn aggregate_field<F: PrimeField>(scalar: &F, nums: Vec<&F>) -> F {
    assert!(!nums.is_empty());

    let mut result = *nums[0];
    let mut power = *scalar;
    for num in nums[1..].iter() {
        result += **num * power;
        power *= scalar;
    }
    result
}
