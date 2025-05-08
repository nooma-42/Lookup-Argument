use crate::backend::lasso::affine_prod::*;
use crate::poly::multilinear::MultilinearPolynomial;
use crate::util::arithmetic::Field;
use halo2_curves::bn256::Fr as Scalar;

pub fn eq_mle_poly(w: &[Scalar]) -> AffinePolynomial {
    let mut prod = Vec::new();
    for i in 0..w.len() {
        let coeff = Scalar::from(2).mul(&w[i]).sub(&Scalar::one());
        let constant = Scalar::one() - w[i];

        let term = AffineTerm::new(coeff, i + 1, constant);
        prod.push(term);
    }

    let product = AffineProduct::new(Scalar::one(), &prod);
    AffinePolynomial::new(vec![product], Scalar::zero())
}

// Need test for this function
pub fn mle_mul<F: Field>(
    p1: &MultilinearPolynomial<F>,
    p2: &MultilinearPolynomial<F>,
) -> MultilinearPolynomial<F> {
    let mut result = Vec::with_capacity(p1.evals().len() * p2.evals().len());

    for &a in p1.evals() {
        for &b in p2.evals() {
            result.push(a * b);
        }
    }
    MultilinearPolynomial::new(result)
}
