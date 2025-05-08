use halo2_curves::bn256::Fr as Scalar;

pub fn log_ceil(n: usize) -> usize {
    // return the smallest int that >= log(n)
    (usize::BITS - (n - 1).leading_zeros()) as usize
}

pub fn hash_tuple(element: (Scalar, Scalar, Scalar), tau: Scalar, gamma: Scalar) -> Scalar {
    element.0 * gamma * gamma + element.1 * gamma + element.2 - tau
}
