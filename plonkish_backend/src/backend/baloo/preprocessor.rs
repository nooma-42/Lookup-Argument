use crate::{
    backend::baloo::{
        key::{BalooProverKey, BalooVerifierKey},
        util::log_2,
    },
    pcs::{
        univariate::{UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam},
        PolynomialCommitmentScheme,
    },
    poly::{univariate::UnivariatePolynomial, Polynomial},
    util::arithmetic::{root_of_unity, Field},
    Error,
};
use halo2_curves::bn256::{Bn256, Fr};
use rand::rngs::OsRng;
use std::cmp::max;

type Pcs = UnivariateKzg<Bn256>;

/// Legacy generic preprocessing function for backward compatibility
pub fn preprocess(
    t: usize,
    m: usize,
) -> Result<
    (
        UnivariateKzgParam<Bn256>,
        UnivariateKzgProverParam<Bn256>,
        UnivariateKzgVerifierParam<Bn256>,
    ),
    Error,
> {
    let mut rng = OsRng;
    let poly_size = max(t.next_power_of_two() * 2, m.next_power_of_two() * 2);
    let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
    let (pp, vp) = Pcs::trim(&param, poly_size, 1).unwrap();

    Ok((param, pp, vp))
}

/// Table-specific preprocessing that performs all O(t log t) computations upfront.
/// This is the core performance optimization that moves expensive operations
/// from the online proving phase to the offline preprocessing phase.
pub fn preprocess_for_table(
    table: &[Fr],
    m_max: usize,
) -> Result<(BalooProverKey, BalooVerifierKey), Error> {
    let t = table.len();
    assert!(t.is_power_of_two(), "Table size must be a power of 2");
    
    let mut rng = OsRng;
    let poly_size = max(t.next_power_of_two() * 2, m_max.next_power_of_two() * 2);
    let param = Pcs::setup(poly_size, 1, &mut rng)?;
    let (pp, vp) = Pcs::trim(&param, poly_size, 1)?;
    
    let d = (1 << pp.k()) - 2;
    
    // Step 1: Interpolate the complete table to get T(X)
    let t_poly = UnivariatePolynomial::lagrange(table.to_vec()).ifft();
    
    // Step 2: Construct the vanishing polynomial Z_H(X) = X^t - 1
    let mut z_h_coeffs = vec![Fr::zero(); t + 1];
    z_h_coeffs[0] = -Fr::one();
    z_h_coeffs[t] = Fr::one();
    let z_h_poly = UnivariatePolynomial::monomial(z_h_coeffs);
    
    // Step 3: Commit to T(X) and Z_H(X)
    let table_comm = Pcs::commit_monomial(&pp, t_poly.coeffs());
    let z_h_comm = Pcs::commit_monomial(&pp, z_h_poly.coeffs());
    
    // Clone the commitments for separate use in prover and verifier keys
    let table_comm_clone = table_comm.clone();
    let z_h_comm_clone = z_h_comm.clone();
    
    // Step 4: Precompute all witness proofs for table elements
    // This is the most expensive part that we're moving to preprocessing
    let log_t = log_2(t);
    let t_root_of_unity = root_of_unity::<Fr>(log_t);
    let t_roots: Vec<Fr> = (0..t)
        .map(|i| t_root_of_unity.pow([i as u64]))
        .collect();
    
    let mut table_element_proofs = Vec::with_capacity(t);
    let mut subgroup_element_proofs = Vec::with_capacity(t);
    
    for i in 0..t {
        let omega_i = t_roots[i];
        let table_value = table[i];
        
        // Compute Q_i(X) = (T(X) - table[i]) / (X - ω^i)
        let divisor = UnivariatePolynomial::monomial(vec![-omega_i, Fr::one()]);
        let numerator = t_poly.clone() + (-table_value);
        let (q_i_poly, remainder) = numerator.div_rem(&divisor);
        
        // Sanity check: remainder should be zero
        assert!(remainder.coeffs().iter().all(|&c| c == Fr::zero()),
                "Polynomial division should have zero remainder");
        
        let q_i_comm = Pcs::commit_monomial(&pp, q_i_poly.coeffs());
        table_element_proofs.push(q_i_comm.to_affine());
        
        // Compute H_i(X) = Z_H(X) / (X - ω^i)
        let (h_i_poly, remainder) = z_h_poly.div_rem(&divisor);
        assert!(remainder.coeffs().iter().all(|&c| c == Fr::zero()),
                "Z_H polynomial division should have zero remainder");
        
        let h_i_comm = Pcs::commit_monomial(&pp, h_i_poly.coeffs());
        subgroup_element_proofs.push(h_i_comm.to_affine());
    }
    
    // Step 5: Construct the prover and verifier keys
    let prover_key = BalooProverKey {
        pp,
        param: param.clone(),
        table_comm: table_comm_clone,
        z_h_comm: z_h_comm_clone,
        table_element_proofs,
        subgroup_element_proofs,
        table: table.to_vec(),
        t,
        d,
    };
    
    let verifier_key = BalooVerifierKey {
        vp,
        param,
        table_comm,
        z_h_comm,
        t,
        d,
    };
    
    Ok((prover_key, verifier_key))
}

#[cfg(test)]
mod tests {
    use crate::backend::baloo::preprocessor::preprocess;

    #[test]
    fn test_preprocess() {
        let (param, pp, vp) = preprocess(10, 10).unwrap();
        println!("param: {:?}", param);
        println!("pp: {:?}", pp);
        println!("vp: {:?}", vp);
    }
}
