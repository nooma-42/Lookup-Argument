use crate::{
    backend::logupgkr::{LogupGkrInfo, LogupGkrProverParam, LogupGkrVerifierParam},
    poly::multilinear::MultilinearPolynomial,
    util::arithmetic::PrimeField,
    Error,
};
use crate::backend::logupgkr::util::{
    generate_binary_combinations, create_multilinear_poly, p, q
};

pub fn preprocess<F>(
    info: &LogupGkrInfo<F>,
) -> Result<(LogupGkrProverParam<F>, LogupGkrVerifierParam<F>), Error>
where
    F: PrimeField,
{
    // Create MultilinearPolynomials from the provided value maps
    let m_poly = create_multilinear_poly(info.m_values.clone());
    let t_poly = create_multilinear_poly(info.t_values.clone());
    
    let mut w_polys = Vec::new();
    for w_map in &info.w_values {
        w_polys.push(create_multilinear_poly(w_map.clone()));
    }
    
    // Determine the number of variables from the first map (assuming all maps have the same dimensions)
    let num_vars = info.m_values.keys().next().unwrap().len();
    
    // Generate all binary inputs for the multilinear polynomials
    let all_inputs = generate_binary_combinations(2 * num_vars as u32);
    
    // Calculate p and q values for all inputs
    let mut p_values = Vec::new();
    let mut q_values = Vec::new();
    
    for input in all_inputs {
        let x = &input[..num_vars];
        let y = &input[num_vars..];
        
        p_values.push(p(x, y, &m_poly));
        q_values.push(q(x, y, &t_poly, &w_polys, info.a));
    }
    
    // Create multilinear polynomials from p and q values
    let ps = MultilinearPolynomial::new(p_values);
    let qs = MultilinearPolynomial::new(q_values);
    
    // Create claims (using None as placeholders)
    let claims = vec![None; 2];
    let (p_0s, q_0s) = claims.split_at(1);
    
    // Create prover and verifier parameters
    let pp = LogupGkrProverParam {
        m_poly,
        t_poly,
        w_polys,
        a: info.a,
        ps,
        qs,
        p_0s: p_0s.to_vec(),
        q_0s: q_0s.to_vec(),
    };
    
    let vp = LogupGkrVerifierParam {
        num_vars,
        p_0s: p_0s.to_vec(),
        q_0s: q_0s.to_vec(),
    };
    
    Ok((pp, vp))
} 