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
    let num_vars_column = info.w_values.len().ilog2() as usize;
    let num_vars_row = info.w_values[0].keys().next().unwrap().len(); // NOTE: [0] is arbitrary, all w_values have the same number of variables
    let num_vars = num_vars_column + num_vars_row;
    
    // Generate all binary inputs for the multilinear polynomials
    let all_inputs = generate_binary_combinations((num_vars_column + num_vars_row) as u32);
    
    // Calculate p and q values for all inputs
    let mut p_values = Vec::new();
    let mut q_values = Vec::new();
    
    for input in all_inputs {
        // row part of multilinear input
        let x = &input[..num_vars_row];
        // column part of multilinear input
        let y = &input[num_vars_row..]; 
        
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