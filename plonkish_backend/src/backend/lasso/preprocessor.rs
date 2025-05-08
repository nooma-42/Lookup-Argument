// plonkish_backend/src/backend/lasso/preprocessor.rs
use super::{LassoInfo, util::{get_poly_and_vec, get_index}}; 
use super::prover::LassoProverParam; 
use super::verifier::LassoVerifierParam; // Import from verifier module
use crate::{
    poly::multilinear::MultilinearPolynomial,
    util::arithmetic::{Field, PrimeField},
    Error,
    pcs::multilinear::MultilinearKzg,
    pcs::PolynomialCommitmentScheme,
};
use serde::{Serialize, de::DeserializeOwned};
use std::{fmt::Debug, hash::Hash, collections::HashMap};

// Helper to calculate necessary vector length (power of 2 >= m)
fn padded_size(m: usize) -> usize {
    if m == 0 { 0 } else { m.next_power_of_two() }
}

pub fn preprocess<F, Pcs>(
    pcs_param_raw: &Pcs::Param, // Use associated type
    info: &LassoInfo<F, Pcs>,
) -> Result<(LassoProverParam<F, Pcs>, LassoVerifierParam<F, Pcs>), Error>
where
    F: PrimeField + Hash + Eq + Send + Sync,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>, // Use base trait
    // Bounds match associated types from PolynomialCommitmentScheme
    Pcs::Commitment: Default + Clone + Debug + Serialize + DeserializeOwned + PartialEq + Send + Sync + AsRef<[Pcs::CommitmentChunk]>,
    Pcs::CommitmentChunk: Clone + Debug + Default + Serialize + DeserializeOwned + Send + Sync,
    Pcs::ProverParam: Clone + Debug + Serialize + DeserializeOwned + Send + Sync,
    Pcs::VerifierParam: Clone + Debug + Serialize + DeserializeOwned + Send + Sync,
    Pcs::Param: Clone + Debug + Serialize + DeserializeOwned + Send + Sync,
{
    // --- Parameters ---
    let l = info.l;
    let c = info.c;
    let k = info.k;
    let alpha = info.subtables.len();
    if alpha != c * k {
        return Err(Error::InvalidSnark("alpha != c * k".to_string()));
    }
    if l == 0 || c == 0 || k == 0 {
         return Err(Error::InvalidSnark("l, c, k must be > 0".to_string()));
    }
    let subtable_size = 1 << l;
    for table in &info.subtables {
        if table.len() != subtable_size {
            return Err(Error::InvalidSnark("Subtable size mismatch".to_string()));
        }
    }

    let witness = &info.witness;
    let m = witness.len();
    let logm = padded_size(m).ilog2() as usize; // Use padded size for MLEs
    let m_padded = 1 << logm;

    // --- Calculate derived values ---
    // 1. Pad witness `a` and get indices
    let mut a_padded = vec![F::ZERO; m_padded];
    let mut indices_decomposed: Vec<Vec<usize>> = Vec::with_capacity(m); // Store decomposed indices [idx_i^(1), ..., idx_i^(c)] for each witness element
    let mut witness_map = HashMap::new(); // Track witness elements for error checking if needed

    for (i, w) in witness.iter().enumerate() {
        a_padded[i] = *w;
        // Find index decomposition for witness element w
        // This requires the inverse of the g function. Assuming simple range table for now.
        let val_u64 = w.to_repr().as_ref()[0]; // Assumes Fr fits in u64 - needs care for larger fields
        let decomposed = get_index(val_u64 as usize, l, c)?; // Get [idx^(1), ..., idx^(c)]
        indices_decomposed.push(decomposed);
        witness_map.insert(w, i); // Track witness element
    }
    // Pad `a` by repeating the last element if m < m_padded
    if m > 0 && m < m_padded {
        let last_val = a_padded[m - 1];
        for i in m..m_padded {
            a_padded[i] = last_val;
        }
    }
    let (a_poly, _) = get_poly_and_vec(a_padded, logm)?;


    // 2. Calculate `dim` polynomials
    let mut dim_poly = Vec::with_capacity(c);
    let mut dim_values_all = Vec::with_capacity(c); // Store values for later use if needed
    for j in 0..c { // For each chunk j
        let mut dim_j_values = vec![0; m_padded];
        for i in 0..m {
            dim_j_values[i] = indices_decomposed[i][j]; // dim_j[i] = idx_i^(j)
        }
        // Pad dim_j if needed
        if m > 0 && m < m_padded {
            let last_idx = dim_j_values[m - 1];
             for i in m..m_padded {
                 dim_j_values[i] = last_idx;
             }
        }
        let (poly, values_f) = get_poly_and_vec(dim_j_values.into_iter().map(|v| F::from(v as u64)).collect(), logm)?;
        dim_poly.push(poly);
        dim_values_all.push(values_f); // Store F representation
    }


    // 3. Calculate `E`, `read`, `write`, `final` polynomials for each subtable
    let mut E_poly = Vec::with_capacity(alpha);
    let mut read_poly = Vec::with_capacity(alpha);
    let mut write_poly = Vec::with_capacity(alpha); // Store temporarily
    let mut final_poly = Vec::with_capacity(alpha);

    let mut E_values_all = Vec::with_capacity(alpha);
    let mut read_values_all = Vec::with_capacity(alpha);
    let mut write_values_all = Vec::with_capacity(alpha);
    let mut final_values_all = Vec::with_capacity(alpha);


    for i in 0..alpha { // For each subtable T_i
        let chunk_index = i / k; // Which dim polynomial (0..c-1) to use for indices
        let subtable = &info.subtables[i];

        let mut E_i_values = vec![F::ZERO; m_padded];
        let mut read_i_values = vec![F::ZERO; m_padded];
        let mut write_i_values = vec![F::ZERO; m_padded];
        let mut final_i_counters = vec![0u64; subtable_size]; // Track counters during simulation

        for access_idx in 0..m {
            let table_idx = dim_values_all[chunk_index][access_idx].to_repr().as_ref()[0] as usize; // idx_i^(chunk_index)
            if table_idx >= subtable_size {
                 return Err(Error::InvalidSnark(format!("Subtable index {} out of bounds {}", table_idx, subtable_size)));
            }
            E_i_values[access_idx] = subtable[table_idx];
            let current_counter = final_i_counters[table_idx];
            read_i_values[access_idx] = F::from(current_counter);
            write_i_values[access_idx] = F::from(current_counter + 1);
            final_i_counters[table_idx] = current_counter + 1;
    }

        // Pad E, read, write if needed
        if m > 0 && m < m_padded {
            let last_E = E_i_values[m - 1];
            let last_read = read_i_values[m - 1];
            let last_write = write_i_values[m - 1];
            for idx in m..m_padded {
                E_i_values[idx] = last_E;
                read_i_values[idx] = last_read;
                write_i_values[idx] = last_write;
            }
        }

        // Convert final counters to field elements
        let final_i_values_f: Vec<F> = final_i_counters.iter().map(|&cnt| F::from(cnt)).collect();

        // Create polynomials
        let (poly_E, vals_E) = get_poly_and_vec(E_i_values, logm)?;
        let (poly_read, vals_read) = get_poly_and_vec(read_i_values, logm)?;
        let (poly_write, vals_write) = get_poly_and_vec(write_i_values, logm)?;
        let (poly_final, vals_final) = get_poly_and_vec(final_i_values_f, l)?; // final has size 2^l (subtable size)

        E_poly.push(poly_E);
        read_poly.push(poly_read);
        write_poly.push(poly_write); // Store temporarily for prover
        final_poly.push(poly_final);

        E_values_all.push(vals_E);
        read_values_all.push(vals_read);
        write_values_all.push(vals_write);
        final_values_all.push(vals_final);
    }


    // 4. Trim PCS parameters
    // Estimate max poly size and batch size
    let max_vars_for_commit = logm.max(l).max(l + 1); // Consider GP polys (l+1 or logm+1)
    let pcs_poly_size = 1 << max_vars_for_commit;
    // Estimate batch size based on number of polys committed/opened - Not used by MultilinearKzg::trim
    // let pcs_batch_size = alpha * 4 + c + 1; // a + c*dim + alpha*E + alpha*read + alpha*final + (alpha*4 for grand products)
    // let (pcs_pp, pcs_vp) = Pcs::trim(pcs_param_raw, pcs_poly_size, pcs_batch_size)?;
    let (pcs_pp, pcs_vp) = Pcs::trim(pcs_param_raw, pcs_poly_size, 0)?; // MultilinearKzg ignores batch_size


    // 5. Commit to polynomials (parallelize?)
    let a_comm = Pcs::commit(&pcs_pp, &a_poly)?;
    let dim_comm_results: Result<Vec<_>, Error> = dim_poly.iter().map(|p| Pcs::commit(&pcs_pp, p)).collect();
    let dim_comm = dim_comm_results?;
    let E_comm_results: Result<Vec<_>, Error> = E_poly.iter().map(|p| Pcs::commit(&pcs_pp, p)).collect();
    let E_comm = E_comm_results?;
    let read_comm_results: Result<Vec<_>, Error> = read_poly.iter().map(|p| Pcs::commit(&pcs_pp, p)).collect();
    let read_comm = read_comm_results?;
    let final_comm_results: Result<Vec<_>, Error> = final_poly.iter().map(|p| Pcs::commit(&pcs_pp, p)).collect();
    let final_comm = final_comm_results?;


    // 6. Construct Params
    let pp = LassoProverParam {
        l, c, k, alpha, logm,
        subtables: info.subtables.clone(),
        a_poly, dim_poly, E_poly, read_poly, write_poly, final_poly,
        a_comm: a_comm.clone(),
        dim_comm: dim_comm.clone(),
        E_comm: E_comm.clone(),
        read_comm: read_comm.clone(),
        final_comm: final_comm.clone(),
        pcs_param: pcs_pp,
        // Initialize temporary storage fields
        S0_poly: None,
        S_poly: None,
        RS_poly: None,
        WS_poly: None,
    };

    // Verifier doesn't get the polynomials directly, only commitments and params
    let vp = LassoVerifierParam {
        l, c, k, alpha, logm,
        subtables: info.subtables.clone(), // Verifier needs tables for hash checks
        a_comm, dim_comm, E_comm, read_comm, final_comm,
        pcs_param: pcs_vp,
        _marker: PhantomData,
    };

    Ok((pp, vp))
}

