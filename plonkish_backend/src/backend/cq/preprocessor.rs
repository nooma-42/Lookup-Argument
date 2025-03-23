use crate::pcs::univariate::{
    UnivariateKzg, UnivariateKzgParam, UnivariateKzgProverParam, UnivariateKzgVerifierParam,
};
use crate::{
    backend::cq::util::{ec_fft, ec_ifft, fft, ifft, log_2},
    pcs::PolynomialCommitmentScheme,
    util::arithmetic::{root_of_unity, Field},
    Error,
};
use halo2_curves::bn256::{Bn256, Fq, Fr, G1Affine, G1};
use rand::rngs::OsRng;
use std::time::Instant;
use std::{cmp::max, iter};

type Pcs = UnivariateKzg<Bn256>;
type Scalar = Fr;

pub fn preprocess(
    t: usize,
    m: usize,
    table: &Vec<Fr>,
) -> Result<
    (
        UnivariateKzgParam<Bn256>,
        UnivariateKzgProverParam<Bn256>,
        UnivariateKzgVerifierParam<Bn256>,
        Vec<G1>,
    ),
    Error,
> {
    let mut rng = OsRng;
    let poly_size = max(t.next_power_of_two() * 2, m.next_power_of_two() * 2);

    // let start = Instant::now();
    let param = Pcs::setup(poly_size, 1, &mut rng).unwrap();
    // let duration1 = start.elapsed();
    // println!("\n ------------Setup: {}ms----------- \n",duration1.as_millis());

    let (pp, vp) = Pcs::trim(&param, poly_size, 1).unwrap();
    let mut powers_of_x = pp.monomial_g1()[..t].to_vec();

    // let start = Instant::now();
    let q_t_comm_poly_coeffs = precompute_with_fk(table, &mut powers_of_x);
    // let duration2 = start.elapsed();
    // println!("\n ------------precompute_with_fk: {}ms----------- \n",duration2.as_millis());

    Ok((param, pp, vp, q_t_comm_poly_coeffs))
}

fn fk(coeffs: &mut Vec<Fr>, powers_of_x: &mut Vec<G1Affine>) -> Vec<G1> {
    let start = Instant::now();
    //println!("\n ***************** Start fk() ****************");
    assert_eq!(coeffs.len(), powers_of_x.len(), "length should be equal");
    let n = coeffs.len();
    assert!(n.is_power_of_two());

    // Get first column of circulant matrix in length of 2 * len(coeffs)
    // For example: coeffs is [1, 2, 3, 4]
    // The first column of circulant matrix should be: [4, 0, 0, 0, 0, 0, 2, 3]
    let mut first_col = coeffs.to_vec();
    // [1, 2, 3, 4] -> [0, 2, 3, 4]
    first_col[0] = Scalar::from(0);

    // get first column of circulant matrix in 2n size
    // 1. padding 0
    // [0, 2, 3, 4] -> [0, 0, 0, 0, 0, 2, 3, 4]
    first_col.splice(0..0, iter::repeat(Scalar::from(0)).take(n));
    // 2. roll by 1 to right
    // [0, 0, 0, 0, 0, 2, 3, 4] -> [4, 0, 0, 0, 0, 0, 2, 3]
    first_col.rotate_right(1);

    // inverse srs: delete last one then inverse
    let mut inv_powers_of_x = powers_of_x[..n - 1].to_vec();
    inv_powers_of_x.reverse();
    let group1_zero = G1Affine {
        x: Fq::zero(),
        y: Fq::zero(),
    };
    inv_powers_of_x.push(group1_zero); // 假设 b.Z1 是 Scalar(0.0)

    // padding n 0s to the end
    let ec_neutral_vals = vec![group1_zero; n];
    inv_powers_of_x.extend(ec_neutral_vals);

    // We have circulant matrix C, C = F_inv * diag(F * first_col) * F
    // F: DFT matrix, F_inv: inverse DFT matrix
    // We want to get Q_T_comm_poly_coeffs = C * x = F_inv * diag(F * first_col) * F * x
    // 1. right hand side: F * x

    // let start1 = Instant::now();
    let rhs = ec_fft(&mut inv_powers_of_x);
    // let duration1 = start1.elapsed();
    // println!("\n ------------ec_fft: {}ms----------- \n",duration1.as_millis());

    // 2. middle hand side: F * first_col
    // let start2 = Instant::now();
    let mhs = fft(&mut first_col);
    // let duration2 = start2.elapsed();
    // println!("\n ------------fft: {}ms----------- \n",duration2.as_millis());

    // middle * right (element wise) to get diagonal: diag(F * first_col) * F * x
    let m_r_hs: Vec<G1Affine> = rhs
        .iter()
        .zip(mhs.iter())
        .map(|(&r, &m)| (r * m).into())
        .collect();

    // 3. ifft
    // let start3 = Instant::now();
    let result = ec_ifft(&mut m_r_hs.clone());
    // let duration3 = start3.elapsed();
    // println!("\n ------------ec_ifft: {}ms----------- \n",duration3.as_millis());

    // 4. return first n values
    let mut q_comm_poly_coeffs = result[..n].to_vec();
    //println!("\n ***************** End fk() ****************");

    // let duration = start.elapsed();
    // println!("\n ------------fk function: {}ms----------- \n",duration.as_millis());

    // let start4 = Instant::now();
    let t = coeffs.len();
    let log_values: usize = log_2(coeffs.len());
    let root_of_unity = root_of_unity::<Fr>(log_values);
    let roots = (0..coeffs.len())
        .map(|i| root_of_unity.pow([i as u64]))
        .collect::<Vec<Fr>>();
    // let duration4 = start4.elapsed();
    // println!("\n ------------root_of_unity: {}ms----------- \n",duration4.as_millis());

    let k_t_commit = ec_fft(&mut q_comm_poly_coeffs);

    let mut q_t_comm = vec![
        G1 {
            x: Fq::zero(),
            y: Fq::zero(),
            z: Fq::zero()
        };
        t
    ];
    // let start5 = Instant::now();
    for i in 0..t {
        let scale = roots[i] * (Scalar::from(t as u64).invert().unwrap());
        q_t_comm[i] = k_t_commit[i] * scale;
    }
    // let duration5 = start5.elapsed();
    // println!("\n ------------k_t_commit: {}ms----------- \n",duration5.as_millis());

    q_t_comm
}

pub fn precompute_with_fk(table: &Vec<Fr>, powers_of_x: &mut Vec<G1Affine>) -> Vec<G1> {
    // Convert table values to Scalars
    let t_values: Vec<Scalar> = table.iter().map(|&val| Scalar::from(val)).collect();

    // Compute t_poly_coeffs using inverse FFT
    let mut t_poly_coeffs = t_values.clone();
    // let start = Instant::now();
    t_poly_coeffs = ifft(&mut t_poly_coeffs);
    // let duration = start.elapsed();
    // println!("\n ------------ifft: {}ms----------- \n",duration.as_millis());

    // Compute h values with fk

    fk(&mut t_poly_coeffs, powers_of_x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::cq::preprocessor::preprocess;

    #[test]
    fn test_preprocess() {
        let table = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let (param, pp, vp, q_t_comm_poly_coeffs) = preprocess(4, 4, &table).unwrap();
        println!("param: {:?}", param);
        println!("pp: {:?}", pp);
        println!("vp: {:?}", vp);
        println!("q_t_comm_poly_coeffs: {:?}", q_t_comm_poly_coeffs);
    }
}
