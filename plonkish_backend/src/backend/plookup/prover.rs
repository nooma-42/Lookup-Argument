use std::ops::Div;

use crate::{
    backend::plookup::{PlookupInfo, PlookupProverParam},
    pcs::PolynomialCommitmentScheme,
    poly::univariate::UnivariatePolynomial,
    util::{arithmetic::PrimeField, transcript::TranscriptWrite},
    halo2_curves::ff::WithSmallOrderMulGroup,
    Error,
};

pub(super) fn prove<
    F: PrimeField + WithSmallOrderMulGroup<3>, 
    Pcs: PolynomialCommitmentScheme<F, Polynomial = UnivariatePolynomial<F>>
>(
    pp: PlookupProverParam<F, Pcs>,
    info: &PlookupInfo<F>,
    transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
) -> Result<(), Error> {
    let order = 1 << info.k;
    // TODO: add another entry to Error for lookup itself?
    if info.table.len() != order {
        return Err(Error::InvalidPcsParam(String::from("unmatched table length and order")))
    }
    if info.lookup.len() >= order {
        return Err(Error::InvalidPcsParam(String::from("lookup length too long")))
    }
    let t = info.table.clone();
    let mut f = info.lookup.clone();

    // round 1
    // pad f to length order-1
    let last = f[f.len()-1];
    while f.len() < order-1 {
        f.push(last);
    }
    let s = sorted_by_table(&t, &f);
    let h1 = s[..order].to_vec();
    let h2 = s[order-1..].to_vec();
    let t_poly = UnivariatePolynomial::lagrange(t.clone()).ifft();
    let mut f_clone = f.clone();
    f_clone.push(last); // pad to length order for polynomial
    let f_poly = UnivariatePolynomial::lagrange(f_clone).ifft();
    let h1_poly = UnivariatePolynomial::lagrange(h1.clone()).ifft();
    let h2_poly = UnivariatePolynomial::lagrange(h2.clone()).ifft();
    // write f_comm, h1_comm, h2_comm
    Pcs::commit_and_write(&pp.pcs, &f_poly, transcript)?;
    Pcs::commit_and_write(&pp.pcs, &h1_poly, transcript)?;
    Pcs::commit_and_write(&pp.pcs, &h2_poly, transcript)?;

    // round 2
    let mut challenges: Vec<F> = Vec::with_capacity(2);
    challenges.extend(transcript.squeeze_challenges(2));
    let beta = &challenges[0];
    let gamma = &challenges[1];
    // TODO: return a vector
    let z = compute_inner_table(beta, gamma, &f, &t, &s);
    // write z_comm
    let z_poly = UnivariatePolynomial::lagrange(z.clone()).ifft();
    Pcs::commit_and_write(&pp.pcs, &z_poly, transcript)?;

    // round 3
    let delta = &transcript.squeeze_challenge();
    let q_poly = compute_quotient_polynomial(
        &pp.g, beta, gamma, delta, &t, &h1, &h2, &z,
        t_poly.clone(), f_poly.clone(), h1_poly.clone(), h2_poly.clone(), z_poly.clone(),
    );
    // write q_comm
    Pcs::commit_and_write(&pp.pcs, &q_poly, transcript).unwrap();

    // round 4
    let zeta = transcript.squeeze_challenge();
    let f_eval = f_poly.evaluate(&zeta);
    let h1_eval = h1_poly.evaluate(&zeta);
    let h2_eval = h2_poly.evaluate(&zeta);
    let z_eval = z_poly.evaluate(&zeta);
    let q_eval = q_poly.evaluate(&zeta);
    transcript.write_field_elements(
        vec![&f_eval, &h1_eval, &h2_eval, &z_eval, &q_eval])?;
    let g_zeta = pp.g * zeta;
    let h1_g_eval = h1_poly.evaluate(&g_zeta);
    let h2_g_eval = h2_poly.evaluate(&g_zeta);
    let z_g_eval = z_poly.evaluate(&g_zeta);
    transcript.write_field_elements(
        vec![&h1_g_eval, &h2_g_eval, &z_g_eval])?;

    // round 5
    let eps = &transcript.squeeze_challenge();
    let _agg_witness_comm = {
        let agg_poly = aggregate_poly(eps, vec![
            &f_poly, &h1_poly, &h2_poly, &z_poly, &q_poly
        ]);
        let agg_eval = aggregate_field(eps, vec![
            &f_eval, &h1_eval, &h2_eval, &z_eval, &q_eval
        ]);
        let numer = agg_poly + -agg_eval;
        let denom = UnivariatePolynomial::monomial(vec![-zeta, F::ONE]);
        let agg_witness = numer.div(&denom);
        Pcs::commit_and_write(&pp.pcs, &agg_witness, transcript).unwrap();
    };
    let _agg_g_witness_comm = {
        let agg_poly = aggregate_poly(eps, vec![
            &h1_poly, &h2_poly, &z_poly
        ]);
        let agg_eval = aggregate_field(eps, vec![
            &h1_g_eval, &h2_g_eval, &z_g_eval
        ]);
        let numer = agg_poly + -agg_eval;
        let denom = UnivariatePolynomial::monomial(vec![-g_zeta, F::ONE]);
        let agg_witness = numer.div(&denom);
        Pcs::commit_and_write(&pp.pcs, &agg_witness, transcript)?;
    };
    Ok(())
}

fn sorted_by_table<F: PrimeField>(
    table: &Vec<F>,
    witness: &Vec<F>,
) -> Vec<F> {
    let n = table.len();
    let count = {
        let mut ret = vec![1; table.len()];
        for w in witness {
            let mut i = 0;
            while i < n {
                if *w == table[i] {
                    ret[i] += 1;
                    break;
                }
                i += 1;
            }
            assert!(i < n); // ensure w in table
        }
        ret
    };
    let sorted = {
        let mut ret: Vec<F> = Vec::with_capacity(n+witness.len());
        let mut i = 0;
        while i < n {
            let mut cnt = count[i];
            while cnt > 0 {
                ret.push(table[i]);
                cnt -= 1;
            }
            i += 1;
        }
        ret
    }; 
    sorted
}

fn compute_inner_table<F: PrimeField>(
    beta: &F,
    gamma: &F,
    f: &Vec<F>,
    t: &Vec<F>,
    s: &Vec<F>,
) -> Vec<F> {
    let n = t.len();
    let mut numer = F::ONE;
    let mut denom = F::ONE;
    let values = {
        let mut ret: Vec<F> = Vec::with_capacity(n);
        ret.push(F::ONE);
        let mut i = 1;
        while i < n {
            numer = numer * (F::ONE+*beta) * (*gamma+f[i-1])
                * (*gamma*(F::ONE+*beta)+t[i-1]+*beta*t[i]);
            denom = denom * (*gamma*(F::ONE+*beta)+s[i-1]+*beta*s[i])
                * (*gamma*(F::ONE+*beta)+s[n+i-2]+*beta*s[n+i-1]);
            i += 1;
            ret.push(numer*F::invert(&denom).unwrap());
        }
        assert!(ret[n-1] == F::ONE);
        ret
    };
    values
}

/// send both table and polynomial (in monomial basis) for efficiency
fn compute_quotient_polynomial<F: PrimeField+WithSmallOrderMulGroup<3>>(
    g: &F,
    beta: &F,
    gamma: &F,
    delta: &F,
    t: &Vec<F>,
    h1: &Vec<F>,
    h2: &Vec<F>,
    z: &Vec<F>,
    t_poly: UnivariatePolynomial<F>,
    f_poly: UnivariatePolynomial<F>,
    h1_poly: UnivariatePolynomial<F>,
    h2_poly: UnivariatePolynomial<F>,
    z_poly: UnivariatePolynomial<F>,
) -> UnivariatePolynomial<F> {
    let n = t.len();
    let l0_poly = {
        let mut values = vec![F::ZERO; n];
        values[0] = F::ONE;
        UnivariatePolynomial::lagrange(values).ifft()
    };
    let ln_poly = {
        let mut values = vec![F::ZERO; n];
        values[n-1] = F::ONE;
        UnivariatePolynomial::lagrange(values).ifft()
    };
    let z_minus_1 = z_poly.clone() + (-F::ONE);
    let poly_a = &l0_poly * &z_minus_1;
    let t_shift_poly = to_shifted_poly(t);
    let z_shift_poly = to_shifted_poly(z);
    let h1_shift_poly = to_shifted_poly(h1);
    let h2_shift_poly = to_shifted_poly(h2);
    let poly_b = {
        let front = UnivariatePolynomial::monomial(vec![-(F::invert(g).unwrap()), F::ONE]); // x-1
        let beta_plus_1 = F::ONE + beta;
        let lhs = &z_poly * &beta_plus_1 * (f_poly + *gamma) *
            (&t_poly + &t_shift_poly * beta + beta_plus_1 * gamma);
        let rhs = z_shift_poly * (&h1_poly + &h1_shift_poly * beta + beta_plus_1 * gamma)
            * (&h2_poly + &h2_shift_poly * beta + beta_plus_1 * gamma);
        front * (&lhs - &rhs)
    };
    let poly_c = &ln_poly * &(&h1_poly - &h2_shift_poly);
    let poly_d = &ln_poly * &z_minus_1;
    let agg = aggregate_poly(delta, vec![&poly_a, &poly_b, &poly_c, &poly_d]);
    let vanish = { // x^n - 1
        let mut coeffs = vec![F::ZERO; n+1];
        coeffs[0] = -F::ONE;
        coeffs[n] = F::ONE;
        UnivariatePolynomial::monomial(coeffs)
    };
    agg.div(&vanish)
}

fn to_shifted_poly<F: PrimeField+WithSmallOrderMulGroup<3>>(table: &Vec<F>) -> UnivariatePolynomial<F> {
    let coeffs = table[1..]
        .iter()
        .chain(table[..1].iter())
        .cloned()
        .collect();
    UnivariatePolynomial::lagrange(coeffs).ifft()
}

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

#[cfg(test)]
mod tests {
    use crate::{
        backend::plookup::{
            preprocessor::{get_root_of_power_of_2_order, preprocess},
            prover::{
                prove,
                compute_inner_table,
                compute_quotient_polynomial,
                sorted_by_table,
            }, PlookupInfo
        },
        halo2_curves::{bn256::{Bn256, Fr}, ff::Field},
        pcs::{univariate::UnivariateKzg, PolynomialCommitmentScheme},
        poly::univariate::UnivariatePolynomial,
        util::{test::std_rng, transcript::{Keccak256Transcript,InMemoryTranscript}}
    };

    type Poly = UnivariatePolynomial<Fr>;
    type Pcs = UnivariateKzg<Bn256>;

    #[test]
    fn test_sorted_by_table() {
        let lookup = vec![Fr::from(1), Fr::from(1), Fr::from(2)];
        let table = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let sorted = sorted_by_table(&table, &lookup);
        let expected = vec![
            Fr::from(1), Fr::from(1), Fr::from(1), Fr::from(2),
            Fr::from(2), Fr::from(3), Fr::from(4)
        ];
        let is_the_same = {
            (sorted.len() == expected.len()) &&
            sorted.iter().zip(expected).all(|(a,b)| (*a==b))
        };
        assert!(is_the_same);
    }

    #[test]
    fn test_compute_inner_table() {
        let lookup = vec![Fr::from(1), Fr::from(1), Fr::from(2)];
        let table = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let sorted = sorted_by_table(&table, &lookup);
        let beta = Fr::from(1);
        let gamma = Fr::from(2);
        let z = compute_inner_table(
            &beta,
            &gamma,
            &lookup,
            &table,
            &sorted,
        );
        let ans = {
            let seven_over_eight = Fr::from(7)*Fr::invert(&Fr::from(8)).unwrap();
            vec![Fr::from(1), seven_over_eight, seven_over_eight, Fr::from(1)]
        };
        let is_the_same = {
            (z.len() == ans.len()) &&
            z.iter().zip(ans).all(|(a,b)| (*a==b))
        };
        assert!(is_the_same);
    }

    #[test]
    fn test_compute_quotient_polynomial() {
        let lookup = vec![Fr::from(1), Fr::from(1), Fr::from(2)];
        let table = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let sorted = sorted_by_table(&table, &lookup);
        let h1 = sorted[..4].to_vec();
        let h2 = sorted[4..].to_vec();
        let beta = Fr::from(1);
        let gamma = Fr::from(2);
        let t_poly = Poly::lagrange(table.clone());
        let f_poly = Poly::lagrange(lookup.clone());
        let h1_poly = Poly::lagrange(h1.clone());
        let h2_poly = Poly::lagrange(h2.clone());
        let z = compute_inner_table(
            &beta,
            &gamma,
            &lookup,
            &table,
            &sorted,
        );
        let z_poly = Poly::lagrange(z.clone());
        let delta = Fr::from(3);
        let g = get_root_of_power_of_2_order(2).unwrap();
        let poly = compute_quotient_polynomial(
            &g, &beta, &gamma, &delta, &table, &h1, &h2, &z,
            t_poly, f_poly, h1_poly, h2_poly, z_poly,
        );
        println!("{:#?}", poly);
    }

    #[test]
    fn test_prove() {
        let mut rng = std_rng();
        let k = 2;
        let n = 1 << k;
        let table = vec![Fr::one(), Fr::from(2), Fr::from(3), Fr::from(4)];
        let lookup = vec![Fr::one(), Fr::one(), Fr::from(2)];
        let param = Pcs::setup(n*4, 1, &mut rng).unwrap();
        let info = PlookupInfo{k, table, lookup};
        let (pp, _) = preprocess::<Fr, Pcs>(&param, &info).unwrap();
        let mut transcript = Keccak256Transcript::new(());
        prove::<Fr,Pcs>(pp, &info, &mut transcript).unwrap()
    }
}
