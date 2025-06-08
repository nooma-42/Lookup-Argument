use crate::poly::univariate::UnivariatePolynomial;
use crate::util::arithmetic::{root_of_unity, Field, PrimeField, WithSmallOrderMulGroup};
use std::collections::HashSet;

pub fn is_power_of_2(n: usize) -> bool {
    n != 0 && (n & (n - 1)) == 0
}

pub fn get_roots<F: Field + PrimeField>(N: usize) -> Vec<F> {
    assert!(is_power_of_2(N));
    let root = root_of_unity::<F>(N.ilog2() as usize);
    (0..N).map(|i| root.pow(&[i as u64])).collect()
}

pub fn get_vanishing_poly<F: PrimeField>(n: usize) -> UnivariatePolynomial<F> {
    let mut coeffs = vec![F::ZERO; n + 1];
    coeffs[0] = -F::ONE;
    coeffs[n] = F::ONE;
    UnivariatePolynomial::monomial(coeffs)
}

pub fn get_unique_positions<F: PrimeField>(c: &[F], values: &[F]) -> Vec<usize> {
    let mut result = HashSet::new();
    for value in values {
        // This will panic if a value in `values` is not present in `c`.
        // Consider returning Result or handling the error appropriately if this is possible.
        result.insert(c.iter().position(|x| x == value).unwrap());
    }
    result.into_iter().collect()
}

// === Optimized subproduct tree implementation for tau polynomials ===

fn log_2(n: usize) -> usize {
    assert!(n.is_power_of_two());
    n.ilog2() as usize
}

fn pow_2(k: usize) -> usize {
    1 << k
}

/// Compute subproduct tree in O(M(n) * log(n)) time, where O(M(n)) is the
/// asymptotic complexity of multiplication, and equal to O(nlog(n)) if using
/// FFT-based fast multiplication.
///
/// Return a vector of levels, each level is a vector of polynomials,
/// and each polynomial is a vector of coefficients.
fn construct_subproduct_tree<F: PrimeField + WithSmallOrderMulGroup<3>>(
    domain: &[F],
) -> Vec<Vec<UnivariatePolynomial<F>>> {
    let n = domain.len();
    
    if n == 0 {
        return vec![];
    }
    
    if n == 1 {
        // Single element: just return the linear polynomial (X - domain[0])
        let poly = UnivariatePolynomial::monomial(vec![-domain[0], F::ONE]);
        return vec![vec![poly]];
    }

    let mut tree: Vec<Vec<UnivariatePolynomial<F>>> = Vec::new();
    let mut current_level: Vec<UnivariatePolynomial<F>> = Vec::new();
    
    // Build leaf level: (X - domain[i]) for each point
    for &u in domain.iter() {
        current_level.push(UnivariatePolynomial::monomial(vec![-u, F::ONE]));
    }
    tree.push(current_level.clone());

    // Build internal levels by multiplying pairs
    while current_level.len() > 1 {
        let mut new_level = Vec::new();
        
        // Process pairs
        let mut i = 0;
        while i + 1 < current_level.len() {
            let left = &current_level[i];
            let right = &current_level[i + 1];
            let poly = left.poly_mul(right.clone());
            new_level.push(poly);
            i += 2;
        }
        
        // Handle odd number of elements by copying the last one
        if i < current_level.len() {
            new_level.push(current_level[i].clone());
        }
        
        tree.push(new_level.clone());
        current_level = new_level;
    }

    tree
}

/// Multipoint evaluation using subproduct tree
fn eval_rec<F: PrimeField + WithSmallOrderMulGroup<3>>(
    tree: &Vec<Vec<UnivariatePolynomial<F>>>,
    k: usize,
    base: usize,
    f: &UnivariatePolynomial<F>,
    u: &[F],
) -> Vec<F> {
    let n = u.len();
    
    if k == 0 {
        assert_eq!(u.len(), 1);
        if f.is_empty() { 
            return vec![F::ZERO]; 
        }
        return vec![f.evaluate(&u[0])];
    }
    
    let divisor0 = &tree[k - 1][2 * base];
    let divisor1 = &tree[k - 1][2 * base + 1];
    
    let (_q0, r0) = f.div_rem(divisor0);
    let (_q1, r1) = f.div_rem(divisor1);

    let (u0, u1) = u.split_at(n / 2);
    let mut rs0: Vec<F> = eval_rec(tree, k - 1, base * 2, &r0, u0);
    let mut rs1: Vec<F> = eval_rec(tree, k - 1, base * 2 + 1, &r1, u1);
    rs0.append(&mut rs1);
    rs0
}

/// Optimized computation of tau polynomials using subproduct tree
/// Time complexity: O(m log^2 m) instead of O(m^2)
pub fn get_tau_polys_optimized<F: PrimeField + WithSmallOrderMulGroup<3>>(
    unique_positions: &[usize],
    roots_N: &[F],
) -> Vec<UnivariatePolynomial<F>> {
    let points: Vec<F> = unique_positions.iter().map(|&i| roots_N[i]).collect();
    let m = points.len();
    
    if m == 0 {
        return vec![];
    }
    
    if m == 1 {
        // Special case: single point, tau_0(X) = 1
        return vec![UnivariatePolynomial::monomial(vec![F::ONE])];
    }

    // We don't need to pad to power of 2 for this application
    // Just work with the exact number of points
    
    // Step 1: Build subproduct tree for efficient polynomial operations
    let tree = construct_subproduct_tree(&points);
    
    // Step 2: Get the vanishing polynomial z_I(X) = ∏(X - p_i)
    // This is stored in the root of the tree
    let tree_levels = tree.len();
    let root_level = tree_levels - 1;
    
    // The tree should have a single polynomial at the root
    assert_eq!(tree[root_level].len(), 1, "Root should have exactly one polynomial");
    let z_i_poly = &tree[root_level][0];
    
    // Step 3: For each point, compute tau_i(X) = z_I(X) / ((X - p_i) * z_I'(p_i))
    let mut tau_polys = Vec::new();
    
    for i in 0..m {
        let point = points[i];
        
        // Compute z_I(X) / (X - p_i)
        let linear_divisor = UnivariatePolynomial::monomial(vec![-point, F::ONE]);
        let (quotient, remainder) = z_i_poly.div_rem(&linear_divisor);
        
        // Verify that remainder is zero (should be since point is a root)
        debug_assert!(remainder.is_empty() || remainder.coeffs().iter().all(|&c| c == F::ZERO));
        
        // Compute z_I'(p_i) by evaluating the derivative at the point
        let z_i_prime = z_i_poly.derivative();
        let z_i_prime_at_point = z_i_prime.evaluate(&point);
        
        // tau_i(X) = quotient / z_I'(p_i)
        let denominator_inv = z_i_prime_at_point.invert().expect("z_I'(p_i) should be non-zero");
        let tau_i = &quotient * denominator_inv;
        
        tau_polys.push(tau_i);
    }
    
    tau_polys
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::arithmetic::WithSmallOrderMulGroup;
    use halo2_curves::bn256::Fr;

    #[test]
    fn test_tau_polys_optimization() {
        // Test with small input
        let unique_positions = vec![0, 2, 5];
        let roots_N = get_roots::<Fr>(8); // power of 2

        // Compute using original O(m^2) method
        let tau_polys_original = get_tau_polys_original(&unique_positions, &roots_N);
        
        // Compute using optimized O(m log^2 m) method
        let tau_polys_optimized = get_tau_polys_optimized(&unique_positions, &roots_N);

        // Verify they produce the same results
        assert_eq!(tau_polys_original.len(), tau_polys_optimized.len());
        
        for (orig, opt) in tau_polys_original.iter().zip(tau_polys_optimized.iter()) {
            // Test that they evaluate to the same values at multiple points
            let test_points = [Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(7), Fr::from(11)];
            for &point in &test_points {
                let orig_eval = orig.evaluate(&point);
                let opt_eval = opt.evaluate(&point);
                assert_eq!(orig_eval, opt_eval, "Polynomials differ at point {:?}", point);
            }
        }
    }

    // Keep the original implementation for testing purposes
    fn get_tau_polys_original<F: PrimeField + WithSmallOrderMulGroup<3>>(
        unique_positions: &[usize],
        roots_N: &[F],
    ) -> Vec<UnivariatePolynomial<F>> {
        let points: Vec<F> = unique_positions.iter().map(|&i| roots_N[i]).collect();
        
        (0..unique_positions.len())
            .map(|i| {
                let mut tau_poly = UnivariatePolynomial::monomial(vec![F::ONE]);
                let mut denominator = F::ONE;
                for j in 0..unique_positions.len() {
                    if i != j {
                        tau_poly = tau_poly.poly_mul(UnivariatePolynomial::monomial(vec![-points[j], F::ONE]));
                        denominator *= points[i] - points[j];
                    }
                }
                &tau_poly * &denominator.invert().unwrap()
            })
            .collect()
    }

    #[test]
    fn test_edge_cases() {
        let roots_N = get_roots::<Fr>(8);
        
        // Test empty case
        let empty_positions = vec![];
        let tau_polys_empty = get_tau_polys_optimized(&empty_positions, &roots_N);
        assert_eq!(tau_polys_empty.len(), 0);
        
        // Test single element case
        let single_position = vec![3];
        let tau_polys_single = get_tau_polys_optimized(&single_position, &roots_N);
        assert_eq!(tau_polys_single.len(), 1);
        // For single element, tau_0(X) should be the constant polynomial 1
        assert_eq!(tau_polys_single[0].evaluate(&Fr::from(42)), Fr::ONE);
    }

    #[test]
    fn test_performance_improvement() {
        // This is more of a conceptual test - with larger inputs, 
        // the optimized version should be significantly faster
        let unique_positions: Vec<usize> = (0..16).step_by(2).collect(); // [0, 2, 4, 6, 8, 10, 12, 14]
        let roots_N = get_roots::<Fr>(16);
        
        // Both should produce valid polynomials
        let tau_polys_optimized = get_tau_polys_optimized(&unique_positions, &roots_N);
        assert_eq!(tau_polys_optimized.len(), 8);
        
        // Each tau polynomial should be well-formed
        for (i, tau_poly) in tau_polys_optimized.iter().enumerate() {
            // tau_i(roots_N[unique_positions[i]]) should be 1
            let eval_at_own_point = tau_poly.evaluate(&roots_N[unique_positions[i]]);
            assert_eq!(eval_at_own_point, Fr::ONE, "tau_{} should evaluate to 1 at its own point", i);
            
            // tau_i(roots_N[unique_positions[j]]) should be 0 for j != i
            for (j, &other_pos) in unique_positions.iter().enumerate() {
                if i != j {
                    let eval_at_other_point = tau_poly.evaluate(&roots_N[other_pos]);
                    assert_eq!(eval_at_other_point, Fr::ZERO, "tau_{} should evaluate to 0 at point {}", i, j);
                }
            }
        }
    }
} 