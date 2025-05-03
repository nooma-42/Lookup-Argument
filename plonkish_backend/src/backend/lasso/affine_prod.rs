use generic_array::typenum::True;
use halo2_curves::bn256::Fr;
use std::cmp;
use std::collections::HashSet;
use std::fmt;
use std::ops::Mul;
use std::sync::Arc;

type Scalar = Fr;

#[derive(Clone, Debug, PartialEq)]
pub struct AffineTerm {
    pub coeff: Scalar,
    pub x_i: usize,
    pub constant: Scalar,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AffineProduct {
    pub coeff: Scalar,
    pub terms: Vec<AffineTerm>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AffinePolynomial {
    pub terms: Vec<AffineProduct>,
    pub constant: Scalar,
}

#[derive(Clone, Debug, PartialEq)]
pub struct UnivariateExpansion {
    pub coeffs: Vec<Scalar>,
    pub degree: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MultivariateExpansion {
    pub terms: Vec<Vec<Scalar>>,
    pub v: usize,
}

// affine term = coeff_i * x_i + const_i
impl AffineTerm {
    pub fn new(coeff: Scalar, i: usize, constant: Scalar) -> Self {
        if i < 1 {
            panic!("i should be greater than 0");
        }
        Self {
            coeff,
            x_i: i,
            constant,
        }
    }

    pub fn eval(&self, x: Scalar) -> Scalar {
        self.coeff * x + self.constant
    }

    pub fn is_constant(&self) -> bool {
        self.coeff == Scalar::zero()
    }

    pub fn convert(&self) -> UnivariateExpansion {
        UnivariateExpansion::new(vec![self.constant, self.coeff], 1)
    }
}

impl Mul<Scalar> for AffineTerm {
    type Output = Self;
    fn mul(self, rhs: Scalar) -> Self {
        Self {
            coeff: self.coeff * rhs,
            x_i: self.x_i,
            constant: self.constant * rhs,
        }
    }
}

impl fmt::Display for AffineTerm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} * x_{} + {})", self.coeff, self.x_i, self.constant)
    }
}

// γ * (term₁ * term₂ * ...) → γ * ((α₁ * x_i₁ + β₁) * (α₂ * x_i₂ + β₂) * ...)
impl AffineProduct {
    pub fn new(coeff: Scalar, terms: Vec<AffineTerm>) -> Self {
        Self { coeff, terms }
    }

    pub fn mult(&mut self, n: Scalar) {
        self.coeff *= n;
    }

    pub fn apply(&self) -> Result<Scalar, Self> {
        let mut res = self.coeff;
        let mut new_terms = vec![];

        if res == Scalar::zero() {
            return Ok(Scalar::zero());
        }

        for t in &self.terms {
            if t.coeff == Scalar::zero() {
                if t.constant == Scalar::zero() {
                    return Ok(Scalar::zero());
                }
                res *= t.constant;
            } else {
                new_terms.push(t.clone());
            }
        }

        if new_terms.is_empty() {
            Ok(res)
        } else {
            Err(AffineProduct::new(res, new_terms))
        }
    }

    pub fn eval_univariate(&self, x: Scalar) -> Scalar {
        let mut res = self.coeff;
        for t in &self.terms {
            let val = t.eval(x);
            if val == Scalar::zero() {
                return Scalar::zero();
            }
            res *= val;
        }
        res
    }
    pub fn get_expansion(&self) -> UnivariateExpansion {
        let mut res = self.terms[0].convert() * self.coeff;
        for t in self.terms.iter().skip(1) {
            res = res * t.clone();
        }
        res
    }

    pub fn get_multi_expansion(&self, v: usize) -> MultivariateExpansion {
        let mut term = [Scalar::zero(); v + 1];
        tern[0] = self.coeff.clone();
        let mut mexp = MultivariateExpansion::new(vec![term], v);
        for t in &self.terms {
            mexp = mexp * t.clone();
        }
        mexp
    }
}

impl Mul for AffineProduct {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let mut combined_terms = self.terms.clone();
        combined_terms.extend(rhs.terms.clone());
        AffineProduct::new(self.coeff * rhs.coeff, combined_terms)
    }
}

impl fmt::Display for AffineProduct {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let terms_str = self
            .terms
            .iter()
            .map(|t| format!("{}", t))
            .collect::<Vec<_>>()
            .join(" * ");
        write!(f, "{} * ({})", self.coeff, terms_str)
    }
}

// P(x) = Σ γᵢ * Π (α * x_j + β) + c
impl AffinePolynomial {
    pub fn new(terms: Vec<AffineProduct>, constant: Scalar) -> Self {
        Self { terms, constant }
    }

    pub fn eval_i(&self, x_i: Scalar, i: usize) -> Self {
        assert!(i > 0, "i should start from 1");
        let mut new_terms = vec![];
        let mut new_const = self.constant;

        for mono in &self.terms {
            let mut reduced_terms = vec![];
            let mut coeff = mono.coeff;

            for term in &mono.terms {
                if term.x_i == i {
                    let val = term.eval(x_i);
                    if val == Scalar::zero() {
                        reduced_terms.clear();
                        coeff = Scalar::zero();
                        break;
                    } else {
                        coeff *= val;
                    }
                } else {
                    reduced_terms.push(term.clone());
                }
            }

            if reduced_terms.is_empty() {
                new_const += coeff;
            } else {
                new_terms.push(AffineProduct::new(coeff, reduced_terms));
            }
        }
        AffinePolynomial::new(new_terms, new_const).apply_all()
    }

    pub fn eval(&self, x: &[Scalar]) -> Scalar {
        let mut poly = self.clone();
        for i in 0..x.len() {
            poly = poly.eval_i(x[i], i + 1);
        }
        poly.constant
    }

    pub fn is_univariate(&self) -> bool {
        let mut seen = None;
        for mono in &self.terms {
            for term in &mono.terms {
                match seen {
                    None => seen = Some(term.x_i),
                    Some(v) if v != term.x_i => return false,
                    _ => {}
                }
            }
        }
        seen.is_some()
    }
    pub fn apply_all(&self) -> Self {
        let mut new_terms = vec![];
        let mut new_const = self.constant;
        for t in &self.terms {
            match t.apply() {
                TermOrScalar::Scalar(c) => new_const += c,
                TermOrScalar::AffineProduct(m) => new_terms.push(m),
            }
        }
        Self::new(new_terms, new_const)
    }

    pub fn eval_univariate(&self, x: Scalar) -> Scalar {
        let mut res = Scalar::zero();
        for t in &self.terms {
            res += t.eval_univariate(x);
        }
        res + self.constant
    }
    pub fn get_highest_degree(&self) -> usize {
        let mut max_degree = 0;
        for term in &self.terms {
            max_degree = cmp(term.terms.len(), max_degree);
        }
        max_degree
    }

    pub fn get_highest_index(&self) -> usize {
        let mut max_index = 0;
        for mono in &self.terms {
            for term in &mono.terms {
                max_index = cmp(term.x_i, max_index);
            }
        }
        max_index
    }

    pub fn get_expansion(&self) -> UnivariateExpansion {
        let mut res = UnivariateExpansion::new(vec![Scalar::zero()], 0);
        for t in &self.terms {
            res = res + t.get_expansion();
        }
        res
    }

    pub fn get_multi_expansion(&self, v: usize) -> MultivariateExpansion {
        let mut term = [Scalar::zero(); v + 1];
        tern[0] = self.coeff.clone();
        let mut mexp = MultivariateExpansion::new(vec![term], v);
        for mono in &self.terms {
            mexp = mexp * mono.get_multi_expansion(v);
        }
        mexp
    }

    pub fn quotient_single_term(&self, value: Scalar, i: usize) -> Self {
        assert!(i > 0, "i should start from 1");
        let mut new_terms_poly = vec![];
        let mut new_const = Scalar::zero();

        for mono in &self.terms {
            let mut terms = mono.terms.clone();
            let mut coeff = mono.coeff;

            while true {
                let mut new_terms = Vec::new();
                let mut this_coeff = coeff;
                let mut next_coeff = coeff;
                let mut has_x_i = false;
                for term in terms {
                    if has_x_i {
                        new_terms.push(term);
                    } else {
                        if term.x_i == i {
                            has_x_i = true;
                            this_coeff *= term.coeff;
                            next_coeff += term.eval(value);
                        } else {
                            new_terms.push(term);
                        }
                    }
                }
                if has_x_i == false {
                    break;
                }
                if new_terms.len() == 0 {
                    new_const += this_coeff;
                    break;
                } else {
                    new_terms_poly.push(AffineProduct::new(this_coeff, new_terms));
                    terms = new_terms;
                    coeff = next_coeff;
                }
            }
        }
        AffinePolynomial::new(new_terms_poly, new_const).apply_all()
    }

    // Todo: Test this function since it is implemented as a copy of one in pylookup.
    pub fn evaluate(&self, args: &[Scalar]) -> Scalar {
        // Evaluate the polynomial at given point(s).
        // - param args: either
        //      - a single Scalar for univariate polynomials, or
        //      - a list of Scalars for multivariate polynomials.
        // - return: The result of the polynomial evaluation.
        let max_x_i = self
            .terms
            .iter()
            .flat_map(|mono| mono.terms.iter().map(|t| t.x_i))
            .max()
            .unwrap_or(0);

        assert!(
            args.len() >= max_x_i,
            "Not enough arguments provided. Expected at least {}.",
            max_x_i
        );

        let mut result = self.constant;
        for mono in &self.terms {
            let mut term_result = mono.coeff;
            for term in &mono.terms {
                let x_val = args[term.x_i - 1].clone();
                term_result *= term.coeff * x_val + term.constant;
            }
            result += term_result;
        }
        result
    }

    pub fn num_var(&self) -> usize {
        let mut variables = HashSet::new();
        for mono in &self.terms {
            for term in &mono.terms {
                variables.insert(term.x_i);
            }
        }
        variables.len()
    }
}

impl fmt::Display for AffinePolynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let terms_str: Vec<String> = self.terms.iter().map(|t| t.to_string()).collect();
        write!(f, "{} + {}", terms_str.join(" + "), self.constant)
    }
}

// PartialEq is derived
impl PartialEq for AffinePolynomial {
    fn eq(&self, other: &Self) -> bool {
        if self.terms.len() != other.terms.len() {
            return false;
        }

        for (a, b) in self.terms.iter().zip(&other.terms) {
            if a != b {
                return false;
            }
        }

        self.constant == other.constant
    }
}

impl UnivariateExpansion {
    pub fn new(coeffs: Vec<Scalar>, degree: usize) -> Self {
        Self { coeffs, degree }
    }
    // ToDo: Implement the following parts of this Struct
}

impl MultivariateExpansion {
    pub fn new(terms: Vec<Vec<Scalar>>, v: usize) -> Self {
        Self { terms, v }
    }
    // ToDo: Implement the following parts of this Struct
}
