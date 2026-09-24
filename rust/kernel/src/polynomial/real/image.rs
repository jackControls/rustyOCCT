//! Exact rational-function images in a square-free polynomial quotient ring.
//! Keep common denominators and carry the complete dependence witness through
//! every integer row operation. See MATHEMATICS.md for the algebraic contract.
use super::{abs, content, gcd_integer, isolate, one, AlgebraicRoot, Budget, IntPolynomial, R};
use crate::{exact, Error, Result};
use num_bigint::{BigInt, Sign};

pub(crate) struct ImageBudget {
    remaining: usize,
    max_bits: u64,
}
impl ImageBudget {
    pub(crate) fn new(updates: usize, max_bits: u64) -> Self {
        Self {
            remaining: updates,
            max_bits,
        }
    }
    fn charge(&mut self, updates: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(updates)
            .ok_or(Error::ComputationLimit(
                "algebraic image coefficient updates",
            ))?;
        Ok(())
    }
    fn check<'a>(&self, values: impl IntoIterator<Item = &'a BigInt>) -> Result<()> {
        if values.into_iter().any(|x| x.bits() > self.max_bits) {
            Err(Error::ComputationLimit("algebraic image coefficient size"))
        } else {
            Ok(())
        }
    }
}

#[derive(Clone)]
struct RationalPolynomial {
    coefficients: Vec<BigInt>,
    denominator: BigInt,
}
fn trim(v: &mut Vec<BigInt>) {
    while v.last().is_some_and(exact::zero) {
        v.pop();
    }
}
impl RationalPolynomial {
    fn new(
        coefficients: Vec<BigInt>,
        denominator: BigInt,
        budget: &mut ImageBudget,
    ) -> Result<Self> {
        budget.charge(coefficients.len() + 1)?;
        budget.check(coefficients.iter().chain([&denominator]))?;
        let mut result = Self {
            coefficients,
            denominator,
        };
        trim(&mut result.coefficients);
        let mut common = content(&result.coefficients, abs(&result.denominator));
        if result.denominator.sign() == Sign::Minus {
            common = -common;
        }
        for c in &mut result.coefficients {
            *c /= &common;
        }
        result.denominator /= common;
        Ok(result)
    }
    fn from_rationals(values: &[R], budget: &mut ImageBudget) -> Result<Self> {
        let mut den = BigInt::from(1);
        for x in values {
            budget.charge(1)?;
            den = (&den / gcd_integer(den.clone(), x.denom().clone())) * x.denom();
            budget.check([&den])?;
        }
        Self::new(
            values
                .iter()
                .map(|x| x.numer() * (&den / x.denom()))
                .collect(),
            den,
            budget,
        )
    }
    fn product(&self, other: &Self, budget: &mut ImageBudget) -> Result<Self> {
        if self.coefficients.is_empty() || other.coefficients.is_empty() {
            return Self::new(vec![], BigInt::from(1), budget);
        }
        budget.charge(self.coefficients.len() * other.coefficients.len())?;
        let mut result =
            vec![BigInt::from(0); self.coefficients.len() + other.coefficients.len() - 1];
        for (i, a) in self.coefficients.iter().enumerate() {
            for (j, b) in other.coefficients.iter().enumerate() {
                result[i + j] += a * b;
            }
            budget.check(&result)?;
        }
        Self::new(result, &self.denominator * &other.denominator, budget)
    }
    fn remainder(mut self, p: &IntPolynomial, budget: &mut ImageBudget) -> Result<Self> {
        let lead = p.0.last().expect("nonzero modulus");
        debug_assert!(lead.sign() == Sign::Plus);
        while self.coefficients.len() >= p.0.len() {
            budget.charge(self.coefficients.len() + p.0.len())?;
            let shift = self.coefficients.len() - p.0.len();
            let scale = self.coefficients.last().unwrap().clone();
            for c in &mut self.coefficients {
                *c *= lead;
            }
            for (i, c) in p.0.iter().enumerate() {
                self.coefficients[i + shift] -= &scale * c;
            }
            self = Self::new(self.coefficients, self.denominator * lead, budget)?;
        }
        Ok(self)
    }
}

// r = witness * input_denominator modulo the defining polynomial. Removing
// content from r must divide its witness by the SAME content.
fn normalize_pair(
    mut r: Vec<BigInt>,
    witness: RationalPolynomial,
    budget: &mut ImageBudget,
) -> Result<(IntPolynomial, RationalPolynomial)> {
    trim(&mut r);
    budget.charge(r.len())?;
    budget.check(&r)?;
    let common = content(&r, BigInt::from(0)).max(BigInt::from(1));
    for c in &mut r {
        *c /= &common;
    }
    let witness =
        RationalPolynomial::new(witness.coefficients, witness.denominator * common, budget)?;
    Ok((IntPolynomial(r), witness))
}

fn extended_gcd(
    d: &RationalPolynomial,
    p: &IntPolynomial,
    budget: &mut ImageBudget,
) -> Result<(IntPolynomial, RationalPolynomial)> {
    let (mut a, mut av) = normalize_pair(
        d.coefficients.clone(),
        RationalPolynomial {
            coefficients: vec![d.denominator.clone()],
            denominator: BigInt::from(1),
        },
        budget,
    )?;
    let mut b = p.clone();
    let mut bv = RationalPolynomial {
        coefficients: vec![],
        denominator: BigInt::from(1),
    };
    while !b.is_zero() {
        let (mut r, mut v) = (a, av);
        while !r.is_zero() && r.0.len() >= b.0.len() {
            let shift = r.0.len() - b.0.len();
            let lead = b.0.last().unwrap();
            let magnitude = abs(lead);
            let factor = if lead.sign() == Sign::Minus {
                -r.0.last().unwrap()
            } else {
                r.0.last().unwrap().clone()
            };
            budget.charge(r.0.len() + b.0.len() + v.coefficients.len() + bv.coefficients.len())?;
            for c in &mut r.0 {
                *c *= &magnitude;
            }
            for (i, c) in b.0.iter().enumerate() {
                r.0[i + shift] -= &factor * c;
            }
            let common = (&v.denominator
                / gcd_integer(v.denominator.clone(), bv.denominator.clone()))
                * &bv.denominator;
            let left = &magnitude * (&common / &v.denominator);
            let right = &factor * (&common / &bv.denominator);
            for c in &mut v.coefficients {
                *c *= &left;
            }
            v.coefficients.resize(
                v.coefficients.len().max(bv.coefficients.len() + shift),
                BigInt::from(0),
            );
            for (i, c) in bv.coefficients.iter().enumerate() {
                v.coefficients[i + shift] -= &right * c;
            }
            (r, v) = normalize_pair(
                r.0,
                RationalPolynomial {
                    coefficients: v.coefficients,
                    denominator: common,
                },
                budget,
            )?;
        }
        (a, av, b, bv) = (b, bv, r, v);
    }
    if a.0.last().unwrap().sign() == Sign::Minus {
        a = a.negate();
        for c in &mut av.coefficients {
            *c = -&*c;
        }
    }
    Ok((a, av))
}

/// All real image roots, not just the selected parameter's image. The caller
/// selects its unique image with exact rational threshold comparisons and may
/// reuse this result for every parameter with the same defining polynomial.
pub(crate) fn image_roots(
    root: &AlgebraicRoot,
    numerator: &[R],
    denominator: &[R],
    roots: &mut Budget,
    budget: &mut ImageBudget,
) -> Result<Vec<AlgebraicRoot>> {
    if root.vanishes_polynomial(&IntPolynomial::from_rationals(denominator)) {
        return Err(Error::Degenerate("algebraic image denominator"));
    }
    let image = image_polynomial(&root.defining.polynomial, numerator, denominator, budget)?;
    let lead = abs(image.0.last().unwrap());
    let bound = image.0[..image.0.len() - 1]
        .iter()
        .map(|c| R::new(abs(c), lead.clone()))
        .max()
        .unwrap()
        + one();
    isolate(&image, -&bound, bound, roots)
}

/// Minimal annihilating polynomial of N/D over all non-pole roots of the
/// square-free defining polynomial. Used independently of isolation in tests.
fn image_polynomial(
    defining: &IntPolynomial,
    numerator: &[R],
    denominator: &[R],
    budget: &mut ImageBudget,
) -> Result<IntPolynomial> {
    let mut p = defining.clone();
    budget.check(&p.0)?;
    let d = RationalPolynomial::from_rationals(denominator, budget)?;
    if d.coefficients.is_empty() {
        return Err(Error::Degenerate("algebraic image denominator"));
    }
    let (common, mut inverse) = extended_gcd(&d, &p, budget)?;
    if !common.is_constant() {
        // Complex or other real poles of P can make a raw resultant vanish,
        // even though the selected root is not a pole. Remove just those factors.
        budget.charge(p.0.len() * common.0.len())?;
        p = p.quotient_exact(&common).positive();
        if p.is_constant() {
            return Err(Error::Degenerate("algebraic image denominator"));
        }
        let (g, v) = extended_gcd(&d, &p, budget)?;
        debug_assert!(g.is_constant());
        inverse = v;
    }
    let n = RationalPolynomial::from_rationals(numerator, budget)?;
    let image = n.product(&inverse, budget)?.remainder(&p, budget)?;
    let degree = p.0.len() - 1;
    let mut powers = RationalPolynomial {
        coefficients: vec![BigInt::from(1)],
        denominator: BigInt::from(1),
    };
    let mut pivots: Vec<Option<Vec<BigInt>>> = vec![None; degree];
    for k in 0..=degree {
        let mut row = powers.coefficients.clone();
        row.resize(2 * degree + 1, BigInt::from(0));
        row[degree + k] = powers.denominator.clone();
        let mut independent = false;
        for j in 0..degree {
            if exact::zero(&row[j]) {
                continue;
            }
            let Some(base) = &pivots[j] else {
                pivots[j] = Some(row.clone());
                independent = true;
                break;
            };
            budget.charge(3 * row.len())?;
            let common = gcd_integer(abs(&row[j]), abs(&base[j]));
            let (left, right) = (&base[j] / &common, &row[j] / &common);
            for (c, b) in row.iter_mut().zip(base) {
                *c = &left * &*c - &right * b;
            }
            budget.check(&row)?;
            let common = content(&row, BigInt::from(0));
            if common > BigInt::from(1) {
                for c in &mut row {
                    *c /= &common;
                }
            }
        }
        if !independent {
            let image = IntPolynomial::new(row[degree..].to_vec()).positive();
            debug_assert_eq!(image.0.len(), k + 1);
            return Ok(image);
        }
        powers = powers.product(&image, budget)?.remainder(&p, budget)?;
    }
    unreachable!("degree+1 quotient-ring powers are linearly dependent")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn every_image_coefficient_matches_independent_sympy_resultants() {
        let mut count = 0;
        for line in include_str!("../../../../fixtures/algebraic-images.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let mut words = line.split_whitespace();
            let name = words.next().unwrap();
            let mut read = || -> Vec<R> {
                let n: usize = words.next().unwrap().parse().unwrap();
                (0..n)
                    .map(|_| R::from_str(words.next().unwrap()).unwrap())
                    .collect()
            };
            let (p, n, d, expected) = (read(), read(), read(), read());
            assert!(words.next().is_none(), "{name}");
            let actual = image_polynomial(
                &IntPolynomial::from_rationals(&p).positive(),
                &n,
                &d,
                &mut ImageBudget::new(2_000_000, 65_536),
            )
            .unwrap_or_else(|e| panic!("{name}: {e}"));
            assert_eq!(
                actual,
                IntPolynomial::from_rationals(&expected).positive(),
                "{name}"
            );
            count += 1;
        }
        assert_eq!(count, 94);
    }

    #[test]
    fn poles_and_image_work_limits_return_errors() {
        let p = IntPolynomial::new(vec![BigInt::from(-2), BigInt::from(0), BigInt::from(1)]);
        let n = [one()];
        let d = [one()];
        for limits in [(0, 65_536), (2_000_000, 0)] {
            assert!(matches!(
                image_polynomial(&p, &n, &d, &mut ImageBudget::new(limits.0, limits.1)),
                Err(Error::ComputationLimit(_))
            ));
        }
        for d in [
            vec![],
            vec![
                R::from_integer(BigInt::from(-2)),
                super::super::zero(),
                one(),
            ],
        ] {
            assert!(matches!(
                image_polynomial(&p, &n, &d, &mut ImageBudget::new(2_000_000, 65_536)),
                Err(Error::Degenerate(_))
            ));
        }
        let root = AlgebraicRoot::rational(one());
        assert!(matches!(
            image_roots(
                &root,
                &n,
                &[-one(), one()],
                &mut Budget::new(Default::default()),
                &mut ImageBudget::new(2_000_000, 65_536)
            ),
            Err(Error::Degenerate(_))
        ));
    }
}
