//! Certified real roots of finite binary64 polynomials of degree at most two.
//!
//! OCCT reference: math_DirectPolynomialRoots::Solve(A,B,C), including degree
//! reduction, multiplicity and stable handling of cancellation. This module
//! deliberately uses exact zero/discriminant decisions, with no coefficient
//! threshold or near-double-root merging. Cubic/quartic roots are not included.
use crate::{algebraic::Surd, exact, math::finite, Result, ScalarInterval};
use num_bigint::{BigInt, Sign};
use std::cmp::Ordering;

/// Retains the exact algebraic root, even when no finite binary64 enclosure
/// exists. A rounded interval is a view of this root, not its identity.
#[derive(Debug, Clone)]
pub struct QuadraticRoot {
    pub(crate) value: Surd,
    multiplicity: u8,
}

impl QuadraticRoot {
    pub fn multiplicity(&self) -> u8 {
        self.multiplicity
    }
    pub fn bounds(&self) -> Result<ScalarInterval> {
        self.value.bounds(0, "polynomial root")
    }
    /// Exact ordering relative to any finite binary64 value.
    pub fn compare(&self, value: f64) -> Result<Ordering> {
        finite(value, "root comparison value")?;
        Ok(self.value.compare(value, 0))
    }
    pub(crate) fn rational(num: BigInt, den: BigInt, multiplicity: u8) -> Self {
        Self {
            value: Surd::rational(num, den),
            multiplicity,
        }
    }
}

#[derive(Debug, Clone)]
pub enum QuadraticRoots {
    None,
    /// The identically zero polynomial: every real number is a root.
    All,
    One(QuadraticRoot),
    /// Two distinct real roots, in exact increasing order. Their binary64
    /// enclosures may overlap; callers must not merge them on that basis.
    Two {
        lower: QuadraticRoot,
        upper: QuadraticRoot,
    },
}

/// Roots of a*x² + b*x + c = 0 for all finite coefficients. Degree reduction
/// uses exact zero. A double root is returned once with multiplicity two.
pub fn quadratic_roots(a: f64, b: f64, c: f64) -> Result<QuadraticRoots> {
    for value in [a, b, c] {
        finite(value, "polynomial coefficient")?;
    }
    Ok(solve_integer(
        exact::integer(a),
        exact::integer(b),
        exact::integer(c),
    ))
}

pub(crate) fn solve_integer(mut a: BigInt, mut b: BigInt, mut c: BigInt) -> QuadraticRoots {
    // Remove only an exactly common power of two. Typical binary64 inputs have
    // a large shared lattice factor; retaining it would waste integer work.
    if let Some(shift) = [&a, &b, &c].iter().filter_map(|x| x.trailing_zeros()).min() {
        a >>= shift as usize;
        b >>= shift as usize;
        c >>= shift as usize;
    }
    if exact::zero(&a) {
        return if !exact::zero(&b) {
            QuadraticRoots::One(QuadraticRoot::rational(-c, b, 1))
        } else if exact::zero(&c) {
            QuadraticRoots::All
        } else {
            QuadraticRoots::None
        };
    }
    if a.sign() == Sign::Minus {
        a = -a;
        b = -b;
        c = -c;
    }
    let discriminant = &b * &b - BigInt::from(4) * &a * c;
    match discriminant.sign() {
        Sign::Minus => QuadraticRoots::None,
        Sign::NoSign => QuadraticRoots::One(QuadraticRoot::rational(-b, a * 2, 2)),
        Sign::Plus => {
            let numerator = -b;
            let denominator = a << 1_usize;
            QuadraticRoots::Two {
                lower: QuadraticRoot {
                    value: Surd::new(
                        numerator.clone(),
                        BigInt::from(-1),
                        discriminant.clone(),
                        denominator.clone(),
                    ),
                    multiplicity: 1,
                },
                upper: QuadraticRoot {
                    value: Surd::new(numerator, BigInt::from(1), discriminant, denominator),
                    multiplicity: 1,
                },
            }
        }
    }
}
