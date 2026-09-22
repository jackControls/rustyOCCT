//! Exact degree-two real expressions: (n + k*sqrt(radicand)) / d, with d>0.
//! No floating square root, rounded discriminant, or convergence threshold.
use crate::{exact, interval, Result, ScalarInterval};
use num_bigint::{BigInt, Sign};
use std::cmp::Ordering;

#[derive(Debug, Clone)]
pub(crate) struct Surd {
    n: BigInt,
    k: BigInt,
    radicand: BigInt,
    d: BigInt,
}

impl Surd {
    pub(crate) fn new(n: BigInt, k: BigInt, radicand: BigInt, d: BigInt) -> Self {
        debug_assert!(radicand.sign() != Sign::Minus && !exact::zero(&d));
        if d.sign() == Sign::Minus {
            Self {
                n: -n,
                k: -k,
                radicand,
                d: -d,
            }
        } else {
            Self { n, k, radicand, d }
        }
    }
    pub(crate) fn rational(n: BigInt, d: BigInt) -> Self {
        Self::new(n, BigInt::from(0), BigInt::from(0), d)
    }
    /// Return offset + multiplier*self, with exact integer affine coefficients.
    pub(crate) fn affine(&self, offset: &BigInt, multiplier: &BigInt) -> Self {
        Self::new(
            offset * &self.d + multiplier * &self.n,
            multiplier * &self.k,
            self.radicand.clone(),
            self.d.clone(),
        )
    }
    pub(crate) fn compare(&self, value: f64, scale: i32) -> Ordering {
        debug_assert!(value.is_finite() && (scale == 0 || scale == -1074));
        let shift = (scale + 1074) as usize;
        let rational = (&self.n << shift) - &self.d * exact::integer(value);
        let radical = &self.k << shift;
        sign_of_sum(&rational, &radical, &self.radicand)
    }
    pub(crate) fn bounds(&self, scale: i32, what: &'static str) -> Result<ScalarInterval> {
        interval::enclose(|x| self.compare(x, scale), what)
    }
}

fn sign(value: &BigInt) -> Ordering {
    value.cmp(&BigInt::from(0))
}

fn sign_of_sum(a: &BigInt, b: &BigInt, radicand: &BigInt) -> Ordering {
    if exact::zero(b) || exact::zero(radicand) {
        return sign(a);
    }
    if exact::zero(a) {
        return sign(b);
    }
    if a.sign() == b.sign() {
        return sign(a);
    }
    // With opposite signs, comparing nonnegative squared magnitudes is
    // equivalent to comparing |a| and |b|*sqrt(radicand), including equality.
    match (a * a).cmp(&(b * b * radicand)) {
        Ordering::Greater => sign(a),
        Ordering::Less => sign(b),
        Ordering::Equal => Ordering::Equal,
    }
}
