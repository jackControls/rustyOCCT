//! Exact real-root isolation and sign evaluation by primitive Sturm sequences.
//!
//! Polynomials are primitive integer vectors; pseudo-division scales by positive
//! leading-coefficient magnitudes so sign variations remain valid. Sturm-Tarski
//! queries decide polynomial signs at an isolated algebraic root, including zero.
//! See MATHEMATICS.md for the proof reference and resource contract.
use crate::{exact, interval, math::finite, Error, Result, ScalarInterval};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational as R;
use std::{cmp::Ordering, sync::Arc};

pub const MAX_POLYNOMIAL_DEGREE: usize = 25;

/// A deterministic work limit, independent of geometric tolerance. Exhaustion
/// returns an error and never publishes an incomplete set of roots.
#[derive(Debug, Clone, Copy)]
pub struct RootIsolationOptions {
    pub max_subdivisions: usize,
}
impl Default for RootIsolationOptions {
    fn default() -> Self {
        Self {
            max_subdivisions: 65_536,
        }
    }
}
pub(crate) struct Budget {
    remaining: usize,
}
impl Budget {
    pub(crate) fn new(options: RootIsolationOptions) -> Self {
        Self {
            remaining: options.max_subdivisions,
        }
    }
    fn split(&mut self) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(1)
            .ok_or(Error::ComputationLimit("polynomial root subdivision"))?;
        Ok(())
    }
}

/// Exact polynomial of the represented coefficients, in ascending power order.
/// At most 26 finite coefficients. Empty input denotes the zero polynomial.
#[derive(Debug, Clone)]
pub struct Polynomial {
    pub(crate) exact: IntPolynomial,
}
impl Polynomial {
    pub fn new(coefficients: &[f64]) -> Result<Self> {
        if coefficients.len() > MAX_POLYNOMIAL_DEGREE + 1 {
            return Err(Error::LimitExceeded("polynomial degree"));
        }
        for &x in coefficients {
            finite(x, "polynomial coefficient")?;
        }
        Ok(Self {
            exact: IntPolynomial::new(coefficients.iter().copied().map(exact::integer).collect()),
        })
    }
    /// None for the zero polynomial.
    pub fn degree(&self) -> Option<usize> {
        self.exact.0.len().checked_sub(1)
    }
    pub fn sign(&self, x: f64) -> Result<Ordering> {
        finite(x, "polynomial argument")?;
        Ok(self.exact.sign_at(&rat(x)))
    }
    pub fn real_roots(&self) -> Result<RealRoots> {
        self.real_roots_with_options(RootIsolationOptions::default())
    }
    pub fn real_roots_with_options(&self, options: RootIsolationOptions) -> Result<RealRoots> {
        if self.exact.is_zero() {
            return Ok(RealRoots::All);
        }
        if self.exact.is_constant() {
            return Ok(RealRoots::Finite(vec![]));
        }
        // Cauchy's strict root bound. It need not be representable by f64.
        let leading = abs(self.exact.0.last().unwrap());
        let bound = self.exact.0[..self.exact.0.len() - 1]
            .iter()
            .map(|c| R::new(abs(c), leading.clone()))
            .max()
            .unwrap()
            + one();
        Ok(RealRoots::Finite(isolate(
            &self.exact,
            -&bound,
            bound,
            &mut Budget::new(options),
        )?))
    }
    /// All roots in a finite closed parameter interval, including endpoints.
    pub fn roots_in(&self, lower: f64, upper: f64) -> Result<RealRoots> {
        self.roots_in_with_options(lower, upper, RootIsolationOptions::default())
    }
    pub fn roots_in_with_options(
        &self,
        lower: f64,
        upper: f64,
        options: RootIsolationOptions,
    ) -> Result<RealRoots> {
        finite(lower, "root interval")?;
        finite(upper, "root interval")?;
        if lower > upper {
            return Err(Error::OutOfDomain("reversed root interval"));
        }
        if self.exact.is_zero() {
            return Ok(RealRoots::All);
        }
        Ok(RealRoots::Finite(isolate(
            &self.exact,
            rat(lower),
            rat(upper),
            &mut Budget::new(options),
        )?))
    }
}

#[derive(Debug, Clone)]
pub enum RealRoots {
    /// Every value in the queried domain is a root of the zero polynomial.
    All,
    /// Distinct roots in exact increasing order, each retaining multiplicity.
    Finite(Vec<AlgebraicRoot>),
}

#[derive(Debug)]
struct RootPolynomial {
    polynomial: IntPolynomial,
    sturm: Vec<IntPolynomial>,
}

/// An exact root identity: square-free polynomial and rational isolating
/// interval, or an exact rational witness. Binary64 bounds are only a view and
/// can overlap those of another, distinct root. Values outside f64 still exist.
#[derive(Debug, Clone)]
pub struct AlgebraicRoot {
    defining: Arc<RootPolynomial>,
    lower: R,
    upper: R,
    multiplicity: usize,
}
impl AlgebraicRoot {
    /// Tighten a single-root interval by a bounded amount before repeated sign
    /// queries. Its square-free polynomial has one simple root here, so opposite
    /// endpoint signs certify each bisection without another Sturm chain. This
    /// Verified rational candidates can collapse the interval exactly. This
    /// is only an exact fast-filter aid: undecided signs still use Sturm-Tarski.
    pub(crate) fn refine_for_signs(&mut self, steps: usize) {
        if self.lower == self.upper || self.recognize_rational() {
            return;
        }
        let defining = self.defining.clone();
        let p = &defining.polynomial;
        let left = p.sign_at(&self.lower);
        debug_assert!(left != Ordering::Equal && left != p.sign_at(&self.upper));
        for i in 0..steps {
            let middle = (&self.lower + &self.upper) / R::from_integer(BigInt::from(2));
            let sign = p.sign_at(&middle);
            if sign == Ordering::Equal {
                self.lower = middle.clone();
                self.upper = middle;
                break;
            }
            if sign == left {
                self.lower = middle;
            } else {
                self.upper = middle;
            }
            if (i + 1) % 16 == 0 && self.recognize_rational() {
                break;
            }
        }
    }
    fn recognize_rational(&mut self) -> bool {
        let value = rational_in_interval(&self.lower, &self.upper);
        // Membership and a zero of the defining polynomial certify that this
        // is the unique isolated root. A guess alone never changes the result.
        if self.lower <= value
            && value <= self.upper
            && self.defining.polynomial.sign_at(&value) == Ordering::Equal
        {
            self.lower = value.clone();
            self.upper = value;
            true
        } else {
            false
        }
    }
    pub fn multiplicity(&self) -> usize {
        self.multiplicity
    }
    pub fn bounds(&self) -> Result<ScalarInterval> {
        interval::enclose(|x| self.compare_rational(&rat(x)), "polynomial root")
    }
    pub fn compare(&self, x: f64) -> Result<Ordering> {
        finite(x, "root comparison value")?;
        Ok(self.compare_rational(&rat(x)))
    }
    /// Exact ordering of two real algebraic roots, including roots of different
    /// polynomials. Multiplicity and overlapping binary64 bounds do not affect
    /// equality. No iterative approximate-equality criterion is used.
    pub fn compare_root(&self, other: &Self) -> Ordering {
        if other.lower == other.upper {
            return self.compare_rational(&other.lower);
        }
        if self.compare_rational(&other.lower) != Ordering::Greater {
            return Ordering::Less;
        }
        if self.compare_rational(&other.upper) != Ordering::Less {
            return Ordering::Greater;
        }
        // self is strictly inside other's isolator. Its square-free defining
        // polynomial has exactly one simple root there and changes sign once.
        // A zero proves equality even when the two defining polynomials differ.
        let sign = self.sign_polynomial(&other.defining.polynomial);
        if sign == Ordering::Equal {
            Ordering::Equal
        } else if sign == other.defining.polynomial.sign_at(&other.lower) {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }
    /// Exact sign of another represented polynomial at this root. In
    /// particular, a zero is decided algebraically, without approximate equality.
    pub fn sign_at(&self, polynomial: &Polynomial) -> Ordering {
        self.sign_polynomial(&polynomial.exact)
    }
    pub(crate) fn compare_rational(&self, x: &R) -> Ordering {
        if self.lower == self.upper {
            return self.lower.cmp(x);
        }
        if x <= &self.lower {
            return Ordering::Greater;
        }
        if x >= &self.upper {
            return Ordering::Less;
        }
        if self.defining.polynomial.sign_at(x) == Ordering::Equal {
            return Ordering::Equal;
        }
        if variations(&self.defining.sturm, &self.lower) > variations(&self.defining.sturm, x) {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    }
    pub(crate) fn sign_polynomial(&self, g: &IntPolynomial) -> Ordering {
        if self.lower == self.upper {
            return g.sign_at(&self.lower);
        }
        // At a root of p, g and its positive pseudo-remainder modulo p have
        // the same sign. Reduce before interval evaluation and Sturm-Tarski:
        // repeated contacts can have a low-degree square-free p while their
        // coordinates and implicit derivatives have much higher degree.
        let reduced;
        let g = if g.0.len() >= self.defining.polynomial.0.len() {
            reduced = g.remainder(&self.defining.polynomial);
            &reduced
        } else {
            g
        };
        if g.is_zero() {
            return Ordering::Equal;
        }
        if g.is_constant() {
            return g.0[0].cmp(&BigInt::from(0));
        }
        // Exact interval Horner evaluation is a sufficient fast test only.
        // If it cannot decide, use the complete Sturm-Tarski query.
        let (lo, hi) = g.range_signs(&self.lower, &self.upper);
        if lo == Ordering::Greater {
            return Ordering::Greater;
        }
        if hi == Ordering::Less {
            return Ordering::Less;
        }
        // A wide isolator can make interval dependency obscure an easy sign.
        // Try a bounded exact refinement before building the query's Sturm
        // chain. This cannot decide a false sign: unresolved/zero values still
        // take the complete algebraic path below.
        let mut refined = self.clone();
        refined.refine_for_signs(64);
        if refined.lower == refined.upper {
            return g.sign_at(&refined.lower);
        }
        let (lo, hi) = g.range_signs(&refined.lower, &refined.upper);
        if lo == Ordering::Greater {
            return Ordering::Greater;
        }
        if hi == Ordering::Less {
            return Ordering::Less;
        }
        let p = &self.defining.polynomial;
        let query = sequence(p.clone(), p.derivative().multiply(g));
        let sign = variations(&query, &refined.lower) - variations(&query, &refined.upper);
        debug_assert!((-1..=1).contains(&sign));
        sign.cmp(&0)
    }
}

/// A rational candidate in a closed interval, using common continued-fraction
/// prefixes. If no integer lies in [a,b], both share floor(a) and have positive
/// fractional parts; subtract that integer and reciprocate, reversing bounds.
/// Reversing these exact maps recovers a candidate in the original interval.
/// The iterative form avoids a call stack proportional to rational bit length.
fn rational_in_interval(a: &R, b: &R) -> R {
    let (mut lower, mut upper) = (a.clone(), b.clone());
    let mut prefixes = Vec::new();
    let mut value = loop {
        let integer = lower.ceil();
        if integer <= upper {
            break integer;
        }
        let floor = lower.floor();
        let next_lower = one() / (&upper - &floor);
        let next_upper = one() / (&lower - &floor);
        prefixes.push(floor);
        (lower, upper) = (next_lower, next_upper);
    };
    for prefix in prefixes.into_iter().rev() {
        value = prefix + one() / value;
    }
    value
}

/// Primitive coefficients, ascending powers. Only positive common factors are
/// removed: normalizing the leading sign here would corrupt Sturm sequences.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct IntPolynomial(pub(crate) Vec<BigInt>);
impl IntPolynomial {
    pub(crate) fn new(mut coefficients: Vec<BigInt>) -> Self {
        while coefficients.last().is_some_and(exact::zero) {
            coefficients.pop();
        }
        let content = coefficients
            .iter()
            .fold(BigInt::from(0), |g, x| gcd_integer(g, abs(x)));
        if content > BigInt::from(1) {
            for c in &mut coefficients {
                *c /= &content;
            }
        }
        Self(coefficients)
    }
    pub(crate) fn from_rationals(coefficients: &[R]) -> Self {
        let den = coefficients.iter().fold(BigInt::from(1), |d, x| {
            (&d / gcd_integer(d.clone(), x.denom().clone())) * x.denom()
        });
        Self::new(
            coefficients
                .iter()
                .map(|c| c.numer() * (&den / c.denom()))
                .collect(),
        )
    }
    pub(crate) fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    fn is_constant(&self) -> bool {
        self.0.len() <= 1
    }
    fn positive(mut self) -> Self {
        if self.0.last().is_some_and(|c| c.sign() == Sign::Minus) {
            self = self.negate();
        }
        self
    }
    fn negate(mut self) -> Self {
        for c in &mut self.0 {
            *c = -&*c;
        }
        self
    }
    pub(crate) fn derivative(&self) -> Self {
        Self::new(
            self.0
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, c)| c * BigInt::from(i))
                .collect(),
        )
    }
    fn multiply(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self(vec![]);
        }
        let mut result = vec![BigInt::from(0); self.0.len() + other.0.len() - 1];
        for (i, a) in self.0.iter().enumerate() {
            for (j, b) in other.0.iter().enumerate() {
                result[i + j] += a * b;
            }
        }
        Self::new(result)
    }
    /// A positive multiple of the true rational remainder. Positive scaling at
    /// EVERY elimination step is essential when the divisor's lead is negative.
    fn remainder(&self, divisor: &Self) -> Self {
        debug_assert!(!divisor.is_zero());
        let mut r = self.clone();
        let lead = divisor.0.last().unwrap();
        let magnitude = abs(lead);
        while !r.is_zero() && r.0.len() >= divisor.0.len() {
            let shift = r.0.len() - divisor.0.len();
            let multiplier = if lead.sign() == Sign::Minus {
                -r.0.last().unwrap()
            } else {
                r.0.last().unwrap().clone()
            };
            for c in &mut r.0 {
                *c *= &magnitude;
            }
            for (i, c) in divisor.0.iter().enumerate() {
                r.0[i + shift] -= &multiplier * c;
            }
            r = Self::new(r.0);
        }
        r
    }
    fn gcd(&self, other: &Self) -> Self {
        let (mut a, mut b) = (self.clone(), other.clone());
        while !b.is_zero() {
            let r = a.remainder(&b);
            a = b;
            b = r;
        }
        a.positive()
    }
    fn quotient_exact(&self, divisor: &Self) -> Self {
        debug_assert!(!divisor.is_zero());
        let mut remainder: Vec<R> = self.0.iter().cloned().map(R::from_integer).collect();
        let mut quotient = vec![zero(); self.0.len() - divisor.0.len() + 1];
        let leading = R::from_integer(divisor.0.last().unwrap().clone());
        while !remainder.is_empty() && remainder.len() >= divisor.0.len() {
            let shift = remainder.len() - divisor.0.len();
            let c = remainder.last().unwrap() / &leading;
            quotient[shift] = c.clone();
            for (i, b) in divisor.0.iter().enumerate() {
                remainder[i + shift] -= &c * R::from_integer(b.clone());
            }
            while remainder.last().is_some_and(|x| x == &zero()) {
                remainder.pop();
            }
        }
        debug_assert!(remainder.is_empty());
        Self::from_rationals(&quotient)
    }
    pub(crate) fn sign_at(&self, x: &R) -> Ordering {
        let Some(last) = self.0.last() else {
            return Ordering::Equal;
        };
        let mut value = last.clone();
        let mut den = BigInt::from(1);
        for c in self.0.iter().rev().skip(1) {
            den *= x.denom();
            value = value * x.numer() + c * &den;
        }
        value.cmp(&BigInt::from(0))
    }
    /// Signs of the exact interval-Horner bounds on [a,b]. Keep both bounds
    /// over the same positive denominator instead of reducing four rational
    /// products at every coefficient. Only their signs are needed by the filter.
    fn range_signs(&self, a: &R, b: &R) -> (Ordering, Ordering) {
        let Some((last, rest)) = self.0.split_last() else {
            return (Ordering::Equal, Ordering::Equal);
        };
        let common = a.denom() / gcd_integer(a.denom().clone(), b.denom().clone()) * b.denom();
        let a = a.numer() * (&common / a.denom());
        let b = b.numer() * (&common / b.denom());
        let (mut low, mut high) = (last.clone(), last.clone());
        let mut denominator = BigInt::from(1);
        for c in rest.iter().rev() {
            let products = [&low * &a, &low * &b, &high * &a, &high * &b];
            denominator *= &common;
            let offset = c * &denominator;
            low = products.iter().min().unwrap() + &offset;
            high = products.iter().max().unwrap() + offset;
        }
        (low.cmp(&BigInt::from(0)), high.cmp(&BigInt::from(0)))
    }
}

fn sequence(first: IntPolynomial, second: IntPolynomial) -> Vec<IntPolynomial> {
    let mut chain = vec![first];
    let mut next = second;
    while !next.is_zero() {
        let remainder = chain.last().unwrap().remainder(&next).negate();
        chain.push(next);
        next = remainder;
    }
    chain
}
fn variations(chain: &[IntPolynomial], x: &R) -> i32 {
    let mut previous = Ordering::Equal;
    let mut count = 0;
    for p in chain {
        let sign = p.sign_at(x);
        if sign != Ordering::Equal {
            if previous != Ordering::Equal && sign != previous {
                count += 1;
            }
            previous = sign;
        }
    }
    count
}

pub(crate) fn isolate(
    p: &IntPolynomial,
    lower: R,
    upper: R,
    budget: &mut Budget,
) -> Result<Vec<AlgebraicRoot>> {
    debug_assert!(lower <= upper && !p.is_zero());
    if p.is_constant() {
        return Ok(vec![]);
    }
    let mut common = p.gcd(&p.derivative());
    let square_free = p.quotient_exact(&common).positive();
    let defining = Arc::new(RootPolynomial {
        sturm: sequence(square_free.clone(), square_free.derivative()),
        polynomial: square_free,
    });
    let root = |lower: R, upper: R| AlgebraicRoot {
        defining: defining.clone(),
        lower,
        upper,
        multiplicity: 1,
    };
    let mut roots = Vec::new();
    let q = &defining.polynomial;
    if q.0.len() == 2 {
        let value = R::new(-&q.0[0], q.0[1].clone());
        if value >= lower && value <= upper {
            roots.push(root(value.clone(), value));
        }
    } else if lower == upper {
        if q.sign_at(&lower) == Ordering::Equal {
            roots.push(root(lower.clone(), lower));
        }
    } else {
        if q.sign_at(&lower) == Ordering::Equal {
            roots.push(root(lower.clone(), lower.clone()));
        }
        // Iterative traversal avoids stack overflow for tightly clustered roots.
        enum Entry {
            Interval(R, R, i32, i32),
            Exact(R),
        }
        let chain = &defining.sturm;
        let mut stack = vec![Entry::Interval(
            lower.clone(),
            upper.clone(),
            variations(chain, &lower),
            variations(chain, &upper),
        )];
        while let Some(entry) = stack.pop() {
            let Entry::Interval(a, b, va, vb) = entry else {
                if let Entry::Exact(x) = entry {
                    roots.push(root(x.clone(), x));
                }
                continue;
            };
            let b_zero = q.sign_at(&b) == Ordering::Equal;
            // V(a)-V(b) counts (a,b]; subtract a root at b for this open interval.
            let count = va - vb - i32::from(b_zero);
            if count == 0 {
                continue;
            }
            debug_assert!(count > 0);
            if count == 1 && !b_zero && q.sign_at(&a) != Ordering::Equal {
                roots.push(root(a, b));
                continue;
            }
            budget.split()?;
            let middle = (&a + &b) / R::from_integer(BigInt::from(2));
            let vm = variations(chain, &middle);
            stack.push(Entry::Interval(middle.clone(), b, vm, vb));
            if q.sign_at(&middle) == Ordering::Equal {
                stack.push(Entry::Exact(middle.clone()));
            }
            stack.push(Entry::Interval(a, middle, va, vm));
        }
        if q.sign_at(&upper) == Ordering::Equal {
            roots.push(root(upper.clone(), upper));
        }
    }
    // Repeated gcd layers identify each root's multiplicity. Each square-free
    // layer is a divisor of q, so a root's isolating endpoints cannot be its roots.
    while !common.is_constant() {
        let next = common.gcd(&common.derivative());
        let layer = common.quotient_exact(&next).positive();
        let chain = sequence(layer.clone(), layer.derivative());
        for r in &mut roots {
            let contains = if r.lower == r.upper {
                layer.sign_at(&r.lower) == Ordering::Equal
            } else {
                variations(&chain, &r.lower) > variations(&chain, &r.upper)
            };
            if contains {
                r.multiplicity += 1;
            }
        }
        common = next;
    }
    Ok(roots)
}
fn gcd_integer(mut a: BigInt, mut b: BigInt) -> BigInt {
    while !exact::zero(&b) {
        let r = &a % &b;
        a = b;
        b = r;
    }
    a
}
fn abs(x: &BigInt) -> BigInt {
    if x.sign() == Sign::Minus {
        -x
    } else {
        x.clone()
    }
}
pub(crate) fn rat(x: f64) -> R {
    R::from_float(x).expect("validated finite polynomial input")
}
fn zero() -> R {
    R::from_integer(BigInt::from(0))
}
fn one() -> R {
    R::from_integer(BigInt::from(1))
}
