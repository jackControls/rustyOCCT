//! Exact real-root isolation and sign evaluation by subresultant Sturm sequences.
//!
//! Polynomials are integer vectors; pseudo-division scales by positive
//! leading-coefficient magnitudes so sign variations remain valid. Sturm-Tarski
//! queries decide polynomial signs at an isolated algebraic root, including zero.
//! See MATHEMATICS.md for the proof reference and resource contract.
use crate::{exact, interval, math::finite, Error, Result, ScalarInterval};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational as R;
use std::{cmp::Ordering, sync::Arc};

pub(crate) mod image;

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
/// Extra bisections allowed to reach the rational-root-theorem candidate.
const LEAD_BOUND_STEPS: i64 = 64;
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
    pub(crate) fn rational(value: R) -> Self {
        let polynomial = IntPolynomial::new(vec![-value.numer(), value.denom().clone()]);
        Self {
            defining: Arc::new(RootPolynomial {
                sturm: vec![polynomial.clone(), polynomial.derivative()],
                polynomial,
            }),
            lower: value.clone(),
            upper: value,
            multiplicity: 1,
        }
    }
    pub(crate) fn rational_value(&self) -> Option<&R> {
        (self.lower == self.upper).then_some(&self.lower)
    }
    pub(crate) fn isolator(&self) -> (&R, &R) {
        (&self.lower, &self.upper)
    }
    /// The square-free polynomial with exactly one root in the isolator.
    pub(crate) fn defining(&self) -> &IntPolynomial {
        &self.defining.polynomial
    }
    /// Zero-only query without constructing a signed Sturm-Tarski chain. The
    /// gcd is a square-free divisor of the defining polynomial. Our isolator
    /// contains at most one of its roots, so a sign change is necessary and
    /// sufficient for membership. Its endpoints cannot be roots of the divisor.
    pub(crate) fn vanishes_polynomial(&self, g: &IntPolynomial) -> bool {
        if self.lower == self.upper {
            return g.sign_at(&self.lower) == Ordering::Equal;
        }
        let common = self.defining.polynomial.gcd(g);
        !common.is_constant() && common.sign_at(&self.lower) != common.sign_at(&self.upper)
    }
    /// Tighten a single-root interval by a bounded amount before repeated sign
    /// queries. Its square-free polynomial has one simple root here, so opposite
    /// endpoint signs certify each bisection without another Sturm chain.
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
        // Bisect the ends' numerators over one positive denominator, which
        // doubles each step: the midpoints are the same rationals as
        // `(lower + upper) / 2`, reduced only where they are read (every
        // 16 steps, for the rational root test, and at the end).
        let mut done = 0;
        while done < steps {
            let stop = steps.min((done / 16 + 1) * 16);
            let (ld, ud) = (self.lower.denom(), self.upper.denom());
            let mut den = ld / gcd_integer(ld.clone(), ud.clone()) * ud;
            let mut lo = self.lower.numer() * (&den / ld);
            let mut hi = self.upper.numer() * (&den / ud);
            while done < stop {
                let middle = &lo + &hi;
                den <<= 1;
                let sign = p.sign_at_fraction(&middle, &den);
                if sign == Ordering::Equal {
                    let middle = R::new(middle, den);
                    self.lower = middle.clone();
                    self.upper = middle;
                    return;
                }
                if sign == left {
                    lo = middle;
                    hi <<= 1;
                } else {
                    hi = middle;
                    lo <<= 1;
                }
                done += 1;
            }
            self.lower = lowest(lo, den.clone());
            self.upper = lowest(hi, den);
            if done % 16 == 0 && self.recognize_rational() {
                break;
            }
        }
    }
    /// Bisect until the isolator is narrower than 1/|lead(p)|, where the rational
    /// root theorem leaves one candidate, then test it. Only a short remaining
    /// distance is attempted: repeated long bisections at irrational roots cost
    /// more than the algebraic fallback they would avoid.
    pub(crate) fn refine_to_lead_bound(&mut self) {
        if self.lower == self.upper {
            return;
        }
        let lead = R::from_integer(abs(self.defining.polynomial.0.last().unwrap()));
        let scaled = (&self.upper - &self.lower) * &lead;
        let steps = scaled.numer().bits() as i64 - scaled.denom().bits() as i64 + 2;
        if !(1..=LEAD_BOUND_STEPS).contains(&steps) {
            return;
        }
        self.refine_for_signs(steps as usize);
        if self.lower != self.upper {
            self.recognize_rational();
        }
    }
    fn recognize_rational(&mut self) -> bool {
        // A rational root a/b in lowest terms has b | lead(p), so lead*root is
        // an integer. Once the isolator is narrower than 1/|lead|, that integer
        // is its only candidate. This recognizes huge-denominator roots long
        // before their simplest continued-fraction prefix becomes unique.
        let lead = abs(self.defining.polynomial.0.last().unwrap());
        let (ln, ld) = (self.lower.numer(), self.lower.denom());
        let (un, ud) = (self.upper.numer(), self.upper.denom());
        // (upper-lower)*lead < 1, cleared over the positive denominators.
        let unique = (un * ld - ln * ud) * &lead < ud * ld;
        let mut candidates = vec![rational_in_interval(&self.lower, &self.upper)];
        if unique {
            candidates.push(R::new(div_ceil(&(ln * &lead), ld), lead));
        }
        for value in candidates {
            // Membership and a zero of the defining polynomial certify that this
            // is the unique isolated root. A guess alone never changes the result.
            if le(&self.lower, &value)
                && le(&value, &self.upper)
                && self.defining.polynomial.sign_at(&value) == Ordering::Equal
            {
                self.lower = value.clone();
                self.upper = value;
                return true;
            }
        }
        false
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
    /// Compare two positive affine images without constructing image resultants.
    /// Callers own normalized offsets and strictly positive widths.
    pub(crate) fn compare_affine(
        &self,
        offset: &R,
        width: &R,
        other: &Self,
        other_offset: &R,
        other_width: &R,
    ) -> Ordering {
        debug_assert!(width > &zero() && other_width > &zero());
        if offset == other_offset && width == other_width {
            return self.compare_root(other);
        }
        let delta = offset - other_offset;
        let threshold = |x: &R| (other_width * x - &delta) / width;
        if let Some(x) = other.rational_value() {
            return self.compare_rational(&threshold(x));
        }
        if self.compare_rational(&threshold(&other.lower)) != Ordering::Greater {
            return Ordering::Less;
        }
        if self.compare_rational(&threshold(&other.upper)) != Ordering::Less {
            return Ordering::Greater;
        }
        // The mapped self root lies strictly inside other's unique-root
        // isolator. Positive denominator clearing preserves this sign query.
        let a = delta / other_width;
        let b = width / other_width;
        let mut coefficients = Vec::<R>::new();
        for c in other.defining.polynomial.0.iter().rev() {
            let mut next = vec![zero(); coefficients.len() + 1];
            for (i, v) in coefficients.into_iter().enumerate() {
                next[i] += &a * &v;
                next[i + 1] += &b * &v;
            }
            next[0] += R::from_integer(c.clone());
            coefficients = next;
        }
        let sign = self.sign_polynomial(&IntPolynomial::from_rationals(&coefficients));
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
        let p = &self.defining.polynomial;
        let sign = p.sign_at(x);
        if sign == Ordering::Equal {
            return Ordering::Equal;
        }
        // A square-free polynomial has exactly one simple root in this
        // isolator. An interior point is left of that root iff its sign agrees
        // with the left endpoint; evaluating the entire Sturm chain is needless.
        if sign == p.sign_at(&self.lower) {
            Ordering::Greater
        } else {
            Ordering::Less
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
        // Bisecting the defining polynomial is far cheaper than a gcd with a
        // large query. Reaching width 1/|lead| lets the rational root theorem
        // recognize any rational root; otherwise the tighter isolator can
        // still decide the interval test before the algebraic paths below.
        refined.refine_to_lead_bound();
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
        // A shared square-free factor can prove zero far more cheaply than
        // constructing P'Q. In particular, geometric coordinates may vanish at
        // many roots of a higher-degree stationary equation. Nonzero values
        // still follow the complete signed query below.
        if refined.vanishes_polynomial(g) {
            return Ordering::Equal;
        }
        let p = &self.defining.polynomial;
        let query = sequence(p.clone(), p.derivative().multiply(g));
        let sign = variations(&query, &refined.lower) - variations(&query, &refined.upper);
        debug_assert!((-1..=1).contains(&sign));
        sign.cmp(&0)
    }
}

/// The Mersenne prime `2^61 - 1`.
const PRIME: u64 = (1 << 61) - 1;

fn mul_mod(a: u64, b: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(PRIME)) as u64
}

fn inverse_mod(a: u64) -> u64 {
    // Fermat: a^(p - 2).
    let (mut base, mut e, mut out) = (a, PRIME - 2, 1u64);
    while e > 0 {
        if e & 1 == 1 {
            out = mul_mod(out, base);
        }
        base = mul_mod(base, base);
        e >>= 1;
    }
    out
}

fn reduce_mod(p: &[BigInt]) -> Vec<u64> {
    let m = BigInt::from(PRIME);
    let mut out: Vec<u64> = p
        .iter()
        .map(|c| {
            let r = num_integer::Integer::mod_floor(c, &m);
            r.to_u64_digits().1.first().copied().unwrap_or(0)
        })
        .collect();
    while out.last() == Some(&0) {
        out.pop();
    }
    out
}

/// Whether two nonzero integer polynomials are certainly coprime: their
/// reductions modulo `PRIME`, which divides neither leading coefficient,
/// have a constant gcd (`false`: undecided).
fn coprime_modulo_prime(a: &[BigInt], b: &[BigInt]) -> bool {
    let (mut x, mut y) = (reduce_mod(a), reduce_mod(b));
    if x.len() != a.len() || y.len() != b.len() || x.is_empty() || y.is_empty() {
        return false;
    }
    if x.len() < y.len() {
        std::mem::swap(&mut x, &mut y);
    }
    // Euclid's algorithm in GF(p).
    while y.len() > 1 {
        let inv = inverse_mod(*y.last().expect("nonzero"));
        while x.len() >= y.len() {
            let c = mul_mod(*x.last().expect("nonempty"), inv);
            let shift = x.len() - y.len();
            for (i, t) in y.iter().enumerate() {
                let s = mul_mod(c, *t);
                x[i + shift] = (x[i + shift] + PRIME - s) % PRIME;
            }
            while x.last() == Some(&0) {
                x.pop();
            }
        }
        std::mem::swap(&mut x, &mut y);
    }
    // A nonzero constant divides both: coprime; zero: `x` the gcd.
    y.len() == 1 || x.len() == 1
}

/// `n / d` (`d > 0`) in lowest terms, as `R::new` gives it, by the crate's
/// Lehmer gcd (`num_rational` reduces by Stein's binary one, quadratic in
/// the operands' length: an isolator's ends after hundreds of bisections).
fn lowest(n: BigInt, d: BigInt) -> R {
    let g = crate::rational::gcd(&n, &d);
    if g == BigInt::from(1) {
        R::new_raw(n, d)
    } else {
        R::new_raw(n / &g, d / &g)
    }
}

/// A rational candidate in a closed interval, using common continued-fraction
/// prefixes. If no integer lies in [a,b], both share floor(a) and have positive
/// fractional parts; subtract that integer and reciprocate, reversing bounds.
/// Reversing these exact maps recovers a candidate in the original interval.
/// Integer numerator/denominator pairs avoid a rational reduction per step,
/// and continued-fraction convergents are already in lowest terms.
pub(crate) fn rational_in_interval(a: &R, b: &R) -> R {
    // lower = ln/ld and upper = un/ud with positive denominators.
    let (mut ln, mut ld) = (a.numer().clone(), a.denom().clone());
    let (mut un, mut ud) = (b.numer().clone(), b.denom().clone());
    let mut prefixes = Vec::new();
    let integer = loop {
        let integer = div_ceil(&ln, &ld);
        if &integer * &ud <= un {
            break integer;
        }
        // lower is not an integer here, so both remainders are positive.
        let floor = &integer - 1;
        let (rl, ru) = (&ln - &floor * &ld, &un - &floor * &ud);
        prefixes.push(floor);
        (ln, ld, un, ud) = (ud, ru, ld, rl);
    };
    let (mut numer, mut denom) = (integer, BigInt::from(1));
    for prefix in prefixes.into_iter().rev() {
        (numer, denom) = (&prefix * &numer + &denom, numer);
    }
    R::new_raw(numer, denom)
}
/// a <= b by one cross-multiplication. num-rational's Ord walks floor-division
/// continued fractions, which is far slower for multi-thousand-bit operands.
fn le(a: &R, b: &R) -> bool {
    a.numer() * b.denom() <= b.numer() * a.denom()
}
fn div_ceil(n: &BigInt, d: &BigInt) -> BigInt {
    num_integer::Integer::div_ceil(n, d)
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
        let content = content(&coefficients, BigInt::from(0));
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
    pub(crate) fn is_constant(&self) -> bool {
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
        Self::new(self.pseudo_remainder(divisor).0)
    }
    /// |lead(divisor)|^(deg difference + 1) times the rational remainder,
    /// without content removal.
    fn pseudo_remainder(&self, divisor: &Self) -> Self {
        debug_assert!(!divisor.is_zero());
        let mut r = self.clone();
        let lead = divisor.0.last().unwrap();
        let magnitude = abs(lead);
        // The classical pseudo-remainder scales by exactly |lead|^(delta+1),
        // even when cancellation skips elimination steps; subresultant
        // divisibility depends on that exact power.
        let mut missing = (self.0.len() + 1).saturating_sub(divisor.0.len());
        while !r.is_zero() && r.0.len() >= divisor.0.len() {
            missing -= 1;
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
            // Content removal only bounds growth; one pass at the end costs
            // far less than a multiprecision gcd sweep after every step.
            while r.0.last().is_some_and(exact::zero) {
                r.0.pop();
            }
        }
        if missing > 0 && !r.is_zero() {
            let scale = magnitude.pow(missing as u32);
            for c in &mut r.0 {
                *c *= &scale;
            }
        }
        r
    }
    fn divide_exact(mut self, divisor: &BigInt) -> Self {
        for c in &mut self.0 {
            debug_assert!(exact::zero(&(&*c % divisor)));
            *c /= divisor;
        }
        self
    }
    /// Whether the two (nonzero) are certainly coprime, by their reductions
    /// modulo a prime (`false`: undecided).
    pub(crate) fn coprime_with(&self, other: &Self) -> bool {
        !self.is_zero() && !other.is_zero() && coprime_modulo_prime(&self.0, &other.0)
    }
    pub(crate) fn gcd(&self, other: &Self) -> Self {
        let (a, b) = if self.0.len() >= other.0.len() {
            (self, other)
        } else {
            (other, self)
        };
        if b.is_zero() {
            return Self::new(a.0.clone()).positive();
        }
        // Coprime modulo a prime dividing neither leading coefficient:
        // coprime over the rationals (the gcd's reduction divides both
        // reductions and keeps its degree), the subresultant chain's
        // constant's one form, `1`.
        if coprime_modulo_prime(&a.0, &b.0) {
            return Self(vec![BigInt::from(1)]);
        }
        let mut chain = subresultants(Self::new(a.0.clone()), Self::new(b.0.clone()), false);
        Self::new(chain.pop().unwrap().0).positive()
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
        self.sign_at_fraction(x.numer(), x.denom())
    }
    /// The sign at `numerator / denominator`, `denominator > 0`, in any
    /// terms: the homogenized value is the polynomial's times a positive
    /// power of the denominator.
    fn sign_at_fraction(&self, numerator: &BigInt, denominator: &BigInt) -> Ordering {
        debug_assert!(denominator.sign() == Sign::Plus);
        let Some(last) = self.0.last() else {
            return Ordering::Equal;
        };
        if let Some(sign) = self.sign_filter(numerator, denominator) {
            return sign;
        }
        let mut value = last.clone();
        let mut den = BigInt::from(1);
        for c in self.0.iter().rev().skip(1) {
            den *= denominator;
            value = value * numerator + c * &den;
        }
        value.cmp(&BigInt::from(0))
    }
    /// The sign at `numerator / denominator` where a floating-point
    /// evaluation certifies it, `None` otherwise (the exact Horner decides).
    /// The homogenized Horner runs in binary64 mantissas with unbounded
    /// exponents (`Xf`): each coefficient and both ends rounded once (with
    /// the digits below their top 128 bits dropped, relative error below
    /// `2u`), each product and sum once more, so every term of the value is
    /// perturbed by at most `5 d + 5` factors `1 + u` (`d` the degree) and
    /// the computed value is within `gamma_{5d+5}` of the terms' absolute
    /// sum, which the same evaluation on absolute values bounds from below
    /// within the same factor (Higham, Horner's rule). A value beyond `(6 d
    /// + 8) 2u` times that sum has the exact value's sign; a smaller one,
    /// zero among them, is left undecided.
    fn sign_filter(&self, numerator: &BigInt, denominator: &BigInt) -> Option<Ordering> {
        let (xn, xd) = (Xf::of(numerator), Xf::of(denominator));
        let an = xn.abs();
        let mut coefficients = self.0.iter().rev();
        let top = Xf::of(coefficients.next()?);
        let (mut value, mut sum, mut power) = (top, top.abs(), Xf::ONE);
        for c in coefficients {
            power = power.mul(xd);
            let term = Xf::of(c).mul(power);
            value = value.mul(xn).add(term);
            sum = sum.mul(an).add(term.abs());
        }
        let k = (6 * self.0.len() + 8) as f64;
        let bound = sum.mul(Xf::new(k * f64::EPSILON, 0));
        if value.m == 0.0 || !value.abs().exceeds(bound) {
            return None;
        }
        Some(if value.m > 0.0 {
            Ordering::Greater
        } else {
            Ordering::Less
        })
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

/// A binary64 mantissa in `[0.5, 1)` (or zero) times `2^e`, `e` unbounded:
/// big integers' leading bits for `IntPolynomial::sign_filter`. A product
/// or a sum rounds once (a sum's smaller operand scaled exactly by a power
/// of two, or dropped where it lies below `2^-899` of the larger).
#[derive(Debug, Clone, Copy)]
struct Xf {
    m: f64,
    e: i64,
}

impl Xf {
    const ONE: Self = Self { m: 0.5, e: 1 };

    fn new(m: f64, e: i64) -> Self {
        if m == 0.0 {
            return Self { m: 0.0, e: 0 };
        }
        let bits = m.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i64;
        debug_assert!(exponent != 0 && exponent != 0x7ff, "a normal mantissa");
        Self {
            m: f64::from_bits((bits & !(0x7ff << 52)) | (1022 << 52)),
            e: e + exponent - 1022,
        }
    }

    /// The integer's top 128 bits rounded to nearest.
    fn of(x: &BigInt) -> Self {
        let mut digits = x.magnitude().iter_u64_digits();
        let n = digits.len() as i64;
        let Some(top) = digits.next_back() else {
            return Self::new(0.0, 0);
        };
        let next = digits.next_back().unwrap_or(0);
        let m = ((u128::from(top) << 64) | u128::from(next)) as f64;
        let e = 64 * (n - 2);
        Self::new(if x.sign() == Sign::Minus { -m } else { m }, e)
    }

    fn abs(self) -> Self {
        Self {
            m: self.m.abs(),
            e: self.e,
        }
    }

    fn mul(self, other: Self) -> Self {
        Self::new(self.m * other.m, self.e + other.e)
    }

    fn add(self, other: Self) -> Self {
        let (big, small) = if self.m == 0.0 || (other.m != 0.0 && other.e > self.e) {
            (other, self)
        } else {
            (self, other)
        };
        if small.m == 0.0 || big.e - small.e > 900 {
            return big;
        }
        let scaled = small.m * f64::powi(2.0, (small.e - big.e) as i32);
        Self::new(big.m + scaled, big.e)
    }

    /// Whether a nonnegative value exceeds another.
    fn exceeds(self, other: Self) -> bool {
        if other.m == 0.0 {
            return self.m > 0.0;
        }
        self.m > 0.0 && (self.e > other.e || (self.e == other.e && self.m > other.m))
    }
}

fn sequence(first: IntPolynomial, second: IntPolynomial) -> Vec<IntPolynomial> {
    subresultants(first, second, true)
}
/// Magnitude subresultant remainder sequence (Collins; Cohen, Algorithm 3.3.1).
/// Each element is a positive multiple of the next Euclidean remainder, negated
/// for a Sturm chain. Dividing by |g| h^delta is exact because these terms
/// differ from the classical subresultants only by sign, so no integer content
/// is ever computed. The last element is a nonzero multiple of the gcd.
fn subresultants(first: IntPolynomial, second: IntPolynomial, sturm: bool) -> Vec<IntPolynomial> {
    let mut chain = vec![first];
    if second.is_zero() {
        return chain;
    }
    chain.push(second);
    let (mut g, mut h) = (BigInt::from(1), BigInt::from(1));
    loop {
        let (a, b) = (&chain[chain.len() - 2], &chain[chain.len() - 1]);
        if b.is_constant() {
            break;
        }
        if a.0.len() < b.0.len() {
            // A Sturm-Tarski query can start with a higher-degree second term.
            // Its remainder is the first term itself; restart the scaling there.
            let r = a.clone();
            (g, h) = (BigInt::from(1), BigInt::from(1));
            chain.push(if sturm { r.negate() } else { r });
            continue;
        }
        let delta = (a.0.len() - b.0.len()) as u32;
        let r = a.pseudo_remainder(b);
        if r.is_zero() {
            break;
        }
        let r = r.divide_exact(&(&g * h.pow(delta)));
        let lead = abs(b.0.last().unwrap());
        h = if delta == 0 {
            h
        } else {
            lead.pow(delta) / h.pow(delta - 1)
        };
        g = lead;
        chain.push(if sturm { r.negate() } else { r });
    }
    chain
}
fn variations(chain: &[IntPolynomial], x: &R) -> i32 {
    variations_at(chain, x.numer(), x.denom())
}
/// Sign variations at `numerator / denominator`, `denominator > 0`, in any
/// terms.
fn variations_at(chain: &[IntPolynomial], numerator: &BigInt, denominator: &BigInt) -> i32 {
    let mut previous = Ordering::Equal;
    let mut count = 0;
    for p in chain {
        let sign = p.sign_at_fraction(numerator, denominator);
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
    isolate_with_gcd(p, p.gcd(&p.derivative()), lower, upper, budget)
}

/// `isolate` with `gcd(p, p')` (as `IntPolynomial::gcd` gives it) already
/// computed by the caller.
pub(crate) fn isolate_with_gcd(
    p: &IntPolynomial,
    mut common: IntPolynomial,
    lower: R,
    upper: R,
    budget: &mut Budget,
) -> Result<Vec<AlgebraicRoot>> {
    debug_assert!(lower <= upper && !p.is_zero());
    if p.is_constant() {
        return Ok(vec![]);
    }
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
        // An interval's ends are numerators over the ends' common
        // denominator times `2^level`: a midpoint is their sum one level
        // down, the same rational as `(a + b) / 2`, reduced only where a
        // root is published (no gcd per bisection).
        enum Entry {
            Interval(BigInt, BigInt, u32, i32, i32),
            Exact(BigInt, u32),
        }
        let chain = &defining.sturm;
        let (ld, ud) = (lower.denom(), upper.denom());
        let base = ld / gcd_integer(ld.clone(), ud.clone()) * ud;
        let den = |level: u32| &base << level as usize;
        let at = |n: BigInt, level: u32| R::new(n, den(level));
        let mut stack = vec![Entry::Interval(
            lower.numer() * (&base / ld),
            upper.numer() * (&base / ud),
            0,
            variations(chain, &lower),
            variations(chain, &upper),
        )];
        while let Some(entry) = stack.pop() {
            let (a, b, level, va, vb) = match entry {
                Entry::Exact(x, level) => {
                    let x = at(x, level);
                    roots.push(root(x.clone(), x));
                    continue;
                }
                Entry::Interval(a, b, level, va, vb) => (a, b, level, va, vb),
            };
            let d = den(level);
            let b_zero = q.sign_at_fraction(&b, &d) == Ordering::Equal;
            // V(a)-V(b) counts (a,b]; subtract a root at b for this open interval.
            let count = va - vb - i32::from(b_zero);
            if count == 0 {
                continue;
            }
            debug_assert!(count > 0);
            if count == 1 && !b_zero && q.sign_at_fraction(&a, &d) != Ordering::Equal {
                roots.push(root(at(a, level), at(b, level)));
                continue;
            }
            budget.split()?;
            let middle = &a + &b;
            let (a, b, level) = (a << 1usize, b << 1usize, level + 1);
            let d = den(level);
            let vm = variations_at(chain, &middle, &d);
            let exact = q.sign_at_fraction(&middle, &d) == Ordering::Equal;
            stack.push(Entry::Interval(middle.clone(), b, level, vm, vb));
            if exact {
                stack.push(Entry::Exact(middle.clone(), level));
            }
            stack.push(Entry::Interval(a, middle, level, va, vm));
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
/// Nonnegative gcd of nonnegative integers. Stein's binary algorithm avoids a
/// multiprecision division per Euclidean step, but it removes only about one
/// bit per step from the larger operand. One initial division balances the
/// operands, and a unit operand needs no work at all.
fn gcd_integer(a: BigInt, b: BigInt) -> BigInt {
    let (big, small) = if a >= b { (a, b) } else { (b, a) };
    if exact::zero(&small) {
        return big;
    }
    if small == BigInt::from(1) {
        return small;
    }
    let rest = &big % &small;
    num_integer::Integer::gcd(&small, &rest)
}
/// Integer content, stopping as soon as it is one.
fn content(v: &[BigInt], initial: BigInt) -> BigInt {
    let mut g = initial;
    for x in v {
        if g == BigInt::from(1) {
            break;
        }
        g = gcd_integer(g, abs(x));
    }
    g
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The modular shortcut answers only coprime pairs, as the subresultant
    /// The floating-point sign filter answers only where the exact Horner
    /// agrees: coefficients of a few to a few thousand bits, points far out,
    /// near roots (a product of known linear factors, at their roots and an
    /// ulp-sized step beside them) and with cancelling terms.
    #[test]
    fn sign_filter_agrees_with_the_exact_horner() {
        let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let big = |bits: u32, next: &mut dyn FnMut() -> u64| -> BigInt {
            let mut x = BigInt::from(0);
            for _ in 0..bits.div_ceil(64) {
                x = (x << 64) + BigInt::from(next());
            }
            x >>= (64 - bits % 64) % 64;
            if next() % 2 == 0 {
                -x
            } else {
                x
            }
        };
        let (mut decided, mut total) = (0, 0);
        for round in 0..400u32 {
            let bits = [3, 40, 70, 200, 1500, 3000][round as usize % 6];
            let degree = 1 + (round as usize % 9);
            // Half the polynomials vanish at known rationals.
            let roots: Vec<(BigInt, BigInt)> = (0..degree)
                .map(|_| {
                    let d = big(bits / 3 + 2, &mut next).magnitude().clone() + 1u32;
                    (big(bits / 3 + 4, &mut next), BigInt::from(d))
                })
                .collect();
            let p = if round % 2 == 0 {
                let mut c = vec![BigInt::from(1)];
                for (n, d) in &roots {
                    // times (d t - n)
                    let mut out = vec![BigInt::from(0); c.len() + 1];
                    for (i, x) in c.iter().enumerate() {
                        out[i + 1] += x * d;
                        out[i] -= x * n;
                    }
                    c = out;
                }
                IntPolynomial::new(c)
            } else {
                IntPolynomial::new((0..=degree).map(|_| big(bits, &mut next)).collect())
            };
            if p.is_zero() {
                continue;
            }
            let mut points: Vec<(BigInt, BigInt)> = roots.clone();
            for (n, d) in &roots {
                let k = BigInt::from(1) << 300;
                points.push((n * &k + 1, d * &k));
                points.push((n * &k - 1, d * &k));
            }
            for _ in 0..4 {
                let d = big(bits.min(400), &mut next).magnitude().clone() + 1u32;
                points.push((big(bits.min(400) + 30, &mut next), BigInt::from(d)));
            }
            for (n, d) in &points {
                let exact = {
                    let mut value = p.0.last().unwrap().clone();
                    let mut den = BigInt::from(1);
                    for c in p.0.iter().rev().skip(1) {
                        den *= d;
                        value = value * n + c * &den;
                    }
                    value.cmp(&BigInt::from(0))
                };
                total += 1;
                if let Some(sign) = p.sign_filter(n, d) {
                    decided += 1;
                    assert_eq!(sign, exact, "{:?} at {n}/{d}", p.0);
                }
                assert_eq!(p.sign_at_fraction(n, d), exact);
            }
        }
        // Most points are decided by the filter; the roots never are.
        assert!(decided * 2 > total, "{decided} of {total}");
    }

    /// The modular shortcut answers only coprime pairs, as the subresultant
    /// chain does; a common factor, or a prime dividing a leading
    /// coefficient, goes the exact way.
    #[test]
    fn modular_coprimality_agrees_with_subresultants() {
        let p = |c: &[i64]| IntPolynomial::new(c.iter().map(|&x| BigInt::from(x)).collect());
        let chain_gcd = |a: &IntPolynomial, b: &IntPolynomial| {
            let mut chain = subresultants(a.clone(), b.clone(), false);
            IntPolynomial::new(chain.pop().unwrap().0).positive()
        };
        let f = p(&[-7, 3, 0, 5, -2, 11]);
        let g = p(&[4, 0, -9, 1]);
        let h = p(&[2, -3, 1]);
        assert!(coprime_modulo_prime(&f.0, &g.0));
        assert_eq!(f.gcd(&g), chain_gcd(&f, &g));
        assert_eq!(f.gcd(&g), p(&[1]));
        let (fh, gh) = (f.multiply(&h), g.multiply(&h));
        assert!(!coprime_modulo_prime(&fh.0, &gh.0));
        assert_eq!(fh.gcd(&gh), h);
        // A leading coefficient the prime divides: undecided.
        let big = BigInt::from(PRIME) * 3;
        let q = IntPolynomial(vec![BigInt::from(1), BigInt::from(2), big]);
        assert!(!coprime_modulo_prime(&q.0, &g.0));
        assert_eq!(q.gcd(&g), chain_gcd(&q, &g));
        // A square-free test of a polynomial with its derivative.
        let sq = h.multiply(&h).multiply(&g);
        assert!(!coprime_modulo_prime(&sq.0, &sq.derivative().0));
        assert_eq!(sq.gcd(&sq.derivative()), h);
    }

    /// The previous rational-arithmetic implementation, kept as a reference.
    fn reference(a: &R, b: &R) -> R {
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

    /// The previous primitive Euclidean Sturm chain, kept as a reference.
    fn primitive_chain(first: IntPolynomial, second: IntPolynomial) -> Vec<IntPolynomial> {
        let mut chain = vec![first];
        let mut next = second;
        while !next.is_zero() {
            let remainder = chain.last().unwrap().remainder(&next).negate();
            chain.push(next);
            next = remainder;
        }
        chain
    }

    #[test]
    fn subresultant_chains_are_positive_multiples_of_primitive_chains() {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let poly = |c: Vec<i64>| IntPolynomial::new(c.into_iter().map(BigInt::from).collect());
        let mut pairs = vec![
            // Degree drops by more than one; the second term has higher degree.
            (poly(vec![1, 0, 0, 0, 1]), poly(vec![0, 0, 1, 0, 0, 1])),
            (poly(vec![-2, 0, 1]), poly(vec![0, 2])),
            (
                poly(vec![1, 0, 0, 0, 0, 0, -1]),
                poly(vec![0, 0, 0, 0, 0, 6]),
            ),
        ];
        for _ in 0..400 {
            let random = |n: usize, next: &mut dyn FnMut() -> u64| {
                poly((0..n).map(|_| (next() % 41) as i64 - 20).collect())
            };
            let a = random(2 + (next() % 12) as usize, &mut next);
            let b = random(1 + (next() % 12) as usize, &mut next);
            if !a.is_zero() && !b.is_zero() {
                pairs.push((a, b));
            }
        }
        for (a, b) in pairs {
            let fast = subresultants(a.clone(), b.clone(), true);
            let slow = primitive_chain(a.clone(), b.clone());
            assert_eq!(fast.len(), slow.len(), "{a:?} {b:?}");
            for (x, y) in fast.iter().zip(&slow) {
                // Same primitive part and the same sign: a positive multiple.
                let px = IntPolynomial::new(x.0.clone());
                assert_eq!(px, *y, "{a:?} {b:?}");
            }
            let g = a.gcd(&b);
            assert_eq!(
                g,
                primitive_chain(a.clone(), b.clone())
                    .pop()
                    .unwrap()
                    .positive()
            );
        }
    }

    #[test]
    fn integer_continued_fraction_matches_rational_reference() {
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let big = |bits: u32, x: u64| BigInt::from(x) << bits as usize;
        let mut cases = vec![
            (R::from_integer(3.into()), R::from_integer(3.into())),
            (R::new((-7).into(), 3.into()), R::new((-2).into(), 1.into())),
            (R::new(1.into(), 3.into()), R::new(1.into(), 3.into())),
            (
                R::new(big(0, 1) + big(1012, 3), big(1012, 3) * 3),
                R::new(big(0, 1) + big(1012, 3) + 1, big(1012, 3) * 3),
            ),
        ];
        for _ in 0..2000 {
            let scale = next() % 700;
            let d1 = big(scale as u32, 1 + next() % 1000);
            let d2 = big((next() % 700) as u32, 1 + next() % 1000);
            let n1 = BigInt::from(next() as i64 >> (next() % 60)) * &d1 / 1000;
            let a = R::new(n1, d1);
            let b = &a + R::new(BigInt::from(1 + next() % 5000), d2);
            cases.push((a, b));
        }
        for (a, b) in cases {
            let actual = rational_in_interval(&a, &b);
            assert_eq!(actual, reference(&a, &b), "{a} {b}");
            assert!(a <= actual && actual <= b);
            // new_raw requires lowest terms and a positive denominator.
            assert_eq!(
                actual,
                R::new(actual.numer().clone(), actual.denom().clone())
            );
        }
    }
}
