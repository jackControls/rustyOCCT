//! Certified real enclosures for geometric validation.
//!
//! An [`Interval`] has rational endpoints and always contains the real value it
//! represents. Endpoints are rounded outward to a fixed dyadic grid after
//! multiplication so sizes stay bounded. `pi` comes from Machin's formula with
//! alternating-series remainder bounds; `cos_sin` reduces by multiples of pi/2
//! and sums Taylor series with explicit remainder terms; `atan` halves its
//! argument before an alternating series. Nothing here uses binary64 libm.
use num_bigint::{BigInt, Sign};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::OnceLock;

/// Outward rounding grid: 2^-GRID absolute (about 1.6e-58).
const GRID: usize = 192;
/// Series terms: |z| <= ~0.2 for atan and |r| <= ~0.8 for sin/cos leave
/// remainders far below the rounding grid.
const ATAN_TERMS: usize = 40;
const TAYLOR_TERMS: usize = 24;
/// Integer square roots are taken on a 2^-SQRT_BITS grid.
const SQRT_BITS: usize = 200;

fn zero() -> R {
    R::from_integer(BigInt::from(0))
}
fn one() -> R {
    R::from_integer(BigInt::from(1))
}
fn int(n: i64) -> R {
    R::from_integer(BigInt::from(n))
}
/// Sign of a normalized rational (its denominator is positive).
fn sgn(x: &R) -> Sign {
    x.numer().sign()
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Interval {
    lo: R,
    hi: R,
}

fn floor_grid(x: &R) -> R {
    if x.denom().bits() as usize <= GRID + 1 {
        return x.clone();
    }
    let scale = BigInt::from(1) << GRID;
    R::new(
        (x * R::from_integer(scale.clone())).floor().to_integer(),
        scale,
    )
}
fn ceil_grid(x: &R) -> R {
    if x.denom().bits() as usize <= GRID + 1 {
        return x.clone();
    }
    let scale = BigInt::from(1) << GRID;
    R::new(
        (x * R::from_integer(scale.clone())).ceil().to_integer(),
        scale,
    )
}

impl Interval {
    pub(crate) fn exact(x: R) -> Self {
        Self {
            lo: x.clone(),
            hi: x,
        }
    }
    pub(crate) fn from_f64(x: f64) -> Self {
        Self::exact(R::from_float(x).expect("finite certified input"))
    }
    pub(crate) fn new(lo: R, hi: R) -> Self {
        debug_assert!(lo <= hi);
        Self { lo, hi }
    }
    pub(crate) fn lo(&self) -> &R {
        &self.lo
    }
    pub(crate) fn hi(&self) -> &R {
        &self.hi
    }
    fn rounded(lo: R, hi: R) -> Self {
        Self {
            lo: floor_grid(&lo),
            hi: ceil_grid(&hi),
        }
    }
    pub(crate) fn add(&self, o: &Self) -> Self {
        Self::rounded(&self.lo + &o.lo, &self.hi + &o.hi)
    }
    pub(crate) fn sub(&self, o: &Self) -> Self {
        Self::rounded(&self.lo - &o.hi, &self.hi - &o.lo)
    }
    pub(crate) fn neg(&self) -> Self {
        Self::new(-&self.hi, -&self.lo)
    }
    pub(crate) fn mul(&self, o: &Self) -> Self {
        let p = [
            &self.lo * &o.lo,
            &self.lo * &o.hi,
            &self.hi * &o.lo,
            &self.hi * &o.hi,
        ];
        let lo = p.iter().min().unwrap().clone();
        let hi = p.iter().max().unwrap().clone();
        Self::rounded(lo, hi)
    }
    pub(crate) fn scale(&self, k: &R) -> Self {
        self.mul(&Self::exact(k.clone()))
    }
    pub(crate) fn square(&self) -> Self {
        let (a, b) = (self.abs_lo(), self.abs_hi());
        Self::rounded(&a * &a, &b * &b)
    }
    /// Smallest absolute value in the interval.
    pub(crate) fn abs_lo(&self) -> R {
        if sgn(&self.lo) != Sign::Minus {
            self.lo.clone()
        } else if sgn(&self.hi) != Sign::Plus {
            -&self.hi
        } else {
            zero()
        }
    }
    pub(crate) fn abs_hi(&self) -> R {
        (-&self.lo).max(self.hi.clone())
    }
    /// Certified sign, or None when the interval contains zero in its interior
    /// or on one side of a nonzero width.
    pub(crate) fn sign(&self) -> Option<Ordering> {
        if sgn(&self.lo) == Sign::Plus {
            Some(Ordering::Greater)
        } else if sgn(&self.hi) == Sign::Minus {
            Some(Ordering::Less)
        } else if sgn(&self.lo) == Sign::NoSign && sgn(&self.hi) == Sign::NoSign {
            Some(Ordering::Equal)
        } else {
            None
        }
    }
    /// Certified comparison with an exact value.
    pub(crate) fn cmp_exact(&self, x: &R) -> Option<Ordering> {
        if &self.lo > x {
            Some(Ordering::Greater)
        } else if &self.hi < x {
            Some(Ordering::Less)
        } else if &self.lo == x && &self.hi == x {
            Some(Ordering::Equal)
        } else {
            None
        }
    }
    /// Enclosure of sqrt over a nonnegative interval (negative parts clamp to zero).
    pub(crate) fn sqrt(&self) -> Self {
        let lo = self.lo.clone().max(zero());
        let hi = self.hi.clone().max(zero());
        Self::new(sqrt_bound(&lo, false), sqrt_bound(&hi, true))
    }
    /// Division by an interval certainly excluding zero.
    pub(crate) fn div(&self, o: &Self) -> Option<Self> {
        o.sign().filter(|s| *s != Ordering::Equal)?;
        let inv = Self::rounded(one() / &o.hi, one() / &o.lo);
        Some(self.mul(&inv))
    }
}

/// sqrt(x) rounded down (up=false) or up on the 2^-SQRT_BITS grid, from integer isqrt.
fn sqrt_bound(x: &R, up: bool) -> R {
    const BITS: usize = SQRT_BITS;
    let scale = BigInt::from(1) << (2 * BITS);
    let scaled = x * R::from_integer(scale);
    let n = if up { scaled.ceil() } else { scaled.floor() }.to_integer();
    let mut root = n.sqrt();
    if up && &root * &root < n {
        root += 1;
    }
    R::new(root, BigInt::from(1) << BITS)
}

/// sum_{k<ATAN_TERMS} (-1)^k z^(2k+1)/(2k+1) for |z| <= 1/4, plus the
/// alternating-series remainder bound.
fn atan_series(z: &Interval) -> Interval {
    let mut total = Interval::exact(zero());
    let z2 = z.square();
    let mut power = z.clone();
    for k in 0..ATAN_TERMS {
        let term = power.scale(&R::new(1.into(), BigInt::from(2 * k + 1)));
        total = if k % 2 == 0 {
            total.add(&term)
        } else {
            total.sub(&term)
        };
        power = power.mul(&z2);
    }
    let bound = power.abs_hi() / int((2 * ATAN_TERMS + 1) as i64);
    total.add(&Interval::new(-&bound, bound))
}

pub(crate) fn pi() -> &'static Interval {
    static PI: OnceLock<Interval> = OnceLock::new();
    PI.get_or_init(|| {
        let a = atan_series(&Interval::exact(R::new(1.into(), 5.into())));
        let b = atan_series(&Interval::exact(R::new(1.into(), 239.into())));
        a.scale(&int(16)).sub(&b.scale(&int(4)))
    })
}

/// atan of an exact rational, certified.
fn atan_exact(x: &R) -> Interval {
    if sgn(x) == Sign::Minus {
        return atan_exact(&-x).neg();
    }
    if x > &one() {
        // atan(x) = pi/2 - atan(1/x)
        return pi()
            .scale(&R::new(1.into(), 2.into()))
            .sub(&atan_exact(&(one() / x)));
    }
    // Two halvings: atan(x) = 2 atan(x / (1 + sqrt(1 + x^2))), interval-wise.
    let mut z = Interval::exact(x.clone());
    for _ in 0..2 {
        let d = Interval::exact(one())
            .add(&z.square())
            .sqrt()
            .add(&Interval::exact(one()));
        z = z.div(&d).expect("positive denominator");
    }
    atan_series(&z).scale(&int(4))
}

/// atan over an interval (monotone).
pub(crate) fn atan(z: &Interval) -> Interval {
    Interval::new(atan_exact(z.lo()).lo.clone(), atan_exact(z.hi()).hi.clone())
}

/// Angle of (x, y) in (-pi, pi], certified; None near the branch cut or origin.
pub(crate) fn atan2(y: &Interval, x: &Interval) -> Option<Interval> {
    let half = pi().scale(&R::new(1.into(), 2.into()));
    if x.sign() == Some(Ordering::Greater) {
        return Some(atan(&y.div(x)?));
    }
    match y.sign() {
        Some(Ordering::Greater) => Some(half.sub(&atan(&x.div(y)?))),
        Some(Ordering::Less) => Some(half.neg().sub(&atan(&x.div(y)?))),
        _ => None,
    }
}

/// (cos x, sin x) for an exact rational x, certified.
pub(crate) fn cos_sin(x: &R) -> (Interval, Interval) {
    if sgn(x) == Sign::NoSign {
        return (Interval::exact(one()), Interval::exact(zero()));
    }
    // k = nearest multiple of pi/2, from a rational midpoint of pi.
    let half_pi = pi().scale(&R::new(1.into(), 2.into()));
    let mid = (half_pi.lo() + half_pi.hi()) / int(2);
    let k = (x / &mid).round().to_integer();
    let r = Interval::exact(x.clone()).sub(&half_pi.scale(&R::from_integer(k.clone())));
    let (c, s) = taylor_cos_sin(&r);
    let quadrant = ((k % 4) + 4) % 4;
    if quadrant == BigInt::from(0) {
        (c, s)
    } else if quadrant == BigInt::from(1) {
        (s.neg(), c)
    } else if quadrant == BigInt::from(2) {
        (c.neg(), s.neg())
    } else {
        (s, c.neg())
    }
}

/// Taylor series for |r| <= ~0.8 with Lagrange remainder bounds.
fn taylor_cos_sin(r: &Interval) -> (Interval, Interval) {
    let r2 = r.square();
    let (mut cos, mut sin) = (Interval::exact(zero()), Interval::exact(zero()));
    let mut term_c = Interval::exact(one()); // r^(2k)/(2k)!
    let mut term_s = r.clone(); // r^(2k+1)/(2k+1)!
    for k in 0..TAYLOR_TERMS {
        if k % 2 == 0 {
            cos = cos.add(&term_c);
            sin = sin.add(&term_s);
        } else {
            cos = cos.sub(&term_c);
            sin = sin.sub(&term_s);
        }
        let a = R::new(1.into(), BigInt::from((2 * k + 1) * (2 * k + 2)));
        let b = R::new(1.into(), BigInt::from((2 * k + 2) * (2 * k + 3)));
        term_c = term_c.mul(&r2).scale(&a);
        term_s = term_s.mul(&r2).scale(&b);
    }
    // Alternating series with decreasing terms (|r| < 1): the remainder is
    // bounded by the next term's magnitude.
    let bc = term_c.abs_hi();
    let bs = term_s.abs_hi();
    (
        cos.add(&Interval::new(-&bc, bc)),
        sin.add(&Interval::new(-&bs, bs)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn width(i: &Interval) -> f64 {
        num_float(&(i.hi() - i.lo()))
    }
    fn num_float(x: &R) -> f64 {
        let (n, d) = (x.numer(), x.denom());
        let shift = (d.bits() as i64 - 60).max(0);
        let n = n >> shift as usize;
        let d = d >> shift as usize;
        let nf = n.to_string().parse::<f64>().unwrap();
        let df = d.to_string().parse::<f64>().unwrap();
        nf / df
    }
    fn contains(i: &Interval, x: f64) -> bool {
        // x is a binary64 approximation: allow its rounding.
        let tol = 4.0 * f64::EPSILON * x.abs().max(1.0);
        num_float(i.lo()) <= x + tol && x - tol <= num_float(i.hi())
    }

    #[test]
    fn pi_is_tight_and_contains_the_binary64_approximation() {
        let p = pi();
        assert!(width(p) < 1e-50);
        assert!(contains(p, std::f64::consts::PI));
        // f64 PI is below the real pi.
        assert_eq!(
            p.cmp_exact(&R::from_float(std::f64::consts::PI).unwrap()),
            Some(Ordering::Greater)
        );
    }

    #[test]
    fn cos_sin_match_libm_and_the_pythagorean_identity() {
        for &x in &[
            0.5,
            -0.5,
            1.0,
            2.0,
            3.0,
            -3.0,
            std::f64::consts::FRAC_PI_2,
            std::f64::consts::PI,
            std::f64::consts::TAU,
            -std::f64::consts::TAU,
            10.0,
            100.25,
            -1234.5,
            1e-30,
        ] {
            let (c, s) = cos_sin(&R::from_float(x).unwrap());
            assert!(width(&c) < 1e-45 && width(&s) < 1e-45, "{x}");
            assert!(contains(&c, x.cos()) && contains(&s, x.sin()), "{x}");
            let one = c.square().add(&s.square());
            assert_eq!(
                one.cmp_exact(&R::from_integer(2.into())),
                Some(Ordering::Less)
            );
            assert!(contains(&one, 1.0));
        }
        // f64 TAU is below 2 pi, so sin(TAU_f64) is certainly negative.
        let (_, s) = cos_sin(&R::from_float(std::f64::consts::TAU).unwrap());
        assert_eq!(s.sign(), Some(Ordering::Less));
    }

    #[test]
    fn atan2_covers_all_quadrants_and_refuses_the_branch_cut() {
        for &(y, x) in &[
            (1.0, 2.0),
            (2.0, -1.0),
            (-2.0, -1.0),
            (-1.0, 3.0),
            (0.0, 1.0),
            (5.0, 0.0),
            (-5.0, 0.0),
            (1e-8, -1.0),
        ] {
            let a = atan2(&Interval::from_f64(y), &Interval::from_f64(x)).unwrap();
            assert!(width(&a) < 1e-40, "{y} {x}");
            assert!(contains(&a, f64::atan2(y, x)), "{y} {x}");
        }
        assert!(atan2(&Interval::from_f64(0.0), &Interval::from_f64(-1.0)).is_none());
    }

    #[test]
    fn sqrt_and_division_enclose() {
        let two = Interval::exact(int(2));
        let root = two.sqrt();
        assert!(contains(&root, std::f64::consts::SQRT_2));
        assert!(width(&root) < 1e-25);
        assert_eq!(root.square().cmp_exact(&int(2)), None);
        assert!(Interval::exact(one())
            .div(&Interval::new(int(-1), int(1)))
            .is_none());
    }
}

// ------------------------------------------------------------------ tiers

/// Arithmetic used by certified geometric checks. Implementations must only
/// return enclosures: every result contains the exact real value.
pub(crate) trait Real: Clone + std::fmt::Debug {
    /// The exact value of a finite binary64 number.
    fn exact_f64(x: f64) -> Self;
    /// An enclosure of an exact rational.
    fn from_r(x: &R) -> Self;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    fn neg(&self) -> Self;
    fn square(&self) -> Self;
    fn sqrt(&self) -> Self;
    fn div(&self, o: &Self) -> Option<Self>;
    fn sign(&self) -> Option<Ordering>;
    /// Certified comparison of two enclosed values.
    fn cmp(&self, o: &Self) -> Option<Ordering> {
        self.sub(o).sign()
    }
    /// (cos x, sin x) over an enclosed angle.
    fn cos_sin(x: &Self) -> (Self, Self);
    /// A rational strictly within reach of the enclosure (its midpoint).
    fn midpoint(&self) -> R;
    /// Half width, as a rational upper bound.
    fn radius(&self) -> R;
    /// Enlarge by +-w.
    fn widen(&self, w: &R) -> Self;
    fn atan2(y: &Self, x: &Self) -> Option<Self>;
}

impl Real for Interval {
    fn exact_f64(x: f64) -> Self {
        Self::from_f64(x)
    }
    fn from_r(x: &R) -> Self {
        Self::exact(x.clone())
    }
    fn add(&self, o: &Self) -> Self {
        Interval::add(self, o)
    }
    fn sub(&self, o: &Self) -> Self {
        Interval::sub(self, o)
    }
    fn mul(&self, o: &Self) -> Self {
        Interval::mul(self, o)
    }
    fn neg(&self) -> Self {
        Interval::neg(self)
    }
    fn square(&self) -> Self {
        Interval::square(self)
    }
    fn sqrt(&self) -> Self {
        Interval::sqrt(self)
    }
    fn div(&self, o: &Self) -> Option<Self> {
        Interval::div(self, o)
    }
    fn sign(&self) -> Option<Ordering> {
        Interval::sign(self)
    }
    fn cos_sin(x: &Self) -> (Self, Self) {
        if x.lo == x.hi {
            return cached_cos_sin(&x.lo);
        }
        // |d cos|, |d sin| <= |d angle| around the midpoint.
        let (c, s) = cached_cos_sin(&Real::midpoint(x));
        let w = Real::radius(x);
        (c.widen(&w), s.widen(&w))
    }
    fn midpoint(&self) -> R {
        (&self.lo + &self.hi) / int(2)
    }
    fn radius(&self) -> R {
        (&self.hi - &self.lo) / int(2)
    }
    fn widen(&self, w: &R) -> Self {
        Self::new(&self.lo - w, &self.hi + w)
    }
    fn atan2(y: &Self, x: &Self) -> Option<Self> {
        atan2(y, x)
    }
}

thread_local! {
    static TRIG: std::cell::RefCell<std::collections::HashMap<R, (Interval, Interval)>> =
        std::cell::RefCell::new(std::collections::HashMap::new());
}

/// Rational cos/sin, memoized per exact angle: models reuse few angles.
fn cached_cos_sin(x: &R) -> (Interval, Interval) {
    if let Some(hit) = TRIG.with(|t| t.borrow().get(x).cloned()) {
        return hit;
    }
    let value = cos_sin(x);
    TRIG.with(|t| {
        let mut t = t.borrow_mut();
        if t.len() > 65_536 {
            t.clear();
        }
        t.insert(x.clone(), value.clone());
    });
    value
}

fn next_up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits(1);
    }
    let bits = x.to_bits();
    f64::from_bits(if x > 0.0 { bits + 1 } else { bits - 1 })
}
fn next_down(x: f64) -> f64 {
    -next_up(-x)
}

/// Binary64 interval with one-ulp outward rounding after every operation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Fast {
    lo: f64,
    hi: f64,
}

impl Fast {
    fn out(lo: f64, hi: f64) -> Self {
        let lo = if lo.is_nan() {
            f64::NEG_INFINITY
        } else {
            next_down(lo)
        };
        let hi = if hi.is_nan() {
            f64::INFINITY
        } else {
            next_up(hi)
        };
        Self { lo, hi }
    }
    /// Outward binary64 enclosure of an exact rational interval.
    fn from_interval(i: &Interval) -> Self {
        let bound = |x: &R, up: bool| -> f64 {
            match crate::interval::enclose(|y| x.cmp(&R::from_float(y).unwrap()), "certified bound")
            {
                Ok(b) => {
                    if up {
                        b.upper()
                    } else {
                        b.lower()
                    }
                }
                Err(_) if up => f64::INFINITY,
                Err(_) => f64::NEG_INFINITY,
            }
        };
        Self {
            lo: bound(&i.lo, false),
            hi: bound(&i.hi, true),
        }
    }
}

impl Real for Fast {
    fn exact_f64(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn from_r(x: &R) -> Self {
        bracket(x).unwrap_or_else(|| Self::from_interval(&Interval::exact(x.clone())))
    }
    fn add(&self, o: &Self) -> Self {
        Self::bounds(sum(self.lo, o.lo), sum(self.hi, o.hi))
    }
    fn sub(&self, o: &Self) -> Self {
        Self::bounds(sum(self.lo, -o.hi), sum(self.hi, -o.lo))
    }
    fn mul(&self, o: &Self) -> Self {
        if self.lo == self.hi && o.lo == o.hi {
            let p = product(self.lo, o.lo);
            return Self::bounds(p, p);
        }
        let p = [
            product(self.lo, o.lo),
            product(self.lo, o.hi),
            product(self.hi, o.lo),
            product(self.hi, o.hi),
        ];
        if p.iter().any(|x| x.0.is_nan()) {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        let lo = p.iter().map(|x| down(*x)).fold(f64::INFINITY, f64::min);
        let hi = p.iter().map(|x| up(*x)).fold(f64::NEG_INFINITY, f64::max);
        Self { lo, hi }
    }
    fn neg(&self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
    fn square(&self) -> Self {
        let a = if self.lo >= 0.0 {
            self.lo
        } else if self.hi <= 0.0 {
            -self.hi
        } else {
            0.0
        };
        let b = (-self.lo).max(self.hi);
        Self {
            lo: down(product(a, a)),
            hi: up(product(b, b)),
        }
        .clamp_nonnegative()
    }
    fn sqrt(&self) -> Self {
        Self::out(self.lo.max(0.0).sqrt(), self.hi.max(0.0).sqrt()).clamp_nonnegative()
    }
    fn div(&self, o: &Self) -> Option<Self> {
        if !(o.lo > 0.0 || o.hi < 0.0) {
            return None;
        }
        let inv = Self::out(1.0 / o.hi, 1.0 / o.lo);
        Some(self.mul(&inv))
    }
    fn sign(&self) -> Option<Ordering> {
        if self.lo > 0.0 {
            Some(Ordering::Greater)
        } else if self.hi < 0.0 {
            Some(Ordering::Less)
        } else if self.lo == 0.0 && self.hi == 0.0 {
            Some(Ordering::Equal)
        } else {
            None
        }
    }
    fn cos_sin(x: &Self) -> (Self, Self) {
        fast_cos_sin(x)
    }
    fn midpoint(&self) -> R {
        let m = 0.5 * self.lo + 0.5 * self.hi;
        R::from_float(if m.is_finite() { m } else { 0.0 }).unwrap()
    }
    fn radius(&self) -> R {
        // Covers the distance from the (rounded) midpoint to both ends.
        let m = self.midpoint();
        match (R::from_float(self.lo), R::from_float(self.hi)) {
            (Some(lo), Some(hi)) => (&m - lo).max(hi - &m),
            _ => R::from_integer(BigInt::from(1) << 1100),
        }
    }
    fn widen(&self, w: &R) -> Self {
        let wf = Self::from_r(w).hi;
        Self::out(self.lo - wf, self.hi + wf)
    }
    fn atan2(y: &Self, x: &Self) -> Option<Self> {
        let to = |f: &Self| -> Option<Interval> {
            Some(Interval::new(R::from_float(f.lo)?, R::from_float(f.hi)?))
        };
        Some(Self::from_interval(&atan2(&to(y)?, &to(x)?)?))
    }
}

impl Fast {
    /// Rounded values with their exact error signs, as lower and upper bounds.
    fn bounds(lo: (f64, i8), hi: (f64, i8)) -> Self {
        if lo.0.is_nan() || hi.0.is_nan() {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        Self {
            lo: down(lo),
            hi: up(hi),
        }
    }
    fn clamp_nonnegative(self) -> Self {
        Self {
            lo: self.lo.max(0.0),
            hi: self.hi,
        }
    }
}

/// Rounded a + b and the sign of the exact error (TwoSum), or 2 when the sum
/// overflowed and the error is unknown.
fn sum(a: f64, b: f64) -> (f64, i8) {
    let s = a + b;
    if !s.is_finite() {
        return (s, if a.is_finite() && b.is_finite() { 2 } else { 0 });
    }
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, sign_of(err))
}

/// Rounded a * b and the sign of the exact error. The fused multiply-add
/// error is exact away from underflow; near it the sign is unknown (2).
fn product(a: f64, b: f64) -> (f64, i8) {
    let p = a * b;
    if a == 0.0 || b == 0.0 || !a.is_finite() || !b.is_finite() {
        return (p, 0);
    }
    if !p.is_finite() || p.abs() < f64::powi(2.0, -960) {
        return (p, 2);
    }
    (p, sign_of(a.mul_add(b, -p)))
}

fn sign_of(err: f64) -> i8 {
    if err > 0.0 {
        1
    } else if err < 0.0 {
        -1
    } else {
        0
    }
}

/// Lower bound of the exact value behind a rounded result.
fn down((x, err): (f64, i8)) -> f64 {
    if err == 0 || err == 1 {
        x
    } else {
        next_down(x)
    }
}

/// Upper bound of the exact value behind a rounded result.
fn up((x, err): (f64, i8)) -> f64 {
    if err == 0 || err == -1 {
        x
    } else {
        next_up(x)
    }
}

/// Leading 64 bits of a magnitude as f64 and the binary exponent dropped.
fn leading(v: &BigInt) -> (f64, i64) {
    let shift = (v.bits() as i64 - 64).max(0);
    let top = v.magnitude() >> shift as usize;
    (
        top.to_u64_digits().first().copied().unwrap_or(0) as f64,
        shift,
    )
}

/// The tightest binary64 bracket of a rational: a point when it is exactly
/// representable, else adjacent numbers. The quotient of leading bits is
/// within a few ulps; exact comparisons then fix the bracket. None outside
/// the comfortably normal range.
fn bracket(x: &R) -> Option<Fast> {
    if x.numer().sign() == Sign::NoSign {
        return Some(Fast::exact_f64(0.0));
    }
    let ((n, sn), (d, sd)) = (leading(x.numer()), leading(x.denom()));
    let e = sn - sd;
    if !(-900..=900).contains(&e) {
        return None;
    }
    let mut f = n / d * f64::powi(2.0, e as i32);
    if x.numer().sign() == Sign::Minus {
        f = -f;
    }
    if !f.is_finite() || f.abs() < 1e-250 || f.abs() > 1e250 {
        return None;
    }
    for _ in 0..8 {
        match R::from_float(f)?.cmp(x) {
            Ordering::Equal => return Some(Fast::exact_f64(f)),
            Ordering::Greater => f = next_down(f),
            Ordering::Less => {
                let up = next_up(f);
                if R::from_float(up)? >= *x {
                    return Some(if R::from_float(up)? == *x {
                        Fast::exact_f64(up)
                    } else {
                        Fast { lo: f, hi: up }
                    });
                }
                f = up;
            }
        }
    }
    None
}

/// pi/2 lies strictly between FRAC_PI_2 and the next binary64 number (checked
/// against the Machin enclosure in tests).
fn fast_half_pi() -> Fast {
    Fast {
        lo: std::f64::consts::FRAC_PI_2,
        hi: next_up(std::f64::consts::FRAC_PI_2),
    }
}

/// Certified binary64 cos/sin: reduce by an enclosure of pi/2, then sum Taylor
/// series in outward-rounded intervals with an explicit remainder term.
fn fast_cos_sin(x: &Fast) -> (Fast, Fast) {
    if x.lo == 0.0 && x.hi == 0.0 {
        return (Fast::exact_f64(1.0), Fast::exact_f64(0.0));
    }
    let whole = Fast { lo: -1.0, hi: 1.0 };
    if !(x.lo.is_finite() && x.hi.is_finite()) || x.lo.abs().max(x.hi.abs()) > 1.0e6 {
        return (whole, whole);
    }
    let middle = 0.5 * x.lo + 0.5 * x.hi;
    let k = (middle / std::f64::consts::FRAC_PI_2).round();
    let r = x.sub(&fast_half_pi().mul(&Fast::exact_f64(k)));
    let reach = r.lo.abs().max(r.hi.abs());
    if reach > 0.9 {
        return (whole, whole);
    }
    let r2 = r.square();
    let one = Fast::exact_f64(1.0);
    let (mut cos, mut sin) = (Fast::exact_f64(0.0), Fast::exact_f64(0.0));
    let (mut tc, mut ts) = (one, r);
    const TERMS: usize = 13;
    for j in 0..TERMS {
        if j % 2 == 0 {
            cos = cos.add(&tc);
            sin = sin.add(&ts);
        } else {
            cos = cos.sub(&tc);
            sin = sin.sub(&ts);
        }
        let a = ((2 * j + 1) * (2 * j + 2)) as f64;
        let b = ((2 * j + 2) * (2 * j + 3)) as f64;
        tc = tc.mul(&r2).div(&Fast::exact_f64(a)).unwrap();
        ts = ts.mul(&r2).div(&Fast::exact_f64(b)).unwrap();
    }
    // Alternating series with decreasing terms (|r| < 1): the remainder is
    // bounded by the next term's magnitude.
    let bound = |t: &Fast| {
        let m = t.lo.abs().max(t.hi.abs());
        Fast::out(-m, m)
    };
    let cos = cos.add(&bound(&tc));
    let sin = sin.add(&bound(&ts));
    let quadrant = (k.rem_euclid(4.0)) as i64;
    match quadrant {
        0 => (cos, sin),
        1 => (sin.neg(), cos),
        2 => (cos.neg(), sin.neg()),
        _ => (sin, cos.neg()),
    }
}

#[cfg(test)]
mod tier_tests {
    use super::*;

    #[test]
    fn binary64_half_pi_interval_contains_pi_over_two() {
        let half = pi().scale(&R::new(1.into(), 2.into()));
        let lo = R::from_float(std::f64::consts::FRAC_PI_2).unwrap();
        let hi = R::from_float(next_up(std::f64::consts::FRAC_PI_2)).unwrap();
        assert!(&lo < half.lo() && half.hi() < &hi);
    }

    #[test]
    fn fast_arithmetic_is_tight_when_exact_and_outward_otherwise() {
        let exact = |f: &Fast| (R::from_float(f.lo).unwrap(), R::from_float(f.hi).unwrap());
        let mut x = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let m = (x >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
            m * f64::powi(2.0, (x % 64) as i32 - 32)
        };
        for _ in 0..4000 {
            let (a, b) = (next(), next());
            let (fa, fb) = (Fast::exact_f64(a), Fast::exact_f64(b));
            let (ra, rb) = (R::from_float(a).unwrap(), R::from_float(b).unwrap());
            for (f, r) in [
                (fa.add(&fb), &ra + &rb),
                (fa.sub(&fb), &ra - &rb),
                (fa.mul(&fb), &ra * &rb),
                (fa.square(), &ra * &ra),
            ] {
                let (lo, hi) = exact(&f);
                // Exact results stay points; inexact ones sit strictly
                // inside adjacent binary64 numbers.
                if lo == hi {
                    assert_eq!(lo, r, "{a} {b}");
                } else {
                    assert!(lo < r && r < hi && f.hi == next_up(f.lo), "{a} {b}");
                }
            }
        }
        let one = Fast::exact_f64(1.0);
        assert_eq!(one.add(&one), Fast::exact_f64(2.0));
        assert_eq!(
            fast_cos_sin(&Fast::exact_f64(0.0)),
            (one, Fast::exact_f64(0.0))
        );
    }

    #[test]
    fn rational_brackets_are_adjacent_binary64_numbers() {
        let mut x = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for _ in 0..2000 {
            let n = (BigInt::from(next() as i64) * BigInt::from(next())) >> (next() % 90) as usize;
            let d = BigInt::from(next() | 1) << (next() % 70) as usize;
            let v = R::new(n, d);
            let b = bracket(&v).unwrap();
            let (lo, hi) = (R::from_float(b.lo).unwrap(), R::from_float(b.hi).unwrap());
            assert!(lo <= v && v <= hi);
            assert!(b.lo == b.hi || b.hi == next_up(b.lo));
            assert_eq!(b.lo == b.hi, lo == v);
        }
        assert_eq!(
            bracket(&R::from_float(0.1).unwrap()),
            Some(Fast::exact_f64(0.1))
        );
        assert_eq!(
            bracket(&R::new(1.into(), 3.into())).map(|b| b.hi == next_up(b.lo)),
            Some(true)
        );
    }

    #[test]
    fn fast_trig_encloses_the_rational_enclosure() {
        // Angles through +-150 radians, off any multiple of pi/2.
        for i in -50..=50 {
            let x = i as f64 * 2.97 + 0.001 * (i as f64).sin();
            let (c, s) = fast_cos_sin(&Fast::exact_f64(x));
            let (ce, se) = cos_sin(&R::from_float(x).unwrap());
            for (fast, exact) in [(c, ce), (s, se)] {
                assert!(R::from_float(fast.lo).unwrap() <= *exact.lo(), "{x}");
                assert!(R::from_float(fast.hi).unwrap() >= *exact.hi(), "{x}");
                assert!(fast.hi - fast.lo < 1e-13, "{x}");
            }
        }
    }
}
