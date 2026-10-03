//! Certified Gauss–Legendre quadrature of spline mass properties (F8 of
//! REVIEW_NOTES.md; the derivation is in MATHEMATICS.md). On a piece
//! `[a, b]` the `n`-point rule's node sum is exact up to `K_n (b - a)^(2n+1)
//! f_(2n)(ξ)`, with the Taylor coefficient `f_(2n)` enclosed over the whole
//! piece by interval Taylor arithmetic (`Series`); on a box the tensor rule
//! adds one such term per axis, each from a univariate series with the
//! other variable held as an interval. Pieces are halved while a remainder
//! is large against the piece's own absolute integral; every accepted piece
//! is an enclosure whatever the acceptance decides.
//!
//! Integrands are written once over `Num`, so the node values run on plain
//! enclosures (`Point`) and the remainders on series. Three integrals use
//! the rule: `∫ M/W^k` along rational pcurve pieces (`quotient_integral`),
//! `∫ F(u, v) du` along spline pcurves (`pcurve_integrals`), and the face
//! integrals of a spline surface by Green's theorem through its patches
//! (`spline_face`).
use super::bernstein::{c, derivative, pcurve_arcs, piece_of, r, ratio, span_arcs, Bern};
use super::Lp;
use crate::certified::{Interval as I, Real};
use crate::surface::ExactBezierSurface3;
use crate::topology::Curve2;
use crate::BSplineSurface3;
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::OnceLock;

/// Nodes of the rule per direction.
const NODES: usize = 8;
/// Series length: every coefficient up to the remainder's, the `2n`-th.
const LENGTH: usize = 2 * NODES + 1;
/// A piece is accepted when each remainder is at most this share (`2^-44`)
/// of the piece's absolute integral, scaled by its relative size.
const RELATIVE: f64 = 1.0 / 17_592_186_044_416.0;
/// Halvings per direction before a piece is accepted as it is.
const DEPTH: u32 = 12;
/// Pieces (or boxes) one integral may examine, and the halvings (summed
/// over the directions) past which a series still undefined means a
/// singularity, not a wide enclosure: then the rule gives up and the
/// first-order route runs (a degenerate patch edge, where `|N|` vanishes,
/// would otherwise be refined down to the depth limit all along it).
const WORK: usize = 2048;
const SINGULAR: u32 = 12;

fn zero<T: Real>() -> T {
    T::exact_f64(0.0)
}

/// The arithmetic the integrands need, shared by plain enclosures (the
/// rule's node values) and truncated Taylor series (its remainders).
pub(super) trait Num<T: Real>: Clone {
    /// The constant `x`, shaped like `self`.
    fn lift(&self, x: &T) -> Self;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn neg(&self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    /// The square, certainly nonnegative.
    fn square(&self) -> Self;
    fn scale(&self, s: &T) -> Self;
    fn shift(&self, s: &T) -> Self;
    /// `None` where `o` may vanish.
    fn div(&self, o: &Self) -> Option<Self>;
    /// `None` unless `self` is certainly positive.
    fn sqrt(&self) -> Option<Self>;
    fn cos_sin(&self) -> (Self, Self);
    fn powi(&self, n: u8) -> Self {
        (0..n).fold(self.lift(&T::exact_f64(1.0)), |acc, _| acc.mul(self))
    }
    /// The middle of the value's enclosure at the base (a binary64 choice
    /// between equal forms, S9f.2b).
    fn mid(&self) -> f64;
    /// Whether the value's enclosure at the base is about a point, within
    /// rounding (S9f.2b: a de Casteljau step's form).
    fn sharp(&self) -> bool;
}

fn middle<T: Real>(x: &T) -> f64 {
    let (lo, hi) = x.bounds_f64();
    0.5 * lo + 0.5 * hi
}

/// An enclosure about a point (a node's or a jet's base, within rounding),
/// not over a range.
fn sharp<T: Real>(x: &T) -> bool {
    let (lo, hi) = x.bounds_f64();
    hi - lo <= 1e-12 * lo.abs().max(hi.abs()).max(1.0)
}

/// A plain enclosure.
#[derive(Clone, Debug)]
pub(super) struct Point<T>(T);

impl<T: Real> Num<T> for Point<T> {
    fn lift(&self, x: &T) -> Self {
        Point(x.clone())
    }
    fn add(&self, o: &Self) -> Self {
        Point(self.0.add(&o.0))
    }
    fn sub(&self, o: &Self) -> Self {
        Point(self.0.sub(&o.0))
    }
    fn neg(&self) -> Self {
        Point(self.0.neg())
    }
    fn mul(&self, o: &Self) -> Self {
        Point(self.0.mul(&o.0))
    }
    fn square(&self) -> Self {
        Point(self.0.square())
    }
    fn scale(&self, s: &T) -> Self {
        Point(self.0.mul(s))
    }
    fn shift(&self, s: &T) -> Self {
        Point(self.0.add(s))
    }
    fn div(&self, o: &Self) -> Option<Self> {
        if !matches!(o.0.sign(), Some(Ordering::Greater | Ordering::Less)) {
            return None;
        }
        Some(Point(self.0.div(&o.0)?))
    }
    fn sqrt(&self) -> Option<Self> {
        (self.0.sign() == Some(Ordering::Greater)).then(|| Point(self.0.sqrt()))
    }
    fn cos_sin(&self) -> (Self, Self) {
        let (co, si) = T::cos_sin(&self.0);
        (Point(co), Point(si))
    }
    fn mid(&self) -> f64 {
        middle(&self.0)
    }
    fn sharp(&self) -> bool {
        sharp(&self.0)
    }
}

/// Taylor jets (`crate::jet`) as the integrands' numbers (S9f.2b: a spline
/// wall's meeting evaluated once for its jets and its quadrature).
impl<T: Real> Num<T> for crate::jet::Jet<T> {
    fn lift(&self, x: &T) -> Self {
        Self::constant(x.clone(), self.order())
    }
    fn add(&self, o: &Self) -> Self {
        crate::jet::Jet::add(self, o)
    }
    fn sub(&self, o: &Self) -> Self {
        crate::jet::Jet::sub(self, o)
    }
    fn neg(&self) -> Self {
        crate::jet::Jet::neg(self)
    }
    fn mul(&self, o: &Self) -> Self {
        crate::jet::Jet::mul(self, o)
    }
    fn square(&self) -> Self {
        crate::jet::Jet::square(self)
    }
    fn scale(&self, s: &T) -> Self {
        crate::jet::Jet::scale(self, s)
    }
    fn shift(&self, s: &T) -> Self {
        self.add_constant(s)
    }
    fn div(&self, o: &Self) -> Option<Self> {
        crate::jet::Jet::div(self, o)
    }
    fn sqrt(&self) -> Option<Self> {
        crate::jet::Jet::sqrt(self)
    }
    fn cos_sin(&self) -> (Self, Self) {
        crate::jet::Jet::cos_sin(self)
    }
    fn mid(&self) -> f64 {
        middle(&self.c[0])
    }
    fn sharp(&self) -> bool {
        sharp(&self.c[0])
    }
}

/// A truncated Taylor series `Σ_(k<len) c_k ε^k` in one variable, standing
/// for a function at `x + ε`; coefficients past `c` are exactly zero (a
/// constant keeps one, a linear argument two).
#[derive(Clone, Debug)]
pub(super) struct Series<T> {
    c: Vec<T>,
    len: usize,
}

impl<T: Real> Series<T> {
    fn constant(x: T, len: usize) -> Self {
        Self { c: vec![x], len }
    }

    fn variable(x: T, len: usize) -> Self {
        let mut c = vec![x];
        if len > 1 {
            c.push(T::exact_f64(1.0));
        }
        Self { c, len }
    }

    fn coefficient(&self, k: usize) -> T {
        self.c.get(k).cloned().unwrap_or_else(zero)
    }
}

impl<T: Real> Num<T> for Series<T> {
    fn lift(&self, x: &T) -> Self {
        Self::constant(x.clone(), self.len)
    }

    fn add(&self, o: &Self) -> Self {
        let n = self.c.len().max(o.c.len());
        let c = (0..n)
            .map(|k| match (self.c.get(k), o.c.get(k)) {
                (Some(a), Some(b)) => a.add(b),
                (Some(a), None) | (None, Some(a)) => a.clone(),
                (None, None) => unreachable!("k is below one length"),
            })
            .collect();
        Self { c, len: self.len }
    }

    fn sub(&self, o: &Self) -> Self {
        self.add(&o.neg())
    }

    fn neg(&self) -> Self {
        Self {
            c: self.c.iter().map(|x| x.neg()).collect(),
            len: self.len,
        }
    }

    fn mul(&self, o: &Self) -> Self {
        let n = (self.c.len() + o.c.len() - 1).min(self.len);
        let c = (0..n)
            .map(|k| {
                let lo = k.saturating_sub(o.c.len() - 1);
                let hi = k.min(self.c.len() - 1);
                (lo..=hi).fold(zero(), |s: T, i| s.add(&self.c[i].mul(&o.c[k - i])))
            })
            .collect();
        Self { c, len: self.len }
    }

    /// With each diagonal product squared (sign-exact, so a sum of squares
    /// stays certainly nonnegative).
    fn square(&self) -> Self {
        let m = self.c.len();
        let n = (2 * m - 1).min(self.len);
        let two = T::exact_f64(2.0);
        let c = (0..n)
            .map(|k| {
                let mut s = zero::<T>();
                let mut i = k.saturating_sub(m - 1);
                while 2 * i < k {
                    s = s.add(&self.c[i].mul(&self.c[k - i]));
                    i += 1;
                }
                s = s.mul(&two);
                if k % 2 == 0 {
                    s = s.add(&self.c[k / 2].square());
                }
                s
            })
            .collect();
        Self { c, len: self.len }
    }

    fn scale(&self, s: &T) -> Self {
        Self {
            c: self.c.iter().map(|x| x.mul(s)).collect(),
            len: self.len,
        }
    }

    fn shift(&self, s: &T) -> Self {
        let mut c = self.c.clone();
        c[0] = c[0].add(s);
        Self { c, len: self.len }
    }

    fn div(&self, o: &Self) -> Option<Self> {
        if !matches!(o.c[0].sign(), Some(Ordering::Greater | Ordering::Less)) {
            return None;
        }
        if o.c.len() == 1 {
            let c = self
                .c
                .iter()
                .map(|x| x.div(&o.c[0]))
                .collect::<Option<_>>()?;
            return Some(Self { c, len: self.len });
        }
        let mut q: Vec<T> = Vec::with_capacity(self.len);
        for k in 0..self.len {
            let mut s = self.coefficient(k);
            for j in 1..=k.min(o.c.len() - 1) {
                s = s.sub(&o.c[j].mul(&q[k - j]));
            }
            q.push(s.div(&o.c[0])?);
        }
        Some(Self {
            c: q,
            len: self.len,
        })
    }

    fn sqrt(&self) -> Option<Self> {
        if self.c[0].sign() != Some(Ordering::Greater) {
            return None;
        }
        let r0 = self.c[0].sqrt();
        if self.c.len() == 1 {
            return Some(Self::constant(r0, self.len));
        }
        let twice = r0.mul(&T::exact_f64(2.0));
        let mut out = vec![r0];
        for k in 1..self.len {
            let mut s = self.coefficient(k);
            for j in 1..k {
                s = s.sub(&out[j].mul(&out[k - j]));
            }
            out.push(s.div(&twice)?);
        }
        Some(Self {
            c: out,
            len: self.len,
        })
    }

    fn cos_sin(&self) -> (Self, Self) {
        let (c0, s0) = T::cos_sin(&self.c[0]);
        if self.c.len() == 1 {
            return (Self::constant(c0, self.len), Self::constant(s0, self.len));
        }
        let (mut co, mut si) = (vec![c0], vec![s0]);
        for k in 1..self.len {
            let (mut a, mut b) = (zero::<T>(), zero::<T>());
            for j in 1..=k.min(self.c.len() - 1) {
                let ja = self.c[j].mul(&T::exact_f64(j as f64));
                a = a.add(&ja.mul(&si[k - j]));
                b = b.add(&ja.mul(&co[k - j]));
            }
            let inverse = c::<T>(&ratio(1, k as i64));
            co.push(a.mul(&inverse).neg());
            si.push(b.mul(&inverse));
        }
        let len = self.len;
        (Self { c: co, len }, Self { c: si, len })
    }

    fn mid(&self) -> f64 {
        middle(&self.c[0])
    }
    fn sharp(&self) -> bool {
        sharp(&self.c[0])
    }
}

/// A Bernstein polynomial on `[0, 1]` with `Num` coefficients at `t`, by
/// de Casteljau (`a + t (b - a)` at every step).
fn casteljau<T: Real, N: Num<T>>(mut row: Vec<N>, t: &N) -> N {
    for last in (1..row.len()).rev() {
        for i in 0..last {
            row[i] = row[i].add(&t.mul(&row[i + 1].sub(&row[i])));
        }
    }
    row.swap_remove(0)
}

fn bernstein<T: Real, N: Num<T>>(b: &[T], t: &N) -> N {
    casteljau(b.iter().map(|x| t.lift(x)).collect(), t)
}

/// The rule's nodes and weights on `[-1, 1]` and `K_n`, as rational
/// enclosures, computed once.
struct ExactRule {
    nodes: Vec<[R; 2]>,
    weights: Vec<[R; 2]>,
    remainder: R,
}

/// `(P_n(x), P_(n-1)(x))` by the three-term recurrence.
fn legendre<T: Real>(n: usize, x: &T) -> (T, T) {
    let (mut previous, mut current) = (T::exact_f64(1.0), x.clone());
    for k in 1..n {
        let next = c::<T>(&ratio(2 * k as i64 + 1, k as i64 + 1))
            .mul(x)
            .mul(&current)
            .sub(&c::<T>(&ratio(k as i64, k as i64 + 1)).mul(&previous));
        previous = current;
        current = next;
    }
    (current, previous)
}

/// `(P_n(x), P_(n-1)(x))` in binary64, for Newton's iteration.
fn legendre_f64(n: usize, x: f64) -> (f64, f64) {
    let (mut previous, mut current) = (1.0, x);
    for k in 1..n {
        let k = k as f64;
        let next = ((2.0 * k + 1.0) * x * current - k * previous) / (k + 1.0);
        previous = current;
        current = next;
    }
    (current, previous)
}

/// The sign of `P_n(x)`, exactly.
fn legendre_sign(n: usize, x: &R) -> Ordering {
    let (mut previous, mut current) = (ratio(1, 1), x.clone());
    for k in 1..n {
        let next = (ratio(2 * k as i64 + 1, 1) * x * &current - ratio(k as i64, 1) * &previous)
            / ratio(k as i64 + 1, 1);
        previous = current;
        current = next;
    }
    current.cmp(&ratio(0, 1))
}

fn exact_rule() -> ExactRule {
    let n = NODES;
    let mut nodes = Vec::with_capacity(n);
    let mut weights = Vec::with_capacity(n);
    for i in 0..n {
        // Newton's iteration in binary64 from the classical guess.
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (p, q) = legendre_f64(n, x);
            let step = p / (n as f64 * (q - x * p) / (1.0 - x * x));
            x -= step;
            if step.abs() < 1e-17 {
                break;
            }
        }
        // Widened until the exact values of P_n at the ends differ in sign,
        // then bisected exactly down to 2^-100, so that the nodes and
        // weights lift into binary64 as adjacent values.
        let mut delta = f64::EPSILON;
        let (mut lo, mut hi, low) = loop {
            let (lo, hi) = (r(x) - r(delta), r(x) + r(delta));
            let (a, b) = (legendre_sign(n, &lo), legendre_sign(n, &hi));
            if a != b && a != Ordering::Equal && b != Ordering::Equal {
                break (lo, hi, a);
            }
            delta *= 2.0;
            assert!(delta < 1e-6, "a Legendre root is bracketed");
        };
        let fine = ratio(1, 1) / R::from_integer(num_bigint::BigInt::from(1) << 100);
        while &hi - &lo > fine {
            let mid = (&lo + &hi) / ratio(2, 1);
            match legendre_sign(n, &mid) {
                s if s == low => lo = mid,
                Ordering::Equal => break,
                _ => hi = mid,
            }
        }
        // w = 2 (1 - x^2) / (n^2 (P_(n-1) - x P_n)^2) over the bracket.
        let x = I::new(lo.clone(), hi.clone());
        let (p, q) = legendre::<I>(n, &x);
        let d = q.sub(&x.mul(&p));
        let w = I::exact(ratio(1, 1))
            .sub(&x.square())
            .mul(&I::exact(ratio(2, 1)))
            .div(&d.square().mul(&I::exact(ratio((n * n) as i64, 1))))
            .expect("a Legendre root is not a root of P_(n-1) - x P_n");
        nodes.push([lo, hi]);
        weights.push([w.lo().clone(), w.hi().clone()]);
    }
    for pair in nodes.windows(2) {
        assert!(pair[0][0] > pair[1][1], "disjoint brackets, one per root");
    }
    let factorial = |k: usize| (1..=k).fold(ratio(1, 1), |a, j| a * ratio(j as i64, 1));
    let (f, g) = (factorial(n), factorial(2 * n));
    let remainder = &f * &f * &f * &f / (ratio(2 * n as i64 + 1, 1) * &g * &g);
    ExactRule {
        nodes,
        weights,
        remainder,
    }
}

/// The rule lifted into a tier.
struct Rule<T> {
    nodes: Vec<T>,
    weights: Vec<T>,
    remainder: T,
}

fn rule<T: Real>() -> Rule<T> {
    static EXACT: OnceLock<ExactRule> = OnceLock::new();
    let exact = EXACT.get_or_init(exact_rule);
    let hull = |[lo, hi]: &[R; 2]| c::<T>(lo).union(&c(hi));
    Rule {
        nodes: exact.nodes.iter().map(hull).collect(),
        weights: exact.weights.iter().map(hull).collect(),
        remainder: c(&exact.remainder),
    }
}

/// The rule's points and weights on `[lo, hi]`.
fn scaled_rule<T: Real>(rule: &Rule<T>, lo: f64, hi: f64) -> (Vec<T>, Vec<T>) {
    let half = T::exact_f64(0.5);
    let mid = T::exact_f64(lo).add(&T::exact_f64(hi)).mul(&half);
    let h = T::exact_f64(hi).sub(&T::exact_f64(lo)).mul(&half);
    (
        rule.nodes.iter().map(|x| mid.add(&h.mul(x))).collect(),
        rule.weights.iter().map(|w| h.mul(w)).collect(),
    )
}

/// The largest magnitude in an enclosure.
fn magnitude<T: Real>(x: &T) -> f64 {
    let (lo, hi) = x.bounds_f64();
    lo.abs().max(hi.abs())
}

fn width<T: Real>(x: &T) -> f64 {
    let (lo, hi) = x.bounds_f64();
    hi - lo
}

fn power<T: Real>(x: f64, k: usize) -> T {
    let x = T::exact_f64(x);
    (0..k).fold(T::exact_f64(1.0), |p, _| p.mul(&x))
}

fn accumulate<T: Real>(total: &mut [T], values: &[T]) {
    for (t, v) in total.iter_mut().zip(values) {
        *t = t.add(v);
    }
}

/// A vector integrand of one variable on `[0, 1]`: `None` where it is not
/// certainly smooth.
trait Integrand1<T: Real> {
    fn at<N: Num<T>>(&self, t: &N) -> Option<Vec<N>>;
}

/// A vector integrand of `(τ, σ)` on the unit square.
trait Integrand2<T: Real> {
    fn at<N: Num<T>>(&self, tau: &N, sigma: &N) -> Option<Vec<N>>;
}

/// The node sum of every component, and the sum of their magnitudes.
fn sum_nodes<T: Real>(
    values: impl Iterator<Item = Option<(Vec<Point<T>>, T)>>,
) -> Option<(Vec<T>, Vec<f64>)> {
    let (mut sum, mut absolute): (Vec<T>, Vec<f64>) = (Vec::new(), Vec::new());
    for item in values {
        let (values, w) = item?;
        if sum.is_empty() {
            sum = vec![zero(); values.len()];
            absolute = vec![0.0; values.len()];
        }
        for ((s, a), v) in sum.iter_mut().zip(absolute.iter_mut()).zip(&values) {
            let term = v.0.mul(&w);
            *a += magnitude(&term);
            *s = s.add(&term);
        }
    }
    Some((sum, absolute))
}

/// Whether every remainder `e` is within its budget: `size` times the
/// larger of `RELATIVE` of the piece's absolute integral and four times the
/// first rule's rounding width.
fn within(e: &[f64], size: f64, scale: &[(f64, f64)]) -> bool {
    e.iter()
        .zip(scale)
        .all(|(e, (share, rounding))| *e <= size * (RELATIVE * share).max(4.0 * rounding))
}

/// `∫_0^1` of every component of `f`, enclosed.
fn integrate_1d<T: Real, F: Integrand1<T>>(f: &F) -> Option<Vec<T>> {
    let rule = rule::<T>();
    let nodes = |a: f64, b: f64| {
        let (points, weights) = scaled_rule(&rule, a, b);
        sum_nodes(
            points
                .into_iter()
                .zip(weights)
                .map(|(t, w)| Some((f.at(&Point(t))?, w))),
        )
    };
    // K_n (b - a)^(2n+1) f_(2n)([a, b]) per component.
    let remainder = |a: f64, b: f64| -> Option<Vec<T>> {
        let x = T::exact_f64(a).union(&T::exact_f64(b));
        let scale = rule.remainder.mul(&power(b - a, 2 * NODES + 1));
        Some(
            f.at(&Series::variable(x, LENGTH))?
                .iter()
                .map(|s| s.coefficient(2 * NODES).mul(&scale))
                .collect(),
        )
    };
    let (first, absolute) = nodes(0.0, 1.0)?;
    let scale: Vec<(f64, f64)> = absolute
        .iter()
        .zip(&first)
        .map(|(s, q)| (*s, width(q)))
        .collect();
    let mut total: Vec<T> = vec![zero(); first.len()];
    let mut stack = vec![(0.0_f64, 1.0_f64, 0_u32)];
    let mut work = 0;
    while let Some((a, b, depth)) = stack.pop() {
        work += 1;
        if work > WORK {
            return None;
        }
        // A series undefined over the whole piece (an enclosure too wide to
        // exclude a zero) is retried on its halves; the node sum is formed
        // only for accepted pieces.
        let rest = remainder(a, b);
        if rest.is_none() && depth >= SINGULAR {
            return None;
        }
        let accepted = rest
            .as_ref()
            .is_some_and(|e| within(&e.iter().map(magnitude).collect::<Vec<_>>(), b - a, &scale));
        if accepted || (depth == DEPTH && rest.is_some()) {
            let sum = if (a, b) == (0.0, 1.0) {
                first.clone()
            } else {
                nodes(a, b)?.0
            };
            for (t, (q, e)) in total.iter_mut().zip(sum.iter().zip(rest?.iter())) {
                *t = t.add(&q.add(e));
            }
            continue;
        }
        if depth == DEPTH {
            return None;
        }
        let m = 0.5 * (a + b);
        stack.push((a, m, depth + 1));
        stack.push((m, b, depth + 1));
    }
    Some(total)
}

/// `∫_0^1 ∫_0^1` of every component of `f(τ, σ)`, enclosed by the tensor
/// rule and its two remainders.
fn integrate_2d<T: Real, F: Integrand2<T>>(f: &F) -> Option<Vec<T>> {
    type Box2 = [[f64; 2]; 2];
    let rule = rule::<T>();
    let nodes = |[[a, b], [c0, d]]: Box2| {
        let (taus, wt) = scaled_rule(&rule, a, b);
        let (sigmas, ws) = scaled_rule(&rule, c0, d);
        let pairs = taus.iter().zip(&wt).flat_map(|(tau, w1)| {
            sigmas.iter().zip(&ws).map(move |(sigma, w2)| {
                Some((
                    f.at(&Point(tau.clone()), &Point(sigma.clone()))?,
                    w1.mul(w2),
                ))
            })
        });
        sum_nodes(pairs)
    };
    // The remainders along τ and along σ, per component.
    let remainders = |[[a, b], [c0, d]]: Box2| -> Option<(Vec<T>, Vec<T>)> {
        let (x, y) = (
            T::exact_f64(a).union(&T::exact_f64(b)),
            T::exact_f64(c0).union(&T::exact_f64(d)),
        );
        let along_tau = rule
            .remainder
            .mul(&power(b - a, 2 * NODES + 1))
            .mul(&T::exact_f64(d - c0));
        let along_sigma = rule
            .remainder
            .mul(&T::exact_f64(b - a))
            .mul(&power(d - c0, 2 * NODES + 1));
        let tau = f.at(
            &Series::variable(x.clone(), LENGTH),
            &Series::constant(y.clone(), LENGTH),
        )?;
        let sigma = f.at(&Series::constant(x, LENGTH), &Series::variable(y, LENGTH))?;
        let last = |s: &Series<T>, k: &T| s.coefficient(2 * NODES).mul(k);
        Some((
            tau.iter().map(|s| last(s, &along_tau)).collect(),
            sigma.iter().map(|s| last(s, &along_sigma)).collect(),
        ))
    };
    let whole: Box2 = [[0.0, 1.0], [0.0, 1.0]];
    let (first, absolute) = nodes(whole)?;
    let scale: Vec<(f64, f64)> = absolute
        .iter()
        .zip(&first)
        .map(|(s, q)| (*s, width(q)))
        .collect();
    let mut total: Vec<T> = vec![zero(); first.len()];
    let mut stack = vec![(whole, [0_u32; 2])];
    let mut work = 0;
    while let Some((bx, depth)) = stack.pop() {
        work += 1;
        if work > WORK {
            return None;
        }
        let area = (bx[0][1] - bx[0][0]) * (bx[1][1] - bx[1][0]);
        // A series undefined over the whole box (an enclosure too wide to
        // exclude a zero) is retried on halves across the less divided
        // side; the node sum is formed only for accepted boxes.
        let rest = remainders(bx);
        if rest.is_none() && depth[0] + depth[1] >= SINGULAR {
            return None;
        }
        let (accepted, axis) = match &rest {
            Some((et, es)) => {
                let t: Vec<f64> = et.iter().map(magnitude).collect();
                let s: Vec<f64> = es.iter().map(magnitude).collect();
                let both: Vec<f64> = t.iter().zip(&s).map(|(a, b)| a + b).collect();
                // Halve across the direction whose remainder is larger
                // against its budget.
                let ratio = |e: &[f64]| {
                    e.iter()
                        .zip(&scale)
                        .map(|(e, (share, _))| e / (RELATIVE * share).max(f64::MIN_POSITIVE))
                        .fold(0.0_f64, f64::max)
                };
                (
                    within(&both, area, &scale),
                    usize::from(ratio(&s) > ratio(&t)),
                )
            }
            None => (false, usize::from(depth[1] < depth[0])),
        };
        let axis = if depth[axis] == DEPTH { 1 - axis } else { axis };
        if accepted || depth[axis] == DEPTH {
            let (et, es) = rest?;
            let sum = if bx == whole {
                first.clone()
            } else {
                nodes(bx)?.0
            };
            for (t, ((q, a), b)) in total.iter_mut().zip(sum.iter().zip(&et).zip(&es)) {
                *t = t.add(&q.add(a).add(b));
            }
            continue;
        }
        let [lo, hi] = bx[axis];
        let m = 0.5 * (lo + hi);
        let (mut left, mut right, mut next) = (bx, bx, depth);
        left[axis] = [lo, m];
        right[axis] = [m, hi];
        next[axis] += 1;
        stack.push((left, next));
        stack.push((right, next));
    }
    Some(total)
}

/// `M/W^k` of Bernstein polynomials.
struct Quotient<'a, T> {
    m: &'a Bern<T>,
    w: &'a Bern<T>,
    k: u8,
}

impl<T: Real> Integrand1<T> for Quotient<'_, T> {
    fn at<N: Num<T>>(&self, t: &N) -> Option<Vec<N>> {
        let denominator = bernstein(self.w, t).powi(self.k);
        Some(vec![bernstein(self.m, t).div(&denominator)?])
    }
}

/// `∫_0^1 M/W^k` of Bernstein polynomials, enclosed (a rational pcurve
/// piece's flux or mass integrand); `None` where the rule cannot run.
pub(super) fn quotient_integral<T: Real>(m: &Bern<T>, w: &Bern<T>, k: i32) -> Option<T> {
    let k = u8::try_from(k).ok()?;
    integrate_1d(&Quotient { m, w, k })?.pop()
}

/// A piece of a pcurve as homogeneous Bernstein coordinates `(U, V, W)` on
/// `[0, 1]` and the derivatives `U'`, `W'`; a uniform weight is kept as a
/// constant.
struct Piece<T> {
    u: Bern<T>,
    v: Bern<T>,
    w: Bern<T>,
    du: Bern<T>,
    dw: Bern<T>,
}

impl<T: Real> Piece<T> {
    fn new(u: Bern<T>, v: Bern<T>, w: Bern<T>) -> Self {
        let (du, dw) = (derivative(&u), derivative(&w));
        Self { u, v, w, du, dw }
    }

    /// `(u, v, u')` at `t`.
    fn at<N: Num<T>>(&self, t: &N) -> Option<[N; 3]> {
        let (u, v, w) = (
            bernstein(&self.u, t),
            bernstein(&self.v, t),
            bernstein(&self.w, t),
        );
        let (du, dw) = (bernstein(&self.du, t), bernstein(&self.dw, t));
        let slope = du.mul(&w).sub(&u.mul(&dw)).div(&w.square())?;
        Some([u.div(&w)?, v.div(&w)?, slope])
    }
}

/// The exact homogeneous coordinates of a pcurve arc's piece as a `Piece`,
/// optionally mapped affinely first (`X ↦ (X - lo W)/(hi - lo)` per axis).
fn exact_piece<T: Real>(h: &[Vec<R>; 4], local: Option<&[[R; 2]; 2]>) -> Piece<T> {
    let map = |k: usize| -> Bern<T> {
        h[k].iter()
            .zip(&h[3])
            .map(|(x, w)| match local {
                Some(d) => c(&((x - &d[k][0] * w) / (&d[k][1] - &d[k][0]))),
                None => c(x),
            })
            .collect()
    };
    let w: Bern<T> = if h[3].iter().all(|x| *x == h[3][0]) {
        vec![c(&h[3][0])]
    } else {
        h[3].iter().map(c).collect()
    };
    Piece::new(map(0), map(1), w)
}

/// A vector function of a pcurve's point `(u, v)` (every `F` whose
/// `∫ F du` is wanted).
pub(super) trait AlongPcurve<T: Real> {
    fn at<N: Num<T>>(&self, u: &N, v: &N) -> Option<Vec<N>>;
}

/// `F(u(τ), v(τ)) u'(τ)` along a piece.
struct OnPiece<'a, T, F> {
    piece: &'a Piece<T>,
    f: &'a F,
}

impl<T: Real, F: AlongPcurve<T>> Integrand1<T> for OnPiece<'_, T, F> {
    fn at<N: Num<T>>(&self, t: &N) -> Option<Vec<N>> {
        let [u, v, slope] = self.piece.at(t)?;
        Some(self.f.at(&u, &v)?.iter().map(|x| x.mul(&slope)).collect())
    }
}

/// `∫_0^1 F(u(τ), v(τ)) u'(τ) dτ` over a spline pcurve translated by
/// `-about` (exactly, on its homogeneous poles), for every component of
/// `F`, enclosed; `None` where the rule cannot run.
pub(super) fn pcurve_integrals<T: Real, F: AlongPcurve<T>>(
    span: &crate::topology::SplineSpan<crate::BSplineCurve2>,
    about: &[R; 2],
    count: usize,
    f: &F,
) -> Option<Vec<T>> {
    let mut total = vec![zero::<T>(); count];
    for (_, _, arc) in &span_arcs(span)? {
        let poles = arc.homogeneous_poles();
        let h: [Vec<R>; 4] = std::array::from_fn(|k| {
            poles
                .iter()
                .map(|p| match k {
                    0 | 1 => &p[k] - &about[k] * &p[3],
                    _ => p[k].clone(),
                })
                .collect()
        });
        // Constant u: nothing to integrate.
        if h[0]
            .iter()
            .zip(&h[3])
            .all(|(x, w)| x * &h[3][0] == &h[0][0] * w)
        {
            continue;
        }
        let piece = exact_piece::<T>(&h, None);
        accumulate(&mut total, &integrate_1d(&OnPiece { piece: &piece, f })?);
    }
    Some(total)
}

/// A Bézier patch translated by `-centre` exactly and lifted: per
/// homogeneous coordinate, the nets of `h`, `h_ū` and `h_v̄` in the patch's
/// local coordinates.
struct Net<T> {
    domain: [[R; 2]; 2],
    nets: Vec<[Vec<Vec<T>>; 3]>,
}

fn net<T: Real>(q: &ExactBezierSurface3, centre: &[R; 3]) -> Net<T> {
    let [du, dv] = q.degrees();
    let poles = q.homogeneous_poles();
    let nets = (0..4)
        .map(|k| {
            let h: Vec<Vec<T>> = (0..=du)
                .map(|i| {
                    (0..=dv)
                        .map(|j| {
                            let p = &poles[i * (dv + 1) + j];
                            if k < 3 {
                                c(&(&p[k] - &centre[k] * &p[3]))
                            } else {
                                c(&p[3])
                            }
                        })
                        .collect()
                })
                .collect();
            let hu: Vec<Vec<T>> = if du == 0 {
                vec![vec![zero(); dv + 1]]
            } else {
                let m = T::exact_f64(du as f64);
                h.windows(2)
                    .map(|w| {
                        w[1].iter()
                            .zip(&w[0])
                            .map(|(a, b)| a.sub(b).mul(&m))
                            .collect()
                    })
                    .collect()
            };
            let hv: Vec<Vec<T>> = h.iter().map(derivative).collect();
            [h, hu, hv]
        })
        .collect();
    Net {
        domain: q.domain().clone(),
        nets,
    }
}

/// A tensor Bernstein net at `(ū, v̄)`: rows in `v̄`, then `ū`.
fn tensor<T: Real, N: Num<T>>(net: &[Vec<T>], u: &N, v: &N) -> N {
    casteljau(net.iter().map(|row| bernstein(row, v)).collect(), u)
}

/// The fourteen mass integrands (or the four with `|N|`) of a patch at
/// `(ū, v̄)`, with `N̄ = S_ū × S_v̄`, relative to the origin (`rest` is the
/// centre the net is translated by, minus the origin).
fn patch_integrands<T: Real, N: Num<T>>(
    net: &Net<T>,
    u: &N,
    v: &N,
    rest: &[T; 3],
    third: Option<&T>,
) -> Option<Vec<N>> {
    let h: Vec<[N; 3]> = net
        .nets
        .iter()
        .map(|n| std::array::from_fn(|d| tensor(&n[d], u, v)))
        .collect();
    let w = &h[3];
    let (mut p, mut su, mut sv) = (Vec::new(), Vec::new(), Vec::new());
    for (k, hk) in h.iter().take(3).enumerate() {
        let s = hk[0].div(&w[0])?;
        su.push(hk[1].sub(&s.mul(&w[1])).div(&w[0])?);
        sv.push(hk[2].sub(&s.mul(&w[2])).div(&w[0])?);
        p.push(s.shift(&rest[k]));
    }
    let n: Vec<N> = (0..3)
        .map(|k| {
            let (a, b) = ((k + 1) % 3, (k + 2) % 3);
            su[a].mul(&sv[b]).sub(&su[b].mul(&sv[a]))
        })
        .collect();
    let norm = n[0]
        .square()
        .add(&n[1].square())
        .add(&n[2].square())
        .sqrt()?;
    let mut out = Vec::with_capacity(14);
    if let Some(third) = third {
        let half = T::exact_f64(0.5);
        let sq: Vec<N> = p.iter().map(|x| x.square()).collect();
        out.push(
            p[0].mul(&n[0])
                .add(&p[1].mul(&n[1]))
                .add(&p[2].mul(&n[2]))
                .scale(third),
        );
        for i in 0..3 {
            out.push(sq[i].mul(&n[i]).scale(&half));
        }
        for i in 0..3 {
            out.push(sq[i].mul(&p[i]).mul(&n[i]).scale(third));
        }
        for (i, j) in [(0, 1), (1, 2), (2, 0)] {
            out.push(sq[i].mul(&p[j]).mul(&n[i]).scale(&half));
        }
    }
    out.push(norm.clone());
    for pk in &p {
        out.push(pk.mul(&norm));
    }
    Some(out)
}

/// A pcurve piece's `(ū, v̄, ū')` in a patch's local coordinates at `τ`.
trait PieceAt<T: Real> {
    fn at<N: Num<T>>(&self, t: &N) -> Option<[N; 3]>;
}

impl<T: Real> PieceAt<T> for Piece<T> {
    fn at<N: Num<T>>(&self, t: &N) -> Option<[N; 3]> {
        Piece::at(self, t)
    }
}

/// A piece of a spline wall's meeting (S9f.2b) on its wall between two
/// knots, in the knot span's patch: the pcurve's fraction `g = ga + τ len`
/// (the edge's `g`, or `1 - g` reversed), the wall's `u` affine in it, the
/// `v` the curve's on that span's polynomial (`wall_meet::eval`).
struct WallPiece<'a> {
    m: &'a crate::topology::WallMeet,
    span: &'a super::wall_meet::Span,
    ga: R,
    len: R,
    reversed: bool,
}

impl WallPiece<'_> {
    /// The piece's local `ū` at its fraction's end `t` (0 or 1), exactly.
    fn u_local(&self, t: i64) -> Option<R> {
        let g = &self.ga + &self.len * ratio(t, 1);
        let f = if self.reversed { ratio(1, 1) - g } else { g };
        let u = R::from_float(self.m.start)? + R::from_float(self.m.sweep)? * f;
        let [u0, u1] = &self.span.u;
        Some((u - u0) / (u1 - u0))
    }
}

impl<T: Real> PieceAt<T> for WallPiece<'_> {
    fn at<N: Num<T>>(&self, t: &N) -> Option<[N; 3]> {
        let g = t.scale(&c(&self.len)).shift(&c(&self.ga));
        let f = if self.reversed {
            g.neg().shift(&T::exact_f64(1.0))
        } else {
            g
        };
        let u = f
            .scale(&T::exact_f64(self.m.sweep))
            .shift(&T::exact_f64(self.m.start));
        let (v, _) = super::wall_meet::eval(self.m, self.span, &u)?;
        let ([u0, u1], [v0, v1]) = (&self.span.u, &self.span.v);
        let ub = u.shift(&c(&-u0)).scale(&c(&(ratio(1, 1) / (u1 - u0))));
        let vb = v.shift(&c(&-v0)).scale(&c(&(ratio(1, 1) / (v1 - v0))));
        let sense = if self.reversed { -1 } else { 1 };
        let slope = R::from_float(self.m.sweep)? * &self.len * ratio(sense, 1) / (u1 - u0);
        Some([ub, vb, t.lift(&c(&slope))])
    }
}

/// Halvings of a spline wall's meeting's piece before its sweep is given up
/// (S9f.2b: its `v` a square root whose series over a wide range near a
/// turning point beyond the piece is undefined).
const WALL_SPLITS: usize = 40;

/// A spline wall's meeting's piece's sweep (`sweep_piece`), its fraction
/// range halved exactly where the rule cannot run on it whole. A piece
/// whose discriminant (binary64, at its ends and middle: a guide, not an
/// enclosure) varies by more than a factor of two is halved without a try:
/// its square root's series about the piece reaches a turning point's
/// distance, and a sweep that fails spends its whole budget first (a
/// meeting ending `10^-7` short of a turning point needs some twenty
/// halvings toward that end).
fn sweep_wall<T: Real>(
    nets: &[Net<T>],
    at: usize,
    piece: &WallPiece<'_>,
    rest: &[T; 3],
    third: Option<&T>,
    count: usize,
    depth: usize,
) -> Option<Vec<T>> {
    let (ua, ub) = (piece.u_local(0)?, piece.u_local(1)?);
    let graded = {
        let f = crate::solid::split::rational_f64;
        let (a, b) = (f(&ua), f(&ub));
        let d = [a, 0.5 * (a + b), b].map(|x| super::wall_meet::discriminant_f64(piece.span, x));
        let (lo, hi) = (d[0].min(d[1]).min(d[2]), d[0].max(d[1]).max(d[2]));
        lo > 0.5 * hi
    };
    if graded || depth >= WALL_SPLITS {
        let ends = [c(&ua), c(&ub)];
        let mut local = vec![zero::<T>(); count];
        if sweep_piece(nets, at, piece, ends, rest, third, &mut local).is_some() {
            return Some(local);
        }
    }
    if depth >= WALL_SPLITS {
        return None;
    }
    let half = &piece.len / ratio(2, 1);
    let mut out = vec![zero::<T>(); count];
    for ga in [piece.ga.clone(), &piece.ga + &half] {
        let part = WallPiece {
            m: piece.m,
            span: piece.span,
            ga,
            len: half.clone(),
            reversed: piece.reversed,
        };
        accumulate(
            &mut out,
            &sweep_wall(nets, at, &part, rest, third, count, depth + 1)?,
        );
    }
    Some(out)
}

/// `-ū'(τ) v̄(τ) f̄(ū(τ), σ v̄(τ))` of a piece in a patch's local
/// coordinates.
struct Sweep<'a, T, P> {
    net: &'a Net<T>,
    piece: &'a P,
    rest: &'a [T; 3],
    /// `1/3` when all fourteen integrands are wanted.
    third: Option<&'a T>,
}

/// A piece's sweep in patch `at` and the same `J` of the line between its
/// end `ū` values (`ends`) in every patch below it in its column, added to
/// `total`.
#[allow(clippy::too_many_arguments)]
fn sweep_piece<T: Real, P: PieceAt<T>>(
    nets: &[Net<T>],
    at: usize,
    piece: &P,
    ends: [T; 2],
    rest: &[T; 3],
    third: Option<&T>,
    total: &mut [T],
) -> Option<()> {
    let patch = &nets[at];
    accumulate(
        total,
        &integrate_2d(&Sweep {
            net: patch,
            piece,
            rest,
            third,
        })?,
    );
    let [us, vs] = &patch.domain;
    let one = || vec![T::exact_f64(1.0)];
    let line = Piece::new(ends.to_vec(), one(), one());
    for below in nets {
        let [bu, bv] = &below.domain;
        if bu == us && bv[1] <= vs[0] {
            accumulate(
                total,
                &integrate_2d(&Sweep {
                    net: below,
                    piece: &line,
                    rest,
                    third,
                })?,
            );
        }
    }
    Some(())
}

impl<T: Real, P: PieceAt<T>> Integrand2<T> for Sweep<'_, T, P> {
    fn at<N: Num<T>>(&self, tau: &N, sigma: &N) -> Option<Vec<N>> {
        let [u, v, slope] = self.piece.at(tau)?;
        let values = patch_integrands(self.net, &u, &sigma.mul(&v), self.rest, self.third)?;
        let factor = slope.mul(&v).neg();
        Some(values.iter().map(|x| x.mul(&factor)).collect())
    }
}

/// The patch of a spline pcurve piece not lying in one patch by its control
/// points (S9f.1), and per integrand a bound of the error of integrating it
/// there: a hull leaving the surface's domain while the curve certainly
/// keeps to it (its coordinate against the bound exactly nonnegative: a
/// crease nearly touching a cap) is boxed by the bound; a hull across a
/// `u` boundary by a sliver at most `2^-20` of the patch's width (a
/// crease's pcurve over a wall's knot line, its rounded identity in `u` a
/// step past it) is integrated on that patch's polynomial, the slivers'
/// error within `2 (deg + 1) δ̄` (the piece's `ū` turns at most `deg + 1`
/// times there) times `|G|` over the sliver's columns, each patch's and its
/// neighbour's (`G` the antiderivative in `v̄` from the domain's start, at
/// most the sum of `|f̄|` over the columns' boxes). `None` otherwise.
fn sliver_net<T: Real>(
    nets: &[Net<T>],
    points: &[(R, R)],
    h: &[Vec<R>; 4],
    rest: &[T; 3],
    all: bool,
) -> Option<(usize, Vec<T>)> {
    let span = |f: &dyn Fn(&(R, R)) -> &R| -> [R; 2] {
        let lo = points.iter().map(f).min().expect("a point").clone();
        let hi = points.iter().map(f).max().expect("a point").clone();
        [lo, hi]
    };
    let mut ranges = [span(&|p| &p.0), span(&|p| &p.1)];
    for (axis, range) in ranges.iter_mut().enumerate() {
        let lo = nets.iter().map(|q| &q.domain[axis][0]).min()?.clone();
        let hi = nets.iter().map(|q| &q.domain[axis][1]).max()?.clone();
        let (x, w) = (&h[axis], &h[3]);
        if range[0] < lo {
            let off: Vec<R> = x.iter().zip(w).map(|(x, w)| x - &lo * w).collect();
            if !super::bernstein::nonnegative(&off) {
                return None;
            }
            range[0] = lo;
        }
        if range[1] > hi {
            let off: Vec<R> = x.iter().zip(w).map(|(x, w)| &hi * w - x).collect();
            if !super::bernstein::nonnegative(&off) {
                return None;
            }
            range[1] = hi;
        }
    }
    let [us, vs] = &ranges;
    let count = if all { 14 } else { 4 };
    let holds_v = |q: &Net<T>| q.domain[1][0] <= vs[0] && vs[1] <= q.domain[1][1];
    let within = |q: &Net<T>| holds_v(q) && q.domain[0][0] <= us[0] && us[1] <= q.domain[0][1];
    if let Some(at) = nets.iter().position(within) {
        return Some((at, vec![zero(); count]));
    }
    let mid = (&us[0] + &us[1]) / ratio(2, 1);
    let at = nets
        .iter()
        .position(|q| holds_v(q) && q.domain[0][0] <= mid && mid <= q.domain[0][1])?;
    let [u0, u1] = &nets[at].domain[0];
    let width = u1 - u0;
    let limit = &width * ratio(1, 1 << 20);
    if u0 - &us[0] > limit || &us[1] - u1 > limit {
        return None;
    }
    let third = c::<T>(&ratio(1, 3));
    // Σ |f̄| over the u box of the patch's column (itself and the patches
    // below it), per integrand.
    let column = |q: &Net<T>, ua: &R, ub: &R| -> Option<Vec<R>> {
        let mut out = vec![ratio(0, 1); count];
        for below in nets
            .iter()
            .filter(|b| b.domain[0] == q.domain[0] && b.domain[1][0] <= q.domain[1][0])
        {
            let [b0, b1] = &below.domain[0];
            let local = |x: &R| c::<T>(&((x - b0) / (b1 - b0)));
            let u = Point(local(ua).union(&local(ub)));
            let v = Point(T::exact_f64(0.0).union(&T::exact_f64(1.0)));
            let values = patch_integrands(below, &u, &v, rest, all.then_some(&third))?;
            for (o, x) in out.iter_mut().zip(values) {
                let (lo, hi) = x.0.bounds_f64();
                let m = lo.abs().max(hi.abs());
                *o += R::from_float(m).filter(|_| m.is_finite())?;
            }
        }
        Some(out)
    };
    let degree = h[0].len().saturating_sub(1) as i64;
    let mut bound = vec![ratio(0, 1); count];
    for (outside, edge, other_edge) in [
        (u0 - &us[0], u0.clone(), 1usize),
        (&us[1] - u1, u1.clone(), 0usize),
    ] {
        if outside <= ratio(0, 1) {
            continue;
        }
        let neighbour = nets
            .iter()
            .find(|q| q.domain[0][other_edge] == edge && holds_v(q))?;
        let (lo, hi) = if other_edge == 1 {
            (us[0].clone(), edge.clone())
        } else {
            (edge.clone(), us[1].clone())
        };
        let mine = column(&nets[at], &lo, &hi)?;
        let theirs = column(neighbour, &lo, &hi)?;
        let wn = &neighbour.domain[0][1] - &neighbour.domain[0][0];
        let turns = ratio(2 * (degree + 1), 1) * &outside;
        for k in 0..count {
            bound[k] += &turns * (&mine[k] / &width + &theirs[k] / &wn);
        }
    }
    Some((
        at,
        bound.iter().map(|e| T::exact_f64(0.0).widen(e)).collect(),
    ))
}

/// The fourteen mass integrals (`all`) or the four with `|N|` of a face on a
/// nonperiodic spline surface relative to `origin`, by Green's theorem
/// through its patches with `G` from the domain's start: a boundary piece
/// in one patch contributes `J = -∫∫ ū'(τ) v̄(τ) f̄(ū(τ), σ v̄(τ)) dσ dτ` in
/// the patch's local coordinates, and every patch below it in its column
/// the same `J` of the line from the piece's end `ū` values at `v̄ = 1` (the
/// patches of a column share their `u` span, so their local `ū` agree).
/// `None` when a piece does not certainly lie in one patch or the rule
/// cannot run (the caller falls back to the strips).
pub(super) fn spline_face<T: Real>(
    surface: &BSplineSurface3,
    loops: &[Lp],
    origin: &[T; 3],
    all: bool,
) -> Option<Vec<T>> {
    if surface.u_knots().is_periodic() || surface.v_knots().is_periodic() {
        return None;
    }
    let centre: [R; 3] = std::array::from_fn(|k| origin[k].midpoint());
    let rest: [T; 3] = std::array::from_fn(|k| c::<T>(&centre[k]).sub(&origin[k]));
    let nets: Vec<Net<T>> = surface
        .bezier_patches()
        .ok()?
        .iter()
        .map(|q| net(q, &centre))
        .collect();
    let mut total = vec![zero::<T>(); if all { 14 } else { 4 }];
    let third = c::<T>(&ratio(1, 3));
    let third = all.then_some(&third);
    // The slivers' error bounds (`sliver_net`).
    let mut extras: Vec<Vec<T>> = Vec::new();
    for lp in loops {
        for fin in &lp.fins {
            // A spline wall's meeting's projection onto this wall (S9f.2b):
            // its pieces between the wall's knots, exactly, each in its
            // knot span's patch.
            if let Curve2::Projection(pr) = &fin.pcurve {
                let crate::topology::Curve3::WallMeet(m) = &pr.curve else {
                    return None;
                };
                if m.wall != *surface {
                    return None;
                }
                let spans = super::wall_meet::spans(m)?;
                let one = ratio(1, 1);
                for (fa, fb, k) in super::wall_meet::pieces(m)? {
                    let (ga, gb) = if pr.reversed {
                        (&one - &fb, &one - &fa)
                    } else {
                        (fa, fb)
                    };
                    let span = spans.spans.get(k)?;
                    let at = nets.iter().position(|q| q.domain[0] == span.u)?;
                    let piece = WallPiece {
                        m,
                        span,
                        len: &gb - &ga,
                        ga,
                        reversed: pr.reversed,
                    };
                    let values = sweep_wall(&nets, at, &piece, &rest, third, total.len(), 0)?;
                    accumulate(&mut total, &values);
                }
                continue;
            }
            let arcs = pcurve_arcs(&fin.pcurve)?;
            // Lines are split where they cross a patch boundary.
            let mut cuts = vec![ratio(0, 1), ratio(1, 1)];
            for (from, to, _) in &arcs {
                cuts.extend([from.clone(), to.clone()]);
            }
            if let Curve2::LineSegment { start, end } = &fin.pcurve {
                let (a, b) = ([r(start.x), r(start.y)], [r(end.x), r(end.y)]);
                for p in &nets {
                    for (axis, [lo, hi]) in p.domain.iter().enumerate() {
                        let d = &b[axis] - &a[axis];
                        if d == ratio(0, 1) {
                            continue;
                        }
                        for k in [lo, hi] {
                            let t = (k - &a[axis]) / &d;
                            if t > ratio(0, 1) && t < ratio(1, 1) {
                                cuts.push(t);
                            }
                        }
                    }
                }
            }
            cuts.sort();
            cuts.dedup();
            for w in cuts.windows(2) {
                let h = piece_of(&arcs, &w[0], &w[1])?;
                let points: Vec<(R, R)> = (0..h[3].len())
                    .map(|i| (&h[0][i] / &h[3][i], &h[1][i] / &h[3][i]))
                    .collect();
                // Constant u: `du = 0` along the piece.
                if points.iter().all(|p| p.0 == points[0].0) {
                    continue;
                }
                let (at, extra) = match nets.iter().position(|p| {
                    let [[u0, u1], [v0, v1]] = &p.domain;
                    points
                        .iter()
                        .all(|(u, v)| u0 <= u && u <= u1 && v0 <= v && v <= v1)
                }) {
                    Some(at) => (at, None),
                    None => {
                        let (at, extra) = sliver_net(&nets, &points, &h, &rest, all)?;
                        (at, Some(extra))
                    }
                };
                let [[u0, u1], _] = &nets[at].domain;
                let local = |u: &R| c::<T>(&((u - u0) / (u1 - u0)));
                let ends = [local(&points[0].0), local(&points[points.len() - 1].0)];
                let piece = exact_piece(&h, Some(&nets[at].domain));
                sweep_piece(&nets, at, &piece, ends, &rest, third, &mut total)?;
                extras.extend(extra);
            }
        }
        // Chords closing the loop's gaps: their ends must certainly share a
        // patch.
        for (a, b) in super::chords::<T>(lp) {
            if a[0].bounds_f64() == b[0].bounds_f64() && width(&a[0]) == 0.0 {
                continue;
            }
            let at = super::chord_patch(nets.iter().map(|p| &p.domain), &a, &b)?;
            let [[u0, u1], [v0, v1]] = &nets[at].domain;
            let local = |x: &T, lo: &R, hi: &R| x.sub(&c(lo)).div(&c(&(hi - lo)));
            let ends = [local(&a[0], u0, u1)?, local(&b[0], u0, u1)?];
            let piece = Piece::new(
                ends.to_vec(),
                vec![local(&a[1], v0, v1)?, local(&b[1], v0, v1)?],
                vec![T::exact_f64(1.0)],
            );
            sweep_piece(&nets, at, &piece, ends, &rest, third, &mut total)?;
        }
    }
    for extra in &extras {
        accumulate(&mut total, extra);
    }
    Some(total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::certified::Fast;

    fn contains<T: Real>(x: &T, want: &R) -> bool {
        let (m, w) = (x.midpoint(), x.radius());
        &m - &w <= *want && *want <= &m + &w
    }

    /// The rule integrates polynomials of degree `2n - 1` exactly: its
    /// weights sum to 2 and the moments of `x^14` and `x^15` are `2/15` and
    /// 0, within tight enclosures; `K_8 = (8!)^4 / (17 (16!)^2)`.
    #[test]
    fn the_rule_is_exact_to_degree_fifteen() {
        let rule = rule::<I>();
        let moment = |k: usize| {
            rule.nodes
                .iter()
                .zip(&rule.weights)
                .fold(I::exact(ratio(0, 1)), |s, (x, w)| {
                    let power = (0..k).fold(I::exact(ratio(1, 1)), |p, _| p.mul(x));
                    s.add(&w.mul(&power))
                })
        };
        assert!(contains(&moment(0), &ratio(2, 1)));
        assert!(contains(&moment(14), &ratio(2, 15)));
        assert!(contains(&moment(15), &ratio(0, 1)));
        assert!(moment(14).radius() < ratio(1, 1_000_000_000_000_000_000));
        let k = rule.remainder.midpoint();
        assert!(k > ratio(3551, 10_000_000_000_000) && k < ratio(3552, 10_000_000_000_000));
    }

    /// The series recurrences against known expansions at 0: `1/(1 + ε)`,
    /// `cos ε`, `sin ε` and `√(1 + ε)`.
    #[test]
    fn series_follow_their_expansions() {
        let x = Series::variable(I::exact(ratio(0, 1)), 8);
        let one = x.lift(&I::exact(ratio(1, 1)));
        let q = one.div(&one.add(&x)).unwrap();
        let (co, si) = x.cos_sin();
        for k in 0..8 {
            let sign = if k % 2 == 0 { 1 } else { -1 };
            assert!(contains(&q.coefficient(k), &ratio(sign, 1)), "1/(1+e) {k}");
            let factorial = (1..=k as i64).product::<i64>();
            let (c, s) = [(1, 0), (0, 1), (-1, 0), (0, -1)][k % 4];
            assert!(
                contains(&co.coefficient(k), &ratio(c, factorial)),
                "cos {k}"
            );
            assert!(
                contains(&si.coefficient(k), &ratio(s, factorial)),
                "sin {k}"
            );
        }
        let root = one.add(&x).sqrt().unwrap();
        let want = [
            ratio(1, 1),
            ratio(1, 2),
            ratio(-1, 8),
            ratio(1, 16),
            ratio(-5, 128),
        ];
        for (k, w) in want.iter().enumerate() {
            assert!(contains(&root.coefficient(k), w), "sqrt {k}");
        }
    }

    /// `√t` (and `√τ` on the square) is not analytic at `0`: the rule gives
    /// up there, after a bounded amount of work, rather than refining down
    /// to the depth limit along the singular side.
    struct Root;

    impl Integrand1<Fast> for Root {
        fn at<N: Num<Fast>>(&self, t: &N) -> Option<Vec<N>> {
            Some(vec![t.sqrt()?])
        }
    }

    impl Integrand2<Fast> for Root {
        fn at<N: Num<Fast>>(&self, tau: &N, sigma: &N) -> Option<Vec<N>> {
            Some(vec![tau.sqrt()?.mul(sigma)])
        }
    }

    #[test]
    fn singular_integrands_give_up() {
        let clock = std::time::Instant::now();
        assert!(integrate_1d(&Root).is_none());
        assert!(integrate_2d(&Root).is_none());
        assert!(clock.elapsed().as_secs_f64() < 5.0);
    }

    /// `∫_0^1 dt/(1 + t²) = π/4` and `∫_0^1 dt/(1 + t²)² = π/8 + 1/4`, with
    /// `1 + t²` in Bernstein form `[1, 1, 2]`, in both tiers.
    #[test]
    fn quotients_enclose_their_closed_forms() {
        fn check<T: Real>() {
            let w: Bern<T> = [1.0, 1.0, 2.0].map(T::exact_f64).to_vec();
            let m: Bern<T> = vec![T::exact_f64(1.0)];
            let pi = crate::certified::pi();
            for (k, want) in [
                (1, [pi.lo() / ratio(4, 1), pi.hi() / ratio(4, 1)]),
                (
                    2,
                    [
                        pi.lo() / ratio(8, 1) + ratio(1, 4),
                        pi.hi() / ratio(8, 1) + ratio(1, 4),
                    ],
                ),
            ] {
                let x = quotient_integral(&m, &w, k).unwrap();
                assert!(contains(&x, &want[0]) && contains(&x, &want[1]), "k = {k}");
                assert!(width(&x) < 1e-14, "k = {k}: {}", width(&x));
            }
        }
        check::<Fast>();
        check::<I>();
    }
}
