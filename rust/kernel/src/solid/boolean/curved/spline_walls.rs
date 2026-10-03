//! S9f.1: spline walls in the curved engine (REVIEW_NOTES.md, "S9f refined,
//! before its code" and "S9f.1 evidence").
//!
//! A spline profile segment (S8b's nonrational B-spline of degree `p` at
//! most 7) is its exact Bézier arcs in rationals, run the way the profile
//! runs: its run parameter `tau` is the curve's own `t`, or `first + last -
//! t` for a span reversed along its curve, so `tau` always grows along the
//! run. Its wall `o + S_x(tau) x + S_y(tau) y + w n` (`Surf::Spline`) meets
//!
//! * a plane `m . (P - p0) = 0` where `a(tau) + b w = 0`, `a = m . (o - p0)
//!   + (m . x) S_x + (m . y) S_y`, `b = m . n`: for `b != 0` in the crease
//!   `w = -a / b` (`Crv::Spline`, an affine image of the profile over the
//!   wall, S8b.3's `plane_image` once rounded), for `b = 0` in generatrices
//!   at the roots of `a` on each arc (degree `p`);
//! * a line where its trace's equation vanishes on an arc (degree `p`);
//!
//! and a cap edge or a crease (a curve `w = h0 + h1 u + h2 v` over the
//! profile) meets a plane at the roots of the plane's function along it.
//! Every root is exact: rational, or the generator `alpha` of `Q(alpha)`,
//! the arc's own parameter `s` in `[0, 1]`, shared by every root of one
//! polynomial wherever it is found (the polynomial's primitive integer
//! coefficients and the root's index key it): a point found twice, by a
//! plane's function along a cap edge and by the generatrix of that plane,
//! is one number of one field. So every vertex on a wall lies in `Q(alpha)`
//! of degree at most `p` at a known parameter, found again from its
//! coordinates (`locate`: the arc whose polynomials at the point's
//! generator give it, or for a rational point the common roots of `S_x -
//! u` and `S_y - v`).
//!
//! Membership in the prism is exact: a point on a spline segment is on the
//! profile's boundary (its side from the run's tangent there), any other
//! rational point counts the `+u` ray's crossings with each arc at the
//! roots of `S_y - v`, right of it by the exact sign of `S_x - u` there. An
//! irrational point off every spline segment (none arises against a
//! polyhedral partner: every irrational point of S9f.1 lies on a spline
//! wall or on an arc's cylinder) is undecided here; the profile classifies
//! it at a rational point of a box about it that no element of the
//! profile meets (S9f.2a, `Prism::rational_proxy`), and a point of
//! another field on a segment is found by the arc's implicit equation
//! (`foreign_param`).
//!
//! Tangencies are `Degenerate`: a plane along a generatrix (a root of even
//! multiplicity), a plane touching the wall at a knot (a root of
//! multiplicity above one on either side of it), a line touching the wall;
//! so is a plane within rounding of the wall's axis (the crease's slope past
//! `10^12`, S9c.1's (c) on the exact models: frames built from one normal
//! round it differently, `TILT` against `TILTX` leaves `n . m` at `-8.9e-17`
//! where `TILT` against `SIDE` leaves it exactly zero, generatrices).
use super::meet::{EdgeMeet, Pos, Section};
use super::model::{Affine, Crv, Prism};
use super::num::*;
use super::spline_parallel::Implicit;
use crate::polynomial::real::{AlgebraicRoot, IntPolynomial};
use crate::solid::split::spline::{power, roots, Span};
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::{Arc, OnceLock};

/// An exact Bézier arc of a spline segment, in the segment's run.
#[derive(Debug, Clone)]
pub(super) struct BArc {
    /// Control points in run order.
    pub(super) cps: Vec<[R; 2]>,
    /// The coordinates in powers of the arc's parameter `s` in `[0, 1]`.
    pub(super) x: Vec<R>,
    pub(super) y: Vec<R>,
    /// Their derivatives in `s`.
    pub(super) dx: Vec<R>,
    pub(super) dy: Vec<R>,
    /// The run parameter's range: `tau = d0 + s (d1 - d0)`.
    pub(super) d: [R; 2],
    /// The control points in binary64 (for views only).
    cf: Vec<[f64; 2]>,
    /// Its implicit equation and inversion in its frame (S9f.2a), made once.
    implicit: OnceLock<Arc<Implicit>>,
}

impl BArc {
    /// The arc's implicit equation and its parameter's inversion (S9f.2a).
    pub(super) fn implicit(&self) -> Arc<Implicit> {
        self.implicit
            .get_or_init(|| Arc::new(Implicit::of(&self.x, &self.y)))
            .clone()
    }
}

/// A spline profile segment, exact.
#[derive(Debug)]
pub(super) struct SplineSeg {
    pub(super) arcs: Vec<BArc>,
    /// The run parameter's range (the curve's domain).
    pub(super) first: R,
    pub(super) last: R,
    /// Whether the segment runs against its curve's parameter.
    pub(super) rev: bool,
    /// The stored segment with each knot of multiplicity `p` removed once
    /// where it is C1 (R4: `topology::c1_reduced_span`), the curve every
    /// result's piece is cut from, so each stays C1 rounded.
    pub(super) span: Span,
    /// Whether it bounds a hole (its material on its right).
    pub(super) hole: bool,
    /// Its prism's operand (0 for the object): two spline walls' crossings
    /// lie in the object's segment's fields (S9f.2a).
    pub(super) op: usize,
}

pub(super) fn trim(mut p: Vec<R>) -> Vec<R> {
    while p.last().is_some_and(|c| *c == zero()) {
        p.pop();
    }
    p
}

pub(super) fn derivative(p: &[R]) -> Vec<R> {
    p.iter()
        .enumerate()
        .skip(1)
        .map(|(i, c)| c * int(i as i64))
        .collect()
}

/// A polynomial's value at a number of any field (Horner).
pub(super) fn peval(p: &[R], s: &Qd) -> Qd {
    p.iter()
        .rev()
        .fold(Qd::rat(zero()), |acc, c| acc.mul(s).add_r(c))
}

pub(super) fn peval_r(p: &[R], s: &R) -> R {
    p.iter().rev().fold(zero(), |acc, c| acc * s + c)
}

/// `a + b X + c Y` of two power series.
pub(super) fn combine(a: &R, b: &R, x: &[R], c: &R, y: &[R]) -> Vec<R> {
    let n = x.len().max(y.len()).max(1);
    let mut out = vec![zero(); n];
    out[0] = a.clone();
    for (i, v) in x.iter().enumerate() {
        out[i] += b * v;
    }
    for (i, v) in y.iter().enumerate() {
        out[i] += c * v;
    }
    trim(out)
}

// ------------------------------------------------------------- generators

/// A generator kept: its polynomial's primitive coefficients, its root's
/// index and its field's generator.
type KeptGen = (Vec<BigInt>, usize, Arc<Gen>);

thread_local! {
    /// Generators by their polynomial (primitive integer coefficients, a
    /// positive leading one) and root index: every root of one polynomial is
    /// one number of one field wherever it is found.
    static GENS: std::cell::RefCell<Vec<KeptGen>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// The generators kept (the oldest dropped first).
const KEPT_GENS: usize = 4096;

/// The number of the `index`-th root of `p` in `[0, 1]`: rational exactly,
/// else the generator of its field.
fn root_value(p: &[R], root: &AlgebraicRoot, index: usize) -> K {
    if let Some(x) = root.rational_value() {
        return K::Rat(x.clone());
    }
    let ip = IntPolynomial::from_rationals(p);
    let mut key = ip.0;
    if key
        .last()
        .is_some_and(|x| x.sign() == num_bigint::Sign::Minus)
    {
        key = key.into_iter().map(|x| -x).collect();
    }
    let g = GENS.with(|g| {
        let mut g = g.borrow_mut();
        if let Some((_, _, gen)) = g.iter().find(|(k, i, _)| *i == index && *k == key) {
            return gen.clone();
        }
        let poly: Vec<R> = key.iter().map(|c| R::from_integer(c.clone())).collect();
        let gen = Arc::new(Gen::new(poly, root.clone()));
        if g.len() >= KEPT_GENS {
            g.drain(..KEPT_GENS / 2);
        }
        g.push((key, index, gen.clone()));
        gen
    });
    K::generator(&g)
}

// ------------------------------------------------------------------ roots

/// A root of a function along a segment: its run parameter, its greatest
/// multiplicity on the arcs holding it, whether it lies at an interior knot
/// or at the segment's start or end.
#[derive(Debug, Clone)]
pub(super) struct Root {
    pub(super) tau: Qd,
    pub(super) mult: usize,
    pub(super) knot: bool,
    pub(super) end: bool,
}

/// The roots of a function given on each arc, or where it vanishes on an
/// arc.
pub(super) enum Roots {
    At(Vec<Root>),
    /// Zero on every arc.
    Along,
    /// Zero on some arcs only.
    Partly,
}

/// Where a profile point lies against a spline segment.
pub(super) enum Against {
    /// On it, at this run parameter.
    On(Qd),
    /// Off it: whether it flips the `+u` ray's parity.
    Flips(bool),
    /// An irrational point off it: its ray's crossings are not decided.
    Undecided,
}

impl SplineSeg {
    pub(super) fn new(span: &Span, hole: bool, op: usize) -> Result<Self> {
        let curve = span.curve().as_curve3();
        if curve.is_rational() || curve.is_periodic() {
            return Err(Error::OutOfDomain(
                "a rational or periodic spline profile segment (S8b)",
            ));
        }
        let (a, b) = curve.domain();
        if span.range() != [a, b] {
            return Err(Error::OutOfDomain(
                "a spline profile segment over part of its curve (S8b)",
            ));
        }
        let raw = crate::solid::split::spline::arcs(span)?;
        let (first, last) = (q(a), q(b));
        let rev = span.is_reversed();
        let mirror = |t: &R| &first + &last - t;
        let mut arcs = Vec::new();
        let order: Vec<usize> = if rev {
            (0..raw.len()).rev().collect()
        } else {
            (0..raw.len()).collect()
        };
        for k in order {
            let arc = &raw[k];
            let (cps, d) = if rev {
                (
                    arc.cps.iter().rev().cloned().collect::<Vec<_>>(),
                    [mirror(&arc.domain[1]), mirror(&arc.domain[0])],
                )
            } else {
                (arc.cps.clone(), arc.domain.clone())
            };
            let x = trim(power(&cps.iter().map(|c| c[0].clone()).collect::<Vec<_>>()));
            let y = trim(power(&cps.iter().map(|c| c[1].clone()).collect::<Vec<_>>()));
            let cf = cps
                .iter()
                .map(|c| [rational_f64(&c[0]), rational_f64(&c[1])])
                .collect();
            arcs.push(BArc {
                dx: derivative(&x),
                dy: derivative(&y),
                x,
                y,
                cps,
                d,
                cf,
                implicit: OnceLock::new(),
            });
        }
        Ok(Self {
            arcs,
            first,
            last,
            rev,
            span: crate::topology::c1_reduced_span(span)?,
            hole,
            op,
        })
    }

    /// The run's start and end points.
    pub(super) fn start(&self) -> &[R; 2] {
        &self.arcs[0].cps[0]
    }

    pub(super) fn end(&self) -> &[R; 2] {
        let a = &self.arcs[self.arcs.len() - 1];
        &a.cps[a.cps.len() - 1]
    }

    /// The stored curve's own parameter of a run parameter.
    pub(super) fn t_of(&self, tau: &Qd) -> Qd {
        if self.rev {
            tau.neg().add_r(&(&self.first + &self.last))
        } else {
            tau.clone()
        }
    }

    /// The run parameter of an arc's parameter.
    pub(super) fn tau_of(&self, k: usize, s: &K) -> Qd {
        let d = &self.arcs[k].d;
        Qd::of(s.scale(&(&d[1] - &d[0]))).add_r(&d[0])
    }

    /// The arc holding a run parameter (the earlier at a knot) and its
    /// parameter there.
    fn arc_at(&self, tau: &Qd) -> (usize, Qd) {
        let n = self.arcs.len();
        let k = (0..n)
            .find(|&k| tau.cmp(&Qd::rat(self.arcs[k].d[1].clone())) != Ordering::Greater)
            .unwrap_or(n - 1);
        let d = &self.arcs[k].d;
        let s = tau.add_r(&-&d[0]).scale(&(int(1) / (&d[1] - &d[0])));
        (k, s)
    }

    /// The profile point at a run parameter.
    pub(super) fn point(&self, tau: &Qd) -> [Qd; 2] {
        let (k, s) = self.arc_at(tau);
        let a = &self.arcs[k];
        [peval(&a.x, &s), peval(&a.y, &s)]
    }

    /// The derivative in the run parameter (C1 inside: either arc's at a
    /// knot).
    pub(super) fn deriv(&self, tau: &Qd) -> [Qd; 2] {
        let (k, s) = self.arc_at(tau);
        let a = &self.arcs[k];
        let inv = int(1) / (&a.d[1] - &a.d[0]);
        [peval(&a.dx, &s).scale(&inv), peval(&a.dy, &s).scale(&inv)]
    }

    /// The run's direction at its start or end (its first or last nonzero
    /// control leg).
    pub(super) fn end_tangent(&self, at_start: bool) -> [R; 2] {
        let mut legs = self
            .arcs
            .iter()
            .flat_map(|a| a.cps.windows(2))
            .map(|w| [&w[1][0] - &w[0][0], &w[1][1] - &w[0][1]])
            .filter(|d| d[0] != zero() || d[1] != zero());
        if at_start {
            legs.next().expect("a nonzero control leg")
        } else {
            legs.next_back().expect("a nonzero control leg")
        }
    }

    /// The roots along the segment of a function given on each arc.
    pub(super) fn roots_of(&self, poly: &dyn Fn(&BArc) -> Vec<R>) -> Result<Roots> {
        let n = self.arcs.len();
        let mut out: Vec<Root> = Vec::new();
        let mut zero_arcs = 0;
        let one = int(1);
        for (k, arc) in self.arcs.iter().enumerate() {
            let p = trim(poly(arc));
            if p.is_empty() {
                zero_arcs += 1;
                continue;
            }
            for (i, root) in roots(&p)?.iter().enumerate() {
                let s = root_value(&p, root, i);
                let at0 = matches!(&s, K::Rat(x) if *x == zero());
                let at1 = matches!(&s, K::Rat(x) if *x == one);
                let tau = self.tau_of(k, &s);
                let mult = root.multiplicity();
                if at0 && k > 0 {
                    // The previous arc's end: one root at the knot.
                    if let Some(last) = out.last_mut() {
                        if last.knot && last.tau.cmp(&tau) == Ordering::Equal {
                            last.mult = last.mult.max(mult);
                            continue;
                        }
                    }
                    out.push(Root {
                        tau,
                        mult,
                        knot: true,
                        end: false,
                    });
                    continue;
                }
                out.push(Root {
                    tau,
                    mult,
                    knot: at1 && k + 1 < n,
                    end: (at0 && k == 0) || (at1 && k + 1 == n),
                });
            }
        }
        Ok(match zero_arcs {
            0 => Roots::At(out),
            z if z == n => Roots::Along,
            _ => Roots::Partly,
        })
    }

    /// The run parameter of a profile point on the segment, if it lies on
    /// it: for a point of `Q(alpha)` the arc whose polynomials at `alpha`
    /// give it (every such point was found at its arc's parameter), for a
    /// rational point the common roots of `S_x - u` and `S_y - v`; a surd's
    /// point is never on a spline (none but a joint's, rational). S9f.2a: a
    /// point of another field (a crossing found on the other input's spline
    /// or cylinder) by each arc's implicit equation and its parameter's
    /// inversion (`foreign_param`).
    pub(super) fn locate(&self, x: &[Qd; 2]) -> Option<Qd> {
        if x.iter().any(|c| c.field().is_some()) {
            return None;
        }
        if let Some(g) = x[0].gen().or(x[1].gen()).cloned() {
            let alpha = Qd::of(K::generator(&g));
            for (k, a) in self.arcs.iter().enumerate() {
                if x[0].cmp(&peval(&a.x, &alpha)) == Ordering::Equal
                    && x[1].cmp(&peval(&a.y, &alpha)) == Ordering::Equal
                {
                    return Some(self.tau_of(k, &alpha.a));
                }
            }
            return self.foreign_param(x);
        }
        let (u, v) = (x[0].rational()?, x[1].rational()?);
        self.rational_param(&[u.clone(), v.clone()])
    }

    /// The run parameter of a point of an algebraic field on the segment
    /// whose generator is not its arcs' own (S9f.2a): on each arc whose
    /// control box holds the point's enclosure, the exact sign of the arc's
    /// implicit equation there and, where it vanishes, the arc's parameter
    /// from its inversion (`Implicit::param`), the arc's point there the
    /// given one exactly.
    fn foreign_param(&self, x: &[Qd; 2]) -> Option<Qd> {
        let (ix, iy) = (x[0].interval(), x[1].interval());
        let one = Qd::rat(int(1));
        let nil = Qd::rat(zero());
        for (k, a) in self.arcs.iter().enumerate() {
            let apart = |i: usize, e: &crate::certified::Interval| {
                a.cps.iter().all(|c| &c[i] < e.lo()) || a.cps.iter().all(|c| &c[i] > e.hi())
            };
            if apart(0, &ix) || apart(1, &iy) {
                continue;
            }
            let imp = a.implicit();
            if imp.value(x).sign() != Ordering::Equal {
                continue;
            }
            let Some(s) = imp.param(x) else {
                continue;
            };
            if s.cmp(&nil) == Ordering::Less || s.cmp(&one) == Ordering::Greater {
                continue;
            }
            if x[0].cmp(&peval(&a.x, &s)) != Ordering::Equal
                || x[1].cmp(&peval(&a.y, &s)) != Ordering::Equal
            {
                continue;
            }
            return Some(self.tau_of(k, &s.a));
        }
        None
    }

    /// A rational point's run parameter on the segment.
    fn rational_param(&self, x: &[R; 2]) -> Option<Qd> {
        for (k, a) in self.arcs.iter().enumerate() {
            // The control hull's box first.
            let outside =
                |i: usize| a.cps.iter().all(|c| c[i] < x[i]) || a.cps.iter().all(|c| c[i] > x[i]);
            if outside(0) || outside(1) {
                continue;
            }
            let px = combine(&-&x[0], &int(1), &a.x, &zero(), &[]);
            let py = combine(&-&x[1], &int(1), &a.y, &zero(), &[]);
            let common = match (px.is_empty(), py.is_empty()) {
                // A constant arc (none: the profile's arcs are not points).
                (true, true) => return Some(Qd::rat(a.d[0].clone())),
                (true, false) => IntPolynomial::from_rationals(&py),
                (false, true) => IntPolynomial::from_rationals(&px),
                (false, false) => {
                    IntPolynomial::from_rationals(&px).gcd(&IntPolynomial::from_rationals(&py))
                }
            };
            if common.is_constant() {
                continue;
            }
            let c: Vec<R> = common
                .0
                .iter()
                .map(|c| R::from_integer(c.clone()))
                .collect();
            let rs = roots(&c).ok()?;
            if let Some((i, r)) = rs.iter().enumerate().next() {
                return Some(self.tau_of(k, &root_value(&c, r, i)));
            }
        }
        None
    }

    /// A profile point against the segment: on it, or whether it flips the
    /// parity of the `+u` ray from it (rational points only).
    pub(super) fn against(&self, x: &[Qd; 2]) -> Result<Against> {
        if let Some(tau) = self.locate(x) {
            return Ok(Against::On(tau));
        }
        let (Some(u), Some(v)) = (x[0].rational(), x[1].rational()) else {
            return Ok(Against::Undecided);
        };
        Ok(Against::Flips(self.ray(&[u.clone(), v.clone()])? % 2 == 1))
    }

    /// How often the segment's points turn between above the point's `v`
    /// (strictly) and not above it right of the point, along the `+u` ray
    /// from a rational point off the segment: a crossing's parity, by the
    /// predicate the chords of lines and arcs count at their ends (a joint
    /// at the point's height not above it on either side).
    fn ray(&self, x: &[R; 2]) -> Result<usize> {
        let mut count = 0usize;
        let one = int(1);
        for a in &self.arcs {
            // Entirely above, below or left of the point: no turn right of
            // it (the control hull holds the arc).
            if a.cps.iter().all(|c| c[1] > x[1])
                || a.cps.iter().all(|c| c[1] < x[1])
                || a.cps.iter().all(|c| c[0] < x[0])
            {
                continue;
            }
            let py = combine(&-&x[1], &int(1), &a.y, &zero(), &[]);
            let px = combine(&-&x[0], &int(1), &a.x, &zero(), &[]);
            if py.is_empty() {
                // At the point's height throughout: never above.
                continue;
            }
            let ipx = IntPolynomial::from_rationals(&px);
            let right_at = |s: &R| peval_r(&px, s) > zero();
            // Where a status changes: a rational parameter, or a root.
            #[derive(Clone)]
            enum At {
                R(R),
                Root(Box<AlgebraicRoot>),
            }
            let above_at = |s: &R| peval_r(&py, s) > zero();
            let mut seq: Vec<(bool, Option<At>)> = vec![(above_at(&zero()), Some(At::R(zero())))];
            // Each gap's status at a rational between its roots' isolators
            // (no root of `py` lies there), each interior root's not above.
            let mut last_hi = zero();
            for r in &roots(&py)? {
                let (lo, hi) = (r.isolator().0.clone(), r.isolator().1.clone());
                if lo == hi && (lo == zero() || lo == one) {
                    continue;
                }
                seq.push((above_at(&((&last_hi + &lo) / int(2))), None));
                seq.push((false, Some(At::Root(Box::new(r.clone())))));
                last_hi = hi;
            }
            seq.push((above_at(&((&last_hi + &one) / int(2))), None));
            seq.push((above_at(&one), Some(At::R(one.clone()))));
            for w in seq.windows(2) {
                if w[0].0 == w[1].0 {
                    continue;
                }
                let at = w[1].1.as_ref().or(w[0].1.as_ref()).expect("a place");
                let right = match at {
                    At::R(s) => right_at(s),
                    At::Root(r) => r.sign_polynomial(&ipx) == Ordering::Greater,
                };
                if right {
                    count += 1;
                }
            }
        }
        Ok(count)
    }

    /// A binary64 view of a run parameter's point.
    pub(super) fn point_f64(&self, tau: f64) -> [f64; 2] {
        let n = self.arcs.len();
        let k = (0..n)
            .find(|&k| tau <= rational_f64(&self.arcs[k].d[1]))
            .unwrap_or(n - 1);
        let a = &self.arcs[k];
        let (d0, d1) = (rational_f64(&a.d[0]), rational_f64(&a.d[1]));
        let s = ((tau - d0) / (d1 - d0)).clamp(0.0, 1.0);
        casteljau(&a.cf, s)
    }

    /// The run parameter of the segment's point nearest a profile point
    /// (binary64: a coarse search, then Newton steps; a view only).
    pub(super) fn tau_near(&self, u: f64, v: f64) -> f64 {
        let mut best = (f64::INFINITY, 0usize, 0.0f64);
        for (k, a) in self.arcs.iter().enumerate() {
            for i in 0..=32 {
                let s = f64::from(i) / 32.0;
                let p = casteljau(&a.cf, s);
                let d = (p[0] - u).powi(2) + (p[1] - v).powi(2);
                if d < best.0 {
                    best = (d, k, s);
                }
            }
        }
        let (_, k, mut s) = best;
        let a = &self.arcs[k];
        for _ in 0..24 {
            let (p, d1, d2) = jets(&a.cf, s);
            let e = [p[0] - u, p[1] - v];
            let g = e[0] * d1[0] + e[1] * d1[1];
            let h = d1[0] * d1[0] + d1[1] * d1[1] + e[0] * d2[0] + e[1] * d2[1];
            if h.abs() < f64::MIN_POSITIVE {
                break;
            }
            let next = (s - g / h).clamp(0.0, 1.0);
            if next == s {
                break;
            }
            s = next;
        }
        let (d0, d1) = (rational_f64(&a.d[0]), rational_f64(&a.d[1]));
        d0 + s * (d1 - d0)
    }

    /// The control points of every arc (their hull holds the segment).
    pub(super) fn hull(&self) -> impl Iterator<Item = &[R; 2]> {
        self.arcs.iter().flat_map(|a| a.cps.iter())
    }
}

/// A point of a Bézier arc (binary64).
fn casteljau(cps: &[[f64; 2]], s: f64) -> [f64; 2] {
    let mut w: Vec<[f64; 2]> = cps.to_vec();
    let n = w.len();
    for r in 1..n {
        for i in 0..n - r {
            w[i] = [
                w[i][0] * (1.0 - s) + w[i + 1][0] * s,
                w[i][1] * (1.0 - s) + w[i + 1][1] * s,
            ];
        }
    }
    w[0]
}

/// A Bézier arc's point and first two derivatives (binary64).
fn jets(cps: &[[f64; 2]], s: f64) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let p = cps.len() - 1;
    let diff = |c: &[[f64; 2]]| -> Vec<[f64; 2]> {
        let m = (c.len() - 1) as f64;
        c.windows(2)
            .map(|w| [m * (w[1][0] - w[0][0]), m * (w[1][1] - w[0][1])])
            .collect()
    };
    let d1 = if p >= 1 { diff(cps) } else { vec![[0.0, 0.0]] };
    let d2 = if p >= 2 { diff(&d1) } else { vec![[0.0, 0.0]] };
    (casteljau(cps, s), casteljau(&d1, s), casteljau(&d2, s))
}

// ------------------------------------------------------------- the curves

/// A curve over a spline segment on its wall: `w = h0 + h1 u + h2 v` above
/// the profile point `(u, v) = S(tau)`, on the model's frame; a cap edge at
/// a constant height or a plane's crease. Placed by the run parameter.
#[derive(Debug, Clone)]
pub(super) struct WallCrv {
    pub(super) seg: Arc<SplineSeg>,
    pub(super) f: Affine,
    pub(super) h: [R; 3],
    /// The model's stored frame (the rounded curves are built on it).
    pub(super) frame: crate::Frame3,
}

impl PartialEq for WallCrv {
    fn eq(&self, o: &Self) -> bool {
        Arc::ptr_eq(&self.seg, &o.seg) && self.h == o.h && self.f == o.f
    }
}

/// A model point `o + u x + v y + w n` of any field.
pub(super) fn qpoint(f: &Affine, u: &Qd, v: &Qd, w: &Qd) -> QV {
    qadd(
        &qadd(&qadd(&qv(&f.o), &qscale(&f.x, u)), &qscale(&f.y, v)),
        &qscale(&f.n, w),
    )
}

impl WallCrv {
    fn height(&self, uv: &[Qd; 2]) -> Qd {
        uv[0]
            .scale(&self.h[1])
            .add(&uv[1].scale(&self.h[2]))
            .add_r(&self.h[0])
    }

    /// The height's exact value at a rational profile point, rounded once.
    pub(super) fn height_f64(&self, p: crate::Point2) -> f64 {
        rational_f64(&(&self.h[0] + &self.h[1] * q(p.x) + &self.h[2] * q(p.y)))
    }

    pub(super) fn point(&self, tau: &Qd) -> QV {
        let uv = self.seg.point(tau);
        let w = self.height(&uv);
        qpoint(&self.f, &uv[0], &uv[1], &w)
    }

    /// The tangent along increasing `tau` (unit-free).
    pub(super) fn tangent(&self, tau: &Qd) -> QV {
        let d = self.seg.deriv(tau);
        let dw = d[0].scale(&self.h[1]).add(&d[1].scale(&self.h[2]));
        qadd(
            &qadd(&qscale(&self.f.x, &d[0]), &qscale(&self.f.y, &d[1])),
            &qscale(&self.f.n, &dw),
        )
    }

    /// The run parameter of a point on the curve, if it is on it.
    pub(super) fn place(&self, x: &QV) -> Option<Qd> {
        let l = self.f.local_q(x);
        let uv = [l[0].clone(), l[1].clone()];
        let tau = self.seg.locate(&uv)?;
        (l[2].cmp(&self.height(&uv)) == Ordering::Equal).then_some(tau)
    }

    pub(super) fn on(&self, x: &QV) -> bool {
        self.place(x).is_some()
    }

    /// Whether it lies at one height (a cap edge's curve).
    pub(super) fn level(&self) -> bool {
        self.h[1] == zero() && self.h[2] == zero()
    }

    /// Points in binary64 from `t0` to `t1` (run parameters).
    pub(super) fn samples(&self, t0: f64, t1: f64, n: usize) -> Vec<[f64; 3]> {
        let fl = |v: &V| v.clone().map(|x| rational_f64(&x));
        let (o, x, y, nn) = (fl(&self.f.o), fl(&self.f.x), fl(&self.f.y), fl(&self.f.n));
        let h = self.h.clone().map(|x| rational_f64(&x));
        (0..=n)
            .map(|i| {
                let t = t0 + (t1 - t0) * i as f64 / n as f64;
                let [u, v] = self.seg.point_f64(t);
                let w = h[0] + h[1] * u + h[2] * v;
                [0, 1, 2].map(|k| o[k] + x[k] * u + y[k] * v + nn[k] * w)
            })
            .collect()
    }
}

// ------------------------------------------------------------ the meetings

fn tangent_plane() -> Error {
    Error::Degenerate("a plane tangent to a spline wall along a generatrix")
}

fn tangent_knot() -> Error {
    Error::Degenerate("a plane tangent to a spline wall along the generatrix of a knot")
}

fn tangent_edge() -> Error {
    Error::Degenerate("an edge of one input tangent to a spline wall of the other")
}

fn tangent_curve() -> Error {
    Error::Degenerate("a spline edge of one input tangent to a face of the other")
}

/// A plane's section of a spline wall of the model: the crease, or the
/// generatrices strictly inside the segment's run.
pub(super) fn plane_wall(model: &Prism, seg: &Arc<SplineSeg>, p: &V, m: &V) -> Result<Section> {
    let f = &model.f;
    let (mx, my, mn) = (dot(m, &f.x), dot(m, &f.y), dot(m, &f.n));
    let k = dot(m, &sub(&f.o, p));
    if mn != zero() {
        // A plane within rounding of the wall's axis: the crease's slope
        // past 10^12 (S9c.1's (c)).
        if &mx * &mx + &my * &my > int(10).pow(24) * &mn * &mn {
            return Err(Error::Degenerate(
                "a plane within rounding of a spline wall's axis",
            ));
        }
        let inv = int(-1) / &mn;
        return Ok(Section::Curves(vec![Crv::Spline(Box::new(WallCrv {
            seg: seg.clone(),
            f: f.clone(),
            h: [&k * &inv, &mx * &inv, &my * &inv],
            frame: model.frame,
        }))]));
    }
    // Parallel to the axis: generatrices at the roots of `a` inside the run.
    let rs = match seg.roots_of(&|a: &BArc| combine(&k, &mx, &a.x, &my, &a.y))? {
        Roots::At(r) => r,
        Roots::Along | Roots::Partly => {
            return Err(Error::OutOfDomain("a spline wall along a plane (S9f)"))
        }
    };
    let mut out = Vec::new();
    for r in rs {
        if r.knot && r.mult > 1 {
            return Err(tangent_knot());
        }
        if r.mult > 1 {
            return Err(tangent_plane());
        }
        if r.end {
            // The wall's vertical edge, on the plane: met there.
            continue;
        }
        let uv = seg.point(&r.tau);
        out.push(Crv::Line {
            p: qpoint(f, &uv[0], &uv[1], &Qd::rat(zero())),
            d: f.n.clone(),
        });
    }
    Ok(Section::Curves(out))
}

/// Where a line `p + t d` meets a spline wall of the model on frame `f`
/// (a model edge: `p` rational).
pub(super) fn line_wall(p: &QV, d: &V, f: &Affine, seg: &SplineSeg) -> Result<EdgeMeet> {
    let l = f.local_q(p);
    let ld = f.local_dir(d);
    let (Some(u0), Some(v0)) = (l[0].rational(), l[1].rational()) else {
        return Err(Error::ComputationLimit(
            "an irrational line against a spline wall",
        ));
    };
    let (u0, v0) = (u0.clone(), v0.clone());
    if ld[0] == zero() && ld[1] == zero() {
        // Along the wall's axis: on it or apart.
        return Ok(if seg.rational_param(&[u0, v0]).is_some() {
            EdgeMeet::Along
        } else {
            EdgeMeet::None
        });
    }
    // The trace's equation: dv (u - u0) - du (v - v0).
    let c = &ld[0] * &v0 - &ld[1] * &u0;
    let (a, b) = (ld[1].clone(), -&ld[0]);
    let rs = match seg.roots_of(&|arc: &BArc| combine(&c, &a, &arc.x, &b, &arc.y))? {
        Roots::At(r) => r,
        Roots::Along => return Ok(EdgeMeet::Along),
        Roots::Partly => {
            return Err(Error::Degenerate(
                "an edge of one input along a spline wall's straight arc",
            ))
        }
    };
    let dd = &ld[0] * &ld[0] + &ld[1] * &ld[1];
    let mut out = Vec::new();
    for r in rs {
        if r.mult > 1 {
            return Err(if r.knot {
                tangent_knot()
            } else {
                tangent_edge()
            });
        }
        let uv = seg.point(&r.tau);
        // Along the line: ((u, v) - (u0, v0)) . (du, dv) / |(du, dv)|^2.
        let t = uv[0]
            .add_r(&-&u0)
            .scale(&ld[0])
            .add(&uv[1].add_r(&-&v0).scale(&ld[1]))
            .scale(&(int(1) / &dd));
        let x = qadd(p, &qscale(d, &t));
        out.push((Pos::T(t), x));
    }
    Ok(EdgeMeet::Points(out))
}

/// Where a curve over a spline segment meets a plane: its points and run
/// parameters.
pub(super) fn wallcrv_plane(c: &WallCrv, p0: &V, m: &V) -> Result<EdgeMeet> {
    let f = &c.f;
    let mn = dot(m, &f.n);
    let c0 = dot(m, &sub(&f.o, p0)) + &mn * &c.h[0];
    let cu = dot(m, &f.x) + &mn * &c.h[1];
    let cv = dot(m, &f.y) + &mn * &c.h[2];
    let rs = match c
        .seg
        .roots_of(&|a: &BArc| combine(&c0, &cu, &a.x, &cv, &a.y))?
    {
        Roots::At(r) => r,
        Roots::Along => return Ok(EdgeMeet::Along),
        Roots::Partly => {
            return Err(Error::Degenerate(
                "an edge of one input on a face of the other",
            ))
        }
    };
    let mut out = Vec::new();
    for r in rs {
        if r.mult > 1 {
            return Err(if r.knot {
                tangent_knot()
            } else {
                tangent_curve()
            });
        }
        let x = c.point(&r.tau);
        out.push((Pos::T(r.tau), x));
    }
    Ok(EdgeMeet::Points(out))
}
