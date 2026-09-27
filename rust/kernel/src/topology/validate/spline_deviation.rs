//! Certified deviation of spline uses by exact composition (S4b of
//! REVIEW_NOTES.md). Where every factor is rational (a line or spline edge,
//! a line or spline pcurve, a plane or spline surface), the use's deviation
//! `D(t) = C(s(t)) - S(P(t))` is a rational function on each common Bézier
//! piece. In Bernstein form over each piece the kernel forms the
//! homogeneous numerator `N = C_xyz Ŝ_w - Ŝ_xyz C_w` and the positive
//! denominator `W = C_w Ŝ_w`, where `Ŝ = S_h(P) P_w^(du+dv)` clears the
//! pcurve's weight. Then `|D| <= max_i |N_i| / min_j W_j` whenever every
//! `W_j > 0`: cancellation is exact, since the pieces are extracted exactly
//! and composed in binary64 intervals, then exactly when those cannot
//! decide. A piece is halved until the bound holds or a depth is reached;
//! a failure is certified by a sample.
use super::bernstein::{
    c, curve_arcs, difference, lift, pcurve_arcs, piece_of, power, product, r, ratio, scaled, sum,
    value, Arcs, Bern,
};
use super::{tiered, Verdict};
use crate::certified::{Fast, Interval as I, Real};
use crate::surface::ExactBezierSurface3;
use crate::topology::{Curve2, Curve3, Surface};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// The composed degree above which the exact path declines (the Taylor path
/// then decides: high-degree products are slow and wide in intervals).
const MAX_DEGREE: usize = 24;
/// Halvings of one piece before the bound is left undecided.
const MAX_DEPTH: usize = 10;

/// One common piece: the edge's and the pcurve's homogeneous Bernstein
/// coordinates, and the surface patch under the pcurve (none on a plane).
struct Piece {
    edge: [Vec<R>; 4],
    uv: [Vec<R>; 4],
    patch: Option<ExactBezierSurface3>,
}

/// A use whose factors are all rational, piece by piece.
pub(super) struct Exact {
    plane: Option<[[R; 3]; 3]>,
    pieces: Vec<Piece>,
}

/// The numerator and denominator of `D` on one piece.
#[derive(Clone)]
struct Quotient<T> {
    n: [Bern<T>; 3],
    w: Bern<T>,
}

impl<T: Real> Quotient<T> {
    fn halves(&self) -> (Self, Self) {
        let split = |a: &Bern<T>| super::bernstein::halves(a);
        let n = self.n.clone().map(|x| split(&x));
        let w = split(&self.w);
        (
            Self {
                n: [n[0].0.clone(), n[1].0.clone(), n[2].0.clone()],
                w: w.0,
            },
            Self {
                n: [n[0].1.clone(), n[1].1.clone(), n[2].1.clone()],
                w: w.1,
            },
        )
    }
    /// `|N(t)|^2 > tol^2 W(t)^2` certainly.
    fn beyond_at(&self, t: &T, tol2: &T) -> bool {
        let w = value(&self.w, t);
        let n2 = self
            .n
            .iter()
            .fold(T::exact_f64(0.0), |s, x| s.add(&value(x, t).square()));
        n2.cmp(&tol2.mul(&w.square())) == Some(Ordering::Greater)
    }
    /// Every `W_j > 0` and `|N_i|^2 <= tol^2 W_j^2` for all `i`, `j`.
    fn within(&self, tol2: &T) -> bool {
        if self.w.iter().any(|w| w.sign() != Some(Ordering::Greater)) {
            return false;
        }
        let n2: Vec<T> = (0..self.n[0].len())
            .map(|i| {
                self.n
                    .iter()
                    .fold(T::exact_f64(0.0), |s, x| s.add(&x[i].square()))
            })
            .collect();
        self.w.iter().all(|w| {
            let limit = tol2.mul(&w.square());
            n2.iter()
                .all(|n| matches!(n.cmp(&limit), Some(Ordering::Less | Ordering::Equal)))
        })
    }
}

impl Exact {
    fn quotient<T: Real>(&self, piece: &Piece) -> Quotient<T> {
        let edge: [Bern<T>; 4] = std::array::from_fn(|k| lift(&piece.edge[k]));
        let uv: [Bern<T>; 4] = std::array::from_fn(|k| lift(&piece.uv[k]));
        let s = match (&self.plane, &piece.patch) {
            (Some(frame), _) => plane_composed(frame, &uv),
            (None, Some(patch)) => patch_composed(patch, &uv),
            (None, None) => unreachable!("a surface for every piece"),
        };
        Quotient {
            n: std::array::from_fn(|k| {
                difference(&product(&edge[k], &s[3]), &product(&s[k], &edge[3]))
            }),
            w: product(&edge[3], &s[3]),
        }
    }

    fn decide_in<T: Real>(&self, tol: f64) -> Verdict {
        let tol2 = T::exact_f64(tol).square();
        let samples: Vec<T> = (0..=4).map(|k| c(&ratio(k, 4))).collect();
        let mut stack = Vec::new();
        for piece in &self.pieces {
            let q = self.quotient::<T>(piece);
            if samples.iter().any(|t| q.beyond_at(t, &tol2)) {
                return Verdict::Beyond;
            }
            stack.push((q, 0));
        }
        let half = T::exact_f64(0.5);
        let mut undecided = false;
        while let Some((q, depth)) = stack.pop() {
            if q.within(&tol2) {
                continue;
            }
            if q.beyond_at(&half, &tol2) {
                return Verdict::Beyond;
            }
            if depth == MAX_DEPTH {
                undecided = true;
                continue;
            }
            let (a, b) = q.halves();
            stack.push((a, depth + 1));
            stack.push((b, depth + 1));
        }
        if undecided {
            Verdict::Unknown
        } else {
            Verdict::Within
        }
    }

    /// A certified upper bound of `|D|` over the use (binary64 intervals,
    /// rounded outward), after up to four halvings of a piece whose
    /// denominator is not yet certainly positive; `None` otherwise.
    pub(super) fn upper_bound(&self) -> Option<f64> {
        let mut worst = 0.0_f64;
        for piece in &self.pieces {
            let mut stack = vec![(self.quotient::<Fast>(piece), 0)];
            while let Some((q, depth)) = stack.pop() {
                if q.w.iter().any(|w| w.sign() != Some(Ordering::Greater)) {
                    if depth == 4 {
                        return None;
                    }
                    let (a, b) = q.halves();
                    stack.push((a, depth + 1));
                    stack.push((b, depth + 1));
                    continue;
                }
                let n = (0..q.n[0].len())
                    .map(|i| {
                        q.n.iter()
                            .fold(Fast::exact_f64(0.0), |s, x| s.add(&x[i].square()))
                            .bounds_f64()
                            .1
                    })
                    .fold(0.0_f64, f64::max);
                let w =
                    q.w.iter()
                        .map(|w| w.bounds_f64().0)
                        .fold(f64::INFINITY, f64::min);
                let bound = Fast::exact_f64(n)
                    .sqrt()
                    .div(&Fast::exact_f64(w))?
                    .bounds_f64()
                    .1;
                if !bound.is_finite() {
                    return None;
                }
                worst = worst.max(bound);
            }
        }
        Some(worst)
    }

    /// Within `tol`, certainly beyond it, or undecided.
    pub(super) fn decide(&self, tol: f64) -> Verdict {
        tiered(
            Verdict::Unknown,
            || self.decide_in::<Fast>(tol),
            || self.decide_in::<I>(tol),
        )
    }
}

fn plane_composed<T: Real>(frame: &[[R; 3]; 3], p: &[Bern<T>; 4]) -> [Bern<T>; 4] {
    let (pu, pv, pw) = (&p[0], &p[1], &p[3]);
    let [o, x, y] = frame;
    let xyz: [Bern<T>; 3] = std::array::from_fn(|k| {
        sum(
            &sum(&scaled(pw, &c(&o[k])), &scaled(pu, &c(&x[k]))),
            &scaled(pv, &c(&y[k])),
        )
    });
    let [a, b, cc] = xyz;
    [a, b, cc, pw.clone()]
}

/// `S_h(P) P_w^(du+dv)` on a Bézier patch, with the local coordinates
/// `U = (P_u - u0 P_w)/(u1 - u0)` and `V` likewise.
fn patch_composed<T: Real>(patch: &ExactBezierSurface3, p: &[Bern<T>; 4]) -> [Bern<T>; 4] {
    let (pu, pv, pw) = (&p[0], &p[1], &p[3]);
    let [du, dv] = patch.degrees();
    let [[u0, u1], [v0, v1]] = patch.domain().clone();
    let local = |x: &Bern<T>, lo: &R, hi: &R| {
        scaled(
            &difference(x, &scaled(pw, &c(lo))),
            &c(&(ratio(1, 1) / (hi - lo))),
        )
    };
    let (big_u, big_v) = (local(pu, &u0, &u1), local(pv, &v0, &v1));
    let (rest_u, rest_v) = (difference(pw, &big_u), difference(pw, &big_v));
    let basis = |a: &Bern<T>, b: &Bern<T>, n: usize, i: usize| {
        let k = (0..i).fold(ratio(1, 1), |acc, j| {
            acc * ratio((n - j) as i64, (j + 1) as i64)
        });
        scaled(&product(&power(a, i), &power(b, n - i)), &c(&k))
    };
    let bu: Vec<Bern<T>> = (0..=du).map(|i| basis(&big_u, &rest_u, du, i)).collect();
    let bv: Vec<Bern<T>> = (0..=dv).map(|j| basis(&big_v, &rest_v, dv, j)).collect();
    let poles = patch.homogeneous_poles();
    let mut out: [Bern<T>; 4] = std::array::from_fn(|_| vec![T::exact_f64(0.0)]);
    for i in 0..=du {
        for j in 0..=dv {
            let b = product(&bu[i], &bv[j]);
            let q = &poles[i * (dv + 1) + j];
            for k in 0..4 {
                out[k] = sum(&out[k], &scaled(&b, &c(&q[k])));
            }
        }
    }
    out
}

/// Fractions where a line pcurve crosses a patch's boundary.
fn crossings(p: &Curve2, patches: &[ExactBezierSurface3]) -> Vec<R> {
    let Curve2::LineSegment { start, end } = p else {
        return Vec::new();
    };
    let (a, b) = ([r(start.x), r(start.y)], [r(end.x), r(end.y)]);
    let (zero, one) = (ratio(0, 1), ratio(1, 1));
    let mut out = Vec::new();
    for q in patches {
        for (axis, [lo, hi]) in q.domain().iter().enumerate() {
            let d = &b[axis] - &a[axis];
            if d == zero {
                continue;
            }
            for k in [lo, hi] {
                let t = (k - &a[axis]) / &d;
                if t > zero && t < one {
                    out.push(t);
                }
            }
        }
    }
    out
}

/// The patch whose closed box holds every control point of a pcurve piece,
/// and so the whole piece by its convex hull.
fn patch_under<'a>(
    patches: &'a [ExactBezierSurface3],
    uv: &[Vec<R>; 4],
) -> Option<&'a ExactBezierSurface3> {
    let points: Vec<(R, R)> = (0..uv[3].len())
        .map(|i| (&uv[0][i] / &uv[3][i], &uv[1][i] / &uv[3][i]))
        .collect();
    patches.iter().find(|q| {
        let [[u0, u1], [v0, v1]] = q.domain();
        points
            .iter()
            .all(|(u, v)| u0 <= u && u <= u1 && v0 <= v && v <= v1)
    })
}

/// The pieces of a use whose factors are all rational, or `None` (an arc,
/// an analytic curved surface, a spline pcurve across a surface's knot
/// line, or too high a composed degree).
pub(super) fn rational_use(
    curve: &Curve3,
    surface: &Surface,
    pcurve: &Curve2,
    forward: bool,
) -> Option<Exact> {
    let edge: Arcs = curve_arcs(curve)?;
    let uv: Arcs = pcurve_arcs(pcurve)?;
    let (plane, patches) = match surface {
        Surface::Plane(f) => (
            Some([
                f.origin().to_array().map(r),
                f.x().to_array().map(r),
                f.y().to_array().map(r),
            ]),
            Vec::new(),
        ),
        Surface::BSpline(s) => (None, s.bezier_patches().ok()?),
        _ => return None,
    };
    let (zero, one) = (ratio(0, 1), ratio(1, 1));
    let flip = |s: &R| &one - s;
    let mut cuts: Vec<R> = vec![zero.clone(), one.clone()];
    for (from, to, _) in &uv {
        cuts.extend([from.clone(), to.clone()]);
    }
    for (from, to, _) in &edge {
        for s in [from, to] {
            cuts.push(if forward { s.clone() } else { flip(s) });
        }
    }
    cuts.extend(crossings(pcurve, &patches));
    cuts.sort();
    cuts.dedup();
    let mut pieces = Vec::new();
    for w in cuts.windows(2) {
        let (t0, t1) = (&w[0], &w[1]);
        let p = piece_of(&uv, t0, t1)?;
        let e = if forward {
            piece_of(&edge, t0, t1)?
        } else {
            let mut e = piece_of(&edge, &flip(t1), &flip(t0))?;
            e.iter_mut().for_each(|x| x.reverse());
            e
        };
        let patch = match &plane {
            Some(_) => None,
            None => {
                let patch = patch_under(&patches, &p)?;
                let [du, dv] = patch.degrees();
                if (du + dv) * (p[3].len() - 1) + e[3].len() > MAX_DEGREE {
                    return None;
                }
                Some(patch.clone())
            }
        };
        pieces.push(Piece {
            edge: e,
            uv: p,
            patch,
        });
    }
    Some(Exact { plane, pieces })
}
