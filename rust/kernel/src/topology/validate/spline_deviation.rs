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
    /// A bound of the surface's departure from the patch's polynomial
    /// where the pcurve piece's hull crosses into a neighbouring patch by a
    /// sliver (zero when it lies in the patch).
    extra: f64,
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
        let samples: Vec<T> = (0..=4).map(|k| c(&ratio(k, 4))).collect();
        let mut stack = Vec::new();
        let mut undecided = false;
        for piece in &self.pieces {
            // A sliver's departure narrows what the polynomial may deviate
            // by within the tolerance and widens it beyond.
            let (inner, outer) = (
                T::exact_f64(tol).sub(&T::exact_f64(piece.extra)),
                T::exact_f64(tol).add(&T::exact_f64(piece.extra)),
            );
            if inner.sign() != Some(Ordering::Greater) {
                undecided = true;
                continue;
            }
            let (inner2, outer2) = (inner.square(), outer.square());
            let q = self.quotient::<T>(piece);
            if samples.iter().any(|t| q.beyond_at(t, &outer2)) {
                return Verdict::Beyond;
            }
            stack.push((q, 0, inner2, outer2));
        }
        let half = T::exact_f64(0.5);
        while let Some((q, depth, inner2, outer2)) = stack.pop() {
            if q.within(&inner2) {
                continue;
            }
            if q.beyond_at(&half, &outer2) {
                return Verdict::Beyond;
            }
            if depth == MAX_DEPTH {
                undecided = true;
                continue;
            }
            let (a, b) = q.halves();
            stack.push((a, depth + 1, inner2.clone(), outer2.clone()));
            stack.push((b, depth + 1, inner2, outer2));
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
                    .add(&Fast::exact_f64(piece.extra))
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

/// A patch's nonrational control net re-expressed over a box of its own
/// parameters (inside or outside its domain: a polynomial extends), by
/// blossoming each row and column; `None` for a rational patch.
fn over_box(patch: &ExactBezierSurface3, b: &[[R; 2]; 2]) -> Option<Vec<[R; 3]>> {
    let [du, dv] = patch.degrees();
    let [[u0, u1], [v0, v1]] = patch.domain().clone();
    let poles = patch.homogeneous_poles();
    let w0 = &poles[0][3];
    if poles.iter().any(|q| &q[3] != w0) {
        return None;
    }
    let net: Vec<[R; 3]> = poles
        .iter()
        .map(|q| std::array::from_fn(|k| &q[k] / w0))
        .collect();
    // The control points of a degree-n polynomial given over [0, 1] by `cs`
    // over [a, b]: its blossom at (a^(n-i), b^i).
    let reexpress = |cs: &[[R; 3]], a: &R, b: &R| -> Vec<[R; 3]> {
        let n = cs.len() - 1;
        (0..=n)
            .map(|i| {
                let mut row = cs.to_vec();
                for step in 0..n {
                    let t = if step < n - i { a } else { b };
                    row = row
                        .windows(2)
                        .map(|w| std::array::from_fn(|k| &w[0][k] + t * (&w[1][k] - &w[0][k])))
                        .collect();
                }
                row[0].clone()
            })
            .collect()
    };
    let local = |x: &R, lo: &R, hi: &R| (x - lo) / (hi - lo);
    let (ua, ub) = (local(&b[0][0], &u0, &u1), local(&b[0][1], &u0, &u1));
    let (va, vb) = (local(&b[1][0], &v0, &v1), local(&b[1][1], &v0, &v1));
    // Rows along v first, then columns along u.
    let mut rows: Vec<Vec<[R; 3]>> = (0..=du)
        .map(|i| reexpress(&net[i * (dv + 1)..(i + 1) * (dv + 1)], &va, &vb))
        .collect();
    for j in 0..=dv {
        let column: Vec<[R; 3]> = rows.iter().map(|r| r[j].clone()).collect();
        for (i, x) in reexpress(&column, &ua, &ub).into_iter().enumerate() {
            rows[i][j] = x;
        }
    }
    Some(rows.into_iter().flatten().collect())
}

/// The patch holding a pcurve piece's control points but for slivers
/// across its boundaries at most `2^-20` of its size wide, with a certified
/// bound of the surface's departure from that patch's polynomial on them:
/// on each neighbouring patch's part of the piece's box, both polynomials
/// re-expressed exactly over it, the largest control point of their
/// difference (the difference lies in its hull). A piece whose control
/// polygon leaves the surface's domain while its curve certainly keeps to
/// it (its homogeneous coordinate against the bound exactly nonnegative: a
/// crease nearly touching a cap, S9f.1, as Green's exact path takes it,
/// S8b.3's (d)) is boxed by that bound; a piece reaching beyond the
/// surface's domain, or across a rational patch, is `None`.
fn patch_near<'a>(
    patches: &'a [ExactBezierSurface3],
    uv: &[Vec<R>; 4],
) -> Option<(&'a ExactBezierSurface3, f64)> {
    if let Some(patch) = patch_under(patches, uv) {
        return Some((patch, 0.0));
    }
    let points: Vec<(R, R)> = (0..uv[3].len())
        .map(|i| (&uv[0][i] / &uv[3][i], &uv[1][i] / &uv[3][i]))
        .collect();
    let lo_hi = |f: &dyn Fn(&(R, R)) -> &R| {
        let lo = points.iter().map(f).min().expect("a point").clone();
        let hi = points.iter().map(f).max().expect("a point").clone();
        [lo, hi]
    };
    let mut bbox = [lo_hi(&|p| &p.0), lo_hi(&|p| &p.1)];
    // Inside the surface's domain: the hull, or the curve where its hull
    // leaves the domain.
    for (axis, [below, above]) in bbox.iter_mut().enumerate() {
        let lo = patches.iter().map(|q| &q.domain()[axis][0]).min()?.clone();
        let hi = patches.iter().map(|q| &q.domain()[axis][1]).max()?.clone();
        let (x, w) = (&uv[axis], &uv[3]);
        if *below < lo {
            let off: Vec<R> = x.iter().zip(w).map(|(x, w)| x - &lo * w).collect();
            if !super::bernstein::nonnegative(&off) {
                return None;
            }
            *below = lo;
        }
        if *above > hi {
            let off: Vec<R> = x.iter().zip(w).map(|(x, w)| &hi * w - x).collect();
            if !super::bernstein::nonnegative(&off) {
                return None;
            }
            *above = hi;
        }
    }
    let two = ratio(2, 1);
    let mid = [
        (&bbox[0][0] + &bbox[0][1]) / &two,
        (&bbox[1][0] + &bbox[1][1]) / &two,
    ];
    let inside = |q: &ExactBezierSurface3, p: &[R; 2]| {
        (0..2).all(|a| q.domain()[a][0] <= p[a] && p[a] <= q.domain()[a][1])
    };
    let patch = patches.iter().find(|q| inside(q, &mid))?;
    let sliver = ratio(1, 1 << 20);
    for (axis, [below, above]) in bbox.iter().enumerate() {
        let [lo, hi] = &patch.domain()[axis];
        let limit = (hi - lo) * &sliver;
        if lo - below > limit || above - hi > limit {
            return None;
        }
    }
    let own_net = |b: &[[R; 2]; 2]| over_box(patch, b);
    let mut extra = 0.0_f64;
    for other in patches {
        if std::ptr::eq(other, patch) {
            continue;
        }
        // The neighbour's part of the box, with room inside it.
        let part: [[R; 2]; 2] = std::array::from_fn(|a| {
            [
                (&bbox[a][0]).max(&other.domain()[a][0]).clone(),
                (&bbox[a][1]).min(&other.domain()[a][1]).clone(),
            ]
        });
        if part.iter().any(|[lo, hi]| lo >= hi) {
            continue;
        }
        let theirs = over_box(other, &part)?;
        let mine = own_net(&part)?;
        for (a, b) in theirs.iter().zip(&mine) {
            let d2: R = (0..3).map(|k| (&a[k] - &b[k]) * (&a[k] - &b[k])).sum();
            let bound = I::exact(d2).sqrt().bounds_f64().1;
            extra = extra.max(bound);
        }
    }
    if !extra.is_finite() {
        return None;
    }
    Some((patch, extra))
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
        let (patch, extra) = match &plane {
            Some(_) => (None, 0.0),
            None => {
                let (patch, extra) = patch_near(&patches, &p)?;
                let [du, dv] = patch.degrees();
                if (du + dv) * (p[3].len() - 1) + e[3].len() > MAX_DEGREE {
                    return None;
                }
                (Some(patch.clone()), extra)
            }
        };
        pieces.push(Piece {
            edge: e,
            uv: p,
            patch,
            extra,
        });
    }
    Some(Exact { plane, pieces })
}
