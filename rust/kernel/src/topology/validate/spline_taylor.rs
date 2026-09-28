//! Certified deviation of spline uses by second-order Taylor enclosures
//! (S4b of REVIEW_NOTES.md), for every use exact composition does not
//! decide: an arc edge or pcurve, an analytic curved surface under spline
//! geometry, or a spline pcurve across a spline surface's knot lines. On a
//! piece of the use, in its local parameter `τ` in `[0, 1]`,
//! `D(τ) = C(τ) - S(P(τ))` obeys
//! `|D| <= |D(½)| + |D'(½)|/2 + sup|D''|/8`, with
//! `D' = C' - S_u P_u' - S_v P_v'` and
//! `D'' = C'' - (S_uu P_u'^2 + 2 S_uv P_u' P_v' + S_vv P_v'^2 + S_u P_u'' + S_v P_v'')`.
//! Values and first derivatives at `½` are tight; second derivatives over
//! the piece are enclosed from the certified trigonometry (lines, arcs,
//! analytic surfaces) or from the hulls of Bézier derivative nets (splines,
//! through the quotient rule). The pieces and patches are extracted exactly
//! once and lifted into each arithmetic tier, where they are halved and
//! trimmed by interval de Casteljau: every enclosure stays sound. A spline
//! surface under a box that crosses knot lines takes the hull over the
//! patches it meets: the surface is C1 there (R4) and piecewise C2, so the
//! bound still holds. A piece that does not meet the tolerance is halved; a
//! midpoint beyond it certifies a failure.
use super::bernstein::{c, derivative, exact, halves, r, value, Bern};
use super::{tiered, Verdict};
use crate::certified::{Fast, Interval as I, Real};
use crate::curve::ExactBezierCurve3;
use crate::surface::ExactBezierSurface3;
use crate::topology::{Curve2, Curve3, Surface};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// Halvings of one piece before the bound is left undecided, per tier.
const DEPTH_FAST: usize = 14;
const DEPTH_EXACT: usize = 6;

/// The smallest interval holding every enclosure.
fn hull<T: Real>(xs: &[T]) -> T {
    xs[1..].iter().fold(xs[0].clone(), |acc, x| acc.union(x))
}

fn norm<T: Real>(v: &[T]) -> T {
    v.iter()
        .fold(T::exact_f64(0.0), |s, x| s.add(&x.square()))
        .sqrt()
}

/// One piece of a curve in its local parameter `τ` in `[0, 1]`, in `N`
/// coordinates, in one arithmetic tier.
#[derive(Clone)]
enum Piece<T, const N: usize> {
    Line([T; N], [T; N]),
    /// `center + radius (cos θ x + sin θ y)`, `θ` from `a0` to `a1`; `x`
    /// and `y` are the frame's axes (the plane's for a pcurve).
    Arc {
        center: [T; N],
        x: [T; N],
        y: [T; N],
        radius: T,
        a0: T,
        a1: T,
    },
    /// Homogeneous Bernstein coordinates (the last is the weight).
    Bezier(Vec<Bern<T>>),
}

/// The same, exact, as extracted.
#[derive(Clone)]
enum Exact<const N: usize> {
    Line([R; N], [R; N]),
    Arc {
        center: [R; N],
        x: [R; N],
        y: [R; N],
        radius: R,
        a0: R,
        a1: R,
    },
    Bezier(Vec<Vec<R>>),
}

impl<const N: usize> Exact<N> {
    fn lift<T: Real>(&self) -> Piece<T, N> {
        let v = |a: &[R; N]| -> [T; N] { std::array::from_fn(|k| c(&a[k])) };
        match self {
            Exact::Line(a, b) => Piece::Line(v(a), v(b)),
            Exact::Arc {
                center,
                x,
                y,
                radius,
                a0,
                a1,
            } => Piece::Arc {
                center: v(center),
                x: v(x),
                y: v(y),
                radius: c(radius),
                a0: c(a0),
                a1: c(a1),
            },
            Exact::Bezier(h) => {
                Piece::Bezier(h.iter().map(|b| b.iter().map(c).collect()).collect())
            }
        }
    }

    fn reversed(&self) -> Self {
        match self {
            Exact::Line(a, b) => Exact::Line(b.clone(), a.clone()),
            Exact::Arc {
                center,
                x,
                y,
                radius,
                a0,
                a1,
            } => Exact::Arc {
                center: center.clone(),
                x: x.clone(),
                y: y.clone(),
                radius: radius.clone(),
                a0: a1.clone(),
                a1: a0.clone(),
            },
            Exact::Bezier(h) => Exact::Bezier(
                h.iter()
                    .map(|x| x.iter().rev().cloned().collect())
                    .collect(),
            ),
        }
    }
}

/// A value, first and second derivative (or their enclosures).
type Jet<T, const N: usize> = [[T; N]; 3];

impl<T: Real, const N: usize> Piece<T, N> {
    fn halves(&self) -> (Self, Self) {
        let half = T::exact_f64(0.5);
        match self {
            Piece::Line(a, b) => {
                let m: [T; N] = std::array::from_fn(|k| a[k].add(&b[k]).mul(&half));
                (Piece::Line(a.clone(), m.clone()), Piece::Line(m, b.clone()))
            }
            Piece::Arc {
                center,
                x,
                y,
                radius,
                a0,
                a1,
            } => {
                let m = a0.add(a1).mul(&half);
                let arc = |a0: &T, a1: &T| Piece::Arc {
                    center: center.clone(),
                    x: x.clone(),
                    y: y.clone(),
                    radius: radius.clone(),
                    a0: a0.clone(),
                    a1: a1.clone(),
                };
                (arc(a0, &m), arc(&m, a1))
            }
            Piece::Bezier(h) => {
                let parts: Vec<(Bern<T>, Bern<T>)> = h.iter().map(halves).collect();
                (
                    Piece::Bezier(parts.iter().map(|p| p.0.clone()).collect()),
                    Piece::Bezier(parts.iter().map(|p| p.1.clone()).collect()),
                )
            }
        }
    }

    /// The jet at `τ = ½`, or its enclosure over the piece.
    fn jet(&self, over: bool) -> Option<Jet<T, N>> {
        let zero_v: [T; N] = std::array::from_fn(|_| T::exact_f64(0.0));
        let half = T::exact_f64(0.5);
        match self {
            Piece::Line(a, b) => {
                let value = std::array::from_fn(|k| {
                    if over {
                        a[k].union(&b[k])
                    } else {
                        a[k].add(&b[k]).mul(&half)
                    }
                });
                Some([value, std::array::from_fn(|k| b[k].sub(&a[k])), zero_v])
            }
            Piece::Arc {
                center,
                x,
                y,
                radius,
                a0,
                a1,
            } => {
                let theta = if over {
                    a0.union(a1)
                } else {
                    a0.add(a1).mul(&half)
                };
                let (co, si) = T::cos_sin(&theta);
                let sweep = a1.sub(a0);
                let at = |k: usize, p: &T, q: &T| x[k].mul(p).add(&y[k].mul(q));
                let value = std::array::from_fn(|k| center[k].add(&radius.mul(&at(k, &co, &si))));
                let d1 = std::array::from_fn(|k| radius.mul(&sweep).mul(&at(k, &si.neg(), &co)));
                let d2 = std::array::from_fn(|k| {
                    radius.mul(&sweep.square()).mul(&at(k, &co, &si)).neg()
                });
                Some([value, d1, d2])
            }
            Piece::Bezier(h) => {
                let hs: Vec<[Bern<T>; 3]> = h
                    .iter()
                    .map(|b| {
                        let d = derivative(b);
                        let dd = derivative(&d);
                        [b.clone(), d, dd]
                    })
                    .collect();
                let at = |b: &Bern<T>| -> T {
                    if over {
                        hull(b)
                    } else {
                        value(b, &half)
                    }
                };
                let w = [at(&hs[N][0]), at(&hs[N][1]), at(&hs[N][2])];
                if w[0].sign() != Some(Ordering::Greater) {
                    return None;
                }
                let mut value_v = zero_v.clone();
                let mut d1 = zero_v.clone();
                let mut d2 = zero_v;
                for k in 0..N {
                    let (hk, dk, ddk) = (at(&hs[k][0]), at(&hs[k][1]), at(&hs[k][2]));
                    let v = hk.div(&w[0])?;
                    let v1 = dk.sub(&v.mul(&w[1])).div(&w[0])?;
                    let v2 = ddk
                        .sub(&v1.mul(&w[1]).mul(&T::exact_f64(2.0)))
                        .sub(&v.mul(&w[2]))
                        .div(&w[0])?;
                    value_v[k] = v;
                    d1[k] = v1;
                    d2[k] = v2;
                }
                Some([value_v, d1, d2])
            }
        }
    }
}

/// The surface's value and five partials, `S, S_u, S_v, S_uu, S_uv,
/// S_vv`, over the box `(u, v)`.
pub(super) type SurfaceJet<T> = [[T; 3]; 6];

fn analytic_jet<T: Real>(s: &Surface, u: &T, v: &T) -> Option<SurfaceJet<T>> {
    let zero = || T::exact_f64(0.0);
    let z3 = || [zero(), zero(), zero()];
    let (f, kind) = match s {
        Surface::Plane(f) => (f, 0),
        Surface::Cylinder { frame, .. } => (frame, 1),
        Surface::Cone { frame, .. } => (frame, 2),
        Surface::Sphere { frame, .. } => (frame, 3),
        Surface::Torus { frame, .. } => (frame, 4),
        Surface::BSpline(_) => return None,
    };
    let fr = super::frame::<T>(f);
    let (co, si) = T::cos_sin(u);
    let e: [T; 3] = std::array::from_fn(|k| fr.x[k].mul(&co).add(&fr.y[k].mul(&si)));
    let e1: [T; 3] = std::array::from_fn(|k| fr.y[k].mul(&co).sub(&fr.x[k].mul(&si)));
    let comb = |terms: &[(&[T; 3], T)]| -> [T; 3] {
        std::array::from_fn(|k| {
            terms
                .iter()
                .fold(zero(), |acc, (vec, s)| acc.add(&vec[k].mul(s)))
        })
    };
    let one = T::exact_f64(1.0);
    Some(match (kind, s) {
        (0, _) => [
            comb(&[(&fr.o, one.clone()), (&fr.x, u.clone()), (&fr.y, v.clone())]),
            fr.x.clone(),
            fr.y.clone(),
            z3(),
            z3(),
            z3(),
        ],
        (1, Surface::Cylinder { radius, .. }) => {
            let r = T::exact_f64(*radius);
            [
                comb(&[(&fr.o, one), (&e, r.clone()), (&fr.n, v.clone())]),
                comb(&[(&e1, r.clone())]),
                fr.n.clone(),
                comb(&[(&e, r.neg())]),
                z3(),
                z3(),
            ]
        }
        (
            2,
            Surface::Cone {
                radius, half_angle, ..
            },
        ) => {
            let (ca, sa) = T::cos_sin(&T::exact_f64(*half_angle));
            let rho = T::exact_f64(*radius).add(&sa.mul(v));
            [
                comb(&[(&fr.o, one), (&e, rho.clone()), (&fr.n, ca.mul(v))]),
                comb(&[(&e1, rho.clone())]),
                comb(&[(&e, sa.clone()), (&fr.n, ca)]),
                comb(&[(&e, rho.neg())]),
                comb(&[(&e1, sa)]),
                z3(),
            ]
        }
        (3, Surface::Sphere { radius, .. }) => {
            let r = T::exact_f64(*radius);
            let (cv, sv) = T::cos_sin(v);
            [
                comb(&[(&fr.o, one), (&e, r.mul(&cv)), (&fr.n, r.mul(&sv))]),
                comb(&[(&e1, r.mul(&cv))]),
                comb(&[(&e, r.mul(&sv).neg()), (&fr.n, r.mul(&cv))]),
                comb(&[(&e, r.mul(&cv).neg())]),
                comb(&[(&e1, r.mul(&sv).neg())]),
                comb(&[(&e, r.mul(&cv).neg()), (&fr.n, r.mul(&sv).neg())]),
            ]
        }
        (4, Surface::Torus { major, minor, .. }) => {
            let (big, small) = (T::exact_f64(*major), T::exact_f64(*minor));
            let (cv, sv) = T::cos_sin(v);
            let rho = big.add(&small.mul(&cv));
            [
                comb(&[(&fr.o, one), (&e, rho.clone()), (&fr.n, small.mul(&sv))]),
                comb(&[(&e1, rho.clone())]),
                comb(&[(&e, small.mul(&sv).neg()), (&fr.n, small.mul(&cv))]),
                comb(&[(&e, rho.neg())]),
                comb(&[(&e1, small.mul(&sv).neg())]),
                comb(&[(&e, small.mul(&cv).neg()), (&fr.n, small.mul(&sv).neg())]),
            ]
        }
        _ => unreachable!("kind matches the surface"),
    })
}

/// A Bézier patch lifted into a tier: its domain and its homogeneous tensor
/// net, `[coordinate][i][j]`.
#[derive(Clone)]
pub(super) struct Patch<T> {
    lifted: [[T; 2]; 2],
    /// Per homogeneous coordinate, the nets of `h, h_u, h_v, h_uu, h_uv,
    /// h_vv` in the patch's local parameters, computed once.
    partials: Vec<[Vec<Vec<T>>; 6]>,
    /// Whether each side (`[axis][low, high]`) is the surface domain's:
    /// there the patch's polynomial extends past it, as OCCT evaluates a
    /// pcurve printed a hair outside the domain; elsewhere a neighbour
    /// patch covers the box.
    open: [[bool; 2]; 2],
}

pub(super) fn lift_patches<T: Real>(patches: &[ExactBezierSurface3]) -> Vec<Patch<T>> {
    lift_patches_about(patches, None)
}

/// `lift_patches` with the surface translated by `-origin` first, exactly
/// (`X - origin w` on the homogeneous poles): a rational patch's jets are
/// quotients of enclosures, whose widths grow with `|S|`, so integrands
/// relative to a point near the body stay at its own scale.
pub(super) fn lift_patches_about<T: Real>(
    patches: &[ExactBezierSurface3],
    origin: Option<&[R; 3]>,
) -> Vec<Patch<T>> {
    let bound = |axis: usize, side: usize| {
        let values = patches.iter().map(|q| q.domain()[axis][side].clone());
        if side == 0 {
            values.min()
        } else {
            values.max()
        }
    };
    let outer = [[bound(0, 0), bound(0, 1)], [bound(1, 0), bound(1, 1)]];
    patches
        .iter()
        .map(|q| {
            let mut p = lift_patch(q, origin);
            p.open = std::array::from_fn(|axis| {
                std::array::from_fn(|side| {
                    Some(&q.domain()[axis][side]) == outer[axis][side].as_ref()
                })
            });
            p
        })
        .collect()
}

fn lift_patch<T: Real>(q: &ExactBezierSurface3, origin: Option<&[R; 3]>) -> Patch<T> {
    let [du, dv] = q.degrees();
    let mut poles = q.homogeneous_poles().to_vec();
    if let Some(o) = origin {
        for p in &mut poles {
            for k in 0..3 {
                p[k] = &p[k] - &o[k] * &p[3];
            }
        }
    }
    let domain = q.domain();
    let d_u = |a: &Vec<Vec<T>>| -> Vec<Vec<T>> {
        if a.len() == 1 {
            return vec![vec![T::exact_f64(0.0); a[0].len()]];
        }
        let m = T::exact_f64((a.len() - 1) as f64);
        a.windows(2)
            .map(|w| {
                w[1].iter()
                    .zip(&w[0])
                    .map(|(p, q)| p.sub(q).mul(&m))
                    .collect()
            })
            .collect()
    };
    let d_v = |a: &Vec<Vec<T>>| -> Vec<Vec<T>> { a.iter().map(derivative).collect() };
    let partials = (0..4)
        .map(|k| {
            let h: Vec<Vec<T>> = (0..=du)
                .map(|i| (0..=dv).map(|j| c(&poles[i * (dv + 1) + j][k])).collect())
                .collect();
            let (hu, hv) = (d_u(&h), d_v(&h));
            let (huu, huv, hvv) = (d_u(&hu), d_v(&hu), d_v(&hv));
            [h, hu, hv, huu, huv, hvv]
        })
        .collect();
    Patch {
        partials,
        open: [[false; 2]; 2],
        lifted: std::array::from_fn(|a| std::array::from_fn(|b| c(&domain[a][b]))),
    }
}

/// The Bernstein coefficients of a polynomial restricted to `[a, b]`, by
/// two de Casteljau splits (valid past `[0, 1]`, the polynomial's
/// extension): at `b` then at `a / b` of the left part, or at `a` then at
/// `(b - a) / (1 - a)` of the right part, whichever divides by the larger
/// magnitude. In interval arithmetic the result encloses the exact one.
fn restricted<T: Real>(coeffs: &[T], a: &T, b: &T) -> Vec<T> {
    let split = |cs: &[T], t: &T| -> (Vec<T>, Vec<T>) {
        let s = T::exact_f64(1.0).sub(t);
        let mut row = cs.to_vec();
        let (mut left, mut right) = (vec![row[0].clone()], vec![row[row.len() - 1].clone()]);
        while row.len() > 1 {
            row = row
                .windows(2)
                .map(|w| s.mul(&w[0]).add(&t.mul(&w[1])))
                .collect();
            left.push(row[0].clone());
            right.push(row[row.len() - 1].clone());
        }
        right.reverse();
        (left, right)
    };
    let one = T::exact_f64(1.0);
    let magnitude = |x: &T| {
        let (lo, hi) = x.bounds_f64();
        lo.abs().min(hi.abs())
    };
    let rest = one.sub(a);
    if magnitude(b) >= magnitude(&rest) {
        let (left, _) = split(coeffs, b);
        match a.div(b) {
            Some(t) => split(&left, &t).1,
            None => left,
        }
    } else {
        let (_, right) = split(coeffs, a);
        match b.sub(a).div(&rest) {
            Some(t) => split(&right, &t).0,
            None => right,
        }
    }
}

/// A spline surface's jet over a box: the hull over the patches it meets,
/// each restricted to the box and enclosed by its derivative nets.
pub(super) fn spline_jet<T: Real>(patches: &[Patch<T>], u: &T, v: &T) -> Option<SurfaceJet<T>> {
    spline_jet_to(patches, u, v, true)
}

/// `S`, `S_u` and `S_v` of a spline surface over a box, as `spline_jet`
/// without the second derivatives (the flux and mass strips need none).
pub(super) fn spline_jet1<T: Real>(patches: &[Patch<T>], u: &T, v: &T) -> Option<[[T; 3]; 3]> {
    let [s, su, sv, ..] = spline_jet_to(patches, u, v, false)?;
    Some([s, su, sv])
}

/// The jet over a box; its second derivatives are zero unless `second`.
fn spline_jet_to<T: Real>(
    patches: &[Patch<T>],
    u: &T,
    v: &T,
    second: bool,
) -> Option<SurfaceJet<T>> {
    let (ua, ub) = u.bounds_f64();
    let (va, vb) = v.bounds_f64();
    let mut out: Option<SurfaceJet<T>> = None;
    for q in patches {
        let [[pu0, pu1], [pv0, pv1]] = &q.lifted;
        // Certainly disjoint patches are skipped (binary64 bounds are
        // outward, so a patch that may meet the box is kept).
        if ub < pu0.bounds_f64().0
            || ua > pu1.bounds_f64().1
            || vb < pv0.bounds_f64().0
            || va > pv1.bounds_f64().1
        {
            continue;
        }
        let jet = patch_jet(q, u, v, second)?;
        out = Some(match out {
            None => jet,
            Some(acc) => {
                std::array::from_fn(|i| std::array::from_fn(|k| acc[i][k].union(&jet[i][k])))
            }
        });
    }
    out
}

/// The enclosure of a patch's jet over the box, in global parameters.
fn patch_jet<T: Real>(q: &Patch<T>, u: &T, v: &T, second: bool) -> Option<SurfaceJet<T>> {
    let [[u0, u1], [v0, v1]] = &q.lifted;
    let (wu, wv) = (u1.sub(u0), v1.sub(v0));
    // The box in local coordinates, clipped to the patch except past the
    // surface domain's own sides (the polynomial extension).
    let local = |x: &T, lo: &T, w: &T, open: [bool; 2]| -> Option<(T, T)> {
        let (a, b) = x.bounds_f64();
        let (a, b) = (
            T::exact_f64(a).sub(lo).div(w)?,
            T::exact_f64(b).sub(lo).div(w)?,
        );
        let low = if open[0] { f64::NEG_INFINITY } else { 0.0 };
        let high = if open[1] { f64::INFINITY } else { 1.0 };
        let a = T::exact_f64(a.bounds_f64().0.clamp(low, high));
        let b = T::exact_f64(b.bounds_f64().1.clamp(low, high));
        Some((a, b))
    };
    let (a_u, b_u) = local(u, u0, &wu, q.open[0])?;
    let (a_v, b_v) = local(v, v0, &wv, q.open[1])?;
    // Restrict every row in v, then every column in u.
    let restrict = |net: &Vec<Vec<T>>| -> Vec<Vec<T>> {
        let rows: Vec<Vec<T>> = net.iter().map(|row| restricted(row, &a_v, &b_v)).collect();
        let n = rows[0].len();
        let cols: Vec<Vec<T>> = (0..n)
            .map(|j| {
                let col: Vec<T> = rows.iter().map(|r| r[j].clone()).collect();
                restricted(&col, &a_u, &b_u)
            })
            .collect();
        (0..cols[0].len())
            .map(|i| cols.iter().map(|c| c[i].clone()).collect())
            .collect()
    };
    // The hull of a derivative net of the whole patch restricted to the box.
    let flat = |a: &Vec<Vec<T>>| -> T {
        let xs: Vec<T> = restrict(a).into_iter().flatten().collect();
        hull(&xs)
    };
    // [h, h_u, h_v, h_uu, h_uv, h_vv] of a coordinate, in global units (the
    // patch's local derivatives divided by its widths).
    let partials = |k: usize| -> Option<[T; 6]> {
        let [h, hu, hv, huu, huv, hvv] = &q.partials[k];
        if !second {
            let zero = T::exact_f64(0.0);
            return Some([
                flat(h),
                flat(hu).div(&wu)?,
                flat(hv).div(&wv)?,
                zero.clone(),
                zero.clone(),
                zero,
            ]);
        }
        Some([
            flat(h),
            flat(hu).div(&wu)?,
            flat(hv).div(&wv)?,
            flat(huu).div(&wu.square())?,
            flat(huv).div(&wu.mul(&wv))?,
            flat(hvv).div(&wv.square())?,
        ])
    };
    let w = partials(3)?;
    if w[0].sign() != Some(Ordering::Greater) {
        return None;
    }
    let two = T::exact_f64(2.0);
    let column = |k: usize| -> Option<[T; 6]> {
        let h = partials(k)?;
        let s = h[0].div(&w[0])?;
        let s_u = h[1].sub(&s.mul(&w[1])).div(&w[0])?;
        let s_v = h[2].sub(&s.mul(&w[2])).div(&w[0])?;
        let s_uu = h[3]
            .sub(&two.mul(&s_u).mul(&w[1]))
            .sub(&s.mul(&w[3]))
            .div(&w[0])?;
        let s_uv = h[4]
            .sub(&s_u.mul(&w[2]))
            .sub(&s_v.mul(&w[1]))
            .sub(&s.mul(&w[4]))
            .div(&w[0])?;
        let s_vv = h[5]
            .sub(&two.mul(&s_v).mul(&w[2]))
            .sub(&s.mul(&w[5]))
            .div(&w[0])?;
        Some([s, s_u, s_v, s_uu, s_uv, s_vv])
    };
    let columns = [column(0)?, column(1)?, column(2)?];
    Some(std::array::from_fn(|i| {
        std::array::from_fn(|k| columns[k][i].clone())
    }))
}

/// A use as exact pieces of its edge and pcurve, over the surface.
pub(super) struct Taylor {
    surface: Surface,
    patches: Vec<ExactBezierSurface3>,
    pieces: Vec<(Exact<3>, Exact<2>)>,
}

/// The bound's three terms on a piece, or `None` when a jet is undefined.
fn terms<T: Real>(
    surface: &Surface,
    patches: &[Patch<T>],
    edge: &Piece<T, 3>,
    uv: &Piece<T, 2>,
) -> Option<(T, T, T)> {
    let [c_mid, c1_mid, _] = edge.jet(false)?;
    let [_, _, c2] = edge.jet(true)?;
    let [p_mid, p1_mid, _] = uv.jet(false)?;
    let [p_over, p1, p2] = uv.jet(true)?;
    let jet = |u: &T, v: &T| match surface {
        Surface::BSpline(_) => spline_jet(patches, u, v),
        s => analytic_jet(s, u, v),
    };
    let s_mid = jet(&p_mid[0], &p_mid[1])?;
    let s_over = jet(&p_over[0], &p_over[1])?;
    let d0: [T; 3] = std::array::from_fn(|k| c_mid[k].sub(&s_mid[0][k]));
    let d1: [T; 3] = std::array::from_fn(|k| {
        c1_mid[k]
            .sub(&s_mid[1][k].mul(&p1_mid[0]))
            .sub(&s_mid[2][k].mul(&p1_mid[1]))
    });
    let two = T::exact_f64(2.0);
    let d2: [T; 3] = std::array::from_fn(|k| {
        let [_, su, sv, suu, suv, svv] = &s_over;
        let chain = suu[k]
            .mul(&p1[0].square())
            .add(&two.mul(&suv[k]).mul(&p1[0]).mul(&p1[1]))
            .add(&svv[k].mul(&p1[1].square()))
            .add(&su[k].mul(&p2[0]))
            .add(&sv[k].mul(&p2[1]));
        c2[k].sub(&chain)
    });
    Some((norm(&d0), norm(&d1), norm(&d2)))
}

/// A lifted edge piece and its pcurve piece over the same fraction range.
type PiecePair<T> = (Piece<T, 3>, Piece<T, 2>);

impl Taylor {
    fn lifted<T: Real>(&self) -> (Vec<Patch<T>>, Vec<PiecePair<T>>) {
        (
            lift_patches(&self.patches),
            self.pieces
                .iter()
                .map(|(e, p)| (e.lift(), p.lift()))
                .collect(),
        )
    }

    fn decide_in<T: Real>(&self, tol: f64, depth_limit: usize) -> Verdict {
        let tol = T::exact_f64(tol);
        let (patches, pieces) = self.lifted::<T>();
        let mut stack: Vec<(Piece<T, 3>, Piece<T, 2>, usize)> =
            pieces.into_iter().map(|(e, p)| (e, p, 0)).collect();
        let mut undecided = false;
        while let Some((edge, uv, depth)) = stack.pop() {
            // An undefined jet (off the surface's domain, a weight not
            // certainly positive) leaves the piece undecided.
            let Some((d0, d1, d2)) = terms(&self.surface, &patches, &edge, &uv) else {
                undecided = true;
                continue;
            };
            if d0.cmp(&tol) == Some(Ordering::Greater) {
                return Verdict::Beyond;
            }
            let bound = d0
                .add(&d1.mul(&T::exact_f64(0.5)))
                .add(&d2.mul(&T::exact_f64(0.125)));
            if matches!(bound.cmp(&tol), Some(Ordering::Less | Ordering::Equal)) {
                continue;
            }
            if depth == depth_limit {
                undecided = true;
                continue;
            }
            let (e0, e1) = edge.halves();
            let (p0, p1) = uv.halves();
            stack.push((e0, p0, depth + 1));
            stack.push((e1, p1, depth + 1));
        }
        if undecided {
            Verdict::Unknown
        } else {
            Verdict::Within
        }
    }

    pub(super) fn decide(&self, tol: f64) -> Verdict {
        tiered(
            Verdict::Unknown,
            || self.decide_in::<Fast>(tol, DEPTH_FAST),
            || self.decide_in::<I>(tol, DEPTH_EXACT),
        )
    }

    /// A certified upper bound of `|D|`, refined up to eight halvings while
    /// the remainder dominates the bound (binary64 intervals).
    pub(super) fn upper_bound(&self) -> Option<f64> {
        let (patches, pieces) = self.lifted::<Fast>();
        let mut worst = 0.0_f64;
        let mut stack: Vec<(Piece<Fast, 3>, Piece<Fast, 2>, usize)> =
            pieces.into_iter().map(|(e, p)| (e, p, 0)).collect();
        while let Some((edge, uv, depth)) = stack.pop() {
            let (d0, d1, d2) = terms(&self.surface, &patches, &edge, &uv)?;
            let point = d0.add(&d1.mul(&Fast::exact_f64(0.5))).bounds_f64().1;
            let rest = d2.mul(&Fast::exact_f64(0.125)).bounds_f64().1;
            if rest > point && depth < 8 {
                let (e0, e1) = edge.halves();
                let (p0, p1) = uv.halves();
                stack.push((e0, p0, depth + 1));
                stack.push((e1, p1, depth + 1));
                continue;
            }
            let bound = Fast::exact_f64(point)
                .add(&Fast::exact_f64(rest))
                .bounds_f64()
                .1;
            if !bound.is_finite() {
                return None;
            }
            worst = worst.max(bound);
        }
        Some(worst)
    }
}

/// A line, arc or spline in `N` coordinates, over its whole fraction.
enum CurveKind<const N: usize> {
    Line([R; N], [R; N]),
    Arc {
        center: [R; N],
        x: [R; N],
        y: [R; N],
        radius: R,
        start: R,
        sweep: R,
    },
    Spline(Vec<(R, R, ExactBezierCurve3)>),
}

impl<const N: usize> CurveKind<N> {
    fn piece(&self, f0: &R, f1: &R) -> Option<Exact<N>> {
        Some(match self {
            CurveKind::Line(a, b) => {
                let at = |f: &R| -> [R; N] { std::array::from_fn(|k| &a[k] + (&b[k] - &a[k]) * f) };
                Exact::Line(at(f0), at(f1))
            }
            CurveKind::Arc {
                center,
                x,
                y,
                radius,
                start,
                sweep,
            } => Exact::Arc {
                center: center.clone(),
                x: x.clone(),
                y: y.clone(),
                radius: radius.clone(),
                a0: start + sweep * f0,
                a1: start + sweep * f1,
            },
            CurveKind::Spline(arcs) => {
                let h = super::bernstein::piece_of(arcs, f0, f1)?;
                // x, y (, z) and the weight.
                let keep: Vec<usize> = if N == 2 {
                    vec![0, 1, 3]
                } else {
                    vec![0, 1, 2, 3]
                };
                Exact::Bezier(keep.into_iter().map(|k| h[k].clone()).collect())
            }
        })
    }

    fn knots(&self) -> Vec<R> {
        match self {
            CurveKind::Spline(arcs) => arcs
                .iter()
                .flat_map(|(a, b, _)| [a.clone(), b.clone()])
                .collect(),
            _ => Vec::new(),
        }
    }
}

fn pieces<const N: usize>(curve: &CurveKind<N>, cuts: &[R]) -> Option<Vec<Exact<N>>> {
    cuts.windows(2).map(|w| curve.piece(&w[0], &w[1])).collect()
}

fn edge_kind(curve: &Curve3) -> Option<CurveKind<3>> {
    Some(match curve {
        Curve3::LineSegment { start, end } => {
            CurveKind::Line(start.to_array().map(r), end.to_array().map(r))
        }
        Curve3::Circle { frame, radius } => CurveKind::Arc {
            center: frame.origin().to_array().map(r),
            x: frame.x().to_array().map(r),
            y: frame.y().to_array().map(r),
            radius: r(*radius),
            start: exact(0),
            sweep: r(std::f64::consts::TAU),
        },
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => CurveKind::Arc {
            center: frame.origin().to_array().map(r),
            x: frame.x().to_array().map(r),
            y: frame.y().to_array().map(r),
            radius: r(*radius),
            start: r(*start_angle),
            sweep: r(*sweep_angle),
        },
        Curve3::BSpline(s) => CurveKind::Spline(super::bernstein::edge_arcs(s)?),
        Curve3::HyperbolaArc { .. } | Curve3::ParabolaArc { .. } | Curve3::Section(_) => {
            return None
        }
        // The unit circle's arc on the semi-axes as axes.
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => CurveKind::Arc {
            center: frame.origin().to_array().map(r),
            x: frame.x().to_array().map(|x| r(x) * r(*major)),
            y: frame.y().to_array().map(|y| r(y) * r(*minor)),
            radius: exact(1),
            start: r(*start_angle),
            sweep: r(*sweep_angle),
        },
    })
}

fn pcurve_kind(p: &Curve2) -> Option<CurveKind<2>> {
    Some(match p {
        Curve2::LineSegment { start, end } => {
            CurveKind::Line([r(start.x), r(start.y)], [r(end.x), r(end.y)])
        }
        Curve2::CircularArc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => CurveKind::Arc {
            center: [r(center.x), r(center.y)],
            x: [exact(1), exact(0)],
            y: [exact(0), exact(1)],
            radius: r(*radius),
            start: r(*start_angle),
            sweep: r(*sweep_angle),
        },
        Curve2::BSpline(s) => CurveKind::Spline(super::bernstein::span_arcs(s)?),
        Curve2::EllipseArc {
            center,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => CurveKind::Arc {
            center: [r(center.x), r(center.y)],
            x: [r(*major), exact(0)],
            y: [exact(0), r(*minor)],
            radius: exact(1),
            start: r(*start_angle),
            sweep: r(*sweep_angle),
        },
        Curve2::Sinusoid { .. } | Curve2::Projection(_) => return None,
    })
}

/// The Taylor form of a use, or `None` on a periodic spline surface.
pub(super) fn taylor_use(
    curve: &Curve3,
    surface: &Surface,
    pcurve: &Curve2,
    forward: bool,
) -> Option<Taylor> {
    let patches = match surface {
        Surface::BSpline(s) => {
            if s.u_knots().is_periodic() || s.v_knots().is_periodic() {
                return None;
            }
            s.bezier_patches().ok()?
        }
        _ => Vec::new(),
    };
    let edge = edge_kind(curve)?;
    let uv = pcurve_kind(pcurve)?;
    let one = exact(1);
    let mut cuts: Vec<R> = vec![exact(0), one.clone()];
    cuts.extend(uv.knots());
    cuts.extend(
        edge.knots()
            .into_iter()
            .map(|s| if forward { s } else { &one - s }),
    );
    cuts.sort();
    cuts.dedup();
    let uv_pieces = pieces(&uv, &cuts)?;
    let edge_pieces: Vec<Exact<3>> = if forward {
        pieces(&edge, &cuts)?
    } else {
        let flipped: Vec<R> = cuts.iter().rev().map(|t| &one - t).collect();
        pieces(&edge, &flipped)?
            .into_iter()
            .rev()
            .map(|p| p.reversed())
            .collect()
    };
    Some(Taylor {
        surface: surface.clone(),
        patches,
        pieces: edge_pieces.into_iter().zip(uv_pieces).collect(),
    })
}
