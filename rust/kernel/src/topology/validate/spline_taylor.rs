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
//! analytic surfaces) or from the hulls of exact Bézier derivative nets
//! (splines, through the quotient rule). A spline surface under a box that
//! crosses knot lines takes the hull over the patches it meets: the surface
//! is C1 there (R4) and piecewise C2, so the bound still holds. A piece that
//! does not meet the tolerance is halved; a midpoint beyond it certifies a
//! failure.
use super::bernstein::{c, derivative, ends, exact, r, ratio, Bern};
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

fn zero() -> R {
    exact(0)
}

fn half(x: &R) -> R {
    x * ratio(1, 2)
}

/// The smallest interval holding every enclosure.
fn hull<T: Real>(xs: &[T]) -> T {
    let bounds: Vec<(R, R)> = xs.iter().map(ends).collect();
    let lo = bounds.iter().map(|b| b.0.clone()).min().unwrap();
    let hi = bounds.iter().map(|b| b.1.clone()).max().unwrap();
    T::from_r(&half(&(&lo + &hi))).widen(&half(&(hi - lo)))
}

/// An enclosure of `[a, b]` (either order).
fn span<T: Real>(a: &R, b: &R) -> T {
    let d = a - b;
    let d = if d < zero() { -d } else { d };
    T::from_r(&half(&(a + b))).widen(&half(&d))
}

fn norm<T: Real>(v: &[T]) -> T {
    v.iter()
        .fold(T::exact_f64(0.0), |s, x| s.add(&x.square()))
        .sqrt()
}

fn exact_halves(a: &[R]) -> (Vec<R>, Vec<R>) {
    let mut row = a.to_vec();
    let (mut left, mut right) = (vec![row[0].clone()], vec![row[row.len() - 1].clone()]);
    while row.len() > 1 {
        row = row.windows(2).map(|w| half(&(&w[0] + &w[1]))).collect();
        left.push(row[0].clone());
        right.push(row[row.len() - 1].clone());
    }
    right.reverse();
    (left, right)
}

/// One piece of a curve in its local parameter `τ` in `[0, 1]`, in `N`
/// coordinates.
#[derive(Clone)]
enum Piece<const N: usize> {
    Line([R; N], [R; N]),
    /// `center + radius (cos θ x + sin θ y)`, `θ` from `a0` to `a1`; `x`
    /// and `y` are the frame's axes (the plane's for a pcurve).
    Arc {
        center: [R; N],
        x: [R; N],
        y: [R; N],
        radius: R,
        a0: R,
        a1: R,
    },
    /// Homogeneous Bernstein coordinates (the last is the weight).
    Bezier(Vec<Vec<R>>),
}

/// A value, first and second derivative (or their enclosures).
type Jet<T, const N: usize> = [[T; N]; 3];

impl<const N: usize> Piece<N> {
    fn halves(&self) -> (Self, Self) {
        match self {
            Piece::Line(a, b) => {
                let m: [R; N] = std::array::from_fn(|k| half(&(&a[k] + &b[k])));
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
                let m = half(&(a0 + a1));
                let arc = |a0: &R, a1: &R| Piece::Arc {
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
                let parts: Vec<(Vec<R>, Vec<R>)> = h.iter().map(|x| exact_halves(x)).collect();
                (
                    Piece::Bezier(parts.iter().map(|p| p.0.clone()).collect()),
                    Piece::Bezier(parts.iter().map(|p| p.1.clone()).collect()),
                )
            }
        }
    }

    fn reversed(&self) -> Self {
        match self {
            Piece::Line(a, b) => Piece::Line(b.clone(), a.clone()),
            Piece::Arc {
                center,
                x,
                y,
                radius,
                a0,
                a1,
            } => Piece::Arc {
                center: center.clone(),
                x: x.clone(),
                y: y.clone(),
                radius: radius.clone(),
                a0: a1.clone(),
                a1: a0.clone(),
            },
            Piece::Bezier(h) => Piece::Bezier(
                h.iter()
                    .map(|x| x.iter().rev().cloned().collect())
                    .collect(),
            ),
        }
    }

    /// The jet at `τ = ½` (`at == None`) or its enclosure over the piece.
    fn jet<T: Real>(&self, over: bool) -> Option<Jet<T, N>> {
        let zero_v: [T; N] = std::array::from_fn(|_| T::exact_f64(0.0));
        match self {
            Piece::Line(a, b) => {
                let value = if over {
                    std::array::from_fn(|k| span(&a[k], &b[k]))
                } else {
                    std::array::from_fn(|k| c(&half(&(&a[k] + &b[k]))))
                };
                Some([value, std::array::from_fn(|k| c(&(&b[k] - &a[k]))), zero_v])
            }
            Piece::Arc {
                center,
                x,
                y,
                radius,
                a0,
                a1,
            } => {
                let theta: T = if over {
                    span(a0, a1)
                } else {
                    c(&half(&(a0 + a1)))
                };
                let (co, si) = T::cos_sin(&theta);
                let (rad, sweep) = (c::<T>(radius), c::<T>(&(a1 - a0)));
                let at = |k: usize, p: &T, q: &T| c::<T>(&x[k]).mul(p).add(&c::<T>(&y[k]).mul(q));
                let value =
                    std::array::from_fn(|k| c::<T>(&center[k]).add(&rad.mul(&at(k, &co, &si))));
                let d1 = std::array::from_fn(|k| rad.mul(&sweep).mul(&at(k, &si.neg(), &co)));
                let d2 =
                    std::array::from_fn(|k| rad.mul(&sweep.square()).mul(&at(k, &co, &si)).neg());
                Some([value, d1, d2])
            }
            Piece::Bezier(h) => {
                let lift = |x: &Vec<R>| -> Bern<T> { x.iter().map(c).collect() };
                let hs: Vec<[Bern<T>; 3]> = h
                    .iter()
                    .map(|x| {
                        let b = lift(x);
                        let d = derivative(&b);
                        let dd = derivative(&d);
                        [b, d, dd]
                    })
                    .collect();
                let at = |b: &Bern<T>| -> T {
                    if over {
                        hull(b)
                    } else {
                        super::bernstein::value(b, &T::exact_f64(0.5))
                    }
                };
                let w = [at(&hs[N][0]), at(&hs[N][1]), at(&hs[N][2])];
                if w[0].sign() != Some(Ordering::Greater) {
                    return None;
                }
                let mut value: [T; N] = zero_v.clone();
                let mut d1: [T; N] = zero_v.clone();
                let mut d2: [T; N] = zero_v;
                for k in 0..N {
                    let (hk, dk, ddk) = (at(&hs[k][0]), at(&hs[k][1]), at(&hs[k][2]));
                    let v = hk.div(&w[0])?;
                    let v1 = dk.sub(&v.mul(&w[1])).div(&w[0])?;
                    let v2 = ddk
                        .sub(&v1.mul(&w[1]).mul(&T::exact_f64(2.0)))
                        .sub(&v.mul(&w[2]))
                        .div(&w[0])?;
                    value[k] = v;
                    d1[k] = v1;
                    d2[k] = v2;
                }
                Some([value, d1, d2])
            }
        }
    }
}

/// The surface's value and five partials, `S, S_u, S_v, S_uu, S_uv,
/// S_vv`, over the box `(u, v)`.
type SurfaceJet<T> = [[T; 3]; 6];

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

/// A spline surface's jet over a box: the hull over the patches it meets,
/// each trimmed exactly to the box and enclosed by its derivative nets.
fn spline_jet<T: Real>(patches: &[ExactBezierSurface3], u: &T, v: &T) -> Option<SurfaceJet<T>> {
    let (ua, ub) = ends(u);
    let (va, vb) = ends(v);
    let mut all: Vec<SurfaceJet<T>> = Vec::new();
    for q in patches {
        let [[pu0, pu1], [pv0, pv1]] = q.domain().clone();
        let lo_u = ua.clone().max(pu0.clone());
        let hi_u = ub.clone().min(pu1.clone());
        let lo_v = va.clone().max(pv0.clone());
        let hi_v = vb.clone().min(pv1.clone());
        if lo_u > hi_u || lo_v > hi_v {
            continue;
        }
        // A degenerate box widens within the patch (a larger box is still
        // sound).
        let widen = |lo: R, hi: R, d0: &R, d1: &R| -> (R, R) {
            if lo < hi {
                return (lo, hi);
            }
            let eps = (d1 - d0) * ratio(1, 1 << 40);
            ((&lo - &eps).max(d0.clone()), (&hi + &eps).min(d1.clone()))
        };
        let (lo_u, hi_u) = widen(lo_u, hi_u, &pu0, &pu1);
        let (lo_v, hi_v) = widen(lo_v, hi_v, &pv0, &pv1);
        let sub = q.trim(&lo_u, &hi_u, &lo_v, &hi_v).ok()?;
        all.push(patch_jet(&sub)?);
    }
    if all.is_empty() {
        return None;
    }
    Some(std::array::from_fn(|i| {
        std::array::from_fn(|k| {
            let xs: Vec<T> = all.iter().map(|j| j[i][k].clone()).collect();
            hull(&xs)
        })
    }))
}

/// The enclosure of a trimmed patch's jet in global parameters.
fn patch_jet<T: Real>(q: &ExactBezierSurface3) -> Option<SurfaceJet<T>> {
    let [du, dv] = q.degrees();
    let [[u0, u1], [v0, v1]] = q.domain().clone();
    let (su, sv) = (
        c::<T>(&(exact(1) / (&u1 - &u0))),
        c::<T>(&(exact(1) / (&v1 - &v0))),
    );
    let poles = q.homogeneous_poles();
    // Hulls of a coordinate's tensor net and of its partial nets.
    let net = |k: usize| -> Vec<Vec<T>> {
        (0..=du)
            .map(|i| (0..=dv).map(|j| c(&poles[i * (dv + 1) + j][k])).collect())
            .collect()
    };
    let d_u = |a: &Vec<Vec<T>>| -> Vec<Vec<T>> {
        let m = T::exact_f64((a.len() - 1) as f64);
        if a.len() == 1 {
            return vec![vec![T::exact_f64(0.0); a[0].len()]];
        }
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
    let flat = |a: &Vec<Vec<T>>| -> T {
        let xs: Vec<T> = a.iter().flatten().cloned().collect();
        hull(&xs)
    };
    // [h, h_u, h_v, h_uu, h_uv, h_vv] for each coordinate, global units.
    let partials = |k: usize| -> [T; 6] {
        let h = net(k);
        let hu = d_u(&h);
        let hv = d_v(&h);
        [
            flat(&h),
            flat(&hu).mul(&su),
            flat(&hv).mul(&sv),
            flat(&d_u(&hu)).mul(&su.square()),
            flat(&d_v(&hu)).mul(&su).mul(&sv),
            flat(&d_v(&hv)).mul(&sv.square()),
        ]
    };
    let w = partials(3);
    if w[0].sign() != Some(Ordering::Greater) {
        return None;
    }
    let two = T::exact_f64(2.0);
    // Per coordinate: S and its five partials, by the quotient rule.
    let column = |k: usize| -> Option<[T; 6]> {
        let h = partials(k);
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
    let out: SurfaceJet<T> =
        std::array::from_fn(|i| std::array::from_fn(|k| columns[k][i].clone()));
    Some(out)
}

/// A use as pieces of its edge and pcurve, over the surface.
pub(super) struct Taylor {
    surface: Surface,
    patches: Vec<ExactBezierSurface3>,
    pieces: Vec<(Piece<3>, Piece<2>)>,
}

/// The bound's three terms on a piece, or `None` when a jet is undefined.
fn terms<T: Real>(t: &Taylor, edge: &Piece<3>, uv: &Piece<2>) -> Option<(T, T, T)> {
    let [c_mid, c1_mid, _] = edge.jet::<T>(false)?;
    let [_, _, c2] = edge.jet::<T>(true)?;
    let [p_mid, p1_mid, _] = uv.jet::<T>(false)?;
    let [p_over, p1, p2] = uv.jet::<T>(true)?;
    let jet = |u: &T, v: &T| match &t.surface {
        Surface::BSpline(_) => spline_jet(&t.patches, u, v),
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

impl Taylor {
    fn decide_in<T: Real>(&self, tol: f64, depth_limit: usize) -> Verdict {
        let tol = T::exact_f64(tol);
        let mut stack: Vec<(Piece<3>, Piece<2>, usize)> = self
            .pieces
            .iter()
            .map(|(e, p)| (e.clone(), p.clone(), 0))
            .collect();
        let mut undecided = false;
        while let Some((edge, uv, depth)) = stack.pop() {
            // An undefined jet (off the surface's domain, a weight not
            // certainly positive) leaves the piece undecided.
            let Some((d0, d1, d2)) = terms::<T>(self, &edge, &uv) else {
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
        let mut worst = 0.0_f64;
        let mut stack: Vec<(Piece<3>, Piece<2>, usize)> = self
            .pieces
            .iter()
            .map(|(e, p)| (e.clone(), p.clone(), 0))
            .collect();
        while let Some((edge, uv, depth)) = stack.pop() {
            let (d0, d1, d2) = terms::<Fast>(self, &edge, &uv)?;
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

/// A curve's pieces between the given fractions (sorted, from 0 to 1).
fn curve_pieces<const N: usize>(curve: &CurveKind<N>, cuts: &[R]) -> Option<Vec<Piece<N>>> {
    cuts.windows(2).map(|w| curve.piece(&w[0], &w[1])).collect()
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
    fn piece(&self, f0: &R, f1: &R) -> Option<Piece<N>> {
        Some(match self {
            CurveKind::Line(a, b) => {
                let at = |f: &R| -> [R; N] { std::array::from_fn(|k| &a[k] + (&b[k] - &a[k]) * f) };
                Piece::Line(at(f0), at(f1))
            }
            CurveKind::Arc {
                center,
                x,
                y,
                radius,
                start,
                sweep,
            } => Piece::Arc {
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
                Piece::Bezier(keep.into_iter().map(|k| h[k].clone()).collect())
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
            start: zero(),
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
        Curve3::BSpline(s) => {
            CurveKind::Spline(super::bernstein::spline_arcs(s.curve(), s.range())?)
        }
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
    let mut cuts: Vec<R> = vec![zero(), one.clone()];
    cuts.extend(uv.knots());
    cuts.extend(
        edge.knots()
            .into_iter()
            .map(|s| if forward { s } else { &one - s }),
    );
    cuts.sort();
    cuts.dedup();
    let uv_pieces = curve_pieces(&uv, &cuts)?;
    let edge_pieces: Vec<Piece<3>> = if forward {
        curve_pieces(&edge, &cuts)?
    } else {
        let flipped: Vec<R> = cuts.iter().rev().map(|t| &one - t).collect();
        curve_pieces(&edge, &flipped)?
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
