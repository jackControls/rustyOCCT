//! S9f.2b: spline walls against cylinder walls on crossing axes
//! (REVIEW_NOTES.md, "S9f.2b refined, before its code" and "S9f.2b.2
//! refined, before its code"); S9f.3a: against spheres, whose function
//! along a ruling is the same quadratic over the world's three rows
//! (`spline_sphere.rs`); S9f.3b: against cones, the same quadratic with the
//! radius row of negative sign (`spline_cone.rs`): every row carries its
//! sign (`terms`).
//!
//! Along a spline wall's ruling at the run parameter `tau`, `X = o + S_x(tau)
//! x + S_y(tau) y + w n` (the spline prism's exact model), the cylinder's
//! function in its own frame's exact rows (`procedural::Other`: `sum_i (g_i
//! . X - e_i)^2 - r^2`, an elliptic cylinder in the world where the stored
//! axes are not exactly orthonormal) is `F = A w^2 + 2 B(tau) w + C(tau)`:
//! with `P_i = g_i . (o + x S_x + y S_y) - e_i` and `q_i = g_i . n`, `A = sum
//! q_i^2` (positive: the axes cross), `B = sum q_i P_i` (degree `p` on each
//! Bézier arc) and `C = sum P_i^2 - r^2` (degree `2 p`). The meeting's
//! branches are `w = (-B +- sqrt(D)) / A`, `D = B^2 - A C` of degree `2 p`,
//! its roots the turning points (the ruling tangent to the cylinder).
//!
//! The meeting of a spline wall and a cylinder (`meeting`, once per pair of
//! faces) is each branch over each maximal range of the segment's run where
//! `D > 0` (`Crv::WallMeet`, a graph over `tau`, S9f.2b.1): a turning point
//! ends a range where it lies outside either face (the pieces beside it are
//! outside too); one on a face's boundary, at an interior knot or a
//! segment's end, or of multiplicity above one (the cylinder tangent to the
//! wall) is `Degenerate`. A turning point inside both faces (a loop,
//! S9f.2b.2) gets a graph over the height instead (`Window`): the run
//! parameter the one root of `F(., w)` in a rational window of the run
//! inside one arc, over the heights `[w-(tau_s), w+(tau_s)]` of a rational
//! switch `tau_s` on the side where `D > 0`, the graphs over `tau` on both
//! branches ending there and the two switch points vertices; each such piece
//! verified exactly before it is kept (`height_piece`). A point is on a
//! graph over `tau` when its profile point lies on the segment within the
//! range, it lies on the cylinder and `A w + B(tau) = sum q_i (g_i . X -
//! e_i)` (half `F`'s derivative in `w`) has the branch's sign, exactly; on
//! a graph over the height when its profile point lies on the segment
//! strictly inside the window, it lies on the cylinder and its height in the
//! range. A graph over `tau`'s midpoint at a rational `tau` lies in
//! `Q(sqrt(D(tau)))`, a graph over the height's at a rational height in the
//! window's root's field, of degree at most `2 p`.
//!
//! Vertices: a curve over a spline segment (a cap edge or a crease, `w = h0
//! + h1 S_x + h2 S_y`) meets the cylinder at the roots of `F` along it,
//! degree `2 p` (`wallcrv_cyl`); a cylinder's cap circle meets the wall at
//! the roots of `F` along its cap plane's crease on the wall, degree `2 p`,
//! at its angle there (`conic_wall`). A cap plane holding the wall's axis
//! direction meets the wall in generatrices whose points with the circle lie
//! in a tower `Q(alpha)(sqrt(delta))`: they are found instead where each
//! arc's implicit equation vanishes at the circle's projection, a
//! polynomial of degree `2 p` in its half-angle tangent `t`, each point in
//! the one field `Q(t)` (a primitive element of the tower, S9f.2b.2,
//! `tower_points`). The cylinder prism's vertical edges meet the wall by
//! S9f.1's `line_wall` (degree `p`), the spline prism's the cylinder at
//! quadratic surds. So every vertex on the wall lies in `Q(alpha)` of degree
//! at most `2 p` or in `Q(sqrt(d))`.
use super::meet::{
    conic_angle, overlap_heights, within_rounding_of_parallel, CylPair, EdgeMeet, Pos,
};
use super::model::{Affine, Crv, FaceKind, Loc, Prism, Seg, P2};
use super::num::*;
use super::procedural::{other_of, Other};
use super::spheres::Mixed;
use super::spline_walls::{
    combine, derivative, peval, peval_r, qpoint, trim, BArc, Roots, SplineSeg, WallCrv,
};
use super::turned::roots_repeated;
use crate::polynomial::real::{isolate, AlgebraicRoot, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::solid::split::{rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::Arc;

/// The quadric a spline wall meets: a cylinder (S9f.2b), a sphere
/// (S9f.3a) or a cone (S9f.3b), for the refusals' labels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Partner {
    Cylinder,
    Sphere,
    Cone,
}

impl Partner {
    fn tangent(self) -> Error {
        Error::Degenerate(match self {
            Partner::Cylinder => "a cylinder tangent to a spline wall",
            Partner::Sphere => "a sphere tangent to a spline wall",
            Partner::Cone => "a cone tangent to a spline wall",
        })
    }

    fn knot(self) -> Error {
        Error::Degenerate(match self {
            Partner::Cylinder => "a spline wall's meeting with a cylinder turning back at a knot",
            Partner::Sphere => "a spline wall's meeting with a sphere turning back at a knot",
            Partner::Cone => "a spline wall's meeting with a cone turning back at a knot",
        })
    }

    /// How far the switch lies from a loop's turning point toward where the
    /// gentler branch's slope has fallen to one: a quarter of the way for a
    /// cylinder (S9f.2b.2), a sixty-fourth for a sphere (S9f.3a: its loops'
    /// graphs over the height, as tall as the sphere's section there, cost
    /// ten times more to integrate at a quarter on the slowest fuzz variants,
    /// the graphs over the run beside them little more); a cone's as a
    /// sphere's (S9f.3b).
    fn switch(self) -> f64 {
        match self {
            Partner::Cylinder => 0.25,
            Partner::Sphere | Partner::Cone => 0.015625,
        }
    }

    pub(super) fn edge(self) -> Error {
        Error::Degenerate(match self {
            Partner::Cylinder => {
                "a spline wall's meeting with a cylinder turning back on a face's boundary"
            }
            Partner::Sphere => {
                "a spline wall's meeting with a sphere turning back on a face's boundary"
            }
            Partner::Cone => {
                "a spline wall's meeting with a cone turning back on a face's boundary"
            }
        })
    }
}

fn tangent_curve() -> Error {
    Error::Degenerate("a spline edge of one input tangent to a face of the other")
}

fn unverified() -> Error {
    Error::ComputationLimit(
        "a spline wall's meeting's graph over the height not verified (S9f.2b.2)",
    )
}

// ------------------------------------------------------------ polynomials

pub(super) fn pmul(a: &[R], b: &[R]) -> Vec<R> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if *x == zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trim(out)
}

pub(super) fn padd(a: &[R], b: &[R]) -> Vec<R> {
    let mut out = vec![zero(); a.len().max(b.len())];
    for (i, x) in a.iter().enumerate() {
        out[i] += x;
    }
    for (i, x) in b.iter().enumerate() {
        out[i] += x;
    }
    trim(out)
}

pub(super) fn pscale(a: &[R], k: &R) -> Vec<R> {
    trim(a.iter().map(|x| x * k).collect())
}

/// The other quadric as signed rows and a constant: `F = c0 + sum_i s_i (g_i
/// . X - e_i)^2` (a cylinder's two rows or a sphere's three, positive, and
/// `c0 = -r^2`; S9f.3b: a cone's two positive and its radius row `t h . X -
/// (t e_h - r)` negative, `c0 = 0`).
pub(super) fn terms(other: &Other) -> (Vec<(V, R, R)>, R) {
    let mut rows: Vec<(V, R, R)> = other
        .g
        .iter()
        .zip(&other.e)
        .map(|(g, e)| (g.clone(), e.clone(), int(1)))
        .collect();
    if other.t == zero() {
        return (rows, -(&other.r * &other.r));
    }
    rows.push((
        scale(&other.h, &other.t),
        &other.t * &other.eh - &other.r,
        int(-1),
    ));
    (rows, zero())
}

/// The other quadric's rows along a curve over an arc, `w = h0 + h1 S_x +
/// h2 S_y` (`h = None`: free, the ruling's `w` coefficient kept apart): per
/// row `P_i + q_i w`, `(P_i(s), q_i, s_i)` with its sign (`terms`).
fn rows(f: &Affine, arc: &BArc, other: &Other, h: Option<&[R; 3]>) -> Vec<(Vec<R>, R, R)> {
    terms(other)
        .0
        .into_iter()
        .map(|(g, e, sign)| {
            let q = dot(&g, &f.n);
            let (gx, gy) = (dot(&g, &f.x), dot(&g, &f.y));
            let c0 = dot(&g, &f.o) - e;
            match h {
                None => (combine(&c0, &gx, &arc.x, &gy, &arc.y), q, sign),
                Some(h) => (
                    combine(
                        &(&c0 + &q * &h[0]),
                        &(&gx + &q * &h[1]),
                        &arc.x,
                        &(&gy + &q * &h[2]),
                        &arc.y,
                    ),
                    zero(),
                    sign,
                ),
            }
        })
        .collect()
}

/// `A`, `B(s)` and `C(s)` of the other quadric's function along the arc's
/// rulings.
fn ruling_coeffs(f: &Affine, arc: &BArc, other: &Other) -> (R, Vec<R>, Vec<R>) {
    let rs = rows(f, arc, other, None);
    let a = rs.iter().fold(zero(), |acc, (_, q, s)| acc + s * q * q);
    let b = rs.iter().fold(Vec::new(), |acc: Vec<R>, (p, q, s)| {
        padd(&acc, &pscale(p, &(q * s)))
    });
    let c = rs
        .iter()
        .fold(vec![terms(other).1], |acc: Vec<R>, (p, _, s)| {
            padd(&acc, &pscale(&pmul(p, p), s))
        });
    (a, b, c)
}

/// The other quadric's function along a curve over an arc (`w` given by
/// `h`).
pub(super) fn along_curve(f: &Affine, arc: &BArc, other: &Other, h: &[R; 3]) -> Vec<R> {
    rows(f, arc, other, Some(h))
        .iter()
        .fold(vec![terms(other).1], |acc: Vec<R>, (p, _, s)| {
            padd(&acc, &pscale(&pmul(p, p), s))
        })
}

/// `A`, `B` and `C` at a rational profile point.
fn coeffs_at(f: &Affine, other: &Other, uv: &[R; 2]) -> (R, R, R) {
    let base = add(&add(&f.o, &scale(&f.x, &uv[0])), &scale(&f.y, &uv[1]));
    let (rows, c0) = terms(other);
    let (mut a, mut b, mut c) = (zero(), zero(), c0);
    for (g, e, s) in &rows {
        let (p, q) = (dot(g, &base) - e, dot(g, &f.n));
        a += s * &q * &q;
        b += s * &q * &p;
        c += s * &p * &p;
    }
    (a, b, c)
}

/// Whether a spline wall and a cylinder whose axes lie within rounding of
/// parallel are certainly apart where both faces' boxes overlap (a meeting
/// lies in both): over each arc, at the overlap's middle height `wm` in the
/// wall's frame, the wall's point in the cylinder's rows `X(s) = P(s) + q
/// wm` lies beyond the cylinder's radius, or within it, by more than `m >=
/// |q| half` (the ruling's drift over half the overlap's heights), exactly:
/// `|X|^2 - (r +- m)^2` of one sign on `[0, 1]` (no root there, isolated).
fn near_parallel_apart(
    sm: &Prism,
    sf: usize,
    seg: &SplineSeg,
    cm: &Prism,
    cf: usize,
    other: &Other,
) -> bool {
    let f = &sm.f;
    let Some((wm, half)) = overlap_heights(f, [&sm.boxes[sf], &cm.boxes[cf]]) else {
        return false;
    };
    let drift2 = other
        .g
        .iter()
        .fold(zero(), |acc, g| acc + dot(g, &f.n) * dot(g, &f.n))
        * &half
        * &half;
    // A rational bound of the drift: its binary64 root, widened, checked.
    let mut m = crate::solid::split::q(rational_f64(&drift2).sqrt() * (1.0 + 1e-6) + 1e-300);
    while &m * &m < drift2 {
        m = &m * int(2);
    }
    let (zero_r, one) = (zero(), int(1));
    // `|X|^2 - k^2` of the sign `want` all over the arc.
    let signed = |arc: &BArc, k: &R, want: Ordering| {
        let g = rows(f, arc, other, None)
            .iter()
            .fold(vec![-(k * k)], |acc: Vec<R>, (p, q, s)| {
                let x = padd(p, &[q * &wm]);
                padd(&acc, &pscale(&pmul(&x, &x), s))
            });
        let g = trim(g);
        matches!(roots_within(&g, &zero_r, &one), Some(r) if r.is_empty())
            && g.first().map_or(zero(), R::clone).cmp(&zero()) == want
    };
    let (beyond, within) = (&other.r + &m, &other.r - &m);
    seg.arcs.iter().all(|arc| {
        signed(arc, &beyond, Ordering::Greater)
            || (within > zero() && signed(arc, &within, Ordering::Less))
    })
}

/// The distinct real roots of a polynomial strictly inside `(lo, hi)`;
/// `None` for the zero polynomial or a root at either end.
fn roots_within(p: &[R], lo: &R, hi: &R) -> Option<Vec<AlgebraicRoot>> {
    let p = trim(p.to_vec());
    if p.is_empty() {
        return None;
    }
    let ip = IntPolynomial::from_rationals(&p);
    if ip.is_constant() {
        return Some(Vec::new());
    }
    if ip.sign_at(lo) == Ordering::Equal || ip.sign_at(hi) == Ordering::Equal {
        return None;
    }
    isolate(
        &ip,
        lo.clone(),
        hi.clone(),
        &mut Budget::new(RootIsolationOptions::default()),
    )
    .ok()
}

fn horner(p: &[f64], s: f64) -> f64 {
    p.iter().rev().fold(0.0, |acc, c| acc * s + c)
}

fn floats(p: &[R]) -> Vec<f64> {
    p.iter().map(rational_f64).collect()
}

// ------------------------------------------------------------ the curve

/// A graph over the height's window (S9f.2b.2): its run parameters, its
/// arc and that arc's own parameters at the window's ends.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Window {
    pub(super) tau: [R; 2],
    pub(super) arc: usize,
    pub(super) s: [R; 2],
}

/// A branch of a spline wall's meeting with a cylinder over a range of the
/// segment's run (S9f.2b.1), placed by the run parameter; or a graph over
/// the height about a turning point (S9f.2b.2, `window`), placed by the
/// height.
#[derive(Debug, Clone)]
pub(super) struct WallMeetCrv {
    pub(super) seg: Arc<SplineSeg>,
    /// The spline prism's frame (its model's).
    pub(super) f: Affine,
    pub(super) other: Other,
    /// The spline prism's operand: its face in the section is the wall.
    pub(super) carrier: usize,
    /// The branch: `A w + B(tau)` positive or negative (over the run).
    pub(super) plus: bool,
    /// Run parameters, or heights for a graph over the height.
    pub(super) range: [Qd; 2],
    pub(super) window: Option<Window>,
}

impl PartialEq for WallMeetCrv {
    fn eq(&self, o: &Self) -> bool {
        // One wall, one cylinder (a cylinder's faces split at seams share
        // it), one branch, one range.
        Arc::ptr_eq(&self.seg, &o.seg)
            && self.carrier == o.carrier
            && self.other.g == o.other.g
            && self.other.e == o.other.e
            && self.other.r == o.other.r
            && self.other.t == o.other.t
            && self.other.h == o.other.h
            && self.other.eh == o.other.eh
            && self.plus == o.plus
            && self.window == o.window
            && self.range[0].cmp(&o.range[0]) == Ordering::Equal
            && self.range[1].cmp(&o.range[1]) == Ordering::Equal
    }
}

impl WallMeetCrv {
    /// Half `F`'s derivative in `w` at a point: `sum s_i q_i (g_i . X -
    /// e_i)`.
    fn slope(&self, x: &QV) -> Qd {
        terms(&self.other)
            .0
            .iter()
            .fold(Qd::rat(zero()), |acc, (g, e, s)| {
                acc.add(
                    &qdot(x, g)
                        .add_r(&-e.clone())
                        .scale(&(dot(g, &self.f.n) * s)),
                )
            })
    }

    fn within(&self, t: &Qd) -> bool {
        t.cmp(&self.range[0]) != Ordering::Less && t.cmp(&self.range[1]) != Ordering::Greater
    }

    /// The curve's place of a point on it (its run parameter, or its
    /// height for a graph over the height), if it is on it.
    pub(super) fn locate(&self, x: &QV) -> Option<Qd> {
        let l = self.f.local_q(x);
        if let Some(win) = &self.window {
            if !self.within(&l[2]) || self.other.value(x).sign() != Ordering::Equal {
                return None;
            }
            let tau = self.seg.locate(&[l[0].clone(), l[1].clone()])?;
            let inside = tau.cmp(&Qd::rat(win.tau[0].clone())) == Ordering::Greater
                && tau.cmp(&Qd::rat(win.tau[1].clone())) == Ordering::Less;
            return inside.then(|| l[2].clone());
        }
        let tau = self.seg.locate(&[l[0].clone(), l[1].clone()])?;
        if !self.within(&tau) || self.other.value(x).sign() != Ordering::Equal {
            return None;
        }
        let want = if self.plus {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        (self.slope(x).sign() == want).then_some(tau)
    }

    pub(super) fn on(&self, x: &QV) -> bool {
        self.locate(x).is_some()
    }

    pub(super) fn place(&self, x: &QV) -> Qd {
        self.locate(x)
            .expect("a point on a spline wall's meeting at a known parameter")
    }

    /// The curve's point at a rational place inside its range: over the
    /// run in `Q(sqrt(D))`, over the height the window's root.
    pub(super) fn at(&self, k: &R) -> Option<QV> {
        if let Some(win) = &self.window {
            return self.height_at(win, k);
        }
        let uv = self.seg.point(&Qd::rat(k.clone()));
        let (Some(u), Some(v)) = (uv[0].rational(), uv[1].rational()) else {
            return None;
        };
        let uv = [u.clone(), v.clone()];
        let (a, b, c) = coeffs_at(&self.f, &self.other, &uv);
        let d = &b * &b - &a * &c;
        if d <= zero() {
            return None;
        }
        let inv = int(1) / &a;
        let s = if self.plus { int(1) } else { int(-1) };
        let w = Qd::new(-&b * &inv, s * &inv, d);
        Some(qpoint(
            &self.f,
            &Qd::rat(uv[0].clone()),
            &Qd::rat(uv[1].clone()),
            &w,
        ))
    }

    /// A graph over the height's point at a rational height: the one root
    /// of `F(., w)` in the window, its field's generator that root (S9f.1's
    /// generators: one per polynomial and root).
    fn height_at(&self, win: &Window, w: &R) -> Option<QV> {
        let arc = &self.seg.arcs[win.arc];
        let (a, b, c) = ruling_coeffs(&self.f, arc, &self.other);
        let p = padd(&padd(&pscale(&b, &(int(2) * w)), &c), &[&a * w * w]);
        let mut rs = roots_within(&p, &win.s[0], &win.s[1])?;
        if rs.len() != 1 {
            return None;
        }
        let root = rs.pop()?;
        let s = match root.rational_value() {
            Some(x) => K::Rat(x.clone()),
            None => K::generator(&Arc::new(Gen::new(p, root))),
        };
        let s = Qd::of(s);
        Some(qpoint(
            &self.f,
            &peval(&arc.x, &s),
            &peval(&arc.y, &s),
            &Qd::rat(w.clone()),
        ))
    }

    /// The tangent at a point (unit-free): `F_w dS - F_tau n`, along
    /// increasing `tau` on the plus branch (`F_w > 0`), its opposite on the
    /// other; over the height along rising `w` (its `w` component `-F_tau`,
    /// never zero on a verified window).
    pub(super) fn tangent(&self, x: &QV, place: &Qd) -> QV {
        let tau = if self.window.is_some() {
            let l = self.f.local_q(x);
            self.seg
                .locate(&[l[0].clone(), l[1].clone()])
                .expect("a point of a spline wall's meeting on its wall")
        } else {
            place.clone()
        };
        let d = self.seg.deriv(&tau);
        let ds = qadd(&qscale(&self.f.x, &d[0]), &qscale(&self.f.y, &d[1]));
        let mut ft = Qd::rat(zero());
        for (g, e, s) in &terms(&self.other).0 {
            ft = ft.add(&qdot(x, g).add_r(&-e.clone()).mul(&qdot(&ds, g)).scale(s));
        }
        let fw = self.slope(x);
        let t: QV = std::array::from_fn(|k| ds[k].mul(&fw).sub(&ft.scale(&self.f.n[k])));
        let flip = if self.window.is_some() {
            ft.sign() == Ordering::Greater
        } else {
            !self.plus
        };
        if flip {
            t.map(|c| c.neg())
        } else {
            t
        }
    }

    /// Points in binary64 from `t0` to `t1` (run parameters, or heights).
    pub(super) fn samples(&self, t0: f64, t1: f64, n: usize) -> Vec<[f64; 3]> {
        let fl = |v: &V| v.clone().map(|x| rational_f64(&x));
        let (o, x, y, nn) = (fl(&self.f.o), fl(&self.f.x), fl(&self.f.y), fl(&self.f.n));
        let (signed, c0) = terms(&self.other);
        let rows: Vec<([f64; 3], f64, f64)> = signed
            .iter()
            .map(|(g, e, s)| (fl(g), rational_f64(e), rational_f64(s)))
            .collect();
        // (`-r^2` in binary64 as before for a cylinder or a sphere.)
        let c0 = if self.other.t == zero() {
            let r = rational_f64(&self.other.r);
            -r * r
        } else {
            rational_f64(&c0)
        };
        // `A`, `B` and `C` at a profile point, and its world point at `w`.
        let abc = |u: f64, v: f64| {
            let base: [f64; 3] = [0, 1, 2].map(|k| o[k] + x[k] * u + y[k] * v);
            let (mut a, mut b, mut c) = (0.0, 0.0, c0);
            for (g, e, s) in &rows {
                let p = g[0] * base[0] + g[1] * base[1] + g[2] * base[2] - e;
                let q = g[0] * nn[0] + g[1] * nn[1] + g[2] * nn[2];
                a += s * q * q;
                b += s * q * p;
                c += s * p * p;
            }
            (base, a, b, c)
        };
        if let Some(win) = &self.window {
            let arc = &self.seg.arcs[win.arc];
            let (px, py) = (floats(&arc.x), floats(&arc.y));
            let [s0, s1] = [rational_f64(&win.s[0]), rational_f64(&win.s[1])];
            let f = |s: f64, w: f64| {
                let (_, a, b, c) = abc(horner(&px, s), horner(&py, s));
                a * w * w + 2.0 * b * w + c
            };
            return (0..=n)
                .map(|i| {
                    let w = t0 + (t1 - t0) * i as f64 / n as f64;
                    let (mut a, mut b) = (s0, s1);
                    let fa = f(a, w);
                    for _ in 0..200 {
                        let mid = 0.5 * a + 0.5 * b;
                        if mid <= a || mid >= b {
                            break;
                        }
                        if f(mid, w).signum() == fa.signum() {
                            a = mid;
                        } else {
                            b = mid;
                        }
                    }
                    let s = 0.5 * a + 0.5 * b;
                    let (base, ..) = abc(horner(&px, s), horner(&py, s));
                    [0, 1, 2].map(|k| base[k] + nn[k] * w)
                })
                .collect();
        }
        (0..=n)
            .map(|i| {
                let t = t0 + (t1 - t0) * i as f64 / n as f64;
                let [u, v] = self.seg.point_f64(t);
                let (base, a, b, c) = abc(u, v);
                let s = if self.plus { 1.0 } else { -1.0 };
                let sq = s * (b * b - a * c).max(0.0).sqrt();
                let (pp, mm) = (-b + sq, -b - sq);
                let w = if pp.abs() >= mm.abs() { pp / a } else { c / mm };
                [0, 1, 2].map(|k| base[k] + nn[k] * w)
            })
            .collect()
    }
}

// ------------------------------------------------------------ the meeting

/// A spline wall's meeting with a cylinder wall of another prism on a
/// crossing axis (`sm`'s face `sf` holding segment `seg`, `cm`'s face `cf`
/// on the cylinder `c`, `r`): each branch over each range of the run where
/// the discriminant is positive, a graph over the height about each turning
/// point inside both faces, and those graphs' switches.
#[allow(clippy::too_many_arguments)]
pub(super) fn meeting(
    sm: &Prism,
    sf: usize,
    seg: &Arc<SplineSeg>,
    cm: &Prism,
    cf: usize,
    c: &P2,
    r: &R,
) -> Result<CylPair> {
    meeting_with(
        sm,
        sf,
        seg,
        cm,
        cf,
        other_of(&cm.f, c, r),
        Partner::Cylinder,
    )
}

/// A spline wall's meeting with the quadric `other` of face `cf` of `cm`
/// (a cylinder's on a crossing axis, S9f.2b; a sphere's, S9f.3a; a cone's,
/// S9f.3b, its `A` nonzero of either sign: for `A < 0` the discriminant
/// vanishes only where a ruling passes the apex, refused before).
pub(super) fn meeting_with(
    sm: &Prism,
    sf: usize,
    seg: &Arc<SplineSeg>,
    cm: &Prism,
    cf: usize,
    other: Other,
    partner: Partner,
) -> Result<CylPair> {
    let f = &sm.f;
    // Axes within rounding of parallel (a frame's normal normalized again,
    // an ulp off the wall's): a meeting within the faces runs along a
    // sliver of the run no binary64 edge holds (its heights' turn rate past
    // 10^12 per unit of the run, the stored axes perhaps exactly parallel),
    // as two such cylinders' does (`meet::cyl_pair`). Apart where certainly
    // apart within the faces' boxes, `Degenerate` otherwise.
    if partner == Partner::Cylinder && within_rounding_of_parallel(&f.n, &cm.f.n) {
        if near_parallel_apart(sm, sf, seg, cm, cf, &other) {
            return Ok(CylPair::Apart);
        }
        return Err(Error::Degenerate(
            "a spline wall's and a cylinder's axes within rounding of parallel",
        ));
    }
    let rs = match seg.roots_of(&|arc: &BArc| {
        let (a, b, cc) = ruling_coeffs(f, arc, &other);
        padd(&pmul(&b, &b), &pscale(&cc, &-a))
    })? {
        Roots::At(rs) => rs,
        Roots::Along | Roots::Partly => return Err(partner.tangent()),
    };
    let first = Qd::rat(seg.first.clone());
    let last = Qd::rat(seg.last.clone());
    // The ranges' bounds: the segment's ends and the turning points, those
    // inside both faces marked (loops, S9f.2b.2).
    let mut bounds: Vec<(Qd, bool)> = vec![(first, false)];
    for root in &rs {
        // The turning point: on the ruling at w = -B / A.
        let uv = seg.point(&root.tau);
        let base = qpoint(f, &uv[0], &uv[1], &Qd::rat(zero()));
        let (mut a, mut b) = (zero(), Qd::rat(zero()));
        for (g, e, s) in &terms(&other).0 {
            let q = dot(g, &f.n);
            a += s * &q * &q;
            b = b.add(&qdot(&base, g).add_r(&-e.clone()).scale(&(&q * s)));
        }
        let w = b.scale(&(int(-1) / a));
        let x = qpoint(f, &uv[0], &uv[1], &w);
        let inside = (sm.in_face(sf, &x), cm.in_face(cf, &x));
        // S9f.3a: on the hemispheres' split, which is no edge of the input:
        // the split is tried at another seam.
        if partner == Partner::Sphere
            && inside.0 != Loc::Out
            && inside.1 == Loc::On
            && super::graph::seam_at(cm, cf, &x)
        {
            return Err(Error::ComputationLimit(super::graph::SEAM));
        }
        let near = !matches!(inside, (Loc::Out, _) | (_, Loc::Out));
        // Outside a face but within the resolution of both: the meeting's
        // end there a rounding away from turning back, its square root's
        // slope unbounded (a frame normalized to an ulp off a tangent
        // ruling on a cap's edge).
        let tol = sm.tolerance.linear().max(cm.tolerance.linear());
        let close = |pr: &Prism, fi: usize, at: Loc| at != Loc::Out || outside(pr, fi, &x) <= tol;
        if !near && close(sm, sf, inside.0) && close(cm, cf, inside.1) {
            return Err(partner.edge());
        }
        if near {
            if root.mult > 1 {
                return Err(partner.tangent());
            }
            if root.knot {
                return Err(partner.knot());
            }
            if root.end || inside != (Loc::In, Loc::In) {
                return Err(partner.edge());
            }
            // S9f.3a: inside both faces but within the resolution of a
            // cap's or a rim's plane: a rounding away from turning back on
            // that edge (a sphere about a turned prism's frame origin), the
            // graph over the height cut there.
            if inner_gap(sm, sf, &x) <= tol || inner_gap(cm, cf, &x) <= tol {
                return Err(partner.edge());
            }
            // A loop's turning point (S9f.2b.2): a graph over the height.
            bounds.push((root.tau.clone(), true));
            continue;
        }
        if !root.end {
            bounds.push((root.tau.clone(), false));
        }
    }
    bounds.push((last, false));
    let mut out = Vec::new();
    let mut switches = Vec::new();
    let mut switch_at: Vec<Option<R>> = vec![None; bounds.len()];
    for i in 0..bounds.len() {
        if !bounds[i].1 {
            continue;
        }
        let (piece, tau_s, xs) = height_piece(seg, f, &other, &bounds, i, partner)?;
        out.push(Crv::WallMeet(Box::new(WallMeetCrv {
            carrier: seg.op,
            ..piece
        })));
        switch_at[i] = Some(tau_s);
        switches.extend(xs);
    }
    for i in 0..bounds.len() - 1 {
        let (lo, hi) = (&bounds[i].0, &bounds[i + 1].0);
        let k = super::graph::rational_between_num(lo, hi)?;
        let uv = seg.point(&Qd::rat(k));
        let (Some(u), Some(v)) = (uv[0].rational(), uv[1].rational()) else {
            unreachable!("a rational run parameter's point is rational")
        };
        let (a, b, cc) = coeffs_at(f, &other, &[u.clone(), v.clone()]);
        if &b * &b - &a * &cc <= zero() {
            continue;
        }
        // A loop's turning point ends the ranges beside it at its switch.
        let lo = switch_at[i].clone().map_or(lo.clone(), Qd::rat);
        let hi = switch_at[i + 1].clone().map_or(hi.clone(), Qd::rat);
        for plus in [true, false] {
            out.push(Crv::WallMeet(Box::new(WallMeetCrv {
                seg: seg.clone(),
                f: f.clone(),
                other: other.clone(),
                carrier: seg.op,
                plus,
                range: [lo.clone(), hi.clone()],
                window: None,
            })));
        }
    }
    Ok(CylPair::Mixed(Box::new(Mixed {
        pieces: out,
        switches,
    })))
}

/// A binary64 view of the other quadric's function along an arc's rulings:
/// each row's `P_i(s)`, `P_i'(s)`, `q_i` and sign, `A`, the constant, the
/// arc's coordinates' derivatives and the frame's `x`, `y` and `|n|`.
struct View {
    rows: Vec<(Vec<f64>, Vec<f64>, f64, f64)>,
    a: f64,
    c0: f64,
    dx: Vec<f64>,
    dy: Vec<f64>,
    fx: [f64; 3],
    fy: [f64; 3],
    nn: f64,
}

impl View {
    fn new(f: &Affine, arc: &BArc, other: &Other) -> Self {
        let rows: Vec<(Vec<f64>, Vec<f64>, f64, f64)> = rows(f, arc, other, None)
            .iter()
            .map(|(p, q, s)| {
                (
                    floats(p),
                    floats(&derivative(p)),
                    rational_f64(q),
                    rational_f64(s),
                )
            })
            .collect();
        let a = rows.iter().map(|(_, _, q, s)| s * q * q).sum();
        let fl = |v: &V| v.clone().map(|x| rational_f64(&x));
        let n = fl(&f.n);
        Self {
            rows,
            a,
            c0: rational_f64(&terms(other).1),
            dx: floats(&arc.dx),
            dy: floats(&arc.dy),
            fx: fl(&f.x),
            fy: fl(&f.y),
            nn: (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt(),
        }
    }

    /// The gentler branch's slope at the arc's parameter `s`: `|dw/ds|`
    /// over the profile's speed (infinite where `D <= 0`).
    fn slope(&self, s: f64) -> f64 {
        let (mut b, mut db, mut c, mut dc) = (0.0, 0.0, self.c0, 0.0);
        for (p, dp, q, sg) in &self.rows {
            let (pv, dpv) = (horner(p, s), horner(dp, s));
            b += sg * q * pv;
            db += sg * q * dpv;
            c += sg * pv * pv;
            dc += 2.0 * sg * pv * dpv;
        }
        let d = b * b - self.a * c;
        if d.is_nan() || d <= 0.0 {
            return f64::INFINITY;
        }
        let dd = 2.0 * b * db - self.a * dc;
        let sq = d.sqrt();
        let (wp, wm) = (
            (-db + dd / (2.0 * sq)) / self.a,
            (-db - dd / (2.0 * sq)) / self.a,
        );
        let (sx, sy) = (horner(&self.dx, s), horner(&self.dy, s));
        let speed = (0..3)
            .map(|k| (self.fx[k] * sx + self.fy[k] * sy).powi(2))
            .sum::<f64>()
            .sqrt();
        wp.abs().min(wm.abs()) * self.nn / speed
    }
}

/// The graph over the height about the turning point `bounds[i]` (inside
/// both faces), its switch's run parameter and its two switch points (the
/// carrier set by the caller). On the side of the turning point `s*` where
/// `D > 0` (the sign of `D'` there, exactly), `d1` is where the gentler
/// branch's slope over the profile's arc length has fallen to one (binary64,
/// at most a third of the way to the next bound and inside the arc); the
/// switch `s_s` a quarter of the way for a cylinder (the slope there about
/// two: it grows as the distance's inverse square root; the graphs over the
/// height cost more per point than those over the run, the run's more per
/// piece the nearer they end to `s*`), a sixty-fourth for a sphere
/// (`Partner::switch`), the window's near end `s0` at `d1`, its far end
/// `s1` past `s*` as far (at most half the way to the bound behind and
/// inside the arc), all rational. The piece is kept when, exactly: (i) `F`'s
/// roots in `w` at `s0` (surds, or none) lie outside the range
/// `[w-(s_s), w+(s_s)]` and `D(s1) < 0`; (ii) at the rational height `-B(s_s)
/// / A` the window holds one root of `F`; (iii) `H = A C'^2 - 4 B B' C' + 4
/// C B'^2` has no root in the closed window (every root in the window is
/// then simple at every height). Otherwise the distances halve, at most
/// twenty times.
fn height_piece(
    seg: &Arc<SplineSeg>,
    f: &Affine,
    other: &Other,
    bounds: &[(Qd, bool)],
    i: usize,
    partner: Partner,
) -> Result<(WallMeetCrv, R, [QV; 2])> {
    let tau = &bounds[i].0;
    // Its arc (not at a knot: refused before).
    let k = (0..seg.arcs.len())
        .find(|&k| tau.cmp(&Qd::rat(seg.arcs[k].d[1].clone())) == Ordering::Less)
        .ok_or_else(unverified)?;
    let arc = &seg.arcs[k];
    let [d0, d1] = &arc.d;
    let span = d1 - d0;
    let s_of = |t: &Qd| t.add_r(&-d0).scale(&(int(1) / &span));
    let star = s_of(tau);
    let (a, b, c) = ruling_coeffs(f, arc, other);
    let dpoly = padd(&pmul(&b, &b), &pscale(&c, &-a.clone()));
    let side = match peval(&derivative(&dpoly), &star).sign() {
        Ordering::Greater => 1.0,
        Ordering::Less => -1.0,
        Ordering::Equal => return Err(partner.tangent()),
    };
    let s_star = star.to_f64();
    let span_f = rational_f64(&span);
    let gap = |j: usize| (bounds[j].0.to_f64() - tau.to_f64()).abs() / span_f;
    let (ahead, behind) = if side > 0.0 {
        (i + 1, i - 1)
    } else {
        (i - 1, i + 1)
    };
    let to_end = |dir: f64| if dir > 0.0 { 1.0 - s_star } else { s_star };
    let d_max = (gap(ahead) / 3.0).min(0.9 * to_end(side));
    let back_max = (gap(behind) / 2.0).min(0.9 * to_end(-side));
    if !(d_max > 0.0 && back_max > 0.0) {
        return Err(unverified());
    }
    let view = View::new(f, arc, other);
    let steep = |d: f64| view.slope(s_star + side * d) > 1.0;
    let mut reach = if steep(d_max) {
        d_max
    } else {
        let (mut lo, mut hi) = (d_max * 1e-12, d_max);
        for _ in 0..100 {
            let mid = (lo * hi).sqrt();
            if !(lo < mid && mid < hi) {
                break;
            }
            if steep(mid) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        hi
    };
    let mut back = reach.min(back_max);
    let dyadic = |x: f64| R::from_float(x).ok_or_else(unverified);
    let one = int(1);
    for _ in 0..20 {
        let ss = dyadic(s_star + side * partner.switch() * reach)?;
        let s0 = dyadic(s_star + side * reach)?;
        let s1 = dyadic(s_star - side * back)?;
        reach *= 0.5;
        back *= 0.5;
        // Strictly ordered about the turning point, inside the arc.
        let q = |x: &R| Qd::rat(x.clone());
        let ordered = if side > 0.0 {
            zero() < s1
                && q(&s1).cmp(&star) == Ordering::Less
                && star.cmp(&q(&ss)) == Ordering::Less
                && ss < s0
                && s0 < one
        } else {
            zero() < s0
                && s0 < ss
                && q(&ss).cmp(&star) == Ordering::Less
                && star.cmp(&q(&s1)) == Ordering::Less
                && s1 < one
        };
        if !ordered {
            continue;
        }
        // The range: the branches' heights at the switch.
        let at = |s: &R| {
            let (bs, cs) = (peval_r(&b, s), peval_r(&c, s));
            let d = &bs * &bs - &a * &cs;
            (bs, d)
        };
        let (bs, ds) = at(&ss);
        if ds <= zero() {
            continue;
        }
        let inv = int(1) / &a;
        let w = [
            Qd::new(-&bs * &inv, -&inv, ds.clone()),
            Qd::new(-&bs * &inv, inv.clone(), ds.clone()),
        ];
        // (i) The window's ends.
        let (b0, d0v) = at(&s0);
        let clear0 = match d0v.cmp(&zero()) {
            Ordering::Less => true,
            Ordering::Equal => false,
            Ordering::Greater => [-&inv, inv.clone()].iter().all(|sg| {
                let root = Qd::new(-&b0 * &inv, sg.clone(), d0v.clone());
                root.cmp(&w[0]) == Ordering::Less || root.cmp(&w[1]) == Ordering::Greater
            }),
        };
        if !clear0 || at(&s1).1 >= zero() {
            continue;
        }
        let (lo, hi) = if s0 < s1 {
            (s0.clone(), s1.clone())
        } else {
            (s1.clone(), s0.clone())
        };
        // (ii) One root at a height inside the range.
        let wp = -&bs * &inv;
        let probe = padd(&padd(&pscale(&b, &(int(2) * &wp)), &c), &[&a * &wp * &wp]);
        if !matches!(roots_within(&probe, &lo, &hi), Some(rs) if rs.len() == 1) {
            continue;
        }
        // (iii) No double root in the window at any height.
        let (db, dc) = (derivative(&b), derivative(&c));
        let h = padd(
            &padd(
                &pscale(&pmul(&dc, &dc), &a),
                &pscale(&pmul(&pmul(&b, &db), &dc), &int(-4)),
            ),
            &pscale(&pmul(&c, &pmul(&db, &db)), &int(4)),
        );
        if !matches!(roots_within(&h, &lo, &hi), Some(rs) if rs.is_empty()) {
            continue;
        }
        let tau_of = |s: &R| d0 + s * &span;
        let u = peval_r(&arc.x, &ss);
        let v = peval_r(&arc.y, &ss);
        let xs = [0, 1].map(|j| qpoint(f, &Qd::rat(u.clone()), &Qd::rat(v.clone()), &w[j]));
        let piece = WallMeetCrv {
            seg: seg.clone(),
            f: f.clone(),
            other: other.clone(),
            carrier: 0,
            plus: true,
            range: w,
            window: Some(Window {
                tau: [tau_of(&lo), tau_of(&hi)],
                arc: k,
                s: [lo, hi],
            }),
        };
        return Ok((piece, tau_of(&ss), xs));
    }
    Err(unverified())
}

/// How far (binary64) a point inside a face lies from the planes of its
/// caps or rims: a prism wall's heights, a sphere's ends, a cone's rims'
/// planes (S9f.3b); infinite for other faces.
fn inner_gap(pr: &Prism, fi: usize, x: &QV) -> f64 {
    match (pr.faces[fi].kind, &pr.ball) {
        (FaceKind::Wall(..), _) if pr.ball.is_none() && pr.funnel.is_none() => {
            let w = pr.f.local_q(x)[2].to_f64();
            (w - rational_f64(&pr.lo))
                .abs()
                .min((rational_f64(&pr.hi) - w).abs())
        }
        (FaceKind::Half(_), Some(ball)) => {
            let nn = rational_f64(&dot(&ball.n, &ball.n));
            let along = qdot(&qsub(x, &qv(&ball.c)), &ball.n).to_f64();
            ball.ends
                .iter()
                .flatten()
                .map(|h| (along - rational_f64(h) * nn).abs() / nn.sqrt())
                .fold(f64::INFINITY, f64::min)
        }
        // S9f.3b: a cone's wall, from its rims' planes.
        (FaceKind::ConeWall, _) => {
            let w = pr.f.local_q(x)[2].to_f64();
            w.abs().min((rational_f64(&pr.hi) - w).abs())
        }
        _ => f64::INFINITY,
    }
}

/// How far (binary64, in the prism's local units) a point on a wall face's
/// surface lies outside the face: past its heights and, where its generatrix
/// leaves the face's line or arc, from the nearer of their ends.
fn outside(pr: &Prism, fi: usize, x: &QV) -> f64 {
    let FaceKind::Wall(b, j) = pr.faces[fi].kind else {
        return f64::INFINITY;
    };
    let l = pr.f.local_q(x);
    let w = if l[2].cmp(&Qd::rat(pr.lo.clone())) == Ordering::Less {
        Qd::rat(pr.lo.clone())
    } else if l[2].cmp(&Qd::rat(pr.hi.clone())) == Ordering::Greater {
        Qd::rat(pr.hi.clone())
    } else {
        l[2].clone()
    };
    let height = (l[2].to_f64() - w.to_f64()).abs();
    let across = if pr.in_face(fi, &qpoint(&pr.f, &l[0], &l[1], &w)) == Loc::Out {
        let at = [l[0].to_f64(), l[1].to_f64()];
        let gap = |e: &P2| (at[0] - rational_f64(&e[0])).hypot(at[1] - rational_f64(&e[1]));
        let seg = &pr.bounds[b].segs[j];
        if matches!(seg, Seg::Spline(_)) {
            0.0
        } else {
            gap(seg.start()).min(gap(seg.end()))
        }
    } else {
        0.0
    };
    height.hypot(across)
}

// ------------------------------------------------------------ the vertices

/// Where a curve over a spline segment (a cap edge or a crease) meets a
/// cylinder wall of another prism on a crossing axis (`c`, `r` on frame
/// `cf`): the roots of the cylinder's function along it.
pub(super) fn wallcrv_cyl(curve: &WallCrv, cf: &Affine, c: &P2, r: &R) -> Result<EdgeMeet> {
    wallcrv_quadric(curve, &other_of(cf, c, r))
}

/// Where a curve over a spline segment meets the quadric `other` (a
/// cylinder's, S9f.2b; a sphere's, S9f.3a): the roots of its function
/// along it, degree `2 p`.
pub(super) fn wallcrv_quadric(curve: &WallCrv, other: &Other) -> Result<EdgeMeet> {
    let rs = match curve
        .seg
        .roots_of(&|arc: &BArc| along_curve(&curve.f, arc, other, &curve.h))?
    {
        Roots::At(rs) => rs,
        Roots::Along => return Ok(EdgeMeet::Along),
        Roots::Partly => return Err(tangent_curve()),
    };
    let mut out = Vec::new();
    for root in rs {
        if root.mult > 1 {
            return Err(tangent_curve());
        }
        let x = curve.point(&root.tau);
        out.push((Pos::T(root.tau), x));
    }
    Ok(EdgeMeet::Points(out))
}

/// Where a cylinder's cap edge (the conic `cc + a cos + b sin` on its own
/// cylinder `c`, `r` of the model `cm`) meets a spline wall (segment `seg` on
/// frame `sf`) on a crossing axis: the cylinder's function along the cap
/// plane's crease on the wall, each root at its angle on the circle; a cap
/// plane holding the wall's axis direction by the arcs' implicit equations
/// along the circle (`tower_points`).
#[allow(clippy::too_many_arguments)]
pub(super) fn conic_wall(
    cm: &Prism,
    c: &P2,
    r: &R,
    cc: &V,
    a: &V,
    b: &V,
    sf: &Affine,
    seg: &SplineSeg,
) -> Result<EdgeMeet> {
    let m = cross(a, b);
    let mn = dot(&m, &sf.n);
    let k = dot(&m, &sub(&sf.o, cc));
    let (mx, my) = (dot(&m, &sf.x), dot(&m, &sf.y));
    if mn == zero() {
        // The cap plane holds the wall's axis direction: its generatrices
        // on the wall, the circle's points there in a tower (S9f.2b.2);
        // none when the plane misses the segment.
        return match seg.roots_of(&|arc: &BArc| combine(&k, &mx, &arc.x, &my, &arc.y))? {
            Roots::At(rs) if rs.is_empty() => Ok(EdgeMeet::None),
            Roots::At(_) => tower_points(cc, a, b, sf, seg, Partner::Cylinder),
            Roots::Along | Roots::Partly => {
                Err(Error::OutOfDomain("a spline wall along a plane (S9f)"))
            }
        };
    }
    conic_crease(&other_of(&cm.f, c, r), cc, a, b, sf, seg)
}

/// Where a conic `cc + a cos + b sin` meets a spline wall on its plane's
/// crease, the plane not holding the wall's axis: the roots of the quadric
/// `other`'s function along the crease (a quadric holding the conic: its
/// cylinder's, S9f.2b; S9f.3b: a cone's rim's own elliptic cylinder),
/// each placed by its angle on the conic.
pub(super) fn conic_crease(
    other: &Other,
    cc: &V,
    a: &V,
    b: &V,
    sf: &Affine,
    seg: &SplineSeg,
) -> Result<EdgeMeet> {
    let m = cross(a, b);
    let mn = dot(&m, &sf.n);
    let k = dot(&m, &sub(&sf.o, cc));
    let (mx, my) = (dot(&m, &sf.x), dot(&m, &sf.y));
    let inv = int(-1) / &mn;
    let h = [&k * &inv, &mx * &inv, &my * &inv];
    let rs = match seg.roots_of(&|arc: &BArc| along_curve(sf, arc, other, &h))? {
        Roots::At(rs) => rs,
        Roots::Along | Roots::Partly => return Err(tangent_curve()),
    };
    let mut out = Vec::new();
    for root in rs {
        if root.mult > 1 {
            return Err(tangent_curve());
        }
        let uv = seg.point(&root.tau);
        let w = uv[0].scale(&h[1]).add(&uv[1].scale(&h[2])).add_r(&h[0]);
        let x = qpoint(sf, &uv[0], &uv[1], &w);
        let cs = conic_angle(cc, a, b, &x);
        out.push((Pos::Ang(cs), x));
    }
    Ok(EdgeMeet::Points(out))
}

/// The rotations tried in turn for a chart of the circle (rational points
/// of the unit circle, S9d.2c's).
fn rotations() -> Vec<[R; 2]> {
    let r = |a: i64, b: i64, c: i64| [R::new(a.into(), c.into()), R::new(b.into(), c.into())];
    vec![
        r(1, 0, 1),
        r(3, 4, 5),
        r(5, 12, 13),
        r(8, 15, 17),
        r(7, 24, 25),
        r(20, 21, 29),
        r(12, 35, 37),
        r(9, 40, 41),
    ]
}

/// S9f.2b.2: where the circle `cc + a cos + b sin`, in a plane holding the
/// wall's axis direction, meets the spline wall (segment `seg` on frame
/// `sf`): its projection along the axis meets the profile where each arc's
/// implicit equation vanishes, `f(u(t), v(t)) (1 + t^2)^d`, a polynomial of
/// degree at most `2 p` in the circle's half-angle tangent `t` in a chart
/// whose antipode is no point of the curve; each real root whose point lies
/// on the segment (`locate`: not another branch of the implicit curve) is a
/// crossing, its point and its `(cos, sin)` in `Q(t)` (a primitive element
/// of the tower `Q(alpha)(sqrt(delta))`); a point found on two arcs (a
/// knot) is kept once. A repeated root on the segment is the circle tangent
/// to a generatrix there: the meeting's turning point on the cap's rim.
pub(super) fn tower_points(
    cc: &V,
    a: &V,
    b: &V,
    sf: &Affine,
    seg: &SplineSeg,
    partner: Partner,
) -> Result<EdgeMeet> {
    let (l0, la, lb) = (sf.local(cc), sf.local_dir(a), sf.local_dir(b));
    let mut out: Vec<(Pos, QV)> = Vec::new();
    for arc in &seg.arcs {
        let imp = arc.implicit();
        let mut done = false;
        for [c0, s0] in rotations() {
            // The antipode `(cos, sin) = -(c0, s0)`: a rational point.
            let anti = [0, 1].map(|i| Qd::rat(&l0[i] - &la[i] * &c0 - &lb[i] * &s0));
            if imp.value(&anti).sign() == Ordering::Equal {
                continue;
            }
            // cos = (c0 (1 - t^2) - 2 s0 t) / q, sin = (s0 (1 - t^2) + 2 c0 t)
            // / q, q = 1 + t^2.
            let q = vec![int(1), zero(), int(1)];
            let cn = vec![c0.clone(), &s0 * int(-2), -c0.clone()];
            let sn = vec![s0.clone(), &c0 * int(2), -s0.clone()];
            let coord = |i: usize| {
                padd(
                    &padd(&pscale(&q, &l0[i]), &pscale(&cn, &la[i])),
                    &pscale(&sn, &lb[i]),
                )
            };
            let poly = imp.homogeneous(&coord(0), &coord(1), &q);
            if poly.is_empty() {
                return Err(Error::OutOfDomain("a spline wall along a plane (S9f)"));
            }
            done = true;
            if poly.len() <= 1 {
                break;
            }
            let (sf_poly, roots) = roots_repeated(&poly)?;
            for (root, repeated) in roots {
                let t = match root.rational_value() {
                    Some(x) => K::Rat(x.clone()),
                    None => K::generator(&Arc::new(Gen::new(sf_poly.clone(), root))),
                };
                let inv = t
                    .mul(&t)
                    .add(&K::Rat(int(1)))
                    .recip()
                    .expect("1 + t^2 is positive");
                let (cn, sn) = (K::Rat(int(1)).sub(&t.mul(&t)), t.scale(&int(2)));
                let cs = [
                    Qd::of(cn.scale(&c0).sub(&sn.scale(&s0)).mul(&inv)),
                    Qd::of(cn.scale(&s0).add(&sn.scale(&c0)).mul(&inv)),
                ];
                let x: QV = std::array::from_fn(|i| {
                    cs[0].scale(&a[i]).add(&cs[1].scale(&b[i])).add_r(&cc[i])
                });
                let l = sf.local_q(&x);
                if seg.locate(&[l[0].clone(), l[1].clone()]).is_none() {
                    continue;
                }
                if repeated {
                    return Err(partner.edge());
                }
                if !out.iter().any(|(_, y)| qv_eq(y, &x)) {
                    out.push((Pos::Ang(cs), x));
                }
            }
            break;
        }
        if !done {
            return Err(Error::ComputationLimit(
                "a cap circle's chart through the spline wall (S9f.2b.2)",
            ));
        }
    }
    Ok(EdgeMeet::Points(out))
}
