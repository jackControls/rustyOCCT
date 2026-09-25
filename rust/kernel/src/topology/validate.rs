//! Complete, typed validation of analytic boundary representations.
//!
//! Combinatorial decisions are exact. Geometric decisions are certified with
//! rational intervals: a check passes only when a certified upper bound is
//! within tolerance, fails only when a certified lower bound exceeds it, and
//! otherwise reports an explicit `Uncertified*` issue. Nothing is sampled to
//! decide validity. Reference: OCCT `BRepCheck_Analyzer` and its per-shape
//! checkers (see SOURCE_MAP.md); OCCT samples curve/pcurve agreement instead.
//!
//! Conventions: every curve and pcurve uses a normalized fraction in [0,1].
//! A use's pcurve follows the oriented face, so a reversed use pairs pcurve
//! fraction t with edge fraction 1-t. Outer loops wind counter-clockwise about
//! the oriented face normal; inner loops wind clockwise.
use super::{Coedge, Curve2, Curve3, Edge, Face, FaceId, Orientation, Surface, Vertex};
use crate::certified::{pi, Fast, Interval as I, Real};
use crate::{Frame3, Tolerance};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;
use std::fmt;

/// Issue classes of the validation contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IssueKind {
    Reference,
    UnusedVertex,
    UnusedEdge,
    FaceWithoutShell,
    FaceReused,
    EmptyFace,
    EmptyLoop,
    EmptyShell,
    OpenLoop,
    FreeEdge,
    NonManifoldEdge,
    SameSenseUses,
    EdgeAcrossShells,
    DisconnectedShell,
    NonManifoldVertex,
    Euler,
    DegenerateVertex,
    DegenerateCurve,
    DegenerateSurface,
    DegeneratePcurve,
    VertexOffCurve,
    PcurveOffEdge,
    UvGap,
    LoopWinding,
    InnerLoopOutside,
    ShellOrientation,
    CavityOutside,
    NestedCavity,
    UncertifiedVertexOffCurve,
    UncertifiedPcurveOffEdge,
    UncertifiedUvGap,
    UncertifiedLoopWinding,
    UncertifiedContainment,
    UncertifiedShellOrientation,
}

impl IssueKind {
    pub fn name(self) -> &'static str {
        use IssueKind::*;
        match self {
            Reference => "reference",
            UnusedVertex => "unused_vertex",
            UnusedEdge => "unused_edge",
            FaceWithoutShell => "face_without_shell",
            FaceReused => "face_reused",
            EmptyFace => "empty_face",
            EmptyLoop => "empty_loop",
            EmptyShell => "empty_shell",
            OpenLoop => "open_loop",
            FreeEdge => "free_edge",
            NonManifoldEdge => "non_manifold_edge",
            SameSenseUses => "same_sense_uses",
            EdgeAcrossShells => "edge_across_shells",
            DisconnectedShell => "disconnected_shell",
            NonManifoldVertex => "non_manifold_vertex",
            Euler => "euler",
            DegenerateVertex => "degenerate_vertex",
            DegenerateCurve => "degenerate_curve",
            DegenerateSurface => "degenerate_surface",
            DegeneratePcurve => "degenerate_pcurve",
            VertexOffCurve => "vertex_off_curve",
            PcurveOffEdge => "pcurve_off_edge",
            UvGap => "uv_gap",
            LoopWinding => "loop_winding",
            InnerLoopOutside => "inner_loop_outside",
            ShellOrientation => "shell_orientation",
            CavityOutside => "cavity_outside",
            NestedCavity => "nested_cavity",
            UncertifiedVertexOffCurve => "uncertified_vertex_off_curve",
            UncertifiedPcurveOffEdge => "uncertified_pcurve_off_edge",
            UncertifiedUvGap => "uncertified_uv_gap",
            UncertifiedLoopWinding => "uncertified_loop_winding",
            UncertifiedContainment => "uncertified_containment",
            UncertifiedShellOrientation => "uncertified_shell_orientation",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EdgeEnd {
    Start,
    End,
}

/// The offending entity; loops and uses are indexed within their face.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Entity {
    Vertex(usize),
    Edge(usize),
    EdgeEnd(usize, EdgeEnd),
    Face(usize),
    Loop(usize, usize),
    Use(usize, usize, usize),
    Shell(usize),
}

impl fmt::Display for Entity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Vertex(v) => write!(f, "vertex {v}"),
            Self::Edge(e) => write!(f, "edge {e}"),
            Self::EdgeEnd(e, EdgeEnd::Start) => write!(f, "edge {e} start"),
            Self::EdgeEnd(e, EdgeEnd::End) => write!(f, "edge {e} end"),
            Self::Face(i) => write!(f, "face {i}"),
            Self::Loop(i, l) => write!(f, "loop {i}.{l}"),
            Self::Use(i, l, u) => write!(f, "use {i}.{l}.{u}"),
            Self::Shell(s) => write!(f, "shell {s}"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Issue {
    pub kind: IssueKind,
    pub entity: Entity,
}

impl fmt::Display for Issue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.kind.name(), self.entity)
    }
}

// ------------------------------------------------------------------ helpers
//
// Geometry is generic over `Real`. Checks first run with outward-rounded
// binary64 intervals (`Fast`); any undecided result is recomputed with exact
// rational intervals. Both tiers only ever return enclosures.

fn r(x: f64) -> R {
    R::from_float(x).expect("finite value checked before exact conversion")
}
fn int(n: i64) -> R {
    R::from_integer(BigInt::from(n))
}
type V3<T> = [T; 3];
type V2<T> = [T; 2];

fn c<T: Real>(x: f64) -> T {
    T::exact_f64(x)
}
fn q<T: Real>(x: &R) -> T {
    T::from_r(x)
}
fn v3<T: Real>(p: [f64; 3]) -> V3<T> {
    p.map(|x| c(x))
}
fn vadd<T: Real>(a: &V3<T>, b: &V3<T>) -> V3<T> {
    std::array::from_fn(|i| a[i].add(&b[i]))
}
fn vsub<T: Real>(a: &V3<T>, b: &V3<T>) -> V3<T> {
    std::array::from_fn(|i| a[i].sub(&b[i]))
}
fn vscale<T: Real>(a: &V3<T>, s: &T) -> V3<T> {
    std::array::from_fn(|i| a[i].mul(s))
}
fn vdot<T: Real>(a: &V3<T>, b: &V3<T>) -> T {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}
fn vcross<T: Real>(a: &V3<T>, b: &V3<T>) -> V3<T> {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}
fn zero3<T: Real>() -> V3<T> {
    std::array::from_fn(|_| c(0.0))
}

/// Stored frame vectors: they define the represented geometry exactly.
struct FrameV<T> {
    o: V3<T>,
    x: V3<T>,
    y: V3<T>,
    n: V3<T>,
}
fn frame<T: Real>(f: &Frame3) -> FrameV<T> {
    FrameV {
        o: v3(f.origin().to_array()),
        x: v3(f.x().to_array()),
        y: v3(f.y().to_array()),
        n: v3(f.normal().to_array()),
    }
}

/// cos and sin over an interval, from the midpoint and |d cos| <= |d angle|.
fn cos_sin_i<T: Real>(u: &T) -> (T, T) {
    T::cos_sin(u)
}

/// (start, sweep) of an arc; a full circle is the arc from 0 with sweep TAU.
fn arc_of(c: &Curve3) -> Option<(&Frame3, f64, f64, f64)> {
    match c {
        Curve3::Circle { frame, radius } => Some((frame, *radius, 0.0, TAU)),
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => Some((frame, *radius, *start_angle, *sweep_angle)),
        Curve3::LineSegment { .. } => None,
    }
}

fn curve_at<T: Real>(curve: &Curve3, t: f64) -> V3<T> {
    match curve {
        Curve3::LineSegment { start, .. } if t == 0.0 => v3::<T>(start.to_array()),
        Curve3::LineSegment { end, .. } if t == 1.0 => v3::<T>(end.to_array()),
        Curve3::LineSegment { start, end } => {
            let (a, b) = (v3::<T>(start.to_array()), v3::<T>(end.to_array()));
            vadd(&a, &vscale(&vsub(&b, &a), &c(t)))
        }
        _ => {
            let (f, radius, start, sweep) = arc_of(curve).unwrap();
            let fr = frame::<T>(f);
            let (co, si) = T::cos_sin(&c::<T>(start).add(&c::<T>(sweep).mul(&c(t))));
            let rad = c::<T>(radius);
            vadd(
                &fr.o,
                &vadd(&vscale(&fr.x, &rad.mul(&co)), &vscale(&fr.y, &rad.mul(&si))),
            )
        }
    }
}

fn pcurve_at<T: Real>(p: &Curve2, t: f64) -> V2<T> {
    match p {
        Curve2::LineSegment { start, .. } if t == 0.0 => [c(start.x), c(start.y)],
        Curve2::LineSegment { end, .. } if t == 1.0 => [c(end.x), c(end.y)],
        Curve2::LineSegment { start, end } => {
            let tt = c::<T>(t);
            [
                c::<T>(start.x).add(&c::<T>(end.x).sub(&c(start.x)).mul(&tt)),
                c::<T>(start.y).add(&c::<T>(end.y).sub(&c(start.y)).mul(&tt)),
            ]
        }
        Curve2::CircularArc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let angle = c::<T>(*start_angle).add(&c::<T>(*sweep_angle).mul(&c(t)));
            let (co, si) = T::cos_sin(&angle);
            let rad = c::<T>(*radius);
            [
                c::<T>(center.x).add(&rad.mul(&co)),
                c::<T>(center.y).add(&rad.mul(&si)),
            ]
        }
    }
}

fn surface_at<T: Real>(s: &Surface, uv: &V2<T>) -> V3<T> {
    match s {
        Surface::Plane(f) => {
            let fr = frame::<T>(f);
            vadd(&fr.o, &vadd(&vscale(&fr.x, &uv[0]), &vscale(&fr.y, &uv[1])))
        }
        Surface::Cylinder { frame: f, radius } => {
            let fr = frame::<T>(f);
            let (co, si) = cos_sin_i(&uv[0]);
            let rad = c::<T>(*radius);
            let radial = vadd(&vscale(&fr.x, &rad.mul(&co)), &vscale(&fr.y, &rad.mul(&si)));
            vadd(&vadd(&fr.o, &radial), &vscale(&fr.n, &uv[1]))
        }
    }
}

/// Certified three-valued comparison of a squared quantity with tol^2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    Within,
    Beyond,
    Unknown,
}
fn within<T: Real>(d2: &T, tol2: &T) -> Verdict {
    match d2.cmp(tol2) {
        Some(Ordering::Greater) => Verdict::Beyond,
        Some(_) => Verdict::Within,
        // Undecided only if tol2 lies strictly inside: within iff hi <= tol2.
        None => Verdict::Unknown,
    }
}

/// Run a certified decision with binary64 intervals, then exactly if needed.
fn tiered<X: PartialEq>(unknown: X, fast: impl FnOnce() -> X, exact: impl FnOnce() -> X) -> X {
    let first = fast();
    if first == unknown {
        exact()
    } else {
        first
    }
}

/// A(t) = a0 + a1 t + sum_w c_w cos(w t) + s_w sin(w t), frequencies exact.
struct Harmonic<T> {
    a0: V3<T>,
    a1: V3<T>,
    terms: Vec<(R, V3<T>, V3<T>)>,
}
impl<T: Real> Harmonic<T> {
    fn new() -> Self {
        Self {
            a0: zero3(),
            a1: zero3(),
            terms: Vec::new(),
        }
    }
    fn affine(&mut self, a0: &V3<T>, a1: &V3<T>, plus: bool) {
        let s = c::<T>(if plus { 1.0 } else { -1.0 });
        self.a0 = vadd(&self.a0, &vscale(a0, &s));
        self.a1 = vadd(&self.a1, &vscale(a1, &s));
    }
    /// plus/minus (cos(alpha + omega t) cv + sin(alpha + omega t) sv).
    fn rotating(&mut self, alpha: &T, omega: &R, cv: &V3<T>, sv: &V3<T>, plus: bool) {
        let (ca, sa) = T::cos_sin(alpha);
        let sign = omega.numer().sign();
        if sign == Sign::NoSign {
            let value = vadd(&vscale(cv, &ca), &vscale(sv, &sa));
            self.affine(&value, &zero3(), plus);
            return;
        }
        let g = c::<T>(if sign == Sign::Plus { 1.0 } else { -1.0 });
        let w = if sign == Sign::Plus {
            omega.clone()
        } else {
            -omega
        };
        let cterm = vadd(&vscale(cv, &ca), &vscale(sv, &sa));
        let sterm = vadd(&vscale(cv, &sa.neg().mul(&g)), &vscale(sv, &ca.mul(&g)));
        let s = c::<T>(if plus { 1.0 } else { -1.0 });
        let (cterm, sterm) = (vscale(&cterm, &s), vscale(&sterm, &s));
        if let Some(term) = self.terms.iter_mut().find(|t| t.0 == w) {
            term.1 = vadd(&term.1, &cterm);
            term.2 = vadd(&term.2, &sterm);
        } else {
            self.terms.push((w, cterm, sterm));
        }
    }
    /// Enclosure of an upper bound of sup_t |A(t)| over [0,1].
    fn upper(&self) -> T {
        let at0 = vdot(&self.a0, &self.a0).sqrt();
        let end = vadd(&self.a0, &self.a1);
        let at1 = vdot(&end, &end).sqrt();
        // max(at0, at1) <= at0 + at1 is looser; use the sum of both ends'
        // squares' root instead: max(a,b) <= sqrt(a^2 + b^2).
        let mut bound = at0.square().add(&at1.square()).sqrt();
        for (_, cv, sv) in &self.terms {
            let (cc, ss, cs) = (vdot(cv, cv), vdot(sv, sv), vdot(cv, sv));
            let root = cc.sub(&ss).square().add(&cs.square().mul(&c(4.0))).sqrt();
            let lambda = cc.add(&ss).add(&root).mul(&c(0.5));
            bound = bound.add(&lambda.sqrt());
        }
        bound
    }
}

fn add_curve<T: Real>(h: &mut Harmonic<T>, curve: &Curve3, forward: bool) {
    match curve {
        Curve3::LineSegment { start, end } => {
            let (a, b) = (v3::<T>(start.to_array()), v3::<T>(end.to_array()));
            if forward {
                h.affine(&a, &vsub(&b, &a), true);
            } else {
                h.affine(&b, &vsub(&a, &b), true);
            }
        }
        _ => {
            let (f, radius, start, sweep) = arc_of(curve).unwrap();
            let fr = frame::<T>(f);
            let (mut alpha, mut omega) = (c::<T>(start), r(sweep));
            if !forward {
                alpha = alpha.add(&c(sweep));
                omega = -omega;
            }
            let rad = c::<T>(radius);
            h.affine(&fr.o, &zero3(), true);
            h.rotating(
                &alpha,
                &omega,
                &vscale(&fr.x, &rad),
                &vscale(&fr.y, &rad),
                true,
            );
        }
    }
}

/// Subtract S(P(t)); false when the composition is not harmonic.
fn sub_use<T: Real>(h: &mut Harmonic<T>, s: &Surface, p: &Curve2) -> bool {
    match (s, p) {
        (Surface::Plane(f), Curve2::LineSegment { start, end }) => {
            let fr = frame::<T>(f);
            let base = vadd(
                &fr.o,
                &vadd(&vscale(&fr.x, &c(start.x)), &vscale(&fr.y, &c(start.y))),
            );
            let du = c::<T>(end.x).sub(&c(start.x));
            let dv = c::<T>(end.y).sub(&c(start.y));
            let slope = vadd(&vscale(&fr.x, &du), &vscale(&fr.y, &dv));
            h.affine(&base, &slope, false);
            true
        }
        (
            Surface::Plane(f),
            Curve2::CircularArc {
                center,
                radius,
                start_angle,
                sweep_angle,
            },
        ) => {
            let fr = frame::<T>(f);
            let base = vadd(
                &fr.o,
                &vadd(&vscale(&fr.x, &c(center.x)), &vscale(&fr.y, &c(center.y))),
            );
            h.affine(&base, &zero3(), false);
            let rad = c::<T>(*radius);
            h.rotating(
                &c(*start_angle),
                &r(*sweep_angle),
                &vscale(&fr.x, &rad),
                &vscale(&fr.y, &rad),
                false,
            );
            true
        }
        (Surface::Cylinder { frame: f, radius }, Curve2::LineSegment { start, end }) => {
            let fr = frame::<T>(f);
            let dv = c::<T>(end.y).sub(&c(start.y));
            h.affine(
                &vadd(&fr.o, &vscale(&fr.n, &c(start.y))),
                &vscale(&fr.n, &dv),
                false,
            );
            let rad = c::<T>(*radius);
            let du = r(end.x) - r(start.x);
            h.rotating(
                &c(start.x),
                &du,
                &vscale(&fr.x, &rad),
                &vscale(&fr.y, &rad),
                false,
            );
            true
        }
        (Surface::Cylinder { .. }, Curve2::CircularArc { .. }) => false,
    }
}

const SAMPLES: i64 = 32;

fn deviation<T: Real>(
    curve: &Curve3,
    s: &Surface,
    p: &Curve2,
    forward: bool,
    tol: &T,
    tol2: &T,
) -> Verdict {
    let mut h = Harmonic::<T>::new();
    add_curve(&mut h, curve, forward);
    if sub_use(&mut h, s, p) && matches!(h.upper().cmp(tol), Some(Ordering::Less | Ordering::Equal))
    {
        return Verdict::Within;
    }
    for k in 0..=SAMPLES {
        let t = k as f64 / SAMPLES as f64;
        let tc = if forward { t } else { 1.0 - t };
        let d = vsub(
            &curve_at::<T>(curve, tc),
            &surface_at(s, &pcurve_at::<T>(p, t)),
        );
        if vdot(&d, &d).cmp(tol2) == Some(Ordering::Greater) {
            return Verdict::Beyond;
        }
    }
    Verdict::Unknown
}

fn finite(values: &[f64]) -> bool {
    values.iter().all(|x| x.is_finite())
}

/// Whether |end - start|^2 > tol2, or None when this tier cannot decide.
fn longer<T: Real>(start: [f64; 3], end: [f64; 3], tol2: &T) -> Option<bool> {
    let d = vsub(&v3::<T>(end), &v3::<T>(start));
    Some(vdot(&d, &d).cmp(tol2)? == Ordering::Greater)
}

fn curve_valid(curve: &Curve3, tol: &R, fast_tol2: &Fast, exact_tol2: &I) -> bool {
    match curve {
        Curve3::LineSegment { start, end } => {
            let (a, b) = (start.to_array(), end.to_array());
            if !finite(&[a, b].concat()) {
                return false;
            }
            tiered(
                None,
                || longer(a, b, fast_tol2),
                || longer(a, b, exact_tol2),
            )
            .unwrap_or_else(|| {
                // Below the interval grid: compare the exact rationals.
                let d: [R; 3] = std::array::from_fn(|i| r(b[i]) - r(a[i]));
                d.iter().map(|x| x * x).sum::<R>() > tol * tol
            })
        }
        _ => {
            let (_, radius, start, sweep) = arc_of(curve).unwrap();
            finite(&[radius, start, sweep])
                && r(radius) > *tol
                && sweep != 0.0
                && sweep.abs() <= TAU
        }
    }
}

fn surface_valid(s: &Surface, tol: &R) -> bool {
    match s {
        Surface::Plane(_) => true,
        Surface::Cylinder { radius, .. } => radius.is_finite() && r(*radius) > *tol,
    }
}

fn pcurve_valid(p: &Curve2) -> bool {
    match p {
        Curve2::LineSegment { start, end } => {
            finite(&[start.x, start.y, end.x, end.y]) && (start.x, start.y) != (end.x, end.y)
        }
        Curve2::CircularArc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            finite(&[center.x, center.y, *radius, *start_angle, *sweep_angle])
                && *radius > 0.0
                && *sweep_angle != 0.0
                && sweep_angle.abs() <= TAU
        }
    }
}

fn use_vertices(edges: &[Edge], u: &Coedge) -> (usize, usize) {
    let edge = &edges[u.edge.0];
    if u.orientation == Orientation::Forward {
        (edge.start.0, edge.end.0)
    } else {
        (edge.end.0, edge.start.0)
    }
}

fn vertex_gap<T: Real>(curve: &Curve3, vertex: [f64; 3], t: f64, tol2: &T) -> Verdict {
    let d = vsub(&curve_at::<T>(curve, t), &v3::<T>(vertex));
    within(&vdot(&d, &d), tol2)
}

fn uv_gap<T: Real>(s: &Surface, p: &Curve2, next: &Curve2, tol2: &T) -> Verdict {
    let (a, b) = (pcurve_at::<T>(p, 1.0), pcurve_at::<T>(next, 0.0));
    let mut du = a[0].sub(&b[0]);
    if let Surface::Cylinder { radius, .. } = s {
        du = du.mul(&c(*radius));
    }
    let dv = a[1].sub(&b[1]);
    within(&du.square().add(&dv.square()), tol2)
}

// ------------------------------------------------------------------ UV geometry

/// Twice the signed area contribution: integral of (u dv - v du).
fn area_term<T: Real>(p: &Curve2) -> T {
    match p {
        Curve2::LineSegment { start, end } => c::<T>(start.x)
            .mul(&c(end.y))
            .sub(&c::<T>(end.x).mul(&c(start.y))),
        Curve2::CircularArc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let a0 = c::<T>(*start_angle);
            let a1 = a0.add(&c(*sweep_angle));
            let ((c0, s0), (c1, s1)) = (T::cos_sin(&a0), T::cos_sin(&a1));
            let rho = c::<T>(*radius);
            let linear = c::<T>(center.x)
                .mul(&s1.sub(&s0))
                .sub(&c::<T>(center.y).mul(&c1.sub(&c0)));
            rho.mul(&linear).add(&rho.square().mul(&c(*sweep_angle)))
        }
    }
}

/// Twice the signed area of a loop closed by the same straight segments
/// between consecutive uses as [`crossings`].
fn loop_area<T: Real>(lp: &[Coedge]) -> T {
    lp.iter().enumerate().fold(c::<T>(0.0), |acc, (k, u)| {
        let (a, b) = (
            pcurve_at::<T>(&u.pcurve, 1.0),
            pcurve_at::<T>(&lp[(k + 1) % lp.len()].pcurve, 0.0),
        );
        let chord = a[0].mul(&b[1]).sub(&b[0].mul(&a[1]));
        acc.add(&area_term::<T>(&u.pcurve)).add(&chord)
    })
}

/// Number of +u ray crossings from p over the given loops; None when a
/// decision is not certified. Each loop is closed exactly by straight
/// segments from each use's end to the next use's start (gaps the uv_gap
/// check certified to be within tolerance), so parity is well defined.
fn crossings<T: Real>(loops: &[Vec<Coedge>], p: &V2<T>) -> Option<u32> {
    let mut count = 0;
    for lp in loops {
        for (k, u) in lp.iter().enumerate() {
            let next = &lp[(k + 1) % lp.len()];
            count += segment_crossing(
                pcurve_at::<T>(&u.pcurve, 1.0),
                pcurve_at::<T>(&next.pcurve, 0.0),
                p,
            )?;
            count += match &u.pcurve {
                Curve2::LineSegment { start, end } => {
                    segment_crossing([c(start.x), c(start.y)], [c(end.x), c(end.y)], p)?
                }
                Curve2::CircularArc {
                    center,
                    radius,
                    start_angle,
                    sweep_angle,
                } => arc_crossings(
                    [r(center.x), r(center.y)],
                    &r(*radius),
                    &r(*start_angle),
                    &r(*sweep_angle),
                    p,
                )?,
            };
        }
    }
    Some(count)
}

/// Half-open rule: a crossing when exactly one end lies strictly above p.v.
fn above<T: Real>(v: &T, pv: &T) -> Option<bool> {
    Some(v.sub(pv).sign()? == Ordering::Greater)
}

fn segment_crossing<T: Real>(a: V2<T>, b: V2<T>, p: &V2<T>) -> Option<u32> {
    let (ua, ub) = (above(&a[1], &p[1])?, above(&b[1], &p[1])?);
    if ua == ub {
        return Some(0);
    }
    // Certainly straddling p.v, so b.v - a.v excludes zero.
    let slope = b[0].sub(&a[0]).div(&b[1].sub(&a[1]))?;
    let x = a[0].add(&p[1].sub(&a[1]).mul(&slope));
    match x.sub(&p[0]).sign()? {
        Ordering::Greater => Some(1),
        Ordering::Less => Some(0),
        Ordering::Equal => None,
    }
}

/// Split at the v-extrema (angles pi/2 + k pi); each piece is v-monotone and
/// stays in one u half of the circle.
fn arc_crossings<T: Real>(center: [R; 2], rho: &R, start: &R, sweep: &R, p: &V2<T>) -> Option<u32> {
    let end = start + sweep;
    let (lo, hi) = if sweep.numer().sign() == Sign::Minus {
        (end.clone(), start.clone())
    } else {
        (start.clone(), end.clone())
    };
    let pi_i = pi();
    let approx = (pi_i.lo() + pi_i.hi()) / int(2);
    let half = &approx / int(2);
    let kmin: BigInt = ((&lo - &half) / &approx).floor().to_integer() - 1;
    let kmax: BigInt = ((&hi - &half) / &approx).ceil().to_integer() + 1;
    // (midpoint angle, extremum sign or None for an endpoint)
    let mut cuts: Vec<(R, Option<i64>)> = vec![(lo.clone(), None)];
    let mut k = kmin;
    while k <= kmax {
        let angle = pi_i
            .scale(&R::new(1.into(), 2.into()))
            .add(&pi_i.scale(&R::from_integer(k.clone())));
        let gt = angle.cmp_exact(&lo)?;
        let lt = angle.cmp_exact(&hi)?;
        if gt == Ordering::Greater && lt == Ordering::Less {
            // sin(pi/2 + k pi) = (-1)^k exactly.
            let s = if (&k % 2i32) == BigInt::from(0) {
                1
            } else {
                -1
            };
            cuts.push(((angle.lo() + angle.hi()) / int(2), Some(s)));
        }
        k += 1;
    }
    cuts.push((hi.clone(), None));
    let (cx, cy, rho_t) = (q::<T>(&center[0]), q::<T>(&center[1]), q::<T>(rho));
    let height = |cut: &(R, Option<i64>)| -> T {
        match cut.1 {
            Some(s) => cy.add(&rho_t.mul(&c(s as f64))),
            None => cy.add(&rho_t.mul(&T::cos_sin(&q(&cut.0)).1)),
        }
    };
    let mut count = 0;
    for pair in cuts.windows(2) {
        if above(&height(&pair[0]), &p[1])? == above(&height(&pair[1]), &p[1])? {
            continue;
        }
        let (cm, _) = T::cos_sin(&q(&((&pair[0].0 + &pair[1].0) / int(2))));
        let side = match cm.sign()? {
            Ordering::Greater => 1.0,
            Ordering::Less => -1.0,
            Ordering::Equal => return None,
        };
        let dy = p[1].sub(&cy);
        let reach = rho_t.square().sub(&dy.square()).sqrt();
        let x = cx.add(&reach.mul(&c(side)));
        match x.sub(&p[0]).sign()? {
            Ordering::Greater => count += 1,
            Ordering::Less => {}
            Ordering::Equal => return None,
        }
    }
    Some(count)
}

// ------------------------------------------------------------------ volumes

/// Integral over one use of G(u,v) dv, where dG/du = S.(S_u x S_v).
fn volume_term<T: Real>(s: &Surface, p: &Curve2) -> Option<T> {
    match s {
        Surface::Plane(f) => {
            let fr = frame::<T>(f);
            let h = vdot(&fr.o, &vcross(&fr.x, &fr.y));
            match p {
                Curve2::LineSegment { start, end } => {
                    let du = c::<T>(end.x).sub(&c(start.x));
                    let dv = c::<T>(end.y).sub(&c(start.y));
                    let mean = c::<T>(start.x).add(&du.mul(&c(0.5)));
                    Some(h.mul(&dv).mul(&mean))
                }
                Curve2::CircularArc {
                    center,
                    radius,
                    start_angle,
                    sweep_angle,
                } => {
                    let a0 = c::<T>(*start_angle);
                    let a1 = a0.add(&c(*sweep_angle));
                    let (_, s0) = T::cos_sin(&a0);
                    let (_, s1) = T::cos_sin(&a1);
                    let (_, t0) = T::cos_sin(&a0.mul(&c(2.0)));
                    let (_, t1) = T::cos_sin(&a1.mul(&c(2.0)));
                    let rho = c::<T>(*radius);
                    let first = c::<T>(center.x).mul(&rho).mul(&s1.sub(&s0));
                    let second = rho.square().mul(
                        &c::<T>(*sweep_angle)
                            .mul(&c(0.5))
                            .add(&t1.sub(&t0).mul(&c(0.25))),
                    );
                    Some(h.mul(&first.add(&second)))
                }
            }
        }
        Surface::Cylinder { frame: f, radius } => {
            let Curve2::LineSegment { start, end } = p else {
                return None;
            };
            let fr = frame::<T>(f);
            let rad = c::<T>(*radius);
            let ay = vdot(&fr.o, &vcross(&fr.y, &fr.n));
            let bx = vdot(&fr.o, &vcross(&fr.x, &fr.n));
            let det = vdot(&fr.x, &vcross(&fr.y, &fr.n));
            let (u0, u1) = (r(start.x), r(end.x));
            let du = &u1 - &u0;
            let dv = c::<T>(end.y).sub(&c(start.y));
            let (c0, s0) = T::cos_sin(&c(start.x));
            let (is, ic) = if du.numer().sign() == Sign::NoSign {
                (s0, c0)
            } else {
                let (c1, s1) = T::cos_sin(&c(end.x));
                let inv = q::<T>(&(int(1) / &du));
                (c0.sub(&c1).mul(&inv), s1.sub(&s0).mul(&inv))
            };
            let mean_u = q::<T>(&(&u0 + &du / int(2)));
            let g = ay
                .mul(&is)
                .add(&bx.mul(&ic))
                .add(&rad.mul(&det).mul(&mean_u));
            Some(rad.mul(&g).mul(&dv))
        }
    }
}

fn shell_volume<T: Real>(faces: &[Face], shell: &[FaceId]) -> Option<T> {
    let mut total = c::<T>(0.0);
    for f in shell {
        for u in faces[f.0].loops.iter().flatten() {
            total = total.add(&volume_term::<T>(&faces[f.0].surface, &u.pcurve)?);
        }
    }
    Some(total)
}

// ------------------------------------------------------------------ 3D containment

const DIRECTIONS: [[i64; 3]; 8] = [
    [3, 5, 7],
    [7, -2, 5],
    [-4, 9, 6],
    [5, 6, -8],
    [-9, 4, -7],
    [2, -11, 3],
    [11, 3, -2],
    [-3, -7, 10],
];

/// Parity of ray hits from an exact point against a shell; None when every
/// ray direction meets an uncertified decision. Adjacent faces only meet
/// within tolerance, so a hit also needs a certified margin `sqrt(margin2)`
/// from its face's boundary. The parity is then the same for any watertight
/// surface within tolerance of the faces.
fn inside<T: Real>(faces: &[Face], shell: &[FaceId], point: &[R; 3], margin2: &T) -> Option<bool> {
    'direction: for d in DIRECTIONS {
        let d: [R; 3] = d.map(int);
        let mut hits = 0;
        for face_id in shell {
            match face_hits::<T>(&faces[face_id.0], point, &d, margin2) {
                Some(n) => hits += n,
                None => continue 'direction,
            }
        }
        return Some(hits % 2 == 1);
    }
    None
}

fn det3(a: &[R; 3], b: &[R; 3], cc: &[R; 3]) -> R {
    &a[0] * (&b[1] * &cc[2] - &b[2] * &cc[1]) - &a[1] * (&b[0] * &cc[2] - &b[2] * &cc[0])
        + &a[2] * (&b[0] * &cc[1] - &b[1] * &cc[0])
}

fn face_hits<T: Real>(face: &Face, p: &[R; 3], d: &[R; 3], margin2: &T) -> Option<u32> {
    let exact = |f: &Frame3| {
        (
            f.origin().to_array().map(r),
            f.x().to_array().map(r),
            f.y().to_array().map(r),
            f.normal().to_array().map(r),
        )
    };
    match &face.surface {
        Surface::Plane(f) => {
            let (o, x, y, _) = exact(f);
            let qv: [R; 3] = std::array::from_fn(|i| &p[i] - &o[i]);
            let nd: [R; 3] = std::array::from_fn(|i| -&d[i]);
            let det = det3(&x, &y, &nd);
            if det.numer().sign() == Sign::NoSign {
                let n = [
                    &x[1] * &y[2] - &x[2] * &y[1],
                    &x[2] * &y[0] - &x[0] * &y[2],
                    &x[0] * &y[1] - &x[1] * &y[0],
                ];
                let side: R = (0..3).map(|i| &qv[i] * &n[i]).sum();
                return if side.numer().sign() == Sign::NoSign {
                    None
                } else {
                    Some(0)
                };
            }
            let u = det3(&qv, &y, &nd) / &det;
            let v = det3(&x, &qv, &nd) / &det;
            let t = det3(&x, &y, &qv) / &det;
            if t.numer().sign() == Sign::Minus {
                return Some(0);
            }
            let hit = [q(&u), q(&v)];
            if !clear_of_boundary(&face.loops, &hit, &c(1.0), margin2)? {
                return None;
            }
            let inside_face = crossings::<T>(&face.loops, &hit)? % 2;
            match t.numer().sign() {
                // The point lies on the face's plane: harmless only when it is
                // certainly outside the face region itself.
                Sign::NoSign => (inside_face == 0).then_some(0),
                _ => Some(inside_face),
            }
        }
        Surface::Cylinder { frame: f, radius } => {
            let (o, x, y, n) = exact(f);
            let det = det3(&x, &y, &n);
            let coords = |v: &[R; 3]| {
                [
                    det3(v, &y, &n) / &det,
                    det3(&x, v, &n) / &det,
                    det3(&x, &y, v) / &det,
                ]
            };
            let q0 = coords(&std::array::from_fn(|i| &p[i] - &o[i]));
            let q1 = coords(d);
            let rad = r(*radius);
            let a = &q1[0] * &q1[0] + &q1[1] * &q1[1];
            let b = int(2) * (&q0[0] * &q1[0] + &q0[1] * &q1[1]);
            let cc = &q0[0] * &q0[0] + &q0[1] * &q0[1] - &rad * &rad;
            if a.numer().sign() == Sign::NoSign {
                return if cc.numer().sign() == Sign::NoSign {
                    None
                } else {
                    Some(0)
                };
            }
            let disc = &b * &b - int(4) * &a * &cc;
            match disc.numer().sign() {
                Sign::Minus => return Some(0),
                Sign::NoSign => return None,
                Sign::Plus => {}
            }
            let root = q::<T>(&disc).sqrt();
            let rad_t = c::<T>(*radius);
            let window = uv_window(&face.loops)?;
            let (cw, sw) = T::cos_sin(&q(&window));
            let mut count = 0;
            for sign in [-1.0, 1.0] {
                let t = q::<T>(&-&b)
                    .add(&root.mul(&c(sign)))
                    .div(&q(&(int(2) * &a)))?;
                match t.sign()? {
                    Ordering::Less => continue,
                    Ordering::Equal => return None,
                    Ordering::Greater => {}
                }
                let alpha = q::<T>(&q0[0]).add(&t.mul(&q(&q1[0])));
                let beta = q::<T>(&q0[1]).add(&t.mul(&q(&q1[1])));
                let v = q::<T>(&q0[2]).add(&t.mul(&q(&q1[2])));
                // Angle relative to the face's UV window center.
                let ar = alpha.mul(&cw).add(&beta.mul(&sw));
                let br = beta.mul(&cw).sub(&alpha.mul(&sw));
                let hit = [q::<T>(&window).add(&T::atan2(&br, &ar)?), v];
                if !clear_of_boundary(&face.loops, &hit, &rad_t, margin2)? {
                    return None;
                }
                count += crossings::<T>(&face.loops, &hit)? % 2;
            }
            Some(count)
        }
    }
}

/// Whether a UV point is certainly farther than `sqrt(margin2)` from every
/// use and closing chord of the loops, with u distances scaled by `su` (the
/// cylinder radius; one on planes). Some(false) when it may be closer.
fn clear_of_boundary<T: Real>(
    loops: &[Vec<Coedge>],
    p: &V2<T>,
    su: &T,
    margin2: &T,
) -> Option<bool> {
    let scaled = |a: &V2<T>| [a[0].mul(su), a[1].clone()];
    let p = scaled(p);
    let far = |d2: &T| d2.cmp(margin2) == Some(Ordering::Greater);
    let dist2 = |a: &V2<T>| {
        let (du, dv) = (a[0].sub(&p[0]), a[1].sub(&p[1]));
        du.square().add(&dv.square())
    };
    // Squared distance from p to segment ab exceeds margin2.
    let segment_clear = |a: V2<T>, b: V2<T>| -> Option<bool> {
        let (a, b) = (scaled(&a), scaled(&b));
        let d = [b[0].sub(&a[0]), b[1].sub(&a[1])];
        let w = [p[0].sub(&a[0]), p[1].sub(&a[1])];
        let len2 = d[0].square().add(&d[1].square());
        if len2.sign()? != Ordering::Greater {
            return Some(far(&dist2(&a)));
        }
        let cross = d[0].mul(&w[1]).sub(&d[1].mul(&w[0]));
        if far(&cross.square().div(&len2)?) {
            return Some(true);
        }
        // Near the line: clear only beyond an end, away from both ends.
        let along = d[0].mul(&w[0]).add(&d[1].mul(&w[1]));
        let outside =
            along.sign()? == Ordering::Less || along.sub(&len2).sign()? == Ordering::Greater;
        Some(outside && far(&dist2(&a)) && far(&dist2(&b)))
    };
    for lp in loops {
        for (k, u) in lp.iter().enumerate() {
            let next = &lp[(k + 1) % lp.len()];
            if !segment_clear(pcurve_at(&u.pcurve, 1.0), pcurve_at(&next.pcurve, 0.0))? {
                return Some(false);
            }
            let clear = match &u.pcurve {
                Curve2::LineSegment { start, end } => {
                    segment_clear([c(start.x), c(start.y)], [c(end.x), c(end.y)])?
                }
                // Plane arcs only: clear of the whole circle.
                Curve2::CircularArc { center, radius, .. } => {
                    let r0 = dist2(&scaled(&[c(center.x), c(center.y)])).sqrt();
                    far(&r0.sub(&c(*radius)).square())
                }
            };
            if !clear {
                return Some(false);
            }
        }
    }
    Some(true)
}

/// Midpoint of the u extent of a face's line pcurve endpoints.
fn uv_window(loops: &[Vec<Coedge>]) -> Option<R> {
    let mut us: Vec<R> = Vec::new();
    for u in loops.iter().flatten() {
        match &u.pcurve {
            Curve2::LineSegment { start, end } => {
                us.push(r(start.x));
                us.push(r(end.x));
            }
            Curve2::CircularArc { .. } => return None,
        }
    }
    let lo = us.iter().min()?.clone();
    let hi = us.iter().max()?.clone();
    Some((lo + hi) / int(2))
}

// ------------------------------------------------------------------ the contract

pub(crate) fn check(
    vertices: &[Vertex],
    edges: &[Edge],
    faces: &[Face],
    shells: &[Vec<FaceId>],
    tolerance: Tolerance,
) -> Vec<Issue> {
    let mut issues: BTreeSet<Issue> = BTreeSet::new();
    let add = |issues: &mut BTreeSet<Issue>, kind, entity| {
        issues.insert(Issue { kind, entity });
    };
    use Entity as En;
    use IssueKind as K;
    let (nv, ne, nf) = (vertices.len(), edges.len(), faces.len());
    for (i, edge) in edges.iter().enumerate() {
        if edge.start.0 >= nv || edge.end.0 >= nv {
            add(&mut issues, K::Reference, En::Edge(i));
        }
    }
    for (fi, face) in faces.iter().enumerate() {
        for (li, lp) in face.loops.iter().enumerate() {
            for (ui, u) in lp.iter().enumerate() {
                if u.edge.0 >= ne {
                    add(&mut issues, K::Reference, En::Use(fi, li, ui));
                }
            }
        }
    }
    for (si, shell) in shells.iter().enumerate() {
        if shell.iter().any(|f| f.0 >= nf) {
            add(&mut issues, K::Reference, En::Shell(si));
        }
    }
    if !issues.is_empty() {
        return issues.into_iter().collect();
    }
    let tol = r(tolerance.linear());
    let tol2 = &tol * &tol;
    let (fast_tol, fast_tol2) = (Fast::from_r(&tol), Fast::from_r(&tol2));
    let (exact_tol, exact_tol2) = (I::from_r(&tol), I::from_r(&tol2));
    // Containment rays keep twice the tolerance from face boundaries.
    let (fast_margin2, exact_margin2) = (fast_tol2.mul(&c(4.0)), exact_tol2.mul(&c(4.0)));

    let mut used_vertex = vec![false; nv];
    for edge in edges {
        used_vertex[edge.start.0] = true;
        used_vertex[edge.end.0] = true;
    }
    for (v, used) in used_vertex.iter().enumerate() {
        if !used {
            add(&mut issues, K::UnusedVertex, En::Vertex(v));
        }
    }
    let mut uses: BTreeMap<usize, Vec<(usize, bool)>> = BTreeMap::new();
    for (fi, face) in faces.iter().enumerate() {
        for lp in &face.loops {
            for u in lp {
                uses.entry(u.edge.0)
                    .or_default()
                    .push((fi, u.orientation == Orientation::Forward));
            }
        }
    }
    for i in 0..ne {
        if !uses.contains_key(&i) {
            add(&mut issues, K::UnusedEdge, En::Edge(i));
        }
    }
    let mut owner: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (si, shell) in shells.iter().enumerate() {
        if shell.is_empty() {
            add(&mut issues, K::EmptyShell, En::Shell(si));
        }
        for f in shell {
            owner.entry(f.0).or_default().push(si);
        }
    }
    let mut structural_faces = BTreeSet::new();
    let mut bad_shells: BTreeSet<usize> = BTreeSet::new();
    for (fi, face) in faces.iter().enumerate() {
        match owner.get(&fi).map(Vec::len) {
            None => add(&mut issues, K::FaceWithoutShell, En::Face(fi)),
            Some(n) if n > 1 => add(&mut issues, K::FaceReused, En::Face(fi)),
            _ => {}
        }
        let mut broken = false;
        if face.loops.is_empty() {
            add(&mut issues, K::EmptyFace, En::Face(fi));
            broken = true;
        }
        for (li, lp) in face.loops.iter().enumerate() {
            if lp.is_empty() {
                add(&mut issues, K::EmptyLoop, En::Loop(fi, li));
                structural_faces.insert(fi);
                broken = true;
                continue;
            }
            let ends: Vec<_> = lp.iter().map(|u| use_vertices(edges, u)).collect();
            if (0..ends.len()).any(|k| ends[k].1 != ends[(k + 1) % ends.len()].0) {
                add(&mut issues, K::OpenLoop, En::Loop(fi, li));
                structural_faces.insert(fi);
                broken = true;
            }
        }
        if broken {
            bad_shells.extend(owner.get(&fi).into_iter().flatten().copied());
        }
    }

    // Edge accounting over faces owned by exactly one shell.
    let shell_of: BTreeMap<usize, usize> = owner
        .iter()
        .filter(|(_, s)| s.len() == 1)
        .map(|(f, s)| (*f, s[0]))
        .collect();
    for (edge, list) in &uses {
        let list: Vec<_> = list
            .iter()
            .filter(|x| shell_of.contains_key(&x.0))
            .collect();
        let set: BTreeSet<usize> = list.iter().map(|x| shell_of[&x.0]).collect();
        let kind = if set.len() > 1 {
            Some(K::EdgeAcrossShells)
        } else if list.len() == 1 {
            Some(K::FreeEdge)
        } else if list.len() > 2 {
            Some(K::NonManifoldEdge)
        } else if list.len() == 2 && list[0].1 == list[1].1 {
            Some(K::SameSenseUses)
        } else {
            None
        };
        if let Some(kind) = kind {
            add(&mut issues, kind, En::Edge(*edge));
            bad_shells.extend(set);
        }
    }
    for (si, shell) in shells.iter().enumerate() {
        if bad_shells.contains(&si)
            || shell.is_empty()
            || shell
                .iter()
                .any(|f| owner.get(&f.0).map_or(0, Vec::len) != 1)
        {
            bad_shells.insert(si);
            continue;
        }
        let members: BTreeSet<usize> = shell.iter().map(|f| f.0).collect();
        let mut adjacency: BTreeMap<usize, BTreeSet<usize>> =
            members.iter().map(|f| (*f, BTreeSet::new())).collect();
        for list in uses.values() {
            let fs: BTreeSet<usize> = list
                .iter()
                .map(|x| x.0)
                .filter(|f| members.contains(f))
                .collect();
            for a in &fs {
                adjacency
                    .get_mut(a)
                    .unwrap()
                    .extend(fs.iter().filter(|b| *b != a));
            }
        }
        if !connected(&adjacency) {
            add(&mut issues, K::DisconnectedShell, En::Shell(si));
            bad_shells.insert(si);
            continue;
        }
        let mut links: BTreeMap<usize, BTreeMap<usize, BTreeSet<usize>>> = BTreeMap::new();
        for f in &members {
            for lp in &faces[*f].loops {
                for (k, u) in lp.iter().enumerate() {
                    let w = &lp[(k + 1) % lp.len()];
                    let v = use_vertices(edges, u).1;
                    let g = links.entry(v).or_default();
                    g.entry(u.edge.0).or_default().insert(w.edge.0);
                    g.entry(w.edge.0).or_default().insert(u.edge.0);
                }
            }
        }
        let mut pinched = false;
        for (v, g) in &links {
            if !connected(g) {
                add(&mut issues, K::NonManifoldVertex, En::Vertex(*v));
                pinched = true;
            }
        }
        if pinched {
            bad_shells.insert(si);
            continue;
        }
        let mut vs = BTreeSet::new();
        let mut es = BTreeSet::new();
        let mut loops = 0i64;
        for f in &members {
            loops += faces[*f].loops.len() as i64;
            for u in faces[*f].loops.iter().flatten() {
                let (a, b) = use_vertices(edges, u);
                vs.insert(a);
                vs.insert(b);
                es.insert(u.edge.0);
            }
        }
        let chi = vs.len() as i64 - es.len() as i64 + 2 * members.len() as i64 - loops;
        if chi % 2 != 0 || chi > 2 {
            add(&mut issues, K::Euler, En::Shell(si));
            bad_shells.insert(si);
        }
    }

    // ------------------------------------------------ certified geometry
    let mut vertex_ok = vec![true; nv];
    for (v, vertex) in vertices.iter().enumerate() {
        if !finite(&vertex.position.to_array()) {
            vertex_ok[v] = false;
            add(&mut issues, K::DegenerateVertex, En::Vertex(v));
        }
    }
    let curve_ok: Vec<bool> = edges
        .iter()
        .map(|edge| curve_valid(&edge.curve, &tol, &fast_tol2, &exact_tol2))
        .collect();
    for (i, ok) in curve_ok.iter().enumerate() {
        if !ok {
            add(&mut issues, K::DegenerateCurve, En::Edge(i));
        }
    }
    let surface_ok: Vec<bool> = faces
        .iter()
        .map(|f| surface_valid(&f.surface, &tol))
        .collect();
    let mut bad_faces: BTreeSet<usize> = BTreeSet::new();
    for (fi, ok) in surface_ok.iter().enumerate() {
        if !ok {
            add(&mut issues, K::DegenerateSurface, En::Face(fi));
            bad_faces.insert(fi);
        }
    }
    for (i, edge) in edges.iter().enumerate() {
        if !curve_ok[i] {
            continue;
        }
        for (end, v, t) in [
            (EdgeEnd::Start, edge.start.0, 0.0),
            (EdgeEnd::End, edge.end.0, 1.0),
        ] {
            if !vertex_ok[v] {
                continue;
            }
            let gap = |curve: &Curve3| {
                tiered(
                    Verdict::Unknown,
                    || vertex_gap::<Fast>(curve, vertices[v].position.to_array(), t, &fast_tol2),
                    || vertex_gap::<I>(curve, vertices[v].position.to_array(), t, &exact_tol2),
                )
            };
            match gap(&edge.curve) {
                Verdict::Within => {}
                Verdict::Beyond => add(&mut issues, K::VertexOffCurve, En::EdgeEnd(i, end)),
                Verdict::Unknown => add(
                    &mut issues,
                    K::UncertifiedVertexOffCurve,
                    En::EdgeEnd(i, end),
                ),
            }
        }
    }
    for (fi, face) in faces.iter().enumerate() {
        for (li, lp) in face.loops.iter().enumerate() {
            for (ui, u) in lp.iter().enumerate() {
                if !pcurve_valid(&u.pcurve) {
                    add(&mut issues, K::DegeneratePcurve, En::Use(fi, li, ui));
                    bad_faces.insert(fi);
                    continue;
                }
                if !surface_ok[fi] || !curve_ok[u.edge.0] {
                    continue;
                }
                let forward = u.orientation == Orientation::Forward;
                let curve = &edges[u.edge.0].curve;
                match tiered(
                    Verdict::Unknown,
                    || {
                        deviation::<Fast>(
                            curve,
                            &face.surface,
                            &u.pcurve,
                            forward,
                            &fast_tol,
                            &fast_tol2,
                        )
                    },
                    || {
                        deviation::<I>(
                            curve,
                            &face.surface,
                            &u.pcurve,
                            forward,
                            &exact_tol,
                            &exact_tol2,
                        )
                    },
                ) {
                    Verdict::Within => {}
                    Verdict::Beyond => {
                        add(&mut issues, K::PcurveOffEdge, En::Use(fi, li, ui));
                        bad_faces.insert(fi);
                    }
                    Verdict::Unknown => {
                        add(
                            &mut issues,
                            K::UncertifiedPcurveOffEdge,
                            En::Use(fi, li, ui),
                        );
                        bad_faces.insert(fi);
                    }
                }
            }
        }
    }
    for (fi, face) in faces.iter().enumerate() {
        if !surface_ok[fi] {
            continue;
        }
        for (li, lp) in face.loops.iter().enumerate() {
            if lp.iter().any(|u| !pcurve_valid(&u.pcurve)) {
                continue;
            }
            for (ui, u) in lp.iter().enumerate() {
                let w = &lp[(ui + 1) % lp.len()];
                match tiered(
                    Verdict::Unknown,
                    || uv_gap::<Fast>(&face.surface, &u.pcurve, &w.pcurve, &fast_tol2),
                    || uv_gap::<I>(&face.surface, &u.pcurve, &w.pcurve, &exact_tol2),
                ) {
                    Verdict::Within => {}
                    Verdict::Beyond => {
                        add(&mut issues, K::UvGap, En::Use(fi, li, ui));
                        bad_faces.insert(fi);
                    }
                    Verdict::Unknown => {
                        add(&mut issues, K::UncertifiedUvGap, En::Use(fi, li, ui));
                        bad_faces.insert(fi);
                    }
                }
            }
        }
    }

    // Loop winding and imbrication on sound faces.
    let mut winding_faces = BTreeSet::new();
    for (fi, face) in faces.iter().enumerate() {
        if bad_faces.contains(&fi) || structural_faces.contains(&fi) || face.loops.is_empty() {
            continue;
        }
        for (li, lp) in face.loops.iter().enumerate() {
            let sign = tiered(
                None,
                || loop_area::<Fast>(lp).sign(),
                || loop_area::<I>(lp).sign(),
            );
            let outer = li == 0;
            let want = if outer == (face.orientation == Orientation::Forward) {
                Ordering::Greater
            } else {
                Ordering::Less
            };
            match sign {
                Some(s) if s == want => {}
                Some(_) => {
                    add(&mut issues, K::LoopWinding, En::Loop(fi, li));
                    winding_faces.insert(fi);
                }
                None => {
                    add(&mut issues, K::UncertifiedLoopWinding, En::Loop(fi, li));
                    winding_faces.insert(fi);
                }
            }
        }
        for li in 1..face.loops.len() {
            let outer = std::slice::from_ref(&face.loops[0]);
            let start = &face.loops[li][0].pcurve;
            match tiered(
                None,
                || crossings::<Fast>(outer, &pcurve_at(start, 0.0)),
                || crossings::<I>(outer, &pcurve_at(start, 0.0)),
            ) {
                Some(n) if n % 2 == 1 => {}
                Some(_) => add(&mut issues, K::InnerLoopOutside, En::Loop(fi, li)),
                None => add(&mut issues, K::UncertifiedContainment, En::Loop(fi, li)),
            }
        }
    }

    // Shell orientation and cavity nesting on fully sound shells.
    bad_faces.extend(structural_faces);
    bad_faces.extend(winding_faces);
    let sound: Vec<usize> = (0..shells.len())
        .filter(|si| {
            !bad_shells.contains(si) && shells[*si].iter().all(|f| !bad_faces.contains(&f.0))
        })
        .collect();
    let mut oriented = Vec::new();
    for &si in &sound {
        let sign = tiered(
            None,
            || shell_volume::<Fast>(faces, &shells[si])?.sign(),
            || shell_volume::<I>(faces, &shells[si])?.sign(),
        );
        let want = if si == 0 {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        match sign {
            Some(s) if s == want => oriented.push(si),
            Some(_) => add(&mut issues, K::ShellOrientation, En::Shell(si)),
            None => add(&mut issues, K::UncertifiedShellOrientation, En::Shell(si)),
        }
    }
    if oriented.first() == Some(&0) {
        let cavities: Vec<usize> = oriented.iter().copied().filter(|s| *s != 0).collect();
        for &si in &cavities {
            let first = &faces[shells[si][0].0].loops[0][0];
            let v = use_vertices(edges, first).0;
            if !vertex_ok[v] {
                continue;
            }
            let point = vertices[v].position.to_array().map(r);
            let contained = |shell: &[FaceId]| {
                tiered(
                    None,
                    || inside::<Fast>(faces, shell, &point, &fast_margin2),
                    || inside::<I>(faces, shell, &point, &exact_margin2),
                )
            };
            match contained(&shells[0]) {
                None => {
                    add(&mut issues, K::UncertifiedContainment, En::Shell(si));
                    continue;
                }
                Some(false) => {
                    add(&mut issues, K::CavityOutside, En::Shell(si));
                    continue;
                }
                Some(true) => {}
            }
            let mut nested = false;
            for &sj in &cavities {
                if sj == si {
                    continue;
                }
                match contained(&shells[sj]) {
                    Some(true) => nested = true,
                    Some(false) => {}
                    None => {
                        add(&mut issues, K::UncertifiedContainment, En::Shell(si));
                        nested = false;
                        break;
                    }
                }
            }
            if nested {
                add(&mut issues, K::NestedCavity, En::Shell(si));
            }
        }
    }
    issues.into_iter().collect()
}

fn connected<T: Ord + Copy>(graph: &BTreeMap<T, BTreeSet<T>>) -> bool {
    let Some(&start) = graph.keys().next() else {
        return true;
    };
    let mut seen = BTreeSet::new();
    let mut stack = vec![start];
    while let Some(x) = stack.pop() {
        if seen.insert(x) {
            stack.extend(graph[&x].iter().copied());
        }
    }
    seen.len() == graph.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::EdgeId;
    use crate::Point2;

    fn square() -> Vec<Vec<Coedge>> {
        let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        vec![(0..4)
            .map(|k| {
                let (a, b) = (corners[k], corners[(k + 1) % 4]);
                Coedge {
                    edge: EdgeId::new(k),
                    orientation: Orientation::Forward,
                    pcurve: Curve2::LineSegment {
                        start: Point2::new(a.0, a.1),
                        end: Point2::new(b.0, b.1),
                    },
                }
            })
            .collect()]
    }

    fn clear(u: f64, v: f64, su: f64, margin: f64) -> bool {
        let fast = clear_of_boundary::<Fast>(
            &square(),
            &[c(u), c(v)],
            &c(su),
            &c::<Fast>(margin).square(),
        );
        let exact =
            clear_of_boundary::<I>(&square(), &[c(u), c(v)], &c(su), &c::<I>(margin).square());
        assert_eq!(fast, exact, "{u} {v}");
        exact.unwrap()
    }

    #[test]
    fn containment_hits_keep_a_margin_from_every_boundary() {
        assert!(clear(0.5, 0.5, 1.0, 0.1));
        // Within the margin of an edge, including a hit on its shared line.
        assert!(!clear(0.5, 1e-9, 1.0, 1e-8));
        assert!(!clear(0.5, 0.0, 1.0, 1e-8));
        // Beyond a segment's end but near its line: judged by the end points.
        assert!(clear(2.0, 0.0, 1.0, 0.5));
        assert!(!clear(1.2, 0.0, 1.0, 0.5));
        // Cylinder u distances are scaled by the radius.
        assert!(clear(0.5, 0.5, 1.0, 0.4));
        assert!(!clear(0.5, 0.5, 0.1, 0.4));
    }
}
