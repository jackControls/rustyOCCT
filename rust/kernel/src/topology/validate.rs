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
use super::{
    Curve2, Curve3, Edge, Enclosure, Face, Fin, Loop, Orientation, Region, RegionKind, Shell,
    ShellId, Side, Surface, Vertex,
};
use crate::certified::{pi, Fast, Interval as I, Real};
use crate::jet::Jet;
use crate::{Frame3, Tolerance};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;
use std::fmt;

mod bernstein;
pub(crate) mod continuity;
mod mass;
mod projection;
mod quadrature;
pub(crate) use projection::{conic_point_fast, projection_range, section_rates};
mod spline_deviation;
mod spline_flux;
mod spline_taylor;
mod wall_meet;
pub(crate) use mass::{edge_length, face_mass, mass, sheet_measure};

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
    WireEdgeWithFins,
    AcornVertexUsed,
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
    // Cell-complex structure (TOPOLOGY_MODEL.md).
    FinWithoutLoop,
    FinReused,
    LoopWithoutFace,
    LoopReused,
    EdgeFinsMismatch,
    RingEdgeWithVertex,
    RingEdgeOpen,
    SideWithoutShell,
    SideInTwoShells,
    SideRegionMismatch,
    WindingMismatch,
    DoubleBounding,
    RegionWithoutShell,
    NoInfiniteRegion,
    RegionShellMismatch,
    SeamEdge,
    RadialOrderInconsistent,
    VertexLoopOffSurface,
    UncertifiedVertexLoop,
    // Enclosures (Contract 5 of IDENTITY_AND_HISTORY.md).
    EnclosureMissing,
    EnclosureExceedsResolution,
    EnclosureUnsound,
    UncertifiedEnclosure,
    // Poles (REVIEW_NOTES.md S3).
    PoleOffApex,
    // Continuity (R4 of REVIEW_NOTES.md).
    EdgeNotC1,
    PcurveNotC1,
    FaceNotC1,
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
            WireEdgeWithFins => "wire_edge_with_fins",
            AcornVertexUsed => "acorn_vertex_used",
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
            FinWithoutLoop => "fin_without_loop",
            FinReused => "fin_reused",
            LoopWithoutFace => "loop_without_face",
            LoopReused => "loop_reused",
            EdgeFinsMismatch => "edge_fins_mismatch",
            RingEdgeWithVertex => "ring_edge_with_vertex",
            RingEdgeOpen => "ring_edge_open",
            SideWithoutShell => "side_without_shell",
            SideInTwoShells => "side_in_two_shells",
            SideRegionMismatch => "side_region_mismatch",
            WindingMismatch => "winding_mismatch",
            DoubleBounding => "double_bounding",
            RegionWithoutShell => "region_without_shell",
            NoInfiniteRegion => "no_infinite_region",
            RegionShellMismatch => "region_shell_mismatch",
            SeamEdge => "seam_edge",
            RadialOrderInconsistent => "radial_order_inconsistent",
            VertexLoopOffSurface => "vertex_loop_off_surface",
            UncertifiedVertexLoop => "uncertified_vertex_loop",
            EnclosureMissing => "enclosure_missing",
            EnclosureExceedsResolution => "enclosure_exceeds_resolution",
            EnclosureUnsound => "enclosure_unsound",
            UncertifiedEnclosure => "uncertified_enclosure",
            PoleOffApex => "pole_off_apex",
            EdgeNotC1 => "edge_not_c1",
            PcurveNotC1 => "pcurve_not_c1",
            FaceNotC1 => "face_not_c1",
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
    Region(usize),
    /// A fin that no face's loop reaches, by arena index.
    FinSlot(usize),
    /// A loop that no face reaches, by arena index.
    LoopSlot(usize),
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
            Self::Region(r) => write!(f, "region {r}"),
            Self::FinSlot(k) => write!(f, "fin {k}"),
            Self::LoopSlot(l) => write!(f, "loop slot {l}"),
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

/// (frame, semi-axes along x and y, start, sweep) of an arc or an ellipse
/// arc; a full circle is the arc from 0 with sweep TAU.
fn arc_of(c: &Curve3) -> Option<(&Frame3, [f64; 2], f64, f64)> {
    match c {
        Curve3::Circle { frame, radius } => Some((frame, [*radius; 2], 0.0, TAU)),
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => Some((frame, [*radius; 2], *start_angle, *sweep_angle)),
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => Some((frame, [*major, *minor], *start_angle, *sweep_angle)),
        Curve3::LineSegment { .. }
        | Curve3::BSpline(_)
        | Curve3::HyperbolaArc { .. }
        | Curve3::ParabolaArc { .. }
        | Curve3::Section(_)
        | Curve3::Meet(_)
        | Curve3::Rise(_)
        | Curve3::Toric(_)
        | Curve3::WallMeet(_) => None,
    }
}

/// A spline's exact point at an edge or pcurve fraction, the parameter
/// `a + t (b - a)` exact (its ends exactly the domain's).
fn spline_point(curve: &crate::BSplineCurve3, range: [f64; 2], t: f64) -> [R; 3] {
    // The ends exactly; elsewhere the curve's exact point at the rounded
    // parameter (a point of the curve, at a fraction within rounding).
    let u = crate::topology::spline_parameter(range.into(), t);
    curve
        .exact_point(u)
        .expect("a spline evaluates exactly in its range")
}

/// `-∫ F du` along a chord from `a` to `b` when its closed form cannot
/// decide (a `du` straddling zero): `-du` times `F` over the chord's box,
/// which holds the chord's mean of `F`.
pub(super) fn chord_enclosure<T: Real>(
    a: &V2<T>,
    b: &V2<T>,
    f: &dyn Fn(&T, &T) -> Option<T>,
) -> Option<T> {
    let du = b[0].sub(&a[0]);
    Some(du.mul(&f(&a[0].union(&b[0]), &a[1].union(&b[1]))?).neg())
}

/// Halvings of a spline pcurve's Bézier pieces when a flux or mass
/// integral along it is enclosed (S4d).
pub(super) const SPLINE_DEPTH: usize = 6;

/// Whether every pole of a spline pcurve has the same u, so that `du` is
/// identically zero along it.
fn constant_u(spline: &super::SplineSpan<crate::BSplineCurve2>) -> bool {
    let poles = spline.curve().poles();
    poles.iter().all(|p| p.x == poles[0].x)
}

/// How a spline use's deviation is decided (S4b).
enum SplineUse {
    Exact(Box<spline_deviation::Exact>),
    Taylor(Box<spline_taylor::Taylor>),
}

/// Whether a use involves spline geometry, which only the endpoint and
/// continuity checks certify before the rest of S4.
fn spline_use(curve: &Curve3, s: &Surface, p: &Curve2) -> bool {
    // A spline wall's meeting's projection onto its own wall (S9f.2b) lies
    // on it by definition: no spline rule decides it.
    if let Curve2::Projection(pr) = p {
        if matches!(curve, Curve3::WallMeet(_)) && pr.curve == *curve && pr.surface == *s {
            return false;
        }
    }
    matches!(curve, Curve3::BSpline(_))
        || matches!(s, Surface::BSpline(_))
        || matches!(p, Curve2::BSpline(_))
}

fn curve_at<T: Real>(curve: &Curve3, t: f64) -> V3<T> {
    match curve {
        Curve3::BSpline(span) => {
            let t = if span.is_reversed() { 1.0 - t } else { t };
            spline_point(span.curve(), span.range(), t).map(|x| q::<T>(&x))
        }
        Curve3::LineSegment { start, .. } if t == 0.0 => v3::<T>(start.to_array()),
        Curve3::LineSegment { end, .. } if t == 1.0 => v3::<T>(end.to_array()),
        Curve3::LineSegment { start, end } => {
            let (a, b) = (v3::<T>(start.to_array()), v3::<T>(end.to_array()));
            vadd(&a, &vscale(&vsub(&b, &a), &c(t)))
        }
        Curve3::HyperbolaArc { .. }
        | Curve3::ParabolaArc { .. }
        | Curve3::Section(_)
        | Curve3::Meet(_)
        | Curve3::Rise(_)
        | Curve3::Toric(_)
        | Curve3::WallMeet(_) => {
            projection::conic_point::<T>(curve, t).expect("a conic or section evaluates")
        }
        _ => {
            let (f, [rx, ry], start, sweep) = arc_of(curve).unwrap();
            let fr = frame::<T>(f);
            let (co, si) = T::cos_sin(&c::<T>(start).add(&c::<T>(sweep).mul(&c(t))));
            vadd(
                &fr.o,
                &vadd(
                    &vscale(&fr.x, &c::<T>(rx).mul(&co)),
                    &vscale(&fr.y, &c::<T>(ry).mul(&si)),
                ),
            )
        }
    }
}

fn pcurve_at<T: Real>(p: &Curve2, t: f64) -> V2<T> {
    match p {
        Curve2::BSpline(spline) => {
            let t = if spline.is_reversed() { 1.0 - t } else { t };
            let [x, y, _] = spline_point(spline.curve().as_curve3(), spline.range(), t);
            [q(&x), q(&y)]
        }
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
        Curve2::EllipseArc {
            center,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => {
            let angle = c::<T>(*start_angle).add(&c::<T>(*sweep_angle).mul(&c(t)));
            let (co, si) = T::cos_sin(&angle);
            [
                c::<T>(center.x).add(&c::<T>(*major).mul(&co)),
                c::<T>(center.y).add(&c::<T>(*minor).mul(&si)),
            ]
        }
        Curve2::Projection(pr) => projection::projection_at::<T>(pr, t).unwrap_or_else(|| {
            let wide = c::<T>(0.0).widen(&int(1 << 30));
            [wide.clone(), wide]
        }),
        Curve2::Sinusoid { start, sweep, a } => {
            let u = c::<T>(*start).add(&c::<T>(*sweep).mul(&c(t)));
            let (co, si) = T::cos_sin(&u);
            let v = c::<T>(a[0])
                .add(&c::<T>(a[1]).mul(&co))
                .add(&c::<T>(a[2]).mul(&si));
            [u, v]
        }
    }
}

fn surface_at<T: Real>(s: &Surface, uv: &V2<T>) -> V3<T> {
    match s {
        Surface::BSpline(_) => unreachable!("spline uses are uncertified before sampling"),
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
        Surface::Cone {
            frame: f,
            radius,
            half_angle,
        } => {
            let fr = frame::<T>(f);
            let (ca, sa) = T::cos_sin(&c(*half_angle));
            let (co, si) = cos_sin_i(&uv[0]);
            let rho = c::<T>(*radius).add(&sa.mul(&uv[1]));
            let radial = vadd(&vscale(&fr.x, &rho.mul(&co)), &vscale(&fr.y, &rho.mul(&si)));
            vadd(&vadd(&fr.o, &radial), &vscale(&fr.n, &ca.mul(&uv[1])))
        }
        Surface::Sphere { frame: f, radius } => {
            let fr = frame::<T>(f);
            let (co, si) = cos_sin_i(&uv[0]);
            let (cv, sv) = T::cos_sin(&uv[1]);
            let rad = c::<T>(*radius);
            let rho = rad.mul(&cv);
            let radial = vadd(&vscale(&fr.x, &rho.mul(&co)), &vscale(&fr.y, &rho.mul(&si)));
            vadd(&vadd(&fr.o, &radial), &vscale(&fr.n, &rad.mul(&sv)))
        }
        Surface::Torus {
            frame: f,
            major,
            minor,
        } => {
            let fr = frame::<T>(f);
            let (co, si) = cos_sin_i(&uv[0]);
            let (cv, sv) = T::cos_sin(&uv[1]);
            let r = c::<T>(*minor);
            let rho = c::<T>(*major).add(&r.mul(&cv));
            let radial = vadd(&vscale(&fr.x, &rho.mul(&co)), &vscale(&fr.y, &rho.mul(&si)));
            vadd(&vadd(&fr.o, &radial), &vscale(&fr.n, &r.mul(&sv)))
        }
    }
}

/// A cone's apex parameter `v = -radius / sin a`; `None` for another surface
/// or when the division is not certain.
fn apex_v<T: Real>(s: &Surface) -> Option<T> {
    let Surface::Cone {
        radius, half_angle, ..
    } = s
    else {
        return None;
    };
    let (_, sa) = T::cos_sin(&c(*half_angle));
    c::<T>(*radius).neg().div(&sa)
}

/// A cone's apex.
fn apex_point<T: Real>(s: &Surface) -> Option<V3<T>> {
    let Surface::Cone {
        frame: f,
        half_angle,
        ..
    } = s
    else {
        return None;
    };
    let fr = frame::<T>(f);
    let (ca, _) = T::cos_sin(&c(*half_angle));
    Some(vadd(&fr.o, &vscale(&fr.n, &ca.mul(&apex_v::<T>(s)?))))
}

/// `pi / 2`, certified.
fn half_pi<T: Real>() -> T {
    T::from_r(&(pi().midpoint() / int(2))).widen(&(pi().radius() / int(2)))
}

/// The parameter `v` of a pole: a cone's apex, or a sphere's north (`pi/2`)
/// or south (`-pi/2`) pole; `None` for another surface.
fn pole_v<T: Real>(s: &Surface, north: bool) -> Option<T> {
    match s {
        Surface::Cone { .. } => apex_v(s),
        Surface::Sphere { .. } if north => Some(half_pi()),
        Surface::Sphere { .. } => Some(half_pi::<T>().neg()),
        _ => None,
    }
}

/// A pole's point (see [`pole_v`]).
fn pole_point<T: Real>(s: &Surface, north: bool) -> Option<V3<T>> {
    match s {
        Surface::Cone { .. } => apex_point(s),
        Surface::Sphere { frame: f, radius } => {
            let fr = frame::<T>(f);
            let r = c::<T>(*radius);
            Some(vadd(
                &fr.o,
                &vscale(&fr.n, &if north { r } else { r.neg() }),
            ))
        }
        _ => None,
    }
}

/// Which pole closes a band whose edge loops wind `turns` in total on a face
/// of `sense`: the material lies left of a loop's traversal, so a band
/// winding `+u` on a forward face closes at the larger `v`.
fn pole_north(turns: i32, sense: Orientation) -> bool {
    (turns > 0) == (sense == Orientation::Forward)
}

/// The edge loops' total winding in `u` of a face.
fn face_turns(face: &Face, loops: &[Loop]) -> i32 {
    face.loops
        .iter()
        .filter_map(|l| match loops.get(l.0) {
            Some(Loop::Edges { winding, .. }) => Some(winding[0]),
            _ => None,
        })
        .sum()
}

/// The length per unit of u at parameter v: a cylinder's radius, a cone's
/// `radius + v sin a`, a sphere's `radius cos v`, a torus's
/// `major + minor cos v`.
fn u_scale<T: Real>(s: &Surface, v: &T) -> Option<T> {
    match s {
        Surface::Plane(_) => None,
        Surface::Cylinder { radius, .. } => Some(c(*radius)),
        Surface::Cone {
            radius, half_angle, ..
        } => {
            let (_, sa) = T::cos_sin(&c(*half_angle));
            Some(c::<T>(*radius).add(&sa.mul(v)))
        }
        Surface::Sphere { radius, .. } => Some(c::<T>(*radius).mul(&T::cos_sin(v).0)),
        Surface::Torus { major, minor, .. } => {
            Some(c::<T>(*major).add(&c::<T>(*minor).mul(&T::cos_sin(v).0)))
        }
        // Not a length scale: a gap on a spline surface is decided only
        // when it is exactly zero (uv_gap).
        Surface::BSpline(_) => None,
    }
}

/// The length per unit of v: 1 along a cylinder's or cone's rulings, a
/// sphere's radius, a torus's minor radius.
fn v_scale(s: &Surface) -> f64 {
    match s {
        Surface::Sphere { radius, .. } => *radius,
        Surface::Torus { minor, .. } => *minor,
        _ => 1.0,
    }
}

/// The position in `face.loops` of a cone's or sphere's pole: the face's
/// first vertex loop, when its edge loops wind once in total, so the pole
/// closes the band (REVIEW_NOTES.md S3). On a sphere it is the pole on the
/// band's material side ([`pole_north`]).
pub(crate) fn pole_position(face: &Face, loops: &[Loop]) -> Option<usize> {
    if !matches!(face.surface, Surface::Cone { .. } | Surface::Sphere { .. }) {
        return None;
    }
    let total: i32 = face
        .loops
        .iter()
        .filter_map(|l| match loops.get(l.0) {
            Some(Loop::Edges { winding, .. }) => Some(winding[0]),
            _ => None,
        })
        .sum();
    if total.abs() != 1 {
        return None;
    }
    face.loops
        .iter()
        .position(|l| matches!(loops.get(l.0), Some(Loop::Vertex(_))))
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

/// A distance threshold (a finite binary64 value) and its square, enclosed
/// in either tier.
struct Threshold(f64);
impl Threshold {
    fn tier<T: Real>(&self) -> (T, T) {
        let t = T::exact_f64(self.0);
        let t2 = t.square();
        (t, t2)
    }
}

/// A decision against a stored enclosure, then against the resolution: the
/// geometric verdict, and the enclosure's verdict when the geometry is within
/// the resolution (`Beyond`: the bound is unsound; `Unknown`: uncertified).
fn bounded(
    bound: Option<f64>,
    tol: f64,
    decide: impl Fn(&Threshold) -> Verdict,
) -> (Verdict, Option<Verdict>) {
    let tol = Threshold(tol);
    let Some(b) = bound else {
        return (decide(&tol), None);
    };
    let at_bound = decide(&Threshold(b));
    if at_bound == Verdict::Within {
        return (Verdict::Within, None);
    }
    let at_tol = decide(&tol);
    let enclosure = (at_tol == Verdict::Within).then_some(at_bound);
    (at_tol, enclosure)
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

/// `tiered` for a sign from integrals along projections (S8d.2): binary64
/// intervals at the loose sign width first (`projection::LOOSE`), which
/// decide every sign not within it of zero, then at the tight one, then
/// exactly. Every tier's enclosure holds the value, so a sign any of them
/// decides is the value's.
fn tiered_integral<X: PartialEq>(unknown: X, fast: impl Fn() -> X, exact: impl FnOnce() -> X) -> X {
    let first = projection::with_sign_width(projection::LOOSE, &fast);
    if first != unknown {
        return first;
    }
    tiered(unknown, fast, exact)
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
        let ellipse = |cv: &V3<T>, sv: &V3<T>| {
            let (cc, ss, cs) = (vdot(cv, cv), vdot(sv, sv), vdot(cv, sv));
            let root = cc.sub(&ss).square().add(&cs.square().mul(&c(4.0))).sqrt();
            cc.add(&ss).add(&root).mul(&c(0.5)).sqrt()
        };
        // Two terms at nearby frequencies w1 < w2 are Re(z1 e^{i w1 t}) and
        // Re(z2 e^{i w2 t}) with z = c - i s; on [0,1] their sum is at most
        // the ellipse bound of z1 + z2 plus |z2| |w2 - w1|, since
        // |e^{i d t} - 1| <= |d| t. An arc and its pcurve whose sweeps differ
        // by a rounding cancel this way; the smaller of the two bounds is used.
        let mut terms: Vec<&(R, V3<T>, V3<T>)> = self.terms.iter().collect();
        terms.sort_by(|a, b| a.0.cmp(&b.0));
        let mut k = 0;
        while k < terms.len() {
            let (w1, c1, s1) = terms[k];
            let single = ellipse(c1, s1);
            if let Some((w2, c2, s2)) = terms.get(k + 1).copied() {
                let apart = single.add(&ellipse(c2, s2));
                let gap = T::from_r(&(w2 - w1));
                let near = vdot(c2, c2).add(&vdot(s2, s2)).sqrt().mul(&gap);
                let paired = ellipse(&vadd(c1, c2), &vadd(s1, s2)).add(&near);
                if paired.cmp(&apart) == Some(Ordering::Less) {
                    bound = bound.add(&paired);
                    k += 2;
                    continue;
                }
            }
            bound = bound.add(&single);
            k += 1;
        }
        bound
    }
}

/// Add C(t); false for a spline, which is not a harmonic sum.
fn add_curve<T: Real>(h: &mut Harmonic<T>, curve: &Curve3, forward: bool) -> bool {
    match curve {
        Curve3::BSpline(_)
        | Curve3::HyperbolaArc { .. }
        | Curve3::ParabolaArc { .. }
        | Curve3::Section(_)
        | Curve3::Meet(_)
        | Curve3::Rise(_)
        | Curve3::Toric(_)
        | Curve3::WallMeet(_) => return false,
        Curve3::LineSegment { start, end } => {
            let (a, b) = (v3::<T>(start.to_array()), v3::<T>(end.to_array()));
            if forward {
                h.affine(&a, &vsub(&b, &a), true);
            } else {
                h.affine(&b, &vsub(&a, &b), true);
            }
        }
        _ => {
            let (f, [rx, ry], start, sweep) = arc_of(curve).unwrap();
            let fr = frame::<T>(f);
            let (mut alpha, mut omega) = (c::<T>(start), r(sweep));
            if !forward {
                alpha = alpha.add(&c(sweep));
                omega = -omega;
            }
            h.affine(&fr.o, &zero3(), true);
            h.rotating(
                &alpha,
                &omega,
                &vscale(&fr.x, &c(rx)),
                &vscale(&fr.y, &c(ry)),
                true,
            );
        }
    }
    true
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
        (
            Surface::Plane(f),
            Curve2::EllipseArc {
                center,
                major,
                minor,
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
            h.rotating(
                &c(*start_angle),
                &r(*sweep_angle),
                &vscale(&fr.x, &c(*major)),
                &vscale(&fr.y, &c(*minor)),
                false,
            );
            true
        }
        // O + a0 n + cos u (R x + a1 n) + sin u (R y + a2 n).
        (Surface::Cylinder { frame: f, radius }, Curve2::Sinusoid { start, sweep, a }) => {
            let fr = frame::<T>(f);
            let rad = c::<T>(*radius);
            h.affine(&vadd(&fr.o, &vscale(&fr.n, &c(a[0]))), &zero3(), false);
            h.rotating(
                &c(*start),
                &r(*sweep),
                &vadd(&vscale(&fr.x, &rad), &vscale(&fr.n, &c(a[1]))),
                &vadd(&vscale(&fr.y, &rad), &vscale(&fr.n, &c(a[2]))),
                false,
            );
            true
        }
        (_, Curve2::EllipseArc { .. } | Curve2::Sinusoid { .. } | Curve2::Projection(_)) => false,
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
        // Harmonic along a ruling (du = 0) or a parallel (dv = 0) only.
        (
            Surface::Cone {
                frame: f,
                radius,
                half_angle,
            },
            Curve2::LineSegment { start, end },
        ) => {
            let fr = frame::<T>(f);
            let (ca, sa) = T::cos_sin(&c(*half_angle));
            let (du, dv) = (r(end.x) - r(start.x), r(end.y) - r(start.y));
            let zero = int(0);
            let v0 = c::<T>(start.y);
            let rho0 = c::<T>(*radius).add(&sa.mul(&v0));
            if du == zero {
                let (co, si) = T::cos_sin(&c(start.x));
                let e = vadd(&vscale(&fr.x, &co), &vscale(&fr.y, &si));
                let base = vadd(
                    &fr.o,
                    &vadd(&vscale(&e, &rho0), &vscale(&fr.n, &ca.mul(&v0))),
                );
                let dvt = c::<T>(end.y).sub(&v0);
                let slope = vadd(&vscale(&e, &sa.mul(&dvt)), &vscale(&fr.n, &ca.mul(&dvt)));
                h.affine(&base, &slope, false);
                true
            } else if dv == zero {
                h.affine(&vadd(&fr.o, &vscale(&fr.n, &ca.mul(&v0))), &zero3(), false);
                h.rotating(
                    &c(start.x),
                    &du,
                    &vscale(&fr.x, &rho0),
                    &vscale(&fr.y, &rho0),
                    false,
                );
                true
            } else {
                false
            }
        }
        (Surface::Cone { .. }, Curve2::CircularArc { .. }) => false,
        // A meridian (du = 0) is a great-circle arc in v; a parallel
        // (dv = 0) a circle in u.
        (Surface::Sphere { frame: f, radius }, Curve2::LineSegment { start, end }) => {
            let fr = frame::<T>(f);
            let rad = c::<T>(*radius);
            let (du, dv) = (r(end.x) - r(start.x), r(end.y) - r(start.y));
            let zero = int(0);
            if du == zero {
                let (co, si) = T::cos_sin(&c(start.x));
                let e = vadd(&vscale(&fr.x, &co), &vscale(&fr.y, &si));
                h.affine(&fr.o, &zero3(), false);
                h.rotating(
                    &c(start.y),
                    &dv,
                    &vscale(&e, &rad),
                    &vscale(&fr.n, &rad),
                    false,
                );
                true
            } else if dv == zero {
                let (cv, sv) = T::cos_sin(&c(start.y));
                let rho = rad.mul(&cv);
                h.affine(&vadd(&fr.o, &vscale(&fr.n, &rad.mul(&sv))), &zero3(), false);
                h.rotating(
                    &c(start.x),
                    &du,
                    &vscale(&fr.x, &rho),
                    &vscale(&fr.y, &rho),
                    false,
                );
                true
            } else {
                false
            }
        }
        (Surface::Sphere { .. }, Curve2::CircularArc { .. }) => false,
        // A meridian (du = 0) is an arc of the tube's circle about
        // O + major e(u0); a parallel (dv = 0) a circle about the axis.
        (
            Surface::Torus {
                frame: f,
                major,
                minor,
            },
            Curve2::LineSegment { start, end },
        ) => {
            let fr = frame::<T>(f);
            let (big, tube) = (c::<T>(*major), c::<T>(*minor));
            let (du, dv) = (r(end.x) - r(start.x), r(end.y) - r(start.y));
            let zero = int(0);
            if du == zero {
                let (co, si) = T::cos_sin(&c(start.x));
                let e = vadd(&vscale(&fr.x, &co), &vscale(&fr.y, &si));
                h.affine(&vadd(&fr.o, &vscale(&e, &big)), &zero3(), false);
                h.rotating(
                    &c(start.y),
                    &dv,
                    &vscale(&e, &tube),
                    &vscale(&fr.n, &tube),
                    false,
                );
                true
            } else if dv == zero {
                let (cv, sv) = T::cos_sin(&c(start.y));
                let rho = big.add(&tube.mul(&cv));
                h.affine(
                    &vadd(&fr.o, &vscale(&fr.n, &tube.mul(&sv))),
                    &zero3(),
                    false,
                );
                h.rotating(
                    &c(start.x),
                    &du,
                    &vscale(&fr.x, &rho),
                    &vscale(&fr.y, &rho),
                    false,
                );
                true
            } else {
                false
            }
        }
        (Surface::Torus { .. }, Curve2::CircularArc { .. }) => false,
        (Surface::BSpline(_), _) | (_, Curve2::BSpline(_)) => false,
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
    if spline_use(curve, s, p) {
        return Verdict::Unknown;
    }
    // A projection of the fin's own edge lies on it by definition (D13).
    if let Curve2::Projection(pr) = p {
        if projection::own(pr, curve, s, forward) {
            return Verdict::Within;
        }
    }
    let mut h = Harmonic::<T>::new();
    if add_curve(&mut h, curve, forward)
        && sub_use(&mut h, s, p)
        && matches!(h.upper().cmp(tol), Some(Ordering::Less | Ordering::Equal))
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
        // Degenerate when every pole lies within tolerance of the first:
        // the whole curve then lies within tolerance of a point.
        Curve3::BSpline(span) => {
            let poles = span.curve().poles();
            let first = poles[0].to_array().map(r);
            poles.iter().any(|p| {
                let d: [R; 3] = std::array::from_fn(|i| r(p.to_array()[i]) - &first[i]);
                d.iter().map(|x| x * x).sum::<R>() > tol * tol
            })
        }
        Curve3::HyperbolaArc {
            major,
            minor,
            start,
            sweep,
            ..
        } => {
            finite(&[*major, *minor, *start, *sweep])
                && r(*major) > *tol
                && r(*minor) > *tol
                && *sweep != 0.0
        }
        Curve3::ParabolaArc {
            focal,
            start,
            sweep,
            ..
        } => finite(&[*focal, *start, *sweep]) && *focal > 0.0 && *sweep != 0.0,
        // A torus section (S8d.3): a ring torus, a unit plane, a range within
        // a turn.
        Curve3::Section(s) => {
            finite(&[s.major, s.minor, s.start, s.sweep, s.sign])
                && finite(&s.plane)
                && r(s.minor) > *tol
                && s.major > s.minor
                && (s.sign == 1.0 || s.sign == -1.0)
                && s.sweep != 0.0
                && s.sweep.abs() <= TAU
                && (s.plane[0].hypot(s.plane[1]).hypot(s.plane[2]) - 1.0).abs() <= 1e-12
        }
        // Two cylinders' meeting (S9c.2): positive radii, axes not parallel,
        // a range within a turn.
        Curve3::Meet(m) => {
            let (n, x2, y2) = (m.frame.normal(), m.other.x(), m.other.y());
            let cone = m.half_angle != 0.0 || m.other_half_angle != 0.0;
            finite(&[
                m.radius,
                m.other_radius,
                m.start,
                m.sweep,
                m.sign,
                m.half_angle,
                m.other_half_angle,
            ]) && (r(m.radius) > *tol || (m.half_angle != 0.0 && m.radius >= 0.0))
                // A cone's radius at its frame's origin: zero at an apex.
                && (r(m.other_radius) > *tol
                    || (m.other_half_angle != 0.0 && m.other_radius >= 0.0))
                && m.half_angle.abs() < std::f64::consts::FRAC_PI_2
                && m.other_half_angle.abs() < std::f64::consts::FRAC_PI_2
                && (m.sign == 1.0 || m.sign == -1.0)
                && m.sweep != 0.0
                && m.sweep.abs() <= TAU
                && (m.other_sphere || cone || n.dot(x2).hypot(n.dot(y2)) > 1e-6)
        }
        // A cylinder's and a sphere's meeting over the height (S9d.2b):
        // positive radii, a sign, a nonempty height range, the sphere's
        // centre off the axis.
        Curve3::Rise(m) => {
            let ([a, b], _) = m.coefficients();
            finite(&[
                m.radius,
                m.sphere_radius,
                m.start,
                m.sweep,
                m.sign,
                m.half_angle,
            ]) && finite(&m.centre.to_array())
                && m.half_angle.abs() < std::f64::consts::FRAC_PI_2
                && (r(m.radius) > *tol || (m.half_angle != 0.0 && m.radius >= 0.0))
                && r(m.sphere_radius) > *tol
                && (m.sign == 1.0 || m.sign == -1.0)
                && m.sweep != 0.0
                && a.hypot(b) > 1e-12
        }
        // A spline wall's meeting with a cylinder (S9f.2b) or a sphere
        // (S9f.3a): a nonrational wall of degree one in `v` over two pole
        // rows, a range inside its `u` domain, a sign, a positive radius, a
        // ruling crossing a cylinder's axis, the discriminant positive at
        // the ends. A graph over `v` (S9f.2b.2): its window inside one knot
        // span of the `u` domain, its range inside the `v` domain, the other
        // surface's function along the ruling of opposite signs at the
        // window's ends at the range's ends and middle.
        Curve3::WallMeet(m) => {
            let s = &m.wall;
            let ((u0, u1), (v0, v1)) = s.domain();
            let (a, b) = (m.start, m.start + m.sweep);
            let inside = |x: f64| u0 <= x && x <= u1;
            // S9f.3a: a sphere meets every ruling's line; a cylinder only
            // one crossing its axis.
            let crossing = |u: f64| {
                let (_, dir) = m.ruling(u);
                let (x2, y2) = (m.other.x(), m.other.y());
                m.other_sphere || dir.dot(x2).hypot(dir.dot(y2)) > 1e-6 * dir.length()
            };
            let common = finite(&[m.other_radius, m.sign, m.start, m.sweep])
                && !s.is_rational()
                && !s.u_knots().is_periodic()
                && s.v_knots().degree() == 1
                && s.v_knots().pole_count() == 2
                && r(m.other_radius) > *tol
                && (m.sign == 1.0 || m.sign == -1.0)
                && m.sweep != 0.0;
            match m.window {
                None => common && inside(a) && inside(b) && crossing(a) && crossing(b),
                Some([wa, wb]) => {
                    let knots = s.u_knots().knots();
                    let one_span = !knots.iter().any(|k| wa < *k && *k < wb);
                    let in_v = |x: f64| v0 <= x && x <= v1;
                    let brackets = [0.0, 0.5, 1.0].iter().all(|f| {
                        let t = m.start + m.sweep * f - v0;
                        m.along(wa, t) * m.along(wb, t) < 0.0
                    });
                    common
                        && finite(&[wa, wb])
                        && wa < wb
                        && inside(wa)
                        && inside(wb)
                        && one_span
                        && in_v(a)
                        && in_v(b)
                        && crossing(0.5 * wa + 0.5 * wb)
                        && brackets
                }
            }
        }
        // A torus's meeting with a quadric (S9d.4b.2): a ring torus, a
        // positive radius (a cone's may be zero at its frame's origin), a
        // window of the other angle under a turn, a range within a turn;
        // with another ring torus (S9d.4b.2b), neither a sphere nor a cone.
        Curve3::Toric(m) => {
            finite(&[
                m.major,
                m.minor,
                m.other_radius,
                m.other_half_angle,
                m.other_minor,
                m.window[0],
                m.window[1],
                m.start,
                m.sweep,
            ]) && r(m.minor) > *tol
                && m.major > m.minor
                && (r(m.other_radius) > *tol
                    || (m.other_half_angle != 0.0 && m.other_radius >= 0.0))
                && m.other_half_angle.abs() < std::f64::consts::FRAC_PI_2
                && (!m.other_sphere || m.other_half_angle == 0.0)
                && (m.other_minor == 0.0
                    || (r(m.other_minor) > *tol
                        && m.other_radius > m.other_minor
                        && !m.other_sphere
                        && m.other_half_angle == 0.0))
                && m.window[0] < m.window[1]
                && m.window[1] - m.window[0] < TAU
                && m.sweep != 0.0
                && m.sweep.abs() <= TAU
        }
        _ => {
            let (_, [rx, ry], start, sweep) = arc_of(curve).unwrap();
            finite(&[rx, ry, start, sweep])
                && r(rx) > *tol
                && r(ry) > *tol
                && sweep != 0.0
                && sweep.abs() <= TAU
        }
    }
}

fn surface_valid(s: &Surface, tol: &R) -> bool {
    match s {
        Surface::Plane(_) => true,
        Surface::Cylinder { radius, .. } => radius.is_finite() && r(*radius) > *tol,
        // A cone may be given at its apex (radius 0); its angle is strictly
        // between 0 and a right angle in magnitude.
        Surface::Cone {
            radius, half_angle, ..
        } => {
            radius.is_finite()
                && *radius >= 0.0
                && half_angle.is_finite()
                && *half_angle != 0.0
                && half_angle.abs() < std::f64::consts::FRAC_PI_2
        }
        Surface::Sphere { radius, .. } => radius.is_finite() && r(*radius) > *tol,
        // A ring torus: the tube's radius above tolerance and clear of the
        // axis.
        Surface::Torus { major, minor, .. } => {
            major.is_finite()
                && minor.is_finite()
                && r(*minor) > *tol
                && r(*major) - r(*minor) > *tol
        }
        // Finite with positive weights by construction.
        Surface::BSpline(_) => true,
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
        // Degenerate when every pole is the first.
        Curve2::BSpline(span) => {
            let poles = span.curve().poles();
            poles.iter().any(|p| (p.x, p.y) != (poles[0].x, poles[0].y))
        }
        Curve2::EllipseArc {
            center,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => {
            finite(&[
                center.x,
                center.y,
                *major,
                *minor,
                *start_angle,
                *sweep_angle,
            ]) && *major > 0.0
                && *minor > 0.0
                && *sweep_angle != 0.0
                && sweep_angle.abs() <= TAU
        }
        Curve2::Projection(pr) => {
            pr.lifts.len() >= 2 && pr.lifts.iter().all(|l| finite(&[l.x, l.y]))
        }
        Curve2::Sinusoid { start, sweep, a } => {
            finite(&[*start, *sweep, a[0], a[1], a[2]]) && *sweep != 0.0 && sweep.abs() <= TAU
        }
    }
}

/// The arenas of one body, as validated.
pub(crate) struct View<'a> {
    pub vertices: &'a [Vertex],
    pub edges: &'a [Edge],
    pub fins: &'a [Fin],
    pub loops: &'a [Loop],
    pub faces: &'a [Face],
    pub shells: &'a [Shell],
    pub regions: &'a [Region],
}

/// A face's edge loop resolved: its fins in order and its windings in u
/// and v.
struct Lp<'a> {
    fins: Vec<&'a Fin>,
    winding: i32,
    winding_v: i32,
}

fn fin_vertices(edges: &[Edge], fin: &Fin) -> (Option<usize>, Option<usize>) {
    let edge = &edges[fin.edge.0];
    let (a, b) = (edge.start.map(|v| v.0), edge.end.map(|v| v.0));
    if fin.sense == Orientation::Forward {
        (a, b)
    } else {
        (b, a)
    }
}

/// Squared distance from a vertex to a curve's point at fraction `t`.
fn vertex_gap2<T: Real>(curve: &Curve3, vertex: [f64; 3], t: f64) -> T {
    let d = vsub(&curve_at::<T>(curve, t), &v3::<T>(vertex));
    vdot(&d, &d)
}

fn vertex_gap<T: Real>(curve: &Curve3, vertex: [f64; 3], t: f64, tol2: &T) -> Verdict {
    within(&vertex_gap2::<T>(curve, vertex, t), tol2)
}

/// The squared gap from the end of `p` to the start of `next` shifted by
/// `shift` in u, angles scaled by the radius.
/// A pcurve's end (`t` 0 or 1) exactly, for a line or spline.
fn pcurve_end_exact(p: &Curve2, t: f64) -> Option<[R; 2]> {
    match p {
        Curve2::LineSegment { start, end } => {
            let q = if t == 0.0 { start } else { end };
            Some([r(q.x), r(q.y)])
        }
        Curve2::BSpline(span) => {
            let t = if span.is_reversed() { 1.0 - t } else { t };
            let [x, y, _] = spline_point(span.curve().as_curve3(), span.range(), t);
            Some([x, y])
        }
        Curve2::CircularArc { .. }
        | Curve2::EllipseArc { .. }
        | Curve2::Sinusoid { .. }
        | Curve2::Projection(_) => None,
    }
}

/// The squared 3D distance between the surface points at two pcurves'
/// meeting ends on a spline surface, exactly (S4); `None` when an end is
/// not exact or off the surface's domain.
fn spline_gap2(s: &crate::BSplineSurface3, p: &Curve2, next: &Curve2) -> Option<R> {
    let (a, b) = (
        spline_point_at(s, pcurve_end_exact(p, 1.0)?)?,
        spline_point_at(s, pcurve_end_exact(next, 0.0)?)?,
    );
    Some((0..3).map(|k| (&a[k] - &b[k]) * (&a[k] - &b[k])).sum())
}

/// A pcurve's end (`t` 0 or 1) on a spline surface in 3D, enclosed: an
/// exact pcurve's on the surface's patch, a spline wall's meeting's
/// projection onto its own wall (S9f.2b) the meeting's own end point.
fn spline_end<T: Real>(s: &crate::BSplineSurface3, p: &Curve2, t: f64) -> Option<V3<T>> {
    if let Curve2::Projection(pr) = p {
        if matches!(&pr.curve, Curve3::WallMeet(m) if m.wall == *s) {
            let f = if pr.reversed { 1.0 - t } else { t };
            return projection::conic_point::<T>(&pr.curve, f);
        }
    }
    Some(spline_point_at(s, pcurve_end_exact(p, t)?)?.map(|x| q::<T>(&x)))
}

/// A spline surface's point at exact parameters, exactly (`None` off its
/// domain).
fn spline_point_at(s: &crate::BSplineSurface3, q: [R; 2]) -> Option<[R; 3]> {
    let patches = s.bezier_patches().ok()?;
    let ((ua, ub), (va, vb)) = s.domain();
    let (ua, ub, va, vb) = (r(ua), r(ub), r(va), r(vb));
    // A point within 1e-9 of the domain's width outside it (a pcurve end
    // printed at fifteen digits) is on the boundary patch's polynomial
    // extension, as OCCT evaluates it; farther out it is off the surface.
    let slack = |lo: &R, hi: &R| (hi - lo) * R::new(BigInt::from(1), BigInt::from(1_000_000_000));
    let (su, sv) = (slack(&ua, &ub), slack(&va, &vb));
    let at = |q: [R; 2]| -> Option<[R; 3]> {
        if q[0] < &ua - &su || q[0] > &ub + &su || q[1] < &va - &sv || q[1] > &vb + &sv {
            return None;
        }
        let clamp = |x: &R, lo: &R, hi: &R| x.clone().max(lo.clone()).min(hi.clone());
        let (cu, cv) = (clamp(&q[0], &ua, &ub), clamp(&q[1], &va, &vb));
        let patch = patches.iter().find(|b| {
            let [[u0, u1], [v0, v1]] = b.domain();
            *u0 <= cu && cu <= *u1 && *v0 <= cv && cv <= *v1
        })?;
        Some(patch_point(patch, &q[0], &q[1]))
    };
    at(q)
}

/// A Bézier patch's point at `(u, v)`, exactly, by de Casteljau: inside its
/// domain or on its polynomial extension.
fn patch_point(patch: &crate::ExactBezierSurface3, u: &R, v: &R) -> [R; 3] {
    let [du, dv] = patch.degrees();
    let [[u0, u1], [v0, v1]] = patch.domain().clone();
    let (s, t) = ((u - &u0) / (&u1 - &u0), (v - &v0) / (&v1 - &v0));
    let poles = patch.homogeneous_poles();
    let casteljau = |mut row: Vec<[R; 4]>, t: &R| -> [R; 4] {
        let one = R::from_integer(BigInt::from(1));
        while row.len() > 1 {
            row = row
                .windows(2)
                .map(|w| std::array::from_fn(|k| (&one - t) * &w[0][k] + t * &w[1][k]))
                .collect();
        }
        row.pop().unwrap()
    };
    let column: Vec<[R; 4]> = (0..=du)
        .map(|i| {
            casteljau(
                (0..=dv).map(|j| poles[i * (dv + 1) + j].clone()).collect(),
                &t,
            )
        })
        .collect();
    let h = casteljau(column, &s);
    [&h[0] / &h[3], &h[1] / &h[3], &h[2] / &h[3]]
}

fn uv_gap2<T: Real>(s: &Surface, p: &Curve2, next: &Curve2, shift: [f64; 2]) -> T {
    if let Surface::BSpline(spline) = s {
        // Measured in 3D, exactly; a spline wall's meeting's end enclosed
        // (S9f.2b); unknown otherwise.
        if let Some(g) = spline_gap2(spline, p, next) {
            return q(&g);
        }
        return match (
            spline_end::<T>(spline, p, 1.0),
            spline_end::<T>(spline, next, 0.0),
        ) {
            (Some(a), Some(b)) => {
                let d = vsub(&a, &b);
                vdot(&d, &d)
            }
            _ => c::<T>(0.0).widen(&R::from_integer(BigInt::from(1) << 1000)),
        };
    }
    let (a, b) = (pcurve_at::<T>(p, 1.0), pcurve_at::<T>(next, 0.0));
    let mut du = a[0].sub(&b[0].add(&c(shift[0])));
    if let Some(scale) = u_scale::<T>(s, &a[1]) {
        du = du.mul(&scale);
    }
    let dv = a[1].sub(&b[1].add(&c(shift[1]))).mul(&c(v_scale(s)));
    du.square().add(&dv.square())
}

fn uv_gap<T: Real>(s: &Surface, p: &Curve2, next: &Curve2, shift: [f64; 2], tol2: &T) -> Verdict {
    within(&uv_gap2::<T>(s, p, next, shift), tol2)
}

/// Squared distances from a point to the surface's pieces; the distance to
/// the surface is the smallest. One piece for a plane or cylinder; for a
/// cone the two generatrix lines of the meridian half-plane (both nappes):
/// `r cos a - R cos a - z sin a` and `r cos a + R cos a + z sin a`.
fn surface_gap2<T: Real>(s: &Surface, p: [f64; 3]) -> Vec<T> {
    // No certified distance to a spline surface yet.
    if matches!(s, Surface::BSpline(_)) {
        return Vec::new();
    }
    let (Surface::Plane(f)
    | Surface::Cylinder { frame: f, .. }
    | Surface::Cone { frame: f, .. }
    | Surface::Sphere { frame: f, .. }
    | Surface::Torus { frame: f, .. }) = s
    else {
        unreachable!("spline surfaces return above")
    };
    let fr = frame::<T>(f);
    let rel = vsub(&v3::<T>(p), &fr.o);
    let axial = vdot(&rel, &fr.n);
    let radial = || {
        let radial = vsub(&rel, &vscale(&fr.n, &axial));
        vdot(&radial, &radial).sqrt()
    };
    match s {
        Surface::Plane(_) => vec![axial.square()],
        Surface::Cylinder { radius, .. } => vec![radial().sub(&c(*radius)).square()],
        Surface::Cone {
            radius, half_angle, ..
        } => {
            let (ca, sa) = T::cos_sin(&c(*half_angle));
            let rc = radial().mul(&ca);
            let shift = c::<T>(*radius).mul(&ca).add(&axial.mul(&sa));
            vec![rc.sub(&shift).square(), rc.add(&shift).square()]
        }
        Surface::Sphere { radius, .. } => {
            let dist = vdot(&rel, &rel).sqrt();
            vec![dist.sub(&c(*radius)).square()]
        }
        // The distance to the tube's core circle, less the minor radius.
        Surface::Torus { major, minor, .. } => {
            let to_core = radial()
                .sub(&c(*major))
                .square()
                .add(&axial.square())
                .sqrt();
            vec![to_core.sub(&c(*minor)).square()]
        }
        Surface::BSpline(_) => unreachable!("returned above"),
    }
}

/// Distance from a point to a surface within tolerance: within when any
/// piece is, beyond when every piece is.
fn on_surface<T: Real>(s: &Surface, p: [f64; 3], tol2: &T) -> Verdict {
    let verdicts: Vec<Verdict> = surface_gap2::<T>(s, p)
        .iter()
        .map(|d2| within(d2, tol2))
        .collect();
    if verdicts.is_empty() {
        Verdict::Unknown
    } else if verdicts.contains(&Verdict::Within) {
        Verdict::Within
    } else if verdicts.iter().all(|v| *v == Verdict::Beyond) {
        Verdict::Beyond
    } else {
        Verdict::Unknown
    }
}

/// Squared distance from a point to a pole (see [`pole_v`]).
fn pole_gap2<T: Real>(s: &Surface, north: bool, p: [f64; 3]) -> Option<T> {
    let d = vsub(&v3::<T>(p), &pole_point::<T>(s, north)?);
    Some(vdot(&d, &d))
}

// ------------------------------------------------------------------ enclosures

/// The smallest stored bound, 2^-80: its square stays far above the
/// rational tier's 2^-192 grid, so an exactly zero gap still decides.
const MIN_BOUND: f64 = 8.271806125530277e-25;

/// A stored bound from a certified upper bound `x >= 0`: the next binary64
/// value above it, so the checker's comparison is strict and decides in the
/// same tier, and at least [`MIN_BOUND`].
fn next_above(x: f64) -> f64 {
    f64::from_bits(x.to_bits() + 1).max(MIN_BOUND)
}

/// A certified upper bound of `sqrt(value)`, from binary64 intervals and then
/// rational ones, or `None` when neither gives a finite bound.
fn root_bound(fast: impl FnOnce() -> Fast, exact: impl FnOnce() -> I) -> Option<f64> {
    let hi = fast().sqrt().bounds_f64().1;
    let hi = if hi.is_finite() {
        hi
    } else {
        exact().sqrt().bounds_f64().1
    };
    (hi.is_finite() && hi >= 0.0).then(|| next_above(hi))
}

/// Measured enclosures: certified upper bounds of every vertex, fin and face
/// gap of well-formed geometry, `None` where none is finite or the
/// deviation is not harmonic (an arc pcurve on a cylinder).
pub(crate) struct Measured {
    pub vertices: Vec<Option<f64>>,
    pub fins: Vec<Option<f64>>,
    pub faces: Vec<Option<f64>>,
}

pub(crate) fn measure(view: &View) -> Measured {
    bernstein::clear_memo();
    let max = |a: Option<f64>, b: Option<f64>| Some(a?.max(b?));
    // Every stored bound is at least MIN_BOUND; a face without gaps (a
    // whole sphere) gets exactly that.
    let mut vertices = vec![Some(MIN_BOUND); view.vertices.len()];
    for edge in view.edges {
        for (v, t) in [(edge.start, 0.0), (edge.end, 1.0)] {
            let Some(v) = v else { continue };
            let at = view.vertices[v.0].position.to_array();
            let bound = root_bound(
                || vertex_gap2::<Fast>(&edge.curve, at, t),
                || vertex_gap2::<I>(&edge.curve, at, t),
            );
            vertices[v.0] = max(vertices[v.0], bound);
        }
    }
    let mut fins = vec![None; view.fins.len()];
    let mut faces = vec![Some(MIN_BOUND); view.faces.len()];
    for (fi, face) in view.faces.iter().enumerate() {
        for l in &face.loops {
            match &view.loops[l.0] {
                Loop::Vertex(v) => {
                    let at = view.vertices[v.0].position.to_array();
                    // The nearest piece of the surface: the smallest bound.
                    let pieces = surface_gap2::<Fast>(&face.surface, at).len();
                    let bound = (0..pieces)
                        .map(|k| {
                            root_bound(
                                || surface_gap2::<Fast>(&face.surface, at).swap_remove(k),
                                || surface_gap2::<I>(&face.surface, at).swap_remove(k),
                            )
                        })
                        .fold(None, |a: Option<f64>, b| match (a, b) {
                            (Some(a), Some(b)) => Some(a.min(b)),
                            (a, b) => a.or(b),
                        });
                    vertices[v.0] = max(vertices[v.0], bound);
                    // A pole also stays at the apex.
                    let pole = pole_position(face, view.loops).map(|p| face.loops[p]);
                    if pole == Some(*l) {
                        let north = pole_north(face_turns(face, view.loops), face.sense);
                        let apex = root_bound(
                            || {
                                pole_gap2::<Fast>(&face.surface, north, at)
                                    .unwrap_or(c(f64::INFINITY))
                            },
                            || pole_gap2::<I>(&face.surface, north, at).unwrap_or(c(f64::MAX)),
                        );
                        vertices[v.0] = max(vertices[v.0], apex);
                    }
                }
                Loop::Edges {
                    fins: list,
                    winding,
                } => {
                    for (ui, k) in list.iter().enumerate() {
                        let u = &view.fins[k.0];
                        let curve = &view.edges[u.edge.0].curve;
                        let forward = u.sense == Orientation::Forward;
                        let harmonic = |h: &mut Harmonic<Fast>| {
                            add_curve(h, curve, forward) && sub_use(h, &face.surface, &u.pcurve)
                        };
                        let mut h = Harmonic::<Fast>::new();
                        // A fin that already has an enclosure (an imported
                        // tolerance) keeps it: nothing to measure.
                        fins[k.0] = if u.enclosure.is_some() {
                            None
                        } else if spline_use(curve, &face.surface, &u.pcurve) {
                            match spline_deviation::rational_use(
                                curve,
                                &face.surface,
                                &u.pcurve,
                                forward,
                            ) {
                                // Undecided exactly (a sliver's bound, a
                                // denominator): the Taylor enclosure.
                                Some(exact) => exact.upper_bound().or_else(|| {
                                    spline_taylor::taylor_use(
                                        curve,
                                        &face.surface,
                                        &u.pcurve,
                                        forward,
                                    )
                                    .and_then(|taylor| taylor.upper_bound())
                                }),
                                None => spline_taylor::taylor_use(
                                    curve,
                                    &face.surface,
                                    &u.pcurve,
                                    forward,
                                )
                                .and_then(|taylor| taylor.upper_bound()),
                            }
                            .map(next_above)
                        } else if matches!(&u.pcurve, Curve2::Projection(pr)
                            if projection::own(pr, curve, &face.surface, forward))
                        {
                            Some(MIN_BOUND)
                        } else if harmonic(&mut h) {
                            let hi = h.upper().bounds_f64().1;
                            let hi = if hi.is_finite() {
                                hi
                            } else {
                                let mut h = Harmonic::<I>::new();
                                add_curve(&mut h, curve, forward);
                                sub_use(&mut h, &face.surface, &u.pcurve);
                                h.upper().bounds_f64().1
                            };
                            (hi.is_finite() && hi >= 0.0).then(|| next_above(hi))
                        } else {
                            None
                        };
                        let w = &view.fins[list[(ui + 1) % list.len()].0];
                        let shift = if face.surface.is_periodic() && ui + 1 == list.len() {
                            [TAU * f64::from(winding[0]), TAU * f64::from(winding[1])]
                        } else {
                            [0.0, 0.0]
                        };
                        let gap = root_bound(
                            || uv_gap2::<Fast>(&face.surface, &u.pcurve, &w.pcurve, shift),
                            || uv_gap2::<I>(&face.surface, &u.pcurve, &w.pcurve, shift),
                        );
                        faces[fi] = max(faces[fi], gap);
                    }
                }
            }
        }
    }
    Measured {
        vertices,
        fins,
        faces,
    }
}

// ------------------------------------------------------------------ UV geometry

/// Twice the signed area contribution relative to the point `o` (exact):
/// integral of ((u - o_u) dv - (v - o_v) du).
fn area_term<T: Real>(p: &Curve2, o: &[R; 2]) -> T {
    let (ou, ov) = (T::from_r(&o[0]), T::from_r(&o[1]));
    match p {
        Curve2::LineSegment { start, end } => {
            let (su, sv) = (c::<T>(start.x).sub(&ou), c::<T>(start.y).sub(&ov));
            let (eu, ev) = (c::<T>(end.x).sub(&ou), c::<T>(end.y).sub(&ov));
            su.mul(&ev).sub(&eu.mul(&sv))
        }
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
                .sub(&ou)
                .mul(&s1.sub(&s0))
                .sub(&c::<T>(center.y).sub(&ov).mul(&c1.sub(&c0)));
            rho.mul(&linear).add(&rho.square().mul(&c(*sweep_angle)))
        }
        // Undecidable (a weight not certainly positive): any sign.
        Curve2::BSpline(spline) => bernstein::twice_area(spline, o)
            .unwrap_or_else(|| T::exact_f64(0.0).widen(&R::from_integer(BigInt::from(1) << 1000))),
        // The circle's terms with the radius along u and along v.
        Curve2::EllipseArc {
            center,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => {
            let a0 = c::<T>(*start_angle);
            let a1 = a0.add(&c(*sweep_angle));
            let ((c0, s0), (c1, s1)) = (T::cos_sin(&a0), T::cos_sin(&a1));
            let (rx, ry) = (c::<T>(*major), c::<T>(*minor));
            ry.mul(&c::<T>(center.x).sub(&ou).mul(&s1.sub(&s0)))
                .sub(&rx.mul(&c::<T>(center.y).sub(&ov).mul(&c1.sub(&c0))))
                .add(&rx.mul(&ry).mul(&c(*sweep_angle)))
        }
        // Along a projection (S8d.2), integrated with jets; unknown (any
        // sign) when the integral cannot be enclosed.
        Curve2::Projection(pr) => projection::integrate_along::<T>(pr, &|u, v, du, dv| {
            Some(
                u.add_constant(&ou.neg())
                    .mul(dv)
                    .sub(&v.add_constant(&ov.neg()).mul(du)),
            )
        })
        .unwrap_or_else(|| T::exact_f64(0.0).widen(&R::from_integer(BigInt::from(1) << 1000))),
        // With w = u - o_u: a1 w cos u + a2 w sin u - 2 a1 sin u + 2 a2 cos u
        // between the ends, less (a0 - o_v) times the sweep.
        Curve2::Sinusoid { start, sweep, a } => {
            let (a1, a2) = (c::<T>(a[1]), c::<T>(a[2]));
            let g = |u: &T| {
                let (co, si) = T::cos_sin(u);
                let w = u.sub(&ou);
                a1.mul(&w)
                    .mul(&co)
                    .add(&a2.mul(&w).mul(&si))
                    .sub(&a1.mul(&si).mul(&c(2.0)))
                    .add(&a2.mul(&co).mul(&c(2.0)))
            };
            let u0 = c::<T>(*start);
            let u1 = u0.add(&c(*sweep));
            g(&u1)
                .sub(&g(&u0))
                .sub(&c::<T>(a[0]).sub(&ov).mul(&c(*sweep)))
        }
    }
}

/// `-∫ v du` and `∫ u dv` along a sinusoid, in closed form.
fn sinusoid_areas<T: Real>(start: f64, sweep: f64, a: &[f64; 3]) -> (T, T) {
    let u0 = c::<T>(start);
    let u1 = u0.add(&c(sweep));
    let ((c0, s0), (c1, s1)) = (T::cos_sin(&u0), T::cos_sin(&u1));
    let (a0, a1, a2) = (c::<T>(a[0]), c::<T>(a[1]), c::<T>(a[2]));
    let minus_v_du = a0
        .mul(&c(sweep))
        .add(&a1.mul(&s1.sub(&s0)))
        .sub(&a2.mul(&c1.sub(&c0)))
        .neg();
    // ∫ u (-a1 sin u + a2 cos u) du = a1 (u cos u - sin u) + a2 (u sin u + cos u).
    let f = |u: &T, co: &T, si: &T| a1.mul(&u.mul(co).sub(si)).add(&a2.mul(&u.mul(si).add(co)));
    let u_dv = f(&u1, &c1, &s1).sub(&f(&u0, &c0, &s0));
    (minus_v_du, u_dv)
}

/// The chords closing a loop exactly: each fin's end to the next fin's
/// start, the last shifted by the winding times the period.
fn chords<T: Real>(lp: &Lp) -> Vec<(V2<T>, V2<T>)> {
    let n = lp.fins.len();
    (0..n)
        .map(|k| {
            let a = pcurve_at::<T>(&lp.fins[k].pcurve, 1.0);
            let mut b = pcurve_at::<T>(&lp.fins[(k + 1) % n].pcurve, 0.0);
            if k == n - 1 && lp.winding != 0 {
                b[0] = b[0].add(&c(TAU * f64::from(lp.winding)));
            }
            if k == n - 1 && lp.winding_v != 0 {
                b[1] = b[1].add(&c(TAU * f64::from(lp.winding_v)));
            }
            (a, b)
        })
        .collect()
}

/// The patch a closing chord's ends certainly share on a spline surface:
/// the first whose domain holds both enclosures, else the first holding
/// them within `2^-40` of its width past its sides inside the surface's
/// domain and anywhere past its sides on the domain's edges (S9f.2b: a
/// spline wall's meeting's end enclosed by its jets about a point on a
/// cap's or a knot's line; a chord of the loop's gap, its integrand on that
/// patch's polynomial or its extension, closes a gap within the tolerance;
/// past the domain's edge that patch's polynomial is the only one there).
fn chord_patch<'a, T: Real>(
    mut domains: impl Iterator<Item = &'a [[R; 2]; 2]> + Clone,
    a: &V2<T>,
    b: &V2<T>,
) -> Option<usize> {
    let mut whole: Option<[[R; 2]; 2]> = None;
    for d in domains.clone() {
        whole = Some(match whole {
            None => d.clone(),
            Some(w) => std::array::from_fn(|k| {
                [
                    w[k][0].clone().min(d[k][0].clone()),
                    w[k][1].clone().max(d[k][1].clone()),
                ]
            }),
        });
    }
    let whole = whole?;
    let holds = |d: &[[R; 2]; 2], x: &V2<T>, slack: bool| {
        d.iter()
            .zip(x)
            .zip(&whole)
            .all(|(([lo, hi], v), [first, last])| {
                let s = if slack {
                    (hi - lo) / R::from_integer(BigInt::from(1u64 << 40))
                } else {
                    R::from_integer(0.into())
                };
                (slack && lo == first
                    || matches!(
                        v.cmp(&q::<T>(&(lo - &s))),
                        Some(Ordering::Greater | Ordering::Equal)
                    ))
                    && (slack && hi == last
                        || matches!(
                            v.cmp(&q::<T>(&(hi + &s))),
                            Some(Ordering::Less | Ordering::Equal)
                        ))
            })
    };
    domains
        .clone()
        .position(|d| holds(d, a, false) && holds(d, b, false))
        .or_else(|| domains.position(|d| holds(d, a, true) && holds(d, b, true)))
}

/// Twice the signed area of an unwound loop closed by chords.
fn loop_area<T: Real>(lp: &Lp) -> T {
    // Relative to a point of the loop (any exact point gives the same
    // closed integral), so a loop far from the UV origin keeps its
    // enclosures at its own scale.
    let o: [R; 2] = match lp.fins.first() {
        Some(u) => pcurve_at::<T>(&u.pcurve, 0.0).map(|x| x.midpoint()),
        None => [R::from_integer(0.into()), R::from_integer(0.into())],
    };
    let (ou, ov) = (T::from_r(&o[0]), T::from_r(&o[1]));
    let mut total = c::<T>(0.0);
    for u in &lp.fins {
        total = total.add(&area_term::<T>(&u.pcurve, &o));
    }
    for (a, b) in chords::<T>(lp) {
        let (au, av) = (a[0].sub(&ou), a[1].sub(&ov));
        let (bu, bv) = (b[0].sub(&ou), b[1].sub(&ov));
        total = total.add(&au.mul(&bv).sub(&bu.mul(&av)));
    }
    total
}

/// -integral of v du over a loop on the universal cover, closed by chords;
/// seam segments (du = 0) would contribute nothing. Lines only.
fn periodic_area<T: Real>(lp: &Lp) -> Option<T> {
    let term = |a: &V2<T>, b: &V2<T>| a[1].add(&b[1]).mul(&c(-0.5)).mul(&b[0].sub(&a[0]));
    let mut total = c::<T>(0.0);
    for u in &lp.fins {
        total = total.add(&match &u.pcurve {
            Curve2::LineSegment { start, end } => {
                term(&[c(start.x), c(start.y)], &[c(end.x), c(end.y)])
            }
            Curve2::BSpline(spline) => bernstein::minus_v_du(spline)?,
            Curve2::Sinusoid { start, sweep, a } => sinusoid_areas::<T>(*start, *sweep, a).0,
            Curve2::Projection(pr) => {
                projection::integrate_along::<T>(pr, &|_, v, du, _| Some(v.mul(du).neg()))?
            }
            Curve2::CircularArc { .. } | Curve2::EllipseArc { .. } => return None,
        });
    }
    for (a, b) in chords::<T>(lp) {
        total = total.add(&term(&a, &b));
    }
    Some(total)
}

/// The integral of u dv over a loop on the universal cover, closed by chords:
/// the periodic area of a loop winding in v (on a torus). Lines only.
fn periodic_area_v<T: Real>(lp: &Lp) -> Option<T> {
    let term = |a: &V2<T>, b: &V2<T>| a[0].add(&b[0]).mul(&c(0.5)).mul(&b[1].sub(&a[1]));
    let mut total = c::<T>(0.0);
    for u in &lp.fins {
        total = total.add(&match &u.pcurve {
            Curve2::LineSegment { start, end } => {
                term(&[c(start.x), c(start.y)], &[c(end.x), c(end.y)])
            }
            Curve2::BSpline(spline) => bernstein::u_dv(spline)?,
            Curve2::Sinusoid { start, sweep, a } => sinusoid_areas::<T>(*start, *sweep, a).1,
            Curve2::Projection(pr) => {
                projection::integrate_along::<T>(pr, &|u, _, _, dv| Some(u.mul(dv)))?
            }
            Curve2::CircularArc { .. } | Curve2::EllipseArc { .. } => return None,
        });
    }
    for (a, b) in chords::<T>(lp) {
        total = total.add(&term(&a, &b));
    }
    Some(total)
}

/// Number of +u ray crossings from p over the given unwound loops (closed by
/// chords); None when a decision is not certified.
fn crossings<T: Real>(loops: &[&Lp], p: &V2<T>) -> Option<u32> {
    let mut count = 0;
    for lp in loops {
        for (a, b) in chords::<T>(lp) {
            count += segment_crossing(a, b, p)?;
        }
        for u in &lp.fins {
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
                Curve2::BSpline(spline) => bernstein::crossing_parity(spline, p)?,
                // An axis-aligned ellipse's crossings are its unit circle's
                // after scaling each axis about the centre.
                Curve2::EllipseArc {
                    center,
                    major,
                    minor,
                    start_angle,
                    sweep_angle,
                } => {
                    let scaled = [
                        p[0].sub(&c(center.x)).div(&c(*major))?,
                        p[1].sub(&c(center.y)).div(&c(*minor))?,
                    ];
                    let zero = R::from_integer(0.into());
                    arc_crossings(
                        [zero.clone(), zero],
                        &int(1),
                        &r(*start_angle),
                        &r(*sweep_angle),
                        &scaled,
                    )?
                }
                Curve2::Projection(pr) => projection::crossings::<T>(pr, p)?,
                Curve2::Sinusoid { start, sweep, a } => sinusoid_crossings(*start, *sweep, a, p)?,
            };
        }
    }
    Some(count)
}

/// Crossings of the +u ray from p with `v = a0 + a1 cos u + a2 sin u` over
/// `u` in the pcurve's range: the roots of `rho cos(u - phi) = v - a0`,
/// `phi +- acos(w) + 2 pi k`, each certainly inside the range and past p.u
/// (a root at an end or at p.u is undecided).
fn sinusoid_crossings<T: Real>(start: f64, sweep: f64, a: &[f64; 3], p: &V2<T>) -> Option<u32> {
    let end = start + sweep;
    let (lo, hi) = if sweep < 0.0 {
        (end, start)
    } else {
        (start, end)
    };
    let (a0, a1, a2) = (c::<T>(a[0]), c::<T>(a[1]), c::<T>(a[2]));
    if a[1] == 0.0 && a[2] == 0.0 {
        // A line at height a0.
        return segment_crossing([c(lo), a0.clone()], [c(hi), a0], p);
    }
    let rho = a1.square().add(&a2.square()).sqrt();
    let w = p[1].sub(&a0).div(&rho)?;
    let one = c::<T>(1.0);
    match (w.cmp(&one)?, w.cmp(&one.neg())?) {
        (Ordering::Greater, _) | (_, Ordering::Less) => return Some(0),
        (Ordering::Less, Ordering::Greater) => {}
        _ => return None,
    }
    let phi = T::atan2(&a2, &a1)?;
    let theta = T::atan2(&one.sub(&w.square()).sqrt(), &w)?;
    let pi_i = pi();
    let two_pi = T::from_r(&(pi_i.midpoint() * int(2))).widen(&(pi_i.radius() * int(2)));
    let (plo, _) = phi.bounds_f64();
    let kmin = ((lo - plo - 4.0) / TAU).floor() as i64 - 1;
    let kmax = ((hi - plo + 4.0) / TAU).ceil() as i64 + 1;
    let (lo_t, hi_t) = (c::<T>(lo), c::<T>(hi));
    let mut count = 0;
    for k in kmin..=kmax {
        let base = phi.add(&two_pi.mul(&c(k as f64)));
        for root in [base.add(&theta), base.sub(&theta)] {
            let inside = match (root.cmp(&lo_t)?, root.cmp(&hi_t)?) {
                (Ordering::Greater, Ordering::Less) => true,
                (Ordering::Less, _) | (_, Ordering::Greater) => false,
                _ => return None,
            };
            if inside && root.cmp(&p[0])? == Ordering::Greater {
                count += 1;
            }
        }
    }
    Some(count)
}

/// Crossings of the +v ray from p with the face's pcurves and chords, over
/// every u alias k*TAU (a cylinder's universal cover); None when a decision
/// is not certified or a pcurve is not a line, a sinusoid or a projection.
fn cover_crossings<T: Real>(loops: &[&Lp], p: &V2<T>) -> Option<u32> {
    Some(cover_hits(loops, p)?.len() as u32)
}

/// Signed crossings of the +v ray from p with the loops' pcurves and
/// chords, over every u alias: +1 where a piece runs in -u, -1 in +u (the
/// region on the left lies below a piece running in -u). None when a
/// decision is not certified or a pcurve is not a line, a sinusoid or a
/// projection.
fn signed_cover_crossings<T: Real>(loops: &[&Lp], p: &V2<T>) -> Option<i64> {
    Some(cover_hits(loops, p)?.iter().sum())
}

/// The +u ray's crossings from `p` (each +1 running in -v, -1 in +v) with
/// the loops' chords and pcurves over every `v` alias: `cover_hits` with
/// the parameters' roles exchanged (S9e.3b: a hole in a torus face wound
/// in `v`); a sinusoid's pcurve (a cylinder's) is not taken.
fn cover_hits_v<T: Real>(loops: &[&Lp], p: &V2<T>) -> Option<Vec<i64>> {
    let swap = |x: &V2<T>| [x[1].clone(), x[0].clone()];
    let q = swap(p);
    let mut out = Vec::new();
    for lp in loops {
        let mut segments: Vec<(V2<T>, V2<T>)> = chords::<T>(lp)
            .into_iter()
            .map(|(a, b)| (swap(&a), swap(&b)))
            .collect();
        for u in &lp.fins {
            match &u.pcurve {
                Curve2::LineSegment { start, end } => {
                    segments.push(([c(start.y), c(start.x)], [c(end.y), c(end.x)]))
                }
                Curve2::Projection(pr) => {
                    out.extend(projection::cover_crossings_on::<T>(pr, &q, true)?)
                }
                _ => return None,
            }
        }
        for (a, b) in segments {
            for k in alias_range(&a[0], &b[0], &q[0])? {
                let u = q[0].sub(&c(TAU * k as f64));
                let (ra, rb) = (above(&a[0], &u)?, above(&b[0], &u)?);
                if ra == rb {
                    continue;
                }
                let slope = b[1].sub(&a[1]).div(&b[0].sub(&a[0]))?;
                let v = a[1].add(&u.sub(&a[0]).mul(&slope));
                match v.sub(&q[1]).sign()? {
                    Ordering::Greater => out.push(if rb { -1 } else { 1 }),
                    Ordering::Less => {}
                    Ordering::Equal => return None,
                }
            }
        }
    }
    Some(out)
}

/// The +v ray's crossings from p (each +1 running in -u, -1 in +u) with
/// the loops' chords and pcurves over every u alias: lines exactly as
/// segments, sinusoids at the alias (S9c.1), projections by certified
/// pieces (S9c.2).
fn cover_hits<T: Real>(loops: &[&Lp], p: &V2<T>) -> Option<Vec<i64>> {
    let mut out = Vec::new();
    for lp in loops {
        let mut segments = chords::<T>(lp);
        for u in &lp.fins {
            match &u.pcurve {
                Curve2::LineSegment { start, end } => {
                    segments.push(([c(start.x), c(start.y)], [c(end.x), c(end.y)]))
                }
                Curve2::Sinusoid { start, sweep, a } => {
                    out.extend(sinusoid_cover(*start, *sweep, a, p)?)
                }
                Curve2::Projection(pr) => out.extend(projection::cover_crossings::<T>(pr, p)?),
                _ => return None,
            }
        }
        for (a, b) in segments {
            for k in alias_range(&a[0], &b[0], &p[0])? {
                let u = p[0].sub(&c(TAU * k as f64));
                // Half-open in u: exactly one end strictly right of u.
                let (ra, rb) = (above(&a[0], &u)?, above(&b[0], &u)?);
                if ra == rb {
                    continue;
                }
                let slope = b[1].sub(&a[1]).div(&b[0].sub(&a[0]))?;
                let v = a[1].add(&u.sub(&a[0]).mul(&slope));
                match v.sub(&p[1]).sign()? {
                    Ordering::Greater => out.push(if rb { -1 } else { 1 }),
                    Ordering::Less => {}
                    Ordering::Equal => return None,
                }
            }
        }
    }
    Some(out)
}

/// The +v ray's crossings from p with `v = a0 + a1 cos u + a2 sin u` over
/// its range, at every u alias (half-open in u, as a segment's).
fn sinusoid_cover<T: Real>(start: f64, sweep: f64, a: &[f64; 3], p: &V2<T>) -> Option<Vec<i64>> {
    let (s, e) = (c::<T>(start), c::<T>(start + sweep));
    let mut out = Vec::new();
    for k in alias_range(&s, &e, &p[0])? {
        let u = p[0].sub(&c(TAU * k as f64));
        let (ra, rb) = (above(&s, &u)?, above(&e, &u)?);
        if ra == rb {
            continue;
        }
        let (co, si) = T::cos_sin(&u);
        let v = c::<T>(a[0])
            .add(&c::<T>(a[1]).mul(&co))
            .add(&c::<T>(a[2]).mul(&si));
        match v.sub(&p[1]).sign()? {
            Ordering::Greater => out.push(if rb { -1 } else { 1 }),
            Ordering::Less => {}
            Ordering::Equal => return None,
        }
    }
    Some(out)
}

/// Every k with p.u - k TAU possibly within the u span of a..b, plus one each
/// side; None for spans too wide to enumerate.
fn alias_range<T: Real>(a: &T, b: &T, pu: &T) -> Option<std::ops::RangeInclusive<i64>> {
    let ((a0, a1), (b0, b1), (p0, p1)) = (a.bounds_f64(), b.bounds_f64(), pu.bounds_f64());
    let (lo, hi) = (a0.min(b0), a1.max(b1));
    let kmin = ((p0 - hi) / TAU).floor() - 1.0;
    let kmax = ((p1 - lo) / TAU).ceil() + 1.0;
    (kmin.is_finite() && kmax.is_finite() && kmax - kmin <= 64.0)
        .then_some(kmin as i64..=kmax as i64)
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

/// -integral of v f(u) du along a line from a to b, where f = A sin u +
/// B cos u + C (A = B = 0 on a plane); exact closed forms, None when a
/// cylinder line's du may be zero without being zero.
fn line_flux<T: Real>(a: &V2<T>, b: &V2<T>, coeffs: &(T, T, T), plane: bool) -> Option<T> {
    let (ca, cb, cc) = coeffs;
    let du = b[0].sub(&a[0]);
    if plane {
        return Some(cc.mul(&du).mul(&a[1].add(&b[1])).mul(&c(-0.5)));
    }
    if du.sign()? == Ordering::Equal {
        return Some(c(0.0));
    }
    let (c0, s0) = T::cos_sin(&a[0]);
    let (c1, s1) = T::cos_sin(&b[0]);
    // I0 = integral of f, I1 = integral of s f, over s in [0, du].
    let i0 = ca
        .mul(&c0.sub(&c1))
        .add(&cb.mul(&s1.sub(&s0)))
        .add(&cc.mul(&du));
    let i1 = ca
        .mul(&du.mul(&c1).neg().add(&s1).sub(&s0))
        .add(&cb.mul(&du.mul(&s1).add(&c1).sub(&c0)))
        .add(&cc.mul(&du.square()).mul(&c(0.5)));
    let dv = b[1].sub(&a[1]);
    Some(a[1].mul(&i0).add(&dv.mul(&i1).div(&du)?).neg())
}

/// On a cone, `S.(S_u x S_v) = rho(v) h(u)` with `rho = R + v sin a` and
/// `h = A cos u + B sin u + C`. Its v-antiderivative is taken from the apex
/// on a face with a pole, `rho^2 / (2 sin a) h(u)`, zero at the pole, and
/// from 0 otherwise, `(R v + v^2 sin a / 2) h(u)`, which needs no division
/// by `sin a` (a nearly cylindrical cone has its apex far away; any constant
/// cancels over a band's two loops). This is `-integral of that du` along a
/// line from a to b, in closed form with the moments `integral of w^k
/// cos(u0 + w)` and `w^k sin(u0 + w)` for k <= 2; `None` when du may be zero
/// without being zero.
fn cone_line_flux<T: Real>(
    a: &V2<T>,
    b: &V2<T>,
    sa: &T,
    radius: &T,
    h: &(T, T, T),
    poled: bool,
) -> Option<T> {
    let (ha, hb, hc) = h;
    let d = b[0].sub(&a[0]);
    if d.sign()? == Ordering::Equal {
        return Some(c(0.0));
    }
    let m = b[1].sub(&a[1]).div(&d)?;
    let rho0 = radius.add(&sa.mul(&a[1]));
    let (c0, s0) = T::cos_sin(&a[0]);
    let (c1, s1) = T::cos_sin(&b[0]);
    // Moments over w in [0, d] of cos(u0 + w) and sin(u0 + w).
    let jc0 = s1.sub(&s0);
    let js0 = c0.sub(&c1);
    let jc1 = d.mul(&s1).add(&c1).sub(&c0);
    let js1 = d.mul(&c1).neg().add(&s1).sub(&s0);
    let jc2 = d.square().mul(&s1).sub(&js1.mul(&c(2.0)));
    let js2 = d.square().mul(&c1).neg().add(&jc1.mul(&c(2.0)));
    let k = |jc: &T, js: &T, power: T| ha.mul(jc).add(&hb.mul(js)).add(&hc.mul(&power));
    let k0 = k(&jc0, &js0, d.clone());
    let k1 = k(&jc1, &js1, d.square().mul(&c(0.5)));
    let k2 = k(&jc2, &js2, d.square().mul(&d).div(&c(3.0))?);
    // F(v0 + m w) = F(v0) + rho0 m w + (s m^2 / 2) w^2.
    let at_start = if poled {
        rho0.square().div(&sa.mul(&c(2.0)))?
    } else {
        radius.mul(&a[1]).add(&sa.mul(&a[1].square()).mul(&c(0.5)))
    };
    let first = at_start.mul(&k0);
    let second = rho0.mul(&m).mul(&k1);
    let third = sa.mul(&m.square()).mul(&c(0.5)).mul(&k2);
    Some(first.add(&second).add(&third).neg())
}

/// The integral over the face of (S - origin).(S_u x S_v) du dv, as -loop
/// integral of v f(u) du (f does not depend on v), loops closed by chords;
/// their orientation carries the face's sense. A shell's faces share the
/// origin, a point of the shell, so a far body's fluxes do not cancel.
fn face_flux<T: Real>(face: &Face, loops: &[Lp], origin: &V3<T>) -> Option<T> {
    if let Surface::BSpline(spline) = &face.surface {
        // Exact on a nonrational surface with pieces inside patches,
        // enclosed otherwise.
        return spline_flux::spline_face_flux(spline, loops, origin)
            .or_else(|| spline_flux::enclosed_face_flux(spline, loops, origin));
    }
    if matches!(face.surface, Surface::Sphere { .. } | Surface::Torus { .. }) {
        return mass::sphere_flux(face, loops, origin);
    }
    if let Surface::Cone {
        frame: f,
        radius,
        half_angle,
    } = &face.surface
    {
        let fr = frame::<T>(f);
        let o = vsub(&fr.o, origin);
        let (ca, sa) = T::cos_sin(&c(*half_angle));
        let rad = c::<T>(*radius);
        let h = (
            ca.mul(&vdot(&o, &fr.x)),
            ca.mul(&vdot(&o, &fr.y)),
            ca.mul(&rad).sub(&sa.mul(&vdot(&o, &fr.n))),
        );
        // A sound face with a pole winds once in total (pole_position).
        let poled = loops.iter().map(|lp| lp.winding).sum::<i32>().abs() == 1;
        let mut total = c::<T>(0.0);
        for lp in loops {
            for u in &lp.fins {
                let (start, end) = match &u.pcurve {
                    Curve2::LineSegment { start, end } => (start, end),
                    // -∫ F(v) h(u) du along a spline, enclosed (S4d).
                    Curve2::BSpline(spline) => {
                        let (ha, hb, hc) = &h;
                        let g = |uu: &T, v: &T| -> Option<T> {
                            let (co, si) = T::cos_sin(uu);
                            let hu = ha.mul(&co).add(&hb.mul(&si)).add(hc);
                            let fv = if poled {
                                rad.add(&sa.mul(v)).square().div(&sa.mul(&c(2.0)))?
                            } else {
                                rad.mul(v).add(&sa.mul(&v.square()).mul(&c(0.5)))
                            };
                            Some(fv.mul(&hu))
                        };
                        total =
                            total.add(&bernstein::green_integral(spline, SPLINE_DEPTH, &g)?.neg());
                        continue;
                    }
                    // -∫ F(v) h(u) du along a projection (S8d.2).
                    Curve2::Projection(pr) => {
                        let (ha, hb, hc) = &h;
                        let term = projection::integrate_along::<T>(pr, &|u, v, du, _| {
                            let (co, si) = u.cos_sin();
                            let hu = co.scale(ha).add(&si.scale(hb)).add_constant(hc);
                            let fv = if poled {
                                v.scale(&sa)
                                    .add_constant(&rad)
                                    .square()
                                    .scale(&c::<T>(1.0).div(&sa.mul(&c(2.0)))?)
                            } else {
                                v.scale(&rad).add(&v.square().scale(&sa.mul(&c(0.5))))
                            };
                            Some(fv.mul(&hu).mul(du).neg())
                        })?;
                        total = total.add(&term);
                        continue;
                    }
                    Curve2::CircularArc { .. }
                    | Curve2::EllipseArc { .. }
                    | Curve2::Sinusoid { .. } => return None,
                };
                let term = cone_line_flux(
                    &[c(start.x), c(start.y)],
                    &[c(end.x), c(end.y)],
                    &sa,
                    &rad,
                    &h,
                    poled,
                )?;
                total = total.add(&term);
            }
            for (a, b) in chords::<T>(lp) {
                let term = match cone_line_flux(&a, &b, &sa, &rad, &h, poled) {
                    Some(term) => term,
                    None => {
                        let (ha, hb, hc) = &h;
                        let g = |uu: &T, v: &T| -> Option<T> {
                            let (co, si) = T::cos_sin(uu);
                            let hu = ha.mul(&co).add(&hb.mul(&si)).add(hc);
                            let fv = if poled {
                                rad.add(&sa.mul(v)).square().div(&sa.mul(&c(2.0)))?
                            } else {
                                rad.mul(v).add(&sa.mul(&v.square()).mul(&c(0.5)))
                            };
                            Some(fv.mul(&hu))
                        };
                        chord_enclosure(&a, &b, &g)?
                    }
                };
                total = total.add(&term);
            }
        }
        return Some(total);
    }
    let (plane, coeffs) = match &face.surface {
        Surface::Plane(f) => {
            let fr = frame::<T>(f);
            let o = vsub(&fr.o, origin);
            (true, (c(0.0), c(0.0), vdot(&o, &vcross(&fr.x, &fr.y))))
        }
        Surface::Cylinder { frame: f, radius } => {
            let fr = frame::<T>(f);
            let o = vsub(&fr.o, origin);
            let rad = c::<T>(*radius);
            let a = vdot(&o, &vcross(&fr.x, &fr.n));
            let b = vdot(&o, &vcross(&fr.y, &fr.n));
            let det = vdot(&fr.x, &vcross(&fr.y, &fr.n));
            (
                false,
                (rad.mul(&a).neg(), rad.mul(&b), rad.square().mul(&det)),
            )
        }
        Surface::Cone { .. }
        | Surface::Sphere { .. }
        | Surface::Torus { .. }
        | Surface::BSpline(_) => unreachable!("handled above"),
    };
    let mut total = c::<T>(0.0);
    for lp in loops {
        // On a plane `∮ du = 0` around the closed loop, so `v` is taken
        // relative to a point of the loop: a loop far from the UV origin
        // keeps its enclosures at its own scale.
        let v0: R = match lp.fins.first() {
            Some(u) if plane => pcurve_at::<T>(&u.pcurve, 0.0)[1].midpoint(),
            _ => R::from_integer(0.into()),
        };
        let shift = T::from_r(&v0);
        for u in &lp.fins {
            let term = match &u.pcurve {
                Curve2::LineSegment { start, end } => line_flux(
                    &[c(start.x), c::<T>(start.y).sub(&shift)],
                    &[c(end.x), c::<T>(end.y).sub(&shift)],
                    &coeffs,
                    plane,
                )?,
                Curve2::CircularArc {
                    center,
                    radius,
                    start_angle,
                    sweep_angle,
                } if plane => {
                    // h integral of (cy + rho sin t) rho sin t dt.
                    let a0 = c::<T>(*start_angle);
                    let a1 = a0.add(&c(*sweep_angle));
                    let (c0, _) = T::cos_sin(&a0);
                    let (c1, _) = T::cos_sin(&a1);
                    let (_, t0) = T::cos_sin(&a0.mul(&c(2.0)));
                    let (_, t1) = T::cos_sin(&a1.mul(&c(2.0)));
                    let rho = c::<T>(*radius);
                    let first = c::<T>(center.y).sub(&shift).mul(&rho).mul(&c0.sub(&c1));
                    let second = rho.square().mul(
                        &c::<T>(*sweep_angle)
                            .mul(&c(0.5))
                            .sub(&t1.sub(&t0).mul(&c(0.25))),
                    );
                    coeffs.2.mul(&first.add(&second))
                }
                // h integral of (cy + ry sin t) rx sin t dt.
                Curve2::EllipseArc {
                    center,
                    major,
                    minor,
                    start_angle,
                    sweep_angle,
                } if plane => {
                    let a0 = c::<T>(*start_angle);
                    let a1 = a0.add(&c(*sweep_angle));
                    let (c0, _) = T::cos_sin(&a0);
                    let (c1, _) = T::cos_sin(&a1);
                    let (_, t0) = T::cos_sin(&a0.mul(&c(2.0)));
                    let (_, t1) = T::cos_sin(&a1.mul(&c(2.0)));
                    let (rx, ry) = (c::<T>(*major), c::<T>(*minor));
                    let first = c::<T>(center.y).sub(&shift).mul(&rx).mul(&c0.sub(&c1));
                    let second = rx.mul(&ry).mul(
                        &c::<T>(*sweep_angle)
                            .mul(&c(0.5))
                            .sub(&t1.sub(&t0).mul(&c(0.25))),
                    );
                    coeffs.2.mul(&first.add(&second))
                }
                // -∫ v (A sin u + B cos u + C) du on a cylinder, in closed form.
                Curve2::Sinusoid { start, sweep, a } if !plane => {
                    sinusoid_flux(*start, *sweep, a, &coeffs)
                }
                // On a plane, `h (-∮ v du)`.
                Curve2::BSpline(spline) if plane => {
                    let about = [R::from_integer(0.into()), v0.clone()];
                    coeffs
                        .2
                        .mul(&bernstein::minus_v_du_about::<T>(spline, &about)?)
                }
                // A spline pcurve of constant u adds nothing to -∮ v f(u) du.
                Curve2::BSpline(spline) if constant_u(spline) => c(0.0),
                // Otherwise -∫ v f(u) du, enclosed (S4d).
                Curve2::BSpline(spline) => {
                    let f = |u: &T, v: &T| -> Option<T> {
                        let (co, si) = T::cos_sin(u);
                        Some(v.mul(&coeffs.0.mul(&si).add(&coeffs.1.mul(&co)).add(&coeffs.2)))
                    };
                    bernstein::green_integral(spline, SPLINE_DEPTH, &f)?.neg()
                }
                // -∫ (v - v0) f(u) du along a projection (S8d.2): `f` the
                // plane's constant or the cylinder's `A sin u + B cos u + C`.
                Curve2::Projection(pr) => projection::integrate_along::<T>(pr, &|u, v, du, _| {
                    let f = if plane {
                        Jet::constant(coeffs.2.clone(), u.order())
                    } else {
                        let (co, si) = u.cos_sin();
                        si.scale(&coeffs.0)
                            .add(&co.scale(&coeffs.1))
                            .add_constant(&coeffs.2)
                    };
                    Some(v.add_constant(&shift.neg()).mul(&f).mul(du).neg())
                })?,
                Curve2::CircularArc { .. }
                | Curve2::EllipseArc { .. }
                | Curve2::Sinusoid { .. } => return None,
            };
            total = total.add(&term);
        }
        for (a, b) in chords::<T>(lp) {
            let (a, b) = (
                [a[0].clone(), a[1].sub(&shift)],
                [b[0].clone(), b[1].sub(&shift)],
            );
            let term = match line_flux(&a, &b, &coeffs, plane) {
                Some(term) => term,
                None => {
                    let g = |u: &T, v: &T| -> Option<T> {
                        let (co, si) = T::cos_sin(u);
                        let f = if plane {
                            coeffs.2.clone()
                        } else {
                            coeffs.0.mul(&si).add(&coeffs.1.mul(&co)).add(&coeffs.2)
                        };
                        Some(v.mul(&f))
                    };
                    chord_enclosure(&a, &b, &g)?
                }
            };
            total = total.add(&term);
        }
    }
    Some(total)
}

/// `-∫ v (A sin u + B cos u + C) du` along `v = a0 + a1 cos u + a2 sin u`
/// from `start` over `sweep`, `(A, B, C)` the cylinder's flux coefficients.
fn sinusoid_flux<T: Real>(start: f64, sweep: f64, a: &[f64; 3], k: &(T, T, T)) -> T {
    let u0 = c::<T>(start);
    let u1 = u0.add(&c(sweep));
    let ((c0, s0), (c1, s1)) = (T::cos_sin(&u0), T::cos_sin(&u1));
    let (_, t0) = T::cos_sin(&u0.mul(&c(2.0)));
    let (_, t1) = T::cos_sin(&u1.mul(&c(2.0)));
    let du = c::<T>(sweep);
    // The integrals of sin, cos, sin cos, cos^2 and sin^2 between the ends.
    let i_s = c0.sub(&c1);
    let i_c = s1.sub(&s0);
    let i_sc = s1.square().sub(&s0.square()).mul(&c(0.5));
    let quarter = t1.sub(&t0).mul(&c(0.25));
    let i_cc = du.mul(&c(0.5)).add(&quarter);
    let i_ss = du.mul(&c(0.5)).sub(&quarter);
    let (a0, a1, a2) = (c::<T>(a[0]), c::<T>(a[1]), c::<T>(a[2]));
    let (ka, kb, kc) = k;
    a0.mul(kc)
        .mul(&du)
        .add(&a0.mul(ka).mul(&i_s))
        .add(&a0.mul(kb).mul(&i_c))
        .add(&a1.mul(kc).mul(&i_c))
        .add(&a1.mul(ka).mul(&i_sc))
        .add(&a1.mul(kb).mul(&i_cc))
        .add(&a2.mul(kc).mul(&i_s))
        .add(&a2.mul(ka).mul(&i_ss))
        .add(&a2.mul(kb).mul(&i_sc))
        .neg()
}

/// The faces a shell lists on both sides: a sheet's (S6). They bound no
/// region of the shell's own and add nothing to its flux.
fn two_sided(shell: &Shell) -> BTreeSet<usize> {
    let sides: BTreeSet<(usize, Side)> = shell.sides.iter().map(|(f, s)| (f.0, *s)).collect();
    shell
        .sides
        .iter()
        .filter(|(f, side)| {
            let other = match side {
                Side::Front => Side::Back,
                Side::Back => Side::Front,
            };
            sides.contains(&(f.0, other))
        })
        .map(|(f, _)| f.0)
        .collect()
}

/// The origin of a shell's fluxes: any point works (the divergence of
/// `S - origin` is 3), and a vertex of the shell keeps every face's
/// integrand at the body's own scale; a shell without vertices takes its
/// first surface's frame origin or first pole.
fn shell_origin(view: &View, resolved: &[Vec<Lp>], shell: &Shell) -> [f64; 3] {
    for (f, _) in &shell.sides {
        for lp in &resolved[f.0] {
            for fin in &lp.fins {
                let edge = &view.edges[fin.edge.0];
                if let Some(v) = edge.start.or(edge.end) {
                    return view.vertices[v.0].position.to_array();
                }
            }
        }
    }
    match shell.sides.first().map(|(f, _)| &view.faces[f.0].surface) {
        Some(
            Surface::Plane(f)
            | Surface::Cylinder { frame: f, .. }
            | Surface::Cone { frame: f, .. }
            | Surface::Sphere { frame: f, .. }
            | Surface::Torus { frame: f, .. },
        ) => f.origin().to_array(),
        Some(Surface::BSpline(s)) => s.poles()[0].to_array(),
        None => [0.0; 3],
    }
}

/// A shell's flux: + for front sides, - for back sides. Positive when the
/// shell encloses its region.
fn shell_flux<T: Real>(
    faces: &[Face],
    resolved: &[Vec<Lp>],
    shell: &Shell,
    origin: [f64; 3],
) -> Option<T> {
    let origin = v3::<T>(origin);
    let mut total = c::<T>(0.0);
    let both = two_sided(shell);
    for (f, side) in &shell.sides {
        if both.contains(&f.0) {
            continue;
        }
        let flux = face_flux::<T>(&faces[f.0], &resolved[f.0], &origin)?;
        total = match side {
            Side::Front => total.add(&flux),
            Side::Back => total.sub(&flux),
        };
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

/// Parity of ray hits from an exact point against a shell's faces; None when
/// every ray direction meets an uncertified decision. Adjacent faces only
/// meet within tolerance, so a hit also needs a certified margin
/// `sqrt(margin2)` from its face's boundary. The parity is then the same for
/// any watertight surface within tolerance of the faces.
fn inside<T: Real>(
    faces: &[Face],
    resolved: &[Vec<Lp>],
    shell: &Shell,
    point: &[R; 3],
    margin2: &T,
) -> Option<bool> {
    'direction: for d in DIRECTIONS {
        let d: [R; 3] = d.map(int);
        let mut hits = 0;
        for (f, _) in &shell.sides {
            match face_hits::<T>(&faces[f.0], &resolved[f.0], point, &d, margin2) {
                Some(n) => hits += n,
                None => continue 'direction,
            }
        }
        return Some(hits % 2 == 1);
    }
    None
}

/// A ray's crossings of a whole torus: the positive real roots of
/// `(|w|^2 + R^2 - r^2)^2 - 4 R^2 (w_x^2 + w_y^2)`, `w = q0 + t q1` in the
/// frame's coordinates (Cramer on its stored axes), counted exactly; `None`
/// for a root at the start or a repeated root (a tangent ray).
fn torus_ray_hits(frame: &Frame3, major: f64, minor: f64, p: &[R; 3], d: &[R; 3]) -> Option<u32> {
    use crate::polynomial::real::{isolate, Budget, IntPolynomial};
    let exact = |v: [f64; 3]| v.map(r);
    let (o, x, y, n) = (
        exact(frame.origin().to_array()),
        exact(frame.x().to_array()),
        exact(frame.y().to_array()),
        exact(frame.normal().to_array()),
    );
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
    let (big, small) = (r(major), r(minor));
    let k = &big * &big - &small * &small;
    let dotr = |a: &[R; 3], b: &[R; 3]| &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2];
    let (s2, s1, s0) = (dotr(&q1, &q1), int(2) * dotr(&q0, &q1), dotr(&q0, &q0) + &k);
    let (p2, p1, p0) = (
        &q1[0] * &q1[0] + &q1[1] * &q1[1],
        int(2) * (&q0[0] * &q1[0] + &q0[1] * &q1[1]),
        &q0[0] * &q0[0] + &q0[1] * &q0[1],
    );
    let b4 = int(4) * &big * &big;
    let poly = [
        &s0 * &s0 - &b4 * &p0,
        int(2) * &s0 * &s1 - &b4 * &p1,
        &s1 * &s1 + int(2) * &s0 * &s2 - &b4 * &p2,
        int(2) * &s1 * &s2,
        &s2 * &s2,
    ];
    if poly[0].numer().sign() == Sign::NoSign {
        return None;
    }
    let ip = IntPolynomial::from_rationals(&poly);
    if !ip.gcd(&ip.derivative()).is_constant() {
        return None;
    }
    let abs = |c: &R| {
        if c.numer().sign() == Sign::Minus {
            -c.clone()
        } else {
            c.clone()
        }
    };
    let lead = abs(&poly[4]);
    let bound = poly[..4]
        .iter()
        .map(|c| abs(c) / &lead)
        .fold(int(0), |m, c| if c > m { c } else { m })
        + int(1);
    let roots = isolate(
        &ip,
        int(0),
        bound,
        &mut Budget::new(crate::polynomial::RootIsolationOptions::default()),
    )
    .ok()?;
    Some(roots.len() as u32)
}

fn det3(a: &[R; 3], b: &[R; 3], cc: &[R; 3]) -> R {
    &a[0] * (&b[1] * &cc[2] - &b[2] * &cc[1]) - &a[1] * (&b[0] * &cc[2] - &b[2] * &cc[0])
        + &a[2] * (&b[0] * &cc[1] - &b[1] * &cc[0])
}

fn face_hits<T: Real>(
    face: &Face,
    loops: &[Lp],
    p: &[R; 3],
    d: &[R; 3],
    margin2: &T,
) -> Option<u32> {
    // A whole sphere (no loops, S9d.1): the ray's hits are the positive
    // roots of `|p + t d - c|^2 = R^2`, exactly (a tangency or a root at
    // the start undecided).
    if let (Surface::Sphere { frame, radius }, true) = (&face.surface, loops.is_empty()) {
        let c = frame.origin().to_array().map(r);
        let rr = r(*radius);
        let w: [R; 3] = std::array::from_fn(|i| &p[i] - &c[i]);
        let dotr = |a: &[R; 3], b: &[R; 3]| &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2];
        let (a2, b2, c2) = (dotr(d, d), int(2) * dotr(&w, d), dotr(&w, &w) - &rr * &rr);
        let disc = &b2 * &b2 - int(4) * &a2 * &c2;
        let zero = int(0);
        if disc < zero {
            return Some(0);
        }
        if disc == zero || c2 == zero {
            return None;
        }
        // Roots' product c2 / a2 and sum -b2 / a2 (a2 > 0).
        return Some(if c2 < zero {
            1
        } else if b2 < zero {
            2
        } else {
            0
        });
    }
    // A whole torus (no loops, S9d.4a): the positive roots of the ray's
    // quartic, exactly (a root at the start or a repeated one undecided).
    if let (
        Surface::Torus {
            frame,
            major,
            minor,
        },
        true,
    ) = (&face.surface, loops.is_empty())
    {
        return torus_ray_hits(frame, *major, *minor, p, d);
    }
    // Rays against spheres with loops and tori with loops are not decided
    // yet: the containment is uncertified, unless the ray misses the whole
    // surface (S9d.4b.1: a cavity in a torus segment's barrel).
    if let Surface::Torus {
        frame,
        major,
        minor,
    } = &face.surface
    {
        return (torus_ray_hits(frame, *major, *minor, p, d)? == 0).then_some(0);
    }
    if let Surface::Sphere { frame, radius } = &face.surface {
        let c = frame.origin().to_array().map(r);
        let rr = r(*radius);
        let w: [R; 3] = std::array::from_fn(|i| &p[i] - &c[i]);
        let dotr = |a: &[R; 3], b: &[R; 3]| &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2];
        let (a2, b2, c2) = (dotr(d, d), int(2) * dotr(&w, d), dotr(&w, &w) - &rr * &rr);
        let zero = int(0);
        // No real root, or both behind the start (their product positive,
        // their sum negative).
        let miss = &b2 * &b2 < int(4) * &a2 * &c2 || (c2 > zero && b2 > zero);
        return miss.then_some(0);
    }
    if matches!(face.surface, Surface::BSpline(_)) {
        return None;
    }
    let refs: Vec<&Lp> = loops.iter().collect();
    let exact = |f: &Frame3| {
        (
            f.origin().to_array().map(r),
            f.x().to_array().map(r),
            f.y().to_array().map(r),
            f.normal().to_array().map(r),
        )
    };
    match &face.surface {
        Surface::Sphere { .. } | Surface::Torus { .. } | Surface::BSpline(_) => None,
        // S9d.3a: `X^2 + Y^2 = (R + Z tan a)^2` in the frame, the hits on
        // the nappe where `R + Z tan a > 0`; each hit's `(u, v)` against the
        // loops by the +v ray's crossings, the apex (where the face closes
        // at it, its loops winding once) one more when it lies up the ray.
        Surface::Cone {
            frame: f,
            radius,
            half_angle,
        } => {
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
            let (ca, sa) = T::cos_sin(&c(*half_angle));
            let tan = sa.div(&ca)?;
            let rad = c::<T>(*radius);
            let [x0, y0, z0] = q0.clone().map(|v| q::<T>(&v));
            let [x1, y1, z1] = q1.clone().map(|v| q::<T>(&v));
            let rho0 = rad.add(&z0.mul(&tan));
            let k1 = z1.mul(&tan);
            let a = x1.square().add(&y1.square()).sub(&k1.square());
            let b = x0
                .mul(&x1)
                .add(&y0.mul(&y1))
                .sub(&rho0.mul(&k1))
                .mul(&c(2.0));
            let cc = x0.square().add(&y0.square()).sub(&rho0.square());
            let a_sign = a.sign()?;
            if a_sign == Ordering::Equal {
                return None;
            }
            let disc = b.square().sub(&a.mul(&cc).mul(&c(4.0)));
            match disc.sign()? {
                Ordering::Less => return Some(0),
                Ordering::Equal => return None,
                Ordering::Greater => {}
            }
            let root = disc.sqrt();
            let apex_up = loops.iter().map(|lp| lp.winding).sum::<i32>().abs() == 1
                && sa.sign()? == Ordering::Less;
            let mut count = 0;
            for sign in [-1.0, 1.0] {
                let t = b.neg().add(&root.mul(&c(sign))).div(&a.mul(&c(2.0)))?;
                match t.sign()? {
                    Ordering::Less => continue,
                    Ordering::Equal => return None,
                    Ordering::Greater => {}
                }
                let (hx, hy, hz) = (
                    x0.add(&t.mul(&x1)),
                    y0.add(&t.mul(&y1)),
                    z0.add(&t.mul(&z1)),
                );
                let rho = rad.add(&hz.mul(&tan));
                match rho.sign()? {
                    // The other nappe: not the surface's.
                    Ordering::Less => continue,
                    Ordering::Equal => return None,
                    Ordering::Greater => {}
                }
                let hit = [T::atan2(&hy, &hx)?, hz.div(&ca)?];
                if !clear_of_boundary(&refs, &hit, &rho, true, margin2)? {
                    return None;
                }
                let crossed = cover_crossings::<T>(&refs, &hit)? + u32::from(apex_up);
                count += crossed % 2;
            }
            Some(count)
        }
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
            if !clear_of_boundary(&refs, &hit, &c(1.0), false, margin2)? {
                return None;
            }
            let inside_face = crossings::<T>(&refs, &hit)? % 2;
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
                // Any representative angle: the loops are tested on every alias.
                let hit = [T::atan2(&beta, &alpha)?, v];
                if !clear_of_boundary(&refs, &hit, &rad_t, true, margin2)? {
                    return None;
                }
                count += cover_crossings::<T>(&refs, &hit)? % 2;
            }
            Some(count)
        }
    }
}

/// Whether a UV point is certainly farther than `sqrt(margin2)` from every
/// fin and closing chord of the loops, with u distances scaled by `su` (the
/// cylinder radius; one on planes) and, on a periodic surface, every u alias
/// of each piece considered. Some(false) when it may be closer.
fn clear_of_boundary<T: Real>(
    loops: &[&Lp],
    p: &V2<T>,
    su: &T,
    periodic: bool,
    margin2: &T,
) -> Option<bool> {
    let scaled = |a: &V2<T>| [a[0].mul(su), a[1].clone()];
    let far = |d2: &T| d2.cmp(margin2) == Some(Ordering::Greater);
    // Squared distance from p to segment ab exceeds margin2.
    let segment_clear = |a: &V2<T>, b: &V2<T>, p: &V2<T>| -> Option<bool> {
        let p = scaled(p);
        let dist2 = |x: &V2<T>| {
            let (du, dv) = (x[0].sub(&p[0]), x[1].sub(&p[1]));
            du.square().add(&dv.square())
        };
        let (a, b) = (scaled(a), scaled(b));
        let d = [b[0].sub(&a[0]), b[1].sub(&a[1])];
        let w = [p[0].sub(&a[0]), p[1].sub(&a[1])];
        let len2 = d[0].square().add(&d[1].square());
        match len2.sign() {
            Some(Ordering::Greater) => {}
            // A chord of no length within rounding (a ring's, closing its
            // one fin): clear when its start is farther than the margin by
            // more than its length.
            _ => {
                let gap = dist2(&a).sqrt().sub(&len2.sqrt());
                return Some(gap.sign()? == Ordering::Greater && far(&gap.square()));
            }
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
    let each_alias = |a: &V2<T>, b: &V2<T>| -> Option<bool> {
        if !periodic {
            return segment_clear(a, b, p);
        }
        for k in alias_range(&a[0], &b[0], &p[0])? {
            let shifted = [p[0].sub(&c(TAU * k as f64)), p[1].clone()];
            if !segment_clear(a, b, &shifted)? {
                return Some(false);
            }
        }
        Some(true)
    };
    for lp in loops {
        for (a, b) in chords::<T>(lp) {
            if !each_alias(&a, &b)? {
                return Some(false);
            }
        }
        for u in &lp.fins {
            let clear = match &u.pcurve {
                Curve2::LineSegment { start, end } => {
                    each_alias(&[c(start.x), c(start.y)], &[c(end.x), c(end.y)])?
                }
                // Plane arcs only: clear of the whole circle.
                Curve2::CircularArc { .. } if periodic => return None,
                Curve2::CircularArc { center, radius, .. } => {
                    let d = [p[0].sub(&c(center.x)), p[1].sub(&c(center.y))];
                    let r0 = d[0].square().add(&d[1].square()).sqrt();
                    far(&r0.sub(&c(*radius)).square())
                }
                // An ellipse is the image of the unit circle under scalings
                // by its semi-axes, which shrink distances by at most the
                // minor one: clear beyond the margin over it.
                Curve2::EllipseArc { .. } if periodic => return None,
                Curve2::EllipseArc {
                    center,
                    major,
                    minor,
                    ..
                } => {
                    let d = [
                        p[0].sub(&c(center.x)).div(&c(*major))?,
                        p[1].sub(&c(center.y)).div(&c(*minor))?,
                    ];
                    let r0 = d[0].square().add(&d[1].square()).sqrt();
                    far(&r0.sub(&c(1.0)).square().mul(&c::<T>(*minor).square()))
                }
                Curve2::BSpline(_) | Curve2::Sinusoid { .. } | Curve2::Projection(_) => {
                    return None
                }
            };
            if !clear {
                return Some(false);
            }
        }
    }
    Some(true)
}
// ------------------------------------------------------------------ the contract

fn closed_curve(curve: &Curve3) -> bool {
    match curve {
        Curve3::Circle { .. } => true,
        Curve3::CircularArc { sweep_angle, .. } | Curve3::EllipseArc { sweep_angle, .. } => {
            sweep_angle.abs() == TAU
        }
        Curve3::HyperbolaArc { .. } | Curve3::ParabolaArc { .. } => false,
        // A torus section over a whole turn of its angle (S8d.3).
        Curve3::Section(s) => s.sweep.abs() == TAU,
        // Two cylinders' meeting over a whole turn: a ring (S9c.2).
        Curve3::Meet(m) => m.sweep.abs() == TAU,
        // A meeting over the height never closes (S9d.2b).
        Curve3::Rise(_) => false,
        // A torus's meeting over a whole turn of its angle (S9d.4b.2).
        Curve3::Toric(m) => m.sweep.abs() == TAU,
        // A spline wall's meeting over its wall's `u` never closes (S9f.2b).
        Curve3::WallMeet(_) => false,
        Curve3::LineSegment { .. } => false,
        // A spline ring edge is a full period; its seam is tested for C1.
        Curve3::BSpline(span) => span.is_closed_period(),
    }
}

pub(crate) fn check(view: &View, tolerance: Tolerance) -> Vec<Issue> {
    bernstein::clear_memo();
    let View {
        vertices,
        edges,
        fins,
        loops,
        faces,
        shells,
        regions,
    } = *view;
    let mut issues: BTreeSet<Issue> = BTreeSet::new();
    let add = |issues: &mut BTreeSet<Issue>, kind, entity| {
        issues.insert(Issue { kind, entity });
    };
    use Entity as En;
    use IssueKind as K;
    let (nv, ne, nfin, nl, nf, ns, nr) = (
        vertices.len(),
        edges.len(),
        fins.len(),
        loops.len(),
        faces.len(),
        shells.len(),
        regions.len(),
    );
    // Fin and loop locations (by first occurrence) name the entities.
    let mut place: BTreeMap<usize, (usize, usize, usize)> = BTreeMap::new();
    let mut loop_place: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
    for (fi, face) in faces.iter().enumerate() {
        for (li, l) in face.loops.iter().enumerate() {
            if l.0 < nl {
                loop_place.entry(l.0).or_insert((fi, li));
                if let Loop::Edges { fins: list, .. } = &loops[l.0] {
                    for (ui, k) in list.iter().enumerate() {
                        if k.0 < nfin {
                            place.entry(k.0).or_insert((fi, li, ui));
                        }
                    }
                }
            }
        }
    }
    let fin_name = |k: usize| match place.get(&k) {
        Some(&(f, l, u)) => En::Use(f, l, u),
        None => En::FinSlot(k),
    };
    let loop_name = |l: usize| match loop_place.get(&l) {
        Some(&(f, i)) => En::Loop(f, i),
        None => En::LoopSlot(l),
    };

    // ------------------------------------------------ references
    for (i, edge) in edges.iter().enumerate() {
        if [edge.start, edge.end].iter().flatten().any(|v| v.0 >= nv)
            || edge.fins.iter().any(|k| k.0 >= nfin)
        {
            add(&mut issues, K::Reference, En::Edge(i));
        }
    }
    for (k, fin) in fins.iter().enumerate() {
        if fin.edge.0 >= ne {
            add(&mut issues, K::Reference, fin_name(k));
        }
    }
    for (l, lp) in loops.iter().enumerate() {
        let bad = match lp {
            Loop::Edges { fins: list, .. } => list.iter().any(|k| k.0 >= nfin),
            Loop::Vertex(v) => v.0 >= nv,
        };
        if bad {
            add(&mut issues, K::Reference, loop_name(l));
        }
    }
    for (fi, face) in faces.iter().enumerate() {
        if face.loops.iter().any(|l| l.0 >= nl) || face.front.0 >= ns || face.back.0 >= ns {
            add(&mut issues, K::Reference, En::Face(fi));
        }
    }
    for (si, shell) in shells.iter().enumerate() {
        if shell.region.0 >= nr
            || shell.sides.iter().any(|(f, _)| f.0 >= nf)
            || shell.wire_edges.iter().any(|e| e.0 >= ne)
            || shell.acorn_vertices.iter().any(|v| v.0 >= nv)
        {
            add(&mut issues, K::Reference, En::Shell(si));
        }
    }
    for (ri, region) in regions.iter().enumerate() {
        if region.shells.iter().any(|s| s.0 >= ns) {
            add(&mut issues, K::Reference, En::Region(ri));
        }
    }
    if !issues.is_empty() {
        return issues.into_iter().collect();
    }
    let tol = r(tolerance.linear());
    let tol2 = &tol * &tol;
    let (fast_tol2, exact_tol2) = (Fast::from_r(&tol2), I::from_r(&tol2));
    // Containment rays keep twice the tolerance from face boundaries.
    let (fast_margin2, exact_margin2) = (fast_tol2.mul(&c(4.0)), exact_tol2.mul(&c(4.0)));

    // ------------------------------------------------ structure (exact)
    let mut fin_count = vec![0usize; nfin];
    for lp in loops {
        if let Loop::Edges { fins: list, .. } = lp {
            for k in list {
                fin_count[k.0] += 1;
            }
        }
    }
    for (k, n) in fin_count.iter().enumerate() {
        match n {
            0 => add(&mut issues, K::FinWithoutLoop, fin_name(k)),
            1 => {}
            _ => add(&mut issues, K::FinReused, fin_name(k)),
        }
    }
    let mut loop_count = vec![0usize; nl];
    for face in faces {
        for l in &face.loops {
            loop_count[l.0] += 1;
        }
    }
    for (l, n) in loop_count.iter().enumerate() {
        match n {
            0 => add(&mut issues, K::LoopWithoutFace, loop_name(l)),
            1 => {}
            _ => add(&mut issues, K::LoopReused, loop_name(l)),
        }
    }
    let mut users: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (k, fin) in fins.iter().enumerate() {
        users.entry(fin.edge.0).or_default().push(k);
    }
    let sorted = |v: &[usize]| {
        let mut v = v.to_vec();
        v.sort_unstable();
        v
    };
    for (i, edge) in edges.iter().enumerate() {
        let listed: Vec<usize> = edge.fins.iter().map(|k| k.0).collect();
        if sorted(&listed) != sorted(users.get(&i).map_or(&[][..], Vec::as_slice)) {
            add(&mut issues, K::EdgeFinsMismatch, En::Edge(i));
        }
    }
    let mut used_vertex = vec![false; nv];
    for edge in edges {
        for v in [edge.start, edge.end].into_iter().flatten() {
            used_vertex[v.0] = true;
        }
    }
    for lp in loops {
        if let Loop::Vertex(v) = lp {
            used_vertex[v.0] = true;
        }
    }
    for shell in shells {
        for v in &shell.acorn_vertices {
            used_vertex[v.0] = true;
        }
    }
    for (v, used) in used_vertex.iter().enumerate() {
        if !used {
            add(&mut issues, K::UnusedVertex, En::Vertex(v));
        }
    }
    let wire: BTreeSet<usize> = shells
        .iter()
        .flat_map(|s| s.wire_edges.iter().map(|e| e.0))
        .collect();
    for (i, edge) in edges.iter().enumerate() {
        if !users.contains_key(&i) && !wire.contains(&i) {
            add(&mut issues, K::UnusedEdge, En::Edge(i));
        }
        if edge.start.is_none() != edge.end.is_none() {
            add(&mut issues, K::RingEdgeWithVertex, En::Edge(i));
        } else if edge.start.is_none() && !closed_curve(&edge.curve) {
            add(&mut issues, K::RingEdgeOpen, En::Edge(i));
        }
    }
    // Face sides against the shells listing them.
    let mut listings: BTreeMap<(usize, Side), Vec<usize>> = BTreeMap::new();
    for (si, shell) in shells.iter().enumerate() {
        if shell.sides.is_empty() && shell.wire_edges.is_empty() && shell.acorn_vertices.is_empty()
        {
            add(&mut issues, K::EmptyShell, En::Shell(si));
        }
        for (f, side) in &shell.sides {
            listings.entry((f.0, *side)).or_default().push(si);
        }
    }
    let mut side_bad: BTreeSet<usize> = BTreeSet::new();
    for (fi, face) in faces.iter().enumerate() {
        let none = Vec::new();
        let fr = listings.get(&(fi, Side::Front)).unwrap_or(&none);
        let bk = listings.get(&(fi, Side::Back)).unwrap_or(&none);
        let repeated = |v: &Vec<usize>| v.iter().any(|s| v.iter().filter(|t| *t == s).count() > 1);
        let kind = if fr.is_empty() && bk.is_empty() {
            Some(K::FaceWithoutShell)
        } else if fr.len() > 1 || bk.len() > 1 {
            Some(if repeated(fr) || repeated(bk) {
                K::FaceReused
            } else {
                K::SideInTwoShells
            })
        } else if fr.is_empty() || bk.is_empty() {
            Some(K::SideWithoutShell)
        } else if fr[0] != face.front.0 || bk[0] != face.back.0 {
            Some(K::SideRegionMismatch)
        } else {
            None
        };
        if let Some(kind) = kind {
            add(&mut issues, kind, En::Face(fi));
            side_bad.insert(fi);
        }
    }
    // Resolved loops per face; winding in u only.
    let resolved: Vec<Vec<Lp>> = faces
        .iter()
        .map(|face| {
            face.loops
                .iter()
                .filter_map(|l| match &loops[l.0] {
                    Loop::Edges {
                        fins: list,
                        winding,
                    } => Some(Lp {
                        fins: list.iter().map(|k| &fins[k.0]).collect(),
                        winding: winding[0],
                        winding_v: winding[1],
                    }),
                    Loop::Vertex(_) => None,
                })
                .collect()
        })
        .collect();
    let mut structural_faces = BTreeSet::new();
    for (fi, face) in faces.iter().enumerate() {
        // A sphere or torus without loops is the whole closed surface.
        if face.loops.is_empty()
            && !matches!(face.surface, Surface::Sphere { .. } | Surface::Torus { .. })
        {
            add(&mut issues, K::EmptyFace, En::Face(fi));
            structural_faces.insert(fi);
        }
        let periodic = face.surface.is_periodic();
        for (li, l) in face.loops.iter().enumerate() {
            let Loop::Edges {
                fins: list,
                winding,
            } = &loops[l.0]
            else {
                continue;
            };
            if list.is_empty() {
                add(&mut issues, K::EmptyLoop, En::Loop(fi, li));
                structural_faces.insert(fi);
                continue;
            }
            if (winding[1] != 0 && !face.surface.is_periodic_v()) || (!periodic && winding[0] != 0)
            {
                add(&mut issues, K::WindingMismatch, En::Loop(fi, li));
                structural_faces.insert(fi);
            }
            let ends: Vec<_> = list
                .iter()
                .map(|k| fin_vertices(edges, &fins[k.0]))
                .collect();
            if ends.len() == 1 && ends[0] == (None, None) {
                continue;
            }
            let open = ends.iter().any(|(a, b)| a.is_none() || b.is_none())
                || (0..ends.len()).any(|k| ends[k].1 != ends[(k + 1) % ends.len()].0);
            if open {
                add(&mut issues, K::OpenLoop, En::Loop(fi, li));
                structural_faces.insert(fi);
            }
        }
        if periodic {
            // Windings balance, except on a cone where one pole (a vertex
            // loop at the apex) closes a band that winds once.
            let wound: Vec<i32> = resolved[fi].iter().map(|lp| lp.winding).collect();
            let total = wound.iter().sum::<i32>();
            let poled = total.abs() == 1 && pole_position(face, loops).is_some();
            let total_v: i32 = resolved[fi].iter().map(|lp| lp.winding_v).sum();
            if (wound.iter().any(|w| *w != 0) && total != 0 && !poled) || total_v != 0 {
                add(&mut issues, K::WindingMismatch, En::Loop(fi, 0));
                structural_faces.insert(fi);
            }
        }
    }
    for (ri, region) in regions.iter().enumerate() {
        let distinct: BTreeSet<usize> = region.shells.iter().map(|s| s.0).collect();
        if distinct.len() != region.shells.len() {
            add(&mut issues, K::DoubleBounding, En::Region(ri));
        }
        if ri > 0 && region.shells.is_empty() {
            add(&mut issues, K::RegionWithoutShell, En::Region(ri));
        }
    }
    if regions.first().map(|r| r.kind) != Some(RegionKind::Void) {
        add(&mut issues, K::NoInfiniteRegion, En::Region(0));
    }
    for (si, shell) in shells.iter().enumerate() {
        if !regions[shell.region.0].shells.contains(&ShellId(si)) {
            add(&mut issues, K::RegionShellMismatch, En::Shell(si));
        }
    }
    // Seams are forbidden: an edge with two fins in one face.
    let fin_face: BTreeMap<usize, usize> = place.iter().map(|(k, p)| (*k, p.0)).collect();
    for (i, list) in &users {
        let faces_of: Vec<Option<&usize>> = list.iter().map(|k| fin_face.get(k)).collect();
        let distinct: BTreeSet<_> = faces_of.iter().collect();
        if distinct.len() != faces_of.len() {
            add(&mut issues, K::SeamEdge, En::Edge(*i));
            structural_faces.extend(list.iter().filter_map(|k| fin_face.get(k).copied()));
        }
    }

    // Edge accounting: shells must alternate around each edge's fins, in the
    // stored radial order.
    let mut bad_shells: BTreeSet<usize> = BTreeSet::new();
    for (i, edge) in edges.iter().enumerate() {
        let listed: Vec<usize> = edge.fins.iter().map(|k| k.0).collect();
        if listed.is_empty()
            || sorted(&listed) != sorted(users.get(&i).map_or(&[][..], Vec::as_slice))
        {
            continue;
        }
        let around: Vec<usize> = listed
            .iter()
            .copied()
            .filter(|k| fin_face.get(k).is_some_and(|f| !side_bad.contains(f)))
            .collect();
        if around.is_empty() {
            continue;
        }
        let face_of = |k: usize| &faces[fin_face[&k]];
        let forward = |k: usize| fins[k].sense == Orientation::Forward;
        let ahead: Vec<usize> = around
            .iter()
            .map(|&k| {
                if forward(k) {
                    face_of(k).front.0
                } else {
                    face_of(k).back.0
                }
            })
            .collect();
        let behind: Vec<usize> = around
            .iter()
            .map(|&k| {
                if forward(k) {
                    face_of(k).back.0
                } else {
                    face_of(k).front.0
                }
            })
            .collect();
        let n = around.len();
        if (0..n).all(|j| ahead[j] == behind[(j + 1) % n]) {
            continue;
        }
        let pairs: BTreeSet<(usize, usize)> = around
            .iter()
            .map(|&k| {
                let f = face_of(k);
                (f.front.0.min(f.back.0), f.front.0.max(f.back.0))
            })
            .collect();
        let kind = if pairs.len() > 1 {
            K::EdgeAcrossShells
        } else if n == 1 {
            K::FreeEdge
        } else if n == 2 {
            if forward(around[0]) == forward(around[1]) {
                K::SameSenseUses
            } else {
                K::RadialOrderInconsistent
            }
        } else {
            K::NonManifoldEdge
        };
        add(&mut issues, kind, En::Edge(i));
        for &k in &around {
            bad_shells.insert(face_of(k).front.0);
            bad_shells.insert(face_of(k).back.0);
        }
    }
    for &fi in &structural_faces {
        bad_shells.insert(faces[fi].front.0);
        bad_shells.insert(faces[fi].back.0);
    }
    // A shell listing exactly the opposite sides of an earlier shell is the
    // same surface: its shell-level checks are not repeated.
    let opposite = |side: Side| match side {
        Side::Front => Side::Back,
        Side::Back => Side::Front,
    };
    let twins: BTreeSet<usize> = (0..ns)
        .filter(|&si| {
            let mut mine: Vec<(usize, Side)> = shells[si]
                .sides
                .iter()
                .map(|(f, s)| (f.0, opposite(*s)))
                .collect();
            mine.sort();
            !mine.is_empty()
                && (0..si).any(|sj| {
                    let mut theirs: Vec<(usize, Side)> =
                        shells[sj].sides.iter().map(|(f, s)| (f.0, *s)).collect();
                    theirs.sort();
                    theirs == mine
                })
        })
        .collect();
    let shell_faces =
        |si: usize| -> Vec<usize> { shells[si].sides.iter().map(|(f, _)| f.0).collect() };
    for (si, shell) in shells.iter().enumerate().take(ns) {
        if twins.contains(&si) {
            continue;
        }
        let members_list = shell_faces(si);
        if bad_shells.contains(&si)
            || members_list.is_empty()
            || members_list.iter().any(|f| side_bad.contains(f))
        {
            bad_shells.insert(si);
            continue;
        }
        let members: BTreeSet<usize> = members_list.iter().copied().collect();
        let mut adjacency: BTreeMap<usize, BTreeSet<usize>> =
            members.iter().map(|f| (*f, BTreeSet::new())).collect();
        for list in users.values() {
            let fs: BTreeSet<usize> = list
                .iter()
                .filter_map(|k| fin_face.get(k).copied())
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
            for lp in &resolved[*f] {
                for (k, u) in lp.fins.iter().enumerate() {
                    let w = lp.fins[(k + 1) % lp.fins.len()];
                    let Some(v) = fin_vertices(edges, u).1 else {
                        continue;
                    };
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
        let mut loop_total = 0i64;
        for f in &members {
            for l in &faces[*f].loops {
                loop_total += 1;
                match &loops[l.0] {
                    Loop::Vertex(v) => {
                        vs.insert(v.0);
                    }
                    Loop::Edges { fins: list, .. } => {
                        for k in list {
                            let e = fins[k.0].edge.0;
                            if let (Some(a), Some(b)) = (edges[e].start, edges[e].end) {
                                vs.insert(a.0);
                                vs.insert(b.0);
                                es.insert(e);
                            }
                        }
                    }
                }
            }
        }
        let chi = vs.len() as i64 - es.len() as i64 + 2 * members.len() as i64 - loop_total;
        // A sheet's two sides make no closed surface of its faces (S6).
        if two_sided(shell).is_empty() && (chi % 2 != 0 || chi > 2) {
            add(&mut issues, K::Euler, En::Shell(si));
            bad_shells.insert(si);
        }
    }

    // Wire edges and acorn vertices (S6): a wire edge has no fins, an
    // acorn vertex bounds nothing else, and a shell of wire edges and
    // acorns alone is connected through their vertices.
    let edge_ends: BTreeSet<usize> = edges
        .iter()
        .flat_map(|e| [e.start, e.end])
        .flatten()
        .map(|v| v.0)
        .collect();
    let loop_vertices: BTreeSet<usize> = loops
        .iter()
        .filter_map(|l| match l {
            Loop::Vertex(v) => Some(v.0),
            Loop::Edges { .. } => None,
        })
        .collect();
    for (si, shell) in shells.iter().enumerate() {
        for e in &shell.wire_edges {
            if !edges[e.0].fins.is_empty() || users.contains_key(&e.0) {
                add(&mut issues, K::WireEdgeWithFins, En::Edge(e.0));
            }
        }
        for v in &shell.acorn_vertices {
            if edge_ends.contains(&v.0) || loop_vertices.contains(&v.0) {
                add(&mut issues, K::AcornVertexUsed, En::Vertex(v.0));
            }
        }
        if shell.sides.is_empty()
            && !(shell.wire_edges.is_empty() && shell.acorn_vertices.is_empty())
        {
            let parts: Vec<BTreeSet<usize>> = shell
                .wire_edges
                .iter()
                .map(|e| {
                    [edges[e.0].start, edges[e.0].end]
                        .into_iter()
                        .flatten()
                        .map(|v| v.0)
                        .collect()
                })
                .chain(shell.acorn_vertices.iter().map(|v| BTreeSet::from([v.0])))
                .collect();
            let mut seen = vec![false; parts.len()];
            let mut stack = vec![0];
            seen[0] = true;
            while let Some(p) = stack.pop() {
                for q in 0..parts.len() {
                    if !seen[q] && !parts[p].is_disjoint(&parts[q]) {
                        seen[q] = true;
                        stack.push(q);
                    }
                }
            }
            if seen.iter().any(|x| !x) {
                add(&mut issues, K::DisconnectedShell, En::Shell(si));
            }
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
    // C1 in each cell's own parameterisation (R4): exact, never uncertified.
    for (i, edge) in edges.iter().enumerate() {
        if let Curve3::BSpline(spline) = &edge.curve {
            if !continuity::curve_c1(spline.curve(), spline.range(), spline.is_closed_period()) {
                add(&mut issues, K::EdgeNotC1, En::Edge(i));
            }
        }
    }
    for (fi, face) in faces.iter().enumerate() {
        if let Surface::BSpline(spline) = &face.surface {
            if !continuity::surface_c1(spline) {
                add(&mut issues, K::FaceNotC1, En::Face(fi));
                bad_faces.insert(fi);
            }
        }
        for (li, l) in face.loops.iter().enumerate() {
            let Loop::Edges { fins: list, .. } = &loops[l.0] else {
                continue;
            };
            for (ui, k) in list.iter().enumerate() {
                if let Curve2::BSpline(spline) = &fins[k.0].pcurve {
                    if !continuity::curve_c1(
                        spline.curve().as_curve3(),
                        spline.range(),
                        spline.is_closed_period(),
                    ) {
                        add(&mut issues, K::PcurveNotC1, En::Use(fi, li, ui));
                        bad_faces.insert(fi);
                    }
                }
            }
        }
    }
    // Enclosures (M5): a usable one lies in [0, resolution]; the geometric
    // checks below decide against it first, then against the resolution.
    let usable = |issues: &mut BTreeSet<Issue>, e: &Option<Enclosure>, entity| match e {
        None => {
            add(issues, K::EnclosureMissing, entity);
            None
        }
        Some(e) if !(e.bound >= 0.0 && e.bound <= tolerance.linear()) => {
            add(issues, K::EnclosureExceedsResolution, entity);
            None
        }
        Some(e) => Some(e.bound),
    };
    let vertex_bound: Vec<Option<f64>> = vertices
        .iter()
        .enumerate()
        .map(|(v, x)| usable(&mut issues, &x.enclosure, En::Vertex(v)))
        .collect();
    let face_bound: Vec<Option<f64>> = faces
        .iter()
        .enumerate()
        .map(|(fi, f)| usable(&mut issues, &f.enclosure, En::Face(fi)))
        .collect();
    let mut fin_bound: Vec<Option<f64>> = vec![None; nfin];
    let mut fin_seen = vec![false; nfin];
    for (fi, face) in faces.iter().enumerate() {
        for (li, l) in face.loops.iter().enumerate() {
            if let Loop::Edges { fins: list, .. } = &loops[l.0] {
                for (ui, k) in list.iter().enumerate() {
                    if !std::mem::replace(&mut fin_seen[k.0], true) {
                        fin_bound[k.0] =
                            usable(&mut issues, &fins[k.0].enclosure, En::Use(fi, li, ui));
                    }
                }
            }
        }
    }
    let enclosure_verdict =
        |issues: &mut BTreeSet<Issue>, bound: Option<Verdict>, entity| match bound {
            Some(Verdict::Beyond) => add(issues, K::EnclosureUnsound, entity),
            Some(Verdict::Unknown) => add(issues, K::UncertifiedEnclosure, entity),
            _ => {}
        };
    for (i, edge) in edges.iter().enumerate() {
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            continue;
        };
        if !curve_ok[i] {
            continue;
        }
        for (which, v, t) in [(EdgeEnd::Start, start.0, 0.0), (EdgeEnd::End, end.0, 1.0)] {
            if !vertex_ok[v] {
                continue;
            }
            let at = vertices[v].position.to_array();
            let decide = |th: &Threshold| {
                tiered(
                    Verdict::Unknown,
                    || vertex_gap::<Fast>(&edge.curve, at, t, &th.tier::<Fast>().1),
                    || vertex_gap::<I>(&edge.curve, at, t, &th.tier::<I>().1),
                )
            };
            let (verdict, bound) = bounded(vertex_bound[v], tolerance.linear(), decide);
            enclosure_verdict(&mut issues, bound, En::Vertex(v));
            match verdict {
                Verdict::Within => {}
                Verdict::Beyond => add(&mut issues, K::VertexOffCurve, En::EdgeEnd(i, which)),
                Verdict::Unknown => add(
                    &mut issues,
                    K::UncertifiedVertexOffCurve,
                    En::EdgeEnd(i, which),
                ),
            }
        }
    }
    for (fi, face) in faces.iter().enumerate() {
        for (li, l) in face.loops.iter().enumerate() {
            let list = match &loops[l.0] {
                Loop::Vertex(v) => {
                    if surface_ok[fi] && vertex_ok[v.0] {
                        let at = vertices[v.0].position.to_array();
                        let decide = |th: &Threshold| {
                            tiered(
                                Verdict::Unknown,
                                || on_surface::<Fast>(&face.surface, at, &th.tier::<Fast>().1),
                                || on_surface::<I>(&face.surface, at, &th.tier::<I>().1),
                            )
                        };
                        let (verdict, bound) =
                            bounded(vertex_bound[v.0], tolerance.linear(), decide);
                        enclosure_verdict(&mut issues, bound, En::Vertex(v.0));
                        match verdict {
                            Verdict::Within if pole_position(face, loops) == Some(li) => {
                                // A pole must also sit at the apex or pole.
                                let north = pole_north(face_turns(face, loops), face.sense);
                                let decide = |th: &Threshold| {
                                    let apex = |t2: Option<Verdict>| t2.unwrap_or(Verdict::Unknown);
                                    tiered(
                                        Verdict::Unknown,
                                        || {
                                            apex(
                                                pole_gap2::<Fast>(&face.surface, north, at)
                                                    .map(|d| within(&d, &th.tier::<Fast>().1)),
                                            )
                                        },
                                        || {
                                            apex(
                                                pole_gap2::<I>(&face.surface, north, at)
                                                    .map(|d| within(&d, &th.tier::<I>().1)),
                                            )
                                        },
                                    )
                                };
                                let (verdict, bound) =
                                    bounded(vertex_bound[v.0], tolerance.linear(), decide);
                                enclosure_verdict(&mut issues, bound, En::Vertex(v.0));
                                match verdict {
                                    Verdict::Within => {}
                                    Verdict::Beyond => {
                                        add(&mut issues, K::PoleOffApex, En::Loop(fi, li));
                                        bad_faces.insert(fi);
                                    }
                                    Verdict::Unknown => {
                                        add(
                                            &mut issues,
                                            K::UncertifiedVertexLoop,
                                            En::Loop(fi, li),
                                        );
                                        bad_faces.insert(fi);
                                    }
                                }
                            }
                            Verdict::Within => {}
                            Verdict::Beyond => {
                                add(&mut issues, K::VertexLoopOffSurface, En::Loop(fi, li));
                                bad_faces.insert(fi);
                            }
                            Verdict::Unknown => {
                                add(&mut issues, K::UncertifiedVertexLoop, En::Loop(fi, li));
                                bad_faces.insert(fi);
                            }
                        }
                    }
                    continue;
                }
                Loop::Edges { fins: list, .. } => list,
            };
            for (ui, k) in list.iter().enumerate() {
                let u = &fins[k.0];
                if !pcurve_valid(&u.pcurve) {
                    add(&mut issues, K::DegeneratePcurve, En::Use(fi, li, ui));
                    bad_faces.insert(fi);
                    continue;
                }
                if !surface_ok[fi] || !curve_ok[u.edge.0] {
                    continue;
                }
                let forward = u.sense == Orientation::Forward;
                let curve = &edges[u.edge.0].curve;
                // Spline uses whose factors are all rational are decided by
                // exact composition (S4b); the others are not yet.
                // Spline uses are decided by exact composition where every
                // factor is rational, by Taylor enclosures otherwise (S4b).
                let spline = if spline_use(curve, &face.surface, &u.pcurve) {
                    match spline_deviation::rational_use(curve, &face.surface, &u.pcurve, forward) {
                        Some(exact) => Some(SplineUse::Exact(Box::new(exact))),
                        None => match spline_taylor::taylor_use(
                            curve,
                            &face.surface,
                            &u.pcurve,
                            forward,
                        ) {
                            Some(taylor) => Some(SplineUse::Taylor(Box::new(taylor))),
                            None => {
                                add(
                                    &mut issues,
                                    K::UncertifiedPcurveOffEdge,
                                    En::Use(fi, li, ui),
                                );
                                bad_faces.insert(fi);
                                continue;
                            }
                        },
                    }
                } else {
                    None
                };
                let decide = |th: &Threshold| match &spline {
                    Some(SplineUse::Exact(exact)) => match exact.decide(th.0) {
                        Verdict::Unknown => {
                            spline_taylor::taylor_use(curve, &face.surface, &u.pcurve, forward)
                                .map_or(Verdict::Unknown, |taylor| taylor.decide(th.0))
                        }
                        verdict => verdict,
                    },
                    Some(SplineUse::Taylor(taylor)) => taylor.decide(th.0),
                    None => tiered(
                        Verdict::Unknown,
                        || {
                            let (t, t2) = th.tier::<Fast>();
                            deviation::<Fast>(curve, &face.surface, &u.pcurve, forward, &t, &t2)
                        },
                        || {
                            let (t, t2) = th.tier::<I>();
                            deviation::<I>(curve, &face.surface, &u.pcurve, forward, &t, &t2)
                        },
                    ),
                };
                let (verdict, bound) = bounded(fin_bound[k.0], tolerance.linear(), decide);
                enclosure_verdict(&mut issues, bound, En::Use(fi, li, ui));
                match verdict {
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
    // UV continuity; the last fin closes on the first shifted by the winding.
    for (fi, face) in faces.iter().enumerate() {
        if !surface_ok[fi] {
            continue;
        }
        let periodic = face.surface.is_periodic();
        for (li, l) in face.loops.iter().enumerate() {
            let Loop::Edges {
                fins: list,
                winding,
            } = &loops[l.0]
            else {
                continue;
            };
            if list.iter().any(|k| !pcurve_valid(&fins[k.0].pcurve)) {
                continue;
            }
            for (ui, k) in list.iter().enumerate() {
                let u = &fins[k.0];
                let w = &fins[list[(ui + 1) % list.len()].0];
                let shift = if ui + 1 == list.len() && periodic {
                    [TAU * f64::from(winding[0]), TAU * f64::from(winding[1])]
                } else {
                    [0.0, 0.0]
                };
                let decide = |th: &Threshold| {
                    tiered(
                        Verdict::Unknown,
                        || {
                            let t2 = th.tier::<Fast>().1;
                            uv_gap::<Fast>(&face.surface, &u.pcurve, &w.pcurve, shift, &t2)
                        },
                        || {
                            let t2 = th.tier::<I>().1;
                            uv_gap::<I>(&face.surface, &u.pcurve, &w.pcurve, shift, &t2)
                        },
                    )
                };
                let (verdict, bound) = bounded(face_bound[fi], tolerance.linear(), decide);
                enclosure_verdict(&mut issues, bound, En::Face(fi));
                match verdict {
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
        let forward = face.sense == Orientation::Forward;
        let want_outer = if forward {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        let want_inner = want_outer.reverse();
        let edge_loops: Vec<(usize, &Lp)> = face
            .loops
            .iter()
            .enumerate()
            .filter(|(_, l)| matches!(loops[l.0], Loop::Edges { .. }))
            .map(|(li, _)| li)
            .zip(resolved[fi].iter())
            .collect();
        let wound = face.surface.is_periodic() && edge_loops.iter().any(|(_, lp)| lp.winding != 0);
        let wound_v =
            face.surface.is_periodic_v() && edge_loops.iter().any(|(_, lp)| lp.winding_v != 0);
        if wound && wound_v {
            // A face wound in both directions (a torus knot): not decided.
            add(&mut issues, K::UncertifiedLoopWinding, En::Loop(fi, 0));
            winding_faces.insert(fi);
            continue;
        }
        if wound_v {
            // The periodic area in v, the integral of u dv, has the face's
            // sign over all loops; each unwound loop the opposite one.
            fn total_v<T: Real>(loops: &[(usize, &Lp)]) -> Option<T> {
                loops.iter().try_fold(c::<T>(0.0), |acc, (_, lp)| {
                    Some(acc.add(&periodic_area_v::<T>(lp)?))
                })
            }
            let total = tiered_integral(
                None,
                || total_v::<Fast>(&edge_loops)?.sign(),
                || total_v::<I>(&edge_loops)?.sign(),
            );
            match total {
                Some(s) if s == want_outer => {}
                Some(_) => {
                    add(&mut issues, K::LoopWinding, En::Loop(fi, 0));
                    winding_faces.insert(fi);
                }
                None => {
                    add(&mut issues, K::UncertifiedLoopWinding, En::Loop(fi, 0));
                    winding_faces.insert(fi);
                }
            }
            for (li, lp) in &edge_loops {
                if lp.winding_v != 0 {
                    continue;
                }
                let sign = tiered_integral(
                    None,
                    || periodic_area_v::<Fast>(lp)?.sign(),
                    || periodic_area_v::<I>(lp)?.sign(),
                );
                match sign {
                    Some(s) if s == want_inner => {
                        // Inside when the signed crossings of the +u ray
                        // from a point of it with the other loops give the
                        // face's sign, the band's right side running in +v
                        // (S9e.3b: the +v ray's rule of a band wound in u
                        // with the parameters exchanged, which reverses the
                        // sign).
                        let others: Vec<&Lp> = edge_loops
                            .iter()
                            .filter(|(l, _)| l != li)
                            .map(|(_, other)| *other)
                            .collect();
                        let want: i64 = if forward { -1 } else { 1 };
                        let start = &lp.fins[0].pcurve;
                        let decided = [0.0, 0.376_953_125, 0.678_710_937_5].iter().find_map(|&t| {
                            tiered(
                                None,
                                || {
                                    Some(
                                        cover_hits_v::<Fast>(&others, &pcurve_at(start, t))?
                                            .iter()
                                            .sum::<i64>(),
                                    )
                                },
                                || {
                                    Some(
                                        cover_hits_v::<I>(&others, &pcurve_at(start, t))?
                                            .iter()
                                            .sum::<i64>(),
                                    )
                                },
                            )
                        });
                        match decided {
                            Some(n) if n == want => {}
                            Some(_) => add(&mut issues, K::InnerLoopOutside, En::Loop(fi, *li)),
                            None => {
                                add(&mut issues, K::UncertifiedContainment, En::Loop(fi, *li));
                            }
                        }
                    }
                    Some(_) => {
                        add(&mut issues, K::LoopWinding, En::Loop(fi, *li));
                        winding_faces.insert(fi);
                    }
                    None => {
                        add(&mut issues, K::UncertifiedLoopWinding, En::Loop(fi, *li));
                        winding_faces.insert(fi);
                    }
                }
            }
            continue;
        }
        if wound {
            // A pole is the line v = v_apex traversed against the band: it
            // adds 2 pi W v_apex, W the edge loops' total winding.
            let pole = pole_position(face, loops).is_some();
            let turns: i32 = edge_loops.iter().map(|(_, lp)| lp.winding).sum();
            let north = pole_north(turns, face.sense);
            fn pole_term<T: Real>(s: &Surface, turns: i32, north: bool) -> Option<T> {
                let two_pi =
                    T::from_r(&(pi().midpoint() * int(2))).widen(&(pi().radius() * int(2)));
                Some(
                    two_pi
                        .mul(&c(f64::from(turns)))
                        .mul(&pole_v::<T>(s, north)?),
                )
            }
            fn total<T: Real>(
                s: &Surface,
                loops: &[(usize, &Lp)],
                pole: bool,
                turns: i32,
                north: bool,
            ) -> Option<T> {
                let mut sum = loops.iter().try_fold(c::<T>(0.0), |acc, (_, lp)| {
                    Some(acc.add(&periodic_area::<T>(lp)?))
                })?;
                if pole {
                    sum = sum.add(&pole_term::<T>(s, turns, north)?);
                }
                Some(sum)
            }
            let total = tiered_integral(
                None,
                || total::<Fast>(&face.surface, &edge_loops, pole, turns, north)?.sign(),
                || total::<I>(&face.surface, &edge_loops, pole, turns, north)?.sign(),
            );
            match total {
                Some(s) if s == want_outer => {}
                Some(_) => {
                    add(&mut issues, K::LoopWinding, En::Loop(fi, 0));
                    winding_faces.insert(fi);
                }
                None => {
                    add(&mut issues, K::UncertifiedLoopWinding, En::Loop(fi, 0));
                    winding_faces.insert(fi);
                }
            }
            for (li, lp) in &edge_loops {
                if lp.winding != 0 {
                    continue;
                }
                let sign = tiered_integral(
                    None,
                    || periodic_area::<Fast>(lp)?.sign(),
                    || periodic_area::<I>(lp)?.sign(),
                );
                match sign {
                    Some(s) if s == want_inner => {
                        // Inside when the signed crossings of the +v ray from
                        // its first point with the other loops, a north pole
                        // counting as one more above, give the face's sign
                        // (S8d.2: a hole in a band).
                        let others: Vec<&Lp> = edge_loops
                            .iter()
                            .filter(|(l, _)| l != li)
                            .map(|(_, other)| *other)
                            .collect();
                        let want: i64 = if forward { 1 } else { -1 };
                        let polar = if pole && north { want } else { 0 };
                        // Its start, else two points along its first fin
                        // (a start on a ray through another loop's vertex is
                        // undecided; S9d.2).
                        let start = &lp.fins[0].pcurve;
                        let decided = [0.0, 0.376_953_125, 0.678_710_937_5].iter().find_map(|&t| {
                            tiered(
                                None,
                                || signed_cover_crossings::<Fast>(&others, &pcurve_at(start, t)),
                                || signed_cover_crossings::<I>(&others, &pcurve_at(start, t)),
                            )
                        });
                        match decided {
                            Some(n) if n + polar == want => {}
                            Some(_) => add(&mut issues, K::InnerLoopOutside, En::Loop(fi, *li)),
                            None => {
                                add(&mut issues, K::UncertifiedContainment, En::Loop(fi, *li));
                            }
                        }
                    }
                    Some(_) => {
                        add(&mut issues, K::LoopWinding, En::Loop(fi, *li));
                        winding_faces.insert(fi);
                    }
                    None => {
                        add(&mut issues, K::UncertifiedLoopWinding, En::Loop(fi, *li));
                        winding_faces.insert(fi);
                    }
                }
            }
            continue;
        }
        // A torus or sphere face whose first loop runs as a hole is the
        // surface less its loops (S8d.3, S9d.1): every loop then has the inner
        // sign, and no outer loop holds them.
        let complement = matches!(face.surface, Surface::Torus { .. } | Surface::Sphere { .. })
            && edge_loops.first().is_some_and(|(_, lp)| {
                tiered_integral(
                    None,
                    || loop_area::<Fast>(lp).sign(),
                    || loop_area::<I>(lp).sign(),
                ) == Some(want_inner)
            });
        for (pos, (li, lp)) in edge_loops.iter().enumerate() {
            let sign = tiered_integral(
                None,
                || loop_area::<Fast>(lp).sign(),
                || loop_area::<I>(lp).sign(),
            );
            let want = if pos == 0 && !complement {
                want_outer
            } else {
                want_inner
            };
            match sign {
                Some(s) if s == want => {}
                Some(_) => {
                    add(&mut issues, K::LoopWinding, En::Loop(fi, *li));
                    winding_faces.insert(fi);
                }
                None => {
                    add(&mut issues, K::UncertifiedLoopWinding, En::Loop(fi, *li));
                    winding_faces.insert(fi);
                }
            }
        }
        if let Some((_, outer)) = edge_loops.first().filter(|_| !complement) {
            for (li, lp) in edge_loops.iter().skip(1) {
                let start = &lp.fins[0].pcurve;
                let outer = [*outer];
                // Any point of the inner loop decides: its start, else two
                // points along its first fin (a start on the ray through
                // another loop's start, as two rings' may be, is undecided).
                let decide = || {
                    [0.0, 0.376_953_125, 0.678_710_937_5].iter().find_map(|&t| {
                        tiered(
                            None,
                            || crossings::<Fast>(&outer, &pcurve_at(start, t)),
                            || crossings::<I>(&outer, &pcurve_at(start, t)),
                        )
                    })
                };
                match decide() {
                    Some(n) if n % 2 == 1 => {}
                    Some(_) => add(&mut issues, K::InnerLoopOutside, En::Loop(fi, *li)),
                    None => add(&mut issues, K::UncertifiedContainment, En::Loop(fi, *li)),
                }
            }
        }
    }

    // Region orientation and cavity nesting on fully sound shells.
    bad_faces.extend(structural_faces);
    bad_faces.extend(winding_faces);
    bad_faces.extend(side_bad);
    let sound: BTreeSet<usize> = (0..ns)
        .filter(|si| {
            !bad_shells.contains(si)
                && !twins.contains(si)
                && !shells[*si].sides.is_empty()
                && shells[*si]
                    .sides
                    .iter()
                    .all(|(f, _)| !bad_faces.contains(&f.0))
                // A sheet shell has no orientation to decide (S6).
                && shells[*si].sides.len() > 2 * two_sided(&shells[*si]).len()
        })
        .collect();
    for (ri, region) in regions.iter().enumerate() {
        let mut oriented = Vec::new();
        for (pos, s) in region.shells.iter().enumerate() {
            let si = s.0;
            if !sound.contains(&si) {
                continue;
            }
            let origin = shell_origin(view, &resolved, &shells[si]);
            let sign = tiered_integral(
                None,
                || shell_flux::<Fast>(faces, &resolved, &shells[si], origin)?.sign(),
                || shell_flux::<I>(faces, &resolved, &shells[si], origin)?.sign(),
            );
            let want = if ri > 0 && pos == 0 {
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
        let Some(outer) = region.shells.first().map(|s| s.0) else {
            continue;
        };
        if ri == 0 || !oriented.contains(&outer) {
            continue;
        }
        let cavities: Vec<usize> = oriented.iter().copied().filter(|s| *s != outer).collect();
        for &si in &cavities {
            let Some(point) = shell_point(view, &resolved, si, &vertex_ok) else {
                continue;
            };
            let contained = |shell: usize| {
                tiered(
                    None,
                    || inside::<Fast>(faces, &resolved, &shells[shell], &point, &fast_margin2),
                    || inside::<I>(faces, &resolved, &shells[shell], &point, &exact_margin2),
                )
            };
            match contained(outer) {
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
                match contained(sj) {
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

/// A point of a shell: its first face's first fin's start vertex, or the
/// surface point at that fin's pcurve start when it has no vertex.
fn shell_point(view: &View, resolved: &[Vec<Lp>], si: usize, vertex_ok: &[bool]) -> Option<[R; 3]> {
    let f = view.shells[si].sides.first()?.0 .0;
    let face = &view.faces[f];
    if let Some(Loop::Vertex(v)) = face.loops.first().map(|l| &view.loops[l.0]) {
        return Some(view.vertices[v.0].position.to_array().map(r));
    }
    // A face without loops (a whole sphere): its point at (0, 0).
    if face.loops.is_empty() {
        let p = face.surface.point(crate::Point2::default());
        return finite(&p.to_array()).then(|| p.to_array().map(r));
    }
    let fin = resolved[f].first()?.fins.first()?;
    match fin_vertices(view.edges, fin).0 {
        Some(v) if vertex_ok[v] => Some(view.vertices[v].position.to_array().map(r)),
        Some(_) => None,
        None => {
            let p = face.surface.point(fin.pcurve.point(0.0));
            finite(&p.to_array()).then(|| p.to_array().map(r))
        }
    }
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

    fn fin(a: (f64, f64), b: (f64, f64)) -> Fin {
        Fin {
            edge: EdgeId::new(0),
            sense: Orientation::Forward,
            pcurve: Curve2::LineSegment {
                start: Point2::new(a.0, a.1),
                end: Point2::new(b.0, b.1),
            },
            enclosure: None,
        }
    }

    fn clear_in(
        fins: &[Fin],
        winding: i32,
        p: (f64, f64),
        su: f64,
        periodic: bool,
        margin: f64,
    ) -> bool {
        let lp = Lp {
            fins: fins.iter().collect(),
            winding,
            winding_v: 0,
        };
        let fast = clear_of_boundary::<Fast>(
            &[&lp],
            &[c(p.0), c(p.1)],
            &c(su),
            periodic,
            &c::<Fast>(margin).square(),
        );
        let exact = clear_of_boundary::<I>(
            &[&lp],
            &[c(p.0), c(p.1)],
            &c(su),
            periodic,
            &c::<I>(margin).square(),
        );
        assert_eq!(fast, exact, "{p:?}");
        exact.unwrap()
    }

    fn clear(u: f64, v: f64, su: f64, margin: f64) -> bool {
        let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let square: Vec<Fin> = (0..4)
            .map(|k| fin(corners[k], corners[(k + 1) % 4]))
            .collect();
        clear_in(&square, 0, (u, v), su, false, margin)
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

    #[test]
    fn periodic_margins_see_every_alias_of_a_ring_loop() {
        // A ring loop at v = 0 winding once; a point just below it on the far
        // side of the period is near an alias, one above the strip is clear.
        let ring = [fin((0.0, 0.0), (TAU, 0.0))];
        assert!(!clear_in(&ring, 1, (TAU + 0.5, 1e-9), 1.0, true, 1e-8));
        assert!(!clear_in(&ring, 1, (-3.0, -1e-9), 1.0, true, 1e-8));
        assert!(clear_in(&ring, 1, (-3.0, 0.5), 1.0, true, 1e-8));
    }
}
