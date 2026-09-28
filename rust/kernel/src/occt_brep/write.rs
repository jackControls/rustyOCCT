//! The writer: a cell topology's solid regions as `.brep` version 1 text,
//! with OCCT's structure inserted by rule. A face that winds around a
//! cylinder or cone gets one seam at the `u` where its ring loops start: its
//! forward use (in the unoriented face) at `u0 + 2 pi` going from the lower
//! ring to the upper one, its reversed use at `u0`; every ring edge gets a
//! vertex where the seam meets it and becomes a closed edge. A cone's pole
//! becomes a degenerated edge at the apex, as `BRepPrim_Cone` builds it: no
//! 3D curve, the pcurve `v = v_apex` over one turn from `u0`, the apex its
//! vertex at both ends, and the seam ends there. Every pcurve shares its
//! edge's parameter, as OCCT's SameParameter edges do. A body without a
//! solid region (S6) is written as its closed shells, its face or open
//! shell, its edge or wire, or its vertex.
use super::BrepError;
use crate::topology::{
    BodyClass, Curve2, Curve3, EdgeId, FaceId, Loop, Orientation, RegionKind, Side, Surface,
    Topology,
};
use crate::{BSplineCurve3, BSplineSurface3, Point2, Point3};
use std::collections::BTreeMap;
use std::f64::consts::TAU;
use std::fmt::Write as _;

fn num(x: f64) -> String {
    // Shortest round-trip text; OCCT reads it with operator>>.
    let s = format!("{x:?}");
    s.strip_suffix(".0").map(str::to_string).unwrap_or(s)
}

fn nums(v: &[f64]) -> String {
    v.iter().map(|x| num(*x)).collect::<Vec<_>>().join(" ")
}

fn unwritable(what: &'static str) -> BrepError {
    BrepError::Unwritable { what }
}

/// A 3D curve record with the edge's parameter range.
struct EdgeGeometry {
    record: String,
    range: [f64; 2],
}

/// The 3D record of a curve and its range, the parameter increasing along
/// the edge.
fn curve_record(c: &Curve3, ring_start: Option<f64>) -> EdgeGeometry {
    match c {
        Curve3::BSpline(span) => EdgeGeometry {
            record: bspline_record(span.curve(), |p| nums(&p.to_array())),
            range: span.range(),
        },
        Curve3::LineSegment { start, end } => {
            let d = *end - *start;
            let l = d.length();
            let u = d * (1.0 / l);
            EdgeGeometry {
                record: format!("1 {} {}", nums(&start.to_array()), nums(&u.to_array())),
                range: [0.0, l],
            }
        }
        Curve3::Circle { frame, radius } => curve_record(
            &Curve3::CircularArc {
                frame: *frame,
                radius: *radius,
                start_angle: 0.0,
                sweep_angle: TAU,
            },
            ring_start,
        ),
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => {
            // A negative sweep is the positive sweep about the flipped axis:
            // t = -angle.
            let flip = *sweep_angle < 0.0;
            let (n, y) = if flip {
                (-frame.normal(), -frame.normal().cross(frame.x()))
            } else {
                (frame.normal(), frame.normal().cross(frame.x()))
            };
            let start = match ring_start {
                Some(s) => s,
                None if flip => -start_angle,
                None => *start_angle,
            };
            EdgeGeometry {
                record: format!(
                    "2 {} {} {} {} {}",
                    nums(&frame.origin().to_array()),
                    nums(&n.to_array()),
                    nums(&frame.x().to_array()),
                    nums(&y.to_array()),
                    num(*radius)
                ),
                range: [start, start + sweep_angle.abs()],
            }
        }
        // `Geom_Hyperbola` (record 5) and `Geom_Parabola` (record 6), their
        // parameters the edge's; a negative sweep about the flipped axis
        // (both are even in `x` and odd in `y` under `t -> -t`).
        Curve3::HyperbolaArc {
            frame,
            major,
            minor,
            start,
            sweep,
        } => {
            let flip = *sweep < 0.0;
            let (n, y) = if flip {
                (-frame.normal(), -frame.y())
            } else {
                (frame.normal(), frame.y())
            };
            let s = if flip { -start } else { *start };
            EdgeGeometry {
                record: format!(
                    "5 {} {} {} {} {} {}",
                    nums(&frame.origin().to_array()),
                    nums(&n.to_array()),
                    nums(&frame.x().to_array()),
                    nums(&y.to_array()),
                    num(*major),
                    num(*minor)
                ),
                range: [s, s + sweep.abs()],
            }
        }
        Curve3::ParabolaArc {
            frame,
            focal,
            start,
            sweep,
        } => {
            let flip = *sweep < 0.0;
            let (n, y) = if flip {
                (-frame.normal(), -frame.y())
            } else {
                (frame.normal(), frame.y())
            };
            let s = if flip { -start } else { *start };
            EdgeGeometry {
                record: format!(
                    "6 {} {} {} {} {}",
                    nums(&frame.origin().to_array()),
                    nums(&n.to_array()),
                    nums(&frame.x().to_array()),
                    nums(&y.to_array()),
                    num(*focal)
                ),
                range: [s, s + sweep.abs()],
            }
        }
        // `Geom_Ellipse` (record 3), its major radius first: an ellipse
        // longer along y is written about its axes turned a quarter turn,
        // `t - pi/2` its angle; a negative sweep about the flipped axis.
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => {
            let (mut x, mut y) = (frame.x(), frame.y());
            let (mut big, mut small, mut shift) = (*major, *minor, 0.0);
            if major < minor {
                (x, y) = (frame.y(), -frame.x());
                (big, small, shift) = (*minor, *major, -std::f64::consts::FRAC_PI_2);
            }
            let flip = *sweep_angle < 0.0;
            let (n, y) = if flip {
                (-frame.normal(), -y)
            } else {
                (frame.normal(), y)
            };
            let angle = start_angle + shift;
            let start = match ring_start {
                Some(s) => s,
                None if flip => -angle,
                None => angle,
            };
            EdgeGeometry {
                record: format!(
                    "3 {} {} {} {} {} {}",
                    nums(&frame.origin().to_array()),
                    nums(&n.to_array()),
                    nums(&x.to_array()),
                    nums(&y.to_array()),
                    num(big),
                    num(small)
                ),
                range: [start, start + sweep_angle.abs()],
            }
        }
    }
}

/// A B-spline curve record (7), 3D or 2D: flags, degree, counts, poles
/// (each with its weight when rational), knots with multiplicities, as
/// `GeomTools_CurveSet` prints it.
fn bspline_record(c: &BSplineCurve3, pole: impl Fn(&Point3) -> String) -> String {
    let rational = c.is_rational();
    let mut text = format!(
        "7 {} {} {} {} {}",
        u8::from(rational),
        u8::from(c.is_periodic()),
        c.degree(),
        c.poles().len(),
        c.knots().len()
    );
    for (p, w) in c.poles().iter().zip(c.weights()) {
        text.push_str("\n ");
        text.push_str(&pole(p));
        if rational {
            text.push(' ');
            text.push_str(&num(*w));
        }
    }
    for (k, m) in c.knots().iter().zip(c.multiplicities()) {
        let _ = write!(text, "\n {} {m}", num(*k));
    }
    text
}

/// A B-spline surface record (9): flags, degrees, counts, poles U-major
/// (each with its weight when rational), then the u and v knots.
fn bspline_surface_record(s: &BSplineSurface3) -> String {
    let (u, v) = (s.u_knots(), s.v_knots());
    let rational = s.is_rational();
    let mut text = format!(
        "9 {r} {r} {} {} {} {} {} {} {} {}",
        u8::from(u.is_periodic()),
        u8::from(v.is_periodic()),
        u.degree(),
        v.degree(),
        u.pole_count(),
        v.pole_count(),
        u.knots().len(),
        v.knots().len(),
        r = u8::from(rational)
    );
    for (p, w) in s.poles().iter().zip(s.weights()) {
        text.push_str("\n ");
        text.push_str(&nums(&p.to_array()));
        if rational {
            text.push(' ');
            text.push_str(&num(*w));
        }
    }
    for axis in [u, v] {
        for (k, m) in axis.knots().iter().zip(axis.multiplicities()) {
            let _ = write!(text, "\n {} {m}", num(*k));
        }
    }
    text
}

/// The point of a curve at its own parameter `t` (as `curve_record` writes it).
fn curve_at(c: &Curve3, t: f64) -> Point3 {
    match c {
        Curve3::BSpline(span) => span.curve().point(t).expect("a parameter in the range"),
        Curve3::LineSegment { start, end } => {
            let d = *end - *start;
            *start + d * (t / d.length())
        }
        Curve3::Circle { frame, radius } => {
            frame.point(Point2::new(radius * t.cos(), radius * t.sin()), 0.0)
        }
        Curve3::CircularArc {
            frame,
            radius,
            sweep_angle,
            ..
        } => {
            let a = if *sweep_angle < 0.0 { -t } else { t };
            frame.point(Point2::new(radius * a.cos(), radius * a.sin()), 0.0)
        }
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            sweep_angle,
            ..
        } => {
            let a = if *sweep_angle < 0.0 { -t } else { t };
            let a = if major < minor {
                a + std::f64::consts::FRAC_PI_2
            } else {
                a
            };
            frame.point(Point2::new(major * a.cos(), minor * a.sin()), 0.0)
        }
        Curve3::HyperbolaArc {
            frame,
            major,
            minor,
            sweep,
            ..
        } => {
            let a = if *sweep < 0.0 { -t } else { t };
            frame.point(Point2::new(major * a.cosh(), minor * a.sinh()), 0.0)
        }
        Curve3::ParabolaArc {
            frame,
            focal,
            sweep,
            ..
        } => {
            let a = if *sweep < 0.0 { -t } else { t };
            frame.point(Point2::new(a * a / (4.0 * focal), a), 0.0)
        }
    }
}

/// The 2D record of a fin's pcurve in the edge's parameter `[t0, t1]`.
fn pcurve_record(
    surface: &Surface,
    pcurve: &Curve2,
    forward: bool,
    range: [f64; 2],
) -> Result<String, BrepError> {
    let [t0, t1] = range;
    Ok(match (surface, pcurve) {
        // A spline pcurve is written in the edge's direction and parameter:
        // unflagged there, over the edge's range to printing precision.
        (_, Curve2::BSpline(span)) => {
            let span = if forward {
                span.clone()
            } else {
                span.reversed()
            };
            let Some(span) = span.unflagged() else {
                return Err(unwritable("a spline pcurve whose mirror rounds"));
            };
            let [a, b] = span.range();
            let near = |x: f64, y: f64| (x - y).abs() <= 1e-12 * (1.0 + (t1 - t0).abs() + y.abs());
            if span.is_reversed() || !near(a, t0) || !near(b, t1) {
                return Err(unwritable("a spline pcurve off its edge's parameter"));
            }
            let c = span.curve().as_curve3();
            bspline_record(c, |p| nums(&[p.x, p.y]))
        }
        (Surface::BSpline(_), Curve2::CircularArc { .. }) => {
            return Err(unwritable("an arc pcurve on a spline surface"));
        }
        (_, Curve2::LineSegment { start, end }) => {
            // In the edge's direction.
            let (a, b) = if forward {
                (*start, *end)
            } else {
                (*end, *start)
            };
            let (du, dv) = (b.x - a.x, b.y - a.y);
            let span = t1 - t0;
            let (dx, dy) = match surface {
                // On a cylinder, cone or sphere a horizontal pcurve runs in u
                // at unit speed against an angle parameter.
                Surface::Cylinder { .. } | Surface::Cone { .. } | Surface::Sphere { .. }
                    if dv == 0.0 =>
                {
                    (du.signum(), 0.0)
                }
                _ => {
                    let l = du.hypot(dv);
                    (du / l, dv / l)
                }
            };
            let _ = span;
            format!(
                "1 {} {}",
                nums(&[a.x - dx * t0, a.y - dy * t0]),
                nums(&[dx, dy])
            )
        }
        (
            Surface::Plane(_),
            Curve2::CircularArc {
                center,
                radius,
                start_angle,
                sweep_angle,
            },
        ) => {
            let (start, sweep) = if forward {
                (*start_angle, *sweep_angle)
            } else {
                (start_angle + sweep_angle, -sweep_angle)
            };
            let sigma = sweep.signum();
            let alpha = start - sigma * t0;
            let (x, y) = (
                (alpha.cos(), alpha.sin()),
                (-sigma * alpha.sin(), sigma * alpha.cos()),
            );
            format!(
                "2 {} {} {} {}",
                nums(&[center.x, center.y]),
                nums(&[x.0, x.1]),
                nums(&[y.0, y.1]),
                num(*radius)
            )
        }
        (
            Surface::Cylinder { .. }
            | Surface::Cone { .. }
            | Surface::Sphere { .. }
            | Surface::Torus { .. },
            Curve2::CircularArc { .. },
        ) => return Err(unwritable("an arc pcurve on a periodic surface")),
        // `Geom2d_Ellipse` (record 3) on its own axes, the edge's parameter
        // its angle: y flipped when the angle decreases along the edge.
        (
            Surface::Plane(_),
            Curve2::EllipseArc {
                center,
                major,
                minor,
                start_angle,
                sweep_angle,
            },
        ) => {
            let (start, sweep) = if forward {
                (*start_angle, *sweep_angle)
            } else {
                (start_angle + sweep_angle, -sweep_angle)
            };
            let sigma = sweep.signum();
            let alpha = start - sigma * t0;
            if major < minor || alpha.abs() > 1e-9 * (1.0 + t0.abs()) {
                return Err(unwritable(
                    "an ellipse pcurve whose angle is not its edge's parameter",
                ));
            }
            format!(
                "3 {} 1 0 0 {} {} {}",
                nums(&[center.x, center.y]),
                num(sigma),
                num(*major),
                num(*minor)
            )
        }
        (_, Curve2::EllipseArc { .. }) => return Err(unwritable("an ellipse pcurve off a plane")),
        // OCCT has no analytic record for a cylinder's plane section.
        (_, Curve2::Sinusoid { .. }) => {
            return Err(unwritable("a sinusoid pcurve (a cylinder's plane section)"))
        }
        // D13's exact projections are written as approximations only.
        (_, Curve2::Projection(_)) => return Err(unwritable("a projection pcurve (D13)")),
    })
}

/// Where a seam meets a wound loop: an existing vertex, the seam vertex of
/// a ring edge, or a pole vertex the topology does not have (a whole
/// sphere's), by index into the extra vertices.
#[derive(Clone, Copy)]
enum Vertex {
    Shared(usize),
    Ring(usize),
    Pole(usize),
}

#[derive(Clone, Copy)]
struct End {
    vertex: Vertex,
    /// `v` of the seam end.
    height: f64,
}

/// A wound loop: its winding in the unoriented face, its ring edge if it is
/// one, and otherwise each vertex with its position in the face (a pole's
/// vertex at any `u`); whether it is a pole.
type WoundLoop = (i32, Option<usize>, Vec<(usize, Point2)>, bool);

/// Whether consecutive fins meeting at `a` and `b` pass through a pole of
/// the face's cone or sphere (a different `u` on the pole's line), where
/// OCCT has a degenerated edge (as `Topology::pole_passes` counts them).
fn pole_pass(surface: &Surface, a: Point2, b: Point2) -> bool {
    let poles: Vec<f64> = match surface {
        Surface::Sphere { .. } => vec![std::f64::consts::FRAC_PI_2, -std::f64::consts::FRAC_PI_2],
        Surface::Cone {
            radius, half_angle, ..
        } => vec![-radius / half_angle.sin()],
        _ => return false,
    };
    let at_pole = poles
        .iter()
        .any(|p| (a.y - p).abs() <= 1e-12 * p.abs().max(1.0));
    let du = (a.x - b.x).rem_euclid(TAU);
    at_pole && du.min(TAU - du) > 1e-9
}

/// A surface of revolution: its frame and meridian.
#[derive(Clone, Copy)]
enum Revolved {
    /// `rho(v) = radius + v sin`, height `v cos` (a cylinder has `sin = 0`).
    Ruled {
        frame: crate::Frame3,
        radius: f64,
        sin: f64,
        cos: f64,
    },
    /// `rho(v) = radius cos v`, height `radius sin v`.
    Sphere { frame: crate::Frame3, radius: f64 },
    /// `rho(v) = major + minor cos v`, height `minor sin v`.
    Torus {
        frame: crate::Frame3,
        major: f64,
        minor: f64,
    },
}

impl Revolved {
    fn of(surface: &Surface) -> Option<Self> {
        match surface {
            Surface::BSpline(_) => None,
            Surface::Cylinder { frame, radius } => Some(Self::Ruled {
                frame: *frame,
                radius: *radius,
                sin: 0.0,
                cos: 1.0,
            }),
            Surface::Cone {
                frame,
                radius,
                half_angle,
            } => Some(Self::Ruled {
                frame: *frame,
                radius: *radius,
                sin: half_angle.sin(),
                cos: half_angle.cos(),
            }),
            Surface::Sphere { frame, radius } => Some(Self::Sphere {
                frame: *frame,
                radius: *radius,
            }),
            Surface::Torus {
                frame,
                major,
                minor,
            } => Some(Self::Torus {
                frame: *frame,
                major: *major,
                minor: *minor,
            }),
            Surface::Plane(_) => None,
        }
    }
    fn frame(&self) -> crate::Frame3 {
        match self {
            Self::Ruled { frame, .. } | Self::Sphere { frame, .. } | Self::Torus { frame, .. } => {
                *frame
            }
        }
    }
    /// `(rho, height)` at `v`.
    fn meridian(&self, v: f64) -> (f64, f64) {
        match *self {
            Self::Ruled {
                radius, sin, cos, ..
            } => (radius + v * sin, v * cos),
            Self::Sphere { radius, .. } => (radius * v.cos(), radius * v.sin()),
            Self::Torus { major, minor, .. } => (major + minor * v.cos(), minor * v.sin()),
        }
    }
    fn point(&self, u: f64, v: f64) -> Point3 {
        let (rho, height) = self.meridian(v);
        self.frame()
            .point(Point2::new(rho * u.cos(), rho * u.sin()), height)
    }
    /// `v` of the pole a band closes: a cone's apex, a sphere's north pole
    /// when the band winds `+u` in the unoriented face, else its south pole.
    fn pole_v(&self, north: bool) -> Option<f64> {
        match *self {
            Self::Ruled { radius, sin, .. } if sin != 0.0 => Some(-radius / sin),
            Self::Ruled { .. } => None,
            Self::Sphere { .. } if north => Some(std::f64::consts::FRAC_PI_2),
            Self::Sphere { .. } => Some(-std::f64::consts::FRAC_PI_2),
            Self::Torus { .. } => None,
        }
    }
}

/// A wound face's seam at `u0`, from the loop winding `+u` (in the
/// unoriented face) to the loop winding `-u`.
struct Seam {
    u0: f64,
    bottom: End,
    top: End,
}

/// A torus face's seam in v: a latitude arc at `v0` over `range` in u, from
/// the `-v` loop's seam vertex to the `+v` loop's.
struct VSeam {
    v0: f64,
    from: Vertex,
    to: Vertex,
    range: [f64; 2],
}

struct Records {
    /// Shape records in file order.
    shapes: Vec<String>,
}

impl Records {
    fn push(&mut self, record: String) -> usize {
        self.shapes.push(record);
        self.shapes.len() - 1
    }
}

/// Version 1 `.brep` text of every solid region of `topology`; one solid is
/// the root, several a compound.
pub fn write(topology: &Topology, tolerance: f64) -> Result<String, BrepError> {
    let t = topology;
    // A spline edge traversed against its curve would need mirrored knots.
    if t.edges()
        .iter()
        .any(|e| matches!(&e.curve, Curve3::BSpline(s) if s.is_reversed()))
    {
        return Err(unwritable("a spline edge against its curve"));
    }
    let tol = num(tolerance);
    let mut curves2d: Vec<String> = Vec::new();
    let mut curves: Vec<String> = Vec::new();
    let mut surfaces: Vec<String> = Vec::new();
    let mut records = Records { shapes: Vec::new() };
    // Seams: one per wound face, at a u0 where both wound loops have a
    // vertex (a ring loop gets its seam vertex there).
    let mut ring_start: BTreeMap<usize, f64> = BTreeMap::new();
    let mut wound: BTreeMap<usize, Seam> = BTreeMap::new();
    // The poles of each wound face: the vertex, `v` of the pole and whether
    // its degenerated edge runs `+u` in the unoriented face.
    let mut poles: BTreeMap<usize, Vec<(Vertex, f64, bool)>> = BTreeMap::new();
    // Pole points the topology has no vertex for (a whole sphere's).
    let mut extra: Vec<Point3> = Vec::new();
    for (fi, face) in t.faces().iter().enumerate() {
        let Some(rev) = Revolved::of(&face.surface) else {
            continue;
        };
        let frame = rev.frame();
        let rho = |v: f64| rev.meridian(v).0.abs();
        let flip = face.sense == Orientation::Reversed;
        let mut loops: Vec<WoundLoop> = Vec::new();
        let mut pole = None;
        for l in &face.loops {
            let Loop::Edges { fins, winding } = &t.loops()[l.index()] else {
                if let (Loop::Vertex(v), Surface::Cone { .. } | Surface::Sphere { .. }, None) =
                    (&t.loops()[l.index()], &face.surface, pole)
                {
                    pole = Some(v.index());
                }
                continue;
            };
            if winding[1] != 0 {
                // A torus wound in v has its seam in v (below).
                if matches!(rev, Revolved::Torus { .. }) && winding[0] == 0 {
                    continue;
                }
                return Err(unwritable("a winding in v and u"));
            }
            if winding[0] == 0 {
                continue;
            }
            if winding[0].abs() != 1 {
                return Err(unwritable("a loop winding more than once"));
            }
            let first = &t.fins()[fins[0].index()];
            let ring = (fins.len() == 1 && t.edges()[first.edge.index()].is_ring())
                .then_some(first.edge.index());
            let mut at = Vec::new();
            for k in fins {
                let f = &t.fins()[k.index()];
                let edge = &t.edges()[f.edge.index()];
                let v = if f.sense == Orientation::Forward {
                    edge.start
                } else {
                    edge.end
                };
                if let Some(v) = v {
                    at.push((v.index(), f.pcurve.point(0.0)));
                }
            }
            if ring.is_none() && at.len() != fins.len() {
                return Err(unwritable("a wound loop mixing rings and vertices"));
            }
            loops.push((if flip { -winding[0] } else { winding[0] }, ring, at, false));
        }
        if loops.is_empty() {
            // A whole sphere (no edge loops at all): a seam at u = 0 from its
            // south pole to its north pole, each with a degenerated edge (as
            // BRepPrim_Sphere).
            let edge_loops = face
                .loops
                .iter()
                .any(|l| matches!(t.loops()[l.index()], Loop::Edges { .. }));
            if let (Revolved::Sphere { .. }, None, false) = (rev, pole, edge_loops) {
                let half = std::f64::consts::FRAC_PI_2;
                let (south, north) = (extra.len(), extra.len() + 1);
                extra.push(rev.point(0.0, -half));
                extra.push(rev.point(0.0, half));
                // In the unoriented face the south pole's line runs +u.
                poles.insert(
                    fi,
                    vec![
                        (Vertex::Pole(south), -half, true),
                        (Vertex::Pole(north), half, false),
                    ],
                );
                wound.insert(
                    fi,
                    Seam {
                        u0: 0.0,
                        bottom: End {
                            vertex: Vertex::Pole(south),
                            height: -half,
                        },
                        top: End {
                            vertex: Vertex::Pole(north),
                            height: half,
                        },
                    },
                );
            }
            continue;
        }
        // The pole closes a band that winds once: it runs the other way.
        if let Some(v) = pole {
            let total: i32 = loops.iter().map(|l| l.0).sum();
            if total.abs() != 1 {
                return Err(unwritable("a pole on a face not wound once"));
            }
            let Some(apex) = rev.pole_v(total > 0) else {
                return Err(unwritable("a pole on a cylinder"));
            };
            loops.push((-total, None, vec![(v, Point2::new(0.0, apex))], true));
            poles.insert(fi, vec![(Vertex::Shared(v), apex, -total > 0)]);
        }
        let (Some(bottom), Some(top), 2) = (
            loops.iter().find(|l| l.0 > 0),
            loops.iter().find(|l| l.0 < 0),
            loops.len(),
        ) else {
            return Err(unwritable("a wound face without one loop each way"));
        };
        let same = |a: Point2, b: Point2| {
            let d = (a.x - b.x).rem_euclid(TAU);
            d.min(TAU - d) * rho(a.y).max(rho(b.y)) <= tolerance
        };
        // The seam position: a vertex both loops have, or any vertex of a
        // loop facing a ring, or where the rings' pcurves start.
        let ring_u = |e: usize| {
            t.face_fins(FaceId::new(fi))
                .into_iter()
                .flatten()
                .find(|f| f.edge.index() == e)
                .map(|f| f.pcurve.point(0.0).x)
        };
        // A pole meets the seam anywhere: the other loop decides.
        let pole_end = |l: &WoundLoop| l.3.then(|| (l.2[0].0, l.2[0].1.y));
        let choice = match (bottom.1, top.1) {
            _ if bottom.3 || top.3 => {
                let other = if bottom.3 { top } else { bottom };
                let at = match other.1 {
                    Some(e) => ring_u(e).map(|u| (u, None)),
                    None => other.2.first().map(|(w, q)| (q.x, Some((*w, q.y)))),
                };
                at.map(|(u, vertex)| {
                    if bottom.3 {
                        (u, pole_end(bottom), vertex)
                    } else {
                        (u, vertex, pole_end(top))
                    }
                })
            }
            (Some(e), Some(_)) => ring_u(e).map(|u| (u, None, None)),
            (Some(_), None) => top.2.first().map(|(w, q)| (q.x, None, Some((*w, q.y)))),
            (None, Some(_)) => bottom.2.first().map(|(v, p)| (p.x, Some((*v, p.y)), None)),
            (None, None) => bottom.2.iter().find_map(|(v, p)| {
                top.2
                    .iter()
                    .find(|(_, q)| same(*p, *q))
                    .map(|(w, q)| (p.x, Some((*v, p.y)), Some((*w, q.y))))
            }),
        };
        let Some((u, bottom_vertex, top_vertex)) = choice else {
            return Err(unwritable(
                "a wound face whose loops share no seam position",
            ));
        };
        let u0 = u.rem_euclid(TAU);
        let mut end =
            |ring: Option<usize>, vertex: Option<(usize, f64)>| -> Result<End, BrepError> {
                if let Some((v, height)) = vertex {
                    return Ok(End {
                        vertex: Vertex::Shared(v),
                        height,
                    });
                }
                let e = ring.expect("a ring loop");
                let fin = t
                    .face_fins(FaceId::new(fi))
                    .into_iter()
                    .flatten()
                    .find(|f| f.edge.index() == e)
                    .expect("the ring's fin on this face");
                let height = fin.pcurve.point(0.0).y;
                // The ring's own parameter at the seam point.
                let (Curve3::Circle {
                    frame: cf,
                    radius: cr,
                }
                | Curve3::CircularArc {
                    frame: cf,
                    radius: cr,
                    ..
                }) = &t.edges()[e].curve
                else {
                    return Err(unwritable("a ring edge that is not a circle"));
                };
                if (cr - rho(height)).abs() > tolerance {
                    return Err(unwritable("a ring edge off its wall's radius"));
                }
                let dir = frame.x() * u0.cos() + frame.normal().cross(frame.x()) * u0.sin();
                let ny = cf.normal().cross(cf.x());
                let mut angle = dir.dot(ny).atan2(dir.dot(cf.x()));
                if let Curve3::CircularArc { sweep_angle, .. } = &t.edges()[e].curve {
                    if *sweep_angle < 0.0 {
                        angle = -angle;
                    }
                }
                ring_start.entry(e).or_insert(angle);
                Ok(End {
                    vertex: Vertex::Ring(e),
                    height,
                })
            };
        let bottom_end = end(bottom.1, bottom_vertex)?;
        let top_end = end(top.1, top_vertex)?;
        wound.insert(
            fi,
            Seam {
                u0,
                bottom: bottom_end,
                top: top_end,
            },
        );
    }
    // Seams in v (a torus): a face wound in v gets a latitude seam at the v0
    // where its ring loops start, from the -v loop's seam vertex to the +v
    // loop's (the rectangle's bottom, its forward use at v0); a whole torus
    // gets both seams through one vertex of its own.
    let mut vseams: BTreeMap<usize, VSeam> = BTreeMap::new();
    let mut whole_tori: BTreeMap<usize, usize> = BTreeMap::new();
    for (fi, face) in t.faces().iter().enumerate() {
        let Surface::Torus {
            frame,
            major,
            minor,
        } = &face.surface
        else {
            continue;
        };
        let flip = face.sense == Orientation::Reversed;
        let (mut plus, mut minus, mut edge_loops) = (None, None, false);
        for l in &face.loops {
            let Loop::Edges { fins, winding } = &t.loops()[l.index()] else {
                continue;
            };
            edge_loops = true;
            if winding[1] == 0 {
                continue;
            }
            let first = &t.fins()[fins[0].index()];
            let e = first.edge.index();
            if winding[1].abs() != 1 || fins.len() != 1 || !t.edges()[e].is_ring() {
                return Err(unwritable("a loop winding in v that is not one ring"));
            }
            let start = first.pcurve.point(0.0);
            let wv = if flip { -winding[1] } else { winding[1] };
            if wv > 0 {
                plus = Some((e, start));
            } else {
                minus = Some((e, start));
            }
        }
        if !edge_loops {
            whole_tori.insert(fi, extra.len());
            extra.push(face.surface.point(Point2::default()));
            continue;
        }
        let (Some(plus), Some(minus)) = (plus, minus) else {
            continue;
        };
        let v0 = plus.1.y.rem_euclid(TAU);
        let (u_from, u_to) = (minus.1.x, plus.1.x);
        if u_to <= u_from {
            return Err(unwritable("a torus wedge running backwards in u"));
        }
        for (e, at) in [plus, minus] {
            let (Curve3::Circle {
                frame: cf,
                radius: cr,
            }
            | Curve3::CircularArc {
                frame: cf,
                radius: cr,
                ..
            }) = &t.edges()[e].curve
            else {
                return Err(unwritable("a ring edge that is not a circle"));
            };
            if (cr - minor).abs() > tolerance {
                return Err(unwritable("a ring edge off its tube's radius"));
            }
            let dir = face.surface.point(Point2::new(at.x, v0)) - cf.origin();
            let ny = cf.normal().cross(cf.x());
            let mut angle = dir.dot(ny).atan2(dir.dot(cf.x()));
            if let Curve3::CircularArc { sweep_angle, .. } = &t.edges()[e].curve {
                if *sweep_angle < 0.0 {
                    angle = -angle;
                }
            }
            ring_start.entry(e).or_insert(angle);
        }
        let _ = (frame, major);
        vseams.insert(
            fi,
            VSeam {
                v0,
                from: Vertex::Ring(minus.0),
                to: Vertex::Ring(plus.0),
                range: [u_from, u_to],
            },
        );
    }
    // A ring edge no seam meets (bounding a plane or in a wire, S6) is
    // closed at the start of its record, as BRepBuilderAPI_MakeEdge closes a
    // circle.
    for (e, edge) in t.edges().iter().enumerate() {
        let start = match &edge.curve {
            Curve3::Circle { .. } => 0.0,
            Curve3::CircularArc {
                start_angle,
                sweep_angle,
                ..
            } if *sweep_angle < 0.0 => -start_angle,
            Curve3::CircularArc { start_angle, .. } => *start_angle,
            _ => continue,
        };
        if edge.is_ring() {
            ring_start.entry(e).or_insert(start);
        }
    }
    // Vertices, then seam vertices of ring edges.
    let mut vertex_record: BTreeMap<usize, usize> = BTreeMap::new();
    let vertex = |p: Point3, records: &mut Records| {
        records.push(format!(
            "Ve\n{tol}\n{}\n0 0\n\n0101101\n*",
            nums(&p.to_array())
        ))
    };
    for (v, x) in t.vertices().iter().enumerate() {
        vertex_record.insert(v, vertex(x.position, &mut records));
    }
    let mut ring_vertex: BTreeMap<usize, usize> = BTreeMap::new();
    for (e, start) in &ring_start {
        let p = curve_at(&t.edges()[*e].curve, *start);
        ring_vertex.insert(*e, vertex(p, &mut records));
    }
    let extra_vertex: Vec<usize> = extra.iter().map(|p| vertex(*p, &mut records)).collect();
    // Surfaces.
    for face in t.faces() {
        surfaces.push(match &face.surface {
            Surface::BSpline(spline) => bspline_surface_record(spline),
            Surface::Plane(f) => format!(
                "1 {} {} {} {}",
                nums(&f.origin().to_array()),
                nums(&f.normal().to_array()),
                nums(&f.x().to_array()),
                nums(&f.normal().cross(f.x()).to_array())
            ),
            Surface::Cylinder { frame: f, radius } => format!(
                "2 {} {} {} {} {}",
                nums(&f.origin().to_array()),
                nums(&f.normal().to_array()),
                nums(&f.x().to_array()),
                nums(&f.normal().cross(f.x()).to_array()),
                num(*radius)
            ),
            Surface::Cone {
                frame: f,
                radius,
                half_angle,
            } => format!(
                "3 {} {} {} {} {} {}",
                nums(&f.origin().to_array()),
                nums(&f.normal().to_array()),
                nums(&f.x().to_array()),
                nums(&f.normal().cross(f.x()).to_array()),
                num(*radius),
                num(*half_angle)
            ),
            Surface::Torus {
                frame: f,
                major,
                minor,
            } => format!(
                "5 {} {} {} {} {} {}",
                nums(&f.origin().to_array()),
                nums(&f.normal().to_array()),
                nums(&f.x().to_array()),
                nums(&f.normal().cross(f.x()).to_array()),
                num(*major),
                num(*minor)
            ),
            Surface::Sphere { frame: f, radius } => format!(
                "4 {} {} {} {} {}",
                nums(&f.origin().to_array()),
                nums(&f.normal().to_array()),
                nums(&f.x().to_array()),
                nums(&f.normal().cross(f.x()).to_array()),
                num(*radius)
            ),
        });
    }
    // Edge geometry and pcurves per face.
    let geometry: Vec<EdgeGeometry> = t
        .edges()
        .iter()
        .enumerate()
        .map(|(e, edge)| curve_record(&edge.curve, ring_start.get(&e).copied()))
        .collect();
    let mut reps: Vec<Vec<String>> = vec![Vec::new(); t.edges().len()];
    for (fi, face) in t.faces().iter().enumerate() {
        for fins in t.face_fins(FaceId::new(fi)) {
            for fin in fins {
                let e = fin.edge.index();
                let range = geometry[e].range;
                let mut pcurve = fin.pcurve.clone();
                // On a wound face every pcurve lies in [u0, u0 + 2 pi].
                if let (Some(seam), Curve2::LineSegment { start, end }) = (wound.get(&fi), &pcurve)
                {
                    let low = start.x.min(end.x);
                    let k = -((low - seam.u0) / TAU + 1e-9).floor() * TAU;
                    pcurve = Curve2::LineSegment {
                        start: Point2::new(start.x + k, start.y),
                        end: Point2::new(end.x + k, end.y),
                    };
                }
                // Likewise in [v0, v0 + 2 pi] on a face wound in v: a ring
                // starting a rounding below its seam's v0 (reduced into
                // [0, 2 pi)) lies a period below it otherwise.
                if let (Some(seam), Curve2::LineSegment { start, end }) = (vseams.get(&fi), &pcurve)
                {
                    let low = start.y.min(end.y);
                    let k = -((low - seam.v0) / TAU + 1e-9).floor() * TAU;
                    pcurve = Curve2::LineSegment {
                        start: Point2::new(start.x, start.y + k),
                        end: Point2::new(end.x, end.y + k),
                    };
                }
                let forward = fin.sense == Orientation::Forward;
                curves2d.push(pcurve_record(&face.surface, &pcurve, forward, range)?);
                reps[e].push(format!(
                    "2 {} {} 0 {}",
                    curves2d.len(),
                    fi + 1,
                    nums(&range)
                ));
            }
        }
    }
    // Seam edges: from the bottom loop's seam vertex to the top one's.
    let record_of = |vertex: Vertex| match vertex {
        Vertex::Shared(v) => vertex_record[&v],
        Vertex::Ring(e) => ring_vertex[&e],
        Vertex::Pole(k) => extra_vertex[k],
    };
    let seam_vertex = |end: &End| record_of(end.vertex);
    let mut seam_of: BTreeMap<usize, usize> = BTreeMap::new();
    for (fi, seam) in &wound {
        let rev = Revolved::of(&t.faces()[*fi].surface).expect("wound faces are revolved");
        let (u0, low, high) = (seam.u0, seam.bottom.height, seam.top.height);
        let (a, b) = (rev.point(u0, low), rev.point(u0, high));
        if b.distance(a) <= tolerance {
            return Err(unwritable("a wound face of no height at its seam"));
        }
        // The seam's 3D curve and range: a ruling from low to high at unit
        // speed, or on a sphere the meridian circle in the angle v itself.
        let range = match rev {
            Revolved::Ruled { .. } => {
                let l = b.distance(a);
                let d = (b - a) * (1.0 / l);
                curves.push(format!("1 {} {}", nums(&a.to_array()), nums(&d.to_array())));
                let dv = (high - low).signum();
                curves2d.push(format!("1 {} 0 {}", nums(&[u0 + TAU, low]), num(dv)));
                curves2d.push(format!("1 {} 0 {}", nums(&[u0, low]), num(dv)));
                [0.0, l]
            }
            Revolved::Sphere { frame, .. } | Revolved::Torus { frame, .. } => {
                if high <= low {
                    return Err(unwritable("a seam running down"));
                }
                let e = frame.x() * u0.cos() + frame.normal().cross(frame.x()) * u0.sin();
                let n = frame.normal();
                // The meridian circle: the sphere's, or the tube's about
                // O + major e(u0); its parameter is v.
                let (centre, radius) = match rev {
                    Revolved::Torus { major, minor, .. } => (frame.origin() + e * major, minor),
                    Revolved::Sphere { radius, .. } => (frame.origin(), radius),
                    Revolved::Ruled { .. } => unreachable!("ruled seams are lines"),
                };
                curves.push(format!(
                    "2 {} {} {} {} {}",
                    nums(&centre.to_array()),
                    nums(&e.cross(n).to_array()),
                    nums(&e.to_array()),
                    nums(&n.to_array()),
                    num(radius)
                ));
                curves2d.push(format!("1 {} 0 1", nums(&[u0 + TAU, 0.0])));
                curves2d.push(format!("1 {} 0 1", nums(&[u0, 0.0])));
                [low, high]
            }
        };
        let c3 = curves.len();
        let (p1, p2) = (curves2d.len() - 1, curves2d.len());
        let edge = records.push(format!(
            "Ed\n {tol} 1 1 0\n1 {c3} 0 {}\n3 {p1} {p2} CN {} 0 {}\n0\n\n0101000\n{{v+{}}} 0 {{v-{}}} 0 *",
            nums(&range),
            fi + 1,
            nums(&range),
            seam_vertex(&seam.bottom),
            seam_vertex(&seam.top),
        ));
        seam_of.insert(*fi, edge);
    }
    // Seams in v and whole tori.
    let latitude = |frame: crate::Frame3, major: f64, minor: f64, v: f64| {
        let centre = frame.origin() + frame.normal() * (minor * v.sin());
        format!(
            "2 {} {} {} {} {}",
            nums(&centre.to_array()),
            nums(&frame.normal().to_array()),
            nums(&frame.x().to_array()),
            nums(&frame.normal().cross(frame.x()).to_array()),
            num(major + minor * v.cos())
        )
    };
    let mut vseam_of: BTreeMap<usize, usize> = BTreeMap::new();
    for (fi, seam) in &vseams {
        let Surface::Torus {
            frame,
            major,
            minor,
        } = &t.faces()[*fi].surface
        else {
            unreachable!("v seams are on tori");
        };
        curves.push(latitude(*frame, *major, *minor, seam.v0));
        let c3 = curves.len();
        curves2d.push(format!("1 {} 1 0", nums(&[0.0, seam.v0])));
        curves2d.push(format!("1 {} 1 0", nums(&[0.0, seam.v0 + TAU])));
        let (p1, p2) = (curves2d.len() - 1, curves2d.len());
        let edge = records.push(format!(
            "Ed\n {tol} 1 1 0\n1 {c3} 0 {}\n3 {p1} {p2} CN {} 0 {}\n0\n\n0101000\n{{v+{}}} 0 {{v-{}}} 0 *",
            nums(&seam.range),
            fi + 1,
            nums(&seam.range),
            record_of(seam.from),
            record_of(seam.to),
        ));
        vseam_of.insert(*fi, edge);
    }
    let mut torus_seams: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
    for (fi, k) in &whole_tori {
        let Surface::Torus {
            frame,
            major,
            minor,
        } = &t.faces()[*fi].surface
        else {
            unreachable!("whole tori");
        };
        let at = extra_vertex[*k];
        let full = nums(&[0.0, TAU]);
        // The meridian circle at u = 0 (parameter v), its forward use at
        // u = 2 pi; the latitude circle at v = 0 (parameter u), its forward
        // use at v = 0 (as BRepPrim_OneAxis builds them).
        let centre = frame.origin() + frame.x() * *major;
        curves.push(format!(
            "2 {} {} {} {} {}",
            nums(&centre.to_array()),
            nums(&frame.x().cross(frame.normal()).to_array()),
            nums(&frame.x().to_array()),
            nums(&frame.normal().to_array()),
            num(*minor)
        ));
        let cu = curves.len();
        curves2d.push(format!("1 {} 0 1", nums(&[TAU, 0.0])));
        curves2d.push("1 0 0 0 1".to_string());
        let (pu1, pu2) = (curves2d.len() - 1, curves2d.len());
        let u_seam = records.push(format!(
            "Ed\n {tol} 1 1 0\n1 {cu} 0 {full}\n3 {pu1} {pu2} CN {} 0 {full}\n0\n\n0101000\n{{v+{at}}} 0 {{v-{at}}} 0 *",
            fi + 1,
        ));
        curves.push(latitude(*frame, *major, *minor, 0.0));
        let cv = curves.len();
        curves2d.push("1 0 0 1 0".to_string());
        curves2d.push(format!("1 {} 1 0", nums(&[0.0, TAU])));
        let (pv1, pv2) = (curves2d.len() - 1, curves2d.len());
        let v_seam = records.push(format!(
            "Ed\n {tol} 1 1 0\n1 {cv} 0 {full}\n3 {pv1} {pv2} CN {} 0 {full}\n0\n\n0101000\n{{v+{at}}} 0 {{v-{at}}} 0 *",
            fi + 1,
        ));
        torus_seams.insert(*fi, (u_seam, v_seam));
    }
    // Degenerated edges: a cone's or sphere's poles, one turn at the pole's
    // v from the seam.
    let mut pole_edge: BTreeMap<usize, Vec<(usize, bool)>> = BTreeMap::new();
    for (fi, list) in &poles {
        let u0 = wound[fi].u0;
        for (v, apex, forward) in list {
            curves2d.push(format!("1 {} 1 0", nums(&[u0, *apex])));
            let apex_vertex = record_of(*v);
            let edge = records.push(format!(
                "Ed\n {tol} 1 1 1\n2 {} {} 0 0 {}\n0\n\n0101000\n{{v+{apex_vertex}}} 0 {{v-{apex_vertex}}} 0 *",
                curves2d.len(),
                fi + 1,
                num(TAU),
            ));
            pole_edge.entry(*fi).or_default().push((edge, *forward));
        }
    }
    // Edges.
    let mut edge_record: BTreeMap<usize, usize> = BTreeMap::new();
    for (e, edge) in t.edges().iter().enumerate() {
        // A wire edge (S6) has no fin: only its 3D curve.
        curves.push(geometry[e].record.clone());
        let mut text = format!(
            "Ed\n {tol} 1 1 0\n1 {} 0 {}\n",
            curves.len(),
            nums(&geometry[e].range)
        );
        for r in &reps[e] {
            text.push_str(r);
            text.push('\n');
        }
        text.push_str("0\n\n0101000\n");
        let (s, en) = match (edge.start, edge.end) {
            (Some(s), Some(en)) => (vertex_record[&s.index()], vertex_record[&en.index()]),
            _ => {
                let v = ring_vertex[&e];
                (v, v)
            }
        };
        let _ = write!(text, "{{v+{s}}} 0 {{v-{en}}} 0 *");
        edge_record.insert(e, records.push(text));
    }
    // Wires and faces, in the unoriented face: a reversed face's loops are
    // traversed backwards.
    let mut face_record: BTreeMap<usize, usize> = BTreeMap::new();
    for (fi, face) in t.faces().iter().enumerate() {
        let flip = face.sense == Orientation::Reversed;
        let mut loops: Vec<(Vec<(usize, bool)>, i32)> = Vec::new();
        for l in &face.loops {
            let Loop::Edges { fins, winding } = &t.loops()[l.index()] else {
                if pole_edge.contains_key(&fi) {
                    continue;
                }
                return Err(unwritable("a vertex loop"));
            };
            // In the oriented face, each fin and, where the loop passes
            // through a pole, OCCT's degenerated edge from one fin's end to
            // the next one's start along the pole's line.
            let mut uses: Vec<(usize, bool)> = Vec::new();
            for (k, fk) in fins.iter().enumerate() {
                let f = &t.fins()[fk.index()];
                uses.push((f.edge.index(), f.sense == Orientation::Forward));
                let next = &t.fins()[fins[(k + 1) % fins.len()].index()];
                let (a, b) = (f.pcurve.point(1.0), next.pcurve.point(0.0));
                if !pole_pass(&face.surface, a, b) {
                    continue;
                }
                if winding[0] != 0 || seam_of.contains_key(&fi) {
                    return Err(unwritable("a loop through a pole on a wound face"));
                }
                let edge = &t.edges()[f.edge.index()];
                let end = if f.sense == Orientation::Forward {
                    edge.end
                } else {
                    edge.start
                };
                let Some(v) = end else {
                    return Err(unwritable("a ring edge through a pole"));
                };
                let du = b.x - a.x;
                curves2d.push(format!("1 {} {} 0", nums(&[a.x, a.y]), num(du.signum())));
                let at = vertex_record[&v.index()];
                let record = records.push(format!(
                    "Ed\n {tol} 1 1 1\n2 {} {} 0 0 {}\n0\n\n0101000\n{{v+{at}}} 0 {{v-{at}}} 0 *",
                    curves2d.len(),
                    fi + 1,
                    num(du.abs()),
                ));
                let key = t.edges().len() + edge_record.len();
                edge_record.insert(key, record);
                uses.push((key, true));
            }
            let mut uses: Vec<(usize, bool)> =
                uses.into_iter().map(|(e, fwd)| (e, fwd != flip)).collect();
            if flip {
                uses.reverse();
            }
            loops.push((uses, if flip { -winding[0] } else { winding[0] }));
        }
        let mut wires = Vec::new();
        let use_text = |(e, fwd): &(usize, bool)| {
            format!("{{{}{}}} 0", if *fwd { "+" } else { "-" }, edge_record[e])
        };
        if let Some((u_seam, v_seam)) = torus_seams.get(&fi) {
            wires.push(records.push(format!(
                "Wi\n\n0101100\n{{+{v_seam}}} 0 {{+{u_seam}}} 0 {{-{v_seam}}} 0 {{-{u_seam}}} 0 *"
            )));
        } else if let Some(seam_edge) = vseam_of.get(&fi) {
            // One wire: the seam at v0 from the -v ring to the +v ring, the +v
            // ring, the seam back at v0 + 2 pi, the -v ring.
            let ring_use = |up: bool| -> Option<String> {
                face.loops.iter().find_map(|l| match &t.loops()[l.index()] {
                    Loop::Edges { fins, winding } if winding[1] != 0 => {
                        let wv = if flip { -winding[1] } else { winding[1] };
                        let f = &t.fins()[fins[0].index()];
                        let fwd = (f.sense == Orientation::Forward) != flip;
                        ((wv > 0) == up).then(|| {
                            format!(
                                "{{{}{}}} 0",
                                if fwd { "+" } else { "-" },
                                edge_record[&f.edge.index()]
                            )
                        })
                    }
                    _ => None,
                })
            };
            let (Some(up), Some(down)) = (ring_use(true), ring_use(false)) else {
                return Err(unwritable("a torus wedge without both rings"));
            };
            wires.push(records.push(format!(
                "Wi\n\n0101100\n{{+{seam_edge}}} 0 {up} {{-{seam_edge}}} 0 {down} *"
            )));
        } else if let Some(seam_edge) = seam_of.get(&fi) {
            // One wire: the bottom loop (+u) from the seam vertex, the seam
            // up, the top loop (-u) from its seam vertex, the seam down.
            let seam = &wound[&fi];
            let starting = |uses: &[(usize, bool)], end: &End| -> Vec<(usize, bool)> {
                let at = match end.vertex {
                    Vertex::Ring(_) | Vertex::Pole(_) => 0,
                    Vertex::Shared(v) => uses
                        .iter()
                        .position(|(e, fwd)| {
                            let edge = &t.edges()[*e];
                            let s = if *fwd { edge.start } else { edge.end };
                            s.map(|s| s.index()) == Some(v)
                        })
                        .unwrap_or(0),
                };
                uses[at..].iter().chain(&uses[..at]).copied().collect()
            };
            let face_poles = pole_edge.get(&fi).cloned().unwrap_or_default();
            let wound_text = |up: bool, end: &End| -> Option<Vec<String>> {
                if let Some((edge, forward)) = face_poles.iter().find(|(_, f)| *f == up) {
                    let sign = if *forward { "+" } else { "-" };
                    return Some(vec![format!("{{{sign}{edge}}} 0")]);
                }
                let (uses, _) = loops.iter().find(|(_, w)| (*w > 0) == up && *w != 0)?;
                Some(starting(uses, end).iter().map(use_text).collect())
            };
            let (Some(b), Some(u)) = (wound_text(true, &seam.bottom), wound_text(false, &seam.top))
            else {
                return Err(unwritable("a wound face without both wound loops"));
            };
            let mut text = b;
            text.push(format!("{{+{seam_edge}}} 0"));
            text.extend(u);
            text.push(format!("{{-{seam_edge}}} 0"));
            wires.push(records.push(format!("Wi\n\n0101100\n{} *", text.join(" "))));
            for (uses, _) in loops.iter().filter(|(_, w)| *w == 0) {
                let text: Vec<String> = uses.iter().map(use_text).collect();
                wires.push(records.push(format!("Wi\n\n0101100\n{} *", text.join(" "))));
            }
        } else {
            for (uses, w) in loops {
                if w != 0 {
                    return Err(unwritable("a wound loop on a face without a seam"));
                }
                let text: Vec<String> = uses
                    .iter()
                    .map(|(e, fwd)| {
                        format!("{{{}{}}} 0", if *fwd { "+" } else { "-" }, edge_record[e])
                    })
                    .collect();
                wires.push(records.push(format!("Wi\n\n0101100\n{} *", text.join(" "))));
            }
        }
        let text: Vec<String> = wires.iter().map(|w| format!("{{+{w}}} 0")).collect();
        face_record.insert(
            fi,
            records.push(format!(
                "Fa\n0 {tol} {} 0\n\n0101000\n{} *",
                fi + 1,
                text.join(" ")
            )),
        );
    }
    // Shells and solids: one solid per solid region, its shells' sides with
    // their orientation against the surface.
    let mut solids = Vec::new();
    for region in t.regions().iter().filter(|r| r.kind == RegionKind::Solid) {
        let mut shells = Vec::new();
        for s in &region.shells {
            let shell = &t.shells()[s.index()];
            let text: Vec<String> = shell
                .sides
                .iter()
                .map(|(f, side)| {
                    let face = &t.faces()[f.index()];
                    let forward = (face.sense == Orientation::Forward) == (*side == Side::Front);
                    format!(
                        "{{{}{}}} 0",
                        if forward { "+" } else { "-" },
                        face_record[&f.index()]
                    )
                })
                .collect();
            shells.push(records.push(format!("Sh\n\n0101100\n{} *", text.join(" "))));
        }
        let text: Vec<String> = shells.iter().map(|s| format!("{{+{s}}} 0")).collect();
        solids.push(records.push(format!("So\n\n0100100\n{} *", text.join(" "))));
    }
    let compound = |items: &[String], records: &mut Records| {
        records.push(format!("Co\n\n1100000\n{} *", items.join(" ")))
    };
    let (root, forward) = match solids.as_slice() {
        [one] => (*one, true),
        [] => free_root(t, &face_record, &edge_record, &vertex_record, &mut records)?,
        many => {
            let text: Vec<String> = many.iter().map(|s| format!("{{+{s}}} 0")).collect();
            (compound(&text, &mut records), true)
        }
    };
    let _ = EdgeId::new(0);
    // Records are numbered backwards: the first written is n, the last 1.
    let n = records.shapes.len();
    let number = |i: usize| n - i;
    let mut out = String::from(
        "DBRep_DrawableShape\n\nCASCADE Topology V1, (c) Matra-Datavision\nLocations 0\n",
    );
    let _ = writeln!(out, "Curve2ds {}", curves2d.len());
    for c in &curves2d {
        let _ = writeln!(out, "{c}");
    }
    let _ = writeln!(out, "Curves {}", curves.len());
    for c in &curves {
        let _ = writeln!(out, "{c}");
    }
    let _ = writeln!(out, "Polygon3D 0\nPolygonOnTriangulations 0");
    let _ = writeln!(out, "Surfaces {}", surfaces.len());
    for s in &surfaces {
        let _ = writeln!(out, "{s}");
    }
    let _ = writeln!(out, "Triangulations 0\n\nTShapes {n}");
    for record in &records.shapes {
        // Resolve {v+i}, {v-i}, {+i} and {-i} placeholders into numbers.
        let mut text = String::new();
        let mut rest = record.as_str();
        while let Some(open) = rest.find('{') {
            text.push_str(&rest[..open]);
            let close = open + rest[open..].find('}').expect("closed placeholder");
            let inner = &rest[open + 1..close];
            let (sign, index) = match inner.strip_prefix('v') {
                Some(v) => (&v[..1], &v[1..]),
                None => (&inner[..1], &inner[1..]),
            };
            let i: usize = index.parse().expect("record index");
            let _ = write!(text, "{sign}{}", number(i));
            rest = &rest[close + 1..];
        }
        text.push_str(rest);
        let _ = writeln!(out, "{text}");
    }
    let _ = writeln!(
        out,
        "\n{}{} 0 \n0",
        if forward { "+" } else { "-" },
        number(root)
    );
    Ok(out)
}

/// The root of a body without a solid region (S6): a closed shell (of a
/// bounded void) or several in a compound; an open sheet's one face, or its
/// faces in one open shell, each by its front side; a wire's one edge, or
/// its edges in one wire, each oriented to continue the last; an acorn's
/// vertex.
fn free_root(
    t: &Topology,
    face_record: &BTreeMap<usize, usize>,
    edge_record: &BTreeMap<usize, usize>,
    vertex_record: &BTreeMap<usize, usize>,
    records: &mut Records,
) -> Result<(usize, bool), BrepError> {
    let sign = |forward: bool| if forward { "+" } else { "-" };
    let side_text = |f: FaceId, side: Side| {
        let face = &t.faces()[f.index()];
        let forward = (face.sense == Orientation::Forward) == (side == Side::Front);
        format!("{{{}{}}} 0", sign(forward), face_record[&f.index()])
    };
    match t.class() {
        BodyClass::Sheet => {
            let closed: Vec<usize> = t.regions()[1..]
                .iter()
                .flat_map(|r| r.shells.iter())
                .map(|s| {
                    let text: Vec<String> = t.shells()[s.index()]
                        .sides
                        .iter()
                        .map(|(f, side)| side_text(*f, *side))
                        .collect();
                    records.push(format!("Sh\n\n0101100\n{} *", text.join(" ")))
                })
                .collect();
            match closed.as_slice() {
                [one] => return Ok((*one, true)),
                [] => {}
                many => {
                    let text: Vec<String> = many.iter().map(|s| format!("{{+{s}}} 0")).collect();
                    return Ok((
                        records.push(format!("Co\n\n1100000\n{} *", text.join(" "))),
                        true,
                    ));
                }
            }
            if t.faces().len() == 1 {
                let forward = t.faces()[0].sense == Orientation::Forward;
                return Ok((face_record[&0], forward));
            }
            let text: Vec<String> = (0..t.faces().len())
                .map(|f| side_text(FaceId::new(f), Side::Front))
                .collect();
            Ok((
                records.push(format!("Sh\n\n0101000\n{} *", text.join(" "))),
                true,
            ))
        }
        BodyClass::Wire => {
            let edges: Vec<usize> = t
                .shells()
                .iter()
                .flat_map(|s| s.wire_edges.iter().map(|e| e.index()))
                .collect();
            if let [one] = edges.as_slice() {
                return Ok((edge_record[one], true));
            }
            let ends = |e: usize| (t.edges()[e].start, t.edges()[e].end);
            let mut at = None;
            let mut text = Vec::new();
            for (k, e) in edges.iter().enumerate() {
                let (a, b) = ends(*e);
                let forward = match at {
                    Some(v) => a == Some(v) || b != Some(v),
                    // The first edge runs towards the second.
                    None => edges.get(k + 1).is_none_or(|n| {
                        let (c, d) = ends(*n);
                        b == c || b == d || !(a == c || a == d)
                    }),
                };
                at = if forward { b } else { a };
                text.push(format!("{{{}{}}} 0", sign(forward), edge_record[e]));
            }
            Ok((
                records.push(format!("Wi\n\n0101100\n{} *", text.join(" "))),
                true,
            ))
        }
        BodyClass::Acorn => {
            let v = t
                .shells()
                .iter()
                .find_map(|s| s.acorn_vertices.first())
                .expect("an acorn has its vertex");
            Ok((vertex_record[&v.index()], true))
        }
        BodyClass::Solid | BodyClass::General => {
            Err(unwritable("a topology without a solid region"))
        }
    }
}
