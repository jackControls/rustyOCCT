//! Deflection-controlled, watertight tessellation (T-a of `REVIEW_NOTES.md`,
//! `TESSELLATION.md`).
//!
//! Every edge is discretized once, uniformly in its fraction; every face
//! that uses the edge takes exactly that polyline as boundary, so faces
//! share nodes and the mesh is watertight by construction. A face is
//! triangulated in a planar chart of its domain on the universal cover: a
//! plane in its frame coordinates, a periodic face whose loops do not wind
//! in the sinusoidal chart `((u - u_c) |S_u|(v), s(v))`, which collapses a
//! pole to a point, a face wound in `u` in an annulus chart `P(v) (cos u,
//! sin u)` and a torus face wound in `v` the same way with `u` and `v`
//! exchanged; a whole sphere or torus is a structured grid. No seam is
//! meshed. A constrained Delaunay triangulation of the chart polygons
//! (`cdt.rs`) is refined by Steiner points until every triangle's certified
//! deflection and normal turn (`bounds.rs`) are within the request.
//!
//! The bound reported for a triangle is certified: under the map
//! `Σ λ_i X_i ↦ S(Σ λ_i p_i)`, `X_i` its nodes and `p_i` their parameter
//! points, no point moves farther than it; for an edge segment the same with
//! the curve. Reference: OCCT `BRepMesh_IncrementalMesh` and its
//! `IMeshTools_Parameters` (`Deflection`, `Angle`), which measure deflection
//! by sampling and mesh seams; see `SOURCE_MAP.md`.
use crate::predicates::{orient2d_finite, Orientation2};
use crate::topology::validate::continuity;
use crate::topology::{
    Curve2, Curve3, EdgeId, Face, FaceId, FinId, Loop, Orientation, RegionKind, Surface, Topology,
};
use crate::{Error, Point2, Point3, Result, Vec3};
use std::collections::{BTreeMap, VecDeque};
use std::f64::consts::{FRAC_PI_2, TAU};

mod bounds;
mod cdt;
mod spline;

use cdt::{Cdt, CORNERS};
use spline::{CurveCells, SurfaceCells};

/// More nodes than this is `ComputationLimit`.
pub const MAX_NODES: usize = 4_000_000;
/// Doublings of a face's curved edges when its chords cross or nest wrongly.
const MAX_ROUNDS: usize = 12;
/// Fractions of the request an edge's own and thin-triangle bounds may use.
const CURVE_SHARE: f64 = 0.9;
const THIN_SHARE: f64 = 0.45;

/// Absolute bounds of a tessellation, as OCCT's `IMeshTools_Parameters`
/// `Deflection` (linear, in the body's length unit) and `Angle` (radians).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parameters {
    deflection: f64,
    angle: f64,
}

impl Parameters {
    /// A finite deflection above zero and an angle in `(0, π/2]`.
    pub fn new(deflection: f64, angle: f64) -> Result<Self> {
        if !deflection.is_finite() || !angle.is_finite() {
            return Err(Error::NonFinite("tessellation parameters"));
        }
        if deflection <= 0.0 || angle <= 0.0 || angle > FRAC_PI_2 {
            return Err(Error::OutOfDomain("tessellation parameters"));
        }
        Ok(Self { deflection, angle })
    }
    pub fn deflection(self) -> f64 {
        self.deflection
    }
    pub fn angle(self) -> f64 {
        self.angle
    }
}

/// A triangle's or segment's certified deflection and normal (or tangent)
/// turn.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bound {
    pub deflection: f64,
    pub angle: f64,
}

/// A face's triangles, a range of `Mesh::triangles`.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceMesh {
    pub face: FaceId,
    pub triangles: std::ops::Range<usize>,
}

/// An edge's polyline: nodes from the curve's start to its end (a ring edge
/// closes on its first node), and each segment's bound.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeMesh {
    pub edge: EdgeId,
    pub nodes: Vec<usize>,
    pub segments: Vec<Bound>,
}

/// A watertight triangle mesh of a topology's faces with its edges'
/// polylines. Every triangle's normal (right-hand rule) leaves the solid
/// region behind its face, or follows the face's oriented normal between
/// void regions.
#[derive(Debug, Clone, PartialEq)]
pub struct Mesh {
    pub nodes: Vec<Point3>,
    pub triangles: Vec<[usize; 3]>,
    pub triangle_bounds: Vec<Bound>,
    /// One per face, in face order.
    pub faces: Vec<FaceMesh>,
    /// One per edge, in edge order.
    pub edges: Vec<EdgeMesh>,
    /// The node of each vertex.
    pub vertices: Vec<usize>,
    /// The largest certified deflection of a triangle or segment.
    pub deflection: f64,
    /// The largest certified normal or tangent turn.
    pub angle: f64,
}

/// Tessellate every face and edge of a topology within `parameters`.
pub fn tessellate(topology: &Topology, parameters: Parameters) -> Result<Mesh> {
    let ranges: Vec<(f64, f64)> = topology
        .faces()
        .iter()
        .map(|face| v_range(topology, face))
        .collect();
    let mut fin_face = vec![usize::MAX; topology.fins().len()];
    for (f, face) in topology.faces().iter().enumerate() {
        for l in &face.loops {
            if let Loop::Edges { fins, .. } = &topology.loops()[l.index()] {
                for fin in fins {
                    fin_face[fin.index()] = f;
                }
            }
        }
    }
    let splines = Splines::new(topology, &fin_face)?;
    let mut multiplier = vec![1usize; topology.edges().len()];
    for _ in 0..MAX_ROUNDS {
        let mut nodes: Vec<Point3> = topology.vertices().iter().map(|v| v.position).collect();
        let mut polylines = Vec::with_capacity(topology.edges().len());
        for (e, edge) in topology.edges().iter().enumerate() {
            let mut n = segment_count(topology, e, &fin_face, &ranges, &splines, parameters)?
                * multiplier[e];
            let cells = splines.edges[e].as_ref();
            let mark = nodes.len();
            let mut polyline = discretize(topology, e, &edge.curve, cells, n, &mut nodes)?;
            // T-b: a spline edge's segments, and the boundary segments of
            // spline faces, checked; each failure doubles the count.
            let mut doublings = 0;
            while splines.checked(topology, e, &fin_face)
                && !splines.edge_ok(topology, e, &polyline, &fin_face, parameters)
            {
                doublings += 1;
                n = n.saturating_mul(2);
                if doublings > MAX_ROUNDS || n > MAX_NODES {
                    return Err(Error::ComputationLimit("tessellation spline edge"));
                }
                nodes.truncate(mark);
                polyline = discretize(topology, e, &edge.curve, cells, n, &mut nodes)?;
            }
            polylines.push(polyline);
            if nodes.len() > MAX_NODES {
                return Err(Error::ComputationLimit("tessellation size"));
            }
        }
        let mut triangles = Vec::new();
        let mut triangle_bounds = Vec::new();
        let mut faces = Vec::new();
        let mut refine = Vec::new();
        for (f, face) in topology.faces().iter().enumerate() {
            let start = triangles.len();
            let mesher = FaceMesher {
                topology,
                face,
                polylines: &polylines,
                parameters,
                outward: outward_sign(topology, face),
                spline: splines.faces[f].as_ref(),
            };
            match mesher.mesh(&mut nodes, &mut triangles, &mut triangle_bounds) {
                Ok(()) => {}
                Err(FaceError::Conflict) => {
                    triangles.truncate(start);
                    triangle_bounds.truncate(start);
                    refine.push(f);
                }
                Err(FaceError::Fatal(error)) => return Err(error),
            }
            faces.push(FaceMesh {
                face: FaceId::new(f),
                triangles: start..triangles.len(),
            });
        }
        if refine.is_empty() {
            let mesh = Mesh {
                deflection: triangle_bounds
                    .iter()
                    .chain(polylines.iter().flat_map(|p| p.segments.iter()))
                    .map(|b| b.deflection)
                    .fold(0.0, f64::max),
                angle: triangle_bounds
                    .iter()
                    .chain(polylines.iter().flat_map(|p| p.segments.iter()))
                    .map(|b| b.angle)
                    .fold(0.0, f64::max),
                nodes,
                triangles,
                triangle_bounds,
                faces,
                vertices: (0..topology.vertices().len()).collect(),
                edges: polylines
                    .into_iter()
                    .enumerate()
                    .map(|(e, p)| EdgeMesh {
                        edge: EdgeId::new(e),
                        nodes: p.nodes,
                        segments: p.segments,
                    })
                    .collect(),
            };
            check_mesh(topology, &mesh, parameters)?;
            return Ok(mesh);
        }
        // Coarse chords crossed or changed nesting: refine every curved
        // edge of those faces and start again.
        let mut any = false;
        for f in refine {
            for e in face_edges(topology, &topology.faces()[f]) {
                if !matches!(topology.edges()[e].curve, Curve3::LineSegment { .. }) {
                    multiplier[e] *= 2;
                    any = true;
                }
            }
        }
        if !any {
            return Err(Error::InvalidTopology("face boundary crosses itself"));
        }
    }
    Err(Error::ComputationLimit("tessellation boundary refinement"))
}

/// The spline cells of a topology (T-b): each spline edge's and each
/// nonperiodic spline face's, their C1 continuity checked (R4), and the
/// speed of each spline pcurve on a curved analytic surface.
struct Splines {
    edges: Vec<Option<CurveCells>>,
    faces: Vec<Option<SurfaceCells>>,
    pcurves: Vec<Option<f64>>,
}

impl Splines {
    fn new(topology: &Topology, fin_face: &[usize]) -> Result<Self> {
        let not_c1 = Error::OutOfDomain("tessellation of a spline that is not C1");
        let mut faces = Vec::with_capacity(topology.faces().len());
        for face in topology.faces() {
            faces.push(match &face.surface {
                Surface::BSpline(s) => {
                    if s.u_knots().is_periodic() || s.v_knots().is_periodic() {
                        return Err(Error::OutOfDomain(
                            "tessellation of a periodic spline surface",
                        ));
                    }
                    if !continuity::surface_c1(s) {
                        return Err(not_c1);
                    }
                    Some(SurfaceCells::new(s)?)
                }
                _ => None,
            });
        }
        let mut edges = Vec::with_capacity(topology.edges().len());
        for edge in topology.edges() {
            edges.push(match &edge.curve {
                Curve3::BSpline(span) => {
                    if !continuity::curve_c1(span.curve(), span.range(), span.is_closed_period()) {
                        return Err(not_c1);
                    }
                    Some(CurveCells::new(span.curve(), span.range())?)
                }
                _ => None,
            });
        }
        let mut pcurves = Vec::with_capacity(topology.fins().len());
        for (k, fin) in topology.fins().iter().enumerate() {
            let curved = fin_face[k] != usize::MAX
                && !matches!(
                    topology.faces()[fin_face[k]].surface,
                    Surface::Plane(_) | Surface::BSpline(_)
                );
            pcurves.push(match &fin.pcurve {
                Curve2::BSpline(span) if curved => {
                    Some(CurveCells::new(span.curve().as_curve3(), span.range())?.speed())
                }
                _ => None,
            });
        }
        Ok(Self {
            edges,
            faces,
            pcurves,
        })
    }

    fn face_of(&self, fin: FinId, fin_face: &[usize]) -> Option<&SurfaceCells> {
        let f = fin_face[fin.index()];
        if f == usize::MAX {
            None
        } else {
            self.faces[f].as_ref()
        }
    }

    /// Whether an edge's polyline needs T-b's checks: a spline edge, or an
    /// edge with a fin on a spline face.
    fn checked(&self, topology: &Topology, e: usize, fin_face: &[usize]) -> bool {
        self.edges[e].is_some()
            || topology.edges()[e]
                .fins
                .iter()
                .any(|f| self.face_of(*f, fin_face).is_some())
    }

    /// A spline edge's segments within the request, and on every spline face
    /// using the edge each boundary segment's thin-triangle condition over
    /// the cells of its chord's box, within a share of the request.
    fn edge_ok(
        &self,
        topology: &Topology,
        e: usize,
        polyline: &Polyline,
        fin_face: &[usize],
        parameters: Parameters,
    ) -> bool {
        let (delta, theta) = (parameters.deflection, parameters.angle);
        if self.edges[e].is_some()
            && polyline
                .segments
                .iter()
                .any(|b| !(b.deflection <= delta && b.angle <= theta))
        {
            return false;
        }
        for id in &topology.edges()[e].fins {
            let Some(cells) = self.face_of(*id, fin_face) else {
                continue;
            };
            let fin = &topology.fins()[id.index()];
            let forward = fin.sense == Orientation::Forward;
            let points: Vec<Point2> = polyline
                .fractions
                .iter()
                .map(|&t| {
                    let s = if forward {
                        t
                    } else if t == 0.0 {
                        1.0
                    } else if t == 1.0 {
                        0.0
                    } else {
                        1.0 - t
                    };
                    cells.clamp(fin.pcurve.point(s))
                })
                .collect();
            for w in points.windows(2) {
                let (deviation, turn, _) = cells.box_bound(
                    w[0].x.min(w[1].x),
                    w[0].x.max(w[1].x),
                    w[0].y.min(w[1].y),
                    w[0].y.max(w[1].y),
                );
                if !(deviation <= THIN_SHARE * delta && turn <= THIN_SHARE * theta) {
                    return false;
                }
            }
        }
        true
    }
}

fn face_edges(topology: &Topology, face: &Face) -> Vec<usize> {
    let mut out = Vec::new();
    for l in &face.loops {
        if let Loop::Edges { fins, .. } = &topology.loops()[l.index()] {
            out.extend(fins.iter().map(|f| topology.fins()[f.index()].edge.index()));
        }
    }
    out
}

/// The sign taking the parametric normal `S_u x S_v` to the triangles'
/// normal: the face's sense, reversed when the solid region lies on its back
/// side only (the triangles then leave the material).
fn outward_sign(topology: &Topology, face: &Face) -> f64 {
    let kind = |s: crate::topology::ShellId| {
        topology.regions()[topology.shells()[s.index()].region.index()].kind
    };
    let material = if kind(face.front) == RegionKind::Void && kind(face.back) == RegionKind::Solid {
        -1.0
    } else {
        1.0
    };
    material * face.sense.sign()
}

// ------------------------------------------------------------------ edges

struct Polyline {
    nodes: Vec<usize>,
    fractions: Vec<f64>,
    segments: Vec<Bound>,
}

/// The poles' `v` on a sphere or a cone.
fn poles(surface: &Surface) -> Vec<f64> {
    match surface {
        Surface::Sphere { .. } => vec![-FRAC_PI_2, FRAC_PI_2],
        Surface::Cone {
            radius, half_angle, ..
        } => vec![-radius / half_angle.sin()],
        _ => Vec::new(),
    }
}

/// The range of `v` a face's pcurves reach (arcs by their extreme points),
/// with its poles.
fn v_range(topology: &Topology, face: &Face) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut add = |v: f64| {
        lo = lo.min(v);
        hi = hi.max(v);
    };
    for l in &face.loops {
        match &topology.loops()[l.index()] {
            Loop::Edges { fins, .. } => {
                for fin in fins {
                    match &topology.fins()[fin.index()].pcurve {
                        Curve2::LineSegment { start, end } => {
                            add(start.y);
                            add(end.y);
                        }
                        Curve2::CircularArc { center, radius, .. } => {
                            add(center.y - radius.abs());
                            add(center.y + radius.abs());
                        }
                        Curve2::EllipseArc { center, minor, .. } => {
                            add(center.y - minor.abs());
                            add(center.y + minor.abs());
                        }
                        Curve2::Sinusoid { a, .. } => {
                            let amplitude = a[1].hypot(a[2]) * (1.0 + 1e-12);
                            add(a[0] - amplitude);
                            add(a[0] + amplitude);
                        }
                        Curve2::BSpline(span) => {
                            // The curve lies in its control hull.
                            for p in span.curve().poles() {
                                add(p.y);
                            }
                        }
                    }
                }
            }
            Loop::Vertex(_) => {
                for p in poles(&face.surface) {
                    add(p);
                }
            }
        }
    }
    match face.surface {
        Surface::Sphere { .. } if face.loops.is_empty() => (-FRAC_PI_2, FRAC_PI_2),
        Surface::Torus { .. } if face.loops.is_empty() => (0.0, TAU),
        _ if lo > hi => (0.0, 0.0),
        Surface::Sphere { .. } => (lo.max(-FRAC_PI_2), hi.min(FRAC_PI_2)),
        _ => (lo, hi),
    }
}

/// An edge's segment count: its curve's deflection and turn within the
/// request, and each curved face's thin-triangle condition (a triangle on a
/// boundary segment is never smaller than the segment) within a share of it.
fn segment_count(
    topology: &Topology,
    e: usize,
    fin_face: &[usize],
    ranges: &[(f64, f64)],
    splines: &Splines,
    parameters: Parameters,
) -> Result<usize> {
    let edge = &topology.edges()[e];
    let (delta, theta) = (parameters.deflection, parameters.angle);
    let mut n = 1.0f64;
    let closed = edge.is_ring() || edge.start == edge.end;
    // The largest semi-axis, the turn per unit angle, the sweep.
    if let Some((frame, radius, ratio, sweep)) = match &edge.curve {
        Curve3::Circle { frame, radius } => Some((frame, radius.abs(), 1.0, TAU)),
        Curve3::CircularArc {
            frame,
            radius,
            sweep_angle,
            ..
        } => Some((frame, radius.abs(), 1.0, *sweep_angle)),
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            sweep_angle,
            ..
        } => {
            let (big, small) = (major.abs().max(minor.abs()), major.abs().min(minor.abs()));
            Some((frame, big, big / small * (1.0 + 1e-12), *sweep_angle))
        }
        Curve3::LineSegment { .. } | Curve3::BSpline(_) => None,
    } {
        let sigma = bounds::frame_norm(frame);
        let phi = sweep.abs() * (1.0 + 1e-12);
        n = n
            .max(phi * ratio / theta)
            .max(phi * (sigma * radius / (8.0 * CURVE_SHARE * delta)).sqrt());
    }
    // T-b: `D2 h² / 8` within the request over the span's parameter length;
    // each segment's turn and bound are checked once discretized.
    if let (Curve3::BSpline(span), Some(cells)) = (&edge.curve, &splines.edges[e]) {
        let [a, b] = span.range();
        let h = bounds::sum_up(b, -a);
        n = n.max(h * (cells.curvature() / (8.0 * CURVE_SHARE * delta)).sqrt());
    }
    if closed {
        n = n.max(3.0);
    }
    for fin in &edge.fins {
        let f = fin_face[fin.index()];
        if f == usize::MAX {
            continue;
        }
        let face = &topology.faces()[f];
        if matches!(face.surface, Surface::Plane(_)) {
            continue;
        }
        let (du, dv) = match &topology.fins()[fin.index()].pcurve {
            Curve2::LineSegment { start, end } => {
                ((end.x - start.x).abs(), (end.y - start.y).abs())
            }
            Curve2::CircularArc {
                radius,
                sweep_angle,
                ..
            } => {
                let l = (radius * sweep_angle).abs();
                (l, l)
            }
            Curve2::EllipseArc {
                major,
                minor,
                sweep_angle,
                ..
            } => {
                let l = (major.abs().max(minor.abs()) * sweep_angle).abs();
                (l, l)
            }
            Curve2::Sinusoid { sweep, a, .. } => (
                sweep.abs(),
                (a[1].hypot(a[2]) * sweep).abs() * (1.0 + 1e-12),
            ),
            // T-b: its parameter length times its speed, in u and in v.
            Curve2::BSpline(span) => match splines.pcurves[fin.index()] {
                Some(speed) => {
                    let [a, b] = span.range();
                    let l = bounds::sum_up(b, -a) * speed * (1.0 + 1e-12);
                    (l, l)
                }
                None => (0.0, 0.0),
            },
        };
        let (v0, v1) = ranges[f];
        let [a, b, c, nu, nv] = bounds::coefficients(&face.surface, v0, v1);
        let q = a * du * du + 2.0 * b * du * dv + c * dv * dv;
        n = n
            .max((q / (8.0 * THIN_SHARE * delta)).sqrt() * (1.0 + 1e-9))
            .max((nu * du + nv * dv) / (THIN_SHARE * theta) * (1.0 + 1e-9));
    }
    if !n.is_finite() || n > MAX_NODES as f64 {
        return Err(Error::ComputationLimit("tessellation size"));
    }
    Ok(n.ceil() as usize)
}

/// An edge's nodes at fractions `k / n` (its vertices' nodes at the ends) and
/// each segment's certified bound.
fn discretize(
    topology: &Topology,
    e: usize,
    curve: &Curve3,
    cells: Option<&CurveCells>,
    n: usize,
    nodes: &mut Vec<Point3>,
) -> Result<Polyline> {
    let edge = &topology.edges()[e];
    let unrepresentable = Error::Unrepresentable("tessellation node");
    let mut ids = Vec::with_capacity(n + 1);
    let mut fractions = Vec::with_capacity(n + 1);
    // Each node's distance from the exact curve point at its fraction.
    let mut errors = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let t = if k == n { 1.0 } else { k as f64 / n as f64 };
        // A spline's node carries its span's parameter at the fraction.
        let exact = match (curve, cells) {
            (Curve3::BSpline(span), Some(cells)) => cells.point(span.parameter(t)),
            _ => bounds::curve_point(curve, t),
        }
        .ok_or(unrepresentable.clone())?;
        let end = match (k, edge.start, edge.end) {
            (0, Some(v), _) => Some(v.index()),
            (k, _, Some(v)) if k == n => Some(v.index()),
            (k, None, None) if k == n => Some(ids[0]),
            _ => None,
        };
        let id = match end {
            Some(id) => id,
            None => {
                let (p, _) = bounds::settle(&exact);
                nodes.push(p);
                nodes.len() - 1
            }
        };
        errors.push(bounds::gap(nodes[id], &exact));
        ids.push(id);
        fractions.push(t);
    }
    let mut segments = Vec::with_capacity(n);
    for k in 0..n {
        let dt = bounds::sum_up(fractions[k + 1], -fractions[k]);
        let (deviation, angle) = match (curve, cells) {
            (Curve3::BSpline(span), Some(cells)) => cells.segment(
                span.parameter(fractions[k]),
                span.parameter(fractions[k + 1]),
            ),
            _ => bounds::segment_bound(curve, dt),
        };
        let deflection = bounds::sum_up(deviation, errors[k].max(errors[k + 1]));
        if !deflection.is_finite() {
            return Err(unrepresentable);
        }
        segments.push(Bound { deflection, angle });
    }
    Ok(Polyline {
        nodes: ids,
        fractions,
        segments,
    })
}

// ------------------------------------------------------------------ charts

/// `|S_u|` up to the frame's norm and the arc length `s` along `v`.
fn metric(surface: &Surface, v: f64) -> (f64, f64) {
    match surface {
        Surface::Cylinder { radius, .. } => (radius.abs(), v),
        Surface::Cone {
            radius, half_angle, ..
        } => ((radius + v * half_angle.sin()).abs(), v),
        Surface::Sphere { radius, .. } => ((radius * v.cos()).max(0.0), radius * v),
        Surface::Torus { major, minor, .. } => (major + minor * v.cos(), minor * v),
        _ => (1.0, v),
    }
}

fn arc_length_inverse(surface: &Surface, s: f64) -> f64 {
    match surface {
        Surface::Sphere { radius, .. } => (s / radius).clamp(-FRAC_PI_2, FRAC_PI_2),
        Surface::Torus { minor, .. } => s / minor,
        _ => s,
    }
}

#[derive(Debug, Clone, Copy)]
enum Chart {
    /// `(u - cu, v - cv)`.
    Plane { cu: f64, cv: f64 },
    /// `((u - uc) g(v), s(v) - sc)`: poles collapse to points.
    Sinusoidal { uc: f64, sc: f64 },
    /// `P(v) (cos u, sin u)`, `P = exp(k (s - sref) / gref)`, or `|s - s_p|`
    /// around a pole on the side `side` (+1: the band above it).
    RingU {
        k: f64,
        sref: f64,
        gref: f64,
        pole: Option<(f64, f64)>,
    },
    /// A torus wound in `v`: `exp(k (u - uref) scale) (cos v, sin v)`.
    RingV { k: f64, uref: f64, scale: f64 },
    /// A spline face (T-b): `(gu (u - cu), gv (v - cv))`.
    Affine { cu: f64, cv: f64, gu: f64, gv: f64 },
}

impl Chart {
    /// +1 when the chart keeps the parameter orientation.
    fn sign(self) -> f64 {
        match self {
            Chart::Plane { .. }
            | Chart::Sinusoidal { .. }
            | Chart::RingV { .. }
            | Chart::Affine { .. } => 1.0,
            Chart::RingU { pole: None, .. } => -1.0,
            Chart::RingU {
                pole: Some((_, side)),
                ..
            } => -side,
        }
    }

    fn forward(self, surface: &Surface, uv: Point2, pole: Option<f64>) -> Point2 {
        match self {
            Chart::Plane { cu, cv } => Point2::new(uv.x - cu, uv.y - cv),
            Chart::Sinusoidal { uc, sc } => match pole {
                Some(vp) => Point2::new(0.0, metric(surface, vp).1 - sc),
                None => {
                    let (g, s) = metric(surface, uv.y);
                    Point2::new((uv.x - uc) * g, s - sc)
                }
            },
            Chart::RingU {
                k,
                sref,
                gref,
                pole: axis,
            } => {
                if pole.is_some() {
                    return Point2::new(0.0, 0.0);
                }
                let s = metric(surface, uv.y).1;
                let rho = match axis {
                    Some((sp, _)) => (s - sp).abs(),
                    None => (k * (s - sref) / gref).exp(),
                };
                let (sin, cos) = uv.x.sin_cos();
                Point2::new(rho * cos, rho * sin)
            }
            Chart::RingV { k, uref, scale } => {
                let rho = (k * (uv.x - uref) * scale).exp();
                let (sin, cos) = uv.y.sin_cos();
                Point2::new(rho * cos, rho * sin)
            }
            Chart::Affine { cu, cv, gu, gv } => Point2::new((uv.x - cu) * gu, (uv.y - cv) * gv),
        }
    }

    fn inverse(self, surface: &Surface, p: Point2) -> Point2 {
        match self {
            Chart::Plane { cu, cv } => Point2::new(p.x + cu, p.y + cv),
            Chart::Sinusoidal { uc, sc } => {
                let v = arc_length_inverse(surface, p.y + sc);
                let g = metric(surface, v).0;
                Point2::new(if g > 0.0 { uc + p.x / g } else { uc }, v)
            }
            Chart::RingU {
                k,
                sref,
                gref,
                pole,
            } => {
                let rho = p.x.hypot(p.y);
                let s = match pole {
                    Some((sp, side)) => sp + side * rho,
                    None => sref + gref * rho.ln() / k,
                };
                Point2::new(p.y.atan2(p.x), arc_length_inverse(surface, s))
            }
            Chart::RingV { k, uref, scale } => {
                let rho = p.x.hypot(p.y);
                Point2::new(uref + rho.ln() / (k * scale), p.y.atan2(p.x))
            }
            Chart::Affine { cu, cv, gu, gv } => Point2::new(p.x / gu + cu, p.y / gv + cv),
        }
    }
}

// ------------------------------------------------------------------ faces

enum FaceError {
    /// Chords that cross, touch or nest wrongly: refine the face's edges.
    Conflict,
    Fatal(Error),
}

impl From<cdt::Conflict> for FaceError {
    fn from(_: cdt::Conflict) -> Self {
        FaceError::Conflict
    }
}

/// A node of a face: its parameter point on the cover (any `u` at a pole).
#[derive(Debug, Clone, Copy)]
struct Local {
    node: usize,
    uv: Point2,
    pole: Option<f64>,
}

struct FaceMesher<'a> {
    topology: &'a Topology,
    face: &'a Face,
    polylines: &'a [Polyline],
    parameters: Parameters,
    outward: f64,
    /// A spline face's cells (T-b).
    spline: Option<&'a SurfaceCells>,
}

impl FaceMesher<'_> {
    fn surface(&self) -> &Surface {
        &self.face.surface
    }

    /// The exact point of the face's surface at a parameter point, enclosed.
    fn surface_point(&self, uv: Point2) -> Option<[crate::certified::Fast; 3]> {
        match self.spline {
            Some(cells) => cells.jet(uv).map(|j| j[0]),
            None => bounds::surface_point(self.surface(), uv),
        }
    }

    fn fatal<T>(what: &'static str) -> std::result::Result<T, FaceError> {
        Err(FaceError::Fatal(Error::OutOfDomain(what)))
    }

    fn mesh(
        &self,
        nodes: &mut Vec<Point3>,
        triangles: &mut Vec<[usize; 3]>,
        bounds_out: &mut Vec<Bound>,
    ) -> std::result::Result<(), FaceError> {
        if self.face.loops.is_empty() {
            return self.grid(nodes, triangles, bounds_out);
        }
        let surface = self.surface();
        let periodic_v = surface.is_periodic_v();
        let mut rings: Vec<Vec<Local>> = Vec::new();
        let mut pole_vertices = Vec::new();
        let (mut wound_u, mut wound_v) = (0, 0);
        let pole_values = poles(surface);
        for l in &self.face.loops {
            match &self.topology.loops()[l.index()] {
                Loop::Edges { fins, winding } => {
                    match *winding {
                        [0, 0] => {}
                        [1 | -1, 0] if surface.is_periodic() => wound_u += 1,
                        [0, 1 | -1] if periodic_v => wound_v += 1,
                        _ => return Self::fatal("tessellation of a loop wound twice or both ways"),
                    }
                    rings.push(self.boundary(fins, &pole_values));
                }
                Loop::Vertex(v) => {
                    let position = self.topology.vertices()[v.index()].position;
                    let nearest = pole_values.iter().copied().min_by(|a, b| {
                        let d = |p: f64| {
                            bounds::surface_point(surface, Point2::new(0.0, p))
                                .map(|q| bounds::settle(&q).0.distance(position))
                                .unwrap_or(f64::INFINITY)
                        };
                        d(*a).total_cmp(&d(*b))
                    });
                    let Some(vp) = nearest else {
                        return Self::fatal("tessellation of a vertex loop that is not a pole");
                    };
                    pole_vertices.push(Local {
                        node: v.index(),
                        uv: Point2::new(0.0, vp),
                        pole: Some(vp),
                    });
                }
            }
        }
        if wound_u > 0 && wound_v > 0 {
            return Self::fatal("tessellation of a loop wound twice or both ways");
        }
        let wound = wound_u + wound_v;
        if (wound > 0 && wound + pole_vertices.len() != 2)
            || (wound == 0 && !pole_vertices.is_empty())
        {
            return Self::fatal("tessellation of this face's loops");
        }
        let chart = self.chart(&rings, &pole_vertices, wound_u > 0, wound_v > 0);
        let mut locals: Vec<Local> = rings.iter().flatten().copied().collect();
        locals.extend(pole_vertices.iter().copied());
        {
            // A node twice at one parameter point: loops that touch. The
            // same vertex a period apart on the cover (a wall pinched to a
            // point where a plane touches its cap's circle, S8a.2) is two
            // points of the chart.
            let mut seen: Vec<(usize, Point2)> = locals.iter().map(|l| (l.node, l.uv)).collect();
            seen.sort_by_key(|x| x.0);
            let touching = seen.windows(2).any(|w| {
                w[0].0 == w[1].0
                    && (w[0].1.x - w[1].1.x).abs() <= 1e-9 * (1.0 + w[0].1.x.abs())
                    && (w[0].1.y - w[1].1.y).abs() <= 1e-9 * (1.0 + w[0].1.y.abs())
            });
            if touching {
                return Self::fatal("tessellation of a face whose loops touch");
            }
        }
        let points: Vec<Point2> = locals
            .iter()
            .map(|l| chart.forward(surface, l.uv, l.pole))
            .collect();
        if points.iter().any(|p| !(p.x.is_finite() && p.y.is_finite())) {
            return Err(FaceError::Fatal(Error::Unrepresentable(
                "tessellation chart",
            )));
        }
        let mut cdt = Cdt::new(&points)?;
        let mut start = 0u32;
        for ring in &rings {
            let n = ring.len() as u32;
            for k in 0..n {
                cdt.constrain(CORNERS + start + k, CORNERS + start + (k + 1) % n)?;
            }
            start += n;
        }
        cdt.make_delaunay();
        cdt.classify()?;
        // Material lies left of each loop's traversal seen along the oriented
        // normal: in the chart, left when its sign and the sense agree.
        let left = chart.sign() * self.face.sense.sign() > 0.0;
        let mut start = 0u32;
        for ring in &rings {
            let n = ring.len() as u32;
            for k in 0..n {
                let (a, b) = (CORNERS + start + k, CORNERS + start + (k + 1) % n);
                let (inner, outer) = if left { (a, b) } else { (b, a) };
                let inside = |x: u32, y: u32| {
                    cdt.directed_edge(x, y)
                        .map(|(t, _)| cdt.alive_inside(t))
                        .unwrap_or(false)
                };
                if !inside(inner, outer) || inside(outer, inner) {
                    return Err(FaceError::Conflict);
                }
            }
            start += n;
        }
        let mut gaps: Vec<f64> = locals
            .iter()
            .map(|l| self.gap(nodes[l.node], l.uv, l.pole))
            .collect();
        if gaps
            .iter()
            .any(|g| g.is_nan() || *g > 0.5 * self.parameters.deflection)
        {
            return Err(FaceError::Fatal(Error::InvalidTopology(
                "a node farther than half the deflection from its face",
            )));
        }
        let curved = !matches!(surface, Surface::Plane(_));
        if curved {
            // A spline face's parameter span: refinement below 2^-30 of it
            // in both directions is a budget failure.
            let span = |f: fn(&Local) -> f64| {
                let (lo, hi) = locals
                    .iter()
                    .map(f)
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
                        (l.min(x), h.max(x))
                    });
                (hi - lo) * 2f64.powi(-30)
            };
            let floor = (span(|l| l.uv.x), span(|l| l.uv.y));
            self.refine(&mut cdt, chart, &mut locals, &mut gaps, nodes, floor)?;
        }
        let flip = chart.sign() * self.outward < 0.0;
        for t in 0..cdt.slots() {
            if !cdt.alive_inside(t) {
                continue;
            }
            let (ok, bound) = self.evaluate(&cdt, t, chart, &locals, &gaps, nodes);
            if curved && !ok {
                return Err(FaceError::Fatal(Error::ComputationLimit(
                    "tessellation refinement",
                )));
            }
            let v = cdt.vertices(t).map(|v| locals[(v - CORNERS) as usize].node);
            triangles.push(if flip { [v[0], v[2], v[1]] } else { v });
            bounds_out.push(bound);
        }
        Ok(())
    }

    /// A loop's nodes in traversal order with their pcurve points: each fin
    /// runs its edge's polyline forwards or backwards, the pcurve at the
    /// matching fraction (a reversed use pairs pcurve fraction `s` with edge
    /// fraction `1 - s`), its last node left to the next fin.
    fn boundary(&self, fins: &[crate::topology::FinId], pole_values: &[f64]) -> Vec<Local> {
        let mut out = Vec::new();
        for fin in fins {
            let fin = &self.topology.fins()[fin.index()];
            let line = &self.polylines[fin.edge.index()];
            let n = line.nodes.len() - 1;
            let forward = fin.sense == Orientation::Forward;
            for k in 0..n {
                let j = if forward { k } else { n - k };
                let s = if forward {
                    line.fractions[j]
                } else if j == n {
                    0.0
                } else {
                    1.0 - line.fractions[j]
                };
                let uv = fin.pcurve.point(s);
                let uv = self.spline.map_or(uv, |cells| cells.clamp(uv));
                let node = line.nodes[j];
                let is_vertex = node < self.topology.vertices().len();
                let pole = pole_values
                    .iter()
                    .copied()
                    .find(|p| is_vertex && (uv.y - p).abs() <= 1e-12 * p.abs().max(1.0));
                out.push(Local { node, uv, pole });
            }
        }
        out
    }

    fn chart(&self, rings: &[Vec<Local>], poles: &[Local], wound_u: bool, wound_v: bool) -> Chart {
        let surface = self.surface();
        let points: Vec<Point2> = rings
            .iter()
            .flatten()
            .filter(|l| l.pole.is_none())
            .map(|l| l.uv)
            .collect();
        let span = |f: fn(&Point2) -> f64| {
            points
                .iter()
                .map(f)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
                    (l.min(x), h.max(x))
                })
        };
        let (u0, u1) = span(|p| p.x);
        let (v0, v1) = span(|p| p.y);
        let (um, vm) = (0.5 * u0 + 0.5 * u1, 0.5 * v0 + 0.5 * v1);
        if matches!(surface, Surface::Plane(_)) {
            return Chart::Plane { cu: um, cv: vm };
        }
        if let Some(cells) = self.spline {
            // Scaled by the partials' lengths at the region's centre, so the
            // chart is about isometric there.
            let (gu, gv) = cells.speeds(cells.clamp(Point2::new(um, vm)));
            return Chart::Affine {
                cu: um,
                cv: vm,
                gu,
                gv,
            };
        }
        if wound_v {
            let (major, minor) = match surface {
                Surface::Torus { major, minor, .. } => (*major, *minor),
                _ => (1.0, 1.0),
            };
            let scale = major / minor;
            let range = (u1 - u0) * scale;
            let k = if range > 4.0 { 4.0 / range } else { 1.0 };
            return Chart::RingV { k, uref: um, scale };
        }
        if !wound_u {
            return Chart::Sinusoidal {
                uc: um,
                sc: metric(surface, vm).1,
            };
        }
        if let Some(pole) = poles.first() {
            let vp = pole.pole.unwrap_or(0.0);
            let side = if vm > vp { 1.0 } else { -1.0 };
            return Chart::RingU {
                k: 1.0,
                sref: 0.0,
                gref: 1.0,
                pole: Some((metric(surface, vp).1, side)),
            };
        }
        let (s0, s1) = (metric(surface, v0).1, metric(surface, v1).1);
        let gref = metric(surface, vm)
            .0
            .max(metric(surface, v0).0)
            .max(metric(surface, v1).0)
            .max(f64::MIN_POSITIVE);
        let range = (s1 - s0) / gref;
        let k = if range > 4.0 { 4.0 / range } else { 1.0 };
        Chart::RingU {
            k,
            sref: metric(surface, vm).1,
            gref,
            pole: None,
        }
    }

    /// The certified distance from a node to the surface at its parameter
    /// point (at a pole, any `u`: the gap is taken at `u = 0` and at the
    /// triangle's own `u` when a triangle is evaluated).
    fn gap(&self, x: Point3, uv: Point2, pole: Option<f64>) -> f64 {
        let uv = pole.map_or(uv, |vp| Point2::new(uv.x, vp));
        match self.surface_point(uv) {
            Some(p) => bounds::gap(x, &p),
            None => f64::INFINITY,
        }
    }

    /// A triangle's parameter points lifted into one period, a pole taking
    /// the mean `u` of the others.
    fn lift(&self, locals: [&Local; 3]) -> [Point2; 3] {
        let surface = self.surface();
        let (pu, pv) = (surface.is_periodic(), surface.is_periodic_v());
        let reference = locals
            .iter()
            .find(|l| l.pole.is_none())
            .map(|l| l.uv)
            .unwrap_or(locals[0].uv);
        let near = |x: f64, r: f64, periodic: bool| {
            if periodic {
                x - TAU * ((x - r) / TAU).round()
            } else {
                x
            }
        };
        let mut out = locals
            .map(|l| Point2::new(near(l.uv.x, reference.x, pu), near(l.uv.y, reference.y, pv)));
        let others: Vec<f64> = locals
            .iter()
            .zip(out.iter())
            .filter(|(l, _)| l.pole.is_none())
            .map(|(_, p)| p.x)
            .collect();
        let mean = if others.is_empty() {
            reference.x
        } else {
            others.iter().sum::<f64>() / others.len() as f64
        };
        for (l, p) in locals.iter().zip(out.iter_mut()) {
            if let Some(vp) = l.pole {
                *p = Point2::new(mean, vp);
            }
        }
        out
    }

    /// Whether a triangle meets every criterion, and its certified bound.
    fn evaluate(
        &self,
        cdt: &Cdt,
        t: u32,
        chart: Chart,
        locals: &[Local],
        gaps: &[f64],
        nodes: &[Point3],
    ) -> (bool, Bound) {
        let v = cdt.vertices(t).map(|v| (v - CORNERS) as usize);
        let l = [&locals[v[0]], &locals[v[1]], &locals[v[2]]];
        let p = self.lift(l);
        let mut gap = 0.0f64;
        for k in 0..3 {
            gap = gap.max(gaps[v[k]]);
            if l[k].pole.is_some() {
                gap = gap.max(self.gap(nodes[l[k].node], p[k], l[k].pole));
            }
        }
        // A triangle with one vertex at a pole is certified by its fan map.
        let poles: Vec<usize> = (0..3).filter(|k| l[*k].pole.is_some()).collect();
        let pole = (poles.len() == 1).then(|| poles[0]);
        // A spline face's bounds and normal come from its cells (T-b).
        let (deviation, turn, spline_normal) = match self.spline {
            Some(cells) => cells.triangle_bound(p),
            None => {
                let (deviation, turn) = bounds::triangle_bound(self.surface(), p, pole);
                (deviation, turn, None)
            }
        };
        let bound = Bound {
            deflection: bounds::sum_up(deviation, gap),
            angle: turn,
        };
        if matches!(self.surface(), Surface::Plane(_)) {
            return (true, bound);
        }
        let expected = if chart.sign() > 0.0 {
            Orientation2::CounterClockwise
        } else {
            Orientation2::Clockwise
        };
        let within =
            bound.deflection <= self.parameters.deflection && bound.angle <= self.parameters.angle;
        if !within || orient2d_finite(p[0], p[1], p[2]) != expected {
            return (false, bound);
        }
        // The flat triangle's normal against the surface's at its middle, in
        // the parameter orientation.
        let x = l.map(|l| nodes[l.node]);
        let mut n = (x[1] - x[0]).cross(x[2] - x[0]);
        if chart.sign() < 0.0 {
            n = n * -1.0;
        }
        let mid = Point2::new(
            (p[0].x + p[1].x + p[2].x) / 3.0,
            (p[0].y + p[1].y + p[2].y) / 3.0,
        );
        let normal = match (self.spline, spline_normal) {
            (Some(_), Some(normal)) => normal,
            (Some(_), None) => return (false, bound),
            (None, _) => self.parametric_normal(mid),
        };
        (n.dot(normal) > 0.0, bound)
    }

    /// The direction of `S_u x S_v` (a cone's on either nappe).
    fn parametric_normal(&self, uv: Point2) -> Vec3 {
        let n = self.surface().normal(uv);
        match self.surface() {
            Surface::Cone {
                radius, half_angle, ..
            } if radius + uv.y * half_angle.sin() < 0.0 => n * -1.0,
            _ => n,
        }
    }

    /// The edge to split, opposite vertex `i`: the one contributing most to
    /// the triangle's bounds relative to the request (its own extents in the
    /// deviation and turn formulas; an edge at a pole has no `u` extent),
    /// the longest in space when none contributes.
    fn split_edge(&self, cdt: &Cdt, t: u32, locals: &[Local], nodes: &[Point3]) -> usize {
        let v = cdt.vertices(t).map(|v| (v - CORNERS) as usize);
        let l = [&locals[v[0]], &locals[v[1]], &locals[v[2]]];
        let p = self.lift(l);
        let (v0, v1) = p
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), q| {
                (a.min(q.y), b.max(q.y))
            });
        let [a, b, c, nu, nv] = match self.spline {
            // T-b: the cells over the triangle's box, the normal's rates
            // relative to |S_u x S_v| at its centre.
            Some(cells) => {
                let (u0, u1) = p
                    .iter()
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), q| {
                        (a.min(q.x), b.max(q.x))
                    });
                let [a, b, c, du, dv] = cells.coefficients(u0, u1, v0, v1);
                let m = cells
                    .normal(Point2::new(0.5 * u0 + 0.5 * u1, 0.5 * v0 + 0.5 * v1))
                    .map_or(0.0, |n| n.length())
                    .max(f64::MIN_POSITIVE);
                [a, b, c, (a * dv + du * b) / m, (b * dv + du * c) / m]
            }
            None => bounds::coefficients(self.surface(), v0, v1),
        };
        let (delta, theta) = (self.parameters.deflection, self.parameters.angle);
        let score = |k: usize| {
            let (m, n) = ((k + 1) % 3, (k + 2) % 3);
            let du = if l[m].pole.is_some() || l[n].pole.is_some() {
                0.0
            } else {
                (p[m].x - p[n].x).abs()
            };
            let dv = (p[m].y - p[n].y).abs();
            let deviation = (a * du * du + 2.0 * b * du * dv + c * dv * dv) / (8.0 * delta);
            deviation.max((nu * du + nv * dv) / theta)
        };
        let length = |k: usize| nodes[l[(k + 1) % 3].node].distance(nodes[l[(k + 2) % 3].node]);
        let best = (0..3)
            .max_by(|x, y| score(*x).total_cmp(&score(*y)))
            .unwrap();
        if score(best) > 0.0 {
            best
        } else {
            (0..3)
                .max_by(|x, y| length(*x).total_cmp(&length(*y)))
                .unwrap()
        }
    }

    fn refine(
        &self,
        cdt: &mut Cdt,
        chart: Chart,
        locals: &mut Vec<Local>,
        gaps: &mut Vec<f64>,
        nodes: &mut Vec<Point3>,
        floor: (f64, f64),
    ) -> std::result::Result<(), FaceError> {
        let surface = self.surface();
        let mut queue: VecDeque<(u32, u32)> = (0..cdt.slots())
            .filter(|t| cdt.alive_inside(*t))
            .map(|t| (t, cdt.stamp(t)))
            .collect();
        while let Some((t, stamp)) = queue.pop_front() {
            if !cdt.alive_inside(t) || cdt.stamp(t) != stamp {
                continue;
            }
            if self.evaluate(cdt, t, chart, locals, gaps, nodes).0 {
                continue;
            }
            if self.spline.is_some() {
                // T-b: a failing triangle this small will not certify (a
                // vanishing normal, a collapsed row).
                let q = self.lift(cdt.vertices(t).map(|v| &locals[(v - CORNERS) as usize]));
                let extent = |f: fn(&Point2) -> f64| {
                    let (lo, hi) = q
                        .iter()
                        .map(f)
                        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
                            (l.min(x), h.max(x))
                        });
                    hi - lo
                };
                if extent(|q| q.x) <= floor.0 && extent(|q| q.y) <= floor.1 {
                    return Err(FaceError::Fatal(Error::ComputationLimit(
                        "tessellation refinement",
                    )));
                }
            }
            let v = cdt.vertices(t);
            let fixed = cdt.fixed(t);
            let i = self.split_edge(cdt, t, locals, nodes);
            let c = v.map(|v| cdt.point(v));
            let centroid = Point2::new(
                (c[0].x + c[1].x + c[2].x) / 3.0,
                (c[0].y + c[1].y + c[2].y) / 3.0,
            );
            let mut candidates = Vec::with_capacity(2);
            if !fixed[i] {
                let (a, b) = (c[(i + 1) % 3], c[(i + 2) % 3]);
                candidates.push(Point2::new(0.5 * a.x + 0.5 * b.x, 0.5 * a.y + 0.5 * b.y));
            }
            candidates.push(centroid);
            let mut inserted = None;
            for p in candidates {
                if let Some(done) = cdt.insert_inside(p, t) {
                    inserted = Some((p, done));
                    break;
                }
            }
            let Some((p, (vertex, touched))) = inserted else {
                return Err(FaceError::Fatal(Error::ComputationLimit(
                    "tessellation refinement",
                )));
            };
            debug_assert_eq!(vertex as usize, CORNERS as usize + locals.len());
            let uv = chart.inverse(surface, p);
            let uv = self.spline.map_or(uv, |cells| cells.clamp(uv));
            let exact = self
                .surface_point(uv)
                .ok_or(FaceError::Fatal(Error::Unrepresentable(
                    "tessellation node",
                )))?;
            let (x, error) = bounds::settle(&exact);
            if !(error.is_finite() && x.x.is_finite() && x.y.is_finite() && x.z.is_finite()) {
                return Err(FaceError::Fatal(Error::Unrepresentable(
                    "tessellation node",
                )));
            }
            nodes.push(x);
            if nodes.len() > MAX_NODES {
                return Err(FaceError::Fatal(Error::ComputationLimit(
                    "tessellation size",
                )));
            }
            locals.push(Local {
                node: nodes.len() - 1,
                uv,
                pole: None,
            });
            gaps.push(error);
            for s in touched {
                if cdt.alive_inside(s) {
                    queue.push_back((s, cdt.stamp(s)));
                }
            }
        }
        Ok(())
    }

    /// A face without loops (the whole sphere or torus): a grid in `(u, v)`
    /// refined until every triangle's bounds hold.
    fn grid(
        &self,
        nodes: &mut Vec<Point3>,
        triangles: &mut Vec<[usize; 3]>,
        bounds_out: &mut Vec<Bound>,
    ) -> std::result::Result<(), FaceError> {
        let surface = self.surface();
        let sphere = match surface {
            Surface::Sphere { .. } => true,
            Surface::Torus { .. } => false,
            _ => return Self::fatal("tessellation of a face without loops"),
        };
        let (delta, theta) = (self.parameters.deflection, self.parameters.angle);
        let (v0, v1) = if sphere {
            (-FRAC_PI_2, FRAC_PI_2)
        } else {
            (0.0, TAU)
        };
        let [a, _, c, nu, nv] = bounds::coefficients(surface, v0, v1);
        // Equal steps h in u and v from (a + c) h^2 / 8 and (nu + nv) h
        // within the request (the mixed term is mostly tangential), shrunk
        // until every triangle's certified bounds hold.
        let mut h = (8.0 * CURVE_SHARE * delta / (a + c))
            .sqrt()
            .min(CURVE_SHARE * theta / (nu + nv));
        for _ in 0..16 {
            let steps = |span: f64, least: f64| ((span / h).ceil()).max(least);
            let (mu, mv) = (
                steps(TAU, 3.0),
                steps(v1 - v0, if sphere { 2.0 } else { 3.0 }),
            );
            if mu * mv > MAX_NODES as f64 {
                return Err(FaceError::Fatal(Error::ComputationLimit(
                    "tessellation size",
                )));
            }
            let (mu, mv) = (mu as usize, mv as usize);
            let base = nodes.len();
            let mut local = Vec::new();
            // Rows of v: a sphere's poles are single nodes.
            let rows: Vec<f64> = (0..=mv)
                .map(|j| {
                    if j == mv {
                        v1
                    } else {
                        v0 + (v1 - v0) * j as f64 / mv as f64
                    }
                })
                .collect();
            let mut index = vec![vec![0usize; mu]; mv + 1];
            let last_row = if sphere { mv } else { mv - 1 };
            for j in 0..=last_row {
                let pole = sphere && (j == 0 || j == mv);
                for i in 0..mu {
                    if pole && i > 0 {
                        index[j][i] = index[j][0];
                        continue;
                    }
                    let uv = Point2::new(
                        if pole {
                            0.0
                        } else {
                            TAU * i as f64 / mu as f64
                        },
                        rows[j],
                    );
                    let exact = bounds::surface_point(surface, uv).ok_or(FaceError::Fatal(
                        Error::Unrepresentable("tessellation node"),
                    ))?;
                    let (x, _) = bounds::settle(&exact);
                    local.push(x);
                    index[j][i] = base + local.len() - 1;
                }
            }
            if !sphere {
                index[mv] = index[0].clone();
            }
            let at = |i: usize, j: usize| index[j][i % mu];
            let position = |id: usize| {
                if id < base {
                    nodes[id]
                } else {
                    local[id - base]
                }
            };
            let mut out = Vec::new();
            let mut ok = true;
            for j in 0..mv {
                for i in 0..mu {
                    let (ua, ub) = (TAU * i as f64 / mu as f64, TAU * (i + 1) as f64 / mu as f64);
                    let (va, vb) = (rows[j], rows[j + 1]);
                    let quad = [
                        (at(i, j), Point2::new(ua, va)),
                        (at(i + 1, j), Point2::new(ub, va)),
                        (at(i + 1, j + 1), Point2::new(ub, vb)),
                        (at(i, j + 1), Point2::new(ua, vb)),
                    ];
                    let mut tris = vec![[quad[0], quad[1], quad[2]], [quad[0], quad[2], quad[3]]];
                    if sphere && j == 0 {
                        tris = vec![[quad[0], quad[2], quad[3]]];
                    } else if sphere && j == mv - 1 {
                        tris = vec![[quad[0], quad[1], quad[2]]];
                    }
                    for tri in tris {
                        let mut p = tri.map(|(_, uv)| uv);
                        // A pole's parameter point: the others' mean u.
                        for k in 0..3 {
                            if sphere && (p[k].y == v0 || p[k].y == v1) {
                                let others: Vec<f64> =
                                    (0..3).filter(|m| *m != k).map(|m| p[m].x).collect();
                                p[k].x = 0.5 * others[0] + 0.5 * others[1];
                            }
                        }
                        let mut gap = 0.0f64;
                        for k in 0..3 {
                            let exact = bounds::surface_point(surface, p[k]).ok_or(
                                FaceError::Fatal(Error::Unrepresentable("tessellation node")),
                            )?;
                            gap = gap.max(bounds::gap(position(tri[k].0), &exact));
                        }
                        let pole = (0..3).find(|k| sphere && (p[*k].y == v0 || p[*k].y == v1));
                        let (deviation, turn) = bounds::triangle_bound(surface, p, pole);
                        let bound = Bound {
                            deflection: bounds::sum_up(deviation, gap),
                            angle: turn,
                        };
                        ok &= bound.deflection <= delta && bound.angle <= theta;
                        out.push((tri.map(|(n, _)| n), bound));
                    }
                }
            }
            if !ok {
                h /= 1.5;
                continue;
            }
            nodes.extend(local);
            if nodes.len() > MAX_NODES {
                return Err(FaceError::Fatal(Error::ComputationLimit(
                    "tessellation size",
                )));
            }
            let flip = self.outward < 0.0;
            for (v, bound) in out {
                triangles.push(if flip { [v[0], v[2], v[1]] } else { v });
                bounds_out.push(bound);
            }
            return Ok(());
        }
        Err(FaceError::Fatal(Error::ComputationLimit(
            "tessellation refinement",
        )))
    }
}

// ------------------------------------------------------------------ checks

/// The contract, checked on every result: every triangle's nodes distinct,
/// every segment of an edge's polyline in as many triangles as the edge has
/// fins, every other mesh edge in exactly two, and a mesh edge in two
/// triangles traversed once each way; every bound within the request.
fn check_mesh(topology: &Topology, mesh: &Mesh, parameters: Parameters) -> Result<()> {
    let broken = Error::ComputationLimit("tessellation contract");
    if mesh.deflection > parameters.deflection || mesh.angle > parameters.angle {
        return Err(broken);
    }
    let mut directed: BTreeMap<(usize, usize), u32> = BTreeMap::new();
    for t in &mesh.triangles {
        if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
            return Err(broken);
        }
        for k in 0..3 {
            *directed.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
        }
    }
    let mut expected: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (e, polyline) in mesh.edges.iter().enumerate() {
        let fins = topology.edges()[e].fins.len();
        for w in polyline.nodes.windows(2) {
            if w[0] == w[1] {
                return Err(broken);
            }
            *expected
                .entry((w[0].min(w[1]), w[0].max(w[1])))
                .or_default() += fins;
        }
    }
    let mut uses: BTreeMap<(usize, usize), (u32, u32)> = BTreeMap::new();
    for (&(a, b), &k) in &directed {
        let entry = uses.entry((a.min(b), a.max(b))).or_default();
        if a < b {
            entry.0 += k;
        } else {
            entry.1 += k;
        }
    }
    for (key, (forward, backward)) in &uses {
        let total = (forward + backward) as usize;
        match expected.get(key) {
            Some(want) if *want != total => return Err(broken),
            None if total != 2 => return Err(broken),
            _ => {}
        }
        if total == 2 && (*forward != 1 || *backward != 1) {
            return Err(broken);
        }
    }
    for (key, want) in &expected {
        if *want > 0 && !uses.contains_key(key) {
            return Err(broken);
        }
    }
    Ok(())
}
