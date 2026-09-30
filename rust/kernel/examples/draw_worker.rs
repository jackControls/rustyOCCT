//! Test-only DRAW command adapter. Production code never depends on Tcl or OCCT.
//! One persistent shape table per test process; no expected values live here.
//!
//! History commands delegate to the kernel's `History`. DRAW's `generated`
//! follows OCCT's `BRepPrimAPI_MakePrism::Generated`: a profile vertex gives
//! its swept (vertical or seam) edge, a profile edge its wall, the profile
//! face the solid. Start and end copies are OCCT's First/Last shapes, which
//! DRAW does not query; the kernel reports them as generated with their roles.
use rusty_occt::history::{History, Relation};
use rusty_occt::identity::{InputLabel, OperationId, Parent, Role};
use rusty_occt::topology::{
    Curve2, Curve3, EdgeId, FaceId, Loop, Orientation, Slot, Surface, Topology,
};
use rusty_occt::{
    Boundary, BoundaryLabels, Frame3, Point2, Point3, Profile, RigidTransform, Segment, Solid,
    Tolerance, Vec3,
};
use std::collections::{BTreeMap, BTreeSet};

#[path = "../tests/support/draw_geometry.rs"]
mod draw_geometry;
use std::f64::consts::TAU;
use std::io::{self, BufRead, Write};

enum Failure {
    Error(String),
    Unsupported(String),
}
impl From<rusty_occt::Error> for Failure {
    fn from(value: rusty_occt::Error) -> Self {
        Self::Error(value.to_string())
    }
}
type Result<T> = std::result::Result<T, Failure>;

/// A closed planar polyline and its labels in input order.
#[derive(Clone)]
struct Polyline {
    points: Vec<Point3>,
    labels: BoundaryLabels,
    /// S5: per segment, `None` for a line or its arc (centre, radius,
    /// counter-clockwise about the plane's normal), with the plane the
    /// `profile` command built it in.
    arcs: Vec<Option<(Point3, f64, bool)>>,
    plane: Option<Frame3>,
}

impl Polyline {
    fn has_arcs(&self) -> bool {
        self.arcs.iter().any(Option::is_some)
    }
}

/// A DRAW surface or curve (S6), in its DRAW placement; since S7 also
/// spheres, cones (the half-angle in radians) and tori, and `intersect`'s
/// curves.
#[derive(Clone)]
enum Geom {
    Plane(Frame3),
    Cylinder(Frame3, f64),
    Line(Point3, Vec3),
    Circle(Frame3, f64),
    Sphere(Frame3, f64),
    Cone(Frame3, f64, f64),
    Torus(Frame3, f64, f64),
    Curve(draw_geometry::Curve),
}

impl Geom {
    /// The analytic surface, for `intersect` and `xdistcs`.
    fn surface(&self) -> Option<Surface> {
        Some(match *self {
            Geom::Plane(frame) => Surface::Plane(frame),
            Geom::Cylinder(frame, radius) => Surface::Cylinder { frame, radius },
            Geom::Sphere(frame, radius) => Surface::Sphere { frame, radius },
            Geom::Cone(frame, half_angle, radius) => Surface::Cone {
                frame,
                radius,
                half_angle,
            },
            Geom::Torus(frame, major, minor) => Surface::Torus {
                frame,
                major,
                minor,
            },
            _ => return None,
        })
    }

    /// The curve, for `bounds`, `xdistcs` and `dump`.
    fn curve(&self) -> Option<draw_geometry::Curve> {
        Some(match self {
            Geom::Line(p, d) => draw_geometry::Curve::Line(*p, *d),
            Geom::Circle(f, r) => draw_geometry::Curve::Circle(*f, *r),
            Geom::Curve(c) => c.clone(),
            _ => return None,
        })
    }
}

#[derive(Clone)]
enum Shape {
    Solid(Box<Solid>),
    Wire(Polyline),
    Face(Polyline),
    ProfileVertex {
        label: InputLabel,
    },
    ProfileEdge {
        label: InputLabel,
        start: Point3,
        end: Point3,
    },
    Sub {
        body: Box<BodyRef>,
        slot: Slot,
    },
    /// A solid restored from a `.brep` file by the T2 converter, with the
    /// resolution its OCCT tolerances give it; since S6 also a restored free
    /// shape or a `mkface`/`mkedge` result, `kind` its OCCT type.
    Body {
        body: Box<BodyRef>,
        resolution: Tolerance,
        kind: &'static str,
        /// The kernel's sheet or wire, when the adapter built it (a plane
        /// `mkface`, a `bsplit` piece): what `bsplit` splits (S8e).
        kernel: Option<Box<rusty_occt::Body>>,
    },
    /// A DRAW geometric object (S6): `plane`, `cylinder`, `line`, `circle`.
    Geometry(Geom),
    Compound(Vec<Shape>),
    /// `bsplit`'s result (S8e): every object's pieces in the kernel's order.
    /// `groups` holds each object's first piece, piece count and whether it
    /// was a wire (OCCT's General Fuse keeps a wire one wire of its split
    /// edges); `alias` maps each entity a piece shares with an earlier piece
    /// of its object (on the plane: a cut face, a chord, a crossing) to that
    /// one's key, since OCCT shares them between the pieces.
    Split {
        pieces: Vec<Shape>,
        groups: Vec<(usize, usize, bool)>,
        alias: BTreeMap<String, String>,
    },
    /// A prism `prism` made without `Copy` (S9a): OCCT builds its end face
    /// as the start face moved by a location, reusing its `TShape`s, so DRAW
    /// counts fewer vertices, edges and faces than the kernel's prism has
    /// (and so does a Boolean of it); volumes, areas and lengths per use
    /// agree. Its counts and history are unsupported.
    Uncopied(Box<Solid>),
    /// A Boolean's result (S9a): OCCT's compound of the result's solids,
    /// none for an empty result, each a prism, a stack (S9a.2) or a
    /// polyhedron (S9b.1), the last two general bodies, a closed cavity
    /// their second shell. The solids share nothing (the kernel refuses
    /// results touching themselves). `uncopied` when an argument was an
    /// uncopied prism, whose shared shapes OCCT's result keeps.
    Boolean {
        solids: Vec<Solid>,
        uncopied: bool,
    },
    /// A native pick with no entity in the cell model (a seam, its vertex),
    /// or none the selector could single out.
    Lost(String),
}

/// A body's topology and the tag its entities are counted under: a prism's
/// body id, or a serial for a restored body (imported ids depend only on
/// entity counts, so two restored bodies may share them).
#[derive(Clone)]
struct BodyRef {
    tag: String,
    topology: Topology,
}

impl BodyRef {
    fn of(solid: &Solid) -> Self {
        Self {
            tag: solid.topology().body_id().to_string(),
            topology: solid.topology().clone(),
        }
    }
}

#[derive(Clone)]
struct Saved {
    history: History,
    output: Solid,
    /// The profile face the operation swept.
    face: Polyline,
}

#[derive(Default)]
struct Session {
    shapes: BTreeMap<String, Shape>,
    histories: BTreeMap<String, Saved>,
    last: Option<Saved>,
    next_label: u64,
    next_operation: u64,
    /// `explode` calls so far; the native selector numbers them alike.
    explodes: usize,
    selector: Option<Vec<Pick>>,
    /// Restored bodies so far, for their tags.
    restored: u64,
    /// DRAW's numeric variables (S7: `bounds`, `dval`).
    numbers: BTreeMap<String, f64>,
    /// The General Fuse arguments (S8e): `baddobjects`, `baddtools`, and
    /// the argument lists' generation `bfillds` last filled.
    objects: Vec<Shape>,
    tools: Vec<Shape>,
    arguments: u64,
    filled: Option<u64>,
    /// The last operation was a split, a Boolean or an uncopied prism,
    /// whose history the adapter does not keep for `savehistory`.
    unkept_history: bool,
    /// The arguments `bop` last prepared (S9a): the object and the tool.
    bop: Option<(Shape, Shape)>,
}

fn unsupported(args: &[String]) -> Failure {
    Failure::Unsupported(format!("unsupported DRAW signature: {}", args.join(" ")))
}
/// Tcl's `string match` for `*` and `?`.
fn glob(pattern: &str, name: &str) -> bool {
    fn go(p: &[u8], n: &[u8]) -> bool {
        match (p.first(), n.first()) {
            (None, None) => true,
            (Some(b'*'), _) => go(&p[1..], n) || (!n.is_empty() && go(p, &n[1..])),
            (Some(b'?'), Some(_)) => go(&p[1..], &n[1..]),
            (Some(a), Some(b)) if a == b => go(&p[1..], &n[1..]),
            _ => false,
        }
    }
    go(pattern.as_bytes(), name.as_bytes())
}
fn error(message: &str) -> Failure {
    Failure::Error(message.into())
}
fn numbers(args: &[String]) -> Result<Vec<f64>> {
    args.iter()
        .map(|s| {
            s.parse::<f64>()
                .ok()
                .filter(|n| n.is_finite())
                .ok_or_else(|| Failure::Error(format!("invalid finite number: {s}")))
        })
        .collect()
}
fn get<'a>(shapes: &'a BTreeMap<String, Shape>, name: &str) -> Result<&'a Shape> {
    match shapes.get(name) {
        Some(Shape::Lost(why)) => return Err(Failure::Unsupported(format!("{name}: {why}"))),
        Some(Shape::Geometry(_)) => {
            return Err(Failure::Unsupported(format!(
                "{name}: a geometric object, not a shape"
            )))
        }
        _ => {}
    }
    shapes
        .get(name)
        .ok_or_else(|| Failure::Error(format!("unknown shape: {name}")))
}
fn solid<'a>(
    shapes: &'a BTreeMap<String, Shape>,
    name: &str,
    args: &[String],
) -> Result<&'a Solid> {
    match get(shapes, name)? {
        Shape::Solid(s) => Ok(s),
        _ => Err(unsupported(args)),
    }
}

// ------------------------------------------------------------------ measures

/// Twice the signed area and twice the first moments of a pcurve.
fn moments(p: &Curve2) -> [f64; 3] {
    match p {
        Curve2::BSpline(_) => unreachable!("no spline geometry reaches the adapter before S4"),
        Curve2::LineSegment { start: a, end: b } => {
            let (du, dv) = (b.x - a.x, b.y - a.y);
            [
                a.x * b.y - b.x * a.y,
                dv * (a.x * a.x + a.x * du + du * du / 3.0),
                -du * (a.y * a.y + a.y * dv + dv * dv / 3.0),
            ]
        }
        Curve2::CircularArc {
            center: c,
            radius: r,
            start_angle: t0,
            sweep_angle: sweep,
        } => {
            let t1 = t0 + sweep;
            let (s0, c0, s1, c1) = (t0.sin(), t0.cos(), t1.sin(), t1.cos());
            let sq_c = |t: f64| t / 2.0 + (2.0 * t).sin() / 4.0;
            let sq_s = |t: f64| t / 2.0 - (2.0 * t).sin() / 4.0;
            [
                r * (c.x * (s1 - s0) - c.y * (c1 - c0)) + r * r * sweep,
                r * (c.x * c.x * (s1 - s0)
                    + 2.0 * c.x * r * (sq_c(t1) - sq_c(*t0))
                    + r * r * ((s1 - s1.powi(3) / 3.0) - (s0 - s0.powi(3) / 3.0))),
                r * (-c.y * c.y * (c1 - c0)
                    + 2.0 * c.y * r * (sq_s(t1) - sq_s(*t0))
                    + r * r * ((-c1 + c1.powi(3) / 3.0) - (-c0 + c0.powi(3) / 3.0))),
            ]
        }
        Curve2::EllipseArc {
            center: c,
            major,
            minor,
            start_angle: t0,
            sweep_angle: sweep,
        } => {
            let (big, small) = (*major, *minor);
            let t1 = t0 + sweep;
            let (s0, c0, s1, c1) = (t0.sin(), t0.cos(), t1.sin(), t1.cos());
            let cube_s = |s: f64| s - s * s * s / 3.0;
            let cube_c = |c: f64| -c + c * c * c / 3.0;
            let sq_c = |t: f64| t / 2.0 + (2.0 * t).sin() / 4.0;
            let sq_s = |t: f64| t / 2.0 - (2.0 * t).sin() / 4.0;
            [
                small * c.x * (s1 - s0) - big * c.y * (c1 - c0) + big * small * sweep,
                small
                    * (c.x * c.x * (s1 - s0)
                        + 2.0 * c.x * big * (sq_c(t1) - sq_c(*t0))
                        + big * big * (cube_s(s1) - cube_s(s0))),
                big * (-c.y * c.y * (c1 - c0)
                    + 2.0 * c.y * small * (sq_s(t1) - sq_s(*t0))
                    + small * small * (cube_c(c1) - cube_c(c0))),
            ]
        }
        Curve2::Sinusoid { .. } => unreachable!("a sinusoid lies on a cylinder"),
        // A projection (S8d.2): Gauss nodes on 64 pieces, the derivative by
        // central differences.
        Curve2::Projection(_) => {
            let nodes: [(f64, f64); 3] = [
                (0.5 - 0.5 * 0.6f64.sqrt(), 5.0 / 18.0),
                (0.5, 8.0 / 18.0),
                (0.5 + 0.5 * 0.6f64.sqrt(), 5.0 / 18.0),
            ];
            let pieces = 256;
            let h = 1e-6;
            let mut out = [0.0; 3];
            for k in 0..pieces {
                for (x, w) in nodes {
                    let t = (k as f64 + x) / pieces as f64;
                    let (a, b) = (p.point((t - h).max(0.0)), p.point((t + h).min(1.0)));
                    let span = (t + h).min(1.0) - (t - h).max(0.0);
                    let (du, dv) = ((b.x - a.x) / span, (b.y - a.y) / span);
                    let q = p.point(t);
                    let w = w / pieces as f64;
                    out[0] += w * (q.x * dv - q.y * du);
                    out[1] += w * q.x * q.x * dv;
                    out[2] -= w * q.y * q.y * du;
                }
            }
            out
        }
    }
}

fn face_area(t: &Topology, face: usize) -> f64 {
    // An exact projection's pcurve (an ellipse on a cap, S9c.1): the
    // kernel's certified area.
    let projected = t
        .face_fins(FaceId::new(face))
        .iter()
        .flatten()
        .any(|u| matches!(u.pcurve, Curve2::Projection(_)));
    if projected {
        return t
            .face_area_and_centre(FaceId::new(face))
            .map_or(f64::NAN, |(area, _)| area);
    }
    let fins = t.face_fins(FaceId::new(face));
    let fins = fins.iter().flatten();
    match t.faces()[face].surface {
        // A restored spline face (S4): the kernel's certified area.
        Surface::BSpline(_) => t
            .face_area_and_centre(FaceId::new(face))
            .map_or(f64::NAN, |(area, _)| area),
        Surface::Plane(_) => 0.5 * fins.map(|u| moments(&u.pcurve)[0]).sum::<f64>().abs(),
        // Wall pcurves are lines on the universal cover: the area is the
        // radius times the periodic area -∮ v du, wound loops included.
        Surface::Cylinder { radius, .. } => {
            let periodic = fins
                .map(|u| match u.pcurve {
                    Curve2::BSpline(_) => {
                        unreachable!("no spline geometry reaches the adapter before S4")
                    }
                    Curve2::LineSegment { start: a, end: b } => -0.5 * (a.y + b.y) * (b.x - a.x),
                    Curve2::Sinusoid { start, sweep, a } => {
                        let (u0, u1) = (start, start + sweep);
                        -(a[0] * sweep + a[1] * (u1.sin() - u0.sin())
                            - a[2] * (u1.cos() - u0.cos()))
                    }
                    Curve2::CircularArc { .. }
                    | Curve2::EllipseArc { .. }
                    | Curve2::Projection(_) => f64::NAN,
                })
                .sum::<f64>();
            radius * periodic.abs()
        }
        // Cone and sphere faces are measured by the kernel's general mass
        // properties.
        Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } => t
            .face_area_and_centre(FaceId::new(face))
            .map_or(f64::NAN, |(area, _)| area),
    }
}

/// A whole sphere: a sphere face without edge loops (OCCT's one wire of a
/// seam and two degenerated edges at its pole vertices).
fn whole_sphere(t: &Topology, face: usize) -> bool {
    let face = &t.faces()[face];
    matches!(face.surface, Surface::Sphere { .. })
        && face
            .loops
            .iter()
            .all(|l| matches!(t.loops()[l.index()], Loop::Vertex(_)))
}

/// The seam OCCT's wound faces carry and the kernel does not: its length (the
/// face's extent across the winding, to the apex or pole of a pole) when the
/// face has one. `v` is arc length along a cylinder's or cone's rulings and
/// the angle on a sphere's meridian.
fn seam_length(t: &Topology, face: usize) -> Option<f64> {
    if let (true, Surface::Sphere { radius, .. }) =
        (whole_sphere(t, face), &t.faces()[face].surface)
    {
        return Some(std::f64::consts::PI * radius);
    }
    // A torus wound in v: its seam is the latitude arc at the rings' v over
    // their u range.
    if let Surface::Torus { major, minor, .. } = t.faces()[face].surface {
        let mut u = (f64::MAX, f64::MIN);
        let mut v0 = None;
        for l in &t.faces()[face].loops {
            if let Loop::Edges { fins, winding } = &t.loops()[l.index()] {
                if winding[1] != 0 {
                    let p = t.fins()[fins[0].index()].pcurve.point(0.0);
                    u = (u.0.min(p.x), u.1.max(p.x));
                    v0 = Some(p.y);
                }
            }
        }
        if let Some(v0) = v0 {
            return Some((u.1 - u.0) * (major + minor * v0.cos()));
        }
    }
    let face = &t.faces()[face];
    let turns: i32 = face
        .loops
        .iter()
        .filter_map(|l| match &t.loops()[l.index()] {
            Loop::Edges { winding, .. } => Some(winding[0]),
            Loop::Vertex(_) => None,
        })
        .sum();
    let mut wound = false;
    let mut v = (f64::MAX, f64::MIN);
    for l in &face.loops {
        match &t.loops()[l.index()] {
            Loop::Edges { fins, winding } => {
                wound |= winding != &[0, 0];
                for k in fins {
                    let y = t.fins()[k.index()].pcurve.point(0.0).y;
                    v = (v.0.min(y), v.1.max(y));
                }
            }
            Loop::Vertex(_) => {
                let pole = match face.surface {
                    Surface::Cone {
                        radius, half_angle, ..
                    } => Some(-radius / half_angle.sin()),
                    // The pole on the band's material side.
                    Surface::Sphere { .. } => {
                        Some(if (turns > 0) == (face.sense == Orientation::Forward) {
                            std::f64::consts::FRAC_PI_2
                        } else {
                            -std::f64::consts::FRAC_PI_2
                        })
                    }
                    _ => None,
                };
                if let Some(p) = pole {
                    v = (v.0.min(p), v.1.max(p));
                }
            }
        }
    }
    let scale = match face.surface {
        Surface::Sphere { radius, .. } => radius,
        Surface::Torus { minor, .. } => minor,
        _ => 1.0,
    };
    wound.then_some(scale * (v.1 - v.0))
}

fn edge_length(curve: &Curve3) -> f64 {
    match curve {
        Curve3::BSpline(_) => unreachable!("no spline geometry reaches the adapter before S4"),
        Curve3::LineSegment { start, end } => start.distance(*end),
        Curve3::Circle { radius, .. } => TAU * radius,
        Curve3::CircularArc {
            radius,
            sweep_angle,
            ..
        } => radius * sweep_angle.abs(),
        // An elliptic integral: Gauss nodes on 64 pieces of the angle.
        Curve3::EllipseArc {
            major,
            minor,
            start_angle,
            sweep_angle,
            ..
        } => {
            let nodes = gauss();
            let pieces = 64;
            let mut total = 0.0;
            for k in 0..pieces {
                for (x, w) in &nodes {
                    let a = start_angle + sweep_angle * (k as f64 + x) / pieces as f64;
                    total += w * (major * a.sin()).hypot(minor * a.cos());
                }
            }
            total * sweep_angle.abs() / pieces as f64
        }
        // Hyperbolas, parabolas, torus sections (S8d.2-3), cylinders'
        // meetings (S9c.2) and a torus's with a quadric (S9d.4b.2): chords
        // of 4096 pieces.
        Curve3::HyperbolaArc { .. }
        | Curve3::ParabolaArc { .. }
        | Curve3::Section(_)
        | Curve3::Meet(_)
        | Curve3::Rise(_)
        | Curve3::Toric(_) => {
            let pieces = 4096;
            (0..pieces)
                .map(|k| {
                    curve
                        .point(k as f64 / pieces as f64)
                        .distance(curve.point((k + 1) as f64 / pieces as f64))
                })
                .sum()
        }
    }
}

/// Distinct vertices, edges, faces and solids of a shape (by body and slot
/// for kernel shapes, by label for profile shapes).
#[derive(Default)]
struct Parts {
    vertices: BTreeSet<String>,
    edges: BTreeSet<String>,
    wires: usize,
    faces: BTreeSet<String>,
    shells: usize,
    solids: usize,
    compounds: usize,
    length: f64,
    area: f64,
    volume: f64,
    /// A `bsplit` piece's prefix, and the entities pieces share (S8e).
    scope: String,
    alias: BTreeMap<String, String>,
}

impl Parts {
    /// An entity's key: its body's tag and its id, within the current
    /// piece's scope; a shared entity's is the first piece's.
    fn key(&self, b: &BodyRef, slot: Slot) -> String {
        let own = format!(
            "{}{}:{}",
            self.scope,
            b.tag,
            b.topology.id_of(slot).unwrap()
        );
        self.alias.get(&own).cloned().unwrap_or(own)
    }
}

fn collect(shape: &Shape, parts: &mut Parts) {
    // Lengths sum per occurrence, as OCCT's lprops explores edges; counts are
    // of distinct shapes, as nbshapes reports them.
    let add_edge = |s: &BodyRef, e: usize, parts: &mut Parts| {
        let edge = &s.topology.edges()[e];
        let key = parts.key(s, Slot::Edge(EdgeId::new(e)));
        parts.length += edge_length(&edge.curve);
        for v in [edge.start, edge.end].into_iter().flatten() {
            let vertex = parts.key(s, Slot::Vertex(v));
            parts.vertices.insert(vertex);
        }
        // OCCT splits a ring edge at a seam vertex.
        if edge.is_ring() {
            parts.vertices.insert(format!("{key}:seam"));
        }
        parts.edges.insert(key);
    };
    let add_face = |s: &BodyRef, f: usize, parts: &mut Parts| {
        let t = &s.topology;
        let face = parts.key(s, Slot::Face(FaceId::new(f)));
        let first = parts.faces.insert(face.clone());
        // A face two split pieces share (S8e) is counted once and measured
        // per use, as OCCT's explorer visits it in each solid.
        let shared = parts.alias.values().any(|k| *k == face);
        if !first && !shared {
            return;
        }
        parts.area += face_area(t, f);
        for fins in t.face_fins(FaceId::new(f)) {
            for u in fins {
                add_edge(s, u.edge.index(), parts);
            }
        }
        if !first {
            return;
        }
        let unwound = t.faces()[f]
            .loops
            .iter()
            .filter(|l| {
                matches!(
                    &t.loops()[l.index()],
                    Loop::Edges {
                        winding: [0, 0],
                        ..
                    }
                )
            })
            .count();
        parts.wires += unwound;
        // A cone's or sphere's pole is OCCT's degenerated edge at the pole
        // vertex, of no length, in the wound face's wire; so is each pass of
        // a loop through a pole, and a whole sphere has two, at vertices of
        // their own.
        for l in &t.faces()[f].loops {
            if let Loop::Vertex(v) = &t.loops()[l.index()] {
                let vertex = parts.key(s, Slot::Vertex(*v));
                parts.vertices.insert(vertex);
                parts.edges.insert(format!("{face}:pole"));
            }
        }
        for k in 0..t.pole_passes(&t.faces()[f]) {
            parts.edges.insert(format!("{face}:pass{k}"));
        }
        if whole_sphere(t, f) {
            for end in ["south", "north"] {
                parts.vertices.insert(format!("{face}:{end}"));
                parts.edges.insert(format!("{face}:{end}"));
            }
        }
        // A whole torus: one vertex and two seams (the meridian and the
        // latitude circles), each used twice in its one wire.
        if let (Surface::Torus { major, minor, .. }, true) = (
            &t.faces()[f].surface,
            t.faces()[f]
                .loops
                .iter()
                .all(|l| matches!(t.loops()[l.index()], Loop::Vertex(_))),
        ) {
            parts.vertices.insert(format!("{face}:corner"));
            parts.edges.insert(format!("{face}:u-seam"));
            parts.edges.insert(format!("{face}:v-seam"));
            parts.wires += 1;
            parts.length += 2.0 * TAU * minor + 2.0 * TAU * (major + minor);
        }
        // A wound face's two ring loops are one OCCT wire, joined by a seam
        // used in both directions.
        if let Some(length) = seam_length(t, f) {
            parts.wires += 1;
            parts.edges.insert(format!("{face}:seam"));
            parts.length += 2.0 * length;
        }
    };
    let polyline = |p: &Polyline, parts: &mut Parts, face: bool| {
        let n = p.points.len();
        let boundary = if p.has_arcs() {
            boundary_of(p, Tolerance::default()).ok().map(|(_, b)| b)
        } else {
            None
        };
        for k in 0..n {
            parts
                .vertices
                .insert(format!("L{}", p.labels.vertices[k].0));
            parts.edges.insert(format!("L{}", p.labels.segments[k].0));
            let chord = p.points[k].distance(p.points[(k + 1) % n]);
            parts.length += match p.arcs[k] {
                None => chord,
                // The arc's length from its sweep (S5).
                Some((center, radius, ccw)) => {
                    let frame = p.plane.expect("arcs come from a planar profile");
                    let local = |q: Point3| {
                        let [x, y, _] = frame.coordinates(q);
                        Point2::new(
                            x - frame.coordinates(center)[0],
                            y - frame.coordinates(center)[1],
                        )
                    };
                    let (a, b) = (local(p.points[k]), local(p.points[(k + 1) % n]));
                    let mut turn = (a.x * b.y - a.y * b.x).atan2(a.x * b.x + a.y * b.y);
                    if turn <= 0.0 {
                        turn += TAU;
                    }
                    radius * if ccw { turn } else { TAU - turn }
                }
            };
        }
        parts.wires += 1;
        if face && parts.faces.insert(format!("L{}", p.labels.boundary.0)) {
            parts.area += boundary.map_or_else(|| polygon_area(&p.points), |b| b.area());
        }
    };
    let add_solid = |s: &Solid, parts: &mut Parts| {
        parts.solids += 1;
        parts.shells += s.topology().occt_counts().shells;
        parts.volume += s.mass_properties().volume;
        let b = BodyRef::of(s);
        for f in 0..b.topology.faces().len() {
            add_face(&b, f, parts);
        }
    };
    match shape {
        Shape::Solid(s) | Shape::Uncopied(s) => add_solid(s, parts),
        // A Boolean's solids share nothing: each is counted in its own
        // scope, since two may carry equal ids (a fuse of disjoint copies).
        Shape::Boolean { solids, .. } => {
            parts.compounds += 1;
            let scope = parts.scope.clone();
            for (k, s) in solids.iter().enumerate() {
                parts.scope = format!("{scope}B{k}/");
                add_solid(s, parts);
            }
            parts.scope = scope;
        }
        Shape::Body { body, .. } => {
            let t = &body.topology;
            let c = t.occt_counts();
            parts.solids += c.solids;
            parts.shells += c.shells;
            // General certified mass properties (REVIEW_NOTES S3, U2).
            if c.solids > 0 {
                parts.volume += t
                    .mass_enclosure()
                    .map_or(f64::NAN, |m| m.midpoints().volume);
            }
            for f in 0..t.faces().len() {
                add_face(body, f, parts);
            }
            // S6: a wire's edges (one wire when several) and an acorn's
            // vertex.
            for shell in t.shells() {
                for e in &shell.wire_edges {
                    add_edge(body, e.index(), parts);
                }
                parts.wires += usize::from(shell.wire_edges.len() > 1);
                for v in &shell.acorn_vertices {
                    let key = parts.key(body, Slot::Vertex(*v));
                    parts.vertices.insert(key);
                }
            }
        }
        Shape::Geometry(_) => unreachable!("geometry is not a shape"),
        Shape::Wire(p) => polyline(p, parts, false),
        Shape::Face(p) => polyline(p, parts, true),
        Shape::ProfileVertex { label, .. } => {
            parts.vertices.insert(format!("L{}", label.0));
        }
        Shape::ProfileEdge { label, start, end } => {
            parts.edges.insert(format!("L{}", label.0));
            parts.length += start.distance(*end);
        }
        Shape::Sub { body, slot } => match slot {
            Slot::Vertex(v) => {
                let vertex = parts.key(body, Slot::Vertex(*v));
                parts.vertices.insert(vertex);
            }
            Slot::Edge(e) => add_edge(body, e.index(), parts),
            Slot::Face(f) => add_face(body, f.index(), parts),
            Slot::Region(_) => {}
        },
        Shape::Compound(items) => {
            parts.compounds += 1;
            for item in items {
                collect(item, parts);
            }
        }
        // OCCT's splitter returns a compound of every object's images,
        // or the one image itself; a wire's image is one wire.
        Shape::Split {
            pieces,
            groups,
            alias,
        } => {
            parts.compounds += usize::from(split_items(groups) > 1);
            let scope = parts.scope.clone();
            parts.alias.extend(alias.clone());
            for &(first, count, wire) in groups {
                let wires = parts.wires;
                for (k, piece) in pieces.iter().enumerate().skip(first).take(count) {
                    parts.scope = format!("{scope}P{k}/");
                    collect(piece, parts);
                }
                if wire {
                    parts.wires = wires + 1;
                }
            }
            parts.scope = scope;
        }
        Shape::Lost(_) => unreachable!("lost picks are never read"),
    }
}

/// The shapes OCCT's splitter result holds: each piece of a solid or a
/// sheet, one wire per split wire.
fn split_items(groups: &[(usize, usize, bool)]) -> usize {
    groups
        .iter()
        .map(|&(_, count, wire)| if wire { 1 } else { count })
        .sum()
}

fn parts(shape: &Shape) -> Parts {
    let mut p = Parts::default();
    collect(shape, &mut p);
    p
}

/// Planar polygon area by Newell's vector.
fn polygon_area(points: &[Point3]) -> f64 {
    let n = points.len();
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for k in 0..n {
        let (a, b) = (
            points[k] - Point3::ORIGIN,
            points[(k + 1) % n] - Point3::ORIGIN,
        );
        sum = sum + a.cross(b);
    }
    0.5 * sum.length()
}

fn type_name(shape: &Shape) -> &'static str {
    match shape {
        Shape::Solid(_) | Shape::Uncopied(_) => "SOLID",
        Shape::Boolean { .. } => "COMPOUND",
        Shape::Body { kind, .. } => kind,
        Shape::Geometry(_) => unreachable!("geometry is not a shape"),
        Shape::Wire(_) => "WIRE",
        Shape::Face(_) => "FACE",
        Shape::ProfileVertex { .. } => "VERTEX",
        Shape::ProfileEdge { .. } => "EDGE",
        Shape::Sub { slot, .. } => match slot {
            Slot::Vertex(_) => "VERTEX",
            Slot::Edge(_) => "EDGE",
            Slot::Face(_) => "FACE",
            Slot::Region(_) => "SOLID",
        },
        Shape::Compound(_) => "COMPOUND",
        Shape::Split { pieces, groups, .. } => match groups.as_slice() {
            [(_, _, true)] => "WIRE",
            [(first, 1, false)] => type_name(&pieces[*first]),
            _ => "COMPOUND",
        },
        Shape::Lost(_) => unreachable!("lost picks are never read"),
    }
}

/// Only the mass is reported: checkprops reads nothing else, and a centre of
/// gravity is not computed for these shapes.
fn props(mass: f64) -> String {
    format!("Mass : {mass:.17e}\n")
}

// ------------------------------------------------------------------ profiles

/// The plane of a closed polyline: Newell normal, origin at the first point,
/// x toward the second.
fn plane_of(p: &Polyline, tolerance: Tolerance) -> Result<(Frame3, Vec<Point2>)> {
    if let Some(frame) = p.plane {
        let local = p
            .points
            .iter()
            .map(|q| {
                let [x, y, _] = frame.coordinates(*q);
                Point2::new(x, y)
            })
            .collect();
        return Ok((frame, local));
    }
    let n = p.points.len();
    let mut normal = Vec3::new(0.0, 0.0, 0.0);
    for k in 0..n {
        let (a, b) = (
            p.points[k] - Point3::ORIGIN,
            p.points[(k + 1) % n] - Point3::ORIGIN,
        );
        normal = normal + a.cross(b);
    }
    let frame = Frame3::new(p.points[0], normal, p.points[1] - p.points[0], tolerance)?;
    let mut local = Vec::new();
    for point in &p.points {
        let [x, y, z] = frame.coordinates(*point);
        if z.abs() > tolerance.linear() {
            return Err(error("wire is not planar"));
        }
        local.push(Point2::new(x, y));
    }
    Ok((frame, local))
}

/// The kernel boundary of a planar wire or face: a polygon, or a path when
/// it has arcs (S5), in its plane's frame, with its labels.
fn boundary_of(p: &Polyline, tolerance: Tolerance) -> Result<(Frame3, Boundary)> {
    let (frame, local) = plane_of(p, tolerance)?;
    // A whole circle (`profile ... C r 360`, S9a): one arc from its point
    // back to it, the kernel's circle.
    if let ([_], [Some((center, radius, _))]) = (local.as_slice(), p.arcs.as_slice()) {
        let [x, y, _] = frame.coordinates(*center);
        let circle = Boundary::circle(Point2::new(x, y), *radius, tolerance)?;
        return Ok((frame, circle.with_labels(p.labels.clone())?));
    }
    let boundary = if p.has_arcs() {
        let segments = p
            .arcs
            .iter()
            .map(|arc| match arc {
                None => Segment::Line,
                Some((center, radius, ccw)) => {
                    let [x, y, _] = frame.coordinates(*center);
                    Segment::Arc {
                        center: Point2::new(x, y),
                        radius: *radius,
                        ccw: *ccw,
                    }
                }
            })
            .collect();
        Boundary::path(local, segments, tolerance)?
    } else {
        Boundary::polygon(local, tolerance)?
    };
    Ok((frame, boundary.with_labels(p.labels.clone())?))
}

fn history_output(saved: &Saved, parent: InputLabel, roles: &[Role]) -> Vec<Shape> {
    let t = saved.output.topology();
    saved
        .history
        .relations
        .iter()
        .filter_map(|r| match r {
            Relation::Generated { from, to, role }
                if from == &[Parent::Label(parent)] && roles.contains(role) =>
            {
                Some(Shape::Sub {
                    body: Box::new(BodyRef::of(&saved.output)),
                    slot: t.slot_of(*to).unwrap(),
                })
            }
            _ => None,
        })
        .collect()
}

// ------------------------------------------------------------------ native selector
//
// OCCT's exploration order is its own. For `explode` of a kernel solid the
// runner first runs the test in native DRAW, which records every pick's
// measure (sprops/lprops mass, printed to 6 significant digits) and centre of
// gravity (Draw variables, full precision). Each pick selects the one kernel
// entity with the same kind, measure and centre; a pick with no such entity
// (a seam) or several is lost, and using it is unsupported.

#[derive(Clone)]
struct Pick {
    explode: usize,
    kind: String,
    mass: f64,
    centre: Point3,
}

fn load_selector() -> Result<Vec<Pick>> {
    let Ok(path) = std::env::var("RUSTY_DRAW_SELECTOR") else {
        return Ok(Vec::new());
    };
    let text = std::fs::read_to_string(path).map_err(|e| error(&format!("selector: {e}")))?;
    let mut picks = Vec::new();
    for line in text.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        let (Some(&"pick"), Some(explode)) = (w.first(), w.get(1)) else {
            return Err(error("malformed selector line"));
        };
        let explode = explode
            .parse()
            .map_err(|_| error("malformed selector line"))?;
        let n = if w.len() == 8 {
            numbers(&w[4..].iter().map(|x| x.to_string()).collect::<Vec<_>>())?
        } else {
            vec![f64::NAN; 4]
        };
        picks.push(Pick {
            explode,
            kind: w.get(3).unwrap_or(&"other").to_string(),
            mass: n[0],
            centre: Point3::new(n[1], n[2], n[3]),
        });
    }
    Ok(picks)
}

/// Gauss-Legendre nodes and weights on [0, 1].
fn gauss() -> Vec<(f64, f64)> {
    const N: usize = 24;
    (1..=N)
        .map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 - 0.25) / (N as f64 + 0.5)).cos();
            let mut d = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=N {
                    let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                    p0 = p1;
                    p1 = p2;
                }
                d = N as f64 * (x * p1 - p0) / (x * x - 1.0);
                let dx = p1 / d;
                x -= dx;
                if dx.abs() < 1e-16 {
                    break;
                }
            }
            ((1.0 - x) / 2.0, 1.0 / ((1.0 - x * x) * d * d))
        })
        .collect()
}

/// A pcurve's point and derivative at `t` in [0, 1].
fn pcurve_at(p: &Curve2, t: f64) -> (Point2, Point2) {
    match p {
        Curve2::BSpline(_) => unreachable!("no spline geometry reaches the adapter before S4"),
        Curve2::LineSegment { start: a, end: b } => (p.point(t), Point2::new(b.x - a.x, b.y - a.y)),
        Curve2::CircularArc {
            radius,
            start_angle,
            sweep_angle,
            ..
        } => {
            let (sin, cos) = (start_angle + sweep_angle * t).sin_cos();
            let k = radius * sweep_angle;
            (p.point(t), Point2::new(-k * sin, k * cos))
        }
        Curve2::EllipseArc {
            major,
            minor,
            start_angle,
            sweep_angle,
            ..
        } => {
            let (sin, cos) = (start_angle + sweep_angle * t).sin_cos();
            (
                p.point(t),
                Point2::new(-major * sweep_angle * sin, minor * sweep_angle * cos),
            )
        }
        Curve2::Sinusoid { start, sweep, a } => {
            let (sin, cos) = (start + sweep * t).sin_cos();
            (
                p.point(t),
                Point2::new(*sweep, sweep * (a[2] * cos - a[1] * sin)),
            )
        }
        Curve2::Projection(_) => {
            let (lo, hi) = ((t - 1e-6).max(0.0), (t + 1e-6).min(1.0));
            let (a, b) = (p.point(lo), p.point(hi));
            (
                p.point(t),
                Point2::new((b.x - a.x) / (hi - lo), (b.y - a.y) / (hi - lo)),
            )
        }
    }
}

/// Area and centre of gravity of a face: -∮ G du over its pcurves for
/// ∂G/∂v = 1, u, v (planes) or 1, cos u, sin u, v (cylinders). Seams are
/// vertical, so wound faces need no closing sides.
fn face_centre(t: &Topology, face: usize) -> (f64, Point3) {
    let nodes = gauss();
    let surface = &t.faces()[face].surface;
    let projected = t
        .face_fins(FaceId::new(face))
        .iter()
        .flatten()
        .any(|u| matches!(u.pcurve, Curve2::Projection(_)));
    // A restored spline face (S4), or one with an exact projection's
    // pcurve (S9c.1): the kernel's certified area and centre.
    if matches!(surface, Surface::BSpline(_)) || projected {
        return t
            .face_area_and_centre(FaceId::new(face))
            .unwrap_or((f64::NAN, Point3::ORIGIN));
    }
    let mut m = [0.0; 4];
    for fin in t.face_fins(FaceId::new(face)).iter().flatten() {
        for (x, w) in &nodes {
            let (q, d) = pcurve_at(&fin.pcurve, *x);
            let g = match surface {
                Surface::BSpline(_) => {
                    unreachable!("no spline geometry reaches the adapter before S4")
                }
                Surface::Plane(_) => [q.y, q.x * q.y, q.y * q.y / 2.0, 0.0],
                Surface::Cylinder { .. } => {
                    [q.y, q.y * q.x.cos(), q.y * q.x.sin(), q.y * q.y / 2.0]
                }
                Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } => [0.0; 4],
            };
            for k in 0..4 {
                m[k] -= w * g[k] * d.x;
            }
        }
    }
    match surface {
        Surface::BSpline(_) => unreachable!("no spline geometry reaches the adapter before S4"),
        Surface::Plane(_) => (
            m[0].abs(),
            surface.point(Point2::new(m[1] / m[0], m[2] / m[0])),
        ),
        Surface::Cylinder { frame, radius } => (
            radius * m[0].abs(),
            frame.point(
                Point2::new(radius * m[1] / m[0], radius * m[2] / m[0]),
                m[3] / m[0],
            ),
        ),
        Surface::Cone { .. } | Surface::Sphere { .. } | Surface::Torus { .. } => t
            .face_area_and_centre(FaceId::new(face))
            .unwrap_or((f64::NAN, Point3::new(f64::NAN, f64::NAN, f64::NAN))),
    }
}

/// Length and centre of gravity of an edge (curves are uniform in their
/// fraction).
fn edge_centre(curve: &Curve3) -> (f64, Point3) {
    let mut c = Vec3::new(0.0, 0.0, 0.0);
    for (x, w) in gauss() {
        c = c + (curve.point(x) - Point3::ORIGIN) * w;
    }
    (edge_length(curve), Point3::ORIGIN + c)
}

fn same_pick(pick: &Pick, mass: f64, centre: Point3) -> bool {
    let scale = 1.0 + pick.centre.distance(Point3::ORIGIN);
    (pick.mass - mass).abs() <= 1e-5 * pick.mass.abs().max(1.0)
        && pick.centre.distance(centre) <= 1e-7 * scale
}

fn select(body: &BodyRef, pick: &Pick) -> Shape {
    let t = &body.topology;
    let found: Vec<Slot> = match pick.kind.as_str() {
        "face" => (0..t.faces().len())
            .filter(|f| {
                let (m, c) = face_centre(t, *f);
                same_pick(pick, m, c)
            })
            .map(|f| Slot::Face(FaceId::new(f)))
            .collect(),
        "edge" => (0..t.edges().len())
            .filter(|e| {
                let (m, c) = edge_centre(&t.edges()[*e].curve);
                same_pick(pick, m, c)
            })
            .map(|e| Slot::Edge(EdgeId::new(e)))
            .collect(),
        _ => Vec::new(),
    };
    match found.as_slice() {
        [slot] => Shape::Sub {
            body: Box::new(body.clone()),
            slot: *slot,
        },
        [] => Shape::Lost(format!(
            "native {} pick has no entity in the cell model",
            pick.kind
        )),
        _ => Shape::Lost(format!(
            "native {} pick matches several entities",
            pick.kind
        )),
    }
}

// ------------------------------------------------------------------ restore

/// A `.brep` file through the T2 reader and converter: its compounds kept
/// as compounds, each solid a body. Constructs the kernel cannot represent
/// make the whole restore unsupported, by name; a solid the converter
/// builds but the validator rejects is an error.
fn restore(path: &str, serial: &mut u64) -> Result<Shape> {
    use rusty_occt::occt_brep::{import, read, read::Kind, Rejected};
    // OCCT reads nothing after the root reference, and some upstream files
    // carry stray bytes there; decode leniently.
    let bytes = std::fs::read(path).map_err(|e| error(&format!("restore: {e}")))?;
    let text = String::from_utf8_lossy(&bytes);
    // DRAW's restore reads any saved Draw object; only shapes are converted.
    let kind = text
        .lines()
        .map(str::trim)
        .find(|l| !l.is_empty())
        .unwrap_or("");
    if kind != "DBRep_DrawableShape" && !kind.starts_with("CASCADE Topology") {
        return Err(Failure::Unsupported(format!("restore: a {kind} object")));
    }
    let doc = read(&text).map_err(|e| error(&format!("restore: {e}")))?;
    let imported = import(&doc);
    if !imported.unsupported.is_empty() {
        let names: Vec<String> = imported
            .unsupported
            .iter()
            .map(|(k, n)| format!("{k}={n}"))
            .collect();
        return Err(Failure::Unsupported(format!(
            "restore: unsupported constructs {}",
            names.join(",")
        )));
    }
    let mut solids = imported.solids.into_iter();
    let mut free = imported.free.into_iter();
    fn build(
        doc: &rusty_occt::occt_brep::Document,
        sub: rusty_occt::occt_brep::read::Sub,
        solids: &mut std::vec::IntoIter<rusty_occt::occt_brep::ImportedSolid>,
        free: &mut std::vec::IntoIter<rusty_occt::occt_brep::ImportedFree>,
        serial: &mut u64,
    ) -> Result<Shape> {
        let shape = &doc.shapes[sub.shape];
        match shape.kind {
            Kind::Compound => Ok(Shape::Compound(
                shape
                    .subs
                    .iter()
                    .map(|s| build(doc, *s, solids, free, serial))
                    .collect::<Result<_>>()?,
            )),
            Kind::Solid => {
                let solid = solids.next().ok_or_else(|| error("restore: solid order"))?;
                match solid.result {
                    Ok(topology) => {
                        *serial += 1;
                        Ok(Shape::Body {
                            body: Box::new(BodyRef {
                                tag: format!("R{serial}"),
                                topology,
                            }),
                            resolution: solid.tolerance,
                            kind: "SOLID",
                            kernel: None,
                        })
                    }
                    Err(Rejected::Invalid { issues, .. }) => Err(error(&format!(
                        "restore: the validator rejects solid record {}: {}",
                        solid.record,
                        issues.first().map(|i| i.to_string()).unwrap_or_default()
                    ))),
                    Err(Rejected::Unsupported(names)) => Err(Failure::Unsupported(format!(
                        "restore: unsupported constructs {}",
                        names.join(",")
                    ))),
                }
            }
            // S6: a free shell, face, wire, edge or vertex.
            Kind::Shell | Kind::Face | Kind::Wire | Kind::Edge | Kind::Vertex => {
                let shape = free
                    .next()
                    .ok_or_else(|| error("restore: free shape order"))?;
                let kind = match shape.kind {
                    Kind::Shell => "SHELL",
                    Kind::Face => "FACE",
                    Kind::Wire => "WIRE",
                    Kind::Edge => "EDGE",
                    _ => "VERTEX",
                };
                match shape.result {
                    Ok(topology) => {
                        *serial += 1;
                        Ok(Shape::Body {
                            body: Box::new(BodyRef {
                                tag: format!("R{serial}"),
                                topology,
                            }),
                            resolution: shape.tolerance,
                            kind,
                            kernel: None,
                        })
                    }
                    Err(Rejected::Invalid { issues, .. }) => Err(error(&format!(
                        "restore: the validator rejects free shape record {}: {}",
                        shape.record,
                        issues.first().map(|i| i.to_string()).unwrap_or_default()
                    ))),
                    Err(Rejected::Unsupported(names)) => Err(Failure::Unsupported(format!(
                        "restore: unsupported constructs {}",
                        names.join(",")
                    ))),
                }
            }
            _ => Err(Failure::Unsupported("restore: a compsolid".into())),
        }
    }
    build(&doc, doc.root, &mut solids, &mut free, serial)
}

// ------------------------------------------------------------------ geometry (S6)

/// A DRAW placement: origin, main direction (default Z) and X direction,
/// chosen as gp_Ax2 chooses it when only the main direction is given.
fn placement(n: &[f64], t: Tolerance) -> Result<Frame3> {
    let origin = Point3::new(n[0], n[1], n[2]);
    let normal = if n.len() >= 6 {
        Vec3::new(n[3], n[4], n[5])
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let x = if n.len() >= 9 {
        Vec3::new(n[6], n[7], n[8])
    } else {
        let (a, b, c) = (normal.x, normal.y, normal.z);
        let (aa, ba, ca) = (a.abs(), b.abs(), c.abs());
        if ba <= aa && ba <= ca {
            if aa > ca {
                Vec3::new(-c, 0.0, a)
            } else {
                Vec3::new(c, 0.0, -a)
            }
        } else if aa <= ba && aa <= ca {
            if ba > ca {
                Vec3::new(0.0, -c, b)
            } else {
                Vec3::new(0.0, c, -b)
            }
        } else if aa > ba {
            Vec3::new(-b, a, 0.0)
        } else {
            Vec3::new(b, -a, 0.0)
        }
    };
    Ok(Frame3::new(origin, normal, x, t)?)
}

fn geometry<'a>(shapes: &'a BTreeMap<String, Shape>, name: &str) -> Result<&'a Geom> {
    match shapes.get(name) {
        Some(Shape::Geometry(g)) => Ok(g),
        _ => Err(Failure::Unsupported(format!(
            "{name}: not a plane, cylinder, line or circle"
        ))),
    }
}

/// A body of parts: validated against the resolution, measured enclosures.
fn body_of(
    parts: rusty_occt::topology::TopologyParts,
    t: Tolerance,
    serial: &mut u64,
    kind: &'static str,
) -> Result<Shape> {
    let topology = Topology::from_parts(parts.with_measured_enclosures(), t).map_err(|issues| {
        error(&format!(
            "invalid construction: {}",
            issues.first().map(|i| i.to_string()).unwrap_or_default()
        ))
    })?;
    *serial += 1;
    Ok(Shape::Body {
        body: Box::new(BodyRef {
            tag: format!("M{serial}"),
            topology,
        }),
        resolution: t,
        kind,
        kernel: None,
    })
}

/// `mkface` on a cylinder: the patch `[u0, u1] x [v0, v1]` of its cover,
/// both sides in the void; one turn is a band bounded by two ring edges.
fn cylinder_patch(
    frame: Frame3,
    radius: f64,
    b: [f64; 4],
    t: Tolerance,
) -> Result<rusty_occt::topology::TopologyParts> {
    use rusty_occt::topology::{
        Edge, Face, Fin, FinId, Loop, LoopId, Region, RegionId, RegionKind, Shell, ShellId, Side,
        TopologyParts, Vertex, VertexId,
    };
    let [u0, u1, v0, v1] = b;
    if !(u0 < u1 && v0 < v1) || u1 - u0 > TAU * (1.0 + 1e-12) {
        return Err(error("mkface: empty or over-wound bounds"));
    }
    let at = |v: f64| {
        Frame3::new(
            frame.point(Point2::default(), v),
            frame.normal(),
            frame.x(),
            t,
        )
    };
    let point = |u: f64, v: f64| {
        at(v).map(|f| f.point(Point2::new(radius * u.cos(), radius * u.sin()), 0.0))
    };
    let mut p = TopologyParts::default();
    let fin = |edge: usize, forward: bool, a: (f64, f64), c: (f64, f64)| Fin {
        edge: EdgeId::new(edge),
        sense: if forward {
            Orientation::Forward
        } else {
            Orientation::Reversed
        },
        pcurve: Curve2::LineSegment {
            start: Point2::new(a.0, a.1),
            end: Point2::new(c.0, c.1),
        },
        enclosure: None,
    };
    let band = (u1 - u0 - TAU).abs() <= 1e-12 * TAU;
    let loops: Vec<(Vec<Fin>, [i32; 2])> = if band {
        for v in [v0, v1] {
            p.edges.push(Edge {
                start: None,
                end: None,
                curve: Curve3::Circle {
                    frame: at(v)?,
                    radius,
                },
                fins: Vec::new(),
            });
        }
        vec![
            (vec![fin(0, true, (u0, v0), (u0 + TAU, v0))], [1, 0]),
            (vec![fin(1, false, (u0 + TAU, v1), (u0, v1))], [-1, 0]),
        ]
    } else {
        for (u, v) in [(u0, v0), (u1, v0), (u1, v1), (u0, v1)] {
            p.vertices.push(Vertex {
                position: point(u, v)?,
                enclosure: None,
            });
        }
        let arc = |v: f64| -> Result<Curve3> {
            Ok(Curve3::CircularArc {
                frame: at(v)?,
                radius,
                start_angle: u0,
                sweep_angle: u1 - u0,
            })
        };
        let line = |a: usize, b: usize, p: &TopologyParts| Curve3::LineSegment {
            start: p.vertices[a].position,
            end: p.vertices[b].position,
        };
        let ends = [
            (0, 1, arc(v0)?),
            (1, 2, line(1, 2, &p)),
            (3, 2, arc(v1)?),
            (0, 3, line(0, 3, &p)),
        ];
        for (a, b, curve) in ends {
            p.edges.push(Edge {
                start: Some(VertexId::new(a)),
                end: Some(VertexId::new(b)),
                curve,
                fins: Vec::new(),
            });
        }
        vec![(
            vec![
                fin(0, true, (u0, v0), (u1, v0)),
                fin(1, true, (u1, v0), (u1, v1)),
                fin(2, false, (u1, v1), (u0, v1)),
                fin(3, false, (u0, v1), (u0, v0)),
            ],
            [0, 0],
        )]
    };
    let mut face_loops = Vec::new();
    for (fins, winding) in loops {
        let ids = fins
            .into_iter()
            .map(|f| {
                p.edges[f.edge.index()].fins.push(FinId::new(p.fins.len()));
                p.fins.push(f);
                FinId::new(p.fins.len() - 1)
            })
            .collect();
        p.loops.push(Loop::Edges { fins: ids, winding });
        face_loops.push(LoopId::new(p.loops.len() - 1));
    }
    p.faces.push(Face {
        surface: Surface::Cylinder { frame, radius },
        sense: Orientation::Forward,
        loops: face_loops,
        front: ShellId::new(0),
        back: ShellId::new(0),
        enclosure: None,
    });
    p.shells.push(Shell {
        region: RegionId::new(0),
        sides: vec![(FaceId::new(0), Side::Front), (FaceId::new(0), Side::Back)],
        wire_edges: Vec::new(),
        acorn_vertices: Vec::new(),
    });
    p.regions.push(Region {
        kind: RegionKind::Void,
        shells: vec![ShellId::new(0)],
    });
    Ok(p)
}

/// `mkedge` on a line or circle: one wire edge in the void, bounded by the
/// parameters; a whole circle is a ring edge.
fn edge_parts(curve: Curve3, ends: Option<[Point3; 2]>) -> rusty_occt::topology::TopologyParts {
    use rusty_occt::topology::{
        Edge, Region, RegionId, RegionKind, Shell, ShellId, TopologyParts, Vertex, VertexId,
    };
    let mut p = TopologyParts::default();
    let (start, end) = match ends {
        Some([a, b]) => {
            p.vertices.push(Vertex {
                position: a,
                enclosure: None,
            });
            p.vertices.push(Vertex {
                position: b,
                enclosure: None,
            });
            (Some(VertexId::new(0)), Some(VertexId::new(1)))
        }
        None => (None, None),
    };
    p.edges.push(Edge {
        start,
        end,
        curve,
        fins: Vec::new(),
    });
    p.shells.push(Shell {
        region: RegionId::new(0),
        sides: Vec::new(),
        wire_edges: vec![EdgeId::new(0)],
        acorn_vertices: Vec::new(),
    });
    p.regions.push(Region {
        kind: RegionKind::Void,
        shells: vec![ShellId::new(0)],
    });
    p
}

// ------------------------------------------------------------------ splitter (S8e)
//
// OCCT's General Fuse splitter (`bsplit`, `bapisplit`) with one planar tool
// face is the kernel's split by the face's plane where the face reaches
// across every object and nothing else interferes. The adapter checks that
// much before it splits, and reports anything else unsupported: several
// tools, a tool that is not a convex planar face containing the objects'
// shadow on its plane, objects that may touch each other (General Fuse
// intersects them too), and a plane through a vertex or tangent to an arc,
// where OCCT's sharing is not confirmed. Ring edges (OCCT's seams) are
// outside the subset.

/// A split's capability limits are unsupported, not failures.
fn split_failure(e: rusty_occt::Error) -> Failure {
    use rusty_occt::Error as E;
    match e {
        E::OutOfDomain(_)
        | E::ComputationLimit(_)
        | E::Degenerate(_)
        | E::LimitExceeded(_)
        | E::PrecisionLoss => Failure::Unsupported(format!("bsplit: {e}")),
        e => e.into(),
    }
}

fn split_unsupported(why: &str) -> Failure {
    Failure::Unsupported(format!("bsplit: {why}"))
}

/// The tool's plane and its face's profile in the plane's frame: a convex
/// polygon without holes.
fn split_tool(tool: &Shape, t: Tolerance) -> Result<(Frame3, Profile)> {
    let (frame, profile) = match tool {
        Shape::Body {
            kernel: Some(body),
            kind: "FACE",
            ..
        } => (
            body.frame(),
            body.profile()
                .ok_or_else(|| split_unsupported("the tool is not a planar face"))?
                .clone(),
        ),
        Shape::Face(p) => {
            let (frame, boundary) = boundary_of(p, t)?;
            (frame, Profile::new(boundary, vec![], t)?)
        }
        _ => return Err(split_unsupported("the tool is not a planar face")),
    };
    let convex = profile.holes().is_empty()
        && profile.outer().polygon_vertices().is_some_and(|p| {
            let n = p.len();
            let turns: Vec<f64> = (0..n)
                .map(|k| {
                    let (a, b, c) = (p[k], p[(k + 1) % n], p[(k + 2) % n]);
                    (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x)
                })
                .collect();
            turns.iter().all(|z| *z > 0.0) || turns.iter().all(|z| *z < 0.0)
        });
    if !convex {
        return Err(split_unsupported("the tool face is not a convex polygon"));
    }
    Ok((frame, profile))
}

/// An object as the kernel splits it: a prism, or a sheet or closed wire.
enum Object {
    Solid(Solid),
    Body(rusty_occt::Body),
}

impl Object {
    fn of(shape: &Shape, operation: OperationId, t: Tolerance) -> Result<Self> {
        Ok(match shape {
            Shape::Solid(s) if s.profile().is_some() => Object::Solid((**s).clone()),
            Shape::Body {
                kernel: Some(body), ..
            } if body.profile().is_some() || body.boundary().is_some() => {
                Object::Body((**body).clone())
            }
            Shape::Face(p) => {
                let (frame, boundary) = boundary_of(p, t)?;
                Object::Body(
                    rusty_occt::Body::face_from_profile_with(
                        operation,
                        Profile::new(boundary, vec![], t)?,
                        frame,
                    )?
                    .0,
                )
            }
            Shape::Wire(p) => {
                let (frame, boundary) = boundary_of(p, t)?;
                Object::Body(
                    rusty_occt::Body::wire_from_boundary_with(operation, boundary, frame, t)?.0,
                )
            }
            _ => {
                return Err(split_unsupported(
                    "an object other than a prism, a planar face or a closed wire",
                ))
            }
        })
    }

    fn topology(&self) -> &Topology {
        match self {
            Object::Solid(s) => s.topology(),
            Object::Body(b) => b.topology(),
        }
    }
}

/// An axis-aligned box about a topology of lines and arcs, after checking
/// that the plane crosses its edges only in their interiors: no vertex
/// within the resolution of the plane, no arc tangent to it.
fn split_extent(topology: &Topology, plane: Frame3, t: Tolerance) -> Result<[Point3; 2]> {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    let mut add = |p: Point3, r: [f64; 3]| {
        for (i, x) in p.to_array().into_iter().enumerate() {
            lo[i] = lo[i].min(x - r[i]);
            hi[i] = hi[i].max(x + r[i]);
        }
    };
    let m = plane.normal();
    for v in topology.vertices() {
        if plane.coordinates(v.position)[2].abs() <= t.linear() {
            return Err(split_unsupported(
                "the tool's plane passes through a vertex",
            ));
        }
        add(v.position, [0.0; 3]);
    }
    for edge in topology.edges() {
        if edge.is_ring() {
            return Err(split_unsupported(
                "a closed edge (OCCT's seam and its vertex)",
            ));
        }
        match edge.curve {
            Curve3::LineSegment { start, end } => {
                add(start, [0.0; 3]);
                add(end, [0.0; 3]);
            }
            Curve3::CircularArc { frame, radius, .. } => {
                // The circle's height over the plane is d0 + r A cos(u - u0).
                let d0 = plane.coordinates(frame.origin())[2];
                let along = m.dot(frame.normal());
                let a = (1.0 - along * along).max(0.0).sqrt();
                if (d0.abs() - radius * a).abs() <= t.linear() {
                    return Err(split_unsupported("the tool's plane is tangent to an arc"));
                }
                let n = frame.normal().to_array();
                add(
                    frame.origin(),
                    n.map(|c| radius * (1.0 - c * c).max(0.0).sqrt()),
                );
            }
            _ => {
                return Err(split_unsupported(
                    "an object with edges other than line segments and arcs",
                ))
            }
        }
    }
    Ok([
        Point3::new(lo[0], lo[1], lo[2]),
        Point3::new(hi[0], hi[1], hi[2]),
    ])
}

fn corners([lo, hi]: [Point3; 2]) -> impl Iterator<Item = Point3> {
    (0..8).map(move |k| {
        Point3::new(
            if k & 1 == 0 { lo.x } else { hi.x },
            if k & 2 == 0 { lo.y } else { hi.y },
            if k & 4 == 0 { lo.z } else { hi.z },
        )
    })
}

/// A piece's body, as `collect` keys its entities.
fn body_ref(shape: &Shape) -> Option<BodyRef> {
    match shape {
        Shape::Solid(s) => Some(BodyRef::of(s)),
        Shape::Body { body, .. } => Some((**body).clone()),
        _ => None,
    }
}

/// The entities one object's pieces share, by geometry: vertices at one
/// point, edges between the same vertices through one midpoint, faces
/// bounded by the same edges. Each maps to the first piece's key.
fn shared_entities(
    pieces: &[Shape],
    first: usize,
    t: Tolerance,
    alias: &mut BTreeMap<String, String>,
) {
    use rusty_occt::topology::VertexId;
    let refs: Vec<(usize, BodyRef)> = pieces
        .iter()
        .enumerate()
        .filter_map(|(k, p)| body_ref(p).map(|b| (first + k, b)))
        .collect();
    let key = |k: usize, b: &BodyRef, slot: Slot| {
        format!("P{k}/{}:{}", b.tag, b.topology.id_of(slot).unwrap())
    };
    let canonical =
        |alias: &BTreeMap<String, String>, own: String| alias.get(&own).cloned().unwrap_or(own);
    let mut points: Vec<(usize, Point3, String)> = Vec::new();
    for (k, b) in &refs {
        for (i, v) in b.topology.vertices().iter().enumerate() {
            let own = key(*k, b, Slot::Vertex(VertexId::new(i)));
            match points
                .iter()
                .find(|(j, p, _)| j != k && p.distance(v.position) <= t.linear())
            {
                Some((_, _, c)) => {
                    alias.insert(own, c.clone());
                }
                None => points.push((*k, v.position, own)),
            }
        }
    }
    type EdgeSignature = (usize, [Option<String>; 2], Point3, String);
    let mut edges: Vec<EdgeSignature> = Vec::new();
    for (k, b) in &refs {
        for (i, e) in b.topology.edges().iter().enumerate() {
            let own = key(*k, b, Slot::Edge(EdgeId::new(i)));
            let mut ends =
                [e.start, e.end].map(|v| v.map(|v| canonical(alias, key(*k, b, Slot::Vertex(v)))));
            ends.sort();
            let middle = e.curve.point(0.5);
            match edges
                .iter()
                .find(|(j, u, p, _)| j != k && *u == ends && p.distance(middle) <= t.linear())
            {
                Some((_, _, _, c)) => {
                    alias.insert(own, c.clone());
                }
                None => edges.push((*k, ends, middle, own)),
            }
        }
    }
    let mut faces: Vec<(usize, Vec<String>, String)> = Vec::new();
    for (k, b) in &refs {
        for f in 0..b.topology.faces().len() {
            let own = key(*k, b, Slot::Face(FaceId::new(f)));
            let mut bounds: Vec<String> = b
                .topology
                .face_fins(FaceId::new(f))
                .iter()
                .flatten()
                .map(|u| canonical(alias, key(*k, b, Slot::Edge(u.edge))))
                .collect();
            bounds.sort();
            match faces.iter().find(|(j, u, _)| j != k && *u == bounds) {
                Some((_, _, c)) => {
                    alias.insert(own, c.clone());
                }
                None => faces.push((*k, bounds, own)),
            }
        }
    }
}

/// `bsplit` and `bapisplit`: every object split by the one tool's plane.
fn split(session: &mut Session, t: Tolerance) -> Result<Shape> {
    let [tool] = session.tools.as_slice() else {
        return Err(split_unsupported("one planar tool face is supported"));
    };
    if session.objects.is_empty() {
        return Err(split_unsupported("no objects"));
    }
    let (plane, face) = split_tool(tool, t)?;
    let mut objects = Vec::new();
    for shape in &session.objects {
        session.next_operation += 1;
        let object = Object::of(shape, OperationId(session.next_operation), t)?;
        let extent = split_extent(object.topology(), plane, t)?;
        objects.push((object, extent));
    }
    // General Fuse also intersects the objects with each other.
    for (i, (_, [a0, a1])) in objects.iter().enumerate() {
        for (_, [b0, b1]) in &objects[..i] {
            let (a0, a1, b0, b1) = (a0.to_array(), a1.to_array(), b0.to_array(), b1.to_array());
            let apart = (0..3).any(|k| a1[k] + t.linear() < b0[k] || b1[k] + t.linear() < a0[k]);
            if !apart {
                return Err(split_unsupported(
                    "objects that may interfere with each other",
                ));
            }
        }
    }
    let (mut pieces, mut groups, mut alias) = (Vec::new(), Vec::new(), BTreeMap::new());
    for (object, extent) in objects {
        let heights: Vec<f64> = corners(extent).map(|c| plane.coordinates(c)[2]).collect();
        let missed =
            heights.iter().all(|h| *h > t.linear()) || heights.iter().all(|h| *h < -t.linear());
        // Where the plane reaches the object, the tool face must cover the
        // object's whole shadow on it.
        if !missed {
            for c in corners(extent) {
                let [x, y, _] = plane.coordinates(c);
                if face.classify(Point2::new(x, y))? != rusty_occt::Location::Inside {
                    return Err(split_unsupported(
                        "the tool face does not reach across an object",
                    ));
                }
            }
        }
        session.next_operation += 1;
        let operation = OperationId(session.next_operation);
        let first = pieces.len();
        let wire = match &object {
            Object::Solid(s) => {
                let (split, _) = s.split_by_plane(operation, plane).map_err(split_failure)?;
                pieces.extend(split.into_iter().map(|(_, p)| Shape::Solid(Box::new(p))));
                false
            }
            Object::Body(b) => {
                let (split, _) = b.split_by_plane(operation, plane).map_err(split_failure)?;
                let wire = b.topology().faces().is_empty();
                for (_, piece) in split {
                    session.restored += 1;
                    pieces.push(Shape::Body {
                        body: Box::new(BodyRef {
                            tag: format!("S{}", session.restored),
                            topology: piece.topology().clone(),
                        }),
                        resolution: piece.resolution(),
                        kind: if wire { "WIRE" } else { "FACE" },
                        kernel: Some(Box::new(piece)),
                    });
                }
                wire
            }
        };
        shared_entities(&pieces[first..], first, t, &mut alias);
        groups.push((first, pieces.len() - first, wire));
    }
    session.last = None;
    session.unkept_history = true;
    Ok(Shape::Split {
        pieces,
        groups,
        alias,
    })
}

// ------------------------------------------------------------------ Booleans (S9a)
//
// OCCT's Boolean commands (`bfuse`, `bcut`, `bcommon`, `btuc`; `bop` and
// `bopfuse`, `bopcut`, `boptuc`, `bopcommon`; `bbop` and `bapibop` on the
// General Fuse arguments) with one object and one tool, each a prism the
// adapter made, through `Solid::fuse`, `cut` and `common`. The result is
// OCCT's compound of the result's solids. Whatever the kernel refuses
// (frames with different axes, an offset that is not binary64 in them,
// S9a.2's stacks, results touching themselves, other solids) is
// unsupported, never an answer.

#[derive(Clone, Copy)]
enum BooleanOp {
    Common,
    Fuse,
    Cut,
    /// OCCT's CUT21: the tool minus the object.
    Tuc,
}

fn boolean_unsupported(why: &str) -> Failure {
    Failure::Unsupported(format!("Boolean: {why}"))
}

/// A Boolean's capability limits are unsupported, not failures.
fn boolean_failure(e: rusty_occt::Error) -> Failure {
    use rusty_occt::Error as E;
    match e {
        E::OutOfDomain(_)
        | E::ComputationLimit(_)
        | E::Degenerate(_)
        | E::LimitExceeded(_)
        | E::PrecisionLoss => Failure::Unsupported(format!("Boolean: {e}")),
        e => e.into(),
    }
}

/// A Boolean's result solid that is a general body, not a prism: S9a.2's
/// stack or S9b's polyhedron (every result solid without a profile).
fn general(s: &Solid) -> bool {
    s.profile().is_none()
}

/// A Boolean argument: one solid the adapter made (a prism, a cone, a
/// sphere or a torus; the kernel decides what it supports), or a Boolean
/// result of one solid (a prism, a stack or a polyhedron: S9b.2 takes the
/// last two on their stored geometry); and whether OCCT's shape reuses
/// its shapes.
fn boolean_argument(shape: &Shape) -> Result<(&Solid, bool)> {
    match shape {
        Shape::Solid(s) => Ok((s, false)),
        Shape::Uncopied(s) => Ok((s, true)),
        Shape::Boolean { solids, uncopied } if solids.len() == 1 => Ok((&solids[0], *uncopied)),
        Shape::Boolean { .. } => Err(boolean_unsupported("an argument of several solids or none")),
        _ => Err(boolean_unsupported(
            "an argument other than a solid the adapter made",
        )),
    }
}

/// The same prism extruded again under an operation of its own: new ids,
/// and a construction the Boolean indexes by profile element.
fn rebuilt(s: &Solid, operation: OperationId) -> Result<Solid> {
    let profile = s
        .profile()
        .ok_or_else(|| boolean_unsupported("a solid other than a prism"))?;
    Ok(Solid::extrude_with(
        operation,
        profile.clone(),
        s.frame(),
        s.start_offset(),
        s.end_offset(),
    )?
    .0)
}

/// `object op tool`, under `operations` (the Boolean's, then those of
/// rebuilt arguments). DRAW's shapes are distinct whatever made them; the
/// kernel's ids are not (every `box` is a cuboid of the unspecified
/// operation, a `copy` keeps its ids), and a Boolean's history needs its
/// inputs' ids apart: a tool sharing ids with the object is built again.
/// So is a Boolean result that is a prism, which the kernel's Boolean
/// indexes by its profile (its entities descend from the inputs); a stack
/// or a polyhedron is taken as it is (S9b.2, on its stored geometry), the
/// other argument built again when they share ids (two general bodies, or
/// a cone, a sphere or a torus, sharing ids are unsupported).
fn boolean(
    object: &Shape,
    tool: &Shape,
    op: BooleanOp,
    operations: [OperationId; 3],
) -> Result<Shape> {
    let (a, ua) = boolean_argument(object)?;
    let (b, ub) = boolean_argument(tool)?;
    let mut a = match object {
        Shape::Boolean { .. } if !general(a) => rebuilt(a, operations[1])?,
        _ => a.clone(),
    };
    let mut b = match tool {
        Shape::Boolean { .. } if !general(b) => rebuilt(b, operations[2])?,
        _ => b.clone(),
    };
    let ids: BTreeSet<_> = a.topology().ids().map(|(id, _)| id).collect();
    if b.topology().ids().any(|(id, _)| ids.contains(&id)) {
        if !general(&b) {
            b = rebuilt(&b, operations[2])?;
        } else if !general(&a) {
            a = rebuilt(&a, operations[1])?;
        } else if matches!(
            (object, tool),
            (Shape::Boolean { .. }, Shape::Boolean { .. })
        ) {
            return Err(boolean_unsupported("two stacks or polyhedra sharing ids"));
        } else {
            // A cone, a sphere or a torus sharing ids with the other
            // argument: nothing to extrude again.
            return Err(boolean_unsupported(
                "a solid other than a prism sharing ids",
            ));
        }
    }
    let (a, b, operation) = (&a, &b, operations[0]);
    let (solids, _) = match op {
        BooleanOp::Common => a.common(operation, b),
        BooleanOp::Fuse => a.fuse(operation, b),
        BooleanOp::Cut => a.cut(operation, b),
        BooleanOp::Tuc => b.cut(operation, a),
    }
    .map_err(boolean_failure)?;
    Ok(Shape::Boolean {
        solids,
        uncopied: ua || ub,
    })
}

/// Whether a Boolean result is what `unifysamedom` would leave: each
/// solid a prism whose boundaries have no two consecutive collinear lines
/// or arcs of one circle (the walls and edges OCCT's unifier merges; the
/// caps are one face each already), or a stack (S9a.2) or a polyhedron
/// (S9b.1), which the kernel builds unified: a stack's walls on one line
/// or circle facing one way joined across slab heights and piece ends, a
/// polyhedron's coplanar fragments facing one way joined into maximal
/// faces, and edges without a vertex where they run straight on between
/// the same two faces.
fn unified(solids: &[Solid]) -> bool {
    let collinear = |a: Point2, b: Point2, c: Point2| {
        let (u, v) = ((b.x - a.x, b.y - a.y), (c.x - b.x, c.y - b.y));
        let cross = u.0 * v.1 - u.1 * v.0;
        let scale = (u.0.hypot(u.1) * v.0.hypot(v.1)).max(f64::MIN_POSITIVE);
        cross.abs() <= 1e-12 * scale && u.0 * v.0 + u.1 * v.1 > 0.0
    };
    solids.iter().all(|s| {
        let Some(profile) = s.profile() else {
            return general(s);
        };
        std::iter::once(profile.outer())
            .chain(profile.holes())
            .all(|b| {
                if let Some(p) = b.polygon_vertices() {
                    let n = p.len();
                    return (0..n).all(|k| !collinear(p[k], p[(k + 1) % n], p[(k + 2) % n]));
                }
                let Some((points, segments)) = b.path_geometry() else {
                    return b.circle_geometry().is_some();
                };
                let n = points.len();
                (0..n).all(|k| match (&segments[k], &segments[(k + 1) % n]) {
                    (Segment::Line, Segment::Line) => {
                        !collinear(points[k], points[(k + 1) % n], points[(k + 2) % n])
                    }
                    (
                        Segment::Arc {
                            center: c0,
                            radius: r0,
                            ..
                        },
                        Segment::Arc {
                            center: c1,
                            radius: r1,
                            ..
                        },
                    ) => c0.x != c1.x || c0.y != c1.y || r0 != r1,
                    _ => true,
                })
            })
    })
}

/// Whether DRAW's counts of a shape differ from the kernel's structure for
/// a reason no mapping covers (an uncopied prism's shared shapes).
fn uncounted(shape: &Shape) -> bool {
    match shape {
        Shape::Uncopied(_) => true,
        Shape::Boolean { uncopied, .. } => *uncopied,
        Shape::Compound(items) | Shape::Split { pieces: items, .. } => items.iter().any(uncounted),
        _ => false,
    }
}

/// A `trotate` that is a whole number of quarter turns about a coordinate
/// axis: its centre and the signed permutation of coordinates it is. The
/// kernel's `RigidTransform::rotation` takes the sine and cosine of the
/// angle in radians, as OCCT's `gp_Trsf` does (the cosine of a quarter
/// turn rounds to 6.1e-17): OCCT's tolerances absorb that, but the
/// kernel's Booleans decide exactly and would find a wall turned onto
/// another's plane tilted off it (S9b.1).
#[derive(Clone, Copy)]
struct QuarterTurn {
    centre: Vec3,
    axis: usize,
    turns: u8,
}

impl QuarterTurn {
    /// The turn of a vector: negations and a permutation, exact.
    fn apply(self, mut v: Vec3) -> Vec3 {
        for _ in 0..self.turns {
            v = match self.axis {
                0 => Vec3::new(v.x, -v.z, v.y),
                1 => Vec3::new(v.z, v.y, -v.x),
                _ => Vec3::new(-v.y, v.x, v.z),
            };
        }
        v
    }
}

/// `trotate`'s numbers (centre, axis, angle in degrees) as quarter turns.
fn quarter_turn(n: &[f64]) -> Option<QuarterTurn> {
    let [cx, cy, cz, ax, ay, az, angle] = *n else {
        return None;
    };
    let axes = [ax, ay, az];
    if axes.iter().filter(|a| **a != 0.0).count() != 1 || angle % 90.0 != 0.0 || angle.abs() > 1e15
    {
        return None;
    }
    let axis = (0..3).find(|&k| axes[k] != 0.0)?;
    // A whole multiple of 90 divides exactly.
    let quarters = (angle / 90.0) as i64 * if axes[axis] > 0.0 { 1 } else { -1 };
    Some(QuarterTurn {
        centre: Vec3::new(cx, cy, cz),
        axis,
        turns: quarters.rem_euclid(4) as u8,
    })
}

/// A solid moved by `ttranslate` or `trotate`: the kernel's rigid motion
/// of it, keeping every id; a prism turned by quarter turns about a
/// coordinate axis is the same prism in its frame turned exactly, its
/// origin by DRAW's location arithmetic (the centre less its image, then
/// the origin's image plus that: one rounding each), with the same ids.
fn moved(s: &Solid, transform: RigidTransform, quarter: Option<QuarterTurn>) -> Result<Solid> {
    if let (Some(q), Some(profile)) = (quarter, s.profile()) {
        let f = s.frame();
        let shift = q.centre - q.apply(q.centre);
        let origin = Point3::ORIGIN + (q.apply(f.origin() - Point3::ORIGIN) + shift);
        let (normal, x, y) = (q.apply(f.normal()), q.apply(f.x()), q.apply(f.y()));
        let frame = Frame3::new(origin, normal, x, s.resolution())?;
        // The frame keeps its turned axes bit for bit, and the prism its ids.
        if frame.origin() == origin && frame.normal() == normal && frame.x() == x && frame.y() == y
        {
            let (turned, _) = Solid::extrude_with(
                s.operation(),
                profile.clone(),
                frame,
                s.start_offset(),
                s.end_offset(),
            )?;
            if turned.topology().body_id() == s.topology().body_id()
                && turned.topology().ids().eq(s.topology().ids())
            {
                return Ok(turned);
            }
        }
    }
    // Other solids (cones, spheres, tori) by the quarter turn's exact
    // matrix: their frames keep exact axes, as the prism's (S9d.3b.1's
    // survey: a cone's rounded axis crossed a cylinder's cap within
    // rounding of parallel).
    if let Some(q) = quarter {
        let origin = Point3::ORIGIN + q.centre;
        let exact = RigidTransform::quarter_turn(origin, q.axis, u32::from(q.turns))?;
        return Ok(s.transform_with(OperationId::UNSPECIFIED, exact)?.0);
    }
    Ok(s.transform_with(OperationId::UNSPECIFIED, transform)?.0)
}

/// A DRAW number (`Draw::Atof`): a literal, or an expression of `dset`
/// variables; one that cannot be evaluated is a capability gap.
fn draw_number(text: &str, numbers: &BTreeMap<String, f64>) -> Result<f64> {
    if let Ok(v) = text.trim_end_matches(',').parse::<f64>() {
        return Ok(v);
    }
    draw_geometry::evaluate(text, &|n| numbers.get(n).copied())
        .map_err(|e| Failure::Unsupported(format!("DRAW expression {text}: {e}")))
}

/// Constructors whose numbers may be DRAW expressions: the arguments from
/// the first index on (up to the second, when given).
fn numeric_arguments(command: &str) -> Option<(usize, usize)> {
    Some(match command {
        "box" | "pcylinder" | "pcone" | "psphere" | "ptorus" | "ttranslate" | "trotate"
        | "polyline" => (2, usize::MAX),
        "prism" => (3, 6),
        _ => return None,
    })
}

/// Whether a `psphere` is placed on a DRAW plane: its third word names
/// one (`psphere name plane R [angle1 angle2]`).
fn on_plane(session: &Session, args: &[String]) -> bool {
    args.first().is_some_and(|c| c == "psphere") && plane_named(&session.shapes, args)
}

/// Whether a command's third word names a DRAW plane.
fn plane_named(shapes: &BTreeMap<String, Shape>, args: &[String]) -> bool {
    args.len() > 2 && matches!(shapes.get(&args[2]), Some(Shape::Geometry(Geom::Plane(_))))
}

fn dispatch(session: &mut Session, args: &[String]) -> Result<String> {
    // DRAW evaluates `dset` variables in numeric arguments (S9a).
    let evaluated: Vec<String>;
    let args = match numeric_arguments(args.first().map(String::as_str).unwrap_or("")) {
        Some((first, last)) => {
            // `psphere name plane R ...`: the plane is a name, not a number.
            let first = first + usize::from(on_plane(session, args));
            evaluated = args
                .iter()
                .enumerate()
                .map(|(k, a)| {
                    if (first..last).contains(&k) && a.trim_end_matches(',').parse::<f64>().is_err()
                    {
                        draw_number(a, &session.numbers).map(|v| v.to_string())
                    } else {
                        Ok(a.clone())
                    }
                })
                .collect::<Result<_>>()?;
            &evaluated[..]
        }
        None => args,
    };
    let t = Tolerance::default();
    let command = args.first().map(String::as_str).unwrap_or("");
    if command == "explode" {
        session.explodes += 1;
    }
    let shapes = &mut session.shapes;
    match command {
        // S6: DRAW's geometric objects and the faces and edges made on them.
        "plane" if matches!(args.len(), 2 | 5 | 8 | 11) => {
            let n = numbers(&args[2..])?;
            let frame = if n.is_empty() {
                Frame3::xy()
            } else {
                placement(&n, t)?
            };
            shapes.insert(args[1].clone(), Shape::Geometry(Geom::Plane(frame)));
            Ok(String::new())
        }
        "cylinder" | "circle"
            if matches!(args.len(), 3 | 6 | 9 | 12)
                && (command == "cylinder" || args.len() > 3) =>
        {
            let n = numbers(&args[2..])?;
            let (radius, axes) = n.split_last().expect("a radius");
            let frame = if axes.is_empty() {
                Frame3::xy()
            } else {
                placement(axes, t)?
            };
            if *radius <= t.linear() {
                return Err(error("radius below the resolution"));
            }
            let g = if command == "cylinder" {
                Geom::Cylinder(frame, *radius)
            } else {
                Geom::Circle(frame, *radius)
            };
            shapes.insert(args[1].clone(), Shape::Geometry(g));
            Ok(String::new())
        }
        "line" if args.len() == 8 => {
            let n = numbers(&args[2..])?;
            let d = Vec3::new(n[3], n[4], n[5]);
            if d.length() <= 0.0 {
                return Err(error("line: null direction"));
            }
            let g = Geom::Line(Point3::new(n[0], n[1], n[2]), d * (1.0 / d.length()));
            shapes.insert(args[1].clone(), Shape::Geometry(g));
            Ok(String::new())
        }
        // S7 (lowalgos): DRAW's analytic surfaces (`anasurface`: a placement,
        // then the radius; a cone's half-angle in degrees, then its radius; a
        // torus's major and minor radii).
        "sphere" | "cone" | "torus"
            if matches!(
                (command, args.len()),
                ("sphere", 3 | 6 | 9 | 12) | ("cone" | "torus", 4 | 7 | 10 | 13)
            ) =>
        {
            let n = numbers(&args[2..])?;
            let params = if command == "sphere" { 1 } else { 2 };
            let (axes, p) = n.split_at(n.len() - params);
            let frame = if axes.is_empty() {
                Frame3::xy()
            } else {
                placement(axes, t)?
            };
            let g = match command {
                "sphere" => Geom::Sphere(frame, p[0]),
                "cone" => Geom::Cone(frame, p[0].to_radians(), p[1]),
                _ => Geom::Torus(frame, p[0], p[1]),
            };
            shapes.insert(args[1].clone(), Shape::Geometry(g));
            Ok(String::new())
        }
        "intersect" if args.len() == 4 || args.len() == 5 => {
            let surface = |name: &str| -> Result<Surface> {
                geometry(shapes, name)?
                    .surface()
                    .ok_or_else(|| Failure::Unsupported(format!("{name}: not an analytic surface")))
            };
            let (a, b) = (surface(&args[2])?, surface(&args[3])?);
            let (curves, points) = match draw_geometry::intersect(&a, &b) {
                Ok(found) => found,
                Err(why) if why.contains("computation budget") || why.contains("not yet") => {
                    return Err(Failure::Unsupported(format!("intersect: {why}")))
                }
                Err(why) => return Err(error(&why)),
            };
            if curves.is_empty() && points.is_empty() {
                return Err(error("No intersections found!"));
            }
            let name = &args[1];
            let mut names = Vec::new();
            for (k, c) in curves.into_iter().enumerate() {
                names.push(if k == 0 && points.is_empty() {
                    name.clone()
                } else {
                    format!("{name}_{}", k + 1)
                });
                shapes.insert(names[k].clone(), Shape::Geometry(Geom::Curve(c)));
            }
            // DRAW names a single curve `res`, several `res_1 .. res_n`.
            if names.len() > 1 && names[0] == *name {
                let c = shapes.remove(name).expect("the first curve");
                names[0] = format!("{name}_1");
                shapes.insert(names[0].clone(), c);
            }
            for (k, p) in points.into_iter().enumerate() {
                let n = format!("{name}_p_{}", k + 1);
                shapes.insert(
                    n.clone(),
                    Shape::Geometry(Geom::Curve(draw_geometry::Curve::Line(
                        p,
                        Vec3::new(1.0, 0.0, 0.0),
                    ))),
                );
                names.push(n);
            }
            Ok(names.iter().map(|n| format!("{n} ")).collect())
        }
        "bounds" if args.len() == 4 => {
            let c = geometry(shapes, &args[1])?
                .curve()
                .ok_or_else(|| Failure::Unsupported(format!("{}: not a curve", args[1])))?;
            let [lo, hi] = c.bounds();
            session.numbers.insert(args[2].clone(), lo);
            session.numbers.insert(args[3].clone(), hi);
            Ok(String::new())
        }
        "dval" if args.len() == 2 => {
            let numbers = &session.numbers;
            let v = draw_geometry::evaluate(&args[1], &|n| numbers.get(n).copied())
                .map_err(|e| error(&e))?;
            Ok(draw_geometry::format_g(v, 17))
        }
        "renamevar" if args.len() >= 3 && args.len() % 2 == 1 => {
            for pair in args[1..].chunks(2) {
                if let Some(shape) = shapes.remove(&pair[0]) {
                    shapes.insert(pair[1].clone(), shape);
                } else if let Some(v) = session.numbers.remove(&pair[0]) {
                    session.numbers.insert(pair[1].clone(), v);
                } else {
                    return Err(error(&format!("{} does not exist", pair[0])));
                }
            }
            Ok(String::new())
        }
        "dump" if args.len() == 2 => {
            let g = geometry(shapes, &args[1])?;
            let body = match (g.curve(), g) {
                (Some(c), _) => c.kind().to_string(),
                (None, Geom::Plane(_)) => "Plane".into(),
                (None, Geom::Cylinder(..)) => "CylindricalSurface".into(),
                (None, Geom::Sphere(..)) => "SphericalSurface".into(),
                (None, Geom::Cone(..)) => "ConicalSurface".into(),
                (None, _) => "ToroidalSurface".into(),
            };
            Ok(format!(
                "\n\n*********** Dump of {} *************\n{body}\n",
                args[1]
            ))
        }
        "xdistcs" if (6..=8).contains(&args.len()) => {
            let numbers = &session.numbers;
            let value = |s: &str| {
                draw_geometry::evaluate(s, &|n| numbers.get(n).copied()).map_err(|e| error(&e))
            };
            let c = match geometry(shapes, &args[1]).ok().and_then(Geom::curve) {
                Some(c) => c,
                None => return Ok(format!("Error: {} is not a curve!\n", args[1])),
            };
            let s = match geometry(shapes, &args[2]).ok().and_then(Geom::surface) {
                Some(s) => s,
                None => return Ok(format!("Error: {} is not a surface!\n", args[2])),
            };
            let (t1, t2) = (value(&args[3])?, value(&args[4])?);
            let n = (value(&args[5])? as usize).max(2);
            let tol = args
                .get(6)
                .map(|a| value(a))
                .transpose()?
                .unwrap_or(f64::MAX);
            let warn = args
                .get(7)
                .map(|a| value(a))
                .transpose()?
                .unwrap_or(f64::MAX);
            let mut out = String::new();
            let (mut max, mut at) = (0.0f64, t1);
            for i in 0..n {
                let t = if i + 1 == n {
                    t2
                } else {
                    t1 + (t2 - t1) / (n - 1) as f64 * i as f64
                };
                let p = c.point(t).map_err(|e| error(&e))?;
                let d = draw_geometry::distance(&s, p);
                if d > tol {
                    out += &format!("Error in {}:", args[1]);
                } else if d > warn {
                    out += "Attention (critical value of tolerance) :";
                }
                out += &format!(
                    " T={}\tD={}\n",
                    draw_geometry::format_g(t, 6),
                    draw_geometry::format_g(d, 6)
                );
                if d > max {
                    (max, at) = (d, t);
                }
            }
            out += &format!(
                "Max distance = {}\nParam = {}\n",
                draw_geometry::format_g(max, 17),
                draw_geometry::format_g(at, 17)
            );
            Ok(out)
        }
        "directory" if args.len() <= 2 => {
            let pattern = args.get(1).map(String::as_str).unwrap_or("*");
            let names: Vec<&String> = shapes
                .keys()
                .chain(session.numbers.keys())
                .filter(|n| glob(pattern, n))
                .collect();
            Ok(names.iter().map(|n| format!("{n} ")).collect())
        }
        "dsetsignal" => Ok(String::new()),
        // DRAW's numeric variables set by `dset` (the boolean group's begin
        // sets SCALE); `protect` only guards a variable against deletion.
        "dset" if args.len() >= 3 && args.len() % 2 == 1 => {
            // An expression the evaluator cannot read (a function other
            // than sqrt) is a capability gap.
            for pair in args[1..].chunks(2) {
                let v = draw_number(&pair[1], &session.numbers)?;
                session.numbers.insert(pair[0].clone(), v);
            }
            Ok(String::new())
        }
        "protect" if args.len() >= 2 => Ok(String::new()),
        // S8e: OCCT's General Fuse splitter with one plane face as the tool.
        "bclearobjects" | "bcleartools" if args.len() == 1 => {
            if command == "bclearobjects" {
                session.objects.clear();
            } else {
                session.tools.clear();
            }
            session.arguments += 1;
            Ok(String::new())
        }
        "baddobjects" | "baddtools" if args.len() >= 2 => {
            let mut added = Vec::new();
            for name in &args[1..] {
                added.push(get(shapes, name)?.clone());
            }
            if command == "baddobjects" {
                session.objects.extend(added);
            } else {
                session.tools.extend(added);
            }
            session.arguments += 1;
            Ok(String::new())
        }
        "bfillds" if args.len() == 1 => {
            if session.objects.is_empty() {
                return Ok("No objects to process\n".into());
            }
            session.filled = Some(session.arguments);
            Ok(String::new())
        }
        "bsplit" | "bapisplit" if args.len() == 2 => {
            if command == "bsplit" {
                match session.filled {
                    None => return Ok("Prepare PaveFiller first\n".into()),
                    // The arguments changed since the filler was prepared.
                    Some(g) if g != session.arguments => return Err(unsupported(args)),
                    Some(_) => {}
                }
            }
            let result = split(session, t)?;
            session.shapes.insert(args[1].clone(), result);
            Ok(String::new())
        }
        // S9a: Booleans of one object and one tool.
        "bop" if args.len() == 3 => {
            let object = get(shapes, &args[1])?.clone();
            let tool = get(shapes, &args[2])?.clone();
            session.bop = Some((object, tool));
            Ok(String::new())
        }
        "bfuse" | "bcut" | "bcommon" | "btuc" | "bopfuse" | "bopcut" | "boptuc" | "bopcommon"
        | "bbop" | "bapibop"
            if args.len()
                == match command {
                    "bbop" | "bapibop" => 3,
                    "bfuse" | "bcut" | "bcommon" | "btuc" => 4,
                    _ => 2,
                } =>
        {
            let op = match command {
                "bfuse" | "bopfuse" => BooleanOp::Fuse,
                "bcut" | "bopcut" => BooleanOp::Cut,
                "bcommon" | "bopcommon" => BooleanOp::Common,
                "btuc" | "boptuc" => BooleanOp::Tuc,
                // BOPAlgo_Operation: COMMON, FUSE, CUT, CUT21 (SECTION is 4).
                _ => match args[2].as_str() {
                    "0" => BooleanOp::Common,
                    "1" => BooleanOp::Fuse,
                    "2" => BooleanOp::Cut,
                    "3" => BooleanOp::Tuc,
                    _ => return Err(unsupported(args)),
                },
            };
            let (object, tool) = match command {
                "bfuse" | "bcut" | "bcommon" | "btuc" => (
                    get(shapes, &args[2])?.clone(),
                    get(shapes, &args[3])?.clone(),
                ),
                "bbop" | "bapibop" => {
                    if command == "bbop" {
                        match session.filled {
                            None => return Ok("Prepare PaveFiller first\n".into()),
                            // The arguments changed since the filler was prepared.
                            Some(g) if g != session.arguments => return Err(unsupported(args)),
                            Some(_) => {}
                        }
                    }
                    let ([object], [tool]) = (session.objects.as_slice(), session.tools.as_slice())
                    else {
                        return Err(boolean_unsupported("one object and one tool are supported"));
                    };
                    (object.clone(), tool.clone())
                }
                _ => session
                    .bop
                    .clone()
                    .ok_or_else(|| error("Prepare BOPAlgo_PaveFiller first"))?,
            };
            let n = session.next_operation;
            session.next_operation += 3;
            let operations = [1, 2, 3].map(|k| OperationId(n + k));
            let result = boolean(&object, &tool, op, operations)?;
            session.last = None;
            session.unkept_history = true;
            session.shapes.insert(args[1].clone(), result);
            Ok(String::new())
        }
        // OCCT's unifier leaves a Boolean result of the kernel as it is,
        // when its walls and edges are merged already (S9a).
        "unifysamedom" if args.len() == 3 => match get(shapes, &args[2])? {
            shape @ Shape::Boolean { solids, .. } if unified(solids) => {
                let shape = shape.clone();
                shapes.insert(args[1].clone(), shape);
                Ok(String::new())
            }
            _ => Err(unsupported(args)),
        },
        "mkface" if args.len() == 7 => {
            let b = numbers(&args[3..])?;
            let shape = match geometry(shapes, &args[2])?.clone() {
                Geom::Plane(frame) => {
                    if !(b[0] < b[1] && b[2] < b[3]) {
                        return Err(error("mkface: empty bounds"));
                    }
                    let rectangle = Boundary::polygon(
                        vec![
                            Point2::new(b[0], b[2]),
                            Point2::new(b[1], b[2]),
                            Point2::new(b[1], b[3]),
                            Point2::new(b[0], b[3]),
                        ],
                        t,
                    )?;
                    session.next_operation += 1;
                    let (body, _) = rusty_occt::Body::face_from_profile_with(
                        OperationId(session.next_operation),
                        Profile::new(rectangle, vec![], t)?,
                        frame,
                    )?;
                    Shape::Body {
                        body: Box::new(BodyRef {
                            tag: body.topology().body_id().to_string(),
                            topology: body.topology().clone(),
                        }),
                        resolution: t,
                        kind: "FACE",
                        kernel: Some(Box::new(body)),
                    }
                }
                Geom::Cylinder(frame, radius) => body_of(
                    cylinder_patch(frame, radius, [b[0], b[1], b[2], b[3]], t)?,
                    t,
                    &mut session.restored,
                    "FACE",
                )?,
                _ => return Err(unsupported(args)),
            };
            session.shapes.insert(args[1].clone(), shape);
            Ok(String::new())
        }
        "mkedge" if args.len() == 3 || args.len() == 5 => {
            let bounds = numbers(&args[3..])?;
            let (curve, ends) = match (geometry(shapes, &args[2])?.clone(), bounds.as_slice()) {
                (Geom::Line(p, d), [f, l]) if f < l => {
                    let (a, b) = (p + d * *f, p + d * *l);
                    (Curve3::LineSegment { start: a, end: b }, Some([a, b]))
                }
                (Geom::Circle(frame, radius), []) => (Curve3::Circle { frame, radius }, None),
                (Geom::Circle(frame, radius), [f, l]) if f < l && l - f < TAU => {
                    let at =
                        |a: f64| frame.point(Point2::new(radius * a.cos(), radius * a.sin()), 0.0);
                    (
                        Curve3::CircularArc {
                            frame,
                            radius,
                            start_angle: *f,
                            sweep_angle: l - f,
                        },
                        Some([at(*f), at(*l)]),
                    )
                }
                _ => return Err(unsupported(args)),
            };
            let shape = body_of(edge_parts(curve, ends), t, &mut session.restored, "EDGE")?;
            session.shapes.insert(args[1].clone(), shape);
            Ok(String::new())
        }
        "box" if args.len() == 5 || args.len() == 8 => {
            let n = numbers(&args[2..])?;
            let (origin, size) = if n.len() == 3 {
                (Point3::ORIGIN, Vec3::new(n[0], n[1], n[2]))
            } else {
                (Point3::new(n[0], n[1], n[2]), Vec3::new(n[3], n[4], n[5]))
            };
            let solid = Solid::box_at(origin, size, t)?;
            shapes.insert(args[1].clone(), Shape::Solid(Box::new(solid)));
            Ok(String::new())
        }
        "pcylinder" if args.len() == 4 => {
            let n = numbers(&args[2..])?;
            let frame = Frame3::new(
                Point3::ORIGIN,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
                t,
            )?;
            let solid = Solid::cylinder(frame, n[0], 0.0, n[1], t)?;
            shapes.insert(args[1].clone(), Shape::Solid(Box::new(solid)));
            Ok(String::new())
        }
        // `pcone name R1 R2 H`: a cone or frustum on the z axis, the seam
        // along x (a partial angle is not supported).
        "pcone" if args.len() == 5 => {
            let n = numbers(&args[2..])?;
            let frame = Frame3::new(
                Point3::ORIGIN,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
                t,
            )?;
            let (solid, _) =
                Solid::cone_with(OperationId::UNSPECIFIED, frame, n[0], n[1], n[2], t)?;
            shapes.insert(args[1].clone(), Shape::Solid(Box::new(solid)));
            Ok(String::new())
        }
        // `psphere name [plane] R [angle1 angle2]`: a sphere or zone about
        // the z axis or the plane's normal, centred at the origin or the
        // plane's origin, latitudes in degrees (a partial longitude is not
        // supported).
        "psphere" if matches!(args.len() - usize::from(plane_named(shapes, args)), 3 | 5) => {
            let (frame, n) = match shapes.get(&args[2]) {
                Some(Shape::Geometry(Geom::Plane(frame))) => (*frame, numbers(&args[3..])?),
                _ => (
                    Frame3::new(
                        Point3::ORIGIN,
                        Vec3::new(0.0, 0.0, 1.0),
                        Vec3::new(1.0, 0.0, 0.0),
                        t,
                    )?,
                    numbers(&args[2..])?,
                ),
            };
            let (low, high) = if n.len() == 3 {
                (n[1].to_radians(), n[2].to_radians())
            } else {
                (-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2)
            };
            let (solid, _) =
                Solid::sphere_with(OperationId::UNSPECIFIED, frame, n[0], low, high, t)?;
            shapes.insert(args[1].clone(), Shape::Solid(Box::new(solid)));
            Ok(String::new())
        }
        // `ptorus name R1 R2 [angle1 angle2] [angle]`: a torus about the z
        // axis, angles in degrees (a segment and a wedge at once are out of
        // the kernel's domain).
        "ptorus" if (4..=7).contains(&args.len()) => {
            let n = numbers(&args[2..])?;
            let frame = Frame3::new(
                Point3::ORIGIN,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
                t,
            )?;
            let (low, high, angle) = match n.len() {
                2 => (0.0, TAU, TAU),
                3 => (0.0, TAU, n[2].to_radians()),
                4 => (n[2].to_radians(), n[3].to_radians(), TAU),
                _ => (n[2].to_radians(), n[3].to_radians(), n[4].to_radians()),
            };
            let (solid, _) = Solid::torus_with(
                OperationId::UNSPECIFIED,
                frame,
                n[0],
                n[1],
                low,
                high,
                angle,
                t,
            )?;
            shapes.insert(args[1].clone(), Shape::Solid(Box::new(solid)));
            Ok(String::new())
        }
        // `tcopy` copies the geometry too, keeping what the shape shares
        // (an uncopied prism stays one): the same shape to the kernel.
        "copy" | "tcopy" if args.len() == 3 => {
            let shape = get(shapes, &args[1])?.clone();
            shapes.insert(args[2].clone(), shape);
            Ok(String::new())
        }
        "ttranslate" | "trotate" if args.len() == (if command == "ttranslate" { 5 } else { 9 }) => {
            let n = numbers(&args[2..])?;
            let transform = if command == "ttranslate" {
                RigidTransform::translation(Vec3::new(n[0], n[1], n[2]))?
            } else {
                RigidTransform::rotation(
                    Point3::new(n[0], n[1], n[2]),
                    Vec3::new(n[3], n[4], n[5]),
                    n[6].to_radians(),
                )?
            };
            let quarter = if command == "trotate" {
                quarter_turn(&n)
            } else {
                None
            };
            // DRAW moves the location and records no history.
            let moved = match get(shapes, &args[1])? {
                Shape::Solid(s) => Shape::Solid(Box::new(moved(s, transform, quarter)?)),
                Shape::Uncopied(s) => Shape::Uncopied(Box::new(moved(s, transform, quarter)?)),
                // A `profile` sketch moved before its prism (S9a): its
                // points, arc centres and plane.
                Shape::Face(p) | Shape::Wire(p) if command == "ttranslate" && p.plane.is_some() => {
                    let v = Vec3::new(n[0], n[1], n[2]);
                    let mut p = p.clone();
                    for q in &mut p.points {
                        *q = *q + v;
                    }
                    for (center, _, _) in p.arcs.iter_mut().flatten() {
                        *center = *center + v;
                    }
                    p.plane = Some(p.plane.expect("checked").transformed(transform, t)?);
                    boundary_of(&p, t)?;
                    match get(shapes, &args[1])? {
                        Shape::Face(_) => Shape::Face(p),
                        _ => Shape::Wire(p),
                    }
                }
                _ => return Err(unsupported(args)),
            };
            shapes.insert(args[1].clone(), moved);
            Ok(String::new())
        }
        "polyline" if args.len() >= 11 && (args.len() - 2) % 3 == 0 => {
            let n = numbers(&args[2..])?;
            let mut points: Vec<Point3> =
                n.chunks(3).map(|c| Point3::new(c[0], c[1], c[2])).collect();
            // Only closed polylines (last point repeats the first) make faces.
            if points.first().unwrap().distance(*points.last().unwrap()) > t.linear() {
                return Err(unsupported(args));
            }
            points.pop();
            let count = points.len() as u64;
            let base = session.next_label;
            session.next_label += 2 * count + 1;
            let labels = BoundaryLabels {
                boundary: InputLabel(base),
                segments: (0..count).map(|k| InputLabel(base + 1 + k)).collect(),
                vertices: (0..count)
                    .map(|k| InputLabel(base + 1 + count + k))
                    .collect(),
            };
            let arcs = vec![None; points.len()];
            shapes.insert(
                args[1].clone(),
                Shape::Wire(Polyline {
                    points,
                    labels,
                    arcs,
                    plane: None,
                }),
            );
            Ok(String::new())
        }
        // Upstream's sketch command (BRepTest_CurveCommands.cxx, `profile`):
        // a moving point and direction in a plane, closed by a line to the
        // first point, a face unless `W`. Lines and arcs (S5); `S` (a face's
        // surface) and open wires (`WW`) are not supported.
        "profile" if args.len() >= 3 => {
            let words = &args[2..];
            let (mut x0, mut y0, mut x, mut y, mut dx, mut dy): (f64, f64, f64, f64, f64, f64) =
                (0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
            let (mut origin, mut normal, mut xdir) = (
                Point3::ORIGIN,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
            );
            let (mut face, mut first, mut stay_first) = (true, true, false);
            let mut points: Vec<Point2> = Vec::new();
            let mut arcs: Vec<Option<(Point2, f64, bool)>> = Vec::new();
            let confusion = 1e-7;
            let mut i = 0;
            // Draw::Atof: a literal or an expression of `dset` variables.
            let variables = &session.numbers;
            let value = |k: usize| -> Result<f64> {
                draw_number(
                    words
                        .get(k)
                        .ok_or_else(|| error("profile : bad number of arguments"))?,
                    variables,
                )
            };
            enum Move {
                Line(f64),
                Circle(f64, f64),
                None,
            }
            loop {
                let closing = i >= words.len();
                let mut step = Move::None;
                if closing {
                    let (ex, ey) = (x0 - x, y0 - y);
                    let length = (ex * ex + ey * ey).sqrt();
                    if length <= confusion {
                        break;
                    }
                    dx = ex / length;
                    dy = ey / length;
                    step = Move::Line(length);
                } else {
                    let code = words[i].to_ascii_uppercase();
                    let second = code.chars().nth(1);
                    match code.chars().next().unwrap_or(' ') {
                        'F' => {
                            if !first {
                                return Err(error(
                                    "profile: The F instruction must precede all moves",
                                ));
                            }
                            (x0, y0) = (value(i + 1)?, value(i + 2)?);
                            (x, y) = (x0, y0);
                            stay_first = true;
                            i += 2;
                        }
                        'O' => {
                            origin = Point3::new(value(i + 1)?, value(i + 2)?, value(i + 3)?);
                            stay_first = true;
                            i += 3;
                        }
                        'P' => {
                            normal = Vec3::new(value(i + 1)?, value(i + 2)?, value(i + 3)?);
                            xdir = Vec3::new(value(i + 4)?, value(i + 5)?, value(i + 6)?);
                            stay_first = true;
                            i += 6;
                        }
                        'X' => {
                            let mut length = value(i + 1)?;
                            if second == Some('X') {
                                length -= x;
                            }
                            (dx, dy) = (1.0, 0.0);
                            step = Move::Line(length);
                            i += 1;
                        }
                        'Y' => {
                            let mut length = value(i + 1)?;
                            if second == Some('Y') {
                                length -= y;
                            }
                            (dx, dy) = (0.0, 1.0);
                            step = Move::Line(length);
                            i += 1;
                        }
                        'L' => {
                            step = Move::Line(value(i + 1)?);
                            i += 1;
                        }
                        'T' => {
                            let (mut vx, mut vy) = (value(i + 1)?, value(i + 2)?);
                            if second == Some('T') {
                                vx -= x;
                                vy -= y;
                            }
                            let length = (vx * vx + vy * vy).sqrt();
                            if length > confusion {
                                (dx, dy) = (vx / length, vy / length);
                                step = Move::Line(length);
                            }
                            i += 2;
                        }
                        'R' => {
                            let angle = value(i + 1)?.to_radians();
                            if second == Some('R') {
                                (dx, dy) = (angle.cos(), angle.sin());
                            } else {
                                let (c, s) = (angle.cos(), angle.sin());
                                (dx, dy) = (c * dx - s * dy, s * dx + c * dy);
                            }
                            i += 1;
                        }
                        'D' => {
                            let (vx, vy) = (value(i + 1)?, value(i + 2)?);
                            let length = (vx * vx + vy * vy).sqrt();
                            if length > confusion {
                                (dx, dy) = (vx / length, vy / length);
                            }
                            i += 2;
                        }
                        'C' => {
                            let radius = value(i + 1)?;
                            if radius.abs() > confusion {
                                step = Move::Circle(radius, value(i + 2)?.to_radians());
                            }
                            i += 2;
                        }
                        'I' => {
                            let target = value(i + 1)?;
                            let (along, from) = match second {
                                Some('X') => (dx, x),
                                Some('Y') => (dy, y),
                                _ => return Err(unsupported(args)),
                            };
                            if along.abs() < confusion {
                                return Err(error("Profile : cannot intersect"));
                            }
                            step = Move::Line((target - from) / along);
                            i += 1;
                        }
                        'W' => {
                            if second == Some('W') {
                                return Err(unsupported(args));
                            }
                            face = false;
                            i = words.len() - 1;
                        }
                        _ => return Err(unsupported(args)),
                    }
                }
                match step {
                    Move::Line(mut length) => {
                        if length < 0.0 {
                            length = -length;
                            (dx, dy) = (-dx, -dy);
                        }
                        points.push(Point2::new(x, y));
                        arcs.push(None);
                        x += length * dx;
                        y += length * dy;
                    }
                    Move::Circle(mut radius, mut angle) => {
                        let mut sense = true;
                        if radius < 0.0 {
                            radius = -radius;
                            sense = !sense;
                            (dx, dy) = (-dx, -dy);
                        }
                        // gp_Ax2d: the centre left of the direction, the
                        // x axis from the centre to the point.
                        let center = Point2::new(x - radius * dy, y + radius * dx);
                        let (ux, uy) = (dy, -dx);
                        if angle < 0.0 {
                            angle = -angle;
                            sense = !sense;
                        }
                        // A direct circle's y axis turns x a quarter left.
                        let (vx, vy) = if sense { (-uy, ux) } else { (uy, -ux) };
                        points.push(Point2::new(x, y));
                        arcs.push(Some((center, radius, sense)));
                        let (c, s) = (angle.cos(), angle.sin());
                        x = center.x + radius * (c * ux + s * vx);
                        y = center.y + radius * (c * uy + s * vy);
                        (dx, dy) = (-s * ux + c * vx, -s * uy + c * vy);
                    }
                    Move::None => {}
                }
                if closing {
                    break;
                }
                first = stay_first;
                stay_first = false;
                i += 1;
            }
            let plane = Frame3::new(origin, normal, xdir, t)?;
            let at = |q: Point2| plane.point(q, 0.0);
            let count = points.len() as u64;
            let base = session.next_label;
            session.next_label += 2 * count + 1;
            let labels = BoundaryLabels {
                boundary: InputLabel(base),
                segments: (0..count).map(|k| InputLabel(base + 1 + k)).collect(),
                vertices: (0..count)
                    .map(|k| InputLabel(base + 1 + count + k))
                    .collect(),
            };
            let polyline = Polyline {
                points: points.iter().map(|q| at(*q)).collect(),
                labels,
                arcs: arcs
                    .iter()
                    .map(|a| a.map(|(c, r, ccw)| (at(c), r, ccw)))
                    .collect(),
                plane: Some(plane),
            };
            // The boundary must be valid now, as OCCT's face is.
            boundary_of(&polyline, t)?;
            let shape = if face {
                Shape::Face(polyline)
            } else {
                Shape::Wire(polyline)
            };
            shapes.insert(args[1].clone(), shape);
            Ok(String::new())
        }
        "mkplane" if args.len() == 3 => match get(shapes, &args[2])? {
            Shape::Wire(p) => {
                plane_of(p, t)?;
                let face = Shape::Face(p.clone());
                shapes.insert(args[1].clone(), face);
                Ok(String::new())
            }
            _ => Err(unsupported(args)),
        },
        // DRAW's Copy builds distinct end entities, as the kernel does;
        // without it OCCT reuses the start shapes under a moved location
        // (an uncopied prism: its counts and history are unsupported).
        "prism"
            if args.len() == 6
                || args.len() == 7 && args[6].to_ascii_lowercase().starts_with('c') =>
        {
            let copy = args.len() == 7;
            let Shape::Face(face) = get(shapes, &args[2])?.clone() else {
                return Err(unsupported(args));
            };
            let n = numbers(&args[3..6])?;
            let vector = Vec3::new(n[0], n[1], n[2]);
            let (frame, boundary) = boundary_of(&face, t)?;
            let height = vector.dot(frame.normal());
            // Only prisms normal to the profile plane are supported.
            if (vector - frame.normal() * height).length() > t.linear() {
                return Err(unsupported(args));
            }
            let profile = Profile::new(boundary, vec![], t)?;
            session.next_operation += 1;
            let operation = OperationId(session.next_operation);
            let (output, history) = Solid::extrude_with(operation, profile, frame, 0.0, height)?;
            if !copy {
                session.last = None;
                session.unkept_history = true;
                session
                    .shapes
                    .insert(args[1].clone(), Shape::Uncopied(Box::new(output)));
                return Ok(String::new());
            }
            session.last = Some(Saved {
                history,
                output: output.clone(),
                face,
            });
            session.unkept_history = false;
            session
                .shapes
                .insert(args[1].clone(), Shape::Solid(Box::new(output)));
            Ok(String::new())
        }
        "explode" if args.len() == 3 => {
            let polyline = match get(shapes, &args[1])? {
                Shape::Wire(p) | Shape::Face(p) => p.clone(),
                shape @ (Shape::Solid(_) | Shape::Body { .. }) => {
                    let s = match shape {
                        Shape::Solid(s) => BodyRef::of(s),
                        Shape::Body { body, .. } => (**body).clone(),
                        _ => unreachable!("matched above"),
                    };
                    // DBRep's explode reads the type from its first letter.
                    let kind = match args[2].bytes().next().map(|b| b.to_ascii_lowercase()) {
                        Some(b'f') => "face",
                        Some(b'e') => "edge",
                        _ => return Err(unsupported(args)),
                    };
                    if session.selector.is_none() {
                        session.selector = Some(load_selector()?);
                    }
                    let picks: Vec<&Pick> = session
                        .selector
                        .as_ref()
                        .unwrap()
                        .iter()
                        .filter(|p| p.explode == session.explodes)
                        .collect();
                    if picks.is_empty() || picks.iter().any(|p| p.kind != kind) {
                        return Err(unsupported(args));
                    }
                    let mut names = Vec::new();
                    for (k, pick) in picks.iter().enumerate() {
                        let name = format!("{}_{}", args[1], k + 1);
                        session.shapes.insert(name.clone(), select(&s, pick));
                        names.push(name);
                    }
                    return Ok(names.join(" "));
                }
                _ => return Err(unsupported(args)),
            };
            // Profile edges are chords: an arc's is not its edge.
            if polyline.has_arcs() && args[2].eq_ignore_ascii_case("e") {
                return Err(unsupported(args));
            }
            let n = polyline.points.len();
            let items: Vec<Shape> = match args[2].to_ascii_lowercase().as_str() {
                "e" => (0..n)
                    .map(|k| Shape::ProfileEdge {
                        label: polyline.labels.segments[k],
                        start: polyline.points[k],
                        end: polyline.points[(k + 1) % n],
                    })
                    .collect(),
                "v" => (0..n)
                    .map(|k| Shape::ProfileVertex {
                        label: polyline.labels.vertices[k],
                    })
                    .collect(),
                _ => return Err(unsupported(args)),
            };
            let mut names = Vec::new();
            for (k, item) in items.into_iter().enumerate() {
                let name = format!("{}_{}", args[1], k + 1);
                shapes.insert(name.clone(), item);
                names.push(name);
            }
            Ok(names.join(" "))
        }
        "savehistory" if args.len() == 2 && session.unkept_history => Err(unsupported(args)),
        "savehistory" if args.len() == 2 => {
            let saved = session
                .last
                .clone()
                .ok_or_else(|| error("No history has been prepared yet."))?;
            session.histories.insert(args[1].clone(), saved);
            Ok(String::new())
        }
        "generated" | "modified" if args.len() == 4 => {
            let saved = session
                .histories
                .get(&args[2])
                .ok_or_else(|| {
                    Failure::Error(format!("History with the name {} does not exist.", args[2]))
                })?
                .clone();
            let found = match (command, get(shapes, &args[3])?) {
                // An extrusion modifies no profile element.
                ("modified", _) => Vec::new(),
                (_, Shape::ProfileEdge { label, .. }) => {
                    history_output(&saved, *label, &[Role::Wall])
                }
                (_, Shape::ProfileVertex { label, .. }) => {
                    history_output(&saved, *label, &[Role::Vertical, Role::Seam])
                }
                (_, Shape::Face(p)) if p.labels == saved.face.labels => {
                    vec![Shape::Solid(Box::new(saved.output.clone()))]
                }
                _ => return Err(unsupported(args)),
            };
            let shape = match found.len() {
                0 if command == "generated" => {
                    return Ok("No shapes were generated from the shape.".into())
                }
                0 => return Ok("The shape has not been modified.".into()),
                1 => found.into_iter().next().unwrap(),
                _ => Shape::Compound(found),
            };
            session.shapes.insert(args[1].clone(), shape);
            Ok(String::new())
        }
        "isdeleted" if args.len() == 3 => {
            session.histories.get(&args[1]).ok_or_else(|| {
                Failure::Error(format!("History with the name {} does not exist.", args[1]))
            })?;
            get(shapes, &args[2])?;
            // Every profile element of an extrusion is generated from, not removed.
            Ok("Not deleted.".into())
        }
        "isdraw" if args.len() == 2 => Ok(if shapes.contains_key(&args[1]) {
            "1"
        } else {
            "0"
        }
        .into()),
        "whatis" if args.len() == 2 => {
            // DRAW's `dtyp`: geometry and numbers as DrawTrSurf and
            // Draw_Number describe themselves; nothing for an unknown name.
            match shapes.get(&args[1]) {
                Some(Shape::Geometry(g)) => {
                    let what = if g.curve().is_some() {
                        " a 3d curve"
                    } else {
                        "a surface"
                    };
                    return Ok(format!("{} is a {what}", args[1]));
                }
                None if session.numbers.contains_key(&args[1]) => {
                    return Ok(format!("{} is a numeric", args[1]))
                }
                None => return Ok(format!("{} is a ", args[1])),
                _ => {}
            }
            let shape = get(shapes, &args[1])?;
            Ok(format!(
                "{} is a shape {} FORWARD Free Modified",
                args[1],
                type_name(shape)
            ))
        }
        "checkshape" if args.len() == 2 => {
            // Every solid of the shape validates against its own resolution.
            fn bodies<'a>(shape: &'a Shape, out: &mut Vec<(&'a Topology, Tolerance)>) -> bool {
                match shape {
                    Shape::Solid(s) | Shape::Uncopied(s) => {
                        out.push((s.topology(), s.resolution()))
                    }
                    // An empty Boolean result is a valid empty compound.
                    Shape::Boolean { solids, .. } => {
                        out.extend(solids.iter().map(|s| (s.topology(), s.resolution())))
                    }
                    Shape::Body {
                        body, resolution, ..
                    } => out.push((&body.topology, *resolution)),
                    Shape::Compound(items) | Shape::Split { pieces: items, .. } => {
                        return items.iter().all(|i| bodies(i, out))
                    }
                    _ => return false,
                }
                true
            }
            // S6: a profile face or wire validates as the face or wire body
            // the kernel builds of it.
            if let Shape::Face(p) | Shape::Wire(p) = get(shapes, &args[1])? {
                let (frame, boundary) = boundary_of(p, t)?;
                let body = if matches!(get(shapes, &args[1])?, Shape::Face(_)) {
                    rusty_occt::Body::face_from_profile_with(
                        OperationId::UNSPECIFIED,
                        Profile::new(boundary, vec![], t)?,
                        frame,
                    )?
                } else {
                    rusty_occt::Body::wire_from_boundary_with(
                        OperationId::UNSPECIFIED,
                        boundary,
                        frame,
                        t,
                    )?
                }
                .0;
                body.topology().validate(t)?;
                return Ok("This shape seems to be valid".into());
            }
            let mut found = Vec::new();
            if !bodies(get(shapes, &args[1])?, &mut found) {
                return Err(unsupported(args));
            }
            for (topology, resolution) in found {
                topology.validate(resolution)?;
            }
            Ok("This shape seems to be valid".into())
        }
        "restore" if args.len() == 3 => {
            let shape = restore(&args[1], &mut session.restored)?;
            session.shapes.insert(args[2].clone(), shape);
            Ok(String::new())
        }
        "nbshapes" if args.len() == 2 => {
            let shape = get(shapes, &args[1])?;
            if uncounted(shape) {
                return Err(Failure::Unsupported(format!(
                    "nbshapes {}: an uncopied prism's shapes OCCT shares",
                    args[1]
                )));
            }
            let p = parts(shape);
            let mut counts = [
                ("VERTEX", p.vertices.len()),
                ("EDGE", p.edges.len()),
                ("WIRE", p.wires),
                ("FACE", p.faces.len()),
                ("SHELL", p.shells),
                ("SOLID", p.solids),
                ("COMPSOLID", 0),
                ("COMPOUND", p.compounds),
            ];
            // A whole solid reports the kernel's synthesized OCCT counts; the
            // per-shape collection must agree with them.
            let whole = match shape {
                Shape::Solid(s) => Some(s.topology()),
                Shape::Body { body, .. } => Some(&body.topology),
                _ => None,
            };
            if let Some(t) = whole {
                let c = t.occt_counts();
                let synthesized = [c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids];
                if counts[..6].iter().map(|(_, n)| *n).ne(synthesized) {
                    return Err(error("shape counts disagree with the count synthesizer"));
                }
                for (k, n) in synthesized.into_iter().enumerate() {
                    counts[k].1 = n;
                }
            }
            let mut result = format!("Number of shapes in {}\n", args[1]);
            for (kind, count) in counts {
                result.push_str(&format!(" {kind:<10}: {count}\n"));
            }
            result.push_str(&format!(
                " SHAPE     : {}\n",
                counts.iter().map(|(_, n)| n).sum::<usize>()
            ));
            Ok(result)
        }
        // A split's volume (S8e): its solids', a wire's being zero; OCCT
        // gives a face a signed volume against the origin, not supported.
        "vprops"
            if (args.len() == 2 || args.len() == 3)
                && matches!(shapes.get(&args[1]), Some(Shape::Split { .. })) =>
        {
            if args.len() == 3 && numbers(&args[2..])?[0] <= 0.0 {
                return Err(unsupported(args));
            }
            let Some(Shape::Split { pieces, .. }) = shapes.get(&args[1]) else {
                unreachable!("matched above")
            };
            let mut volume = 0.0;
            for piece in pieces {
                match piece {
                    Shape::Solid(s) => volume += s.mass_properties().volume,
                    Shape::Body { kind: "WIRE", .. } => {}
                    _ => return Err(unsupported(args)),
                }
            }
            Ok(props(volume))
        }
        // A Boolean's volume and centre (S9a): its solids', an empty
        // result's zero.
        "vprops"
            if (args.len() == 2 || args.len() == 3)
                && matches!(shapes.get(&args[1]), Some(Shape::Boolean { .. })) =>
        {
            if args.len() == 3 && numbers(&args[2..])?[0] <= 0.0 {
                return Err(unsupported(args));
            }
            let Some(Shape::Boolean { solids, .. }) = shapes.get(&args[1]) else {
                unreachable!("matched above")
            };
            let (mut volume, mut moment) = (0.0, Vec3::new(0.0, 0.0, 0.0));
            for s in solids {
                let m = s.mass_properties();
                volume += m.volume;
                moment = moment + (m.centroid - Point3::ORIGIN) * m.volume;
            }
            let c = if volume > 0.0 {
                moment * (1.0 / volume)
            } else {
                moment
            };
            Ok(format!(
                "{}\nCenter of gravity :\nX = {:.17e}\nY = {:.17e}\nZ = {:.17e}\n",
                props(volume),
                c.x,
                c.y,
                c.z
            ))
        }
        "vprops"
            if (args.len() == 2 || args.len() == 3)
                && matches!(
                    get(shapes, &args[1])?,
                    Shape::Solid(_) | Shape::Uncopied(_) | Shape::Body { kind: "SOLID", .. }
                ) =>
        {
            if args.len() == 3 && numbers(&args[2..])?[0] <= 0.0 {
                return Err(unsupported(args));
            }
            // Epsilon controls OCCT quadrature. The kernel's properties are
            // exact or certified enclosures, reported at their midpoints.
            let m = match get(shapes, &args[1])? {
                Shape::Body { body, .. } => body
                    .topology
                    .mass_enclosure()
                    .ok_or_else(|| error("mass properties not certified"))?
                    .midpoints(),
                Shape::Solid(s) | Shape::Uncopied(s) => s.mass_properties(),
                _ => unreachable!("matched above"),
            };
            Ok(format!("Mass : {:.17e}\n\nCenter of gravity :\nX = {:.17e}\nY = {:.17e}\nZ = {:.17e}\nMatrix of Inertia :\n{:.17e} {:.17e} {:.17e}\n{:.17e} {:.17e} {:.17e}\n{:.17e} {:.17e} {:.17e}\n",
                m.volume, m.centroid.x, m.centroid.y, m.centroid.z,
                m.inertia[0][0], m.inertia[0][1], m.inertia[0][2],
                m.inertia[1][0], m.inertia[1][1], m.inertia[1][2],
                m.inertia[2][0], m.inertia[2][1], m.inertia[2][2]))
        }
        "sprops" | "lprops" if args.len() == 2 || args.len() == 3 => {
            if args.len() == 3 && numbers(&args[2..])?[0] <= 0.0 {
                return Err(unsupported(args));
            }
            let shape = get(shapes, &args[1])?;
            let p = parts(shape);
            let mass = if command == "sprops" {
                p.area
            } else {
                p.length
            };
            // A selected face or edge also reports its centre of gravity.
            let centre = match (command, shape) {
                (
                    "sprops",
                    Shape::Sub {
                        body,
                        slot: Slot::Face(f),
                    },
                ) => Some(face_centre(&body.topology, f.index()).1),
                (
                    "lprops",
                    Shape::Sub {
                        body,
                        slot: Slot::Edge(e),
                    },
                ) => Some(edge_centre(&body.topology.edges()[e.index()].curve).1),
                _ => None,
            };
            Ok(match centre {
                Some(c) => format!(
                    "{}\nCenter of gravity : \nX = {:.17e}\nY = {:.17e}\nZ = {:.17e}\n",
                    props(mass),
                    c.x,
                    c.y,
                    c.z
                ),
                None => props(mass),
            })
        }
        "isbbinterf" if args.len() == 3 => {
            let a = solid(shapes, &args[1], args)?.bounds();
            let b = solid(shapes, &args[2], args)?.bounds();
            Ok(if a.intersects(b, t)? {
                "The shapes are interfered by AABB.\n"
            } else {
                "The shapes are NOT interfered by AABB.\n"
            }
            .into())
        }
        _ => Err(unsupported(args)),
    }
}

// Hex-encoded UTF-8 tokens keep Tcl names, spaces and newlines out of framing.
fn unhex(s: &str) -> Result<String> {
    let bytes = s.as_bytes();
    if bytes.len() % 2 != 0 || !bytes.is_ascii() {
        return Err(Failure::Error("malformed wire token".into()));
    }
    let decoded = (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|_| Failure::Error("malformed wire hex".into()))?;
    String::from_utf8(decoded).map_err(|_| Failure::Error("malformed wire UTF-8".into()))
}
fn hex(s: &str) -> String {
    s.as_bytes().iter().map(|b| format!("{b:02x}")).collect()
}
fn main() -> io::Result<()> {
    let mut session = Session::default();
    let mut output = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let result = line?
            .split(' ')
            .map(unhex)
            .collect::<Result<Vec<_>>>()
            .and_then(|args| dispatch(&mut session, &args));
        let (status, message) = match result {
            Ok(value) => ("OK", value),
            Err(Failure::Error(value)) => ("ERROR", value),
            Err(Failure::Unsupported(value)) => ("UNSUPPORTED", value),
        };
        writeln!(output, "{status} {}", hex(&message))?;
        output.flush()?;
    }
    Ok(())
}
