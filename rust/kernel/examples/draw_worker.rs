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
use rusty_occt::topology::{Curve2, Curve3, EdgeId, FaceId, Loop, Slot, Surface, Topology};
use rusty_occt::{
    Boundary, BoundaryLabels, Frame3, Point2, Point3, Profile, RigidTransform, Solid, Tolerance,
    Vec3,
};
use std::collections::{BTreeMap, BTreeSet};
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
    /// resolution its OCCT tolerances give it.
    Body {
        body: Box<BodyRef>,
        resolution: Tolerance,
    },
    Compound(Vec<Shape>),
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
}

fn unsupported(args: &[String]) -> Failure {
    Failure::Unsupported(format!("unsupported DRAW signature: {}", args.join(" ")))
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
    if let Some(Shape::Lost(why)) = shapes.get(name) {
        return Err(Failure::Unsupported(format!("{name}: {why}")));
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
    }
}

fn face_area(t: &Topology, face: usize) -> f64 {
    let fins = t.face_fins(FaceId::new(face));
    let fins = fins.iter().flatten();
    match t.faces()[face].surface {
        Surface::Plane(_) => 0.5 * fins.map(|u| moments(&u.pcurve)[0]).sum::<f64>().abs(),
        // Wall pcurves are lines on the universal cover: the area is the
        // radius times the periodic area -∮ v du, wound loops included.
        Surface::Cylinder { radius, .. } => {
            let periodic = fins
                .map(|u| match u.pcurve {
                    Curve2::LineSegment { start: a, end: b } => -0.5 * (a.y + b.y) * (b.x - a.x),
                    Curve2::CircularArc { .. } => f64::NAN,
                })
                .sum::<f64>();
            radius * periodic.abs()
        }
    }
}

/// The seam OCCT's wound faces carry and the kernel does not: its length (the
/// face's extent across the winding) when the face has one.
fn seam_length(t: &Topology, face: usize) -> Option<f64> {
    let face = &t.faces()[face];
    let mut wound = false;
    let mut v = (f64::MAX, f64::MIN);
    for l in &face.loops {
        if let Loop::Edges { fins, winding } = &t.loops()[l.index()] {
            wound |= winding != &[0, 0];
            for k in fins {
                let y = t.fins()[k.index()].pcurve.point(0.0).y;
                v = (v.0.min(y), v.1.max(y));
            }
        }
    }
    wound.then_some(v.1 - v.0)
}

fn edge_length(curve: &Curve3) -> f64 {
    match curve {
        Curve3::LineSegment { start, end } => start.distance(*end),
        Curve3::Circle { radius, .. } => TAU * radius,
        Curve3::CircularArc {
            radius,
            sweep_angle,
            ..
        } => radius * sweep_angle.abs(),
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
    solids: usize,
    compounds: usize,
    length: f64,
    area: f64,
    volume: f64,
}

fn collect(shape: &Shape, parts: &mut Parts) {
    let key = |b: &BodyRef, slot: Slot| format!("{}:{}", b.tag, b.topology.id_of(slot).unwrap());
    // Lengths sum per occurrence, as OCCT's lprops explores edges; counts are
    // of distinct shapes, as nbshapes reports them.
    let add_edge = |s: &BodyRef, e: usize, parts: &mut Parts| {
        let edge = &s.topology.edges()[e];
        parts.edges.insert(key(s, Slot::Edge(EdgeId::new(e))));
        parts.length += edge_length(&edge.curve);
        for v in [edge.start, edge.end].into_iter().flatten() {
            parts.vertices.insert(key(s, Slot::Vertex(v)));
        }
        // OCCT splits a ring edge at a seam vertex.
        if edge.is_ring() {
            parts
                .vertices
                .insert(format!("{}:seam", key(s, Slot::Edge(EdgeId::new(e)))));
        }
    };
    let add_face = |s: &BodyRef, f: usize, parts: &mut Parts| {
        let t = &s.topology;
        let face = key(s, Slot::Face(FaceId::new(f)));
        if !parts.faces.insert(face.clone()) {
            return;
        }
        parts.area += face_area(t, f);
        for fins in t.face_fins(FaceId::new(f)) {
            for u in fins {
                add_edge(s, u.edge.index(), parts);
            }
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
        for k in 0..n {
            parts
                .vertices
                .insert(format!("L{}", p.labels.vertices[k].0));
            parts.edges.insert(format!("L{}", p.labels.segments[k].0));
            parts.length += p.points[k].distance(p.points[(k + 1) % n]);
        }
        parts.wires += 1;
        if face && parts.faces.insert(format!("L{}", p.labels.boundary.0)) {
            parts.area += polygon_area(&p.points);
        }
    };
    match shape {
        Shape::Solid(s) => {
            parts.solids += 1;
            parts.volume += s.mass_properties().volume;
            let b = BodyRef::of(s);
            for f in 0..b.topology.faces().len() {
                add_face(&b, f, parts);
            }
        }
        Shape::Body { body, .. } => {
            parts.solids += 1;
            // General mass properties arrive with the cone (REVIEW_NOTES S3).
            parts.volume = f64::NAN;
            for f in 0..body.topology.faces().len() {
                add_face(body, f, parts);
            }
        }
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
                parts.vertices.insert(key(body, Slot::Vertex(*v)));
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
        Shape::Lost(_) => unreachable!("lost picks are never read"),
    }
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
        Shape::Solid(_) | Shape::Body { .. } => "SOLID",
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
    }
}

/// Area and centre of gravity of a face: -∮ G du over its pcurves for
/// ∂G/∂v = 1, u, v (planes) or 1, cos u, sin u, v (cylinders). Seams are
/// vertical, so wound faces need no closing sides.
fn face_centre(t: &Topology, face: usize) -> (f64, Point3) {
    let nodes = gauss();
    let surface = &t.faces()[face].surface;
    let mut m = [0.0; 4];
    for fin in t.face_fins(FaceId::new(face)).iter().flatten() {
        for (x, w) in &nodes {
            let (q, d) = pcurve_at(&fin.pcurve, *x);
            let g = match surface {
                Surface::Plane(_) => [q.y, q.x * q.y, q.y * q.y / 2.0, 0.0],
                Surface::Cylinder { .. } => {
                    [q.y, q.y * q.x.cos(), q.y * q.x.sin(), q.y * q.y / 2.0]
                }
            };
            for k in 0..4 {
                m[k] -= w * g[k] * d.x;
            }
        }
    }
    match surface {
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
    fn build(
        doc: &rusty_occt::occt_brep::Document,
        sub: rusty_occt::occt_brep::read::Sub,
        solids: &mut std::vec::IntoIter<rusty_occt::occt_brep::ImportedSolid>,
        serial: &mut u64,
    ) -> Result<Shape> {
        let shape = &doc.shapes[sub.shape];
        match shape.kind {
            Kind::Compound => Ok(Shape::Compound(
                shape
                    .subs
                    .iter()
                    .map(|s| build(doc, *s, solids, serial))
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
            _ => Err(Failure::Unsupported("restore: a free shape".into())),
        }
    }
    build(&doc, doc.root, &mut solids, serial)
}

fn dispatch(session: &mut Session, args: &[String]) -> Result<String> {
    let t = Tolerance::default();
    let command = args.first().map(String::as_str).unwrap_or("");
    if command == "explode" {
        session.explodes += 1;
    }
    let shapes = &mut session.shapes;
    match command {
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
        "copy" if args.len() == 3 => {
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
            // DRAW moves the location and records no history.
            let moved = solid(shapes, &args[1], args)?
                .transform_with(OperationId::UNSPECIFIED, transform)?
                .0;
            shapes.insert(args[1].clone(), Shape::Solid(Box::new(moved)));
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
            shapes.insert(args[1].clone(), Shape::Wire(Polyline { points, labels }));
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
        // without it OCCT reuses the start shapes under a moved location.
        "prism" if args.len() == 7 && args[6].to_ascii_lowercase().starts_with('c') => {
            let Shape::Face(face) = get(shapes, &args[2])?.clone() else {
                return Err(unsupported(args));
            };
            let n = numbers(&args[3..6])?;
            let vector = Vec3::new(n[0], n[1], n[2]);
            let (frame, local) = plane_of(&face, t)?;
            let height = vector.dot(frame.normal());
            // Only prisms normal to the profile plane are supported.
            if (vector - frame.normal() * height).length() > t.linear() {
                return Err(unsupported(args));
            }
            let boundary = Boundary::polygon(local, t)?.with_labels(face.labels.clone())?;
            let profile = Profile::new(boundary, vec![], t)?;
            session.next_operation += 1;
            let operation = OperationId(session.next_operation);
            let (output, history) = Solid::extrude_with(operation, profile, frame, 0.0, height)?;
            session.last = Some(Saved {
                history,
                output: output.clone(),
                face,
            });
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
                    Shape::Solid(s) => out.push((s.topology(), s.profile().tolerance())),
                    Shape::Body { body, resolution } => out.push((&body.topology, *resolution)),
                    Shape::Compound(items) => return items.iter().all(|i| bodies(i, out)),
                    _ => return false,
                }
                true
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
            let p = parts(shape);
            let mut counts = [
                ("VERTEX", p.vertices.len()),
                ("EDGE", p.edges.len()),
                ("WIRE", p.wires),
                ("FACE", p.faces.len()),
                ("SHELL", p.solids),
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
        "vprops"
            if (args.len() == 2 || args.len() == 3)
                && matches!(get(shapes, &args[1])?, Shape::Solid(_)) =>
        {
            if args.len() == 3 && numbers(&args[2..])?[0] <= 0.0 {
                return Err(unsupported(args));
            }
            // Epsilon controls OCCT quadrature. These box properties are analytic.
            let m = solid(shapes, &args[1], args)?.mass_properties();
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
