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
use rusty_occt::topology::{Curve2, Curve3, Slot, Surface, Topology};
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
        solid: Box<Solid>,
        slot: Slot,
    },
    Compound(Vec<Shape>),
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
    let face = &t.faces()[face];
    match face.surface {
        Surface::Plane(_) => {
            0.5 * face
                .loops
                .iter()
                .flatten()
                .map(|u| moments(&u.pcurve)[0])
                .sum::<f64>()
                .abs()
        }
        Surface::Cylinder { radius, .. } => {
            let ends: Vec<Point2> = face.loops[0].iter().map(|u| u.pcurve.point(0.0)).collect();
            let span = |f: fn(&Point2) -> f64| {
                let v: Vec<f64> = ends.iter().map(f).collect();
                v.iter().cloned().fold(f64::MIN, f64::max)
                    - v.iter().cloned().fold(f64::MAX, f64::min)
            };
            radius * span(|p| p.x) * span(|p| p.y)
        }
    }
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
    let key = |s: &Solid, slot: Slot| {
        format!(
            "{}:{}",
            s.topology().body_id(),
            s.topology().id_of(slot).unwrap()
        )
    };
    let add_edge = |s: &Solid, e: usize, parts: &mut Parts| {
        let t = s.topology();
        let edge = &t.edges()[e];
        if parts
            .edges
            .insert(key(s, Slot::Edge(rusty_occt::topology::EdgeId::new(e))))
        {
            parts.length += edge_length(&edge.curve);
        }
        for v in [edge.start, edge.end] {
            parts.vertices.insert(key(s, Slot::Vertex(v)));
        }
    };
    let add_face = |s: &Solid, f: usize, parts: &mut Parts| {
        let t = s.topology();
        if parts
            .faces
            .insert(key(s, Slot::Face(rusty_occt::topology::FaceId::new(f))))
        {
            parts.area += face_area(t, f);
            parts.wires += t.faces()[f].loops.len();
            for u in t.faces()[f].loops.iter().flatten() {
                add_edge(s, u.edge.index(), parts);
            }
        }
    };
    let polyline = |p: &Polyline, parts: &mut Parts, face: bool| {
        let n = p.points.len();
        for k in 0..n {
            parts
                .vertices
                .insert(format!("L{}", p.labels.vertices[k].0));
            if parts.edges.insert(format!("L{}", p.labels.segments[k].0)) {
                parts.length += p.points[k].distance(p.points[(k + 1) % n]);
            }
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
            for f in 0..s.topology().faces().len() {
                add_face(s, f, parts);
            }
        }
        Shape::Wire(p) => polyline(p, parts, false),
        Shape::Face(p) => polyline(p, parts, true),
        Shape::ProfileVertex { label, .. } => {
            parts.vertices.insert(format!("L{}", label.0));
        }
        Shape::ProfileEdge { label, start, end } => {
            if parts.edges.insert(format!("L{}", label.0)) {
                parts.length += start.distance(*end);
            }
        }
        Shape::Sub { solid, slot } => match slot {
            Slot::Vertex(v) => {
                parts.vertices.insert(key(solid, Slot::Vertex(*v)));
            }
            Slot::Edge(e) => add_edge(solid, e.index(), parts),
            Slot::Face(f) => add_face(solid, f.index(), parts),
        },
        Shape::Compound(items) => {
            parts.compounds += 1;
            for item in items {
                collect(item, parts);
            }
        }
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
        Shape::Solid(_) => "SOLID",
        Shape::Wire(_) => "WIRE",
        Shape::Face(_) => "FACE",
        Shape::ProfileVertex { .. } => "VERTEX",
        Shape::ProfileEdge { .. } => "EDGE",
        Shape::Sub { slot, .. } => match slot {
            Slot::Vertex(_) => "VERTEX",
            Slot::Edge(_) => "EDGE",
            Slot::Face(_) => "FACE",
        },
        Shape::Compound(_) => "COMPOUND",
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
                    solid: Box::new(saved.output.clone()),
                    slot: t.slot_of(*to).unwrap(),
                })
            }
            _ => None,
        })
        .collect()
}

fn dispatch(session: &mut Session, args: &[String]) -> Result<String> {
    let t = Tolerance::default();
    let command = args.first().map(String::as_str).unwrap_or("");
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
            solid(shapes, &args[1], args)?.topology().validate(t)?;
            Ok("This shape seems to be valid".into())
        }
        "nbshapes" if args.len() == 2 => {
            let shape = get(shapes, &args[1])?;
            let p = parts(shape);
            let counts = [
                ("VERTEX", p.vertices.len()),
                ("EDGE", p.edges.len()),
                ("WIRE", p.wires),
                ("FACE", p.faces.len()),
                ("SHELL", p.solids),
                ("SOLID", p.solids),
                ("COMPSOLID", 0),
                ("COMPOUND", p.compounds),
            ];
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
            Ok(props(if command == "sprops" {
                p.area
            } else {
                p.length
            }))
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
