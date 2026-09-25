//! Case protocol written by rust/tools/identity_reference.py::encode_case, and
//! entity rows with structural locators found from geometry alone.
use rusty_occt::history::{History, Relation};
use rusty_occt::identity::{InputLabel, OperationId, Parent, ProfileElement, Role};
use rusty_occt::topology::{Curve3, Slot, Surface};
use rusty_occt::{
    Boundary, BoundaryLabels, Frame3, Point2, Point3, Profile, RigidTransform, Solid, Tolerance,
    Vec3,
};

pub struct CaseSpec {
    pub name: String,
    pub tolerance: Tolerance,
    pub operation: OperationId,
    pub frame: [f64; 9],
    pub start: f64,
    pub end: f64,
    pub boundaries: Vec<Boundary>,
    pub box_at: Option<([f64; 3], [f64; 3])>,
    pub transforms: Vec<RigidTransform>,
}

fn labels(words: &[&str]) -> BoundaryLabels {
    let n = |w: &&str| InputLabel(w.parse().unwrap());
    let bar = words.iter().position(|w| *w == "|").unwrap();
    BoundaryLabels {
        boundary: n(&words[0]),
        segments: words[1..bar].iter().map(n).collect(),
        vertices: words[bar + 1..].iter().map(n).collect(),
    }
}

/// Parse one block (without its trailing `end`).
pub fn parse(block: &str) -> CaseSpec {
    let mut spec = CaseSpec {
        name: String::new(),
        tolerance: Tolerance::default(),
        operation: OperationId(0),
        frame: [0.0; 9],
        start: 0.0,
        end: 0.0,
        boundaries: Vec::new(),
        box_at: None,
        transforms: Vec::new(),
    };
    for line in block.lines().filter(|l| !l.trim().is_empty()) {
        let w: Vec<&str> = line.split_whitespace().collect();
        let f = |i: usize| -> f64 { w[i].parse().unwrap() };
        match w[0] {
            "case" => {
                spec.name = w[1].to_string();
                spec.tolerance = Tolerance::new(f(2), 1e-12).unwrap();
            }
            "op" => spec.operation = OperationId(w[1].parse().unwrap()),
            "frame" => spec.frame = std::array::from_fn(|i| f(i + 1)),
            "offsets" => (spec.start, spec.end) = (f(1), f(2)),
            "box" => spec.box_at = Some(([f(1), f(2), f(3)], [f(4), f(5), f(6)])),
            "boundary" => {
                let (boundary, rest) = if w[1] == "C" {
                    let b =
                        Boundary::circle(Point2::new(f(2), f(3)), f(4), spec.tolerance).unwrap();
                    (b, 5)
                } else {
                    let n: usize = w[2].parse().unwrap();
                    let points = (0..n)
                        .map(|k| Point2::new(f(3 + 2 * k), f(4 + 2 * k)))
                        .collect();
                    (
                        Boundary::polygon(points, spec.tolerance).unwrap(),
                        3 + 2 * n,
                    )
                };
                let boundary = if w.len() > rest && w[rest] == "labels" {
                    boundary.with_labels(labels(&w[rest + 1..])).unwrap()
                } else {
                    boundary
                };
                spec.boundaries.push(boundary);
            }
            "transform" => spec.transforms.push(if w[1] == "T" {
                RigidTransform::translation(Vec3::new(f(2), f(3), f(4))).unwrap()
            } else {
                RigidTransform::rotation(
                    Point3::new(f(2), f(3), f(4)),
                    Vec3::new(f(5), f(6), f(7)),
                    f(8),
                )
                .unwrap()
            }),
            other => panic!("unknown row {other}"),
        }
    }
    spec
}

/// Every case block of a protocol file.
pub fn cases(text: &str) -> Vec<CaseSpec> {
    text.split("\nend")
        .filter(|b| !b.trim().is_empty())
        .map(|b| parse(b.trim()))
        .collect()
}

/// The solid before any transform.
pub fn build(spec: &CaseSpec) -> Solid {
    build_tracked(spec).0
}

/// The solid before any transform, with its construction history.
pub fn build_tracked(spec: &CaseSpec) -> (Solid, History) {
    if let Some((origin, size)) = spec.box_at {
        let o = Point3::new(origin[0], origin[1], origin[2]);
        let size = Vec3::new(size[0], size[1], size[2]);
        return Solid::box_at_with(spec.operation, o, size, spec.tolerance).unwrap();
    }
    let f = spec.frame;
    let frame = Frame3::new(
        Point3::new(f[0], f[1], f[2]),
        Vec3::new(f[3], f[4], f[5]),
        Vec3::new(f[6], f[7], f[8]),
        spec.tolerance,
    )
    .unwrap();
    let profile = Profile::new(
        spec.boundaries[0].clone(),
        spec.boundaries[1..].to_vec(),
        spec.tolerance,
    )
    .unwrap();
    Solid::extrude_with(spec.operation, profile, frame, spec.start, spec.end).unwrap()
}

pub fn role_name(role: Role) -> &'static str {
    match role {
        Role::StartCap => "start_cap",
        Role::EndCap => "end_cap",
        Role::Wall => "wall",
        Role::BottomEdge => "bottom_edge",
        Role::TopEdge => "top_edge",
        Role::Vertical => "vertical",
        Role::Seam => "seam",
        Role::BottomVertex => "bottom_vertex",
        Role::TopVertex => "top_vertex",
        Role::SeamVertex => "seam_vertex",
        Role::Body => "body",
        Role::External => "external",
    }
}

pub fn parent_text(p: &Parent) -> String {
    match p {
        Parent::Label(l) => format!("L{}", l.0),
        Parent::Profile { boundary, element } => match element {
            ProfileElement::Boundary => format!("P{boundary}.b0"),
            ProfileElement::Segment(i) => format!("P{boundary}.s{i}"),
            ProfileElement::Vertex(i) => format!("P{boundary}.v{i}"),
        },
        Parent::Entity(id) => format!("E{id}"),
    }
}

/// Structural locator of every vertex, edge and face, found from geometry:
/// vertices by position against the profile's stored points at the start or
/// end offset, edges by their end vertices, faces by their edges.
pub fn rows(solid: &Solid) -> Vec<String> {
    let t = solid.topology();
    let frame = solid.frame();
    let profile = solid.profile();
    let tolerance = profile.tolerance();
    let sides = [("start", solid.start_offset()), ("end", solid.end_offset())];
    let boundaries: Vec<&Boundary> = std::iter::once(profile.outer())
        .chain(profile.holes())
        .collect();
    let scale = solid.bounds().max.distance(solid.bounds().min).max(1.0);
    let budget = 64.0 * f64::EPSILON * scale + 4.0 * tolerance.linear().min(scale * 1e-9);
    let locate_vertex = |p: Point3| -> (usize, usize, &'static str) {
        let mut found = Vec::new();
        for (b, boundary) in boundaries.iter().enumerate() {
            let points: Vec<Point2> = match boundary.polygon_vertices() {
                Some(points) => points.to_vec(),
                None => {
                    let (c, r) = boundary.circle_geometry().unwrap();
                    vec![Point2::new(c.x + r, c.y)]
                }
            };
            for (j, q) in points.iter().enumerate() {
                for (side, offset) in sides {
                    if frame.point(*q, offset).distance(p) <= budget {
                        found.push((b, j, side));
                    }
                }
            }
        }
        assert_eq!(found.len(), 1, "vertex {p:?} matches {found:?}");
        found[0]
    };
    let vertex_loc: Vec<_> = t
        .vertices()
        .iter()
        .map(|v| locate_vertex(v.position))
        .collect();
    let mut edge_loc = Vec::new();
    for edge in t.edges() {
        let (a, b) = (vertex_loc[edge.start.index()], vertex_loc[edge.end.index()]);
        let loc = match &edge.curve {
            Curve3::LineSegment { .. } if a.0 == b.0 && a.1 == b.1 && a.2 != b.2 => {
                format!("{} vertex {} both", a.0, a.1)
            }
            Curve3::LineSegment { .. } => {
                assert_eq!((a.0, a.2), (b.0, b.2), "edge between boundaries or sides");
                let n = boundaries[a.0].polygon_vertices().unwrap().len();
                let j = if (a.1 + 1) % n == b.1 { a.1 } else { b.1 };
                assert!((a.1 + 1) % n == b.1 || (b.1 + 1) % n == a.1);
                format!("{} segment {} {}", a.0, j, a.2)
            }
            _ => {
                assert_eq!(a, b);
                format!("{} segment 0 {}", a.0, a.2)
            }
        };
        edge_loc.push(loc);
    }
    let mut face_loc = Vec::new();
    for face in t.faces() {
        let edges: Vec<&String> = face
            .loops
            .iter()
            .flatten()
            .map(|u| &edge_loc[u.edge.index()])
            .collect();
        let walls: Vec<&String> = edges
            .iter()
            .copied()
            .filter(|l| l.contains(" segment "))
            .collect();
        let loc = if edges.iter().all(|l| l.contains(" segment "))
            && face.loops.len() == boundaries.len()
        {
            let side = walls[0].split(' ').next_back().unwrap();
            assert!(walls.iter().all(|l| l.ends_with(side)));
            assert!(matches!(face.surface, Surface::Plane(_)));
            format!("cap {side}")
        } else {
            let w: Vec<&str> = walls[0].split(' ').collect();
            format!("{} segment {} both", w[0], w[2])
        };
        face_loc.push(loc);
    }
    let mut out = Vec::new();
    for (id, slot) in t.ids() {
        let d = t.derivation(id).unwrap();
        assert_eq!(d.id(), id, "derivation must re-hash to its id");
        assert_eq!(t.slot_of(id), Some(slot));
        assert_eq!(t.id_of(slot), Some(id));
        let (kind, loc) = match slot {
            Slot::Vertex(v) => {
                let (b, j, side) = vertex_loc[v.index()];
                ("vertex", format!("{b} vertex {j} {side}"))
            }
            Slot::Edge(e) => ("edge", edge_loc[e.index()].clone()),
            Slot::Face(f) => ("face", face_loc[f.index()].clone()),
        };
        let parents: Vec<String> = d.parents.iter().map(parent_text).collect();
        out.push(format!(
            "{id} {kind} {} {} {} {loc}",
            role_name(d.role),
            d.ordinal,
            parents.join(",")
        ));
    }
    out.sort();
    out
}

/// One relation as fixture text (identity_reference.py::relation_text).
pub fn relation_text(r: &Relation) -> String {
    let ids = |v: &[rusty_occt::identity::EntityId]| {
        v.iter()
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(",")
    };
    match r {
        Relation::Generated { from, to, role } => {
            let parents: Vec<String> = from.iter().map(parent_text).collect();
            format!("generated {} {to} {}", parents.join(","), role_name(*role))
        }
        Relation::Modified { from, to } => format!("modified {from} {to}"),
        Relation::Unchanged { id } => format!("unchanged {id}"),
        Relation::Deleted { id } => format!("deleted {id}"),
        Relation::Split { from, into } => format!("split {from} {}", ids(into)),
        Relation::Merged { from, into } => format!("merged {} {into}", ids(from)),
    }
}
