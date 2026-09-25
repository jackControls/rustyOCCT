//! Valid prisms (star outlines, round or square holes, box cavities) mutated
//! into specific invalid shells. Each mutation predicts its issues from what
//! it changed: the exact report where the contract is local, or a required
//! issue on the mutated entity. Cavities and inversions are assembled here,
//! not by a kernel builder, so shell merging is independent of production.
use crate::byte;
use rusty_occt::identity::OperationId;
use rusty_occt::topology::{
    Curve2, Curve3, Edge, EdgeId, FaceId, Orientation, Topology, TopologyParts, Vertex, VertexId,
};
use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
use std::f64::consts::TAU;

struct Bytes<'a>(&'a [u8], usize);
impl Bytes<'_> {
    fn next(&mut self) -> u8 {
        self.1 += 1;
        byte(self.0, self.1 - 1)
    }
    /// In [0, 1].
    fn unit(&mut self) -> f64 {
        f64::from(self.next()) / 255.0
    }
    /// In [-1, 1].
    fn signed(&mut self) -> f64 {
        2.0 * self.unit() - 1.0
    }
    fn pick(&mut self, n: usize) -> usize {
        usize::from(self.next()) % n
    }
}

fn flip(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
    }
}

fn reversed(p: &Curve2) -> Curve2 {
    match p {
        Curve2::LineSegment { start, end } => Curve2::LineSegment {
            start: *end,
            end: *start,
        },
        Curve2::CircularArc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => Curve2::CircularArc {
            center: *center,
            radius: *radius,
            start_angle: start_angle + sweep_angle,
            sweep_angle: -sweep_angle,
        },
    }
}

/// Turn a shell inside out: flip each face, and traverse each loop backwards
/// with every use and pcurve reversed.
fn invert(parts: &mut TopologyParts, shell: usize) {
    for id in parts.shells[shell].clone() {
        let face = &mut parts.faces[id.index()];
        face.orientation = flip(face.orientation);
        for lp in &mut face.loops {
            lp.reverse();
            for u in lp.iter_mut() {
                u.orientation = flip(u.orientation);
                u.pcurve = reversed(&u.pcurve);
            }
        }
    }
}

fn parts_of(t: &Topology) -> TopologyParts {
    TopologyParts {
        vertices: t.vertices().to_vec(),
        edges: t.edges().to_vec(),
        faces: t.faces().to_vec(),
        shells: t.shells().to_vec(),
    }
}

/// Append `other` as further shells, renumbering its entities.
fn merge(parts: &mut TopologyParts, other: TopologyParts) {
    let (nv, ne, nf) = (parts.vertices.len(), parts.edges.len(), parts.faces.len());
    parts.vertices.extend(other.vertices);
    parts.edges.extend(other.edges.into_iter().map(|e| Edge {
        start: VertexId::new(e.start.index() + nv),
        end: VertexId::new(e.end.index() + nv),
        curve: e.curve,
    }));
    parts.faces.extend(other.faces.into_iter().map(|mut f| {
        for u in f.loops.iter_mut().flatten() {
            u.edge = EdgeId::new(u.edge.index() + ne);
        }
        f
    }));
    parts.shells.extend(
        other
            .shells
            .into_iter()
            .map(|s| s.into_iter().map(|f| FaceId::new(f.index() + nf)).collect()),
    );
}

fn shift_v(p: &Curve2, dv: f64) -> Curve2 {
    let up = |q: Point2| Point2::new(q.x, q.y + dv);
    match p {
        Curve2::LineSegment { start, end } => Curve2::LineSegment {
            start: up(*start),
            end: up(*end),
        },
        Curve2::CircularArc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => Curve2::CircularArc {
            center: up(*center),
            radius: *radius,
            start_angle: *start_angle,
            sweep_angle: *sweep_angle,
        },
    }
}

enum Expect {
    Valid,
    Exactly(Vec<String>),
    /// Every listed issue, plus at least one from each alternative group.
    Contains(Vec<String>, Vec<Vec<String>>),
}

fn report(parts: &TopologyParts, tolerance: Tolerance) -> Vec<String> {
    let issues = parts.check(tolerance);
    // Deterministic, duplicate-free, and consistent with the constructor.
    assert_eq!(issues, parts.check(tolerance));
    let mut names: Vec<String> = issues.iter().map(|i| i.to_string()).collect();
    names.sort();
    let count = names.len();
    names.dedup();
    assert_eq!(count, names.len(), "duplicate issues {names:?}");
    match Topology::from_parts(parts.clone(), tolerance) {
        Ok(_) => assert!(names.is_empty()),
        Err(e) => assert_eq!(e, issues),
    }
    names
}

fn use_at(parts: &TopologyParts, b: &mut Bytes) -> (usize, usize, usize) {
    let fi = b.pick(parts.faces.len());
    let li = b.pick(parts.faces[fi].loops.len());
    let ui = b.pick(parts.faces[fi].loops[li].len());
    (fi, li, ui)
}

pub fn check_brep_validation(data: &[u8]) {
    let mut b = Bytes(data, 0);
    let mutation = b.next() % 16;
    let count = 3 + usize::from(b.next() % 10);
    let feature = b.next() % 4;
    let scale = 2.0_f64.powi(i32::from(b.next() % 21) - 10);
    let tolerance = Tolerance::new(1e-9 * scale, 1e-12).unwrap();
    let tau = tolerance.linear();
    // A star-shaped outline around the origin always is a simple polygon.
    let outline: Vec<_> = (0..count)
        .map(|i| {
            let angle = TAU * i as f64 / count as f64;
            let radius = scale * (0.8 + 0.6 * b.unit());
            Point2::new(radius * angle.cos(), radius * angle.sin())
        })
        .collect();
    let square = |h: f64| {
        Boundary::polygon(
            vec![
                Point2::new(-h, -h),
                Point2::new(h, -h),
                Point2::new(h, h),
                Point2::new(-h, h),
            ],
            tolerance,
        )
    };
    let holes = match feature {
        1 => vec![Boundary::circle(Point2::default(), 0.2 * scale, tolerance).unwrap()],
        2 => vec![square(0.15 * scale).unwrap()],
        _ => vec![],
    };
    let Ok(outer) = Boundary::polygon(outline, tolerance) else {
        return;
    };
    let Ok(profile) = Profile::new(outer, holes, tolerance) else {
        return;
    };
    let normal = Vec3::new(b.signed(), b.signed(), 0.5 + b.unit());
    let origin = Point3::new(
        10.0 * scale * b.signed(),
        10.0 * scale * b.signed(),
        10.0 * scale * b.signed(),
    );
    let Ok(frame) = Frame3::new(origin, normal, Vec3::X, tolerance) else {
        return;
    };
    let low = -scale * (0.5 + b.unit());
    let high = scale * (0.5 + b.unit());
    let solid = Solid::extrude_with(OperationId::UNSPECIFIED, profile, frame, low, high)
        .map(|(s, _)| s)
        .expect("valid prism");
    let mut parts = parts_of(solid.topology());
    let with_cavity = feature == 3 || mutation == 11;
    if with_cavity {
        let inner = Profile::new(square(0.1 * scale).unwrap(), vec![], tolerance).unwrap();
        let cavity = Solid::extrude_with(
            OperationId::UNSPECIFIED,
            inner,
            frame,
            0.5 * low,
            0.5 * high,
        )
        .map(|(s, _)| s)
        .expect("valid box");
        let mut cavity = parts_of(cavity.topology());
        if mutation != 11 {
            invert(&mut cavity, 0);
        }
        merge(&mut parts, cavity);
    }
    if mutation != 11 {
        assert_eq!(report(&parts, tolerance), Vec::<String>::new());
    }
    let one = |s: String| vec![s];
    let expect = match mutation {
        0 => Expect::Valid,
        1 => {
            let n = parts.vertices.len();
            parts.vertices.push(Vertex {
                position: Point3::new(b.signed(), b.signed(), b.signed()),
            });
            Expect::Exactly(one(format!("unused_vertex:vertex {n}")))
        }
        2 => {
            let (a, c) = (b.pick(parts.vertices.len()), b.pick(parts.vertices.len()));
            let (p, q) = (parts.vertices[a].position, parts.vertices[c].position);
            if a == c || (q - p).length() <= 2.0 * tau {
                return;
            }
            let n = parts.edges.len();
            parts.edges.push(Edge {
                start: VertexId::new(a),
                end: VertexId::new(c),
                curve: Curve3::LineSegment { start: p, end: q },
            });
            Expect::Exactly(one(format!("unused_edge:edge {n}")))
        }
        3 => {
            parts.shells.push(vec![]);
            Expect::Exactly(one(format!("empty_shell:shell {}", parts.shells.len() - 1)))
        }
        4 => {
            let fi = b.pick(parts.faces.len());
            parts.faces[fi].loops.push(vec![]);
            let li = parts.faces[fi].loops.len() - 1;
            Expect::Exactly(one(format!("empty_loop:loop {fi}.{li}")))
        }
        5 => {
            let (fi, li, ui) = use_at(&parts, &mut b);
            let bad = parts.edges.len() + usize::from(b.next());
            parts.faces[fi].loops[li][ui].edge = EdgeId::new(bad);
            Expect::Exactly(one(format!("reference:use {fi}.{li}.{ui}")))
        }
        6 => {
            let si = b.pick(parts.shells.len());
            let k = b.pick(parts.shells[si].len());
            let fi = parts.shells[si].remove(k).index();
            Expect::Contains(one(format!("face_without_shell:face {fi}")), vec![])
        }
        7 => {
            let (fi, li, ui) = use_at(&parts, &mut b);
            let u = &mut parts.faces[fi].loops[li][ui];
            u.orientation = flip(u.orientation);
            let e = u.edge.index();
            Expect::Contains(one(format!("same_sense_uses:edge {e}")), vec![])
        }
        8 => {
            let v = b.pick(parts.vertices.len());
            let p = parts.vertices[v].position;
            parts.vertices[v].position = Point3::new(p.x + 1000.0 * tau, p.y, p.z);
            // Some incident edge end must now miss the moved vertex.
            let ends: Vec<String> = parts
                .edges
                .iter()
                .enumerate()
                .flat_map(|(i, e)| {
                    let mut out = vec![];
                    if e.start.index() == v {
                        out.push(format!("vertex_off_curve:edge {i} start"));
                    }
                    if e.end.index() == v {
                        out.push(format!("vertex_off_curve:edge {i} end"));
                    }
                    out
                })
                .collect();
            Expect::Contains(vec![], vec![ends])
        }
        9 => {
            // The v axis has unit speed on planes and cylinders alike.
            let (fi, li, ui) = use_at(&parts, &mut b);
            let u = &mut parts.faces[fi].loops[li][ui];
            u.pcurve = shift_v(&u.pcurve, 1000.0 * tau);
            Expect::Contains(one(format!("pcurve_off_edge:use {fi}.{li}.{ui}")), vec![])
        }
        10 => {
            invert(&mut parts, 0);
            Expect::Exactly(one("shell_orientation:shell 0".into()))
        }
        11 => Expect::Exactly(one("shell_orientation:shell 1".into())),
        12 => {
            let fi = b.pick(parts.faces.len());
            let f = &mut parts.faces[fi];
            f.orientation = flip(f.orientation);
            Expect::Contains(one(format!("loop_winding:loop {fi}.0")), vec![])
        }
        13 => {
            // A shift the looser tolerance absorbs and the original rejects.
            let (fi, li, ui) = use_at(&parts, &mut b);
            let u = &mut parts.faces[fi].loops[li][ui];
            u.pcurve = shift_v(&u.pcurve, 10.0 * tau);
            let loose = Tolerance::new(100.0 * tau, 1e-12).unwrap();
            assert_eq!(report(&parts, loose), Vec::<String>::new());
            Expect::Contains(one(format!("pcurve_off_edge:use {fi}.{li}.{ui}")), vec![])
        }
        14 => {
            let si = b.pick(parts.shells.len());
            let face = parts.shells[si][b.pick(parts.shells[si].len())];
            parts.shells[si].push(face);
            Expect::Contains(one(format!("face_reused:face {}", face.index())), vec![])
        }
        _ => {
            // Swap a face's outer loop with a hole: both wind the wrong way.
            let Some(fi) = parts.faces.iter().position(|f| f.loops.len() > 1) else {
                return;
            };
            parts.faces[fi].loops.swap(0, 1);
            Expect::Contains(
                vec![],
                vec![vec![
                    format!("loop_winding:loop {fi}.0"),
                    format!("loop_winding:loop {fi}.1"),
                ]],
            )
        }
    };
    let got = report(&parts, tolerance);
    match expect {
        Expect::Valid => assert_eq!(got, Vec::<String>::new()),
        Expect::Exactly(mut want) => {
            want.sort();
            assert_eq!(got, want);
        }
        Expect::Contains(all, any) => {
            assert!(!got.is_empty());
            for w in &all {
                assert!(got.contains(w), "missing {w} in {got:?}");
            }
            for group in &any {
                assert!(
                    group.iter().any(|w| got.contains(w)),
                    "none of {group:?} in {got:?}"
                );
            }
        }
    }
}
