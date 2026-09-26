//! Valid prisms (star outlines, round or square holes, box cavities) mutated
//! into specific invalid cell complexes. Each mutation predicts its issues
//! from what it changed: the exact report where the contract is local, or a
//! required issue on the mutated entity. Cavities and inversions are
//! assembled here, not by a kernel builder, so shell and region merging is
//! independent of production. The cell model's own failure modes
//! (TOPOLOGY_MODEL.md) are mutations 16-21; the period shift needs a loop of
//! several fins on a cylinder, which prisms of polygons and circles do not
//! have, so it mutates the valid stadium fixtures. Enclosures (M5) are
//! mutations 24-26: a missing bound, a bound outside [0, resolution], and a
//! vertex moved within the resolution but past its measured bound; mutation
//! 27 places a circle prism far out at a tolerance near the coordinates'
//! resolution, where construction must fail or enclose within it. Mutation
//! 28 builds a cone or frustum (S3 of REVIEW_NOTES.md) and moves its pole
//! along a ruling or off the surface, drops the pole, makes the surface
//! degenerate or shifts a ring's pcurve, as the cone fixtures do.
use crate::byte;
use rusty_occt::identity::OperationId;
use rusty_occt::topology::{
    Curve2, Curve3, Edge, EdgeId, Enclosure, FaceId, FinId, Loop, LoopId, Orientation, Region,
    RegionId, RegionKind, Shell, ShellId, Surface, Topology, TopologyParts, Vertex, VertexId,
};
use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
use std::f64::consts::TAU;

#[path = "../../kernel/tests/support/brep_protocol.rs"]
mod brep_protocol;

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

/// Turn a shell's faces inside out: flip each face, and traverse each loop
/// backwards with every fin and pcurve reversed and its winding negated.
/// Sides and regions are left as they were. A reversed arc's start angle
/// rounds, so the changed fins and faces are measured again, as a builder
/// would measure its output.
fn invert(parts: &mut TopologyParts, shell: usize) {
    let faces: Vec<FaceId> = parts.shells[shell].sides.iter().map(|(f, _)| *f).collect();
    for id in faces {
        let face = &mut parts.faces[id.index()];
        face.sense = flip(face.sense);
        face.enclosure = None;
        for l in face.loops.clone() {
            if let Loop::Edges { fins, winding } = &mut parts.loops[l.index()] {
                fins.reverse();
                winding[0] = -winding[0];
                for k in fins.iter() {
                    let fin = &mut parts.fins[k.index()];
                    fin.sense = flip(fin.sense);
                    fin.pcurve = reversed(&fin.pcurve);
                    fin.enclosure = None;
                }
            }
        }
    }
    *parts = std::mem::take(parts).with_measured_enclosures();
}

fn parts_of(t: &Topology) -> TopologyParts {
    TopologyParts {
        vertices: t.vertices().to_vec(),
        edges: t.edges().to_vec(),
        fins: t.fins().to_vec(),
        loops: t.loops().to_vec(),
        faces: t.faces().to_vec(),
        shells: t.shells().to_vec(),
        regions: t.regions().to_vec(),
    }
}

/// Append a prism's complex as a cavity: its material shell (0) bounds the
/// body's solid region as a further shell, its twin (1) a new bounded void.
fn merge_cavity(parts: &mut TopologyParts, other: TopologyParts) {
    let (nv, ne, nk) = (parts.vertices.len(), parts.edges.len(), parts.fins.len());
    let (nl, nf, ns) = (parts.loops.len(), parts.faces.len(), parts.shells.len());
    let nr = parts.regions.len();
    let v = |v: VertexId| VertexId::new(v.index() + nv);
    let k = |k: &FinId| FinId::new(k.index() + nk);
    parts.vertices.extend(other.vertices);
    parts.edges.extend(other.edges.into_iter().map(|e| Edge {
        start: e.start.map(v),
        end: e.end.map(v),
        curve: e.curve,
        fins: e.fins.iter().map(k).collect(),
    }));
    parts.fins.extend(other.fins.into_iter().map(|mut f| {
        f.edge = EdgeId::new(f.edge.index() + ne);
        f
    }));
    parts.loops.extend(other.loops.into_iter().map(|l| match l {
        Loop::Edges { fins, winding } => Loop::Edges {
            fins: fins.iter().map(k).collect(),
            winding,
        },
        Loop::Vertex(x) => Loop::Vertex(v(x)),
    }));
    parts.faces.extend(other.faces.into_iter().map(|mut f| {
        f.loops = f
            .loops
            .iter()
            .map(|l| LoopId::new(l.index() + nl))
            .collect();
        f.front = ShellId::new(f.front.index() + ns);
        f.back = ShellId::new(f.back.index() + ns);
        f
    }));
    for (i, shell) in other.shells.into_iter().enumerate() {
        parts.shells.push(Shell {
            region: RegionId::new(if i == 0 { 1 } else { nr }),
            sides: shell
                .sides
                .iter()
                .map(|(f, side)| (FaceId::new(f.index() + nf), *side))
                .collect(),
            wire_edges: shell
                .wire_edges
                .iter()
                .map(|e| EdgeId::new(e.index() + ne))
                .collect(),
            acorn_vertices: shell.acorn_vertices.iter().map(|x| v(*x)).collect(),
        });
    }
    parts.regions[1].shells.push(ShellId::new(ns));
    parts.regions.push(Region {
        kind: RegionKind::Void,
        shells: vec![ShellId::new(ns + 1)],
    });
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

/// A face, a loop position in it, a fin position in that loop, and the fin.
fn fin_at(parts: &TopologyParts, b: &mut Bytes) -> Option<(usize, usize, usize, FinId)> {
    let fi = b.pick(parts.faces.len());
    let li = b.pick(parts.faces[fi].loops.len());
    let Loop::Edges { fins, .. } = &parts.loops[parts.faces[fi].loops[li].index()] else {
        return None;
    };
    let ui = b.pick(fins.len());
    Some((fi, li, ui, fins[ui]))
}

/// The stadium fixtures: prisms whose walls are cylinder patches bounded by
/// loops of four fins.
fn stadium(b: &mut Bytes) -> (TopologyParts, Tolerance) {
    let text = include_str!("../../fixtures/brep-cases.txt");
    let blocks: Vec<&str> = text
        .split("\nend")
        .filter(|x| {
            let name = x.trim().lines().next().unwrap_or("");
            name == "case stadium" || name == "case stadium_rotated"
        })
        .collect();
    let (_, tolerance, parts) = brep_protocol::parse(blocks[b.pick(blocks.len())].trim());
    (parts, Tolerance::new(tolerance, 1e-12).unwrap())
}

/// A circle prism far from the origin, at a tolerance a few ulps of its
/// coordinates: the builder either fails or returns a valid body whose every
/// enclosure fits the resolution. It never widens a bound to succeed.
fn far_prism(b: &mut Bytes) {
    let far = 2.0_f64.powi(i32::from(b.next() % 48));
    let ulps = f64::from(1 + b.next() % 64);
    let tolerance = match Tolerance::new(far * f64::EPSILON * ulps, 1e-12) {
        Ok(t) => t,
        Err(_) => return,
    };
    let radius = 1.0 + b.unit();
    let Ok(circle) = Boundary::circle(Point2::new(b.signed(), b.signed()), radius, tolerance)
    else {
        return;
    };
    let Ok(profile) = Profile::new(circle, vec![], tolerance) else {
        return;
    };
    let origin = Point3::new(far * (1.0 + b.unit()), far * b.signed(), far * b.signed());
    let normal = Vec3::new(b.signed(), b.signed(), 0.5 + b.unit());
    let Ok(frame) = Frame3::new(origin, normal, Vec3::X, tolerance) else {
        return;
    };
    let Ok((solid, _)) = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile,
        frame,
        0.0,
        1.0 + b.unit(),
    ) else {
        return;
    };
    let t = solid.topology();
    assert!(t.check(tolerance).is_empty());
    let bounds = t
        .vertices()
        .iter()
        .map(|v| v.enclosure)
        .chain(t.fins().iter().map(|f| f.enclosure))
        .chain(t.faces().iter().map(|f| f.enclosure));
    for e in bounds {
        let e = e.expect("a builder encloses every entity");
        assert!(e.bound >= 0.0 && e.bound <= tolerance.linear(), "{e:?}");
    }
}

/// A cone or frustum from `Solid::cone_with`, valid with enclosures within
/// the resolution and certified mass properties, then one mutation with its
/// predicted issues.
fn cone(b: &mut Bytes) {
    let scale = 2.0_f64.powi(i32::from(b.next() % 21) - 10);
    let tolerance = Tolerance::new(1e-9 * scale, 1e-12).unwrap();
    let tau = tolerance.linear();
    let radius = |b: &mut Bytes| {
        if b.next() % 3 == 0 {
            0.0
        } else {
            scale * (0.1 + b.unit())
        }
    };
    let (bottom, top, height) = (radius(b), radius(b), scale * (0.1 + b.unit()));
    let normal = Vec3::new(b.signed(), b.signed(), 0.5 + b.unit());
    let origin = Point3::new(
        10.0 * scale * b.signed(),
        10.0 * scale * b.signed(),
        10.0 * scale * b.signed(),
    );
    let Ok(frame) = Frame3::new(origin, normal, Vec3::X, tolerance) else {
        return;
    };
    let Ok((solid, _)) = Solid::cone_with(
        OperationId::UNSPECIFIED,
        frame,
        bottom,
        top,
        height,
        tolerance,
    ) else {
        // Equal radii are a cylinder, two apices no solid.
        assert!(bottom == top, "{bottom} {top} {height}");
        return;
    };
    let t = solid.topology();
    assert!(t.check(tolerance).is_empty());
    for e in t
        .vertices()
        .iter()
        .map(|v| v.enclosure)
        .chain(t.fins().iter().map(|f| f.enclosure))
        .chain(t.faces().iter().map(|f| f.enclosure))
    {
        let e = e.expect("a builder encloses every entity");
        assert!(e.bound >= 0.0 && e.bound <= tau, "{e:?}");
    }
    // The certified volume holds the frustum's, pi h (R^2 + R r + r^2) / 3,
    // up to the rounding of the stored slant and angle.
    let m = t.mass_enclosure().expect("certified cone mass properties");
    let exact = std::f64::consts::PI * height * (bottom * bottom + bottom * top + top * top) / 3.0;
    assert!(m.volume[0] <= m.volume[1]);
    assert!(
        (0.5 * (m.volume[0] + m.volume[1]) - exact).abs() <= 1e-9 * exact,
        "{m:?} {exact}"
    );
    let apices = usize::from(bottom == 0.0) + usize::from(top == 0.0);
    let c = t.occt_counts();
    assert_eq!(
        (c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids),
        (2, 3, 3 - apices, 3 - apices, 1, 1)
    );
    let mut parts = parts_of(t);
    let fi = (0..parts.faces.len())
        .find(|f| matches!(parts.faces[*f].surface, Surface::Cone { .. }))
        .expect("a wall");
    let loops = parts.faces[fi].loops.clone();
    let pole = loops
        .iter()
        .position(|l| matches!(parts.loops[l.index()], Loop::Vertex(_)));
    let ring = loops
        .iter()
        .position(|l| matches!(parts.loops[l.index()], Loop::Edges { .. }))
        .expect("a ring loop");
    let axis = frame.normal();
    let on_axis = |z: f64| frame.point(Point2::default(), z);
    let expect = match (b.next() % 6, pole) {
        (1, Some(li)) => {
            // Along a ruling towards the other end: on the surface, off
            // the apex.
            let Loop::Vertex(v) = parts.loops[loops[li].index()] else {
                unreachable!("the pole is a vertex loop")
            };
            let at = parts.vertices[v.index()].position;
            let (rim, z) = if at.distance(on_axis(0.0)) < at.distance(on_axis(height)) {
                (top, height)
            } else {
                (bottom, 0.0)
            };
            let angle = TAU * b.unit();
            let target =
                on_axis(z) + (frame.x() * angle.cos() + axis.cross(frame.x()) * angle.sin()) * rim;
            let Ok(d) = (target - at).normalized() else {
                return;
            };
            // Declared at the resolution, as the fixtures declare a moved
            // entity: the surface gap is rounding, the apex 1000 away.
            parts.vertices[v.index()].position = at + d * (1000.0 * tau);
            parts.vertices[v.index()].enclosure = Some(Enclosure::computed(tau));
            Expect::Exactly(vec![format!("pole_off_apex:loop {fi}.{li}")])
        }
        (2, Some(li)) => {
            // Along the axis: off the surface.
            let Loop::Vertex(v) = parts.loops[loops[li].index()] else {
                unreachable!("the pole is a vertex loop")
            };
            let vertex = &mut parts.vertices[v.index()];
            vertex.position = vertex.position + axis * (1000.0 * tau);
            vertex.enclosure = Some(Enclosure::computed(tau));
            Expect::Contains(
                vec![format!("vertex_loop_off_surface:loop {fi}.{li}")],
                vec![],
            )
        }
        (3, Some(li)) => {
            let slot = loops[li].index();
            parts.faces[fi].loops.remove(li);
            let ring = if ring > li { ring - 1 } else { ring };
            Expect::Contains(
                vec![
                    format!("loop_without_face:loop slot {slot}"),
                    format!("winding_mismatch:loop {fi}.{ring}"),
                ],
                vec![],
            )
        }
        (4, _) => {
            let Surface::Cone { half_angle, .. } = &mut parts.faces[fi].surface else {
                unreachable!("the wall is a cone")
            };
            *half_angle = std::f64::consts::FRAC_PI_2;
            Expect::Contains(vec![format!("degenerate_surface:face {fi}")], vec![])
        }
        (5, _) => {
            let Loop::Edges { fins, .. } = &parts.loops[loops[ring].index()] else {
                unreachable!("a ring loop")
            };
            let fin = &mut parts.fins[fins[0].index()];
            fin.pcurve = shift_v(&fin.pcurve, 1000.0 * tau);
            Expect::Contains(vec![format!("pcurve_off_edge:use {fi}.{ring}.0")], vec![])
        }
        _ => Expect::Valid,
    };
    verify(report(&parts, tolerance), expect);
}

fn verify(got: Vec<String>, expect: Expect) {
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

pub fn check_brep_validation(data: &[u8]) {
    let mut b = Bytes(data, 0);
    let mutation = b.next() % 29;
    if mutation == 27 {
        far_prism(&mut b);
        return;
    }
    if mutation == 28 {
        cone(&mut b);
        return;
    }
    if mutation == 16 {
        // A period shift of one fin of a multi-fin loop on a cylinder opens
        // the loop at both of its ends.
        let (mut parts, tolerance) = stadium(&mut b);
        assert_eq!(report(&parts, tolerance), Vec::<String>::new());
        let walls: Vec<usize> = (0..parts.faces.len())
            .filter(|f| matches!(parts.faces[*f].surface, Surface::Cylinder { .. }))
            .collect();
        let fi = walls[b.pick(walls.len())];
        let Loop::Edges { fins, .. } = &parts.loops[parts.faces[fi].loops[0].index()] else {
            unreachable!("stadium walls have edge loops");
        };
        let n = fins.len();
        let ui = b.pick(n);
        let turns = f64::from(1 + b.next() % 3) * if b.next() % 2 == 0 { 1.0 } else { -1.0 };
        let fin = &mut parts.fins[fins[ui].index()];
        let Curve2::LineSegment { start, end } = fin.pcurve else {
            unreachable!("stadium walls have line pcurves");
        };
        let du = turns * TAU;
        fin.pcurve = Curve2::LineSegment {
            start: Point2::new(start.x + du, start.y),
            end: Point2::new(end.x + du, end.y),
        };
        let got = report(&parts, tolerance);
        for want in [
            format!("uv_gap:use {fi}.0.{ui}"),
            format!("uv_gap:use {fi}.0.{}", (ui + n - 1) % n),
        ] {
            assert!(got.contains(&want), "missing {want} in {got:?}");
        }
        return;
    }
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
        merge_cavity(&mut parts, cavity);
    }
    if mutation != 11 {
        assert_eq!(report(&parts, tolerance), Vec::<String>::new());
    }
    let one = |s: String| vec![s];
    // A point of the profile clear of every hole and of the outline.
    let clear = {
        let angle = TAU * b.unit();
        Point2::new(0.3 * scale * angle.cos(), 0.3 * scale * angle.sin())
    };
    let expect = match mutation {
        0 => Expect::Valid,
        1 => {
            let n = parts.vertices.len();
            parts.vertices.push(Vertex {
                position: Point3::new(b.signed(), b.signed(), b.signed()),
                enclosure: Some(Enclosure::computed(tau)),
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
                start: Some(VertexId::new(a)),
                end: Some(VertexId::new(c)),
                curve: Curve3::LineSegment { start: p, end: q },
                fins: vec![],
            });
            Expect::Exactly(one(format!("unused_edge:edge {n}")))
        }
        3 => {
            let n = parts.shells.len();
            parts.shells.push(Shell {
                region: RegionId::new(1),
                sides: vec![],
                wire_edges: vec![],
                acorn_vertices: vec![],
            });
            parts.regions[1].shells.push(ShellId::new(n));
            Expect::Exactly(one(format!("empty_shell:shell {n}")))
        }
        4 => {
            let fi = b.pick(parts.faces.len());
            parts.loops.push(Loop::Edges {
                fins: vec![],
                winding: [0, 0],
            });
            parts.faces[fi]
                .loops
                .push(LoopId::new(parts.loops.len() - 1));
            let li = parts.faces[fi].loops.len() - 1;
            Expect::Exactly(one(format!("empty_loop:loop {fi}.{li}")))
        }
        5 => {
            let Some((fi, li, ui, k)) = fin_at(&parts, &mut b) else {
                return;
            };
            let bad = parts.edges.len() + usize::from(b.next());
            parts.fins[k.index()].edge = EdgeId::new(bad);
            Expect::Exactly(one(format!("reference:use {fi}.{li}.{ui}")))
        }
        6 => {
            // Both sides of a face leave their shells.
            let fi = b.pick(parts.faces.len());
            for shell in &mut parts.shells {
                shell.sides.retain(|(f, _)| f.index() != fi);
            }
            Expect::Contains(one(format!("face_without_shell:face {fi}")), vec![])
        }
        7 => {
            let Some((_, _, _, k)) = fin_at(&parts, &mut b) else {
                return;
            };
            let fin = &mut parts.fins[k.index()];
            fin.sense = flip(fin.sense);
            let e = fin.edge.index();
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
                    if e.start == Some(VertexId::new(v)) {
                        out.push(format!("vertex_off_curve:edge {i} start"));
                    }
                    if e.end == Some(VertexId::new(v)) {
                        out.push(format!("vertex_off_curve:edge {i} end"));
                    }
                    out
                })
                .collect();
            Expect::Contains(vec![], vec![ends])
        }
        9 => {
            // The v axis has unit speed on planes and cylinders alike.
            let Some((fi, li, ui, k)) = fin_at(&parts, &mut b) else {
                return;
            };
            let fin = &mut parts.fins[k.index()];
            fin.pcurve = shift_v(&fin.pcurve, 1000.0 * tau);
            Expect::Contains(one(format!("pcurve_off_edge:use {fi}.{li}.{ui}")), vec![])
        }
        10 => {
            invert(&mut parts, 0);
            Expect::Exactly(one("shell_orientation:shell 0".into()))
        }
        // The cavity's material shell was merged without turning it inside out.
        11 => Expect::Exactly(one("shell_orientation:shell 2".into())),
        12 => {
            let fi = b.pick(parts.faces.len());
            let f = &mut parts.faces[fi];
            f.sense = flip(f.sense);
            Expect::Contains(one(format!("loop_winding:loop {fi}.0")), vec![])
        }
        13 => {
            // A shift the looser tolerance absorbs and the original rejects.
            let Some((fi, li, ui, k)) = fin_at(&parts, &mut b) else {
                return;
            };
            let fin = &mut parts.fins[k.index()];
            fin.pcurve = shift_v(&fin.pcurve, 10.0 * tau);
            // The moved fin and its face get bounds measured on the moved
            // geometry; under the original tolerance they exceed it.
            fin.enclosure = None;
            parts.faces[fi].enclosure = None;
            parts = parts.with_measured_enclosures();
            let loose = Tolerance::new(100.0 * tau, 1e-12).unwrap();
            assert_eq!(report(&parts, loose), Vec::<String>::new());
            Expect::Contains(one(format!("pcurve_off_edge:use {fi}.{li}.{ui}")), vec![])
        }
        14 => {
            let si = b.pick(parts.shells.len());
            let side = parts.shells[si].sides[b.pick(parts.shells[si].sides.len())];
            parts.shells[si].sides.push(side);
            Expect::Contains(one(format!("face_reused:face {}", side.0.index())), vec![])
        }
        15 => {
            // Swap a planar face's outer loop with a hole: both wind the wrong way.
            let Some(fi) = parts
                .faces
                .iter()
                .position(|f| f.loops.len() > 1 && matches!(f.surface, Surface::Plane(_)))
            else {
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
        17 => {
            // Flip the winding of a ring loop: it no longer closes on the cover.
            let Some(fi) = parts
                .faces
                .iter()
                .position(|f| matches!(f.surface, Surface::Cylinder { .. }))
            else {
                return;
            };
            let li = b.pick(parts.faces[fi].loops.len());
            let Loop::Edges { winding, .. } = &mut parts.loops[parts.faces[fi].loops[li].index()]
            else {
                unreachable!("walls have edge loops");
            };
            winding[0] = -winding[0];
            Expect::Contains(
                vec![
                    format!("uv_gap:use {fi}.{li}.0"),
                    format!("winding_mismatch:loop {fi}.0"),
                ],
                vec![],
            )
        }
        18 => {
            // Two edges exchange their fin lists: neither lists its own fins.
            let (a, c) = (b.pick(parts.edges.len()), b.pick(parts.edges.len()));
            if a == c {
                return;
            }
            let fins = std::mem::take(&mut parts.edges[a].fins);
            parts.edges[a].fins = std::mem::replace(&mut parts.edges[c].fins, fins);
            Expect::Contains(
                vec![
                    format!("edge_fins_mismatch:edge {a}"),
                    format!("edge_fins_mismatch:edge {c}"),
                ],
                vec![],
            )
        }
        19 => {
            // A face's front side is recorded at its back shell.
            let fi = b.pick(parts.faces.len());
            let f = &mut parts.faces[fi];
            f.front = f.back;
            Expect::Contains(one(format!("side_region_mismatch:face {fi}")), vec![])
        }
        20 | 22 => {
            // A vertex loop on a cap, or off it along the cap normal.
            let fi = b.pick(2);
            let Surface::Plane(plane) = parts.faces[fi].surface else {
                unreachable!("faces 0 and 1 are caps");
            };
            let height = frame.coordinates(plane.origin())[2];
            let lift = if mutation == 20 { 1000.0 * tau } else { 0.0 };
            let n = parts.vertices.len();
            // One resolution: sound on the cap, and not what an off-surface
            // vertex loop is reported for.
            parts.vertices.push(Vertex {
                position: frame.point(clear, height) + plane.normal() * lift,
                enclosure: Some(Enclosure::computed(tau)),
            });
            parts.loops.push(Loop::Vertex(VertexId::new(n)));
            parts.faces[fi]
                .loops
                .push(LoopId::new(parts.loops.len() - 1));
            let li = parts.faces[fi].loops.len() - 1;
            if mutation == 22 {
                Expect::Valid
            } else {
                Expect::Exactly(one(format!("vertex_loop_off_surface:loop {fi}.{li}")))
            }
        }
        21 => {
            // An open face loses its only loop.
            let Some(fi) = (0..parts.faces.len()).find(|f| parts.faces[*f].loops.len() == 1) else {
                return;
            };
            let l = parts.faces[fi].loops.pop().unwrap();
            Expect::Contains(
                vec![
                    format!("empty_face:face {fi}"),
                    format!("loop_without_face:loop slot {}", l.index()),
                ],
                vec![],
            )
        }
        24 | 25 => {
            // Remove one bound, or set it outside [0, resolution].
            let bad = [2.0 * tau, -tau, f64::NAN, f64::INFINITY][b.pick(4)];
            let (slot, entity) = match b.pick(3) {
                0 => {
                    let v = b.pick(parts.vertices.len());
                    (&mut parts.vertices[v].enclosure, format!("vertex {v}"))
                }
                1 => {
                    let f = b.pick(parts.faces.len());
                    (&mut parts.faces[f].enclosure, format!("face {f}"))
                }
                _ => {
                    let Some((fi, li, ui, k)) = fin_at(&parts, &mut b) else {
                        return;
                    };
                    (
                        &mut parts.fins[k.index()].enclosure,
                        format!("use {fi}.{li}.{ui}"),
                    )
                }
            };
            if mutation == 24 {
                *slot = None;
                Expect::Exactly(one(format!("enclosure_missing:{entity}")))
            } else {
                *slot = Some(Enclosure::computed(bad));
                Expect::Exactly(one(format!("enclosure_exceeds_resolution:{entity}")))
            }
        }
        26 => {
            // A vertex moved half the resolution: still on its curves' ends
            // within the resolution, but past its measured bound.
            let v = b.pick(parts.vertices.len());
            let d = Vec3::new(b.signed(), b.signed(), b.signed());
            let Ok(d) = d.normalized() else {
                return;
            };
            parts.vertices[v].position = parts.vertices[v].position + d * (0.5 * tau);
            Expect::Exactly(one(format!("enclosure_unsound:vertex {v}")))
        }
        _ => {
            // Radial order of a two-fin edge is cyclic: reversing it is no change.
            let e = b.pick(parts.edges.len());
            parts.edges[e].fins.reverse();
            Expect::Valid
        }
    };
    verify(report(&parts, tolerance), expect);
}
