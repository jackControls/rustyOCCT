//! Test-only protocol for the source-pinned history comparison. Reads case
//! blocks (tools/identity_reference.py::encode_case) and prints each
//! construction relation with its parent's profile locator and a geometric
//! signature of its target, then every transform step's before/after pairs,
//! in the signature format of rust/tools/occt_history_oracle.cpp. A region
//! is signed as its body; `C` rows are the synthesized OCCT counts.
#[path = "../tests/support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use identity_protocol::{cases, role_name};
use rusty_occt::history::Relation;
use rusty_occt::identity::{OperationId, Parent, ProfileElement};
use rusty_occt::topology::{Curve2, Curve3, Slot, Surface};
use rusty_occt::{Point2, Point3, RigidTransform, Solid, Vec3};
use std::io::Read;

fn p3(p: Point3) -> String {
    format!("{:?} {:?} {:?}", p.x, p.y, p.z)
}

/// Twice the signed area and twice the first moments of one pcurve,
/// from the closed forms of 1/2 ∮ (u dv - v du), 1/2 ∮ u² dv, -1/2 ∮ v² du.
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
            let cube_s = |s: f64| s - s * s * s / 3.0;
            let cube_c = |c: f64| -c + c * c * c / 3.0;
            let sq_c = |t: f64| t / 2.0 + (2.0 * t).sin() / 4.0;
            let sq_s = |t: f64| t / 2.0 - (2.0 * t).sin() / 4.0;
            [
                r * (c.x * (s1 - s0) - c.y * (c1 - c0)) + r * r * sweep,
                r * (c.x * c.x * (s1 - s0)
                    + 2.0 * c.x * r * (sq_c(t1) - sq_c(*t0))
                    + r * r * (cube_s(s1) - cube_s(s0))),
                r * (-c.y * c.y * (c1 - c0)
                    + 2.0 * c.y * r * (sq_s(t1) - sq_s(*t0))
                    + r * r * (cube_c(c1) - cube_c(c0))),
            ]
        }
    }
}

fn signature(s: &Solid, slot: Slot) -> String {
    let t = s.topology();
    match slot {
        Slot::Vertex(v) => format!("V {}", p3(t.vertices()[v.index()].position)),
        Slot::Edge(e) => {
            let c = &t.edges()[e.index()].curve;
            let kind = if matches!(c, Curve3::LineSegment { .. }) {
                "line"
            } else {
                "circle"
            };
            format!(
                "E {kind} {} {} {}",
                p3(c.point(0.0)),
                p3(c.point(0.5)),
                p3(c.point(1.0))
            )
        }
        Slot::Face(f) => {
            let face = &t.faces()[f.index()];
            let fins = t.face_fins(f);
            let fins = fins.iter().flatten();
            match face.surface {
                Surface::Plane(frame) => {
                    let mut m = [0.0; 3];
                    for u in fins {
                        let x = moments(&u.pcurve);
                        (0..3).for_each(|i| m[i] += x[i]);
                    }
                    let (area, cu, cv) = (0.5 * m[0], m[1] / m[0], m[2] / m[0]);
                    let c = frame.point(Point2::new(cu, cv), 0.0);
                    format!("F plane {:?} {}", area.abs(), p3(c))
                }
                Surface::Cylinder { frame, radius } => {
                    // A full-turn wall of two ring loops: its area is the
                    // radius times the periodic area -∮ v du, its centroid on
                    // the axis at mid height.
                    let (mut periodic, mut v) = (0.0, (f64::MAX, f64::MIN));
                    for u in fins {
                        let Curve2::LineSegment { start: a, end: b } = u.pcurve else {
                            unreachable!("prism walls have line pcurves");
                        };
                        periodic -= 0.5 * (a.y + b.y) * (b.x - a.x);
                        v = (v.0.min(a.y), v.1.max(a.y));
                    }
                    let c = frame.point(Point2::default(), 0.5 * (v.0 + v.1));
                    format!("F cylinder {:?} {}", radius * periodic.abs(), p3(c))
                }
            }
        }
        // A prism has one solid region: the body.
        Slot::Region(_) => mass(s),
    }
}

fn mass(s: &Solid) -> String {
    let m = s.mass_properties();
    format!("S {:?} {}", m.volume, p3(m.centroid))
}

fn body(s: &Solid) -> String {
    format!("B {}", mass(s))
}

/// The counts OCCT reports for the same body, seams synthesized.
fn counts(s: &Solid) -> String {
    let c = s.topology().occt_counts();
    format!(
        "C {} {} {} {} {} {}",
        c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids
    )
}

fn locator(s: &Solid, parents: &[Parent]) -> String {
    let boundaries: Vec<_> = std::iter::once(s.profile().outer())
        .chain(s.profile().holes())
        .collect();
    if parents.len() > 1
        || matches!(
            parents[0],
            Parent::Profile {
                element: ProfileElement::Boundary,
                ..
            }
        )
    {
        return "face 0 0".into();
    }
    let (b, element) = match parents[0] {
        Parent::Profile { boundary, element } => (boundary as usize, element),
        Parent::Label(l) => boundaries
            .iter()
            .enumerate()
            .find_map(|(b, x)| {
                let labels = x.labels()?;
                if labels.boundary == l {
                    return Some((b, ProfileElement::Boundary));
                }
                if let Some(j) = labels.segments.iter().position(|y| *y == l) {
                    return Some((b, ProfileElement::Segment(j as u32)));
                }
                let j = labels.vertices.iter().position(|y| *y == l)?;
                Some((b, ProfileElement::Vertex(j as u32)))
            })
            .unwrap(),
        Parent::Entity(_) => unreachable!("constructions have no entity parents"),
    };
    match element {
        ProfileElement::Boundary => "face 0 0".into(),
        ProfileElement::Segment(j) => format!("{b} edge {j}"),
        ProfileElement::Vertex(j) => format!("{b} vertex {j}"),
    }
}

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for spec in cases(&input) {
        // box_at is a cuboid construction then a translation; report both
        // steps, as the native probe does.
        let (mut solid, construct, mut transforms) = match spec.box_at {
            Some((origin, size)) => {
                let (s, h) = Solid::cuboid_with(
                    spec.operation,
                    size[0].abs(),
                    size[1].abs(),
                    size[2].abs(),
                    spec.tolerance,
                )
                .unwrap();
                let corner = Vec3::new(
                    origin[0] + size[0].min(0.0),
                    origin[1] + size[1].min(0.0),
                    origin[2] + size[2].min(0.0),
                );
                (s, h, vec![RigidTransform::translation(corner).unwrap()])
            }
            None => {
                let (s, h) = identity_protocol::build_tracked(&spec);
                (s, h, Vec::new())
            }
        };
        transforms.extend(spec.transforms.iter().copied());
        println!("R {}", spec.name);
        let t = solid.topology();
        for r in &construct.relations {
            let Relation::Generated { from, to, role } = r else {
                unreachable!("a construction only generates");
            };
            let d = t.derivation(*to).unwrap();
            println!(
                "G {} {} {} | {}",
                role_name(*role),
                d.ordinal,
                locator(&solid, from),
                signature(&solid, t.slot_of(*to).unwrap())
            );
        }
        println!("{}", body(&solid));
        println!("{}", counts(&solid));
        for (k, transform) in transforms.iter().enumerate() {
            let (next, h) = solid
                .transform_with(OperationId(k as u64), *transform)
                .unwrap();
            for r in &h.relations {
                let Relation::Modified { from, to } = r else {
                    unreachable!("a rigid motion only modifies");
                };
                let before = signature(&solid, solid.topology().slot_of(*from).unwrap());
                let after = signature(&next, next.topology().slot_of(*to).unwrap());
                println!("T {k} {before} -> {after}");
            }
            println!("{}", body(&next));
            solid = next;
        }
        println!("end");
    }
}
