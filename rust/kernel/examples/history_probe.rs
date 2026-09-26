//! Test-only protocol for the source-pinned history comparison. Reads case
//! blocks (tools/identity_reference.py::encode_case) and prints each
//! construction relation with its parent's profile locator and a geometric
//! signature of its target, then every transform step's before/after pairs,
//! in the signature format of rust/tools/occt_history_oracle.cpp. A region
//! is signed as its body; `C` rows are the synthesized OCCT counts.
#[path = "../tests/support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
#[path = "../tests/support/signatures.rs"]
mod signatures;
use identity_protocol::{cases, role_name};
use rusty_occt::history::Relation;
use rusty_occt::identity::{OperationId, Parent, ProfileElement};
use rusty_occt::{RigidTransform, Solid, Vec3};
use signatures::{body, counts, signature};
use std::io::Read;

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
