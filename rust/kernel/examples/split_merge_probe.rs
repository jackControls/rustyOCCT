//! Test-only protocol for the source-pinned split/fuse comparison. Reads
//! scenario blocks (tools/generate_split_merge_fixtures.py::encode) that the
//! native oracle represents, runs their split and fuse steps and prints, per
//! stage (`split`, `fuse`, and `composed` for a split whose pieces are fused
//! back), every relation with the geometric signatures of its sources and
//! targets, and every output body with its synthesized OCCT counts and its
//! entities' signatures (format of rust/tools/occt_split_merge_oracle.cpp).
#[path = "../tests/support/signatures.rs"]
#[allow(dead_code)]
mod signatures;
#[path = "../tests/support/split_merge_protocol.rs"]
#[allow(dead_code)]
mod split_merge_protocol;
use rusty_occt::history::{History, Relation};
use rusty_occt::identity::{EntityId, Parent};
use rusty_occt::Solid;
use signatures::{counts, mass, signature};
use split_merge_protocol::{scenarios, Step};
use std::collections::BTreeMap;
use std::io::Read;

fn sig(bodies: &[&Solid], id: EntityId) -> String {
    for b in bodies {
        if let Some(slot) = b.topology().slot_of(id) {
            return signature(b, slot);
        }
    }
    panic!("unknown id {id}");
}

fn stage(name: &str, inputs: &[&Solid], outputs: &[&Solid], h: &History) {
    let all: Vec<&Solid> = inputs.iter().chain(outputs).copied().collect();
    let sigs = |ids: &[EntityId]| {
        ids.iter()
            .map(|i| sig(&all, *i))
            .collect::<Vec<_>>()
            .join(" ; ")
    };
    for r in &h.relations {
        let (kind, sources, targets) = match r {
            Relation::Unchanged { id } => ("unchanged", vec![*id], vec![*id]),
            Relation::Modified { from, to } => ("modified", vec![*from], vec![*to]),
            Relation::Split { from, into } => ("split", vec![*from], into.clone()),
            Relation::Merged { from, into } => ("merged", from.clone(), vec![*into]),
            Relation::Deleted { id } => ("deleted", vec![*id], vec![]),
            Relation::Generated { from, to, role } => {
                let parents: Vec<EntityId> = from
                    .iter()
                    .map(|p| match p {
                        Parent::Entity(e) => *e,
                        other => panic!("unexpected parent {other:?}"),
                    })
                    .collect();
                println!(
                    "G {name} {role:?} | {} -> {}",
                    sigs(&parents),
                    sig(&all, *to)
                );
                continue;
            }
        };
        for s in sources {
            println!("Q {name} {kind} | {} -> {}", sig(&all, s), sigs(&targets));
        }
    }
    for (k, b) in outputs.iter().enumerate() {
        println!("B {name} {k} {}", mass(b));
        println!("{}", counts(b).replacen("C", &format!("C {name} {k}"), 1));
        for (id, _) in b.topology().ids() {
            println!("O {name} {}", sig(&[b], id));
        }
    }
}

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for s in scenarios(&input) {
        println!("R {}", s.name);
        let mut bodies: BTreeMap<String, Solid> = s.bodies.iter().cloned().collect();
        let mut histories = Vec::new();
        for step in &s.steps {
            match step {
                Step::Split {
                    body,
                    height,
                    operation,
                    lower,
                    upper,
                } => {
                    let ([l, u], h) = bodies[body].split_at_height(*operation, *height).unwrap();
                    stage("split", &[&bodies[body]], &[&l, &u], &h);
                    bodies.insert(lower.clone(), l);
                    bodies.insert(upper.clone(), u);
                    histories.push((vec![body.clone()], vec![lower.clone(), upper.clone()], h));
                }
                Step::Fuse {
                    a,
                    b,
                    operation,
                    out,
                } => {
                    let (f, h) = bodies[a].fuse_stacked(&bodies[b], *operation).unwrap();
                    stage("fuse", &[&bodies[a], &bodies[b]], &[&f], &h);
                    bodies.insert(out.clone(), f);
                    histories.push((vec![a.clone(), b.clone()], vec![out.clone()], h));
                }
                Step::Compose(i, j) => {
                    let composed = histories[*i].2.then(&histories[*j].2).unwrap();
                    let inputs: Vec<&Solid> = histories[*i].0.iter().map(|n| &bodies[n]).collect();
                    let outputs: Vec<&Solid> = histories[*j].1.iter().map(|n| &bodies[n]).collect();
                    stage("composed", &inputs, &outputs, &composed);
                    histories.push((histories[*i].0.clone(), histories[*j].1.clone(), composed));
                }
            }
        }
        println!("end");
    }
}
