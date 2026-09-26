//! Scenario protocol written by rust/tools/generate_split_merge_fixtures.py:
//! named prisms in the identity case protocol, then `split`, `fuse` and
//! `compose` steps. `evaluate` runs a scenario and returns its rows in the
//! fixture's order and format.
#[path = "identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
pub use identity_protocol::{relation_text, rows};
use rusty_occt::history::{check, History};
use rusty_occt::identity::OperationId;
use rusty_occt::Solid;
use std::collections::BTreeMap;

pub enum Step {
    Split {
        body: String,
        height: f64,
        operation: OperationId,
        lower: String,
        upper: String,
    },
    Fuse {
        a: String,
        b: String,
        operation: OperationId,
        out: String,
    },
    Compose(usize, usize),
}

pub struct Scenario {
    pub name: String,
    pub bodies: Vec<(String, Solid)>,
    pub steps: Vec<Step>,
}

pub fn parse(block: &str) -> Scenario {
    let mut lines = block.lines();
    let head: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    let (name, tolerance) = (head[1].to_string(), head[2]);
    let mut sections: Vec<(String, Vec<String>)> = Vec::new();
    let mut steps = Vec::new();
    for line in lines {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w[0] {
            "body" => sections.push((w[1].to_string(), Vec::new())),
            "split" => steps.push(Step::Split {
                body: w[1].into(),
                height: w[2].parse().unwrap(),
                operation: OperationId(w[3].parse().unwrap()),
                lower: w[4].into(),
                upper: w[5].into(),
            }),
            "fuse" => steps.push(Step::Fuse {
                a: w[1].into(),
                b: w[2].into(),
                operation: OperationId(w[3].parse().unwrap()),
                out: w[4].into(),
            }),
            "compose" => steps.push(Step::Compose(w[1].parse().unwrap(), w[2].parse().unwrap())),
            _ => sections.last_mut().unwrap().1.push(line.to_string()),
        }
    }
    let bodies = sections
        .into_iter()
        .map(|(body, rows)| {
            let text = format!("case {body} {tolerance}\n{}", rows.join("\n"));
            let spec = identity_protocol::parse(&text);
            let mut solid = identity_protocol::build(&spec);
            for t in &spec.transforms {
                solid = solid
                    .transform_with(OperationId::UNSPECIFIED, *t)
                    .unwrap()
                    .0;
            }
            (body, solid)
        })
        .collect();
    Scenario {
        name,
        bodies,
        steps,
    }
}

pub fn scenarios(text: &str) -> Vec<Scenario> {
    text.split("\nend")
        .filter(|b| !b.trim().is_empty())
        .map(|b| parse(b.trim()))
        .collect()
}

/// (section, key, text) rows, checking every history and output on the way.
pub fn evaluate(s: &Scenario) -> Vec<(String, String, String)> {
    let mut bodies: BTreeMap<String, Solid> = s.bodies.iter().cloned().collect();
    let mut histories: Vec<Option<History>> = Vec::new();
    let mut out = Vec::new();
    let mut produced = Vec::new();
    let set = |b: &Solid| b.topology().entity_set(b.profile().tolerance());
    for (k, step) in s.steps.iter().enumerate() {
        let result = match step {
            Step::Split {
                body,
                height,
                operation,
                lower,
                upper,
            } => bodies[body]
                .split_at_height(*operation, *height)
                .map(|([l, u], h)| {
                    assert_eq!(
                        check(&[set(&bodies[body])], &[set(&l), set(&u)], &h),
                        vec![],
                        "{}",
                        s.name
                    );
                    bodies.insert(lower.clone(), l);
                    bodies.insert(upper.clone(), u);
                    produced.extend([lower.clone(), upper.clone()]);
                    h
                }),
            Step::Fuse {
                a,
                b,
                operation,
                out: o,
            } => bodies[a]
                .fuse_stacked(&bodies[b], *operation)
                .map(|(f, h)| {
                    assert_eq!(
                        check(&[set(&bodies[a]), set(&bodies[b])], &[set(&f)], &h),
                        vec![],
                        "{}",
                        s.name
                    );
                    bodies.insert(o.clone(), f);
                    produced.push(o.clone());
                    h
                }),
            Step::Compose(i, j) => {
                let (first, second) = (histories[*i].as_ref(), histories[*j].as_ref());
                let composed = first.unwrap().then(second.unwrap()).unwrap();
                out.extend(
                    composed
                        .relations
                        .iter()
                        .map(|r| ("compose".into(), format!("{i} {j}"), relation_text(r))),
                );
                histories.push(Some(composed));
                continue;
            }
        };
        match result {
            Ok(h) => {
                out.extend(
                    h.relations
                        .iter()
                        .map(|r| ("step".into(), k.to_string(), relation_text(r))),
                );
                histories.push(Some(h));
            }
            Err(_) => {
                out.push(("step".into(), k.to_string(), "error".into()));
                histories.push(None);
            }
        }
    }
    for name in produced {
        let b = &bodies[&name];
        assert!(b.topology().check(b.profile().tolerance()).is_empty());
        out.push((
            "body".into(),
            name.clone(),
            format!("id {}", b.topology().body_id()),
        ));
        out.extend(
            rows(b)
                .into_iter()
                .map(|r| ("body".into(), name.clone(), r)),
        );
    }
    out
}
