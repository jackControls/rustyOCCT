//! Scenario protocol written by rust/tools/generate_split_merge_fixtures.py
//! and generate_attribute_fixtures.py: named prisms in the identity case
//! protocol, attribute `policy` rows and `attribute` rows (an entity by its
//! structural locator), then `transform`, `split`, `fuse` and `compose`
//! steps. `evaluate` runs a scenario and returns its rows in the fixtures'
//! order and format.
#[path = "identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
pub use identity_protocol::{relation_text, rows};
use rusty_occt::attributes::{
    Attribute, AttributeKey, OnMerge, OnModify, OnSplit, OnTransform, OutcomeResult, Policy,
    PolicyTable,
};
use rusty_occt::history::{check, History};
use rusty_occt::identity::{EntityId, OperationId};
use rusty_occt::{Context, RigidTransform, Solid, Vec3};
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
    Transform {
        body: String,
        operation: OperationId,
        by: [f64; 3],
        out: String,
    },
}

pub struct Scenario {
    pub name: String,
    pub bodies: Vec<(String, Solid)>,
    pub steps: Vec<Step>,
    pub policies: PolicyTable,
}

fn policy(w: &[&str]) -> Policy {
    Policy {
        on_modify: if w[0] == "keep" {
            OnModify::Keep
        } else {
            OnModify::Drop
        },
        on_transform: match w[1] {
            "keep" => OnTransform::Keep,
            "drop" => OnTransform::Drop,
            _ => OnTransform::Recompute,
        },
        on_split: if w[2] == "copy" {
            OnSplit::CopyToAll
        } else {
            OnSplit::Drop
        },
        on_merge: if w[3] == "keep_if_equal" {
            OnMerge::KeepIfAllEqual
        } else {
            OnMerge::Drop
        },
    }
}

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
        .collect()
}

/// The fixtures' recompute callback: the old value followed by 0x2a.
pub fn recompute(_: AttributeKey, _: EntityId, value: &[u8]) -> Vec<u8> {
    let mut v = value.to_vec();
    v.push(0x2a);
    v
}

fn result_name(r: OutcomeResult) -> &'static str {
    match r {
        OutcomeResult::Kept => "kept",
        OutcomeResult::Copied => "copied",
        OutcomeResult::Recomputed => "recomputed",
        OutcomeResult::Dropped => "dropped",
        OutcomeResult::Conflict => "conflict",
    }
}

pub fn parse(block: &str) -> Scenario {
    let mut lines = block.lines();
    let head: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    let (name, tolerance) = (head[1].to_string(), head[2]);
    let mut sections: Vec<(String, Vec<String>)> = Vec::new();
    let mut steps = Vec::new();
    let mut policies = PolicyTable::default();
    let mut attributes: Vec<(String, u64, Vec<u8>, String)> = Vec::new();
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
            "transform" if w.len() == 7 => steps.push(Step::Transform {
                body: w[1].into(),
                operation: OperationId(w[2].parse().unwrap()),
                by: [w[3], w[4], w[5]].map(|x| x.parse().unwrap()),
                out: w[6].into(),
            }),
            "policy" => {
                policies
                    .0
                    .insert(AttributeKey(w[1].parse().unwrap()), policy(&w[2..]));
            }
            "attribute" => attributes.push((
                w[1].into(),
                w[2].parse().unwrap(),
                unhex(w[3]),
                w[4..].join(" "),
            )),
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
        .collect::<Vec<(String, Solid)>>();
    // Attributes land on the entity each locator names.
    let bodies = bodies
        .into_iter()
        .map(|(body, mut solid)| {
            let located: BTreeMap<String, EntityId> = rows(&solid)
                .into_iter()
                .map(|r| {
                    let w: Vec<&str> = r.splitn(6, ' ').collect();
                    (w[5].to_string(), w[0].parse().unwrap())
                })
                .collect();
            for (b, key, value, loc) in &attributes {
                if *b == body {
                    let attribute = Attribute {
                        key: AttributeKey(*key),
                        value: value.clone(),
                    };
                    solid = solid.with_attribute(located[loc], attribute).unwrap();
                }
            }
            (body, solid)
        })
        .collect();
    Scenario {
        name,
        policies,
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
    let context = |op: OperationId| {
        Context::new(op)
            .with_policies(&s.policies)
            .with_recompute(&recompute)
    };
    for (k, step) in s.steps.iter().enumerate() {
        let result = match step {
            Step::Transform {
                body,
                operation,
                by,
                out: o,
            } => bodies[body]
                .transform_in(
                    &context(*operation),
                    RigidTransform::translation(Vec3::new(by[0], by[1], by[2])).unwrap(),
                )
                .map(|(m, h)| {
                    assert_eq!(
                        check(&[set(&bodies[body])], &[set(&m)], &h),
                        vec![],
                        "{}",
                        s.name
                    );
                    bodies.insert(o.clone(), m);
                    produced.push(o.clone());
                    h
                }),
            Step::Split {
                body,
                height,
                operation,
                lower,
                upper,
            } => bodies[body]
                .split_at_height_in(&context(*operation), *height)
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
                .fuse_stacked_in(&context(*operation), &bodies[b])
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
                out.extend(h.attributes.iter().map(|o| {
                    let to: Vec<String> = o.to.iter().map(|t| t.to_string()).collect();
                    let to = if to.is_empty() {
                        "-".into()
                    } else {
                        to.join(",")
                    };
                    let text = format!(
                        "outcome {} {} {to} {}",
                        o.key.0,
                        o.from,
                        result_name(o.result)
                    );
                    ("step".into(), k.to_string(), text)
                }));
                histories.push(Some(h));
            }
            Err(_) => {
                out.push(("step".into(), k.to_string(), "error".into()));
                histories.push(None);
            }
        }
    }
    for name in &produced {
        let mut lines: Vec<String> = bodies[name]
            .attributes()
            .iter()
            .flat_map(|(id, list)| {
                list.iter().map(move |a| {
                    let hex: String = a.value.iter().map(|b| format!("{b:02x}")).collect();
                    format!("{id} {} {hex}", a.key.0)
                })
            })
            .collect();
        lines.sort();
        out.extend(lines.into_iter().map(|l| ("attr".into(), name.clone(), l)));
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
