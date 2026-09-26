//! Attributes through operations (M4 of IDENTITY_AND_HISTORY.md, contract
//! 4). Random policies and random values on the `identity` target's prisms,
//! through a transform, a split at a random height, value changes on one
//! piece, and the fuse of the pieces. Every step's outcomes pass the
//! independent `attributes::check`; a merge conflicts exactly when its
//! parents' values differ or one lacks the key; no attribute survives on a
//! deleted id; re-running a step gives the same bodies and history; a key
//! with no policy is refused.
use crate::identity::{build, spec};
use libfuzzer_sys::arbitrary::Unstructured;
use rusty_occt::attributes::{
    check, Attribute, AttributeKey, AttributeMap, OnMerge, OnModify, OnSplit, OnTransform,
    OutcomeResult, Policy, PolicyTable,
};
use rusty_occt::history::{History, Relation};
use rusty_occt::identity::{EntityId, OperationId};
use rusty_occt::{Context, Error, RigidTransform, Solid, Vec3};

const KEYS: u64 = 6;

fn append(_: AttributeKey, _: EntityId, value: &[u8]) -> Vec<u8> {
    let mut v = value.to_vec();
    v.push(0x2a);
    v
}

fn union(bodies: &[&Solid]) -> AttributeMap {
    bodies.iter().flat_map(|b| b.attributes().clone()).collect()
}

fn value(map: &AttributeMap, id: EntityId, key: AttributeKey) -> Option<&Vec<u8>> {
    map.get(&id)?
        .iter()
        .find(|a| a.key == key)
        .map(|a| &a.value)
}

/// The checker is clean, conflicts are exactly the unequal merges, and no
/// attribute lands on a deleted id.
fn verify(table: &PolicyTable, inputs: &[&Solid], outputs: &[&Solid], h: &History) {
    let (before, after) = (union(inputs), union(outputs));
    assert_eq!(check(&before, &after, h, table), vec![]);
    for r in &h.relations {
        match r {
            Relation::Deleted { id } => assert!(!after.contains_key(id)),
            Relation::Merged { from, .. } => {
                for key in (1..=KEYS).map(AttributeKey) {
                    if table.0[&key].on_merge != OnMerge::KeepIfAllEqual {
                        continue;
                    }
                    let values: Vec<_> = from.iter().map(|p| value(&before, *p, key)).collect();
                    let equal = values.iter().all(|v| v.is_some() && *v == values[0]);
                    for p in from {
                        let outcome = h.attributes.iter().find(|o| o.from == *p && o.key == key);
                        match outcome {
                            Some(o) => assert_eq!(
                                o.result == OutcomeResult::Conflict,
                                !equal,
                                "{values:?}"
                            ),
                            None => assert!(value(&before, *p, key).is_none()),
                        }
                    }
                }
            }
            _ => {}
        }
    }
}

pub fn check_attributes(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = spec(&mut u) else {
        return;
    };
    let Some(prism) = build(&s, s.labels.as_deref(), 1.0, false) else {
        return;
    };
    let byte = |u: &mut Unstructured| u.arbitrary::<u8>().unwrap_or(0);
    let mut table = PolicyTable::default();
    for key in 1..=KEYS {
        let b = byte(&mut u);
        table.0.insert(
            AttributeKey(key),
            Policy {
                on_modify: if b & 1 == 0 {
                    OnModify::Keep
                } else {
                    OnModify::Drop
                },
                on_transform: [OnTransform::Keep, OnTransform::Drop, OnTransform::Recompute]
                    [usize::from(b >> 1) % 3],
                on_split: if b & 8 == 0 {
                    OnSplit::CopyToAll
                } else {
                    OnSplit::Drop
                },
                on_merge: if b & 16 == 0 {
                    OnMerge::KeepIfAllEqual
                } else {
                    OnMerge::Drop
                },
            },
        );
    }
    // Random keys and small values, so equal and unequal merges both occur.
    let mut tagged = prism.clone();
    for (id, _) in prism.topology().ids() {
        let mask = byte(&mut u);
        for key in 1..=KEYS {
            if mask & (1 << (key - 1)) != 0 {
                let attribute = Attribute {
                    key: AttributeKey(key),
                    value: vec![byte(&mut u) % 3],
                };
                tagged = tagged.with_attribute(id, attribute).unwrap();
            }
        }
    }
    let context = |op| {
        Context::new(OperationId(op))
            .with_policies(&table)
            .with_recompute(&append)
    };
    // A key without a policy is refused, never guessed.
    if let Some((id, _)) = prism.topology().ids().next() {
        let stray = tagged
            .with_attribute(
                id,
                Attribute {
                    key: AttributeKey(KEYS + 1),
                    value: vec![0],
                },
            )
            .unwrap();
        let t = RigidTransform::translation(Vec3::new(1.0, 0.0, 0.0)).unwrap();
        assert_eq!(
            stray.transform_in(&context(1), t).unwrap_err(),
            Error::MissingAttributePolicy(KEYS + 1)
        );
    }
    let scale = s.tolerance.linear() / 1e-9;
    let t = RigidTransform::translation(Vec3::new(scale, -2.0 * scale, 0.5 * scale)).unwrap();
    let Ok((moved, h0)) = tagged.transform_in(&context(2), t) else {
        return;
    };
    verify(&table, &[&tagged], &[&moved], &h0);
    assert_eq!(
        tagged.transform_in(&context(2), t).unwrap(),
        (moved.clone(), h0)
    );
    let (low, high) = (
        moved.start_offset().min(moved.end_offset()),
        moved.start_offset().max(moved.end_offset()),
    );
    let fraction = f64::from(u.arbitrary::<u16>().unwrap_or(32768)) / 65536.0;
    let Ok(([lower, upper], h1)) =
        moved.split_at_height_in(&context(3), low + (high - low) * fraction)
    else {
        return;
    };
    verify(&table, &[&moved], &[&lower, &upper], &h1);
    // Change some values on the upper piece so merges can disagree.
    let mut changed = upper.clone();
    for (id, list) in upper.attributes() {
        for a in list {
            if byte(&mut u) % 4 == 0 {
                let attribute = Attribute {
                    key: a.key,
                    value: vec![9],
                };
                changed = changed.with_attribute(*id, attribute).unwrap();
            }
        }
    }
    let (fused, h2) = lower.fuse_stacked_in(&context(4), &changed).unwrap();
    verify(&table, &[&lower, &changed], &[&fused], &h2);
    assert_eq!(
        lower.fuse_stacked_in(&context(4), &changed).unwrap(),
        (fused, h2)
    );
}
