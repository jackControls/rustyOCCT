//! Attributes through extrusion, transform, split and fuse (M4) against the
//! independent enumeration (fixtures/attribute-*.txt|tsv from
//! tools/generate_attribute_fixtures.py): every outcome and every output
//! attribute for the full policy matrix, each history also passing the
//! independent attribute checker in debug builds.
#[path = "support/split_merge_protocol.rs"]
mod split_merge_protocol;
use rusty_occt::identity::{fnv1a128, EntityId};
use split_merge_protocol::{evaluate, scenarios};
use std::collections::BTreeMap;

fn expected() -> BTreeMap<String, Vec<String>> {
    let mut map: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in include_str!("../../fixtures/attribute-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (case, rest) = line.split_once('\t').unwrap();
        map.entry(case.to_string())
            .or_default()
            .push(rest.to_string());
    }
    map
}

#[test]
fn every_outcome_and_attribute_matches_the_independent_enumeration() {
    let expected = expected();
    let all = scenarios(include_str!("../../fixtures/attribute-cases.txt"));
    assert_eq!(all.len(), expected.len());
    let mut failures = Vec::new();
    for s in &all {
        // Entity rows are the split/fuse fixtures' business.
        let got: Vec<(String, String, String)> =
            evaluate(s).into_iter().filter(|r| r.0 != "body").collect();
        let want = &expected[&s.name];
        let matches = match want[0].strip_prefix("digest\t") {
            Some(digest) => {
                let (count, hash) = digest.split_once('\t').unwrap();
                let text: Vec<String> =
                    got.iter().map(|(a, b, c)| format!("{a} {b} {c}")).collect();
                count.parse::<usize>().unwrap() == text.len()
                    && EntityId(fnv1a128(text.join("\n").as_bytes())).to_string() == hash
            }
            None => {
                let text: Vec<String> = got
                    .iter()
                    .map(|(a, b, c)| format!("{a}\t{b}\t{c}"))
                    .collect();
                text == *want
            }
        };
        if !matches {
            failures.push(s.name.clone());
        }
    }
    assert!(failures.is_empty(), "{failures:?}");
}

mod laws {
    use super::split_merge_protocol::scenarios;
    use rusty_occt::attributes::{
        Attribute, AttributeKey, OnMerge, OnModify, OnSplit, OnTransform, Policy, PolicyTable,
    };
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Context, Error, RigidTransform, Solid, Vec3};

    fn prism() -> Solid {
        scenarios(include_str!("../../fixtures/split-merge-cases.txt"))
            .into_iter()
            .find(|s| s.name == "split_holes_3")
            .unwrap()
            .bodies[0]
            .1
            .clone()
    }

    fn policy(on_transform: OnTransform) -> Policy {
        Policy {
            on_modify: OnModify::Keep,
            on_transform,
            on_split: OnSplit::CopyToAll,
            on_merge: OnMerge::KeepIfAllEqual,
        }
    }

    #[test]
    fn a_key_without_a_policy_or_callback_is_an_operation_error() {
        let p = prism();
        let (id, _) = p.topology().ids().next().unwrap();
        let tagged = p
            .with_attribute(
                id,
                Attribute {
                    key: AttributeKey(7),
                    value: vec![1],
                },
            )
            .unwrap();
        let t = RigidTransform::translation(Vec3::new(1.0, 0.0, 0.0)).unwrap();
        // No table: the operation fails rather than guess.
        assert_eq!(
            tagged.transform_with(OperationId(1), t).unwrap_err(),
            Error::MissingAttributePolicy(7)
        );
        assert_eq!(
            tagged.split_at_height(OperationId(1), 1.0).unwrap_err(),
            Error::MissingAttributePolicy(7)
        );
        let mut table = PolicyTable::default();
        table
            .0
            .insert(AttributeKey(7), policy(OnTransform::Recompute));
        let context = Context::new(OperationId(1)).with_policies(&table);
        assert_eq!(
            tagged.transform_in(&context, t).unwrap_err(),
            Error::MissingRecompute(7)
        );
        // An attribute on an id the body does not have is refused.
        let (other, _) =
            Solid::cuboid_with(OperationId(99), 1.0, 1.0, 1.0, p.profile().tolerance()).unwrap();
        let foreign = other.topology().ids().next().unwrap().0;
        assert!(p
            .with_attribute(
                foreign,
                Attribute {
                    key: AttributeKey(7),
                    value: vec![1]
                }
            )
            .is_err());
    }

    #[test]
    fn outcomes_are_deterministic_and_nothing_survives_on_a_deleted_id() {
        let p = prism();
        let mut table = PolicyTable::default();
        table.0.insert(AttributeKey(3), policy(OnTransform::Keep));
        let mut tagged = p.clone();
        for (id, _) in p.topology().ids() {
            tagged = tagged
                .with_attribute(
                    id,
                    Attribute {
                        key: AttributeKey(3),
                        value: vec![5, 5],
                    },
                )
                .unwrap();
        }
        let context = Context::new(OperationId(4)).with_policies(&table);
        let ([lower, upper], split) = tagged.split_at_height_in(&context, 1.25).unwrap();
        let again = tagged.split_at_height_in(&context, 1.25).unwrap();
        assert_eq!(again.1, split);
        assert_eq!(again.0, [lower.clone(), upper.clone()]);
        let (fused, fuse) = lower.fuse_stacked_in(&context, &upper).unwrap();
        for id in fuse.relations.iter().filter_map(|r| match r {
            rusty_occt::history::Relation::Deleted { id } => Some(*id),
            _ => None,
        }) {
            assert!(!fused.attributes().contains_key(&id));
        }
        // Equal copies merge back: every entity of the fused body keeps the value.
        assert_eq!(fused.attributes().len(), fused.topology().ids().count());
    }
}
