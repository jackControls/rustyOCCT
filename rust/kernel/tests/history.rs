//! Complete histories of extrusions and rigid transforms (M2) against the
//! independent enumeration (fixtures/history-expected.tsv), each checked by
//! the history checker, composed, and resolved.
#[path = "support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use identity_protocol::{build_tracked, cases, relation_text};
use rusty_occt::history::{check, History, Resolution};
use rusty_occt::identity::{fnv1a128, EntityId, OperationId, OperationKind};
use rusty_occt::topology::Slot;
use rusty_occt::Solid;
use std::collections::BTreeMap;

fn expected() -> BTreeMap<String, Vec<String>> {
    let mut map: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in include_str!("../../fixtures/history-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (case, rest) = line.split_once('\t').unwrap();
        let (step, relation) = rest.split_once('\t').unwrap();
        map.entry(case.to_string())
            .or_default()
            .push(format!("{step} {relation}"));
    }
    map
}

fn text(step: &str, h: &History) -> Vec<String> {
    h.relations
        .iter()
        .map(|r| format!("{step} {}", relation_text(r)))
        .collect()
}

#[test]
fn every_history_matches_the_independent_enumeration_and_checks_clean() {
    let expected = expected();
    let specs = cases(include_str!("../../fixtures/identity-cases.txt"));
    assert_eq!(specs.len(), expected.len());
    let mut relations = 0;
    let mut failures = Vec::new();
    for spec in &specs {
        let (solid, construct) = build_tracked(spec);
        let tolerance = spec.tolerance;
        let first = solid.topology().entity_set(tolerance);
        assert_eq!(
            check(&[], std::slice::from_ref(&first), &construct),
            vec![],
            "{}",
            spec.name
        );
        let mut got = text("construct", &construct);
        let mut current: Solid = solid;
        let mut composed = construct.clone();
        for (k, transform) in spec.transforms.iter().enumerate() {
            let (next, h) = current.transform_with(OperationId(77), *transform).unwrap();
            let before = current.topology().entity_set(tolerance);
            let after = next.topology().entity_set(tolerance);
            assert_eq!(check(&[before], &[after], &h), vec![], "{}", spec.name);
            assert_eq!(h.kind, OperationKind::Transform);
            for (id, _) in current.topology().ids() {
                assert_eq!(h.resolve(id), Resolution::Same(id));
            }
            got.extend(text(&format!("transform{k}"), &h));
            composed = composed.then(&h).unwrap();
            current = next;
        }
        if !spec.transforms.is_empty() {
            let last = current.topology().entity_set(tolerance);
            assert_eq!(check(&[], &[last], &composed), vec![], "{}", spec.name);
            got.extend(text("composed", &composed));
        }
        relations += got.len();
        let want = &expected[&spec.name];
        let matches = match want[0].strip_prefix("digest ") {
            Some(digest) => {
                // Region relations are listed after the digest, which covers the rest.
                let regions: Vec<String> = current
                    .topology()
                    .ids()
                    .filter(|(_, slot)| matches!(slot, Slot::Region(_)))
                    .map(|(id, _)| id.to_string())
                    .collect();
                let (regions, rest): (Vec<_>, Vec<_>) = got
                    .iter()
                    .cloned()
                    .partition(|r| r.split(' ').any(|w| regions.iter().any(|id| id == w)));
                let (count, hash) = digest.split_once(' ').unwrap();
                count.parse::<usize>().unwrap() == rest.len()
                    && EntityId(fnv1a128(rest.join("\n").as_bytes())).to_string() == hash
                    && regions == want[1..]
            }
            None => got == *want,
        };
        if !matches {
            failures.push(format!(
                "{}: {:?}",
                spec.name,
                got.iter().take(3).collect::<Vec<_>>()
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(relations, 114_270);
}

#[test]
#[allow(deprecated)]
fn deprecated_wrappers_discard_the_same_history() {
    let spec = &cases(include_str!("../../fixtures/identity-cases.txt"))[0];
    let (solid, _) = build_tracked(spec);
    let legacy = Solid::extrude(
        solid.profile().clone(),
        solid.frame(),
        solid.start_offset(),
        solid.end_offset(),
    )
    .unwrap();
    assert!(legacy.topology().ids().eq(Solid::extrude_with(
        OperationId::UNSPECIFIED,
        solid.profile().clone(),
        solid.frame(),
        solid.start_offset(),
        solid.end_offset(),
    )
    .unwrap()
    .0
    .topology()
    .ids()));
    let transform =
        rusty_occt::RigidTransform::translation(rusty_occt::Vec3::new(1.0, 2.0, 3.0)).unwrap();
    assert_eq!(
        solid.transformed(transform).unwrap(),
        solid.transform_with(OperationId(5), transform).unwrap().0
    );
}

#[test]
fn a_construction_cannot_follow_its_own_output() {
    let spec = &cases(include_str!("../../fixtures/identity-cases.txt"))[0];
    let (solid, construct) = build_tracked(spec);
    let (_, moved) = solid
        .transform_with(OperationId(1), rusty_occt::RigidTransform::identity())
        .unwrap();
    // A construction has no input bodies: it cannot follow a transform.
    assert!(moved.then(&construct).is_err());
    assert!(construct.then(&moved).is_ok());
}

#[test]
fn operations_record_their_level_and_replay_at_it() {
    use rusty_occt::identity::AlgorithmLevel;
    use rusty_occt::{Error, RigidTransform, Vec3};
    let spec = &cases(include_str!("../../fixtures/identity-cases.txt"))[0];
    let (solid, construct) = build_tracked(spec);
    assert_eq!(construct.level, AlgorithmLevel::CURRENT);
    let replayed = Solid::extrude_at(
        construct.level,
        solid.operation(),
        solid.profile().clone(),
        solid.frame(),
        solid.start_offset(),
        solid.end_offset(),
    )
    .unwrap();
    assert_eq!(replayed, (solid.clone(), construct));
    let transform = RigidTransform::translation(Vec3::new(1.0, 2.0, 3.0)).unwrap();
    let (moved, h) = solid.transform_with(OperationId(4), transform).unwrap();
    assert_eq!(h.level, AlgorithmLevel::CURRENT);
    assert_eq!(
        solid
            .transform_at(h.level, OperationId(4), transform)
            .unwrap(),
        (moved, h)
    );
    // No level before the first, and none this build does not provide.
    for level in [
        AlgorithmLevel(0),
        AlgorithmLevel(2),
        AlgorithmLevel(u32::MAX),
    ] {
        assert_eq!(
            solid
                .transform_at(level, OperationId(4), transform)
                .unwrap_err(),
            Error::UnknownAlgorithmLevel(level.0)
        );
        assert_eq!(
            Solid::extrude_at(
                level,
                solid.operation(),
                solid.profile().clone(),
                solid.frame(),
                solid.start_offset(),
                solid.end_offset(),
            )
            .unwrap_err(),
            Error::UnknownAlgorithmLevel(level.0)
        );
    }
}
