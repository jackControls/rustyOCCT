//! Value ids on extrusions and rigid transforms against the independent
//! reference (fixtures/identity-*.txt|tsv from tools/generate_identity_fixtures.py).
#[path = "support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use identity_protocol::{build, cases, rows};
use rusty_occt::identity::{fnv1a128, EntityId, InputLabel, OperationId};
use rusty_occt::{Boundary, BoundaryLabels, Error, Frame3, Point2, Profile, Solid, Tolerance};
use std::collections::BTreeMap;

fn expected() -> BTreeMap<String, Vec<String>> {
    let mut map: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in include_str!("../../fixtures/identity-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (case, row) = line.split_once('\t').unwrap();
        map.entry(case.to_string())
            .or_default()
            .push(row.to_string());
    }
    map
}

#[test]
fn every_case_matches_the_independent_ids_and_transforms_keep_them() {
    let expected = expected();
    let specs = cases(include_str!("../../fixtures/identity-cases.txt"));
    assert_eq!(specs.len(), expected.len());
    let mut failures = Vec::new();
    let mut entities = 0;
    for spec in &specs {
        let solid = build(spec);
        let got = rows(&solid);
        entities += got.len();
        let want = &expected[&spec.name];
        let matches = match want[0].strip_prefix("digest ") {
            Some(digest) => {
                // Region rows are listed after the digest, which covers the rest.
                let (regions, rest): (Vec<_>, Vec<_>) = got
                    .iter()
                    .cloned()
                    .partition(|r| r.contains(" region region "));
                let (count, hash) = digest.split_once(' ').unwrap();
                count.parse::<usize>().unwrap() == rest.len()
                    && EntityId(fnv1a128(rest.join("\n").as_bytes())).to_string() == hash
                    && regions == want[1..]
            }
            None => got == *want,
        };
        if !matches {
            failures.push(format!(
                "{}:\n  got  {:?}",
                spec.name,
                got.iter().take(4).collect::<Vec<_>>()
            ));
        }
        // Every rigid motion keeps every id, derivation and locator.
        let mut moved = solid.clone();
        for transform in &spec.transforms {
            moved = moved
                .transform_with(OperationId::UNSPECIFIED, *transform)
                .map(|(s, _)| s)
                .unwrap();
            assert_eq!(
                moved.topology().body_id(),
                solid.topology().body_id(),
                "{}",
                spec.name
            );
            assert_eq!(rows(&moved), got, "{}: transform", spec.name);
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(specs.len(), 546);
    assert_eq!(entities, 46494);
}

#[test]
fn labels_decide_ids_and_indices_stand_in_without_them() {
    let specs = cases(include_str!("../../fixtures/identity-cases.txt"));
    let by_name: BTreeMap<&str, Vec<String>> = specs
        .iter()
        .map(|s| (s.name.as_str(), rows(&build(s))))
        .collect();
    let ids = |name: &str| -> Vec<String> {
        let mut v: Vec<String> = by_name[name].iter().map(|r| r[..32].to_string()).collect();
        v.sort();
        v
    };
    // Moving labelled points keeps every id; permuting label values permutes
    // them (the same set of ids on different entities).
    assert_eq!(ids("labels_moved"), ids("labels_present"));
    assert_eq!(ids("labels_permuted"), ids("labels_present"));
    assert_ne!(by_name["labels_permuted"], by_name["labels_present"]);
    assert!(ids("labels_absent")
        .iter()
        .all(|id| !ids("labels_present").contains(id)));
    // Inserting a vertex adds its entities and keeps every other id.
    let before = ids("labels_present");
    let after = ids("labels_inserted_vertex");
    assert!(before.iter().all(|id| after.contains(id)));
    assert_eq!(after.len() - before.len(), 6);
    // Reversing the extrusion direction keeps every id.
    assert_eq!(ids("holes_3_reversed_direction"), ids("holes_3_labelled"));
    // Every box sign combination has the same ids (the cuboid's).
    for tag in ["ppn", "pnp", "pnn", "npp", "npn", "nnp", "nnn"] {
        assert_eq!(ids(&format!("box_{tag}")), ids("box_ppp"));
    }
}

#[test]
fn invalid_labels_are_rejected() {
    let t = Tolerance::default();
    let square = || {
        Boundary::polygon(
            vec![
                Point2::new(0., 0.),
                Point2::new(4., 0.),
                Point2::new(4., 4.),
                Point2::new(0., 4.),
            ],
            t,
        )
        .unwrap()
    };
    let labels = |b: u64, s: [u64; 4], v: [u64; 4]| BoundaryLabels {
        boundary: InputLabel(b),
        segments: s.iter().map(|x| InputLabel(*x)).collect(),
        vertices: v.iter().map(|x| InputLabel(*x)).collect(),
    };
    assert!(square()
        .with_labels(labels(0, [1, 2, 3, 4], [5, 6, 7, 8]))
        .is_ok());
    let short = BoundaryLabels {
        boundary: InputLabel(0),
        segments: vec![InputLabel(1)],
        vertices: vec![InputLabel(2)],
    };
    assert!(matches!(
        square().with_labels(short),
        Err(Error::InvalidLabel(_))
    ));
    assert!(matches!(
        square().with_labels(labels(0, [1, 2, 3, 4], [5, 6, 7, 1])),
        Err(Error::InvalidLabel(_))
    ));
    // Labels must also be distinct across the boundaries of one profile.
    let outer = Boundary::polygon(
        vec![
            Point2::new(-9., -9.),
            Point2::new(9., -9.),
            Point2::new(9., 9.),
            Point2::new(-9., 9.),
        ],
        t,
    )
    .unwrap()
    .with_labels(labels(10, [11, 12, 13, 14], [15, 16, 17, 18]))
    .unwrap();
    let hole = square()
        .with_labels(labels(20, [21, 22, 23, 24], [25, 26, 27, 18]))
        .unwrap();
    assert!(matches!(
        Profile::new(outer.clone(), vec![hole], t),
        Err(Error::InvalidLabel(_))
    ));
    let hole = square()
        .with_labels(labels(20, [21, 22, 23, 24], [25, 26, 27, 28]))
        .unwrap();
    let profile = Profile::new(outer, vec![hole], t).unwrap();
    let a = Solid::extrude_with(OperationId(1), profile.clone(), Frame3::xy(), 0., 1.)
        .unwrap()
        .0;
    let b = Solid::extrude_with(OperationId(2), profile, Frame3::xy(), 0., 1.)
        .unwrap()
        .0;
    // A different operation id changes every id.
    assert!(a
        .topology()
        .ids()
        .all(|(id, _)| b.topology().slot_of(id).is_none()));
}
