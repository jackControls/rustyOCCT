//! Complete histories of extrusions and rigid transforms (M2 of
//! IDENTITY_AND_HISTORY.md). Every produced history must check clean, resolve
//! every input and compose associatively; one mutation per input must be
//! reported with its predicted issue on the mutated relation.
use crate::identity::{build_tracked, spec};
use libfuzzer_sys::arbitrary::Unstructured;
use rusty_occt::history::{
    check, EntityInfo, EntitySet, Geometry, History, HistoryIssueKind as K, Relation, Resolution,
};
use rusty_occt::identity::{
    AlgorithmLevel, EntityId, EntityKind, OperationId, OperationKind, Role,
};
use rusty_occt::topology::Curve3;
use rusty_occt::Solid;

fn kinds(issues: &[rusty_occt::history::HistoryIssue]) -> Vec<(K, EntityId)> {
    issues.iter().map(|i| (i.kind, i.id)).collect()
}

fn has(issues: &[(K, EntityId)], kind: K, id: EntityId) -> bool {
    issues.contains(&(kind, id))
}

/// Replace one line edge of `set` by its two halves, as a split would.
fn split_edge(set: &EntitySet, pick: usize) -> Option<(EntityId, EntitySet, [EntityId; 2])> {
    let lines: Vec<(&EntityId, &EntityInfo)> = set
        .entities
        .iter()
        .filter(|(_, e)| matches!(e.geometry, Geometry::Curve(Curve3::LineSegment { .. })))
        .collect();
    let (id, info) = lines.get(pick % lines.len().max(1))?;
    let Geometry::Curve(Curve3::LineSegment { start, end }) = info.geometry else {
        return None;
    };
    let mid = start + (end - start) * 0.5;
    let mut out = set.clone();
    out.entities.remove(id);
    let children = [EntityId([0xa5; 16]), EntityId([0x5a; 16])];
    for (k, (a, b)) in [(start, mid), (mid, end)].into_iter().enumerate() {
        out.entities.insert(
            children[k],
            EntityInfo {
                kind: EntityKind::Edge,
                ordinal: k as u32,
                geometry: Geometry::Curve(Curve3::LineSegment { start: a, end: b }),
                structure: Vec::new(),
            },
        );
    }
    Some((**id, out, children))
}

pub fn check_history(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = spec(&mut u) else {
        return;
    };
    let mutation = u.arbitrary::<u8>().unwrap_or(0) % 7;
    let pick = usize::from(u.arbitrary::<u16>().unwrap_or(0));
    // Either extrusion direction, circles included (seamless ring edges).
    let reverse = u.arbitrary::<bool>().unwrap_or(false);
    let Some((solid, construct)) = build_tracked(&s, s.labels.as_deref(), 1.0, reverse) else {
        return;
    };
    let tol = s.tolerance;
    let set = |x: &Solid| x.topology().entity_set(tol);
    // The construction and every transform check clean and compose.
    assert_eq!(check(&[], &[set(&solid)], &construct), vec![]);
    assert_eq!(construct.kind, OperationKind::Extrude);
    // H8: the construction replays identically at its recorded level.
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
    assert_eq!(replayed, (solid.clone(), construct.clone()));
    // The body's one region is generated from every profile boundary.
    let regions = construct
        .relations
        .iter()
        .filter(|r| {
            matches!(
                r,
                Relation::Generated {
                    role: Role::Region,
                    ..
                }
            )
        })
        .count();
    assert_eq!(regions, 1);
    let mut steps = vec![construct.clone()];
    let mut sets = vec![set(&solid)];
    let mut current = solid.clone();
    for (k, transform) in s.transforms.iter().enumerate() {
        let Ok((next, h)) = current.transform_with(OperationId(k as u64), *transform) else {
            break;
        };
        assert_eq!(check(&[set(&current)], &[set(&next)], &h), vec![]);
        assert_eq!(
            current.transform_at(h.level, OperationId(k as u64), *transform),
            Ok((next.clone(), h.clone()))
        );
        for (id, _) in current.topology().ids() {
            assert_eq!(h.resolve(id), Resolution::Same(id));
        }
        sets.push(set(&next));
        steps.push(h);
        current = next;
    }
    let mut composed = steps[0].clone();
    for h in &steps[1..] {
        composed = composed.then(h).unwrap();
    }
    assert_eq!(
        check(&[], &[sets.last().unwrap().clone()], &composed),
        vec![]
    );
    if steps.len() == 3 {
        let right = steps[1].then(&steps[2]).unwrap();
        assert_eq!(steps[0].then(&right).unwrap(), composed, "associativity");
    }
    let (last, last_set) = (steps.last().unwrap().clone(), sets.last().unwrap().clone());
    let inputs: Vec<EntitySet> = if steps.len() > 1 {
        vec![sets[sets.len() - 2].clone()]
    } else {
        vec![]
    };
    let n = last.relations.len();
    let target = |h: &History, i: usize| h.relations[i].targets()[0];
    let run = |h: &History| kinds(&check(&inputs, std::slice::from_ref(&last_set), h));
    match mutation {
        0 => {
            // Dropping a relation loses its target (and source, if any).
            let mut h = last.clone();
            let gone = h.relations.remove(pick % n);
            let issues = run(&h);
            assert!(
                has(&issues, K::MissingTarget, gone.targets()[0]),
                "{issues:?}"
            );
            for source in gone.sources() {
                assert!(has(&issues, K::MissingSource, source));
            }
        }
        1 => {
            let mut h = last.clone();
            let copy = h.relations[pick % n].clone();
            h.relations.insert(pick % n, copy.clone());
            let issues = run(&h);
            let id = copy.sources().first().copied().unwrap_or(copy.targets()[0]);
            assert!(has(&issues, K::RelationOrder, id), "{issues:?}");
            for source in copy.sources() {
                assert!(has(&issues, K::DuplicateSource, source));
            }
        }
        2 => {
            // Changing a kind: a generated target becomes a phantom source;
            // a rigid motion reported as Unchanged differs in geometry.
            let mut h = last.clone();
            let i = pick % n;
            let id = target(&h, i);
            h.relations[i] = match &h.relations[i] {
                Relation::Generated { to, .. } => Relation::Modified { from: *to, to: *to },
                _ => Relation::Unchanged { id },
            };
            let issues = run(&h);
            if steps.len() == 1 {
                assert!(has(&issues, K::DanglingId, id), "{issues:?}");
            } else if inputs[0].entities[&id].geometry != last_set.entities[&id].geometry {
                assert!(has(&issues, K::UnchangedGeometryDiffers, id), "{issues:?}");
            }
        }
        3 => {
            let mut h = last.clone();
            let i = pick % n;
            let old = target(&h, i);
            let phantom = EntityId([0xee; 16]);
            match &mut h.relations[i] {
                Relation::Generated { to, .. } | Relation::Modified { to, .. } => *to = phantom,
                _ => unreachable!("extrusions and transforms report no other kind"),
            }
            let issues = run(&h);
            assert!(has(&issues, K::DanglingId, phantom), "{issues:?}");
            assert!(has(&issues, K::MissingTarget, old));
        }
        4 | 5 => {
            // A split of a real edge: clean in canonical order, reported when
            // the children are swapped.
            let Some((parent, out, [a, b])) = split_edge(&last_set, pick) else {
                return;
            };
            let relations = |into: Vec<EntityId>| {
                last_set
                    .entities
                    .keys()
                    .filter(|id| **id != parent)
                    .map(|id| Relation::Unchanged { id: *id })
                    .chain([Relation::Split { from: parent, into }])
                    .collect()
            };
            let mut src = last_set.clone();
            src.body = EntityId([1; 16]);
            let mut dst = out;
            dst.body = EntityId([2; 16]);
            let history = |into| {
                History::new(
                    OperationId(9),
                    OperationKind::Extrude,
                    vec![src.body],
                    vec![dst.body],
                    relations(into),
                    vec![],
                )
            };
            let clean = check(&[src.clone()], &[dst.clone()], &history(vec![a, b]));
            assert_eq!(clean, vec![]);
            let swapped = kinds(&check(&[src.clone()], &[dst.clone()], &history(vec![b, a])));
            let mut want = vec![(K::OrdinalNotCanonical, a), (K::OrdinalNotCanonical, b)];
            want.sort();
            assert_eq!(swapped, want);
        }
        _ => {
            // Out of order: a construction cannot follow its own output.
            if steps.len() > 1 {
                assert_eq!(
                    steps[1].then(&steps[0]).unwrap_err().kind,
                    K::CompositionInvalid
                );
            }
        }
    }
}
