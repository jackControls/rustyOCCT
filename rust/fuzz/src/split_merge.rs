//! Height split and stacked fuse (M3 of IDENTITY_AND_HISTORY.md) on the
//! structure-aware prisms of the `identity` target. A random split height,
//! then the fuse of the pieces, and a fuse of two independently extruded
//! stacked prisms. Every history checks clean; the split followed by the fuse
//! resolves every input to `Same`; volumes add up within the rounding of
//! area x height; every input id resolves, no output id does. One mutation
//! per input must be reported with its predicted issue. Enclosures (M5, T5)
//! never fall: every continued entity's bound is at least its parents', and
//! every bound lies within the resolution.
use crate::identity::{build, spec};
use libfuzzer_sys::arbitrary::Unstructured;
use rusty_occt::history::{check, EntitySet, History, HistoryIssueKind as K, Relation, Resolution};
use rusty_occt::identity::{EntityId, OperationId};
use rusty_occt::topology::Slot;
use rusty_occt::Solid;

fn set(s: &Solid) -> EntitySet {
    s.topology().entity_set(s.resolution())
}

fn has(issues: &[rusty_occt::history::HistoryIssue], kind: K, id: EntityId) -> bool {
    issues.iter().any(|i| i.kind == kind && i.id == id)
}

fn resolves_exhaustively(h: &History, inputs: &[&Solid], outputs: &[&Solid]) {
    for s in inputs {
        for (id, _) in s.topology().ids() {
            assert_ne!(h.resolve(id), Resolution::Unknown);
        }
    }
    for s in outputs {
        for (id, _) in s.topology().ids() {
            if inputs.iter().all(|i| i.topology().slot_of(id).is_none()) {
                assert_eq!(h.resolve(id), Resolution::Unknown);
            }
        }
    }
}

/// An entity's enclosure bound: its own, or an edge's largest over its fins.
fn bound(s: &Solid, id: EntityId) -> Option<f64> {
    let t = s.topology();
    match t.slot_of(id)? {
        Slot::Vertex(v) => t.vertices()[v.index()].enclosure.map(|e| e.bound),
        Slot::Face(f) => t.faces()[f.index()].enclosure.map(|e| e.bound),
        Slot::Edge(e) => t.edges()[e.index()]
            .fins
            .iter()
            .filter_map(|k| t.fins()[k.index()].enclosure.map(|e| e.bound))
            .reduce(f64::max),
        Slot::Region(_) => None,
    }
}

/// T5: no continued entity's bound falls below a parent's; all bounds exist
/// and fit the resolution.
fn enclosures_carry(h: &History, inputs: &[&Solid], outputs: &[&Solid]) {
    let of = |id: EntityId, bodies: &[&Solid]| bodies.iter().find_map(|s| bound(s, id));
    for r in &h.relations {
        let pairs: Vec<(EntityId, EntityId)> = match r {
            Relation::Unchanged { id } => vec![(*id, *id)],
            Relation::Modified { from, to } => vec![(*from, *to)],
            Relation::Split { from, into } => into.iter().map(|t| (*from, *t)).collect(),
            Relation::Merged { from, into } => from.iter().map(|f| (*f, *into)).collect(),
            _ => Vec::new(),
        };
        for (from, to) in pairs {
            if let (Some(a), Some(b)) = (of(from, inputs), of(to, outputs)) {
                assert!(b >= a, "{from:?} {a} -> {to:?} {b}");
            }
        }
    }
    for s in outputs {
        let tol = s.resolution().linear();
        let t = s.topology();
        let all = t
            .vertices()
            .iter()
            .map(|v| v.enclosure)
            .chain(t.fins().iter().map(|f| f.enclosure))
            .chain(t.faces().iter().map(|f| f.enclosure));
        for e in all {
            let e = e.expect("every entity is enclosed");
            assert!(e.bound >= 0.0 && e.bound <= tol, "{e:?}");
        }
    }
}

pub fn check_split_merge(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = spec(&mut u) else {
        return;
    };
    let mutation = u.arbitrary::<u8>().unwrap_or(0) % 6;
    let t = f64::from(u.arbitrary::<u16>().unwrap_or(1)) / 65536.0;
    let pick = usize::from(u.arbitrary::<u16>().unwrap_or(0));
    let Some(mut parent) = build(&s, s.labels.as_deref(), 1.0, false) else {
        return;
    };
    for transform in &s.transforms {
        if let Ok((moved, h)) = parent.transform_with(OperationId::UNSPECIFIED, *transform) {
            enclosures_carry(&h, &[&parent], &[&moved]);
            parent = moved;
        }
    }
    let (low, high) = (
        parent.start_offset().min(parent.end_offset()),
        parent.start_offset().max(parent.end_offset()),
    );
    let h = low + (high - low) * t;
    let tol = s.tolerance.linear();
    let result = parent.split_at_height(OperationId(71), h);
    // A split fails exactly when the height is outside or a piece is thinner
    // than the tolerance.
    let thin = !(low < h && h < high) || !(h - low > tol) || !(high - h > tol);
    let Ok(([lower, upper], split)) = result else {
        assert!(thin, "split at {h} of [{low}, {high}] failed");
        return;
    };
    assert!(!thin);
    assert_eq!(
        check(&[set(&parent)], &[set(&lower), set(&upper)], &split),
        vec![]
    );
    resolves_exhaustively(&split, &[&parent], &[&lower, &upper]);
    enclosures_carry(&split, &[&parent], &[&lower, &upper]);
    let v = parent.mass_properties().volume;
    let sum = lower.mass_properties().volume + upper.mass_properties().volume;
    assert!((sum - v).abs() <= 8.0 * f64::EPSILON * v, "{sum} {v}");
    for piece in [&lower, &upper] {
        assert!(piece.topology().check(s.tolerance).is_empty());
    }
    // The pieces fuse back to the prism; every input resolves to one entity.
    let (fused, fuse) = lower.fuse_stacked(&upper, OperationId(72)).unwrap();
    assert_eq!(
        check(&[set(&lower), set(&upper)], &[set(&fused)], &fuse),
        vec![]
    );
    resolves_exhaustively(&fuse, &[&lower, &upper], &[&fused]);
    enclosures_carry(&fuse, &[&lower, &upper], &[&fused]);
    assert_eq!(fused.mass_properties(), parent.mass_properties());
    let both = split.then(&fuse).unwrap();
    assert_eq!(check(&[set(&parent)], &[set(&fused)], &both), vec![]);
    for (id, _) in parent.topology().ids() {
        assert!(matches!(both.resolve(id), Resolution::Same(_)));
    }
    // Fusing a piece with itself or with the other piece twice is refused.
    assert!(lower.fuse_stacked(&lower, OperationId(73)).is_err());
    assert!(fused.fuse_stacked(&upper, OperationId(73)).is_err());
    // One mutation of a produced history must be reported.
    let run_split = |h: &History| check(&[set(&parent)], &[set(&lower), set(&upper)], h);
    let run_fuse = |h: &History| check(&[set(&lower), set(&upper)], &[set(&fused)], h);
    let splits: Vec<usize> = (0..split.relations.len())
        .filter(|k| matches!(split.relations[*k], Relation::Split { .. }))
        .collect();
    let merges: Vec<usize> = (0..fuse.relations.len())
        .filter(|k| matches!(fuse.relations[*k], Relation::Merged { .. }))
        .collect();
    match mutation {
        0 => {
            // Swapped children: both ordinals are no longer canonical.
            let mut m = split.clone();
            let Relation::Split { into, .. } = &mut m.relations[splits[pick % splits.len()]] else {
                unreachable!();
            };
            into.swap(0, 1);
            let children = into.clone();
            let issues = run_split(&m);
            for child in &children {
                assert!(has(&issues, K::OrdinalNotCanonical, *child), "{issues:?}");
            }
        }
        1 => {
            // A deleted shared cap entity left out: its source is missing.
            let mut m = fuse.clone();
            let deleted: Vec<usize> = (0..m.relations.len())
                .filter(|k| matches!(m.relations[*k], Relation::Deleted { .. }))
                .collect();
            let Relation::Deleted { id } = m.relations.remove(deleted[pick % deleted.len()]) else {
                unreachable!();
            };
            assert!(has(&run_fuse(&m), K::MissingSource, id));
        }
        2 => {
            // A deleted entity reported unchanged: it is not in the output.
            let mut m = fuse.clone();
            let k = (0..m.relations.len())
                .filter(|k| matches!(m.relations[*k], Relation::Deleted { .. }))
                .nth(pick % 3)
                .unwrap();
            let Relation::Deleted { id } = m.relations[k] else {
                unreachable!();
            };
            m.relations[k] = Relation::Unchanged { id };
            assert!(has(&run_fuse(&m), K::DanglingId, id));
        }
        3 => {
            // A split parent given another parent's children.
            if splits.len() < 2 {
                return;
            }
            let (a, b) = (
                splits[pick % splits.len()],
                splits[(pick + 1) % splits.len()],
            );
            let mut m = split.clone();
            let Relation::Split { into: other, .. } = split.relations[b].clone() else {
                unreachable!();
            };
            let Relation::Split { from, into } = &mut m.relations[a] else {
                unreachable!();
            };
            let from = *from;
            *into = other.clone();
            let issues = run_split(&m);
            // Different support, or a different kind of entity altogether.
            assert!(
                other
                    .iter()
                    .any(|c| has(&issues, K::SplitSupportDiffers, *c)
                        || has(&issues, K::DimensionMismatch, *c)),
                "{from} {issues:?}"
            );
        }
        4 => {
            // A merge that names a parent twice loses the other parent.
            let mut m = fuse.clone();
            let Relation::Merged { from, .. } = &mut m.relations[merges[pick % merges.len()]]
            else {
                unreachable!();
            };
            let (kept, lost) = (from[0], from[1]);
            from[1] = kept;
            let issues = run_fuse(&m);
            assert!(has(&issues, K::MissingSource, lost), "{issues:?}");
            assert!(has(&issues, K::DuplicateSource, kept));
        }
        _ => {
            // Out of order: the fuse cannot precede the split.
            assert_eq!(fuse.then(&split).unwrap_err().kind, K::CompositionInvalid);
        }
    }
    // An independent construction stacked on the prism's end fuses with it.
    if s.operation == OperationId(74) {
        return;
    }
    let (start, end) = (parent.start_offset(), parent.end_offset());
    let Ok((above, _)) = Solid::extrude_with(
        OperationId(74),
        parent.profile().expect("a prism").clone(),
        parent.frame(),
        end,
        end + (end - start),
    ) else {
        return;
    };
    let (stacked, h2) = parent.fuse_stacked(&above, OperationId(75)).unwrap();
    assert_eq!(
        check(&[set(&parent), set(&above)], &[set(&stacked)], &h2),
        vec![]
    );
    resolves_exhaustively(&h2, &[&parent, &above], &[&stacked]);
    let (va, vb) = (v, above.mass_properties().volume);
    let vs = stacked.mass_properties().volume;
    assert!(
        (va + vb - vs).abs() <= 8.0 * f64::EPSILON * vs,
        "{va} {vb} {vs}"
    );
    // Parents are in axial order: the call order changes no id.
    let (swapped, h3) = above.fuse_stacked(&parent, OperationId(75)).unwrap();
    assert_eq!(swapped, stacked);
    assert_eq!(h3.relations, h2.relations);
    assert_eq!(
        h3.input_bodies,
        vec![h2.input_bodies[1], h2.input_bodies[0]]
    );
}
