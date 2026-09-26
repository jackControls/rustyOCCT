//! Enclosures (M5, Contract 5 of IDENTITY_AND_HISTORY.md) on every prism of
//! the identity corpus: every entity is enclosed within the resolution and
//! the body validates against its stored bounds; through a rigid transform,
//! a split at mid-height and the fuse of the pieces, no continued entity's
//! bound falls below its parents' (T5).
#[path = "support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use rusty_occt::history::{History, Relation};
use rusty_occt::identity::{EntityId, OperationId};
use rusty_occt::topology::{Provenance, Slot};
use rusty_occt::{Point3, RigidTransform, Solid, Vec3};

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

fn enclosed(s: &Solid) {
    let tol = s.profile().tolerance();
    let t = s.topology();
    let all = t
        .vertices()
        .iter()
        .map(|v| v.enclosure)
        .chain(t.fins().iter().map(|f| f.enclosure))
        .chain(t.faces().iter().map(|f| f.enclosure));
    for e in all {
        let e = e.expect("every entity is enclosed");
        assert_eq!(e.provenance, Provenance::Computed);
        assert!(e.bound > 0.0 && e.bound <= tol.linear(), "{e:?}");
    }
    assert!(t.check(tol).is_empty());
}

fn carried(h: &History, inputs: &[&Solid], outputs: &[&Solid]) -> usize {
    let of = |id: EntityId, bodies: &[&Solid]| bodies.iter().find_map(|s| bound(s, id));
    let mut n = 0;
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
                n += 1;
            }
        }
    }
    n
}

#[test]
fn bounds_fit_the_resolution_and_never_fall_through_operations() {
    let specs = identity_protocol::cases(include_str!("../../fixtures/identity-cases.txt"));
    // Every prism in release builds; every eighth in debug builds, where the
    // certified arithmetic is slow.
    let step = if cfg!(debug_assertions) { 8 } else { 1 };
    let mut pairs = 0;
    for spec in specs.iter().step_by(step) {
        let solid = identity_protocol::build(spec);
        enclosed(&solid);
        let turn =
            RigidTransform::rotation(Point3::new(1.0, -2.0, 0.5), Vec3::new(0.3, 1.0, -0.2), 0.7)
                .unwrap();
        let (moved, h) = solid.transform_with(OperationId(900), turn).unwrap();
        enclosed(&moved);
        pairs += carried(&h, &[&solid], &[&moved]);
        let (low, high) = (
            moved.start_offset().min(moved.end_offset()),
            moved.start_offset().max(moved.end_offset()),
        );
        let ([lower, upper], split) = moved
            .split_at_height(OperationId(901), 0.5 * (low + high))
            .unwrap();
        enclosed(&lower);
        enclosed(&upper);
        pairs += carried(&split, &[&moved], &[&lower, &upper]);
        let (fused, fuse) = lower.fuse_stacked(&upper, OperationId(902)).unwrap();
        enclosed(&fused);
        pairs += carried(&fuse, &[&lower, &upper], &[&fused]);
    }
    assert!(pairs > 100_000 / step, "{pairs}");
}
