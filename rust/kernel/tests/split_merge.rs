//! Height split and stacked fuse (M3) against the independent enumeration
//! (fixtures/split-merge-*.txt|tsv from tools/generate_split_merge_fixtures.py):
//! every relation, composition, output id and entity row, each history
//! checked and each output valid.
#[path = "support/split_merge_protocol.rs"]
mod split_merge_protocol;
use rusty_occt::identity::{fnv1a128, EntityId};
use split_merge_protocol::{evaluate, scenarios};
use std::collections::BTreeMap;

fn expected() -> BTreeMap<String, Vec<String>> {
    let mut map: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in include_str!("../../fixtures/split-merge-expected.tsv")
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
fn every_scenario_matches_the_independent_enumeration() {
    let expected = expected();
    let all = scenarios(include_str!("../../fixtures/split-merge-cases.txt"));
    assert_eq!(all.len(), expected.len());
    let mut failures = Vec::new();
    let mut rows = 0;
    for s in &all {
        let got = evaluate(s);
        rows += got.len();
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
    assert!(rows > 0);
}

mod laws {
    use super::split_merge_protocol::scenarios;
    use rusty_occt::history::{check, Resolution};
    use rusty_occt::identity::{AlgorithmLevel, OperationId};
    use rusty_occt::{Error, RigidTransform, Solid, Vec3};

    fn bases() -> Vec<Solid> {
        scenarios(include_str!("../../fixtures/split-merge-cases.txt"))
            .into_iter()
            .filter(|s| s.name.starts_with("split_") && !s.name.contains("rejected"))
            .map(|s| s.bodies[0].1.clone())
            .collect()
    }

    fn set(s: &Solid) -> rusty_occt::history::EntitySet {
        s.topology().entity_set(s.profile().tolerance())
    }

    #[test]
    fn split_then_fuse_restores_the_prism_and_resolves_every_input() {
        for parent in bases() {
            let (low, high) = (
                parent.start_offset().min(parent.end_offset()),
                parent.start_offset().max(parent.end_offset()),
            );
            let h = low + (high - low) * 0.625;
            let ([lower, upper], split) = parent.split_at_height(OperationId(5), h).unwrap();
            // Volumes add up to the whole within the rounding of area x height.
            let v = parent.mass_properties().volume;
            let sum = lower.mass_properties().volume + upper.mass_properties().volume;
            assert!((sum - v).abs() <= 8.0 * f64::EPSILON * v, "{sum} {v}");
            // Every input id resolves; an id that was no input is unknown.
            for (id, _) in parent.topology().ids() {
                assert_ne!(split.resolve(id), Resolution::Unknown);
            }
            assert_eq!(
                split.resolve(lower.topology().body_id()),
                Resolution::Unknown
            );
            let (fused, fuse) = lower.fuse_stacked(&upper, OperationId(6)).unwrap();
            // Parents are in axial order: the call order changes no id.
            let (again, other) = upper.fuse_stacked(&lower, OperationId(6)).unwrap();
            assert_eq!(again, fused);
            assert_eq!(other.relations, fuse.relations);
            assert_eq!(fused.mass_properties(), parent.mass_properties());
            assert_eq!(fused.bounds(), parent.bounds());
            let both = split.then(&fuse).unwrap();
            assert_eq!(check(&[set(&parent)], &[set(&fused)], &both), vec![]);
            for (id, _) in parent.topology().ids() {
                match both.resolve(id) {
                    Resolution::Same(to) => assert!(fused.topology().slot_of(to).is_some()),
                    other => panic!("{other:?}"),
                }
            }
        }
    }

    #[test]
    fn pieces_keep_their_ids_through_rigid_motion_and_replay_at_their_level() {
        let parent = &bases()[3];
        let h = 0.5 * (parent.start_offset() + parent.end_offset());
        let ([lower, _], split) = parent.split_at_height(OperationId(5), h).unwrap();
        assert_eq!(split.level(), Some(AlgorithmLevel::CURRENT));
        assert_eq!(
            parent
                .split_at_height_at(split.level().unwrap(), OperationId(5), h)
                .unwrap()
                .0[0],
            lower
        );
        let t = RigidTransform::rotation(rusty_occt::Point3::ORIGIN, Vec3::new(0.0, 0.6, 0.8), 1.3)
            .unwrap();
        let (moved, h1) = lower.transform_with(OperationId(9), t).unwrap();
        assert_eq!(check(&[set(&lower)], &[set(&moved)], &h1), vec![]);
        assert!(lower.topology().ids().eq(moved.topology().ids()));
        for level in [AlgorithmLevel(0), AlgorithmLevel(2)] {
            assert_eq!(
                parent
                    .split_at_height_at(level, OperationId(5), h)
                    .unwrap_err(),
                Error::UnknownAlgorithmLevel(level.0)
            );
            let (a, b) = (&lower, &moved);
            assert_eq!(
                a.fuse_stacked_at(level, b, OperationId(6)).unwrap_err(),
                Error::UnknownAlgorithmLevel(level.0)
            );
        }
    }

    #[test]
    fn a_split_height_that_is_not_a_number_is_outside_the_prism() {
        let parent = &bases()[0];
        for h in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(parent.split_at_height(OperationId(5), h).is_err());
        }
    }
}

mod support {
    use super::split_merge_protocol::scenarios;
    use rusty_occt::history::{check, EntitySet, Geometry, HistoryIssueKind as K};
    use rusty_occt::identity::{OperationId, Role};
    use rusty_occt::topology::Surface;
    use rusty_occt::{Frame3, Solid};

    fn base(name: &str) -> Solid {
        scenarios(include_str!("../../fixtures/split-merge-cases.txt"))
            .into_iter()
            .find(|s| s.name == name)
            .unwrap()
            .bodies[0]
            .1
            .clone()
    }

    /// The split history's issues after `change` edits the upper piece's first
    /// wall surface.
    fn perturbed(name: &str, change: impl Fn(&Surface, f64) -> Surface) -> Vec<K> {
        let parent = base(name);
        let tol = parent.profile().tolerance();
        let ([lower, upper], h) = parent.split_at_height(OperationId(5), 1.25).unwrap();
        let set = |s: &Solid| s.topology().entity_set(tol);
        let mut up: EntitySet = set(&upper);
        let (id, info) = up
            .entities
            .iter_mut()
            .find(|(id, _)| upper.topology().derivation(**id).unwrap().role == Role::Wall)
            .unwrap();
        let Geometry::Surface { surface, .. } = &mut info.geometry else {
            unreachable!();
        };
        *surface = change(surface, tol.linear());
        let id = *id;
        check(&[set(&parent)], &[set(&lower), up], &h)
            .into_iter()
            .filter(|i| i.id == id)
            .map(|i| i.kind)
            .collect()
    }

    fn moved(frame: &Frame3, along: f64, across: f64, tol: rusty_occt::Tolerance) -> Frame3 {
        Frame3::new(
            frame.origin() + frame.normal() * along + frame.x() * across,
            frame.normal(),
            frame.x(),
            tol,
        )
        .unwrap()
    }

    #[test]
    fn a_plane_piece_must_lie_on_the_whole_plane_and_face_its_way() {
        let t = rusty_occt::Tolerance::default();
        // Sliding within the plane keeps the support; leaving it does not.
        assert_eq!(
            perturbed("split_polygon_5", |s, _| match s {
                Surface::Plane(f) => Surface::Plane(moved(f, 0.0, 3.0, t)),
                _ => unreachable!(),
            }),
            vec![]
        );
        assert_eq!(
            perturbed("split_polygon_5", |s, tol| match s {
                Surface::Plane(f) => Surface::Plane(moved(f, 2.0 * tol, 0.0, t)),
                _ => unreachable!(),
            }),
            vec![K::SplitSupportDiffers]
        );
        // The boundary is what counts: a tilted plane through the same
        // origin still passes if the piece's edges stay on the whole's plane,
        // but a flipped normal never does.
        assert_eq!(
            perturbed("split_polygon_5", |s, _| match s {
                Surface::Plane(f) => {
                    Surface::Plane(Frame3::new(f.origin(), -f.normal(), f.x(), t).unwrap())
                }
                _ => unreachable!(),
            }),
            vec![K::SplitSupportDiffers]
        );
    }

    #[test]
    fn a_cylinder_piece_may_differ_by_axis_distance_plus_radius_up_to_the_resolution() {
        let t = rusty_occt::Tolerance::default();
        let shifted = |across: f64, grow: f64| {
            move |s: &Surface, tol: f64| match s {
                Surface::Cylinder { frame, radius } => Surface::Cylinder {
                    frame: moved(frame, 0.0, across * tol, t),
                    radius: radius + grow * tol,
                },
                _ => unreachable!(),
            }
        };
        assert_eq!(perturbed("split_circle", shifted(0.4, 0.4)), vec![]);
        assert_eq!(
            perturbed("split_circle", shifted(0.6, 0.6)),
            vec![K::SplitSupportDiffers]
        );
        assert_eq!(
            perturbed("split_circle", shifted(0.0, 1.5)),
            vec![K::SplitSupportDiffers]
        );
    }
}
