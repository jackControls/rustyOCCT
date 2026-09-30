//! S9e.2: a stack with an arc wall (S9a.2), an S9b.1 result of line prisms
//! given with arcs, and one solid of a first result of several, given to
//! another Boolean, against the independent reference
//! (`fixtures/boolean-given-*` from
//! `tools/generate_given_boolean_fixtures.py`): each case's first Boolean
//! and its second on the given solid (the first result's one solid, or the
//! one holding the case's pick point) and a third prism, as object or tool
//! (`swapped`). Each case runs once (on a few threads) for the checks that
//! read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::Error;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!("../../fixtures/boolean-given-cases.txt"))
}

/// `f` over the cases on at most six threads, in the cases' order.
fn each<T: Send>(f: impl Fn(&protocol::Case) -> T + Sync) -> Vec<(String, T)> {
    let cases = cases();
    let next = AtomicUsize::new(0);
    let workers = std::thread::available_parallelism()
        .map_or(2, |n| n.get())
        .clamp(1, 6);
    let mut out: Vec<(usize, String, T)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    let mut mine = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::SeqCst);
                        let Some(case) = cases.get(i) else { break };
                        mine.push((i, case.name.clone(), f(case)));
                    }
                    mine
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("a worker"))
            .collect()
    });
    out.sort_by_key(|x| x.0);
    out.into_iter().map(|(_, n, t)| (n, t)).collect()
}

/// Every case's run (the second Boolean's inputs, results and history),
/// once.
fn runs() -> &'static [(String, Result<protocol::Run, Error>)] {
    static RUNS: OnceLock<Vec<(String, Result<protocol::Run, Error>)>> = OnceLock::new();
    RUNS.get_or_init(|| each(protocol::run))
}

/// The reference's kind and `(solids, [volume, area])` per case.
type Expected = std::collections::BTreeMap<String, (String, Option<(usize, [f64; 2])>)>;

fn expected() -> Expected {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-given-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let entry = expect
            .entry(name.to_string())
            .or_insert((String::new(), None));
        match w[0] {
            "expect" => entry.0 = w[1].to_string(),
            "result" => {
                let v: Vec<f64> = w[2..4].iter().map(|x| x.parse().unwrap()).collect();
                entry.1 = Some((w[1].parse::<usize>().unwrap(), [v[0], v[1]]));
            }
            _ => {}
        }
    }
    expect
}

/// A result's volume and area enclosures within the reference's, each at
/// most `1e-9` wide (relative).
fn within(out: &[rusty_occt::Solid], count: usize, v: [f64; 2]) -> Result<(), String> {
    let mut sums = [0.0f64; 4];
    for s in out {
        let m = s
            .topology()
            .mass_enclosure()
            .ok_or("a result's mass undecided")?;
        for (k, x) in [
            m.volume[0],
            m.volume[1],
            m.surface_area[0],
            m.surface_area[1],
        ]
        .into_iter()
        .enumerate()
        {
            sums[k] += x;
        }
    }
    if out.len() != count {
        return Err(format!("{} solids for {count}", out.len()));
    }
    let near = |x: f64, lo: f64, hi: f64| {
        let slack = 1e-9 * x.abs().max(1.0);
        lo - slack <= x && x <= hi + slack
    };
    let narrow = |lo: f64, hi: f64| hi - lo <= 1e-9 * lo.abs().max(hi.abs()).max(1.0);
    if !narrow(sums[0], sums[1]) || !narrow(sums[2], sums[3]) {
        return Err(format!("wide enclosures {sums:?}"));
    }
    if !near(v[0], sums[0], sums[1]) || !near(v[1], sums[2], sums[3]) {
        return Err(format!("volume and area {v:?} against {sums:?}"));
    }
    Ok(())
}

#[test]
fn every_case_matches_the_reference() {
    let expect = expected();
    let mut failures = Vec::new();
    for (name, run) in runs() {
        let (kind, want) = &expect[name];
        match (kind.as_str(), run) {
            ("degenerate", Err(Error::Degenerate(_))) => continue,
            ("degenerate", other) => {
                failures.push(format!(
                    "{name}: {:?} not refused",
                    other.as_ref().map(|_| ())
                ));
                continue;
            }
            (_, Err(e)) => {
                failures.push(format!("{name}: {e}"));
                continue;
            }
            ("empty", Ok((_, _, out, _))) => {
                if !out.is_empty() {
                    failures.push(format!("{name}: {} solids, not empty", out.len()));
                }
                continue;
            }
            _ => {}
        }
        let Ok((_, _, out, _)) = run else {
            unreachable!("refused above")
        };
        let (count, v) = want.expect("a result");
        if let Err(e) = within(out, count, v) {
            failures.push(format!("{name}: {e}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures: {failures:#?}",
        failures.len()
    );
}

/// The second Boolean's history is complete over its inputs, the given
/// solid's entities among them: it chains onto the first Boolean's, every
/// entity of the given solid resolved by it.
#[test]
fn given_histories_are_complete() {
    for (name, run) in runs() {
        let Ok((a, b, out, h)) = run else {
            continue;
        };
        let ins = [
            a.topology().entity_set(a.resolution()),
            b.topology().entity_set(b.resolution()),
        ];
        let outs: Vec<_> = out
            .iter()
            .map(|s| s.topology().entity_set(s.resolution()))
            .collect();
        let issues = history::check(&ins, &outs, h);
        assert!(issues.is_empty(), "{name}: {issues:?}");
        for s in [a, b] {
            for (id, _) in s.topology().ids() {
                assert!(
                    !matches!(h.resolve(id), Resolution::Unknown)
                        || h.relations.iter().any(|r| r.sources().contains(&id)),
                    "{name}: an input entity the history does not name"
                );
            }
        }
    }
}

/// The given solid's ids are the second Boolean's inputs: the history of the
/// first Boolean names every entity of every solid of its result, the
/// second reads the given one's.
#[test]
fn the_given_solid_is_the_second_input() {
    for case in cases() {
        let Ok((a, b, first, h1)) = protocol::run_first(&case) else {
            panic!("{}: the first Boolean", case.name)
        };
        let then = case.then.as_ref().unwrap();
        assert_eq!(
            first.len() > 1,
            then.pick.is_some(),
            "{}: several solids exactly where a solid is picked",
            case.name
        );
        let ins = [
            a.topology().entity_set(a.resolution()),
            b.topology().entity_set(b.resolution()),
        ];
        let outs: Vec<_> = first
            .iter()
            .map(|s| s.topology().entity_set(s.resolution()))
            .collect();
        assert!(history::check(&ins, &outs, &h1).is_empty(), "{}", case.name);
        let picked = protocol::given(&case, first);
        let Some(Ok((x, y, _, _))) = runs()
            .iter()
            .find(|(n, _)| *n == case.name)
            .map(|(_, r)| r.as_ref())
        else {
            continue;
        };
        let given = if then.swapped { y } else { x };
        let ids = |s: &rusty_occt::Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
        assert_eq!(ids(given), ids(&picked), "{}", case.name);
    }
}

#[test]
fn results_are_deterministic_and_move_rigidly() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Location, Point3, RigidTransform, Vec3};
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ids = |s: &rusty_occt::Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    let first: std::collections::BTreeMap<&str, _> =
        runs().iter().map(|(n, r)| (n.as_str(), r)).collect();
    let failures: Vec<String> = each(|case| {
        let Ok((_, _, out, h)) = first[case.name.as_str()] else {
            return Vec::new();
        };
        let mut bad = Vec::new();
        let (_, _, again, h2) = protocol::run(case).unwrap();
        if *h != h2 {
            bad.push(format!("{}: histories differ", case.name));
        }
        for (x, y) in out.iter().zip(&again) {
            if ids(x) != ids(y) || x.topology().edges() != y.topology().edges() {
                bad.push(format!("{}: results differ", case.name));
            }
        }
        for s in out {
            let (moved, _) = s
                .transform_with(OperationId(900), motion)
                .unwrap_or_else(|e| panic!("{}: {e}", case.name));
            if ids(&moved) != ids(s) {
                bad.push(format!("{}: moved ids", case.name));
            }
            let (v0, v1) = (s.mass_properties().volume, moved.mass_properties().volume);
            if (v0 - v1).abs() > 1e-9 * v0.abs().max(1.0) {
                bad.push(format!("{}: moved volume {v0} {v1}", case.name));
            }
            for v in moved.topology().vertices() {
                if moved.classify(v.position).unwrap() != Location::Boundary {
                    bad.push(format!("{}: a moved vertex off the boundary", case.name));
                }
            }
        }
        bad
    })
    .into_iter()
    .flat_map(|(_, b)| b)
    .collect();
    assert!(failures.is_empty(), "{failures:#?}");
}

/// A given solid moved rigidly (a stack moved on its frame, an S9b.1
/// result with its inputs, one solid of several), given to the second
/// Boolean with the third prism moved alike, keeps the reference's volumes:
/// its model built again from its moved construction and matched to its
/// moved topology. The motion is a translation by binary64 steps, as
/// S9e.1's: a turned motion rounds the frames apart.
#[test]
fn moved_results_keep_their_volumes() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{RigidTransform, Vec3};
    let motion = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let expect = expected();
    for name in [
        "rollex_turned_cut",
        "rollex_flat_common",
        "boss_sliced_common",
        "boss_tool_cut",
        "turned_bored_cut",
        "slanted_bored_fuse",
        "severed_drilled_cut",
        "severed_tool_fuse",
        "split_bored_common",
    ] {
        let case = cases().into_iter().find(|c| c.name == name).unwrap();
        let (_, _, first, _) = protocol::run_first(&case).unwrap();
        let then = case.then.as_ref().unwrap();
        let given = protocol::given(&case, first);
        let r = given.transform_with(OperationId(801), motion).unwrap().0;
        let c = protocol::build(&then.third)
            .transform_with(OperationId(802), motion)
            .unwrap()
            .0;
        let (x, y) = if then.swapped { (&c, &r) } else { (&r, &c) };
        let out = match then.op.as_str() {
            "fuse" => x.fuse(then.operation, y),
            "cut" => x.cut(then.operation, y),
            _ => x.common(then.operation, y),
        }
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .0;
        let (count, v) = expect[name].1.unwrap();
        within(&out, count, v).unwrap_or_else(|e| panic!("{name}: {e}"));
    }
}

/// The given solid of several is the second input alone: the other solid's
/// entities are neither the second result's nor in its history.
#[test]
fn the_other_solid_is_not_an_input() {
    for case in cases()
        .iter()
        .filter(|c| c.then.as_ref().unwrap().pick.is_some())
    {
        let Some(Ok((_, _, out, h))) = runs()
            .iter()
            .find(|(n, _)| *n == case.name)
            .map(|(_, r)| r.as_ref())
        else {
            continue;
        };
        let (_, _, first, _) = protocol::run_first(case).unwrap();
        let picked = protocol::given(case, first.clone());
        let mine: std::collections::BTreeSet<_> =
            picked.topology().ids().map(|(id, _)| id).collect();
        let others: Vec<_> = first
            .iter()
            .flat_map(|s| s.topology().ids().map(|(id, _)| id))
            .filter(|id| !mine.contains(id))
            .collect();
        assert!(!others.is_empty(), "{}", case.name);
        for id in others {
            assert!(
                out.iter().all(|s| s.topology().slot_of(id).is_none()),
                "{}: the other solid's entity in the result",
                case.name
            );
            assert!(
                !h.relations.iter().any(|r| r.sources().contains(&id)),
                "{}: the other solid's entity in the history",
                case.name
            );
        }
    }
}

/// A third prism strictly inside one solid of several meets none of its
/// faces: the adjacency that sorts the second arrangement's pieces by solid
/// has nothing to go by (the decisions' `ComputationLimit`).
#[test]
fn a_partner_inside_one_of_several_solids_is_a_limit() {
    use rusty_occt::identity::OperationId;
    let case = cases()
        .into_iter()
        .find(|c| c.name == "severed_drilled_cut")
        .unwrap();
    let (_, _, first, _) = protocol::run_first(&case).unwrap();
    let given = protocol::given(&case, first);
    let inner = protocol::build(&protocol::parse(
        "case inner 1e-07\nop 97\nframe 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0\noffsets 1.0 2.0\n\
         boundary P 4 4.0 0.25 6.0 0.25 6.0 1.0 4.0 1.0",
    ));
    for run in [
        given.cut(OperationId(98), &inner),
        given.fuse(OperationId(98), &inner),
        inner.common(OperationId(98), &given),
    ] {
        assert!(
            matches!(run, Err(Error::ComputationLimit(_))),
            "{:?}",
            run.map(|r| r.0.len())
        );
    }
}

/// Declared degenerate classes are refused: a tangency between the stack's
/// rim and the third cylinder, and (the conservative refusal) a tangency
/// between the third box and the solid not given.
#[test]
fn tangencies_are_degenerate() {
    for (name, run) in runs() {
        if name.starts_with("rollex_tangent") || name.starts_with("severed_tangent") {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}: {run:?}");
        }
    }
}

/// The kernel stores the frames the reference sliced, bit for bit.
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-given-frames.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let w: Vec<&str> = row.split('\t').collect();
        let case = cases.iter().find(|c| c.name == w[0]).unwrap();
        let (which, axis) = w[1].split_once(' ').unwrap();
        let spec = match which {
            "obj" => &case.object,
            "tool" => &case.tool,
            _ => &case.then.as_ref().unwrap().third,
        };
        let frame = protocol::build(spec).frame();
        let v = match axis {
            "n" => frame.normal(),
            "x" => frame.x(),
            _ => frame.y(),
        };
        let got = [v.x, v.y, v.z].map(hex).join(" ");
        assert_eq!(got, w[2], "{} {}", w[0], w[1]);
    }
}
