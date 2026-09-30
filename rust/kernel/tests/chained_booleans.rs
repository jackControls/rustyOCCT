//! S9e.1: a Boolean's result given to another Boolean, against the
//! independent reference (`fixtures/boolean-chained-*` from
//! `tools/generate_chained_boolean_fixtures.py`): each case's first Boolean
//! (S9c.1's, one solid) and its second on that result and a third prism, as
//! object or tool (`swapped`). Each case runs once (on a few threads) for
//! the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::Error;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!("../../fixtures/boolean-chained-cases.txt"))
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
    for line in include_str!("../../fixtures/boolean-chained-expected.tsv")
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

/// The second Boolean's history is complete over its inputs, the first
/// result's entities among them: it chains onto the first Boolean's, every
/// entity of the first result resolved by it.
#[test]
fn chained_histories_are_complete() {
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

/// The first result's ids are the second Boolean's inputs: the history of
/// the first Boolean names every entity the second one reads.
#[test]
fn the_first_result_is_the_second_input() {
    for case in cases() {
        let Ok((a, b, first, h1)) = protocol::run_first(&case) else {
            panic!("{}: the first Boolean", case.name)
        };
        assert_eq!(first.len(), 1, "{}", case.name);
        let ins = [
            a.topology().entity_set(a.resolution()),
            b.topology().entity_set(b.resolution()),
        ];
        let outs = [first[0].topology().entity_set(first[0].resolution())];
        assert!(history::check(&ins, &outs, &h1).is_empty(), "{}", case.name);
        let Some(Ok((x, y, _, _))) = runs()
            .iter()
            .find(|(n, _)| *n == case.name)
            .map(|(_, r)| r.as_ref())
        else {
            continue;
        };
        let given = if case.then.as_ref().unwrap().swapped {
            y
        } else {
            x
        };
        let ids = |s: &rusty_occt::Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
        assert_eq!(ids(given), ids(&first[0]), "{}", case.name);
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

/// A first result moved rigidly, given to the second Boolean with the third
/// prism moved alike, keeps the reference's volumes (the result's model
/// built again from its moved inputs, checked against its moved vertices).
/// The motion is a translation by binary64 steps: a turned motion rounds
/// the frames apart, and the fixtures' planes parallel to another frame's
/// cylinder (a box's wall along the tilted hole) then lie within rounding
/// of its direction, S9c.1's refusal.
#[test]
fn moved_results_keep_their_volumes() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{RigidTransform, Vec3};
    let motion = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let expect = expected();
    for name in [
        "hole_halved_common",
        "hole_tool_cut",
        "groove_drilled_cut",
        "quarter_bored_fuse",
        "pin_slab_cut",
        "hole_capped_fuse",
    ] {
        let case = cases().into_iter().find(|c| c.name == name).unwrap();
        let (_, _, first, _) = protocol::run_first(&case).unwrap();
        let then = case.then.as_ref().unwrap();
        let r = first[0].transform_with(OperationId(801), motion).unwrap().0;
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

/// The groove leaves the box's top face in two result faces on one input
/// face; the drilled hole's faces continue the strip it goes through, the
/// other strip is kept whole.
#[test]
fn a_face_in_parts_names_each_part() {
    let (_, _, first, _) = protocol::run_first(
        &cases()
            .into_iter()
            .find(|c| c.name == "groove_drilled_cut")
            .unwrap(),
    )
    .unwrap();
    let t = first[0].topology();
    let top: Vec<_> = t
        .ids()
        .filter_map(|(id, slot)| match slot {
            rusty_occt::topology::Slot::Face(f) => match &t.faces()[f.index()].surface {
                rusty_occt::topology::Surface::Plane(frame)
                    if frame.normal().z.abs() == 1.0 && frame.origin().z == 4.0 =>
                {
                    Some(id)
                }
                _ => None,
            },
            _ => None,
        })
        .collect();
    assert_eq!(top.len(), 2, "the grooved top in two faces");
    let (_, _, out, h) = runs()
        .iter()
        .find(|(n, _)| n == "groove_drilled_cut")
        .unwrap()
        .1
        .as_ref()
        .unwrap();
    let kept: Vec<_> = top
        .iter()
        .filter(|id| out[0].topology().ids().any(|(x, _)| x == **id))
        .collect();
    // One strip unchanged, the other drilled: modified in place.
    assert_eq!(
        kept.len(),
        2,
        "both strips keep their ids: {:?}",
        h.relations
    );
}

/// Declared degenerate classes are refused.
#[test]
fn tangencies_and_edges_on_faces_are_degenerate() {
    for (name, run) in runs() {
        if name.starts_with("groove_tangent") || name.starts_with("groove_edge") {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}: {run:?}");
        }
    }
}

/// The kernel stores the frames the reference sliced, bit for bit.
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-chained-frames.tsv")
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
