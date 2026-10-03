//! S9e.3a: a Boolean's result of spheres, cones and tori (or with procedural
//! edges the partner does not reach), deeper chains, and given results
//! against a sphere, a cone or a torus, against the independent reference
//! (`fixtures/boolean-given-curved-*` from
//! `tools/generate_given_curved_boolean_fixtures.py`): each case's Booleans
//! in turn, each on the previous one's result's one solid and the next
//! solid, as object or tool (`swapped`). Each case runs once (on a few
//! threads) for the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Solid};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-given-curved-cases.txt"
    ))
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

/// Every case's run (the last Boolean's inputs, results and history), once.
fn runs() -> &'static [(String, Result<protocol::Run, Error>)] {
    static RUNS: OnceLock<Vec<(String, Result<protocol::Run, Error>)>> = OnceLock::new();
    RUNS.get_or_init(|| each(protocol::run))
}

/// The reference's kind and `(solids, [volume, area])` per case.
type Expected = std::collections::BTreeMap<String, (String, Option<(usize, [f64; 2])>)>;

fn expected() -> Expected {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-given-curved-expected.tsv")
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
fn within(out: &[Solid], count: usize, v: [f64; 2]) -> Result<(), String> {
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

/// Every Boolean of a case in turn, each on the previous result's one solid:
/// its inputs, results and history.
fn stages(case: &protocol::Case) -> Result<Vec<protocol::Run>, Error> {
    let first = protocol::run_first(case)?;
    let mut out = vec![first];
    for then in case.then.iter().chain(&case.more) {
        let prev = out.last().unwrap().2.clone();
        let given = protocol::given_by(&case.name, then, prev);
        let next = protocol::build(&then.third);
        let (a, b) = if then.swapped {
            (next, given)
        } else {
            (given, next)
        };
        let (res, h) = match then.op.as_str() {
            "fuse" => a.fuse(then.operation, &b),
            "cut" => a.cut(then.operation, &b),
            _ => a.common(then.operation, &b),
        }?;
        out.push((a, b, res, h));
    }
    Ok(out)
}

/// Each Boolean's history is complete over its inputs, the previous
/// result's solid among them: the chain's histories chain stage by stage,
/// every entity of every given solid resolved by the next.
#[test]
fn chained_histories_are_complete() {
    let failures: Vec<String> = each(|case| {
        let Ok(all) = stages(case) else {
            return Vec::new();
        };
        let mut bad = Vec::new();
        for (k, (a, b, out, h)) in all.iter().enumerate() {
            let ins = [
                a.topology().entity_set(a.resolution()),
                b.topology().entity_set(b.resolution()),
            ];
            let outs: Vec<_> = out
                .iter()
                .map(|s| s.topology().entity_set(s.resolution()))
                .collect();
            let issues = history::check(&ins, &outs, h);
            if !issues.is_empty() {
                bad.push(format!("{} stage {k}: {issues:?}", case.name));
            }
            for s in [a, b] {
                for (id, _) in s.topology().ids() {
                    if matches!(h.resolve(id), Resolution::Unknown)
                        && !h.relations.iter().any(|r| r.sources().contains(&id))
                    {
                        bad.push(format!("{} stage {k}: an input entity unnamed", case.name));
                    }
                }
            }
            // The next stage's given solid is this stage's result's.
            if let Some((x, y, _, _)) = all.get(k + 1) {
                let then = if k == 0 {
                    case.then.as_ref().unwrap()
                } else {
                    &case.more[k - 1]
                };
                let given = if then.swapped { y } else { x };
                let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
                if !out.iter().any(|s| ids(s) == ids(given)) {
                    bad.push(format!("{} stage {k}: the given solid", case.name));
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

#[test]
fn results_are_deterministic_and_move_rigidly() {
    use rusty_occt::{Location, Point3, RigidTransform, Vec3};
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
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

/// A given result moved rigidly, given to the last Boolean with the last
/// solid moved alike, keeps the reference's volumes: its model built again
/// from its moved construction. The motion is a translation by binary64
/// steps, as S9e.1's: a turned motion rounds the frames apart.
#[test]
fn moved_results_keep_their_volumes() {
    use rusty_occt::{RigidTransform, Vec3};
    let motion = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let expect = expected();
    for name in [
        "dome_tilt_cut",
        "g9_clear_cut",
        "countersink_drill_common",
        "groove_drill_fuse",
        "holes_tilt_common",
        "dome_deep_cut",
        "holed_ball_common",
        "holed_torus_cut",
    ] {
        let case = cases().into_iter().find(|c| c.name == name).unwrap();
        let all = stages(&case).unwrap_or_else(|e| panic!("{name}: {e}"));
        let last = case.more.last().or(case.then.as_ref()).unwrap();
        let (x0, y0, _, _) = &all[all.len() - 1];
        let given = if last.swapped { y0 } else { x0 };
        let r = given.transform_with(OperationId(801), motion).unwrap().0;
        let c = protocol::build(&last.third)
            .transform_with(OperationId(802), motion)
            .unwrap()
            .0;
        let (x, y) = if last.swapped { (&c, &r) } else { (&r, &c) };
        let out = match last.op.as_str() {
            "fuse" => x.fuse(last.operation, y),
            "cut" => x.cut(last.operation, y),
            _ => x.common(last.operation, y),
        }
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .0;
        let (count, v) = expect[name].1.unwrap();
        within(&out, count, v).unwrap_or_else(|e| panic!("{name}: {e}"));
    }
}

/// Declared degenerate cases are refused: DRAW's G9 (the third cylinder
/// touching the frustum's top circle) and a box touching the dome's top.
#[test]
fn tangencies_are_degenerate() {
    for (name, run) in runs() {
        let g9 = name.strip_prefix("g9_").is_some_and(|op| !op.contains('_'));
        if g9 || name.starts_with("dome_touch") {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}: {run:?}");
        }
    }
}

fn solid(rows: &str) -> Solid {
    protocol::build(&protocol::parse(&format!("case extra 1e-07\n{rows}")))
}

fn boolean(a: &Solid, op: &str, id: u64, b: &Solid) -> Result<Vec<Solid>, Error> {
    let id = OperationId(id);
    match op {
        "fuse" => a.fuse(id, b),
        "cut" => a.cut(id, b),
        _ => a.common(id, b),
    }
    .map(|r| r.0)
}

/// A construction tree deeper than three Booleans is a limit: the domed box
/// drilled twice more, then cut (`chain::MAX_DEPTH`).
#[test]
fn deeper_trees_are_a_limit() {
    let box_ = solid(
        "op 91\nframe 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0\noffsets 0.0 4.0\n\
         boundary P 4 0.0 0.0 10.0 0.0 10.0 10.0 0.0 10.0",
    );
    let ball = solid("op 92\nframe 5.0 5.0 4.0 0.0 0.0 1.0 1.0 0.0 0.0\nsphere 3.0 -1.5707963267948966 1.5707963267948966");
    let drill = |x: f64, op: u64| {
        solid(&format!(
            "op {op}\nframe 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0\noffsets -1.0 8.0\nboundary C {x} 5.0 0.5"
        ))
    };
    // A slab clear of the dome's faces, across the first drill.
    let slab = solid(
        "op 99\nframe 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0\noffsets 1.5 2.5\n\
         boundary P 4 0.5 0.5 1.6 0.5 1.6 9.5 0.5 9.5",
    );
    let one = |v: Vec<Solid>| <[Solid; 1]>::try_from(v).ok().unwrap()[0].clone();
    let domed = one(boolean(&box_, "fuse", 93, &ball).unwrap());
    let d1 = one(boolean(&domed, "cut", 95, &drill(1.25, 94)).unwrap());
    let d2 = one(boolean(&d1, "cut", 97, &drill(9.0, 96)).unwrap());
    // Three Booleans: given once more.
    let kept = boolean(&d2, "common", 98, &slab);
    assert!(kept.is_ok(), "{:?}", kept.map(|r| r.len()));
    let d3 = one(boolean(&d2, "cut", 101, &drill(3.0, 100)).unwrap());
    // Four: a limit.
    assert!(matches!(
        boolean(&d3, "common", 102, &slab),
        Err(Error::ComputationLimit(_))
    ));
}

/// A partner face meeting a given meeting of two curved faces (S9e.3b):
/// the sphere fused with the peg, its `Rise` crossed by a box.
#[test]
fn a_meeting_of_curved_faces_met_is_s9e3b() {
    let case = cases()
        .into_iter()
        .find(|c| c.name == "peg_clear_cut")
        .unwrap();
    let (_, _, first, _) = protocol::run_first(&case).unwrap();
    let given = protocol::given(&case, first);
    let cross = solid(
        "op 97\nframe 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0\noffsets 2.5 7.0\n\
         boundary P 4 -4.0 -4.0 4.0 -4.0 4.0 4.0 -4.0 4.0",
    );
    match boolean(&given, "cut", 98, &cross) {
        Err(Error::OutOfDomain(m)) => assert!(m.contains("S9e.3b"), "{m}"),
        other => panic!("{:?}", other.map(|r| r.len())),
    }
}

/// The kernel stores the frames the reference swept, bit for bit.
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-given-curved-frames.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let w: Vec<&str> = row.split('\t').collect();
        let case = cases.iter().find(|c| c.name == w[0]).unwrap();
        let (which, axis) = w[1].split_once(' ').unwrap();
        let k: usize = which.trim_start_matches("solid").parse().unwrap();
        let spec = match k {
            0 => &case.object,
            1 => &case.tool,
            2 => &case.then.as_ref().unwrap().third,
            _ => &case.more[k - 3].third,
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
