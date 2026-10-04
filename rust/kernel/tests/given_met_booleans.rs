//! S9e.3b: a Boolean's result whose edges on meetings of two curved faces
//! (`Rise`, `Meet`, `Toric`), or on a cone's or a torus's general plane
//! section, a face of the other input meets, against the independent
//! reference (`fixtures/boolean-given-met-*` from
//! `tools/generate_given_met_boolean_fixtures.py`): each case's Booleans in
//! turn, the second on the first's one solid and the third solid, as object
//! or tool (`swapped`). Each case runs once (on a few threads) for the checks
//! that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Solid};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!("../../fixtures/boolean-given-met-cases.txt"))
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
    for line in include_str!("../../fixtures/boolean-given-met-expected.tsv")
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
/// from its moved construction (a translation by binary64 steps, as
/// S9e.3a's).
#[test]
fn moved_results_keep_their_volumes() {
    use rusty_occt::{RigidTransform, Vec3};
    let motion = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let expect = expected();
    for name in [
        "peg_tilt_cut",
        "peg_wall_common",
        "cross_ball_fuse",
        "torus_rod_wall_fuse",
        "cone_cut_ball_cut",
        "torus_cut_pipe_common",
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

/// Declared degenerate cases are refused: the box touching the peg's
/// meeting with the sphere at its lowest point and the ball through that
/// point with dependent normals (the partner tangent to the given meeting),
/// and the rod whose top cap's plane holds the frustum's apex (S9d.3a's).
#[test]
fn tangencies_are_degenerate() {
    for (name, run) in runs() {
        if name.starts_with("peg_touch") || name.starts_with("peg_kiss") {
            match run {
                Err(Error::Degenerate(m)) => assert!(m.contains("S9e.3b"), "{name}: {m}"),
                other => panic!("{name}: {:?}", other.as_ref().map(|_| ())),
            }
        }
        if name.starts_with("cone_cut_rod") {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}: {run:?}");
        }
    }
}

/// The kernel stores the frames the reference swept, bit for bit.
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-given-met-frames.tsv")
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

fn volume(out: &[Solid]) -> f64 {
    out.iter().map(|s| s.mass_properties().volume).sum()
}

/// The given result of a case's first Boolean.
fn given_of(name: &str) -> Solid {
    let case = cases().into_iter().find(|c| c.name == name).unwrap();
    let (_, _, first, _) = protocol::run_first(&case).unwrap();
    protocol::given(&case, first)
}

/// The three operations with a partner keep the pair identities by the
/// kernel's own measures: `V(X u C) + V(X n C) = V(X) + V(C)` and `V(X - C)
/// = V(X) - V(X n C)`.
fn identities(given: &Solid, partner: &Solid) -> Result<(), String> {
    let f = boolean(given, "fuse", 97, partner).map_err(|e| e.to_string())?;
    let c = boolean(given, "cut", 97, partner).map_err(|e| e.to_string())?;
    let m = boolean(given, "common", 97, partner).map_err(|e| e.to_string())?;
    let (vx, vc) = (
        given.mass_properties().volume,
        partner.mass_properties().volume,
    );
    let near = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0);
    if !near(volume(&f) + volume(&m), vx + vc) || !near(volume(&c), vx - volume(&m)) {
        return Err(format!(
            "fuse {} cut {} common {} given {vx} partner {vc}",
            volume(&f),
            volume(&c),
            volume(&m)
        ));
    }
    Ok(())
}

/// A cone's elliptic section met by a rod (the declared fixture's rod
/// lowered clear of the frustum's apex): a conic against a conic in the
/// section's plane, the three operations consistent.
#[test]
fn a_cone_section_met_by_a_rod() {
    let given = given_of("cone_cut_rod_cut");
    let rod = solid(
        "op 96\nframe 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0\noffsets -1.0 5.5\nboundary C 0.1 2.5 0.4",
    );
    identities(&given, &rod).unwrap();
}

/// Two crossing rods' meeting met by a wall through the thick rod's axis
/// (parallel to its rulings, the carrier's fibre holding two points): the
/// projection retried, the three operations consistent.
#[test]
fn a_meeting_met_by_a_wall_along_the_carriers_rulings() {
    let given = given_of("cross_wall_fuse");
    let wall = solid(
        "op 96\nframe 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0\noffsets -8.0 8.0\n\
         boundary P 4 -8.0 0.25 8.0 0.25 8.0 8.0 -8.0 8.0",
    );
    identities(&given, &wall).unwrap();
}

/// A torus among three curved surfaces is refused by cost (S9e.3b): the
/// torus and the rod's `Toric` met by a ball.
#[test]
fn a_torus_among_three_curved_surfaces_is_refused() {
    let given = given_of("torus_rod_wall_fuse");
    let ball = solid(
        "op 96\nframe 3.49 0.09 0.87 0.0 0.0 1.0 1.0 0.0 0.0\nsphere 0.3 -1.5707963267948966 1.5707963267948966",
    );
    match boolean(&given, "cut", 97, &ball) {
        Err(Error::OutOfDomain(m)) => assert!(m.contains("S9e.3b"), "{m}"),
        other => panic!("{:?}", other.map(|r| r.len())),
    }
}

/// The fuzz replay's chained stage (`fuzz/regressions/README.md`): a bar in
/// the tilted frame less a torus band across it, given to a ball about the
/// middle of its first torus section, which lies on the bar's wall, so the
/// wall holds the ball's axis and its section is a meridian through the
/// ball's pole inside the bar (`uv_gap` at `86b1d834`; `sphere_booleans.rs`'s
/// box and ball alike): the three operations consistent.
#[test]
fn a_ball_about_a_section_on_a_wall_through_its_axis() {
    use rusty_occt::topology::Curve3;
    use rusty_occt::{Boundary, Frame3, Point3, Profile, Tolerance, Vec3};
    let tol = Tolerance::default();
    let frame = |o: Point3| Frame3::new(o, Vec3::new(0.0, 3.0, 4.0), Vec3::new(1.0, 0.0, 0.0), tol);
    let fa = frame(Point3::new(1.0, -2.0, 0.5)).unwrap();
    let fb = frame(Point3::new(
        1.0 + 0.625 / 3.0,
        -2.0 - 0.875 / 3.0,
        0.5 + 2.25 / 5.0,
    ))
    .unwrap();
    let bar = Profile::new(Boundary::rectangle(4.25, 0.5, tol).unwrap(), vec![], tol).unwrap();
    let (a, _) = Solid::extrude_with(OperationId(1), bar, fa, 0.0, 2.25).unwrap();
    let turn = std::f64::consts::TAU;
    let (band, _) =
        Solid::torus_with(OperationId(2), fb, 1.3125, 0.65625, 0.5, 2.25, turn, tol).unwrap();
    let given = boolean(&a, "cut", 4, &band).unwrap().remove(0);
    let at = given
        .topology()
        .edges()
        .iter()
        .find_map(|e| {
            use Curve3::{Meet, Rise, Section, Toric};
            matches!(e.curve, Meet(_) | Rise(_) | Toric(_) | Section(_)).then(|| e.curve.point(0.5))
        })
        .unwrap();
    let f = Frame3::new(at, fa.normal(), fa.x() * 3.0 + fa.y() * 4.0, tol).unwrap();
    let half = std::f64::consts::FRAC_PI_2;
    let (ball, _) = Solid::sphere_with(OperationId(7), f, 1.25, -half, half, tol).unwrap();
    identities(&given, &ball).unwrap();
}
