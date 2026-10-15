//! S9e.4b.4d.2a: imported turned bodies of cylinders, cones and planes about
//! one axis meeting at corners only (a point, `bug28626`'s kind, about the
//! world's `z` and `x`; a waist; steps in and out) given to Booleans,
//! against the independent reference (`fixtures/boolean-turned-profiles-*`
//! from `tools/generate_turned_profiles_boolean_fixtures.py`, the bodies
//! under `fixtures/imported/profile_*.brep`). The kernel's own turned bodies
//! are no evidence: its writer refuses their coaxial meetings (a cone's
//! section stored as a procedural edge, D13), as it does the fuzz targets'.
//! Each case runs once (on a few threads) for the checks that read its
//! result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Location, Point3, RigidTransform, Solid, Vec3};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-turned-profiles-cases.txt"
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

/// The reference's kind, reason and `(solids, [volume, area])` per case.
type Expected = BTreeMap<String, (String, String, Option<(usize, [f64; 2])>)>;

fn expected() -> Expected {
    let mut expect = BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-turned-profiles-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let entry = expect
            .entry(name.to_string())
            .or_insert((String::new(), String::new(), None));
        match w[0] {
            "expect" => entry.0 = w[1].to_string(),
            "reason" => entry.1 = row["reason ".len()..].to_string(),
            "result" => {
                let v: Vec<f64> = w[2..4].iter().map(|x| x.parse().unwrap()).collect();
                entry.2 = Some((w[1].parse::<usize>().unwrap(), [v[0], v[1]]));
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

/// Every case as declared: solid within the reference's enclosures, the
/// degenerate and unsupported refused for their declared reasons.
#[test]
fn every_case_matches_the_reference() {
    let expect = expected();
    let mut failures = Vec::new();
    for (name, run) in runs() {
        let (kind, reason, want) = &expect[name];
        match (kind.as_str(), run) {
            ("degenerate", Err(Error::Degenerate(m)))
            | ("unsupported", Err(Error::OutOfDomain(m))) => {
                if m != reason {
                    failures.push(format!("{name}: refused as {m}, not {reason}"));
                }
                continue;
            }
            ("degenerate" | "unsupported", other) => {
                failures.push(format!(
                    "{name}: {:?} not refused as {kind}",
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

/// Each Boolean's history is complete over its inputs, the imported ones'
/// stored ids among them: every entity of every input resolved, and no
/// relation from an id outside the inputs.
#[test]
fn histories_are_complete_over_the_imported_ids() {
    let failures: Vec<String> = runs()
        .iter()
        .filter_map(|(name, run)| run.as_ref().ok().map(|r| (name, r)))
        .flat_map(|(name, (a, b, out, h))| {
            let mut bad = Vec::new();
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
                bad.push(format!("{name}: {issues:?}"));
            }
            let inputs: std::collections::BTreeSet<_> = [a, b]
                .iter()
                .flat_map(|s| s.topology().ids().map(|(id, _)| id))
                .collect();
            for s in [a, b] {
                for (id, _) in s.topology().ids() {
                    if matches!(h.resolve(id), Resolution::Unknown)
                        && !h.relations.iter().any(|r| r.sources().contains(&id))
                    {
                        bad.push(format!("{name}: an input entity unnamed"));
                    }
                }
            }
            for r in &h.relations {
                if r.sources().iter().any(|id| !inputs.contains(id)) {
                    bad.push(format!("{name}: a relation from outside the inputs"));
                }
            }
            bad
        })
        .collect();
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
fn results_are_deterministic_and_move_rigidly() {
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    let first: BTreeMap<&str, _> = runs().iter().map(|(n, r)| (n.as_str(), r)).collect();
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

/// Both inputs moved by a dyadic translation and by a rotation that rounds
/// their frames keep the reference's volumes (a turn moving an exact
/// incidence within rounding is refused as `Degenerate`, a given meeting's
/// eliminant over the turned frames' data past its budget as
/// `ComputationLimit`, S9e.3b's): the imported bodies read again off their
/// moved stored topologies, their cones on their moved axes exactly. On a
/// turned frame a coaxial cylinder's and cone's meeting is a quartic graph
/// (`Curve3::Meet`, the frame's stored axes orthonormal within rounding),
/// and a partner's turned plane meeting it is S9e.3b's triple past its
/// budget: the ball's case evaluates turned, the others' shift alone.
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    let mut turned = 0;
    for name in [
        "point_box_cut",
        "waist_rod_common",
        "steps_ball_fuse",
        "box_crown_cut",
        "point_tilt_cut",
    ] {
        let case = all.iter().find(|c| c.name == name).unwrap();
        for (k, motion) in [shift, turn].into_iter().enumerate() {
            let (a, b, _, _) = protocol::run_first(case).unwrap();
            let a = a.transform_with(OperationId(801), motion).unwrap().0;
            let b = b.transform_with(OperationId(802), motion).unwrap().0;
            let run = match case.op.as_str() {
                "fuse" => a.fuse(case.operation, &b),
                "cut" => a.cut(case.operation, &b),
                _ => a.common(case.operation, &b),
            };
            let (count, [v, _]) = expect[name].2.unwrap();
            match run {
                Err(Error::Degenerate(_) | Error::ComputationLimit(_)) if k == 1 => continue,
                Err(e) => panic!("{name} {k}: {e}"),
                Ok((out, _)) => {
                    let got: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
                    assert_eq!(out.len(), count, "{name} {k}");
                    assert!((got - v).abs() <= 1e-9 * v, "{name} {k}: {got} for {v}");
                    turned += k;
                }
            }
        }
    }
    assert!(turned >= 1, "{turned} turned cases evaluated");
}

fn body(name: &str, op: u64) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(op), topology, resolution).map(|(s, _)| s)
}

/// A profile's volume by Pappus: its frusta's `pi h (a^2 + a b + b^2) / 3`.
fn pappus(points: &[(f64, f64)]) -> f64 {
    points
        .windows(2)
        .map(|p| {
            let ((a, w0), (b, w1)) = (p[0], p[1]);
            std::f64::consts::PI * (w1 - w0) * (a * a + a * b + b * b) / 3.0
        })
        .sum()
}

/// The bodies' profiles `(r, w)` from the axis out and back to it, their
/// axes' origins and whether along the world's `x` (else `z`).
type Turned = (&'static str, Vec<(f64, f64)>, [f64; 3], bool);

fn profiles() -> Vec<Turned> {
    let point = vec![
        (0.0, 0.0),
        (3.0, 3.0),
        (3.0, 6.0),
        (4.0, 8.5),
        (4.0, 16.0),
        (0.0, 16.0),
    ];
    vec![
        ("profile_point", point.clone(), [0.0, 0.0, 0.0], false),
        ("profile_point_x", point, [-1.5, 4.125, 0.375], true),
        (
            "profile_waist",
            vec![
                (0.0, 0.0),
                (8.0, 0.0),
                (8.0, 10.0),
                (5.0, 14.0),
                (5.0, 20.0),
                (0.0, 20.0),
            ],
            [0.0, 0.0, 0.0],
            false,
        ),
        (
            "profile_steps",
            vec![
                (0.0, 0.0),
                (6.0, 0.0),
                (6.0, 4.0),
                (5.0, 5.0),
                (3.0, 5.0),
                (3.0, 9.0),
                (4.0, 11.0),
                (4.0, 13.0),
                (0.0, 13.0),
            ],
            [0.0, 0.0, 0.0],
            false,
        ),
        (
            "profile_crown",
            vec![
                (0.0, 0.0),
                (3.0, 0.0),
                (3.0, 4.0),
                (6.0, 4.0),
                (6.0, 7.0),
                (4.0, 9.0),
                (4.0, 11.0),
                (0.0, 11.0),
            ],
            [1.0, 2.0, 0.0],
            false,
        ),
    ]
}

/// A body's point at radius `r` along a direction `turn` radians about its
/// axis from the frame's `x` and height `w`.
fn at((_, _, o, along_x): &Turned, r: f64, turn: f64, w: f64) -> Point3 {
    let (c, s) = (turn.cos(), turn.sin());
    if *along_x {
        Point3::new(o[0] + w, o[1] + r * c, o[2] + r * s)
    } else {
        Point3::new(o[0] + r * c, o[1] + r * s, o[2] + w)
    }
}

/// Every body imports as the lattice of its primitives (its model matched
/// on import), its volume its profile's by Pappus, its stored vertices on
/// its boundary, points within its profile inside and points beyond it
/// outside; moved by a turn, its volume and its vertices' places alike.
#[test]
fn imported_bodies_are_their_profiles() {
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    for t in profiles() {
        let (name, points, _, _) = &t;
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let volume = pappus(points);
        let m = s.mass_properties();
        assert!(
            (m.volume - volume).abs() <= 1e-9 * volume,
            "{name}: {} for {volume}",
            m.volume
        );
        for p in points.windows(2) {
            let ((a, w0), (b, w1)) = (p[0], p[1]);
            if w1 <= w0 {
                continue;
            }
            let w = 0.5 * (w0 + w1);
            let r = 0.5 * (a + b);
            for (k, angle) in [0.3, 2.1, 4.4].into_iter().enumerate() {
                let inside = at(&t, r * 0.75, angle, w);
                let outside = at(&t, r + 0.25, angle, w);
                assert_eq!(s.classify(inside).unwrap(), Location::Inside, "{name} {k}");
                assert_eq!(
                    s.classify(outside).unwrap(),
                    Location::Outside,
                    "{name} {k}"
                );
            }
        }
        for v in s.topology().vertices() {
            assert_eq!(
                s.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
        let (moved, _) = s.transform_with(OperationId(900), turn).unwrap();
        let v = moved.mass_properties().volume;
        assert!((v - volume).abs() <= 1e-9 * volume, "{name}: moved {v}");
        for v in moved.topology().vertices() {
            assert_eq!(
                moved.classify(v.position).unwrap(),
                Location::Boundary,
                "{name} moved"
            );
        }
    }
}

/// Profiles this step does not take keep their refusals: a cylinder OCCT
/// stores as two faces on one line (`turned_split`, a vertex on its wall),
/// a torus at corners (`turned_corner`); a profile S9e.4b.4b.1's chains take
/// (`turned_step`) keeps its chain.
#[test]
fn other_profiles_keep_their_outcomes() {
    const NO_CHAIN: &str =
        "an imported body of several primitives other than a Boolean chain of them (S9e.4b.4c)";
    const GENERAL: &str = "an imported solid other than a prism, a sphere, a cone, a torus or a \
                           Boolean of its primitives (S9e.4b)";
    match body("turned_split", 91) {
        Err(Error::OutOfDomain(m)) => assert_eq!(m, NO_CHAIN),
        other => panic!("turned_split: {:?}", other.map(|_| ())),
    }
    match body("turned_corner", 91) {
        Err(Error::OutOfDomain(m)) => assert_eq!(m, GENERAL),
        other => panic!("turned_corner: {:?}", other.map(|_| ())),
    }
    let s = body("turned_step", 91).unwrap();
    let volume = std::f64::consts::PI * (180.0 + 0.001 * 76.0 / 3.0 + 16.0 * 6.999);
    let v = s.mass_properties().volume;
    assert!((v - volume).abs() <= 1e-9 * volume, "{v} for {volume}");
}
