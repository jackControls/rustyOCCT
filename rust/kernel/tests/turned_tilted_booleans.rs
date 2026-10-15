//! S9e.4b.4d.2b's moved turned bodies: imported turned bodies about an axis
//! off the world's (`(1, 2, 2) / 3`, OCCT storing every surface's frame
//! rounded: a waist of corners and a pin whose bottom and step are cones
//! within the resolution of a plane) given to Booleans against an imported
//! box OCCT wrote on a `SKEW` frame and a ball, against the independent
//! reference (`fixtures/boolean-turned-tilted-*` from
//! `tools/generate_turned_tilted_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/tilted_*.brep`): a coaxial cylinder's and cone's
//! meeting on such a frame a plane's section (`cones.rs`'s `coaxial_plane`),
//! where its rulings' quartic over the frame's rounded data met the imported
//! box's rounded planes past S9e.3b's eliminant's budget. Each case runs once
//! (on a few threads) for the checks that read its result.
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
        "../../fixtures/boolean-turned-tilted-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-turned-tilted-expected.tsv")
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
/// their frames again keep the reference's volumes, the imported bodies read
/// again off their moved stored topologies: on the turned frames a coaxial
/// meeting stays a plane's section, and every case evaluates.
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    let mut turned = 0;
    for name in ["waist_box_cut", "pin_box_common", "waist_ball_fuse"] {
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
    assert_eq!(turned, 3, "{turned} turned cases evaluated");
}

fn body(name: &str, op: u64) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(op), topology, resolution).map(|(s, _)| s)
}

/// A profile's volume by Pappus, its heights three times its `w`: its
/// frusta's `pi h (a^2 + a b + b^2) / 3`.
fn pappus(points: &[(f64, f64)]) -> f64 {
    points
        .windows(2)
        .map(|p| {
            let ((a, w0), (b, w1)) = (p[0], p[1]);
            std::f64::consts::PI * (w1 - w0) * (a * a + a * b + b * b)
        })
        .sum()
}

/// The bodies' model profiles `(r, w)` from the axis out and back to it, as
/// the generator's (the pin's bottom and step level).
fn profiles() -> Vec<(&'static str, Vec<(f64, f64)>)> {
    vec![
        (
            "tilted_waist",
            vec![
                (0.0, 0.0),
                (6.0, 0.0),
                (6.0, 1.0),
                (4.0, 2.0),
                (4.0, 3.0),
                (0.0, 3.0),
            ],
        ),
        (
            "tilted_pin",
            vec![
                (0.0, 0.0),
                (6.0, 0.0),
                (6.0, 1.0),
                (4.0, 1.0),
                (3.0, 3.0),
                (0.0, 3.0),
            ],
        ),
    ]
}

/// The point at radius `r`, `turn` radians about the axis from `x` = `(2,
/// 1, -2) / 3` towards `y` = `(-2, 2, -1) / 3`, and height `w` (`3 w` along
/// the axis `(1, 2, 2) / 3` from `(1/2, -1/4, 0)`).
fn at(r: f64, turn: f64, w: f64) -> Point3 {
    let (c, s) = (turn.cos(), turn.sin());
    let x = [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0];
    let y = [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0];
    let a = [1.0, 2.0, 2.0];
    let o = [0.5, -0.25, 0.0];
    let p: [f64; 3] = std::array::from_fn(|i| o[i] + r * (c * x[i] + s * y[i]) + w * a[i]);
    Point3::new(p[0], p[1], p[2])
}

/// A box about the world's axes holding every body.
fn holder(op: u64) -> Solid {
    use rusty_occt::{Boundary, Frame3, Point2, Profile, Tolerance};
    let tol = Tolerance::default();
    let square = [(-16.0, -16.0), (16.0, -16.0), (16.0, 16.0), (-16.0, 16.0)]
        .map(|(x, y)| Point2::new(x, y))
        .to_vec();
    let frame = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let profile = Profile::new(Boundary::polygon(square, tol).unwrap(), Vec::new(), tol).unwrap();
    Solid::extrude_with(OperationId(op), profile, frame, -16.0, 16.0)
        .unwrap()
        .0
}

/// Every body imports as the lattice of its model profile on its tilted
/// axis: its volume and its common with a box holding it its model profile's
/// by Pappus, points within the profile inside and beyond it outside, its
/// stored vertices on its boundary; moved by a turn, alike.
#[test]
fn imported_bodies_are_their_model_profiles() {
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    for (name, model) in profiles() {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let want = pappus(&model);
        let v = s.mass_properties().volume;
        assert!((v - want).abs() <= 1e-9 * want, "{name}: {v} for {want}");
        let (out, _) = s.common(OperationId(95), &holder(94)).unwrap();
        let got: f64 = out.iter().map(|x| x.mass_properties().volume).sum();
        assert_eq!(out.len(), 1, "{name}");
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "{name}: {got} for {want}"
        );
        for p in model.windows(2) {
            let ((a, w0), (b, w1)) = (p[0], p[1]);
            if w1 <= w0 {
                continue;
            }
            let (w, r) = (0.5 * (w0 + w1), 0.5 * (a + b));
            for (k, angle) in [0.3f64, 2.1, 4.4].into_iter().enumerate() {
                let inside = at(0.75 * r, angle, w);
                let outside = at(r + 0.25, angle, w);
                assert_eq!(s.classify(inside).unwrap(), Location::Inside, "{name} {k}");
                assert_eq!(
                    s.classify(outside).unwrap(),
                    Location::Outside,
                    "{name} {k}"
                );
            }
        }
        let (moved, _) = s.transform_with(OperationId(900), turn).unwrap();
        let v = moved.mass_properties().volume;
        assert!((v - want).abs() <= 1e-9 * want, "{name}: moved {v}");
        for x in [&s, &moved] {
            for v in x.topology().vertices() {
                assert_eq!(
                    x.classify(v.position).unwrap(),
                    Location::Boundary,
                    "{name}"
                );
            }
        }
    }
}
