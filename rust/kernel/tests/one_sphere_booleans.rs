//! S9e.4b.3c.1: two pieces of one sphere given to a Boolean (faces of both
//! inputs on one sphere, in general position: S9c.1's faces on one surface,
//! the circles of both crossing on it), and an imported sphere piece whose
//! rim OCCT divided into two arcs (the DRAW survey's `so1` and `so4`),
//! against the independent reference (`fixtures/boolean-one-sphere-*` from
//! `tools/generate_one_sphere_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`). Each case runs once (on a few threads) for the
//! checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Location, Point2, Point3, Profile, RigidTransform, Solid, Tolerance,
    Vec3,
};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!("../../fixtures/boolean-one-sphere-cases.txt"))
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
    for line in include_str!("../../fixtures/boolean-one-sphere-expected.tsv")
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
            ("unsupported", Err(Error::OutOfDomain(_))) => continue,
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

/// The declared refusals name their reasons: a ball within the resolution
/// of the cap's sphere as S9's. The hemisphere and the octant on one frame,
/// S9e.4b.3c.2's exact incidences, evaluate since its kernel.
#[test]
fn refusals_name_their_reasons() {
    for (name, run) in runs() {
        let reason = match run {
            Err(Error::Degenerate(m) | Error::OutOfDomain(m)) => *m,
            _ => "",
        };
        if name.starts_with("ball_near") {
            assert_eq!(
                reason, "two spheres within the resolution of one sphere",
                "{name}"
            );
        }
        if name.starts_with("hemi_octant") {
            assert!(run.is_ok(), "{name}: {reason}");
        }
    }
}

/// Each Boolean's history is complete over its inputs, the imported ones'
/// stored ids among them: every entity of every input resolved, and no
/// relation from an id outside the inputs (a primitive's or a hull's).
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

/// Both inputs moved rigidly keep the reference's volumes, translated and
/// turned: each piece read again off its moved stored topology.
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    for name in [
        "hemi_cap_cut",
        "cap_octant_common",
        "octant_cap_cut",
        "ball_cap_cut",
        "wedge_cap_fuse",
    ] {
        let case = all.iter().find(|c| c.name == name).unwrap();
        for (k, motion) in [shift, turn].into_iter().enumerate() {
            let (a, b, _, _) = protocol::run_first(case).unwrap();
            let a = a.transform_with(OperationId(801), motion).unwrap().0;
            let b = b.transform_with(OperationId(802), motion).unwrap().0;
            let out = match case.op.as_str() {
                "fuse" => a.fuse(case.operation, &b),
                "cut" => a.cut(case.operation, &b),
                _ => a.common(case.operation, &b),
            }
            .unwrap_or_else(|e| panic!("{name} {k}: {e}"))
            .0;
            let (count, [v, _]) = expect[name].1.unwrap();
            let got: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
            assert_eq!(out.len(), count, "{name} {k}");
            assert!((got - v).abs() <= 1e-9 * v, "{name} {k}: {got} for {v}");
        }
    }
}

/// A body file's imported solid under an operation.
fn body(name: &str, op: u64) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(op), topology, resolution).map(|(s, _)| s)
}

/// The divided rims: the hemisphere and the cap import as pieces whose rim
/// is the stored two arcs (two stored vertices on the rim, a third at the
/// stored pole), their volumes their closed forms, every stored vertex on
/// the boundary.
#[test]
fn divided_rims_import_as_pieces() {
    use std::f64::consts::PI;
    for (name, volume) in [
        ("sphere_hemi", 2.0 * PI * 125.0 / 3.0),
        ("sphere_cap", PI * 3.5 * 3.5 * (15.0 - 3.5) / 3.0),
    ] {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let t = s.topology();
        assert_eq!(
            (t.faces().len(), t.edges().len(), t.vertices().len()),
            (2, 2, 3),
            "{name}"
        );
        let m = s.mass_properties();
        assert!(
            (m.volume - volume).abs() <= 1e-9 * volume,
            "{name}: {} for {volume}",
            m.volume
        );
        for v in t.vertices() {
            assert_eq!(
                s.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
        assert_eq!(
            s.classify(Point3::new(5.0, 5.0, 4.0) + Vec3::new(4.0, 1.0, 8.0) * (4.0 / 9.0))
                .unwrap(),
            Location::Inside,
            "{name}"
        );
    }
}

/// A stored vertex on a section circle's axis (the stored pole above a
/// hemisphere's rim on the world's axes, exactly) has no place on the
/// circle: the divided rim's split passes it by (it divided by zero, a
/// panic on import, before S9e.4b.3c.2's evidence found it).
#[test]
fn a_stored_pole_above_its_rim_imports() {
    use std::f64::consts::PI;
    let s = body("incidence_hemi", 91).unwrap();
    let t = s.topology();
    assert_eq!(
        (t.faces().len(), t.edges().len(), t.vertices().len()),
        (2, 2, 3)
    );
    let volume = 2.0 * PI * 125.0 / 3.0;
    let m = s.mass_properties();
    assert!((m.volume - volume).abs() <= 1e-9 * volume, "{}", m.volume);
    assert_eq!(
        s.classify(Point3::new(5.0, 5.0, 9.0)).unwrap(),
        Location::Boundary
    );
}

/// A block less a leaning rod across its back face (a groove: its one
/// curved face's material outside the cylinder) is no piece of its
/// primitive common its planes: refused on import as S9e.4b.3c's before its
/// model was built until S9e.4b.3c.3a, its hull less its primitive since (its
/// volume the kernel's own notch's).
#[test]
fn a_notch_is_its_hull_less_its_primitive() {
    let tol = Tolerance::default();
    let frame = Frame3::xy();
    let (block, _) = Solid::extrude_with(
        OperationId(1),
        Profile::new(
            Boundary::polygon(
                vec![
                    Point2::new(0.0, 0.0),
                    Point2::new(10.0, 0.0),
                    Point2::new(10.0, 7.0),
                    Point2::new(0.0, 7.0),
                ],
                tol,
            )
            .unwrap(),
            Vec::new(),
            tol,
        )
        .unwrap(),
        frame,
        0.0,
        5.0,
    )
    .unwrap();
    // A rod leaning across the block's back face: the groove's wall is no
    // prism wall of the block's direction (not S9e.4a's prism).
    let lean = Frame3::new(
        Point3::new(5.0, 6.5, 2.5),
        Vec3::new(1.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        tol,
    )
    .unwrap();
    let (rod, _) = Solid::extrude_with(
        OperationId(2),
        Profile::new(
            Boundary::circle(Point2::new(0.0, 0.0), 1.0, tol).unwrap(),
            Vec::new(),
            tol,
        )
        .unwrap(),
        lean,
        -6.0,
        6.0,
    )
    .unwrap();
    let (out, _) = block.cut(OperationId(3), &rod).unwrap();
    let [notch] = <[Solid; 1]>::try_from(out).unwrap();
    let (imported, _) =
        Solid::imported_with(OperationId(4), notch.topology().clone(), notch.resolution()).unwrap();
    let (v0, v1) = (
        notch.mass_properties().volume,
        imported.mass_properties().volume,
    );
    assert!((v0 - v1).abs() <= 1e-9 * v0, "{v1} for {v0}");
}

/// Two whole balls of one sphere on different frames (their splits' great
/// circles only on the sphere): their fuse and common the ball, their cut
/// empty; a ball and an imported cap of its sphere, the cap's ids kept.
#[test]
fn two_balls_of_one_sphere() {
    let tol = Tolerance::default();
    let half = std::f64::consts::FRAC_PI_2;
    let c = Point3::new(5.0, 5.0, 4.0);
    let upright = Frame3::new(c, Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0), tol).unwrap();
    let (a, _) = Solid::sphere_with(OperationId(1), upright, 5.0, -half, half, tol).unwrap();
    let turned = Frame3::new(c, Vec3::new(4.0, 1.0, 8.0), Vec3::new(-7.0, -4.0, 4.0), tol).unwrap();
    let (b, _) = Solid::sphere_with(OperationId(2), turned, 5.0, -half, half, tol).unwrap();
    let ball = 4.0 * std::f64::consts::PI * 125.0 / 3.0;
    for (op, want) in [("fuse", Some(ball)), ("common", Some(ball)), ("cut", None)] {
        let out = match op {
            "fuse" => a.fuse(OperationId(5), &b),
            "cut" => a.cut(OperationId(5), &b),
            _ => a.common(OperationId(5), &b),
        }
        .unwrap_or_else(|e| panic!("{op}: {e}"))
        .0;
        match want {
            None => assert!(out.is_empty(), "{op}"),
            Some(v) => {
                assert_eq!(out.len(), 1, "{op}");
                let got = out[0].mass_properties().volume;
                assert!((got - v).abs() <= 1e-9 * v, "{op}: {got}");
            }
        }
    }
}
