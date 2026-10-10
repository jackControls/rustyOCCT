//! S9e.4b.2: imported polyhedra other than prisms (bodies of plane faces and
//! line edges OCCT wrote to `.brep` files: wedges, polyhedra of given points,
//! results of boxes), given to Booleans, against the independent reference
//! (`fixtures/boolean-imported-polyhedra-*` from
//! `tools/generate_imported_polyhedra_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`): each case's Boolean (or chain) with its imported
//! inputs made by `Solid::imported_with`, a polyhedron decided on its stored
//! vertices (S9b.2's stored model). Each case runs once (on a few threads)
//! for the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Location, Point3, Solid};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-imported-polyhedra-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-imported-polyhedra-expected.tsv")
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

/// The declared refusal names its reason: the fuse of the slab on the
/// pyramid's base's plane as S9's `Degenerate` (the cavity's cases solid since
/// S9e.4b.4c.1).
#[test]
fn refusals_name_their_reasons() {
    for (name, run) in runs() {
        let reason = match run {
            Err(Error::Degenerate(m) | Error::OutOfDomain(m)) => *m,
            _ => "",
        };
        if name == "pyramid_flush_fuse" {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}: {reason}");
        }
    }
}

/// Each Boolean's history is complete over its inputs, the imported ones'
/// stored ids among them: every entity of every input resolved, and no
/// relation from outside the inputs.
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
    use rusty_occt::{RigidTransform, Vec3};
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

/// Both inputs moved rigidly keep the reference's volumes: translated, and
/// turned where no exact contact between them is broken by the turn's
/// rounding (the moved stored vertices are roundings of the turned ones).
#[test]
fn moved_inputs_keep_their_volumes() {
    use rusty_occt::{RigidTransform, Vec3};
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    for name in [
        "tetra_box_cut",
        "box_tetra_common",
        "octa_slab_cut",
        "pyramid_box_fuse",
        "truncated_box_common",
        "wedge_slab_cut",
        "notched_rod_cut",
        "ell_slab_fuse",
        "tetra_octa_fuse",
        "tetra_octa_cut",
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

/// Every body imports (a polyhedron on its stored vertices; `ridge` as
/// S9e.4a's prism; `hollow` with its cavity since S9e.4b.4c.1), its mass
/// its construction's (the reference's exact volumes), its stored vertices
/// on its boundary and points off it classified.
#[test]
fn imported_bodies_are_their_constructions() {
    let bodies = [
        ("tetra", 56.0, [3.6896, 4.5736, 1.0]),
        ("octa", 256.0 / 3.0, [5.0, 5.0, 4.0]),
        ("pyramid", 128.0, [5.0, 4.6, 2.3]),
        ("truncated", 112.0, [3.8748, 4.0207, 4.0]),
        ("wedge", 368.0 / 3.0, [4.6756, 4.9047, 3.587]),
        ("notched", 394.48066666666665, [4.9369, 4.9531, 1.9886]),
        ("ell", 166.31088913245537, [3.6327, 3.8469, 2.2835]),
        ("steps_low", 632.0 / 3.0, [108.0, -36.0, 6.9304]),
        ("steps_high", 182.0 / 3.0, [108.0, -36.0, 8.4396]),
        ("pedestal", 4072.306558201019, [-1.3912, 1.4243, 6.0012]),
        ("draft", 335.0, [1.99, 6.3657, 3.595]),
        ("ridge", 470.25, [7.0, 6.2081, 3.3134]),
        ("vane_up", 49.0 / 6.0, [5.0, 5.0, 5.2092]),
        ("vane_down", 49.0, [5.0, 5.0, 3.7449]),
        ("hollow", 936.0, [1.5, 5.0, 5.0]),
    ];
    for (name, volume, inside) in bodies {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let m = s.mass_properties();
        assert!(
            (m.volume - volume).abs() <= 1e-9 * volume,
            "{name}: {} for {volume}",
            m.volume
        );
        let [x, y, z] = inside;
        assert_eq!(
            s.classify(Point3::new(x, y, z)).unwrap(),
            Location::Inside,
            "{name}"
        );
        assert_eq!(
            s.classify(Point3::new(x + 500.0, y, z)).unwrap(),
            Location::Outside,
            "{name}"
        );
        assert_eq!(
            s.classify(Point3::new(x, y, z + 50.0)).unwrap(),
            Location::Outside,
            "{name}"
        );
        for v in s.topology().vertices() {
            assert_eq!(
                s.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
    }
    // The hollow box's cavity is outside it.
    let hollow = body("hollow", 91).unwrap();
    assert_eq!(
        hollow.classify(Point3::new(5.0, 5.0, 5.0)).unwrap(),
        Location::Outside
    );
}

/// An imported polyhedron against a body with curved faces, and a result of
/// one given to a Boolean with curved faces, evaluate since S9e.4b.4c.1 (the
/// polyhedron's stored triangles a leaf of the curved engine): the pair
/// identities hold.
#[test]
fn curved_partners_evaluate() {
    use rusty_occt::{Boundary, Frame3, Point2, Profile, Tolerance, Vec3};
    let tolerance = Tolerance::default();
    let rod = |op: u64| {
        let frame = Frame3::new(Point3::new(0.0, 0.0, -1.0), Vec3::Z, Vec3::X, tolerance).unwrap();
        let disc = Boundary::circle(Point2::new(5.0, 4.0), 1.5, tolerance).unwrap();
        let profile = Profile::new(disc, vec![], tolerance).unwrap();
        Solid::extrude_with(OperationId(op), profile, frame, 0.0, 9.0)
            .unwrap()
            .0
    };
    let volume = |out: &[Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let pyramid = body("pyramid", 91).unwrap();
    let (r, p) = (rod(92), pyramid.mass_properties().volume);
    let vr = r.mass_properties().volume;
    let fuse = volume(&r.fuse(OperationId(93), &pyramid).unwrap().0);
    let cut = volume(&pyramid.cut(OperationId(93), &r).unwrap().0);
    let common = volume(&pyramid.common(OperationId(93), &r).unwrap().0);
    assert!(
        (fuse + common - p - vr).abs() <= 1e-9 * p,
        "{fuse} {common}"
    );
    assert!((cut + common - p).abs() <= 1e-9 * p, "{cut} {common}");
    let frame = Frame3::new(Point3::new(0.0, 0.0, 3.0), Vec3::Z, Vec3::X, tolerance).unwrap();
    let square = Boundary::polygon(
        vec![
            Point2::new(2.0, 4.5),
            Point2::new(8.0, 4.5),
            Point2::new(8.0, 10.5),
            Point2::new(2.0, 10.5),
        ],
        tolerance,
    )
    .unwrap();
    let box_ = Solid::extrude_with(
        OperationId(94),
        Profile::new(square, vec![], tolerance).unwrap(),
        frame,
        0.0,
        4.0,
    )
    .unwrap()
    .0;
    let (cut, _) = pyramid.cut(OperationId(95), &box_).unwrap();
    let [first] = <[Solid; 1]>::try_from(cut).unwrap();
    let v = first.mass_properties().volume;
    let cut = volume(&first.cut(OperationId(96), &rod(97)).unwrap().0);
    let common = volume(&first.common(OperationId(96), &rod(97)).unwrap().0);
    assert!((cut + common - v).abs() <= 1e-9 * v, "{cut} {common}");
}

/// The kernel stores the frames the reference swept, bit for bit (the
/// solids built from rows).
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-imported-polyhedra-frames.tsv")
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
