//! S9e.4b.4a: imported prisms whose arcs of two circles meet at a joint
//! (bodies OCCT wrote to `.brep` files: four discs' common, a pointed arch,
//! an arc tangent inside another, two circles crossing at a small angle),
//! given to Booleans, against the independent reference
//! (`fixtures/boolean-imported-joints-*` from
//! `tools/generate_imported_joints_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`): each case's Boolean (or chain) with its imported
//! inputs made by `Solid::imported_with`, decided on S9e.4a's construction
//! with each arc ending at such a joint taken through its two ends. Each
//! case runs once (on a few threads) for the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Location, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3,
};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-imported-joints-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-imported-joints-expected.tsv")
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

/// The declared refusals name their reasons: the split lens's joint of two
/// circles each holding other arcs as S9e.4b.4's, the rod on the quad's
/// stored circle as two surfaces within rounding of one (S9's).
#[test]
fn refusals_name_their_reasons() {
    for (name, run) in runs() {
        let reason = match run {
            Err(Error::Degenerate(m) | Error::OutOfDomain(m)) => *m,
            _ => "",
        };
        if name.starts_with("split_lens") {
            assert!(
                matches!(run, Err(Error::OutOfDomain(_)))
                    && reason
                        == "an imported prism's joint of two circles each holding several arcs (S9e.4b.4)",
                "{name}: {reason}"
            );
        }
        if name.starts_with("quad_seat") {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}: {reason}");
        }
    }
}

/// Each Boolean's history is complete over its inputs, the imported ones'
/// stored ids among them: every entity of every input resolved, and no
/// relation from a construction's own ids.
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
    use rusty_occt::RigidTransform;
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
/// turned where no face of one input is exactly parallel to the other's
/// cylinder (a turn rounds such a pair into S9's `Degenerate` "a plane
/// within rounding of a cylinder's direction"); an imported prism's arcs
/// are taken through their ends again in the moved frame.
#[test]
fn moved_inputs_keep_their_volumes() {
    use rusty_occt::RigidTransform;
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    for (name, turned) in [
        ("quad_box_cut", false),
        ("box_quad_common", false),
        ("quad_rod_cut", true),
        ("quad_ball_fuse", true),
        ("arch_box_common", false),
        ("arch_slab_cut", false),
        ("cam_box_fuse", false),
        ("cam_rod_cut", true),
        ("blade_box_cut", false),
        ("both_fuse", true),
        ("both_cut", true),
    ] {
        let case = all.iter().find(|c| c.name == name).unwrap();
        let motions = if turned {
            vec![shift, turn]
        } else {
            vec![shift]
        };
        for (k, motion) in motions.into_iter().enumerate() {
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

/// A circular segment's area beyond its chord: `r^2 (t - sin t) / 2` for
/// the arc's angle `t`.
fn segment(r: f64, t: f64) -> f64 {
    r * r * (t - t.sin()) / 2.0
}

/// Every body imports as S9e.4a's construction (its joints taken through
/// only in a Boolean's exact model), its mass its profile's closed form (a
/// polygon and its arcs' segments) times its height, its stored vertices
/// on its boundary and points off it classified.
#[test]
fn imported_bodies_are_their_constructions() {
    use std::f64::consts::{FRAC_PI_2, PI};
    let theta = 2.0 * 0.6f64.asin();
    let blade_a = 2.0 * 0.28f64.asin();
    let blade_b = FRAC_PI_2 - (2.5f64 / 6.0).atan();
    let bodies = [
        // Four discs' common: the square and four segments of angle 2 asin(3/5).
        (
            "quad",
            4.0 * (36.0 + 4.0 * segment(5.0, theta)),
            [5.0, 5.0, 2.0],
        ),
        // A triangle of base 4 and height 4, two segments of angle atan(4/3).
        (
            "arch",
            3.0 * (8.0 + 2.0 * segment(5.0, (4.0f64 / 3.0).atan())),
            [4.25, 4.3, 1.0],
        ),
        (
            "cam",
            3.0 * (35.5 + segment(5.0, PI - (4.0f64 / 3.0).atan()) + segment(3.0, FRAC_PI_2)),
            [5.0, 5.0, 1.5],
        ),
        (
            "blade",
            3.0 * (120.0 + segment(12.5, blade_a) + segment(6.5, blade_b)),
            [5.0, 3.0, 1.5],
        ),
        (
            "split_lens",
            4.0 * (25.0 * 2.0 * 0.9272952180016122 - 24.0),
            [5.0, 3.0, 2.5],
        ),
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
            s.classify(Point3::new(x + 50.0, y, z)).unwrap(),
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
}

/// The kernel's own prisms whose arcs of two circles meet at rational
/// points (four discs' common, a lens, an S curve, an arc tangent inside
/// another), each in a turned frame, written by the kernel's writer, read
/// back and imported (their joints rounded off both circles by the frames
/// normalized again), give each Boolean with a box and a rod the kernel's
/// own results' volumes within `1e-9`.
#[test]
fn kernel_written_joints_import_as_themselves() {
    use rusty_occt::occt_brep::{import, read, write};
    let tol = Tolerance::new(1e-7, 1e-9).unwrap();
    let arc = |cx: f64, cy: f64, r: f64, ccw: bool| Segment::Arc {
        center: Point2::new(cx, cy),
        radius: r,
        ccw,
    };
    let frame = |o: (f64, f64, f64), x: (f64, f64)| {
        Frame3::new(
            Point3::new(o.0, o.1, o.2),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(x.0, x.1, 0.0),
            tol,
        )
        .unwrap()
    };
    let prism = |pts: &[(f64, f64)], segs: Vec<Segment>, f: Frame3, h: f64| {
        let pts = pts.iter().map(|&(x, y)| Point2::new(x, y)).collect();
        let b = Boundary::path(pts, segs, tol).unwrap();
        Solid::extrude_with(
            OperationId(1),
            Profile::new(b, vec![], tol).unwrap(),
            f,
            0.0,
            h,
        )
        .unwrap()
        .0
    };
    let upright = |outline: Boundary, z0: f64, z1: f64| {
        let f = frame((0.0, 0.0, 0.0), (1.0, 0.0));
        Solid::extrude_with(
            OperationId(2),
            Profile::new(outline, vec![], tol).unwrap(),
            f,
            z0,
            z1,
        )
        .unwrap()
        .0
    };
    let square = |x0: f64, y0: f64, x1: f64, y1: f64| {
        Boundary::polygon(
            vec![
                Point2::new(x0, y0),
                Point2::new(x1, y0),
                Point2::new(x1, y1),
                Point2::new(x0, y1),
            ],
            tol,
        )
        .unwrap()
    };
    let (t30, r125, tilt) = ((0.8660254037844386, 0.5), (12.0, 5.0), (0.6, 0.8));
    let bodies = [
        prism(
            &[(3.0, -3.0), (3.0, 3.0), (-3.0, 3.0), (-3.0, -3.0)],
            vec![
                arc(-1.0, 0.0, 5.0, true),
                arc(0.0, -1.0, 5.0, true),
                arc(1.0, 0.0, 5.0, true),
                arc(0.0, 1.0, 5.0, true),
            ],
            frame((5.0, 5.0, 0.0), r125),
            4.0,
        ),
        prism(
            &[(0.0, -4.0), (0.0, 4.0)],
            vec![arc(-3.0, 0.0, 5.0, true), arc(3.0, 0.0, 5.0, true)],
            frame((5.0, 5.0, 0.0), tilt),
            4.0,
        ),
        prism(
            &[
                (-4.0, -3.0),
                (4.0, -3.0),
                (4.0, 1.0),
                (2.0, 3.0),
                (0.0, 5.0),
                (-4.0, 5.0),
            ],
            vec![
                Segment::Line,
                Segment::Line,
                arc(2.0, 1.0, 2.0, true),
                arc(2.0, 5.0, 2.0, false),
                Segment::Line,
                Segment::Line,
            ],
            frame((5.0, 3.0, 0.0), t30),
            3.0,
        ),
        prism(
            &[(-3.0, -4.0), (5.0, 0.0), (2.0, 3.0), (-3.0, 3.0)],
            vec![
                arc(0.0, 0.0, 5.0, true),
                arc(2.0, 0.0, 3.0, true),
                Segment::Line,
                Segment::Line,
            ],
            frame((5.0, 5.0, 0.0), r125),
            3.0,
        ),
    ];
    let tools = [
        upright(square(6.0, 6.0, 12.0, 12.0), 1.0, 2.5),
        upright(
            Boundary::circle(Point2::new(7.5, 4.75), 1.0, tol).unwrap(),
            -3.0,
            8.0,
        ),
    ];
    let volume = |r: Result<(Vec<Solid>, rusty_occt::history::History), Error>| {
        r.map(|(out, _)| out.iter().map(|s| s.mass_properties().volume).sum::<f64>())
    };
    for (k, solid) in bodies.iter().enumerate() {
        let text = write(solid.topology(), solid.resolution().linear()).unwrap();
        let doc = read(&text).unwrap();
        let [stored] = <[_; 1]>::try_from(import(&doc).solids).unwrap();
        let (again, _) =
            Solid::imported_with(OperationId(3), stored.result.unwrap(), stored.tolerance).unwrap();
        for tool in &tools {
            for op in 0..3 {
                let run = |s: &Solid| match op {
                    0 => volume(s.fuse(OperationId(5), tool)),
                    1 => volume(s.cut(OperationId(5), tool)),
                    _ => volume(s.common(OperationId(5), tool)),
                };
                let (direct, via) = (run(solid), run(&again));
                match (&direct, &via) {
                    (Ok(x), Ok(y)) => assert!(
                        (x - y).abs() <= 1e-9 * x.abs().max(1.0),
                        "body {k} op {op}: {x} {y}"
                    ),
                    _ => panic!("body {k} op {op}: {direct:?} {via:?}"),
                }
            }
        }
    }
}

/// The kernel stores the frames the reference swept, bit for bit (the
/// solids built from rows; the imported ones' are the converter's).
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-imported-joints-frames.tsv")
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
