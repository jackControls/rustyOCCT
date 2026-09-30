//! S9d.4b.2a: Booleans of a whole torus against prisms with arcs, spheres
//! and cones, against the independent reference
//! (`fixtures/boolean-torus-curved-*` from
//! `tools/generate_torus_curved_boolean_fixtures.py`). Each case runs once
//! (on a few threads) for the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history;
use rusty_occt::Error;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

/// Two tori meet in fields of degree eight about different axes: S9d.4b.2b's
/// (every `tori_*` case).
const LATER: &[&str] = &[
    "tori_coax_fuse",
    "tori_coax_common",
    "tori_link_fuse",
    "tori_link_common",
    "tori_ring_fuse",
    "tori_ring_cut",
    "tori_ring_common",
    "tori_side_cut",
    "tori_side_common",
    "tori_cross_cut",
    "tori_cross_common",
    "tori_kiss_fuse",
];

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-torus-curved-cases.txt"
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

/// Every case's run, once.
fn runs() -> &'static [(String, Result<protocol::Run, Error>)] {
    static RUNS: OnceLock<Vec<(String, Result<protocol::Run, Error>)>> = OnceLock::new();
    RUNS.get_or_init(|| each(protocol::run))
}

/// The reference's kind and `(solids, [volume, area])` per case.
type Expected = std::collections::BTreeMap<String, (String, Option<(usize, [f64; 2])>)>;

fn expected() -> Expected {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-torus-curved-expected.tsv")
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

#[test]
fn every_case_matches_the_reference() {
    let expect = expected();
    // The later sub-step's cases are exactly the two tori.
    let names: Vec<String> = cases().into_iter().map(|c| c.name).collect();
    for name in &names {
        assert_eq!(
            LATER.contains(&name.as_str()),
            name.starts_with("tori_"),
            "{name}"
        );
    }
    let mut failures = Vec::new();
    for (name, run) in runs() {
        let (kind, want) = &expect[name];
        if LATER.contains(&name.as_str()) {
            match run {
                Err(Error::OutOfDomain(m)) if m.contains("S9d.4b.2b") => {}
                other => failures.push(format!(
                    "{name}: {:?} not S9d.4b.2b's",
                    other.as_ref().map(|_| ())
                )),
            }
            continue;
        }
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
        let mut sums = [0.0f64; 4];
        for s in out {
            let Some(m) = s.topology().mass_enclosure() else {
                failures.push(format!("{name}: a result's mass undecided"));
                continue;
            };
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
            failures.push(format!("{name}: {} solids for {count}", out.len()));
            continue;
        }
        let near = |x: f64, lo: f64, hi: f64| {
            let slack = 1e-9 * x.abs().max(1.0);
            lo - slack <= x && x <= hi + slack
        };
        // The enclosures narrow: a wide one (its midpoint the reported
        // value) would hold the reference yet report another.
        let narrow = |lo: f64, hi: f64| hi - lo <= 1e-9 * lo.abs().max(hi.abs()).max(1.0);
        if !narrow(sums[0], sums[1]) || !narrow(sums[2], sums[3]) {
            failures.push(format!("{name}: wide enclosures {sums:?}"));
        }
        if !near(v[0], sums[0], sums[1]) || !near(v[1], sums[2], sums[3]) {
            failures.push(format!("{name}: volume and area {v:?} against {sums:?}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures: {failures:#?}",
        failures.len()
    );
}

#[test]
fn fixture_histories_are_complete() {
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

/// Coaxial pairs meet in circles about the axis (`Curve3::Circle` or arcs),
/// others in `Curve3::Toric` graphs; no result edge is a plane's spiric
/// section where no plane meets the torus.
#[test]
fn coaxial_pairs_meet_in_circles() {
    use rusty_occt::topology::Curve3;
    let coaxial = [
        "pipe_hole_fuse",
        "sleeve_common",
        "ball_coax_cut",
        "ball_axis_common",
        "cone_coax_common",
        "pipe_tilt_common",
    ];
    for (name, run) in runs() {
        let Ok((_, _, out, _)) = run else { continue };
        let torics = out
            .iter()
            .flat_map(|s| s.topology().edges().iter())
            .filter(|e| matches!(e.curve, Curve3::Toric(_)))
            .count();
        if coaxial.contains(&name.as_str()) {
            assert_eq!(torics, 0, "{name}");
        } else if [
            "rod_cut",
            "bore_cut",
            "ball_top_cut",
            "spike_cut",
            "stadium_cut",
        ]
        .contains(&name.as_str())
        {
            assert!(torics > 0, "{name}");
        }
    }
}

/// The fuzz replay's variant (S9d.4b.2a): a torus about `x` inside a ball of
/// radius 3 but for its top, one section passing within 0.003 of the
/// ball's stored pole; its projection on the sphere needs more anchors than
/// 16 to pin its lift. The common is the torus less the cut.
#[test]
fn a_section_near_a_spheres_pole_projects() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Frame3, Point3, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let half = std::f64::consts::FRAC_PI_2;
    let (a, _) = Solid::sphere_with(OperationId(1), Frame3::xy(), 3.0, -half, half, tol).unwrap();
    let fb = Frame3::new(
        Point3::new(-0.625, -0.625, 0.5),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        tol,
    )
    .unwrap();
    let turn = std::f64::consts::TAU;
    let (t, _) =
        Solid::torus_with(OperationId(2), fb, 1.875, 0.9375, 0.0, turn, turn, tol).unwrap();
    let volume = |out: &[Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let common = volume(&t.common(OperationId(3), &a).unwrap().0);
    let cut = volume(&t.cut(OperationId(4), &a).unwrap().0);
    let vt = t.mass_properties().volume;
    assert!(
        (common + cut - vt).abs() <= 1e-9 * vt,
        "{common} {cut} {vt}"
    );
    assert!(cut > 0.0 && common > cut);
}

/// The fuzz replay's other variant (S9d.4b.2a): a torus in a turned frame
/// through a square slab's round hole, the hole's wall meeting it in curves
/// wound about both the axis and the tube (a torus knot's face, whose loop
/// winding the validator leaves undecided): `ComputationLimit`, not an
/// invalid body.
#[test]
fn a_torus_face_wound_both_ways_is_undecided() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let s = 4.75;
    let square = Boundary::polygon(
        [(-s, -s), (s, -s), (s, s), (-s, s)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec(),
        tol,
    )
    .unwrap();
    let hole = Boundary::circle(Point2::new(0.0, 0.0), 2.375, tol).unwrap();
    let profile = Profile::new(square, vec![hole], tol).unwrap();
    let (a, _) = Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 2.25).unwrap();
    let fb = Frame3::new(
        Point3::new(1.0, 1.0, 1.125),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let turn = std::f64::consts::TAU;
    let (t, _) =
        Solid::torus_with(OperationId(2), fb, 3.5625, 1.78125, 0.0, turn, turn, tol).unwrap();
    assert!(matches!(
        a.common(OperationId(5), &t),
        Err(Error::ComputationLimit(m)) if m.contains("loop winding")
    ));
}

/// Only a tangency of the surfaces is refused: a coaxial cylinder along
/// the outer equator and a sphere touching the equator at a point.
#[test]
fn tangencies_are_degenerate() {
    for (name, run) in runs() {
        if ["pipe_equator_fuse", "ball_touch_fuse"].contains(&name.as_str()) {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}");
        }
    }
}

/// A torus segment or wedge against a curved face stays S9d.4b's
/// `OutOfDomain`; two tori name S9d.4b.2b.
#[test]
fn parts_and_tori_pairs_stay_later() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Frame3, Solid, Tolerance};
    use std::f64::consts::{FRAC_PI_2, TAU};
    let tol = Tolerance::default();
    let torus = |op: u64, low: f64, high: f64, angle: f64| {
        Solid::torus_with(
            OperationId(op),
            Frame3::xy(),
            2.5,
            1.0,
            low,
            high,
            angle,
            tol,
        )
        .unwrap()
        .0
    };
    let (ball, _) = Solid::sphere_with(
        OperationId(9),
        Frame3::xy(),
        1.25,
        -FRAC_PI_2,
        FRAC_PI_2,
        tol,
    )
    .unwrap();
    let half = torus(1, -FRAC_PI_2, FRAC_PI_2, TAU);
    match half.common(OperationId(3), &ball) {
        Err(Error::OutOfDomain(m)) => assert!(m.contains("segment or wedge"), "{m}"),
        other => panic!("{:?}", other.map(|_| ())),
    }
    let (a, b) = (torus(4, 0.0, TAU, TAU), torus(5, 0.0, TAU, TAU));
    let b = b
        .transform_with(
            OperationId(6),
            rusty_occt::RigidTransform::translation(rusty_occt::Vec3::new(0.5, 0.0, 0.0)).unwrap(),
        )
        .unwrap()
        .0;
    match a.fuse(OperationId(7), &b) {
        Err(Error::OutOfDomain(m)) => assert!(m.contains("S9d.4b.2b"), "{m}"),
        other => panic!("{:?}", other.map(|_| ())),
    }
}
