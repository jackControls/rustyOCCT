//! S9d.2c: Booleans of a sphere against cylinders where a frame is turned
//! (a turned cap's circles of unequal axes against a cylinder, a sphere
//! meeting a turned cylinder in a loop), against the independent reference
//! (`fixtures/boolean-spheres-turned-*` from
//! `tools/generate_spheres_turned_boolean_fixtures.py`). Each case runs once
//! (on a few threads) for the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history;
use rusty_occt::Error;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-spheres-turned-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-spheres-turned-expected.tsv")
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

/// A turned cylinder's loop keeps `Curve3::Rise`'s closed form for its
/// pieces over the height: on the stored axes its binary64 points lie on the
/// cylinder and within rounding of the sphere (the exact meeting's function
/// on the stored ellipse is of degree two; the closed form reads it as the
/// circle it rounds).
#[test]
fn a_turned_loops_height_pieces_lie_on_both_surfaces() {
    use rusty_occt::topology::Curve3;
    let loops = [
        "bite_lean_fuse",
        "bite_lean_cut",
        "bite_lean_common",
        "bite_tilt_fuse",
        "bite_tilt_common",
        "graze_tiltx_cut",
        "graze_tiltx_common",
        "bite_r125_cut",
        "dome_bite_turned_cut",
        "dome_bite_turned_common",
    ];
    for (name, run) in runs() {
        if !loops.contains(&name.as_str()) {
            continue;
        }
        let Ok((_, _, out, _)) = run else {
            panic!("{name}: refused")
        };
        let mut rises = 0;
        for e in out.iter().flat_map(|s| s.topology().edges().iter()) {
            let Curve3::Rise(m) = &e.curve else { continue };
            rises += 1;
            for k in 0..=16 {
                let p = m.point(f64::from(k) / 16.0);
                let off = (p - m.centre).length() - m.sphere_radius;
                assert!(off.abs() < 1e-12, "{name}: {off} off the sphere");
                let [x, y, _] = m.frame.coordinates(p);
                let off = x.hypot(y) - m.radius;
                assert!(off.abs() < 1e-12, "{name}: {off} off the cylinder");
            }
        }
        assert!(rises > 0, "{name}: no piece over the height");
    }
}

/// A turned cap's circles meet the cylinder at points of `Q(alpha)`: the
/// hemispheres' results hold vertices on the rim's plane (their binary64
/// views within rounding of it) where the rim crosses the cylinder (not in
/// `dome_tiltx_rings`, whose rings lie on either side of the equator and
/// whose split's great circle crosses the rod).
#[test]
fn a_turned_caps_rim_meets_the_cylinder_on_its_plane() {
    use rusty_occt::Vec3;
    let caps = [
        ("dome_tilt_pipe_cut", Vec3::new(0.0, 3.0, 4.0)),
        ("dome_tilt_pipe_common", Vec3::new(0.0, 3.0, 4.0)),
        ("dome_lean_bite_common", Vec3::new(3.0, 0.0, 4.0)),
    ];
    for (name, run) in runs() {
        let Some((_, n)) = caps.iter().find(|(c, _)| c == name) else {
            continue;
        };
        let Ok((_, _, out, _)) = run else {
            panic!("{name}: refused")
        };
        let n = *n * (1.0 / n.length());
        let on_rim = out
            .iter()
            .flat_map(|s| s.topology().vertices().iter())
            .filter(|v| {
                (v.position - rusty_occt::Point3::new(0.0, 0.0, 0.0))
                    .dot(n)
                    .abs()
                    < 1e-12
            })
            .count();
        assert!(on_rim >= 2, "{name}: {on_rim} vertices on the rim's plane");
    }
}

/// A turned cone's loop against a sphere (S9d.3b's `ball_tilt`) and a turned
/// cap's circle against a cone stay S9d.3c's.
#[test]
fn turned_cones_stay_later() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Frame3, Point3, Solid, Tolerance, Vec3};
    let later = |r: Result<(Vec<rusty_occt::Solid>, history::History), Error>| matches!(r, Err(Error::OutOfDomain(m)) if m.contains("S9d.3c"));
    let ball_tilt = protocol::cases(include_str!("../../fixtures/boolean-cones-cases.txt"))
        .into_iter()
        .find(|c| c.name == "ball_tilt_common")
        .expect("S9d.3b's fixture");
    assert!(later(
        protocol::run(&ball_tilt).map(|(_, _, out, h)| (out, h))
    ));
    let tol = Tolerance::default();
    let tilt = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let half = std::f64::consts::FRAC_PI_2;
    let (dome, _) = Solid::sphere_with(OperationId(1), tilt, 2.0, 0.0, half, tol).unwrap();
    let base = Frame3::new(
        Point3::new(0.25, 0.0, -3.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let (cone, _) = Solid::cone_with(OperationId(2), base, 1.75, 1.5, 6.0, tol).unwrap();
    assert!(later(dome.common(OperationId(3), &cone)));
}
