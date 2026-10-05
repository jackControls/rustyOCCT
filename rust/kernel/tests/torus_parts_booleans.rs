//! S9d.4c: Booleans of a sphere's cap or zone against a whole torus, and of
//! torus v-segments and wedges against prisms with arcs, spheres, cones and
//! whole tori, against the independent reference
//! (`fixtures/boolean-torus-parts-*` from
//! `tools/generate_torus_parts_boolean_fixtures.py`). Each case runs once (on
//! a few threads) for the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history;
use rusty_occt::Error;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!("../../fixtures/boolean-torus-parts-cases.txt"))
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
    for line in include_str!("../../fixtures/boolean-torus-parts-expected.tsv")
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
    // The enclosures narrow: a wide one (its midpoint the reported value)
    // would hold the reference yet report another.
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

/// A cap's rim within the resolution of tangency to the torus and a half's
/// rim tangent to a pipe are refused.
#[test]
fn tangencies_are_degenerate() {
    for (name, run) in runs() {
        if ["rim_touch_fuse", "oh_rim_touch_common"].contains(&name.as_str()) {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}");
        }
    }
}

/// Coaxial pairs meet in circles about the axis, others in `Curve3::Toric`
/// graphs; the pairs across a part's rim of surd radius, or a cap's, have a
/// vertex on it.
#[test]
fn coaxial_pairs_meet_in_circles() {
    use rusty_occt::topology::Curve3;
    let coaxial = [
        "zone_coax_fuse",
        "cap_coax_common",
        "oh_pipe_common",
        "ih_pipe_common",
        "ih_cone_common",
        "ih_ball_common",
        "oh_dome_common",
        "band_tilt_pipe_common",
    ];
    let meeting = [
        "cap_top_cut",
        "zone_side_cut",
        "dome_lean_cut",
        "oh_ball_cut",
        "oh_torus_cut",
        "band_ball_cut",
        "band_cone_cut",
        "band_torus_cut",
        "qw_pipe_cut",
        "qw_torus_cut",
        "hw_ball_cut",
        "tw_cone_cut",
        "qw_lean_ball_cut",
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
        } else if meeting.contains(&name.as_str()) {
            assert!(torics > 0, "{name}");
        }
    }
}

/// Parts and caps in a turned frame (S9d.4c): a band, a quarter wedge and a
/// cap moved as the results are above, against their tools moved alike,
/// keep the reference's volumes (their rims' crossings in fields of the
/// turned frame's rounded axes).
#[test]
fn turned_parts_keep_their_volumes() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Point3, RigidTransform, Vec3};
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let expect = expected();
    for name in ["band_ball_common", "qw_pipe_common", "cap_top_common"] {
        let case = cases().into_iter().find(|c| c.name == name).unwrap();
        let (a, b) = (protocol::build(&case.object), protocol::build(&case.tool));
        let a = a.transform_with(OperationId(801), motion).unwrap().0;
        let b = b.transform_with(OperationId(802), motion).unwrap().0;
        let (out, _) = a.common(case.operation, &b).unwrap();
        let (count, v) = expect[name].1.unwrap();
        within(&out, count, v).unwrap_or_else(|e| panic!("{name}: {e}"));
    }
}

/// The fuzz replay's first variant (S9d.4c): a box over a band's hole, just
/// above its upper end disc, fused with it: the void between the disc, the
/// band's inner wall and the box's bottom is a cavity of the fuse (a shell of
/// both inputs' faces turned inward), whose containment the validator's
/// rays leave undecided: `ComputationLimit`, not an invalid solid; the cut
/// and the common evaluate.
#[test]
fn a_void_under_a_bands_hole_is_a_cavity() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let square = Boundary::polygon(
        [(-3.0, -3.0), (3.0, -3.0), (3.0, 3.0), (-3.0, 3.0)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec(),
        tol,
    )
    .unwrap();
    let profile = Profile::new(square, vec![], tol).unwrap();
    let (b, _) = Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.5567, 2.0).unwrap();
    let f = Frame3::new(
        Point3::new(0.625, 1.125, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let turn = std::f64::consts::TAU;
    let (t, _) =
        Solid::torus_with(OperationId(2), f, 1.125, 0.703125, 0.5, 2.25, turn, tol).unwrap();
    assert!(matches!(
        b.fuse(OperationId(3), &t),
        Err(Error::ComputationLimit(m)) if m.contains("cavity")
    ));
    let volume =
        |out: &[rusty_occt::Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let cut = volume(&b.cut(OperationId(4), &t).unwrap().0);
    let common = volume(&b.common(OperationId(5), &t).unwrap().0);
    let vb = b.mass_properties().volume;
    assert!(
        (cut + common - vb).abs() <= 1e-9 * vb,
        "{cut} {common} {vb}"
    );
}

/// The fuzz replay's second variant (S9d.4c): a quarter wedge stood on its
/// side, its end half-plane on the rounded direction of its turn within
/// rounding of a sphere's poles; the sphere face's loop through the pole
/// winds by its pcurves' lifted ends (half a turn at the pole), and the
/// three operations keep the volumes' identities.
#[test]
fn a_loop_through_a_spheres_pole_winds_by_its_ends() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Frame3, Point3, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let half = std::f64::consts::FRAC_PI_2;
    let (a, _) =
        Solid::sphere_with(OperationId(1), Frame3::xy(), 3.1875, -half, half, tol).unwrap();
    let f = Frame3::new(
        Point3::new(0.875, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        tol,
    )
    .unwrap();
    let turn = std::f64::consts::TAU;
    let (t, _) =
        Solid::torus_with(OperationId(2), f, 2.0625, 1.2890625, 0.0, turn, half, tol).unwrap();
    let volume =
        |out: &[rusty_occt::Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let fuse = volume(&a.fuse(OperationId(3), &t).unwrap().0);
    let cut = volume(&a.cut(OperationId(4), &t).unwrap().0);
    let common = volume(&a.common(OperationId(5), &t).unwrap().0);
    let (va, vt) = (a.mass_properties().volume, t.mass_properties().volume);
    assert!((fuse - (va + vt - common)).abs() <= 1e-9 * (va + vt));
    assert!((cut - (va - common)).abs() <= 1e-9 * va);
    assert!(common > 0.0 && common < vt);
}

/// A frustum against a torus band in the tilted frame and its offset (a
/// `TURNED_PARTS` variant of the fuzz replay): on the frustum's wall a hole
/// whose point no piece's 24-sided polygon holds, near one of them, is
/// refused as a sliver within the resolution, as before the wall's loops
/// were sampled again where a hole's point lies near another loop's (finer
/// polygons nested it in the wall's piece, and the result was open).
#[test]
fn a_hole_no_piece_holds_is_refused() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Frame3, Point3, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let (n, x) = (Vec3::new(0.0, 3.0, 4.0), Vec3::new(1.0, 0.0, 0.0));
    let f = Frame3::new(Point3::new(1.0, -2.0, 0.5), n, x, tol).unwrap();
    let g = Frame3::new(
        Point3::new(1.0 + 1.0 / 3.0, -2.0 + 1.0 / 3.0, 0.5 + 2.25 / 5.0),
        n,
        x,
        tol,
    )
    .unwrap();
    let (a, _) = Solid::cone_with(OperationId(1), f, 3.5625, 1.78125, 2.25, tol).unwrap();
    let turn = std::f64::consts::TAU;
    let (b, _) =
        Solid::torus_with(OperationId(2), g, 3.5625, 2.2265625, -2.5, -0.25, turn, tol).unwrap();
    for r in [
        a.fuse(OperationId(3), &b),
        a.cut(OperationId(4), &b),
        a.common(OperationId(5), &b),
    ] {
        match r {
            Err(Error::Degenerate(m)) => assert!(m.contains("thinner"), "{m}"),
            other => panic!("{:?}", other.map(|r| r.0.len())),
        }
    }
}
