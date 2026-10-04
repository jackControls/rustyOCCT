//! S9d.4b.2a and S9d.4b.2b: Booleans of a whole torus against prisms with
//! arcs, spheres, cones and whole tori, against the independent reference
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
        "tori_coax_fuse",
        "tori_coax_common",
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
            "tori_ring_cut",
            "tori_side_cut",
            "tori_cross_cut",
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
/// the outer equator, a sphere touching the equator at a point and a
/// coaxial torus touching along a circle.
#[test]
fn tangencies_are_degenerate() {
    for (name, run) in runs() {
        if ["pipe_equator_fuse", "ball_touch_fuse", "tori_kiss_fuse"].contains(&name.as_str()) {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}");
        }
    }
}

/// A whole torus of radii 2.5 and 0.75 about an axis parallel to `T`'s
/// (radii 2.5 and 1), its centre on `T`'s core circle (S9d.4b.2b; the
/// evidence's correction (c)): with the top circles at one height (its
/// centre a quarter up, or `T` moved half a unit along x: both of radius
/// 2.5, crossing) or its bottom on `T`'s top the surfaces are tangent where
/// those circles cross (both normals along the axis), `Degenerate`; with
/// the equators at one height (its centre on `T`'s plane) they cross
/// transversally (the normals there along the two radii, 43 degrees apart)
/// and the meeting only turns: the common two solids of 7.74729, as
/// native DRAW's (`bcommon`, `vprops`: 7.74729, the fuse 69.359 and the cut
/// 41.6007, each valid).
#[test]
fn parallel_tori_at_one_height() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Frame3, Point3, RigidTransform, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let turn = std::f64::consts::TAU;
    let (t, _) =
        Solid::torus_with(OperationId(1), Frame3::xy(), 2.5, 1.0, 0.0, turn, turn, tol).unwrap();
    let at = |z: f64| {
        let f = Frame3::new(
            Point3::new(2.5, 0.0, z),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            tol,
        )
        .unwrap();
        Solid::torus_with(OperationId(2), f, 2.5, 0.75, 0.0, turn, turn, tol)
            .unwrap()
            .0
    };
    let moved = Solid::torus_with(OperationId(7), Frame3::xy(), 2.5, 1.0, 0.0, turn, turn, tol)
        .unwrap()
        .0
        .transform_with(
            OperationId(6),
            RigidTransform::translation(Vec3::new(0.5, 0.0, 0.0)).unwrap(),
        )
        .unwrap()
        .0;
    for other in [at(0.25), at(1.75), moved] {
        assert!(
            matches!(t.fuse(OperationId(3), &other), Err(Error::Degenerate(_))),
            "a tangency"
        );
    }
    let (out, h) = t.common(OperationId(5), &at(0.0)).unwrap();
    assert_eq!(out.len(), 2);
    let v: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
    assert!((v - 7.74729).abs() < 5e-6, "{v}");
    let ins = [
        t.topology().entity_set(t.resolution()),
        at(0.0).topology().entity_set(tol),
    ];
    let outs: Vec<_> = out
        .iter()
        .map(|s| s.topology().entity_set(s.resolution()))
        .collect();
    assert!(history::check(&ins, &outs, &h).is_empty());
}

/// Two tori in a turned frame (S9d.4b.2b): the stored axes not exactly
/// orthonormal, the other torus's function is of degree four in each of the
/// first's angles, its points algebraic of degree eight and its turning
/// points found by subdivision; the torus ringing the tube, moved as the
/// results are above, keeps the reference's volume.
#[test]
fn turned_tori_meet_in_fields_of_degree_eight() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::topology::Curve3;
    use rusty_occt::{Point3, RigidTransform, Vec3};
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let case = cases()
        .into_iter()
        .find(|c| c.name == "tori_ring_common")
        .unwrap();
    let (a, b) = (protocol::build(&case.object), protocol::build(&case.tool));
    let a = a.transform_with(OperationId(801), motion).unwrap().0;
    let b = b.transform_with(OperationId(802), motion).unwrap().0;
    let (out, _) = a.common(case.operation, &b).unwrap();
    let (_, want) = &expected()["tori_ring_common"];
    let (count, v) = want.unwrap();
    assert_eq!(out.len(), count);
    let m = out[0].topology().mass_enclosure().unwrap();
    let slack = 1e-9 * v[0];
    assert!(
        m.volume[0] - slack <= v[0] && v[0] <= m.volume[1] + slack,
        "{m:?}"
    );
    assert!(out[0]
        .topology()
        .edges()
        .iter()
        .any(|e| { matches!(&e.curve, Curve3::Toric(t) if t.other_minor > 0.0) }));
}

/// Two tori of equal radii about one centre, the second's frame built from
/// the first's normal normalized again (an ulp off it) or with its axes
/// turned (rounded otherwise): one surface within rounding, `Degenerate`
/// (its meeting's projections were left unpinned, `PrecisionLoss`, before);
/// a thinner torus there, nested within it, evaluates as about the first's
/// own frame.
#[test]
fn tori_within_rounding_of_one_surface_are_degenerate() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Frame3, Point2, Point3, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let turn = std::f64::consts::TAU;
    let tilt = Frame3::new(
        Point3::new(1.0, -2.0, 0.5),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let centre = tilt.point(Point2::new(0.0, 0.0), 1.5);
    let ring = |id: u64, f: Frame3, small: f64| {
        Solid::torus_with(OperationId(id), f, 2.0, small, 0.0, turn, turn, tol)
            .unwrap()
            .0
    };
    let own = Frame3::new(
        centre,
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    assert_eq!(own.normal(), tilt.normal());
    let first = ring(1, own, 0.75);
    // The normal normalized again, an ulp off (where this platform's
    // `hypot` keeps it, turned by an ulp first).
    let hint = tilt.x() * 3.0 + tilt.y() * 4.0;
    let up = |x: f64, k: i64| f64::from_bits(x.to_bits().wrapping_add_signed(k));
    let n = tilt.normal();
    let again = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)]
        .into_iter()
        .map(|(ky, kz)| Frame3::new(centre, Vec3::new(n.x, up(n.y, ky), up(n.z, kz)), hint, tol))
        .find_map(|f| f.ok().filter(|f| f.normal() != n))
        .expect("a frame an ulp off");
    let turned = Frame3::new(centre, Vec3::new(0.0, 3.0, 4.0), tilt.y(), tol).unwrap();
    for f in [again, turned] {
        let other = ring(2, f, 0.75);
        for r in [
            first.fuse(OperationId(3), &other),
            first.cut(OperationId(4), &other),
            first.common(OperationId(5), &other),
        ] {
            assert!(
                matches!(&r, Err(Error::Degenerate(m)) if m.contains("one surface")),
                "{:?}",
                r.map(|x| x.0.len())
            );
        }
    }
    let (near, own) = (ring(2, again, 0.5), ring(2, own, 0.5));
    let volume = |r: Vec<Solid>| r.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let (a, b) = (
        volume(first.cut(OperationId(4), &near).unwrap().0),
        volume(first.cut(OperationId(4), &own).unwrap().0),
    );
    // The tube's ring between the radii: 2 pi^2 R (r1^2 - r2^2).
    let want = 2.0 * std::f64::consts::PI.powi(2) * 2.0 * (0.5625 - 0.25);
    assert!(
        (a - b).abs() <= 1e-9 * a && (a - want).abs() <= 1e-9 * want,
        "{a} {b} {want}"
    );
}
