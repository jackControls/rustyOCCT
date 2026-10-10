//! S9d.1: Booleans of a sphere against polyhedral prisms, against the
//! independent reference (`fixtures/boolean-sphere-*` from
//! `tools/generate_sphere_boolean_fixtures.py`).
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

/// None left to a later sub-step.
const LATER: &[&str] = &[];

#[test]
fn every_case_matches_the_reference() {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-sphere-expected.tsv")
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
    let mut failures = Vec::new();
    for case in protocol::cases(include_str!("../../fixtures/boolean-sphere-cases.txt")) {
        let (kind, want) = &expect[&case.name];
        let rows = match protocol::rows(&case) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{}: {e}", case.name));
                continue;
            }
        };
        if LATER.contains(&case.name.as_str()) {
            if rows != ["unsupported"] {
                failures.push(format!("{}: {rows:?} not unsupported", case.name));
            }
            continue;
        }
        match kind.as_str() {
            "degenerate" => {
                if rows != ["refused"] {
                    failures.push(format!("{}: {rows:?} not refused", case.name));
                }
                continue;
            }
            "empty" => {
                if rows != ["empty"] {
                    failures.push(format!("{}: {rows:?} not empty", case.name));
                }
                continue;
            }
            _ => {}
        }
        let (count, v) = want.expect("a result");
        let solids: Vec<Vec<f64>> = rows
            .iter()
            .map(|r| {
                r.split(' ')
                    .skip(1)
                    .take(4)
                    .map(|x| x.parse().unwrap_or(f64::NAN))
                    .collect()
            })
            .collect();
        if solids.len() != count || solids.iter().any(|s| s.len() < 4) {
            failures.push(format!("{}: {rows:?} for {count} solids", case.name));
            continue;
        }
        let sum = |i: usize| solids.iter().map(|s| s[i]).sum::<f64>();
        let near = |x: f64, lo: f64, hi: f64| {
            let slack = 1e-9 * x.abs().max(1.0);
            lo - slack <= x && x <= hi + slack
        };
        // The enclosures narrow: a wide one (its midpoint the reported
        // value) would hold the reference yet report another.
        let narrow = |lo: f64, hi: f64| hi - lo <= 1e-9 * lo.abs().max(hi.abs()).max(1.0);
        if !narrow(sum(0), sum(1)) || !narrow(sum(2), sum(3)) {
            failures.push(format!(
                "{}: wide enclosures [{}, {}], [{}, {}]",
                case.name,
                sum(0),
                sum(1),
                sum(2),
                sum(3)
            ));
        }
        if !near(v[0], sum(0), sum(1)) || !near(v[1], sum(2), sum(3)) {
            failures.push(format!(
                "{}: volume {v:?} against [{}, {}], area [{}, {}]",
                case.name,
                sum(0),
                sum(1),
                sum(2),
                sum(3)
            ));
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-sphere-cases.txt")) {
        let Ok((a, b, out, h)) = protocol::run(&case) else {
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
        let issues = history::check(&ins, &outs, &h);
        assert!(issues.is_empty(), "{}: {issues:?}", case.name);
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-sphere-cases.txt")) {
        let Ok((_, _, out, h)) = protocol::run(&case) else {
            continue;
        };
        let (_, _, again, h2) = protocol::run(&case).unwrap();
        assert_eq!(h, h2, "{}", case.name);
        for (x, y) in out.iter().zip(&again) {
            assert_eq!(ids(x), ids(y), "{}", case.name);
        }
        for s in &out {
            let (moved, _) = s
                .transform_with(OperationId(900), motion)
                .unwrap_or_else(|e| panic!("{}: {e}", case.name));
            assert_eq!(ids(&moved), ids(s), "{}", case.name);
            let (v0, v1) = (s.mass_properties().volume, moved.mass_properties().volume);
            assert!(
                (v0 - v1).abs() <= 1e-9 * v0.abs().max(1.0),
                "{}: {v0} {v1}",
                case.name
            );
            for v in moved.topology().vertices() {
                assert_eq!(
                    moved.classify(v.position).unwrap(),
                    Location::Boundary,
                    "{}",
                    case.name
                );
            }
        }
    }
}

/// A tilted U whose corner and bottom cap pass through a whole sphere's
/// centre (the `boolean` target's): the fused sphere face is the sphere
/// less a patch clear of its poles, one loop running as a hole (the
/// surface less its loops, as a torus's); the volumes agree.
#[test]
fn a_sphere_less_a_patch_clear_of_its_poles() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let (s, t) = (2.75, 0.625);
    let outer = Boundary::polygon(
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(3.0 * s, 0.0),
            Point2::new(3.0 * s, 2.0 * t),
            Point2::new(2.0 * s, 2.0 * t),
            Point2::new(2.0 * s, t),
            Point2::new(s, t),
            Point2::new(s, 2.0 * t),
            Point2::new(0.0, 2.0 * t),
        ],
        tol,
    )
    .unwrap();
    let profile = Profile::new(outer, vec![], tol).unwrap();
    let frame = Frame3::new(
        Point3::new(1.0, -2.0, 0.5),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let (a, _) = Solid::extrude_with(OperationId(1), profile, frame, 0.0, 1.0).unwrap();
    let half = std::f64::consts::FRAC_PI_2;
    let (b, _) = Solid::sphere_with(OperationId(2), frame, 1.125, -half, half, tol).unwrap();
    let volume = |out: &[Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let fused = a.fuse(OperationId(3), &b).unwrap().0;
    let common = a.common(OperationId(4), &b).unwrap().0;
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let (vf, vc) = (volume(&fused), volume(&common));
    assert!(
        (vf - (va + vb - vc)).abs() <= 1e-9 * vf,
        "{vf} {va} {vb} {vc}"
    );
}

/// A box's wall through an upright whole sphere's centre, its section a
/// great circle through the stored sphere's poles (the DRAW grids'
/// `ZI4`-`ZI7` alike): the section has vertices at the poles and meridians'
/// pcurves; the volumes are the half balls', histories complete, results
/// moved rigidly. A plane through a pole off the axis is `OutOfDomain`.
#[test]
fn sections_through_a_spheres_poles() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{
        history, Boundary, Error, Frame3, Location, Point2, Point3, Profile, RigidTransform, Solid,
        Tolerance, Vec3,
    };
    let tol = Tolerance::default();
    let prism = |frame: Frame3, pts: &[(f64, f64)], h: f64| {
        let outer =
            Boundary::polygon(pts.iter().map(|p| Point2::new(p.0, p.1)).collect(), tol).unwrap();
        let profile = Profile::new(outer, vec![], tol).unwrap();
        Solid::extrude_with(OperationId(1), profile, frame, 0.0, h)
            .unwrap()
            .0
    };
    let up = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let a = prism(
        up,
        &[(0.0, -4.0), (4.0, -4.0), (4.0, 4.0), (0.0, 4.0)],
        16.0,
    );
    let centre = Frame3::new(
        Point3::new(0.0, 0.0, 8.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let half = std::f64::consts::FRAC_PI_2;
    let (b, _) = Solid::sphere_with(OperationId(2), centre, 2.0, -half, half, tol).unwrap();
    let ball = 32.0 * std::f64::consts::PI / 3.0;
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ins = [
        a.topology().entity_set(a.resolution()),
        b.topology().entity_set(b.resolution()),
    ];
    for (r, want) in [
        (a.fuse(OperationId(3), &b), 512.0 + ball / 2.0),
        (a.cut(OperationId(3), &b), 512.0 - ball / 2.0),
        (a.common(OperationId(3), &b), ball / 2.0),
    ] {
        let (out, h) = r.unwrap();
        let v: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
        assert!((v - want).abs() <= 1e-9 * want, "{v} {want}");
        let outs: Vec<_> = out
            .iter()
            .map(|s| s.topology().entity_set(s.resolution()))
            .collect();
        assert!(history::check(&ins, &outs, &h).is_empty());
        for s in &out {
            let (moved, _) = s.transform_with(OperationId(900), motion).unwrap();
            for v in moved.topology().vertices() {
                assert_eq!(moved.classify(v.position).unwrap(), Location::Boundary);
            }
        }
    }
    // A wall through the north pole, off the axis.
    let side = Frame3::new(
        Point3::new(0.0, -4.0, 10.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let t = prism(
        side,
        &[(2.0, -8.0), (6.0, -8.0), (6.0, 16.0), (-4.0, 16.0)],
        8.0,
    );
    assert!(matches!(
        t.fuse(OperationId(3), &b),
        Err(Error::OutOfDomain(_))
    ));
}

/// A box's wall through a sphere's axis, one pole inside the box and the
/// other beyond a cap's plane (the fuzz replay's chained ball about a
/// meeting on a prism's wall, `fuzz/regressions/README.md`): the sphere
/// face's loop runs up its meridian through the pole and down the other
/// side and closes along the cap's circle, winding none. Its pcurves from
/// the pole on are lifted by a turn so the loop closes on its first fin
/// (S9e.4b.3a set its winding to none without lifting them, `uv_gap`
/// unless the loop happened to start at the pole: 76 of these 80 at
/// `86b1d834`). The volumes are the closed forms' and the histories
/// complete, for either pole and whatever the sphere's reference direction.
#[test]
fn a_meridian_loop_through_one_pole_closes_on_its_first_fin() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{history, Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    let half = std::f64::consts::FRAC_PI_2;
    let outer = Boundary::polygon(
        [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec(),
        tol,
    )
    .unwrap();
    let profile = Profile::new(outer, vec![], tol).unwrap();
    let (a, _) = Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 2.0).unwrap();
    let r = 1.25;
    // Half the ball less its cap of height 0.75 beyond z = 0 or z = 2.
    let ball = 4.0 * pi * r * r * r / 3.0;
    let common = (ball - pi * 0.75 * 0.75 * (3.0 * r - 0.75) / 3.0) / 2.0;
    for z in [0.5, 1.5] {
        for x in [(1.0, 0.0), (3.0, 4.0), (0.0, 1.0), (-3.0, 4.0), (-1.0, 0.0)] {
            for up in [1.0, -1.0] {
                let f = Frame3::new(
                    Point3::new(2.0, 0.0, z),
                    Vec3::new(0.0, 0.0, up),
                    Vec3::new(x.0, x.1, 0.0),
                    tol,
                )
                .unwrap();
                let (b, _) = Solid::sphere_with(OperationId(2), f, r, -half, half, tol).unwrap();
                let ins = [
                    a.topology().entity_set(a.resolution()),
                    b.topology().entity_set(b.resolution()),
                ];
                for (out, want) in [
                    (a.fuse(OperationId(3), &b), 16.0 + ball - common),
                    (a.cut(OperationId(4), &b), 16.0 - common),
                    (a.common(OperationId(5), &b), common),
                    (b.cut(OperationId(6), &a), ball - common),
                ] {
                    let (out, h) = out.unwrap_or_else(|e| panic!("{z} {x:?} {up}: {e}"));
                    let v: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
                    assert!(
                        (v - want).abs() <= 1e-9 * want,
                        "{z} {x:?} {up}: {v} {want}"
                    );
                    let outs: Vec<_> = out
                        .iter()
                        .map(|s| s.topology().entity_set(s.resolution()))
                        .collect();
                    assert!(history::check(&ins, &outs, &h).is_empty());
                }
            }
        }
    }
}

/// A box's vertical edge along an upright sphere's axis (through its
/// stored poles), whatever the sphere's reference direction, and the DRAW
/// grids' `ZI5` (a sphere turned a quarter turn about x, then about y, its
/// section through its poles a lune across its seam): the volumes are the
/// closed forms' and their enclosures narrow (a closing chord at a pole,
/// its ends' `v` apart by rounding, was enclosed over its `u` hull, the
/// midpoint far off).
#[test]
fn wedges_through_a_spheres_poles() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{
        Boundary, Frame3, Point2, Point3, Profile, RigidTransform, Solid, Tolerance, Vec3,
    };
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    let half = std::f64::consts::FRAC_PI_2;
    let ball = 32.0 * pi / 3.0;
    let prism = |pts: &[(f64, f64)], h: f64| {
        let outer =
            Boundary::polygon(pts.iter().map(|p| Point2::new(p.0, p.1)).collect(), tol).unwrap();
        let profile = Profile::new(outer, vec![], tol).unwrap();
        let up = Frame3::new(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            tol,
        )
        .unwrap();
        Solid::extrude_with(OperationId(1), profile, up, 0.0, h)
            .unwrap()
            .0
    };
    let check = |out: Vec<Solid>, want: f64| {
        let mut v = 0.0;
        for s in &out {
            let e = s.topology().mass_enclosure().unwrap();
            assert!(e.volume[1] - e.volume[0] <= 1e-9 * want, "{:?}", e.volume);
            v += s.mass_properties().volume;
        }
        assert!((v - want).abs() <= 1e-9 * want, "{v} {want}");
    };
    let quadrant = prism(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)], 16.0);
    for x in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(-1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 1.0, 0.0),
    ] {
        let c = Frame3::new(Point3::new(0.0, 0.0, 8.0), Vec3::new(0.0, 0.0, 1.0), x, tol).unwrap();
        let (b, _) = Solid::sphere_with(OperationId(2), c, 2.0, -half, half, tol).unwrap();
        check(quadrant.common(OperationId(3), &b).unwrap().0, ball / 4.0);
        check(b.cut(OperationId(3), &quadrant).unwrap().0, 0.75 * ball);
    }
    let block = prism(&[(-4.0, -4.0), (4.0, -4.0), (4.0, 4.0), (-4.0, 4.0)], 8.0);
    let c = Frame3::new(
        Point3::new(0.0, 0.0, 8.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let (b, _) = Solid::sphere_with(OperationId(2), c, 2.0, -half, half, tol).unwrap();
    let centre = Point3::new(0.0, 0.0, 8.0);
    let rx = RigidTransform::rotation(centre, Vec3::new(1.0, 0.0, 0.0), half).unwrap();
    let ry = RigidTransform::rotation(centre, Vec3::new(0.0, 1.0, 0.0), half).unwrap();
    let (b, _) = b.transform_with(OperationId(5), rx).unwrap();
    let (b, _) = b.transform_with(OperationId(6), ry).unwrap();
    check(block.common(OperationId(3), &b).unwrap().0, ball / 2.0);
    check(
        block.fuse(OperationId(3), &b).unwrap().0,
        512.0 + ball / 2.0,
    );
    check(block.cut(OperationId(3), &b).unwrap().0, 512.0 - ball / 2.0);
    check(b.cut(OperationId(3), &block).unwrap().0, ball / 2.0);
}

/// A plane crossing a sphere within the resolution of tangency is
/// `Degenerate` (the `boolean` target's
/// `fuzz/regressions/boolean/crash-d238291d9edbe60570ae31479609845872d763c0.bin`):
/// the tilted stadium's flat wall, 1.25 from the centre of a sphere of
/// radius 1.25 about the turned box's frame but for rounding (4e-16 inside
/// it), met it in a circle of radius 1e-8 the validator refused as a
/// degenerate curve, for the stadium and for its cut given again. A box's
/// wall `2^-24` inside an upright sphere is refused alike, and `2^-24`
/// outside it too (a gap within the resolution, the ball inside the box:
/// `a_ball_within_the_resolution_of_a_plane_face_is_degenerate`); `2^-16`
/// inside or outside it evaluate.
#[test]
fn a_wall_crossing_within_the_resolution_of_tangency_is_degenerate() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{
        Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3,
    };
    let tol = Tolerance::default();
    let half = std::f64::consts::FRAC_PI_2;
    let stadium = |s: f64, t: f64| {
        let outer = Boundary::path(
            vec![
                Point2::new(0.0, -t),
                Point2::new(s, -t),
                Point2::new(s, t),
                Point2::new(0.0, t),
            ],
            vec![
                Segment::Line,
                Segment::Arc {
                    center: Point2::new(s, 0.0),
                    radius: t,
                    ccw: true,
                },
                Segment::Line,
                Segment::Arc {
                    center: Point2::new(0.0, 0.0),
                    radius: t,
                    ccw: true,
                },
            ],
            tol,
        )
        .unwrap();
        Profile::new(outer, vec![], tol).unwrap()
    };
    let degenerate = |stage: &str, a: &Solid, b: &Solid| {
        for (op, r) in [
            ("fuse", a.fuse(OperationId(3), b)),
            ("cut", a.cut(OperationId(4), b)),
            ("common", a.common(OperationId(5), b)),
        ] {
            assert!(
                matches!(&r, Err(Error::Degenerate(m)) if m.contains("tangency")),
                "{stage} {op}: {:?}",
                r.map(|x| x.0.len())
            );
        }
    };
    let tilt = Frame3::new(
        Point3::new(1.0, -2.0, 0.5),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let h = 1.75;
    let (a, _) = Solid::extrude_with(OperationId(1), stadium(2.75, 1.0), tilt, 0.0, h).unwrap();
    let (stack, _) =
        Solid::extrude_with(OperationId(2), stadium(1.25, 2.25), tilt, h, h + 1.0).unwrap();
    let turned = Frame3::new(
        tilt.point(Point2::new(0.5, 0.25), h / 3.0),
        tilt.normal(),
        tilt.x() * 3.0 + tilt.y() * 4.0,
        tol,
    )
    .unwrap();
    let (ball, _) = Solid::sphere_with(OperationId(7), turned, 1.25, -half, half, tol).unwrap();
    degenerate("stadium", &a, &ball);
    // The fuzz input's chained stage: the stadium less the touching stack.
    let given = a.cut(OperationId(4), &stack).unwrap().0;
    assert_eq!(given.len(), 1);
    degenerate("given", &given[0], &ball);

    let up = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let square = Boundary::polygon(
        vec![
            Point2::new(-4.0, -4.0),
            Point2::new(4.0, -4.0),
            Point2::new(4.0, 4.0),
            Point2::new(-4.0, 4.0),
        ],
        tol,
    )
    .unwrap();
    let (block, _) = Solid::extrude_with(
        OperationId(1),
        Profile::new(square, vec![], tol).unwrap(),
        up,
        0.0,
        8.0,
    )
    .unwrap();
    let r = 2.0;
    let sphere = |e: f64| {
        let c = Frame3::new(
            Point3::new(4.0 - r + e, 0.0, 4.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            tol,
        )
        .unwrap();
        Solid::sphere_with(OperationId(2), c, r, -half, half, tol)
            .unwrap()
            .0
    };
    let tiny = 2f64.powi(-24);
    degenerate("box", &block, &sphere(tiny));
    degenerate("box gap", &block, &sphere(-tiny));
    let ball = 4.0 * std::f64::consts::PI * r * r * r / 3.0;
    for e in [2f64.powi(-16), -(2f64.powi(-16))] {
        // The common: the ball less the cap of height `e` past the wall.
        let k = e.max(0.0);
        let cap = std::f64::consts::PI * k * k * (3.0 * r - k) / 3.0;
        let common = block.common(OperationId(5), &sphere(e)).unwrap().0;
        let v: f64 = common.iter().map(|s| s.mass_properties().volume).sum();
        assert!((v - (ball - cap)).abs() <= 1e-9 * ball, "{e}: {v}");
    }
}

/// A ball within the resolution of tangency to a plane face is `Degenerate`
/// on either side of it (S9d.1's sections): missing the face's plane by a
/// gap no wider than the resolution it was a miss, and a ball resting on a
/// turned face but for rounding (an imported octahedron's face missed by
/// 4.6e-17) was fused with it into two solids, a ball inside a box cut
/// from it as a cavity behind a wall thinner than the resolution. Each
/// ball is placed by an exact search over its centre's last bits (the
/// frames' axes are rounded by each host's `hypot`): its distance from the
/// face's exact plane (the model's: a cap's through `o + h n` normal to
/// `x * y`, a wall's through its run's start normal to `d * n`) within a
/// factor of two of `r (1 + d)`, for `d` from `1e-17` to `1e-8` on either
/// side, the ball outside the body or inside it: on a turned box's top
/// (its foot inside the face and on its edge), a triangular prism's
/// slanted wall in a frame turned about two axes, that prism less a box (a
/// polyhedral result given), a level box's top (its faces' boxes, padded
/// far less than the resolution, apart by the gap) and two balls apart or
/// nested (their nearest points). Every operation is refused with the
/// existing rule's reason; `1e-6` of the radius on either side evaluates,
/// in the pair identities, and so do nearest points outside either face
/// (the foot past the top's edge, a zone whose sphere the top's plane
/// misses above its rims) and a gap inside both inputs (a slab's floor
/// beneath a dimple's sphere; a block below the dimple, the gap outside
/// it, is refused).
#[test]
fn a_ball_within_the_resolution_of_a_plane_face_is_degenerate() {
    use num_rational::BigRational as R;
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    type V = [R; 3];
    let tol = Tolerance::default();
    let q = |x: f64| R::from_float(x).unwrap();
    let v = |a: [f64; 3]| a.map(q);
    let add = |a: &V, b: &V| [&a[0] + &b[0], &a[1] + &b[1], &a[2] + &b[2]];
    let sub = |a: &V, b: &V| [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]];
    let scale = |a: &V, s: &R| [&a[0] * s, &a[1] * s, &a[2] * s];
    let dot = |a: &V, b: &V| &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2];
    let cross = |a: &V, b: &V| {
        [
            &a[1] * &b[2] - &a[2] * &b[1],
            &a[2] * &b[0] - &a[0] * &b[2],
            &a[0] * &b[1] - &a[1] * &b[0],
        ]
    };
    let frame = |o: [f64; 3], n: [f64; 3], x: [f64; 3]| {
        Frame3::new(
            Point3::new(o[0], o[1], o[2]),
            Vec3::new(n[0], n[1], n[2]),
            Vec3::new(x[0], x[1], x[2]),
            tol,
        )
        .unwrap()
    };
    let axes = |f: &Frame3| {
        [
            f.origin().to_array(),
            f.x().to_array(),
            f.y().to_array(),
            f.normal().to_array(),
        ]
        .map(v)
    };
    let prism = |f: Frame3, pts: &[(f64, f64)], hi: f64| {
        let b =
            Boundary::polygon(pts.iter().map(|&(x, y)| Point2::new(x, y)).collect(), tol).unwrap();
        Solid::extrude_with(
            OperationId(1),
            Profile::new(b, vec![], tol).unwrap(),
            f,
            0.0,
            hi,
        )
        .unwrap()
        .0
    };
    let ball = |c: [f64; 3], r: f64, op: u64| {
        let half = std::f64::consts::FRAC_PI_2;
        let f = frame(c, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        Solid::sphere_with(OperationId(op), f, r, -half, half, tol)
            .unwrap()
            .0
    };
    // A face's exact plane and its binary64 unit normal: a cap at `h`, the
    // wall over the run `a -> b`.
    let cap = |f: &Frame3, h: f64| {
        let a = axes(f);
        let p = add(&a[0], &scale(&a[3], &q(h)));
        (p, cross(&a[1], &a[2]), f.normal().to_array())
    };
    let wall = |f: &Frame3, a: (f64, f64), b: (f64, f64)| {
        let x = axes(f);
        let p = add(&x[0], &add(&scale(&x[1], &q(a.0)), &scale(&x[2], &q(a.1))));
        let d = add(&scale(&x[1], &q(b.0 - a.0)), &scale(&x[2], &q(b.1 - a.1)));
        let u = (f.x() * (b.0 - a.0) + f.y() * (b.1 - a.1))
            .cross(f.normal())
            .normalized()
            .unwrap();
        (p, cross(&d, &x[3]), u.to_array())
    };
    // `k` ulps up (away from zero for a negative value).
    let step = |x: f64, k: i64| f64::from_bits((x.to_bits() as i64 + k) as u64);
    // A centre near `c0` whose squared distance `dist2` lies within
    // `[(want + lo)^2, (want + hi)^2]` for the band about `gap`, over its
    // coordinates' last bits.
    let search = |c0: [f64; 3], dist2: &dyn Fn([f64; 3]) -> R, want: &R, gap: f64| {
        let (lo, hi) = if gap > 0.0 {
            (gap / 2.0, gap * 2.0)
        } else {
            (gap * 2.0, gap / 2.0)
        };
        let (a, b) = (want + q(lo), want + q(hi));
        let (a, b) = (&a * &a, &b * &b);
        for k in 0..=12i64 {
            for i in -k..=k {
                for j in -k..=k {
                    for l in -k..=k {
                        if i.abs().max(j.abs()).max(l.abs()) != k {
                            continue;
                        }
                        let c = [step(c0[0], i), step(c0[1], j), step(c0[2], l)];
                        let d2 = dist2(c);
                        if a <= d2 && d2 <= b {
                            return c;
                        }
                    }
                }
            }
        }
        panic!("no centre near {c0:?} for a gap of {gap:e}");
    };
    // A ball's centre at `foot + side r (1 + d) u`, its distance from the
    // plane `(p, m)` in the band about `r (1 + d)`.
    let place = |foot: [f64; 3], r: f64, d: f64, side: f64, face: &(V, V, [f64; 3])| {
        let (p, m, u) = face;
        let c0 = [0, 1, 2].map(|k| foot[k] + side * r * (1.0 + d) * u[k]);
        let mm = dot(m, m);
        let dist2 = |c: [f64; 3]| {
            let off = dot(m, &sub(&v(c), p));
            &off * &off / &mm
        };
        search(c0, &dist2, &q(r), d * r)
    };
    let refused = |what: &str, body: &Solid, b: &Solid| {
        for (op, out) in [
            ("fuse", body.fuse(OperationId(3), b)),
            ("cut", body.cut(OperationId(4), b)),
            ("cut back", b.cut(OperationId(4), body)),
            ("common", body.common(OperationId(5), b)),
        ] {
            assert!(
                matches!(&out, Err(Error::Degenerate(m))
                    if *m == "a plane crossing a sphere within the resolution of tangency (S9d.1)"),
                "{what} {op}: {:?}",
                out.map(|x| x.0.len())
            );
        }
    };
    let volume = |s: &[Solid]| s.iter().map(|x| x.mass_properties().volume).sum::<f64>();
    // Beyond the resolution every operation evaluates, the fuse `fused`
    // solids, in the pair identities.
    let evaluates = |what: &str, body: &Solid, b: &Solid, fused: usize| {
        let f = body.fuse(OperationId(3), b).unwrap().0;
        let c = body.cut(OperationId(4), b).unwrap().0;
        let m = body.common(OperationId(5), b).unwrap().0;
        assert_eq!(f.len(), fused, "{what}");
        let (va, vb) = (body.mass_properties().volume, b.mass_properties().volume);
        let near = |x: f64, y: f64| (x - y).abs() <= 1e-9 * (va + vb);
        assert!(near(volume(&f), va + vb - volume(&m)), "{what} fuse");
        assert!(near(volume(&c), va - volume(&m)), "{what} cut");
    };

    // A box turned off the world's axes (the boolean target's tilted
    // frame), its top: the foot inside the face, on its edge and past it.
    let tilt = frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let tilted = prism(
        tilt,
        &[(-3.0, -2.0), (3.0, -2.0), (3.0, 2.0), (-3.0, 2.0)],
        4.0,
    );
    let top = cap(&tilt, 4.0);
    let foot = tilt.point(Point2::new(0.5, 0.25), 4.0).to_array();
    for side in [1.0, -1.0] {
        for d in [1e-17, 1e-15, -1e-15, 1e-12, -1e-12, 1e-8] {
            let b = ball(place(foot, 1.5, d, side, &top), 1.5, 9);
            refused(&format!("turned box {side} {d:e}"), &tilted, &b);
        }
    }
    // On the top's edge, then past it (the ball 3.3e-5 from the edge: no
    // contact, evaluating as before).
    let edge = tilt.point(Point2::new(3.0, 0.25), 4.0).to_array();
    let b = ball(place(edge, 1.5, 1e-12, 1.0, &top), 1.5, 9);
    refused("turned box's edge", &tilted, &b);
    let past = tilt.point(Point2::new(3.01, 0.25), 4.0).to_array();
    let b = ball(place(past, 1.5, 1e-12, 1.0, &top), 1.5, 9);
    evaluates("past the turned box's edge", &tilted, &b, 2);
    // A zone of a sphere the top's plane misses within the resolution
    // above its rims (the boolean target's zone under a prism's top): no
    // contact, evaluating; a whole ball there is refused.
    let c = place(foot, 1.5, 1e-12, -1.0, &top);
    let about = Frame3::new(Point3::new(c[0], c[1], c[2]), tilt.normal(), tilt.x(), tol).unwrap();
    let (zone, _) = Solid::sphere_with(OperationId(9), about, 1.5, -0.5, 0.75, tol).unwrap();
    evaluates("a zone under the turned box's top", &tilted, &zone, 1);
    refused(
        "a ball under the turned box's top",
        &tilted,
        &ball(c, 1.5, 9),
    );
    for (d, side, fused) in [(1e-6, 1.0, 2), (-1e-6, 1.0, 1), (1e-6, -1.0, 1)] {
        let b = ball(place(foot, 1.5, d, side, &top), 1.5, 9);
        evaluates(&format!("turned box {side} {d:e}"), &tilted, &b, fused);
    }
    // The box less a ball through its top (a dimple, given) against a slab
    // across the top whose floor its sphere misses within the resolution:
    // the gap inside both inputs (their materials overlapping there, the
    // boolean target's and S9e.4b.4b.2's chained dimple), evaluating; a
    // block below the dimple with its top on that plane leaves the gap
    // outside the block, a wall thinner than the resolution: refused.
    let floor = Frame3::new(
        tilt.point(Point2::new(0.0, 0.0), 3.0),
        tilt.normal(),
        tilt.x(),
        tol,
    )
    .unwrap();
    let square = [(-8.0, -8.0), (8.0, -8.0), (8.0, 8.0), (-8.0, 8.0)];
    let slab = |lo: f64, hi: f64| {
        let b = Boundary::polygon(
            square.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
            tol,
        )
        .unwrap();
        let p = Profile::new(b, vec![], tol).unwrap();
        Solid::extrude_with(OperationId(21), p, floor, lo, hi)
            .unwrap()
            .0
    };
    let foot = floor.point(Point2::new(0.5, 0.25), 0.0).to_array();
    let c = place(foot, 1.5, 1e-12, 1.0, &cap(&floor, 0.0));
    let dimple = tilted.cut(OperationId(20), &ball(c, 1.5, 9)).unwrap().0;
    assert_eq!(dimple.len(), 1);
    evaluates("a dimple under a slab", &dimple[0], &slab(0.0, 2.0), 1);
    refused("a dimple over a block", &dimple[0], &slab(-2.0, 0.0));

    // A triangular prism in a frame turned about two axes, its slanted
    // wall; then less a box across its top (a polyhedral result given).
    let turn = frame([0.25, 0.5, -1.0], [1.0, 2.0, 2.0], [2.0, 1.0, -2.0]);
    let tri = prism(turn, &[(0.0, 0.0), (4.0, 0.0), (0.0, 3.0)], 5.0);
    let slant = wall(&turn, (4.0, 0.0), (0.0, 3.0));
    let foot = turn.point(Point2::new(2.0, 1.5), 2.5).to_array();
    let lid = Solid::box_at(
        Point3::new(-20.0, -20.0, foot[2] + 1.0),
        Vec3::new(40.0, 40.0, 40.0),
        tol,
    )
    .unwrap();
    let given = tri.cut(OperationId(2), &lid).unwrap().0;
    assert_eq!(given.len(), 1);
    for side in [1.0, -1.0] {
        for d in [1e-16, -1e-16, 1e-10] {
            let b = ball(place(foot, 0.75, d, side, &slant), 0.75, 9);
            refused(&format!("slanted wall {side} {d:e}"), &tri, &b);
            refused(&format!("given wall {side} {d:e}"), &given[0], &b);
        }
    }
    for (d, fused) in [(1e-6, 2), (-1e-6, 1)] {
        let b = ball(place(foot, 0.75, d, 1.0, &slant), 0.75, 9);
        evaluates(&format!("slanted wall {d:e}"), &tri, &b, fused);
        evaluates(&format!("given wall {d:e}"), &given[0], &b, fused);
    }

    // A level box's top: a gap wider than its faces' boxes' padding.
    let up = frame([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let level = prism(
        up,
        &[(-4.0, -4.0), (4.0, -4.0), (4.0, 4.0), (-4.0, 4.0)],
        8.0,
    );
    let top = cap(&up, 8.0);
    for side in [1.0, -1.0] {
        for d in [1e-8, -1e-8, 1e-12] {
            let b = ball(place([0.5, 0.25, 8.0], 1.5, d, side, &top), 1.5, 9);
            refused(&format!("level box {side} {d:e}"), &level, &b);
        }
    }
    let b = ball(place([0.5, 0.25, 8.0], 1.5, 1e-6, 1.0, &top), 1.5, 9);
    evaluates("level box", &level, &b, 2);

    // Two balls apart or nested within the resolution of tangency.
    let (c1, r1, r2) = ([0.25, -0.5, 0.125], 1.5, 0.75);
    let first = ball(c1, r1, 8);
    let u = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let dist2 = |c: [f64; 3]| {
        let w = sub(&v(c), &v(c1));
        dot(&w, &w)
    };
    for (what, side) in [("apart", 1.0), ("nested", -1.0)] {
        let want = q(r1 + side * r2);
        for d in [1e-12, -1e-12] {
            let c0 = [0, 1, 2].map(|k| c1[k] + (r1 + side * r2 * (1.0 + d)) * u[k]);
            let c = search(c0, &dist2, &want, side * d * r2);
            refused(&format!("balls {what} {d:e}"), &first, &ball(c, r2, 9));
        }
    }
}

/// A ball within the resolution of an edge or a vertex of the other input
/// is `Degenerate` on either side of it (S9d.1, as a plane face within it
/// is): missing a face's plane past its edge it was a miss, and a ball
/// resting on a box's edge or corner but for rounding was fused with it
/// into two solids, a ball inside an L's reflex edge cut from it behind a
/// wall thinner than the resolution; crossing it within the resolution
/// evaluated, was refused by the validator or as a piece thinner than the
/// resolution, by turns. Each ball's distance from the edge's exact line,
/// its exact circle or the vertex lies within a factor of two of `r (1 +
/// d)` (an exact search over its centre's last bits, the frames' axes
/// rounded by each host's `hypot`), for `d` from `1e-15` to `1e-9` on
/// either side: a turned box's top edge (the ball's foot on the top's plane
/// `1e-3 r` past the edge, and the ball on the faces' bisector) and corner,
/// an L's reflex edge from inside, a stadium's arc's rim (an exact circle;
/// in a turned frame, a conic within rounding of one), an imported
/// octahedron's edge and vertex and a hemisphere's rim circle. Every
/// operation is refused; `1e-6` of the radius evaluates in the pair
/// identities, and so does a gap inside both inputs (a cavity's sphere
/// within the resolution of the L's reflex edge inside its material).
#[test]
fn a_ball_within_the_resolution_of_an_edge_or_vertex_is_degenerate() {
    use num_rational::BigRational as R;
    use rusty_occt::identity::OperationId;
    use rusty_occt::{
        Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3,
    };
    type V = [R; 3];
    const EDGE: &str = "a sphere within the resolution of tangency to an edge (S9d.1)";
    const VERTEX: &str = "a sphere within the resolution of a vertex (S9d.1)";
    let tol = Tolerance::default();
    let half = std::f64::consts::FRAC_PI_2;
    let q = |x: f64| R::from_float(x).unwrap();
    let v = |a: [f64; 3]| a.map(q);
    let add = |a: &V, b: &V| [&a[0] + &b[0], &a[1] + &b[1], &a[2] + &b[2]];
    let sub = |a: &V, b: &V| [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]];
    let scale = |a: &V, s: &R| [&a[0] * s, &a[1] * s, &a[2] * s];
    let dot = |a: &V, b: &V| &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2];
    let unit = |a: [f64; 3]| {
        let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        a.map(|x| x / n)
    };
    let frame = |o: [f64; 3], n: [f64; 3], x: [f64; 3]| {
        Frame3::new(
            Point3::new(o[0], o[1], o[2]),
            Vec3::new(n[0], n[1], n[2]),
            Vec3::new(x[0], x[1], x[2]),
            tol,
        )
        .unwrap()
    };
    // A frame's exact point `(u, v, w)`, as its model takes it.
    let exact = |f: &Frame3, u: f64, w: f64, h: f64| {
        let a = [
            f.origin().to_array(),
            f.x().to_array(),
            f.y().to_array(),
            f.normal().to_array(),
        ]
        .map(v);
        let p = add(&a[0], &scale(&a[1], &q(u)));
        add(&p, &add(&scale(&a[2], &q(w)), &scale(&a[3], &q(h))))
    };
    let prism = |f: Frame3, pts: &[(f64, f64)], lo: f64, hi: f64, op: u64| {
        let b =
            Boundary::polygon(pts.iter().map(|&(x, y)| Point2::new(x, y)).collect(), tol).unwrap();
        Solid::extrude_with(
            OperationId(op),
            Profile::new(b, vec![], tol).unwrap(),
            f,
            lo,
            hi,
        )
        .unwrap()
        .0
    };
    let ball = |c: [f64; 3], r: f64| {
        let f = frame(c, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        Solid::sphere_with(OperationId(9), f, r, -half, half, tol)
            .unwrap()
            .0
    };
    // `k` ulps up (away from zero for a negative value).
    let step = |x: f64, k: i64| f64::from_bits((x.to_bits() as i64 + k) as u64);
    // A centre near `c0` (its coordinates `held` kept) whose squared
    // distance lies within `[(r + lo)^2, (r + hi)^2]` for the band about
    // the gap `r d`, over its coordinates' last bits.
    let search = |c0: [f64; 3], held: [bool; 3], dist2: &dyn Fn([f64; 3]) -> R, r: f64, d: f64| {
        let gap = r * d;
        let (lo, hi) = if gap > 0.0 {
            (gap / 2.0, gap * 2.0)
        } else {
            (gap * 2.0, gap / 2.0)
        };
        let (a, b) = (q(r) + q(lo), q(r) + q(hi));
        let (a, b) = (&a * &a, &b * &b);
        let span = |k: i64, i: usize| if held[i] { 0..=0 } else { -k..=k };
        for k in 0..=24i64 {
            for i in span(k, 0) {
                for j in span(k, 1) {
                    for l in span(k, 2) {
                        if i.abs().max(j.abs()).max(l.abs()) != k {
                            continue;
                        }
                        let c = [step(c0[0], i), step(c0[1], j), step(c0[2], l)];
                        let d2 = dist2(c);
                        if a <= d2 && d2 <= b {
                            return c;
                        }
                    }
                }
            }
        }
        panic!("no centre near {c0:?} for {d:e}");
    };
    let line = |p0: V, p1: V| {
        let e = sub(&p1, &p0);
        let ee = dot(&e, &e);
        move |c: [f64; 3]| {
            let w = sub(&v(c), &p0);
            let t = dot(&w, &e);
            dot(&w, &w) - &t * &t / &ee
        }
    };
    let point = |p: V| {
        move |c: [f64; 3]| {
            let w = sub(&v(c), &p);
            dot(&w, &w)
        }
    };
    // A ball of radius `r` about `at + r (1 + d) u`, placed.
    let near = |at: [f64; 3], u: [f64; 3], r: f64, d: f64, dist2: &dyn Fn([f64; 3]) -> R| {
        let c0 = [0, 1, 2].map(|k| at[k] + r * (1.0 + d) * u[k]);
        ball(search(c0, [false; 3], dist2, r, d), r)
    };
    let refused = |what: &str, body: &Solid, b: &Solid, why: &str| {
        for (op, out) in [
            ("fuse", body.fuse(OperationId(3), b)),
            ("cut", body.cut(OperationId(4), b)),
            ("cut back", b.cut(OperationId(4), body)),
            ("common", body.common(OperationId(5), b)),
        ] {
            assert!(
                matches!(&out, Err(Error::Degenerate(m)) if *m == why),
                "{what} {op}: {:?}",
                out.map(|x| x.0.len())
            );
        }
    };
    let volume = |s: &[Solid]| s.iter().map(|x| x.mass_properties().volume).sum::<f64>();
    let evaluates = |what: &str, body: &Solid, b: &Solid, fused: usize| {
        let f = body.fuse(OperationId(3), b).unwrap().0;
        let c = body.cut(OperationId(4), b).unwrap().0;
        let m = body.common(OperationId(5), b).unwrap().0;
        assert_eq!(f.len(), fused, "{what}");
        let (va, vb) = (body.mass_properties().volume, b.mass_properties().volume);
        let close = |x: f64, y: f64| (x - y).abs() <= 1e-9 * (va + vb);
        assert!(close(volume(&f), va + vb - volume(&m)), "{what} fuse");
        assert!(close(volume(&c), va - volume(&m)), "{what} cut");
    };
    let within = [1e-15, -1e-15, 1e-12, -1e-12, 1e-9, -1e-9];

    // A box turned off the world's axes: its top's edge at `u = 3` and its
    // corner at `(3, 2, 4)`.
    let tilt = frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let tilted = prism(
        tilt,
        &[(-3.0, -2.0), (3.0, -2.0), (3.0, 2.0), (-3.0, 2.0)],
        0.0,
        4.0,
        1,
    );
    let (n, x, y) = (tilt.normal(), tilt.x(), tilt.y());
    let mix = |a: f64, b: f64, c: f64| unit((n * a + x * b + y * c).to_array());
    let edge = line(exact(&tilt, 3.0, -2.0, 4.0), exact(&tilt, 3.0, 2.0, 4.0));
    let foot = tilt.point(Point2::new(3.0, 0.25), 4.0).to_array();
    for (dir, u) in [
        ("past", mix(1.0, 1e-3, 0.0)),
        ("bisector", mix(1.0, 1.0, 0.0)),
    ] {
        for d in within {
            let b = near(foot, u, 1.5, d, &edge);
            refused(&format!("turned box's edge {dir} {d:e}"), &tilted, &b, EDGE);
        }
    }
    for (d, fused) in [(1e-6, 2), (-1e-6, 1)] {
        let b = near(foot, mix(1.0, 1.0, 0.0), 1.5, d, &edge);
        evaluates(&format!("turned box's edge {d:e}"), &tilted, &b, fused);
    }
    let corner = point(exact(&tilt, 3.0, 2.0, 4.0));
    let at = tilt.point(Point2::new(3.0, 2.0), 4.0).to_array();
    for d in [1e-12, -1e-12, 1e-9] {
        let b = near(at, mix(1.0, 1.0, 1.0), 1.5, d, &corner);
        refused(&format!("turned box's corner {d:e}"), &tilted, &b, VERTEX);
    }
    let b = near(at, mix(1.0, 1.0, 1.0), 1.5, 1e-6, &corner);
    evaluates("turned box's corner", &tilted, &b, 2);

    // An L's reflex edge, a ball inside it (a wall thinner than the
    // resolution in the cut); a cavity's sphere there inside the L's
    // material (the gap inside both inputs) evaluating.
    let up = frame([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let ell = prism(
        up,
        &[
            (0.0, 0.0),
            (8.0, 0.0),
            (8.0, 4.0),
            (4.0, 4.0),
            (4.0, 8.0),
            (0.0, 8.0),
        ],
        0.0,
        8.0,
        1,
    );
    let reflex = line(exact(&up, 4.0, 4.0, 0.0), exact(&up, 4.0, 4.0, 8.0));
    let inward = unit([-1.0, -1.0, 0.0]);
    for d in [1e-12, -1e-12] {
        let b = near([4.0, 4.0, 4.0], inward, 1.5, d, &reflex);
        refused(&format!("L's reflex edge {d:e}"), &ell, &b, EDGE);
    }
    let b = near([4.0, 4.0, 4.0], inward, 1.5, 1e-6, &reflex);
    evaluates("L's reflex edge", &ell, &b, 1);
    let block = prism(
        up,
        &[(-2.0, -2.0), (10.0, -2.0), (10.0, 10.0), (-2.0, 10.0)],
        -2.0,
        10.0,
        7,
    );
    let b = near([4.0, 4.0, 4.0], inward, 1.5, 1e-12, &reflex);
    let hollow = block.cut(OperationId(2), &b).unwrap().0;
    assert_eq!(hollow.len(), 1);
    let fuse = ell.fuse(OperationId(3), &hollow[0]).unwrap().0;
    let common = ell.common(OperationId(5), &hollow[0]).unwrap().0;
    let back = hollow[0].cut(OperationId(4), &ell).unwrap().0;
    assert_eq!((fuse.len(), common.len(), back.len()), (1, 1, 1));
    let (vl, vh) = (
        ell.mass_properties().volume,
        hollow[0].mass_properties().volume,
    );
    assert!((volume(&fuse) - (vl + vh - volume(&common))).abs() <= 1e-9 * (vl + vh));
    assert!((volume(&back) - (vh - volume(&common))).abs() <= 1e-9 * (vl + vh));

    // A stadium's arc's rim about `(s, 0)` at its point `(s + t, 0, h)`:
    // the ball on its plane `y = 0`, its distance from the exact circle
    // `(|c_x - s| - t)^2 + (c_z - h)^2`.
    let (s, t, h) = (2.75, 1.0, 1.75);
    let stadium = |f: Frame3| {
        let outer = Boundary::path(
            vec![
                Point2::new(0.0, -t),
                Point2::new(s, -t),
                Point2::new(s, t),
                Point2::new(0.0, t),
            ],
            vec![
                Segment::Line,
                Segment::Arc {
                    center: Point2::new(s, 0.0),
                    radius: t,
                    ccw: true,
                },
                Segment::Line,
                Segment::Arc {
                    center: Point2::new(0.0, 0.0),
                    radius: t,
                    ccw: true,
                },
            ],
            tol,
        )
        .unwrap();
        Solid::extrude_with(
            OperationId(1),
            Profile::new(outer, vec![], tol).unwrap(),
            f,
            0.0,
            h,
        )
        .unwrap()
        .0
    };
    let level = stadium(up);
    let rim = |c: [f64; 3]| {
        let c = v(c);
        let (e, z) = (&c[0] - q(s + t), &c[2] - q(h));
        &e * &e + &z * &z
    };
    let r = 0.75;
    let rim_ball = |th: f64, d: f64| {
        let c0 = [
            s + t + r * (1.0 + d) * th.sin(),
            0.0,
            h + r * (1.0 + d) * th.cos(),
        ];
        ball(search(c0, [false, true, false], &rim, r, d), r)
    };
    for th in [1e-3, std::f64::consts::FRAC_PI_4] {
        for d in within {
            refused(
                &format!("stadium's rim {th} {d:e}"),
                &level,
                &rim_ball(th, d),
                EDGE,
            );
        }
    }
    evaluates(
        "stadium's rim",
        &level,
        &rim_ball(std::f64::consts::FRAC_PI_4, 1e-6),
        2,
    );
    // In the turned frame (a conic within rounding of a circle; gaps far
    // wider than its rounding, placed in binary64).
    let turned = stadium(tilt);
    for (d, fused) in [(1e-9, 0), (-1e-9, 0), (1e-6, 2)] {
        let k = r * (1.0 + d) * std::f64::consts::FRAC_1_SQRT_2;
        let c = (tilt.x() * (s + t + k) + tilt.normal() * (h + k)).to_array();
        let o = tilt.origin().to_array();
        let b = ball([0, 1, 2].map(|i| o[i] + c[i]), r);
        if fused == 0 {
            refused(&format!("turned stadium's rim {d:e}"), &turned, &b, EDGE);
        } else {
            evaluates("turned stadium's rim", &turned, &b, fused);
        }
    }

    // An imported octahedron (a polyhedron on its stored vertices): its edge
    // from `(5, 5, 8)` to `(9, 5, 4)` and that vertex.
    let (topology, resolution) = protocol::imported_topology("imported/octa.brep");
    let octa = Solid::imported_with(OperationId(1), topology, resolution)
        .unwrap()
        .0;
    let ridge = line(v([5.0, 5.0, 8.0]), v([9.0, 5.0, 4.0]));
    let apex = point(v([5.0, 5.0, 8.0]));
    let out = unit([1.0, 0.0, 1.0]);
    for d in [1e-12, -1e-12, 1e-9] {
        let b = near([7.0, 5.0, 6.0], out, 1.5, d, &ridge);
        refused(&format!("octahedron's edge {d:e}"), &octa, &b, EDGE);
        let b = near([5.0, 5.0, 8.0], unit([0.1, 0.2, 1.0]), 1.5, d, &apex);
        refused(&format!("octahedron's vertex {d:e}"), &octa, &b, VERTEX);
    }
    for (d, fused) in [(1e-6, 2), (-1e-6, 1)] {
        let b = near([7.0, 5.0, 6.0], out, 1.5, d, &ridge);
        evaluates(&format!("octahedron's edge {d:e}"), &octa, &b, fused);
    }

    // A hemisphere below `z = 0`: its rim circle at `(2, 0, 0)`.
    let (hemi, _) = Solid::sphere_with(OperationId(1), up, 2.0, -half, 0.0, tol).unwrap();
    let circle = |c: [f64; 3]| {
        let c = v(c);
        let e = &c[0] - q(2.0);
        &e * &e + &c[2] * &c[2]
    };
    let hemi_ball = |d: f64| {
        let k = r * (1.0 + d) * std::f64::consts::FRAC_1_SQRT_2;
        ball(
            search([2.0 + k, 0.0, k], [false, true, false], &circle, r, d),
            r,
        )
    };
    for d in [1e-12, -1e-12, 1e-9] {
        refused(
            &format!("hemisphere's rim {d:e}"),
            &hemi,
            &hemi_ball(d),
            EDGE,
        );
    }
    evaluates("hemisphere's rim", &hemi, &hemi_ball(1e-6), 2);
}
