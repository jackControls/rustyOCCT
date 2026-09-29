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
