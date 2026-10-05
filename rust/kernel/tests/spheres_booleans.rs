//! S9d.2: Booleans of spheres against prisms with arcs and of two spheres,
//! against the independent reference (`fixtures/boolean-spheres-*` from
//! `tools/generate_spheres_boolean_fixtures.py`).
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

/// The cases left to a later sub-step (`OutOfDomain`): none since S9d.2c
/// took the loop in a turned frame (`bite_tilt_cut`).
const LATER: &[&str] = &[];

#[test]
fn every_case_matches_the_reference() {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-spheres-expected.tsv")
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-spheres-cases.txt")) {
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-spheres-cases.txt")) {
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-spheres-cases.txt")) {
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

/// A sphere centred on a cylinder's cap with the cylinder's radius touches
/// its wall all round the equator (the DRAW grids' `ZH5`): a tangency,
/// `Degenerate`, whatever the operation.
#[test]
fn a_sphere_tangent_to_a_cylinder_all_round_is_degenerate() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let outer = Boundary::circle(Point2::new(0.0, 0.0), 4.0, tol).unwrap();
    let profile = Profile::new(outer, vec![], tol).unwrap();
    let frame = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let (a, _) = Solid::extrude_with(OperationId(1), profile, frame, 0.0, 8.0).unwrap();
    let centre = Frame3::new(
        Point3::new(0.0, 0.0, 8.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let half = std::f64::consts::FRAC_PI_2;
    let (b, _) = Solid::sphere_with(OperationId(2), centre, 4.0, -half, half, tol).unwrap();
    for r in [
        a.fuse(OperationId(3), &b).map(|_| ()),
        a.cut(OperationId(3), &b).map(|_| ()),
        a.common(OperationId(3), &b).map(|_| ()),
    ] {
        assert!(matches!(r, Err(Error::Degenerate(_))), "{r:?}");
    }
}

/// The boolean target's `crash-43d93718` (the chained stage's sphere about
/// the middle of a holed prism's meeting with a sphere, `GIVEN_MET`), its
/// direct form: a ball whose axis is parallel to a round hole's and whose
/// centre lies on the hole's cylinder meets it in a loop through the
/// ball's north pole, inside the wall. Exactly on it, a section through a
/// sphere's pole off its meridians is refused (S9d.1); a centre within
/// rounding of it (the fuzz input's, an ulp off an exact one) passes too
/// near the pole for the pcurve's 256 anchors to pin its lift on the
/// sphere, `ComputationLimit` (it was `PrecisionLoss`, which a Boolean
/// does not document); a centre `0.0125` off evaluates.
#[test]
fn a_section_within_rounding_of_a_spheres_pole_is_a_computation_limit() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let s = 1.25;
    let outer = Boundary::polygon(
        [(-s, -s), (s, -s), (s, s), (-s, s)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec(),
        tol,
    )
    .unwrap();
    let hole = Boundary::circle(Point2::new(0.0, 0.0), 0.9375, tol).unwrap();
    let profile = Profile::new(outer, vec![hole], tol).unwrap();
    let (a, _) = Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 2.25).unwrap();
    let half = std::f64::consts::FRAC_PI_2;
    let ball = |x: f64, y: f64| {
        let f = Frame3::new(
            Point3::new(x, y, 0.5625),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(3.0, 4.0, 0.0),
            tol,
        )
        .unwrap();
        Solid::sphere_with(OperationId(2), f, 1.25, -half, half, tol)
            .unwrap()
            .0
    };
    let all = |b: &Solid| {
        [
            a.fuse(OperationId(3), b),
            a.cut(OperationId(4), b),
            a.common(OperationId(5), b),
        ]
    };
    // On the hole's cylinder exactly: 0.5625^2 + 0.75^2 = 0.9375^2.
    for r in all(&ball(0.5625, 0.75)) {
        assert!(
            matches!(&r, Err(Error::OutOfDomain(m)) if m.contains("pole")),
            "{:?}",
            r.map(|_| ())
        );
    }
    // Within rounding of it: an ulp off that point, and the fuzz input's
    // centre (the middle of a rounded meeting on the hole's wall).
    let ulp = f64::from_bits(0.75f64.to_bits() + 1);
    for (x, y) in [(0.5625, ulp), (-0.7184401491032446, -0.6022873086463865)] {
        for r in all(&ball(x, y)) {
            assert!(
                matches!(&r, Err(Error::ComputationLimit(m)) if m.contains("unpinned")),
                "{:?}",
                r.map(|_| ())
            );
        }
    }
    // 0.0125 off the cylinder (radius 0.95): fuse, cut and common evaluate.
    let b = ball(0.57, 0.76);
    let [f, c, m] = all(&b).map(|r| {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum::<f64>()
    });
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    assert!((f + m - va - vb).abs() <= 1e-9 * (va + vb), "{f} {m}");
    assert!((c + m - va).abs() <= 1e-9 * va, "{c} {m}");
}
