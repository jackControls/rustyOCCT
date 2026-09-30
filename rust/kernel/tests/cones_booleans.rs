//! S9d.3b: Booleans of a cone or frustum against prisms with arcs, spheres
//! and cones, against the independent reference (`fixtures/boolean-cones-*`
//! from `tools/generate_cones_boolean_fixtures.py`).
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

/// The cases left to a later sub-step (`OutOfDomain`): none since S9d.3c
/// (`ball_tilt_common`, a loop in a turned frame, its height graph S9d.2c's
/// with the cone's radius).
const LATER: &[&str] = &[];

#[test]
fn every_case_matches_the_reference() {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-cones-expected.tsv")
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-cones-cases.txt")) {
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-cones-cases.txt")) {
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-cones-cases.txt")) {
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

/// A cone with its apex at its base (radius 0 at its frame's origin)
/// crossed by a rod: rings over the rod's angle or the cone's, the volumes'
/// identities and narrow enclosures, histories complete (the fuzz target's
/// replay: a `Curve3::Meet` on a cone of radius zero at its origin).
#[test]
fn a_cone_with_its_apex_at_its_base_carries_a_ring() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{history, Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let up = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let (cone, _) = Solid::cone_with(OperationId(1), up, 0.0, 2.0, 3.0, tol).unwrap();
    let side = Frame3::new(
        Point3::new(-4.0, 0.25, 2.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        tol,
    )
    .unwrap();
    let rod = Profile::new(
        Boundary::circle(Point2::new(0.0, 0.0), 0.5, tol).unwrap(),
        vec![],
        tol,
    )
    .unwrap();
    let (rod, _) = Solid::extrude_with(OperationId(2), rod, side, 0.0, 8.0).unwrap();
    let ins = [
        cone.topology().entity_set(cone.resolution()),
        rod.topology().entity_set(rod.resolution()),
    ];
    let volume = |out: &[Solid]| {
        out.iter()
            .map(|s| {
                let e = s.topology().mass_enclosure().unwrap();
                assert!(e.volume[1] - e.volume[0] <= 1e-9 * e.volume[1].abs().max(1.0));
                s.mass_properties().volume
            })
            .sum::<f64>()
    };
    let mut v = Vec::new();
    for r in [
        cone.fuse(OperationId(3), &rod),
        cone.cut(OperationId(3), &rod),
        cone.common(OperationId(3), &rod),
    ] {
        let (out, h) = r.unwrap();
        let outs: Vec<_> = out
            .iter()
            .map(|s| s.topology().entity_set(s.resolution()))
            .collect();
        assert!(history::check(&ins, &outs, &h).is_empty());
        v.push(volume(&out));
    }
    let (va, vb) = (cone.mass_properties().volume, rod.mass_properties().volume);
    assert!((v[0] - (va + vb - v[2])).abs() <= 1e-9 * v[0]);
    assert!((v[1] - (va - v[2])).abs() <= 1e-9 * va);
    assert!(v[2] > 0.0);
}

/// A sphere off a frustum's axis crossing its wall (`ball_side` with the
/// sphere's radius 0.875, clear of the end planes): a loop over the cone's
/// angle and height (S9d.3b.2); volumes, areas and centres the reference's
/// (`cones_boolean_reference.py` on that pair: fuse 16.8219782646248678,
/// cut 14.0158160766527348, common 0.64494964009963367), enclosures narrow,
/// histories complete.
#[test]
#[allow(clippy::excessive_precision)] // the reference's rows as printed
fn a_sphere_off_a_cones_axis_meets_it_in_a_loop() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{history, Frame3, Point3, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let up = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let (f, _) = Solid::cone_with(OperationId(1), up, 2.0, 1.0, 2.0, tol).unwrap();
    let at = Frame3::new(
        Point3::new(1.75, 0.5, 1.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let half = std::f64::consts::FRAC_PI_2;
    let (s, _) = Solid::sphere_with(OperationId(2), at, 0.875, -half, half, tol).unwrap();
    let ins = [
        f.topology().entity_set(f.resolution()),
        s.topology().entity_set(s.resolution()),
    ];
    // (volume, area, centre) of each operation's one solid.
    let want = [
        (
            16.821_978_264_624_868,
            41.647_480_310_556_612,
            [
                0.241_946_110_141_536_17,
                0.069_127_460_040_438_906,
                0.822_327_814_205_085_56,
            ],
        ),
        (
            14.015_816_076_652_735,
            37.532_610_200_928_205,
            [
                -0.059_987_347_032_922_15,
                -0.017_139_242_009_406_328,
                0.786_755_360_421_065_96,
            ],
        ),
        (
            0.644_949_640_099_633_67,
            4.756_054_652_133_274_7,
            [
                1.303_623_679_532_392_9,
                0.372_463_908_437_826_55,
                0.763_090_020_647_585_45,
            ],
        ),
    ];
    for (r, (v, a, c)) in [
        f.fuse(OperationId(3), &s),
        f.cut(OperationId(3), &s),
        f.common(OperationId(3), &s),
    ]
    .into_iter()
    .zip(want)
    {
        let (out, h) = r.unwrap();
        assert_eq!(out.len(), 1);
        let e = out[0].topology().mass_enclosure().unwrap();
        let slack = |x: f64| 1e-12 * x.abs().max(1.0);
        assert!(
            e.volume[0] - slack(v) <= v && v <= e.volume[1] + slack(v),
            "{:?} {v}",
            e.volume
        );
        assert!(e.surface_area[0] - slack(a) <= a && a <= e.surface_area[1] + slack(a));
        assert!(e.volume[1] - e.volume[0] <= 1e-9 * v);
        let m = out[0].mass_properties().centroid.to_array();
        assert!((0..3).all(|i| (m[i] - c[i]).abs() <= 1e-9), "{m:?} {c:?}");
        let outs = [out[0].topology().entity_set(out[0].resolution())];
        assert!(history::check(&ins, &outs, &h).is_empty());
    }
}
