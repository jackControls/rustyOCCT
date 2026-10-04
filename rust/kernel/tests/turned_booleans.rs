//! S9c.2b.1: Booleans of prisms whose cylinders in turned frames meet in
//! quartics, against the independent reference (`fixtures/boolean-turned-*`
//! from `tools/generate_turned_boolean_fixtures.py`).
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

/// None left to a later sub-step.
const LATER: &[&str] = &[];

#[test]
fn every_case_matches_the_reference() {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-turned-expected.tsv")
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-turned-cases.txt")) {
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-turned-cases.txt")) {
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-turned-cases.txt")) {
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

/// A boolean fuzz variant: a square with a round hole in the tilted frame
/// fused with a torus band, its first solid (the holed prism) given to the
/// chained stage's turned cylinder, whose frame's normal is the tilted
/// frame's normalized again, an ulp off it. The cylinders' models crossed
/// within rounding of parallel, their meeting a piece of an ulp's sweep on
/// stored axes exactly parallel, which the validator could not evaluate (a
/// panic). Meeting, such axes are `Degenerate`; certainly apart within the
/// faces (one within the other, or beyond it), the pair evaluates.
#[test]
fn cylinders_within_rounding_of_parallel_are_degenerate_unless_apart() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let tilt = Frame3::new(
        Point3::new(1.0, -2.0, 0.5),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let disc = |r: f64| {
        Profile::new(
            Boundary::circle(Point2::new(0.0, 0.0), r, tol).unwrap(),
            vec![],
            tol,
        )
        .unwrap()
    };
    let s = 4.25;
    let square = Boundary::polygon(
        vec![
            Point2::new(-s, -s),
            Point2::new(s, -s),
            Point2::new(s, s),
            Point2::new(-s, s),
        ],
        tol,
    )
    .unwrap();
    let hole = Boundary::circle(Point2::new(0.0, 0.0), 1.375, tol).unwrap();
    let holed = Profile::new(square, vec![hole], tol).unwrap();
    let h = 2.25;
    let (prism, _) = Solid::extrude_with(OperationId(1), holed, tilt, 0.0, h).unwrap();
    let (rod, _) = Solid::extrude_with(OperationId(1), disc(1.375), tilt, 0.0, h).unwrap();
    let big = 0.75 * 2.75;
    let (band, _) = Solid::torus_with(
        OperationId(2),
        tilt,
        big,
        big * 0.5,
        -2.5,
        -0.25,
        std::f64::consts::TAU,
        tol,
    )
    .unwrap();
    let (fused, _) = prism.fuse(OperationId(3), &band).unwrap();
    // The chained stage's partner: the tilted frame's normal normalized
    // again (an ulp off), its x turned.
    let partner = |r: f64, at: Point2| {
        let f = Frame3::new(
            tilt.point(at, h / 3.0),
            tilt.normal(),
            tilt.x() * 3.0 + tilt.y() * 4.0,
            tol,
        )
        .unwrap();
        assert_ne!(f.normal(), tilt.normal(), "an ulp off");
        Solid::extrude_with(OperationId(7), disc(r), f, 0.0, h)
            .unwrap()
            .0
    };
    let meeting = partner(1.25, Point2::new(0.5, 0.25));
    for (stage, a) in [("given", &fused[0]), ("holed", &prism), ("rod", &rod)] {
        for (op, r) in [
            ("fuse", a.fuse(OperationId(9), &meeting)),
            ("cut", a.cut(OperationId(9), &meeting)),
            ("common", a.common(OperationId(9), &meeting)),
        ] {
            assert!(
                matches!(&r, Err(Error::Degenerate(m)) if m.contains("within rounding of parallel")),
                "{stage} {op}: {:?}",
                r.map(|x| x.0.len())
            );
        }
    }
    // Within the hole or the rod, holding the rod, and beyond them.
    let volume = |out: &[Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let near = |x: f64, y: f64| (x - y).abs() <= 1e-9 * x.abs().max(y.abs()).max(1.0);
    let third = h - h / 3.0;
    for (a, r, at, common) in [
        (&prism, 0.5, Point2::new(0.25, 0.0), 0.0),
        (&prism, 1.0, Point2::new(0.125, 0.0), 0.0),
        (&rod, 0.5, Point2::new(0.25, 0.0), 0.25 * third),
        (&rod, 2.0, Point2::new(0.25, 0.0), 1.890625 * third),
        (&rod, 0.25, Point2::new(3.0, 0.0), 0.0),
    ] {
        let b = partner(r, at);
        let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
        let (f, _) = a.fuse(OperationId(9), &b).unwrap();
        let (c, _) = a.cut(OperationId(9), &b).unwrap();
        let (m, _) = a.common(OperationId(9), &b).unwrap();
        let (f, c, m) = (volume(&f), volume(&c), volume(&m));
        assert!(near(f, va + vb - m), "{r}: fuse {f} for {va} + {vb} - {m}");
        assert!(near(c, va - m), "{r}: cut {c} for {va} - {m}");
        // The partner's part, or the rod's, over the heights h / 3..h.
        let want = std::f64::consts::PI * common;
        assert!(near(m, want), "{r}: common {m} for {want}");
    }
}
