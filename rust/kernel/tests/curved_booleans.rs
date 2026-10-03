//! S9c.1: Booleans of prisms with arcs in any relative position, against
//! the independent reference (`fixtures/boolean-curved-*` from
//! `tools/generate_curved_boolean_fixtures.py`).
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

/// None left to a later sub-step.
const LATER: &[&str] = &[];

/// Equal cylinders with meeting axes in stored turned frames: the models'
/// extents across the common perpendicular are equal exactly, a double
/// tangency (S9c.2b's decisions: `Degenerate`).
const NODES: &[&str] = &[
    "steinmetz_oblique_fuse",
    "steinmetz_oblique_common",
    "steinmetz_tilted_common",
];

#[test]
fn every_case_matches_the_reference() {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-curved-expected.tsv")
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-curved-cases.txt")) {
        let (kind, want) = &expect[&case.name];
        let rows = match protocol::rows(&case) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{}: {e}", case.name));
                continue;
            }
        };
        if NODES.contains(&case.name.as_str()) {
            if rows != ["refused"] {
                failures.push(format!("{}: {rows:?} not refused", case.name));
            }
            continue;
        }
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-curved-cases.txt")) {
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-curved-cases.txt")) {
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

/// Parallel circular cylinders (exact frames) where the second input's cap
/// circle crosses the first's wall (found by S9e.1): the second's arc meets
/// the first's cylinder where its own circle meets the first's, both in its
/// own frame (read in the first's frame, its crossings were missed and the
/// result left open). The common is the circles' lens over the heights both
/// hold, whichever input comes first.
#[test]
fn a_second_inputs_arc_across_a_parallel_cylinder() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let disc = |cx: f64, cy: f64, r: f64| {
        Profile::new(
            Boundary::circle(Point2::new(cx, cy), r, tol).unwrap(),
            vec![],
            tol,
        )
        .unwrap()
    };
    let (a, _) =
        Solid::extrude_with(OperationId(1), disc(0.0, 0.0, 3.0), Frame3::xy(), 0.0, 5.0).unwrap();
    let down = Frame3::new(
        Point3::new(0.0, 0.0, 6.0),
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    // Over z in [2, 6], its circle about (4, 0.5) of radius 2.
    let (b, _) = Solid::extrude_with(OperationId(2), disc(4.0, -0.5, 2.0), down, 0.0, 4.0).unwrap();
    let (r1, r2, d) = (3.0f64, 2.0f64, 16.25f64.sqrt());
    let xc = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let (a1, a2) = ((xc / r1).acos(), ((d - xc) / r2).acos());
    let lens = r1 * r1 * (a1 - a1.sin() * a1.cos()) + r2 * r2 * (a2 - a2.sin() * a2.cos());
    let volume = |out: Vec<Solid>| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    for common in [
        volume(a.common(OperationId(3), &b).unwrap().0),
        volume(b.common(OperationId(4), &a).unwrap().0),
    ] {
        assert!((common - 3.0 * lens).abs() <= 1e-9 * lens, "{common}");
    }
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let fuse = volume(a.fuse(OperationId(5), &b).unwrap().0);
    let cut = volume(a.cut(OperationId(6), &b).unwrap().0);
    assert!((fuse - (va + vb - 3.0 * lens)).abs() <= 1e-9 * fuse);
    assert!((cut - (va - 3.0 * lens)).abs() <= 1e-9 * cut);
}

/// A circle as two arcs (two faces on one cylinder, a periodic face split
/// at a seam) against a box crossing the arcs' joint, in frames with equal
/// axes and an offset that rounds (the curved engine): the same solids as
/// the circle's. A point on the arcs' shared chord, run either way, was
/// taken inside both circular segments and the box's cap piece inside the
/// circle kept, the fuse left open (S9e.4a found it; latent since S9c.1).
#[test]
fn two_arcs_of_one_circle_are_its_circle() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let at = |x: f64, z: f64| {
        Frame3::new(
            Point3::new(x, 0.0, z),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            tol,
        )
        .unwrap()
    };
    let arc = Segment::Arc {
        center: Point2::new(5.0, 5.0),
        radius: 2.5,
        ccw: true,
    };
    let halves = Boundary::path(
        vec![Point2::new(7.5, 5.0), Point2::new(2.5, 5.0)],
        vec![arc.clone(), arc],
        tol,
    )
    .unwrap();
    let circle = Boundary::circle(Point2::new(5.0, 5.0), 2.5, tol).unwrap();
    let square = Boundary::polygon(
        vec![
            Point2::new(5.5, 2.0),
            Point2::new(9.0, 2.0),
            Point2::new(9.0, 8.0),
            Point2::new(5.5, 8.0),
        ],
        tol,
    )
    .unwrap();
    let prism = |b: Boundary, op: u64, f: Frame3, h: [f64; 2]| {
        let profile = Profile::new(b, vec![], tol).unwrap();
        Solid::extrude_with(OperationId(op), profile, f, h[0], h[1])
            .unwrap()
            .0
    };
    let tool = prism(square, 2, at(0.0, 1.0), [0.0, 3.0]);
    let volume = |out: Vec<Solid>| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    for (k, profile) in [halves, circle].into_iter().enumerate() {
        let object = prism(profile, 1, at(0.1, 0.0), [0.0, 5.0]);
        let f = volume(object.fuse(OperationId(3), &tool).unwrap().0);
        let c = volume(object.cut(OperationId(4), &tool).unwrap().0);
        let m = volume(object.common(OperationId(5), &tool).unwrap().0);
        let (va, vb) = (
            object.mass_properties().volume,
            tool.mass_properties().volume,
        );
        assert!((f - (va + vb - m)).abs() <= 1e-9 * f, "{k}: fuse {f}");
        assert!((c - (va - m)).abs() <= 1e-9 * va, "{k}: cut {c}");
        assert!(
            (m - 23.478130341535252).abs() <= 1e-9 * m,
            "{k}: common {m}"
        );
    }
}
