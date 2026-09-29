//! S9d.4b.1: Booleans of a torus v-segment or wedge against polyhedral
//! prisms, against the independent reference
//! (`fixtures/boolean-torus-segment-*` from
//! `tools/generate_torus_segment_boolean_fixtures.py`).
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

/// None left to a later sub-step.
const LATER: &[&str] = &[];

#[test]
fn every_case_matches_the_reference() {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-torus-segment-expected.tsv")
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
    for case in protocol::cases(include_str!(
        "../../fixtures/boolean-torus-segment-cases.txt"
    )) {
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
    for case in protocol::cases(include_str!(
        "../../fixtures/boolean-torus-segment-cases.txt"
    )) {
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
    for case in protocol::cases(include_str!(
        "../../fixtures/boolean-torus-segment-cases.txt"
    )) {
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

#[test]
fn segments_and_wedges_classify_exactly() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Frame3, Location, Point3, Solid, Tolerance};
    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    let part = |low: f64, high: f64, angle: f64| {
        let tol = Tolerance::default();
        Solid::torus_with(
            OperationId(1),
            Frame3::xy(),
            2.5,
            1.5,
            low,
            high,
            angle,
            tol,
        )
        .unwrap()
        .0
    };
    let (i, b, o) = (Location::Inside, Location::Boundary, Location::Outside);
    type Points<'a> = &'a [([f64; 3], Location)];
    let cases: [(Solid, Points); 5] = [
        // The outer half: a barrel from the axis to the tube's outside.
        (
            part(-FRAC_PI_2, FRAC_PI_2, TAU),
            &[
                ([0.0, 0.0, 0.0], i),
                ([3.5, 0.0, 0.0], i),
                ([1.0, 0.0, 1.4999], i),
                ([4.5, 0.0, 0.0], o),
                ([0.0, 0.0, 1.6], o),
                ([2.5, 0.0, 1.5], b),
                ([1.0, 1.0, -1.5], b),
                ([0.0, 4.0, 0.0], b),
            ],
        ),
        // The inner half: a spool, its waist of radius 1.
        (
            part(FRAC_PI_2, 3.0 * FRAC_PI_2, TAU),
            &[
                ([0.0, 0.0, 0.0], i),
                ([1.5, 0.0, 0.0], o),
                ([0.0, 1.0, 0.0], b),
                ([1.9, 0.0, 1.4], i),
                ([2.0, 0.0, 1.6], o),
            ],
        ),
        // An upper band: the tube's cap above its high end's disc.
        (
            part(0.5, 2.25, TAU),
            &[
                ([2.5, 0.0, 1.4], i),
                ([0.0, 0.0, 1.0], i),
                ([0.0, 0.0, 1.3], o),
                ([0.0, 0.0, 0.5], o),
                ([3.0, 0.0, 0.2], o),
            ],
        ),
        // A quarter turn of the whole tube.
        (
            part(0.0, TAU, FRAC_PI_2),
            &[
                ([2.5, 0.5, 0.0], i),
                ([0.5, 2.5, 0.0], i),
                ([2.5, -0.5, 0.0], o),
                ([-0.5, 2.5, 0.0], o),
                ([2.5, 0.0, 0.0], b),
                ([2.5, 2.5, 0.0], i),
            ],
        ),
        // Three quarters: two half-planes' union.
        (
            part(0.0, TAU, 1.5 * PI),
            &[
                ([-2.5, 0.5, 0.0], i),
                ([-0.5, -2.5, 0.0], i),
                ([2.5, -0.5, 0.0], o),
                ([0.0, -2.5, 0.0], b),
            ],
        ),
    ];
    for (k, (s, points)) in cases.iter().enumerate() {
        for (p, want) in points.iter() {
            assert_eq!(
                s.classify(Point3::new(p[0], p[1], p[2])).unwrap(),
                *want,
                "part {k} at {p:?}"
            );
        }
    }
}

/// The fuzz replay's three-quarter wedge (`fuzz/regressions/README.md`):
/// coaxial with a tilted square frame with a square hole, its fused wall a
/// patch over more than half a turn whose holes' pcurves, lifted from the
/// surface's principal angles, lay a turn away from its outer loop.
#[test]
fn a_three_quarter_wedge_fused_holds_its_holes() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let frame = Frame3::new(
        Point3::new(1.0, -2.0, 0.5),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let square = |s: f64| {
        let pts = [(-s, -s), (s, -s), (s, s), (-s, s)];
        Boundary::polygon(pts.map(|(x, y)| Point2::new(x, y)).to_vec(), tol).unwrap()
    };
    let profile = Profile::new(square(3.25), vec![square(1.625)], tol).unwrap();
    let (a, _) = Solid::extrude_with(OperationId(1), profile, frame, 0.0, 0.75).unwrap();
    let turn = std::f64::consts::TAU;
    let (b, _) = Solid::torus_with(
        OperationId(2),
        frame,
        3.1875,
        0.796875,
        0.0,
        turn,
        3.0 * std::f64::consts::FRAC_PI_2,
        tol,
    )
    .unwrap();
    let volume = |out: &[Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let fused = volume(&a.fuse(OperationId(3), &b).unwrap().0);
    let common = volume(&a.common(OperationId(4), &b).unwrap().0);
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    assert!(
        (fused - (va + vb - common)).abs() <= 1e-9 * fused,
        "{fused} {va} {vb} {common}"
    );
}

/// A cavity in a segment's wall region: the validator's rays against a
/// torus face with loops are not decided, so the cut is `ComputationLimit`
/// (not an invalid body).
#[test]
fn a_cavity_in_a_segment_is_undecided() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Error, Frame3, Point2, Profile, Solid, Tolerance};
    use std::f64::consts::{FRAC_PI_2, TAU};
    let tol = Tolerance::default();
    let (half, _) = Solid::torus_with(
        OperationId(1),
        Frame3::xy(),
        2.5,
        1.5,
        -FRAC_PI_2,
        FRAC_PI_2,
        TAU,
        tol,
    )
    .unwrap();
    let pts = [(3.0, -0.25), (3.5, -0.25), (3.5, 0.25), (3.0, 0.25)];
    let square = Boundary::polygon(pts.map(|(x, y)| Point2::new(x, y)).to_vec(), tol).unwrap();
    let profile = Profile::new(square, vec![], tol).unwrap();
    let (block, _) =
        Solid::extrude_with(OperationId(2), profile, Frame3::xy(), -0.25, 0.25).unwrap();
    assert!(matches!(
        half.cut(OperationId(3), &block),
        Err(Error::ComputationLimit(_))
    ));
}
