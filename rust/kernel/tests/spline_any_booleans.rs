//! S9f.1: Booleans of spline prisms against polyhedral prisms in any
//! relative position, against the independent reference
//! (`fixtures/boolean-spline-any-*` from
//! `tools/generate_spline_any_boolean_fixtures.py`), and R4 under rounding:
//! a spline profile with an interior knot of multiplicity `p`, exactly C1,
//! extruded in a turned frame.
use rusty_occt::identity::OperationId;
use rusty_occt::topology::{Curve3, SplineSpan, Surface};
use rusty_occt::{
    BSplineCurve2, Boundary, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// `TILT` of the fixtures: the normal `(0, 3, 4) / 5`, `x` along the
/// world's `x` (stored bit for bit, `boolean-spline-any-frames.tsv`).
fn tilt() -> Frame3 {
    Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .expect("TILT")
}

/// R4's profile of the evidence: a quadratic from `(10, 2)` to `(0, 6)`
/// whose interior knot at `(5, 4)` has multiplicity two (the degree), C1
/// exactly (its pole the midpoint of `(7, 3)` and `(3, 5)`, equal spans),
/// closed by lines.
fn knot_profile() -> Profile {
    let points = vec![
        Point2::new(0.0, 0.0),
        Point2::new(10.0, 0.0),
        Point2::new(10.0, 2.0),
        Point2::new(0.0, 6.0),
    ];
    let poles = [(10.0, 2.0), (7.0, 3.0), (5.0, 4.0), (3.0, 5.0), (0.0, 6.0)]
        .iter()
        .map(|(x, y)| Point2::new(*x, *y))
        .collect();
    let curve = BSplineCurve2::new(2, poles, None, vec![0.0, 1.0, 2.0], vec![3, 2, 3])
        .expect("a quadratic");
    let segments = vec![
        Segment::Line,
        Segment::Line,
        Segment::Spline(SplineSpan::whole(curve)),
        Segment::Line,
    ];
    let outer = Boundary::path(points, segments, tol()).expect("a C1 spline path");
    Profile::new(outer, vec![], tol()).expect("a profile")
}

/// R4 under rounding: the knot's lifted pole lands 2^-53 off its
/// neighbours' midpoint in `TILT`, which refused the prism
/// (`InvalidTopology("edge_not_c1")`) before S9f.1. The lifting removes one
/// copy of the knot exactly (the profile is C1 there), so the walls and cap
/// edges keep the profile's poles less the knot's, C1 by construction, and
/// the prism validates with the profile's exact volume.
#[test]
fn a_knot_of_multiplicity_p_extrudes_in_a_turned_frame() {
    let profile = knot_profile();
    let area = profile.area();
    let (solid, _) =
        Solid::extrude_with(OperationId(1), profile, tilt(), 0.0, 5.0).expect("R4's prism in TILT");
    let t = solid.topology();
    let mut reduced = 0;
    for e in t.edges() {
        if let Curve3::BSpline(span) = &e.curve {
            assert_eq!(span.curve().multiplicities(), &[3, 1, 3]);
            reduced += 1;
        }
    }
    for f in t.faces() {
        if let Surface::BSpline(s) = &f.surface {
            assert_eq!(s.u_knots().multiplicities(), &[3, 1, 3]);
        }
    }
    assert_eq!(reduced, 2, "both cap edges");
    let m = t.mass_enclosure().expect("certified mass");
    let volume = area * 5.0;
    assert!(
        m.volume[0] <= volume * (1.0 + 1e-12) && volume * (1.0 - 1e-12) <= m.volume[1],
        "{:?} against {volume}",
        m.volume
    );
    // The same profile in the world's frame, its lift exact: the same
    // reduced curve, the profile's poles less the knot's.
    let (flat, _) = Solid::extrude_with(OperationId(2), knot_profile(), Frame3::xy(), 0.0, 5.0)
        .expect("R4's prism in XY");
    for e in flat.topology().edges() {
        if let Curve3::BSpline(span) = &e.curve {
            assert_eq!(span.curve().multiplicities(), &[3, 1, 3]);
            let poles: Vec<[f64; 2]> = span.curve().poles().iter().map(|p| [p.x, p.y]).collect();
            assert_eq!(
                poles,
                [[10.0, 2.0], [7.0, 3.0], [3.0, 5.0], [0.0, 6.0]],
                "the profile's poles less the knot's"
            );
        }
    }
}

#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

const CASES: &str = include_str!("../../fixtures/boolean-spline-any-cases.txt");
const EXPECTED: &str = include_str!("../../fixtures/boolean-spline-any-expected.tsv");
/// The fixtures added with the kernel (after the native capture): R4's
/// rounded knot in `TILT` and the blob in `TILT` against a `SIDE` box.
const R4_CASES: &str = include_str!("../../fixtures/boolean-spline-any-r4-cases.txt");
const R4_EXPECTED: &str = include_str!("../../fixtures/boolean-spline-any-r4-expected.tsv");

/// Every case of both sets.
fn all_cases() -> Vec<protocol::Case> {
    let mut out = protocol::cases(CASES);
    out.extend(protocol::cases(R4_CASES));
    out
}

/// Each case's kind and, for a result, its solids, volume, area and centre.
type Expect = (String, Option<(usize, [f64; 5])>);

fn expected(text: &str) -> std::collections::BTreeMap<String, Expect> {
    let mut expect = std::collections::BTreeMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let entry = expect
            .entry(name.to_string())
            .or_insert((String::new(), None));
        match w[0] {
            "expect" => entry.0 = w[1].to_string(),
            "result" => {
                let v: Vec<f64> = w[2..7].iter().map(|x| x.parse().unwrap()).collect();
                entry.1 = Some((
                    w[1].parse::<usize>().unwrap(),
                    [v[0], v[1], v[2], v[3], v[4]],
                ));
            }
            _ => {}
        }
    }
    expect
}

/// Every case of a fixture set as the reference declares it: degenerate
/// ones refused, empty ones empty, the others with the reference's solids
/// and their volume, area and centre within enclosures narrower than 1e-9
/// relative (a wide one would hold the reference yet report another).
fn check_set(cases: &str, expected_text: &str) {
    let expect = expected(expected_text);
    let mut failures = Vec::new();
    for case in protocol::cases(cases) {
        let (kind, want) = &expect[&case.name];
        let rows = match protocol::rows(&case) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{}: {e}", case.name));
                continue;
            }
        };
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
                    .take(10)
                    .map(|x| x.parse().unwrap_or(f64::NAN))
                    .collect()
            })
            .collect();
        if solids.len() != count || solids.iter().any(|s| s.len() < 10) {
            failures.push(format!("{}: {rows:?} for {count} solids", case.name));
            continue;
        }
        let sum = |i: usize| solids.iter().map(|s| s[i]).sum::<f64>();
        let near = |x: f64, lo: f64, hi: f64| {
            let slack = 1e-9 * x.abs().max(1.0);
            lo - slack <= x && x <= hi + slack
        };
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
                "{}: volume {} area {} against [{}, {}], [{}, {}]",
                case.name,
                v[0],
                v[1],
                sum(0),
                sum(1),
                sum(2),
                sum(3)
            ));
        }
        // The centre: the solids' first moments over the total volume.
        for i in 0..3 {
            let moment = |pick: fn(f64, f64) -> f64| {
                solids
                    .iter()
                    .map(|s| {
                        let (vl, vh, cl, ch) = (s[0], s[1], s[4 + 2 * i], s[5 + 2 * i]);
                        pick(pick(vl * cl, vl * ch), pick(vh * cl, vh * ch))
                    })
                    .sum::<f64>()
            };
            let (lo, hi) = (moment(f64::min), moment(f64::max));
            let m = v[0] * v[2 + i];
            let slack = 1e-9 * m.abs().max(v[0]).max(1.0);
            if !(lo - slack <= m && m <= hi + slack) || hi - lo > slack {
                failures.push(format!(
                    "{}: centre {i} {} against [{}, {}]",
                    case.name,
                    v[2 + i],
                    lo / v[0],
                    hi / v[0]
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures: {failures:#?}",
        failures.len()
    );
}

#[test]
fn every_case_matches_the_reference() {
    check_set(CASES, EXPECTED);
}

/// R4's rounded knot as a Boolean fixture (deferred by the evidence until
/// the extrusion took it), and generatrices of planes exactly parallel to a
/// turned wall's axis (`TILT` against `XY`'s and `SIDE`'s planes: `n . m`
/// exactly zero, where `TILTX`'s -8.9e-17 is `blob_rounding`'s refusal).
#[test]
fn the_rounded_knot_and_exact_parallels_match_the_reference() {
    check_set(R4_CASES, R4_EXPECTED);
}

#[test]
#[ignore]
fn debug_errors() {
    let want = std::env::var("CASE").unwrap_or_default();
    for case in all_cases() {
        if !case.name.contains(&want) {
            continue;
        }
        let t = std::time::Instant::now();
        match protocol::run(&case) {
            Ok((_, _, out, _)) => println!(
                "{} ({:?}): {} solids {:?}",
                case.name,
                t.elapsed(),
                out.len(),
                out.iter()
                    .map(|s| s.mass_properties().volume)
                    .collect::<Vec<_>>()
            ),
            Err(e) => println!("{} ({:?}): {e:?}", case.name, t.elapsed()),
        }
    }
}

#[test]
fn fixture_histories_are_complete() {
    for case in all_cases() {
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
    use rusty_occt::{Location, RigidTransform};
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    for case in all_cases() {
        let Ok((_, _, out, h)) = protocol::run(&case) else {
            continue;
        };
        let (_, _, again, h2) = protocol::run(&case).unwrap();
        assert_eq!(h, h2, "{}", case.name);
        for (x, y) in out.iter().zip(&again) {
            assert_eq!(ids(x), ids(y), "{}", case.name);
            assert_eq!(
                format!("{:?}", x.topology().vertices()),
                format!("{:?}", y.topology().vertices()),
                "{}",
                case.name
            );
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

fn spline_segment(poles: &[(f64, f64)], knots: Vec<f64>, mults: Vec<usize>) -> Segment {
    let degree = mults[0] - 1;
    let poles = poles.iter().map(|(x, y)| Point2::new(*x, *y)).collect();
    let curve = BSplineCurve2::new(degree, poles, None, knots, mults).expect("a spline");
    Segment::Spline(SplineSpan::whole(curve))
}

fn polygon(points: &[(f64, f64)]) -> Boundary {
    Boundary::polygon(
        points.iter().map(|(x, y)| Point2::new(*x, *y)).collect(),
        tol(),
    )
    .expect("a polygon")
}

/// The boolean target's replays of spline variants (S9f.1, its decode): a
/// box apart from a spline prism, the `+u` ray from the box's points
/// through the prism's joints at their height. The spline's crossings
/// count by the chords' predicate (an end above the point strictly), so a
/// joint at the point's height is not above it on either side; counted the
/// other way round (below) the parity flipped there and the box's faces
/// were kept inside the prism (an open result).
#[test]
fn a_ray_through_a_spline_joint_counts_as_the_chords_do() {
    let bulge = Boundary::path(
        vec![
            Point2::new(0.0, -0.5),
            Point2::new(1.0, -0.5),
            Point2::new(1.0, 0.5),
            Point2::new(0.0, 0.5),
        ],
        vec![
            Segment::Line,
            spline_segment(
                &[(1.0, -0.5), (1.5, 0.0), (1.0, 0.5)],
                vec![0.0, 1.0],
                vec![3, 3],
            ),
            Segment::Line,
            Segment::Line,
        ],
        tol(),
    )
    .expect("a bulge");
    let profile = Profile::new(bulge, vec![], tol()).expect("a profile");
    let (a, _) = Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 0.5).unwrap();
    let lean = Frame3::new(
        Point3::new(-2.0, 0.0, 0.25),
        Vec3::new(3.0, 0.0, 4.0),
        Vec3::new(0.0, 1.0, 0.0),
        tol(),
    )
    .unwrap();
    let rect = Profile::new(
        polygon(&[(0.0, 0.0), (1.0, 0.0), (1.0, 0.5), (0.0, 0.5)]),
        vec![],
        tol(),
    )
    .unwrap();
    let (b, _) = Solid::extrude_with(OperationId(2), rect, lean, 0.0, 0.5).unwrap();
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let (fused, _) = a.fuse(OperationId(3), &b).expect("apart: both");
    assert_eq!(fused.len(), 2);
    let total: f64 = fused.iter().map(|s| s.mass_properties().volume).sum();
    assert!((total - va - vb).abs() <= 1e-12 * (va + vb));
    let (cut, _) = a.cut(OperationId(4), &b).expect("apart: the object");
    assert_eq!(cut.len(), 1);
    assert!((cut[0].mass_properties().volume - va).abs() <= 1e-12 * va);
    let (common, _) = a.common(OperationId(5), &b).expect("apart: empty");
    assert!(common.is_empty());
}

/// The boolean target's replays of spline variants (S9f.1): a tilted box's
/// hole wall in the plane of a lens hole's corner (a joint of its two
/// cubics, the lens on one side of the plane) touches the prism along that
/// joint's vertical edge: a contact of the inputs, refused (`Degenerate`),
/// where it had been taken as the plane crossing there (a non-manifold
/// result).
#[test]
fn a_plane_touching_a_spline_joint_is_refused() {
    let (a, h) = (1.375, 1.375);
    let lower = [(-a, 0.0), (-a / 2.0, -h), (a / 2.0, -h), (a, 0.0)];
    let upper = [(a, 0.0), (a / 2.0, h), (-a / 2.0, h), (-a, 0.0)];
    let rev = |p: [(f64, f64); 4]| [p[3], p[2], p[1], p[0]];
    let cubic = |p: &[(f64, f64)]| spline_segment(p, vec![0.0, 1.0], vec![4, 4]);
    let hole = Boundary::path(
        vec![Point2::new(-a, 0.0), Point2::new(a, 0.0)],
        vec![cubic(&rev(upper)), cubic(&rev(lower))],
        tol(),
    )
    .expect("a lens");
    let s = 2.75;
    let outer = polygon(&[(-s, -s), (s, -s), (s, s), (-s, s)]);
    let profile = Profile::new(outer, vec![hole], tol()).unwrap();
    let (object, _) =
        Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 1.75).unwrap();
    let tilt = Frame3::new(
        Point3::new(-0.75, -1.25, 0.875),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    let (t, q) = (1.25, 0.625);
    let tool = Profile::new(
        polygon(&[(-t, -t), (t, -t), (t, t), (-t, t)]),
        vec![polygon(&[(-q, -q), (q, -q), (q, q), (-q, q)])],
        tol(),
    )
    .unwrap();
    let (tool, _) = Solid::extrude_with(OperationId(2), tool, tilt, 0.4375, 1.3125).unwrap();
    for r in [
        object.fuse(OperationId(3), &tool),
        object.cut(OperationId(4), &tool),
        object.common(OperationId(5), &tool),
    ] {
        assert!(
            matches!(r, Err(rusty_occt::Error::Degenerate(_))),
            "{:?}",
            r.map(|x| x.0.len())
        );
    }
}

/// R4 under rounding in one frame (S9a.2's profile Booleans, found by the
/// boolean target's corpus decoding R4's knot): a piece of a spline cut
/// across its knot of multiplicity `p` had its poles next to the knot
/// rounded off C1, and the result's profile was refused
/// (`InvalidCurve("a spline profile segment not C1 inside")`). Pieces are
/// cut from the curve with the knot removed once (`spline::piece`), C1 by
/// construction.
#[test]
fn a_piece_across_a_knot_of_multiplicity_p_stays_c1() {
    let (a, _) = Solid::extrude_with(OperationId(1), knot_profile(), Frame3::xy(), 0.0, 5.0)
        .expect("R4's prism");
    let square = Profile::new(
        polygon(&[(4.25, 3.0), (6.0, 3.0), (6.0, 7.0), (4.25, 7.0)]),
        vec![],
        tol(),
    )
    .unwrap();
    let (b, _) = Solid::extrude_with(OperationId(2), square, Frame3::xy(), 0.0, 5.0).unwrap();
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let volume = |r: Vec<Solid>| r.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let fused = volume(a.fuse(OperationId(3), &b).expect("a fuse").0);
    let common = volume(a.common(OperationId(5), &b).expect("a common").0);
    let cut = volume(a.cut(OperationId(4), &b).expect("a cut").0);
    assert!((fused - (va + vb - common)).abs() <= 1e-9 * fused);
    assert!((cut - (va - common)).abs() <= 1e-9 * va);
    assert!(common > 0.0 && common < vb);
}
