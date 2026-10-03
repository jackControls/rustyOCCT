//! S9f.2b.1: Booleans of spline prisms against prisms with arc or circle
//! walls whose axes cross, against the independent reference
//! (`fixtures/boolean-spline-crossing-*` from
//! `tools/generate_spline_crossing_boolean_fixtures.py`): every meeting a
//! graph over the spline's parameter inside the faces as the reference,
//! S9f.2b.2's loops refused (`OutOfDomain`), tangencies and turning points
//! at knots `Degenerate`; and the decisions' other refusals: a cylinder's
//! cap circle on a spline wall in a plane holding the wall's axis (a tower
//! field), spline walls against spline walls on crossing axes.
use rusty_occt::identity::OperationId;
use rusty_occt::topology::SplineSpan;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance,
    Vec3,
};

#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

const CASES: &str = include_str!("../../fixtures/boolean-spline-crossing-cases.txt");
const EXPECTED: &str = include_str!("../../fixtures/boolean-spline-crossing-expected.tsv");

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Each case's kind, its sub-step and, for a result, its solids, volume,
/// area and centre.
type Expect = (String, String, Option<(usize, [f64; 5])>);

fn expected(text: &str) -> std::collections::BTreeMap<String, Expect> {
    let mut expect = std::collections::BTreeMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let entry = expect
            .entry(name.to_string())
            .or_insert((String::new(), String::new(), None));
        match w[0] {
            "expect" => {
                entry.0 = w[1].to_string();
                entry.1 = w[2].to_string();
            }
            "result" => {
                let v: Vec<f64> = w[2..7].iter().map(|x| x.parse().unwrap()).collect();
                entry.2 = Some((
                    w[1].parse::<usize>().unwrap(),
                    [v[0], v[1], v[2], v[3], v[4]],
                ));
            }
            _ => {}
        }
    }
    expect
}

/// Every case as the reference declares it: degenerate ones refused,
/// S9f.2b.2's loops out of the domain, the others with the reference's
/// solids and their volume, area and centre within enclosures narrower
/// than 1e-9 relative (a wide one would hold the reference yet report
/// another).
#[test]
fn every_case_matches_the_reference() {
    let expect = expected(EXPECTED);
    let mut failures = Vec::new();
    for case in protocol::cases(CASES) {
        let (kind, step, want) = &expect[&case.name];
        if step == "S9f.2b.2" {
            match protocol::run(&case) {
                Err(Error::OutOfDomain(m)) if m.contains("S9f.2b.2") => {}
                r => failures.push(format!(
                    "{}: {:?} not S9f.2b.2's",
                    case.name,
                    r.map(|x| x.2.len())
                )),
            }
            continue;
        }
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
#[ignore]
fn debug_errors() {
    let want = std::env::var("CASE").unwrap_or_default();
    for case in protocol::cases(CASES) {
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
    for case in protocol::cases(CASES) {
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
    for case in protocol::cases(CASES) {
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

// ------------------------------------------------------------ the decisions' refusals

fn spline_segment(poles: &[(f64, f64)], knots: Vec<f64>, mults: Vec<usize>) -> Segment {
    let degree = mults[0] - 1;
    let poles = poles.iter().map(|(x, y)| Point2::new(*x, *y)).collect();
    let curve = BSplineCurve2::new(degree, poles, None, knots, mults).expect("a spline");
    Segment::Spline(SplineSpan::whole(curve))
}

/// The dome `y = x (4 - x) / 2` over `[0, 4]`, closed by its base.
fn dome() -> Profile {
    let outer = Boundary::path(
        vec![Point2::new(0.0, 0.0), Point2::new(4.0, 0.0)],
        vec![
            Segment::Line,
            spline_segment(
                &[(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)],
                vec![0.0, 1.0],
                vec![3, 3],
            ),
        ],
        tol(),
    )
    .expect("the dome");
    Profile::new(outer, vec![], tol()).expect("a profile")
}

fn frame(origin: (f64, f64, f64), normal: (f64, f64, f64), x: (f64, f64, f64)) -> Frame3 {
    Frame3::new(
        Point3::new(origin.0, origin.1, origin.2),
        Vec3::new(normal.0, normal.1, normal.2),
        Vec3::new(x.0, x.1, x.2),
        tol(),
    )
    .expect("a frame")
}

fn prism(id: u64, profile: Profile, f: Frame3, lo: f64, hi: f64) -> Solid {
    Solid::extrude_with(OperationId(id), profile, f, lo, hi)
        .expect("a prism")
        .0
}

fn disc(cx: f64, cy: f64, r: f64) -> Profile {
    Profile::new(
        Boundary::circle(Point2::new(cx, cy), r, tol()).expect("a circle"),
        vec![],
        tol(),
    )
    .expect("a disc")
}

/// A rod along `x` (its axis exactly perpendicular to the dome's) ending at
/// `x = 2` inside the dome's prism: its cap's plane holds the dome's axis,
/// meeting the wall in the apex's generatrix, where the circle's points lie
/// in a tower field (S9f.2b.2's); the same rod through the whole prism, its
/// caps outside, is S9f.2b.1's (`dome_side`).
#[test]
fn a_cap_circle_on_a_spline_wall_along_its_axis_is_s9f2b2s() {
    let a = prism(1, dome(), Frame3::xy(), 0.0, 2.0);
    let side = frame((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0));
    // The rod about (y, z) = (1, -1) of radius 2 to x = 2: its cap circle
    // meets the arch's generatrix at x = 2 (y = 2) at z = sqrt(3) - 1.
    let short = prism(2, disc(1.0, -1.0, 2.0), side, 0.0, 3.0);
    let r = a.common(OperationId(3), &short);
    assert!(
        matches!(&r, Err(Error::OutOfDomain(m)) if m.contains("S9f.2b.2")),
        "{:?}",
        r.map(|x| x.0.len())
    );
    // Through the whole prism: its upper branch over the arch from end to
    // end.
    let long = prism(4, disc(1.0, -1.0, 2.0), side, 0.0, 6.0);
    let (va, vb) = (a.mass_properties().volume, long.mass_properties().volume);
    let volume = |r: Vec<Solid>| r.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let fused = volume(a.fuse(OperationId(5), &long).expect("a fuse").0);
    let common = volume(a.common(OperationId(6), &long).expect("a common").0);
    let cut = volume(a.cut(OperationId(7), &long).expect("a cut").0);
    assert!((fused - (va + vb - common)).abs() <= 1e-9 * fused);
    assert!((cut - (va - common)).abs() <= 1e-9 * va);
    assert!(common > 0.0 && common < va.min(vb));
}

/// Spline walls against spline walls on crossing axes stay refused (by
/// design, "S9f refined").
#[test]
fn crossing_spline_walls_stay_refused() {
    let a = prism(1, dome(), Frame3::xy(), 0.0, 5.0);
    let side = frame((-1.0, 1.0, 2.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0));
    let c = prism(4, dome(), side, 0.0, 6.0);
    let r = a.fuse(OperationId(5), &c);
    assert!(
        matches!(&r, Err(Error::OutOfDomain(m)) if m.contains("crossing axes")),
        "{:?}",
        r.map(|x| x.0.len())
    );
}

/// The square `[-s, s]^2`.
fn square(s: f64) -> Boundary {
    Boundary::polygon(
        vec![
            Point2::new(-s, -s),
            Point2::new(s, -s),
            Point2::new(s, s),
            Point2::new(-s, s),
        ],
        tol(),
    )
    .expect("a square")
}

/// Three operations' inclusion and exclusion, each result validated as it
/// is built.
fn agree(a: &Solid, b: &Solid) {
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let volume = |r: Vec<Solid>| r.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let fused = volume(a.fuse(OperationId(5), b).expect("a fuse").0);
    let common = volume(a.common(OperationId(6), b).expect("a common").0);
    let cut = volume(a.cut(OperationId(7), b).expect("a cut").0);
    assert!(
        (fused - (va + vb - common)).abs() <= 1e-9 * fused,
        "{fused} {va} {vb} {common}"
    );
    assert!(
        (cut - (va - common)).abs() <= 1e-9 * va,
        "{cut} {va} {common}"
    );
}

/// A slab with the fuzz target's lens hole (two cubics, clockwise) and a
/// leaning slab with a round hole: the hole's cylinder meets the lens's
/// lower wall in a curve ending on the top cap `6.7e-7` in the wall's `u`
/// short of a turning point `0.0049` above the cap, where its `v` is a
/// square root of `d` near zero. Its ends' enclosures (`d` from its exact
/// polynomial, not from its enclosed factors: 3.5e-12 wide, not 2.2e-10),
/// the closing chords past the wall's domain and the integrals' halvings
/// toward the turning point (some twenty) certify it (a boolean fuzz
/// variant: `uncertified_shell_orientation` before).
#[test]
fn a_meeting_ending_near_a_turning_point_validates() {
    let lens = Boundary::path(
        vec![Point2::new(-1.5, 0.0), Point2::new(1.5, 0.0)],
        vec![
            spline_segment(
                &[(-1.5, 0.0), (-0.75, 1.5), (0.75, 1.5), (1.5, 0.0)],
                vec![0.0, 1.0],
                vec![4, 4],
            ),
            spline_segment(
                &[(1.5, 0.0), (0.75, -1.5), (-0.75, -1.5), (-1.5, 0.0)],
                vec![0.0, 1.0],
                vec![4, 4],
            ),
        ],
        tol(),
    )
    .expect("the lens");
    let object = Profile::new(square(3.0), vec![lens], tol()).expect("a profile");
    let hole = Boundary::circle(Point2::new(0.0, 0.0), 1.5, tol()).expect("a circle");
    let tool = Profile::new(square(3.0), vec![hole], tol()).expect("a profile");
    let a = prism(1, object, Frame3::xy(), 0.0, 0.5);
    let lean = frame((-1.625, -1.625, 0.25), (3.0, 0.0, 4.0), (0.0, 1.0, 0.0));
    let b = prism(2, tool, lean, 0.0, 0.5);
    agree(&a, &b);
}

/// R4's knot profile (a quadratic of two spans from `(2, 0.5)` to `(0, 1)`)
/// against a leaning stadium: the ruling at `u = 1/2`, `(1.625, 0.75)`, is
/// tangent to the stadium's arc's cylinder on the object's top cap but for
/// the frame's rounding (its normal `(3, 0, 4) / 5`), its turning point a
/// rounding away from the cap's edge: `Degenerate` (a boolean fuzz variant:
/// its meeting's end left its loop winding undecided).
#[test]
fn a_turning_point_within_rounding_of_a_cap_is_degenerate() {
    let outer = Boundary::path(
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 0.5),
            Point2::new(0.0, 1.0),
        ],
        vec![
            Segment::Line,
            Segment::Line,
            spline_segment(
                &[
                    (2.0, 0.5),
                    (1.75, 0.75),
                    (1.0, 1.0),
                    (0.25, 1.25),
                    (0.0, 1.0),
                ],
                vec![0.0, 1.0, 2.0],
                vec![3, 2, 3],
            ),
            Segment::Line,
        ],
        tol(),
    )
    .expect("the knot profile");
    let knot = Profile::new(outer, vec![], tol()).expect("a profile");
    let (s, t) = (2.0, 0.625);
    let stadium = Boundary::path(
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
        tol(),
    )
    .expect("a stadium");
    let stadium = Profile::new(stadium, vec![], tol()).expect("a profile");
    let a = prism(1, knot, Frame3::xy(), 0.0, 2.0);
    let lean = frame((0.875, 1.375, 1.0), (3.0, 0.0, 4.0), (0.0, 1.0, 0.0));
    let b = prism(2, stadium, lean, 0.5, 1.5);
    for r in [
        a.fuse(OperationId(3), &b),
        a.cut(OperationId(4), &b),
        a.common(OperationId(5), &b),
    ] {
        assert!(
            matches!(&r, Err(Error::Degenerate(m)) if m.contains("turning back on a face's boundary")),
            "{:?}",
            r.map(|x| x.0.len())
        );
    }
}
