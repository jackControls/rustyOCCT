//! S9f.2a: Booleans of spline prisms against prisms with arc, circle or
//! spline walls whose axes are exactly parallel, against the independent
//! reference (`fixtures/boolean-spline-parallel-*` from
//! `tools/generate_spline_parallel_boolean_fixtures.py`), and the decisions'
//! refusals: coincident spline walls, crossing axes, arcs whose degrees'
//! product exceeds 16.
use rusty_occt::identity::OperationId;
use rusty_occt::topology::SplineSpan;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance,
    Vec3,
};

#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

const CASES: &str = include_str!("../../fixtures/boolean-spline-parallel-cases.txt");
const EXPECTED: &str = include_str!("../../fixtures/boolean-spline-parallel-expected.tsv");

fn tol() -> Tolerance {
    Tolerance::default()
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

/// Every case as the reference declares it: degenerate ones refused, the
/// others with the reference's solids and their volume, area and centre
/// within enclosures narrower than 1e-9 relative (a wide one would hold the
/// reference yet report another).
#[test]
fn every_case_matches_the_reference() {
    let expect = expected(EXPECTED);
    let mut failures = Vec::new();
    for case in protocol::cases(CASES) {
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

/// A region under one Bezier arch from `(w, 0)` back to the origin through
/// `poles` (closed by the base line).
fn arch(w: f64, poles: &[(f64, f64)]) -> Profile {
    let degree = poles.len() - 1;
    let outer = Boundary::path(
        vec![Point2::new(0.0, 0.0), Point2::new(w, 0.0)],
        vec![
            Segment::Line,
            spline_segment(poles, vec![0.0, 1.0], vec![degree + 1, degree + 1]),
        ],
        tol(),
    )
    .expect("an arch");
    Profile::new(outer, vec![], tol()).expect("a profile")
}

/// The dome `y = x (4 - x) / 2` over `[0, 4]`.
fn dome() -> Profile {
    arch(4.0, &[(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)])
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

fn disc(r: f64) -> Profile {
    Profile::new(
        Boundary::circle(Point2::new(0.0, 0.0), r, tol()).expect("a circle"),
        vec![],
        tol(),
    )
    .expect("a disc")
}

/// The dome against itself mirrored about `x = 2` (the axis reversed, a
/// frame whose `x` runs along `-x`): its arch on the other's, one surface
/// in two frames, refused (S9f.2a's (5)); here first as a vertex of one on
/// the other's face, the arches' ends meeting (the module's tests reach the
/// coincident arcs themselves).
#[test]
fn coincident_spline_walls_are_refused() {
    let a = prism(1, dome(), Frame3::xy(), 0.0, 5.0);
    let b = prism(
        2,
        dome(),
        frame((4.0, 0.0, 6.0), (0.0, 0.0, -1.0), (-1.0, 0.0, 0.0)),
        0.0,
        4.0,
    );
    for r in [
        a.fuse(OperationId(3), &b),
        a.cut(OperationId(4), &b),
        a.common(OperationId(5), &b),
    ] {
        assert!(
            matches!(&r, Err(Error::OutOfDomain(m)) if m.contains("one surface"))
                || matches!(&r, Err(Error::Degenerate(_))),
            "{:?}",
            r.map(|x| x.0.len())
        );
    }
}

/// Crossing axes stay refused: a cylinder's (S9f.2b's) and another spline
/// wall's (by design).
#[test]
fn crossing_axes_stay_refused() {
    let a = prism(1, dome(), Frame3::xy(), 0.0, 5.0);
    let side = frame((-1.0, 1.0, 2.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0));
    let b = prism(2, disc(1.0), side, 0.0, 6.0);
    let r = a.fuse(OperationId(3), &b);
    assert!(
        matches!(&r, Err(Error::OutOfDomain(m)) if m.contains("S9f.2b")),
        "{:?}",
        r.map(|x| x.0.len())
    );
    let c = prism(4, dome(), side, 0.0, 6.0);
    let r = a.fuse(OperationId(5), &c);
    assert!(
        matches!(&r, Err(Error::OutOfDomain(m)) if m.contains("crossing axes")),
        "{:?}",
        r.map(|x| x.0.len())
    );
}

/// Two crossing arches of degrees five and four (a field of degree 20):
/// `ComputationLimit` (S9f.2a's (4)); of degrees four and four (16) taken,
/// its operations consistent.
#[test]
fn arcs_whose_degrees_multiply_past_16_are_a_limit() {
    let quintic = arch(
        5.0,
        &[
            (5.0, 0.0),
            (4.0, 3.0),
            (3.0, 1.0),
            (2.0, 4.0),
            (1.0, 2.0),
            (0.0, 0.0),
        ],
    );
    let quartic = || {
        arch(
            4.0,
            &[(4.0, 0.0), (3.0, 3.0), (2.0, 1.0), (1.0, 3.0), (0.0, 0.0)],
        )
    };
    let turn = frame((3.75, -0.5, -1.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0));
    let a = prism(1, quintic, Frame3::xy(), 0.0, 3.0);
    let b = prism(2, quartic(), turn, 0.0, 5.0);
    let r = a.common(OperationId(3), &b);
    assert!(
        matches!(&r, Err(Error::ComputationLimit(m)) if m.contains("16")),
        "{:?}",
        r.map(|x| x.0.len())
    );
    let a = prism(4, quartic(), Frame3::xy(), 0.0, 3.0);
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let volume = |r: Vec<Solid>| r.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let fused = volume(a.fuse(OperationId(5), &b).expect("a fuse").0);
    let common = volume(a.common(OperationId(6), &b).expect("a common").0);
    let cut = volume(a.cut(OperationId(7), &b).expect("a cut").0);
    assert!((fused - (va + vb - common)).abs() <= 1e-9 * fused);
    assert!((cut - (va - common)).abs() <= 1e-9 * va);
    assert!(common > 0.0 && common < va.min(vb));
}

/// A disc crossing the wave at its interior knot `(6.5, 6.5)` (its circle
/// through the knot across the wave's legs there), a quarter turn of the
/// disc's frame against the wave's `XY`: the same Boolean as the disc in
/// `XY` (a circle is unchanged by the turn), which S9a.2's one-frame
/// profiles decide, within 1e-9.
#[test]
fn a_crossing_at_a_knot_is_the_one_frame_booleans() {
    let wave = || {
        let outer = Boundary::path(
            vec![
                Point2::new(0.0, 0.0),
                Point2::new(10.0, 0.0),
                Point2::new(10.0, 6.0),
                Point2::new(0.0, 6.0),
            ],
            vec![
                Segment::Line,
                Segment::Line,
                spline_segment(
                    &[(10.0, 6.0), (8.0, 8.0), (5.0, 5.0), (2.0, 8.0), (0.0, 6.0)],
                    vec![0.0, 1.0, 2.0, 3.0],
                    vec![3, 1, 1, 3],
                ),
                Segment::Line,
            ],
            tol(),
        )
        .expect("the wave");
        Profile::new(outer, vec![], tol()).expect("a profile")
    };
    let a = prism(1, wave(), Frame3::xy(), 0.0, 5.0);
    let turned = prism(
        2,
        disc(1.0),
        frame((6.5, 7.5, -1.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0)),
        0.0,
        4.0,
    );
    let flat = prism(
        3,
        disc(1.0),
        frame((6.5, 7.5, -1.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
        0.0,
        4.0,
    );
    let volume = |r: Vec<Solid>| r.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    for (k, op) in ["fuse", "cut", "common"].iter().enumerate() {
        let id = |n: u64| OperationId(10 + 2 * k as u64 + n);
        let run = |b: &Solid, n: u64| match *op {
            "fuse" => a.fuse(id(n), b),
            "cut" => a.cut(id(n), b),
            _ => a.common(id(n), b),
        };
        let x = volume(run(&turned, 0).expect("S9f.2a").0);
        let y = volume(run(&flat, 1).expect("S9a.2").0);
        assert!((x - y).abs() <= 1e-9 * x.max(1.0), "{op}: {x} against {y}");
    }
}
