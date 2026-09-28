//! S7d.1: pairs of lines, circles, ellipses and hyperbolas' branches,
//! against the independent reference (`fixtures/curve-curve-*.txt|tsv` from
//! `tools/generate_curve_curve_fixtures.py`).
#[path = "support/curve_curve_protocol.rs"]
mod protocol;
use protocol::{cases, rows, Case};
use rusty_occt::intersection::{curve_curve, AnalyticCurve, Conic, CurveCurveIntersection};
use rusty_occt::topology::{Curve3, SplineSpan};
use rusty_occt::{BSplineCurve3, Error, Frame3, Point3, Tolerance, Vec3};
use std::collections::BTreeMap;
use std::f64::consts::TAU;

fn all() -> Vec<Case> {
    cases(include_str!("../../fixtures/curve-curve-cases.txt"))
}

fn expected() -> BTreeMap<String, Vec<Vec<String>>> {
    let mut out: BTreeMap<String, Vec<Vec<String>>> = BTreeMap::new();
    for line in include_str!("../../fixtures/curve-curve-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        out.entry(name.to_string())
            .or_default()
            .push(row.split(' ').map(str::to_string).collect());
    }
    out
}

fn frame_of(c: &AnalyticCurve) -> Option<Frame3> {
    match c {
        AnalyticCurve::Edge(Curve3::Circle { frame, .. })
        | AnalyticCurve::Conic(Conic::Ellipse { frame, .. })
        | AnalyticCurve::Conic(Conic::Hyperbola { frame, .. }) => Some(*frame),
        _ => None,
    }
}

/// The kernel stores the frames the reference intersected, bit for bit.
#[test]
fn stored_frames_are_the_reference_inputs() {
    let all: BTreeMap<String, Case> = all().into_iter().map(|c| (c.name.clone(), c)).collect();
    let mut checked = 0;
    for row in include_str!("../../fixtures/curve-curve-frames.tsv")
        .lines()
        .skip(1)
    {
        let w: Vec<&str> = row.split('\t').collect();
        let frame = frame_of(&all[w[0]].curves[w[1].parse::<usize>().unwrap()]).unwrap();
        let got = match w[2] {
            "n" => frame.normal(),
            "x" => frame.x(),
            _ => frame.y(),
        };
        let want: Vec<u64> = w[3]
            .split(' ')
            .map(|h| u64::from_str_radix(h, 16).unwrap())
            .collect();
        assert_eq!(got.to_array().map(f64::to_bits).to_vec(), want, "{}", w[0]);
        checked += 1;
    }
    let framed: usize = all
        .values()
        .map(|c| c.curves.iter().filter(|x| frame_of(x).is_some()).count())
        .sum();
    assert_eq!(checked, 3 * framed);
}

fn inside(x: f64, [lo, hi]: [f64; 2], turn: bool) -> bool {
    let slack = 1e-25 * x.abs();
    let shifts: &[f64] = if turn { &[0.0, TAU, -TAU] } else { &[0.0] };
    shifts
        .iter()
        .any(|s| lo - slack <= x + s && x + s <= hi + slack)
}

fn periodic(c: &AnalyticCurve) -> bool {
    matches!(
        c,
        AnalyticCurve::Edge(Curve3::Circle { .. }) | AnalyticCurve::Conic(Conic::Ellipse { .. })
    )
}

#[test]
fn every_case_matches_the_exact_reference() {
    let want = expected();
    let all = all();
    assert_eq!(all.len(), want.len());
    let mut failures = Vec::new();
    for case in &all {
        let got = rows(case).unwrap_or_else(|e| panic!("{}: {e}", case.name));
        let w = &want[&case.name];
        let ok = got.len() == w.len()
            && got.iter().zip(w).all(|(g, w)| {
                let g: Vec<&str> = g.split(' ').collect();
                if g[0] != w[0] {
                    return false;
                }
                if w[0] != "point" {
                    return true;
                }
                let bound = |k: usize| -> [f64; 2] {
                    [g[1 + 2 * k].parse().unwrap(), g[2 + 2 * k].parse().unwrap()]
                };
                g[11] == w[6]
                    && (0..5).all(|k| {
                        let turn = k < 2 && periodic(&case.curves[k]);
                        inside(w[1 + k].parse().unwrap(), bound(k), turn)
                    })
            });
        if !ok {
            failures.push(format!("{}: {got:?} for {w:?}", case.name));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The order of the curves swaps the parameters and keeps the rest (two
/// conics in one plane are computed through the first: the enclosures may
/// differ, the values not).
#[test]
fn swapping_the_curves_swaps_the_parameters() {
    for case in all() {
        let [a, b] = &case.curves;
        let (x, y) = (curve_curve(a, b).unwrap(), curve_curve(b, a).unwrap());
        match (&x, &y) {
            (CurveCurveIntersection::Points(p), CurveCurveIntersection::Points(q)) => {
                assert_eq!(p.len(), q.len(), "{}", case.name);
                for p in p {
                    let hit = q.iter().any(|q| {
                        let overlap = |[a, b]: [f64; 2], [c, d]: [f64; 2]| a <= d && c <= b;
                        q.tangent == p.tangent
                            && overlap(q.parameters[0], p.parameters[1])
                            && overlap(q.parameters[1], p.parameters[0])
                            && (0..3).all(|k| {
                                q.point[k][0] <= p.point[k][1] && p.point[k][0] <= q.point[k][1]
                            })
                    });
                    assert!(hit, "{}: {p:?} in {q:?}", case.name);
                }
            }
            _ => assert_eq!(x, y, "{}", case.name),
        }
    }
}

#[test]
fn splines_and_degenerate_lines_are_refused() {
    let frame = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Tolerance::default(),
    )
    .unwrap();
    let circle = AnalyticCurve::Edge(Curve3::Circle { frame, radius: 1.0 });
    let p = Point3::new(1.0, 2.0, 3.0);
    let point = AnalyticCurve::Edge(Curve3::LineSegment { start: p, end: p });
    assert!(matches!(
        curve_curve(&point, &circle),
        Err(Error::Degenerate(_))
    ));
    let poles = vec![Point3::new(0.0, 0.0, -1.0), Point3::new(1.0, 0.0, 1.0)];
    let spline = BSplineCurve3::new(1, poles, None, vec![0.0, 1.0], vec![2, 2]).unwrap();
    let edge = AnalyticCurve::Edge(Curve3::BSpline(SplineSpan::whole(spline)));
    assert!(matches!(
        curve_curve(&edge, &circle),
        Err(Error::OutOfDomain(_))
    ));
}

/// Circles of a sphere through its poles (great circles in planes through
/// its axis at rational angles) meet at both poles, crossing.
#[test]
fn great_circles_through_the_poles_cross_there() {
    let o = Point3::new(0.5, -0.25, 1.0);
    let polar = |n: Vec3| {
        AnalyticCurve::Edge(Curve3::Circle {
            frame: Frame3::new(o, n, Vec3::new(0.0, 0.0, 1.0), Tolerance::default()).unwrap(),
            radius: 2.0,
        })
    };
    let normals = [
        (1.0, 0.0),
        (0.0, 1.0),
        (3.0, 4.0),
        (-4.0, 3.0),
        (1.0, 1.0),
        (5.0, -12.0),
    ];
    for (i, &(a, b)) in normals.iter().enumerate() {
        for &(c, d) in &normals[i + 1..] {
            let got =
                curve_curve(&polar(Vec3::new(a, b, 0.0)), &polar(Vec3::new(c, d, 0.0))).unwrap();
            let CurveCurveIntersection::Points(points) = got else {
                panic!("{a} {b} {c} {d}: {got:?}");
            };
            assert_eq!(points.len(), 2);
            for p in points {
                assert!(!p.tangent);
                let z = 0.5 * (p.point[2][0] + p.point[2][1]);
                assert!((z - 1.0).abs() > 1.99, "{z}");
            }
        }
    }
}
