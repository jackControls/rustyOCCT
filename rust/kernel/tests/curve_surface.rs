//! S7c.1: lines and circles against planes, cylinders, cones, spheres and
//! tori, against the independent reference (`fixtures/curve-surface-*.txt|tsv`
//! from `tools/generate_curve_surface_fixtures.py`).
use rusty_occt::intersection::{curve_surface, CurvePoint, CurveSurfaceIntersection};
use rusty_occt::topology::{Curve3, SplineSpan, Surface};
use rusty_occt::{BSplineCurve3, Error, Frame3, Point3, Tolerance, Vec3};
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

fn frame(v: &[f64]) -> Frame3 {
    Frame3::new(
        Point3::new(v[0], v[1], v[2]),
        Vec3::new(v[3], v[4], v[5]),
        Vec3::new(v[6], v[7], v[8]),
        Tolerance::default(),
    )
    .unwrap()
}

fn numbers(words: &[&str]) -> Vec<f64> {
    words.iter().map(|w| w.parse().unwrap()).collect()
}

fn curve(words: &[&str]) -> Curve3 {
    let v = numbers(&words[1..]);
    match words[0] {
        "line" => Curve3::LineSegment {
            start: Point3::new(v[0], v[1], v[2]),
            end: Point3::new(v[3], v[4], v[5]),
        },
        _ => Curve3::Circle {
            frame: frame(&v),
            radius: v[9],
        },
    }
}

fn surface(words: &[&str]) -> Surface {
    let v = numbers(&words[1..]);
    let frame = frame(&v);
    match words[0] {
        "plane" => Surface::Plane(frame),
        "sphere" => Surface::Sphere {
            frame,
            radius: v[9],
        },
        "cylinder" => Surface::Cylinder {
            frame,
            radius: v[9],
        },
        "torus" => Surface::Torus {
            frame,
            major: v[9],
            minor: v[10],
        },
        _ => Surface::Cone {
            frame,
            radius: v[9],
            half_angle: v[10],
        },
    }
}

fn cases() -> Vec<(String, Curve3, Surface)> {
    include_str!("../../fixtures/curve-surface-cases.txt")
        .split("\nend")
        .filter(|b| !b.trim().is_empty())
        .map(|block| {
            let lines: Vec<Vec<&str>> = block
                .trim()
                .lines()
                .map(|l| l.split_whitespace().collect())
                .collect();
            (
                lines[0][1].to_string(),
                curve(&lines[1][1..]),
                surface(&lines[2][1..]),
            )
        })
        .collect()
}

fn expected() -> BTreeMap<String, Vec<Vec<String>>> {
    let mut out: BTreeMap<String, Vec<Vec<String>>> = BTreeMap::new();
    for line in include_str!("../../fixtures/curve-surface-expected.tsv")
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

fn stored(c: &Curve3) -> Option<&Frame3> {
    match c {
        Curve3::Circle { frame, .. } | Curve3::CircularArc { frame, .. } => Some(frame),
        _ => None,
    }
}

fn surface_frame(s: &Surface) -> &Frame3 {
    match s {
        Surface::Plane(f)
        | Surface::Cylinder { frame: f, .. }
        | Surface::Cone { frame: f, .. }
        | Surface::Sphere { frame: f, .. }
        | Surface::Torus { frame: f, .. } => f,
        Surface::BSpline(_) => unreachable!(),
    }
}

/// The kernel stores the frames the reference intersected, bit for bit.
#[test]
fn stored_normals_are_the_reference_inputs() {
    let all: BTreeMap<String, (Curve3, Surface)> =
        cases().into_iter().map(|(n, c, s)| (n, (c, s))).collect();
    let rows: Vec<&str> = include_str!("../../fixtures/curve-surface-frames.tsv")
        .lines()
        .skip(1)
        .collect();
    let circles = all.values().filter(|(c, _)| stored(c).is_some()).count();
    assert_eq!(rows.len(), all.len() + circles);
    for row in rows {
        let words: Vec<&str> = row.split('\t').collect();
        let (c, s) = &all[words[0]];
        // A circle's frame, then the surface's.
        let framed: Vec<&Frame3> = stored(c).into_iter().chain([surface_frame(s)]).collect();
        let frame = framed[words[1].parse::<usize>().unwrap()];
        let want: Vec<u64> = words[2]
            .split(' ')
            .map(|h| u64::from_str_radix(h, 16).unwrap())
            .collect();
        let got: Vec<u64> = frame.normal().to_array().map(f64::to_bits).to_vec();
        assert_eq!(got, want, "{}", words[0]);
    }
}

fn inside(x: f64, [lo, hi]: [f64; 2], turn: bool) -> bool {
    let slack = 1e-25 * x.abs();
    let shifts: &[f64] = if turn { &[0.0, TAU, -TAU] } else { &[0.0] };
    shifts
        .iter()
        .any(|s| lo - slack <= x + s && x + s <= hi + slack)
}

fn check(
    name: &str,
    circle: bool,
    got: &CurveSurfaceIntersection,
    rows: &[Vec<String>],
) -> Result<(), String> {
    match got {
        CurveSurfaceIntersection::Empty if rows[0][0] == "empty" => Ok(()),
        CurveSurfaceIntersection::Contained if rows[0][0] == "contained" => Ok(()),
        CurveSurfaceIntersection::Points(points) if points.len() == rows.len() => {
            for (p, row) in points.iter().zip(rows) {
                let want: Vec<f64> = row[1..5].iter().map(|w| w.parse().unwrap()).collect();
                let contact = if p.tangent { "tangent" } else { "crossing" };
                if row[0] != "point" || row[5] != contact {
                    return Err(format!("{name}: {contact} for {row:?}"));
                }
                if !inside(want[0], p.parameter, circle)
                    || (0..3).any(|k| !inside(want[1 + k], p.point[k], false))
                {
                    return Err(format!("{name}: {p:?} misses {row:?}"));
                }
            }
            Ok(())
        }
        other => Err(format!("{name}: {other:?} for {rows:?}")),
    }
}

#[test]
fn every_case_matches_the_exact_reference() {
    let want = expected();
    let all = cases();
    assert_eq!(all.len(), want.len());
    let mut failures = Vec::new();
    for (name, c, s) in &all {
        let got = curve_surface(c, s).unwrap_or_else(|e| panic!("{name}: {e}"));
        if let Err(why) = check(name, stored(c).is_some(), &got, &want[name]) {
            failures.push(why);
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// An arc is its whole circle, and every point lies on the curve.
#[test]
fn arcs_are_their_circles_and_points_lie_on_the_curve() {
    for (name, c, s) in cases() {
        let got = curve_surface(&c, &s).unwrap();
        if let Curve3::Circle { frame, radius } = c {
            let arc = Curve3::CircularArc {
                frame,
                radius,
                start_angle: 0.3,
                sweep_angle: 1.0,
            };
            assert_eq!(curve_surface(&arc, &s).unwrap(), got, "{name}");
        }
        let CurveSurfaceIntersection::Points(points) = got else {
            continue;
        };
        for w in points.windows(2) {
            assert!(w[0].parameter[1] <= w[1].parameter[0], "{name}: unsorted");
        }
        for CurvePoint {
            parameter, point, ..
        } in points
        {
            let t = 0.5 * parameter[0] + 0.5 * parameter[1];
            let at = match &c {
                Curve3::LineSegment { start, end } => *start + (*end - *start) * t,
                Curve3::Circle { frame, radius } => {
                    frame.origin() + (frame.x() * t.cos() + frame.y() * t.sin()) * *radius
                }
                _ => unreachable!(),
            };
            let mid = point.map(|[lo, hi]| 0.5 * lo + 0.5 * hi);
            let gap = (0..3)
                .map(|k| (at.to_array()[k] - mid[k]).abs())
                .fold(0.0, f64::max);
            assert!(gap < 1e-9, "{name}: {gap}");
        }
    }
}

#[test]
fn unsupported_and_degenerate_inputs_are_errors() {
    let plane = Surface::Plane(frame(&[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]));
    let point = Point3::new(1.0, 2.0, 3.0);
    let degenerate = Curve3::LineSegment {
        start: point,
        end: point,
    };
    assert!(matches!(
        curve_surface(&degenerate, &plane),
        Err(Error::Degenerate(_))
    ));
    let poles = vec![Point3::new(0.0, 0.0, -1.0), Point3::new(1.0, 0.0, 1.0)];
    let spline = BSplineCurve3::new(1, poles, None, vec![0.0, 1.0], vec![2, 2]).unwrap();
    let edge = Curve3::BSpline(SplineSpan::whole(spline));
    assert!(matches!(
        curve_surface(&edge, &plane),
        Err(Error::OutOfDomain(_))
    ));
}

/// Lines through a sphere's centre at many rational directions: two
/// crossings at distance `r` each side, symmetric in the parameter.
#[test]
fn lines_through_a_sphere_cross_it_twice() {
    let sphere = Surface::Sphere {
        frame: frame(&[1.0, -2.0, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]),
        radius: 1.5,
    };
    for k in 0..40 {
        let a = PI * f64::from(k) / 40.0;
        let d = [a.cos(), a.sin(), 0.25 * f64::from(k % 7) - 0.75];
        let start = Point3::new(1.0 - d[0], -2.0 - d[1], 0.5 - d[2]);
        let end = Point3::new(1.0 + d[0], -2.0 + d[1], 0.5 + d[2]);
        let got = curve_surface(&Curve3::LineSegment { start, end }, &sphere).unwrap();
        let CurveSurfaceIntersection::Points(points) = got else {
            panic!("{k}: {got:?}");
        };
        assert_eq!(points.len(), 2, "{k}");
        assert!(points.iter().all(|p| !p.tangent));
        let s = [0, 1].map(|i| 0.5 * points[i].parameter[0] + 0.5 * points[i].parameter[1]);
        assert!((s[0] + s[1] - 1.0).abs() < 1e-12, "{k}: {s:?}");
    }
}
