//! S7c: lines and circles (S7c.1), ellipses, hyperbolas and splines (S7c.2)
//! against planes, cylinders, cones, spheres and tori, against the
//! independent reference (`fixtures/curve-surface-*.txt|tsv` from
//! `tools/generate_curve_surface_fixtures.py`).
#[path = "support/curve_surface_protocol.rs"]
mod protocol;
use protocol::{cases, rows, Case, Curve};
use rusty_occt::intersection::{
    conic_surface, curve_surface, spline_cone, spline_torus, Conic, CurveSurfaceIntersection,
};
use rusty_occt::topology::{Curve3, SplineSpan, Surface};
use rusty_occt::{BSplineCurve3, Error, Frame3, Point3, Tolerance, Vec3};
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

fn all() -> Vec<Case> {
    cases(include_str!("../../fixtures/curve-surface-cases.txt"))
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

fn curve_frame(c: &Curve) -> Option<Frame3> {
    match c {
        Curve::Edge(Curve3::Circle { frame, .. })
        | Curve::Conic(Conic::Ellipse { frame, .. })
        | Curve::Conic(Conic::Hyperbola { frame, .. }) => Some(*frame),
        _ => None,
    }
}

fn surface_frame(s: &Surface) -> Frame3 {
    match s {
        Surface::Plane(f)
        | Surface::Cylinder { frame: f, .. }
        | Surface::Cone { frame: f, .. }
        | Surface::Sphere { frame: f, .. }
        | Surface::Torus { frame: f, .. } => *f,
        Surface::BSpline(_) => unreachable!(),
    }
}

/// The kernel stores the frames the reference intersected, bit for bit: the
/// normals, and a conic's axes.
#[test]
fn stored_frames_are_the_reference_inputs() {
    let all: BTreeMap<String, Case> = all().into_iter().map(|c| (c.name.clone(), c)).collect();
    let rows: Vec<&str> = include_str!("../../fixtures/curve-surface-frames.tsv")
        .lines()
        .skip(1)
        .collect();
    let mut checked = 0;
    for row in rows {
        let words: Vec<&str> = row.split('\t').collect();
        let case = &all[words[0]];
        let own = curve_frame(&case.curve);
        let got = match words[1] {
            "x" => own.unwrap().x(),
            "y" => own.unwrap().y(),
            k => {
                let framed: Vec<Frame3> = own
                    .into_iter()
                    .chain([surface_frame(&case.surface)])
                    .collect();
                framed[k.parse::<usize>().unwrap()].normal()
            }
        };
        let want: Vec<u64> = words[2]
            .split(' ')
            .map(|h| u64::from_str_radix(h, 16).unwrap())
            .collect();
        assert_eq!(
            got.to_array().map(f64::to_bits).to_vec(),
            want,
            "{}",
            words[0]
        );
        checked += 1;
    }
    let conics = all
        .values()
        .filter(|c| matches!(c.curve, Curve::Conic(_)))
        .count();
    let framed = all
        .values()
        .filter(|c| curve_frame(&c.curve).is_some())
        .count();
    assert_eq!(checked, all.len() + framed + 2 * conics);
}

fn inside(x: f64, [lo, hi]: [f64; 2], turn: bool) -> bool {
    let slack = 1e-25 * x.abs();
    let shifts: &[f64] = if turn { &[0.0, TAU, -TAU] } else { &[0.0] };
    shifts
        .iter()
        .any(|s| lo - slack <= x + s && x + s <= hi + slack)
}

fn check(case: &Case, got: &[String], want: &[Vec<String>]) -> Result<(), String> {
    let name = &case.name;
    let turn = matches!(
        case.curve,
        Curve::Edge(Curve3::Circle { .. }) | Curve::Conic(Conic::Ellipse { .. })
    );
    if got.len() != want.len() {
        return Err(format!("{name}: {got:?} for {want:?}"));
    }
    for (g, w) in got.iter().zip(want) {
        let g: Vec<&str> = g.split(' ').collect();
        if g[0] != w[0] {
            return Err(format!("{name}: {g:?} for {w:?}"));
        }
        match w[0].as_str() {
            "empty" | "contained" => {}
            "overlap" => {
                let (a, b): (f64, f64) = (w[1].parse().unwrap(), w[2].parse().unwrap());
                if g[1].parse::<f64>().unwrap() != a || g[2].parse::<f64>().unwrap() != b {
                    return Err(format!("{name}: {g:?} for {w:?}"));
                }
            }
            _ => {
                let numbers: Vec<f64> = w[1..5].iter().map(|x| x.parse().unwrap()).collect();
                let bound = |k: usize| -> [f64; 2] {
                    [g[1 + 2 * k].parse().unwrap(), g[2 + 2 * k].parse().unwrap()]
                };
                if g[9] != w[5] {
                    return Err(format!("{name}: contact {g:?} for {w:?}"));
                }
                if !inside(numbers[0], bound(0), turn)
                    || (1..4).any(|k| !inside(numbers[k], bound(k), false))
                {
                    return Err(format!("{name}: {g:?} misses {w:?}"));
                }
            }
        }
    }
    Ok(())
}

#[test]
fn every_case_matches_the_exact_reference() {
    let want = expected();
    let all = all();
    assert_eq!(all.len(), want.len());
    let mut failures = Vec::new();
    for case in &all {
        let got = rows(case).unwrap_or_else(|e| panic!("{}: {e}", case.name));
        if let Err(why) = check(case, &got, &want[&case.name]) {
            failures.push(why);
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// A point's curve position at its parameter, in binary64.
fn at(c: &Curve, t: f64) -> Option<Point3> {
    Some(match c {
        Curve::Edge(Curve3::LineSegment { start, end }) => *start + (*end - *start) * t,
        Curve::Edge(Curve3::Circle { frame, radius }) => {
            frame.origin() + (frame.x() * t.cos() + frame.y() * t.sin()) * *radius
        }
        Curve::Conic(Conic::Ellipse {
            frame,
            major,
            minor,
        }) => frame.origin() + frame.x() * (major * t.cos()) + frame.y() * (minor * t.sin()),
        Curve::Conic(Conic::Hyperbola {
            frame,
            major,
            minor,
        }) => frame.origin() + frame.x() * (major * t.cosh()) + frame.y() * (minor * t.sinh()),
        Curve::Spline(s) => s.point(t).ok()?,
        _ => return None,
    })
}

/// An arc is its whole circle; every point is its curve's point at its
/// parameter; points are sorted.
#[test]
fn arcs_are_their_circles_and_points_lie_on_the_curve() {
    for case in all() {
        let got = rows(&case).unwrap();
        if let Curve::Edge(Curve3::Circle { frame, radius }) = case.curve {
            let arc = Curve3::CircularArc {
                frame,
                radius,
                start_angle: 0.3,
                sweep_angle: 1.0,
            };
            let circle = curve_surface(&Curve3::Circle { frame, radius }, &case.surface).unwrap();
            assert_eq!(
                curve_surface(&arc, &case.surface).unwrap(),
                circle,
                "{}",
                case.name
            );
        }
        let mut last = f64::NEG_INFINITY;
        for row in got {
            let w: Vec<&str> = row.split(' ').collect();
            if w[0] != "point" {
                continue;
            }
            let n: Vec<f64> = w[1..9].iter().map(|x| x.parse().unwrap()).collect();
            assert!(n[0] >= last, "{}: unsorted", case.name);
            last = n[1];
            let t = 0.5 * n[0] + 0.5 * n[1];
            let Some(p) = at(&case.curve, t) else {
                continue;
            };
            let mid = [
                0.5 * (n[2] + n[3]),
                0.5 * (n[4] + n[5]),
                0.5 * (n[6] + n[7]),
            ];
            let gap = (0..3)
                .map(|k| (p.to_array()[k] - mid[k]).abs())
                .fold(0.0, f64::max);
            assert!(gap < 1e-9, "{}: {gap}", case.name);
        }
    }
}

#[test]
fn unsupported_and_degenerate_inputs_are_errors() {
    let frame = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Tolerance::default(),
    )
    .unwrap();
    let plane = Surface::Plane(frame);
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
    let edge = Curve3::BSpline(SplineSpan::whole(spline.clone()));
    assert!(matches!(
        curve_surface(&edge, &plane),
        Err(Error::OutOfDomain(_))
    ));
    assert!(matches!(
        spline_torus(&spline, &plane),
        Err(Error::OutOfDomain(_))
    ));
    assert!(matches!(
        spline_cone(&spline, &plane),
        Err(Error::OutOfDomain(_))
    ));
    // An ellipse's major semi-axis is the larger, as OCCT's.
    for (major, minor) in [(1.0, 2.0), (0.0, 0.0), (-1.0, 1.0)] {
        let e = Conic::Ellipse {
            frame,
            major,
            minor,
        };
        assert!(matches!(
            conic_surface(&e, &plane),
            Err(Error::OutOfDomain(_))
        ));
    }
    let nan = Conic::Hyperbola {
        frame,
        major: f64::NAN,
        minor: 1.0,
    };
    assert!(matches!(
        conic_surface(&nan, &plane),
        Err(Error::NonFinite(_))
    ));
}

/// Lines through a sphere's centre at many rational directions: two
/// crossings at distance `r` each side, symmetric in the parameter.
#[test]
fn lines_through_a_sphere_cross_it_twice() {
    let sphere = Surface::Sphere {
        frame: Frame3::new(
            Point3::new(1.0, -2.0, 0.5),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Tolerance::default(),
        )
        .unwrap(),
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

/// An ellipse about a sphere's centre in a plane through it, with semi-axes
/// either side of the radius: four crossings, symmetric in the angle; a
/// hyperbola's branch meets it twice where its vertex is inside.
#[test]
fn conics_about_a_sphere_centre_cross_it_symmetrically() {
    let o = Point3::new(0.5, 0.25, -1.0);
    let frame = Frame3::new(
        o,
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        Tolerance::default(),
    )
    .unwrap();
    let sphere = Surface::Sphere { frame, radius: 1.5 };
    for k in 1..12 {
        let minor = 0.125 * f64::from(k);
        let e = Conic::Ellipse {
            frame,
            major: 2.0,
            minor,
        };
        let got = conic_surface(&e, &sphere).unwrap();
        let CurveSurfaceIntersection::Points(points) = got else {
            panic!("{k}: {got:?}");
        };
        assert_eq!(points.len(), 4, "{k}");
        let h = Conic::Hyperbola {
            frame,
            major: 1.0,
            minor,
        };
        let CurveSurfaceIntersection::Points(points) = conic_surface(&h, &sphere).unwrap() else {
            panic!("{k}");
        };
        assert_eq!(points.len(), 2, "{k}");
        let t = [0, 1].map(|i| 0.5 * points[i].parameter[0] + 0.5 * points[i].parameter[1]);
        assert!((t[0] + t[1]).abs() < 1e-12, "{k}: {t:?}");
    }
}
