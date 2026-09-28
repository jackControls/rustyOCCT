//! Pairs of lines, circles, ellipses and hyperbolas (S7d.1 of
//! REVIEW_NOTES.md) on dyadic data, independent or with an exact degeneracy
//! made on purpose: the second curve in the first's plane; the same curve in
//! another representation (a line through two other of its points, a circle
//! with its normal reversed, an ellipse or a hyperbola with its axis
//! reversed, which is the other branch); two circles touching in their plane;
//! a line along a conic's axis. The intersection never panics and fails only
//! by a certified comparison it cannot decide; swapping the curves swaps the
//! parameters; every point is each curve's point at its parameter; the
//! points of a coincident pair lie on both curves.
use crate::analytic_intersections::{frame, Bytes};
use rusty_occt::intersection::{
    curve_curve, spline_curve, AnalyticCurve, Conic, CurveCurveIntersection,
};
use rusty_occt::topology::Curve3;
use rusty_occt::BSplineCurve3;
use rusty_occt::{Error, Frame3, Point3, Vec3};

fn mid([lo, hi]: [f64; 2]) -> f64 {
    0.5 * lo + 0.5 * hi
}

fn at(c: &AnalyticCurve, t: f64) -> Point3 {
    match c {
        AnalyticCurve::Edge(Curve3::LineSegment { start, end }) => *start + (*end - *start) * t,
        AnalyticCurve::Edge(Curve3::Circle { frame, radius }) => {
            frame.origin() + (frame.x() * t.cos() + frame.y() * t.sin()) * *radius
        }
        AnalyticCurve::Conic(Conic::Ellipse {
            frame,
            major,
            minor,
        }) => frame.origin() + frame.x() * (major * t.cos()) + frame.y() * (minor * t.sin()),
        AnalyticCurve::Conic(Conic::Hyperbola {
            frame,
            major,
            minor,
        }) => frame.origin() + frame.x() * (major * t.cosh()) + frame.y() * (minor * t.sinh()),
        _ => unreachable!(),
    }
}

/// A binary64 residual of a point against a curve (zero on it).
fn residual(c: &AnalyticCurve, p: Point3) -> f64 {
    match c {
        AnalyticCurve::Edge(Curve3::LineSegment { start, end }) => {
            let d = *end - *start;
            (p - *start).cross(d).length() / d.length()
        }
        AnalyticCurve::Edge(Curve3::Circle { frame, radius }) => {
            let w = p - frame.origin();
            let h = w.dot(frame.normal());
            h.abs() + ((w - frame.normal() * h).length() - radius).abs()
        }
        AnalyticCurve::Conic(conic) => {
            let (f, a, b, sign): (Frame3, f64, f64, f64) = match *conic {
                Conic::Ellipse {
                    frame,
                    major,
                    minor,
                } => (frame, major, minor, 1.0),
                Conic::Hyperbola {
                    frame,
                    major,
                    minor,
                } => (frame, major, minor, -1.0),
            };
            let w = p - f.origin();
            let (xi, eta) = (w.dot(f.x()), w.dot(f.y()));
            w.dot(f.normal()).abs()
                + (xi * xi / (a * a) + sign * eta * eta / (b * b) - 1.0).abs() * a.min(b)
        }
        _ => unreachable!(),
    }
}

fn make(kind: u8, f: Frame3, a: f64, b: f64) -> AnalyticCurve {
    match kind % 4 {
        0 => AnalyticCurve::Edge(Curve3::LineSegment {
            start: f.origin(),
            end: f.origin() + f.x() * a,
        }),
        1 => AnalyticCurve::Edge(Curve3::Circle {
            frame: f,
            radius: a,
        }),
        2 => AnalyticCurve::Conic(Conic::Ellipse {
            frame: f,
            major: a.max(b),
            minor: a.min(b),
        }),
        _ => AnalyticCurve::Conic(Conic::Hyperbola {
            frame: f,
            major: a,
            minor: b,
        }),
    }
}

pub fn check_curve_curve(data: &[u8]) {
    let mut b = Bytes(data, 0);
    let (k1, k2, mode) = (b.next(), b.next(), b.next() % 5);
    let o = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
    let n = Vec3::new(b.small(), b.small(), b.small());
    let (a1, b1) = (
        0.25 + f64::from(b.next() % 32) / 8.0,
        0.25 + f64::from(b.next() % 32) / 8.0,
    );
    let Some(f1) = frame(o, n) else { return };
    if k1 >= 128 {
        return check_spline(&mut b, k2, mode, make(k2 % 3 + 1, f1, a1, b1));
    }
    let c1 = make(k1, f1, a1, b1);
    let (a2, b2) = (
        0.25 + f64::from(b.next() % 32) / 8.0,
        0.25 + f64::from(b.next() % 32) / 8.0,
    );
    let c2 = match mode {
        0 => {
            let o2 = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
            let n2 = Vec3::new(b.small(), b.small(), b.small());
            let Some(f2) = frame(o2, n2) else { return };
            make(k2, f2, a2, b2)
        }
        // In the first's plane (a circle's: normal to its stored normal).
        1 => {
            let o2 = o + f1.x() * b.dyadic() + f1.y() * b.dyadic();
            let Some(f2) = frame(o2, f1.normal()) else {
                return;
            };
            make(k2, f2, a2, b2)
        }
        // The same curve in another representation.
        2 => match &c1 {
            AnalyticCurve::Edge(Curve3::LineSegment { start, end }) => {
                let d = *end - *start;
                AnalyticCurve::Edge(Curve3::LineSegment {
                    start: *start + d * 3.0,
                    end: *start + d * -0.5,
                })
            }
            AnalyticCurve::Edge(Curve3::Circle { radius, .. }) => {
                let Some(f2) = frame(o, -f1.normal()) else {
                    return;
                };
                AnalyticCurve::Edge(Curve3::Circle {
                    frame: f2,
                    radius: *radius,
                })
            }
            AnalyticCurve::Conic(c) => {
                let (major, minor, hyperbola) = match *c {
                    Conic::Ellipse { major, minor, .. } => (major, minor, false),
                    Conic::Hyperbola { major, minor, .. } => (major, minor, true),
                };
                let Ok(f2) = Frame3::new(o, f1.normal(), -f1.x(), rusty_occt::Tolerance::default())
                else {
                    return;
                };
                if hyperbola {
                    AnalyticCurve::Conic(Conic::Hyperbola {
                        frame: f2,
                        major,
                        minor,
                    })
                } else {
                    AnalyticCurve::Conic(Conic::Ellipse {
                        frame: f2,
                        major,
                        minor,
                    })
                }
            }
            _ => return,
        },
        // Two circles touching in the first's plane.
        3 => {
            let Some(f2) = frame(o + f1.x() * (a1 + a2), f1.normal()) else {
                return;
            };
            make(1, f2, a2, b2)
        }
        // A line along the first curve's x axis.
        _ => AnalyticCurve::Edge(Curve3::LineSegment {
            start: o + f1.x() * b.dyadic(),
            end: o + f1.x() * (b.dyadic() + 1.0),
        }),
    };
    let result = match curve_curve(&c1, &c2) {
        Ok(r) => r,
        Err(Error::ComputationLimit(_)) => return,
        Err(Error::Degenerate(_)) if matches!(c2, AnalyticCurve::Edge(Curve3::LineSegment { start, end }) if start == end) => {
            return
        }
        Err(e) => panic!("unexpected error {e}"),
    };
    let swapped = curve_curve(&c2, &c1).expect("the swapped pair decides too");
    let scale = |p: Point3| (p - Point3::ORIGIN).length().max(1.0) * 8.0;
    match (&result, &swapped) {
        (CurveCurveIntersection::Points(p), CurveCurveIntersection::Points(q)) => {
            assert_eq!(p.len(), q.len(), "swapped");
            for w in p.windows(2) {
                assert!(mid(w[0].parameters[0]) <= mid(w[1].parameters[0]), "sorted");
            }
            for p in p {
                // Two conics in one plane are computed through the first:
                // the enclosures differ, the values and contacts do not.
                let overlap = |[a, b]: [f64; 2], [c, d]: [f64; 2]| a <= d && c <= b;
                assert!(
                    q.iter().any(|q| q.tangent == p.tangent
                        && overlap(q.parameters[0], p.parameters[1])
                        && overlap(q.parameters[1], p.parameters[0])),
                    "swapped parameters"
                );
                let x = Point3::new(mid(p.point[0]), mid(p.point[1]), mid(p.point[2]));
                for (k, c) in [&c1, &c2].into_iter().enumerate() {
                    let on = at(c, mid(p.parameters[k]));
                    assert!(
                        (on - x).length() <= 1e-9 * scale(x),
                        "{p:?} is not curve {k}'s point"
                    );
                }
            }
        }
        (CurveCurveIntersection::Coincident, CurveCurveIntersection::Coincident) => {
            for k in 0..5 {
                let p = at(&c1, -1.0 + 0.5 * f64::from(k));
                let r = residual(&c2, p);
                assert!(r <= 1e-9 * scale(p), "{r}: coincident {c1:?} and {c2:?}");
            }
        }
        _ => assert_eq!(result, swapped, "swapped"),
    }
}

/// S7d.2: a rational B-spline (a first kind byte of 128 or more) against a
/// circle, an ellipse or a hyperbola: independent poles, poles in the conic's
/// plane, or the exact rational quarter circle on a circle (then independent
/// poles).
fn check_spline(b: &mut Bytes, pick: u8, mode: u8, conic: AnalyticCurve) {
    let f = match &conic {
        AnalyticCurve::Edge(Curve3::Circle { frame, .. }) => *frame,
        AnalyticCurve::Conic(Conic::Ellipse { frame, .. } | Conic::Hyperbola { frame, .. }) => {
            *frame
        }
        _ => unreachable!(),
    };
    let (o, x, y) = (f.origin(), f.x(), f.y());
    let degree = 1 + usize::from(pick / 4 % 3);
    let n = degree + 1 + usize::from(b.next() % 3);
    let mut poles: Vec<Point3> = match (mode % 3, &conic) {
        (2, AnalyticCurve::Edge(Curve3::Circle { radius, .. })) => {
            vec![o + x * *radius, o + (x + y) * *radius, o + y * *radius]
        }
        (1, _) => Vec::new(),
        _ => (0..n)
            .map(|_| Point3::new(b.dyadic(), b.dyadic(), b.dyadic()))
            .collect(),
    };
    let quarter = poles.len() == 3 && mode % 3 == 2;
    while poles.len() < n.max(3) {
        poles.push(if mode % 3 == 1 {
            o + x * b.dyadic() + y * b.dyadic()
        } else {
            Point3::new(b.dyadic(), b.dyadic(), b.dyadic())
        });
    }
    let n = poles.len();
    let degree = if quarter { 2 } else { degree.min(n - 1) };
    let weights: Vec<f64> = (0..n)
        .map(|i| match (quarter, i) {
            (true, 2) => 2.0,
            (true, 0 | 1) => 1.0,
            _ => [1.0, 2.0, 0.5, 1.5][usize::from(b.next() % 4)],
        })
        .collect();
    let spans = n - degree;
    let knots: Vec<f64> = (0..=spans).map(|k| k as f64).collect();
    let mut mults = vec![1; spans + 1];
    mults[0] = degree + 1;
    mults[spans] = degree + 1;
    if quarter && spans > 1 {
        mults[1] = degree;
    }
    if mults.iter().sum::<usize>() != n + degree + 1 {
        return;
    }
    let Ok(c) = BSplineCurve3::new(degree, poles, Some(weights), knots, mults) else {
        return;
    };
    let found = match spline_curve(&c, &conic) {
        Ok(r) => r,
        Err(Error::ComputationLimit(_)) => return,
        Err(e) => panic!("unexpected error {e}"),
    };
    let scale = |p: Point3| (p - Point3::ORIGIN).length().max(1.0) * 8.0;
    for [lo, hi] in &found.overlaps {
        assert!(lo < hi, "an overlap");
        for k in 0..=4 {
            let p = c.point(lo + (hi - lo) * f64::from(k) / 4.0).unwrap();
            let r = residual(&conic, p);
            assert!(r <= 1e-9 * scale(p), "{r}: overlap off {conic:?}");
        }
    }
    for w in found.points.windows(2) {
        assert!(mid(w[0].parameters[0]) <= mid(w[1].parameters[0]), "sorted");
    }
    for p in &found.points {
        let x = Point3::new(mid(p.point[0]), mid(p.point[1]), mid(p.point[2]));
        let on = c.point(mid(p.parameters[0])).unwrap();
        assert!(
            (on - x).length() <= 1e-9 * scale(x),
            "{p:?} is not the spline's point"
        );
        let on = at(&conic, mid(p.parameters[1]));
        assert!(
            (on - x).length() <= 1e-9 * scale(x),
            "{p:?} is not the conic's point"
        );
    }
}
