//! Lines and circles against analytic surfaces (S7c.1 of REVIEW_NOTES.md).
//! A plane, cylinder, cone, sphere or torus (the bytes of
//! `analytic_intersections`) and a line or a circle on dyadic data,
//! independent or with an exact degeneracy made on purpose: a line along the
//! surface's axis at its radius (a cylinder's ruling), or across it at that
//! offset (tangent to a cylinder or a sphere), or through the origin; a
//! circle coaxial with the surface at its radius, beside it at the sum of
//! the radii, or in a plane through the axis. The intersection never panics
//! and fails only by a certified comparison it cannot decide; points are
//! sorted with ordered enclosures, each on the curve at its parameter and on
//! the surface; a contained curve lies on the surface; an arc is its whole
//! circle; and an exact dyadic translation keeps the result wherever the
//! moved frames keep their stored axes.
use crate::analytic_intersections::{distance, frame, make, translated, Bytes};
use rusty_occt::intersection::{curve_surface, CurveSurfaceIntersection};
use rusty_occt::topology::{Curve3, Surface};
use rusty_occt::{Error, Point3, RigidTransform, Vec3};

fn mid([lo, hi]: [f64; 2]) -> f64 {
    0.5 * lo + 0.5 * hi
}

fn at(c: &Curve3, t: f64) -> Point3 {
    match c {
        Curve3::LineSegment { start, end } => *start + (*end - *start) * t,
        Curve3::Circle { frame, radius } => {
            frame.origin() + (frame.x() * t.cos() + frame.y() * t.sin()) * *radius
        }
        _ => unreachable!(),
    }
}

fn moved(c: &Curve3, t: RigidTransform) -> Curve3 {
    match c {
        Curve3::LineSegment { start, end } => Curve3::LineSegment {
            start: t.point(*start),
            end: t.point(*end),
        },
        Curve3::Circle { frame, radius } => Curve3::Circle {
            frame: frame
                .transformed(t, rusty_occt::Tolerance::default())
                .unwrap(),
            radius: *radius,
        },
        _ => unreachable!(),
    }
}

pub fn check_curve_surface(data: &[u8]) {
    let mut b = Bytes(data, 0);
    let (kind, shape, mode) = (b.next(), b.next(), b.next() % 4);
    let o = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
    let n = Vec3::new(b.small(), b.small(), b.small());
    let r = 0.25 + f64::from(b.next() % 32) / 8.0;
    let a = 0.05 + f64::from(b.next() % 64) / 48.0;
    let Some(f) = frame(o, n) else { return };
    let s = make(kind, f, r, a);
    // The surface's radius where a degeneracy is exact.
    let reach = match &s {
        Surface::Torus { major, minor, .. } => major + minor,
        Surface::Cone { radius, .. } => *radius,
        _ => r,
    };
    let (x, y, axis) = (f.x(), f.y(), f.normal());
    let curve = if shape % 2 == 0 {
        let (p0, d) = match mode {
            0 => (
                Point3::new(b.dyadic(), b.dyadic(), b.dyadic()),
                Vec3::new(b.small(), b.small(), b.small()),
            ),
            1 => (o + x * reach, axis),
            2 => (o + x * reach, y),
            _ => (o, Vec3::new(b.small(), b.small(), b.small())),
        };
        if d.length() == 0.0 {
            return;
        }
        Curve3::LineSegment {
            start: p0,
            end: p0 + d,
        }
    } else {
        let rc = 0.25 + f64::from(b.next() % 16) / 8.0;
        let (c, cn, rc) = match mode {
            0 => (
                Point3::new(b.dyadic(), b.dyadic(), b.dyadic()),
                Vec3::new(b.small(), b.small(), b.small()),
                rc,
            ),
            1 => (o + axis * b.dyadic(), axis, reach),
            2 => (o + x * (reach + rc), axis, rc),
            _ => (o, x, reach),
        };
        let Some(cf) = frame(c, cn) else { return };
        Curve3::Circle {
            frame: cf,
            radius: rc,
        }
    };
    let result = match curve_surface(&curve, &s) {
        Ok(r) => r,
        Err(Error::ComputationLimit(_)) => return,
        Err(e) => panic!("unexpected error {e}"),
    };
    let scale = |p: Point3| (p - Point3::ORIGIN).length().max(1.0) * 8.0;
    match &result {
        CurveSurfaceIntersection::Empty => {}
        CurveSurfaceIntersection::Contained => {
            for k in 0..6 {
                let p = at(&curve, -2.0 + 0.7 * f64::from(k));
                let gap = distance(&s, p - Point3::ORIGIN);
                assert!(
                    gap <= 1e-9 * scale(p),
                    "{gap}: contained {curve:?} off {s:?}"
                );
            }
        }
        CurveSurfaceIntersection::Points(points) => {
            assert!(!points.is_empty());
            for w in points.windows(2) {
                assert!(mid(w[0].parameter) <= mid(w[1].parameter), "sorted");
            }
            for p in points {
                assert!(p.parameter[0] <= p.parameter[1], "ordered enclosure");
                for [lo, hi] in p.point {
                    assert!(lo <= hi, "ordered enclosure");
                }
                let q = Point3::new(mid(p.point[0]), mid(p.point[1]), mid(p.point[2]));
                let on = at(&curve, mid(p.parameter));
                assert!(
                    (on - q).length() <= 1e-9 * scale(q),
                    "{p:?} is not the curve's point at its parameter"
                );
                let gap = distance(&s, q - Point3::ORIGIN);
                assert!(gap <= 1e-9 * scale(q), "{gap}: {p:?} off {s:?}");
            }
        }
    }
    if let Curve3::Circle { frame, radius } = &curve {
        let arc = Curve3::CircularArc {
            frame: *frame,
            radius: *radius,
            start_angle: 0.5,
            sweep_angle: 2.0,
        };
        assert_eq!(
            curve_surface(&arc, &s).as_ref().ok(),
            Some(&result),
            "an arc is its circle"
        );
    }
    // An exact dyadic translation (of dyadic points only, so it moves them
    // exactly): the same classes and contacts, the parameters overlapping
    // (lines keep theirs; circles their angles).
    let dyadic = |p: Point3| p.to_array().iter().all(|x| (x * 16.0).fract() == 0.0);
    let points = match &curve {
        Curve3::LineSegment { start, end } => vec![*start, *end],
        Curve3::Circle { frame, .. } => vec![frame.origin()],
        _ => unreachable!(),
    };
    if !points.into_iter().chain([o]).all(dyadic) {
        return;
    }
    let shift = RigidTransform::translation(Vec3::new(b.dyadic(), b.dyadic(), b.dyadic())).unwrap();
    let (curve2, s2) = (moved(&curve, shift), translated(&s, shift));
    // A moved frame is normalised again, which may change its stored axes
    // and so the exact sets: compare only frames that keep them.
    let axes = |f: rusty_occt::Frame3| [f.normal(), f.x(), f.y()].map(|v| v.to_array());
    let same = |a: &Curve3, b: &Curve3| match (a, b) {
        (Curve3::Circle { frame: f, .. }, Curve3::Circle { frame: g, .. }) => axes(*f) == axes(*g),
        _ => true,
    };
    if axes(crate::analytic_intersections::frame_of(&s))
        != axes(crate::analytic_intersections::frame_of(&s2))
        || !same(&curve, &curve2)
    {
        return;
    }
    let Ok(t) = curve_surface(&curve2, &s2) else {
        return;
    };
    match (&result, &t) {
        (CurveSurfaceIntersection::Points(p), CurveSurfaceIntersection::Points(q)) => {
            assert_eq!(p.len(), q.len(), "translated");
            for (p, q) in p.iter().zip(q) {
                assert_eq!(p.tangent, q.tangent, "translated");
                let ([a, b], [c, d]) = (p.parameter, q.parameter);
                assert!(a <= d && c <= b, "translated parameters overlap");
            }
        }
        _ => assert_eq!(
            std::mem::discriminant(&result),
            std::mem::discriminant(&t),
            "translated"
        ),
    }
}
