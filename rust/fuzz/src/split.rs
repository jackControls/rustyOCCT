//! Prisms split by planes (S8 of REVIEW_NOTES.md): a rectangle, a regular
//! polygon, a stadium, a U or a rectangle with a round or square hole on
//! dyadic sizes, split by a plane chosen with an exact degeneracy on purpose:
//! parallel to the axis at a dyadic offset, through a profile vertex, along
//! a profile edge, tangent to an arc or the hole; normal to the axis at a
//! dyadic height or in a cap; or oblique (S8a.2), through a cap at a dyadic
//! point or touching a cap's edge. The split never panics and fails only as
//! documented (a certified comparison it cannot decide, a degenerate piece);
//! every piece validates as it is built and its history passes the
//! independent check (in the split, debug builds); the pieces' volumes add
//! up to the solid's, every piece lies on its side of the plane and moves
//! rigidly with its ids; a missing or touching plane returns the solid.
use crate::analytic_intersections::Bytes;
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Side, Solid, Tolerance, Vec3,
};

fn profile(kind: u8, s: f64, t: f64) -> Option<Profile> {
    let tol = Tolerance::default();
    let (outer, holes) = match kind % 6 {
        0 => (Boundary::rectangle(s, t, tol).ok()?, vec![]),
        1 => {
            let n = 3 + (s as usize) % 6;
            let pts = (0..n)
                .map(|k| {
                    let a = std::f64::consts::TAU * k as f64 / n as f64;
                    Point2::new(t * a.cos(), t * a.sin())
                })
                .collect();
            (Boundary::polygon(pts, tol).ok()?, vec![])
        }
        2 => (
            Boundary::path(
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
                tol,
            )
            .ok()?,
            vec![],
        ),
        3 => (
            Boundary::polygon(
                vec![
                    Point2::new(0.0, 0.0),
                    Point2::new(3.0 * s, 0.0),
                    Point2::new(3.0 * s, 2.0 * t),
                    Point2::new(2.0 * s, 2.0 * t),
                    Point2::new(2.0 * s, t),
                    Point2::new(s, t),
                    Point2::new(s, 2.0 * t),
                    Point2::new(0.0, 2.0 * t),
                ],
                tol,
            )
            .ok()?,
            vec![],
        ),
        4 => (
            Boundary::polygon(
                vec![
                    Point2::new(-s, -s),
                    Point2::new(s, -s),
                    Point2::new(s, s),
                    Point2::new(-s, s),
                ],
                tol,
            )
            .ok()?,
            vec![Boundary::circle(Point2::new(0.0, 0.0), (t).min(s * 0.75), tol).ok()?],
        ),
        _ => (
            Boundary::polygon(
                vec![
                    Point2::new(-s, -s),
                    Point2::new(s, -s),
                    Point2::new(s, s),
                    Point2::new(-s, s),
                ],
                tol,
            )
            .ok()?,
            vec![Boundary::polygon(
                vec![
                    Point2::new(-s / 2.0, -s / 2.0),
                    Point2::new(s / 2.0, -s / 2.0),
                    Point2::new(s / 2.0, s / 2.0),
                    Point2::new(-s / 2.0, s / 2.0),
                ],
                tol,
            )
            .ok()?],
        ),
    };
    Profile::new(outer, holes, tol).ok()
}

pub fn check_split(data: &[u8]) {
    if data.first().is_some_and(|k| *k >= 224) {
        return check_revolved(&data[1..]);
    }
    let mut b = Bytes(data, 0);
    let (kind, mode) = (b.next(), b.next() % 9);
    let s = 1.0 + f64::from(b.next() % 32) / 4.0;
    let t = 0.5 + f64::from(b.next() % 32) / 8.0;
    let h = 0.5 + f64::from(b.next() % 16) / 4.0;
    let Some(p) = profile(kind, s, t) else { return };
    let o = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
    let tilted = b.next() % 2 == 1;
    let frame = if !tilted {
        Frame3::xy()
    } else {
        match Frame3::new(
            o,
            Vec3::new(0.0, 3.0, 4.0),
            Vec3::new(1.0, 0.0, 0.0),
            Tolerance::default(),
        ) {
            Ok(f) => f,
            Err(_) => return,
        }
    };
    let Ok((solid, _)) = Solid::extrude_with(OperationId(1), p.clone(), frame, 0.0, h) else {
        return;
    };
    let at = |u: f64, v: f64, w: f64| frame.point(Point2::new(u, v), w);
    let dir = |u: f64, v: f64, w: f64| frame.x() * u + frame.y() * v + frame.normal() * w;
    let (a, c) = (b.small(), b.small());
    // A profile vertex and an edge for the degenerate modes.
    let vertex = match p
        .outer()
        .polygon_vertices()
        .or(p.outer().path_geometry().map(|g| g.0))
    {
        Some(v) if !v.is_empty() => v[b.next() as usize % v.len()],
        _ => Point2::new(0.0, 0.0),
    };
    let (point, normal) = match mode {
        0 => (at(b.dyadic(), b.dyadic(), 0.0), dir(a, c, 0.0)),
        1 => (at(vertex.x, vertex.y, 0.0), dir(a, c, 0.0)),
        2 => match p.outer().polygon_vertices() {
            Some(v) if v.len() > 1 => {
                let k = b.next() as usize % v.len();
                let (p0, p1) = (v[k], v[(k + 1) % v.len()]);
                (at(p0.x, p0.y, 0.0), dir(p0.y - p1.y, p1.x - p0.x, 0.0))
            }
            _ => return,
        },
        // Tangent to the stadium's right arc, or to the hole.
        3 => match kind % 6 {
            2 => (at(s + t, 0.0, 0.0), dir(1.0, 0.0, 0.0)),
            4 => (at((t).min(s * 0.75), 0.0, 0.0), dir(1.0, 0.0, 0.0)),
            _ => return,
        },
        4 => (
            at(0.0, 0.0, f64::from(b.next() % 16) / 16.0 * h),
            dir(0.0, 0.0, 1.0),
        ),
        // Oblique through a profile vertex at a cap, or touching a cap's
        // circle or arc where the stadium's or the hole's is.
        7 => (at(vertex.x, vertex.y, 0.0), dir(a, c, 1.0)),
        8 => match kind % 6 {
            2 => (at(s + t, 0.0, 0.0), dir(1.0, 0.0, 1.0)),
            4 => (at((t).min(s * 0.75), 0.0, h), dir(1.0, 0.0, -1.0)),
            _ => return,
        },
        5 => (at(0.0, 0.0, h), dir(0.0, 0.0, 1.0)),
        _ => (at(b.dyadic(), b.dyadic(), h / 2.0), dir(a, c, 1.0)),
    };
    if normal.length() == 0.0 {
        return;
    }
    let hint = if normal.x.abs() < 0.5 * normal.length() {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let Ok(plane) = Frame3::new(point, normal, hint, Tolerance::default()) else {
        return;
    };
    let (pieces, _history) = match solid.split_by_plane(OperationId(2), plane) {
        Ok(r) => r,
        Err(Error::ComputationLimit(_)) => return,
        Err(Error::Degenerate(_)) => return,
        Err(e) => panic!("unexpected error {e}"),
    };
    let total: f64 = pieces.iter().map(|(_, p)| p.mass_properties().volume).sum();
    let whole = solid.mass_properties().volume;
    assert!(
        (total - whole).abs() <= 1e-9 * whole.max(1.0),
        "{total} for {whole}"
    );
    if pieces.len() == 1 {
        assert_eq!(
            pieces[0].1.topology().body_id(),
            solid.topology().body_id(),
            "unchanged"
        );
        return;
    }
    let n = plane.normal();
    let q = plane.origin();
    let motion =
        rusty_occt::RigidTransform::rotation(o, Vec3::new(1.0, 2.0, 2.0), 0.5).expect("a rotation");
    for (side, piece) in &pieces {
        let c = piece.mass_properties().centroid;
        let g = (c - q).dot(n);
        match side {
            Side::Below => assert!(g < 0.0, "{g}: a piece below its plane"),
            Side::Above => assert!(g > 0.0, "{g}: a piece above its plane"),
        }
        let (moved, _) = piece
            .transform_with(OperationId(3), motion)
            .expect("a piece moves rigidly");
        let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
        assert_eq!(ids(piece), ids(&moved), "a moved piece keeps its ids");
        let (v0, v1) = (
            piece.mass_properties().volume,
            moved.mass_properties().volume,
        );
        assert!((v0 - v1).abs() <= 1e-9 * v0.max(1.0), "{v0} moved to {v1}");
    }
}

/// S8c and S8d: a cone, frustum, apex cone, sphere, zone or whole torus on
/// dyadic sizes, in an
/// axis-aligned or a tilted frame, split by a plane normal to its axis at a
/// dyadic height, through an apex, a pole or a cap, a whole sphere by any
/// plane through a dyadic point, a plane containing the axis, or one oblique
/// to it (a cone's conic or a zone's circle waits for S8d: `OutOfDomain`
/// then). The split never
/// panics and fails only as documented; volumes add up, pieces lie on their
/// sides and move rigidly with their ids.
fn check_revolved(data: &[u8]) {
    use std::f64::consts::FRAC_PI_2;
    let mut b = Bytes(data, 0);
    let (kind, mode) = (b.next() % 5, b.next() % 6);
    let s = 0.5 + f64::from(b.next() % 16) / 4.0;
    let t = f64::from(b.next() % 16) / 4.0;
    let h = 0.5 + f64::from(b.next() % 16) / 4.0;
    let o = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
    let tilted = b.next() % 2 == 1;
    let tolerance = Tolerance::default();
    let axis = if tilted {
        Vec3::new(0.0, 3.0, 4.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let Ok(frame) = Frame3::new(o, axis, Vec3::new(1.0, 0.0, 0.0), tolerance) else {
        return;
    };
    // A cone (t may be zero: an apex), a whole sphere, a zone or a cap.
    let (made, low_w, high_w, whole) = match kind {
        0 => (
            Solid::cone_with(OperationId(1), frame, s, t, h, tolerance),
            0.0,
            h,
            false,
        ),
        1 => (
            Solid::sphere_with(OperationId(1), frame, s, -FRAC_PI_2, FRAC_PI_2, tolerance),
            -s,
            s,
            true,
        ),
        2 => {
            let (lo, hi) = (
                -0.5 - f64::from(b.next() % 8) / 8.0,
                0.25 + f64::from(b.next() % 8) / 8.0,
            );
            (
                Solid::sphere_with(OperationId(1), frame, s, lo, hi, tolerance),
                s * lo.sin(),
                s * hi.sin(),
                false,
            )
        }
        3 => (
            Solid::sphere_with(OperationId(1), frame, s, 0.25, FRAC_PI_2, tolerance),
            s * 0.25f64.sin(),
            s,
            false,
        ),
        // S8d.1: a whole torus (its tube's radius from t, the gap to the
        // axis from s).
        _ => {
            let minor = 0.25 + t / 4.0;
            (
                Solid::torus_with(
                    OperationId(1),
                    frame,
                    minor + s,
                    minor,
                    -std::f64::consts::PI,
                    std::f64::consts::PI,
                    std::f64::consts::TAU,
                    tolerance,
                ),
                -minor,
                minor,
                false,
            )
        }
    };
    let Ok((solid, _)) = made else { return };
    let at = |w: f64| frame.point(rusty_occt::Point2::new(0.0, 0.0), w);
    let (point, normal, may_refuse) = match mode {
        // Normal to the axis at a dyadic height, at an end, or beyond.
        0 => (
            at(low_w + (high_w - low_w) * f64::from(b.next() % 17) / 16.0),
            frame.normal(),
            false,
        ),
        1 => (at(high_w), frame.normal() * -1.0, false),
        // Through a dyadic point, tilted: a whole sphere always splits.
        2 | 3 => {
            let p = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
            let n = Vec3::new(b.small(), b.small(), b.small() + 0.5);
            (o + (p - Point3::ORIGIN) * 0.25, n, !whole)
        }
        // Containing the axis (S8c.2): always split.
        4 => (
            o,
            frame.x() * b.small() + frame.y() * (b.small() + 0.5),
            false,
        ),
        _ => (
            at(0.5 * (low_w + high_w)),
            frame.x() + frame.normal(),
            !whole,
        ),
    };
    if normal.length() == 0.0 {
        return;
    }
    let hint = if normal.x.abs() < 0.5 * normal.length() {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let Ok(plane) = Frame3::new(point, normal, hint, tolerance) else {
        return;
    };
    let pieces = match solid.split_by_plane(OperationId(2), plane) {
        Ok((pieces, _)) => pieces,
        Err(Error::ComputationLimit(_) | Error::Degenerate(_)) => return,
        Err(Error::OutOfDomain(_)) if may_refuse => return,
        Err(e) => panic!("unexpected error {e}"),
    };
    let total: f64 = pieces.iter().map(|(_, p)| p.mass_properties().volume).sum();
    let whole_volume = solid.mass_properties().volume;
    assert!(
        (total - whole_volume).abs() <= 1e-9 * whole_volume.max(1.0),
        "{total} for {whole_volume}"
    );
    if pieces.len() == 1 {
        assert_eq!(pieces[0].1.topology().body_id(), solid.topology().body_id());
        return;
    }
    let motion =
        rusty_occt::RigidTransform::rotation(o, Vec3::new(1.0, 2.0, 2.0), 0.5).expect("a rotation");
    for (side, piece) in &pieces {
        let g = (piece.mass_properties().centroid - plane.origin()).dot(plane.normal());
        match side {
            Side::Below => assert!(g < 0.0, "{g}: a piece below its plane"),
            Side::Above => assert!(g > 0.0, "{g}: a piece above its plane"),
        }
        let (moved, _) = piece
            .transform_with(OperationId(3), motion)
            .expect("a piece moves rigidly");
        let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
        assert_eq!(ids(piece), ids(&moved), "a moved piece keeps its ids");
    }
}
