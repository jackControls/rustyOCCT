//! Prisms split by planes (S8 of REVIEW_NOTES.md): a rectangle, a regular
//! polygon, a stadium, a U or a rectangle with a round or square hole on
//! dyadic sizes, split by a plane chosen with an exact degeneracy on purpose:
//! parallel to the axis at a dyadic offset, through a profile vertex, along
//! a profile edge, tangent to an arc or the hole; normal to the axis at a
//! dyadic height or in a cap; or oblique (S8a.2). The split never panics
//! and fails only as documented (an oblique plane `OutOfDomain`, a
//! certified comparison it cannot decide, a piece within binary64 of a cap);
//! its history passes the independent check (in the split, debug builds);
//! the pieces' volumes add up to the solid's and every piece lies on its
//! side of the plane; a missing or touching plane returns the solid.
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
    let mut b = Bytes(data, 0);
    let (kind, mode) = (b.next(), b.next() % 7);
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
        // Oblique: mode 6, or any plane in the tilted frame (its stored
        // axes are not exactly orthogonal, so its "parallel" planes are not).
        Err(Error::OutOfDomain(_)) if mode == 6 || tilted => return,
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
    for (side, piece) in &pieces {
        let c = piece.mass_properties().centroid;
        let g = (c - q).dot(n);
        match side {
            Side::Below => assert!(g < 0.0, "{g}: a piece below its plane"),
            Side::Above => assert!(g > 0.0, "{g}: a piece above its plane"),
        }
    }
}
