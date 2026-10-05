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
//! S8b.3: a leading byte in `192..224` takes a spline profile instead (a
//! rectangle with a quadratic bulge, a wave of three quadratic spans, a
//! square with a lens hole of two cubics given either way round), with the
//! tangent modes touching the spline at its apex. S8e: a leading byte in
//! `160..192` takes a sheet or a closed wire of any of those profiles (a
//! wire the outer boundary), split by a plane through a dyadic point,
//! through a vertex, along an edge, tangent, parallel to the body or in its
//! plane, or at an angle: the measures (areas or lengths) add up, each piece
//! lies on its side and moves rigidly with its ids. S9e.4b.3b: each oblique
//! piece of a prism with arcs and each piece of a cone, zone or cap by a
//! plane not normal to its axis is given to Booleans (`PIECE_BOOLEANS`): its
//! common with a box holding it is itself, and the first piece's fuse, cut
//! and common with a box turned about its centroid obey the pair identities
//! when all three evaluate; a refusal names a later step's or S9's rule,
//! never S9e.4's. S9e.4b.3c.1: a zone's or cap's first piece also meets a
//! whole ball of its sphere in a turned frame: their common is the piece,
//! their fuse the ball.
use crate::analytic_intersections::Bytes;
use rusty_occt::identity::OperationId;
use rusty_occt::topology::SplineSpan;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Side, Solid,
    Tolerance, Vec3,
};

/// Whether split pieces with curved faces are given to Booleans (S9e.4b.3b:
/// the given model of the split's primitive common its planes'
/// half-spaces): the prisms' oblique pieces and the cones', zones' and
/// caps' pieces by planes not normal to their axes. On; a torus's pieces
/// are not given (their arrangements take seconds, too long under the
/// sanitizer).
const PIECE_BOOLEANS: bool = true;

/// S9e.4b.3b: a split piece with a curved face given to Booleans: its common
/// with a box holding it is itself; with a box turned about its centroid
/// (`turned`: the first piece only, for time), the pair identities where
/// fuse, cut and common evaluate. A refusal is
/// documented (`Degenerate`, `ComputationLimit`, or `OutOfDomain` of a later
/// step), never S9e.4's refusal of a plane's piece.
fn piece_booleans(piece: &Solid, turned: bool) {
    use rusty_occt::history::History;
    let curved = piece
        .topology()
        .faces()
        .iter()
        .any(|f| !matches!(f.surface, rusty_occt::topology::Surface::Plane(_)));
    if !PIECE_BOOLEANS || !curved {
        return;
    }
    let tolerance = Tolerance::default();
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> Option<f64> {
        match r {
            Ok((out, _)) => Some(out.iter().map(|s| s.mass_properties().volume).sum()),
            Err(Error::Degenerate(_) | Error::ComputationLimit(_)) => None,
            Err(Error::OutOfDomain(m)) => {
                assert!(!m.ends_with("(S9e.4)"), "a split piece refused: {m}");
                None
            }
            Err(e) => panic!("unexpected error {e}"),
        }
    };
    let own = piece.mass_properties().volume;
    let bounds = piece.bounds();
    let reach = (bounds.max - bounds.min).length().max(1.0);
    let corner = bounds.min + Vec3::new(-reach, -reach, -reach);
    let Ok(holder) = Solid::box_at(
        corner,
        Vec3::new(3.0 * reach, 3.0 * reach, 3.0 * reach),
        tolerance,
    ) else {
        return;
    };
    if let Some(v) = volume(piece.common(OperationId(4), &holder)) {
        assert!((v - own).abs() <= 1e-9 * own.max(1.0), "{v} for {own}");
    }
    if !turned {
        return;
    }
    // S9e.4b.3c.1: a zone's or cap's piece and a whole ball of its sphere in
    // a turned frame (faces of both inputs on one sphere): their common is
    // the piece, their fuse the ball.
    let sphere = piece
        .topology()
        .faces()
        .iter()
        .find_map(|f| match f.surface {
            rusty_occt::topology::Surface::Sphere { frame, radius } => {
                Some((frame.origin(), radius))
            }
            _ => None,
        });
    if let Some((centre, radius)) = sphere {
        let half = std::f64::consts::FRAC_PI_2;
        if let Ok((ball, _)) = Frame3::new(
            centre,
            Vec3::new(4.0, 1.0, 8.0),
            Vec3::new(-7.0, -4.0, 4.0),
            tolerance,
        )
        .and_then(|f| Solid::sphere_with(OperationId(7), f, radius, -half, half, tolerance))
        {
            let whole = ball.mass_properties().volume;
            if let Some(v) = volume(piece.common(OperationId(8), &ball)) {
                assert!((v - own).abs() <= 1e-9 * whole, "{v} for {own}");
            }
            if let Some(v) = volume(piece.fuse(OperationId(8), &ball)) {
                assert!((v - whole).abs() <= 1e-9 * whole, "{v} for {whole}");
            }
        }
    }
    let c = piece.mass_properties().centroid;
    let Ok(turned) = Frame3::new(
        c,
        Vec3::new(1.0, 2.0, 2.0),
        Vec3::new(2.0, 1.0, -2.0),
        tolerance,
    ) else {
        return;
    };
    let side = 0.5 * reach;
    let square =
        [(0.0, 0.0), (side, 0.0), (side, side), (0.0, side)].map(|(x, y)| Point2::new(x, y));
    let Ok(profile) = Boundary::polygon(square.to_vec(), tolerance)
        .and_then(|b| Profile::new(b, vec![], tolerance))
    else {
        return;
    };
    let Ok((other, _)) = Solid::extrude_with(OperationId(5), profile, turned, 0.0, side) else {
        return;
    };
    let vb = other.mass_properties().volume;
    let (f, k, m) = (
        volume(piece.fuse(OperationId(6), &other)),
        volume(piece.cut(OperationId(6), &other)),
        volume(piece.common(OperationId(6), &other)),
    );
    if let (Some(f), Some(k), Some(m)) = (f, k, m) {
        let size = own + vb;
        assert!(
            (f + m - own - vb).abs() <= 1e-9 * size,
            "{f} + {m} for {own} + {vb}"
        );
        assert!((k - (own - m)).abs() <= 1e-9 * size, "{k} for {own} - {m}");
    }
}

/// A nonrational spline segment over its whole domain.
fn spline(
    degree: usize,
    poles: &[(f64, f64)],
    knots: Vec<f64>,
    mults: Vec<usize>,
) -> Option<Segment> {
    let poles = poles.iter().map(|(x, y)| Point2::new(*x, *y)).collect();
    let curve = BSplineCurve2::new(degree, poles, None, knots, mults).ok()?;
    Some(Segment::Spline(SplineSpan::whole(curve)))
}

/// S8b.3's spline profiles on dyadic sizes.
pub(crate) fn spline_profile(kind: u8, s: f64, t: f64) -> Option<Profile> {
    let tol = Tolerance::default();
    let (outer, holes) = match kind % 3 {
        // A rectangle whose right side bulges to x = s + t / 2.
        0 => (
            Boundary::path(
                vec![
                    Point2::new(0.0, -t),
                    Point2::new(s, -t),
                    Point2::new(s, t),
                    Point2::new(0.0, t),
                ],
                vec![
                    Segment::Line,
                    spline(
                        2,
                        &[(s, -t), (s + t, 0.0), (s, t)],
                        vec![0.0, 1.0],
                        vec![3, 3],
                    )?,
                    Segment::Line,
                    Segment::Line,
                ],
                tol,
            )
            .ok()?,
            vec![],
        ),
        // A rectangle under a wave of three quadratic spans.
        1 => (
            Boundary::path(
                vec![
                    Point2::new(0.0, 0.0),
                    Point2::new(3.0 * s, 0.0),
                    Point2::new(3.0 * s, 2.0 * t),
                    Point2::new(0.0, 2.0 * t),
                ],
                vec![
                    Segment::Line,
                    Segment::Line,
                    spline(
                        2,
                        &[
                            (3.0 * s, 2.0 * t),
                            (2.5 * s, 3.0 * t),
                            (1.5 * s, t),
                            (0.5 * s, 3.0 * t),
                            (0.0, 2.0 * t),
                        ],
                        vec![0.0, 1.0, 2.0, 3.0],
                        vec![3, 1, 1, 3],
                    )?,
                    Segment::Line,
                ],
                tol,
            )
            .ok()?,
            vec![],
        ),
        // A square with a lens hole reaching y = +-3s/8, given
        // counter-clockwise or clockwise.
        _ => {
            let (a, h) = (s / 2.0, s / 2.0);
            let lower = [(-a, 0.0), (-a / 2.0, -h), (a / 2.0, -h), (a, 0.0)];
            let upper = [(a, 0.0), (a / 2.0, h), (-a / 2.0, h), (-a, 0.0)];
            let cubic = |p: &[(f64, f64)]| spline(3, p, vec![0.0, 1.0], vec![4, 4]);
            let hole = if t > 2.0 {
                Boundary::path(
                    vec![Point2::new(-a, 0.0), Point2::new(a, 0.0)],
                    vec![cubic(&lower)?, cubic(&upper)?],
                    tol,
                )
            } else {
                let rev = |p: [(f64, f64); 4]| [p[3], p[2], p[1], p[0]];
                Boundary::path(
                    vec![Point2::new(-a, 0.0), Point2::new(a, 0.0)],
                    vec![cubic(&rev(upper))?, cubic(&rev(lower))?],
                    tol,
                )
            }
            .ok()?;
            (
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
                vec![hole],
            )
        }
    };
    Profile::new(outer, holes, tol).ok()
}

pub(crate) fn profile(kind: u8, s: f64, t: f64) -> Option<Profile> {
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
    if data.first().is_some_and(|k| *k >= 192) {
        return check_prism(&data[1..], true);
    }
    if data.first().is_some_and(|k| *k >= 160) {
        return check_body(&data[1..]);
    }
    check_prism(data, false);
}

/// S8e: a sheet or a closed wire split by a plane.
fn check_body(data: &[u8]) {
    use rusty_occt::Body;
    let mut b = Bytes(data, 0);
    let (kind, mode, flags) = (b.next(), b.next() % 8, b.next());
    let (splined, wire, tilted) = (flags & 1 == 1, flags & 2 == 2, flags & 4 == 4);
    let s = 1.0 + f64::from(b.next() % 32) / 4.0;
    let t = 0.5 + f64::from(b.next() % 32) / 8.0;
    let made = if splined {
        spline_profile(kind, s, t)
    } else {
        profile(kind, s, t)
    };
    let Some(p) = made else { return };
    let o = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
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
    let made = if wire {
        Body::wire_from_boundary_with(OperationId(1), p.outer().clone(), frame, p.tolerance())
    } else {
        Body::face_from_profile_with(OperationId(1), p.clone(), frame)
    };
    let Ok((body, _)) = made else { return };
    let at = |u: f64, v: f64, w: f64| frame.point(Point2::new(u, v), w);
    let dir = |u: f64, v: f64, w: f64| frame.x() * u + frame.y() * v + frame.normal() * w;
    let (a, c) = (b.small(), b.small());
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
        // Tangent to the stadium's arc, the hole, the bulge or the lens.
        3 => match (splined, kind % if splined { 3 } else { 6 }) {
            (false, 2) => (at(s + t, 0.0, 0.0), dir(1.0, 0.0, 0.0)),
            (false, 4) => (at(t.min(s * 0.75), 0.0, 0.0), dir(1.0, 0.0, 0.0)),
            (true, 0) => (at(s + t / 2.0, 0.0, 0.0), dir(1.0, 0.0, 0.0)),
            (true, 2) => (at(0.0, 3.0 * s / 8.0, 0.0), dir(0.0, 1.0, 0.0)),
            _ => return,
        },
        // Parallel to the body, off it or in its plane.
        4 => (
            at(0.0, 0.0, f64::from(b.next() % 3) - 1.0),
            dir(0.0, 0.0, 1.0),
        ),
        // Through a vertex at an angle.
        5 => (at(vertex.x, vertex.y, 0.0), dir(a, c, 1.0)),
        _ => (at(b.dyadic(), b.dyadic(), 0.0), dir(a, c, b.small())),
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
    let (pieces, _history) = match body.split_by_plane(OperationId(2), plane) {
        Ok(r) => r,
        Err(Error::ComputationLimit(_) | Error::Degenerate(_)) => return,
        Err(e) => panic!("unexpected error {e}"),
    };
    let measure = |x: &Body| {
        let m = x.measure().expect("a certified measure");
        (0.5 * m.measure[0] + 0.5 * m.measure[1], m.centre)
    };
    let total: f64 = pieces.iter().map(|(_, x)| measure(x).0).sum();
    let whole = measure(&body).0;
    assert!(
        (total - whole).abs() <= 1e-9 * whole.max(1.0),
        "{total} for {whole}"
    );
    if pieces.len() == 1 {
        assert_eq!(
            pieces[0].1.topology().body_id(),
            body.topology().body_id(),
            "unchanged"
        );
        return;
    }
    let motion =
        rusty_occt::RigidTransform::rotation(o, Vec3::new(1.0, 2.0, 2.0), 0.5).expect("a rotation");
    for (side, piece) in &pieces {
        let (_, centre) = measure(piece);
        let c = Point3::new(
            0.5 * centre[0][0] + 0.5 * centre[0][1],
            0.5 * centre[1][0] + 0.5 * centre[1][1],
            0.5 * centre[2][0] + 0.5 * centre[2][1],
        );
        let g = (c - plane.origin()).dot(plane.normal());
        match side {
            Side::Below => assert!(g < 0.0, "{g}: a piece below its plane"),
            Side::Above => assert!(g > 0.0, "{g}: a piece above its plane"),
        }
        let (moved, _) = piece
            .transform_with(OperationId(3), motion)
            .expect("a piece moves rigidly");
        let ids = |x: &Body| x.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
        assert_eq!(ids(piece), ids(&moved), "a moved piece keeps its ids");
        let (m0, m1) = (measure(piece).0, measure(&moved).0);
        assert!((m0 - m1).abs() <= 1e-9 * m0.max(1.0), "{m0} moved to {m1}");
    }
}

fn check_prism(data: &[u8], splined: bool) {
    let mut b = Bytes(data, 0);
    let (kind, mode) = (b.next(), b.next() % 9);
    let s = 1.0 + f64::from(b.next() % 32) / 4.0;
    let t = 0.5 + f64::from(b.next() % 32) / 8.0;
    let h = 0.5 + f64::from(b.next() % 16) / 4.0;
    let made = if splined {
        spline_profile(kind, s, t)
    } else {
        profile(kind, s, t)
    };
    let Some(p) = made else { return };
    // Where the tangent modes touch: the stadium's right arc or the round
    // hole, the bulge's apex or the lens's top.
    let touch = match (splined, kind % if splined { 3 } else { 6 }) {
        (false, 2) => Some((Point2::new(s + t, 0.0), (1.0, 0.0))),
        (false, 4) => Some((Point2::new(t.min(s * 0.75), 0.0), (1.0, 0.0))),
        (true, 0) => Some((Point2::new(s + t / 2.0, 0.0), (1.0, 0.0))),
        (true, 2) => Some((Point2::new(0.0, 3.0 * s / 8.0), (0.0, 1.0))),
        _ => None,
    };
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
        // Tangent to the stadium's right arc, the hole or a spline.
        3 => match touch {
            Some((q, (nx, ny))) => (at(q.x, q.y, 0.0), dir(nx, ny, 0.0)),
            None => return,
        },
        4 => (
            at(0.0, 0.0, f64::from(b.next() % 16) / 16.0 * h),
            dir(0.0, 0.0, 1.0),
        ),
        // Oblique through a profile vertex at a cap, or touching a cap's
        // circle or arc where the stadium's or the hole's is.
        7 => (at(vertex.x, vertex.y, 0.0), dir(a, c, 1.0)),
        // (at the low cap for the stadium, the high one for the hole, either
        // for a spline).
        8 => match touch {
            Some((q, (nx, ny)))
                if (!splined && kind % 6 == 2) || (splined && b.next() % 2 == 0) =>
            {
                (at(q.x, q.y, 0.0), dir(nx, ny, 1.0))
            }
            Some((q, (nx, ny))) => (at(q.x, q.y, h), dir(nx, ny, -1.0)),
            None => return,
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
        // S8a.2's oblique pieces (the other modes' pieces are prisms).
        if mode >= 6 {
            piece_booleans(piece, std::ptr::eq(piece, &pieces[0].1));
        }
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
    let (kind, mode) = (b.next() % 5, b.next() % 9);
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
    let (made, low_w, high_w) = match kind {
        0 => (
            Solid::cone_with(OperationId(1), frame, s, t, h, tolerance),
            0.0,
            h,
        ),
        1 => (
            Solid::sphere_with(OperationId(1), frame, s, -FRAC_PI_2, FRAC_PI_2, tolerance),
            -s,
            s,
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
            )
        }
        3 => (
            Solid::sphere_with(OperationId(1), frame, s, 0.25, FRAC_PI_2, tolerance),
            s * 0.25f64.sin(),
            s,
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
            )
        }
    };
    let Ok((solid, _)) = made else { return };
    let at = |w: f64| frame.point(rusty_occt::Point2::new(0.0, 0.0), w);
    // Every plane splits since S8d.3 (a torus's two caps aside).
    let spiric = false;
    let (point, normal, may_refuse) = match mode {
        // Normal to the axis at a dyadic height, at an end, or beyond.
        0 => (
            at(low_w + (high_w - low_w) * f64::from(b.next() % 17) / 16.0),
            frame.normal(),
            false,
        ),
        1 => (at(high_w), frame.normal() * -1.0, false),
        // Through a dyadic point, tilted.
        2 | 3 => {
            let p = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
            let n = Vec3::new(b.small(), b.small(), b.small() + 0.5);
            (o + (p - Point3::ORIGIN) * 0.25, n, spiric)
        }
        // Containing the axis (S8c.2): always split.
        4 => (
            o,
            frame.x() * b.small() + frame.y() * (b.small() + 0.5),
            false,
        ),
        5 => (
            at(0.5 * (low_w + high_w)),
            frame.x() + frame.normal(),
            spiric,
        ),
        // S8d.2: through a point of a cone's end circle, containing its
        // tangent there (touching the rim, or crossing it when tilted).
        6 => {
            let (w, r) = if b.next() % 2 == 0 { (0.0, s) } else { (h, t) };
            let (w, r) = if kind == 0 {
                (w, r)
            } else {
                (low_w, (s * s - low_w * low_w).max(0.0).sqrt())
            };
            (
                frame.point(rusty_occt::Point2::new(r, 0.0), w),
                frame.x() * b.small() + frame.normal() * (b.small() + 0.5),
                spiric,
            )
        }
        // Parallel to a cone's ruling at u = 0 (a parabola).
        7 => (
            at(low_w + (high_w - low_w) * f64::from(b.next() % 17) / 16.0)
                + frame.x() * (b.small() * 0.25),
            frame.x() * h + frame.normal() * (s - t),
            spiric,
        ),
        // Through a frustum's virtual apex (its rulings), or a sphere's
        // centre.
        _ => {
            let w = if kind == 0 && s != t {
                -s * h / (t - s)
            } else {
                0.0
            };
            (
                at(w),
                frame.x() + frame.y() * (b.small() * 0.25) + frame.normal() * (b.small() * 0.25),
                spiric,
            )
        }
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
        // A plane through an apex or pole off the axis (S8d.2's domain), a
        // torus cut in two caps (S8d.3's).
        Err(Error::OutOfDomain(m)) if m.contains("apex or pole") || m.contains("two caps") => {
            return
        }
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
        // S8c.2's and S8d.2's halves of a cone, zone or cap (a whole
        // sphere's pieces and those normal to the axis are caps, zones and
        // frusta; a torus's are slow).
        if matches!(kind, 0 | 2 | 3) && mode >= 2 {
            piece_booleans(piece, std::ptr::eq(piece, &pieces[0].1));
        }
    }
}
