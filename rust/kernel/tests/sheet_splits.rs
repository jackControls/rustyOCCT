//! S8e: sheets and wires split by a plane.
use rusty_occt::history;
use rusty_occt::identity::OperationId;
use rusty_occt::{Body, Boundary, Frame3, Point2, Point3, Profile, Segment, Side, Tolerance, Vec3};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn plane(origin: [f64; 3], normal: [f64; 3]) -> Frame3 {
    let n = Vec3::new(normal[0], normal[1], normal[2]);
    let hint = if n.x.abs() < 0.5 * n.length() {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    Frame3::new(Point3::new(origin[0], origin[1], origin[2]), n, hint, tol()).unwrap()
}

fn check(body: &Body, pieces: &[(Side, Body)], h: &rusty_occt::history::History) {
    let before = [body.topology().entity_set(body.resolution())];
    let after: Vec<_> = pieces
        .iter()
        .map(|(_, p)| p.topology().entity_set(p.resolution()))
        .collect();
    let issues = history::check(&before, &after, h);
    assert!(issues.is_empty(), "{issues:?}");
}

fn measure(b: &Body) -> f64 {
    let m = b.measure().expect("a certified measure");
    0.5 * m.measure[0] + 0.5 * m.measure[1]
}

#[test]
fn a_square_sheet_splits_in_two() {
    let square = Boundary::rectangle(4.0, 2.0, tol()).unwrap();
    let profile = Profile::new(square, vec![], tol()).unwrap();
    let (sheet, _) = Body::face_from_profile_with(OperationId(1), profile, Frame3::xy()).unwrap();
    let (pieces, h) = sheet
        .split_by_plane(OperationId(2), plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.3]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    assert_eq!(pieces[0].0, Side::Below);
    let areas: Vec<f64> = pieces.iter().map(|(_, p)| measure(p)).collect();
    assert!(
        (areas[0] - 2.0).abs() < 1e-12 && (areas[1] - 6.0).abs() < 1e-12,
        "{areas:?}"
    );
    check(&sheet, &pieces, &h);
}

#[test]
fn a_plane_parallel_to_a_sheet_returns_it() {
    let square = Boundary::rectangle(4.0, 2.0, tol()).unwrap();
    let profile = Profile::new(square, vec![], tol()).unwrap();
    let (sheet, _) = Body::face_from_profile_with(OperationId(1), profile, Frame3::xy()).unwrap();
    for (z, side) in [(1.0, Side::Below), (-1.0, Side::Above), (0.0, Side::Below)] {
        let (pieces, h) = sheet
            .split_by_plane(OperationId(2), plane([0.0, 0.0, z], [0.0, 0.0, 1.0]))
            .unwrap();
        assert_eq!(pieces.len(), 1);
        assert_eq!(pieces[0].0, side);
        assert_eq!(pieces[0].1.topology().body_id(), sheet.topology().body_id());
        check(&sheet, &pieces, &h);
    }
}

#[test]
fn a_holed_sheet_splits_through_its_hole() {
    let outer = Boundary::rectangle(10.0, 10.0, tol()).unwrap();
    let hole = Boundary::circle(Point2::new(5.0, 5.0), 2.0, tol()).unwrap();
    let profile = Profile::new(outer, vec![hole], tol()).unwrap();
    let (sheet, _) = Body::face_from_profile_with(OperationId(1), profile, Frame3::xy()).unwrap();
    let (pieces, h) = sheet
        .split_by_plane(OperationId(2), plane([5.0, 0.0, 0.0], [1.0, 0.0, 0.0]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    let area = 100.0 - std::f64::consts::PI * 4.0;
    let total: f64 = pieces.iter().map(|(_, p)| measure(p)).sum();
    assert!((total - area).abs() < 1e-9, "{total}");
    check(&sheet, &pieces, &h);
}

#[test]
fn a_rectangle_wire_splits_into_two_open_wires() {
    let square = Boundary::rectangle(4.0, 2.0, tol()).unwrap();
    let (wire, _) =
        Body::wire_from_boundary_with(OperationId(1), square, Frame3::xy(), tol()).unwrap();
    let (pieces, h) = wire
        .split_by_plane(OperationId(2), plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    let lengths: Vec<f64> = pieces.iter().map(|(_, p)| measure(p)).collect();
    assert!((lengths[0] - 4.0).abs() < 1e-12, "{lengths:?}");
    assert!((lengths[1] - 8.0).abs() < 1e-12, "{lengths:?}");
    for (_, p) in &pieces {
        assert!(p.path().is_some());
    }
    check(&wire, &pieces, &h);
}

#[test]
fn a_wire_through_its_vertices_splits_them() {
    // A square cut along its diagonal through two stored vertices.
    let square = Boundary::rectangle(2.0, 2.0, tol()).unwrap();
    let (wire, _) =
        Body::wire_from_boundary_with(OperationId(1), square, Frame3::xy(), tol()).unwrap();
    let (pieces, h) = wire
        .split_by_plane(OperationId(2), plane([0.0, 0.0, 0.0], [1.0, -1.0, 0.0]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    for (_, p) in &pieces {
        assert!((measure(p) - 4.0).abs() < 1e-12);
    }
    check(&wire, &pieces, &h);
}

#[test]
fn a_circle_wire_splits_into_two_arcs() {
    let circle = Boundary::circle(Point2::new(0.0, 0.0), 1.0, tol()).unwrap();
    let (wire, _) =
        Body::wire_from_boundary_with(OperationId(1), circle, Frame3::xy(), tol()).unwrap();
    let (pieces, h) = wire
        .split_by_plane(OperationId(2), plane([0.5, 0.0, 0.0], [1.0, 0.0, 0.0]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    let total: f64 = pieces.iter().map(|(_, p)| measure(p)).sum();
    assert!((total - std::f64::consts::TAU).abs() < 1e-9, "{total}");
    check(&wire, &pieces, &h);
}

#[test]
fn a_spline_sheet_splits_through_its_spline() {
    use rusty_occt::topology::SplineSpan;
    let curve = rusty_occt::BSplineCurve2::new(
        2,
        vec![
            Point2::new(4.0, 2.0),
            Point2::new(2.0, 3.0),
            Point2::new(0.0, 2.0),
        ],
        None,
        vec![0.0, 1.0],
        vec![3, 3],
    )
    .unwrap();
    let boundary = Boundary::path(
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 2.0),
            Point2::new(0.0, 2.0),
        ],
        vec![
            Segment::Line,
            Segment::Line,
            Segment::Spline(SplineSpan::whole(curve)),
            Segment::Line,
        ],
        tol(),
    )
    .unwrap();
    let profile = Profile::new(boundary.clone(), vec![], tol()).unwrap();
    let (sheet, _) = Body::face_from_profile_with(OperationId(1), profile, Frame3::xy()).unwrap();
    let (pieces, h) = sheet
        .split_by_plane(OperationId(2), plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    let total: f64 = pieces.iter().map(|(_, p)| measure(p)).sum();
    assert!((total - (8.0 + 4.0 / 3.0)).abs() < 1e-9, "{total}");
    check(&sheet, &pieces, &h);
    let (wire, _) =
        Body::wire_from_boundary_with(OperationId(1), boundary, Frame3::xy(), tol()).unwrap();
    // The wire's length: 8 of lines and the parabola's
    // ∫ sqrt(16 + (2 - 4t)^2) dt over [0, 1] (mpmath, 40 digits).
    let length = 8.0 + 4.160_915_277_738_203;
    let m = wire.measure().expect("a certified length");
    assert!(
        m.measure[0] <= length && length <= m.measure[1],
        "{:?}",
        m.measure
    );
    assert!(
        m.measure[1] - m.measure[0] <= 1e-10 * length,
        "{:?}",
        m.measure
    );
    let (pieces, h) = wire
        .split_by_plane(OperationId(2), plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    check(&wire, &pieces, &h);
}

/// The split target's `crash-93175910`: a square wire `17.5` across in the
/// tilted frame split by a plane in its plane, its normal the frame's
/// normalized again (glibc's `hypot` turns it an ulp): the trace `a u + b v
/// + d = 0` of an ulp's tilt crossed the square, and each piece lay within
/// 1e-15 of the plane, its centre on its side only exactly, not by the
/// target's binary64 test. A plane within the resolution of the body's
/// plane over the whole body that would split it is now `Degenerate`, a
/// sheet's and a wire's alike, on every normal an ulp or two off the
/// frame's (a plane exactly parallel, or missing the body, returns it); a
/// plane tilted `1e-6` across the body still splits it.
#[test]
fn a_plane_within_the_resolution_of_the_bodys_plane_is_degenerate() {
    let s = 8.75;
    let square = || {
        Boundary::polygon(
            [(-s, -s), (s, -s), (s, s), (-s, s)]
                .map(|(x, y)| Point2::new(x, y))
                .to_vec(),
            tol(),
        )
        .unwrap()
    };
    let frame = Frame3::new(
        Point3::new(7.9375, 7.9375, 7.9375),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    let (wire, _) = Body::wire_from_boundary_with(OperationId(1), square(), frame, tol()).unwrap();
    let profile = Profile::new(square(), vec![], tol()).unwrap();
    let (sheet, _) = Body::face_from_profile_with(OperationId(1), profile, frame).unwrap();
    let n = frame.normal();
    let mut refused = 0;
    for body in [&wire, &sheet] {
        for k in -2i64..=2 {
            for (i, x) in [n.y, n.z].into_iter().enumerate() {
                let x = f64::from_bits((x.to_bits() as i64 + k) as u64);
                let m = if i == 0 {
                    Vec3::new(n.x, x, n.z)
                } else {
                    Vec3::new(n.x, n.y, x)
                };
                // Through the body's centre: a trace of any tilt crosses it.
                let at = Frame3::new(frame.origin(), m, frame.x(), tol()).unwrap();
                match body.split_by_plane(OperationId(2), at) {
                    Err(rusty_occt::Error::Degenerate(m)) => {
                        assert!(m.contains("body's plane"), "{m}");
                        refused += 1;
                    }
                    Ok((pieces, h)) => {
                        assert_eq!(pieces.len(), 1, "{k} {i}");
                        assert_eq!(pieces[0].1.topology().body_id(), body.topology().body_id());
                        check(body, &pieces, &h);
                    }
                    Err(e) => panic!("{e}"),
                }
                // 1e-3 above the body's plane (beyond the resolution): the
                // body below it, whole.
                let off = Frame3::new(frame.origin() + n * 1e-3, m, frame.x(), tol()).unwrap();
                let (pieces, _) = body.split_by_plane(OperationId(2), off).unwrap();
                assert_eq!(pieces.len(), 1);
                assert_eq!(pieces[0].0, Side::Below);
            }
        }
        // Tilted 1e-6 about the frame's y through its centre: the trace is
        // the y axis, the halves 8.75e-6 off the plane at the square's sides.
        let m = n + frame.x() * 1e-6;
        let at = Frame3::new(frame.origin(), m, frame.y(), tol()).unwrap();
        let (pieces, h) = body.split_by_plane(OperationId(2), at).unwrap();
        assert_eq!(pieces.len(), 2);
        assert_eq!(pieces[0].0, Side::Below);
        check(body, &pieces, &h);
        let total: f64 = pieces.iter().map(|(_, p)| measure(p)).sum();
        assert!((total - measure(body)).abs() <= 1e-9 * measure(body));
    }
    assert!(refused > 0);
}
