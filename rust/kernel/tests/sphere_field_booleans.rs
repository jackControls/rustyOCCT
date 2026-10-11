//! Lines of one quadratic field against quadrics meeting them at points of
//! another: the kernel's own stack of two prisms (a stadium plate fused
//! with a square ring over its arc end, the `boolean` fuzz target's `RIMS`
//! well on the world's axes as first drafted) cut by a ball, a rod and a
//! cone. The ring's side wall `y = 1/4` meets the plate's arc wall (radius
//! 3/8 about `(1/2, 0)`) along the vertical line `x = 1/2 + sqrt(5)/8`,
//! `y = 1/4`: a section's line through points of `Q(sqrt 5)`. A quadric
//! centred on `x = 1/2` (or with its axis along `x`) meets it in a rational
//! quadratic, its roots surds `sqrt E` of their own: the point
//! `p + t d` has its `x` in `Q(sqrt 5)` and its `z` in `Q(sqrt E)`, which
//! no point of one field holds, and a hemisphere's membership of it (a dot
//! product of its coordinates) panicked with "surds of one field". Such a
//! point is now written in the line's surd where `E / 5` is a rational
//! square, no meeting where it lies certainly outside the edge's and the
//! face's boxes, and refused (`ComputationLimit`, the line against that
//! quadric) elsewhere (`meet::line_point`).
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3,
};

fn tolerance() -> Tolerance {
    Tolerance::new(1e-7, 1e-12).unwrap()
}

fn frame(o: [f64; 3], n: [f64; 3], x: [f64; 3]) -> Frame3 {
    Frame3::new(
        Point3::new(o[0], o[1], o[2]),
        Vec3::new(n[0], n[1], n[2]),
        Vec3::new(x[0], x[1], x[2]),
        tolerance(),
    )
    .unwrap()
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Boundary {
    let points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        .map(|(x, y)| Point2::new(x, y))
        .to_vec();
    Boundary::polygon(points, tolerance()).unwrap()
}

/// The stadium plate of radius `3k` about `(0, 0)` and `(4k, 0)` over
/// `[0, h]` fused with the square `[2k, 9k] x [-2k, 2k]` less the square
/// `[7k/2, 11k/2] x [-k, k]` over `[h/2, 3h/2]`, on the world's axes, for
/// `k = 1/8`, `h = 7/4`.
fn well() -> Solid {
    let (k, h) = (0.125, 1.75);
    let tol = tolerance();
    let f = Frame3::xy();
    let r = 3.0 * k;
    let arc = |x: f64| Segment::Arc {
        center: Point2::new(x, 0.0),
        radius: r,
        ccw: true,
    };
    let stadium = Boundary::path(
        vec![
            Point2::new(0.0, -r),
            Point2::new(4.0 * k, -r),
            Point2::new(4.0 * k, r),
            Point2::new(0.0, r),
        ],
        vec![Segment::Line, arc(4.0 * k), Segment::Line, arc(0.0)],
        tol,
    )
    .unwrap();
    let prism = |profile: Profile, h0: f64, h1: f64, op: u64| {
        Solid::extrude_with(OperationId(op), profile, f, h0, h1)
            .unwrap()
            .0
    };
    let plate = prism(Profile::new(stadium, vec![], tol).unwrap(), 0.0, h, 31);
    let ring = prism(
        Profile::new(
            square(2.0 * k, -2.0 * k, 9.0 * k, 2.0 * k),
            vec![square(3.5 * k, -k, 5.5 * k, k)],
            tol,
        )
        .unwrap(),
        h / 2.0,
        1.5 * h,
        32,
    );
    let (out, _) = plate.fuse(OperationId(33), &ring).unwrap();
    assert_eq!(out.len(), 1);
    out.into_iter().next().unwrap()
}

/// A ball of radius `r` about `c`, on the world's `z`, its seam turned
/// `(3, 4, 0)`.
fn ball(c: [f64; 3], r: f64) -> Solid {
    let half = std::f64::consts::FRAC_PI_2;
    Solid::sphere_with(
        OperationId(7),
        frame(c, [0.0, 0.0, 1.0], [3.0, 4.0, 0.0]),
        r,
        -half,
        half,
        tolerance(),
    )
    .unwrap()
    .0
}

fn volume(r: rusty_occt::Result<(Vec<Solid>, rusty_occt::history::History)>) -> f64 {
    let (out, _) = r.unwrap();
    out.iter().map(|s| s.mass_properties().volume).sum()
}

/// Cut, common and fuse each evaluate, deterministically, with one volume:
/// the cut and the common the body's, the fuse the body's and the tool's
/// less their common.
fn evaluates_alike(body: &Solid, tool: &Solid) -> (f64, f64) {
    let vb = body.mass_properties().volume;
    let vt = tool.mass_properties().volume;
    let cut = volume(body.cut(OperationId(34), tool));
    let common = volume(body.common(OperationId(35), tool));
    let fuse = volume(body.fuse(OperationId(36), tool));
    assert!(cut > 0.0 && common > 0.0, "{cut} {common}");
    assert!(
        (cut + common - vb).abs() <= 1e-9 * vb,
        "{cut} + {common} for {vb}"
    );
    assert!(
        (fuse - (vb + vt - common)).abs() <= 1e-9 * (vb + vt),
        "{fuse} for {vb} + {vt} - {common}"
    );
    assert_eq!(cut, volume(body.cut(OperationId(34), tool)));
    (cut, common)
}

/// The line's points of two fields within the edge: each refused as the
/// line against its quadric (each panicked before).
#[test]
fn points_of_two_fields_are_refused() {
    let body = well();
    let refused = |tool: &Solid, text: &str| match body.cut(OperationId(34), tool) {
        Err(Error::ComputationLimit(t)) => assert_eq!(t, text),
        other => panic!("{:?}", other.map(|(o, _)| o.len())),
    };
    // A ball about (1/2, 1/4, 7/12): z = 7/12 +- sqrt(95)/8 on the line.
    refused(
        &ball([0.5, 0.25, 1.75 / 3.0], 1.25),
        "an irrational line against a sphere",
    );
    // A rod along x about y = 0, z = 21/16 of radius 1/2: z = 21/16 +-
    // sqrt(3)/4.
    let rod = Solid::cylinder_with(
        OperationId(7),
        frame([-1.0, 0.0, 1.3125], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        0.5,
        0.0,
        3.0,
        tolerance(),
    )
    .unwrap()
    .0;
    refused(&rod, "an irrational line against a cylinder");
    // A cone about the axis through (1/2, 1/8), meeting the line where its
    // radius is sqrt(6)/8.
    let cone = Solid::cone_with(
        OperationId(7),
        frame([0.5, 0.125, 0.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        0.1,
        0.6,
        1.5,
        tolerance(),
    )
    .unwrap()
    .0;
    refused(&cone, "an irrational line against a cone");
}

/// A ball whose roots on the line are surds of its own field: about
/// `(1/2, 7/16, 5/4)` of radius 7/16, `(z - 5/4)^2 = 5/64`, written in the
/// line's surd, each operation evaluating (it panicked before).
#[test]
fn points_of_the_line_s_field_evaluate() {
    let (cut, common) = evaluates_alike(&well(), &ball([0.5, 0.4375, 1.25], 0.4375));
    assert!(common < 0.2 && cut > 1.7, "{cut} {common}");
}

/// Balls meeting the line at points of two fields beyond its ends (above
/// and below the ring's wall): no meeting, each operation evaluating (they
/// panicked before).
#[test]
fn points_of_two_fields_beyond_the_edge_are_no_meeting() {
    for (cz, r) in [(1.3125, 1.0), (1.25, 1.125), (1.375, 1.0625)] {
        evaluates_alike(&well(), &ball([0.5, 0.25, cz], r));
    }
}
