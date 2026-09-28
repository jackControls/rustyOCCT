//! S8b.1-2: profiles with spline segments and their prisms.
use rusty_occt::identity::OperationId;
use rusty_occt::topology::SplineSpan;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Profile, Segment, Solid, Tolerance,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn spline(poles: &[(f64, f64)], knots: Vec<f64>, mults: Vec<usize>) -> Segment {
    let degree = mults[0] - 1;
    let poles = poles.iter().map(|(x, y)| Point2::new(*x, *y)).collect();
    let curve = BSplineCurve2::new(degree, poles, None, knots, mults).expect("a spline");
    Segment::Spline(SplineSpan::whole(curve))
}

/// A 4 x 2 rectangle whose top side bulges as a quadratic spline through
/// the control point (2, 3): its area is 8 plus two thirds of the control
/// triangle's (4/3).
fn bulge() -> Boundary {
    let points = vec![
        Point2::new(0.0, 0.0),
        Point2::new(4.0, 0.0),
        Point2::new(4.0, 2.0),
        Point2::new(0.0, 2.0),
    ];
    let segments = vec![
        Segment::Line,
        Segment::Line,
        spline(
            &[(4.0, 2.0), (2.0, 3.0), (0.0, 2.0)],
            vec![0.0, 1.0],
            vec![3, 3],
        ),
        Segment::Line,
    ];
    Boundary::path(points, segments, tol()).expect("a valid spline path")
}

#[test]
fn a_spline_bulge_has_its_exact_area_and_centre() {
    let b = bulge();
    let area = 8.0 + 4.0 / 3.0;
    assert!((b.area() - area).abs() <= 1e-14 * area, "{}", b.area());
    // The bulge's centroid: the rectangle's (2, 1) and the parabolic
    // segment's two fifths up its height 1/2, weighted by their areas.
    let cy = (8.0 * 1.0 + 4.0 / 3.0 * (2.0 + 0.2)) / area;
    let c = b.centroid();
    assert!(
        (c.x - 2.0).abs() <= 1e-13 && (c.y - cy).abs() <= 1e-13,
        "{c:?}"
    );
}

#[test]
fn a_spline_prism_validates_with_its_certified_mass() {
    let profile = Profile::new(bulge(), vec![], tol()).expect("a profile");
    let (solid, _) = Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 3.0)
        .expect("a spline prism");
    let m = solid
        .topology()
        .mass_enclosure()
        .expect("certified mass properties");
    let volume = (8.0 + 4.0 / 3.0) * 3.0;
    assert!(
        m.volume[0] <= volume * (1.0 + 1e-14) && volume * (1.0 - 1e-14) <= m.volume[1],
        "{:?}",
        m.volume
    );
    assert!(
        m.volume[1] - m.volume[0] <= 1e-10 * volume,
        "{:?}",
        m.volume
    );
    // Written to .brep and read back as one certified solid.
    let text = rusty_occt::occt_brep::write(solid.topology(), 1e-7).expect("writable");
    let doc = rusty_occt::occt_brep::read(&text).expect("readable");
    let back = rusty_occt::occt_brep::import(&doc);
    assert_eq!(back.solids.len(), 1);
    assert!(
        back.solids[0].result.is_ok(),
        "{:?}",
        back.solids[0].result.as_ref().err()
    );
}

#[test]
fn a_spline_turning_a_quarter_turn_is_out_of_domain() {
    // A quadratic whose legs turn through 90 degrees.
    let points = vec![
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(0.0, 2.0),
    ];
    let segments = vec![
        Segment::Line,
        spline(
            &[(2.0, 0.0), (2.0, 2.0), (0.0, 2.0)],
            vec![0.0, 1.0],
            vec![3, 3],
        ),
        Segment::Line,
    ];
    assert!(matches!(
        Boundary::path(points, segments, tol()),
        Err(Error::OutOfDomain(_))
    ));
}

#[test]
fn a_spline_crossing_its_path_is_refused() {
    // The bulge pushed below the rectangle's bottom side: it crosses it.
    let points = vec![
        Point2::new(0.0, 0.0),
        Point2::new(4.0, 0.0),
        Point2::new(4.0, 2.0),
        Point2::new(0.0, 2.0),
    ];
    let segments = vec![
        Segment::Line,
        Segment::Line,
        spline(
            &[(4.0, 2.0), (3.0, -3.0), (1.0, -3.0), (0.0, 2.0)],
            vec![0.0, 1.0],
            vec![4, 4],
        ),
        Segment::Line,
    ];
    assert!(Boundary::path(points, segments, tol()).is_err());
}
