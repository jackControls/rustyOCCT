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
fn a_spline_turning_through_a_half_turn_is_valid() {
    // A 4 x 2 rectangle whose left side is a cubic end bulging to x = -9/8
    // and turning through a half-turn: the end adds 9/5.
    let points = vec![
        Point2::new(0.0, 0.0),
        Point2::new(4.0, 0.0),
        Point2::new(4.0, 2.0),
        Point2::new(0.0, 2.0),
    ];
    let segments = vec![
        Segment::Line,
        Segment::Line,
        Segment::Line,
        spline(
            &[(0.0, 2.0), (-1.5, 2.0), (-1.5, 0.0), (0.0, 0.0)],
            vec![0.0, 1.0],
            vec![4, 4],
        ),
    ];
    let b = Boundary::path(points, segments, tol()).expect("a valid half-turn end");
    let area = 8.0 + 9.0 / 5.0;
    assert!((b.area() - area).abs() <= 1e-14 * area, "{}", b.area());
    let profile = Profile::new(b, vec![], tol()).expect("a profile");
    let (solid, _) = Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 1.0)
        .expect("a half-turn spline prism");
    let m = solid
        .topology()
        .mass_enclosure()
        .expect("certified mass properties");
    assert!(m.volume[0] <= area * (1.0 + 1e-14) && area * (1.0 - 1e-14) <= m.volume[1]);
}

#[test]
fn a_spline_with_a_cusp_is_refused() {
    // A cubic whose derivative vanishes at its middle, closed by a line:
    // its two halves double back at the cusp.
    let points = vec![Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)];
    let segments = vec![
        spline(
            &[(0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0)],
            vec![0.0, 1.0],
            vec![4, 4],
        ),
        Segment::Line,
    ];
    assert!(matches!(
        Boundary::path(points, segments, tol()),
        Err(Error::SelfIntersection)
    ));
}

#[test]
fn a_spline_looping_over_itself_is_refused() {
    // A cubic crossing itself before it returns to the line's far end.
    let points = vec![Point2::new(0.0, 0.0), Point2::new(4.0, 0.0)];
    let segments = vec![
        Segment::Line,
        spline(
            &[(4.0, 0.0), (-2.0, 4.0), (6.0, 4.0), (0.0, 0.0)],
            vec![0.0, 1.0],
            vec![4, 4],
        ),
    ];
    assert!(matches!(
        Boundary::path(points, segments, tol()),
        Err(Error::SelfIntersection)
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

#[test]
fn a_prism_with_a_clockwise_spline_hole_writes_and_reads_back() {
    // A 10 x 10 square with a lens hole of two cubics given clockwise: its
    // stored spans run against their curves (the prism's cap edges too).
    let outer = Boundary::rectangle(10.0, 10.0, tol()).expect("a square");
    let points = vec![Point2::new(3.0, 5.0), Point2::new(7.0, 5.0)];
    let segments = vec![
        spline(
            &[(3.0, 5.0), (4.0, 7.0), (6.0, 7.0), (7.0, 5.0)],
            vec![0.0, 1.0],
            vec![4, 4],
        ),
        spline(
            &[(7.0, 5.0), (6.0, 3.0), (4.0, 3.0), (3.0, 5.0)],
            vec![0.0, 1.0],
            vec![4, 4],
        ),
    ];
    let hole = Boundary::path(points, segments, tol()).expect("a lens");
    let profile = Profile::new(outer, vec![hole], tol()).expect("a profile");
    let (solid, _) = Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 2.0)
        .expect("a prism with a spline hole");
    let text = rusty_occt::occt_brep::write(solid.topology(), 1e-7).expect("writable");
    let doc = rusty_occt::occt_brep::read(&text).expect("readable");
    let back = rusty_occt::occt_brep::import(&doc);
    assert_eq!(back.solids.len(), 1);
    let read = back.solids[0].result.as_ref().expect("a certified solid");
    let (a, b) = (
        solid.mass_properties().volume,
        read.mass_enclosure().expect("certified mass").volume,
    );
    assert!(b[0] - 1e-9 * a <= a && a <= b[1] + 1e-9 * a, "{a} {b:?}");
}
