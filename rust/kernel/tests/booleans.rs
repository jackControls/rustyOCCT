//! S9a: Booleans of prisms in one frame.
use rusty_occt::history::{self, History};
use rusty_occt::identity::OperationId;
use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Boundary {
    Boundary::polygon(
        vec![
            Point2::new(x0, y0),
            Point2::new(x1, y0),
            Point2::new(x1, y1),
            Point2::new(x0, y1),
        ],
        tol(),
    )
    .unwrap()
}

fn prism(outer: Boundary, holes: Vec<Boundary>, frame: Frame3, lo: f64, hi: f64, op: u64) -> Solid {
    let profile = Profile::new(outer, holes, tol()).unwrap();
    Solid::extrude_with(OperationId(op), profile, frame, lo, hi)
        .unwrap()
        .0
}

fn check(a: &Solid, b: &Solid, out: &[Solid], h: &History) {
    let ins = [
        a.topology().entity_set(a.resolution()),
        b.topology().entity_set(b.resolution()),
    ];
    let outs: Vec<_> = out
        .iter()
        .map(|s| s.topology().entity_set(s.resolution()))
        .collect();
    let issues = history::check(&ins, &outs, h);
    assert!(issues.is_empty(), "{issues:?}");
}

fn volume(out: &[Solid]) -> f64 {
    out.iter().map(|s| s.mass_properties().volume).sum()
}

#[test]
fn overlapping_boxes() {
    let a = prism(rect(0.0, 0.0, 4.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(rect(2.0, 1.0, 6.0, 3.0), vec![], Frame3::xy(), 0.0, 1.0, 2);
    let (f, h) = a.fuse(OperationId(3), &b).unwrap();
    assert_eq!(f.len(), 1);
    assert!((volume(&f) - 14.0).abs() < 1e-12, "{}", volume(&f));
    check(&a, &b, &f, &h);
    let (c, h) = a.cut(OperationId(4), &b).unwrap();
    assert_eq!(c.len(), 1);
    assert!((volume(&c) - 6.0).abs() < 1e-12, "{}", volume(&c));
    check(&a, &b, &c, &h);
    let (m, h) = a.common(OperationId(5), &b).unwrap();
    assert_eq!(m.len(), 1);
    assert!((volume(&m) - 2.0).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
}

#[test]
fn a_box_and_a_cylinder() {
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 2.0, 1);
    let circle = Boundary::circle(Point2::new(4.0, 2.0), 1.0, tol()).unwrap();
    let b = prism(circle, vec![], Frame3::xy(), -1.0, 3.0, 2);
    let half = std::f64::consts::PI / 2.0;
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!((volume(&m) - 2.0 * half).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
    let (c, h) = a.cut(OperationId(4), &b).unwrap();
    assert!(
        (volume(&c) - (32.0 - 2.0 * half)).abs() < 1e-12,
        "{}",
        volume(&c)
    );
    check(&a, &b, &c, &h);
}

#[test]
fn a_lens_of_two_cylinders() {
    let c1 = Boundary::circle(Point2::new(0.0, 0.0), 1.0, tol()).unwrap();
    let c2 = Boundary::circle(Point2::new(1.0, 0.0), 1.0, tol()).unwrap();
    let a = prism(c1, vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(c2, vec![], Frame3::xy(), 0.0, 1.0, 2);
    // Two unit circles one apart: the lens is 2π/3 - √3/2.
    let lens = 2.0 * std::f64::consts::PI / 3.0 - 3f64.sqrt() / 2.0;
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!((volume(&m) - lens).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
    let (f, h) = a.fuse(OperationId(4), &b).unwrap();
    assert!((volume(&f) - (2.0 * std::f64::consts::PI - lens)).abs() < 1e-12);
    check(&a, &b, &f, &h);
}

#[test]
fn a_hole_through_a_box() {
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(
        Boundary::circle(Point2::new(2.0, 2.0), 1.0, tol()).unwrap(),
        vec![],
        Frame3::xy(),
        -1.0,
        2.0,
        2,
    );
    let (c, h) = a.cut(OperationId(3), &b).unwrap();
    assert_eq!(c.len(), 1);
    assert_eq!(c[0].profile().unwrap().holes().len(), 1);
    assert!((volume(&c) - (16.0 - std::f64::consts::PI)).abs() < 1e-12);
    check(&a, &b, &c, &h);
}

#[test]
fn identical_and_disjoint_and_inside() {
    let a = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let same = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 2);
    let (f, h) = a.fuse(OperationId(3), &same).unwrap();
    assert!((volume(&f) - 4.0).abs() < 1e-12);
    check(&a, &same, &f, &h);
    let (c, h) = a.cut(OperationId(4), &same).unwrap();
    assert!(c.is_empty());
    check(&a, &same, &c, &h);
    let far = prism(rect(5.0, 0.0, 6.0, 1.0), vec![], Frame3::xy(), 0.0, 1.0, 5);
    let (f, h) = a.fuse(OperationId(6), &far).unwrap();
    assert_eq!(f.len(), 2);
    check(&a, &far, &f, &h);
    let (m, h) = a.common(OperationId(7), &far).unwrap();
    assert!(m.is_empty());
    check(&a, &far, &m, &h);
    let inner = prism(rect(0.5, 0.5, 1.5, 1.5), vec![], Frame3::xy(), 0.0, 1.0, 8);
    let (f, h) = a.fuse(OperationId(9), &inner).unwrap();
    assert!((volume(&f) - 4.0).abs() < 1e-12);
    check(&a, &inner, &f, &h);
    let (m, h) = a.common(OperationId(10), &inner).unwrap();
    assert!((volume(&m) - 1.0).abs() < 1e-12);
    check(&a, &inner, &m, &h);
}

#[test]
fn boxes_sharing_part_of_a_wall() {
    let a = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(rect(2.0, 1.0, 4.0, 3.0), vec![], Frame3::xy(), 0.0, 1.0, 2);
    let (f, h) = a.fuse(OperationId(3), &b).unwrap();
    assert_eq!(f.len(), 1);
    assert!((volume(&f) - 8.0).abs() < 1e-12);
    check(&a, &b, &f, &h);
}

#[test]
fn frames_with_equal_axes() {
    // An offset origin whose coordinates are binary64 in the axes.
    let other = Frame3::new(
        Point3::new(1.0, 0.5, 0.25),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    let a = prism(rect(0.0, 0.0, 4.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], other, -0.25, 0.75, 2);
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!((volume(&m) - 2.0 * 1.5).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
    // The tilted frame, one origin.
    let tilted = Frame3::new(
        Point3::new(1.0, -2.0, 0.5),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    let a = prism(rect(0.0, 0.0, 4.0, 2.0), vec![], tilted, 0.0, 1.0, 4);
    let b = prism(rect(1.0, 0.5, 3.0, 2.5), vec![], tilted, 0.0, 1.0, 5);
    let (f, h) = a.fuse(OperationId(6), &b).unwrap();
    assert!(
        (volume(&f) - (8.0 + 4.0 - 3.0)).abs() < 1e-9,
        "{}",
        volume(&f)
    );
    check(&a, &b, &f, &h);
    // An offset that rounds in the tilted axes: S9b's.
    let shifted = Frame3::new(
        tilted.point(Point2::new(1.0, 0.5), 0.25),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    let c = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], shifted, 0.0, 1.0, 7);
    assert!(matches!(
        a.common(OperationId(8), &c),
        Err(rusty_occt::Error::OutOfDomain(_))
    ));
}

#[test]
fn an_arc_profile_cut_by_a_box() {
    let stadium = Boundary::path(
        vec![
            Point2::new(0.0, -1.0),
            Point2::new(3.0, -1.0),
            Point2::new(3.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        vec![
            Segment::Line,
            Segment::Arc {
                center: Point2::new(3.0, 0.0),
                radius: 1.0,
                ccw: true,
            },
            Segment::Line,
            Segment::Arc {
                center: Point2::new(0.0, 0.0),
                radius: 1.0,
                ccw: true,
            },
        ],
        tol(),
    )
    .unwrap();
    let a = prism(stadium, vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(
        rect(1.0, -2.0, 2.0, 2.0),
        vec![],
        Frame3::xy(),
        -1.0,
        2.0,
        2,
    );
    let (c, h) = a.cut(OperationId(3), &b).unwrap();
    assert_eq!(c.len(), 2);
    let area = 6.0 + std::f64::consts::PI;
    assert!((volume(&c) - (area - 2.0)).abs() < 1e-12, "{}", volume(&c));
    check(&a, &b, &c, &h);
}

fn circle(x: f64, y: f64, r: f64) -> Boundary {
    Boundary::circle(Point2::new(x, y), r, tol()).unwrap()
}

#[test]
fn tangent_cylinders() {
    use rusty_occt::Error;
    // Outside each other, touching at (1, 0): a fuse would touch itself
    // along a line, the common is empty.
    let a = prism(circle(0.0, 0.0, 1.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(circle(2.0, 0.0, 1.0), vec![], Frame3::xy(), 0.0, 1.0, 2);
    assert!(matches!(
        a.fuse(OperationId(3), &b),
        Err(Error::Degenerate(_))
    ));
    let (m, h) = a.common(OperationId(4), &b).unwrap();
    assert!(m.is_empty());
    check(&a, &b, &m, &h);
    let (c, h) = a.cut(OperationId(5), &b).unwrap();
    assert!((volume(&c) - std::f64::consts::PI).abs() < 1e-12);
    check(&a, &b, &c, &h);
    // Inside, touching at (2, 0): the fuse is the big one, the common the
    // small one, the cut would touch itself.
    let big = prism(circle(0.0, 0.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 6);
    let small = prism(circle(1.0, 0.0, 1.0), vec![], Frame3::xy(), 0.0, 1.0, 7);
    let (f, h) = big.fuse(OperationId(8), &small).unwrap();
    assert!((volume(&f) - 4.0 * std::f64::consts::PI).abs() < 1e-12);
    check(&big, &small, &f, &h);
    let (m, h) = big.common(OperationId(9), &small).unwrap();
    assert!((volume(&m) - std::f64::consts::PI).abs() < 1e-12);
    check(&big, &small, &m, &h);
    assert!(matches!(
        big.cut(OperationId(10), &small),
        Err(Error::Degenerate(_))
    ));
}

#[test]
fn arcs_of_one_circle() {
    // Two half-discs of the unit circle overlapping in a quarter.
    let half = |a0: f64, op: u64| {
        let (c, s) = (a0.cos(), a0.sin());
        let p = Point2::new(c, s);
        let q = Point2::new(-c, -s);
        let boundary = Boundary::path(
            vec![p, q],
            vec![
                Segment::Arc {
                    center: Point2::new(0.0, 0.0),
                    radius: 1.0,
                    ccw: true,
                },
                Segment::Line,
            ],
            tol(),
        )
        .unwrap();
        prism(boundary, vec![], Frame3::xy(), 0.0, 1.0, op)
    };
    let a = half(0.0, 1);
    let b = half(std::f64::consts::FRAC_PI_2, 2);
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!(
        (volume(&m) - std::f64::consts::FRAC_PI_4).abs() < 1e-9,
        "{}",
        volume(&m)
    );
    check(&a, &b, &m, &h);
    let (f, h) = a.fuse(OperationId(4), &b).unwrap();
    assert!(
        (volume(&f) - 3.0 * std::f64::consts::FRAC_PI_4).abs() < 1e-9,
        "{}",
        volume(&f)
    );
    check(&a, &b, &f, &h);
}

#[test]
fn a_vertex_on_the_others_edge() {
    // A triangle whose apex lies on the square's top edge.
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let tri = Boundary::polygon(
        vec![
            Point2::new(1.0, 6.0),
            Point2::new(3.0, 6.0),
            Point2::new(2.0, 4.0),
        ],
        tol(),
    )
    .unwrap();
    let b = prism(tri, vec![], Frame3::xy(), 0.0, 1.0, 2);
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!(m.is_empty());
    check(&a, &b, &m, &h);
    let tri = Boundary::polygon(
        vec![
            Point2::new(1.0, 2.0),
            Point2::new(3.0, 2.0),
            Point2::new(2.0, 4.0),
        ],
        tol(),
    )
    .unwrap();
    let inner = prism(tri, vec![], Frame3::xy(), 0.0, 1.0, 4);
    let (m, h) = a.common(OperationId(5), &inner).unwrap();
    assert!((volume(&m) - 2.0).abs() < 1e-12);
    check(&a, &inner, &m, &h);
}

#[test]
fn a_holed_box_and_a_tool() {
    let a = prism(
        rect(0.0, 0.0, 6.0, 6.0),
        vec![circle(3.0, 3.0, 1.0)],
        Frame3::xy(),
        0.0,
        1.0,
        1,
    );
    // Through the hole: the common is the box's part around it, in two.
    let b = prism(
        rect(2.5, -1.0, 3.5, 7.0),
        vec![],
        Frame3::xy(),
        -1.0,
        2.0,
        2,
    );
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert_eq!(m.len(), 2);
    let band = 0.75f64.sqrt() + std::f64::consts::PI / 3.0;
    let expected = 6.0 - band;
    assert!((volume(&m) - expected).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
    let (c, h) = a.cut(OperationId(4), &b).unwrap();
    assert_eq!(c.len(), 2);
    assert!((volume(&c) - (36.0 - std::f64::consts::PI - expected)).abs() < 1e-12);
    check(&a, &b, &c, &h);
}

#[test]
fn stacks_wait_for_s9a2() {
    use rusty_occt::Error;
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(rect(1.0, 1.0, 3.0, 3.0), vec![], Frame3::xy(), 0.5, 2.0, 2);
    assert!(matches!(
        a.fuse(OperationId(3), &b),
        Err(Error::OutOfDomain(_))
    ));
    assert!(matches!(
        a.cut(OperationId(4), &b),
        Err(Error::OutOfDomain(_))
    ));
    let (m, h) = a.common(OperationId(5), &b).unwrap();
    assert!((volume(&m) - 2.0).abs() < 1e-12);
    check(&a, &b, &m, &h);
}
