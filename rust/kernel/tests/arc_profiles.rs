//! S5: profiles with circular-arc segments and their prisms. Ids and
//! histories are checked against the independent reference in `identity.rs`
//! and natively by `compare_history.py --family arc`; this file covers the
//! path's own contract: validity, orientation, labels, classification, mass
//! properties and the height split and stacked fuse.
#[path = "support/brep_protocol.rs"]
#[allow(dead_code)]
mod brep_protocol;
use rusty_occt::identity::{InputLabel, OperationId};
use rusty_occt::topology::Topology;
use rusty_occt::{
    Boundary, BoundaryLabels, Error, Frame3, Location, Point2, Point3, Profile, Segment, Solid,
    Tolerance, Vec3,
};

fn tol() -> Tolerance {
    Tolerance::new(1e-7, 1e-12).unwrap()
}

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn arc(cx: f64, cy: f64, r: f64, ccw: bool) -> Segment {
    Segment::Arc {
        center: p(cx, cy),
        radius: r,
        ccw,
    }
}

fn stadium() -> (Vec<Point2>, Vec<Segment>) {
    (
        vec![p(0., -1.), p(3., -1.), p(3., 1.), p(0., 1.)],
        vec![
            Segment::Line,
            arc(3., 0., 1., true),
            Segment::Line,
            arc(0., 0., 1., true),
        ],
    )
}

fn notch() -> (Vec<Point2>, Vec<Segment>) {
    (
        vec![
            p(0., 0.),
            p(10., 0.),
            p(10., 6.),
            p(6., 6.),
            p(4., 6.),
            p(0., 6.),
        ],
        vec![
            Segment::Line,
            Segment::Line,
            Segment::Line,
            arc(5., 6., 1., false),
            Segment::Line,
            Segment::Line,
        ],
    )
}

fn frame() -> Frame3 {
    Frame3::new(
        Point3::new(1.0, -2.0, 3.0),
        Vec3::new(0.3, -0.4, 0.8),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap()
}

fn prism(points: Vec<Point2>, segments: Vec<Segment>) -> Solid {
    let profile = Profile::new(
        Boundary::path(points, segments, tol()).unwrap(),
        vec![],
        tol(),
    )
    .unwrap();
    Solid::extrude_with(OperationId(5), profile, frame(), 0.0, 2.0)
        .unwrap()
        .0
}

#[test]
fn a_path_of_lines_is_the_polygon() {
    let points = vec![p(0., 0.), p(4., 0.), p(4., 3.), p(0., 3.)];
    let path = Boundary::path(points.clone(), vec![Segment::Line; 4], tol()).unwrap();
    assert_eq!(path, Boundary::polygon(points, tol()).unwrap());
    assert!(path.path_geometry().is_none());
}

#[test]
fn clockwise_input_is_stored_counter_clockwise_with_its_labels() {
    let (points, segments) = stadium();
    let n = points.len();
    // The same stadium entered clockwise from the same first point.
    let mut cw_points = vec![points[0]];
    cw_points.extend(points[1..].iter().rev());
    let cw_segments: Vec<Segment> = (0..n)
        .map(|j| match segments[n - 1 - j] {
            Segment::Arc { center, radius, .. } => Segment::Arc {
                center,
                radius,
                ccw: false,
            },
            s => s,
        })
        .collect();
    let label = |k: u64| InputLabel(k);
    let labels = BoundaryLabels {
        boundary: label(1),
        segments: (0..n as u64).map(|j| label(10 + j)).collect(),
        vertices: (0..n as u64).map(|j| label(20 + j)).collect(),
    };
    let ccw = Boundary::path(points, segments, tol()).unwrap();
    let cw = Boundary::path(cw_points, cw_segments, tol())
        .unwrap()
        .with_labels(labels)
        .unwrap();
    assert_eq!(cw.path_geometry(), ccw.path_geometry());
    let stored = cw.labels().unwrap();
    // Stored segment j is input segment n-1-j; stored vertex j is input
    // vertex (n-j) mod n.
    assert_eq!(
        stored.segments,
        vec![label(13), label(12), label(11), label(10)]
    );
    assert_eq!(
        stored.vertices,
        vec![label(20), label(23), label(22), label(21)]
    );
    assert!((ccw.area() - (6.0 + std::f64::consts::PI)).abs() < 1e-12);
}

#[test]
fn invalid_paths_are_rejected() {
    let t = tol();
    let (points, mut segments) = stadium();
    // An arc end off its circle.
    segments[1] = arc(3., 0., 1.01, true);
    assert_eq!(
        Boundary::path(points.clone(), segments, t),
        Err(Error::InvalidCurve("arc point off its circle"))
    );
    // A figure eight: two lens-like loops crossing.
    let crossing = Boundary::path(
        vec![p(0., 0.), p(4., 4.), p(4., 0.), p(0., 4.)],
        vec![Segment::Line, Segment::Line, Segment::Line, Segment::Line],
        t,
    );
    assert_eq!(crossing, Err(Error::SelfIntersection));
    // An arc doubling back along its neighbour: the line (0,0)-(2,0) then
    // an arc leaving (2,0) back towards (0,0) tangentially.
    let back = Boundary::path(
        vec![p(0., 0.), p(2., 0.), p(1., 1.)],
        vec![Segment::Line, arc(2., 1., 1., false), Segment::Line],
        t,
    );
    assert!(back.is_err(), "{back:?}");
    // A line crossing its neighbouring arc a second time: the upper half
    // circle from (0,0) over (2,2) to (4,0), then a line to (1,3) through
    // its top.
    let again = Boundary::path(
        vec![p(0., 0.), p(4., 0.), p(1., 3.)],
        vec![arc(2., 0., 2., false), Segment::Line, Segment::Line],
        t,
    );
    assert_eq!(again, Err(Error::SelfIntersection));
    // Two arcs of the path touching away from their joins: a rectangle
    // whose top and bottom bulge inwards until they meet.
    let pinched = Boundary::path(
        vec![p(0., 0.), p(4., 0.), p(4., 2.), p(0., 2.)],
        vec![
            arc(2., 1.5, 2.5, false),
            Segment::Line,
            arc(2., 0.5, 2.5, false),
            Segment::Line,
        ],
        t,
    );
    assert_eq!(pinched, Err(Error::SelfIntersection));
    // A zero-length segment.
    let short = Boundary::path(
        vec![p(0., 0.), p(0., 0.), p(1., 1.)],
        vec![
            Segment::Line,
            arc(0.5, 0.5, std::f64::consts::FRAC_1_SQRT_2, true),
            Segment::Line,
        ],
        t,
    );
    assert_eq!(short, Err(Error::Degenerate("path segment")));
    // A hole touching the outer path's arc.
    let (points, segments) = stadium();
    let outer = Boundary::path(points, segments, t).unwrap();
    let hole = Boundary::circle(p(3.4, 0.0), 0.6, t).unwrap();
    assert!(Profile::new(outer, vec![hole], t).is_err());
}

#[test]
fn points_classify_against_arcs() {
    let (points, segments) = notch();
    let solid = prism(points, segments);
    let at = |x: f64, y: f64, z: f64| solid.frame().point(p(x, y), z);
    let cases = [
        (at(5.0, 3.0, 1.0), Location::Inside),
        // Inside the notch's circle: outside the material.
        (at(5.0, 5.5, 1.0), Location::Outside),
        (at(5.0, 5.0, 1.0), Location::Boundary),
        (
            at(5.7, 6.0 - 0.714_142_842_854_285, 1.0),
            Location::Boundary,
        ),
        (at(3.5, 5.9, 1.0), Location::Inside),
        (at(11.0, 3.0, 1.0), Location::Outside),
        (at(5.0, 3.0, 2.0), Location::Boundary),
    ];
    for (point, want) in cases {
        assert_eq!(solid.classify(point).unwrap(), want, "{point:?}");
    }
    let (points, segments) = stadium();
    let solid = prism(points, segments);
    let at = |x: f64, y: f64| solid.frame().point(p(x, y), 1.0);
    for (x, y, want) in [
        (3.9, 0.0, Location::Inside),
        (4.0, 0.0, Location::Boundary),
        (4.1, 0.0, Location::Outside),
        (-0.9, 0.3, Location::Inside),
        (-0.9, 0.5, Location::Outside),
        (1.5, 1.0, Location::Boundary),
        (1.5, 0.99, Location::Inside),
    ] {
        assert_eq!(solid.classify(at(x, y)).unwrap(), want, "({x}, {y})");
    }
}

#[test]
fn closed_form_mass_lies_in_the_certified_enclosure() {
    for (points, segments) in [stadium(), notch()] {
        let solid = prism(points, segments);
        let m = solid.mass_properties();
        let e = solid.topology().mass_enclosure().unwrap();
        let inside = |x: f64, [lo, hi]: [f64; 2]| {
            let slack = 1e-12 * x.abs().max(1.0);
            lo - slack <= x && x <= hi + slack
        };
        assert!(inside(m.volume, e.volume));
        assert!(inside(m.surface_area, e.surface_area));
        for (c, r) in m.centroid.to_array().into_iter().zip(e.centroid) {
            assert!(inside(c, r), "{c} {r:?}");
        }
        for a in 0..3 {
            for b in 0..3 {
                assert!(inside(m.inertia[a][b], e.inertia[a][b]), "{a}{b}");
            }
        }
    }
}

#[test]
fn arc_prisms_split_and_fuse_by_height() {
    let (points, segments) = stadium();
    let solid = prism(points, segments);
    let ([low, high], _) = solid.split_at_height(OperationId(6), 0.5).unwrap();
    for piece in [&low, &high] {
        assert_eq!(piece.topology().check(piece.resolution()), Vec::new());
    }
    let (fused, _) = low.fuse_stacked(&high, OperationId(7)).unwrap();
    assert_eq!(fused.topology().check(fused.resolution()), Vec::new());
    assert!((fused.mass_properties().volume - solid.mass_properties().volume).abs() < 1e-12);
}

/// The neutral generator's arc prisms (`generate_brep_fixtures.prism`, the
/// independent seamed builder behind `brep-cases.txt`) are what the kernel's
/// builder makes of the same profiles: the same vertices, the same counts,
/// overlapping certified volumes.
#[test]
fn the_neutral_arc_prisms_are_builder_prisms() {
    let fixture = |name: &str| {
        let block = include_str!("../../fixtures/brep-cases.txt")
            .split("\nend")
            .find(|b| b.trim().lines().next() == Some(&format!("case {name}")))
            .unwrap();
        let (_, tolerance, parts) = brep_protocol::parse(block.trim());
        Topology::from_parts(parts, Tolerance::new(tolerance, 1e-12).unwrap()).unwrap()
    };
    let xy = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    type Case = (&'static str, Vec<Point2>, Vec<Segment>, f64, f64);
    let cases: [Case; 3] = [
        (
            "stadium",
            vec![p(0., -1.), p(3., -1.), p(3., 1.), p(0., 1.)],
            vec![
                Segment::Line,
                arc(3., 0., 1., true),
                Segment::Line,
                arc(0., 0., 1., true),
            ],
            0.0,
            1.0,
        ),
        (
            "concave_notch",
            vec![
                p(0., 0.),
                p(4., 0.),
                p(4., 2.),
                p(2.75, 2.),
                p(1.25, 2.),
                p(0., 2.),
            ],
            vec![
                Segment::Line,
                Segment::Line,
                Segment::Line,
                arc(2., 2., 0.75, false),
                Segment::Line,
                Segment::Line,
            ],
            0.0,
            1.5,
        ),
        (
            "half_disc",
            vec![p(-1., 0.), p(1., 0.)],
            vec![Segment::Line, arc(0., 0., 1., true)],
            0.0,
            0.4,
        ),
    ];
    for (name, points, segments, low, high) in cases {
        let profile = Profile::new(
            Boundary::path(points, segments, tol()).unwrap(),
            vec![],
            tol(),
        )
        .unwrap();
        let built = Solid::extrude_with(OperationId(1), profile, xy, low, high)
            .unwrap()
            .0;
        let neutral = fixture(name);
        let t = built.topology();
        assert_eq!(t.check(tol()), Vec::new(), "{name}");
        assert_eq!(t.occt_counts(), neutral.occt_counts(), "{name}");
        let key = |t: &Topology| {
            let mut v: Vec<[u64; 3]> = t
                .vertices()
                .iter()
                .map(|v| v.position.to_array().map(f64::to_bits))
                .collect();
            v.sort();
            v
        };
        assert_eq!(key(t), key(&neutral), "{name}");
        let (a, b) = (
            t.mass_enclosure().unwrap().volume,
            neutral.mass_enclosure().unwrap().volume,
        );
        assert!(a[0] <= b[1] && b[0] <= a[1], "{name}: {a:?} {b:?}");
    }
}
