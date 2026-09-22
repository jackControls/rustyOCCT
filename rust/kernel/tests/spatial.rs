use rusty_occt::intersection::{
    line_plane, segment_plane, segment_triangle, IntersectionPoint, LinePlaneIntersection as L,
    SegmentPlaneIntersection as P, SegmentTriangleIntersection as T,
};
use rusty_occt::predicates::{in_sphere, orient3d, Orientation3, SphereLocation};
use rusty_occt::{Error, Plane3, Point3, Triangle3};

fn number(word: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(word, 16).unwrap())
}
fn points(words: &[&str]) -> Vec<Point3> {
    words
        .chunks_exact(3)
        .map(|p| Point3::new(number(p[0]), number(p[1]), number(p[2])))
        .collect()
}
fn orientation_sign(o: Orientation3) -> i8 {
    match o {
        Orientation3::Negative => -1,
        Orientation3::Coplanar => 0,
        Orientation3::Positive => 1,
    }
}
fn sphere_sign(o: SphereLocation) -> i8 {
    match o {
        SphereLocation::Outside => -1,
        SphereLocation::Boundary => 0,
        SphereLocation::Inside => 1,
    }
}

#[test]
fn spatial_predicates_match_independent_fraction_matrix_oracle() {
    let mut count = 0;
    for row in include_str!("../../fixtures/predicates3d.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let w: Vec<_> = row.split_whitespace().collect();
        let p = points(&w[2..w.len() - 1]);
        let expected = w[w.len() - 1];
        if w[0] == "o" {
            // All 24 permutations, with independently counted inversion parity.
            for a in 0..4 {
                for b in 0..4 {
                    for c in 0..4 {
                        for d in 0..4 {
                            let order = [a, b, c, d];
                            if (0..4).any(|i| (i + 1..4).any(|j| order[i] == order[j])) {
                                continue;
                            }
                            let inversions = (0..4)
                                .map(|i| (i + 1..4).filter(|&j| order[i] > order[j]).count())
                                .sum::<usize>();
                            let parity = if inversions % 2 == 0 { 1 } else { -1 };
                            assert_eq!(
                                orientation_sign(orient3d(p[a], p[b], p[c], p[d]).unwrap()),
                                parity * expected.parse::<i8>().unwrap(),
                                "{} {order:?}",
                                w[1]
                            );
                        }
                    }
                }
            }
        } else {
            for [a, b, c, d] in [[0, 1, 2, 3], [1, 0, 2, 3], [1, 2, 3, 0], [3, 2, 1, 0]] {
                let result = in_sphere(p[a], p[b], p[c], p[d], p[4]);
                if expected == "D" {
                    assert!(matches!(result, Err(Error::Degenerate(_))), "{}", w[1]);
                } else {
                    assert_eq!(
                        sphere_sign(result.unwrap()),
                        expected.parse::<i8>().unwrap(),
                        "{}",
                        w[1]
                    );
                }
            }
        }
        count += 1;
    }
    assert_eq!(count, 1648);
}

fn check_point(point: IntersectionPoint, expected: &[&str], label: &str) {
    let bounds = point.bounds();
    let low = bounds.min.to_array();
    let high = bounds.max.to_array();
    let actual = [
        point.parameter().lower(),
        point.parameter().upper(),
        low[0],
        high[0],
        low[1],
        high[1],
        low[2],
        high[2],
    ];
    for (i, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
        assert_eq!(actual, number(expected), "{label}: component {i}");
    }
    for (i, coordinate) in point.position().to_array().into_iter().enumerate() {
        assert!(
            coordinate.is_finite() && low[i] <= coordinate && coordinate <= high[i],
            "{label}"
        );
    }
}

#[test]
fn intersections_match_independent_barycentric_oracle_and_minimal_enclosures() {
    let mut count = 0;
    for row in include_str!("../../fixtures/intersections.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        count += 1;
        let w: Vec<_> = row.split_whitespace().collect();
        let p = points(&w[2..17]);
        let plane = Plane3::through_points(p[0], p[1], p[2]);
        if let Err(e) = plane {
            assert!(matches!(e, Error::Degenerate(_)));
            assert_eq!(w[17], "D");
            continue;
        }
        let plane = plane.unwrap();
        let result =
            match w[0] {
                "l" => line_plane(p[3], p[4], &plane).map(|r| match r {
                    L::Disjoint => ("N", vec![]),
                    L::Contained => ("C", vec![]),
                    L::Point(p) => ("P", vec![p]),
                }),
                "p" => segment_plane(p[3], p[4], &plane).map(|r| match r {
                    P::Disjoint => ("N", vec![]),
                    P::Contained => ("C", vec![]),
                    P::Point(p) => ("P", vec![p]),
                }),
                "t" => segment_triangle(p[3], p[4], &Triangle3::new(p[0], p[1], p[2]).unwrap())
                    .map(|r| match r {
                        T::Disjoint => ("N", vec![]),
                        T::Point(p) => ("P", vec![p]),
                        T::Overlap { start, end } => ("S", vec![start, end]),
                    }),
                _ => panic!("unknown fixture kind"),
            };
        match result {
            Err(Error::Unrepresentable(_)) => assert_eq!(w[17], "U", "{}", w[1]),
            Err(Error::Degenerate(_)) => assert_eq!(w[17], "D", "{}", w[1]),
            Ok((status, points)) => {
                assert_eq!(status, w[17], "{} {}", w[0], w[1]);
                assert_eq!(w.len(), 18 + 8 * points.len());
                for (point, expected) in points.into_iter().zip(w[18..].chunks_exact(8)) {
                    check_point(point, expected, w[1]);
                }
            }
            Err(e) => panic!("{}: {e}", w[1]),
        }
    }
    assert_eq!(count, 963);
}

#[test]
fn canonical_sign_and_sphere_convention_are_explicit() {
    let [a, b, c, d] = [
        Point3::ORIGIN,
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
        Point3::new(0., 0., 1.),
    ];
    assert_eq!(orient3d(a, b, c, d).unwrap(), Orientation3::Positive);
    assert_eq!(
        in_sphere(a, b, c, d, Point3::new(0.5, 0.5, 0.5)).unwrap(),
        SphereLocation::Inside
    );
    assert_eq!(
        in_sphere(a, b, c, d, Point3::new(1., 1., 1.)).unwrap(),
        SphereLocation::Boundary
    );
    assert_eq!(
        in_sphere(a, b, c, d, Point3::new(2., 2., 2.)).unwrap(),
        SphereLocation::Outside
    );
}

#[test]
fn spatial_apis_reject_every_nonfinite_coordinate() {
    let original = [
        Point3::ORIGIN,
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
        Point3::new(0., 0., -1.),
        Point3::new(0., 0., 1.),
    ];
    for index in 0..15 {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut values = original.map(Point3::to_array);
            values[index / 3][index % 3] = invalid;
            let p = values.map(|[x, y, z]| Point3::new(x, y, z));
            assert!(matches!(
                in_sphere(p[0], p[1], p[2], p[3], p[4]),
                Err(Error::NonFinite(_))
            ));
            if index < 12 {
                assert!(matches!(
                    orient3d(p[0], p[1], p[2], p[3]),
                    Err(Error::NonFinite(_))
                ));
            }
            if index < 9 {
                assert!(matches!(
                    Plane3::through_points(p[0], p[1], p[2]),
                    Err(Error::NonFinite(_))
                ));
                assert!(matches!(
                    Triangle3::new(p[0], p[1], p[2]),
                    Err(Error::NonFinite(_))
                ));
            } else {
                let triangle = Triangle3::new(p[0], p[1], p[2]).unwrap();
                assert!(matches!(
                    line_plane(p[3], p[4], triangle.plane()),
                    Err(Error::NonFinite(_))
                ));
                assert!(matches!(
                    segment_plane(p[3], p[4], triangle.plane()),
                    Err(Error::NonFinite(_))
                ));
                assert!(matches!(
                    segment_triangle(p[3], p[4], &triangle),
                    Err(Error::NonFinite(_))
                ));
            }
        }
    }
}

#[test]
fn coplanar_clipping_respects_winding_endpoints_and_tangency() {
    let [a, b, c] = [
        Point3::ORIGIN,
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
    ];
    let p = Point3::new(-1., 0.25, 0.);
    let q = Point3::new(2., 0.25, 0.);
    for [a, b, c] in [
        [a, b, c],
        [b, c, a],
        [c, a, b],
        [a, c, b],
        [b, a, c],
        [c, b, a],
    ] {
        let t = Triangle3::new(a, b, c).unwrap();
        let T::Overlap { start, end } = segment_triangle(p, q, &t).unwrap() else {
            panic!("overlap missing")
        };
        let T::Overlap {
            start: reverse_start,
            end: reverse_end,
        } = segment_triangle(q, p, &t).unwrap()
        else {
            panic!("reverse overlap missing")
        };
        assert_eq!(start.bounds(), reverse_end.bounds());
        assert_eq!(end.bounds(), reverse_start.bounds());
        assert_eq!(start.position(), Point3::new(0., 0.25, 0.));
        assert_eq!(end.position(), Point3::new(0.75, 0.25, 0.));
        let T::Point(touch) =
            segment_triangle(Point3::new(-1., 1., 0.), Point3::new(1., -1., 0.), &t).unwrap()
        else {
            panic!("tangent missing")
        };
        assert_eq!(touch.position(), Point3::ORIGIN);
    }
}
