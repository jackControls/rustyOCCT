use num_rational::BigRational as R;
use rusty_occt::intersection::{spline_cylinder_in, spline_sphere_in, Cylinder3, Sphere3};
use rusty_occt::intersection::{
    spline_plane, spline_plane_in, spline_plane_in_with_options, spline_plane_with_options, Plane3,
    SplinePlaneContact, SplinePlaneOptions,
};
use rusty_occt::polynomial::RootIsolationOptions;
use rusty_occt::{BSplineCurve3, Error, Point3, Vec3};
use std::cmp::Ordering::{Equal, Greater, Less};
#[path = "support/exact_spline_edits.rs"]
mod exact_edits;
#[allow(dead_code)]
#[path = "support/knot_reference.rs"]
mod knot_identity;

fn bits(word: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(word, 16).unwrap())
}

#[test]
fn complete_results_match_independent_basis_and_continued_fraction_oracle() {
    check_fixtures(include_str!("../../fixtures/spline-plane.tsv"), 284);
}

#[test]
fn sphere_and_cylinder_results_match_independent_exact_oracle() {
    check_fixtures(include_str!("../../fixtures/spline-quadric.tsv"), 118);
}

#[test]
fn quadric_work_limits_fail_atomically_and_domains_are_explicit() {
    use rusty_occt::intersection::{
        spline_cylinder_with_options, spline_sphere_with_options, SplineSurfaceOptions,
    };
    let curve = BSplineCurve3::new(
        1,
        vec![Point3::new(-2., 0., 0.), Point3::new(2., 0., 0.)],
        None,
        vec![0., 1.],
        vec![2, 2],
    )
    .unwrap();
    let sphere = Sphere3::new(Point3::new(0., 0., 0.), 1.).unwrap();
    let cylinder = Cylinder3::new(Point3::new(0., 0., 0.), Vec3::new(0., 0., 1.), 1.).unwrap();
    for options in [
        SplineSurfaceOptions {
            max_spans: 0,
            ..Default::default()
        },
        SplineSurfaceOptions {
            root_isolation: RootIsolationOptions {
                max_subdivisions: 0,
            },
            ..Default::default()
        },
    ] {
        assert!(matches!(
            spline_sphere_with_options(&curve, &sphere, options),
            Err(Error::ComputationLimit(_))
        ));
        assert!(matches!(
            spline_cylinder_with_options(&curve, &cylinder, options),
            Err(Error::ComputationLimit(_))
        ));
    }
    for (first, last) in [(-1., 1.), (0., 2.), (1., 0.), (0.5, 0.5)] {
        assert!(matches!(
            spline_sphere_in(&curve, &sphere, first, last),
            Err(Error::OutOfDomain(_))
        ));
        assert!(matches!(
            spline_cylinder_in(&curve, &cylinder, first, last),
            Err(Error::OutOfDomain(_))
        ));
    }
    assert!(matches!(
        spline_sphere_in(&curve, &sphere, f64::NAN, 1.),
        Err(Error::NonFinite(_))
    ));
    assert!(matches!(
        spline_cylinder_in(&curve, &cylinder, 0., f64::INFINITY),
        Err(Error::NonFinite(_))
    ));
}

fn check_fixtures(text: &str, expected_count: usize) {
    let mut count = 0;
    for row in text.lines().filter(|s| !s.starts_with('#')) {
        let (input, expected) = row.split_once(" | ").unwrap();
        let mut w = input.split_whitespace();
        let name = w.next().unwrap();
        let kind = w.next().unwrap();
        let (surface_kind, kind) = kind.split_once(':').unwrap_or(("plane", kind));
        let degree = w.next().unwrap().parse().unwrap();
        let np: usize = w.next().unwrap().parse().unwrap();
        let nk: usize = w.next().unwrap().parse().unwrap();
        let mut number = || w.next().unwrap().parse().unwrap();
        let plane: [Point3; 3] = std::array::from_fn(|_| Point3::new(number(), number(), number()));
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for _ in 0..np {
            poles.push(Point3::new(number(), number(), number()));
            weights.push(number());
        }
        let mut knots = Vec::new();
        let mut mults = Vec::new();
        for _ in 0..nk {
            knots.push(w.next().unwrap().parse().unwrap());
            mults.push(w.next().unwrap().parse().unwrap());
        }
        let range = if kind.ends_with('T') {
            Some((
                w.next().unwrap().parse().unwrap(),
                w.next().unwrap().parse().unwrap(),
            ))
        } else {
            None
        };
        assert!(w.next().is_none());
        let curve = if kind.starts_with('P') {
            BSplineCurve3::new_periodic(degree, poles, Some(weights), knots, mults)
        } else {
            BSplineCurve3::new(degree, poles, Some(weights), knots, mults)
        }
        .unwrap();
        let (a, b) = range.unwrap_or_else(|| curve.domain());
        let result = match surface_kind {
            "sphere" => {
                spline_sphere_in(&curve, &Sphere3::new(plane[0], plane[2].x).unwrap(), a, b)
            }
            "cylinder" => spline_cylinder_in(
                &curve,
                &Cylinder3::new(
                    plane[0],
                    Vec3::new(plane[1].x, plane[1].y, plane[1].z),
                    plane[2].x,
                )
                .unwrap(),
                a,
                b,
            ),
            "plane" => {
                let plane = Plane3::through_points(plane[0], plane[1], plane[2]).unwrap();
                if range.is_some() {
                    spline_plane_in(&curve, &plane, a, b)
                } else {
                    spline_plane(&curve, &plane)
                }
            }
            _ => panic!("invalid primitive"),
        }
        .unwrap_or_else(|e| panic!("{name}: {e}"));
        let mut e = expected.split_whitespace();
        assert_eq!(
            result.points().len(),
            e.next().unwrap().parse::<usize>().unwrap(),
            "{name}"
        );
        assert_eq!(
            result.overlaps().len(),
            e.next().unwrap().parse::<usize>().unwrap(),
            "{name}"
        );
        assert_eq!(
            result.is_disjoint(),
            result.points().is_empty() && result.overlaps().is_empty()
        );
        for point in result.points() {
            for (c, bounds) in std::iter::once(point.parameter())
                .chain(point.coordinate_bounds())
                .enumerate()
            {
                let lo = bits(e.next().unwrap());
                let hi = bits(e.next().unwrap());
                assert_eq!(
                    (bounds.lower(), bounds.upper()),
                    (lo, hi),
                    "{name}, component {c}"
                );
                let compare = |x| {
                    if c == 0 {
                        point.compare_parameter(x)
                    } else {
                        point.compare_coordinate(c - 1, x)
                    }
                    .unwrap()
                };
                assert_eq!(
                    compare(lo),
                    if lo == hi { Equal } else { Greater },
                    "{name}"
                );
                assert_eq!(compare(hi), if lo == hi { Equal } else { Less }, "{name}");
            }
            let contact = match e.next().unwrap() {
                "B" => SplinePlaneContact::Boundary,
                "C" => SplinePlaneContact::Crossing,
                "T" => SplinePlaneContact::Tangent,
                _ => panic!(),
            };
            assert_eq!(point.contact(), contact, "{name}");
            let orders = std::array::from_fn(|_| {
                let n = e.next().unwrap().parse::<usize>().unwrap();
                (n != 0).then_some(n)
            });
            assert_eq!(point.multiplicities(), orders, "{name}");
            assert!(matches!(
                point.compare_coordinate(3, 0.),
                Err(Error::OutOfDomain(_))
            ));
            assert!(matches!(
                point.compare_parameter(f64::NAN),
                Err(Error::NonFinite(_))
            ));
        }
        for overlap in result.overlaps() {
            for (i, bounds) in overlap.parameter_bounds().into_iter().enumerate() {
                let lo = bits(e.next().unwrap());
                let hi = if range.is_some() {
                    bits(e.next().unwrap())
                } else {
                    lo
                };
                assert_eq!((bounds.lower(), bounds.upper()), (lo, hi), "{name}");
                assert_eq!(
                    overlap.compare_parameter(i, lo).unwrap(),
                    if lo == hi { Equal } else { Greater },
                    "{name}"
                );
                assert_eq!(
                    overlap.compare_parameter(i, hi).unwrap(),
                    if lo == hi { Equal } else { Less },
                    "{name}"
                );
            }
            assert!(matches!(
                overlap.compare_parameter(2, 0.),
                Err(Error::OutOfDomain(_))
            ));
        }
        assert!(e.next().is_none(), "{name}");
        // Every independent fixture also covers exact conversion, rational
        // refinement, and exact removal. The complete function identity is
        // checked by a Cox coefficient oracle independent of production.
        let original = curve.to_exact();
        let (fine, cuts) = exact_edits::refined(&original);
        assert!(knot_identity::equal(&original, &fine), "{name}");
        let restored = exact_edits::restored(&fine, &cuts);
        assert_eq!(restored, original, "{name}");
        for exact in [&original, &fine, &restored] {
            use rusty_occt::intersection::{
                exact_spline_cylinder_in, exact_spline_plane_in, exact_spline_sphere_in,
            };
            let (a, b) = (R::from_float(a).unwrap(), R::from_float(b).unwrap());
            let exact_result = match surface_kind {
                "plane" => exact_spline_plane_in(
                    exact,
                    &Plane3::through_points(plane[0], plane[1], plane[2]).unwrap(),
                    &a,
                    &b,
                ),
                "sphere" => exact_spline_sphere_in(
                    exact,
                    &Sphere3::new(plane[0], plane[2].x).unwrap(),
                    &a,
                    &b,
                ),
                "cylinder" => exact_spline_cylinder_in(
                    exact,
                    &Cylinder3::new(
                        plane[0],
                        Vec3::new(plane[1].x, plane[1].y, plane[1].z),
                        plane[2].x,
                    )
                    .unwrap(),
                    &a,
                    &b,
                ),
                _ => unreachable!(),
            }
            .unwrap_or_else(|e| panic!("{name}: exact intersection: {e}"));
            assert_eq!(exact_result.is_disjoint(), result.is_disjoint(), "{name}");
            let enclosed = exact_result.enclosed().unwrap();
            assert_eq!(enclosed.points().len(), result.points().len(), "{name}");
            for (i, (x, y)) in enclosed.points().iter().zip(result.points()).enumerate() {
                assert_eq!(x.parameter(), y.parameter(), "{name}");
                assert_eq!(x.coordinate_bounds(), y.coordinate_bounds(), "{name}");
                assert_eq!(x.contact(), y.contact(), "{name}");
                assert_eq!(x.multiplicities(), y.multiplicities(), "{name}");
                for value in [y.parameter().lower(), y.parameter().upper()] {
                    assert_eq!(
                        exact_result.points()[i]
                            .compare_parameter(&R::from_float(value).unwrap())
                            .unwrap(),
                        y.compare_parameter(value).unwrap(),
                        "{name}"
                    );
                }
            }
            assert_eq!(enclosed.overlaps(), result.overlaps(), "{name}");
        }
        count += 1;
    }
    assert_eq!(count, expected_count);
}

#[test]
fn overlap_touch_becomes_a_boundary_point_when_the_interval_is_outside_the_overlap() {
    let curve = BSplineCurve3::new(
        1,
        vec![
            Point3::new(0., 0., -1.),
            Point3::new(1., 0., 0.),
            Point3::new(2., 0., 0.),
            Point3::new(3., 0., 1.),
        ],
        None,
        vec![0., 1., 2., 3.],
        vec![2, 1, 1, 2],
    )
    .unwrap();
    let plane = Plane3::through_points(
        Point3::new(0., 0., 0.),
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
    )
    .unwrap();
    for (a, b, u, orders) in [
        (0.5, 1., 1., [Some(1), None]),
        (2., 2.5, 2., [None, Some(1)]),
    ] {
        let hit = spline_plane_in(&curve, &plane, a, b).unwrap();
        assert!(hit.overlaps().is_empty());
        assert_eq!(hit.points().len(), 1);
        assert_eq!(hit.points()[0].compare_parameter(u).unwrap(), Equal);
        assert_eq!(hit.points()[0].contact(), SplinePlaneContact::Boundary);
        assert_eq!(hit.points()[0].multiplicities(), orders);
    }
    for (a, b, expected) in [(0.5, 2.5, (1., 2.)), (1.25, 1.75, (1.25, 1.75))] {
        let hit = spline_plane_in(&curve, &plane, a, b).unwrap();
        assert!(hit.points().is_empty());
        assert_eq!(hit.overlaps().len(), 1);
        assert_eq!(hit.overlaps()[0].parameters(), expected);
    }
}

#[test]
fn periodic_ranges_keep_exact_parameters_and_preflight_traversal_limits() {
    let curve = BSplineCurve3::new_periodic(
        1,
        vec![Point3::new(0., 0., -1.), Point3::new(1., 0., 1.)],
        None,
        vec![0., 0.25, 0.5],
        vec![1; 3],
    )
    .unwrap();
    let plane = Plane3::through_points(
        Point3::new(0., 0., 0.),
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
    )
    .unwrap();
    let first = 2f64.powi(53);
    let last = first + 2.;
    let hits = spline_plane_in(&curve, &plane, first, last).unwrap();
    assert_eq!(hits.points().len(), 8);
    assert!(hits
        .points()
        .windows(2)
        .all(|p| p[0].parameter() == p[1].parameter()));
    for p in hits.points() {
        assert_eq!(p.compare_parameter(first).unwrap(), Greater);
        assert_eq!(p.compare_parameter(last).unwrap(), Less);
    }
    let options = SplinePlaneOptions {
        max_spans: 7,
        ..SplinePlaneOptions::default()
    };
    assert!(matches!(
        spline_plane_in_with_options(&curve, &plane, first, last, options),
        Err(Error::ComputationLimit(_))
    ));
    assert_eq!(
        spline_plane_in_with_options(
            &curve,
            &plane,
            first,
            last,
            SplinePlaneOptions {
                max_spans: 8,
                ..options
            }
        )
        .unwrap()
        .points()
        .len(),
        8
    );
    assert!(matches!(
        spline_plane_in(&curve, &plane, -f64::MAX, f64::MAX),
        Err(Error::ComputationLimit(_))
    ));
    for (a, b) in [(0., 0.), (1., 0.), (f64::NAN, 1.), (0., f64::INFINITY)] {
        assert!(spline_plane_in(&curve, &plane, a, b).is_err());
    }
    let nonperiodic =
        BSplineCurve3::new(1, curve.poles().to_vec(), None, vec![0., 1.], vec![2; 2]).unwrap();
    assert!(matches!(
        spline_plane_in(&nonperiodic, &plane, -0.25, 1.),
        Err(Error::OutOfDomain(_))
    ));
    assert!(matches!(
        spline_plane_in(&nonperiodic, &plane, 0., 1.25),
        Err(Error::OutOfDomain(_))
    ));
}

#[test]
fn subdivision_exhaustion_never_publishes_a_partial_intersection() {
    // Two distinct roots: z = 16(t-1/4)(t-3/4).
    let curve = BSplineCurve3::new(
        2,
        vec![
            Point3::new(0., 0., 3.),
            Point3::new(1., 0., -5.),
            Point3::new(2., 0., 3.),
        ],
        None,
        vec![0., 1.],
        vec![3, 3],
    )
    .unwrap();
    let plane = Plane3::through_points(
        Point3::new(0., 0., 0.),
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
    )
    .unwrap();
    assert!(matches!(
        spline_plane_with_options(
            &curve,
            &plane,
            RootIsolationOptions {
                max_subdivisions: 0
            }
        ),
        Err(Error::ComputationLimit(_))
    ));
    assert_eq!(spline_plane(&curve, &plane).unwrap().points().len(), 2);
}
