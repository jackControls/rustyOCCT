use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::proximity::{
    closest_points_on_exact_spline_in, closest_points_on_spline, closest_points_on_spline_in,
    SplineProximityOptions,
};
use rusty_occt::{BSplineCurve3, Error, ExactBSplineCurve3, ExactKnotVector, Point3};
use std::cmp::Ordering::{Equal, Greater, Less};

fn rational(n: i32, d: i32) -> R {
    R::new(BigInt::from(n), BigInt::from(d))
}

#[test]
fn zero_minima_merge_shared_knots_and_intervals_without_losing_later_hits() {
    // The first span has strictly positive minimum distance. Later spans hit
    // the query at shared knots, stay there for two spans, leave, and return.
    let controls = [
        [-3., -1., 0.],
        [-2., 1., 0.],
        [0., 0., 0.],
        [1., 1., 0.],
        [0., 0., 0.],
        [0., 0., 0.],
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 0., 0.],
    ];
    let curve = BSplineCurve3::new(
        1,
        controls
            .into_iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect(),
        None,
        (0..=8).map(f64::from).collect(),
        vec![2, 1, 1, 1, 1, 1, 1, 1, 2],
    )
    .unwrap();
    let result = closest_points_on_spline(&curve, Point3::ORIGIN).unwrap();
    assert_eq!(result.compare_squared_distance(0.).unwrap(), Equal);
    assert_eq!(result.points().len(), 2);
    for (point, u) in result.points().iter().zip([2, 8]) {
        assert_eq!(point.rational_parameter(), Some(rational(u, 1)));
        for axis in 0..3 {
            assert_eq!(
                point.compare_coordinate(axis, &rational(0, 1)).unwrap(),
                Equal
            );
        }
    }
    assert_eq!(result.intervals().len(), 1);
    assert_eq!(
        result.intervals()[0].parameters(),
        &[rational(4, 1), rational(6, 1)]
    );

    let result =
        closest_points_on_spline_in(&curve, Point3::ORIGIN, 4.25, 5.75, Default::default())
            .unwrap();
    assert!(result.points().is_empty());
    assert_eq!(result.intervals().len(), 1);
    assert_eq!(
        result.intervals()[0].parameters(),
        &[rational(17, 4), rational(23, 4)]
    );

    // A common polynomial can have no real roots in the queried range.
    // C(t)=(t^2+1,0,0) on [-1,1] has its positive minimum at t=0.
    let away = BSplineCurve3::new(
        2,
        vec![
            Point3::new(2., 0., 0.),
            Point3::ORIGIN,
            Point3::new(2., 0., 0.),
        ],
        None,
        vec![-1., 1.],
        vec![3, 3],
    )
    .unwrap();
    let result = closest_points_on_spline(&away, Point3::ORIGIN).unwrap();
    assert_eq!(result.points().len(), 1);
    assert_eq!(
        result.points()[0].rational_parameter(),
        Some(rational(0, 1))
    );
    assert_eq!(result.compare_squared_distance(1.).unwrap(), Equal);
}

#[test]
fn exact_minima_survive_unrepresentable_parameter_coordinate_and_distance_views() {
    let huge = R::from_integer(BigInt::from(1) << 2048usize);
    let half = rational(1, 2);
    let zero = rational(0, 1);
    let one = rational(1, 1);
    let controls = vec![
        [zero.clone(), zero.clone(), zero.clone(), one.clone()],
        [one.clone(), zero.clone(), zero.clone(), one.clone()],
    ];
    let curve = ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new(1, vec![huge.clone(), &huge + &one], vec![2, 2]).unwrap(),
        controls,
    )
    .unwrap();
    let result = closest_points_on_exact_spline_in(
        &curve,
        &[half.clone(), huge.clone(), zero.clone()],
        &huge,
        &(&huge + &one),
        Default::default(),
    )
    .unwrap();
    assert_eq!(result.points().len(), 1);
    assert_eq!(result.points()[0].rational_parameter(), Some(&huge + &half));
    assert!(matches!(
        result.points()[0].parameter_bounds(),
        Err(Error::Unrepresentable(_))
    ));
    assert_eq!(
        result.points()[0].compare_coordinate(0, &half).unwrap(),
        Equal
    );
    assert_eq!(
        result
            .compare_squared_distance_exact(&(&huge * &huge))
            .unwrap(),
        Equal
    );
    assert!(matches!(
        result.squared_distance_bounds(),
        Err(Error::Unrepresentable(_))
    ));
    assert!(matches!(
        result.distance_bounds(),
        Err(Error::Unrepresentable(_))
    ));

    let curve = ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new(1, vec![zero.clone(), one.clone()], vec![2, 2]).unwrap(),
        vec![
            [huge.clone(), zero.clone(), zero.clone(), one.clone()],
            [&huge + &one, zero.clone(), zero.clone(), one.clone()],
        ],
    )
    .unwrap();
    let result = closest_points_on_exact_spline_in(
        &curve,
        &[&huge + &half, one.clone(), zero.clone()],
        &zero,
        &one,
        Default::default(),
    )
    .unwrap();
    assert_eq!(result.points()[0].rational_parameter(), Some(half));
    assert_eq!(result.compare_squared_distance(1.).unwrap(), Equal);
    assert!(matches!(
        result.points()[0].point_bounds(),
        Err(Error::Unrepresentable(_))
    ));
    assert_eq!(
        result.points()[0]
            .compare_coordinate(0, &(&huge + rational(1, 2)))
            .unwrap(),
        Equal
    );
}
fn parabola() -> BSplineCurve3 {
    BSplineCurve3::new(
        2,
        vec![
            Point3::new(-2., 4., 0.),
            Point3::new(0., -4., 0.),
            Point3::new(2., 4., 0.),
        ],
        None,
        vec![-2., 2.],
        vec![3, 3],
    )
    .unwrap()
}

#[test]
fn boundaries_singletons_corners_and_unclamped_domains_are_complete() {
    let line = BSplineCurve3::new(
        1,
        vec![Point3::ORIGIN, Point3::new(1., 0., 0.)],
        None,
        vec![0., 1.],
        vec![2, 2],
    )
    .unwrap();
    for (point, u, d2) in [
        (Point3::new(-1., 2., 0.), 0., 5.),
        (Point3::new(2., 2., 0.), 1., 5.),
        (Point3::new(0.5, 2., 0.), 0.5, 4.),
    ] {
        let result = closest_points_on_spline(&line, point).unwrap();
        assert_eq!(result.points().len(), 1);
        assert!(result.intervals().is_empty());
        assert_eq!(
            result.points()[0]
                .compare_parameter(&R::from_float(u).unwrap())
                .unwrap(),
            Equal
        );
        assert_eq!(result.compare_squared_distance(d2).unwrap(), Equal);
    }
    let result = closest_points_on_spline_in(
        &line,
        Point3::new(0.5, 2., 0.),
        0.125,
        0.25,
        SplineProximityOptions::default(),
    )
    .unwrap();
    assert_eq!(
        result.points()[0].rational_parameter(),
        Some(rational(1, 4))
    );
    assert_eq!(result.compare_squared_distance(4.0625).unwrap(), Equal);
    let result = closest_points_on_spline_in(
        &line,
        Point3::new(0., 2., 0.),
        0.25,
        0.25,
        SplineProximityOptions::default(),
    )
    .unwrap();
    assert_eq!(result.points().len(), 1);
    assert_eq!(
        result.points()[0].rational_parameter(),
        Some(rational(1, 4))
    );
    assert_eq!(result.compare_squared_distance(4.0625).unwrap(), Equal);
    let corner = BSplineCurve3::new(
        1,
        vec![
            Point3::new(0., 1., 0.),
            Point3::ORIGIN,
            Point3::new(1., 0., 0.),
        ],
        None,
        vec![0., 1., 2.],
        vec![2, 1, 2],
    )
    .unwrap();
    let result = closest_points_on_spline(&corner, Point3::new(-1., -1., 0.)).unwrap();
    assert_eq!(result.points().len(), 1);
    assert_eq!(
        result.points()[0].rational_parameter(),
        Some(rational(1, 1))
    );
    assert_eq!(result.compare_squared_distance(2.).unwrap(), Equal);
    let unclamped = BSplineCurve3::new(
        2,
        (0..5).map(|i| Point3::new(f64::from(i), 0., 0.)).collect(),
        None,
        (0..8).map(f64::from).collect(),
        vec![1; 8],
    )
    .unwrap();
    let result = closest_points_on_spline(&unclamped, Point3::new(1.5, 1., 0.)).unwrap();
    assert_eq!(
        result.points()[0].rational_parameter(),
        Some(rational(3, 1))
    );
    assert_eq!(result.compare_squared_distance(1.).unwrap(), Equal);
}

#[test]
fn all_equal_minima_and_sub_float_ties_are_distinguished_exactly() {
    let curve = parabola();
    let symmetric = closest_points_on_spline(&curve, Point3::new(0., 1., 0.)).unwrap();
    assert_eq!(symmetric.points().len(), 2);
    assert_eq!(symmetric.compare_squared_distance(0.75).unwrap(), Equal);
    assert_eq!(
        symmetric.points()[0]
            .compare_parameter(&rational(0, 1))
            .unwrap(),
        Less
    );
    assert_eq!(
        symmetric.points()[1]
            .compare_parameter(&rational(0, 1))
            .unwrap(),
        Greater
    );
    let shifted = closest_points_on_spline(&curve, Point3::new(2f64.powi(-70), 1., 0.)).unwrap();
    assert_eq!(shifted.points().len(), 1);
    assert_eq!(
        shifted.points()[0]
            .compare_parameter(&rational(0, 1))
            .unwrap(),
        Greater
    );
    assert_eq!(shifted.compare_squared_distance(0.75).unwrap(), Less);
    assert_eq!(
        shifted
            .distance_cmp(&symmetric, SplineProximityOptions::default())
            .unwrap(),
        Less
    );
    let degenerate = closest_points_on_spline(&curve, Point3::new(0., 0.5, 0.)).unwrap();
    assert_eq!(degenerate.points().len(), 1);
    assert_eq!(
        degenerate.points()[0].rational_parameter(),
        Some(rational(0, 1))
    );
    assert_eq!(degenerate.compare_squared_distance(0.25).unwrap(), Equal);
}

#[test]
fn whole_constant_distance_intervals_and_periodic_parameter_aliases_survive() {
    let circle = BSplineCurve3::new(
        2,
        vec![
            Point3::new(1., 0., 0.),
            Point3::new(1., 1., 0.),
            Point3::new(0., 1., 0.),
        ],
        Some(vec![1., 1., 2.]),
        vec![0., 1.],
        vec![3, 3],
    )
    .unwrap();
    for (z, d2) in [(0., 1.), (2., 5.)] {
        let result = closest_points_on_spline(&circle, Point3::new(0., 0., z)).unwrap();
        assert!(result.points().is_empty());
        assert_eq!(result.intervals().len(), 1);
        assert_eq!(
            result.intervals()[0].parameters(),
            &[rational(0, 1), rational(1, 1)]
        );
        assert_eq!(result.compare_squared_distance(d2).unwrap(), Equal);
    }
    let curve = BSplineCurve3::new_periodic(
        1,
        vec![
            Point3::ORIGIN,
            Point3::new(1., 0., 0.),
            Point3::new(0., 1., 0.),
        ],
        None,
        vec![0., 1., 2., 3.],
        vec![1; 4],
    )
    .unwrap();
    let result = closest_points_on_spline_in(
        &curve,
        Point3::new(-1., -1., 0.),
        -3.,
        6.,
        SplineProximityOptions::default(),
    )
    .unwrap();
    assert_eq!(result.points().len(), 4);
    assert!(result.intervals().is_empty());
    for (p, u) in result.points().iter().zip([-3, 0, 3, 6]) {
        assert_eq!(p.rational_parameter(), Some(rational(u, 1)));
        for axis in 0..3 {
            assert_eq!(p.compare_coordinate(axis, &rational(0, 1)).unwrap(), Equal);
        }
    }
    assert_eq!(result.compare_squared_distance(2.).unwrap(), Equal);
}

#[test]
fn computation_limits_are_atomic_and_invalid_queries_are_rejected() {
    let curve = parabola();
    let point = Point3::new(0., 1., 0.);
    for options in [
        SplineProximityOptions {
            max_spans: 0,
            ..Default::default()
        },
        SplineProximityOptions {
            max_candidates: 0,
            ..Default::default()
        },
        SplineProximityOptions {
            max_image_coefficient_updates: 0,
            ..Default::default()
        },
        SplineProximityOptions {
            max_image_coefficient_bits: 0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            closest_points_on_spline_in(&curve, point, -2., 2., options),
            Err(Error::ComputationLimit(_))
        ));
    }
    let result = closest_points_on_spline_in(
        &curve,
        point,
        0.,
        2.,
        SplineProximityOptions {
            max_image_coefficient_updates: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.points().len(), 1);
    for (a, b) in [(2., -2.), (-3., 2.), (0., 3.), (3., 3.)] {
        assert!(matches!(
            closest_points_on_spline_in(&curve, point, a, b, Default::default()),
            Err(Error::OutOfDomain(_))
        ));
    }
    assert!(matches!(
        closest_points_on_spline(&curve, Point3::new(f64::NAN, 0., 0.)),
        Err(Error::NonFinite(_))
    ));
    let bad = R::new_raw(BigInt::from(1), BigInt::from(0));
    assert!(matches!(
        closest_points_on_exact_spline_in(
            &curve.to_exact(),
            &[bad, rational(0, 1), rational(0, 1)],
            &rational(-2, 1),
            &rational(2, 1),
            Default::default()
        ),
        Err(Error::InvalidSpline(_))
    ));
}

#[test]
fn unused_complex_poles_do_not_destroy_exact_distance_images() {
    // Independently expanded Bernstein coefficients on [-2,2], degree 6:
    // X=t^4+t^2, W=t^2+1, Y=Z=0. Every Bernstein weight is positive.
    // Thus C(t)=(t^2,0,0), despite removable complex poles at +/-i.
    // About Q=(2,0,1), D=(t^2-2)^2+1 has two equal minima at +/-sqrt(2).
    let x = [
        rational(20, 1),
        rational(-4, 1),
        rational(-4, 3),
        rational(12, 5),
        rational(-4, 3),
        rational(-4, 1),
        rational(20, 1),
    ];
    let w = [
        rational(5, 1),
        rational(7, 3),
        rational(11, 15),
        rational(1, 5),
        rational(11, 15),
        rational(7, 3),
        rational(5, 1),
    ];
    let basis = ExactKnotVector::new(6, vec![rational(-2, 1), rational(2, 1)], vec![7, 7]).unwrap();
    let curve = ExactBSplineCurve3::from_homogeneous(
        basis,
        x.into_iter()
            .zip(w)
            .map(|(x, w)| [x, rational(0, 1), rational(0, 1), w])
            .collect(),
    )
    .unwrap();
    let result = closest_points_on_exact_spline_in(
        &curve,
        &[rational(2, 1), rational(0, 1), rational(1, 1)],
        &rational(-2, 1),
        &rational(2, 1),
        Default::default(),
    )
    .unwrap();
    assert_eq!(result.points().len(), 2);
    assert_eq!(result.compare_squared_distance(1.).unwrap(), Equal);
    assert_eq!(
        result.points()[0]
            .compare_parameter(&rational(-1, 1))
            .unwrap(),
        Less
    );
    assert_eq!(
        result.points()[1]
            .compare_parameter(&rational(1, 1))
            .unwrap(),
        Greater
    );
    for point in result.points() {
        assert_eq!(point.compare_coordinate(0, &rational(2, 1)).unwrap(), Equal);
    }
}
