use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::intersection::*;
use rusty_occt::{
    BSplineCurve3, Error, ExactBSplineCurve3, ExactKnotVector, Point3, ScalarInterval,
};
use std::cmp::Ordering::{self, Equal, Greater, Less};

fn r(n: i64) -> R {
    R::from_integer(n.into())
}
fn q(n: i64, d: i64) -> R {
    R::new(n.into(), d.into())
}
fn xyz(x: R, y: R, z: R) -> [R; 3] {
    [x, y, z]
}
fn exact(controls: Vec<[R; 4]>, domain: [R; 2]) -> ExactBSplineCurve3 {
    let n = controls.len();
    ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new(n - 1, domain.to_vec(), vec![n; 2]).unwrap(),
        controls,
    )
    .unwrap()
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
        vec![3; 2],
    )
    .unwrap()
}
fn backwards() -> ExactBSplineCurve3 {
    exact(
        [4, -4, 4].map(|x| [r(x), r(0), r(0), r(1)]).to_vec(),
        [r(-2), r(2)],
    )
}
fn interval_parameters(interval: &SplineLinearOverlap, a: R, b: R) {
    for (p, expected) in interval.endpoints().iter().zip([a, b]) {
        assert_eq!(p.compare_parameter(&expected).unwrap(), Equal);
    }
}
fn bounds(interval: ScalarInterval, compare: impl Fn(&R) -> Ordering) {
    let a = R::from_float(interval.lower()).unwrap();
    let b = R::from_float(interval.upper()).unwrap();
    if a == b {
        assert_eq!(compare(&a), Equal);
    } else {
        assert_eq!(
            interval
                .lower()
                .to_bits()
                .abs_diff(interval.upper().to_bits()),
            1
        );
        assert_eq!(compare(&a), Greater);
        assert_eq!(compare(&b), Less);
    }
}

#[test]
fn complete_secants_tangencies_and_spatial_misses() {
    let curve = parabola();
    let line = spline_line(&curve, Point3::new(2., 1., 0.), Point3::new(3., 1., 0.)).unwrap();
    assert!(line.overlaps().is_empty());
    assert_eq!(line.points().len(), 2);
    for (point, t) in line.points().iter().zip([-1, 1]) {
        assert_eq!(point.compare_parameter(&r(t)).unwrap(), Equal);
        assert_eq!(point.compare_coordinate(0, &r(t)).unwrap(), Equal);
        assert_eq!(point.compare_coordinate(1, &r(1)).unwrap(), Equal);
        assert_eq!(point.compare_coordinate(2, &r(0)).unwrap(), Equal);
        assert_eq!(point.compare_linear_parameter(&r(t - 2)).unwrap(), Equal);
    }
    assert!(
        spline_segment(&curve, Point3::new(2., 1., 0.), Point3::new(3., 1., 0.))
            .unwrap()
            .is_disjoint()
    );
    let a = Point3::new(-2., 1., 0.);
    let b = Point3::new(2., 1., 0.);
    let forward = spline_segment(&curve, a, b).unwrap();
    let reverse = spline_segment(&curve, b, a).unwrap();
    assert_eq!(forward.points().len(), 2);
    assert_eq!(reverse.points().len(), 2);
    for (i, (x, y)) in forward.points().iter().zip(reverse.points()).enumerate() {
        let t = q(1 + 2 * i as i64, 4);
        assert_eq!(x.parameter_cmp(y), Equal);
        assert_eq!(x.compare_linear_parameter(&t).unwrap(), Equal);
        assert_eq!(y.compare_linear_parameter(&(r(1) - &t)).unwrap(), Equal);
        bounds(x.linear_parameter_bounds().unwrap(), |v| t.cmp(v));
    }
    let tangent =
        spline_segment(&curve, Point3::new(-1., 0., 0.), Point3::new(1., 0., 0.)).unwrap();
    assert_eq!(tangent.points().len(), 1);
    assert_eq!(tangent.points()[0].compare_parameter(&r(0)).unwrap(), Equal);
    assert!(
        spline_line(&curve, Point3::new(-2., 1., 1.), Point3::new(2., 1., 1.))
            .unwrap()
            .is_disjoint()
    );
    // C(t)=(3t,3t^2,t^3) has independently expanded integer Bernstein controls.
    let space = exact(
        [(-6, 12, -8), (-2, -4, 8), (2, -4, -8), (6, 12, 8)]
            .map(|(x, y, z)| [r(x), r(y), r(z), r(1)])
            .to_vec(),
        [r(-2), r(2)],
    );
    let hits =
        exact_spline_segment(&space, &xyz(r(-3), r(3), r(-1)), &xyz(r(3), r(3), r(1))).unwrap();
    assert_eq!(hits.points().len(), 2);
    for (point, t) in hits.points().iter().zip([-1, 1]) {
        assert_eq!(point.compare_parameter(&r(t)).unwrap(), Equal);
    }
    let touch = exact_spline_line(&space, &xyz(r(-1), r(0), r(0)), &xyz(r(1), r(0), r(0))).unwrap();
    assert_eq!(touch.points().len(), 1);
    assert_eq!(touch.points()[0].compare_parameter(&r(0)).unwrap(), Equal);
}

#[test]
fn backtracking_has_complete_algebraic_overlap_intervals() {
    let curve = backwards();
    let a = xyz(r(1), r(0), r(0));
    let b = xyz(r(2), r(0), r(0));
    let line = exact_spline_line(&curve, &a, &b).unwrap();
    assert!(line.points().is_empty());
    assert_eq!(line.overlaps().len(), 1);
    interval_parameters(&line.overlaps()[0], r(-2), r(2));
    let result = exact_spline_segment(&curve, &a, &b).unwrap();
    assert!(result.points().is_empty());
    assert_eq!(result.overlaps().len(), 2);
    for (interval, positive) in result.overlaps().iter().zip([false, true]) {
        let [lo, hi] = interval.endpoints();
        let (radical, integer) = if positive { (hi, lo) } else { (lo, hi) };
        assert_eq!(
            integer
                .compare_parameter(&r(if positive { 1 } else { -1 }))
                .unwrap(),
            Equal
        );
        bounds(radical.parameter_bounds().unwrap(), |v| {
            if positive {
                if v < &r(0) {
                    Greater
                } else {
                    r(2).cmp(&(v * v))
                }
            } else if v > &r(0) {
                Less
            } else {
                (v * v).cmp(&r(2))
            }
        });
        assert_eq!(radical.compare_coordinate(0, &r(2)).unwrap(), Equal);
        assert_eq!(radical.compare_linear_parameter(&r(1)).unwrap(), Equal);
        assert_eq!(integer.compare_coordinate(0, &r(1)).unwrap(), Equal);
        assert_eq!(integer.compare_linear_parameter(&r(0)).unwrap(), Equal);
    }
    let joined = exact_spline_segment(&curve, &xyz(r(0), r(0), r(0)), &a).unwrap();
    assert!(joined.points().is_empty());
    assert_eq!(joined.overlaps().len(), 1);
    interval_parameters(&joined.overlaps()[0], r(-1), r(1));
    let touch =
        exact_spline_segment(&curve, &xyz(r(-1), r(0), r(0)), &xyz(r(0), r(0), r(0))).unwrap();
    assert_eq!(touch.points().len(), 1);
    assert!(touch.overlaps().is_empty());
    assert_eq!(touch.points()[0].compare_parameter(&r(0)).unwrap(), Equal);
    let point = exact_spline_segment(&curve, &a, &a).unwrap();
    assert_eq!(point.points().len(), 2);
    for hit in point.points() {
        assert_eq!(hit.compare_linear_parameter(&r(0)).unwrap(), Equal);
    }
}

#[test]
fn constant_curves_and_singletons_keep_the_parameter_contract() {
    let at = xyz(q(1, 3), r(2), r(-4));
    let curve = exact(
        [1, 3, 2]
            .map(|w| [&at[0] * r(w), &at[1] * r(w), &at[2] * r(w), r(w)])
            .to_vec(),
        [q(1, 7), q(10, 7)],
    );
    let result = exact_spline_segment(&curve, &at, &at).unwrap();
    assert!(result.points().is_empty());
    assert_eq!(result.overlaps().len(), 1);
    interval_parameters(&result.overlaps()[0], q(1, 7), q(10, 7));
    for endpoint in result.overlaps()[0].endpoints() {
        assert_eq!(endpoint.compare_linear_parameter(&r(0)).unwrap(), Equal);
        for (i, value) in at.iter().enumerate() {
            assert_eq!(endpoint.compare_coordinate(i, value).unwrap(), Equal);
            bounds(endpoint.coordinate_bounds().unwrap()[i], |x| value.cmp(x));
        }
    }
    let singleton = exact_spline_segment_in(&curve, &at, &at, &q(3, 7), &q(3, 7)).unwrap();
    assert_eq!(singleton.points().len(), 1);
    assert!(singleton.overlaps().is_empty());
    assert_eq!(
        singleton.points()[0].compare_parameter(&q(3, 7)).unwrap(),
        Equal
    );
    assert!(
        exact_spline_segment(&curve, &xyz(r(0), r(0), r(0)), &xyz(r(0), r(0), r(0)))
            .unwrap()
            .is_disjoint()
    );
}

#[test]
fn shared_knots_periodic_aliases_and_trimmed_overlaps_are_complete() {
    let curve = ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new_periodic(1, (0..=3).map(r).collect(), vec![1; 4]).unwrap(),
        vec![
            [r(0), r(0), r(0), r(1)],
            [r(1), r(0), r(0), r(1)],
            [r(0), r(1), r(0), r(1)],
        ],
    )
    .unwrap();
    let a = xyz(r(-1), r(0), r(0));
    let b = xyz(r(2), r(0), r(0));
    let result = exact_spline_segment_in(&curve, &a, &b, &r(-3), &r(6)).unwrap();
    assert_eq!(result.points().len(), 1);
    assert_eq!(result.points()[0].compare_parameter(&r(6)).unwrap(), Equal);
    assert_eq!(result.overlaps().len(), 3);
    for (interval, (lo, hi)) in result.overlaps().iter().zip([(-3, -2), (0, 1), (3, 4)]) {
        interval_parameters(interval, r(lo), r(hi));
    }
    let trim = exact_spline_segment_in(&curve, &a, &b, &q(-5, 2), &q(7, 2)).unwrap();
    assert!(trim.points().is_empty());
    assert_eq!(trim.overlaps().len(), 3);
    for (interval, (lo, hi)) in
        trim.overlaps()
            .iter()
            .zip([(q(-5, 2), r(-2)), (r(0), r(1)), (r(3), q(7, 2))])
    {
        interval_parameters(interval, lo, hi);
    }
    let origin = xyz(r(0), r(0), r(0));
    let aliases = exact_spline_segment_in(&curve, &origin, &origin, &r(-3), &r(6)).unwrap();
    assert!(aliases.overlaps().is_empty());
    assert_eq!(aliases.points().len(), 4);
    for (point, t) in aliases.points().iter().zip([-3, 0, 3, 6]) {
        assert_eq!(point.compare_parameter(&r(t)).unwrap(), Equal);
        assert_eq!(point.compare_coordinate(0, &r(0)).unwrap(), Equal);
    }
    let huge = R::from_integer(BigInt::from(1) << 2048usize);
    assert!(matches!(
        exact_spline_line_in_with_options(
            &curve,
            &a,
            &b,
            &r(0),
            &huge,
            SplineLinearOptions {
                max_spans: 1,
                ..Default::default()
            }
        ),
        Err(Error::ComputationLimit(_))
    ));
}

#[test]
fn exact_knot_and_degree_edits_preserve_rational_intersections() {
    // H=(t+t^2, t^2-t+2/9, 0, 1+t), t in [0,1], with original
    // curve parameter u=-5/7+2t. The axis hits are exactly t=1/3,2/3.
    let curve = exact(
        vec![
            [r(0), q(2, 9), r(0), r(1)],
            [q(1, 3), q(-1, 9), r(0), q(4, 3)],
            [r(1), q(-1, 9), r(0), q(5, 3)],
            [r(2), q(2, 9), r(0), r(2)],
        ],
        [q(-5, 7), q(9, 7)],
    );
    let edited = curve
        .insert_knot(&r(0), 3)
        .unwrap()
        .elevated(8)
        .unwrap()
        .insert_knot(&q(1, 7), 5)
        .unwrap();
    let a = xyz(r(0), r(0), r(0));
    let b = xyz(r(1), r(0), r(0));
    let original = exact_spline_line(&curve, &a, &b).unwrap();
    let result = exact_spline_line(&edited, &a, &b).unwrap();
    assert_eq!(original.points().len(), 2);
    assert_eq!(result.points().len(), 2);
    for (i, (x, y)) in original.points().iter().zip(result.points()).enumerate() {
        let t = q(1 + i as i64, 3);
        let parameter = q(-5, 7) + r(2) * &t;
        assert_eq!(x.parameter_cmp(y), Equal);
        assert_eq!(y.compare_parameter(&parameter).unwrap(), Equal);
        assert_eq!(y.compare_coordinate(0, &t).unwrap(), Equal);
        assert_eq!(y.compare_linear_parameter(&t).unwrap(), Equal);
    }
    let clipped = exact_spline_segment(&edited, &xyz(q(1, 2), r(0), r(0)), &b).unwrap();
    assert_eq!(clipped.points().len(), 1);
    assert_eq!(
        clipped.points()[0].compare_parameter(&q(13, 21)).unwrap(),
        Equal
    );
    let unclamped = ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new(2, (-2..=4).map(r).collect(), vec![1; 7]).unwrap(),
        [(-3, 1), (-1, 2), (1, 3), (3, 4)]
            .map(|(x, w)| [r(x * w), r(0), r(0), r(w)])
            .to_vec(),
    )
    .unwrap();
    let entire = exact_spline_line(&unclamped, &a, &b).unwrap();
    assert_eq!(entire.overlaps().len(), 1);
    interval_parameters(&entire.overlaps()[0], r(0), r(2));
}

#[test]
fn sub_float_gaps_and_coincident_point_aliases_are_not_merged() {
    let epsilon = R::new(1.into(), BigInt::from(1) << 180usize);
    // C(t)=((t-1)(t-1-epsilon),0,0), t in [0,2].
    let curve = exact(
        [r(1) + &epsilon, r(-1), r(1) - &epsilon]
            .map(|x| [x, r(0), r(0), r(1)])
            .to_vec(),
        [r(0), r(2)],
    );
    let zero = xyz(r(0), r(0), r(0));
    let result = exact_spline_segment(&curve, &zero, &xyz(r(2), r(0), r(0))).unwrap();
    assert!(result.points().is_empty());
    assert_eq!(result.overlaps().len(), 2);
    interval_parameters(&result.overlaps()[0], r(0), r(1));
    interval_parameters(&result.overlaps()[1], r(1) + &epsilon, r(2));
    let a = &result.overlaps()[0].endpoints()[1];
    let b = &result.overlaps()[1].endpoints()[0];
    assert_eq!(a.parameter_cmp(b), Less);
    assert_eq!(
        a.parameter_bounds().unwrap().upper(),
        b.parameter_bounds().unwrap().lower()
    );
    let points = exact_spline_segment(&curve, &zero, &zero).unwrap();
    assert_eq!(points.points().len(), 2);
    assert!(points.overlaps().is_empty());
    for point in points.points() {
        for component in 0..3 {
            assert_eq!(point.compare_coordinate(component, &r(0)).unwrap(), Equal);
        }
    }
    assert_eq!(points.points()[0].parameter_cmp(&points.points()[1]), Less);
}

#[test]
fn exact_results_survive_unrepresentable_coordinates_and_both_parameters() {
    let huge = R::from_integer(BigInt::from(1) << 2048usize);
    let tiny = r(1) / &huge;
    let curve = exact(
        vec![[huge.clone(), r(0), r(0), r(1)]; 2],
        [huge.clone(), &huge + r(1)],
    );
    let result = exact_spline_line(&curve, &xyz(r(0), r(0), r(0)), &xyz(tiny, r(0), r(0))).unwrap();
    assert!(result.points().is_empty());
    assert_eq!(result.overlaps().len(), 1);
    let point = &result.overlaps()[0].endpoints()[0];
    assert_eq!(point.compare_parameter(&huge).unwrap(), Equal);
    assert_eq!(point.compare_coordinate(0, &huge).unwrap(), Equal);
    assert_eq!(
        point.compare_linear_parameter(&(&huge * &huge)).unwrap(),
        Equal
    );
    assert!(matches!(
        point.parameter_bounds(),
        Err(Error::Unrepresentable(_))
    ));
    assert!(matches!(
        point.coordinate_bounds(),
        Err(Error::Unrepresentable(_))
    ));
    assert!(matches!(
        point.linear_parameter_bounds(),
        Err(Error::Unrepresentable(_))
    ));
    let segment = exact_spline_segment(
        &curve,
        &xyz(&huge - r(1), r(0), r(0)),
        &xyz(&huge + r(1), r(0), r(0)),
    )
    .unwrap();
    assert_eq!(segment.overlaps().len(), 1);
    assert_eq!(
        segment.overlaps()[0].endpoints()[0]
            .compare_linear_parameter(&q(1, 2))
            .unwrap(),
        Equal
    );
}

#[test]
fn invalid_inputs_and_exhausted_budgets_return_atomic_errors() {
    let curve = parabola();
    let a = Point3::new(-2., 1., 0.);
    let b = Point3::new(2., 1., 0.);
    for options in [
        SplineLinearOptions {
            max_spans: 0,
            ..Default::default()
        },
        SplineLinearOptions {
            max_candidates: 0,
            ..Default::default()
        },
        SplineLinearOptions {
            root_isolation: rusty_occt::polynomial::RootIsolationOptions {
                max_subdivisions: 0,
            },
            ..Default::default()
        },
    ] {
        assert!(matches!(
            spline_line_in_with_options(&curve, a, b, -2., 2., options),
            Err(Error::ComputationLimit(_))
        ));
    }
    assert!(matches!(
        spline_line(&curve, a, a),
        Err(Error::Degenerate(_))
    ));
    assert!(spline_line(&curve, Point3::new(f64::NAN, 0., 0.), b).is_err());
    assert!(spline_segment_in(&curve, a, b, f64::NEG_INFINITY, 2.).is_err());
    assert!(spline_segment_in(&curve, a, b, 2., -2.).is_err());
    assert!(spline_segment_in(&curve, a, b, -3., -3.).is_err());
    let exact = curve.to_exact();
    let bad = R::new_raw(1.into(), 0.into());
    assert!(exact_spline_segment(
        &exact,
        &xyz(bad.clone(), r(0), r(0)),
        &xyz(r(1), r(0), r(0))
    )
    .is_err());
    let result = spline_line(&curve, a, b).unwrap();
    assert!(result.points()[0].compare_parameter(&bad).is_err());
    assert!(result.points()[0].compare_linear_parameter(&bad).is_err());
    assert!(result.points()[0].compare_coordinate(3, &r(0)).is_err());
}

#[test]
fn rational_circle_keeps_irrational_parameters_and_exact_endpoint_contacts() {
    let circle = exact(
        vec![
            [r(1), r(0), r(0), r(1)],
            [r(1), r(1), r(0), r(1)],
            [r(0), r(2), r(0), r(2)],
        ],
        [r(0), r(1)],
    );
    let hit =
        exact_spline_segment(&circle, &xyz(r(0), r(0), r(0)), &xyz(r(1), r(1), r(0))).unwrap();
    assert_eq!(hit.points().len(), 1);
    let point = &hit.points()[0];
    bounds(point.parameter_bounds().unwrap(), |x| {
        let y = x + r(1);
        if y < r(0) {
            Greater
        } else {
            r(2).cmp(&(&y * &y))
        }
    });
    let compare = |x: &R| {
        if x < &r(0) {
            Greater
        } else {
            q(1, 2).cmp(&(x * x))
        }
    };
    let [x, y, z] = point.coordinate_bounds().unwrap();
    bounds(x, compare);
    bounds(y, compare);
    bounds(z, |x| r(0).cmp(x));
    bounds(point.linear_parameter_bounds().unwrap(), compare);
    let chord =
        exact_spline_segment(&circle, &xyz(r(1), r(0), r(0)), &xyz(r(0), r(1), r(0))).unwrap();
    assert_eq!(chord.points().len(), 2);
    for (point, t) in chord.points().iter().zip([0, 1]) {
        assert_eq!(point.compare_parameter(&r(t)).unwrap(), Equal);
    }
}
