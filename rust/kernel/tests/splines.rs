use rusty_occt::curve::{BSplineCurve3, BezierCurve3, DerivativeOrder as D, KnotSide as S};
use rusty_occt::{Error, Point3};

fn number(word: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(word, 16).unwrap())
}

#[test]
fn positions_and_derivatives_match_independent_exact_basis_oracle() {
    let mut count = 0;
    for row in include_str!("../../fixtures/splines.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let mut w = row.split_whitespace();
        let name = w.next().unwrap();
        let kind = w.next().unwrap();
        let degree = w.next().unwrap().parse().unwrap();
        let np: usize = w.next().unwrap().parse().unwrap();
        let nk: usize = w.next().unwrap().parse().unwrap();
        let u = number(w.next().unwrap());
        let side = match w.next().unwrap() {
            "L" => S::Left,
            "R" => S::Right,
            "A" => S::Automatic,
            _ => panic!(),
        };
        let order = match w.next().unwrap() {
            "0" => D::Position,
            "1" => D::First,
            "2" => D::Second,
            _ => panic!(),
        };
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for _ in 0..np {
            poles.push(Point3::new(
                number(w.next().unwrap()),
                number(w.next().unwrap()),
                number(w.next().unwrap()),
            ));
            weights.push(number(w.next().unwrap()));
        }
        let mut knots = Vec::new();
        let mut mults = Vec::new();
        for _ in 0..nk {
            knots.push(number(w.next().unwrap()));
            mults.push(w.next().unwrap().parse().unwrap());
        }
        let result = if kind == "B" {
            BezierCurve3::new(poles, Some(weights))
                .unwrap()
                .evaluate(u, order)
        } else {
            BSplineCurve3::new(degree, poles, Some(weights), knots, mults)
                .unwrap()
                .evaluate(u, order, side)
        };
        let status = w.next().unwrap();
        match result {
            Err(Error::Unrepresentable(_)) => assert_eq!(status, "U", "{name}"),
            Err(Error::OutOfDomain(_)) => assert_eq!(status, "O", "{name}"),
            Err(Error::DiscontinuousDerivative) => assert_eq!(status, "D", "{name}"),
            Err(e) => panic!("{name}: {e}"),
            Ok(value) => {
                assert_eq!(status, "P", "{name}");
                for v in [
                    Some(value.position_bounds()),
                    value.derivative_bounds(1),
                    value.derivative_bounds(2),
                ]
                .into_iter()
                .flatten()
                {
                    for component in v {
                        assert_eq!(component.lower(), number(w.next().unwrap()), "{name}");
                        assert_eq!(component.upper(), number(w.next().unwrap()), "{name}");
                        assert!(component.representative().is_finite());
                    }
                }
            }
        }
        assert!(w.next().is_none(), "{name}");
        count += 1;
    }
    assert_eq!(count, 385);
}

#[test]
fn repeated_knots_require_sides_only_when_the_requested_derivative_disagrees() {
    let curve = BSplineCurve3::new(
        1,
        vec![
            Point3::ORIGIN,
            Point3::new(1., 0., 0.),
            Point3::new(1., 1., 0.),
        ],
        None,
        vec![0., 1., 2.],
        vec![2, 1, 2],
    )
    .unwrap();
    assert_eq!(
        curve
            .evaluate(1., D::Position, S::Automatic)
            .unwrap()
            .position(),
        Point3::new(1., 0., 0.)
    );
    assert!(matches!(
        curve.evaluate(1., D::First, S::Automatic),
        Err(Error::DiscontinuousDerivative)
    ));
    let left = curve.evaluate(1., D::Second, S::Left).unwrap();
    let right = curve.evaluate(1., D::Second, S::Right).unwrap();
    assert_eq!(
        left.derivative_bounds(1)
            .unwrap()
            .map(|x| x.representative()),
        [1., 0., 0.]
    );
    assert_eq!(
        right
            .derivative_bounds(1)
            .unwrap()
            .map(|x| x.representative()),
        [0., 1., 0.]
    );
    let smooth = BSplineCurve3::new(
        1,
        vec![
            Point3::ORIGIN,
            Point3::new(1., 1., 0.),
            Point3::new(2., 2., 0.),
        ],
        None,
        vec![0., 1., 2.],
        vec![2, 1, 2],
    )
    .unwrap();
    let value = smooth.evaluate(1., D::Second, S::Automatic).unwrap();
    assert_eq!(
        value
            .derivative_bounds(2)
            .unwrap()
            .map(|x| x.representative()),
        [0.; 3]
    );
}

#[test]
fn overflow_in_derivatives_does_not_prevent_position_queries() {
    let curve = BezierCurve3::new(
        vec![
            Point3::new(-f64::MAX, 0., 0.),
            Point3::new(f64::MAX, 0., 0.),
        ],
        None,
    )
    .unwrap();
    let value = curve.evaluate(0.5, D::Position).unwrap();
    assert_eq!(value.position(), Point3::ORIGIN);
    assert_eq!(value.derivative_bounds(1), None);
    assert!(matches!(
        curve.evaluate(0.5, D::First),
        Err(Error::Unrepresentable(_))
    ));
}

#[test]
fn common_weight_scaling_and_coordinate_permutation_are_exact_invariants() {
    let poles = vec![
        Point3::new(1., 2., 3.),
        Point3::new(-5., 4., 7.),
        Point3::new(9., -3., 2.),
        Point3::new(-1., 6., 4.),
    ];
    let weights = [1., 2., 0.5, 4.];
    let base = BezierCurve3::new(poles.clone(), Some(weights.to_vec())).unwrap();
    let permuted = BezierCurve3::new(
        poles.iter().map(|p| Point3::new(p.z, p.x, p.y)).collect(),
        Some(weights.to_vec()),
    )
    .unwrap();
    for u in [0., 0.125, 0.5, 0.875, 1.] {
        let a = base.evaluate(u, D::Second).unwrap();
        let p = permuted.evaluate(u, D::Second).unwrap();
        for (x, y) in [
            (a.position_bounds(), p.position_bounds()),
            (
                a.derivative_bounds(1).unwrap(),
                p.derivative_bounds(1).unwrap(),
            ),
            (
                a.derivative_bounds(2).unwrap(),
                p.derivative_bounds(2).unwrap(),
            ),
        ] {
            assert_eq!([x[2], x[0], x[1]], y);
        }
        for power in [-1000, -100, 100, 1000] {
            let scaled = BezierCurve3::new(
                poles.clone(),
                Some(weights.map(|w| w * 2.0_f64.powi(power)).to_vec()),
            )
            .unwrap();
            assert_eq!(a, scaled.evaluate(u, D::Second).unwrap());
        }
    }
}

#[test]
fn original_parameter_units_determine_derivatives_and_closed_domain() {
    let curve = BSplineCurve3::new(
        1,
        vec![Point3::ORIGIN, Point3::new(8., 4., 2.)],
        None,
        vec![3., 5.],
        vec![2, 2],
    )
    .unwrap();
    assert_eq!(curve.domain(), (3., 5.));
    let value = curve.evaluate(4., D::Second, S::Automatic).unwrap();
    assert_eq!(value.position(), Point3::new(4., 2., 1.));
    assert_eq!(
        value
            .derivative_bounds(1)
            .unwrap()
            .map(|x| x.representative()),
        [4., 2., 1.]
    );
    for (u, side) in [
        (2., S::Automatic),
        (6., S::Automatic),
        (3., S::Left),
        (5., S::Right),
    ] {
        assert!(matches!(
            curve.evaluate(u, D::Position, side),
            Err(Error::OutOfDomain(_))
        ));
    }
    for u in [3., 5.] {
        assert!(curve.evaluate(u, D::Second, S::Automatic).is_ok());
    }
}

#[test]
fn construction_rejects_invalid_data_and_nonfinite_values() {
    let poles = vec![Point3::ORIGIN, Point3::new(1., 1., 1.)];
    for degree in [0, 26, usize::MAX] {
        assert!(BSplineCurve3::new(degree, poles.clone(), None, vec![0., 1.], vec![2, 2]).is_err());
    }
    for (knots, mults) in [
        (vec![], vec![]),
        (vec![0.], vec![2]),
        (vec![0., 0.], vec![2, 2]),
        (vec![1., 0.], vec![2, 2]),
        (vec![0., 1.], vec![0, 4]),
        (vec![0., 1.], vec![usize::MAX, 2]),
        (vec![0., 1.], vec![1, 2]),
        (vec![0., 1.], vec![2]),
    ] {
        assert!(BSplineCurve3::new(1, poles.clone(), None, knots, mults).is_err());
    }
    for weights in [vec![], vec![1.], vec![0., 1.], vec![-1., 1.]] {
        assert!(BezierCurve3::new(poles.clone(), Some(weights)).is_err());
    }
    assert!(BezierCurve3::new(vec![], None).is_err());
    assert!(BezierCurve3::new(vec![Point3::ORIGIN], None).is_err());
    assert!(BezierCurve3::new(vec![Point3::ORIGIN; 27], None).is_err());
    assert!(matches!(
        BSplineCurve3::new(
            1,
            vec![Point3::ORIGIN; 4097],
            None,
            vec![0., 1.],
            vec![2, 2]
        ),
        Err(Error::LimitExceeded(_))
    ));
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for coord in 0..3 {
            let mut xyz = [0.; 3];
            xyz[coord] = x;
            assert!(matches!(
                BezierCurve3::new(
                    vec![Point3::new(xyz[0], xyz[1], xyz[2]), Point3::ORIGIN],
                    None
                ),
                Err(Error::NonFinite(_))
            ));
        }
        assert!(matches!(
            BezierCurve3::new(poles.clone(), Some(vec![x, 1.])),
            Err(Error::NonFinite(_))
        ));
        assert!(matches!(
            BSplineCurve3::new(1, poles.clone(), None, vec![x, 1.], vec![2, 2]),
            Err(Error::NonFinite(_))
        ));
        let curve = BezierCurve3::new(poles.clone(), None).unwrap();
        assert!(matches!(
            curve.evaluate(x, D::Position),
            Err(Error::NonFinite(_))
        ));
    }
}
