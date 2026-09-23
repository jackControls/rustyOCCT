#[path = "support/bezier_protocol.rs"]
mod protocol;
#[path = "support/bezier_reference.rs"]
mod reference;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use reference::integer;
use rusty_occt::curve::{BezierExtractionOptions, DerivativeOrder as D, KnotSide};
use rusty_occt::{BSplineCurve3, BezierCurve3, Error, Point3};

#[test]
fn exact_controls_match_independent_polynomials_and_original_parameter_jets() {
    let mut count = 0;
    let all_oracles = std::env::var_os("RUSTY_VERIFY_ALL_BEZIER_ORACLES").is_some();
    for row in include_str!("../../fixtures/bezier-editing.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let (input, expected) = row.split_once(" | ").unwrap();
        let input = protocol::parse(input);
        let actual = protocol::execute(&input);
        assert_eq!(
            protocol::encode(&input.name, &actual),
            expected,
            "{}",
            input.name
        );
        // Every case has an independent Python coefficient oracle. Exercise a
        // second, Rust coefficient oracle on all moderate degrees and selected
        // high degrees; fuzzing uses it for every generated case.
        if all_oracles || input.curve.degree() <= 8 || count % 13 == 0 {
            let references = reference::extract(&input.curve, input.first, input.last)
                .into_iter()
                .flat_map(|c| reference::apply(c, input.op, input.elevation))
                .collect::<Vec<_>>();
            assert_eq!(actual.len(), references.len());
            for (a, b) in actual.iter().zip(&references) {
                reference::check(a, b, &R::new(1.into(), 3.into()));
            }
        }
        count += 1;
    }
    assert_eq!(count, 636);
}

#[test]
fn rational_cuts_survive_rounding_and_edit_composition() {
    let original = BezierCurve3::new(
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(1., 2., 3.),
            Point3::new(-2., 1., 4.),
        ],
        Some(vec![1., 3., 2.]),
    )
    .unwrap();
    let curve = original.to_exact();
    assert_eq!(curve, original.as_bspline().bezier_arcs().unwrap()[0]);
    let third = R::new(1.into(), 3.into());
    let tiny = R::new(1.into(), BigInt::from(1) << 2048);
    let next = &third + &tiny;
    let trimmed = curve.trim(&third, &next).unwrap();
    assert_eq!(trimmed.domain(), &[third.clone(), next.clone()]);
    assert_eq!(
        trimmed.exact_evaluate(&third, D::Second).unwrap(),
        curve.exact_evaluate(&third, D::Second).unwrap()
    );
    assert_eq!(
        trimmed.exact_evaluate(&next, D::Second).unwrap(),
        curve.exact_evaluate(&next, D::Second).unwrap()
    );
    assert_eq!(
        curve.elevated(25).unwrap().trim(&third, &next).unwrap(),
        trimmed.elevated(25).unwrap()
    );
    assert_eq!(trimmed.reversed().reversed(), trimmed);
    let [left, right] = curve.split_at(&third).unwrap();
    assert_eq!(
        left.homogeneous_poles().last(),
        right.homogeneous_poles().first()
    );
    assert_eq!(
        left.exact_evaluate(&third, D::Second).unwrap(),
        right.exact_evaluate(&third, D::Second).unwrap()
    );
    assert_eq!(
        curve.elevated(6).unwrap().split_at(&third).unwrap(),
        [left.elevated(6).unwrap(), right.elevated(6).unwrap()]
    );
    let reversed = curve.reversed();
    let a = curve.exact_evaluate(&third, D::Second).unwrap();
    let b = reversed
        .exact_evaluate(&(integer(1) - &third), D::Second)
        .unwrap();
    assert_eq!(a.position(), b.position());
    assert_eq!(
        a.derivative(1).unwrap().clone().map(|x| -x),
        *b.derivative(1).unwrap()
    );
    assert_eq!(a.derivative(2), b.derivative(2));
}

#[test]
fn validation_limits_and_exact_evaluation_survive_float_overflow() {
    let poles = vec![Point3::new(0., 0., 0.), Point3::new(1., 0., 0.)];
    let curve = BezierCurve3::new(poles.clone(), Some(vec![1., 2.]))
        .unwrap()
        .to_exact();
    let bad = R::new_raw(1.into(), 0.into());
    for u in [&bad, &-integer(1), &integer(2)] {
        assert!(curve.split_at(u).is_err());
        assert!(curve.exact_evaluate(u, D::Second).is_err());
        assert!(curve.trim(u, &integer(1)).is_err());
    }
    for u in [integer(0), integer(1)] {
        assert!(curve.split_at(&u).is_err());
    }
    assert!(curve.trim(&integer(1), &integer(0)).is_err());
    assert!(curve.trim(&integer(0), &integer(0)).is_err());
    assert!(curve.trim(&integer(0), &bad).is_err());
    assert!(curve.elevated(0).is_err());
    assert!(curve.elevated(26).is_err());
    assert_eq!(curve.elevated(1).unwrap(), curve);
    assert_eq!(curve.trim(&integer(0), &integer(1)).unwrap(), curve);
    let raw = R::new_raw((-2).into(), (-6).into());
    assert_eq!(
        curve.split_at(&raw).unwrap(),
        curve.split_at(&R::new(1.into(), 3.into())).unwrap()
    );
    for u in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            curve.evaluate(u, D::Position),
            Err(Error::NonFinite(_))
        ));
    }
    let tiny = f64::from_bits(1);
    let spline = BSplineCurve3::new(1, poles.clone(), None, vec![0., tiny], vec![2, 2]).unwrap();
    let arc = spline.bezier_arcs().unwrap().remove(0);
    let exact = arc.exact_evaluate(&integer(0), D::Second).unwrap();
    assert_eq!(
        exact.derivative(1).unwrap()[0],
        R::new(BigInt::from(1) << 1074, 1.into())
    );
    assert!(matches!(exact.enclosed(), Err(Error::Unrepresentable(_))));
    assert!(arc.evaluate(0., D::Position).is_ok());
    assert!(matches!(
        spline.bezier_arcs_with_options(0., tiny, BezierExtractionOptions { max_arcs: 0 }),
        Err(Error::ComputationLimit(_))
    ));
    for (a, b) in [(0., 0.), (tiny, 0.), (-tiny, tiny), (0., f64::NAN)] {
        assert!(spline.bezier_arcs_in(a, b).is_err());
    }
    let periodic =
        BSplineCurve3::new_periodic(1, poles, None, vec![0., tiny, 2. * tiny], vec![1, 1, 1])
            .unwrap();
    assert!(matches!(
        periodic.bezier_arcs_in(-f64::MAX, f64::MAX),
        Err(Error::ComputationLimit(_))
    ));
    assert_eq!(
        periodic
            .bezier_arcs_with_options(0., 4. * tiny, BezierExtractionOptions { max_arcs: 4 })
            .unwrap()
            .len(),
        4
    );
    assert!(matches!(
        periodic.bezier_arcs_with_options(0., 4. * tiny, BezierExtractionOptions { max_arcs: 3 }),
        Err(Error::ComputationLimit(_))
    ));
}

#[test]
fn extracted_endpoints_keep_distinct_one_sided_derivatives() {
    let curve = BSplineCurve3::new(
        1,
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(1., 0., 0.),
            Point3::new(1., 1., 0.),
        ],
        None,
        vec![0., 1., 2.],
        vec![2, 1, 2],
    )
    .unwrap();
    let arcs = curve.bezier_arcs().unwrap();
    let a = arcs[0].exact_evaluate(&integer(1), D::First).unwrap();
    let b = arcs[1].exact_evaluate(&integer(1), D::First).unwrap();
    assert_eq!(a.position(), b.position());
    assert_ne!(a.derivative(1), b.derivative(1));
    assert!(matches!(
        curve.evaluate(1., D::First, KnotSide::Automatic),
        Err(Error::DiscontinuousDerivative)
    ));
    assert_eq!(
        arcs[0].evaluate(1., D::First).unwrap(),
        curve.evaluate(1., D::First, KnotSide::Left).unwrap()
    );
    assert_eq!(
        arcs[1].evaluate(1., D::First).unwrap(),
        curve.evaluate(1., D::First, KnotSide::Right).unwrap()
    );
}
