#[allow(dead_code)]
#[path = "support/bezier_reference.rs"]
mod bernstein;
#[path = "support/knot_protocol.rs"]
mod protocol;
#[path = "support/knot_reference.rs"]
mod reference;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use reference::integer as r;
use rusty_occt::curve::{BezierExtractionOptions, DerivativeOrder as D, KnotSide as S};
use rusty_occt::{BSplineCurve3, BezierCurve3, Error, ExactBSplineCurve3, ExactKnotVector, Point3};

fn cubic() -> ExactBSplineCurve3 {
    BezierCurve3::new(
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(1., 2., 3.),
            Point3::new(-1., 4., 2.),
            Point3::new(5., 1., -2.),
        ],
        Some(vec![1., 2., 3., 1.]),
    )
    .unwrap()
    .as_bspline()
    .to_exact()
}

#[test]
fn exact_controls_and_removal_decisions_match_independent_equations() {
    let mut count = 0;
    for row in include_str!("../../fixtures/knot-editing.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let (input, expected) = row.split_once(" | ").unwrap();
        let input = protocol::parse(input);
        let (flags, curve) = protocol::execute(&input);
        assert_eq!(
            protocol::encode(&input.name, &flags, &curve),
            expected,
            "{}",
            input.name
        );
        if input.curve.degree() <= 5
            || count % 47 == 0
            || std::env::var_os("RUSTY_VERIFY_ALL_KNOT_ORACLES").is_some()
        {
            assert!(reference::equal(&input.curve, &curve), "{}", input.name);
            let mut current = input.curve.clone();
            for (op, u, m) in &input.operations {
                if *op == 'I' {
                    current = current.insert_knot(u, *m).unwrap();
                } else {
                    let actual = current.remove_knot(u, *m).unwrap();
                    assert_eq!(
                        actual,
                        reference::removed(&current, u, *m),
                        "{}",
                        input.name
                    );
                    if let Some(next) = actual {
                        current = next;
                    }
                }
            }
        }
        count += 1;
    }
    assert_eq!(count, 723);
}

#[test]
fn rational_edits_batch_aliases_and_unclamped_inactive_controls() {
    let original = cubic();
    let third = R::new(1.into(), 3.into());
    let next = &third + R::new(1.into(), BigInt::from(1) << 2048);
    let requests = vec![(next.clone(), 2), (third.clone(), 3), (third.clone(), 1)];
    let edited = original.refined(&requests).unwrap();
    assert!(reference::equal(&original, &edited));
    assert_eq!(
        edited,
        original
            .insert_knot(&third, 3)
            .unwrap()
            .insert_knot(&next, 2)
            .unwrap()
    );
    assert_eq!(
        edited
            .remove_knot(&next, 0)
            .unwrap()
            .unwrap()
            .remove_knot(&third, 0)
            .unwrap(),
        Some(original.clone())
    );
    assert_eq!(original.refined(&[]).unwrap(), original);
    assert_eq!(original.insert_knot(&third, 0).unwrap(), original);
    for u in [&r(0), &third, &next, &r(1)] {
        assert_eq!(
            original.exact_evaluate(u, D::Second, S::Automatic),
            edited.exact_evaluate(u, D::Second, S::Automatic)
        );
    }
    let inactive = BSplineCurve3::new(
        2,
        vec![
            Point3::new(97., 41., -93.),
            Point3::new(0., 0., 0.),
            Point3::new(1., 2., 3.),
            Point3::new(4., 3., 2.),
        ],
        Some(vec![3., 1., 2., 1.]),
        vec![-1., 0., 1., 2.],
        vec![2, 2, 1, 2],
    )
    .unwrap()
    .to_exact();
    assert_eq!(inactive.domain(), &[r(0), r(1)]);
    let refined = inactive.refined(&[(third, 2), (r(1), 2)]).unwrap();
    assert_eq!(
        inactive.homogeneous_poles()[0],
        refined.homogeneous_poles()[0]
    );
    assert!(reference::equal(&inactive, &refined));
}

#[test]
fn exact_evaluation_and_extraction_keep_original_units_and_one_sided_jets() {
    let original = cubic();
    let third = R::new(1.into(), 3.into());
    let edited = original
        .refined(&[(third.clone(), 3), (R::new(2.into(), 3.into()), 2)])
        .unwrap();
    for arc in edited.bezier_arcs().unwrap() {
        let [a, b] = arc.domain();
        let coefficients = reference::polynomial(&edited, a, b);
        assert_eq!(bernstein::coefficients(&arc), coefficients);
        let polynomial = bernstein::Arc {
            degree: edited.degree(),
            domain: [a.clone(), b.clone()],
            coefficients,
        };
        for u in [a.clone(), (a + b) / r(2), b.clone()] {
            let expected = polynomial.jet(&u);
            let jet = edited.exact_evaluate(&u, D::Second, S::Automatic).unwrap();
            assert_eq!(jet.position().coordinates(), &expected[0]);
            assert_eq!(jet.derivative(1).unwrap(), &expected[1]);
            assert_eq!(jet.derivative(2).unwrap(), &expected[2]);
            assert_eq!(arc.exact_evaluate(&u, D::Second).unwrap(), jet);
        }
        assert_eq!(arc.to_bspline().bezier_arcs().unwrap(), vec![arc]);
    }
    let kink = BSplineCurve3::new(
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
    .unwrap()
    .to_exact();
    assert!(matches!(
        kink.exact_evaluate(&r(1), D::First, S::Automatic),
        Err(Error::DiscontinuousDerivative)
    ));
    assert_ne!(
        kink.exact_evaluate(&r(1), D::First, S::Left)
            .unwrap()
            .derivative(1),
        kink.exact_evaluate(&r(1), D::First, S::Right)
            .unwrap()
            .derivative(1)
    );
    assert!(kink.remove_knot(&r(1), 0).unwrap().is_none());
    assert!(kink.exact_evaluate(&r(0), D::Position, S::Left).is_err());
    assert!(kink.exact_evaluate(&r(2), D::Position, S::Right).is_err());
}

#[test]
fn periodic_seam_aliases_origin_shift_and_exact_multiturn_extraction() {
    let basis =
        ExactKnotVector::new_periodic(2, vec![r(0), r(1), r(2), r(3), r(4)], vec![1; 5]).unwrap();
    let constant =
        ExactBSplineCurve3::from_homogeneous(basis, vec![[r(2), -r(3), r(5), r(1)]; 4]).unwrap();
    let raised = constant.refined(&[(r(4), 2), (r(0), 1)]).unwrap();
    assert_eq!(raised, constant.insert_knot(&r(0), 2).unwrap());
    assert_eq!(
        raised.remove_knot(&r(4), 1).unwrap(),
        Some(constant.clone())
    );
    let shifted = constant.remove_knot(&r(0), 0).unwrap().unwrap();
    assert_eq!(shifted.domain(), &[r(1), r(5)]);
    assert_eq!(
        Some(shifted.clone()),
        constant.remove_knot(&r(4), 0).unwrap()
    );
    assert!(reference::equal(&constant, &shifted));
    // Reinsert the old origin as a knot of the shifted domain. Its cyclic
    // representation need not have the old origin, but the function is equal.
    assert!(reference::equal(
        &constant,
        &shifted.insert_knot(&r(4), 1).unwrap()
    ));
    let offset = R::from_integer(BigInt::from(1) << 4096);
    let arcs = raised.bezier_arcs_in(&offset, &(&offset + r(8))).unwrap();
    assert_eq!(arcs.len(), 8);
    assert_eq!(arcs[0].domain()[0], offset);
    for u in [-r(9), r(0), r(4), r(13), offset] {
        assert_eq!(
            raised.exact_evaluate(&u, D::Second, S::Left).unwrap(),
            shifted.exact_evaluate(&u, D::Second, S::Right).unwrap()
        );
    }
    assert!(matches!(
        raised.bezier_arcs_with_options(&r(0), &r(8), BezierExtractionOptions { max_arcs: 7 }),
        Err(Error::ComputationLimit(_))
    ));
}

#[test]
fn nonconstant_periodic_origin_removal_preserves_canonical_control_order() {
    let original = BSplineCurve3::new_periodic(
        2,
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(2., 1., 0.),
            Point3::new(1., 3., 2.),
            Point3::new(-1., 1., 1.),
        ],
        None,
        vec![0., 1., 2., 3., 4.],
        vec![1; 5],
    )
    .unwrap()
    .to_exact();
    let shifted_basis = ExactKnotVector::new_periodic(
        2,
        vec![
            R::new(1.into(), 2.into()),
            r(1),
            r(2),
            r(3),
            r(4),
            R::new(9.into(), 2.into()),
        ],
        vec![1; 6],
    )
    .unwrap();
    let independent_controls = reference::recover(&original, &shifted_basis).unwrap();
    let fine = ExactBSplineCurve3::from_homogeneous(shifted_basis, independent_controls).unwrap();
    assert!(reference::equal(&original, &fine));
    let expected = reference::removed(&fine, &fine.domain()[0], 0).unwrap();
    assert_eq!(
        fine.remove_knot(&fine.domain()[0], 0).unwrap(),
        Some(expected.clone())
    );
    assert_eq!(
        fine.remove_knot(&fine.domain()[1], 0).unwrap(),
        Some(expected.clone())
    );
    assert_eq!(expected.domain(), &[r(1), r(5)]);
    let mut controls = original.homogeneous_poles().to_vec();
    controls.rotate_left(1);
    assert_eq!(expected.homogeneous_poles(), controls);
    for u in [-r(1), r(0), r(1), R::new(5.into(), 2.into()), r(5)] {
        assert_eq!(
            original.exact_evaluate(&u, D::First, S::Automatic),
            expected.exact_evaluate(&u, D::First, S::Automatic)
        );
    }
}

#[test]
fn retained_degree25_removal_exercises_complete_integer_equation_solver() {
    let bytes =
        include_bytes!("../../fuzz/regressions/knot_editing/constant-periodic-d25-removal.bin");
    assert_eq!(&bytes[..5], &[3, 24, 1, 1, 4]);
    let row = include_str!("../../fixtures/knot-editing.tsv")
        .lines()
        .find(|r| r.starts_with("constant_periodic_d25_remove_existing "))
        .unwrap();
    let (input, expected) = row.split_once(" | ").unwrap();
    let input = protocol::parse(input);
    let actual = reference::removed(&input.curve, &r(1), 0).unwrap();
    assert_eq!(
        input.curve.remove_knot(&r(1), 0).unwrap(),
        Some(actual.clone())
    );
    assert!(reference::equal(&input.curve, &actual));
    assert_eq!(protocol::encode(&input.name, &[true], &actual), expected);
}

#[test]
fn validation_resource_preflight_and_nonpositive_coarse_weights() {
    let curve = cubic();
    let bad = R::new_raw(1.into(), 0.into());
    for u in [&bad, &-r(1), &r(2)] {
        assert!(curve.insert_knot(u, 0).is_err());
        assert!(curve.exact_evaluate(u, D::Second, S::Automatic).is_err());
    }
    assert!(curve.refined(&[(r(0), 5)]).is_err());
    assert!(curve.insert_knot(&R::new(1.into(), 2.into()), 4).is_err());
    assert!(curve.remove_knot(&r(0), 0).is_err());
    assert!(curve.remove_knot(&R::new(1.into(), 2.into()), 0).is_err());
    assert!(curve.refined(&vec![(r(0), 0); 4097]).is_err());
    assert!(curve.bezier_arcs_in(&r(0), &bad).is_err());
    assert!(curve.bezier_arcs_in(&r(1), &r(0)).is_err());
    let raw = R::new_raw((-2).into(), (-6).into());
    assert_eq!(
        curve.insert_knot(&raw, 3).unwrap(),
        curve.insert_knot(&R::new(1.into(), 3.into()), 3).unwrap()
    );
    let basis = ExactKnotVector::new(1, vec![r(0), r(1)], vec![2, 2]).unwrap();
    assert!(ExactBSplineCurve3::from_homogeneous(
        basis.clone(),
        vec![[bad.clone(), r(0), r(0), r(1)]; 2]
    )
    .is_err());
    assert!(
        ExactBSplineCurve3::from_homogeneous(basis.clone(), vec![[r(0), r(0), r(0), r(0)]; 2])
            .is_err()
    );
    let normalized = ExactBSplineCurve3::from_homogeneous(
        basis,
        vec![
            [
                raw.clone(),
                r(0),
                r(0),
                R::new_raw((-2).into(), (-2).into())
            ];
            2
        ],
    )
    .unwrap();
    assert_eq!(
        normalized.homogeneous_poles()[0][0],
        R::new(1.into(), 3.into())
    );
    for (degree, knots, mults) in [
        (0, vec![r(0), r(1)], vec![1, 1]),
        (26, vec![r(0), r(1)], vec![27, 27]),
        (1, vec![r(1), r(0)], vec![2, 2]),
        (1, vec![r(0), bad], vec![2, 2]),
        (2, vec![r(0), r(1)], vec![1, 1]),
    ] {
        assert!(ExactKnotVector::new(degree, knots, mults).is_err());
    }
    let basis = ExactKnotVector::new(
        1,
        (0..4096).map(r).collect(),
        std::iter::once(2)
            .chain(std::iter::repeat_n(1, 4094))
            .chain(std::iter::once(2))
            .collect(),
    )
    .unwrap();
    let maximum =
        ExactBSplineCurve3::from_homogeneous(basis, vec![[r(0), r(0), r(0), r(1)]; 4096]).unwrap();
    assert!(matches!(
        maximum.insert_knot(&R::new(1.into(), 2.into()), 1),
        Err(Error::LimitExceeded(_))
    ));
    // Refining the signed quadratic weight row [1,-1/4,1] at 1/2 gives
    // positive weights [1,3/8,3/8,1]. Exact removal is outside our family.
    let basis = ExactKnotVector::new(
        2,
        vec![r(0), R::new(1.into(), 2.into()), r(1)],
        vec![3, 1, 3],
    )
    .unwrap();
    let weights = [
        r(1),
        R::new(3.into(), 8.into()),
        R::new(3.into(), 8.into()),
        r(1),
    ];
    let fine = ExactBSplineCurve3::from_homogeneous(
        basis,
        weights.into_iter().map(|w| [r(0), r(0), r(0), w]).collect(),
    )
    .unwrap();
    assert!(fine
        .remove_knot(&R::new(1.into(), 2.into()), 0)
        .unwrap()
        .is_none());
    assert!(reference::removed(&fine, &R::new(1.into(), 2.into()), 0).is_none());
    let tiny = R::new(1.into(), BigInt::from(1) << 2048);
    let basis = ExactKnotVector::new(1, vec![r(0), tiny.clone()], vec![2, 2]).unwrap();
    let huge_jet = ExactBSplineCurve3::from_homogeneous(
        basis,
        vec![[r(0), r(0), r(0), r(1)], [r(1), r(0), r(0), r(1)]],
    )
    .unwrap();
    assert_eq!(
        huge_jet
            .exact_evaluate(&tiny, D::Second, S::Automatic)
            .unwrap()
            .derivative(1)
            .unwrap()[0],
        r(1) / tiny
    );
    assert!(matches!(
        huge_jet.evaluate(0., D::First, S::Automatic),
        Err(Error::Unrepresentable(_))
    ));
}
