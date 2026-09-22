#[path = "support/bezier_reference.rs"]
#[allow(dead_code)]
mod bezier_reference;
#[path = "support/surface_editing_protocol.rs"]
mod protocol;
#[path = "support/surface_editing_reference.rs"]
mod reference;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::curve::{DerivativeOrder as D, KnotSide};
use rusty_occt::surface::BezierPatchExtractionOptions as Options;
use rusty_occt::{BSplineSurface3, BezierSurface3, Error, KnotVector, Point3};
fn q(n: i32, d: i32) -> R {
    R::new(n.into(), d.into())
}

#[test]
fn complete_tensor_controls_and_independent_quotient_jets() {
    let mut count = 0;
    for row in include_str!("../../fixtures/surface-editing.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let (input, expected) = row.split_once(" | ").unwrap();
        let input = protocol::parse(input);
        let actual = protocol::execute(&input);
        assert_eq!(
            protocol::encode_fixture(&input.name, &actual),
            expected,
            "{}",
            input.name
        );
        // Python checks every complete result. A second independent Rust
        // coefficient/jet oracle checks selected cases here and every fuzz input.
        if count % 17 == 0 {
            let references: Vec<_> = reference::extract(&input.surface, input.rectangle)
                .into_iter()
                .flat_map(|p| reference::apply(p, input.op, input.elevation))
                .collect();
            assert_eq!(actual.len(), references.len());
            for (a, b) in actual.iter().zip(&references) {
                match (a, b) {
                    (protocol::Item::Patch(a), reference::Item::Patch(b)) => {
                        reference::check(a, b, [&q(1, 3), &q(2, 7)])
                    }
                    (protocol::Item::Curve(a), reference::Item::Curve(b)) => {
                        reference::curves::check(a, b, &q(2, 7))
                    }
                    _ => panic!("inconsistent result types"),
                }
            }
        }
        count += 1;
    }
    assert_eq!(count, 739);
}

#[test]
fn exact_boundaries_rational_restrictions_and_axis_commutation() {
    let poles = (0..3)
        .flat_map(|i| (0..4).map(move |j| Point3::new(i as f64, j as f64, (i * j) as f64)))
        .collect();
    let original = BezierSurface3::new(
        2,
        3,
        poles,
        Some((0..12).map(|i| (1 + i % 5) as f64).collect()),
    )
    .unwrap();
    let p = original.to_exact();
    assert_eq!(p, original.as_bspline().bezier_patches().unwrap()[0]);
    let (u, v) = (q(1, 3), q(2, 7));
    let [left, right] = p.split_u_at(&u).unwrap();
    let [bottom, top] = p.split_v_at(&v).unwrap();
    assert_eq!(left.u_iso(&u).unwrap(), right.u_iso(&u).unwrap());
    assert_eq!(bottom.v_iso(&v).unwrap(), top.v_iso(&v).unwrap());
    assert_eq!(
        left.split_v_at(&v).unwrap()[0],
        bottom.split_u_at(&u).unwrap()[0]
    );
    assert_eq!(p.exchanged_uv().exchanged_uv(), p);
    assert_eq!(p.u_reversed().u_reversed(), p);
    assert_eq!(p.v_reversed().v_reversed(), p);
    assert_eq!(p.u_reversed().v_reversed(), p.v_reversed().u_reversed());
    assert_eq!(p.exchanged_uv().u_reversed(), p.v_reversed().exchanged_uv());
    assert_eq!(
        p.elevated(5, 3).unwrap().elevated(5, 7).unwrap(),
        p.elevated(2, 7).unwrap().elevated(5, 7).unwrap()
    );
    assert_eq!(
        p.elevated(5, 7).unwrap().split_u_at(&u).unwrap(),
        [left.elevated(5, 7).unwrap(), right.elevated(5, 7).unwrap()]
    );
    assert_eq!(
        p.elevated(5, 7).unwrap().u_iso(&u).unwrap(),
        p.u_iso(&u).unwrap().elevated(7).unwrap()
    );
    assert_eq!(
        p.elevated(5, 7).unwrap().v_iso(&v).unwrap(),
        p.v_iso(&v).unwrap().elevated(5).unwrap()
    );
    let jet = p.exact_evaluate(&u, &v, D::Second).unwrap();
    assert_eq!(jet, left.exact_evaluate(&u, &v, D::Second).unwrap());
    assert_eq!(jet, right.exact_evaluate(&u, &v, D::Second).unwrap());
    let swapped = p.exchanged_uv().exact_evaluate(&v, &u, D::Second).unwrap();
    let reversed = p
        .u_reversed()
        .exact_evaluate(&(q(1, 1) - &u), &v, D::Second)
        .unwrap();
    assert_eq!(jet.position(), swapped.position());
    assert_eq!(jet.position(), reversed.position());
    for (du, dv) in [(1, 0), (0, 1), (2, 0), (0, 2), (1, 1)] {
        assert_eq!(jet.derivative(du, dv), swapped.derivative(dv, du));
        let expected =
            jet.derivative(du, dv)
                .unwrap()
                .clone()
                .map(|x| if du % 2 == 1 { -x } else { x });
        assert_eq!(reversed.derivative(du, dv), Some(&expected));
    }
    let iso = p.u_iso(&u).unwrap().exact_evaluate(&v, D::Second).unwrap();
    assert_eq!(iso.position(), jet.position());
    assert_eq!(iso.derivative(1), jet.derivative(0, 1));
    assert_eq!(iso.derivative(2), jet.derivative(0, 2));
    let tiny = R::new(1.into(), BigInt::from(1) << 512);
    let end = &u + &tiny;
    let trim = p.trim(&u, &end, &v, &(v.clone() + &tiny)).unwrap();
    assert_eq!(trim.exact_evaluate(&u, &v, D::Second).unwrap(), jet);
    assert_eq!(
        p.elevated(4, 5)
            .unwrap()
            .trim(&u, &end, &v, &(v.clone() + &tiny))
            .unwrap(),
        trim.elevated(4, 5).unwrap()
    );
    assert_eq!(trim.domain(), &[[u.clone(), end], [v.clone(), v + tiny]]);
    for axis in 0..2 {
        for end in 0..2 {
            let iso = if axis == 0 {
                p.u_iso(&p.domain()[0][end])
            } else {
                p.v_iso(&p.domain()[1][end])
            }
            .unwrap();
            let controls: Vec<_> = (0..=p.degrees()[1 - axis])
                .map(|i| {
                    let (u, v) = if axis == 0 {
                        (end * p.degrees()[0], i)
                    } else {
                        (i, end * p.degrees()[1])
                    };
                    p.homogeneous_poles()[u * (p.degrees()[1] + 1) + v].clone()
                })
                .collect();
            assert_eq!(iso.homogeneous_poles(), controls);
        }
    }
}

#[test]
fn invalid_parameters_and_unrepresentable_partials_keep_exact_data() {
    let tiny = f64::from_bits(1);
    let axis = KnotVector::new(1, vec![0., tiny], vec![2, 2]).unwrap();
    let s = BSplineSurface3::new(
        axis.clone(),
        axis,
        vec![
            Point3::new(0., 0., 0.),
            Point3::new(0., 1., 0.),
            Point3::new(1., 0., 0.),
            Point3::new(1., 1., 1.),
        ],
        None,
    )
    .unwrap();
    let p = s.bezier_patches().unwrap().remove(0);
    let zero = q(0, 1);
    let end = R::from_float(tiny).unwrap();
    let exact = p.exact_evaluate(&zero, &zero, D::Second).unwrap();
    assert_eq!(
        exact.derivative(1, 1).unwrap()[2],
        R::from_integer(BigInt::from(1) << 2148)
    );
    assert!(matches!(exact.enclosed(), Err(Error::Unrepresentable(_))));
    assert!(p.evaluate(0., 0., D::Position).is_ok());
    assert_eq!(
        p.exact_evaluate(&zero, &zero, D::Position)
            .unwrap()
            .derivative(1, 0),
        None
    );
    assert_eq!(exact.derivative(2, 1), None);
    assert_eq!(exact.derivative(0, 0), None);
    assert_eq!(p.trim(&zero, &end, &zero, &end).unwrap(), p);
    assert_eq!(p.elevated(1, 1).unwrap(), p);
    let bad = R::new_raw(1.into(), 0.into());
    for x in [&bad, &-q(1, 1), &q(1, 1)] {
        assert!(p.split_u_at(x).is_err());
        assert!(p.split_v_at(x).is_err());
        assert!(p.u_iso(x).is_err());
        assert!(p.v_iso(x).is_err());
        assert!(p.exact_evaluate(x, &zero, D::Position).is_err());
        assert!(p.exact_evaluate(&zero, x, D::Position).is_err());
        for i in 0..4 {
            let mut range = [&zero, &end, &zero, &end];
            range[i] = x;
            assert!(p.trim(range[0], range[1], range[2], range[3]).is_err());
        }
    }
    assert!(p.split_u_at(&zero).is_err());
    assert!(p.split_v_at(&end).is_err());
    assert!(p.trim(&zero, &zero, &zero, &end).is_err());
    assert!(p.trim(&zero, &end, &end, &zero).is_err());
    for (u, v) in [(0, 1), (1, 0), (26, 1), (1, 26)] {
        assert!(matches!(p.elevated(u, v), Err(Error::InvalidSurface(_))));
    }
    let raw = R::new_raw(-end.numer(), -end.denom() * 2);
    assert_eq!(
        p.split_u_at(&raw).unwrap(),
        p.split_u_at(&(&end / q(2, 1))).unwrap()
    );
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            p.evaluate(x, 0., D::Position),
            Err(Error::NonFinite(_))
        ));
        assert!(matches!(
            p.evaluate(0., x, D::Position),
            Err(Error::NonFinite(_))
        ));
    }
}

#[test]
fn cartesian_patch_and_control_limits_precede_geometry_allocation() {
    let tiny = f64::from_bits(1);
    let axis = KnotVector::new_periodic(1, vec![0., tiny, 2. * tiny], vec![1, 1, 1]).unwrap();
    let s =
        BSplineSurface3::new(axis.clone(), axis, vec![Point3::new(0., 0., 0.); 4], None).unwrap();
    assert_eq!(
        s.bezier_patches_with_options(
            0.,
            4. * tiny,
            0.,
            4. * tiny,
            Options {
                max_patches: 16,
                max_controls: 64
            }
        )
        .unwrap()
        .len(),
        16
    );
    for options in [
        Options {
            max_patches: 15,
            max_controls: 64,
        },
        Options {
            max_patches: 16,
            max_controls: 63,
        },
        Options {
            max_patches: 0,
            max_controls: usize::MAX,
        },
        Options {
            max_patches: usize::MAX,
            max_controls: 0,
        },
    ] {
        assert!(matches!(
            s.bezier_patches_with_options(0., 4. * tiny, 0., 4. * tiny, options),
            Err(Error::ComputationLimit(_))
        ));
    }
    for rect in [
        [-f64::MAX, f64::MAX, 0., tiny],
        [0., tiny, -f64::MAX, f64::MAX],
        [-f64::MAX, f64::MAX, -f64::MAX, f64::MAX],
    ] {
        assert!(matches!(
            s.bezier_patches_in(rect[0], rect[1], rect[2], rect[3]),
            Err(Error::ComputationLimit(_))
        ));
    }
    for rect in [
        [0., 0., 0., tiny],
        [0., tiny, tiny, 0.],
        [0., f64::NAN, 0., tiny],
        [0., tiny, 0., f64::INFINITY],
    ] {
        assert!(s
            .bezier_patches_in(rect[0], rect[1], rect[2], rect[3])
            .is_err());
    }
    let p = s.bezier_patches().unwrap();
    assert_eq!(
        p[0].u_iso(&R::from_float(tiny).unwrap()).unwrap(),
        p[2].u_iso(&R::from_float(tiny).unwrap()).unwrap()
    );
}

#[test]
fn knot_rectangles_keep_distinct_one_sided_partial_derivatives() {
    let axis = KnotVector::new(1, vec![0., 1., 2.], vec![2, 1, 2]).unwrap();
    let poles = (0..3)
        .flat_map(|i| {
            (0..3).map(move |j| Point3::new((i * i) as f64, (j * j) as f64, (i * j) as f64))
        })
        .collect();
    let s = BSplineSurface3::new(axis.clone(), axis, poles, None).unwrap();
    let patches = s.bezier_patches().unwrap();
    assert!(matches!(
        s.evaluate(1., 1., D::First, [KnotSide::Automatic; 2]),
        Err(Error::DiscontinuousDerivative)
    ));
    for (i, p) in patches.iter().enumerate() {
        let sides = [
            if i / 2 == 0 {
                KnotSide::Left
            } else {
                KnotSide::Right
            },
            if i % 2 == 0 {
                KnotSide::Left
            } else {
                KnotSide::Right
            },
        ];
        assert_eq!(
            p.evaluate(1., 1., D::Second).unwrap(),
            s.evaluate(1., 1., D::Second, sides).unwrap()
        );
    }
    let a = patches[0]
        .exact_evaluate(&q(1, 1), &q(1, 1), D::Second)
        .unwrap();
    let b = patches[3]
        .exact_evaluate(&q(1, 1), &q(1, 1), D::Second)
        .unwrap();
    assert_eq!(a.position(), b.position());
    assert_ne!(a.derivative(1, 0), b.derivative(1, 0));
    assert_ne!(a.derivative(0, 1), b.derivative(0, 1));
}
