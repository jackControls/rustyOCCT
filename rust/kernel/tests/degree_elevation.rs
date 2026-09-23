#[path = "support/knot_protocol.rs"]
mod curves;
#[allow(dead_code)]
#[path = "support/knot_reference.rs"]
mod knot_reference;
#[path = "support/degree_reference.rs"]
mod reference;
#[path = "support/surface_knot_protocol.rs"]
mod surfaces;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::curve::{DerivativeOrder as D, KnotSide as S};
use rusty_occt::{Error, ExactBSplineCurve3, ExactBSplineSurface3, ExactKnotVector};

fn r(n: i32) -> R {
    R::from_integer(n.into())
}

#[test]
fn independent_degree_equations_reject_unrepresentable_higher_coefficients() {
    // The first span is identically zero, so a bad solver that stops at full
    // rank can miss the quadratic coefficient appearing only in the last span.
    let old = ExactKnotVector::new(2, vec![r(0), r(1), r(2)], vec![3, 2, 3]).unwrap();
    let linear = ExactKnotVector::new(1, vec![r(0), r(1), r(2)], vec![2, 1, 2]).unwrap();
    let mut controls = vec![[r(0), r(0), r(0), r(1)]; 5];
    assert!(knot_reference::recover_controls(&old, &[&controls], &linear).is_some());
    controls[4][0] = R::new(1.into(), BigInt::from(1) << 4096);
    assert!(knot_reference::recover_controls(&old, &[&controls], &linear).is_none());
    controls[4][0] = R::from_integer(BigInt::from(1) << 4096);
    assert!(knot_reference::recover_controls(&old, &[&controls], &linear).is_none());
}

#[test]
fn retained_periodic_tensor_timeout_reconstructs_every_control() {
    let data =
        include_bytes!("../../fuzz/regressions/degree_elevation/uniform-periodic-both-24-25.bin");
    let axis = ExactKnotVector::new_periodic(24, (0..28).map(r).collect(), vec![1; 28]).unwrap();
    let controls = (0..27 * 27)
        .map(|i| {
            std::array::from_fn(|c| {
                let byte = data[16 + 4 * i + c];
                r(if c == 3 {
                    1 + i32::from(byte % 8)
                } else {
                    byte as i8 as i32
                })
            })
        })
        .collect();
    let surface = ExactBSplineSurface3::from_homogeneous(axis.clone(), axis, controls).unwrap();
    let result = surface.elevated(25, 25).unwrap();
    assert_eq!(result, reference::elevated_surface(&surface, [25, 25]));
    assert_eq!(result.domain(), surface.domain());
    assert_eq!(result.pole_counts(), [54, 54]);
    assert_eq!(
        result.exchanged_uv(),
        surface.exchanged_uv().elevated(25, 25).unwrap()
    );
}

#[test]
fn minimal_periodic_support_reconstructs_every_source_unit_column() {
    // With only p+1 cyclic poles, the support reaches furthest into the
    // neighboring periods. Pack every source unit column into x/y/z fields
    // of a transverse tensor row; positive unit weights keep valid geometry.
    // The global periodic Cox solver checks all map coefficients without
    // using the production working extension or its clamping/cropping.
    for p in 1usize..25 {
        let n = p + 1;
        let width = n.div_ceil(3).max(2);
        let mut v_mults = vec![1; width];
        v_mults[0] = 2;
        v_mults[width - 1] = 2;
        let v =
            ExactKnotVector::new(1, (0..width).map(|i| r(i as i32)).collect(), v_mults).unwrap();
        let mut seams = vec![1, p.div_ceil(2), p];
        seams.sort_unstable();
        seams.dedup();
        for seam in seams {
            let u = ExactKnotVector::new_periodic(
                p,
                vec![r(0), R::new(1.into(), 257.into()), r(3)],
                vec![seam, n - seam, seam],
            )
            .unwrap();
            let controls = (0..n)
                .flat_map(|i| {
                    (0..width).map(move |j| {
                        std::array::from_fn(|c| r(i32::from(c == 3 || i == 3 * j + c)))
                    })
                })
                .collect();
            let surface = ExactBSplineSurface3::from_homogeneous(u, v.clone(), controls).unwrap();
            let mut targets = vec![p + 1, 25];
            targets.dedup();
            for q in targets {
                assert_eq!(
                    surface.elevated(q, 1).unwrap(),
                    reference::elevated_surface(&surface, [q, 1]),
                    "degree {p} to {q}, seam multiplicity {seam}"
                );
            }
        }
    }
}

#[test]
fn every_curve_control_matches_independent_coefficient_equations() {
    let mut count = 0;
    for row in include_str!("../../fixtures/degree-elevation-curves.tsv").lines() {
        let (input, expected) = row.split_once(" | ").unwrap();
        let input = curves::parse(input);
        let (flags, result) = curves::execute(&input);
        assert_eq!(
            curves::encode(&input.name, &flags, &result),
            expected,
            "{}",
            input.name
        );
        assert_eq!(result.domain(), input.curve.domain());
        assert_eq!(result.is_periodic(), input.curve.is_periodic());
        if (input.curve.degree() <= 3 && result.degree() <= 4)
            || std::env::var_os("RUSTY_VERIFY_ALL_DEGREE_ORACLES").is_some()
        {
            assert_eq!(
                result,
                reference::elevated(&input.curve, result.degree()),
                "{}",
                input.name
            );
        }
        if input.curve.degree() <= 3 && result.degree() <= 4 {
            for numerator in [0, 1, 3, 4] {
                let u = &result.domain()[0]
                    + (&result.domain()[1] - &result.domain()[0]) * r(numerator) / r(4);
                for side in [S::Left, S::Automatic, S::Right] {
                    assert_eq!(
                        input.curve.exact_evaluate(&u, D::Second, side),
                        result.exact_evaluate(&u, D::Second, side),
                        "{}",
                        input.name
                    );
                }
            }
        }
        count += 1;
    }
    assert_eq!(count, 142);
}

#[test]
fn every_tensor_control_matches_independent_coefficient_equations() {
    let mut count = 0;
    for row in include_str!("../../fixtures/degree-elevation-surfaces.tsv").lines() {
        let (input, expected) = row.split_once(" | ").unwrap();
        let input = surfaces::parse(input);
        let (flags, result) = surfaces::execute(&input);
        assert_eq!(
            surfaces::encode(&input.name, &flags, &result, true),
            expected,
            "{}",
            input.name
        );
        assert_eq!(result.domain(), input.surface.domain());
        if std::env::var_os("RUSTY_VERIFY_ALL_DEGREE_ORACLES").is_some() {
            assert_eq!(
                result,
                reference::elevated_surface(&input.surface, result.degrees()),
                "{}",
                input.name
            );
        }
        if result.degrees().iter().all(|&p| p <= 4) {
            let [u, v] = result.degrees();
            let reversed_order = input
                .surface
                .elevated(input.surface.degrees()[0], v)
                .unwrap()
                .elevated(u, v)
                .unwrap();
            assert_eq!(result, reversed_order, "{}", input.name);
            assert_eq!(
                result.exchanged_uv(),
                input.surface.exchanged_uv().elevated(v, u).unwrap()
            );
        }
        count += 1;
    }
    assert_eq!(count, 125);
}

#[test]
fn exact_extremes_and_inactive_unclamped_support_keep_the_domain() {
    let huge = R::from_integer(BigInt::from(1) << 4096);
    let tiny = R::new(1.into(), BigInt::from(1) << 1100);
    let start = &huge + r(1) / r(3);
    let axis = ExactKnotVector::new(
        2,
        (0..8).map(|i| &start + r(i) * &tiny).collect(),
        vec![1; 8],
    )
    .unwrap();
    let controls = (0..5)
        .map(|i| {
            [
                r(i),
                r(-3 * i),
                huge.clone(),
                if i % 2 == 0 { tiny.clone() } else { r(1) },
            ]
        })
        .collect();
    let original = ExactBSplineCurve3::from_homogeneous(axis, controls).unwrap();
    let elevated = original.elevated(3).unwrap();
    assert_eq!(elevated.domain(), original.domain());
    assert_eq!(
        elevated.knots(),
        &(1..7).map(|i| &start + r(i) * &tiny).collect::<Vec<_>>()
    );
    assert_eq!(elevated.multiplicities(), &[2; 6]);
    assert_eq!(elevated.homogeneous_poles().len(), 8);
    let u = &start + r(7) * &tiny / r(2);
    assert_eq!(
        original.exact_evaluate(&u, D::Second, S::Automatic),
        elevated.exact_evaluate(&u, D::Second, S::Automatic)
    );
    assert_eq!(original.elevated(4).unwrap(), elevated.elevated(4).unwrap());
    let clamped = ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new(2, vec![start.clone(), &start + tiny.clone()], vec![3, 3]).unwrap(),
        vec![
            [huge.clone(), r(0), r(1), tiny.clone()],
            [r(-1), huge.clone(), r(3), r(2)],
            [r(3), r(2), huge.clone(), huge],
        ],
    )
    .unwrap();
    let bezier = clamped.bezier_arcs().unwrap().remove(0);
    assert_eq!(
        bezier.to_bspline().elevated(5).unwrap(),
        bezier.elevated(5).unwrap().to_bspline()
    );
}

fn linear_axis(count: usize) -> ExactKnotVector {
    let mut mults = vec![1; count];
    mults[0] = 2;
    mults[count - 1] = 2;
    ExactKnotVector::new(1, (0..count).map(|i| r(i as i32)).collect(), mults).unwrap()
}

#[test]
fn degree_and_cartesian_limits_reject_atomically() {
    let pole = [r(0), r(0), r(0), r(1)];
    let curve =
        ExactBSplineCurve3::from_homogeneous(linear_axis(4096), vec![pole.clone(); 4096]).unwrap();
    let saved = curve.clone();
    assert!(matches!(curve.elevated(2), Err(Error::LimitExceeded(_))));
    for degree in [0, 26, usize::MAX] {
        assert!(matches!(
            curve.elevated(degree),
            Err(Error::InvalidSpline(_))
        ));
    }
    assert_eq!(curve.elevated(1).unwrap(), saved);
    assert_eq!(curve, saved);
    let surface =
        ExactBSplineSurface3::from_homogeneous(linear_axis(64), linear_axis(64), vec![pole; 4096])
            .unwrap();
    let saved = surface.clone();
    for degrees in [(2, 1), (1, 2), (2, 2)] {
        assert!(matches!(
            surface.elevated(degrees.0, degrees.1),
            Err(Error::LimitExceeded(_))
        ));
    }
    for degrees in [(0, 1), (1, 0), (26, 1), (1, 26), (usize::MAX, usize::MAX)] {
        assert!(matches!(
            surface.elevated(degrees.0, degrees.1),
            Err(Error::InvalidSpline(_))
        ));
    }
    assert_eq!(surface.elevated(1, 1).unwrap(), saved);
    assert_eq!(surface, saved);
}
