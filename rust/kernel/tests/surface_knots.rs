#[allow(dead_code)]
#[path = "support/knot_reference.rs"]
mod knot_reference;
#[path = "support/surface_knot_protocol.rs"]
mod protocol;
#[path = "support/surface_knot_reference.rs"]
mod reference;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::curve::{DerivativeOrder as D, KnotSide as S};
use rusty_occt::surface::BezierPatchExtractionOptions;
use rusty_occt::{
    BSplineSurface3, BezierSurface3, Error, ExactBSplineSurface3, ExactKnotVector, KnotVector,
    Point3,
};
fn q(n: i32, d: i32) -> R {
    R::new(n.into(), d.into())
}
fn patch() -> ExactBSplineSurface3 {
    BezierSurface3::new(
        2,
        3,
        (0..3)
            .flat_map(|i| (0..4).map(move |j| Point3::new(i as f64, j as f64, (i * j) as f64)))
            .collect(),
        Some((0..12).map(|i| (1 + i % 5) as f64).collect()),
    )
    .unwrap()
    .to_exact()
    .to_bspline()
}
#[test]
fn full_grids_and_removal_decisions_match_independent_equations() {
    let mut count = 0;
    for row in include_str!("../../fixtures/surface-knots.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let (input, expected) = row.split_once(" | ").unwrap();
        let input = protocol::parse(input);
        let (flags, surface) = protocol::execute(&input);
        assert_eq!(
            protocol::encode(&input.name, &flags, &surface, true),
            expected,
            "{}",
            input.name
        );
        // All fixtures come from the complete Python coefficient equations.
        // Selected ordinary cases and every fuzz edit also use a separately
        // implemented fraction-free Rust polynomial/feasibility oracle.
        if (count % 97 == 0 && input.surface.degrees().iter().all(|&d| d <= 8))
            || std::env::var_os("RUSTY_VERIFY_ALL_SURFACE_KNOT_ORACLES").is_some()
        {
            let mut current = input.surface.clone();
            for (op, axis, u, m) in &input.operations {
                if *op == 'I' {
                    let next = if *axis == 0 {
                        current.insert_u_knot(u, *m)
                    } else {
                        current.insert_v_knot(u, *m)
                    }
                    .unwrap();
                    reference::equal_axis(&current, &next, *axis);
                    current = next;
                } else if let Some(next) = reference::check_removal(&current, *axis, u, *m) {
                    current = next;
                }
            }
            assert_eq!(current, surface);
        }
        count += 1;
    }
    assert_eq!(count, 1756);
}
#[test]
fn rational_batches_commute_and_preserve_complete_control_grid() {
    let original = patch();
    let u = q(1, 3);
    let v = q(2, 7);
    let close = &u + R::new(1.into(), BigInt::from(1) << 2048);
    let ur = [(u.clone(), 2), (close.clone(), 1), (u.clone(), 0)];
    let vr = [(v.clone(), 3), (v.clone(), 1)];
    let first = original.refined_u(&ur).unwrap();
    reference::equal_axis(&original, &first, 0);
    let edited = first.refined_v(&vr).unwrap();
    reference::equal_axis(&first, &edited, 1);
    assert_eq!(edited, original.refined(&ur, &vr).unwrap());
    assert_eq!(
        edited,
        original.refined_v(&vr).unwrap().refined_u(&ur).unwrap()
    );
    assert_eq!(
        edited.exchanged_uv(),
        original.exchanged_uv().refined(&vr, &ur).unwrap()
    );
    assert_eq!(edited.exchanged_uv().exchanged_uv(), edited);
    let restored = edited
        .remove_v_knot(&v, 0)
        .unwrap()
        .unwrap()
        .remove_u_knot(&close, 0)
        .unwrap()
        .unwrap()
        .remove_u_knot(&u, 0)
        .unwrap()
        .unwrap();
    assert_eq!(restored, original);
    for (a, b) in [
        (q(0, 1), v.clone()),
        (u.clone(), v.clone()),
        (close, v),
        (q(1, 1), q(1, 1)),
    ] {
        assert_eq!(
            original.exact_evaluate(&a, &b, D::Second, [S::Automatic; 2]),
            edited.exact_evaluate(&a, &b, D::Second, [S::Automatic; 2])
        );
    }
}
#[test]
fn exact_patches_isocurves_and_partials_keep_original_parameter_units() {
    let original = patch();
    let u = q(1, 3);
    let v = q(2, 7);
    let edited = original
        .refined(&[(u.clone(), 2)], &[(v.clone(), 3)])
        .unwrap();
    for p in edited.bezier_patches().unwrap() {
        let [[a, b], [c, d]] = p.domain();
        let expected = original.bezier_patches_in(a, b, c, d).unwrap().remove(0);
        assert_eq!(p, expected);
        assert_eq!(p.to_bspline().bezier_patches().unwrap(), vec![p.clone()]);
        let (x, y) = ((a + b) / q(2, 1), (c + d) / q(2, 1));
        let jet = edited
            .exact_evaluate(&x, &y, D::Second, [S::Automatic; 2])
            .unwrap();
        assert_eq!(jet, p.exact_evaluate(&x, &y, D::Second).unwrap());
        let iso = edited.u_iso(&x).unwrap();
        let along = iso.exact_evaluate(&y, D::Second, S::Automatic).unwrap();
        assert_eq!(along.position(), jet.position());
        assert_eq!(along.derivative(1), jet.derivative(0, 1));
        assert_eq!(along.derivative(2), jet.derivative(0, 2));
        assert_eq!(
            edited
                .v_iso(&y)
                .unwrap()
                .exact_evaluate(&x, D::Second, S::Automatic)
                .unwrap()
                .derivative(2),
            jet.derivative(2, 0)
        );
        assert_eq!(iso.bezier_arcs_in(c, d).unwrap()[0], p.u_iso(&x).unwrap());
    }
    let end = &u + R::new(1.into(), BigInt::from(1) << 1100);
    let strip = edited.bezier_patches_in(&u, &end, &q(0, 1), &v).unwrap();
    assert_eq!(strip[0].domain(), &[[u.clone(), end], [q(0, 1), v.clone()]]);
    assert_eq!(
        strip[0].exact_evaluate(&u, &v, D::Second).unwrap(),
        edited
            .exact_evaluate(&u, &v, D::Second, [S::Automatic; 2])
            .unwrap()
    );
}
#[test]
fn periodic_origin_changes_are_independent_and_support_rational_multiturns() {
    let axis =
        ExactKnotVector::new_periodic(2, (0..5).map(|n| q(n, 1)).collect(), vec![1; 5]).unwrap();
    let other = ExactKnotVector::new(2, vec![q(0, 1), q(1, 1)], vec![3, 3]).unwrap();
    let controls = (0..4)
        .flat_map(|_| {
            (0..3).map(|j| {
                [
                    q(j * (j + 1), 1),
                    q(j * j * (j + 1), 1),
                    q(0, 1),
                    q(j + 1, 1),
                ]
            })
        })
        .collect();
    let original = ExactBSplineSurface3::from_homogeneous(axis, other, controls).unwrap();
    let seam = original.refined_u(&[(q(4, 1), 2), (q(0, 1), 1)]).unwrap();
    assert_eq!(seam, original.insert_u_knot(&q(0, 1), 2).unwrap());
    assert_eq!(
        seam.remove_u_knot(&q(4, 1), 1).unwrap(),
        Some(original.clone())
    );
    let shifted = reference::check_removal(&original, 0, &q(0, 1), 0).unwrap();
    assert_eq!(shifted.domain()[0], [q(1, 1), q(5, 1)]);
    assert_eq!(shifted.v_knots(), original.v_knots());
    assert_eq!(
        original.remove_u_knot(&q(4, 1), 0).unwrap(),
        Some(shifted.clone())
    );
    assert_eq!(
        original.exchanged_uv().remove_v_knot(&q(0, 1), 0).unwrap(),
        Some(shifted.exchanged_uv())
    );
    let far = R::from_integer(BigInt::from(1) << 2048);
    let pieces = shifted
        .bezier_patches_in(&far, &(&far + q(8, 1)), &q(0, 1), &q(1, 1))
        .unwrap();
    assert_eq!(pieces.len(), 7);
    assert_eq!(pieces[0].domain()[0][0], far);
    for p in pieces {
        let u = (&p.domain()[0][0] + &p.domain()[0][1]) / q(2, 1);
        let v = q(1, 3);
        assert_eq!(
            original
                .exact_evaluate(&u, &v, D::Second, [S::Automatic; 2])
                .unwrap(),
            p.exact_evaluate(&u, &v, D::Second).unwrap()
        );
    }
    assert!(matches!(
        shifted.bezier_patches_with_options(
            &q(1, 1),
            &q(9, 1),
            &q(0, 1),
            &q(1, 1),
            BezierPatchExtractionOptions {
                max_patches: 5,
                max_controls: 4096
            }
        ),
        Err(Error::ComputationLimit(_))
    ));
}
#[test]
fn one_nonremovable_transverse_row_rejects_the_whole_surface() {
    let u = ExactKnotVector::new(1, vec![q(0, 1), q(1, 1), q(2, 1)], vec![2, 1, 2]).unwrap();
    let v = ExactKnotVector::new(1, vec![q(0, 1), q(1, 1)], vec![2, 2]).unwrap();
    let original = ExactBSplineSurface3::from_homogeneous(
        u,
        v,
        (0..3)
            .flat_map(|i| {
                (0..2).map(move |j| [q(i, 1), q(j, 1), q(i32::from(i == 1 && j == 1), 1), q(1, 1)])
            })
            .collect(),
    )
    .unwrap();
    assert!(reference::check_removal(&original, 0, &q(1, 1), 0).is_none());
    assert!(original
        .exact_evaluate(&q(1, 1), &q(1, 1), D::First, [S::Automatic; 2])
        .is_err());
    let left = original
        .exact_evaluate(&q(1, 1), &q(1, 1), D::First, [S::Left, S::Left])
        .unwrap();
    let right = original
        .exact_evaluate(&q(1, 1), &q(1, 1), D::First, [S::Right, S::Left])
        .unwrap();
    assert_eq!(left.derivative(1, 0).unwrap()[2], q(1, 1));
    assert_eq!(right.derivative(1, 0).unwrap()[2], q(-1, 1));
    let u = ExactKnotVector::new(2, vec![q(0, 1), q(1, 2), q(1, 1)], vec![3, 1, 3]).unwrap();
    let v = original.v_knots().clone();
    let invalid_inverse = ExactBSplineSurface3::from_homogeneous(
        u,
        v,
        [q(1, 1), q(3, 8), q(3, 8), q(1, 1)]
            .iter()
            .flat_map(|w| (0..2).map(move |_| [q(0, 1), q(0, 1), q(0, 1), w.clone()]))
            .collect(),
    )
    .unwrap();
    assert!(reference::check_removal(&invalid_inverse, 0, &q(1, 2), 0).is_none());
}
#[test]
fn invalid_batches_and_cartesian_limits_are_atomic() {
    let original = patch();
    let bad = R::new_raw(1.into(), 0.into());
    for (a, b) in [
        (vec![], vec![(bad.clone(), 0)]),
        (vec![(q(1, 3), 1)], vec![(q(2, 1), 0)]),
        (vec![(q(1, 3), 3)], vec![]),
    ] {
        assert!(original.refined(&a, &b).is_err());
    }
    assert!(original
        .refined(&vec![(q(1, 3), 0); 2048], &vec![(q(1, 3), 0); 2049])
        .is_err());
    assert!(original.remove_u_knot(&bad, 0).is_err());
    assert!(original.u_iso(&bad).is_err());
    assert!(original
        .exact_evaluate(&bad, &q(0, 1), D::Position, [S::Automatic; 2])
        .is_err());
    assert!(original
        .evaluate(f64::NAN, 0., D::Position, [S::Automatic; 2])
        .is_err());
    assert!(original
        .bezier_patches_in(&bad, &q(1, 1), &q(0, 1), &q(1, 1))
        .is_err());
    assert!(original
        .bezier_patches_with_options(
            &q(0, 1),
            &q(1, 1),
            &q(0, 1),
            &q(1, 1),
            BezierPatchExtractionOptions {
                max_patches: 1,
                max_controls: 11
            }
        )
        .is_err());
    let axis = ExactKnotVector::new(
        1,
        (0..64).map(|i| q(i, 1)).collect(),
        std::iter::once(2)
            .chain(std::iter::repeat_n(1, 62))
            .chain(std::iter::once(2))
            .collect(),
    )
    .unwrap();
    let large = ExactBSplineSurface3::from_homogeneous(
        axis.clone(),
        axis,
        vec![[q(0, 1), q(0, 1), q(0, 1), q(1, 1)]; 4096],
    )
    .unwrap();
    assert!(matches!(
        large.insert_u_knot(&q(1, 2), 1),
        Err(Error::LimitExceeded(_))
    ));
    let mut controls = original.homogeneous_poles().to_vec();
    controls[0][3] = bad;
    assert!(ExactBSplineSurface3::from_homogeneous(
        original.u_knots().clone(),
        original.v_knots().clone(),
        controls
    )
    .is_err());
    assert_eq!(original, patch());
}
#[test]
fn binary64_conversion_preserves_unrepresentable_exact_partials() {
    let tiny = f64::from_bits(1);
    let axis = KnotVector::new(1, vec![0., tiny], vec![2, 2]).unwrap();
    let source = BSplineSurface3::new(
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
    let exact = source.to_exact();
    let u = R::from_float(tiny).unwrap() / q(2, 1);
    let edited = exact.refined(&[(u.clone(), 1)], &[(u.clone(), 1)]).unwrap();
    let jet = edited
        .exact_evaluate(&u, &u, D::Second, [S::Automatic; 2])
        .unwrap();
    assert_eq!(jet.position().coordinates(), &[q(1, 2), q(1, 2), q(1, 4)]);
    assert!(jet.enclosed().is_err());
    assert_eq!(
        exact
            .exact_evaluate(&u, &u, D::Second, [S::Automatic; 2])
            .unwrap(),
        jet
    );
}

#[allow(dead_code)]
#[path = "support/bezier_reference.rs"]
mod bezier_reference;
#[allow(dead_code)]
#[path = "support/surface_editing_reference.rs"]
mod tensor_reference;

#[test]
fn retained_degree25_periodic_surface_timeout_keeps_complete_oracles() {
    // Minimized structural input: a degree-25 periodic tensor, redundant in
    // one direction, with nonconstant weighted geometry in the other.
    let axis =
        ExactKnotVector::new_periodic(25, (0..29).map(|i| q(i, 1)).collect(), vec![1; 29]).unwrap();
    let controls = (0..28)
        .flat_map(|_| {
            (0..28).map(|j| {
                let w = q(1 + j % 5, 1);
                [q(j, 1) * &w, q(j * j, 1) * &w, w.clone(), w]
            })
        })
        .collect();
    let original = ExactBSplineSurface3::from_homogeneous(axis.clone(), axis, controls).unwrap();
    let edited = reference::check_removal(&original, 0, &q(0, 1), 0).unwrap();
    assert_eq!(
        original.exchanged_uv().remove_v_knot(&q(0, 1), 0).unwrap(),
        Some(edited.exchanged_uv())
    );
    let domain = [[q(1, 1), q(2, 1)], [q(0, 1), q(1, 1)]];
    let expected = tensor_reference::exact_patch(
        &original,
        domain.clone(),
        knot_reference::basis_coefficients,
    );
    let p = edited
        .bezier_patches_in(&domain[0][0], &domain[0][1], &domain[1][0], &domain[1][1])
        .unwrap()
        .remove(0);
    tensor_reference::check(&p, &expected, [&q(85, 257), &q(85, 257)]);
    let (u, v) = (q(1, 1) + q(85, 257), q(85, 257));
    let jet = edited
        .exact_evaluate(&u, &v, D::Second, [S::Automatic; 2])
        .unwrap();
    assert_eq!(jet, p.exact_evaluate(&u, &v, D::Second).unwrap());
    let offset = R::from_integer(BigInt::from(1) << 2048) * q(28, 1);
    assert_eq!(
        jet,
        edited
            .exact_evaluate(&(u + &offset), &(v + offset), D::Second, [S::Automatic; 2])
            .unwrap()
    );
}

#[test]
fn edited_surface_isocurves_feed_certified_analytic_intersections() {
    use rusty_occt::intersection::{
        exact_spline_cylinder, exact_spline_plane, exact_spline_sphere, Cylinder3, Plane3, Sphere3,
    };
    use rusty_occt::Vec3;
    use std::cmp::Ordering::Equal;
    let original = BezierSurface3::new(
        1,
        1,
        vec![
            Point3::new(-2., -1., 0.),
            Point3::new(-2., 1., 0.),
            Point3::new(2., -1., 0.),
            Point3::new(2., 1., 0.),
        ],
        None,
    )
    .unwrap()
    .to_exact()
    .to_bspline();
    let edited = original.refined(&[(q(1, 3), 1)], &[(q(2, 7), 1)]).unwrap();
    let curve = edited.v_iso(&q(1, 2)).unwrap();
    let plane = Plane3::through_points(
        Point3::ORIGIN,
        Point3::new(0., 1., 0.),
        Point3::new(0., 0., 1.),
    )
    .unwrap();
    let hit = exact_spline_plane(&curve, &plane).unwrap();
    assert_eq!(hit.points().len(), 1);
    assert!(hit.overlaps().is_empty());
    assert_eq!(hit.points()[0].compare_parameter(&q(1, 2)).unwrap(), Equal);
    let sphere = Sphere3::new(Point3::ORIGIN, 1.).unwrap();
    let cylinder = Cylinder3::new(Point3::ORIGIN, Vec3::new(0., 0., 1.), 1.).unwrap();
    for hits in [
        exact_spline_sphere(&curve, &sphere).unwrap(),
        exact_spline_cylinder(&curve, &cylinder).unwrap(),
    ] {
        assert_eq!(hits.points().len(), 2);
        assert!(hits.overlaps().is_empty());
        for (p, t) in hits.points().iter().zip([q(1, 4), q(3, 4)]) {
            assert_eq!(p.compare_parameter(&t).unwrap(), Equal);
        }
    }
}

#[test]
fn retained_degree25_unclamped_batch_timeout_checks_every_raw_support_coefficient() {
    for u_target in [2, 25] {
        check_unclamped_batch(u_target);
    }
}

fn check_unclamped_batch(u_target: usize) {
    let axis = ExactKnotVector::new(25, (0..54).map(|i| q(i, 1)).collect(), vec![1; 54]).unwrap();
    let controls = (0..784)
        .map(|i| {
            let byte = |c: usize| ((37 * (4 * i + c) + 1) % 256) as u8;
            let w = q(1 + i32::from(byte(3) % 8), 1);
            [
                q(byte(0) as i8 as i32, 1) * &w,
                q(byte(1) as i8 as i32, 1) * &w,
                q(byte(2) as i8 as i32, 1) * &w,
                w,
            ]
        })
        .collect();
    let original = ExactBSplineSurface3::from_homogeneous(axis.clone(), axis, controls).unwrap();
    let cut = q(25, 1) + q(85, 257);
    let first = original.insert_u_knot(&cut, u_target).unwrap();
    reference::equal_axis(&original, &first, 0);
    let both = first.insert_v_knot(&cut, 25).unwrap();
    reference::equal_axis(&first, &both, 1);
    assert_eq!(
        both,
        original
            .refined(&[(cut.clone(), u_target)], &[(cut.clone(), 25)])
            .unwrap()
    );
    assert_eq!(
        both.exchanged_uv(),
        original
            .exchanged_uv()
            .refined(&[(cut.clone(), 25)], &[(cut.clone(), u_target)])
            .unwrap()
    );
    let domain = [[q(25, 1), cut.clone()], [q(25, 1), cut.clone()]];
    let expected =
        tensor_reference::exact_patch(&original, domain, knot_reference::basis_coefficients);
    let p = both
        .bezier_patches_in(&q(25, 1), &cut, &q(25, 1), &cut)
        .unwrap()
        .remove(0);
    tensor_reference::check(&p, &expected, [&q(85, 257), &q(85, 257)]);
}

#[test]
fn coefficient_checker_rejects_every_control_field_and_large_signed_residuals() {
    let original = patch();
    let edited = original.insert_u_knot(&q(1, 3), 2).unwrap();
    let before = reference::curves(&original, 0);
    let after = reference::curves(&edited, 0);
    assert!(knot_reference::equal_many(&before, &after));
    for pole in 0..edited.homogeneous_poles().len() {
        for c in 0..4 {
            let mut controls = edited.homogeneous_poles().to_vec();
            controls[pole][c] += q(1, 7);
            let bad = ExactBSplineSurface3::from_homogeneous(
                edited.u_knots().clone(),
                edited.v_knots().clone(),
                controls,
            )
            .unwrap();
            assert!(
                !knot_reference::equal_many(&before, &reference::curves(&bad, 0)),
                "undetected grid field {pole}/{c}"
            );
        }
    }
    let huge = R::from_integer(BigInt::from(1) << 4096);
    let tiny = R::new(1.into(), BigInt::from(1) << 4096);
    // Large and tiny signed residuals must not cancel across distinct fields.
    for perturb in [huge, tiny] {
        let mut controls = edited.homogeneous_poles().to_vec();
        controls[1][0] += &perturb;
        controls[1][1] -= &perturb;
        let bad = ExactBSplineSurface3::from_homogeneous(
            edited.u_knots().clone(),
            edited.v_knots().clone(),
            controls,
        )
        .unwrap();
        assert!(!knot_reference::equal_many(
            &before,
            &reference::curves(&bad, 0)
        ));
    }
}
