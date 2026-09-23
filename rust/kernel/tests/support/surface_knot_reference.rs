//! Independent full-support coefficient checks for every transverse control.
//! The coefficient oracle uses Cox polynomials and fraction-free elimination,
//! independently of the kernel's local Boehm/inverse-Boehm implementation.
use super::knot_reference as reference;
use num_rational::BigRational as R;
use rusty_occt::{ExactBSplineCurve3, ExactBSplineSurface3};

pub fn curves(surface: &ExactBSplineSurface3, axis: usize) -> Vec<ExactBSplineCurve3> {
    let [nu, nv] = surface.pole_counts();
    let basis = if axis == 0 {
        surface.u_knots()
    } else {
        surface.v_knots()
    };
    (0..[nv, nu][axis])
        .map(|j| {
            ExactBSplineCurve3::from_homogeneous(
                basis.clone(),
                (0..[nu, nv][axis])
                    .map(|i| {
                        surface.homogeneous_poles()[if axis == 0 { i * nv + j } else { j * nv + i }]
                            .clone()
                    })
                    .collect(),
            )
            .unwrap()
        })
        .collect()
}

pub fn equal_axis(original: &ExactBSplineSurface3, candidate: &ExactBSplineSurface3, axis: usize) {
    assert_eq!(
        if axis == 0 {
            original.v_knots()
        } else {
            original.u_knots()
        },
        if axis == 0 {
            candidate.v_knots()
        } else {
            candidate.u_knots()
        }
    );
    let before = curves(original, axis);
    let after = curves(candidate, axis);
    assert_eq!(before.len(), after.len());
    assert!(
        reference::equal_many(&before, &after),
        "complete tensor coefficient identity"
    );
}

pub fn check_removal(
    original: &ExactBSplineSurface3,
    axis: usize,
    u: &R,
    target: usize,
) -> Option<ExactBSplineSurface3> {
    let result = if axis == 0 {
        original.remove_u_knot(u, target)
    } else {
        original.remove_v_knot(u, target)
    }
    .unwrap();
    let source = curves(original, axis);
    let basis = reference::removal_basis(&source[0], u, target);
    let answers = reference::recover_many(&source, &basis).and_then(|rows| {
        rows.into_iter()
            .map(|c| {
                if c.iter().any(|p| p[3] <= R::from_integer(0.into())) {
                    None
                } else {
                    Some(ExactBSplineCurve3::from_homogeneous(basis.clone(), c).unwrap())
                }
            })
            .collect::<Option<Vec<_>>>()
    });
    assert_eq!(
        result.is_some(),
        answers.is_some(),
        "exact surface removal feasibility"
    );
    if let (Some(candidate), Some(answers)) = (&result, answers) {
        assert_eq!(curves(candidate, axis), answers);
    }
    result
}
