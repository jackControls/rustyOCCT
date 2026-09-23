//! Independent complete degree-elevation control reconstruction. Fraction-free
//! Cox equations determine every working control, including inactive fields;
//! no production averaging, insertion, extraction or evaluation is called.
use super::knot_reference;
use num_rational::BigRational as R;
use rusty_occt::{ExactBSplineCurve3, ExactBSplineSurface3, ExactKnotVector};

fn integer(n: usize) -> R {
    R::from_integer(n.into())
}

fn grouped(flat: &[R]) -> (Vec<R>, Vec<usize>) {
    let mut knots = Vec::new();
    let mut mults = Vec::new();
    for k in flat {
        if knots.last() == Some(k) {
            *mults.last_mut().unwrap() += 1;
        } else {
            knots.push(k.clone());
            mults.push(1);
        }
    }
    (knots, mults)
}

/// Working axes must fit the test helper's validated 4096-pole basis type.
/// Production preflights only its exposed output; this is an oracle limit.
pub fn elevated_rows(original: &[ExactBSplineCurve3], degree: usize) -> Vec<ExactBSplineCurve3> {
    assert!(!original.is_empty());
    let old = original[0].knot_vector();
    assert!(original.iter().all(|r| r.knot_vector() == old));
    assert!((old.degree()..=25).contains(&degree));
    if old.degree() == degree {
        return original.to_vec();
    }
    let delta = degree - old.degree();
    // When there are many transverse fields, solve all source unit columns
    // once. Their complete coefficient identities prove the linear map for
    // every field, independent of its coordinate size. This does not use the
    // production control map or infer equality from selected geometric rows.
    let symbolic = original.len() > 1 && 4 * original.len() > old.pole_count();
    let mut controls: Vec<Vec<[R; 4]>> = if symbolic {
        (0..old.pole_count().div_ceil(4))
            .map(|block| {
                (0..old.pole_count())
                    .map(|i| std::array::from_fn(|c| integer(usize::from(i == 4 * block + c))))
                    .collect()
            })
            .collect()
    } else {
        original
            .iter()
            .map(|r| r.homogeneous_poles().to_vec())
            .collect()
    };
    let (working, target, final_axis, skip) = if old.is_periodic() {
        let target = ExactKnotVector::new_periodic(
            degree,
            old.knots().to_vec(),
            old.multiplicities().iter().map(|m| m + delta).collect(),
        )
        .unwrap();
        (old.clone(), target.clone(), target, 0)
    } else {
        let p = old.degree();
        let first = old
            .knots()
            .iter()
            .position(|k| k == &old.domain()[0])
            .unwrap();
        let last = old
            .knots()
            .iter()
            .position(|k| k == &old.domain()[1])
            .unwrap();
        let (prefix, suffix) = (
            p + 1 - old.multiplicities()[0],
            p + 1 - old.multiplicities().last().unwrap(),
        );
        let mut mults = old.multiplicities().to_vec();
        mults[0] = p + 1;
        *mults.last_mut().unwrap() = p + 1;
        let working = ExactKnotVector::new(p, old.knots().to_vec(), mults.clone()).unwrap();
        for m in &mut mults {
            *m += delta;
        }
        let target = ExactKnotVector::new(degree, old.knots().to_vec(), mults.clone()).unwrap();
        let flat: Vec<_> = old
            .knots()
            .iter()
            .zip(&mults)
            .flat_map(|(k, &m)| std::iter::repeat_n(k.clone(), m))
            .collect();
        let skip = prefix + delta * first;
        let tail = suffix + delta * (old.knots().len() - 1 - last);
        let (knots, mults) = grouped(&flat[skip..flat.len() - tail]);
        let final_axis = ExactKnotVector::new(degree, knots, mults).unwrap();
        assert!((0..skip).all(|i| flat[i + degree + 1] <= old.domain()[0]));
        assert!(
            (target.pole_count() - tail..target.pole_count()).all(|i| flat[i] >= old.domain()[1])
        );
        let zero: [R; 4] = std::array::from_fn(|_| integer(0));
        for row in &mut controls {
            row.splice(0..0, std::iter::repeat_n(zero.clone(), prefix));
            row.extend(std::iter::repeat_n(zero.clone(), suffix));
        }
        (working, target, final_axis, skip)
    };
    assert_eq!(old.domain(), final_axis.domain());
    let recovered = knot_reference::recover_controls(
        &working,
        &controls.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        &target,
    )
    .expect("complete degree elevation coefficient system is consistent");
    let recovered: Vec<Vec<[R; 4]>> = if symbolic {
        let transform: Vec<Vec<_>> = (skip..skip + final_axis.pole_count())
            .map(|j| {
                (0..old.pole_count())
                    .map(|i| (i, &recovered[i / 4][j][i % 4]))
                    .filter(|(_, x)| **x != integer(0))
                    .collect()
            })
            .collect();
        original
            .iter()
            .map(|curve| {
                transform
                    .iter()
                    .map(|row| {
                        std::array::from_fn(|c| {
                            row.iter()
                                .map(|(i, x)| *x * &curve.homogeneous_poles()[*i][c])
                                .sum()
                        })
                    })
                    .collect()
            })
            .collect()
    } else {
        recovered
            .into_iter()
            .map(|row| row[skip..skip + final_axis.pole_count()].to_vec())
            .collect()
    };
    recovered
        .into_iter()
        .map(|row| ExactBSplineCurve3::from_homogeneous(final_axis.clone(), row).unwrap())
        .collect()
}

pub fn elevated(curve: &ExactBSplineCurve3, degree: usize) -> ExactBSplineCurve3 {
    elevated_rows(std::slice::from_ref(curve), degree).remove(0)
}

pub fn elevated_surface(
    surface: &ExactBSplineSurface3,
    degrees: [usize; 2],
) -> ExactBSplineSurface3 {
    let mut result = surface.clone();
    for axis in 0..2 {
        let counts = result.pole_counts();
        let mut axes = [result.u_knots().clone(), result.v_knots().clone()];
        if degrees[axis] == axes[axis].degree() {
            continue;
        }
        let rows: Vec<_> = (0..counts[1 - axis])
            .map(|fixed| {
                ExactBSplineCurve3::from_homogeneous(
                    axes[axis].clone(),
                    (0..counts[axis])
                        .map(|i| {
                            result.homogeneous_poles()[if axis == 0 {
                                i * counts[1] + fixed
                            } else {
                                fixed * counts[1] + i
                            }]
                            .clone()
                        })
                        .collect(),
                )
                .unwrap()
            })
            .collect();
        let rows = elevated_rows(&rows, degrees[axis]);
        axes[axis] = rows[0].knot_vector().clone();
        let controls = (0..axes[0].pole_count())
            .flat_map(|i| {
                let rows = &rows;
                (0..axes[1].pole_count()).map(move |j| {
                    if axis == 0 {
                        rows[j].homogeneous_poles()[i].clone()
                    } else {
                        rows[i].homogeneous_poles()[j].clone()
                    }
                })
            })
            .collect();
        let [u, v] = axes;
        result = ExactBSplineSurface3::from_homogeneous(u, v, controls).unwrap();
    }
    result
}
