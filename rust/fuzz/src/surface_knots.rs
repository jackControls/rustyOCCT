//! Complete coefficient identities and exact removal decisions for tensor
//! surfaces, independently periodic axes and rational operation sequences.
use crate::{byte, knot_reference};
#[path = "../../kernel/tests/support/surface_knot_reference.rs"]
mod reference;
use crate::surface_reference as tensor;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::curve::{DerivativeOrder as D, KnotSide as S};
use rusty_occt::{BSplineSurface3, Error, ExactBSplineSurface3, KnotVector, Point3};
pub fn check_surface_knots(data: &[u8]) {
    check(data, |_| {});
}
pub fn profile_surface_knots(
    data: &[u8],
) -> std::collections::BTreeMap<&'static str, std::time::Duration> {
    let mut times = std::collections::BTreeMap::new();
    let mut at = std::time::Instant::now();
    check(data, |label| {
        *times.entry(label).or_default() += at.elapsed();
        at = std::time::Instant::now();
    });
    times
}
fn insert(
    s: &ExactBSplineSurface3,
    axis: usize,
    u: &R,
    m: usize,
    mark: &mut impl FnMut(&'static str),
) -> ExactBSplineSurface3 {
    let next = if axis == 0 {
        s.insert_u_knot(u, m)
    } else {
        s.insert_v_knot(u, m)
    }
    .unwrap();
    mark("kernel refinement");
    reference::equal_axis(s, &next, axis);
    mark("complete coefficient identity");
    next
}
fn check(data: &[u8], mut mark: impl FnMut(&'static str)) {
    let mode = byte(data, 0) % 4;
    let bound = match mode {
        0 => 2,
        1 => 5,
        _ => 25,
    };
    let degrees: [usize; 2] = std::array::from_fn(|i| 1 + byte(data, 1 + i) as usize % bound);
    let kinds: [u8; 2] = std::array::from_fn(|i| if mode == 3 { 1 } else { byte(data, 3 + i) % 3 });
    let axis = byte(data, 5) as usize % 2;
    let op = byte(data, 6) % 5;
    let axes: [KnotVector; 2] = std::array::from_fn(|i| {
        let d = degrees[i];
        let scale = match byte(data, 10 + i) % 8 {
            0 => f64::from_bits(1),
            1 => 2_f64.powi(-500),
            2 => 2_f64.powi(500),
            _ => 1.,
        };
        let (k, m) = if mode == 3 {
            (
                (0..d + 4).map(|i| i as f64).collect::<Vec<_>>(),
                vec![1; d + 4],
            )
        } else {
            match kinds[i] {
                0 => (vec![0., 1., 3.], vec![d + 1, 1, d + 1]),
                1 => (vec![0., 1., 3.], vec![1, d, 1]),
                _ => (
                    (0..2 * d + 4).map(|i| i as f64).collect(),
                    vec![1; 2 * d + 4],
                ),
            }
        };
        let k = k.into_iter().map(|x| x * scale).collect();
        if kinds[i] == 1 {
            KnotVector::new_periodic(d, k, m)
        } else {
            KnotVector::new(d, k, m)
        }
        .unwrap()
    });
    let [nu, nv] = axes.each_ref().map(|a| a.pole_count());
    let scale = 2_f64.powi((byte(data, 9) as i32 - 128) * 7);
    let scalar = |i: usize, j: usize, c: usize| {
        if mode == 0 {
            f64::from_bits(u64::from_le_bytes(std::array::from_fn(|k| {
                byte(data, 16 + 8 * (4 * (i * nv + j) + c) + k)
            })))
        } else if mode == 3 {
            let t = if axis == 0 { j } else { i };
            [t as f64, (t * t) as f64, 1., 1. + (t % 5) as f64][c]
        } else if c == 3 {
            1. + (byte(data, 16 + 4 * (i * nv + j) + c) % 8) as f64
        } else {
            (byte(data, 16 + 4 * (i * nv + j) + c) as i8 as f64)
                * if mode == 1 { scale } else { 1. }
        }
    };
    let poles: Vec<_> = (0..nu)
        .flat_map(|i| {
            (0..nv).map(move |j| Point3::new(scalar(i, j, 0), scalar(i, j, 1), scalar(i, j, 2)))
        })
        .collect();
    let weights: Vec<_> = (0..nu)
        .flat_map(|i| (0..nv).map(move |j| scalar(i, j, 3)))
        .collect();
    let invalid = poles
        .iter()
        .any(|p| p.to_array().iter().any(|x| !x.is_finite()))
        || weights.iter().any(|w| !w.is_finite() || *w <= 0.);
    let [u, v] = axes;
    let original = BSplineSurface3::new(u, v, poles, Some(weights));
    if invalid {
        assert!(matches!(
            original,
            Err(Error::NonFinite(_) | Error::InvalidSurface(_))
        ));
        return;
    }
    let original = original.unwrap().to_exact();
    let fraction = R::new((1 + u16::from(byte(data, 7))).into(), 257.into());
    let cuts: [R; 2] = std::array::from_fn(|i| {
        let a = if i == 0 {
            original.u_knots()
        } else {
            original.v_knots()
        };
        let lo = &a.domain()[0];
        let hi = a.knots().iter().find(|k| *k > lo).unwrap();
        lo + (hi - lo) * &fraction
    });
    let target = 1 + byte(data, 8) as usize % degrees[axis];
    mark("input construction");
    let edited = if op == 3 && mode == 3 {
        reference::check_removal(&original, axis, &original.domain()[axis][0], 0).unwrap()
    } else if op == 4 {
        let a = if axis == 0 {
            original.u_knots()
        } else {
            original.v_knots()
        };
        let k = a
            .knots()
            .iter()
            .find(|k| *k > &a.domain()[0] && *k < &a.domain()[1])
            .unwrap();
        let idx = a.knots().iter().position(|x| x == k).unwrap();
        if a.pole_count() > a.degree() + a.multiplicities()[idx] {
            reference::check_removal(&original, axis, k, 0).unwrap_or_else(|| original.clone())
        } else {
            insert(&original, axis, &cuts[axis], target, &mut mark)
        }
    } else {
        let first = insert(&original, axis, &cuts[axis], target, &mut mark);
        if op == 1 {
            reference::check_removal(&first, axis, &cuts[axis], 0).unwrap()
        } else if op == 2 {
            let other = insert(
                &first,
                1 - axis,
                &cuts[1 - axis],
                degrees[1 - axis],
                &mut mark,
            );
            let ur = [(cuts[0].clone(), if axis == 0 { target } else { degrees[0] })];
            let vr = [(cuts[1].clone(), if axis == 1 { target } else { degrees[1] })];
            assert_eq!(other, original.refined(&ur, &vr).unwrap());
            assert_eq!(
                other.exchanged_uv(),
                original.exchanged_uv().refined(&vr, &ur).unwrap()
            );
            other
        } else {
            first
        }
    };
    mark("remaining edit and removal equations");
    let domain: [[R; 2]; 2] = std::array::from_fn(|i| {
        let a = if i == 0 {
            edited.u_knots()
        } else {
            edited.v_knots()
        };
        let lo = &a.domain()[0];
        let hi = a.knots().iter().find(|k| *k > lo).unwrap();
        [lo.clone(), hi.clone()]
    });
    let [[ua, ub], [va, vb]] = &domain;
    let p = edited.bezier_patches_in(ua, ub, va, vb).unwrap().remove(0);
    mark("kernel patch extraction");
    // A refined cell lies inside an original span. Successful removal proves
    // equality across the old interior knots, so its polynomial also extends
    // across the merged cell, including a shifted periodic origin.
    let expected = tensor::exact_patch(&original, domain, knot_reference::basis_coefficients);
    mark("independent tensor coefficient construction");
    tensor::check(&p, &expected, [&fraction, &fraction]);
    mark("complete patch coefficients and quotient jets");
    let u = &p.domain()[0][0] + (&p.domain()[0][1] - &p.domain()[0][0]) * &fraction;
    let v = &p.domain()[1][0] + (&p.domain()[1][1] - &p.domain()[1][0]) * &fraction;
    let jet = edited
        .exact_evaluate(&u, &v, D::Second, [S::Automatic; 2])
        .unwrap();
    assert_eq!(jet, p.exact_evaluate(&u, &v, D::Second).unwrap());
    mark("kernel exact surface jets");
    let iso = edited
        .u_iso(&u)
        .unwrap()
        .bezier_arcs_in(&p.domain()[1][0], &p.domain()[1][1])
        .unwrap()
        .remove(0);
    crate::bezier_reference::check(&iso, &expected.iso(0, &fraction), &fraction);
    mark("isocurve coefficients and quotient jets");
    for i in 0..2 {
        let a = if i == 0 {
            edited.u_knots()
        } else {
            edited.v_knots()
        };
        if a.is_periodic() {
            let mut uv = [u.clone(), v.clone()];
            uv[i] += (&a.domain()[1] - &a.domain()[0]) * R::from_integer(BigInt::from(1) << 2048);
            assert_eq!(
                jet,
                edited
                    .exact_evaluate(&uv[0], &uv[1], D::Second, [S::Automatic; 2])
                    .unwrap()
            );
        }
    }
    mark("independent patches isocurves jets and enclosures");
    let bad = R::new_raw(1.into(), 0.into());
    assert!(edited.refined(&[], &[(bad.clone(), 0)]).is_err());
    assert!(edited.u_iso(&bad).is_err());
    assert!(edited.v_iso(&bad).is_err());
    mark("rejection checks");
}
