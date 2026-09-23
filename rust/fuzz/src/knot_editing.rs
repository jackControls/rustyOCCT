//! Full-support polynomial identities, independent removal feasibility, exact
//! jets, periodic seam edits and rational operation sequences.
use crate::knot_reference as reference;
use crate::{bezier_reference as bernstein, byte};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use reference::integer as r;
use rusty_occt::curve::{DerivativeOrder as D, KnotSide as S};
use rusty_occt::{BSplineCurve3, Error, ExactBSplineCurve3, Point3};

pub fn check_knot_editing(data: &[u8]) {
    check(data, |_| {});
}
pub fn profile_knot_editing(
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

fn remove(
    curve: &ExactBSplineCurve3,
    u: &R,
    target: usize,
    mark: &mut impl FnMut(&'static str),
) -> ExactBSplineCurve3 {
    let actual = curve.remove_knot(u, target).unwrap();
    mark("kernel removal");
    let expected = reference::removed(curve, u, target);
    mark("oracle removal equations");
    assert_eq!(actual, expected);
    actual.unwrap_or_else(|| curve.clone())
}

fn check(data: &[u8], mut mark: impl FnMut(&'static str)) {
    let mode = byte(data, 0) % 4;
    let degree = 1 + byte(data, 1) as usize
        % match mode {
            0 => 3,
            1 => 8,
            _ => 25,
        };
    let kind = if mode == 3 { 1 } else { byte(data, 2) % 3 };
    let mult = 1 + byte(data, 3) as usize % degree;
    let op = byte(data, 4) % 6;
    let scale = 2_f64.powi((byte(data, 5) as i32 - 128) * 7);
    let parameter_scale = match byte(data, 6) % 8 {
        0 => f64::from_bits(1),
        1 => 2_f64.powi(-500),
        2 => 2_f64.powi(500),
        _ => 1.,
    };
    let (knots, mults) = if mode == 3 {
        (
            (0..degree + 4).map(|i| i as f64).collect(),
            vec![1; degree + 4],
        )
    } else {
        match kind {
            0 => (vec![0., 1., 4.], vec![degree + 1, mult, degree + 1]),
            1 => (vec![0., 1., 3., 4.], vec![mult, degree, 1, mult]),
            _ => (
                vec![-2., -1., 0., 1., 2., 3., 4.],
                vec![degree, degree, 1, degree, 1, degree, degree],
            ),
        }
    };
    let np = mults.iter().sum::<usize>() - if kind == 1 { mults[0] } else { degree + 1 };
    let scalar = |i| {
        if mode == 0 {
            f64::from_bits(u64::from_le_bytes(std::array::from_fn(|j| {
                byte(data, 16 + 8 * i + j)
            })))
        } else if mode == 3 {
            [2., -3., 5., 1.][i % 4]
        } else if i % 4 == 3 {
            1. + f64::from(byte(data, 16 + i) % 8)
        } else {
            (byte(data, 16 + i) as i8 as f64) * if mode == 1 { scale } else { 1. }
        }
    };
    let poles: Vec<_> = (0..np)
        .map(|i| Point3::new(scalar(4 * i), scalar(4 * i + 1), scalar(4 * i + 2)))
        .collect();
    let weights: Vec<_> = (0..np).map(|i| scalar(4 * i + 3)).collect();
    let invalid = poles
        .iter()
        .any(|p| p.to_array().iter().any(|x| !x.is_finite()))
        || weights.iter().any(|w| !w.is_finite() || *w <= 0.);
    let build = if kind == 1 {
        BSplineCurve3::new_periodic
    } else {
        BSplineCurve3::new
    };
    let curve = build(
        degree,
        poles,
        Some(weights),
        knots.into_iter().map(|k| k * parameter_scale).collect(),
        mults,
    );
    if invalid {
        assert!(matches!(
            curve,
            Err(Error::NonFinite(_) | Error::InvalidCurve(_))
        ));
        return;
    }
    let original = curve.unwrap().to_exact();
    let a = &original.domain()[0];
    let b = original.knots().iter().find(|k| *k > a).unwrap();
    let fraction = R::new((1 + u16::from(byte(data, 7))).into(), 257.into());
    let u = a + (b - a) * &fraction;
    let v = (&u + b) / r(2);
    let target = 1 + byte(data, 8) as usize % degree;
    mark("input construction");
    let edited = match op {
        0 => original.insert_knot(&u, target).unwrap(),
        1 => {
            let actual = original
                .refined(&[
                    (v.clone(), degree),
                    (u.clone(), target),
                    (v.clone(), 1),
                    (u.clone(), 0),
                ])
                .unwrap();
            assert_eq!(
                actual,
                original
                    .insert_knot(&u, target)
                    .unwrap()
                    .insert_knot(&v, degree)
                    .unwrap()
            );
            actual
        }
        2 => {
            let fine = original.insert_knot(&u, target).unwrap();
            let coarse = remove(&fine, &u, 0, &mut mark);
            assert_eq!(coarse, original);
            coarse
        }
        3 if kind == 1 => {
            let fine = original
                .refined(&[(original.domain()[1].clone(), degree), (a.clone(), 1)])
                .unwrap();
            assert_eq!(fine, original.insert_knot(a, degree).unwrap());
            let coarse = remove(&fine, a, original.multiplicities()[0], &mut mark);
            assert_eq!(coarse, original);
            fine
        }
        4 | 5 if np > degree + 1 => {
            let index = if kind == 1 && (op == 5 || byte(data, 9) % 2 == 0) {
                0
            } else {
                original
                    .knots()
                    .iter()
                    .position(|k| k > a && k < &original.domain()[1])
                    .unwrap()
            };
            let knot = &original.knots()[index];
            remove(
                &original,
                knot,
                original.multiplicities()[index] - 1,
                &mut mark,
            )
        }
        _ => original.refined(&[(u.clone(), degree), (v, 1)]).unwrap(),
    };
    mark("kernel edits");
    assert!(reference::equal(&original, &edited));
    mark("oracle complete identity");
    assert_eq!(edited.degree(), degree);
    assert_eq!(edited.is_periodic(), kind == 1);
    assert_eq!(
        &edited.domain()[1] - &edited.domain()[0],
        &original.domain()[1] - &original.domain()[0]
    );
    if op == 0 {
        let mut expected: std::collections::BTreeMap<_, _> = original
            .knots()
            .iter()
            .cloned()
            .zip(original.multiplicities().iter().copied())
            .collect();
        expected.insert(u.clone(), target);
        assert_eq!(
            edited
                .knots()
                .iter()
                .cloned()
                .zip(edited.multiplicities().iter().copied())
                .collect::<std::collections::BTreeMap<_, _>>(),
            expected
        );
    }
    // Check an independently selected full polynomial span, its exact
    // Bernstein controls, jets and minimal binary64 enclosures.
    let spans: Vec<_> = edited
        .knots()
        .windows(2)
        .filter(|k| k[0] >= edited.domain()[0] && k[1] <= edited.domain()[1])
        .collect();
    let span = spans[byte(data, 10) as usize % spans.len()];
    let (lo, hi) = (&span[0], &span[1]);
    let arc = edited.bezier_arcs_in(lo, hi).unwrap().remove(0);
    let polynomial = bernstein::Arc {
        degree,
        domain: [lo.clone(), hi.clone()],
        coefficients: reference::polynomial(&edited, lo, hi),
    };
    bernstein::check(&arc, &polynomial, &fraction);
    let probe = lo + (hi - lo) * &fraction;
    let jet = edited
        .exact_evaluate(&probe, D::Second, S::Automatic)
        .unwrap();
    assert_eq!(jet, arc.exact_evaluate(&probe, D::Second).unwrap());
    if kind == 1 {
        let period = &edited.domain()[1] - &edited.domain()[0];
        let translated = probe + period * R::from_integer(BigInt::from(1) << 2048);
        assert_eq!(
            jet,
            edited
                .exact_evaluate(&translated, D::Second, S::Automatic)
                .unwrap()
        );
    }
    mark("exact extraction jets and bounds");
    let bad = R::new_raw(1.into(), 0.into());
    assert!(edited.insert_knot(&bad, 0).is_err());
    assert!(edited
        .refined(&[(edited.domain()[0].clone(), degree + 2)])
        .is_err());
    assert!(edited
        .exact_evaluate(&bad, D::Position, S::Automatic)
        .is_err());
    mark("rejection checks");
}
