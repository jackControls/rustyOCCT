//! Complete independent control reconstruction, including inactive support.
//! Binary64 atoms, rational extremes, all degrees, both tensor directions,
//! periodic origins, operation composition, and atomic resource rejection.
use crate::{byte, knot_reference};
#[path = "../../kernel/tests/support/degree_reference.rs"]
mod reference;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::curve::{DerivativeOrder as D, KnotSide as S};
use rusty_occt::{Error, ExactBSplineCurve3, ExactBSplineSurface3, ExactKnotVector};

fn r(n: i32) -> R {
    R::from_integer(n.into())
}

pub fn check_degree_elevation(data: &[u8]) {
    check(data, |_| {});
}

pub fn profile_degree_elevation(
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

fn target(p: usize, code: u8) -> usize {
    match code % 32 {
        0 => p,
        1..=25 => p + (code as usize % 32 - 1) % (26 - p),
        26 => p - 1,
        27 => 26,
        28 => usize::MAX,
        _ => 25,
    }
}

fn count(axis: &ExactKnotVector, degree: usize) -> Option<usize> {
    if !(axis.degree()..=25).contains(&degree) {
        return None;
    }
    let active = axis
        .knots()
        .iter()
        .filter(|k| *k >= &axis.domain()[0] && *k <= &axis.domain()[1])
        .count();
    Some(axis.pole_count() + (degree - axis.degree()) * (active - 1))
}

fn pow2(exponent: i32) -> R {
    if exponent >= 0 {
        R::from_integer(BigInt::from(1) << exponent as usize)
    } else {
        R::new(1.into(), BigInt::from(1) << (-exponent) as usize)
    }
}

fn axis(data: &[u8], index: usize, p: usize) -> ExactKnotVector {
    let kind = byte(data, 4 + index) % 5;
    let m = 1 + byte(data, 12 + index) as usize % p;
    let (knots, mults) = match kind {
        0 => (vec![0, 1, 4], vec![p + 1, m, p + 1]),
        1 => (
            (0..2 * p + 4).map(|i| i as i32).collect(),
            vec![1; 2 * p + 4],
        ),
        2 => (vec![-2, -1, 0, 1, 2, 3, 4], vec![p, p, 1, p, 1, p, p]),
        3 => (vec![0, 1, 3], vec![m, p, m]),
        _ => ((0..p + 4).map(|i| i as i32).collect(), vec![1; p + 4]),
    };
    let scale = match byte(data, 10 + index) % 8 {
        0 => pow2(-1100),
        1 => pow2(4096),
        2 => R::new(1.into(), 257.into()),
        _ => r(1),
    };
    let knots = knots.into_iter().map(|k| r(k) * &scale).collect();
    if kind >= 3 {
        ExactKnotVector::new_periodic(p, knots, mults)
    } else {
        ExactKnotVector::new(p, knots, mults)
    }
    .unwrap()
}

fn controls(data: &[u8], n: usize) -> Vec<[R; 4]> {
    let mode = byte(data, 1) % 4;
    let scale = match mode {
        1 => pow2((i32::from(byte(data, 9)) - 128) * 7),
        3 => match byte(data, 9) % 3 {
            0 => pow2(-1100),
            1 => pow2(4096),
            _ => R::new(1.into(), 257.into()),
        },
        _ => r(1),
    };
    (0..n)
        .map(|i| {
            std::array::from_fn(|c| {
                if mode == 0 {
                    let x = f64::from_bits(u64::from_le_bytes(std::array::from_fn(|j| {
                        byte(data, 16 + 8 * (4 * i + c) + j)
                    })));
                    // Invalid binary64 atoms become explicitly invalid rational atoms;
                    // construction must reject them before any editing is attempted.
                    R::from_float(x).unwrap_or_else(|| R::new_raw(1.into(), 0.into()))
                } else if c == 3 {
                    r(1 + i32::from(byte(data, 16 + 4 * i + c) % 8))
                } else {
                    r(byte(data, 16 + 4 * i + c) as i8 as i32) * &scale
                }
            })
        })
        .collect()
}

fn valid(controls: &[[R; 4]]) -> bool {
    controls
        .iter()
        .all(|p| p.iter().all(|c| *c.denom() > BigInt::from(0)) && p[3] > r(0))
}

fn interior(axis: &ExactKnotVector, fraction: &R) -> R {
    let a = &axis.domain()[0];
    let b = axis.knots().iter().find(|k| *k > a).unwrap();
    a + (b - a) * fraction
}

fn rejected<T>(result: rusty_occt::Result<T>, expected_count: Option<usize>) {
    let error = result
        .err()
        .expect("invalid degree or output count must fail");
    if expected_count.is_none() {
        assert!(matches!(error, Error::InvalidSpline(_)));
    } else {
        assert!(matches!(error, Error::LimitExceeded(_)));
    }
}

fn check(data: &[u8], mut mark: impl FnMut(&'static str)) {
    let family = byte(data, 0) % 3;
    let degrees: [usize; 2] = std::array::from_fn(|i| 1 + byte(data, 2 + i) as usize % 25);
    let desired = std::array::from_fn(|i| target(degrees[i], byte(data, 6 + i)));
    if family == 2 {
        check_limit(data, degrees, desired);
        mark("atomic count and degree limits");
        return;
    }
    let u = axis(data, 0, degrees[0]);
    let fraction = R::new((1 + u16::from(byte(data, 14))).into(), 257.into());
    let op = byte(data, 8) % 3;
    if family == 0 {
        let grid = controls(data, u.pole_count());
        let okay = valid(&grid);
        let original = ExactBSplineCurve3::from_homogeneous(u, grid);
        if !okay {
            assert!(original.is_err());
            return;
        }
        let original = original.unwrap();
        let expected_count = count(original.knot_vector(), desired[0]);
        mark("input construction");
        let result = original.elevated(desired[0]);
        if expected_count.is_none_or(|n| n > 4096) {
            rejected(result, expected_count);
            return;
        }
        let result = result.unwrap();
        mark("kernel degree elevation");
        assert_eq!(result, reference::elevated(&original, desired[0]));
        mark("complete independent coefficient reconstruction");
        assert_eq!(result.domain(), original.domain());
        assert_eq!(result.homogeneous_poles().len(), expected_count.unwrap());
        assert_eq!(result.is_periodic(), original.is_periodic());
        let u = interior(original.knot_vector(), &fraction);
        if op == 0 {
            let mid = (degrees[0] + desired[0]) / 2;
            assert_eq!(
                original
                    .elevated(mid)
                    .unwrap()
                    .elevated(desired[0])
                    .unwrap(),
                result
            );
        } else if op == 1 {
            let m = 1 + byte(data, 15) as usize % degrees[0];
            assert_eq!(
                original
                    .insert_knot(&u, m)
                    .unwrap()
                    .elevated(desired[0])
                    .unwrap(),
                result.insert_knot(&u, m + desired[0] - degrees[0]).unwrap()
            );
        } else {
            for k in original.domain() {
                for side in [S::Left, S::Automatic, S::Right] {
                    assert_eq!(
                        original.exact_evaluate(k, D::Second, side),
                        result.exact_evaluate(k, D::Second, side)
                    );
                }
            }
        }
        assert_eq!(
            original.exact_evaluate(&u, D::Second, S::Automatic),
            result.exact_evaluate(&u, D::Second, S::Automatic)
        );
    } else {
        let v = axis(data, 1, degrees[1]);
        let n = u.pole_count() * v.pole_count();
        let expected_count = count(&u, desired[0])
            .zip(count(&v, desired[1]))
            .map(|(a, b)| a * b);
        let grid = controls(data, n);
        let okay = valid(&grid) && n <= 4096;
        let original = ExactBSplineSurface3::from_homogeneous(u, v, grid);
        if !okay {
            assert!(original.is_err());
            return;
        }
        let original = original.unwrap();
        mark("input construction");
        let result = original.elevated(desired[0], desired[1]);
        if expected_count.is_none_or(|n| n > 4096) {
            rejected(result, expected_count);
            return;
        }
        let result = result.unwrap();
        mark("kernel degree elevation");
        assert_eq!(result, reference::elevated_surface(&original, desired));
        mark("complete independent coefficient reconstruction");
        assert_eq!(result.domain(), original.domain());
        assert_eq!(result.homogeneous_poles().len(), expected_count.unwrap());
        let u = interior(original.u_knots(), &fraction);
        let v = interior(original.v_knots(), &fraction);
        if op == 0 {
            assert_eq!(
                original
                    .elevated(degrees[0], desired[1])
                    .unwrap()
                    .elevated(desired[0], desired[1])
                    .unwrap(),
                result
            );
        } else if op == 1 {
            assert_eq!(
                original
                    .exchanged_uv()
                    .elevated(desired[1], desired[0])
                    .unwrap(),
                result.exchanged_uv()
            );
        } else {
            let which = byte(data, 15) as usize % 2;
            let m = 1 + byte(data, 15) as usize % degrees[which];
            let before = if which == 0 {
                original.insert_u_knot(&u, m)
            } else {
                original.insert_v_knot(&v, m)
            };
            let after = if which == 0 {
                result.insert_u_knot(&u, m + desired[0] - degrees[0])
            } else {
                result.insert_v_knot(&v, m + desired[1] - degrees[1])
            };
            // Composition itself can cross the exposed grid-size limit. Its
            // success/failure must agree in both operation orders as well.
            match (
                before.and_then(|s| s.elevated(desired[0], desired[1])),
                after,
            ) {
                (Ok(a), Ok(b)) => assert_eq!(a, b),
                (Err(Error::LimitExceeded(_)), Err(Error::LimitExceeded(_))) => (),
                outcomes => panic!("refinement/elevation composition disagrees: {outcomes:?}"),
            }
        }
        assert_eq!(
            original.exact_evaluate(&u, &v, D::Second, [S::Automatic; 2]),
            result.exact_evaluate(&u, &v, D::Second, [S::Automatic; 2])
        );
    }
    mark("composition and exact jets");
}

fn uniform_axis(p: usize, n: usize) -> ExactKnotVector {
    let knots = (0..=n - p).map(|k| r(k as i32)).collect::<Vec<_>>();
    let mut mults = vec![1; knots.len()];
    mults[0] = p + 1;
    *mults.last_mut().unwrap() = p + 1;
    ExactKnotVector::new(p, knots, mults).unwrap()
}

fn check_limit(data: &[u8], degrees: [usize; 2], desired: [usize; 2]) {
    let control = [r(1), r(-2), r(3), r(1)];
    if byte(data, 8) % 2 == 0 {
        let axis = uniform_axis(degrees[0], 4096);
        let expected = count(&axis, desired[0]);
        let original = ExactBSplineCurve3::from_homogeneous(axis, vec![control; 4096]).unwrap();
        let result = original.elevated(desired[0]);
        if expected.is_some_and(|n| n <= 4096) {
            assert_eq!(result.unwrap(), original);
        } else {
            rejected(result, expected);
        }
    } else {
        let u = uniform_axis(degrees[0], 64);
        let v = uniform_axis(degrees[1], 64);
        let expected = count(&u, desired[0])
            .zip(count(&v, desired[1]))
            .map(|(a, b)| a * b);
        let original = ExactBSplineSurface3::from_homogeneous(u, v, vec![control; 4096]).unwrap();
        let result = original.elevated(desired[0], desired[1]);
        if expected.is_some_and(|n| n <= 4096) {
            assert_eq!(result.unwrap(), original);
        } else {
            rejected(result, expected);
        }
    }
}
