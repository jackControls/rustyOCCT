//! Independent Cox-de Boor basis derivatives plus closed quotient formulas.
//! Production differentiates homogeneous pole interpolation instead.
use super::{byte, one, zero};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::curve::{BSplineCurve3, DerivativeOrder as D, KnotSide as S};
use rusty_occt::{Error, Point3, ScalarInterval};

pub(crate) fn rat(x: f64) -> R {
    R::from_float(x).unwrap()
}
pub(crate) fn basis_jet(degree: usize, knots: &[R], u: &R, span: usize) -> Vec<[R; 3]> {
    let mut basis: Vec<[R; 3]> = (0..knots.len() - 1)
        .map(|i| [if i == span { one() } else { zero() }, zero(), zero()])
        .collect();
    for p in 1..=degree {
        let previous = basis;
        basis = (0..previous.len() - 1)
            .map(|i| {
                let a = &knots[i + p] - &knots[i];
                let b = &knots[i + p + 1] - &knots[i + 1];
                std::array::from_fn(|r| {
                    if r == 0 {
                        (if a == zero() {
                            zero()
                        } else {
                            (u - &knots[i]) * &previous[i][0] / &a
                        }) + (if b == zero() {
                            zero()
                        } else {
                            (&knots[i + p + 1] - u) * &previous[i + 1][0] / &b
                        })
                    } else {
                        R::from_integer(BigInt::from(p))
                            * ((if a == zero() {
                                zero()
                            } else {
                                &previous[i][r - 1] / &a
                            }) - (if b == zero() {
                                zero()
                            } else {
                                &previous[i + 1][r - 1] / &b
                            }))
                    }
                })
            })
            .collect();
    }
    basis
}
fn reference(
    degree: usize,
    poles: &[Point3],
    weights: &[f64],
    flat: &[R],
    u: &R,
    span: usize,
) -> [[R; 3]; 3] {
    let basis = basis_jet(degree, flat, u, span);
    let h: [[R; 4]; 3] = std::array::from_fn(|r| {
        std::array::from_fn(|c| {
            basis
                .iter()
                .enumerate()
                .map(|(i, basis)| {
                    let p = poles[i % poles.len()];
                    let w = weights[i % weights.len()];
                    &basis[r] * rat(w) * if c == 3 { one() } else { rat(p.to_array()[c]) }
                })
                .sum()
        })
    });
    let w = &h[0][3];
    let two = R::from_integer(BigInt::from(2));
    std::array::from_fn(|r| {
        std::array::from_fn(|c| match r {
            0 => &h[0][c] / w,
            1 => (&h[1][c] * w - &h[0][c] * &h[1][3]) / (w * w),
            _ => {
                (&h[2][c] * w * w - &h[0][c] * &h[2][3] * w - &two * &h[1][c] * w * &h[1][3]
                    + &two * &h[0][c] * &h[1][3] * &h[1][3])
                    / (w * w * w)
            }
        })
    })
}
pub(crate) fn bounds(b: ScalarInterval, expected: &R) {
    let (lo, hi) = (b.lower(), b.upper());
    assert!(lo.is_finite() && hi.is_finite());
    if lo == hi {
        assert_eq!(&rat(lo), expected);
    } else {
        assert_eq!(lo.next_up(), hi);
        assert!(rat(lo) < *expected && *expected < rat(hi));
    }
    assert!(lo <= b.representative() && b.representative() <= hi);
}

// Independently index an infinite periodic knot sequence by Euclidean division.
pub(crate) fn axis(
    degree: usize,
    knots: &[f64],
    mults: &[usize],
    periodic: bool,
) -> (Vec<R>, R, R) {
    let base: Vec<_> = knots
        .iter()
        .zip(mults)
        .flat_map(|(&k, &m)| std::iter::repeat_n(rat(k), m))
        .collect();
    if !periodic {
        return (
            base.clone(),
            base[degree].clone(),
            base[base.len() - degree - 1].clone(),
        );
    }
    let (start, end) = (rat(knots[0]), rat(knots[knots.len() - 1]));
    let period = &end - &start;
    let n = (base.len() - mults[0]) as isize;
    let pad = (degree + 1 - mults[0]) as isize;
    let flat = (-pad..n + mults[0] as isize + pad)
        .map(|i| {
            &base[i.rem_euclid(n) as usize]
                + R::from_integer(BigInt::from(i.div_euclid(n))) * &period
        })
        .collect();
    (flat, start, end)
}
pub(crate) fn parameters(
    flat: &[R],
    start: &R,
    end: &R,
    u: f64,
    side: S,
    periodic: bool,
) -> Vec<(R, usize)> {
    let mut u = rat(u);
    if periodic {
        let period = end - start;
        let mut remainder = (&u - start) % &period;
        if remainder < zero() {
            remainder += &period;
        }
        u = start + remainder;
    }
    let at = |value: R, left: bool| {
        let span = (0..flat.len())
            .rfind(|&i| {
                if left {
                    flat[i] < value
                } else {
                    flat[i] <= value
                }
            })
            .unwrap();
        (value, span)
    };
    let seam = periodic && &u == start;
    let mut choices = vec![if seam && side == S::Left {
        at(end.clone(), true)
    } else {
        at(u.clone(), side == S::Left || &u == end)
    }];
    if side == S::Automatic && (seam || (&u > start && &u < end && flat.contains(&u))) {
        choices.push(at(if seam { end.clone() } else { u }, true));
    }
    choices
}

pub fn check_splines(data: &[u8]) {
    let mode = byte(data, 0) % 4;
    let periodic = byte(data, 7) & 1 != 0;
    let degree = 1 + byte(data, 1) as usize
        % if mode < 2 {
            3
        } else if mode == 2 {
            8
        } else {
            25
        };
    let mult = 1 + byte(data, 2) as usize % degree;
    let count = degree + 1 + mult;
    let side = match byte(data, 3) % 3 {
        0 => S::Automatic,
        1 => S::Left,
        _ => S::Right,
    };
    let order = byte(data, 4) % 3;
    let scale = 2.0_f64.powi(byte(data, 5) as i32 - 128);
    let scalar = |i: usize| {
        if mode < 2 {
            f64::from_bits(u64::from_le_bytes(std::array::from_fn(|j| {
                byte(data, 8 + 8 * i + j)
            })))
        } else {
            (byte(data, 8 + i) as i8 as f64) * if mode == 2 { scale } else { 1. }
        }
    };
    let poles: Vec<_> = (0..count)
        .map(|i| Point3::new(scalar(4 * i), scalar(4 * i + 1), scalar(4 * i + 2)))
        .collect();
    let weights: Vec<_> = (0..count)
        .map(|i| {
            let x = scalar(4 * i + 3);
            match mode {
                0 => x,
                1 => x.abs(),
                _ => 1. + x.abs(),
            }
        })
        .collect();
    let knot_count = if periodic { 4 } else { 3 };
    let knots = if mode == 0 {
        (0..knot_count)
            .map(|i| scalar(4 * count + i))
            .collect::<Vec<_>>()
    } else {
        let locations = if periodic {
            vec![0., 1., 2., 3.]
        } else {
            vec![0., 1., 3.]
        };
        locations
            .into_iter()
            .map(|k| k * if mode == 2 { scale } else { 1. })
            .collect()
    };
    let u = if mode == 0 {
        scalar(4 * count + knot_count)
    } else if periodic && byte(data, 6) % 8 == 7 {
        f64::MAX
    } else {
        knots[1]
            * [0., 0.25, 1., 1.5, 3., -6., 9.]
                [byte(data, 6) as usize % if periodic { 7 } else { 5 }]
    };
    let mults = if periodic {
        vec![mult, degree, 1, mult]
    } else {
        vec![degree + 1, mult, degree + 1]
    };
    let curve = if periodic {
        BSplineCurve3::new_periodic(
            degree,
            poles.clone(),
            Some(weights.clone()),
            knots.clone(),
            mults.clone(),
        )
    } else {
        BSplineCurve3::new(
            degree,
            poles.clone(),
            Some(weights.clone()),
            knots.clone(),
            mults.clone(),
        )
    };
    let valid = poles.iter().flat_map(|p| p.to_array()).all(f64::is_finite)
        && weights.iter().all(|w| w.is_finite() && *w > 0.)
        && knots.iter().all(|k| k.is_finite())
        && knots.windows(2).all(|w| w[0] < w[1]);
    if !valid {
        assert!(curve.is_err());
        return;
    }
    let curve = curve.unwrap();
    let requested = match order {
        0 => D::Position,
        1 => D::First,
        _ => D::Second,
    };
    let actual = curve.evaluate(u, requested, side);
    if !u.is_finite() {
        assert!(matches!(actual, Err(Error::NonFinite(_))));
        return;
    }
    if !periodic
        && (u < knots[0]
            || u > knots[2]
            || (u == knots[0] && side == S::Left)
            || (u == knots[2] && side == S::Right))
    {
        assert!(matches!(actual, Err(Error::OutOfDomain(_))));
        return;
    }
    let (flat, start, end) = axis(degree, &knots, &mults, periodic);
    let choices = parameters(&flat, &start, &end, u, side, periodic);
    let expected = reference(degree, &poles, &weights, &flat, &choices[0].0, choices[0].1);
    for (value, span) in &choices[1..] {
        let opposite = reference(degree, &poles, &weights, &flat, value, *span);
        if expected[..=order as usize] != opposite[..=order as usize] {
            assert!(matches!(actual, Err(Error::DiscontinuousDerivative)));
            return;
        }
    }
    let max = rat(f64::MAX);
    if expected[..=order as usize]
        .iter()
        .flatten()
        .any(|x| x > &max || x < &-max.clone())
    {
        assert!(matches!(actual, Err(Error::Unrepresentable(_))));
        return;
    }
    let actual = actual.unwrap();
    for (i, row) in expected.iter().enumerate().take(order as usize + 1) {
        let components = if i == 0 {
            actual.position_bounds()
        } else {
            actual.derivative_bounds(i).unwrap()
        };
        for (b, r) in components.into_iter().zip(row) {
            bounds(b, r);
        }
    }
    for i in (order as usize + 1).max(1)..=2 {
        assert!(actual.derivative_bounds(i).is_none());
    }
}
