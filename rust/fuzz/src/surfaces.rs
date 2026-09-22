//! Independent tensor basis derivatives and closed bivariate quotient rules.
use super::{
    byte, one,
    splines::{axis, basis_jet, bounds, parameters, rat},
    zero,
};
use num_rational::BigRational as R;
use rusty_occt::curve::{DerivativeOrder as D, KnotSide as S};
use rusty_occt::{BSplineSurface3, Error, KnotVector, Point3};

const PARTIALS: [(usize, usize); 6] = [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1)];
fn reference(
    poles: &[Point3],
    weights: &[f64],
    nv: usize,
    ub: &[[R; 3]],
    vb: &[[R; 3]],
) -> [[R; 3]; 6] {
    let nu = poles.len() / nv;
    let h: [[R; 4]; 6] = std::array::from_fn(|r| {
        let (a, b) = PARTIALS[r];
        std::array::from_fn(|c| {
            let mut value = zero();
            for (i, u) in ub.iter().enumerate() {
                for (j, v) in vb.iter().enumerate() {
                    let k = (i % nu) * nv + j % nv;
                    value += &u[a]
                        * &v[b]
                        * rat(weights[k])
                        * if c == 3 {
                            one()
                        } else {
                            rat(poles[k].to_array()[c])
                        };
                }
            }
            value
        })
    });
    let w = &h[0][3];
    let two = rat(2.);
    std::array::from_fn(|r| {
        std::array::from_fn(|c| match r {
            0 => &h[0][c] / w,
            1 | 2 => (&h[r][c] * w - &h[0][c] * &h[r][3]) / (w * w),
            3 | 4 => {
                let first = r - 2;
                (&h[r][c] * w * w
                    - &h[0][c] * &h[r][3] * w
                    - &two * &h[first][c] * w * &h[first][3]
                    + &two * &h[0][c] * &h[first][3] * &h[first][3])
                    / (w * w * w)
            }
            _ => {
                (&h[5][c] * w * w
                    - &h[1][c] * w * &h[2][3]
                    - &h[2][c] * w * &h[1][3]
                    - &h[0][c] * w * &h[5][3]
                    + &two * &h[0][c] * &h[1][3] * &h[2][3])
                    / (w * w * w)
            }
        })
    })
}

pub fn check_surfaces(data: &[u8]) {
    let mode = byte(data, 0) % 4;
    let limit = if mode < 2 {
        2
    } else if mode == 2 {
        5
    } else {
        25
    };
    let mut degrees = [
        1 + byte(data, 1) as usize % limit,
        1 + byte(data, 2) as usize % if mode == 3 { 2 } else { limit },
    ];
    if mode == 3 && byte(data, 3) & 4 != 0 {
        degrees.swap(0, 1);
    }
    let periodic = [byte(data, 3) & 1 != 0, byte(data, 3) & 2 != 0];
    let sides = [byte(data, 4), byte(data, 5)].map(|x| match x % 3 {
        0 => S::Automatic,
        1 => S::Left,
        _ => S::Right,
    });
    let order = byte(data, 6) % 3;
    let count = [1, 3, 6][order as usize];
    let scale = 2.0_f64.powi(byte(data, 7) as i32 - 128);
    let scalar = |i: usize| {
        if mode < 2 {
            f64::from_bits(u64::from_le_bytes(std::array::from_fn(|j| {
                byte(data, 12 + 8 * i + j)
            })))
        } else {
            (byte(data, 12 + i) as i8 as f64) * if mode == 2 { scale } else { 1. }
        }
    };
    let mults = std::array::from_fn::<_, 2, _>(|i| {
        if periodic[i] {
            vec![1, degrees[i], 1]
        } else {
            vec![degrees[i] + 1, 1, degrees[i] + 1]
        }
    });
    let sizes = std::array::from_fn::<_, 2, _>(|i| {
        if periodic[i] {
            degrees[i] + 1
        } else {
            degrees[i] + 2
        }
    });
    let knots = std::array::from_fn::<_, 2, _>(|i| {
        if mode == 0 {
            (0..3)
                .map(|j| scalar(4 * sizes[0] * sizes[1] + 4 * i + j))
                .collect::<Vec<_>>()
        } else {
            vec![
                0.,
                if mode == 2 { scale } else { 1. },
                if mode == 2 { 3. * scale } else { 3. },
            ]
        }
    });
    let poles: Vec<_> = (0..sizes[0] * sizes[1])
        .map(|i| Point3::new(scalar(4 * i), scalar(4 * i + 1), scalar(4 * i + 2)))
        .collect();
    let weights: Vec<_> = (0..poles.len())
        .map(|i| {
            let x = scalar(4 * i + 3);
            match mode {
                0 => x,
                1 => {
                    if x == 0. {
                        1.
                    } else {
                        x.abs()
                    }
                }
                _ => 1. + x.abs(),
            }
        })
        .collect();
    let valid = poles.iter().flat_map(|p| p.to_array()).all(f64::is_finite)
        && weights.iter().all(|w| w.is_finite() && *w > 0.)
        && knots
            .iter()
            .all(|k| k.iter().all(|x| x.is_finite()) && k.windows(2).all(|w| w[0] < w[1]));
    let axes = std::array::from_fn::<_, 2, _>(|i| {
        if periodic[i] {
            KnotVector::new_periodic(degrees[i], knots[i].clone(), mults[i].clone())
        } else {
            KnotVector::new(degrees[i], knots[i].clone(), mults[i].clone())
        }
    });
    let surface = match axes {
        [Ok(u), Ok(v)] => BSplineSurface3::new(u, v, poles.clone(), Some(weights.clone())),
        _ => {
            assert!(!valid);
            return;
        }
    };
    if !valid {
        assert!(surface.is_err());
        return;
    }
    let surface = surface.unwrap();
    let query = std::array::from_fn::<_, 2, _>(|i| {
        if mode == 0 {
            scalar(4 * poles.len() + 4 * i + 3)
        } else {
            let q = byte(data, 8 + i);
            if periodic[i] && q % 9 == 8 {
                if i == 0 {
                    f64::MAX
                } else {
                    -f64::MAX
                }
            } else {
                knots[i][1] * [0., 0.25, 1., 1.5, 3., -3., 6., -0.5][q as usize % 8]
            }
        }
    });
    let requested = match order {
        0 => D::Position,
        1 => D::First,
        _ => D::Second,
    };
    let actual = surface.evaluate(query[0], query[1], requested, sides);
    if !query.iter().all(|q| q.is_finite()) {
        assert!(actual.is_err());
        return;
    }
    for i in 0..2 {
        if !periodic[i]
            && (query[i] < knots[i][0]
                || query[i] > knots[i][2]
                || (query[i] == knots[i][0] && sides[i] == S::Left)
                || (query[i] == knots[i][2] && sides[i] == S::Right))
        {
            assert!(matches!(actual, Err(Error::OutOfDomain(_))));
            return;
        }
    }
    let bases = std::array::from_fn::<_, 2, _>(|i| {
        let (flat, start, end) = axis(degrees[i], &knots[i], &mults[i], periodic[i]);
        parameters(&flat, &start, &end, query[i], sides[i], periodic[i])
            .iter()
            .map(|(u, span)| basis_jet(degrees[i], &flat, u, *span))
            .collect::<Vec<_>>()
    });
    let expected = reference(&poles, &weights, sizes[1], &bases[0][0], &bases[1][0]);
    for (i, u) in bases[0].iter().enumerate() {
        for (j, v) in bases[1].iter().enumerate() {
            if (i != 0 || j != 0)
                && expected[..count] != reference(&poles, &weights, sizes[1], u, v)[..count]
            {
                assert!(matches!(actual, Err(Error::DiscontinuousDerivative)));
                return;
            }
        }
    }
    let max = rat(f64::MAX);
    if expected[..count]
        .iter()
        .flatten()
        .any(|x| x > &max || x < &-max.clone())
    {
        assert!(matches!(actual, Err(Error::Unrepresentable(_))));
        return;
    }
    let actual = actual.unwrap();
    for (i, (a, b)) in PARTIALS.into_iter().enumerate() {
        if i == 0 {
            for (b, r) in actual.position_bounds().into_iter().zip(&expected[0]) {
                bounds(b, r);
            }
        } else if i < count {
            for (b, r) in actual
                .derivative_bounds(a, b)
                .unwrap()
                .into_iter()
                .zip(&expected[i])
            {
                bounds(b, r);
            }
        } else {
            assert!(actual.derivative_bounds(a, b).is_none());
        }
    }
}
