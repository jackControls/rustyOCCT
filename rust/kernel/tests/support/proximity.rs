//! Test-only global convex optimality checks. No active-face search or linear
//! system solve: certify the returned witnesses against the original sets.
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::proximity::{ClosestPair, LinearPrimitive3 as P};
use rusty_occt::{Error, Point3, Result, ScalarInterval};
use std::cmp::Ordering;

fn zero() -> R {
    R::from_integer(BigInt::from(0))
}
fn one() -> R {
    R::from_integer(BigInt::from(1))
}
fn rational(x: f64) -> R {
    R::from_float(x).unwrap()
}
fn sub(a: &[R; 3], b: &[R; 3]) -> [R; 3] {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn dot(a: &[R; 3], b: &[R; 3]) -> R {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

pub fn points(p: &P) -> Vec<Point3> {
    match p {
        P::Point(p) => vec![*p],
        P::Line(p) | P::Segment(p) => p.to_vec(),
        P::Plane(p) => p.defining_points().to_vec(),
        P::Triangle(p) => p.vertices().to_vec(),
    }
}

fn adjacent(lo: f64, hi: f64) {
    assert!(lo <= hi && lo.is_finite() && hi.is_finite());
    if lo != hi {
        assert_eq!(lo.to_bits().abs_diff(hi.to_bits()), 1);
    }
}

fn bound(value: &R, b: Result<ScalarInterval>, squared: bool) {
    let maximum = rational(f64::MAX);
    if *value
        > if squared {
            &maximum * &maximum
        } else {
            maximum.clone()
        }
        || (!squared && *value < -maximum)
    {
        assert!(matches!(b, Err(Error::Unrepresentable(_))));
        return;
    }
    let b = b.unwrap();
    adjacent(b.lower(), b.upper());
    let [lo, hi] = [b.lower(), b.upper()].map(|x| {
        let r = rational(x);
        if squared {
            &r * &r
        } else {
            r
        }
    });
    assert!(lo <= *value && *value <= hi);
    if lo != hi {
        assert!(lo < *value && *value < hi);
    }
}

pub fn certify(a: &P, b: &P, pair: &ClosestPair) {
    let delta = sub(&pair.exact_points()[0], &pair.exact_points()[1]);
    assert_eq!(pair.exact_squared_distance(), &dot(&delta, &delta));
    bound(
        pair.exact_squared_distance(),
        pair.squared_distance_bounds(),
        false,
    );
    bound(pair.exact_squared_distance(), pair.distance_bounds(), true);
    for (operand, shape) in [a, b].into_iter().enumerate() {
        let vertices: Vec<_> = points(shape)
            .into_iter()
            .map(|p| p.to_array().map(rational))
            .collect();
        let p = &pair.exact_points()[operand];
        let parameters = &pair.exact_parameters()[operand];
        assert_eq!(parameters.len(), vertices.len() - 1);
        for k in 0..3 {
            let reconstructed = &vertices[0][k]
                + parameters
                    .iter()
                    .zip(&vertices[1..])
                    .map(|(t, v)| t * (&v[k] - &vertices[0][k]))
                    .sum::<R>();
            assert_eq!(p[k], reconstructed);
            for x in [-f64::MAX, -0., f64::from_bits(1), f64::MAX] {
                assert_eq!(
                    pair.compare_coordinate(operand, k, x).unwrap(),
                    p[k].cmp(&rational(x))
                );
            }
        }
        if matches!(shape, P::Point(_) | P::Segment(_) | P::Triangle(_)) {
            assert!(
                parameters.iter().all(|t| *t >= zero()) && parameters.iter().sum::<R>() <= one()
            );
            for v in &vertices {
                let support = dot(&delta, &sub(v, p));
                assert!(if operand == 0 {
                    support >= zero()
                } else {
                    support <= zero()
                });
            }
        } else {
            for v in &vertices[1..] {
                assert_eq!(dot(&delta, &sub(v, &vertices[0])), zero());
            }
        }
        let maximum = rational(f64::MAX);
        if p.iter().any(|x| *x > maximum || *x < -&maximum) {
            assert!(matches!(
                pair.point_bounds(operand),
                Err(Error::Unrepresentable(_))
            ));
            assert!(matches!(
                pair.point(operand),
                Err(Error::Unrepresentable(_))
            ));
        } else {
            let b = pair.point_bounds(operand).unwrap();
            let representative = pair.point(operand).unwrap().to_array();
            for k in 0..3 {
                let (lo, hi) = (b.min.to_array()[k], b.max.to_array()[k]);
                adjacent(lo, hi);
                let (l, h) = (rational(lo), rational(hi));
                assert!(l <= p[k] && p[k] <= h);
                if lo != hi {
                    assert!(l < p[k] && p[k] < h);
                }
                assert!(lo <= representative[k] && representative[k] <= hi);
            }
        }
        if parameters.iter().any(|x| *x > maximum || *x < -&maximum) {
            assert!(matches!(
                pair.parameter_bounds(operand),
                Err(Error::Unrepresentable(_))
            ));
        } else {
            let bounds = pair.parameter_bounds(operand).unwrap();
            assert_eq!(bounds.len(), parameters.len());
            for (i, (t, b)) in parameters.iter().zip(bounds).enumerate() {
                bound(t, Ok(b), false);
                for x in [b.lower(), b.upper()] {
                    assert_eq!(
                        pair.compare_parameter(operand, i, x).unwrap(),
                        t.cmp(&rational(x))
                    );
                }
            }
        }
    }
    for x in [-f64::MAX, -0., f64::from_bits(1), 1., f64::MAX] {
        let q = rational(x);
        let expected = if x < 0. {
            Ordering::Greater
        } else {
            pair.exact_squared_distance().cmp(&(&q * &q))
        };
        assert_eq!(pair.compare_distance(x).unwrap(), expected);
        assert_eq!(
            pair.compare_squared_distance(x).unwrap(),
            pair.exact_squared_distance().cmp(&q)
        );
    }
}
