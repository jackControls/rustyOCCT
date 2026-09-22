//! Complete known-factor Bezier contacts and independent piecewise rational
//! line/plane intersections, including periodic seams and contained runs.
use super::{byte, one, zero};
use crate::splines::{bounds, rat};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::intersection::{spline_plane, Plane3, SplinePlaneContact as C};
use rusty_occt::{BSplineCurve3, Point3};
use std::cmp::Ordering;
use std::collections::BTreeMap;

fn integer(n: i64) -> R {
    R::from_integer(BigInt::from(n))
}
fn choose(n: usize, k: usize) -> u64 {
    (0..k).fold(1, |a, i| a * (n - i) as u64 / (i + 1) as u64)
}
fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
fn frame(data: &[u8]) -> ([[f64; 3]; 3], f64) {
    (
        if byte(data, 3) & 1 == 0 {
            [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        } else {
            [[2., 1., 0.], [0., 2., 1.], [1., -2., 4.]]
        },
        2f64.powi((i32::from(byte(data, 4)) - 128) * 4),
    )
}
fn transform(p: [f64; 3], frame: [[f64; 3]; 3], scale: f64) -> Point3 {
    let [x, y, z] =
        std::array::from_fn(|c| (0..3).map(|i| p[i] * frame[i][c]).sum::<f64>() * scale);
    Point3::new(x, y, z)
}
fn plane(frame: [[f64; 3]; 3], scale: f64) -> Plane3 {
    Plane3::through_points(
        Point3::new(0., 0., 0.),
        transform([1., 0., 0.], frame, scale),
        transform([0., 1., 0.], frame, scale),
    )
    .unwrap()
}

fn bezier(data: &[u8]) {
    let degree = 1 + usize::from(byte(data, 1) % 8);
    let mut polynomial = vec![one()];
    let mut expected = BTreeMap::<R, usize>::new();
    for i in 0..degree {
        let root = R::new(
            BigInt::from(i32::from(byte(data, 10 + i) % 7) - 1),
            BigInt::from(4),
        );
        *expected.entry(root.clone()).or_default() += 1;
        let mut next = vec![zero(); polynomial.len() + 1];
        for (j, a) in polynomial.iter().enumerate() {
            next[j] -= a * &root;
            next[j + 1] += a;
        }
        polynomial = next;
    }
    expected.retain(|r, _| r >= &zero() && r <= &one());
    let denominator = (0..=degree)
        .map(|k| choose(degree, k))
        .fold(1, |a, b| a / gcd(a, b) * b)
        * 4u64.pow(degree as u32);
    let (frame, scale) = frame(data);
    let mut poles = Vec::new();
    let mut weights = Vec::new();
    for i in 0..=degree {
        let b: R = (0..=i)
            .map(|k| {
                &polynomial[k] * R::new(BigInt::from(choose(i, k)), BigInt::from(choose(degree, k)))
            })
            .sum::<R>()
            * R::from_integer(BigInt::from(denominator));
        assert_eq!(b.denom(), &BigInt::from(1));
        let z = b.numer().to_string().parse::<f64>().unwrap();
        assert_eq!(rat(z), b);
        let w = 2f64.powi(i32::from(byte(data, 30 + i) % 9) - 4);
        let p = [
            f64::from(byte(data, 50 + 2 * i) as i8),
            f64::from(byte(data, 51 + 2 * i) as i8),
            z / w,
        ];
        poles.push(transform(p, frame, scale));
        weights.push(w);
    }
    let parameter_scale = 2f64.powi((i32::from(byte(data, 2)) - 128) * 4);
    let curve = BSplineCurve3::new(
        degree,
        poles.clone(),
        Some(weights.clone()),
        vec![-parameter_scale, 3. * parameter_scale],
        vec![degree + 1; 2],
    )
    .unwrap();
    let result = spline_plane(&curve, &plane(frame, scale)).unwrap();
    assert!(result.overlaps().is_empty());
    assert_eq!(result.points().len(), expected.len());
    for (point, (root, order)) in result.points().iter().zip(expected) {
        let u = (-one() + integer(4) * &root) * rat(parameter_scale);
        bounds(point.parameter(), &u);
        let basis: Vec<R> = (0..=degree)
            .map(|i| {
                R::from_integer(BigInt::from(choose(degree, i)))
                    * root.pow(i as i32)
                    * (one() - &root).pow((degree - i) as i32)
                    * rat(weights[i])
            })
            .collect();
        let w: R = basis.iter().sum();
        for (c, b) in point.coordinate_bounds().into_iter().enumerate() {
            let expected: R = poles
                .iter()
                .zip(&basis)
                .map(|(p, b)| rat(p.to_array()[c]) * b)
                .sum::<R>()
                / &w;
            bounds(b, &expected);
        }
        assert_eq!(
            point.multiplicities(),
            [
                (root > zero()).then_some(order),
                (root < one()).then_some(order)
            ]
        );
        assert_eq!(
            point.contact(),
            if root == zero() || root == one() {
                C::Boundary
            } else if order % 2 == 0 {
                C::Tangent
            } else {
                C::Crossing
            }
        );
    }
}

struct Hit {
    position: [R; 3],
    signs: [Option<Ordering>; 2],
}
fn polygon(data: &[u8]) {
    let periodic = byte(data, 0) & 2 != 0;
    let spans = 2 + usize::from(byte(data, 1) % 7);
    let np = spans + usize::from(!periodic);
    let (frame, scale) = frame(data);
    let mut poles = Vec::new();
    let mut weights = Vec::new();
    let mut heights = Vec::new();
    for i in 0..np {
        let z = i32::from(byte(data, 10 + i) % 5) - 2;
        poles.push(transform(
            [i as f64, f64::from(byte(data, 30 + i) as i8), f64::from(z)],
            frame,
            scale,
        ));
        weights.push(2f64.powi(i32::from(byte(data, 50 + i) % 9) - 4));
        heights.push(integer(i64::from(z)) * rat(weights[i]));
    }
    let mut knots = vec![0.];
    for i in 0..spans {
        knots.push(knots[i] + 1. + f64::from(byte(data, 70 + i) % 4));
    }
    let mut mults = vec![1; spans + 1];
    if !periodic {
        mults[0] = 2;
        mults[spans] = 2;
    }
    let curve = if periodic {
        BSplineCurve3::new_periodic(
            1,
            poles.clone(),
            Some(weights.clone()),
            knots.clone(),
            mults,
        )
    } else {
        BSplineCurve3::new(
            1,
            poles.clone(),
            Some(weights.clone()),
            knots.clone(),
            mults,
        )
    }
    .unwrap();
    let mut hits = BTreeMap::<R, Hit>::new();
    let mut overlaps: Vec<(f64, f64)> = Vec::new();
    for i in 0..spans {
        let j = (i + 1) % np;
        let (a, b) = (&heights[i], &heights[j]);
        if a == &zero() && b == &zero() {
            if let Some(last) = overlaps.last_mut().filter(|last| last.1 == knots[i]) {
                last.1 = knots[i + 1];
            } else {
                overlaps.push((knots[i], knots[i + 1]));
            }
            continue;
        }
        if a == b {
            continue;
        }
        let root = -a / (b - a);
        if root < zero() || root > one() {
            continue;
        }
        let u = rat(knots[i]) + (rat(knots[i + 1]) - rat(knots[i])) * &root;
        let w = (one() - &root) * rat(weights[i]) + &root * rat(weights[j]);
        let position = std::array::from_fn(|c| {
            ((one() - &root) * rat(weights[i]) * rat(poles[i].to_array()[c])
                + &root * rat(weights[j]) * rat(poles[j].to_array()[c]))
                / &w
        });
        let signs = [
            (root > zero()).then(|| a.cmp(&zero())),
            (root < one()).then(|| b.cmp(&zero())),
        ];
        if let Some(previous) = hits.get_mut(&u) {
            assert_eq!(previous.position, position);
            previous.signs[1] = signs[1];
        } else {
            hits.insert(u, Hit { position, signs });
        }
    }
    hits.retain(|u, _| !overlaps.iter().any(|&(a, b)| &rat(a) <= u && u <= &rat(b)));
    let result = spline_plane(&curve, &plane(frame, scale)).unwrap();
    assert_eq!(
        result
            .overlaps()
            .iter()
            .map(|o| o.parameters())
            .collect::<Vec<_>>(),
        overlaps
    );
    assert_eq!(result.points().len(), hits.len());
    for (point, (u, hit)) in result.points().iter().zip(hits) {
        bounds(point.parameter(), &u);
        for (b, x) in point.coordinate_bounds().into_iter().zip(hit.position) {
            bounds(b, &x);
        }
        assert_eq!(point.multiplicities(), hit.signs.map(|s| s.map(|_| 1)));
        assert_eq!(
            point.contact(),
            match hit.signs {
                [Some(a), Some(b)] =>
                    if a == b {
                        C::Tangent
                    } else {
                        C::Crossing
                    },
                _ => C::Boundary,
            }
        );
    }
}
pub fn check_spline_intersections(data: &[u8]) {
    if byte(data, 0) & 1 == 0 {
        bezier(data)
    } else {
        polygon(data)
    }
}
