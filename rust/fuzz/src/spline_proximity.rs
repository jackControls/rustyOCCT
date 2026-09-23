//! Complete minimum sets from known factors, closed-form parabolas, circles,
//! and exhaustive rational segment projection. No production root isolator or
//! spline evaluation supplies the expected answer.
use crate::{byte, roots::interval, splines::bounds};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::proximity::{closest_points_on_exact_spline_in, SplineClosestSet};
use rusty_occt::{ExactBSplineCurve3, ExactKnotVector};
use std::cmp::Ordering::{self, Equal, Greater, Less};

fn r(n: i64) -> R {
    R::from_integer(n.into())
}
fn q(n: i64, d: i64) -> R {
    R::new(n.into(), d.into())
}
fn pow2(n: i32) -> R {
    let value = BigInt::from(1) << n.unsigned_abs() as usize;
    if n < 0 {
        R::new(1.into(), value)
    } else {
        R::from_integer(value)
    }
}
fn multiply(a: &[R], b: &[R]) -> Vec<R> {
    let mut out = vec![r(0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}
fn choose(n: usize, k: usize) -> R {
    R::from_integer((0..k).fold(BigInt::from(1), |a, i| a * (n - i) / (i + 1)))
}
fn bernstein(p: &[R], degree: usize) -> Vec<R> {
    // t=-3+6u; direct binomial power-to-Bernstein identity, without de Boor.
    let mut power = vec![r(0); degree + 1];
    for (i, c) in p.iter().enumerate() {
        for (j, v) in power.iter_mut().enumerate().take(i + 1) {
            *v += c * choose(i, j) * r(-3).pow((i - j) as i32) * r(6).pow(j as i32);
        }
    }
    (0..=degree)
        .map(|i| {
            (0..=i)
                .map(|j| &power[j] * choose(i, j) / choose(degree, j))
                .sum()
        })
        .collect()
}
fn curve(
    p: usize,
    controls: Vec<[R; 4]>,
    knots: Vec<R>,
    mults: Vec<usize>,
    periodic: bool,
) -> ExactBSplineCurve3 {
    let axis = if periodic {
        ExactKnotVector::new_periodic(p, knots, mults)
    } else {
        ExactKnotVector::new(p, knots, mults)
    }
    .unwrap();
    ExactBSplineCurve3::from_homogeneous(axis, controls).unwrap()
}
fn distance(result: &SplineClosestSet, expected: &R) {
    assert_eq!(
        result.compare_squared_distance_exact(expected).unwrap(),
        Equal
    );
    bounds(result.squared_distance_bounds().unwrap(), expected);
    interval(result.distance_bounds().unwrap(), |x| {
        if x < &r(0) {
            Greater
        } else {
            expected.cmp(&(x * x))
        }
    });
    assert_eq!(result.compare_distance(-1.).unwrap(), Greater);
}

// 0..10 are sorted: -sqrt(5), -2, -sqrt(3), -sqrt(2), -1, 0,
// 1, sqrt(2), sqrt(3), 2, sqrt(5).
fn compare_known(index: usize, x: &R) -> Ordering {
    let (square, positive) = match index {
        0 => (5, false),
        2 => (3, false),
        3 => (2, false),
        7 => (2, true),
        8 => (3, true),
        10 => (5, true),
        _ => {
            return r(match index {
                1 => -2,
                4 => -1,
                5 => 0,
                6 => 1,
                9 => 2,
                _ => unreachable!(),
            })
            .cmp(x)
        }
    };
    if positive {
        if x < &r(0) {
            Greater
        } else {
            r(square).cmp(&(x * x))
        }
    } else if x > &r(0) {
        Less
    } else {
        (x * x).cmp(&r(square))
    }
}

fn factors(data: &[u8]) {
    let power_degree = 1 + usize::from(byte(data, 1) % 21);
    let degree = power_degree + 4;
    let mut p = vec![r(0), r(1)];
    let mut present = [false; 11];
    present[5] = true;
    while p.len() <= power_degree {
        let v = byte(data, 3 + p.len());
        let factor = if p.len() < power_degree && v & 1 == 1 {
            let (square, indices) = match (v / 2) % 3 {
                0 => (2, [3, 7]),
                1 => (3, [2, 8]),
                _ => (5, [0, 10]),
            };
            for i in indices {
                present[i] = true;
            }
            vec![-r(square), r(0), r(1)]
        } else {
            let index = [1, 4, 5, 6, 9][usize::from(v % 5)];
            present[index] = true;
            let value = [-2, -1, 0, 1, 2][usize::from(v % 5)];
            vec![-r(value), r(1)]
        };
        p = multiply(&p, &factor);
    }
    assert_eq!(p.len(), power_degree + 1);
    let base = [
        p.clone(),
        multiply(&p, &[r(0), r(0), r(1)]),
        multiply(&p, &[r(0), r(0), r(0), r(0), r(1)]),
    ];
    let components = base.map(|c| bernstein(&c, degree));
    let query: [R; 3] = std::array::from_fn(|i| q(i64::from(byte(data, 28 + i) as i8), 16));
    let shift = q(i64::from(byte(data, 31) as i8), 8);
    let width = pow2(i32::from(byte(data, 32) % 17) - 8);
    let scale = pow2(i32::from(byte(data, 33) % 33) - 16);
    let axis = usize::from(byte(data, 34) % 3);
    let controls = (0..=degree)
        .map(|i| {
            // Positive Bernstein weights, often full degree. Squared-distance
            // stationarity can reach 3p-2 although the complete zero set is known.
            let w = r(1 + i64::from(byte(data, 40 + i) % 5));
            let xyz: [R; 3] =
                std::array::from_fn(|c| &query[c] * &w + &scale * &components[(c + axis) % 3][i]);
            [xyz[0].clone(), xyz[1].clone(), xyz[2].clone(), w]
        })
        .collect();
    let mut c = curve(
        degree,
        controls,
        vec![&shift - r(3) * &width, &shift + r(3) * &width],
        vec![degree + 1; 2],
        false,
    );
    if degree <= 8 && byte(data, 35) & 1 != 0 {
        c = c.insert_knot(&(&shift + &width / r(3)), degree).unwrap();
        if byte(data, 35) & 2 != 0 {
            c = c.elevated(degree + 1).unwrap();
        }
    }
    let (lo, hi) = match byte(data, 36) % 4 {
        0 => (-3, 3),
        1 => (0, 3),
        2 => (-3, 0),
        _ => (0, 0),
    };
    let first = &shift + r(lo) * &width;
    let last = &shift + r(hi) * &width;
    let result =
        closest_points_on_exact_spline_in(&c, &query, &first, &last, Default::default()).unwrap();
    let expected: Vec<_> = (0..11)
        .filter(|&i| {
            present[i] && compare_known(i, &r(lo)) != Less && compare_known(i, &r(hi)) != Greater
        })
        .collect();
    assert_eq!(result.points().len(), expected.len());
    assert!(result.intervals().is_empty());
    distance(&result, &r(0));
    for (actual, i) in result.points().iter().zip(expected) {
        interval(actual.parameter_bounds().unwrap(), |x| {
            compare_known(i, &((x - &shift) / &width))
        });
        for (axis, value) in query.iter().enumerate() {
            assert_eq!(actual.compare_coordinate(axis, value).unwrap(), Equal);
        }
        for (b, value) in actual.coordinate_bounds().unwrap().into_iter().zip(&query) {
            bounds(b, value);
        }
    }
}

fn parabola(data: &[u8]) {
    let c = curve(
        2,
        vec![
            [r(-2), r(4), r(0), r(1)],
            [r(0), r(-4), r(0), r(1)],
            [r(2), r(4), r(0), r(1)],
        ],
        vec![r(-2), r(2)],
        vec![3; 2],
        false,
    );
    let axial = q(i64::from(byte(data, 2) % 9), 4);
    if byte(data, 1) & 1 == 0 {
        let square = q(i64::from(byte(data, 3) % 9), 4);
        let y = &square + q(1, 2);
        let result = closest_points_on_exact_spline_in(
            &c,
            &[r(0), y.clone(), axial.clone()],
            &r(-2),
            &r(2),
            Default::default(),
        )
        .unwrap();
        assert!(result.intervals().is_empty());
        assert_eq!(result.points().len(), if square == r(0) { 1 } else { 2 });
        distance(&result, &(y - q(1, 4) + &axial * &axial));
        for (i, point) in result.points().iter().enumerate() {
            let positive = i == 1;
            let compare = |x: &R| {
                if square == r(0) {
                    r(0).cmp(x)
                } else if positive {
                    if x < &r(0) {
                        Greater
                    } else {
                        square.cmp(&(x * x))
                    }
                } else if x > &r(0) {
                    Less
                } else {
                    (x * x).cmp(&square)
                }
            };
            interval(point.parameter_bounds().unwrap(), compare);
            let [x, y, z] = point.coordinate_bounds().unwrap();
            interval(x, compare);
            bounds(y, &square);
            bounds(z, &r(0));
        }
    } else {
        let epsilon = pow2(-8 - i32::from(byte(data, 3) % 128));
        let result = closest_points_on_exact_spline_in(
            &c,
            &[epsilon.clone(), r(1), axial.clone()],
            &r(-2),
            &r(2),
            Default::default(),
        )
        .unwrap();
        assert_eq!(result.points().len(), 1);
        assert!(result.intervals().is_empty());
        // Reflection makes the positive minimizer strictly better. On [1/2,1]
        // 2t^3-t-epsilon is strictly increasing and has exactly one root.
        let compare = |x: &R| {
            if x < &q(1, 2) {
                Greater
            } else if x > &r(1) {
                Less
            } else {
                (r(2) * x * x * x - x - &epsilon).cmp(&r(0)).reverse()
            }
        };
        let point = &result.points()[0];
        interval(point.parameter_bounds().unwrap(), compare);
        let [x, y, z] = point.coordinate_bounds().unwrap();
        interval(x, compare);
        bounds(z, &r(0));
        interval(y, |v| {
            // 2t(t^2-v) = (1-2v)t+epsilon at the positive root.
            let a = r(1) - r(2) * v;
            if a == r(0) {
                Greater
            } else {
                let order = compare(&(-&epsilon / &a));
                if a < r(0) {
                    order.reverse()
                } else {
                    order
                }
            }
        });
        let symmetric_distance = q(3, 4) + &axial * &axial;
        assert_eq!(
            result
                .compare_squared_distance_exact(&symmetric_distance)
                .unwrap(),
            Less
        );
        assert_eq!(
            result
                .compare_squared_distance_exact(&(&symmetric_distance - r(4) * &epsilon))
                .unwrap(),
            Greater
        );
        // Independent interval sign evaluation at the monotone cubic root.
        // Its rational root theorem rules out roots in (1/2,1) for these
        // epsilon powers, so a nonzero quadratic remainder cannot vanish.
        interval(result.squared_distance_bounds().unwrap(), |v| {
            let coefficients = [
                r(1) + &epsilon * &epsilon + &axial * &axial - v,
                -q(3, 2) * &epsilon,
                -q(1, 2),
            ];
            let (mut lo, mut hi) = (q(1, 2), r(1));
            for _ in 0..2048 {
                let (mut a, mut b) = (r(0), r(0));
                for c in coefficients.iter().rev() {
                    let products = [&a * &lo, &a * &hi, &b * &lo, &b * &hi];
                    a = products.iter().min().unwrap() + c;
                    b = products.iter().max().unwrap() + c;
                }
                if a > r(0) {
                    return Greater;
                }
                if b < r(0) {
                    return Less;
                }
                let mid = (&lo + &hi) / r(2);
                if compare(&mid) == Greater {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            panic!("independent cubic distance sign did not separate");
        });
    }
}

fn circles(data: &[u8]) {
    let periodic = byte(data, 1) & 1 == 1;
    let controls = if periodic {
        vec![
            [0, -2, 0, 2],
            [1, -1, 0, 1],
            [1, 0, 0, 1],
            [1, 1, 0, 1],
            [0, 2, 0, 2],
            [-1, 1, 0, 1],
            [-1, 0, 0, 1],
            [-1, -1, 0, 1],
        ]
    } else {
        vec![[1, 0, 0, 1], [1, 1, 0, 1], [0, 2, 0, 2]]
    };
    let c = curve(
        2,
        controls.into_iter().map(|c| c.map(r)).collect(),
        if periodic {
            (0..=4).map(r).collect()
        } else {
            vec![r(0), r(1)]
        },
        if periodic { vec![2; 5] } else { vec![3; 2] },
        periodic,
    );
    let (a, b) = if periodic {
        (
            q(i64::from(byte(data, 2) as i8), 8),
            r(1 + i64::from(byte(data, 3) % 12)),
        )
    } else {
        (q(i64::from(byte(data, 2) % 4), 8), q(1, 2))
    };
    let last = &a + &b;
    let z = q(i64::from(byte(data, 4) as i8), 16);
    let result = closest_points_on_exact_spline_in(
        &c,
        &[r(0), r(0), z.clone()],
        &a,
        &last,
        Default::default(),
    )
    .unwrap();
    assert!(result.points().is_empty());
    assert_eq!(result.intervals().len(), 1);
    assert_eq!(result.intervals()[0].parameters(), &[a, last]);
    distance(&result, &(r(1) + &z * &z));
}

fn dot(a: &[R; 3], b: &[R; 3]) -> R {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn polyline(data: &[u8]) {
    let n = 3 + usize::from(byte(data, 1) % 4);
    let periodic = byte(data, 2) & 1 != 0;
    let points: Vec<[R; 3]> = (0..n)
        .map(|i| std::array::from_fn(|j| q(i64::from(byte(data, 8 + 3 * i + j) as i8), 8)))
        .collect();
    let query: [R; 3] = std::array::from_fn(|j| q(i64::from(byte(data, 3 + j) as i8), 8));
    let end = if periodic { n } else { n - 1 };
    let knots = (0..=end).map(|i| r(i as i64)).collect();
    let mut mults = vec![1; end + 1];
    if !periodic {
        mults[0] = 2;
        mults[end] = 2;
    }
    let c = curve(
        1,
        points
            .iter()
            .map(|p| [p[0].clone(), p[1].clone(), p[2].clone(), r(1)])
            .collect(),
        knots,
        mults,
        periodic,
    );
    let first = if periodic {
        -r(end as i64) + q(1, 3)
    } else {
        q(1, 3)
    };
    let last = if periodic {
        r(2 * end as i64) - q(1, 7)
    } else {
        r(end as i64) - q(1, 7)
    };
    let mut minimum: Option<R> = None;
    let mut expected: Vec<(R, [R; 3])> = vec![];
    let mut intervals: Vec<[R; 2]> = vec![];
    for j in -(end as i64)..2 * end as i64 {
        if !periodic && (j < 0 || j >= end as i64) {
            continue;
        }
        let lo = first.clone().max(r(j));
        let hi = last.clone().min(r(j + 1));
        if lo >= hi {
            continue;
        }
        let index = j.rem_euclid(n as i64) as usize;
        let a = &points[index];
        let b = &points[(index + 1) % n];
        let direction: [R; 3] = std::array::from_fn(|i| &b[i] - &a[i]);
        let delta: [R; 3] = std::array::from_fn(|i| &query[i] - &a[i]);
        let denominator = dot(&direction, &direction);
        let u = if denominator == r(0) {
            &lo - r(j)
        } else {
            (dot(&delta, &direction) / &denominator)
                .max(&lo - r(j))
                .min(&hi - r(j))
        };
        let at: [R; 3] = std::array::from_fn(|i| &a[i] + &u * &direction[i]);
        let error: [R; 3] = std::array::from_fn(|i| &at[i] - &query[i]);
        let d = dot(&error, &error);
        let order = minimum.as_ref().map_or(Less, |v| d.cmp(v));
        if order == Less {
            minimum = Some(d);
            expected.clear();
            intervals.clear();
        }
        if order != Greater {
            if denominator == r(0) {
                intervals.push([lo, hi]);
            } else {
                expected.push((r(j) + u, at));
            }
        }
    }
    let mut merged: Vec<[R; 2]> = vec![];
    for interval in intervals {
        if let Some(last) = merged.last_mut().filter(|s| s[1] == interval[0]) {
            last[1] = interval[1].clone();
        } else {
            merged.push(interval);
        }
    }
    expected.sort_by(|a, b| a.0.cmp(&b.0));
    expected.dedup_by(|a, b| a.0 == b.0);
    expected.retain(|(u, _)| !merged.iter().any(|s| s[0] <= *u && *u <= s[1]));
    let result =
        closest_points_on_exact_spline_in(&c, &query, &first, &last, Default::default()).unwrap();
    assert_eq!(result.points().len(), expected.len());
    assert_eq!(result.intervals().len(), merged.len());
    distance(&result, &minimum.unwrap());
    for (got, (u, xyz)) in result.points().iter().zip(expected) {
        assert_eq!(got.compare_parameter(&u).unwrap(), Equal);
        bounds(got.parameter_bounds().unwrap(), &u);
        for (i, x) in xyz.iter().enumerate() {
            assert_eq!(got.compare_coordinate(i, x).unwrap(), Equal);
        }
        for (b, x) in got.coordinate_bounds().unwrap().into_iter().zip(xyz) {
            bounds(b, &x);
        }
    }
    for (got, range) in result.intervals().iter().zip(merged) {
        assert_eq!(got.parameters(), &range);
    }
}

pub fn check_spline_proximity(data: &[u8]) {
    match byte(data, 0) % 4 {
        0 => factors(data),
        1 => parabola(data),
        2 => circles(data),
        _ => polyline(data),
    }
}
