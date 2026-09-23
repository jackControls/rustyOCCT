//! Complete known-factor, rational secant and periodic certificates for exact
//! edited curves. Expected roots are constructed algebraically, not recovered
//! with production de Boor interpolation or its root isolator.
use crate::byte;
use crate::knot_reference as identity;
use crate::splines::{bounds, rat};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::intersection::*;
use rusty_occt::{Error, ExactBSplineCurve3, ExactKnotVector, Point3, ScalarInterval, Vec3};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

fn r(n: i64) -> R {
    R::from_integer(n.into())
}
fn q(n: i64, d: i64) -> R {
    R::new(n.into(), d.into())
}
fn pow2(n: i32) -> R {
    let x = BigInt::from(1) << n.unsigned_abs() as usize;
    if n < 0 {
        R::new(1.into(), x)
    } else {
        R::from_integer(x)
    }
}
fn choose(n: usize, k: usize) -> R {
    R::from_integer((0..k).fold(BigInt::from(1), |a, i| a * (n - i) / (i + 1)))
}
fn rotate<T: Clone>(p: &[T; 3], axis: usize) -> [T; 3] {
    std::array::from_fn(|c| p[(c + axis) % 3].clone())
}
fn point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}
fn primitive_plane(axis: usize) -> Plane3 {
    Plane3::through_points(
        Point3::new(0., 0., 0.),
        point(rotate(&[1., 0., 0.], axis)),
        point(rotate(&[0., 1., 0.], axis)),
    )
    .unwrap()
}
fn primitive_cylinder(axis: usize) -> Cylinder3 {
    let [x, y, z] = rotate(&[0., 0., 1.], axis);
    Cylinder3::new(Point3::new(0., 0., 0.), Vec3::new(x, y, z), 1.).unwrap()
}
fn intersect(
    curve: &ExactBSplineCurve3,
    family: usize,
    axis: usize,
    first: &R,
    last: &R,
    options: SplineSurfaceOptions,
) -> rusty_occt::Result<ExactSplineSurfaceIntersection> {
    match family {
        0 => {
            exact_spline_plane_in_with_options(curve, &primitive_plane(axis), first, last, options)
        }
        1 => exact_spline_sphere_in_with_options(
            curve,
            &Sphere3::new(Point3::new(0., 0., 0.), 1.).unwrap(),
            first,
            last,
            options,
        ),
        2 => exact_spline_cylinder_in_with_options(
            curve,
            &primitive_cylinder(axis),
            first,
            last,
            options,
        ),
        _ => unreachable!(),
    }
}
fn enclosure(actual: rusty_occt::Result<ScalarInterval>, value: &R) -> bool {
    if value < &-rat(f64::MAX) || value > &rat(f64::MAX) {
        assert!(matches!(actual, Err(Error::Unrepresentable(_))));
        false
    } else {
        bounds(actual.unwrap(), value);
        true
    }
}

pub fn check_exact_spline_intersections(data: &[u8]) {
    check(data, |_| {});
}
pub fn profile_exact_spline_intersections(
    data: &[u8],
) -> BTreeMap<&'static str, std::time::Duration> {
    let mut times = BTreeMap::new();
    let mut at = std::time::Instant::now();
    check(data, |label| {
        *times.entry(label).or_default() += at.elapsed();
        at = std::time::Instant::now();
    });
    times
}

fn check(data: &[u8], mut mark: impl FnMut(&'static str)) {
    let family = usize::from(byte(data, 0) % 3);
    let mode = byte(data, 2) % 6;
    let degree = if mode == 5 {
        2
    } else {
        1 + usize::from(byte(data, 1) % 25)
    };
    let axis = usize::from(byte(data, 4) % 3);
    let (a, width) = match byte(data, 3) % 8 {
        0 => (q(1, 3), q(5, 7)),
        1 => (pow2(2048), pow2(-2048)),
        2 => (-pow2(2048), q(7, 3)),
        3 => (-pow2(-2048), pow2(-2047)),
        4 => (
            q(i64::from(byte(data, 4) as i8), 257),
            pow2(i32::from(byte(data, 5) as i8)),
        ),
        5 => (r(0), r(1)),
        6 => (pow2(1024), pow2(1024)),
        _ => (-pow2(1024), pow2(1025)),
    };
    let periodic = mode == 3;
    let (lo, hi) = if periodic {
        let lo = q(i64::from(byte(data, 7) % 9) - 4, 2);
        let period = if family == 0 { 3 } else { 4 };
        (
            lo.clone(),
            lo + r(period * i64::from(1 + byte(data, 9) % 3)),
        )
    } else if byte(data, 7) & 8 != 0 {
        (q(1, 7), q(2, 3))
    } else {
        let lo = q(i64::from(byte(data, 7) % 8), 8);
        let hi = &lo + (r(1) - &lo) * q(i64::from(1 + byte(data, 8) % 8), 8);
        (lo, hi)
    };
    let first = &a + &width * &lo;
    let last = &a + &width * &hi;
    let mut expected = BTreeMap::<R, ([R; 3], usize)>::new();
    let mut overlap = false;
    let (controls, knots, mults, p) = if periodic {
        if family == 0 {
            for k in -4..=10 {
                let t = q(3 * k, 2);
                if lo <= t && t <= hi {
                    expected.insert(&a + &width * &t, ([r(0), r(0), r(0)], 1));
                }
            }
            (
                vec![
                    [r(0), r(0), r(0), r(1)],
                    [r(1), r(0), r(1), r(1)],
                    [-r(1), r(0), -r(1), r(1)],
                ],
                (0..=3).map(|i| &a + &width * r(i)).collect(),
                vec![1; 4],
                1,
            )
        } else {
            overlap = true;
            // Four rational circle quadrants with shared positive endpoint
            // weights; every homogeneous span satisfies X^2+Y^2=W^2.
            let controls = [
                [0, -2, 2],
                [1, -1, 1],
                [1, 0, 1],
                [1, 1, 1],
                [0, 2, 2],
                [-1, 1, 1],
                [-1, 0, 1],
                [-1, -1, 1],
            ]
            .map(|[x, y, w]| [r(x), r(y), r(0), r(w)])
            .to_vec();
            (
                controls,
                (0..=4).map(|i| &a + &width * r(i)).collect(),
                vec![2; 5],
                2,
            )
        }
    } else if mode == 2 {
        overlap = true;
        let controls = if family == 0 {
            vec![[r(0), r(0), r(0), r(1)]; 2]
        } else {
            vec![
                [r(1), r(0), r(0), r(1)],
                [r(1), r(1), r(0), r(1)],
                [r(0), r(2), r(0), r(2)],
            ]
        };
        let p = controls.len() - 1;
        (controls, vec![a.clone(), &a + &width], vec![p + 1; 2], p)
    } else if mode == 4 {
        let scale = r(i64::from(1 + byte(data, 5) % 8));
        let axial = if byte(data, 6) & 1 != 0 {
            pow2(2048)
        } else {
            r(1)
        };
        let roots = if family == 0 {
            vec![q(1, 2)]
        } else {
            vec![
                (&scale - r(1)) / (r(2) * &scale + r(1)),
                (&scale + r(1)) / (r(2) * &scale - r(1)),
            ]
        };
        for t in roots {
            if lo <= t && t <= hi {
                let w = r(1) + &t;
                let x = &scale * (r(2) * &t - r(1)) / &w;
                let coordinates = if family == 0 {
                    [&t / &w, r(0), r(0)]
                } else {
                    [x, r(0), if family == 1 { r(0) } else { &t * &axial / &w }]
                };
                expected.insert(&a + &width * &t, (coordinates, 1));
            }
        }
        let controls = (0..=degree)
            .map(|i| {
                let t = R::new(i.into(), degree.into());
                let w = r(1) + &t;
                let n = &scale * (r(2) * &t - r(1));
                if family == 0 {
                    [t, r(0), n, w]
                } else {
                    [n, r(0), if family == 1 { r(0) } else { t * &axial }, w]
                }
            })
            .collect();
        (
            controls,
            vec![a.clone(), &a + &width],
            vec![degree + 1; 2],
            degree,
        )
    } else {
        let denominator = i64::from(2 + byte(data, 10) % 8);
        let roots: Vec<_> = if mode == 5 {
            vec![
                q(1, 3),
                q(1, 3) + pow2(-(1 + i32::from(byte(data, 11))) * 8),
            ]
        } else {
            (0..degree)
                .map(|i| {
                    q(
                        i64::from(
                            byte(data, 16 + if mode == 1 { 0 } else { i })
                                % ((denominator + 3) as u8),
                        ) - 1,
                        denominator,
                    )
                })
                .collect()
        };
        let mut polynomial = vec![r(1)];
        for root in &roots {
            let mut next = vec![r(0); polynomial.len() + 1];
            for (i, c) in polynomial.iter().enumerate() {
                next[i] -= c * root;
                next[i + 1] += c;
            }
            polynomial = next;
        }
        let weight_scale = if byte(data, 6) & 1 != 0 {
            pow2(-2048)
        } else {
            r(1)
        };
        let weights: Vec<_> = (0..=degree)
            .map(|i| q(i64::from(1 + byte(data, 64 + i) % 31), 11) * &weight_scale)
            .collect();
        let distinct: BTreeSet<_> = roots.iter().cloned().collect();
        for t in distinct {
            if lo <= t && t <= hi {
                let w: R = weights
                    .iter()
                    .enumerate()
                    .map(|(i, w)| {
                        w * choose(degree, i)
                            * t.pow(i as i32)
                            * (r(1) - &t).pow((degree - i) as i32)
                    })
                    .sum();
                let coordinates = if family == 0 {
                    [&t / &w, r(0), r(0)]
                } else {
                    [r(0), r(1), if family == 1 { r(0) } else { &t / &w }]
                };
                let order =
                    roots.iter().filter(|x| **x == t).count() * if family == 0 { 1 } else { 2 };
                expected.insert(&a + &width * &t, (coordinates, order));
            }
        }
        let controls = (0..=degree)
            .map(|i| {
                let n: R = (0..=i)
                    .map(|j| &polynomial[j] * choose(i, j) / choose(degree, j))
                    .sum();
                let t = R::new(i.into(), degree.into());
                let w = weights[i].clone();
                if family == 0 {
                    [t, r(0), n, w]
                } else {
                    [n, w.clone(), if family == 1 { r(0) } else { t }, w]
                }
            })
            .collect();
        (
            controls,
            vec![a.clone(), &a + &width],
            vec![degree + 1; 2],
            degree,
        )
    };
    let common = pow2(i32::from(byte(data, 12) as i8) * 8);
    let controls = controls
        .into_iter()
        .map(|p: [R; 4]| {
            let [x, y, z] = rotate(&[p[0].clone(), p[1].clone(), p[2].clone()], axis);
            [x * &common, y * &common, z * &common, &p[3] * &common]
        })
        .collect();
    let basis = if periodic {
        ExactKnotVector::new_periodic(p, knots, mults)
    } else {
        ExactKnotVector::new(p, knots, mults)
    }
    .unwrap();
    let original = ExactBSplineCurve3::from_homogeneous(basis, controls).unwrap();
    mark("input construction");
    let b = original
        .knots()
        .iter()
        .find(|k| *k > original.domain().first().unwrap())
        .unwrap();
    let cut = &a + (b - &a) * q(1 + i64::from(byte(data, 13)), 257);
    let next = (&cut + b) / r(2);
    let edited = match byte(data, 8) % 4 {
        0 => original.clone(),
        1 => original
            .refined(&[(cut.clone(), p), (next.clone(), 1), (cut.clone(), 0)])
            .unwrap(),
        2 => {
            let fine = original.insert_knot(&cut, p).unwrap();
            let restored = fine.remove_knot(&cut, 0).unwrap().unwrap();
            assert_eq!(restored, original);
            restored
        }
        _ if periodic => {
            let removed = original
                .remove_knot(&a, original.multiplicities()[0] - 1)
                .unwrap();
            assert_eq!(
                removed,
                identity::removed(&original, &a, original.multiplicities()[0] - 1)
            );
            removed.unwrap_or_else(|| original.clone())
        }
        _ => original
            .refined(&[(cut.clone(), p), (next, 1)])
            .unwrap()
            .remove_knot(&cut, 0)
            .unwrap()
            .unwrap(),
    };
    mark("kernel edits");
    assert!(identity::equal(&original, &edited));
    mark("oracle complete identity");
    let hits = intersect(
        &edited,
        family,
        axis,
        &first,
        &last,
        SplineSurfaceOptions::default(),
    )
    .unwrap();
    mark("kernel intersections");
    assert_eq!(hits.points().len(), expected.len());
    assert_eq!(hits.overlaps().len(), usize::from(overlap));
    assert_eq!(hits.is_disjoint(), expected.is_empty() && !overlap);
    for (p, (parameter, (coordinates, order))) in hits.points().iter().zip(&expected) {
        assert_eq!(p.compare_parameter(parameter).unwrap(), Ordering::Equal);
        assert_eq!(
            p.compare_parameter(&(parameter - q(1, 257))).unwrap(),
            Ordering::Greater
        );
        assert_eq!(
            p.compare_parameter(&(parameter + q(1, 257))).unwrap(),
            Ordering::Less
        );
        mark("kernel parameter comparisons");
        let boundary = parameter == &first || parameter == &last;
        assert_eq!(
            p.contact(),
            if boundary {
                SplineSurfaceContact::Boundary
            } else if order % 2 == 1 {
                SplineSurfaceContact::Crossing
            } else {
                SplineSurfaceContact::Tangent
            }
        );
        assert_eq!(
            p.multiplicities(),
            [
                (parameter != &first).then_some(*order),
                (parameter != &last).then_some(*order)
            ]
        );
        let coordinates = rotate(coordinates, axis);
        enclosure(p.parameter_bounds(), parameter);
        mark("parameter bounds and contact orders");
        for (i, value) in coordinates.iter().enumerate() {
            assert_eq!(p.compare_coordinate(i, value).unwrap(), Ordering::Equal);
            mark("kernel coordinate comparisons");
            enclosure(p.coordinate_bound(i), value);
            mark("coordinate bounds and oracle");
        }
    }
    if overlap {
        let interval = &hits.overlaps()[0];
        assert_eq!(interval.parameters(), &[first.clone(), last.clone()]);
        for (i, value) in [&first, &last].into_iter().enumerate() {
            assert_eq!(
                interval.compare_parameter(i, value).unwrap(),
                Ordering::Equal
            );
            assert_eq!(
                interval.compare_parameter(i, &(value - q(1, 257))).unwrap(),
                Ordering::Greater
            );
        }
        if [&first, &last]
            .into_iter()
            .any(|v| v < &-rat(f64::MAX) || v > &rat(f64::MAX))
        {
            assert!(matches!(
                interval.parameter_bounds(),
                Err(Error::Unrepresentable(_))
            ));
            assert!(matches!(hits.enclosed(), Err(Error::Unrepresentable(_))));
        } else {
            for (bound, value) in interval
                .parameter_bounds()
                .unwrap()
                .into_iter()
                .zip([&first, &last])
            {
                bounds(bound, value);
            }
            let enclosed = hits.enclosed().unwrap();
            assert!(enclosed.points().is_empty());
            assert_eq!(enclosed.overlaps().len(), 1);
        }
    }
    mark("exact contacts and bounds");
    let bad = R::new_raw(1.into(), 0.into());
    assert!(intersect(
        &edited,
        family,
        axis,
        &bad,
        &last,
        SplineSurfaceOptions::default()
    )
    .is_err());
    assert!(intersect(
        &edited,
        family,
        axis,
        &first,
        &first,
        SplineSurfaceOptions::default()
    )
    .is_err());
    assert!(matches!(
        intersect(
            &edited,
            family,
            axis,
            &first,
            &last,
            SplineSurfaceOptions {
                max_spans: 0,
                ..Default::default()
            }
        ),
        Err(Error::ComputationLimit(_))
    ));
    if let Some(p) = hits.points().first() {
        assert!(p.compare_parameter(&bad).is_err());
        assert!(p.coordinate_bound(3).is_err());
    }
    if let Some(interval) = hits.overlaps().first() {
        assert!(interval.compare_parameter(0, &bad).is_err());
        assert!(matches!(
            interval.compare_parameter(2, &first),
            Err(Error::OutOfDomain(_))
        ));
    }
    mark("rejection checks");
}
