//! Complete known-factor Bezier contacts and independent piecewise rational
//! line/plane and weighted-line/quadric intersections, including periodic seams
//! and contained runs. Quadric tangency orders come from squared known factors.
use super::{byte, one, zero};
use crate::splines::{bounds, rat};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::intersection::{
    spline_cylinder_in, spline_plane, spline_plane_in, spline_sphere_in, Cylinder3, Plane3,
    Sphere3, SplinePlaneContact as C,
};
use rusty_occt::{BSplineCurve3, Point3, Vec3};
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
    let quadric = byte(data, 0) & 16 != 0;
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
    let (mut frame, scale) = frame(data);
    if quadric {
        frame = if byte(data, 3) & 1 == 0 {
            [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        } else {
            [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]]
        };
    }
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
        let mut p = [
            f64::from(byte(data, 50 + 2 * i) as i8),
            f64::from(byte(data, 51 + 2 * i) as i8),
            z / w,
        ];
        if quadric {
            p = [z / w, 1., 0.];
        }
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
    let clipped = byte(data, 0) & 4 != 0;
    let (mut first, last) = if clipped {
        let a = byte(data, 90) % 4;
        let b = a + 1 + byte(data, 91) % (4 - a);
        (
            (f64::from(a) - 1.) * parameter_scale,
            (f64::from(b) - 1.) * parameter_scale,
        )
    } else {
        curve.domain()
    };
    if clipped && byte(data, 0) & 8 != 0 {
        first = first.next_up();
    }
    let lower = (rat(first) + rat(parameter_scale)) / rat(4. * parameter_scale);
    let upper = (rat(last) + rat(parameter_scale)) / rat(4. * parameter_scale);
    expected.retain(|r, _| r >= &lower && r <= &upper);
    let result = if quadric {
        if byte(data, 3) & 2 == 0 {
            spline_sphere_in(
                &curve,
                &Sphere3::new(Point3::new(0., 0., 0.), scale).unwrap(),
                first,
                last,
            )
        } else {
            let [x, y, z] = frame[2];
            spline_cylinder_in(
                &curve,
                &Cylinder3::new(Point3::new(0., 0., 0.), Vec3::new(x, y, z), scale).unwrap(),
                first,
                last,
            )
        }
    } else if clipped {
        spline_plane_in(&curve, &plane(frame, scale), first, last)
    } else {
        spline_plane(&curve, &plane(frame, scale))
    }
    .unwrap();
    assert!(result.overlaps().is_empty());
    assert_eq!(result.points().len(), expected.len());
    for (point, (root, order)) in result.points().iter().zip(expected) {
        let order = if quadric { 2 * order } else { order };
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
                (root > lower).then_some(order),
                (root < upper).then_some(order)
            ]
        );
        assert_eq!(
            point.contact(),
            if root == lower || root == upper {
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
    let clipped = byte(data, 0) & 4 != 0;
    let period = knots[spans];
    let (first, last) = if clipped {
        let offset = if periodic {
            if byte(data, 0) & 8 != 0 {
                2f64.powi(53) * period
            } else {
                (f64::from(byte(data, 90) % 5) - 2.) * period
            }
        } else {
            0.
        };
        let a = f64::from(byte(data, 91) % 4) * period / 4.;
        let length = if periodic {
            f64::from(1 + byte(data, 92) % 8) * period / 4.
        } else {
            period - a
        };
        let first = offset + a;
        (first, (first + length).max(first.next_up()))
    } else {
        curve.domain()
    };
    let (lower, upper) = (rat(first), rat(last));
    let start_turn = if periodic {
        (&lower / rat(period)).floor()
    } else {
        zero()
    };
    let turns = if periodic {
        ((&upper / rat(period)).ceil() - &start_turn)
            .to_integer()
            .to_string()
            .parse::<usize>()
            .unwrap()
    } else {
        1
    };
    assert!(turns <= 16);
    let mut hits = BTreeMap::<R, Hit>::new();
    let mut overlaps: Vec<(R, R)> = Vec::new();
    for (turn, i) in (0..turns).flat_map(|t| (0..spans).map(move |i| (t, i))) {
        let offset = (&start_turn + integer(turn as i64)) * rat(period);
        let start = rat(knots[i]) + &offset;
        let end = rat(knots[i + 1]) + offset;
        let low = start.clone().max(lower.clone());
        let high = end.clone().min(upper.clone());
        if low >= high {
            continue;
        }
        let j = (i + 1) % np;
        let (a, b) = (&heights[i], &heights[j]);
        if a == &zero() && b == &zero() {
            if let Some(last) = overlaps.last_mut().filter(|last| last.1 == low) {
                last.1 = high;
            } else {
                overlaps.push((low, high));
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
        let u = &start + (&end - &start) * &root;
        if u < low || u > high {
            continue;
        }
        let w = (one() - &root) * rat(weights[i]) + &root * rat(weights[j]);
        let position = std::array::from_fn(|c| {
            ((one() - &root) * rat(weights[i]) * rat(poles[i].to_array()[c])
                + &root * rat(weights[j]) * rat(poles[j].to_array()[c]))
                / &w
        });
        let signs = [
            (u > low).then(|| a.cmp(&zero())),
            (u < high).then(|| b.cmp(&zero())),
        ];
        if let Some(previous) = hits.get_mut(&u) {
            assert_eq!(previous.position, position);
            previous.signs[1] = signs[1];
        } else {
            hits.insert(u, Hit { position, signs });
        }
    }
    hits.retain(|u, _| !overlaps.iter().any(|(a, b)| a <= u && u <= b));
    let result = if clipped {
        spline_plane_in(&curve, &plane(frame, scale), first, last)
    } else {
        spline_plane(&curve, &plane(frame, scale))
    }
    .unwrap();
    assert_eq!(result.overlaps().len(), overlaps.len());
    for (o, (a, b)) in result.overlaps().iter().zip(overlaps) {
        let [lo, hi] = o.parameter_bounds();
        bounds(lo, &a);
        bounds(hi, &b);
    }
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
// Geometric line roots come from a rational projection/quadratic-sign oracle.
// Convert the physical affine parameter through the independent rational
// degree-one weight map; production instead solves its homogeneous polynomial.
fn quadric_line(data: &[u8]) {
    use crate::curved::{affine_compare, expected_surface, interval};
    let scale = 2f64.powi((i32::from(byte(data, 4)) - 128) * 4);
    let center = Point3::new(0., 0., 0.);
    let p = Point3::new(
        f64::from(byte(data, 10) as i8) * scale,
        f64::from(byte(data, 11) as i8) * scale,
        f64::from(byte(data, 12) as i8) * scale,
    );
    let q = Point3::new(
        f64::from(byte(data, 13) as i8) * scale,
        f64::from(byte(data, 14) as i8) * scale,
        f64::from(byte(data, 15) as i8) * scale,
    );
    let axis = Vec3::new(
        1.,
        f64::from(byte(data, 16) as i8),
        f64::from(byte(data, 17) as i8),
    );
    let radius = f64::from(1 + byte(data, 18) % 128) * scale;
    let weights = [
        2f64.powi(i32::from(byte(data, 19) % 17) - 8),
        2f64.powi(i32::from(byte(data, 20) % 17) - 8),
    ];
    let curve = BSplineCurve3::new(
        1,
        vec![p, q],
        Some(weights.to_vec()),
        vec![0., 1.],
        vec![2, 2],
    )
    .unwrap();
    let cylinder = byte(data, 3) & 2 != 0;
    let (first, last) = if byte(data, 0) & 4 != 0 {
        (0.25, 0.75)
    } else {
        (0., 1.)
    };
    let result = if cylinder {
        spline_cylinder_in(
            &curve,
            &Cylinder3::new(center, axis, radius).unwrap(),
            first,
            last,
        )
    } else {
        spline_sphere_in(&curve, &Sphere3::new(center, radius).unwrap(), first, last)
    }
    .unwrap();
    let mut input = Vec::new();
    input.extend(center.to_array().map(rat));
    input.extend(axis.to_array().map(rat));
    input.push(rat(radius));
    input.extend(p.to_array().map(rat));
    input.extend(q.to_array().map(rat));
    let w = weights.map(rat);
    let physical = |t: &R| t * &w[1] / ((one() - t) * &w[0] + t * &w[1]);
    let value = |p: Point3| {
        let v = p.to_array().map(rat);
        let a = axis.to_array().map(rat);
        let dot: R = v.iter().zip(&a).map(|(x, y)| x * y).sum();
        v.iter().map(|x| x * x).sum::<R>()
            - if cylinder {
                &dot * &dot / a.iter().map(|x| x * x).sum::<R>()
            } else {
                zero()
            }
            - rat(radius) * rat(radius)
    };
    // Constant geometric curves on the surface occupy their entire parameter
    // interval; a zero-length segment API has a different point contract.
    if p == q && value(p) == zero() {
        assert!(result.points().is_empty());
        assert_eq!(result.overlaps().len(), 1);
        assert_eq!(result.overlaps()[0].parameters(), (first, last));
        return;
    }
    let expected = expected_surface(if cylinder { 2 } else { 1 }, &input, true);
    let Some(expected) = expected else {
        assert!(result.points().is_empty());
        assert_eq!(result.overlaps().len(), 1);
        assert_eq!(result.overlaps()[0].parameters(), (first, last));
        return;
    };
    let low = physical(&rat(first));
    let high = physical(&rat(last));
    let expected: Vec<_> = expected
        .into_iter()
        .filter(|(r, _)| r.compare(&low) != Ordering::Less && r.compare(&high) != Ordering::Greater)
        .collect();
    assert!(result.overlaps().is_empty());
    assert_eq!(result.points().len(), expected.len());
    for (point, (root, _)) in result.points().iter().zip(expected) {
        interval(point.parameter(), |t| {
            if t < &zero() {
                Ordering::Greater
            } else if t > &one() {
                Ordering::Less
            } else {
                root.compare(&physical(t))
            }
        });
        for c in 0..3 {
            let offset = rat(p.to_array()[c]);
            let direction = rat(q.to_array()[c]) - &offset;
            interval(point.coordinate_bounds()[c], |x| {
                affine_compare(&root, &offset, &direction, x)
            });
        }
        let at_start = root.compare(&low) == Ordering::Equal;
        let at_end = root.compare(&high) == Ordering::Equal;
        let order = usize::from(root.multiplicity());
        assert_eq!(
            point.multiplicities(),
            [(!at_start).then_some(order), (!at_end).then_some(order)]
        );
        assert_eq!(
            point.contact(),
            if at_start || at_end {
                C::Boundary
            } else if order == 2 {
                C::Tangent
            } else {
                C::Crossing
            }
        );
    }
}

pub fn check_spline_intersections(data: &[u8]) {
    if byte(data, 0) & 48 == 48 {
        power_quadric(data)
    } else if byte(data, 0) & 1 == 0 {
        bezier(data)
    } else if byte(data, 0) & 16 != 0 {
        quadric_line(data)
    } else {
        polygon(data)
    }
}

// Degree-50 equations have one known positive root: 2*t^degree=1. Monotonic
// rational powers certify it independently of the production Sturm isolator.
fn power_quadric(data: &[u8]) {
    use crate::curved::interval;
    let degree = 1 + usize::from(byte(data, 1) % 25);
    let scale = 2f64.powi((i32::from(byte(data, 4)) - 128) * 4);
    let mut poles = vec![Point3::new(0., 0., 0.); degree + 1];
    poles[degree] = Point3::new(2. * scale, 0., 0.);
    let curve = BSplineCurve3::new(degree, poles, None, vec![0., 1.], vec![degree + 1; 2]).unwrap();
    let first = if byte(data, 0) & 4 == 0 { 0. } else { 0.75 };
    let last = if byte(data, 0) & 8 == 0 { 1. } else { 0.875 };
    let result = if byte(data, 3) & 2 == 0 {
        spline_sphere_in(
            &curve,
            &Sphere3::new(Point3::new(0., 0., 0.), scale).unwrap(),
            first,
            last,
        )
    } else {
        spline_cylinder_in(
            &curve,
            &Cylinder3::new(Point3::new(0., 0., 0.), Vec3::new(0., 0., 1.), scale).unwrap(),
            first,
            last,
        )
    }
    .unwrap();
    let compare = |t: &R| {
        if t < &zero() {
            Ordering::Greater
        } else if t > &one() {
            Ordering::Less
        } else {
            (integer(2) * t.pow(degree as i32) - one())
                .cmp(&zero())
                .reverse()
        }
    };
    assert!(result.overlaps().is_empty());
    if compare(&rat(first)) == Ordering::Less || compare(&rat(last)) == Ordering::Greater {
        assert!(result.points().is_empty());
        return;
    }
    assert_eq!(result.points().len(), 1);
    let point = &result.points()[0];
    interval(point.parameter(), compare);
    for (b, x) in point
        .coordinate_bounds()
        .into_iter()
        .zip([rat(scale), zero(), zero()])
    {
        bounds(b, &x);
    }
    let start = compare(&rat(first)) == Ordering::Equal;
    let end = compare(&rat(last)) == Ordering::Equal;
    assert_eq!(
        point.multiplicities(),
        [(!start).then_some(1), (!end).then_some(1)]
    );
    assert_eq!(
        point.contact(),
        if start || end {
            C::Boundary
        } else {
            C::Crossing
        }
    );
}
