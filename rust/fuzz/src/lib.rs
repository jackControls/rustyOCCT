//! Fuzz-only, independently arranged rational matrix / barycentric oracles.
//! This crate is never linked into the kernel or ordinary application builds.
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::intersection::{
    line_plane, segment_plane, segment_triangle, IntersectionPoint, LinePlaneIntersection as L,
    SegmentPlaneIntersection as P, SegmentTriangleIntersection as T,
};
use rusty_occt::predicates::{
    in_sphere, orient2d, orient3d, Orientation2, Orientation3, SphereLocation,
};
use rusty_occt::{Error, Point2, Point3, Triangle3};
use std::cmp::Ordering;

mod curved;
pub use curved::check_curved;
mod splines;
pub use splines::check_splines;
mod surfaces;
pub use surfaces::check_surfaces;

pub fn byte(data: &[u8], i: usize) -> u8 {
    data.get(i).copied().unwrap_or(0)
}

/// Raw binary64, scaled integers, coplanar and sphere-boundary modes. Short
/// inputs are zero-padded, so minimization can remove irrelevant suffixes.
pub fn decode(data: &[u8]) -> [Point3; 5] {
    let mode = byte(data, 0) % 4;
    let mut p: [[f64; 3]; 5] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let offset = 1 + (i * 3 + j) * 8;
            let word = u64::from_le_bytes(std::array::from_fn(|k| byte(data, offset + k)));
            if mode == 0 {
                f64::from_bits(word)
            } else {
                (word as i16 as f64) * 2.0_f64.powi(byte(data, 120) as i32 - 128)
            }
        })
    });
    if mode == 2 {
        for p in &mut p {
            p[2] = 0.;
        }
    }
    if mode == 3 {
        p[4] = p[0];
    }
    p.map(|[x, y, z]| Point3::new(x, y, z))
}

fn zero() -> R {
    R::from_integer(BigInt::from(0))
}
fn one() -> R {
    R::from_integer(BigInt::from(1))
}
fn rational(p: Point3) -> [R; 3] {
    p.to_array().map(|x| R::from_float(x).unwrap())
}

/// Rational Gaussian elimination, distinct from the production determinant
/// expansions and their fixed binary64 integer lattice.
fn determinant(mut a: Vec<Vec<R>>) -> R {
    let mut det = one();
    for i in 0..a.len() {
        let Some(pivot) = (i..a.len()).find(|&j| a[j][i] != zero()) else {
            return zero();
        };
        if pivot != i {
            a.swap(i, pivot);
            det = -det;
        }
        det *= &a[i][i];
        for j in i + 1..a.len() {
            let (prior, current) = a.split_at_mut(j);
            let (pivot, row) = (&prior[i], &mut current[0]);
            let ratio = &row[i] / &pivot[i];
            for (cell, value) in row.iter_mut().zip(pivot).skip(i + 1) {
                *cell -= &ratio * value;
            }
        }
    }
    det
}

fn orientation(p: &[[R; 3]]) -> R {
    determinant(
        p.iter()
            .map(|p| [vec![one()], p.to_vec()].concat())
            .collect(),
    )
}

pub fn check_predicates(data: &[u8]) {
    let p = decode(data);
    let o = orient3d(p[0], p[1], p[2], p[3]);
    let s = in_sphere(p[0], p[1], p[2], p[3], p[4]);
    if !p.iter().flat_map(|p| p.to_array()).all(f64::is_finite) {
        assert!(matches!(s, Err(Error::NonFinite(_))));
        if !p[..4].iter().flat_map(|p| p.to_array()).all(f64::is_finite) {
            assert!(matches!(o, Err(Error::NonFinite(_))));
        }
        return;
    }
    let r = p.map(rational);
    let d = orientation(&r[..4]);
    assert_eq!(
        o.unwrap(),
        match d.cmp(&zero()) {
            Ordering::Less => Orientation3::Negative,
            Ordering::Equal => Orientation3::Coplanar,
            Ordering::Greater => Orientation3::Positive,
        }
    );
    if d == zero() {
        assert!(matches!(s, Err(Error::Degenerate(_))));
    } else {
        let sphere = determinant(
            r.iter()
                .map(|p| [vec![one()], p.to_vec(), vec![p.iter().map(|x| x * x).sum()]].concat())
                .collect(),
        );
        let expected = if sphere == zero() {
            SphereLocation::Boundary
        } else if sphere.cmp(&zero()) == d.cmp(&zero()) {
            SphereLocation::Outside
        } else {
            SphereLocation::Inside
        };
        assert_eq!(s.unwrap(), expected);
        assert_eq!(in_sphere(p[1], p[0], p[2], p[3], p[4]).unwrap(), expected);
    }
    let d2 = (&r[1][0] - &r[0][0]) * (&r[2][1] - &r[0][1])
        - (&r[1][1] - &r[0][1]) * (&r[2][0] - &r[0][0]);
    let expected = match d2.cmp(&zero()) {
        Ordering::Less => Orientation2::Clockwise,
        Ordering::Equal => Orientation2::Collinear,
        Ordering::Greater => Orientation2::CounterClockwise,
    };
    assert_eq!(
        orient2d(
            Point2::new(p[0].x, p[0].y),
            Point2::new(p[1].x, p[1].y),
            Point2::new(p[2].x, p[2].y)
        )
        .unwrap(),
        expected
    );
}

fn sub(a: &[R; 3], b: &[R; 3]) -> [R; 3] {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn cross(a: &[R; 3], b: &[R; 3]) -> [R; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

fn encloses(low: f64, high: f64, value: &R) {
    assert!(low.is_finite() && high.is_finite());
    assert!(R::from_float(low).unwrap() <= *value && *value <= R::from_float(high).unwrap());
    assert!(
        low == high || low.next_up() == high,
        "enclosure must be minimal"
    );
}

fn check_point(actual: IntersectionPoint, p: &[R; 3], q: &[R; 3], t: &R) {
    encloses(actual.parameter().lower(), actual.parameter().upper(), t);
    let bounds = actual.bounds();
    for i in 0..3 {
        encloses(
            bounds.min.to_array()[i],
            bounds.max.to_array()[i],
            &(&p[i] + t * (&q[i] - &p[i])),
        );
        assert!(
            bounds.min.to_array()[i] <= actual.position().to_array()[i]
                && actual.position().to_array()[i] <= bounds.max.to_array()[i]
        );
    }
}

fn check_result(
    result: rusty_occt::Result<(u8, Vec<IntersectionPoint>)>,
    status: u8,
    parameters: Vec<R>,
    p: &[R; 3],
    q: &[R; 3],
) {
    let maximum = R::from_float(f64::MAX).unwrap();
    let outside = |r: R| r > maximum || r < -&maximum;
    let unrepresentable = parameters
        .iter()
        .any(|t| outside(t.clone()) || (0..3).any(|i| outside(&p[i] + t * (&q[i] - &p[i]))));
    if unrepresentable {
        assert!(matches!(result, Err(Error::Unrepresentable(_))));
        return;
    }
    let (actual, points) = result.unwrap();
    assert_eq!(actual, status);
    assert_eq!(points.len(), parameters.len());
    for (point, t) in points.into_iter().zip(parameters) {
        check_point(point, p, q, &t);
    }
}

pub fn check_intersections(data: &[u8]) {
    let p = decode(data);
    let triangle = Triangle3::new(p[0], p[1], p[2]);
    if !p[..3].iter().flat_map(|p| p.to_array()).all(f64::is_finite) {
        assert!(matches!(triangle, Err(Error::NonFinite(_))));
        return;
    }
    let a = rational(p[0]);
    let b = rational(p[1]);
    let c = rational(p[2]);
    let u = sub(&b, &a);
    let v = sub(&c, &a);
    let normal = cross(&u, &v);
    if normal.iter().all(|x| *x == zero()) {
        assert!(matches!(triangle, Err(Error::Degenerate(_))));
        return;
    }
    let triangle = triangle.unwrap();
    let line = line_plane(p[3], p[4], triangle.plane()).map(|r| match r {
        L::Disjoint => (b'N', vec![]),
        L::Contained => (b'C', vec![]),
        L::Point(x) => (b'P', vec![x]),
    });
    let segment = segment_plane(p[3], p[4], triangle.plane()).map(|r| match r {
        P::Disjoint => (b'N', vec![]),
        P::Contained => (b'C', vec![]),
        P::Point(x) => (b'P', vec![x]),
    });
    let tri = segment_triangle(p[3], p[4], &triangle).map(|r| match r {
        T::Disjoint => (b'N', vec![]),
        T::Point(x) => (b'P', vec![x]),
        T::Overlap { start, end } => (b'S', vec![start, end]),
    });
    if !p[3..].iter().flat_map(|p| p.to_array()).all(f64::is_finite) {
        for result in [line, segment, tri] {
            assert!(matches!(result, Err(Error::NonFinite(_))));
        }
        return;
    }
    let (rp, rq) = (rational(p[3]), rational(p[4]));
    let distance = |p: &[R; 3]| -> R { normal.iter().zip(sub(p, &a)).map(|(n, x)| n * x).sum() };
    let (dp, dq) = (distance(&rp), distance(&rq));
    if rp == rq {
        assert!(matches!(line, Err(Error::Degenerate(_))));
    } else if dp == dq {
        check_result(
            line,
            if dp == zero() { b'C' } else { b'N' },
            vec![],
            &rp,
            &rq,
        );
    } else {
        check_result(line, b'P', vec![&dp / (&dp - &dq)], &rp, &rq);
    }
    if dp == zero() && dq == zero() {
        if rp == rq {
            check_result(segment, b'P', vec![zero()], &rp, &rq);
        } else {
            check_result(segment, b'C', vec![], &rp, &rq);
        }
    } else if dp.cmp(&zero()) == dq.cmp(&zero()) {
        check_result(segment, b'N', vec![], &rp, &rq);
    } else {
        check_result(segment, b'P', vec![&dp / (&dp - &dq)], &rp, &rq);
    }

    let drop = normal.iter().position(|x| *x != zero()).unwrap();
    let axes: Vec<_> = (0..3).filter(|&i| i != drop).collect();
    let (i, j) = (axes[0], axes[1]);
    let det = &u[i] * &v[j] - &u[j] * &v[i];
    let barycentric = |p: &[R; 3]| {
        let w = sub(p, &a);
        let beta = (&w[i] * &v[j] - &w[j] * &v[i]) / &det;
        let gamma = (&u[i] * &w[j] - &u[j] * &w[i]) / &det;
        [one() - &beta - &gamma, beta, gamma]
    };
    let (bp, bq) = (barycentric(&rp), barycentric(&rq));
    let (status, params) = if dp != zero() || dq != zero() {
        if dp.cmp(&zero()) == dq.cmp(&zero()) {
            (b'N', vec![])
        } else {
            let t = &dp / (&dp - &dq);
            if (0..3).any(|i| &bp[i] + &t * (&bq[i] - &bp[i]) < zero()) {
                (b'N', vec![])
            } else {
                (b'P', vec![t])
            }
        }
    } else {
        let (mut low, mut high) = (zero(), one());
        for i in 0..3 {
            if bp[i] < zero() && bq[i] < zero() {
                low = one();
                high = zero();
                break;
            }
            match bq[i].cmp(&bp[i]) {
                Ordering::Greater => low = low.max(-&bp[i] / (&bq[i] - &bp[i])),
                Ordering::Less => high = high.min(-&bp[i] / (&bq[i] - &bp[i])),
                Ordering::Equal => (),
            }
        }
        if low > high {
            (b'N', vec![])
        } else if rp == rq {
            (b'P', vec![zero()])
        } else if low == high {
            (b'P', vec![low])
        } else {
            (b'S', vec![low, high])
        }
    };
    check_result(tri, status, params, &rp, &rq);
}
