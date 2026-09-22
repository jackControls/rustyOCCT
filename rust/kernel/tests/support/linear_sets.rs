//! Independent exact boundary oracle: closed-form line intersections, triangle
//! edges, and a 3D gift-wrapping hull. No production affine elimination/halfspace
//! vertex enumeration. Used by deterministic tests and instrumented fuzzing.
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::intersection::{ExactPoint3, LinearIntersection as I, LinearPrimitive3 as P};
use rusty_occt::{Error, Point3};
type V = [R; 3];
fn z() -> R {
    R::from_integer(BigInt::from(0))
}
fn one() -> R {
    R::from_integer(BigInt::from(1))
}
fn sub(a: &V, b: &V) -> V {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn add(a: &V, b: &V) -> V {
    std::array::from_fn(|i| &a[i] + &b[i])
}
fn mul(a: &V, t: &R) -> V {
    std::array::from_fn(|i| &a[i] * t)
}
fn dot(a: &V, b: &V) -> R {
    (0..3).map(|i| &a[i] * &b[i]).sum()
}
fn cross(a: &V, b: &V) -> V {
    std::array::from_fn(|i| &a[(i + 1) % 3] * &b[(i + 2) % 3] - &a[(i + 2) % 3] * &b[(i + 1) % 3])
}
fn norm(a: &V) -> R {
    dot(a, a)
}
fn normal(p: &[V]) -> V {
    cross(&sub(&p[1], &p[0]), &sub(&p[2], &p[0]))
}
fn point(p: Point3) -> V {
    p.to_array().map(|x| R::from_float(x).unwrap())
}
fn exact(p: V) -> ExactPoint3 {
    ExactPoint3::from_coordinates(p)
}
#[derive(Clone)]
struct Shape {
    kind: u8,
    p: Vec<V>,
}
fn shape(p: &P) -> Shape {
    let (kind, p) = match p {
        P::Point(p) => (b'P', vec![*p]),
        P::Line(p) => (b'L', p.to_vec()),
        P::Segment(p) => (b'S', p.to_vec()),
        P::Plane(p) => (b'F', p.defining_points().to_vec()),
        P::Triangle(p) => (b'T', p.vertices().to_vec()),
    };
    let mut result = Shape {
        kind,
        p: p.into_iter().map(point).collect(),
    };
    if kind == b'S' && result.p[0] == result.p[1] {
        result.kind = b'P';
        result.p.truncate(1);
    }
    result
}
fn contains(a: &Shape, p: &V) -> bool {
    let w = sub(p, &a.p[0]);
    match a.kind {
        b'P' => p == &a.p[0],
        b'L' | b'S' => {
            let d = sub(&a.p[1], &a.p[0]);
            norm(&cross(&d, &w)) == z()
                && (a.kind == b'L' || (dot(&d, &w) >= z() && dot(&d, &w) <= norm(&d)))
        }
        b'F' | b'T' => {
            if dot(&normal(&a.p), &w) != z() {
                return false;
            }
            if a.kind == b'F' {
                return true;
            }
            // Independent Gram barycentric membership, not oriented halfspaces.
            let (d, e) = (sub(&a.p[1], &a.p[0]), sub(&a.p[2], &a.p[0]));
            let (dd, de, ee, dw, ew) = (norm(&d), dot(&d, &e), norm(&e), dot(&d, &w), dot(&e, &w));
            let det = &dd * &ee - &de * &de;
            let u = &dw * &ee - &ew * &de;
            let v = &ew * &dd - &dw * &de;
            u >= z() && v >= z() && u + v <= det
        }
        _ => unreachable!(),
    }
}
fn line(p: &V, d: &V) -> I {
    let axis = d.iter().position(|x| *x != z()).unwrap();
    let direction = mul(d, &(one() / &d[axis]));
    I::Line {
        origin: exact(sub(p, &mul(&direction, &p[axis]))),
        direction,
    }
}
fn plane(p: &[V]) -> I {
    let n = normal(p);
    let k = n.iter().find(|x| **x != z()).unwrap();
    I::Plane {
        normal: mul(&n, &(one() / k)),
        offset: dot(&n, &p[0]) / k,
    }
}
fn hull(mut p: Vec<V>) -> I {
    p.sort();
    p.dedup();
    if p.is_empty() {
        return I::Empty;
    }
    if p.len() == 1 {
        return I::Point(exact(p.remove(0)));
    }
    let d = sub(&p[1], &p[0]);
    let n = p[2..]
        .iter()
        .map(|q| cross(&d, &sub(q, &p[0])))
        .find(|n| norm(n) != z());
    let Some(n) = n else {
        return I::Segment([exact(p[0].clone()), exact(p.last().unwrap().clone())]);
    };
    let mut result = vec![p[0].clone()];
    loop {
        let a = result.last().unwrap();
        let mut b = p.iter().find(|p| *p != a).unwrap();
        for c in &p {
            let turn = dot(&n, &cross(&sub(b, a), &sub(c, a)));
            if turn < z() || (turn == z() && norm(&sub(c, a)) > norm(&sub(b, a))) {
                b = c;
            }
        }
        if b == &result[0] {
            break;
        }
        assert!(!result.contains(b), "reference hull must advance");
        result.push(b.clone());
    }
    if result[1] > result[result.len() - 1] {
        result[1..].reverse();
    }
    I::Polygon(result.into_iter().map(exact).collect())
}
fn line_line(a: &Shape, b: &Shape) -> I {
    let (u, v, w) = (
        sub(&a.p[1], &a.p[0]),
        sub(&b.p[1], &b.p[0]),
        sub(&b.p[0], &a.p[0]),
    );
    let n = cross(&u, &v);
    if norm(&n) != z() {
        if dot(&w, &n) != z() {
            return I::Empty;
        }
        let t = dot(&cross(&w, &v), &n) / norm(&n);
        let s = dot(&cross(&w, &u), &n) / norm(&n);
        if (a.kind == b'S' && (t < z() || t > one())) || (b.kind == b'S' && (s < z() || s > one()))
        {
            return I::Empty;
        }
        return I::Point(exact(add(&a.p[0], &mul(&u, &t))));
    }
    if norm(&cross(&w, &u)) != z() {
        return I::Empty;
    }
    if a.kind == b'L' && b.kind == b'L' {
        return line(&a.p[0], &u);
    }
    let mut points = Vec::new();
    for (x, y) in [(a, b), (b, a)] {
        if x.kind == b'S' {
            points.extend(x.p.iter().filter(|p| contains(y, p)).cloned());
        }
    }
    hull(points)
}
fn line_plane(a: &Shape, b: &Shape) -> I {
    let (d, n) = (sub(&a.p[1], &a.p[0]), normal(&b.p));
    let (den, num) = (dot(&d, &n), dot(&sub(&b.p[0], &a.p[0]), &n));
    if den == z() {
        return if num != z() {
            I::Empty
        } else if a.kind == b'S' {
            hull(a.p.clone())
        } else {
            line(&a.p[0], &d)
        };
    }
    let t = num / den;
    if a.kind == b'S' && (t < z() || t > one()) {
        I::Empty
    } else {
        I::Point(exact(add(&a.p[0], &mul(&d, &t))))
    }
}
fn edges(a: &Shape) -> Vec<Shape> {
    (0..3)
        .map(|i| Shape {
            kind: b'S',
            p: vec![a.p[i].clone(), a.p[(i + 1) % 3].clone()],
        })
        .collect()
}
fn append(result: I, points: &mut Vec<V>) {
    assert!(matches!(result, I::Empty | I::Point(_) | I::Segment(_)));
    points.extend(result.vertices().iter().map(|p| p.coordinates().clone()));
}
fn intersect(a: &Shape, b: &Shape) -> I {
    if a.kind == b'P' {
        return if contains(b, &a.p[0]) {
            I::Point(exact(a.p[0].clone()))
        } else {
            I::Empty
        };
    }
    if b.kind == b'P' {
        return intersect(b, a);
    }
    if a.kind == b'L' || a.kind == b'S' {
        if b.kind == b'L' || b.kind == b'S' {
            return line_line(a, b);
        }
        let hit = line_plane(a, b);
        if b.kind == b'F' || matches!(hit, I::Empty) {
            return hit;
        }
        if let I::Point(ref p) = hit {
            return if contains(b, p.coordinates()) {
                hit
            } else {
                I::Empty
            };
        }
        let mut points = if a.kind == b'S' {
            a.p.iter().filter(|p| contains(b, p)).cloned().collect()
        } else {
            Vec::new()
        };
        for edge in edges(b) {
            append(intersect(a, &edge), &mut points);
        }
        return hull(points);
    }
    if b.kind == b'L' || b.kind == b'S' {
        return intersect(b, a);
    }
    if a.kind == b'F' && b.kind == b'F' {
        let (u, v) = (normal(&a.p), normal(&b.p));
        let d = cross(&u, &v);
        if norm(&d) == z() {
            return if contains(a, &b.p[0]) {
                plane(&a.p)
            } else {
                I::Empty
            };
        }
        let p = mul(
            &add(
                &mul(&cross(&v, &d), &dot(&u, &a.p[0])),
                &mul(&cross(&d, &u), &dot(&v, &b.p[0])),
            ),
            &(one() / norm(&d)),
        );
        return line(&p, &d);
    }
    if a.kind == b'F' {
        return intersect(b, a);
    }
    let mut points = Vec::new();
    for edge in edges(a) {
        append(intersect(&edge, b), &mut points);
    }
    if b.kind == b'T' {
        for edge in edges(b) {
            append(intersect(&edge, a), &mut points);
        }
    }
    hull(points)
}

pub fn reference(a: &P, b: &P) -> I {
    intersect(&shape(a), &shape(b))
}

pub fn check_bounds(result: &I) {
    let mut points: Vec<_> = result.vertices().iter().collect();
    if let I::Line { origin, .. } = result {
        points.push(origin);
    }
    let max = R::from_float(f64::MAX).unwrap();
    for p in points {
        let mut representable = true;
        for i in 0..3 {
            let r = &p.coordinates()[i];
            let bounds = p.coordinate_bounds(i);
            if r > &max || r < &(-max.clone()) {
                representable = false;
                assert!(matches!(bounds, Err(Error::Unrepresentable(_))));
                continue;
            }
            let b = bounds.unwrap();
            let (lo, hi) = (b.lower(), b.upper());
            assert!(R::from_float(lo).unwrap() <= *r && *r <= R::from_float(hi).unwrap());
            if lo == hi {
                assert_eq!(*r, R::from_float(lo).unwrap());
            } else {
                let next = if lo == 0. {
                    f64::from_bits(1)
                } else if lo < 0. {
                    f64::from_bits(lo.to_bits() - 1)
                } else {
                    f64::from_bits(lo.to_bits() + 1)
                };
                assert_eq!(next, hi, "minimal enclosure");
                assert!(R::from_float(lo).unwrap() < *r && *r < R::from_float(hi).unwrap());
            }
        }
        if representable {
            let b = p.bounds().unwrap();
            let q = p.position().unwrap();
            for i in 0..3 {
                assert!(
                    b.min.to_array()[i] <= q.to_array()[i]
                        && q.to_array()[i] <= b.max.to_array()[i]
                );
            }
        } else {
            assert!(matches!(p.bounds(), Err(Error::Unrepresentable(_))));
            assert!(matches!(p.position(), Err(Error::Unrepresentable(_))));
        }
    }
}
