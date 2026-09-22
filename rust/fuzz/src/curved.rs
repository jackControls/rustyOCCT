//! Rational polynomial-sign and axial-projection oracle. It does not use a
//! radical expression, a floating square root, or production root coefficients.
use super::{byte, one, zero};
use num_rational::BigRational as R;
use rusty_occt::intersection::*;
use rusty_occt::polynomial::{quadratic_roots, QuadraticRoot, QuadraticRoots};
use rusty_occt::{Error, Point3, Result, ScalarInterval, Vec3};
use std::cmp::Ordering;

#[derive(Clone)]
pub(crate) enum RefRoot {
    Rational(R, u8),
    Quadratic { a: R, b: R, c: R, upper: bool },
}
impl RefRoot {
    pub(crate) fn compare(&self, x: &R) -> Ordering {
        match self {
            Self::Rational(value, _) => value.cmp(x),
            Self::Quadratic { a, b, c, upper } => {
                let vertex = -b / (a + a);
                let polynomial = a * x * x + b * x + c;
                if *upper {
                    if *x < vertex {
                        Ordering::Greater
                    } else {
                        polynomial.cmp(&zero()).reverse()
                    }
                } else if *x > vertex {
                    Ordering::Less
                } else {
                    polynomial.cmp(&zero())
                }
            }
        }
    }
    pub(crate) fn multiplicity(&self) -> u8 {
        match self {
            Self::Rational(_, n) => *n,
            _ => 1,
        }
    }
}
fn roots(mut a: R, mut b: R, mut c: R) -> Option<Vec<RefRoot>> {
    if a == zero() {
        return if b != zero() {
            Some(vec![RefRoot::Rational(-c / b, 1)])
        } else if c == zero() {
            None
        } else {
            Some(vec![])
        };
    }
    if a < zero() {
        a = -a;
        b = -b;
        c = -c;
    }
    let d = &b * &b - (&a + &a) * (&c + &c);
    if d < zero() {
        return Some(vec![]);
    }
    if d == zero() {
        return Some(vec![RefRoot::Rational(-b / (&a + &a), 2)]);
    }
    Some(vec![
        RefRoot::Quadratic {
            a: a.clone(),
            b: b.clone(),
            c: c.clone(),
            upper: false,
        },
        RefRoot::Quadratic {
            a,
            b,
            c,
            upper: true,
        },
    ])
}

fn unrepresentable(compare: &impl Fn(&R) -> Ordering) -> bool {
    compare(&R::from_float(-f64::MAX).unwrap()) == Ordering::Less
        || compare(&R::from_float(f64::MAX).unwrap()) == Ordering::Greater
}
pub(crate) fn interval(bounds: ScalarInterval, compare: impl Fn(&R) -> Ordering) {
    let (lo, hi) = (bounds.lower(), bounds.upper());
    assert!(lo.is_finite() && hi.is_finite() && lo <= hi);
    let (low, high) = (
        compare(&R::from_float(lo).unwrap()),
        compare(&R::from_float(hi).unwrap()),
    );
    if lo == hi {
        assert_eq!(low, Ordering::Equal);
        assert_eq!(high, Ordering::Equal);
    } else {
        assert_eq!(lo.next_up(), hi);
        assert_eq!(low, Ordering::Greater);
        assert_eq!(high, Ordering::Less);
    }
    assert!(lo <= bounds.representative() && bounds.representative() <= hi);
}
fn root_check(root: QuadraticRoot, expected: &RefRoot, probe: f64) {
    assert_eq!(root.multiplicity(), expected.multiplicity());
    let compare = |x: &R| expected.compare(x);
    if unrepresentable(&compare) {
        assert!(matches!(root.bounds(), Err(Error::Unrepresentable(_))));
    } else {
        interval(root.bounds().unwrap(), compare);
    }
    if let Some(probe_r) = R::from_float(probe) {
        assert_eq!(root.compare(probe).unwrap(), expected.compare(&probe_r));
    } else {
        assert!(matches!(root.compare(probe), Err(Error::NonFinite(_))));
    }
}
fn polynomial(v: &[f64; 13]) {
    let actual = quadratic_roots(v[0], v[1], v[2]);
    if !v[..3].iter().all(|x| x.is_finite()) {
        assert!(matches!(actual, Err(Error::NonFinite(_))));
        return;
    }
    let r = v[..3]
        .iter()
        .map(|&x| R::from_float(x).unwrap())
        .collect::<Vec<_>>();
    let expected = roots(r[0].clone(), r[1].clone(), r[2].clone());
    let actual = match actual.unwrap() {
        QuadraticRoots::All => {
            assert!(expected.is_none());
            return;
        }
        QuadraticRoots::None => vec![],
        QuadraticRoots::One(root) => vec![root],
        QuadraticRoots::Two { lower, upper } => vec![lower, upper],
    };
    let expected = expected.unwrap();
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.into_iter().zip(expected) {
        root_check(actual, &expected, v[3]);
    }
}

fn dot(a: &[R], b: &[R]) -> R {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn sub(a: &[R], b: &[R]) -> Vec<R> {
    a.iter().zip(b).map(|(a, b)| a - b).collect()
}
fn evaluate(a: &R, b: &R, c: &R, t: &R) -> R {
    a * t * t + b * t + c
}

pub(crate) fn expected_surface(
    kind: u8,
    v: &[R],
    segment: bool,
) -> Option<Vec<(RefRoot, ContactKind)>> {
    let (center, axis, radius, p, q) = (&v[..3], &v[3..6], &v[6], &v[7..10], &v[10..13]);
    let (w, d) = (sub(p, center), sub(q, p));
    let (a, b, c) = if kind == 2 {
        let length2 = dot(axis, axis);
        let wd = dot(&w, axis);
        let dd = dot(&d, axis);
        let a = dot(&d, &d) - &dd * &dd / &length2;
        let b = dot(&w, &d) - &wd * &dd / &length2;
        let c = dot(&w, &w) - &wd * &wd / &length2 - radius * radius;
        (a, &b + &b, c)
    } else {
        let b = dot(&w, &d);
        (dot(&d, &d), &b + &b, dot(&w, &w) - radius * radius)
    };
    let candidates = if kind == 3 {
        let (initial, direction) = (dot(&w, axis), dot(&d, axis));
        if direction == zero() {
            if initial != zero() {
                return Some(vec![]);
            }
            roots(a, b, c)
        } else {
            let t = -initial / direction;
            Some(if evaluate(&a, &b, &c, &t) == zero() {
                vec![RefRoot::Rational(t, 1)]
            } else {
                vec![]
            })
        }
    } else {
        roots(a, b, c)
    };
    match candidates {
        None if p == q => Some(vec![(
            RefRoot::Rational(zero(), 1),
            ContactKind::DegeneratePoint,
        )]),
        None => None,
        Some(roots) => Some(
            roots
                .into_iter()
                .filter(|root| {
                    !segment
                        || (root.compare(&zero()) != Ordering::Less
                            && root.compare(&one()) != Ordering::Greater)
                })
                .map(|root| {
                    let kind = if root.multiplicity() == 2 {
                        ContactKind::Tangent
                    } else {
                        ContactKind::Crossing
                    };
                    (root, kind)
                })
                .collect(),
        ),
    }
}

pub(crate) fn affine_compare(root: &RefRoot, offset: &R, direction: &R, x: &R) -> Ordering {
    if *direction == zero() {
        return offset.cmp(x);
    }
    let order = root.compare(&((x - offset) / direction));
    if *direction < zero() {
        order.reverse()
    } else {
        order
    }
}
fn surface_result(
    actual: Result<CurvedIntersection>,
    expected: Option<Vec<(RefRoot, ContactKind)>>,
    v: &[R],
) {
    let Some(expected) = expected else {
        assert!(matches!(actual, Ok(CurvedIntersection::Contained)));
        return;
    };
    let (p, q) = (&v[7..10], &v[10..13]);
    let d = sub(q, p);
    let outside = expected.iter().any(|(root, _)| {
        unrepresentable(&|x| root.compare(x))
            || (0..3).any(|i| unrepresentable(&|x| affine_compare(root, &p[i], &d[i], x)))
    });
    if outside {
        assert!(matches!(actual, Err(Error::Unrepresentable(_))));
        return;
    }
    let hits = match actual.unwrap() {
        CurvedIntersection::Disjoint => vec![],
        CurvedIntersection::One(hit) => vec![hit],
        CurvedIntersection::Two { first, second } => vec![first, second],
        CurvedIntersection::Contained => panic!("unexpected containment"),
    };
    assert_eq!(hits.len(), expected.len());
    for (hit, (root, kind)) in hits.into_iter().zip(expected) {
        assert_eq!(hit.kind, kind);
        interval(hit.point.parameter(), |x| root.compare(x));
        let bounds = hit.point.bounds();
        for i in 0..3 {
            let (lo, hi) = (bounds.min.to_array()[i], bounds.max.to_array()[i]);
            let compare = |x: &R| affine_compare(&root, &p[i], &d[i], x);
            assert!(lo.is_finite() && hi.is_finite());
            if lo == hi {
                assert_eq!(compare(&R::from_float(lo).unwrap()), Ordering::Equal);
            } else {
                assert_eq!(lo.next_up(), hi);
                assert_eq!(compare(&R::from_float(lo).unwrap()), Ordering::Greater);
                assert_eq!(compare(&R::from_float(hi).unwrap()), Ordering::Less);
            }
            assert!(
                lo <= hit.point.position().to_array()[i]
                    && hit.point.position().to_array()[i] <= hi
            );
        }
    }
}

pub fn check_curved(data: &[u8]) {
    let mode = byte(data, 0) % 4;
    let kind = byte(data, 1) % 4;
    let mut v: [f64; 13] = std::array::from_fn(|i| {
        let word = u64::from_le_bytes(std::array::from_fn(|j| byte(data, 2 + 8 * i + j)));
        if mode == 0 {
            f64::from_bits(word)
        } else {
            (word as i16 as f64) * 2.0_f64.powi(byte(data, 106) as i32 - 128)
        }
    });
    if kind == 0 {
        if mode == 2 {
            v[0] = 1.;
            v[1] = -2.;
            v[2] = f64::from_bits(1.0_f64.to_bits() - 1 + (byte(data, 107) % 3) as u64);
        }
        if mode == 3 {
            v[0] = f64::from_bits(1);
            v[1] = 1.;
            v[2] = byte(data, 107) as f64;
        }
        polynomial(&v);
        return;
    }
    if mode != 0 {
        v[6] = v[6].abs() + 1.;
        if v[3..6].iter().all(|&x| x == 0.) {
            v[5] = 1.;
        }
    }
    if mode == 2 || mode == 3 {
        let radius = 2.0_f64.powi(byte(data, 106) as i32 - 128);
        v = [
            0.,
            0.,
            0.,
            0.,
            0.,
            1.,
            radius,
            -2. * radius,
            radius,
            0.,
            2. * radius,
            radius,
            0.,
        ];
        if mode == 2 {
            let y = f64::from_bits(radius.to_bits() - 1 + (byte(data, 107) % 3) as u64);
            v[8] = y;
            v[11] = y;
        } else {
            v[7] = radius;
            v[8] = 0.;
            v[9] = -2. * radius;
            v[10] = radius;
            v[11] = 0.;
            v[12] = 2. * radius;
            if byte(data, 107).is_multiple_of(2) {
                v[9] = 0.;
                v[12] = 0.;
            }
        }
    }
    let center = Point3::new(v[0], v[1], v[2]);
    let axis = Vec3::new(v[3], v[4], v[5]);
    let p = Point3::new(v[7], v[8], v[9]);
    let q = Point3::new(v[10], v[11], v[12]);
    let valid = v[..3].iter().all(|x| x.is_finite())
        && v[6].is_finite()
        && v[6] > 0.
        && (kind == 1
            || (v[3..6].iter().all(|x| x.is_finite()) && v[3..6].iter().any(|&x| x != 0.)));
    let pair = match kind {
        1 => Sphere3::new(center, v[6]).map(|s| (line_sphere(p, q, &s), segment_sphere(p, q, &s))),
        2 => Cylinder3::new(center, axis, v[6])
            .map(|s| (line_cylinder(p, q, &s), segment_cylinder(p, q, &s))),
        3 => Circle3::new(center, axis, v[6])
            .map(|s| (line_circle(p, q, &s), segment_circle(p, q, &s))),
        _ => unreachable!(),
    };
    if !valid {
        assert!(pair.is_err());
        return;
    }
    let (line, segment) = pair.unwrap();
    if !v[7..].iter().all(|x| x.is_finite()) {
        assert!(matches!(line, Err(Error::NonFinite(_))));
        assert!(matches!(segment, Err(Error::NonFinite(_))));
        return;
    }
    // A sphere does not consume the axis fields, including NaNs in those bytes.
    if kind == 1 {
        v[3..6].fill(0.);
    }
    let rational = v.map(|x| R::from_float(x).unwrap());
    if p == q {
        assert!(matches!(line, Err(Error::Degenerate(_))));
    } else {
        surface_result(line, expected_surface(kind, &rational, false), &rational);
    }
    surface_result(segment, expected_surface(kind, &rational, true), &rational);
}
