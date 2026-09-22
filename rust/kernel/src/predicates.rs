//! Exact decisions about the supplied finite binary64 coordinates.
//!
//! Predicates do not apply modeling tolerance and do not construct new geometry.
//! A point can be strictly on one side of a line and also inside its tolerance
//! band. Keep those questions separate. See `rust/MATHEMATICS.md`.
use crate::math::finite;
use crate::{exact as integer, Error, Point2, Point3, Result};
use num_bigint::Sign;
use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation2 {
    Clockwise,
    Collinear,
    CounterClockwise,
}

impl Orientation2 {
    fn from_ordering(ordering: Ordering) -> Self {
        match ordering {
            Ordering::Less => Self::Clockwise,
            Ordering::Equal => Self::Collinear,
            Ordering::Greater => Self::CounterClockwise,
        }
    }
}

/// The exact sign of `(b - a) x (c - a)` for the represented inputs.
///
/// Handles every finite `f64`, including subnormals and coordinates for which
/// floating-point differences/products would overflow. NaN/infinity are errors.
/// This gives no guarantee about coordinates rounded *before* this call, or
/// about the accuracy of a subsequently constructed intersection or distance.
pub fn orient2d(a: Point2, b: Point2, c: Point2) -> Result<Orientation2> {
    for value in [a.x, a.y, b.x, b.y, c.x, c.y] {
        finite(value, "orientation coordinate")?;
    }
    Ok(orient2d_finite(a, b, c))
}

pub(crate) fn orient2d_finite(a: Point2, b: Point2, c: Point2) -> Orientation2 {
    let coordinates = [a.x, a.y, b.x, b.y, c.x, c.y];
    debug_assert!(coordinates.iter().all(|v| v.is_finite()));
    if let Some(result) = filtered(a, b, c) {
        return result;
    }
    exact(a, b, c)
}

fn filtered(a: Point2, b: Point2, c: Point2) -> Option<Orientation2> {
    // Shewchuk, predicates.c (public domain, 1996): orient2d / ccwerrboundA.
    // The bound assumes no floating-point overflow/underflow. Conservatively
    // restrict nonzero coordinates to [2^-400, 2^400]; their nonzero differences
    // are at least 2^-452, keeping products AND the error bound normal.
    const SMALL: f64 = f64::from_bits((1023 - 400) << 52);
    const LARGE: f64 = f64::from_bits((1023 + 400) << 52);
    if ![a.x, a.y, b.x, b.y, c.x, c.y]
        .iter()
        .all(|v| *v == 0.0 || (SMALL..=LARGE).contains(&v.abs()))
    {
        return None;
    }
    let left = (a.x - c.x) * (b.y - c.y);
    let right = (a.y - c.y) * (b.x - c.x);
    let determinant = left - right;
    const UNIT_ROUNDOFF: f64 = f64::EPSILON * 0.5;
    const ERROR_BOUND: f64 = (3.0 + 16.0 * UNIT_ROUNDOFF) * UNIT_ROUNDOFF;
    let bound = ERROR_BOUND * (left.abs() + right.abs());
    if determinant.abs() > bound {
        Some(if determinant > 0.0 {
            Orientation2::CounterClockwise
        } else {
            Orientation2::Clockwise
        })
    } else {
        None
    }
}

// Every finite f64 is sign * mantissa * 2^exponent, where mantissa has at most
// 53 bits and exponent is in [-1074, 971]. Express each exact product in units
// of 2^-2148. Its shift is in [0, 4090], and its integer occupies <=4196 bits.
// Even six products need <4199 bits. 66 limbs (4224 bits) suffice, with no heap,
// rounded subtraction, floating-point multiplication, or variable-size input.
const LIMBS: usize = 66;

fn parts(value: f64) -> (bool, u64, i32) {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1_u64 << 52) - 1);
    let (mantissa, power) = if exponent == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1_u64 << 52), exponent - 1075)
    };
    (bits >> 63 != 0, mantissa, power)
}

fn add_word(sum: &mut [u64; LIMBS], mut index: usize, mut word: u64) {
    while word != 0 {
        let (value, carry) = sum[index].overflowing_add(word);
        sum[index] = value;
        word = u64::from(carry);
        index += 1;
    }
}

fn add_product(sum: &mut [u64; LIMBS], product: u128, shift: usize) {
    let index = shift / 64;
    let offset = shift % 64;
    let low = product as u64;
    let high = (product >> 64) as u64;
    add_word(sum, index, low << offset);
    if offset == 0 {
        add_word(sum, index + 1, high);
    } else {
        add_word(sum, index + 1, (low >> (64 - offset)) | (high << offset));
        add_word(sum, index + 2, high >> (64 - offset));
    }
}

fn exact(a: Point2, b: Point2, c: Point2) -> Orientation2 {
    let mut positive = [0_u64; LIMBS];
    let mut negative = [0_u64; LIMBS];
    // Expanded determinant cancels the a.x*a.y terms algebraically, avoiding
    // any rounded coordinate differences: ax*by + bx*cy + cx*ay - ay*bx - by*cx - cy*ax.
    for (x, y, subtract) in [
        (a.x, b.y, false),
        (b.x, c.y, false),
        (c.x, a.y, false),
        (a.y, b.x, true),
        (b.y, c.x, true),
        (c.y, a.x, true),
    ] {
        let (sx, mx, ex) = parts(x);
        let (sy, my, ey) = parts(y);
        let product = u128::from(mx) * u128::from(my);
        let sum = if sx ^ sy ^ subtract {
            &mut negative
        } else {
            &mut positive
        };
        add_product(sum, product, (ex + ey + 2148) as usize);
    }
    Orientation2::from_ordering(positive.iter().rev().cmp(negative.iter().rev()))
}

/// Sign of `((b - a) x (c - a)) dot (d - a)`.
/// Positive is above the oriented plane (opposite Shewchuk's orient3d sign).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation3 {
    Negative,
    Coplanar,
    Positive,
}

/// Exact sidedness for all finite binary64 coordinates, without a tolerance.
/// Degenerate defining planes return `Coplanar`; NaN/infinity return an error.
pub fn orient3d(a: Point3, b: Point3, c: Point3, d: Point3) -> Result<Orientation3> {
    for p in [a, b, c, d] {
        for x in p.to_array() {
            finite(x, "orientation coordinate")?;
        }
    }
    if let Some(sign) = filtered3d(a, b, c, d) {
        return Ok(sign);
    }
    let value = integer::orientation(
        &integer::point(a)?,
        &integer::point(b)?,
        &integer::point(c)?,
        &integer::point(d)?,
    );
    Ok(match value.sign() {
        Sign::Minus => Orientation3::Negative,
        Sign::NoSign => Orientation3::Coplanar,
        Sign::Plus => Orientation3::Positive,
    })
}

fn filtered3d(a: Point3, b: Point3, c: Point3, d: Point3) -> Option<Orientation3> {
    // Shewchuk predicates.c: orient3d / o3derrboundA (public domain, 1996).
    // Nonzero differences >=2^-252, triple products >=2^-756. Products,
    // subtraction residuals, and the error bound stay normal on this domain.
    const SMALL: f64 = f64::from_bits((1023 - 200) << 52);
    const LARGE: f64 = f64::from_bits((1023 + 200) << 52);
    if ![a, b, c, d]
        .iter()
        .flat_map(|p| p.to_array())
        .all(|x| x == 0.0 || (SMALL..=LARGE).contains(&x.abs()))
    {
        return None;
    }
    let [ax, ay, az] = (a - d).to_array();
    let [bx, by, bz] = (b - d).to_array();
    let [cx, cy, cz] = (c - d).to_array();
    let [bc, cb, ca, ac, ab, ba] = [bx * cy, cx * by, cx * ay, ax * cy, ax * by, bx * ay];
    let determinant = az * (bc - cb) + bz * (ca - ac) + cz * (ab - ba);
    let permanent = (bc.abs() + cb.abs()) * az.abs()
        + (ca.abs() + ac.abs()) * bz.abs()
        + (ab.abs() + ba.abs()) * cz.abs();
    const U: f64 = f64::EPSILON * 0.5;
    const BOUND: f64 = (7.0 + 56.0 * U) * U;
    if determinant.abs() > BOUND * permanent {
        Some(if determinant < 0.0 {
            Orientation3::Positive
        } else {
            Orientation3::Negative
        })
    } else {
        None
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SphereLocation {
    Inside,
    Boundary,
    Outside,
}

/// Exact position relative to the sphere through four finite points.
/// Independent of defining-point order. A coplanar defining tetrahedron has
/// no unique sphere and returns `Error::Degenerate`. No center is constructed.
pub fn in_sphere(
    a: Point3,
    b: Point3,
    c: Point3,
    d: Point3,
    query: Point3,
) -> Result<SphereLocation> {
    let [a, b, c, d, q] = [
        integer::point(a)?,
        integer::point(b)?,
        integer::point(c)?,
        integer::point(d)?,
        integer::point(query)?,
    ];
    let orientation = integer::orientation(&a, &b, &c, &d);
    if integer::zero(&orientation) {
        return Err(Error::Degenerate("sphere defining tetrahedron"));
    }
    let [a, b, c, d] = [&a, &b, &c, &d].map(|p| integer::sub(p, &q));
    // Expand the [dx, dy, dz, squared distance] determinant along its last
    // column. Its interior sign is opposite our above-plane convention.
    let determinant = -integer::dot(&a, &a) * integer::determinant(&b, &c, &d)
        + integer::dot(&b, &b) * integer::determinant(&a, &c, &d)
        - integer::dot(&c, &c) * integer::determinant(&a, &b, &d)
        + integer::dot(&d, &d) * integer::determinant(&a, &b, &c);
    Ok(if integer::zero(&determinant) {
        SphereLocation::Boundary
    } else if determinant.sign() == orientation.sign() {
        SphereLocation::Outside
    } else {
        SphereLocation::Inside
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orientation3_filter_and_fallback_domains() {
        let a = Point3::ORIGIN;
        let b = Point3::new(1., 0., 0.);
        let c = Point3::new(0., 1., 0.);
        assert_eq!(
            filtered3d(a, b, c, Point3::new(0., 0., 1.)),
            Some(Orientation3::Positive)
        );
        assert_eq!(filtered3d(a, b, c, Point3::new(0.5, 0.5, 0.)), None);
        assert_eq!(
            filtered3d(a, b, c, Point3::new(0., 0., f64::from_bits(1))),
            None
        );
        assert_eq!(filtered3d(a, b, c, Point3::new(0., 0., f64::MAX)), None);
    }

    #[test]
    fn filter_certifies_only_well_separated_normal_arithmetic() {
        assert_eq!(
            filtered(
                Point2::default(),
                Point2::new(1.0, 0.0),
                Point2::new(0.0, 1.0)
            ),
            Some(Orientation2::CounterClockwise)
        );
        let u = 134_217_728.0;
        assert_eq!(
            filtered(
                Point2::new(u, u - 1.0),
                Point2::new(u + 1.0, u),
                Point2::default()
            ),
            None
        );
        assert_eq!(
            filtered(
                Point2::default(),
                Point2::new(f64::MAX, 0.0),
                Point2::new(0.0, f64::MAX)
            ),
            None
        );
        assert_eq!(
            filtered(
                Point2::default(),
                Point2::new(f64::from_bits(1), 0.0),
                Point2::new(0.0, f64::from_bits(1))
            ),
            None
        );
    }

    #[test]
    fn exact_accumulator_can_carry_six_maximum_products() {
        let mantissa = (1_u128 << 53) - 1;
        let mut repeated = [0; LIMBS];
        for _ in 0..6 {
            add_product(&mut repeated, mantissa * mantissa, 4090);
        }
        let mut combined = [0; LIMBS];
        add_product(&mut combined, 6 * mantissa * mantissa, 4090);
        assert_eq!(repeated, combined);
        assert_ne!(repeated[LIMBS - 1], 0);
    }
}
