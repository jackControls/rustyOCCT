//! Exact decisions about the supplied finite binary64 coordinates.
//!
//! Predicates do not apply modeling tolerance and do not construct new geometry.
//! A point can be strictly on one side of a line and also inside its tolerance
//! band. Keep those questions separate. See `rust/MATHEMATICS.md`.
use crate::math::finite;
use crate::{Point2, Result};
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

#[cfg(test)]
mod tests {
    use super::*;

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
