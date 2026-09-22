//! Finite binary64 values expressed as integers in units of 2^-1074.
//! Input size and polynomial degree bound every allocation; no user precision
//! parameter or convergence loop enters the exact geometry algorithms.
use crate::{math::finite, Point3, Result};
use num_bigint::{BigInt, Sign};

pub(crate) type Vector = [BigInt; 3];

pub(crate) fn integer(value: f64) -> BigInt {
    debug_assert!(value.is_finite());
    let bits = value.to_bits();
    let exponent = (bits >> 52) & 0x7ff;
    let fraction = bits & ((1_u64 << 52) - 1);
    let (mantissa, shift) = if exponent == 0 {
        (fraction, 0)
    } else {
        (fraction | (1_u64 << 52), exponent - 1)
    };
    let magnitude = BigInt::from(mantissa) << shift as usize;
    if bits >> 63 == 0 {
        magnitude
    } else {
        -magnitude
    }
}

pub(crate) fn point(p: Point3) -> Result<Vector> {
    for value in p.to_array() {
        finite(value, "exact geometry coordinate")?;
    }
    Ok(p.to_array().map(integer))
}

pub(crate) fn sub(a: &Vector, b: &Vector) -> Vector {
    std::array::from_fn(|i| &a[i] - &b[i])
}

pub(crate) fn cross(a: &Vector, b: &Vector) -> Vector {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

pub(crate) fn dot(a: &Vector, b: &Vector) -> BigInt {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

pub(crate) fn determinant(a: &Vector, b: &Vector, c: &Vector) -> BigInt {
    dot(&cross(a, b), c)
}

pub(crate) fn orientation(a: &Vector, b: &Vector, c: &Vector, d: &Vector) -> BigInt {
    determinant(&sub(b, a), &sub(c, a), &sub(d, a))
}

pub(crate) fn zero(value: &BigInt) -> bool {
    value.sign() == Sign::NoSign
}
