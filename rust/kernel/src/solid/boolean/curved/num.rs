//! S9c.1's numbers: rational vectors, quadratic surds `a + b sqrt(d)` over
//! the rationals and their exact signs, one surd or two apart, and their
//! enclosures.
use crate::certified::Interval as I;
use crate::solid::split::{rational_f64, zero};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::cmp::Ordering;

pub(super) type V = [R; 3];

pub(super) fn int(n: i64) -> R {
    R::from_integer(BigInt::from(n))
}

pub(super) fn add(a: &V, b: &V) -> V {
    [&a[0] + &b[0], &a[1] + &b[1], &a[2] + &b[2]]
}

pub(super) fn sub(a: &V, b: &V) -> V {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}

pub(super) fn scale(a: &V, k: &R) -> V {
    [&a[0] * k, &a[1] * k, &a[2] * k]
}

pub(super) fn dot(a: &V, b: &V) -> R {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

pub(super) fn cross(a: &V, b: &V) -> V {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

pub(super) fn neg(a: &V) -> V {
    [-&a[0], -&a[1], -&a[2]]
}

pub(super) fn is_zero(a: &V) -> bool {
    a.iter().all(|x| *x == zero())
}

pub(super) fn sign(x: &R) -> Ordering {
    x.cmp(&zero())
}

/// The exact square root of a rational square, if it is one.
pub(super) fn rational_sqrt(x: &R) -> Option<R> {
    if *x < zero() {
        return None;
    }
    let root = |n: &BigInt| {
        let r = n.sqrt();
        (&r * &r == *n).then_some(r)
    };
    Some(R::new(root(x.numer())?, root(x.denom())?))
}

/// `a + b sqrt(d)`, `d >= 0`; a rational has `b = d = 0`. A surd whose `d`
/// is a rational square is folded into its rational part.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Qd {
    pub(super) a: R,
    pub(super) b: R,
    pub(super) d: R,
}

impl Qd {
    pub(super) fn rat(a: R) -> Self {
        Self {
            a,
            b: zero(),
            d: zero(),
        }
    }

    pub(super) fn new(a: R, b: R, d: R) -> Self {
        debug_assert!(d >= zero());
        if b == zero() || d == zero() {
            return Self::rat(a);
        }
        if let Some(s) = rational_sqrt(&d) {
            return Self::rat(a + b * s);
        }
        Self { a, b, d }
    }

    pub(super) fn is_rational(&self) -> bool {
        self.b == zero()
    }

    /// The field's `d` of the surd part, when there is one.
    pub(super) fn field(&self) -> Option<&R> {
        (!self.is_rational()).then_some(&self.d)
    }

    /// The common field of two numbers, or `None` when both are surds of
    /// different `d`.
    fn common(&self, o: &Self) -> Option<R> {
        match (self.field(), o.field()) {
            (None, None) => Some(zero()),
            (Some(d), None) | (None, Some(d)) => Some(d.clone()),
            (Some(d), Some(e)) => (d == e).then(|| d.clone()),
        }
    }

    pub(super) fn add(&self, o: &Self) -> Self {
        let d = self.common(o).expect("surds of one field");
        Self::new(&self.a + &o.a, &self.b + &o.b, d)
    }

    pub(super) fn sub(&self, o: &Self) -> Self {
        let d = self.common(o).expect("surds of one field");
        Self::new(&self.a - &o.a, &self.b - &o.b, d)
    }

    pub(super) fn mul(&self, o: &Self) -> Self {
        let d = self.common(o).expect("surds of one field");
        Self::new(
            &self.a * &o.a + &self.b * &o.b * &d,
            &self.a * &o.b + &self.b * &o.a,
            d,
        )
    }

    pub(super) fn scale(&self, k: &R) -> Self {
        Self::new(&self.a * k, &self.b * k, self.d.clone())
    }

    pub(super) fn neg(&self) -> Self {
        Self::new(-&self.a, -&self.b, self.d.clone())
    }

    pub(super) fn add_r(&self, k: &R) -> Self {
        Self::new(&self.a + k, self.b.clone(), self.d.clone())
    }

    /// The exact sign.
    pub(super) fn sign(&self) -> Ordering {
        let sa = sign(&self.a);
        let sb = sign(&self.b);
        if sb == Ordering::Equal {
            return sa;
        }
        if sa == Ordering::Equal || sa == sb {
            return sb;
        }
        // Opposite signs: the larger magnitude wins.
        match (&self.a * &self.a).cmp(&(&self.b * &self.b * &self.d)) {
            Ordering::Greater => sa,
            Ordering::Less => sb,
            Ordering::Equal => Ordering::Equal,
        }
    }

    /// The exact order of two numbers of any fields.
    pub(super) fn cmp(&self, o: &Self) -> Ordering {
        if let Some(d) = self.common(o) {
            return Self::new(&self.a - &o.a, &self.b - &o.b, d).sign();
        }
        // (a - a' + b sqrt d) + (-b') sqrt d'.
        tower_sign(
            &Self::new(&self.a - &o.a, self.b.clone(), self.d.clone()),
            &Self::rat(-&o.b),
            &o.d,
        )
    }

    pub(super) fn interval(&self) -> I {
        let a = I::exact(self.a.clone());
        if self.is_rational() {
            return a;
        }
        a.add(&I::exact(self.b.clone()).mul(&I::exact(self.d.clone()).sqrt()))
    }

    /// The value rounded to binary64 (from a tight enclosure).
    pub(super) fn to_f64(&self) -> f64 {
        if self.is_rational() {
            return rational_f64(&self.a);
        }
        let i = self.interval();
        rational_f64(&((i.lo() + i.hi()) / int(2)))
    }
}

/// The exact sign of `x + y sqrt(e)`, `x` and `y` of one field, `e >= 0`.
pub(super) fn tower_sign(x: &Qd, y: &Qd, e: &R) -> Ordering {
    let sy = if *e == zero() {
        Ordering::Equal
    } else {
        y.sign()
    };
    let sx = x.sign();
    if sy == Ordering::Equal {
        return sx;
    }
    if sx == Ordering::Equal || sx == sy {
        return sy;
    }
    match x.mul(x).sub(&y.mul(y).scale(e)).sign() {
        Ordering::Greater => sx,
        Ordering::Less => sy,
        Ordering::Equal => Ordering::Equal,
    }
}

/// The exact sign of `sum a_i b_i`, the `a_i` of one field and the `b_i`
/// of another (or the same): `sum a_i p_i + (sum a_i q_i) sqrt e` for
/// `b_i = p_i + q_i sqrt e`.
pub(super) fn mixed_dot_sign(a: &[Qd], b: &[Qd]) -> Ordering {
    let e = b
        .iter()
        .find_map(|x| x.field().cloned())
        .unwrap_or_else(zero);
    let field_a = a
        .iter()
        .find_map(|x| x.field().cloned())
        .unwrap_or_else(zero);
    let mut x = Qd::new(zero(), zero(), field_a.clone());
    let mut y = Qd::new(zero(), zero(), field_a);
    for (ai, bi) in a.iter().zip(b) {
        x = x.add(&ai.scale(&bi.a));
        y = y.add(&ai.scale(&bi.b));
    }
    tower_sign(&x, &y, &e)
}

/// A point or vector of surds of one field.
pub(super) type QV = [Qd; 3];

pub(super) fn qv(p: &V) -> QV {
    [
        Qd::rat(p[0].clone()),
        Qd::rat(p[1].clone()),
        Qd::rat(p[2].clone()),
    ]
}

pub(super) fn qadd(a: &QV, b: &QV) -> QV {
    [a[0].add(&b[0]), a[1].add(&b[1]), a[2].add(&b[2])]
}

pub(super) fn qsub(a: &QV, b: &QV) -> QV {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}

pub(super) fn qscale(a: &V, k: &Qd) -> QV {
    [k.scale(&a[0]), k.scale(&a[1]), k.scale(&a[2])]
}

pub(super) fn qdot(a: &QV, b: &V) -> Qd {
    a[0].scale(&b[0])
        .add(&a[1].scale(&b[1]))
        .add(&a[2].scale(&b[2]))
}

pub(super) fn qqdot(a: &QV, b: &QV) -> Qd {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

pub(super) fn qcross(a: &QV, b: &QV) -> QV {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}

pub(super) fn qv_eq(a: &QV, b: &QV) -> bool {
    a.iter().zip(b).all(|(x, y)| x.cmp(y) == Ordering::Equal)
}

pub(super) fn qv_f64(a: &QV) -> [f64; 3] {
    [a[0].to_f64(), a[1].to_f64(), a[2].to_f64()]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(a: i64, b: i64) -> R {
        R::new(BigInt::from(a), BigInt::from(b))
    }

    #[test]
    fn signs_of_surds() {
        // 1 - sqrt(2) < 0, 3 - 2 sqrt(2) > 0, 3 - sqrt(9) folds to 0.
        assert_eq!(Qd::new(int(1), int(-1), int(2)).sign(), Ordering::Less);
        assert_eq!(Qd::new(int(3), int(-2), int(2)).sign(), Ordering::Greater);
        assert_eq!(Qd::new(int(3), int(-1), int(9)).sign(), Ordering::Equal);
        // sqrt(2) vs sqrt(3) / 1.2 (= 1.4434) and 1.4142.
        let a = Qd::new(zero(), int(1), int(2));
        let b = Qd::new(zero(), r(5, 6), int(3));
        assert_eq!(a.cmp(&b), Ordering::Less);
        assert_eq!(b.cmp(&a), Ordering::Greater);
        // sqrt(2) + sqrt(3) vs 3.146...: (sqrt 2 - 3.15) + sqrt 3 < 0.
        let x = Qd::new(r(-315, 100), int(1), int(2));
        assert_eq!(tower_sign(&x, &Qd::rat(int(1)), &int(3)), Ordering::Less);
        let x = Qd::new(r(-314, 100), int(1), int(2));
        assert_eq!(tower_sign(&x, &Qd::rat(int(1)), &int(3)), Ordering::Greater);
    }

    #[test]
    fn mixed_products() {
        // sqrt2 * sqrt3 - 2.449 > 0, - 2.45 < 0.
        let a = [Qd::new(zero(), int(1), int(2)), Qd::rat(r(-2449, 1000))];
        let b = [Qd::new(zero(), int(1), int(3)), Qd::rat(int(1))];
        assert_eq!(mixed_dot_sign(&a, &b), Ordering::Greater);
        let a = [Qd::new(zero(), int(1), int(2)), Qd::rat(r(-245, 100))];
        assert_eq!(mixed_dot_sign(&a, &b), Ordering::Less);
    }
}
