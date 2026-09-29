//! Rational arithmetic in lowest terms, the same canonical values as
//! `num_rational`'s operators, with a faster gcd.
//!
//! `num_rational` reduces every result by `num_integer`'s gcd of `BigUint`,
//! Stein's binary algorithm, which removes about one bit per multiprecision
//! step and so takes time quadratic in the operands' length, however small
//! the other operand is. Here one Euclidean division balances the operands
//! and Lehmer's algorithm (Knuth, TAOCP 4.5.2, Algorithm L) takes the
//! quotients of their leading 63 bits, a word of multiprecision work for
//! each ~60 bits removed. Products and sums reduce as in Knuth 4.5.1: the
//! gcds of cross terms only, with the result's numerator and denominator
//! coprime and its denominator positive (zero as `0/1`), exactly what
//! `BigRational::new` would give.
use num_bigint::{BigInt, BigUint, Sign};
use num_rational::BigRational as R;

fn is_zero(x: &BigUint) -> bool {
    x.bits() == 0
}

fn is_one(x: &BigInt) -> bool {
    x.sign() == Sign::Plus && x.bits() == 1
}

fn low_u128(x: &BigUint) -> u128 {
    let d = x.to_u64_digits();
    u128::from(d.first().copied().unwrap_or(0)) | (u128::from(d.get(1).copied().unwrap_or(0)) << 64)
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// The gcd of two naturals.
fn gcd_nat(a: &BigUint, b: &BigUint) -> BigUint {
    let (mut u, mut v) = if a >= b {
        (a.clone(), b.clone())
    } else {
        (b.clone(), a.clone())
    };
    loop {
        if is_zero(&v) {
            return u;
        }
        if u.bits() <= 128 {
            return BigUint::from(gcd_u128(low_u128(&u), low_u128(&v)));
        }
        // The leading 63 bits of u and v at one shift; Euclid on them with
        // cofactors while both bounds agree on the quotient (Knuth's L2).
        let shift = u.bits() - 63;
        let top = |x: &BigUint| low_u128(&(x >> shift)) as i128;
        let (mut uh, mut vh) = (top(&u), top(&v));
        let (mut a, mut b, mut c, mut d) = (1i128, 0i128, 0i128, 1i128);
        loop {
            let (n0, d0, n1, d1) = (uh + a, vh + c, uh + b, vh + d);
            if d0 <= 0 || d1 <= 0 || n0 < 0 || n1 < 0 {
                break;
            }
            let q = n0 / d0;
            if q != n1 / d1 {
                break;
            }
            (a, c) = (c, a - q * c);
            (b, d) = (d, b - q * d);
            (uh, vh) = (vh, uh - q * vh);
        }
        if b == 0 {
            let r = &u % &v;
            (u, v) = (v, r);
        } else {
            // The cofactors' matrix is unimodular: the gcd is kept, and the
            // agreed quotients make both images the true remainders.
            let (su, sv) = (BigInt::from(u), BigInt::from(v));
            let nu = &su * a + &sv * b;
            let nv = &su * c + &sv * d;
            let (nu, nv) = (nu.into_parts().1, nv.into_parts().1);
            (u, v) = if nu >= nv { (nu, nv) } else { (nv, nu) };
        }
    }
}

/// The nonnegative gcd of two integers.
pub(crate) fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    BigInt::from(gcd_nat(a.magnitude(), b.magnitude()))
}

fn zero() -> R {
    R::from_integer(BigInt::from(0))
}

fn is_zero_r(x: &R) -> bool {
    x.numer().sign() == Sign::NoSign
}

/// `a b`.
pub(crate) fn mul(a: &R, b: &R) -> R {
    if is_zero_r(a) || is_zero_r(b) {
        return zero();
    }
    let g1 = gcd(a.numer(), b.denom());
    let g2 = gcd(b.numer(), a.denom());
    let part = |x: &BigInt, g: &BigInt| if is_one(g) { x.clone() } else { x / g };
    R::new_raw(
        part(a.numer(), &g1) * part(b.numer(), &g2),
        part(a.denom(), &g2) * part(b.denom(), &g1),
    )
}

/// `a / b`, `b` nonzero.
pub(crate) fn div(a: &R, b: &R) -> R {
    assert!(!is_zero_r(b), "division by zero");
    let inv = if b.numer().sign() == Sign::Minus {
        R::new_raw(-b.denom(), -b.numer())
    } else {
        R::new_raw(b.denom().clone(), b.numer().clone())
    };
    mul(a, &inv)
}

/// `a + b`.
pub(crate) fn add(a: &R, b: &R) -> R {
    if is_zero_r(a) {
        return b.clone();
    }
    if is_zero_r(b) {
        return a.clone();
    }
    let d1 = gcd(a.denom(), b.denom());
    if is_one(&d1) {
        return R::new_raw(
            a.numer() * b.denom() + b.numer() * a.denom(),
            a.denom() * b.denom(),
        );
    }
    let (ad, bd) = (a.denom() / &d1, b.denom() / &d1);
    let t = a.numer() * &bd + b.numer() * &ad;
    if t.sign() == Sign::NoSign {
        return zero();
    }
    let d2 = gcd(&t, &d1);
    if is_one(&d2) {
        R::new_raw(t, ad * b.denom())
    } else {
        R::new_raw(t / &d2, ad * (b.denom() / d2))
    }
}

/// `a - b`.
pub(crate) fn sub(a: &R, b: &R) -> R {
    add(a, &-b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_integer::Integer;

    /// A deterministic stream of integers of assorted lengths and shapes.
    fn samples() -> Vec<BigInt> {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut out = vec![BigInt::from(0), BigInt::from(1), BigInt::from(-1)];
        for k in 0..400 {
            let words = 1 + (next() % 24) as usize;
            let mut x = BigInt::from(0);
            for _ in 0..words {
                x = (x << 64) + BigInt::from(next());
            }
            // Shared factors, powers of two and small values among them.
            match k % 5 {
                0 => x = (x * BigInt::from(3u64.pow(20))) << 37,
                1 => x = BigInt::from(next() % 1000),
                2 => x = -x,
                _ => {}
            }
            out.push(x);
        }
        out
    }

    #[test]
    fn gcds_agree_with_steins() {
        let xs = samples();
        for (i, a) in xs.iter().enumerate() {
            for b in xs.iter().skip(i).step_by(7) {
                let want = a.gcd(b);
                assert_eq!(gcd(a, b), want, "gcd({a}, {b})");
                // Common factors: gcd(a c, b c) = gcd(a, b) c.
                let c = &xs[(i * 31) % xs.len()];
                if c.sign() != Sign::NoSign {
                    assert_eq!(
                        gcd(&(a * c), &(b * c)),
                        want * BigInt::from(c.magnitude().clone())
                    );
                }
            }
        }
    }

    #[test]
    fn operations_agree_with_num_rational() {
        let xs = samples();
        let rs: Vec<R> = xs
            .iter()
            .zip(xs.iter().rev())
            .filter(|(_, d)| d.sign() != Sign::NoSign)
            .map(|(n, d)| R::new(n.clone(), d.clone()))
            .collect();
        for (i, a) in rs.iter().enumerate() {
            for b in rs.iter().skip(i).step_by(5).chain([a, &-a]) {
                let (x, y) = (mul(a, b), a * b);
                assert!(x.numer() == y.numer() && x.denom() == y.denom());
                let (x, y) = (add(a, b), a + b);
                assert!(x.numer() == y.numer() && x.denom() == y.denom());
                let (x, y) = (sub(a, b), a - b);
                assert!(x.numer() == y.numer() && x.denom() == y.denom());
                if b.numer().sign() != Sign::NoSign {
                    let (x, y) = (div(a, b), a / b);
                    assert!(x.numer() == y.numer() && x.denom() == y.denom());
                }
            }
        }
    }
}
