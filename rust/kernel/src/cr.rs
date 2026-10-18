//! Correctly rounded elementary functions: `sin`, `cos`, `tan`, `atan`,
//! `atan2`, `asin`, `acos`, `hypot`, `exp`, `ln`, `log2`, `sinh`, `cosh`,
//! `asinh` and `atanh`, each the exact value rounded once to the nearest binary64.
//! No transcendental value at a binary64 argument is a tie (other than the
//! exact cases each function returns first); `hypot` decides its ties
//! exactly, to even.
//!
//! Platform libms differ in the last place: macOS's `hypot` and `tan`
//! against glibc's, for two (the review of 18 October 2026, F11 and X1).
//! A kernel construction that called them built different bits on
//! different hosts, and a Boolean's history differed between them. R13: no
//! kernel construction calls the platform's libm; every one calls these
//! functions, and every host builds the same bits from the same inputs.
//! `sqrt` and the four operations are correctly rounded by IEEE 754
//! everywhere and stay as they are.
//!
//! Two tiers, as the certified decisions'. The `Fast` tier evaluates in
//! double-double with an explicit error bound, and decides the rounding
//! when the bound keeps the value off both of its rounding boundaries. Its
//! error-free products use the fused multiply-add instruction where the
//! target has one and Dekker's product elsewhere (the same exact errors;
//! never a libm `fma` call). Its passes go quickest first: sine and cosine
//! from tables of `j/256` with their leading terms exact (all arguments but
//! about one in 250), then the tables in double-double (one in a thousand
//! of those), then series in double-double (one in `2^30`); the
//! arctangent's from tables, then series. Otherwise the rational tier
//! (`certified.rs`: Machin's pi, Taylor sine and cosine, `atan`, `ln` and
//! the integer square root on a `2^-200` grid; an exponential series here)
//! encloses the value and produces the bits; `hypot` decides on the exact
//! squares, in 128-bit integers for the ties of Pythagorean data. Sines and
//! cosines past `2^20 pi/2`, or too near a multiple of `pi/2` for the
//! first reduction, are reduced again with 1344 bits of `pi/2` before
//! either tier. Where neither tier decides, which needs a value within
//! `2^-120` of a rounding boundary (beyond the hardest cases known for
//! binary64) or a result in the subnormal range, the result is the
//! enclosure's nearest, still the same on every host.
//!
//! The functions are free functions and, through [`Rn`], methods named
//! after `f64`'s with `_rn` (round to nearest), as the fixture generators'
//! `cos_rn`, `sin_rn`, `atan2_rn` and `hypot_rn` are.
use crate::certified::{self, Interval};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::OnceLock;

// ------------------------------------------------------------- binary64

/// `2^k` for `k` in the normal range.
const fn pow2(k: i32) -> f64 {
    f64::from_bits(((k + 1023) as u64) << 52)
}

/// Relative error budget of a double-double evaluation. Every chain below
/// stays under `2^-95` (each double-double operation under `2^-102`,
/// Joldes, Muller and Popescu 2017, and at most about twenty in a chain);
/// the budget leaves a factor of 128 more.
const FAST: f64 = pow2(-88);

/// Relative error budget of the tabled evaluations (sine, cosine and
/// arctangent from tables of `j/256`, a short series in binary64 beyond
/// the leading term): their analyses stay under `2^-70`; the budget
/// leaves a factor of 64. They decide all but about one argument in a
/// thousand; the series' evaluation then decides.
const TABLED: f64 = pow2(-64);

/// Relative error budget of the quick sine and cosine (the tables'
/// leading terms exact, every other term in binary64): their analysis stays
/// under `2^-66`; the budget leaves a factor of 16. They decide all but
/// about one argument in 250; the tabled evaluation then decides.
const QUICK: f64 = pow2(-62);

fn next_up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits(1);
    }
    let b = x.to_bits();
    f64::from_bits(if x > 0.0 { b + 1 } else { b - 1 })
}

fn next_down(x: f64) -> f64 {
    -next_up(-x)
}

/// `x 2^k` in two steps, each factor representable: exact while the
/// result is normal.
fn ldexp(x: f64, k: i32) -> f64 {
    let k = k.clamp(-2200, 2200);
    let step = |mut v: f64, mut e: i32| {
        while e > 1000 {
            v *= pow2(1000);
            e -= 1000;
        }
        while e < -1000 {
            v *= pow2(-1000);
            e += 1000;
        }
        v * pow2(e)
    };
    let half = k / 2;
    step(step(x, half), k - half)
}

/// The binary exponent of a positive finite number: `x` in `[2^e, 2^(e+1))`.
fn exponent(x: f64) -> i32 {
    if x < f64::MIN_POSITIVE {
        return exponent(x * pow2(54)) - 54;
    }
    ((x.to_bits() >> 52) & 0x7ff) as i32 - 1023
}

/// `1/n!` in binary64, by successive divisions (the tails' coefficients,
/// whose own rounding is far below the tails' weight).
const fn inv_fact(n: u32) -> f64 {
    let mut p = 1.0;
    let mut m = 1;
    while m <= n {
        p /= m as f64;
        m += 1;
    }
    p
}

// ---------------------------------------------- error-free transformations

/// `a + b` and its exact error (Knuth).
#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

/// `a + b` and its exact error for `|a| >= |b|` (Dekker).
#[inline]
fn fast_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

/// Veltkamp's split: `a = hi + lo`, each of at most 26 bits.
#[inline]
#[cfg(not(any(
    target_arch = "aarch64",
    all(target_arch = "x86_64", target_feature = "fma")
)))]
fn split(a: f64) -> (f64, f64) {
    let c = 134_217_729.0 * a;
    let hi = c - (c - a);
    (hi, a - hi)
}

/// `a b` and its exact error, away from overflow and underflow: by the
/// fused multiply-add instruction where the target has one, else by
/// Dekker's product. Both give the exact error, so the bits are the same
/// either way (no libm `fma` call: Rust lowers `mul_add` to the
/// instruction on these targets).
#[inline]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    #[cfg(any(
        target_arch = "aarch64",
        all(target_arch = "x86_64", target_feature = "fma")
    ))]
    {
        (p, a.mul_add(b, -p))
    }
    #[cfg(not(any(
        target_arch = "aarch64",
        all(target_arch = "x86_64", target_feature = "fma")
    )))]
    {
        let (ah, al) = split(a);
        let (bh, bl) = split(b);
        (p, ((ah * bh - p) + ah * bl + al * bh) + al * bl)
    }
}

/// A double-double `hi + lo`, `|lo|` at most half an ulp of `hi`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Dd {
    hi: f64,
    lo: f64,
}

const fn dd(hi: u64, lo: u64) -> Dd {
    Dd {
        hi: f64::from_bits(hi),
        lo: f64::from_bits(lo),
    }
}

impl Dd {
    const fn f(x: f64) -> Self {
        Self { hi: x, lo: 0.0 }
    }
    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
    /// AccurateDWPlusDW: relative error below `3 u^2`.
    fn add(self, o: Self) -> Self {
        let (sh, sl) = two_sum(self.hi, o.hi);
        let (th, tl) = two_sum(self.lo, o.lo);
        let (vh, vl) = fast_two_sum(sh, sl + th);
        let (hi, lo) = fast_two_sum(vh, tl + vl);
        Self { hi, lo }
    }
    fn sub(self, o: Self) -> Self {
        self.add(o.neg())
    }
    /// DWPlusFP: relative error below `2 u^2`.
    fn add_f(self, b: f64) -> Self {
        let (sh, sl) = two_sum(self.hi, b);
        let (hi, lo) = fast_two_sum(sh, self.lo + sl);
        Self { hi, lo }
    }
    /// DWTimesDW1: relative error below `7 u^2`.
    fn mul(self, o: Self) -> Self {
        let (ch, cl) = two_prod(self.hi, o.hi);
        let t = self.hi * o.lo + self.lo * o.hi;
        let (hi, lo) = fast_two_sum(ch, cl + t);
        Self { hi, lo }
    }
    /// DWTimesFP1: relative error below `2 u^2`.
    fn mul_f(self, b: f64) -> Self {
        let (ch, cl1) = two_prod(self.hi, b);
        let (th, tl) = fast_two_sum(ch, self.lo * b);
        let (hi, lo) = fast_two_sum(th, tl + cl1);
        Self { hi, lo }
    }
    /// DWDivDW2: relative error below `15 u^2 + 56 u^3`.
    fn div(self, o: Self) -> Self {
        let th = self.hi / o.hi;
        let r = o.mul_f(th);
        let d = (self.hi - r.hi) + (self.lo - r.lo);
        let (hi, lo) = fast_two_sum(th, d / o.hi);
        Self { hi, lo }
    }
    /// One Newton step from the binary64 root: relative error below `4 u^2`.
    fn sqrt(self) -> Self {
        let s = self.hi.sqrt();
        if s == 0.0 {
            return Self::f(0.0);
        }
        let (p, e) = two_prod(s, s);
        let d = (((self.hi - p) - e) + self.lo) / (2.0 * s);
        let (hi, lo) = fast_two_sum(s, d);
        Self { hi, lo }
    }
    /// `self 2^k`, exact away from underflow.
    fn scale(self, k: i32) -> Self {
        Self {
            hi: ldexp(self.hi, k),
            lo: ldexp(self.lo, k),
        }
    }
}

/// The nearest binary64 to every value within `err` of `x`, if they share
/// one: the rounding boundaries, half an ulp each side of `x.hi + x.lo`,
/// lie beyond the bound and the comparison's own rounding.
fn decide(x: Dd, err: f64) -> Option<f64> {
    let c = x.hi + x.lo;
    if !c.is_finite() || c == 0.0 || err.is_nan() {
        return None;
    }
    // The gaps to the neighbours away from zero and towards it, exactly
    // (half as large towards zero at a power of two): from the exponent
    // field, the subnormals' a fixed `2^-1074`.
    let bits = c.to_bits() & !(1u64 << 63);
    let field = (bits >> 52) as i32;
    // 2^(field - 1075): normal from field 53 on, subnormal below.
    let away = if field >= 53 {
        f64::from_bits(((field - 52) as u64) << 52)
    } else {
        f64::from_bits(1u64 << (field - 1).max(0))
    };
    let toward = if field > 1 && bits & ((1 << 52) - 1) == 0 {
        0.5 * away
    } else {
        away
    };
    // `x - c` (`x.hi - c` exactly: two neighbouring numbers), positive
    // away from zero.
    let d = (x.hi - c) + x.lo;
    let u = if c < 0.0 { -d } else { d };
    let slack = 2.0 * (err + u.abs() * pow2(-50) + f64::from_bits(4));
    (u + slack < 0.5 * away && u - slack > -0.5 * toward).then_some(c)
}

/// Horner in binary64 over `c`: `sum c_n z^n`, a series' tail.
fn tail(z: f64, c: &[f64]) -> f64 {
    c.iter().rev().fold(0.0, |p, &k| p * z + k)
}

/// Horner in double-double over the head `c`, from the tail's value `t`
/// (its terms divided by `z^len`): `sum_{n < len} c_n z^n + z^len t`.
fn head(z: Dd, c: &[Dd], t: f64) -> Dd {
    c.iter().rev().fold(Dd::f(t), |p, &k| p.mul(z).add(k))
}

// ------------------------------------------------------------ constants

/// `1/n!`, `n = 0..=13`.
const INV_FACT: [Dd; 14] = [
    dd(0x3ff0000000000000, 0x0000000000000000),
    dd(0x3ff0000000000000, 0x0000000000000000),
    dd(0x3fe0000000000000, 0x0000000000000000),
    dd(0x3fc5555555555555, 0x3c65555555555555),
    dd(0x3fa5555555555555, 0x3c45555555555555),
    dd(0x3f81111111111111, 0x3c01111111111111),
    dd(0x3f56c16c16c16c17, 0xbbef49f49f49f49f),
    dd(0x3f2a01a01a01a01a, 0x3b6a01a01a01a01a),
    dd(0x3efa01a01a01a01a, 0x3b3a01a01a01a01a),
    dd(0x3ec71de3a556c734, 0xbb6c154f8ddc6c00),
    dd(0x3e927e4fb7789f5c, 0x3b3cbbc05b4fa99a),
    dd(0x3e5ae64567f544e4, 0xbafc062e06d1f209),
    dd(0x3e21eed8eff8d898, 0xbac2aec959e14c06),
    dd(0x3de6124613a86d09, 0x3a8f28e0cc748ebe),
];

/// `1/(2k+1)`, `k = 0..=6`.
const INV_ODD: [Dd; 7] = [
    dd(0x3ff0000000000000, 0x0000000000000000),
    dd(0x3fd5555555555555, 0x3c75555555555555),
    dd(0x3fc999999999999a, 0xbc6999999999999a),
    dd(0x3fc2492492492492, 0x3c62492492492492),
    dd(0x3fbc71c71c71c71c, 0x3c5c71c71c71c71c),
    dd(0x3fb745d1745d1746, 0xbc4745d1745d1746),
    dd(0x3fb3b13b13b13b14, 0xbc53b13b13b13b14),
];

/// `pi/2` in four parts of 33, 33, 33 and 53 bits: `k` times each of the
/// first three is exact for `|k| < 2^20`; the rest is below `2^-159`.
const HALF_PI_PARTS: [f64; 4] = [
    f64::from_bits(0x3ff921fb54400000),
    f64::from_bits(0x3dd0b4611a600000),
    f64::from_bits(0x3ba3198a2e000000),
    f64::from_bits(0x397b839a252049c1),
];
const PI: Dd = dd(0x400921fb54442d18, 0x3ca1a62633145c07);
const HALF_PI: Dd = dd(0x3ff921fb54442d18, 0x3c91a62633145c07);

/// `atan(j/16)`, `j = 0..=16`.
const ATAN_TABLE: [Dd; 17] = [
    dd(0x0000000000000000, 0x0000000000000000),
    dd(0x3faff55bb72cfdea, 0xbc3c934d86d23f1d),
    dd(0x3fbfd5ba9aac2f6e, 0xbc4cd37686760c17),
    dd(0x3fc7b97b4bce5b02, 0x3c5347b0b4f881ca),
    dd(0x3fcf5b75f92c80dd, 0x3c68ab6e3cf7afbd),
    dd(0x3fd362773707ebcc, 0xbc6963a544b672d8),
    dd(0x3fd6f61941e4def1, 0xbc7c63aae6f6e918),
    dd(0x3fda64eec3cc23fd, 0xbc724dec1b50b7ff),
    dd(0x3fddac670561bb4f, 0x3c7a2b7f222f65e2),
    dd(0x3fe0657e94db30d0, 0xbc7d5b495f6349e6),
    dd(0x3fe1e00babdefeb4, 0xbc5928df287a668f),
    dd(0x3fe345f01cce37bb, 0x3c81021137c71102),
    dd(0x3fe4978fa3269ee1, 0x3c72419a87f2a458),
    dd(0x3fe5d58987169b18, 0x3c60028e4bc5e7ca),
    dd(0x3fe700a7c5784634, 0xbc78c34d25aadef6),
    dd(0x3fe819d0b7158a4d, 0xbc7bf76229d3b917),
    dd(0x3fe921fb54442d18, 0x3c81a62633145c07),
];

/// `ln 2` in three parts of 42, 42 and 53 bits: `k` times each of the
/// first two is exact for `|k| < 2^11`; the rest is below `2^-139`.
const LN2_PARTS: [f64; 3] = [
    f64::from_bits(0x3fe62e42fefa3800),
    f64::from_bits(0x3d2ef35793c76000),
    f64::from_bits(0x3a8cc01f97b57a08),
];
const LN2: Dd = dd(0x3fe62e42fefa39ef, 0x3c7abc9e3b39803f);

/// `(-1)^n / (2n+1)!` and `(-1)^n / (2n)!`, `n = 7..=14`: the sine's and
/// cosine's tails.
const SIN_TAIL: [f64; 8] = {
    let mut a = [0.0; 8];
    let mut i = 0;
    while i < 8 {
        let n = (i + 7) as u32;
        let v = inv_fact(2 * n + 1);
        a[i] = if n % 2 == 0 { v } else { -v };
        i += 1;
    }
    a
};
const COS_TAIL: [f64; 8] = {
    let mut a = [0.0; 8];
    let mut i = 0;
    while i < 8 {
        let n = (i + 7) as u32;
        let v = inv_fact(2 * n);
        a[i] = if n % 2 == 0 { v } else { -v };
        i += 1;
    }
    a
};
/// `(-1)^n / (2n+1)`, `n = 4..=10`: the arctangent's tail.
const ATAN_TAIL: [f64; 7] = {
    let mut a = [0.0; 7];
    let mut i = 0;
    while i < 7 {
        let n = i + 4;
        let v = 1.0 / (2 * n + 1) as f64;
        a[i] = if n % 2 == 0 { v } else { -v };
        i += 1;
    }
    a
};
/// `1 / (2n+1)`, `n = 7..=20`: the logarithm's tail.
const LN_TAIL: [f64; 14] = {
    let mut a = [0.0; 14];
    let mut i = 0;
    while i < 14 {
        a[i] = 1.0 / (2 * (i + 7) + 1) as f64;
        i += 1;
    }
    a
};
/// `1/n!`, `n = 11..=22`: the exponential's tail.
const EXP_TAIL: [f64; 12] = {
    let mut a = [0.0; 12];
    let mut i = 0;
    while i < 12 {
        a[i] = inv_fact(i as u32 + 11);
        i += 1;
    }
    a
};
/// `1/(2n+1)!` and `1/(2n)!`, `n = 6..=12`: the hyperbolic sine's and
/// cosine's tails near zero.
const SINH_TAIL: [f64; 7] = {
    let mut a = [0.0; 7];
    let mut i = 0;
    while i < 7 {
        a[i] = inv_fact(2 * (i as u32 + 6) + 1);
        i += 1;
    }
    a
};
const COSH_TAIL: [f64; 7] = {
    let mut a = [0.0; 7];
    let mut i = 0;
    while i < 7 {
        a[i] = inv_fact(2 * (i as u32 + 6));
        i += 1;
    }
    a
};

// ---------------------------------------------------- the rational tier

fn exact(x: f64) -> R {
    R::from_float(x).expect("a finite binary64")
}

fn int(n: i64) -> R {
    R::from_integer(BigInt::from(n))
}

fn pow2_r(k: i64) -> R {
    if k >= 0 {
        R::from_integer(BigInt::from(1) << k as usize)
    } else {
        R::new(BigInt::from(1), BigInt::from(1) << (-k) as usize)
    }
}

/// A binary64 within a few ulps of a nonzero rational in `[2^-1000, 2^1000]`.
fn approx(q: &R) -> Option<f64> {
    let lead = |v: &BigInt| -> (f64, i64) {
        let shift = (v.bits() as i64 - 64).max(0);
        let top = v.magnitude() >> shift as usize;
        (
            top.to_u64_digits().first().copied().unwrap_or(0) as f64,
            shift,
        )
    };
    if q.numer().sign() == Sign::NoSign {
        return None;
    }
    let ((n, sn), (d, sd)) = (lead(q.numer()), lead(q.denom()));
    let e = sn - sd;
    if !(-1100..=1100).contains(&e) {
        return None;
    }
    let f = ldexp(n / d, e as i32);
    let f = if q.numer().sign() == Sign::Minus {
        -f
    } else {
        f
    };
    (f.is_finite() && f.abs() >= pow2(-1000) && f.abs() <= pow2(1000)).then_some(f)
}

/// The binary64 nearest a rational, `None` at a tie or outside the range.
fn nearest(q: &R) -> Option<f64> {
    let mid = |a: f64, b: f64| (exact(a) + exact(b)) / int(2);
    let mut c = approx(q)?;
    for _ in 0..16 {
        let (down, up) = (next_down(c), next_up(c));
        match q.cmp(&mid(c, up)) {
            Ordering::Greater => {
                c = up;
                continue;
            }
            Ordering::Equal => return None,
            Ordering::Less => {}
        }
        match q.cmp(&mid(down, c)) {
            Ordering::Less => {
                c = down;
                continue;
            }
            Ordering::Equal => return None,
            Ordering::Greater => {}
        }
        return Some(c);
    }
    None
}

/// The binary64 every value of the enclosure rounds to, if one.
fn round_enclosure(i: &Interval) -> Option<f64> {
    let c = nearest(i.lo())?;
    (nearest(i.hi())? == c).then_some(c)
}

/// The rational tier's binary64: the enclosure's if it decides, else its
/// midpoint's nearest (exact arithmetic still, the same on every host),
/// else `fallback`.
fn rational(i: Option<Interval>, fallback: f64) -> f64 {
    let Some(i) = i else {
        return fallback;
    };
    round_enclosure(&i)
        .or_else(|| approx(&((i.lo() + i.hi()) / int(2))))
        .unwrap_or(fallback)
}

/// `i` widened by `e` each way.
fn widened(i: &Interval, e: f64) -> Interval {
    let e = exact(e);
    Interval::new(i.lo() - &e, i.hi() + &e)
}

/// `exp(x) = 2^k m` with an enclosure of `m` (`|x - k ln 2| <= ~0.35`), by
/// a Taylor series and its remainder bound.
fn exp_rational(x: f64) -> (i32, Interval) {
    let ln2 = certified::ln2();
    let mid = (ln2.lo() + ln2.hi()) / int(2);
    let k = (exact(x) / mid).round().to_integer();
    let ki: i32 = i32::try_from(&k).unwrap_or(0);
    let r = Interval::from_f64(x).sub(&ln2.scale(&R::from_integer(k)));
    let mut sum = Interval::exact(int(1));
    let mut term = sum.clone();
    const TERMS: i64 = 40;
    for n in 1..TERMS {
        term = term
            .mul(&r)
            .scale(&R::new(BigInt::from(1), BigInt::from(n)));
        sum = sum.add(&term);
    }
    // The remainder: |r|^N / N! times a geometric factor below 2.
    let rem = term.mul(&r).abs_hi() * R::new(BigInt::from(2), BigInt::from(TERMS));
    (ki, sum.add(&Interval::new(-rem.clone(), rem)))
}

// ------------------------------------------------------------- hypot

/// `sqrt(x^2 + y^2)` correctly rounded (IEEE 754's `hypot`: infinite when
/// either is, else NaN when either is).
pub fn hypot(x: f64, y: f64) -> f64 {
    if x.is_infinite() || y.is_infinite() {
        return f64::INFINITY;
    }
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    let (a, b) = if x.abs() >= y.abs() {
        (x.abs(), y.abs())
    } else {
        (y.abs(), x.abs())
    };
    if b == 0.0 {
        return a;
    }
    // b < 2^-27 a: the excess, a (b/a)^2 / 2, is under a quarter ulp.
    if b < a * pow2(-27) {
        return a;
    }
    if !(pow2(-400)..pow2(400)).contains(&a) {
        return hypot_exact(a, b);
    }
    // Unscaled: the squares and their errors stay normal.
    let (p, pe) = two_prod(a, a);
    let (q, qe) = two_prod(b, b);
    let r = Dd { hi: p, lo: pe }.add(Dd { hi: q, lo: qe }).sqrt();
    if let Some(c) = decide(r, FAST * r.hi) {
        return c;
    }
    let e = exponent(a);
    let (a1, b1) = (ldexp(a, -e), ldexp(b, -e));
    if b1 >= pow2(-9) {
        ldexp(hypot_scaled(a1, b1, ldexp(r.hi + r.lo, -e)), e)
    } else {
        hypot_exact(a, b)
    }
}

/// `hypot` of `a` in `[1, 2)` and `b` in `[2^-9, a]` from their exact
/// squares in 128-bit integers (the ties and near ties the double-double
/// leaves; Pythagorean data make them common): on the `2^-61` grid `a`,
/// `b` and twice every candidate's midpoints are integers below `2^64`,
/// their squares below `2^128`. Ties to even.
fn hypot_scaled(a: f64, b: f64, guess: f64) -> f64 {
    let grid = |x: f64| (x * pow2(61)) as u128;
    let s = grid(a) * grid(a) + grid(b) * grid(b);
    // The midpoint between c and its neighbour d, squared, on the grid.
    let mid2 = |c: f64, d: f64| {
        let m = grid(c) + grid(d);
        // (m / 2)^2 = m^2 / 4: compare 4 s with m^2 instead.
        m * m
    };
    let s4 = s << 2;
    let mut c = guess;
    for _ in 0..8 {
        let (up, down) = (next_up(c), next_down(c));
        match s4.cmp(&mid2(c, up)) {
            Ordering::Greater => {
                c = up;
                continue;
            }
            Ordering::Equal => return if c.to_bits() & 1 == 0 { c } else { up },
            Ordering::Less => {}
        }
        match s4.cmp(&mid2(down, c)) {
            Ordering::Less => {
                c = down;
                continue;
            }
            Ordering::Equal => return if c.to_bits() & 1 == 0 { c } else { down },
            Ordering::Greater => {}
        }
        return c;
    }
    c
}

/// A finite nonzero binary64's integer significand and exponent: `|x| = m 2^e`.
fn parts(x: f64) -> (BigInt, i64) {
    let b = x.abs().to_bits();
    let (field, frac) = ((b >> 52) as i64, b & ((1 << 52) - 1));
    if field == 0 {
        (BigInt::from(frac), -1074)
    } else {
        (BigInt::from(frac | (1 << 52)), field - 1075)
    }
}

/// `hypot` from the exact squares in integers: the binary64 whose rounding
/// cell's squared ends bracket `a^2 + b^2`, ties to even.
fn hypot_exact(a: f64, b: f64) -> f64 {
    let ((ma, ea), (mb, eb)) = (parts(a), parts(b));
    // a^2 + b^2 = s 2^es exactly.
    let es = 2 * ea.min(eb);
    let s = ((&ma * &ma) << (2 * ea - es) as usize) + ((&mb * &mb) << (2 * eb - es) as usize);
    // sign(m^2 - a^2 - b^2) for a midpoint m = n 2^em.
    let cmp = |n: BigInt, em: i64| -> Ordering {
        let e = (2 * em).min(es);
        let lhs = (&n * &n) << (2 * em - e) as usize;
        lhs.cmp(&(&s << (es - e) as usize))
    };
    let mut c = (a * (1.0 + (b / a) * (b / a)).sqrt()).max(a);
    if !c.is_finite() {
        c = f64::MAX;
    }
    for _ in 0..64 {
        let (mc, ec) = parts(c);
        let up = next_up(c);
        // The upper midpoint (the overflow threshold past the largest).
        match cmp(&mc * 2 + 1, ec - 1) {
            Ordering::Less => {
                if up.is_infinite() {
                    return up;
                }
                c = up;
                continue;
            }
            Ordering::Equal => {
                return if c.to_bits() & 1 == 0 { c } else { up };
            }
            Ordering::Greater => {}
        }
        let down = next_down(c);
        // The lower midpoint, a quarter ulp below a power of two.
        let lower = if mc == BigInt::from(1u64 << 52) && ec > -1074 {
            cmp(&mc * 4 - 1, ec - 2)
        } else {
            cmp(&mc * 2 - 1, ec - 1)
        };
        match lower {
            Ordering::Greater => {
                c = down;
                continue;
            }
            Ordering::Equal => return if c.to_bits() & 1 == 0 { c } else { down },
            Ordering::Less => {}
        }
        return c;
    }
    c
}

// ------------------------------------------------------ sine and cosine

/// `x - k pi/2` in double-double for `|k| < 2^20` (Cody and Waite), and
/// its absolute error.
fn reduce(x: f64, k: f64) -> (Dd, f64) {
    if k == 0.0 {
        return (Dd::f(x), 0.0);
    }
    let [p1, p2, p3, p4] = HALF_PI_PARTS;
    // x - k p1 is exact (Sterbenz); k p1, k p2 and k p3 are exact.
    let t = x - k * p1;
    let (h, l) = two_sum(t, -(k * p2));
    let r = Dd { hi: h, lo: l }.add_f(-(k * p3)).add_f(-(k * p4));
    // The rest of pi/2 (2^-159 |k|), k p4's rounding (2^-150 |k|), the
    // double-double sums' (4 u^2 |r|).
    let err = k.abs() * pow2(-148) + r.hi.abs() * pow2(-103);
    (r, err)
}

/// `x - k pi/2` for the quick pass: the last two parts' products (below
/// `2^-48`) summed in binary64 into the low part, an absolute error under
/// `2^-98` more.
fn reduce_quick(x: f64, k: f64) -> (Dd, f64) {
    if k == 0.0 {
        return (Dd::f(x), 0.0);
    }
    let [p1, p2, p3, p4] = HALF_PI_PARTS;
    let t = x - k * p1;
    let (h, l) = two_sum(t, -(k * p2));
    let (h, l) = two_sum(h, l - (k * p3 + k * p4));
    (Dd { hi: h, lo: l }, pow2(-98))
}

/// Bits of `pi/2` for the second reduction.
const BIG: usize = 1344;

/// `pi/2` as an integer over `2^BIG`, within two units: Machin's formula,
/// `pi/2 = 8 atan(1/5) - 2 atan(1/239)`, in integers with 16 guard bits
/// (each series' terms truncate, a unit each, some 300 of them).
fn half_pi_fixed() -> &'static BigInt {
    static P: OnceLock<BigInt> = OnceLock::new();
    P.get_or_init(|| {
        let one = BigInt::from(1) << (BIG + 16);
        let atan_inv = |n: i64| {
            let n2 = BigInt::from(n * n);
            let mut term = &one / BigInt::from(n);
            let mut sum = BigInt::from(0);
            let mut k: i64 = 0;
            while term.sign() != Sign::NoSign {
                let t = &term / BigInt::from(2 * k + 1);
                if k % 2 == 0 {
                    sum += t;
                } else {
                    sum -= t;
                }
                term /= &n2;
                k += 1;
            }
            sum
        };
        (atan_inv(5) * 8 - atan_inv(239) * 2) >> 16
    })
}

/// `x` reduced with `BIG` bits of `pi/2`: `x - k pi/2` exactly as a
/// rational of those bits, `k mod 4`, the bound on the rational's error
/// (two units of `pi/2`'s last bit, `|k|` times) and its double-double.
struct Reduced {
    r: R,
    quadrant: i64,
    err: f64,
    dd: Dd,
}

fn reduce_precise(x: f64) -> Reduced {
    let p = half_pi_fixed();
    let xi = (exact(x) * pow2_r(BIG as i64)).to_integer();
    let two_p = p * BigInt::from(2);
    let k = (&xi * BigInt::from(2) + p).div_floor(&two_p);
    let r = R::new(&xi - &k * p, BigInt::from(1) << BIG);
    let quadrant: i64 = i64::try_from(k.mod_floor(&BigInt::from(4))).unwrap_or(0);
    let e = k.bits() as i64 + 2 - BIG as i64;
    let err = pow2(e.max(-1000) as i32);
    let hi = approx(&r).unwrap_or(0.0);
    let lo = if hi == 0.0 {
        0.0
    } else {
        approx(&(&r - exact(hi))).unwrap_or(0.0)
    };
    let (h, l) = two_sum(hi, lo);
    Reduced {
        r,
        quadrant,
        err,
        dd: Dd { hi: h, lo: l },
    }
}

/// Sine and cosine of a reduced argument `|r| < 0.8`, each to `2^-96`
/// relative: Taylor in `z = r^2`, heads to `z^6` in double-double, tails
/// to `z^14` (`z^7/15! < 2^-44`) in binary64.
fn sin_cos_reduced(r: Dd) -> (Dd, Dd) {
    let z = r.mul(r);
    let signed = |c: Dd, n: usize| if n % 2 == 0 { c } else { c.neg() };
    let sin_head: [Dd; 7] = std::array::from_fn(|n| signed(INV_FACT[2 * n + 1], n));
    let cos_head: [Dd; 7] = std::array::from_fn(|n| signed(INV_FACT[2 * n], n));
    let s = head(z, &sin_head, tail(z.hi, &SIN_TAIL));
    let c = head(z, &cos_head, tail(z.hi, &COS_TAIL));
    (r.mul(s), c)
}

/// The tables: `sin` and `cos` of `j/256` for `j = 0..=205` (every
/// reduced argument's), `atan(j/256)` for `j = 0..=256`, each from the
/// series in double-double (`2^-96`), computed once.
struct Tables {
    sin_cos: Vec<(Dd, Dd)>,
    atan: Vec<Dd>,
}

fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| Tables {
        sin_cos: (0..=205)
            .map(|j| sin_cos_reduced(Dd::f(f64::from(j) / 256.0)))
            .collect(),
        atan: (0..=256)
            .map(|j| atan_unit(Dd::f(f64::from(j) / 256.0)))
            .collect(),
    })
}

/// Sine and cosine of a reduced argument `|r| < 0.8` from the tables:
/// `r = a + b`, `a = j/256`, `|b| <= 1/512`; `sin b = b + d` and
/// `cos b = 1 - b^2/2 + e` with `d` (`<= 2^-20 |b|`) and `e` in binary64
/// (`b`'s low part through `d`'s first term), `b^2` exact; then the sum
/// formulas in double-double. Each within `2^-70` relative.
fn sin_cos_tabled(r: Dd) -> (Dd, Dd) {
    let j = (r.hi * 256.0).round().clamp(-205.0, 205.0);
    let b = r.add_f(-j / 256.0);
    let (bh, bl) = (b.hi, b.lo);
    let z = bh * bh;
    let d = -(bh * z) * (1.0 / 6.0 - z * (1.0 / 120.0 - z * (1.0 / 5040.0))) - 0.5 * z * bl;
    let sb = b.add_f(d);
    let (p, pe) = two_prod(bh, bh);
    let (h, l) = fast_two_sum(p, pe + 2.0 * bh * bl);
    let e = z * z * (1.0 / 24.0 - z * (1.0 / 720.0 - z * (1.0 / 40320.0)));
    let cb = Dd::f(1.0).sub(Dd { hi: h, lo: l }.scale(-1)).add_f(e);
    if j == 0.0 {
        return (sb, cb);
    }
    let (sa, ca) = tables().sin_cos[j.abs() as usize];
    let sa = if j < 0.0 { sa.neg() } else { sa };
    (sa.mul(cb).add(ca.mul(sb)), ca.mul(cb).sub(sa.mul(sb)))
}

/// Sine and cosine of a reduced argument `|r| < 0.8` from the tables, the
/// leading terms exact and the rest in binary64: `r = a + b`, `a = j/256`,
/// `|b| <= 1/512`, `sin(a + b) = Sa + Ca b + (Sa (cos b - 1) + Ca (sin b -
/// b))` with `Sa + Ca b` exact in double-double (`Ca b` by Dekker's product)
/// and the rest, below `2^-18` of it, in binary64; the cosine alike. Each
/// within `2^-66` relative.
fn sin_cos_quick(r: Dd) -> (Dd, Dd) {
    let j = (r.hi * 256.0).round().clamp(-205.0, 205.0);
    let b = r.add_f(-j / 256.0);
    let (bh, bl) = (b.hi, b.lo);
    let z = bh * bh;
    // cos b - 1 and sin b - b, their high parts' terms past the first.
    let cm1 = -0.5 * z + z * z * (1.0 / 24.0 - z * (1.0 / 720.0));
    let sm1 = -(bh * z) * (1.0 / 6.0 - z * (1.0 / 120.0)) - 0.5 * z * bl;
    if j == 0.0 {
        let (h, l) = two_sum(1.0, cm1);
        return (b.add_f(sm1), Dd { hi: h, lo: l });
    }
    let (sa, ca) = tables().sin_cos[j.abs() as usize];
    let sa = if j < 0.0 { sa.neg() } else { sa };
    // sin: Sa + Ca b, then the small terms.
    let (p, pe) = two_prod(ca.hi, bh);
    let (s0, s1) = two_sum(sa.hi, p);
    let ts = sa.lo + pe + ca.hi * bl + ca.lo * bh + sa.hi * cm1 + ca.hi * sm1;
    let (sh, sl) = fast_two_sum(s0, s1 + ts);
    // cos: Ca - Sa b, then the small terms.
    let (q, qe) = two_prod(sa.hi, bh);
    let (c0, c1) = two_sum(ca.hi, -q);
    let tc = ca.lo - qe - sa.hi * bl - sa.lo * bh + ca.hi * cm1 - sa.hi * sm1;
    let (ch, cl) = fast_two_sum(c0, c1 + tc);
    (Dd { hi: sh, lo: sl }, Dd { hi: ch, lo: cl })
}

/// `(sin, cos)` of `x` from those of its reduced argument and `k mod 4`.
fn quadrant<T>(q: i64, (s, c): (T, T), neg: impl Fn(T) -> T) -> (T, T) {
    match q.rem_euclid(4) {
        0 => (s, c),
        1 => (c, neg(s)),
        2 => (neg(s), neg(c)),
        _ => (neg(c), s),
    }
}

/// The evaluations of sine and cosine, quickest first.
#[derive(Clone, Copy)]
enum Pass {
    Quick,
    Tabled,
    Series,
}

/// Double-double sine and cosine with their absolute errors: by Cody and
/// Waite for `|k| < 2^20`, else (`second`) from [`reduce_precise`]; by
/// the pass's evaluation.
fn sin_cos_dd(x: f64, second: Option<&Reduced>, pass: Pass) -> Option<((Dd, f64), (Dd, f64))> {
    let (r, er, q) = match second {
        Some(red) => (red.dd, red.err + red.dd.hi.abs() * pow2(-104), red.quadrant),
        None => {
            let k = (x * std::f64::consts::FRAC_2_PI).round();
            if k.abs() >= pow2(20) {
                return None;
            }
            let (r, er) = match pass {
                Pass::Quick => reduce_quick(x, k),
                _ => reduce(x, k),
            };
            (r, er, k as i64)
        }
    };
    let ((s, c), rel) = match pass {
        Pass::Quick => (sin_cos_quick(r), QUICK),
        Pass::Tabled => (sin_cos_tabled(r), TABLED),
        Pass::Series => (sin_cos_reduced(r), FAST),
    };
    let s = (s, rel * s.hi.abs() + er);
    let c = (c, rel * c.hi.abs() + er);
    Some(quadrant(q, (s, c), |(v, e): (Dd, f64)| (v.neg(), e)))
}

/// The rational tier's sine and cosine of a reduced argument, widened by
/// its error, in `x`'s quadrant.
fn sin_cos_rational(red: &Reduced) -> (Interval, Interval) {
    let (c, s) = certified::cos_sin(&red.r);
    let (s, c) = (widened(&s, red.err), widened(&c, red.err));
    quadrant(red.quadrant, (s, c), |v: Interval| v.neg())
}

/// `sin x` correctly rounded.
pub fn sin(x: f64) -> f64 {
    sin_cos(x).0
}

/// `cos x` correctly rounded.
pub fn cos(x: f64) -> f64 {
    sin_cos(x).1
}

/// `(sin x, cos x)`, each correctly rounded. The last argument's pair is
/// kept per thread: callers taking `sin` and `cos` of one angle one after
/// the other (most of them) pay for one evaluation.
pub fn sin_cos(x: f64) -> (f64, f64) {
    use std::cell::Cell;
    thread_local! {
        static LAST: Cell<(u64, f64, f64)> = const { Cell::new((0, 0.0, 1.0)) };
    }
    let key = x.to_bits();
    if let Ok((k, s, c)) = LAST.try_with(Cell::get) {
        if k == key {
            return (s, c);
        }
    }
    let (s, c) = sin_cos_uncached(x);
    let _ = LAST.try_with(|l| l.set((key, s, c)));
    (s, c)
}

fn sin_cos_uncached(x: f64) -> (f64, f64) {
    if !x.is_finite() {
        return (f64::NAN, f64::NAN);
    }
    // |x| < 2^-27: sin x within a quarter ulp of x, cos x of 1.
    if x.abs() < pow2(-27) {
        return (x, 1.0);
    }
    // The tables, the series, the series after the second reduction:
    // each fills what the ones before left undecided.
    type Pair = Option<((Dd, f64), (Dd, f64))>;
    fn fill(v: Pair, sc: (Option<f64>, Option<f64>)) -> (Option<f64>, Option<f64>) {
        match v {
            Some(((sv, es), (cv, ec))) => (
                sc.0.or_else(|| decide(sv, es)),
                sc.1.or_else(|| decide(cv, ec)),
            ),
            None => sc,
        }
    }
    let mut sc = (None, None);
    for pass in [Pass::Quick, Pass::Tabled, Pass::Series] {
        sc = fill(sin_cos_dd(x, None, pass), sc);
        if sc.0.is_some() && sc.1.is_some() {
            break;
        }
    }
    if let (Some(s), Some(c)) = sc {
        return (s, c);
    }
    let red = reduce_precise(x);
    let fast = sin_cos_dd(x, Some(&red), Pass::Series);
    let (s, c) = fill(fast, sc);
    if let (Some(s), Some(c)) = (s, c) {
        return (s, c);
    }
    let (si, ci) = sin_cos_rational(&red);
    let (sf, cf) = fast.map_or((0.0, 1.0), |((s, _), (c, _))| (s.hi, c.hi));
    let unit = |v: f64| v.clamp(-1.0, 1.0);
    (
        s.unwrap_or_else(|| unit(rational(Some(si), sf))),
        c.unwrap_or_else(|| unit(rational(Some(ci), cf))),
    )
}

/// `tan x` correctly rounded.
pub fn tan(x: f64) -> f64 {
    if !x.is_finite() {
        return f64::NAN;
    }
    // |x| < 2^-27: tan x within a quarter ulp of x.
    if x.abs() < pow2(-27) {
        return x;
    }
    let quotient = |v: Option<((Dd, f64), (Dd, f64))>| {
        let ((s, es), (c, ec)) = v?;
        let t = s.div(c);
        let rel = es / s.hi.abs() + ec / c.hi.abs() + pow2(-100);
        Some((t, rel * t.hi.abs()))
    };
    for pass in [Pass::Quick, Pass::Tabled, Pass::Series] {
        if let Some(v) = quotient(sin_cos_dd(x, None, pass)).and_then(|(t, e)| decide(t, e)) {
            return v;
        }
    }
    let red = reduce_precise(x);
    let fast = quotient(sin_cos_dd(x, Some(&red), Pass::Series));
    if let Some(v) = fast.and_then(|(t, e)| decide(t, e)) {
        return v;
    }
    let (si, ci) = sin_cos_rational(&red);
    rational(si.div(&ci), fast.map_or(0.0, |(t, _)| t.hi))
}

// ------------------------------------------------- arctangent and kin

/// `atan t` for a double-double `0 <= t <= 1`, to `2^-97` relative:
/// `atan(j/16) + atan(u)`, `u = (t - j/16) / (1 + t j/16)`, `|u| <= 1/32`,
/// heads to `u^7` in double-double, tails to `u^21`.
fn atan_unit(t: Dd) -> Dd {
    let j = (t.hi * 16.0).round().clamp(0.0, 16.0);
    let c = j / 16.0;
    let u = if j == 0.0 {
        t
    } else {
        t.add_f(-c).div(Dd::f(1.0).add(t.mul_f(c)))
    };
    let w = u.mul(u);
    let heads = [INV_ODD[0], INV_ODD[1].neg(), INV_ODD[2], INV_ODD[3].neg()];
    let v = u.mul(head(w, &heads, tail(w.hi, &ATAN_TAIL)));
    if j == 0.0 {
        v
    } else {
        ATAN_TABLE[j as usize].add(v)
    }
}

/// `atan(y/x)` for `0 <= y <= x` from the tables: `atan(j/256) + atan u`,
/// `u = (y - c x) / (x + c y)`, `c = j/256`, `|u| <= 1/512`, `atan u = u +
/// d` with `d` (`<= 2^-19 |u|`) in binary64. Within `2^-70` relative.
fn atan_tabled(y: Dd, x: Dd) -> Dd {
    let j = (y.hi / x.hi * 256.0).round().clamp(0.0, 256.0);
    let u = if j == 0.0 {
        y.div(x)
    } else {
        let c = j / 256.0;
        y.sub(x.mul_f(c)).div(x.add(y.mul_f(c)))
    };
    let (uh, ul) = (u.hi, u.lo);
    let w = uh * uh;
    let d = -(uh * w)
        * (1.0 / 3.0 - w * (1.0 / 5.0 - w * (1.0 / 7.0 - w * (1.0 / 9.0 - w * (1.0 / 11.0)))))
        - w * ul;
    let v = u.add_f(d);
    if j == 0.0 {
        v
    } else {
        tables().atan[j as usize].add(v)
    }
}

/// The angle of `(x, y)` for `x, y >= 0`, not both zero, in double-double:
/// from the tables (`tabled`) or the series.
fn angle(y: Dd, x: Dd, tabled: bool) -> Dd {
    let at = |y: Dd, x: Dd| {
        if tabled {
            atan_tabled(y, x)
        } else {
            atan_unit(y.div(x))
        }
    };
    if y.hi <= x.hi {
        at(y, x)
    } else {
        HALF_PI.sub(at(x, y))
    }
}

/// The Fast tier's binary64 for a value built by `f` from the tables
/// (`f(true)`, error `TABLED`), else from the series (`FAST`); the series'
/// value either way.
fn two_pass(f: impl Fn(bool) -> Dd) -> (Option<f64>, Dd) {
    let t = f(true);
    if let Some(v) = decide(t, TABLED * t.hi.abs()) {
        return (Some(v), t);
    }
    let t = f(false);
    (decide(t, FAST * t.hi.abs()), t)
}

/// `atan x` correctly rounded.
pub fn atan(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        return std::f64::consts::FRAC_PI_2.copysign(x);
    }
    // |x| < 2^-28: atan x within a quarter ulp of x.
    if x.abs() < pow2(-28) {
        return x;
    }
    let a = x.abs();
    let (v, t) = two_pass(|tabled| {
        let v = angle(Dd::f(a), Dd::f(1.0), tabled);
        if x < 0.0 {
            v.neg()
        } else {
            v
        }
    });
    v.unwrap_or_else(|| rational(Some(certified::atan(&Interval::from_f64(x))), t.hi))
}

/// `atan2(y, x)` correctly rounded, with C99's special values. A quotient
/// of the smaller magnitude by the larger below `2^-900` takes the IEEE
/// quotient (within an ulp; no kernel construction comes near one).
pub fn atan2(y: f64, x: f64) -> f64 {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI as PI64};
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    if y == 0.0 {
        return if x.is_sign_negative() {
            PI64.copysign(y)
        } else {
            y
        };
    }
    if x == 0.0 {
        return FRAC_PI_2.copysign(y);
    }
    if x.is_infinite() {
        let three_quarters = f64::from_bits(0x4002d97c7f3321d2);
        return match (x > 0.0, y.is_infinite()) {
            (true, false) => 0.0f64.copysign(y),
            (false, false) => PI64.copysign(y),
            (true, true) => FRAC_PI_4.copysign(y),
            (false, true) => three_quarters.copysign(y),
        };
    }
    if y.is_infinite() {
        return FRAC_PI_2.copysign(y);
    }
    let (a, b) = (y.abs(), x.abs());
    if a.min(b) < a.max(b) * pow2(-900) {
        let q = if a < b { a / b } else { b / a };
        let v = match (a < b, x > 0.0) {
            (true, true) => q,
            (true, false) => PI64 - q,
            (false, _) => FRAC_PI_2,
        };
        return v.copysign(y);
    }
    // Both scaled alike (exactly) into range.
    let e = exponent(a.max(b));
    let (v, t) = two_pass(|tabled| {
        let t = angle(Dd::f(ldexp(a, -e)), Dd::f(ldexp(b, -e)), tabled);
        let t = if x < 0.0 { PI.sub(t) } else { t };
        if y < 0.0 {
            t.neg()
        } else {
            t
        }
    });
    v.unwrap_or_else(|| {
        let i = certified::atan2(&Interval::from_f64(y), &Interval::from_f64(x));
        rational(i, t.hi)
    })
}

/// `sqrt(1 - x^2)` in double-double for `|x| < 1`: `(1 - |x|)(1 + |x|)`,
/// each factor exact in double-double.
fn cosine_of(x: f64) -> Dd {
    let a = x.abs();
    let (h, l) = two_sum(1.0, -a);
    let (g, m) = two_sum(1.0, a);
    Dd { hi: h, lo: l }.mul(Dd { hi: g, lo: m }).sqrt()
}

/// The rational tier's `sqrt(1 - x^2)`.
fn cosine_interval(x: f64) -> Interval {
    let q = exact(x);
    Interval::exact(int(1) - &q * &q).sqrt()
}

/// `asin x` correctly rounded.
pub fn asin(x: f64) -> f64 {
    if x.is_nan() || x.abs() > 1.0 {
        return f64::NAN;
    }
    if x.abs() == 1.0 {
        return std::f64::consts::FRAC_PI_2.copysign(x);
    }
    // |x| < 2^-27: asin x within a quarter ulp of x.
    if x.abs() < pow2(-27) {
        return x;
    }
    let c = cosine_of(x);
    let (v, t) = two_pass(|tabled| {
        let t = angle(Dd::f(x.abs()), c, tabled);
        if x < 0.0 {
            t.neg()
        } else {
            t
        }
    });
    v.unwrap_or_else(|| {
        let i = certified::atan2(&Interval::from_f64(x), &cosine_interval(x));
        rational(i, t.hi)
    })
}

/// `acos x` correctly rounded.
pub fn acos(x: f64) -> f64 {
    if x.is_nan() || x.abs() > 1.0 {
        return f64::NAN;
    }
    if x == 1.0 {
        return 0.0;
    }
    if x == -1.0 {
        return std::f64::consts::PI;
    }
    // |x| < 2^-55: acos x within a quarter ulp of pi/2's nearest.
    if x.abs() < pow2(-55) {
        return std::f64::consts::FRAC_PI_2;
    }
    let c = cosine_of(x);
    let (v, t) = two_pass(|tabled| {
        let t = angle(c, Dd::f(x.abs()), tabled);
        if x < 0.0 {
            PI.sub(t)
        } else {
            t
        }
    });
    v.unwrap_or_else(|| {
        let i = certified::atan2(&cosine_interval(x), &Interval::from_f64(x));
        rational(i, t.hi)
    })
}

// ------------------------------------------- exponential and logarithm

/// `exp(x) = 2^k m` with `m` in double-double to `2^-96` relative, for
/// `|x| <= 746`: `r = x - k ln 2`, `|r| <= 0.35`, by Taylor, heads to
/// `r^10` in double-double, tails to `r^22`.
fn exp_fast(x: f64) -> (i32, Dd) {
    let k = (x * std::f64::consts::LOG2_E).round();
    let [l1, l2, l3] = LN2_PARTS;
    let r = if k == 0.0 {
        Dd::f(x)
    } else {
        // x - k l1 is exact (Sterbenz); k l1 and k l2 are exact.
        let (h, l) = two_sum(x - k * l1, -(k * l2));
        Dd { hi: h, lo: l }.add_f(-(k * l3))
    };
    (k as i32, head(r, &INV_FACT[..11], tail(r.hi, &EXP_TAIL)))
}

/// `exp x` correctly rounded (in the subnormal range, the double-double's
/// nearest).
pub fn exp(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x > 710.0 {
        return f64::INFINITY;
    }
    if x < -746.0 {
        return 0.0;
    }
    // |x| < 2^-54: exp x within a quarter ulp of 1.
    if x.abs() < pow2(-54) {
        return 1.0;
    }
    let (k, m) = exp_fast(x);
    // Rounding in the scaled domain is rounding in the normal range, and
    // its overflow to infinity is round to nearest's.
    if k > -1021 {
        if let Some(v) = decide(m, FAST * m.hi) {
            return ldexp(v, k);
        }
        let (kr, i) = exp_rational(x);
        if let Some(v) = round_enclosure(&i) {
            return ldexp(v, kr);
        }
    }
    ldexp(m.hi + m.lo, k)
}

/// `atanh s` for a double-double `|s| <= 0.172`, to `2^-97` relative:
/// `s sum s^2n / (2n+1)`, heads to `s^13` in double-double, tails to
/// `s^41`.
fn atanh_series(s: Dd) -> Dd {
    let w = s.mul(s);
    s.mul(head(w, &INV_ODD, tail(w.hi, &LN_TAIL)))
}

/// `ln y` in double-double to `2^-96` relative, for a positive finite
/// double-double `y` away from 1 or exactly 1 apart: `y = 2^e m`, `m` in
/// `[0.7, 1.42)`, `ln m = 2 atanh s`, `s = (m - 1)/(m + 1)`, `|s| <= 0.172`.
fn ln_dd(y: Dd) -> Dd {
    let mut e = exponent(y.hi);
    let mut m = y.scale(-e);
    if m.hi > std::f64::consts::SQRT_2 {
        m = m.scale(-1);
        e += 1;
    }
    // m - 1 is exact in its high part (Sterbenz).
    let s = m.add_f(-1.0).div(m.add_f(1.0));
    let lm = atanh_series(s).scale(1);
    if e == 0 {
        return lm;
    }
    let ef = f64::from(e);
    let [l1, l2, l3] = LN2_PARTS;
    // e l1 and e l2 are exact (|e| < 2^11).
    let (h, l) = two_sum(ef * l1, ef * l2);
    Dd { hi: h, lo: l }.add_f(ef * l3).add(lm)
}

fn ln_fast(x: f64) -> Dd {
    ln_dd(Dd::f(x))
}

/// `ln x` correctly rounded.
pub fn ln(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x.is_infinite() || x == 1.0 {
        return if x == 1.0 { 0.0 } else { x };
    }
    let v = ln_fast(x);
    decide(v, FAST * v.hi.abs())
        .unwrap_or_else(|| rational(certified::ln(&Interval::from_f64(x)), v.hi))
}

/// `log2 x` correctly rounded (a power of two's exponent exactly).
pub fn log2(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x.is_infinite() {
        return x;
    }
    let e = exponent(x);
    if ldexp(x, -e) == 1.0 {
        return f64::from(e);
    }
    let v = ln_fast(x).div(LN2);
    decide(v, FAST * v.hi.abs()).unwrap_or_else(|| {
        let i = certified::ln(&Interval::from_f64(x)).and_then(|l| l.div(certified::ln2()));
        rational(i, v.hi)
    })
}

/// `asinh x` correctly rounded.
pub fn asinh(x: f64) -> f64 {
    if !x.is_finite() || x == 0.0 {
        return x;
    }
    let a = x.abs();
    // |x| < 2^-27: asinh x within a quarter ulp of x.
    if a < pow2(-27) {
        return x;
    }
    let v = if a > pow2(26) {
        // ln(2a) + 1/(4a^2), the next term below 2^-110 of it.
        ln_fast(a).add(LN2).add_f(0.25 / a / a)
    } else {
        // y - 1 = a + a^2 / (1 + sqrt(1 + a^2)), without cancellation.
        let d = Dd::f(a);
        let a2 = d.mul(d);
        let ym1 = d.add(a2.div(a2.add_f(1.0).sqrt().add_f(1.0)));
        if ym1.hi <= 0.4 {
            atanh_series(ym1.div(ym1.add_f(2.0))).scale(1)
        } else {
            ln_dd(ym1.add_f(1.0))
        }
    };
    decide(v, FAST * v.hi).map_or_else(
        || {
            rational(
                Some(certified::asinh(&Interval::from_f64(x))),
                v.hi.copysign(x),
            )
        },
        |r| r.copysign(x),
    )
}

/// `atanh x` correctly rounded.
pub fn atanh(x: f64) -> f64 {
    if x.is_nan() || x.abs() > 1.0 {
        return f64::NAN;
    }
    if x.abs() == 1.0 {
        return f64::INFINITY.copysign(x);
    }
    let a = x.abs();
    // |x| < 2^-27: atanh x within a quarter ulp of x.
    if a < pow2(-27) {
        return x;
    }
    let v = if a <= 0.17 {
        atanh_series(Dd::f(a))
    } else {
        // ln((1 + a)/(1 - a)) / 2, each factor exact in double-double.
        let (p, pl) = two_sum(1.0, a);
        let (m, ml) = two_sum(1.0, -a);
        ln_dd(Dd { hi: p, lo: pl }.div(Dd { hi: m, lo: ml })).scale(-1)
    };
    decide(v, FAST * v.hi).map_or_else(
        || {
            let q = exact(x);
            let arg = (int(1) + &q) / (int(1) - &q);
            let i =
                certified::ln(&Interval::exact(arg)).map(|l| l.scale(&R::new(1.into(), 2.into())));
            rational(i, v.hi.copysign(x))
        },
        |r| r.copysign(x),
    )
}

/// The rational tier's `cosh` (`even`) or `sinh` of `|x| <= 40`.
fn hyperbolic_rational(a: f64, even: bool, fallback: f64) -> f64 {
    let (k, i) = exp_rational(a);
    let e = i.scale(&pow2_r(i64::from(k)));
    let i = Interval::exact(int(1))
        .div(&e)
        .map(|inv| if even { e.add(&inv) } else { e.sub(&inv) }.scale(&R::new(1.into(), 2.into())));
    rational(i, fallback)
}

/// `cosh` (`even`) or `sinh` of `|x|`.
fn hyperbolic(a: f64, even: bool) -> f64 {
    if a > 711.0 {
        return f64::INFINITY;
    }
    if a < 0.5 {
        // sum x^2n / (2n)! or x sum x^2n / (2n+1)!, heads to n = 5 in
        // double-double, tails to n = 12.
        let x = Dd::f(a);
        let w = x.mul(x);
        let v = if even {
            let heads: [Dd; 6] = std::array::from_fn(|n| INV_FACT[2 * n]);
            head(w, &heads, tail(w.hi, &COSH_TAIL))
        } else {
            let heads: [Dd; 6] = std::array::from_fn(|n| INV_FACT[2 * n + 1]);
            x.mul(head(w, &heads, tail(w.hi, &SINH_TAIL)))
        };
        return decide(v, FAST * v.hi).unwrap_or_else(|| hyperbolic_rational(a, even, v.hi));
    }
    let (k, m) = exp_fast(a);
    if a > 40.0 {
        // e^-a is below 2^-115 of e^a: half of e^a, widened by that.
        let err = FAST * m.hi + pow2(-110) * m.hi;
        if let Some(v) = decide(m, err) {
            return ldexp(v, k - 1);
        }
        let (kr, i) = exp_rational(a);
        let tiny = i.hi() * pow2_r(-110);
        let i = Interval::new(i.lo() - &tiny, i.hi() + &tiny);
        return ldexp(round_enclosure(&i).unwrap_or(m.hi + m.lo), kr - 1);
    }
    let e = m.scale(k);
    let inv = Dd::f(1.0).div(e);
    let v = if even { e.add(inv) } else { e.sub(inv) }.scale(-1);
    // e^a's and e^-a's errors, the difference's cancellation at most 2.2.
    decide(v, 4.0 * FAST * v.hi).unwrap_or_else(|| hyperbolic_rational(a, even, v.hi))
}

/// `sinh x` correctly rounded.
pub fn sinh(x: f64) -> f64 {
    if !x.is_finite() || x == 0.0 {
        return x;
    }
    // |x| < 2^-27: sinh x within a quarter ulp of x.
    if x.abs() < pow2(-27) {
        return x;
    }
    hyperbolic(x.abs(), false).copysign(x)
}

/// `cosh x` correctly rounded.
pub fn cosh(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    // |x| < 2^-27: cosh x within a quarter ulp of 1.
    if x.abs() < pow2(-27) {
        return 1.0;
    }
    hyperbolic(x.abs(), true)
}

// --------------------------------------------------------------- methods

/// `f64`'s elementary functions correctly rounded, by `f64`'s names with
/// `_rn`: `x.sin_rn()` is [`sin`]`(x)`, `y.atan2_rn(x)` is [`atan2`]`(y,
/// x)`, `x.hypot_rn(y)` is [`hypot`]`(x, y)`.
pub trait Rn {
    fn sin_rn(self) -> f64;
    fn cos_rn(self) -> f64;
    fn sin_cos_rn(self) -> (f64, f64);
    fn tan_rn(self) -> f64;
    fn atan_rn(self) -> f64;
    fn atan2_rn(self, x: f64) -> f64;
    fn asin_rn(self) -> f64;
    fn acos_rn(self) -> f64;
    fn hypot_rn(self, y: f64) -> f64;
    fn exp_rn(self) -> f64;
    fn ln_rn(self) -> f64;
    fn log2_rn(self) -> f64;
    fn sinh_rn(self) -> f64;
    fn cosh_rn(self) -> f64;
    fn asinh_rn(self) -> f64;
    fn atanh_rn(self) -> f64;
}

impl Rn for f64 {
    fn sin_rn(self) -> f64 {
        sin(self)
    }
    fn cos_rn(self) -> f64 {
        cos(self)
    }
    fn sin_cos_rn(self) -> (f64, f64) {
        sin_cos(self)
    }
    fn tan_rn(self) -> f64 {
        tan(self)
    }
    fn atan_rn(self) -> f64 {
        atan(self)
    }
    fn atan2_rn(self, x: f64) -> f64 {
        atan2(self, x)
    }
    fn asin_rn(self) -> f64 {
        asin(self)
    }
    fn acos_rn(self) -> f64 {
        acos(self)
    }
    fn hypot_rn(self, y: f64) -> f64 {
        hypot(self, y)
    }
    fn exp_rn(self) -> f64 {
        exp(self)
    }
    fn ln_rn(self) -> f64 {
        ln(self)
    }
    fn log2_rn(self) -> f64 {
        log2(self)
    }
    fn sinh_rn(self) -> f64 {
        sinh(self)
    }
    fn cosh_rn(self) -> f64 {
        cosh(self)
    }
    fn asinh_rn(self) -> f64 {
        asinh(self)
    }
    fn atanh_rn(self) -> f64 {
        atanh(self)
    }
}

#[cfg(test)]
mod tests;
