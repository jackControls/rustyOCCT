//! S9c.2b.1: two cylinders in turned frames (affine models on their stored
//! axes, elliptic in the world) whose axes cross, meeting in a quartic
//! (REVIEW_NOTES.md, S9c.2b refined). Each cylinder's discriminant `D_K`
//! over its angle is a quartic in the half-angle tangent `t` of a chart,
//! its roots (the turning points of `K`'s graph) isolated exactly. `D_K > 0`
//! all round gives two rings over `K`; otherwise every interval of
//! `D_A >= 0` holds a loop, its turning points of both kinds ordered along
//! it and a rational point of `A`'s angle placed between each adjacent pair
//! of different kinds; the pieces are then verified exactly (no root of the
//! carrier's discriminant in a piece's range). A double root, or an
//! extremum of `D_K` within the resolution of zero (two branches that
//! close), is `Degenerate`. A cap's circle meets the other cylinder where a
//! quartic in its own `t` vanishes: roots on the edge's arc are S9c.2b.2's.
use super::graph::between_ccw;
use super::meet::{tangency, CylPair};
use super::model::*;
use super::num::*;
use super::procedural::{other_of, MeetCrv, Other, Quartic, Ruled};
use crate::polynomial::real::{isolate, isolate_with_gcd, AlgebraicRoot, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;

fn limit(what: &'static str) -> Error {
    Error::ComputationLimit(what)
}

// ------------------------------------------------------------ polynomials

/// A polynomial with rational coefficients, ascending powers.
pub(super) type Poly = Vec<R>;

pub(super) fn trim(mut p: Poly) -> Poly {
    while p.last().is_some_and(|c| *c == zero()) {
        p.pop();
    }
    p
}

pub(super) fn padd(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| match (a.get(i), b.get(i)) {
                (Some(x), Some(y)) => crate::rational::add(x, y),
                (Some(x), None) | (None, Some(x)) => x.clone(),
                (None, None) => zero(),
            })
            .collect(),
    )
}

pub(super) fn pmul(a: &Poly, b: &Poly) -> Poly {
    super::num::product(a, b)
}

pub(super) fn pscale(a: &Poly, k: &R) -> Poly {
    trim(a.iter().map(|x| crate::rational::mul(x, k)).collect())
}

pub(super) fn pderiv(a: &Poly) -> Poly {
    trim(
        a.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| crate::rational::mul(c, &int(i as i64)))
            .collect(),
    )
}

/// The remainder of `a` by `b` (`b` nonzero).
fn prem(a: &Poly, b: &Poly) -> Poly {
    let mut r = a.clone();
    let lead = b.last().expect("a nonzero divisor").clone();
    while r.len() >= b.len() && !r.is_empty() {
        let shift = r.len() - b.len();
        let c = crate::rational::div(r.last().expect("nonempty"), &lead);
        for (i, x) in b.iter().enumerate() {
            r[i + shift] = crate::rational::sub(&r[i + shift], &crate::rational::mul(&c, x));
        }
        r = trim(r);
    }
    r
}

/// The Sturm chain of a square-free polynomial.
pub(super) fn sturm(p: &Poly) -> Vec<Poly> {
    let mut chain = vec![p.clone(), pderiv(p)];
    while !chain.last().expect("a chain").is_empty() {
        let n = chain.len();
        let r = prem(&chain[n - 2], &chain[n - 1]);
        if r.is_empty() {
            break;
        }
        chain.push(pscale(&r, &int(-1)));
    }
    chain.retain(|x| !x.is_empty());
    chain
}

/// A polynomial with integer coefficients, ascending powers.
pub(super) type IPoly = Vec<BigInt>;

/// `sturm` in integers: each member a positive multiple of `sturm(p)`'s
/// (positively scaled pseudo-remainders, contents removed), so the sign
/// changes at every point are the same.
pub(super) fn sturm_int(p: &Poly) -> Vec<IPoly> {
    let primitive = |mut v: IPoly| -> IPoly {
        while v.last().is_some_and(|c| c.sign() == Sign::NoSign) {
            v.pop();
        }
        let g = v
            .iter()
            .fold(BigInt::from(0), |g, c| crate::rational::gcd(&g, c));
        if g > BigInt::from(1) {
            for c in &mut v {
                *c /= &g;
            }
        }
        v
    };
    let first = primitive(IntPolynomial::from_rationals(p).0);
    let second = primitive(
        first
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| c * BigInt::from(i))
            .collect(),
    );
    let mut chain = vec![first, second];
    while !chain.last().expect("a chain").is_empty() {
        let n = chain.len();
        let (a, b) = (&chain[n - 2], &chain[n - 1]);
        // `|lead|^k a` reduced by `b`: a positive multiple of the remainder.
        let lead = b.last().expect("a nonzero divisor");
        let (magnitude, negative) = (
            BigInt::from(lead.magnitude().clone()),
            lead.sign() == Sign::Minus,
        );
        let mut r = a.clone();
        while r.len() >= b.len() && !r.is_empty() {
            let shift = r.len() - b.len();
            let top = r.last().expect("nonempty").clone();
            let m = if negative { -top } else { top };
            for c in &mut r {
                *c *= &magnitude;
            }
            for (i, c) in b.iter().enumerate() {
                r[i + shift] -= &m * c;
            }
            while r.last().is_some_and(|c| c.sign() == Sign::NoSign) {
                r.pop();
            }
        }
        if r.is_empty() {
            break;
        }
        let r = primitive(r);
        chain.push(r.into_iter().map(|c| -c).collect());
    }
    chain.retain(|x| !x.is_empty());
    chain
}

/// Sign changes of an integer chain at a rational (zeros skipped): the
/// values homogenized by a positive power of its denominator.
pub(super) fn changes_int(chain: &[IPoly], x: &R) -> usize {
    let (n, d) = (x.numer(), x.denom());
    let signs: Vec<Sign> = chain
        .iter()
        .map(|p| {
            let mut it = p.iter().rev();
            let Some(last) = it.next() else {
                return Sign::NoSign;
            };
            let mut value = last.clone();
            let mut den = BigInt::from(1);
            for c in it {
                den *= d;
                value = value * n + c * &den;
            }
            value.sign()
        })
        .filter(|s| *s != Sign::NoSign)
        .collect();
    signs.windows(2).filter(|w| w[0] != w[1]).count()
}

/// A polynomial's exact value at a surd.
fn eval(p: &Poly, x: &Qd) -> Qd {
    let mut acc = Qd::rat(zero());
    for c in p.iter().rev() {
        acc = acc.mul(x).add_r(c);
    }
    acc
}

/// Sign changes of a Sturm chain at a surd (zeros skipped).
pub(super) fn changes(chain: &[Poly], x: &Qd) -> usize {
    let signs: Vec<Ordering> = chain
        .iter()
        .map(|p| eval(p, x).sign())
        .filter(|s| *s != Ordering::Equal)
        .collect();
    signs.windows(2).filter(|w| w[0] != w[1]).count()
}

fn int_poly(p: &Poly) -> IntPolynomial {
    IntPolynomial::from_rationals(p)
}

/// The real roots of a nonzero polynomial, exactly isolated in increasing
/// order; a repeated root is a tangency.
pub(super) fn roots(p: &Poly) -> Result<Vec<AlgebraicRoot>> {
    let ip = int_poly(p);
    if ip.is_zero() {
        return Err(tangency());
    }
    if ip.is_constant() {
        return Ok(Vec::new());
    }
    // A repeated real root is a tangency (complex ones do not matter: a
    // constant discriminant times (1 + t^2)^2, S9d.2).
    let g = ip.gcd(&ip.derivative());
    if !g.is_constant() {
        let lead = g.0.last().expect("nonzero").clone();
        let bound = g.0[..g.0.len() - 1]
            .iter()
            .map(|c| R::new(c.clone(), lead.clone()))
            .map(|x| if x < zero() { -x } else { x })
            .fold(zero(), |m, x| if x > m { x } else { m })
            + int(1);
        if !isolate(
            &g,
            -bound.clone(),
            bound,
            &mut Budget::new(RootIsolationOptions::default()),
        )?
        .is_empty()
        {
            return Err(tangency());
        }
    }
    let abs = |x: &R| if *x < zero() { -x.clone() } else { x.clone() };
    let lead = abs(p.last().expect("nonzero"));
    let bound = p[..p.len() - 1]
        .iter()
        .map(|c| abs(c) / &lead)
        .fold(zero(), |m, x| if x > m { x } else { m })
        + int(1);
    // `isolate` with the gcd above (it would compute it again).
    isolate_with_gcd(
        &ip,
        g,
        -bound.clone(),
        bound,
        &mut Budget::new(RootIsolationOptions::default()),
    )
}

/// `gcd(p, p')` where `p = (1 + t^2)^m q` (a chart's factor, `m >= 1`) and
/// `q` is certainly coprime with its derivative: then `(1 + t^2)^(m - 1)`,
/// primitive with a positive leading coefficient as `IntPolynomial::gcd`
/// gives it (`p' = (1 + t^2)^(m - 1) (2 m t q + (1 + t^2) q')`, and the
/// irreducible `1 + t^2` divides neither `q` nor `2 m t q`). `None`
/// otherwise: the subresultant chain decides.
fn chart_gcd(p: &IntPolynomial) -> Option<IntPolynomial> {
    use num_bigint::BigInt;
    // Exact division by `t^2 + 1`, where it divides.
    let divide = |p: &[BigInt]| -> Option<Vec<BigInt>> {
        if p.len() < 3 {
            return None;
        }
        let mut r = p.to_vec();
        let mut q = vec![BigInt::from(0); p.len() - 2];
        for k in (0..q.len()).rev() {
            let c = r[k + 2].clone();
            r[k] -= &c;
            r[k + 2] = BigInt::from(0);
            q[k] = c;
        }
        (r[0].sign() == num_bigint::Sign::NoSign && r[1].sign() == num_bigint::Sign::NoSign)
            .then_some(q)
    };
    let mut q = p.0.clone();
    let mut m = 0;
    while let Some(next) = divide(&q) {
        q = next;
        m += 1;
    }
    if m == 0 {
        return None;
    }
    let q = IntPolynomial::new(q);
    if q.is_constant() || !q.coprime_with(&q.derivative()) {
        return None;
    }
    // (1 + t^2)^(m - 1) by its binomial coefficients.
    let k = m - 1;
    let mut out = vec![BigInt::from(0); 2 * k + 1];
    let mut binomial = BigInt::from(1);
    for j in 0..=k {
        out[2 * j] = binomial.clone();
        binomial = binomial * BigInt::from(k - j) / BigInt::from(j + 1);
    }
    Some(IntPolynomial::new(out))
}

/// The distinct real roots of a nonzero, nonconstant polynomial, each with
/// whether it is repeated, and its square-free part (S9d.4b.1: a tangency
/// off a part's rim is no contact).
pub(super) fn roots_repeated(p: &Poly) -> Result<(Poly, Vec<(AlgebraicRoot, bool)>)> {
    let ip = int_poly(p);
    if ip.is_zero() || ip.is_constant() {
        return Err(tangency());
    }
    // `g` and `ip` primitive, so `h = ip / g` is in integers (Gauss's
    // lemma): exactly, and `p / g` is `h` times `p`'s positive multiple of
    // `ip`, each coefficient reduced once (the same as in rationals).
    let ig = chart_gcd(&ip).unwrap_or_else(|| ip.gcd(&ip.derivative()));
    let g = &ig.0;
    let lead = g.last().expect("a nonzero gcd");
    let mut r = ip.0.clone();
    let mut h = vec![BigInt::from(0); r.len().saturating_sub(g.len()) + 1];
    while r.len() >= g.len() && !r.is_empty() {
        let shift = r.len() - g.len();
        let (c, rest) = r.last().expect("nonempty").div_rem(lead);
        debug_assert!(rest.sign() == Sign::NoSign, "an exact quotient");
        for (i, x) in g.iter().enumerate() {
            r[i + shift] -= &c * x;
        }
        h[shift] = c;
        while r.last().is_some_and(|c| c.sign() == Sign::NoSign) {
            r.pop();
        }
    }
    while h.last().is_some_and(|c| c.sign() == Sign::NoSign) {
        h.pop();
    }
    let top = p.iter().rev().find(|c| **c != zero()).expect("nonzero");
    let scale = crate::rational::div(top, &R::from_integer(ip.0.last().expect("nonzero").clone()));
    let sf: Poly = h
        .iter()
        .map(|c| crate::rational::mul(&scale, &R::from_integer(c.clone())))
        .collect();
    let lead = h.last().expect("nonzero").magnitude().clone();
    let bound = h[..h.len() - 1]
        .iter()
        .map(|c| {
            reduced(
                BigInt::from(c.magnitude().clone()),
                &BigInt::from(lead.clone()),
            )
        })
        .fold(zero(), |m, x| if x > m { x } else { m })
        + int(1);
    let rs = isolate(
        &IntPolynomial::new(h),
        -bound.clone(),
        bound,
        &mut Budget::new(RootIsolationOptions::default()),
    )?;
    let out = rs
        .into_iter()
        .map(|r| {
            let repeated = !ig.is_constant() && r.vanishes_polynomial(&ig);
            (r, repeated)
        })
        .collect();
    Ok((sf, out))
}

/// `n / d` in lowest terms (`d` positive), by the crate's gcd.
pub(super) fn reduced(n: BigInt, d: &BigInt) -> R {
    if n.sign() == Sign::NoSign {
        return zero();
    }
    let g = crate::rational::gcd(&n, d);
    if g == BigInt::from(1) {
        R::new_raw(n, d.clone())
    } else {
        R::new_raw(n / &g, d / &g)
    }
}

/// The least common multiple of two positive integers, by the crate's gcd.
pub(super) fn lcm(a: &BigInt, b: &BigInt) -> BigInt {
    a / crate::rational::gcd(a, b) * b
}

/// A root's midpoint, rational.
pub(super) fn middle(r: &AlgebraicRoot) -> R {
    let (a, b) = r.isolator();
    (a + b) / int(2)
}

// ------------------------------------------------------------ charts

/// A chart of the circle: `(cos, sin)` the base `(c0, s0)` turned by
/// `2 atan t`, counter-clockwise with `t`; the base's antipode at infinity.
#[derive(Debug, Clone)]
pub(super) struct Chart {
    pub(super) c0: R,
    pub(super) s0: R,
}

impl Chart {
    /// The numerators of `cos` and `sin` over `1 + t^2`.
    pub(super) fn numerators(&self) -> [Poly; 2] {
        let (c0, s0) = (&self.c0, &self.s0);
        [
            trim(vec![c0.clone(), -(int(2) * s0), -c0.clone()]),
            trim(vec![s0.clone(), int(2) * c0, -s0.clone()]),
        ]
    }

    /// The direction at `t`: over one denominator in integers, reduced once
    /// (the same rationals as each operation reduced in turn).
    pub(super) fn at(&self, t: &R) -> [R; 2] {
        let (n, d) = (t.numer(), t.denom());
        let (nn, dd) = (n * n, d * d);
        // `(cos, sin)` of the turn: `(d^2 - n^2, 2 n d) / (d^2 + n^2)`.
        let (c, s, w) = (&dd - &nn, (n * d) << 1usize, &dd + &nn);
        // The base over its own denominator.
        let lcm = lcm(self.c0.denom(), self.s0.denom());
        let a = self.c0.numer() * (&lcm / self.c0.denom());
        let b = self.s0.numer() * (&lcm / self.s0.denom());
        let den = w * lcm;
        [
            reduced(&a * &c - &b * &s, &den),
            reduced(&b * &c + &a * &s, &den),
        ]
    }

    /// The chart's `t` of a direction (`None` at the antipode).
    pub(super) fn t_of(&self, cs: &[Qd; 2]) -> Option<Qd> {
        let c = cs[0].scale(&self.c0).add(&cs[1].scale(&self.s0));
        let s = cs[1].scale(&self.c0).sub(&cs[0].scale(&self.s0));
        let den = c.add_r(&int(1));
        if den.sign() == Ordering::Equal {
            return None;
        }
        Some(s.mul(&den.recip()?))
    }
}

/// A polynomial form in `(cos, sin)` of degree at most `deg`: its
/// coefficients by exponents. A chart's polynomial is the form times `(1 +
/// t^2)^deg` (a quadratic form's a quartic in `t`; S9d.3b's discriminants
/// over a cone's angle are quartic forms, octics in `t`).
#[derive(Debug, Clone, Default)]
pub(super) struct Form {
    terms: BTreeMap<(u32, u32), R>,
    deg: u32,
}

/// A linear form `l0 + l1 c + l2 s`.
pub(super) type Lin = [R; 3];

impl Form {
    /// A constant of the given degree.
    pub(super) fn constant(k: R, deg: u32) -> Self {
        let mut f = Self {
            terms: BTreeMap::new(),
            deg,
        };
        f.add_const(&k);
        f
    }

    /// A form from its coefficients by exponents `(cos, sin)`, of the given
    /// degree (S9d.4b.2).
    pub(super) fn from_terms(terms: BTreeMap<(u32, u32), R>, deg: u32) -> Self {
        Self {
            terms: terms.into_iter().filter(|(_, x)| *x != zero()).collect(),
            deg,
        }
    }

    /// A linear form, of degree 1.
    pub(super) fn lin(l: &Lin) -> Self {
        let mut terms = BTreeMap::new();
        for (e, x) in [((0, 0), &l[0]), ((1, 0), &l[1]), ((0, 1), &l[2])] {
            if *x != zero() {
                terms.insert(e, x.clone());
            }
        }
        Self { terms, deg: 1 }
    }

    pub(super) fn add_const(&mut self, k: &R) {
        let e = self.terms.entry((0, 0)).or_insert_with(zero);
        *e = crate::rational::add(e, k);
        if *e == zero() {
            self.terms.remove(&(0, 0));
        }
    }

    pub(super) fn scaled(&self, a: &R) -> Self {
        Self {
            terms: self
                .terms
                .iter()
                .filter(|_| *a != zero())
                .map(|(e, x)| (*e, crate::rational::mul(x, a)))
                .collect(),
            deg: self.deg,
        }
    }

    pub(super) fn add(&self, o: &Self) -> Self {
        let mut terms = self.terms.clone();
        for (e, x) in &o.terms {
            let t = terms.entry(*e).or_insert_with(zero);
            *t = crate::rational::add(t, x);
            if *t == zero() {
                terms.remove(e);
            }
        }
        Self {
            terms,
            deg: self.deg.max(o.deg),
        }
    }

    pub(super) fn sub(&self, o: &Self) -> Self {
        self.add(&o.scaled(&int(-1)))
    }

    pub(super) fn mul(&self, o: &Self) -> Self {
        let mut terms: BTreeMap<(u32, u32), R> = BTreeMap::new();
        for (ea, x) in &self.terms {
            for (eb, y) in &o.terms {
                let t = terms.entry((ea.0 + eb.0, ea.1 + eb.1)).or_insert_with(zero);
                *t = crate::rational::add(t, &crate::rational::mul(x, y));
            }
        }
        terms.retain(|_, x| *x != zero());
        Self {
            terms,
            deg: self.deg + o.deg,
        }
    }

    /// The largest coefficient's size.
    pub(super) fn max_coefficient(&self) -> R {
        self.terms
            .values()
            .map(|x| if *x < zero() { -x.clone() } else { x.clone() })
            .fold(zero(), |m, x| if x > m { x } else { m })
    }

    /// The form's constant value, when it has no other term.
    pub(super) fn as_constant(&self) -> Option<R> {
        match self.terms.len() {
            0 => Some(zero()),
            1 => self.terms.get(&(0, 0)).cloned(),
            _ => None,
        }
    }

    pub(super) fn degree(&self) -> u32 {
        self.deg
    }

    pub(super) fn value(&self, cs: &[R; 2]) -> R {
        use crate::rational::{add, mul};
        let pow = |x: &R, n: u32| (0..n).fold(int(1), |acc, _| mul(&acc, x));
        self.terms.iter().fold(zero(), |acc, ((i, j), x)| {
            add(&acc, &mul(&mul(x, &pow(&cs[0], *i)), &pow(&cs[1], *j)))
        })
    }

    /// The form's value at any `(cos, sin)` of a field.
    pub(super) fn value_q(&self, cs: &[Qd; 2]) -> Qd {
        let pow = |x: &Qd, n: u32| (0..n).fold(Qd::rat(int(1)), |acc, _| acc.mul(x));
        self.terms.iter().fold(Qd::rat(zero()), |acc, ((i, j), x)| {
            acc.add(&pow(&cs[0], *i).mul(&pow(&cs[1], *j)).scale(x))
        })
    }

    /// Times `(1 + t^2)^deg` in a chart: a polynomial in `t`. In integers
    /// over one denominator, each coefficient reduced once (the same
    /// polynomial as in rationals reduced at every product).
    pub(super) fn poly(&self, chart: &Chart) -> Poly {
        let mul = |p: &IPoly, q: &IPoly| -> IPoly {
            if p.is_empty() || q.is_empty() {
                return Vec::new();
            }
            let mut out = vec![BigInt::from(0); p.len() + q.len() - 1];
            for (i, x) in p.iter().enumerate() {
                for (j, y) in q.iter().enumerate() {
                    out[i + j] += x * y;
                }
            }
            out
        };
        let pw = |p: &IPoly, n: u32| (0..n).fold(vec![BigInt::from(1)], |acc, _| mul(&acc, p));
        // The chart's numerators over the base's denominator `l`, each
        // power of a numerator over `l` once more: `cn^i sn^j w^rest` over
        // `l^(i + j)`, brought to the largest such power.
        let (c0, s0) = (&chart.c0, &chart.s0);
        let l = lcm(c0.denom(), s0.denom());
        let (a, b) = (
            c0.numer() * (&l / c0.denom()),
            s0.numer() * (&l / s0.denom()),
        );
        let cn: IPoly = vec![a.clone(), -(&b << 1usize), -a.clone()];
        let sn: IPoly = vec![b.clone(), &a << 1usize, -b];
        let w: IPoly = vec![BigInt::from(1), BigInt::from(0), BigInt::from(1)];
        let den = self
            .terms
            .values()
            .fold(BigInt::from(1), |m, x| lcm(&m, x.denom()));
        let top = self
            .terms
            .keys()
            .map(|(i, j)| i + j)
            .fold(self.deg, u32::max);
        let mut out: IPoly = Vec::new();
        for ((i, j), x) in &self.terms {
            let rest = self.deg.saturating_sub(i + j);
            let k = x.numer() * (&den / x.denom()) * l.pow(top - i - j);
            let term = mul(&mul(&pw(&cn, *i), &pw(&sn, *j)), &pw(&w, rest));
            if out.len() < term.len() {
                out.resize(term.len(), BigInt::from(0));
            }
            for (o, t) in out.iter_mut().zip(term) {
                *o += t * &k;
            }
        }
        let total = den * l.pow(top);
        trim(out.into_iter().map(|c| reduced(c, &total)).collect())
    }
}

pub(super) fn square_sum(ls: &[Lin]) -> Form {
    let mut f = Form::constant(zero(), 2);
    for l in ls {
        let x = Form::lin(l);
        f = f.add(&x.mul(&x));
    }
    f.deg = 2;
    f
}

/// A cylinder of an operand: its frame, circle centre and radius.
type Cyl<'a> = (&'a Affine, &'a P2, &'a R);

/// A ruling's quadratic against the other cylinder, `A w^2 + 2 B w + C`:
/// `A`, and the discriminant `B^2 - A C` as a form in the carrier's
/// `(cos, sin)`.
pub(super) fn discriminant(k: Cyl, o: &Other) -> (R, Form) {
    let slope = zero();
    let (a, d) = ruled_discriminant(
        Ruled {
            f: k.0,
            c: k.1,
            r: k.2,
            k: &slope,
        },
        o,
    );
    (a.as_constant().expect("a cylinder's constant A"), d)
}

/// A ruled carrier's ruling `o + r e + w (n + k e)` (`e = cos x + sin y`)
/// against the other quadric: `A` and `D = B^2 - A C` as forms in its
/// `(cos, sin)` (quadratic and quartic on a cone, S9d.3b).
pub(super) fn ruled_discriminant(k: Ruled, o: &Other) -> (Form, Form) {
    let (a, b, cc) = ruled_quadratic(k, o);
    let d = b.mul(&b).sub(&a.mul(&cc));
    (a, d)
}

/// A ruled carrier's ruling against the other quadric: `A`, `B` and `C` of
/// `A w^2 + 2 B w + C` as quadratic forms in its `(cos, sin)`.
pub(super) fn ruled_quadratic(k: Ruled, o: &Other) -> (Form, Form, Form) {
    let Ruled { f, c, r, k: slope } = k;
    let base = f.point(&c[0], &c[1], &zero());
    let at_base = |g: &V, e: &R| -> Lin { [dot(g, &base) - e, r * dot(g, &f.x), r * dot(g, &f.y)] };
    let along = |g: &V| -> Lin { [dot(g, &f.n), slope * dot(g, &f.x), slope * dot(g, &f.y)] };
    let mut a = Form::constant(zero(), 2);
    let mut b = Form::constant(zero(), 2);
    let mut cc = Form::constant(zero(), 2);
    for (g, e) in o.g.iter().zip(&o.e) {
        let (s, n) = (Form::lin(&at_base(g, e)), Form::lin(&along(g)));
        a = a.add(&n.mul(&n));
        b = b.add(&s.mul(&n));
        cc = cc.add(&s.mul(&s));
    }
    // The other's radius term: `r0 + w rd`.
    let r0 = Form::lin(&[
        &o.r + &o.t * (dot(&o.h, &base) - &o.eh),
        &o.t * r * dot(&o.h, &f.x),
        &o.t * r * dot(&o.h, &f.y),
    ]);
    let rd = Form::lin(&along(&o.h).map(|x| &o.t * x));
    a = a.sub(&rd.mul(&rd));
    b = b.sub(&r0.mul(&rd));
    cc = cc.sub(&r0.mul(&r0));
    (a, b, cc)
}

/// A chart whose antipode has `D < 0` (a loop's piece never reaches it),
/// or `None` when `D >= 0` on the four axis points and every root gap.
pub(super) fn negative_chart(d: &Form) -> Result<Option<Chart>> {
    for (c0, s0) in [(1, 0), (0, 1), (-1, 0), (0, -1)] {
        if d.value(&[int(-c0), int(-s0)]) < zero() {
            return Ok(Some(Chart {
                c0: int(c0),
                s0: int(s0),
            }));
        }
    }
    // A negative point between two roots of a chart with a nonzero lead.
    for (c0, s0) in [(1, 0), (0, 1)] {
        let chart = Chart {
            c0: int(c0),
            s0: int(s0),
        };
        if d.value(&[int(-c0), int(-s0)]) == zero() {
            continue;
        }
        let rs = roots(&d.poly(&chart))?;
        for w in rs.windows(2) {
            let t = (w[0].isolator().1 + w[1].isolator().0) / int(2);
            let p = chart.at(&t);
            if d.value(&p) < zero() {
                return Ok(Some(Chart {
                    c0: -p[0].clone(),
                    s0: -p[1].clone(),
                }));
            }
        }
        return Ok(None);
    }
    Err(limit(
        "a cylinders' discriminant vanishing at every axis point",
    ))
}

/// Two branches within the resolution: an extremum of `D` (a critical
/// point of its chart's quartic) where `|D| < (A res / 2)^2`, or a
/// repeated root (a tangency, `roots`).
pub(super) fn near_node(a: &R, d: &Form, chart: &Chart, res: f64) -> Result<()> {
    let p = d.poly(chart);
    let _ = roots(&p)?;
    let delta = {
        let h = a * q(res) / int(2);
        &h * &h
    };
    let w2 = (0..d.degree()).fold(vec![int(1)], |acc, _| {
        pmul(&acc, &vec![int(1), zero(), int(1)])
    });
    let lo = int_poly(&padd(&p, &pscale(&w2, &delta)));
    let hi = int_poly(&padd(&p, &pscale(&w2, &-delta.clone())));
    let dp = pderiv(&p);
    if dp.is_empty() {
        return Ok(());
    }
    let crit = match roots(&dp) {
        Ok(r) => r,
        // A repeated critical point is an inflection, not a node.
        Err(Error::Degenerate(_)) => return Ok(()),
        Err(e) => return Err(e),
    };
    for r in crit {
        if r.sign_polynomial(&lo) == Ordering::Greater && r.sign_polynomial(&hi) == Ordering::Less {
            return Err(Error::Degenerate(
                "two cylinders' section within the resolution of a node",
            ));
        }
    }
    Ok(())
}

/// `near_node`'s test over a chart's closed range `t0..t1` alone (S9d.4b.1:
/// a node off a segment's or wedge's wall does not matter): whether an
/// extremum of `D` there lies within the resolution of zero.
pub(super) fn near_node_within(a: &R, d: &Form, chart: &Chart, res: f64, t0: &R, t1: &R) -> bool {
    let p = d.poly(chart);
    let delta = {
        let h = a * q(res) / int(2);
        &h * &h
    };
    let w2 = (0..d.degree()).fold(vec![int(1)], |acc, _| {
        pmul(&acc, &vec![int(1), zero(), int(1)])
    });
    let lo = int_poly(&padd(&p, &pscale(&w2, &delta)));
    let hi = int_poly(&padd(&p, &pscale(&w2, &-delta.clone())));
    let dp = pderiv(&p);
    if dp.is_empty() {
        return false;
    }
    // A repeated critical point (an inflection) or none decides nothing.
    let Ok(crit) = roots(&dp) else {
        return false;
    };
    crit.into_iter().any(|r| {
        r.compare_rational(t0) != Ordering::Less
            && r.compare_rational(t1) != Ordering::Greater
            && r.sign_polynomial(&lo) == Ordering::Greater
            && r.sign_polynomial(&hi) == Ordering::Less
    })
}

// ------------------------------------------------------------ the pair

/// The meeting of two cylinders in turned frames whose axes cross (`x` of
/// operand 0, `y` of operand 1).
pub(super) fn crossing(x: Cyl, y: Cyl, res: f64) -> Result<CylPair> {
    let cyl = [x, y];
    let others = [other_of(y.0, y.1, y.2), other_of(x.0, x.1, x.2)];
    let disc = [
        discriminant(cyl[0], &others[0]),
        discriminant(cyl[1], &others[1]),
    ];
    let piece = |k: usize, plus: bool, range: Option<[[Qd; 2]; 2]>| {
        MeetCrv::new(k, cyl[k], &others[k], plus, range)
    };
    // Charts with a negative antipode, each discriminant's roots, and the
    // near-node test.
    let mut charts: Vec<Option<Chart>> = Vec::new();
    for (a, d) in &disc {
        let chart = negative_chart(d)?;
        let test = chart.clone().unwrap_or(Chart {
            c0: int(1),
            s0: zero(),
        });
        near_node(a, d, &test, res)?;
        charts.push(chart);
    }
    // D_K >= 0 all round (no negative point): rings over K.
    for (k, chart) in charts.iter().enumerate() {
        if chart.is_none() {
            return Ok(CylPair::Quartic(Box::new(Quartic {
                pieces: vec![piece(k, true, None), piece(k, false, None)],
                switches: Vec::new(),
            })));
        }
    }
    let chart_a = charts[0].clone().expect("a chart");
    let chart_b = charts[1].clone().expect("a chart");
    let pa = disc[0].1.poly(&chart_a);
    let pb = disc[1].1.poly(&chart_b);
    let (mut ra, mut rb) = (roots(&pa)?, roots(&pb)?);
    // Isolators narrowed far below binary64, so their midpoints place the
    // switches.
    for r in ra.iter_mut().chain(rb.iter_mut()) {
        r.refine_for_signs(160);
    }
    if ra.is_empty() || rb.is_empty() {
        // D < 0 all round for one cylinder: no ruling meets the other.
        return Ok(CylPair::Apart);
    }
    if ra.len() % 2 != 0 || rb.len() % 2 != 0 {
        return Err(limit("an odd count of a cylinders' turning points"));
    }
    // Turning points of B's graph: where B's ruling touches A, in A's
    // chart with A's branch (binary64 views; verified below).
    let probe_a = piece(0, true, None);
    let probe_b = piece(1, true, None);
    let mut b_turns: Vec<(f64, bool)> = Vec::new();
    for r in &rb {
        let cs = chart_b.at(&middle(r));
        let (f, c, rr) = cyl[1];
        let base = f.point(&(&c[0] + rr * &cs[0]), &(&c[1] + rr * &cs[1]), &zero());
        // The double root's height: w = -B / A.
        let o = &others[1];
        let n: Vec<R> = (0..o.g.len()).map(|i| dot(&o.g[i], &f.n)).collect();
        let s: Vec<R> = (0..o.g.len())
            .map(|i| dot(&o.g[i], &base) - &o.e[i])
            .collect();
        let a2 = n.iter().fold(zero(), |acc, x| acc + x * x);
        let w = -n.iter().zip(&s).fold(zero(), |acc, (x, y)| acc + x * y) / a2;
        let pt = qv(&add(&base, &scale(&f.n, &w)));
        let place = probe_a.place(&pt);
        let t = chart_a
            .t_of(&place)
            .ok_or(limit("a turning point at a chart's antipode"))?;
        let branch = match probe_a.branch_sign(&pt) {
            Ordering::Greater => true,
            Ordering::Less => false,
            Ordering::Equal => return Err(limit("a turning point of both graphs")),
        };
        b_turns.push((t.to_f64(), branch));
    }
    let mut pieces = Vec::new();
    let mut switches = Vec::new();
    for w in ra.chunks(2) {
        let (a, b) = (&w[0], &w[1]);
        let (af, bf) = (rational_f64(&middle(a)), rational_f64(&middle(b)));
        // Events along the loop: lambda in [0, 2), + branch from a to b,
        // then - branch back.
        let mut events: Vec<(f64, bool)> = vec![(0.0, true), (1.0, true)];
        for &(t, plus) in &b_turns {
            if t <= af || t >= bf {
                continue;
            }
            let l = (t - af) / (bf - af);
            events.push((if plus { l } else { 2.0 - l }, false));
        }
        events.sort_by(|x, y| x.0.total_cmp(&y.0));
        if events.len() < 3 {
            return Err(limit("a loop without turning points of the other graph"));
        }
        let n = events.len();
        for i in 0..n {
            let (l0, l1) = (
                events[i].0,
                events[(i + 1) % n].0 + if i + 1 == n { 2.0 } else { 0.0 },
            );
            if l1 - l0 < 1e-9 {
                return Err(Error::Degenerate("two turning points within rounding"));
            }
        }
        // (lambda) -> (t, branch).
        let at = |l: f64| -> (R, bool) {
            let l = l.rem_euclid(2.0);
            if l < 1.0 {
                (q(af + l * (bf - af)), true)
            } else {
                (q(bf - (l - 1.0) * (bf - af)), false)
            }
        };
        let point = |l: f64| -> Result<(QV, R, bool)> {
            let (t, plus) = at(l);
            if a.compare_rational(&t) != Ordering::Less
                || b.compare_rational(&t) != Ordering::Greater
            {
                return Err(limit("a switch point off its loop"));
            }
            let p = piece(0, plus, None)
                .at(&chart_a.at(&t))
                .ok_or(limit("a switch point off its piece"))?;
            Ok((p, t, plus))
        };
        // Switches between adjacent events of different kinds.
        let mut sw: Vec<(f64, QV, R, bool)> = Vec::new();
        for i in 0..n {
            let (e0, e1) = (events[i], events[(i + 1) % n]);
            if e0.1 == e1.1 {
                continue;
            }
            let l1 = e1.0 + if i + 1 == n { 2.0 } else { 0.0 };
            let l = 0.5 * (e0.0 + l1);
            let (p, t, plus) = point(l)?;
            sw.push((l.rem_euclid(2.0), p, t, plus));
        }
        sw.sort_by(|x, y| x.0.total_cmp(&y.0));
        let m = sw.len();
        if m < 2 {
            return Err(limit("a loop with one switch"));
        }
        for i in 0..m {
            let (s0, s1) = (&sw[i], &sw[(i + 1) % m]);
            let l1 = s1.0 + if i + 1 == m { 2.0 } else { 0.0 };
            // The run's events.
            let inner: Vec<&(f64, bool)> = events
                .iter()
                .filter(|e| {
                    let l = if e.0 < s0.0 { e.0 + 2.0 } else { e.0 };
                    l > s0.0 && l < l1
                })
                .collect();
            let Some(first) = inner.first() else {
                return Err(limit("a run without turning points"));
            };
            if first.1 {
                // A's turning points: a graph over B's angle.
                let lm = 0.5
                    * (s0.0
                        + if first.0 < s0.0 {
                            first.0 + 2.0
                        } else {
                            first.0
                        });
                let (mid, _, _) = point(lm)?;
                let sign = probe_b.branch_sign(&mid);
                if sign == Ordering::Equal
                    || probe_b.branch_sign(&s0.1) != sign
                    || probe_b.branch_sign(&s1.1) != sign
                {
                    return Err(limit("a run over B's angle changing branch"));
                }
                let (p0, p1, pm) = (
                    probe_b.place(&s0.1),
                    probe_b.place(&s1.1),
                    probe_b.place(&mid),
                );
                let range = if between_ccw(&p0, &pm, &p1) {
                    [p0, p1]
                } else {
                    [p1, p0]
                };
                verify(&pb, &chart_b, &range)?;
                pieces.push(piece(1, sign == Ordering::Greater, Some(range)));
            } else {
                // B's turning points: a graph over A's angle, one branch.
                if s0.3 != s1.3 {
                    return Err(limit("a run over A's angle changing branch"));
                }
                let (lo, hi) = if s0.2 < s1.2 {
                    (&s0.2, &s1.2)
                } else {
                    (&s1.2, &s0.2)
                };
                if ra.iter().any(|r| {
                    r.compare_rational(lo) == Ordering::Greater
                        && r.compare_rational(hi) == Ordering::Less
                }) {
                    return Err(limit("a turning point inside a piece"));
                }
                let rng = [lo, hi].map(|t| {
                    let cs = chart_a.at(t);
                    [Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]
                });
                pieces.push(piece(0, s0.3, Some(rng)));
            }
        }
        switches.extend(sw.into_iter().map(|s| s.1));
    }
    Ok(CylPair::Quartic(Box::new(Quartic { pieces, switches })))
}

/// No root of a carrier's discriminant within a counter-clockwise range
/// (the range clear of the chart's antipode): exact Sturm counts at its
/// ends.
fn verify(p: &Poly, chart: &Chart, range: &[[Qd; 2]; 2]) -> Result<()> {
    let (Some(t0), Some(t1)) = (chart.t_of(&range[0]), chart.t_of(&range[1])) else {
        return Err(limit("a piece through a chart's antipode"));
    };
    if t0.cmp(&t1) != Ordering::Less {
        return Err(limit("a piece through a chart's antipode"));
    }
    let chain = sturm(&trim(p.clone()));
    if changes(&chain, &t0) != changes(&chain, &t1) {
        return Err(limit("a turning point inside a piece"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    fn ip(c: &[i64]) -> IntPolynomial {
        IntPolynomial::new(c.iter().map(|&x| BigInt::from(x)).collect())
    }

    fn times(a: &IntPolynomial, b: &IntPolynomial) -> IntPolynomial {
        let mut out = vec![BigInt::from(0); a.0.len() + b.0.len() - 1];
        for (i, x) in a.0.iter().enumerate() {
            for (j, y) in b.0.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        IntPolynomial::new(out)
    }

    /// A chart's factor `(1 + t^2)^m` times a square-free part: the gcd
    /// with the derivative is `(1 + t^2)^(m - 1)`, the subresultant
    /// chain's; a part not square-free goes the exact way.
    #[test]
    fn a_charts_factor_gives_the_gcd_with_the_derivative() {
        let w = ip(&[1, 0, 1]);
        let q = ip(&[-6, 1, 4, -3, 7]);
        for m in 1..5 {
            let p = (0..m).fold(q.clone(), |acc, _| times(&acc, &w));
            assert_eq!(chart_gcd(&p), Some(p.gcd(&p.derivative())));
        }
        let square = times(&times(&q, &ip(&[2, -1])), &ip(&[2, -1]));
        let p = times(&times(&square, &w), &w);
        assert_eq!(chart_gcd(&p), None);
        assert_eq!(chart_gcd(&q), None);
    }

    /// A small linear congruential stream.
    fn stream(mut seed: u64) -> impl FnMut() -> i64 {
        move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 33) as i64
        }
    }

    /// `Form::poly` product by product in rationals.
    fn poly_rational(f: &Form, chart: &Chart) -> Poly {
        let [cn, sn] = chart.numerators();
        let w = vec![int(1), zero(), int(1)];
        let pw = |p: &Poly, n: u32| (0..n).fold(vec![int(1)], |acc, _| pmul(&acc, p));
        let mut out: Poly = Vec::new();
        for ((i, j), x) in &f.terms {
            let rest = f.deg.saturating_sub(i + j);
            let term = pmul(&pmul(&pw(&cn, *i), &pw(&sn, *j)), &pw(&w, rest));
            out = padd(&out, &pscale(&term, x));
        }
        out
    }

    /// `roots_repeated` with its quotient in rationals.
    fn roots_repeated_rational(p: &Poly) -> (Poly, Vec<(AlgebraicRoot, bool)>) {
        let ip = int_poly(p);
        let g = ip.gcd(&ip.derivative());
        let g: Poly = trim(g.0.iter().map(|c| R::from_integer(c.clone())).collect());
        let mut r = p.clone();
        let lead = g.last().expect("a nonzero gcd").clone();
        let mut quot = vec![zero(); r.len().saturating_sub(g.len()) + 1];
        while r.len() >= g.len() && !r.is_empty() {
            let shift = r.len() - g.len();
            let c = r.last().expect("nonempty") / &lead;
            for (i, x) in g.iter().enumerate() {
                r[i + shift] -= &c * x;
            }
            quot[shift] = c;
            r = trim(r);
        }
        let sf = trim(quot);
        let ig = int_poly(&g);
        let abs = |x: &R| if *x < zero() { -x.clone() } else { x.clone() };
        let lead = abs(sf.last().expect("nonzero"));
        let bound = sf[..sf.len() - 1]
            .iter()
            .map(|c| abs(c) / &lead)
            .fold(zero(), |m, x| if x > m { x } else { m })
            + int(1);
        let rs = isolate(
            &int_poly(&sf),
            -bound.clone(),
            bound,
            &mut Budget::new(RootIsolationOptions::default()),
        )
        .expect("isolated");
        let out = rs
            .into_iter()
            .map(|r| {
                let repeated = !ig.is_constant() && r.vanishes_polynomial(&ig);
                (r, repeated)
            })
            .collect();
        (sf, out)
    }

    #[test]
    fn integer_square_free_part_is_the_rational_one() {
        let mut next = stream(0x1405_7b7e_f767_814f);
        for case in 0..60 {
            let mut r = |m: i64| int(next() % (2 * m + 1) - m) / int(1 + next() % 61);
            // Products of linear and quadratic factors, some repeated, over
            // a rational multiple.
            let mut p: Poly = vec![r(9) + int(1) / int(7)];
            for k in 0..(1 + case % 4) {
                let f = if k % 2 == 0 {
                    vec![r(20), int(1)]
                } else {
                    vec![r(20), r(9), r(9) + int(10)]
                };
                p = pmul(&p, &f);
                if (case + k) % 3 == 0 {
                    p = pmul(&p, &f);
                }
            }
            let p = trim(p);
            if p.len() < 2 {
                continue;
            }
            let (sf, rs) = roots_repeated(&p).expect("roots");
            let (sf2, rs2) = roots_repeated_rational(&p);
            assert_eq!(sf, sf2);
            assert_eq!(rs.len(), rs2.len());
            for ((a, x), (b, y)) in rs.iter().zip(&rs2) {
                assert_eq!((a.isolator(), x), (b.isolator(), y));
            }
        }
    }

    #[test]
    fn integer_chart_arithmetic_is_the_rational_one() {
        let mut next = stream(0x5851_f42d_4c95_7f2d);
        for case in 0..60 {
            let mut r = |m: i64| int(next() % (2 * m + 1) - m) / int(1 + next() % 89);
            let base = super::super::model::circle_point(&[zero(), zero()], &int(1), &r(40));
            let chart = Chart {
                c0: base[0].clone(),
                s0: base[1].clone(),
            };
            let deg = 1 + case % 4;
            let mut terms = BTreeMap::new();
            for i in 0..=deg {
                for j in 0..=(deg - i) {
                    if (i + j + case) % 3 != 1 {
                        terms.insert((i, j), r(500));
                    }
                }
            }
            let f = Form::from_terms(terms, deg);
            assert_eq!(f.poly(&chart), poly_rational(&f, &chart));
            let t = r(30);
            let den = int(1) + &t * &t;
            let (c, s) = ((int(1) - &t * &t) / &den, int(2) * &t / &den);
            assert_eq!(
                chart.at(&t),
                [
                    &chart.c0 * &c - &chart.s0 * &s,
                    &chart.s0 * &c + &chart.c0 * &s
                ]
            );
            // The chains' sign changes at rational points.
            let p = f.poly(&chart);
            if p.len() > 1 {
                let (a, b) = (sturm(&p), sturm_int(&p));
                assert_eq!(a.len(), b.len());
                for _ in 0..8 {
                    let x = r(5);
                    assert_eq!(changes(&a, &Qd::rat(x.clone())), changes_int(&b, &x));
                }
            }
        }
    }
}
