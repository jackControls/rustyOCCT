//! S9e.3b: a given result's edge on a meeting of two curved faces (`Meet`,
//! `Rise`, `Toric`) or on a plane's general section of a cone or a torus,
//! met by a face of the other input: three surfaces, two or three of them
//! curved (REVIEW_NOTES.md, S9e.3b refined).
//!
//! Every meeting point is found in one `Q(alpha)`. With a plane among the
//! three, the other two are restricted to the plane's rational affine
//! coordinates `(s, t)` (a quadric's conic, a torus's spiric quartic); with
//! three quadrics, the two others are restricted to a ruled one's rulings,
//! the ruling's angle by the half-angle tangent `x` of a rational chart and
//! the place `y` along it. The eliminant is the two restrictions' resultant
//! in the second parameter (degree 4 for two conics, 8 for a conic and a
//! quartic or for three quadrics), its real roots isolated exactly; at each
//! the second parameter is the fibre's gcd over `Q(alpha)` of degree one, so
//! the point is a rational map of `(alpha, y)`. A fibre holding two points
//! (a plane parallel to a cylinder's rulings, a partner parallel to the
//! carrier) retries the projection: the parameters exchanged, sheared, or
//! another ruled surface as carrier. Each point is verified on all three
//! surfaces exactly, kept where the given curve's own test holds it (its
//! branch, window and range), and refused `Degenerate` where the three
//! gradients are dependent (the partner tangent to the curve) or two
//! meetings lie within the resolution. A torus among three curved surfaces,
//! or two tori in one plane, is refused by cost.
use super::meet::{EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use super::procedural::Other;
use super::turned::{padd, pmul, pscale, roots_repeated, trim, Chart, Poly};
use crate::certified::{Fast, Real};
use crate::polynomial::real::AlgebraicRoot;
use crate::solid::split::zero;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::Arc;
use std::sync::OnceLock;

/// The partner's surface tangent to a given meeting at a point on it.
pub(super) const TANGENT: &str =
    "a given meeting of two faces tangent to a face of the other input (S9e.3b)";
/// Two meetings of one given edge within the resolution.
const NEAR: &str = "two meetings of a given meeting within the resolution (S9e.3b)";
/// A torus among three curved surfaces, or two tori in one plane.
pub(super) const COSTLY: &str =
    "a given result's meeting met by a third surface of degree sixteen \
     or more: a torus among three curved surfaces, or two tori in one plane (S9e.3b)";
const UNSEPARATED: &str = "a given meeting's points not separated by a projection (S9e.3b)";
/// The eliminant's budget: its coefficients' bits (numerator and
/// denominator), past which the fields' arithmetic is refused.
const BUDGET: u64 = 4096;
const COSTLY_FIELD: &str =
    "a given meeting's eliminant past its budget of coefficient bits (S9e.3b)";

/// A ruled surface: `o + (r + k w) (cos x + sin y) + w n` (a cylinder's
/// `k` zero, a cone's its slope).
#[derive(Debug, Clone)]
pub(super) struct Ruling {
    pub(super) o: V,
    pub(super) x: V,
    pub(super) y: V,
    pub(super) n: V,
    pub(super) r: R,
    pub(super) k: R,
}

/// A surface's implicit function in world coordinates, exactly.
#[derive(Debug, Clone)]
enum Imp {
    /// `n . p - d`.
    Plane { n: V, d: R },
    /// A cylinder, sphere or cone (`Other`'s form).
    Quad(Other),
    /// `(|l|^2 + R^2 - r^2)^2 - 4 R^2 (l_u^2 + l_v^2)`, `l` the frame's
    /// local coordinates.
    Torus { f: Box<Affine>, big: R, small: R },
}

// ------------------------------------------------------------ bivariate

/// A polynomial in `y` whose coefficients are polynomials in `x`
/// (ascending in both).
type Bv = Vec<Poly>;

fn bv_trim(mut p: Bv) -> Bv {
    for c in p.iter_mut() {
        *c = trim(std::mem::take(c));
    }
    while p.last().is_some_and(|c| c.is_empty()) {
        p.pop();
    }
    p
}

fn bv_add(a: &Bv, b: &Bv) -> Bv {
    let n = a.len().max(b.len());
    let none = Vec::new();
    bv_trim(
        (0..n)
            .map(|j| padd(a.get(j).unwrap_or(&none), b.get(j).unwrap_or(&none)))
            .collect(),
    )
}

fn bv_scale(a: &Bv, k: &R) -> Bv {
    bv_trim(a.iter().map(|c| pscale(c, k)).collect())
}

fn bv_sub(a: &Bv, b: &Bv) -> Bv {
    bv_add(a, &bv_scale(b, &int(-1)))
}

fn bv_mul(a: &Bv, b: &Bv) -> Bv {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out: Bv = vec![Vec::new(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = padd(&out[i + j], &pmul(x, y));
        }
    }
    bv_trim(out)
}

/// `sum c_i p_i`.
fn bv_dot(c: &V, p: &[Bv; 3]) -> Bv {
    (0..3).fold(Vec::new(), |acc, i| bv_add(&acc, &bv_scale(&p[i], &c[i])))
}

/// The polynomial with the roles of `x` and `y` exchanged.
fn bv_transpose(a: &Bv) -> Bv {
    let n = a.iter().map(|c| c.len()).max().unwrap_or(0);
    bv_trim(
        (0..n)
            .map(|i| {
                a.iter()
                    .map(|c| c.get(i).cloned().unwrap_or_else(zero))
                    .collect()
            })
            .collect(),
    )
}

fn bv_deg_x(a: &Bv) -> usize {
    a.iter()
        .map(|c| c.len().saturating_sub(1))
        .max()
        .unwrap_or(0)
}

fn peval(p: &Poly, x: &R) -> R {
    p.iter().rev().fold(zero(), |acc, c| acc * x + c)
}

fn peval_k(p: &Poly, x: &K) -> K {
    p.iter()
        .rev()
        .fold(K::Rat(zero()), |acc, c| acc.mul(x).add(&K::Rat(c.clone())))
}

/// A bivariate polynomial at `x` (a polynomial in `y` over `Q(alpha)`).
fn bv_at_k(a: &Bv, x: &K) -> Vec<K> {
    a.iter().map(|c| peval_k(c, x)).collect()
}

fn keval(p: &[K], y: &K) -> K {
    p.iter()
        .rev()
        .fold(K::Rat(zero()), |acc, c| acc.mul(y).add(c))
}

// ------------------------------------------------------------ resultants

/// The determinant of a square rational matrix (Gaussian elimination).
fn det(mut m: Vec<Vec<R>>) -> R {
    let n = m.len();
    let mut out = int(1);
    for c in 0..n {
        let Some(p) = (c..n).find(|&r| m[r][c] != zero()) else {
            return zero();
        };
        if p != c {
            m.swap(p, c);
            out = -out;
        }
        let pivot = m[c][c].clone();
        out *= &pivot;
        let (top, rest) = m.split_at_mut(c + 1);
        let row = &top[c];
        for other in rest.iter_mut() {
            if other[c] == zero() {
                continue;
            }
            let f = &other[c] / &pivot;
            for (x, y) in other.iter_mut().zip(row).skip(c) {
                *x -= &f * y;
            }
        }
    }
    out
}

/// The Sylvester resultant of `a` and `b` of formal degrees `da` and `db`
/// (ascending rational coefficients).
fn sylvester(a: &[R], da: usize, b: &[R], db: usize) -> R {
    let n = da + db;
    let coef = |p: &[R], i: usize| p.get(i).cloned().unwrap_or_else(zero);
    let mut m = vec![vec![zero(); n]; n];
    for r in 0..db {
        for i in 0..=da {
            m[r][r + da - i] = coef(a, i);
        }
    }
    for r in 0..da {
        for i in 0..=db {
            m[db + r][r + db - i] = coef(b, i);
        }
    }
    det(m)
}

/// The resultant in `y` of two bivariate polynomials, a polynomial in `x`,
/// by its values at `deg + 1` integers and Newton's interpolation.
fn resultant_y(a: &Bv, b: &Bv) -> Poly {
    let (da, db) = (a.len() - 1, b.len() - 1);
    if da == 0 {
        return (0..db).fold(vec![int(1)], |acc, _| pmul(&acc, &a[0]));
    }
    if db == 0 {
        return (0..da).fold(vec![int(1)], |acc, _| pmul(&acc, &b[0]));
    }
    let deg = da * bv_deg_x(b) + db * bv_deg_x(a);
    let xs: Vec<R> = (0..=deg as i64)
        .map(|j| int(j - (deg as i64) / 2))
        .collect();
    let vals: Vec<R> = xs
        .iter()
        .map(|x| {
            let pa: Vec<R> = a.iter().map(|c| peval(c, x)).collect();
            let pb: Vec<R> = b.iter().map(|c| peval(c, x)).collect();
            sylvester(&pa, da, &pb, db)
        })
        .collect();
    interpolate(&xs, &vals)
}

/// The polynomial through `(x_j, v_j)` (Newton's divided differences).
fn interpolate(xs: &[R], vs: &[R]) -> Poly {
    let n = xs.len();
    let mut c = vs.to_vec();
    for k in 1..n {
        for j in (k..n).rev() {
            c[j] = (&c[j] - &c[j - 1]) / (&xs[j] - &xs[j - k]);
        }
    }
    let mut out: Poly = vec![c[n - 1].clone()];
    for j in (0..n - 1).rev() {
        // out = out (x - x_j) + c_j.
        let shifted = pmul(&out, &vec![-xs[j].clone(), int(1)]);
        out = padd(&shifted, &vec![c[j].clone()]);
    }
    trim(out)
}

/// `p / d` where `d` divides `p` exactly, else `None`.
fn pdiv_exact(p: &Poly, d: &Poly) -> Option<Poly> {
    let d = trim(d.clone());
    if d.len() < 2 {
        return None;
    }
    let mut r = trim(p.clone());
    if r.len() < d.len() {
        return None;
    }
    let lead = d.last().expect("nonzero").clone();
    let mut quot = vec![zero(); r.len() - d.len() + 1];
    while r.len() >= d.len() {
        let shift = r.len() - d.len();
        let c = r.last().expect("nonempty") / &lead;
        for (i, x) in d.iter().enumerate() {
            r[i + shift] -= &c * x;
        }
        quot[shift] = c;
        r = trim(r);
    }
    r.is_empty().then(|| trim(quot))
}

// ------------------------------------------------------------ fibres

fn ktrim(mut p: Vec<K>) -> Vec<K> {
    while p.last().is_some_and(|c| c.sign() == Ordering::Equal) {
        p.pop();
    }
    p
}

/// The gcd over `Q(alpha)` of two polynomials (Euclid's algorithm, each
/// leading coefficient's sign exact at `alpha`).
fn kgcd(a: Vec<K>, b: Vec<K>) -> Result<Vec<K>> {
    let (mut a, mut b) = (ktrim(a), ktrim(b));
    while !b.is_empty() {
        let inv = b
            .last()
            .expect("nonempty")
            .recip()
            .ok_or(Error::ComputationLimit(UNSEPARATED))?;
        let mut r = a;
        while r.len() >= b.len() {
            let c = r.last().expect("nonempty").mul(&inv);
            let shift = r.len() - b.len();
            for (i, x) in b.iter().enumerate() {
                r[i + shift] = r[i + shift].sub(&c.mul(x));
            }
            r.pop();
            r = ktrim(r);
        }
        a = b;
        b = r;
    }
    Ok(a)
}

// ------------------------------------------------------------ surfaces

impl Imp {
    /// The function on a parameterization `p / q` times `q^degree`.
    fn restrict(&self, p: &[Bv; 3], q: &Bv) -> Bv {
        let lin = |g: &V, e: &R| bv_sub(&bv_dot(g, p), &bv_scale(q, e));
        match self {
            Imp::Plane { n, d } => lin(n, d),
            Imp::Quad(o) => {
                let mut acc = Vec::new();
                for (g, e) in o.g.iter().zip(&o.e) {
                    let l = lin(g, e);
                    acc = bv_add(&acc, &bv_mul(&l, &l));
                }
                let rad = bv_add(&bv_scale(q, &o.r), &bv_scale(&lin(&o.h, &o.eh), &o.t));
                bv_sub(&acc, &bv_mul(&rad, &rad))
            }
            Imp::Torus { f, big, small } => {
                let l: Vec<Bv> = (0..3)
                    .map(|k| lin(f.row(k), &dot(f.row(k), &f.o)))
                    .collect();
                let q2 = bv_mul(q, q);
                let s = l
                    .iter()
                    .fold(Vec::new(), |acc, x| bv_add(&acc, &bv_mul(x, x)));
                let t = bv_add(&s, &bv_scale(&q2, &(big * big - small * small)));
                let rho2 = bv_add(&bv_mul(&l[0], &l[0]), &bv_mul(&l[1], &l[1]));
                bv_sub(
                    &bv_mul(&t, &t),
                    &bv_scale(&bv_mul(&q2, &rho2), &(int(4) * big * big)),
                )
            }
        }
    }

    fn value(&self, x: &QV) -> Qd {
        match self {
            Imp::Plane { n, d } => qdot(x, n).add_r(&-d.clone()),
            Imp::Quad(o) => o.value(x),
            Imp::Torus { f, big, small } => {
                let l = f.local_q(x);
                let p2 = l[0].mul(&l[0]).add(&l[1].mul(&l[1]));
                let t = p2.add(&l[2].mul(&l[2])).add_r(&(big * big - small * small));
                t.mul(&t).sub(&p2.scale(&(int(4) * big * big)))
            }
        }
    }

    /// The gradient (a positive multiple).
    fn gradient(&self, x: &QV) -> QV {
        match self {
            Imp::Plane { n, .. } => qv(n),
            Imp::Quad(o) => o.gradient(x),
            Imp::Torus { f, big, small } => {
                let l = f.local_q(x);
                let p2 = l[0].mul(&l[0]).add(&l[1].mul(&l[1]));
                let t = p2.add(&l[2].mul(&l[2])).add_r(&(big * big - small * small));
                let t2 = t.add_r(&-(int(2) * big * big));
                let g = [l[0].mul(&t2), l[1].mul(&t2), l[2].mul(&t)];
                (0..3).fold(qv(&[zero(), zero(), zero()]), |acc, k| {
                    qadd(&acc, &qscale(f.row(k), &g[k]))
                })
            }
        }
    }
}

/// A parameterization `p / q` of a plane or of a ruled surface's rulings,
/// bivariate in `(x, y)`.
struct Param {
    p: [Bv; 3],
    q: Bv,
}

impl Param {
    /// The plane `n . p = d` in affine coordinates along `a` and `b`.
    fn plane(n: &V, d: &R, a: &V, b: &V) -> Self {
        // A point of the plane on the axis of `n`'s largest component.
        let abs = |x: &R| if *x < zero() { -x.clone() } else { x.clone() };
        let k = (0..3)
            .max_by(|&i, &j| abs(&n[i]).cmp(&abs(&n[j])))
            .expect("an axis");
        let mut p0 = [zero(), zero(), zero()];
        p0[k] = d / &n[k];
        let p =
            [0, 1, 2].map(|i| bv_trim(vec![vec![p0[i].clone(), a[i].clone()], vec![b[i].clone()]]));
        Self {
            p,
            q: vec![vec![int(1)]],
        }
    }

    /// A ruled surface's rulings over a chart's half-angle tangent: a
    /// cylinder's circle point `Q o + r (C x + S y)` and `y` along its axis
    /// (`y` the height times `Q`), a cone's apex and `y` its radius along
    /// `C x + S y + Q n / k`, over `Q = 1 + x^2`.
    fn ruled(s: &Ruling, chart: &Chart) -> Self {
        let [c, sn] = chart.numerators();
        let qq: Poly = vec![int(1), zero(), int(1)];
        let p = [0, 1, 2].map(|i| {
            let radial = padd(&pscale(&c, &s.x[i]), &pscale(&sn, &s.y[i]));
            if s.k == zero() {
                let base = padd(&pscale(&qq, &s.o[i]), &pscale(&radial, &s.r));
                bv_trim(vec![base, vec![s.n[i].clone()]])
            } else {
                let apex = &s.o[i] - &s.r / &s.k * &s.n[i];
                let dir = padd(&radial, &pscale(&qq, &(&s.n[i] / &s.k)));
                bv_trim(vec![pscale(&qq, &apex), dir])
            }
        });
        Self { p, q: vec![qq] }
    }

    /// The ruling at a rational `(cos, sin)`, over `y` alone.
    fn ruling_at(s: &Ruling, cs: &[R; 2]) -> Self {
        let p = [0, 1, 2].map(|i| {
            let radial = &cs[0] * &s.x[i] + &cs[1] * &s.y[i];
            if s.k == zero() {
                bv_trim(vec![vec![&s.o[i] + &s.r * &radial], vec![s.n[i].clone()]])
            } else {
                let apex = &s.o[i] - &s.r / &s.k * &s.n[i];
                bv_trim(vec![vec![apex], vec![radial + &s.n[i] / &s.k]])
            }
        });
        Self {
            p,
            q: vec![vec![int(1)]],
        }
    }

    /// A ruled parameterization (linear in the place `y` along the ruling)
    /// in `(x, w)`, `w = y + c x`: each point's ruling and place mixed, so
    /// two points on one ruling, or at one place on two, have different
    /// `w`.
    fn sheared(&self, c: &R) -> Self {
        let p = [0, 1, 2].map(|i| {
            let none = Vec::new();
            let (base, dir) = (
                self.p[i].first().unwrap_or(&none),
                self.p[i].get(1).unwrap_or(&none),
            );
            let shift = pscale(&pmul(&vec![zero(), int(1)], dir), &-c.clone());
            bv_trim(vec![padd(base, &shift), dir.clone()])
        });
        Self {
            p,
            q: self.q.clone(),
        }
    }

    fn transposed(&self) -> Self {
        Self {
            p: [0, 1, 2].map(|i| bv_transpose(&self.p[i])),
            q: bv_transpose(&self.q),
        }
    }
}

/// A box `(lo, hi)` of world coordinates.
type Clip = ([f64; 3], [f64; 3]);

/// A bivariate polynomial's coefficients as binary64 intervals, all scaled
/// by one power of two so the largest is near one (the same roots; `None`
/// where one is out of binary64's range even so).
type Fv = Vec<Vec<Fast>>;

/// The exponent of two that brings a set of rationals' largest magnitude
/// near one.
fn scale_exponent<'a>(c: impl Iterator<Item = &'a R>) -> i64 {
    c.filter(|x| **x != zero())
        .map(|x| x.numer().bits() as i64 - x.denom().bits() as i64)
        .max()
        .unwrap_or(0)
}

/// A rational's binary64 enclosure (`Fast::near_r`), a tiny one's the
/// interval about zero holding it.
fn fast_of(x: &R) -> Option<Fast> {
    Fast::near_r(x).or_else(|| {
        (x.numer().bits() < x.denom().bits())
            .then(|| Fast::exact_f64(-1e-250).union(&Fast::exact_f64(1e-250)))
    })
}

fn scaled(a: &Bv, e: i64) -> Option<Fv> {
    let k = if e >= 0 {
        R::new(1.into(), num_bigint::BigInt::from(1) << e as usize)
    } else {
        R::from_integer(num_bigint::BigInt::from(1) << (-e) as usize)
    };
    a.iter()
        .map(|c| c.iter().map(|x| fast_of(&(x * &k))).collect())
        .collect()
}

/// The binary64 enclosures' data of a projection: the two restrictions,
/// each scaled on its own, and the parameterization's numerators and
/// denominator scaled together.
struct Encl {
    g2: Fv,
    g3: Fv,
    p: [Fv; 3],
    q: Fv,
}

impl Encl {
    fn new(g2: &Bv, g3: &Bv, par: &Param) -> Option<Self> {
        let e2 = scale_exponent(g2.iter().flatten());
        let e3 = scale_exponent(g3.iter().flatten());
        let ep = scale_exponent(
            par.p
                .iter()
                .flatten()
                .flatten()
                .chain(par.q.iter().flatten()),
        );
        Some(Self {
            g2: scaled(g2, e2)?,
            g3: scaled(g3, e3)?,
            p: [
                scaled(&par.p[0], ep)?,
                scaled(&par.p[1], ep)?,
                scaled(&par.p[2], ep)?,
            ],
            q: scaled(&par.q, ep)?,
        })
    }
}

/// A polynomial at a binary64 interval.
fn feval(p: &[Fast], x: &Fast) -> Fast {
    p.iter()
        .rev()
        .fold(Fast::exact_f64(0.0), |acc, c| acc.mul(x).add(c))
}

/// An enclosure of the fibre's common point at a root `x` (an interval),
/// from binary64 intervals: Euclid's algorithm on the two polynomials in
/// `y` down to a linear remainder whose root is `y`, then the point `p / q`
/// there. None where a divisor's leading interval holds zero (the exact
/// path decides).
fn enclose_point(e: &Encl, x: &Fast) -> Option<[Fast; 3]> {
    let at = |g: &Fv| -> Vec<Fast> { g.iter().map(|c| feval(c, x)).collect() };
    let (a, b) = (at(&e.g2), at(&e.g3));
    let (mut hi, mut lo) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    // Euclid's algorithm in intervals down to a linear remainder, each
    // divisor's leading coefficient certainly nonzero.
    while lo.len() > 2 {
        let lead = *lo.last().expect("nonempty");
        if !matches!(lead.sign(), Some(Ordering::Less | Ordering::Greater)) {
            return None;
        }
        while hi.len() >= lo.len() {
            let c = hi.last().expect("nonempty").div(&lead)?;
            let shift = hi.len() - lo.len();
            for (i, l) in lo.iter().enumerate() {
                hi[i + shift] = hi[i + shift].sub(&c.mul(l));
            }
            hi.pop();
        }
        (hi, lo) = (lo, hi);
    }
    if lo.len() < 2 {
        return None;
    }
    let (r0, r1) = (lo[0], lo[1]);
    let y = r0.neg().div(&r1)?;
    let val = |g: &Fv| -> Fast {
        g.iter()
            .rev()
            .fold(Fast::exact_f64(0.0), |acc, c| acc.mul(&y).add(&feval(c, x)))
    };
    let q = val(&e.q);
    Some([
        val(&e.p[0]).div(&q)?,
        val(&e.p[1]).div(&q)?,
        val(&e.p[2]).div(&q)?,
    ])
}

/// Whether an enclosure lies clear of a box (padded): no point there.
fn clear_of(enc: &[Fast; 3], clip: &Clip) -> bool {
    let size = (0..3)
        .map(|i| clip.0[i].abs().max(clip.1[i].abs()))
        .fold(1.0f64, f64::max);
    let pad = 1e-6 * size;
    (0..3).any(|i| {
        let (lo, hi) = enc[i].bounds_f64();
        hi < clip.0[i] - pad || lo > clip.1[i] + pad
    })
}

/// A projection's outcome.
enum Solved {
    Elim(Arc<Elim>),
    /// A fibre held two points or more: another projection.
    Retry,
    /// The two restrictions share a component.
    Shared,
}

/// A projection that separates the points: the two restrictions, the
/// parameterization, the eliminant's square-free part, its real roots, and
/// each root's point once computed (none there: a fibre without a common
/// root).
struct Elim {
    g2: Bv,
    g3: Bv,
    par: Param,
    encl: Option<Encl>,
    sf: Poly,
    roots: Vec<AlgebraicRoot>,
    points: Vec<OnceLock<Option<QV>>>,
}

/// The integers of a rational vector over its denominators' lcm, divided by
/// their gcd (the same polynomial up to a positive factor).
fn primitive(c: &[R]) -> Vec<R> {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let den = c.iter().fold(BigInt::from(1), |l, x| l.lcm(x.denom()));
    let num: Vec<BigInt> = c.iter().map(|x| x.numer() * (&den / x.denom())).collect();
    let g = num.iter().fold(BigInt::from(0), |g, x| g.gcd(x));
    if g == BigInt::from(0) {
        return c.to_vec();
    }
    num.into_iter().map(|x| R::from_integer(x / &g)).collect()
}

/// A bivariate polynomial made primitive (one factor for all coefficients).
fn bv_primitive(a: &Bv) -> Bv {
    let flat: Vec<R> = a.iter().flatten().cloned().collect();
    let p = primitive(&flat);
    let mut k = 0;
    a.iter()
        .map(|c| {
            let out = p[k..k + c.len()].to_vec();
            k += c.len();
            out
        })
        .collect()
}

/// The fibre's common point at a root, exactly: `Err(())` where the fibre
/// holds no single point to take (its gcd zero or of degree two or more).
fn fibre_point(e: &Elim, i: usize) -> Result<std::result::Result<Option<QV>, ()>> {
    let g = Arc::new(Gen::new(e.sf.clone(), e.roots[i].clone()));
    let a = K::generator(&g);
    let (fa, fb) = (ktrim(bv_at_k(&e.g2, &a)), ktrim(bv_at_k(&e.g3, &a)));
    if let Some(found) = fibre_low(e, &a, &fa, &fb) {
        return Ok(found);
    }
    // Higher degrees (a torus's quartic): Euclid's algorithm over `Q(alpha)`.
    let gcd = kgcd(fa, fb)?;
    match gcd.len() {
        1 => return Ok(Ok(None)),
        2 => {}
        _ => return Ok(Err(())),
    }
    let y = gcd[0]
        .neg()
        .mul(&gcd[1].recip().ok_or(Error::ComputationLimit(UNSEPARATED))?);
    let den = keval(&bv_at_k(&e.par.q, &a), &y);
    if den.sign() == Ordering::Equal {
        return Ok(Ok(None));
    }
    let inv = den.recip().ok_or(Error::ComputationLimit(UNSEPARATED))?;
    Ok(Ok(Some([0, 1, 2].map(|k| {
        Qd::of(keval(&bv_at_k(&e.par.p[k], &a), &y).mul(&inv))
    }))))
}

/// `sum p_i n^i m^(d - i)`: a polynomial in `y` at `y = n / m` times
/// `m^d` (`d` at least its degree).
fn homog(p: &[K], n: &K, m: &K, d: usize) -> K {
    let mut np = vec![K::Rat(int(1))];
    let mut mp = vec![K::Rat(int(1))];
    for j in 1..=d {
        np.push(np[j - 1].mul(n));
        mp.push(mp[j - 1].mul(m));
    }
    p.iter().enumerate().fold(K::Rat(zero()), |acc, (i, c)| {
        acc.add(&c.mul(&np[i]).mul(&mp[d - i]))
    })
}

/// `fibre_point` for two polynomials of degree one or two in the second
/// parameter (the restrictions of quadrics): their combination without the
/// leading terms (`b_2 a - a_2 b`, the linear one itself, or the two
/// linear ones' resultant), whose root `-r_0 / r_1` is the common one where
/// the other vanishes there, and the point `p / q` homogenized at `(-r_0,
/// r_1)`: one inverse where Euclid's algorithm over `Q(alpha)` takes four,
/// the same gcd and the same point. `None` for other degrees.
fn fibre_low(e: &Elim, a: &K, fa: &[K], fb: &[K]) -> Option<std::result::Result<Option<QV>, ()>> {
    let (la, lb) = (fa.len(), fb.len());
    if !(2..=3).contains(&la) || !(2..=3).contains(&lb) {
        return None;
    }
    let (r, other): (Vec<K>, &[K]) = match (la, lb) {
        (3, 3) => {
            let r = (0..2)
                .map(|j| fb[2].mul(&fa[j]).sub(&fa[2].mul(&fb[j])))
                .collect();
            (ktrim(r), fb)
        }
        (3, 2) => (fb.to_vec(), fa),
        (2, 3) => (fa.to_vec(), fb),
        _ => {
            let c = fa[1].mul(&fb[0]).sub(&fb[1].mul(&fa[0]));
            if c.sign() != Ordering::Equal {
                return Some(Ok(None));
            }
            (fa.to_vec(), fb)
        }
    };
    match r.len() {
        // The two quadratics proportional: a gcd of degree two.
        0 => return Some(Err(())),
        // A nonzero constant: no common root.
        1 => return Some(Ok(None)),
        _ => {}
    }
    let (n, m) = (r[0].neg(), r[1].clone());
    if homog(other, &n, &m, other.len() - 1).sign() != Ordering::Equal {
        return Some(Ok(None));
    }
    let par: Vec<Vec<K>> = e
        .par
        .p
        .iter()
        .chain(std::iter::once(&e.par.q))
        .map(|c| bv_at_k(c, a))
        .collect();
    let d = par
        .iter()
        .map(|c| c.len().saturating_sub(1))
        .max()
        .unwrap_or(0);
    let den = homog(&par[3], &n, &m, d);
    if den.sign() == Ordering::Equal {
        return Some(Ok(None));
    }
    let inv = den.recip()?;
    Some(Ok(Some(
        [0, 1, 2].map(|k| Qd::of(homog(&par[k], &n, &m, d).mul(&inv))),
    )))
}

/// A root's isolator as a binary64 interval, narrowed by 60 bisections.
fn root_box(root: &AlgebraicRoot) -> Option<Fast> {
    let mut r = root.clone();
    r.refine_for_signs(60);
    let (lo, hi) = r.isolator();
    Some(fast_of(lo)?.union(&fast_of(hi)?))
}

/// The eliminant of `g2 = g3 = 0` on a parameterization and its roots;
/// `Retry` where a fibre holds two points (decided by binary64 intervals
/// where the fibre's remainder is certainly linear, exactly otherwise).
fn solve(g2: &Bv, g3: &Bv, par: Param) -> Result<Solved> {
    if g2.is_empty() || g3.is_empty() {
        return Ok(Solved::Shared);
    }
    let (g2, g3) = (bv_primitive(g2), bv_primitive(g3));
    // The budget, before the resultant: its coefficients' size bound (each
    // restriction's degree in `y` times the other's largest coefficient)
    // past `BUDGET` bits is a limit.
    let bits = |g: &Bv| {
        g.iter()
            .flatten()
            .map(|c| c.numer().bits() + c.denom().bits())
            .max()
            .unwrap_or(0)
    };
    let bound = (g2.len() as u64 - 1) * bits(&g3) + (g3.len() as u64 - 1) * bits(&g2);
    if bound > BUDGET {
        return Err(Error::ComputationLimit(COSTLY_FIELD));
    }
    let mut e = trim(resultant_y(&g2, &g3));
    if e.is_empty() {
        return Ok(Solved::Shared);
    }
    // The chart's `1 + x^2` (no real root) where it divides.
    if par.q.len() == 1 && par.q[0].len() > 1 {
        while let Some(d) = pdiv_exact(&e, &par.q[0]) {
            e = d;
        }
    }
    let (sf, roots, repeated): (Poly, Vec<AlgebraicRoot>, Vec<bool>) = if e.len() < 2 {
        (e, Vec::new(), Vec::new())
    } else {
        let (sf, rs) = roots_repeated(&primitive(&e))?;
        let (roots, repeated) = rs.into_iter().unzip();
        (sf, roots, repeated)
    };
    // The budget: an eliminant whose coefficients pass `BUDGET` bits (the
    // fields' every later product as large: a torus's section in turned
    // frames) is a limit.
    if sf
        .iter()
        .any(|c| c.numer().bits() + c.denom().bits() > BUDGET)
    {
        return Err(Error::ComputationLimit(COSTLY_FIELD));
    }
    let points = roots.iter().map(|_| OnceLock::new()).collect();
    let encl = Encl::new(&g2, &g3, &par);
    let elim = Elim {
        g2,
        g3,
        par,
        encl,
        sf,
        roots,
        points,
    };
    for (i, &rep) in repeated.iter().enumerate().take(elim.roots.len()) {
        // A simple root of the eliminant holds at most one common point in
        // its fibre (the resultant's multiplicity counts them): only a
        // repeated one may need another projection.
        if !rep {
            continue;
        }
        let linear = root_box(&elim.roots[i])
            .zip(elim.encl.as_ref())
            .and_then(|(x, encl)| enclose_point(encl, &x))
            .is_some();
        if linear {
            continue;
        }
        match fibre_point(&elim, i)? {
            Ok(p) => {
                let _ = elim.points[i].set(p);
            }
            Err(()) => return Ok(Solved::Retry),
        }
    }
    Ok(Solved::Elim(Arc::new(elim)))
}

/// The eliminant's points, but those certainly outside `clip` (each
/// computed once).
fn points_of(e: &Elim, clip: Option<&Clip>) -> Result<Vec<QV>> {
    let mut out = Vec::new();
    for i in 0..e.roots.len() {
        if let Some(p) = e.points[i].get() {
            out.extend(p.clone());
            continue;
        }
        // A point certainly outside the edge's and the face's boxes is no
        // meeting: no exact work there.
        if let Some(clip) = clip {
            if root_box(&e.roots[i])
                .zip(e.encl.as_ref())
                .and_then(|(x, encl)| enclose_point(encl, &x))
                .is_some_and(|enc| clear_of(&enc, clip))
            {
                continue;
            }
        }
        let p = fibre_point(e, i)?.map_err(|_| Error::ComputationLimit(UNSEPARATED))?;
        let _ = e.points[i].set(p.clone());
        out.extend(p);
    }
    Ok(out)
}

/// A meeting's eliminant by its three surfaces (`None`: the curve shares a
/// component with the third), kept for the next edge of the same curve
/// (every piece of one meeting against one face).
type Found = Option<Arc<Elim>>;

thread_local! {
    static FOUND: std::cell::RefCell<Vec<(String, Found)>> = const { std::cell::RefCell::new(Vec::new()) };
}

fn cached(imps: &[Imp; 3], f: impl FnOnce() -> Result<Found>) -> Result<Found> {
    let key = format!("{imps:?}");
    if let Some(hit) = FOUND.with(|c| {
        c.borrow()
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| v.clone())
    }) {
        return Ok(hit);
    }
    let v = f()?;
    FOUND.with(|c| {
        let mut c = c.borrow_mut();
        if c.len() >= 32 {
            c.remove(0);
        }
        c.push((key, v.clone()));
    });
    Ok(v)
}

/// The partner's face as a surface, with its rulings where it is ruled.
fn partner(py: &Prism, fy: usize) -> Result<(Imp, Option<Ruling>)> {
    let f = &py.f;
    let ruled = |o: V, r: &R, k: &R| Ruling {
        o,
        x: f.x.clone(),
        y: f.y.clone(),
        n: f.n.clone(),
        r: r.clone(),
        k: k.clone(),
    };
    Ok(match &py.faces[fy].surf {
        Surf::Plane { p, m } => (
            Imp::Plane {
                n: m.clone(),
                d: dot(m, p),
            },
            None,
        ),
        Surf::Cyl { c, r, .. } => (
            Imp::Quad(super::cones::other_face(py, fy)),
            Some(ruled(f.point(&c[0], &c[1], &zero()), r, &zero())),
        ),
        Surf::Cone { b, k } => (
            Imp::Quad(super::cones::other_face(py, fy)),
            Some(ruled(f.o.clone(), b, k)),
        ),
        Surf::Sphere { .. } => (Imp::Quad(super::cones::other_face(py, fy)), None),
        Surf::Torus => {
            let ring = py.ring.as_ref().expect("a torus");
            (
                Imp::Torus {
                    f: Box::new(f.clone()),
                    big: ring.big.clone(),
                    small: ring.small.clone(),
                },
                None,
            )
        }
        Surf::Spline(_) => {
            return Err(Error::OutOfDomain(
                "a spline wall against a curved face in any position (S9f.2)",
            ))
        }
    })
}

/// A given curve's two surfaces, and its carrier's rulings where ruled.
fn curve_surfaces(curve: &Crv) -> Result<(Imp, Imp, Option<Ruling>)> {
    Ok(match curve {
        Crv::Meet(m) => {
            let (ruling, carrier, other) = m.surfaces();
            (Imp::Quad(carrier), Imp::Quad(other), Some(ruling))
        }
        Crv::Rise(c) => {
            let (ruling, carrier, sphere) = c.surfaces();
            (Imp::Quad(carrier), Imp::Quad(sphere), Some(ruling))
        }
        Crv::Toric(c) => {
            let m = &c.m;
            let torus = Imp::Torus {
                f: Box::new(m.f.clone()),
                big: m.big.clone(),
                small: m.small.clone(),
            };
            let other = match &m.other {
                super::torus_curved::Far::Quadric(o) => Imp::Quad((**o).clone()),
                super::torus_curved::Far::Torus { f, big, small } => Imp::Torus {
                    f: f.clone(),
                    big: big.clone(),
                    small: small.clone(),
                },
            };
            (torus, other, None)
        }
        Crv::Cone(c) => {
            let ((n, d), cone) = c.surfaces();
            (Imp::Plane { n, d }, Imp::Quad(cone), None)
        }
        Crv::Torus(c) => {
            let f = &c.f;
            let [a, b, mu, kappa] = &c.plane;
            let n = add(
                &add(&scale(f.row(0), a), &scale(f.row(1), b)),
                &scale(f.row(2), mu),
            );
            let d = dot(&n, &f.o) - kappa;
            (
                Imp::Plane { n, d },
                Imp::Torus {
                    f: Box::new(f.clone()),
                    big: c.big.clone(),
                    small: c.small.clone(),
                },
                None,
            )
        }
        _ => unreachable!("a meeting of two surfaces"),
    })
}

/// Whether the given curve's piece holds a point (its own exact test).
fn holds(curve: &Crv, x: &QV) -> bool {
    match curve {
        Crv::Meet(m) => m.on(x),
        Crv::Rise(c) => c.on(x),
        Crv::Toric(c) => c.on(x),
        Crv::Cone(c) => c.on(x),
        Crv::Torus(c) => c.on(x),
        _ => false,
    }
}

/// Two rational vectors spanning the plane normal to `n`.
fn plane_basis(n: &V) -> [V; 2] {
    let abs = |x: &R| if *x < zero() { -x.clone() } else { x.clone() };
    let k = (0..3)
        .min_by(|&i, &j| abs(&n[i]).cmp(&abs(&n[j])))
        .expect("an axis");
    let mut axis = [zero(), zero(), zero()];
    axis[k] = int(1);
    let a = primitive(&cross(n, &axis));
    let a = [a[0].clone(), a[1].clone(), a[2].clone()];
    let b = primitive(&cross(n, &a));
    [a, [b[0].clone(), b[1].clone(), b[2].clone()]]
}

/// The points of three surfaces, a plane among them: the other two on the
/// plane's coordinates, projected along `s` (sheared where a fibre holds
/// two points).
fn with_plane(plane: (&V, &R), a: &Imp, b: &Imp) -> Result<Found> {
    let [e1, e2] = plane_basis(plane.0);
    let shear = |k: i64| add(&e2, &scale(&e1, &int(k)));
    let half = scale(&e1, &(int(1) / int(2)));
    let variants = [
        (e1.clone(), e2.clone()),
        (e2.clone(), e1.clone()),
        (e1.clone(), shear(1)),
        (e1.clone(), shear(-1)),
        (e1.clone(), shear(2)),
        (e1.clone(), add(&e2, &half)),
    ];
    for (u, v) in &variants {
        let par = Param::plane(plane.0, plane.1, u, v);
        let (ra, rb) = (a.restrict(&par.p, &par.q), b.restrict(&par.p, &par.q));
        match solve(&ra, &rb, par)? {
            Solved::Elim(e) => return Ok(Some(e)),
            Solved::Shared => return Ok(None),
            Solved::Retry => continue,
        }
    }
    Err(Error::ComputationLimit(UNSEPARATED))
}

/// The points of three quadrics along a ruled one's rulings (the curve's
/// carrier, then the partner), over charts whose antipode holds none.
fn ruled(imps: &[Imp; 3], rulings: &[(usize, Ruling)]) -> Result<Found> {
    let bases = [
        [int(1), zero()],
        [R::new(3.into(), 5.into()), R::new(4.into(), 5.into())],
        [R::new((-4).into(), 5.into()), R::new(3.into(), 5.into())],
    ];
    // A projection over each chart whose antipode holds no common point:
    // `None` where a fibre holds two points (the next projection).
    let project =
        |ci: usize, s: &Ruling, make: &dyn Fn(&Chart) -> Param| -> Result<Option<Found>> {
            let [a, b] = match ci {
                0 => [&imps[1], &imps[2]],
                1 => [&imps[0], &imps[2]],
                _ => [&imps[0], &imps[1]],
            };
            for base in &bases {
                let chart = Chart {
                    c0: base[0].clone(),
                    s0: base[1].clone(),
                };
                // The chart's antipode: a common point of the two there (or a
                // common complex one) moves the chart.
                let anti = Param::ruling_at(s, &[-base[0].clone(), -base[1].clone()]);
                let (ra, rb) = (a.restrict(&anti.p, &anti.q), b.restrict(&anti.p, &anti.q));
                if ra.is_empty() || rb.is_empty() || resultant_y(&ra, &rb).is_empty() {
                    continue;
                }
                let par = make(&chart);
                let (ra, rb) = (a.restrict(&par.p, &par.q), b.restrict(&par.p, &par.q));
                return Ok(match solve(&ra, &rb, par)? {
                    Solved::Elim(e) => Some(Some(e)),
                    Solved::Shared => Some(None),
                    Solved::Retry => None,
                });
            }
            Ok(None)
        };
    for (ci, s) in rulings {
        // Over the rulings' angle, then over the place along them.
        for swap in [false, true] {
            let made = project(*ci, s, &|chart| {
                let par = Param::ruled(s, chart);
                if swap {
                    par.transposed()
                } else {
                    par
                }
            })?;
            if let Some(found) = made {
                return Ok(found);
            }
        }
    }
    // Over the place sheared by the angle's chart, `y + c x`: points on one
    // ruling at two places and at one place on two rulings both apart (two
    // cylinders with parallel axes meet in two rulings, which a sphere may
    // cross at equal heights: every point then shares its ruling with one
    // and its height with another).
    for (ci, s) in rulings {
        for c in [int(1), int(-1), int(2), int(1) / int(2)] {
            let made = project(*ci, s, &|chart| {
                Param::ruled(s, chart).sheared(&c).transposed()
            })?;
            if let Some(found) = made {
                return Ok(found);
            }
        }
    }
    Err(Error::ComputationLimit(UNSEPARATED))
}

/// Where a given edge's meeting of two surfaces (or a cone's or a torus's
/// section) meets face `fy` of `py`: the points on the curve's piece, each
/// with its place; a point certainly outside `clip` (the edge's and the
/// face's boxes) is not constructed.
pub(super) fn meet(curve: &Crv, py: &Prism, fy: usize, clip: Option<&Clip>) -> Result<EdgeMeet> {
    let (c1, c2, carrier) = curve_surfaces(curve)?;
    let (c3, partner_ruling) = partner(py, fy)?;
    // A closed meeting all round within the resolution of the partner's
    // plane, and not on it (a coaxial pair's circle in a plane normal to
    // their axis, their frames rounded apart): a sliver, as a tangency.
    if let (Crv::Meet(m), Imp::Plane { n, d }) = (curve, &c3) {
        let f = crate::solid::split::rational_f64;
        let (n, d) = (n.clone().map(|x| f(&x)), f(d));
        let size = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let res = py.tolerance.linear();
        let near = |p: &[f64; 3]| (n[0] * p[0] + n[1] * p[1] + n[2] * p[2] - d).abs() <= res * size;
        if m.range.is_none()
            && size > 0.0
            && m.samples(0.0, std::f64::consts::TAU, 32).iter().all(near)
        {
            return Err(Error::Degenerate(TANGENT));
        }
    }
    let imps = [c1, c2, c3];
    let planes: Vec<usize> = (0..3)
        .filter(|&i| matches!(imps[i], Imp::Plane { .. }))
        .collect();
    let tori = imps
        .iter()
        .filter(|i| matches!(i, Imp::Torus { .. }))
        .count();
    let found = match planes.as_slice() {
        [p] => {
            if tori == 2 {
                return Err(Error::OutOfDomain(COSTLY));
            }
            let Imp::Plane { n, d } = &imps[*p] else {
                unreachable!("a plane")
            };
            let others: Vec<&Imp> = (0..3).filter(|i| i != p).map(|i| &imps[i]).collect();
            cached(&imps, || with_plane((n, d), others[0], others[1]))?
        }
        [] => {
            if tori > 0 {
                return Err(Error::OutOfDomain(COSTLY));
            }
            let mut rulings = Vec::new();
            if let Some(s) = carrier {
                rulings.push((0, s));
            }
            if let Some(s) = partner_ruling {
                rulings.push((2, s));
            }
            cached(&imps, || ruled(&imps, &rulings))?
        }
        _ => {
            return Err(Error::ComputationLimit(
                "a given section met by a plane outside S9e.3a's meetings",
            ))
        }
    };
    let Some(elim) = found else {
        // The curve's surfaces share a component with the partner's: the
        // curve on its surface.
        return Ok(EdgeMeet::Along);
    };
    let points = points_of(&elim, clip)?;
    let mut out: Vec<(Pos, QV)> = Vec::new();
    let mut views: Vec<[f64; 3]> = Vec::new();
    let res = py.tolerance.linear();
    for x in points {
        if imps.iter().any(|s| s.value(&x).sign() != Ordering::Equal) {
            return Err(Error::ComputationLimit(
                "a given meeting's point off its surfaces",
            ));
        }
        if !holds(curve, &x) {
            continue;
        }
        let g: Vec<QV> = imps.iter().map(|s| s.gradient(&x)).collect();
        if qqdot(&g[0], &qcross(&g[1], &g[2])).sign() == Ordering::Equal {
            return Err(Error::Degenerate(TANGENT));
        }
        let view = qv_f64(&x);
        for w in &views {
            let d = (0..3).map(|i| (view[i] - w[i]).powi(2)).sum::<f64>().sqrt();
            if d <= res {
                return Err(Error::Degenerate(NEAR));
            }
        }
        views.push(view);
        out.push((super::graph::place(curve, &x), x));
    }
    Ok(EdgeMeet::Points(out))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: i64, y: i64, z: i64) -> V {
        [int(x), int(y), int(z)]
    }

    /// A cylinder of radius 5 about the `z` axis through `(cx, 0)`.
    fn cylinder(cx: i64) -> (Imp, Ruling) {
        let other = Other {
            g: vec![v(1, 0, 0), v(0, 1, 0)],
            e: vec![int(cx), zero()],
            r: int(5),
            h: v(0, 0, 0),
            eh: zero(),
            t: zero(),
        };
        let ruling = Ruling {
            o: v(cx, 0, 0),
            x: v(1, 0, 0),
            y: v(0, 1, 0),
            n: v(0, 0, 1),
            r: int(5),
            k: zero(),
        };
        (Imp::Quad(other), ruling)
    }

    /// Two cylinders with parallel axes meet in two rulings, `x = 3, y =
    /// +-4`, which a ball about `(3, 0, 0)` of radius 5 crosses at the same
    /// heights `z = +-3`: each of the four points shares its ruling with one
    /// and its height with another, so neither surface's rulings nor their
    /// heights separate them; the sheared projection does.
    #[test]
    fn points_on_two_rulings_at_equal_heights_are_separated() {
        let (carrier, along) = cylinder(0);
        let (partner, across) = cylinder(6);
        let ball = Imp::Quad(Other {
            g: vec![v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)],
            e: vec![int(3), zero(), zero()],
            r: int(5),
            h: v(0, 0, 0),
            eh: zero(),
            t: zero(),
        });
        let imps = [carrier, ball, partner];
        let found = ruled(&imps, &[(0, along), (2, across)])
            .expect("separated")
            .expect("isolated points");
        let points = points_of(&found, None).expect("the points");
        let mut views: Vec<[f64; 3]> = points.iter().map(qv_f64).collect();
        views.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
        assert_eq!(
            views,
            vec![
                [3.0, -4.0, -3.0],
                [3.0, -4.0, 3.0],
                [3.0, 4.0, -3.0],
                [3.0, 4.0, 3.0]
            ]
        );
        for p in &points {
            for s in &imps {
                assert_eq!(s.value(p).sign(), Ordering::Equal);
            }
        }
    }
}
