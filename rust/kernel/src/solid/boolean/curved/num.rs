//! S9c.1's numbers: rational vectors, quadratic surds `a + b sqrt(d)` over
//! the rationals and their exact signs, one surd or two apart, and their
//! enclosures.
use crate::certified::{Fast, Interval as I, Real};
use crate::polynomial::real::{AlgebraicRoot, IntPolynomial};
use crate::rational::{add as radd, div as rdiv, mul as rmul, sub as rsub};
use crate::solid::split::{rational_f64, zero};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::{Arc, Mutex, OnceLock};

pub(super) type V = [R; 3];

pub(super) fn int(n: i64) -> R {
    R::from_integer(BigInt::from(n))
}

pub(super) fn add(a: &V, b: &V) -> V {
    [radd(&a[0], &b[0]), radd(&a[1], &b[1]), radd(&a[2], &b[2])]
}

pub(super) fn sub(a: &V, b: &V) -> V {
    [rsub(&a[0], &b[0]), rsub(&a[1], &b[1]), rsub(&a[2], &b[2])]
}

pub(super) fn scale(a: &V, k: &R) -> V {
    [rmul(&a[0], k), rmul(&a[1], k), rmul(&a[2], k)]
}

pub(super) fn dot(a: &V, b: &V) -> R {
    radd(
        &radd(&rmul(&a[0], &b[0]), &rmul(&a[1], &b[1])),
        &rmul(&a[2], &b[2]),
    )
}

pub(super) fn cross(a: &V, b: &V) -> V {
    let m = |x: &R, y: &R, z: &R, w: &R| rsub(&rmul(x, y), &rmul(z, w));
    [
        m(&a[1], &b[2], &a[2], &b[1]),
        m(&a[2], &b[0], &a[0], &b[2]),
        m(&a[0], &b[1], &a[1], &b[0]),
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

// ------------------------------------------------------------ fields

/// A real algebraic generator `alpha` (S9c.2b.2): a root of `poly`
/// (rational coefficients, ascending), isolated by `root`.
#[derive(Debug)]
pub(super) struct Gen {
    pub(super) poly: Vec<R>,
    pub(super) root: AlgebraicRoot,
    /// The isolator narrowed by a count of bisections, kept by that count:
    /// every enclosure of the field's numbers narrows the same root the
    /// same way, and a sign is exact from any isolator.
    narrowed: Mutex<Vec<(usize, AlgebraicRoot)>>,
    /// Exact signs decided (Sturm-Tarski) by polynomial, kept: a vertex's
    /// coordinates are asked the same questions again and again (S9e.3b:
    /// its surfaces' values on every section through it).
    exact: Mutex<std::collections::HashMap<Vec<R>, Ordering>>,
    /// A factor of `poly` known to vanish at the root (the gcd with every
    /// polynomial found zero there): a polynomial it divides is zero at
    /// the root without another Sturm-Tarski count.
    divisor: Mutex<Option<Vec<R>>>,
    /// `x^j mod poly` for `j` from its degree `n` to `2 n - 2`, numerators
    /// over one common denominator (S9d.4b.2b: a product reduced without a
    /// rational operation per coefficient), computed once.
    powers: OnceLock<(Vec<Vec<BigInt>>, BigInt)>,
}

/// The bisections of the isolator a sign query starts from (those of
/// `Qd::interval`'s enclosure, so the two share it).
const SIGN_STEPS: usize = 96;

impl Gen {
    pub(super) fn new(poly: Vec<R>, root: AlgebraicRoot) -> Self {
        Self {
            poly,
            root,
            narrowed: Mutex::new(Vec::new()),
            exact: Mutex::new(std::collections::HashMap::new()),
            divisor: Mutex::new(None),
            powers: OnceLock::new(),
        }
    }

    /// `x^j mod poly`, `j` in `n..=2 n - 2`, over a common denominator.
    fn powers(&self) -> &(Vec<Vec<BigInt>>, BigInt) {
        self.powers.get_or_init(|| {
            let m = &self.poly;
            let n = m.len() - 1;
            let lead = m[n].clone();
            // x^n = -(m_0 + ... + m_(n-1) x^(n-1)) / m_n.
            let first: Vec<R> = m[..n].iter().map(|c| -rdiv(c, &lead)).collect();
            let mut rows = vec![first.clone()];
            for _ in n + 1..=2 * n - 2 {
                let cur = rows.last().expect("a row");
                let top = cur[n - 1].clone();
                let next: Vec<R> = (0..n)
                    .map(|i| {
                        let shifted = if i == 0 { zero() } else { cur[i - 1].clone() };
                        radd(&shifted, &rmul(&top, &first[i]))
                    })
                    .collect();
                rows.push(next);
            }
            let den = rows
                .iter()
                .flatten()
                .fold(BigInt::from(1), |l, c| lcm(&l, c.denom()));
            let table = rows
                .iter()
                .map(|row| row.iter().map(|c| c.numer() * (&den / c.denom())).collect())
                .collect();
            (table, den)
        })
    }

    /// `p q mod poly` (both reduced), exactly: integer products over the
    /// operands' common denominators, the high terms by `powers`, one
    /// reduction per coefficient (the same value as `pmod(pmul(p, q))`).
    fn mul_mod(&self, p: &[R], q: &[R]) -> Vec<R> {
        if p.is_empty() || q.is_empty() {
            return Vec::new();
        }
        let n = self.poly.len() - 1;
        debug_assert!(p.len() <= n && q.len() <= n, "reduced operands");
        let (pi, dp) = integral(p);
        let (qi, dq) = integral(q);
        let mut prod = vec![BigInt::from(0); pi.len() + qi.len() - 1];
        for (i, x) in pi.iter().enumerate() {
            if x.sign() == num_bigint::Sign::NoSign {
                continue;
            }
            for (j, y) in qi.iter().enumerate() {
                prod[i + j] += x * y;
            }
        }
        let mut den = dp * dq;
        let out = if prod.len() <= n {
            prod
        } else {
            let (table, d) = self.powers();
            let mut out: Vec<BigInt> = prod[..n].iter().map(|c| c * d).collect();
            for (j, c) in prod.iter().enumerate().skip(n) {
                if c.sign() == num_bigint::Sign::NoSign {
                    continue;
                }
                for (o, t) in out.iter_mut().zip(&table[j - n]) {
                    *o += c * t;
                }
            }
            den *= d;
            out
        };
        ptrim(out.into_iter().map(|c| ratio(c, &den)).collect())
    }
}

fn lcm(a: &BigInt, b: &BigInt) -> BigInt {
    let g = crate::rational::gcd(a, b);
    a / g * b
}

/// Integer numerators over the least common denominator.
fn integral(p: &[R]) -> (Vec<BigInt>, BigInt) {
    let den = p.iter().fold(BigInt::from(1), |l, c| lcm(&l, c.denom()));
    (
        p.iter().map(|c| c.numer() * (&den / c.denom())).collect(),
        den,
    )
}

/// `n / d` reduced (`d` positive).
fn ratio(n: BigInt, d: &BigInt) -> R {
    if n.sign() == num_bigint::Sign::NoSign {
        return zero();
    }
    let g = crate::rational::gcd(&n, d);
    R::new_raw(n / &g, d / &g)
}

impl Gen {
    /// The root refined by `steps` bisections (`refine_for_signs`) from
    /// its isolation, computed once.
    fn narrowed(&self, steps: usize) -> AlgebraicRoot {
        let mut kept = self.narrowed.lock().unwrap_or_else(|e| e.into_inner());
        if let Some((_, root)) = kept.iter().find(|(n, _)| *n == steps) {
            return root.clone();
        }
        let mut root = self.root.clone();
        root.refine_for_signs(steps);
        kept.push((steps, root.clone()));
        root
    }

    /// The exact sign of a polynomial at the root (Sturm-Tarski on a
    /// narrowed isolator: the same root, so the same sign), a binary64
    /// enclosure over the isolator first (S9d.4b.2b: certain where it
    /// excludes zero).
    fn sign_of(&self, p: &[R]) -> Ordering {
        let root = self.narrowed(SIGN_STEPS);
        let (lo, hi) = root.isolator();
        // A few ulps' enclosures of the rationals (no exact comparisons).
        let v = Fast::near_r(lo).zip(Fast::near_r(hi)).and_then(|(lo, hi)| {
            let x = lo.union(&hi);
            p.iter().rev().try_fold(Fast::exact_f64(0.0), |acc, c| {
                Some(acc.mul(&x).add(&Fast::near_r(c)?))
            })
        });
        match v.and_then(|v| v.sign()) {
            Some(s @ (Ordering::Less | Ordering::Greater)) => s,
            _ => {
                let divisor = self
                    .divisor
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .clone();
                if let Some(d) = &divisor {
                    if pmod(p, d).is_empty() {
                        return Ordering::Equal;
                    }
                }
                let known = self
                    .exact
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .get(p)
                    .copied();
                known.unwrap_or_else(|| {
                    let ip = IntPolynomial::from_rationals(p);
                    let s = root.sign_polynomial(&ip);
                    self.exact
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .insert(p.to_vec(), s);
                    if s == Ordering::Equal && !ip.is_constant() {
                        // The gcd with the known factor (or the generator's
                        // polynomial) vanishes at the root too.
                        let base = divisor.unwrap_or_else(|| self.poly.clone());
                        let g = IntPolynomial::from_rationals(&base).gcd(&ip);
                        if !g.is_constant() {
                            let g: Vec<R> =
                                g.0.iter().map(|c| R::from_integer(c.clone())).collect();
                            *self.divisor.lock().unwrap_or_else(|e| e.into_inner()) = Some(g);
                        }
                    }
                    s
                })
            }
        }
    }
}

/// An element of `Q` or of one `Q(alpha)`: a polynomial in `alpha`
/// (reduced by its generator's polynomial; a constant is `Rat`).
#[derive(Debug, Clone)]
pub(super) enum K {
    Rat(R),
    Alg(Arc<Gen>, Vec<R>),
}

impl PartialEq for K {
    fn eq(&self, o: &Self) -> bool {
        match (self, o) {
            (K::Rat(a), K::Rat(b)) => a == b,
            (K::Alg(g, p), K::Alg(h, q)) => Arc::ptr_eq(g, h) && p == q,
            _ => false,
        }
    }
}

fn ptrim(mut p: Vec<R>) -> Vec<R> {
    while p.last().is_some_and(|c| *c == zero()) {
        p.pop();
    }
    p
}

/// `a mod m` over the rationals (`m` nonconstant).
fn pmod(a: &[R], m: &[R]) -> Vec<R> {
    let mut r = ptrim(a.to_vec());
    let lead = m.last().expect("a divisor").clone();
    while r.len() >= m.len() {
        let shift = r.len() - m.len();
        let c = rdiv(r.last().expect("nonempty"), &lead);
        for (i, x) in m.iter().enumerate() {
            r[i + shift] = rsub(&r[i + shift], &rmul(&c, x));
        }
        r = ptrim(r);
    }
    r
}

/// `(quotient, remainder)` of `a` by `b` (`b` nonzero).
fn pdivmod(a: &[R], b: &[R]) -> (Vec<R>, Vec<R>) {
    let mut r = ptrim(a.to_vec());
    let lead = b.last().expect("a divisor").clone();
    let mut q = vec![zero(); r.len().saturating_sub(b.len()) + 1];
    while r.len() >= b.len() && !r.is_empty() {
        let shift = r.len() - b.len();
        let c = rdiv(r.last().expect("nonempty"), &lead);
        for (i, x) in b.iter().enumerate() {
            r[i + shift] = rsub(&r[i + shift], &rmul(&c, x));
        }
        q[shift] = c;
        r = ptrim(r);
    }
    (ptrim(q), r)
}

fn pmul(a: &[R], b: &[R]) -> Vec<R> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = radd(&out[i + j], &rmul(x, y));
        }
    }
    ptrim(out)
}

fn psub(a: &[R], b: &[R]) -> Vec<R> {
    let n = a.len().max(b.len());
    ptrim(
        (0..n)
            .map(|i| match (a.get(i), b.get(i)) {
                (Some(x), Some(y)) => rsub(x, y),
                (Some(x), None) => x.clone(),
                (None, Some(y)) => -y,
                (None, None) => zero(),
            })
            .collect(),
    )
}

fn two_fields() -> ! {
    panic!("two algebraic fields in one exact expression (S9c.2b.2)")
}

impl K {
    fn norm(g: &Arc<Gen>, p: Vec<R>) -> Self {
        let p = pmod(&p, &g.poly);
        match p.len() {
            0 => K::Rat(zero()),
            1 => K::Rat(p[0].clone()),
            _ => K::Alg(g.clone(), p),
        }
    }

    /// `alpha` itself.
    pub(super) fn generator(g: &Arc<Gen>) -> Self {
        Self::norm(g, vec![zero(), int(1)])
    }

    pub(super) fn gen(&self) -> Option<&Arc<Gen>> {
        match self {
            K::Rat(_) => None,
            K::Alg(g, _) => Some(g),
        }
    }

    fn poly(&self) -> Vec<R> {
        match self {
            K::Rat(a) => ptrim(vec![a.clone()]),
            K::Alg(_, p) => p.clone(),
        }
    }

    /// The common generator of two elements (`None` inside for `Q`), or
    /// `None` for two different fields.
    fn common<'a>(&'a self, o: &'a Self) -> Option<Option<&'a Arc<Gen>>> {
        match (self.gen(), o.gen()) {
            (None, None) => Some(None),
            (Some(g), None) | (None, Some(g)) => Some(Some(g)),
            (Some(g), Some(h)) => Arc::ptr_eq(g, h).then_some(Some(g)),
        }
    }

    fn lift(&self, o: &Self, f: impl Fn(&[R], &[R]) -> Vec<R>, r: impl Fn(&R, &R) -> R) -> Self {
        if let (K::Rat(a), K::Rat(b)) = (self, o) {
            return K::Rat(r(a, b));
        }
        match self.common(o) {
            None => two_fields(),
            Some(None) => {
                let (K::Rat(a), K::Rat(b)) = (self, o) else {
                    unreachable!("rationals")
                };
                K::Rat(r(a, b))
            }
            Some(Some(g)) => {
                let g = g.clone();
                K::norm(&g, f(&self.poly(), &o.poly()))
            }
        }
    }

    pub(super) fn add(&self, o: &Self) -> Self {
        self.lift(
            o,
            |a, b| psub(a, &b.iter().map(|x| -x).collect::<Vec<_>>()),
            radd,
        )
    }

    pub(super) fn sub(&self, o: &Self) -> Self {
        self.lift(o, psub, rsub)
    }

    pub(super) fn mul(&self, o: &Self) -> Self {
        if let (K::Alg(g, p), K::Alg(h, q)) = (self, o) {
            if Arc::ptr_eq(g, h) {
                return K::norm(g, g.mul_mod(p, q));
            }
        }
        self.lift(o, pmul, rmul)
    }

    pub(super) fn scale(&self, k: &R) -> Self {
        match self {
            K::Rat(a) => K::Rat(rmul(a, k)),
            K::Alg(g, p) => K::norm(g, p.iter().map(|x| rmul(x, k)).collect()),
        }
    }

    pub(super) fn neg(&self) -> Self {
        match self {
            K::Rat(a) => K::Rat(-a),
            K::Alg(g, p) => K::Alg(g.clone(), p.iter().map(|x| -x).collect()),
        }
    }

    pub(super) fn is_zero(&self) -> bool {
        matches!(self, K::Rat(a) if *a == zero())
    }

    /// The exact sign (Sturm-Tarski at the generator).
    pub(super) fn sign(&self) -> Ordering {
        match self {
            K::Rat(a) => sign(a),
            K::Alg(g, p) => g.sign_of(p),
        }
    }

    /// The inverse, or `None` for zero.
    pub(super) fn recip(&self) -> Option<Self> {
        match self {
            K::Rat(a) => (*a != zero()).then(|| K::Rat(int(1) / a)),
            K::Alg(g, p) => {
                let mut m = g.poly.clone();
                loop {
                    // Extended Euclid: s p = gcd (mod m).
                    let (mut r0, mut r1) = (m.clone(), pmod(p, &m));
                    let (mut s0, mut s1): (Vec<R>, Vec<R>) = (Vec::new(), vec![int(1)]);
                    while !r1.is_empty() {
                        let (quo, rem) = pdivmod(&r0, &r1);
                        let s2 = psub(&s0, &pmul(&quo, &s1));
                        (r0, r1) = (r1, rem);
                        (s0, s1) = (s1, s2);
                    }
                    if r0.len() == 1 {
                        let k = int(1) / &r0[0];
                        return Some(K::norm(g, s0.iter().map(|x| x * &k).collect()));
                    }
                    // A common factor: zero at alpha, or alpha is a root of
                    // the cofactor.
                    if g.sign_of(&r0) == Ordering::Equal {
                        return None;
                    }
                    m = pdivmod(&m, &r0).0;
                }
            }
        }
    }

    /// Its binary64 value at the middle of the generator's isolator narrowed
    /// by `steps` bisections, evaluated exactly (a view, not an enclosure).
    fn value_near(&self, steps: usize) -> f64 {
        match self {
            K::Rat(a) => rational_f64(a),
            K::Alg(g, p) => {
                let root = g.narrowed(steps);
                let (lo, hi) = root.isolator();
                let m = (lo + hi) / int(2);
                // Horner in integers over one denominator, reduced once:
                // `sum n_i a^i b^(d - i) / (D b^d)` for `m = a / b` and the
                // coefficients `n_i / D` (the same value).
                let (n, den) = integral(p);
                let (a, b) = (m.numer(), m.denom());
                let mut b_power = BigInt::from(1);
                let mut acc = BigInt::from(0);
                for (i, c) in n.iter().enumerate().rev() {
                    if i + 1 < n.len() {
                        b_power *= b;
                    }
                    acc = acc * a + c * &b_power;
                }
                rational_f64(&R::new(acc, den * b_power))
            }
        }
    }

    /// A binary64 enclosure over the generator's isolator narrowed by
    /// `steps` bisections, from a few ulps' enclosures of its rationals
    /// (`None` where one is out of binary64's comfortable range).
    fn enclose_fast(&self, steps: usize) -> Option<Fast> {
        match self {
            K::Rat(a) => Fast::near_r(a),
            K::Alg(g, p) => {
                let root = g.narrowed(steps);
                let (lo, hi) = root.isolator();
                let x = Fast::near_r(lo)?.union(&Fast::near_r(hi)?);
                p.iter().rev().try_fold(Fast::exact_f64(0.0), |acc, c| {
                    Some(acc.mul(&x).add(&Fast::near_r(c)?))
                })
            }
        }
    }

    /// An enclosure, the generator's isolator narrowed by `steps`
    /// bisections.
    pub(super) fn enclose(&self, steps: usize) -> I {
        match self {
            K::Rat(a) => I::exact(a.clone()),
            K::Alg(g, p) => {
                let root = g.narrowed(steps);
                let (lo, hi) = root.isolator();
                let x = I::new(lo.clone(), hi.clone());
                p.iter().rev().fold(I::exact(zero()), |acc, c| {
                    acc.mul(&x).add(&I::exact(c.clone()))
                })
            }
        }
    }
}

/// The sign of an enclosed expression refined until it excludes zero;
/// equal when still undecided below `1e-40` of its magnitude (two numbers
/// of different fields that close count as equal: S9c.2b.2).
pub(super) fn approx_sign(f: impl Fn(usize) -> I) -> Ordering {
    let tiny = R::new(BigInt::from(1), BigInt::from(10).pow(40));
    for steps in [64, 160, 320, 480] {
        let i = f(steps);
        if let Some(s) = i.sign() {
            return s;
        }
        let width = i.hi() - i.lo();
        let size = i.abs_hi().max(int(1));
        if width < &tiny * size {
            return Ordering::Equal;
        }
    }
    Ordering::Equal
}

// ------------------------------------------------------------ surds

/// `a + b sqrt(d)`, `d >= 0`, `a` and `b` in `Q` or one `Q(alpha)`; a
/// number of the base field has `b = d = 0`. A surd whose `d` is a
/// rational square is folded into its base part.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Qd {
    pub(super) a: K,
    pub(super) b: K,
    pub(super) d: R,
}

impl Qd {
    pub(super) fn rat(a: R) -> Self {
        Self::of(K::Rat(a))
    }

    /// A number of the base field.
    pub(super) fn of(a: K) -> Self {
        Self {
            a,
            b: K::Rat(zero()),
            d: zero(),
        }
    }

    pub(super) fn new(a: R, b: R, d: R) -> Self {
        Self::parts(K::Rat(a), K::Rat(b), d)
    }

    pub(super) fn parts(a: K, b: K, d: R) -> Self {
        debug_assert!(d >= zero());
        if b.is_zero() || d == zero() {
            return Self::of(a);
        }
        if let Some(s) = rational_sqrt(&d) {
            return Self::of(a.add(&b.scale(&s)));
        }
        Self { a, b, d }
    }

    /// `parts` for a `d` of an existing number's (zero, or no rational
    /// square: its surd was folded when made), whose square root need not
    /// be tried again.
    fn known(a: K, b: K, d: R) -> Self {
        if b.is_zero() || d == zero() {
            return Self::of(a);
        }
        Self { a, b, d }
    }

    /// A rational, when it is one.
    pub(super) fn rational(&self) -> Option<&R> {
        match (&self.a, self.b.is_zero()) {
            (K::Rat(a), true) => Some(a),
            _ => None,
        }
    }

    pub(super) fn is_rational(&self) -> bool {
        self.rational().is_some()
    }

    /// The surd's `d`, when there is one.
    pub(super) fn field(&self) -> Option<&R> {
        (!self.b.is_zero()).then_some(&self.d)
    }

    /// The base field's generator, when algebraic.
    pub(super) fn gen(&self) -> Option<&Arc<Gen>> {
        self.a.gen().or(self.b.gen())
    }

    /// Whether two numbers' base fields agree.
    fn same_base(&self, o: &Self) -> bool {
        match (self.gen(), o.gen()) {
            (Some(g), Some(h)) => Arc::ptr_eq(g, h),
            _ => true,
        }
    }

    /// The common surd of two numbers of one base field, or `None` when
    /// both are surds of different `d`.
    fn common(&self, o: &Self) -> Option<R> {
        match (self.field(), o.field()) {
            (None, None) => Some(zero()),
            (Some(d), None) | (None, Some(d)) => Some(d.clone()),
            (Some(d), Some(e)) => (d == e).then(|| d.clone()),
        }
    }

    pub(super) fn add(&self, o: &Self) -> Self {
        let d = self.common(o).expect("surds of one field");
        Self::known(self.a.add(&o.a), self.b.add(&o.b), d)
    }

    pub(super) fn sub(&self, o: &Self) -> Self {
        let d = self.common(o).expect("surds of one field");
        Self::known(self.a.sub(&o.a), self.b.sub(&o.b), d)
    }

    pub(super) fn mul(&self, o: &Self) -> Self {
        let d = self.common(o).expect("surds of one field");
        Self::known(
            self.a.mul(&o.a).add(&self.b.mul(&o.b).scale(&d)),
            self.a.mul(&o.b).add(&self.b.mul(&o.a)),
            d,
        )
    }

    pub(super) fn scale(&self, k: &R) -> Self {
        Self::known(self.a.scale(k), self.b.scale(k), self.d.clone())
    }

    /// Times an element of the base field.
    pub(super) fn scale_k(&self, k: &K) -> Self {
        Self::known(self.a.mul(k), self.b.mul(k), self.d.clone())
    }

    pub(super) fn neg(&self) -> Self {
        Self::known(self.a.neg(), self.b.neg(), self.d.clone())
    }

    pub(super) fn add_r(&self, k: &R) -> Self {
        Self::known(
            self.a.add(&K::Rat(k.clone())),
            self.b.clone(),
            self.d.clone(),
        )
    }

    /// The inverse (`None` for zero): `(a - b sqrt d) / (a^2 - b^2 d)`.
    pub(super) fn recip(&self) -> Option<Self> {
        let norm = self.a.mul(&self.a).sub(&self.b.mul(&self.b).scale(&self.d));
        let inv = norm.recip()?;
        Some(Self::known(
            self.a.mul(&inv),
            self.b.neg().mul(&inv),
            self.d.clone(),
        ))
    }

    /// The exact sign (a binary64 enclosure first where the number is a
    /// surd over an algebraic field, S9d.4c: certain where it excludes
    /// zero; the exact `a^2 - b^2 d` a product in the field).
    pub(super) fn sign(&self) -> Ordering {
        if self.d != zero() && !self.b.is_zero() && self.gen().is_some() {
            if let Some(s @ (Ordering::Less | Ordering::Greater)) =
                self.enclose_fast(SIGN_STEPS).and_then(|v| v.sign())
            {
                return s;
            }
        }
        let sa = self.a.sign();
        let sb = if self.d == zero() {
            Ordering::Equal
        } else {
            self.b.sign()
        };
        if sb == Ordering::Equal {
            return sa;
        }
        if sa == Ordering::Equal || sa == sb {
            return sb;
        }
        // Opposite signs: the larger magnitude wins.
        let big = self.a.mul(&self.a).sub(&self.b.mul(&self.b).scale(&self.d));
        match big.sign() {
            Ordering::Greater => sa,
            Ordering::Less => sb,
            Ordering::Equal => Ordering::Equal,
        }
    }

    /// The exact order of two numbers of any fields (enclosures for two
    /// different base fields).
    pub(super) fn cmp(&self, o: &Self) -> Ordering {
        if !self.same_base(o) {
            // Binary64 enclosures first: certain where they exclude zero.
            let fast = self
                .enclose_fast(64)
                .zip(o.enclose_fast(64))
                .and_then(|(x, y)| x.sub(&y).sign());
            if let Some(s @ (Ordering::Less | Ordering::Greater)) = fast {
                return s;
            }
            return approx_sign(|n| self.enclose(n).sub(&o.enclose(n)));
        }
        if let Some(d) = self.common(o) {
            return Self::known(self.a.sub(&o.a), self.b.sub(&o.b), d).sign();
        }
        // (a - a' + b sqrt d) + (-b') sqrt d'.
        tower_sign(
            &Self::known(self.a.sub(&o.a), self.b.clone(), self.d.clone()),
            &Self::of(o.b.neg()),
            &o.d,
        )
    }

    /// `enclose` in binary64 (`K::enclose_fast`).
    fn enclose_fast(&self, n: usize) -> Option<Fast> {
        let a = self.a.enclose_fast(n)?;
        if self.b.is_zero() {
            return Some(a);
        }
        Some(a.add(&self.b.enclose_fast(n)?.mul(&Fast::near_r(&self.d)?.sqrt())))
    }

    fn enclose(&self, n: usize) -> I {
        let a = self.a.enclose(n);
        if self.b.is_zero() {
            return a;
        }
        a.add(&self.b.enclose(n).mul(&I::exact(self.d.clone()).sqrt()))
    }

    pub(super) fn interval(&self) -> I {
        self.enclose(96)
    }

    /// The value rounded to binary64 (from a tight enclosure: narrowed
    /// further where a field's large coefficients leave the first one wide,
    /// S9d.4b.2).
    pub(super) fn to_f64(&self) -> f64 {
        if let Some(a) = self.rational() {
            return rational_f64(a);
        }
        let tiny = R::new(BigInt::from(1), BigInt::from(1u64 << 62));
        let i = self.interval();
        let size = i.abs_hi().max(int(1));
        if i.hi() - i.lo() <= &tiny * size {
            return rational_f64(&((i.lo() + i.hi()) / int(2)));
        }
        // Wide: the value at the narrowed isolator's middle, exactly (a
        // rational interval's products of large coefficients are slow), at
        // more bisections until two agree.
        let mut last = None;
        for steps in [160, 320, 640, 1280] {
            let a = self.a.value_near(steps);
            let b = if self.b.is_zero() {
                0.0
            } else {
                self.b.value_near(steps) * rational_f64(&self.d).sqrt()
            };
            let v = a + b;
            if last == Some(v) {
                return v;
            }
            last = Some(v);
        }
        last.unwrap_or(f64::NAN)
    }
}

/// The exact sign of `x + y sqrt(e)`, `x` and `y` of one field, `e >= 0`.
pub(super) fn tower_sign(x: &Qd, y: &Qd, e: &R) -> Ordering {
    // A binary64 enclosure first (S9d.4c): certain where it excludes zero.
    if *e != zero() {
        let fast = x
            .enclose_fast(SIGN_STEPS)
            .zip(y.enclose_fast(SIGN_STEPS))
            .and_then(|(a, b)| Some(a.add(&b.mul(&Fast::near_r(e)?.sqrt()))));
        if let Some(s @ (Ordering::Less | Ordering::Greater)) = fast.and_then(|v| v.sign()) {
            return s;
        }
    }
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
/// `b_i = p_i + q_i sqrt e`; by enclosures when their base fields are
/// different algebraic ones.
pub(super) fn mixed_dot_sign(a: &[Qd], b: &[Qd]) -> Ordering {
    let ga = a.iter().find_map(|x| x.gen());
    let gb = b.iter().find_map(|x| x.gen());
    // Binary64 enclosures first where a field is algebraic: certain where
    // they exclude zero (S9d.4c: before the exact products in one field
    // too, a tower's the costliest).
    if ga.is_some() || gb.is_some() {
        let steps = if ga.zip(gb).is_some_and(|(g, h)| !Arc::ptr_eq(g, h)) {
            64
        } else {
            SIGN_STEPS
        };
        let fast = a
            .iter()
            .zip(b)
            .try_fold(Fast::exact_f64(0.0), |acc, (x, y)| {
                Some(acc.add(&x.enclose_fast(steps)?.mul(&y.enclose_fast(steps)?)))
            });
        if let Some(s @ (Ordering::Less | Ordering::Greater)) = fast.and_then(|v| v.sign()) {
            return s;
        }
    }
    if let (Some(g), Some(h)) = (ga, gb) {
        if !Arc::ptr_eq(g, h) {
            return approx_sign(|n| {
                a.iter().zip(b).fold(I::exact(zero()), |acc, (x, y)| {
                    acc.add(&x.enclose(n).mul(&y.enclose(n)))
                })
            });
        }
    }
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
        x = x.add(&ai.scale_k(&bi.a));
        y = y.add(&ai.scale_k(&bi.b));
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

    /// Products reduced through the powers of the generator (S9d.4b.2b)
    /// are the rational reduction's, and signs from the binary64
    /// enclosure are the exact ones.
    #[test]
    fn products_in_a_field_reduce_as_rationals() {
        // A root of 3 x^5 - 7/2 x^3 + x - 5/9 near 1.
        let poly = vec![r(-5, 9), int(1), zero(), r(-7, 2), zero(), int(3)];
        let ip = IntPolynomial::from_rationals(&poly);
        let roots = crate::polynomial::real::isolate(
            &ip,
            int(0),
            int(2),
            &mut crate::polynomial::real::Budget::new(
                crate::polynomial::RootIsolationOptions::default(),
            ),
        )
        .unwrap();
        let g = Arc::new(Gen::new(poly.clone(), roots.last().unwrap().clone()));
        let a = vec![r(1, 3), r(-2, 5), r(7, 11), r(3, 2), r(-1, 7)];
        let b = vec![r(-4, 9), int(2), r(1, 13), zero(), r(5, 3)];
        assert_eq!(g.mul_mod(&a, &b), pmod(&pmul(&a, &b), &poly));
        assert_eq!(
            g.mul_mod(&a[..2], &b[..3]),
            pmod(&pmul(&a[..2], &b[..3]), &poly)
        );
        let (x, y) = (K::norm(&g, a.clone()), K::norm(&g, b.clone()));
        let xy = x.mul(&y);
        assert_eq!(xy, K::norm(&g, pmod(&pmul(&a, &b), &poly)));
        let exact = g
            .narrowed(SIGN_STEPS)
            .sign_polynomial(&IntPolynomial::from_rationals(&xy.poly()));
        assert_eq!(xy.sign(), exact);
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
