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
    exact: Mutex<std::collections::HashMap<Ip, Ordering>>,
    /// A factor of `poly` known to vanish at the root (the gcd with every
    /// polynomial found zero there): a polynomial it divides is zero at
    /// the root without another Sturm-Tarski count.
    divisor: Mutex<Option<Vec<R>>>,
    /// What `poly` gives every root of it alike, shared by their
    /// generators (`Shared`).
    shared: Arc<Shared>,
    /// A dyadic point `x / 2^bits` within `2^-bits` of the root, certified
    /// by its defining polynomial's signs at `(x -+ 1) / 2^bits` inside the
    /// isolator; the finest found is kept.
    tight: Mutex<Option<(u64, BigInt)>>,
    /// The precisions `tight` failed to reach from no point kept (Newton's
    /// iteration from the isolator's narrowings never certified): the same
    /// failure again, not its iteration.
    failed: Mutex<Vec<u64>>,
    /// A dyadic point for signs alone (`sign_near`) where `tight` finds
    /// none: Newton's iteration from narrower isolators (roots so close
    /// that `tight`'s starts lie outside its convergence), the finest kept.
    sign_point: Mutex<Option<(u64, BigInt)>>,
    /// Inverses found (`K::recip`), by the element's coefficients: a
    /// point's coordinates are divided by the same numbers again and again
    /// (S9e.3b: a meeting's place on its curve at every test of a vertex).
    inverses: Mutex<std::collections::HashMap<Ip, Option<Ip>>>,
    /// The isolator's width's bits above one (`filter_steps`), once.
    above: OnceLock<usize>,
}

/// What a generator's polynomial gives every root of it alike, kept by the
/// polynomial and shared by those roots' generators (S9d.4c: a circle's
/// crossings with a torus are the roots of one resultant, each root its
/// own generator, and each crossing's place divides by one polynomial's
/// value there): `x^j mod poly` for `j` from its degree `n` to `2 n - 2`,
/// numerators over one common denominator (S9d.4b.2b: a product reduced
/// without a rational operation per coefficient), computed once; and the
/// inverses of the elements prime to `poly` (an inverse modulo `poly`, one
/// polynomial whatever the root).
#[derive(Debug, Default)]
struct Shared {
    powers: OnceLock<(Vec<Vec<BigInt>>, BigInt)>,
    inverses: Mutex<std::collections::HashMap<Ip, Ip>>,
}

/// `poly`'s `Shared`, kept for the last few hundred polynomials.
fn shared(poly: &[R]) -> Arc<Shared> {
    static KEPT: OnceLock<Mutex<std::collections::HashMap<Ip, Arc<Shared>>>> = OnceLock::new();
    let key = Ip::from_rats(poly);
    let mut kept = KEPT
        .get_or_init(Default::default)
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    if kept.len() >= 256 && !kept.contains_key(&key) {
        kept.clear();
    }
    kept.entry(key).or_default().clone()
}

/// A polynomial with rational coefficients as integer numerators over one
/// positive denominator, trimmed. A field's products and sums are integer
/// ones, not a gcd per coefficient and operation: a product is reduced to
/// the one form (the denominator sharing no factor with all numerators),
/// a sum is not (the next product reduces it), so `same` compares two by
/// value and `canonical` gives the one form, which keys the kept signs and
/// inverses by their integers (`num_rational`'s hash expands every
/// rational into its continued fraction).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(super) struct Ip {
    /// The numerators, then the denominator (one vector: an element of a
    /// field as small as before beside a rational).
    v: Vec<BigInt>,
}

impl Ip {
    fn raw(mut num: Vec<BigInt>, den: BigInt) -> Self {
        num.push(den);
        Ip { v: num }
    }

    fn num(&self) -> &[BigInt] {
        &self.v[..self.v.len() - 1]
    }

    fn den(&self) -> &BigInt {
        self.v.last().expect("a denominator")
    }

    /// `num / den` (`den` nonzero) in its one form.
    fn new(mut num: Vec<BigInt>, mut den: BigInt) -> Self {
        while num
            .last()
            .is_some_and(|c| c.sign() == num_bigint::Sign::NoSign)
        {
            num.pop();
        }
        if num.is_empty() {
            return Ip::raw(num, BigInt::from(1));
        }
        if den.sign() == num_bigint::Sign::Minus {
            den = -den;
            for c in &mut num {
                *c = -&*c;
            }
        }
        let one = BigInt::from(1);
        let mut g = den.clone();
        for c in &num {
            if g == one {
                break;
            }
            if c.sign() != num_bigint::Sign::NoSign {
                g = crate::rational::gcd(&g, c);
            }
        }
        if g != one {
            for c in &mut num {
                *c /= &g;
            }
            den /= &g;
        }
        Ip::raw(num, den)
    }

    fn from_rats(p: &[R]) -> Self {
        let (num, den) = integral(p);
        Ip::new(num, den)
    }

    /// The rational coefficients.
    fn rats(&self) -> Vec<R> {
        self.num()
            .iter()
            .map(|c| ratio(c.clone(), self.den()))
            .collect()
    }

    fn len(&self) -> usize {
        self.num().len()
    }

    /// `a + s b`, `s` one or minus one.
    fn add(&self, o: &Self, minus: bool) -> Self {
        let n = self.num().len().max(o.num().len());
        let zero = BigInt::from(0);
        if self.den() == o.den() {
            let num = (0..n)
                .map(|i| {
                    let (x, y) = (
                        self.num().get(i).unwrap_or(&zero),
                        o.num().get(i).unwrap_or(&zero),
                    );
                    if minus {
                        x - y
                    } else {
                        x + y
                    }
                })
                .collect();
            return Ip::loose(num, self.den().clone());
        }
        let g = crate::rational::gcd(self.den(), o.den());
        let (ka, kb) = (o.den() / &g, self.den() / &g);
        let num = (0..n)
            .map(|i| {
                let x = self.num().get(i).unwrap_or(&zero) * &ka;
                let y = o.num().get(i).unwrap_or(&zero) * &kb;
                if minus {
                    x - y
                } else {
                    x + y
                }
            })
            .collect();
        Ip::loose(num, self.den() * &ka)
    }

    /// Times a rational `a / b`: in its one form where it was, by gcds of
    /// `a` with the denominator and of `b` with the numerators (each with
    /// a factor as small as the rational's).
    fn scale(&self, k: &R) -> Self {
        let (a, b) = (k.numer(), k.denom());
        if a.sign() == num_bigint::Sign::NoSign || self.num().is_empty() {
            return Ip::zero();
        }
        let one = BigInt::from(1);
        let ga = crate::rational::gcd(a, self.den());
        let mut gb = b.clone();
        for c in self.num() {
            if gb == one {
                break;
            }
            gb = crate::rational::gcd(&gb, c);
        }
        let (a, b) = (a / &ga, b / &gb);
        Ip::raw(
            self.num().iter().map(|c| c / &gb * &a).collect(),
            self.den() / &ga * &b,
        )
    }

    /// A sum's numerators over a denominator, trimmed, the denominator
    /// positive, without the content's gcd (a product reduces it): equal
    /// polynomials may differ in this form, compared by `same` and keyed
    /// by `canonical`.
    fn loose(mut num: Vec<BigInt>, mut den: BigInt) -> Self {
        while num
            .last()
            .is_some_and(|c| c.sign() == num_bigint::Sign::NoSign)
        {
            num.pop();
        }
        if num.is_empty() {
            return Ip::zero();
        }
        if den.sign() == num_bigint::Sign::Minus {
            den = -den;
            for c in &mut num {
                *c = -&*c;
            }
        }
        Ip::raw(num, den)
    }

    /// The one form.
    fn canonical(&self) -> Self {
        Ip::new(self.num().to_vec(), self.den().clone())
    }

    /// Whether the two are one polynomial (numerators cross-multiplied).
    fn same(&self, o: &Self) -> bool {
        self.num().len() == o.num().len()
            && self
                .num()
                .iter()
                .zip(o.num())
                .all(|(x, y)| x * o.den() == y * self.den())
    }

    fn neg(&self) -> Self {
        Ip::raw(self.num().iter().map(|c| -c).collect(), self.den().clone())
    }

    fn zero() -> Self {
        Ip::new(Vec::new(), BigInt::from(1))
    }

    /// The product as polynomials (no reduction by a field's polynomial).
    fn mul(&self, o: &Self) -> Self {
        if self.num().is_empty() || o.num().is_empty() {
            return Ip::zero();
        }
        let mut out = vec![BigInt::from(0); self.num().len() + o.num().len() - 1];
        for (i, x) in self.num().iter().enumerate() {
            if x.sign() == num_bigint::Sign::NoSign {
                continue;
            }
            for (j, y) in o.num().iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        Ip::new(out, self.den() * o.den())
    }

    /// Quotient and remainder by `b` (nonzero), by pseudo-division of the
    /// numerators: `L^e A = Q B + R` for `B`'s leading `L`, so `a = (Q d_b /
    /// (L^e d_a)) b + R / (L^e d_a)`.
    fn divmod(&self, b: &Self) -> (Self, Self) {
        let bn = b.num();
        debug_assert!(!bn.is_empty(), "a nonzero divisor");
        if self.num().len() < bn.len() {
            return (Ip::zero(), self.clone());
        }
        let lead = bn.last().expect("a nonzero divisor");
        let mut r = self.num().to_vec();
        let mut q = vec![BigInt::from(0); r.len() - bn.len() + 1];
        let mut scale = BigInt::from(1);
        while r.len() >= bn.len() {
            let shift = r.len() - bn.len();
            let t = r.last().expect("nonempty").clone();
            for c in r.iter_mut() {
                *c *= lead;
            }
            for c in q.iter_mut() {
                *c *= lead;
            }
            scale *= lead;
            q[shift] += &t;
            for (i, x) in bn.iter().enumerate() {
                r[i + shift] -= &t * x;
            }
            while r
                .last()
                .is_some_and(|c| c.sign() == num_bigint::Sign::NoSign)
            {
                r.pop();
            }
        }
        let den = &scale * self.den();
        (
            Ip::new(q.into_iter().map(|c| c * b.den()).collect(), den.clone()),
            Ip::new(r, den),
        )
    }
}

/// The bisections of the isolator a sign query starts from (those of
/// `Qd::interval`'s enclosure, so the two share it).
pub(super) const SIGN_STEPS: usize = 96;

impl Gen {
    pub(super) fn new(poly: Vec<R>, root: AlgebraicRoot) -> Self {
        Self {
            shared: shared(&poly),
            poly,
            root,
            narrowed: Mutex::new(Vec::new()),
            exact: Mutex::new(std::collections::HashMap::new()),
            divisor: Mutex::new(None),
            tight: Mutex::new(None),
            failed: Mutex::new(Vec::new()),
            sign_point: Mutex::new(None),
            inverses: Mutex::new(std::collections::HashMap::new()),
            above: OnceLock::new(),
        }
    }

    /// `x^j mod poly`, `j` in `n..=2 n - 2`, over a common denominator.
    fn powers(&self) -> &(Vec<Vec<BigInt>>, BigInt) {
        self.shared.powers.get_or_init(|| {
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
    /// operands' denominators, the high terms by `powers`, one content
    /// reduction (the same polynomial as `pmod(pmul(p, q))`).
    fn mul_mod(&self, p: &Ip, q: &Ip) -> Ip {
        if p.num().is_empty() || q.num().is_empty() {
            return Ip::new(Vec::new(), BigInt::from(1));
        }
        let n = self.poly.len() - 1;
        debug_assert!(p.len() <= n && q.len() <= n, "reduced operands");
        let (pi, dp) = (p.num(), p.den());
        let (qi, dq) = (q.num(), q.den());
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
        Ip::new(out, den)
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

    /// `steps` bisections more than the isolator's width has bits above one
    /// (an eliminant's root isolated over its root bound, `1e66` wide, is
    /// that wide still after `SIGN_STEPS`): an isolator narrower than
    /// `2^-steps` for the filters, whose decisions it does not change.
    fn filter_steps(&self, steps: usize) -> usize {
        steps
            + *self.above.get_or_init(|| {
                let (lo, hi) = self.root.isolator();
                let w = hi - lo;
                (w.numer().bits() as i64 - w.denom().bits() as i64 + 1).max(0) as usize
            })
    }

    /// The exact sign of a polynomial at the root (Sturm-Tarski on a
    /// narrowed isolator: the same root, so the same sign), a binary64
    /// enclosure over the isolator first (S9d.4b.2b: certain where it
    /// excludes zero).
    fn sign_of(&self, p: &Ip) -> Ordering {
        // A few ulps' enclosures of the rationals (no exact comparisons).
        let v = self.fast_root(SIGN_STEPS).and_then(|x| {
            p.num()
                .iter()
                .rev()
                .try_fold(Fast::exact_f64(0.0), |acc, c| {
                    Some(acc.mul(&x).add(&Fast::near_parts(c, p.den())?))
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
                    if pmod(&p.rats(), d).is_empty() {
                        return Ordering::Equal;
                    }
                }
                let known = self
                    .exact
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .get(&p.canonical())
                    .copied();
                if let Some(s) = known {
                    return s;
                }
                // A dyadic point near the root and the polynomial's
                // derivative bounded about it: certain where the value there
                // exceeds the bound.
                if let Some(s) = self.sign_near(p) {
                    return s;
                }
                {
                    let root = self.narrowed(self.filter_steps(SIGN_STEPS));
                    let ip = IntPolynomial::new(p.num().to_vec());
                    let s = root.sign_polynomial(&ip);
                    self.exact
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .insert(p.canonical(), s);
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
                }
            }
        }
    }

    /// A dyadic point `x / 2^b` with `|alpha - x / 2^b| <= 2^-b`, `b` at
    /// least `want`: Newton's iteration on the root's defining polynomial
    /// in integers, the precision doubling each step, from the finest point
    /// found before or the middle of the isolator narrowed to `2^-24` (to
    /// `2^-96` where that fails: roots that close), then certified by the
    /// polynomial's opposite nonzero signs at `(x -+ 1) / 2^b`, both inside
    /// the isolator (where it has one root). `None` where that fails too (a
    /// rational root): the exact path decides.
    fn tight(&self, want: u64) -> Option<(BigInt, u64)> {
        let mut kept = self.tight.lock().unwrap_or_else(|e| e.into_inner());
        if let Some((b, x)) = kept.as_ref() {
            if *b >= want {
                return Some((x.clone(), *b));
            }
        }
        if self.root.rational_value().is_some() || self.root.defining().0.len() < 2 {
            return None;
        }
        if kept.is_none()
            && self
                .failed
                .lock()
                .unwrap_or_else(|e| e.into_inner())
                .contains(&want)
        {
            return None;
        }
        let start = |steps: usize| {
            let root = self.narrowed(self.filter_steps(steps));
            let (lo, hi) = root.isolator();
            let w = hi - lo;
            let b = (w.denom().bits() as i64 - w.numer().bits() as i64).clamp(16, 4096) as u64;
            let m = (lo + hi) * R::from_integer(BigInt::from(1) << (b - 1) as usize);
            (m.floor().to_integer(), b)
        };
        let found = match kept.as_ref() {
            Some((b, x)) => self.newton(x.clone(), *b, want + 16),
            None => {
                let (x, b) = start(24);
                let found = self.newton(x, b, want + 16).or_else(|| {
                    let (x, b) = start(SIGN_STEPS);
                    self.newton(x, b, want + 16)
                });
                if found.is_none() {
                    self.failed
                        .lock()
                        .unwrap_or_else(|e| e.into_inner())
                        .push(want);
                }
                found
            }
        }?;
        *kept = Some((found.1, found.0.clone()));
        Some(found)
    }

    /// A point for signs alone where `tight` has none (`sign_point`):
    /// Newton's iteration from the isolator narrowed by 192, 384, 768 and
    /// 1,536 bisections more than its width's bits, each certified as
    /// `tight`'s. Only `sign_near` takes it, and `tight` is not asked again
    /// for it: `tight`'s point (the binary64 views', from the isolator where
    /// it has none) is the one it was.
    fn sign_point(&self, want: u64) -> Option<(BigInt, u64)> {
        if self.root.rational_value().is_some() || self.root.defining().0.len() < 2 {
            return None;
        }
        let mut kept = self.sign_point.lock().unwrap_or_else(|e| e.into_inner());
        if let Some((b, x)) = kept.as_ref() {
            if *b >= want {
                return Some((x.clone(), *b));
            }
        }
        let start = |steps: usize| {
            let root = self.narrowed(self.filter_steps(steps));
            let (lo, hi) = root.isolator();
            let w = hi - lo;
            let b = (w.denom().bits() as i64 - w.numer().bits() as i64).clamp(16, 1 << 20) as u64;
            let m = (lo + hi) * R::from_integer(BigInt::from(1) << (b - 1) as usize);
            (m.floor().to_integer(), b)
        };
        let found = match kept.as_ref() {
            Some((b, x)) => self.newton(x.clone(), *b, want + 16),
            None => [192, 384, 768, 1536].iter().find_map(|&steps| {
                let (x, b) = start(steps);
                self.newton(x, b, (want + 16).max(b))
            }),
        }?;
        *kept = Some((found.1, found.0.clone()));
        Some(found)
    }

    /// `tight`'s iteration from `x / 2^b` to `target` bits and its
    /// certificate.
    fn newton(&self, mut x: BigInt, mut b: u64, target: u64) -> Option<(BigInt, u64)> {
        let f = &self.root.defining().0;
        let d = f.len() - 1;
        let mut steps = 0;
        while b < target {
            steps += 1;
            if steps > 40 {
                return None;
            }
            let next = (2 * b).min(target);
            // `F = f(x) 2^(b d)` and `G = f'(x) 2^(b (d - 1))`.
            let mut fx = f[d].clone();
            for i in (0..d).rev() {
                fx = fx * &x + (&f[i] << (b as usize * (d - i)));
            }
            let mut gx = &f[d] * BigInt::from(d);
            for j in (0..d - 1).rev() {
                gx = gx * &x + ((&f[j + 1] * BigInt::from(j + 1)) << (b as usize * (d - 1 - j)));
            }
            if gx.sign() == num_bigint::Sign::NoSign {
                return None;
            }
            // The step `f / f' = F / (G 2^b)` at the next precision.
            let num = fx << (next - b) as usize;
            x = (x << (next - b) as usize) - num_integer::Integer::div_floor(&num, &gx);
            b = next;
        }
        let den = BigInt::from(1) << b as usize;
        let lo = R::new(&x - 1, den.clone());
        let hi = R::new(&x + 1, den);
        let (ilo, ihi) = self.root.isolator();
        if lo < *ilo || hi > *ihi {
            return None;
        }
        let poly = self.root.defining();
        let (sl, sh) = (poly.sign_at(&lo), poly.sign_at(&hi));
        if sl == Ordering::Equal || sh == Ordering::Equal || sl == sh {
            return None;
        }
        Some((x, b))
    }

    /// A binary64 interval holding the root: about the dyadic point at
    /// `steps` bits, else the isolator narrowed by `steps` bisections more
    /// than its width's bits.
    fn fast_root(&self, steps: usize) -> Option<Fast> {
        if let Some((x, b)) = self.tight(steps as u64) {
            let den = BigInt::from(1) << b as usize;
            return Some(
                Fast::near_parts(&(&x - 1), &den)?.union(&Fast::near_parts(&(&x + 1), &den)?),
            );
        }
        let root = self.narrowed(self.filter_steps(steps));
        let (lo, hi) = root.isolator();
        Some(Fast::near_r(lo)?.union(&Fast::near_r(hi)?))
    }

    /// The sign of a polynomial at the root from a dyadic point `m = x /
    /// 2^b` within `2^-b` of it: `p(m)` exactly (integers over `2^(b k)`),
    /// and `|p(alpha) - p(m)| <= 2^-b max |p'|` over `|t| <= (|x| + 1) /
    /// 2^b`, bounded by the absolute coefficients; certain where `|p(m)|`
    /// exceeds that bound. At three precisions from the coefficients'
    /// size; `None` where none decides (a zero, or a value too small).
    fn sign_near(&self, p: &Ip) -> Option<Ordering> {
        if p.len() < 2 {
            return None;
        }
        // `tight`'s point where it has one; else (roots so close its
        // Newton starts fail) the point for signs alone.
        let tight = self.size_bits(p.num(), p.den());
        let point = |want: u64| {
            if tight.is_some() {
                self.tight(want)
            } else {
                self.sign_point(want)
            }
        };
        let t = match tight {
            Some(t) => t,
            None => self.size_bits_at(p.num(), p.den(), point)?,
        };
        for want in [t + 96, 2 * t + 384, 4 * t + 1536] {
            let (v, bound, _) = self.near_value_at(p.num(), want, point)?;
            if v.magnitude() > bound.magnitude() {
                return Some(if v.sign() == num_bigint::Sign::Minus {
                    Ordering::Less
                } else {
                    Ordering::Greater
                });
            }
        }
        None
    }

    /// The bits above one of a polynomial's largest term near the root
    /// (`num / den` its coefficients): its derivative's size there, which
    /// the dyadic point's precision must pass.
    fn size_bits(&self, num: &[BigInt], den: &BigInt) -> Option<u64> {
        self.size_bits_at(num, den, |want| self.tight(want))
    }

    /// `size_bits` with the dyadic point `point` gives.
    fn size_bits_at(
        &self,
        num: &[BigInt],
        den: &BigInt,
        point: impl Fn(u64) -> Option<(BigInt, u64)>,
    ) -> Option<u64> {
        let top = num.iter().map(|c| c.bits()).max().unwrap_or(0) as i64 - den.bits() as i64;
        let (x, b) = point(128)?;
        let mag = (x.bits() as i64 - b as i64).max(0) * (num.len() as i64 - 1);
        Some((top + mag).max(0) as u64)
    }

    /// The integer polynomial `num` (degree `k` at least one) near the
    /// root: `(v, bound, b)` with `|num(alpha) 2^(b k) - v| <= bound`, from
    /// the dyadic point `m = x / 2^b` within `2^-b` of it (`b >= want`):
    /// `v = num(m) 2^(b k)` exactly, `bound` the absolute coefficients'
    /// derivative over `|t| <= (|x| + 1) / 2^b` times `2^(b (k - 1))`.
    fn near_value(&self, num: &[BigInt], want: u64) -> Option<(BigInt, BigInt, u64)> {
        self.near_value_at(num, want, |want| self.tight(want))
    }

    /// `near_value` at the dyadic point `point` gives.
    fn near_value_at(
        &self,
        num: &[BigInt],
        want: u64,
        point: impl Fn(u64) -> Option<(BigInt, u64)>,
    ) -> Option<(BigInt, BigInt, u64)> {
        let k = num.len() - 1;
        let (x, b) = point(want)?;
        let shift = |c: &BigInt, n: usize| c << (b as usize * n);
        let mut v = num[k].clone();
        for i in (0..k).rev() {
            v = v * &x + shift(&num[i], k - i);
        }
        let y = BigInt::from(x.magnitude() + 1u32);
        let mut bound = BigInt::from(num[k].magnitude().clone()) * BigInt::from(k);
        for j in (0..k - 1).rev() {
            let c = BigInt::from(num[j + 1].magnitude().clone()) * BigInt::from(j + 1);
            bound = bound * &y + shift(&c, k - 1 - j);
        }
        Some((v, bound, b))
    }
}

/// A dyadic enclosure `[lo, hi] / 2^e` of a number of the fields.
#[derive(Debug, Clone)]
struct Dy {
    lo: BigInt,
    hi: BigInt,
    e: usize,
}

/// `floor(x / d)` and `ceil(x / d)`, `d > 0`.
fn floor_div(x: &BigInt, d: &BigInt) -> BigInt {
    num_integer::Integer::div_floor(x, d)
}

fn ceil_div(x: &BigInt, d: &BigInt) -> BigInt {
    -num_integer::Integer::div_floor(&-x, d)
}

impl Dy {
    /// `n / d` (`d > 0`) at `e` bits.
    fn ratio(n: &BigInt, d: &BigInt, e: usize) -> Self {
        let s = n << e;
        Dy {
            lo: floor_div(&s, d),
            hi: ceil_div(&s, d),
            e,
        }
    }

    fn add(&self, o: &Self) -> Self {
        debug_assert_eq!(self.e, o.e);
        Dy {
            lo: &self.lo + &o.lo,
            hi: &self.hi + &o.hi,
            e: self.e,
        }
    }

    fn mul(&self, o: &Self) -> Self {
        debug_assert_eq!(self.e, o.e);
        let p = [
            &self.lo * &o.lo,
            &self.lo * &o.hi,
            &self.hi * &o.lo,
            &self.hi * &o.hi,
        ];
        let lo = p.iter().min().expect("four products");
        let hi = p.iter().max().expect("four products");
        Dy {
            lo: lo >> self.e,
            hi: -((-hi) >> self.e),
            e: self.e,
        }
    }

    /// `sqrt(n / d)`, `n >= 0`, `d > 0`, at `e` bits.
    fn sqrt(n: &BigInt, d: &BigInt, e: usize) -> Self {
        let s = n << (2 * e);
        let (lo, hi) = (floor_div(&s, d), ceil_div(&s, d));
        let r = hi.sqrt();
        Dy {
            lo: lo.sqrt(),
            hi: if &r * &r == hi { r } else { r + 1 },
            e,
        }
    }

    /// The sign, where the enclosure excludes zero.
    fn sign(&self) -> Option<Ordering> {
        if self.lo.sign() == num_bigint::Sign::Plus {
            Some(Ordering::Greater)
        } else if self.hi.sign() == num_bigint::Sign::Minus {
            Some(Ordering::Less)
        } else {
            None
        }
    }

    /// The binary64 number both ends round to (to nearest), which the
    /// enclosed number rounds to as well; `None` where they differ or lie
    /// outside binary64's normal range.
    fn rounded(&self) -> Option<f64> {
        let lo = rounded_ratio(&self.lo, self.e)?;
        (rounded_ratio(&self.hi, self.e)? == lo).then_some(lo)
    }
}

/// `n / 2^e` rounded to the nearest binary64 number (rounded to odd at 64
/// bits, then to nearest by the conversion of a `u64`), or `None` outside
/// `1e-290 ..= 1e290` in magnitude (zero is zero).
fn rounded_ratio(n: &BigInt, e: usize) -> Option<f64> {
    if n.sign() == num_bigint::Sign::NoSign {
        return Some(0.0);
    }
    let m = n.magnitude();
    // The top 64 bits, the last set where lower bits were dropped.
    let drop = (m.bits() as i64 - 64).max(0) as usize;
    let mut q = (m >> drop).to_u64_digits().first().copied().unwrap_or(0);
    if drop > 0 && m.trailing_zeros().unwrap_or(0) < drop as u64 {
        q |= 1;
    }
    let f = q as f64;
    let exp = drop as i64 - e as i64;
    if !(-1000..=1000).contains(&exp) {
        return None;
    }
    let v = f * f64::powi(2.0, exp as i32);
    if !v.is_finite() || !(1e-290..=1e290).contains(&v.abs()) {
        return None;
    }
    Some(if n.sign() == num_bigint::Sign::Minus {
        -v
    } else {
        v
    })
}

/// An element of `Q` or of one `Q(alpha)`: a polynomial in `alpha`
/// (reduced by its generator's polynomial; a constant is `Rat`).
#[derive(Debug, Clone)]
pub(super) enum K {
    Rat(R),
    Alg(Arc<Gen>, Ip),
}

impl PartialEq for K {
    fn eq(&self, o: &Self) -> bool {
        match (self, o) {
            (K::Rat(a), K::Rat(b)) => a == b,
            (K::Alg(g, p), K::Alg(h, q)) => Arc::ptr_eq(g, h) && p.same(q),
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

#[cfg(test)]
fn pmul(a: &[R], b: &[R]) -> Vec<R> {
    product(a, b)
}

/// The product of two polynomials with rational coefficients: integer
/// products over the operands' common denominators, each coefficient
/// reduced once (the same rationals as a product of rationals term by
/// term, without a gcd per term).
pub(super) fn product(a: &[R], b: &[R]) -> Vec<R> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    if a.len() == 1 || b.len() == 1 {
        let (k, p) = if a.len() == 1 { (&a[0], b) } else { (&b[0], a) };
        return ptrim(p.iter().map(|c| rmul(c, k)).collect());
    }
    let (ai, da) = integral(a);
    let (bi, db) = integral(b);
    // Denominators far apart make the common one far larger than each:
    // then term by term.
    let widest = |p: &[R]| p.iter().map(|c| c.denom().bits()).max().unwrap_or(0);
    if da.bits() > widest(a) + 64 || db.bits() > widest(b) + 64 {
        let mut out = vec![zero(); a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                out[i + j] = radd(&out[i + j], &rmul(x, y));
            }
        }
        return ptrim(out);
    }
    let mut out = vec![BigInt::from(0); ai.len() + bi.len() - 1];
    for (i, x) in ai.iter().enumerate() {
        if x.sign() == num_bigint::Sign::NoSign {
            continue;
        }
        for (j, y) in bi.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    let den = da * db;
    ptrim(out.into_iter().map(|c| ratio(c, &den)).collect())
}

#[cfg(test)]
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
            _ => K::Alg(g.clone(), Ip::from_rats(&p)),
        }
    }

    /// A reduced polynomial's element (a constant is `Rat`).
    fn of_ip(g: &Arc<Gen>, p: Ip) -> Self {
        match p.len() {
            0 => K::Rat(zero()),
            1 => K::Rat(ratio(p.num()[0].clone(), p.den())),
            _ => K::Alg(g.clone(), p),
        }
    }

    /// The element's polynomial in the integer form (a rational's
    /// constant).
    fn ip(&self) -> Ip {
        match self {
            K::Rat(a) => Ip::new(vec![a.numer().clone()], a.denom().clone()),
            K::Alg(_, p) => p.clone(),
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

    /// The common generator of two elements (`None` inside for `Q`), or
    /// `None` for two different fields.
    fn common<'a>(&'a self, o: &'a Self) -> Option<Option<&'a Arc<Gen>>> {
        match (self.gen(), o.gen()) {
            (None, None) => Some(None),
            (Some(g), None) | (None, Some(g)) => Some(Some(g)),
            (Some(g), Some(h)) => Arc::ptr_eq(g, h).then_some(Some(g)),
        }
    }

    /// `a + b` or `a - b`.
    fn sum(&self, o: &Self, minus: bool) -> Self {
        if let (K::Rat(a), K::Rat(b)) = (self, o) {
            return K::Rat(if minus { rsub(a, b) } else { radd(a, b) });
        }
        match self.common(o) {
            None => two_fields(),
            Some(None) => unreachable!("rationals"),
            Some(Some(g)) => {
                let g = g.clone();
                K::of_ip(&g, self.ip().add(&o.ip(), minus))
            }
        }
    }

    pub(super) fn add(&self, o: &Self) -> Self {
        self.sum(o, false)
    }

    pub(super) fn sub(&self, o: &Self) -> Self {
        self.sum(o, true)
    }

    pub(super) fn mul(&self, o: &Self) -> Self {
        match (self, o) {
            (K::Rat(a), K::Rat(b)) => K::Rat(rmul(a, b)),
            (K::Alg(g, p), K::Alg(h, q)) => {
                if !Arc::ptr_eq(g, h) {
                    two_fields()
                }
                K::of_ip(g, g.mul_mod(p, q))
            }
            (K::Alg(_, _), K::Rat(k)) => self.scale(k),
            (K::Rat(k), K::Alg(_, _)) => o.scale(k),
        }
    }

    pub(super) fn scale(&self, k: &R) -> Self {
        match self {
            K::Rat(a) => K::Rat(rmul(a, k)),
            K::Alg(g, p) => K::of_ip(g, p.scale(k)),
        }
    }

    pub(super) fn neg(&self) -> Self {
        match self {
            K::Rat(a) => K::Rat(-a),
            K::Alg(g, p) => K::Alg(g.clone(), p.neg()),
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
                let key = p.canonical();
                if let Some(known) = g
                    .inverses
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .get(&key)
                {
                    return known.as_ref().map(|q| K::Alg(g.clone(), q.clone()));
                }
                // Prime to the polynomial: another root's inverse is this
                // root's too.
                if let Some(q) = g
                    .shared
                    .inverses
                    .lock()
                    .unwrap_or_else(|e| e.into_inner())
                    .get(&key)
                {
                    return Some(K::of_ip(g, q.clone()));
                }
                let inv = Self::inverse(g, p);
                let mut kept = g.inverses.lock().unwrap_or_else(|e| e.into_inner());
                if kept.len() >= 4096 {
                    kept.clear();
                }
                match &inv {
                    None => {
                        kept.insert(key, None);
                    }
                    Some(K::Alg(_, q)) => {
                        kept.insert(key, Some(q.clone()));
                    }
                    Some(K::Rat(_)) => {}
                }
                inv
            }
        }
    }

    /// `recip`'s extended Euclid in `Q[x]` modulo the generator's
    /// polynomial (or a factor of it vanishing at the root), in the
    /// integer form: pseudo-divisions, one content reduction per quotient,
    /// remainder and cofactor (the inverse modulo `m` is one polynomial of
    /// lower degree, whatever the algorithm).
    fn inverse(g: &Arc<Gen>, p: &Ip) -> Option<Self> {
        let mut m = Ip::from_rats(&g.poly);
        let mut whole = true;
        loop {
            // Extended Euclid: s p = gcd (mod m).
            let (mut r0, mut r1) = (m.clone(), p.divmod(&m).1);
            let (mut s0, mut s1) = (Ip::zero(), Ip::new(vec![BigInt::from(1)], BigInt::from(1)));
            while r1.len() > 0 {
                let (quo, rem) = r0.divmod(&r1);
                let s2 = s0.add(&quo.mul(&s1), true);
                (r0, r1) = (r1, rem);
                (s0, s1) = (s1, s2);
            }
            if r0.len() == 1 {
                // `s0 / r0`.
                let k = R::new(r0.den().clone(), r0.num()[0].clone());
                let q = s0.scale(&k).divmod(&Ip::from_rats(&g.poly)).1;
                // Prime to the whole polynomial (no factor divided out):
                // every root's.
                if whole {
                    let mut kept = g.shared.inverses.lock().unwrap_or_else(|e| e.into_inner());
                    if kept.len() >= 4096 {
                        kept.clear();
                    }
                    kept.insert(p.canonical(), q.clone());
                }
                return Some(K::of_ip(g, q));
            }
            // A common factor: zero at alpha, or alpha is a root of the
            // cofactor.
            if g.sign_of(&r0) == Ordering::Equal {
                return None;
            }
            m = m.divmod(&r0).0;
            whole = false;
        }
    }

    /// A dyadic enclosure at `e` bits after the point, from the
    /// generator's dyadic point (`None` where that fails).
    fn dy(&self, e: usize) -> Option<Dy> {
        match self {
            K::Rat(a) => Some(Dy::ratio(a.numer(), a.denom(), e)),
            K::Alg(g, p) => {
                let (num, den) = (p.num(), p.den());
                let t = g.size_bits(num, den)?;
                let (v, bound, b) = g.near_value(num, e as u64 + t + 16)?;
                let s = den << (b as usize * (num.len() - 1));
                Some(Dy {
                    lo: floor_div(&((&v - &bound) << e), &s),
                    hi: ceil_div(&((v + bound) << e), &s),
                    e,
                })
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
                let (n, den) = (p.num(), p.den().clone());
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
                let x = g.fast_root(steps)?;
                p.num()
                    .iter()
                    .rev()
                    .try_fold(Fast::exact_f64(0.0), |acc, c| {
                        Some(acc.mul(&x).add(&Fast::near_parts(c, p.den())?))
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
                p.rats().iter().rev().fold(I::exact(zero()), |acc, c| {
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
            // Then a dyadic enclosure, before the product `a^2 - b^2 d`.
            if let Some(s) = self.dy(128).and_then(|v| v.sign()) {
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

    /// A dyadic enclosure at `e` bits after the point (`K::dy`; the surd's
    /// root by integer square roots).
    fn dy(&self, e: usize) -> Option<Dy> {
        let a = self.a.dy(e)?;
        if self.b.is_zero() {
            return Some(a);
        }
        let b = self.b.dy(e)?;
        Some(a.add(&b.mul(&Dy::sqrt(self.d.numer(), self.d.denom(), e))))
    }

    /// `enclose` in binary64 (`K::enclose_fast`).
    pub(super) fn enclose_fast(&self, n: usize) -> Option<Fast> {
        let a = self.a.enclose_fast(n)?;
        if self.b.is_zero() {
            return Some(a);
        }
        Some(a.add(&self.b.enclose_fast(n)?.mul(&Fast::near_r(&self.d)?.sqrt())))
    }

    /// An enclosure, its generator's isolator narrowed by `n` bisections.
    pub(super) fn enclose_at(&self, n: usize) -> I {
        self.enclose(n)
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
        // Dyadic enclosures from the generator's dyadic point: the number
        // rounded to nearest where both ends round alike.
        for e in [96, 384, 1536] {
            if let Some(v) = self.dy(e).and_then(|d| d.rounded()) {
                return v;
            }
        }
        self.to_f64_isolator()
    }

    /// `to_f64` from the isolator where the dyadic point is not found (a
    /// rational root) or the value is zero or tiny.
    fn to_f64_isolator(&self) -> f64 {
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
        // Then dyadic enclosures, before the products.
        if x.gen().is_some() || y.gen().is_some() {
            let v = x
                .dy(128)
                .zip(y.dy(128))
                .map(|(a, b)| a.add(&b.mul(&Dy::sqrt(e.numer(), e.denom(), 128))));
            if let Some(s) = v.and_then(|v| v.sign()) {
                return s;
            }
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
            return approx_dot_sign(a, b);
        }
    }
    // Dyadic enclosures from the generator's dyadic point next: certain
    // where they exclude zero, before the exact products of a field of
    // degree eight.
    if ga.is_some() || gb.is_some() {
        if let Some(s) = dy_dot_sign(a, b) {
            return s;
        }
    }
    mixed_dot_exact(a, b)
}

/// A part of a number in a key of the kept approximate signs: a rational
/// by its integers, an element of a field by its generator's identity (the
/// memo keeps the generator alive while the key is kept) and its integers.
#[derive(Debug, PartialEq, Eq, Hash)]
enum Part {
    Rat(BigInt, BigInt),
    Alg(usize, Ip),
}

/// Approximate signs kept before the memo starts again.
const APPROX_LIMIT: usize = 1024;

/// The kept approximate signs, with the generators their keys name.
type Approx = std::collections::HashMap<Vec<Part>, (Ordering, Vec<Arc<Gen>>)>;

thread_local! {
    /// The signs of dot products of numbers of two different fields by
    /// enclosures (`approx_dot_sign`), by their numbers: the arrangement
    /// asks one direction's order against another's again and again (a
    /// vertex's every edge, every face). A hit is the value the enclosures
    /// give (the same generators narrowed the same way).
    static APPROX: std::cell::RefCell<Approx> = std::cell::RefCell::new(Approx::new());
}

/// The sign of `sum a_i b_i` over two different fields, by enclosures
/// (`approx_sign`), kept (`APPROX`).
fn approx_dot_sign(a: &[Qd], b: &[Qd]) -> Ordering {
    let mut key = vec![Part::Rat(BigInt::from(a.len()), BigInt::from(b.len()))];
    let mut gens = Vec::new();
    for x in a.iter().chain(b) {
        for k in [&x.a, &x.b] {
            key.push(match k {
                K::Rat(r) => Part::Rat(r.numer().clone(), r.denom().clone()),
                K::Alg(g, p) => {
                    gens.push(g.clone());
                    Part::Alg(Arc::as_ptr(g) as usize, p.clone())
                }
            });
        }
        key.push(Part::Rat(x.d.numer().clone(), x.d.denom().clone()));
    }
    if let Some(s) = APPROX.with(|m| m.borrow().get(&key).map(|(s, _)| *s)) {
        return s;
    }
    let s = approx_sign(|n| {
        a.iter().zip(b).fold(I::exact(zero()), |acc, (x, y)| {
            acc.add(&x.enclose(n).mul(&y.enclose(n)))
        })
    });
    APPROX.with(|m| {
        let mut m = m.borrow_mut();
        if m.len() >= APPROX_LIMIT {
            m.clear();
        }
        m.insert(key, (s, gens));
    });
    s
}

/// The sign of `sum a_i b_i` from dyadic enclosures 128 bits after the
/// point (`None`: undecided, or no dyadic point).
fn dy_dot_sign(a: &[Qd], b: &[Qd]) -> Option<Ordering> {
    let e = 128;
    let mut acc: Option<Dy> = None;
    for (x, y) in a.iter().zip(b) {
        let p = x.dy(e)?.mul(&y.dy(e)?);
        acc = Some(match acc {
            None => p,
            Some(s) => s.add(&p),
        });
    }
    acc?.sign()
}

/// `mixed_dot_sign`'s exact products in one field.
fn mixed_dot_exact(a: &[Qd], b: &[Qd]) -> Ordering {
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
        let ip = |p: &[R]| Ip::from_rats(p);
        assert_eq!(
            g.mul_mod(&ip(&a), &ip(&b)).rats(),
            pmod(&pmul(&a, &b), &poly)
        );
        assert_eq!(
            g.mul_mod(&ip(&a[..2]), &ip(&b[..3])).rats(),
            pmod(&pmul(&a[..2], &b[..3]), &poly)
        );
        let (x, y) = (K::norm(&g, a.clone()), K::norm(&g, b.clone()));
        let xy = x.mul(&y);
        assert_eq!(xy, K::norm(&g, pmod(&pmul(&a, &b), &poly)));
        // Sums, differences and multiples in the integer form are the
        // rational ones, and one polynomial has one form.
        assert_eq!(
            x.add(&y),
            K::norm(&g, psub(&a, &b.iter().map(|c| -c).collect::<Vec<_>>()))
        );
        assert_eq!(x.sub(&y), K::norm(&g, psub(&a, &b)));
        assert_eq!(
            x.scale(&r(-6, 35)),
            K::norm(&g, a.iter().map(|c| c * r(-6, 35)).collect())
        );
        assert_eq!(
            Ip::new(
                vec![int(6).to_integer(), int(-4).to_integer()],
                int(-10).to_integer()
            ),
            Ip::from_rats(&[r(-3, 5), r(2, 5)])
        );
        let exact = g
            .narrowed(SIGN_STEPS)
            .sign_polynomial(&IntPolynomial::new(xy.ip().num().to_vec()));
        assert_eq!(xy.sign(), exact);
        // The filters agree with Sturm-Tarski on the elements' signs, and
        // the dyadic enclosure holds the value the isolator does.
        for k in [&x, &y, &xy, &x.sub(&y), &xy.add(&K::Rat(r(-1, 3)))] {
            let p = k.ip();
            let exact = g
                .narrowed(SIGN_STEPS)
                .sign_polynomial(&IntPolynomial::new(p.num().to_vec()));
            assert_eq!(g.sign_near(&p), Some(exact));
            assert_eq!(k.sign(), exact);
            let d = k.dy(80).unwrap();
            let i = k.enclose(200);
            let scale = R::from_integer(BigInt::from(1) << 80);
            assert!(R::from_integer(d.lo.clone()) <= i.lo() * &scale);
            assert!(R::from_integer(d.hi.clone()) >= i.hi() * &scale);
            let inv = k.recip().unwrap();
            assert_eq!(k.mul(&inv).sub(&K::Rat(int(1))).sign(), Ordering::Equal);
            assert_eq!(inv, k.recip().unwrap());
        }
    }

    /// `sqrt 2` beside a complex pair `c +- 2^-66 i`, `c` within `1e-19`
    /// of it: the real root's isolator need not be narrow, `tight`'s Newton
    /// starts lie outside its convergence and it finds no point (nor is it
    /// asked again), while the point for signs alone, from narrower
    /// isolators, decides the signs Sturm-Tarski gives; `tight` keeps none.
    #[test]
    fn close_roots_decide_signs_without_tight_points() {
        let c = R::new(
            BigInt::from(14_142_135_623_730_950_488u128),
            BigInt::from(10u64).pow(19),
        );
        let eps = R::new(BigInt::from(1), BigInt::from(1) << 132usize);
        let two = int(2);
        // (x^2 - 2)(x^2 - 2 c x + c^2 + 2^-132).
        let quad = [&c * &c + &eps, -(int(2) * &c), int(1)];
        let poly: Vec<R> = (0..5)
            .map(|k| {
                let mut sum = zero();
                for (i, q) in quad.iter().enumerate() {
                    if k >= i && k - i <= 2 {
                        let f = [-two.clone(), zero(), int(1)][k - i].clone();
                        sum += q * f;
                    }
                }
                sum
            })
            .collect();
        let ip = IntPolynomial::from_rationals(&poly);
        let roots = crate::polynomial::real::isolate(
            &ip,
            int(1),
            int(2),
            &mut crate::polynomial::real::Budget::new(
                crate::polynomial::RootIsolationOptions::default(),
            ),
        )
        .unwrap();
        assert_eq!(roots.len(), 1);
        for root in roots {
            let g = Gen::new(poly.clone(), root);
            assert!(g.tight(128).is_none());
            // Values near the roots and their gap's size, and products
            // with coefficients of a few hundred bits.
            let near = R::new(
                BigInt::from(1_414_213_562_373_095_049u64),
                BigInt::from(10u64).pow(18),
            );
            let cases: Vec<Vec<R>> = vec![
                vec![-near.clone(), int(1)],
                vec![-(&near * &near), zero(), int(1)],
                vec![-(&c * &c), zero(), int(1)],
                vec![r(-3, 7), int(1), r(5, 11)],
            ];
            for c in cases {
                let p = Ip::from_rats(&c);
                let big = g.mul_mod(&g.mul_mod(&p, &p), &g.mul_mod(&p, &p));
                for p in [p, big] {
                    let exact = g
                        .narrowed(SIGN_STEPS)
                        .sign_polynomial(&IntPolynomial::new(p.num().to_vec()));
                    assert_ne!(exact, Ordering::Equal);
                    assert_eq!(g.sign_near(&p), Some(exact));
                }
            }
            assert!(g.tight.lock().unwrap().is_none());
            assert!(g.tight(128).is_none());
        }
    }

    /// A ratio of integers over a power of two rounds to nearest: both
    /// ends of a narrow enclosure agree only where the value rounds alike.
    #[test]
    fn dyadic_views_round_to_nearest() {
        // 1/3 at 70 bits after the point, its ends one unit apart.
        let e: usize = 70;
        let lo: BigInt = (BigInt::from(1) << e) / BigInt::from(3);
        let d = Dy {
            lo: lo.clone(),
            hi: &lo + 1,
            e,
        };
        assert_eq!(d.rounded(), Some(1.0 / 3.0));
        assert_eq!(rounded_ratio(&-lo, e), Some(-(1.0 / 3.0)));
        // An end on each side of a rounding boundary: none.
        let one = BigInt::from(1) << 60;
        let half_ulp = BigInt::from(1) << 7;
        let d = Dy {
            lo: &one + &half_ulp - 1,
            hi: &one + &half_ulp + 1,
            e: 60,
        };
        assert_eq!(d.rounded(), None);
        // Square roots enclose: sqrt(2) at 64 bits.
        let s = Dy::sqrt(&BigInt::from(2), &BigInt::from(1), 64);
        assert_eq!(s.rounded(), Some(std::f64::consts::SQRT_2));
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
