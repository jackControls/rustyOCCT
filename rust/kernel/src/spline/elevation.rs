//! Exact degree elevation, referencing BSplCLib::IncreaseDegree and its
//! Prautzsch rank-curve averaging. Sparse control maps are shared across all
//! transverse fields. An independent Cox coefficient solver checks the maps.
//!
//! Unclamped cropping is derived from the original active knot indices, rather
//! than OCCT's output-copy loop, which can change domains or emit bad counts.
//! Periodic data is extended over three periods and cropped by the complete
//! canonical extended knot sequence, preserving the original cyclic origin.
//!
//! The map's coefficients are exact rationals in lowest terms. They are
//! first computed in machine words (`Word`), with the knots moved by one
//! positive affine map onto integers: every insertion ratio
//! `(u - t_i) / (t_{i+p} - t_i)` is invariant under it, so every coefficient
//! is the same rational. Any overflow abandons the words and the map is
//! computed again with `BigRational` (a fuzz timeout's degree-19 periodic
//! axis, `fuzz/regressions/README.md`: the same map, 4.7 s to 0.36 s).
use super::refinement::Sparse;
use super::{integer, ExactKnotVector, KnotRefinementTransform};
use crate::rational;
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::collections::btree_map::Entry;
use std::collections::BTreeMap;

/// The exact arithmetic of the control map. Every operation gives the
/// reduced rational, or `None` when it does not fit the representation.
trait Field: Clone + Sized {
    type Knot: Clone + Ord;
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    /// `(u - a) / (b - a)`, `a < b`.
    fn ratio(u: &Self::Knot, a: &Self::Knot, b: &Self::Knot) -> Option<Self>;
    fn add(&self, other: &Self) -> Option<Self>;
    fn sub(&self, other: &Self) -> Option<Self>;
    fn mul(&self, other: &Self) -> Option<Self>;
    fn div(&self, n: usize) -> Option<Self>;
    fn in_unit(&self) -> bool;
    fn rational(&self) -> R;
}

impl Field for R {
    type Knot = R;
    fn zero() -> Self {
        integer(0)
    }
    fn one() -> Self {
        integer(1)
    }
    fn is_zero(&self) -> bool {
        self == &integer(0)
    }
    fn is_one(&self) -> bool {
        self == &integer(1)
    }
    fn ratio(u: &R, a: &R, b: &R) -> Option<Self> {
        Some(rational::div(&rational::sub(u, a), &rational::sub(b, a)))
    }
    fn add(&self, other: &Self) -> Option<Self> {
        Some(rational::add(self, other))
    }
    fn sub(&self, other: &Self) -> Option<Self> {
        Some(rational::sub(self, other))
    }
    fn mul(&self, other: &Self) -> Option<Self> {
        Some(rational::mul(self, other))
    }
    fn div(&self, n: usize) -> Option<Self> {
        Some(rational::div(self, &integer(n)))
    }
    fn in_unit(&self) -> bool {
        self >= &integer(0) && self <= &integer(1)
    }
    fn rational(&self) -> R {
        self.clone()
    }
}

/// Euclid's gcd, in 64-bit words (one machine division a step) once both
/// operands fit them.
fn gcd_words(mut a: u128, mut b: u128) -> u128 {
    while a > u128::from(u64::MAX) || b > u128::from(u64::MAX) {
        if b == 0 {
            return a;
        }
        (a, b) = (b, a % b);
    }
    let (mut a, mut b) = (a as u64, b as u64);
    while b != 0 {
        (a, b) = (b, a % b);
    }
    u128::from(a)
}

/// A rational `n / d` in lowest terms with `d > 0` (zero as `0/1`), both
/// within `i128` and neither `i128::MIN`: `BigRational`'s value, in words.
/// Products and sums reduce as in Knuth 4.5.1, as `crate::rational` does.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Word {
    n: i128,
    d: i128,
}

impl Word {
    fn new(n: i128, d: i128) -> Option<Self> {
        if n == i128::MIN || d == i128::MIN || d == 0 {
            return None;
        }
        let (n, d) = if d < 0 { (-n, -d) } else { (n, d) };
        let g = gcd_words(n.unsigned_abs(), d.unsigned_abs()) as i128;
        Some(Self { n: n / g, d: d / g })
    }
    fn gcd(a: i128, b: i128) -> i128 {
        gcd_words(a.unsigned_abs(), b.unsigned_abs()) as i128
    }
    fn checked(n: Option<i128>, d: Option<i128>) -> Option<Self> {
        let (n, d) = (n?, d?);
        (n != i128::MIN && d != i128::MIN).then_some(Self { n, d })
    }
}

impl Field for Word {
    type Knot = i128;
    fn zero() -> Self {
        Self { n: 0, d: 1 }
    }
    fn one() -> Self {
        Self { n: 1, d: 1 }
    }
    fn is_zero(&self) -> bool {
        self.n == 0
    }
    fn is_one(&self) -> bool {
        self.n == 1 && self.d == 1
    }
    fn ratio(u: &i128, a: &i128, b: &i128) -> Option<Self> {
        Self::new(u.checked_sub(*a)?, b.checked_sub(*a)?)
    }
    fn add(&self, other: &Self) -> Option<Self> {
        if self.n == 0 {
            return Some(*other);
        }
        if other.n == 0 {
            return Some(*self);
        }
        let g = Self::gcd(self.d, other.d);
        let (ad, bd) = (self.d / g, other.d / g);
        let t = self
            .n
            .checked_mul(bd)?
            .checked_add(other.n.checked_mul(ad)?)?;
        if t == 0 {
            return Some(Self::zero());
        }
        let h = Self::gcd(t, g);
        Self::checked(Some(t / h), ad.checked_mul(other.d / h))
    }
    fn sub(&self, other: &Self) -> Option<Self> {
        self.add(&Self {
            n: other.n.checked_neg()?,
            d: other.d,
        })
    }
    fn mul(&self, other: &Self) -> Option<Self> {
        if self.n == 0 || other.n == 0 {
            return Some(Self::zero());
        }
        let g1 = Self::gcd(self.n, other.d);
        let g2 = Self::gcd(other.n, self.d);
        Self::checked(
            (self.n / g1).checked_mul(other.n / g2),
            (self.d / g2).checked_mul(other.d / g1),
        )
    }
    fn div(&self, n: usize) -> Option<Self> {
        Self::new(self.n, self.d.checked_mul(i128::try_from(n).ok()?)?)
    }
    fn in_unit(&self) -> bool {
        0 <= self.n && self.n <= self.d
    }
    fn rational(&self) -> R {
        R::new(self.n.into(), self.d.into())
    }
}

/// The knots moved by `k -> (k - k_0) / g` onto coprime integers, `g` the
/// rational gcd of their differences from the first, when they fit words.
fn integer_knots(knots: &[R]) -> Option<Vec<i128>> {
    let differences: Vec<R> = knots.iter().map(|k| rational::sub(k, &knots[0])).collect();
    let mut denominator = BigInt::from(1);
    for x in &differences {
        denominator = &denominator / rational::gcd(&denominator, x.denom()) * x.denom();
    }
    let integers: Vec<BigInt> = differences
        .iter()
        .map(|x| x.numer() * (&denominator / x.denom()))
        .collect();
    let g = integers
        .iter()
        .fold(BigInt::from(0), |g, x| rational::gcd(&g, x));
    if g == BigInt::from(0) {
        return None;
    }
    integers
        .iter()
        .map(|x| i128::try_from(&(x / &g)).ok())
        .collect()
}

type Row<F> = Vec<(usize, F)>;

/// `(1 - alpha) left + alpha right`, zero coefficients dropped; rows are
/// sorted by control index.
fn blend<F: Field>(left: &Row<F>, right: &Row<F>, alpha: &F) -> Option<Row<F>> {
    if alpha.is_zero() {
        return Some(left.clone());
    }
    if alpha.is_one() {
        return Some(right.clone());
    }
    let complement = F::one().sub(alpha)?;
    let mut result = Vec::with_capacity(left.len() + right.len());
    let (mut a, mut b) = (left.iter().peekable(), right.iter().peekable());
    loop {
        let (i, x) = match (a.peek(), b.peek()) {
            (Some((i, x)), Some((j, y))) if i == j => {
                let sum = complement.mul(x)?.add(&alpha.mul(y)?)?;
                a.next();
                b.next();
                (*i, sum)
            }
            (Some((i, x)), Some((j, _))) if i < j => {
                a.next();
                (*i, complement.mul(x)?)
            }
            (Some((i, x)), None) => {
                a.next();
                (*i, complement.mul(x)?)
            }
            (_, Some((j, y))) => {
                b.next();
                (*j, alpha.mul(y)?)
            }
            (None, None) => break,
        };
        if !x.is_zero() {
            result.push((i, x));
        }
    }
    Some(result)
}

fn expand<K: Clone>(knots: &[K], mults: &[usize]) -> Vec<K> {
    knots
        .iter()
        .zip(mults)
        .flat_map(|(k, &m)| std::iter::repeat_n(k.clone(), m))
        .collect()
}

fn insert<F: Field>(
    p: usize,
    flat: &mut Vec<F::Knot>,
    rows: &mut Vec<Row<F>>,
    u: &F::Knot,
) -> Option<()> {
    let k = flat.partition_point(|x| x <= u) - 1;
    let mult = k + 1 - flat.partition_point(|x| x < u);
    let (first, last) = (k - p + 1, k - mult);
    let changed = (first..=last)
        .map(|i| {
            let alpha = F::ratio(u, &flat[i], &flat[i + p])?;
            debug_assert!(alpha.in_unit());
            blend(&rows[i - 1], &rows[i], &alpha)
        })
        .collect::<Option<Vec<_>>>()?;
    rows.insert(last + 1, rows[last].clone());
    rows[first..=last].clone_from_slice(&changed);
    flat.insert(k + 1, u.clone());
    Some(())
}

fn increment<F: Field>(
    p: usize,
    knots: &[F::Knot],
    mults: &mut [usize],
    rows: &[Row<F>],
) -> Option<Vec<Row<F>>> {
    let mut sums = vec![BTreeMap::<usize, F>::new(); rows.len() + knots.len() - 1];
    for rank in 0..=p {
        let mut duplicated: Vec<_> = rows
            .iter()
            .enumerate()
            .flat_map(|(i, r)| std::iter::repeat_n(r.clone(), 1 + usize::from(i % (p + 1) == rank)))
            .collect();
        let mut cumulative = 0;
        let mut next = rank + 1;
        let mut missing = Vec::new();
        let selected: Vec<_> = knots
            .iter()
            .zip(mults.iter())
            .map(|(k, &m)| {
                cumulative += m;
                if cumulative >= next {
                    next += p + 1;
                    m + 1
                } else {
                    missing.push(k);
                    m
                }
            })
            .collect();
        let mut flat = expand(knots, &selected);
        debug_assert_eq!(flat.len(), duplicated.len() + p + 2);
        for k in missing {
            insert(p + 1, &mut flat, &mut duplicated, k)?;
        }
        debug_assert_eq!(duplicated.len(), sums.len());
        for (sum, row) in sums.iter_mut().zip(duplicated) {
            for (i, x) in row {
                match sum.entry(i) {
                    Entry::Occupied(mut total) => {
                        let next = total.get().add(&x)?;
                        total.insert(next);
                    }
                    Entry::Vacant(slot) => {
                        slot.insert(x);
                    }
                }
            }
        }
    }
    for m in mults {
        *m += 1;
    }
    sums.into_iter()
        .map(|sum| {
            sum.into_iter()
                .map(|(i, x)| Some((i, x.div(p + 1)?)))
                .collect()
        })
        .collect()
}

/// Elevate the unit rows `units` (a control index, or a zero control) on the
/// clamped working knots from degree `p` to `q`.
fn elevate<F: Field>(
    p: usize,
    q: usize,
    knots: &[F::Knot],
    mults: &[usize],
    units: &[Option<usize>],
) -> Option<Vec<Sparse>> {
    let mut mults = mults.to_vec();
    let mut rows: Vec<Row<F>> = units
        .iter()
        .map(|unit| unit.iter().map(|&i| (i, F::one())).collect())
        .collect();
    for degree in p..q {
        rows = increment(degree, knots, &mut mults, &rows)?;
    }
    Some(
        rows.iter()
            .map(|row| row.iter().map(|(i, x)| (*i, x.rational())).collect())
            .collect(),
    )
}

pub(crate) struct DegreeElevationTransform(KnotRefinementTransform);

impl DegreeElevationTransform {
    /// Both bases and the complete output grid have already been preflighted.
    pub(crate) fn new(old: &ExactKnotVector, new: &ExactKnotVector) -> Self {
        let p = old.degree();
        let n = old.pole_count();
        let (knots, mut mults, mut units, prefix, suffix) = if old.is_periodic() {
            let period = &old.domain()[1] - &old.domain()[0];
            let e = p + 1 - old.multiplicities()[0];
            let mut knots = Vec::new();
            let mut mults = Vec::new();
            // n > p puts the active interval of these three periods strictly
            // around the entire fundamental period. Zero clamping at the two
            // outer ends therefore cannot affect it. After elevation, n' > q
            // and e = q + 1 - m' is unchanged: the canonical knot extension
            // lies strictly inside the neighboring periods, away from either
            // clamped end. Matching that full window preserves the cyclic
            // origin. Extra outer periods add work but no supporting basis.
            for turn in -1..=1 {
                let shift = R::from_integer(turn.into()) * &period;
                for (k, &m) in old
                    .knots()
                    .iter()
                    .zip(old.multiplicities())
                    .take(old.knots().len() - 1)
                {
                    knots.push(k + &shift);
                    mults.push(m);
                }
            }
            knots.push(&old.knots()[0] + integer(2) * &period);
            mults.push(old.multiplicities()[0]);
            let units: Vec<_> = (0..3 * n - e).map(|i| Some((i + e) % n)).collect();
            (knots, mults, units, e, e)
        } else {
            (
                old.knots().to_vec(),
                old.multiplicities().to_vec(),
                (0..n).map(Some).collect(),
                p + 1 - old.multiplicities()[0],
                p + 1 - old.multiplicities().last().unwrap(),
            )
        };
        // Zero units are temporary zero homogeneous controls, never exposed
        // as positive-weight curve or surface data.
        units.splice(0..0, std::iter::repeat_n(None, prefix));
        units.extend(std::iter::repeat_n(None, suffix));
        mults[0] = p + 1;
        *mults.last_mut().unwrap() = p + 1;
        let q = new.degree();
        let rows = integer_knots(&knots)
            .and_then(|words| elevate::<Word>(p, q, &words, &mults, &units))
            .or_else(|| elevate::<R>(p, q, &knots, &mults, &units))
            .expect("rational arithmetic does not overflow");
        for m in &mut mults {
            *m += q - p;
        }
        let offset = if old.is_periodic() {
            expand(&knots, &mults)
                .windows(new.flat.len())
                .position(|w| w == new.flat)
                .expect("elevated periodic support contains canonical extension")
        } else {
            let first = old.knots().binary_search(&old.domain()[0]).unwrap();
            let last = old.knots().binary_search(&old.domain()[1]).unwrap();
            let delta = q - p;
            let offset = prefix + delta * first;
            let tail = suffix + delta * (old.knots().len() - 1 - last);
            let flat = expand(&knots, &mults);
            debug_assert_eq!(&flat[offset..flat.len() - tail], new.flat);
            offset
        };
        let rows = &rows[offset..offset + new.pole_count()];
        debug_assert!(rows.iter().all(|r| !r.is_empty()));
        Self(KnotRefinementTransform::from_rows(rows))
    }

    pub(crate) fn apply(&self, controls: &[[R; 4]]) -> Vec<[R; 4]> {
        self.0.apply(controls)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words() -> Vec<Word> {
        let big = [1i128 << 62, (1 << 100) + 7, 3i128.pow(70), i128::MAX];
        let mut out = Vec::new();
        for n in (-7..=7).chain(big).chain(big.map(|x| -x)) {
            for d in [
                1,
                2,
                3,
                12,
                257,
                (1i128 << 61) - 1,
                10i128.pow(18),
                i128::MAX,
            ] {
                out.extend(Word::new(n, d));
            }
        }
        out
    }

    fn same(w: Option<Word>, r: R) {
        if let Some(w) = w {
            assert!(&BigInt::from(w.n) == r.numer() && &BigInt::from(w.d) == r.denom());
        }
    }

    /// Word operations give `num_rational`'s reduced values or overflow.
    #[test]
    fn words_agree_with_num_rational() {
        let xs = words();
        for (i, a) in xs.iter().enumerate() {
            let ra = R::new(a.n.into(), a.d.into());
            same(Some(*a), ra.clone());
            for b in xs.iter().skip(i % 11).step_by(11) {
                let rb = R::new(b.n.into(), b.d.into());
                same(a.add(b), &ra + &rb);
                same(a.sub(b), &ra - &rb);
                same(a.mul(b), &ra * &rb);
            }
            for n in [1, 2, 25, 26] {
                same(a.div(n), &ra / integer(n));
            }
        }
        assert_eq!(Word::new(1, 0), None);
        assert_eq!(Word::new(i128::MIN, 1), None);
        assert_eq!(Word { n: i128::MAX, d: 1 }.add(&Word::one()), None);
        assert_eq!(Word { n: i128::MAX, d: 1 }.mul(&Word { n: 2, d: 1 }), None);
    }

    /// The words and `BigRational` give the same map on clamped working
    /// axes, scaled far outside binary64 or by `1/257`; knots too far apart
    /// for words leave them to `BigRational`.
    #[test]
    fn word_maps_are_the_rational_maps() {
        let scales = [
            integer(1),
            R::from_integer(BigInt::from(1) << 4096),
            R::new(1.into(), BigInt::from(1) << 1100),
            R::new(1.into(), 257.into()),
        ];
        let uniform = std::iter::once(6)
            .chain([1; 12])
            .chain([6])
            .collect::<Vec<usize>>();
        for scale in &scales {
            for (p, q, knots, mults) in [
                (2, 5, vec![0, 1, 4], vec![3, 1, 3]),
                (3, 7, vec![-2, -1, 0, 1, 2, 3, 4], vec![4, 3, 1, 3, 1, 3, 4]),
                (5, 9, (0..14).collect(), uniform.clone()),
            ] {
                let knots: Vec<R> = knots
                    .into_iter()
                    .map(|k: i32| R::from_integer(k.into()) * scale)
                    .collect();
                let count = mults.iter().sum::<usize>() - p - 1;
                let units: Vec<_> = (0..count).map(Some).collect();
                let words = integer_knots(&knots).expect("integer knots");
                assert!(words.iter().all(|k| k.unsigned_abs() < 64));
                let fast = elevate::<Word>(p, q, &words, &mults, &units).expect("words");
                let exact = elevate::<R>(p, q, &knots, &mults, &units).unwrap();
                assert_eq!(fast, exact);
            }
        }
        let far = [
            integer(0),
            integer(1),
            R::from_integer(BigInt::from(1) << 200),
        ];
        assert_eq!(integer_knots(&far), None);
    }
}
