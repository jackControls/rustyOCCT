//! Validated knot vectors shared by curves and tensor-product surfaces.
//!
//! OCCT references: BSplCLib::NbPoles, KnotSequence, PoleIndex; curve/surface
//! PrepareEval and PLib/BSplSLib::RationalDerivative. Extended periodic knots
//! and parameter reduction are exact, including outside binary64's range.
use crate::{curve::KnotSide, exact, interval, math::finite, Error, Result, ScalarInterval};
use num_bigint::BigInt;
use num_rational::BigRational as R;

mod refinement;
pub(crate) use refinement::KnotRefinementTransform;
mod rational;
pub(crate) use rational::normalize;
pub use rational::ExactKnotVector;

pub const MAX_DEGREE: usize = 25;
/// Also the maximum total number of control points in a surface.
pub const MAX_POLES: usize = 4096;

/// Immutable knot vector. Pole count follows from the degree and multiplicities.
/// Both periodic and nonperiodic vectors require more poles than their degree.
/// Periodic vectors use OCCT's pole ordering and equal end multiplicities.
#[derive(Debug, Clone)]
pub struct KnotVector {
    degree: usize,
    pole_count: usize,
    periodic: bool,
    knots: Vec<f64>,
    multiplicities: Vec<usize>,
    flat: Vec<R>,
    domain: (f64, f64),
}

#[derive(Debug, Clone)]
pub(crate) struct Parameter {
    pub(crate) value: R,
    pub(crate) span: usize,
}

/// One nonempty piece of an explicit query interval. Shifted periodic knots
/// need not be representable as f64; keep both span and clipping bounds exact.
pub(crate) struct KnotSpan {
    pub index: usize,
    pub start: R,
    pub end: R,
    pub lower: R,
    pub upper: R,
}

impl KnotVector {
    /// Distinct, finite increasing knots; degree 1..=25. End multiplicities
    /// are 1..=degree+1, interior multiplicities 1..=degree. No snapping.
    pub fn new(degree: usize, knots: Vec<f64>, multiplicities: Vec<usize>) -> Result<Self> {
        Self::build(degree, knots, multiplicities, false)
    }

    /// Periodic vectors have equal end multiplicities in 1..=degree. Their
    /// pole count is sum(multiplicities) minus one end multiplicity. Finite
    /// query parameters are reduced exactly modulo the knot period.
    pub fn new_periodic(
        degree: usize,
        knots: Vec<f64>,
        multiplicities: Vec<usize>,
    ) -> Result<Self> {
        Self::build(degree, knots, multiplicities, true)
    }

    fn build(
        degree: usize,
        knots: Vec<f64>,
        multiplicities: Vec<usize>,
        periodic: bool,
    ) -> Result<Self> {
        if degree == 0 || degree > MAX_DEGREE {
            return Err(Error::InvalidSpline("degree must be in 1..=25"));
        }
        if knots.len() > MAX_POLES + MAX_DEGREE + 1 {
            return Err(Error::LimitExceeded("spline knot data"));
        }
        if knots.len() < 2 || knots.len() != multiplicities.len() {
            return Err(Error::InvalidSpline("knot/multiplicity counts"));
        }
        for &k in &knots {
            finite(k, "spline knot")?;
        }
        if knots.windows(2).any(|k| k[0] >= k[1]) {
            return Err(Error::InvalidSpline("knots must be strictly increasing"));
        }
        let mut total = 0;
        for (i, &m) in multiplicities.iter().enumerate() {
            let end = i == 0 || i + 1 == knots.len();
            if m == 0 || m > degree + usize::from(end && !periodic) {
                return Err(Error::InvalidSpline("knot multiplicity"));
            }
            total += m; // Degree and knot-count bounds prevent overflow.
        }
        if periodic && multiplicities[0] != multiplicities[knots.len() - 1] {
            return Err(Error::InvalidSpline(
                "periodic end multiplicities must agree",
            ));
        }
        let pole_count = total
            .checked_sub(if periodic {
                multiplicities[0]
            } else {
                degree + 1
            })
            .ok_or(Error::InvalidSpline("sum of knot multiplicities"))?;
        if pole_count > MAX_POLES {
            return Err(Error::LimitExceeded("spline control data"));
        }
        if pole_count <= degree {
            return Err(Error::InvalidSpline("more poles than degree required"));
        }
        let mut flat: Vec<_> = knots
            .iter()
            .zip(&multiplicities)
            .flat_map(|(&k, &m)| std::iter::repeat_n(rational(k), m))
            .collect();
        let domain = if periodic {
            let domain = (knots[0], knots[knots.len() - 1]);
            let period = rational(domain.1) - rational(domain.0);
            let extension = degree + 1 - multiplicities[0];
            let before: Vec<_> = flat[pole_count - extension..pole_count]
                .iter()
                .map(|k| k - &period)
                .collect();
            let after: Vec<_> = flat[multiplicities[0]..multiplicities[0] + extension]
                .iter()
                .map(|k| k + &period)
                .collect();
            flat = before.into_iter().chain(flat).chain(after).collect();
            domain
        } else {
            // Original f64 atoms, recovered without rounding any arithmetic.
            let expanded: Vec<_> = knots
                .iter()
                .zip(&multiplicities)
                .flat_map(|(&k, &m)| std::iter::repeat_n(k, m))
                .collect();
            (expanded[degree], expanded[pole_count])
        };
        if domain.0 >= domain.1 {
            return Err(Error::InvalidSpline("empty parameter domain"));
        }
        Ok(Self {
            degree,
            pole_count,
            periodic,
            knots,
            multiplicities,
            flat,
            domain,
        })
    }

    pub fn degree(&self) -> usize {
        self.degree
    }
    pub fn pole_count(&self) -> usize {
        self.pole_count
    }
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }
    pub fn multiplicities(&self) -> &[usize] {
        &self.multiplicities
    }
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }
    /// Fundamental domain for periodic vectors; [U[degree], U[pole_count]] otherwise.
    pub fn domain(&self) -> (f64, f64) {
        self.domain
    }

    pub(crate) fn pole_index(&self, index: usize) -> usize {
        if self.periodic {
            index % self.pole_count
        } else {
            index
        }
    }

    pub(crate) fn spans(&self) -> impl Iterator<Item = (f64, f64, usize)> + '_ {
        self.knots.windows(2).filter_map(|k| {
            if k[0] < self.domain.0 || k[1] > self.domain.1 {
                return None;
            }
            let start = rational(k[0]);
            Some((k[0], k[1], self.flat.partition_point(|x| x <= &start) - 1))
        })
    }

    pub(crate) fn spans_in(
        &self,
        first: f64,
        last: f64,
        max_spans: usize,
    ) -> Result<Vec<KnotSpan>> {
        finite(first, "curve interval")?;
        finite(last, "curve interval")?;
        if first >= last {
            return Err(Error::OutOfDomain(
                "curve interval must have positive length",
            ));
        }
        if !self.periodic && (first < self.domain.0 || last > self.domain.1) {
            return Err(Error::OutOfDomain("curve interval"));
        }
        let (lower, upper) = (rational(first), rational(last));
        let base: Vec<_> = self
            .spans()
            .map(|(a, b, i)| (rational(a), rational(b), i))
            .collect();
        let period = rational(self.domain.1) - rational(self.domain.0);
        // Count before enumeration. A tiny period and a huge query must return
        // a resource error without iterating an astronomical number of turns.
        let mut count = BigInt::from(0);
        for (a, b, _) in &base {
            count += if self.periodic {
                (((&upper - a) / &period).ceil() - ((&lower - b) / &period).floor() - integer(1))
                    .to_integer()
            } else {
                BigInt::from(usize::from(a < &upper && b > &lower))
            };
            if count > BigInt::from(max_spans) {
                return Err(Error::ComputationLimit("spline span traversal"));
            }
        }
        let mut offset = if self.periodic {
            ((&lower - rational(self.domain.0)) / &period).floor() * &period
        } else {
            zero()
        };
        let mut result = Vec::new();
        loop {
            for (a, b, index) in &base {
                let (start, end) = (a + &offset, b + &offset);
                if start >= upper {
                    break;
                }
                if end <= lower {
                    continue;
                }
                result.push(KnotSpan {
                    index: *index,
                    lower: start.clone().max(lower.clone()),
                    upper: end.clone().min(upper.clone()),
                    start,
                    end,
                });
            }
            if !self.periodic {
                break;
            }
            offset += &period;
            if rational(self.domain.0) + &offset >= upper {
                break;
            }
        }
        debug_assert_eq!(BigInt::from(result.len()), count);
        Ok(result)
    }

    /// One selected side, and (only if continuity is not guaranteed by knot
    /// multiplicity) the other side needed for an automatic derivative check.
    pub(crate) fn locate(&self, u: f64, side: KnotSide, order: usize) -> Result<Vec<Parameter>> {
        finite(u, "spline parameter")?;
        let (start, end) = self.domain;
        if !self.periodic
            && (u < start
                || u > end
                || (u == start && side == KnotSide::Left)
                || (u == end && side == KnotSide::Right))
        {
            return Err(Error::OutOfDomain("spline parameter or requested side"));
        }
        let start = rational(start);
        let end = rational(end);
        let mut u = rational(u);
        if self.periodic {
            let period = &end - &start;
            let turns = ((&u - &start) / &period).floor();
            u -= turns * period;
            debug_assert!(u >= start && u < end);
        }
        let seam = self.periodic && u == start;
        let parameter = |value: R, left: bool| {
            let span = self
                .flat
                .partition_point(|k| if left { k < &value } else { k <= &value })
                - 1;
            Parameter { value, span }
        };
        let mut sides = vec![if seam && side == KnotSide::Left {
            parameter(end.clone(), true)
        } else {
            parameter(
                u.clone(),
                side == KnotSide::Left || (!self.periodic && u == end),
            )
        }];
        if side == KnotSide::Automatic && order > 0 {
            let multiplicity = if seam {
                self.multiplicities[0]
            } else if u > start && u < end {
                self.flat.partition_point(|k| k <= &u) - self.flat.partition_point(|k| k < &u)
            } else {
                0
            };
            if multiplicity > 0 && order > self.degree - multiplicity {
                sides.push(parameter(if seam { end } else { u }, true));
            }
        }
        Ok(sides)
    }
}

/// Exact homogeneous power coefficients on one span, with local parameter
/// t=(u-U[span])/(U[span+1]-U[span]). Interpolate polynomials through de Boor;
/// do not recover coefficients from rounded point samples.
pub(crate) fn span_polynomial<const N: usize>(
    axis: &KnotVector,
    span: usize,
    poles: Vec<[R; N]>,
) -> [Vec<R>; N] {
    span_polynomial_from_knots(axis.degree, &axis.flat, span, poles)
}

pub(crate) fn span_polynomial_from_knots<const N: usize>(
    p: usize,
    knots: &[R],
    span: usize,
    poles: Vec<[R; N]>,
) -> [Vec<R>; N] {
    // One positive denominator per vector of polynomials avoids reducing a
    // rational number at every coefficient product and sum. Every de Boor
    // stage still performs the same exact linear-polynomial interpolation.
    let denominator = common_denominator(poles.iter().flatten());
    let mut d: Vec<_> = poles
        .into_iter()
        .map(|p| {
            (
                p.map(|c| vec![c.numer() * (&denominator / c.denom())]),
                denominator.clone(),
            )
        })
        .collect();
    let start = &knots[span];
    let length = &knots[span + 1] - start;
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = span - p + j;
            let width = &knots[i + p - r + 1] - &knots[i];
            let a = (start - &knots[i]) / &width;
            let b = &length / width;
            let (left, left_den) = &d[j - 1];
            let (right, right_den) = &d[j];
            let mix_den = integer_lcm(a.denom(), b.denom());
            let a_num = a.numer() * (&mix_den / a.denom());
            let b_num = b.numer() * (&mix_den / b.denom());
            let common = integer_lcm(left_den, right_den);
            let (ls, rs) = (&common / left_den, &common / right_den);
            let (l0, r0) = ((&mix_den - &a_num) * &ls, &a_num * &rs);
            let (l1, r1) = (-&b_num * &ls, &b_num * &rs);
            let mut coefficients: [Vec<BigInt>; N] = std::array::from_fn(|c| {
                let mut value = vec![BigInt::from(0); r + 1];
                for (k, entry) in value.iter_mut().enumerate() {
                    if k < r {
                        *entry += &l0 * &left[c][k] + &r0 * &right[c][k];
                    }
                    if k > 0 {
                        *entry += &l1 * &left[c][k - 1] + &r1 * &right[c][k - 1];
                    }
                }
                value
            });
            let mut denominator = common * mix_den;
            let mut content = denominator.clone();
            for value in coefficients.iter().flatten() {
                if content == BigInt::from(1) {
                    break;
                }
                content = integer_gcd(content, value.clone());
            }
            if content > BigInt::from(1) {
                denominator /= &content;
                for value in coefficients.iter_mut().flatten() {
                    *value /= &content;
                }
            }
            d[j] = (coefficients, denominator);
        }
    }
    let (coefficients, denominator) = d.pop().unwrap();
    coefficients.map(|v| {
        v.into_iter()
            .map(|x| R::new(x, denominator.clone()))
            .collect()
    })
}

/// Evaluate the spline blossom at `p-i` lower and `i` upper arguments to
/// obtain Bernstein controls on a clipped span. The arguments are fractions
/// of the *base* span, so periodic translations do not round the knot data.
pub(crate) fn bezier_controls(
    axis: &KnotVector,
    span: usize,
    mut poles: Vec<[R; 4]>,
    lower: &R,
    upper: &R,
) -> Vec<[R; 4]> {
    let p = axis.degree;
    let start = &axis.flat[span];
    let width = &axis.flat[span + 1] - start;
    let a = start + lower * &width;
    let b = start + upper * &width;
    let stage = |row: &mut [[R; 4]], r: usize, u: &R| {
        for j in (r..=p).rev() {
            let i = span - p + j;
            let width = &axis.flat[i + p - r + 1] - &axis.flat[i];
            debug_assert!(width > zero());
            let alpha = (u - &axis.flat[i]) / width;
            debug_assert!(alpha >= zero() && alpha <= integer(1));
            if alpha == zero() {
                row[j] = row[j - 1].clone();
            } else if alpha != integer(1) {
                let complement = integer(1) - &alpha;
                row[j] =
                    std::array::from_fn(|c| &complement * &row[j - 1][c] + &alpha * &row[j][c]);
            }
        }
    };
    let mut result = Vec::with_capacity(p + 1);
    // Share the prefixes containing lower arguments between blossom queries.
    for lower_count in 0..=p {
        let mut row = poles.clone();
        for r in lower_count + 1..=p {
            stage(&mut row, r, &b);
        }
        result.push(row.pop().unwrap());
        if lower_count < p {
            stage(&mut poles, lower_count + 1, &a);
        }
    }
    result.reverse();
    result
}

/// Reusable extraction map for a tensor axis. Exact knot insertion raises the
/// two clipped boundary multiplicities on unit controls, then the active
/// Bernstein block supplies the linear map for every homogeneous grid row.
/// OCCT reference: BSplCLib::InsertKnots and its local insertion schema.
/// Clearing denominators before each dot product avoids a rational reduction
/// at every multiply/add; no coefficient or output rounding is introduced.
pub(crate) struct BezierSpanTransform {
    rows: Vec<(Vec<BigInt>, BigInt)>,
}
impl BezierSpanTransform {
    pub(crate) fn new(axis: &KnotVector, span: usize, lower: &R, upper: &R) -> Self {
        Self::from_knots(axis.degree, &axis.flat, span, lower, upper)
    }

    pub(crate) fn from_knots(p: usize, flat: &[R], span: usize, lower: &R, upper: &R) -> Self {
        let n = p + 1;
        // Exactly the local p+1 controls and their 2p+2 surrounding knots.
        // At an unclamped domain end, use the last equal knot index, which
        // can exceed the last pole index; clamped FindSpan conventions do not
        // apply to this local representation.
        let mut knots = flat[span - p..=span + p + 1].to_vec();
        let width = &flat[span + 1] - &flat[span];
        let a = &flat[span] + lower * &width;
        let b = &flat[span] + upper * width;
        let mut matrix: Vec<Vec<R>> = (0..n)
            .map(|i| (0..n).map(|j| integer(usize::from(i == j))).collect())
            .collect();
        for u in [&a, &b] {
            let mut mult = knots.iter().filter(|k| *k == u).count();
            while mult < p {
                let k = knots.partition_point(|x| x <= u) - 1;
                matrix = (0..=matrix.len())
                    .map(|i| {
                        if i <= k - p {
                            return matrix[i].clone();
                        }
                        if i > k - mult {
                            return matrix[i - 1].clone();
                        }
                        let alpha = (u - &knots[i]) / (&knots[i + p] - &knots[i]);
                        debug_assert!(alpha >= zero() && alpha <= integer(1));
                        if alpha == zero() {
                            return matrix[i - 1].clone();
                        }
                        if alpha == integer(1) {
                            return matrix[i].clone();
                        }
                        let complement = integer(1) - &alpha;
                        matrix[i - 1]
                            .iter()
                            .zip(&matrix[i])
                            .map(|(a, b)| {
                                if a == b {
                                    a.clone()
                                } else {
                                    &complement * a + &alpha * b
                                }
                            })
                            .collect()
                    })
                    .collect();
                knots.insert(k + 1, u.clone());
                mult += 1;
            }
        }
        let k = knots.partition_point(|x| x <= &a) - 1;
        Self {
            rows: matrix[k - p..=k]
                .iter()
                .map(|row| {
                    let denominator = common_denominator(row.iter());
                    let numerators = row
                        .iter()
                        .map(|x| x.numer() * (&denominator / x.denom()))
                        .collect();
                    (numerators, denominator)
                })
                .collect(),
        }
    }

    pub(crate) fn apply(&self, controls: &[[R; 4]]) -> Vec<[R; 4]> {
        let denominator = common_denominator(controls.iter().flatten());
        let numerators: Vec<[BigInt; 4]> = controls
            .iter()
            .map(|p| std::array::from_fn(|c| p[c].numer() * (&denominator / p[c].denom())))
            .collect();
        self.rows
            .iter()
            .map(|(row, d)| {
                let den = &denominator * d;
                std::array::from_fn(|c| {
                    R::new(
                        row.iter().zip(&numerators).map(|(a, b)| a * &b[c]).sum(),
                        den.clone(),
                    )
                })
            })
            .collect()
    }
}

/// A reusable differentiated de Boor map for one parameter. Propagating the
/// final interpolation weight backwards through the triangular de Boor graph
/// constructs all local pole weights and their first/second derivatives once.
/// Tensor rows then use exact integer dot products with shared denominators.
/// This is algebraically the same interpolation/derivative rule as de_boor_exact;
/// no sampling, finite differencing or floating-point accumulation is used.
pub(crate) struct SplineJetTransform(BezierSpanTransform);
impl SplineJetTransform {
    pub(crate) fn new(p: usize, flat: &[R], at: &Parameter, order: usize) -> Self {
        debug_assert!(order <= 2);
        let mut weights: Vec<[R; 3]> = (0..=p).map(|_| std::array::from_fn(|_| zero())).collect();
        weights[p][0] = integer(1);
        for r in (1..=p).rev() {
            for j in r..=p {
                let i = at.span - p + j;
                let width = &flat[i + p - r + 1] - &flat[i];
                debug_assert!(width > zero());
                let alpha = (&at.value - &flat[i]) / &width;
                let complement = integer(1) - &alpha;
                let incoming = weights[j].clone();
                for n in 0..=order {
                    let correction = if n == 0 {
                        zero()
                    } else {
                        integer(n) * &incoming[n - 1] / &width
                    };
                    weights[j][n] = &alpha * &incoming[n] + &correction;
                    weights[j - 1][n] += &complement * &incoming[n] - correction;
                }
            }
        }
        Self(BezierSpanTransform {
            rows: (0..=order)
                .map(|n| {
                    let denominator = common_denominator(weights.iter().map(|w| &w[n]));
                    let numerators = weights
                        .iter()
                        .map(|w| w[n].numer() * (&denominator / w[n].denom()))
                        .collect();
                    (numerators, denominator)
                })
                .collect(),
        })
    }
    pub(crate) fn apply(&self, controls: &[[R; 4]]) -> Vec<[R; 4]> {
        self.0.apply(controls)
    }
}

/// Exact homogeneous affine blend. Clear the eight control denominators once
/// and reduce only the four completed coordinates, rather than every product.
pub(crate) fn blend_homogeneous(left: &[R; 4], right: &[R; 4], alpha: &R) -> [R; 4] {
    if alpha == &zero() {
        return left.clone();
    }
    if alpha == &integer(1) {
        return right.clone();
    }
    let denominator = common_denominator(left.iter().chain(right));
    let a = alpha.numer();
    let b = alpha.denom() - a;
    let output_den = &denominator * alpha.denom();
    std::array::from_fn(|c| {
        R::new(
            &b * left[c].numer() * (&denominator / left[c].denom())
                + a * right[c].numer() * (&denominator / right[c].denom()),
            output_den.clone(),
        )
    })
}

pub(crate) fn common_denominator<'a>(values: impl Iterator<Item = &'a R>) -> BigInt {
    let mut denominator = BigInt::from(1);
    for x in values {
        denominator = integer_lcm(&denominator, x.denom());
    }
    denominator
}

fn integer_lcm(a: &BigInt, b: &BigInt) -> BigInt {
    a / integer_gcd(a.clone(), b.clone()) * b
}

fn integer_gcd(mut a: BigInt, mut b: BigInt) -> BigInt {
    while b != BigInt::from(0) {
        (a, b) = (b.clone(), a % b);
    }
    if a < BigInt::from(0) {
        -a
    } else {
        a
    }
}

/// Differentiate de Boor's homogeneous pole interpolation, along one axis.
/// `previous[n]` supplies (lower jet index, derivative multiplier) for alpha'.
pub(crate) fn de_boor<const N: usize>(
    axis: &KnotVector,
    at: &Parameter,
    poles: Vec<[[R; N]; 4]>,
    previous: [Option<(usize, usize)>; N],
    count: usize,
) -> [[R; N]; 4] {
    de_boor_exact(axis.degree, &axis.flat, at, poles, previous, count)
}

pub(crate) fn de_boor_exact<const N: usize>(
    p: usize,
    flat: &[R],
    at: &Parameter,
    mut poles: Vec<[[R; N]; 4]>,
    previous: [Option<(usize, usize)>; N],
    count: usize,
) -> [[R; N]; 4] {
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = at.span - p + j;
            let low = &flat[i];
            let width = &flat[i + p - r + 1] - low;
            debug_assert!(width > zero());
            let alpha = (&at.value - low) / &width;
            let complement = integer(1) - &alpha;
            poles[j] = std::array::from_fn(|c| {
                std::array::from_fn(|n| {
                    if n >= count {
                        return zero();
                    }
                    let mut value = &complement * &poles[j - 1][c][n] + &alpha * &poles[j][c][n];
                    if let Some((lower, multiplier)) = previous[n] {
                        value += integer(multiplier)
                            * (&poles[j][c][lower] - &poles[j - 1][c][lower])
                            / &width;
                    }
                    value
                })
            });
        }
    }
    poles.pop().unwrap()
}

/// Recursive multivariate quotient rule from X = W * P. Partial indices are
/// topologically ordered and contain every lower derivative through order two.
pub(crate) fn rationalize<const N: usize>(
    h: [[R; N]; 4],
    partials: [(usize, usize); N],
    count: usize,
) -> Vec<[R; 3]> {
    debug_assert!(h[3][0] > zero());
    let mut result: Vec<[R; 3]> = Vec::with_capacity(count);
    for n in 0..count {
        let (u, v) = partials[n];
        result.push(std::array::from_fn(|c| {
            let mut numerator = h[c][n].clone();
            for (k, &(a, b)) in partials.iter().enumerate().take(n + 1).skip(1) {
                if a <= u && b <= v {
                    let lower = partials.iter().position(|&x| x == (u - a, v - b)).unwrap();
                    let binomial = if (u == 2 && a == 1) || (v == 2 && b == 1) {
                        2
                    } else {
                        1
                    };
                    numerator -= integer(binomial) * &h[3][k] * &result[lower][c];
                }
            }
            numerator / &h[3][0]
        }));
    }
    result
}

pub(crate) fn homogeneous<const N: usize>(point: crate::Point3, weight: f64) -> [[R; N]; 4] {
    let w = rational(weight);
    std::array::from_fn(|c| {
        std::array::from_fn(|n| {
            if n != 0 {
                zero()
            } else if c == 3 {
                w.clone()
            } else {
                rational(point.to_array()[c]) * &w
            }
        })
    })
}
pub(crate) fn rational(x: f64) -> R {
    R::from_float(x).expect("validated finite spline value")
}
fn integer(x: usize) -> R {
    R::from_integer(BigInt::from(x))
}
fn zero() -> R {
    integer(0)
}
pub(crate) fn bounds(v: &[R; 3]) -> Result<[ScalarInterval; 3]> {
    let enclose = |value: &R| {
        let numerator = value.numer() << 1074_usize;
        interval::enclose(
            |x| numerator.cmp(&(value.denom() * exact::integer(x))),
            "spline evaluation",
        )
    };
    Ok([enclose(&v[0])?, enclose(&v[1])?, enclose(&v[2])?])
}
