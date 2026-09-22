//! Validated knot vectors shared by curves and tensor-product surfaces.
//!
//! OCCT references: BSplCLib::NbPoles, KnotSequence, PoleIndex; curve/surface
//! PrepareEval and PLib/BSplSLib::RationalDerivative. Extended periodic knots
//! and parameter reduction are exact, including outside binary64's range.
use crate::{curve::KnotSide, exact, interval, math::finite, Error, Result, ScalarInterval};
use num_bigint::BigInt;
use num_rational::BigRational as R;

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
    let mut d: Vec<[Vec<R>; N]> = poles.into_iter().map(|p| p.map(|c| vec![c])).collect();
    let p = axis.degree;
    let start = &axis.flat[span];
    let length = &axis.flat[span + 1] - start;
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = span - p + j;
            let width = &axis.flat[i + p - r + 1] - &axis.flat[i];
            let a = (start - &axis.flat[i]) / &width;
            let b = &length / width;
            d[j] = std::array::from_fn(|c| {
                let mut value = vec![zero(); r + 1];
                for (k, entry) in value.iter_mut().enumerate() {
                    if k < r {
                        *entry += (integer(1) - &a) * &d[j - 1][c][k] + &a * &d[j][c][k];
                    }
                    if k > 0 {
                        *entry += &b * (&d[j][c][k - 1] - &d[j - 1][c][k - 1]);
                    }
                }
                value
            });
        }
    }
    d.pop().unwrap()
}

/// Differentiate de Boor's homogeneous pole interpolation, along one axis.
/// `previous[n]` supplies (lower jet index, derivative multiplier) for alpha'.
pub(crate) fn de_boor<const N: usize>(
    axis: &KnotVector,
    at: &Parameter,
    mut poles: Vec<[[R; N]; 4]>,
    previous: [Option<(usize, usize)>; N],
    count: usize,
) -> [[R; N]; 4] {
    let p = axis.degree;
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = at.span - p + j;
            let low = &axis.flat[i];
            let width = &axis.flat[i + p - r + 1] - low;
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
