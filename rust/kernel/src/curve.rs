//! Certified nonperiodic rational B-spline and Bezier curve evaluation.
//!
//! OCCT references: Geom_BSplineCurve::CheckCurveData, BSplCLib::Eval/Bohm,
//! BSplCLib_D0/D1/D2, and PLib::RationalDerivative. Exact de Boor interpolation
//! with differentiated recurrences replaces floating evaluation. Every output
//! component has its smallest finite binary64 enclosure. See MATHEMATICS.md.
use crate::{exact, interval, math::finite, Error, Point3, Result, ScalarInterval};
use num_bigint::BigInt;
use num_rational::BigRational as R;

pub const MAX_DEGREE: usize = 25;
pub const MAX_POLES: usize = 4096;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivativeOrder {
    Position,
    First,
    Second,
}
impl DerivativeOrder {
    fn count(self) -> usize {
        match self {
            Self::Position => 0,
            Self::First => 1,
            Self::Second => 2,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KnotSide {
    /// At an interior knot, require the requested derivatives to agree exactly
    /// from both sides. At a domain endpoint, use the only interior side.
    Automatic,
    /// Left limit within the curve's domain; invalid at the first parameter.
    Left,
    /// Right limit within the curve's domain; invalid at the last parameter.
    Right,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveEvaluation {
    position: [ScalarInterval; 3],
    derivatives: [Option<[ScalarInterval; 3]>; 2],
}
impl CurveEvaluation {
    pub fn position_bounds(self) -> [ScalarInterval; 3] {
        self.position
    }
    /// Rounded representative; retain the bounds for subsequent decisions.
    pub fn position(self) -> Point3 {
        let p = self.position.map(ScalarInterval::representative);
        Point3::new(p[0], p[1], p[2])
    }
    /// Bounds for derivative 1 or 2, or None if it was not requested.
    /// These are derivatives with respect to the original knot parameter.
    pub fn derivative_bounds(self, order: usize) -> Option<[ScalarInterval; 3]> {
        match order {
            1 | 2 => self.derivatives[order - 1],
            _ => None,
        }
    }
}

/// Immutable nonperiodic curve, with positive weights and distinct ordered
/// knots plus multiplicities. Clamped and unclamped knot vectors are accepted.
#[derive(Debug, Clone)]
pub struct BSplineCurve3 {
    degree: usize,
    poles: Vec<Point3>,
    weights: Vec<f64>,
    knots: Vec<f64>,
    multiplicities: Vec<usize>,
    flat_knots: Vec<f64>,
}
impl BSplineCurve3 {
    /// Degree 1..=25; at most 4096 poles. Interior multiplicities are 1..=degree,
    /// end multiplicities 1..=degree+1, and sum(mults)=poles.len()+degree+1.
    /// All coordinates/knots are finite; every supplied weight is positive and
    /// finite. No knot/weight is snapped, normalized or discarded by tolerance.
    pub fn new(
        degree: usize,
        poles: Vec<Point3>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<usize>,
    ) -> Result<Self> {
        if degree == 0 || degree > MAX_DEGREE {
            return Err(Error::InvalidCurve("degree must be in 1..=25"));
        }
        if poles.len() > MAX_POLES || knots.len() > MAX_POLES + MAX_DEGREE + 1 {
            return Err(Error::LimitExceeded("spline control data"));
        }
        if poles.len() <= degree || knots.len() < 2 || knots.len() != multiplicities.len() {
            return Err(Error::InvalidCurve("pole/knot/multiplicity counts"));
        }
        for p in &poles {
            exact::point(*p)?;
        }
        for &k in &knots {
            finite(k, "spline knot")?;
        }
        if knots.windows(2).any(|k| k[0] >= k[1]) {
            return Err(Error::InvalidCurve("knots must be strictly increasing"));
        }
        let mut total = 0;
        for (i, &m) in multiplicities.iter().enumerate() {
            let end = i == 0 || i + 1 == knots.len();
            if m == 0 || m > degree + usize::from(end) {
                return Err(Error::InvalidCurve("knot multiplicity"));
            }
            total += m; // Bounded counts and validated degree prevent overflow.
        }
        if total != poles.len() + degree + 1 {
            return Err(Error::InvalidCurve("sum of knot multiplicities"));
        }
        let weights = weights.unwrap_or_else(|| vec![1.; poles.len()]);
        if weights.len() != poles.len() {
            return Err(Error::InvalidCurve("weight count"));
        }
        for &w in &weights {
            finite(w, "spline weight")?;
            if w <= 0. {
                return Err(Error::InvalidCurve("weights must be positive"));
            }
        }
        let flat_knots: Vec<_> = knots
            .iter()
            .zip(&multiplicities)
            .flat_map(|(&k, &m)| std::iter::repeat_n(k, m))
            .collect();
        if flat_knots[degree] >= flat_knots[poles.len()] {
            return Err(Error::InvalidCurve("empty parameter domain"));
        }
        Ok(Self {
            degree,
            poles,
            weights,
            knots,
            multiplicities,
            flat_knots,
        })
    }
    pub fn degree(&self) -> usize {
        self.degree
    }
    pub fn poles(&self) -> &[Point3] {
        &self.poles
    }
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }
    pub fn multiplicities(&self) -> &[usize] {
        &self.multiplicities
    }
    pub fn is_rational(&self) -> bool {
        self.weights.windows(2).any(|w| w[0] != w[1])
    }
    pub fn domain(&self) -> (f64, f64) {
        (
            self.flat_knots[self.degree],
            self.flat_knots[self.poles.len()],
        )
    }

    /// Evaluate within the closed domain, without extrapolation. An automatic
    /// derivative at a knot is rejected if the exact one-sided jets disagree.
    /// Unrepresentable requested components fail atomically; requesting only a
    /// position does not attempt potentially overflowing derivatives.
    pub fn evaluate(
        &self,
        u: f64,
        order: DerivativeOrder,
        side: KnotSide,
    ) -> Result<CurveEvaluation> {
        finite(u, "spline parameter")?;
        let (start, end) = self.domain();
        if u < start
            || u > end
            || (u == start && side == KnotSide::Left)
            || (u == end && side == KnotSide::Right)
        {
            return Err(Error::OutOfDomain("spline parameter or requested side"));
        }
        let count = order.count();
        let left = side == KnotSide::Left || u == end;
        let span = self
            .flat_knots
            .partition_point(|&k| if left { k < u } else { k <= u })
            - 1;
        let values = self.exact_jet(u, span, count);
        if side == KnotSide::Automatic && count > 0 && u > start && u < end {
            if let Ok(k) = self.knots.binary_search_by(|k| k.partial_cmp(&u).unwrap()) {
                if count > self.degree - self.multiplicities[k] {
                    let other = self.flat_knots.partition_point(|&k| k < u) - 1;
                    if values != self.exact_jet(u, other, count) {
                        return Err(Error::DiscontinuousDerivative);
                    }
                }
            }
        }
        let bounds = |v: &[R; 3]| -> Result<[ScalarInterval; 3]> {
            Ok([enclose(&v[0])?, enclose(&v[1])?, enclose(&v[2])?])
        };
        Ok(CurveEvaluation {
            position: bounds(&values[0])?,
            derivatives: [
                if count >= 1 {
                    Some(bounds(&values[1])?)
                } else {
                    None
                },
                if count >= 2 {
                    Some(bounds(&values[2])?)
                } else {
                    None
                },
            ],
        })
    }

    fn exact_jet(&self, u: f64, span: usize, order: usize) -> Vec<[R; 3]> {
        let p = self.degree;
        let u = rational(u);
        // Each local homogeneous pole carries value, first and second derivative.
        let mut d: Vec<[[R; 3]; 4]> = (span - p..=span)
            .map(|i| {
                let w = rational(self.weights[i]);
                let xyz = self.poles[i].to_array();
                std::array::from_fn(|c| {
                    [
                        if c == 3 {
                            w.clone()
                        } else {
                            rational(xyz[c]) * &w
                        },
                        zero(),
                        zero(),
                    ]
                })
            })
            .collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = span - p + j;
                let low = rational(self.flat_knots[i]);
                let width = rational(self.flat_knots[i + p - r + 1]) - &low;
                // This interval contains the selected positive-width span.
                debug_assert!(width > zero());
                let alpha = (&u - low) / &width;
                let complement = R::from_integer(BigInt::from(1)) - &alpha;
                d[j] = std::array::from_fn(|c| {
                    std::array::from_fn(|n| {
                        if n > order {
                            return zero();
                        }
                        let mut value = &complement * &d[j - 1][c][n] + &alpha * &d[j][c][n];
                        if n > 0 {
                            value += R::from_integer(BigInt::from(n))
                                * (&d[j][c][n - 1] - &d[j - 1][c][n - 1])
                                / &width;
                        }
                        value
                    })
                });
            }
        }
        let h = &d[p];
        debug_assert!(h[3][0] > zero());
        let mut result: Vec<[R; 3]> = Vec::with_capacity(order + 1);
        for n in 0..=order {
            result.push(std::array::from_fn(|c| {
                let mut numerator = h[c][n].clone();
                for k in 1..=n {
                    let binomial = if n == 2 && k == 1 { 2 } else { 1 };
                    numerator -=
                        R::from_integer(BigInt::from(binomial)) * &h[3][k] * &result[n - k][c];
                }
                numerator / &h[3][0]
            }));
        }
        result
    }
}

/// Rational Bezier curve on [0,1], represented exactly as a clamped B-spline.
#[derive(Debug, Clone)]
pub struct BezierCurve3(BSplineCurve3);
impl BezierCurve3 {
    pub fn new(poles: Vec<Point3>, weights: Option<Vec<f64>>) -> Result<Self> {
        let degree = poles
            .len()
            .checked_sub(1)
            .ok_or(Error::InvalidCurve("Bezier pole count"))?;
        if degree == 0 || degree > MAX_DEGREE {
            return Err(Error::InvalidCurve("Bezier degree must be in 1..=25"));
        }
        Ok(Self(BSplineCurve3::new(
            degree,
            poles,
            weights,
            vec![0., 1.],
            vec![degree + 1; 2],
        )?))
    }
    pub fn as_bspline(&self) -> &BSplineCurve3 {
        &self.0
    }
    pub fn evaluate(&self, u: f64, order: DerivativeOrder) -> Result<CurveEvaluation> {
        self.0.evaluate(u, order, KnotSide::Automatic)
    }
}

fn rational(x: f64) -> R {
    R::from_float(x).expect("validated finite spline value")
}
fn zero() -> R {
    R::from_integer(BigInt::from(0))
}
fn enclose(value: &R) -> Result<ScalarInterval> {
    let numerator = value.numer() << 1074_usize;
    interval::enclose(
        |x| numerator.cmp(&(value.denom() * exact::integer(x))),
        "spline evaluation",
    )
}
