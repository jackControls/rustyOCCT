//! Exact homogeneous spline extraction and immutable Bernstein editing.
use super::{BSplineCurve3, BezierCurve3, CurveEvaluation, DerivativeOrder, MAX_DEGREE, MAX_POLES};
use crate::{intersection::ExactPoint3, spline, Error, Result};
use num_bigint::BigInt;
use num_rational::BigRational as R;

fn integer(n: usize) -> R {
    R::from_integer(n.into())
}

fn parameter(value: &R) -> Result<R> {
    if value.denom() == &0.into() {
        return Err(Error::InvalidCurve(
            "rational parameter has zero denominator",
        ));
    }
    Ok(R::new(value.numer().clone(), value.denom().clone()))
}

/// Preflight limit on the number of extracted arcs, including periodic turns.
/// Raising this limit raises the permitted allocation/work; it is not a time
/// limit on exact rational arithmetic. Zero permits no arcs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BezierExtractionOptions {
    pub max_arcs: usize,
}
impl Default for BezierExtractionOptions {
    fn default() -> Self {
        Self {
            max_arcs: MAX_POLES,
        }
    }
}

impl BSplineCurve3 {
    /// Extract the fundamental domain without rounding any control data.
    pub fn bezier_arcs(&self) -> Result<Vec<ExactBezierCurve3>> {
        let (first, last) = self.domain();
        self.bezier_arcs_in(first, last)
    }

    /// Extract a positive-length closed interval. Periodic intervals may cross
    /// seams or cover multiple turns. The default limit is 4096 arcs.
    pub fn bezier_arcs_in(&self, first: f64, last: f64) -> Result<Vec<ExactBezierCurve3>> {
        self.bezier_arcs_with_options(first, last, BezierExtractionOptions::default())
    }

    pub fn bezier_arcs_with_options(
        &self,
        first: f64,
        last: f64,
        options: BezierExtractionOptions,
    ) -> Result<Vec<ExactBezierCurve3>> {
        let spans = self.basis.spans_in(first, last, options.max_arcs)?;
        Ok(spans
            .into_iter()
            .map(|span| {
                let controls = (span.index - self.degree()..=span.index)
                    .map(|i| {
                        let i = self.basis.pole_index(i);
                        spline::homogeneous::<1>(self.poles[i], self.weights[i]).map(|[x]| x)
                    })
                    .collect();
                let length = &span.end - &span.start;
                let lower = (&span.lower - &span.start) / &length;
                let upper = (&span.upper - &span.start) / length;
                ExactBezierCurve3 {
                    controls: spline::bezier_controls(
                        &self.basis,
                        span.index,
                        controls,
                        &lower,
                        &upper,
                    ),
                    domain: [span.lower, span.upper],
                }
            })
            .collect())
    }
}

impl BezierCurve3 {
    /// Preserve the input binary64 values as exact homogeneous rationals.
    pub fn to_exact(&self) -> ExactBezierCurve3 {
        ExactBezierCurve3 {
            controls: self
                .0
                .poles
                .iter()
                .zip(&self.0.weights)
                .map(|(&p, &w)| spline::homogeneous::<1>(p, w).map(|[x]| x))
                .collect(),
            domain: [integer(0), integer(1)],
        }
    }
}

/// Positive-weight rational Bézier geometry with an exact increasing domain.
/// Controls use `(w*x,w*y,w*z,w)` and local `t=(u-a)/(b-a)`. All public
/// parameters and derivatives use the original parameter `u`. Endpoint
/// evaluations are the limits of this arc, independent of adjacent arcs.
///
/// ```
/// use rusty_occt::{BezierCurve3, Point3};
/// use num_rational::BigRational;
/// let curve = BezierCurve3::new(vec![Point3::new(0.,0.,0.), Point3::new(1.,2.,0.)], None)?;
/// let third = BigRational::new(1.into(), 3.into());
/// let [left, right] = curve.to_exact().split_at(&third)?;
/// assert_eq!(left.domain()[1], right.domain()[0]);
/// assert_eq!(left.homogeneous_poles().last(), right.homogeneous_poles().first());
/// # Ok::<(), rusty_occt::Error>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactBezierCurve3 {
    controls: Vec<[R; 4]>,
    domain: [R; 2],
}

impl ExactBezierCurve3 {
    /// Internal construction from already validated positive homogeneous data.
    pub(crate) fn from_homogeneous(controls: Vec<[R; 4]>, domain: [R; 2]) -> Self {
        debug_assert!((2..=MAX_DEGREE + 1).contains(&controls.len()));
        debug_assert!(controls.iter().all(|p| p[3] > integer(0)));
        debug_assert!(domain[0] < domain[1]);
        Self { controls, domain }
    }

    pub fn degree(&self) -> usize {
        self.controls.len() - 1
    }
    pub fn domain(&self) -> &[R; 2] {
        &self.domain
    }
    pub fn homogeneous_poles(&self) -> &[[R; 4]] {
        &self.controls
    }

    /// Subdivide at a strictly interior parameter, with no float conversion.
    pub fn split_at(&self, u: &R) -> Result<[Self; 2]> {
        let u = parameter(u)?;
        if u <= self.domain[0] || u >= self.domain[1] {
            return Err(Error::OutOfDomain("Bezier split parameter"));
        }
        let t = (&u - &self.domain[0]) / (&self.domain[1] - &self.domain[0]);
        let mut row = IntegerRow::new(&self.controls);
        let mut left = vec![self.controls[0].clone()];
        let mut right = vec![self.controls.last().unwrap().clone()];
        while row.poles.len() > 1 {
            row.interpolate(&t);
            left.push(row.rational(0));
            right.push(row.rational(row.poles.len() - 1));
        }
        right.reverse();
        Ok([
            Self {
                controls: left,
                domain: [self.domain[0].clone(), u.clone()],
            },
            Self {
                controls: right,
                domain: [u, self.domain[1].clone()],
            },
        ])
    }

    /// Retain a positive-length closed subinterval, including the full domain.
    pub fn trim(&self, first: &R, last: &R) -> Result<Self> {
        let first = parameter(first)?;
        let last = parameter(last)?;
        if first < self.domain[0] || last > self.domain[1] || first >= last {
            return Err(Error::OutOfDomain("Bezier trim interval"));
        }
        let result = if last < self.domain[1] {
            let [left, _] = self.split_at(&last)?;
            left
        } else {
            self.clone()
        };
        Ok(if first > result.domain[0] {
            let [_, right] = result.split_at(&first)?;
            right
        } else {
            result
        })
    }

    /// `reversed(u) = original(a+b-u)` on the same increasing domain.
    pub fn reversed(&self) -> Self {
        let mut result = self.clone();
        result.controls.reverse();
        result
    }

    /// Raise the Bernstein degree without altering the rational function.
    /// Equal degree returns the identity; reduction and degree >25 are errors.
    pub fn elevated(&self, degree: usize) -> Result<Self> {
        if degree < self.degree() || degree > MAX_DEGREE {
            return Err(Error::InvalidCurve("Bezier elevation degree"));
        }
        if degree == self.degree() {
            return Ok(self.clone());
        }
        let mut row = IntegerRow::new(&self.controls);
        while row.poles.len() <= degree {
            let n = row.poles.len();
            let mut controls = Vec::with_capacity(n + 1);
            controls.push(std::array::from_fn(|c| &row.poles[0][c] * n));
            for i in 1..n {
                controls.push(std::array::from_fn(|c| {
                    &row.poles[i - 1][c] * i + &row.poles[i][c] * (n - i)
                }));
            }
            controls.push(std::array::from_fn(|c| &row.poles[n - 1][c] * n));
            row.poles = controls;
            row.denominator *= n;
        }
        Ok(Self {
            controls: (0..row.poles.len()).map(|i| row.rational(i)).collect(),
            domain: self.domain.clone(),
        })
    }

    /// Exact point and optional derivatives, retained even when their finite
    /// binary64 enclosures cannot be represented. Normalize rational inputs;
    /// reject a zero denominator or parameters outside the closed domain.
    pub fn exact_evaluate(&self, u: &R, order: DerivativeOrder) -> Result<ExactCurveEvaluation> {
        let u = parameter(u)?;
        if u < self.domain[0] || u > self.domain[1] {
            return Err(Error::OutOfDomain("Bezier evaluation parameter"));
        }
        let count = order.count();
        let length = &self.domain[1] - &self.domain[0];
        let t = (&u - &self.domain[0]) / &length;
        let jets = homogeneous_jet(&self.controls, &t, &length, count);
        let h = std::array::from_fn(|c| std::array::from_fn(|n| jets[n][c].clone()));
        let values = spline::rationalize(h, [(0, 0), (1, 0), (2, 0)], count + 1);
        Ok(ExactCurveEvaluation {
            position: ExactPoint3::from_coordinates(values[0].clone()),
            derivatives: [values.get(1).cloned(), values.get(2).cloned()],
        })
    }

    /// Convenience conversion with minimal finite binary64 enclosures.
    pub fn evaluate(&self, u: f64, order: DerivativeOrder) -> Result<CurveEvaluation> {
        crate::math::finite(u, "Bezier parameter")?;
        self.exact_evaluate(&spline::rational(u), order)?.enclosed()
    }
}

/// An exact rational jet. Conversion failure does not discard its values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactCurveEvaluation {
    position: ExactPoint3,
    derivatives: [Option<[R; 3]>; 2],
}
impl ExactCurveEvaluation {
    pub(crate) fn from_values(values: Vec<[R; 3]>) -> Self {
        Self {
            position: ExactPoint3::from_coordinates(values[0].clone()),
            derivatives: [values.get(1).cloned(), values.get(2).cloned()],
        }
    }
    pub fn position(&self) -> &ExactPoint3 {
        &self.position
    }
    pub fn derivative(&self, order: usize) -> Option<&[R; 3]> {
        match order {
            1 | 2 => self.derivatives[order - 1].as_ref(),
            _ => None,
        }
    }
    pub fn enclosed(&self) -> Result<CurveEvaluation> {
        Ok(CurveEvaluation {
            position: spline::bounds(self.position.coordinates())?,
            derivatives: [
                self.derivatives[0]
                    .as_ref()
                    .map(spline::bounds)
                    .transpose()?,
                self.derivatives[1]
                    .as_ref()
                    .map(spline::bounds)
                    .transpose()?,
            ],
        })
    }
}

/// A shared positive denominator lets each de Casteljau stage use integer
/// multiply/add operations. At stage r the denominator is D*t.denom()^r;
/// normalize only the requested boundary controls or final value. The exact
/// Bernstein recurrence and published reduced rational controls are unchanged.
struct IntegerRow {
    poles: Vec<[BigInt; 4]>,
    denominator: BigInt,
}
impl IntegerRow {
    fn new(controls: &[[R; 4]]) -> Self {
        let mut denominator = BigInt::from(1);
        for x in controls.iter().flatten() {
            let (mut a, mut b) = (denominator.clone(), x.denom().clone());
            while b != BigInt::from(0) {
                (a, b) = (b.clone(), a % b);
            }
            denominator = denominator / a * x.denom();
        }
        let poles = controls
            .iter()
            .map(|p| std::array::from_fn(|c| p[c].numer() * (&denominator / p[c].denom())))
            .collect();
        Self { poles, denominator }
    }
    fn rational(&self, i: usize) -> [R; 4] {
        std::array::from_fn(|c| R::new(self.poles[i][c].clone(), self.denominator.clone()))
    }
    fn interpolate(&mut self, t: &R) {
        let complement = t.denom() - t.numer();
        self.poles = self
            .poles
            .windows(2)
            .map(|p| std::array::from_fn(|c| &complement * &p[0][c] + t.numer() * &p[1][c]))
            .collect();
        self.denominator *= t.denom();
    }
}
pub(crate) fn homogeneous_value(controls: &[[R; 4]], t: &R) -> [R; 4] {
    if t == &integer(0) {
        return controls[0].clone();
    }
    if t == &integer(1) {
        return controls.last().unwrap().clone();
    }
    let mut row = IntegerRow::new(controls);
    while row.poles.len() > 1 {
        row.interpolate(t);
    }
    row.rational(0)
}

/// Homogeneous jets of a Bernstein row; the parameter is local but derivative
/// units use `length`. Also accepts signed/zero derivative control rows.
pub(crate) fn homogeneous_jet(controls: &[[R; 4]], t: &R, length: &R, order: usize) -> [[R; 4]; 3] {
    let mut jets = std::array::from_fn(|_| std::array::from_fn(|_| integer(0)));
    if order == 0 {
        jets[0] = homogeneous_value(controls, t);
        return jets;
    }
    let degree = controls.len() - 1;
    // The final two/three rows of one de Casteljau triangle give D1/D2.
    // This avoids separately forming and evaluating rational derivative grids.
    // At an endpoint those rows are just the first/last controls.
    let endpoint_count = controls.len().min(order + 1);
    let selected = if t == &integer(0) {
        &controls[..endpoint_count]
    } else if t == &integer(1) {
        &controls[controls.len() - endpoint_count..]
    } else {
        controls
    };
    let mut row = IntegerRow::new(selected);
    while row.poles.len() > 1 {
        if row.poles.len() == 3 && order == 2 {
            let denominator = &row.denominator * length.numer() * length.numer();
            let factor = length.denom() * length.denom() * degree * (degree - 1);
            jets[2] = std::array::from_fn(|c| {
                R::new(
                    (&row.poles[2][c] - &row.poles[1][c] * 2 + &row.poles[0][c]) * &factor,
                    denominator.clone(),
                )
            });
        }
        if row.poles.len() == 2 {
            let denominator = &row.denominator * length.numer();
            let factor = length.denom() * degree;
            jets[1] = std::array::from_fn(|c| {
                R::new(
                    (&row.poles[1][c] - &row.poles[0][c]) * &factor,
                    denominator.clone(),
                )
            });
        }
        row.interpolate(t);
    }
    jets[0] = row.rational(0);
    jets
}
