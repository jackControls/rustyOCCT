//! Exact homogeneous spline extraction and immutable Bernstein editing.
use super::{BSplineCurve3, BezierCurve3, CurveEvaluation, DerivativeOrder, MAX_DEGREE, MAX_POLES};
use crate::{intersection::ExactPoint3, spline, Error, Result};
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
        let mut row = self.controls.clone();
        let mut left = vec![row[0].clone()];
        let mut right = vec![row.last().unwrap().clone()];
        while row.len() > 1 {
            row = interpolate(&row, &t);
            left.push(row[0].clone());
            right.push(row.last().unwrap().clone());
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
            self.split_at(&last)?[0].clone()
        } else {
            self.clone()
        };
        Ok(if first > result.domain[0] {
            result.split_at(&first)?[1].clone()
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
        let mut result = self.clone();
        while result.degree() < degree {
            let n = result.controls.len();
            let mut controls = Vec::with_capacity(n + 1);
            controls.push(result.controls[0].clone());
            for i in 1..n {
                let alpha = integer(i) / integer(n);
                controls.push(blend(&result.controls[i], &result.controls[i - 1], &alpha));
            }
            controls.push(result.controls.last().unwrap().clone());
            result.controls = controls;
        }
        Ok(result)
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
        let mut controls = self.controls.clone();
        let mut jets: [[R; 4]; 3] = std::array::from_fn(|_| std::array::from_fn(|_| integer(0)));
        for (n, jet) in jets.iter_mut().enumerate().take(count + 1) {
            *jet = casteljau(&controls, &t);
            if n < count {
                let factor = integer(controls.len() - 1) / &length;
                controls = if controls.len() == 1 {
                    vec![std::array::from_fn(|_| integer(0))]
                } else {
                    controls
                        .windows(2)
                        .map(|p| std::array::from_fn(|c| &factor * (&p[1][c] - &p[0][c])))
                        .collect()
                };
            }
        }
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

fn blend(a: &[R; 4], b: &[R; 4], t: &R) -> [R; 4] {
    let complement = integer(1) - t;
    std::array::from_fn(|c| &complement * &a[c] + t * &b[c])
}
fn interpolate(row: &[[R; 4]], t: &R) -> Vec<[R; 4]> {
    row.windows(2).map(|p| blend(&p[0], &p[1], t)).collect()
}
fn casteljau(controls: &[[R; 4]], t: &R) -> [R; 4] {
    if t == &integer(0) {
        return controls[0].clone();
    }
    if t == &integer(1) {
        return controls.last().unwrap().clone();
    }
    let mut row = controls.to_vec();
    while row.len() > 1 {
        row = interpolate(&row, t);
    }
    row.pop().unwrap()
}
