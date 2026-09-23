//! Certified rational B-spline and Bezier curve evaluation.
//!
//! OCCT references: Geom_BSplineCurve::CheckCurveData, BSplCLib::Eval/Bohm,
//! BSplCLib_D0/D1/D2, and PLib::RationalDerivative. Exact de Boor interpolation
//! with differentiated recurrences replaces floating evaluation. Every output
//! component has its smallest finite binary64 enclosure. See MATHEMATICS.md.
use crate::{exact, math::finite, spline, Error, KnotVector, Point3, Result, ScalarInterval};
use num_rational::BigRational as R;

mod editing;
mod knot_editing;
pub(crate) use editing::{homogeneous_jet, homogeneous_value};
pub use editing::{BezierExtractionOptions, ExactBezierCurve3, ExactCurveEvaluation};
pub use knot_editing::ExactBSplineCurve3;

pub use crate::spline::{MAX_DEGREE, MAX_POLES};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivativeOrder {
    Position,
    First,
    Second,
}
impl DerivativeOrder {
    pub(crate) fn count(self) -> usize {
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
    /// from both sides, including a periodic seam. At a nonperiodic endpoint,
    /// use the only interior side.
    Automatic,
    /// Left limit. At a periodic seam, use the end of the preceding period.
    /// Invalid at the first parameter of a nonperiodic domain.
    Left,
    /// Right limit. At a periodic seam, use the start of the following period.
    /// Invalid at the last parameter of a nonperiodic domain.
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

/// Immutable rational curve. Nonperiodic domains may be clamped or unclamped;
/// periodic curves use OCCT's cyclic control-point order.
#[derive(Debug, Clone)]
pub struct BSplineCurve3 {
    basis: KnotVector,
    poles: Vec<Point3>,
    weights: Vec<f64>,
}
impl BSplineCurve3 {
    /// Degree 1..=25; degree < pole count <= 4096. Interior multiplicities are
    /// 1..=degree, end multiplicities 1..=degree+1, sum(mults)=poles+degree+1.
    /// All coordinates/knots are finite; every supplied weight is positive and
    /// finite. No knot/weight is snapped, normalized or discarded by tolerance.
    pub fn new(
        degree: usize,
        poles: Vec<Point3>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<usize>,
    ) -> Result<Self> {
        Self::build(
            KnotVector::new(degree, knots, multiplicities),
            poles,
            weights,
        )
    }
    /// Same finite-data limits as `new`. Periodic end multiplicities must agree
    /// and are at most the degree. The pole count is sum(mults) minus one end
    /// multiplicity. Every finite parameter is accepted and wrapped exactly.
    pub fn new_periodic(
        degree: usize,
        poles: Vec<Point3>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<usize>,
    ) -> Result<Self> {
        Self::build(
            KnotVector::new_periodic(degree, knots, multiplicities),
            poles,
            weights,
        )
    }
    fn build(
        basis: Result<KnotVector>,
        poles: Vec<Point3>,
        weights: Option<Vec<f64>>,
    ) -> Result<Self> {
        if poles.len() > MAX_POLES {
            return Err(Error::LimitExceeded("spline control data"));
        }
        let basis = basis.map_err(|e| match e {
            Error::InvalidSpline(what) => Error::InvalidCurve(what),
            other => other,
        })?;
        if basis.pole_count() != poles.len() {
            return Err(Error::InvalidCurve("pole/knot/multiplicity counts"));
        }
        for p in &poles {
            exact::point(*p)?;
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
        Ok(Self {
            basis,
            poles,
            weights,
        })
    }
    pub fn degree(&self) -> usize {
        self.basis.degree()
    }
    pub fn poles(&self) -> &[Point3] {
        &self.poles
    }
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    pub fn knots(&self) -> &[f64] {
        self.basis.knots()
    }
    pub fn multiplicities(&self) -> &[usize] {
        self.basis.multiplicities()
    }
    pub fn knot_vector(&self) -> &KnotVector {
        &self.basis
    }
    pub fn is_rational(&self) -> bool {
        self.weights.windows(2).any(|w| w[0] != w[1])
    }
    pub fn is_periodic(&self) -> bool {
        self.basis.is_periodic()
    }
    pub fn domain(&self) -> (f64, f64) {
        self.basis.domain()
    }

    /// Evaluate within a nonperiodic closed domain without extrapolation, or at
    /// any finite periodic parameter. An automatic derivative at a knot/seam
    /// is rejected if the exact one-sided jets disagree. Unrepresentable
    /// requested components fail atomically; position-only queries do not
    /// compute potentially overflowing derivatives.
    pub fn evaluate(
        &self,
        u: f64,
        order: DerivativeOrder,
        side: KnotSide,
    ) -> Result<CurveEvaluation> {
        let count = order.count();
        let at = self.basis.locate(u, side, count)?;
        let values = self.exact_jet(&at[0], count);
        for other in &at[1..] {
            if values != self.exact_jet(other, count) {
                return Err(Error::DiscontinuousDerivative);
            }
        }
        Ok(CurveEvaluation {
            position: spline::bounds(&values[0])?,
            derivatives: [
                if count >= 1 {
                    Some(spline::bounds(&values[1])?)
                } else {
                    None
                },
                if count >= 2 {
                    Some(spline::bounds(&values[2])?)
                } else {
                    None
                },
            ],
        })
    }
    fn exact_jet(&self, at: &spline::Parameter, order: usize) -> Vec<[R; 3]> {
        let local = (at.span - self.degree()..=at.span)
            .map(|i| {
                let i = self.basis.pole_index(i);
                spline::homogeneous::<3>(self.poles[i], self.weights[i])
            })
            .collect();
        let h = spline::de_boor(
            &self.basis,
            at,
            local,
            [None, Some((0, 1)), Some((1, 2))],
            order + 1,
        );
        spline::rationalize(h, [(0, 0), (1, 0), (2, 0)], order + 1)
    }

    pub(crate) fn span_polynomial(&self, span: usize) -> [Vec<R>; 4] {
        let poles = (span - self.degree()..=span)
            .map(|i| {
                let i = self.basis.pole_index(i);
                spline::homogeneous::<1>(self.poles[i], self.weights[i]).map(|[x]| x)
            })
            .collect();
        spline::span_polynomial(&self.basis, span, poles)
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
