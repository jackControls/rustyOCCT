//! Certified rational tensor-product Bezier and B-spline surface evaluation.
//!
//! OCCT references: Geom_BSplineSurface::CheckSurfaceData and LocalD0/D1/D2,
//! BSplSLib::PrepareEval, D0/D1/D2 and RationalDerivative. Exact differentiated
//! homogeneous de Boor interpolation is applied in each parameter direction.
use crate::curve::{DerivativeOrder, KnotSide};
use crate::spline::{self, Parameter, MAX_DEGREE, MAX_POLES};
use crate::{exact, math::finite, Error, KnotVector, Point3, Result, ScalarInterval};
use num_rational::BigRational as R;

mod editing;
pub use editing::{BezierPatchExtractionOptions, ExactBezierSurface3, ExactSurfaceEvaluation};

// Ordered by total degree. Every lower partial precedes its dependents.
const PARTIALS: [(usize, usize); 6] = [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1)];
const U_PREVIOUS: [Option<(usize, usize)>; 6] =
    [None, Some((0, 1)), None, Some((1, 2)), None, Some((2, 1))];
const V_PREVIOUS: [Option<(usize, usize)>; 6] =
    [None, None, Some((0, 1)), None, Some((2, 2)), Some((1, 1))];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceEvaluation {
    position: [ScalarInterval; 3],
    derivatives: [Option<[ScalarInterval; 3]>; 5],
}
impl SurfaceEvaluation {
    pub fn position_bounds(self) -> [ScalarInterval; 3] {
        self.position
    }
    /// Rounded representative; keep the bounds for subsequent decisions.
    pub fn position(self) -> Point3 {
        let p = self.position.map(ScalarInterval::representative);
        Point3::new(p[0], p[1], p[2])
    }
    /// Bounds for Du, Dv, Duu, Dvv or Duv. Returns None for other orders or if
    /// the derivative was not requested. Units use the original U/V parameters.
    pub fn derivative_bounds(self, u_order: usize, v_order: usize) -> Option<[ScalarInterval; 3]> {
        PARTIALS
            .iter()
            .skip(1)
            .position(|&x| x == (u_order, v_order))
            .and_then(|i| self.derivatives[i])
    }
}

/// Immutable positive-weight tensor-product surface. U and V may independently
/// be periodic. Control points and weights are U-major: index = u * v_count + v.
/// At most 4096 total poles, degrees 1..=25; no trimming or topology is implied.
#[derive(Debug, Clone)]
pub struct BSplineSurface3 {
    u: KnotVector,
    v: KnotVector,
    poles: Vec<Point3>,
    weights: Vec<f64>,
}
impl BSplineSurface3 {
    pub fn new(
        u: KnotVector,
        v: KnotVector,
        poles: Vec<Point3>,
        weights: Option<Vec<f64>>,
    ) -> Result<Self> {
        let count = u
            .pole_count()
            .checked_mul(v.pole_count())
            .ok_or(Error::LimitExceeded("spline surface control data"))?;
        if count > MAX_POLES || poles.len() > MAX_POLES {
            return Err(Error::LimitExceeded("spline surface control data"));
        }
        if poles.len() != count {
            return Err(Error::InvalidSurface("control grid dimensions"));
        }
        for p in &poles {
            exact::point(*p)?;
        }
        let weights = weights.unwrap_or_else(|| vec![1.; count]);
        if weights.len() != count {
            return Err(Error::InvalidSurface("weight grid dimensions"));
        }
        for &w in &weights {
            finite(w, "spline surface weight")?;
            if w <= 0. {
                return Err(Error::InvalidSurface("weights must be positive"));
            }
        }
        Ok(Self {
            u,
            v,
            poles,
            weights,
        })
    }
    pub fn u_knots(&self) -> &KnotVector {
        &self.u
    }
    pub fn v_knots(&self) -> &KnotVector {
        &self.v
    }
    pub fn poles(&self) -> &[Point3] {
        &self.poles
    }
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
    pub fn is_rational(&self) -> bool {
        self.weights.windows(2).any(|w| w[0] != w[1])
    }
    pub fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (self.u.domain(), self.v.domain())
    }

    /// Position, or position and every partial derivative through the requested
    /// total order (at most two). Independent U/V sides select a quadrant at
    /// intersecting knots. Automatic sides must agree exactly across all
    /// relevant quadrants, including periodic seams. No nonperiodic extrapolation.
    /// Any unrepresentable requested component fails the complete query.
    pub fn evaluate(
        &self,
        u: f64,
        v: f64,
        order: DerivativeOrder,
        sides: [KnotSide; 2],
    ) -> Result<SurfaceEvaluation> {
        let order = order.count();
        let count = [1, 3, 6][order];
        let us = self.u.locate(u, sides[0], order)?;
        let vs = self.v.locate(v, sides[1], order)?;
        let values = self.exact_jet(&us[0], &vs[0], count);
        for (i, at_u) in us.iter().enumerate() {
            for (j, at_v) in vs.iter().enumerate() {
                if (i != 0 || j != 0) && values != self.exact_jet(at_u, at_v, count) {
                    return Err(Error::DiscontinuousDerivative);
                }
            }
        }
        let mut derivatives = [None; 5];
        for i in 1..count {
            derivatives[i - 1] = Some(spline::bounds(&values[i])?);
        }
        Ok(SurfaceEvaluation {
            position: spline::bounds(&values[0])?,
            derivatives,
        })
    }
    fn exact_jet(&self, u: &Parameter, v: &Parameter, count: usize) -> Vec<[R; 3]> {
        let local = (u.span - self.u.degree()..=u.span)
            .map(|i| {
                let i = self.u.pole_index(i);
                let row = (v.span - self.v.degree()..=v.span)
                    .map(|j| {
                        let index = i * self.v.pole_count() + self.v.pole_index(j);
                        spline::homogeneous::<6>(self.poles[index], self.weights[index])
                    })
                    .collect();
                spline::de_boor(&self.v, v, row, V_PREVIOUS, count)
            })
            .collect();
        let h = spline::de_boor(&self.u, u, local, U_PREVIOUS, count);
        spline::rationalize(h, PARTIALS, count)
    }
}

/// Rational Bezier patch on [0,1] x [0,1], with U-major control data.
#[derive(Debug, Clone)]
pub struct BezierSurface3(BSplineSurface3);
impl BezierSurface3 {
    /// `(u_degree+1)*(v_degree+1)` poles and optional positive weights.
    pub fn new(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Point3>,
        weights: Option<Vec<f64>>,
    ) -> Result<Self> {
        if !(1..=MAX_DEGREE).contains(&u_degree) || !(1..=MAX_DEGREE).contains(&v_degree) {
            return Err(Error::InvalidSurface("Bezier degrees must be in 1..=25"));
        }
        Ok(Self(BSplineSurface3::new(
            KnotVector::new(u_degree, vec![0., 1.], vec![u_degree + 1; 2])?,
            KnotVector::new(v_degree, vec![0., 1.], vec![v_degree + 1; 2])?,
            poles,
            weights,
        )?))
    }
    pub fn as_bspline(&self) -> &BSplineSurface3 {
        &self.0
    }
    pub fn evaluate(&self, u: f64, v: f64, order: DerivativeOrder) -> Result<SurfaceEvaluation> {
        self.0.evaluate(u, v, order, [KnotSide::Automatic; 2])
    }
}
