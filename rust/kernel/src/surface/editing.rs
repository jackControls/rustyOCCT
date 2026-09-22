//! Exact tensor Bernstein editing. OCCT reference: GeomConvert surface-to-
//! Bezier conversion, Geom_BezierSurface editing, BSplSLib Iso/IncreaseDegree.
//! Shared exact curve operations act independently on homogeneous grid axes.
use super::{BSplineSurface3, BezierSurface3, SurfaceEvaluation, MAX_DEGREE, PARTIALS};
use crate::curve::{homogeneous_jet, homogeneous_value, DerivativeOrder, ExactBezierCurve3};
use crate::{intersection::ExactPoint3, spline, Error, Result};
use num_rational::BigRational as R;

fn integer(n: usize) -> R {
    R::from_integer(n.into())
}
fn parameter(value: &R) -> Result<R> {
    if value.denom() == &0.into() {
        return Err(Error::InvalidSurface(
            "rational parameter has zero denominator",
        ));
    }
    Ok(R::new(value.numer().clone(), value.denom().clone()))
}

/// Preflight limits for tensor extraction, including all periodic turns.
/// Both patch count and the sum of their control counts must fit. Zero permits
/// no results. These are allocation/work limits, not exact-arithmetic deadlines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BezierPatchExtractionOptions {
    pub max_patches: usize,
    pub max_controls: usize,
}
impl Default for BezierPatchExtractionOptions {
    fn default() -> Self {
        Self {
            max_patches: 4096,
            max_controls: 1_048_576,
        }
    }
}

impl BSplineSurface3 {
    /// Exact patches over the fundamental domain, ordered by U then V.
    pub fn bezier_patches(&self) -> Result<Vec<ExactBezierSurface3>> {
        let ((ua, ub), (va, vb)) = self.domain();
        self.bezier_patches_in(ua, ub, va, vb)
    }

    /// Extract a positive closed rectangle. Each periodic axis may cover
    /// multiple turns; nonperiodic axes must remain inside their domain.
    pub fn bezier_patches_in(
        &self,
        ua: f64,
        ub: f64,
        va: f64,
        vb: f64,
    ) -> Result<Vec<ExactBezierSurface3>> {
        self.bezier_patches_with_options(ua, ub, va, vb, BezierPatchExtractionOptions::default())
    }

    pub fn bezier_patches_with_options(
        &self,
        ua: f64,
        ub: f64,
        va: f64,
        vb: f64,
        options: BezierPatchExtractionOptions,
    ) -> Result<Vec<ExactBezierSurface3>> {
        let (du, dv) = (self.u.degree(), self.v.degree());
        let controls_per_patch = (du + 1) * (dv + 1);
        let limit = options
            .max_patches
            .min(options.max_controls / controls_per_patch);
        let us = self.u.spans_in(ua, ub, limit)?;
        // Every valid positive interval contains at least one span. Restrict
        // the second traversal before enumeration to bound the Cartesian product.
        let vs = self.v.spans_in(va, vb, limit / us.len())?;
        let count = us
            .len()
            .checked_mul(vs.len())
            .ok_or(Error::ComputationLimit("Bezier patch count"))?;
        let mut result = Vec::with_capacity(count);
        let v_maps: Vec<_> = vs
            .iter()
            .map(|v| {
                let length = &v.end - &v.start;
                spline::BezierSpanTransform::new(
                    &self.v,
                    v.index,
                    &((&v.lower - &v.start) / &length),
                    &((&v.upper - &v.start) / length),
                )
            })
            .collect();
        for u in &us {
            let ulength = &u.end - &u.start;
            let ulow = (&u.lower - &u.start) / &ulength;
            let uhigh = (&u.upper - &u.start) / ulength;
            let u_map = spline::BezierSpanTransform::new(&self.u, u.index, &ulow, &uhigh);
            for (v, v_map) in vs.iter().zip(&v_maps) {
                let columns: Vec<_> = (v.index - dv..=v.index)
                    .map(|j| {
                        let poles: Vec<_> = (u.index - du..=u.index)
                            .map(|i| {
                                let index = self.u.pole_index(i) * self.v.pole_count()
                                    + self.v.pole_index(j);
                                spline::homogeneous::<1>(self.poles[index], self.weights[index])
                                    .map(|[x]| x)
                            })
                            .collect();
                        u_map.apply(&poles)
                    })
                    .collect();
                let controls = (0..=du)
                    .flat_map(|i| {
                        let row: Vec<_> = columns.iter().map(|c| c[i].clone()).collect();
                        v_map.apply(&row)
                    })
                    .collect();
                result.push(ExactBezierSurface3 {
                    degrees: [du, dv],
                    controls,
                    domain: [
                        [u.lower.clone(), u.upper.clone()],
                        [v.lower.clone(), v.upper.clone()],
                    ],
                });
            }
        }
        Ok(result)
    }
}

impl BezierSurface3 {
    /// Exact homogeneous input data on [0,1]², without rounding products.
    pub fn to_exact(&self) -> ExactBezierSurface3 {
        ExactBezierSurface3 {
            degrees: [self.0.u.degree(), self.0.v.degree()],
            controls: self
                .0
                .poles
                .iter()
                .zip(&self.0.weights)
                .map(|(&p, &w)| spline::homogeneous::<1>(p, w).map(|[x]| x))
                .collect(),
            domain: [[integer(0), integer(1)], [integer(0), integer(1)]],
        }
    }
}

/// Positive-weight rational tensor Bézier patch. U-major controls are
/// `(wx,wy,wz,w)`; parameters and derivatives use the retained original U/V
/// domains. Each patch owns its one-sided boundary limits. This geometry does
/// not imply face topology, arbitrary trim loops or continuity with neighbors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactBezierSurface3 {
    degrees: [usize; 2],
    controls: Vec<[R; 4]>,
    domain: [[R; 2]; 2],
}
impl ExactBezierSurface3 {
    pub fn degrees(&self) -> [usize; 2] {
        self.degrees
    }
    pub fn domain(&self) -> &[[R; 2]; 2] {
        &self.domain
    }
    pub fn homogeneous_poles(&self) -> &[[R; 4]] {
        &self.controls
    }

    fn rows(&self, axis: usize) -> Vec<ExactBezierCurve3> {
        let other = 1 - axis;
        (0..=self.degrees[other])
            .map(|fixed| {
                let controls = (0..=self.degrees[axis])
                    .map(|i| {
                        let index = if axis == 0 {
                            i * (self.degrees[1] + 1) + fixed
                        } else {
                            fixed * (self.degrees[1] + 1) + i
                        };
                        self.controls[index].clone()
                    })
                    .collect();
                ExactBezierCurve3::from_homogeneous(controls, self.domain[axis].clone())
            })
            .collect()
    }
    fn with_rows(&self, axis: usize, rows: &[ExactBezierCurve3]) -> Self {
        let mut degrees = self.degrees;
        let mut domain = self.domain.clone();
        degrees[axis] = rows[0].degree();
        domain[axis] = rows[0].domain().clone();
        let controls = (0..=degrees[0])
            .flat_map(|i| {
                (0..=degrees[1]).map(move |j| {
                    let (row, col) = if axis == 0 { (j, i) } else { (i, j) };
                    rows[row].homogeneous_poles()[col].clone()
                })
            })
            .collect();
        Self {
            degrees,
            domain,
            controls,
        }
    }
    fn split(&self, axis: usize, value: &R) -> Result<[Self; 2]> {
        let value = parameter(value)?;
        if value <= self.domain[axis][0] || value >= self.domain[axis][1] {
            return Err(Error::OutOfDomain("Bezier surface split parameter"));
        }
        let mut left = Vec::new();
        let mut right = Vec::new();
        for row in self.rows(axis) {
            let [a, b] = row.split_at(&value)?;
            left.push(a);
            right.push(b);
        }
        Ok([self.with_rows(axis, &left), self.with_rows(axis, &right)])
    }
    pub fn split_u_at(&self, u: &R) -> Result<[Self; 2]> {
        self.split(0, u)
    }
    pub fn split_v_at(&self, v: &R) -> Result<[Self; 2]> {
        self.split(1, v)
    }

    /// Retain a positive closed rectangle, including the full domain.
    pub fn trim(&self, ua: &R, ub: &R, va: &R, vb: &R) -> Result<Self> {
        let ranges = [
            [parameter(ua)?, parameter(ub)?],
            [parameter(va)?, parameter(vb)?],
        ];
        for (range, domain) in ranges.iter().zip(&self.domain) {
            if range[0] < domain[0] || range[1] > domain[1] || range[0] >= range[1] {
                return Err(Error::OutOfDomain("Bezier surface trim rectangle"));
            }
        }
        let mut result = self.clone();
        for (axis, range) in ranges.iter().enumerate() {
            if range == &result.domain[axis] {
                continue;
            }
            let rows = result
                .rows(axis)
                .iter()
                .map(|r| r.trim(&range[0], &range[1]))
                .collect::<Result<Vec<_>>>()?;
            result = result.with_rows(axis, &rows);
        }
        Ok(result)
    }

    fn reversed(&self, axis: usize) -> Self {
        self.with_rows(
            axis,
            &self
                .rows(axis)
                .iter()
                .map(ExactBezierCurve3::reversed)
                .collect::<Vec<_>>(),
        )
    }
    /// `result(u,v) = original(a+b-u,v)` on the same increasing U domain.
    pub fn u_reversed(&self) -> Self {
        self.reversed(0)
    }
    /// `result(u,v) = original(u,c+d-v)` on the same increasing V domain.
    pub fn v_reversed(&self) -> Self {
        self.reversed(1)
    }

    /// Transpose controls, degrees and domains: `result(v,u) = original(u,v)`.
    /// The corresponding parameter-space orientation is reversed.
    pub fn exchanged_uv(&self) -> Self {
        Self {
            degrees: [self.degrees[1], self.degrees[0]],
            domain: [self.domain[1].clone(), self.domain[0].clone()],
            controls: (0..=self.degrees[1])
                .flat_map(|j| {
                    (0..=self.degrees[0])
                        .map(move |i| self.controls[i * (self.degrees[1] + 1) + j].clone())
                })
                .collect(),
        }
    }

    /// Independently raise both degrees without changing the homogeneous
    /// tensor polynomial. Reduction or a degree above 25 is rejected.
    pub fn elevated(&self, u_degree: usize, v_degree: usize) -> Result<Self> {
        let target = [u_degree, v_degree];
        if target
            .iter()
            .zip(self.degrees)
            .any(|(&new, old)| new < old || new > MAX_DEGREE)
        {
            return Err(Error::InvalidSurface("Bezier surface elevation degrees"));
        }
        let mut result = self.clone();
        for (axis, &degree) in target.iter().enumerate() {
            if degree == result.degrees[axis] {
                continue;
            }
            let rows = result
                .rows(axis)
                .iter()
                .map(|r| r.elevated(degree))
                .collect::<Result<Vec<_>>>()?;
            result = result.with_rows(axis, &rows);
        }
        Ok(result)
    }

    fn iso(&self, axis: usize, value: &R) -> Result<ExactBezierCurve3> {
        let value = parameter(value)?;
        let [a, b] = &self.domain[axis];
        if &value < a || &value > b {
            return Err(Error::OutOfDomain("Bezier surface isoparameter"));
        }
        let t = (&value - a) / (b - a);
        let controls = self
            .rows(axis)
            .iter()
            .map(|r| homogeneous_value(r.homogeneous_poles(), &t))
            .collect();
        Ok(ExactBezierCurve3::from_homogeneous(
            controls,
            self.domain[1 - axis].clone(),
        ))
    }
    /// Constant U, varying V in the original V domain, including boundaries.
    pub fn u_iso(&self, u: &R) -> Result<ExactBezierCurve3> {
        self.iso(0, u)
    }
    /// Constant V, varying U in the original U domain, including boundaries.
    pub fn v_iso(&self, v: &R) -> Result<ExactBezierCurve3> {
        self.iso(1, v)
    }

    /// Exact position and every partial through the requested total order.
    /// Conversion failure does not discard these retained rational values.
    pub fn exact_evaluate(
        &self,
        u: &R,
        v: &R,
        order: DerivativeOrder,
    ) -> Result<ExactSurfaceEvaluation> {
        let parameters = [parameter(u)?, parameter(v)?];
        for (p, d) in parameters.iter().zip(&self.domain) {
            if p < &d[0] || p > &d[1] {
                return Err(Error::OutOfDomain("Bezier surface evaluation parameter"));
            }
        }
        let lengths: [R; 2] = std::array::from_fn(|i| &self.domain[i][1] - &self.domain[i][0]);
        let t: [R; 2] =
            std::array::from_fn(|i| (&parameters[i] - &self.domain[i][0]) / &lengths[i]);
        let order = order.count();
        let count = [1, 3, 6][order];
        let v_jets: Vec<_> = self
            .controls
            .chunks(self.degrees[1] + 1)
            .map(|row| homogeneous_jet(row, &t[1], &lengths[1], order))
            .collect();
        let mut jets: [[R; 4]; 6] = std::array::from_fn(|_| std::array::from_fn(|_| integer(0)));
        for v_order in 0..=order {
            let controls: Vec<_> = v_jets.iter().map(|row| row[v_order].clone()).collect();
            let u_jets = homogeneous_jet(&controls, &t[0], &lengths[0], order - v_order);
            for (i, &(u, v)) in PARTIALS.iter().enumerate().take(count) {
                if v == v_order {
                    jets[i] = u_jets[u].clone();
                }
            }
        }
        let h = std::array::from_fn(|c| std::array::from_fn(|i| jets[i][c].clone()));
        let values = spline::rationalize(h, PARTIALS, count);
        Ok(ExactSurfaceEvaluation {
            position: ExactPoint3::from_coordinates(values[0].clone()),
            derivatives: std::array::from_fn(|i| values.get(i + 1).cloned()),
        })
    }

    /// Convenience binary64 enclosures, failing atomically on overflow.
    pub fn evaluate(&self, u: f64, v: f64, order: DerivativeOrder) -> Result<SurfaceEvaluation> {
        crate::math::finite(u, "Bezier surface U parameter")?;
        crate::math::finite(v, "Bezier surface V parameter")?;
        self.exact_evaluate(&spline::rational(u), &spline::rational(v), order)?
            .enclosed()
    }
}

/// A retained exact surface jet. Only requested derivatives are populated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactSurfaceEvaluation {
    position: ExactPoint3,
    derivatives: [Option<[R; 3]>; 5],
}
impl ExactSurfaceEvaluation {
    pub fn position(&self) -> &ExactPoint3 {
        &self.position
    }
    pub fn derivative(&self, u_order: usize, v_order: usize) -> Option<&[R; 3]> {
        PARTIALS
            .iter()
            .skip(1)
            .position(|&p| p == (u_order, v_order))
            .and_then(|i| self.derivatives[i].as_ref())
    }
    pub fn enclosed(&self) -> Result<SurfaceEvaluation> {
        let mut derivatives = [None; 5];
        for (actual, exact) in derivatives.iter_mut().zip(&self.derivatives) {
            *actual = exact.as_ref().map(spline::bounds).transpose()?;
        }
        Ok(SurfaceEvaluation {
            position: spline::bounds(self.position.coordinates())?,
            derivatives,
        })
    }
}
