//! Exact tensor B-spline editing. OCCT references: Geom_BSplineSurface U/V
//! insertion/removal and BSplSLib's homogeneous transverse-curve packing.
//! The existing exact curve algorithms preserve each complete tensor axis.
use super::{
    BSplineSurface3, BezierPatchExtractionOptions, ExactBezierSurface3, ExactSurfaceEvaluation,
    SurfaceEvaluation, PARTIALS,
};
use crate::curve::{DerivativeOrder, KnotSide};
use crate::spline::{self, Parameter, MAX_POLES};
use crate::{Error, ExactBSplineCurve3, ExactKnotVector, Result};
use num_rational::BigRational as R;

fn integer(n: usize) -> R {
    R::from_integer(n.into())
}

fn control_count(axes: &[ExactKnotVector; 2]) -> Result<usize> {
    let count = axes[0]
        .pole_count()
        .checked_mul(axes[1].pole_count())
        .ok_or(Error::LimitExceeded("spline surface control data"))?;
    if count > MAX_POLES {
        return Err(Error::LimitExceeded("spline surface control data"));
    }
    Ok(count)
}

/// Exact rational tensor B-spline with independent U/V periodicity. Controls
/// `(wx,wy,wz,w)` are U-major and weights are strictly positive. No face trim
/// loops or B-rep topology are implied. At most 4096 complete grid controls.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactBSplineSurface3 {
    axes: [ExactKnotVector; 2],
    controls: Vec<[R; 4]>,
}

impl BSplineSurface3 {
    /// Preserve every represented input atom, without rounding homogeneous products.
    pub fn to_exact(&self) -> ExactBSplineSurface3 {
        ExactBSplineSurface3 {
            axes: [self.u.to_exact(), self.v.to_exact()],
            controls: self
                .poles
                .iter()
                .zip(&self.weights)
                .map(|(&p, &w)| spline::homogeneous::<1>(p, w).map(|[x]| x))
                .collect(),
        }
    }
}

impl ExactBezierSurface3 {
    /// An exact clamped tensor B-spline on the same rational rectangle.
    pub fn to_bspline(&self) -> ExactBSplineSurface3 {
        let degrees = self.degrees();
        ExactBSplineSurface3 {
            axes: std::array::from_fn(|i| {
                ExactKnotVector::new(
                    degrees[i],
                    self.domain()[i].to_vec(),
                    vec![degrees[i] + 1; 2],
                )
                .expect("validated Bezier surface domain")
            }),
            controls: self.homogeneous_poles().to_vec(),
        }
    }
}

impl ExactBSplineSurface3 {
    /// Normalize all homogeneous rationals; reject zero denominators, nonpositive
    /// weights, grid mismatches and Cartesian control counts above 4096.
    pub fn from_homogeneous(
        u: ExactKnotVector,
        v: ExactKnotVector,
        controls: Vec<[R; 4]>,
    ) -> Result<Self> {
        let axes = [u, v];
        if controls.len() != control_count(&axes)? {
            return Err(Error::InvalidSurface("control grid dimensions"));
        }
        let controls: Vec<_> = controls
            .into_iter()
            .map(|[x, y, z, w]| {
                Ok([
                    spline::normalize(&x)?,
                    spline::normalize(&y)?,
                    spline::normalize(&z)?,
                    spline::normalize(&w)?,
                ])
            })
            .collect::<Result<_>>()?;
        if controls.iter().any(|p| p[3] <= integer(0)) {
            return Err(Error::InvalidSurface("weights must be positive"));
        }
        Ok(Self { axes, controls })
    }
    pub fn u_knots(&self) -> &ExactKnotVector {
        &self.axes[0]
    }
    pub fn v_knots(&self) -> &ExactKnotVector {
        &self.axes[1]
    }
    pub fn degrees(&self) -> [usize; 2] {
        self.axes.each_ref().map(|a| a.degree())
    }
    pub fn pole_counts(&self) -> [usize; 2] {
        self.axes.each_ref().map(|a| a.pole_count())
    }
    pub fn domain(&self) -> [[R; 2]; 2] {
        self.axes.each_ref().map(|a| a.domain().clone())
    }
    pub fn homogeneous_poles(&self) -> &[[R; 4]] {
        &self.controls
    }

    fn rows(&self, axis: usize) -> Vec<ExactBSplineCurve3> {
        let counts = self.pole_counts();
        (0..counts[1 - axis])
            .map(|fixed| {
                let controls = (0..counts[axis])
                    .map(|i| {
                        self.controls[if axis == 0 {
                            i * counts[1] + fixed
                        } else {
                            fixed * counts[1] + i
                        }]
                        .clone()
                    })
                    .collect();
                ExactBSplineCurve3::from_homogeneous(self.axes[axis].clone(), controls)
                    .expect("validated homogeneous grid row")
            })
            .collect()
    }
    fn with_rows(&self, axis: usize, rows: &[ExactBSplineCurve3]) -> Self {
        let mut axes = self.axes.clone();
        axes[axis] = rows[0].knot_vector().clone();
        debug_assert!(rows.iter().all(|r| r.knot_vector() == &axes[axis]));
        debug_assert!(control_count(&axes).is_ok());
        let controls = (0..axes[0].pole_count())
            .flat_map(|i| {
                (0..axes[1].pole_count()).map(move |j| {
                    let (row, col) = if axis == 0 { (j, i) } else { (i, j) };
                    rows[row].homogeneous_poles()[col].clone()
                })
            })
            .collect();
        Self { axes, controls }
    }

    /// Atomic U/V refinement by total multiplicities. Both request lists and
    /// the final Cartesian grid are checked before any control arithmetic.
    /// At most 4096 combined requests; no parameter snapping or approximation.
    pub fn refined(&self, u_requests: &[(R, usize)], v_requests: &[(R, usize)]) -> Result<Self> {
        if u_requests.len().saturating_add(v_requests.len()) > MAX_POLES {
            return Err(Error::LimitExceeded("surface knot refinement requests"));
        }
        let plans = [
            self.axes[0].refinement_plan(u_requests)?,
            self.axes[1].refinement_plan(v_requests)?,
        ];
        control_count(&[plans[0].0.clone(), plans[1].0.clone()])?;
        let mut result = self.clone();
        for axis in 0..2 {
            if plans[axis].1.is_empty() {
                continue;
            }
            let transform = spline::KnotRefinementTransform::new(
                &result.axes[axis],
                &plans[axis].0,
                &plans[axis].1,
            );
            let counts = result.pole_counts();
            let rows: Vec<_> = (0..counts[1 - axis])
                .map(|fixed| {
                    let controls: Vec<_> = (0..counts[axis])
                        .map(|i| {
                            result.controls[if axis == 0 {
                                i * counts[1] + fixed
                            } else {
                                fixed * counts[1] + i
                            }]
                            .clone()
                        })
                        .collect();
                    transform.apply(&controls)
                })
                .collect();
            let mut axes = result.axes.clone();
            axes[axis] = plans[axis].0.clone();
            let controls = (0..axes[0].pole_count())
                .flat_map(|i| {
                    let rows = &rows;
                    (0..axes[1].pole_count()).map(move |j| {
                        if axis == 0 {
                            rows[j][i].clone()
                        } else {
                            rows[i][j].clone()
                        }
                    })
                })
                .collect();
            result = Self { axes, controls };
        }
        Ok(result)
    }
    pub fn refined_u(&self, requests: &[(R, usize)]) -> Result<Self> {
        self.refined(requests, &[])
    }
    pub fn refined_v(&self, requests: &[(R, usize)]) -> Result<Self> {
        self.refined(&[], requests)
    }
    pub fn insert_u_knot(&self, u: &R, target: usize) -> Result<Self> {
        self.refined_u(&[(u.clone(), target)])
    }
    pub fn insert_v_knot(&self, v: &R, target: usize) -> Result<Self> {
        self.refined_v(&[(v.clone(), target)])
    }

    fn removed(&self, axis: usize, parameter: &R, target: usize) -> Result<Option<Self>> {
        let mut result = Vec::new();
        for row in self.rows(axis) {
            match row.remove_knot(parameter, target)? {
                Some(row) => result.push(row),
                None => return Ok(None),
            }
        }
        Ok(Some(self.with_rows(axis, &result)))
    }
    /// Preserve every homogeneous row exactly or return None. Invalid requests
    /// are errors. Periodic seam deletion can advance the U origin.
    pub fn remove_u_knot(&self, u: &R, target: usize) -> Result<Option<Self>> {
        self.removed(0, u, target)
    }
    /// Preserve every homogeneous column exactly or return None. Invalid requests
    /// are errors. Periodic seam deletion can advance the V origin.
    pub fn remove_v_knot(&self, v: &R, target: usize) -> Result<Option<Self>> {
        self.removed(1, v, target)
    }

    fn jet(&self, u: &Parameter, v: &Parameter, count: usize) -> Vec<[R; 3]> {
        let [du, dv] = self.degrees();
        let order = match count {
            1 => 0,
            3 => 1,
            _ => 2,
        };
        let um = spline::SplineJetTransform::new(du, &self.axes[0].flat, u, order);
        let vm = spline::SplineJetTransform::new(dv, &self.axes[1].flat, v, order);
        let v_jets: Vec<_> = (u.span - du..=u.span)
            .map(|i| {
                let poles: Vec<_> = (v.span - dv..=v.span)
                    .map(|j| {
                        self.controls[self.axes[0].pole_index(i) * self.axes[1].pole_count()
                            + self.axes[1].pole_index(j)]
                        .clone()
                    })
                    .collect();
                vm.apply(&poles)
            })
            .collect();
        let mut h: [[R; 6]; 4] = std::array::from_fn(|_| std::array::from_fn(|_| integer(0)));
        for v_order in 0..=order {
            let poles: Vec<_> = v_jets.iter().map(|row| row[v_order].clone()).collect();
            let jets = um.apply(&poles);
            for (n, &(a, b)) in PARTIALS.iter().take(count).enumerate() {
                if b == v_order {
                    for c in 0..4 {
                        h[c][n] = jets[a][c].clone();
                    }
                }
            }
        }
        spline::rationalize(h, PARTIALS, count)
    }
    /// Exact position and all requested partials in original parameter units.
    /// Automatic sides require equality across every relevant knot quadrant.
    pub fn exact_evaluate(
        &self,
        u: &R,
        v: &R,
        order: DerivativeOrder,
        sides: [KnotSide; 2],
    ) -> Result<ExactSurfaceEvaluation> {
        let order = order.count();
        let count = [1, 3, 6][order];
        let us = self.axes[0].locate(u, sides[0], order)?;
        let vs = self.axes[1].locate(v, sides[1], order)?;
        let values = self.jet(&us[0], &vs[0], count);
        for (i, at_u) in us.iter().enumerate() {
            for (j, at_v) in vs.iter().enumerate() {
                if (i != 0 || j != 0) && values != self.jet(at_u, at_v, count) {
                    return Err(Error::DiscontinuousDerivative);
                }
            }
        }
        Ok(ExactSurfaceEvaluation::from_values(values))
    }
    pub fn evaluate(
        &self,
        u: f64,
        v: f64,
        order: DerivativeOrder,
        sides: [KnotSide; 2],
    ) -> Result<SurfaceEvaluation> {
        crate::math::finite(u, "surface U parameter")?;
        crate::math::finite(v, "surface V parameter")?;
        self.exact_evaluate(&spline::rational(u), &spline::rational(v), order, sides)?
            .enclosed()
    }

    fn iso(&self, axis: usize, parameter: &R) -> Result<ExactBSplineCurve3> {
        let at = self.axes[axis]
            .locate(parameter, KnotSide::Automatic, 0)?
            .remove(0);
        let counts = self.pole_counts();
        let degree = self.axes[axis].degree();
        let map = spline::SplineJetTransform::new(degree, &self.axes[axis].flat, &at, 0);
        let controls = (0..counts[1 - axis])
            .map(|fixed| {
                let poles = (at.span - degree..=at.span)
                    .map(|i| {
                        let i = self.axes[axis].pole_index(i);
                        let index = if axis == 0 {
                            i * counts[1] + fixed
                        } else {
                            fixed * counts[1] + i
                        };
                        self.controls[index].clone()
                    })
                    .collect::<Vec<_>>();
                map.apply(&poles).remove(0)
            })
            .collect();
        ExactBSplineCurve3::from_homogeneous(self.axes[1 - axis].clone(), controls)
    }
    /// Fix U, retaining the complete exact V B-spline and its parameterization.
    pub fn u_iso(&self, u: &R) -> Result<ExactBSplineCurve3> {
        self.iso(0, u)
    }
    /// Fix V, retaining the complete exact U B-spline and its parameterization.
    pub fn v_iso(&self, v: &R) -> Result<ExactBSplineCurve3> {
        self.iso(1, v)
    }

    /// Exchange axes, transposing the complete U-major homogeneous grid.
    pub fn exchanged_uv(&self) -> Self {
        let [nu, nv] = self.pole_counts();
        Self {
            axes: [self.axes[1].clone(), self.axes[0].clone()],
            controls: (0..nv)
                .flat_map(|j| (0..nu).map(move |i| self.controls[i * nv + j].clone()))
                .collect(),
        }
    }

    pub fn bezier_patches(&self) -> Result<Vec<ExactBezierSurface3>> {
        let [[ua, ub], [va, vb]] = self.domain();
        self.bezier_patches_in(&ua, &ub, &va, &vb)
    }
    pub fn bezier_patches_in(
        &self,
        ua: &R,
        ub: &R,
        va: &R,
        vb: &R,
    ) -> Result<Vec<ExactBezierSurface3>> {
        self.bezier_patches_with_options(ua, ub, va, vb, BezierPatchExtractionOptions::default())
    }
    /// Complete rational rectangle extraction. Cartesian patch/control limits
    /// are checked before constructing any tensor transformation or patch.
    pub fn bezier_patches_with_options(
        &self,
        ua: &R,
        ub: &R,
        va: &R,
        vb: &R,
        options: BezierPatchExtractionOptions,
    ) -> Result<Vec<ExactBezierSurface3>> {
        let [du, dv] = self.degrees();
        let per_patch = (du + 1) * (dv + 1);
        let limit = options.max_patches.min(options.max_controls / per_patch);
        let us = self.axes[0].spans_in(ua, ub, limit)?;
        let vs = self.axes[1].spans_in(va, vb, limit / us.len())?;
        let count = us
            .len()
            .checked_mul(vs.len())
            .ok_or(Error::ComputationLimit("Bezier patch count"))?;
        let v_maps: Vec<_> = vs
            .iter()
            .map(|v| {
                let length = &v.end - &v.start;
                spline::BezierSpanTransform::from_knots(
                    dv,
                    &self.axes[1].flat,
                    v.index,
                    &((&v.lower - &v.start) / &length),
                    &((&v.upper - &v.start) / length),
                )
            })
            .collect();
        let mut result = Vec::with_capacity(count);
        for u in &us {
            let length = &u.end - &u.start;
            let u_map = spline::BezierSpanTransform::from_knots(
                du,
                &self.axes[0].flat,
                u.index,
                &((&u.lower - &u.start) / &length),
                &((&u.upper - &u.start) / length),
            );
            for (v, v_map) in vs.iter().zip(&v_maps) {
                let columns: Vec<_> = (v.index - dv..=v.index)
                    .map(|j| {
                        let poles: Vec<_> = (u.index - du..=u.index)
                            .map(|i| {
                                self.controls[self.axes[0].pole_index(i)
                                    * self.axes[1].pole_count()
                                    + self.axes[1].pole_index(j)]
                                .clone()
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
                result.push(ExactBezierSurface3::from_homogeneous(
                    [du, dv],
                    controls,
                    [
                        [u.lower.clone(), u.upper.clone()],
                        [v.lower.clone(), v.upper.clone()],
                    ],
                ));
            }
        }
        Ok(result)
    }
}
