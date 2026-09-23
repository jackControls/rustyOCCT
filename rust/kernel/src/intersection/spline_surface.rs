//! Shared certified rational spline/implicit-surface intersection engine.
//! OCCT references: GeomAPI_IntCS, quadric/curve scalar substitution in
//! IntCurveSurface, and math_FunctionAllRoots' point/interval result model.
//! Exact homogeneous span polynomials and Sturm isolation replace sampling and
//! tolerance-based zero intervals. All calculations refer to represented inputs.
use crate::polynomial::{
    real::{self, Budget, IntPolynomial},
    AlgebraicRoot, RootIsolationOptions,
};
use crate::{
    interval, math::finite, spline, BSplineCurve3, Bounds3, Error, ExactBSplineCurve3, Point3,
    Result, ScalarInterval,
};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::{cmp::Ordering, sync::Arc};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplineSurfaceContact {
    /// The curve changes implicit surface sign, including at a nonsmooth knot.
    Crossing,
    /// An isolated contact with the same implicit surface sign on either side.
    Tangent,
    /// The hit is at an endpoint of the queried closed parameter domain.
    Boundary,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SplineSurfaceOverlap {
    exact: [R; 2],
    bounds: [ScalarInterval; 2],
}
impl SplineSurfaceOverlap {
    /// Representatives of the maximal closed overlap endpoints. For shifted
    /// periodic knots these may be rounded; retain `parameter_bounds` or compare
    /// against the exact endpoint. Distinct endpoints may share an enclosure.
    pub fn parameters(&self) -> (f64, f64) {
        (
            self.bounds[0].representative(),
            self.bounds[1].representative(),
        )
    }
    pub fn parameter_bounds(&self) -> [ScalarInterval; 2] {
        self.bounds
    }
    pub fn compare_parameter(&self, endpoint: usize, value: f64) -> Result<Ordering> {
        finite(value, "overlap endpoint comparison")?;
        let exact = self
            .exact
            .get(endpoint)
            .ok_or(Error::OutOfDomain("overlap endpoint"))?;
        Ok(exact.cmp(&real::rat(value)))
    }
}

/// Independent limits on visited spans and polynomial subdivisions. Neither
/// limit is a floating geometric tolerance or a hard CPU/allocation deadline.
#[derive(Debug, Clone, Copy)]
pub struct SplineSurfaceOptions {
    pub root_isolation: RootIsolationOptions,
    pub max_spans: usize,
}
impl Default for SplineSurfaceOptions {
    fn default() -> Self {
        Self {
            root_isolation: RootIsolationOptions::default(),
            max_spans: crate::spline::MAX_POLES,
        }
    }
}

/// Exact maximal closed interval contained in the surface. Endpoints can be
/// outside binary64 range or arbitrarily close together.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactSplineSurfaceOverlap {
    parameters: [R; 2],
}
impl ExactSplineSurfaceOverlap {
    pub fn parameters(&self) -> &[R; 2] {
        &self.parameters
    }
    pub fn compare_parameter(&self, endpoint: usize, value: &R) -> Result<Ordering> {
        let value = spline::normalize(value)?;
        Ok(self
            .parameters
            .get(endpoint)
            .ok_or(Error::OutOfDomain("overlap endpoint"))?
            .cmp(&value))
    }
    pub fn parameter_bounds(&self) -> Result<[ScalarInterval; 2]> {
        let bound = |i: usize| {
            interval::enclose(
                |x| self.parameters[i].cmp(&real::rat(x)),
                "spline overlap endpoint",
            )
        };
        Ok([bound(0)?, bound(1)?])
    }
    pub fn enclosed(&self) -> Result<SplineSurfaceOverlap> {
        Ok(SplineSurfaceOverlap {
            exact: self.parameters.clone(),
            bounds: self.parameter_bounds()?,
        })
    }
}

/// Exact algebraic contact in the original rational parameter units. No
/// floating enclosure is required to retain or compare this result.
#[derive(Debug, Clone)]
pub struct ExactSplineSurfacePoint {
    span: Arc<Span>,
    root: AlgebraicRoot,
    multiplicities: [Option<usize>; 2],
    contact: SplineSurfaceContact,
}
impl ExactSplineSurfacePoint {
    pub fn contact(&self) -> SplineSurfaceContact {
        self.contact
    }
    pub fn multiplicities(&self) -> [Option<usize>; 2] {
        self.multiplicities
    }
    pub fn compare_parameter(&self, value: &R) -> Result<Ordering> {
        Ok(self
            .span
            .compare_rational_parameter(&self.root, &spline::normalize(value)?))
    }
    pub fn compare_coordinate(&self, component: usize, value: &R) -> Result<Ordering> {
        if component >= 3 {
            return Err(Error::OutOfDomain("coordinate component"));
        }
        Ok(self
            .span
            .compare_rational_coordinate(&self.root, component, &spline::normalize(value)?))
    }
    pub fn parameter_bounds(&self) -> Result<ScalarInterval> {
        interval::enclose(
            |x| self.span.compare_parameter(&self.root, x),
            "spline intersection parameter",
        )
    }
    pub fn coordinate_bound(&self, component: usize) -> Result<ScalarInterval> {
        if component >= 3 {
            return Err(Error::OutOfDomain("coordinate component"));
        }
        interval::enclose(
            |x| self.span.compare_coordinate(&self.root, component, x),
            "spline intersection coordinate",
        )
    }
    pub fn coordinate_bounds(&self) -> Result<[ScalarInterval; 3]> {
        Ok([
            self.coordinate_bound(0)?,
            self.coordinate_bound(1)?,
            self.coordinate_bound(2)?,
        ])
    }
    /// All finite enclosures, or a typed error. The exact contact is unchanged.
    pub fn enclosed(&self) -> Result<SplineSurfacePoint> {
        self.clone().into_enclosed()
    }
    fn into_enclosed(self) -> Result<SplineSurfacePoint> {
        Ok(SplineSurfacePoint {
            parameter: self.parameter_bounds()?,
            coordinates: self.coordinate_bounds()?,
            exact: self,
        })
    }
}

/// Complete exact contacts and overlaps. Positive weights permit denominator
/// clearing without adding roots. Failure never returns a partial result.
#[derive(Debug, Clone)]
pub struct ExactSplineSurfaceIntersection {
    points: Vec<ExactSplineSurfacePoint>,
    overlaps: Vec<ExactSplineSurfaceOverlap>,
}
impl ExactSplineSurfaceIntersection {
    pub fn points(&self) -> &[ExactSplineSurfacePoint] {
        &self.points
    }
    pub fn overlaps(&self) -> &[ExactSplineSurfaceOverlap] {
        &self.overlaps
    }
    pub fn is_disjoint(&self) -> bool {
        self.points.is_empty() && self.overlaps.is_empty()
    }
    /// Enclose every parameter and coordinate atomically. Exact results remain
    /// usable if any value is outside finite binary64 range.
    pub fn enclosed(&self) -> Result<SplineSurfaceIntersection> {
        self.clone().into_enclosed()
    }
    fn into_enclosed(self) -> Result<SplineSurfaceIntersection> {
        Ok(SplineSurfaceIntersection {
            points: self
                .points
                .into_iter()
                .map(ExactSplineSurfacePoint::into_enclosed)
                .collect::<Result<_>>()?,
            overlaps: self
                .overlaps
                .iter()
                .map(ExactSplineSurfaceOverlap::enclosed)
                .collect::<Result<_>>()?,
        })
    }
}

#[derive(Debug)]
enum ControlBounds {
    Binary64([(f64, f64); 3]),
    Rational(Box<[(R, R); 3]>),
}
impl ControlBounds {
    fn compare_rational(&self, component: usize, value: &R) -> Option<Ordering> {
        let compare = |lo: &R, hi: &R| {
            if value < lo {
                Some(Ordering::Greater)
            } else if value > hi {
                Some(Ordering::Less)
            } else if lo == hi {
                Some(Ordering::Equal)
            } else {
                None
            }
        };
        match self {
            Self::Binary64(bounds) => {
                let (lo, hi) = bounds[component];
                compare(&real::rat(lo), &real::rat(hi))
            }
            Self::Rational(bounds) => compare(&bounds[component].0, &bounds[component].1),
        }
    }
    fn compare_float(&self, component: usize, value: f64) -> Option<Ordering> {
        if let Self::Binary64(bounds) = self {
            let (lo, hi) = bounds[component];
            if value < lo {
                Some(Ordering::Greater)
            } else if value > hi {
                Some(Ordering::Less)
            } else if lo == hi {
                Some(Ordering::Equal)
            } else {
                None
            }
        } else {
            self.compare_rational(component, &real::rat(value))
        }
    }
}

#[derive(Debug)]
struct Span {
    start: R,
    end: R,
    // All four coordinates share one positive scale. Clearing denominators
    // once preserves H_c/H_w and avoids rational reductions at every query.
    homogeneous: [Vec<BigInt>; 4],
    control_bounds: ControlBounds,
}
impl Span {
    fn compare_coordinate(&self, root: &AlgebraicRoot, c: usize, x: f64) -> Ordering {
        self.control_bounds
            .compare_float(c, x)
            .unwrap_or_else(|| self.coordinate_sign(root, c, &real::rat(x)))
    }
    fn compare_rational_coordinate(&self, root: &AlgebraicRoot, c: usize, x: &R) -> Ordering {
        self.control_bounds
            .compare_rational(c, x)
            .unwrap_or_else(|| self.coordinate_sign(root, c, x))
    }
    fn coordinate_sign(&self, root: &AlgebraicRoot, c: usize, x: &R) -> Ordering {
        let coefficients: Vec<_> = self.homogeneous[c]
            .iter()
            .zip(&self.homogeneous[3])
            .map(|(a, w)| a * x.denom() - x.numer() * w)
            .collect();
        root.sign_polynomial(&IntPolynomial::new(coefficients))
    }
    fn compare_parameter(&self, root: &AlgebraicRoot, x: f64) -> Ordering {
        self.compare_rational_parameter(root, &real::rat(x))
    }
    fn compare_rational_parameter(&self, root: &AlgebraicRoot, x: &R) -> Ordering {
        root.compare_rational(&((x - &self.start) / (&self.end - &self.start)))
    }
}

/// Retains the exact local algebraic parameter and homogeneous coordinates.
/// Distinct points in the result have distinct parameters even when their
/// binary64 enclosures coincide or the curve visits the same position twice.
#[derive(Debug, Clone)]
pub struct SplineSurfacePoint {
    exact: ExactSplineSurfacePoint,
    parameter: ScalarInterval,
    coordinates: [ScalarInterval; 3],
}
impl SplineSurfacePoint {
    pub fn exact(&self) -> &ExactSplineSurfacePoint {
        &self.exact
    }
    pub fn parameter(&self) -> ScalarInterval {
        self.parameter
    }
    pub fn position(&self) -> Point3 {
        let [x, y, z] = self.coordinates.map(ScalarInterval::representative);
        Point3::new(x, y, z)
    }
    pub fn coordinate_bounds(&self) -> [ScalarInterval; 3] {
        self.coordinates
    }
    pub fn bounds(&self) -> Bounds3 {
        let [x, y, z] = self.coordinates.map(ScalarInterval::lower);
        let [a, b, c] = self.coordinates.map(ScalarInterval::upper);
        Bounds3 {
            min: Point3::new(x, y, z),
            max: Point3::new(a, b, c),
        }
    }
    pub fn contact(&self) -> SplineSurfaceContact {
        self.exact.contact
    }
    /// Contact orders in the left and right nonzero span polynomials. None
    /// denotes the outside of the queried parameter domain. At a knot, orders
    /// can differ; there is no invented single multiplicity across a corner.
    pub fn multiplicities(&self) -> [Option<usize>; 2] {
        self.exact.multiplicities
    }
    pub fn compare_parameter(&self, value: f64) -> Result<Ordering> {
        finite(value, "intersection parameter comparison")?;
        Ok(self.exact.span.compare_parameter(&self.exact.root, value))
    }
    pub fn compare_coordinate(&self, component: usize, value: f64) -> Result<Ordering> {
        finite(value, "intersection coordinate comparison")?;
        if component >= 3 {
            return Err(Error::OutOfDomain("coordinate component"));
        }
        Ok(self
            .exact
            .span
            .compare_coordinate(&self.exact.root, component, value))
    }
}

#[derive(Debug, Clone)]
pub struct SplineSurfaceIntersection {
    points: Vec<SplineSurfacePoint>,
    overlaps: Vec<SplineSurfaceOverlap>,
}
impl SplineSurfaceIntersection {
    /// Isolated points in exact parameter order. Overlap endpoints are excluded.
    pub fn points(&self) -> &[SplineSurfacePoint] {
        &self.points
    }
    pub fn overlaps(&self) -> &[SplineSurfaceOverlap] {
        &self.overlaps
    }
    pub fn is_disjoint(&self) -> bool {
        self.points.is_empty() && self.overlaps.is_empty()
    }
}

struct Candidate {
    span: Arc<Span>,
    root: AlgebraicRoot,
    boundary: Option<R>,
    multiplicities: [Option<usize>; 2],
    signs: [Option<Ordering>; 2],
}

trait CurveInput {
    fn polynomial(&self, span: usize) -> [Vec<R>; 4];
    fn control_bounds(&self, span: usize) -> ControlBounds;
}
impl CurveInput for BSplineCurve3 {
    fn polynomial(&self, span: usize) -> [Vec<R>; 4] {
        self.span_polynomial(span)
    }
    fn control_bounds(&self, span: usize) -> ControlBounds {
        ControlBounds::Binary64(std::array::from_fn(|c| {
            let mut values = (span - self.degree()..=span)
                .map(|i| self.poles()[self.knot_vector().pole_index(i)].to_array()[c]);
            let first = values.next().unwrap();
            values.fold((first, first), |(lo, hi), x| (lo.min(x), hi.max(x)))
        }))
    }
}
impl CurveInput for ExactBSplineCurve3 {
    fn polynomial(&self, span: usize) -> [Vec<R>; 4] {
        self.span_polynomial(span)
    }
    fn control_bounds(&self, span: usize) -> ControlBounds {
        ControlBounds::Rational(Box::new(std::array::from_fn(|c| {
            let mut values = (span - self.degree()..=span).map(|i| {
                let p = &self.homogeneous_poles()[self.knot_vector().pole_index(i)];
                &p[c] / &p[3]
            });
            let first = values.next().unwrap();
            values.fold((first.clone(), first), |(lo, hi), x| {
                (lo.min(x.clone()), hi.max(x))
            })
        })))
    }
}

pub(crate) fn intersect(
    curve: &BSplineCurve3,
    first: f64,
    last: f64,
    options: SplineSurfaceOptions,
    refinement_steps: usize,
    polynomial: impl Fn(&[Vec<R>; 4]) -> IntPolynomial,
) -> Result<SplineSurfaceIntersection> {
    let spans = curve
        .knot_vector()
        .spans_in(first, last, options.max_spans)?;
    collect(curve, spans, options, refinement_steps, polynomial)?.into_enclosed()
}

pub(crate) fn intersect_exact(
    curve: &ExactBSplineCurve3,
    first: &R,
    last: &R,
    options: SplineSurfaceOptions,
    refinement_steps: usize,
    polynomial: impl Fn(&[Vec<R>; 4]) -> IntPolynomial,
) -> Result<ExactSplineSurfaceIntersection> {
    let spans = curve
        .knot_vector()
        .spans_in(first, last, options.max_spans)?;
    collect(curve, spans, options, refinement_steps, polynomial)
}

fn collect(
    curve: &impl CurveInput,
    spans: Vec<spline::KnotSpan>,
    options: SplineSurfaceOptions,
    refinement_steps: usize,
    polynomial: impl Fn(&[Vec<R>; 4]) -> IntPolynomial,
) -> Result<ExactSplineSurfaceIntersection> {
    let mut budget = Budget::new(options.root_isolation);
    let mut candidates: Vec<Candidate> = Vec::new();
    let mut overlaps: Vec<[R; 2]> = Vec::new();
    for at in spans {
        let (start, end, index) = (at.start, at.end, at.index);
        let length = &end - &start;
        let lower = (&at.lower - &start) / &length;
        let upper = (&at.upper - &start) / &length;
        let homogeneous = curve.polynomial(index);
        let polynomial = polynomial(&homogeneous);
        if polynomial.is_zero() {
            if let Some(last) = overlaps.last_mut().filter(|o| o[1] == at.lower) {
                last[1] = at.upper;
            } else {
                overlaps.push([at.lower, at.upper]);
            }
            continue;
        }
        let roots = real::isolate(&polynomial, lower.clone(), upper.clone(), &mut budget)?;
        if roots.is_empty() {
            continue;
        }
        let control_bounds = curve.control_bounds(index);
        let denominator = spline::common_denominator(homogeneous.iter().flatten());
        let homogeneous = homogeneous.map(|v| {
            v.into_iter()
                .map(|c| c.numer() * (&denominator / c.denom()))
                .collect()
        });
        let span = Arc::new(Span {
            start,
            end,
            homogeneous,
            control_bounds,
        });
        for mut root in roots {
            // Bounded refinement also recognizes exactly verified rational
            // roots. Quadrics request more steps for their higher-degree signs.
            root.refine_for_signs(refinement_steps);
            let at_start = root.compare_rational(&lower) == Ordering::Equal;
            let at_end = root.compare_rational(&upper) == Ordering::Equal;
            let m = root.multiplicity();
            let boundary = if at_start {
                Some(at.lower.clone())
            } else if at_end {
                Some(at.upper.clone())
            } else {
                None
            };
            let mut derivative = polynomial.clone();
            for _ in 0..m {
                derivative = derivative.derivative();
            }
            let right_sign = root.sign_polynomial(&derivative);
            debug_assert!(right_sign != Ordering::Equal);
            let left_sign = if m % 2 == 0 {
                right_sign
            } else {
                right_sign.reverse()
            };
            let item = Candidate {
                span: span.clone(),
                root,
                boundary: boundary.clone(),
                multiplicities: [(!at_start).then_some(m), (!at_end).then_some(m)],
                signs: [
                    (!at_start).then_some(left_sign),
                    (!at_end).then_some(right_sign),
                ],
            };
            if let Some(last) = candidates
                .last_mut()
                .filter(|last| boundary.is_some() && last.boundary == boundary)
            {
                debug_assert!(last.multiplicities[0].is_some() && item.multiplicities[1].is_some());
                last.multiplicities[1] = item.multiplicities[1];
                last.signs[1] = item.signs[1];
            } else {
                candidates.push(item);
            }
        }
    }
    let mut points = Vec::new();
    for item in candidates {
        if item
            .boundary
            .as_ref()
            .is_some_and(|k| overlaps.iter().any(|o| &o[0] <= k && k <= &o[1]))
        {
            continue;
        }
        let contact = match item.signs {
            [Some(a), Some(b)] if a != b => SplineSurfaceContact::Crossing,
            [Some(_), Some(_)] => SplineSurfaceContact::Tangent,
            _ => SplineSurfaceContact::Boundary,
        };
        points.push(ExactSplineSurfacePoint {
            span: item.span,
            root: item.root,
            multiplicities: item.multiplicities,
            contact,
        });
    }
    let overlaps = overlaps
        .into_iter()
        .map(|parameters| ExactSplineSurfaceOverlap { parameters })
        .collect();
    Ok(ExactSplineSurfaceIntersection { points, overlaps })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersection::{spline_plane_in, Plane3};
    use num_bigint::BigInt;

    #[test]
    fn exact_parameters_remain_distinct_when_enclosures_coincide() {
        let first = 2f64.powi(53);
        let plane = Plane3::through_points(
            Point3::new(0., 0., 0.),
            Point3::new(1., 0., 0.),
            Point3::new(0., 1., 0.),
        )
        .unwrap();
        let curve = BSplineCurve3::new_periodic(
            1,
            vec![Point3::new(0., 0., -1.), Point3::new(1., 0., 1.)],
            None,
            vec![0., 0.25, 0.5],
            vec![1; 3],
        )
        .unwrap();
        let hits = spline_plane_in(&curve, &plane, first, first + 2.).unwrap();
        assert_eq!(hits.points.len(), 8);
        for (i, point) in hits.points.iter().enumerate() {
            // Independently known affine crossings: 1/8, 3/8, ... 15/8.
            let expected = real::rat(first) + R::new(BigInt::from(2 * i + 1), BigInt::from(8));
            assert_eq!(
                point.exact().compare_parameter(&expected).unwrap(),
                Ordering::Equal
            );
        }
        let curve = BSplineCurve3::new_periodic(
            1,
            vec![
                Point3::new(0., 0., 0.),
                Point3::new(1., 1., 0.),
                Point3::new(2., 0., 1.),
            ],
            None,
            vec![0., 0.125, 0.375, 0.5],
            vec![1; 4],
        )
        .unwrap();
        let hits = spline_plane_in(&curve, &plane, first, first + 2.).unwrap();
        assert_eq!(hits.overlaps.len(), 4);
        for (i, overlap) in hits.overlaps.iter().enumerate() {
            let start = real::rat(first) + R::new(BigInt::from(i), BigInt::from(2));
            let end = &start + R::new(BigInt::from(1), BigInt::from(8));
            assert_eq!(overlap.exact, [start, end]);
        }
    }
}
