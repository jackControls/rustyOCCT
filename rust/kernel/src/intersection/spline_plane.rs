//! Certified rational spline/plane intersections, including zero spans.
//! OCCT references: GeomAPI_IntCS, quadric/curve scalar substitution in
//! IntCurveSurface, and math_FunctionAllRoots' point/interval result model.
//! Exact homogeneous span polynomials and Sturm isolation replace sampling and
//! tolerance-based zero intervals. All calculations refer to represented inputs.
use super::Plane3;
use crate::polynomial::{
    real::{self, Budget, IntPolynomial},
    AlgebraicRoot, RootIsolationOptions,
};
use crate::{
    interval, math::finite, BSplineCurve3, Bounds3, Error, Point3, Result, ScalarInterval,
};
use num_rational::BigRational as R;
use std::{cmp::Ordering, sync::Arc};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplinePlaneContact {
    /// The curve changes plane side, including at a nonsmooth knot.
    Crossing,
    /// An isolated contact with the same plane side on either side.
    Tangent,
    /// The hit is at an endpoint of the queried closed parameter domain.
    Boundary,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SplinePlaneOverlap {
    exact: [R; 2],
    bounds: [ScalarInterval; 2],
}
impl SplinePlaneOverlap {
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
pub struct SplinePlaneOptions {
    pub root_isolation: RootIsolationOptions,
    pub max_spans: usize,
}
impl Default for SplinePlaneOptions {
    fn default() -> Self {
        Self {
            root_isolation: RootIsolationOptions::default(),
            max_spans: crate::spline::MAX_POLES,
        }
    }
}

#[derive(Debug)]
struct Span {
    start: R,
    end: R,
    homogeneous: [Vec<R>; 4],
    control_bounds: [(f64, f64); 3],
}
impl Span {
    fn compare_coordinate(&self, root: &AlgebraicRoot, c: usize, x: f64) -> Ordering {
        let (lo, hi) = self.control_bounds[c];
        if x < lo {
            return Ordering::Greater;
        }
        if x > hi {
            return Ordering::Less;
        }
        if lo == hi {
            return Ordering::Equal;
        }
        let x = real::rat(x);
        let coefficients: Vec<_> = self.homogeneous[c]
            .iter()
            .zip(&self.homogeneous[3])
            .map(|(a, w)| a - &x * w)
            .collect();
        root.sign_polynomial(&IntPolynomial::from_rationals(&coefficients))
    }
    fn compare_parameter(&self, root: &AlgebraicRoot, x: f64) -> Ordering {
        root.compare_rational(&((real::rat(x) - &self.start) / (&self.end - &self.start)))
    }
}

/// Retains the exact local algebraic parameter and homogeneous coordinates.
/// Distinct points in the result have distinct parameters even when their
/// binary64 enclosures coincide or the curve visits the same position twice.
#[derive(Debug, Clone)]
pub struct SplinePlanePoint {
    span: Arc<Span>,
    root: AlgebraicRoot,
    parameter: ScalarInterval,
    coordinates: [ScalarInterval; 3],
    multiplicities: [Option<usize>; 2],
    contact: SplinePlaneContact,
}
impl SplinePlanePoint {
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
    pub fn contact(&self) -> SplinePlaneContact {
        self.contact
    }
    /// Contact orders in the left and right nonzero span polynomials. None
    /// denotes the outside of the queried parameter domain. At a knot, orders
    /// can differ; there is no invented single multiplicity across a corner.
    pub fn multiplicities(&self) -> [Option<usize>; 2] {
        self.multiplicities
    }
    pub fn compare_parameter(&self, value: f64) -> Result<Ordering> {
        finite(value, "intersection parameter comparison")?;
        Ok(self.span.compare_parameter(&self.root, value))
    }
    pub fn compare_coordinate(&self, component: usize, value: f64) -> Result<Ordering> {
        finite(value, "intersection coordinate comparison")?;
        if component >= 3 {
            return Err(Error::OutOfDomain("coordinate component"));
        }
        Ok(self.span.compare_coordinate(&self.root, component, value))
    }
}

#[derive(Debug, Clone)]
pub struct SplinePlaneIntersection {
    points: Vec<SplinePlanePoint>,
    overlaps: Vec<SplinePlaneOverlap>,
}
impl SplinePlaneIntersection {
    /// Isolated points in exact parameter order. Overlap endpoints are excluded.
    pub fn points(&self) -> &[SplinePlanePoint] {
        &self.points
    }
    pub fn overlaps(&self) -> &[SplinePlaneOverlap] {
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

/// Intersect the complete closed fundamental domain with an infinite plane.
/// Periodic start and end are distinct parameter-boundary events even when
/// they represent the same position. Adjacent knot hits are merged exactly;
/// no proximity threshold merges distinct roots. Contained spans become maximal
/// closed overlap intervals. The function never returns a partial success.
pub fn spline_plane(curve: &BSplineCurve3, plane: &Plane3) -> Result<SplinePlaneIntersection> {
    spline_plane_with_options(curve, plane, RootIsolationOptions::default())
}
pub fn spline_plane_with_options(
    curve: &BSplineCurve3,
    plane: &Plane3,
    options: RootIsolationOptions,
) -> Result<SplinePlaneIntersection> {
    let (first, last) = curve.domain();
    spline_plane_in_with_options(
        curve,
        plane,
        first,
        last,
        SplinePlaneOptions {
            root_isolation: options,
            ..SplinePlaneOptions::default()
        },
    )
}

/// Intersect a finite closed parameter interval of positive length. Nonperiodic
/// bounds must lie in the curve domain. Periodic bounds remain in the caller's
/// parameter units, may cross any seam and may span multiple turns. No snapping,
/// automatic period adjustment, extrapolation or reversed-parameter sense is
/// applied. Endpoints are included and reported as boundary contacts.
pub fn spline_plane_in(
    curve: &BSplineCurve3,
    plane: &Plane3,
    first: f64,
    last: f64,
) -> Result<SplinePlaneIntersection> {
    spline_plane_in_with_options(curve, plane, first, last, SplinePlaneOptions::default())
}
pub fn spline_plane_in_with_options(
    curve: &BSplineCurve3,
    plane: &Plane3,
    first: f64,
    last: f64,
    options: SplinePlaneOptions,
) -> Result<SplinePlaneIntersection> {
    let spans = curve
        .knot_vector()
        .spans_in(first, last, options.max_spans)?;
    let mut budget = Budget::new(options.root_isolation);
    let mut candidates: Vec<Candidate> = Vec::new();
    let mut overlaps: Vec<[R; 2]> = Vec::new();
    let normal = plane.normal.clone().map(R::from_integer);
    let anchor = plane.vertices[0].to_array().map(real::rat);
    for at in spans {
        let (start, end, index) = (at.start, at.end, at.index);
        let length = &end - &start;
        let lower = (&at.lower - &start) / &length;
        let upper = (&at.upper - &start) / &length;
        let homogeneous = curve.span_polynomial(index);
        let coefficients: Vec<R> = (0..=curve.degree())
            .map(|i| {
                (0..3)
                    .map(|c| &normal[c] * (&homogeneous[c][i] - &anchor[c] * &homogeneous[3][i]))
                    .sum()
            })
            .collect();
        let polynomial = IntPolynomial::from_rationals(&coefficients);
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
        let control_bounds = std::array::from_fn(|c| {
            let mut values = (index - curve.degree()..=index)
                .map(|i| curve.poles()[curve.knot_vector().pole_index(i)].to_array()[c]);
            let first = values.next().unwrap();
            values.fold((first, first), |(lo, hi), x| (lo.min(x), hi.max(x)))
        });
        let span = Arc::new(Span {
            start,
            end,
            homogeneous,
            control_bounds,
        });
        for root in roots {
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
            [Some(a), Some(b)] if a != b => SplinePlaneContact::Crossing,
            [Some(_), Some(_)] => SplinePlaneContact::Tangent,
            _ => SplinePlaneContact::Boundary,
        };
        let coordinate = |c| {
            interval::enclose(
                |x| item.span.compare_coordinate(&item.root, c, x),
                "spline intersection coordinate",
            )
        };
        points.push(SplinePlanePoint {
            parameter: interval::enclose(
                |x| item.span.compare_parameter(&item.root, x),
                "spline intersection parameter",
            )?,
            coordinates: [coordinate(0)?, coordinate(1)?, coordinate(2)?],
            span: item.span,
            root: item.root,
            multiplicities: item.multiplicities,
            contact,
        });
    }
    let overlaps = overlaps
        .into_iter()
        .map(|exact| {
            let bound = |i: usize| {
                interval::enclose(|x| exact[i].cmp(&real::rat(x)), "spline overlap endpoint")
            };
            Ok(SplinePlaneOverlap {
                bounds: [bound(0)?, bound(1)?],
                exact,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(SplinePlaneIntersection { points, overlaps })
}

#[cfg(test)]
mod tests {
    use super::*;
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
            let local = (&expected - &point.span.start) / (&point.span.end - &point.span.start);
            assert_eq!(point.root.compare_rational(&local), Ordering::Equal);
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
