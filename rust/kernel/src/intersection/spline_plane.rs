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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SplinePlaneOverlap {
    start: f64,
    end: f64,
}
impl SplinePlaneOverlap {
    /// Maximal closed interval entirely on the plane. Bounds are original knots.
    pub fn parameters(self) -> (f64, f64) {
        (self.start, self.end)
    }
}

#[derive(Debug)]
struct Span {
    start: f64,
    end: f64,
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
        root.compare_rational(
            &((real::rat(x) - real::rat(self.start))
                / (real::rat(self.end) - real::rat(self.start))),
        )
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
    knot: Option<f64>,
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
    let mut budget = Budget::new(options);
    let mut candidates: Vec<Candidate> = Vec::new();
    let mut overlaps: Vec<SplinePlaneOverlap> = Vec::new();
    let normal = plane.normal.clone().map(R::from_integer);
    let anchor = plane.vertices[0].to_array().map(real::rat);
    for (start, end, index) in curve.knot_vector().spans() {
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
            if let Some(last) = overlaps.last_mut().filter(|o| o.end == start) {
                last.end = end;
            } else {
                overlaps.push(SplinePlaneOverlap { start, end });
            }
            continue;
        }
        let roots = real::isolate(&polynomial, real::rat(0.), real::rat(1.), &mut budget)?;
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
            let at_start = root.compare_rational(&real::rat(0.)) == Ordering::Equal;
            let at_end = root.compare_rational(&real::rat(1.)) == Ordering::Equal;
            let m = root.multiplicity();
            let knot = if at_start {
                Some(start)
            } else if at_end {
                Some(end)
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
                knot,
                multiplicities: [(!at_start).then_some(m), (!at_end).then_some(m)],
                signs: [
                    (!at_start).then_some(left_sign),
                    (!at_end).then_some(right_sign),
                ],
            };
            if let Some(last) = candidates
                .last_mut()
                .filter(|last| knot.is_some() && last.knot == knot)
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
            .knot
            .is_some_and(|k| overlaps.iter().any(|o| o.start <= k && k <= o.end))
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
    Ok(SplinePlaneIntersection { points, overlaps })
}
