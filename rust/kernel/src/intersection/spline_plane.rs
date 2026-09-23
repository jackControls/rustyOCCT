//! Certified rational spline/plane intersections. See SOURCE_MAP.md.
use super::Plane3;
use crate::polynomial::{
    real::{self, IntPolynomial},
    RootIsolationOptions,
};
use crate::{BSplineCurve3, Result};
use num_rational::BigRational as R;
// Preserve the original plane API while sharing exact result identities.
pub use super::spline_surface::{
    SplineSurfaceContact as SplinePlaneContact,
    SplineSurfaceIntersection as SplinePlaneIntersection,
    SplineSurfaceOptions as SplinePlaneOptions, SplineSurfaceOverlap as SplinePlaneOverlap,
    SplineSurfacePoint as SplinePlanePoint,
};

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
    super::spline_surface::intersect(curve, first, last, options, 32, polynomial(plane))
}

pub(super) fn polynomial(plane: &Plane3) -> impl Fn(&[Vec<R>; 4]) -> IntPolynomial {
    let normal = plane.normal.clone().map(R::from_integer);
    let anchor = plane.vertices[0].to_array().map(real::rat);
    move |homogeneous| {
        let coefficients: Vec<R> = (0..homogeneous[3].len())
            .map(|i| {
                (0..3)
                    .map(|c| &normal[c] * (&homogeneous[c][i] - &anchor[c] * &homogeneous[3][i]))
                    .sum()
            })
            .collect();
        IntPolynomial::from_rationals(&coefficients)
    }
}
