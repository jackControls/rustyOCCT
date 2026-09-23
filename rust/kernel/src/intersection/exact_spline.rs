//! Exact edited B-spline intersections with retained algebraic results.
//! The same implicit equations, root isolation and contact merging power the
//! binary64 curve entry points. No knot, pole, trim or result is rounded here.
use super::{Cylinder3, ExactSplineSurfaceIntersection, Plane3, Sphere3, SplineSurfaceOptions};
use crate::polynomial::real::rat;
use crate::{ExactBSplineCurve3, Result};
use num_rational::BigRational as R;

/// Intersect the complete closed fundamental domain with the infinite plane.
/// Exact contacts and overlaps remain available outside binary64 range.
pub fn exact_spline_plane(
    curve: &ExactBSplineCurve3,
    plane: &Plane3,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_plane_with_options(curve, plane, SplineSurfaceOptions::default())
}
/// Complete fundamental domain with explicit span/root work limits.
pub fn exact_spline_plane_with_options(
    curve: &ExactBSplineCurve3,
    plane: &Plane3,
    options: SplineSurfaceOptions,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_plane_in_with_options(
        curve,
        plane,
        &curve.domain()[0],
        &curve.domain()[1],
        options,
    )
}
/// Closed positive-length rational interval, without snapping or extrapolation.
/// Periodic intervals preserve the caller's parameters across seams and turns.
pub fn exact_spline_plane_in(
    curve: &ExactBSplineCurve3,
    plane: &Plane3,
    first: &R,
    last: &R,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_plane_in_with_options(curve, plane, first, last, SplineSurfaceOptions::default())
}
/// Rational interval with explicit span/root limits. Exhaustion returns an
/// error, never a partial result. Enclosures can be requested on the result.
pub fn exact_spline_plane_in_with_options(
    curve: &ExactBSplineCurve3,
    plane: &Plane3,
    first: &R,
    last: &R,
    options: SplineSurfaceOptions,
) -> Result<ExactSplineSurfaceIntersection> {
    super::spline_surface::intersect_exact(
        curve,
        first,
        last,
        options,
        0,
        super::spline_plane::polynomial(plane),
    )
}

/// Intersect the complete closed fundamental domain with the sphere surface.
/// Exact contacts and overlaps remain available outside binary64 range.
pub fn exact_spline_sphere(
    curve: &ExactBSplineCurve3,
    sphere: &Sphere3,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_sphere_with_options(curve, sphere, SplineSurfaceOptions::default())
}
/// Complete fundamental domain with explicit span/root work limits.
pub fn exact_spline_sphere_with_options(
    curve: &ExactBSplineCurve3,
    sphere: &Sphere3,
    options: SplineSurfaceOptions,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_sphere_in_with_options(
        curve,
        sphere,
        &curve.domain()[0],
        &curve.domain()[1],
        options,
    )
}
/// Closed positive-length rational interval, without snapping or extrapolation.
/// Periodic intervals preserve the caller's parameters across seams and turns.
pub fn exact_spline_sphere_in(
    curve: &ExactBSplineCurve3,
    sphere: &Sphere3,
    first: &R,
    last: &R,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_sphere_in_with_options(curve, sphere, first, last, SplineSurfaceOptions::default())
}
/// Rational interval with explicit span/root limits. Exhaustion returns an
/// error, never a partial result. Enclosures can be requested on the result.
pub fn exact_spline_sphere_in_with_options(
    curve: &ExactBSplineCurve3,
    sphere: &Sphere3,
    first: &R,
    last: &R,
    options: SplineSurfaceOptions,
) -> Result<ExactSplineSurfaceIntersection> {
    super::spline_surface::intersect_exact(
        curve,
        first,
        last,
        options,
        128,
        super::spline_quadric::polynomial(
            sphere.center().to_array().map(rat),
            rat(sphere.radius()),
            None,
        ),
    )
}

/// Intersect the complete closed fundamental domain with the infinite cylinder surface, without end caps.
/// Exact contacts and overlaps remain available outside binary64 range.
pub fn exact_spline_cylinder(
    curve: &ExactBSplineCurve3,
    cylinder: &Cylinder3,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_cylinder_with_options(curve, cylinder, SplineSurfaceOptions::default())
}
/// Complete fundamental domain with explicit span/root work limits.
pub fn exact_spline_cylinder_with_options(
    curve: &ExactBSplineCurve3,
    cylinder: &Cylinder3,
    options: SplineSurfaceOptions,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_cylinder_in_with_options(
        curve,
        cylinder,
        &curve.domain()[0],
        &curve.domain()[1],
        options,
    )
}
/// Closed positive-length rational interval, without snapping or extrapolation.
/// Periodic intervals preserve the caller's parameters across seams and turns.
pub fn exact_spline_cylinder_in(
    curve: &ExactBSplineCurve3,
    cylinder: &Cylinder3,
    first: &R,
    last: &R,
) -> Result<ExactSplineSurfaceIntersection> {
    exact_spline_cylinder_in_with_options(
        curve,
        cylinder,
        first,
        last,
        SplineSurfaceOptions::default(),
    )
}
/// Rational interval with explicit span/root limits. Exhaustion returns an
/// error, never a partial result. Enclosures can be requested on the result.
pub fn exact_spline_cylinder_in_with_options(
    curve: &ExactBSplineCurve3,
    cylinder: &Cylinder3,
    first: &R,
    last: &R,
    options: SplineSurfaceOptions,
) -> Result<ExactSplineSurfaceIntersection> {
    super::spline_surface::intersect_exact(
        curve,
        first,
        last,
        options,
        128,
        super::spline_quadric::polynomial(
            cylinder.origin().to_array().map(rat),
            rat(cylinder.radius()),
            Some(cylinder.axis().to_array().map(rat)),
        ),
    )
}
