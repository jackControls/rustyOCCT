//! Exact homogeneous spline substitution into sphere/cylinder equations.
//! OCCT references: IntCurveSurface_QuadricCurveExactInterUtils::PerformIntersection
//! and IntSurf_Quadric::Distance. Positive weights let us clear denominators
//! without introducing roots. Equations have degree at most twice the curve's
//! degree (50); the shared exact isolator retains every root and zero span.
use super::{Cylinder3, Sphere3, SplineSurfaceIntersection, SplineSurfaceOptions};
use crate::polynomial::real::{rat, IntPolynomial};
use crate::{BSplineCurve3, Result};
use num_rational::BigRational as R;

/// Complete closed fundamental domain against the sphere surface. Being inside
/// the ball is not an intersection. Identically contained spans are overlaps.
pub fn spline_sphere(curve: &BSplineCurve3, sphere: &Sphere3) -> Result<SplineSurfaceIntersection> {
    spline_sphere_with_options(curve, sphere, SplineSurfaceOptions::default())
}
pub fn spline_sphere_with_options(
    curve: &BSplineCurve3,
    sphere: &Sphere3,
    options: SplineSurfaceOptions,
) -> Result<SplineSurfaceIntersection> {
    let (first, last) = curve.domain();
    spline_sphere_in_with_options(curve, sphere, first, last, options)
}
/// Closed positive-length parameter interval. Periodic parameters are preserved
/// across any number of turns, subject to the same work limits as plane queries.
pub fn spline_sphere_in(
    curve: &BSplineCurve3,
    sphere: &Sphere3,
    first: f64,
    last: f64,
) -> Result<SplineSurfaceIntersection> {
    spline_sphere_in_with_options(curve, sphere, first, last, SplineSurfaceOptions::default())
}
pub fn spline_sphere_in_with_options(
    curve: &BSplineCurve3,
    sphere: &Sphere3,
    first: f64,
    last: f64,
    options: SplineSurfaceOptions,
) -> Result<SplineSurfaceIntersection> {
    intersect(
        curve,
        sphere.center().to_array().map(rat),
        rat(sphere.radius()),
        None,
        first,
        last,
        options,
    )
}
/// Complete closed fundamental domain against an infinite cylinder surface;
/// there are no end caps. The represented axis is never rounded to a unit vector.
pub fn spline_cylinder(
    curve: &BSplineCurve3,
    cylinder: &Cylinder3,
) -> Result<SplineSurfaceIntersection> {
    spline_cylinder_with_options(curve, cylinder, SplineSurfaceOptions::default())
}
pub fn spline_cylinder_with_options(
    curve: &BSplineCurve3,
    cylinder: &Cylinder3,
    options: SplineSurfaceOptions,
) -> Result<SplineSurfaceIntersection> {
    let (first, last) = curve.domain();
    spline_cylinder_in_with_options(curve, cylinder, first, last, options)
}
/// Closed interval against the infinite cylinder. Query endpoints are included;
/// no snapping, extrapolation or reversal of curve sense is applied.
pub fn spline_cylinder_in(
    curve: &BSplineCurve3,
    cylinder: &Cylinder3,
    first: f64,
    last: f64,
) -> Result<SplineSurfaceIntersection> {
    spline_cylinder_in_with_options(
        curve,
        cylinder,
        first,
        last,
        SplineSurfaceOptions::default(),
    )
}
pub fn spline_cylinder_in_with_options(
    curve: &BSplineCurve3,
    cylinder: &Cylinder3,
    first: f64,
    last: f64,
    options: SplineSurfaceOptions,
) -> Result<SplineSurfaceIntersection> {
    intersect(
        curve,
        cylinder.origin().to_array().map(rat),
        rat(cylinder.radius()),
        Some(cylinder.axis().to_array().map(rat)),
        first,
        last,
        options,
    )
}

fn intersect(
    curve: &BSplineCurve3,
    center: [R; 3],
    radius: R,
    axis: Option<[R; 3]>,
    first: f64,
    last: f64,
    options: SplineSurfaceOptions,
) -> Result<SplineSurfaceIntersection> {
    super::spline_surface::intersect(
        curve,
        first,
        last,
        options,
        128,
        polynomial(center, radius, axis),
    )
}

pub(super) fn polynomial(
    center: [R; 3],
    radius: R,
    axis: Option<[R; 3]>,
) -> impl Fn(&[Vec<R>; 4]) -> IntPolynomial {
    let radius2 = &radius * &radius;
    let norm: R = axis
        .as_ref()
        .map_or_else(|| rat(1.), |v| v.iter().map(|x| x * x).sum());
    move |h| {
        let n = h[3].len();
        let delta: [Vec<R>; 3] = std::array::from_fn(|c| {
            h[c].iter()
                .zip(&h[3])
                .map(|(p, w)| p - &center[c] * w)
                .collect()
        });
        let projection: Vec<R> = (0..n)
            .map(|i| {
                axis.as_ref()
                    .map_or_else(|| rat(0.), |v| (0..3).map(|c| &v[c] * &delta[c][i]).sum())
            })
            .collect();
        let mut f = vec![rat(0.); 2 * n - 1];
        for i in 0..n {
            for j in 0..n {
                let dot: R = (0..3).map(|c| &delta[c][i] * &delta[c][j]).sum();
                f[i + j] += &norm * (dot - &radius2 * &h[3][i] * &h[3][j])
                    - &projection[i] * &projection[j];
            }
        }
        IntPolynomial::from_rationals(&f)
    }
}
