//! Certified tolerance decisions (Contract 5, T4, of IDENTITY_AND_HISTORY.md).
//!
//! Every comparison of a computed length against a linear tolerance goes
//! through here: a distance, a sum or difference of coordinates, or an area
//! against the tolerance times a perimeter. Inputs are binary64 values taken
//! exactly. Polynomial decisions run first in binary64 intervals and, if
//! those cannot decide, in exact rationals, so they always decide. The
//! area/perimeter test involves square roots and `π`: it runs in rational
//! intervals and, in the rare case they cannot decide, answers
//! "degenerate", the conservative result for a validity screen.
pub(crate) mod arcs;
pub(crate) mod splines;

use crate::certified::{pi, Fast, Interval, Real};
use crate::Point2;
use num_rational::BigRational as R;
use std::cmp::Ordering;

fn zero() -> R {
    R::from_integer(0.into())
}

fn q(x: f64) -> R {
    R::from_float(x).expect("tolerance decisions take finite values")
}

fn sum_r(terms: &[f64]) -> R {
    terms.iter().map(|x| q(*x)).fold(zero(), |a, b| a + b)
}

fn sum_t<T: Real>(terms: &[f64]) -> T {
    terms
        .iter()
        .fold(T::exact_f64(0.0), |a, x| a.add(&T::exact_f64(*x)))
}

/// `a <= b` if certain.
fn le<T: Real>(a: &T, b: &T) -> Option<bool> {
    Some(a.cmp(b)? != Ordering::Greater)
}

/// The exact sum of `lhs` is at most the exact sum of `rhs`.
pub(crate) fn sum_le(lhs: &[f64], rhs: &[f64]) -> bool {
    sum_r(lhs) <= sum_r(rhs)
}

/// The exact sum of `lhs` exceeds the exact sum of `rhs`.
pub(crate) fn sum_gt(lhs: &[f64], rhs: &[f64]) -> bool {
    !sum_le(lhs, rhs)
}

fn dist2_t<T: Real>(a: Point2, b: Point2) -> T {
    let dx = T::exact_f64(a.x).sub(&T::exact_f64(b.x));
    let dy = T::exact_f64(a.y).sub(&T::exact_f64(b.y));
    dx.square().add(&dy.square())
}

fn dist2_r(a: Point2, b: Point2) -> R {
    let (dx, dy) = (q(a.x) - q(b.x), q(a.y) - q(b.y));
    &dx * &dx + &dy * &dy
}

/// Squared distance from `p` to the segment `ab`, if the binary64 intervals
/// place the foot of the perpendicular certainly.
fn segment_dist2_t<T: Real>(p: Point2, a: Point2, b: Point2) -> Option<T> {
    let e = |x: f64| T::exact_f64(x);
    let (dx, dy) = (e(b.x).sub(&e(a.x)), e(b.y).sub(&e(a.y)));
    let (wx, wy) = (e(p.x).sub(&e(a.x)), e(p.y).sub(&e(a.y)));
    let len2 = dx.square().add(&dy.square());
    let t = wx.mul(&dx).add(&wy.mul(&dy));
    if t.sign()? != Ordering::Greater {
        return Some(wx.square().add(&wy.square()));
    }
    if t.cmp(&len2)? != Ordering::Less {
        return Some(dist2_t(p, b));
    }
    let cross = wx.mul(&dy).sub(&wy.mul(&dx));
    cross.square().div(&len2)
}

fn segment_dist2_r(p: Point2, a: Point2, b: Point2) -> R {
    let (dx, dy) = (q(b.x) - q(a.x), q(b.y) - q(a.y));
    let (wx, wy) = (q(p.x) - q(a.x), q(p.y) - q(a.y));
    let len2 = &dx * &dx + &dy * &dy;
    let t = &wx * &dx + &wy * &dy;
    if len2 == zero() || t <= zero() {
        return &wx * &wx + &wy * &wy;
    }
    if t >= len2 {
        return dist2_r(p, b);
    }
    let cross = &wx * &dy - &wy * &dx;
    &cross * &cross / len2
}

/// `sqrt(d2) <= Σ threshold`, from both tiers.
fn root_le(fast: Option<Fast>, exact: impl FnOnce() -> R, threshold: &[f64]) -> bool {
    let th = sum_r(threshold);
    if th < zero() {
        return false;
    }
    if let Some(d2) = fast {
        let t = sum_t::<Fast>(threshold);
        if let Some(answer) = le(&d2, &t.square()) {
            return answer;
        }
    }
    exact() <= &th * &th
}

/// `sqrt(d2) >= Σ threshold`, from both tiers.
fn root_ge(fast: Option<Fast>, exact: impl FnOnce() -> R, threshold: &[f64]) -> bool {
    let th = sum_r(threshold);
    if th <= zero() {
        return true;
    }
    if let Some(d2) = fast {
        let t = sum_t::<Fast>(threshold).square();
        if let Some(answer) = le(&t, &d2) {
            return answer;
        }
    }
    exact() >= &th * &th
}

/// `|a - b| <= Σ threshold`.
pub(crate) fn distance_le(a: Point2, b: Point2, threshold: &[f64]) -> bool {
    root_le(Some(dist2_t(a, b)), || dist2_r(a, b), threshold)
}

/// `|a - b| >= Σ threshold`.
pub(crate) fn distance_ge(a: Point2, b: Point2, threshold: &[f64]) -> bool {
    root_ge(Some(dist2_t(a, b)), || dist2_r(a, b), threshold)
}

/// The distance from `p` to the segment `ab` is at most `Σ threshold`.
pub(crate) fn segment_distance_le(p: Point2, a: Point2, b: Point2, threshold: &[f64]) -> bool {
    root_le(
        segment_dist2_t(p, a, b),
        || segment_dist2_r(p, a, b),
        threshold,
    )
}

/// One boundary of a profile, for the area screen.
pub(crate) enum Outline<'a> {
    Polygon(&'a [Point2]),
    Circle(f64),
    /// Points and the segment from each to the next (S5).
    Path(&'a [Point2], &'a [crate::profile::Segment]),
}

/// `π` in a tier.
fn pi_t<T: Real>() -> T {
    T::from_r(&pi().midpoint()).widen(&pi().radius())
}

/// An outline's area and perimeter in a tier; `None` when a quotient is not
/// defined there.
fn area_perimeter<T: Real>(outline: &Outline) -> Option<(T, T)> {
    let abs = |x: T| -> T {
        match x.sign() {
            Some(Ordering::Less) => x.neg(),
            Some(_) => x,
            None => T::exact_f64(0.0).union(&x.neg()).union(&x),
        }
    };
    match outline {
        Outline::Polygon(points) => {
            let n = points.len();
            let mut twice = zero();
            let mut perimeter = T::exact_f64(0.0);
            for k in 0..n {
                let (a, b) = (points[k], points[(k + 1) % n]);
                twice += q(a.x) * q(b.y) - q(b.x) * q(a.y);
                perimeter = perimeter.add(&T::from_r(&dist2_r(a, b)).sqrt());
            }
            let area = twice / R::from_integer(2.into());
            let area = if area < zero() { -area } else { area };
            Some((T::from_r(&area), perimeter))
        }
        Outline::Circle(radius) => {
            let r = T::exact_f64(*radius);
            Some((
                pi_t::<T>().mul(&r.square()),
                pi_t::<T>().mul(&r).mul(&T::exact_f64(2.0)),
            ))
        }
        Outline::Path(points, segments) => {
            let n = points.len();
            let mut twice = T::exact_f64(0.0);
            let mut perimeter = T::exact_f64(0.0);
            for (k, segment) in segments.iter().enumerate() {
                let (a, b) = (points[k], points[(k + 1) % n]);
                twice = twice.add(&T::from_r(&(q(a.x) * q(b.y) - q(b.x) * q(a.y))));
                match segment {
                    crate::profile::Segment::Line => {
                        perimeter = perimeter.add(&T::from_r(&dist2_r(a, b)).sqrt());
                    }
                    // S8b: the spline's exact twice-area in place of its
                    // chord's, its length between the chord and the control
                    // polygon.
                    crate::profile::Segment::Spline(span) => {
                        let parts = splines::spline(span)?;
                        let origin = [R::from_integer(0.into()), R::from_integer(0.into())];
                        let g = splines::green(&parts, &origin);
                        let chord = q(a.x) * q(b.y) - q(b.x) * q(a.y);
                        twice = twice.add(&T::from_r(&(&g[0] * R::from_integer(2.into()) - chord)));
                        perimeter = perimeter.add(&splines::length_bounds::<T>(&parts));
                    }
                    crate::profile::Segment::Arc {
                        center,
                        radius,
                        ccw,
                    } => {
                        // The circular segment r^2 (φ - sin φ) and the arc
                        // length r |φ|, φ the signed sweep.
                        let (ux, uy) = (q(a.x) - q(center.x), q(a.y) - q(center.y));
                        let (vx, vy) = (q(b.x) - q(center.x), q(b.y) - q(center.y));
                        let cross = &ux * &vy - &uy * &vx;
                        let dot = &ux * &vx + &uy * &vy;
                        let two_pi = pi_t::<T>().mul(&T::exact_f64(2.0));
                        let turn = if cross == zero() {
                            pi_t::<T>()
                        } else {
                            let angle = T::atan2(&T::from_r(&cross), &T::from_r(&dot))?;
                            if cross < zero() {
                                angle.add(&two_pi)
                            } else {
                                angle
                            }
                        };
                        let phi = if *ccw { turn } else { turn.sub(&two_pi) };
                        let lengths =
                            T::from_r(&((&ux * &ux + &uy * &uy) * (&vx * &vx + &vy * &vy))).sqrt();
                        let sine = T::from_r(&cross).div(&lengths)?;
                        let r = T::exact_f64(*radius);
                        twice = twice.add(&r.square().mul(&phi.sub(&sine)));
                        perimeter = perimeter.add(&abs(r.mul(&phi)));
                    }
                }
            }
            Some((abs(twice.mul(&T::exact_f64(0.5))), perimeter))
        }
    }
}

/// The screen in one tier: `None` when it cannot decide.
fn degenerate_in<T: Real>(outlines: &[Outline], tolerance: f64) -> Option<bool> {
    let mut area = T::exact_f64(0.0);
    let mut perimeter = T::exact_f64(0.0);
    for (i, outline) in outlines.iter().enumerate() {
        let (a, p) = area_perimeter::<T>(outline)?;
        area = if i == 0 { area.add(&a) } else { area.sub(&a) };
        perimeter = perimeter.add(&p);
    }
    let limit = perimeter.mul(&T::from_r(&(q(tolerance) * R::new(1.into(), 2.into()))));
    le(&area, &limit)
}

/// The material area (the first outline's minus the others') is at most the
/// tolerance times half the total perimeter: thinner than the tolerance on
/// average. Binary64 intervals first, rational intervals when those cannot
/// decide; undecidable in rational intervals counts as degenerate.
pub(crate) fn area_is_degenerate(outlines: &[Outline], tolerance: f64) -> bool {
    degenerate_in::<Fast>(outlines, tolerance)
        .or_else(|| degenerate_in::<Interval>(outlines, tolerance))
        .unwrap_or(true)
}

/// Where local coordinates `(x, y, z)` lie against the right circular cone or
/// frustum of radius `bottom` at z = 0 and `top` at z = height: 0 inside, 1
/// on the boundary within `tolerance`, 2 outside. Exact: the radius at z is
/// rational in the inputs, and radial distances compare squared. The band
/// on the lateral face is measured radially.
pub(crate) fn cone_location(
    point: [f64; 3],
    bottom: f64,
    top: f64,
    height: f64,
    tolerance: f64,
) -> u8 {
    let [x, y, z] = point.map(q);
    let (b, t, h, tol) = (q(bottom), q(top), q(height), q(tolerance));
    if z < -tol.clone() || z > &h + &tol {
        return 2;
    }
    let radius = &b + (&t - &b) * &z / &h;
    let radial2 = &x * &x + &y * &y;
    let outer = &radius + &tol;
    if outer < zero() || radial2 > &outer * &outer {
        return 2;
    }
    let near_end = |end: &R| {
        let d = &z - end;
        d <= tol.clone() && -d <= tol.clone()
    };
    let inner = &radius - &tol;
    if near_end(&zero()) || near_end(&h) || inner <= zero() || radial2 >= &inner * &inner {
        return 1;
    }
    0
}

/// Where a point (in the sphere's frame) lies against a sphere or zone of
/// `radius` between the axial heights `low` and `high` (at a pole, `-radius`
/// or `radius`): 0 inside, 1 on the boundary within `tolerance`, 2 outside.
/// Exact: distances compare squared.
pub(crate) fn sphere_location(
    point: [f64; 3],
    radius: f64,
    low: f64,
    high: f64,
    tolerance: f64,
) -> u8 {
    let [x, y, z] = point.map(q);
    let (r, lo, hi, tol) = (q(radius), q(low), q(high), q(tolerance));
    let d2 = &x * &x + &y * &y + &z * &z;
    let outer = &r + &tol;
    if d2 > &outer * &outer || z < &lo - &tol || z > &hi + &tol {
        return 2;
    }
    let near = |end: &R| {
        let d = &z - end;
        d <= tol.clone() && -d <= tol.clone()
    };
    let inner = &r - &tol;
    if near(&lo) || near(&hi) || inner <= zero() || d2 >= &inner * &inner {
        return 1;
    }
    0
}

/// Where a point (in the torus's frame) lies against a whole ring torus:
/// 0 inside, 1 on the boundary within `tolerance`, 2 outside. Exact: with
/// `rho` the distance from the axis, the distance to the core circle is at
/// most `d` exactly when `rho^2 + R^2 + z^2 - d^2 <= 2 R rho`, decided by
/// squaring (`d = r - tol`, `r + tol`).
pub(crate) fn torus_location(point: [f64; 3], major: f64, minor: f64, tolerance: f64) -> u8 {
    let [x, y, z] = point.map(q);
    let (big, r, tol) = (q(major), q(minor), q(tolerance));
    let rho2 = &x * &x + &y * &y;
    let base = &rho2 + &big * &big + &z * &z;
    // (rho - R)^2 + z^2 <= d^2  <=>  base - d^2 <= 2 R rho. On a ring torus
    // R > r + tol >= d, so base - d^2 > 0 and both sides square.
    let inside = |d: &R| -> bool {
        let lhs = &base - d * d;
        &lhs * &lhs <= R::from_integer(4.into()) * &big * &big * &rho2
    };
    let outer = &r + &tol;
    if !inside(&outer) {
        return 2;
    }
    let inner = &r - &tol;
    if inner <= zero() || !inside(&inner) {
        return 1;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ties_decide_exactly() {
        let (a, b) = (Point2::new(0.0, 0.0), Point2::new(3.0, 4.0));
        assert!(distance_le(a, b, &[5.0]));
        assert!(!distance_le(a, b, &[5.0 - f64::EPSILON * 4.0]));
        assert!(distance_ge(a, b, &[5.0]));
        assert!(!distance_ge(a, b, &[5.0 + f64::EPSILON * 8.0]));
        // 0.1 + 0.2 rounds up; the exact sum does not reach 0.3's value.
        assert!(sum_le(&[0.1, 0.2], &[0.30000000000000004]));
        assert!(sum_gt(&[0.1, 0.2], &[0.3]));
        let (p, s0, s1) = (
            Point2::new(1.0, 1e-9),
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
        );
        assert!(segment_distance_le(p, s0, s1, &[1e-9]));
        assert!(!segment_distance_le(p, s0, s1, &[1e-9 * (1.0 - 1e-15)]));
        assert!(segment_distance_le(Point2::new(3.0, 0.0), s0, s1, &[1.0]));
        assert!(!segment_distance_le(
            Point2::new(3.0, 0.0),
            s0,
            s1,
            &[0.5, 0.49]
        ));
    }

    #[test]
    fn thin_material_is_degenerate() {
        let sliver = [
            Point2::new(0.0, 0.0),
            Point2::new(10.0, 0.0),
            Point2::new(10.0, 1e-8),
            Point2::new(0.0, 1e-8),
        ];
        assert!(area_is_degenerate(&[Outline::Polygon(&sliver)], 1e-7));
        assert!(!area_is_degenerate(&[Outline::Polygon(&sliver)], 1e-9));
        assert!(!area_is_degenerate(
            &[Outline::Circle(2.0), Outline::Circle(1.0)],
            1e-7
        ));
    }
}
