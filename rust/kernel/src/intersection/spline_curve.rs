//! Rational B-splines against circles, ellipses and hyperbolas' branches
//! (S7d.2 of `REVIEW_NOTES.md`).
//!
//! On each span, with the homogeneous polynomials `(X, W)`, the conic's
//! plane gives `L = n . (X - o W)` (degree `p`) and its equation `Q = W^2 E(X
//! / W)` (degree `2 p`), both exact. A span lies on the conic exactly when `L
//! = Q = 0`; otherwise the points are the real roots of `gcd(L, Q)` in the
//! span, a hyperbola's branch kept by a certified sign. A root's order in the
//! gcd is its least order in `L` and `Q`: at least two puts the spline's
//! tangent in the plane and along the conic (a tangency).
use super::curve_curve::{AnalyticCurve, CurveCurvePoint, Exact};
use super::procedural::{bounds, bounds3, limit, E};
use super::tangency::P;
use crate::certified::{Interval as I, Real};
use crate::polynomial::real::{isolate, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::{BSplineCurve3, Error, Result};
use num_rational::BigRational as R;

/// A rational B-spline against a conic: isolated points (the spline's
/// parameter first) and overlaps.
#[derive(Debug, Clone, PartialEq)]
pub struct SplineCurveIntersection {
    /// Sorted by the spline's parameter, none inside an overlap.
    pub points: Vec<CurveCurvePoint>,
    /// Maximal runs of whole spans on the conic, between exact knots.
    pub overlaps: Vec<[f64; 2]>,
}

/// The whole closed domain of a rational B-spline against a circle (an
/// arc's whole circle), an ellipse or a hyperbola's branch. A line is
/// `spline_line`'s.
pub fn spline_curve(
    curve: &BSplineCurve3,
    other: &AnalyticCurve,
) -> Result<SplineCurveIntersection> {
    let conic = match Exact::of(other)? {
        Exact::Plane(c) => c,
        Exact::Line(_) => return Err(Error::OutOfDomain("a spline against a line (spline_line)")),
    };
    let (first, last) = curve.domain();
    let spans = curve
        .knot_vector()
        .spans_in(first, last, crate::spline::MAX_POLES)?;
    let mut points: Vec<CurveCurvePoint> = Vec::new();
    let mut overlaps: Vec<[R; 2]> = Vec::new();
    let mut budget = Budget::new(RootIsolationOptions::default());
    for at in spans {
        let length = &at.end - &at.start;
        let lower = (&at.lower - &at.start) / &length;
        let upper = (&at.upper - &at.start) / &length;
        let h = curve.span_polynomial(at.index);
        let (l, q) = conic.homogeneous(&h);
        let (il, iq) = (
            IntPolynomial::from_rationals(&l),
            IntPolynomial::from_rationals(&q),
        );
        if il.is_zero() && iq.is_zero() {
            match overlaps.last_mut().filter(|o| o[1] == at.lower) {
                Some(last) => last[1] = at.upper.clone(),
                None => overlaps.push([at.lower.clone(), at.upper.clone()]),
            }
            continue;
        }
        let common = il.gcd(&iq);
        if common.is_constant() {
            continue;
        }
        for mut root in isolate(&common, lower.clone(), upper.clone(), &mut budget)? {
            let tangent = root.multiplicity() > 1;
            root.refine_for_signs(300);
            let (lo, hi) = root.isolator();
            let x = I::new(lo.clone(), hi.clone());
            let point = point_at(&h, &x)?;
            if !conic.on_branch(&point)? {
                continue;
            }
            let s = x
                .mul(&I::exact(length.clone()))
                .add(&I::exact(at.start.clone()));
            points.push(CurveCurvePoint {
                parameters: [bounds(&s), conic.parameter(&point)?],
                point: bounds3(&point),
                tangent,
            });
        }
    }
    let overlaps: Vec<[f64; 2]> = overlaps
        .iter()
        .map(|[a, b]| {
            [
                bounds(&I::exact(a.clone()))[0],
                bounds(&I::exact(b.clone()))[1],
            ]
        })
        .collect();
    // A root at a knot comes from both spans (exact, so equal enclosures);
    // points at an overlap's ends belong to it.
    points.sort_by(|a, b| a.parameters[0][0].total_cmp(&b.parameters[0][0]));
    points.dedup_by(|b, a| a.parameters[0] == b.parameters[0]);
    points.retain(|p| {
        !overlaps
            .iter()
            .any(|[lo, hi]| *lo <= p.parameters[0][0] && p.parameters[0][1] <= *hi)
    });
    Ok(SplineCurveIntersection { points, overlaps })
}

/// `X(x) / W(x)` of a span's power-basis polynomials, enclosed.
fn point_at(h: &[P; 4], x: &I) -> Result<E<I>> {
    let v: [I; 4] = std::array::from_fn(|c| {
        h[c].iter().rev().fold(I::exact_f64(0.0), |acc, k| {
            acc.mul(x).add(&I::exact(k.clone()))
        })
    });
    let mut p: E<I> = std::array::from_fn(|_| I::exact_f64(0.0));
    for c in 0..3 {
        p[c] = v[c].div(&v[3]).ok_or(limit("a spline's point"))?;
    }
    Ok(p)
}
