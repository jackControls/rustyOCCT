//! C1 continuity of spline cells (R4 of REVIEW_NOTES.md). A binary64
//! B-spline has interior multiplicities at most its degree `p`, so it is C1
//! by construction at a knot of lower multiplicity; at a knot of
//! multiplicity `p` (a periodic seam included) it is C1 exactly when one
//! exact homogeneous removal of that knot has zero residual. That is C1 of
//! the homogeneous curve, which implies C1 of the rational one. A surface is
//! C1 across a knot line when every control row across it is.
use crate::curve::ExactBSplineCurve3;
use crate::{BSplineCurve3, BSplineSurface3};
use num_rational::BigRational as R;

/// The knots of multiplicity equal to the degree that a C1 test must
/// remove: those strictly inside the range `[first, last]` (an unclamped
/// basis has knots beyond its domain, a trimmed edge beyond its range), and
/// on a closed period also the knot at the range's ends, where the ring
/// joins. A periodic knot counts at every translate by the period, reduced
/// to the base period for the removal.
fn tested(curve: &ExactBSplineCurve3, first: &R, last: &R, closed: bool) -> Vec<R> {
    let (knots, mults) = (curve.knots(), curve.multiplicities());
    let [a, b] = curve.domain();
    let n = knots.len();
    let mut out = Vec::new();
    for i in 0..n {
        if mults[i] < curve.degree() {
            continue;
        }
        let k = &knots[i];
        if !curve.is_periodic() {
            if k > first && k < last {
                out.push(k.clone());
            }
            continue;
        }
        if i + 1 == n {
            continue; // the seam's copy
        }
        let period = b - a;
        // The translates of k in [first, last].
        let mut t = k - ((k - first) / &period).floor() * &period;
        while t <= *last {
            if (t > *first && t < *last) || (closed && (t == *first || t == *last)) {
                out.push(k.clone());
                break;
            }
            t += &period;
        }
    }
    out
}

/// One exact removal at `u`, down to `degree - 1`. A small periodic basis
/// may not allow the removal's pole count: a simple knot inserted at the
/// middle of the longest span (the curve unchanged, `u` untouched) makes
/// room first.
fn c1_at(curve: &ExactBSplineCurve3, u: &R) -> bool {
    let target = curve.degree() - 1;
    if let Ok(removable) = curve.removable(u, target) {
        return removable;
    }
    let knots = curve.knots();
    let Some((a, b)) = knots
        .windows(2)
        .map(|w| (&w[0], &w[1]))
        .max_by(|x, y| (x.1 - x.0).cmp(&(y.1 - y.0)))
    else {
        return false;
    };
    let middle = (a + b) / R::from_integer(2.into());
    curve
        .insert_knot(&middle, 1)
        .and_then(|refined| refined.removable(u, target))
        .unwrap_or(false)
}

/// C1 over the range, in the curve's own parameter.
pub(crate) fn curve_c1(curve: &BSplineCurve3, range: [f64; 2], closed: bool) -> bool {
    let exact = curve.to_exact();
    let first = R::from_float(range[0]).expect("finite range");
    let last = R::from_float(range[1]).expect("finite range");
    tested(&exact, &first, &last, closed)
        .iter()
        .all(|u| c1_at(&exact, u))
}

/// C1 across every knot line of the surface, in both parameters.
pub(crate) fn surface_c1(surface: &BSplineSurface3) -> bool {
    let exact = surface.to_exact();
    (0..2).all(|axis| {
        let rows = exact.rows(axis);
        let Some(first) = rows.first() else {
            return true;
        };
        // Every knot of the whole surface.
        let [lo, hi] = first.domain().clone();
        tested(first, &lo, &hi, first.is_periodic())
            .iter()
            .all(|u| rows.iter().all(|row| c1_at(row, u)))
    })
}
