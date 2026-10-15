//! C1 continuity of spline cells (R4 of REVIEW_NOTES.md, with the rational
//! test of U10 (b)). A binary64 B-spline has interior multiplicities at
//! most its degree `p`, so it is C1 by construction at a knot of lower
//! multiplicity; at a knot of multiplicity `p` (a periodic seam included)
//! the rational curve `P = A / w` must be C1 there. One exact homogeneous
//! removal of that knot with zero residual (C1 of the homogeneous curve
//! `H = (A, w)`) decides it first, since it implies C1 of the rational
//! curve; where the removal leaves a residual the exact one-sided first
//! derivatives of the rational curve, `(A' w - A w') / w^2`, must be equal
//! (the jump of `H'` a multiple of `H`: the usual NURBS circle's quarter
//! arcs at double knots). Equal weights make the two tests one. A surface is
//! C1 across a knot line when every control row across it is homogeneously
//! C1, or else when its rational partial across the line agrees from both
//! sides along the whole line: the jump's numerator `ΔA_u w - A Δw_u` is a
//! polynomial of degree at most `2 q` on each span of the other parameter
//! (`q` its degree), so exact equality of the two one-sided partials at
//! `2 q + 1` distinct points inside each span decides it. Rows each C1 as
//! rational curves are not enough: the rows' blend mixes their multiples.
use crate::curve::{DerivativeOrder, ExactBSplineCurve3, KnotSide};
use crate::surface::ExactBSplineSurface3;
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

/// C1 of the homogeneous curve at `u`: one exact removal, down to
/// `degree - 1`, with zero residual. A small periodic basis may not allow
/// the removal's pole count: a simple knot inserted at the middle of the
/// longest span (the curve unchanged, `u` untouched) makes room first.
fn homogeneous_c1_at(curve: &ExactBSplineCurve3, u: &R) -> bool {
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

/// Whether every weight is the same: the rational test is then the
/// homogeneous one (`w' = 0` on both sides).
fn equal_weights(curve: &ExactBSplineCurve3) -> bool {
    let poles = curve.homogeneous_poles();
    poles.iter().all(|h| h[3] == poles[0][3])
}

/// C1 of the rational curve at `u` (U10 (b)): its exact first derivatives
/// from the left and from the right are equal (a periodic seam's left side
/// the end of the preceding period). Weights are positive, so `w > 0`.
fn rational_c1_at(curve: &ExactBSplineCurve3, u: &R) -> bool {
    let side = |side| curve.exact_evaluate(u, DerivativeOrder::First, side).ok();
    match (side(KnotSide::Left), side(KnotSide::Right)) {
        (Some(left), Some(right)) => left == right,
        _ => false,
    }
}

/// C1 at a tested knot: homogeneously, or else as a rational curve.
fn c1_at(curve: &ExactBSplineCurve3, u: &R) -> bool {
    homogeneous_c1_at(curve, u) || (!equal_weights(curve) && rational_c1_at(curve, u))
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

/// The rational surface's partial across the knot line at `k` of `axis`
/// agrees from both sides at `2 q + 1` distinct points inside each span of
/// the other parameter's domain (`q` its degree): exactly when the jump's
/// numerator, a polynomial of degree at most `2 q` there, vanishes.
fn rational_line_c1(surface: &ExactBSplineSurface3, axis: usize, k: &R) -> bool {
    let other = if axis == 0 {
        surface.v_knots()
    } else {
        surface.u_knots()
    };
    let [lo, hi] = other.domain();
    let mut ends: Vec<R> = other
        .knots()
        .iter()
        .filter(|t| *t > lo && *t < hi)
        .cloned()
        .collect();
    ends.insert(0, lo.clone());
    ends.push(hi.clone());
    let points = 2 * other.degree() + 1;
    let steps = R::from_integer((points + 1).into());
    let order = DerivativeOrder::First;
    let at = |t: &R, side: KnotSide| {
        let (u, v, sides) = if axis == 0 {
            (k, t, [side, KnotSide::Automatic])
        } else {
            (t, k, [KnotSide::Automatic, side])
        };
        surface.exact_evaluate(u, v, order, sides).ok()
    };
    ends.windows(2).all(|w| {
        (1..=points).all(|i| {
            let t = &w[0] + (&w[1] - &w[0]) * R::from_integer(i.into()) / &steps;
            matches!(
                (at(&t, KnotSide::Left), at(&t, KnotSide::Right)),
                (Some(left), Some(right)) if left == right
            )
        })
    })
}

/// C1 across every knot line of the surface, in both parameters.
pub(crate) fn surface_c1(surface: &BSplineSurface3) -> bool {
    let exact = surface.to_exact();
    let rational = surface.is_rational();
    (0..2).all(|axis| {
        let rows = exact.rows(axis);
        let Some(first) = rows.first() else {
            return true;
        };
        // Every knot of the whole surface.
        let [lo, hi] = first.domain().clone();
        tested(first, &lo, &hi, first.is_periodic())
            .iter()
            .all(|u| {
                rows.iter().all(|row| homogeneous_c1_at(row, u))
                    || (rational && rational_line_c1(&exact, axis, u))
            })
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{KnotVector, Point3};

    const W: f64 = std::f64::consts::FRAC_1_SQRT_2;

    fn circle(knots: Vec<f64>) -> BSplineCurve3 {
        let r = 1.5;
        let on = [[r, 0.0], [0.0, r], [-r, 0.0], [0.0, -r]];
        let corners = [[r, r], [-r, r], [-r, -r], [r, -r]];
        let poles = (0..4)
            .flat_map(|i| [on[i], corners[i]])
            .map(|[x, y]| Point3::new(x, y, 0.0))
            .collect();
        let mults = vec![2; knots.len()];
        BSplineCurve3::new_periodic(2, poles, Some([1.0, W].repeat(4)), knots, mults).unwrap()
    }

    fn line(weights: [f64; 5]) -> BSplineCurve3 {
        let poles = [0.0, 0.75, 1.5, 2.25, 3.0]
            .map(|x| Point3::new(x, 0.0, 0.0))
            .to_vec();
        BSplineCurve3::new(
            2,
            poles,
            Some(weights.to_vec()),
            vec![0.0, 0.5, 1.0],
            vec![3, 2, 3],
        )
        .unwrap()
    }

    fn whole(c: &BSplineCurve3) -> bool {
        let (a, b) = c.domain();
        curve_c1(c, [a, b], c.is_periodic())
    }

    #[test]
    fn the_nurbs_circle_is_c1_rationally_only() {
        let c = circle(vec![0.0, 1.0, 2.0, 3.0, 4.0]);
        assert!(whole(&c));
        let exact = c.to_exact();
        let seam = &exact.domain()[0];
        assert!(!homogeneous_c1_at(&exact, seam));
        assert!(rational_c1_at(&exact, seam));
        // Unequal spans: the speed halves at 3 and doubles at the seam (G1).
        assert!(!whole(&circle(vec![0.0, 1.0, 2.0, 3.0, 5.0])));
    }

    #[test]
    fn a_line_whose_weights_alone_jump_is_c1_and_a_speed_jump_is_not() {
        assert!(whole(&line([1.0, 1.0, 2.0, 1.0, 1.0])));
        assert!(!whole(&line([1.0, 1.0, 2.0, 3.0, 1.0])));
        // Equal weights other than one: the homogeneous verdict.
        assert!(whole(&line([2.0; 5])));
        let mut kinked = line([1.0; 5]);
        kinked = kinked
            .with_poles(
                kinked
                    .poles()
                    .iter()
                    .enumerate()
                    .map(|(i, p)| {
                        if i == 1 {
                            Point3::new(1.0, 0.0, 0.0)
                        } else {
                            *p
                        }
                    })
                    .collect(),
            )
            .unwrap();
        assert!(!whole(&kinked));
    }

    fn cap(rows: [[f64; 5]; 3]) -> BSplineSurface3 {
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for (i, fu) in [0.0, 0.25, 0.5, 0.75, 1.0].iter().enumerate() {
            for (j, fv) in [0.0, 0.5, 1.0].iter().enumerate() {
                poles.push(Point3::new(3.0 * fu, 2.0 * fv, 1.0));
                weights.push(rows[j][i]);
            }
        }
        BSplineSurface3::new(
            KnotVector::new(2, vec![0.0, 1.5, 3.0], vec![3, 2, 3]).unwrap(),
            KnotVector::new(2, vec![0.0, 2.0], vec![3, 3]).unwrap(),
            poles,
            Some(weights),
        )
        .unwrap()
    }

    #[test]
    fn a_surface_is_c1_by_its_rational_partial_along_the_whole_line() {
        let c1 = [1.0, 1.0, 2.0, 1.0, 1.0];
        let g1 = [1.0, 1.0, 2.0, 3.0, 1.0];
        assert!(surface_c1(&cap([c1; 3])));
        assert!(!surface_c1(&cap([g1; 3])));
        // Each row C1 as a rational curve, the surface not.
        let rows = cap([c1, [1.0; 5], c1]);
        let exact = rows.to_exact();
        let u = R::from_float(1.5).unwrap();
        assert!(exact.rows(0).iter().all(|row| c1_at(row, &u)));
        assert!(!surface_c1(&rows));
        // Product weights across v: the rows' common factor cancels.
        let scaled = |k: f64| c1.map(|w| w * k);
        assert!(surface_c1(&cap([scaled(1.0), scaled(3.0), scaled(0.5)])));
    }
}
