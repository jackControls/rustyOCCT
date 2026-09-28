//! S9a.2's spline profiles: where a nonrational spline segment meets a
//! circle, another spline, or passes near a point (REVIEW_NOTES' S9a.2
//! spline decisions). A line's meetings are S8b.3's
//! (`solid::split::spline::meets`).
//!
//! On each exact Bézier arc of a spline, a circle's `(x - c_x)^2 + (y -
//! c_y)^2 - r^2` is a polynomial of degree `2 p` in the arc's parameter; its
//! roots come from the same isolation as a line's (S8b.3's `meets_with`).
//! Two arcs meet where each one's parameter is a root of the resultant of
//! the other's implicit equation on it: the implicit equation of an arc
//! `(x(t), y(t))` is `Res_t(x(t) - X, y(t) - Y)` (a Sylvester determinant),
//! its value on the other arc's points a polynomial of degree at most `p q`
//! in that arc's parameter, found exactly by exact evaluation at `p q + 1`
//! rational parameters and Newton interpolation. Roots of both resultants
//! in their arcs are paired by their certified points' boxes (refined until
//! each root has one partner or none); a root with no partner lies on the
//! other arc's algebraic curve but not on the arc.
use crate::certified::Interval as I;
use crate::decide::splines::to_f64;
use crate::solid::split::spline::{arcs, eval_interval, meets_with, power, roots, Meeting, Span};
use crate::solid::split::{q, zero};
use crate::{Error, Point2, Result};
use num_rational::BigRational as R;

fn one() -> R {
    R::from_integer(1.into())
}

/// The product of two polynomials in ascending powers.
fn mul(a: &[R], b: &[R]) -> Vec<R> {
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (k, y) in b.iter().enumerate() {
            out[i + k] += x * y;
        }
    }
    out
}

fn add(a: &[R], b: &[R]) -> Vec<R> {
    let mut out = vec![zero(); a.len().max(b.len())];
    for (i, x) in a.iter().enumerate() {
        out[i] += x;
    }
    for (i, x) in b.iter().enumerate() {
        out[i] += x;
    }
    out
}

/// A spline segment's meetings with a circle (`meets_with`'s crossings,
/// touches and sides).
pub(super) fn with_circle(span: &Span, center: Point2, radius: f64) -> Result<Meeting> {
    let (cx, cy, r) = (q(center.x), q(center.y), q(radius));
    let condition = move |cps: &[[R; 2]]| -> Vec<R> {
        let mut x = power(&cps.iter().map(|c| c[0].clone()).collect::<Vec<_>>());
        let mut y = power(&cps.iter().map(|c| c[1].clone()).collect::<Vec<_>>());
        x[0] -= &cx;
        y[0] -= &cy;
        let mut out = add(&mul(&x, &x), &mul(&y, &y));
        out[0] -= &r * &r;
        out
    };
    meets_with(span, &condition)
}

/// Whether two spline segments are one curve: equal degree, knots and
/// poles over their ranges, or the one's poles reversed and knots mirrored
/// (exactly) in the other's.
pub(super) fn same_curve(a: &Span, b: &Span) -> bool {
    let (ca, cb) = (a.curve().as_curve3(), b.curve().as_curve3());
    if ca.degree() != cb.degree() {
        return false;
    }
    let (pa, pb) = (a.curve().poles(), b.curve().poles());
    if ca.knots() == cb.knots()
        && ca.multiplicities() == cb.multiplicities()
        && pa == pb
        && a.range() == b.range()
    {
        return true;
    }
    let (lo, hi) = ca.domain();
    let mirror = |k: f64| q(lo) + q(hi) - q(k);
    let rev_poles: Vec<Point2> = pa.iter().rev().copied().collect();
    let mirrored_knots = ca
        .knots()
        .iter()
        .rev()
        .map(|k| mirror(*k))
        .eq(cb.knots().iter().map(|k| q(*k)));
    let mirrored_mults = ca
        .multiplicities()
        .iter()
        .rev()
        .eq(cb.multiplicities().iter());
    let [ra0, ra1] = a.range();
    let [rb0, rb1] = b.range();
    rev_poles == pb
        && mirrored_knots
        && mirrored_mults
        && mirror(ra1) == q(rb0)
        && mirror(ra0) == q(rb1)
}

/// A polynomial's degree (its last nonzero coefficient), none for zero.
fn degree(p: &[R]) -> Option<usize> {
    p.iter().rposition(|c| *c != zero())
}

/// The determinant of a square rational matrix (fraction-free Gaussian
/// elimination with row swaps).
fn det(mut m: Vec<Vec<R>>) -> R {
    let n = m.len();
    let mut sign = one();
    for col in 0..n {
        let Some(pivot) = (col..n).find(|&r| m[r][col] != zero()) else {
            return zero();
        };
        if pivot != col {
            m.swap(pivot, col);
            sign = -sign;
        }
        let pivot_row = m[col].clone();
        for row in m.iter_mut().skip(col + 1) {
            if row[col] == zero() {
                continue;
            }
            let f = &row[col] / &pivot_row[col];
            for (x, p) in row.iter_mut().zip(&pivot_row).skip(col) {
                *x -= &f * p;
            }
        }
    }
    let mut out = sign;
    for (i, row) in m.iter().enumerate() {
        out *= &row[i];
    }
    out
}

/// `Res_t(f, g)` of two polynomials in ascending powers, by their actual
/// degrees.
fn resultant(f: &[R], g: &[R]) -> R {
    let (Some(df), Some(dg)) = (degree(f), degree(g)) else {
        return zero();
    };
    if df == 0 {
        return (0..dg).fold(one(), |acc, _| acc * &f[0]);
    }
    if dg == 0 {
        return (0..df).fold(one(), |acc, _| acc * &g[0]);
    }
    let n = df + dg;
    let mut m = vec![vec![zero(); n]; n];
    for r in 0..dg {
        for k in 0..=df {
            m[r][r + k] = f[df - k].clone();
        }
    }
    for r in 0..df {
        for k in 0..=dg {
            m[dg + r][r + k] = g[dg - k].clone();
        }
    }
    det(m)
}

/// Arc `a`'s implicit equation on arc `b`'s points, a polynomial in `b`'s
/// parameter (ascending powers), by exact evaluation and interpolation.
fn implicit_on(a: &[[R; 2]], b: &[[R; 2]]) -> Vec<R> {
    let (ax, ay) = (
        power(&a.iter().map(|c| c[0].clone()).collect::<Vec<_>>()),
        power(&a.iter().map(|c| c[1].clone()).collect::<Vec<_>>()),
    );
    let (bx, by) = (
        power(&b.iter().map(|c| c[0].clone()).collect::<Vec<_>>()),
        power(&b.iter().map(|c| c[1].clone()).collect::<Vec<_>>()),
    );
    let da = degree(&ax)
        .unwrap_or(0)
        .max(degree(&ay).unwrap_or(0))
        .max(1);
    let db = degree(&bx)
        .unwrap_or(0)
        .max(degree(&by).unwrap_or(0))
        .max(1);
    let count = da * db + 1;
    let n = R::from_integer(((count - 1) as i64).into());
    let xs: Vec<R> = (0..count)
        .map(|k| R::from_integer((k as i64).into()) / &n)
        .collect();
    let value = |p: &[R], s: &R| p.iter().rev().fold(zero(), |acc, c| acc * s + c);
    let ys: Vec<R> = xs
        .iter()
        .map(|s| {
            let (x, y) = (value(&bx, s), value(&by, s));
            let mut f = ax.clone();
            let mut g = ay.clone();
            f[0] -= x;
            g[0] -= y;
            resultant(&f, &g)
        })
        .collect();
    // Newton's divided differences, then the power basis.
    let mut coef = ys;
    for j in 1..count {
        for i in (j..count).rev() {
            coef[i] = (&coef[i] - &coef[i - 1]) / (&xs[i] - &xs[i - j]);
        }
    }
    let mut out = vec![zero(); count];
    for i in (0..count).rev() {
        // out = out * (s - x_i) + coef_i
        let mut next = vec![zero(); count];
        for (k, c) in out.iter().enumerate() {
            if *c == zero() {
                continue;
            }
            if k + 1 < count {
                next[k + 1] += c;
            }
            next[k] -= c * &xs[i];
        }
        next[0] += &coef[i];
        out = next;
    }
    out
}

/// A crossing of two spline segments: each one's curve parameter (rounded
/// to binary64) and the first's exact point there.
#[derive(Debug, Clone)]
pub(super) struct SplineCrossing {
    pub(super) t: [f64; 2],
    pub(super) exact: [R; 2],
}

/// Where two spline segments cross strictly inside both (their own ends
/// are the vertices' business).
pub(super) fn with_spline(a: &Span, b: &Span) -> Result<Vec<SplineCrossing>> {
    let limit = || Error::ComputationLimit("a crossing of two splines, paired");
    let (arcs_a, arcs_b) = (arcs(a)?, arcs(b)?);
    let (ca, cb) = (a.curve(), b.curve());
    let (da, db) = (ca.as_curve3().domain(), cb.as_curve3().domain());
    let mut out: Vec<SplineCrossing> = Vec::new();
    for x in &arcs_a {
        for y in &arcs_b {
            // Screen: the control boxes apart.
            let bx = |c: &[[R; 2]], i: usize| {
                let lo = c.iter().map(|p| p[i].clone()).min().expect("a pole");
                let hi = c.iter().map(|p| p[i].clone()).max().expect("a pole");
                (lo, hi)
            };
            if (0..2).any(|i| {
                let ((l1, h1), (l2, h2)) = (bx(&x.cps, i), bx(&y.cps, i));
                h1 < l2 || h2 < l1
            }) {
                continue;
            }
            let on_y = implicit_on(&x.cps, &y.cps);
            let on_x = implicit_on(&y.cps, &x.cps);
            if on_y.iter().all(|c| *c == zero()) || on_x.iter().all(|c| *c == zero()) {
                return Err(Error::OutOfDomain(
                    "two spline segments along one curve of different forms (S9a.2)",
                ));
            }
            let mut sy = roots(&on_y)?;
            let mut sx = roots(&on_x)?;
            let point = |cps: &[[R; 2]], s: &I| -> [I; 2] {
                let px = power(&cps.iter().map(|c| c[0].clone()).collect::<Vec<_>>());
                let py = power(&cps.iter().map(|c| c[1].clone()).collect::<Vec<_>>());
                [eval_interval(&px, s), eval_interval(&py, s)]
            };
            let overlaps = |u: &[I; 2], v: &[I; 2]| {
                (0..2).all(|i| u[i].lo() <= v[i].hi() && v[i].lo() <= u[i].hi())
            };
            // Pair each root on `x` with the one on `y` at the same point.
            let mut pairs: Vec<(usize, usize)> = Vec::new();
            for (i, rx) in sx.iter_mut().enumerate() {
                let mut found = None;
                for _ in 0..40 {
                    let bx_ = point(
                        &x.cps,
                        &I::new(rx.isolator().0.clone(), rx.isolator().1.clone()),
                    );
                    let candidates: Vec<usize> = sy
                        .iter()
                        .enumerate()
                        .filter(|(_, ry)| {
                            let by_ = point(
                                &y.cps,
                                &I::new(ry.isolator().0.clone(), ry.isolator().1.clone()),
                            );
                            overlaps(&bx_, &by_)
                        })
                        .map(|(k, _)| k)
                        .collect();
                    match candidates.len() {
                        0 => {
                            found = Some(None);
                            break;
                        }
                        1 => {
                            // Unique while the boxes are tight: refine both
                            // until the box is small and the pairing stays.
                            let k = candidates[0];
                            let width = |r: &crate::polynomial::real::AlgebraicRoot| {
                                r.isolator().1 - r.isolator().0
                            };
                            let tight = R::new(1.into(), (1u64 << 40).into());
                            if width(rx) <= tight && width(&sy[k]) <= tight {
                                found = Some(Some(k));
                                break;
                            }
                            rx.refine_for_signs(8);
                            sy[k].refine_for_signs(8);
                        }
                        _ => {
                            rx.refine_for_signs(8);
                            for k in candidates {
                                sy[k].refine_for_signs(8);
                            }
                        }
                    }
                }
                match found {
                    Some(Some(k)) => pairs.push((i, k)),
                    Some(None) => {}
                    None => return Err(limit()),
                }
            }
            if pairs
                .iter()
                .map(|p| p.1)
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                != pairs.len()
            {
                return Err(limit());
            }
            for (i, k) in pairs {
                let ta = crate::solid::split::spline::rounded_parameter(&mut sx[i], &x.domain);
                let tb = crate::solid::split::spline::rounded_parameter(&mut sy[k], &y.domain);
                // The segments' own ends are vertices.
                let inside = |t: f64, d: (f64, f64), s: &Span| {
                    let [lo, hi] = s.range();
                    t > lo && t < hi && t > d.0 && t < d.1
                };
                if !inside(ta, da, a) || !inside(tb, db, b) {
                    continue;
                }
                // A tangency of two splines (S9a.2's decisions leave its
                // side undecided): refused.
                if sx[i].multiplicity() != 1 || sy[k].multiplicity() != 1 {
                    return Err(Error::Degenerate("a spline tangent to another"));
                }
                if out.iter().any(|c| c.t == [ta, tb]) {
                    // At an interior knot, from both arcs beside it.
                    continue;
                }
                let e = ca.as_curve3().exact_point(ta)?;
                out.push(SplineCrossing {
                    t: [ta, tb],
                    exact: [e[0].clone(), e[1].clone()],
                });
            }
        }
    }
    Ok(out)
}

/// The curve parameter (rounded) of a spline segment's point nearest `p`
/// within its range, and whether that point lies within `tol` of `p`
/// (decided exactly at the rounded parameter).
pub(super) fn nearest(span: &Span, p: Point2, tol: f64) -> Result<Option<f64>> {
    let (px, py) = (q(p.x), q(p.y));
    let [lo, hi] = span.range();
    let mut best: Option<(R, f64)> = None;
    for arc in arcs(span)? {
        let x = power(&arc.cps.iter().map(|c| c[0].clone()).collect::<Vec<_>>());
        let y = power(&arc.cps.iter().map(|c| c[1].clone()).collect::<Vec<_>>());
        let deriv = |p: &[R]| -> Vec<R> {
            p.iter()
                .enumerate()
                .skip(1)
                .map(|(k, c)| c * R::from_integer((k as i64).into()))
                .collect()
        };
        let (mut dx, mut dy) = (x.clone(), y.clone());
        dx[0] -= &px;
        dy[0] -= &py;
        // d/ds |B(s) - p|^2 / 2 = (B - p) . B'.
        let (vx, vy) = (deriv(&x), deriv(&y));
        let nothing = vec![zero()];
        let condition = if vx.is_empty() && vy.is_empty() {
            nothing.clone()
        } else {
            add(
                &mul(&dx, if vx.is_empty() { &nothing } else { &vx }),
                &mul(&dy, if vy.is_empty() { &nothing } else { &vy }),
            )
        };
        let mut candidates: Vec<f64> = Vec::new();
        if condition.iter().any(|c| *c != zero()) {
            for mut r in roots(&condition)? {
                candidates.push(crate::solid::split::spline::rounded_parameter(
                    &mut r,
                    &arc.domain,
                ));
            }
        }
        candidates.push(to_f64(&arc.domain[0]));
        candidates.push(to_f64(&arc.domain[1]));
        for t in candidates {
            if !(lo <= t && t <= hi) {
                continue;
            }
            let e = span.curve().as_curve3().exact_point(t)?;
            let d2 = (&e[0] - &px) * (&e[0] - &px) + (&e[1] - &py) * (&e[1] - &py);
            if best.as_ref().is_none_or(|b| d2 < b.0) {
                best = Some((d2, t));
            }
        }
    }
    let t2 = q(tol) * q(tol);
    Ok(best.filter(|b| b.0 <= t2).map(|b| b.1))
}
