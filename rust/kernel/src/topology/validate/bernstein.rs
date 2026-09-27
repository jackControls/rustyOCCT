//! Bernstein polynomials on `[0, 1]` over either arithmetic tier, and exact
//! Bézier pieces of the curves and surfaces of a topology (S4 of
//! REVIEW_NOTES.md). Pieces are extracted exactly; the arithmetic on them
//! runs in binary64 intervals first and exactly when those cannot decide.
use crate::certified::Real;
use crate::curve::ExactBezierCurve3;
use crate::topology::{Curve2, Curve3};
use num_bigint::BigInt;
use num_rational::BigRational as R;

pub(super) type Bern<T> = Vec<T>;

pub(super) fn exact(x: i64) -> R {
    R::from_integer(BigInt::from(x))
}

pub(super) fn ratio(a: i64, b: i64) -> R {
    R::new(BigInt::from(a), BigInt::from(b))
}

pub(super) fn r(x: f64) -> R {
    R::from_float(x).expect("finite geometry")
}

pub(super) fn c<T: Real>(x: &R) -> T {
    T::from_r(x)
}

fn binomial(n: usize, k: usize) -> BigInt {
    let mut b = BigInt::from(1);
    for i in 0..k {
        b = b * BigInt::from(n - i) / BigInt::from(i + 1);
    }
    b
}

pub(super) fn lift<T: Real>(a: &[R]) -> Bern<T> {
    a.iter().map(c).collect()
}

/// `C(n, k)` in the tier, exactly: binomials up to degree 56 are binary64
/// integers, larger ones pass through rationals.
fn binomial_t<T: Real>(n: usize, k: usize) -> T {
    if n <= 56 {
        let mut x = 1.0_f64;
        for i in 0..k {
            x = x * (n - i) as f64 / (i + 1) as f64;
        }
        T::exact_f64(x.round())
    } else {
        c(&R::from_integer(binomial(n, k)))
    }
}

/// The product in the scaled basis `C(n, i) t^i (1 - t)^(n - i)`, where it
/// is a convolution, converted back by `1/C(m + n, k)`.
pub(super) fn product<T: Real>(a: &Bern<T>, b: &Bern<T>) -> Bern<T> {
    let (m, n) = (a.len() - 1, b.len() - 1);
    let sa: Vec<T> = a
        .iter()
        .enumerate()
        .map(|(i, x)| x.mul(&binomial_t(m, i)))
        .collect();
    let sb: Vec<T> = b
        .iter()
        .enumerate()
        .map(|(j, x)| x.mul(&binomial_t(n, j)))
        .collect();
    (0..=m + n)
        .map(|k| {
            let mut sum = T::exact_f64(0.0);
            for i in k.saturating_sub(n)..=k.min(m) {
                sum = sum.add(&sa[i].mul(&sb[k - i]));
            }
            sum.div(&binomial_t(m + n, k)).expect("a positive binomial")
        })
        .collect()
}

pub(super) fn elevate<T: Real>(a: &Bern<T>, degree: usize) -> Bern<T> {
    let mut a = a.clone();
    while a.len() - 1 < degree {
        let n = a.len() as i64;
        let mut next = Vec::with_capacity(a.len() + 1);
        next.push(a[0].clone());
        for i in 1..a.len() {
            let f = c::<T>(&ratio(i as i64, n));
            let g = c::<T>(&ratio(n - i as i64, n));
            next.push(f.mul(&a[i - 1]).add(&g.mul(&a[i])));
        }
        next.push(a[a.len() - 1].clone());
        a = next;
    }
    a
}

pub(super) fn sum<T: Real>(a: &Bern<T>, b: &Bern<T>) -> Bern<T> {
    let n = a.len().max(b.len()) - 1;
    let (a, b) = (elevate(a, n), elevate(b, n));
    a.iter().zip(&b).map(|(x, y)| x.add(y)).collect()
}

pub(super) fn scaled<T: Real>(a: &Bern<T>, s: &T) -> Bern<T> {
    a.iter().map(|x| x.mul(s)).collect()
}

pub(super) fn difference<T: Real>(a: &Bern<T>, b: &Bern<T>) -> Bern<T> {
    sum(a, &b.iter().map(|x| x.neg()).collect())
}

pub(super) fn power<T: Real>(a: &Bern<T>, k: usize) -> Bern<T> {
    (0..k).fold(vec![T::exact_f64(1.0)], |acc, _| product(&acc, a))
}

/// The derivative, of one degree less (a constant's is zero).
pub(super) fn derivative<T: Real>(a: &Bern<T>) -> Bern<T> {
    let n = a.len() - 1;
    if n == 0 {
        return vec![T::exact_f64(0.0)];
    }
    let scale = T::exact_f64(n as f64);
    a.windows(2).map(|w| w[1].sub(&w[0]).mul(&scale)).collect()
}

/// The integral over `[0, 1]`: the mean of the coefficients.
pub(super) fn integral<T: Real>(a: &Bern<T>) -> T {
    let total = a.iter().fold(T::exact_f64(0.0), |s, x| s.add(x));
    total.mul(&c(&ratio(1, a.len() as i64)))
}

/// De Casteljau halves.
pub(super) fn halves<T: Real>(a: &Bern<T>) -> (Bern<T>, Bern<T>) {
    let half = T::exact_f64(0.5);
    let mut row = a.clone();
    let (mut left, mut right) = (vec![row[0].clone()], vec![row[row.len() - 1].clone()]);
    while row.len() > 1 {
        row = row.windows(2).map(|w| w[0].add(&w[1]).mul(&half)).collect();
        left.push(row[0].clone());
        right.push(row[row.len() - 1].clone());
    }
    right.reverse();
    (left, right)
}

pub(super) fn value<T: Real>(a: &Bern<T>, t: &T) -> T {
    let one = T::exact_f64(1.0);
    let s = one.sub(t);
    let mut row = a.clone();
    while row.len() > 1 {
        row = row
            .windows(2)
            .map(|w| s.mul(&w[0]).add(&t.mul(&w[1])))
            .collect();
    }
    row.pop().unwrap()
}

/// A curve as exact Bézier arcs over `[0, 1]` of its fraction: `(from, to,
/// arc)` with the arc on its own domain.
pub(super) type Arcs = Vec<(R, R, ExactBezierCurve3)>;

fn line_arc(a: [f64; 3], b: [f64; 3]) -> Arcs {
    let controls = vec![
        [r(a[0]), r(a[1]), r(a[2]), exact(1)],
        [r(b[0]), r(b[1]), r(b[2]), exact(1)],
    ];
    let arc = ExactBezierCurve3::from_homogeneous(controls, [exact(0), exact(1)]);
    vec![(exact(0), exact(1), arc)]
}

thread_local! {
    /// Arcs of the splines of the topology being validated, by the curve's
    /// address and range: a topology is borrowed immutably while it is
    /// validated, so addresses are stable; `clear_memo` runs at every entry.
    static ARCS: std::cell::RefCell<std::collections::HashMap<(usize, u64, u64), Arcs>> =
        std::cell::RefCell::new(std::collections::HashMap::new());
}

/// Forget the memoized arcs (at the start of a validation, measurement or
/// mass computation).
pub(super) fn clear_memo() {
    ARCS.with(|m| m.borrow_mut().clear());
}

pub(super) fn spline_arcs(curve: &crate::BSplineCurve3, range: [f64; 2]) -> Option<Arcs> {
    let key = (
        curve as *const crate::BSplineCurve3 as usize,
        range[0].to_bits(),
        range[1].to_bits(),
    );
    if let Some(arcs) = ARCS.with(|m| m.borrow().get(&key).cloned()) {
        return Some(arcs);
    }
    let arcs = compute_arcs(curve, range)?;
    ARCS.with(|m| m.borrow_mut().insert(key, arcs.clone()));
    Some(arcs)
}

fn compute_arcs(curve: &crate::BSplineCurve3, range: [f64; 2]) -> Option<Arcs> {
    let exact = curve.to_exact();
    let (a, b) = (r(range[0]), r(range[1]));
    let span = &b - &a;
    let arcs = exact.bezier_arcs_in(&a, &b).ok()?;
    Some(
        arcs.into_iter()
            .map(|arc| {
                let [u0, u1] = arc.domain().clone();
                ((&u0 - &a) / &span, (&u1 - &a) / &span, arc)
            })
            .collect(),
    )
}

/// A spline span's exact arcs over its range, in its direction: a reversed
/// span's arcs reversed, in reverse order, on mirrored fractions.
fn directed(arcs: Arcs, reversed: bool) -> Arcs {
    if !reversed {
        return arcs;
    }
    let one = exact(1);
    arcs.into_iter()
        .rev()
        .map(|(a, b, arc)| (&one - b, &one - a, arc.reversed()))
        .collect()
}

/// A spline pcurve's exact arcs over its range, in its direction.
pub(super) fn span_arcs(span: &crate::topology::SplineSpan<crate::BSplineCurve2>) -> Option<Arcs> {
    Some(directed(
        spline_arcs(span.curve().as_curve3(), span.range())?,
        span.is_reversed(),
    ))
}

/// A spline edge's exact arcs over its range, in its direction.
pub(super) fn edge_arcs(span: &crate::topology::SplineSpan<crate::BSplineCurve3>) -> Option<Arcs> {
    Some(directed(
        spline_arcs(span.curve(), span.range())?,
        span.is_reversed(),
    ))
}

pub(super) fn curve_arcs(curve: &Curve3) -> Option<Arcs> {
    match curve {
        Curve3::LineSegment { start, end } => Some(line_arc(start.to_array(), end.to_array())),
        Curve3::BSpline(span) => edge_arcs(span),
        _ => None,
    }
}

pub(super) fn pcurve_arcs(p: &Curve2) -> Option<Arcs> {
    match p {
        Curve2::LineSegment { start, end } => {
            Some(line_arc([start.x, start.y, 0.0], [end.x, end.y, 0.0]))
        }
        Curve2::BSpline(span) => span_arcs(span),
        Curve2::CircularArc { .. } => None,
    }
}

/// The homogeneous Bernstein coordinates of the arc covering fractions
/// `[f0, f1]`, reparameterized onto `[0, 1]`.
pub(super) fn piece_of(arcs: &Arcs, f0: &R, f1: &R) -> Option<[Vec<R>; 4]> {
    let (from, to, arc) = arcs.iter().find(|(from, to, _)| from <= f0 && f1 <= to)?;
    let [u0, u1] = arc.domain().clone();
    let at = |f: &R| &u0 + (f - from) / (to - from) * (&u1 - &u0);
    let trimmed = arc.trim(&at(f0), &at(f1)).ok()?;
    let poles = trimmed.homogeneous_poles();
    Some(std::array::from_fn(|k| {
        poles.iter().map(|p| p[k].clone()).collect()
    }))
}

/// The lower and upper bounds of an enclosure, as rationals.
pub(super) fn ends<T: Real>(x: &T) -> (R, R) {
    let (m, w) = (x.midpoint(), x.radius());
    (&m - &w, m + w)
}

/// `∫ M/W^k` over `[0, 1]`, enclosed: exact where every coefficient of `W`
/// is the same exact value; elsewhere, on pieces halved six times, with
/// `1/W^k` in `[a, b]`, within `a ∫M ± (b - a) max |M_i|`. `None` when `W`
/// is not certainly positive.
pub(super) fn quotient_integral<T: Real>(
    m: &Bern<T>,
    w: &Bern<T>,
    k: i32,
    uniform: bool,
) -> Option<T> {
    let unit = T::exact_f64(0.0).widen(&exact(1));
    let mut total = T::exact_f64(0.0);
    let mut stack = vec![(m.clone(), w.clone(), if uniform { 0 } else { 6 })];
    while let Some((m, w, depth)) = stack.pop() {
        if depth > 0 {
            let (m, w) = (halves(&m), halves(&w));
            // Each half covers half the parameter.
            let half = T::exact_f64(0.5);
            stack.push((scaled(&m.0, &half), w.0, depth - 1));
            stack.push((scaled(&m.1, &half), w.1, depth - 1));
            continue;
        }
        let inverse: Vec<T> = w
            .iter()
            .map(|x| {
                let power = (0..k).fold(T::exact_f64(1.0), |acc, _| acc.mul(x));
                (x.sign() == Some(std::cmp::Ordering::Greater))
                    .then(|| T::exact_f64(1.0).div(&power))
                    .flatten()
            })
            .collect::<Option<_>>()?;
        let bounds: Vec<(R, R)> = inverse.iter().map(ends).collect();
        let a = bounds.iter().map(|b| b.0.clone()).min()?;
        let b = bounds.iter().map(|b| b.1.clone()).max()?;
        let largest = m
            .iter()
            .map(|x| {
                let (lo, hi) = ends(x);
                if -&lo > hi {
                    -lo
                } else {
                    hi
                }
            })
            .max()?;
        let spread = c::<T>(&((b - &a) * largest)).mul(&unit);
        total = total.add(&integral(&m).mul(&c(&a)).add(&spread));
    }
    Some(total)
}

/// `∫ M/W^k` over a spline pcurve, `M` built from the homogeneous Bernstein
/// coordinates `(U, V, W)` of each Bézier piece.
pub(super) fn rational_integral<T: Real>(
    curve: &crate::topology::SplineSpan<crate::BSplineCurve2>,
    k: i32,
    integrand: impl Fn(&Bern<T>, &Bern<T>, &Bern<T>) -> Bern<T>,
) -> Option<T> {
    let arcs = span_arcs(curve)?;
    let mut total = T::exact_f64(0.0);
    for (_, _, arc) in &arcs {
        let poles = arc.homogeneous_poles();
        let uniform = poles.iter().all(|p| p[3] == poles[0][3]);
        let [u, v, w]: [Bern<T>; 3] = [0, 1, 3].map(|k| poles.iter().map(|p| c(&p[k])).collect());
        total = total.add(&quotient_integral(&integrand(&u, &v, &w), &w, k, uniform)?);
    }
    Some(total)
}

/// Twice the signed area `∮ (u dv - v du)` of a spline pcurve: on a piece
/// `u = U/W`, `v = V/W`, the integrand is `(U V' - V U')/W^2`.
pub(super) fn twice_area<T: Real>(
    curve: &crate::topology::SplineSpan<crate::BSplineCurve2>,
) -> Option<T> {
    rational_integral(curve, 2, |u, v, _| {
        difference(&product(u, &derivative(v)), &product(v, &derivative(u)))
    })
}

/// `∮ u dv` of a spline pcurve: `U (V' W - V W') / W^3`.
pub(super) fn u_dv<T: Real>(
    curve: &crate::topology::SplineSpan<crate::BSplineCurve2>,
) -> Option<T> {
    rational_integral(curve, 3, |u: &Bern<T>, v: &Bern<T>, w: &Bern<T>| {
        let dv = difference(&product(&derivative(v), w), &product(v, &derivative(w)));
        product(u, &dv)
    })
}

/// `-∮ v du` of a spline pcurve: `-V (U' W - U W') / W^3`.
pub(super) fn minus_v_du<T: Real>(
    curve: &crate::topology::SplineSpan<crate::BSplineCurve2>,
) -> Option<T> {
    rational_integral(curve, 3, |u: &Bern<T>, v: &Bern<T>, w: &Bern<T>| {
        let du = difference(&product(&derivative(u), w), &product(u, &derivative(w)));
        product(v, &du).iter().map(|x| x.neg()).collect()
    })
}

/// The parity of the crossings of the ray `u > p.u`, `v = p.v` with a spline
/// pcurve, counted half-open at `v = p.v` as for segments (`above` is
/// `v > p.v`). A piece certainly right of the point contributes whether its
/// ends lie on different sides, whatever its shape; a piece certainly left,
/// above or below contributes nothing; any other piece is halved exactly,
/// up to twelve times. `None` when a piece stays undecided (the pcurve
/// passes through or next to the point).
pub(super) fn crossing_parity<T: Real>(
    curve: &crate::topology::SplineSpan<crate::BSplineCurve2>,
    p: &[T; 2],
) -> Option<u32> {
    use std::cmp::Ordering::{Greater, Less};
    let arcs = span_arcs(curve)?;
    let above = |v: &R| -> Option<bool> { Some(c::<T>(v).sub(&p[1]).sign()? == Greater) };
    let mut parity = 0;
    for (_, _, arc) in &arcs {
        let h = arc.homogeneous_poles();
        let piece: [Vec<R>; 3] = [0, 1, 3].map(|k| h.iter().map(|q| q[k].clone()).collect());
        let mut stack = vec![(piece, 0)];
        while let Some(([u, v, w], depth)) = stack.pop() {
            let points: Vec<(R, R)> = (0..w.len())
                .map(|i| (&u[i] / &w[i], &v[i] / &w[i]))
                .collect();
            // Certain sides of the control hull, and so of the piece.
            let all = |f: &dyn Fn(&(R, R)) -> Option<bool>| -> Option<bool> {
                let mut out = true;
                for q in &points {
                    out &= f(q)?;
                }
                Some(out)
            };
            let side_u =
                |q: &(R, R), want: std::cmp::Ordering| c::<T>(&q.0).cmp(&p[0]).map(|o| o == want);
            let side_v =
                |q: &(R, R), want: std::cmp::Ordering| c::<T>(&q.1).cmp(&p[1]).map(|o| o == want);
            let left = all(&|q| side_u(q, Less)).unwrap_or(false);
            let over = all(&|q| side_v(q, Greater)).unwrap_or(false);
            let under = all(&|q| side_v(q, Less)).unwrap_or(false);
            if left || over || under {
                continue;
            }
            if all(&|q| side_u(q, Greater)).unwrap_or(false) {
                let first = above(&points[0].1)?;
                let last = above(&points[points.len() - 1].1)?;
                parity ^= u32::from(first != last);
                continue;
            }
            if depth == 12 {
                return None;
            }
            let halves3 = |x: &Vec<R>| {
                let mut row = x.clone();
                let two = exact(2);
                let (mut l, mut r) = (vec![row[0].clone()], vec![row[row.len() - 1].clone()]);
                while row.len() > 1 {
                    row = row.windows(2).map(|w| (&w[0] + &w[1]) / &two).collect();
                    l.push(row[0].clone());
                    r.push(row[row.len() - 1].clone());
                }
                r.reverse();
                (l, r)
            };
            let (ua, ub) = halves3(&u);
            let (va, vb) = halves3(&v);
            let (wa, wb) = halves3(&w);
            stack.push(([ua, va, wa], depth + 1));
            stack.push(([ub, vb, wb], depth + 1));
        }
    }
    Some(parity)
}

/// `∫ F(u, v) du` along a spline pcurve, enclosed, for a function `F` the
/// caller evaluates over boxes (the antiderivative in `v` of a surface's
/// flux or mass integrand). Each Bézier piece is halved `depth` times; on a
/// piece in its own parameter `τ` in `[0, 1]`, the integrand `F(P) P_u'`
/// is enclosed over the box of the piece's control points (its convex
/// hull) with `P_u' = (U' W - U W')/W^2`, so its integral lies in that
/// enclosure. The width is `O(h^2)` per piece, `O(h)` in total.
pub(super) fn green_integral<T: Real>(
    curve: &crate::topology::SplineSpan<crate::BSplineCurve2>,
    depth: usize,
    f: &dyn Fn(&T, &T) -> Option<T>,
) -> Option<T> {
    let hull = |xs: &[T]| xs[1..].iter().fold(xs[0].clone(), |acc, x| acc.union(x));
    let mut total = T::exact_f64(0.0);
    for (_, _, arc) in &span_arcs(curve)? {
        let poles = arc.homogeneous_poles();
        let coords: [Bern<T>; 3] = [0, 1, 3].map(|k| poles.iter().map(|p| c(&p[k])).collect());
        let mut stack = vec![(coords, depth)];
        while let Some(([u, v, w], d)) = stack.pop() {
            if d > 0 {
                let (u, v, w) = (halves(&u), halves(&v), halves(&w));
                // Each half covers half the piece's parameter.
                stack.push(([u.0, v.0, w.0], d - 1));
                stack.push(([u.1, v.1, w.1], d - 1));
                continue;
            }
            let points: Vec<[T; 2]> = (0..w.len())
                .map(|i| Some([u[i].div(&w[i])?, v[i].div(&w[i])?]))
                .collect::<Option<_>>()?;
            let bu = hull(&points.iter().map(|p| p[0].clone()).collect::<Vec<_>>());
            let bv = hull(&points.iter().map(|p| p[1].clone()).collect::<Vec<_>>());
            let wb = hull(&w);
            if wb.sign() != Some(std::cmp::Ordering::Greater) {
                return None;
            }
            let du = difference(&product(&derivative(&u), &w), &product(&u, &derivative(&w)));
            let du = hull(&du).div(&wb.square())?;
            // The piece is 2^-(depth - d) of the arc; P_u' is in the piece's
            // own parameter, so the enclosure is the piece's integral.
            total = total.add(&f(&bu, &bv)?.mul(&du));
        }
    }
    Some(total)
}
