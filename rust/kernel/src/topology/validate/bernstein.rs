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

pub(super) fn product<T: Real>(a: &Bern<T>, b: &Bern<T>) -> Bern<T> {
    let (m, n) = (a.len() - 1, b.len() - 1);
    (0..=m + n)
        .map(|k| {
            let total = R::from_integer(binomial(m + n, k));
            let mut sum = T::exact_f64(0.0);
            for i in k.saturating_sub(n)..=k.min(m) {
                let w = R::from_integer(binomial(m, i) * binomial(n, k - i)) / &total;
                sum = sum.add(&a[i].mul(&b[k - i]).mul(&c(&w)));
            }
            sum
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

pub(super) fn spline_arcs(curve: &crate::BSplineCurve3) -> Option<Arcs> {
    let exact = curve.to_exact();
    let [a, b] = exact.domain().clone();
    let span = &b - &a;
    let arcs = exact.bezier_arcs().ok()?;
    Some(
        arcs.into_iter()
            .map(|arc| {
                let [u0, u1] = arc.domain().clone();
                ((&u0 - &a) / &span, (&u1 - &a) / &span, arc)
            })
            .collect(),
    )
}

pub(super) fn curve_arcs(curve: &Curve3) -> Option<Arcs> {
    match curve {
        Curve3::LineSegment { start, end } => Some(line_arc(start.to_array(), end.to_array())),
        Curve3::BSpline(spline) => spline_arcs(spline),
        _ => None,
    }
}

pub(super) fn pcurve_arcs(p: &Curve2) -> Option<Arcs> {
    match p {
        Curve2::LineSegment { start, end } => {
            Some(line_arc([start.x, start.y, 0.0], [end.x, end.y, 0.0]))
        }
        Curve2::BSpline(spline) => spline_arcs(spline.as_curve3()),
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
    curve: &crate::BSplineCurve2,
    k: i32,
    integrand: impl Fn(&Bern<T>, &Bern<T>, &Bern<T>) -> Bern<T>,
) -> Option<T> {
    let arcs = spline_arcs(curve.as_curve3())?;
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
pub(super) fn twice_area<T: Real>(curve: &crate::BSplineCurve2) -> Option<T> {
    rational_integral(curve, 2, |u, v, _| {
        difference(&product(u, &derivative(v)), &product(v, &derivative(u)))
    })
}

/// `-∮ v du` of a spline pcurve: `-V (U' W - U W') / W^3`.
pub(super) fn minus_v_du<T: Real>(curve: &crate::BSplineCurve2) -> Option<T> {
    rational_integral(curve, 3, |u: &Bern<T>, v: &Bern<T>, w: &Bern<T>| {
        let du = difference(&product(&derivative(u), w), &product(u, &derivative(w)));
        product(v, &du).iter().map(|x| x.neg()).collect()
    })
}
