//! S8b.3: a spline profile segment against a plane's trace.
//!
//! The line's function `a x + b y + d` on each exact Bézier arc of a
//! nonrational spline is a polynomial with rational coefficients in the
//! arc's parameter; its real roots are isolated exactly with their
//! multiplicities (`polynomial::real::isolate`). A simple root is a
//! crossing, at the curve parameter rounded to binary64 from its isolator;
//! an even one a touch; any other, and a tangency at a knot, is a plane
//! tangent to the profile where it crosses it (`Degenerate`). Sides come
//! from exact signs at rational parameters between distinct roots, never
//! from the rounded crossings. A piece of the segment is the exact
//! restriction of its curve to the rounded parameters (Boehm's knot
//! insertion in rationals), its poles rounded: its ends are the rounded
//! crossing points exactly.
use super::{q, zero};
use crate::certified::Interval as I;
use crate::decide::splines::to_f64;
use crate::polynomial::real::{isolate, AlgebraicRoot, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::topology::SplineSpan;
use crate::{BSplineCurve2, Error, Point2, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

pub(crate) type Span = SplineSpan<BSplineCurve2>;

fn one() -> R {
    R::from_integer(1.into())
}

/// An exact Bézier arc of a spline: its control points and its range of
/// the curve's parameter.
pub(crate) struct Arc {
    pub(crate) cps: Vec<[R; 2]>,
    pub(crate) domain: [R; 2],
}

/// The curve's arcs in increasing parameter.
pub(crate) fn arcs(span: &Span) -> Result<Vec<Arc>> {
    Ok(span
        .curve()
        .as_curve3()
        .bezier_arcs()?
        .iter()
        .map(|a| Arc {
            cps: a
                .homogeneous_poles()
                .iter()
                .map(|h| [&h[0] / &h[3], &h[1] / &h[3]])
                .collect(),
            domain: a.domain().clone(),
        })
        .collect())
}

fn binomial(n: usize, k: usize) -> R {
    let mut c = one();
    for i in 0..k {
        c = c * R::from_integer((n - i).into()) / R::from_integer((i + 1).into());
    }
    c
}

/// Bernstein coefficients over `[0, 1]` in ascending powers.
pub(crate) fn power(c: &[R]) -> Vec<R> {
    let n = c.len() - 1;
    let mut out = vec![zero(); n + 1];
    for (i, ci) in c.iter().enumerate() {
        for (k, o) in out.iter_mut().enumerate().skip(i) {
            let term = ci * binomial(n, i) * binomial(n - i, k - i);
            if (k - i) % 2 == 0 {
                *o += term;
            } else {
                *o -= term;
            }
        }
    }
    out
}

pub(crate) fn eval(p: &[R], s: &R) -> R {
    p.iter().rev().fold(zero(), |acc, c| acc * s + c)
}

pub(crate) fn eval_interval(p: &[R], s: &I) -> I {
    p.iter().rev().fold(I::exact(zero()), |acc, c| {
        acc.mul(s).add(&I::exact(c.clone()))
    })
}

/// `a x + b y + k` on an arc, in ascending powers of its parameter.
fn line_poly(cps: &[[R; 2]], line: [&R; 3]) -> Vec<R> {
    let [a, b, k] = line;
    let g: Vec<R> = cps.iter().map(|p| a * &p[0] + b * &p[1] + k).collect();
    power(&g)
}

fn budget() -> Budget {
    Budget::new(RootIsolationOptions::default())
}

/// The distinct roots of `p` in `[0, 1]`, increasing.
pub(crate) fn roots(p: &[R]) -> Result<Vec<AlgebraicRoot>> {
    let mut roots = isolate(
        &IntPolynomial::from_rationals(p),
        zero(),
        one(),
        &mut budget(),
    )?;
    roots.sort_by(|x, y| x.compare_root(y));
    Ok(roots)
}

/// Rational parameters in every gap between consecutive distinct roots in
/// `[0, 1]` (and before the first and after the last), each with its place
/// among the roots: sample `i` lies before root `i`.
fn gaps(roots: &[AlgebraicRoot]) -> Vec<R> {
    let two = R::from_integer(2.into());
    let mut out = Vec::new();
    let mut last = zero();
    for r in roots {
        let (lo, _) = r.isolator();
        if lo > &last {
            out.push((&last + lo) / &two);
        } else {
            // No room before this root (it is at 0, or at the previous
            // isolator's end): a sample that is never taken.
            out.push(last.clone());
        }
        last = r.isolator().1.clone();
    }
    out.push(if last < one() {
        (&last + one()) / &two
    } else {
        one()
    });
    out
}

/// Whether `a x + b y + k` takes the sign `want` (strictly) somewhere on
/// the spline: its sign at each arc's ends and between the arc's roots.
pub(super) fn takes(span: &Span, line: [&R; 3], want: Ordering) -> Result<bool> {
    for arc in arcs(span)? {
        let p = line_poly(&arc.cps, line);
        if p.iter().all(|c| c == &zero()) {
            continue;
        }
        let rs = roots(&p)?;
        let mut samples = gaps(&rs);
        samples.extend([zero(), one()]);
        if samples.iter().any(|s| eval(&p, s).cmp(&zero()) == want) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// A crossing of the line: its rounded curve parameter, the curve's exact
/// point there and that point rounded.
#[derive(Debug, Clone)]
pub(crate) struct Crossing {
    pub(crate) t: f64,
    pub(crate) exact: [R; 2],
    pub(crate) rounded: Point2,
}

/// Where a spline segment meets the line strictly inside it.
#[derive(Debug, Clone, Default)]
pub(crate) struct Meeting {
    /// Crossings along the segment's stored direction.
    pub(crate) crossings: Vec<Crossing>,
    /// Enclosures of the points where it touches the line.
    pub(crate) touches: Vec<[I; 2]>,
    /// The line function's sign on each piece between crossings (one more
    /// than the crossings), along the stored direction.
    pub(crate) sides: Vec<i8>,
}

enum Event {
    Sign(i8),
    Cross(Crossing),
    Touch([I; 2]),
}

/// The binary64 nearest a root's curve parameter `u0 + s (u1 - u0)`.
pub(crate) fn rounded_parameter(root: &mut AlgebraicRoot, domain: &[R; 2]) -> f64 {
    let map = |s: &R| &domain[0] + s * (&domain[1] - &domain[0]);
    for _ in 0..80 {
        let (lo, hi) = root.isolator();
        let (l, h) = (to_f64(&map(lo)), to_f64(&map(hi)));
        if l == h {
            return l;
        }
        root.refine_for_signs(16);
    }
    let (lo, hi) = root.isolator();
    to_f64(&map(&((lo + hi) / R::from_integer(2.into()))))
}

/// The segment's crossings, touches and sides with `a x + b y + d = 0`.
pub(crate) fn meets(span: &Span, line: [&R; 3]) -> Result<Meeting> {
    meets_with(span, &|cps: &[[R; 2]]| line_poly(cps, line))
}

/// The segment's crossings, touches and sides with the zero set of a
/// polynomial condition on its points, given per Bézier arc in the arc's
/// parameter (ascending powers): a line's `a x + b y + d` (degree `p`), a
/// circle's `(x - c_x)^2 + (y - c_y)^2 - r^2` (degree `2 p`).
pub(crate) fn meets_with(span: &Span, condition: &dyn Fn(&[[R; 2]]) -> Vec<R>) -> Result<Meeting> {
    let tangent =
        || Error::Degenerate("a plane tangent to a spline profile segment where it crosses it");
    let curve = span.curve().as_curve3();
    let (first, last) = curve.domain();
    let arcs = arcs(span)?;
    let count = arcs.len();
    let mut events: Vec<Event> = Vec::new();
    for (k, arc) in arcs.iter().enumerate() {
        let p = condition(&arc.cps);
        if p.iter().all(|c| c == &zero()) {
            return Err(Error::OutOfDomain(
                "a spline profile segment along the plane",
            ));
        }
        let (px, py) = (
            power(&arc.cps.iter().map(|c| c[0].clone()).collect::<Vec<_>>()),
            power(&arc.cps.iter().map(|c| c[1].clone()).collect::<Vec<_>>()),
        );
        let mut rs = roots(&p)?;
        let samples = gaps(&rs);
        for (i, root) in rs.iter_mut().enumerate() {
            let sign = eval(&p, &samples[i]).cmp(&zero()) as i8;
            if sign != 0 {
                events.push(Event::Sign(sign));
            }
            let (lo, hi) = (root.isolator().0.clone(), root.isolator().1.clone());
            // The segment's own ends are its path points; an arc's start
            // is the previous arc's end.
            if lo == zero() && hi == zero() {
                continue;
            }
            if lo == one() && hi == one() {
                if k + 1 == count {
                    continue;
                }
                // At an interior knot: a crossing when simple.
                if root.multiplicity() != 1 {
                    return Err(tangent());
                }
                let t = to_f64(&arc.domain[1]);
                let exact = arc.cps[arc.cps.len() - 1].clone();
                let rounded = Point2::new(to_f64(&exact[0]), to_f64(&exact[1]));
                events.push(Event::Cross(Crossing { t, exact, rounded }));
                continue;
            }
            match root.multiplicity() {
                1 => {
                    let t = rounded_parameter(root, &arc.domain);
                    if !(first < t && t < last) {
                        return Err(Error::Degenerate(
                            "a split within the resolution of a profile vertex",
                        ));
                    }
                    let e = curve.exact_point(t)?;
                    let exact = [e[0].clone(), e[1].clone()];
                    let rounded = Point2::new(to_f64(&exact[0]), to_f64(&exact[1]));
                    events.push(Event::Cross(Crossing { t, exact, rounded }));
                }
                m if m % 2 == 0 => {
                    let s = I::new(lo, hi);
                    events.push(Event::Touch([
                        eval_interval(&px, &s),
                        eval_interval(&py, &s),
                    ]));
                }
                _ => return Err(tangent()),
            }
        }
        let sign = eval(&p, &samples[rs.len()]).cmp(&zero()) as i8;
        if sign != 0 {
            events.push(Event::Sign(sign));
        }
    }
    let mut out = Meeting {
        sides: vec![0],
        ..Default::default()
    };
    for e in events {
        match e {
            Event::Sign(s) => {
                let side = out.sides.last_mut().expect("a side");
                if *side == 0 {
                    *side = s;
                } else if *side != s {
                    return Err(Error::ComputationLimit(
                        "a spline's side between its crossings",
                    ));
                }
            }
            Event::Cross(c) => {
                if out.crossings.last().is_some_and(|l| l.t >= c.t) {
                    return Err(Error::Degenerate(
                        "two crossings of a spline profile segment within rounding",
                    ));
                }
                out.crossings.push(c);
                out.sides.push(0);
            }
            Event::Touch(t) => out.touches.push(t),
        }
    }
    if out.sides.contains(&0) {
        return Err(Error::ComputationLimit(
            "a spline's side between its crossings",
        ));
    }
    if span.is_reversed() {
        out.crossings.reverse();
        out.touches.reverse();
        out.sides.reverse();
    }
    Ok(out)
}

/// Inserts `u` into a clamped knot sequence until its multiplicity is the
/// degree, exactly; the index of its first copy.
fn insert(knots: &mut Vec<R>, poles: &mut Vec<[R; 2]>, p: usize, u: &R) -> usize {
    let existing = knots.iter().filter(|k| *k == u).count();
    for _ in existing..p {
        let n = poles.len();
        let k = (p..n)
            .rev()
            .find(|&k| &knots[k] <= u)
            .expect("an interior parameter");
        let mut new = Vec::with_capacity(n + 1);
        for i in 0..=n {
            if i + p <= k {
                new.push(poles[i].clone());
            } else if i > k {
                new.push(poles[i - 1].clone());
            } else {
                let alpha = (u - &knots[i]) / (&knots[i + p] - &knots[i]);
                let beta = one() - &alpha;
                new.push([
                    &alpha * &poles[i][0] + &beta * &poles[i - 1][0],
                    &alpha * &poles[i][1] + &beta * &poles[i - 1][1],
                ]);
            }
        }
        knots.insert(k + 1, u.clone());
        *poles = new;
    }
    knots.iter().position(|k| k == u).expect("an inserted knot")
}

/// The spline's curve restricted to `[t0, t1]` of its parameter, exactly by
/// knot insertion, its poles rounded; traversed in the span's direction.
pub(crate) fn restrict(span: &Span, t0: f64, t1: f64) -> Result<Span> {
    let curve = span.curve().as_curve3();
    let p = curve.degree();
    let (first, last) = curve.domain();
    if (t0, t1) == (first, last) {
        return Ok(span.clone());
    }
    let mut knots: Vec<R> = Vec::new();
    for (k, m) in curve.knots().iter().zip(curve.multiplicities()) {
        knots.extend(std::iter::repeat_n(q(*k), *m));
    }
    let mut poles: Vec<[R; 2]> = curve.poles().iter().map(|c| [q(c.x), q(c.y)]).collect();
    if t1 < last {
        let u = q(t1);
        let j = insert(&mut knots, &mut poles, p, &u);
        knots.truncate(j + p);
        knots.push(u);
        poles.truncate(j);
    }
    if t0 > first {
        let u = q(t0);
        let j = insert(&mut knots, &mut poles, p, &u);
        knots.drain(..j);
        knots.insert(0, u);
        poles.drain(..j - 1);
    }
    let mut values: Vec<f64> = Vec::new();
    let mut mults: Vec<usize> = Vec::new();
    for k in &knots {
        let v = to_f64(k);
        if values.last() == Some(&v) {
            *mults.last_mut().expect("a knot") += 1;
        } else {
            values.push(v);
            mults.push(1);
        }
    }
    let poles: Vec<Point2> = poles
        .iter()
        .map(|c| Point2::new(to_f64(&c[0]), to_f64(&c[1])))
        .collect();
    let piece = SplineSpan::whole(BSplineCurve2::new(p, poles, None, values, mults)?);
    Ok(if span.is_reversed() {
        piece.reversed()
    } else {
        piece
    })
}

/// A piece's direction at its start (`start`) or end, along its traversal:
/// its first or last nonzero control leg.
pub(crate) fn tangent(span: &Span, start: bool) -> (f64, f64) {
    let poles = span.curve().poles();
    let forward = start != span.is_reversed();
    let legs: Vec<(f64, f64)> = poles
        .windows(2)
        .map(|w| (w[1].x - w[0].x, w[1].y - w[0].y))
        .filter(|d| *d != (0.0, 0.0))
        .collect();
    let d = if forward {
        legs[0]
    } else {
        legs[legs.len() - 1]
    };
    if span.is_reversed() {
        (-d.0, -d.1)
    } else {
        d
    }
}

/// Twice the signed area between a piece and its chord, `∮ x dy - y dx`
/// along the piece less the chord's, in binary64 from its Bézier arcs'
/// control points in the power basis (for a cycle's orientation only).
pub(crate) fn twice_area_beyond_chord(span: &Span) -> f64 {
    let Ok(arcs) = arcs(span) else {
        return 0.0;
    };
    let power_f64 = |c: Vec<f64>| -> Vec<f64> {
        let n = c.len() - 1;
        let binomial =
            |n: usize, k: usize| (0..k).fold(1.0, |acc, m| acc * (n - m) as f64 / (m + 1) as f64);
        let mut out = vec![0.0; n + 1];
        for (i, ci) in c.iter().enumerate() {
            for (k, o) in out.iter_mut().enumerate().skip(i) {
                let term = ci * binomial(n, i) * binomial(n - i, k - i);
                *o += if (k - i) % 2 == 0 { term } else { -term };
            }
        }
        out
    };
    // ∫_0^1 a(s) b'(s) ds of two power series.
    let cross_integral = |a: &[f64], b: &[f64]| {
        let mut total = 0.0;
        for (i, ai) in a.iter().enumerate() {
            for (k, bk) in b.iter().enumerate().skip(1) {
                total += ai * bk * k as f64 / (i + k) as f64;
            }
        }
        total
    };
    let mut twice = 0.0;
    for arc in &arcs {
        let x = power_f64(arc.cps.iter().map(|c| to_f64(&c[0])).collect());
        let y = power_f64(arc.cps.iter().map(|c| to_f64(&c[1])).collect());
        twice += cross_integral(&x, &y) - cross_integral(&y, &x);
    }
    let (first, last) = (
        &arcs[0].cps[0],
        &arcs[arcs.len() - 1].cps[arcs[arcs.len() - 1].cps.len() - 1],
    );
    let (a, b) = (
        Point2::new(to_f64(&first[0]), to_f64(&first[1])),
        Point2::new(to_f64(&last[0]), to_f64(&last[1])),
    );
    let chord = a.x * b.y - b.x * a.y;
    if span.is_reversed() {
        -(twice - chord)
    } else {
        twice - chord
    }
}

/// The part of a spline's curve between parameters `u0` and `u1`, traversed
/// from `u0` (the same restriction a section of the curve makes).
pub(crate) fn piece(span: &Span, u0: f64, u1: f64) -> Result<Span> {
    let whole = SplineSpan::whole(span.curve().clone());
    let part = restrict(&whole, u0.min(u1), u0.max(u1))?;
    Ok(if u0 > u1 { part.reversed() } else { part })
}

/// A piece's image on a plane over the profile: each pole `p` at the height
/// `height(p)`, exact for a height affine in `p` up to rounding.
pub(super) fn plane_image(
    piece: &Span,
    frame: crate::Frame3,
    height: &dyn Fn(Point2) -> f64,
) -> Result<crate::topology::Curve3> {
    let c = piece.curve().as_curve3();
    let poles = c
        .poles()
        .iter()
        .map(|p| {
            let p = Point2::new(p.x, p.y);
            frame.point(p, height(p))
        })
        .collect();
    let curve = crate::BSplineCurve3::new(
        c.degree(),
        poles,
        None,
        c.knots().to_vec(),
        c.multiplicities().to_vec(),
    )?;
    let whole = SplineSpan::whole(curve);
    Ok(crate::topology::Curve3::BSpline(if piece.is_reversed() {
        whole.reversed()
    } else {
        whole
    }))
}

/// A piece's pcurve on its wall (`u` the curve's parameter, `v` the height
/// above the wall's low end): `u` itself, as the spline with its Greville
/// abscissae for poles, and `v(p)` at each pole.
pub(super) fn wall_pcurve(
    piece: &Span,
    v: &dyn Fn(Point2) -> f64,
) -> Result<crate::topology::Curve2> {
    let c = piece.curve().as_curve3();
    let p = c.degree();
    let mut knots: Vec<R> = Vec::new();
    for (k, m) in c.knots().iter().zip(c.multiplicities()) {
        knots.extend(std::iter::repeat_n(q(*k), *m));
    }
    let degree = R::from_integer((p as i64).into());
    let poles: Vec<Point2> = c
        .poles()
        .iter()
        .enumerate()
        .map(|(i, pole)| {
            let greville: R = knots[i + 1..=i + p].iter().sum::<R>() / &degree;
            Point2::new(to_f64(&greville), v(Point2::new(pole.x, pole.y)))
        })
        .collect();
    let curve = BSplineCurve2::new(
        p,
        poles,
        None,
        c.knots().to_vec(),
        c.multiplicities().to_vec(),
    )?;
    let whole = SplineSpan::whole(curve);
    Ok(crate::topology::Curve2::BSpline(if piece.is_reversed() {
        whole.reversed()
    } else {
        whole
    }))
}
