//! Rational B-splines against tori and cones (S7c.2 of `REVIEW_NOTES.md`).
//!
//! * A torus: its quartic in homogeneous form,
//!   `(|D|^2 + (R^2 - r^2) W^2)^2 - 4 R^2 W^2 (|D|^2 - (D . a)^2 / |a|^2)` with `D = X - o W`, of degree
//!   `4 p` on each span, through the exact isolator shared with
//!   `spline_sphere` and `spline_cylinder`: every root with its contact
//!   orders, every span on the torus an overlap.
//! * A cone, the set `rho = |r + h tan a|` about its stored axis: on a span
//!   `W^2 F = Q0 + tau Q1 + tau^2 Q2` with exact `Q0 = |D|^2 - (D . a)^2 /
//!   |a|^2 - r^2 W^2`, `Q1 = -2 r (D . a) W`, `Q2 = -(D . a)^2` and one
//!   irrational `tau = tan a / |a|`, transcendental (`a` a nonzero rational,
//!   Lindemann-Weierstrass). So a span lies on the cone exactly when all
//!   three vanish, and a root common to all three (their exact gcd: a
//!   rational apex) is an exact point, a tangency (the apex is singular).
//!   The other roots are isolated by certified subdivision with `tau` in
//!   intervals (binary64 first, rational when that cannot decide); a root
//!   that cannot be certified is `ComputationLimit`.
use super::analytic::Enclosure;
use super::conic_surface::positive_roots;
use super::curve_surface::CurvePoint;
use super::procedural::{bounds, bounds3, limit, q, E, X};
use super::spline_surface::{intersect, SplineSurfaceIntersection, SplineSurfaceOptions};
use super::tangency::{padd, pmul, pscale, P};
use crate::certified::{Fast, Interval as I, Real};
use crate::polynomial::real::{isolate, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::topology::Surface;
use crate::{BSplineCurve3, Error, Result};
use num_rational::BigRational as R;

fn frame_data(surface: &Surface) -> (X, X) {
    let f = match surface {
        Surface::Plane(f)
        | Surface::Cylinder { frame: f, .. }
        | Surface::Cone { frame: f, .. }
        | Surface::Sphere { frame: f, .. }
        | Surface::Torus { frame: f, .. } => f,
        Surface::BSpline(_) => unreachable!("an analytic surface"),
    };
    (f.origin().to_array().map(q), f.normal().to_array().map(q))
}

/// `D = X - o W`, `|D|^2`, `D . a` and `W` of a span's homogeneous
/// polynomials `(X, W)`.
fn relative(h: &[Vec<R>; 4], o: &X, a: &X) -> (P, P, P) {
    let d: [P; 3] = std::array::from_fn(|c| padd(&h[c], &pscale(&h[3], &-&o[c])));
    let dd = (0..3).fold(Vec::new(), |acc, c| padd(&acc, &pmul(&d[c], &d[c])));
    let da = (0..3).fold(Vec::new(), |acc, c| padd(&acc, &pscale(&d[c], &a[c])));
    (dd, da, h[3].clone())
}

/// The whole closed domain of a rational B-spline against a torus.
pub fn spline_torus(curve: &BSplineCurve3, torus: &Surface) -> Result<SplineSurfaceIntersection> {
    let Surface::Torus { major, minor, .. } = torus else {
        return Err(Error::OutOfDomain("spline_torus takes a torus"));
    };
    let (o, a) = frame_data(torus);
    let (big2, k) = (
        q(*major) * q(*major),
        q(*major) * q(*major) - q(*minor) * q(*minor),
    );
    let aa: R = a.iter().map(|x| x * x).sum();
    let (first, last) = curve.domain();
    intersect(
        curve,
        first,
        last,
        SplineSurfaceOptions::default(),
        128,
        move |h| {
            let (dd, da, w) = relative(h, &o, &a);
            let w2 = pmul(&w, &w);
            let s = padd(&dd, &pscale(&w2, &k));
            let m = padd(
                &dd,
                &pscale(&pmul(&da, &da), &-(R::from_integer(1.into()) / &aa)),
            );
            let f = padd(
                &pmul(&s, &s),
                &pscale(&pmul(&w2, &m), &-(R::from_integer(4.into()) * &big2)),
            );
            IntPolynomial::from_rationals(&f)
        },
    )
}

/// A rational B-spline against a cone: isolated points and overlaps.
#[derive(Debug, Clone, PartialEq)]
pub struct SplineConeIntersection {
    /// Sorted by parameter, none inside an overlap.
    pub points: Vec<CurvePoint>,
    /// Maximal runs of whole spans on the cone, between exact knots.
    pub overlaps: Vec<[f64; 2]>,
}

/// The whole closed domain of a rational B-spline against a cone.
pub fn spline_cone(curve: &BSplineCurve3, cone: &Surface) -> Result<SplineConeIntersection> {
    let Surface::Cone {
        radius, half_angle, ..
    } = cone
    else {
        return Err(Error::OutOfDomain("spline_cone takes a cone"));
    };
    let (o, a) = frame_data(cone);
    let r = q(*radius);
    let aa: R = a.iter().map(|x| x * x).sum();
    let (first, last) = curve.domain();
    let spans = curve
        .knot_vector()
        .spans_in(first, last, crate::spline::MAX_POLES)?;
    let mut points: Vec<CurvePoint> = Vec::new();
    let mut overlaps: Vec<[R; 2]> = Vec::new();
    let mut budget = Budget::new(RootIsolationOptions::default());
    for at in spans {
        let length = &at.end - &at.start;
        let lower = (&at.lower - &at.start) / &length;
        let upper = (&at.upper - &at.start) / &length;
        let h = curve.span_polynomial(at.index);
        let (dd, da, w) = relative(&h, &o, &a);
        let q0 = padd(
            &padd(
                &dd,
                &pscale(&pmul(&da, &da), &-(R::from_integer(1.into()) / &aa)),
            ),
            &pscale(&pmul(&w, &w), &-(&r * &r)),
        );
        let q1 = pscale(&pmul(&da, &w), &-(R::from_integer(2.into()) * &r));
        let q2 = pscale(&pmul(&da, &da), &R::from_integer((-1).into()));
        // Integer images only for exact zero tests and the gcd: each is a
        // different positive multiple, so the tau combination stays rational.
        let qp = [q0, q1, q2];
        let qs = qp.clone().map(|p| IntPolynomial::from_rationals(&p));
        if qs.iter().all(IntPolynomial::is_zero) {
            match overlaps.last_mut().filter(|o| o[1] == at.lower) {
                Some(last) => last[1] = at.upper.clone(),
                None => overlaps.push([at.lower.clone(), at.upper.clone()]),
            }
            continue;
        }
        let to_global = |x: &I| {
            x.mul(&I::exact(length.clone()))
                .add(&I::exact(at.start.clone()))
        };
        let point_at = |x: &I| -> Result<E<I>> {
            let hv: [I; 4] = std::array::from_fn(|c| {
                h[c].iter().rev().fold(I::exact_f64(0.0), |acc, k| {
                    acc.mul(x).add(&I::exact(k.clone()))
                })
            });
            let mut p: E<I> = std::array::from_fn(|_| I::exact_f64(0.0));
            for c in 0..3 {
                p[c] = hv[c].div(&hv[3]).ok_or(limit("a spline's point"))?;
            }
            Ok(p)
        };
        // Exact common roots (a rational apex): tau being transcendental, a
        // root's order in F is its least order in the Q's, its order in
        // their gcd.
        let common = qs[0].gcd(&qs[1]).gcd(&qs[2]);
        if !common.is_constant() {
            for mut root in isolate(&common, lower.clone(), upper.clone(), &mut budget)? {
                let tangent = root.multiplicity() > 1;
                root.refine_for_signs(128);
                let (lo, hi) = root.isolator();
                let x = I::new(lo.clone(), hi.clone());
                points.push(CurvePoint {
                    parameter: bounds(&to_global(&x)),
                    point: bounds3(&point_at(&x)?),
                    tangent,
                });
            }
        }
        // The rest: F / common in intervals of tau.
        let common: P = common.0.iter().cloned().map(R::from_integer).collect();
        let rest: Vec<P> = qp.iter().map(|p| divide(p, &common)).collect();
        let found = match span_roots::<Fast>(&rest, *half_angle, &aa, &lower, &upper) {
            Err(Error::ComputationLimit(_)) => {
                span_roots::<I>(&rest, *half_angle, &aa, &lower, &upper)?
            }
            other => other?,
        };
        for x in found {
            let x = I::new(q(x[0]), q(x[1]));
            points.push(CurvePoint {
                parameter: bounds(&to_global(&x)),
                point: bounds3(&point_at(&x)?),
                tangent: false,
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
    // Exact apexes at a knot come from both spans; points at an overlap's
    // ends belong to it.
    points.sort_by(|a, b| a.parameter[0].total_cmp(&b.parameter[0]));
    points.dedup_by(|b, a| a.parameter == b.parameter && a.tangent && b.tangent);
    points.retain(|p| {
        !overlaps
            .iter()
            .any(|[lo, hi]| *lo <= p.parameter[0] && p.parameter[1] <= *hi)
    });
    Ok(SplineConeIntersection { points, overlaps })
}

/// The roots in `[lower, upper]` (the span's local parameter) of
/// `Q0 + tau Q1 + tau^2 Q2` in one tier.
fn span_roots<T: Real>(
    qs: &[P],
    half_angle: f64,
    aa: &R,
    lower: &R,
    upper: &R,
) -> Result<Vec<Enclosure>> {
    let (c, s) = T::cos_sin(&T::exact_f64(half_angle));
    let tau = s
        .div(&c)
        .and_then(|t| t.div(&T::from_r(aa).sqrt()))
        .ok_or(limit("a cone's tangent"))?;
    let n = qs.iter().map(Vec::len).max().unwrap_or(0);
    let mut coeffs: Vec<T> = (0..n).map(|_| T::exact_f64(0.0)).collect();
    let mut power = T::exact_f64(1.0);
    for p in qs {
        for (k, c) in p.iter().enumerate() {
            coeffs[k] = coeffs[k].add(&power.mul(&T::from_r(c)));
        }
        power = power.mul(&tau);
    }
    // Substitute x = lower + (upper - lower) z / (1 + z), z in (0, oo):
    // positive roots of (1 + z)^n F(x(z)) are the span's interior roots.
    let (lo, hi) = (T::from_r(lower), T::from_r(upper));
    let width = hi.sub(&lo);
    // F(lower + width y), y in [0, 1], then y = z / (1 + z).
    let mut shifted = coeffs.clone();
    taylor_shift(&mut shifted, &lo);
    let mut w = T::exact_f64(1.0);
    for c in shifted.iter_mut() {
        *c = c.mul(&w);
        w = w.mul(&width);
    }
    // (1 + z)^(n-1) sum_k g_k (z / (1 + z))^k = sum_k g_k z^k (1 + z)^(n-1-k).
    let deg = shifted.len().saturating_sub(1);
    let mut out: Vec<T> = (0..=deg).map(|_| T::exact_f64(0.0)).collect();
    for (k, g) in shifted.iter().enumerate() {
        // z^k (1 + z)^(deg - k): binomial coefficients.
        let mut binom = 1.0f64;
        for j in 0..=deg - k {
            out[k + j] = out[k + j].add(&g.mul(&T::exact_f64(binom)));
            binom = binom * ((deg - k - j) as f64) / ((j + 1) as f64);
        }
    }
    if deg == 0 {
        return Ok(Vec::new());
    }
    let roots = positive_roots(&out)?;
    // Back to the span's local parameter: y = z / (1 + z).
    Ok(roots
        .into_iter()
        .map(|[a, b]| {
            let y = |z: f64| {
                let v = T::exact_f64(z);
                v.div(&T::exact_f64(1.0).add(&v)).expect("1 + z > 0")
            };
            let (ya, yb) = (y(a), y(b));
            let x = |v: &T| lo.add(&width.mul(v));
            [x(&ya).bounds_f64().0, x(&yb).bounds_f64().1]
        })
        .collect())
}

/// The exact quotient `p / d` (`d` divides `p`; the zero polynomial stays).
fn divide(p: &P, d: &P) -> P {
    let mut rem: P = p.clone();
    while rem.last().is_some_and(|x| *x == R::from_integer(0.into())) {
        rem.pop();
    }
    if rem.is_empty() {
        return rem;
    }
    let mut quotient = vec![R::from_integer(0.into()); rem.len() + 1 - d.len()];
    let lead = d.last().expect("a nonzero divisor");
    while rem.len() >= d.len() {
        let shift = rem.len() - d.len();
        let c = rem.last().unwrap() / lead;
        for (i, b) in d.iter().enumerate() {
            rem[i + shift] -= &c * b;
        }
        quotient[shift] = c;
        rem.pop();
        while rem.last().is_some_and(|x| *x == R::from_integer(0.into())) {
            rem.pop();
        }
    }
    debug_assert!(rem.is_empty());
    quotient
}

/// `c(x) -> c(x + s)` in place.
fn taylor_shift<T: Real>(c: &mut [T], s: &T) {
    let n = c.len();
    for i in 0..n {
        for j in (i..n - 1).rev() {
            let t = c[j + 1].mul(s);
            c[j] = c[j].add(&t);
        }
    }
}
