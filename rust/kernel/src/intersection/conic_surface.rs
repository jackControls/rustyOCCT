//! Ellipses and hyperbolas against analytic surfaces (S7c.2 of
//! `REVIEW_NOTES.md`).
//!
//! A conic has stored binary64 data, as OCCT's `Geom_Ellipse` and
//! `Geom_Hyperbola`: a `Frame3` and semi-axes `major`, `minor`, the exact
//! point sets `o + major cos t x + minor sin t y` and `o + major cosh t x +
//! minor sinh t y` in the stored axes (one branch of the hyperbola).
//!
//! * Against a plane, cylinder, sphere or torus: a point `o + xi x + eta y`
//!   is on the curve exactly when `minor^2 xi^2 +- major^2 eta^2 = major^2
//!   minor^2` (and `xi > 0` on the hyperbola's branch), so the circle's
//!   resultant method applies (`curve_surface::section`): exact
//!   multiplicities (tangencies) and containment. The parameter is a
//!   certified `atan2(eta / minor, xi / major)` or `asinh(eta / minor)`.
//! * Against a cone: an ellipse is the circle's trigonometric polynomial of
//!   degree two; a hyperbola's function times `4 z^2` with `z = e^t` is a
//!   quartic in `z > 0` with interval coefficients, its roots isolated by
//!   certified subdivision between Cauchy bounds of it and its reciprocal,
//!   `t = ln z`. A tangency there is `ComputationLimit`.
use super::analytic::Enclosure;
use super::curve_surface::{
    section, trig_cone, CurvePoint, CurveSurfaceIntersection, Implicit, Section,
};
use super::procedural::{bounds, bounds3, e3, edot, esub, limit, pi, q, E, X};
use super::tangency::{badd, bmul, bscale, B};
use super::torus_curves::angle;
use crate::certified::{Fast, Interval as I, Real};
use crate::topology::Surface;
use crate::{Error, Frame3, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// An ellipse or one branch of a hyperbola with stored binary64 data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Conic {
    /// `o + major cos t x + minor sin t y`, `major >= minor > 0`, `t` in
    /// `(-pi, pi]`.
    Ellipse {
        frame: Frame3,
        major: f64,
        minor: f64,
    },
    /// `o + major cosh t x + minor sinh t y`, both positive, `t` real.
    Hyperbola {
        frame: Frame3,
        major: f64,
        minor: f64,
    },
}

/// The intersection of an ellipse or a hyperbola's branch with an analytic
/// surface: `Empty`, `Contained`, or points sorted by the conic's parameter.
pub fn conic_surface(conic: &Conic, surface: &Surface) -> Result<CurveSurfaceIntersection> {
    let (frame, major, minor, hyperbola) = match *conic {
        Conic::Ellipse {
            frame,
            major,
            minor,
        } => (frame, major, minor, false),
        Conic::Hyperbola {
            frame,
            major,
            minor,
        } => (frame, major, minor, true),
    };
    if !(major.is_finite() && minor.is_finite()) {
        return Err(Error::NonFinite("a conic's semi-axes"));
    }
    if !(minor > 0.0 && major > 0.0) || (!hyperbola && major < minor) {
        return Err(Error::OutOfDomain(
            "a conic's semi-axes (positive, an ellipse's major the larger)",
        ));
    }
    let s = Implicit::of(surface)?;
    let o = frame.origin().to_array().map(q);
    let (x, y) = (frame.x().to_array().map(q), frame.y().to_array().map(q));
    let (a, b) = (q(major), q(minor));
    let mut points = if s.cone.is_some() {
        if hyperbola {
            hyperbola_cone(&o, &x, &y, &a, &b, &s)?
        } else {
            let ax = x.clone().map(|v| v * &a);
            let ay = y.clone().map(|v| v * &b);
            trig_cone(&o, &ax, &ay, &s)?
                .into_iter()
                .map(|(parameter, point)| CurvePoint {
                    parameter,
                    point,
                    tangent: false,
                })
                .collect()
        }
    } else {
        // minor^2 xi^2 +- major^2 eta^2 - major^2 minor^2.
        let (a2, b2) = (&a * &a, &b * &b);
        let sign = if hyperbola { -1 } else { 1 };
        let conic = |xi: &B, eta: &B| {
            badd(
                &badd(
                    &bscale(&bmul(xi, xi), &b2),
                    &bscale(&bmul(eta, eta), &(&a2 * R::from_integer(sign.into()))),
                ),
                &vec![vec![-(&a2 * &b2)]],
            )
        };
        match section(&o, &x, &y, &conic, &|p| s.at(p))? {
            Section::Contained => return Ok(CurveSurfaceIntersection::Contained),
            Section::Points(found) => {
                let mut out = Vec::new();
                for p in found {
                    let u = p.xi.div(&I::from_r(&a)).ok_or(limit("a conic's point"))?;
                    let v = p.eta.div(&I::from_r(&b)).ok_or(limit("a conic's point"))?;
                    let t = if hyperbola {
                        // The other branch (xi < 0) is not this curve.
                        match u.sign() {
                            Some(Ordering::Greater) => {}
                            Some(Ordering::Less) => continue,
                            _ => return Err(limit("a hyperbola's branch")),
                        }
                        crate::certified::asinh(&v)
                    } else {
                        turn_angle(&v, &u)?
                    };
                    out.push(CurvePoint {
                        parameter: bounds(&t),
                        point: bounds3(&p.point),
                        tangent: p.tangent,
                    });
                }
                out
            }
        }
    };
    points.sort_by(|a, b| {
        (a.parameter[0] + a.parameter[1]).total_cmp(&(b.parameter[0] + b.parameter[1]))
    });
    Ok(if points.is_empty() {
        CurveSurfaceIntersection::Empty
    } else {
        CurveSurfaceIntersection::Points(points)
    })
}

/// `atan2(y, x)` in `(-pi, pi]` unless it straddles `pi`.
pub(super) fn turn_angle(y: &I, x: &I) -> Result<I> {
    let t = angle(y, x)?;
    Ok(
        if y.sign() == Some(Ordering::Less) && x.sign() == Some(Ordering::Less) {
            t.sub(&pi::<I>().mul(&I::exact_f64(2.0)))
        } else {
            t
        },
    )
}

/// A hyperbola's branch against a cone: with `z = e^t`,
/// `z (p - V) = c0 + c1 z + c2 z^2` (`c0 = (A x - B y) / 2`, `c1 = o - V`,
/// `c2 = (A x + B y) / 2`), so `z^2 F` is the quartic
/// `sum_{i, j} Phi(ci, cj) z^(i + j)` with the cone's bilinear form `Phi(u,
/// w) = cos^2 h u . w - (u . a)(w . a)`. Binary64 intervals first, rational
/// ones when those cannot decide.
fn hyperbola_cone(o: &X, x: &X, y: &X, a: &R, b: &R, s: &Implicit) -> Result<Vec<CurvePoint>> {
    let found = match quartic_roots::<Fast>(o, x, y, a, b, s) {
        Err(Error::ComputationLimit(_)) => quartic_roots::<I>(o, x, y, a, b, s)?,
        other => other?,
    };
    let mut out = Vec::new();
    for z in found {
        let zz = I::new(q(z[0]), q(z[1]));
        let t = crate::certified::ln(&zz).ok_or(limit("a hyperbola's parameter"))?;
        let inv = I::exact_f64(1.0)
            .div(&zz)
            .ok_or(limit("a hyperbola's point"))?;
        let (ch, sh) = (
            zz.add(&inv).mul(&I::exact_f64(0.5)),
            zz.sub(&inv).mul(&I::exact_f64(0.5)),
        );
        let (ax, by) = (I::from_r(a), I::from_r(b));
        let p: E<I> = std::array::from_fn(|i| {
            I::from_r(&o[i])
                .add(&ax.mul(&ch).mul(&I::from_r(&x[i])))
                .add(&by.mul(&sh).mul(&I::from_r(&y[i])))
        });
        out.push(CurvePoint {
            parameter: bounds(&t),
            point: bounds3(&p),
            tangent: false,
        });
    }
    Ok(out)
}

/// The positive roots of the hyperbola's quartic in `z`, in one tier.
fn quartic_roots<T: Real>(
    o: &X,
    x: &X,
    y: &X,
    a: &R,
    b: &R,
    s: &Implicit,
) -> Result<Vec<Enclosure>> {
    let (v, axis, c2) = s.cone::<T>()?;
    let half = R::new(1.into(), 2.into());
    let ax: X = x.clone().map(|c| c * a * &half);
    let by: X = y.clone().map(|c| c * b * &half);
    let c: [E<T>; 3] = [
        std::array::from_fn(|i| T::from_r(&(&ax[i] - &by[i]))),
        esub(&e3::<T>(o), &v),
        std::array::from_fn(|i| T::from_r(&(&ax[i] + &by[i]))),
    ];
    let phi = |u: &E<T>, w: &E<T>| {
        edot(u, w)
            .mul(&c2)
            .sub(&edot(u, &axis).mul(&edot(w, &axis)))
    };
    let mut coeffs: Vec<T> = (0..5).map(|_| T::exact_f64(0.0)).collect();
    for i in 0..3 {
        for j in 0..3 {
            coeffs[i + j] = coeffs[i + j].add(&phi(&c[i], &c[j]));
        }
    }
    positive_roots(&coeffs)
}

/// The positive real roots of `sum coeffs[k] z^k` (interval coefficients,
/// both ends certainly nonzero): isolated between the Cauchy bounds of the
/// polynomial and its reciprocal by certified subdivision (the mean-value
/// form excluding a box; a root certified by a sign change and a derivative
/// of one sign, narrowed by bisection on certain signs).
pub(super) fn positive_roots<T: Real>(coeffs: &[T]) -> Result<Vec<Enclosure>> {
    let n = coeffs.len() - 1;
    let certain = |x: &T| matches!(x.sign(), Some(Ordering::Less | Ordering::Greater));
    let (lead, last) = (&coeffs[n], &coeffs[0]);
    if !certain(lead) || !certain(last) {
        return Err(limit("a polynomial's end coefficients"));
    }
    // 1 + max |c_k / c_n| bounds every root; its reciprocal's bound gives a
    // lower bound of the positive roots.
    let ratio = |k: usize, d: &T| -> Result<f64> {
        let r = coeffs[k].div(d).ok_or(limit("a Cauchy bound"))?;
        let (lo, hi) = r.bounds_f64();
        Ok(lo.abs().max(hi.abs()))
    };
    let upper = 1.0
        + (0..n)
            .map(|k| ratio(k, lead))
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .fold(0.0, f64::max);
    let lower = 1.0
        / (1.0
            + (1..=n)
                .map(|k| ratio(k, last))
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(0.0, f64::max));
    // Rounded outward.
    let (lower, upper) = (lower * (1.0 - 1e-12), upper * (1.0 + 1e-12));
    let derivative: Vec<T> = (1..=n)
        .map(|k| coeffs[k].mul(&T::exact_f64(k as f64)))
        .collect();
    let eval = |c: &[T], z: &T| {
        c.iter()
            .rev()
            .fold(T::exact_f64(0.0), |acc, k| acc.mul(z).add(k))
    };
    let sign = |z: f64| {
        eval(coeffs, &T::exact_f64(z))
            .sign()
            .filter(|o| *o != Ordering::Equal)
    };
    let mut out = Vec::new();
    let mut pending = vec![(lower, upper)];
    let mut budget = 20_000usize;
    while let Some((lo, hi)) = pending.pop() {
        budget = budget
            .checked_sub(1)
            .ok_or(limit("a polynomial's positive roots"))?;
        let zb = super::procedural::span::<T>(lo, hi);
        if certain(&eval(coeffs, &zb)) {
            continue;
        }
        let m = 0.5 * lo + 0.5 * hi;
        let db = eval(&derivative, &zb);
        let mv = eval(coeffs, &T::exact_f64(m)).add(&db.mul(&zb.sub(&T::exact_f64(m))));
        if certain(&mv) {
            continue;
        }
        let (sa, sb) = (sign(lo), sign(hi));
        if sa.is_some() && sb.is_some() && sa != sb && certain(&db) {
            let (mut a, mut b) = (lo, hi);
            for _ in 0..80 {
                let mid = 0.5 * a + 0.5 * b;
                if !(a < mid && mid < b) {
                    break;
                }
                match sign(mid) {
                    Some(x) if Some(x) == sa => a = mid,
                    Some(_) => b = mid,
                    None => break,
                }
            }
            out.push([a, b]);
            continue;
        }
        // Split off the centre (a root at a round value never on an end).
        let split = lo + (hi - lo) * 0.4453125;
        if !(lo < split && split < hi) {
            return Err(limit("a tangency in intervals"));
        }
        pending.push((split, hi));
        pending.push((lo, split));
    }
    out.sort_by(|a, b| a[0].total_cmp(&b[0]));
    Ok(out)
}
