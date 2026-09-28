//! Lines and circles against analytic surfaces (S7c.1 of `REVIEW_NOTES.md`).
//!
//! A line is the exact set through its two rational points, parameterised
//! `p0 + s (p1 - p0)`; a circle the exact point set of its stored frame (the
//! plane through the origin normal to the stored normal) and radius,
//! parameterised by the angle from the stored x axis. A surface is the exact
//! point set of S7a.
//!
//! * A line against a plane, cylinder, sphere or torus: the surface's
//!   function along it is a polynomial of degree one, two or four with
//!   rational coefficients, its real roots isolated exactly, a multiple root
//!   a tangency, the zero polynomial containment.
//! * A circle against them: in a rational basis `(U, V)` of its plane
//!   (`V` sheared until distinct points have distinct `u`) the surface's
//!   function reduced modulo the circle's conic `a2 v^2 + a1 v + a0` is
//!   `c0(u) + c1(u) v`; the points are the real roots of the resultant
//!   `a2 c0^2 - a1 c0 c1 + a0 c1^2` with `v = -c0 / c1`, a multiple root a
//!   tangency, `c0 = c1 = 0` containment.
//! * A cone's `cos` and `sin` are transcendental: a line's quadratic and a
//!   circle's trigonometric polynomial of degree two in certified intervals,
//!   a root certified by a sign change and a derivative of one sign, a
//!   tangency (undecidable) `ComputationLimit`, except a line through a
//!   rational apex (the apex, a double root).
use super::analytic::{Enclosure, Enclosure3};
use super::procedural::{
    bounds, bounds3, cross, dot, e3, edot, esub, limit, pi, q, span, sub, unit, zero, E, X,
};
use super::tangency::{
    badd, bdot, bl, bmul, bscale, cauchy, enclose, padd, pmul, pscale, reduce, rz, B, P,
};
use super::torus_curves::angle;
use crate::certified::{Fast, Interval as I, Real};
use crate::polynomial::real::{isolate, Budget, IntPolynomial};
use crate::polynomial::RootIsolationOptions;
use crate::topology::{Curve3, Surface};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// A point where a curve meets a surface.
#[derive(Debug, Clone, PartialEq)]
pub struct CurvePoint {
    /// The curve's parameter: a line's `s`, a circle's angle.
    pub parameter: Enclosure,
    pub point: Enclosure3,
    /// A tangency (a multiple root), else a crossing.
    pub tangent: bool,
}

/// The intersection of a curve and a surface.
#[derive(Debug, Clone, PartialEq)]
pub enum CurveSurfaceIntersection {
    Empty,
    /// The curve lies on the surface.
    Contained,
    /// Points sorted by parameter.
    Points(Vec<CurvePoint>),
}

/// The intersection of a line or circle (a `Curve3` line segment's whole
/// line, an arc's whole circle) with an analytic surface.
pub fn curve_surface(curve: &Curve3, surface: &Surface) -> Result<CurveSurfaceIntersection> {
    let s = Implicit::of(surface)?;
    let mut out = match curve {
        Curve3::LineSegment { start, end } => {
            let p0 = start.to_array().map(q);
            let d = sub(&end.to_array().map(q), &p0);
            if zero(&dot(&d, &d)) {
                return Err(Error::Degenerate("a line through one point"));
            }
            if s.cone.is_some() {
                line_cone(&p0, &d, &s)?
            } else {
                line_exact(&p0, &d, &s)?
            }
        }
        Curve3::Circle { frame, radius } | Curve3::CircularArc { frame, radius, .. } => {
            let c = Circle {
                o: frame.origin().to_array().map(q),
                n: frame.normal().to_array().map(q),
                x: frame.x().to_array().map(q),
                y: frame.y().to_array().map(q),
                r: q(*radius),
            };
            if s.cone.is_some() {
                circle_cone(&c, &s)?
            } else {
                circle_exact(&c, &s)?
            }
        }
        Curve3::BSpline(_) => {
            return Err(Error::OutOfDomain(
                "a spline edge against an analytic surface (spline_plane and exact_spline)",
            ))
        }
    };
    if let CurveSurfaceIntersection::Points(points) = &mut out {
        points.sort_by(|a, b| {
            (a.parameter[0] + a.parameter[1]).total_cmp(&(b.parameter[0] + b.parameter[1]))
        });
        if points.is_empty() {
            return Ok(CurveSurfaceIntersection::Empty);
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------ surfaces

/// A surface's implicit function: `n . (p - o)` (a plane), `|p - o|^2 -
/// ((p - o) . a)^2 / |a|^2 - r^2` (a cylinder), `|p - o|^2 - r^2` (a
/// sphere), `(|w|^2 + R^2 - r^2)^2 - 4 R^2 (|w|^2 - (w . a)^2 / |a|^2)` (a
/// torus), all rational; a cone's `cos^2 h |p - V|^2 - ((p - V) . a)^2` in
/// intervals.
struct Implicit {
    kind: Kind,
    o: X,
    a: X,
    r: R,
    minor: R,
    /// A cone's half-angle and radius at its origin.
    cone: Option<(f64, R)>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Plane,
    Cylinder,
    Sphere,
    Torus,
    Cone,
}

impl Implicit {
    fn of(s: &Surface) -> Result<Self> {
        let (kind, f, r, minor, cone) = match s {
            Surface::Plane(f) => (Kind::Plane, f, 0.0, 0.0, None),
            Surface::Cylinder { frame, radius } => (Kind::Cylinder, frame, *radius, 0.0, None),
            Surface::Sphere { frame, radius } => (Kind::Sphere, frame, *radius, 0.0, None),
            Surface::Torus {
                frame,
                major,
                minor,
            } => (Kind::Torus, frame, *major, *minor, None),
            Surface::Cone {
                frame,
                radius,
                half_angle,
            } => (
                Kind::Cone,
                frame,
                *radius,
                0.0,
                Some((*half_angle, q(*radius))),
            ),
            Surface::BSpline(_) => {
                return Err(Error::OutOfDomain("a curve against a spline surface"))
            }
        };
        Ok(Self {
            kind,
            o: f.origin().to_array().map(q),
            a: f.normal().to_array().map(q),
            r: q(r),
            minor: q(minor),
            cone,
        })
    }
    /// The function at a point given by polynomials (in `(u, v)`: `B`).
    fn at(&self, p: &[B; 3]) -> B {
        let w: [B; 3] = std::array::from_fn(|i| badd(&p[i], &vec![vec![-self.o[i].clone()]]));
        let a: [B; 3] = std::array::from_fn(|i| vec![vec![self.a[i].clone()]]);
        let aa = dot(&self.a, &self.a);
        let ww = bdot(&w, &w);
        let wa = bdot(&w, &a);
        let r2 = &self.r * &self.r;
        match self.kind {
            Kind::Plane => wa,
            Kind::Sphere => badd(&ww, &vec![vec![-r2]]),
            Kind::Cylinder => badd(
                &badd(
                    &ww,
                    &bscale(&bmul(&wa, &wa), &(-(R::from_integer(1.into()) / &aa))),
                ),
                &vec![vec![-r2]],
            ),
            Kind::Torus => {
                let k = &r2 - &self.minor * &self.minor;
                let s = badd(&ww, &vec![vec![k]]);
                let m = badd(
                    &ww,
                    &bscale(&bmul(&wa, &wa), &(-(R::from_integer(1.into()) / &aa))),
                );
                badd(
                    &bmul(&s, &s),
                    &bscale(&m, &(-(R::from_integer(4.into()) * &r2))),
                )
            }
            Kind::Cone => unreachable!("a cone is in intervals"),
        }
    }
    /// A cone's apex, axis and `cos^2` in a tier.
    fn cone<T: Real>(&self) -> Result<(E<T>, E<T>, T)> {
        let (h, r) = self.cone.clone().expect("a cone");
        let (c, s) = T::cos_sin(&T::exact_f64(h));
        let a = unit::<T>(&self.a)?;
        let v = if zero(&r) {
            e3(&self.o)
        } else {
            let back = T::from_r(&r)
                .mul(&c)
                .div(&s)
                .ok_or(limit("a cone's apex"))?;
            esub(
                &e3(&self.o),
                &[a[0].mul(&back), a[1].mul(&back), a[2].mul(&back)],
            )
        };
        Ok((v, a, c.square()))
    }
}

/// A root's enclosure as an interval, narrowed.
fn root_interval(root: &mut crate::polynomial::real::AlgebraicRoot) -> I {
    root.refine_for_signs(300);
    let (lo, hi) = root.isolator();
    I::new(lo.clone(), hi.clone())
}

// ------------------------------------------------------------------ lines

/// A line against a plane, cylinder, sphere or torus: the exact polynomial
/// along it.
fn line_exact(p0: &X, d: &X, s: &Implicit) -> Result<CurveSurfaceIntersection> {
    // p(s) = p0 + s d as polynomials in s (the v-degree 0 slot of B).
    let p: [B; 3] = std::array::from_fn(|i| vec![vec![p0[i].clone(), d[i].clone()]]);
    let f = s.at(&p);
    let poly: P = f.first().cloned().unwrap_or_default();
    let int = IntPolynomial::from_rationals(&poly);
    if int.is_zero() {
        return Ok(CurveSurfaceIntersection::Contained);
    }
    let bound = cauchy(&poly);
    let mut budget = Budget::new(RootIsolationOptions::default());
    let roots = isolate(&int, -bound.clone(), bound, &mut budget)?;
    let mut out = Vec::new();
    for mut root in roots {
        let tangent = root.multiplicity() > 1;
        let sv = root_interval(&mut root);
        let point: E<I> =
            std::array::from_fn(|i| I::exact(p0[i].clone()).add(&sv.mul(&I::exact(d[i].clone()))));
        out.push(CurvePoint {
            parameter: bounds(&sv),
            point: bounds3(&point),
            tangent,
        });
    }
    Ok(CurveSurfaceIntersection::Points(out))
}

/// A line against a cone (both nappes): `A s^2 + 2 B s + C` in intervals,
/// binary64 first.
fn line_cone(p0: &X, d: &X, s: &Implicit) -> Result<CurveSurfaceIntersection> {
    // A line through a rational apex meets the cone there only (a double
    // root): its generatrices are never rational.
    if let Some((_, r)) = &s.cone {
        if zero(r) && cross(&sub(&s.o, p0), d).iter().all(zero) {
            let sv = I::exact(dot(&sub(&s.o, p0), d) / dot(d, d));
            return Ok(CurveSurfaceIntersection::Points(vec![CurvePoint {
                parameter: bounds(&sv),
                point: bounds3(&e3::<I>(&s.o)),
                tangent: true,
            }]));
        }
    }
    match line_cone_in::<Fast>(p0, d, s) {
        Err(Error::ComputationLimit(_)) => line_cone_in::<I>(p0, d, s),
        other => other,
    }
}

fn line_cone_in<T: Real>(p0: &X, d: &X, s: &Implicit) -> Result<CurveSurfaceIntersection> {
    let (v, a, c2) = s.cone::<T>()?;
    let form = |x: &E<T>, y: &E<T>| edot(x, y).mul(&c2).sub(&edot(x, &a).mul(&edot(y, &a)));
    let dd = e3::<T>(d);
    let w = esub(&e3(p0), &v);
    let (qa, qb, qc) = (form(&dd, &dd), form(&w, &dd), form(&w, &w));
    let disc = qb.square().sub(&qa.mul(&qc));
    match disc.sign() {
        Some(Ordering::Less) => return Ok(CurveSurfaceIntersection::Empty),
        Some(Ordering::Greater) => {}
        _ => return Err(limit("a line tangent to a cone")),
    }
    let root = disc.sqrt();
    let mut out = Vec::new();
    for sgn in [-1.0, 1.0] {
        let sv = qb
            .neg()
            .add(&root.mul(&T::exact_f64(sgn)))
            .div(&qa)
            .ok_or(limit("a line along a cone's generatrix"))?;
        let point: E<T> =
            std::array::from_fn(|i| T::from_r(&p0[i]).add(&sv.mul(&T::from_r(&d[i]))));
        out.push(CurvePoint {
            parameter: bounds(&sv),
            point: bounds3(&point),
            tangent: false,
        });
    }
    Ok(CurveSurfaceIntersection::Points(out))
}

// ------------------------------------------------------------------ circles

/// A circle's exact data: origin, stored normal and axes, radius.
struct Circle {
    o: X,
    n: X,
    x: X,
    y: X,
    r: R,
}

impl Circle {
    /// The angle of a point of the circle in its stored frame, in
    /// `(-pi, pi]` unless it straddles `pi`.
    fn angle(&self, p: &E<I>) -> Result<Enclosure> {
        let w = esub(p, &e3(&self.o));
        let (y, x) = (edot(&w, &e3(&self.y)), edot(&w, &e3(&self.x)));
        let mut t = angle(&y, &x)?;
        if y.sign() == Some(Ordering::Less) && x.sign() == Some(Ordering::Less) {
            t = t.sub(&pi::<I>().mul(&I::exact_f64(2.0)));
        }
        Ok(bounds(&t))
    }
}

/// A circle against a plane, cylinder, sphere or torus, exactly.
fn circle_exact(c: &Circle, s: &Implicit) -> Result<CurveSurfaceIntersection> {
    let u0 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        .iter()
        .map(|e| cross(&c.n, &e.map(|k| R::from_integer(k.into()))))
        .find(|u| !u.iter().all(zero))
        .expect("a nonzero normal");
    let v0 = cross(&c.n, &u0);
    let two = R::from_integer(2.into());
    for shear in [0i64, 1, -2, 3, -5, 7, 11, -13] {
        let sh = R::new(shear.into(), 7.into());
        let u = u0.clone();
        let v: X = std::array::from_fn(|i| &v0[i] + &u0[i] * &sh);
        // E0 = a2 v^2 + a1 v + a0.
        let a0: P = vec![-(&c.r * &c.r), rz(), dot(&u, &u)];
        let a1: P = vec![rz(), &two * dot(&u, &v)];
        let a2 = dot(&v, &v);
        let p: [B; 3] = std::array::from_fn(|i| bl(c.o[i].clone(), u[i].clone(), v[i].clone()));
        let (c0, c1) = reduce(&s.at(&p), &a0, &a1, &a2);
        let (i0, i1) = (
            IntPolynomial::from_rationals(&c0),
            IntPolynomial::from_rationals(&c1),
        );
        if i0.is_zero() && i1.is_zero() {
            return Ok(CurveSurfaceIntersection::Contained);
        }
        let res = padd(
            &padd(
                &pscale(&pmul(&c0, &c0), &a2),
                &pscale(&pmul(&pmul(&a1, &c0), &c1), &R::from_integer((-1).into())),
            ),
            &pmul(&a0, &pmul(&c1, &c1)),
        );
        let resultant = IntPolynomial::from_rationals(&res);
        if resultant.is_zero() {
            return Err(limit("a circle's resultant"));
        }
        let bound = cauchy(&res);
        let mut budget = Budget::new(RootIsolationOptions::default());
        let roots = isolate(&resultant, -bound.clone(), bound, &mut budget)?;
        if roots
            .iter()
            .any(|x| i1.is_zero() || x.vanishes_polynomial(&i1))
        {
            continue;
        }
        let mut out = Vec::new();
        for mut root in roots {
            // Each real root carries one point (`c1 != 0` there), so its
            // multiplicity is that point's intersection multiplicity.
            let tangent = root.multiplicity() > 1;
            let uu = root_interval(&mut root);
            let vv = enclose(&c0, &uu)
                .div(&enclose(&c1, &uu))
                .ok_or(limit("a circle's point"))?
                .neg();
            let point: E<I> = std::array::from_fn(|i| {
                I::exact(c.o[i].clone())
                    .add(&uu.mul(&I::exact(u[i].clone())))
                    .add(&vv.mul(&I::exact(v[i].clone())))
            });
            out.push(CurvePoint {
                parameter: c.angle(&point)?,
                point: bounds3(&point),
                tangent,
            });
        }
        return Ok(CurveSurfaceIntersection::Points(out));
    }
    Err(limit("a circle's points"))
}

/// A circle against a cone: `G(theta) = f(o + r (cos theta x + sin theta y))`
/// of degree two in `cos theta`, `sin theta`. Its roots are isolated by
/// certified subdivision over one turn cut where `G` is certainly nonzero
/// (the mean-value form excluding a box), each certified by a sign change
/// and a derivative of one sign and narrowed by bisection on certain signs.
fn circle_cone(c: &Circle, s: &Implicit) -> Result<CurveSurfaceIntersection> {
    let fast = Ring::<Fast>::new(c, s)?;
    let cut = match cut(&fast) {
        Err(Error::ComputationLimit(_)) => cut(&Ring::<I>::new(c, s)?)?,
        other => other?,
    };
    // Binary64 intervals first; the boxes they leave undecided (merged)
    // again in exact arithmetic, under a small budget.
    let mut budget = 20_000;
    let (mut found, left) = search(&fast, vec![cut], &mut budget, 1e-10)?;
    let mut clusters: Vec<(f64, f64)> = Vec::new();
    for (lo, hi) in left {
        match clusters.last_mut() {
            Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
            _ => clusters.push((lo, hi)),
        }
    }
    let mut budget = 16;
    let exact = match clusters.is_empty() {
        true => None,
        false => Some(Ring::<I>::new(c, s)?),
    };
    for cluster in clusters {
        let ring = exact.as_ref().expect("an exact ring for the clusters");
        let (roots, left) = search(ring, vec![cluster], &mut budget, 0.0)?;
        if !left.is_empty() {
            return Err(limit("a circle tangent to a cone"));
        }
        found.extend(roots);
    }
    found.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let mut out = Vec::new();
    for t in found {
        // Into (-pi, pi] where the root lies wholly beyond.
        let mut tt = span::<I>(t[0], t[1]);
        let tau = pi::<I>().mul(&I::exact_f64(2.0));
        if tt.cmp(&pi::<I>()) == Some(Ordering::Greater) {
            tt = tt.sub(&tau);
        } else if tt.cmp(&pi::<I>().neg()) != Some(Ordering::Greater)
            && tt.cmp(&pi::<I>().neg()).is_some()
        {
            tt = tt.add(&tau);
        }
        let (p, _) = fast.at(&span(t[0], t[1]));
        out.push(CurvePoint {
            parameter: bounds(&tt),
            point: bounds3(&p),
            tangent: false,
        });
    }
    Ok(CurveSurfaceIntersection::Points(out))
}

/// A circle and a cone in one tier, the cone's apex, axis and `cos^2`
/// computed once.
struct Ring<T> {
    o: E<T>,
    x: E<T>,
    y: E<T>,
    r: T,
    v: E<T>,
    a: E<T>,
    c2: T,
}

impl<T: Real> Ring<T> {
    fn new(c: &Circle, s: &Implicit) -> Result<Self> {
        let (v, a, c2) = s.cone::<T>()?;
        Ok(Self {
            o: e3(&c.o),
            x: e3(&c.x),
            y: e3(&c.y),
            r: T::from_r(&c.r),
            v,
            a,
            c2,
        })
    }
    /// The circle's point at `t` and `(G, G')` there.
    fn at(&self, t: &T) -> (E<T>, (T, T)) {
        let (ct, st) = T::cos_sin(t);
        let (o, x, y, r) = (&self.o, &self.x, &self.y, &self.r);
        let p: E<T> = std::array::from_fn(|i| o[i].add(&r.mul(&x[i].mul(&ct).add(&y[i].mul(&st)))));
        let dp: E<T> = std::array::from_fn(|i| r.mul(&y[i].mul(&ct).sub(&x[i].mul(&st))));
        let w = esub(&p, &self.v);
        let form = |x: &E<T>, y: &E<T>| {
            edot(x, y)
                .mul(&self.c2)
                .sub(&edot(x, &self.a).mul(&edot(y, &self.a)))
        };
        let g = (form(&w, &w), form(&w, &dp).mul(&T::exact_f64(2.0)));
        (p, g)
    }
}

/// One turn `[theta0 - pi, theta0 + pi]` whose ends lie inside a box around
/// `theta0 - pi` (mod a turn) where `G` is certainly nonzero, so no root is
/// found twice.
fn cut<T: Real>(ring: &Ring<T>) -> Result<(f64, f64)> {
    for t0 in [0.0, 1.0, 2.0, -1.0, 0.5, -2.5] {
        let lo = bounds(&I::exact_f64(t0).sub(&pi::<I>()));
        let hi = bounds(&I::exact_f64(t0).add(&pi::<I>()));
        let margin = 1e-9 * (1.0 + f64::abs(t0));
        let zone = span::<T>(lo[0] - margin, lo[1] + margin);
        let (_, (g, _)) = ring.at(&zone);
        if matches!(g.sign(), Some(Ordering::Less | Ordering::Greater)) {
            return Ok((lo[0], hi[1]));
        }
    }
    Err(limit("a circle against a cone"))
}

/// Certified roots, and the boxes left undecided.
type Search = (Vec<Enclosure>, Vec<(f64, f64)>);

/// The roots of `G` in `boxes`: those certified (a sign change and a
/// derivative of one sign, narrowed by bisection on certain signs), and the
/// boxes left undecided at `min_width` relative or at binary64 resolution.
fn search<T: Real>(
    ring: &Ring<T>,
    boxes: Vec<(f64, f64)>,
    budget: &mut usize,
    min_width: f64,
) -> Result<Search> {
    let certain = |x: &T| matches!(x.sign(), Some(Ordering::Less | Ordering::Greater));
    let sign = |x: &T| x.sign().filter(|o| *o != Ordering::Equal);
    let value = |t: f64| Ok::<T, Error>(ring.at(&T::exact_f64(t)).1 .0);
    let (mut roots, mut left) = (Vec::new(), Vec::new());
    let mut pending = boxes;
    while let Some((lo, hi)) = pending.pop() {
        *budget = budget
            .checked_sub(1)
            .ok_or(limit("a circle against a cone"))?;
        let m = 0.5 * lo + 0.5 * hi;
        let (_, (gb, db)) = ring.at(&span(lo, hi));
        if certain(&gb) {
            continue;
        }
        let mv = value(m)?.add(&db.mul(&span::<T>(lo, hi).sub(&T::exact_f64(m))));
        if certain(&mv) {
            continue;
        }
        let (sa, sb) = (sign(&value(lo)?), sign(&value(hi)?));
        if sa.is_some() && sb.is_some() && sa != sb && certain(&db) {
            let (mut a, mut b) = (lo, hi);
            for _ in 0..80 {
                let mid = 0.5 * a + 0.5 * b;
                if !(a < mid && mid < b) {
                    break;
                }
                match sign(&value(mid)?) {
                    Some(x) if Some(x) == sa => a = mid,
                    Some(_) => b = mid,
                    None => break,
                }
            }
            roots.push([a, b]);
            continue;
        }
        if hi - lo <= min_width * (1.0 + m.abs()) || !(lo < m && m < hi) {
            left.push((lo, hi));
            continue;
        }
        // Split off the centre, so that a root at a round angle (a turn's
        // middle, a quarter) never lies on a box's end.
        let split = lo + (hi - lo) * 0.4453125;
        let split = if lo < split && split < hi { split } else { m };
        pending.push((split, hi));
        pending.push((lo, split));
    }
    left.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok((roots, left))
}
