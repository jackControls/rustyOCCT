//! Intersections with a torus (S7b.3a of `REVIEW_NOTES.md`): a torus and a
//! plane or a sphere in any position, and every coaxial torus pair.
//!
//! A torus (axis `a` normalised exactly, major radius `R`, minor `r < R`) is
//! `C(phi) + r (cos t e + sin t a)`, `C = o + R e`, `e(phi) = cos phi x +
//! sin phi y` in an exactly orthonormal frame, `x` along the component normal
//! to the axis of the plane's stored normal or of the direction from the
//! torus's origin to the sphere's centre. On a meridian circle the plane or
//! sphere is `f0 + alpha cos t + beta sin t`, so each `phi` gives the points
//! of two branches, `(cos t, sin t) = (-alpha f0 -+ beta sqrt(D), -beta f0 +-
//! alpha sqrt(D)) / (alpha^2 + beta^2)`, `D = alpha^2 + beta^2 - f0^2`,
//! joined where `D` vanishes. `D` is a quadratic `P(c)` in `c = m cos phi`
//! with rational coefficients (`m^2` rational) and a negative leading one:
//! the classes come from the exact signs of `P(m)`, `P(-m)` and the vertex
//! against `+-m`, and a loop's ends are `arccos` of a root of `P` over `m`.
//!
//! Special cases, decided exactly: a plane normal to the axis, a sphere
//! centred on it and a coaxial cylinder, cone or torus meet the torus in
//! circles about its axis, where their meridians meet in a half-plane; a
//! plane containing the axis in two meridian circles. A sphere containing a
//! meridian circle, and a cylinder, cone or torus off the axis (S7b.3b), are
//! `NotConic`.
use super::analytic::{self, AnalyticItem, Enclosure, Enclosure3, SurfaceIntersection};
use super::procedural::{
    arccos, bounds, bounds3, cross, dot, e3, eadd, ecross, edot, escale, esub, i, limit,
    perpendicular, pi, q, span, sub, unit, zero, Branch, Component, ProceduralCurve, E, X,
};
use crate::certified::{Fast, Interval as I, Real};
use crate::topology::Surface;
use crate::{Error, Result};
use num_bigint::Sign;
use num_rational::BigRational as R;
use std::cmp::Ordering;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
}

/// A surface's exact data: origin, (unnormalised) axis or normal, radius (a
/// torus's major), a cone's half-angle and a torus's minor radius.
struct Data {
    kind: Kind,
    o: X,
    a: X,
    r: R,
    angle: f64,
    minor: R,
}

impl Data {
    fn of(s: &Surface) -> Option<Self> {
        let (kind, f, r, angle, minor) = match s {
            Surface::Plane(f) => (Kind::Plane, f, 0.0, 0.0, 0.0),
            Surface::Cylinder { frame, radius } => (Kind::Cylinder, frame, *radius, 0.0, 0.0),
            Surface::Cone {
                frame,
                radius,
                half_angle,
            } => (Kind::Cone, frame, *radius, *half_angle, 0.0),
            Surface::Sphere { frame, radius } => (Kind::Sphere, frame, *radius, 0.0, 0.0),
            Surface::Torus {
                frame,
                major,
                minor,
            } => (Kind::Torus, frame, *major, 0.0, *minor),
            Surface::BSpline(_) => return None,
        };
        Some(Self {
            kind,
            o: f.origin().to_array().map(q),
            a: f.normal().to_array().map(q),
            r: q(r),
            angle,
            minor: q(minor),
        })
    }
}

fn int(k: i64) -> R {
    R::from_integer(k.into())
}
fn negative(x: &R) -> bool {
    x.numer().sign() == Sign::Minus
}
fn sign_of(x: &R) -> Ordering {
    match x.numer().sign() {
        Sign::Minus => Ordering::Less,
        Sign::NoSign => Ordering::Equal,
        Sign::Plus => Ordering::Greater,
    }
}
fn all_zero(v: &X) -> bool {
    v.iter().all(zero)
}

/// The sign of `u + v sqrt(q)` (`q >= 0`), exactly.
fn surd_sign(u: &R, v: &R, q: &R) -> Ordering {
    let (su, sv) = (
        sign_of(u),
        if zero(q) { Ordering::Equal } else { sign_of(v) },
    );
    if sv == Ordering::Equal || su == sv {
        return if su == Ordering::Equal { sv } else { su };
    }
    if su == Ordering::Equal {
        return sv;
    }
    match (u * u).cmp(&(v * v * q)) {
        Ordering::Greater => su,
        Ordering::Less => sv,
        Ordering::Equal => Ordering::Equal,
    }
}

// ------------------------------------------------------------------ entry

/// The intersection of a torus and any analytic surface.
pub(super) fn intersect(a: &Surface, b: &Surface) -> Result<SurfaceIntersection> {
    let (Some(da), Some(db)) = (Data::of(a), Data::of(b)) else {
        return Err(Error::OutOfDomain(
            "analytic intersections of spline surfaces",
        ));
    };
    // The torus; of two, the first by its stored data (order independence).
    let first = match (da.kind, db.kind) {
        (Kind::Torus, Kind::Torus) => key(a) <= key(b),
        (Kind::Torus, _) => true,
        _ => false,
    };
    let (ts, os, t, o) = if first {
        (a, b, da, db)
    } else {
        (b, a, db, da)
    };
    let rel = sub(&o.o, &t.o);
    let parallel = all_zero(&cross(&o.a, &t.a));
    let on_axis = all_zero(&cross(&rel, &t.a));
    match o.kind {
        Kind::Plane => {
            if parallel {
                return normal_plane(&t, &o);
            }
            if zero(&dot(&o.a, &t.a)) && zero(&dot(&o.a, &rel)) {
                return meridian_circles(&t, &o.a);
            }
            general(ts, os, &t, &o)
        }
        Kind::Sphere => {
            if on_axis {
                let (za, d2) = axial(&t, &o.o);
                return circles_of(
                    &t,
                    circle_circle(
                        (i(&t.r), I::exact_f64(0.0)),
                        &t.minor,
                        (i(&int(0)), za),
                        &o.r,
                        &(&t.r * &t.r + d2),
                    ),
                );
            }
            general(ts, os, &t, &o)
        }
        Kind::Cylinder | Kind::Cone if !(parallel && on_axis) => {
            let curve = super::torus_curves::intersect(ts, os)?;
            Ok(if curve.tracks().is_empty() && curve.nodes().is_empty() {
                SurfaceIntersection::Empty
            } else {
                SurfaceIntersection::Traced(Box::new(curve))
            })
        }
        _ if !(parallel && on_axis) => Ok(SurfaceIntersection::NotConic),
        Kind::Cylinder => {
            let k = &o.r - &t.r;
            let r2 = &t.minor * &t.minor;
            let h = match (&k * &k).cmp(&r2) {
                Ordering::Greater => return Ok(SurfaceIntersection::Empty),
                Ordering::Equal => I::exact_f64(0.0),
                Ordering::Less => i::<I>(&(r2 - &k * &k)).sqrt(),
            };
            let mut heights = vec![h.clone()];
            if h.sign() != Some(Ordering::Equal) {
                heights.push(h.neg());
            }
            circles_of(&t, heights.into_iter().map(|z| (i(&o.r), z)).collect())
        }
        Kind::Torus => {
            let (za, za2) = axial(&t, &o.o);
            if o.r == t.r && zero(&za2) {
                return Ok(if o.minor == t.minor {
                    SurfaceIntersection::Same
                } else {
                    SurfaceIntersection::Empty
                });
            }
            let d = &o.r - &t.r;
            circles_of(
                &t,
                circle_circle(
                    (i(&t.r), I::exact_f64(0.0)),
                    &t.minor,
                    (i(&o.r), za),
                    &o.minor,
                    &(&d * &d + za2),
                ),
            )
        }
        _ => coaxial_cone(&t, &o),
    }
}

/// The comparison key of a torus: stored normal, origin, radii.
fn key(s: &Surface) -> [f64; 8] {
    let Surface::Torus {
        frame,
        major,
        minor,
    } = s
    else {
        unreachable!("tori only")
    };
    let [x, y, z] = frame.normal().to_array();
    let [u, v, w] = frame.origin().to_array();
    [x, y, z, u, v, w, *major, *minor]
}

/// A point's axial coordinate from the torus's origin, enclosed, and its
/// square exactly.
fn axial(t: &Data, p: &X) -> (I, R) {
    let h = dot(&sub(p, &t.o), &t.a);
    let aa = dot(&t.a, &t.a);
    let za = i::<I>(&h).div(&i::<I>(&aa).sqrt()).expect("a nonzero axis");
    (za, &h * &h / aa)
}

/// Circles about the torus's axis through meridian points `(rho, z)`.
fn circles_of(t: &Data, points: Vec<(I, I)>) -> Result<SurfaceIntersection> {
    let axis = unit::<I>(&t.a)?;
    let o = e3::<I>(&t.o);
    let items = points
        .iter()
        .map(|(rho, z)| analytic::circle(&eadd(&o, &escale(&axis, z)), &t.a, rho))
        .collect::<Result<Vec<_>>>()?;
    Ok(analytic::sorted(items))
}

/// Two circles in a meridian half-plane: centres `(rho, z)` enclosed, radii
/// and the squared distance of the centres exact (nonzero); their meeting
/// points by exact classes (tangent: one).
fn circle_circle(c1: (I, I), r1: &R, c2: (I, I), r2: &R, d2: &R) -> Vec<(I, I)> {
    let (s, f) = ((r1 + r2) * (r1 + r2), (r1 - r2) * (r1 - r2));
    if d2 > &s || d2 < &f {
        return Vec::new();
    }
    let d = i::<I>(d2).sqrt();
    let two = I::exact_f64(2.0);
    let div = |x: &I| x.div(&d).expect("distinct centres");
    let along = div(&i(&(d2 + r1 * r1 - r2 * r2))).div(&two).expect("two");
    let u = (div(&c2.0.sub(&c1.0)), div(&c2.1.sub(&c1.1)));
    let foot = (c1.0.add(&along.mul(&u.0)), c1.1.add(&along.mul(&u.1)));
    if d2 == &s || d2 == &f {
        return vec![foot];
    }
    let h = i::<I>(&(r1 * r1)).sub(&along.square()).sqrt();
    let n = (u.1.neg(), u.0);
    vec![
        (foot.0.add(&h.mul(&n.0)), foot.1.add(&h.mul(&n.1))),
        (foot.0.sub(&h.mul(&n.0)), foot.1.sub(&h.mul(&n.1))),
    ]
}

/// A plane normal to the axis at height `z`: `sin t = z / r`, circles of
/// radius `R +- sqrt(r^2 - z^2)`, one of radius `R` when tangent.
fn normal_plane(t: &Data, p: &Data) -> Result<SurfaceIntersection> {
    let n = &p.a;
    // The plane meets the axis at o + k a.
    let k = dot(n, &sub(&p.o, &t.o)) / dot(n, &t.a);
    let z2 = &k * &k * dot(&t.a, &t.a);
    let r2 = &t.minor * &t.minor;
    let z = i::<I>(&k).mul(&i::<I>(&dot(&t.a, &t.a)).sqrt());
    Ok(match z2.cmp(&r2) {
        Ordering::Greater => SurfaceIntersection::Empty,
        Ordering::Equal => return circles_of(t, vec![(i(&t.r), z)]),
        Ordering::Less => {
            let h = i::<I>(&(r2 - z2)).sqrt();
            return circles_of(
                t,
                vec![(i::<I>(&t.r).add(&h), z.clone()), (i::<I>(&t.r).sub(&h), z)],
            );
        }
    })
}

/// A plane containing the axis: the meridian circles at `o +- R u`,
/// `u = a x n` unit.
fn meridian_circles(t: &Data, n: &X) -> Result<SurfaceIntersection> {
    let u = unit::<I>(&cross(&t.a, n))?;
    let o = e3::<I>(&t.o);
    let items = [1.0, -1.0]
        .iter()
        .map(|s| {
            let centre = eadd(&o, &escale(&u, &i::<I>(&t.r).mul(&I::exact_f64(*s))));
            analytic::circle(&centre, n, &i(&t.minor))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(analytic::sorted(items))
}

/// A coaxial cone, both nappes: in a meridian half-plane the lines
/// `rho = s (rc + g (z - za) tan h)` (`g = +-1` as the cone's axis runs along
/// the torus's or against it) meet the meridian circle where a quadratic in
/// `z` has certified roots; a cone is never exactly tangent.
fn coaxial_cone(t: &Data, k: &Data) -> Result<SurfaceIntersection> {
    let (za, _) = axial(t, &k.o);
    let (c, s) = I::cos_sin(&I::exact_f64(k.angle));
    let tn = s.div(&c).ok_or(limit("a cone's slope"))?;
    let g = if negative(&dot(&k.a, &t.a)) {
        -1.0
    } else {
        1.0
    };
    let (big, small) = (i::<I>(&t.r), i::<I>(&t.minor));
    let mut points = Vec::new();
    for nappe in [1.0, -1.0] {
        let sn = I::exact_f64(nappe);
        let k0 = sn.mul(&i::<I>(&k.r).sub(&I::exact_f64(g).mul(&za).mul(&tn)));
        let k1 = sn.mul(&I::exact_f64(g)).mul(&tn);
        let a2 = k1.square().add(&I::exact_f64(1.0));
        let off = k0.sub(&big);
        let a1 = k1.mul(&off).mul(&I::exact_f64(2.0));
        let a0 = off.square().sub(&small.square());
        let disc = a1.square().sub(&a2.mul(&a0).mul(&I::exact_f64(4.0)));
        match disc.sign() {
            Some(Ordering::Less) => continue,
            Some(Ordering::Greater) => {}
            _ => return Err(limit("a coaxial cone's meridian")),
        }
        let root = disc.sqrt();
        for z in [a1.neg().add(&root), a1.neg().sub(&root)] {
            let z = z
                .div(&a2.mul(&I::exact_f64(2.0)))
                .ok_or(limit("a coaxial cone's meridian"))?;
            points.push((k0.add(&k1.mul(&z)), z));
        }
    }
    circles_of(t, points)
}

// ------------------------------------------------------------------ curves

/// `D` as an exact quadratic `(P2, P1, P0)` in `c = m cos phi` with `m^2 =
/// q`, and the frame's rational `x` direction.
fn quadratic(t: &Data, o: &Data) -> ([R; 3], R, X) {
    let aa = dot(&t.a, &t.a);
    let (big, small) = (&t.r, &t.minor);
    let r2 = small * small;
    match o.kind {
        Kind::Plane => {
            let n = &o.a;
            let na = dot(n, &t.a);
            let d0 = dot(n, &sub(&t.o, &o.o));
            let beta2 = &r2 * &na * &na / &aa;
            let q = dot(n, n) - &na * &na / &aa;
            // D = r^2 c^2 + beta^2 - (d0 + R c)^2.
            (
                [&r2 - big * big, -(int(2) * big * &d0), beta2 - &d0 * &d0],
                q,
                perpendicular(n, &t.a),
            )
        }
        _ => {
            let w = sub(&o.o, &t.o);
            let wa = dot(&w, &t.a);
            let ww = dot(&w, &w);
            let wa2 = &wa * &wa / &aa;
            let q = &ww - &wa2;
            // D = 4 r^2 (R - c)^2 + 4 r^2 wa^2 - (K - 2 R c)^2 with
            // K = R^2 + |w|^2 + r^2 - rho^2.
            let k = big * big + &ww + &r2 - &o.r * &o.r;
            let four = int(4);
            (
                [
                    &four * (&r2 - big * big),
                    &four * big * (&k - int(2) * &r2),
                    &four * &r2 * big * big + &four * &r2 * &wa2 - &k * &k,
                ],
                q,
                perpendicular(&w, &t.a),
            )
        }
    }
}

/// A torus and a plane or a sphere off the special cases.
fn general(ts: &Surface, os: &Surface, t: &Data, o: &Data) -> Result<SurfaceIntersection> {
    let ([p2, p1, p0], qv, _) = quadratic(t, o);
    let disc = &p1 * &p1 - int(4) * &p2 * &p0;
    let v = -&p1 / (int(2) * &p2);
    if zero(&disc) {
        // P = P2 (c - R)^2 (a sphere centred in the equatorial plane with
        // K = 2 R^2): negative on [-m, m] when R > m, else a whole meridian
        // circle lies on the sphere.
        return Ok(if &v * &v > qv {
            SurfaceIntersection::Empty
        } else {
            SurfaceIntersection::NotConic
        });
    }
    debug_assert!(negative(&p2) && !negative(&disc));
    let at = |s: i64| surd_sign(&(&p2 * &qv + &p0), &(&p1 * int(s)), &qv);
    let (plus, minus) = (at(1), at(-1));
    // The vertex against m and -m.
    let (v_plus, v_minus) = (surd_sign(&v, &int(-1), &qv), surd_sign(&v, &int(1), &qv));
    // The roots of P over m, and their arccos.
    let m = i::<I>(&qv).sqrt();
    let sq = i::<I>(&disc).sqrt();
    let den = i::<I>(&(int(2) * &p2));
    let angle = |root: I| -> Result<Enclosure> {
        let c = root
            .div(&den)
            .and_then(|x| x.div(&m))
            .ok_or(limit("a loop's cosine"))?;
        arccos(&c)
    };
    let lo = || angle(i::<I>(&-&p1).add(&sq));
    let hi = || angle(i::<I>(&-&p1).sub(&sq));
    let neg = |t: Enclosure| [-t[1], -t[0]];
    let tau = pi::<I>().mul(&I::exact_f64(2.0));
    let far = |t: Enclosure| {
        [
            bounds(&tau.sub(&I::exact_f64(t[1])))[0],
            bounds(&tau.sub(&I::exact_f64(t[0])))[1],
        ]
    };
    let (zero_, half, full) = ([0.0, 0.0], bounds(&pi::<I>()), bounds(&tau));
    let lp = |u0: Enclosure, u1: Enclosure| Component::Loop { u: [u0, u1] };
    use Ordering::{Equal as Z, Greater as P, Less as N};
    let components = match (plus, minus) {
        (P, P) => vec![
            Component::Ring {
                branch: Branch::Plus,
            },
            Component::Ring {
                branch: Branch::Minus,
            },
        ],
        (P, N) => {
            let t = lo()?;
            vec![lp(neg(t), t)]
        }
        (N, P) => {
            let t = hi()?;
            vec![lp(t, far(t))]
        }
        (P, Z) => vec![Component::FigureEight { node: half }],
        (Z, P) => vec![Component::FigureEight { node: zero_ }],
        (Z, Z) => vec![lp(zero_, half), lp(half, full)],
        (N, N) => {
            if v_plus == N && v_minus == P {
                let (a, b) = (lo()?, hi()?);
                vec![lp(neg(a), neg(b)), lp(b, a)]
            } else {
                return Ok(SurfaceIntersection::Empty);
            }
        }
        (Z, N) => {
            if v_plus == N {
                let t = lo()?;
                vec![lp(neg(t), zero_), lp(zero_, t)]
            } else {
                return tangent_point(ts, os, zero_);
            }
        }
        (N, Z) => {
            if v_minus == P {
                let t = hi()?;
                vec![lp(t, half), lp(half, far(t))]
            } else {
                return tangent_point(ts, os, half);
            }
        }
    };
    Ok(SurfaceIntersection::Procedural(Box::new(
        ProceduralCurve::new(ts.clone(), os.clone(), components),
    )))
}

fn tangent_point(ts: &Surface, os: &Surface, u: Enclosure) -> Result<SurfaceIntersection> {
    let p = point_at(ts, os, u, Branch::Plus)?;
    Ok(SurfaceIntersection::Items(vec![AnalyticItem::Point(p)]))
}

// ------------------------------------------------------------------ points

enum Other<T> {
    /// `n . (p - o)`.
    Plane { o: E<T>, n: E<T> },
    /// `|p - c|^2 - rho^2`.
    Sphere { c: E<T>, r2: T },
}

/// The meridian parameterisation in a certified tier.
struct Meridians<T> {
    o: E<T>,
    a: E<T>,
    x: E<T>,
    y: E<T>,
    big: T,
    small: T,
    other: Other<T>,
}

impl<T: Real> Meridians<T> {
    fn of(torus: &Surface, other: &Surface) -> Result<Self> {
        let (Some(t), Some(o)) = (Data::of(torus), Data::of(other)) else {
            return Err(Error::OutOfDomain("a torus's meridians"));
        };
        if t.kind != Kind::Torus || !matches!(o.kind, Kind::Plane | Kind::Sphere) {
            return Err(Error::OutOfDomain("a torus's meridians"));
        }
        let (_, _, toward) = quadratic(&t, &o);
        let a = unit::<T>(&t.a)?;
        let x = unit::<T>(&toward)?;
        let y = ecross(&a, &x);
        let other = match o.kind {
            Kind::Plane => Other::Plane {
                o: e3(&o.o),
                n: e3(&o.a),
            },
            _ => Other::Sphere {
                c: e3(&o.o),
                r2: i(&(&o.r * &o.r)),
            },
        };
        Ok(Self {
            o: e3(&t.o),
            a,
            x,
            y,
            big: i(&t.r),
            small: i(&t.minor),
            other,
        })
    }
    fn point(&self, u: &T, branch: Branch) -> Result<E<T>> {
        let (c, s) = T::cos_sin(u);
        let e = eadd(&escale(&self.x, &c), &escale(&self.y, &s));
        let centre = eadd(&self.o, &escale(&e, &self.big));
        let two = T::exact_f64(2.0);
        let (f0, alpha, beta) = match &self.other {
            Other::Plane { o, n } => (
                edot(n, &esub(&centre, o)),
                edot(n, &e).mul(&self.small),
                edot(n, &self.a).mul(&self.small),
            ),
            Other::Sphere { c, r2 } => {
                let rel = esub(&centre, c);
                (
                    edot(&rel, &rel).add(&self.small.square()).sub(r2),
                    edot(&rel, &e).mul(&self.small).mul(&two),
                    edot(&rel, &self.a).mul(&self.small).mul(&two),
                )
            }
        };
        let rho2 = alpha.square().add(&beta.square());
        let d = rho2.sub(&f0.square());
        if d.sign() == Some(Ordering::Less) {
            return Err(Error::OutOfDomain("a parameter off the curve"));
        }
        // sqrt clamps a negative part (the enclosure of a root) to zero.
        let root = match branch {
            Branch::Plus => d.sqrt(),
            Branch::Minus => d.sqrt().neg(),
        };
        let div = |x: T| x.div(&rho2).ok_or(limit("a meridian's point"));
        let cos_t = div(alpha.mul(&f0).add(&beta.mul(&root)).neg())?;
        let sin_t = div(alpha.mul(&root).sub(&beta.mul(&f0)))?;
        let tube = eadd(&escale(&e, &cos_t), &escale(&self.a, &sin_t));
        Ok(eadd(&centre, &escale(&tube, &self.small)))
    }
}

/// A certified enclosure of the points of a torus curve on `branch` at every
/// meridian angle in `u`: binary64 intervals, rational ones when those
/// cannot decide or are wider than `1e-12` relative at an exact parameter.
pub(super) fn point_at(
    torus: &Surface,
    other: &Surface,
    u: Enclosure,
    branch: Branch,
) -> Result<Enclosure3> {
    let [lo, hi] = u;
    let fast = Meridians::<Fast>::of(torus, other)
        .and_then(|s| s.point(&span(lo, hi), branch).map(|p| bounds3(&p)));
    let wide = |p: &Enclosure3| {
        lo == hi
            && p.iter()
                .any(|[a, b]| b - a > 1e-12 * a.abs().max(b.abs()).max(1.0))
    };
    match fast {
        Err(Error::ComputationLimit(_)) => {}
        Ok(p) if wide(&p) => {}
        other => return other,
    }
    let s = Meridians::<I>::of(torus, other)?;
    Ok(bounds3(&s.point(&span(lo, hi), branch)?))
}
