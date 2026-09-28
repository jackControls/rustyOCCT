//! Intersections of analytic surfaces whose result is empty, the same
//! surface, points, lines or a conic (S7a of `REVIEW_NOTES.md`).
//!
//! A surface is the exact point set its stored binary64 data define: a
//! plane through the stored origin with the stored normal; a cylinder,
//! cone or sphere about the line through the stored origin along the stored
//! normal, normalised exactly, with the stored radius and half-angle. A
//! configuration is degenerate (parallel, perpendicular, coaxial, tangent,
//! coincident, through an apex) only when it is exactly: every such decision
//! is an exact rational predicate on the stored data. A cone's `cos` and
//! `sin` are transcendental, so a plane is never exactly parallel to one of
//! its generatrices and no parabola occurs; its comparisons are certified in
//! rational intervals and report [`Error::ComputationLimit`] if they cannot
//! be separated. Every returned number is enclosed: binary64 bounds
//! `[lo, hi]` contain the exact value of each canonical parameter.
//!
//! Canonical forms: a line by its point nearest the global origin and a unit
//! direction; a circle by its centre, unit normal and radius; an ellipse by
//! its centre, unit normal, unit major direction, semi-major and semi-minor
//! axes; a hyperbola (both branches) by its centre, unit normal, unit
//! transverse direction, semi-transverse and semi-conjugate axes. A unit
//! direction's first nonzero coordinate is positive. Items are sorted by
//! kind, then by their numbers.
use crate::certified::{cos_sin, Interval as I, Real};
use crate::topology::Surface;
use crate::{Error, Result};
use num_bigint::Sign;
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// Binary64 bounds `[lo, hi]` of an exact real.
pub type Enclosure = [f64; 2];
/// Bounds of each coordinate of an exact point or vector.
pub type Enclosure3 = [Enclosure; 3];

/// One piece of an intersection, with enclosed canonical parameters.
#[derive(Debug, Clone, PartialEq)]
pub enum AnalyticItem {
    Point(Enclosure3),
    Line {
        point: Enclosure3,
        direction: Enclosure3,
    },
    Circle {
        centre: Enclosure3,
        normal: Enclosure3,
        radius: Enclosure,
    },
    Ellipse {
        centre: Enclosure3,
        normal: Enclosure3,
        major: Enclosure3,
        semi_major: Enclosure,
        semi_minor: Enclosure,
    },
    /// Both branches.
    Hyperbola {
        centre: Enclosure3,
        normal: Enclosure3,
        transverse: Enclosure3,
        semi_transverse: Enclosure,
        semi_conjugate: Enclosure,
    },
}

impl AnalyticItem {
    pub fn kind(&self) -> &'static str {
        match self {
            Self::Point(_) => "point",
            Self::Line { .. } => "line",
            Self::Circle { .. } => "circle",
            Self::Ellipse { .. } => "ellipse",
            Self::Hyperbola { .. } => "hyperbola",
        }
    }
    /// Every enclosed number in canonical order.
    pub fn values(&self) -> Vec<Enclosure> {
        let mut out = Vec::new();
        let mut v3 = |v: &Enclosure3| out.extend(v.iter().copied());
        match self {
            Self::Point(p) => v3(p),
            Self::Line { point, direction } => {
                v3(point);
                v3(direction);
            }
            Self::Circle {
                centre,
                normal,
                radius,
            } => {
                v3(centre);
                v3(normal);
                out.push(*radius);
            }
            Self::Ellipse {
                centre,
                normal,
                major,
                semi_major,
                semi_minor,
            } => {
                v3(centre);
                v3(normal);
                v3(major);
                out.extend([*semi_major, *semi_minor]);
            }
            Self::Hyperbola {
                centre,
                normal,
                transverse,
                semi_transverse,
                semi_conjugate,
            } => {
                v3(centre);
                v3(normal);
                v3(transverse);
                out.extend([*semi_transverse, *semi_conjugate]);
            }
        }
        out
    }
}

/// The intersection of two analytic surfaces.
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceIntersection {
    Empty,
    /// The two surfaces are the same point set.
    Same,
    /// Points, lines and conics, canonically ordered.
    Items(Vec<AnalyticItem>),
    /// A procedural curve (D13, S7b): the other quadric pairs, and a torus
    /// with a plane or a sphere.
    Procedural(Box<super::procedural::ProceduralCurve>),
    /// A traced curve (D13, S7b.3b): a torus and a cylinder or a cone off its
    /// axis, as a graph of certified tracks, folds and tangencies.
    Traced(Box<super::torus_curves::TracedCurve>),
    /// A curve that is not a conic, of a pair not yet parameterised (two
    /// cones, a cone's apex on a sphere, a sphere containing a torus's
    /// meridian circle, two tori off a common axis).
    NotConic,
}

// ------------------------------------------------------------------ exact

type X = [R; 3];

/// Exact zero and sign tests of a rational.
trait Exact {
    fn is_zero(&self) -> bool;
    fn is_positive(&self) -> bool;
}
impl Exact for R {
    fn is_zero(&self) -> bool {
        self.numer().sign() == Sign::NoSign
    }
    fn is_positive(&self) -> bool {
        self.numer().sign() == Sign::Plus
    }
}
fn r0() -> R {
    R::from_integer(0.into())
}

fn q(x: f64) -> R {
    R::from_float(x).expect("finite surface data")
}
fn x3(v: [f64; 3]) -> X {
    v.map(q)
}
fn dot(a: &X, b: &X) -> R {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}
fn cross(a: &X, b: &X) -> X {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
fn sub(a: &X, b: &X) -> X {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}
fn add(a: &X, b: &X) -> X {
    [&a[0] + &b[0], &a[1] + &b[1], &a[2] + &b[2]]
}
fn scale(a: &X, k: &R) -> X {
    [&a[0] * k, &a[1] * k, &a[2] * k]
}
fn is_zero(v: &X) -> bool {
    v.iter().all(Exact::is_zero)
}

// ------------------------------------------------------------------ enclosed

type E = [I; 3];

fn i(x: &R) -> I {
    I::exact(x.clone())
}
fn e3(v: &X) -> E {
    [i(&v[0]), i(&v[1]), i(&v[2])]
}
fn edot(a: &E, b: &E) -> I {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}
fn eadd(a: &E, b: &E) -> E {
    [a[0].add(&b[0]), a[1].add(&b[1]), a[2].add(&b[2])]
}
fn esub(a: &E, b: &E) -> E {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}
fn escale(a: &E, k: &I) -> E {
    [a[0].mul(k), a[1].mul(k), a[2].mul(k)]
}
fn ediv(a: &E, k: &I) -> Result<E> {
    let d = |x: &I| {
        x.div(k)
            .ok_or(Error::ComputationLimit("certified division"))
    };
    Ok([d(&a[0])?, d(&a[1])?, d(&a[2])?])
}
fn bounds(x: &I) -> Enclosure {
    let (lo, hi) = x.bounds_f64();
    [lo, hi]
}
fn bounds3(v: &E) -> Enclosure3 {
    [bounds(&v[0]), bounds(&v[1]), bounds(&v[2])]
}
fn sign(x: &I) -> Result<Ordering> {
    x.sign()
        .ok_or(Error::ComputationLimit("certified comparison"))
}

/// A unit vector along an exact direction, its first nonzero coordinate
/// positive.
fn unit(v: &X) -> Result<E> {
    let lead = v
        .iter()
        .find(|x| !x.is_zero())
        .ok_or(Error::Degenerate("direction"))?;
    let norm = i(&dot(v, v)).sqrt();
    let s = if lead.is_positive() {
        I::exact(R::from_integer(1.into()))
    } else {
        I::exact(R::from_integer((-1).into()))
    };
    ediv(&escale(&e3(v), &s), &norm)
}

/// A unit vector along an enclosed direction, its first coordinate that
/// certainly differs from zero made positive (the coordinates before it are
/// exactly zero where the direction is built from them).
fn unit_enclosed(v: &E) -> Result<E> {
    let lead = v
        .iter()
        .find(|x| x.sign().is_some_and(|s| s != Ordering::Equal))
        .ok_or(Error::ComputationLimit("direction sign"))?;
    let norm = edot(v, v).sqrt();
    let norm = if sign(lead)? == Ordering::Less {
        norm.neg()
    } else {
        norm
    };
    ediv(v, &norm)
}

/// The point of a line nearest the global origin.
fn nearest(p: &E, d: &E) -> Result<E> {
    let k = edot(p, d)
        .div(&edot(d, d))
        .ok_or(Error::ComputationLimit("certified division"))?;
    Ok(esub(p, &escale(d, &k)))
}

pub(super) fn line(p: &E, d: &E) -> Result<AnalyticItem> {
    Ok(AnalyticItem::Line {
        point: bounds3(&nearest(p, d)?),
        direction: bounds3(&unit_enclosed(d)?),
    })
}
fn exact_line(p: &X, d: &X) -> Result<AnalyticItem> {
    Ok(AnalyticItem::Line {
        point: bounds3(&nearest(&e3(p), &e3(d))?),
        direction: bounds3(&unit(d)?),
    })
}
pub(super) fn circle(centre: &E, normal: &X, radius: &I) -> Result<AnalyticItem> {
    Ok(AnalyticItem::Circle {
        centre: bounds3(centre),
        normal: bounds3(&unit(normal)?),
        radius: bounds(radius),
    })
}

// ------------------------------------------------------------------ surfaces

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Kind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
}

/// A surface's exact data: origin, (unnormalised) axis or normal, radius
/// and a cone's half-angle.
struct Q {
    kind: Kind,
    o: X,
    a: X,
    r: R,
    angle: f64,
}

impl Q {
    fn of(s: &Surface) -> Option<Self> {
        let (kind, f, r, angle) = match s {
            Surface::Plane(f) => (Kind::Plane, f, 0.0, 0.0),
            Surface::Cylinder { frame, radius } => (Kind::Cylinder, frame, *radius, 0.0),
            Surface::Cone {
                frame,
                radius,
                half_angle,
            } => (Kind::Cone, frame, *radius, *half_angle),
            Surface::Sphere { frame, radius } => (Kind::Sphere, frame, *radius, 0.0),
            Surface::Torus { .. } | Surface::BSpline(_) => return None,
        };
        Some(Self {
            kind,
            o: x3(f.origin().to_array()),
            a: x3(f.normal().to_array()),
            r: q(r),
            angle,
        })
    }
    /// (cos, sin) of a cone's half-angle, enclosed.
    fn cos_sin(&self) -> (I, I) {
        cos_sin(&q(self.angle))
    }
}

/// The intersection of two analytic surfaces: points, lines and conics
/// (S7a), procedural curves (S7b), and [`SurfaceIntersection::NotConic`]
/// for the pairs not yet parameterised; a spline surface is out of domain.
pub fn surface_surface(a: &Surface, b: &Surface) -> Result<SurfaceIntersection> {
    if matches!(a, Surface::Torus { .. }) || matches!(b, Surface::Torus { .. }) {
        return super::toroidal::intersect(a, b);
    }
    let (sa, sb) = (a, b);
    let (Some(mut a), Some(mut b)) = (Q::of(a), Q::of(b)) else {
        if matches!(a, Surface::BSpline(_)) || matches!(b, Surface::BSpline(_)) {
            return Err(Error::OutOfDomain(
                "analytic intersections of spline surfaces",
            ));
        }
        return Ok(SurfaceIntersection::NotConic);
    };
    if a.kind > b.kind {
        std::mem::swap(&mut a, &mut b);
    }
    let items = match (a.kind, b.kind) {
        (Kind::Plane, Kind::Plane) => plane_plane(&a, &b)?,
        (Kind::Plane, Kind::Sphere) => plane_sphere(&a, &b)?,
        (Kind::Plane, Kind::Cylinder) => plane_cylinder(&a, &b)?,
        (Kind::Plane, Kind::Cone) => plane_cone(&a, &b)?,
        (Kind::Sphere, Kind::Sphere) => sphere_sphere(&a, &b)?,
        _ => {
            let on_axis = |p: &X, q: &Q| is_zero(&cross(&sub(p, &q.o), &q.a));
            let coaxial = if b.kind == Kind::Sphere {
                on_axis(&b.o, &a)
            } else {
                is_zero(&cross(&a.a, &b.a)) && on_axis(&b.o, &a)
            };
            if coaxial {
                coaxial_pair(&a, &b)?
            } else if (a.kind, b.kind) == (Kind::Cylinder, Kind::Cylinder) {
                cylinder_cylinder(&a, &b)?
            } else {
                Out::NotConic
            }
        }
    };
    if let Out::NotConic = items {
        if let Some(found) = super::ruled_curves::intersect(sa, sb)? {
            return Ok(found);
        }
        match super::procedural::intersect(sa, sb)? {
            Some(super::procedural::Found::Empty) => return Ok(SurfaceIntersection::Empty),
            Some(super::procedural::Found::Point(p)) => {
                return Ok(SurfaceIntersection::Items(vec![AnalyticItem::Point(p)]))
            }
            Some(super::procedural::Found::Curve(c)) => {
                return Ok(SurfaceIntersection::Procedural(c))
            }
            None => {}
        }
    }
    Ok(match items {
        Out::Items(items) => sorted(items),
        Out::Same => SurfaceIntersection::Same,
        Out::NotConic => SurfaceIntersection::NotConic,
    })
}

/// Items in canonical order (by kind, then by their numbers), or `Empty`.
pub(super) fn sorted(mut items: Vec<AnalyticItem>) -> SurfaceIntersection {
    let key = |x: &AnalyticItem| {
        (
            x.kind(),
            x.values()
                .iter()
                .map(|[lo, hi]| 0.5 * lo + 0.5 * hi)
                .collect::<Vec<_>>(),
        )
    };
    items.sort_by(|x, y| {
        let (kx, vx) = key(x);
        let (ky, vy) = key(y);
        kx.cmp(ky).then_with(|| {
            vx.iter()
                .zip(&vy)
                .map(|(a, b)| a.total_cmp(b))
                .find(|o| *o != Ordering::Equal)
                .unwrap_or(Ordering::Equal)
        })
    });
    if items.is_empty() {
        SurfaceIntersection::Empty
    } else {
        SurfaceIntersection::Items(items)
    }
}

/// A pair's outcome before sorting.
enum Out {
    Same,
    Items(Vec<AnalyticItem>),
    NotConic,
}
type Found = Result<Out>;

fn plane_plane(p1: &Q, p2: &Q) -> Found {
    let d = cross(&p1.a, &p2.a);
    if is_zero(&d) {
        return Ok(if dot(&sub(&p2.o, &p1.o), &p1.a).is_zero() {
            Out::Same
        } else {
            Out::Items(Vec::new())
        });
    }
    // The point x n1 + y n2 on both planes.
    let (a11, a12, a22) = (dot(&p1.a, &p1.a), dot(&p1.a, &p2.a), dot(&p2.a, &p2.a));
    let (b1, b2) = (dot(&p1.o, &p1.a), dot(&p2.o, &p2.a));
    let det = &a11 * &a22 - &a12 * &a12;
    let x = (&b1 * &a22 - &b2 * &a12) / &det;
    let y = (&a11 * &b2 - &a12 * &b1) / &det;
    let p = add(&scale(&p1.a, &x), &scale(&p2.a, &y));
    Ok(Out::Items(vec![exact_line(&p, &d)?]))
}

fn plane_sphere(p: &Q, s: &Q) -> Found {
    let n = &p.a;
    let off = dot(&sub(&s.o, &p.o), n);
    let nn = dot(n, n);
    let d2 = &off * &off / &nn;
    let r2 = &s.r * &s.r;
    let foot = sub(&s.o, &scale(n, &(&off / &nn)));
    Ok(Out::Items(match d2.cmp(&r2) {
        Ordering::Greater => Vec::new(),
        Ordering::Equal => vec![AnalyticItem::Point(bounds3(&e3(&foot)))],
        Ordering::Less => vec![circle(&e3(&foot), n, &i(&(r2 - d2)).sqrt())?],
    }))
}

fn plane_cylinder(p: &Q, c: &Q) -> Found {
    let (n, a) = (&p.a, &c.a);
    let t = dot(a, n);
    if t.is_zero() {
        // The axis is parallel to the plane: lines along it, offset within
        // the plane perpendicular to it.
        let off = dot(&sub(&c.o, &p.o), n);
        let nn = dot(n, n);
        let d2 = &off * &off / &nn;
        let r2 = &c.r * &c.r;
        let foot = sub(&c.o, &scale(n, &(&off / &nn)));
        return Ok(Out::Items(match d2.cmp(&r2) {
            Ordering::Greater => Vec::new(),
            Ordering::Equal => vec![exact_line(&foot, a)?],
            Ordering::Less => {
                let w = cross(n, a);
                let h = i(&((r2 - d2) / dot(&w, &w))).sqrt();
                let mut out = Vec::new();
                for k in [h.clone(), h.neg()] {
                    let point = eadd(&e3(&foot), &escale(&e3(&w), &k));
                    out.push(line(&point, &e3(a))?);
                }
                out
            }
        }));
    }
    // Where the axis meets the plane.
    let lambda = dot(&sub(&p.o, &c.o), n) / &t;
    let centre = e3(&add(&c.o, &scale(a, &lambda)));
    let w = cross(a, n);
    if is_zero(&w) {
        return Ok(Out::Items(vec![circle(&centre, n, &i(&c.r))?]));
    }
    // Semi-minor r along a x n; semi-major r |a| |n| / |a.n| along the
    // axis's projection n x (a x n).
    let major = cross(n, &w);
    let semi_major = i(&c.r).mul(&i(&(dot(a, a) * dot(n, n) / (&t * &t))).sqrt());
    Ok(Out::Items(vec![AnalyticItem::Ellipse {
        centre: bounds3(&centre),
        normal: bounds3(&unit(n)?),
        major: bounds3(&unit(&major)?),
        semi_major: bounds(&semi_major),
        semi_minor: bounds(&i(&c.r)),
    }]))
}

/// A plane and a double cone (both nappes): `rho cos a = +-(r cos a + h sin
/// a)`, `h` along the unit axis from the cone's origin.
fn plane_cone(p: &Q, k: &Q) -> Found {
    let (n, a) = (&p.a, &k.a);
    let (ca, sa) = k.cos_sin();
    let aa = dot(a, a);
    let nn = dot(n, n);
    let an = dot(a, n);
    let la = i(&aa).sqrt();
    let unit_a = ediv(&e3(a), &la)?;
    let rel = sub(&k.o, &p.o);
    // Axis normal to the plane: a circle, or the apex when the plane passes
    // through a rational apex.
    if is_zero(&cross(a, n)) {
        // h of the plane from the cone's origin along the unit axis.
        let hq = dot(&sub(&p.o, &k.o), a);
        if k.r.is_zero() && hq.is_zero() {
            return Ok(Out::Items(vec![AnalyticItem::Point(bounds3(&e3(&k.o)))]));
        }
        let centre = add(&k.o, &scale(a, &(&hq / &aa)));
        let h = i(&hq).div(&la).expect("a nonzero axis");
        let tan = sa
            .div(&ca)
            .ok_or(Error::ComputationLimit("cone half-angle"))?;
        let rho = i(&k.r).add(&h.mul(&tan));
        let rho = if sign(&rho)? == Ordering::Less {
            rho.neg()
        } else {
            rho
        };
        return Ok(Out::Items(vec![circle(&e3(&centre), n, &rho)?]));
    }
    // The apex, rational when the radius at the origin is zero.
    let tan = sa
        .div(&ca)
        .ok_or(Error::ComputationLimit("cone half-angle"))?;
    let apex = if k.r.is_zero() {
        e3(&k.o)
    } else {
        let back = i(&k.r)
            .div(&tan)
            .ok_or(Error::ComputationLimit("cone half-angle"))?;
        esub(&e3(&k.o), &escale(&unit_a, &back))
    };
    let through_apex = dot(&rel, n).is_zero() && (k.r.is_zero() || an.is_zero());
    // cos b = n^.a^ (signed), sin b = |e1| / |a| with e1 the axis's
    // projection on the plane, e2 = n x e1.
    let e1 = sub(a, &scale(n, &(&an / &nn)));
    let e2 = cross(n, &e1);
    let ln = i(&nn).sqrt();
    let cb = i(&an).div(&la.mul(&ln)).expect("nonzero axes");
    let l1 = i(&dot(&e1, &e1)).sqrt();
    let u1 = ediv(&e3(&e1), &l1)?;
    let u2 = ediv(&e3(&e2), &i(&dot(&e2, &e2)).sqrt())?;
    let sb = l1.div(&la).expect("a nonzero axis");
    let s2 = sa.square();
    // A = cos^2 b - sin^2 a: positive for an ellipse (or the apex), negative
    // for a hyperbola (or two lines); zero (a parabola) is impossible.
    let big_a = cb.square().sub(&s2);
    let kind = sign(&big_a)?;
    if kind == Ordering::Equal {
        return Err(Error::ComputationLimit("a parabolic section"));
    }
    if through_apex {
        if kind == Ordering::Greater {
            return Ok(Out::Items(vec![AnalyticItem::Point(bounds3(&apex))]));
        }
        // Two generatrices in the plane at angle a from the axis: d = cos p
        // u1 +- sin p u2 with cos p = cos a / sin b.
        let cp = ca.div(&sb).ok_or(Error::ComputationLimit("cone section"))?;
        let cp = if sign(&cp)? == Ordering::Less {
            cp.neg()
        } else {
            cp
        };
        let sp = I::exact(R::from_integer(1.into())).sub(&cp.square()).sqrt();
        let mut out = Vec::new();
        for s in [sp.clone(), sp.neg()] {
            let d = eadd(&escale(&u1, &cp), &escale(&u2, &s));
            out.push(line(&apex, &d)?);
        }
        return Ok(Out::Items(out));
    }
    // The signed distance of the apex from the plane, and the foot of the
    // apex on it.
    let unit_n = ediv(&e3(n), &ln)?;
    let dist = edot(&esub(&apex, &e3(&p.o)), &unit_n);
    let foot = esub(&apex, &escale(&unit_n, &dist));
    // A (s - s0)^2 + cos^2 a t^2 + K = 0, s0 = -D cos b sin b / A,
    // K = -D^2 cos^2 a sin^2 a / A.
    let d2 = dist.square();
    let c2 = ca.square();
    let s0 = dist
        .mul(&cb)
        .mul(&sb)
        .neg()
        .div(&big_a)
        .expect("a nonzero A");
    let centre = eadd(&foot, &escale(&u1, &s0));
    let abs = |x: I| -> Result<I> {
        Ok(if sign(&x)? == Ordering::Less {
            x.neg()
        } else {
            x
        })
    };
    // Along e1: |D| |cos a sin a| / |A|; along e2: |D| |sin a| / sqrt|A|.
    let along1 = abs(d2.mul(&c2).mul(&s2))?
        .sqrt()
        .div(&abs(big_a.clone())?)
        .expect("a nonzero A");
    let along2 = abs(d2.mul(&s2))?
        .sqrt()
        .div(&abs(big_a)?.sqrt())
        .expect("a nonzero A");
    let normal = bounds3(&unit(n)?);
    let dir = bounds3(&unit(&e1)?);
    Ok(Out::Items(vec![if kind == Ordering::Greater {
        AnalyticItem::Ellipse {
            centre: bounds3(&centre),
            normal,
            major: dir,
            semi_major: bounds(&along1),
            semi_minor: bounds(&along2),
        }
    } else {
        AnalyticItem::Hyperbola {
            centre: bounds3(&centre),
            normal,
            transverse: dir,
            semi_transverse: bounds(&along1),
            semi_conjugate: bounds(&along2),
        }
    }]))
}

fn sphere_sphere(s1: &Q, s2: &Q) -> Found {
    let d = sub(&s2.o, &s1.o);
    let dd = dot(&d, &d);
    let (r1, r2) = (&s1.r, &s2.r);
    if dd.is_zero() {
        return Ok(if r1 == r2 {
            Out::Same
        } else {
            Out::Items(Vec::new())
        });
    }
    let sum = r1 + r2;
    let diff = r1 - r2;
    if dd > &sum * &sum || dd < &diff * &diff {
        return Ok(Out::Items(Vec::new()));
    }
    // The radical plane at c1 + d k.
    let k = (&dd + r1 * r1 - r2 * r2) / (R::from_integer(2.into()) * &dd);
    let centre = add(&s1.o, &scale(&d, &k));
    let rad2 = r1 * r1 - &k * &k * &dd;
    Ok(Out::Items(if rad2.is_zero() {
        vec![AnalyticItem::Point(bounds3(&e3(&centre)))]
    } else {
        vec![circle(&e3(&centre), &d, &i(&rad2).sqrt())?]
    }))
}

fn cylinder_cylinder(c1: &Q, c2: &Q) -> Found {
    let (a1, a2) = (&c1.a, &c2.a);
    if is_zero(&cross(a1, a2)) {
        // Parallel, not coaxial: as two circles in a cross-section.
        let w = sub(&c2.o, &c1.o);
        let w = sub(&w, &scale(a1, &(dot(&w, a1) / dot(a1, a1))));
        let dd = dot(&w, &w);
        let (r1, r2) = (&c1.r, &c2.r);
        let sum = r1 + r2;
        let diff = r1 - r2;
        if dd > &sum * &sum || dd < &diff * &diff {
            return Ok(Out::Items(Vec::new()));
        }
        let k = (&dd + r1 * r1 - r2 * r2) / (R::from_integer(2.into()) * &dd);
        let foot = add(&c1.o, &scale(&w, &k));
        let h2 = r1 * r1 - &k * &k * &dd;
        if h2.is_zero() {
            return Ok(Out::Items(vec![exact_line(&foot, a1)?]));
        }
        let v = cross(a1, &w);
        let h = i(&(h2 / dot(&v, &v))).sqrt();
        let mut out = Vec::new();
        for s in [h.clone(), h.neg()] {
            out.push(line(&eadd(&e3(&foot), &escale(&e3(&v), &s)), &e3(a1))?);
        }
        return Ok(Out::Items(out));
    }
    // Crossing axes of equal radii: two ellipses in the bisecting planes.
    let m = cross(a1, a2);
    let w = sub(&c2.o, &c1.o);
    if !dot(&w, &m).is_zero() || c1.r != c2.r {
        return Ok(Out::NotConic);
    }
    let u = dot(&cross(&w, a2), &m) / dot(&m, &m);
    let x = e3(&add(&c1.o, &scale(a1, &u)));
    let u1 = ediv(&e3(a1), &i(&dot(a1, a1)).sqrt())?;
    let u2 = ediv(&e3(a2), &i(&dot(a2, a2)).sqrt())?;
    let r = i(&c1.r);
    let mut out = Vec::new();
    for plus in [true, false] {
        let (b, normal) = if plus {
            (eadd(&u1, &u2), esub(&u1, &u2))
        } else {
            (esub(&u1, &u2), eadd(&u1, &u2))
        };
        // Semi-major r over the cosine between an axis and the plane's
        // normal, along the bisector; semi-minor r along the common normal.
        let cos = edot(&u1, &normal).div(&edot(&normal, &normal).sqrt());
        let cos = cos.ok_or(Error::ComputationLimit("cylinder section"))?;
        let cos = if sign(&cos)? == Ordering::Less {
            cos.neg()
        } else {
            cos
        };
        out.push(AnalyticItem::Ellipse {
            centre: bounds3(&x),
            normal: bounds3(&unit_enclosed(&normal)?),
            major: bounds3(&unit_enclosed(&b)?),
            semi_major: bounds(
                &r.div(&cos)
                    .ok_or(Error::ComputationLimit("cylinder section"))?,
            ),
            semi_minor: bounds(&r),
        });
    }
    Ok(Out::Items(out))
}

/// Coaxial cylinders and cones, and a sphere centred on the other's axis:
/// circles where their radius functions of the axial coordinate agree.
fn coaxial_pair(q1: &Q, q2: &Q) -> Found {
    let axis = if q1.kind == Kind::Sphere { q2 } else { q1 };
    let a = &axis.a;
    let aa = dot(a, a);
    let la = i(&aa).sqrt();
    // h: signed distance along the unit axis from the axis surface's origin.
    let h_of = |p: &X| {
        i(&dot(&sub(p, &axis.o), a))
            .div(&la)
            .expect("a nonzero axis")
    };
    let at = |h: &I| {
        eadd(
            &e3(&axis.o),
            &escale(&e3(a), &h.div(&la).expect("a nonzero axis")),
        )
    };
    // Rho as a linear function r0 + t h (cylinders and cones), or a sphere.
    enum F {
        Lin(I, I, bool),
        Sph(I, I),
    }
    let f = |s: &Q| -> Result<F> {
        let h0 = h_of(&s.o);
        Ok(match s.kind {
            Kind::Cylinder => F::Lin(i(&s.r), I::exact(r0()), true),
            Kind::Sphere => F::Sph(h0, i(&s.r)),
            _ => {
                let (c, sn) = s.cos_sin();
                let t = sn
                    .div(&c)
                    .ok_or(Error::ComputationLimit("cone half-angle"))?;
                F::Lin(i(&s.r).sub(&h0.mul(&t)), t, false)
            }
        })
    };
    let (f1, f2) = (f(q1)?, f(q2)?);
    let mut roots: Vec<I> = Vec::new();
    match (&f1, &f2) {
        (F::Lin(c1, t1, cyl1), F::Lin(c2, t2, cyl2)) => {
            // rho1 = +-rho2: (t1 -+ t2) h = -(c1 -+ c2).
            let same_slope = |plus: bool| -> bool {
                if *cyl1 && *cyl2 {
                    return true;
                }
                if *cyl1 || *cyl2 {
                    return false;
                }
                let (x, y) = (q(q1.angle), q(q2.angle));
                if plus {
                    x == y
                } else {
                    x == -y
                }
            };
            for plus in [true, false] {
                let (dt, dc) = if plus {
                    (t1.sub(t2), c1.sub(c2))
                } else {
                    (t1.add(t2), c1.add(c2))
                };
                if same_slope(plus) {
                    // Equal slopes: the same surface when the constants
                    // agree exactly, else no root on this branch.
                    if exact_constants_equal(q1, q2, plus) {
                        return Ok(Out::Same);
                    }
                    continue;
                }
                roots.push(
                    dc.neg()
                        .div(&dt)
                        .ok_or(Error::ComputationLimit("coaxial root"))?,
                );
            }
        }
        (F::Lin(c, t, _), F::Sph(h0, big_r)) | (F::Sph(h0, big_r), F::Lin(c, t, _)) => {
            // (c + t h)^2 = R^2 - (h - h0)^2: (1 + t^2) h^2 + 2 (c t - h0) h
            // + c^2 + h0^2 - R^2 = 0.
            let one = I::exact(R::from_integer(1.into()));
            let qa = one.add(&t.square());
            let qb = c.mul(t).sub(h0);
            let qc = c.square().add(&h0.square()).sub(&big_r.square());
            let disc = qb.square().sub(&qa.mul(&qc));
            match exact_discriminant(q1, q2) {
                Some(Ordering::Less) => {}
                Some(Ordering::Equal) => {
                    roots.push(qb.neg().div(&qa).expect("a positive coefficient"));
                }
                Some(Ordering::Greater) | None => match sign(&disc)? {
                    Ordering::Less => {}
                    Ordering::Equal => {
                        return Err(Error::ComputationLimit("a tangent coaxial sphere"))
                    }
                    Ordering::Greater => {
                        let root = disc.sqrt();
                        for s in [root.clone(), root.neg()] {
                            roots.push(qb.neg().add(&s).div(&qa).expect("a positive coefficient"));
                        }
                    }
                },
            }
        }
        (F::Sph(..), F::Sph(..)) => unreachable!("sphere pairs are handled above"),
    }
    let rho_at = |h: &I| -> I {
        match &f1 {
            F::Lin(c, t, _) => c.add(&t.mul(h)),
            F::Sph(h0, r) => r.square().sub(&h.sub(h0).square()).sqrt(),
        }
    };
    let mut out = Vec::new();
    let mut apex_found = false;
    for h in roots {
        let rho = rho_at(&h);
        let centre = at(&h);
        if apex_root(q1, q2, &h) {
            // Two branches through a common apex meet there once.
            if !apex_found {
                out.push(AnalyticItem::Point(bounds3(&centre)));
                apex_found = true;
            }
            continue;
        }
        let rho = if sign(&rho)? == Ordering::Less {
            rho.neg()
        } else {
            rho
        };
        out.push(circle(&centre, a, &rho)?);
    }
    Ok(Out::Items(out))
}

/// Whether two coaxial linear radius functions with equal slopes (the same
/// or opposite half-angle, or two cylinders) have equal constants exactly:
/// then they are the same surface.
fn exact_constants_equal(q1: &Q, q2: &Q, plus: bool) -> bool {
    let a = &q1.a;
    match (q1.kind, q2.kind) {
        (Kind::Cylinder, Kind::Cylinder) => q1.r == q2.r,
        _ => {
            // r1 - h01 t = +-(r2 - h02 t) with t transcendental: exact only
            // when the rational parts and the h parts agree separately.
            let h01 = dot(&sub(&q1.o, &q1.o), a);
            let h02 = dot(&sub(&q2.o, &q1.o), a);
            if plus {
                q1.r == q2.r && h01 == h02
            } else {
                (&q1.r + &q2.r).is_zero() && h01 == h02
            }
        }
    }
}

/// An exact decision of a coaxial sphere pair's discriminant where it is
/// rational: a cylinder and a sphere (R^2 against r^2).
fn exact_discriminant(q1: &Q, q2: &Q) -> Option<Ordering> {
    let (cyl, sph) = match (q1.kind, q2.kind) {
        (Kind::Cylinder, Kind::Sphere) => (q1, q2),
        (Kind::Sphere, Kind::Cylinder) => (q2, q1),
        _ => return None,
    };
    Some((&sph.r * &sph.r).cmp(&(&cyl.r * &cyl.r)))
}

/// Whether a root is a cone's rational apex lying on the other surface
/// (then the section there is a point).
fn apex_root(q1: &Q, q2: &Q, h: &I) -> bool {
    for (k, other) in [(q1, q2), (q2, q1)] {
        if k.kind != Kind::Cone || !k.r.is_zero() {
            continue;
        }
        let apex = &k.o;
        // On a sphere exactly, or the other cone's rational apex too (a
        // cylinder never passes through a point of its own axis).
        let on_other = match other.kind {
            Kind::Sphere => dot(&sub(apex, &other.o), &sub(apex, &other.o)) == &other.r * &other.r,
            Kind::Cone => other.r.is_zero() && &other.o == apex,
            _ => false,
        };
        if !on_other {
            continue;
        }
        // The root at the apex: h of the apex along the axis surface.
        let axis = if q1.kind == Kind::Sphere { q2 } else { q1 };
        let ha = i(&dot(&sub(apex, &axis.o), &axis.a))
            .div(&i(&dot(&axis.a, &axis.a)).sqrt())
            .expect("a nonzero axis");
        // The exact root is the apex when its enclosure contains the
        // apex's coordinate (the other root is certainly apart).
        if matches!(h.sub(&ha).sign(), None | Some(Ordering::Equal)) {
            return true;
        }
    }
    false
}
