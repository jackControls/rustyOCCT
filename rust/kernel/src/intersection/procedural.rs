//! Procedural intersection curves (D13 of `TOPOLOGY_MODEL.md`, S7b.1 of
//! `REVIEW_NOTES.md`): two cylinders with crossing axes, and a cylinder and
//! a sphere off its axis.
//!
//! The curve is parameterised on the ruled surface (the thinner cylinder;
//! of two equal ones the one whose stored normal and origin come first; the
//! cylinder of a cylinder/sphere pair) by the angle `u` of its ruling in an
//! exactly orthonormal frame: the unit axis `a`, `x` the unit common normal
//! of the axes pointing towards the other axis (towards the sphere's centre),
//! `y = a x x`. A ruling `o + r (cos u x + sin u y) + v a` meets the other
//! quadric where `A v^2 + 2 B(u) v + C(u) = 0`; the curve is
//! `v = (-B +- sqrt(D)) / A` with `D = B^2 - A C`, two branches joined where
//! `D` vanishes. Its class is decided by exact rational predicates (the axes'
//! distance against the radii). In this frame `D` has a closed form: for two
//! cylinders, by Lagrange's identity, `D = A (r2^2 - (r cos u - d)^2)` with
//! `d` the axes' distance; for a sphere whose centre is `e` from the axis,
//! `D = R^2 - e^2 - r^2 + 2 e r cos u`. A loop is therefore `[-t, t]`,
//! `t = arccos c` with `c = (d - r2) / r`, respectively
//! `(e^2 + r^2 - R^2) / (2 e r)`, enclosed by a certified arctangent; rings
//! are `c < -1`, a figure-eight `c = -1` (the node at `u = pi`), a tangent
//! point `c = 1`. Points are evaluated in certified intervals (binary64
//! intervals, rational ones when those cannot decide), so every returned
//! number is an enclosure.
use super::analytic::{Enclosure, Enclosure3};
use crate::certified::{Fast, Interval as I, Real};
use crate::topology::Surface;
use crate::{Error, Result};
use num_bigint::Sign;
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// A branch of the curve: the `+` or the `-` root of the quadratic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Branch {
    Plus,
    Minus,
}

/// A connected piece of a procedural curve.
#[derive(Debug, Clone, PartialEq)]
pub enum Component {
    /// Both branches over `[u0, u1]`, joined at its ends: one closed loop.
    /// `u1 > u0`; either may lie beyond `pi`.
    Loop { u: [Enclosure; 2] },
    /// One branch over a whole turn: a closed ring.
    Ring { branch: Branch },
    /// Both branches over a whole turn, touching at the node (at `u = pi`,
    /// enclosed).
    FigureEight { node: Enclosure },
}

/// D13's procedural curve: two surfaces, the ruled one's parameterisation
/// and the components.
#[derive(Debug, Clone, PartialEq)]
pub struct ProceduralCurve {
    ruled: Surface,
    other: Surface,
    components: Vec<Component>,
}

impl ProceduralCurve {
    /// The surface the curve is parameterised on.
    pub fn ruled(&self) -> &Surface {
        &self.ruled
    }
    pub fn other(&self) -> &Surface {
        &self.other
    }
    pub fn components(&self) -> &[Component] {
        &self.components
    }
    /// A certified enclosure of the curve's points on `branch` at every
    /// parameter in `u = [lo, hi]`; an error where the ruling certainly misses
    /// the other surface.
    pub fn point_at(&self, u: Enclosure, branch: Branch) -> Result<Enclosure3> {
        let [lo, hi] = u;
        let fast = Setup::<Fast>::of(&self.ruled, &self.other)
            .and_then(|s| s.point(&span(lo, hi), branch).map(|p| bounds3(&p)));
        // Near a loop's end D is small and its square root widens binary64
        // intervals: rational intervals then, for an exact parameter.
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
        let s = Setup::<I>::of(&self.ruled, &self.other)?;
        Ok(bounds3(&s.point(&span(lo, hi), branch)?))
    }
}

// ------------------------------------------------------------------ helpers

type X = [R; 3];
type E<T> = [T; 3];

fn q(x: f64) -> R {
    R::from_float(x).expect("finite surface data")
}
fn zero(x: &R) -> bool {
    x.numer().sign() == Sign::NoSign
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
fn scale(a: &X, k: &R) -> X {
    [&a[0] * k, &a[1] * k, &a[2] * k]
}
fn i<T: Real>(x: &R) -> T {
    T::from_r(x)
}
fn e3<T: Real>(v: &X) -> E<T> {
    [i(&v[0]), i(&v[1]), i(&v[2])]
}
fn edot<T: Real>(a: &E<T>, b: &E<T>) -> T {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}
fn ecross<T: Real>(a: &E<T>, b: &E<T>) -> E<T> {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}
fn eadd<T: Real>(a: &E<T>, b: &E<T>) -> E<T> {
    [a[0].add(&b[0]), a[1].add(&b[1]), a[2].add(&b[2])]
}
fn esub<T: Real>(a: &E<T>, b: &E<T>) -> E<T> {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}
fn escale<T: Real>(a: &E<T>, k: &T) -> E<T> {
    [a[0].mul(k), a[1].mul(k), a[2].mul(k)]
}
fn limit(what: &'static str) -> Error {
    Error::ComputationLimit(what)
}
fn unit<T: Real>(v: &X) -> Result<E<T>> {
    let n = i::<T>(&dot(v, v)).sqrt();
    let d = |x: &R| i::<T>(x).div(&n).ok_or(limit("a unit direction"));
    Ok([d(&v[0])?, d(&v[1])?, d(&v[2])?])
}
fn bounds<T: Real>(x: &T) -> Enclosure {
    let (lo, hi) = x.bounds_f64();
    [lo, hi]
}
fn bounds3<T: Real>(v: &E<T>) -> Enclosure3 {
    [bounds(&v[0]), bounds(&v[1]), bounds(&v[2])]
}
/// The enclosure `[lo, hi]` (binary64 ends) in a tier.
fn span<T: Real>(lo: f64, hi: f64) -> T {
    if lo == hi {
        T::exact_f64(lo)
    } else {
        T::exact_f64(lo).union(&T::exact_f64(hi))
    }
}
/// Pi in a tier.
fn pi<T: Real>() -> T {
    let p = crate::certified::pi();
    T::from_r(p.lo()).union(&T::from_r(p.hi()))
}

/// A cylinder's or sphere's exact data.
struct Quad {
    o: X,
    a: X,
    r: R,
    cylinder: bool,
}

impl Quad {
    fn of(s: &Surface) -> Option<Self> {
        let (f, r, cylinder) = match s {
            Surface::Cylinder { frame, radius } => (frame, radius, true),
            Surface::Sphere { frame, radius } => (frame, radius, false),
            _ => return None,
        };
        Some(Self {
            o: f.origin().to_array().map(q),
            a: f.normal().to_array().map(q),
            r: q(*r),
            cylinder,
        })
    }
    /// Whether `a` comes first by its stored normal, then origin (a tie rule
    /// independent of argument order).
    fn first(a: &Surface, b: &Surface) -> bool {
        let key = |s: &Surface| {
            let f = match s {
                Surface::Cylinder { frame, .. } | Surface::Sphere { frame, .. } => *frame,
                _ => unreachable!("cylinders and spheres only"),
            };
            let [x, y, z] = f.normal().to_array();
            let [u, v, w] = f.origin().to_array();
            [x, y, z, u, v, w]
        };
        let (ka, kb) = (key(a), key(b));
        ka.iter()
            .zip(&kb)
            .map(|(x, y)| x.total_cmp(y))
            .find(|o| *o != Ordering::Equal)
            .unwrap_or(Ordering::Equal)
            != Ordering::Greater
    }
}

/// The parameterisation in a certified tier: the ruled frame and the other
/// quadric.
struct Setup<T> {
    o: E<T>,
    a: E<T>,
    x: E<T>,
    y: E<T>,
    r: T,
    /// The other quadric: its origin, its axis and squared length (a
    /// cylinder) and its radius squared.
    oo: E<T>,
    oa: Option<(E<T>, T)>,
    rr: T,
    big_a: T,
}

impl<T: Real> Setup<T> {
    fn of(ruled: &Surface, other: &Surface) -> Result<Self> {
        let (Some(p), Some(o)) = (Quad::of(ruled), Quad::of(other)) else {
            return Err(Error::OutOfDomain("procedural curves of these surfaces"));
        };
        let toward = if o.cylinder {
            let m = cross(&p.a, &o.a);
            if dot(&sub(&o.o, &p.o), &m).numer().sign() == Sign::Minus {
                scale(&m, &R::from_integer((-1).into()))
            } else {
                m
            }
        } else {
            let w = sub(&o.o, &p.o);
            sub(&w, &scale(&p.a, &(dot(&w, &p.a) / dot(&p.a, &p.a))))
        };
        let a = unit::<T>(&p.a)?;
        let x = unit::<T>(&toward)?;
        let y = ecross(&a, &x);
        let big_a = if o.cylinder {
            // 1 - (a . a2)^2 / (|a|^2 |a2|^2), exact.
            let k = dot(&p.a, &o.a);
            i(&(R::from_integer(1.into()) - &k * &k / (dot(&p.a, &p.a) * dot(&o.a, &o.a))))
        } else {
            T::exact_f64(1.0)
        };
        Ok(Self {
            o: e3(&p.o),
            a,
            x,
            y,
            r: i(&p.r),
            oo: e3(&o.o),
            oa: o.cylinder.then(|| (e3(&o.a), i(&dot(&o.a, &o.a)))),
            rr: i(&(&o.r * &o.r)),
            big_a,
        })
    }
    /// The other quadric's bilinear form: `M = 1 - a2 a2^T / |a2|^2` for a
    /// cylinder, the identity for a sphere.
    fn form(&self, p: &E<T>, q: &E<T>) -> Result<T> {
        let base = edot(p, q);
        Ok(match &self.oa {
            Some((a2, aa)) => {
                let k = edot(p, a2).mul(&edot(q, a2)).div(aa);
                base.sub(&k.ok_or(limit("an axis projector"))?)
            }
            None => base,
        })
    }
    /// The ruling's foot at `u` and its derivative in `u`.
    fn ruling(&self, u: &T) -> (E<T>, E<T>) {
        let (c, s) = T::cos_sin(u);
        let radial = eadd(&escale(&self.x, &c), &escale(&self.y, &s));
        let derivative = esub(&escale(&self.y, &c), &escale(&self.x, &s));
        (
            eadd(&self.o, &escale(&radial, &self.r)),
            escale(&derivative, &self.r),
        )
    }
    fn point(&self, u: &T, branch: Branch) -> Result<E<T>> {
        let (p, _) = self.ruling(u);
        let rel = esub(&p, &self.oo);
        let b = self.form(&rel, &self.a)?;
        let c = self.form(&rel, &rel)?.sub(&self.rr);
        let d = b.square().sub(&self.big_a.mul(&c));
        if d.sign() == Some(Ordering::Less) {
            return Err(Error::OutOfDomain("a parameter off the curve"));
        }
        // sqrt clamps a negative part (the enclosure of a root) to zero.
        let root = d.sqrt();
        let root = match branch {
            Branch::Plus => root,
            Branch::Minus => root.neg(),
        };
        let v = b
            .neg()
            .add(&root)
            .div(&self.big_a)
            .ok_or(limit("procedural point"))?;
        Ok(eadd(&p, &escale(&self.a, &v)))
    }
}

// ------------------------------------------------------------------ classes

enum Class {
    Empty,
    Point,
    Loop,
    FigureEight,
    Rings,
}

fn classify_cylinders(p: &Quad, o: &Quad) -> Class {
    let m = cross(&p.a, &o.a);
    let d = dot(&sub(&o.o, &p.o), &m);
    let d2 = &d * &d / dot(&m, &m);
    let sum = &p.r + &o.r;
    let diff = &p.r - &o.r;
    let (s2, f2) = (&sum * &sum, &diff * &diff);
    match d2.cmp(&s2) {
        Ordering::Greater => Class::Empty,
        Ordering::Equal => Class::Point,
        Ordering::Less => match d2.cmp(&f2) {
            Ordering::Equal if !zero(&diff) => Class::FigureEight,
            Ordering::Less => Class::Rings,
            _ => Class::Loop,
        },
    }
}

fn classify_sphere(c: &Quad, s: &Quad) -> Class {
    let w = sub(&s.o, &c.o);
    let w = sub(&w, &scale(&c.a, &(dot(&w, &c.a) / dot(&c.a, &c.a))));
    let e2 = dot(&w, &w);
    let (r, big_r) = (&c.r, &s.r);
    // Sign of R - |e - r|: R^2 against e^2 - 2 e r + r^2, i.e. 2 e r against
    // e^2 + r^2 - R^2.
    let rhs = &e2 + r * r - big_r * big_r;
    let low = if rhs.numer().sign() == Sign::Minus {
        Ordering::Greater
    } else {
        let four = R::from_integer(4.into());
        (&four * &e2 * r * r).cmp(&(&rhs * &rhs))
    };
    match low {
        Ordering::Less => return Class::Empty,
        Ordering::Equal => return Class::Point,
        Ordering::Greater => {}
    }
    // Sign of R - (e + r): (R - r) against e.
    let t = big_r - r;
    let high = if t.numer().sign() == Sign::Minus {
        Ordering::Less
    } else {
        (&t * &t).cmp(&e2)
    };
    match high {
        Ordering::Less => Class::Loop,
        Ordering::Equal => Class::FigureEight,
        Ordering::Greater => Class::Rings,
    }
}

// ------------------------------------------------------------------ loops

/// The half-range `t = arccos c` of a loop `[-t, t]`: `c = (d - r2) / r` for
/// two cylinders (`d` the axes' distance), `c = (e^2 + r^2 - R^2) / (2 e r)`
/// for a cylinder and a sphere (`e` the centre's distance from the axis),
/// both in `(-1, 1)` for a loop.
fn half_range(ruled: &Quad, other: &Quad) -> Result<Enclosure> {
    let c = if other.cylinder {
        let m = cross(&ruled.a, &other.a);
        let d = dot(&sub(&other.o, &ruled.o), &m);
        let distance = I::exact(&d * &d / dot(&m, &m)).sqrt();
        distance
            .sub(&I::exact(other.r.clone()))
            .div(&I::exact(ruled.r.clone()))
    } else {
        let w = sub(&other.o, &ruled.o);
        let w = sub(
            &w,
            &scale(&ruled.a, &(dot(&w, &ruled.a) / dot(&ruled.a, &ruled.a))),
        );
        let e2 = dot(&w, &w);
        let top = I::exact(&e2 + &ruled.r * &ruled.r - &other.r * &other.r);
        let bottom = I::exact(&ruled.r * R::from_integer(2.into())).mul(&I::exact(e2).sqrt());
        top.div(&bottom)
    }
    .ok_or(limit("a loop's cosine"))?;
    let sine = I::exact(R::from_integer(1.into())).sub(&c.square()).sqrt();
    let t = crate::certified::atan2(&sine, &c).ok_or(limit("a loop's half-range"))?;
    Ok(bounds(&t))
}

// ------------------------------------------------------------------ entry

/// The outcome for a cylinder pair with crossing axes (not equal radii
/// through a common point, S7a's ellipses) or a cylinder and a sphere off its
/// axis.
pub(crate) enum Found {
    Empty,
    Point(Enclosure3),
    Curve(Box<ProceduralCurve>),
}

/// `None` for pairs that are not two cylinders with crossing axes or a
/// cylinder and a sphere off its axis.
pub(crate) fn intersect(a: &Surface, b: &Surface) -> Result<Option<Found>> {
    let (Some(qa), Some(qb)) = (Quad::of(a), Quad::of(b)) else {
        return Ok(None);
    };
    let (ruled, other, class) = match (qa.cylinder, qb.cylinder) {
        (true, true) => {
            if cross(&qa.a, &qb.a).iter().all(zero) {
                return Ok(None);
            }
            // The thinner cylinder, else the one first by its stored data.
            let first = match qa.r.cmp(&qb.r) {
                Ordering::Less => true,
                Ordering::Greater => false,
                Ordering::Equal => Quad::first(a, b),
            };
            let (p, o, ps, os) = if first {
                (&qa, &qb, a, b)
            } else {
                (&qb, &qa, b, a)
            };
            (ps, os, classify_cylinders(p, o))
        }
        (true, false) | (false, true) => {
            let (c, s, cs, ss) = if qa.cylinder {
                (&qa, &qb, a, b)
            } else {
                (&qb, &qa, b, a)
            };
            let w = sub(&s.o, &c.o);
            let w = sub(&w, &scale(&c.a, &(dot(&w, &c.a) / dot(&c.a, &c.a))));
            if w.iter().all(zero) {
                return Ok(None);
            }
            (cs, ss, classify_sphere(c, s))
        }
        (false, false) => return Ok(None),
    };
    let fast = Setup::<Fast>::of(ruled, other)?;
    let exact = || Setup::<I>::of(ruled, other);
    let curve = |components| {
        Box::new(ProceduralCurve {
            ruled: ruled.clone(),
            other: other.clone(),
            components,
        })
    };
    Ok(Some(match class {
        Class::Empty => Found::Empty,
        // Tangent at u = 0, towards the other surface.
        Class::Point => Found::Point(match fast.point(&Fast::exact_f64(0.0), Branch::Plus) {
            Err(Error::ComputationLimit(_)) => {
                bounds3(&exact()?.point(&I::exact_f64(0.0), Branch::Plus)?)
            }
            other => bounds3(&other?),
        }),
        Class::Rings => Found::Curve(curve(vec![
            Component::Ring {
                branch: Branch::Plus,
            },
            Component::Ring {
                branch: Branch::Minus,
            },
        ])),
        // The node is on the far side, u = pi: the thinner cylinder touches
        // the thicker from inside away from its axis, Viviani's sphere away
        // from its centre.
        Class::FigureEight => Found::Curve(curve(vec![Component::FigureEight {
            node: bounds(&pi::<I>()),
        }])),
        Class::Loop => {
            let (p, o) = (Quad::of(ruled), Quad::of(other));
            let t = half_range(&p.expect("a cylinder"), &o.expect("a quadric"))?;
            Found::Curve(curve(vec![Component::Loop {
                u: [[-t[1], -t[0]], t],
            }]))
        }
    }))
}
