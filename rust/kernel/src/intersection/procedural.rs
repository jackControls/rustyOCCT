//! Procedural intersection curves (D13 of `TOPOLOGY_MODEL.md`, S7b of
//! `REVIEW_NOTES.md`): two cylinders with crossing axes, a cylinder and a
//! sphere off its axis (S7b.1), a cylinder and a cone whose axes are not
//! coaxial, and a sphere and a cone with the centre off the axis (S7b.2).
//!
//! The curve is parameterised on a ruled surface by the angle `u` of its
//! ruling in an exactly orthonormal frame (unit axis `a`, `x` towards the
//! other surface, `y = a x x`): a cylinder's ruling `P(u) + v a`,
//! `P(u) = o + r (cos u x + sin u y)`; a cone's ruling through its apex
//! `V + v d(u)`, `d(u) = cos h a + sin h (cos u x + sin u y)` (`h` the
//! half-angle). The other quadric is `f(p) = Q(p - q, p - q) - k` (a sphere:
//! `Q` the identity, `k = R^2`; a cylinder: `Q = 1 - a2 a2^T / |a2|^2`,
//! `k = r2^2`; a cone: `q` its apex, `Q = cos^2 h2 - a2 a2^T / |a2|^2`,
//! `k = 0`), so the ruling meets it where `A v^2 + 2 B v + C = 0` with
//! `A = Q(dir, dir)`, `B = Q(foot - q, dir)`, `C = f(foot)`, and the curve is
//! `v = (-B +- sqrt(D)) / A`, `D = B^2 - A C`, two branches joined where `D`
//! vanishes.
//!
//! * Two cylinders (on the thinner; of equal ones the first by stored
//!   normal, then origin) and a cylinder and a sphere (on the cylinder):
//!   exact classes, and `D` in closed form (Lagrange's identity; the
//!   centre's distance), so a loop is `[-t, t]`, `t = arccos c`.
//! * A sphere and a cone (on the cone, `A = 1`): `B(u) = b0 + b1 cos u`, `C`
//!   constant; the apex inside the sphere gives two rings, outside it loops
//!   where `B > sqrt(C)` or `B < -sqrt(C)`, each end an `arccos`.
//! * A cylinder and a cone (on the cylinder, `A` constant): the roots of `D`
//!   are isolated by certified subdivision with a mean-value enclosure and a
//!   work budget; loops lie between consecutive roots where `D > 0`.
//!
//! Points are evaluated in certified intervals (binary64 intervals, rational
//! ones when those cannot decide), so every returned number is an enclosure.
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
    /// `u1 > u0`; either may lie beyond `pi`. The ends' enclosures are
    /// disjoint, `u[0][1] < u[1][0]`, so both branches are defined at every
    /// parameter between them (a loop narrower than that is
    /// `ComputationLimit`).
    Loop { u: [Enclosure; 2] },
    /// One branch over a whole turn: a closed ring.
    Ring { branch: Branch },
    /// Both branches over a whole turn, touching at the node (enclosed: at
    /// `u = pi` for the ruled surfaces, `0` or `pi` on a torus).
    FigureEight { node: Enclosure },
}

/// D13's procedural curve: two surfaces, the parameterisation on the first
/// (a ruled surface's rulings, a torus's meridians) and the components.
#[derive(Debug, Clone, PartialEq)]
pub struct ProceduralCurve {
    carrier: Surface,
    other: Surface,
    components: Vec<Component>,
}

impl ProceduralCurve {
    /// The curve, or `ComputationLimit` for a loop whose ends' enclosures
    /// meet or cross: a loop within an ulp of its parameter (a plane within
    /// rounding of a torus's axis, turned by an ulp on another host's
    /// `hypot`) contains no binary64 parameter to evaluate it at.
    pub(super) fn new(
        carrier: Surface,
        other: Surface,
        components: Vec<Component>,
    ) -> Result<Self> {
        for c in &components {
            if let Component::Loop { u } = c {
                if u[0][1].partial_cmp(&u[1][0]) != Some(Ordering::Less) {
                    return Err(limit("a loop narrower than its ends' enclosures"));
                }
            }
        }
        Ok(Self {
            carrier,
            other,
            components,
        })
    }
    /// The surface the curve is parameterised on: a cylinder or cone by the
    /// angle of its ruling, a torus by the angle of its meridian.
    pub fn carrier(&self) -> &Surface {
        &self.carrier
    }
    pub fn other(&self) -> &Surface {
        &self.other
    }
    pub fn components(&self) -> &[Component] {
        &self.components
    }
    /// A certified enclosure of the curve's points on `branch` at every
    /// parameter in `u = [lo, hi]`; an error where the ruling or meridian
    /// certainly misses the other surface.
    pub fn point_at(&self, u: Enclosure, branch: Branch) -> Result<Enclosure3> {
        if matches!(self.carrier, Surface::Torus { .. }) {
            return super::toroidal::point_at(&self.carrier, &self.other, u, branch);
        }
        let [lo, hi] = u;
        let fast = Setup::<Fast>::of(&self.carrier, &self.other)
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
        let s = Setup::<I>::of(&self.carrier, &self.other)?;
        Ok(bounds3(&s.point(&span(lo, hi), branch)?))
    }
}

// ------------------------------------------------------------------ helpers

pub(super) type X = [R; 3];
pub(super) type E<T> = [T; 3];

pub(super) fn q(x: f64) -> R {
    R::from_float(x).expect("finite surface data")
}
pub(super) fn zero(x: &R) -> bool {
    x.numer().sign() == Sign::NoSign
}
pub(super) fn dot(a: &X, b: &X) -> R {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}
pub(super) fn cross(a: &X, b: &X) -> X {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
pub(super) fn sub(a: &X, b: &X) -> X {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}
pub(super) fn scale(a: &X, k: &R) -> X {
    [&a[0] * k, &a[1] * k, &a[2] * k]
}
/// The component of `w` perpendicular to `a`.
pub(super) fn perpendicular(w: &X, a: &X) -> X {
    sub(w, &scale(a, &(dot(w, a) / dot(a, a))))
}
pub(super) fn i<T: Real>(x: &R) -> T {
    T::from_r(x)
}
pub(super) fn e3<T: Real>(v: &X) -> E<T> {
    [i(&v[0]), i(&v[1]), i(&v[2])]
}
pub(super) fn edot<T: Real>(a: &E<T>, b: &E<T>) -> T {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}
pub(super) fn ecross<T: Real>(a: &E<T>, b: &E<T>) -> E<T> {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}
pub(super) fn eadd<T: Real>(a: &E<T>, b: &E<T>) -> E<T> {
    [a[0].add(&b[0]), a[1].add(&b[1]), a[2].add(&b[2])]
}
pub(super) fn esub<T: Real>(a: &E<T>, b: &E<T>) -> E<T> {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}
pub(super) fn escale<T: Real>(a: &E<T>, k: &T) -> E<T> {
    [a[0].mul(k), a[1].mul(k), a[2].mul(k)]
}
pub(super) fn limit(what: &'static str) -> Error {
    Error::ComputationLimit(what)
}
pub(super) fn unit<T: Real>(v: &X) -> Result<E<T>> {
    let n = i::<T>(&dot(v, v)).sqrt();
    let d = |x: &R| i::<T>(x).div(&n).ok_or(limit("a unit direction"));
    Ok([d(&v[0])?, d(&v[1])?, d(&v[2])?])
}
pub(super) fn bounds<T: Real>(x: &T) -> Enclosure {
    let (lo, hi) = x.bounds_f64();
    [lo, hi]
}
pub(super) fn bounds3<T: Real>(v: &E<T>) -> Enclosure3 {
    [bounds(&v[0]), bounds(&v[1]), bounds(&v[2])]
}
/// The enclosure `[lo, hi]` (binary64 ends) in a tier.
pub(super) fn span<T: Real>(lo: f64, hi: f64) -> T {
    if lo == hi {
        T::exact_f64(lo)
    } else {
        T::exact_f64(lo).union(&T::exact_f64(hi))
    }
}
/// Pi in a tier.
pub(super) fn pi<T: Real>() -> T {
    let p = crate::certified::pi();
    T::from_r(p.lo()).union(&T::from_r(p.hi()))
}
pub(super) fn positive<T: Real>(x: &T) -> Result<bool> {
    match x.sign() {
        Some(Ordering::Greater) => Ok(true),
        Some(Ordering::Less) => Ok(false),
        _ => Err(limit("a certified sign")),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Cylinder,
    Cone,
    Sphere,
}

/// A cylinder's, cone's or sphere's exact data.
struct Quad {
    kind: Kind,
    o: X,
    a: X,
    r: R,
    angle: f64,
}

impl Quad {
    fn of(s: &Surface) -> Option<Self> {
        let (kind, f, r, angle) = match s {
            Surface::Cylinder { frame, radius } => (Kind::Cylinder, frame, radius, 0.0),
            Surface::Cone {
                frame,
                radius,
                half_angle,
            } => (Kind::Cone, frame, radius, *half_angle),
            Surface::Sphere { frame, radius } => (Kind::Sphere, frame, radius, 0.0),
            _ => return None,
        };
        Some(Self {
            kind,
            o: f.origin().to_array().map(q),
            a: f.normal().to_array().map(q),
            r: q(*r),
            angle,
        })
    }
    /// Whether `a` comes first by its stored normal, then origin (a tie rule
    /// independent of argument order).
    fn first(a: &Surface, b: &Surface) -> bool {
        let key = |s: &Surface| {
            let f = match s {
                Surface::Cylinder { frame, .. }
                | Surface::Cone { frame, .. }
                | Surface::Sphere { frame, .. } => *frame,
                _ => unreachable!("quadrics only"),
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
    /// (cos, sin) of a cone's half-angle in a tier.
    fn cos_sin<T: Real>(&self) -> (T, T) {
        T::cos_sin(&T::exact_f64(self.angle))
    }
    /// A cone's apex `o - a r / (|a| tan h)` in a tier.
    fn apex<T: Real>(&self) -> Result<E<T>> {
        if zero(&self.r) {
            return Ok(e3(&self.o));
        }
        let (c, s) = self.cos_sin::<T>();
        let back = i::<T>(&self.r)
            .mul(&c)
            .div(&s)
            .ok_or(limit("a cone's apex"))?;
        Ok(esub(&e3(&self.o), &escale(&unit::<T>(&self.a)?, &back)))
    }
}

/// The ruled surface's frame in a tier.
enum Ruled<T> {
    /// `P(u) + v a`.
    Cylinder { o: E<T>, r: T },
    /// `V + v d(u)`, `d(u) = cos h a + sin h (cos u x + sin u y)`.
    Cone { apex: E<T>, c: T, s: T },
}

/// The parameterisation in a certified tier.
struct Setup<T> {
    ruled: Ruled<T>,
    a: E<T>,
    x: E<T>,
    y: E<T>,
    /// The other quadric: `Q(p - q, p - q) - k`, `Q(u, w) = c2 (u . w) -
    /// (u . a2)(w . a2) / |a2|^2`.
    q: E<T>,
    axis: Option<(E<T>, T)>,
    c2: T,
    k: T,
}

impl<T: Real> Setup<T> {
    fn of(ruled: &Surface, other: &Surface) -> Result<Self> {
        let (Some(p), Some(o)) = (Quad::of(ruled), Quad::of(other)) else {
            return Err(Error::OutOfDomain("procedural curves of these surfaces"));
        };
        // x: towards the other surface, along the common normal of the axes,
        // or perpendicular to the ruled axis for parallel axes and spheres.
        let toward = match o.kind {
            Kind::Sphere => perpendicular(&sub(&o.o, &p.o), &p.a),
            _ => {
                let m = cross(&p.a, &o.a);
                if m.iter().all(zero) {
                    perpendicular(&sub(&o.o, &p.o), &p.a)
                } else if dot(&sub(&o.o, &p.o), &m).numer().sign() == Sign::Minus {
                    scale(&m, &R::from_integer((-1).into()))
                } else {
                    m
                }
            }
        };
        let a = unit::<T>(&p.a)?;
        let x = unit::<T>(&toward)?;
        let y = ecross(&a, &x);
        let ruled = match p.kind {
            Kind::Cylinder => Ruled::Cylinder {
                o: e3(&p.o),
                r: i(&p.r),
            },
            Kind::Cone => {
                let (c, s) = p.cos_sin::<T>();
                Ruled::Cone {
                    apex: p.apex()?,
                    c,
                    s,
                }
            }
            Kind::Sphere => return Err(Error::OutOfDomain("a sphere is not ruled")),
        };
        let one = T::exact_f64(1.0);
        let (q, axis, c2, k) = match o.kind {
            Kind::Sphere => (e3(&o.o), None, one, i(&(&o.r * &o.r))),
            Kind::Cylinder => (
                e3(&o.o),
                Some((e3(&o.a), i(&dot(&o.a, &o.a)))),
                one,
                i(&(&o.r * &o.r)),
            ),
            Kind::Cone => {
                let (c, _) = o.cos_sin::<T>();
                (
                    o.apex()?,
                    Some((e3(&o.a), i(&dot(&o.a, &o.a)))),
                    c.square(),
                    T::exact_f64(0.0),
                )
            }
        };
        Ok(Self {
            ruled,
            a,
            x,
            y,
            q,
            axis,
            c2,
            k,
        })
    }
    fn form(&self, p: &E<T>, w: &E<T>) -> Result<T> {
        let base = edot(p, w).mul(&self.c2);
        Ok(match &self.axis {
            Some((a2, aa)) => {
                let k = edot(p, a2).mul(&edot(w, a2)).div(aa);
                base.sub(&k.ok_or(limit("an axis projector"))?)
            }
            None => base,
        })
    }
    /// The ruling at `u`: its foot, its direction and their derivatives in
    /// `u`.
    fn ruling(&self, u: &T) -> (E<T>, E<T>, E<T>, E<T>) {
        let (c, s) = T::cos_sin(u);
        let radial = eadd(&escale(&self.x, &c), &escale(&self.y, &s));
        let turn = esub(&escale(&self.y, &c), &escale(&self.x, &s));
        let zero = T::exact_f64(0.0);
        let still = [zero.clone(), zero.clone(), zero];
        match &self.ruled {
            Ruled::Cylinder { o, r } => (
                eadd(o, &escale(&radial, r)),
                self.a.clone(),
                escale(&turn, r),
                still,
            ),
            Ruled::Cone { apex, c: ch, s: sh } => (
                apex.clone(),
                eadd(&escale(&self.a, ch), &escale(&radial, sh)),
                still,
                escale(&turn, sh),
            ),
        }
    }
    /// (A, B, C) of the ruling at `u`.
    fn coefficients(&self, u: &T) -> Result<(T, T, T)> {
        let (foot, dir, _, _) = self.ruling(u);
        let rel = esub(&foot, &self.q);
        Ok((
            self.form(&dir, &dir)?,
            self.form(&rel, &dir)?,
            self.form(&rel, &rel)?.sub(&self.k),
        ))
    }
    fn discriminant(&self, u: &T) -> Result<T> {
        let (a, b, c) = self.coefficients(u)?;
        Ok(b.square().sub(&a.mul(&c)))
    }
    /// dD/du for a cylinder's ruling (A constant): 2 B B' - A C', with
    /// B' = Q(P', a) and C' = 2 Q(P', P - q).
    fn derivative(&self, u: &T) -> Result<T> {
        let (foot, dir, dfoot, _) = self.ruling(u);
        let rel = esub(&foot, &self.q);
        let (a, b) = (self.form(&dir, &dir)?, self.form(&rel, &dir)?);
        let db = self.form(&dfoot, &dir)?;
        let two = T::exact_f64(2.0);
        let dc = self.form(&dfoot, &rel)?.mul(&two);
        Ok(b.mul(&db).mul(&two).sub(&a.mul(&dc)))
    }
    fn point(&self, u: &T, branch: Branch) -> Result<E<T>> {
        let (a, b, c) = self.coefficients(u)?;
        let d = b.square().sub(&a.mul(&c));
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
            .div(&a)
            .ok_or(limit("procedural point"))?;
        let (foot, dir, _, _) = self.ruling(u);
        Ok(eadd(&foot, &escale(&dir, &v)))
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
    let w = perpendicular(&sub(&s.o, &c.o), &c.a);
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

/// `arccos c` for `c` in `(-1, 1)`, enclosed.
pub(super) fn arccos(c: &I) -> Result<Enclosure> {
    let sine = I::exact_f64(1.0).sub(&c.square()).sqrt();
    let t = crate::certified::atan2(&sine, c).ok_or(limit("a loop's half-range"))?;
    Ok(bounds(&t))
}

/// The half-range `t = arccos c` of a loop `[-t, t]` of two cylinders
/// (`c = (d - r2) / r`, `d` the axes' distance) or a cylinder and a sphere
/// (`c = (e^2 + r^2 - R^2) / (2 e r)`, `e` the centre's distance).
fn half_range(ruled: &Quad, other: &Quad) -> Result<Enclosure> {
    let c = if other.kind == Kind::Cylinder {
        let m = cross(&ruled.a, &other.a);
        let d = dot(&sub(&other.o, &ruled.o), &m);
        let distance = I::exact(&d * &d / dot(&m, &m)).sqrt();
        distance
            .sub(&I::exact(other.r.clone()))
            .div(&I::exact(ruled.r.clone()))
    } else {
        let w = perpendicular(&sub(&other.o, &ruled.o), &ruled.a);
        let e2 = dot(&w, &w);
        let top = I::exact(&e2 + &ruled.r * &ruled.r - &other.r * &other.r);
        let bottom = I::exact(&ruled.r * R::from_integer(2.into())).mul(&I::exact(e2).sqrt());
        top.div(&bottom)
    }
    .ok_or(limit("a loop's cosine"))?;
    arccos(&c)
}

/// Where `cos u > k` (`above`) or `cos u < k`: nothing, the whole turn, or
/// a loop around 0 (`[-t, t]`) or around pi (`[t, 2 pi - t]`).
enum Region {
    None,
    All,
    Loop([Enclosure; 2]),
}

fn region(k: &I, above: bool) -> Result<Region> {
    let one = I::exact_f64(1.0);
    let beyond = |x: &I, y: &I| positive(&x.sub(y));
    if above {
        if !beyond(&one, k)? {
            return Ok(Region::None);
        }
        if beyond(&one.neg(), k)? {
            return Ok(Region::All);
        }
        let t = arccos(k)?;
        Ok(Region::Loop([[-t[1], -t[0]], t]))
    } else {
        if !beyond(k, &one.neg())? {
            return Ok(Region::None);
        }
        if beyond(k, &one)? {
            return Ok(Region::All);
        }
        let t = arccos(k)?;
        let tau = pi::<I>().mul(&I::exact_f64(2.0));
        let far = [
            bounds(&tau.sub(&I::exact_f64(t[1]))),
            bounds(&tau.sub(&I::exact_f64(t[0]))),
        ];
        Ok(Region::Loop([t, [far[0][0], far[1][1]]]))
    }
}

/// A sphere seen from a cone's apex: `B(u) = b0 + b1 cos u`, `C` constant.
/// `None` when the apex is on the sphere.
fn cone_sphere(cone: &Quad, sphere: &Quad) -> Result<Option<Vec<Component>>> {
    let apex = cone.apex::<I>()?;
    let rel = esub(&apex, &e3(&sphere.o));
    // The apex on the sphere exactly: a rational apex only.
    if zero(&cone.r) {
        let w = sub(&cone.o, &sphere.o);
        if dot(&w, &w) == &sphere.r * &sphere.r {
            return Ok(None);
        }
    }
    let c = edot(&rel, &rel).sub(&I::exact(&sphere.r * &sphere.r));
    let rings = || {
        vec![
            Component::Ring {
                branch: Branch::Plus,
            },
            Component::Ring {
                branch: Branch::Minus,
            },
        ]
    };
    if !positive(&c)? {
        return Ok(Some(rings()));
    }
    // b0 = cos h (V - c) . a, b1 = sin h (V - c) . x.
    let (ch, sh) = cone.cos_sin::<I>();
    let a = unit::<I>(&cone.a)?;
    let x = unit::<I>(&perpendicular(&sub(&sphere.o, &cone.o), &cone.a))?;
    let b0 = ch.mul(&edot(&rel, &a));
    let b1 = sh.mul(&edot(&rel, &x));
    let s = c.sqrt();
    // B > s and B < -s as regions of cos u.
    let mut out = Vec::new();
    let mut all = false;
    for (target, above_if_positive) in [(s.clone(), true), (s.neg(), false)] {
        let k = target.sub(&b0).div(&b1).ok_or(limit("a cone's sinusoid"))?;
        // b1 cos u > t - b0: cos u > k when b1 > 0, < k when b1 < 0; the
        // region below -s reverses.
        let above = positive(&b1)? == above_if_positive;
        match region(&k, above)? {
            Region::None => {}
            Region::All => all = true,
            Region::Loop(u) => out.push(Component::Loop { u }),
        }
    }
    if all {
        return Ok(Some(rings()));
    }
    out.sort_by(|p, q| match (p, q) {
        (Component::Loop { u: a }, Component::Loop { u: b }) => a[0][0].total_cmp(&b[0][0]),
        _ => Ordering::Equal,
    });
    Ok(Some(out))
}

// ------------------------------------------------------------------ roots

/// Certified simple roots of `D` over `[lo, hi]` (binary64 ends), by
/// subdivision: a piece is dropped when `D` has a certain sign on it (its
/// mean-value enclosure `D(m) + D'(piece) (piece - m)`, tight where `D` dips
/// towards zero), kept as a root when `D` changes sign between its ends and
/// `D'` has a certain sign on it, and halved otherwise (up to `depth`
/// halvings); `budget` bounds the pieces a tier examines. Iterative, with an
/// explicit stack of pieces in increasing order: allocation stacks do not
/// grow with the depth (AddressSanitizer records each distinct one).
fn isolate<T: Real>(
    s: &Setup<T>,
    lo: f64,
    hi: f64,
    depth: usize,
    budget: &mut usize,
    out: &mut Vec<Enclosure>,
) -> Result<()> {
    let mut pending = vec![(lo, hi, depth)];
    while let Some((lo, hi, depth)) = pending.pop() {
        *budget = budget
            .checked_sub(1)
            .ok_or(limit("isolating a loop's ends"))?;
        let piece = span::<T>(lo, hi);
        let mid = 0.5 * lo + 0.5 * hi;
        let slope = s.derivative(&piece)?;
        let enclosure = s
            .discriminant(&T::exact_f64(mid))?
            .add(&slope.mul(&T::exact_f64(lo - mid).union(&T::exact_f64(hi - mid))));
        if matches!(enclosure.sign(), Some(Ordering::Less | Ordering::Greater)) {
            continue;
        }
        let (a, b) = (
            s.discriminant(&T::exact_f64(lo))?.sign(),
            s.discriminant(&T::exact_f64(hi))?.sign(),
        );
        let change = matches!(
            (a, b),
            (Some(Ordering::Less), Some(Ordering::Greater))
                | (Some(Ordering::Greater), Some(Ordering::Less))
        );
        if change && matches!(slope.sign(), Some(Ordering::Less | Ordering::Greater)) {
            out.push(refine(s, lo, hi, a));
            continue;
        }
        if depth == 0 || !(lo < mid && mid < hi) {
            return Err(limit("isolating a loop's ends"));
        }
        // The upper half is pushed first so the lower is examined first.
        pending.push((mid, hi, depth - 1));
        pending.push((lo, mid, depth - 1));
    }
    Ok(())
}

/// A certified simple root's enclosure `[lo, hi]` narrowed by bisection on
/// the certain sign of `D` at the middle, until the sign is undecided or
/// the ends are adjacent binary64 values.
fn refine<T: Real>(s: &Setup<T>, mut lo: f64, mut hi: f64, at_lo: Option<Ordering>) -> Enclosure {
    for _ in 0..80 {
        let mid = 0.5 * lo + 0.5 * hi;
        if !(lo < mid && mid < hi) {
            break;
        }
        match s
            .discriminant(&T::exact_f64(mid))
            .ok()
            .and_then(|d| d.sign())
        {
            Some(Ordering::Equal) => return [mid, mid],
            Some(sign) if Some(sign) == at_lo => lo = mid,
            Some(_) => hi = mid,
            None => break,
        }
    }
    [lo, hi]
}

/// The roots of `D` over one turn, `u` in `[-pi, pi]` (binary64 bounds of
/// pi), each once: each of 64 pieces in binary64 intervals, a piece they
/// cannot settle again in rational intervals.
fn roots(fast: &Setup<Fast>, exact: &Setup<I>) -> Result<Vec<Enclosure>> {
    // The binary64 value just above pi (next_up needs Rust 1.86).
    let pi_hi = f64::from_bits(std::f64::consts::PI.to_bits() + 1);
    let lo = -pi_hi;
    let pieces = 64;
    let (mut fast_budget, mut exact_budget) = (20_000, 4_000);
    let mut out = Vec::new();
    for k in 0..pieces {
        let a = lo + (2.0 * pi_hi) * (k as f64) / (pieces as f64);
        let b = if k + 1 == pieces {
            pi_hi
        } else {
            lo + (2.0 * pi_hi) * ((k + 1) as f64) / (pieces as f64)
        };
        let mut found = Vec::new();
        match isolate(fast, a, b, 60, &mut fast_budget, &mut found) {
            Err(Error::ComputationLimit(_)) => {
                found.clear();
                isolate(exact, a, b, 60, &mut exact_budget, &mut found)?;
            }
            other => other?,
        }
        out.extend(found);
    }
    // A root in the sliver beyond pi is the one near -pi, found twice.
    out.sort_by(|x, y| x[0].total_cmp(&y[0]));
    let tau = 2.0 * std::f64::consts::PI;
    if out.len() >= 2 {
        let (first, last) = (out[0], out[out.len() - 1]);
        if last[0] - tau <= first[1] + 1e-9 && first[0] + tau <= last[1] + 1e-9 {
            out.pop();
        }
    }
    Ok(out)
}

/// A cylinder and a cone: loops between consecutive roots of `D` where it is
/// positive, two rings where it is positive throughout, nothing where it is
/// negative throughout.
fn by_roots(ruled: &Surface, other: &Surface) -> Result<Vec<Component>> {
    let fast = Setup::<Fast>::of(ruled, other)?;
    let exact = Setup::<I>::of(ruled, other)?;
    let r = roots(&fast, &exact)?;
    let sign_at = |u: f64| positive(&exact.discriminant(&I::exact_f64(u))?);
    if r.is_empty() {
        return Ok(if sign_at(0.0)? {
            vec![
                Component::Ring {
                    branch: Branch::Plus,
                },
                Component::Ring {
                    branch: Branch::Minus,
                },
            ]
        } else {
            Vec::new()
        });
    }
    let tau = pi::<I>().mul(&I::exact_f64(2.0));
    let n = r.len();
    let mut out = Vec::new();
    for k in 0..n {
        let (u0, mut u1) = (r[k], r[(k + 1) % n]);
        if k + 1 == n {
            let shift = |x: f64| bounds(&I::exact_f64(x).add(&tau));
            u1 = [shift(u1[0])[0], shift(u1[1])[1]];
        }
        let mid = 0.25 * (u0[0] + u0[1] + u1[0] + u1[1]);
        if sign_at(mid)? {
            out.push(Component::Loop { u: [u0, u1] });
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------ entry

/// The outcome for the pairs parameterised here.
pub(crate) enum Found {
    Empty,
    Point(Enclosure3),
    Curve(Box<ProceduralCurve>),
}

/// `None` for pairs not parameterised here (coaxial pairs are S7a's; two
/// cones and a cone whose rational apex lies on a sphere are not yet).
pub(crate) fn intersect(a: &Surface, b: &Surface) -> Result<Option<Found>> {
    let (Some(qa), Some(qb)) = (Quad::of(a), Quad::of(b)) else {
        return Ok(None);
    };
    let curve = |ruled: &Surface, other: &Surface, components: Vec<Component>| {
        Ok(if components.is_empty() {
            Found::Empty
        } else {
            Found::Curve(Box::new(ProceduralCurve::new(
                ruled.clone(),
                other.clone(),
                components,
            )?))
        })
    };
    let off_axis = |axis: &Quad, p: &X| !perpendicular(&sub(p, &axis.o), &axis.a).iter().all(zero);
    let parallel = |p: &Quad, o: &Quad| cross(&p.a, &o.a).iter().all(zero);
    match (qa.kind, qb.kind) {
        (Kind::Cylinder, Kind::Cylinder) => {
            if parallel(&qa, &qb) {
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
            closed(ps, os, classify_cylinders(p, o), p, o).map(Some)
        }
        (Kind::Cylinder, Kind::Sphere) | (Kind::Sphere, Kind::Cylinder) => {
            let (c, s, cs, ss) = if qa.kind == Kind::Cylinder {
                (&qa, &qb, a, b)
            } else {
                (&qb, &qa, b, a)
            };
            if !off_axis(c, &s.o) {
                return Ok(None);
            }
            closed(cs, ss, classify_sphere(c, s), c, s).map(Some)
        }
        (Kind::Cone, Kind::Sphere) | (Kind::Sphere, Kind::Cone) => {
            let (k, s, ks, ss) = if qa.kind == Kind::Cone {
                (&qa, &qb, a, b)
            } else {
                (&qb, &qa, b, a)
            };
            if !off_axis(k, &s.o) {
                return Ok(None);
            }
            cone_sphere(k, s)?
                .map(|components| curve(ks, ss, components))
                .transpose()
        }
        (Kind::Cylinder, Kind::Cone) | (Kind::Cone, Kind::Cylinder) => {
            let (c, k, cs, ks) = if qa.kind == Kind::Cylinder {
                (&qa, &qb, a, b)
            } else {
                (&qb, &qa, b, a)
            };
            if parallel(c, k) && !off_axis(c, &k.o) {
                return Ok(None);
            }
            curve(cs, ks, by_roots(cs, ks)?).map(Some)
        }
        _ => Ok(None),
    }
}

/// Two cylinders or a cylinder and a sphere: their exact class and, for a
/// loop, its closed-form range.
fn closed(ruled: &Surface, other: &Surface, class: Class, p: &Quad, o: &Quad) -> Result<Found> {
    let curve = |components| -> Result<Found> {
        Ok(Found::Curve(Box::new(ProceduralCurve::new(
            ruled.clone(),
            other.clone(),
            components,
        )?)))
    };
    Ok(match class {
        Class::Empty => Found::Empty,
        // Tangent at u = 0, towards the other surface.
        Class::Point => {
            let fast = Setup::<Fast>::of(ruled, other)?;
            Found::Point(match fast.point(&Fast::exact_f64(0.0), Branch::Plus) {
                Err(Error::ComputationLimit(_)) => {
                    bounds3(&Setup::<I>::of(ruled, other)?.point(&I::exact_f64(0.0), Branch::Plus)?)
                }
                other => bounds3(&other?),
            })
        }
        Class::Rings => curve(vec![
            Component::Ring {
                branch: Branch::Plus,
            },
            Component::Ring {
                branch: Branch::Minus,
            },
        ])?,
        // The node is on the far side, u = pi: the thinner cylinder touches
        // the thicker from inside away from its axis, Viviani's sphere away
        // from its centre.
        Class::FigureEight => curve(vec![Component::FigureEight {
            node: bounds(&pi::<I>()),
        }])?,
        Class::Loop => {
            let t = half_range(p, o)?;
            curve(vec![Component::Loop {
                u: [[-t[1], -t[0]], t],
            }])?
        }
    })
}
