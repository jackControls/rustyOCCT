//! Two cones, and a cone whose rational apex lies on a sphere or a cylinder
//! (S7b.4 of `REVIEW_NOTES.md`): traced curves on the cone's rulings.
//!
//! On the first cone's rulings through its apex `V` (of two cones the first
//! by stored data), `V + v d(u)`, `d(u) = cos h a + sin h (cos u x + sin u y)`
//! with `x` the other axis's component normal to `a` (the direction to the
//! other's origin or centre when parallel), the other quadric
//! `(p - q)^T Q (p - q) - k` is `A(u) v^2 + 2 B(u) v + C`: `A = Q(d, d)`,
//! `B = Q(V - q, d)`, `C = Q(V - q, V - q) - k`. With `v = tan(t / 2)`,
//! `G(u, t) = ((A + C) + (C - A) cos t) / 2 + B sin t` is smooth on the torus
//! `(u, t)` mod `2 pi` and its zero set is the curve with its points at
//! infinity (`t = pi`): the traced graph of `torus_curves.rs` applies, and a
//! component's crossings of `t = pi` are counted (an unbounded component).
//!
//! Factors, decided exactly: the apex on the other surface (`C = 0`, a
//! rational apex on a sphere or a cylinder) leaves `A sin psi + 2 B cos psi`
//! (`psi = t / 2`, period `pi`), one point per ruling, through the apex where
//! `B` vanishes; parallel cones of equal half-angles (`A = 0` identically)
//! leave `2 B sin psi + C cos psi`, the conic in their radical plane. Two
//! cones with one rational apex meet in their common generatrices, the
//! certified roots of `A(u)`. A tangency elsewhere needs a symmetric
//! configuration (a cone's `cos` and `sin` are transcendental), such as two
//! congruent cones crossing symmetrically, which meet in two conics crossing
//! at nodes: those are not certified (`ComputationLimit`).
use super::analytic::{self, AnalyticItem, SurfaceIntersection};
use super::procedural::{
    bounds3, cross, dot, e3, eadd, ecross, edot, escale, esub, i, limit, perpendicular, q, span,
    sub, unit, zero, E, X,
};
use super::torus_curves::{
    angle, graph, roots_along, Chart, ChartKind, Jet, Node, TracedComponent, TracedCurve,
};
use crate::certified::{Fast, Interval as I, Real};
use crate::topology::Surface;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

const TAU: f64 = 2.0 * std::f64::consts::PI;

/// A cone's, cylinder's or sphere's exact data.
struct Data {
    o: X,
    a: X,
    r: R,
    angle: f64,
    kind: Kind,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Cone,
    Cylinder,
    Sphere,
}

fn data(s: &Surface) -> Option<Data> {
    let (f, r, angle, kind) = match s {
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => (frame, *radius, *half_angle, Kind::Cone),
        Surface::Cylinder { frame, radius } => (frame, *radius, 0.0, Kind::Cylinder),
        Surface::Sphere { frame, radius } => (frame, *radius, 0.0, Kind::Sphere),
        _ => return None,
    };
    Some(Data {
        o: f.origin().to_array().map(q),
        a: f.normal().to_array().map(q),
        r: q(r),
        angle,
        kind,
    })
}

/// Of two cones the first by stored normal, then origin.
fn first(a: &Surface, b: &Surface) -> bool {
    let key = |s: &Surface| {
        let Surface::Cone { frame, .. } = s else {
            unreachable!("cones only")
        };
        let [x, y, z] = frame.normal().to_array();
        let [u, v, w] = frame.origin().to_array();
        [x, y, z, u, v, w]
    };
    key(a)
        .iter()
        .zip(&key(b))
        .map(|(x, y)| x.total_cmp(y))
        .find(|o| *o != Ordering::Equal)
        .unwrap_or(Ordering::Equal)
        != Ordering::Greater
}

/// The rational apex exactly on a sphere or a cylinder.
fn apex_on(cone: &Data, other: &Data) -> bool {
    if !zero(&cone.r) || other.kind == Kind::Cone {
        return false;
    }
    let w = sub(&cone.o, &other.o);
    let ww = dot(&w, &w);
    let lhs = match other.kind {
        Kind::Sphere => ww,
        _ => {
            let h = dot(&w, &other.a);
            ww - &h * &h / dot(&other.a, &other.a)
        }
    };
    lhs == &other.r * &other.r
}

/// Parallel cones of equal half-angles.
fn twins(cone: &Data, other: &Data) -> bool {
    other.kind == Kind::Cone
        && cone.angle == other.angle
        && cross(&cone.a, &other.a).iter().all(zero)
}

/// Two cones with one rational apex.
fn common_apex(cone: &Data, other: &Data) -> bool {
    other.kind == Kind::Cone && zero(&cone.r) && zero(&other.r) && cone.o == other.o
}

/// A cone's apex in a tier.
fn apex<T: Real>(d: &Data) -> Result<E<T>> {
    if zero(&d.r) {
        return Ok(e3(&d.o));
    }
    let (c, s) = T::cos_sin(&T::exact_f64(d.angle));
    let back = i::<T>(&d.r).mul(&c).div(&s).ok_or(limit("a cone's apex"))?;
    Ok(esub(&e3(&d.o), &escale(&unit::<T>(&d.a)?, &back)))
}

/// The cone's rulings and the other quadric, in a tier.
struct Ruling<T> {
    v: E<T>,
    a: E<T>,
    x: E<T>,
    y: E<T>,
    ch: T,
    sh: T,
    /// The other quadric `(p - q)^T Q (p - q) - k`, `Q = c2 - u u^T` (`u` zero
    /// for a sphere).
    q: E<T>,
    u: E<T>,
    c2: T,
    k: T,
    kind: ChartKind,
}

impl<T: Real> Ruling<T> {
    fn of(cone: &Surface, other: &Surface, kind: ChartKind) -> Result<Self> {
        let (Some(c), Some(o)) = (data(cone), data(other)) else {
            return Err(Error::OutOfDomain("a cone's traced curve"));
        };
        let toward = {
            let w = if o.kind == Kind::Sphere {
                [q(0.0), q(0.0), q(0.0)]
            } else {
                perpendicular(&o.a, &c.a)
            };
            if w.iter().all(zero) {
                perpendicular(&sub(&o.o, &c.o), &c.a)
            } else {
                w
            }
        };
        let a = unit::<T>(&c.a)?;
        let x = unit::<T>(&toward)?;
        let y = ecross(&a, &x);
        let (ch, sh) = T::cos_sin(&T::exact_f64(c.angle));
        let zero3 = || [T::exact_f64(0.0), T::exact_f64(0.0), T::exact_f64(0.0)];
        let (qq, u, c2, k) = match o.kind {
            Kind::Sphere => (e3(&o.o), zero3(), T::exact_f64(1.0), i(&(&o.r * &o.r))),
            Kind::Cylinder => (
                e3(&o.o),
                unit::<T>(&o.a)?,
                T::exact_f64(1.0),
                i(&(&o.r * &o.r)),
            ),
            Kind::Cone => {
                let (c, _) = T::cos_sin(&T::exact_f64(o.angle));
                (
                    apex::<T>(&o)?,
                    unit::<T>(&o.a)?,
                    c.square(),
                    T::exact_f64(0.0),
                )
            }
        };
        Ok(Self {
            v: apex::<T>(&c)?,
            a,
            x,
            y,
            ch,
            sh,
            q: qq,
            u,
            c2,
            k,
            kind,
        })
    }
    fn form(&self, v: &E<T>, w: &E<T>) -> T {
        edot(v, w)
            .mul(&self.c2)
            .sub(&edot(v, &self.u).mul(&edot(w, &self.u)))
    }
    /// `d(u)` and its first two derivatives.
    fn d(&self, u: &T) -> [E<T>; 3] {
        let (c, s) = T::cos_sin(u);
        let radial = eadd(&escale(&self.x, &c), &escale(&self.y, &s));
        let turn = esub(&escale(&self.y, &c), &escale(&self.x, &s));
        [
            eadd(&escale(&self.a, &self.ch), &escale(&radial, &self.sh)),
            escale(&turn, &self.sh),
            escale(&radial, &self.sh.neg()),
        ]
    }
    /// `(A, A_u, A_uu)`, `(B, B_u, B_uu)` and `C`.
    fn abc(&self, u: &T) -> ([T; 3], [T; 3], T) {
        let [d, du, duu] = self.d(u);
        let w = esub(&self.v, &self.q);
        let two = T::exact_f64(2.0);
        (
            [
                self.form(&d, &d),
                self.form(&d, &du).mul(&two),
                self.form(&du, &du).add(&self.form(&d, &duu)).mul(&two),
            ],
            [self.form(&w, &d), self.form(&w, &du), self.form(&w, &duu)],
            self.form(&w, &w).sub(&self.k),
        )
    }
    /// The factor charts' `(alpha, beta)` and their `u`-derivatives:
    /// `alpha sin psi + beta cos psi`.
    fn factor(&self, u: &T) -> ([T; 3], [T; 3]) {
        let (a, b, c) = self.abc(u);
        let two = T::exact_f64(2.0);
        let zero = T::exact_f64(0.0);
        let b2 = [b[0].mul(&two), b[1].mul(&two), b[2].mul(&two)];
        match self.kind {
            ChartKind::Apex => (a, b2),
            _ => (b2, [c, zero.clone(), zero]),
        }
    }
}

impl<T: Real> Chart<T> for Ruling<T> {
    fn value(&self, u: &T, t: &T) -> T {
        self.jet(u, t).g
    }
    fn jet(&self, u: &T, t: &T) -> Jet<T> {
        let (ct, st) = T::cos_sin(t);
        let half = T::exact_f64(0.5);
        if self.kind == ChartKind::Rulings {
            let (a, b, c) = self.abc(u);
            // G = (A + C)/2 + (C - A)/2 cos t + B sin t.
            let one_minus = T::exact_f64(1.0).sub(&ct).mul(&half);
            let d = c.sub(&a[0]).mul(&half);
            return Jet {
                g: a[0].add(&c).mul(&half).add(&d.mul(&ct)).add(&b[0].mul(&st)),
                gp: a[1].mul(&one_minus).add(&b[1].mul(&st)),
                gt: b[0].mul(&ct).sub(&d.mul(&st)),
                gpp: a[2].mul(&one_minus).add(&b[2].mul(&st)),
                gpt: a[1].mul(&st).mul(&half).add(&b[1].mul(&ct)),
                gtt: d.mul(&ct).add(&b[0].mul(&st)).neg(),
            };
        }
        // alpha sin psi + beta cos psi, psi = t.
        let (al, be) = self.factor(u);
        let h = |x: &[T; 3], y: &[T; 3], k: usize| x[k].mul(&st).add(&y[k].mul(&ct));
        Jet {
            g: h(&al, &be, 0),
            gp: h(&al, &be, 1),
            gt: al[0].mul(&ct).sub(&be[0].mul(&st)),
            gpp: h(&al, &be, 2),
            gpt: al[1].mul(&ct).sub(&be[1].mul(&st)),
            gtt: h(&al, &be, 0).neg(),
        }
    }
    fn at(&self, u: &T, t: &T) -> Result<E<T>> {
        let (ct, st) = T::cos_sin(t);
        // v = tan(t / 2) = sin t / (1 + cos t), or tan psi on a factor chart.
        let den = if self.kind == ChartKind::Rulings {
            T::exact_f64(1.0).add(&ct)
        } else {
            ct
        };
        let v = st.div(&den).ok_or(limit("a point at infinity"))?;
        let [d, ..] = self.d(u);
        Ok(eadd(&self.v, &escale(&d, &v)))
    }
}

/// A traced curve's chart in a tier, for `TracedCurve::point_at`.
pub(super) fn chart<T: Real + 'static>(
    kind: ChartKind,
    carrier: &Surface,
    other: &Surface,
) -> Result<Box<dyn Chart<T>>> {
    Ok(Box::new(Ruling::<T>::of(carrier, other, kind)?))
}

/// The intersection of two cones not coaxial, or of a cone whose rational
/// apex lies on a sphere or a cylinder; `None` for other pairs.
pub(super) fn intersect(a: &Surface, b: &Surface) -> Result<Option<SurfaceIntersection>> {
    let (Some(da), Some(db)) = (data(a), data(b)) else {
        return Ok(None);
    };
    let (cs, os, c, o) = match (da.kind, db.kind) {
        (Kind::Cone, Kind::Cone) => {
            if first(a, b) {
                (a, b, da, db)
            } else {
                (b, a, db, da)
            }
        }
        (Kind::Cone, _) => (a, b, da, db),
        (_, Kind::Cone) => (b, a, db, da),
        _ => return Ok(None),
    };
    if o.kind == Kind::Cone {
        if common_apex(&c, &o) {
            return lines(cs, os).map(Some);
        }
    } else if !apex_on(&c, &o) {
        return Ok(None);
    }
    let kind = if apex_on(&c, &o) {
        ChartKind::Apex
    } else if twins(&c, &o) {
        ChartKind::Twins
    } else {
        ChartKind::Rulings
    };
    let fast = Ruling::<Fast>::of(cs, os, kind)?;
    let exact = Ruling::<I>::of(cs, os, kind)?;
    let period = if kind == ChartKind::Rulings {
        TAU
    } else {
        std::f64::consts::PI
    };
    let (folds, tracks, components) = graph(&fast, &exact, period, Vec::new(), &[])?;
    let mut curve = TracedCurve {
        carrier: cs.clone(),
        other: os.clone(),
        chart: kind,
        period,
        folds,
        nodes: Vec::new(),
        tracks,
        components,
    };
    count_infinity(&mut curve, &fast, &exact)?;
    if kind == ChartKind::Apex {
        attach_apex(&mut curve, &exact)?;
    }
    Ok(Some(if curve.tracks.is_empty() && curve.nodes.is_empty() {
        SurfaceIntersection::Empty
    } else {
        SurfaceIntersection::Traced(Box::new(curve))
    }))
}

/// Each component's crossings of the points at infinity: the roots of
/// `G(u, t_inf)` along `u`, each on the track through it.
fn count_infinity(curve: &mut TracedCurve, fast: &Ruling<Fast>, exact: &Ruling<I>) -> Result<()> {
    let t_inf = if curve.chart == ChartKind::Rulings {
        std::f64::consts::PI
    } else {
        0.5 * std::f64::consts::PI
    };
    // t_inf as its enclosure's middle binary64: G there is of one certain
    // sign away from the crossings, which roots_along isolates in u.
    for root in roots_along(fast, exact, t_inf)? {
        let track = curve
            .track_through(root, t_inf)
            .ok_or(limit("a crossing of infinity inside a fold's box"))?;
        for comp in curve.components.iter_mut() {
            match comp {
                TracedComponent::Smooth {
                    tracks, infinite, ..
                }
                | TracedComponent::Crossing {
                    tracks, infinite, ..
                } if tracks.contains(&track) => *infinite += 1,
                _ => {}
            }
        }
    }
    Ok(())
}

/// The apex on the other surface: crossed where `B` vanishes (two simple
/// roots: a crossing), or isolated (none); `B(u) = b0 + b1 cos u + b2 sin u`.
fn attach_apex(curve: &mut TracedCurve, exact: &Ruling<I>) -> Result<()> {
    let w = esub(&exact.v, &exact.q);
    let b0 = exact.form(&w, &escale(&exact.a, &exact.ch));
    let b1 = exact.form(&w, &escale(&exact.x, &exact.sh));
    let b2 = exact.form(&w, &escale(&exact.y, &exact.sh));
    let disc = b1.square().add(&b2.square()).sub(&b0.square());
    let crossing = match disc.sign() {
        Some(Ordering::Greater) => true,
        Some(Ordering::Less) => false,
        _ => return Err(limit("a ring tangent at the apex")),
    };
    // The first root's angle, for the node's parameters.
    let phi = if crossing {
        let r = disc.sqrt();
        let base = angle(&b2, &b1)?;
        let off = angle(&r, &b0.neg())?;
        let (lo, hi) = base.sub(&off).bounds_f64();
        [lo, hi]
    } else {
        [0.0, 0.0]
    };
    curve.nodes.push(Node {
        phi,
        t: [0.0, 0.0],
        point: bounds3(&exact.v),
        crossing,
    });
    if crossing {
        let all: Vec<usize> = (0..curve.tracks.len()).collect();
        let infinite = curve
            .components
            .iter()
            .map(|c| match c {
                TracedComponent::Smooth { infinite, .. } => *infinite,
                _ => 0,
            })
            .sum();
        curve.components = vec![TracedComponent::Crossing {
            tracks: all,
            folds: 0,
            nodes: vec![0],
            infinite,
        }];
    } else {
        curve.components.push(TracedComponent::Isolated { node: 0 });
    }
    Ok(())
}

/// Two cones with one rational apex: the lines along the rulings where `A`
/// vanishes, or the apex alone.
fn lines(cone: &Surface, other: &Surface) -> Result<SurfaceIntersection> {
    let fast = Ruling::<Fast>::of(cone, other, ChartKind::Rulings)?;
    let exact = Ruling::<I>::of(cone, other, ChartKind::Rulings)?;
    // A(u) = G(u, pi) with C = 0.
    let roots = roots_along(&fast, &exact, std::f64::consts::PI)?;
    if roots.is_empty() {
        return Ok(analytic::sorted(vec![AnalyticItem::Point(bounds3(
            &exact.v,
        ))]));
    }
    let items = roots
        .iter()
        .map(|u| {
            let [d, ..] = exact.d(&span(u[0], u[1]));
            analytic::line(&exact.v, &d)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(analytic::sorted(items))
}
