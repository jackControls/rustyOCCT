//! A torus and a cylinder, a cone or another torus off its axis (S7b.3b of
//! `REVIEW_NOTES.md`): traced curves on the torus's meridians.
//!
//! The torus (axis `a` normalised exactly, major radius `R`, minor `r`) is
//! `p(phi, t) = o + (R + r cos t) e(phi) + r sin t a`, `e(phi) = cos phi x +
//! sin phi y`, `x` the other axis's component normal to `a` (the direction to
//! the other's origin when the axes are parallel). The curve is the zero set
//! of `G(phi, t) = f(p(phi, t))` on the flat parameter torus, `f` the other
//! surface's implicit function `(p - q)^T Q (p - q) - k`, `Q = c2 - u u^T`
//! (a cylinder: `q` its origin, `c2 = 1`, `k = r_c^2`; a cone, both nappes:
//! `q` its apex, `c2 = cos^2 h`, `k = 0`; `u` the unit axis), or another
//! torus's quartic `(|w|^2 + R2^2 - r2^2)^2 - 4 R2^2 |w - (w . u) u|^2`,
//! `w = p - q` (S7b.3b.2; of two tori the carrier is the first by stored
//! data). Along a meridian circle `|w|^2` is affine in `cos t`, `sin t`, so
//! `G` is of degree two in them for every pair.
//!
//! The curve is a graph. Its vertices are the tangencies of the surfaces,
//! decided exactly from the pipes' spines and axes (`tangency.rs`); its folds, where a component turns in
//! `phi` (`G = G_t = 0`), are found by subdivision of the parameter torus
//! with mean-value exclusion and certified by the Krawczyk operator. Each fold
//! and tangency gets a box whose top and bottom edges carry no zero of `G` and
//! whose sides carry the certified simple roots of its local picture (two and
//! none for a fold, two and two for a crossing, none for an isolated point).
//! Between them the branches are graphs `t(phi)`, followed by chains of
//! certified windows (`G` changes sign across the window for every `phi` in
//! the step, `G_t` keeps a sign on it); the windows fix how branches continue
//! and give every point evaluation its bracket. Components are the connected
//! pieces; a closed smooth one keeps its winding numbers on the torus.
use super::analytic::{Enclosure, Enclosure3};
use super::procedural::{
    bounds, bounds3, e3, eadd, ecross, edot, escale, esub, i, limit, perpendicular, q, span, sub,
    unit, zero, E, X,
};
use super::tangency;
use crate::certified::{Fast, Interval as I, Real};
use crate::topology::Surface;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// The binary64 value just above pi (next_up needs Rust 1.86).
const PI_HI: f64 = 3.1415926535897936;

/// The binary64 value just above a positive `x`.
fn above(x: f64) -> f64 {
    f64::from_bits(x.to_bits() + 1)
}
const TAU: f64 = 2.0 * std::f64::consts::PI;

// ------------------------------------------------------------------ public

/// A fold: where a component turns in `phi` (`G = G_t = 0`).
#[derive(Debug, Clone, PartialEq)]
pub struct Fold {
    pub phi: Enclosure,
    pub t: Enclosure,
    pub point: Enclosure3,
}

/// A tangency of the two surfaces: a singular point of the curve.
#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    pub phi: Enclosure,
    pub t: Enclosure,
    pub point: Enclosure3,
    /// Two branches cross here; otherwise the point is isolated.
    pub crossing: bool,
}

/// A connected component of a traced curve.
#[derive(Debug, Clone, PartialEq)]
pub enum TracedComponent {
    /// A closed smooth curve: its tracks, its folds and its winding numbers
    /// `[w_phi, w_t]` on the torus (`w_phi > 0`, or `w_t >= 0` when
    /// `w_phi = 0`).
    Smooth {
        tracks: Vec<usize>,
        folds: usize,
        winding: [i64; 2],
        /// Its crossings of the points at infinity (a cone's rulings: an
        /// unbounded component crosses them, and may return); zero on a
        /// torus.
        infinite: usize,
    },
    /// Branches through crossings: their tracks, folds and nodes.
    Crossing {
        tracks: Vec<usize>,
        folds: usize,
        nodes: Vec<usize>,
        infinite: usize,
    },
    /// An isolated tangency point.
    Isolated { node: usize },
}

/// A branch `t(phi)` over an unwrapped `phi`-range, by its certified windows.
#[derive(Debug, Clone, PartialEq)]
pub struct Track {
    /// The unwrapped `phi` range, increasing.
    pub phi: [f64; 2],
    steps: Vec<Step>,
}

#[derive(Debug, Clone, PartialEq)]
struct Step {
    phi: [f64; 2],
    /// Unwrapped `t`: exactly one root of `G(phi, .)` for every `phi` in the
    /// step.
    window: [f64; 2],
}

/// Which chart a traced curve lives on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ChartKind {
    /// A torus's meridians (S7b.3b).
    Meridians,
    /// A cone's rulings with `v = tan(t / 2)` (S7b.4).
    Rulings,
    /// A cone's rulings through its apex on the other surface: the factor
    /// `A sin psi + 2 B cos psi` (S7b.4).
    Apex,
    /// Parallel cones of equal half-angles: the factor `2 B sin psi +
    /// C cos psi` (S7b.4).
    Twins,
}

/// A torus's curve with a cylinder, a cone or another torus, or a cone's
/// with another quadric, as a graph of tracks, folds and tangencies (D13's
/// procedural curve for S7b.3b and S7b.4).
#[derive(Debug, Clone, PartialEq)]
pub struct TracedCurve {
    pub(super) carrier: Surface,
    pub(super) other: Surface,
    pub(super) chart: ChartKind,
    /// The chart's period in `t`: `2 pi`, or `pi` for a factor's `psi`.
    pub(super) period: f64,
    pub(super) folds: Vec<Fold>,
    pub(super) nodes: Vec<Node>,
    pub(super) tracks: Vec<Track>,
    pub(super) components: Vec<TracedComponent>,
}

impl TracedCurve {
    /// The torus the curve is parameterised on.
    pub fn carrier(&self) -> &Surface {
        &self.carrier
    }
    pub fn other(&self) -> &Surface {
        &self.other
    }
    pub fn folds(&self) -> &[Fold] {
        &self.folds
    }
    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }
    pub fn tracks(&self) -> &[Track] {
        &self.tracks
    }
    pub fn components(&self) -> &[TracedComponent] {
        &self.components
    }
    /// The period of the chart's second parameter: `2 pi` (a torus's tube
    /// angle, or `t` with `v = tan(t / 2)` on a cone's rulings), or `pi` for
    /// the factor `psi` (`v = tan psi`) of a cone through its apex on the
    /// other surface or of parallel twin cones.
    pub fn period(&self) -> f64 {
        self.period
    }
    /// The curve's point on `track` at meridian angle `phi` (any turn of it
    /// inside the track's range), enclosed; an error outside the range.
    pub fn point_at(&self, track: usize, phi: f64) -> Result<Enclosure3> {
        Ok(self.locate(track, phi)?.1)
    }
    /// The angle `t` round the tube (unwrapped along the track) of the
    /// curve's point on `track` at `phi`, enclosed.
    pub fn t_at(&self, track: usize, phi: f64) -> Result<Enclosure> {
        Ok(self.locate(track, phi)?.0)
    }
    fn locate(&self, track: usize, phi: f64) -> Result<(Enclosure, Enclosure3)> {
        let tr = self
            .tracks
            .get(track)
            .ok_or(Error::OutOfDomain("a traced curve's track"))?;
        let k = ((tr.phi[0] - phi) / TAU).ceil();
        // The range's end up to rounding of the caller's arithmetic.
        let at = (phi + k * TAU).min(tr.phi[1]).max(tr.phi[0]);
        if (phi + k * TAU - at).abs() > 1e-9 * (1.0 + at.abs()) {
            return Err(Error::OutOfDomain("a meridian outside the track"));
        }
        let step = tr
            .steps
            .iter()
            .find(|s| s.phi[0] <= at && at <= s.phi[1])
            .ok_or(Error::OutOfDomain("a meridian outside the track"))?;
        let fast = self.chart::<Fast>()?;
        let fast = fast.as_ref();
        let narrow = |p: &Enclosure3| {
            p.iter()
                .all(|[a, b]| b - a <= 1e-12 * a.abs().max(b.abs()).max(1.0))
        };
        // Binary64 intervals resolve t only to about their rounding over
        // |G_t|; rational ones narrow further from their bracket.
        let bracket = match root_in(fast, at, step.window) {
            Ok(t) => match point(fast, at, t) {
                Ok(p) if narrow(&p) => return Ok((t, p)),
                _ => t,
            },
            Err(_) => step.window,
        };
        let exact = self.chart::<I>()?;
        let t = newton(exact.as_ref(), at, bracket)?;
        Ok((t, point(exact.as_ref(), at, t)?))
    }
    /// The chart in a certified tier.
    fn chart<T: Real + 'static>(&self) -> Result<Box<dyn Chart<T>>> {
        Ok(match self.chart {
            ChartKind::Meridians => Box::new(Field::<T>::of(&self.carrier, &self.other)?),
            kind => super::ruled_curves::chart::<T>(kind, &self.carrier, &self.other)?,
        })
    }
}

// ------------------------------------------------------------------ field

/// A surface's exact data.
struct Data {
    o: X,
    a: X,
    r: R,
    minor: R,
    angle: f64,
    cone: bool,
    torus: bool,
}

fn data(s: &Surface) -> Option<Data> {
    let (f, r, minor, angle, cone, torus) = match s {
        Surface::Torus {
            frame,
            major,
            minor,
        } => (frame, *major, *minor, 0.0, false, true),
        Surface::Cylinder { frame, radius } => (frame, *radius, 0.0, 0.0, false, false),
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => (frame, *radius, 0.0, *half_angle, true, false),
        _ => return None,
    };
    Some(Data {
        o: f.origin().to_array().map(q),
        a: f.normal().to_array().map(q),
        r: q(r),
        minor: q(minor),
        angle,
        cone,
        torus,
    })
}

/// `G` on the parameter torus in a certified tier.
struct Field<T> {
    o: E<T>,
    a: E<T>,
    x: E<T>,
    y: E<T>,
    big: T,
    small: T,
    /// `f(p) = (p - q)^T Q (p - q) - k`, `Q = c2 - u u^T`; for another
    /// torus (`torus_pair`) `c2 = 1` and `f = (|w|^2 + K)^2 - 4 R2^2
    /// (p - q)^T Q (p - q)`, `w = p - q`, with `(R2^2, K = R2^2 - r2^2)`.
    q: E<T>,
    u: E<T>,
    c2: T,
    k: T,
    quartic: Option<(T, T)>,
}

/// `G` and its derivatives up to the second order.
pub(super) struct Jet<T> {
    pub(super) g: T,
    pub(super) gp: T,
    pub(super) gt: T,
    pub(super) gpp: T,
    pub(super) gpt: T,
    pub(super) gtt: T,
}

/// A function `G(u, t)` on a parameter torus (`u` mod `2 pi`, `t` mod the
/// chart's period) whose zero set is the curve, with the point of the curve
/// at a parameter: the torus's meridians here, a cone's rulings in
/// `ruled_curves.rs`.
pub(super) trait Chart<T: Real> {
    fn value(&self, u: &T, t: &T) -> T;
    fn jet(&self, u: &T, t: &T) -> Jet<T>;
    /// The curve's point at `(u, t)`; an error at infinity.
    fn at(&self, u: &T, t: &T) -> Result<E<T>>;
}

impl<T: Real> Chart<T> for Field<T> {
    fn value(&self, u: &T, t: &T) -> T {
        Field::value(self, u, t)
    }
    fn jet(&self, u: &T, t: &T) -> Jet<T> {
        Field::jet(self, u, t)
    }
    fn at(&self, u: &T, t: &T) -> Result<E<T>> {
        let [p, ..] = self.frame(u, t);
        Ok(p)
    }
}

/// The rational `x` direction of the torus's frame.
fn toward(t: &Data, o: &Data) -> X {
    let w = perpendicular(&o.a, &t.a);
    if w.iter().all(zero) {
        perpendicular(&sub(&o.o, &t.o), &t.a)
    } else {
        w
    }
}

impl<T: Real> Field<T> {
    fn of(torus: &Surface, other: &Surface) -> Result<Self> {
        let (Some(t), Some(o)) = (data(torus), data(other)) else {
            return Err(Error::OutOfDomain("a traced torus curve"));
        };
        let a = unit::<T>(&t.a)?;
        let x = unit::<T>(&toward(&t, &o))?;
        let y = ecross(&a, &x);
        let u = unit::<T>(&o.a)?;
        let quartic = o.torus.then(|| {
            (
                i::<T>(&(&o.r * &o.r)),
                i::<T>(&(&o.r * &o.r - &o.minor * &o.minor)),
            )
        });
        let (q, c2, k) = if o.cone {
            let (c, s) = T::cos_sin(&T::exact_f64(o.angle));
            let apex = if zero(&o.r) {
                e3(&o.o)
            } else {
                let back = i::<T>(&o.r).mul(&c).div(&s).ok_or(limit("a cone's apex"))?;
                esub(&e3(&o.o), &escale(&u, &back))
            };
            (apex, c.square(), T::exact_f64(0.0))
        } else {
            (e3(&o.o), T::exact_f64(1.0), i(&(&o.r * &o.r)))
        };
        Ok(Self {
            o: e3(&t.o),
            a,
            x,
            y,
            big: i(&t.r),
            small: i(&t.minor),
            q,
            u,
            c2,
            k,
            quartic,
        })
    }
    fn qv(&self, w: &E<T>) -> E<T> {
        esub(&escale(w, &self.c2), &escale(&self.u, &edot(w, &self.u)))
    }
    fn form(&self, v: &E<T>, w: &E<T>) -> T {
        edot(v, w)
            .mul(&self.c2)
            .sub(&edot(v, &self.u).mul(&edot(w, &self.u)))
    }
    /// The point and its derivatives: p, p_phi, p_t, p_pp, p_pt, p_tt.
    fn frame(&self, phi: &T, t: &T) -> [E<T>; 6] {
        let (cp, sp) = T::cos_sin(phi);
        let (ct, st) = T::cos_sin(t);
        let e = eadd(&escale(&self.x, &cp), &escale(&self.y, &sp));
        let de = esub(&escale(&self.y, &cp), &escale(&self.x, &sp));
        let rho = self.big.add(&self.small.mul(&ct));
        let rs = self.small.mul(&st);
        let rc = self.small.mul(&ct);
        let p = eadd(&eadd(&self.o, &escale(&e, &rho)), &escale(&self.a, &rs));
        let pp = escale(&de, &rho);
        let pt = esub(&escale(&self.a, &rc), &escale(&e, &rs));
        let ppp = escale(&e, &rho.neg());
        let ppt = escale(&de, &rs.neg());
        let ptt = escale(
            &eadd(&escale(&e, &rc), &escale(&self.a, &rs)),
            &T::exact_f64(-1.0),
        );
        [p, pp, pt, ppp, ppt, ptt]
    }
    fn value(&self, phi: &T, t: &T) -> T {
        let [p, ..] = self.frame(phi, t);
        let w = esub(&p, &self.q);
        match &self.quartic {
            None => self.form(&w, &w).sub(&self.k),
            Some((r22, kk)) => edot(&w, &w)
                .add(kk)
                .square()
                .sub(&self.form(&w, &w).mul(r22).mul(&T::exact_f64(4.0))),
        }
    }
    fn jet(&self, phi: &T, t: &T) -> Jet<T> {
        let [p, pp, pt, ppp, ppt, ptt] = self.frame(phi, t);
        let w = esub(&p, &self.q);
        let qw = self.qv(&w);
        let two = T::exact_f64(2.0);
        // The gradient of f and its Hessian as a form.
        let (g, grad, s, r22) = match &self.quartic {
            None => (
                edot(&w, &qw).sub(&self.k),
                escale(&qw, &two),
                None,
                T::exact_f64(0.0),
            ),
            Some((r22, kk)) => {
                let s = edot(&w, &w).add(kk);
                let four = T::exact_f64(4.0);
                let g = s.square().sub(&edot(&w, &qw).mul(r22).mul(&four));
                let grad = esub(
                    &escale(&w, &s.mul(&four)),
                    &escale(&qw, &r22.mul(&T::exact_f64(8.0))),
                );
                (g, grad, Some(s), r22.clone())
            }
        };
        let d = |v: &E<T>| edot(&grad, v);
        let h = |v: &E<T>, z: &E<T>| match &s {
            None => self.form(v, z).mul(&two),
            Some(s) => edot(v, z)
                .mul(s)
                .mul(&T::exact_f64(4.0))
                .add(&edot(&w, v).mul(&edot(&w, z)).mul(&T::exact_f64(8.0)))
                .sub(&self.form(v, z).mul(&r22).mul(&T::exact_f64(8.0))),
        };
        Jet {
            g,
            gp: d(&pp),
            gt: d(&pt),
            gpp: h(&pp, &pp).add(&d(&ppp)),
            gpt: h(&pp, &pt).add(&d(&ppt)),
            gtt: h(&pt, &pt).add(&d(&ptt)),
        }
    }
}

fn point<T: Real>(f: &dyn Chart<T>, phi: f64, t: [f64; 2]) -> Result<Enclosure3> {
    Ok(bounds3(&f.at(&T::exact_f64(phi), &span(t[0], t[1]))?))
}

fn certain<T: Real>(x: &T) -> bool {
    matches!(x.sign(), Some(Ordering::Less | Ordering::Greater))
}
fn sign<T: Real>(x: &T) -> Option<Ordering> {
    x.sign().filter(|s| *s != Ordering::Equal)
}
fn mid(e: [f64; 2]) -> f64 {
    0.5 * e[0] + 0.5 * e[1]
}
/// `[lo, hi] - m` in a tier, exactly.
fn offset<T: Real>(e: [f64; 2], m: f64) -> T {
    span::<T>(e[0], e[1]).sub(&T::exact_f64(m))
}

// ------------------------------------------------------------------ 1D

/// The certified simple roots of `G(phi, .)` over `[lo, hi]` (unwrapped `t`,
/// `phi` exact), each narrowed by bisection: subdivision with the mean-value
/// enclosure `G(m) + G_t(piece) (piece - m)`, a piece kept as a root when `G`
/// changes sign between its ends and `G_t` keeps a sign on it.
fn roots_on<T: Real>(
    f: &dyn Chart<T>,
    phi: f64,
    lo: f64,
    hi: f64,
    budget: &mut usize,
) -> Result<Vec<Enclosure>> {
    let p = T::exact_f64(phi);
    let at = |t: f64| f.value(&p, &T::exact_f64(t));
    let mut out = Vec::new();
    let mut pending = vec![(lo, hi, 60usize)];
    while let Some((a, b, depth)) = pending.pop() {
        *budget = budget.checked_sub(1).ok_or(limit("a meridian's roots"))?;
        let m = 0.5 * a + 0.5 * b;
        let j = f.jet(&p, &span(a, b));
        if certain(&j.g) {
            continue;
        }
        let mv = at(m).add(&j.gt.mul(&offset::<T>([a, b], m)));
        if certain(&mv) {
            continue;
        }
        let (sa, sb) = (sign(&at(a)), sign(&at(b)));
        if certain(&j.gt) {
            // Monotone on the piece: a root where the ends' signs differ,
            // or exactly at an end (a subdivision point can hit a symmetric
            // root), else none.
            let exact = |x: f64| at(x).sign() == Some(Ordering::Equal);
            if sa.is_some() && sb.is_some() {
                if sa != sb {
                    out.push(newton(f, phi, [a, b])?);
                }
                continue;
            }
            if sa.is_none() && exact(a) || sb.is_none() && exact(b) {
                let x = if sa.is_none() && exact(a) { a } else { b };
                if !out.iter().any(|r: &Enclosure| r[0] == x && r[1] == x) {
                    out.push([x, x]);
                }
                continue;
            }
        }
        if depth == 0 || !(a < m && m < b) {
            return Err(limit("a meridian's roots"));
        }
        pending.push((m, b, depth - 1));
        pending.push((a, m, depth - 1));
    }
    out.sort_by(|x, y| x[0].total_cmp(&y[0]));
    Ok(out)
}

/// A simple root's bracket `[a, b]` narrowed by bisection on certain signs.
fn narrow(
    s: impl Fn(f64) -> Option<Ordering>,
    mut a: f64,
    mut b: f64,
    at_a: Option<Ordering>,
) -> Enclosure {
    for _ in 0..80 {
        let m = 0.5 * a + 0.5 * b;
        if !(a < m && m < b) {
            break;
        }
        match s(m) {
            Some(x) if Some(x) == at_a => a = m,
            Some(_) => b = m,
            None => break,
        }
    }
    [a, b]
}

/// The root of `G(phi, .)` in a certified window.
fn root_in<T: Real>(f: &dyn Chart<T>, phi: f64, w: [f64; 2]) -> Result<Enclosure> {
    let p = T::exact_f64(phi);
    let s = |t: f64| sign(&f.value(&p, &T::exact_f64(t)));
    let (a, b) = (s(w[0]), s(w[1]));
    if a.is_none() || b.is_none() || a == b {
        return Err(limit("a track's root"));
    }
    Ok(narrow(s, w[0], w[1], a))
}

/// The simple root of `G(phi, .)` in a certified bracket (one root, `G_t` of
/// one sign) narrowed by the interval Newton operator `m - G(m) / G_t(T)`,
/// which keeps every root of the bracket: quadratically, a few evaluations
/// where bisection would need dozens.
fn newton<T: Real>(f: &dyn Chart<T>, phi: f64, w: [f64; 2]) -> Result<Enclosure> {
    let p = T::exact_f64(phi);
    let mut b = w;
    for _ in 0..12 {
        let m = mid(b);
        if !(b[0] < m && m < b[1]) {
            break;
        }
        let gm = f.value(&p, &T::exact_f64(m));
        let gt = f.jet(&p, &span(b[0], b[1])).gt;
        let Some(q) = gm.div(&gt) else {
            break;
        };
        let (lo, hi) = T::exact_f64(m).sub(&q).bounds_f64();
        let next = [lo.max(b[0]), hi.min(b[1])];
        if next[0] > next[1] {
            return Err(limit("a track's root"));
        }
        let stalled = next[1] - next[0] >= 0.5 * (b[1] - b[0]);
        b = next;
        if stalled {
            break;
        }
    }
    Ok(b)
}

/// The simple root of `G(., t)` in a certified bracket (`G_u` of one sign)
/// narrowed by the interval Newton operator `m - G(m) / G_u(bracket)`,
/// which keeps the root.
fn newton_along<T: Real>(f: &dyn Chart<T>, t: &T, w: [f64; 2]) -> [f64; 2] {
    let mut b = w;
    for _ in 0..4 {
        let m = mid(b);
        if !(b[0] < m && m < b[1]) {
            break;
        }
        let gm = f.value(&T::exact_f64(m), t);
        let Some(q) = gm.div(&f.jet(&span(b[0], b[1]), t).gp) else {
            break;
        };
        let (lo, hi) = T::exact_f64(m).sub(&q).bounds_f64();
        let next = [lo.max(b[0]), hi.min(b[1])];
        if next[0] > next[1] {
            break;
        }
        let stalled = next[1] - next[0] >= 0.5 * (b[1] - b[0]);
        b = next;
        if stalled {
            break;
        }
    }
    b
}

/// Whether a piece's jet shows `G` free of zeros on it: `G`, or its
/// mean-value form `G(m) + G_u(piece) (piece - m)`, of one certain sign.
fn excluded<T: Real>(j: &Jet<T>, gm: &T, piece: [f64; 2], m: f64) -> bool {
    certain(&j.g) || certain(&gm.add(&j.gp.mul(&offset::<T>(piece, m))))
}

/// The certified simple roots of `G(., t)` over `u` in `[-pi, pi]` (`t`
/// exact): where the curve crosses the line `t` (a cone's points at
/// infinity). A piece on which `G_u` keeps a sign holds one root if `G`'s
/// certain signs at its ends differ and none if they agree. Binary64
/// intervals fail where a subdivision point lies within their rounding of a
/// root (a symmetric pair's crossing an ulp or two off `+-pi/2`); the
/// rational tier then takes what binary64 settles (a piece free of zeros,
/// a point's sign, a root's bracket down to binary64's rounding) and
/// decides the pieces holding roots itself, as before, so it evaluates
/// rational cosines only on those and within binary64's rounding of a
/// root. Its roots are bisected to adjacent binary64 values, as before.
pub(super) fn roots_along(
    fast: &dyn Chart<Fast>,
    exact: &dyn Chart<I>,
    t: f64,
) -> Result<Vec<Enclosure>> {
    fn run<T: Real>(
        f: &dyn Chart<T>,
        quick: Option<&dyn Chart<Fast>>,
        t: f64,
        budget: &mut usize,
    ) -> Result<Vec<Enclosure>> {
        let tt = T::exact_f64(t);
        let ft = Fast::exact_f64(t);
        let at = |u: f64| f.value(&T::exact_f64(u), &tt);
        let s = |u: f64| {
            quick
                .and_then(|q| sign(&q.value(&Fast::exact_f64(u), &ft)))
                .or_else(|| sign(&at(u)))
        };
        let mut out = Vec::new();
        let mut pending = vec![(-PI_HI, PI_HI, 60usize)];
        while let Some((a, b, depth)) = pending.pop() {
            *budget = budget
                .checked_sub(1)
                .ok_or(limit("a curve's crossings of infinity"))?;
            let m = 0.5 * a + 0.5 * b;
            let mut ends = None;
            if let Some(q) = quick {
                let j = q.jet(&span(a, b), &ft);
                if excluded(&j, &q.value(&Fast::exact_f64(m), &ft), [a, b], m) {
                    continue;
                }
                if certain(&j.gp) {
                    let (sa, sb) = (s(a), s(b));
                    if sa.is_some() && sa == sb {
                        continue;
                    }
                    ends = Some((sa, sb));
                }
            }
            let j = f.jet(&span(a, b), &tt);
            if excluded(&j, &at(m), [a, b], m) {
                continue;
            }
            let (sa, sb) = ends.unwrap_or_else(|| (s(a), s(b)));
            if sa.is_some() && sb.is_some() && certain(&j.gp) {
                // Monotone: a root between the ends if their signs differ,
                // else none (halving such a piece beside a root at its end
                // would descend to the floor).
                if sa != sb {
                    out.push(match quick {
                        // Binary64 bisection to within its rounding of the
                        // root, then interval Newton steps: few rational
                        // cosines before the last bisections, whose bracket
                        // is the one bisection alone would reach.
                        Some(q) => {
                            let w = narrow(|u| sign(&q.value(&Fast::exact_f64(u), &ft)), a, b, sa);
                            let w = newton_along(f, &tt, w);
                            narrow(s, w[0], w[1], sa)
                        }
                        None => narrow(s, a, b, sa),
                    });
                }
                continue;
            }
            if depth == 0 || !(a < m && m < b) {
                return Err(limit("a curve's crossings of infinity"));
            }
            pending.push((m, b, depth - 1));
            pending.push((a, m, depth - 1));
        }
        out.sort_by(|x, y| x[0].total_cmp(&y[0]));
        Ok(dedupe_seam(out, TAU))
    }
    let mut budget = 20_000;
    match run(fast, None, t, &mut budget) {
        Err(Error::ComputationLimit(_)) => run(exact, Some(fast), t, &mut budget),
        other => other,
    }
}

impl TracedCurve {
    /// The track through `(u, t)` (some turn of each), if any.
    pub(super) fn track_through(&self, u: Enclosure, t: f64) -> Option<usize> {
        let m = mid(u);
        self.tracks.iter().position(|tr| {
            let k = ((tr.phi[0] - m) / TAU).ceil();
            let at = m + k * TAU;
            tr.steps.iter().any(|s| {
                s.phi[0] <= at
                    && at <= s.phi[1]
                    && [-2.0, -1.0, 0.0, 1.0, 2.0].iter().any(|j| {
                        s.window[0] <= t + j * self.period && t + j * self.period <= s.window[1]
                    })
            })
        })
    }
}

/// Whether `G` vanishes on a box's edge `[p0, p1] x {t}`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Edge {
    /// No zero: every piece excluded.
    Free,
    /// A zero, certainly: `G` has opposite certain signs at two exact
    /// points of the edge, so no tier and no budget can free it.
    Crossed,
    /// A piece the subdivision could not settle.
    Unsettled,
}

/// Whether `G` vanishes on `[p0, p1] x {t}`: subdivision with the mean-value
/// enclosure in `phi`. A certain sign change between exact points (the ends,
/// or a piece's end and its midpoint) is a zero, reported at once: the
/// subdivision would otherwise descend to its floor there, which in rational
/// intervals costs seconds of transcendental evaluations beside a nearly
/// singular point.
fn edge_free<T: Real>(
    f: &dyn Chart<T>,
    p0: f64,
    p1: f64,
    t: f64,
    budget: &mut usize,
) -> Result<Edge> {
    let tt = T::exact_f64(t);
    let at = |x: f64| f.value(&T::exact_f64(x), &tt);
    let changes = |x: Option<Ordering>, y: Option<Ordering>| x.is_some() && y.is_some() && x != y;
    let (s0, s1) = (sign(&at(p0)), sign(&at(p1)));
    if changes(s0, s1) {
        return Ok(Edge::Crossed);
    }
    let mut pending = vec![(p0, p1, s0, s1, 40usize)];
    while let Some((a, b, sa, sb, depth)) = pending.pop() {
        *budget = budget.checked_sub(1).ok_or(limit("a box's edge"))?;
        let j = f.jet(&span(a, b), &tt);
        if certain(&j.g) {
            continue;
        }
        let m = 0.5 * a + 0.5 * b;
        let gm = at(m);
        let mv = gm.add(&j.gp.mul(&offset::<T>([a, b], m)));
        if certain(&mv) {
            continue;
        }
        if depth == 0 || !(a < m && m < b) {
            return Ok(Edge::Unsettled);
        }
        let sm = sign(&gm);
        if changes(sa, sm) || changes(sm, sb) {
            return Ok(Edge::Crossed);
        }
        pending.push((m, b, sm, sb, depth - 1));
        pending.push((a, m, sa, sm, depth - 1));
    }
    Ok(Edge::Free)
}

// ------------------------------------------------------------------ 2D

type Bx = [[f64; 2]; 2];

fn inflate(b: &Bx, k: f64) -> Bx {
    b.map(|[lo, hi]| {
        let w = (hi - lo) * k;
        [lo - w, hi + w]
    })
}

enum Kraw {
    Unique(Bx),
    Empty,
    Unknown,
}

/// The Krawczyk operator of `(G, G_t)` (a fold) or of `(G_phi, G_t)` (a
/// critical point of `G`) on a box: `Unique` when it maps the box into its
/// interior (one regular zero there, inside the returned box).
fn krawczyk<T: Real>(f: &dyn Chart<T>, b: &Bx, critical: bool) -> Kraw {
    let m = [mid(b[0]), mid(b[1])];
    let (mp, mt) = (T::exact_f64(m[0]), T::exact_f64(m[1]));
    let jm = f.jet(&mp, &mt);
    let jb = f.jet(&span(b[0][0], b[0][1]), &span(b[1][0], b[1][1]));
    let (fm, jj) = if critical {
        ([jm.gp, jm.gt], [[jb.gpp, jb.gpt.clone()], [jb.gpt, jb.gtt]])
    } else {
        ([jm.g, jm.gt], [[jb.gp, jb.gt], [jb.gpt, jb.gtt]])
    };
    // Y, the inverse of the Jacobian's midpoint.
    let c = |x: &T| {
        let (lo, hi) = x.bounds_f64();
        0.5 * lo + 0.5 * hi
    };
    let jmid = [[c(&jj[0][0]), c(&jj[0][1])], [c(&jj[1][0]), c(&jj[1][1])]];
    let det = jmid[0][0] * jmid[1][1] - jmid[0][1] * jmid[1][0];
    if !(det.is_finite() && det != 0.0) {
        return Kraw::Unknown;
    }
    let y = [
        [jmid[1][1] / det, -jmid[0][1] / det],
        [-jmid[1][0] / det, jmid[0][0] / det],
    ];
    if !y.iter().flatten().all(|v| v.is_finite()) {
        return Kraw::Unknown;
    }
    let d = [offset::<T>(b[0], m[0]), offset::<T>(b[1], m[1])];
    let mut k = [[0.0; 2]; 2];
    let mut empty = false;
    for r in 0..2 {
        let yr = [T::exact_f64(y[r][0]), T::exact_f64(y[r][1])];
        let mut acc = T::exact_f64(m[r]).sub(&yr[0].mul(&fm[0]).add(&yr[1].mul(&fm[1])));
        for (col, dc) in d.iter().enumerate() {
            let yj = yr[0].mul(&jj[0][col]).add(&yr[1].mul(&jj[1][col]));
            let delta = T::exact_f64(if r == col { 1.0 } else { 0.0 });
            acc = acc.add(&delta.sub(&yj).mul(dc));
        }
        let (lo, hi) = acc.bounds_f64();
        if !(lo.is_finite() && hi.is_finite()) {
            return Kraw::Unknown;
        }
        if hi < b[r][0] || lo > b[r][1] {
            empty = true;
        }
        k[r] = [lo, hi];
    }
    if empty {
        return Kraw::Empty;
    }
    if (0..2).all(|r| b[r][0] < k[r][0] && k[r][1] < b[r][1]) {
        Kraw::Unique(k)
    } else {
        Kraw::Unknown
    }
}

/// A fold's box refined by repeated Krawczyk steps.
fn refine<T: Real>(f: &dyn Chart<T>, mut b: Bx, critical: bool) -> Bx {
    for _ in 0..40 {
        match krawczyk(f, &inflate(&b, 1e-3), critical) {
            Kraw::Unique(k)
                if k[0][1] - k[0][0] < b[0][1] - b[0][0]
                    || k[1][1] - k[1][0] < b[1][1] - b[1][0] =>
            {
                b = [
                    [k[0][0].max(b[0][0]), k[0][1].min(b[0][1])],
                    [k[1][0].max(b[1][0]), k[1][1].min(b[1][1])],
                ];
                if b[0][0] > b[0][1] || b[1][0] > b[1][1] {
                    return k;
                }
            }
            _ => break,
        }
    }
    b
}

enum Verdict {
    Excluded,
    Fold(Bx),
    Split,
}

fn verdict<T: Real>(f: &dyn Chart<T>, b: &Bx) -> Verdict {
    let (p, t) = (span::<T>(b[0][0], b[0][1]), span::<T>(b[1][0], b[1][1]));
    let jb = f.jet(&p, &t);
    if certain(&jb.g) || certain(&jb.gt) {
        return Verdict::Excluded;
    }
    let m = [mid(b[0]), mid(b[1])];
    let jm = f.jet(&T::exact_f64(m[0]), &T::exact_f64(m[1]));
    let (dp, dt) = (offset::<T>(b[0], m[0]), offset::<T>(b[1], m[1]));
    let g = jm.g.add(&jb.gp.mul(&dp)).add(&jb.gt.mul(&dt));
    let gt = jm.gt.add(&jb.gpt.mul(&dp)).add(&jb.gtt.mul(&dt));
    if certain(&g) || certain(&gt) {
        return Verdict::Excluded;
    }
    match krawczyk(f, &inflate(b, 0.25), false) {
        Kraw::Unique(k) => Verdict::Fold(k),
        Kraw::Empty => Verdict::Excluded,
        Kraw::Unknown => Verdict::Split,
    }
}

fn inside(b: &Bx, z: &Bx) -> bool {
    (0..2).all(|r| z[r][0] <= b[r][0] && b[r][1] <= z[r][1])
}

/// Every fold on the parameter torus, outside the tangencies' boxes: a
/// subdivision of `[-pi, pi]^2` (binary64 intervals, a box they cannot settle
/// at a small size again in rational intervals) within a work budget.
fn folds(
    fast: &dyn Chart<Fast>,
    exact: &dyn Chart<I>,
    period: f64,
    zones: &[Bx],
) -> Result<Vec<Bx>> {
    let mut out: Vec<Bx> = Vec::new();
    let half = above(0.5 * period);
    let mut pending = vec![([[-PI_HI, PI_HI], [-half, half]], 0usize)];
    // The fixtures need at most about two thousand boxes and no rational
    // one (each costs tens of milliseconds): binary64 intervals fail only
    // beside a nearly singular point, where a few rational boxes rarely
    // help.
    let mut budget = 50_000usize;
    let mut exact_budget = 4usize;
    while let Some((b, depth)) = pending.pop() {
        budget = budget
            .checked_sub(1)
            .ok_or(limit("a torus curve's folds"))?;
        let shifted =
            |z: &Bx, k: f64, l: f64| [[z[0][0] + k, z[0][1] + k], [z[1][0] + l, z[1][1] + l]];
        if zones.iter().any(|z| {
            [-TAU, 0.0, TAU].iter().any(|k| {
                [-period, 0.0, period]
                    .iter()
                    .any(|l| inside(&b, &shifted(z, *k, *l)))
            })
        }) {
            continue;
        }
        let small = b.iter().all(|[lo, hi]| hi - lo < 1e-7);
        let mut v = verdict(fast, &b);
        if matches!(v, Verdict::Split) && small {
            exact_budget = exact_budget
                .checked_sub(1)
                .ok_or(limit("a torus curve's folds"))?;
            v = verdict(exact, &b);
        }
        match v {
            Verdict::Excluded => {}
            Verdict::Fold(k) => {
                let k = refine(fast, k, false);
                let same = |o: &Bx| {
                    (0..2).all(|r| {
                        let p = if r == 0 { TAU } else { period };
                        [-p, 0.0, p]
                            .iter()
                            .any(|s| o[r][0] <= k[r][1] + s && k[r][0] + s <= o[r][1])
                    })
                };
                if !out.iter().any(same) {
                    out.push(k);
                }
            }
            Verdict::Split => {
                if depth > 64 {
                    return Err(limit("a torus curve's folds"));
                }
                let r = if b[0][1] - b[0][0] >= b[1][1] - b[1][0] {
                    0
                } else {
                    1
                };
                let m = mid(b[r]);
                if !(b[r][0] < m && m < b[r][1]) {
                    return Err(limit("a torus curve's folds"));
                }
                let (mut lo, mut hi) = (b, b);
                lo[r][1] = m;
                hi[r][0] = m;
                pending.push((hi, depth + 1));
                pending.push((lo, depth + 1));
            }
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------ boxes

/// A fold's or a tangency's box: no zero of `G` on its top and bottom edges,
/// the certified roots on its sides.
#[derive(Clone)]
pub(super) struct Local {
    b: Bx,
    left: Vec<Enclosure>,
    right: Vec<Enclosure>,
    kind: LocalKind,
}

#[derive(Clone, Copy, PartialEq)]
enum LocalKind {
    Fold,
    Node(usize),
}

type Sides = (Vec<Enclosure>, Vec<Enclosure>);

/// A box's certified side roots, if its top and bottom are free of the
/// curve: `Err(true)` when an edge certainly meets it (another tier cannot
/// help), `Err(false)` when this tier could not settle an edge or a side.
fn local_edges<T: Real>(f: &dyn Chart<T>, b: &Bx) -> Result<std::result::Result<Sides, bool>> {
    let mut budget = 20_000;
    for t in [b[1][0], b[1][1]] {
        match edge_free(f, b[0][0], b[0][1], t, &mut budget)? {
            Edge::Free => {}
            Edge::Crossed => return Ok(Err(true)),
            Edge::Unsettled => return Ok(Err(false)),
        }
    }
    let left = roots_on(f, b[0][0], b[1][0], b[1][1], &mut budget);
    let right = roots_on(f, b[0][1], b[1][0], b[1][1], &mut budget);
    match (left, right) {
        (Ok(l), Ok(r)) => Ok(Ok((l, r))),
        (Err(Error::ComputationLimit(_)), _) | (_, Err(Error::ComputationLimit(_))) => {
            Ok(Err(false))
        }
        (Err(e), _) | (_, Err(e)) => Err(e),
    }
}

/// A box's side roots in binary64 intervals, else in rational ones unless
/// the binary64 tier proved the curve crosses its top or bottom.
fn box_edges(fast: &dyn Chart<Fast>, exact: &dyn Chart<I>, b: &Bx) -> Result<Option<Sides>> {
    Ok(match local_edges(fast, b)? {
        Ok(e) => Some(e),
        Err(true) => None,
        Err(false) => local_edges(exact, b)?.ok(),
    })
}

/// The box of a fold at `k`: `G_phi` of one sign on it, no zero on its top
/// and bottom, two roots on one side and none on the other.
fn fold_box(fast: &dyn Chart<Fast>, exact: &dyn Chart<I>, k: &Bx) -> Result<Local> {
    let c = [mid(k[0]), mid(k[1])];
    let j = fast.jet(&Fast::exact_f64(c[0]), &Fast::exact_f64(c[1]));
    let ratio = {
        let (gp, gtt) = (j.gp.bounds_f64(), j.gtt.bounds_f64());
        (mid([gp.0, gp.1]) / mid([gtt.0, gtt.1])).abs()
    };
    let mut hp = 1e-3f64;
    for _ in 0..12 {
        // The two branches are about sqrt(2 hp |G_phi / G_tt|) from the fold
        // at hp from it.
        let spread = (2.0 * hp * ratio).sqrt();
        let ht = 3.0 * spread + 1e-9;
        let b = [[c[0] - hp, c[0] + hp], [c[1] - ht, c[1] + ht]];
        let over = |f: &dyn Fn(&I, &I) -> bool| f(&span(b[0][0], b[0][1]), &span(b[1][0], b[1][1]));
        let gp_fast = certain(
            &fast
                .jet(&span(b[0][0], b[0][1]), &span(b[1][0], b[1][1]))
                .gp,
        );
        let gp = gp_fast || over(&|p, t| certain(&exact.jet(p, t).gp));
        if (0..2).all(|r| b[r][0] < k[r][0] && k[r][1] < b[r][1]) && gp {
            let edges = box_edges(fast, exact, &b)?;
            if let Some((left, right)) = edges {
                let counts = (left.len(), right.len());
                if counts == (2, 0) || counts == (0, 2) {
                    return Ok(Local {
                        b,
                        left,
                        right,
                        kind: LocalKind::Fold,
                    });
                }
            }
        }
        hp /= 4.0;
    }
    Err(limit("a fold's box"))
}

/// The box of a tangency: a unique critical point of `G` in it, no zero on
/// its top and bottom, two roots on each side for a crossing, none for an
/// isolated point.
fn node_box(
    fast: &dyn Chart<Fast>,
    exact: &dyn Chart<I>,
    at: &Bx,
    crossing: bool,
    index: usize,
) -> Result<Local> {
    let c = [mid(at[0]), mid(at[1])];
    let j = exact.jet(&span(at[0][0], at[0][1]), &span(at[1][0], at[1][1]));
    let m = |x: &I| {
        let (lo, hi) = x.bounds_f64();
        mid([lo, hi])
    };
    // The branches' slopes dt/dphi solve gtt s^2 + 2 gpt s + gpp = 0.
    let (a, bq, cq) = (m(&j.gtt), m(&j.gpt), m(&j.gpp));
    let slope = if crossing {
        let disc = (bq * bq - a * cq).max(0.0).sqrt();
        if a == 0.0 {
            return Err(limit("a crossing along a meridian"));
        }
        ((-bq + disc) / a).abs().max(((-bq - disc) / a).abs())
    } else {
        1.0
    };
    let mut h = 1e-2f64;
    for _ in 0..12 {
        let ht = h * (2.0 * slope + 0.5);
        let b = [[c[0] - h, c[0] + h], [c[1] - ht, c[1] + ht]];
        let unique = matches!(krawczyk(fast, &b, true), Kraw::Unique(_))
            || matches!(krawczyk(exact, &b, true), Kraw::Unique(_));
        if unique && (0..2).all(|r| b[r][0] < at[r][0] && at[r][1] < b[r][1]) {
            let edges = box_edges(fast, exact, &b)?;
            if let Some((left, right)) = edges {
                let want = if crossing { (2, 2) } else { (0, 0) };
                if (left.len(), right.len()) == want {
                    return Ok(Local {
                        b,
                        left,
                        right,
                        kind: LocalKind::Node(index),
                    });
                }
            }
        }
        h /= 4.0;
    }
    Err(limit("a tangency's box"))
}

// ------------------------------------------------------------------ tracks

/// Where a track begins or ends.
#[derive(Clone, Copy, Debug, PartialEq)]
enum End {
    /// The sweep's start or end meridian, root `j`.
    Seam(usize),
    /// A box's side: box, side (0 left, 1 right), root.
    Side(usize, usize, usize),
}

struct Traced {
    track: Track,
    ends: [End; 2],
    /// Unwrapped `t` at the ends.
    t: [f64; 2],
}

/// One certified step of a branch from `phi0` (its root in `w0`) to `phi1`:
/// a window where `G` changes sign across the ends for every `phi` of the step
/// and `G_t` keeps a sign; the root at `phi1` narrowed in it.
fn step<T: Real>(
    f: &dyn Chart<T>,
    phi0: f64,
    w0: Enclosure,
    phi1: f64,
) -> Option<(Enclosure, Enclosure)> {
    let tm = mid(w0);
    let j = f.jet(&T::exact_f64(phi0), &T::exact_f64(tm));
    let c = |x: &T| {
        let (lo, hi) = x.bounds_f64();
        mid([lo, hi])
    };
    // The branch's slope and curvature from implicit differentiation:
    // t' = -G_phi / G_t, t'' = -(G_pp + 2 G_pt t' + G_tt t'^2) / G_t.
    let (gt, slope) = (c(&j.gt), -c(&j.gp) / c(&j.gt));
    let bend = -(c(&j.gpp) + 2.0 * c(&j.gpt) * slope + c(&j.gtt) * slope * slope) / gt;
    if !(slope.is_finite() && bend.is_finite()) {
        return None;
    }
    let h = phi1 - phi0;
    let dt = slope * h + 0.5 * bend * h * h;
    let margin = 0.5 * (slope * h).abs()
        + bend.abs() * h * h
        + 4.0 * (w0[1] - w0[0])
        + 1e-12 * (1.0 + tm.abs());
    let w = [
        w0[0].min(w0[0] + dt) - margin,
        w0[1].max(w0[1] + dt) + margin,
    ];
    let ps = span::<T>(phi0, phi1);
    let jb = f.jet(&ps, &span(w[0], w[1]));
    if !certain(&jb.gt) {
        return None;
    }
    // G on the window's ends over the step, by the mean-value form in phi.
    let pm = 0.5 * phi0 + 0.5 * phi1;
    let edge = |t: f64| {
        let tt = T::exact_f64(t);
        let j = f.jet(&ps, &tt);
        sign(
            &f.value(&T::exact_f64(pm), &tt)
                .add(&j.gp.mul(&offset::<T>([phi0, phi1], pm))),
        )
    };
    let (lo, hi) = (edge(w[0]), edge(w[1]));
    if lo.is_none() || hi.is_none() || lo == hi {
        return None;
    }
    // The window is certified at phi1 too: its root narrowed by Newton.
    let root = newton(f, phi1, w).ok()?;
    Some((w, root))
}

/// A branch followed from `phi0` (root `w0`) to `phi1`, adaptively.
fn follow(
    fast: &dyn Chart<Fast>,
    exact: &dyn Chart<I>,
    phi0: f64,
    w0: Enclosure,
    phi1: f64,
    budget: &mut usize,
) -> Result<(Vec<Step>, Enclosure)> {
    let mut steps = Vec::new();
    let (mut at, mut w) = (phi0, w0);
    let mut h = ((phi1 - phi0) / 8.0).max(1e-9);
    while at < phi1 {
        *budget = budget
            .checked_sub(1)
            .ok_or(limit("a torus curve's branches"))?;
        let next = if at + h >= phi1 { phi1 } else { at + h };
        let done = step(fast, at, w, next).or_else(|| {
            if h < 1e-9 {
                step(exact, at, w, next)
            } else {
                None
            }
        });
        match done {
            Some((window, root)) => {
                steps.push(Step {
                    phi: [at, next],
                    window,
                });
                at = next;
                w = root;
                h *= 2.0;
            }
            None => {
                h *= 0.25;
                if h < 1e-13 {
                    return Err(limit("a torus curve's branches"));
                }
            }
        }
    }
    Ok((steps, w))
}

fn overlaps_mod(a: Enclosure, b: Enclosure, period: f64) -> Option<f64> {
    for k in [-2.0, -1.0, 0.0, 1.0, 2.0] {
        let s = k * period;
        if a[0] <= b[1] + s && b[0] + s <= a[1] {
            return Some(s);
        }
    }
    None
}

// ------------------------------------------------------------------ entry

/// `atan2(y, x)` enclosed, continued across `pi` when `x` is certainly
/// negative (an angle at `pi` is `pi + atan2(-y, -x)`).
pub(super) fn angle(y: &I, x: &I) -> Result<I> {
    if x.sign() == Some(Ordering::Less) {
        let a = crate::certified::atan2(&y.neg(), &x.neg()).ok_or(limit("a tangency's angle"))?;
        let p = crate::certified::pi();
        return Ok(a.add(&I::new(p.lo().clone(), p.hi().clone())));
    }
    crate::certified::atan2(y, x).ok_or(limit("a tangency's angle"))
}

/// The intersection of a torus and a cylinder, a cone or another torus
/// (`torus_pair`, the carrier first by stored data) off its axis.
pub(super) fn intersect(torus: &Surface, other: &Surface) -> Result<TracedCurve> {
    let fast = Field::<Fast>::of(torus, other)?;
    let exact = Field::<I>::of(torus, other)?;
    let (Some(td), Some(od)) = (data(torus), data(other)) else {
        return Err(Error::OutOfDomain("a traced torus curve"));
    };
    // Tangencies, exactly (a cone is never exactly tangent).
    let mut nodes = Vec::new();
    let mut zones = Vec::new();
    if !od.cone {
        // A torus pair's tangencies (torus_pair), else a cylinder's.
        let contacts = if od.torus {
            tangency::torus_torus(
                &td.o, &td.a, &td.r, &td.minor, &od.o, &od.a, &od.r, &od.minor,
            )?
        } else {
            tangency::torus_cylinder(&td.o, &td.a, &td.r, &td.minor, &od.o, &od.a, &od.r)?
        };
        for c in contacts {
            let rel = esub(&c.spine, &exact.o);
            let phi = angle(&edot(&rel, &exact.y), &edot(&rel, &exact.x))?;
            let (cp, sp) = I::cos_sin(&phi);
            let e = eadd(&escale(&exact.x, &cp), &escale(&exact.y, &sp));
            let n = esub(&c.point, &c.spine);
            let t = angle(&edot(&n, &exact.a), &edot(&n, &e))?;
            let at = [bounds(&phi), bounds(&t)];
            let j = exact.jet(&span(at[0][0], at[0][1]), &span(at[1][0], at[1][1]));
            let det = j.gpp.mul(&j.gtt).sub(&j.gpt.square());
            let crossing = match det.sign() {
                Some(Ordering::Less) => true,
                Some(Ordering::Greater) => false,
                _ => return Err(limit("a degenerate contact")),
            };
            let local = node_box(&fast, &exact, &at, crossing, nodes.len())?;
            zones.push(local.b);
            nodes.push((
                Node {
                    phi: at[0],
                    t: at[1],
                    point: bounds3(&c.point),
                    crossing,
                },
                local,
            ));
        }
    }
    let locals: Vec<Local> = nodes.iter().map(|(_, l)| l.clone()).collect();
    let nodes: Vec<Node> = nodes.into_iter().map(|(n, _)| n).collect();
    let (folds, tracks, components) = graph(&fast, &exact, TAU, locals, &nodes)?;
    Ok(TracedCurve {
        carrier: torus.clone(),
        other: other.clone(),
        chart: ChartKind::Meridians,
        period: TAU,
        folds,
        nodes,
        tracks,
        components,
    })
}

/// The traced graph on a chart with period `period` in `t`, given the
/// tangencies' boxes: folds and their boxes, tracks, components.
pub(super) fn graph(
    fast: &dyn Chart<Fast>,
    exact: &dyn Chart<I>,
    period: f64,
    mut locals: Vec<Local>,
    nodes: &[Node],
) -> Result<(Vec<Fold>, Vec<Track>, Vec<TracedComponent>)> {
    let zones: Vec<Bx> = locals.iter().map(|l| l.b).collect();
    let found = folds(fast, exact, period, &zones)?;
    let mut fold_rows = Vec::new();
    for k in &found {
        locals.push(fold_box(fast, exact, k)?);
        let p = exact.at(&span(k[0][0], k[0][1]), &span(k[1][0], k[1][1]))?;
        fold_rows.push(Fold {
            phi: k[0],
            t: k[1],
            point: bounds3(&p),
        });
    }
    let (tracks, ends) = sweep(fast, exact, period, &locals)?;
    let components = assemble(&tracks, &ends, &locals, nodes, period);
    Ok((fold_rows, tracks, components))
}

/// A live branch of the sweep.
struct Live {
    from: f64,
    begin: End,
    t0: f64,
    root: Enclosure,
    steps: Vec<Step>,
}

/// Every branch followed once round the torus from a meridian outside every
/// box, stopping at boxes' sides and starting again beyond them.
fn sweep(
    fast: &dyn Chart<Fast>,
    exact: &dyn Chart<I>,
    period: f64,
    locals: &[Local],
) -> Result<(Vec<Track>, Vec<Traced>)> {
    // The start: the middle of the widest gap between the boxes' ranges.
    let mut ranges: Vec<[f64; 2]> = locals.iter().map(|l| l.b[0]).collect();
    ranges.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let start = if ranges.is_empty() {
        0.0
    } else {
        let mut best = (0.0f64, 0.0f64);
        for k in 0..ranges.len() {
            let end = ranges[..=k]
                .iter()
                .fold(f64::NEG_INFINITY, |m, r| m.max(r[1]));
            let next = if k + 1 < ranges.len() {
                ranges[k + 1][0]
            } else {
                ranges[0][0] + TAU
            };
            if next - end > best.0 {
                best = (next - end, 0.5 * end + 0.5 * next);
            }
        }
        if best.0 <= 0.0 {
            return Err(limit("boxes round the whole torus"));
        }
        best.1
    };
    // Box sides as events in [start, start + 2 pi).
    let shift = |x: f64| x + TAU * ((start - x) / TAU).ceil();
    let mut events: Vec<(f64, usize, usize)> = Vec::new();
    for (b, l) in locals.iter().enumerate() {
        let lo = shift(l.b[0][0]);
        let hi = lo + (l.b[0][1] - l.b[0][0]);
        events.push((lo, b, 0));
        events.push((hi, b, 1));
    }
    events.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut budget = 200_000usize;
    let first = {
        let mut rb = 100_000usize;
        let half = above(0.5 * period);
        let roots = match roots_on(fast, start, -half, half, &mut rb) {
            Ok(r) => r,
            Err(Error::ComputationLimit(_)) => roots_on(exact, start, -half, half, &mut rb)?,
            Err(e) => return Err(e),
        };
        dedupe_seam(roots, period)
    };
    let mut live: Vec<Live> = first
        .iter()
        .enumerate()
        .map(|(j, r)| Live {
            from: start,
            begin: End::Seam(j),
            t0: mid(*r),
            root: *r,
            steps: Vec::new(),
        })
        .collect();
    let mut done: Vec<Traced> = Vec::new();
    let mut at = start;
    let finish = |mut l: Live, to: f64, end: End, done: &mut Vec<Traced>| {
        // A branch between two boxes' sides at one meridian: its root there.
        if l.steps.is_empty() {
            l.steps.push(Step {
                phi: [l.from, to],
                window: l.root,
            });
        }
        done.push(Traced {
            track: Track {
                phi: [l.from, to],
                steps: l.steps,
            },
            ends: [l.begin, end],
            t: [l.t0, mid(l.root)],
        });
    };
    for (phi, b, side) in
        events
            .iter()
            .copied()
            .chain(std::iter::once((start + TAU, usize::MAX, 0)))
    {
        // Advance every branch to the event.
        for l in live.iter_mut() {
            if phi > at {
                let (steps, root) = follow(fast, exact, at, l.root, phi, &mut budget)?;
                l.steps.extend(steps);
                l.root = root;
            }
        }
        at = phi;
        if b == usize::MAX {
            break;
        }
        let local = &locals[b];
        let (roots, window) = (
            if side == 0 { &local.left } else { &local.right },
            local.b[1],
        );
        if side == 0 {
            // Branches reaching the box's left side end there.
            for (r, root) in roots.iter().enumerate() {
                let hits: Vec<usize> = live
                    .iter()
                    .enumerate()
                    .filter(|(_, l)| overlaps_mod(l.root, *root, period).is_some())
                    .map(|(k, _)| k)
                    .collect();
                let [k] = hits[..] else {
                    return Err(limit("a branch entering a box"));
                };
                let l = live.remove(k);
                finish(l, phi, End::Side(b, 0, r), &mut done);
            }
        } else {
            for (r, root) in roots.iter().enumerate() {
                live.push(Live {
                    from: phi,
                    begin: End::Side(b, 1, r),
                    t0: mid(*root),
                    root: *root,
                    steps: Vec::new(),
                });
            }
        }
        // No other branch passes through the box's window: its top and
        // bottom carry no zero, so one outside stays outside.
        let own = if side == 0 { 0 } else { roots.len() };
        let within = live
            .iter()
            .filter(|l| overlaps_mod(l.root, window, period).is_some())
            .count();
        if within != own {
            return Err(limit("a branch inside a box"));
        }
    }
    // The branches reaching the end meet the start's roots.
    let end = start + TAU;
    for l in std::mem::take(&mut live) {
        let hits: Vec<usize> = first
            .iter()
            .enumerate()
            .filter(|(_, r)| overlaps_mod(l.root, **r, period).is_some())
            .map(|(j, _)| j)
            .collect();
        let [j] = hits[..] else {
            return Err(limit("a branch closing round the torus"));
        };
        finish(l, end, End::Seam(j), &mut done);
    }
    let tracks = done.iter().map(|d| d.track.clone()).collect();
    Ok((tracks, done))
}

/// Roots over `[-pi, pi]` (binary64 bounds of pi): a root in the sliver
/// beyond pi is the one near -pi, found twice.
fn dedupe_seam(mut r: Vec<Enclosure>, period: f64) -> Vec<Enclosure> {
    if r.len() >= 2 {
        let (first, last) = (r[0], r[r.len() - 1]);
        if last[0] - period <= first[1] + 1e-12 && first[0] + period <= last[1] + 1e-12 {
            r.pop();
        }
    }
    r
}

/// Components from the tracks' ends: a track's end joins another's at a
/// seam, inside a fold's box (its two side roots), or at a tangency.
fn assemble(
    tracks: &[Track],
    ends: &[Traced],
    locals: &[Local],
    nodes: &[Node],
    period: f64,
) -> Vec<TracedComponent> {
    let n = tracks.len();
    let mut partner: Vec<[Option<(usize, usize)>; 2]> = vec![[None, None]; n];
    let mut at_node: Vec<[Option<usize>; 2]> = vec![[None, None]; n];
    let find = |e: End| -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for (k, d) in ends.iter().enumerate() {
            for s in 0..2 {
                if d.ends[s] == e {
                    out.push((k, s));
                }
            }
        }
        out
    };
    for (k, d) in ends.iter().enumerate() {
        for s in 0..2 {
            match d.ends[s] {
                End::Seam(j) => {
                    let other = 1 - s;
                    if let Some(&(m, ms)) = find(End::Seam(j))
                        .iter()
                        .find(|(m, ms)| (*m, *ms) != (k, s) && *ms == other)
                    {
                        partner[k][s] = Some((m, ms));
                    }
                }
                End::Side(b, side, r) => match locals[b].kind {
                    LocalKind::Fold => {
                        let e = End::Side(b, side, 1 - r);
                        partner[k][s] = find(e).first().copied();
                    }
                    LocalKind::Node(node) => at_node[k][s] = Some(node),
                },
            }
        }
    }
    let mut out = Vec::new();
    for (k, node) in nodes.iter().enumerate() {
        if !node.crossing {
            out.push(TracedComponent::Isolated { node: k });
        }
    }
    // Union through partners and nodes.
    let mut group: Vec<usize> = (0..n).collect();
    fn root(g: &mut [usize], mut x: usize) -> usize {
        while g[x] != x {
            g[x] = g[g[x]];
            x = g[x];
        }
        x
    }
    for (k, pair) in partner.iter().enumerate() {
        for (m, _) in pair.iter().flatten() {
            let (a, b) = (root(&mut group, k), root(&mut group, *m));
            group[a] = b;
        }
    }
    let mut node_owner: Vec<Option<usize>> = vec![None; nodes.len()];
    for (k, pair) in at_node.iter().enumerate() {
        for v in pair.iter().flatten().copied() {
            {
                match node_owner[v] {
                    Some(m) => {
                        let (a, b) = (root(&mut group, k), root(&mut group, m));
                        group[a] = b;
                    }
                    None => node_owner[v] = Some(k),
                }
            }
        }
    }
    let mut seen = vec![false; n];
    for k in 0..n {
        if seen[k] {
            continue;
        }
        let g = root(&mut group, k);
        let members: Vec<usize> = (0..n).filter(|m| root(&mut group, *m) == g).collect();
        for m in &members {
            seen[*m] = true;
        }
        let folds = members
            .iter()
            .map(|m| {
                (0..2)
                    .filter(|s| {
                        matches!(ends[*m].ends[*s], End::Side(b, _, _)
                            if matches!(locals[b].kind, LocalKind::Fold))
                    })
                    .count()
            })
            .sum::<usize>()
            / 2;
        let node_set: Vec<usize> = {
            let mut v: Vec<usize> = members
                .iter()
                .flat_map(|m| at_node[*m].iter().flatten().copied().collect::<Vec<_>>())
                .collect();
            v.sort_unstable();
            v.dedup();
            v
        };
        if !node_set.is_empty() {
            out.push(TracedComponent::Crossing {
                tracks: members,
                folds,
                nodes: node_set,
                infinite: 0,
            });
            continue;
        }
        // Winding: walk the cycle, the joints adding only their residue.
        let mut wp = 0.0;
        let mut wt = 0.0;
        let mut order = Vec::new();
        let (mut cur, mut fwd) = (members[0], true);
        loop {
            order.push(cur);
            let d = &ends[cur];
            let (a, b) = if fwd { (0, 1) } else { (1, 0) };
            wp += if fwd { 1.0 } else { -1.0 } * (d.track.phi[1] - d.track.phi[0]);
            wt += d.t[b] - d.t[a];
            let Some((next, ns)) = partner[cur][b] else {
                break;
            };
            let jump = ends[next].t[ns] - d.t[b];
            wt += jump - period * (jump / period).round();
            cur = next;
            fwd = ns == 0;
            if (cur == members[0] && fwd) || order.len() > members.len() {
                break;
            }
        }
        let mut w = [(wp / TAU).round() as i64, (wt / period).round() as i64];
        if w[0] < 0 || (w[0] == 0 && w[1] < 0) {
            w = [-w[0], -w[1]];
        }
        out.push(TracedComponent::Smooth {
            tracks: order,
            folds,
            winding: w,
            infinite: 0,
        });
    }
    out
}
