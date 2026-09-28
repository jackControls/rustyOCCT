//! A torus and a cylinder or a cone off its axis (S7b.3b.1 of
//! `REVIEW_NOTES.md`): traced curves on the torus's meridians.
//!
//! The torus (axis `a` normalised exactly, major radius `R`, minor `r`) is
//! `p(phi, t) = o + (R + r cos t) e(phi) + r sin t a`, `e(phi) = cos phi x +
//! sin phi y`, `x` the other axis's component normal to `a` (the direction to
//! the other's origin when the axes are parallel). The curve is the zero set
//! of `G(phi, t) = f(p(phi, t))` on the flat parameter torus, `f` the other
//! surface's implicit function `(p - q)^T Q (p - q) - k`, `Q = c2 - u u^T`
//! (a cylinder: `q` its origin, `c2 = 1`, `k = r_c^2`; a cone, both nappes:
//! `q` its apex, `c2 = cos^2 h`, `k = 0`; `u` the unit axis).
//!
//! The curve is a graph. Its vertices are the tangencies of the surfaces,
//! decided exactly (`tangency.rs`); its folds, where a component turns in
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
    },
    /// Branches through crossings: their tracks, folds and nodes.
    Crossing {
        tracks: Vec<usize>,
        folds: usize,
        nodes: Vec<usize>,
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

/// A torus's curve with a cylinder or a cone, as a graph of tracks, folds
/// and tangencies (D13's procedural curve for S7b.3b).
#[derive(Debug, Clone, PartialEq)]
pub struct TracedCurve {
    carrier: Surface,
    other: Surface,
    folds: Vec<Fold>,
    nodes: Vec<Node>,
    tracks: Vec<Track>,
    components: Vec<TracedComponent>,
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
        let fast = Field::<Fast>::of(&self.carrier, &self.other)?;
        let narrow = |p: &Enclosure3| {
            p.iter()
                .all(|[a, b]| b - a <= 1e-12 * a.abs().max(b.abs()).max(1.0))
        };
        // Binary64 intervals resolve t only to about their rounding over
        // |G_t|; rational ones narrow further from their bracket.
        let bracket = match root_in(&fast, at, step.window) {
            Ok(t) => match point(&fast, at, t) {
                Ok(p) if narrow(&p) => return Ok((t, p)),
                _ => t,
            },
            Err(_) => step.window,
        };
        let exact = Field::<I>::of(&self.carrier, &self.other)?;
        let t = newton(&exact, at, bracket)?;
        Ok((t, point(&exact, at, t)?))
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
}

fn data(s: &Surface) -> Option<Data> {
    let (f, r, minor, angle, cone) = match s {
        Surface::Torus {
            frame,
            major,
            minor,
        } => (frame, *major, *minor, 0.0, false),
        Surface::Cylinder { frame, radius } => (frame, *radius, 0.0, 0.0, false),
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => (frame, *radius, 0.0, *half_angle, true),
        _ => return None,
    };
    Some(Data {
        o: f.origin().to_array().map(q),
        a: f.normal().to_array().map(q),
        r: q(r),
        minor: q(minor),
        angle,
        cone,
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
    /// `f(p) = (p - q)^T Q (p - q) - k`, `Q = c2 - u u^T`.
    q: E<T>,
    u: E<T>,
    c2: T,
    k: T,
}

/// `G` and its derivatives up to the second order.
struct Jet<T> {
    g: T,
    gp: T,
    gt: T,
    gpp: T,
    gpt: T,
    gtt: T,
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
        self.form(&w, &w).sub(&self.k)
    }
    fn jet(&self, phi: &T, t: &T) -> Jet<T> {
        let [p, pp, pt, ppp, ppt, ptt] = self.frame(phi, t);
        let w = esub(&p, &self.q);
        let qw = self.qv(&w);
        let two = T::exact_f64(2.0);
        let d = |v: &E<T>| edot(&qw, v).mul(&two);
        let h = |v: &E<T>, z: &E<T>| self.form(v, z).mul(&two);
        Jet {
            g: edot(&w, &qw).sub(&self.k),
            gp: d(&pp),
            gt: d(&pt),
            gpp: h(&pp, &pp).add(&d(&ppp)),
            gpt: h(&pp, &pt).add(&d(&ppt)),
            gtt: h(&pt, &pt).add(&d(&ptt)),
        }
    }
}

fn point<T: Real>(f: &Field<T>, phi: f64, t: [f64; 2]) -> Result<Enclosure3> {
    let [p, ..] = f.frame(&T::exact_f64(phi), &span(t[0], t[1]));
    Ok(bounds3(&p))
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
    f: &Field<T>,
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
        if sa.is_some() && sb.is_some() && sa != sb && certain(&j.gt) {
            out.push(newton(f, phi, [a, b])?);
            continue;
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
fn root_in<T: Real>(f: &Field<T>, phi: f64, w: [f64; 2]) -> Result<Enclosure> {
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
fn newton<T: Real>(f: &Field<T>, phi: f64, w: [f64; 2]) -> Result<Enclosure> {
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

/// No zero of `G` on `[p0, p1] x {t}`: subdivision with the mean-value
/// enclosure in `phi`.
fn edge_free<T: Real>(f: &Field<T>, p0: f64, p1: f64, t: f64, budget: &mut usize) -> Result<bool> {
    let tt = T::exact_f64(t);
    let mut pending = vec![(p0, p1, 40usize)];
    while let Some((a, b, depth)) = pending.pop() {
        *budget = budget.checked_sub(1).ok_or(limit("a box's edge"))?;
        let j = f.jet(&span(a, b), &tt);
        if certain(&j.g) {
            continue;
        }
        let m = 0.5 * a + 0.5 * b;
        let mv = f
            .value(&T::exact_f64(m), &tt)
            .add(&j.gp.mul(&offset::<T>([a, b], m)));
        if certain(&mv) {
            continue;
        }
        if depth == 0 || !(a < m && m < b) {
            return Ok(false);
        }
        pending.push((m, b, depth - 1));
        pending.push((a, m, depth - 1));
    }
    Ok(true)
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
fn krawczyk<T: Real>(f: &Field<T>, b: &Bx, critical: bool) -> Kraw {
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
fn refine<T: Real>(f: &Field<T>, mut b: Bx, critical: bool) -> Bx {
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

fn verdict<T: Real>(f: &Field<T>, b: &Bx) -> Verdict {
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
fn folds(fast: &Field<Fast>, exact: &Field<I>, zones: &[Bx]) -> Result<Vec<Bx>> {
    let mut out: Vec<Bx> = Vec::new();
    let mut pending = vec![([[-PI_HI, PI_HI], [-PI_HI, PI_HI]], 0usize)];
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
                [-TAU, 0.0, TAU]
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
                        [-TAU, 0.0, TAU]
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
struct Local {
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

fn local_edges<T: Real>(f: &Field<T>, b: &Bx) -> Result<Option<Sides>> {
    let mut budget = 20_000;
    for t in [b[1][0], b[1][1]] {
        if !edge_free(f, b[0][0], b[0][1], t, &mut budget)? {
            return Ok(None);
        }
    }
    let left = roots_on(f, b[0][0], b[1][0], b[1][1], &mut budget);
    let right = roots_on(f, b[0][1], b[1][0], b[1][1], &mut budget);
    match (left, right) {
        (Ok(l), Ok(r)) => Ok(Some((l, r))),
        (Err(Error::ComputationLimit(_)), _) | (_, Err(Error::ComputationLimit(_))) => Ok(None),
        (Err(e), _) | (_, Err(e)) => Err(e),
    }
}

/// The box of a fold at `k`: `G_phi` of one sign on it, no zero on its top
/// and bottom, two roots on one side and none on the other.
fn fold_box(fast: &Field<Fast>, exact: &Field<I>, k: &Bx) -> Result<Local> {
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
            let edges = match local_edges(fast, &b)? {
                Some(e) => Some(e),
                None => local_edges(exact, &b)?,
            };
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
    fast: &Field<Fast>,
    exact: &Field<I>,
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
            let edges = match local_edges(fast, &b)? {
                Some(e) => Some(e),
                None => local_edges(exact, &b)?,
            };
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
    f: &Field<T>,
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
    let slope = -c(&j.gp) / c(&j.gt);
    if !slope.is_finite() {
        return None;
    }
    let dt = slope * (phi1 - phi0);
    let margin = 1.5 * dt.abs() + 4.0 * (w0[1] - w0[0]) + 1e-12 * (1.0 + tm.abs());
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
    fast: &Field<Fast>,
    exact: &Field<I>,
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
            if h < 1e-6 {
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

fn overlaps_mod(a: Enclosure, b: Enclosure) -> Option<f64> {
    for k in [-2.0, -1.0, 0.0, 1.0, 2.0] {
        let s = k * TAU;
        if a[0] <= b[1] + s && b[0] + s <= a[1] {
            return Some(s);
        }
    }
    None
}

// ------------------------------------------------------------------ entry

/// `atan2(y, x)` enclosed, continued across `pi` when `x` is certainly
/// negative (an angle at `pi` is `pi + atan2(-y, -x)`).
fn angle(y: &I, x: &I) -> Result<I> {
    if x.sign() == Some(Ordering::Less) {
        let a = crate::certified::atan2(&y.neg(), &x.neg()).ok_or(limit("a tangency's angle"))?;
        let p = crate::certified::pi();
        return Ok(a.add(&I::new(p.lo().clone(), p.hi().clone())));
    }
    crate::certified::atan2(y, x).ok_or(limit("a tangency's angle"))
}

/// The intersection of a torus and a cylinder or cone off its axis.
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
        for c in tangency::torus_cylinder(&td.o, &td.a, &td.r, &td.minor, &od.o, &od.a, &od.r)? {
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
    let found = folds(&fast, &exact, &zones)?;
    let mut locals: Vec<Local> = nodes.iter().map(|(_, l)| l.clone()).collect();
    let mut fold_rows = Vec::new();
    for k in &found {
        locals.push(fold_box(&fast, &exact, k)?);
        let [p, ..] = exact.frame(&span(k[0][0], k[0][1]), &span(k[1][0], k[1][1]));
        fold_rows.push(Fold {
            phi: k[0],
            t: k[1],
            point: bounds3(&p),
        });
    }
    let nodes: Vec<Node> = nodes.into_iter().map(|(n, _)| n).collect();
    let (tracks, ends) = sweep(&fast, &exact, &locals)?;
    let components = assemble(&tracks, &ends, &locals, &nodes);
    Ok(TracedCurve {
        carrier: torus.clone(),
        other: other.clone(),
        folds: fold_rows,
        nodes,
        tracks,
        components,
    })
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
    fast: &Field<Fast>,
    exact: &Field<I>,
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
        let roots = match roots_on(fast, start, -PI_HI, PI_HI, &mut rb) {
            Ok(r) => r,
            Err(Error::ComputationLimit(_)) => roots_on(exact, start, -PI_HI, PI_HI, &mut rb)?,
            Err(e) => return Err(e),
        };
        dedupe_seam(roots)
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
                    .filter(|(_, l)| overlaps_mod(l.root, *root).is_some())
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
            .filter(|l| overlaps_mod(l.root, window).is_some())
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
            .filter(|(_, r)| overlaps_mod(l.root, **r).is_some())
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
fn dedupe_seam(mut r: Vec<Enclosure>) -> Vec<Enclosure> {
    if r.len() >= 2 {
        let (first, last) = (r[0], r[r.len() - 1]);
        if last[0] - TAU <= first[1] + 1e-12 && first[0] + TAU <= last[1] + 1e-12 {
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
            wt += jump - TAU * (jump / TAU).round();
            cur = next;
            fwd = ns == 0;
            if (cur == members[0] && fwd) || order.len() > members.len() {
                break;
            }
        }
        let mut w = [(wp / TAU).round() as i64, (wt / TAU).round() as i64];
        if w[0] < 0 || (w[0] == 0 && w[1] < 0) {
            w = [-w[0], -w[1]];
        }
        out.push(TracedComponent::Smooth {
            tracks: order,
            folds,
            winding: w,
        });
    }
    out
}
