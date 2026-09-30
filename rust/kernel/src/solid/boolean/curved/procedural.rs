//! S9c.2a: two cylinders meeting in quartics (REVIEW_NOTES.md, S9c.2):
//! perpendicular axes in exact frames, of any radii and offset. The curve
//! is cut into graphs over one cylinder's angle (S7b's ruled
//! parameterisation, D13), each clear of its turning points: two rings
//! over the cylinder whose extent across the common perpendicular lies
//! inside the other's, or one loop in four graphs switched at rational
//! points of the thinner cylinder's angle. Where such a section crosses a
//! cap's circle the vertex is algebraic (`algebraic.rs`, S9c.2b.2).
use super::graph::{between_ccw, rational_between, same_dir};
use super::meet::{tangency, CylPair};
use super::model::*;
use super::num::*;
use crate::solid::split::{rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// A piece of two quadrics' meeting: the carrier's ruling at `(cos, sin)`,
/// `o + r (cos x + sin y) + w (n + k (cos x + sin y))` (a cylinder's `k`
/// zero, a cone's its slope: S9d.3b), meets the other quadric `sum_i (g_i
/// . p - e_i)^2 = (r2 + t (h . p - e_h))^2` at `w = (-B + s sqrt(B^2 - A
/// C)) / A`, `s` the branch (`plus`); an open piece runs counter-clockwise
/// over `range`, a ring (`None`) over a whole turn.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct MeetCrv {
    /// The carrier's operand.
    pub(super) carrier: usize,
    o: V,
    x: V,
    y: V,
    n: V,
    /// The carrier's local rows for `u`, `v` and `w`.
    k: [V; 3],
    r: R,
    /// The carrier's slope (zero for a cylinder).
    slope: R,
    g: Vec<V>,
    e: Vec<R>,
    r2: R,
    /// The other's cone term: `h`, `e_h` and `t` (zero for a cylinder or a
    /// sphere).
    h: V,
    eh: R,
    t: R,
    pub(super) plus: bool,
    pub(super) range: Option<[[Qd; 2]; 2]>,
}

/// The other quadric of a carrier, as its local rows and offsets: `sum_i
/// (g_i . p - e_i)^2 - (r + t (h . p - e_h))^2`, zero on it (a cone's `t`
/// its slope, zero for a cylinder or a sphere).
#[derive(Debug, Clone)]
pub(super) struct Other {
    pub(super) g: Vec<V>,
    pub(super) e: Vec<R>,
    pub(super) r: R,
    pub(super) h: V,
    pub(super) eh: R,
    pub(super) t: R,
}

impl Other {
    /// The radius term `r + t (h . p - e_h)` at a point.
    pub(super) fn radius_at(&self, x: &QV) -> Qd {
        qdot(x, &self.h)
            .add_r(&-self.eh.clone())
            .scale(&self.t)
            .add_r(&self.r)
    }

    /// The quadric's value at a point.
    pub(super) fn value(&self, x: &QV) -> Qd {
        let rho = self.radius_at(x);
        self.g
            .iter()
            .zip(&self.e)
            .fold(Qd::rat(zero()), |acc, (g, e)| {
                let s = qdot(x, g).add_r(&-e.clone());
                acc.add(&s.mul(&s))
            })
            .sub(&rho.mul(&rho))
    }

    /// The quadric's gradient at a point (halved).
    pub(super) fn gradient(&self, x: &QV) -> QV {
        let mut out = qscale(&self.h, &self.radius_at(x).scale(&-self.t.clone()));
        for (g, e) in self.g.iter().zip(&self.e) {
            out = qadd(&out, &qscale(g, &qdot(x, g).add_r(&-e.clone())));
        }
        out
    }

    /// The cone term as a linear form along `c + a cos + b sin` (its
    /// square subtracts).
    pub(super) fn radius_lin(&self, c: &V, a: &V, b: &V) -> [R; 3] {
        [
            &self.r + &self.t * (dot(&self.h, c) - &self.eh),
            &self.t * dot(&self.h, a),
            &self.t * dot(&self.h, b),
        ]
    }
}

pub(super) fn other_of(f: &Affine, c: &P2, r: &R) -> Other {
    let g = vec![f.row(0).clone(), f.row(1).clone()];
    let e = vec![&dot(&g[0], &f.o) + &c[0], &dot(&g[1], &f.o) + &c[1]];
    Other {
        g,
        e,
        r: r.clone(),
        h: [zero(), zero(), zero()],
        eh: zero(),
        t: zero(),
    }
}

/// A sphere as the other quadric (S9d.2): `|p - c|^2 = r^2` over the world's
/// rows.
pub(super) fn other_sphere(c: &V, r: &R) -> Other {
    let rows = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(|e| e.map(int));
    Other {
        g: rows.to_vec(),
        e: c.to_vec(),
        r: r.clone(),
        h: [zero(), zero(), zero()],
        eh: zero(),
        t: zero(),
    }
}

/// A cone as the other quadric (S9d.3b): `u^2 + v^2 = (b + k w)^2` in its
/// frame's coordinates.
pub(super) fn other_cone(f: &Affine, b: &R, k: &R) -> Other {
    let g = vec![f.row(0).clone(), f.row(1).clone()];
    let e = vec![dot(&g[0], &f.o), dot(&g[1], &f.o)];
    Other {
        g,
        e,
        r: b.clone(),
        h: f.row(2).clone(),
        eh: dot(f.row(2), &f.o),
        t: k.clone(),
    }
}

/// A ruled carrier: a cylinder's frame, circle centre and radius (slope
/// zero), or a cone's frame, its axis at `(0, 0)`, its radius at `w = 0`
/// and its slope (S9d.3b).
#[derive(Debug, Clone, Copy)]
pub(super) struct Ruled<'a> {
    pub(super) f: &'a Affine,
    pub(super) c: &'a P2,
    pub(super) r: &'a R,
    pub(super) k: &'a R,
}

impl MeetCrv {
    pub(super) fn new(
        carrier: usize,
        (f, c, r): (&Affine, &P2, &R),
        other: &Other,
        plus: bool,
        range: Option<[[Qd; 2]; 2]>,
    ) -> Self {
        let k = zero();
        Self::ruled(carrier, Ruled { f, c, r, k: &k }, other, plus, range)
    }

    /// A piece over a ruled carrier (a cylinder or a cone, S9d.3b).
    pub(super) fn ruled(
        carrier: usize,
        ruled: Ruled,
        other: &Other,
        plus: bool,
        range: Option<[[Qd; 2]; 2]>,
    ) -> Self {
        let Ruled { f, c, r, k } = ruled;
        Self {
            carrier,
            o: f.point(&c[0], &c[1], &zero()),
            x: f.x.clone(),
            y: f.y.clone(),
            n: f.n.clone(),
            k: [f.row(0).clone(), f.row(1).clone(), f.row(2).clone()],
            r: r.clone(),
            slope: k.clone(),
            g: other.g.clone(),
            e: other.e.clone(),
            r2: other.r.clone(),
            h: other.h.clone(),
            eh: other.eh.clone(),
            t: other.t.clone(),
            plus,
            range,
        }
    }

    fn other(&self) -> Other {
        Other {
            g: self.g.clone(),
            e: self.e.clone(),
            r: self.r2.clone(),
            h: self.h.clone(),
            eh: self.eh.clone(),
            t: self.t.clone(),
        }
    }

    /// The ruling's direction at a `(cos, sin)`.
    fn dir(&self, cs: &[Qd; 2]) -> QV {
        let radial = qadd(&qscale(&self.x, &cs[0]), &qscale(&self.y, &cs[1]));
        qadd(&qv(&self.n), &radial.map(|c| c.scale(&self.slope)))
    }

    /// The carrier's circle point (height 0) at a rational `(cos, sin)`.
    fn base(&self, cs: &[R; 2]) -> V {
        add(
            &self.o,
            &add(
                &scale(&self.x, &(&self.r * &cs[0])),
                &scale(&self.y, &(&self.r * &cs[1])),
            ),
        )
    }

    /// The piece's point at a rational `(cos, sin)` (`None` where the ruling
    /// misses or touches the other cylinder).
    pub(super) fn at(&self, cs: &[R; 2]) -> Option<QV> {
        let base = self.base(cs);
        let dir = add(
            &self.n,
            &scale(
                &add(&scale(&self.x, &cs[0]), &scale(&self.y, &cs[1])),
                &self.slope,
            ),
        );
        let gn: Vec<R> = self.g.iter().map(|g| dot(g, &dir)).collect();
        let s: Vec<R> = self
            .g
            .iter()
            .zip(&self.e)
            .map(|(g, e)| dot(g, &base) - e)
            .collect();
        // The other's radius along the ruling: `r0 + w rd`.
        let r0 = &self.r2 + &self.t * (dot(&self.h, &base) - &self.eh);
        let rd = &self.t * dot(&self.h, &dir);
        let a = gn.iter().fold(zero(), |acc, x| acc + x * x) - &rd * &rd;
        let b = s.iter().zip(&gn).fold(zero(), |acc, (x, y)| acc + x * y) - &r0 * &rd;
        let c = s.iter().fold(zero(), |acc, x| acc + x * x) - &r0 * &r0;
        let d = &b * &b - &a * &c;
        if d <= zero() {
            return None;
        }
        let k = if self.plus { int(1) } else { int(-1) };
        let w = if a == zero() {
            // A ruling along the other's asymptotic direction (S9d.3b.2):
            // its one finite root `-C / 2B`, on the branch of `B`'s sign.
            if sign(&b) != sign(&k) {
                return None;
            }
            Qd::rat(-&c / (int(2) * &b))
        } else {
            Qd::new(-&b / &a, k / &a, d)
        };
        Some(qadd(&qv(&base), &qscale(&dir, &w)))
    }

    /// The carrier's `(cos, sin)` of a point on it: its local `(u, v)` over
    /// the ruling's radius there (a cone's `r + k w`).
    pub(super) fn place(&self, x: &QV) -> [Qd; 2] {
        let d = qsub(x, &qv(&self.o));
        let (u, v) = (qdot(&d, &self.k[0]), qdot(&d, &self.k[1]));
        if self.slope == zero() {
            let inv = int(1) / &self.r;
            return [u.scale(&inv), v.scale(&inv)];
        }
        let rho = qdot(&d, &self.k[2]).scale(&self.slope).add_r(&self.r);
        let inv = rho.recip().expect("a meeting's point off the apex");
        [u.mul(&inv), v.mul(&inv)]
    }

    fn grad_other(&self, x: &QV) -> QV {
        self.other().gradient(x)
    }

    /// The branch a point of both quadrics lies on: the sign of the other's
    /// gradient along the carrier's ruling there.
    pub(super) fn branch_sign(&self, x: &QV) -> Ordering {
        let dir = self.dir(&self.place(x));
        qqdot(&self.grad_other(x), &dir).sign()
    }

    /// Whether a point lies on the piece (on both cylinders, on its branch
    /// and within its run, ends included).
    pub(super) fn on(&self, x: &QV) -> bool {
        // A cone carrier's apex height, where its circle is the apex alone
        // (never a meeting's point: an apex on the other is refused before;
        // S9d.3c's replays: a sphere's pole there).
        if self.slope != zero() {
            let d = qsub(x, &qv(&self.o));
            let rho = qdot(&d, &self.k[2]).scale(&self.slope).add_r(&self.r);
            if rho.sign() == Ordering::Equal {
                return false;
            }
        }
        let cs = self.place(x);
        let unit = cs[0].mul(&cs[0]).add(&cs[1].mul(&cs[1])).add_r(&int(-1));
        if unit.sign() != Ordering::Equal {
            return false;
        }
        if self.other().value(x).sign() != Ordering::Equal {
            return false;
        }
        let branch = qqdot(&self.grad_other(x), &self.dir(&cs)).sign();
        let want = if self.plus {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        if branch != want {
            return false;
        }
        match &self.range {
            None => true,
            Some([lo, hi]) => same_dir(&cs, lo) || same_dir(&cs, hi) || between_ccw(lo, &cs, hi),
        }
    }

    /// The unit-free tangent at a point of the piece, running with its
    /// angle: the two gradients' cross product, turned to the carrier's
    /// counter-clockwise run.
    pub(super) fn tangent(&self, cs: &[Qd; 2], x: &QV) -> QV {
        // The carrier's normal: `cos u_row + sin v_row - k w_row` (a cone's
        // gradient over its radius there).
        let gk = qadd(
            &qadd(&qscale(&self.k[0], &cs[0]), &qscale(&self.k[1], &cs[1])),
            &qv(&scale(&self.k[2], &-self.slope.clone())),
        );
        let t = qcross(&gk, &self.grad_other(x));
        let (du, dv) = (qdot(&t, &self.k[0]), qdot(&t, &self.k[1]));
        let run = cs[1].neg().mul(&du).add(&cs[0].mul(&dv)).sign();
        if run == Ordering::Less {
            t.map(|c| c.neg())
        } else {
            t
        }
    }

    /// Binary64 points over the angle from `t0` turning `sweep`.
    pub(super) fn samples(&self, t0: f64, sweep: f64, n: usize) -> Vec<[f64; 3]> {
        let f = |v: &V| v.clone().map(|y| rational_f64(&y));
        let (o, x, y, nn) = (f(&self.o), f(&self.x), f(&self.y), f(&self.n));
        let g: Vec<[f64; 3]> = self.g.iter().map(f).collect();
        let e: Vec<f64> = self.e.iter().map(rational_f64).collect();
        let (r, r2) = (rational_f64(&self.r), rational_f64(&self.r2));
        let (slope, h, eh, tt) = (
            rational_f64(&self.slope),
            f(&self.h),
            rational_f64(&self.eh),
            rational_f64(&self.t),
        );
        let d3 = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let sign = if self.plus { 1.0 } else { -1.0 };
        (0..=n)
            .map(|i| {
                let t = t0 + sweep * i as f64 / n as f64;
                let radial = [0, 1, 2].map(|j| t.cos() * x[j] + t.sin() * y[j]);
                let base = [0, 1, 2].map(|j| o[j] + r * radial[j]);
                let dir = [0, 1, 2].map(|j| nn[j] + slope * radial[j]);
                let gn: Vec<f64> = g.iter().map(|gi| d3(gi, &dir)).collect();
                let s: Vec<f64> = g
                    .iter()
                    .zip(&e)
                    .map(|(gi, ei)| d3(gi, &base) - ei)
                    .collect();
                let (r0, rd) = (r2 + tt * (d3(&h, &base) - eh), tt * d3(&h, &dir));
                let a: f64 = gn.iter().map(|x| x * x).sum::<f64>() - rd * rd;
                let b: f64 = s.iter().zip(&gn).map(|(x, y)| x * y).sum::<f64>() - r0 * rd;
                let c: f64 = s.iter().map(|x| x * x).sum::<f64>() - r0 * r0;
                let sq = sign * (b * b - a * c).max(0.0).sqrt();
                let (p, m) = (-b + sq, -b - sq);
                let w = if p.abs() >= m.abs() { p / a } else { c / m };
                [0, 1, 2].map(|j| base[j] + w * dir[j])
            })
            .collect()
    }
}

/// Two cylinders' quartic meeting: its pieces and the points where a loop
/// switches between them (perpendicular in exact frames, S9c.2a, or in
/// turned frames, S9c.2b.1).
#[derive(Debug, Clone)]
pub(super) struct Quartic {
    pub(super) pieces: Vec<MeetCrv>,
    pub(super) switches: Vec<QV>,
}

/// A cylinder of an operand: its frame, circle centre and radius.
type Cyl<'a> = (&'a Affine, &'a P2, &'a R);

/// The meeting of two circular cylinders with perpendicular axes (`x` of
/// operand 0, `y` of operand 1), classified exactly by their extents
/// across the common perpendicular.
pub(super) fn perpendicular(x: Cyl, y: Cyl) -> Result<CylPair> {
    let cyl = [x, y];
    let e = cross(&x.0.n, &y.0.n);
    let axis = |z: Cyl| z.0.point(&z.1[0], &z.1[1], &zero());
    let mid = [dot(&axis(x), &e), dot(&axis(y), &e)];
    let rad = [x.2.clone(), y.2.clone()];
    let lo = [&mid[0] - &rad[0], &mid[1] - &rad[1]];
    let hi = [&mid[0] + &rad[0], &mid[1] + &rad[1]];
    // Apart, or touching outside.
    match (lo[1].cmp(&hi[0]), lo[0].cmp(&hi[1])) {
        (Ordering::Greater, _) | (_, Ordering::Greater) => return Ok(CylPair::Apart),
        (Ordering::Equal, _) | (_, Ordering::Equal) => return Err(tangency()),
        _ => {}
    }
    if lo[0] == lo[1] || hi[0] == hi[1] {
        // A node (an inside tangency), or S9c.1's ellipses when both ends
        // agree (handled before).
        return Err(tangency());
    }
    let others = [other_of(y.0, y.1, y.2), other_of(x.0, x.1, x.2)];
    let piece = |k: usize, plus: bool, range: Option<[[Qd; 2]; 2]>| {
        MeetCrv::new(k, cyl[k], &others[k], plus, range)
    };
    // Rings about the cylinder whose extent lies inside the other's.
    for k in 0..2 {
        if lo[k] > lo[1 - k] && hi[k] < hi[1 - k] {
            return Ok(CylPair::Quartic(Box::new(Quartic {
                pieces: vec![piece(k, true, None), piece(k, false, None)],
                switches: Vec::new(),
            })));
        }
    }
    // One loop: at the lower end of the overlap cylinder `zl`'s extreme,
    // at the upper cylinder `zu`'s.
    let zl = if lo[0] > lo[1] { 0 } else { 1 };
    let zu = 1 - zl;
    let ends = [(zl, -1i64), (zu, 1i64)];
    // The direction of cylinder `k`'s extreme (`side` = +-1 along e) in
    // its own (cos, sin).
    let dir = |k: usize, side: i64| -> [Qd; 2] {
        let f = cyl[k].0;
        let s = int(side);
        [
            Qd::rat(dot(f.row(0), &e) * &s),
            Qd::rat(dot(f.row(1), &e) * &s),
        ]
    };
    // The turning point of cylinder `w`'s end on the other `k`'s circle,
    // on the side `sw` of `k`'s coordinate along `w`'s axis.
    let turn_on = |k: usize, w: usize, side: i64, sw: bool| -> Result<[Qd; 2]> {
        let eta = if side < 0 { &lo[w] } else { &hi[w] };
        let off = eta - &mid[k];
        let q2 = &rad[k] * &rad[k] - &off * &off;
        if q2 <= zero() {
            return Err(tangency());
        }
        let nw = &cyl[w].0.n;
        let sgn = if sw { int(1) } else { int(-1) };
        let q = Qd::new(zero(), sgn, q2);
        // v = q nw + off e in k's local directions over r.
        let f = cyl[k].0;
        let inv = int(1) / &rad[k];
        let comp = |row: &V| {
            q.scale(&dot(row, nw))
                .add_r(&(&off * dot(row, &e)))
                .scale(&inv)
        };
        Ok([comp(f.row(0)), comp(f.row(1))])
    };
    // The switch points on the thinner cylinder (the upper one when equal).
    let kk = if rad[0] < rad[1] { 0 } else { 1 };
    let ww = 1 - kk;
    let (kside, wside) = (
        ends.iter().find(|x| x.0 == kk).expect("an end").1,
        ends.iter().find(|x| x.0 == ww).expect("an end").1,
    );
    let ext = dir(kk, kside);
    let opposite = dir(kk, -kside);
    // switch[(sk, sw)]: on kk's circle between its extreme and the turning
    // point of ww's end on side sw, on kk's branch sk.
    let mut sw_pts: Vec<((bool, bool), QV, [Qd; 2])> = Vec::new();
    for sw in [false, true] {
        let t = turn_on(kk, ww, wside, sw)?;
        let ccw = !between_ccw(&t, &opposite, &ext);
        let cs = rational_between(&t, &ext, ccw)?;
        let place = [Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())];
        for sk in [false, true] {
            let p = piece(kk, sk, None)
                .at(&cs)
                .ok_or(Error::ComputationLimit("a switch point off its piece"))?;
            sw_pts.push(((sk, sw), p, place.clone()));
        }
    }
    let find = |sk: bool, sw: bool| {
        sw_pts
            .iter()
            .find(|x| x.0 == (sk, sw))
            .expect("a switch point")
    };
    let run_through = |a: [Qd; 2], b: [Qd; 2], through: &[Qd; 2]| -> [[Qd; 2]; 2] {
        if between_ccw(&a, through, &b) {
            [a, b]
        } else {
            [b, a]
        }
    };
    let mut pieces = Vec::new();
    // kk's pieces about its extreme, one per branch sk.
    for sk in [false, true] {
        let range = run_through(find(sk, false).2.clone(), find(sk, true).2.clone(), &ext);
        pieces.push(piece(kk, sk, Some(range)));
    }
    // ww's pieces about its extreme, one per branch sw.
    let wext = dir(ww, wside);
    for sw in [false, true] {
        let probe = piece(ww, sw, None);
        let (a, b) = (
            probe.place(&find(false, sw).1),
            probe.place(&find(true, sw).1),
        );
        let range = run_through(a, b, &wext);
        pieces.push(piece(ww, sw, Some(range)));
    }
    let switches = sw_pts.into_iter().map(|x| x.1).collect();
    Ok(CylPair::Quartic(Box::new(Quartic { pieces, switches })))
}
