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

/// A piece of two cylinders' meeting: the carrier's ruling at `(cos, sin)`,
/// `o + r (cos x + sin y) + w n`, meets the other cylinder
/// `sum_i (g_i . p - e_i)^2 = r2^2` at `w = (-B + s sqrt(B^2 - A C)) / A`,
/// `s` the branch (`plus`); an open piece runs counter-clockwise over
/// `range`, a ring (`None`) over a whole turn.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct MeetCrv {
    /// The carrier's operand.
    pub(super) carrier: usize,
    o: V,
    x: V,
    y: V,
    n: V,
    /// The carrier's local rows for `u` and `v`.
    k: [V; 2],
    r: R,
    g: Vec<V>,
    e: Vec<R>,
    r2: R,
    pub(super) plus: bool,
    pub(super) range: Option<[[Qd; 2]; 2]>,
}

/// The other cylinder of a carrier, as its local rows and offsets.
pub(super) struct Other {
    pub(super) g: Vec<V>,
    pub(super) e: Vec<R>,
    pub(super) r: R,
}

pub(super) fn other_of(f: &Affine, c: &P2, r: &R) -> Other {
    let g = vec![f.row(0).clone(), f.row(1).clone()];
    let e = vec![&dot(&g[0], &f.o) + &c[0], &dot(&g[1], &f.o) + &c[1]];
    Other { g, e, r: r.clone() }
}

/// A sphere as the other quadric (S9d.2): `|p - c|^2 = r^2` over the world's
/// rows.
pub(super) fn other_sphere(c: &V, r: &R) -> Other {
    let rows = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(|e| e.map(int));
    Other {
        g: rows.to_vec(),
        e: c.to_vec(),
        r: r.clone(),
    }
}

impl MeetCrv {
    pub(super) fn new(
        carrier: usize,
        (f, c, r): (&Affine, &P2, &R),
        other: &Other,
        plus: bool,
        range: Option<[[Qd; 2]; 2]>,
    ) -> Self {
        Self {
            carrier,
            o: f.point(&c[0], &c[1], &zero()),
            x: f.x.clone(),
            y: f.y.clone(),
            n: f.n.clone(),
            k: [f.row(0).clone(), f.row(1).clone()],
            r: r.clone(),
            g: other.g.clone(),
            e: other.e.clone(),
            r2: other.r.clone(),
            plus,
            range,
        }
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
        let gn: Vec<R> = self.g.iter().map(|g| dot(g, &self.n)).collect();
        let s: Vec<R> = self
            .g
            .iter()
            .zip(&self.e)
            .map(|(g, e)| dot(g, &base) - e)
            .collect();
        let a = gn.iter().fold(zero(), |acc, x| acc + x * x);
        let b = s.iter().zip(&gn).fold(zero(), |acc, (x, y)| acc + x * y);
        let c = s.iter().fold(zero(), |acc, x| acc + x * x) - &self.r2 * &self.r2;
        let d = &b * &b - &a * &c;
        if d <= zero() || a == zero() {
            return None;
        }
        let k = if self.plus { int(1) } else { int(-1) };
        let w = Qd::new(-&b / &a, k / &a, d);
        Some(qadd(&qv(&base), &qscale(&self.n, &w)))
    }

    /// The carrier's `(cos, sin)` of a point on it.
    pub(super) fn place(&self, x: &QV) -> [Qd; 2] {
        let d = qsub(x, &qv(&self.o));
        let inv = int(1) / &self.r;
        [
            qdot(&d, &self.k[0]).scale(&inv),
            qdot(&d, &self.k[1]).scale(&inv),
        ]
    }

    /// The other cylinder's gradient at a point (halved).
    fn offsets(&self, x: &QV) -> Vec<Qd> {
        self.g
            .iter()
            .zip(&self.e)
            .map(|(g, e)| qdot(x, g).add_r(&-e.clone()))
            .collect()
    }

    fn grad_other(&self, x: &QV) -> QV {
        let s = self.offsets(x);
        self.g
            .iter()
            .zip(&s)
            .fold(qv(&[zero(), zero(), zero()]), |acc, (g, si)| {
                qadd(&acc, &qscale(g, si))
            })
    }

    /// The branch a point of both cylinders lies on: the sign of the other
    /// cylinder's gradient along the carrier's axis.
    pub(super) fn branch_sign(&self, x: &QV) -> Ordering {
        qdot(&self.grad_other(x), &self.n).sign()
    }

    /// Whether a point lies on the piece (on both cylinders, on its branch
    /// and within its run, ends included).
    pub(super) fn on(&self, x: &QV) -> bool {
        let cs = self.place(x);
        let unit = cs[0].mul(&cs[0]).add(&cs[1].mul(&cs[1])).add_r(&int(-1));
        if unit.sign() != Ordering::Equal {
            return false;
        }
        let g = self
            .offsets(x)
            .iter()
            .fold(Qd::rat(zero()), |acc, si| acc.add(&si.mul(si)))
            .add_r(&-(&self.r2 * &self.r2));
        if g.sign() != Ordering::Equal {
            return false;
        }
        let branch = qdot(&self.grad_other(x), &self.n).sign();
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
        let gk = qadd(&qscale(&self.k[0], &cs[0]), &qscale(&self.k[1], &cs[1]));
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
        let d3 = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let sign = if self.plus { 1.0 } else { -1.0 };
        (0..=n)
            .map(|i| {
                let t = t0 + sweep * i as f64 / n as f64;
                let base = [0, 1, 2].map(|j| o[j] + r * (t.cos() * x[j] + t.sin() * y[j]));
                let gn: Vec<f64> = g.iter().map(|gi| d3(gi, &nn)).collect();
                let s: Vec<f64> = g
                    .iter()
                    .zip(&e)
                    .map(|(gi, ei)| d3(gi, &base) - ei)
                    .collect();
                let a: f64 = gn.iter().map(|x| x * x).sum();
                let b: f64 = s.iter().zip(&gn).map(|(x, y)| x * y).sum();
                let c: f64 = s.iter().map(|x| x * x).sum::<f64>() - r2 * r2;
                let w = (-b + sign * (b * b - a * c).max(0.0).sqrt()) / a;
                [0, 1, 2].map(|j| base[j] + w * nn[j])
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
