//! S9d.4a: a whole torus against polyhedral prisms (REVIEW_NOTES.md, S9d.4
//! refined). Its exact model is `(|l|^2 + R^2 - r^2)^2 <= 4 R^2 (l_u^2 +
//! l_v^2)` in the stored frame's coordinates as rationals. On it the
//! distance from the axis is `rho = (|l|^2 + R^2 - r^2) / 2R`, rational in
//! a point's coordinates, so both its angles have exact places: `u` along
//! `(l_u, l_v)`, `v` along `(rho - R, l_w)`. Its wall is traced in four
//! patches of `(u, v)`, cut at a meridian through a rational direction `e`
//! (and its opposite) and at the parallels of a rational angle `v0` and
//! `v0 + pi`, every seam a circle with rational centre and axes.
//!
//! A plane meets the wall in a spiric section: over `u`, with `A = alpha
//! cos u + beta sin u`, the tube's circle meets it where `r A cos v + r mu
//! sin v = -(R A + kappa)` (`trig`: one quadratic surd at a rational `u`),
//! existing where `D_u = r^2 (A^2 + mu^2) - (R A + kappa)^2 >= 0`; over `v`,
//! `alpha cos u + beta sin u = -(mu r sin v + kappa) / (R + r cos v)`, where
//! `D_v = (alpha^2 + beta^2) (R + r cos v)^2 - (mu r sin v + kappa)^2 >= 0`.
//! `D_u` positive all round gives two rings over `u` (loops about the axis),
//! `D_v` positive all round two over `v` (loops about the tube); otherwise
//! each interval of `D_u >= 0` holds a loop, graphs over `v` about its `u`
//! turning points and over `u` about its `v` ones, switched at rational `u`
//! between turning points of different kinds (S9d.2b's rule), each piece
//! verified exactly.
use super::graph::between_ccw;
use super::meet::{tangency, trig, CylPair, EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use super::turned::{middle, near_node, negative_chart, roots, square_sum, Chart, Form};
use crate::identity::Role;
use crate::profile::boolean::Operand;
use crate::solid::split::{q, rational_f64, zero};
use crate::solid::{Construction, Solid};
use crate::topology::{FaceId, Slot};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::Arc;

/// A torus meeting a plane in a way S9d.4a does not take.
fn later(what: &'static str) -> Error {
    Error::OutOfDomain(what)
}

fn loc_of(sides: &[Ordering]) -> Loc {
    if sides.contains(&Ordering::Less) {
        Loc::Out
    } else if sides.contains(&Ordering::Equal) {
        Loc::On
    } else {
        Loc::In
    }
}

/// The cross product's sign of two directions.
fn cross_sign(a: &[Qd; 2], b: &[Qd; 2]) -> Ordering {
    mixed_dot_sign(&[a[0].clone(), a[1].neg()], &[b[1].clone(), b[0].clone()])
}

/// A torus model's own data.
#[derive(Debug, Clone)]
pub(super) struct Ring {
    pub(super) big: R,
    pub(super) small: R,
    /// The meridian seam's direction of `(l_u, l_v)` (unit, rational).
    pub(super) e: P2,
    /// The parallel seam's `(cos v0, sin v0)` (unit, rational).
    pub(super) v0: P2,
}

impl Ring {
    fn k(&self) -> R {
        &self.big * &self.big - &self.small * &self.small
    }

    /// `|l|^2`.
    fn s(l: &QV) -> Qd {
        l[0].mul(&l[0]).add(&l[1].mul(&l[1])).add(&l[2].mul(&l[2]))
    }

    /// The distance from the axis of a point on the torus.
    pub(super) fn rho(&self, l: &QV) -> Qd {
        Self::s(l)
            .add_r(&self.k())
            .scale(&(int(1) / (int(2) * &self.big)))
    }

    /// The model's function, negative inside.
    pub(super) fn value(&self, l: &QV) -> Qd {
        let t = Self::s(l).add_r(&self.k());
        let p2 = l[0].mul(&l[0]).add(&l[1].mul(&l[1]));
        t.mul(&t).sub(&p2.scale(&(int(4) * &self.big * &self.big)))
    }

    /// Its gradient in local coordinates, over four.
    pub(super) fn gradient(&self, l: &QV) -> QV {
        let t = Self::s(l).add_r(&self.k());
        let b2 = int(2) * &self.big * &self.big;
        [
            l[0].mul(&t).sub(&l[0].scale(&b2)),
            l[1].mul(&t).sub(&l[1].scale(&b2)),
            l[2].mul(&t),
        ]
    }

    /// A point's direction of `u`.
    pub(super) fn u_dir(l: &QV) -> [Qd; 2] {
        [l[0].clone(), l[1].clone()]
    }

    /// A point's direction of `v` (on the torus).
    pub(super) fn v_dir(&self, l: &QV) -> [Qd; 2] {
        [self.rho(l).add_r(&-self.big.clone()), l[2].clone()]
    }

    /// Where a point lies in the torus, pushed along directions in turn.
    pub(super) fn member(&self, f: &Affine, p: &QV, dirs: &[QV]) -> Loc {
        let l = f.local_q(p);
        let mut s = self.value(&l).sign();
        if s == Ordering::Equal {
            let g = self.gradient(&l);
            for d in dirs {
                let ld = f.local_dir_q(d);
                s = g[0]
                    .mul(&ld[0])
                    .add(&g[1].mul(&ld[1]))
                    .add(&g[2].mul(&ld[2]))
                    .sign();
                if s != Ordering::Equal {
                    break;
                }
            }
        }
        loc_of(&[s.reverse()])
    }

    /// The side of the meridian seam (`plus`: counter-clockwise from `e`)
    /// and of the parallel seam (`upper`: from `v0` towards `v0 + pi`).
    fn sides(&self, l: &QV) -> (Ordering, Ordering) {
        let e = [Qd::rat(self.e[0].clone()), Qd::rat(self.e[1].clone())];
        let v0 = [Qd::rat(self.v0[0].clone()), Qd::rat(self.v0[1].clone())];
        (
            cross_sign(&e, &Self::u_dir(l)),
            cross_sign(&v0, &self.v_dir(l)),
        )
    }

    /// Where a point on the torus lies in a patch.
    pub(super) fn in_face(&self, f: &Affine, kind: FaceKind, p: &QV) -> Loc {
        let FaceKind::Patch(upper, plus) = kind else {
            unreachable!("a torus's faces are patches")
        };
        let (su, sv) = self.sides(&f.local_q(p));
        let want = |b: bool| if b { Ordering::Greater } else { Ordering::Less };
        let side = |s: Ordering, b: bool| {
            if s == Ordering::Equal {
                Ordering::Equal
            } else if s == want(b) {
                Ordering::Greater
            } else {
                Ordering::Less
            }
        };
        loc_of(&[side(su, plus), side(sv, upper)])
    }

    /// Whether a meeting may be at a seam (on a meridian or parallel seam).
    pub(super) fn on_seam(&self, f: &Affine, p: &QV) -> bool {
        let (su, sv) = self.sides(&f.local_q(p));
        su == Ordering::Equal || sv == Ordering::Equal
    }

    /// A point's binary64 parameters on a patch: its angles from the seams
    /// (the minus and lower patches' shifted by a turn below), `u` scaled
    /// by `R`. The model's numbers are rounded once, for every point.
    pub(super) fn params(&self, f: &Affine) -> impl Fn([f64; 3], FaceKind) -> [f64; 2] {
        let fl = |v: &V| v.clone().map(|y| rational_f64(&y));
        let (o, x, y, n) = (fl(&f.o), fl(&f.x), fl(&f.y), fl(&f.n));
        let (big, small) = (rational_f64(&self.big), rational_f64(&self.small));
        let fl2 = |v: &P2| [rational_f64(&v[0]), rational_f64(&v[1])];
        let (e, v0) = (fl2(&self.e), fl2(&self.v0));
        move |p, kind| {
            let FaceKind::Patch(upper, plus) = kind else {
                unreachable!("a patch")
            };
            let d = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
            let l = super::graph::solve3(&x, &y, &n, &d);
            let rho = l[0].hypot(l[1]);
            let angle = |a: f64, b: f64, base: &[f64; 2], first: bool| {
                let (c0, s0) = (base[0], base[1]);
                let rel = (b * c0 - a * s0).atan2(a * c0 + b * s0);
                // Within the half's turn: a point on a seam may round across.
                if first {
                    if rel < -std::f64::consts::FRAC_PI_2 {
                        rel + std::f64::consts::TAU
                    } else {
                        rel
                    }
                } else if rel > std::f64::consts::FRAC_PI_2 {
                    rel - std::f64::consts::TAU
                } else {
                    rel
                }
            };
            let u = angle(l[0], l[1], &e, plus);
            let v = angle(rho - big, l[2], &v0, upper);
            [u * big, v * small]
        }
    }
}

/// A plane's section of a torus as a graph over `u` (or over `v`): the
/// plane `alpha l_u + beta l_v + mu l_w + kappa = 0` in the torus's
/// coordinates, the branch `plus` (`trig`'s second root), over the open
/// counter-clockwise `range` of its parameter's direction (a ring: none).
#[derive(Debug, Clone, PartialEq)]
pub(super) struct TorusSec {
    /// The operand whose torus it lies on.
    pub(super) carrier: usize,
    f: Affine,
    big: R,
    small: R,
    pub(super) plane: [R; 4],
    pub(super) over_v: bool,
    pub(super) plus: bool,
    pub(super) range: Option<[[Qd; 2]; 2]>,
}

impl TorusSec {
    fn ring(&self) -> Ring {
        Ring {
            big: self.big.clone(),
            small: self.small.clone(),
            e: [int(1), zero()],
            v0: [int(1), zero()],
        }
    }

    /// The point at a rational `(cos, sin)` of its parameter, on its
    /// branch (none where the section has no point there).
    pub(super) fn at(&self, cs: &[R; 2]) -> Option<QV> {
        let [a, b, m, k] = &self.plane;
        let (big, small) = (&self.big, &self.small);
        let sols = if self.over_v {
            let rv = big + small * &cs[0];
            if rv == zero() {
                return None;
            }
            let g = -(m * small * &cs[1] + k) / &rv;
            trig(a, b, &g).ok()??
        } else {
            let aa = a * &cs[0] + b * &cs[1];
            trig(&(small * &aa), &(small * m), &-(big * &aa + k)).ok()??
        };
        let other = sols.get(usize::from(self.plus))?.clone();
        let (u, v) = if self.over_v {
            (other, [Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())])
        } else {
            ([Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())], other)
        };
        let rho = v[0].scale(small).add_r(big);
        let l = [rho.mul(&u[0]), rho.mul(&u[1]), v[1].scale(small)];
        Some(qadd(
            &qv(&self.f.o),
            &qadd(
                &qadd(&qscale(&self.f.x, &l[0]), &qscale(&self.f.y, &l[1])),
                &qscale(&self.f.n, &l[2]),
            ),
        ))
    }

    /// A point's place: the `(cos, sin)` of the graph's parameter (the
    /// directions over their lengths, `r` and `rho`, exact).
    pub(super) fn place(&self, x: &QV) -> [Qd; 2] {
        let l = self.f.local_q(x);
        let ring = self.ring();
        if self.over_v {
            let inv = int(1) / &self.small;
            ring.v_dir(&l).map(|c| c.scale(&inv))
        } else {
            let inv = ring.rho(&l).recip().expect("a point off the axis");
            Ring::u_dir(&l).map(|c| c.mul(&inv))
        }
    }

    /// The branch a point of the section lies on: the sign of the other
    /// angle's turn from the solution's middle direction.
    fn branch(&self, l: &QV) -> Ordering {
        let [a, b, m, _] = &self.plane;
        let ring = self.ring();
        if self.over_v {
            // cross((alpha, beta), u) = alpha l_v - beta l_u.
            l[1].scale(a).sub(&l[0].scale(b)).sign()
        } else {
            // cross((r A, r mu), v) with A = (alpha l_u + beta l_v) / rho:
            // sign of (alpha l_u + beta l_v) l_w - mu rho (rho - R).
            let rho = ring.rho(l);
            l[0].scale(a)
                .add(&l[1].scale(b))
                .mul(&l[2])
                .sub(&rho.mul(&rho.add_r(&-self.big.clone())).scale(m))
                .sign()
        }
    }

    /// Whether a point lies on the piece (the torus, the plane, its branch
    /// and its range, ends included).
    pub(super) fn on(&self, x: &QV) -> bool {
        let l = self.f.local_q(x);
        let [a, b, m, k] = &self.plane;
        let in_plane = l[0]
            .scale(a)
            .add(&l[1].scale(b))
            .add(&l[2].scale(m))
            .add_r(k)
            .sign()
            == Ordering::Equal;
        if !in_plane || self.ring().value(&l).sign() != Ordering::Equal {
            return false;
        }
        let want = if self.plus {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        if self.branch(&l) != want {
            return false;
        }
        match &self.range {
            None => true,
            Some([lo, hi]) => {
                let p = self.place(x);
                super::graph::same_dir(&p, lo)
                    || super::graph::same_dir(&p, hi)
                    || between_ccw(lo, &p, hi)
            }
        }
    }

    /// The unit-free tangent at a point, running with the parameter: the
    /// torus's and the plane's normals crossed.
    pub(super) fn tangent(&self, x: &QV) -> QV {
        let l = self.f.local_q(x);
        let ring = self.ring();
        let g = ring.gradient(&l);
        let row = |k: usize| self.f.row(k).clone();
        let gw = qadd(
            &qadd(&qscale(&row(0), &g[0]), &qscale(&row(1), &g[1])),
            &qscale(&row(2), &g[2]),
        );
        let [a, b, m, _] = &self.plane;
        let mw = add(
            &add(&scale(&row(0), a), &scale(&row(1), b)),
            &scale(&row(2), m),
        );
        let t = qcross(&gw, &qv(&mw));
        let lt = self.f.local_dir_q(&t);
        // The parameter's turn along t.
        let run = if self.over_v {
            // d(v_dir): (d rho, d w), d rho = (l_u t_u + l_v t_v) / rho.
            let rho = ring.rho(&l);
            let drho = l[0].mul(&lt[0]).add(&l[1].mul(&lt[1]));
            rho.add_r(&-self.big.clone())
                .mul(&rho)
                .mul(&lt[2])
                .sub(&l[2].mul(&drho))
                .sign()
        } else {
            l[0].mul(&lt[1]).sub(&l[1].mul(&lt[0])).sign()
        };
        if run == Ordering::Less {
            t.map(|c| c.neg())
        } else {
            t
        }
    }

    /// Binary64 points from the parameter's angle `t0` turning `sweep`.
    pub(super) fn samples(&self, t0: f64, sweep: f64, n: usize) -> Vec<[f64; 3]> {
        let fl = |v: &V| v.clone().map(|y| rational_f64(&y));
        let (o, x, y, nn) = (fl(&self.f.o), fl(&self.f.x), fl(&self.f.y), fl(&self.f.n));
        let [a, b, m, k] = self.plane.clone().map(|z| rational_f64(&z));
        let (big, small) = (rational_f64(&self.big), rational_f64(&self.small));
        let sign = if self.plus { 1.0 } else { -1.0 };
        (0..=n)
            .map(|i| {
                let t = t0 + sweep * i as f64 / n as f64;
                let (u, v) = if self.over_v {
                    let q = -(m * small * t.sin() + k) / (big + small * t.cos());
                    let u = b.atan2(a) + sign * (q / a.hypot(b)).clamp(-1.0, 1.0).acos();
                    (u, t)
                } else {
                    let aa = a * t.cos() + b * t.sin();
                    let (px, py) = (small * aa, small * m);
                    let ratio = -(big * aa + k) / px.hypot(py);
                    (t, py.atan2(px) + sign * ratio.clamp(-1.0, 1.0).acos())
                };
                let rho = big + small * v.cos();
                let l = [rho * u.cos(), rho * u.sin(), small * v.sin()];
                [0, 1, 2].map(|j| o[j] + l[0] * x[j] + l[1] * y[j] + l[2] * nn[j])
            })
            .collect()
    }
}

/// Where a line meets the torus: the roots of its quartic, algebraic
/// (`Q(alpha)`); a repeated root is a tangency.
pub(super) fn line_torus(p: &QV, d: &V, f: &Affine, ring: &Ring) -> Result<EdgeMeet> {
    let lp = f.local_q(p);
    let (Some(l0), Some(l1), Some(l2)) = (lp[0].rational(), lp[1].rational(), lp[2].rational())
    else {
        return Err(Error::ComputationLimit(
            "an irrational line against a torus",
        ));
    };
    let l = [l0.clone(), l1.clone(), l2.clone()];
    let ld = f.local_dir(d);
    // |l + t d|^2 + K = s2 t^2 + s1 t + s0; rho^2 part likewise.
    let s2 = dot(&ld, &ld);
    let s1 = int(2) * dot(&l, &ld);
    let s0 = dot(&l, &l) + ring.k();
    let p2 = &ld[0] * &ld[0] + &ld[1] * &ld[1];
    let p1 = int(2) * (&l[0] * &ld[0] + &l[1] * &ld[1]);
    let p0 = &l[0] * &l[0] + &l[1] * &l[1];
    let b4 = int(4) * &ring.big * &ring.big;
    let poly = vec![
        &s0 * &s0 - &b4 * &p0,
        int(2) * &s0 * &s1 - &b4 * &p1,
        &s1 * &s1 + int(2) * &s0 * &s2 - &b4 * &p2,
        int(2) * &s1 * &s2,
        &s2 * &s2,
    ];
    let poly = super::turned::trim(poly);
    let mut out = Vec::new();
    for root in roots(&poly)? {
        let g = Arc::new(Gen::new(poly.clone(), root));
        let t = Qd::of(K::generator(&g));
        let x = qadd(p, &qscale(d, &t));
        out.push((Pos::T(t), x));
    }
    Ok(EdgeMeet::Points(out))
}

/// The plane's coefficients in the torus's coordinates.
fn plane_in(f: &Affine, p0: &V, m: &V) -> [R; 4] {
    [
        dot(m, &f.x),
        dot(m, &f.y),
        dot(m, &f.n),
        dot(m, &sub(&f.o, p0)),
    ]
}

/// `D_u` and `D_v` as quadratic forms in `(cos, sin)` of `u` and of `v`.
fn discriminants(ring: &Ring, plane: &[R; 4]) -> (Form, Form) {
    let [a, b, m, k] = plane;
    let (big, small) = (&ring.big, &ring.small);
    // D_u = r^2 (A^2 + mu^2) - (R A + kappa)^2, A = [0, alpha, beta].
    let du = square_sum(&[[zero(), small * a, small * b], [small * m, zero(), zero()]])
        .sub(&square_sum(&[[k.clone(), big * a, big * b]]));
    // D_v = (alpha^2 + beta^2) (R + r cos v)^2 - (mu r sin v + kappa)^2.
    let ab = a * a + b * b;
    let dv = square_sum(&[[big.clone(), small.clone(), zero()]])
        .scaled(&ab)
        .sub(&square_sum(&[[k.clone(), zero(), m * small]]));
    (du, dv)
}

/// Whether a form vanishes on the whole circle.
fn vanishes(a: &Form) -> bool {
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    super::turned::trim(a.poly(&chart)).is_empty() && a.value(&[int(-1), zero()]) == zero()
}

/// A plane's section of a torus (the torus operand `k`): rings over `u` or
/// over `v`, loops of both kinds of graph with their switches, or apart.
pub(super) fn plane_torus(
    k: usize,
    f: &Affine,
    ring: &Ring,
    p0: &V,
    m: &V,
    res: f64,
) -> Result<CylPair> {
    let plane = plane_in(f, p0, m);
    let (du, dv) = discriminants(ring, &plane);
    if vanishes(&du) || vanishes(&dv) {
        return Err(tangency());
    }
    let piece = |over_v: bool, plus: bool, range: Option<[[Qd; 2]; 2]>| TorusSec {
        carrier: k,
        f: f.clone(),
        big: ring.big.clone(),
        small: ring.small.clone(),
        plane: plane.clone(),
        over_v,
        plus,
        range,
    };
    let unit = Chart {
        c0: int(1),
        s0: zero(),
    };
    let bound = int(1);
    let mixed = |pieces: Vec<TorusSec>, switches: Vec<QV>| {
        Ok(CylPair::Mixed(Box::new(super::spheres::Mixed {
            pieces: pieces
                .into_iter()
                .map(|p| Crv::Torus(Box::new(p)))
                .collect(),
            switches,
        })))
    };
    // Rings about the axis (over u) or about the tube (over v), decided
    // before any root is isolated (the other graph's discriminant may be a
    // square's negative, its roots double: a plane normal to the axis, or
    // through it).
    let chart_u = negative_chart(&du)?;
    if chart_u.is_none() {
        near_node(&bound, &du, &unit, res)?;
        return mixed(
            vec![piece(false, true, None), piece(false, false, None)],
            Vec::new(),
        );
    }
    let chart_v = negative_chart(&dv)?;
    if chart_v.is_none() {
        near_node(&bound, &dv, &unit, res)?;
        return mixed(
            vec![piece(true, true, None), piece(true, false, None)],
            Vec::new(),
        );
    }
    let chart_u = chart_u.expect("a chart");
    near_node(&bound, &du, &chart_u, res)?;
    near_node(&bound, &dv, chart_v.as_ref().expect("a chart"), res)?;
    let pu = du.poly(&chart_u);
    let mut ru = roots(&pu)?;
    if ru.is_empty() {
        return Ok(CylPair::Apart);
    }
    if ru.len() % 2 != 0 {
        return Err(Error::ComputationLimit(
            "an odd count of a torus section's turning points",
        ));
    }
    for r in ru.iter_mut() {
        r.refine_for_signs(160);
    }
    // The v graph's turning points: where D_v vanishes, the u there
    // unique (`phi` or `phi + pi` by the sign of the right side).
    let chart_v = chart_v.expect("a chart");
    let pv = dv.poly(&chart_v);
    let mut rv = roots(&pv)?;
    for r in rv.iter_mut() {
        r.refine_for_signs(160);
    }
    let [a, b, mu, kappa] = plane.clone().map(|x| rational_f64(&x));
    let (big, small) = (rational_f64(&ring.big), rational_f64(&ring.small));
    let (cu0, su0) = (rational_f64(&chart_u.c0), rational_f64(&chart_u.s0));
    let (cv0, sv0) = (rational_f64(&chart_v.c0), rational_f64(&chart_v.s0));
    let chart_angle = |c0: f64, s0: f64, t: f64| s0.atan2(c0) + 2.0 * t.atan();
    let to_chart = |c0: f64, s0: f64, ang: f64| {
        let (c, s) = (ang.cos(), ang.sin());
        let (cr, sr) = (c0 * c + s0 * s, c0 * s - s0 * c);
        sr / (1.0 + cr)
    };
    // (chart t of u, the u graph's branch there).
    let mut v_turns: Vec<(f64, bool)> = Vec::new();
    for root in &rv {
        let v = chart_angle(cv0, sv0, rational_f64(&middle(root)));
        let g = -(mu * small * v.sin() + kappa) / (big + small * v.cos());
        let u = b.atan2(a) + if g > 0.0 { 0.0 } else { std::f64::consts::PI };
        // The u graph's branch: sign of (alpha l_u + beta l_v) l_w - mu rho (rho - R).
        let rho = big + small * v.cos();
        let (lu, lv, lw) = (rho * u.cos(), rho * u.sin(), small * v.sin());
        let br = (a * lu + b * lv) * lw - mu * rho * (rho - big);
        if br.abs() < 1e-12 * (1.0 + big * big) {
            return Err(Error::ComputationLimit("a turning point of both graphs"));
        }
        v_turns.push((to_chart(cu0, su0, u), br > 0.0));
    }
    let u_graph = |plus: bool| piece(false, plus, None);
    let v_probe = piece(true, true, None);
    let mut pieces = Vec::new();
    let mut switches = Vec::new();
    for pair in ru.chunks(2) {
        let (ra, rb) = (&pair[0], &pair[1]);
        let (a64, b64) = (rational_f64(&middle(ra)), rational_f64(&middle(rb)));
        // Events along the loop: lambda in [0, 2), the + branch over u from
        // a to b, the - branch back; true for u's turning points.
        let mut events: Vec<(f64, bool)> = vec![(0.0, true), (1.0, true)];
        for &(t, plus) in &v_turns {
            if t <= a64 || t >= b64 {
                continue;
            }
            let l = (t - a64) / (b64 - a64);
            events.push((if plus { l } else { 2.0 - l }, false));
        }
        events.sort_by(|x, y| x.0.total_cmp(&y.0));
        let n = events.len();
        if n < 3 {
            return Err(Error::ComputationLimit(
                "a loop without turning points of the other graph",
            ));
        }
        for i in 0..n {
            let l1 = events[(i + 1) % n].0 + if i + 1 == n { 2.0 } else { 0.0 };
            if l1 - events[i].0 < 1e-9 {
                return Err(Error::Degenerate("two turning points within rounding"));
            }
        }
        let at = |l: f64| -> (R, bool) {
            let l = l.rem_euclid(2.0);
            if l < 1.0 {
                (q(a64 + l * (b64 - a64)), true)
            } else {
                (q(b64 - (l - 1.0) * (b64 - a64)), false)
            }
        };
        let point = |l: f64| -> Result<(QV, R, bool)> {
            let (t, plus) = at(l);
            if ra.compare_rational(&t) != Ordering::Less
                || rb.compare_rational(&t) != Ordering::Greater
            {
                return Err(Error::ComputationLimit("a switch point off its loop"));
            }
            let p = u_graph(plus)
                .at(&chart_u.at(&t))
                .ok_or(Error::ComputationLimit("a switch point off its piece"))?;
            Ok((p, t, plus))
        };
        let mut sw: Vec<(f64, QV, R, bool)> = Vec::new();
        for i in 0..n {
            let (e0, e1) = (events[i], events[(i + 1) % n]);
            if e0.1 == e1.1 {
                continue;
            }
            let l1 = e1.0 + if i + 1 == n { 2.0 } else { 0.0 };
            let l = 0.5 * (e0.0 + l1);
            let (p, t, plus) = point(l)?;
            sw.push((l.rem_euclid(2.0), p, t, plus));
        }
        sw.sort_by(|x, y| x.0.total_cmp(&y.0));
        let m = sw.len();
        if m < 2 {
            return Err(Error::ComputationLimit("a loop with one switch"));
        }
        for i in 0..m {
            let (s0, s1) = (&sw[i], &sw[(i + 1) % m]);
            let l1 = s1.0 + if i + 1 == m { 2.0 } else { 0.0 };
            let first = events
                .iter()
                .find(|e| {
                    let l = if e.0 < s0.0 { e.0 + 2.0 } else { e.0 };
                    l > s0.0 && l < l1
                })
                .ok_or(Error::ComputationLimit("a run without turning points"))?;
            if first.1 {
                // u's turning points: a graph over v.
                let lm = 0.5
                    * (s0.0
                        + if first.0 < s0.0 {
                            first.0 + 2.0
                        } else {
                            first.0
                        });
                let (mid, _, _) = point(lm)?;
                let lv = |x: &QV| f.local_q(x);
                let sign = v_probe.branch(&lv(&mid));
                if sign == Ordering::Equal
                    || v_probe.branch(&lv(&s0.1)) != sign
                    || v_probe.branch(&lv(&s1.1)) != sign
                {
                    return Err(Error::ComputationLimit("a run over v changing branch"));
                }
                let (p0, p1, pm) = (
                    v_probe.place(&s0.1),
                    v_probe.place(&s1.1),
                    v_probe.place(&mid),
                );
                let range = if between_ccw(&p0, &pm, &p1) {
                    [p0, p1]
                } else {
                    [p1, p0]
                };
                verify(&pv, &chart_v, &range)?;
                pieces.push(piece(true, sign == Ordering::Greater, Some(range)));
            } else {
                // v's turning points: a graph over u, one branch.
                if s0.3 != s1.3 {
                    return Err(Error::ComputationLimit("a run over u changing branch"));
                }
                let (lo, hi) = if s0.2 < s1.2 {
                    (&s0.2, &s1.2)
                } else {
                    (&s1.2, &s0.2)
                };
                if ru.iter().any(|r| {
                    r.compare_rational(lo) == Ordering::Greater
                        && r.compare_rational(hi) == Ordering::Less
                }) {
                    return Err(Error::ComputationLimit("a turning point inside a piece"));
                }
                let rng = [lo, hi].map(|t| {
                    let cs = chart_u.at(t);
                    [Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]
                });
                pieces.push(piece(false, s0.3, Some(rng)));
            }
        }
        switches.extend(sw.into_iter().map(|s| s.1));
    }
    mixed(pieces, switches)
}

/// No root of a graph's discriminant within a counter-clockwise range
/// (clear of the chart's antipode): exact Sturm counts at its ends.
fn verify(p: &super::turned::Poly, chart: &Chart, range: &[[Qd; 2]; 2]) -> Result<()> {
    let (Some(t0), Some(t1)) = (chart.t_of(&range[0]), chart.t_of(&range[1])) else {
        return Err(Error::ComputationLimit(
            "a piece through a chart's antipode",
        ));
    };
    if t0.cmp(&t1) != Ordering::Less {
        return Err(Error::ComputationLimit(
            "a piece through a chart's antipode",
        ));
    }
    let chain = super::turned::sturm(&super::turned::trim(p.clone()));
    if super::turned::changes(&chain, &t0) != super::turned::changes(&chain, &t1) {
        return Err(Error::ComputationLimit("a turning point inside a piece"));
    }
    Ok(())
}

/// The exact model of a whole torus; `seam` picks the rational directions
/// of its meridian and parallel seams.
pub(super) fn model(solid: &Solid, op: Operand, seam: &R) -> Result<Prism> {
    let Construction::Torus {
        major,
        minor,
        low,
        high,
        angle,
        ..
    } = &solid.construction
    else {
        unreachable!("a torus");
    };
    if high - low != std::f64::consts::TAU || *angle != std::f64::consts::TAU {
        return Err(later("a Boolean of a torus segment or wedge (S9d.4b)"));
    }
    let f = Affine::new(&solid.frame)?;
    let t = &solid.topology;
    let id = |slot: Slot| {
        t.id_of(slot)
            .ok_or(Error::InvalidTopology("an unnamed slot"))
    };
    let (big, small) = (q(*major), q(*minor));
    let e = circle_point(&[zero(), zero()], &int(1), seam);
    let v0 = circle_point(&[zero(), zero()], &int(1), &(seam + &(int(1) / int(7))));
    let ring = Ring {
        big: big.clone(),
        small: small.clone(),
        e: e.clone(),
        v0: v0.clone(),
    };
    // The wall's four patches: (upper, plus), (upper, minus), (lower,
    // plus), (lower, minus).
    let wall = FaceId(0);
    let kinds = [(true, true), (true, false), (false, true), (false, false)];
    let face_of = |upper: bool, plus: bool| {
        kinds
            .iter()
            .position(|&k| k == (upper, plus))
            .expect("a patch")
    };
    let mut faces = Vec::new();
    for &(upper, plus) in &kinds {
        faces.push(MFace {
            kind: FaceKind::Patch(upper, plus),
            surf: Surf::Torus,
            id: id(Slot::Face(wall))?,
            stored: t.faces()[wall.0].surface.clone(),
            sense: t.faces()[wall.0].sense,
        });
    }
    // The seams' frame: x_e along e, y_e a quarter turn on, n the axis.
    let xe = f.vector(&e[0], &e[1], &zero());
    let ye = f.vector(&-e[1].clone(), &e[0], &zero());
    let n = f.n.clone();
    let point = |su: i64, sv: i64| {
        let (cv, sv_) = (&v0[0] * int(sv), &v0[1] * int(sv));
        let rho = &big + &small * &cv;
        add(
            &f.o,
            &add(&scale(&xe, &(&rho * int(su))), &scale(&n, &(&small * &sv_))),
        )
    };
    let mut verts = Vec::new();
    // Vertices: (u side, v side) = (+e, +v0), (-e, +v0), (+e, -v0), (-e, -v0).
    for (su, sv) in [(1, 1), (-1, 1), (1, -1), (-1, -1)] {
        verts.push(MVert {
            p: qv(&point(su, sv)),
            id: None,
        });
    }
    let vx = |su: i64, sv: i64| match (su, sv) {
        (1, 1) => 0,
        (-1, 1) => 1,
        (1, -1) => 2,
        _ => 3,
    };
    let east = [Qd::rat(int(1)), Qd::rat(zero())];
    let west = [Qd::rat(int(-1)), Qd::rat(zero())];
    let mut edges = Vec::new();
    let mut seam_edge = 0;
    let mut push = |edges: &mut Vec<MEdge>,
                    curve: Crv,
                    from: [Qd; 2],
                    to: [Qd; 2],
                    ends: [usize; 2],
                    faces: [usize; 2]| {
        edges.push(MEdge {
            kind: EdgeKind::Split(seam_edge),
            curve,
            arc: Some((from, to, true)),
            start: ends[0],
            end: ends[1],
            faces,
            id: None,
        });
        seam_edge += 1;
    };
    // The parallels at v0 (upper on the left) and v0 + pi (lower on the
    // left), counter-clockwise about the axis in the (x_e, y_e) basis.
    for (sv, left_upper) in [(1i64, true), (-1, false)] {
        let rho = &big + &small * &v0[0] * int(sv);
        let c = add(&f.o, &scale(&n, &(&small * &v0[1] * int(sv))));
        let curve = || Crv::Conic {
            c: c.clone(),
            a: scale(&xe, &rho),
            b: scale(&ye, &rho),
        };
        for (plus, from, to, ends) in [
            (true, east.clone(), west.clone(), [vx(1, sv), vx(-1, sv)]),
            (false, west.clone(), east.clone(), [vx(-1, sv), vx(1, sv)]),
        ] {
            let faces = [face_of(left_upper, plus), face_of(!left_upper, plus)];
            push(&mut edges, curve(), from, to, ends, faces);
        }
    }
    // The meridians at e (the minus half on the left) and at -e (the plus
    // half), running with v in the (x_e or -x_e, n) basis.
    let (vu, vl) = (
        [Qd::rat(v0[0].clone()), Qd::rat(v0[1].clone())],
        [Qd::rat(-v0[0].clone()), Qd::rat(-v0[1].clone())],
    );
    for (su, left_plus) in [(1i64, false), (-1, true)] {
        let c = add(&f.o, &scale(&xe, &(&big * int(su))));
        let curve = || Crv::Conic {
            c: c.clone(),
            a: scale(&xe, &(&small * int(su))),
            b: scale(&n, &small),
        };
        for (upper, from, to, ends) in [
            (true, vu.clone(), vl.clone(), [vx(su, 1), vx(su, -1)]),
            (false, vl.clone(), vu.clone(), [vx(su, -1), vx(su, 1)]),
        ] {
            let faces = [face_of(upper, left_plus), face_of(upper, !left_plus)];
            push(&mut edges, curve(), from, to, ends, faces);
        }
    }
    let mut info = BTreeMap::new();
    for (eid, _) in t.ids() {
        let role = t.derivation(eid).map_or(Role::External, |d| d.role);
        info.insert(eid, (op, role));
    }
    // One box for every patch: the torus's in its frame.
    let reach = &big + &small;
    let fl = |p: &V| p.clone().map(|x| rational_f64(&x));
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for (du, dv, dw) in [
        (1, 1, 1),
        (1, 1, -1),
        (1, -1, 1),
        (1, -1, -1),
        (-1, 1, 1),
        (-1, 1, -1),
        (-1, -1, 1),
        (-1, -1, -1),
    ] {
        let p = fl(&f.point(
            &(&reach * int(du)),
            &(&reach * int(dv)),
            &(&small * int(dw)),
        ));
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    for k in 0..3 {
        let m = 1e-9 * (1.0 + lo[k].abs().max(hi[k].abs()));
        lo[k] -= m;
        hi[k] += m;
    }
    Ok(Prism {
        frame: solid.frame,
        tolerance: solid.resolution(),
        region: id(Slot::Region(crate::topology::RegionId(1)))?,
        f,
        lo: zero(),
        hi: zero(),
        bounds: Vec::new(),
        faces,
        edges,
        verts,
        info,
        boxes: vec![(lo, hi); 4],
        ball: None,
        funnel: None,
        ring: Some(ring),
    })
}
