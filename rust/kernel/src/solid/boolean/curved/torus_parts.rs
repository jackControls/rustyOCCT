//! S9d.4c: a sphere's cap or zone against a torus, and torus v-segments and
//! wedges against curved faces (REVIEW_NOTES.md, S9d.4c refined).
//!
//! A sphere's own circle `c + dx x + dy y`, `xx dx^2 + 2 xy dx dy + yy dy^2 =
//! r2` (a cap's or zone's rim: a surd radius, and unequal axes on a turned
//! frame's stored axis) meets a torus where the torus's quartic on the
//! circle's plane, a quartic in `(dx, dy)`, vanishes on the circle: in
//! coordinates `(s, t)` turned by a rational rotation the circle is a
//! quadratic in `t` of constant leading coefficient, the quartic reduced
//! modulo it is `r1(s) t + r0(s)`, and their resultant, of degree at most
//! eight in `s`, has one simple real root per crossing, `t = -r0 / r1` there
//! (S9d.2c's rotations, retried where a root repeats or `r1` vanishes). A
//! crossing within the resolution of a tangency is `Degenerate`: the counts
//! with the torus's function offset by its gradient's bound times the
//! resolution must be its own.
//!
//! A part's rim (a segment's end plane's ring over `u`, a wedge's end
//! half-plane's ring over `v`) fixes one of the torus's angles at a
//! quadratic surd: the other surface's function on the angles (S9d.4b.2's
//! `G`) is there a form in the free angle `A + sqrt(d) B`, whose chart
//! polynomial's roots are among its norm's, `A^2 - d B^2`; the rim's
//! crossings are the norm's real roots at which `A + sqrt(d) B` vanishes,
//! decided exactly in `Q(alpha)(sqrt(d))` (the conjugate's belong to the
//! other ring, or the opposite half-plane's circle), a double root a
//! tangency.
use super::meet::{tangency, EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use super::sphere::Circ;
use super::torus::{Ring, TorusSec};
use super::torus_curved::{meeting_form, Far};
use super::turned::{pderiv, roots, roots_repeated, trim, Chart, Form, Poly};
use crate::polynomial::real::AlgebraicRoot;
use crate::solid::split::{q, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::Arc;

fn limit(what: &'static str) -> Error {
    Error::ComputationLimit(what)
}

// ------------------------------------------------------------ bivariate

/// A polynomial in `t` whose coefficients are polynomials in `s`
/// (ascending in both).
type Bv = Vec<Poly>;

fn bv_trim(mut p: Bv) -> Bv {
    for c in p.iter_mut() {
        *c = trim(std::mem::take(c));
    }
    while p.last().is_some_and(|c| c.is_empty()) {
        p.pop();
    }
    p
}

fn padd(a: &Poly, b: &Poly) -> Poly {
    super::turned::padd(a, b)
}

fn pmul(a: &Poly, b: &Poly) -> Poly {
    super::turned::pmul(a, b)
}

fn pscale(a: &Poly, k: &R) -> Poly {
    super::turned::pscale(a, k)
}

fn bv_add(a: &Bv, b: &Bv) -> Bv {
    let n = a.len().max(b.len());
    bv_trim(
        (0..n)
            .map(|j| {
                padd(
                    a.get(j).unwrap_or(&Vec::new()),
                    b.get(j).unwrap_or(&Vec::new()),
                )
            })
            .collect(),
    )
}

fn bv_mul(a: &Bv, b: &Bv) -> Bv {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out: Bv = vec![Vec::new(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = padd(&out[i + j], &pmul(x, y));
        }
    }
    bv_trim(out)
}

fn bv_scale(a: &Bv, k: &R) -> Bv {
    bv_trim(a.iter().map(|c| pscale(c, k)).collect())
}

/// `a + b s + c t`.
fn bv_lin(a: &R, b: &R, c: &R) -> Bv {
    bv_trim(vec![vec![a.clone(), b.clone()], vec![c.clone()]])
}

fn peval_k(p: &Poly, x: &K) -> K {
    p.iter()
        .rev()
        .fold(K::Rat(zero()), |acc, c| acc.mul(x).add(&K::Rat(c.clone())))
}

/// The rotations tried in turn (rational points of the unit circle).
fn rotations() -> Vec<[R; 2]> {
    let r = |a: i64, b: i64, c: i64| [R::new(a.into(), c.into()), R::new(b.into(), c.into())];
    vec![
        r(1, 0, 1),
        r(3, 4, 5),
        r(5, 12, 13),
        r(8, 15, 17),
        r(7, 24, 25),
        r(20, 21, 29),
        r(12, 35, 37),
        r(9, 40, 41),
    ]
}

/// The count of distinct real roots, `None` where one is repeated.
fn count(p: &Poly) -> Result<Option<usize>> {
    if p.len() <= 1 {
        return Ok(Some(0));
    }
    match roots(p) {
        Ok(rs) => Ok(Some(rs.len())),
        Err(Error::Degenerate(_)) => Ok(None),
        Err(e) => Err(e),
    }
}

// ------------------------------------------------------------ a sphere's circle

/// The torus's function in a rotated circle's coordinates `(s, t)` (its
/// local point `l0 + dx lx + dy ly`, `dx = cr s - sr t`, `dy = sr s + cr
/// t`), offset by `k`.
fn torus_in_st(l: &[[R; 3]; 3], big: &R, small: &R, cr: &R, sr: &R, k: &R) -> Bv {
    let [l0, lx, ly] = l;
    let lin: Vec<Bv> = (0..3)
        .map(|i| {
            bv_lin(
                &l0[i],
                &(cr * &lx[i] + sr * &ly[i]),
                &(cr * &ly[i] - sr * &lx[i]),
            )
        })
        .collect();
    let p2 = bv_add(&bv_mul(&lin[0], &lin[0]), &bv_mul(&lin[1], &lin[1]));
    let mut s = bv_add(&p2, &bv_mul(&lin[2], &lin[2]));
    s = bv_add(&s, &vec![vec![big * big - small * small]]);
    let mut t = bv_add(&bv_mul(&s, &s), &bv_scale(&p2, &(int(-4) * big * big)));
    t = bv_add(&t, &vec![vec![k.clone()]]);
    t
}

/// The circle `xx dx^2 + 2 xy dx dy + yy dy^2 - r2` in `(s, t)`.
fn circle_in_st(circ: &Circ, cr: &R, sr: &R) -> Bv {
    let (xx, yy, xy) = (
        dot(&circ.x, &circ.x),
        dot(&circ.y, &circ.y),
        dot(&circ.x, &circ.y),
    );
    let dx = bv_lin(&zero(), cr, &-sr.clone());
    let dy = bv_lin(&zero(), sr, cr);
    let mut e = bv_add(
        &bv_scale(&bv_mul(&dx, &dx), &xx),
        &bv_scale(&bv_mul(&dx, &dy), &(int(2) * xy)),
    );
    e = bv_add(&e, &bv_scale(&bv_mul(&dy, &dy), &yy));
    bv_add(&e, &vec![vec![-circ.r2.clone()]])
}

/// `p` reduced modulo the quadratic `e` in `t` (its leading coefficient a
/// nonzero constant): `[r0, r1]`.
fn reduce(p: &Bv, e: &Bv) -> [Poly; 2] {
    let lead = e[2][0].clone();
    let mut p = p.clone();
    while p.len() > 2 {
        let j = p.len() - 1;
        let c = pscale(&p[j], &(int(1) / &lead));
        let shifted: Bv = (0..j - 2)
            .map(|_| Vec::new())
            .chain(e.iter().map(|x| pmul(x, &c)))
            .collect();
        p = bv_add(&p, &bv_scale(&shifted, &int(-1)));
        // The top coefficient cancels exactly.
        debug_assert!(p.len() <= j);
    }
    let get = |j: usize| p.get(j).cloned().unwrap_or_default();
    [get(0), get(1)]
}

/// The resultant in `t` of the quadratic `e` and `r1 t + r0`.
fn resultant(e: &Bv, r: &[Poly; 2]) -> Poly {
    let get = |j: usize| e.get(j).cloned().unwrap_or_default();
    let [r0, r1] = r;
    trim(padd(
        &padd(
            &pmul(&get(2), &pmul(r0, r0)),
            &pscale(&pmul(&get(1), &pmul(r0, r1)), &int(-1)),
        ),
        &pmul(&get(0), &pmul(r1, r1)),
    ))
}

/// Where a sphere's circle of a surd radius, or on a basis of unequal
/// lengths, meets a torus (model `t`): each crossing's `(dx, dy)` place and
/// point in `Q(alpha)`; a tangency, or a crossing within the resolution
/// `res` of one, `Degenerate`.
pub(super) fn circ_torus(circ: &Circ, t: &Prism, res: f64) -> Result<EdgeMeet> {
    let ring = t.ring.as_ref().expect("a torus");
    let f = &t.f;
    let l = [f.local(&circ.c), f.local_dir(&circ.x), f.local_dir(&circ.y)];
    let (big, small) = (&ring.big, &ring.small);
    // The torus's gradient in local coordinates is `8 R r rho` on it, `rho
    // <= R + r`: the band of its function a displacement of `res` spans.
    let res_q = q(res);
    let delta = int(8) * big * small * (big + small + &res_q) * &res_q;
    let mut tangent = false;
    for [cr, sr] in rotations() {
        let e = circle_in_st(circ, &cr, &sr);
        let r = reduce(&torus_in_st(&l, big, small, &cr, &sr, &zero()), &e);
        let res0 = resultant(&e, &r);
        if res0.is_empty() {
            // The circle on the torus.
            return Ok(EdgeMeet::Along);
        }
        let Some(n) = count(&res0)? else {
            tangent = true;
            continue;
        };
        for k in [delta.clone(), -delta.clone()] {
            let rk = resultant(&e, &reduce(&torus_in_st(&l, big, small, &cr, &sr, &k), &e));
            if count(&rk)? != Some(n) {
                return Err(Error::Degenerate(
                    "a sphere's circle within the resolution of tangency to a torus (S9d.4c)",
                ));
            }
        }
        if n == 0 {
            return Ok(EdgeMeet::Points(Vec::new()));
        }
        let mut out = Vec::new();
        let mut ok = true;
        for root in roots(&res0)? {
            let g = Arc::new(Gen::new(res0.clone(), root));
            let s = K::generator(&g);
            let Some(inv) = peval_k(&r[1], &s).recip() else {
                ok = false;
                break;
            };
            let tt = peval_k(&r[0], &s).neg().mul(&inv);
            let dx = s.scale(&cr).sub(&tt.scale(&sr));
            let dy = s.scale(&sr).add(&tt.scale(&cr));
            let p: QV = [0, 1, 2].map(|j| {
                Qd::of(
                    K::Rat(circ.c[j].clone())
                        .add(&dx.scale(&circ.x[j]))
                        .add(&dy.scale(&circ.y[j])),
                )
            });
            out.push((Pos::Ang([Qd::of(dx), Qd::of(dy)]), p));
        }
        if ok {
            return Ok(EdgeMeet::Points(out));
        }
        tangent = true;
    }
    debug_assert!(tangent);
    Err(tangency())
}

// ------------------------------------------------------------ a part's rim

/// A tangency between a part's rim and the other surface.
fn rim_tangent() -> Error {
    Error::Degenerate("a torus part's rim tangent to the other input's surface (S9d.4c)")
}

/// A chart of the free angle whose antipode is no root of `A + sqrt(d) B`.
fn clear_chart(a: &Form, b: &Form, d: &R) -> Result<Chart> {
    let at = |c: &[R; 2]| Qd::new(a.value(c), b.value(c), d.clone()).sign();
    let mut dirs: Vec<[R; 2]> = [(1, 0), (0, 1), (-1, 0), (0, -1)]
        .iter()
        .map(|&(c, s)| [int(c), int(s)])
        .collect();
    for k in 1..40i64 {
        dirs.push(circle_point(&[zero(), zero()], &int(1), &(int(k) / int(7))));
    }
    for c in dirs {
        if at(&[-c[0].clone(), -c[1].clone()]) != Ordering::Equal {
            return Ok(Chart {
                c0: c[0].clone(),
                s0: c[1].clone(),
            });
        }
    }
    Err(rim_tangent())
}

/// Where a part's rim meets the other input's curved surface `far`: the
/// crossings of `G`'s form at the rim's fixed angle (see the module's
/// note), each with its place on the rim.
pub(super) fn rim_far(c: &TorusSec, far: &Far) -> Result<EdgeMeet> {
    let ring = c.ring();
    // The fixed angle's exact `(cos, sin)`, from the rim's point at the free
    // angle's `(1, 0)`: a wedge's `u` along `(l_u, l_v) / rho`, a segment's
    // `v` along `(rho - R, l_w) / r`.
    let p0 = c
        .at(&[int(1), zero()])
        .ok_or(limit("a torus part's rim without its point"))?;
    let l0 = c.f.local_q(&p0);
    let fixed: [Qd; 2] = if c.over_v {
        let inv = ring
            .rho(&l0)
            .recip()
            .ok_or(limit("a rim's point on the axis"))?;
        Ring::u_dir(&l0).map(|x| x.mul(&inv))
    } else {
        let inv = int(1) / &c.small;
        ring.v_dir(&l0).map(|x| x.scale(&inv))
    };
    let g = meeting_form(&c.f, &c.big, &c.small, far);
    let (fa, fb, d) = g.at_dir(!c.over_v, &fixed);
    let chart = clear_chart(&fa, &fb, &d)?;
    let (pa, pb) = (trim(fa.poly(&chart)), trim(fb.poly(&chart)));
    let surd = d != zero() && !pb.is_empty();
    if pa.is_empty() && !surd {
        return Ok(EdgeMeet::Along);
    }
    let norm = if surd {
        trim(padd(&pmul(&pa, &pa), &pscale(&pmul(&pb, &pb), &-d.clone())))
    } else {
        pa.clone()
    };
    if norm.len() <= 1 {
        return Ok(EdgeMeet::Points(Vec::new()));
    }
    let (sf, rs) = roots_repeated(&norm)?;
    let (da, db) = (pderiv(&pa), pderiv(&pb));
    let value = |a: &Poly, b: &Poly, t: &K| {
        Qd::parts(
            peval_k(a, t),
            if surd { peval_k(b, t) } else { K::Rat(zero()) },
            d.clone(),
        )
    };
    let mut out = Vec::new();
    for (root, repeated) in rs {
        let root: AlgebraicRoot = root;
        let gen = Arc::new(Gen::new(sf.clone(), root));
        let t = K::generator(&gen);
        if value(&pa, &pb, &t).sign() != Ordering::Equal {
            // The conjugate's crossing: the other ring's, or the opposite
            // half-plane's circle's.
            continue;
        }
        if repeated && value(&da, &db, &t).sign() == Ordering::Equal {
            return Err(rim_tangent());
        }
        // The free angle's `(cos, sin)`: `(1 - t^2, 2 t) / (1 + t^2)` turned
        // by the chart's base.
        let inv = t
            .mul(&t)
            .add(&K::Rat(int(1)))
            .recip()
            .ok_or(limit("a rim's crossing at infinity"))?;
        let (cc, ss) = (K::Rat(int(1)).sub(&t.mul(&t)), t.scale(&int(2)));
        let free = [
            Qd::of(cc.scale(&chart.c0).sub(&ss.scale(&chart.s0)).mul(&inv)),
            Qd::of(cc.scale(&chart.s0).add(&ss.scale(&chart.c0)).mul(&inv)),
        ];
        let (cu, cv) = if c.over_v {
            (&fixed, &free)
        } else {
            (&free, &fixed)
        };
        let rho = cv[0].scale(&c.small).add_r(&c.big);
        let l = [rho.mul(&cu[0]), rho.mul(&cu[1]), cv[1].scale(&c.small)];
        let x = qadd(
            &qv(&c.f.o),
            &qadd(
                &qadd(&qscale(&c.f.x, &l[0]), &qscale(&c.f.y, &l[1])),
                &qscale(&c.f.n, &l[2]),
            ),
        );
        out.push((Pos::Ang(c.place(&x)), x));
    }
    Ok(EdgeMeet::Points(out))
}
