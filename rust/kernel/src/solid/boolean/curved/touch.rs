//! S9e.4b.4d.2c.1a: an edge of a face that touches the other input's face
//! along a circle (REVIEW_NOTES.md, "S9e.4b.4d.2c refined" and "S9e.4b.4d.2c.1
//! refined"). In a turned body's own tree a whole torus touches a coaxial
//! plane, cylinder or cone of the other input along a circle that is an edge
//! of the torus's input (`graph.rs`'s `along_touch`); every other edge of the
//! two faces meets the other's surface where it crosses that circle, a double
//! root of its meeting with the other surface there (a cylinder's or a cone's
//! ruling against the torus, a torus's meridian against the plane, cylinder
//! or cone), and wherever else it crosses it transversally. A double root
//! off every touching circle is the tangency it is.
use super::meet::{tangency, EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use super::procedural::Other;
use super::torus::TouchRing;
use super::turned::{roots_repeated, square_sum, trim, Chart, Lin};
use crate::solid::split::zero;
use crate::{Error, Result};
use std::sync::Arc;

/// Where an edge's curve on a face touching face `fy` of `py` along the
/// circles `rings` meets that face's surface, or `None` where the curve is
/// none of those it decides (the meeting is then found as any edge's).
pub(super) fn edge_touching(
    curve: &Crv,
    py: &Prism,
    fy: usize,
    rings: &[&TouchRing],
) -> Option<Result<EdgeMeet>> {
    let on_ring = |x: &QV| rings.iter().any(|r| r.on(x));
    match (curve, &py.faces[fy].surf) {
        // A ruling against the torus it touches.
        (Crv::Line { p, d }, Surf::Torus) => Some((|| {
            let ring = py.ring.as_ref().expect("a torus");
            let poly = super::torus::line_quartic(p, d, &py.f, ring)?;
            let mut out = Vec::new();
            if poly.len() < 2 {
                return Ok(EdgeMeet::Points(out));
            }
            let (sf, rs) = roots_repeated(&poly)?;
            for (root, repeated) in rs {
                let g = Arc::new(Gen::new(sf.clone(), root));
                let t = Qd::of(K::generator(&g));
                let x = qadd(p, &qscale(d, &t));
                if repeated && !on_ring(&x) {
                    return Err(tangency());
                }
                out.push((Pos::T(t), x));
            }
            Ok(EdgeMeet::Points(out))
        })()),
        // A torus's meridian or parallel against the plane it touches.
        (Crv::Conic { c, a, b }, Surf::Plane { p: p0, m }) => {
            let alpha = dot(m, a);
            let beta = dot(m, b);
            let gamma = dot(m, &sub(p0, c));
            let ab = &alpha * &alpha + &beta * &beta;
            if ab == zero() || ab != &gamma * &gamma {
                return None;
            }
            // Tangent: the one point `(alpha, beta) gamma / (alpha^2 +
            // beta^2)`.
            let cs = [
                Qd::rat(&alpha * &gamma / &ab),
                Qd::rat(&beta * &gamma / &ab),
            ];
            let x = super::meet::conic_point(c, a, b, &cs);
            Some(if on_ring(&x) {
                Ok(EdgeMeet::Points(vec![(Pos::Ang(cs), x)]))
            } else {
                Err(tangency())
            })
        }
        // A torus's meridian or parallel against the cylinder or cone it
        // touches.
        (Crv::Conic { c, a, b }, Surf::Cyl { c: cy, r: ry, .. }) => Some(circle_quadric(
            c,
            a,
            b,
            &super::procedural::other_of(&py.f, cy, ry),
            &on_ring,
        )),
        (Crv::Conic { c, a, b }, Surf::Cone { b: rb, k }) => Some(circle_quadric(
            c,
            a,
            b,
            &super::procedural::other_cone(&py.f, rb, k),
            &on_ring,
        )),
        _ => None,
    }
}

/// `algebraic::circle_quadric` with each double root on a touching circle
/// taken as a point of the meeting (a tangency elsewhere).
fn circle_quadric(
    c: &V,
    a: &V,
    b: &V,
    o: &Other,
    on_ring: &dyn Fn(&QV) -> bool,
) -> Result<EdgeMeet> {
    let lin: Vec<Lin> = (0..o.g.len())
        .map(|i| [dot(&o.g[i], c) - &o.e[i], dot(&o.g[i], a), dot(&o.g[i], b)])
        .collect();
    let form = square_sum(&lin).sub(&square_sum(&[o.radius_lin(c, a, b)]));
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let p = trim(form.poly(&chart));
    if p.is_empty() {
        return Err(tangency());
    }
    let point = |cs: &[Qd; 2]| -> QV {
        [0, 1, 2].map(|j| {
            Qd::rat(c[j].clone())
                .add(&cs[0].scale(&a[j]))
                .add(&cs[1].scale(&b[j]))
        })
    };
    let mut out = Vec::new();
    if form.value(&[int(-1), zero()]) == zero() {
        let cs = [Qd::rat(int(-1)), Qd::rat(zero())];
        // A root at the chart's antipode: repeated where the polynomial's
        // degree falls by two or more there.
        let x = point(&cs);
        if p.len() <= 3 && !on_ring(&x) {
            return Err(tangency());
        }
        out.push((Pos::Ang(cs), x));
    }
    if p.len() >= 2 {
        let (sf, rs) = roots_repeated(&p)?;
        for (root, repeated) in rs {
            let g = Arc::new(Gen::new(sf.clone(), root));
            let t = K::generator(&g);
            let den = t.mul(&t).add(&K::Rat(int(1)));
            let inv = den
                .recip()
                .ok_or(Error::ComputationLimit("a circle's crossing at infinity"))?;
            let cos = K::Rat(int(1)).sub(&t.mul(&t)).mul(&inv);
            let sin = t.scale(&int(2)).mul(&inv);
            let cs = [Qd::of(cos), Qd::of(sin)];
            let x = point(&cs);
            if repeated && !on_ring(&x) {
                return Err(tangency());
            }
            out.push((Pos::Ang(cs), x));
        }
    }
    Ok(EdgeMeet::Points(out))
}
