//! S9d.2: spheres against prisms with arcs, and two spheres (REVIEW_NOTES.md,
//! S9d.2 refined). Two spheres meet in their radical plane's section. A
//! cylinder and a sphere meet where the cylinder's rulings meet the sphere:
//! over the cylinder's angle, the ruling's quadratic has a discriminant
//! positive all round (two rings over the angle, graphs of `Curve3::Meet`
//! with the sphere as the other quadric), negative all round (apart), or
//! changing sign (a loop, whose turning points need graphs over the height:
//! S9d.2b's, `OutOfDomain`). A prism's cap circle meets a sphere at the
//! roots of a quartic in its half-angle tangent (`algebraic.rs`); a
//! sphere's own circle (a rim, the split) meets a cylinder where `F0 + s F1
//! = 0`, `s` its radius over its basis's length, at algebraic points with
//! one surd, and another sphere on the spheres' radical plane.
use super::meet::{tangency, CylPair, EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use super::procedural::{other_of, other_sphere, MeetCrv, Other, Quartic};
use super::sphere::{plane_section, Circ};
use super::turned::{roots, square_sum, trim, Chart, Lin};
use crate::polynomial::real::IntPolynomial;
use crate::solid::split::zero;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::sync::Arc;

/// A sphere meeting a cylinder in a loop: S9d.2b's.
fn loop_later() -> Error {
    Error::OutOfDomain("a sphere meeting a cylinder in a loop (S9d.2b)")
}

/// A cylinder of an operand: its frame, circle centre and radius.
type Cyl<'a> = (&'a Affine, &'a P2, &'a R);

/// How a cylinder (of operand `k`) and a sphere meet: rings over the
/// cylinder's angle, or apart.
pub(super) fn sphere_cyl(k: usize, cyl: Cyl, c: &V, r: &R, res: f64) -> Result<CylPair> {
    let other = other_sphere(c, r);
    let (a, d) = super::turned::discriminant(cyl, &other);
    let chart = super::turned::negative_chart(&d)?;
    let test = chart.clone().unwrap_or(Chart {
        c0: int(1),
        s0: zero(),
    });
    super::turned::near_node(&a, &d, &test, res)?;
    let Some(chart) = chart else {
        // D >= 0 all round (and no double root): two rings.
        let piece = |plus| MeetCrv::new(k, cyl, &other, plus, None);
        return Ok(CylPair::Quartic(Box::new(Quartic {
            pieces: vec![piece(true), piece(false)],
            switches: Vec::new(),
        })));
    };
    if roots(&d.poly(&chart))?.is_empty() {
        return Ok(CylPair::Apart);
    }
    Err(loop_later())
}

/// Two spheres' meeting: their radical plane's section (`None` apart), the
/// same sphere refused.
pub(super) fn sphere_sphere(c1: &V, r1: &R, c2: &V, r2: &R) -> Result<Option<Circ>> {
    let m = sub(c2, c1);
    if is_zero(&m) {
        return Err(Error::Degenerate("two spheres about one centre"));
    }
    let (p0, m) = radical(c1, r1, c2, r2);
    plane_section(c1, r1, &p0, &m)
}

/// The radical plane of two spheres: a point and its normal.
pub(super) fn radical(c1: &V, r1: &R, c2: &V, r2: &R) -> (V, V) {
    let m = sub(c2, c1);
    let rhs = dot(c2, c2) - dot(c1, c1) + r1 * r1 - r2 * r2;
    let p0 = scale(&m, &(rhs / (int(2) * dot(&m, &m))));
    (p0, m)
}

/// Where a sphere's own circle meets a cylinder (`f`, `cy`, `ry`): points
/// `c + s (cos x + sin y)`, `s^2 = r2 / |x|^2` (a basis of equal lengths),
/// where `F0 + s F1 = 0`; S9d.2b's for another basis.
pub(super) fn circ_cylinder(circ: &Circ, f: &Affine, cy: &P2, ry: &R) -> Result<EdgeMeet> {
    circ_quadric(circ, &other_of(f, cy, ry))
}

fn circ_quadric(circ: &Circ, o: &Other) -> Result<EdgeMeet> {
    let (xx, yy) = (dot(&circ.x, &circ.x), dot(&circ.y, &circ.y));
    if xx != yy || dot(&circ.x, &circ.y) != zero() {
        return Err(Error::OutOfDomain(
            "a sphere's circle of unequal axes against a cylinder (S9d.2b)",
        ));
    }
    let sigma = &circ.r2 / &xx;
    // a_i + s L_i(cos, sin): F0 = sum a_i^2 + sigma sum L_i^2 - R^2, F1 = 2
    // sum a_i L_i.
    let a: Vec<R> = (0..o.g.len())
        .map(|i| dot(&o.g[i], &circ.c) - &o.e[i])
        .collect();
    let l: Vec<Lin> = (0..o.g.len())
        .map(|i| [zero(), dot(&o.g[i], &circ.x), dot(&o.g[i], &circ.y)])
        .collect();
    let mut f0 = square_sum(&l).scaled(&sigma);
    f0.k += a.iter().fold(zero(), |acc, x| acc + x * x) - &o.r * &o.r;
    let f1: Lin = [0, 1, 2].map(|j| {
        a.iter()
            .zip(&l)
            .fold(zero(), |acc, (ai, li)| acc + int(2) * ai * &li[j])
    });
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let s = rational_sqrt(&sigma);
    // F0 times (1 + t^2)^2, F1 times (1 + t^2).
    let p0 = trim(f0.poly(&chart));
    let [cn, sn] = chart.numerators();
    let w = vec![int(1), zero(), int(1)];
    let p1 = trim(padd(
        &padd(&pscale(&w, &f1[0]), &pscale(&cn, &f1[1])),
        &pscale(&sn, &f1[2]),
    ));
    // The equation's polynomial: F0 + s F1 (s rational), or F0^2 - sigma
    // F1^2 (its roots of the right signs kept).
    let eq = match &s {
        Some(s) => trim(padd(&p0, &pmul(&pscale(&p1, s), &w))),
        None => trim(padd(
            &pmul(&p0, &p0),
            &pscale(&pmul(&pmul(&p1, &p1), &pmul(&w, &w)), &-sigma.clone()),
        )),
    };
    if eq.is_empty() {
        return Err(tangency());
    }
    let point = |cs: &[K; 2]| -> QV {
        let dir = [0, 1, 2].map(|j| cs[0].scale(&circ.x[j]).add(&cs[1].scale(&circ.y[j])));
        [0, 1, 2].map(|j| match &s {
            Some(s) => Qd::of(K::Rat(circ.c[j].clone()).add(&dir[j].scale(s))),
            None => Qd::parts(K::Rat(circ.c[j].clone()), dir[j].clone(), sigma.clone()),
        })
    };
    let place = |cs: &[K; 2]| -> [Qd; 2] {
        match &s {
            Some(s) => [Qd::of(cs[0].scale(s)), Qd::of(cs[1].scale(s))],
            None => [
                Qd::parts(K::Rat(zero()), cs[0].clone(), sigma.clone()),
                Qd::parts(K::Rat(zero()), cs[1].clone(), sigma.clone()),
            ],
        }
    };
    let mut out = Vec::new();
    // The chart's antipode (-1, 0): where the equation drops degrees.
    let anti = [K::Rat(int(-1)), K::Rat(zero())];
    let f0a = f0.value(&[int(-1), zero()]);
    let f1a = &f1[0] - &f1[1];
    let at_anti = match &s {
        Some(s) => f0a.clone() + s * &f1a == zero(),
        None => &f0a * &f0a == &sigma * &f1a * &f1a && (f0a == zero() || sign(&f0a) != sign(&f1a)),
    };
    if at_anti {
        out.push((Pos::Ang(place(&anti)), point(&anti)));
    }
    let deg_ok = |r: &crate::polynomial::real::AlgebraicRoot| -> bool {
        match &s {
            Some(_) => true,
            None => {
                let (g0, g1) = (
                    r.sign_polynomial(&IntPolynomial::from_rationals(&p0)),
                    r.sign_polynomial(&IntPolynomial::from_rationals(&p1)),
                );
                g0 == Ordering::Equal || (g1 != Ordering::Equal && g0 != g1)
            }
        }
    };
    for root in roots(&eq)? {
        if !deg_ok(&root) {
            continue;
        }
        let g = Arc::new(Gen {
            poly: eq.clone(),
            root,
        });
        let t = K::generator(&g);
        let den = t.mul(&t).add(&K::Rat(int(1)));
        let inv = den
            .recip()
            .ok_or(Error::ComputationLimit("a circle's crossing at infinity"))?;
        let cs = [
            K::Rat(int(1)).sub(&t.mul(&t)).mul(&inv),
            t.scale(&int(2)).mul(&inv),
        ];
        out.push((Pos::Ang(place(&cs)), point(&cs)));
    }
    Ok(EdgeMeet::Points(out))
}

fn padd(a: &[R], b: &[R]) -> Vec<R> {
    let n = a.len().max(b.len());
    trim(
        (0..n)
            .map(|i| {
                a.get(i).cloned().unwrap_or_else(zero) + b.get(i).cloned().unwrap_or_else(zero)
            })
            .collect(),
    )
}

fn pmul(a: &[R], b: &[R]) -> Vec<R> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trim(out)
}

fn pscale(a: &[R], k: &R) -> Vec<R> {
    trim(a.iter().map(|x| x * k).collect())
}

/// Where a sphere's own circle (on the sphere `(c1, r1)`) meets another
/// sphere: on the two spheres' radical plane.
pub(super) fn circ_sphere(circ: &Circ, c1: &V, r1: &R, c2: &V, r2: &R) -> Result<EdgeMeet> {
    if is_zero(&sub(c2, c1)) {
        return Err(Error::Degenerate("two spheres about one centre"));
    }
    let (p0, m) = radical(c1, r1, c2, r2);
    Ok(match circ.meet_plane(&p0, &m)? {
        None => EdgeMeet::Along,
        Some(cs) => EdgeMeet::Points(
            cs.into_iter()
                .map(|e| {
                    let x = qadd(
                        &qv(&circ.c),
                        &qadd(&qscale(&circ.x, &e[0]), &qscale(&circ.y, &e[1])),
                    );
                    (Pos::Ang(e), x)
                })
                .collect(),
        ),
    })
}
