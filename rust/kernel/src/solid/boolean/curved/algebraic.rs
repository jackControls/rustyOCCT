//! S9c.2b.2: vertices at algebraic points (REVIEW_NOTES.md, S9c.2b.2).
//! Where two cylinders' section crosses a cap's circle the vertex lies on
//! that circle at a half-angle tangent `alpha`, a root of the circle's
//! quartic against the other cylinder: its coordinates lie in `Q(alpha)`
//! (`num.rs`'s surds over a generator). Parallel cylinders not circular in
//! a common measure meet in generatrices through the first's circle's
//! crossings with the other cylinder, algebraic the same way.
use super::meet::{tangency, CylPair};
use super::model::*;
use super::num::*;
use super::procedural::other_of;
use super::turned::{roots, square_sum, trim, Chart, Lin};
use crate::solid::split::zero;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::sync::Arc;

/// Where the circle `c + a cos + b sin` meets a cylinder: each crossing's
/// `(cos, sin)` place and point, exact (`Q(alpha)` for a root `alpha` of
/// the circle's quartic in its half-angle tangent; the antipode of `(1,
/// 0)` rational). A tangency is `Degenerate`.
pub(super) fn circle_points(
    c: &V,
    a: &V,
    b: &V,
    f: &Affine,
    cy: &P2,
    ry: &R,
) -> Result<Vec<(super::meet::Pos, QV)>> {
    circle_quadric(c, a, b, &other_of(f, cy, ry))
}

/// Where the circle `c + a cos + b sin` meets a quadric (a cylinder, or a
/// sphere: S9d.2).
pub(super) fn circle_quadric(
    c: &V,
    a: &V,
    b: &V,
    o: &super::procedural::Other,
) -> Result<Vec<(super::meet::Pos, QV)>> {
    let lin: Vec<Lin> = (0..o.g.len())
        .map(|i| [dot(&o.g[i], c) - &o.e[i], dot(&o.g[i], a), dot(&o.g[i], b)])
        .collect();
    // Less the other's radius term squared (a cone's varies, S9d.3b).
    let form = square_sum(&lin).sub(&square_sum(&[o.radius_lin(c, a, b)]));
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let p = trim(form.poly(&chart));
    if p.is_empty() {
        // The circle on the cylinder: taken as S9c.1 takes an edge along
        // a surface, never for a cap's circle against a crossing cylinder.
        return Err(tangency());
    }
    let point = |cs: [Qd; 2]| -> QV {
        [0, 1, 2].map(|j| {
            Qd::rat(c[j].clone())
                .add(&cs[0].scale(&a[j]))
                .add(&cs[1].scale(&b[j]))
        })
    };
    let mut out = Vec::new();
    if form.value(&[int(-1), zero()]) == zero() {
        let cs = [Qd::rat(int(-1)), Qd::rat(zero())];
        out.push((super::meet::Pos::Ang(cs.clone()), point(cs)));
    }
    for root in roots(&p)? {
        let g = Arc::new(Gen {
            poly: p.clone(),
            root,
        });
        let t = K::generator(&g);
        let den = t.mul(&t).add(&K::Rat(int(1)));
        let inv = den
            .recip()
            .ok_or(Error::ComputationLimit("a circle's crossing at infinity"))?;
        let cos = K::Rat(int(1)).sub(&t.mul(&t)).mul(&inv);
        let sin = t.scale(&int(2)).mul(&inv);
        let cs = [Qd::of(cos), Qd::of(sin)];
        out.push((super::meet::Pos::Ang(cs.clone()), point(cs)));
    }
    Ok(out)
}

/// A cylinder of an operand: its frame, circle centre and radius.
type Cyl<'a> = (&'a Affine, &'a P2, &'a R);

/// Two parallel cylinders (`x` of operand 0, `y` of operand 1) not
/// circular in a common measure: the generatrices of `x` through its
/// circle's crossings with `y` (both cylinders hold them: their axes are
/// parallel exactly).
pub(super) fn parallel(x: Cyl, y: Cyl) -> Result<CylPair> {
    let (f, c, r) = x;
    let centre = f.point(&c[0], &c[1], &zero());
    let points = circle_points(&centre, &scale(&f.x, r), &scale(&f.y, r), y.0, y.1, y.2)?;
    if points.is_empty() {
        return Ok(CylPair::Apart);
    }
    Ok(CylPair::Lines(
        points.into_iter().map(|(_, p)| p).collect(),
        f.n.clone(),
    ))
}
