//! S9f.3b: spline walls against cones and frustums in any position
//! (REVIEW_NOTES.md, "S9f.3b refined, before its code").
//!
//! Along a spline wall's ruling at the run parameter `tau`, `X = o +
//! S_x(tau) x + S_y(tau) y + w n`, the cone's function on its exact model
//! (`procedural::other_cone`: `u^2 + v^2 - (b + k w_c)^2` in its frame's
//! exact rows) is `A w^2 + 2 B(tau) w + C(tau)`: S9f.2b's quadratic with the
//! radius row of negative sign (`spline_crossing::terms`), `A = q_u^2 +
//! q_v^2 - k^2 q_w^2` one constant for the pair (`q` the prism's axis in the
//! cone's rows). For `A > 0` the meeting is S9f.3a's (branches over the run, loops'
//! graphs over the height about turning points inside both faces, their
//! switches); for `A < 0` every ruling meets the double cone once on each
//! nappe, the discriminant positive but where a ruling passes the apex, and
//! the two branches run over the whole run, the other nappe's beyond the
//! apex and outside the cone's face; `A = 0` (the prism's axis along a
//! generatrix direction) is `Degenerate`, and so is the apex, real or
//! virtual, on the spline wall's surface, both checked first.
//!
//! Vertices: a curve over a spline segment (a cap edge or a crease) meets
//! the cone at the roots of its function along it, degree `2 p`
//! (`wallcrv_cone`); the spline prism's vertical edges meet it by S9d.3a's
//! `line_cone`; a rim (the conic `c + a cos + b sin` of rational radius in
//! the cone's end plane) meets a spline wall along its plane's crease in
//! the rim's own elliptic cylinder, degree `2 p` (`conic_wall`), or, in a
//! plane holding the wall's axis direction, at S9f.2b.2's points in the
//! rim's half-angle chart (`spline_crossing::tower_points`).
use super::meet::{CylPair, EdgeMeet};
use super::model::{Affine, Prism};
use super::num::*;
use super::procedural::{other_cone, Other};
use super::spline_crossing::{conic_crease, meeting_with, terms, tower_points, Partner};
use super::spline_walls::{combine, BArc, Roots, SplineSeg, WallCrv};
use crate::solid::split::zero;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::sync::Arc;

/// A spline wall's meeting with a cone's wall (`sm`'s face `sf` holding
/// segment `seg`, `km`'s cone wall `kf`, its radius `b + k w` on its frame):
/// S9f.2b's meeting over the cone's signed rows, after the decisions'
/// refusals of `A = 0` and of the apex on the wall's surface.
#[allow(clippy::too_many_arguments)]
pub(super) fn meeting(
    sm: &Prism,
    sf: usize,
    seg: &Arc<SplineSeg>,
    km: &Prism,
    kf: usize,
    b: &R,
    k: &R,
) -> Result<CylPair> {
    let other = other_cone(&km.f, b, k);
    let (rows, _) = terms(&other);
    let a = rows.iter().fold(zero(), |acc, (g, _, s)| {
        let q = dot(g, &sm.f.n);
        acc + s * &q * &q
    });
    if a == zero() {
        return Err(Error::Degenerate("a spline wall along a cone's ruling"));
    }
    // The apex (a frustum's virtual one) on the wall's surface: its profile
    // point on the segment.
    let apex = sm.f.local(&km.f.point(&zero(), &zero(), &(-b / k)));
    if seg
        .locate(&[Qd::rat(apex[0].clone()), Qd::rat(apex[1].clone())])
        .is_some()
    {
        return Err(Error::Degenerate(
            "a cone's apex on the other input's surface",
        ));
    }
    meeting_with(sm, sf, seg, km, kf, other, Partner::Cone)
}

/// Where a curve over a spline segment (a cap edge or a crease) meets a
/// cone: the roots of its function along it, degree `2 p`.
pub(super) fn wallcrv_cone(curve: &WallCrv, kf: &Affine, b: &R, k: &R) -> Result<EdgeMeet> {
    super::spline_crossing::wallcrv_quadric(curve, &other_cone(kf, b, k))
}

/// Where a cone's rim, the conic `cc + a cos + b sin` (a circle of rational
/// radius on the cone's stored axes, not exactly orthogonal), meets a spline
/// wall (segment `seg` on frame `sf`): along its plane's crease in the
/// rim's own elliptic cylinder (the rows dual to `a` and `b` in the plane,
/// its function `alpha^2 + beta^2 - 1`), or in a plane holding the wall's
/// axis direction by S9f.2b.2's half-angle chart.
pub(super) fn conic_wall(cc: &V, a: &V, b: &V, sf: &Affine, seg: &SplineSeg) -> Result<EdgeMeet> {
    let m = cross(a, b);
    let mn = dot(&m, &sf.n);
    if mn == zero() {
        let k = dot(&m, &sub(&sf.o, cc));
        let (mx, my) = (dot(&m, &sf.x), dot(&m, &sf.y));
        return match seg.roots_of(&|arc: &BArc| combine(&k, &mx, &arc.x, &my, &arc.y))? {
            Roots::At(rs) if rs.is_empty() => Ok(EdgeMeet::None),
            Roots::At(_) => tower_points(cc, a, b, sf, seg, Partner::Cone),
            Roots::Along | Roots::Partly => {
                Err(Error::OutOfDomain("a spline wall along a plane (S9f)"))
            }
        };
    }
    // `alpha = g1 . (X - cc)`, `beta = g2 . (X - cc)`: `g1 . a = g2 . b = 1`,
    // `g1 . b = g2 . a = 0`, both normal to `m`.
    let det = dot(a, &cross(b, &m));
    let inv = int(1) / &det;
    let g1 = scale(&cross(b, &m), &inv);
    let g2 = scale(&cross(&m, a), &inv);
    let own = Other {
        e: vec![dot(&g1, cc), dot(&g2, cc)],
        g: vec![g1, g2],
        r: int(1),
        h: [zero(), zero(), zero()],
        eh: zero(),
        t: zero(),
    };
    conic_crease(&own, cc, a, b, sf, seg)
}
