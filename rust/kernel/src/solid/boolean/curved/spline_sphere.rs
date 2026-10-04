//! S9f.3a: spline walls against spheres, caps and zones in any position
//! (REVIEW_NOTES.md, "S9f.3 refined, before its code").
//!
//! Along a spline wall's ruling at the run parameter `tau`, `X = o +
//! S_x(tau) x + S_y(tau) y + w n`, the sphere's function `|X - c|^2 - r^2`
//! (its stored centre and radius: `procedural::other_sphere`, the world's
//! three rows) is `A w^2 + 2 B(tau) w + C(tau)` with `A = n . n` (positive:
//! every ruling), `B` of degree `p` and `C` of degree `2 p` on each Bézier
//! arc, and `D = B^2 - A C = A r^2 - |(P - c) x n|^2`: S9f.2b's quadratic with
//! three rows instead of a cylinder's two. The meeting is therefore S9f.2b's
//! (`spline_crossing::meeting_with`): branches over the run between the
//! roots of `D`, graphs over the height about each turning point inside both
//! faces, their switches vertices; a turning point on the hemispheres' split
//! (no edge of the input) tries another seam.
//!
//! Vertices: a curve over a spline segment (a cap edge or a crease) meets
//! the sphere at the roots of its function along it, degree `2 p`
//! (`wallcrv_sphere`); the spline prism's vertical edges meet it at
//! quadratic surds (S9d.1's `line_sphere`); a sphere's circle (`Circ`: a
//! rim of a cap or a zone, the split's great circle, each on its own sphere)
//! meets a spline wall along its plane's crease on the wall, `w = h0 + h1
//! S_x + h2 S_y` in the sphere's function, degree `2 p`, its place on the
//! circle read off rationally (`circ_wall`); in a plane holding the wall's
//! axis direction, the plane's generatrices meet the circle in a tower
//! `Q(alpha)(sqrt(delta))`, found instead by a primitive element: each arc's
//! implicit equation at the circle's projection `l0 + dx lx + dy ly`,
//! reduced by the circle's equation `dx^2 |x|^2 + dy^2 |y|^2 = r2` to `E(dx)
//! + dy O(dx)`, vanishes where `E^2 - (r2 - |x|^2 dx^2) / |y|^2 O^2` does
//! (degree at most `2 p` in `dx`), and there `dy = -E / O`, in the one field
//! `Q(dx)`; where `O` vanishes at a root (two of the circle's points sharing
//! `dx`, as a great circle's whose `x` is the wall's axis) the roles of `dx`
//! and `dy` are swapped (`tower_points`).
use super::meet::{CylPair, EdgeMeet, Pos};
use super::model::{Affine, Prism};
use super::num::*;
use super::procedural::other_sphere;
use super::sphere::{Ball, Circ};
use super::spline_crossing::{along_curve, meeting_with, padd, pmul, pscale, Partner};
use super::spline_walls::{combine, qpoint, trim, BArc, Roots, SplineSeg, WallCrv};
use super::turned::roots_repeated;
use crate::solid::split::zero;
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::sync::Arc;

fn tangent_curve() -> Error {
    Error::Degenerate("a spline edge of one input tangent to a face of the other")
}

/// A spline wall's meeting with a sphere (`sm`'s face `sf` holding segment
/// `seg`, `bm`'s hemisphere `bf` on the sphere `c`, `r`): S9f.2b's meeting
/// over the sphere's rows.
#[allow(clippy::too_many_arguments)]
pub(super) fn meeting(
    sm: &Prism,
    sf: usize,
    seg: &Arc<SplineSeg>,
    bm: &Prism,
    bf: usize,
    c: &V,
    r: &R,
) -> Result<CylPair> {
    meeting_with(sm, sf, seg, bm, bf, other_sphere(c, r), Partner::Sphere)
}

/// Where a curve over a spline segment (a cap edge or a crease) meets a
/// sphere: the roots of its function along it, degree `2 p`.
pub(super) fn wallcrv_sphere(curve: &WallCrv, c: &V, r: &R) -> Result<EdgeMeet> {
    super::spline_crossing::wallcrv_quadric(curve, &other_sphere(c, r))
}

/// Where a sphere's circle (on the sphere `ball`) meets a spline wall
/// (segment `seg` on frame `sf`): along its plane's crease on the wall, or
/// in a plane holding the wall's axis direction by a primitive element of
/// the tower (`tower_points`).
pub(super) fn circ_wall(
    circ: &Circ,
    ball: &Ball,
    sf: &Affine,
    seg: &SplineSeg,
) -> Result<EdgeMeet> {
    let m = circ.normal();
    let mn = dot(&m, &sf.n);
    let k = dot(&m, &sub(&sf.o, &circ.c));
    let (mx, my) = (dot(&m, &sf.x), dot(&m, &sf.y));
    if mn == zero() {
        // The plane holds the wall's axis direction: its generatrices on the
        // wall; none when the plane misses the segment.
        return match seg.roots_of(&|arc: &BArc| combine(&k, &mx, &arc.x, &my, &arc.y))? {
            Roots::At(rs) if rs.is_empty() => Ok(EdgeMeet::None),
            Roots::At(_) => tower_points(circ, sf, seg),
            Roots::Along | Roots::Partly => {
                Err(Error::OutOfDomain("a spline wall along a plane (S9f)"))
            }
        };
    }
    let inv = int(-1) / &mn;
    let h = [&k * &inv, &mx * &inv, &my * &inv];
    let other = other_sphere(&ball.c, &ball.r);
    let rs = match seg.roots_of(&|arc: &BArc| along_curve(sf, arc, &other, &h))? {
        Roots::At(rs) => rs,
        Roots::Along | Roots::Partly => return Err(tangent_curve()),
    };
    let mut out = Vec::new();
    for root in rs {
        if root.mult > 1 {
            return Err(tangent_curve());
        }
        let uv = seg.point(&root.tau);
        let w = uv[0].scale(&h[1]).add(&uv[1].scale(&h[2])).add_r(&h[0]);
        let x = qpoint(sf, &uv[0], &uv[1], &w);
        out.push((Pos::Ang(circ.place(&x)), x));
    }
    Ok(EdgeMeet::Points(out))
}

/// A bivariate polynomial: `p[a][b]` the coefficient of `X^a Y^b`.
type Bi = Vec<Vec<R>>;

fn bmul(p: &Bi, q: &Bi) -> Bi {
    let (na, nb) = (p.len() + q.len() - 1, p[0].len() + q[0].len() - 1);
    let mut out = vec![vec![zero(); nb]; na];
    for (a1, row1) in p.iter().enumerate() {
        for (b1, x) in row1.iter().enumerate() {
            if *x == zero() {
                continue;
            }
            for (a2, row2) in q.iter().enumerate() {
                for (b2, y) in row2.iter().enumerate() {
                    out[a1 + a2][b1 + b2] += x * y;
                }
            }
        }
    }
    out
}

/// `c0 + c1 X + c2 Y`.
fn blinear(c0: &R, c1: &R, c2: &R) -> Bi {
    vec![vec![c0.clone(), c2.clone()], vec![c1.clone(), zero()]]
}

/// An arc's implicit equation `f(u, v)` at `(u, v) = l0 + X lx + Y ly`.
fn substituted(terms: &[(usize, usize, R)], l0: &V, lx: &V, ly: &V) -> Bi {
    let u = blinear(&l0[0], &lx[0], &ly[0]);
    let v = blinear(&l0[1], &lx[1], &ly[1]);
    let d = terms
        .iter()
        .map(|(i, j, _)| i.max(j))
        .max()
        .copied()
        .unwrap_or(0);
    let powers = |p: &Bi| {
        let mut out: Vec<Bi> = vec![vec![vec![int(1)]]];
        for k in 1..=d {
            let next = bmul(&out[k - 1], p);
            out.push(next);
        }
        out
    };
    let (pu, pv) = (powers(&u), powers(&v));
    // Every term's product is square, of side its total degree plus one.
    let n = terms.iter().map(|(i, j, _)| i + j).max().unwrap_or(0) + 1;
    let mut out: Bi = vec![vec![zero(); n]; n];
    for (i, j, c) in terms {
        let t = bmul(&pu[*i], &pv[*j]);
        for (a, row) in t.iter().enumerate() {
            for (b, x) in row.iter().enumerate() {
                out[a][b] += x * c;
            }
        }
    }
    out
}

fn keval(p: &[R], t: &K) -> K {
    p.iter()
        .rev()
        .fold(K::Rat(zero()), |acc, c| acc.mul(t).add(&K::Rat(c.clone())))
}

/// A circle's points on an arc, `(dx, dy)`, each with whether it is
/// repeated.
type Points = Vec<([K; 2], bool)>;

/// The circle's points on one arc with `main` (0: `dx`, 1: `dy`) the
/// primitive element: `[(dx, dy), repeated]`, or `None` where two of the
/// circle's points share a root (`O` vanishing there).
fn eliminate(g: &Bi, main: usize, xx: &R, yy: &R, r2: &R) -> Result<Option<Points>> {
    // `G_k(s)`: the coefficients of the other variable's `k`-th power, as
    // polynomials in the main one `s`.
    let (na, nb) = (g.len(), g[0].len());
    let (nm, no) = if main == 0 { (na, nb) } else { (nb, na) };
    let coef = |m: usize, o: usize| if main == 0 { &g[m][o] } else { &g[o][m] };
    let parts: Vec<Vec<R>> = (0..no)
        .map(|o| trim((0..nm).map(|m| coef(m, o).clone()).collect()))
        .collect();
    // The other variable's square `rho(s)`.
    let (own, theirs) = if main == 0 { (xx, yy) } else { (yy, xx) };
    let rho = trim(vec![r2 / theirs, zero(), -(own / theirs)]);
    let mut pow = vec![int(1)];
    let (mut even, mut odd) = (Vec::new(), Vec::new());
    for (o, p) in parts.iter().enumerate() {
        if o >= 2 && o % 2 == 0 {
            pow = pmul(&pow, &rho);
        }
        let term = pmul(p, &pow);
        if o % 2 == 0 {
            even = padd(&even, &term);
        } else {
            odd = padd(&odd, &term);
        }
    }
    let res = padd(
        &pmul(&even, &even),
        &pscale(&pmul(&rho, &pmul(&odd, &odd)), &int(-1)),
    );
    if res.is_empty() {
        return Ok(None);
    }
    if res.len() == 1 {
        return Ok(Some(Vec::new()));
    }
    let (sf_poly, roots) = roots_repeated(&res)?;
    let mut out = Vec::new();
    for (root, repeated) in roots {
        let s = match root.rational_value() {
            Some(x) => K::Rat(x.clone()),
            None => K::generator(&Arc::new(Gen::new(sf_poly.clone(), root))),
        };
        let (e, o) = (keval(&even, &s), keval(&odd, &s));
        let Some(inv) = o.recip() else {
            return Ok(None);
        };
        let t = e.mul(&inv).neg();
        out.push((if main == 0 { [s, t] } else { [t, s] }, repeated));
    }
    Ok(Some(out))
}

/// S9f.3a: where a sphere's circle in a plane holding the wall's axis
/// direction meets the spline wall (segment `seg` on frame `sf`): each
/// arc's implicit equation at the circle's projection, reduced by the
/// circle's equation, its roots in one field (see the module's notes); a
/// root whose point lies off the segment (another branch of the implicit
/// curve) is dropped, one found on two arcs (a knot) kept once; a repeated
/// root on the segment is the circle tangent to a generatrix there.
fn tower_points(circ: &Circ, sf: &Affine, seg: &SplineSeg) -> Result<EdgeMeet> {
    let (l0, lx, ly) = (
        sf.local(&circ.c),
        sf.local_dir(&circ.x),
        sf.local_dir(&circ.y),
    );
    let (xx, yy) = (dot(&circ.x, &circ.x), dot(&circ.y, &circ.y));
    let mut out: Vec<(Pos, QV)> = Vec::new();
    for arc in &seg.arcs {
        let g = substituted(&arc.implicit().terms(), &l0, &lx, &ly);
        let mut found = None;
        for main in [0, 1] {
            if let Some(points) = eliminate(&g, main, &xx, &yy, &circ.r2)? {
                found = Some(points);
                break;
            }
        }
        let Some(points) = found else {
            return Err(Error::ComputationLimit(
                "a sphere's circle's points on a spline wall sharing both coordinates (S9f.3a)",
            ));
        };
        for ([dx, dy], repeated) in points {
            let place = [Qd::of(dx), Qd::of(dy)];
            let x: QV = std::array::from_fn(|i| {
                place[0]
                    .scale(&circ.x[i])
                    .add(&place[1].scale(&circ.y[i]))
                    .add_r(&circ.c[i])
            });
            let l = sf.local_q(&x);
            if seg.locate(&[l[0].clone(), l[1].clone()]).is_none() {
                continue;
            }
            if repeated {
                return Err(Partner::Sphere.edge());
            }
            if !out.iter().any(|(_, y)| qv_eq(y, &x)) {
                out.push((Pos::Ang(place), x));
            }
        }
    }
    Ok(EdgeMeet::Points(out))
}
