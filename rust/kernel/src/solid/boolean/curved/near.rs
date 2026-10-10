//! S9d's near misses at cylinders, cones and planes (REVIEW_NOTES.md, "A
//! ball within the resolution of a cylinder or cone, and edges within it of
//! a face"): `graph.rs`'s `near_miss` holds a sphere against a plane or a
//! sphere, and `edge_near_misses` an edge or a vertex against a sphere;
//! here a sphere against a cylinder or a cone face (`sphere_quadrics`), an
//! edge's line, conic, circle or other curve against a plane, cylinder or
//! cone face (`edge_faces`), and the points of other curves nearest a
//! sphere for `edge_near_misses` (`nearest_on_run`), by the same rule:
//! within the resolution of tangency, crossing is `Degenerate`, and missing
//! where the gap lies outside either input; a gap inside both inputs
//! evaluates.
use super::graph::{
    boxes_meet, edge_places, intersect, one_surface, place, sphere_toward, strictly_within,
};
use super::meet::Pos;
use super::model::*;
use super::num::*;
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::f64::consts::{PI, TAU};

/// A sphere within the resolution of tangency to a cylinder or cone face.
pub(super) const NEAR_QUADRIC: &str =
    "a sphere within the resolution of tangency to a cylinder or cone (S9d.2)";
/// An edge within the resolution of tangency to a face of the other input.
pub(super) const NEAR_FACE: &str = "an edge within the resolution of tangency to a face (S9d.2)";

/// A face's surface: a plane, or a cylinder or a cone on its model's frame
/// (its local coordinates `(u, v, w)`).
enum Kind<'a> {
    Plane { p: &'a V, m: &'a V },
    Cyl { f: &'a Affine, c: &'a P2, r: &'a R },
    Cone { f: &'a Affine, b: &'a R, k: &'a R },
}

/// A face's surface and its data in binary64 (for searches only): a
/// plane's point and normal, or the frame's map to local coordinates and
/// a cylinder's centre and radius or a cone's `b` and `k`.
struct Quad<'a> {
    kind: Kind<'a>,
    o: [f64; 3],
    rows: [[f64; 3]; 3],
    data: [f64; 3],
}

/// A rational point of the unit circle at about the angle `t` (its half
/// angle's tangent a rational, about the opposite direction near a half
/// turn).
fn unit_at(t: f64) -> [R; 2] {
    let flip = t.cos() < 0.0;
    let s = q((if flip { t - PI } else { t } / 2.0).tan());
    let den = int(1) + &s * &s;
    let (cs, sn) = ((int(1) - &s * &s) / &den, int(2) * &s / &den);
    if flip {
        [-cs, -sn]
    } else {
        [cs, sn]
    }
}

/// The world point of exact local coordinates.
fn world(f: &Affine, l: &[Qd; 3]) -> QV {
    qadd(
        &qadd(&qv(&f.o), &qscale(&f.x, &l[0])),
        &qadd(&qscale(&f.y, &l[1]), &qscale(&f.n, &l[2])),
    )
}

impl<'a> Quad<'a> {
    fn of(m: &'a Prism, fi: usize) -> Option<Self> {
        let (vm, vi) = m.view(fi);
        let g = |v: &V| [0, 1, 2].map(|k| rational_f64(&v[k]));
        let frame = |f: &Affine| (g(&f.o), [g(f.row(0)), g(f.row(1)), g(f.row(2))]);
        let ((o, rows), kind, data) = match &vm.faces[vi].surf {
            Surf::Plane { p, m } => (
                (g(p), [g(m), [0.0; 3], [0.0; 3]]),
                Kind::Plane { p, m },
                [0.0; 3],
            ),
            Surf::Cyl { c, r, .. } => (
                frame(&vm.f),
                Kind::Cyl { f: &vm.f, c, r },
                [rational_f64(&c[0]), rational_f64(&c[1]), rational_f64(r)],
            ),
            Surf::Cone { b, k } => (
                frame(&vm.f),
                Kind::Cone { f: &vm.f, b, k },
                [rational_f64(b), rational_f64(k), 0.0],
            ),
            _ => return None,
        };
        Some(Self {
            kind,
            o,
            rows,
            data,
        })
    }

    /// Local coordinates in binary64.
    fn local64(&self, x: [f64; 3]) -> [f64; 3] {
        let d = [0, 1, 2].map(|k| x[k] - self.o[k]);
        self.rows.map(|r| r[0] * d[0] + r[1] * d[1] + r[2] * d[2])
    }

    /// The surface's function at an exact point: off a plane along its
    /// normal, outside a cylinder's or a cone's circle positive.
    fn value(&self, x: &QV) -> Qd {
        match &self.kind {
            Kind::Plane { p, m } => qdot(&qsub(x, &qv(p)), m),
            Kind::Cyl { f, c, r } => {
                let l = f.local_q(x);
                let (du, dv) = (l[0].add_r(&-&c[0]), l[1].add_r(&-&c[1]));
                du.mul(&du).add(&dv.mul(&dv)).add_r(&-(*r * *r))
            }
            Kind::Cone { f, b, k } => {
                let l = f.local_q(x);
                let rad = l[2].scale(k).add_r(b);
                l[0].mul(&l[0]).add(&l[1].mul(&l[1])).sub(&rad.mul(&rad))
            }
        }
    }

    /// The function's gradient at an exact point.
    fn gradient(&self, x: &QV) -> QV {
        match &self.kind {
            Kind::Plane { m, .. } => qv(m),
            Kind::Cyl { f, c, .. } => {
                let l = f.local_q(x);
                let (du, dv) = (l[0].add_r(&-&c[0]), l[1].add_r(&-&c[1]));
                [0, 1, 2].map(|k| du.scale(&f.row(0)[k]).add(&dv.scale(&f.row(1)[k])))
            }
            Kind::Cone { f, b, k } => {
                let l = f.local_q(x);
                let rk = l[2].scale(k).add_r(b).scale(k);
                [0, 1, 2].map(|j| {
                    l[0].scale(&f.row(0)[j])
                        .add(&l[1].scale(&f.row(1)[j]))
                        .sub(&rk.scale(&f.row(2)[j]))
                })
            }
        }
    }

    /// The function and its gradient in binary64 (for searches).
    fn value64(&self, x: [f64; 3]) -> (f64, [f64; 3]) {
        let rows = &self.rows;
        match &self.kind {
            Kind::Plane { .. } => {
                let (p, m) = (self.o, rows[0]);
                let v = (0..3).map(|k| m[k] * (x[k] - p[k])).sum();
                (v, m)
            }
            Kind::Cyl { .. } => {
                let l = self.local64(x);
                let [cu, cv, r] = self.data;
                let (du, dv) = (l[0] - cu, l[1] - cv);
                let g = [0, 1, 2].map(|k| 2.0 * (du * rows[0][k] + dv * rows[1][k]));
                (du * du + dv * dv - r * r, g)
            }
            Kind::Cone { .. } => {
                let l = self.local64(x);
                let [b, k, _] = self.data;
                let rad = b + k * l[2];
                let g = [0, 1, 2]
                    .map(|j| 2.0 * (l[0] * rows[0][j] + l[1] * rows[1][j] - rad * k * rows[2][j]));
                (l[0] * l[0] + l[1] * l[1] - rad * rad, g)
            }
        }
    }

    /// A point of the surface near the foot of `x` on it: a plane's foot
    /// exactly; on a cylinder or a cone a rational point of it toward the
    /// foot's binary64 place (`None` off a cone's nappe).
    fn foot(&self, x: &QV) -> Option<QV> {
        match &self.kind {
            Kind::Plane { p, m } => {
                let k = qdot(&qsub(x, &qv(p)), m).scale(&(int(1) / dot(m, m)));
                Some(qsub(x, &qscale(m, &k)))
            }
            Kind::Cyl { f, c, r } => {
                let l = self.local64(qv_f64(x));
                let t = (l[1] - self.data[1]).atan2(l[0] - self.data[0]);
                let [cs, sn] = unit_at(t);
                let loc = [
                    Qd::rat(&c[0] + *r * &cs),
                    Qd::rat(&c[1] + *r * &sn),
                    Qd::rat(q(l[2])),
                ];
                Some(world(f, &loc))
            }
            Kind::Cone { f, b, k } => {
                let l = self.local64(qv_f64(x));
                let [bf, kf, _] = self.data;
                let rho = l[0].hypot(l[1]);
                let h = q((l[2] + kf * (rho - bf)) / (1.0 + kf * kf));
                let rad = *b + *k * &h;
                if sign(&rad) != Ordering::Greater {
                    return None;
                }
                let [cs, sn] = unit_at(l[1].atan2(l[0]));
                let loc = [Qd::rat(&rad * &cs), Qd::rat(&rad * &sn), Qd::rat(h)];
                Some(world(f, &loc))
            }
        }
    }
}

// ------------------------------------------------------- spheres

/// A sphere within the resolution of tangency to a cylinder or a cone face
/// of the other input (S9d.2, S9d.3), its nearest points in both faces, is
/// `Degenerate` crossing it and missing it where the gap lies outside
/// either input, as `graph::near_miss` takes a plane: a ball missing a rod's
/// wall by `1e-12` of its radius was fused with it into two solids, a ball
/// inside a rod cut from it behind a wall thinner than the resolution, and
/// one crossing it by that depth ran the validator's integrals for hours.
/// A gap inside both inputs (a block's cavity inside a rod by its wall)
/// evaluates. Each pair of faces whose boxes meet once one is widened by the
/// resolution is tried, before the meetings are found.
pub(super) fn sphere_quadrics(models: &[Prism; 2], res: f64) -> Result<()> {
    let wide = |b: &([f64; 3], [f64; 3])| (b.0.map(|x| x - res), b.1.map(|x| x + res));
    for s in 0..2 {
        let (sm, qm) = (&models[s], &models[1 - s]);
        for fs in 0..sm.faces.len() {
            let (vm, vi) = sm.view(fs);
            let Surf::Sphere { c, r } = &vm.faces[vi].surf else {
                continue;
            };
            for fq in 0..qm.faces.len() {
                let Some(quad) = Quad::of(qm, fq) else {
                    continue;
                };
                if matches!(quad.kind, Kind::Plane { .. })
                    || !boxes_meet(&wide(&sm.boxes[fs]), &qm.boxes[fq])
                {
                    continue;
                }
                sphere_quadric(sm, fs, c, r, qm, fq, &quad, res)?;
            }
        }
    }
    Ok(())
}

/// The cylinder's or cone's point nearest a point `c`, exact in the
/// quadratic field of its distance from the axis (in the local coordinates;
/// in a turned frame the nearest within rounding): on a cylinder its
/// circle's point toward it, on a cone the foot on its generatrix in the
/// plane through the axis (`None` off the cone's nappe); a point on the
/// axis taken toward the local `u`.
fn nearest_to(quad: &Quad, c: &V) -> Option<QV> {
    let (f, centre) = match &quad.kind {
        Kind::Cyl { f, c: cc, .. } => (*f, [cc[0].clone(), cc[1].clone()]),
        Kind::Cone { f, .. } => (*f, [zero(), zero()]),
        Kind::Plane { .. } => return None,
    };
    let l = f.local(c);
    let (du, dv) = (&l[0] - &centre[0], &l[1] - &centre[1]);
    let ww = &du * &du + &dv * &dv;
    // The unit direction from the axis and the distance from it.
    let (dir, rho) = if ww == zero() {
        ([Qd::rat(int(1)), Qd::rat(zero())], Qd::rat(zero()))
    } else {
        let inv = int(1) / &ww;
        (
            [
                Qd::new(zero(), &du * &inv, ww.clone()),
                Qd::new(zero(), &dv * &inv, ww.clone()),
            ],
            Qd::new(zero(), int(1), ww.clone()),
        )
    };
    let (rad, h) = match &quad.kind {
        Kind::Cyl { r, .. } => (Qd::rat((*r).clone()), Qd::rat(l[2].clone())),
        Kind::Cone { b, k, .. } => {
            // The foot of `(rho, h)` on `rho = b + k w`.
            let h = rho
                .add_r(&-(*b).clone())
                .scale(k)
                .add_r(&l[2])
                .scale(&(int(1) / (int(1) + *k * *k)));
            let rad = h.scale(k).add_r(b);
            if rad.sign() != Ordering::Greater {
                return None;
            }
            (rad, h)
        }
        Kind::Plane { .. } => unreachable!("taken above"),
    };
    let loc = [
        dir[0].mul(&rad).add_r(&centre[0]),
        dir[1].mul(&rad).add_r(&centre[1]),
        h,
    ];
    Some(world(f, &loc))
}

/// `sphere_quadrics`' rule for the sphere `(c, r)` on face `fs` of `sm`
/// against face `fq` of `qm`.
#[allow(clippy::too_many_arguments)]
fn sphere_quadric(
    sm: &Prism,
    fs: usize,
    c: &V,
    r: &R,
    qm: &Prism,
    fq: usize,
    quad: &Quad,
    res: f64,
) -> Result<()> {
    let Some(s) = nearest_to(quad, c) else {
        return Ok(());
    };
    let e = qsub(&s, &qv(c));
    let d2 = qqdot(&e, &e);
    let side = d2.add_r(&-(r * r)).sign();
    if side == Ordering::Equal || d2.sign() == Ordering::Equal {
        // Tangent exactly (the tangency's own rules), or about its centre.
        return Ok(());
    }
    let hi = r + q(res);
    if d2.add_r(&-(&hi * &hi)).sign() == Ordering::Greater {
        return Ok(());
    }
    let lo = r - q(res);
    if sign(&lo) == Ordering::Greater && d2.add_r(&-(&lo * &lo)).sign() == Ordering::Less {
        return Ok(());
    }
    // The sphere's point toward the surface's (a rational point of it to
    // rounding: its face's membership is decided there).
    let p = qv(&sphere_toward(c, r, qv_f64(&e)));
    if qm.in_face(fq, &s) == Loc::Out || sm.in_face(fs, &p) == Loc::Out {
        return Ok(());
    }
    if side == Ordering::Less {
        return Err(Error::Degenerate(NEAR_QUADRIC));
    }
    // The gap lies inside an input where its material holds its point
    // pushed across it.
    let gap = qsub(&p, &s);
    let inside = [
        qm.member(&s, std::slice::from_ref(&gap)) == Loc::In,
        sm.member(&p, &[gap.map(|x| x.neg())]) == Loc::In,
    ];
    if inside == [true, true] {
        return Ok(());
    }
    Err(Error::Degenerate(NEAR_QUADRIC))
}

// ------------------------------------------------------- edges

/// An input edge's line, conic or circle within the resolution of tangency
/// to a plane, cylinder or cone face of the other input (S9d.2), at a
/// point strictly inside the edge where the surface's function along the
/// edge is least or greatest and the surface's point nearest it in the
/// face, is `Degenerate` crossing the surface and missing it where the gap
/// lies outside either input, as a sphere is (`graph::edge_near_misses`):
/// a dome's rim circle missing a box's top within the resolution, or a
/// turned rod's cap circle, was fused with it into two solids, a rod's cap
/// circle on another rod's wall likewise. The extrema are exact where the
/// function is linear or quadratic along the curve (a conic or a circle
/// against a plane, its points in the quadratic field of the plane's
/// normal's reach along it; a line against a cylinder or a cone, a rational
/// parameter), and against a cylinder or a cone along a conic or a circle
/// a rational point of it at each extremum found in binary64 (the
/// derivative's root bisected to rounding). A given result's other curves
/// (meetings of two curved faces, cone, torus and spline curves) are taken
/// alike (`Run`): rational parameters at extrema found in binary64, settled
/// by bisection on the derivative's exact sign, a gap below `ROUNDED` an
/// incidence. A function constant along the edge (a circle parallel to a
/// plane, a line along a cylinder's axis, a circle about it) has no such
/// points: the faces' own rules hold those. Edges between faces on one
/// surface and seams' edges are none, as at a sphere; an edge one of whose
/// faces heads toward the surface, or runs along it, is no contact
/// (`leaves_away`).
pub(super) fn edge_faces(models: &[Prism; 2], res: f64) -> Result<()> {
    let wide = |b: &([f64; 3], [f64; 3])| (b.0.map(|x| x - res), b.1.map(|x| x + res));
    for o in 0..2 {
        let (em, qm) = (&models[o], &models[1 - o]);
        for (ei, e) in em.edges.iter().enumerate() {
            if e.id.is_none() || one_surface(em, e.faces[0], e.faces[1]) {
                continue;
            }
            let ebox = intersect(&em.boxes[e.faces[0]], &em.boxes[e.faces[1]]);
            let (from, to, with) = edge_places(em, ei);
            let run = Run::new(em, ei);
            for fq in 0..qm.faces.len() {
                let Some(quad) = Quad::of(qm, fq) else {
                    continue;
                };
                if !boxes_meet(&wide(&qm.boxes[fq]), &ebox) {
                    continue;
                }
                let found = match &run {
                    Some(run) => {
                        let slope = |x: &QV, t: &QV| qqdot(&quad.gradient(x), t).sign();
                        run.extrema(&|x| quad.value64(x).0, true)
                            .into_iter()
                            .filter_map(|(s, kind)| {
                                // Settled only within twice the resolution
                                // in binary64 (exact points cost).
                                let x = run.points(s, s, 1).first().copied()?;
                                let (v, g) = quad.value64(x);
                                let gg = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                                if v.is_nan() || v.abs() > 2.0 * res * gg {
                                    return None;
                                }
                                Some((run.settle(s, &slope)?, kind))
                            })
                            .collect()
                    }
                    None => extrema(&e.curve, &quad),
                };
                for (a, extremum) in found {
                    if strictly_within(&place(&e.curve, &a), &from, &to, with) != Some(true) {
                        continue;
                    }
                    let near = Near {
                        extremum,
                        rounded: run.is_some(),
                    };
                    edge_face(em, ei, qm, fq, &quad, &a, near, res)?;
                }
            }
        }
    }
    Ok(())
}

/// A point of an edge where the surface's function along it is least
/// (`Less`) or greatest (`Greater`), found in binary64 on another curve
/// (`rounded`: a gap below `ROUNDED` there stands for an incidence).
struct Near {
    extremum: Ordering,
    rounded: bool,
}

/// `edge_faces`' rule at a point `a` of edge `ei` of `em`.
#[allow(clippy::too_many_arguments)]
fn edge_face(
    em: &Prism,
    ei: usize,
    qm: &Prism,
    fq: usize,
    quad: &Quad,
    a: &QV,
    near: Near,
    res: f64,
) -> Result<()> {
    let extremum = near.extremum;
    let side = quad.value(a).sign();
    if side == Ordering::Equal {
        // On the surface: the incidences' own rules.
        return Ok(());
    }
    let Some(s) = quad.foot(a) else {
        return Ok(());
    };
    let gap = qsub(&s, a);
    let tol = q(res);
    let gg = qqdot(&gap, &gap);
    if gg.add_r(&-(&tol * &tol)).sign() == Ordering::Greater {
        return Ok(());
    }
    let tiny = q(ROUNDED);
    if near.rounded && gg.add_r(&-(&tiny * &tiny)).sign() != Ordering::Greater {
        return Ok(());
    }
    if qm.in_face(fq, &s) == Loc::Out {
        return Ok(());
    }
    // The inputs touch there only where both of the edge's faces leave it
    // away from the surface (the function growing into them from a least
    // value, falling from a greatest): a face heading toward the surface
    // crosses it there, the inputs overlapping beyond the resolution, and
    // its section near the edge is the arrangement's own (S9f.2b.2's
    // `lens_tilt_loop`: a rod's top rim on a lens's top plane but for the
    // turned frame's rounding, its wall through the plane). A face running
    // along the surface there (a rod's wall on a box's face on its tangent
    // plane, the cap circle within rounding of the plane) is a tangency or
    // an incidence of the two faces, the faces' own rules' (S9c.1's turned
    // box on a rod's tangent plane evaluating as before).
    if !leaves_away(em, ei, quad, a, extremum) {
        return Ok(());
    }
    // A least value below the surface, or a greatest above it: the edge
    // crosses it.
    if side == extremum {
        return Err(Error::Degenerate(NEAR_FACE));
    }
    let inside = [
        em.member(a, std::slice::from_ref(&gap)) == Loc::In,
        qm.member(&s, &[gap.map(|x| x.neg())]) == Loc::In,
    ];
    if inside == [true, true] {
        return Ok(());
    }
    Err(Error::Degenerate(NEAR_FACE))
}

/// Whether both faces of edge `ei` of `em` leave its point `a` away from
/// the surface: the function's gradient against each face's direction into
/// it from the edge (its outward normal across the edge's running tangent,
/// the face on the left `n x t`, on the right `t x n`) positive at a least
/// value, negative at a greatest. A face whose direction lies in the
/// surface's tangent plane there, exactly or within rounding (`along`),
/// runs along the surface: no contact of the edge's own.
fn leaves_away(em: &Prism, ei: usize, quad: &Quad, a: &QV, extremum: Ordering) -> bool {
    let e = &em.edges[ei];
    let pos = place(&e.curve, a);
    let mut t = super::meet::tangent(&e.curve, &pos, a);
    // The running direction: a line's from its start's key to its end's,
    // other curves' with their parameter or against it.
    let (from, to, with) = edge_places(em, ei);
    let forward = match (&e.curve, &from, &to) {
        (Crv::Line { .. }, Pos::T(s0), Pos::T(s1)) => s0.cmp(s1) == Ordering::Less,
        _ => with,
    };
    if !forward {
        t = t.map(|x| x.neg());
    }
    let g = quad.gradient(a);
    let away = [0, 1].map(|k| {
        let n = em.normal_at(e.faces[k], a);
        let inward = if k == 0 {
            qcross(&n, &t)
        } else {
            qcross(&t, &n)
        };
        if along(&g, &inward) {
            return Ordering::Equal;
        }
        qqdot(&g, &inward).sign()
    });
    away.iter().all(|s| *s != extremum && *s != Ordering::Equal)
}

/// Whether a direction `w` lies in the tangent plane of a surface whose
/// gradient is `g`, exactly or within rounding: the sine of its angle to
/// the plane at most `10^-12` (`(g . w)^2 10^24 <= |g|^2 |w|^2`, decided
/// exactly), the band in which the faces' own rules take a plane or an
/// axis along a cylinder's direction (`meet::within_rounding_of_parallel`,
/// a plane within rounding of a cylinder's direction): a direction whose
/// sign there an ulp of a turned frame decides is theirs.
fn along(g: &QV, w: &QV) -> bool {
    let d = qqdot(g, w);
    let lhs = d.mul(&d).scale(&int(10).pow(24));
    lhs.sub(&qqdot(g, g).mul(&qqdot(w, w))).sign() != Ordering::Greater
}

/// The points of a curve where the surface's function along it is least
/// (`Less`) or greatest (`Greater`), exact on the curve.
fn extrema(crv: &Crv, quad: &Quad) -> Vec<(QV, Ordering)> {
    match (crv, &quad.kind) {
        (Crv::Line { .. }, Kind::Plane { .. }) => Vec::new(),
        (Crv::Line { p, d }, Kind::Cyl { f, .. } | Kind::Cone { f, .. }) => {
            // The function is quadratic in the parameter: `A t^2 + B t + C`.
            let (lp, ld) = (f.local_q(p), f.local_dir(d));
            let (a, b) = match &quad.kind {
                Kind::Cyl { c, .. } => {
                    let (u, v) = (lp[0].add_r(&-&c[0]), lp[1].add_r(&-&c[1]));
                    (
                        &ld[0] * &ld[0] + &ld[1] * &ld[1],
                        u.scale(&ld[0]).add(&v.scale(&ld[1])).scale(&int(2)),
                    )
                }
                Kind::Cone { b: cb, k, .. } => {
                    let rad = lp[2].scale(k).add_r(cb);
                    (
                        &ld[0] * &ld[0] + &ld[1] * &ld[1] - *k * *k * &ld[2] * &ld[2],
                        lp[0]
                            .scale(&ld[0])
                            .add(&lp[1].scale(&ld[1]))
                            .sub(&rad.scale(&(*k * &ld[2])))
                            .scale(&int(2)),
                    )
                }
                Kind::Plane { .. } => unreachable!("taken above"),
            };
            if a == zero() {
                return Vec::new();
            }
            let t = b.scale(&(int(-1) / (int(2) * &a)));
            let kind = if a > zero() {
                Ordering::Less
            } else {
                Ordering::Greater
            };
            vec![(qadd(p, &qscale(d, &t)), kind)]
        }
        (Crv::Conic { c, a, b }, Kind::Plane { m, .. }) => {
            // `c + (a (m.a) + b (m.b)) / sqrt(K)`, `K = (m.a)^2 + (m.b)^2`,
            // greatest; less it, least.
            let (ma, mb) = (dot(m, a), dot(m, b));
            let kk = &ma * &ma + &mb * &mb;
            if kk == zero() {
                return Vec::new();
            }
            let dir = add(&scale(a, &ma), &scale(b, &mb));
            [(int(1), Ordering::Greater), (int(-1), Ordering::Less)]
                .into_iter()
                .map(|(s, kind)| {
                    let k = &s / &kk;
                    let x = [0, 1, 2].map(|i| Qd::new(c[i].clone(), &k * &dir[i], kk.clone()));
                    (x, kind)
                })
                .collect()
        }
        (Crv::Circle(circ), Kind::Plane { m, .. }) => {
            // The circle's point along the normal's part in its plane.
            let (xx, yy) = (dot(&circ.x, &circ.x), dot(&circ.y, &circ.y));
            let (dx, dy) = (dot(m, &circ.x) / xx, dot(m, &circ.y) / yy);
            if dx == zero() && dy == zero() {
                return Vec::new();
            }
            vec![
                (circ.at(&[dx.clone(), dy.clone()]), Ordering::Greater),
                (circ.at(&[-dx, -dy]), Ordering::Less),
            ]
        }
        (Crv::Conic { c, a, b }, _) => {
            let g = |v: &V| v.clone().map(|x| rational_f64(&x));
            let (fc, fa, fb) = (g(c), g(a), g(b));
            let at = |t: f64| {
                let (cs, sn) = (t.cos(), t.sin());
                (
                    [0, 1, 2].map(|k| fc[k] + fa[k] * cs + fb[k] * sn),
                    [0, 1, 2].map(|k| -fa[k] * sn + fb[k] * cs),
                )
            };
            search(&at, quad)
                .into_iter()
                .map(|(t, kind)| {
                    let [cs, sn] = unit_at(t);
                    let x = super::meet::conic_point(c, a, b, &[Qd::rat(cs), Qd::rat(sn)]);
                    (x, kind)
                })
                .collect()
        }
        (Crv::Circle(circ), _) => {
            let g = |v: &V| v.clone().map(|x| rational_f64(&x));
            let (fc, fx, fy) = (g(&circ.c), g(&circ.x), g(&circ.y));
            let (xx, yy) = (dot(&circ.x, &circ.x), dot(&circ.y, &circ.y));
            let (lx, ly) = (rational_f64(&xx).sqrt(), rational_f64(&yy).sqrt());
            let rad = rational_f64(&circ.r2).sqrt();
            let at = |t: f64| {
                let (cs, sn) = (t.cos() * rad / lx, t.sin() * rad / ly);
                let (dc, ds) = (-t.sin() * rad / lx, t.cos() * rad / ly);
                (
                    [0, 1, 2].map(|k| fc[k] + fx[k] * cs + fy[k] * sn),
                    [0, 1, 2].map(|k| fx[k] * dc + fy[k] * ds),
                )
            };
            search(&at, quad)
                .into_iter()
                .map(|(t, kind)| {
                    let x = circ.at(&[q(t.cos() / lx), q(t.sin() / ly)]);
                    (x, kind)
                })
                .collect()
        }
        _ => Vec::new(),
    }
}

/// The angles of a closed curve `at(t) = (point, velocity)` where the
/// surface's function along it is least or greatest: each change of its
/// derivative's sign among 64 samples bisected to rounding.
fn search(at: &dyn Fn(f64) -> ([f64; 3], [f64; 3]), quad: &Quad) -> Vec<(f64, Ordering)> {
    let slope = |t: f64| {
        let (x, d) = at(t);
        let (_, g) = quad.value64(x);
        g[0] * d[0] + g[1] * d[1] + g[2] * d[2]
    };
    const N: usize = 64;
    let step = TAU / N as f64;
    let ts: Vec<f64> = (0..N).map(|i| i as f64 * step).collect();
    let ds: Vec<f64> = ts.iter().map(|&t| slope(t)).collect();
    // A function constant along the curve to rounding has none.
    let size = ds.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    if size.is_nan() || size <= 0.0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for i in 0..N {
        let (d0, d1) = (ds[i], ds[(i + 1) % N]);
        let kind = if d0 < 0.0 && d1 >= 0.0 {
            Ordering::Less
        } else if d0 > 0.0 && d1 <= 0.0 {
            Ordering::Greater
        } else {
            continue;
        };
        let (mut lo, mut hi) = (ts[i], ts[i] + step);
        for _ in 0..80 {
            let mid = 0.5 * lo + 0.5 * hi;
            if !(lo < mid && mid < hi) {
                break;
            }
            let s = slope(mid);
            if (s < 0.0) == (kind == Ordering::Less) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        out.push((0.5 * lo + 0.5 * hi, kind));
    }
    out
}

// ------------------------------------------------------- other curves

/// An edge's curve other than a line, a conic or a circle (a meeting of
/// two curved faces, a cone's, a torus's section or meeting: a given
/// result's edges; a spline curve: a spline prism's cap edges too) over
/// its places in binary64, by a
/// fraction `s` of its run: its parameter `t0 + span s` (a height or a run
/// parameter, or an angle turning `span`).
struct Run<'a> {
    crv: &'a Crv,
    t0: f64,
    span: f64,
}

impl<'a> Run<'a> {
    fn new(m: &'a Prism, ei: usize) -> Option<Self> {
        let crv = &m.edges[ei].curve;
        let (from, to, with) = edge_places(m, ei);
        match (crv, &from, &to) {
            (Crv::Rise(_) | Crv::WallMeet(_) | Crv::Spline(_), Pos::T(a), Pos::T(b)) => {
                let (t0, t1) = (a.to_f64(), b.to_f64());
                Some(Self {
                    crv,
                    t0,
                    span: t1 - t0,
                })
            }
            (
                Crv::Meet(_) | Crv::Cone(_) | Crv::Torus(_) | Crv::Toric(_),
                Pos::Ang(a),
                Pos::Ang(b),
            ) => {
                let angle = |p: &[Qd; 2]| p[1].to_f64().atan2(p[0].to_f64());
                let (t0, t1) = (angle(a), angle(b));
                let mut sweep = if with { t1 - t0 } else { t0 - t1 };
                sweep = sweep.rem_euclid(TAU);
                if sweep == 0.0 {
                    sweep = TAU;
                }
                Some(Self {
                    crv,
                    t0,
                    span: if with { sweep } else { -sweep },
                })
            }
            _ => None,
        }
    }

    /// `n + 1` points from fraction `s0` to `s1`.
    fn points(&self, s0: f64, s1: f64, n: usize) -> Vec<[f64; 3]> {
        let (a, b) = (self.t0 + self.span * s0, self.t0 + self.span * s1);
        match self.crv {
            Crv::Rise(c) => c.samples(a, b, n),
            Crv::WallMeet(c) => c.samples(a, b, n),
            Crv::Spline(c) => c.samples(a, b, n),
            Crv::Meet(c) => c.samples(a, b - a, n),
            Crv::Cone(c) => c.samples(a, b - a, n),
            Crv::Torus(c) => c.samples(a, b - a, n),
            Crv::Toric(c) => c.samples(a, b - a, n),
            _ => Vec::new(),
        }
    }

    /// The curve's exact point at a rational parameter near the fraction's
    /// (`None` off its piece).
    fn exact(&self, s: f64) -> Option<QV> {
        self.exact_at(self.t0 + self.span * s)
    }

    /// The curve's exact point at the parameter `t` rounded to a rational.
    fn exact_at(&self, t: f64) -> Option<QV> {
        match self.crv {
            Crv::Rise(c) => c.at(&q(t)).ok().flatten(),
            Crv::WallMeet(c) => c.at(&q(t)),
            Crv::Spline(c) => Some(c.point(&Qd::rat(q(t)))),
            Crv::Meet(c) => c.at(&unit_at(t)),
            Crv::Cone(c) => c.at(&unit_at(t)),
            Crv::Torus(c) => c.at(&unit_at(t)),
            Crv::Toric(c) => c.at(&unit_at(t)),
            _ => None,
        }
    }

    /// The exact point at an extremum found near the fraction `s`, its
    /// parameter narrowed to rounding by bisection on the exact sign of the
    /// function's derivative along the curve, `slope(x, t)` at a point `x`
    /// with its tangent `t` (the binary64 search's place, flat there, is
    /// only within about `1e-8` of it, the function's value off its
    /// extremum by about `1e-16` of its size: an incidence's point not told
    /// from a near miss).
    fn settle(&self, s: f64, slope: &dyn Fn(&QV, &QV) -> Ordering) -> Option<QV> {
        let sign_at = |t: f64| {
            let x = self.exact_at(t)?;
            let tangent = super::meet::tangent(self.crv, &place(self.crv, &x), &x);
            Some(slope(&x, &tangent))
        };
        let h = 1e-6;
        let (mut a, mut b) = (
            self.t0 + self.span * (s - h).max(0.0),
            self.t0 + self.span * (s + h).min(1.0),
        );
        let (sa, sb) = (sign_at(a), sign_at(b));
        let opposite = matches!(
            (sa, sb),
            (Some(Ordering::Less), Some(Ordering::Greater))
                | (Some(Ordering::Greater), Some(Ordering::Less))
        );
        if !opposite {
            return self.exact(s);
        }
        for _ in 0..80 {
            let m = 0.5 * a + 0.5 * b;
            if m == a || m == b {
                break;
            }
            match sign_at(m) {
                Some(Ordering::Equal) => return self.exact_at(m),
                Some(x) if Some(x) == sa => a = m,
                Some(_) => b = m,
                None => break,
            }
        }
        self.exact_at(0.5 * a + 0.5 * b)
    }

    /// The fractions strictly inside the run where `g` is least (`Less`)
    /// or, with `maxima`, greatest (`Greater`): each among 32 samples,
    /// narrowed by eight samples about the best at a time to rounding.
    fn extrema(&self, g: &dyn Fn([f64; 3]) -> f64, maxima: bool) -> Vec<(f64, Ordering)> {
        const N: usize = 32;
        const M: usize = 8;
        let vals: Vec<f64> = self.points(0.0, 1.0, N).into_iter().map(g).collect();
        if vals.len() != N + 1 || vals.iter().any(|v| !v.is_finite()) {
            return Vec::new();
        }
        let mut out = Vec::new();
        for i in 1..N {
            let kind = if vals[i] <= vals[i - 1] && vals[i] < vals[i + 1] {
                Ordering::Less
            } else if maxima && vals[i] >= vals[i - 1] && vals[i] > vals[i + 1] {
                Ordering::Greater
            } else {
                continue;
            };
            let better = |a: f64, b: f64| match kind {
                Ordering::Less => a < b,
                _ => a > b,
            };
            let step = 1.0 / N as f64;
            let (mut lo, mut hi) = ((i - 1) as f64 * step, (i + 1) as f64 * step);
            for _ in 0..64 {
                let vs: Vec<f64> = self.points(lo, hi, M).into_iter().map(g).collect();
                if vs.len() != M + 1 {
                    break;
                }
                let mut j = 0;
                for (k, v) in vs.iter().enumerate() {
                    if better(*v, vs[j]) {
                        j = k;
                    }
                }
                let w = (hi - lo) / M as f64;
                let (nlo, nhi) = (
                    lo + w * j.saturating_sub(1) as f64,
                    lo + w * (j + 1).min(M) as f64,
                );
                if nhi - nlo >= hi - lo {
                    break;
                }
                (lo, hi) = (nlo, nhi);
            }
            out.push((0.5 * lo + 0.5 * hi, kind));
        }
        out
    }
}

/// The points of an edge's curve other than a line, a conic or a circle
/// nearest the centre `c` of a sphere of radius `r` (S9d.1's rule at a
/// sphere, `graph::edge_nearest`): rational parameters at its distance's
/// least values found in binary64, its distance there exceeding the least
/// by far less than the resolution. A point whose distance from the sphere
/// is below `ROUNDED` is no near miss: such a point stands for one the
/// sphere passes through exactly but for its own rounding (S9e.3b's
/// `peg_kiss`, a ball through a given meeting's point tangent to it), the
/// incidences' own rules.
pub(super) fn nearest_on_run(m: &Prism, ei: usize, c: &V, r: &R, res: f64) -> Vec<QV> {
    let Some(run) = Run::new(m, ei) else {
        return Vec::new();
    };
    let fc = [0, 1, 2].map(|k| rational_f64(&c[k]));
    let g = |x: [f64; 3]| (0..3).map(|k| (x[k] - fc[k]).powi(2)).sum::<f64>();
    // `|dd - r^2| = |d - r| (d + r)`, about `2 r` times the gap.
    let bound = r * q(2.0 * ROUNDED);
    let rf = rational_f64(r);
    let slope = |x: &QV, t: &QV| qqdot(&qsub(x, &qv(c)), t).sign();
    run.extrema(&g, false)
        .into_iter()
        .filter_map(|(s, _)| {
            // Settled only within twice the resolution in binary64.
            let x = run.points(s, s, 1).first().copied()?;
            let off = (g(x).sqrt() - rf).abs();
            if off.is_nan() || off > 2.0 * res {
                return None;
            }
            run.settle(s, &slope)
        })
        .filter(|a| {
            let w = qsub(a, &qv(c));
            let off = qqdot(&w, &w).add_r(&-(r * r));
            off.mul(&off).add_r(&-(&bound * &bound)).sign() == Ordering::Greater
        })
        .collect()
}

/// A gap below which a point found in binary64 on another curve stands for
/// an incidence (far below its rounding's effect on the distance, the
/// square of binary64's precision): `2^-80`.
const ROUNDED: f64 = 1.0 / (1u128 << 80) as f64;
