//! S9d's near misses at cylinders, cones and planes (REVIEW_NOTES.md, "A
//! ball within the resolution of a cylinder or cone, and edges within it of
//! a face"): `graph.rs`'s `near_miss` holds a sphere against a plane or a
//! sphere, and `edge_near_misses` an edge or a vertex against a sphere;
//! here a sphere against a cylinder or a cone face (`sphere_quadrics`), an
//! edge's line, conic, circle or other curve against a plane, cylinder or
//! cone face (`edge_faces`), a vertex against such a face (`vertex_faces`,
//! "A vertex within the resolution of a face"; a cone's apex too), the
//! inputs' vertices and edges against each other's (`edges_and_vertices`,
//! "Edges and vertices within the resolution of each other"; two conic or
//! circle edges, `conics_closest`, "Two conic edges within the resolution
//! of each other"), and the
//! points of other curves nearest a sphere for `edge_near_misses`
//! (`nearest_on_run`), by the same rule:
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
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};
use std::sync::Arc;

/// A sphere within the resolution of tangency to a cylinder or cone face.
pub(super) const NEAR_QUADRIC: &str =
    "a sphere within the resolution of tangency to a cylinder or cone (S9d.2)";
/// An edge within the resolution of tangency to a face of the other input.
pub(super) const NEAR_FACE: &str = "an edge within the resolution of tangency to a face (S9d.2)";
/// A vertex within the resolution of a face of the other input.
pub(super) const NEAR_VERTEX_FACE: &str = "a vertex within the resolution of a face (S9d.2)";
/// A cone's apex within the resolution of a face of the other input.
pub(super) const NEAR_APEX: &str = "a cone's apex within the resolution of a face (S9d.3)";
/// A vertex within the resolution of a vertex of the other input.
pub(super) const NEAR_VERTICES: &str =
    "a vertex within the resolution of a vertex of the other input (S9d.2)";
/// A vertex within the resolution of an edge of the other input.
pub(super) const NEAR_VERTEX_EDGE: &str =
    "a vertex within the resolution of an edge of the other input (S9d.2)";
/// Two edges of the inputs within the resolution of each other.
pub(super) const NEAR_EDGES: &str = "two edges within the resolution of each other (S9d.2)";
/// A near miss whose least distance, bracketed, cannot be told from the
/// resolution itself.
pub(super) const NEAR_LIMIT: &str =
    "a near miss's least distance within rounding of the resolution";

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

    /// Whether a point `x` (binary64) lies off a plane or a cylinder by more
    /// than the resolution beyond doubt: its function over its gradient's
    /// length above `4 res`. That ratio is a plane's distance, and within
    /// `7 / 6` of a cylinder's distance `d <= res` where its radius exceeds
    /// `4 res` (`|rho^2 - r^2| / 2 rho = d (rho + r) / 2 rho`, `rho >= 3r /
    /// 4`), far beyond binary64's error. A cone's ratio near its apex is not
    /// so bounded: none.
    fn clearly_off(&self, x: [f64; 3], res: f64) -> bool {
        let (v, g) = self.value64(x);
        let gg = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        match self.kind {
            Kind::Plane { .. } => v.abs() > 4.0 * res * gg,
            Kind::Cyl { .. } => self.data[2] > 4.0 * res && v.abs() > 4.0 * res * gg,
            Kind::Cone { .. } => false,
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
                    // A plane's or a cylinder's exact extremum binary64
                    // puts well beyond the resolution is none (`edge_face`
                    // takes only gaps within it; `Run`'s are settled so).
                    if run.is_none() && quad.clearly_off(qv_f64(&a), res) {
                        continue;
                    }
                    // Strictly inside the edge, or at an end it runs on
                    // through smoothly (a closed edge's vertex, a circle's
                    // seam).
                    match strictly_within(&place(&e.curve, &a), &from, &to, with) {
                        Some(true) => {}
                        None if smooth_end(em, ei, &a).is_some() => {}
                        _ => continue,
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
    let g = quad.gradient(a);
    let away = edge_dirs(em, ei, a).map(|inward| {
        if along(&g, &inward) {
            return Ordering::Equal;
        }
        qqdot(&g, &inward).sign()
    });
    away.iter().all(|s| *s != extremum && *s != Ordering::Equal)
}

/// The directions into the faces of edge `ei` of `em` from its point `a`:
/// each face's outward normal across the edge's running tangent (the face
/// on the left `n x t`, on the right `t x n`).
fn edge_dirs(em: &Prism, ei: usize, a: &QV) -> [QV; 2] {
    let e = &em.edges[ei];
    let t = running(em, ei, a);
    [0, 1].map(|k| {
        let n = em.normal_at(e.faces[k], a);
        if k == 0 {
            qcross(&n, &t)
        } else {
            qcross(&t, &n)
        }
    })
}

/// Whether a point at an end of edge `ei` of `m` lies inside its input
/// edge, the edge's curve running on through its vertex smoothly: the
/// vertex holds exactly two edge ends, of one closed edge or of two model
/// edges of one input edge (a circle's halves at their seam), but for
/// seams' edges (without an input's id, between faces on one surface: a
/// circle's wall's seam). The vertex, where it does.
pub(super) fn smooth_end(m: &Prism, ei: usize, x: &QV) -> Option<usize> {
    let e = &m.edges[ei];
    let id = e.id?;
    let v = [e.start, e.end]
        .into_iter()
        .find(|&v| qv_eq(&m.verts[v].p, x))?;
    let mut ends = 0;
    for f in &m.edges {
        let here = usize::from(f.start == v) + usize::from(f.end == v);
        if here == 0 || (f.id.is_none() && one_surface(m, f.faces[0], f.faces[1])) {
            continue;
        }
        if f.id != Some(id) {
            return None;
        }
        ends += here;
    }
    (ends == 2).then_some(v)
}

/// Whether a direction `w` lies along the plane normal to `g` within the
/// faces' band and a further allowance `rho2` (the square of a length):
/// `(g . w)^2 <= (10^-24 |g|^2 + rho2) |w|^2`, exactly. With `rho2` zero it
/// is `along`.
fn along_by(g: &QV, w: &QV, rho2: &R) -> bool {
    let d = qqdot(g, w);
    let band = qqdot(g, g).scale(&(int(1) / int(10).pow(24))).add_r(rho2);
    d.mul(&d).sub(&band.mul(&qqdot(w, w))).sign() != Ordering::Greater
}

/// The sign every direction takes against `g`, strictly and outside the
/// band (`along`, or `along_by` with an allowance); `None` where one lies
/// within it, they take both signs, or there are none.
fn signs(dirs: &[QV], g: &QV, rho2: Option<&R>) -> Option<Ordering> {
    let mut sign = None;
    for w in dirs {
        let flat = match rho2 {
            None => along(g, w),
            Some(r) => along_by(g, w, r),
        };
        if flat {
            return None;
        }
        let s = qqdot(g, w).sign();
        if sign.is_some_and(|x| x != s) {
            return None;
        }
        sign = Some(s);
    }
    sign
}

/// The sign every ruling of a cone's wall from its apex (`|k| (cos u X + sin
/// u Y) + s N` on the frame's axes, `s` the side the cone opens to, the
/// sign of `k`) takes against `g`, strictly and outside the band: with
/// `(A, B, C)` the components of `g` along `X`, `Y` and `N`, every ruling's
/// is `|k| (A cos u + B sin u) + s C`, all of the sign of `s C` where `C^2 >
/// k^2 (A^2 + B^2)`; the ruling nearest the plane normal to `g` is within
/// the band where `|C| - |k| sqrt(A^2 + B^2) <= E`, `E^2 = (10^-24 |g|^2 +
/// rho2) (1 + k^2)` (the rulings' lengths in the frame's coordinates),
/// decided by squaring. `None` for a model without an apex.
fn apex_side(m: &Prism, g: &QV, rho2: &R) -> Option<Ordering> {
    let fun = m.funnel.as_ref()?;
    fun.apex?;
    let f = &m.f;
    let (ca, cb, cn) = (qdot(g, &f.x), qdot(g, &f.y), qdot(g, &f.n));
    let k = &fun.k;
    let c = if sign(k) == Ordering::Greater {
        cn
    } else {
        cn.neg()
    };
    let kk = k * k;
    let p2 = ca.mul(&ca).add(&cb.mul(&cb)).scale(&kk);
    let c2 = c.mul(&c);
    if c.sign() == Ordering::Equal || c2.sub(&p2).sign() != Ordering::Greater {
        return None;
    }
    let e2 = qqdot(g, g)
        .scale(&(int(1) / int(10).pow(24)))
        .add_r(rho2)
        .scale(&(int(1) + &kk));
    // `|C| <= E + P` (`P^2 = p2`): `C^2 - E^2 - P^2 <= 2 E P`.
    let lhs = c2.sub(&e2).sub(&p2);
    if lhs.sign() != Ordering::Greater
        || lhs.mul(&lhs).sub(&e2.mul(&p2).scale(&int(4))).sign() != Ordering::Greater
    {
        return None;
    }
    Some(c.sign())
}

/// An edge's tangent at its point `a` along its run from its start to its
/// end: a line's from its start's key to its end's, other curves' with
/// their parameter or against it.
fn running(em: &Prism, ei: usize, a: &QV) -> QV {
    let e = &em.edges[ei];
    let t = super::meet::tangent(&e.curve, &place(&e.curve, a), a);
    let (from, to, with) = edge_places(em, ei);
    let forward = match (&e.curve, &from, &to) {
        (Crv::Line { .. }, Pos::T(s0), Pos::T(s1)) => s0.cmp(s1) == Ordering::Less,
        _ => with,
    };
    if forward {
        t
    } else {
        t.map(|x| x.neg())
    }
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

// ------------------------------------------------------- vertices

/// An input vertex within the resolution of a plane, cylinder or cone face
/// of the other input (S9d.2, REVIEW_NOTES.md's "A vertex within the
/// resolution of a face"), the surface's point near its foot in the face
/// and its input's boundary leaving it strictly on one side of the
/// surface's level through it (`corner`), is `Degenerate` crossing the
/// surface and missing it where the gap lies outside either input, as an
/// edge is (`edge_faces`): a cube's corner `1e-12` off a rod's wall was
/// fused with it into two solids, a cube inside a rod cut from it as a
/// cavity behind a wall thinner than the resolution. A vertex on the
/// surface exactly is the incidences'; one whose boundary runs along the
/// surface (within the faces' parallel band) or crosses its level there is
/// no contact of its own. Seams' and poles' vertices (without an input's
/// id) and those among faces of one surface are none, as at a sphere. A
/// cone's apex, which has no edges, is taken with its rulings for its
/// tangent cone (`apex_side`; REVIEW_NOTES.md's "Edges and vertices within
/// the resolution of each other"): a cone standing `1e-12` under a slab was
/// fused with it into two solids.
pub(super) fn vertex_faces(models: &[Prism; 2], res: f64) -> Result<()> {
    let wide = |b: &([f64; 3], [f64; 3])| (b.0.map(|x| x - res), b.1.map(|x| x + res));
    for o in 0..2 {
        let (vm, qm) = (&models[o], &models[1 - o]);
        let ends = edge_ends(vm);
        for (v, mv) in vm.verts.iter().enumerate() {
            let Some(spot) = spot_of(vm, &ends, v) else {
                continue;
            };
            let x = qv_f64(&mv.p);
            for fq in 0..qm.faces.len() {
                let b = wide(&qm.boxes[fq]);
                if !(0..3).all(|k| b.0[k] <= x[k] && x[k] <= b.1[k]) {
                    continue;
                }
                let Some(quad) = Quad::of(qm, fq) else {
                    continue;
                };
                vertex_face(vm, &spot, &mv.p, qm, fq, &quad, res)?;
            }
        }
    }
    Ok(())
}

/// Each vertex's edge ends: the edge and whether it starts there.
fn edge_ends(m: &Prism) -> Vec<Vec<(usize, bool)>> {
    let mut ends: Vec<Vec<(usize, bool)>> = vec![Vec::new(); m.verts.len()];
    for (ei, e) in m.edges.iter().enumerate() {
        ends[e.start].push((ei, true));
        ends[e.end].push((ei, false));
    }
    ends
}

/// A point of an input's boundary whose directions there are decided
/// against a gradient or a gap: a vertex (its edge ends), a cone's apex
/// (its rulings), or a point of an edge.
enum Spot<'a> {
    Vertex(&'a [(usize, bool)]),
    Apex,
    Edge(usize),
}

/// Vertex `v` of `m` as a spot the near-miss rules take: an input vertex
/// (with an input's id) whose edges' faces lie on more than one surface,
/// or a cone's apex.
fn spot_of<'a>(m: &Prism, ends: &'a [Vec<(usize, bool)>], v: usize) -> Option<Spot<'a>> {
    m.verts[v].id?;
    if m.funnel.as_ref().and_then(|f| f.apex) == Some(v) {
        return Some(Spot::Apex);
    }
    let faces: Vec<usize> = ends[v]
        .iter()
        .flat_map(|&(ei, _)| m.edges[ei].faces)
        .collect();
    if faces.is_empty() || faces.iter().all(|&f| one_surface(m, faces[0], f)) {
        return None;
    }
    Some(Spot::Vertex(&ends[v]))
}

/// The sign every direction of `m`'s boundary at the spot `a` takes against
/// `g`, strictly and outside the band (`signs`, `apex_side`): `None` where
/// one lies within it, they take both signs, or the spot is undecided (a
/// reflex or half-turn sector, a face with other than two edge ends).
fn side(m: &Prism, spot: &Spot, a: &QV, g: &QV, rho2: Option<&R>) -> Option<Ordering> {
    match spot {
        Spot::Apex => apex_side(m, g, rho2.unwrap_or(&zero())),
        Spot::Vertex(ends) => signs(&corner_dirs(m, ends, a)?, g, rho2),
        Spot::Edge(ei) => signs(&edge_dirs(m, *ei, a), g, rho2),
    }
}

/// `vertex_faces`' rule at a vertex `a` of `vm` (a spot of it).
fn vertex_face(
    vm: &Prism,
    spot: &Spot,
    a: &QV,
    qm: &Prism,
    fq: usize,
    quad: &Quad,
    res: f64,
) -> Result<()> {
    let why = match spot {
        Spot::Apex => NEAR_APEX,
        _ => NEAR_VERTEX_FACE,
    };
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
    if qqdot(&gap, &gap).add_r(&-(&tol * &tol)).sign() == Ordering::Greater {
        return Ok(());
    }
    if qm.in_face(fq, &s) == Loc::Out {
        return Ok(());
    }
    // Rising along every direction: least there.
    let rises = match spot {
        Spot::Vertex(ends) => corner(vm, ends, quad, a),
        _ => self::side(vm, spot, a, &quad.gradient(a), None),
    };
    let Some(extremum) = rises.map(Ordering::reverse) else {
        return Ok(());
    };
    // A least value below the surface, or a greatest above it: the corner
    // crosses it.
    if side == extremum {
        return Err(Error::Degenerate(why));
    }
    let inside = [
        vm.member(a, std::slice::from_ref(&gap)) == Loc::In,
        qm.member(&s, &[gap.map(|x| x.neg())]) == Loc::In,
    ];
    if inside == [true, true] {
        return Ok(());
    }
    Err(Error::Degenerate(why))
}

/// The sign every direction of a vertex `a` of `m` takes against a
/// surface's gradient there, to first order and strictly (`Greater` where
/// the surface's function rises along every one, least there): every
/// edge's direction away from the vertex and every face's sector between its
/// two edge ends there on one side of the surface's tangent plane, decided
/// exactly against its gradient at the vertex (`corner_dirs`). `None` (no
/// contact of the vertex's own) where a sector is reflex or a half turn or
/// a face has other than two edge ends there (an edge closed at the
/// vertex, a face touching itself), where a direction lies within the
/// faces' parallel band of the tangent plane (`along`: the faces' own
/// tangency or incidence, U1's box on a rod's tangent plane), or where the
/// directions take both signs (an edge or a face heading toward the
/// surface: its own sections).
fn corner(m: &Prism, ends: &[(usize, bool)], quad: &Quad, a: &QV) -> Option<Ordering> {
    let dirs = corner_dirs(m, ends, a)?;
    signs(&dirs, &quad.gradient(a), None)
}

/// A vertex's tangent cone: every edge's direction away from the vertex,
/// where every face's sector between its two edge ends there is convex
/// (its second edge's direction on the face's side of its first: the
/// face's direction into it from the first edge, its outward normal across
/// the edge's running tangent), and so spanned by the two (the faces'
/// directions into them from their edges lie outside an acute sector: they
/// are not taken); `None` where a sector is reflex or a half turn or a
/// face has other than two edge ends there.
fn corner_dirs(m: &Prism, ends: &[(usize, bool)], a: &QV) -> Option<Vec<QV>> {
    // Each face's edge ends: the direction away from the vertex, the
    // edge's running tangent and the face's side of it.
    let mut at_face: BTreeMap<usize, Vec<(QV, QV, usize)>> = BTreeMap::new();
    let mut dirs: Vec<QV> = Vec::new();
    for &(ei, start) in ends {
        let t = running(m, ei, a);
        let away = if start {
            t.clone()
        } else {
            t.clone().map(|x| x.neg())
        };
        for (k, &f) in m.edges[ei].faces.iter().enumerate() {
            at_face
                .entry(f)
                .or_default()
                .push((away.clone(), t.clone(), k));
        }
        dirs.push(away);
    }
    for (&f, list) in &at_face {
        let [(_, t, k), (second, ..)] = list.as_slice() else {
            return None;
        };
        let n = m.normal_at(f, a);
        let into = if *k == 0 {
            qcross(&n, t)
        } else {
            qcross(t, &n)
        };
        if qqdot(&into, second).sign() != Ordering::Greater {
            return None;
        }
    }
    Some(dirs)
}

// ------------------------------------------------ edges and vertices

/// The inputs' vertices and edges within the resolution of each other's
/// (REVIEW_NOTES.md's "Edges and vertices within the resolution of each
/// other"): a vertex (or a cone's apex) and a vertex, a vertex and an
/// edge's point nearest it strictly inside the edge, and two edges' closest
/// points strictly inside both (lines, a line and a conic or a circle, two
/// conics or a conic and a circle),
/// each point at an end where its edge runs on through it smoothly too
/// (`smooth_end`). Each pair whose gap is within the resolution and not
/// zero (an incidence or an exact crossing, the existing rules') is decided
/// by each input's boundary at its point against the gap (`pair`): both
/// leaving it away from the other, `Degenerate` unless the gap lies inside
/// both inputs; both heading toward the other, `Degenerate`; otherwise no
/// contact of the pair's own. A cube's corner `1e-12` off a box's top edge
/// or corner, a box's edge `1e-12` across another's, a corner by a rod's
/// rim were fused with it into two solids. Each pair whose boxes meet once
/// one is widened by the resolution is tried, after `vertex_faces`.
pub(super) fn edges_and_vertices(models: &[Prism; 2], res: f64) -> Result<()> {
    let ends = [edge_ends(&models[0]), edge_ends(&models[1])];
    let spots: [Vec<Option<Spot>>; 2] = [0, 1].map(|o| {
        (0..models[o].verts.len())
            .map(|v| spot_of(&models[o], &ends[o], v))
            .collect()
    });
    let views: [Vec<[f64; 3]>; 2] =
        [0, 1].map(|o| models[o].verts.iter().map(|v| qv_f64(&v.p)).collect());
    let near = |x: &[f64; 3], y: &[f64; 3]| (0..3).all(|k| (x[k] - y[k]).abs() <= 2.0 * res);
    // Vertex and vertex.
    for (va, sa) in spots[0].iter().enumerate() {
        let Some(sa) = sa else { continue };
        for (vb, sb) in spots[1].iter().enumerate() {
            let Some(sb) = sb else { continue };
            if !near(&views[0][va], &views[1][vb]) {
                continue;
            }
            let a = End {
                m: &models[0],
                spot: sa,
                p: models[0].verts[va].p.clone(),
            };
            let b = End {
                m: &models[1],
                spot: sb,
                p: models[1].verts[vb].p.clone(),
            };
            pair(&a, &b, &Bracket::Exact, res, NEAR_VERTICES)?;
        }
    }
    let wide = |b: &([f64; 3], [f64; 3])| (b.0.map(|x| x - 2.0 * res), b.1.map(|x| x + 2.0 * res));
    let edge_box = |m: &Prism, e: &MEdge| intersect(&m.boxes[e.faces[0]], &m.boxes[e.faces[1]]);
    // Vertex and edge, each input's vertices against the other's edges.
    for o in 0..2 {
        let (vm, em) = (&models[o], &models[1 - o]);
        for (ei, e) in em.edges.iter().enumerate() {
            if !takes_edge(em, ei) {
                continue;
            }
            let b = wide(&edge_box(em, e));
            for (v, spot) in spots[o].iter().enumerate() {
                let Some(spot) = spot else { continue };
                let x = &views[o][v];
                if !(0..3).all(|k| b.0[k] <= x[k] && x[k] <= b.1[k]) {
                    continue;
                }
                let p = &vm.verts[v].p;
                for (foot, bracket) in nearest_to_point(em, ei, p, res) {
                    let (vertex, point) = (
                        End {
                            m: vm,
                            spot,
                            p: p.clone(),
                        },
                        End {
                            m: em,
                            spot: &Spot::Edge(ei),
                            p: foot,
                        },
                    );
                    // The gap from the first input's point to the second's.
                    let (a, b) = if o == 0 {
                        (&vertex, &point)
                    } else {
                        (&point, &vertex)
                    };
                    pair(a, b, &bracket, res, NEAR_VERTEX_EDGE)?;
                }
            }
        }
    }
    // Edge and edge.
    for (ea, e) in models[0].edges.iter().enumerate() {
        if !takes_edge(&models[0], ea) {
            continue;
        }
        let ba = wide(&edge_box(&models[0], e));
        for (eb, f) in models[1].edges.iter().enumerate() {
            if !takes_edge(&models[1], eb) || !boxes_meet(&ba, &edge_box(&models[1], f)) {
                continue;
            }
            for (x, y, bracket) in closest_on_edges(&models[0], ea, &models[1], eb, res) {
                let a = End {
                    m: &models[0],
                    spot: &Spot::Edge(ea),
                    p: x,
                };
                let b = End {
                    m: &models[1],
                    spot: &Spot::Edge(eb),
                    p: y,
                };
                pair(&a, &b, &bracket, res, NEAR_EDGES)?;
            }
        }
    }
    Ok(())
}

/// An input edge the rules take: with an input's id, not between faces on
/// one surface, its curve a line, a conic or a circle.
fn takes_edge(m: &Prism, ei: usize) -> bool {
    let e = &m.edges[ei];
    e.id.is_some()
        && !one_surface(m, e.faces[0], e.faces[1])
        && matches!(
            e.curve,
            Crv::Line { .. } | Crv::Conic { .. } | Crv::Circle(_)
        )
}

/// A point of an input's boundary for `pair`.
struct End<'a> {
    m: &'a Prism,
    spot: &'a Spot<'a>,
    p: QV,
}

/// How a pair's points were found: exactly, at a bracketed least of the
/// distance whose curve moves at most the bound over the bracket, or at a
/// binary64 least the exact signs did not bracket.
enum Bracket {
    Exact,
    Within(R),
    Unbracketed,
}

/// Whether the numbers of some points lie in one field (rational, or one
/// surd `sqrt(d)` over one base field), so that they combine.
fn one_field(points: &[&QV]) -> bool {
    let mut surd: Option<&R> = None;
    let mut base = None;
    for x in points.iter().flat_map(|p| p.iter()) {
        if let Some(d) = x.field() {
            if surd.is_some_and(|e| e != d) {
                return false;
            }
            surd = Some(d);
        }
        if let Some(g) = x.gen() {
            if base.is_some_and(|h| !Arc::ptr_eq(h, g)) {
                return false;
            }
            base = Some(g);
        }
    }
    true
}

/// The gap's own rounding: `rho^2`, `rho = 2^(e - 48)`, `2^e` the least
/// power of two at least one and every coordinate of the two points.
fn gap_rounding(a: &QV, b: &QV) -> R {
    let size = qv_f64(a)
        .iter()
        .chain(qv_f64(b).iter())
        .fold(1.0f64, |m, x| m.max(x.abs()));
    let e = if size.is_finite() {
        size.log2().ceil() as i32
    } else {
        0
    };
    R::from_integer(BigInt::from(2)).pow(2 * (e - 48))
}

/// The pair rule at the first input's point `a` and the second's `b`: their
/// gap `n = b - a` within the resolution and not zero; each input's
/// boundary at its point against `n` (`side`, the band `along_by` with the
/// gap's own rounding): the first's all against `n` and the second's all
/// with it (each leaving the other away: missing), `Degenerate` unless the
/// gap lies inside both inputs (each one's membership at its point pushed
/// across it); the reverse (each heading toward the other: crossing),
/// `Degenerate`; otherwise no contact. A bracketed distance above the
/// resolution but within its bound of it, or an unbracketed one within
/// twice it, is a `ComputationLimit` where the pair would be a contact.
fn pair(a: &End, b: &End, bracket: &Bracket, res: f64, why: &'static str) -> Result<()> {
    let Some((n, within)) = gap_within(&a.p, &b.p, bracket, res) else {
        return Ok(());
    };
    let unsure = !within;
    let Some(missing) = contact(a, b, &n) else {
        return Ok(());
    };
    if missing {
        let inside = [
            a.m.member(&a.p, std::slice::from_ref(&n)) == Loc::In,
            b.m.member(&b.p, &[n.clone().map(|x| x.neg())]) == Loc::In,
        ];
        if inside == [true, true] {
            return Ok(());
        }
    }
    if unsure {
        return Err(Error::ComputationLimit(NEAR_LIMIT));
    }
    Err(Error::Degenerate(why))
}

/// Whether the pair at `a` and `b` is a contact of its own against the gap
/// `n = b - a` (`pair`'s test): `Some(true)` where each input's boundary
/// leaves its point strictly away from the other's (missing), `Some(false)`
/// where each heads strictly toward it (crossing), `None` otherwise.
fn contact(a: &End, b: &End, n: &QV) -> Option<bool> {
    let rho2 = gap_rounding(&a.p, &b.p);
    let sa = side(a.m, a.spot, &a.p, n, Some(&rho2));
    let sb = side(b.m, b.spot, &b.p, n, Some(&rho2));
    match (sa, sb) {
        (Some(Ordering::Less), Some(Ordering::Greater)) => Some(true),
        (Some(Ordering::Greater), Some(Ordering::Less)) => Some(false),
        _ => None,
    }
}

/// `pair`'s gap from `a` to `b`, where it may be a near miss: `None` for
/// points of two fields, an incidence or an exact crossing (the existing
/// rules'), and a gap beyond the resolution (beyond it and the bracket's
/// reach); else the gap and whether it lies within the resolution (else
/// only within the bracket's reach of it, unsure).
fn gap_within(a: &QV, b: &QV, bracket: &Bracket, res: f64) -> Option<(QV, bool)> {
    if !one_field(&[a, b]) {
        return None;
    }
    let n = qsub(b, a);
    let nn = qqdot(&n, &n);
    if nn.sign() == Ordering::Equal {
        return None;
    }
    let tol = q(res);
    let above = |reach: &R| nn.add_r(&-(reach * reach)).sign() == Ordering::Greater;
    let within = !above(&tol);
    let unsure = !within
        && match bracket {
            Bracket::Exact => false,
            Bracket::Within(d) => !above(&(&tol + d)),
            Bracket::Unbracketed => !above(&(int(2) * &tol)),
        };
    (within || unsure).then_some((n, within))
}

/// Whether a point lies strictly inside edge `ei` of `m`, or at an end it
/// runs on through smoothly.
fn inside_edge(m: &Prism, ei: usize, x: &QV) -> bool {
    let e = &m.edges[ei];
    let (from, to, with) = edge_places(m, ei);
    match strictly_within(&place(&e.curve, x), &from, &to, with) {
        Some(true) => true,
        Some(false) => false,
        None => smooth_end(m, ei, x).is_some(),
    }
}

/// The points of edge `ei` of `m` nearest a point `p` inside the edge
/// (`inside_edge`): a line's foot, a circle's point toward `p` (exact, in
/// the quadratic field of the projection's length; `p` rational), and on a
/// conic that is no circle exactly each least of the distance bracketed
/// (`bracketed`). None for points of two fields.
fn nearest_to_point(m: &Prism, ei: usize, p: &QV, res: f64) -> Vec<(QV, Bracket)> {
    let e = &m.edges[ei];
    let found: Vec<(QV, Bracket)> = match &e.curve {
        Crv::Line { p: o, d } => {
            if !one_field(&[p, o]) {
                return Vec::new();
            }
            let t = qdot(&qsub(p, o), d).scale(&(int(1) / dot(d, d)));
            vec![(qadd(o, &qscale(d, &t)), Bracket::Exact)]
        }
        Crv::Circle(circ) => match rational_point(p) {
            Some(c) => super::graph::circle_nearest(&circ.c, &circ.x, &circ.y, &circ.r2, &c)
                .map(|x| (x, Bracket::Exact))
                .into_iter()
                .collect(),
            None => Vec::new(),
        },
        Crv::Conic { c, a, b } => {
            let circle = dot(a, b) == zero() && dot(a, a) == dot(b, b);
            match (circle, rational_point(p)) {
                (true, Some(x)) => super::graph::circle_nearest(c, a, b, &dot(a, a), &x)
                    .map(|x| (x, Bracket::Exact))
                    .into_iter()
                    .collect(),
                _ => {
                    let target = qv_f64(p);
                    let far = |x: [f64; 3]| (0..3).map(|k| (x[k] - target[k]).powi(2)).sum();
                    let slope = |x: &QV, t: &QV| qqdot(&qsub(x, p), t).sign();
                    bracketed(&e.curve, &far, &slope, res)
                }
            }
        }
        _ => Vec::new(),
    };
    found
        .into_iter()
        .filter(|(x, _)| inside_edge(m, ei, x))
        .collect()
}

/// A point's rational coordinates, where they are.
fn rational_point(p: &QV) -> Option<V> {
    let [x, y, z] = p;
    Some([
        x.rational()?.clone(),
        y.rational()?.clone(),
        z.rational()?.clone(),
    ])
}

/// The closest points of two edges inside both (`inside_edge`): of two
/// lines (skew; parallel lines have none), exactly; of a line and a conic
/// or a circle, each least of the distance from the curve to the line
/// bracketed on the curve (`bracketed`), with its exact foot on the line;
/// of two conics or circles, at each least found in binary64, decided
/// exactly there or bracketed (`conics_closest`).
fn closest_on_edges(
    ma: &Prism,
    ea: usize,
    mb: &Prism,
    eb: usize,
    res: f64,
) -> Vec<(QV, QV, Bracket)> {
    let (ca, cb) = (&ma.edges[ea].curve, &mb.edges[eb].curve);
    let found: Vec<(QV, QV, Bracket)> = match (ca, cb) {
        (Crv::Line { p: p1, d: d1 }, Crv::Line { p: p2, d: d2 }) => {
            if !one_field(&[p1, p2]) || lines_apart(p1, d1, p2, d2, res) {
                return Vec::new();
            }
            let (a, b, c) = (dot(d1, d1), dot(d1, d2), dot(d2, d2));
            let den = &a * &c - &b * &b;
            if den == zero() {
                return Vec::new();
            }
            let w = qsub(p1, p2);
            let (d, e) = (qdot(&w, d1), qdot(&w, d2));
            let inv = int(1) / &den;
            let s = e.scale(&b).sub(&d.scale(&c)).scale(&inv);
            let t = e.scale(&a).sub(&d.scale(&b)).scale(&inv);
            vec![(
                qadd(p1, &qscale(d1, &s)),
                qadd(p2, &qscale(d2, &t)),
                Bracket::Exact,
            )]
        }
        (Crv::Line { .. }, Crv::Conic { .. } | Crv::Circle(_)) => line_and_curve(cb, ca, res)
            .into_iter()
            .map(|(on_curve, foot, k)| (foot, on_curve, k))
            .collect(),
        (Crv::Conic { .. } | Crv::Circle(_), Crv::Line { .. }) => line_and_curve(ca, cb, res),
        (Crv::Conic { .. } | Crv::Circle(_), Crv::Conic { .. } | Crv::Circle(_)) => {
            conics_closest(ma, ea, mb, eb, res)
        }
        _ => Vec::new(),
    };
    // Only gaps `pair` takes (`gap_within`), told before the points' places
    // along the edges (which cost more).
    found
        .into_iter()
        .filter(|(x, y, k)| {
            gap_within(x, y, k, res).is_some() && inside_edge(ma, ea, x) && inside_edge(mb, eb, y)
        })
        .collect()
}

/// Whether two lines (points `p`, directions `d`) lie farther apart than
/// twice the resolution beyond doubt in binary64: `|w . (d1 x d2)|` above
/// `2 res |d1 x d2|` and a bound of binary64's error in it, the directions
/// not parallel to rounding (their closest points' gap, `pair`'s, then
/// exceeds the resolution).
fn lines_apart(p1: &QV, d1: &V, p2: &QV, d2: &V, res: f64) -> bool {
    let (a, b) = (qv_f64(p1), qv_f64(p2));
    let (u, v) = (
        d1.clone().map(|x| rational_f64(&x)),
        d2.clone().map(|x| rational_f64(&x)),
    );
    let c = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let w = [0, 1, 2].map(|k| a[k] - b[k]);
    let norm = |x: &[f64; 3]| (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    let (cc, uv) = (norm(&c), norm(&u) * norm(&v));
    let reach = norm(&a) + norm(&b) + norm(&w);
    let num = (w[0] * c[0] + w[1] * c[1] + w[2] * c[2]).abs();
    cc > 1e-12 * uv && num > 2.0 * res * cc + 1e-12 * reach * uv
}

/// The points of a conic or a circle `crv` at the leasts of its distance
/// from a line, each with its exact foot on the line.
fn line_and_curve(crv: &Crv, line: &Crv, res: f64) -> Vec<(QV, QV, Bracket)> {
    let Crv::Line { p, d } = line else {
        return Vec::new();
    };
    // A circle's points lie in fields of their own: its line's point
    // rational.
    if matches!(crv, Crv::Circle(_)) && rational_point(p).is_none() {
        return Vec::new();
    }
    let dd = dot(d, d);
    let foot = |x: &QV| {
        let t = qdot(&qsub(x, p), d).scale(&(int(1) / &dd));
        qadd(p, &qscale(d, &t))
    };
    let (fp, fd) = (qv_f64(p), d.clone().map(|x| rational_f64(&x)));
    let fdd = fd.iter().map(|x| x * x).sum::<f64>();
    // `|w x d|^2 / |d|^2`, without the cancellation of `|w|^2 - (w . d)^2 /
    // |d|^2` near the line.
    let far = |x: [f64; 3]| {
        let w = [0, 1, 2].map(|k| x[k] - fp[k]);
        let c = [
            w[1] * fd[2] - w[2] * fd[1],
            w[2] * fd[0] - w[0] * fd[2],
            w[0] * fd[1] - w[1] * fd[0],
        ];
        c.iter().map(|x| x * x).sum::<f64>() / fdd
    };
    let slope = |x: &QV, t: &QV| qqdot(&qsub(x, &foot(x)), t).sign();
    bracketed(crv, &far, &slope, res)
        .into_iter()
        .map(|(x, k)| {
            let y = foot(&x);
            (x, y, k)
        })
        .collect()
}

/// A conic's or a circle's point at the rational `(cos, sin)` of the half
/// angle's tangent `s` (about the opposite direction where `flip`), with
/// its tangent along the increasing angle.
fn curve_at(crv: &Crv, s: &R, flip: bool) -> Option<(QV, QV)> {
    let den = int(1) + s * s;
    let (cs, sn) = ((int(1) - s * s) / &den, int(2) * s / &den);
    let (cs, sn) = if flip { (-cs, -sn) } else { (cs, sn) };
    match crv {
        Crv::Conic { c, a, b } => {
            let u = [Qd::rat(cs), Qd::rat(sn)];
            Some((
                super::meet::conic_point(c, a, b, &u),
                super::meet::conic_tangent(a, b, &u),
            ))
        }
        Crv::Circle(circ) => {
            let x = circ.at(&[cs, sn]);
            let t = circ.tangent(&x);
            Some((x, t))
        }
        _ => None,
    }
}

/// A conic's or a circle's points in binary64 by the angle `t`, its exact
/// data rounded once (`bracketed` takes a few hundred points of each curve;
/// rounding the data at each point was most of the near-miss rule's time).
fn curve64(crv: &Crv) -> Option<impl Fn(f64) -> [f64; 3]> {
    let g = |v: &V| v.clone().map(|x| rational_f64(&x));
    // A conic's centre and axes, or a circle's centre, basis and radius.
    let (fc, fa, fb, r) = match crv {
        Crv::Conic { c, a, b } => (g(c), g(a), g(b), None),
        Crv::Circle(circ) => (
            g(&circ.c),
            g(&circ.x),
            g(&circ.y),
            Some(rational_f64(&circ.r2).sqrt()),
        ),
        _ => return None,
    };
    Some(move |t: f64| match r {
        None => [0, 1, 2].map(|k| fc[k] + fa[k] * t.cos() + fb[k] * t.sin()),
        Some(r) => {
            let dir = [0, 1, 2].map(|k| fa[k] * t.cos() + fb[k] * t.sin());
            let len = dir.iter().map(|x| x * x).sum::<f64>().sqrt();
            [0, 1, 2].map(|k| fc[k] + r * dir[k] / len)
        }
    })
}

/// A rational upper bound of `sqrt(x)`, `x >= 0`.
fn sqrt_up(x: &R) -> R {
    let mut r = q(rational_f64(x).sqrt() * (1.0 + 1e-9) + 1e-300);
    while &r * &r < *x {
        r = &r * int(2);
    }
    r
}

/// A rational lower bound of `sqrt(x)`, `x >= 0`.
fn sqrt_down(x: &R) -> R {
    let mut r = q(rational_f64(x).sqrt() * (1.0 - 1e-9));
    while sign(&r) == Ordering::Greater && &r * &r > *x {
        r = &r / int(2);
    }
    if sign(&r) == Ordering::Less {
        zero()
    } else {
        r
    }
}

/// A bound of a conic's or a circle's speed in its half angle's tangent:
/// `2 (|a| + |b|)` for `c + a cos + b sin` (the angle's rate in the tangent
/// at most two), and for a circle of radius `r` on its basis `(x, y)`
/// (orthogonal) `2 r max(|x|, |y|) / min(|x|, |y|)`.
fn speed(crv: &Crv) -> Option<R> {
    match crv {
        Crv::Conic { a, b, .. } => Some(int(2) * (sqrt_up(&dot(a, a)) + sqrt_up(&dot(b, b)))),
        Crv::Circle(circ) => {
            if dot(&circ.x, &circ.y) != zero() {
                return None;
            }
            let (xx, yy) = (dot(&circ.x, &circ.x), dot(&circ.y, &circ.y));
            let (hi, lo) = if xx > yy { (xx, yy) } else { (yy, xx) };
            let low = sqrt_down(&lo);
            if sign(&low) != Ordering::Greater {
                return None;
            }
            Some(int(2) * sqrt_up(&circ.r2) * sqrt_up(&hi) / low)
        }
        _ => None,
    }
}

/// The leasts of a function `far` (a squared distance) along a conic or a
/// circle: each sampled least among 64 angles narrowed in binary64 by a
/// golden-section search, then bracketed between two rational points of the
/// curve (dyadic tangents of the half angle, about the opposite direction
/// past a quarter turn) where the function's derivative along the curve,
/// `slope` at a point with its tangent, takes opposite exact signs, the
/// bracket narrowed by bisection on that sign. Each is the bracket's middle
/// point with the curve's movement over the bracket (`speed` times its
/// width: the least over the bracket lies within it of that point's); or,
/// where no bracket is found, the binary64 least's rational point; each
/// only where binary64 puts the least within twice the resolution.
fn bracketed(
    crv: &Crv,
    far: &dyn Fn([f64; 3]) -> f64,
    slope: &dyn Fn(&QV, &QV) -> Ordering,
    res: f64,
) -> Vec<(QV, Bracket)> {
    const N: usize = 64;
    let Some(at) = curve64(crv) else {
        return Vec::new();
    };
    let f = |t: f64| far(at(t));
    let step = TAU / N as f64;
    let vals: Vec<f64> = (0..N).map(|i| f(i as f64 * step)).collect();
    if vals.iter().any(|v| !v.is_finite()) {
        return Vec::new();
    }
    // The curve's speed, found at the first near miss (none without it).
    let mut speed_of: Option<R> = None;
    let mut out = Vec::new();
    for i in 0..N {
        let (prev, next) = (vals[(i + N - 1) % N], vals[(i + 1) % N]);
        if !(vals[i] <= prev && vals[i] < next) {
            continue;
        }
        let (mut lo, mut hi) = ((i as f64 - 1.0) * step, (i as f64 + 1.0) * step);
        let g = (5f64.sqrt() - 1.0) / 2.0;
        for _ in 0..80 {
            let (x1, x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
            if f(x1) < f(x2) {
                hi = x2;
            } else {
                lo = x1;
            }
        }
        let t = 0.5 * lo + 0.5 * hi;
        // Only near misses: binary64 within twice the resolution.
        let gap = f(t).max(0.0).sqrt();
        if gap.is_nan() || gap > 2.0 * res {
            continue;
        }
        if speed_of.is_none() {
            let Some(v) = speed(crv) else {
                return Vec::new();
            };
            speed_of = Some(v);
        }
        let Some(speed) = &speed_of else {
            return Vec::new();
        };
        let flip = t.cos() < 0.0;
        let s = (if flip { t - PI } else { t } / 2.0).tan();
        if let Some(found) = bisect(crv, s, flip, slope, speed) {
            out.push(found);
        } else if let Some((x, _)) = curve_at(crv, &q(s), flip) {
            out.push((x, Bracket::Unbracketed));
        }
    }
    out
}

/// The bracket about the half angle's tangent `s` (`bracketed`): opposite
/// exact signs of `slope` at `s - h` and `s + h` (falling, then rising),
/// `h` from `2^-20` widened up to `2^-4`, narrowed by 64 bisections; the
/// point at a middle where the slope vanishes, or at the last middle, with
/// the curve's movement over the bracket (any point of it bounds the least
/// inside it so).
fn bisect(
    crv: &Crv,
    s: f64,
    flip: bool,
    slope: &dyn Fn(&QV, &QV) -> Ordering,
    speed: &R,
) -> Option<(QV, Bracket)> {
    let at = |x: &R| curve_at(crv, x, flip).map(|(p, t)| (slope(&p, &t), p));
    let centre = q(s);
    let mut h = R::new(BigInt::from(1), BigInt::from(1u64 << 20));
    let (mut lo, mut hi) = loop {
        let (lo, hi) = (&centre - &h, &centre + &h);
        let (a, _) = at(&lo)?;
        let (b, _) = at(&hi)?;
        if (a, b) == (Ordering::Less, Ordering::Greater) {
            break (lo, hi);
        }
        h = &h * int(16);
        if h > R::new(BigInt::from(1), BigInt::from(16)) {
            return None;
        }
    };
    for _ in 0..64 {
        let mid = (&lo + &hi) / int(2);
        match at(&mid)? {
            (Ordering::Equal, p) => return Some((p, Bracket::Within(speed * (&hi - &lo)))),
            (Ordering::Less, _) => lo = mid,
            (Ordering::Greater, _) => hi = mid,
        }
    }
    let mid = (&lo + &hi) / int(2);
    let (_, p) = at(&mid)?;
    Some((p, Bracket::Within(speed * (&hi - &lo))))
}

// --------------------------------------------------------- two conics

/// A conic's or a circle's binary64 form `c + a cos + b sin` for the
/// searches (a circle's on its basis orthonormalized and scaled by its
/// radius: `Crv::Circle` places a point by the direction `x cos + y sin`
/// scaled onto it), with the map from its angle to `curve_at`'s.
struct Conic64 {
    c: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    /// A circle's basis lengths `(|x|, |y|)`.
    basis: Option<(f64, f64)>,
}

fn dot64(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub64(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

impl Conic64 {
    fn of(crv: &Crv) -> Option<Self> {
        let g = |v: &V| v.clone().map(|x| rational_f64(&x));
        match crv {
            Crv::Conic { c, a, b } => Some(Self {
                c: g(c),
                a: g(a),
                b: g(b),
                basis: None,
            }),
            Crv::Circle(circ) => {
                let (x, y) = (g(&circ.x), g(&circ.y));
                let (lx, ly) = (dot64(&x, &x).sqrt(), dot64(&y, &y).sqrt());
                let r = rational_f64(&circ.r2).sqrt();
                if !(lx > 0.0 && ly > 0.0 && r.is_finite()) {
                    return None;
                }
                Some(Self {
                    c: g(&circ.c),
                    a: x.map(|v| r * v / lx),
                    b: y.map(|v| r * v / ly),
                    basis: Some((lx, ly)),
                })
            }
            _ => None,
        }
    }

    /// The point and its first and second derivatives at the angle `t`.
    fn jet(&self, t: f64) -> [[f64; 3]; 3] {
        let (cs, sn) = (t.cos(), t.sin());
        let p = [0, 1, 2].map(|k| self.a[k] * cs + self.b[k] * sn);
        [
            [0, 1, 2].map(|k| self.c[k] + p[k]),
            [0, 1, 2].map(|k| self.b[k] * cs - self.a[k] * sn),
            p.map(|x| -x),
        ]
    }

    /// A bound of the speed in the angle, `|a| + |b|`.
    fn reach(&self) -> f64 {
        dot64(&self.a, &self.a).sqrt() + dot64(&self.b, &self.b).sqrt()
    }

    /// `curve_at`'s half-angle tangent and flip at the angle `t`.
    fn param(&self, t: f64) -> (f64, bool) {
        let t = match self.basis {
            None => t,
            Some((lx, ly)) => (t.sin() / ly).atan2(t.cos() / lx),
        };
        let flip = t.cos() < 0.0;
        let u = if flip { t - PI } else { t };
        ((u.sin().atan2(u.cos()) / 2.0).tan(), flip)
    }

    /// The angle of the curve's point nearest `p`, from `s` by Newton's
    /// method on the squared distance's derivative.
    fn nearest(&self, p: &[f64; 3], mut s: f64) -> f64 {
        for _ in 0..64 {
            let [y, d1, d2] = self.jet(s);
            let g = sub64(&y, p);
            let (f1, f2) = (dot64(&g, &d1), dot64(&d1, &d1) + dot64(&g, &d2));
            if f2.is_nan() || f2 <= 0.0 {
                break;
            }
            let ds = (-f1 / f2).clamp(-0.25, 0.25);
            if !ds.is_finite() {
                break;
            }
            s += ds;
            if ds.abs() <= 1e-16 * (1.0 + s.abs()) {
                break;
            }
        }
        s
    }
}

/// The pairs of angles at the local leasts of the distance between two
/// conics (their binary64 forms) within twice the resolution in binary64:
/// over a grid of 64 angles of each, each node at most its eight
/// neighbours (the angles wrapping) whose distance less the curves'
/// movement over a step (`reach` times it, each curve) is within twice the
/// resolution, refined by Newton's method on the gradient (`refine`), kept
/// where its distance is within twice the resolution, once. Every pair of
/// angles lies within half a step of a node in each, so the pair is passed
/// over whole where no node is within that reach.
fn conic_leasts(fa: &Conic64, fb: &Conic64, res: f64) -> Vec<(f64, f64)> {
    const N: usize = 64;
    let step = TAU / N as f64;
    let pa: Vec<[f64; 3]> = (0..N).map(|i| fa.jet(i as f64 * step)[0]).collect();
    let pb: Vec<[f64; 3]> = (0..N).map(|j| fb.jet(j as f64 * step)[0]).collect();
    let reach = 2.0 * res + (fa.reach() + fb.reach()) * step;
    let reach2 = reach * reach;
    let mut d2 = vec![0.0; N * N];
    let mut any = false;
    for i in 0..N {
        for j in 0..N {
            let w = sub64(&pa[i], &pb[j]);
            let v = dot64(&w, &w);
            if !v.is_finite() {
                return Vec::new();
            }
            any |= v <= reach2;
            d2[i * N + j] = v;
        }
    }
    let mut out: Vec<(f64, f64)> = Vec::new();
    if !any {
        return out;
    }
    let close = |p: f64, q: f64| {
        let d = (p - q).rem_euclid(TAU);
        d.min(TAU - d) <= 1e-9
    };
    for i in 0..N {
        for j in 0..N {
            let v = d2[i * N + j];
            if v > reach2 {
                continue;
            }
            let least = [N - 1, 0, 1].iter().all(|di| {
                [N - 1, 0, 1]
                    .iter()
                    .all(|dj| v <= d2[((i + di) % N) * N + (j + dj) % N])
            });
            if !least {
                continue;
            }
            let (t, s) = refine(fa, fb, i as f64 * step, j as f64 * step);
            let w = sub64(&fa.jet(t)[0], &fb.jet(s)[0]);
            let gap = dot64(&w, &w).sqrt();
            if gap.is_nan() || gap > 2.0 * res {
                continue;
            }
            if out.iter().any(|&(t0, s0)| close(t, t0) && close(s, s0)) {
                continue;
            }
            out.push((t, s));
        }
    }
    out
}

/// A least of two conics' distance about the angles `(t, s)`, by Newton's
/// method on the gradient of `|A(t) - B(s)|^2`; one angle's step at a time
/// where the Hessian is not definite (tangent or coincident curves).
fn refine(fa: &Conic64, fb: &Conic64, mut t: f64, mut s: f64) -> (f64, f64) {
    for it in 0..100 {
        let [x, xa, xaa] = fa.jet(t);
        let [y, yb, ybb] = fb.jet(s);
        let g = sub64(&x, &y);
        let (gt, gs) = (dot64(&g, &xa), -dot64(&g, &yb));
        let h11 = dot64(&xa, &xa) + dot64(&g, &xaa);
        let h22 = dot64(&yb, &yb) - dot64(&g, &ybb);
        let h12 = -dot64(&xa, &yb);
        let det = h11 * h22 - h12 * h12;
        let (dt, ds) = if h11 > 0.0 && h22 > 0.0 && det > 1e-12 * h11 * h22 {
            ((h12 * gs - h22 * gt) / det, (h12 * gt - h11 * gs) / det)
        } else if it % 2 == 0 && h11 > 0.0 {
            (-gt / h11, 0.0)
        } else if h22 > 0.0 {
            (0.0, -gs / h22)
        } else {
            break;
        };
        if !(dt.is_finite() && ds.is_finite()) {
            break;
        }
        let (dt, ds) = (dt.clamp(-0.25, 0.25), ds.clamp(-0.25, 0.25));
        t += dt;
        s += ds;
        if dt.abs() <= 1e-16 * (1.0 + t.abs()) && ds.abs() <= 1e-16 * (1.0 + s.abs()) {
            break;
        }
    }
    (t, s)
}

/// The closest points of two conic or circle edges of the inputs
/// (REVIEW_NOTES.md's "Two conic edges within the resolution of each
/// other"; two `Crv::Circle`s not: their points lie in fields of their
/// own): at each least found in binary64 (`conic_leasts`) the curves' exact
/// points at its angles' rational half-angle tangents (`curve_at`), none
/// where their gap is zero or of rounding's size (`|n| <= rho`,
/// `gap_rounding`: every direction then along the gap's plane, no contact);
/// within the resolution, those points (the least at most their gap);
/// beyond it, where the pair is a contact at those points inside both
/// edges, the least bracketed (`conics_bracketed`), or those points
/// unbracketed.
fn conics_closest(
    ma: &Prism,
    ea: usize,
    mb: &Prism,
    eb: usize,
    res: f64,
) -> Vec<(QV, QV, Bracket)> {
    let (ca, cb) = (&ma.edges[ea].curve, &mb.edges[eb].curve);
    if matches!((ca, cb), (Crv::Circle(_), Crv::Circle(_))) {
        return Vec::new();
    }
    let (Some(fa), Some(fb)) = (Conic64::of(ca), Conic64::of(cb)) else {
        return Vec::new();
    };
    let tol = q(res);
    let mut out = Vec::new();
    for (t, s) in conic_leasts(&fa, &fb, res) {
        let ((ua, fla), (ub, flb)) = (fa.param(t), fb.param(s));
        let (Some((x, _)), Some((y, _))) = (curve_at(ca, &q(ua), fla), curve_at(cb, &q(ub), flb))
        else {
            continue;
        };
        if !one_field(&[&x, &y]) {
            continue;
        }
        let n = qsub(&y, &x);
        let nn = qqdot(&n, &n);
        if nn.add_r(&-gap_rounding(&x, &y)).sign() != Ordering::Greater {
            continue;
        }
        if nn.add_r(&-(&tol * &tol)).sign() != Ordering::Greater {
            out.push((x, y, Bracket::Exact));
            continue;
        }
        let (a, b) = (
            End {
                m: ma,
                spot: &Spot::Edge(ea),
                p: x.clone(),
            },
            End {
                m: mb,
                spot: &Spot::Edge(eb),
                p: y.clone(),
            },
        );
        if contact(&a, &b, &n).is_none() || !inside_edge(ma, ea, &x) || !inside_edge(mb, eb, &y) {
            continue;
        }
        out.push(conics_bracketed((ca, &fa, t), (cb, &fb, s)).unwrap_or((
            x,
            y,
            Bracket::Unbracketed,
        )));
    }
    out
}

/// The inner curve of `conics_bracketed`: a circle exactly (its centre,
/// basis and squared radius), or a conic bracketed about a binary64 angle
/// (with its speed bound).
enum Inner<'a> {
    Round(&'a V, &'a V, &'a V, R),
    Conic(&'a Crv, &'a Conic64, f64, R),
}

impl Inner<'_> {
    /// The inner curve's point nearest the rational point `x` and a bound
    /// of its distance from the nearest's: exactly on a circle
    /// (`circle_nearest`), else bracketed on exact signs of `(y - x) . t`
    /// (`bisect`) about the binary64 nearest angle.
    fn nearest(&self, x: &QV) -> Option<(QV, R)> {
        let xr = rational_point(x)?;
        match self {
            Inner::Round(c, a, b, r2) => {
                super::graph::circle_nearest(c, a, b, r2, &xr).map(|y| (y, zero()))
            }
            Inner::Conic(crv, f, s, speed) => {
                let s = f.nearest(&qv_f64(x), *s);
                let (u, flip) = f.param(s);
                let slope = |y: &QV, t: &QV| qqdot(&qsub(y, x), t).sign();
                match bisect(crv, u, flip, &slope, speed)? {
                    (y, Bracket::Within(r)) => Some((y, r)),
                    _ => None,
                }
            }
        }
    }
}

/// Two conics' (or a conic's and a circle's) closest points about the
/// angles found in binary64, bracketed: on the outer curve (a conic:
/// rational points `x`), the exact distance to the inner one, its
/// derivative's sign that of `(x - y) . t_A` at the inner curve's point
/// `y` nearest `x` (`Inner::nearest`; the inner a circle exactly where
/// either is), taken as zero within `|t_A| r_B` of it (`r_B` the bound of
/// `y`'s distance from the nearest's), the outer bracket narrowed by
/// `bisect` on it; the points (the first curve's first) and the bound
/// `D_A + r_B` (the distance 1-Lipschitz in each point: the least over the
/// bracket at least their gap less it).
fn conics_bracketed(
    first: (&Crv, &Conic64, f64),
    second: (&Crv, &Conic64, f64),
) -> Option<(QV, QV, Bracket)> {
    let round = |c: &Crv| match c {
        Crv::Circle(_) => true,
        Crv::Conic { a, b, .. } => dot(a, b) == zero() && dot(a, a) == dot(b, b),
        _ => false,
    };
    let swap = match (first.0, second.0) {
        (Crv::Circle(_), _) => true,
        (_, Crv::Circle(_)) => false,
        (a, b) => round(a) && !round(b),
    };
    let ((co, fo, to), (ci, fi, si)) = if swap {
        (second, first)
    } else {
        (first, second)
    };
    let speed_o = speed(co)?;
    let inner = match ci {
        Crv::Circle(circ) => Inner::Round(&circ.c, &circ.x, &circ.y, circ.r2.clone()),
        Crv::Conic { c, a, b } if round(ci) => Inner::Round(c, a, b, dot(a, a)),
        _ => Inner::Conic(ci, fi, si, speed(ci)?),
    };
    let slope = |x: &QV, ta: &QV| {
        let Some((y, r)) = inner.nearest(x) else {
            return Ordering::Equal;
        };
        let phi = qqdot(&qsub(x, &y), ta);
        let lim = qqdot(ta, ta).scale(&(&r * &r));
        if phi.mul(&phi).sub(&lim).sign() != Ordering::Greater {
            return Ordering::Equal;
        }
        phi.sign()
    };
    let (u, flip) = fo.param(to);
    let (x, Bracket::Within(da)) = bisect(co, u, flip, &slope, &speed_o)? else {
        return None;
    };
    let (y, r) = inner.nearest(&x)?;
    let bound = Bracket::Within(da + r);
    Some(if swap { (y, x, bound) } else { (x, y, bound) })
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::identity::OperationId;
    use crate::profile::boolean::Operand;
    use crate::{Frame3, Point3, Solid, Tolerance, Vec3};

    fn level() -> Frame3 {
        Frame3::new(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Tolerance::default(),
        )
        .unwrap()
    }

    fn rational(p: &QV) -> V {
        rational_point(p).unwrap()
    }

    /// A rod's rim runs on through its seam (two arcs of one input edge,
    /// the wall's seam edge between them) and a cone's rim through its
    /// vertex (a closed edge); a box's edge does not run on through its
    /// corner, and a point inside an edge is no end.
    #[test]
    fn edges_run_on_through_seams_and_closed_vertices() {
        let tol = Tolerance::default();
        let seam = R::new(BigInt::from(2), BigInt::from(7));
        let rod = Solid::cylinder_with(OperationId(1), level(), 1.0, -3.0, 0.0, tol)
            .unwrap()
            .0;
        let cone = Solid::cone_with(OperationId(1), level(), 1.0, 0.0, 1.0, tol)
            .unwrap()
            .0;
        let models = [
            Prism::new(&rod, Operand::A, &seam).unwrap(),
            super::super::cone::model(&cone, Operand::A, &seam).unwrap(),
        ];
        for m in &models {
            let rims: Vec<usize> = (0..m.edges.len())
                .filter(|&ei| {
                    m.edges[ei].id.is_some() && matches!(m.edges[ei].curve, Crv::Conic { .. })
                })
                .collect();
            assert!(!rims.is_empty());
            for ei in rims {
                let e = &m.edges[ei];
                for v in [e.start, e.end] {
                    assert_eq!(smooth_end(m, ei, &m.verts[v].p), Some(v));
                }
            }
        }
        // The cone's apex has no edges.
        let apex = models[1].funnel.as_ref().and_then(|f| f.apex).unwrap();
        assert!(models[1]
            .edges
            .iter()
            .all(|e| e.start != apex && e.end != apex));
        let boxed =
            Solid::box_at(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0), tol).unwrap();
        let m = Prism::new(&boxed, Operand::B, &seam).unwrap();
        for (ei, e) in m.edges.iter().enumerate() {
            assert_eq!(smooth_end(&m, ei, &m.verts[e.start].p), None);
            let (s, t) = (rational(&m.verts[e.start].p), rational(&m.verts[e.end].p));
            let mid = qv(&scale(&add(&s, &t), &(int(1) / int(2))));
            assert_eq!(smooth_end(&m, ei, &mid), None);
        }
    }

    /// A least distance from a point to a conic of a turned frame's kind (no
    /// circle exactly) is bracketed: the exact distance at the bracket's
    /// point is the least to far below the resolution, the bracket's bound
    /// too.
    #[test]
    fn a_conics_least_distance_is_bracketed() {
        let crv = Crv::Conic {
            c: [zero(), zero(), zero()],
            a: [q(1.0 + 2f64.powi(-50)), zero(), q(2f64.powi(-52))],
            b: [zero(), q(1.0), zero()],
        };
        let p = [q(0.6), q(0.8), q(1e-12)];
        let target = p.clone().map(|x| rational_f64(&x));
        let far = |x: [f64; 3]| (0..3).map(|k| (x[k] - target[k]).powi(2)).sum::<f64>();
        let pq = qv(&p);
        let slope = |x: &QV, t: &QV| qqdot(&qsub(x, &pq), t).sign();
        let found = bracketed(&crv, &far, &slope, 1e-7);
        assert_eq!(found.len(), 1);
        let (x, Bracket::Within(d)) = &found[0] else {
            panic!("unbracketed")
        };
        assert!(rational_f64(d) < 1e-15, "{}", rational_f64(d));
        let w = qsub(x, &pq);
        let u = qqdot(&w, &w).to_f64().sqrt();
        // `(0.6, 0.8)` lies on the circle but for the conic's ulps: its
        // distance is about the point's height.
        assert!((u - 1e-12).abs() < 1e-15, "{u}");
        assert!(sqrt_up(&q(2.0)) * sqrt_up(&q(2.0)) >= q(2.0));
        assert!(sqrt_down(&q(2.0)) * sqrt_down(&q(2.0)) <= q(2.0));
    }

    /// Two conics' least distance is found in binary64 and bracketed: a
    /// turned frame's unit circle (no circle exactly) and a circle of
    /// radius 0.5 standing over its point `(0.6, 0.8, 0)` across its
    /// tangent, `1e-9` above it (a circle exactly, the inner curve's point
    /// exact; and within rounding of one, bracketed too): one least each,
    /// the bracket's points `1e-9` apart and its bound far below the
    /// resolution; `0.1` apart, none within twice the resolution.
    #[test]
    fn two_conics_least_distance_is_bracketed() {
        let res = Tolerance::default().linear();
        let a = Crv::Conic {
            c: [zero(), zero(), zero()],
            a: [q(1.0 + 2f64.powi(-50)), zero(), q(2f64.powi(-52))],
            b: [zero(), q(1.0), zero()],
        };
        let other = |d: f64, tilt: f64| Crv::Conic {
            c: [q(0.6), q(0.8), q(0.5) + q(d)],
            a: [zero(), zero(), q(-0.5 - tilt)],
            b: [q(0.3), q(0.4), zero()],
        };
        let fa = Conic64::of(&a).unwrap();
        for tilt in [0.0, 2f64.powi(-50)] {
            let b = other(1e-9, tilt);
            let fb = Conic64::of(&b).unwrap();
            let found = conic_leasts(&fa, &fb, res);
            assert_eq!(found.len(), 1, "{tilt}");
            let (t, s) = found[0];
            let (x, y, Bracket::Within(d)) = conics_bracketed((&a, &fa, t), (&b, &fb, s)).unwrap()
            else {
                panic!("unbracketed")
            };
            assert!(rational_f64(&d) < 1e-15, "{}", rational_f64(&d));
            let w = qsub(&y, &x);
            let u = qqdot(&w, &w).to_f64().sqrt();
            assert!((u - 1e-9).abs() < 1e-14, "{u}");
            let fb = Conic64::of(&other(0.1, tilt)).unwrap();
            assert!(conic_leasts(&fa, &fb, res).is_empty());
        }
    }

    /// `closest_on_edges` passes over two lines only farther apart than
    /// twice the resolution (`lines_apart`): skew lines five resolutions
    /// apart at right angles or at `1e-3` rad, a thousand from the origin
    /// too at right angles; at the resolution, within it, or parallel to
    /// rounding, never; nor where binary64's error bound reaches the gap
    /// (`1e-6` rad a thousand from the origin).
    #[test]
    fn lines_apart_beyond_twice_the_resolution() {
        let res = Tolerance::default().linear();
        let v = |x: [f64; 3]| x.map(q);
        for (far, turn, beyond) in [
            (0.0, 1.0, true),
            (0.0, 1e-3, true),
            (1000.0, 1.0, true),
            (1000.0, 1e-6, false),
        ] {
            for (gap, apart) in [(5.0, beyond), (1.0, false), (0.5, false), (0.0, false)] {
                let p1 = qv(&v([far, 0.25, 0.0]));
                let p2 = qv(&v([far + 0.5, -0.75, gap * res]));
                let (d1, d2) = (v([1.0, 0.0, 0.0]), v([1.0, turn, 0.0]));
                assert_eq!(
                    lines_apart(&p1, &d1, &p2, &d2, res),
                    apart,
                    "{far} {turn} {gap}"
                );
            }
        }
        let p1 = qv(&v([0.0, 0.0, 0.0]));
        let p2 = qv(&v([0.0, 1.0, 1.0]));
        let d = v([1.0, 2.0, 3.0]);
        let tilted = v([1.0, 2.0, 3.0 + 1e-15]);
        assert!(!lines_apart(&p1, &d, &p2, &d, res));
        assert!(!lines_apart(&p1, &d, &p2, &tilted, res));
    }

    /// `edge_faces` passes over an exact extremum binary64 puts clearly off
    /// a plane or a cylinder (`Quad::clearly_off`) only beyond the
    /// resolution: a point within it or at it of a rod's wall or top, on
    /// either side, never; one five resolutions off, always; off a cone's
    /// wall, never.
    #[test]
    fn points_clearly_off_lie_beyond_the_resolution() {
        let tol = Tolerance::default();
        let res = tol.linear();
        let seam = R::new(BigInt::from(2), BigInt::from(7));
        let rod = Solid::cylinder_with(OperationId(1), level(), 1.0, -3.0, 0.0, tol)
            .unwrap()
            .0;
        let cone = Solid::cone_with(OperationId(1), level(), 1.0, 0.0, 1.0, tol)
            .unwrap()
            .0;
        let m = Prism::new(&rod, Operand::A, &seam).unwrap();
        let (mut walls, mut tops) = (0, 0);
        for fi in 0..m.faces.len() {
            let Some(quad) = Quad::of(&m, fi) else {
                continue;
            };
            for (d, off) in [
                (0.5, false),
                (1.0, false),
                (-0.5, false),
                (-1.0, false),
                (5.0, true),
                (-5.0, true),
            ] {
                let x = match quad.kind {
                    Kind::Cyl { .. } => [-(1.0 + d * res) * 0.6, (1.0 + d * res) * 0.8, -1.5],
                    Kind::Plane { .. } if quad.rows[0][2] > 0.0 => [0.25, -0.5, d * res],
                    _ => continue,
                };
                assert_eq!(quad.clearly_off(x, res), off, "{fi} {d}");
            }
            match quad.kind {
                Kind::Cyl { .. } => walls += 1,
                Kind::Plane { .. } if quad.rows[0][2] > 0.0 => tops += 1,
                _ => {}
            }
        }
        assert!(walls > 0 && tops > 0);
        let m = super::super::cone::model(&cone, Operand::A, &seam).unwrap();
        for fi in 0..m.faces.len() {
            let Some(quad) = Quad::of(&m, fi) else {
                continue;
            };
            if matches!(quad.kind, Kind::Cone { .. }) {
                assert!(!quad.clearly_off([0.0, 0.0, 5.0], res));
            }
        }
    }
}
