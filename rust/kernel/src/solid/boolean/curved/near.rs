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
//! of each other"; edges of other curves, points of two fields and a cone's
//! apex against a spline wall or a torus, "Other curves, points of two
//! fields and an apex against a spline or torus face"), and the
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
use crate::certified::{Fast, Interval as I, Real};
use crate::cr::Rn as _;
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
/// A near miss between points of two fields whose gap lies along a face's
/// tangent plane at its point: its membership not decided.
pub(super) const NEAR_FIELDS: &str = "a near miss's gap along a face between points of two fields";

/// A face's surface: a plane, or a cylinder or a cone on its model's frame
/// (its local coordinates `(u, v, w)`); for the edges' and vertices' rules
/// against tori and spline walls (REVIEW_NOTES.md, "Faces with no near-miss
/// rule of their own"), a torus or a spline wall too, their functions their
/// signed distances.
enum Kind<'a> {
    Plane {
        p: &'a V,
        m: &'a V,
    },
    Cyl {
        f: &'a Affine,
        c: &'a P2,
        r: &'a R,
    },
    Cone {
        f: &'a Affine,
        b: &'a R,
        k: &'a R,
    },
    Torus {
        f: &'a Affine,
        ring: &'a super::torus::Ring,
    },
    /// A spline wall, face `fi` of `m` (its outward normal `normal_at`'s),
    /// its arcs in binary64 for its searches.
    Spline {
        f: &'a Affine,
        seg: &'a super::spline_walls::SplineSeg,
        m: &'a Prism,
        fi: usize,
        arcs: Vec<PArc<f64>>,
    },
}

/// A face's surface and its data in binary64 (for searches only): a
/// plane's point and normal, or the frame's map to local coordinates and
/// a cylinder's centre and radius, a cone's `b` and `k` or a torus's radii.
struct Quad<'a> {
    kind: Kind<'a>,
    o: [f64; 3],
    rows: [[f64; 3]; 3],
    data: [f64; 3],
}

/// `unit_at`'s half angle's tangent at `t` and whether it is taken about
/// the opposite direction.
fn unit_param(t: f64) -> (R, bool) {
    let flip = t.cos_rn() < 0.0;
    (q((if flip { t - PI } else { t } / 2.0).tan_rn()), flip)
}

/// A rational point of the unit circle at about the angle `t` (its half
/// angle's tangent a rational, about the opposite direction near a half
/// turn).
fn unit_at(t: f64) -> [R; 2] {
    let (s, flip) = unit_param(t);
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

    /// A torus's or a spline wall's face (the edges' and vertices' rules'
    /// second pass); `None` for the others.
    fn of_curved(m: &'a Prism, fi: usize) -> Option<Self> {
        let (vm, vi) = m.view(fi);
        let g = |v: &V| [0, 1, 2].map(|k| rational_f64(&v[k]));
        let o = g(&vm.f.o);
        let rows = [g(vm.f.row(0)), g(vm.f.row(1)), g(vm.f.row(2))];
        let (kind, data) = match &vm.faces[vi].surf {
            Surf::Torus => {
                let ring = vm.ring.as_ref()?;
                (
                    Kind::Torus { f: &vm.f, ring },
                    [rational_f64(&ring.big), rational_f64(&ring.small), 0.0],
                )
            }
            Surf::Spline(seg) => (
                Kind::Spline {
                    f: &vm.f,
                    seg,
                    m,
                    fi,
                    arcs: parcs(seg),
                },
                [0.0; 3],
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
            // The torus's quartic, negative inside.
            Kind::Torus { f, ring } => ring.value(&f.local_q(x)),
            // The offset along the outward normal at the wall's point near
            // the foot (zero, on the surface, where none is found).
            Kind::Spline { m, fi, .. } => match self.foot(x) {
                Some(s) => qqdot(&qsub(x, &s), &m.normal_at(*fi, &s)),
                None => Qd::rat(zero()),
            },
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
            Kind::Torus { f, ring } => {
                let lg = ring.gradient(&f.local_q(x));
                [0, 1, 2].map(|j| {
                    lg[0]
                        .scale(&f.row(0)[j])
                        .add(&lg[1].scale(&f.row(1)[j]))
                        .add(&lg[2].scale(&f.row(2)[j]))
                })
            }
            // The outward normal at the wall's point near the foot.
            Kind::Spline { m, fi, .. } => match self.foot(x) {
                Some(s) => m.normal_at(*fi, &s),
                None => qv(&[zero(), zero(), zero()]),
            },
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
            // Signed distances and their unit gradients.
            Kind::Torus { .. } => {
                let l = self.local64(x);
                let [big, small, _] = self.data;
                let rho = l[0].hypot_rn(l[1]);
                let (e, z) = (rho - big, l[2]);
                let d = e.hypot_rn(z);
                let lg = [e / d * l[0] / rho, e / d * l[1] / rho, z / d];
                let g =
                    [0, 1, 2].map(|j| lg[0] * rows[0][j] + lg[1] * rows[1][j] + lg[2] * rows[2][j]);
                (d - small, g)
            }
            Kind::Spline { seg, arcs, .. } => {
                let l = self.local64(x);
                let (d, n) = profile_distance(arcs, seg.hole, [l[0], l[1]]);
                let g = [0, 1, 2].map(|j| n[0] * rows[0][j] + n[1] * rows[1][j]);
                (d, g)
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
            // A distance.
            Kind::Torus { .. } | Kind::Spline { .. } => v.abs() > 4.0 * res,
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
                let t = (l[1] - self.data[1]).atan2_rn(l[0] - self.data[0]);
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
                let rho = l[0].hypot_rn(l[1]);
                let h = q((l[2] + kf * (rho - bf)) / (1.0 + kf * kf));
                let rad = *b + *k * &h;
                if sign(&rad) != Ordering::Greater {
                    return None;
                }
                let [cs, sn] = unit_at(l[1].atan2_rn(l[0]));
                let loc = [Qd::rat(&rad * &cs), Qd::rat(&rad * &sn), Qd::rat(h)];
                Some(world(f, &loc))
            }
            // A torus's exact point at the rational half-angle tangents of
            // the binary64 foot's angles (as `apex_face`'s).
            Kind::Torus { f, ring } => {
                let l = self.local64(qv_f64(x));
                let rho = l[0].hypot_rn(l[1]);
                let (u, v) = (l[1].atan2_rn(l[0]), l[2].atan2_rn(rho - self.data[0]));
                if !(u.is_finite() && v.is_finite() && rho > 0.0) {
                    return None;
                }
                let ([cu, su], [cv, sv]) = (unit_at(u), unit_at(v));
                let rad = &ring.big + &ring.small * &cv;
                Some(qv(&f.point(
                    &(&rad * &cu),
                    &(&rad * &su),
                    &(&ring.small * &sv),
                )))
            }
            // A spline wall's exact point at a rational run parameter near the
            // binary64 nearest profile point, at the point's height.
            Kind::Spline { f, seg, arcs, .. } => {
                let l = self.local64(qv_f64(x));
                let tau = near_run(arcs, [l[0], l[1]])?;
                let uv = seg.point(&Qd::rat(run_rational(seg, tau)));
                let h = f.local_q(x);
                Some(super::spline_walls::qpoint(f, &uv[0], &uv[1], &h[2]))
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
        _ => return None,
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
        _ => unreachable!("taken above"),
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
    edge_faces_with(models, res, false)
}

/// `edge_faces` against plane, cylinder and cone faces, or (`curved`)
/// against torus and spline faces (`curved_edges_and_vertices`).
fn edge_faces_with(models: &[Prism; 2], res: f64, curved: bool) -> Result<()> {
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
                let quad = if curved {
                    Quad::of_curved(qm, fq)
                } else {
                    Quad::of(qm, fq)
                };
                let Some(quad) = quad else {
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
                    None if curved && matches!(e.curve, Crv::Line { .. }) => {
                        line_extrema(em, ei, &quad)
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
        // Over the edge's range (`line_extrema`).
        (Crv::Line { .. }, Kind::Torus { .. } | Kind::Spline { .. }) => Vec::new(),
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
                _ => unreachable!("taken above"),
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
                let (cs, sn) = (t.cos_rn(), t.sin_rn());
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
                let (cs, sn) = (t.cos_rn() * rad / lx, t.sin_rn() * rad / ly);
                let (dc, ds) = (-t.sin_rn() * rad / lx, t.cos_rn() * rad / ly);
                (
                    [0, 1, 2].map(|k| fc[k] + fx[k] * cs + fy[k] * sn),
                    [0, 1, 2].map(|k| fx[k] * dc + fy[k] * ds),
                )
            };
            search(&at, quad)
                .into_iter()
                .map(|(t, kind)| {
                    let x = circ.at(&[q(t.cos_rn() / lx), q(t.sin_rn() / ly)]);
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
    vertex_faces_with(models, res, false)
}

/// `vertex_faces` against plane, cylinder and cone faces (and a cone's apex
/// against the others, `apex_face`), or (`curved`) the vertices but a
/// cone's apex against torus and spline faces
/// (`curved_edges_and_vertices`).
fn vertex_faces_with(models: &[Prism; 2], res: f64, curved: bool) -> Result<()> {
    let wide = |b: &([f64; 3], [f64; 3])| (b.0.map(|x| x - res), b.1.map(|x| x + res));
    for o in 0..2 {
        let (vm, qm) = (&models[o], &models[1 - o]);
        let ends = edge_ends(vm);
        for (v, mv) in vm.verts.iter().enumerate() {
            let Some(spot) = spot_of(vm, &ends, v) else {
                continue;
            };
            if curved && matches!(spot, Spot::Apex) {
                continue;
            }
            let x = qv_f64(&mv.p);
            for fq in 0..qm.faces.len() {
                let b = wide(&qm.boxes[fq]);
                if !(0..3).all(|k| b.0[k] <= x[k] && x[k] <= b.1[k]) {
                    continue;
                }
                if curved {
                    if let Some(quad) = Quad::of_curved(qm, fq) {
                        vertex_face(vm, &spot, &mv.p, qm, fq, &quad, res)?;
                    }
                    continue;
                }
                let Some(quad) = Quad::of(qm, fq) else {
                    if matches!(spot, Spot::Apex) {
                        apex_face(vm, &mv.p, qm, fq, res)?;
                    }
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

/// A cone's apex `a` of `vm` within the resolution of a torus or a spline
/// wall face `fq` of `qm` (REVIEW_NOTES.md, "Other curves, points of two
/// fields and an apex against a spline or torus face"): the surface's
/// exact point near the apex's foot (a torus's at the rational half-angle
/// tangents of the binary64 foot's angles; a spline wall's at a rational
/// run parameter near the binary64 nearest profile point, at the apex's
/// height), within the resolution and in the face; the surface's function
/// and gradient at the apex (the torus's quartic; a spline wall's offset
/// along its outward normal at that point, the apex on the wall exactly an
/// incidence); then `vertex_face`'s rule with the apex's rulings.
fn apex_face(vm: &Prism, a: &QV, qm: &Prism, fq: usize, res: f64) -> Result<()> {
    let (view, vi) = qm.view(fq);
    let f = &view.f;
    let l = f.local_q(a);
    let l64 = qv_f64(&l);
    let (s, value, g) = match &view.faces[vi].surf {
        Surf::Torus => {
            let Some(ring) = view.ring.as_ref() else {
                return Ok(());
            };
            let (big, small) = (rational_f64(&ring.big), rational_f64(&ring.small));
            let u = l64[1].atan2_rn(l64[0]);
            let v = l64[2].atan2_rn(l64[0].hypot_rn(l64[1]) - big);
            if !(u.is_finite() && v.is_finite() && small > 0.0) {
                return Ok(());
            }
            let ([cu, su], [cv, sv]) = (unit_at(u), unit_at(v));
            let rad = &ring.big + &ring.small * &cv;
            let s = qv(&f.point(&(&rad * &cu), &(&rad * &su), &(&ring.small * &sv)));
            let lg = ring.gradient(&l);
            let g = [0, 1, 2].map(|j| {
                lg[0]
                    .scale(&f.row(0)[j])
                    .add(&lg[1].scale(&f.row(1)[j]))
                    .add(&lg[2].scale(&f.row(2)[j]))
            });
            (s, ring.value(&l).sign(), g)
        }
        Surf::Spline(seg) => {
            if seg.locate(&[l[0].clone(), l[1].clone()]).is_some() {
                // On the wall exactly: the incidences' own rules.
                return Ok(());
            }
            let tau = Qd::rat(q(seg.tau_near(l64[0], l64[1])));
            let uv = seg.point(&tau);
            let s = super::spline_walls::qpoint(f, &uv[0], &uv[1], &l[2]);
            let m = qm.normal_at(fq, &s);
            let value = qqdot(&qsub(a, &s), &m).sign();
            (s, value, m)
        }
        _ => return Ok(()),
    };
    if value == Ordering::Equal {
        return Ok(());
    }
    let gap = qsub(&s, a);
    let tol = q(res);
    if qqdot(&gap, &gap).add_r(&-(&tol * &tol)).sign() == Ordering::Greater {
        return Ok(());
    }
    if qm.in_face(fq, &s) == Loc::Out {
        return Ok(());
    }
    let Some(extremum) = apex_side(vm, &g, &zero()).map(Ordering::reverse) else {
        return Ok(());
    };
    if value == extremum {
        return Err(Error::Degenerate(NEAR_APEX));
    }
    let inside = [
        vm.member(a, std::slice::from_ref(&gap)) == Loc::In,
        qm.member(&s, &[gap.map(|x| x.neg())]) == Loc::In,
    ];
    if inside == [true, true] {
        return Ok(());
    }
    Err(Error::Degenerate(NEAR_APEX))
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
/// one is widened by the resolution is tried, after `vertex_faces`. Edges of
/// other curves (meetings of two curved faces, a cone's or a torus's
/// sections and meetings, a spline prism's cap edges) and points of two
/// fields are taken as REVIEW_NOTES.md's "Other curves, points of two fields
/// and an apex against a spline or torus face" decides: a corner `1e-12`
/// off a rod's hole's rim, two balls' zones' rims that near each other were
/// fused into two solids.
pub(super) fn edges_and_vertices(models: &[Prism; 2], res: f64) -> Result<()> {
    COARSE.with(|c| c.borrow_mut().clear());
    let out = edges_and_vertices_kept(models, res);
    COARSE.with(|c| c.borrow_mut().clear());
    out
}

/// A run's coarse samples kept through one `edges_and_vertices`
/// (`Run::coarse`).
type Coarse = std::collections::HashMap<(usize, u64, u64), std::rc::Rc<(Vec<[f64; 3]>, f64)>>;

thread_local! {
    static COARSE: std::cell::RefCell<Coarse> = std::cell::RefCell::new(Coarse::new());
}

/// `edges_and_vertices` with its runs' coarse samples kept.
fn edges_and_vertices_kept(models: &[Prism; 2], res: f64) -> Result<()> {
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
            let a = End::at(&models[0], sa, models[0].verts[va].p.clone());
            let b = End::at(&models[1], sb, models[1].verts[vb].p.clone());
            pair(&a, &b, &Bracket::Exact, res, NEAR_VERTICES)?;
        }
    }
    let wide = |b: &([f64; 3], [f64; 3])| (b.0.map(|x| x - 2.0 * res), b.1.map(|x| x + 2.0 * res));
    let edge_box = |m: &Prism, e: &MEdge| intersect(&m.boxes[e.faces[0]], &m.boxes[e.faces[1]]);
    // Vertex and edge, each input's vertices against the other's edges.
    let takes: [Vec<bool>; 2] = [0, 1].map(|o| {
        (0..models[o].edges.len())
            .map(|ei| takes_edge(&models[o], ei))
            .collect()
    });
    for o in 0..2 {
        let (vm, em) = (&models[o], &models[1 - o]);
        for (ei, e) in em.edges.iter().enumerate() {
            if !takes[1 - o][ei] {
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
                // The vertex's side of the pair (the first input's `A`).
                let side = if o == 0 { Side::A } else { Side::B };
                for (foot, bracket) in nearest_to_point(em, ei, p, side, res) {
                    let (vertex, point) = (
                        End::at(vm, spot, p.clone()),
                        End {
                            m: em,
                            spot: &Spot::Edge(ei),
                            p: foot.p,
                            exact: foot.exact,
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
        if !takes[0][ea] {
            continue;
        }
        let ba = wide(&edge_box(&models[0], e));
        for (eb, f) in models[1].edges.iter().enumerate() {
            if !takes[1][eb] || !boxes_meet(&ba, &edge_box(&models[1], f)) {
                continue;
            }
            for (x, y, bracket) in closest_on_edges(&models[0], ea, &models[1], eb, res) {
                let a = End {
                    m: &models[0],
                    spot: &Spot::Edge(ea),
                    p: x.p,
                    exact: x.exact,
                };
                let b = End {
                    m: &models[1],
                    spot: &Spot::Edge(eb),
                    p: y.p,
                    exact: y.exact,
                };
                pair(&a, &b, &bracket, res, NEAR_EDGES)?;
            }
        }
    }
    Ok(())
}

/// An input edge the rules take: with an input's id, not between faces on
/// one surface, its curve a line, a conic, a circle or another curve
/// (`Run`).
fn takes_edge(m: &Prism, ei: usize) -> bool {
    let e = &m.edges[ei];
    e.id.is_some()
        && !one_surface(m, e.faces[0], e.faces[1])
        && match e.curve {
            Crv::Line { .. } | Crv::Conic { .. } | Crv::Circle(_) => true,
            _ => Run::new(m, ei).is_some(),
        }
}

/// A point of a pair found on an edge: its own field's point, and for a
/// point of two fields (a line's foot from a point of another field) the
/// point itself (`End`).
struct Pt {
    p: QV,
    exact: Option<MixV>,
}

impl Pt {
    fn of(p: QV) -> Self {
        Self { p, exact: None }
    }
}

/// A point of an input's boundary for `pair`.
struct End<'a> {
    m: &'a Prism,
    spot: &'a Spot<'a>,
    /// The point; for a point of two fields (a line's foot from a point of
    /// another field, `exact`), a point of the line's own field beside it,
    /// where its directions and its membership are taken (a line's faces'
    /// directions are the same all along it).
    p: QV,
    exact: Option<MixV>,
}

impl<'a> End<'a> {
    fn at(m: &'a Prism, spot: &'a Spot<'a>, p: QV) -> Self {
        Self {
            m,
            spot,
            p,
            exact: None,
        }
    }
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
        size.log2_rn().ceil() as i32
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
    if a.exact.is_some() || b.exact.is_some() || !one_field(&[&a.p, &b.p]) {
        // Points of two fields: one where the surds' ratio is a rational
        // square, else sums over both (`pair_mixed`).
        if a.exact.is_none() && b.exact.is_none() {
            if let Some(p) = rebased(&b.p, &a.p) {
                let b = End::at(b.m, b.spot, p);
                return pair(a, &b, bracket, res, why);
            }
        }
        return pair_mixed(a, b, bracket, res, why);
    }
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

// ---------------------------------------------------- points of two fields

/// Which point of a pair a number belongs with: the first input's (`A`) or
/// the second's (`B`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Side {
    A,
    B,
}

impl Side {
    fn other(self) -> Self {
        match self {
            Side::A => Side::B,
            Side::B => Side::A,
        }
    }
}

/// A number of a pair of points of two fields (REVIEW_NOTES.md, "Other
/// curves, points of two fields and an apex against a spline or torus
/// face"): `a + b + sum a_i b_i`, `a` and each `a_i` of the first point's
/// field, `b` and each `b_i` of the second's. Sums and products of such
/// numbers are such numbers; the sign is `mixed_dot_sign`'s (exact over
/// the rationals and one algebraic base, by enclosures over two).
#[derive(Clone, Debug)]
struct Mix {
    a: Qd,
    b: Qd,
    terms: Vec<(Qd, Qd)>,
}

/// A point or a vector of `Mix`.
type MixV = [Mix; 3];

impl Mix {
    fn zero() -> Self {
        Self {
            a: Qd::rat(zero()),
            b: Qd::rat(zero()),
            terms: Vec::new(),
        }
    }

    fn r(x: &R) -> Self {
        Self::of(&Qd::rat(x.clone()), Side::A)
    }

    /// A number of one side's field.
    fn of(x: &Qd, side: Side) -> Self {
        let mut m = Self::zero();
        match side {
            Side::A => m.a = x.clone(),
            Side::B => m.b = x.clone(),
        }
        m
    }

    /// Adds the product of `x` (the first field's) and `y` (the second's),
    /// into a side's own part where either is rational, beside a term of the
    /// same `y` where there is one.
    fn push(&mut self, x: Qd, y: Qd) {
        if let Some(k) = y.rational() {
            self.a = self.a.add(&x.scale(k));
        } else if let Some(k) = x.rational() {
            self.b = self.b.add(&y.scale(k));
        } else if let Some(t) = self.terms.iter_mut().find(|t| t.1 == y) {
            t.0 = t.0.add(&x);
        } else {
            self.terms.push((x, y));
        }
    }

    fn add(&self, o: &Self) -> Self {
        let mut m = Self {
            a: self.a.add(&o.a),
            b: self.b.add(&o.b),
            terms: self.terms.clone(),
        };
        for (x, y) in &o.terms {
            m.push(x.clone(), y.clone());
        }
        m
    }

    fn neg(&self) -> Self {
        Self {
            a: self.a.neg(),
            b: self.b.neg(),
            terms: self
                .terms
                .iter()
                .map(|(x, y)| (x.neg(), y.clone()))
                .collect(),
        }
    }

    fn sub(&self, o: &Self) -> Self {
        self.add(&o.neg())
    }

    fn scale(&self, k: &R) -> Self {
        Self {
            a: self.a.scale(k),
            b: self.b.scale(k),
            terms: self
                .terms
                .iter()
                .map(|(x, y)| (x.scale(k), y.clone()))
                .collect(),
        }
    }

    fn mul(&self, o: &Self) -> Self {
        let mut m = Self {
            a: self.a.mul(&o.a),
            b: self.b.mul(&o.b),
            terms: Vec::new(),
        };
        m.push(self.a.clone(), o.b.clone());
        m.push(o.a.clone(), self.b.clone());
        for (x, y) in &o.terms {
            m.push(self.a.mul(x), y.clone());
            m.push(x.clone(), self.b.mul(y));
        }
        for (x, y) in &self.terms {
            m.push(x.mul(&o.a), y.clone());
            m.push(x.clone(), y.mul(&o.b));
            for (u, v) in &o.terms {
                m.push(x.mul(u), y.mul(v));
            }
        }
        m
    }

    /// A binary64 interval enclosure (`None` out of binary64's range).
    fn fast(&self) -> Option<Fast> {
        let mut acc = self
            .a
            .enclose_fast(SIGN_STEPS)?
            .add(&self.b.enclose_fast(SIGN_STEPS)?);
        for (x, y) in &self.terms {
            acc = acc.add(
                &x.enclose_fast(SIGN_STEPS)?
                    .mul(&y.enclose_fast(SIGN_STEPS)?),
            );
        }
        Some(acc)
    }

    fn sign(&self) -> Ordering {
        // Certain where the enclosure excludes zero.
        if let Some(s @ (Ordering::Less | Ordering::Greater)) = self.fast().and_then(|v| v.sign()) {
            return s;
        }
        let mut a = vec![self.a.clone(), Qd::rat(int(1))];
        let mut b = vec![Qd::rat(int(1)), self.b.clone()];
        for (x, y) in &self.terms {
            a.push(x.clone());
            b.push(y.clone());
        }
        mixed_dot_sign(&a, &b)
    }

    /// A rational enclosure.
    fn interval(&self) -> I {
        self.terms
            .iter()
            .fold(self.a.interval().add(&self.b.interval()), |acc, (x, y)| {
                acc.add(&x.interval().mul(&y.interval()))
            })
    }

    /// A rational within the enclosure's width of the number.
    fn near(&self) -> R {
        let i = self.interval();
        (i.lo() + i.hi()) / int(2)
    }
}

fn mixv(p: &QV, side: Side) -> MixV {
    [0, 1, 2].map(|k| Mix::of(&p[k], side))
}

fn mv_sub(a: &MixV, b: &MixV) -> MixV {
    [0, 1, 2].map(|k| a[k].sub(&b[k]))
}

fn mv_dot(a: &MixV, b: &MixV) -> Mix {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

/// `n . w`, `w` of one side's field.
fn mv_dot_q(n: &MixV, w: &QV, side: Side) -> Mix {
    mv_dot(n, &mixv(w, side))
}

/// `n . v`, `v` rational.
fn mv_dot_r(n: &MixV, v: &V) -> Mix {
    n[0].scale(&v[0])
        .add(&n[1].scale(&v[1]))
        .add(&n[2].scale(&v[2]))
}

/// A point of `Q(sqrt E)` written in `Q(sqrt D)`, the field of `like`'s
/// numbers, both over the rationals, where `E / D` is a rational square
/// (`sqrt E = r sqrt D`): one field. `None` otherwise.
fn rebased(p: &QV, like: &QV) -> Option<QV> {
    if p.iter().chain(like.iter()).any(|x| x.gen().is_some()) {
        return None;
    }
    let d = like.iter().find_map(|x| x.field())?.clone();
    let e = p.iter().find_map(|x| x.field())?.clone();
    let r = rational_sqrt(&(&e / &d))?;
    let mut out = Vec::new();
    for x in p {
        let (K::Rat(a), K::Rat(b)) = (&x.a, &x.b) else {
            return None;
        };
        out.push(if x.field().is_some() {
            Qd::new(a.clone(), b * &r, d.clone())
        } else {
            x.clone()
        });
    }
    let out: QV = out.try_into().ok()?;
    one_field(&[&out, like]).then_some(out)
}

/// A binary64 interval enclosure of a point's coordinates.
fn fast_point(p: &QV) -> Option<[Fast; 3]> {
    Some([
        p[0].enclose_fast(SIGN_STEPS)?,
        p[1].enclose_fast(SIGN_STEPS)?,
        p[2].enclose_fast(SIGN_STEPS)?,
    ])
}

/// The sign of `(x - y) . t`, `t` of `x`'s field: in one field where `x`
/// and `y` share it, else over both.
fn rel_sign(x: &QV, y: &QV, t: &QV) -> Ordering {
    if one_field(&[x, y]) {
        return qqdot(&qsub(x, y), t).sign();
    }
    // Binary64 enclosures first: certain where they exclude zero.
    let f = |p: &QV| -> Option<[Fast; 3]> {
        Some([
            p[0].enclose_fast(SIGN_STEPS)?,
            p[1].enclose_fast(SIGN_STEPS)?,
            p[2].enclose_fast(SIGN_STEPS)?,
        ])
    };
    if let (Some(fx), Some(fy), Some(ft)) = (f(x), f(y), f(t)) {
        let s = (0..3).fold(Fast::exact_f64(0.0), |acc, k| {
            acc.add(&fx[k].sub(&fy[k]).mul(&ft[k]))
        });
        if let Some(o @ (Ordering::Less | Ordering::Greater)) = s.sign() {
            return o;
        }
    }
    mv_dot_q(&mv_sub(&mixv(x, Side::A), &mixv(y, Side::B)), t, Side::A).sign()
}

/// The squared gap's binary64 view and the gap's own rounding (`rho^2`,
/// `gap_rounding`'s, from the points' binary64 views) of two points of a
/// pair in two fields.
fn mixed_rounding(a: &MixV, b: &MixV) -> R {
    let view = |p: &MixV| p.clone().map(|x| rational_f64(&x.near()));
    let size = view(a)
        .iter()
        .chain(view(b).iter())
        .fold(1.0f64, |m, x| m.max(x.abs()));
    let e = if size.is_finite() {
        size.log2_rn().ceil() as i32
    } else {
        0
    };
    R::from_integer(BigInt::from(2)).pow(2 * (e - 48))
}

/// `pair` for points of two fields: the same rule, each number a `Mix`
/// (`gap_within`'s gap and bound, `contact`'s directions and band, an
/// apex's rulings). A missing pair's gap inside both inputs is each one's
/// membership at its point along a direction of its own field on the gap's
/// side of every tangent plane of its faces there (`member_mixed`); where
/// that is undecided (the gap along such a plane), `ComputationLimit`.
fn pair_mixed(a: &End, b: &End, bracket: &Bracket, res: f64, why: &'static str) -> Result<()> {
    let xa = a.exact.clone().unwrap_or_else(|| mixv(&a.p, Side::A));
    let xb = b.exact.clone().unwrap_or_else(|| mixv(&b.p, Side::B));
    let n = mv_sub(&xb, &xa);
    let nn = mv_dot(&n, &n);
    if nn.sign() == Ordering::Equal {
        return Ok(());
    }
    let tol = q(res);
    let above = |reach: &R| nn.sub(&Mix::r(&(reach * reach))).sign() == Ordering::Greater;
    let within = !above(&tol);
    let unsure = !within
        && match bracket {
            Bracket::Exact => false,
            Bracket::Within(d) => !above(&(&tol + d)),
            Bracket::Unbracketed => !above(&(int(2) * &tol)),
        };
    if !(within || unsure) {
        return Ok(());
    }
    let rho2 = mixed_rounding(&xa, &xb);
    let sa = side_mixed(a, Side::A, &n, &nn, &rho2);
    let sb = side_mixed(b, Side::B, &n, &nn, &rho2);
    let missing = match (sa, sb) {
        (Some(Ordering::Less), Some(Ordering::Greater)) => true,
        (Some(Ordering::Greater), Some(Ordering::Less)) => false,
        _ => return Ok(()),
    };
    if missing {
        let back = n.clone().map(|x| x.neg());
        let ia = member_mixed(a, Side::A, &n);
        let ib = member_mixed(b, Side::B, &back);
        let out = |l: &Option<Loc>| matches!(l, Some(x) if *x != Loc::In);
        if !(out(&ia) || out(&ib)) {
            if ia.is_none() || ib.is_none() {
                return Err(Error::ComputationLimit(NEAR_FIELDS));
            }
            return Ok(());
        }
    }
    if unsure {
        return Err(Error::ComputationLimit(NEAR_LIMIT));
    }
    Err(Error::Degenerate(why))
}

/// `side` for a point of a pair in two fields against the gap `n` (a
/// `Mix`): its directions (of its own field, on `side`) each strictly of
/// one sign outside the band `(n . w)^2 <= (10^-24 |n|^2 + rho^2) |w|^2`.
fn side_mixed(e: &End, side: Side, n: &MixV, nn: &Mix, rho2: &R) -> Option<Ordering> {
    let dirs = match e.spot {
        Spot::Apex => return apex_side_mixed(e.m, n, nn, rho2),
        Spot::Vertex(ends) => corner_dirs(e.m, ends, &e.p)?,
        Spot::Edge(ei) => edge_dirs(e.m, *ei, &e.p).to_vec(),
    };
    let band = nn.scale(&(int(1) / int(10).pow(24))).add(&Mix::r(rho2));
    let mut sign = None;
    for w in &dirs {
        let d = mv_dot_q(n, w, side);
        let ww = Mix::of(&qqdot(w, w), side);
        if d.mul(&d).sub(&band.mul(&ww)).sign() != Ordering::Greater {
            return None;
        }
        let s = d.sign();
        if sign.is_some_and(|x| x != s) {
            return None;
        }
        sign = Some(s);
    }
    sign
}

/// `apex_side` against a gap `n` of two fields.
fn apex_side_mixed(m: &Prism, n: &MixV, nn: &Mix, rho2: &R) -> Option<Ordering> {
    let fun = m.funnel.as_ref()?;
    fun.apex?;
    let f = &m.f;
    let (ca, cb, cn) = (mv_dot_r(n, &f.x), mv_dot_r(n, &f.y), mv_dot_r(n, &f.n));
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
    let e2 = nn
        .scale(&(int(1) / int(10).pow(24)))
        .add(&Mix::r(rho2))
        .scale(&(int(1) + &kk));
    let lhs = c2.sub(&e2).sub(&p2);
    if lhs.sign() != Ordering::Greater
        || lhs.mul(&lhs).sub(&e2.mul(&p2).scale(&int(4))).sign() != Ordering::Greater
    {
        return None;
    }
    Some(c.sign())
}

/// The membership of a pair's point pushed along a gap `n` of two fields:
/// `member` at the point (of its own field) along a rational direction
/// within the enclosures' width of `n`, taken where it lies on `n`'s side
/// of the tangent plane of every face at the point (the material's tangent
/// cone there is bounded by those planes, so both have one membership); a
/// cone's apex's by its rulings exactly (`n` strictly inside the cone,
/// `A^2 + B^2 < k^2 C^2` with `s C > 0`). `None` where undecided: `n` along
/// such a plane, or on the cone.
fn member_mixed(e: &End, side: Side, n: &MixV) -> Option<Loc> {
    let faces: Vec<usize> = match e.spot {
        Spot::Apex => {
            let fun = e.m.funnel.as_ref()?;
            let f = &e.m.f;
            let (ca, cb, cn) = (mv_dot_r(n, &f.x), mv_dot_r(n, &f.y), mv_dot_r(n, &f.n));
            let k = &fun.k;
            let c = if sign(k) == Ordering::Greater {
                cn
            } else {
                cn.neg()
            };
            let q2 = c
                .mul(&c)
                .scale(&(k * k))
                .sub(&ca.mul(&ca).add(&cb.mul(&cb)));
            return match (c.sign(), q2.sign()) {
                (Ordering::Greater, Ordering::Greater) => Some(Loc::In),
                (_, Ordering::Less) | (Ordering::Less, _) => Some(Loc::Out),
                _ => None,
            };
        }
        Spot::Vertex(ends) => {
            let mut fs: Vec<usize> = ends
                .iter()
                .flat_map(|&(ei, _)| e.m.edges[ei].faces)
                .collect();
            fs.sort_unstable();
            fs.dedup();
            fs
        }
        Spot::Edge(ei) => e.m.edges[*ei].faces.to_vec(),
    };
    let approx: V = n.clone().map(|x| x.near());
    let dir = qv(&approx);
    for f in faces {
        let m = e.m.normal_at(f, &e.p);
        let exact = mv_dot_q(n, &m, side).sign();
        if exact == Ordering::Equal || qqdot(&dir, &m).sign() != exact {
            return None;
        }
    }
    Some(e.m.member(&e.p, &[dir]))
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

/// The points of edge `ei` of `m` nearest a point `p` (on the pair's
/// `side`) inside the edge (`inside_edge`): a line's foot (of two fields
/// where the line's points and `p` lie in two, `line_foot`), a circle's
/// point toward `p` (exact, in the quadratic field of the projection's
/// length, where `p` is rational), on a conic that is no circle exactly, or
/// a circle against a point of another field, each least of the distance
/// bracketed (`bracketed`), and on another curve each least found in
/// binary64, its exact point decided there or bracketed (`run_nearest`).
fn nearest_to_point(m: &Prism, ei: usize, p: &QV, side: Side, res: f64) -> Vec<(Pt, Bracket)> {
    let e = &m.edges[ei];
    let target = qv_f64(p);
    let far = |x: [f64; 3]| (0..3).map(|k| (x[k] - target[k]).powi(2)).sum::<f64>();
    let slope = |x: &QV, t: &QV| rel_sign(x, p, t);
    let found: Vec<(Pt, Bracket)> = match &e.curve {
        Crv::Line { p: o, d } => {
            if !one_field(&[p, o]) {
                // Beyond twice the resolution in binary64 by far (its error
                // a few ulps of the points' size): none.
                let (fo, fd) = (qv_f64(o), d.clone().map(|x| rational_f64(&x)));
                let w = sub64(&target, &fo);
                let c = [
                    w[1] * fd[2] - w[2] * fd[1],
                    w[2] * fd[0] - w[0] * fd[2],
                    w[0] * fd[1] - w[1] * fd[0],
                ];
                let size = dot64(&target, &target).sqrt() + dot64(&fo, &fo).sqrt() + 1.0;
                if (dot64(&c, &c) / dot64(&fd, &fd)).sqrt() > 4.0 * res + 1e-12 * size {
                    return Vec::new();
                }
                return line_foot(m, ei, p, side)
                    .map(|(_, foot)| (foot, Bracket::Exact))
                    .into_iter()
                    .collect();
            }
            let t = qdot(&qsub(p, o), d).scale(&(int(1) / dot(d, d)));
            vec![(Pt::of(qadd(o, &qscale(d, &t))), Bracket::Exact)]
        }
        Crv::Circle(circ) => match rational_point(p) {
            Some(c) => super::graph::circle_nearest(&circ.c, &circ.x, &circ.y, &circ.r2, &c)
                .map(|x| (Pt::of(x), Bracket::Exact))
                .into_iter()
                .collect(),
            None => bracketed(&e.curve, &far, &slope, res)
                .into_iter()
                .map(|(x, k)| (Pt::of(x), k))
                .collect(),
        },
        Crv::Conic { c, a, b } => {
            let circle = dot(a, b) == zero() && dot(a, a) == dot(b, b);
            match (circle, rational_point(p)) {
                (true, Some(x)) => super::graph::circle_nearest(c, a, b, &dot(a, a), &x)
                    .map(|x| (Pt::of(x), Bracket::Exact))
                    .into_iter()
                    .collect(),
                _ => bracketed(&e.curve, &far, &slope, res)
                    .into_iter()
                    .map(|(x, k)| (Pt::of(x), k))
                    .collect(),
            }
        }
        _ => match Run::new(m, ei) {
            Some(run) => {
                let pf = fast_point(p);
                let far_fast = |c: &[Fast; 3]| {
                    let pf = pf.as_ref()?;
                    Some((0..3).fold(Fast::exact_f64(0.0), |acc, k| {
                        acc.add(&c[k].sub(&pf[k]).square())
                    }))
                };
                let other = |_: &QV| Some(Pt::of(p.clone()));
                run_nearest(&run, &far, &far_fast, &other, side.other(), res)
                    .into_iter()
                    .map(|(x, _, k)| (Pt::of(x), k))
                    .collect()
            }
            None => Vec::new(),
        },
    };
    found
        .into_iter()
        .filter(|(x, _)| x.exact.is_some() || inside_edge(m, ei, &x.p))
        .collect()
}

/// The foot on line edge `ei` of `m` of a point `x` of another field (on
/// the pair's `side`, the line's point on the other): its parameter `t =
/// (x - o) . d / |d|^2` and the foot `o + t d` (`Mix`), strictly inside the
/// edge by exact comparisons of `t` with its ends' parameters, with a point
/// of the line's own field beside it (`o + t' d`, `t'` a rational within the
/// enclosures' width of `t`, inside the edge too). `None` off the edge, at
/// an end, or where an end lies in a third field.
fn line_foot(m: &Prism, ei: usize, x: &QV, side: Side) -> Option<(Mix, Pt)> {
    let Crv::Line { p: o, d } = &m.edges[ei].curve else {
        return None;
    };
    let (t, foot) = foot_mix(x, side, o, d);
    let proxy = line_inside(m, ei, &t, side.other())?;
    Some((
        t,
        Pt {
            p: proxy,
            exact: Some(foot),
        },
    ))
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
) -> Vec<(Pt, Pt, Bracket)> {
    let (ca, cb) = (&ma.edges[ea].curve, &mb.edges[eb].curve);
    let both = |v: Vec<(QV, QV, Bracket)>| -> Vec<(Pt, Pt, Bracket)> {
        v.into_iter()
            .map(|(x, y, k)| (Pt::of(x), Pt::of(y), k))
            .collect()
    };
    let found: Vec<(Pt, Pt, Bracket)> = match (ca, cb) {
        (Crv::Line { p: p1, d: d1 }, Crv::Line { p: p2, d: d2 }) => {
            if lines_apart(p1, d1, p2, d2, res) {
                return Vec::new();
            }
            if !one_field(&[p1, p2]) {
                return lines_mixed(ma, ea, mb, eb)
                    .into_iter()
                    .filter(|(x, y, k)| takes_gap(x, y, k, res))
                    .collect();
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
            both(vec![(
                qadd(p1, &qscale(d1, &s)),
                qadd(p2, &qscale(d2, &t)),
                Bracket::Exact,
            )])
        }
        (Crv::Line { .. }, Crv::Conic { .. } | Crv::Circle(_)) => {
            line_and_curve((mb, eb), (ma, ea), Side::A, res)
                .into_iter()
                .map(|(on_curve, foot, k)| (foot, Pt::of(on_curve), k))
                .collect()
        }
        (Crv::Conic { .. } | Crv::Circle(_), Crv::Line { .. }) => {
            line_and_curve((ma, ea), (mb, eb), Side::B, res)
                .into_iter()
                .map(|(on_curve, foot, k)| (Pt::of(on_curve), foot, k))
                .collect()
        }
        (Crv::Conic { .. } | Crv::Circle(_), Crv::Conic { .. } | Crv::Circle(_)) => {
            both(conics_closest(ma, ea, mb, eb, res))
        }
        _ => run_pairs(ma, ea, mb, eb, res),
    };
    // Only gaps `pair` takes (`gap_within`), told before the points' places
    // along the edges (which cost more); a point of two fields' place is
    // decided where it is found.
    found
        .into_iter()
        .filter(|(x, y, k)| {
            takes_gap(x, y, k, res)
                && (x.exact.is_some() || inside_edge(ma, ea, &x.p))
                && (y.exact.is_some() || inside_edge(mb, eb, &y.p))
        })
        .collect()
}

/// Whether `pair` would take a pair's points (`gap_within`; over two fields
/// alike).
fn takes_gap(x: &Pt, y: &Pt, bracket: &Bracket, res: f64) -> bool {
    if x.exact.is_none() && y.exact.is_none() && one_field(&[&x.p, &y.p]) {
        return gap_within(&x.p, &y.p, bracket, res).is_some();
    }
    let xa = x.exact.clone().unwrap_or_else(|| mixv(&x.p, Side::A));
    let xb = y.exact.clone().unwrap_or_else(|| mixv(&y.p, Side::B));
    let n = mv_sub(&xb, &xa);
    let nn = mv_dot(&n, &n);
    if nn.sign() == Ordering::Equal {
        return false;
    }
    let tol = q(res);
    let reach = match bracket {
        Bracket::Exact => tol,
        Bracket::Within(d) => &tol + d,
        Bracket::Unbracketed => int(2) * &tol,
    };
    nn.sub(&Mix::r(&(&reach * &reach))).sign() != Ordering::Greater
}

/// The closest points of two skew line edges whose points lie in two fields
/// (`p1` the first's, `p2` the second's): `s` and `t` of `Mix`, each
/// strictly inside its edge (`line_inside`), with points of each line's own
/// field beside them.
fn lines_mixed(ma: &Prism, ea: usize, mb: &Prism, eb: usize) -> Option<(Pt, Pt, Bracket)> {
    let (Crv::Line { p: p1, d: d1 }, Crv::Line { p: p2, d: d2 }) =
        (&ma.edges[ea].curve, &mb.edges[eb].curve)
    else {
        return None;
    };
    let (a, b, c) = (dot(d1, d1), dot(d1, d2), dot(d2, d2));
    let den = &a * &c - &b * &b;
    if den == zero() {
        return None;
    }
    let (o1, o2) = (mixv(p1, Side::A), mixv(p2, Side::B));
    let w = mv_sub(&o1, &o2);
    let (d, e) = (mv_dot_r(&w, d1), mv_dot_r(&w, d2));
    let inv = int(1) / &den;
    let s = e.scale(&b).sub(&d.scale(&c)).scale(&inv);
    let t = e.scale(&a).sub(&d.scale(&b)).scale(&inv);
    let x: MixV = [0, 1, 2].map(|k| o1[k].add(&s.scale(&d1[k])));
    let y: MixV = [0, 1, 2].map(|k| o2[k].add(&t.scale(&d2[k])));
    let px = line_inside(ma, ea, &s, Side::A)?;
    let py = line_inside(mb, eb, &t, Side::B)?;
    Some((
        Pt {
            p: px,
            exact: Some(x),
        },
        Pt {
            p: py,
            exact: Some(y),
        },
        Bracket::Exact,
    ))
}

/// Whether the point `o + t d` of line edge `ei` of `m` (its parameter `t`
/// a `Mix`, the line's point on `side`) lies strictly inside the edge: `t`
/// strictly between its ends' parameters exactly. The line's own field's
/// point beside it (`t'` a rational within the enclosures' width of `t`,
/// inside the edge too); `None` off the edge, at an end, or where an end
/// lies in another field than the line's point.
fn line_inside(m: &Prism, ei: usize, t: &Mix, side: Side) -> Option<QV> {
    let e = &m.edges[ei];
    let Crv::Line { p: o, d } = &e.curve else {
        return None;
    };
    let dd = dot(d, d);
    let mut ends = Vec::new();
    for v in [e.start, e.end] {
        let p = &m.verts[v].p;
        if !one_field(&[p, o]) {
            return None;
        }
        let u = qdot(&qsub(p, o), d).scale(&(int(1) / &dd));
        ends.push(t.sub(&Mix::of(&u, side)).sign());
    }
    if ends[0] == Ordering::Equal || ends[1] == Ordering::Equal || ends[0] == ends[1] {
        return None;
    }
    let proxy = qadd(o, &qscale(d, &Qd::rat(t.near())));
    inside_edge(m, ei, &proxy).then_some(proxy)
}

/// The foot `o + t d` of a point `x` (on `side`) on the line `o + t d` (its
/// point `o` on the other side), of `Mix`, and its parameter.
fn foot_mix(x: &QV, side: Side, o: &QV, d: &V) -> (Mix, MixV) {
    let line = side.other();
    let o = mixv(o, line);
    let t = mv_dot_r(&mv_sub(&mixv(x, side), &o), d).scale(&(int(1) / dot(d, d)));
    let foot = [0, 1, 2].map(|k| o[k].add(&t.scale(&d[k])));
    (t, foot)
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

/// The points of a conic's or a circle's edge `(mc, ec)` at the leasts of
/// its distance from a line edge `(ml, el)` (the line on the pair's
/// `line_side`), each with its exact foot on the line: of one field, or of
/// two where a circle's points and the line's lie in two (`line_foot`).
fn line_and_curve(
    (mc, ec): (&Prism, usize),
    (ml, el): (&Prism, usize),
    line_side: Side,
    res: f64,
) -> Vec<(QV, Pt, Bracket)> {
    let (crv, line) = (&mc.edges[ec].curve, &ml.edges[el].curve);
    let Crv::Line { p, d } = line else {
        return Vec::new();
    };
    let side = line_side.other();
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
    let slope = |x: &QV, t: &QV| {
        if one_field(&[x, p]) {
            return qqdot(&qsub(x, &foot(x)), t).sign();
        }
        let (_, f) = foot_mix(x, side, p, d);
        mv_dot_q(&mv_sub(&mixv(x, side), &f), t, side).sign()
    };
    bracketed(crv, &far, &slope, res)
        .into_iter()
        .filter_map(|(x, k)| {
            let y = if one_field(&[&x, p]) {
                Pt::of(foot(&x))
            } else {
                line_foot(ml, el, &x, side)?.1
            };
            Some((x, y, k))
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
        None => [0, 1, 2].map(|k| fc[k] + fa[k] * t.cos_rn() + fb[k] * t.sin_rn()),
        Some(r) => {
            let dir = [0, 1, 2].map(|k| fa[k] * t.cos_rn() + fb[k] * t.sin_rn());
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
        let flip = t.cos_rn() < 0.0;
        let s = (if flip { t - PI } else { t } / 2.0).tan_rn();
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
        let (cs, sn) = (t.cos_rn(), t.sin_rn());
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
            Some((lx, ly)) => (t.sin_rn() / ly).atan2_rn(t.cos_rn() / lx),
        };
        let flip = t.cos_rn() < 0.0;
        let u = if flip { t - PI } else { t };
        ((u.sin_rn().atan2_rn(u.cos_rn()) / 2.0).tan_rn(), flip)
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

/// The sign of `|y - x|^2 - r2` for two points of any fields.
fn sq_cmp(x: &QV, y: &QV, r2: &R) -> Ordering {
    if one_field(&[x, y]) {
        let n = qsub(y, x);
        return qqdot(&n, &n).add_r(&-r2.clone()).sign();
    }
    let n = mv_sub(&mixv(y, Side::B), &mixv(x, Side::A));
    mv_dot(&n, &n).sub(&Mix::r(r2)).sign()
}

/// The gap's own rounding (`rho^2`) of two points of any fields.
fn rounding_of(x: &QV, y: &QV) -> R {
    if one_field(&[x, y]) {
        gap_rounding(x, y)
    } else {
        mixed_rounding(&mixv(x, Side::A), &mixv(y, Side::B))
    }
}

/// `contact` for two points of any fields (the first input's first).
fn contact_any(a: &End, b: &End) -> Option<bool> {
    if one_field(&[&a.p, &b.p]) {
        return contact(a, b, &qsub(&b.p, &a.p));
    }
    let (xa, xb) = (mixv(&a.p, Side::A), mixv(&b.p, Side::B));
    let n = mv_sub(&xb, &xa);
    let nn = mv_dot(&n, &n);
    let rho2 = mixed_rounding(&xa, &xb);
    match (
        side_mixed(a, Side::A, &n, &nn, &rho2),
        side_mixed(b, Side::B, &n, &nn, &rho2),
    ) {
        (Some(Ordering::Less), Some(Ordering::Greater)) => Some(true),
        (Some(Ordering::Greater), Some(Ordering::Less)) => Some(false),
        _ => None,
    }
}

/// The closest points of two conic or circle edges of the inputs
/// (REVIEW_NOTES.md's "Two conic edges within the resolution of each
/// other"; two `Crv::Circle`s too, "Other curves, points of two fields
/// ...": their points lie in fields of their own): at each least found in
/// binary64 (`conic_leasts`) the curves' exact
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
    let (Some(fa), Some(fb)) = (Conic64::of(ca), Conic64::of(cb)) else {
        return Vec::new();
    };
    let tol = q(res);
    let mut out = Vec::new();
    for (t, s) in conic_leasts(&fa, &fb, res) {
        let ((ua, fla), (ub, flb)) = (fa.param(t), fb.param(s));
        // Of rounding's size beyond doubt in binary64 (coincident rims: a
        // least at every node along them): no contact, no exact point.
        if rounding_gap(&fa, (ua, fla), &fb, (ub, flb)) {
            continue;
        }
        let (Some((x, _)), Some((y, _))) = (curve_at(ca, &q(ua), fla), curve_at(cb, &q(ub), flb))
        else {
            continue;
        };
        if sq_cmp(&x, &y, &rounding_of(&x, &y)) != Ordering::Greater {
            continue;
        }
        if sq_cmp(&x, &y, &(&tol * &tol)) != Ordering::Greater {
            out.push((x, y, Bracket::Exact));
            continue;
        }
        let (sa, sb) = (Spot::Edge(ea), Spot::Edge(eb));
        let (a, b) = (End::at(ma, &sa, x.clone()), End::at(mb, &sb, y.clone()));
        if contact_any(&a, &b).is_none() || !inside_edge(ma, ea, &x) || !inside_edge(mb, eb, &y) {
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

/// Whether binary64 shows the gap between two conics' exact points at the
/// rational half-angle tangents `(u, flip)` of each (`curve_at`'s) to be of
/// rounding's size, `|n| <= rho` (`gap_rounding`; `conics_closest` then
/// passes over the pair, every direction along the gap's plane). Each point
/// is evaluated in binary64 from its curve's data rounded once and the same
/// rational parameter: `c + a cs + b sn`, `cs = (1 - u^2) / (1 + u^2)`
/// within `4 eps` and `sn = 2 u / (1 + u^2)` within `3 eps` (`|u| <= 1`,
/// `eps = 2^-53`), so each coordinate within `(2 |c| + 8 |a| + 7 |b|) eps`
/// of the exact point's to first order (the data's rounding, two products,
/// two sums), bounded here by `10 (|c| + |a| + |b|) eps`; the gap's bound is
/// the binary64 points' distance (enlarged for its own rounding) plus both
/// points' bounds, against `rho` of a size no larger than the exact points'
/// (the largest coordinate less twice the bound). A circle's points
/// (`Crv::Circle`, scaled by a square root) are not taken.
fn rounding_gap(fa: &Conic64, pa: (f64, bool), fb: &Conic64, pb: (f64, bool)) -> bool {
    const EPS: f64 = 1.0 / (1u64 << 53) as f64;
    let point = |f: &Conic64, (u, flip): (f64, bool)| {
        let w = u * u;
        let (cs, sn) = ((1.0 - w) / (1.0 + w), 2.0 * u / (1.0 + w));
        let (cs, sn) = if flip { (-cs, -sn) } else { (cs, sn) };
        let x = [0, 1, 2].map(|k| f.c[k] + (f.a[k] * cs + f.b[k] * sn));
        let e = [0, 1, 2].map(|k| 10.0 * EPS * (f.c[k].abs() + f.a[k].abs() + f.b[k].abs()));
        (x, dot64(&e, &e).sqrt() * (1.0 + 1e-6))
    };
    if fa.basis.is_some() || fb.basis.is_some() || pa.0.abs() > 1.0 || pb.0.abs() > 1.0 {
        return false;
    }
    let ((x, ex), (y, ey)) = (point(fa, pa), point(fb, pb));
    let w = sub64(&y, &x);
    let gap = dot64(&w, &w).sqrt() * (1.0 + 16.0 * EPS) + ex + ey;
    let size = x.iter().chain(y.iter()).fold(0.0f64, |m, v| m.max(v.abs()));
    let low = ((size - 2.0 * ex.max(ey)) * (1.0 - 16.0 * EPS)).max(1.0);
    if !(gap.is_finite() && low.is_finite()) {
        return false;
    }
    let rho = 2f64.powi(low.log2_rn().ceil() as i32 - 48);
    gap <= rho
}

/// The inner curve of `conics_bracketed`: a circle exactly (its centre,
/// basis and squared radius, and the curve for a point `x` not rational),
/// or a conic bracketed about a binary64 angle (with its speed bound).
enum Inner<'a> {
    Round(&'a V, &'a V, &'a V, R, &'a Crv, &'a Conic64, f64),
    Conic(&'a Crv, &'a Conic64, f64, R),
}

impl Inner<'_> {
    /// The inner curve's point nearest the point `x` and a bound of its
    /// distance from the nearest's: exactly on a circle (`circle_nearest`)
    /// where `x` is rational, else bracketed on exact signs of `(y - x) . t`
    /// (`bisect`, over two fields where `x` lies in another) about the
    /// binary64 nearest angle.
    fn nearest(&self, x: &QV) -> Option<(QV, R)> {
        let walk = |crv: &Crv, f: &Conic64, s: f64, speed: &R| {
            let s = f.nearest(&qv_f64(x), s);
            let (u, flip) = f.param(s);
            let slope = |y: &QV, t: &QV| rel_sign(y, x, t);
            match bisect(crv, u, flip, &slope, speed)? {
                (y, Bracket::Within(r)) => Some((y, r)),
                _ => None,
            }
        };
        match self {
            Inner::Round(c, a, b, r2, crv, f, s) => match rational_point(x) {
                Some(xr) => super::graph::circle_nearest(c, a, b, r2, &xr).map(|y| (y, zero())),
                None => walk(crv, f, *s, &speed(crv)?),
            },
            Inner::Conic(crv, f, s, speed) => walk(crv, f, *s, speed),
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
        Crv::Circle(circ) => Inner::Round(&circ.c, &circ.x, &circ.y, circ.r2.clone(), ci, fi, si),
        Crv::Conic { c, a, b } if round(ci) => Inner::Round(c, a, b, dot(a, a), ci, fi, si),
        _ => Inner::Conic(ci, fi, si, speed(ci)?),
    };
    let slope = |x: &QV, ta: &QV| {
        let Some((y, r)) = inner.nearest(x) else {
            return Ordering::Equal;
        };
        beyond_sign(x, &y, ta, &r)
    };
    let (u, flip) = fo.param(to);
    let (x, Bracket::Within(da)) = bisect(co, u, flip, &slope, &speed_o)? else {
        return None;
    };
    let (y, r) = inner.nearest(&x)?;
    let bound = Bracket::Within(da + r);
    Some(if swap { (y, x, bound) } else { (x, y, bound) })
}

/// The sign of `phi = (x - y) . t` (`t` of `x`'s field) where `phi^2`
/// exceeds `|t|^2 r^2`, else equal: a nested bracket's outer slope, the
/// inner point `y` known to within `r` of the nearest's (any two fields).
fn beyond_sign(x: &QV, y: &QV, t: &QV, r: &R) -> Ordering {
    let lim = qqdot(t, t).scale(&(r * r));
    if one_field(&[x, y]) {
        let phi = qqdot(&qsub(x, y), t);
        if phi.mul(&phi).sub(&lim).sign() != Ordering::Greater {
            return Ordering::Equal;
        }
        return phi.sign();
    }
    // Binary64 enclosures first: certain where they decide.
    let f = |p: &QV| -> Option<[Fast; 3]> {
        Some([
            p[0].enclose_fast(SIGN_STEPS)?,
            p[1].enclose_fast(SIGN_STEPS)?,
            p[2].enclose_fast(SIGN_STEPS)?,
        ])
    };
    if let (Some(fx), Some(fy), Some(ft), Some(fl)) =
        (f(x), f(y), f(t), lim.enclose_fast(SIGN_STEPS))
    {
        let phi = (0..3).fold(Fast::exact_f64(0.0), |acc, k| {
            acc.add(&fx[k].sub(&fy[k]).mul(&ft[k]))
        });
        match phi.square().sub(&fl).sign() {
            Some(Ordering::Less) => return Ordering::Equal,
            Some(Ordering::Greater) => {
                if let Some(o @ (Ordering::Less | Ordering::Greater)) = phi.sign() {
                    return o;
                }
            }
            _ => {}
        }
    }
    let phi = mv_dot_q(&mv_sub(&mixv(x, Side::A), &mixv(y, Side::B)), t, Side::A);
    if phi.mul(&phi).sub(&Mix::of(&lim, Side::A)).sign() != Ordering::Greater {
        return Ordering::Equal;
    }
    phi.sign()
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
                let angle = |p: &[Qd; 2]| p[1].to_f64().atan2_rn(p[0].to_f64());
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

    /// The run's 33 binary64 samples and their greatest spacing, kept for
    /// one `edges_and_vertices` (`COARSE`, by the curve's address and run:
    /// the models live through it).
    fn coarse(&self) -> std::rc::Rc<(Vec<[f64; 3]>, f64)> {
        let key = (
            self.crv as *const Crv as usize,
            self.t0.to_bits(),
            self.span.to_bits(),
        );
        if let Some(hit) = COARSE.with(|c| c.borrow().get(&key).cloned()) {
            return hit;
        }
        let pts = self.points(0.0, 1.0, 32);
        let spacing = pts
            .windows(2)
            .map(|w| dot64(&sub64(&w[0], &w[1]), &sub64(&w[0], &w[1])).sqrt())
            .fold(0.0f64, f64::max);
        let v = std::rc::Rc::new((pts, spacing));
        COARSE.with(|c| c.borrow_mut().insert(key, v.clone()));
        v
    }

    /// Whether a parameter lies within the run (its ends included).
    fn holds(&self, t: f64) -> bool {
        let (a, b) = (self.t0, self.t0 + self.span);
        a.min(b) <= t && t <= a.max(b)
    }

    /// A bracket about the fraction `s` between two exact points where
    /// `slope` (at a point with its tangent along the run) takes opposite
    /// exact signs, its half width from `2^-20` of the run widened up to
    /// `2^-4`, narrowed by bisection to the parameters' rounding: the exact
    /// point at its middle and its parameters' ends (both the middle's where
    /// the slope vanishes there).
    fn bracket(&self, s: f64, slope: &dyn Fn(&QV, &QV) -> Ordering) -> Option<(QV, [f64; 2])> {
        let sign_at = |t: f64| {
            let x = self.exact_at(t)?;
            let tangent = super::meet::tangent(self.crv, &place(self.crv, &x), &x);
            Some(slope(&x, &tangent))
        };
        let strict = |o: Ordering| o != Ordering::Equal;
        let tm = self.t0 + self.span * s;
        let mut h = self.span.abs() * 2f64.powi(-20);
        let (mut a, mut b, sa) = loop {
            let (a, b) = (tm - h, tm + h);
            if self.holds(a) && self.holds(b) {
                let (sa, sb) = (sign_at(a)?, sign_at(b)?);
                if strict(sa) && strict(sb) && sa != sb {
                    break (a, b, sa);
                }
            }
            h *= 16.0;
            if h > self.span.abs() / 16.0 {
                return None;
            }
        };
        for _ in 0..80 {
            let m = 0.5 * a + 0.5 * b;
            if m == a || m == b {
                break;
            }
            match sign_at(m)? {
                Ordering::Equal => return Some((self.exact_at(m)?, [m, m])),
                x if x == sa => a = m,
                _ => b = m,
            }
        }
        let m = 0.5 * a + 0.5 * b;
        Some((self.exact_at(m)?, [a, b]))
    }

    /// A certified bound of the curve's movement over a bracket's
    /// parameters `[ta, tb]` (`bracket`) from its point `x` there: the
    /// largest distance from `x`'s enclosure to a binary64 interval
    /// enclosure of the curve's points over the exact parameters the bracket
    /// took (an angle's rational half-angle tangents, `unit_param`; a height
    /// or a run parameter), from the curve's own formula (`MeetCrv`,
    /// `ConeSec`, `TorusSec`, `RiseCrv::enclose`; a spline cap edge's arcs by
    /// Horner's scheme). `None` for a torus's meeting with a quadric, a
    /// turned frame's rise, or a bracket whose enclosure is not certain.
    fn movement(&self, ta: f64, tb: f64, x: &QV) -> Option<R> {
        let cells = self.enclose(ta, tb)?;
        let mut sum = zero();
        for (k, cell) in cells.iter().enumerate() {
            let (bl, bh) = cell.bounds_f64();
            let (xl, xh) = x[k].enclose_fast(64)?.bounds_f64();
            let (Some(bl), Some(bh), Some(xl), Some(xh)) = (
                R::from_float(bl),
                R::from_float(bh),
                R::from_float(xl),
                R::from_float(xh),
            ) else {
                return None;
            };
            let dev = (&bh - &xl).max(&xh - &bl);
            sum += &dev * &dev;
        }
        Some(sqrt_up(&sum))
    }

    /// A binary64 interval enclosure of the curve's exact points at the
    /// parameters `ta`, `tb`, their middle, and every one between (the
    /// exact parameters `exact_at` takes there), from the curve's own
    /// formula; `ta == tb` one point's.
    fn enclose(&self, ta: f64, tb: f64) -> Option<[Fast; 3]> {
        let tm = 0.5 * ta + 0.5 * tb;
        let hull = |v: [R; 3]| {
            let lo = v.iter().min().expect("three").clone();
            let hi = v.iter().max().expect("three").clone();
            Fast::from_r(&lo).union(&Fast::from_r(&hi))
        };
        let cells = match self.crv {
            Crv::Meet(_) | Crv::Cone(_) | Crv::Torus(_) => {
                // The unit circle's arc between the exact `(cos, sin)` the
                // bracket took (`unit_at`; the short arc, its parameters'
                // rounding apart) lies in their box widened by its sagitta,
                // at most a quarter of the squared chord (bounded by the
                // box's squared diagonal).
                let [pa, pb, pm] = [ta, tb, tm].map(unit_at);
                let lo = [0, 1].map(|k| pa[k].clone().min(pb[k].clone()).min(pm[k].clone()));
                let hi = [0, 1].map(|k| pa[k].clone().max(pb[k].clone()).max(pm[k].clone()));
                let w = (0..2).fold(zero(), |acc, k| {
                    let d = &hi[k] - &lo[k];
                    acc + &d * &d
                });
                let cs = [0, 1]
                    .map(|k| Fast::from_r(&(&lo[k] - &w)).union(&Fast::from_r(&(&hi[k] + &w))));
                match self.crv {
                    Crv::Meet(c) => c.enclose(&cs)?,
                    Crv::Cone(c) => c.enclose(&cs)?,
                    Crv::Torus(c) => c.enclose(&cs)?,
                    _ => return None,
                }
            }
            Crv::Rise(c) => c.enclose(&hull([q(ta), q(tb), q(tm)]))?,
            Crv::Spline(c) => wall_enclose(c, &hull([q(ta), q(tb), q(tm)]))?,
            _ => return None,
        };
        Some(cells)
    }

    /// The fraction of the run's point nearest `p` in binary64 about the
    /// fraction `s` (a golden-section search within `1/16` of the run).
    fn nearest64(&self, p: &[f64; 3], s: f64) -> f64 {
        let f = |s: f64| match self.points(s, s, 1).first() {
            Some(x) => (0..3).map(|k| (x[k] - p[k]).powi(2)).sum::<f64>(),
            None => f64::INFINITY,
        };
        let (mut lo, mut hi) = ((s - 1.0 / 16.0).max(0.0), (s + 1.0 / 16.0).min(1.0));
        let g = (5f64.sqrt() - 1.0) / 2.0;
        for _ in 0..80 {
            let (x1, x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
            if f(x1) < f(x2) {
                hi = x2;
            } else {
                lo = x1;
            }
        }
        0.5 * lo + 0.5 * hi
    }
}

/// A binary64 interval enclosure of a spline cap edge's points over an
/// enclosure of its run parameter: each Bézier arc it meets over its part
/// of it by Horner's scheme, `w` from its plane's heights.
fn wall_enclose(c: &super::spline_walls::WallCrv, tau: &Fast) -> Option<[Fast; 3]> {
    let (lo, hi) = tau.bounds_f64();
    let r = |x: &R| Fast::from_r(x);
    let horner = |p: &[R], s: &Fast| {
        p.iter()
            .rev()
            .fold(Fast::exact_f64(0.0), |acc, k| acc.mul(s).add(&r(k)))
    };
    let mut out: Option<[Fast; 3]> = None;
    for arc in &c.seg.arcs {
        let (d0, d1) = (rational_f64(&arc.d[0]), rational_f64(&arc.d[1]));
        if hi < d0.min(d1) || lo > d0.max(d1) {
            continue;
        }
        let s = tau
            .sub(&r(&arc.d[0]))
            .div(&r(&arc.d[1]).sub(&r(&arc.d[0])))?;
        let (u, v) = (horner(&arc.x, &s), horner(&arc.y, &s));
        let w = r(&c.h[0]).add(&r(&c.h[1]).mul(&u)).add(&r(&c.h[2]).mul(&v));
        let f = &c.f;
        let p = [0, 1, 2].map(|j| {
            r(&f.o[j])
                .add(&r(&f.x[j]).mul(&u))
                .add(&r(&f.y[j]).mul(&v))
                .add(&r(&f.n[j]).mul(&w))
        });
        out = Some(match out {
            None => p,
            Some(o) => [0, 1, 2].map(|j| o[j].union(&p[j])),
        });
    }
    out
}

/// The leasts of an edge's other curve's (`run`) distance to another
/// input's point `other(x)` (a vertex, or a line's foot of the run's point
/// `x`), `far` its binary64 square: each least among the run's samples
/// within twice the resolution in binary64, settled on the exact sign of
/// `(x - y) . t` (`Run::bracket`: the run's exact point at the bracket's
/// middle, or at the binary64 place where no bracket is found) and the
/// other's point; none where their gap is of rounding's size (an
/// incidence: no contact); within the resolution those points; beyond it
/// with the curve's movement over the bracket (`Run::movement`), or
/// unbracketed. `side`: the run's side of the pair.
fn run_nearest(
    run: &Run,
    far: &dyn Fn([f64; 3]) -> f64,
    far_fast: &dyn Fn(&[Fast; 3]) -> Option<Fast>,
    other: &dyn Fn(&QV) -> Option<Pt>,
    side: Side,
    res: f64,
) -> Vec<(QV, Pt, Bracket)> {
    // Passed over whole where every sample lies beyond twice the
    // resolution and twice the samples' spacing (the curve's movement
    // between them): the cost of most pairs whose boxes meet.
    let coarse = run.coarse();
    let (pts, spacing) = (&coarse.0, coarse.1);
    let reach = 2.0 * res + 2.0 * spacing;
    if pts.iter().all(|x| far(*x) > reach * reach) {
        return Vec::new();
    }
    let tol = q(res);
    let slope = |x: &QV, t: &QV| match other(x) {
        Some(y) => pt_sign(x, &y, t, side),
        None => Ordering::Equal,
    };
    let mut out = Vec::new();
    for (s, _) in run.extrema(far, false) {
        let Some(x64) = run.points(s, s, 1).first().copied() else {
            continue;
        };
        let gap = far(x64).max(0.0).sqrt();
        if gap.is_nan() || gap > 2.0 * res {
            continue;
        }
        // Of rounding's size beyond doubt in binary64 (an incidence: the
        // run's exact point at the binary64 place enclosed by its formula,
        // its squared distance to the other's enclosed below `rho^2`, `rho`
        // of a size no larger than the points'): no contact, no exact point.
        if let Some(cells) = run.enclose(run.t0 + run.span * s, run.t0 + run.span * s) {
            if let Some(d2) = far_fast(&cells) {
                let size = cells
                    .iter()
                    .map(|c| {
                        let (lo, hi) = c.bounds_f64();
                        if lo > 0.0 {
                            lo
                        } else if hi < 0.0 {
                            -hi
                        } else {
                            0.0
                        }
                    })
                    .fold(1.0f64, f64::max);
                let rho = 2f64.powi(size.log2_rn().ceil() as i32 - 48);
                if d2.bounds_f64().1 <= rho * rho {
                    continue;
                }
            }
        }
        // Settled: the bracket's middle (an incidence then of rounding's
        // size), or the binary64 place's exact point where none is found.
        let (x, span) = match run.bracket(s, &slope) {
            Some((x, span)) => (x, Some(span)),
            None => match run.exact(s) {
                Some(x) => (x, None),
                None => continue,
            },
        };
        let Some(y) = other(&x) else {
            continue;
        };
        match pt_cmp(&x, &y, side, None) {
            Some(Ordering::Greater) => {}
            _ => continue,
        }
        if pt_cmp(&x, &y, side, Some(&(&tol * &tol))) != Some(Ordering::Greater) {
            out.push((x, y, Bracket::Exact));
            continue;
        }
        let bound = match span.and_then(|[ta, tb]| run.movement(ta, tb, &x)) {
            Some(d) => Bracket::Within(d),
            None => Bracket::Unbracketed,
        };
        out.push((x, y, bound));
    }
    out
}

/// The sign of `(x - y) . t` for a point `x` (on `side`, `t` of its field)
/// and a pair's point `y` of any field.
fn pt_sign(x: &QV, y: &Pt, t: &QV, side: Side) -> Ordering {
    match &y.exact {
        None => rel_sign(x, &y.p, t),
        Some(e) => mv_dot_q(&mv_sub(&mixv(x, side), e), t, side).sign(),
    }
}

/// The sign of `|y - x|^2 - r2` for a point `x` (on `side`) and a pair's
/// point `y`, `r2` the gap's own rounding where `None` (`rho^2`).
fn pt_cmp(x: &QV, y: &Pt, side: Side, r2: Option<&R>) -> Option<Ordering> {
    match &y.exact {
        None => Some(sq_cmp(
            x,
            &y.p,
            &r2.cloned().unwrap_or_else(|| rounding_of(x, &y.p)),
        )),
        Some(e) => {
            let xm = mixv(x, side);
            let n = mv_sub(e, &xm);
            let r2 = match r2 {
                Some(r) => r.clone(),
                None => mixed_rounding(&xm, e),
            };
            Some(mv_dot(&n, &n).sub(&Mix::r(&r2)).sign())
        }
    }
}

/// A curve walked by the nested brackets of two curves: a conic or a
/// circle by its angle (`curve_at`, `bisect` with its speed bound), or
/// another curve by its run's fraction (`Run::bracket`, `Run::movement`).
enum Walk<'a> {
    Arc(&'a Crv, Conic64, R),
    Run(Run<'a>),
}

impl Walk<'_> {
    /// The exact point near a binary64 place (an angle, or a fraction).
    fn exact(&self, t: f64) -> Option<QV> {
        match self {
            Walk::Arc(crv, f, _) => {
                let (u, flip) = f.param(t);
                Some(curve_at(crv, &q(u), flip)?.0)
            }
            Walk::Run(run) => run.exact(t),
        }
    }

    /// The binary64 place of the point nearest `p` about the place `t`.
    fn nearest64(&self, p: &[f64; 3], t: f64) -> f64 {
        match self {
            Walk::Arc(_, f, _) => f.nearest(p, t),
            Walk::Run(run) => run.nearest64(p, t),
        }
    }

    /// The bracket about the place `t` on `slope`'s exact signs: the
    /// middle's point and the curve's movement over the bracket.
    fn bracket(&self, t: f64, slope: &dyn Fn(&QV, &QV) -> Ordering) -> Option<(QV, R)> {
        match self {
            Walk::Arc(crv, f, speed) => {
                let (u, flip) = f.param(t);
                match bisect(crv, u, flip, slope, speed)? {
                    (x, Bracket::Within(r)) => Some((x, r)),
                    _ => None,
                }
            }
            Walk::Run(run) => {
                let (x, [ta, tb]) = run.bracket(t, slope)?;
                let d = run.movement(ta, tb, &x)?;
                Some((x, d))
            }
        }
    }

    /// Binary64 samples over the whole curve and a bound of its movement
    /// between neighbours (a conic's `reach` times the step, a run's
    /// greatest gap between samples, doubled).
    fn samples(&self) -> (Vec<[f64; 3]>, f64) {
        match self {
            Walk::Arc(_, f, _) => {
                const N: usize = 64;
                let step = TAU / N as f64;
                (
                    (0..N).map(|i| f.jet(i as f64 * step)[0]).collect(),
                    f.reach() * step,
                )
            }
            Walk::Run(run) => {
                let coarse = run.coarse();
                (coarse.0.clone(), 2.0 * coarse.1)
            }
        }
    }
}

/// Edge `ei`'s walk (`crv` its curve): a conic's or circle's (its binary64
/// form and its speed bound), or another curve's run.
fn walk_of<'a>(m: &'a Prism, ei: usize, crv: &'a Crv) -> Option<Walk<'a>> {
    match crv {
        Crv::Conic { .. } | Crv::Circle(_) => Some(Walk::Arc(crv, Conic64::of(crv)?, speed(crv)?)),
        _ => Some(Walk::Run(Run::new(m, ei)?)),
    }
}

/// The closest points of two edges of which one at least is another
/// curve's (`Run`): a line and another curve, the run's leasts of its
/// distance to the line (`run_nearest`, the foot exact); a conic or a
/// circle or another curve against another curve, the leasts of their
/// distance found in binary64 (the other curve's samples, each one's
/// nearest point on the first), decided at the exact points there within
/// the resolution, beyond it bracketed nested (`nested`) where the pair is
/// a contact at those points inside both edges. `ma`'s point first.
fn run_pairs(ma: &Prism, ea: usize, mb: &Prism, eb: usize, res: f64) -> Vec<(Pt, Pt, Bracket)> {
    let (ca, cb) = (&ma.edges[ea].curve, &mb.edges[eb].curve);
    // A line against another curve.
    for (o, (ml, el), (mr, er)) in [(0, (ma, ea), (mb, eb)), (1, (mb, eb), (ma, ea))] {
        let Crv::Line { p, d } = &ml.edges[el].curve else {
            continue;
        };
        let Some(run) = Run::new(mr, er) else {
            return Vec::new();
        };
        // The run's side of the pair.
        let side = if o == 0 { Side::B } else { Side::A };
        let (fp, fd) = (qv_f64(p), d.clone().map(|x| rational_f64(&x)));
        let fdd = dot64(&fd, &fd);
        let far = |x: [f64; 3]| {
            let w = sub64(&x, &fp);
            let c = [
                w[1] * fd[2] - w[2] * fd[1],
                w[2] * fd[0] - w[0] * fd[2],
                w[0] * fd[1] - w[1] * fd[0],
            ];
            dot64(&c, &c) / fdd
        };
        let foot = |x: &QV| -> Option<Pt> {
            if one_field(&[x, p]) {
                let t = qdot(&qsub(x, p), d).scale(&(int(1) / dot(d, d)));
                return Some(Pt::of(qadd(p, &qscale(d, &t))));
            }
            Some(line_foot(ml, el, x, side)?.1)
        };
        let (of, df) = (fast_point(p), d.clone().map(|x| Fast::from_r(&x)));
        let dd = Fast::from_r(&dot(d, d));
        let far_fast = |c: &[Fast; 3]| {
            let of = of.as_ref()?;
            let w = [0, 1, 2].map(|k| c[k].sub(&of[k]));
            let x = [
                w[1].mul(&df[2]).sub(&w[2].mul(&df[1])),
                w[2].mul(&df[0]).sub(&w[0].mul(&df[2])),
                w[0].mul(&df[1]).sub(&w[1].mul(&df[0])),
            ];
            x[0].square()
                .add(&x[1].square())
                .add(&x[2].square())
                .div(&dd)
        };
        return run_nearest(&run, &far, &far_fast, &foot, side, res)
            .into_iter()
            .map(|(x, y, k)| {
                if o == 0 {
                    (y, Pt::of(x), k)
                } else {
                    (Pt::of(x), y, k)
                }
            })
            .collect();
    }
    let (Some(wa), Some(wb)) = (walk_of(ma, ea, ca), walk_of(mb, eb, cb)) else {
        return Vec::new();
    };
    let tol = q(res);
    let mut out = Vec::new();
    for (t, s) in walk_leasts(&wa, &wb, res) {
        let (Some(mut x), Some(mut y)) = (wa.exact(t), wb.exact(s)) else {
            continue;
        };
        let r2 = rounding_of(&x, &y);
        if sq_cmp(&x, &y, &r2) != Ordering::Greater {
            continue;
        }
        // Within 64 times the gap's rounding: settled (`nested`), an
        // incidence's binary64 places' rounding then gone.
        if sq_cmp(&x, &y, &(&r2 * int(4096))) != Ordering::Greater {
            if let Some((sx, sy, _)) = nested(&wa, t, &wb, s) {
                (x, y) = (sx, sy);
                if sq_cmp(&x, &y, &rounding_of(&x, &y)) != Ordering::Greater {
                    continue;
                }
            }
        }
        if sq_cmp(&x, &y, &(&tol * &tol)) != Ordering::Greater {
            out.push((Pt::of(x), Pt::of(y), Bracket::Exact));
            continue;
        }
        let (sa, sb) = (Spot::Edge(ea), Spot::Edge(eb));
        let (a, b) = (End::at(ma, &sa, x.clone()), End::at(mb, &sb, y.clone()));
        if contact_any(&a, &b).is_none() || !inside_edge(ma, ea, &x) || !inside_edge(mb, eb, &y) {
            continue;
        }
        out.push(match nested(&wa, t, &wb, s) {
            Some((x, y, d)) => (Pt::of(x), Pt::of(y), Bracket::Within(d)),
            None => (Pt::of(x), Pt::of(y), Bracket::Unbracketed),
        });
    }
    out
}

/// The places of the leasts of two walks' distance within twice the
/// resolution in binary64 (as `conic_leasts`): over the grid of both
/// curves' samples, each node at most its neighbours' (an arc's angles
/// wrapping) whose distance less both curves' movements between samples is
/// within twice the resolution, refined by alternating golden-section
/// searches of each curve's place within a step of its node, kept where
/// its distance is within twice the resolution, once.
fn walk_leasts(wa: &Walk, wb: &Walk, res: f64) -> Vec<(f64, f64)> {
    let ((pa, ma), (pb, mb)) = (wa.samples(), wb.samples());
    let (na, nb) = (pa.len(), pb.len());
    let reach = 2.0 * res + ma + mb;
    let mut d2 = vec![0.0; na * nb];
    let mut any = false;
    for (i, x) in pa.iter().enumerate() {
        for (j, y) in pb.iter().enumerate() {
            let w = sub64(x, y);
            let v = dot64(&w, &w);
            if !v.is_finite() {
                return Vec::new();
            }
            any |= v <= reach * reach;
            d2[i * nb + j] = v;
        }
    }
    if !any {
        return Vec::new();
    }
    // A walk's place of sample `i`, its step, and its neighbours (an arc's
    // wrapping).
    let geometry = |w: &Walk, n: usize| match w {
        Walk::Arc(..) => (TAU / n as f64, true),
        Walk::Run(_) => (1.0 / (n - 1) as f64, false),
    };
    let ((sa, wrap_a), (sb, wrap_b)) = (geometry(wa, na), geometry(wb, nb));
    let near = |i: usize, n: usize, wrap: bool| -> Vec<usize> {
        let mut v = Vec::new();
        for d in [-1i64, 0, 1] {
            let k = i as i64 + d;
            if wrap {
                v.push(k.rem_euclid(n as i64) as usize);
            } else if (0..n as i64).contains(&k) {
                v.push(k as usize);
            }
        }
        v
    };
    let dist = |t: f64, s: f64| match (wa.exact64(t), wb.exact64(s)) {
        (Some(x), Some(y)) => {
            let w = sub64(&x, &y);
            dot64(&w, &w)
        }
        _ => f64::INFINITY,
    };
    // The place on one curve, the other's held, nearest within `h`, by
    // `k` golden-section steps.
    let golden = |f: &dyn Fn(f64) -> f64, c: f64, h: f64, lo: f64, hi: f64, k: usize| {
        let (mut a, mut b) = ((c - h).max(lo), (c + h).min(hi));
        let g = (5f64.sqrt() - 1.0) / 2.0;
        for _ in 0..k {
            let (x1, x2) = (b - g * (b - a), a + g * (b - a));
            if f(x1) < f(x2) {
                b = x2;
            } else {
                a = x1;
            }
        }
        0.5 * a + 0.5 * b
    };
    let range = |wrap: bool| {
        if wrap {
            (f64::NEG_INFINITY, f64::INFINITY)
        } else {
            (0.0, 1.0)
        }
    };
    let ((la, ha), (lb, hb)) = (range(wrap_a), range(wrap_b));
    let mut out: Vec<(f64, f64)> = Vec::new();
    let close = |p: f64, q: f64| (p - q).abs() <= 1e-9;
    for i in 0..na {
        for j in 0..nb {
            let v = d2[i * nb + j];
            if v > reach * reach {
                continue;
            }
            let least = near(i, na, wrap_a)
                .iter()
                .all(|&k| near(j, nb, wrap_b).iter().all(|&l| v <= d2[k * nb + l]));
            if !least {
                continue;
            }
            let (mut t, mut s) = (i as f64 * sa, j as f64 * sb);
            // Coarsely first (each place then within `0.618^16` of a step
            // of the least's, the distance off it by about the curves'
            // movements over that): most nodes are passed over.
            for _ in 0..2 {
                t = golden(&|t| dist(t, s), t, sa, la, ha, 16);
                s = golden(&|s| dist(t, s), s, sb, lb, hb, 16);
            }
            let slack = 0.618f64.powi(16) * (ma + mb);
            if dist(t, s).sqrt() > 2.0 * res + 4.0 * slack {
                continue;
            }
            // Then Gauss-Newton steps on `A(t) - B(s)` (central differences
            // of the curves' binary64 points for their tangents), within a
            // step of the node: quadratic at a crossing or an incidence,
            // where alternating searches crawl.
            for _ in 0..40 {
                let (Some(x), Some(y)) = (wa.exact64(t), wb.exact64(s)) else {
                    break;
                };
                let tangent = |w: &Walk, u: f64, h: f64| -> Option<[f64; 3]> {
                    let (p, q) = (w.exact64(u + h)?, w.exact64(u - h)?);
                    Some([0, 1, 2].map(|k| (p[k] - q[k]) / (2.0 * h)))
                };
                let (ha_, hb_) = (sa * 1e-6, sb * 1e-6);
                let (Some(ta), Some(tb)) = (tangent(wa, t, ha_), tangent(wb, s, hb_)) else {
                    break;
                };
                let r = sub64(&x, &y);
                let (aa, ab, bb) = (dot64(&ta, &ta), dot64(&ta, &tb), dot64(&tb, &tb));
                let (ga, gb) = (dot64(&ta, &r), -dot64(&tb, &r));
                // [aa, -ab; -ab, bb] (dt, ds) = -(ga, gb).
                let det = aa * bb - ab * ab;
                if !(det > 1e-24 * aa * bb && det.is_finite()) {
                    break;
                }
                let dt = (-ga * bb - ab * gb) / det;
                let ds = (-gb * aa - ab * ga) / det;
                if !(dt.is_finite() && ds.is_finite()) {
                    break;
                }
                let (nt, ns) = (
                    (t + dt.clamp(-sa, sa)).clamp(la, ha),
                    (s + ds.clamp(-sb, sb)).clamp(lb, hb),
                );
                if dist(nt, ns) > dist(t, s) {
                    break;
                }
                let done = nt == t && ns == s;
                (t, s) = (nt, ns);
                if done {
                    break;
                }
            }
            let gap = dist(t, s).sqrt();
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

impl Walk<'_> {
    /// The binary64 point at a place.
    fn exact64(&self, t: f64) -> Option<[f64; 3]> {
        match self {
            Walk::Arc(_, f, _) => Some(f.jet(t)[0]),
            Walk::Run(run) => run.points(t, t, 1).first().copied(),
        }
    }
}

/// Two curves' closest points about the places found in binary64,
/// bracketed nested (as `conics_bracketed`): on the outer curve the exact
/// distance to the inner one, its derivative's sign that of `(x - y) .
/// t_A` at the inner curve's point `y` nearest `x` (bracketed on it with its
/// bound `r_B`), taken as zero within `|t_A| r_B`; the points (the first's
/// first) and the bound `D_A + r_B`.
fn nested(outer: &Walk, to: f64, inner: &Walk, ti: f64) -> Option<(QV, QV, R)> {
    let near = |x: &QV| -> Option<(QV, R)> {
        let t = inner.nearest64(&qv_f64(x), ti);
        inner.bracket(t, &|y: &QV, tb: &QV| rel_sign(y, x, tb))
    };
    let slope = |x: &QV, ta: &QV| {
        let Some((y, r)) = near(x) else {
            return Ordering::Equal;
        };
        beyond_sign(x, &y, ta, &r)
    };
    let (x, da) = outer.bracket(to, &slope)?;
    let (y, r) = near(&x)?;
    Some((x, y, da + r))
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

// ------------------------------------------------------------ faces

/// Two cylinders' or cones' walls, or a plane and one, within the
/// resolution of tangency (REVIEW_NOTES.md, "Faces with no near-miss rule of
/// their own").
pub(super) const NEAR_WALLS: &str =
    "a cylinder's or cone's wall within the resolution of tangency to the other input's (S9c)";
/// A torus within the resolution of tangency to the other input's face.
pub(super) const NEAR_TORUS: &str =
    "a torus within the resolution of tangency to the other input's surface (S9d.4)";
/// A spline wall within the resolution of tangency to the other input's
/// face.
pub(super) const NEAR_SPLINE: &str =
    "a spline wall within the resolution of tangency to the other input's surface (S9f)";

/// The arithmetic of the face pairs' searches: binary64 (`f64`), forward
/// jets of first and second derivatives in two parameters (`J`), and
/// binary64 intervals (`Iv`) for the certified bounds beyond the resolution.
trait Num: Clone {
    /// A binary64 constant.
    fn k(x: f64) -> Self;
    /// A rational datum (enclosed in the interval arithmetic).
    fn of_r(x: &R) -> Self;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    fn div(&self, o: &Self) -> Self;
    fn sqrt(&self) -> Self;
    fn cos_sin(&self) -> (Self, Self);
    /// A binary64 view of the value (an interval's midpoint).
    fn mid(&self) -> f64;
    /// The square (an interval's never below zero).
    fn sq(&self) -> Self {
        self.mul(self)
    }
}

impl Num for f64 {
    fn k(x: f64) -> Self {
        x
    }
    fn of_r(x: &R) -> Self {
        rational_f64(x)
    }
    fn add(&self, o: &Self) -> Self {
        self + o
    }
    fn sub(&self, o: &Self) -> Self {
        self - o
    }
    fn mul(&self, o: &Self) -> Self {
        self * o
    }
    fn div(&self, o: &Self) -> Self {
        self / o
    }
    fn sqrt(&self) -> Self {
        f64::sqrt(*self)
    }
    fn cos_sin(&self) -> (Self, Self) {
        (self.cos_rn(), self.sin_rn())
    }
    fn mid(&self) -> f64 {
        *self
    }
}

/// A binary64 interval (`Fast`, outward rounded); a division by an interval
/// holding zero gives NaN bounds, never certain.
#[derive(Clone, Copy, Debug)]
struct Iv(Fast);

impl Iv {
    fn span(lo: f64, hi: f64) -> Self {
        Iv(Fast::exact_f64(lo).union(&Fast::exact_f64(hi)))
    }

    fn bounds(&self) -> (f64, f64) {
        self.0.bounds_f64()
    }

    /// Whether it lies strictly within `(lo, hi)`.
    fn inside(&self, lo: f64, hi: f64) -> bool {
        let (a, b) = self.bounds();
        lo < a && b < hi
    }
}

impl Num for Iv {
    fn k(x: f64) -> Self {
        Iv(Fast::exact_f64(x))
    }
    fn of_r(x: &R) -> Self {
        Iv(Fast::from_r(x))
    }
    fn add(&self, o: &Self) -> Self {
        Iv(self.0.add(&o.0))
    }
    fn sub(&self, o: &Self) -> Self {
        Iv(self.0.sub(&o.0))
    }
    fn mul(&self, o: &Self) -> Self {
        Iv(self.0.mul(&o.0))
    }
    fn div(&self, o: &Self) -> Self {
        match self.0.div(&o.0) {
            Some(x) => Iv(x),
            None => Iv(Fast::exact_f64(f64::NAN)),
        }
    }
    fn sqrt(&self) -> Self {
        let (lo, _) = self.bounds();
        if lo.is_nan() || lo < 0.0 {
            return Iv(Fast::exact_f64(f64::NAN));
        }
        Iv(self.0.sqrt())
    }
    fn cos_sin(&self) -> (Self, Self) {
        let (c, s) = Fast::cos_sin(&self.0);
        (Iv(c), Iv(s))
    }
    fn mid(&self) -> f64 {
        let (lo, hi) = self.bounds();
        0.5 * lo + 0.5 * hi
    }
    fn sq(&self) -> Self {
        Iv(self.0.square())
    }
}

/// A forward jet in two parameters: a value, its gradient and its Hessian
/// (`h = [d00, d01, d11]`).
#[derive(Clone, Debug)]
struct J<B> {
    v: B,
    d: [B; 2],
    h: [B; 3],
}

impl<B: Num> J<B> {
    fn con(v: B) -> Self {
        let z = B::k(0.0);
        J {
            v,
            d: [z.clone(), z.clone()],
            h: [z.clone(), z.clone(), z],
        }
    }

    /// The parameter `i` at the value `v`.
    fn var(v: B, i: usize) -> Self {
        let mut j = Self::con(v);
        j.d[i] = B::k(1.0);
        j
    }

    /// `f` of this jet from `f`, `f'` and `f''` at its value.
    fn chain(&self, f: B, f1: B, f2: B) -> Self {
        let d = [self.d[0].mul(&f1), self.d[1].mul(&f1)];
        let h = [
            f2.mul(&self.d[0].mul(&self.d[0])).add(&f1.mul(&self.h[0])),
            f2.mul(&self.d[0].mul(&self.d[1])).add(&f1.mul(&self.h[1])),
            f2.mul(&self.d[1].mul(&self.d[1])).add(&f1.mul(&self.h[2])),
        ];
        J { v: f, d, h }
    }
}

impl<B: Num> Num for J<B> {
    fn k(x: f64) -> Self {
        Self::con(B::k(x))
    }
    fn of_r(x: &R) -> Self {
        Self::con(B::of_r(x))
    }
    fn add(&self, o: &Self) -> Self {
        J {
            v: self.v.add(&o.v),
            d: [self.d[0].add(&o.d[0]), self.d[1].add(&o.d[1])],
            h: [
                self.h[0].add(&o.h[0]),
                self.h[1].add(&o.h[1]),
                self.h[2].add(&o.h[2]),
            ],
        }
    }
    fn sub(&self, o: &Self) -> Self {
        J {
            v: self.v.sub(&o.v),
            d: [self.d[0].sub(&o.d[0]), self.d[1].sub(&o.d[1])],
            h: [
                self.h[0].sub(&o.h[0]),
                self.h[1].sub(&o.h[1]),
                self.h[2].sub(&o.h[2]),
            ],
        }
    }
    fn mul(&self, o: &Self) -> Self {
        let (a, b) = (self, o);
        let cross = |i: usize, j: usize| a.d[i].mul(&b.d[j]).add(&a.d[j].mul(&b.d[i]));
        J {
            v: a.v.mul(&b.v),
            d: [
                a.d[0].mul(&b.v).add(&a.v.mul(&b.d[0])),
                a.d[1].mul(&b.v).add(&a.v.mul(&b.d[1])),
            ],
            h: [
                a.h[0].mul(&b.v).add(&cross(0, 0)).add(&a.v.mul(&b.h[0])),
                a.h[1].mul(&b.v).add(&cross(0, 1)).add(&a.v.mul(&b.h[1])),
                a.h[2].mul(&b.v).add(&cross(1, 1)).add(&a.v.mul(&b.h[2])),
            ],
        }
    }
    fn div(&self, o: &Self) -> Self {
        let one = B::k(1.0);
        let r = one.div(&o.v);
        let r2 = r.mul(&r);
        let recip = o.chain(r.clone(), B::k(0.0).sub(&r2), B::k(2.0).mul(&r2.mul(&r)));
        self.mul(&recip)
    }
    fn sqrt(&self) -> Self {
        let s = self.v.sqrt();
        let f1 = B::k(0.5).div(&s);
        let f2 = B::k(-0.25).div(&s.mul(&self.v));
        self.chain(s, f1, f2)
    }
    fn cos_sin(&self) -> (Self, Self) {
        let (c, s) = self.v.cos_sin();
        let neg = |x: &B| B::k(0.0).sub(x);
        (
            self.chain(c.clone(), neg(&s), neg(&c)),
            self.chain(s.clone(), c, neg(&s)),
        )
    }
    fn mid(&self) -> f64 {
        self.v.mid()
    }
    fn sq(&self) -> Self {
        let two = B::k(2.0);
        self.chain(self.v.sq(), two.mul(&self.v), two)
    }
}

/// Whether `a > b` (neither NaN).
fn above(a: f64, b: f64) -> bool {
    a.partial_cmp(&b) == Some(Ordering::Greater)
}

/// Whether `a <= b` (neither NaN).
fn at_most(a: f64, b: f64) -> bool {
    matches!(a.partial_cmp(&b), Some(Ordering::Less | Ordering::Equal))
}

fn dot3<T: Num>(a: &[T; 3], b: &[T; 3]) -> T {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

fn horner<T: Num>(c: &[T], s: &T) -> T {
    c.iter().rev().fold(T::k(0.0), |acc, ci| acc.mul(s).add(ci))
}

/// A spline wall's arc in any arithmetic: its run range (`d0`, `d1` in
/// binary64 to pick it, `d0t` and `inv = 1 / (d1 - d0)` in its arithmetic)
/// and its coordinates' coefficients in powers of its own parameter `s`
/// (`tau = d0 + s (d1 - d0)`), with their first and second derivatives'.
#[derive(Clone, Debug)]
struct PArc<T> {
    d0: f64,
    d1: f64,
    d0t: T,
    inv: T,
    x: [Vec<T>; 3],
    y: [Vec<T>; 3],
}

fn parcs_of<T: Num>(seg: &super::spline_walls::SplineSeg) -> Vec<PArc<T>> {
    let deriv = |c: &[R]| -> Vec<R> {
        c.iter()
            .enumerate()
            .skip(1)
            .map(|(i, x)| x * int(i as i64))
            .collect()
    };
    let lift = |c: &[R]| c.iter().map(T::of_r).collect::<Vec<T>>();
    let three = |c: &[R]| {
        let d1 = deriv(c);
        let d2 = deriv(&d1);
        [lift(c), lift(&d1), lift(&d2)]
    };
    seg.arcs
        .iter()
        .map(|a| PArc {
            d0: rational_f64(&a.d[0]),
            d1: rational_f64(&a.d[1]),
            d0t: T::of_r(&a.d[0]),
            inv: T::of_r(&(int(1) / (&a.d[1] - &a.d[0]))),
            x: three(&a.x),
            y: three(&a.y),
        })
        .collect()
}

/// A spline wall's arcs in binary64.
fn parcs(seg: &super::spline_walls::SplineSeg) -> Vec<PArc<f64>> {
    parcs_of(seg)
}

/// The arc holding a run parameter (the first whose range reaches it).
fn arc_at<T>(arcs: &[PArc<T>], tau: f64) -> &PArc<T> {
    arcs.iter()
        .find(|a| tau <= a.d1)
        .unwrap_or(&arcs[arcs.len() - 1])
}

/// A profile's point and its first and second derivatives in the run
/// parameter.
fn profile_jet<T: Num>(arcs: &[PArc<T>], tau: &T) -> [[T; 2]; 3] {
    let a = arc_at(arcs, tau.mid());
    let s = tau.sub(&a.d0t).mul(&a.inv);
    let i2 = a.inv.mul(&a.inv);
    [
        [horner(&a.x[0], &s), horner(&a.y[0], &s)],
        [
            horner(&a.x[1], &s).mul(&a.inv),
            horner(&a.y[1], &s).mul(&a.inv),
        ],
        [horner(&a.x[2], &s).mul(&i2), horner(&a.y[2], &s).mul(&i2)],
    ]
}

/// The run parameter of a profile's point nearest a point (binary64: 32
/// samples an arc, then Newton's method; a search only).
fn near_run(arcs: &[PArc<f64>], p: [f64; 2]) -> Option<f64> {
    let mut best = (f64::INFINITY, 0.0f64);
    for a in arcs {
        for i in 0..=32 {
            let tau = a.d0 + (a.d1 - a.d0) * f64::from(i) / 32.0;
            let [c, ..] = profile_jet(arcs, &tau);
            let d = (c[0] - p[0]).powi(2) + (c[1] - p[1]).powi(2);
            if d < best.0 {
                best = (d, tau);
            }
        }
    }
    if !best.0.is_finite() {
        return None;
    }
    let (lo, hi) = (arcs[0].d0, arcs[arcs.len() - 1].d1);
    let mut tau = best.1;
    for _ in 0..32 {
        let [c, d1, d2] = profile_jet(arcs, &tau);
        let e = [c[0] - p[0], c[1] - p[1]];
        let g = e[0] * d1[0] + e[1] * d1[1];
        let h = d1[0] * d1[0] + d1[1] * d1[1] + e[0] * d2[0] + e[1] * d2[1];
        if !above(h, 0.0) {
            break;
        }
        let next = (tau - g / h).clamp(lo, hi);
        if next == tau {
            break;
        }
        tau = next;
    }
    tau.is_finite().then_some(tau)
}

/// A rational run parameter at a binary64 one, within the segment's run.
fn run_rational(seg: &super::spline_walls::SplineSeg, tau: f64) -> R {
    let t = q(tau);
    if t < seg.first {
        seg.first.clone()
    } else if t > seg.last {
        seg.last.clone()
    } else {
        t
    }
}

/// A profile point's signed distance from a spline segment (positive on its
/// outward side: the run's right on an outer boundary, its left on a
/// hole's) and the unit direction of its growth (binary64).
fn profile_distance(arcs: &[PArc<f64>], hole: bool, p: [f64; 2]) -> (f64, [f64; 2]) {
    let Some(tau) = near_run(arcs, p) else {
        return (f64::NAN, [0.0; 2]);
    };
    let [c, d1, _] = profile_jet(arcs, &tau);
    let len = d1[0].hypot_rn(d1[1]);
    let s = if hole { -1.0 } else { 1.0 };
    let m = [s * d1[1] / len, -s * d1[0] / len];
    let d = (p[0] - c[0]) * m[0] + (p[1] - c[1]) * m[1];
    (d, m)
}

/// The kind of a face sampled in its parameters (`P`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PKind {
    Cyl,
    Cone,
    Torus,
    Spline,
}

/// A face's surface as a map of two parameters: a cylinder's or a cone's
/// angle and height, a torus's angles about its axis and its tube, a spline
/// wall's run parameter and height; its frame and shape in any arithmetic.
#[derive(Clone, Debug)]
struct Patch<T> {
    kind: PKind,
    o: [T; 3],
    axes: [[T; 3]; 3],
    /// A cylinder's `(cu, cv, r)`, a cone's `(b, k, 0)`, a torus's `(R, r,
    /// 0)`.
    c: [T; 3],
    arcs: Vec<PArc<T>>,
}

impl<T: Num> Patch<T> {
    fn of(m: &Prism, fi: usize) -> Option<Self> {
        let (vm, vi) = m.view(fi);
        let f = &vm.f;
        let v3 = |v: &V| [T::of_r(&v[0]), T::of_r(&v[1]), T::of_r(&v[2])];
        let z = || T::k(0.0);
        let (kind, c, arcs) = match &vm.faces[vi].surf {
            Surf::Cyl { c, r, .. } => (
                PKind::Cyl,
                [T::of_r(&c[0]), T::of_r(&c[1]), T::of_r(r)],
                Vec::new(),
            ),
            Surf::Cone { b, k } => (PKind::Cone, [T::of_r(b), T::of_r(k), z()], Vec::new()),
            Surf::Torus => {
                let ring = vm.ring.as_ref()?;
                (
                    PKind::Torus,
                    [T::of_r(&ring.big), T::of_r(&ring.small), z()],
                    Vec::new(),
                )
            }
            Surf::Spline(seg) => (PKind::Spline, [z(), z(), z()], parcs_of(seg)),
            _ => return None,
        };
        Some(Self {
            kind,
            o: v3(&f.o),
            axes: [v3(&f.x), v3(&f.y), v3(&f.n)],
            c,
            arcs,
        })
    }

    /// The local coordinates at parameters.
    fn local(&self, p: &[T; 2]) -> [T; 3] {
        match self.kind {
            PKind::Cyl => {
                let (cs, sn) = p[0].cos_sin();
                [
                    self.c[0].add(&self.c[2].mul(&cs)),
                    self.c[1].add(&self.c[2].mul(&sn)),
                    p[1].clone(),
                ]
            }
            PKind::Cone => {
                let (cs, sn) = p[0].cos_sin();
                let rad = self.c[0].add(&self.c[1].mul(&p[1]));
                [rad.mul(&cs), rad.mul(&sn), p[1].clone()]
            }
            PKind::Torus => {
                let (cu, su) = p[0].cos_sin();
                let (cv, sv) = p[1].cos_sin();
                let rad = self.c[0].add(&self.c[1].mul(&cv));
                [rad.mul(&cu), rad.mul(&su), self.c[1].mul(&sv)]
            }
            PKind::Spline => {
                let [c, ..] = profile_jet(&self.arcs, &p[0]);
                let [u, v] = c;
                [u, v, p[1].clone()]
            }
        }
    }

    fn point(&self, p: &[T; 2]) -> [T; 3] {
        let l = self.local(p);
        [0, 1, 2].map(|k| {
            self.o[k]
                .add(&self.axes[0][k].mul(&l[0]))
                .add(&self.axes[1][k].mul(&l[1]))
                .add(&self.axes[2][k].mul(&l[2]))
        })
    }
}

/// The kind of a face held by its signed distance (`Q`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum QKind {
    Plane,
    Sphere,
    Cyl,
    Cone,
    Torus,
    Spline,
}

/// A face's surface held by its signed distance, positive where its
/// function (`Quad::value`'s, the sphere's `|x - c|^2 - r^2`, the torus's
/// quartic, a spline wall's outward offset) grows, as a function of `P`'s
/// local coordinates: its own coordinates `y = t + a l` of `P`'s local `l`
/// (composed exactly, so that a coordinate along `P`'s axis enters only
/// through its coefficient): a plane's offset along its normal, a sphere's
/// point less its centre, the others' local coordinates in their frame.
#[derive(Clone, Debug)]
struct Field<T> {
    kind: QKind,
    a: [[T; 3]; 3],
    t: [T; 3],
    /// A plane's `|m|`, a sphere's `r`, a cylinder's `(cu, cv, r)`, a
    /// cone's `(b, k, sqrt(1 + k^2))`, a torus's `(R, r)`.
    c: [T; 3],
    arcs: Vec<PArc<T>>,
    arcs64: Vec<PArc<f64>>,
    hole: bool,
}

impl<T: Num> Field<T> {
    /// Face `fi` of `m` against `P`'s frame `pf`.
    fn of(m: &Prism, fi: usize, pf: &Affine) -> Option<Self> {
        let (vm, vi) = m.view(fi);
        let f = &vm.f;
        let z = || T::k(0.0);
        let axes = [&pf.x, &pf.y, &pf.n];
        // `y = t + a l` from rows `r` applied to `x - o`.
        let compose = |rows: [V; 3], o: &V| -> ([[T; 3]; 3], [T; 3]) {
            let d = sub(&pf.o, o);
            (
                rows.clone()
                    .map(|r| [0, 1, 2].map(|i| T::of_r(&dot(&r, axes[i])))),
                rows.map(|r| T::of_r(&dot(&r, &d))),
            )
        };
        let frame = || compose([f.row(0).clone(), f.row(1).clone(), f.row(2).clone()], &f.o);
        let none = (Vec::new(), Vec::new(), false);
        let zero3 = || [zero(), zero(), zero()];
        let ((a, t), kind, c, (arcs, arcs64, hole)) = match &vm.faces[vi].surf {
            Surf::Plane { p, m } => (
                compose([m.clone(), zero3(), zero3()], p),
                QKind::Plane,
                [T::of_r(&dot(m, m)).sqrt(), z(), z()],
                none,
            ),
            Surf::Sphere { c, r } => {
                let unit = |k: usize| {
                    let mut e = zero3();
                    e[k] = int(1);
                    e
                };
                (
                    compose([unit(0), unit(1), unit(2)], c),
                    QKind::Sphere,
                    [T::of_r(r), z(), z()],
                    none,
                )
            }
            Surf::Cyl { c, r, .. } => (
                frame(),
                QKind::Cyl,
                [T::of_r(&c[0]), T::of_r(&c[1]), T::of_r(r)],
                none,
            ),
            Surf::Cone { b, k } => (
                frame(),
                QKind::Cone,
                [T::of_r(b), T::of_r(k), T::of_r(&(int(1) + k * k)).sqrt()],
                none,
            ),
            Surf::Torus => {
                let ring = vm.ring.as_ref()?;
                (
                    frame(),
                    QKind::Torus,
                    [T::of_r(&ring.big), T::of_r(&ring.small), z()],
                    none,
                )
            }
            Surf::Spline(seg) => (
                frame(),
                QKind::Spline,
                [z(), z(), z()],
                (parcs_of(seg), parcs(seg), seg.hole),
            ),
        };
        Some(Self {
            kind,
            a,
            t,
            c,
            arcs,
            arcs64,
            hole,
        })
    }

    /// The signed distance at `P`'s local coordinates `l`.
    fn dist(&self, l: &[T; 3]) -> T {
        let y = [0, 1, 2].map(|k| self.t[k].add(&dot3(&self.a[k], l)));
        match self.kind {
            QKind::Plane => y[0].div(&self.c[0]),
            QKind::Sphere => y[0]
                .sq()
                .add(&y[1].sq())
                .add(&y[2].sq())
                .sqrt()
                .sub(&self.c[0]),
            QKind::Cyl => {
                let (du, dv) = (y[0].sub(&self.c[0]), y[1].sub(&self.c[1]));
                du.sq().add(&dv.sq()).sqrt().sub(&self.c[2])
            }
            QKind::Cone => {
                let rho = y[0].sq().add(&y[1].sq()).sqrt();
                let rad = self.c[0].add(&self.c[1].mul(&y[2]));
                rho.sub(&rad).div(&self.c[2])
            }
            QKind::Torus => {
                let rho = y[0].sq().add(&y[1].sq()).sqrt();
                let e = rho.sub(&self.c[0]);
                e.sq().add(&y[2].sq()).sqrt().sub(&self.c[1])
            }
            QKind::Spline => {
                let Some(t0) = near_run(&self.arcs64, [y[0].mid(), y[1].mid()]) else {
                    return T::k(f64::NAN);
                };
                // Two Newton steps in this arithmetic from the binary64
                // foot carry its derivatives.
                let mut tau = T::k(t0);
                for _ in 0..2 {
                    let [c, d1, d2] = profile_jet(&self.arcs, &tau);
                    let e = [y[0].sub(&c[0]), y[1].sub(&c[1])];
                    let g = e[0].mul(&d1[0]).add(&e[1].mul(&d1[1]));
                    let h = e[0]
                        .mul(&d2[0])
                        .add(&e[1].mul(&d2[1]))
                        .sub(&d1[0].sq())
                        .sub(&d1[1].sq());
                    tau = tau.sub(&g.div(&h));
                }
                let [c, d1, _] = profile_jet(&self.arcs, &tau);
                let e = [y[0].sub(&c[0]), y[1].sub(&c[1])];
                let len = d1[0].sq().add(&d1[1].sq()).sqrt();
                let right = e[0].mul(&d1[1]).sub(&e[1].mul(&d1[0])).div(&len);
                if self.hole {
                    T::k(0.0).sub(&right)
                } else {
                    right
                }
            }
        }
    }
}

/// A face's rank as `Q`: a plane, a sphere, a cylinder, a cone, a torus, a
/// spline wall.
fn rank(s: &Surf) -> usize {
    match s {
        Surf::Plane { .. } => 0,
        Surf::Sphere { .. } => 1,
        Surf::Cyl { .. } => 2,
        Surf::Cone { .. } => 3,
        Surf::Torus => 4,
        Surf::Spline(_) => 5,
    }
}

/// Whether the face pairs' rule takes two surfaces of these ranks, and its
/// reason: two cylinders or cones, or a plane and one; a torus and anything
/// but a spline wall; a spline wall and anything but a torus.
fn taken(a: usize, b: usize) -> Option<&'static str> {
    let (lo, hi) = (a.min(b), a.max(b));
    match hi {
        2 | 3 if lo != 1 => Some(NEAR_WALLS),
        4 => Some(NEAR_TORUS),
        5 if lo != 4 => Some(NEAR_SPLINE),
        _ => None,
    }
}

/// Faces of the two inputs within the resolution of touching at a point,
/// along a line or along a circle (REVIEW_NOTES.md, "Faces with no near-miss
/// rule of their own"): two cylinders' or cones' walls, a plane and one, a
/// torus or a spline wall against another surface. One face (`P`) is
/// sampled in its parameters, the other (`Q`) held by its signed distance
/// `phi`; each local least or greatest of `phi` over `P` within twice the
/// resolution in binary64 (Newton's method on forward jets) is decided at
/// exact points of both faces (`face_pair`): within the resolution (their
/// exact gap at most it, above the gap's own rounding), both in their
/// faces, a crossing or a miss with the gap outside either input is
/// `Degenerate`; beyond it, apart where Krawczyk's operator certifies it,
/// else a `ComputationLimit`. Two rods crossing whose walls missed by
/// `1e-12`, a ball by a torus, a rod by a spline wall were fused into two
/// solids.
pub(super) fn face_pairs(models: &[Prism; 2], res: f64) -> Result<()> {
    let wide = |b: &([f64; 3], [f64; 3])| (b.0.map(|x| x - 2.0 * res), b.1.map(|x| x + 2.0 * res));
    for fa in 0..models[0].faces.len() {
        let (va, ia) = models[0].view(fa);
        let ra = rank(&va.faces[ia].surf);
        for fb in 0..models[1].faces.len() {
            let (vb, ib) = models[1].view(fb);
            let rb = rank(&vb.faces[ib].surf);
            let Some(why) = taken(ra, rb) else {
                continue;
            };
            if !boxes_meet(&wide(&models[0].boxes[fa]), &models[1].boxes[fb]) {
                continue;
            }
            // `P` the higher rank, the first input's on a tie.
            let (ip, fp, fq) = if ra >= rb { (0, fa, fb) } else { (1, fb, fa) };
            face_pair(models, ip, fp, fq, res, why)?;
        }
    }
    Ok(())
}

/// A critical point of `phi` over `P`'s parameters (binary64): its
/// parameters, value, gradient and Hessian, and the point's derivatives.
#[derive(Clone, Debug)]
struct Crit {
    p: [f64; 2],
    phi: f64,
    h: [f64; 3],
    dx: [[f64; 3]; 2],
}

/// The kind of a critical point: a point least or greatest (`Less`,
/// `Greater`), or one along a curve of them (the null direction in the
/// parameters); `None` a saddle or coincident surfaces.
#[derive(Clone, Copy, Debug)]
enum Contact {
    Point(Ordering),
    Curve(Ordering, [f64; 2]),
}

/// The parameter grid of a face pair: `P`'s parameters' ranges and counts,
/// whether each is periodic, and `P`'s speed bounds in each.
struct Grid {
    lo: [f64; 2],
    hi: [f64; 2],
    n: [usize; 2],
    periodic: [bool; 2],
    speed: [f64; 2],
}

impl Grid {
    fn step(&self, k: usize) -> f64 {
        (self.hi[k] - self.lo[k]) / self.n[k] as f64
    }

    fn at(&self, i: usize, j: usize) -> [f64; 2] {
        [
            self.lo[0] + self.step(0) * i as f64,
            self.lo[1] + self.step(1) * j as f64,
        ]
    }

    /// Whether parameters lie within the grid's domain (a periodic one
    /// always), widened by a cell.
    fn holds(&self, p: &[f64; 2]) -> bool {
        (0..2).all(|k| {
            self.periodic[k]
                || (self.lo[k] - self.step(k) <= p[k] && p[k] <= self.hi[k] + self.step(k))
        })
    }

    fn wrap(&self, p: [f64; 2]) -> [f64; 2] {
        [0, 1].map(|k| {
            if self.periodic[k] {
                let w = (p[k] - self.lo[k]).rem_euclid(TAU);
                self.lo[k] + w
            } else {
                p[k]
            }
        })
    }
}

/// Whether `Q`'s surface lies beyond `reach` of every point of a box
/// (binary64, with a rounding allowance of `1e-12` of the box's size): a
/// plane's offsets at the box's corners all of one sign beyond it (an
/// affine function's extremes); a sphere's distances from its centre to
/// the box against its radius; a cylinder's, a cone's or a torus's
/// distance over the box of its local coordinates holding the box's image
/// (its corners' extents). A spline wall is never passed over so.
fn box_apart(m: &Prism, fi: usize, b: &([f64; 3], [f64; 3]), reach: f64) -> bool {
    let (vm, vi) = m.view(fi);
    let f = &vm.f;
    let corners: Vec<[f64; 3]> = (0..8)
        .map(|c| std::array::from_fn(|k| if c >> k & 1 == 1 { b.1[k] } else { b.0[k] }))
        .collect();
    let size = corners
        .iter()
        .flatten()
        .fold(1.0f64, |acc, x| acc.max(x.abs()));
    let reach = reach + 1e-12 * size;
    let g = |v: &V| [0, 1, 2].map(|k| rational_f64(&v[k]));
    // Beyond `reach` either way of an enclosed signed distance.
    let apart = |lo: f64, hi: f64| lo > reach || hi < -reach;
    // The box of the corners' local coordinates.
    let local = || {
        let (o, rows) = (g(&f.o), [g(f.row(0)), g(f.row(1)), g(f.row(2))]);
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for c in &corners {
            for k in 0..3 {
                let l: f64 = (0..3).map(|j| rows[k][j] * (c[j] - o[j])).sum();
                lo[k] = lo[k].min(l);
                hi[k] = hi[k].max(l);
            }
        }
        (lo, hi)
    };
    // The least and greatest distance from the origin over a 2D box.
    let radial = |lo: [f64; 2], hi: [f64; 2]| {
        let near = [0, 1].map(|k| 0.0f64.clamp(lo[k], hi[k]));
        let far = [0, 1].map(|k| lo[k].abs().max(hi[k].abs()));
        (near[0].hypot_rn(near[1]), far[0].hypot_rn(far[1]))
    };
    match &vm.faces[vi].surf {
        Surf::Plane { p, m } => {
            let (p, m) = (g(p), g(m));
            let len = m[0].hypot_rn(m[1]).hypot_rn(m[2]);
            let ds: Vec<f64> = corners
                .iter()
                .map(|c| (0..3).map(|k| m[k] * (c[k] - p[k])).sum::<f64>() / len)
                .collect();
            let (lo, hi) = ds
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
                    (a.min(x), b.max(x))
                });
            apart(lo, hi)
        }
        Surf::Sphere { c, r } => {
            let (c, r) = (g(c), rational_f64(r));
            let near: f64 = (0..3)
                .map(|k| (c[k].clamp(b.0[k], b.1[k]) - c[k]).powi(2))
                .sum::<f64>()
                .sqrt();
            let far: f64 = (0..3)
                .map(|k| (c[k] - b.0[k]).abs().max((c[k] - b.1[k]).abs()).powi(2))
                .sum::<f64>()
                .sqrt();
            apart(near - r, far - r)
        }
        Surf::Cyl { c, r, .. } => {
            let (lo, hi) = local();
            let (cu, cv, r) = (rational_f64(&c[0]), rational_f64(&c[1]), rational_f64(r));
            let (near, far) = radial([lo[0] - cu, lo[1] - cv], [hi[0] - cu, hi[1] - cv]);
            apart(near - r, far - r)
        }
        Surf::Cone { b: cb, k } => {
            let (lo, hi) = local();
            let (bf, kf) = (rational_f64(cb), rational_f64(k));
            let (near, far) = radial([lo[0], lo[1]], [hi[0], hi[1]]);
            let rads = [bf + kf * lo[2], bf + kf * hi[2]];
            let (rl, rh) = (rads[0].min(rads[1]), rads[0].max(rads[1]));
            let s = (1.0 + kf * kf).sqrt();
            apart((near - rh) / s, (far - rl) / s)
        }
        Surf::Torus => {
            let Some(ring) = vm.ring.as_ref() else {
                return false;
            };
            let (lo, hi) = local();
            let (big, small) = (rational_f64(&ring.big), rational_f64(&ring.small));
            let (rn, rf) = radial([lo[0], lo[1]], [hi[0], hi[1]]);
            let (near, far) = radial([rn - big, lo[2]], [rf - big, hi[2]]);
            apart(near - small, far - small)
        }
        Surf::Spline(_) => false,
    }
}

/// The local heights in `P`'s frame over a box (its eight corners).
fn heights(o: &[f64; 3], rows: &[[f64; 3]; 3], b: &([f64; 3], [f64; 3])) -> [f64; 2] {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for m in 0..8 {
        let c: [f64; 3] = std::array::from_fn(|k| if m >> k & 1 == 1 { b.1[k] } else { b.0[k] });
        let w: f64 = (0..3).map(|k| rows[2][k] * (c[k] - o[k])).sum();
        lo = lo.min(w);
        hi = hi.max(w);
    }
    [lo, hi]
}

/// `P`'s grid over its face's box met by `Q`'s widened by twice the
/// resolution: an angle's turn in 32 cells, heights and runs in cells of
/// about the angle's cells' length (at most 32, a height range a cell
/// wider each way), and `P`'s speed bounds (its radii or its control
/// polygon, a ten-thousandth for its frame's axes' lengths).
fn grid_of(m: &Prism, fi: usize, b: &([f64; 3], [f64; 3])) -> Option<Grid> {
    let (vm, vi) = m.view(fi);
    let f = &vm.f;
    let g = |v: &V| [0, 1, 2].map(|k| rational_f64(&v[k]));
    let (o, rows) = (g(&f.o), [g(f.row(0)), g(f.row(1)), g(f.row(2))]);
    let slack = 1.0001;
    let turn = TAU / 32.0;
    let span = |lo: f64, hi: f64, cell: f64| -> (f64, f64, usize) {
        let n = (((hi - lo) / cell.max(1e-300)).ceil() as usize).clamp(1, 30);
        let c = (hi - lo) / n as f64;
        (lo - c, hi + c, n + 2)
    };
    let full = (-PI, PI);
    match &vm.faces[vi].surf {
        Surf::Cyl { r, .. } => {
            let r = rational_f64(r).abs();
            let [wl, wh] = heights(&o, &rows, b);
            let (lo, hi, n) = span(wl, wh, r * turn);
            Some(Grid {
                lo: [full.0, lo],
                hi: [full.1, hi],
                n: [32, n],
                periodic: [true, false],
                speed: [r * slack, slack],
            })
        }
        Surf::Cone { b: cb, k } => {
            let (bf, kf) = (rational_f64(cb), rational_f64(k));
            let [mut wl, mut wh] = heights(&o, &rows, b);
            // On the nappe's side of the apex, clear of it.
            if kf != 0.0 {
                let apex = -bf / kf;
                let clear = 1e-9 * (1.0 + apex.abs());
                if kf > 0.0 {
                    wl = wl.max(apex + clear);
                } else {
                    wh = wh.min(apex - clear);
                }
            }
            if !above(wh, wl) {
                return None;
            }
            let rad = (bf + kf * wl).abs().max((bf + kf * wh).abs());
            let (lo, hi, n) = span(wl, wh, rad * turn);
            // The heights a cell wider stay off the apex.
            let (lo, hi) = if kf > 0.0 {
                (lo.max(-bf / kf + 1e-9 * (1.0 + (bf / kf).abs())), hi)
            } else if kf < 0.0 {
                (lo, hi.min(-bf / kf - 1e-9 * (1.0 + (bf / kf).abs())))
            } else {
                (lo, hi)
            };
            let rad = rad.max((bf + kf * lo).abs()).max((bf + kf * hi).abs());
            Some(Grid {
                lo: [full.0, lo],
                hi: [full.1, hi],
                n: [32, n],
                periodic: [true, false],
                speed: [rad * slack, (1.0 + kf * kf).sqrt() * slack],
            })
        }
        Surf::Torus => {
            let ring = vm.ring.as_ref()?;
            let (big, small) = (rational_f64(&ring.big), rational_f64(&ring.small));
            Some(Grid {
                lo: [full.0, full.0],
                hi: [full.1, full.1],
                n: [32, 32],
                periodic: [true, true],
                speed: [(big + small.abs()) * slack, small.abs() * slack],
            })
        }
        Surf::Spline(seg) => {
            // The run's speed from each arc's control polygon (a Bézier
            // derivative's control points hold it).
            let mut speed = 0.0f64;
            for a in &seg.arcs {
                let deg = (a.cps.len() - 1) as f64;
                let len = rational_f64(&(&a.d[1] - &a.d[0]));
                for w in a.cps.windows(2) {
                    let du = rational_f64(&(&w[1][0] - &w[0][0]));
                    let dv = rational_f64(&(&w[1][1] - &w[0][1]));
                    speed = speed.max(deg * du.hypot_rn(dv) / len);
                }
            }
            let (first, last) = (rational_f64(&seg.first), rational_f64(&seg.last));
            let nt = (16 * seg.arcs.len()).min(64);
            let cell = speed * (last - first) / nt as f64;
            let [wl, wh] = heights(&o, &rows, b);
            let (lo, hi, n) = span(wl, wh, cell);
            Some(Grid {
                lo: [first, lo],
                hi: [last, hi],
                n: [nt, n],
                periodic: [false, false],
                speed: [speed * slack, slack],
            })
        }
        _ => None,
    }
}

/// `phi` and its derivatives at parameters (binary64 jets), with the
/// point's first derivatives.
fn jet_at(pp: &Patch<J<f64>>, qf: &Field<J<f64>>, p: [f64; 2]) -> (J<f64>, [[f64; 3]; 2]) {
    let vars = [J::var(p[0], 0), J::var(p[1], 1)];
    let x = pp.point(&vars);
    let dx = [0, 1].map(|i| [x[0].d[i], x[1].d[i], x[2].d[i]]);
    (qf.dist(&pp.local(&vars)), dx)
}

/// Newton's method on `phi`'s gradient from parameters `p` (steps within a
/// cell, directions of an eigenvalue below `10^-9` of the larger left
/// alone): the critical point where it converges within the grid.
fn newton(pp: &Patch<J<f64>>, qf: &Field<J<f64>>, grid: &Grid, p: [f64; 2]) -> Option<Crit> {
    let norm = |j: &J<f64>| j.d[0] * j.d[0] + j.d[1] * j.d[1];
    let crit = |p: [f64; 2], (j, dx): (J<f64>, [[f64; 3]; 2])| Crit {
        p,
        phi: j.v,
        h: [j.h[0], j.h[1], j.h[2]],
        dx,
    };
    let mut p = p;
    let mut cur = jet_at(pp, qf, p);
    for _ in 0..60 {
        let j = &cur.0;
        let (g, h) = ([j.d[0], j.d[1]], [j.h[0], j.h[1], j.h[2]]);
        if !(j.v.is_finite() && g.iter().chain(&h).all(|x| x.is_finite())) {
            return None;
        }
        let (vals, vecs) = sym_eigen(h);
        let big = vals[0].abs().max(vals[1].abs());
        let mut step = [0.0f64; 2];
        for k in 0..2 {
            if vals[k].abs() <= 1e-9 * big || big == 0.0 {
                continue;
            }
            let c = -(vecs[k][0] * g[0] + vecs[k][1] * g[1]) / vals[k];
            step[0] += c * vecs[k][0];
            step[1] += c * vecs[k][1];
        }
        // Within a cell each way.
        let scale = (0..2)
            .map(|k| step[k].abs() / grid.step(k))
            .fold(0.0f64, f64::max);
        if scale > 1.0 {
            step = step.map(|s| s / scale);
        }
        let tiny = (0..2).all(|k| step[k].abs() <= 1e-14 * (1.0 + p[k].abs()));
        if tiny {
            return Some(crit(p, cur));
        }
        // Halved until the gradient shrinks; at rounding's level, where no
        // step shrinks it, the point is taken as found when the step is
        // well within a cell.
        let gn = norm(j);
        let mut t = 1.0;
        let mut next = None;
        for _ in 0..40 {
            let q = grid.wrap([p[0] + t * step[0], p[1] + t * step[1]]);
            if grid.holds(&q) {
                let nj = jet_at(pp, qf, q);
                if norm(&nj.0) < gn {
                    next = Some((q, nj));
                    break;
                }
            }
            t *= 0.5;
        }
        match next {
            Some((q, nj)) => {
                p = q;
                cur = nj;
            }
            None => {
                let small = (0..2).all(|k| step[k].abs() <= 1e-6 * grid.step(k));
                return small.then(|| crit(p, cur));
            }
        }
    }
    None
}

/// A symmetric 2x2 matrix's eigenvalues and unit eigenvectors.
fn sym_eigen(h: [f64; 3]) -> ([f64; 2], [[f64; 2]; 2]) {
    let (a, b, c) = (h[0], h[1], h[2]);
    let m = 0.5 * (a + c);
    let r = (0.5 * (a - c)).hypot_rn(b);
    let vals = [m - r, m + r];
    let vec = |l: f64| -> [f64; 2] {
        let (x, y) = if (a - l).abs() + b.abs() >= (c - l).abs() + b.abs() {
            (-b, a - l)
        } else {
            (c - l, -b)
        };
        let n = x.hypot_rn(y);
        if n > 0.0 {
            [x / n, y / n]
        } else {
            [1.0, 0.0]
        }
    };
    let v0 = if b == 0.0 {
        if a <= c {
            [1.0, 0.0]
        } else {
            [0.0, 1.0]
        }
    } else {
        vec(vals[0])
    };
    (vals, [v0, [-v0[1], v0[0]]])
}

/// The kind of a critical point by its Hessian in `P`'s tangent metric
/// (its generalized eigenvalues, the faces' curvatures' differences across
/// the gap): both of one sign beyond `10^-9` of the larger, a point least
/// or greatest; one within `10^-6` of the larger, a curve of them along its
/// eigenvector; of both signs, a saddle (`None`); both nil, coincident
/// surfaces (`None`).
fn contact_of(c: &Crit, size: f64) -> Option<Contact> {
    let g = [
        dot64(&c.dx[0], &c.dx[0]),
        dot64(&c.dx[0], &c.dx[1]),
        dot64(&c.dx[1], &c.dx[1]),
    ];
    let det = g[0] * g[2] - g[1] * g[1];
    if !above(det, 0.0) {
        return None;
    }
    // det(H - l G) = 0: A l^2 - B l + C = 0.
    let (h0, h1, h2) = (c.h[0], c.h[1], c.h[2]);
    let qa = det;
    let qb = h0 * g[2] + h2 * g[0] - 2.0 * h1 * g[1];
    let qc = h0 * h2 - h1 * h1;
    let disc = (qb * qb - 4.0 * qa * qc).max(0.0).sqrt();
    let (l1, l2) = ((qb - disc) / (2.0 * qa), (qb + disc) / (2.0 * qa));
    let (small, large) = if l1.abs() <= l2.abs() {
        (l1, l2)
    } else {
        (l2, l1)
    };
    if !above(large.abs(), 1e-9 / size.max(1e-300)) {
        return None;
    }
    let kind = |l: f64| {
        if l > 0.0 {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    };
    if small.abs() <= 1e-6 * large.abs() {
        // The null direction: (H - small G) v = 0.
        let (a, b) = (h0 - small * g[0], h1 - small * g[1]);
        let (cc, d) = (h1 - small * g[1], h2 - small * g[2]);
        let v = if a.abs() + b.abs() >= cc.abs() + d.abs() {
            [-b, a]
        } else {
            [-d, cc]
        };
        let n = v[0].hypot_rn(v[1]);
        if !above(n, 0.0) {
            return None;
        }
        return Some(Contact::Curve(kind(large), [v[0] / n, v[1] / n]));
    }
    if small.abs() <= 1e-9 * large.abs() || (small > 0.0) != (large > 0.0) {
        return None;
    }
    Some(Contact::Point(kind(large)))
}

/// The exact point of `P` at rational parameters near `p` (half-angle
/// tangents of its angles, rational heights and run parameters).
fn patch_point(m: &Prism, fi: usize, p: [f64; 2]) -> Option<QV> {
    let (vm, vi) = m.view(fi);
    let f = &vm.f;
    let x = match &vm.faces[vi].surf {
        Surf::Cyl { c, r, .. } => {
            let [cs, sn] = unit_at(p[0]);
            f.point(&(&c[0] + r * &cs), &(&c[1] + r * &sn), &q(p[1]))
        }
        Surf::Cone { b, k } => {
            let [cs, sn] = unit_at(p[0]);
            let w = q(p[1]);
            let rad = b + k * &w;
            if sign(&rad) != Ordering::Greater {
                return None;
            }
            f.point(&(&rad * &cs), &(&rad * &sn), &w)
        }
        Surf::Torus => {
            let ring = vm.ring.as_ref()?;
            let ([cu, su], [cv, sv]) = (unit_at(p[0]), unit_at(p[1]));
            let rad = &ring.big + &ring.small * &cv;
            f.point(&(&rad * &cu), &(&rad * &su), &(&ring.small * &sv))
        }
        Surf::Spline(seg) => {
            let uv = seg.point(&Qd::rat(run_rational(seg, p[0])));
            return Some(super::spline_walls::qpoint(
                f,
                &uv[0],
                &uv[1],
                &Qd::rat(q(p[1])),
            ));
        }
        _ => return None,
    };
    Some(qv(&x))
}

/// `Q`'s exact point near the foot of `a` and `Q`'s function's sign at `a`
/// (`Quad`'s for a plane, a cylinder, a cone, a torus or a spline wall; a
/// rational point of a sphere toward `a`, `|a - c|^2 - r^2`).
fn field_foot(m: &Prism, fi: usize, a: &QV) -> Option<(QV, Ordering)> {
    let (vm, vi) = m.view(fi);
    if let Surf::Sphere { c, r } = &vm.faces[vi].surf {
        let w = qsub(a, &qv(c));
        let s = qv(&sphere_toward(c, r, qv_f64(&w)));
        return Some((s, qqdot(&w, &w).add_r(&-(r * r)).sign()));
    }
    let quad = Quad::of(m, fi).or_else(|| Quad::of_curved(m, fi))?;
    let s = quad.foot(a)?;
    let side = match &quad.kind {
        // Its outward offset at `s` (`value` would find `s` again).
        Kind::Spline { m, fi, .. } => qqdot(&qsub(a, &s), &m.normal_at(*fi, &s)).sign(),
        _ => quad.value(a).sign(),
    };
    Some((s, side))
}

/// What the exact points at a critical point decide: no contact of the
/// rule's own, or a contact beyond the resolution (within twice it in
/// binary64: to be certified apart); within it, the refusal.
enum Decided {
    None,
    Beyond,
}

/// The decision at `P`'s parameters `p` of a critical point of the kind
/// `kind`: `a` exactly on `P` near them, `s` on `Q` near `a`'s foot; a gap
/// of rounding's size, a point outside its face or `a` on `Q` is no
/// contact; a least below `Q` or a greatest above it a crossing, otherwise
/// a miss, a contact unless the gap lies inside both inputs; a contact
/// within the resolution `Degenerate`. `Ok(None)` (with `out` set) where a
/// point lies outside its face.
#[allow(clippy::too_many_arguments)]
fn decide(
    models: &[Prism; 2],
    ip: usize,
    fp: usize,
    fq: usize,
    p: [f64; 2],
    kind: Ordering,
    res: f64,
    why: &'static str,
    out: &mut bool,
) -> Result<Decided> {
    let (mp, mq) = (&models[ip], &models[1 - ip]);
    let Some(a) = patch_point(mp, fp, p) else {
        return Ok(Decided::None);
    };
    let Some((s, side)) = field_foot(mq, fq, &a) else {
        return Ok(Decided::None);
    };
    let gap = qsub(&s, &a);
    let gg = qqdot(&gap, &gap);
    if side == Ordering::Equal || gg.sub(&Qd::rat(gap_rounding(&a, &s))).sign() != Ordering::Greater
    {
        // On `Q`, or a gap of rounding's size: the faces' own rules.
        return Ok(Decided::None);
    }
    if mp.in_face(fp, &a) == Loc::Out || mq.in_face(fq, &s) == Loc::Out {
        *out = true;
        return Ok(Decided::None);
    }
    if side != kind {
        // Missing: no contact where the gap lies inside both inputs.
        let inside = [
            mp.member(&a, std::slice::from_ref(&gap)) == Loc::In,
            mq.member(&s, &[gap.map(|x| x.neg())]) == Loc::In,
        ];
        if inside == [true, true] {
            return Ok(Decided::None);
        }
    }
    let tol = q(res);
    if gg.add_r(&-(&tol * &tol)).sign() != Ordering::Greater {
        return Err(Error::Degenerate(why));
    }
    Ok(Decided::Beyond)
}

/// The face pairs' rule for `P` (face `fp` of input `ip`) and `Q` (face
/// `fq` of the other).
fn face_pair(
    models: &[Prism; 2],
    ip: usize,
    fp: usize,
    fq: usize,
    res: f64,
    why: &'static str,
) -> Result<()> {
    let (mp, mq) = (&models[ip], &models[1 - ip]);
    let bq = &mq.boxes[fq];
    let bp = &mp.boxes[fp];
    let lo: [f64; 3] = std::array::from_fn(|k| bp.0[k].max(bq.0[k] - 2.0 * res));
    let hi: [f64; 3] = std::array::from_fn(|k| bp.1[k].min(bq.1[k] + 2.0 * res));
    if (0..3).any(|k| lo[k] > hi[k]) {
        return Ok(());
    }
    let b = (lo, hi);
    if box_apart(mq, fq, &b, 2.0 * res) {
        return Ok(());
    }
    let Some(grid) = grid_of(mp, fp, &b) else {
        return Ok(());
    };
    let (Some(pp), Some(qf)) = (
        Patch::<f64>::of(mp, fp),
        Field::<f64>::of(mq, fq, &mp.view(fp).0.f),
    ) else {
        return Ok(());
    };
    // The grid's nodes: `phi` and the point (a periodic parameter's last
    // node its first).
    let count = |k: usize| grid.n[k] + usize::from(!grid.periodic[k]);
    let (n0, n1) = (count(0), count(1));
    let mut phi = vec![0.0f64; n0 * n1];
    let mut far = vec![false; n0 * n1];
    let size = b.1.iter().chain(&b.0).fold(1.0f64, |m, x| m.max(x.abs()));
    let cell = 0.5 * (grid.speed[0] * grid.step(0) + grid.speed[1] * grid.step(1));
    // A rounding allowance for the binary64 distances.
    let reach = 2.0 * res + cell + 1e-12 * size;
    let checked = qf.kind != QKind::Spline;
    let mut all_far = true;
    for i in 0..n0 {
        for j in 0..n1 {
            let at = grid.at(i, j);
            let x = pp.point(&at);
            let v = qf.dist(&pp.local(&at));
            let outside = (0..3).any(|k| x[k] < lo[k] - cell || x[k] > hi[k] + cell);
            let f = checked && (outside || v.abs() > reach);
            phi[i * n1 + j] = v;
            far[i * n1 + j] = f;
            all_far &= f;
        }
    }
    if all_far {
        // No point of `P` within twice the resolution of `Q` in the box.
        return Ok(());
    }
    // Seeds: interior nodes, not far, each a least or a greatest of `phi`
    // among its eight neighbours.
    let index = |i: isize, k: usize, n: usize| -> Option<usize> {
        if grid.periodic[k] {
            Some(i.rem_euclid(n as isize) as usize)
        } else if i < 0 || i >= n as isize {
            None
        } else {
            Some(i as usize)
        }
    };
    let at =
        |i: isize, j: isize| -> Option<f64> { Some(phi[index(i, 0, n0)? * n1 + index(j, 1, n1)?]) };
    let inner = |i: usize, k: usize, n: usize| grid.periodic[k] || (i > 0 && i + 1 < n);
    let mut seeds: Vec<(f64, [f64; 2])> = Vec::new();
    for i in 0..n0 {
        if !inner(i, 0, n0) {
            continue;
        }
        for j in 0..n1 {
            if !inner(j, 1, n1) {
                continue;
            }
            let v = phi[i * n1 + j];
            if far[i * n1 + j] || !v.is_finite() {
                continue;
            }
            let (mut least, mut most) = (true, true);
            for di in -1..=1isize {
                for dj in -1..=1isize {
                    if di == 0 && dj == 0 {
                        continue;
                    }
                    if let Some(w) = at(i as isize + di, j as isize + dj) {
                        least &= v <= w;
                        most &= v >= w;
                    }
                }
            }
            if least || most {
                seeds.push((v.abs(), grid.at(i, j)));
            }
        }
    }
    if seeds.is_empty() {
        return Ok(());
    }
    seeds.sort_by(|a, b| a.0.total_cmp(&b.0));
    seeds.truncate(16);
    let (Some(pj), Some(qj)) = (
        Patch::<J<f64>>::of(mp, fp),
        Field::<J<f64>>::of(mq, fq, &mp.view(fp).0.f),
    ) else {
        return Ok(());
    };
    let mut crits: Vec<Crit> = Vec::new();
    for (_, s) in seeds {
        let Some(c) = newton(&pj, &qj, &grid, s) else {
            continue;
        };
        if !at_most(c.phi.abs(), 2.0 * res + 1e-12 * size) {
            continue;
        }
        let near = |d: &Crit| {
            (0..2).all(|k| {
                let mut e = (d.p[k] - c.p[k]).abs();
                if grid.periodic[k] {
                    e = e.min(TAU - e);
                }
                e <= 0.5 * grid.step(k)
            })
        };
        if crits.iter().any(near) {
            continue;
        }
        crits.push(c);
        if crits.len() == 8 {
            break;
        }
    }
    // A curve of critical points is decided once.
    let mut curves: Vec<f64> = Vec::new();
    for c in crits {
        let Some(contact) = contact_of(&c, size) else {
            continue;
        };
        let (kind, null) = match contact {
            Contact::Point(k) => (k, None),
            Contact::Curve(k, v) => {
                if curves.iter().any(|f| (f - c.phi).abs() <= 1e-12 * size) {
                    continue;
                }
                curves.push(c.phi);
                (k, Some(v))
            }
        };
        let mut out = false;
        let mut decided = decide(models, ip, fp, fq, c.p, kind, res, why, &mut out)?;
        let mut at = c.clone();
        if out {
            // Along a curve of critical points: another of its points in
            // both faces (Newton's method again from steps along it).
            let Some(v) = null else {
                continue;
            };
            let Some(next) =
                along_curve(models, ip, fp, fq, &pj, &qj, &grid, &c, v, kind, res, why)?
            else {
                continue;
            };
            (decided, at) = next;
        }
        if let Decided::Beyond = decided {
            if !certified_apart(mp, fp, mq, fq, &grid, &at, null, res) {
                return Err(Error::ComputationLimit(NEAR_LIMIT));
            }
        }
    }
    Ok(())
}

/// A curve of critical points' point inside both faces: Newton's method
/// from steps along its null direction (at most seven, the middle of the
/// grid's range along it first), each decided as `decide`; `None` where no
/// point of it lies in both faces.
#[allow(clippy::too_many_arguments)]
fn along_curve(
    models: &[Prism; 2],
    ip: usize,
    fp: usize,
    fq: usize,
    pj: &Patch<J<f64>>,
    qj: &Field<J<f64>>,
    grid: &Grid,
    c: &Crit,
    v: [f64; 2],
    kind: Ordering,
    res: f64,
    why: &'static str,
) -> Result<Option<(Decided, Crit)>> {
    // The range of steps keeping the parameters in the grid's domain.
    let mut t = [-f64::INFINITY, f64::INFINITY];
    for (k, &vk) in v.iter().enumerate() {
        if grid.periodic[k] || vk.abs() < 1e-12 {
            if grid.periodic[k] && vk.abs() >= 1e-12 {
                let reach = PI / vk.abs();
                t = [t[0].max(-reach), t[1].min(reach)];
            }
            continue;
        }
        let (a, b) = ((grid.lo[k] - c.p[k]) / vk, (grid.hi[k] - c.p[k]) / vk);
        t = [t[0].max(a.min(b)), t[1].min(a.max(b))];
    }
    if !above(t[1], t[0]) || !t[0].is_finite() || !t[1].is_finite() {
        return Ok(None);
    }
    for frac in [0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875] {
        let s = t[0] + (t[1] - t[0]) * frac;
        let start = grid.wrap([c.p[0] + s * v[0], c.p[1] + s * v[1]]);
        let Some(d) = newton(pj, qj, grid, start) else {
            continue;
        };
        if !at_most(d.phi.abs(), 2.0 * res + 1e-12) {
            continue;
        }
        let mut out = false;
        let decided = decide(models, ip, fp, fq, d.p, kind, res, why, &mut out)?;
        if !out {
            return Ok(Some((decided, d)));
        }
    }
    Ok(None)
}

/// `phi` at a box of parameters in binary64 intervals: its jet.
fn jet_over(pp: &Patch<J<Iv>>, qf: &Field<J<Iv>>, p: [Iv; 2]) -> J<Iv> {
    qf.dist(&pp.local(&[J::var(p[0], 0), J::var(p[1], 1)]))
}

/// Whether a contact beyond the resolution is certified apart: Krawczyk's
/// operator on `phi`'s gradient over a box `N` of parameters about the
/// critical point (both parameters for a point; the transverse one for a
/// curve whose null direction is a parameter's, the other the grid's whole
/// range) maps `N` inside itself (one critical point in it, for each value
/// of the other parameter), and `|phi|` over `N` by the mean-value form
/// exceeds the resolution over `1 - g` (`g` the Frobenius norm of `M^T M -
/// I` for `Q`'s frame's axes `M`: its local distance against the world's).
/// A spline wall as `Q` is not certified.
#[allow(clippy::too_many_arguments)]
fn certified_apart(
    mp: &Prism,
    fp: usize,
    mq: &Prism,
    fq: usize,
    grid: &Grid,
    c: &Crit,
    null: Option<[f64; 2]>,
    res: f64,
) -> bool {
    let (Some(pp), Some(qf)) = (
        Patch::<J<Iv>>::of(mp, fp),
        Field::<J<Iv>>::of(mq, fq, &mp.view(fp).0.f),
    ) else {
        return false;
    };
    if qf.kind == QKind::Spline {
        return false;
    }
    let g = match qf.kind {
        QKind::Plane | QKind::Sphere => 0.0,
        _ => frame_defect(&mq.view(fq).0.f),
    };
    if !above(0.5, g) {
        return false;
    }
    let thr = res / (1.0 - g) * (1.0 + 1e-12);
    // The parameters' half widths: about `1e-5` in length each.
    let delta = [0, 1].map(|k| 1e-5 / grid.speed[k].max(1e-300));
    let around = |k: usize| Iv::span(c.p[k] - delta[k], c.p[k] + delta[k]);
    let beyond = |v: &Iv| {
        let (lo, hi) = v.bounds();
        lo > thr || hi < -thr
    };
    match null {
        None => {
            // Both parameters.
            let n = [around(0), around(1)];
            let at = jet_over(&pp, &qf, [Iv::k(c.p[0]), Iv::k(c.p[1])]);
            let over = jet_over(&pp, &qf, n);
            // Y, the inverse of the binary64 Hessian.
            let det = c.h[0] * c.h[2] - c.h[1] * c.h[1];
            if !above(det.abs(), 0.0) {
                return false;
            }
            let y = [[c.h[2] / det, -c.h[1] / det], [-c.h[1] / det, c.h[0] / det]];
            let hn = [[over.h[0], over.h[1]], [over.h[1], over.h[2]]];
            let off = [n[0].sub(&Iv::k(c.p[0])), n[1].sub(&Iv::k(c.p[1]))];
            for i in 0..2 {
                let mut k = Iv::k(c.p[i]);
                for j in 0..2 {
                    k = k.sub(&Iv::k(y[i][j]).mul(&at.d[j]));
                    let mut m = Iv::k(if i == j { 1.0 } else { 0.0 });
                    for l in 0..2 {
                        m = m.sub(&Iv::k(y[i][l]).mul(&hn[l][j]));
                    }
                    k = k.add(&m.mul(&off[j]));
                }
                if !k.inside(c.p[i] - delta[i], c.p[i] + delta[i]) {
                    return false;
                }
            }
            let bound =
                at.v.add(&over.d[0].mul(&off[0]))
                    .add(&over.d[1].mul(&off[1]));
            beyond(&bound)
        }
        Some(v) => {
            // The null direction a parameter's: the other transverse.
            let (k, j) = if v[1].abs() >= 1e6 * v[0].abs() {
                (1, 0)
            } else if v[0].abs() >= 1e6 * v[1].abs() {
                (0, 1)
            } else {
                return false;
            };
            let whole = Iv::span(grid.lo[k] - grid.step(k), grid.hi[k] + grid.step(k));
            let mut point = [Iv::k(c.p[0]), Iv::k(c.p[1])];
            point[k] = whole;
            let mut boxed = point;
            boxed[j] = around(j);
            let at = jet_over(&pp, &qf, point);
            let over = jet_over(&pp, &qf, boxed);
            let hj = c.h[if j == 0 { 0 } else { 2 }];
            if !above(hj.abs(), 0.0) {
                return false;
            }
            let y = 1.0 / hj;
            let off = boxed[j].sub(&Iv::k(c.p[j]));
            let hh = &over.h[if j == 0 { 0 } else { 2 }];
            let kk = Iv::k(c.p[j])
                .sub(&Iv::k(y).mul(&at.d[j]))
                .add(&Iv::k(1.0).sub(&Iv::k(y).mul(hh)).mul(&off));
            if !kk.inside(c.p[j] - delta[j], c.p[j] + delta[j]) {
                return false;
            }
            beyond(&at.v.add(&over.d[j].mul(&off)))
        }
    }
}

/// The Frobenius norm of `M^T M - I` for a frame's axes `M` (an upper bound
/// in binary64 of the exact value).
fn frame_defect(f: &Affine) -> f64 {
    let ax = [&f.x, &f.y, &f.n];
    let mut sum = zero();
    for i in 0..3 {
        for j in 0..3 {
            let e = dot(ax[i], ax[j]) - if i == j { int(1) } else { zero() };
            sum += &e * &e;
        }
    }
    rational_f64(&sum).sqrt() * (1.0 + 1e-12) + 1e-300
}

// ------------------------------------------- edges and vertices, curved

/// The edges' and vertices' rules (`edge_faces`, `vertex_faces`) against a
/// torus's or a spline wall's faces (REVIEW_NOTES.md, "Faces with no
/// near-miss rule of their own"), after the earlier rules: a box's edge or
/// a cube's corner `1e-12` off a torus's equator or a spline dome's crest
/// was fused with it into two solids.
pub(super) fn curved_edges_and_vertices(models: &[Prism; 2], res: f64) -> Result<()> {
    edge_faces_with(models, res, true)?;
    vertex_faces_with(models, res, true)
}

/// The points of a line edge where a torus's or a spline wall's function
/// along it is least or greatest: each change of its derivative's sign
/// among 64 samples over the edge bisected to rounding, at a rational
/// parameter of the line.
fn line_extrema(em: &Prism, ei: usize, quad: &Quad) -> Vec<(QV, Ordering)> {
    let e = &em.edges[ei];
    let Crv::Line { p, d } = &e.curve else {
        return Vec::new();
    };
    let (pf, df) = (qv_f64(p), d.clone().map(|x| rational_f64(&x)));
    let dd = dot64(&df, &df);
    if !above(dd, 0.0) {
        return Vec::new();
    }
    let t_of = |v: usize| {
        let x = qv_f64(&em.verts[v].p);
        dot64(&sub64(&x, &pf), &df) / dd
    };
    let (t0, t1) = (t_of(e.start), t_of(e.end));
    let (t0, t1) = (t0.min(t1), t0.max(t1));
    let at = |t: f64| [0, 1, 2].map(|k| pf[k] + t * df[k]);
    let slope = |t: f64| {
        let (_, g) = quad.value64(at(t));
        dot64(&g, &df)
    };
    const N: usize = 64;
    let ts: Vec<f64> = (0..=N)
        .map(|i| t0 + (t1 - t0) * i as f64 / N as f64)
        .collect();
    let ds: Vec<f64> = ts.iter().map(|&t| slope(t)).collect();
    let mut out = Vec::new();
    for i in 0..N {
        let (d0, d1) = (ds[i], ds[i + 1]);
        let kind = if d0 < 0.0 && d1 >= 0.0 {
            Ordering::Less
        } else if d0 > 0.0 && d1 <= 0.0 {
            Ordering::Greater
        } else {
            continue;
        };
        let (mut lo, mut hi) = (ts[i], ts[i + 1]);
        for _ in 0..80 {
            let mid = 0.5 * lo + 0.5 * hi;
            if !(above(mid, lo) && above(hi, mid)) {
                break;
            }
            if (slope(mid) < 0.0) == (kind == Ordering::Less) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let t = Qd::rat(q(0.5 * lo + 0.5 * hi));
        out.push((qadd(p, &qscale(d, &t)), kind));
    }
    out
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

    /// `rounding_gap` never passes over a gap its exact points exceed `rho`
    /// with: conics of turned frames' kinds (centres up to 3 from the
    /// origin, axes of lengths 0.5 to 2 in turned planes), one shifted by
    /// `k rho` along a direction (`k` from 0 to 2) and its data rounded
    /// again, at the same rational parameters; each binary64 verdict of
    /// rounding's size is one the exact points' gap confirms, and most of the
    /// zero gaps (the same conic) are passed over in binary64 (the others
    /// exactly).
    #[test]
    fn rounding_gaps_in_binary64_hold_exactly() {
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let (mut passed, mut tiny) = (0, 0);
        for i in 0..400 {
            let c = [0, 1, 2].map(|_| 6.0 * next() - 3.0);
            let (th, ph, ps) = (TAU * next(), PI * next(), TAU * next());
            let n = [
                ph.sin_rn() * th.cos_rn(),
                ph.sin_rn() * th.sin_rn(),
                ph.cos_rn(),
            ];
            let e1 = {
                let h = if n[0].abs() < 0.9 {
                    [1.0, 0.0, 0.0]
                } else {
                    [0.0, 1.0, 0.0]
                };
                let d = dot64(&h, &n);
                let v = [0, 1, 2].map(|k| h[k] - d * n[k]);
                let l = dot64(&v, &v).sqrt();
                v.map(|x| x / l)
            };
            let e2 = [
                n[1] * e1[2] - n[2] * e1[1],
                n[2] * e1[0] - n[0] * e1[2],
                n[0] * e1[1] - n[1] * e1[0],
            ];
            let (ra, rb) = (0.5 + 1.5 * next(), 0.5 + 1.5 * next());
            let a = [0, 1, 2].map(|k| ra * (ps.cos_rn() * e1[k] + ps.sin_rn() * e2[k]));
            let b = [0, 1, 2].map(|k| rb * (ps.cos_rn() * e2[k] - ps.sin_rn() * e1[k]));
            let k = if i % 4 == 0 { 0.0 } else { 2.0 * next() };
            let dir = [0, 1, 2].map(|_| next() - 0.5);
            let shift = 2f64.powi(-48) * 4.0 * k / dot64(&dir, &dir).sqrt();
            let c2 = [0, 1, 2].map(|j| c[j] + shift * dir[j]);
            let v = |x: [f64; 3]| x.map(q);
            let ca = Crv::Conic {
                c: v(c),
                a: v(a),
                b: v(b),
            };
            let cb = Crv::Conic {
                c: v(c2),
                a: v(a),
                b: v(b),
            };
            let (fa, fb) = (Conic64::of(&ca).unwrap(), Conic64::of(&cb).unwrap());
            let t = TAU * next();
            let pa = fa.param(t);
            if rounding_gap(&fa, pa, &fb, pa) {
                passed += 1;
                let x = curve_at(&ca, &q(pa.0), pa.1).unwrap().0;
                let y = curve_at(&cb, &q(pa.0), pa.1).unwrap().0;
                let w = qsub(&y, &x);
                let nn = qqdot(&w, &w);
                assert!(
                    nn.add_r(&-gap_rounding(&x, &y)).sign() != Ordering::Greater,
                    "{i}"
                );
                tiny += usize::from(k == 0.0);
            }
        }
        assert!(tiny >= 50 && passed > tiny, "{tiny} {passed}");
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

    /// Points of two fields combine exactly (`Mix`): `x = (sqrt 2, 1, 0)`
    /// and `y = (sqrt 3, 1, 1/2)`, their squared gap `5 - 2 sqrt 6 + 1/4`
    /// (`0.3510205...`) told from rationals `1e-12` either side, a direction's
    /// sign, the gap's square; a point of `Q(sqrt 8)` written in `Q(sqrt 2)`
    /// (`rebased`) and one of `Q(sqrt 3)` not.
    #[test]
    fn numbers_of_two_fields_are_exact() {
        let s = |a: i64, b: i64, d: i64| Qd::new(int(a), int(b), int(d));
        let x = [s(0, 1, 2), Qd::rat(int(1)), Qd::rat(zero())];
        let y = [
            s(0, 1, 3),
            Qd::rat(int(1)),
            Qd::rat(R::new(1.into(), 2.into())),
        ];
        let exact = 5.25 - 2.0 * 6f64.sqrt();
        for (k, want) in [(-1e-12, Ordering::Greater), (1e-12, Ordering::Less)] {
            assert_eq!(sq_cmp(&x, &y, &q(exact + k)), want, "{k}");
        }
        assert!(!one_field(&[&x, &y]));
        // (y - x) . (1, 0, 0) = sqrt 3 - sqrt 2 > 0; against (-1, 0, 2):
        // sqrt 2 - sqrt 3 + 1 > 0; against (-3, 0, -2): 3 sqrt 2 - 3 sqrt 3
        // - 1 < 0.
        let w = |v: [i64; 3]| qv(&v.map(int));
        assert_eq!(rel_sign(&y, &x, &w([1, 0, 0])), Ordering::Greater);
        assert_eq!(rel_sign(&y, &x, &w([-1, 0, 2])), Ordering::Greater);
        assert_eq!(rel_sign(&y, &x, &w([-3, 0, -2])), Ordering::Less);
        let n = mv_sub(&mixv(&y, Side::B), &mixv(&x, Side::A));
        let nn = mv_dot(&n, &n);
        let i = nn.interval();
        assert!(rational_f64(&(i.hi() - i.lo())) < 1e-30);
        assert!((rational_f64(&nn.near()) - exact).abs() < 1e-15);
        assert_eq!(
            nn.mul(&nn)
                .sub(&Mix::r(&q(exact * exact * (1.0 + 1e-12))))
                .sign(),
            Ordering::Less
        );
        let eight = [s(1, 1, 8), Qd::rat(int(2)), Qd::rat(zero())];
        let rebased = rebased(&eight, &x).expect("one field");
        assert!(one_field(&[&rebased, &x]));
        assert!(rebased[0].cmp(&s(1, 2, 2)) == Ordering::Equal);
        assert!(super::rebased(&y, &x).is_none());
    }

    /// Every exact point of another curve over a bracket's parameters lies
    /// within `Run::movement`'s bound of the bracket's middle point: a
    /// rod's hole's rim (a meeting of two cylinders), a dimple's rim (a
    /// cylinder's meeting with a sphere, over the angle and the height), a
    /// cone's and a torus's oblique sections (given results' models) and a
    /// spline dome's cap edges, over brackets from `1e-12` to `1e-3` of their
    /// runs in a level frame and a turned one.
    #[test]
    fn movement_bounds_hold_on_other_curves() {
        let tol = Tolerance::default();
        let pi = PI;
        let seam = R::new(BigInt::from(2), BigInt::from(7));
        let turned = Frame3::new(
            Point3::new(1.0, -2.0, 0.5),
            Vec3::new(0.0, 3.0, 4.0),
            Vec3::new(1.0, 0.0, 0.0),
            tol,
        )
        .unwrap();
        let (mut checked, mut kinds) = (0, std::collections::HashSet::new());
        for f in [level(), turned] {
            let at = |l: [f64; 3]| f.point(crate::Point2::new(l[0], l[1]), l[2]);
            let dir = |l: [f64; 3]| {
                let (x, y, n) = (f.x(), f.y(), f.normal());
                Vec3::new(
                    x.x * l[0] + y.x * l[1] + n.x * l[2],
                    x.y * l[0] + y.y * l[1] + n.y * l[2],
                    x.z * l[0] + y.z * l[1] + n.z * l[2],
                )
            };
            let rod = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 3.0, tol)
                .unwrap()
                .0;
            let across =
                Frame3::new(at([-3.0, 0.0, 0.0]), dir([1.0, 0.0, 0.0]), f.y(), tol).unwrap();
            let bar = Solid::cylinder_with(OperationId(2), across, 0.5, 0.0, 6.0, tol)
                .unwrap()
                .0;
            let bite_frame =
                Frame3::new(at([1.0, 0.0, 0.0]), dir([1.0, 0.0, 0.0]), f.y(), tol).unwrap();
            let bite =
                Solid::sphere_with(OperationId(2), bite_frame, 0.5, -pi / 2.0, pi / 2.0, tol)
                    .unwrap()
                    .0;
            let tilt = |o: f64| {
                let g = Frame3::new(at([0.0, 0.0, o]), dir([0.0, 0.6, 0.8]), f.x(), tol).unwrap();
                let square = crate::Boundary::polygon(
                    [(-6.0, -6.0), (6.0, -6.0), (6.0, 6.0), (-6.0, 6.0)]
                        .map(|(x, y)| crate::Point2::new(x, y))
                        .to_vec(),
                    tol,
                )
                .unwrap();
                Solid::extrude_with(
                    OperationId(2),
                    crate::Profile::new(square, vec![], tol).unwrap(),
                    g,
                    0.0,
                    3.0,
                )
                .unwrap()
                .0
            };
            let cone = Solid::cone_with(OperationId(1), f, 1.0, 0.0, 1.0, tol)
                .unwrap()
                .0;
            let torus = Solid::torus_with(OperationId(1), f, 3.0, 1.0, -pi, pi, 2.0 * pi, tol)
                .unwrap()
                .0;
            let bodies = [
                rod.cut(OperationId(6), &bar).unwrap().0.remove(0),
                rod.cut(OperationId(6), &bite).unwrap().0.remove(0),
                cone.cut(OperationId(6), &tilt(0.4)).unwrap().0.remove(0),
                torus.cut(OperationId(6), &tilt(0.5)).unwrap().0.remove(0),
            ];
            let mut models: Vec<Prism> = bodies
                .iter()
                .map(|b| super::super::given::model(b, Operand::A).unwrap())
                .collect();
            let curve = crate::BSplineCurve2::new(
                2,
                vec![
                    crate::Point2::new(4.0, 0.0),
                    crate::Point2::new(2.0, 4.0),
                    crate::Point2::new(0.0, 0.0),
                ],
                None,
                vec![0.0, 1.0],
                vec![3, 3],
            )
            .unwrap();
            let outer = crate::Boundary::path(
                vec![crate::Point2::new(0.0, 0.0), crate::Point2::new(4.0, 0.0)],
                vec![
                    crate::Segment::Line,
                    crate::Segment::Spline(crate::topology::SplineSpan::whole(curve)),
                ],
                tol,
            )
            .unwrap();
            let dome = Solid::extrude_with(
                OperationId(1),
                crate::Profile::new(outer, vec![], tol).unwrap(),
                f,
                0.0,
                2.0,
            )
            .unwrap()
            .0;
            models.push(Prism::new(&dome, Operand::A, &seam).unwrap());
            for m in &models {
                for ei in 0..m.edges.len() {
                    let Some(run) = Run::new(m, ei) else {
                        continue;
                    };
                    for (i, h) in [1e-12, 1e-9, 1e-6, 1e-3].into_iter().enumerate() {
                        let s = 0.13 + 0.21 * i as f64;
                        let ta = run.t0 + run.span * s;
                        let tb = run.t0 + run.span * (s + h);
                        let tm = 0.5 * ta + 0.5 * tb;
                        let Some(x) = run.exact_at(tm) else {
                            continue;
                        };
                        let Some(d) = run.movement(ta, tb, &x) else {
                            continue;
                        };
                        assert!(
                            rational_f64(&d) < 1e-6 + 10.0 * h,
                            "{h} {}",
                            rational_f64(&d)
                        );
                        for j in 0..=8 {
                            let t = ta + (tb - ta) * j as f64 / 8.0;
                            let Some(y) = run.exact_at(t) else {
                                continue;
                            };
                            assert_ne!(sq_cmp(&x, &y, &(&d * &d)), Ordering::Greater, "{h} {j}");
                            checked += 1;
                        }
                        kinds.insert(std::mem::discriminant(run.crv));
                    }
                }
            }
        }
        assert!(
            checked > 200 && kinds.len() >= 4,
            "{checked} {}",
            kinds.len()
        );
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

    /// A level frame's dome prism, a whole torus, a cone, a rod and a ball:
    /// the face pairs' surfaces.
    fn face_models() -> Vec<Prism> {
        let tol = Tolerance::default();
        let seam = R::new(BigInt::from(2), BigInt::from(7));
        let pi = std::f64::consts::PI;
        let turned = Frame3::new(
            Point3::new(0.3, -0.2, 0.1),
            Vec3::new(0.0, 3.0, 4.0),
            Vec3::new(1.0, 0.0, 0.0),
            tol,
        )
        .unwrap();
        let rod = Solid::cylinder_with(OperationId(1), turned, 1.0, -3.0, 3.0, tol)
            .unwrap()
            .0;
        let cone = Solid::cone_with(OperationId(1), level(), 1.0, 0.0, 1.0, tol)
            .unwrap()
            .0;
        let torus = Solid::torus_with(OperationId(1), turned, 3.0, 1.0, -pi, pi, 2.0 * pi, tol)
            .unwrap()
            .0;
        let ball = Solid::sphere_with(OperationId(1), level(), 0.8, -pi / 2.0, pi / 2.0, tol)
            .unwrap()
            .0;
        let curve = crate::BSplineCurve2::new(
            2,
            vec![
                crate::Point2::new(4.0, 0.0),
                crate::Point2::new(2.0, 4.0),
                crate::Point2::new(0.0, 0.0),
            ],
            None,
            vec![0.0, 1.0],
            vec![3, 3],
        )
        .unwrap();
        let outer = crate::Boundary::path(
            vec![crate::Point2::new(0.0, 0.0), crate::Point2::new(4.0, 0.0)],
            vec![
                crate::Segment::Line,
                crate::Segment::Spline(crate::topology::SplineSpan::whole(curve)),
            ],
            tol,
        )
        .unwrap();
        let dome = Solid::extrude_with(
            OperationId(1),
            crate::Profile::new(outer, vec![], tol).unwrap(),
            turned,
            0.0,
            2.0,
        )
        .unwrap()
        .0;
        vec![
            Prism::new(&rod, Operand::A, &seam).unwrap(),
            super::super::cone::model(&cone, Operand::A, &seam).unwrap(),
            super::super::torus::model(&torus, Operand::A, &seam).unwrap(),
            super::super::sphere::model(&ball, Operand::A, &seam).unwrap(),
            Prism::new(&dome, Operand::A, &seam).unwrap(),
        ]
    }

    /// The face pairs' forward jets agree with central differences of their
    /// binary64 values, and their binary64 interval jets hold them, for
    /// every patch (a cylinder, a cone, a torus, a spline wall) against
    /// every field (a plane, a sphere, a cylinder, a cone, a torus, a spline
    /// wall), at several parameters each.
    #[test]
    fn face_jets_hold_their_differences() {
        let models = face_models();
        let faces = |kinds: &dyn Fn(&Surf) -> bool| -> Vec<(usize, usize)> {
            let mut out = Vec::new();
            for (k, m) in models.iter().enumerate() {
                for fi in 0..m.faces.len() {
                    if kinds(&m.faces[fi].surf) && !out.iter().any(|&(j, _)| j == k) {
                        out.push((k, fi));
                    }
                }
            }
            out
        };
        let patches = faces(&|s| {
            matches!(
                s,
                Surf::Cyl { .. } | Surf::Cone { .. } | Surf::Torus | Surf::Spline(_)
            )
        });
        let fields = faces(&|s| !matches!(s, Surf::Plane { .. }));
        let planes = faces(&|s| matches!(s, Surf::Plane { .. }));
        assert_eq!(patches.len(), 4);
        let mut checked = 0;
        for &(pk, pf) in &patches {
            let mp = &models[pk];
            let frame = &mp.view(pf).0.f;
            let p64 = Patch::<f64>::of(mp, pf).unwrap();
            let pj = Patch::<J<f64>>::of(mp, pf).unwrap();
            let pi = Patch::<J<Iv>>::of(mp, pf).unwrap();
            for &(qk, qf) in fields.iter().chain(planes.iter().take(1)) {
                if qk == pk {
                    continue;
                }
                let mq = &models[qk];
                let q64 = Field::<f64>::of(mq, qf, frame).unwrap();
                let qj = Field::<J<f64>>::of(mq, qf, frame).unwrap();
                let qi = Field::<J<Iv>>::of(mq, qf, frame).unwrap();
                let first = if p64.kind == PKind::Spline { 0.2 } else { -2.1 };
                for p in [[first, 0.4], [first + 0.5, 1.3], [first + 0.9, 0.7]] {
                    let phi = |p: [f64; 2]| q64.dist(&p64.local(&p));
                    let (j, _) = jet_at(&pj, &qj, p);
                    let h = 1e-5;
                    let e = |k: usize| {
                        let mut d = [0.0; 2];
                        d[k] = h;
                        d
                    };
                    let at = |d: [f64; 2]| phi([p[0] + d[0], p[1] + d[1]]);
                    let close = |x: f64, y: f64| (x - y).abs() <= 1e-5 * (1.0 + y.abs());
                    assert!(close(j.v, phi(p)));
                    for k in 0..2 {
                        let d1 = (at(e(k)) - at(e(k).map(|x| -x))) / (2.0 * h);
                        assert!(close(j.d[k], d1), "{pk}/{qk} d{k} {} {d1}", j.d[k]);
                    }
                    let d00 = (at(e(0)) - 2.0 * phi(p) + at(e(0).map(|x| -x))) / (h * h);
                    let d11 = (at(e(1)) - 2.0 * phi(p) + at(e(1).map(|x| -x))) / (h * h);
                    assert!((j.h[0] - d00).abs() <= 1e-3 * (1.0 + d00.abs()));
                    assert!((j.h[2] - d11).abs() <= 1e-3 * (1.0 + d11.abs()));
                    if q64.kind != QKind::Spline {
                        let ji = jet_over(&pi, &qi, [Iv::k(p[0]), Iv::k(p[1])]);
                        let holds = |v: &Iv, x: f64| {
                            let (lo, hi) = v.bounds();
                            lo <= x && x <= hi
                        };
                        assert!(holds(&ji.v, j.v) && holds(&ji.d[0], j.d[0]));
                        assert!(holds(&ji.h[2], j.h[2]));
                    }
                    checked += 1;
                }
            }
        }
        assert!(checked >= 54, "{checked}");
    }

    /// Two rods crossing at right angles with their walls `d` apart: within
    /// the resolution `Degenerate`, just beyond it certified apart by
    /// Krawczyk's operator, within rounding of it a `ComputationLimit`.
    #[test]
    fn face_pairs_decide_about_the_resolution() {
        let tol = Tolerance::default();
        let seam = R::new(BigInt::from(2), BigInt::from(7));
        let res = tol.linear();
        let rod = Solid::cylinder_with(OperationId(1), level(), 1.0, -2.0, 2.0, tol)
            .unwrap()
            .0;
        let across = |d: f64| {
            let f = Frame3::new(
                Point3::new(-3.0, 1.5 + d, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                tol,
            )
            .unwrap();
            let bar = Solid::cylinder_with(OperationId(2), f, 0.5, 0.0, 6.0, tol)
                .unwrap()
                .0;
            [
                Prism::new(&rod, Operand::A, &seam).unwrap(),
                Prism::new(&bar, Operand::B, &seam).unwrap(),
            ]
        };
        assert!(matches!(
            face_pairs(&across(0.9 * res), res),
            Err(Error::Degenerate(m)) if m == NEAR_WALLS
        ));
        assert!(face_pairs(&across(1.5 * res), res).is_ok());
        assert!(face_pairs(&across(1e-3), res).is_ok());
        assert!(matches!(
            face_pairs(&across(res * (1.0 + 1e-14)), res),
            Err(Error::ComputationLimit(m)) if m == NEAR_LIMIT
        ));
    }
}
