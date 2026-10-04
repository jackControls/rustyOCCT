//! S9d.3a: a cone or frustum in S9c's arrangement (REVIEW_NOTES.md, S9d.3).
//! Its exact model is `u^2 + v^2 <= rho(w)^2`, `0 <= w <= h`, `rho(w) = b +
//! k w`, in the stored frame's coordinates as rationals (affine where its
//! axes are not orthonormal). Its wall is a graph over the plane of `(u,
//! v)`, one face traced in that projection (so no seam splits it: the apex
//! is a point inside it, the vertex loop S3's cones close at); each rim is
//! one closed edge with a vertex at a rational point of it (a seam: a
//! meeting there is tried again at another). A plane meets the wall in a
//! graph over the cone's angle (`ConeSec`): the ruling at `(cos u, sin u)`
//! meets it at `rho(u) = G / (mu + k (alpha cos + beta sin))`, rational at a
//! rational angle, on the solid's nappe where `rho > 0`; a hyperbola's or
//! parabola's piece runs between the directions where the plane is
//! parallel to a ruling (its ends at infinity beyond the end planes).
use super::graph::{between_ccw, rational_between};
use super::meet::{quadratic, tangency, trig, EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use crate::identity::Role;
use crate::profile::boolean::Operand;
use crate::solid::split::{q, rational_f64, zero};
use crate::solid::{Construction, Solid};
use crate::topology::{EdgeId, FaceId, Slot, VertexId};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// A cone model's own data.
#[derive(Debug, Clone)]
pub(super) struct Funnel {
    /// The radii at the bottom (`w = 0`) and the top (`w = h`).
    pub(super) b: R,
    pub(super) t: R,
    /// The radius's rate along the axis, `(t - b) / h` (never zero).
    pub(super) k: R,
    /// The rims' seam: a rational direction of `(u, v)` (unit).
    pub(super) e: P2,
    /// The apex's model vertex, when one radius is zero.
    pub(super) apex: Option<usize>,
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

impl Funnel {
    fn rho(&self, w: &Qd) -> Qd {
        w.scale(&self.k).add_r(&self.b)
    }

    /// Where a point lies in the cone, pushed along directions in turn (at
    /// first order; a push along a ruling or the wall's tangent plane keeps
    /// to the wall).
    pub(super) fn member(&self, f: &Affine, h: &R, p: &QV, dirs: &[QV]) -> Loc {
        let l = f.local_q(p);
        let ld: Vec<QV> = dirs.iter().map(|d| f.local_dir_q(d)).collect();
        let rho = self.rho(&l[2]);
        // rho^2 - u^2 - v^2, positive inside; its first-order change along
        // d is 2 (rho k dw - u du - v dv).
        let mut s = rho
            .mul(&rho)
            .sub(&l[0].mul(&l[0]))
            .sub(&l[1].mul(&l[1]))
            .sign();
        let mut i = 0;
        while s == Ordering::Equal && i < ld.len() {
            s = rho
                .mul(&ld[i][2].scale(&self.k))
                .sub(&l[0].mul(&ld[i][0]))
                .sub(&l[1].mul(&ld[i][1]))
                .sign();
            i += 1;
        }
        let mut sides = vec![s];
        for (bound, above) in [(zero(), true), (h.clone(), false)] {
            let mut s = l[2].add_r(&-bound).sign();
            let mut i = 0;
            while s == Ordering::Equal && i < ld.len() {
                s = ld[i][2].sign();
                i += 1;
            }
            sides.push(if above { s } else { s.reverse() });
        }
        loc_of(&sides)
    }

    /// Where a point on a face's surface lies in the face.
    pub(super) fn in_face(&self, f: &Affine, h: &R, kind: FaceKind, p: &QV) -> Loc {
        let l = f.local_q(p);
        match kind {
            FaceKind::Cap(high) => {
                let r = if high { &self.t } else { &self.b };
                let g = l[0].mul(&l[0]).add(&l[1].mul(&l[1])).add_r(&-(r * r));
                loc_of(&[g.sign().reverse()])
            }
            FaceKind::ConeWall => loc_of(&[l[2].sign(), l[2].add_r(&-h).sign().reverse()]),
            FaceKind::Wall(..) | FaceKind::Half(_) | FaceKind::Patch(..) | FaceKind::Facet(_) => {
                unreachable!("a cone's faces")
            }
        }
    }

    /// Whether a point's `(u, v)` runs along the seam's direction: a
    /// meeting there may be the rims' vertices'.
    pub(super) fn on_seam(&self, f: &Affine, p: &QV) -> bool {
        let l = f.local_q(p);
        let across = l[1].scale(&self.e[0]).sub(&l[0].scale(&self.e[1])).sign();
        let along = l[0].scale(&self.e[0]).add(&l[1].scale(&self.e[1])).sign();
        across == Ordering::Equal && along == Ordering::Greater
    }
}

/// Where a line meets a cone's wall (its quadric's points on either nappe;
/// the face's region decides). A double root is a tangency where the wall
/// is (a line through the apex among them), nothing beyond its ends.
pub(super) fn line_cone(p: &QV, d: &V, f: &Affine, fun: &Funnel, h: &R) -> Result<EdgeMeet> {
    let l = f.local_q(p);
    let ld = f.local_dir(d);
    let rho0 = fun.rho(&l[2]);
    let kd = &fun.k * &ld[2];
    // (u0 + t du)^2 + (v0 + t dv)^2 - (rho0 + t k dw)^2.
    let a = &ld[0] * &ld[0] + &ld[1] * &ld[1] - &kd * &kd;
    let b = l[0]
        .scale(&ld[0])
        .add(&l[1].scale(&ld[1]))
        .sub(&rho0.scale(&kd))
        .scale(&int(2));
    let c = l[0].mul(&l[0]).add(&l[1].mul(&l[1])).sub(&rho0.mul(&rho0));
    let (Some(b), Some(c)) = (b.rational(), c.rational()) else {
        return Err(Error::ComputationLimit("an irrational line against a cone"));
    };
    let point = |t: Qd| {
        let x = qadd(p, &qscale(d, &t));
        (Pos::T(t), x)
    };
    if a == zero() {
        // Along a ruling's direction: one root, or none, or on the wall.
        if *b == zero() {
            return Ok(if *c == zero() {
                EdgeMeet::Along
            } else {
                EdgeMeet::None
            });
        }
        return Ok(EdgeMeet::Points(vec![point(Qd::rat(-c / b))]));
    }
    let disc = b * b - int(4) * &a * c;
    if disc == zero() {
        // A double root: a tangency where the wall is, else none.
        let t = -b / (int(2) * &a);
        let w = l[2].add(&Qd::rat(&t * &ld[2]));
        let beyond = w.sign() == Ordering::Less || w.cmp(&Qd::rat(h.clone())) == Ordering::Greater;
        return if beyond {
            Ok(EdgeMeet::None)
        } else {
            Err(tangency())
        };
    }
    Ok(EdgeMeet::Points(
        quadratic(&a, b, c)?.into_iter().map(point).collect(),
    ))
}

/// A plane's section of a cone's wall over the cone's angle: the plane
/// `alpha u + beta v + mu w + kappa = 0` in the cone's coordinates, its
/// points `o + rho (cos x + sin y) + w n` with `rho = G / D`, `D = mu + k
/// (alpha cos + beta sin)`, `G = b mu - k kappa` (zero through the apex),
/// `w = (rho - b) / k`; `range` the open turn a hyperbola's or parabola's
/// branch runs over (counter-clockwise from its first end), none for a
/// closed curve.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct ConeSec {
    /// The operand whose cone it lies on.
    pub(super) carrier: usize,
    f: Affine,
    b: R,
    k: R,
    pub(super) plane: [R; 4],
    g: R,
    pub(super) range: Option<[[Qd; 2]; 2]>,
}

impl ConeSec {
    fn along(&self, cs: &[Qd; 2]) -> Qd {
        cs[0]
            .scale(&self.plane[0])
            .add(&cs[1].scale(&self.plane[1]))
    }

    fn denom(&self, cs: &[Qd; 2]) -> Qd {
        self.along(cs).scale(&self.k).add_r(&self.plane[2])
    }

    /// The point at a rational `(cos, sin)` (none off the solid's nappe).
    pub(super) fn at(&self, cs: &[R; 2]) -> Option<QV> {
        let d = &self.plane[2] + &self.k * (&cs[0] * &self.plane[0] + &cs[1] * &self.plane[1]);
        if d == zero() {
            return None;
        }
        let rho = &self.g / d;
        if sign(&rho) != Ordering::Greater {
            return None;
        }
        let w = (&rho - &self.b) / &self.k;
        Some(qv(&self.f.point(&(&rho * &cs[0]), &(&rho * &cs[1]), &w)))
    }

    /// Where the section meets a plane `m . (x - p0) = 0` (S9e.3a: a given
    /// result's section as an edge): the line of the two planes against the
    /// cone, its roots on the section kept (on its nappe and branch) and
    /// placed by the cone's angle; the planes parallel, none or `Along`; a
    /// double root a tangency.
    pub(super) fn meet_plane(&self, p0: &V, m: &V) -> Result<EdgeMeet> {
        let f = &self.f;
        let [a, b, mu, kappa] = &self.plane;
        // The section's plane in world terms: nc . (x - o) + kappa = 0.
        let nc = add(
            &add(&scale(f.row(0), a), &scale(f.row(1), b)),
            &scale(f.row(2), mu),
        );
        let d = cross(&nc, m);
        if is_zero(&d) {
            return Ok(if dot(&nc, &sub(p0, &f.o)) + kappa == zero() {
                EdgeMeet::Along
            } else {
                EdgeMeet::None
            });
        }
        let (k1, k2) = (dot(&nc, &f.o) - kappa, dot(m, p0));
        let dd = dot(&d, &d);
        let p = scale(
            &add(&scale(&cross(m, &d), &k1), &scale(&cross(&d, &nc), &k2)),
            &(int(1) / dd),
        );
        // (u0 + t du)^2 + (v0 + t dv)^2 - (b + k (w0 + t dw))^2 = 0.
        let (l, ld) = (f.local(&p), f.local_dir(&d));
        let rho0 = &self.b + &self.k * &l[2];
        let kd = &self.k * &ld[2];
        let qa = &ld[0] * &ld[0] + &ld[1] * &ld[1] - &kd * &kd;
        let qb = int(2) * (&l[0] * &ld[0] + &l[1] * &ld[1] - &rho0 * &kd);
        let qc = &l[0] * &l[0] + &l[1] * &l[1] - &rho0 * &rho0;
        let roots = if qa == zero() {
            if qb == zero() {
                return Ok(if qc == zero() {
                    EdgeMeet::Along
                } else {
                    EdgeMeet::None
                });
            }
            vec![Qd::rat(-&qc / &qb)]
        } else {
            quadratic(&qa, &qb, &qc)?
        };
        let p = qv(&p);
        Ok(EdgeMeet::Points(
            roots
                .into_iter()
                .map(|t| qadd(&p, &qscale(&d, &t)))
                .filter(|x| self.on(x))
                .map(|x| (Pos::Ang(self.place(&x)), x))
                .collect(),
        ))
    }

    /// The section's plane `n . p = d` in world terms and the cone as a
    /// quadric (S9e.3b: the curve's two surfaces).
    pub(super) fn surfaces(&self) -> ((V, R), super::procedural::Other) {
        let f = &self.f;
        let [a, b, mu, kappa] = &self.plane;
        let n = add(
            &add(&scale(f.row(0), a), &scale(f.row(1), b)),
            &scale(f.row(2), mu),
        );
        let d = dot(&n, &f.o) - kappa;
        ((n, d), super::procedural::other_cone(f, &self.b, &self.k))
    }

    /// A point's `(cos, sin)` on the cone.
    pub(super) fn place(&self, x: &QV) -> [Qd; 2] {
        let l = self.f.local_q(x);
        let rho = l[2].scale(&self.k).add_r(&self.b);
        let inv = rho.recip().expect("a section's point off the apex");
        [l[0].mul(&inv), l[1].mul(&inv)]
    }

    /// Whether a point lies on the piece (exactly).
    pub(super) fn on(&self, x: &QV) -> bool {
        let l = self.f.local_q(x);
        let [a, b, m, kappa] = &self.plane;
        let in_plane = l[0]
            .scale(a)
            .add(&l[1].scale(b))
            .add(&l[2].scale(m))
            .add_r(kappa)
            .sign()
            == Ordering::Equal;
        if !in_plane {
            return false;
        }
        let rho = l[2].scale(&self.k).add_r(&self.b);
        if rho.sign() != Ordering::Greater {
            return false;
        }
        let on_cone = l[0]
            .mul(&l[0])
            .add(&l[1].mul(&l[1]))
            .sub(&rho.mul(&rho))
            .sign()
            == Ordering::Equal;
        if !on_cone {
            return false;
        }
        match &self.range {
            None => true,
            Some([lo, hi]) => between_ccw(lo, &self.place(x), hi),
        }
    }

    /// The unit-free tangent at a place, with the angle: `sign(G) (-k A'
    /// (cos x + sin y) - A' n + D (-sin x + cos y))`, `A'` the derivative
    /// of `alpha cos + beta sin` (the derivative times `D^2 / |G|`).
    pub(super) fn tangent(&self, cs: &[Qd; 2]) -> QV {
        let d = self.denom(cs);
        let ap = cs[1]
            .scale(&-self.plane[0].clone())
            .add(&cs[0].scale(&self.plane[1]));
        let radial = qadd(&qscale(&self.f.x, &cs[0]), &qscale(&self.f.y, &cs[1]));
        let turn = qadd(&qscale(&self.f.x, &cs[1].neg()), &qscale(&self.f.y, &cs[0]));
        let t = qadd(
            &qadd(
                &radial.map(|r| r.mul(&ap.scale(&-self.k.clone()))),
                &qscale(&self.f.n, &ap.neg()),
            ),
            &turn.map(|r| r.mul(&d)),
        );
        if sign(&self.g) == Ordering::Less {
            t.map(|x| x.neg())
        } else {
            t
        }
    }

    /// Whether the plane is normal to the axis (its section a circle).
    pub(super) fn normal_to_axis(&self) -> bool {
        self.plane[0] == zero() && self.plane[1] == zero()
    }

    /// A section normal to the axis as the conic it is (S9e.3a: a given
    /// result's edge): `c + a cos + b sin` on the cone's frame, its angle
    /// the cone's (the same places).
    pub(super) fn as_conic(&self) -> Option<Crv> {
        if !self.normal_to_axis() {
            return None;
        }
        let (rho, w) = self.circle();
        Some(Crv::Conic {
            c: self.f.point(&zero(), &zero(), &w),
            a: scale(&self.f.x, &rho),
            b: scale(&self.f.y, &rho),
        })
    }

    /// A plane normal to the axis: its circle's radius and height.
    pub(super) fn circle(&self) -> (R, R) {
        let rho = &self.g / &self.plane[2];
        let w = (&rho - &self.b) / &self.k;
        (rho, w)
    }

    /// Binary64 points from the angle `t0` turning `sweep`.
    pub(super) fn samples(&self, t0: f64, sweep: f64, n: usize) -> Vec<[f64; 3]> {
        let fl = |v: &V| v.clone().map(|y| rational_f64(&y));
        let (o, x, y, nn) = (fl(&self.f.o), fl(&self.f.x), fl(&self.f.y), fl(&self.f.n));
        let [a, b, m, _] = self.plane.clone().map(|y| rational_f64(&y));
        let (g, k, b0) = (
            rational_f64(&self.g),
            rational_f64(&self.k),
            rational_f64(&self.b),
        );
        (0..=n)
            .map(|i| {
                let t = t0 + sweep * i as f64 / n as f64;
                let (s, c) = t.sin_cos();
                let rho = g / (m + k * (a * c + b * s));
                let w = (rho - b0) / k;
                [0, 1, 2].map(|j| o[j] + rho * (c * x[j] + s * y[j]) + w * nn[j])
            })
            .collect()
    }
}

/// A plane's section of a cone's wall (the cone of operand `carrier`): a
/// closed curve, a branch, or nothing on the solid's nappe. A plane
/// through the apex (virtual for a frustum: two rulings) is `Degenerate`.
pub(super) fn plane_cone(
    carrier: usize,
    f: &Affine,
    fun: &Funnel,
    p0: &V,
    m: &V,
) -> Result<Vec<Crv>> {
    let (alpha, beta, mu) = (dot(m, &f.x), dot(m, &f.y), dot(m, &f.n));
    let kappa = dot(m, &sub(&f.o, p0));
    let g = &fun.b * &mu - &fun.k * &kappa;
    if g == zero() {
        return Err(Error::Degenerate("a plane through a cone's apex"));
    }
    let mut sec = ConeSec {
        carrier,
        f: f.clone(),
        b: fun.b.clone(),
        k: fun.k.clone(),
        plane: [alpha.clone(), beta.clone(), mu.clone(), kappa],
        g: g.clone(),
        range: None,
    };
    let want = sign(&g);
    let d_at = |cs: &[R; 2]| sign(&(&mu + &fun.k * (&alpha * &cs[0] + &beta * &cs[1])));
    let rho2 = &alpha * &alpha + &beta * &beta;
    if rho2 == zero() {
        // Normal to the axis: a circle about it.
        return Ok(if sign(&mu) == want {
            vec![Crv::Cone(Box::new(sec))]
        } else {
            Vec::new()
        });
    }
    // The directions of the rulings parallel to the plane: `alpha cos +
    // beta sin = gamma`.
    let gamma = -&mu / &fun.k;
    let disc = &rho2 - &gamma * &gamma;
    match sign(&disc) {
        Ordering::Less => {
            // An ellipse, on one nappe.
            if d_at(&[int(1), zero()]) != want {
                return Ok(Vec::new());
            }
        }
        Ordering::Equal => {
            // A parabola: one direction; its branch on the solid's nappe
            // when the rest is.
            let r = [&alpha * &gamma / &rho2, &beta * &gamma / &rho2];
            if d_at(&[-r[0].clone(), -r[1].clone()]) != want {
                return Ok(Vec::new());
            }
            let r = [Qd::rat(r[0].clone()), Qd::rat(r[1].clone())];
            sec.range = Some([r.clone(), r]);
        }
        Ordering::Greater => {
            // A hyperbola: the turn between its asymptotes' directions on
            // the solid's nappe.
            let Some(roots) = trig(&alpha, &beta, &gamma)? else {
                unreachable!("a plane not normal to the axis")
            };
            let (r0, r1) = (roots[0].clone(), roots[1].clone());
            let mid = rational_between(&r0, &r1, true)?;
            sec.range = Some(if d_at(&mid) == want {
                [r0, r1]
            } else {
                [r1, r0]
            });
        }
    }
    Ok(vec![Crv::Cone(Box::new(sec))])
}

/// The exact model of a cone or frustum; `seam` picks the rational point
/// of the rims their vertices lie at.
pub(super) fn model(solid: &Solid, op: Operand, seam: &R) -> Result<Prism> {
    let Construction::Cone { bottom, top, .. } = &solid.construction else {
        unreachable!("a cone");
    };
    let f = Affine::new(&solid.frame)?;
    let t = &solid.topology;
    let id = |slot: Slot| {
        t.id_of(slot)
            .ok_or(Error::InvalidTopology("an unnamed slot"))
    };
    let h = q(solid.end);
    let radii = [q(*bottom), q(*top)];
    let heights = [zero(), h.clone()];
    let k = (&radii[1] - &radii[0]) / &h;
    let e = circle_point(&[zero(), zero()], &int(1), seam);
    let mut fun = Funnel {
        b: radii[0].clone(),
        t: radii[1].clone(),
        k: k.clone(),
        e: e.clone(),
        apex: None,
    };
    // The ends' discs' planes: `(u, v)` spans them, `x * y` their normal
    // (the stored axes are not exactly orthogonal to `n`).
    let up = cross(&f.x, &f.y);
    let mut faces = Vec::new();
    let mut disc_of = [None, None];
    let mut next = 0;
    for end in 0..2 {
        if radii[end] == zero() {
            continue;
        }
        let fid = FaceId(next);
        next += 1;
        disc_of[end] = Some(faces.len());
        faces.push(MFace {
            kind: FaceKind::Cap(end == 1),
            surf: Surf::Plane {
                p: f.point(&zero(), &zero(), &heights[end]),
                m: if end == 1 { up.clone() } else { neg(&up) },
            },
            id: id(Slot::Face(fid))?,
            stored: t.faces()[fid.0].surface.clone(),
            sense: t.faces()[fid.0].sense,
        });
    }
    let wall_id = FaceId(next);
    let wall = faces.len();
    faces.push(MFace {
        kind: FaceKind::ConeWall,
        surf: Surf::Cone {
            b: radii[0].clone(),
            k,
        },
        id: id(Slot::Face(wall_id))?,
        stored: t.faces()[wall_id.0].surface.clone(),
        sense: t.faces()[wall_id.0].sense,
    });
    let mut verts: Vec<MVert> = Vec::new();
    let mut edges: Vec<MEdge> = Vec::new();
    let place = [Qd::rat(e[0].clone()), Qd::rat(e[1].clone())];
    let mut edge_slot = 0;
    for end in 0..2 {
        let r = &radii[end];
        if *r == zero() {
            fun.apex = Some(verts.len());
            verts.push(MVert {
                p: qv(&f.point(&zero(), &zero(), &heights[end])),
                id: Some(id(Slot::Vertex(VertexId(0)))?),
            });
            continue;
        }
        let v = verts.len();
        verts.push(MVert {
            p: qv(&f.point(&(r * &e[0]), &(r * &e[1]), &heights[end])),
            id: None,
        });
        let disc = disc_of[end].expect("a disc");
        // Counter-clockwise about the axis: a top rim runs forward in its
        // disc (on its left seen from outside), a bottom rim backward.
        edges.push(MEdge {
            kind: EdgeKind::Rim(end == 1, 0),
            curve: Crv::Conic {
                c: f.point(&zero(), &zero(), &heights[end]),
                a: scale(&f.x, r),
                b: scale(&f.y, r),
            },
            arc: Some((place.clone(), place.clone(), true)),
            start: v,
            end: v,
            faces: if end == 1 { [disc, wall] } else { [wall, disc] },
            id: Some(id(Slot::Edge(EdgeId(edge_slot)))?),
        });
        edge_slot += 1;
    }
    let mut info = BTreeMap::new();
    for (eid, _) in t.ids() {
        let role = t.derivation(eid).map_or(Role::External, |d| d.role);
        info.insert(eid, (op, role));
    }
    // Boxes: each end circle's square in the frame (its image holds the
    // circle's), the wall's both.
    let fl = |p: &V| p.clone().map(|x| rational_f64(&x));
    let square = |end: usize| -> Vec<[f64; 3]> {
        let r = &radii[end];
        [(1, 1), (-1, 1), (1, -1), (-1, -1)]
            .iter()
            .map(|&(dx, dy)| fl(&f.point(&(r * int(dx)), &(r * int(dy)), &heights[end])))
            .collect()
    };
    let bound = |pts: Vec<[f64; 3]>| {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &pts {
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
        (lo, hi)
    };
    let boxes = faces
        .iter()
        .map(|face| match face.kind {
            FaceKind::Cap(high) => bound(square(usize::from(high))),
            _ => bound([square(0), square(1)].concat()),
        })
        .collect();
    Ok(Prism {
        frame: solid.frame,
        tolerance: solid.resolution(),
        region: id(Slot::Region(crate::topology::RegionId(1)))?,
        f,
        lo: zero(),
        hi: h,
        bounds: Vec::new(),
        faces,
        edges,
        verts,
        info,
        boxes,
        ball: None,
        funnel: Some(fun),
        ring: None,
        given: None,
        hull: None,
    })
}
