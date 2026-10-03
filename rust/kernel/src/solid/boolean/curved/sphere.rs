//! S9d.1: a sphere (whole, a cap or a zone) in S9c's arrangement
//! (REVIEW_NOTES.md, S9d). Its exact model is `|p - c|^2 <= r^2` of the
//! stored centre and radius, a zone's end planes through `c + h n` normal
//! to the stored axis `n` at the stored heights. Its wall is split into two
//! hemispheres by a plane through the axis at a rational point of the
//! equator (each input its own, tried again at another when a meeting
//! falls on it), each traced in its orthographic projection onto that
//! plane; the poles lie on the split. Circles (plane sections, the rims,
//! the split's great circle) are `Circ`: points `c + dx x + dy y` over a
//! rational orthogonal basis of their plane with `dx^2 |x|^2 + dy^2 |y|^2 =
//! r2`, placed by `(dx, dy)`, so a rational direction gives a point of one
//! quadratic field.
use super::meet::{quadratic, tangency, EdgeMeet, Pos};
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

/// A circle over a rational orthogonal basis of its plane.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Circ {
    pub(super) c: V,
    pub(super) x: V,
    pub(super) y: V,
    pub(super) r2: R,
}

impl Circ {
    /// The normal `x * y`: the circle runs counter-clockwise about it.
    pub(super) fn normal(&self) -> V {
        cross(&self.x, &self.y)
    }

    fn norms(&self) -> (R, R) {
        (dot(&self.x, &self.x), dot(&self.y, &self.y))
    }

    /// A point's place `(dx, dy)`.
    pub(super) fn place(&self, p: &QV) -> [Qd; 2] {
        let d = qsub(p, &qv(&self.c));
        let (xx, yy) = self.norms();
        [
            qdot(&d, &self.x).scale(&(int(1) / xx)),
            qdot(&d, &self.y).scale(&(int(1) / yy)),
        ]
    }

    /// Whether a point lies on the circle (exactly).
    pub(super) fn on(&self, p: &QV) -> bool {
        let d = qsub(p, &qv(&self.c));
        qdot(&d, &self.normal()).sign() == Ordering::Equal
            && qqdot(&d, &d).add_r(&-self.r2.clone()).sign() == Ordering::Equal
    }

    /// The circle's point in a rational direction `(dx, dy)` (scaled onto
    /// it): `c + lambda (dx x + dy y)`, `lambda = sqrt(r2 / (dx^2 |x|^2 +
    /// dy^2 |y|^2))`.
    pub(super) fn at(&self, cs: &[R; 2]) -> QV {
        let (xx, yy) = self.norms();
        let s = &cs[0] * &cs[0] * xx + &cs[1] * &cs[1] * yy;
        let lambda = Qd::new(zero(), int(1), &self.r2 / s);
        let dir = add(&scale(&self.x, &cs[0]), &scale(&self.y, &cs[1]));
        qadd(&qv(&self.c), &qscale(&dir, &lambda))
    }

    /// The unit-free tangent at a point, counter-clockwise.
    pub(super) fn tangent(&self, p: &QV) -> QV {
        let n = qv(&self.normal());
        qcross(&n, &qsub(p, &qv(&self.c)))
    }

    /// A place's binary64 angle in the circle's orthonormalized frame.
    pub(super) fn angle(&self, pos: &[Qd; 2]) -> f64 {
        let (xx, yy) = self.norms();
        let (lx, ly) = (rational_f64(&xx).sqrt(), rational_f64(&yy).sqrt());
        (pos[1].to_f64() * ly).atan2(pos[0].to_f64() * lx)
    }

    /// Binary64 points from the angle `t0` turning `sweep`.
    pub(super) fn samples(&self, t0: f64, sweep: f64, n: usize) -> Vec<[f64; 3]> {
        let f = |v: &V| v.clone().map(|y| rational_f64(&y));
        let (c, x, y) = (f(&self.c), f(&self.x), f(&self.y));
        let (lx, ly) = (
            (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt(),
            (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]).sqrt(),
        );
        let rho = rational_f64(&self.r2).sqrt();
        (0..=n)
            .map(|i| {
                let t = t0 + sweep * i as f64 / n as f64;
                [0, 1, 2].map(|k| c[k] + rho * (t.cos() * x[k] / lx + t.sin() * y[k] / ly))
            })
            .collect()
    }

    /// Where the circle meets the plane through `p0` with normal `m`:
    /// `None` when it lies in it, else its places (a tangency refused).
    pub(super) fn meet_plane(&self, p0: &V, m: &V) -> Result<Option<Vec<[Qd; 2]>>> {
        let (alpha, beta) = (dot(m, &self.x), dot(m, &self.y));
        let gamma = dot(m, &sub(p0, &self.c));
        let (xx, yy) = self.norms();
        if alpha == zero() && beta == zero() {
            return Ok(if gamma == zero() {
                None
            } else {
                Some(Vec::new())
            });
        }
        // alpha dx + beta dy = gamma on dx^2 xx + dy^2 yy = r2.
        let mut out = Vec::new();
        if beta != zero() {
            // dy = (gamma - alpha dx) / beta.
            let qa = &xx + &alpha * &alpha * &yy / (&beta * &beta);
            let qb = int(-2) * &alpha * &gamma * &yy / (&beta * &beta);
            let qc = &gamma * &gamma * &yy / (&beta * &beta) - &self.r2;
            for dx in quadratic(&qa, &qb, &qc)? {
                let dy = dx
                    .scale(&-alpha.clone())
                    .add_r(&gamma)
                    .scale(&(int(1) / &beta));
                out.push([dx, dy]);
            }
        } else {
            let dx = &gamma / &alpha;
            let rest = (&self.r2 - &dx * &dx * &xx) / &yy;
            match sign(&rest) {
                Ordering::Less => {}
                Ordering::Equal => return Err(tangency()),
                Ordering::Greater => {
                    for k in [-1, 1] {
                        out.push([Qd::rat(dx.clone()), Qd::new(zero(), int(k), rest.clone())]);
                    }
                }
            }
        }
        Ok(Some(out))
    }
}

/// A plane's section of a sphere: `None` when they miss (a tangency
/// refused). A plane crossing the sphere within the resolution `res` of
/// tangency is `Degenerate` too, as S9d.2c's circles are: a turned frame's
/// rounding never leaves an exact tangency, and the plane would cut a cap
/// no higher than the resolution, its circle as small as `1e-8`.
pub(super) fn plane_section(c: &V, r: &R, p0: &V, m: &V, res: f64) -> Result<Option<Circ>> {
    let mm = dot(m, m);
    let off = dot(m, &sub(p0, c));
    // The centre's projection and the section's radius squared.
    let centre = add(c, &scale(m, &(&off / &mm)));
    let dd = &off * &off / &mm;
    let r2 = r * r - &dd;
    match sign(&r2) {
        Ordering::Less => return Ok(None),
        Ordering::Equal => return Err(tangency()),
        Ordering::Greater => {}
    }
    // The centre's distance `d < r` at least `r - res`: `(r - res)^2 <=
    // d^2`, always when `r <= res`.
    let lo = r - q(res);
    if sign(&lo) != Ordering::Greater || dd >= &lo * &lo {
        return Err(Error::Degenerate(
            "a plane crossing a sphere within the resolution of tangency (S9d.1)",
        ));
    }
    let x = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        .iter()
        .map(|e| cross(m, &e.map(int)))
        .find(|u| !is_zero(u))
        .expect("a nonzero normal");
    let y = cross(m, &x);
    Ok(Some(Circ {
        c: centre,
        x,
        y,
        r2,
    }))
}

/// Where a line meets a sphere.
pub(super) fn line_sphere(p: &QV, d: &V, c: &V, r: &R) -> Result<EdgeMeet> {
    let a = dot(d, d);
    let w = qsub(p, &qv(c));
    let b = qdot(&w, d).scale(&int(2));
    let cc = qqdot(&w, &w).add_r(&-(r * r));
    let (Some(b), Some(cc)) = (b.rational(), cc.rational()) else {
        return Err(Error::ComputationLimit(
            "an irrational line against a sphere",
        ));
    };
    let roots = quadratic(&a, b, cc)?;
    Ok(EdgeMeet::Points(
        roots
            .into_iter()
            .map(|t| {
                let x = qadd(p, &qscale(d, &t));
                (Pos::T(t), x)
            })
            .collect(),
    ))
}

/// A sphere model's own data.
#[derive(Debug, Clone)]
pub(super) struct Ball {
    pub(super) c: V,
    pub(super) r: R,
    /// The stored axis (not unit).
    pub(super) n: V,
    /// The ends' heights `h` (planes `(p - c) . n = h |n|^2`), `None` a pole.
    pub(super) ends: [Option<R>; 2],
    /// The split plane's normal (through `c` and the axis).
    pub(super) split: V,
    /// The split plane's second axis, orthogonal to `n` (of `n`'s length
    /// where it can be: a whole sphere's, an exact frame's).
    pub(super) w: V,
}

impl Ball {
    fn height(&self, p: &QV) -> Qd {
        qdot(&qsub(p, &qv(&self.c)), &self.n)
    }

    /// The side of an end's plane a point lies on: `Greater` inside.
    fn end_side(&self, k: usize, p: &QV, dirs: &[QV]) -> Option<Ordering> {
        let h = self.ends[k].as_ref()?;
        let nn = dot(&self.n, &self.n);
        let mut s = self.height(p).add_r(&-(h * &nn)).sign();
        let mut i = 0;
        while s == Ordering::Equal && i < dirs.len() {
            s = qdot(&dirs[i], &self.n).sign();
            i += 1;
        }
        Some(if k == 0 { s } else { s.reverse() })
    }

    /// Where a point lies in the ball, pushed along directions in turn (at
    /// first order; a push along the sphere's tangent plane keeps to it).
    pub(super) fn member(&self, p: &QV, dirs: &[QV]) -> Loc {
        let w = qsub(p, &qv(&self.c));
        let mut s = qqdot(&w, &w).add_r(&-(&self.r * &self.r)).sign();
        let mut i = 0;
        while s == Ordering::Equal && i < dirs.len() {
            s = qqdot(&w, &dirs[i]).sign();
            i += 1;
        }
        let mut sides = vec![s.reverse()];
        for k in 0..2 {
            if let Some(e) = self.end_side(k, p, dirs) {
                sides.push(e);
            }
        }
        if sides.contains(&Ordering::Less) {
            Loc::Out
        } else if sides.contains(&Ordering::Equal) {
            Loc::On
        } else {
            Loc::In
        }
    }

    /// Where a point on the sphere lies in a hemisphere face (`side` the
    /// split plane's positive side).
    pub(super) fn in_half(&self, side: bool, p: &QV) -> Loc {
        let s = qdot(&qsub(p, &qv(&self.c)), &self.split).sign();
        let want = if side {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        let mut on = false;
        match s {
            Ordering::Equal => on = true,
            x if x != want => return Loc::Out,
            _ => {}
        }
        for k in 0..2 {
            match self.end_side(k, p, &[]) {
                Some(Ordering::Less) => return Loc::Out,
                Some(Ordering::Equal) => on = true,
                _ => {}
            }
        }
        if on {
            Loc::On
        } else {
            Loc::In
        }
    }

    /// An end's rim circle.
    pub(super) fn rim(&self, k: usize) -> Option<Circ> {
        let h = self.ends[k].as_ref()?;
        let nn = dot(&self.n, &self.n);
        let centre = add(&self.c, &scale(&self.n, h));
        let x0 = cross(&self.split, &self.n);
        let x = if dot(&x0, &x0) == zero() {
            unreachable!("the split holds the axis")
        } else {
            x0
        };
        let y = cross(&self.n, &x);
        Some(Circ {
            c: centre,
            x,
            y,
            r2: &self.r * &self.r - h * h * nn,
        })
    }

    /// The split's great circle, counter-clockwise about the split's normal:
    /// from the top pole through `n x split`'s side.
    pub(super) fn great(&self) -> Circ {
        Circ {
            c: self.c.clone(),
            x: self.n.clone(),
            y: self.w.clone(),
            r2: &self.r * &self.r,
        }
    }

    /// A point's binary64 parameters on a hemisphere: its stereographic
    /// projection from the opposite pole of the split onto the split plane
    /// (`(n, split x n)` coordinates, counter-clockwise about the split's
    /// normal; conformal, so the split's rim is not compressed).
    pub(super) fn params(&self, p: [f64; 3], side: bool) -> [f64; 2] {
        let f = |v: &V| v.clone().map(|y| rational_f64(&y));
        let (c, n, w, s) = (
            f(&self.c),
            f(&self.n),
            f(&cross(&self.split, &self.n)),
            f(&self.split),
        );
        let r = rational_f64(&self.r);
        let d = [(p[0] - c[0]) / r, (p[1] - c[1]) / r, (p[2] - c[2]) / r];
        let dotf = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let up = dotf(&d, &s) / dotf(&s, &s).sqrt();
        let den = if side { 1.0 + up } else { 1.0 - up };
        [
            dotf(&d, &n) / dotf(&n, &n).sqrt() / den,
            dotf(&d, &w) / dotf(&w, &w).sqrt() / den,
        ]
    }
}

/// A rim while the model is built: its end, circle, vertices and places.
type Rim = (usize, Circ, [usize; 2], [[Qd; 2]; 2]);

/// The exact model of a sphere, a cap or a zone; `seam` picks the rational
/// point of the equator the split plane holds.
pub(super) fn model(solid: &Solid, op: Operand, seam: &R) -> Result<Prism> {
    let Construction::Sphere {
        radius, low, high, ..
    } = &solid.construction
    else {
        unreachable!("a sphere");
    };
    let f = Affine::new(&solid.frame)?;
    let t = &solid.topology;
    let id = |slot: Slot| {
        t.id_of(slot)
            .ok_or(Error::InvalidTopology("an unnamed slot"))
    };
    let half = std::f64::consts::FRAC_PI_2;
    let r = q(*radius);
    let ends = [
        (*low != -half).then(|| q(crate::math::scaled_sin(*radius, *low))),
        (*high != half).then(|| q(crate::math::scaled_sin(*radius, *high))),
    ];
    // The split plane: the axis and the equator's rational point. A whole
    // sphere's axis is free: two rows of a rational rotation (world axes,
    // of one rational length, S9d.2), picked by the seam, so its poles
    // leave special points and its great circle's basis is of equal axes.
    let whole = ends.iter().all(|x| x.is_none());
    let (n, u) = if whole {
        const ROTATIONS: [[[i64; 3]; 2]; 4] = [
            [[2, 3, 6], [3, -6, 2]],
            [[1, 4, 8], [4, 7, -4]],
            [[2, 6, 9], [6, 7, -6]],
            [[1, 2, 2], [2, 1, -2]],
        ];
        let k = usize::try_from(seam.numer() % num_bigint::BigInt::from(4))
            .unwrap_or(0)
            .min(3);
        let [a, b] = ROTATIONS[k];
        (a.map(int), b.map(int))
    } else {
        let e = circle_point(&[zero(), zero()], &int(1), seam);
        (f.n.clone(), f.vector(&e[0], &e[1], &zero()))
    };
    let split = cross(&n, &u);
    if is_zero(&split) {
        return Err(Error::ComputationLimit(super::graph::SEAM));
    }
    let w = if whole { u.clone() } else { cross(&split, &n) };
    let ball = Ball {
        c: f.o.clone(),
        r: r.clone(),
        n: n.clone(),
        ends: ends.clone(),
        split,
        w,
    };
    // Faces: the discs of the ends that are not poles (in the input's
    // order), then the two hemispheres (the wall's halves).
    let mut faces = Vec::new();
    let mut disc_of = [None, None];
    let mut next = 0;
    for (k, h) in ends.iter().enumerate() {
        let Some(h) = h else { continue };
        let fid = FaceId(next);
        next += 1;
        let m = if k == 1 { n.clone() } else { neg(&n) };
        disc_of[k] = Some(faces.len());
        faces.push(MFace {
            kind: FaceKind::Cap(k == 1),
            surf: Surf::Plane {
                p: add(&ball.c, &scale(&n, h)),
                m,
            },
            id: id(Slot::Face(fid))?,
            stored: t.faces()[fid.0].surface.clone(),
            sense: t.faces()[fid.0].sense,
        });
    }
    let wall = FaceId(next);
    let halves = [faces.len(), faces.len() + 1];
    for side in [true, false] {
        faces.push(MFace {
            kind: FaceKind::Half(side),
            surf: Surf::Sphere {
                c: ball.c.clone(),
                r: r.clone(),
            },
            id: id(Slot::Face(wall))?,
            stored: t.faces()[wall.0].surface.clone(),
            sense: t.faces()[wall.0].sense,
        });
    }
    // The hemisphere a point off the split lies in.
    let half_of = |p: &QV| -> Result<usize> {
        match qdot(&qsub(p, &qv(&ball.c)), &ball.split).sign() {
            Ordering::Greater => Ok(halves[0]),
            Ordering::Less => Ok(halves[1]),
            Ordering::Equal => Err(Error::ComputationLimit(super::graph::SEAM)),
        }
    };
    let mut verts: Vec<MVert> = Vec::new();
    let mut edges: Vec<MEdge> = Vec::new();
    // The split's great circle's stops: each end's rim (its two points on
    // the split) or its pole.
    let great = ball.great();
    let pole_id = if ends.iter().filter(|e| e.is_none()).count() == 1 {
        Some(id(Slot::Vertex(VertexId(0)))?)
    } else {
        None
    };
    // stops[k] = the great circle's places at end k: [+w side, -w side].
    let mut stops: Vec<[(usize, [Qd; 2]); 2]> = Vec::new();
    let mut rims: Vec<Rim> = Vec::new();
    for (k, end) in ends.iter().enumerate() {
        match end {
            None => {
                // The pole: `c +- r n / |n|`, one vertex on both sides.
                let nn = dot(&n, &n);
                let s = Qd::new(zero(), int(if k == 1 { 1 } else { -1 }), &r * &r / nn);
                let p = qadd(&qv(&ball.c), &qscale(&n, &s));
                let v = verts.len();
                verts.push(MVert {
                    p: p.clone(),
                    id: pole_id,
                });
                let pos = great.place(&p);
                stops.push([(v, pos.clone()), (v, pos)]);
            }
            Some(h) => {
                let rim = ball.rim(k).expect("a rim");
                let pts = rim
                    .meet_plane(&ball.c, &ball.split)?
                    .ok_or(Error::Degenerate("a rim in the split plane"))?;
                if pts.len() != 2 {
                    return Err(Error::Degenerate("a zone's end at a pole"));
                }
                let mut vs = [0usize; 2];
                let places = [pts[0].clone(), pts[1].clone()];
                let mut side: [Option<(usize, [Qd; 2])>; 2] = [None, None];
                let w = cross(&ball.split, &ball.n);
                for (i, cs) in pts.iter().enumerate() {
                    let p = qadd(
                        &qv(&rim.c),
                        &qadd(&qscale(&rim.x, &cs[0]), &qscale(&rim.y, &cs[1])),
                    );
                    let v = verts.len();
                    verts.push(MVert {
                        p: p.clone(),
                        id: None,
                    });
                    vs[i] = v;
                    // Which side of the axis in the split plane: `w`'s sign.
                    let plus = qdot(&qsub(&p, &qv(&ball.c)), &w).sign() == Ordering::Greater;
                    side[if plus { 0 } else { 1 }] = Some((v, great.place(&p)));
                }
                let side = [
                    side[0]
                        .clone()
                        .ok_or(Error::Degenerate("a rim on one side of the axis"))?,
                    side[1]
                        .clone()
                        .ok_or(Error::Degenerate("a rim on one side of the axis"))?,
                ];
                let _ = h;
                stops.push(side);
                rims.push((k, rim, vs, places));
            }
        }
    }
    // Rim arcs: each rim's two halves between its points on the split.
    for (edge_slot, (k, rim, vs, places)) in rims.into_iter().enumerate() {
        let eid = id(Slot::Edge(EdgeId(edge_slot)))?;
        let disc = disc_of[k].expect("a disc");
        for (a, b) in [(0usize, 1usize), (1, 0)] {
            // The arc from a counter-clockwise (about n) to b, and the
            // hemisphere its middle lies in.
            let mid = super::graph::rational_between(&places[a], &places[b], true)?;
            let hemi = half_of(&rim.at(&mid))?;
            // A top rim runs forward in its disc (which lies on its left
            // seen from outside), a bottom rim backward.
            let faces_lr = if k == 1 { [disc, hemi] } else { [hemi, disc] };
            edges.push(MEdge {
                kind: EdgeKind::Rim(k == 1, a),
                curve: Crv::Circle(Box::new(rim.clone())),
                arc: Some((places[a].clone(), places[b].clone(), true)),
                start: vs[a],
                end: vs[b],
                faces: faces_lr,
                id: Some(eid),
            });
        }
    }
    // The split's arcs, counter-clockwise about the split's normal: on the
    // `+w` side from the top stop to the bottom, on the `-w` side back up;
    // the positive hemisphere on their left.
    let (bottom, top) = (&stops[0], &stops[1]);
    for (from, to, j) in [(&top[0], &bottom[0], 0usize), (&bottom[1], &top[1], 1)] {
        edges.push(MEdge {
            kind: EdgeKind::Split(j),
            curve: Crv::Circle(Box::new(great.clone())),
            arc: Some((from.1.clone(), to.1.clone(), true)),
            start: from.0,
            end: to.0,
            faces: [halves[0], halves[1]],
            id: None,
        });
    }
    let mut info = BTreeMap::new();
    for (eid, _) in t.ids() {
        let role = t.derivation(eid).map_or(Role::External, |d| d.role);
        info.insert(eid, (op, role));
    }
    let (lo, hi) = (int(0), int(0));
    let mut out = Prism {
        frame: solid.frame,
        tolerance: solid.resolution(),
        region: id(Slot::Region(crate::topology::RegionId(1)))?,
        f,
        lo,
        hi,
        bounds: Vec::new(),
        faces,
        edges,
        verts,
        info,
        boxes: Vec::new(),
        ball: Some(ball),
        funnel: None,
        ring: None,
        given: None,
    };
    let reach = rational_f64(&r) * (1.0 + 1e-9) + 1e-9;
    let c = out.f.o.clone().map(|x| rational_f64(&x));
    out.boxes = (0..out.faces.len())
        .map(|_| {
            (
                [c[0] - reach, c[1] - reach, c[2] - reach],
                [c[0] + reach, c[1] + reach, c[2] + reach],
            )
        })
        .collect();
    Ok(out)
}
