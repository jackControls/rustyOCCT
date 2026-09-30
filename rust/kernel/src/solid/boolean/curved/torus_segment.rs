//! S9d.4b.1: a torus v-segment or wedge against polyhedral prisms
//! (REVIEW_NOTES.md, S9d.4b.1 evidence). Its exact model is S9d.4a's torus
//! with its ends. A v-segment (latitudes `low < high`, a full turn) is the
//! region between the tube's arc and the axis, revolved: its end discs lie
//! in the planes `w = z` at the stored heights `r sin(latitude)`, each from
//! the axis to the arc's end, where the plane meets the tube on the side of
//! the latitude's cosine (`rho = R +- sqrt(r^2 - z^2)`: a ring over `u` of
//! the plane's spiric section, a circle of radius `R` where the plane is
//! tangent to the torus along it, `z = +-r`). In the meridian half-plane a
//! point lies in it when an odd count of the arc's points at its height lie
//! beyond it: between the critical heights (the ends and `+-r`) the tube's
//! outer and inner points at a height lie on the arc or not all along, so
//! the section there is the tube's disc (both), `t < R + q` (the outer),
//! `t < R - q` (the inner) or nothing, decided exactly by `t^2 - R^2` and
//! the torus's function. A wedge (the whole tube, `0 < angle < 2 pi`) is
//! the torus in the sector from the half-plane of `x` to that of the chart
//! direction `(cos angle, sin angle)` rounded, its end discs the tube's in
//! those half-planes, its rims the half-planes' sections of the torus (rings
//! over `v`).
//!
//! The wall is traced as S9d.4a's in two patches, cut at the ends instead
//! of a seam where an end lies: a segment's at its meridian seam (`e` and
//! `-e`), a wedge's at its parallel seam (`v0` and `v0 + pi`). A plane's
//! section of the wall is S9d.4a's; where the torus's section has a node (a
//! plane tangent to the torus) or comes within the resolution of one off
//! the wall, the wall's part is its two branches over the wall's range (over
//! `v` for a segment, over `u` for a wedge), taken when the graph's
//! discriminant is positive all over a range holding the wall's exactly.
use super::graph::{between_ccw, same_dir};
use super::meet::{conic_point, tangency, CylPair, EdgeMeet, Pos};
use super::model::*;
use super::num::*;
use super::torus::{discriminants, line_quartic, loc_of, verify, Ring, TorusSec};
use super::turned::{near_node_within, roots_repeated, trim, Chart};
use crate::identity::Role;
use crate::profile::boolean::Operand;
use crate::solid::split::{q, rational_f64, zero};
use crate::solid::{Construction, Solid};
use crate::topology::{EdgeId, FaceId, Orientation, Slot};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::f64::consts::TAU;
use std::sync::Arc;

/// The part of the torus a solid is.
#[derive(Debug, Clone)]
pub(super) enum Span {
    /// The whole torus (S9d.4a).
    Whole,
    Band(Box<Band>),
    Wedge(Box<Wedge>),
}

/// A v-segment: whether each end's (low, high) arc point lies on the tube's
/// outer side, `r^2 - z^2` of its height `z`, the ends' directions of `v`
/// (`(rho - R, w)`, unnormalized), the critical heights (sorted) and, per
/// band between them, whether the tube's outer and inner points at its
/// heights lie on the arc; the latitudes and their middle (binary64, for
/// parameters).
#[derive(Debug, Clone)]
pub(super) struct Band {
    outer: [bool; 2],
    q2: [R; 2],
    dirs: [[Qd; 2]; 2],
    crit: Vec<R>,
    cases: Vec<(bool, bool)>,
    lats: [f64; 2],
    pub(super) mid: f64,
}

/// A wedge: its end's chart direction `(cos angle, sin angle)` rounded, that
/// direction as a unit, the turn and its middle.
#[derive(Debug, Clone)]
pub(super) struct Wedge {
    end: P2,
    unit: [Qd; 2],
    angle: f64,
    pub(super) mid: f64,
}

/// A value's sign pushed along directions in turn (at first order: its
/// gradient `g` along each).
fn pushed(v: Qd, g: &QV, ld: &[QV]) -> Ordering {
    let mut s = v.sign();
    for d in ld {
        if s != Ordering::Equal {
            break;
        }
        s = g[0]
            .mul(&d[0])
            .add(&g[1].mul(&d[1]))
            .add(&g[2].mul(&d[2]))
            .sign();
    }
    s
}

/// Both, as `Loc`: out when either is (`Less`), on when either is.
fn both(a: Ordering, b: Ordering) -> Ordering {
    if a == Ordering::Less || b == Ordering::Less {
        Ordering::Less
    } else if a == Ordering::Equal || b == Ordering::Equal {
        Ordering::Equal
    } else {
        Ordering::Greater
    }
}

/// Either: in when either is (`Greater`), out when both are out.
fn either(a: Ordering, b: Ordering) -> Ordering {
    if a == Ordering::Greater || b == Ordering::Greater {
        Ordering::Greater
    } else if a == Ordering::Less && b == Ordering::Less {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

fn loc(s: Ordering) -> Loc {
    loc_of(&[s])
}

fn rat2(d: &P2) -> [Qd; 2] {
    [Qd::rat(d[0].clone()), Qd::rat(d[1].clone())]
}

/// A rational unit direction near a binary64 angle (the half angle's
/// tangent, or its cotangent where that stays smaller).
fn rational_dir(angle: f64) -> P2 {
    let half = angle / 2.0;
    let one = int(1);
    if half.cos().abs() >= 0.5 {
        let s = q(half.tan());
        let den = &one + &s * &s;
        [(&one - &s * &s) / &den, int(2) * &s / &den]
    } else {
        let c = q(half.cos() / half.sin());
        let den = &c * &c + &one;
        [(&c * &c - &one) / &den, int(2) * &c / &den]
    }
}

impl Ring {
    /// Where a point lies in a segment or wedge, pushed along directions in
    /// turn (at first order, as the whole torus's).
    pub(super) fn part_member(&self, f: &Affine, p: &QV, dirs: &[QV]) -> Loc {
        let l = f.local_q(p);
        let ld: Vec<QV> = dirs.iter().map(|d| f.local_dir_q(d)).collect();
        let zero_q = || Qd::rat(zero());
        // Inside the tube, and inside the cylinder `t < R`: `Greater`.
        let tube = pushed(
            self.value(&l).neg(),
            &self.gradient(&l).map(|x| x.neg()),
            &ld,
        );
        match &self.span {
            Span::Whole => unreachable!("a part"),
            Span::Band(band_) => {
                let (crit, cases) = (&band_.crit, &band_.cases);
                let t2 = l[0].mul(&l[0]).add(&l[1].mul(&l[1]));
                let cyl = pushed(
                    t2.neg().add_r(&(&self.big * &self.big)),
                    &[l[0].neg(), l[1].neg(), zero_q()],
                    &ld,
                );
                let band = |j: usize| -> Ordering {
                    match cases[j] {
                        (false, false) => Ordering::Less,
                        (true, true) => tube,
                        (true, false) => either(tube, cyl),
                        (false, true) => both(cyl, tube.reverse()),
                    }
                };
                let up = [zero_q(), zero_q(), Qd::rat(int(1))];
                let hs: Vec<Ordering> = crit
                    .iter()
                    .map(|z| pushed(l[2].add_r(&-z.clone()), &up, &ld))
                    .collect();
                match hs.iter().position(|s| *s == Ordering::Equal) {
                    // On a critical plane after every push: in or out when
                    // the bands on both sides agree.
                    Some(i) => {
                        let (a, b) = (band(i), band(i + 1));
                        if a == b {
                            loc(a)
                        } else {
                            Loc::On
                        }
                    }
                    None => loc(band(hs.iter().filter(|s| **s == Ordering::Greater).count())),
                }
            }
            Span::Wedge(wedge) => {
                let end = &wedge.end;
                let s0 = pushed(l[1].clone(), &[zero_q(), Qd::rat(int(1)), zero_q()], &ld);
                let s1 = pushed(
                    l[0].scale(&end[1]).sub(&l[1].scale(&end[0])),
                    &[Qd::rat(end[1].clone()), Qd::rat(-end[0].clone()), zero_q()],
                    &ld,
                );
                // A turn to a half at most is two half-planes' common, a
                // larger one their union.
                let convex = end[1] > zero() || (end[1] == zero() && end[0] < zero());
                let sector = if convex { both(s0, s1) } else { either(s0, s1) };
                loc(both(tube, sector))
            }
        }
    }

    /// Where a point on a face's surface lies in a segment's or wedge's face
    /// (a wall patch or an end disc).
    pub(super) fn part_in_face(&self, f: &Affine, kind: FaceKind, p: &QV) -> Loc {
        let l = f.local_q(p);
        let side = |s: Ordering, b: bool| {
            if s == Ordering::Equal {
                Ordering::Equal
            } else if (s == Ordering::Greater) == b {
                Ordering::Greater
            } else {
                Ordering::Less
            }
        };
        let within = |d: &[Qd; 2], a: &[Qd; 2], b: &[Qd; 2]| {
            if same_dir(d, a) || same_dir(d, b) {
                Ordering::Equal
            } else if between_ccw(a, d, b) {
                Ordering::Greater
            } else {
                Ordering::Less
            }
        };
        let start = [Qd::rat(int(1)), Qd::rat(zero())];
        match (&self.span, kind) {
            (Span::Band(b), FaceKind::Patch(_, plus)) => {
                let (su, _) = self.sides(&l);
                loc_of(&[
                    side(su, plus),
                    within(&self.v_dir(&l), &b.dirs[0], &b.dirs[1]),
                ])
            }
            (Span::Band(b), FaceKind::Cap(high)) => {
                let (outer, q2) = (&b.outer, &b.q2);
                let k = usize::from(high);
                let tube = self.value(&l).sign();
                let t2 = l[0].mul(&l[0]).add(&l[1].mul(&l[1]));
                let cyl = t2.neg().add_r(&(&self.big * &self.big)).sign();
                let (g, e, lt) = (Ordering::Greater, Ordering::Equal, Ordering::Less);
                // The disc to `R + q` (outer), `R - q` (inner) or `R`.
                let (inside, on) = if q2[k] == zero() {
                    (cyl == g, cyl == e)
                } else if outer[k] {
                    (cyl == g || tube == lt, tube == e)
                } else {
                    (cyl == g && tube == g, cyl == g && tube == e)
                };
                if inside {
                    Loc::In
                } else if on {
                    Loc::On
                } else {
                    Loc::Out
                }
            }
            (Span::Wedge(w), FaceKind::Patch(upper, _)) => {
                let (_, sv) = self.sides(&l);
                loc_of(&[side(sv, upper), within(&Ring::u_dir(&l), &start, &w.unit)])
            }
            (Span::Wedge(w), FaceKind::Cap(high)) => {
                let dir = if high { &w.unit } else { &start };
                let half = mixed_dot_sign(&Ring::u_dir(&l), dir);
                if half != Ordering::Greater {
                    return Loc::Out;
                }
                loc(self.value(&l).sign().reverse())
            }
            _ => unreachable!("a segment's or wedge's patches and discs"),
        }
    }
}

/// Where a segment's or wedge's rim (a ring of a plane's section) meets a
/// plane: the line of the two planes against the torus, the points on the
/// rim's branch (`None` apart, `Along` on the plane).
pub(super) fn rim_plane(c: &TorusSec, p0: &V, m: &V) -> Result<EdgeMeet> {
    let f = &c.f;
    let [a, b, mu, k] = &c.plane;
    // The rim's plane in world terms: nc . (x - o) + k = 0.
    let nc = add(
        &add(&scale(f.row(0), a), &scale(f.row(1), b)),
        &scale(f.row(2), mu),
    );
    let d = cross(&nc, m);
    if is_zero(&d) {
        return Ok(if dot(&nc, &sub(p0, &f.o)) + k == zero() {
            EdgeMeet::Along
        } else {
            EdgeMeet::None
        });
    }
    let (k1, k2) = (dot(&nc, &f.o) - k, dot(m, p0));
    let dd = dot(&d, &d);
    let p = scale(
        &add(&scale(&cross(m, &d), &k1), &scale(&cross(&d, &nc), &k2)),
        &(int(1) / dd),
    );
    // The line meets the plane's two rings (or circles): a repeated root on
    // the other one is no contact with the rim.
    let p = qv(&p);
    let poly = line_quartic(&p, &d, f, &c.ring())?;
    let (sf, rs) = roots_repeated(&poly)?;
    let mut out = Vec::new();
    for (root, repeated) in rs {
        let g = Arc::new(Gen::new(sf.clone(), root));
        let x = qadd(&p, &qscale(&d, &Qd::of(K::generator(&g))));
        if !c.on(&x) {
            continue;
        }
        if repeated {
            return Err(tangency());
        }
        out.push((Pos::Ang(c.place(&x)), x));
    }
    Ok(EdgeMeet::Points(out))
}

/// A segment's or wedge's plane section where the torus's has a node, or
/// comes within the resolution of one: its two branches over a range of the
/// wall's parameter (`v` for a segment, `u` for a wedge) holding the wall's
/// exactly, when the graph's discriminant is positive all over it and no
/// extremum of it there lies within the resolution of zero; `None` else.
pub(super) fn window(
    k: usize,
    f: &Affine,
    ring: &Ring,
    plane: &[R; 4],
    res: f64,
) -> Result<Option<CylPair>> {
    let (du, dv) = discriminants(ring, plane);
    let start = [Qd::rat(int(1)), Qd::rat(zero())];
    let (over_v, form, lo, hi, ends) = match &ring.span {
        Span::Band(b) => (true, dv, b.lats[0], b.lats[1], b.dirs.clone()),
        Span::Wedge(w) => (false, du, 0.0, w.angle, [start, w.unit.clone()]),
        Span::Whole => return Ok(None),
    };
    let delta = ((TAU - (hi - lo)) / 4.0).min(0.05);
    let range = [
        rat2(&rational_dir(lo - delta)),
        rat2(&rational_dir(hi + delta)),
    ];
    if !ends.iter().all(|e| between_ccw(&range[0], e, &range[1])) {
        return Ok(None);
    }
    let mid = rational_dir((lo + hi) / 2.0);
    let chart = Chart {
        c0: mid[0].clone(),
        s0: mid[1].clone(),
    };
    let (Some(t0), Some(t1)) = (chart.t_of(&range[0]), chart.t_of(&range[1])) else {
        return Ok(None);
    };
    let (Some(t0), Some(t1)) = (t0.rational().cloned(), t1.rational().cloned()) else {
        return Ok(None);
    };
    let sign_at = |d: &[Qd; 2]| match [d[0].rational(), d[1].rational()] {
        [Some(c), Some(s)] => sign(&form.value(&[c.clone(), s.clone()])),
        _ => Ordering::Equal,
    };
    let s = sign(&form.value(&mid));
    if s == Ordering::Equal || sign_at(&range[0]) != s || sign_at(&range[1]) != s {
        return Ok(None);
    }
    let p = trim(form.poly(&chart));
    if p.is_empty() || verify(&p, &chart, &range).is_err() {
        return Ok(None);
    }
    if near_node_within(&int(1), &form, &chart, res, &t0, &t1) {
        return Ok(None);
    }
    // Negative all over: the plane misses the wall.
    if s == Ordering::Less {
        return Ok(Some(CylPair::Apart));
    }
    let piece = |plus: bool| {
        Crv::Torus(Box::new(TorusSec {
            carrier: k,
            f: f.clone(),
            big: ring.big.clone(),
            small: ring.small.clone(),
            plane: plane.clone(),
            over_v,
            plus,
            range: Some(range.clone()),
        }))
    };
    Ok(Some(CylPair::Mixed(Box::new(super::spheres::Mixed {
        pieces: vec![piece(true), piece(false)],
        switches: Vec::new(),
    }))))
}

/// The torus's box in its frame (every face's, for filtering).
fn torus_box(f: &Affine, big: &R, small: &R) -> ([f64; 3], [f64; 3]) {
    let reach = big + small;
    let fl = |p: &V| p.clone().map(|x| rational_f64(&x));
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for du in [1, -1] {
        for dv in [1, -1] {
            for dw in [1, -1] {
                let p = fl(&f.point(&(&reach * int(du)), &(&reach * int(dv)), &(small * int(dw))));
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
        }
    }
    for k in 0..3 {
        let m = 1e-9 * (1.0 + lo[k].abs().max(hi[k].abs()));
        lo[k] -= m;
        hi[k] += m;
    }
    (lo, hi)
}

/// A model edge's arc on a closed curve from one place to another,
/// counter-clockwise with its parameter.
#[allow(clippy::too_many_arguments)]
fn arc(
    kind: EdgeKind,
    curve: Crv,
    from: [Qd; 2],
    to: [Qd; 2],
    ends: [usize; 2],
    faces: [usize; 2],
    id: Option<crate::identity::EntityId>,
) -> MEdge {
    MEdge {
        kind,
        curve,
        arc: Some((from, to, true)),
        start: ends[0],
        end: ends[1],
        faces,
        id,
    }
}

/// The exact model of a torus v-segment or wedge (S3's slots: its end discs
/// `FaceId(0)` and `FaceId(1)`, its wall `FaceId(2)`, its rims `EdgeId(0)`
/// and `EdgeId(1)`, no vertices); `seam` picks the rational direction of the
/// meridian seam (a segment's) or the parallel seam (a wedge's).
pub(super) fn model(solid: &Solid, op: Operand, seam: &R) -> Result<Prism> {
    let Construction::Torus {
        major,
        minor,
        low,
        high,
        angle,
        ..
    } = &solid.construction
    else {
        unreachable!("a torus");
    };
    let f = Affine::new(&solid.frame)?;
    let t = &solid.topology;
    let id = |slot: Slot| {
        t.id_of(slot)
            .ok_or(Error::InvalidTopology("an unnamed slot"))
    };
    if t.faces().len() != 3 || t.edges().len() != 2 {
        return Err(Error::InvalidTopology("a torus segment's or wedge's slots"));
    }
    let (big, small) = (q(*major), q(*minor));
    let carrier = usize::from(op == Operand::B);
    let wedge = high - low == TAU;
    let wall = FaceId(2);
    let reversed = t.faces()[wall.0].sense == Orientation::Reversed;
    let rims_stored = [0, 1].map(|k| Some(t.edges()[k].curve.clone()));
    let one = int(1);
    let r2 = &small * &small;
    let sec = |plane: [R; 4], over_v: bool, plus: bool| TorusSec {
        carrier,
        f: f.clone(),
        big: big.clone(),
        small: small.clone(),
        plane,
        over_v,
        plus,
        range: None,
    };
    // Faces: the discs (low or start, high or end), the wall's two patches.
    let disc = |k: usize, m: V, p: V| -> Result<MFace> {
        Ok(MFace {
            kind: FaceKind::Cap(k == 1),
            surf: Surf::Plane { p, m },
            id: id(Slot::Face(FaceId(k)))?,
            stored: t.faces()[k].surface.clone(),
            sense: t.faces()[k].sense,
        })
    };
    let patch = |kind: FaceKind| -> Result<MFace> {
        Ok(MFace {
            kind,
            surf: Surf::Torus,
            id: id(Slot::Face(wall))?,
            stored: t.faces()[wall.0].surface.clone(),
            sense: t.faces()[wall.0].sense,
        })
    };
    let rim_id = |k: usize| id(Slot::Edge(EdgeId(k))).map(Some);
    let mut faces = Vec::new();
    let mut verts: Vec<MVert> = Vec::new();
    let mut edges: Vec<MEdge> = Vec::new();
    let (span, e, v0) = if !wedge {
        let lats = [*low, *high];
        let z = lats.map(|x| q(crate::math::scaled_sin(*minor, x)));
        let outer = lats.map(|x| crate::math::scaled_cos(*minor, x) > 0.0);
        let q2 = [&r2 - &z[0] * &z[0], &r2 - &z[1] * &z[1]];
        if q2.iter().any(|x| *x < zero()) {
            return Err(Error::Degenerate("a torus segment's end off its tube"));
        }
        let sgn = |b: bool| if b { int(1) } else { int(-1) };
        let dirs = [0, 1].map(|k| {
            [
                Qd::new(zero(), sgn(outer[k]), q2[k].clone()),
                Qd::rat(z[k].clone()),
            ]
        });
        let mut crit = vec![-small.clone(), small.clone(), z[0].clone(), z[1].clone()];
        crit.sort();
        crit.dedup();
        let cases = (0..=crit.len())
            .map(|j| {
                if j == 0 || j == crit.len() {
                    return (false, false);
                }
                let h = (&crit[j - 1] + &crit[j]) / int(2);
                let qh = &r2 - &h * &h;
                let at = |s: i64| [Qd::new(zero(), int(s), qh.clone()), Qd::rat(h.clone())];
                (
                    between_ccw(&dirs[0], &at(1), &dirs[1]),
                    between_ccw(&dirs[0], &at(-1), &dirs[1]),
                )
            })
            .collect();
        // A disc's outward normal along the axis: the material toward the
        // other end.
        let up = |k: usize| (k == 1) == (z[1] > z[0]);
        for (k, zk) in z.iter().enumerate() {
            let row = f.row(2).clone();
            faces.push(disc(
                k,
                if up(k) { row } else { neg(&row) },
                f.point(&zero(), &zero(), zk),
            )?);
        }
        faces.push(patch(FaceKind::Patch(true, true))?);
        faces.push(patch(FaceKind::Patch(true, false))?);
        let patch_of = |plus: bool| if plus { 2 } else { 3 };
        let e = circle_point(&[zero(), zero()], &one, seam);
        let (xe, n) = (f.vector(&e[0], &e[1], &zero()), f.n.clone());
        let east = [Qd::rat(e[0].clone()), Qd::rat(e[1].clone())];
        let west = [Qd::rat(-e[0].clone()), Qd::rat(-e[1].clone())];
        // Rims: a plane's ring over `u` (a circle where tangent), each split
        // at the seam's two meridians; vertex `2 k + (west)`.
        for k in 0..2 {
            let curve = if q2[k] == zero() {
                Crv::Conic {
                    c: f.point(&zero(), &zero(), &z[k]),
                    a: scale(&f.x, &big),
                    b: scale(&f.y, &big),
                }
            } else {
                Crv::Torus(Box::new(sec(
                    [zero(), zero(), one.clone(), -z[k].clone()],
                    false,
                    !outer[k],
                )))
            };
            for dir in [&east, &west] {
                let p = match &curve {
                    Crv::Conic { c, a, b } => conic_point(c, a, b, dir),
                    Crv::Torus(s) => s
                        .at(&[
                            dir[0].rational().expect("rational").clone(),
                            dir[1].rational().expect("rational").clone(),
                        ])
                        .ok_or(Error::Degenerate("a torus segment's rim"))?,
                    _ => unreachable!("a rim"),
                };
                verts.push(MVert { p, id: None });
            }
            for (j, (from, to, plus)) in [(&east, &west, true), (&west, &east, false)]
                .into_iter()
                .enumerate()
            {
                let ends = if plus {
                    [2 * k, 2 * k + 1]
                } else {
                    [2 * k + 1, 2 * k]
                };
                // Running with `u`, the disc on the left when its normal is
                // the axis's (seen from outside, as a prism's top cap).
                let faces = if up(k) {
                    [k, patch_of(plus)]
                } else {
                    [patch_of(plus), k]
                };
                edges.push(arc(
                    EdgeKind::Rim(k == 1, j),
                    curve.clone(),
                    from.clone(),
                    to.clone(),
                    ends,
                    faces,
                    rim_id(k)?,
                ));
            }
        }
        // The meridians at `e` and `-e`, running with `v` from the low end
        // to the high: on a forward wall the minus half on the left at `e`,
        // the plus half at `-e` (S9d.4a), the other way inside out.
        let places = [0, 1].map(|k| {
            [
                Qd::new(zero(), sgn(outer[k]) / &small, q2[k].clone()),
                Qd::rat(&z[k] / &small),
            ]
        });
        for (j, (su, left_plus)) in [(1i64, false), (-1, true)].into_iter().enumerate() {
            let left_plus = left_plus != reversed;
            let w = if su == 1 { 0 } else { 1 };
            edges.push(arc(
                EdgeKind::Split(j),
                Crv::Conic {
                    c: add(&f.o, &scale(&xe, &(&big * int(su)))),
                    a: scale(&xe, &(&small * int(su))),
                    b: scale(&n, &small),
                },
                places[0].clone(),
                places[1].clone(),
                [w, 2 + w],
                [patch_of(left_plus), patch_of(!left_plus)],
                None,
            ));
        }
        let span = Span::Band(Box::new(Band {
            outer,
            q2,
            dirs,
            crit,
            cases,
            lats,
            mid: (low + high) / 2.0,
        }));
        (span, e, [one.clone(), zero()])
    } else {
        let end = [q(angle.cos()), q(angle.sin())];
        let n2 = &end[0] * &end[0] + &end[1] * &end[1];
        let unit = [
            Qd::new(zero(), &end[0] / &n2, n2.clone()),
            Qd::new(zero(), &end[1] / &n2, n2.clone()),
        ];
        // The discs: the start's outward normal `-y`'s covector, the end's
        // past its half-plane.
        faces.push(disc(0, neg(f.row(1)), f.o.clone())?);
        faces.push(disc(
            1,
            sub(&scale(f.row(1), &end[0]), &scale(f.row(0), &end[1])),
            f.o.clone(),
        )?);
        faces.push(patch(FaceKind::Patch(true, true))?);
        faces.push(patch(FaceKind::Patch(false, true))?);
        let patch_of = |upper: bool| if upper { 2 } else { 3 };
        let v0 = circle_point(&[zero(), zero()], &one, &(seam + &(int(1) / int(7))));
        let above = [Qd::rat(v0[0].clone()), Qd::rat(v0[1].clone())];
        let below = [Qd::rat(-v0[0].clone()), Qd::rat(-v0[1].clone())];
        // Rims: the half-planes' rings over `v`, each split at the seam's
        // two parallels; vertex `2 k + (below)`.
        for k in 0..2 {
            let plane = if k == 0 {
                [zero(), one.clone(), zero(), zero()]
            } else {
                [-end[1].clone(), end[0].clone(), zero(), zero()]
            };
            let s = sec(plane, true, false);
            for dir in [&v0, &[-v0[0].clone(), -v0[1].clone()]] {
                let p = s.at(dir).ok_or(Error::Degenerate("a torus wedge's rim"))?;
                verts.push(MVert { p, id: None });
            }
            let curve = Crv::Torus(Box::new(s));
            for (j, (from, to, upper)) in [(&above, &below, true), (&below, &above, false)]
                .into_iter()
                .enumerate()
            {
                let ends = if upper {
                    [2 * k, 2 * k + 1]
                } else {
                    [2 * k + 1, 2 * k]
                };
                // Running with `v`: the start's disc on the left, the end's
                // on the right (their outward normals against and along `u`).
                let faces = if k == 0 {
                    [0, patch_of(upper)]
                } else {
                    [patch_of(upper), 1]
                };
                edges.push(arc(
                    EdgeKind::Rim(k == 1, j),
                    curve.clone(),
                    from.clone(),
                    to.clone(),
                    ends,
                    faces,
                    rim_id(k)?,
                ));
            }
        }
        // The parallels at `v0` (upper on the left) and `v0 + pi` (lower on
        // the left), running with `u` from the start to the end.
        for (j, (sv, left_upper)) in [(1i64, true), (-1, false)].into_iter().enumerate() {
            let left_upper = left_upper != reversed;
            let rho = &big + &small * &v0[0] * int(sv);
            let w = if sv == 1 { 0 } else { 1 };
            edges.push(arc(
                EdgeKind::Split(j),
                Crv::Conic {
                    c: add(&f.o, &scale(&f.n, &(&small * &v0[1] * int(sv)))),
                    a: scale(&f.x, &rho),
                    b: scale(&f.y, &rho),
                },
                [Qd::rat(one.clone()), Qd::rat(zero())],
                unit.clone(),
                [w, 2 + w],
                [patch_of(left_upper), patch_of(!left_upper)],
                None,
            ));
        }
        let span = Span::Wedge(Box::new(Wedge {
            end,
            unit,
            angle: *angle,
            mid: angle / 2.0,
        }));
        (span, [one.clone(), zero()], v0)
    };
    let ring = Ring {
        big: big.clone(),
        small: small.clone(),
        e,
        v0,
        span,
        reversed,
        rims: rims_stored,
    };
    let mut info = BTreeMap::new();
    for (eid, _) in t.ids() {
        let role = t.derivation(eid).map_or(Role::External, |d| d.role);
        info.insert(eid, (op, role));
    }
    let bx = torus_box(&f, &big, &small);
    Ok(Prism {
        frame: solid.frame,
        tolerance: solid.resolution(),
        region: id(Slot::Region(crate::topology::RegionId(1)))?,
        f,
        lo: zero(),
        hi: zero(),
        bounds: Vec::new(),
        boxes: vec![bx; faces.len()],
        faces,
        edges,
        verts,
        info,
        ball: None,
        funnel: None,
        ring: Some(ring),
        given: None,
    })
}
