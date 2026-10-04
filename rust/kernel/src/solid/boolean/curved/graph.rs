//! S9c.1's arrangement: the vertices where an edge of one solid meets a
//! face of the other (and where two equal cylinders' ellipses cross), each
//! edge split at the vertices on it, the section edges (the parts of two
//! faces' meeting curves inside both faces, between vertices), and each
//! face's pieces traced from them: at a vertex the next edge is the first
//! clockwise from the one arrived by (exact angles about the face's
//! normal), each loop's orientation and nesting from its binary64 image in
//! the face's parameters. Each piece is classified at a point of one of
//! its edges pushed into it and then off the face either way (the other
//! solid's exact membership), and kept as S9b.1 keeps fragments.
use super::meet::*;
use super::model::*;
use super::num::*;
use crate::profile::boolean::Op2;
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

/// A seam conflict: a meeting at a full circle's seam, retried at another.
pub(super) const SEAM: &str = "a meeting at a circle's seam";

fn seam() -> Error {
    Error::ComputationLimit(SEAM)
}

/// What a vertex is.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum VKey {
    /// A model vertex of an operand.
    Input(usize, usize),
    /// Edge `e` of operand `o` through face `f` of the other: root `k`.
    Pierce(usize, usize, usize, usize),
    /// Where two equal cylinders' ellipses cross: faces `fa` (A), `fb` (B).
    Cross(usize, usize, usize),
    /// A ring section's own vertex: section and branch (S9d.1).
    Ring(usize, usize),
    /// A stored sphere's pole on a section: section and pole (S9d.1's
    /// follow-up: a section's pcurves turn half a turn there).
    Pole(usize, usize),
}

#[derive(Debug, Clone)]
pub(super) struct Vx {
    pub(super) p: QV,
    pub(super) key: VKey,
    /// The faces it lies on: (operand, face).
    pub(super) faces: BTreeSet<(usize, usize)>,
}

/// The curve an arrangement edge lies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum CurveRef {
    /// A model edge: operand and edge.
    Edge(usize, usize),
    /// A section's branch: section and branch.
    Section(usize, usize),
}

/// An arrangement edge: a curve's piece between two vertices.
#[derive(Debug, Clone)]
pub(super) struct GEdge {
    pub(super) curve: CurveRef,
    pub(super) crv: Crv,
    pub(super) ends: [usize; 2],
    pub(super) pos: [Pos; 2],
    /// Whether the edge runs with the curve's parameter.
    pub(super) with: bool,
    /// A point strictly inside it, and its place.
    pub(super) mid: QV,
    pub(super) mid_pos: Pos,
}

/// Two faces' meeting: A's face `fa`, B's face `fb`, its branches.
#[derive(Debug, Clone)]
pub(super) struct Sec {
    pub(super) fa: usize,
    pub(super) fb: usize,
}

/// A loop of half-edges `(edge, forward)`.
pub(super) type HLoop = Vec<(usize, bool)>;

/// A face's piece: its loops (the outer first) of half-edges `(edge,
/// forward)` running with the face's own normal, whether the other solid
/// holds its front and its back (`sides`), and for the operation (`for_op`)
/// whether the result keeps it and whether its material lies behind the
/// input face (`behind`: the face keeps its orientation).
#[derive(Debug, Clone)]
pub(super) struct Piece {
    pub(super) op: usize,
    pub(super) face: usize,
    pub(super) loops: Vec<Vec<(usize, bool)>>,
    pub(super) sides: (bool, bool),
    pub(super) keep: bool,
    pub(super) behind: bool,
}

#[derive(Debug, Clone)]
pub(super) struct Arr {
    pub(super) models: [Prism; 2],
    pub(super) vx: Vec<Vx>,
    pub(super) edges: Vec<GEdge>,
    pub(super) secs: Vec<Sec>,
    pub(super) pieces: Vec<Piece>,
    /// Faces on one surface: (A's, B's).
    pub(super) coinc: BTreeSet<(usize, usize)>,
}

fn boxes_meet(a: &([f64; 3], [f64; 3]), b: &([f64; 3], [f64; 3])) -> bool {
    (0..3).all(|k| a.0[k] <= b.1[k] && b.0[k] <= a.1[k])
}

/// Whether two faces of one model lie on one quadric (a circle's halves, a
/// sphere's hemispheres).
fn same_quadric(a: &Surf, b: &Surf) -> bool {
    match (a, b) {
        (Surf::Cyl { c, r, .. }, Surf::Cyl { c: c2, r: r2, .. }) => c == c2 && r == r2,
        (Surf::Sphere { c, r }, Surf::Sphere { c: c2, r: r2 }) => c == c2 && r == r2,
        (Surf::Cone { b, k }, Surf::Cone { b: b2, k: k2 }) => b == b2 && k == k2,
        _ => false,
    }
}

/// A model edge's places at its start and end and whether it runs with its
/// curve's parameter: a line's keys (it runs from its start), an arc's
/// angles, a given edge's first arrangement's places (S9e.3a: of every
/// curve).
pub(super) fn edge_places(m: &Prism, ei: usize) -> (Pos, Pos, bool) {
    let e = &m.edges[ei];
    if let Crv::Line { d, .. } = &e.curve {
        return (
            Pos::T(line_key(&m.verts[e.start].p, d)),
            Pos::T(line_key(&m.verts[e.end].p, d)),
            true,
        );
    }
    if let Some(Some(([a, b], with))) = m.given.as_ref().map(|g| &g.places[ei]) {
        return (a.clone(), b.clone(), *with);
    }
    // A spline cap edge over its whole run (S9f.1).
    if let Crv::Spline(c) = &e.curve {
        return (
            Pos::T(Qd::rat(c.seg.first.clone())),
            Pos::T(Qd::rat(c.seg.last.clone())),
            true,
        );
    }
    match &e.arc {
        Some((a, b, ccw)) => (Pos::Ang(a.clone()), Pos::Ang(b.clone()), *ccw),
        None => unreachable!("an arc edge has its ends' angles"),
    }
}

/// Whether a place lies strictly between two others along an edge running
/// from `a` to `b` (with or against its parameter, `ccw`), `None` at an end.
fn strictly_within(pos: &Pos, a: &Pos, b: &Pos, ccw: bool) -> Option<bool> {
    match (pos, a, b) {
        (Pos::T(t), Pos::T(a), Pos::T(b)) => {
            let (lo, hi) = if a.cmp(b) == Ordering::Greater {
                (b, a)
            } else {
                (a, b)
            };
            match (t.cmp(lo), t.cmp(hi)) {
                (Ordering::Greater, Ordering::Less) => Some(true),
                (Ordering::Equal, _) | (_, Ordering::Equal) => None,
                _ => Some(false),
            }
        }
        (Pos::Ang(x2), Pos::Ang(a), Pos::Ang(b)) => {
            if same_dir(x2, a) || same_dir(x2, b) {
                None
            } else {
                Some(between_run(a, x2, b, ccw))
            }
        }
        _ => unreachable!("places of one kind"),
    }
}

/// The line key of a point: `x . d / |d|^2` (its parameter less a constant
/// of the line).
fn line_key(x: &QV, d: &V) -> Qd {
    qdot(x, d).scale(&(int(1) / dot(d, d)))
}

/// The place of a point on a curve.
pub(super) fn place(crv: &Crv, x: &QV) -> Pos {
    match crv {
        Crv::Line { d, .. } => Pos::T(line_key(x, d)),
        Crv::Conic { c, a, b } => Pos::Ang(conic_angle(c, a, b, x)),
        Crv::Circle(c) => Pos::Ang(c.place(x)),
        Crv::Rise(c) => Pos::T(c.height(x)),
        Crv::Meet(m) => Pos::Ang(m.place(x)),
        Crv::Cone(c) => Pos::Ang(c.place(x)),
        Crv::Torus(c) => Pos::Ang(c.place(x)),
        Crv::Toric(c) => Pos::Ang(c.place(x)),
        // S9f.1: a point found on the spline at its arc's parameter.
        Crv::Spline(c) => Pos::T(
            c.place(x)
                .expect("a point on a spline curve at a known parameter"),
        ),
        // S9f.2b: by the segment's run parameter.
        Crv::WallMeet(c) => Pos::T(c.place(x)),
    }
}

/// The cross product's sign of two directions of any fields.
fn cross2(u: &[Qd; 2], v: &[Qd; 2]) -> Ordering {
    mixed_dot_sign(&[u[0].clone(), u[1].neg()], &[v[1].clone(), v[0].clone()])
}

fn dot2(u: &[Qd; 2], v: &[Qd; 2]) -> Ordering {
    mixed_dot_sign(&[u[0].clone(), u[1].clone()], &[v[0].clone(), v[1].clone()])
}

pub(super) fn same_dir(u: &[Qd; 2], v: &[Qd; 2]) -> bool {
    cross2(u, v) == Ordering::Equal && dot2(u, v) == Ordering::Greater
}

/// Whether `x` lies strictly within the counter-clockwise turn from `a` to
/// `b` (a full turn when `a` and `b` agree).
pub(super) fn between_ccw(a: &[Qd; 2], x: &[Qd; 2], b: &[Qd; 2]) -> bool {
    if same_dir(a, x) || same_dir(x, b) {
        return false;
    }
    match cross2(a, b) {
        Ordering::Greater => cross2(a, x) == Ordering::Greater && cross2(x, b) == Ordering::Greater,
        Ordering::Less => !(cross2(b, x) != Ordering::Less && cross2(x, a) != Ordering::Less),
        Ordering::Equal => {
            if dot2(a, b) == Ordering::Greater {
                true
            } else {
                cross2(a, x) == Ordering::Greater
            }
        }
    }
}

/// Whether `x` comes strictly between `a` and `b` running from `a` with
/// the angle (`ccw`) or against it.
fn between_run(a: &[Qd; 2], x: &[Qd; 2], b: &[Qd; 2], ccw: bool) -> bool {
    if ccw {
        between_ccw(a, x, b)
    } else {
        between_ccw(b, x, a)
    }
}

/// The binary64 angle of a `(cos, sin)`.
fn angle_f64(cs: &[Qd; 2]) -> f64 {
    cs[1].to_f64().atan2(cs[0].to_f64())
}

/// A rational `(cos, sin)` strictly between `a` and `b` running with the
/// angle or against it (`b = a`: a full turn).
pub(super) fn rational_between(a: &[Qd; 2], b: &[Qd; 2], ccw: bool) -> Result<[R; 2]> {
    let (t0, t1) = (angle_f64(a), angle_f64(b));
    let mut sweep = if ccw { t1 - t0 } else { t0 - t1 };
    sweep = sweep.rem_euclid(TAU);
    if sweep == 0.0 {
        sweep = TAU;
    }
    for frac in [0.5, 0.25, 0.75, 0.125, 0.875] {
        let t = if ccw {
            t0 + sweep * frac
        } else {
            t0 - sweep * frac
        };
        let half = t / 2.0;
        // The chart whose parameter stays small.
        let cs = if half.cos().abs() >= 0.5 {
            let s = q(half.tan());
            let one = int(1);
            let den = &one + &s * &s;
            [(&one - &s * &s) / &den, int(2) * &s / &den]
        } else {
            let c = q(half.cos() / half.sin());
            let one = int(1);
            let den = &c * &c + &one;
            [(&c * &c - &one) / &den, int(2) * &c / &den]
        };
        let x = [Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())];
        if between_run(a, &x, b, ccw) {
            return Ok(cs);
        }
    }
    Err(Error::Degenerate(
        "two meetings within rounding along an arc",
    ))
}

/// A rational strictly between two numbers (`a < b`).
pub(super) fn rational_between_num(a: &Qd, b: &Qd) -> Result<R> {
    let (ia, ib) = (a.interval(), b.interval());
    let m = (ia.hi() + ib.lo()) / int(2);
    let x = Qd::rat(m.clone());
    if a.cmp(&x) == Ordering::Less && x.cmp(b) == Ordering::Less {
        return Ok(m);
    }
    Err(Error::Degenerate(
        "two meetings within rounding along a line",
    ))
}

/// Builds the arrangement of two models, the one every operation shares:
/// its pieces' sides decided, not yet kept or dropped (`Arr::for_op`).
pub(super) fn arrange_shared(models: [Prism; 2]) -> Result<Arr> {
    // Cylinder pairs, and faces on one surface (A's, B's).
    let mut pairs: BTreeMap<(usize, usize), CylPair> = BTreeMap::new();
    let mut coinc: BTreeSet<(usize, usize)> = BTreeSet::new();
    // A given result's faces (S9e.1) reach their surfaces' frames and data
    // through their inputs' models (`Prism::view`).
    for fa in 0..models[0].faces.len() {
        let (va, ia) = models[0].view(fa);
        let Surf::Cyl { c: ca, r: ra, .. } = &va.faces[ia].surf else {
            continue;
        };
        for fb in 0..models[1].faces.len() {
            let (vb, ib) = models[1].view(fb);
            let Surf::Cyl { c: cb, r: rb, .. } = &vb.faces[ib].surf else {
                continue;
            };
            let bx = [&models[0].boxes[fa], &models[1].boxes[fb]];
            let pair = cyl_pair(va, ca, ra, bx, vb, cb, rb)?;
            if matches!(pair, CylPair::Same) && boxes_meet(bx[0], bx[1]) {
                coinc.insert((fa, fb));
            }
            pairs.insert((fa, fb), pair);
        }
    }
    // A spline wall and a cylinder on crossing axes (S9f.2b) or a sphere
    // (S9f.3a): the meeting's graphs over the run and over the height and
    // their switches (S9f.2b.2).
    for fa in 0..models[0].faces.len() {
        let (va, ia) = models[0].view(fa);
        for fb in 0..models[1].faces.len() {
            let (vb, ib) = models[1].view(fb);
            let crossing =
                |s: &Prism, c: &Prism| super::spline_parallel::map2(&s.f, &c.f).is_none();
            let meets = boxes_meet(&models[0].boxes[fa], &models[1].boxes[fb]);
            let pair = match (&va.faces[ia].surf, &vb.faces[ib].surf) {
                (Surf::Spline(s), Surf::Cyl { c, r, .. }) if crossing(va, vb) => {
                    if meets {
                        super::spline_crossing::meeting(va, ia, s, vb, ib, c, r)?
                    } else {
                        CylPair::Apart
                    }
                }
                (Surf::Cyl { c, r, .. }, Surf::Spline(s)) if crossing(vb, va) => {
                    if meets {
                        super::spline_crossing::meeting(vb, ib, s, va, ia, c, r)?
                    } else {
                        CylPair::Apart
                    }
                }
                // S9f.3a: a spline wall and a sphere, the same way.
                (Surf::Spline(s), Surf::Sphere { c, r }) => {
                    if meets {
                        super::spline_sphere::meeting(va, ia, s, vb, ib, c, r)?
                    } else {
                        CylPair::Apart
                    }
                }
                (Surf::Sphere { c, r }, Surf::Spline(s)) => {
                    if meets {
                        super::spline_sphere::meeting(vb, ib, s, va, ia, c, r)?
                    } else {
                        CylPair::Apart
                    }
                }
                _ => continue,
            };
            pairs.insert((fa, fb), pair);
        }
    }
    // A cylinder and a sphere (S9d.2): rings over the cylinder's angle, or
    // apart.
    let res = models[0].tolerance.linear();
    for fa in 0..models[0].faces.len() {
        let (va, ia) = models[0].view(fa);
        for fb in 0..models[1].faces.len() {
            let (vb, ib) = models[1].view(fb);
            let pair = match (&va.faces[ia].surf, &vb.faces[ib].surf) {
                (Surf::Cyl { c: ca, r: ra, .. }, Surf::Sphere { c, r }) => {
                    if !boxes_meet(&models[0].boxes[fa], &models[1].boxes[fb]) {
                        CylPair::Apart
                    } else {
                        super::spheres::sphere_cyl(0, (&va.f, ca, ra), c, r, res)?
                    }
                }
                (Surf::Sphere { c, r }, Surf::Cyl { c: cb, r: rb, .. }) => {
                    if !boxes_meet(&models[0].boxes[fa], &models[1].boxes[fb]) {
                        CylPair::Apart
                    } else {
                        super::spheres::sphere_cyl(1, (&vb.f, cb, rb), c, r, res)?
                    }
                }
                _ => continue,
            };
            pairs.insert((fa, fb), pair);
        }
    }
    // A cone and a curved face (S9d.3b): rings over a carrier, a plane, or
    // apart.
    for fa in 0..models[0].faces.len() {
        for fb in 0..models[1].faces.len() {
            let (sa, sb) = (&models[0].faces[fa].surf, &models[1].faces[fb].surf);
            let cone = matches!(sa, Surf::Cone { .. }) || matches!(sb, Surf::Cone { .. });
            let curved = |s: &Surf| !matches!(s, Surf::Plane { .. } | Surf::Torus);
            if !cone || !curved(sa) || !curved(sb) {
                continue;
            }
            let pair = if boxes_meet(&models[0].boxes[fa], &models[1].boxes[fb]) {
                let ((va, ia), (vb, ib)) = (models[0].view(fa), models[1].view(fb));
                super::cones::cone_pair([va, vb], [ia, ib], res)?
            } else {
                CylPair::Apart
            };
            pairs.insert((fa, fb), pair);
        }
    }
    // A plane and a torus (S9d.4a): its spiric section, rings or loops,
    // the same for every patch (found once for each plane).
    let mut spiric: BTreeMap<(usize, usize), CylPair> = BTreeMap::new();
    for fa in 0..models[0].faces.len() {
        for fb in 0..models[1].faces.len() {
            let (sa, sb) = (&models[0].faces[fa].surf, &models[1].faces[fb].surf);
            let (k, plane) = match (sa, sb) {
                (Surf::Torus, Surf::Plane { p, m }) => (0, (p, m)),
                (Surf::Plane { p, m }, Surf::Torus) => (1, (p, m)),
                _ => continue,
            };
            let pair = if boxes_meet(&models[0].boxes[fa], &models[1].boxes[fb]) {
                let key = (k, if k == 0 { fb } else { fa });
                if let Some(pair) = spiric.get(&key) {
                    pair.clone()
                } else {
                    let t = models[k].view(if k == 0 { fa } else { fb }).0;
                    let pair = super::torus::plane_torus(
                        k,
                        &t.f,
                        t.ring.as_ref().expect("a torus"),
                        plane.0,
                        plane.1,
                        res,
                    )?;
                    spiric.insert(key, pair.clone());
                    pair
                }
            } else {
                CylPair::Apart
            };
            pairs.insert((fa, fb), pair);
        }
    }
    // A torus and a quadric face (S9d.4b.2): its meeting's pieces, the same
    // for every patch and for every face on one quadric (found once); two
    // tori (S9d.4b.2b) over the first's angles, the same for every pair of
    // patches.
    let mut toric: Vec<((usize, usize), CylPair)> = Vec::new();
    for fa in 0..models[0].faces.len() {
        for fb in 0..models[1].faces.len() {
            let (sa, sb) = (&models[0].faces[fa].surf, &models[1].faces[fb].surf);
            let quadric = |s: &Surf| {
                matches!(
                    s,
                    Surf::Cyl { .. } | Surf::Sphere { .. } | Surf::Cone { .. }
                )
            };
            let (k, g) = match (sa, sb) {
                (Surf::Torus, Surf::Torus) => (0, fb),
                (Surf::Torus, s) if quadric(s) => (0, fb),
                (s, Surf::Torus) if quadric(s) => (1, fa),
                _ => continue,
            };
            let pair = if boxes_meet(&models[0].boxes[fa], &models[1].boxes[fb]) {
                let other = &models[1 - k];
                let tori = matches!(other.faces[g].surf, Surf::Torus);
                // The first face on the same quadric (or torus).
                let first = if tori {
                    0
                } else {
                    (0..=g)
                        .find(|&h| {
                            std::ptr::eq(other.view(h).0, other.view(g).0)
                                && same_quadric(&other.faces[h].surf, &other.faces[g].surf)
                        })
                        .unwrap_or(g)
                };
                if let Some((_, pair)) = toric.iter().find(|(key, _)| *key == (k, first)) {
                    pair.clone()
                } else {
                    let own = models[k].view(if k == 0 { fa } else { fb }).0;
                    let (vo, io) = other.view(g);
                    let pair = if tori {
                        super::torus_curved::torus_torus(k, own, vo)?
                    } else {
                        super::torus_curved::torus_quadric(
                            k,
                            own,
                            &super::cones::other_face(vo, io),
                        )?
                    };
                    toric.push(((k, first), pair.clone()));
                    pair
                }
            } else {
                CylPair::Apart
            };
            pairs.insert((fa, fb), pair);
        }
    }
    // Coincident planes; planes within the resolution of each other over
    // their faces' boxes' overlap are one plane within it (a sliver between
    // them), parallel or not: a frame normalized again turns its normal by
    // the platform's `hypot`, so a plane rounded parallel to another on one
    // host is off parallel by an ulp on another.
    let tol = q(res);
    for (fa, a) in models[0].faces.iter().enumerate() {
        let Surf::Plane { p: pa, m: ma } = &a.surf else {
            continue;
        };
        for (fb, b) in models[1].faces.iter().enumerate() {
            let Surf::Plane { p: pb, m: mb } = &b.surf else {
                continue;
            };
            if !boxes_meet(&models[0].boxes[fa], &models[1].boxes[fb]) {
                continue;
            }
            let turn = cross(ma, mb);
            let parallel = is_zero(&turn);
            let (lo, hi) = intersect(&models[0].boxes[fa], &models[1].boxes[fb]);
            let (lo, hi) = (lo.map(q), hi.map(q));
            let span = sub(&hi, &lo);
            let (aa, bb) = (dot(ma, ma), dot(mb, mb));
            // Turned apart by more than the resolution across the overlap.
            if !parallel && dot(&turn, &turn) * dot(&span, &span) > &tol * &tol * &aa * &bb {
                continue;
            }
            if parallel && dot(ma, &sub(pb, pa)) == zero() {
                coinc.insert((fa, fb));
                continue;
            }
            // Plane a's offset from the point of plane b nearest the overlap's
            // centre (plane b's own offset when they are parallel).
            let centre = scale(&add(&lo, &hi), &(int(1) / int(2)));
            let foot = sub(&centre, &scale(mb, &(dot(mb, &sub(&centre, pb)) / &bb)));
            let off = dot(ma, &sub(&foot, pa));
            if &off * &off <= &tol * &tol * &aa {
                return Err(Error::Degenerate(
                    "two faces within the resolution of one plane",
                ));
            }
        }
    }
    let in_coinc: [BTreeSet<usize>; 2] = [
        coinc.iter().map(|x| x.0).collect(),
        coinc.iter().map(|x| x.1).collect(),
    ];
    let coincident_with = |o: usize, own: usize, other: usize| -> bool {
        if o == 0 {
            coinc.contains(&(own, other))
        } else {
            coinc.contains(&(other, own))
        }
    };
    let pair_of = |o: usize, own: usize, other: usize| -> Option<&CylPair> {
        if o == 0 {
            pairs.get(&(own, other))
        } else {
            pairs.get(&(other, own))
        }
    };
    // S9f.1: a spline prism's vertical edge between two spline walls (a
    // joint, smooth or not) lying in a plane of the other input, across
    // that plane's face: the plane crosses the solid's boundary along it,
    // so the face holds the edge's parts inside it and meets the walls
    // there, not in generatrices of its own: (operand, edge, other's face).
    let (along, touching) = joint_edges_on_planes(&models)?;
    // Vertices: the inputs'.
    let mut vx: Vec<Vx> = Vec::new();
    // Vertices' binary64 views where computed (`qv_f64`; a vertex's point
    // never changes).
    let mut views: Vec<Option<[f64; 3]>> = Vec::new();
    let mut on_edge: BTreeMap<(usize, usize), Vec<(usize, Pos)>> = BTreeMap::new();
    let mut input_vx: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
    for (o, m) in models.iter().enumerate() {
        for (i, v) in m.verts.iter().enumerate() {
            input_vx[o].push(vx.len());
            vx.push(Vx {
                p: v.p.clone(),
                key: VKey::Input(o, i),
                faces: BTreeSet::new(),
            });
        }
        for (ei, e) in m.edges.iter().enumerate() {
            let (s, t) = (input_vx[o][e.start], input_vx[o][e.end]);
            for f in e.faces {
                vx[s].faces.insert((o, f));
                vx[t].faces.insert((o, f));
            }
            let (ps, pt, _) = edge_places(m, ei);
            on_edge
                .entry((o, ei))
                .or_default()
                .extend([(s, ps), (t, pt)]);
        }
    }
    // Pierces: every edge against every face of the other.
    for o in 0..2 {
        let (me, other) = (&models[o], &models[1 - o]);
        for (ei, e) in me.edges.iter().enumerate() {
            let ebox = intersect(&me.boxes[e.faces[0]], &me.boxes[e.faces[1]]);
            // An arc edge's own cylinder (its wall; a given result's conic's
            // carrier, S9e.1), on its input's model.
            let wall = match &me.given {
                Some(g) => g.own[ei],
                None => e
                    .faces
                    .iter()
                    .copied()
                    .find(|&f| matches!(me.faces[f].surf, Surf::Cyl { .. })),
            };
            let own = match (&e.curve, wall.map(|w| me.view(w))) {
                (Crv::Conic { .. }, Some((vm, vw))) => match &vm.faces[vw].surf {
                    Surf::Cyl { c, r, .. } => Some((vm, c, r)),
                    Surf::Plane { .. }
                    | Surf::Sphere { .. }
                    | Surf::Cone { .. }
                    | Surf::Torus
                    | Surf::Spline(_) => None,
                },
                _ => None,
            };
            let virtual_edge = e.id.is_none();
            // A given edge of a curve whose meetings are S9e.3b's (three
            // surfaces, `triple.rs`): none where one of its faces' surfaces
            // is apart from the other's face (no exact work there).
            let unmet = me.given.is_some()
                && match &e.curve {
                    Crv::Meet(_) | Crv::Rise(_) | Crv::Toric(_) | Crv::Cone(_) => true,
                    Crv::Torus(c) => !c.fixed_angle(),
                    _ => false,
                };
            for (g, gface) in other.faces.iter().enumerate() {
                if !boxes_meet(&ebox, &other.boxes[g]) {
                    continue;
                }
                if unmet
                    && !(matches!(e.curve, Crv::Cone(_) | Crv::Torus(_))
                        && matches!(
                            other.view(g).0.faces[other.view(g).1].surf,
                            Surf::Plane { .. }
                        ))
                    && e.faces.iter().any(|&f| {
                        matches!(pair_of(o, f, g), Some(CylPair::Apart))
                            || super::chain::planes_apart(me, f, other, g)
                    })
                {
                    continue;
                }
                let pair = match (&gface.surf, wall, &e.curve) {
                    (Surf::Cyl { .. }, Some(w), Crv::Conic { .. }) => pair_of(o, w, g),
                    _ => None,
                };
                let (vo, vg) = other.view(g);
                // Parallel cylinders' circles are in the first input's frame:
                // an arc of the second input's cylinder meets the first's
                // where its own circle meets the first's circle in its own
                // frame (found by S9e.1's chains: the relation read as the
                // first input's placed the second's arc at wrong angles).
                let reversed;
                let pair = match (pair, own, &vo.faces[vg].surf) {
                    (
                        Some(CylPair::Parallel { .. }),
                        Some((vm, cx, rx)),
                        Surf::Cyl { c, r, .. },
                    ) if o == 1 => {
                        let wall = wall.expect("an arc's own cylinder");
                        let bx = [&me.boxes[wall], &other.boxes[g]];
                        reversed = cyl_pair(vm, cx, rx, bx, vo, c, r)?;
                        Some(&reversed)
                    }
                    (p, ..) => p,
                };
                // A sphere's circle's own sphere: a given edge's through its
                // faces' views (S9e.3a).
                let own_ball = match &me.given {
                    Some(_) => super::chain::own_model(me, ei, |v| v.ball.is_some())
                        .and_then(|v| v.ball.as_ref()),
                    None => me.ball.as_ref(),
                };
                // The edge's and the face's boxes (S9e.3b: a three
                // surfaces' point certainly outside them not constructed).
                let clip = intersect(&ebox, &other.boxes[g]);
                let meet = match edge_surface(&e.curve, own, own_ball, vo, vg, pair, Some(&clip)) {
                    // A seam's tangency with a torus, or a torus seam's
                    // (S9d.4b.2): another seam is tried.
                    Err(Error::Degenerate(_))
                        if virtual_edge && (me.ring.is_some() || other.ring.is_some()) =>
                    {
                        return Err(seam())
                    }
                    r => r?,
                };
                let points = match meet {
                    EdgeMeet::None => continue,
                    EdgeMeet::Along => {
                        // The edge on the face's surface: taken with the
                        // faces on one surface when one of its faces is,
                        // or a spline joint's edge across a plane's face.
                        if e.faces.iter().any(|&f| coincident_with(o, f, g))
                            || along.contains(&(o, ei, g))
                        {
                            continue;
                        }
                        return Err(if virtual_edge {
                            seam()
                        } else {
                            Error::Degenerate("an edge of one input on a face of the other")
                        });
                    }
                    EdgeMeet::Points(p) => p,
                };
                for (k, (_, x)) in points.into_iter().enumerate() {
                    let pos = place(&e.curve, &x);
                    let (ps, pt) = {
                        let l = &on_edge[&(o, ei)];
                        (l[0].1.clone(), l[1].1.clone())
                    };
                    let inside = strictly_within(&pos, &ps, &pt, edge_places(me, ei).2);
                    let region = other.in_face(g, &x);
                    let end_vertex = if same_end(&pos, &ps) {
                        input_vx[o][e.start]
                    } else {
                        input_vx[o][e.end]
                    };
                    match (inside, region) {
                        (_, Loc::Out) | (Some(false), _) => continue,
                        (None, Loc::In)
                            if vx[end_vertex]
                                .faces
                                .iter()
                                .any(|&(o2, f)| o2 == o && coincident_with(o, f, g)) =>
                        {
                            // A vertex inside a face on its own face's
                            // surface: a vertex of that face's pieces too.
                            vx[end_vertex].faces.insert((1 - o, g));
                            continue;
                        }
                        (Some(true), Loc::On) => {
                            // On an edge of the face whose other face lies on
                            // one of this edge's faces' surface: two edges
                            // crossing on that surface.
                            let seamy = virtual_edge || seam_at(other, g, &x);
                            let meeting = || {
                                if seamy {
                                    seam()
                                } else {
                                    Error::Degenerate(
                                        "an edge of one input meeting an edge of the other",
                                    )
                                }
                            };
                            let Some((f, fpos)) = edge_at(other, g, &x).map_err(|_| meeting())?
                            else {
                                return Err(meeting());
                            };
                            let fe = &other.edges[f];
                            // Or the two edges in one plane, one of them a
                            // spline joint's edge across it (S9f.1).
                            let on_surface =
                                fe.faces.iter().any(|&h| {
                                    h != g && e.faces.iter().any(|&ef| coincident_with(o, ef, h))
                                }) || e.faces.iter().any(|&ef| along.contains(&(1 - o, f, ef)))
                                    || fe
                                        .faces
                                        .iter()
                                        .any(|&h| h != g && along.contains(&(o, ei, h)));
                            if !on_surface {
                                return Err(meeting());
                            }
                            let (id, pos, fpos) = match vx.iter().position(|v| qv_eq(&v.p, &x)) {
                                // Places from the vertex's own point: its
                                // field, where another root found it
                                // (S9f.1).
                                Some(id) => (
                                    id,
                                    place(&e.curve, &vx[id].p),
                                    place(&other.edges[f].curve, &vx[id].p),
                                ),
                                None => {
                                    vx.push(Vx {
                                        p: x.clone(),
                                        key: VKey::Pierce(o, ei, g, k),
                                        faces: BTreeSet::new(),
                                    });
                                    (vx.len() - 1, pos, fpos)
                                }
                            };
                            for &ef in &e.faces {
                                vx[id].faces.insert((o, ef));
                            }
                            for &h in &fe.faces {
                                vx[id].faces.insert((1 - o, h));
                            }
                            let list = on_edge.entry((o, ei)).or_default();
                            if !list.iter().any(|(v, _)| *v == id) {
                                list.push((id, pos));
                            }
                            let list = on_edge.entry((1 - o, f)).or_default();
                            if !list.iter().any(|(v, _)| *v == id) {
                                list.push((id, fpos));
                            }
                            continue;
                        }
                        (None, _) => {
                            // An input vertex on the other's face.
                            let seamy = virtual_edge
                                || seam_at(other, g, &x)
                                || me.verts[if same_end(&pos, &ps) { e.start } else { e.end }]
                                    .id
                                    .is_none();
                            return Err(if seamy {
                                seam()
                            } else {
                                Error::Degenerate("a vertex of one input on the other's face")
                            });
                        }
                        (Some(true), Loc::In) => {}
                    }
                    let id = vx.len();
                    let mut faces: BTreeSet<(usize, usize)> =
                        e.faces.iter().map(|&f| (o, f)).collect();
                    faces.insert((1 - o, g));
                    vx.push(Vx {
                        p: x,
                        key: VKey::Pierce(o, ei, g, k),
                        faces,
                    });
                    on_edge.entry((o, ei)).or_default().push((id, pos));
                }
            }
        }
    }
    // The crossings of equal cylinders' ellipses.
    for ((fa, fb), pair) in &pairs {
        // Crossings of equal cylinders' ellipses, or a loop's switches
        // between its graphs (S9c.2a).
        let points = match pair {
            CylPair::Crossing(cross) => &cross.points,
            CylPair::Quartic(p) => &p.switches,
            CylPair::Mixed(p) => &p.switches,
            _ => continue,
        };
        for (k, x) in points.iter().enumerate() {
            match (models[0].in_face(*fa, x), models[1].in_face(*fb, x)) {
                (Loc::In, Loc::In) => {
                    vx.push(Vx {
                        p: x.clone(),
                        key: VKey::Cross(*fa, *fb, k),
                        faces: BTreeSet::from([(0, *fa), (1, *fb)]),
                    });
                }
                (Loc::Out, _) | (_, Loc::Out) => {}
                _ => return Err(seam()),
            }
        }
    }
    // Model edges split at their vertices.
    let mut edges: Vec<GEdge> = Vec::new();
    let mut parts_of: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    let mut half: BTreeMap<(usize, usize), Vec<(usize, bool)>> = BTreeMap::new();
    for ((o, ei), list) in &on_edge {
        let e = &models[*o].edges[*ei];
        let mut list = list.clone();
        let ccw = edge_places(&models[*o], *ei).2;
        let start = list[0].1.clone();
        // A height running against the edge (a given `Rise`, S9e.3a) orders
        // its places the other way.
        let cmp = |x: &Pos, y: &Pos| -> Ordering {
            let o = order_on(x, y, &start, ccw);
            if matches!(start, Pos::T(_)) && !ccw {
                o.reverse()
            } else {
                o
            }
        };
        let (first, last) = (list.remove(0), list.remove(0));
        list.sort_by(|x, y| cmp(&x.1, &y.1));
        for w in list.windows(2) {
            if cmp(&w[0].1, &w[1].1) == Ordering::Equal {
                return Err(Error::Degenerate("two meetings at one point of an edge"));
            }
        }
        let mut chain = vec![first];
        chain.extend(list);
        chain.push(last);
        for w in chain.windows(2) {
            let (mid, mid_pos) = midpoint(&e.curve, &w[0].1, &w[1].1, ccw)?;
            let gid = edges.len();
            parts_of.entry((*o, *ei)).or_default().push(gid);
            edges.push(GEdge {
                curve: CurveRef::Edge(*o, *ei),
                crv: e.curve.clone(),
                ends: [w[0].0, w[1].0],
                pos: [w[0].1.clone(), w[1].1.clone()],
                with: ccw,
                mid,
                mid_pos,
            });
            half.entry((*o, e.faces[0])).or_default().push((gid, true));
            half.entry((*o, e.faces[1])).or_default().push((gid, false));
        }
    }
    // A spline joint's edge across a plane's face (S9f.1): the face holds
    // its parts inside it.
    for &(o, ei, g) in &along {
        for &gid in parts_of.get(&(o, ei)).map_or(&[][..], |x| x) {
            match models[1 - o].in_face(g, &edges[gid].mid) {
                Loc::In if touching.contains(&(o, ei, g)) => {
                    return Err(Error::Degenerate(
                        "a plane touching a spline prism along a joint's edge",
                    ))
                }
                Loc::In => {
                    let h = half.entry((1 - o, g)).or_default();
                    h.push((gid, true));
                    h.push((gid, false));
                }
                Loc::On => return Err(Error::Degenerate("edges of both inputs overlapping")),
                Loc::Out => {}
            }
        }
    }
    // Faces on one surface: each holds the other's edges within it.
    for &(fa, fb) in &coinc {
        for (o, own, other) in [(0, fa, fb), (1, fb, fa)] {
            let (me, them) = (&models[o], &models[1 - o]);
            for (ei, e) in them.edges.iter().enumerate() {
                if !e.faces.contains(&other) {
                    continue;
                }
                for &gid in parts_of.get(&(1 - o, ei)).map_or(&[][..], |x| x) {
                    match me.in_face(own, &edges[gid].mid) {
                        Loc::In => {
                            // Once, though it bounds several faces on the
                            // surface (a circle's halves).
                            let h = half.entry((o, own)).or_default();
                            if !h.contains(&(gid, true)) {
                                h.push((gid, true));
                                h.push((gid, false));
                            }
                        }
                        Loc::On => {
                            return Err(Error::Degenerate("edges of both inputs overlapping"))
                        }
                        Loc::Out => {}
                    }
                }
            }
        }
    }
    // Sections.
    let mut secs: Vec<Sec> = Vec::new();
    for fa in 0..models[0].faces.len() {
        for fb in 0..models[1].faces.len() {
            if !boxes_meet(&models[0].boxes[fa], &models[1].boxes[fb]) {
                continue;
            }
            if coinc.contains(&(fa, fb)) {
                continue;
            }
            let pair = pairs.get(&(fa, fb));
            let ((va, ia), (vb, ib)) = (models[0].view(fa), models[1].view(fb));
            let curves = match section(va, ia, vb, ib, pair)? {
                Section::Same => continue,
                Section::Curves(c) => c,
            };
            if curves.is_empty() {
                continue;
            }
            let si = secs.len();
            // A stored sphere's poles on a section, or within the
            // resolution of it, inside both faces give it a vertex (the
            // section's point nearest the pole, exactly on it): a section
            // through a pole turns half a turn there on the stored sphere's
            // parameters, its pcurves meridians' lines on either side. A
            // section through a pole off the axis's planes is refused.
            for (o, f) in [(0, fa), (1, fb)] {
                let (vm, vf) = models[o].view(f);
                let Some(poles) = stored_poles(vm, vf) else {
                    continue;
                };
                for (k, pole) in poles.into_iter().enumerate() {
                    for crv in &curves {
                        let at = match crv {
                            _ if on_curve(crv, &pole) => pole.clone(),
                            Crv::Circle(c) => match near_pole(c, &pole, res) {
                                Some(x) => x,
                                None => continue,
                            },
                            _ => continue,
                        };
                        if models[0].in_face(fa, &at) != Loc::In
                            || models[1].in_face(fb, &at) != Loc::In
                        {
                            continue;
                        }
                        // Only a meridian (a circle in a plane holding the
                        // axis, within rounding) runs through a pole on
                        // lines of its pcurves.
                        let Crv::Circle(c) = crv else {
                            return Err(Error::OutOfDomain(
                                "a section through a sphere's pole off its meridians (S9d.1)",
                            ));
                        };
                        if !meridian(&c.normal(), &vm.f.n) {
                            return Err(Error::OutOfDomain(
                                "a section through a sphere's pole off its meridians (S9d.1)",
                            ));
                        }
                        let id = match vx.iter().position(|v| qv_eq(&v.p, &at)) {
                            Some(id) => id,
                            None => {
                                vx.push(Vx {
                                    p: at,
                                    key: VKey::Pole(si, k),
                                    faces: BTreeSet::new(),
                                });
                                vx.len() - 1
                            }
                        };
                        vx[id].faces.insert((0, fa));
                        vx[id].faces.insert((1, fb));
                    }
                }
            }
            // Vertices on both faces.
            let on: Vec<usize> = (0..vx.len())
                .filter(|&v| vx[v].faces.contains(&(0, fa)) && vx[v].faces.contains(&(1, fb)))
                .collect();
            // Their binary64 views, once per vertex (each is tested against
            // every piece of a meeting).
            if curves.iter().any(|c| matches!(c, Crv::Toric(_))) {
                for &v in &on {
                    if views.len() <= v {
                        views.resize(v + 1, None);
                    }
                    if views[v].is_none() {
                        views[v] = Some(qv_f64(&vx[v].p));
                    }
                }
            }
            for (bi, crv) in curves.iter().enumerate() {
                let mut list: Vec<(usize, Pos)> = on
                    .iter()
                    .filter(|&&v| match crv {
                        // On both faces, so on both surfaces: the piece's
                        // window and range decide (S9d.4b.2).
                        Crv::Toric(c) => match views.get(v).copied().flatten() {
                            Some(view) => c.holds_at(&vx[v].p, view),
                            None => c.holds(&vx[v].p),
                        },
                        _ => on_curve(crv, &vx[v].p),
                    })
                    .map(|&v| (v, place(crv, &vx[v].p)))
                    .collect();
                // Conics and rings are closed; an open piece of a meeting
                // runs from its range's start.
                let (closed, zero_dir) = match crv {
                    Crv::Meet(m) => match &m.range {
                        Some([lo, _]) => (false, Pos::Ang(lo.clone())),
                        None => (true, Pos::Ang([Qd::rat(int(1)), Qd::rat(zero())])),
                    },
                    Crv::Cone(c) => match &c.range {
                        Some([lo, _]) => (false, Pos::Ang(lo.clone())),
                        None => (true, Pos::Ang([Qd::rat(int(1)), Qd::rat(zero())])),
                    },
                    Crv::Torus(c) => match &c.range {
                        Some([lo, _]) => (false, Pos::Ang(lo.clone())),
                        None => (true, Pos::Ang([Qd::rat(int(1)), Qd::rat(zero())])),
                    },
                    Crv::Toric(c) => match &c.range {
                        Some([lo, _]) => (false, Pos::Ang(lo.clone())),
                        None => (true, Pos::Ang([Qd::rat(int(1)), Qd::rat(zero())])),
                    },
                    Crv::Conic { .. } | Crv::Circle(_) => {
                        (true, Pos::Ang([Qd::rat(int(1)), Qd::rat(zero())]))
                    }
                    Crv::Line { .. } | Crv::Rise(_) | Crv::Spline(_) | Crv::WallMeet(_) => {
                        (false, Pos::Ang([Qd::rat(int(1)), Qd::rat(zero())]))
                    }
                };
                list.sort_by(|x, y| order_on(&x.1, &y.1, &zero_dir, true));
                for w in list.windows(2) {
                    if order_on(&w[0].1, &w[1].1, &zero_dir, true) == Ordering::Equal {
                        return Err(Error::Degenerate("two meetings at one point of a section"));
                    }
                }
                let n = list.len();
                let segments: Vec<(usize, usize)> = if closed {
                    if n == 0 {
                        // A closed section with no vertex: never inside both
                        // faces for cylinders (their seams cross their
                        // sections); a sphere's circle inside both faces is
                        // a ring, given one vertex at a rational point
                        // (S9d.1).
                        let x = match crv {
                            Crv::Meet(m) => m.at(&[int(1), zero()]).ok_or(
                                Error::ComputationLimit("a closed section without vertices"),
                            )?,
                            Crv::Circle(c) => c.at(&[int(1), zero()]),
                            Crv::Cone(c) => c.at(&[int(1), zero()]).ok_or(
                                Error::ComputationLimit("a closed section without vertices"),
                            )?,
                            Crv::Torus(c) => c.at(&[int(1), zero()]).ok_or(
                                Error::ComputationLimit("a closed section without vertices"),
                            )?,
                            Crv::Toric(c) => c.at(&[int(1), zero()]).ok_or(
                                Error::ComputationLimit("a closed section without vertices"),
                            )?,
                            _ => conic_point_r(crv, &[int(1), zero()]),
                        };
                        let inside = models[0].in_face(fa, &x) == Loc::In
                            && models[1].in_face(fb, &x) == Loc::In;
                        match (inside, crv) {
                            (false, _) => Vec::new(),
                            // A ring on a cone's wall (no seam crosses it,
                            // S9d.3b) likewise.
                            (
                                true,
                                Crv::Circle(_)
                                | Crv::Cone(_)
                                | Crv::Meet(_)
                                | Crv::Torus(_)
                                | Crv::Toric(_),
                            ) => {
                                let v = vx.len();
                                vx.push(Vx {
                                    p: x.clone(),
                                    key: VKey::Ring(si, bi),
                                    faces: BTreeSet::from([(0, fa), (1, fb)]),
                                });
                                list.push((v, place(crv, &x)));
                                vec![(0, 0)]
                            }
                            (true, _) => {
                                return Err(Error::ComputationLimit(
                                    "a closed section without vertices",
                                ))
                            }
                        }
                    } else {
                        (0..n).map(|i| (i, (i + 1) % n)).collect()
                    }
                } else {
                    (0..n.saturating_sub(1)).map(|i| (i, i + 1)).collect()
                };
                for (i, j) in segments {
                    let (a, b) = (&list[i], &list[j]);
                    let (mid, mid_pos) = midpoint(crv, &a.1, &b.1, true)?;
                    let la = models[0].in_face(fa, &mid);
                    let lb = models[1].in_face(fb, &mid);
                    match (la, lb) {
                        (Loc::In, Loc::In) => {}
                        (Loc::Out, _) | (_, Loc::Out) => continue,
                        // Along an edge of a face on the other's surface:
                        // that face holds it.
                        _ if in_coinc[0].contains(&fa) || in_coinc[1].contains(&fb) => continue,
                        _ => return Err(seam()),
                    }
                    let gid = edges.len();
                    edges.push(GEdge {
                        curve: CurveRef::Section(si, bi),
                        crv: crv.clone(),
                        ends: [a.0, b.0],
                        pos: [a.1.clone(), b.1.clone()],
                        with: true,
                        mid,
                        mid_pos,
                    });
                    for key in [(0, fa), (1, fb)] {
                        let h = half.entry(key).or_default();
                        h.push((gid, true));
                        h.push((gid, false));
                    }
                }
            }
            secs.push(Sec { fa, fb });
        }
    }
    let mut arr = Arr {
        models,
        vx,
        edges,
        secs,
        pieces: Vec::new(),
        coinc,
    };
    // Pieces.
    for ((o, f), hs) in &half {
        let loops = arr.trace(*o, *f, hs)?;
        let pieces = arr.group(*o, *f, loops)?;
        for loops in pieces {
            let sides = arr.sides(*o, *f, &loops[0])?;
            arr.pieces.push(Piece {
                op: *o,
                face: *f,
                loops,
                sides,
                keep: false,
                behind: false,
            });
        }
    }
    Ok(arr)
}

/// Spline joints' edges on planes of the other input, and those touching.
type Joints = (
    BTreeSet<(usize, usize, usize)>,
    BTreeSet<(usize, usize, usize)>,
);

/// S9f.1: each spline prism's vertical edge between two spline walls (its
/// joint) lying in the plane of a face of the other input, the plane's
/// face meeting it: `(operand, edge, face)`. Its ends must lie off that
/// face (a vertex of one input on the other's face is `Degenerate`).
///
/// The plane holds the axis, so its trace in the profile is a line through
/// the joint: it crosses the boundary there when the segments arriving and
/// leaving lie on its two sides (their end tangents strictly on one side of
/// its normal, exactly), and touches it otherwise (`touching`: refused
/// where the edge meets the face, a contact of the inputs along an edge).
fn joint_edges_on_planes(models: &[Prism; 2]) -> Result<Joints> {
    let mut out = BTreeSet::new();
    let mut touching = BTreeSet::new();
    for o in 0..2 {
        let (me, other) = (&models[o], &models[1 - o]);
        if me.given.is_some() {
            continue;
        }
        for (ei, e) in me.edges.iter().enumerate() {
            let (EdgeKind::Vertical(b, j), Crv::Line { p, d }) = (e.kind, &e.curve) else {
                continue;
            };
            if !e
                .faces
                .iter()
                .all(|&f| matches!(me.faces[f].surf, Surf::Spline(_)))
            {
                continue;
            }
            let ebox = intersect(&me.boxes[e.faces[0]], &me.boxes[e.faces[1]]);
            for g in 0..other.faces.len() {
                let (vo, vg) = other.view(g);
                let Surf::Plane { p: p0, m } = &vo.faces[vg].surf else {
                    continue;
                };
                if !boxes_meet(&ebox, &other.boxes[g]) || dot(m, d) != zero() {
                    continue;
                }
                if qdot(&qsub(&qv(p0), p), m).sign() != Ordering::Equal {
                    continue;
                }
                // Its ends off the face.
                for v in [e.start, e.end] {
                    if other.in_face(g, &me.verts[v].p) != Loc::Out {
                        return Err(Error::Degenerate(
                            "a vertex of one input on the other's face",
                        ));
                    }
                }
                // The segments' sides of the plane's trace at the joint.
                let segs = &me.bounds[b].segs;
                let (Seg::Spline(arriving), Seg::Spline(leaving)) =
                    (&segs[(j + segs.len() - 1) % segs.len()], &segs[j])
                else {
                    unreachable!("a joint of two spline walls")
                };
                let trace = [dot(m, &me.f.x), dot(m, &me.f.y)];
                let side = |t: [R; 2]| sign(&(&trace[0] * &t[0] + &trace[1] * &t[1]));
                let (a, l) = (
                    side(arriving.end_tangent(false)),
                    side(leaving.end_tangent(true)),
                );
                if a == Ordering::Equal || l == Ordering::Equal || a != l {
                    touching.insert((o, ei, g));
                }
                out.insert((o, ei, g));
            }
        }
    }
    Ok((out, touching))
}

/// The poles of a sphere face's stored surface (`c -+ r n / |n|` on the
/// stored frame's axis), exactly; none for other faces.
fn stored_poles(m: &Prism, f: usize) -> Option<Vec<QV>> {
    let (Surf::Sphere { c, r }, Some(_)) = (&m.faces[f].surf, &m.ball) else {
        return None;
    };
    let n = &m.f.n;
    let nn = dot(n, n);
    Some(
        [-1, 1]
            .map(|k| {
                let s = Qd::new(zero(), int(k), r * r / &nn);
                qadd(&qv(c), &qscale(n, &s))
            })
            .to_vec(),
    )
}

/// A circle's point nearest a pole within the resolution of it, exactly on
/// the circle (its rational direction nearest the pole's), or none.
fn near_pole(c: &super::sphere::Circ, pole: &QV, res: f64) -> Option<QV> {
    let fl = |v: &V| v.clone().map(|x| crate::solid::split::rational_f64(&x));
    let (cc, n) = (fl(&c.c), fl(&c.normal()));
    let p = qv_f64(pole);
    let d = [p[0] - cc[0], p[1] - cc[1], p[2] - cc[2]];
    let nn = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    let h = (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]) / nn;
    let inplane = [0, 1, 2].map(|i| d[i] - h * n[i] / nn);
    let rad = (inplane[0].powi(2) + inplane[1].powi(2) + inplane[2].powi(2)).sqrt();
    let r = crate::solid::split::rational_f64(&c.r2).sqrt();
    if (h * h + (rad - r) * (rad - r)).sqrt() > res {
        return None;
    }
    let place = c.place(pole);
    Some(c.at(&[q(place[0].to_f64()), q(place[1].to_f64())]))
}

/// Whether a plane's normal is orthogonal to an axis within rounding (the
/// plane holds the axis's direction).
fn meridian(m: &V, n: &V) -> bool {
    let fl = |v: &V| v.clone().map(|x| crate::solid::split::rational_f64(&x));
    let (m, n) = (fl(m), fl(n));
    let d = m[0] * n[0] + m[1] * n[1] + m[2] * n[2];
    let l = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt()
        * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    d.abs() <= 1e-12 * l
}

/// The edge of face `g` holding `x` strictly inside it, and its place
/// there (`Degenerate` at an edge's end).
fn edge_at(m: &Prism, g: usize, x: &QV) -> Result<Option<(usize, Pos)>> {
    for (fi, f) in m.edges.iter().enumerate() {
        if !f.faces.contains(&g) || !on_curve(&f.curve, x) {
            continue;
        }
        let pos = place(&f.curve, x);
        let (ps, pt, ccw) = edge_places(m, fi);
        if same_end(&pos, &ps) || same_end(&pos, &pt) {
            return Err(Error::Degenerate(
                "a vertex of one input on an edge of the other",
            ));
        }
        if strictly_within(&pos, &ps, &pt, ccw) == Some(true) {
            return Ok(Some((fi, pos)));
        }
    }
    Ok(None)
}

fn intersect(a: &([f64; 3], [f64; 3]), b: &([f64; 3], [f64; 3])) -> ([f64; 3], [f64; 3]) {
    let mut out = *a;
    for k in 0..3 {
        out.0[k] = a.0[k].max(b.0[k]);
        out.1[k] = a.1[k].min(b.1[k]);
    }
    out
}

/// Whether a meeting at `x` on face `g` may be a seam's: on a full circle's
/// half wall anywhere, on a sphere's hemisphere only on its split (S9d.2).
pub(super) fn seam_at(m: &Prism, g: usize, x: &QV) -> bool {
    let (m, g) = m.view(g);
    if let Some(fun) = &m.funnel {
        // A cone's faces: on the rims' seam direction (S9d.3a).
        return fun.on_seam(&m.f, x);
    }
    if let Some(ring) = &m.ring {
        // A torus's patches: on a meridian or parallel seam (S9d.4a).
        return ring.on_seam(&m.f, x);
    }
    match (m.faces[g].kind, &m.ball) {
        (FaceKind::Half(_), Some(ball)) => {
            qdot(&qsub(x, &qv(&ball.c)), &ball.split).sign() == Ordering::Equal
        }
        _ => on_circle(m, g),
    }
}

/// Whether a face is a full circle's half wall (its sides are a seam).
fn on_circle(m: &Prism, f: usize) -> bool {
    match m.faces[f].kind {
        FaceKind::Wall(b, _) => m.bounds[b].circle,
        FaceKind::Cap(_) | FaceKind::ConeWall => false,
        // A hemisphere's sides are the split, a torus patch's its seams.
        FaceKind::Patch(..) => true,
        FaceKind::Half(_) => true,
    }
}

fn same_end(pos: &Pos, end: &Pos) -> bool {
    match (pos, end) {
        (Pos::T(a), Pos::T(b)) => a.cmp(b) == Ordering::Equal,
        (Pos::Ang(a), Pos::Ang(b)) => same_dir(a, b),
        _ => false,
    }
}

/// The order of two places along a curve from `start` (lines by key,
/// conics by the angle turned from `start` with or against the angle).
fn order_on(x: &Pos, y: &Pos, start: &Pos, ccw: bool) -> Ordering {
    match (x, y, start) {
        (Pos::T(a), Pos::T(b), _) => a.cmp(b),
        (Pos::Ang(a), Pos::Ang(b), Pos::Ang(s)) => {
            if same_dir(a, b) {
                return Ordering::Equal;
            }
            if same_dir(a, s) {
                return Ordering::Less;
            }
            if same_dir(b, s) {
                return Ordering::Greater;
            }
            // a before b when a lies between the start and b.
            if between_run(s, a, b, ccw) {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        }
        _ => unreachable!("places of one kind"),
    }
}

fn conic_point_r(crv: &Crv, cs: &[R; 2]) -> QV {
    let Crv::Conic { c, a, b } = crv else {
        unreachable!("a conic")
    };
    conic_point(c, a, b, &[Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())])
}

/// A point strictly between two places on a curve, and its place.
fn midpoint(crv: &Crv, a: &Pos, b: &Pos, ccw: bool) -> Result<(QV, Pos)> {
    match (crv, a, b) {
        (Crv::Line { p, d }, Pos::T(ta), Pos::T(tb)) => {
            let (lo, hi) = if ta.cmp(tb) == Ordering::Less {
                (ta, tb)
            } else {
                (tb, ta)
            };
            let k = rational_between_num(lo, hi)?;
            // p + (k - key(p)) d.
            let t = Qd::rat(k.clone()).sub(&line_key(p, d));
            let x = qadd(p, &qscale(d, &t));
            Ok((x, Pos::T(Qd::rat(k))))
        }
        // S9f.1: a rational run parameter, so a rational point.
        (Crv::Spline(c), Pos::T(ta), Pos::T(tb)) => {
            let (lo, hi) = if ta.cmp(tb) == Ordering::Less {
                (ta, tb)
            } else {
                (tb, ta)
            };
            let k = Qd::rat(rational_between_num(lo, hi)?);
            Ok((c.point(&k), Pos::T(k)))
        }
        (Crv::Rise(c), Pos::T(ta), Pos::T(tb)) => {
            let (lo, hi) = if ta.cmp(tb) == Ordering::Less {
                (ta, tb)
            } else {
                (tb, ta)
            };
            let k = rational_between_num(lo, hi)?;
            let x = c
                .at(&k)?
                .ok_or(Error::ComputationLimit("a meeting's point off its piece"))?;
            Ok((x, Pos::T(Qd::rat(k))))
        }
        // S9f.2b: a rational run parameter, its point in `Q(sqrt(D))`.
        (Crv::WallMeet(c), Pos::T(ta), Pos::T(tb)) => {
            let (lo, hi) = if ta.cmp(tb) == Ordering::Less {
                (ta, tb)
            } else {
                (tb, ta)
            };
            let k = rational_between_num(lo, hi)?;
            let x = c
                .at(&k)
                .ok_or(Error::ComputationLimit("a meeting's point off its piece"))?;
            Ok((x, Pos::T(Qd::rat(k))))
        }
        (Crv::Circle(c), Pos::Ang(sa), Pos::Ang(sb)) => {
            let cs = rational_between(sa, sb, ccw)?;
            let x = c.at(&cs);
            Ok((
                x,
                Pos::Ang([Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]),
            ))
        }
        (Crv::Torus(m), Pos::Ang(sa), Pos::Ang(sb)) => {
            let cs = rational_between(sa, sb, ccw)?;
            let x = m
                .at(&cs)
                .ok_or(Error::ComputationLimit("a torus section's point off it"))?;
            Ok((
                x,
                Pos::Ang([Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]),
            ))
        }
        (Crv::Toric(m), Pos::Ang(sa), Pos::Ang(sb)) => {
            let cs = rational_between(sa, sb, ccw)?;
            let x = m
                .at(&cs)
                .ok_or(Error::ComputationLimit("a torus meeting's point off it"))?;
            Ok((
                x,
                Pos::Ang([Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]),
            ))
        }
        (Crv::Cone(m), Pos::Ang(sa), Pos::Ang(sb)) => {
            let cs = rational_between(sa, sb, ccw)?;
            let x = m
                .at(&cs)
                .ok_or(Error::ComputationLimit("a cone section's point off it"))?;
            Ok((
                x,
                Pos::Ang([Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]),
            ))
        }
        (Crv::Meet(m), Pos::Ang(sa), Pos::Ang(sb)) => {
            let cs = rational_between(sa, sb, ccw)?;
            let x = m
                .at(&cs)
                .ok_or(Error::ComputationLimit("a meeting's point off its piece"))?;
            Ok((
                x,
                Pos::Ang([Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]),
            ))
        }
        (Crv::Conic { .. }, Pos::Ang(sa), Pos::Ang(sb)) => {
            let cs = rational_between(sa, sb, ccw)?;
            let x = conic_point_r(crv, &cs);
            Ok((
                x,
                Pos::Ang([Qd::rat(cs[0].clone()), Qd::rat(cs[1].clone())]),
            ))
        }
        _ => unreachable!("places of the curve's kind"),
    }
}

/// Whether a point lies on the line `p + t d` (exactly).
pub(super) fn on_line(p: &QV, d: &V, x: &QV) -> bool {
    on_curve(
        &Crv::Line {
            p: p.clone(),
            d: d.clone(),
        },
        x,
    )
}

/// Whether a point lies on a curve (exactly).
fn on_curve(crv: &Crv, x: &QV) -> bool {
    match crv {
        Crv::Line { p, d } => {
            // (x - p) x d = 0, compared term by term (fields may differ).
            let xc = [
                x[1].scale(&d[2]).sub(&x[2].scale(&d[1])),
                x[2].scale(&d[0]).sub(&x[0].scale(&d[2])),
                x[0].scale(&d[1]).sub(&x[1].scale(&d[0])),
            ];
            let pc = [
                p[1].scale(&d[2]).sub(&p[2].scale(&d[1])),
                p[2].scale(&d[0]).sub(&p[0].scale(&d[2])),
                p[0].scale(&d[1]).sub(&p[1].scale(&d[0])),
            ];
            xc.iter().zip(&pc).all(|(a, b)| a.cmp(b) == Ordering::Equal)
        }
        Crv::Meet(m) => m.on(x),
        Crv::Cone(c) => c.on(x),
        Crv::Torus(c) => c.on(x),
        Crv::Toric(c) => c.on(x),
        Crv::Circle(c) => c.on(x),
        Crv::Rise(c) => c.on(x),
        Crv::Spline(c) => c.on(x),
        Crv::WallMeet(c) => c.on(x),
        Crv::Conic { c, a, b } => {
            let cs = conic_angle(c, a, b, x);
            let back = conic_point(c, a, b, &cs);
            let unit = cs[0].mul(&cs[0]).add(&cs[1].mul(&cs[1])).add_r(&int(-1));
            qv_eq(&back, x) && unit.sign() == Ordering::Equal
        }
    }
}

impl Arr {
    /// A half-edge's direction of travel at a place (unit-free).
    fn travel(&self, gid: usize, fwd: bool, pos: &Pos, x: &QV) -> QV {
        let e = &self.edges[gid];
        let t = tangent(&e.crv, pos, x);
        if e.with == fwd {
            t
        } else {
            t.map(|x| x.neg())
        }
    }

    fn start_of(&self, h: (usize, bool)) -> usize {
        let e = &self.edges[h.0];
        if h.1 {
            e.ends[0]
        } else {
            e.ends[1]
        }
    }

    fn end_of(&self, h: (usize, bool)) -> usize {
        self.start_of((h.0, !h.1))
    }

    fn pos_at_start(&self, h: (usize, bool)) -> &Pos {
        &self.edges[h.0].pos[if h.1 { 0 } else { 1 }]
    }

    /// The face's loops: from each half-edge, the next is the first
    /// clockwise (about the face's outward normal) from the way back.
    fn trace(&self, o: usize, f: usize, hs: &[(usize, bool)]) -> Result<Vec<Vec<(usize, bool)>>> {
        let mut out_of: BTreeMap<usize, Vec<(usize, bool)>> = BTreeMap::new();
        for &h in hs {
            out_of.entry(self.start_of(h)).or_default().push(h);
        }
        let mut used: BTreeSet<(usize, bool)> = BTreeSet::new();
        let mut loops = Vec::new();
        for &h0 in hs {
            if used.contains(&h0) {
                continue;
            }
            let mut lp = vec![h0];
            used.insert(h0);
            let mut h = h0;
            loop {
                let v = self.end_of(h);
                let next = self.next_on(o, f, v, h, out_of.get(&v).map_or(&[][..], |x| x))?;
                if next == h0 {
                    break;
                }
                if !used.insert(next) {
                    return Err(Error::InvalidTopology("a face's pieces do not close"));
                }
                lp.push(next);
                h = next;
            }
            loops.push(lp);
        }
        Ok(loops)
    }

    /// The outgoing half-edge after arriving at `v` along `h`.
    pub(super) fn next_on(
        &self,
        o: usize,
        f: usize,
        v: usize,
        h: (usize, bool),
        outs: &[(usize, bool)],
    ) -> Result<(usize, bool)> {
        let back = (h.0, !h.1);
        let cands: Vec<(usize, bool)> = outs.iter().copied().filter(|&c| c != back).collect();
        if cands.is_empty() {
            return if outs.contains(&back) {
                Ok(back)
            } else {
                Err(Error::InvalidTopology("a face's piece ends at a vertex"))
            };
        }
        if cands.len() == 1 {
            return Ok(cands[0]);
        }
        let p = &self.vx[v].p;
        let n = self.models[o].normal_at(f, p);
        let r = self.travel(back.0, back.1, self.pos_at_start(back), p);
        let nr = qcross(&n, &r);
        let coords = |c: (usize, bool)| -> [Qd; 2] {
            let t = self.travel(c.0, c.1, self.pos_at_start(c), p);
            [qqdot(&r, &t), qqdot(&nr, &t)]
        };
        let zero_dir = [Qd::rat(int(1)), Qd::rat(zero())];
        let mut best: Option<((usize, bool), [Qd; 2])> = None;
        for c in cands {
            let a = coords(c);
            if same_dir(&a, &zero_dir) {
                return Err(Error::Degenerate("two curves tangent at a vertex"));
            }
            best = match best {
                None => Some((c, a)),
                Some((bc, ba)) => match angle_cmp(&a, &ba) {
                    Ordering::Greater => Some((c, a)),
                    Ordering::Less => Some((bc, ba)),
                    Ordering::Equal => {
                        return Err(Error::Degenerate("two curves tangent at a vertex"))
                    }
                },
            };
        }
        Ok(best.expect("a candidate").0)
    }

    /// Points of a half-edge in binary64, along its run.
    pub(super) fn samples(&self, h: (usize, bool)) -> Vec<[f64; 3]> {
        let e = &self.edges[h.0];
        match &e.crv {
            Crv::Line { .. } => {
                let (a, b) = (qv_f64(&self.vx[e.ends[0]].p), qv_f64(&self.vx[e.ends[1]].p));
                if h.1 {
                    vec![a, b]
                } else {
                    vec![b, a]
                }
            }
            Crv::Spline(c) => {
                let (Pos::T(t0), Pos::T(t1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a spline curve's places")
                };
                let mut pts = c.samples(t0.to_f64(), t1.to_f64(), 32);
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
            Crv::WallMeet(c) => {
                let (Pos::T(t0), Pos::T(t1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a spline wall's meeting's places")
                };
                let mut pts = c.samples(t0.to_f64(), t1.to_f64(), 48);
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
            Crv::Rise(c) => {
                let (Pos::T(w0), Pos::T(w1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a rise's places")
                };
                let mut pts = c.samples(w0.to_f64(), w1.to_f64(), 48);
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
            Crv::Circle(c) => {
                let (Pos::Ang(p0), Pos::Ang(p1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a circle's places")
                };
                let (t0, t1) = (c.angle(p0), c.angle(p1));
                let mut sweep = if e.with { t1 - t0 } else { t0 - t1 };
                sweep = sweep.rem_euclid(TAU);
                if sweep == 0.0 {
                    sweep = TAU;
                }
                let sweep = if e.with { sweep } else { -sweep };
                let mut pts = c.samples(t0, sweep, 96);
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
            Crv::Torus(m) => {
                let (Pos::Ang(p0), Pos::Ang(p1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a torus section's places")
                };
                let (t0, t1) = (angle_f64(p0), angle_f64(p1));
                let mut sweep = if e.with { t1 - t0 } else { t0 - t1 };
                sweep = sweep.rem_euclid(TAU);
                if sweep == 0.0 {
                    sweep = TAU;
                }
                let sweep = if e.with { sweep } else { -sweep };
                let mut pts = m.samples(t0, sweep, 64);
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
            Crv::Toric(m) => {
                let (Pos::Ang(p0), Pos::Ang(p1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a torus meeting's places")
                };
                let (t0, t1) = (angle_f64(p0), angle_f64(p1));
                let mut sweep = if e.with { t1 - t0 } else { t0 - t1 };
                sweep = sweep.rem_euclid(TAU);
                if sweep == 0.0 {
                    sweep = TAU;
                }
                let sweep = if e.with { sweep } else { -sweep };
                let mut pts = m.samples(t0, sweep, 64);
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
            Crv::Cone(m) => {
                let (Pos::Ang(p0), Pos::Ang(p1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a cone section's places")
                };
                let (t0, t1) = (angle_f64(p0), angle_f64(p1));
                let mut sweep = if e.with { t1 - t0 } else { t0 - t1 };
                sweep = sweep.rem_euclid(TAU);
                if sweep == 0.0 {
                    sweep = TAU;
                }
                let sweep = if e.with { sweep } else { -sweep };
                let mut pts = m.samples(t0, sweep, 48);
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
            Crv::Meet(m) => {
                let (Pos::Ang(p0), Pos::Ang(p1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a meeting's places")
                };
                let (t0, t1) = (angle_f64(p0), angle_f64(p1));
                let mut sweep = if e.with { t1 - t0 } else { t0 - t1 };
                sweep = sweep.rem_euclid(TAU);
                if sweep == 0.0 {
                    sweep = TAU;
                }
                let sweep = if e.with { sweep } else { -sweep };
                let mut pts = m.samples(t0, sweep, 24);
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
            Crv::Conic { c, a, b } => {
                let (Pos::Ang(p0), Pos::Ang(p1)) = (&e.pos[0], &e.pos[1]) else {
                    unreachable!("a conic's places")
                };
                let (t0, t1) = (angle_f64(p0), angle_f64(p1));
                let mut sweep = if e.with { t1 - t0 } else { t0 - t1 };
                sweep = sweep.rem_euclid(TAU);
                if sweep == 0.0 {
                    sweep = TAU;
                }
                let sweep = if e.with { sweep } else { -sweep };
                let f = |x: &V| x.clone().map(|y| crate::solid::split::rational_f64(&y));
                let (cf, af, bf) = (f(c), f(a), f(b));
                let n = 24;
                let mut pts: Vec<[f64; 3]> = (0..=n)
                    .map(|i| {
                        let t = t0 + sweep * i as f64 / n as f64;
                        [0, 1, 2].map(|k| cf[k] + af[k] * t.cos() + bf[k] * t.sin())
                    })
                    .collect();
                if !h.1 {
                    pts.reverse();
                }
                pts
            }
        }
    }

    /// A face's binary64 parameters of a point: the plane's frame
    /// coordinates, or a cylinder's turn from its arc's start and height;
    /// and the orientation of those parameters against the outward normal.
    pub(super) fn params(&self, o: usize, f: usize) -> impl Fn([f64; 3]) -> [f64; 2] + '_ {
        let (m, lf) = self.models[o].view(f);
        let face = &m.faces[lf];
        let fl = |x: &V| x.clone().map(|y| crate::solid::split::rational_f64(&y));
        let (of, xf, yf, nf) = (fl(&m.f.o), fl(&m.f.x), fl(&m.f.y), fl(&m.f.n));
        let surf = face.surf.clone();
        let kind = face.kind;
        let ball = m.ball.clone();
        let ring = m.ring.as_ref().map(|r| r.params(&m.f));
        let lo = rational_f64(&m.lo);
        move |p: [f64; 3]| -> [f64; 2] {
            match &surf {
                // A spline wall (S9f.1): its stored surface's parameters,
                // the curve's own parameter and the height above the low
                // cap.
                Surf::Spline(s) => {
                    let d = [p[0] - of[0], p[1] - of[1], p[2] - of[2]];
                    let l = solve3(&xf, &yf, &nf, &d);
                    let tau = s.tau_near(l[0], l[1]);
                    let t = if s.rev {
                        rational_f64(&s.first) + rational_f64(&s.last) - tau
                    } else {
                        tau
                    };
                    [t, l[2] - lo]
                }
                // A torus's patch: its angles from its seams (S9d.4a).
                Surf::Torus => ring.as_ref().expect("a torus's ring")(p, kind),
                // A cone's wall: its projection on the plane of `(u, v)`,
                // one to one (S9d.3a).
                Surf::Cone { .. } => {
                    let d = [p[0] - of[0], p[1] - of[1], p[2] - of[2]];
                    let l = solve3(&xf, &yf, &nf, &d);
                    [l[0], l[1]]
                }
                Surf::Sphere { .. } => {
                    let FaceKind::Half(side) = kind else {
                        unreachable!("a hemisphere")
                    };
                    ball.as_ref().expect("a sphere's ball").params(p, side)
                }
                Surf::Plane { .. } => {
                    let crate::topology::Surface::Plane(frame) = &face.stored else {
                        unreachable!("a plane face")
                    };
                    let c = frame.coordinates(crate::Point3::new(p[0], p[1], p[2]));
                    [c[0], c[1]]
                }
                Surf::Cyl { c, r, .. } => {
                    // Local coordinates by the (nearly orthonormal) axes'
                    // Gram solve.
                    let d = [p[0] - of[0], p[1] - of[1], p[2] - of[2]];
                    let l = solve3(&xf, &yf, &nf, &d);
                    let (cu, cv) = (
                        crate::solid::split::rational_f64(&c[0]),
                        crate::solid::split::rational_f64(&c[1]),
                    );
                    let rf = crate::solid::split::rational_f64(r);
                    let theta = (l[1] - cv).atan2(l[0] - cu);
                    let FaceKind::Wall(b, j) = kind else {
                        unreachable!("a wall")
                    };
                    let Seg::Arc {
                        c: ac,
                        p: ap,
                        q: aq,
                        ccw,
                        ..
                    } = &m.bounds[b].segs[j]
                    else {
                        unreachable!("an arc wall")
                    };
                    let fv = crate::solid::split::rational_f64;
                    let sweep = crate::profile::arc_sweep(
                        crate::Point2::new(fv(&ac[0]), fv(&ac[1])),
                        crate::Point2::new(fv(&ap[0]), fv(&ap[1])),
                        crate::Point2::new(fv(&aq[0]), fv(&aq[1])),
                        *ccw,
                    )
                    .abs();
                    let a0 = (crate::solid::split::rational_f64(&ap[1])
                        - crate::solid::split::rational_f64(&ac[1]))
                    .atan2(
                        crate::solid::split::rational_f64(&ap[0])
                            - crate::solid::split::rational_f64(&ac[0]),
                    );
                    let turn = if *ccw { theta - a0 } else { a0 - theta };
                    // Within the arc's sweep: a point at its start may turn
                    // a rounding below zero.
                    let mut t = turn.rem_euclid(TAU);
                    if t > (TAU + sweep) / 2.0 {
                        t -= TAU;
                    }
                    [t * rf, l[2]]
                }
            }
        }
    }

    /// The sign turning the parameters' area into the area about the face's
    /// outward normal.
    pub(super) fn param_sign(&self, o: usize, f: usize) -> f64 {
        let (m, lf) = self.models[o].view(f);
        // The primitive face's orientation, its frame and data, reversed
        // where a given result reverses the face (S9e.3a: through every
        // level, a sphere's, cone's or torus's face too).
        let face = &m.faces[lf];
        let sign = if self.models[o].reversed(f) {
            -1.0
        } else {
            1.0
        };
        sign * match &face.surf {
            // A spline wall's `(t, v)` normal is `S_t x n`, the run's right
            // when the run follows `t`: outward on an outer boundary run
            // with `t` or a hole's run against it (S8b.2's walls).
            Surf::Spline(s) => {
                if s.rev == s.hole {
                    1.0
                } else {
                    -1.0
                }
            }
            // `(u, v)` runs counter-clockwise about the axis: the wall's
            // outward normal leans along it where the cone narrows upward.
            // `(u, v)` runs with the torus's outward normal (against an
            // inside-out segment's).
            Surf::Torus => {
                if m.ring.as_ref().is_some_and(|r| r.reversed) {
                    -1.0
                } else {
                    1.0
                }
            }
            Surf::Cone { k, .. } => {
                if *k < zero() {
                    1.0
                } else {
                    -1.0
                }
            }
            // The projection's coordinates run counter-clockwise about the
            // split's normal, a hemisphere's outward normal on its side.
            Surf::Sphere { .. } => {
                if face.kind == FaceKind::Half(true) {
                    1.0
                } else {
                    -1.0
                }
            }
            Surf::Plane { m: out, .. } => {
                let crate::topology::Surface::Plane(frame) = &face.stored else {
                    unreachable!("a plane face")
                };
                let n = frame.normal().to_array();
                let o = out.clone().map(|y| crate::solid::split::rational_f64(&y));
                if n[0] * o[0] + n[1] * o[1] + n[2] * o[2] > 0.0 {
                    1.0
                } else {
                    -1.0
                }
            }
            Surf::Cyl { inside, .. } => {
                let FaceKind::Wall(b, j) = face.kind else {
                    unreachable!("a wall")
                };
                let Seg::Arc { ccw, .. } = &m.bounds[b].segs[j] else {
                    unreachable!("an arc wall")
                };
                // (turn, height) runs with (angle, height) on a
                // counter-clockwise arc, whose frame's area is about the
                // radial outward normal.
                let s = if *ccw { 1.0 } else { -1.0 };
                if *inside {
                    s
                } else {
                    -s
                }
            }
        }
    }

    /// A loop's binary64 polygon in the face's parameters.
    pub(super) fn polygon(&self, o: usize, f: usize, lp: &[(usize, bool)]) -> Vec<[f64; 2]> {
        let par = self.params(o, f);
        let mut out = Vec::new();
        for &h in lp {
            let s = self.samples(h);
            for p in &s[..s.len() - 1] {
                out.push(par(*p));
            }
        }
        out
    }

    /// Loops grouped into pieces: each counter-clockwise loop (about the
    /// outward normal) with the clockwise ones it holds nearest.
    fn group(&self, o: usize, f: usize, loops: Vec<Vec<(usize, bool)>>) -> Result<Vec<Vec<HLoop>>> {
        let sign = self.param_sign(o, f);
        let polys: Vec<Vec<[f64; 2]>> = loops.iter().map(|l| self.polygon(o, f, l)).collect();
        let areas: Vec<f64> = polys.iter().map(|p| sign * area(p)).collect();
        let scale = polys
            .iter()
            .flatten()
            .fold(1.0f64, |m, p| m.max(p[0].abs()).max(p[1].abs()));
        let tiny = 1e-20 * scale * scale;
        let mut outers: Vec<usize> = Vec::new();
        let mut holes: Vec<usize> = Vec::new();
        for (i, a) in areas.iter().enumerate() {
            if a.abs() <= tiny {
                return Err(Error::Degenerate("a piece thinner than the resolution"));
            }
            if *a > 0.0 {
                outers.push(i);
            } else {
                holes.push(i);
            }
        }
        let mut pieces: Vec<Vec<Vec<(usize, bool)>>> =
            outers.iter().map(|&i| vec![loops[i].clone()]).collect();
        for h in holes {
            let probe = polys[h][0];
            let mut best: Option<(usize, f64)> = None;
            let hole_edges: BTreeSet<usize> = loops[h].iter().map(|x| x.0).collect();
            for (k, &i) in outers.iter().enumerate() {
                // The piece the hole bounds from inside shares its edges.
                if loops[i].iter().any(|x| hole_edges.contains(&x.0)) {
                    continue;
                }
                let (inside, dist) = winding(&polys[i], probe);
                if dist <= 1e-9 * scale {
                    return Err(Error::Degenerate("a hole touching a piece's boundary"));
                }
                if inside && best.is_none_or(|(_, a)| areas[i] < a) {
                    best = Some((k, areas[i]));
                }
            }
            // A hole no piece holds: a loop whose binary64 image turned
            // the wrong way (a sliver within the resolution).
            let (k, _) = best.ok_or(Error::Degenerate("a piece thinner than the resolution"))?;
            pieces[k].push(loops[h].clone());
        }
        Ok(pieces)
    }

    /// Every piece's keeping and side for an operation (`classify`).
    pub(super) fn for_op(&mut self, op: Op2) {
        for piece in &mut self.pieces {
            (piece.keep, piece.behind) = classify(piece.op, piece.sides, op);
        }
    }

    /// Whether the other solid holds a piece's front and its back: its
    /// membership at a point of the piece's first edge pushed into it, then
    /// off the face either way.
    fn sides(&self, o: usize, f: usize, outer: &[(usize, bool)]) -> Result<(bool, bool)> {
        let h = outer[0];
        let e = &self.edges[h.0];
        let x = &e.mid;
        let t = self.travel(h.0, h.1, &e.mid_pos, x);
        let n = self.models[o].normal_at(f, x);
        let inward = qcross(&n, &t);
        let neg_n = n.clone().map(|y| y.neg());
        let other = &self.models[1 - o];
        let front = other.member(x, &[inward.clone(), n]);
        let back = other.member(x, &[inward, neg_n]);
        let (front, back) = match (front, back) {
            (Loc::On, _) | (_, Loc::On) => {
                return Err(Error::Degenerate("a piece on the other's boundary"))
            }
            (a, b) => (a == Loc::In, b == Loc::In),
        };
        Ok((front, back))
    }
}

/// Whether the result keeps a piece of operand `o` whose front and back
/// the other solid holds as `(front, back)`, and whether its material lies
/// behind the face.
pub(super) fn classify(o: usize, (front, back): (bool, bool), op: Op2) -> (bool, bool) {
    if o == 0 {
        (
            holds(op, false, front) != holds(op, true, back),
            holds(op, true, back),
        )
    } else {
        (
            front == back && holds(op, front, false) != holds(op, back, true),
            holds(op, back, true),
        )
    }
}

pub(super) fn holds(op: Op2, a: bool, b: bool) -> bool {
    match op {
        Op2::Fuse => a || b,
        Op2::Cut => a && !b,
        Op2::Common => a && b,
    }
}

/// Twice... the signed area of a polygon.
pub(super) fn area(p: &[[f64; 2]]) -> f64 {
    let n = p.len();
    (0..n)
        .map(|i| {
            let (a, b) = (p[i], p[(i + 1) % n]);
            a[0] * b[1] - a[1] * b[0]
        })
        .sum::<f64>()
        / 2.0
}

/// Whether a point lies inside a polygon (winding), and its distance from
/// the polygon.
pub(super) fn winding(poly: &[[f64; 2]], x: [f64; 2]) -> (bool, f64) {
    let n = poly.len();
    let mut wind = 0i32;
    let mut dist = f64::INFINITY;
    for i in 0..n {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let len2 = dx * dx + dy * dy;
        let t = if len2 > 0.0 {
            (((x[0] - a[0]) * dx + (x[1] - a[1]) * dy) / len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let (px, py) = (a[0] + t * dx - x[0], a[1] + t * dy - x[1]);
        dist = dist.min((px * px + py * py).sqrt());
        let cross = dx * (x[1] - a[1]) - dy * (x[0] - a[0]);
        if a[1] <= x[1] {
            if b[1] > x[1] && cross > 0.0 {
                wind += 1;
            }
        } else if b[1] <= x[1] && cross < 0.0 {
            wind -= 1;
        }
    }
    (wind != 0, dist)
}

/// Solves `a u + b v + c w = d` (binary64, Cramer).
pub(super) fn solve3(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3], d: &[f64; 3]) -> [f64; 3] {
    let det = |x: &[f64; 3], y: &[f64; 3], z: &[f64; 3]| {
        x[0] * (y[1] * z[2] - y[2] * z[1]) - x[1] * (y[0] * z[2] - y[2] * z[0])
            + x[2] * (y[0] * z[1] - y[1] * z[0])
    };
    let dd = det(a, b, c);
    [det(d, b, c) / dd, det(a, d, c) / dd, det(a, b, d) / dd]
}
