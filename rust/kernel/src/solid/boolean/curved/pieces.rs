//! S9e.4b.3a: imported plane pieces of a sphere, a cylinder or a cone
//! (REVIEW_NOTES.md, "S9e.4b.3 refined").
//!
//! An imported body of one sphere, cylinder or cone face and plane faces is
//! its primitive (the whole sphere, or a cylinder or cone on the curved
//! face's stored frame reaching past the body's ends) common the half-spaces
//! of its plane faces' stored planes: the given model (S9e.1) of that first
//! Boolean, its second input the convex hull of the planes (`Hull`, bounded
//! by a cube about the primitive), its arrangement and assembly matched to
//! the stored topology (S9e.2's match), whose ids name the model, so a
//! Boolean's history is over the stored ids directly. Each plane is its
//! stored frame's `o + u x + v y` (the stored axes as rationals, its normal
//! `x * y`), the side of the material by the face's region; the hull's
//! vertices are three planes' common points (rationals), its edges lines.
//! A body the model does not match (not convex in its planes, a stored
//! vertex splitting a ring) and two inputs on one sphere are S9e.4b.3c's.
use super::assemble::{assemble_made, Made};
use super::graph::{arrange_shared, Arr};
use super::model::*;
use super::num::*;
use crate::identity::{Derivation, EntityId, EntityKind, OperationId, OperationKind, Role};
use crate::profile::boolean::{Op2, Operand};
use crate::solid::boolean::polyhedra::Component;
use crate::solid::imported::Piece;
use crate::solid::split::{q, rational_f64, zero};
use crate::solid::{Construction, Solid};
use crate::topology::Surface;
use crate::{Error, Frame3, Point3, Result, Vec3};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// A convex body of half-spaces: each face's plane, a point and its
/// outward normal (the material where `m . (x - p) <= 0`).
#[derive(Debug, Clone)]
pub(super) struct Hull {
    planes: Vec<(V, V)>,
    /// Each model face's plane.
    face_plane: Vec<usize>,
}

impl Hull {
    fn side(&self, k: usize, p: &QV) -> Qd {
        let (p0, m) = &self.planes[k];
        qdot(&qsub(p, &qv(p0)), m)
    }

    /// Where a point lies in the hull, pushed along directions in turn.
    pub(super) fn member(&self, p: &QV, dirs: &[QV]) -> Loc {
        let mut on = false;
        for (k, (_, m)) in self.planes.iter().enumerate() {
            let mut s = self.side(k, p).sign();
            let mut i = 0;
            while s == Ordering::Equal && i < dirs.len() {
                s = qdot(&dirs[i], m).sign();
                i += 1;
            }
            match s {
                Ordering::Greater => return Loc::Out,
                Ordering::Equal => on = true,
                Ordering::Less => {}
            }
        }
        if on {
            Loc::On
        } else {
            Loc::In
        }
    }

    /// Where a point on a face's plane lies in the face: every other
    /// half-space's side.
    pub(super) fn in_face(&self, fi: usize, p: &QV) -> Loc {
        let own = self.face_plane[fi];
        let mut on = false;
        for k in 0..self.planes.len() {
            if k == own {
                continue;
            }
            match self.side(k, p).sign() {
                Ordering::Greater => return Loc::Out,
                Ordering::Equal => on = true,
                Ordering::Less => {}
            }
        }
        if on {
            Loc::On
        } else {
            Loc::In
        }
    }
}

/// A plane's exact data from a frame: its origin and `x * y` (exact), turned
/// to leave the material where `flip`.
fn exact_plane(frame: &Frame3, flip: bool) -> (V, V) {
    let v = |a: Vec3| a.to_array().map(q);
    let p = (frame.origin() - Point3::ORIGIN).to_array().map(q);
    let m = cross(&v(frame.x()), &v(frame.y()));
    (p, if flip { neg(&m) } else { m })
}

/// Three planes' common point, if they meet in one.
fn corner(a: &(V, V), b: &(V, V), c: &(V, V)) -> Option<V> {
    let (ma, mb, mc) = (&a.1, &b.1, &c.1);
    let bc = cross(mb, mc);
    let det = dot(ma, &bc);
    if det == zero() {
        return None;
    }
    let (da, db, dc) = (dot(ma, &a.0), dot(mb, &b.0), dot(mc, &c.0));
    let x = add(
        &scale(&bc, &da),
        &add(&scale(&cross(mc, ma), &db), &scale(&cross(ma, mb), &dc)),
    );
    Some(scale(&x, &(int(1) / det)))
}

fn derived(operation: OperationId, entity: EntityKind, ordinal: usize) -> EntityId {
    Derivation {
        operation,
        kind: OperationKind::External,
        entity,
        role: Role::External,
        ordinal: ordinal as u32,
        parents: Vec::new(),
    }
    .id()
}

/// The pseudo-angle order of two rational directions counter-clockwise
/// from `+u`.
fn angle_order(a: &[R; 2], b: &[R; 2]) -> Ordering {
    let lower = |v: &[R; 2]| {
        let s = sign(&v[1]);
        s == Ordering::Less || (s == Ordering::Equal && sign(&v[0]) == Ordering::Less)
    };
    match (lower(a), lower(b)) {
        (false, true) => Ordering::Less,
        (true, false) => Ordering::Greater,
        _ => sign(&(&b[0] * &a[1] - &a[0] * &b[1])),
    }
}

/// The hull of the planes (and a cube about `centre` of half side `half`)
/// as a model: each plane's polygon of the corners on it, counter-clockwise
/// about its outward normal; each edge between two corners of two faces.
/// `frames` are the planes' stored frames (the faces' parameters).
pub(super) fn hull_model(
    planes: &[(Frame3, bool)],
    centre: [f64; 3],
    half: f64,
    op: Operand,
    operation: OperationId,
    tolerance: crate::Tolerance,
) -> Result<Prism> {
    let exact: Vec<((V, V), Frame3)> = planes
        .iter()
        .map(|(f, s)| (exact_plane(f, *s), *f))
        .collect();
    hull_of(&exact, centre, half, op, operation, tolerance)
}

/// [`hull_model`] of exact planes (a point and the outward normal), each
/// with a frame on it for its face's parameters (S9e.4b.3b: a split's
/// plane in its solid's frame).
pub(super) fn hull_of(
    planes: &[((V, V), Frame3)],
    centre: [f64; 3],
    half: f64,
    op: Operand,
    operation: OperationId,
    tolerance: crate::Tolerance,
) -> Result<Prism> {
    let mut all: Vec<(V, V)> = planes.iter().map(|(p, _)| p.clone()).collect();
    let mut frames: Vec<Frame3> = planes.iter().map(|(_, f)| *f).collect();
    for k in 0..3 {
        for s in [-1.0, 1.0] {
            let mut o = centre;
            o[k] += s * half;
            let mut n = [0.0; 3];
            n[k] = s;
            let x = [n[1].abs(), n[2].abs(), n[0].abs()];
            let origin = Point3::new(o[0], o[1], o[2]);
            let frame = Frame3::new(
                origin,
                Vec3::new(n[0], n[1], n[2]),
                Vec3::new(x[0], x[1], x[2]),
                tolerance,
            )?;
            all.push(exact_plane(&frame, false));
            frames.push(frame);
        }
    }
    // Two planes alike (one plane's two faces) are not a hull's.
    for i in 0..all.len() {
        for j in i + 1..all.len() {
            if is_zero(&cross(&all[i].1, &all[j].1))
                && dot(&sub(&all[j].0, &all[i].0), &all[i].1) == zero()
            {
                return Err(Error::OutOfDomain(
                    "an imported plane piece with two faces on one plane (S9e.4b.3c)",
                ));
            }
        }
    }
    let n = all.len();
    let inside = |x: &V| {
        all.iter()
            .all(|(p, m)| sign(&dot(&sub(x, p), m)) != Ordering::Greater)
    };
    let mut corners: Vec<V> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let Some(x) = corner(&all[i], &all[j], &all[k]) else {
                    continue;
                };
                if inside(&x) && !corners.contains(&x) {
                    corners.push(x);
                }
            }
        }
    }
    let on = |k: usize, x: &V| dot(&sub(x, &all[k].0), &all[k].1) == zero();
    // Each plane's polygon, counter-clockwise about its outward normal.
    let mut polygons: Vec<(usize, Vec<usize>)> = Vec::new();
    for (k, (_, m)) in all.iter().enumerate() {
        let mut pts: Vec<usize> = (0..corners.len()).filter(|&c| on(k, &corners[c])).collect();
        if pts.len() < 3 {
            continue;
        }
        // The coordinate plane of the normal's largest component, cyclic.
        let a = (0..3)
            .max_by(|&i, &j| (&m[i] * &m[i]).cmp(&(&m[j] * &m[j])))
            .expect("three");
        let (b, c) = ((a + 1) % 3, (a + 2) % 3);
        let count = int(pts.len() as i64);
        let mid: V =
            [0, 1, 2].map(|i| pts.iter().fold(zero(), |s, &p| s + &corners[p][i]) / &count);
        let dir = |p: usize| [&corners[p][b] - &mid[b], &corners[p][c] - &mid[c]];
        pts.sort_by(|&x, &y| angle_order(&dir(x), &dir(y)));
        if sign(&m[a]) == Ordering::Less {
            pts.reverse();
        }
        // Three corners on a line: a corner on another face's edge.
        for w in 0..pts.len() {
            let (p0, p1, p2) = (
                &corners[pts[w]],
                &corners[pts[(w + 1) % pts.len()]],
                &corners[pts[(w + 2) % pts.len()]],
            );
            if is_zero(&cross(&sub(p1, p0), &sub(p2, p1))) {
                return Err(Error::Degenerate("a hull's corner on an edge"));
            }
        }
        polygons.push((k, pts));
    }
    // Faces, edges and vertices.
    let face_of: BTreeMap<usize, usize> = polygons
        .iter()
        .enumerate()
        .map(|(fi, (k, _))| (*k, fi))
        .collect();
    let mut faces = Vec::new();
    let mut boxes = Vec::new();
    for (fi, (k, pts)) in polygons.iter().enumerate() {
        let (p, m) = &all[*k];
        faces.push(MFace {
            kind: FaceKind::Facet(*k),
            surf: Surf::Plane {
                p: p.clone(),
                m: m.clone(),
            },
            id: derived(operation, EntityKind::Face, fi),
            stored: Surface::Plane(frames[*k]),
            sense: crate::topology::Orientation::Forward,
        });
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &c in pts {
            for i in 0..3 {
                let x = rational_f64(&corners[c][i]);
                lo[i] = lo[i].min(x);
                hi[i] = hi[i].max(x);
            }
        }
        for i in 0..3 {
            let e = 1e-9 * (1.0 + lo[i].abs().max(hi[i].abs()));
            lo[i] -= e;
            hi[i] += e;
        }
        boxes.push((lo, hi));
    }
    let mut used: Vec<Option<usize>> = vec![None; corners.len()];
    let mut verts: Vec<MVert> = Vec::new();
    let mut vertex = |c: usize, verts: &mut Vec<MVert>| -> usize {
        *used[c].get_or_insert_with(|| {
            verts.push(MVert {
                p: qv(&corners[c]),
                id: Some(derived(operation, EntityKind::Vertex, verts.len())),
            });
            verts.len() - 1
        })
    };
    let mut edges: Vec<MEdge> = Vec::new();
    let mut seen: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (fi, (k, pts)) in polygons.iter().enumerate() {
        for w in 0..pts.len() {
            let (s, t) = (pts[w], pts[(w + 1) % pts.len()]);
            let key = (s.min(t), s.max(t));
            if seen.contains_key(&key) {
                continue;
            }
            // The other face holding both corners.
            let other = polygons
                .iter()
                .find(|(j, q)| j != k && q.contains(&s) && q.contains(&t))
                .map(|(j, _)| face_of[j])
                .ok_or(Error::InvalidTopology("a hull's edge on one face"))?;
            let (a, b) = (vertex(s, &mut verts), vertex(t, &mut verts));
            seen.insert(key, edges.len());
            edges.push(MEdge {
                kind: EdgeKind::Facet(edges.len()),
                curve: Crv::Line {
                    p: qv(&corners[s]),
                    d: sub(&corners[t], &corners[s]),
                },
                arc: None,
                start: a,
                end: b,
                faces: [fi, other],
                id: Some(derived(operation, EntityKind::Edge, edges.len())),
            });
        }
    }
    let mut info = BTreeMap::new();
    for f in &faces {
        info.insert(f.id, (op, Role::External));
    }
    for e in &edges {
        info.insert(e.id.expect("an id"), (op, Role::External));
    }
    for v in &verts {
        info.insert(v.id.expect("an id"), (op, Role::External));
    }
    let region = derived(operation, EntityKind::Region, 1);
    info.insert(region, (op, Role::External));
    let face_plane = polygons.iter().map(|(k, _)| *k).collect();
    Ok(Prism {
        frame: Frame3::xy(),
        tolerance,
        region,
        f: Affine::new(&Frame3::xy())?,
        lo: zero(),
        hi: zero(),
        bounds: Vec::new(),
        faces,
        edges,
        verts,
        info,
        boxes,
        ball: None,
        funnel: None,
        ring: None,
        given: None,
        hull: Some(Box::new(Hull {
            planes: all,
            face_plane,
        })),
    })
}

/// The hull's entities' operation: the primitive's, its bits turned.
pub(super) fn hull_operation(primitive: &Solid) -> OperationId {
    OperationId(!primitive.operation.0 ^ 0x5a5a)
}

/// The primitive common the hull, arranged and assembled, a full circle's
/// seam (and a whole sphere's split) tried at several rational points.
fn arranged(piece: &Piece, tolerance: crate::Tolerance) -> Result<(Arr, Vec<(Component, Made)>)> {
    let p = &piece.primitive;
    let (centre, half) = cube(p);
    let hull = hull_model(
        &piece.planes,
        centre,
        half,
        Operand::B,
        hull_operation(p),
        tolerance,
    )?;
    common(p, &hull)
}

/// The cube about a primitive the hull is bounded by: its centre and half
/// side.
pub(super) fn cube(p: &Solid) -> ([f64; 3], f64) {
    let b = p.bounds;
    let (lo, hi) = (b.min.to_array(), b.max.to_array());
    let centre = [0, 1, 2].map(|i| (lo[i] + hi[i]) / 2.0);
    let reach = (0..3).map(|i| hi[i] - lo[i]).fold(1.0, f64::max);
    (centre, 2.0 * reach)
}

/// A primitive common a hull, arranged and assembled, a full circle's seam
/// (and a whole sphere's split) tried at several rational points.
pub(super) fn common(p: &Solid, hull: &Prism) -> Result<(Arr, Vec<(Component, Made)>)> {
    let seams = super::SEAMS;
    let r = |(n, d): (i64, i64)| R::new(n.into(), d.into());
    for &seam in &seams[..seams.len() - 1] {
        let attempt = || -> Result<(Arr, Vec<(Component, Made)>)> {
            let a = super::model_of(p, Operand::A, &r(seam))?;
            let mut arr = arrange_shared([a, hull.clone()])?;
            arr.for_op(Op2::Common);
            let out = assemble_made(&arr, Op2::Common)?;
            Ok((arr, out))
        };
        match attempt() {
            Err(Error::ComputationLimit(m)) if m == super::graph::SEAM => continue,
            r => return r,
        }
    }
    Err(Error::Degenerate("a meeting at every seam tried"))
}

thread_local! {
    /// The last pieces' arrangements, by their content (the piece's and the
    /// stored topology's `Debug` text).
    static ARRANGED: std::cell::RefCell<Vec<(String, Result<Prism>)>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// An imported plane piece's model (S9e.4b.3a): the given model of its
/// primitive common its hull, matched to its stored topology.
pub(super) fn model(s: &Solid, op: Operand) -> Result<Prism> {
    let Construction::Imported(i) = &s.construction else {
        unreachable!("an imported piece")
    };
    let crate::solid::imported::Recognized::Piece(piece) = &i.recognized else {
        unreachable!("an imported piece")
    };
    let key = format!("{piece:?}\n{:?}\n{op:?}", s.topology);
    if let Some(hit) = ARRANGED.with(|k| {
        k.borrow()
            .iter()
            .find(|(known, _)| *known == key)
            .map(|(_, m)| m.clone())
    }) {
        return hit;
    }
    let built = (|| {
        let (arr, out) = arranged(piece, s.resolution())?;
        super::given::built(s, op, arr, out, Op2::Common, None).map_err(|e| match e {
            Error::ComputationLimit(m) if m.contains("rebuilt differently") => Error::OutOfDomain(
                "an imported plane piece other than its primitive common its planes' half-spaces \
                 (S9e.4b.3c)",
            ),
            e => e,
        })
    })();
    ARRANGED.with(|k| {
        let mut k = k.borrow_mut();
        if k.len() >= 2 {
            k.drain(..1);
        }
        k.push((key, built.clone()));
    });
    built
}

/// A solid's sphere, its centre and radius: a sphere's, cap's or zone's, or
/// an imported piece's primitive's.
fn ball_of(s: &Solid) -> Option<(Point3, f64)> {
    if let Some(ball) = super::splits::ball_of(s) {
        return Some(ball);
    }
    let s = match &s.construction {
        Construction::Imported(i) => match &i.recognized {
            crate::solid::imported::Recognized::Piece(p) => &p.primitive,
            _ => return None,
        },
        _ => s,
    };
    match &s.construction {
        Construction::Sphere { radius, .. } => Some((s.frame.origin(), *radius)),
        _ => None,
    }
}

/// Faces of both inputs on one sphere where either is an imported piece
/// (S9e.4b.3b: or a split zone): S9e.4b.3c's (two pieces of one sphere, the
/// DRAW survey's `so1` to `so7`).
pub(super) fn one_sphere(a: &Solid, b: &Solid) -> Result<()> {
    let split = |s: &Solid| super::splits::ball_of(s).is_some();
    if !(is_piece(a) || is_piece(b) || split(a) || split(b)) {
        return Ok(());
    }
    match (ball_of(a), ball_of(b)) {
        (Some(x), Some(y)) if x == y => Err(Error::OutOfDomain(
            "faces of both inputs on one sphere (S9e.4b.3c)",
        )),
        _ => Ok(()),
    }
}

/// Whether a solid is an imported plane piece (S9e.4b.3a).
pub(crate) fn is_piece(s: &Solid) -> bool {
    matches!(&s.construction, Construction::Imported(i) if i.piece())
}

/// An imported piece's model builds and matches its stored topology:
/// checked once on import.
pub(crate) fn check(s: &Solid) -> Result<()> {
    model(s, Operand::A).map(|_| ())
}
