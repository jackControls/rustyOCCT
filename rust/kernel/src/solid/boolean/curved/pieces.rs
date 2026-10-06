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
//! S9e.4b.3c.3: a body may be another Boolean of its primitive and the hull
//! of its other planes (`imported::Form`): the hull less the primitive (a
//! groove, slot or dimple: its curved face's material outside its quadric;
//! the hull the first input), the primitive less the hull of its planes
//! turned over (a bite) or their fuse (a boss), the primitive ending at its
//! caps; two faces on one plane facing one way are one plane of the hull.
//! S9e.4b.3c.3b: a body no form matches is a Boolean tree of its primitive
//! and several hulls (`imported::Tree`), each inner Boolean arranged and
//! assembled in turn and given to the next as its given model (S9e.1's, its
//! own assembly its slots: `Tree::node`), the root's matched to the stored
//! topology. S9e.4b.4b.1: a body of several sphere, cylinder and cone faces
//! is such a tree of its primitives alone (`Tree::Primitive(i)`, a chain of
//! them), each primitive's model at a seam of its own; S9e.4b.4b.2a: its first
//! leaf a prism (S9e.4a's, read off a cap of the body). S9e.4b.3c.1: a stored
//! vertex within the resolution of the model's ring (a rim OCCT split at its
//! sphere's seam; S9e.4b.4b.1: a cylinder's or a cone's meeting with a
//! sphere too) splits it in the first arrangement (`Arr::split_at`), so the
//! match takes the stored arcs;
//! under a partner on its sphere a piece's sphere is split at the second
//! arrangement's seam (`common`'s `first`), its faces then on one surface
//! with the partner's (`graph.rs`).
use super::assemble::{assemble_made, Made};
use super::graph::{arrange_shared, Arr};
use super::model::*;
use super::num::*;
use crate::identity::{Derivation, EntityId, EntityKind, OperationId, OperationKind, Role};
use crate::profile::boolean::{Op2, Operand};
use crate::solid::boolean::polyhedra::Component;
use crate::solid::imported::{Form, Piece, Tree as PieceTree};
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
    // Two faces on one plane facing one way (S9e.4b.3c.3: a groove across
    // a face) are one plane of the hull; facing apart they bound nothing.
    let mut j = 0;
    while j < all.len() {
        let alike = (0..j).find(|&i| {
            is_zero(&cross(&all[i].1, &all[j].1))
                && dot(&sub(&all[j].0, &all[i].0), &all[i].1) == zero()
        });
        match alike {
            Some(i) if sign(&dot(&all[i].1, &all[j].1)) == Ordering::Greater => {
                all.remove(j);
                frames.remove(j);
            }
            Some(_) => {
                return Err(Error::Degenerate(
                    "an imported plane piece's two faces on one plane facing apart",
                ))
            }
            None => j += 1,
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
fn arranged(
    piece: &Piece,
    tolerance: crate::Tolerance,
    first: Option<&R>,
    stored: &[[f64; 3]],
) -> Result<(Arr, Vec<(Component, Made)>)> {
    let p = &piece.primitive;
    if let Some(tree) = &piece.tree {
        // S9e.4b.4b.1: every primitive of a body of several.
        let prims: Vec<&Solid> = std::iter::once(&**p).chain(&piece.others).collect();
        return of_tree(&prims, tree, tolerance, first, stored);
    }
    let (centre, half) = cube(p);
    let hull = hull_model(
        &piece.planes,
        centre,
        half,
        if piece.form == Form::Groove {
            Operand::A
        } else {
            Operand::B
        },
        hull_operation(p),
        tolerance,
    )?;
    of_form(p, &hull, piece.form, first, stored)
}

/// S9e.4b.3c.3: a form's first Boolean's operation; S9e.4b.3c.3b: a tree's
/// root Boolean's.
fn first_boolean(piece: &Piece) -> Op2 {
    match (piece.form, &piece.tree) {
        (Form::Common, _) => Op2::Common,
        (Form::Groove | Form::Bite, _) => Op2::Cut,
        (Form::Boss, _) => Op2::Fuse,
        (Form::Tree, Some(PieceTree::Op(op, _, _))) => *op,
        (Form::Tree, _) => unreachable!("a tree of one Boolean or more"),
    }
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
/// (and a whole sphere's split) tried at several rational points, from
/// `first` where given (S9e.4b.3c.1: the second arrangement's seam, so two
/// pieces of one sphere split it apart and move their splits with its
/// seams).
pub(super) fn common(
    p: &Solid,
    hull: &Prism,
    first: Option<&R>,
    stored: &[[f64; 3]],
) -> Result<(Arr, Vec<(Component, Made)>)> {
    of_form(p, hull, Form::Common, first, stored)
}

/// [`common`] for any form (S9e.4b.3c.3): the hull less the primitive with
/// the hull the first input.
fn of_form(
    p: &Solid,
    hull: &Prism,
    form: Form,
    first: Option<&R>,
    stored: &[[f64; 3]],
) -> Result<(Arr, Vec<(Component, Made)>)> {
    let (op, swapped) = match form {
        Form::Common => (Op2::Common, false),
        Form::Groove => (Op2::Cut, true),
        Form::Bite => (Op2::Cut, false),
        Form::Boss => (Op2::Fuse, false),
        Form::Tree => unreachable!("a tree's piece is arranged by `of_tree`"),
    };
    let r = |(n, d): (i64, i64)| R::new(n.into(), d.into());
    let seams = super::SEAMS;
    let tried: Vec<R> = first
        .cloned()
        .into_iter()
        .chain(seams[..seams.len() - 1].iter().map(|&s| r(s)))
        .collect();
    for (k, seam) in tried.iter().enumerate() {
        let attempt = || -> Result<(Arr, Vec<(Component, Made)>)> {
            let mut arr = if swapped {
                let b = super::model_of(p, Operand::B, seam, false)?;
                arrange_shared([hull.clone(), b])?
            } else {
                let a = super::model_of(p, Operand::A, seam, false)?;
                arrange_shared([a, hull.clone()])?
            };
            arr.split_at(stored, p.resolution().linear())?;
            arr.for_op(op);
            let out = assemble_made(&arr, op)?;
            Ok((arr, out))
        };
        match attempt() {
            // At the second arrangement's seam (S9e.4b.3c.1), a split that
            // fails there is tried at the others in turn: a piece degenerate
            // at every split keeps its own refusal, and a split alike the
            // partner's is that arrangement's seam conflict.
            Err(Error::Degenerate(_)) if first.is_some() && k == 0 => continue,
            Err(Error::ComputationLimit(m)) if m == super::graph::SEAM => continue,
            r => return r,
        }
    }
    Err(Error::Degenerate("a meeting at every seam tried"))
}

/// S9e.4b.3c.3b: what a piece's tree is evaluated with: its primitive, the
/// seam tried, the hulls' cube and the stored vertices that split rims.
struct Tree<'a> {
    p: &'a Solid,
    /// S9e.4b.4b.1: every primitive (the piece's first) and its seam.
    prims: &'a [&'a Solid],
    seams: Vec<&'a R>,
    cube: ([f64; 3], f64),
    stored: &'a [[f64; 3]],
    tolerance: crate::Tolerance,
}

impl Tree<'_> {
    /// A tree's model as an input `operand`: the primitive's model, a hull
    /// leaf model (its entities under an operation of its own, `count` the
    /// hulls and Booleans so far), or an inner Boolean's given model.
    fn node(&self, tree: &PieceTree, operand: Operand, count: &mut u64) -> Result<Prism> {
        *count += 1;
        let own = OperationId(hull_operation(self.p).0 ^ (*count << 40));
        match tree {
            PieceTree::Primitive(i) => {
                super::model_of(self.prims[*i], operand, self.seams[*i], false)
            }
            PieceTree::Hull(planes) => {
                let (centre, half) = self.cube;
                hull_model(planes, centre, half, operand, own, self.tolerance)
            }
            PieceTree::Op(op, a, b) => {
                let (arr, out) = self.boolean(*op, a, b, count)?;
                // One solid, its own assembly naming its given model (under
                // ids of its own operation).
                if out.len() != 1 {
                    return Err(crate::solid::imported::not_its_primitive());
                }
                let parts = out[0].0.parts.clone();
                let topology = inner_topology(parts, own)?;
                let stored = super::given::Stored {
                    topology: &topology,
                    resolution: self.tolerance,
                    frame: self.p.frame,
                    poles: false,
                };
                super::given::built_on(&stored, operand, arr, out, *op, Some(0))
            }
        }
    }

    /// A Boolean of two trees, arranged (the stored vertices splitting its
    /// rims) and assembled.
    fn boolean(
        &self,
        op: Op2,
        a: &PieceTree,
        b: &PieceTree,
        count: &mut u64,
    ) -> Result<(Arr, Vec<(Component, Made)>)> {
        let ma = self.node(a, Operand::A, count)?;
        let mb = self.node(b, Operand::B, count)?;
        let mut arr = arrange_shared([ma, mb])?;
        arr.split_at(self.stored, self.p.resolution().linear())?;
        arr.for_op(op);
        let out = assemble_made(&arr, op)?;
        Ok((arr, out))
    }
}

/// An inner Boolean's assembly as a topology under ids of its own operation
/// (unchecked: it only names that Boolean's given model).
fn inner_topology(
    parts: crate::topology::TopologyParts,
    operation: OperationId,
) -> Result<crate::topology::Topology> {
    use crate::topology::{EdgeId, FaceId, RegionId, Slot, VertexId};
    let topology = crate::topology::Topology::from_parts_unchecked(parts);
    let d = |entity, ordinal: usize| Derivation {
        operation,
        kind: OperationKind::External,
        entity,
        role: Role::External,
        ordinal: ordinal as u32,
        parents: Vec::new(),
    };
    let slots = (0..topology.vertices().len())
        .map(|i| (Slot::Vertex(VertexId(i)), d(EntityKind::Vertex, i)))
        .chain((0..topology.edges().len()).map(|i| (Slot::Edge(EdgeId(i)), d(EntityKind::Edge, i))))
        .chain((0..topology.faces().len()).map(|i| (Slot::Face(FaceId(i)), d(EntityKind::Face, i))))
        .chain(
            (1..topology.regions().len())
                .map(|i| (Slot::Region(RegionId(i)), d(EntityKind::Region, i))),
        )
        .collect();
    topology.renamed(d(EntityKind::Body, 0), slots)
}

/// S9e.4b.3c.3b: the cube a tree's hulls are bounded by: about its
/// primitive, its half side past the stored vertices too.
fn tree_cube(p: &Solid, stored: &[[f64; 3]]) -> ([f64; 3], f64) {
    let (centre, half) = cube(p);
    let reach = stored
        .iter()
        .flat_map(|v| (0..3).map(move |i| (v[i] - centre[i]).abs()))
        .fold(half, f64::max);
    (centre, 2.0 * reach + 0.375)
}

/// S9e.4b.3c.3b: a piece's tree's root Boolean arranged and assembled, its
/// inner Booleans' results their given models, a full circle's seam (and a
/// whole sphere's split) tried at several rational points.
fn of_tree(
    prims: &[&Solid],
    tree: &PieceTree,
    tolerance: crate::Tolerance,
    first: Option<&R>,
    stored: &[[f64; 3]],
) -> Result<(Arr, Vec<(Component, Made)>)> {
    let PieceTree::Op(op, a, b) = tree else {
        return Err(crate::solid::imported::not_its_primitive());
    };
    let r = |(n, d): (i64, i64)| R::new(n.into(), d.into());
    let seams = super::SEAMS;
    let tried: Vec<R> = first
        .cloned()
        .into_iter()
        .chain(seams[..seams.len() - 1].iter().map(|&s| r(s)))
        .collect();
    let p = prims[0];
    let cube = tree_cube(p, stored);
    for k in 0..tried.len() {
        // Each primitive's seam apart from the others' (one surface's seams
        // on one line otherwise; S9e.4b.4b.1).
        let seams = (0..prims.len())
            .map(|i| &tried[(k + i) % tried.len()])
            .collect();
        let ctx = Tree {
            p,
            prims,
            seams,
            cube,
            stored,
            tolerance,
        };
        match ctx.boolean(*op, a, b, &mut 0) {
            Err(Error::Degenerate(_)) if first.is_some() && k == 0 => continue,
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
pub(super) fn model(s: &Solid, op: Operand, first: Option<&R>) -> Result<Prism> {
    let Construction::Imported(i) = &s.construction else {
        unreachable!("an imported piece")
    };
    let crate::solid::imported::Recognized::Piece(piece) = &i.recognized else {
        unreachable!("an imported piece")
    };
    let key = format!("{piece:?}\n{:?}\n{op:?}\n{first:?}", s.topology);
    if let Some(hit) = ARRANGED.with(|k| {
        k.borrow()
            .iter()
            .find(|(known, _)| *known == key)
            .map(|(_, m)| m.clone())
    }) {
        return hit;
    }
    let built = (|| {
        let stored: Vec<[f64; 3]> = s
            .topology
            .vertices()
            .iter()
            .map(|v| v.position.to_array())
            .collect();
        // A tangency in the piece's own arrangement is its curved face's
        // with its own planes (S9e.4b.3c.3: a notch whose wall is tangent to
        // the faces it meets), not between a Boolean's inputs.
        // S9e.4b.4b.1: a body of several primitives' own tangency and
        // mismatch are named for it.
        let several = !piece.others.is_empty();
        let (arr, out) = arranged(piece, s.resolution(), first, &stored).map_err(|e| match e {
            Error::Degenerate(m) if m.contains("tangen") && several => Error::Degenerate(
                "an imported body of several primitives whose faces are tangent along an edge",
            ),
            Error::Degenerate(m) if m.contains("tangen") => Error::Degenerate(
                "an imported plane piece whose curved face is tangent to its plane faces",
            ),
            Error::OutOfDomain(_)
                if several && e == crate::solid::imported::not_its_primitive() =>
            {
                crate::solid::imported::not_their_chain()
            }
            e => e,
        })?;
        let op2 = first_boolean(piece);
        super::given::built(s, op, arr, out, op2, None).map_err(|e| match e {
            Error::ComputationLimit(m) if m.contains("rebuilt differently") && several => {
                crate::solid::imported::not_their_chain()
            }
            Error::ComputationLimit(m) if m.contains("rebuilt differently") => {
                crate::solid::imported::not_its_primitive()
            }
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
            // S9e.4b.4b.1: a body of several primitives' sphere.
            crate::solid::imported::Recognized::Piece(p) => std::iter::once(&*p.primitive)
                .chain(&p.others)
                .find(|x| matches!(x.construction, Construction::Sphere { .. }))
                .unwrap_or(&p.primitive),
            _ => return None,
        },
        _ => s,
    };
    match &s.construction {
        Construction::Sphere { radius, .. } => Some((s.frame.origin(), *radius)),
        _ => None,
    }
}

/// Whether two inputs' spheres are one (S9e.4b.3c.1): a sphere's, cap's or
/// zone's, an imported piece's primitive or a split zone's, equal centres
/// and radii (their exact models' alike); the arrangement then takes their
/// faces as faces on one surface, each piece's sphere split at its seam.
pub(super) fn on_one_sphere(a: &Solid, b: &Solid) -> bool {
    matches!((ball_of(a), ball_of(b)), (Some(x), Some(y)) if x == y)
}

/// Whether a solid is an imported plane piece (S9e.4b.3a).
pub(crate) fn is_piece(s: &Solid) -> bool {
    matches!(&s.construction, Construction::Imported(i) if i.piece())
}

/// An imported piece's model builds and matches its stored topology:
/// checked once on import.
pub(crate) fn check(s: &Solid) -> Result<()> {
    model(s, Operand::A, None).map(|_| ())
}
