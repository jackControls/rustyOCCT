//! S9e.1: a Boolean's result given to another Boolean (REVIEW_NOTES.md,
//! S9e's refined decisions).
//!
//! What the result stores is rounded (its vertices, its edges' curves, its
//! faces' surfaces on frames at rounded points), so it is decided on its
//! construction's exact model instead: the arrangement of its inputs'
//! exact models that built it, run again (the seams tried in the same
//! order), whose kept pieces are the body. The model's faces are the
//! inputs' faces holding kept pieces, each on its input's exact surface in
//! its input's frame (`view`), its region where its input's `in_face`
//! holds and the other input's sides at the point are those the first
//! Boolean keeps, its orientation the kept pieces'; its edges are the
//! arrangement's edges on the kept pieces' boundaries (exact lines and
//! conics), its vertices the arrangement's exact vertices; its membership
//! is the first Boolean's set function over both inputs' exact membership.
//! The stored topology only names it: the re-run's assembly gives the
//! result's slots in the stored order, checked against the stored
//! vertices, and each face, edge and vertex of the model carries the ids of
//! the result's entities it lies in.
use super::assemble::{assemble_made, Made};
use super::graph::{classify, holds, Arr, CurveRef};
use super::model::*;
use super::num::*;
use crate::identity::{EntityId, Role};
use crate::profile::boolean::{Op2, Operand};
use crate::solid::boolean::polyhedra::{Component, Polyhedron};
use crate::solid::split::zero;
use crate::solid::{Construction, Solid};
use crate::topology::{Curve3, EdgeId, FaceId, RegionId, Slot, VertexId};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// A given result's model data (S9e.1).
#[derive(Debug, Clone)]
pub(super) struct Given {
    /// The first Boolean's inputs' models, with its seams.
    leaves: [Prism; 2],
    op: Op2,
    /// The inputs' faces on one surface (A's, B's).
    coinc: BTreeSet<(usize, usize)>,
    /// Each model face's input and input face.
    leaf: Vec<(usize, usize)>,
    /// Each model face keeps its input face's orientation.
    behind: Vec<bool>,
    /// Each model face's result faces' ids, the first its own (several
    /// where the result left its input face in parts).
    pub(super) ids: Vec<Vec<EntityId>>,
    /// The result face's id on a model face's side of each of its edges.
    pub(super) sides: BTreeMap<(usize, usize), EntityId>,
    /// Each model edge's conic's carrier: the model face whose cylinder's
    /// angle places it (its conic's parameter).
    pub(super) own: Vec<Option<usize>>,
    /// Each model edge's result edge's stored curve (none inside a result
    /// face).
    pub(super) curves: Vec<Option<Curve3>>,
}

/// Whether S9e.1's model applies to an input: a Boolean's result from the
/// curved arrangement.
pub(super) fn applies(s: &Solid) -> bool {
    match &s.construction {
        Construction::Polyhedron(p) => super::applies(p),
        _ => false,
    }
}

fn out_of_domain(what: &'static str) -> Error {
    Error::OutOfDomain(what)
}

/// The first Boolean's arrangement and solids with their provenance, the
/// seams tried in `build`'s order.
fn rerun(poly: &Polyhedron) -> Result<(Arr, Vec<(Component, Made)>)> {
    let seams = super::SEAMS;
    let r = |(n, d): (i64, i64)| num_rational::BigRational::new(n.into(), d.into());
    for k in 0..seams.len() - 1 {
        let attempt = || -> Result<(Arr, Vec<(Component, Made)>)> {
            let mut arr = super::shared(poly, &r(seams[k]), &r(seams[k + 1]))?;
            arr.for_op(poly.op);
            let out = assemble_made(&arr, poly.op)?;
            Ok((arr, out))
        };
        match attempt() {
            Err(Error::ComputationLimit(m)) if m == super::graph::SEAM => continue,
            r => return r,
        }
    }
    Err(Error::Degenerate("a meeting at every seam tried"))
}

/// The given model of a Boolean's result.
pub(super) fn model(s: &Solid, op: Operand) -> Result<Prism> {
    let Construction::Polyhedron(poly) = &s.construction else {
        return Err(out_of_domain(
            "a solid other than a Boolean's result given to a Boolean of curved faces (S9e.2)",
        ));
    };
    // S9e.1: a result of prisms (planes and cylinders).
    if [&poly.a, &poly.b]
        .iter()
        .any(|x| !matches!(x.construction, Construction::Prism(_)))
    {
        return Err(out_of_domain(
            "a Boolean's result of solids other than prisms given to another Boolean (S9e.3)",
        ));
    }
    let (arr, mut out) = rerun(poly)?;
    if out.len() != 1 || poly.index != 0 {
        return Err(out_of_domain(
            "a Boolean's result of several solids given to another Boolean (S9e.2)",
        ));
    }
    let (component, made) = out.swap_remove(0);
    let t = &s.topology;
    // The re-run's slots are the stored ones.
    let tol = s.resolution().linear();
    let differ = || Error::ComputationLimit("a given result rebuilt differently");
    let p = &component.parts;
    if p.vertices.len() != t.vertices().len()
        || p.edges.len() != t.edges().len()
        || p.faces.len() != t.faces().len()
    {
        return Err(differ());
    }
    for (a, b) in p.vertices.iter().zip(t.vertices()) {
        let (a, b) = (a.position.to_array(), b.position.to_array());
        if (0..3).any(|k| (a[k] - b[k]).abs() > tol) {
            return Err(differ());
        }
    }
    let id = |slot: Slot| {
        t.id_of(slot)
            .ok_or(Error::InvalidTopology("an unnamed slot"))
    };
    // Every edge on a kept piece a line or a conic (S9c.1's).
    for parts in &made.edges {
        for &(g, _) in parts {
            if !matches!(arr.edges[g].crv, Crv::Line { .. } | Crv::Conic { .. }) {
                return Err(out_of_domain(
                    "a Boolean's result with procedural edges given to another Boolean (S9e.3)",
                ));
            }
        }
    }
    let leaves = arr.models.clone();
    // Faces: the input faces holding kept pieces.
    let mut face_of: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut faces: Vec<MFace> = Vec::new();
    let mut leaf: Vec<(usize, usize)> = Vec::new();
    let mut behind: Vec<bool> = Vec::new();
    let mut ids: Vec<Vec<EntityId>> = Vec::new();
    let mut piece_face: BTreeMap<usize, (usize, EntityId)> = BTreeMap::new();
    for (k, pieces) in made.faces.iter().enumerate() {
        let rid = id(Slot::Face(FaceId(k)))?;
        let stored = &t.faces()[k];
        for &pi in pieces {
            let piece = &arr.pieces[pi];
            let key = (piece.op, piece.face);
            let mf = match face_of.get(&key) {
                Some(&mf) => mf,
                None => {
                    let lf = &leaves[piece.op].faces[piece.face];
                    let surf = match (&lf.surf, piece.behind) {
                        (s, true) => s.clone(),
                        (Surf::Plane { p, m }, false) => Surf::Plane {
                            p: p.clone(),
                            m: neg(m),
                        },
                        (Surf::Cyl { c, r, inside }, false) => Surf::Cyl {
                            c: c.clone(),
                            r: r.clone(),
                            inside: !inside,
                        },
                        _ => {
                            return Err(out_of_domain(
                                "a Boolean's result of solids other than prisms given to another Boolean (S9e.3)",
                            ))
                        }
                    };
                    faces.push(MFace {
                        kind: lf.kind,
                        surf,
                        id: rid,
                        stored: stored.surface.clone(),
                        sense: stored.sense,
                    });
                    leaf.push(key);
                    behind.push(piece.behind);
                    ids.push(Vec::new());
                    face_of.insert(key, faces.len() - 1);
                    faces.len() - 1
                }
            };
            if !ids[mf].contains(&rid) {
                ids[mf].push(rid);
            }
            piece_face.insert(pi, (mf, rid));
        }
    }
    // Edges: the arrangement's edges on the kept pieces' boundaries, each
    // run the way its faces' loops about their own (the model's) normals
    // leave its left face.
    let mut edge_rid: BTreeMap<usize, (EntityId, usize)> = BTreeMap::new();
    for (i, parts) in made.edges.iter().enumerate() {
        let rid = id(Slot::Edge(EdgeId(i)))?;
        for &(g, _) in parts {
            edge_rid.insert(g, (rid, i));
        }
    }
    let mut vertex_rid: BTreeMap<usize, EntityId> = BTreeMap::new();
    for (j, v) in made.vertices.iter().enumerate() {
        if let Some(v) = v {
            vertex_rid.insert(*v, id(Slot::Vertex(VertexId(j)))?);
        }
    }
    // (left, right) model faces of each arrangement edge, and the result
    // face on each side.
    let mut sides_of: BTreeMap<usize, [Option<(usize, EntityId)>; 2]> = BTreeMap::new();
    for (&pi, &(mf, rid)) in &piece_face {
        let piece = &arr.pieces[pi];
        for lp in &piece.loops {
            for &(g, d) in lp {
                let forward = d == piece.behind;
                let slot = &mut sides_of.entry(g).or_insert([None, None])[usize::from(!forward)];
                if slot.is_some() {
                    return Err(Error::InvalidTopology(
                        "a given result's edge used twice one way",
                    ));
                }
                *slot = Some((mf, rid));
            }
        }
    }
    let mut verts: Vec<MVert> = Vec::new();
    let mut vert_of: BTreeMap<usize, usize> = BTreeMap::new();
    let mut edges: Vec<MEdge> = Vec::new();
    let mut own: Vec<Option<usize>> = Vec::new();
    let mut curves: Vec<Option<Curve3>> = Vec::new();
    let mut sides: BTreeMap<(usize, usize), EntityId> = BTreeMap::new();
    for (&g, fs) in &sides_of {
        let [Some((left, lrid)), Some((right, rrid))] = *fs else {
            return Err(Error::InvalidTopology("a given result's open edge"));
        };
        let ge = &arr.edges[g];
        let mut vertex = |v: usize| -> usize {
            *vert_of.entry(v).or_insert_with(|| {
                verts.push(MVert {
                    p: arr.vx[v].p.clone(),
                    id: vertex_rid.get(&v).copied(),
                });
                verts.len() - 1
            })
        };
        let (start, end) = (vertex(ge.ends[0]), vertex(ge.ends[1]));
        let (curve, arc) = match (&ge.crv, &ge.pos) {
            // A line runs from its start to its end along its direction.
            (Crv::Line { p, d }, _) => (
                Crv::Line {
                    p: p.clone(),
                    d: if ge.with { d.clone() } else { neg(d) },
                },
                None,
            ),
            (Crv::Conic { .. }, [super::meet::Pos::Ang(a), super::meet::Pos::Ang(b)]) => {
                (ge.crv.clone(), Some((a.clone(), b.clone(), ge.with)))
            }
            _ => {
                return Err(out_of_domain(
                    "a Boolean's result with procedural edges given to another Boolean (S9e.3)",
                ))
            }
        };
        let ei = edges.len();
        // The conic's carrier: the cylinder whose angle is its parameter.
        let carrier = match (&ge.crv, ge.curve) {
            (Crv::Conic { .. }, CurveRef::Section(si, _)) => {
                let sec = &arr.secs[si];
                if matches!(leaves[0].faces[sec.fa].surf, Surf::Cyl { .. }) {
                    Some((0, sec.fa))
                } else {
                    Some((1, sec.fb))
                }
            }
            (Crv::Conic { .. }, CurveRef::Edge(o, e)) => leaves[o].edges[e]
                .faces
                .iter()
                .copied()
                .find(|&f| matches!(leaves[o].faces[f].surf, Surf::Cyl { .. }))
                .map(|f| (o, f)),
            _ => None,
        };
        own.push(carrier.and_then(|k| face_of.get(&k).copied()));
        let rid = edge_rid.get(&g);
        curves.push(rid.map(|&(_, i)| t.edges()[i].curve.clone()));
        edges.push(MEdge {
            kind: EdgeKind::Given(ei),
            curve,
            arc,
            start,
            end,
            faces: [left, right],
            id: rid.map(|x| x.0),
        });
        sides.insert((ei, left), lrid);
        sides.insert((ei, right), rrid);
    }
    let mut info = BTreeMap::new();
    for (eid, _) in t.ids() {
        let role = t.derivation(eid).map_or(Role::External, |d| d.role);
        info.insert(eid, (op, role));
    }
    let boxes = leaf.iter().map(|&(o, f)| leaves[o].boxes[f]).collect();
    let f = leaves[0].f.clone();
    Ok(Prism {
        frame: s.frame,
        tolerance: s.resolution(),
        region: id(Slot::Region(RegionId(1)))?,
        f,
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
        given: Some(Box::new(Given {
            leaves,
            op: poly.op,
            coinc: arr.coinc.clone(),
            leaf,
            behind,
            ids,
            sides,
            own,
            curves,
        })),
    })
}

impl Given {
    /// A model face's input model and input face.
    pub(super) fn view(&self, fi: usize) -> (&Prism, usize) {
        let (o, f) = self.leaf[fi];
        (&self.leaves[o], f)
    }

    /// The outward normal of a model face: its input face's, reversed
    /// where the result reverses it (a cut's tool's faces).
    pub(super) fn normal_at(&self, fi: usize, p: &QV) -> QV {
        let (m, f) = self.view(fi);
        let n = m.normal_at(f, p);
        if self.behind[fi] {
            n
        } else {
            n.map(|x| x.neg())
        }
    }

    /// Where a point lies in the result: the first Boolean's set function
    /// of both inputs' membership, pushed alike; a point on either's
    /// boundary after every push is `On` unless the other decides it.
    pub(super) fn member(&self, p: &QV, dirs: &[QV]) -> Loc {
        let a = self.leaves[0].member(p, dirs);
        let b = self.leaves[1].member(p, dirs);
        combine(self.op, a, b)
    }

    /// Where a point on a model face's surface lies in its region: in its
    /// input face, where the other input's sides at the point (pushed off
    /// the face both ways where it lies on a face of the other on the same
    /// surface) are those the first Boolean keeps; on the other's boundary
    /// elsewhere (the first arrangement's sections), `On`.
    pub(super) fn in_face(&self, fi: usize, p: &QV) -> Loc {
        let (o, lf) = self.leaf[fi];
        let own = self.leaves[o].in_face(lf, p);
        if own == Loc::Out {
            return Loc::Out;
        }
        let other = &self.leaves[1 - o];
        let kept = |front: bool, back: bool| {
            if classify(o, (front, back), self.op).0 {
                own
            } else {
                Loc::Out
            }
        };
        match other.member(p, &[]) {
            Loc::In => kept(true, true),
            Loc::Out => kept(false, false),
            Loc::On => {
                for g in 0..other.faces.len() {
                    let key = if o == 0 { (lf, g) } else { (g, lf) };
                    if !self.coinc.contains(&key) {
                        continue;
                    }
                    match other.in_face(g, p) {
                        Loc::In => {
                            let n = self.leaves[o].normal_at(lf, p);
                            let back_dir = n.clone().map(|x| x.neg());
                            let front = other.member(p, &[n]);
                            let back = other.member(p, &[back_dir]);
                            return match (front, back) {
                                (Loc::On, _) | (_, Loc::On) => Loc::On,
                                (f, b) => kept(f == Loc::In, b == Loc::In),
                            };
                        }
                        Loc::On => return Loc::On,
                        Loc::Out => {}
                    }
                }
                Loc::On
            }
        }
    }
}

/// A set function of two memberships, `On` standing for either value: the
/// value when every choice gives it, else `On`.
fn combine(op: Op2, a: Loc, b: Loc) -> Loc {
    let values = |l: Loc| match l {
        Loc::In => vec![true],
        Loc::Out => vec![false],
        Loc::On => vec![false, true],
    };
    let mut seen = BTreeSet::new();
    for x in values(a) {
        for y in values(b) {
            seen.insert(holds(op, x, y));
        }
    }
    match (seen.contains(&false), seen.contains(&true)) {
        (false, true) => Loc::In,
        (true, false) => Loc::Out,
        _ => Loc::On,
    }
}
