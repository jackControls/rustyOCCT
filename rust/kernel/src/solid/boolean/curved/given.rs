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
//!
//! S9e.2: a stack (its two profiles' prisms on its frame) and an S9b.1
//! result of line prisms are decided the same way on the curved
//! arrangement of their construction, their stored slots (another
//! assembly's) matched geometrically to the re-run's (`matched.rs`); a
//! result of several solids is given as its whole construction, every
//! solid's faces and edges in the model, and the second arrangement keeps
//! the given solid's (`matched::keep_solid`).
//!
//! S9e.3a: the construction's inputs may be spheres, cones, tori and given
//! results themselves (`chain.rs`): every face keeps its leaf's surface (a
//! reversed one through its orientation alone), every edge its first
//! arrangement's curve and places, and the views go down to the primitive
//! model.
use super::assemble::{assemble_made, Made};
use super::graph::{classify, holds, Arr, CurveRef};
use super::matched::{matched, Match};
use super::meet::Pos;
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
    /// where the result left its input face in parts; none for a face of
    /// another solid of the result, S9e.2).
    pub(super) ids: Vec<Vec<EntityId>>,
    /// The result face's id on a model face's side of each of its edges.
    pub(super) sides: BTreeMap<(usize, usize), EntityId>,
    /// Each model edge's conic's carrier: the model face whose cylinder's
    /// angle places it (its conic's parameter).
    pub(super) own: Vec<Option<usize>>,
    /// Each model edge's result edge's stored curve (none inside a result
    /// face).
    pub(super) curves: Vec<Option<Curve3>>,
    /// Each model edge's stored curve turning against its conic's
    /// parameter (a matched stored frame's normal against the conic's, S9e.2).
    pub(super) flip: Vec<bool>,
    /// Each model edge's solid among the construction's (S9e.2), and the
    /// given solid's, where the construction has several.
    solids: Vec<usize>,
    solid: Option<usize>,
    /// Each model edge's places at its start and end and whether it runs
    /// with its curve's parameter (S9e.3a: a line's none, its direction
    /// runs from its start).
    pub(super) places: Vec<Option<([Pos; 2], bool)>>,
}

/// Whether a solid has an arc (a prism) or a face other than a plane, or is
/// a result of the curved arrangement (S9e.2's condition on a stack or an
/// S9b.1 result and its partner).
fn curved(s: &Solid) -> bool {
    super::applies_arcs(s)
        || matches!(&s.construction, Construction::Polyhedron(p) if super::applies(p))
        || s.topology
            .faces()
            .iter()
            .any(|f| !matches!(f.surface, crate::topology::Surface::Plane(_)))
}

/// Whether a given model applies to an input `x` with the partner `y`: a
/// Boolean's result from the curved arrangement (S9e.1, of several solids
/// S9e.2), or a stack or polyhedral result where either has an arc or a
/// curved face (S9e.2; else S9b.2's stored model decides it).
pub(super) fn applies(x: &Solid, y: &Solid) -> bool {
    match &x.construction {
        Construction::Polyhedron(p) if super::applies(p) => true,
        Construction::Polyhedron(_) | Construction::Stack(_) => curved(x) || curved(y),
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

/// The construction a given solid's model re-runs: the first Boolean's two
/// prisms and operation (a stack's profiles as prisms on its frame over
/// their heights, S9e.2), and whether its stored slots are the re-run's
/// (S9e.1: the curved arrangement built it) rather than matched.
fn construction(s: &Solid) -> Result<(Polyhedron, bool)> {
    match &s.construction {
        Construction::Polyhedron(poly) => {
            // S9e.1 and S9e.2: a result of prisms (planes and cylinders);
            // S9e.3a: of spheres, cones, tori and given results too, to the
            // depth limit.
            if [&poly.a, &poly.b].iter().any(|x| {
                matches!(
                    x.construction,
                    Construction::Clipped(_) | Construction::Half(_)
                )
            }) {
                return Err(out_of_domain(
                    "a Boolean's result of a plane's piece given to another Boolean (S9e.4)",
                ));
            }
            super::chain::check_depth(s)?;
            Ok((poly.as_ref().clone(), super::applies(poly)))
        }
        Construction::Stack(st) => {
            // The stack's two prisms, entity ids of their own (the model's
            // names are the stored result's).
            let prism = |op: u64, profile: &crate::Profile, h: [f64; 2]| {
                Solid::build(
                    crate::identity::OperationId(u64::MAX - op),
                    profile.clone(),
                    s.frame,
                    h[0],
                    h[1],
                )
                .map(Box::new)
            };
            Ok((
                Polyhedron {
                    a: prism(2, &st.a, st.ha)?,
                    b: prism(1, &st.b, st.hb)?,
                    op: st.op,
                    index: st.index,
                },
                false,
            ))
        }
        _ => Err(out_of_domain(
            "a solid other than a Boolean's result given to a Boolean of curved faces (S9e.4)",
        )),
    }
}

/// The given model of a Boolean's result.
#[allow(clippy::too_many_lines)]
pub(super) fn model(s: &Solid, op: Operand) -> Result<Prism> {
    let (poly, direct) = construction(s)?;
    let (arr, out) = rerun(&poly)?;
    let t = &s.topology;
    let tol = s.resolution().linear();
    let differ = || Error::ComputationLimit("a given result rebuilt differently");
    // The given solid among the re-run's and its slots' stored ones.
    let (k, slots) = if direct {
        // S9e.1: the re-run's slots are the stored ones.
        let k = poly.index;
        let p = &out.get(k).ok_or_else(differ)?.0.parts;
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
        let identity = |n: usize| (0..n).collect::<Vec<_>>();
        (
            k,
            Match {
                faces: identity(p.faces.len()),
                edges: identity(p.edges.len()),
                vertices: identity(p.vertices.len()),
            },
        )
    } else {
        // S9e.2: another assembly's slots, matched geometrically; the one
        // re-run solid that matches.
        let mut found = None;
        for (c, (component, _)) in out.iter().enumerate() {
            if let Some(m) = matched(&component.parts, t, tol) {
                if found.is_some() {
                    return Err(differ());
                }
                found = Some((c, m));
            }
        }
        found.ok_or_else(differ)?
    };
    let id = |slot: Slot| {
        t.id_of(slot)
            .ok_or(Error::InvalidTopology("an unnamed slot"))
    };
    let leaves = arr.models.clone();
    // The given solid's region names the model's faces of other solids
    // (their pieces are dropped from the second arrangement: `keep_solid`).
    let region = id(Slot::Region(RegionId(1)))?;
    // Faces: the input faces holding kept pieces, the given solid's first.
    let order: Vec<usize> = std::iter::once(k)
        .chain((0..out.len()).filter(|&c| c != k))
        .collect();
    let mut face_of: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut faces: Vec<MFace> = Vec::new();
    let mut leaf: Vec<(usize, usize)> = Vec::new();
    let mut behind: Vec<bool> = Vec::new();
    let mut ids: Vec<Vec<EntityId>> = Vec::new();
    // Each kept piece's model face, result face (the given solid's) and
    // solid.
    let mut piece_face: BTreeMap<usize, (usize, Option<EntityId>, usize)> = BTreeMap::new();
    for &c in &order {
        let (component, made) = &out[c];
        for (fk, pieces) in made.faces.iter().enumerate() {
            let (rid, stored) = if c == k {
                let j = slots.faces[fk];
                (Some(id(Slot::Face(FaceId(j)))?), &t.faces()[j])
            } else {
                (None, &component.parts.faces[fk])
            };
            for &pi in pieces {
                let piece = &arr.pieces[pi];
                let key = (piece.op, piece.face);
                let mf = match face_of.get(&key) {
                    Some(&mf) => mf,
                    None => {
                        let lf = &leaves[piece.op].faces[piece.face];
                        // A plane's or a cylinder's side reversed with the
                        // face; a sphere's, cone's or torus's surface as its
                        // leaf's, the face reversed through `behind` alone
                        // (S9e.3a).
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
                            (s, false) => s.clone(),
                        };
                        faces.push(MFace {
                            kind: lf.kind,
                            surf,
                            id: rid.unwrap_or(region),
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
                if let Some(rid) = rid {
                    if !ids[mf].contains(&rid) {
                        ids[mf].push(rid);
                    }
                }
                piece_face.insert(pi, (mf, rid, c));
            }
        }
    }
    // Edges: the arrangement's edges on the kept pieces' boundaries, each
    // run the way its faces' loops about their own (the model's) normals
    // leave its left face.
    let (_, made) = &out[k];
    let mut edge_rid: BTreeMap<usize, (EntityId, usize)> = BTreeMap::new();
    for (i, parts) in made.edges.iter().enumerate() {
        let j = slots.edges[i];
        let rid = id(Slot::Edge(EdgeId(j)))?;
        for &(g, _) in parts {
            edge_rid.insert(g, (rid, j));
        }
    }
    let mut vertex_rid: BTreeMap<usize, EntityId> = BTreeMap::new();
    for (j, v) in made.vertices.iter().enumerate() {
        if let Some(v) = v {
            vertex_rid.insert(*v, id(Slot::Vertex(VertexId(slots.vertices[j])))?);
        }
    }
    // (left, right) model faces of each arrangement edge, the result face
    // on each side and the edge's solid.
    type Side = (usize, Option<EntityId>);
    let mut sides_of: BTreeMap<usize, ([Option<Side>; 2], usize)> = BTreeMap::new();
    for (&pi, &(mf, rid, c)) in &piece_face {
        let piece = &arr.pieces[pi];
        for lp in &piece.loops {
            for &(g, d) in lp {
                let forward = d == piece.behind;
                let entry = sides_of.entry(g).or_insert(([None, None], c));
                if entry.1 != c {
                    return Err(Error::InvalidTopology(
                        "a given result's edge on two of its solids",
                    ));
                }
                let slot = &mut entry.0[usize::from(!forward)];
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
    let mut flip: Vec<bool> = Vec::new();
    let mut solids: Vec<usize> = Vec::new();
    let mut places: Vec<Option<([Pos; 2], bool)>> = Vec::new();
    let mut sides: BTreeMap<(usize, usize), EntityId> = BTreeMap::new();
    for (&g, &(fs, c)) in &sides_of {
        let [Some((left, lrid)), Some((right, rrid))] = fs else {
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
        // A cone's section normal to its axis is the circle it is (S9e.3a).
        let crv = match &ge.crv {
            Crv::Cone(c) => c.as_conic().unwrap_or_else(|| ge.crv.clone()),
            x => x.clone(),
        };
        let (curve, arc) = match (&crv, &ge.pos) {
            // A line runs from its start to its end along its direction.
            (Crv::Line { p, d }, _) => (
                Crv::Line {
                    p: p.clone(),
                    d: if ge.with { d.clone() } else { neg(d) },
                },
                None,
            ),
            // A conic, a sphere's circle, a cone's or a torus's section, a
            // meeting placed by an angle (S9e.3a: every curve the first
            // arrangement makes).
            (_, [Pos::Ang(a), Pos::Ang(b)]) => (crv.clone(), Some((a.clone(), b.clone(), ge.with))),
            // A meeting placed by a height (`Rise`).
            (_, [Pos::T(_), Pos::T(_)]) => (crv.clone(), None),
            _ => {
                return Err(Error::InvalidTopology(
                    "a given result's edge of mixed places",
                ))
            }
        };
        places.push(match &crv {
            Crv::Line { .. } => None,
            _ => Some((ge.pos.clone(), ge.with)),
        });
        let ei = edges.len();
        // The conic's carrier: the cylinder whose angle is its parameter.
        let carrier = match (&crv, ge.curve) {
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
        let stored = rid.map(|&(_, j)| t.edges()[j].curve.clone());
        // A stored circle or ellipse turning against the conic's
        // parameter (its frame's normal against the conic's `a x b`).
        flip.push(match (&crv, &stored) {
            (
                Crv::Conic { a, b, .. },
                Some(
                    Curve3::Circle { frame, .. }
                    | Curve3::CircularArc { frame, .. }
                    | Curve3::EllipseArc { frame, .. },
                ),
            ) => {
                let n = cross(a, b).map(|x| crate::solid::split::rational_f64(&x));
                let m = frame.normal();
                n[0] * m.x + n[1] * m.y + n[2] * m.z < 0.0
            }
            _ => false,
        });
        curves.push(stored);
        solids.push(c);
        edges.push(MEdge {
            kind: EdgeKind::Given(ei),
            curve,
            arc,
            start,
            end,
            faces: [left, right],
            id: rid.map(|x| x.0),
        });
        if let Some(lrid) = lrid {
            sides.insert((ei, left), lrid);
        }
        if let Some(rrid) = rrid {
            sides.insert((ei, right), rrid);
        }
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
        region,
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
            flip,
            solids,
            solid: (out.len() > 1).then_some(k),
            places,
        })),
    })
}

impl Given {
    /// Each model edge's solid and the given solid's, where the given
    /// result is one solid of several (S9e.2).
    pub(super) fn several(&self) -> Option<(&[usize], usize)> {
        self.solid.map(|k| (self.solids.as_slice(), k))
    }

    /// A model face's primitive model and its face: its input's, down
    /// through every given level (S9e.3a).
    pub(super) fn view(&self, fi: usize) -> (&Prism, usize) {
        let (o, f) = self.leaf[fi];
        self.leaves[o].view(f)
    }

    /// Whether a model face runs against its primitive model's face: its
    /// orientation composed through every given level (S9e.3a).
    pub(super) fn reversed(&self, fi: usize) -> bool {
        let (o, f) = self.leaf[fi];
        !self.behind[fi] ^ self.leaves[o].reversed(f)
    }

    /// The outward normal of a model face: its input face's (its own
    /// given model's, S9e.3a), reversed where the result reverses it (a
    /// cut's tool's faces).
    pub(super) fn normal_at(&self, fi: usize, p: &QV) -> QV {
        let (o, f) = self.leaf[fi];
        let n = self.leaves[o].normal_at(f, p);
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
