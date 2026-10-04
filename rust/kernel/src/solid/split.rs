//! Splitting a solid by a plane (S8 of `REVIEW_NOTES.md`).
//!
//! S8a.1: a prism split by a plane normal to its axis is M3's height split;
//! by a plane parallel to its axis, each piece is the prism of one piece of
//! the profile's section by the plane's line, its ids renamed from the input
//! by provenance. Which side of the line every profile vertex lies on, where
//! a line segment crosses it, whether a circle meets it (twice, tangentially
//! or not) are exact decisions on the stored data; an arc's crossings are
//! decided within its range by certified intervals. New points are rounded
//! to binary64; the pieces' prisms validate as every prism does.
use super::attrs::{debug_check_attributes, propagate};
use super::{replayable, Construction, Context, Solid};
use crate::certified::{Interval as I, Real};
use crate::history::{History, Relation};
use crate::identity::{
    Derivation, EntityId, EntityKind, OperationKind, Parent, ProfileElement, Role,
};
use crate::profile::{BoundaryKind, Location, Segment};
use crate::{Boundary, Error, Frame3, Point2, Profile, Result, Tolerance};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

pub(crate) mod conic;
mod meridian;
mod oblique;
mod revolved;
mod spiric;
pub(crate) mod spline;
mod torus;
pub(super) use meridian::{Half, Primitive};
pub(crate) use oblique::edge_bounds;
pub(super) use oblique::Clipped;

pub(crate) fn zero() -> R {
    R::from_integer(0.into())
}

pub(crate) fn q(x: f64) -> R {
    R::from_float(x).expect("finite")
}

/// Which side of the plane a solid or piece lies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Side {
    /// Against the plane's normal.
    Below,
    /// Along the plane's normal.
    Above,
}

impl Solid {
    /// Split by the plane through `plane.origin()` normal to `plane.normal()`
    /// (S8). The pieces come with their sides, those below the plane (against
    /// its normal) first; a plane missing the solid, touching it or lying in
    /// one of its faces returns the solid itself, unchanged.
    ///
    /// S8a.1 splits prisms by planes normal to their axis (M3's height split,
    /// its history `HeightSplit`) or parallel to it (history `PlaneSplit`:
    /// walls, cap edges, caps, vertical edges and the region crossed are
    /// `Split`, one child per part; a cut wall is `Generated` from the caps
    /// and the walls it meets, cut edges from the face they cut, cut vertices
    /// and vertical edges from the edge or wall they cut). Other planes and
    /// solids are `OutOfDomain` until S8a.2 and S8c.
    pub fn split_by_plane(
        &self,
        operation: crate::identity::OperationId,
        plane: Frame3,
    ) -> Result<(Vec<(Side, Solid)>, History)> {
        self.split_by_plane_in(&Context::new(operation), plane)
    }

    /// [`Solid::split_by_plane`] in an operation context.
    pub fn split_by_plane_in(
        &self,
        context: &Context,
        plane: Frame3,
    ) -> Result<(Vec<(Side, Solid)>, History)> {
        replayable(context.level)?;
        // The plane in the solid's frame coordinates: a u + b v + c w + d.
        let dot = |u: [f64; 3], v: [f64; 3]| -> R { (0..3).map(|i| q(u[i]) * q(v[i])).sum() };
        let m = plane.normal().to_array();
        let (x, y, n) = (
            self.frame.x().to_array(),
            self.frame.y().to_array(),
            self.frame.normal().to_array(),
        );
        let (o, p0) = (self.frame.origin().to_array(), plane.origin().to_array());
        let (a, b, c) = (dot(m, x), dot(m, y), dot(m, n));
        let d: R = (0..3).map(|i| q(m[i]) * (q(o[i]) - q(p0[i]))).sum();
        let profile = match &self.construction {
            Construction::Prism(profile) => (**profile).clone(),
            Construction::Cone { .. } | Construction::Sphere { .. } => {
                return Ok(self
                    .split_revolved(context, &plane, [a, b, c, d])?
                    .expect("a cone or sphere"))
            }
            Construction::Torus { .. } => {
                return Ok(self
                    .split_torus(context, &plane, [a, b, c, d])?
                    .expect("a torus"))
            }
            _ => {
                return Err(Error::OutOfDomain(
                    "split_by_plane: a split piece split again (S8b on)",
                ))
            }
        };
        let (low, high) = (q(self.start.min(self.end)), q(self.start.max(self.end)));
        if on_one_side(&profile, [&a, &b, &c, &d], &low, &high)? {
            return Ok(self.unchanged(context));
        }
        if (a == zero()) && (b == zero()) {
            // Normal to the axis: the cut at w = -d / c.
            let h = -d / &c;
            if h <= low || h >= high {
                return Ok(self.unchanged(context));
            }
            let height = rational_f64(&h);
            if !(self.start.min(self.end) < height && height < self.start.max(self.end)) {
                return Err(Error::Degenerate("a split within binary64 of a cap"));
            }
            let ([lower, upper], history) = self.split_at_height_in(context, height)?;
            // The lower piece has w < h: g = c (w - h) < 0 there when c > 0.
            let (sl, su) = if c > zero() {
                (Side::Below, Side::Above)
            } else {
                (Side::Above, Side::Below)
            };
            let mut pieces = vec![(sl, lower), (su, upper)];
            pieces.sort_by_key(|p| p.0);
            return Ok((pieces, history));
        }
        if c != zero() {
            return self.split_oblique(context, &profile, [a, b, c, d]);
        }
        let section = Section::new(&profile, [a, b, d])?;
        if section
            .pieces
            .iter()
            .all(|p| p.side == section.pieces[0].side)
        {
            return Ok(self.unchanged(context));
        }
        self.rebuild(context, &profile, section)
    }

    /// The solid itself, every entity `Unchanged`.
    fn unchanged(&self, context: &Context) -> (Vec<(Side, Solid)>, History) {
        let relations = self
            .topology
            .ids()
            .map(|(id, _)| Relation::Unchanged { id })
            .collect();
        let body = self.topology.body_id();
        let history = History::new(
            context.operation,
            OperationKind::PlaneSplit,
            vec![body],
            vec![body],
            relations,
            Vec::new(),
        )
        .at_level(context.level);
        (vec![(Side::Below, self.clone())], history)
    }

    /// The prisms of the section's pieces, renamed from this prism.
    fn rebuild(
        &self,
        context: &Context,
        profile: &Profile,
        section: Section,
    ) -> Result<(Vec<(Side, Solid)>, History)> {
        let operation = context.operation;
        let parent = &self.topology;
        // The input's entities by what they derive from: (entity, role,
        // profile element) of the input's profile.
        let label_of = |b: usize, e: ProfileElement| -> Parent {
            let boundary = profile.boundaries().nth(b).expect("a boundary");
            match boundary.labels() {
                Some(l) => Parent::Label(match e {
                    ProfileElement::Boundary => l.boundary,
                    ProfileElement::Segment(j) => l.segments[j as usize],
                    ProfileElement::Vertex(j) => l.vertices[j as usize],
                }),
                None => Parent::Profile {
                    boundary: b as u32,
                    element: e,
                },
            }
        };
        let mut by_origin: BTreeMap<(EntityKind, Role, Vec<Parent>), EntityId> = BTreeMap::new();
        for (id, _) in parent.ids() {
            let d = parent.derivation(id).expect("a derivation");
            by_origin.insert((d.entity, d.role, d.parents.clone()), id);
        }
        let original = |entity: EntityKind, role: Role, b: usize, e: ProfileElement| {
            by_origin
                .get(&(entity, role, vec![label_of(b, e)]))
                .copied()
                .ok_or(Error::InvalidTopology(
                    "a prism entity without its profile element",
                ))
        };
        let caps: Vec<EntityId> = [Role::StartCap, Role::EndCap]
            .iter()
            .map(|role| {
                parent
                    .ids()
                    .find(|(id, _)| parent.derivation(*id).map(|d| d.role) == Some(*role))
                    .map(|(id, _)| id)
                    .ok_or(Error::InvalidTopology("a prism without its caps"))
            })
            .collect::<Result<_>>()?;
        let region = parent
            .ids()
            .find(|(id, _)| parent.derivation(*id).map(|d| d.entity) == Some(EntityKind::Region))
            .map(|(id, _)| id)
            .ok_or(Error::InvalidTopology("a prism without its region"))?;
        let body = parent.body_id();
        let mut relations = Vec::new();
        let mut kept: BTreeMap<EntityId, usize> = BTreeMap::new();
        // Split children per parent, keyed (piece, local); ordinals follow
        // the keys' order (the history check's canonical ordinals).
        let mut split_keys: BTreeMap<EntityId, Vec<(usize, usize)>> = BTreeMap::new();
        let mut plans: Vec<Vec<(crate::topology::Slot, Plan)>> = Vec::new();
        let mut pieces = Vec::new();
        let mut order: Vec<usize> = (0..section.pieces.len()).collect();
        order.sort_by_key(|&i| (section.pieces[i].side, section.pieces[i].first_origin()));
        for (k, &i) in order.iter().enumerate() {
            let piece2 = &section.pieces[i];
            let mut solid = Self::build(
                operation,
                piece2.profile.clone(),
                self.frame,
                self.start,
                self.end,
            )?;
            let derive = |entity, role, ordinal, parents| Derivation {
                operation,
                kind: OperationKind::PlaneSplit,
                entity,
                role,
                ordinal,
                parents,
            };
            let mut named: Vec<(crate::topology::Slot, Plan)> = Vec::new();
            let mut child = |parent: EntityId, local: usize, entity, role| {
                split_keys.entry(parent).or_default().push((k, local));
                Plan::Child {
                    parent,
                    key: (k, local),
                    entity,
                    role,
                }
            };
            let ids: Vec<(EntityId, crate::topology::Slot)> = solid.topology.ids().collect();
            for (id, slot) in ids {
                let d = solid.topology.derivation(id).expect("a derivation").clone();
                let element = match d.parents.first() {
                    Some(Parent::Profile { boundary, element }) => {
                        Some((*boundary as usize, *element))
                    }
                    _ => None,
                };
                let new = match (d.role, element) {
                    (Role::StartCap | Role::EndCap, _) => {
                        let cap = caps[if d.role == Role::StartCap { 0 } else { 1 }];
                        child(cap, 0, d.entity, d.role)
                    }
                    (Role::Region, _) => child(region, 0, d.entity, d.role),
                    (role, Some((pb, ProfileElement::Segment(pj)))) => {
                        match piece2.segments[pb][pj as usize] {
                            SegOrigin::Whole(b, j) => {
                                let id =
                                    original(d.entity, role, b, ProfileElement::Segment(j as u32))?;
                                *kept.entry(id).or_default() += 1;
                                Plan::Fixed(parent.derivation(id).expect("a derivation").clone())
                            }
                            SegOrigin::Part(b, j, part) => {
                                let from =
                                    original(d.entity, role, b, ProfileElement::Segment(j as u32))?;
                                child(from, part, d.entity, role)
                            }
                            SegOrigin::Chord(chord) => {
                                // A cut wall, or its edge on a cap.
                                let (a, b) = section.chords[chord];
                                let walls = section
                                    .walls_at(a)
                                    .into_iter()
                                    .chain(section.walls_at(b))
                                    .map(|(b, j)| {
                                        original(
                                            EntityKind::Face,
                                            Role::Wall,
                                            b,
                                            ProfileElement::Segment(j as u32),
                                        )
                                    })
                                    .collect::<Result<BTreeSet<_>>>()?;
                                let (role, parents): (Role, Vec<EntityId>) =
                                    if d.entity == EntityKind::Face {
                                        (Role::CutFace, caps.iter().copied().chain(walls).collect())
                                    } else {
                                        let cap = if role == Role::BottomEdge { 0 } else { 1 };
                                        // Which cap is the bottom depends on the direction.
                                        let cap = if self.start < self.end { cap } else { 1 - cap };
                                        (Role::CutEdge, vec![caps[cap]])
                                    };
                                let parents: Vec<Parent> =
                                    parents.into_iter().map(Parent::Entity).collect();
                                let n = derive(
                                    d.entity,
                                    role,
                                    piece_ordinal(
                                        k,
                                        chord * 2 + (d.role == Role::TopEdge) as usize,
                                    ),
                                    parents,
                                );
                                relations.push(Relation::Generated {
                                    from: n.parents.clone(),
                                    to: n.id(),
                                    role,
                                });
                                Plan::Fixed(n)
                            }
                        }
                    }
                    (role, Some((pb, ProfileElement::Vertex(pj)))) => {
                        match piece2.vertices[pb][pj as usize] {
                            PointId::Vertex(b, j) => {
                                let id =
                                    original(d.entity, role, b, ProfileElement::Vertex(j as u32))?;
                                if section.shared.contains(&PointId::Vertex(b, j)) {
                                    child(id, 0, d.entity, role)
                                } else {
                                    *kept.entry(id).or_default() += 1;
                                    Plan::Fixed(
                                        parent.derivation(id).expect("a derivation").clone(),
                                    )
                                }
                            }
                            PointId::Cross(b, j, which) => {
                                // A new vertical edge from the wall it cuts; a
                                // new cap vertex from the cap edge it cuts.
                                let (from, role) = if d.entity == EntityKind::Edge {
                                    (
                                        original(
                                            EntityKind::Face,
                                            Role::Wall,
                                            b,
                                            ProfileElement::Segment(j as u32),
                                        )?,
                                        Role::CutEdge,
                                    )
                                } else {
                                    let edge_role = if role == Role::BottomVertex {
                                        Role::BottomEdge
                                    } else {
                                        Role::TopEdge
                                    };
                                    (
                                        original(
                                            EntityKind::Edge,
                                            edge_role,
                                            b,
                                            ProfileElement::Segment(j as u32),
                                        )?,
                                        Role::CutVertex,
                                    )
                                };
                                let n = derive(
                                    d.entity,
                                    role,
                                    piece_ordinal(k, j * 2 + which),
                                    vec![Parent::Entity(from)],
                                );
                                if !relations.iter().any(|r| matches!(r, Relation::Generated { to, .. } if *to == n.id())) {
                                    relations.push(Relation::Generated { from: n.parents.clone(), to: n.id(), role });
                                }
                                Plan::Fixed(n)
                            }
                        }
                    }
                    _ => return Err(Error::InvalidTopology("a prism slot of unknown provenance")),
                };
                named.push((slot, new));
            }
            plans.push(named);
            solid.operation = operation;
            pieces.push((piece2.side, solid));
        }
        // Phase 2: ordinals by the keys' order, then the renamed pieces. An
        // input entity kept whole keeps its id: `Unchanged` if its stored
        // geometry and bounding ids are the input's, otherwise `Modified`
        // (a bounding vertical edge split in two, a hole's wall now on an
        // outer boundary with its frame reversed).
        for keys in split_keys.values_mut() {
            keys.sort();
            keys.dedup();
        }
        let bases: Vec<crate::topology::Topology> =
            pieces.iter().map(|(_, s)| s.topology.clone()).collect();
        let before = parent.entity_set(self.resolution());
        let mut modified: BTreeSet<EntityId> = BTreeSet::new();
        let children = {
            let mut children: BTreeMap<EntityId, Vec<(u32, EntityId)>> = BTreeMap::new();
            for (k, (named, base)) in plans.iter().zip(&bases).enumerate() {
                let mut derivations = Vec::new();
                for (slot, plan) in named {
                    let d = match plan {
                        Plan::Fixed(d) => d.clone(),
                        Plan::Child {
                            parent: from,
                            key,
                            entity,
                            role,
                        } => {
                            let ordinal =
                                split_keys[from].binary_search(key).expect("a key") as u32;
                            let d = Derivation {
                                operation,
                                kind: OperationKind::PlaneSplit,
                                entity: *entity,
                                role: *role,
                                ordinal,
                                parents: vec![Parent::Entity(*from)],
                            };
                            children.entry(*from).or_default().push((ordinal, d.id()));
                            d
                        }
                    };
                    derivations.push((*slot, d));
                }
                let body_derivation = Derivation {
                    operation,
                    kind: OperationKind::PlaneSplit,
                    entity: EntityKind::Body,
                    role: Role::Body,
                    ordinal: k as u32,
                    parents: vec![Parent::Entity(body)],
                };
                pieces[k].1.topology = base.clone().renamed(body_derivation, derivations)?;
            }
            for (_, piece) in &pieces {
                let after = piece.topology.entity_set(piece.resolution());
                for (id, info) in &after.entities {
                    if !kept.contains_key(id) {
                        continue;
                    }
                    let old = &before.entities[id];
                    if old.geometry != info.geometry || old.structure != info.structure {
                        modified.insert(*id);
                    }
                }
            }
            children
        };
        for (from, mut into) in children {
            // In ordinal order, each once.
            into.sort();
            into.dedup();
            debug_assert_eq!(into.len(), split_keys[&from].len());
            relations.push(Relation::Split {
                from,
                into: into.into_iter().map(|(_, id)| id).collect(),
            });
        }
        for (id, _) in parent.ids() {
            if !kept.contains_key(&id) {
                continue;
            }
            if modified.contains(&id) {
                relations.push(Relation::Modified { from: id, to: id });
            } else {
                relations.push(Relation::Unchanged { id });
            }
        }
        relations.sort_by_cached_key(Relation::sort_key);
        relations.dedup();
        let outputs: Vec<EntityId> = pieces.iter().map(|(_, s)| s.topology.body_id()).collect();
        let history = History::new(
            operation,
            OperationKind::PlaneSplit,
            vec![body],
            outputs,
            relations,
            Vec::new(),
        );
        let mut history = history.at_level(context.level);
        {
            let mut outs: Vec<&mut Solid> = pieces.iter_mut().map(|(_, s)| s).collect();
            super::enclose::carry(&[self], &history, &mut outs);
        }
        let outs: Vec<&Solid> = pieces.iter().map(|(_, s)| s).collect();
        let maps = propagate(context, &[self], &mut history, &outs)?;
        for ((_, s), map) in pieces.iter_mut().zip(maps) {
            s.topology.set_attributes(map);
        }
        let outs: Vec<&Solid> = pieces.iter().map(|(_, s)| s).collect();
        super::stack::debug_check(&[self], &outs, &history);
        debug_check_attributes(context, &[self], &outs, &history);
        Ok((pieces, history))
    }
}

/// A slot's derivation: fixed, or a split child whose ordinal comes from
/// its key's place among its parent's children.
enum Plan {
    Fixed(Derivation),
    Child {
        parent: EntityId,
        key: (usize, usize),
        entity: EntityKind,
        role: Role,
    },
}

/// An ordinal for an entity of piece `k`: one cut in two pieces gives two
/// entities.
pub(crate) fn piece_ordinal(k: usize, local: usize) -> u32 {
    ((k as u32) << 16) | local as u32
}

/// A rational's nearest binary64.
pub(crate) fn rational_f64(x: &R) -> f64 {
    let lo = I::exact(x.clone()).bounds_f64();
    if lo.0 == lo.1 {
        return lo.0;
    }
    // The nearer of the two bounds.
    let (a, b) = (q(lo.0), q(lo.1));
    if x - &a <= &b - x {
        lo.0
    } else {
        lo.1
    }
}

// ------------------------------------------------------------------ section

/// A point of the section's arrangement: a stored profile vertex, or a
/// crossing of the line with stored segment `j` of boundary `b` (the first
/// or second along the segment's stored direction).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum PointId {
    Vertex(usize, usize),
    Cross(usize, usize, usize),
}

/// Where a piece's boundary segment comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum SegOrigin {
    /// A whole stored segment of the input.
    Whole(usize, usize),
    /// A part of one, numbered along the segment's stored direction.
    Part(usize, usize, usize),
    /// A chord of the line inside the profile.
    Chord(usize),
}

#[derive(Debug, Clone)]
pub(crate) enum Kind {
    Line,
    Arc {
        center: Point2,
        radius: f64,
        ccw: bool,
    },
    /// A spline piece, in its traversal's direction.
    Spline(Box<spline::Span>),
}

/// A directed piece of the arrangement.
#[derive(Debug, Clone)]
pub(crate) struct Edge2 {
    pub(crate) from: PointId,
    pub(crate) to: PointId,
    pub(crate) kind: Kind,
    pub(crate) origin: SegOrigin,
    /// -1 below the line, 1 above, 0 along it.
    pub(crate) side: i8,
}

pub(crate) struct SectionPiece {
    pub(crate) side: Side,
    pub(crate) profile: Profile,
    /// Per boundary of the piece's profile (outer first), per stored segment.
    pub(crate) segments: Vec<Vec<SegOrigin>>,
    pub(crate) vertices: Vec<Vec<PointId>>,
}

impl SectionPiece {
    pub(crate) fn first_origin(&self) -> SegOrigin {
        self.segments
            .iter()
            .flatten()
            .min()
            .copied()
            .expect("a segment")
    }
}

pub(crate) struct Section {
    pub(crate) pieces: Vec<SectionPiece>,
    pub(crate) chords: Vec<(PointId, PointId)>,
    /// Stored vertices on the line reached by both sides.
    pub(crate) shared: BTreeSet<PointId>,
    /// The stored segments meeting each point.
    incident: BTreeMap<PointId, Vec<(usize, usize)>>,
    /// Every point's binary64 position.
    pub(crate) positions: BTreeMap<PointId, Point2>,
    /// A spline crossing's curve parameter (S8b.3).
    params: BTreeMap<PointId, f64>,
    /// The boundary's pieces with their sides (a hole's reversed).
    pub(crate) edges: Vec<Edge2>,
}

impl Section {
    fn walls_at(&self, p: PointId) -> Vec<(usize, usize)> {
        self.incident.get(&p).cloned().unwrap_or_default()
    }

    /// The profile's pieces on each side of `a u + b v + d = 0`.
    pub(crate) fn new(profile: &Profile, line: [R; 3]) -> Result<Self> {
        Self::build(profile, line, true)
    }

    /// The boundaries' pieces and their sides only (S8e's wires): no chords,
    /// no pinch checks, no pieces.
    pub(crate) fn arrangement(profile: &Profile, line: [R; 3]) -> Result<Self> {
        Self::build(profile, line, false)
    }

    fn build(profile: &Profile, line: [R; 3], trace_pieces: bool) -> Result<Self> {
        let tolerance = profile.tolerance();
        let [a, b, d] = line;
        let f = |p: Point2| &a * q(p.x) + &b * q(p.y) + &d;
        let sign = |x: &R| x.cmp(&zero()) as i8;
        let mut positions: BTreeMap<PointId, Point2> = BTreeMap::new();
        let mut params: BTreeMap<PointId, f64> = BTreeMap::new();
        let mut on_line: BTreeMap<PointId, I> = BTreeMap::new();
        // The line's parameter: position along (-b, a).
        let along = |p: &[I; 2]| {
            p[0].mul(&I::exact(-b.clone()))
                .add(&p[1].mul(&I::exact(a.clone())))
        };
        let exact2 = |p: Point2| [I::exact(q(p.x)), I::exact(q(p.y))];
        let mut edges: Vec<Edge2> = Vec::new();
        let mut incident: BTreeMap<PointId, Vec<(usize, usize)>> = BTreeMap::new();
        // Where a circle or an arc touches the line: (position along it).
        let mut tangents: Vec<I> = Vec::new();
        let foot_along = |center: Point2| -> I {
            let ab2 = &a * &a + &b * &b;
            let fc = &a * q(center.x) + &b * q(center.y) + &d;
            let foot = [q(center.x) - &fc * &a / &ab2, q(center.y) - &fc * &b / &ab2];
            along(&[I::exact(foot[0].clone()), I::exact(foot[1].clone())])
        };
        for (bi, boundary) in profile.boundaries().enumerate() {
            let (points, segments): (Vec<Point2>, Vec<Segment>) = match &boundary.kind {
                BoundaryKind::Polygon(p) => (p.clone(), vec![Segment::Line; p.len()]),
                BoundaryKind::Path { points, segments } => (points.clone(), segments.clone()),
                BoundaryKind::Circle { center, radius } => {
                    // A circle: split at its crossings, or whole.
                    let fc = f(*center);
                    let ab2 = &a * &a + &b * &b;
                    let r2 = q(*radius) * q(*radius);
                    let delta = &fc * &fc - &r2 * &ab2;
                    if delta == zero() {
                        tangents.push(foot_along(*center));
                    }
                    if delta >= zero() {
                        // Missing or tangent: whole, on the centre's side;
                        // its placeholder vertex a point of it.
                        positions.insert(
                            PointId::Vertex(bi, 0),
                            Point2::new(center.x + radius, center.y),
                        );
                        edges.push(Edge2 {
                            from: PointId::Vertex(bi, 0),
                            to: PointId::Vertex(bi, 0),
                            kind: Kind::Arc {
                                center: *center,
                                radius: *radius,
                                ccw: bi == 0,
                            },
                            origin: SegOrigin::Whole(bi, 0),
                            side: sign(&fc),
                        });
                        continue;
                    }
                    let cross =
                        circle_points(*center, *radius, &[a.clone(), b.clone(), d.clone()])?;
                    let mut ids = Vec::new();
                    for (k, p) in cross.iter().enumerate() {
                        let id = PointId::Cross(bi, 0, k);
                        positions.insert(id, rounded(p));
                        on_line.insert(id, along(p));
                        incident.entry(id).or_default().push((bi, 0));
                        ids.push(id);
                    }
                    // Two arcs, each on one side (region-left: a hole's
                    // circle runs clockwise).
                    let ccw = bi == 0;
                    for (k, (s, e)) in [(ids[0], ids[1]), (ids[1], ids[0])].into_iter().enumerate()
                    {
                        let (from, to) = if ccw { (s, e) } else { (e, s) };
                        let side = arc_side(
                            *center,
                            *radius,
                            positions[&from],
                            positions[&to],
                            ccw,
                            &[a.clone(), b.clone(), d.clone()],
                        )?;
                        edges.push(Edge2 {
                            from,
                            to,
                            kind: Kind::Arc {
                                center: *center,
                                radius: *radius,
                                ccw,
                            },
                            origin: SegOrigin::Part(bi, 0, k),
                            side,
                        });
                    }
                    continue;
                }
            };
            let n = points.len();
            let reach = q(tolerance.linear());
            for (j, p) in points.iter().enumerate() {
                positions.insert(PointId::Vertex(bi, j), *p);
                let fp = f(*p);
                if fp == zero() {
                    on_line.insert(PointId::Vertex(bi, j), along(&exact2(*p)));
                } else if &fp * &fp <= &reach * &reach * (&a * &a + &b * &b) {
                    // Off the line by less than the resolution: its crossings
                    // and pieces would be thinner than it.
                    return Err(Error::Degenerate(
                        "a split within the resolution of a profile vertex",
                    ));
                }
                incident
                    .entry(PointId::Vertex(bi, j))
                    .or_default()
                    .push((bi, j));
                incident
                    .entry(PointId::Vertex(bi, (j + 1) % n))
                    .or_default()
                    .push((bi, j));
            }
            let mut boundary_edges = Vec::new();
            for j in 0..n {
                let (p, e) = (points[j], points[(j + 1) % n]);
                let (fp, fe) = (f(p), f(e));
                let mut cuts: Vec<PointId> = Vec::new();
                // A spline's pieces and their sides, along the segment.
                let mut spline_parts: Option<(Vec<Kind>, Vec<i8>)> = None;
                let kind = match segments[j].clone() {
                    Segment::Spline(span) => {
                        let meeting = spline::meets(&span, [&a, &b, &d])?;
                        for t in &meeting.touches {
                            tangents.push(along(t));
                        }
                        let (first, last) = span.curve().as_curve3().domain();
                        // The pieces' curve parameters, along the segment.
                        let mut ts: Vec<f64> = meeting.crossings.iter().map(|c| c.t).collect();
                        if span.is_reversed() {
                            ts.insert(0, last);
                            ts.push(first);
                        } else {
                            ts.insert(0, first);
                            ts.push(last);
                        }
                        let mut kinds = Vec::new();
                        for w in ts.windows(2) {
                            let (t0, t1) = (w[0].min(w[1]), w[0].max(w[1]));
                            kinds.push(Kind::Spline(Box::new(spline::restrict(&span, t0, t1)?)));
                        }
                        for (k, c) in meeting.crossings.iter().enumerate() {
                            let id = PointId::Cross(bi, j, k);
                            positions.insert(id, c.rounded);
                            params.insert(id, c.t);
                            on_line.insert(
                                id,
                                along(&[
                                    I::exact(c.exact[0].clone()),
                                    I::exact(c.exact[1].clone()),
                                ]),
                            );
                            incident.entry(id).or_default().push((bi, j));
                            cuts.push(id);
                        }
                        spline_parts = Some((kinds, meeting.sides));
                        Kind::Line
                    }
                    Segment::Line => {
                        if sign(&fp) * sign(&fe) < 0 {
                            let t = &fp / (&fp - &fe);
                            let x = I::exact(q(p.x))
                                .add(&I::exact(t.clone()).mul(&I::exact(q(e.x) - q(p.x))));
                            let y =
                                I::exact(q(p.y)).add(&I::exact(t).mul(&I::exact(q(e.y) - q(p.y))));
                            let id = PointId::Cross(bi, j, 0);
                            positions.insert(id, rounded(&[x.clone(), y.clone()]));
                            on_line.insert(id, along(&[x, y]));
                            incident.entry(id).or_default().push((bi, j));
                            cuts.push(id);
                        }
                        Kind::Line
                    }
                    Segment::Arc {
                        center,
                        radius,
                        ccw,
                    } => {
                        let fc = f(center);
                        let ab2 = &a * &a + &b * &b;
                        let delta = &fc * &fc - q(radius) * q(radius) * &ab2;
                        if delta == zero() {
                            let foot =
                                [q(center.x) - &fc * &a / &ab2, q(center.y) - &fc * &b / &ab2];
                            let t = [I::exact(foot[0].clone()), I::exact(foot[1].clone())];
                            let at_end = same_point(&t, p)? || same_point(&t, e)?;
                            if !at_end && within_arc(center, p, e, ccw, &t)?.is_some() {
                                tangents.push(along(&t));
                            }
                        }
                        if delta < zero() {
                            let pts =
                                circle_points(center, radius, &[a.clone(), b.clone(), d.clone()])?;
                            let mut inside: Vec<(ArcPos, [I; 2])> = Vec::new();
                            for pt in pts {
                                // Endpoints on the line are the path's vertices.
                                if ((fp == zero()) && same_point(&pt, p)?)
                                    || ((fe == zero()) && same_point(&pt, e)?)
                                {
                                    continue;
                                }
                                if let Some(t) = within_arc(center, p, e, ccw, &pt)? {
                                    inside.push((t, pt));
                                }
                            }
                            inside.sort_by(|x, y| x.0.compare(&y.0).unwrap_or(Ordering::Equal));
                            for (k, (_, pt)) in inside.into_iter().enumerate() {
                                let id = PointId::Cross(bi, j, k);
                                positions.insert(id, rounded(&pt));
                                on_line.insert(id, along(&pt));
                                incident.entry(id).or_default().push((bi, j));
                                cuts.push(id);
                            }
                        }
                        Kind::Arc {
                            center,
                            radius,
                            ccw,
                        }
                    }
                };
                // The stored direction's parts, then (for a hole) reversed.
                let mut chain = vec![PointId::Vertex(bi, j)];
                chain.extend(cuts);
                chain.push(PointId::Vertex(bi, (j + 1) % n));
                let parts = chain.len() - 1;
                for k in 0..parts {
                    let (s, e) = (chain[k], chain[k + 1]);
                    let origin = if parts == 1 {
                        SegOrigin::Whole(bi, j)
                    } else {
                        SegOrigin::Part(bi, j, k)
                    };
                    if let Some((kinds, sides)) = &spline_parts {
                        boundary_edges.push(Edge2 {
                            from: s,
                            to: e,
                            kind: kinds[k].clone(),
                            origin,
                            side: sides[k],
                        });
                        continue;
                    }
                    let side = match kind.clone() {
                        Kind::Line => {
                            let mid = [
                                positions[&s].x * 0.5 + positions[&e].x * 0.5,
                                positions[&s].y * 0.5 + positions[&e].y * 0.5,
                            ];
                            // Ends on the line, or strictly on one side: the
                            // exact ends' signs decide unless both are cuts.
                            let fs = exact_sign(&s, &positions, &on_line, &f);
                            let fe2 = exact_sign(&e, &positions, &on_line, &f);
                            if fs != 0 {
                                fs
                            } else if fe2 != 0 {
                                fe2
                            } else if matches!((s, e), (PointId::Vertex(..), PointId::Vertex(..))) {
                                0
                            } else {
                                sign(&f(Point2::new(mid[0], mid[1])))
                            }
                        }
                        Kind::Arc {
                            center,
                            radius,
                            ccw,
                        } => arc_side(
                            center,
                            radius,
                            positions[&s],
                            positions[&e],
                            ccw,
                            &[a.clone(), b.clone(), d.clone()],
                        )?,
                        Kind::Spline(_) => unreachable!("a spline's pieces are pushed above"),
                    };
                    boundary_edges.push(Edge2 {
                        from: s,
                        to: e,
                        kind: kind.clone(),
                        origin,
                        side,
                    });
                }
            }
            if bi > 0 {
                // A hole runs clockwise with the region on its left.
                boundary_edges.reverse();
                for e in &mut boundary_edges {
                    std::mem::swap(&mut e.from, &mut e.to);
                    match &mut e.kind {
                        Kind::Arc { ccw, .. } => *ccw = !*ccw,
                        Kind::Spline(span) => **span = span.reversed(),
                        Kind::Line => {}
                    }
                }
            }
            edges.extend(boundary_edges);
        }
        // A crossing within the resolution of its segment's ends, or two
        // points on the line within it of each other, would leave a piece
        // or an edge thinner than the resolution: not representable.
        let tol = tolerance.linear();
        for (id, p) in &positions {
            if let PointId::Cross(bi, j, _) = id {
                for (other, q) in &positions {
                    let same_segment = match other {
                        PointId::Vertex(ob, oj) => {
                            *ob == *bi
                                && (*oj == *j
                                    || (incident
                                        .get(other)
                                        .is_some_and(|s| s.contains(&(*bi, *j)))))
                        }
                        PointId::Cross(ob, oj, _) => ob == bi && oj == j && other != id,
                    };
                    if same_segment && p.distance(*q) <= tol {
                        return Err(Error::Degenerate(
                            "a split within the resolution of a profile vertex",
                        ));
                    }
                }
            }
        }
        if !trace_pieces {
            return Ok(Section {
                pieces: Vec::new(),
                chords: Vec::new(),
                shared: BTreeSet::new(),
                incident,
                positions,
                params,
                edges,
            });
        }
        // Chords: between consecutive points on the line, inside the profile.
        let on_line_at: BTreeMap<PointId, I> = on_line.clone();
        let mut along_line: Vec<(PointId, I)> = on_line.into_iter().collect();
        along_line.sort_by(|x, y| {
            x.1.cmp(&y.1)
                .unwrap_or_else(|| x.1.midpoint().cmp(&y.1.midpoint()))
        });
        let mut chords = Vec::new();
        for w in along_line.windows(2) {
            let (s, e) = (w[0].0, w[1].0);
            let (ps, pe) = (positions[&s], positions[&e]);
            if ps != pe && ps.distance(pe) <= tol {
                return Err(Error::Degenerate("section points within the resolution"));
            }
            // Off the chord's centre: a symmetric touch point sits there.
            let f0 = 0.4453125;
            let mid = Point2::new(ps.x + (pe.x - ps.x) * f0, ps.y + (pe.y - ps.y) * f0);
            if ps == pe {
                continue;
            }
            if profile.classify(mid)? == Location::Inside {
                chords.push((s, e));
            }
        }
        // A point where the boundary only touches the line (a tangency, or a
        // vertex whose pieces lie on one side) with the profile's inside
        // along the line on both of its sides pinches a piece: its hole
        // touches its boundary there, which a profile cannot hold.
        for t in &tangents {
            for (s, e) in &chords {
                let (ts, te) = (&on_line_at[s], &on_line_at[e]);
                let inside = |x: &I, y: &I| x.cmp(y) == Some(Ordering::Less);
                if (inside(ts, t) && inside(t, te)) || (inside(te, t) && inside(t, ts)) {
                    return Err(Error::Degenerate("a plane tangent to the solid inside it"));
                }
            }
        }
        for (p, _) in &along_line {
            if !matches!(p, PointId::Vertex(..)) {
                continue;
            }
            let sides: BTreeSet<i8> = edges
                .iter()
                .filter(|e| e.from == *p || e.to == *p)
                .map(|e| e.side)
                .collect();
            let touch = sides.len() == 1 && !sides.contains(&0);
            let chorded = chords.iter().filter(|(s, e)| s == p || e == p).count();
            if touch && chorded >= 2 {
                return Err(Error::Degenerate(
                    "a plane touching the solid's boundary inside it",
                ));
            }
        }
        // Stored vertices on the line reached from both sides.
        let mut sides_at: BTreeMap<PointId, BTreeSet<i8>> = BTreeMap::new();
        for e in &edges {
            for p in [e.from, e.to] {
                sides_at.entry(p).or_default().insert(e.side);
            }
        }
        for (s, e) in &chords {
            for p in [s, e] {
                sides_at.entry(*p).or_default().extend([-1, 1]);
            }
        }
        let shared: BTreeSet<PointId> = sides_at
            .iter()
            .filter(|(p, s)| matches!(p, PointId::Vertex(..)) && s.contains(&-1) && s.contains(&1))
            .map(|(p, _)| *p)
            .collect();
        let pieces = trace(&edges, &chords, &positions, &[a, b], tolerance)?;
        Ok(Section {
            pieces,
            chords,
            shared,
            incident,
            positions,
            params,
            edges,
        })
    }
}

/// The sign of the line's function at a point of the arrangement: a
/// crossing is on the line; a stored vertex exactly.
fn exact_sign(
    p: &PointId,
    positions: &BTreeMap<PointId, Point2>,
    on_line: &BTreeMap<PointId, I>,
    f: &dyn Fn(Point2) -> R,
) -> i8 {
    if on_line.contains_key(p) || matches!(p, PointId::Cross(..)) {
        return 0;
    }
    f(positions[p]).cmp(&zero()) as i8
}

pub(crate) fn rounded(p: &[I; 2]) -> Point2 {
    let (x0, x1) = p[0].bounds_f64();
    let (y0, y1) = p[1].bounds_f64();
    Point2::new(0.5 * x0 + 0.5 * x1, 0.5 * y0 + 0.5 * y1)
}

/// Whether an enclosed circle point is the exact vertex `v` (a vertex on
/// the line is one of the circle's two points; the other is apart).
pub(crate) fn same_point(p: &[I; 2], v: Point2) -> Result<bool> {
    let contains = |iv: &I, x: f64| iv.lo() <= &q(x) && &q(x) <= iv.hi();
    Ok(contains(&p[0], v.x) && contains(&p[1], v.y))
}

/// The two points of a circle on the line (the line certainly crossing it).
pub(crate) fn circle_points(center: Point2, radius: f64, line: &[R; 3]) -> Result<Vec<[I; 2]>> {
    let [a, b, d] = line;
    let ab2 = a * a + b * b;
    let fc = a * q(center.x) + b * q(center.y) + d;
    // Foot: c - fc (a, b) / |ab|^2; offset along (-b, a) by
    // sqrt(r^2 |ab|^2 - fc^2) / |ab|^2.
    let foot = [q(center.x) - &fc * a / &ab2, q(center.y) - &fc * b / &ab2];
    let h2 = q(radius) * q(radius) * &ab2 - &fc * &fc;
    let h = I::exact(h2)
        .sqrt()
        .div(&I::exact(ab2.clone()))
        .ok_or(Error::Degenerate("a line"))?;
    let mut out = Vec::new();
    for s in [1i64, -1] {
        let s = I::exact(R::from_integer(s.into()));
        out.push([
            I::exact(foot[0].clone()).add(&s.mul(&h).mul(&I::exact(-b.clone()))),
            I::exact(foot[1].clone()).add(&s.mul(&h).mul(&I::exact(a.clone()))),
        ]);
    }
    Ok(out)
}

/// A point's place along an arc from its start: the half-turn it lies in
/// (0 the start's direction, 1 the first half-turn, 2 the opposite
/// direction, 3 the second half-turn) and its vector from the centre, in
/// the arc's turning sense. Places compare by orientation signs alone.
#[derive(Debug, Clone)]
pub(crate) struct ArcPos {
    half: u8,
    v: [I; 2],
    sense: i8,
}

impl ArcPos {
    pub(crate) fn of(r: &[I; 2], v: [I; 2], sense: i8) -> Result<Self> {
        let limit = || Error::ComputationLimit("an arc's crossing at its end");
        let cross = r[0].mul(&v[1]).sub(&r[1].mul(&v[0]));
        let cross = if sense > 0 { cross } else { cross.neg() };
        let half = match cross.sign().ok_or_else(limit)? {
            Ordering::Greater => 1,
            Ordering::Less => 3,
            Ordering::Equal => {
                let dot = r[0].mul(&v[0]).add(&r[1].mul(&v[1]));
                match dot.sign().ok_or_else(limit)? {
                    Ordering::Less => 2,
                    _ => 0,
                }
            }
        };
        Ok(Self { half, v, sense })
    }
    /// Along the arc: the earlier first (within a half-turn, the one the
    /// other lies counter to the sense from).
    pub(crate) fn compare(&self, o: &Self) -> Option<Ordering> {
        if self.half != o.half {
            return Some(self.half.cmp(&o.half));
        }
        if self.half % 2 == 0 {
            return Some(Ordering::Equal);
        }
        let cross = self.v[0].mul(&o.v[1]).sub(&self.v[1].mul(&o.v[0]));
        let cross = if self.sense > 0 { cross } else { cross.neg() };
        Some(cross.sign()?.reverse())
    }
}

/// Whether a point of an arc's circle lies strictly inside the arc from `p`
/// to `e`; its place along the arc if it does. Exact orientation signs, no
/// angles.
pub(crate) fn within_arc(
    center: Point2,
    p: Point2,
    e: Point2,
    ccw: bool,
    x: &[I; 2],
) -> Result<Option<ArcPos>> {
    let c = [I::exact(q(center.x)), I::exact(q(center.y))];
    let rel = |v: [I; 2]| [v[0].sub(&c[0]), v[1].sub(&c[1])];
    let (vp, ve, vx) = (
        rel([I::exact(q(p.x)), I::exact(q(p.y))]),
        rel([I::exact(q(e.x)), I::exact(q(e.y))]),
        rel(x.clone()),
    );
    let sense = if ccw { 1 } else { -1 };
    let px = ArcPos::of(&vp, vx, sense)?;
    if px.half == 0 {
        return Ok(None);
    }
    // A full circle from a point back to it: all but the point.
    if p == e {
        return Ok(Some(px));
    }
    let pe = ArcPos::of(&vp, ve, sense)?;
    match px.compare(&pe) {
        Some(Ordering::Less) => Ok(Some(px)),
        Some(_) => Ok(None),
        None => Err(Error::ComputationLimit("an arc's crossing at its end")),
    }
}

/// Whether the prism lies (weakly) on one side of `a u + b v + c w + d`:
/// the function's extremes over its caps' outer boundaries (vertices
/// exactly; an arc's support point, `f(centre) +- r |(a, b)|`, by the exact
/// sign of a quadratic surd when it lies inside the arc) do not change sign.
fn on_one_side(profile: &Profile, plane: [&R; 4], low: &R, high: &R) -> Result<bool> {
    let [a, b, c, d] = plane;
    let (w_min, w_max) = if c >= &zero() {
        (c * low, c * high)
    } else {
        (c * high, c * low)
    };
    let s = a * a + b * b;
    // The sign of K + t R, R = r sqrt(s) >= 0, t = +-1.
    let surd_sign = |k: &R, r: &R, t: i8| -> Ordering {
        let rr = r * r * &s;
        let kk = k * k;
        match (k.cmp(&zero()), t) {
            (Ordering::Equal, _) => (t as i32).cmp(&0),
            (Ordering::Greater, 1) => Ordering::Greater,
            (Ordering::Less, -1) => Ordering::Less,
            // K - R.
            (Ordering::Greater, _) => kk.cmp(&rr),
            // R - |K|.
            (_, _) => rr.cmp(&kk),
        }
    };
    let f = |p: Point2| a * q(p.x) + b * q(p.y) + d;
    // (sign at the top of the extent, sign at the bottom) per candidate.
    let mut signs: Vec<(Ordering, Ordering)> = Vec::new();
    let vertex = |p: Point2| ((f(p) + &w_max).cmp(&zero()), (f(p) + &w_min).cmp(&zero()));
    let arc = |center: Point2,
               radius: f64,
               ends: Option<(Point2, Point2, bool)>|
     -> Result<Vec<(Ordering, Ordering)>> {
        let mut out = Vec::new();
        if s == zero() {
            return Ok(out);
        }
        let k = f(center);
        for t in [1i8, -1] {
            if let Some((p, e, ccw)) = ends {
                // The support point: centre + t r (a, b) / |(a, b)|.
                let norm = I::exact(s.clone()).sqrt();
                let dir = [I::exact(a.clone()), I::exact(b.clone())].map(|x| {
                    x.mul(&I::exact_f64(radius * t as f64))
                        .div(&norm)
                        .expect("a nonzero normal")
                });
                let pt = [
                    I::exact(q(center.x)).add(&dir[0]),
                    I::exact(q(center.y)).add(&dir[1]),
                ];
                let at_end = same_point(&pt, p)? || same_point(&pt, e)?;
                if at_end || within_arc(center, p, e, ccw, &pt)?.is_none() {
                    continue;
                }
            }
            out.push((
                surd_sign(&(&k + &w_max), &q(radius), t),
                surd_sign(&(&k + &w_min), &q(radius), t),
            ));
        }
        Ok(out)
    };
    match &profile.outer().kind {
        BoundaryKind::Circle { center, radius } => {
            signs.push(vertex(*center));
            signs.extend(arc(*center, *radius, None)?);
        }
        BoundaryKind::Polygon(points) => signs.extend(points.iter().map(|p| vertex(*p))),
        BoundaryKind::Path { points, segments } => {
            signs.extend(points.iter().map(|p| vertex(*p)));
            let n = points.len();
            for (j, seg) in segments.iter().enumerate() {
                match seg {
                    Segment::Arc {
                        center,
                        radius,
                        ccw,
                    } => signs.extend(arc(
                        *center,
                        *radius,
                        Some((points[j], points[(j + 1) % n], *ccw)),
                    )?),
                    // A spline: whether the function takes each sign on it,
                    // exactly between its roots.
                    Segment::Spline(span) => {
                        let (top, bottom) = (d + &w_max, d + &w_min);
                        let above = spline::takes(span, [a, b, &top], Ordering::Greater)?;
                        let below = spline::takes(span, [a, b, &bottom], Ordering::Less)?;
                        signs.push((
                            if above {
                                Ordering::Greater
                            } else {
                                Ordering::Equal
                            },
                            if below {
                                Ordering::Less
                            } else {
                                Ordering::Equal
                            },
                        ));
                    }
                    Segment::Line => {}
                }
            }
        }
    }
    let above = signs.iter().any(|x| x.0 == Ordering::Greater);
    let below = signs.iter().any(|x| x.1 == Ordering::Less);
    Ok(!(above && below))
}

/// The side of an arc piece from `s` to `e` (no crossing inside it): the
/// line's function at the piece's angular midpoint, certainly signed.
fn arc_side(
    center: Point2,
    radius: f64,
    s: Point2,
    e: Point2,
    ccw: bool,
    line: &[R; 3],
) -> Result<i8> {
    let angle = |p: Point2| (p.y - center.y).atan2(p.x - center.x);
    let (a0, a1) = (angle(s), angle(e));
    let mut sweep = if ccw { a1 - a0 } else { a0 - a1 };
    if sweep <= 0.0 {
        sweep += std::f64::consts::TAU;
    }
    let mid = if ccw {
        a0 + sweep / 2.0
    } else {
        a0 - sweep / 2.0
    };
    let (sn, cs) = mid.sin_cos();
    let m = Point2::new(center.x + radius * cs, center.y + radius * sn);
    let [a, b, d] = line;
    let v = a * q(m.x) + b * q(m.y) + d;
    // The midpoint is far from the ends, so far from the line unless the arc
    // touches it tangentially there: then its side is its centre's opposite.
    Ok(match v.cmp(&zero()) {
        // Touching the line there, the arc lies on its centre's side.
        Ordering::Equal => (a * q(center.x) + b * q(center.y) + d).cmp(&zero()) as i8,
        o => o as i8,
    })
}

/// A traced boundary with its segments' and vertices' provenance.
type Traced = (Boundary, Vec<SegOrigin>, Vec<PointId>);

/// Trace each side's cycles into pieces.
fn trace(
    edges: &[Edge2],
    chords: &[(PointId, PointId)],
    positions: &BTreeMap<PointId, Point2>,
    normal: &[R; 2],
    tolerance: Tolerance,
) -> Result<Vec<SectionPiece>> {
    let (na, nb) = (normal[0].clone(), normal[1].clone());
    let mut out = Vec::new();
    for (side, sgn) in [(Side::Below, -1i8), (Side::Above, 1i8)] {
        let mut set: Vec<Edge2> = Vec::new();
        for e in edges {
            let include = match e.side {
                0 => {
                    // Along the line: its region (left) side.
                    let (ps, pe) = (positions[&e.from], positions[&e.to]);
                    let left = [q(ps.y) - q(pe.y), q(pe.x) - q(ps.x)];
                    let s = (&na * &left[0] + &nb * &left[1]).cmp(&zero()) as i8;
                    s == sgn
                }
                s => s == sgn,
            };
            if include {
                set.push(e.clone());
            }
        }
        for (k, &(s, e)) in chords.iter().enumerate() {
            // Directed with this side on its left.
            let (ps, pe) = (positions[&s], positions[&e]);
            let left = [q(ps.y) - q(pe.y), q(pe.x) - q(ps.x)];
            let forward = (&na * &left[0] + &nb * &left[1]).cmp(&zero()) as i8 == sgn;
            let (from, to) = if forward { (s, e) } else { (e, s) };
            set.push(Edge2 {
                from,
                to,
                kind: Kind::Line,
                origin: SegOrigin::Chord(k),
                side: sgn,
            });
        }
        let mut used = vec![false; set.len()];
        let mut cycles: Vec<Vec<usize>> = Vec::new();
        for start in 0..set.len() {
            if used[start] {
                continue;
            }
            let mut cycle = vec![start];
            used[start] = true;
            let mut current = start;
            loop {
                let at = set[current].to;
                if at == set[start].from {
                    break;
                }
                let incoming = tangent(&set[current], positions, false);
                let back = incoming.1.atan2(incoming.0) + std::f64::consts::PI;
                let mut best: Option<(f64, usize)> = None;
                for (i, e) in set.iter().enumerate() {
                    if e.from != at || used[i] {
                        continue;
                    }
                    let d = tangent(e, positions, true);
                    let mut cw = back - d.1.atan2(d.0);
                    while cw <= 1e-12 {
                        cw += std::f64::consts::TAU;
                    }
                    while cw > std::f64::consts::TAU + 1e-12 {
                        cw -= std::f64::consts::TAU;
                    }
                    if best.is_none_or(|b| cw < b.0) {
                        best = Some((cw, i));
                    }
                }
                let Some((_, next)) = best else {
                    return Err(Error::InvalidTopology("an open section cycle"));
                };
                used[next] = true;
                cycle.push(next);
                current = next;
                if cycle.len() > set.len() {
                    return Err(Error::InvalidTopology("a section cycle does not close"));
                }
            }
            cycles.push(cycle);
        }
        // Outer cycles (counter-clockwise) and holes (clockwise).
        let mut outers = Vec::new();
        let mut holes = Vec::new();
        for c in cycles {
            let area = signed_area(
                &c.iter().map(|&i| set[i].clone()).collect::<Vec<_>>(),
                positions,
            );
            if area > 0.0 {
                outers.push(c);
            } else {
                holes.push(c);
            }
        }
        let mut built: Vec<(Traced, Vec<Traced>)> = Vec::new();
        for c in &outers {
            let (bd, segs, verts) = boundary_of(
                c.iter().map(|&i| &set[i]).collect(),
                positions,
                tolerance,
                false,
            )?;
            built.push(((bd, segs, verts), Vec::new()));
        }
        for c in &holes {
            let hole = boundary_of(
                c.iter().map(|&i| &set[i]).collect(),
                positions,
                tolerance,
                true,
            )?;
            let sample = positions[&set[c[0]].from];
            let owner = match built
                .iter()
                .position(|((bd, ..), _)| bd.locate(sample, tolerance) == Location::Inside)
            {
                Some(owner) => owner,
                // On a piece's boundary within the resolution: the hole
                // touches it there (a line grazing a hole's circle).
                None if built
                    .iter()
                    .any(|((bd, ..), _)| bd.locate(sample, tolerance) == Location::Boundary) =>
                {
                    return Err(Error::Degenerate(
                        "a piece's hole within the resolution of its boundary",
                    ));
                }
                None => return Err(Error::InvalidTopology("a section hole outside every piece")),
            };
            built[owner].1.push(hole);
        }
        for ((outer, segs, verts), hs) in built {
            let mut segments = vec![segs];
            let mut vertices = vec![verts];
            let mut hole_boundaries = Vec::new();
            for (h, s, v) in hs {
                hole_boundaries.push(h);
                segments.push(s);
                vertices.push(v);
            }
            // A hole within the resolution of its piece's boundary (a plane
            // passing that close to a hole's vertex) is not representable.
            let profile = Profile::new(outer, hole_boundaries, tolerance).map_err(|e| match e {
                Error::InvalidHole(_) => {
                    Error::Degenerate("a piece's hole within the resolution of its boundary")
                }
                e => e,
            })?;
            out.push(SectionPiece {
                side,
                profile,
                segments,
                vertices,
            });
        }
    }
    Ok(out)
}

/// A piece's direction at its start (`start`) or end, as a vector.
fn tangent(e: &Edge2, positions: &BTreeMap<PointId, Point2>, start: bool) -> (f64, f64) {
    let (ps, pe) = (positions[&e.from], positions[&e.to]);
    match &e.kind {
        Kind::Line => (pe.x - ps.x, pe.y - ps.y),
        Kind::Spline(span) => spline::tangent(span, start),
        &Kind::Arc { center, ccw, .. } => {
            let p = if start { ps } else { pe };
            let (rx, ry) = (p.x - center.x, p.y - center.y);
            if ccw {
                (-ry, rx)
            } else {
                (ry, -rx)
            }
        }
    }
}

fn signed_area(cycle: &[Edge2], positions: &BTreeMap<PointId, Point2>) -> f64 {
    let mut twice = 0.0;
    for e in cycle {
        let (p, r) = (positions[&e.from], positions[&e.to]);
        twice += p.x * r.y - r.x * p.y;
        if let Kind::Spline(span) = &e.kind {
            twice += spline::twice_area_beyond_chord(span);
        }
        if let Kind::Arc {
            center,
            radius,
            ccw,
        } = e.kind
        {
            let a0 = (p.y - center.y).atan2(p.x - center.x);
            let a1 = (r.y - center.y).atan2(r.x - center.x);
            let mut phi = if ccw { a1 - a0 } else { a0 - a1 };
            if phi <= 0.0 {
                phi += std::f64::consts::TAU;
            }
            let phi = if ccw { phi } else { -phi };
            twice += radius * radius * (phi - phi.sin());
        }
    }
    twice / 2.0
}

/// A cycle as a stored counter-clockwise boundary (a hole's reversed), with
/// its segments' and vertices' provenance in stored order.
fn boundary_of(
    cycle: Vec<&Edge2>,
    positions: &BTreeMap<PointId, Point2>,
    tolerance: Tolerance,
    hole: bool,
) -> Result<(Boundary, Vec<SegOrigin>, Vec<PointId>)> {
    // A single whole circle.
    if cycle.len() == 1 && cycle[0].from == cycle[0].to {
        if let Kind::Arc { center, radius, .. } = cycle[0].kind.clone() {
            return Ok((
                Boundary::circle(center, radius, tolerance)?,
                vec![cycle[0].origin],
                vec![],
            ));
        }
    }
    let mut items: Vec<(PointId, Segment, SegOrigin)> = cycle
        .iter()
        .map(|e| {
            let seg = match e.kind.clone() {
                Kind::Line => Segment::Line,
                Kind::Spline(span) => Segment::Spline(*span),
                Kind::Arc {
                    center,
                    radius,
                    ccw,
                } => Segment::Arc {
                    center,
                    radius,
                    ccw,
                },
            };
            (e.from, seg, e.origin)
        })
        .collect();
    if hole {
        // Stored counter-clockwise: reverse the traversal.
        let n = items.len();
        let mut rev = Vec::with_capacity(n);
        for k in 0..n {
            // Segment from point k+1 back to k, reversed.
            let (_, seg, origin) = items[(n - 1 - k) % n].clone();
            let from = items[(n - k) % n].0;
            let seg = match seg {
                Segment::Spline(span) => Segment::Spline(span.reversed()),
                Segment::Line => Segment::Line,
                Segment::Arc {
                    center,
                    radius,
                    ccw,
                } => Segment::Arc {
                    center,
                    radius,
                    ccw: !ccw,
                },
            };
            rev.push((from, seg, origin));
        }
        items = rev;
    }
    let points: Vec<Point2> = items.iter().map(|i| positions[&i.0]).collect();
    let segments: Vec<Segment> = items.iter().map(|i| i.1.clone()).collect();
    // A piece whose boundary touches itself within the resolution (a line
    // grazing a circle, cutting a sliver off it) is thinner than it.
    let boundary = Boundary::path(points.clone(), segments, tolerance).map_err(|e| match e {
        Error::SelfIntersection => Error::Degenerate("a piece thinner than the resolution"),
        e => e,
    })?;
    let stored: Vec<Point2> = match &boundary.kind {
        BoundaryKind::Polygon(p) => p.clone(),
        BoundaryKind::Path { points, .. } => points.clone(),
        BoundaryKind::Circle { .. } => Vec::new(),
    };
    if stored != points {
        return Err(Error::InvalidTopology(
            "a section boundary stored in another order",
        ));
    }
    Ok((
        boundary,
        items.iter().map(|i| i.2).collect(),
        items.iter().map(|i| i.0).collect(),
    ))
}
