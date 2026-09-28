//! S8e: sheets and wires split by a plane (`SPLIT.md`).
//!
//! The plane's trace on the body's plane is `a u + b v + d = 0` in the
//! frame's coordinates, its coefficients exact rationals of the stored data.
//! A sheet's pieces are its profile's section by the trace (S8a.1's and
//! S8b.3's exact `Section`), each a planar sheet renamed from the input by
//! provenance as a prism's caps are; a closed wire's pieces are the maximal
//! runs of its boundary's pieces on one side (an edge along the trace joins
//! the run it continues in stored order), each an open wire.
use super::{Body, Construction};
use crate::history::{self, History, Relation};
use crate::identity::{
    Derivation, EntityId, EntityKind, OperationId, OperationKind, Parent, ProfileElement, Role,
};
use crate::profile::Segment;
use crate::solid::split::{piece_ordinal, q, zero, Kind, PointId, Section, SegOrigin};
use crate::solid::{replayable, Context};
use crate::topology::{Slot, Topology};
use crate::{Boundary, Error, Frame3, Point2, Profile, Result, Side, Tolerance};
use num_rational::BigRational as R;
use std::collections::{BTreeMap, BTreeSet};

/// A slot's derivation: fixed, or a split child whose ordinal comes from its
/// key's place among its parent's children.
enum Plan {
    Fixed(Derivation),
    Child {
        parent: EntityId,
        key: (usize, usize),
        entity: EntityKind,
        role: Role,
    },
}

/// A piece before naming: its side, its body and each slot's plan.
type Unnamed = (Side, Body, Vec<(Slot, Plan)>);

/// The input's entities by (entity, role, the profile element they derive
/// from), and each boundary's parents.
struct Origins {
    by_origin: BTreeMap<(EntityKind, Role, Vec<Parent>), EntityId>,
    labels: Vec<Option<crate::BoundaryLabels>>,
}

impl Origins {
    fn new(topology: &Topology, boundaries: Vec<&Boundary>) -> Self {
        let mut by_origin = BTreeMap::new();
        for (id, _) in topology.ids() {
            let d = topology.derivation(id).expect("a derivation");
            by_origin.insert((d.entity, d.role, d.parents.clone()), id);
        }
        Self {
            by_origin,
            labels: boundaries.iter().map(|b| b.labels().cloned()).collect(),
        }
    }
    fn parent(&self, b: usize, e: ProfileElement) -> Parent {
        match &self.labels[b] {
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
    }
    fn get(&self, entity: EntityKind, role: Role, b: usize, e: ProfileElement) -> Result<EntityId> {
        self.by_origin
            .get(&(entity, role, vec![self.parent(b, e)]))
            .copied()
            .ok_or(Error::InvalidTopology(
                "a body entity without its profile element",
            ))
    }
}

/// Where a piece's slot came from, before naming.
struct Naming {
    operation: OperationId,
    relations: Vec<Relation>,
    /// Input entities kept whole in a piece.
    kept: BTreeSet<EntityId>,
    /// Split children per parent, keyed (piece, local).
    split_keys: BTreeMap<EntityId, Vec<(usize, usize)>>,
}

impl Naming {
    fn derive(
        &self,
        entity: EntityKind,
        role: Role,
        ordinal: u32,
        parents: Vec<Parent>,
    ) -> Derivation {
        Derivation {
            operation: self.operation,
            kind: OperationKind::PlaneSplit,
            entity,
            role,
            ordinal,
            parents,
        }
    }
    fn child(
        &mut self,
        parent: EntityId,
        key: (usize, usize),
        entity: EntityKind,
        role: Role,
    ) -> Plan {
        self.split_keys.entry(parent).or_default().push(key);
        Plan::Child {
            parent,
            key,
            entity,
            role,
        }
    }
    fn keep(&mut self, input: &Topology, id: EntityId) -> Plan {
        self.kept.insert(id);
        Plan::Fixed(input.derivation(id).expect("a derivation").clone())
    }
    /// A new entity generated from `from`.
    fn generated(&mut self, entity: EntityKind, role: Role, ordinal: u32, from: EntityId) -> Plan {
        let d = self.derive(entity, role, ordinal, vec![Parent::Entity(from)]);
        if !self
            .relations
            .iter()
            .any(|r| matches!(r, Relation::Generated { to, .. } if *to == d.id()))
        {
            self.relations.push(Relation::Generated {
                from: d.parents.clone(),
                to: d.id(),
                role,
            });
        }
        Plan::Fixed(d)
    }
}

/// The piece's slot's origin in its own profile: `(boundary, element)`.
fn element_of(topology: &Topology, id: EntityId) -> Option<(usize, ProfileElement)> {
    match topology.derivation(id)?.parents.first() {
        Some(Parent::Profile { boundary, element }) => Some((*boundary as usize, *element)),
        _ => None,
    }
}

fn segment_of(kind: &Kind) -> Segment {
    match kind {
        Kind::Line => Segment::Line,
        Kind::Arc {
            center,
            radius,
            ccw,
        } => Segment::Arc {
            center: *center,
            radius: *radius,
            ccw: *ccw,
        },
        Kind::Spline(span) => Segment::Spline((**span).clone()),
    }
}

fn side_of(sign: i8) -> Side {
    if sign < 0 {
        Side::Below
    } else {
        Side::Above
    }
}

impl Body {
    /// The body's pieces on each side of `plane` (those against its normal
    /// first) and the operation's history (S8e): a sheet's pieces are
    /// sheets of its profile's pieces, a closed wire's its runs on each side
    /// as open wires. A plane parallel to the body, or missing or touching
    /// it, returns the body itself (on the side it lies on; `Below` in the
    /// plane), every entity `Unchanged`. Errors as `Solid::split_by_plane`:
    /// `Degenerate` for a split leaving an edge or a piece thinner than the
    /// resolution or pinching a sheet, `ComputationLimit` for an undecided
    /// comparison; an open wire (a split's piece) split again is
    /// `OutOfDomain`.
    pub fn split_by_plane(
        &self,
        operation: OperationId,
        plane: Frame3,
    ) -> Result<(Vec<(Side, Body)>, History)> {
        self.split_by_plane_in(&Context::new(operation), plane)
    }

    /// [`Body::split_by_plane`] in an operation context.
    pub fn split_by_plane_in(
        &self,
        context: &Context,
        plane: Frame3,
    ) -> Result<(Vec<(Side, Body)>, History)> {
        replayable(context.level)?;
        // The plane in the frame's coordinates at w = 0: a u + b v + d.
        let dot = |u: [f64; 3], v: [f64; 3]| -> R { (0..3).map(|i| q(u[i]) * q(v[i])).sum() };
        let m = plane.normal().to_array();
        let (x, y) = (self.frame.x().to_array(), self.frame.y().to_array());
        let (o, p0) = (self.frame.origin().to_array(), plane.origin().to_array());
        let (a, b) = (dot(m, x), dot(m, y));
        let d: R = (0..3).map(|i| q(m[i]) * (q(o[i]) - q(p0[i]))).sum();
        if a == zero() && b == zero() {
            let side = if d > zero() { Side::Above } else { Side::Below };
            return Ok(self.unchanged(context, side));
        }
        match &self.construction {
            Construction::Face(profile) => {
                let section = Section::new(profile, [a, b, d])?;
                let first = section.pieces[0].side;
                if section.pieces.iter().all(|p| p.side == first) {
                    return Ok(self.unchanged(context, first));
                }
                self.split_sheet(context, profile, &section)
            }
            Construction::Wire {
                boundary,
                tolerance,
            } => {
                let profile = Profile::new(boundary.clone(), Vec::new(), *tolerance)?;
                let section = Section::arrangement(&profile, [a, b, d])?;
                self.split_wire(context, boundary, *tolerance, &section)
            }
            Construction::Path { .. } => Err(Error::OutOfDomain(
                "split_by_plane: an open wire (a split's piece) split again",
            )),
        }
    }

    /// The body itself on `side`, every entity `Unchanged`.
    fn unchanged(&self, context: &Context, side: Side) -> (Vec<(Side, Body)>, History) {
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
        (vec![(side, self.clone())], history)
    }

    /// The sheets of the section's pieces, renamed from this sheet.
    fn split_sheet(
        &self,
        context: &Context,
        profile: &Profile,
        section: &Section,
    ) -> Result<(Vec<(Side, Body)>, History)> {
        let operation = context.operation;
        let input = &self.topology;
        let origins = Origins::new(input, profile.boundaries().collect());
        let face = input
            .ids()
            .find(|(id, _)| input.derivation(*id).map(|d| d.entity) == Some(EntityKind::Face))
            .map(|(id, _)| id)
            .ok_or(Error::InvalidTopology("a sheet without its face"))?;
        let mut naming = Naming {
            operation,
            relations: Vec::new(),
            kept: BTreeSet::new(),
            split_keys: BTreeMap::new(),
        };
        let mut order: Vec<usize> = (0..section.pieces.len()).collect();
        order.sort_by_key(|&i| (section.pieces[i].side, section.pieces[i].first_origin()));
        let mut pieces = Vec::new();
        for (k, &i) in order.iter().enumerate() {
            let piece = &section.pieces[i];
            let built = Body::build(
                operation,
                Construction::Face(Box::new(piece.profile.clone())),
                self.frame,
            )?;
            let mut plans = Vec::new();
            let ids: Vec<(EntityId, Slot)> = built.topology.ids().collect();
            for (id, slot) in ids {
                let d = built.topology.derivation(id).expect("a derivation");
                let plan = match (d.entity, element_of(&built.topology, id)) {
                    (EntityKind::Face, _) => {
                        naming.child(face, (k, 0), EntityKind::Face, Role::Face)
                    }
                    (EntityKind::Edge, Some((pb, ProfileElement::Segment(pj)))) => {
                        match piece.segments[pb][pj as usize] {
                            SegOrigin::Whole(b, j) => {
                                let id = origins.get(
                                    EntityKind::Edge,
                                    Role::Edge,
                                    b,
                                    ProfileElement::Segment(j as u32),
                                )?;
                                naming.keep(input, id)
                            }
                            SegOrigin::Part(b, j, part) => {
                                let from = origins.get(
                                    EntityKind::Edge,
                                    Role::Edge,
                                    b,
                                    ProfileElement::Segment(j as u32),
                                )?;
                                naming.child(from, (k, part), EntityKind::Edge, Role::Edge)
                            }
                            SegOrigin::Chord(c) => naming.generated(
                                EntityKind::Edge,
                                Role::CutEdge,
                                piece_ordinal(k, c),
                                face,
                            ),
                        }
                    }
                    (EntityKind::Vertex, Some((pb, ProfileElement::Vertex(pj)))) => {
                        match piece.vertices[pb][pj as usize] {
                            PointId::Vertex(b, j) => {
                                let id = origins.get(
                                    EntityKind::Vertex,
                                    Role::Vertex,
                                    b,
                                    ProfileElement::Vertex(j as u32),
                                )?;
                                if section.shared.contains(&PointId::Vertex(b, j)) {
                                    naming.child(id, (k, 0), EntityKind::Vertex, Role::Vertex)
                                } else {
                                    naming.keep(input, id)
                                }
                            }
                            PointId::Cross(b, j, which) => {
                                let from = origins.get(
                                    EntityKind::Edge,
                                    Role::Edge,
                                    b,
                                    ProfileElement::Segment(j as u32),
                                )?;
                                naming.generated(
                                    EntityKind::Vertex,
                                    Role::CutVertex,
                                    piece_ordinal(k, j * 2 + which),
                                    from,
                                )
                            }
                        }
                    }
                    _ => return Err(Error::InvalidTopology("a sheet slot of unknown provenance")),
                };
                plans.push((slot, plan));
            }
            pieces.push((piece.side, built, plans));
        }
        self.finish(context, pieces, naming)
    }

    /// The open wires of the boundary's runs on each side, renamed from this
    /// wire.
    fn split_wire(
        &self,
        context: &Context,
        boundary: &Boundary,
        tolerance: Tolerance,
        section: &Section,
    ) -> Result<(Vec<(Side, Body)>, History)> {
        let operation = context.operation;
        let input = &self.topology;
        let edges = &section.edges;
        let n = edges.len();
        // Each piece's side; one along the trace takes the side of the run
        // it continues in stored order.
        let Some(mut current) = edges.iter().rev().find(|e| e.side != 0).map(|e| e.side) else {
            // The whole wire along the trace: in the plane.
            return Ok(self.unchanged(context, Side::Below));
        };
        let mut sides = vec![0i8; n];
        for (i, e) in edges.iter().enumerate() {
            if e.side != 0 {
                current = e.side;
            }
            sides[i] = current;
        }
        if sides.iter().all(|s| *s == sides[0]) {
            return Ok(self.unchanged(context, side_of(sides[0])));
        }
        // Maximal runs, from a change of side round the boundary.
        let start = (0..n)
            .find(|&i| sides[i] != sides[(i + n - 1) % n])
            .expect("a change of side");
        let mut runs: Vec<(i8, Vec<usize>)> = Vec::new();
        for step in 0..n {
            let i = (start + step) % n;
            match runs.last_mut() {
                Some((s, list)) if *s == sides[i] => list.push(i),
                _ => runs.push((sides[i], vec![i])),
            }
        }
        runs.sort_by_key(|(s, list)| (side_of(*s), edges[list[0]].origin));
        let origins = Origins::new(input, vec![boundary]);
        let mut naming = Naming {
            operation,
            relations: Vec::new(),
            kept: BTreeSet::new(),
            split_keys: BTreeMap::new(),
        };
        let mut pieces = Vec::new();
        for (k, (s, list)) in runs.iter().enumerate() {
            let mut ids = vec![edges[list[0]].from];
            ids.extend(list.iter().map(|&i| edges[i].to));
            let points: Vec<Point2> = ids.iter().map(|p| section.positions[p]).collect();
            let segments: Vec<Segment> = list.iter().map(|&i| segment_of(&edges[i].kind)).collect();
            let built = Body::build(
                operation,
                Construction::Path {
                    points,
                    segments,
                    tolerance,
                },
                self.frame,
            )?;
            let last = ids.len() - 1;
            let mut plans = Vec::new();
            let slots: Vec<(EntityId, Slot)> = built.topology.ids().collect();
            for (id, slot) in slots {
                let d = built.topology.derivation(id).expect("a derivation");
                let plan = match (d.entity, element_of(&built.topology, id)) {
                    (EntityKind::Edge, Some((_, ProfileElement::Segment(e)))) => {
                        match edges[list[e as usize]].origin {
                            SegOrigin::Whole(b, j) => {
                                let id = origins.get(
                                    EntityKind::Edge,
                                    Role::Edge,
                                    b,
                                    ProfileElement::Segment(j as u32),
                                )?;
                                naming.keep(input, id)
                            }
                            SegOrigin::Part(b, j, part) => {
                                let from = origins.get(
                                    EntityKind::Edge,
                                    Role::Edge,
                                    b,
                                    ProfileElement::Segment(j as u32),
                                )?;
                                naming.child(from, (k, part), EntityKind::Edge, Role::Edge)
                            }
                            SegOrigin::Chord(_) => {
                                return Err(Error::InvalidTopology("a wire run with a chord"))
                            }
                        }
                    }
                    (EntityKind::Vertex, Some((_, ProfileElement::Vertex(v)))) => {
                        let v = v as usize;
                        match ids[v] {
                            PointId::Vertex(b, j) => {
                                let id = origins.get(
                                    EntityKind::Vertex,
                                    Role::Vertex,
                                    b,
                                    ProfileElement::Vertex(j as u32),
                                )?;
                                // A run's end at a stored vertex: one copy
                                // per run it ends.
                                if v == 0 || v == last {
                                    naming.child(
                                        id,
                                        (k, usize::from(v == last)),
                                        EntityKind::Vertex,
                                        Role::Vertex,
                                    )
                                } else {
                                    naming.keep(input, id)
                                }
                            }
                            PointId::Cross(b, j, which) => {
                                let from = origins.get(
                                    EntityKind::Edge,
                                    Role::Edge,
                                    b,
                                    ProfileElement::Segment(j as u32),
                                )?;
                                naming.generated(
                                    EntityKind::Vertex,
                                    Role::CutVertex,
                                    piece_ordinal(k, j * 2 + which),
                                    from,
                                )
                            }
                        }
                    }
                    _ => return Err(Error::InvalidTopology("a wire slot of unknown provenance")),
                };
                plans.push((slot, plan));
            }
            pieces.push((side_of(*s), built, plans));
        }
        self.finish(context, pieces, naming)
    }

    /// Ordinals by the keys' order, the renamed pieces, their enclosures
    /// raised to their parents' and the history.
    fn finish(
        &self,
        context: &Context,
        pieces: Vec<Unnamed>,
        mut naming: Naming,
    ) -> Result<(Vec<(Side, Body)>, History)> {
        let operation = context.operation;
        let input = &self.topology;
        let body = input.body_id();
        for keys in naming.split_keys.values_mut() {
            keys.sort();
            keys.dedup();
        }
        let before = input.entity_set(self.resolution());
        let mut children: BTreeMap<EntityId, Vec<(u32, EntityId)>> = BTreeMap::new();
        let mut parent_of: BTreeMap<EntityId, EntityId> = BTreeMap::new();
        let mut out = Vec::new();
        for (k, (side, mut built, plans)) in pieces.into_iter().enumerate() {
            let mut derivations = Vec::new();
            for (slot, plan) in plans {
                let d = match plan {
                    Plan::Fixed(d) => d,
                    Plan::Child {
                        parent,
                        key,
                        entity,
                        role,
                    } => {
                        let ordinal = naming.split_keys[&parent]
                            .binary_search(&key)
                            .expect("a key") as u32;
                        let d = naming.derive(entity, role, ordinal, vec![Parent::Entity(parent)]);
                        children.entry(parent).or_default().push((ordinal, d.id()));
                        parent_of.insert(d.id(), parent);
                        d
                    }
                };
                derivations.push((slot, d));
            }
            let body_derivation = naming.derive(
                EntityKind::Body,
                Role::Body,
                k as u32,
                vec![Parent::Entity(body)],
            );
            built.topology = built.topology.renamed(body_derivation, derivations)?;
            // No output's enclosure falls below its parent's (T5).
            built.topology.raise_enclosures(|id| {
                let from = if naming.kept.contains(&id) {
                    id
                } else {
                    *parent_of.get(&id)?
                };
                input.enclosure_bound(from)
            });
            built.operation = operation;
            out.push((side, built));
        }
        let mut relations = std::mem::take(&mut naming.relations);
        for (from, mut into) in children {
            into.sort();
            into.dedup();
            relations.push(Relation::Split {
                from,
                into: into.into_iter().map(|(_, id)| id).collect(),
            });
        }
        let mut modified: BTreeSet<EntityId> = BTreeSet::new();
        for (_, piece) in &out {
            let after = piece.topology.entity_set(piece.resolution());
            for (id, info) in &after.entities {
                if !naming.kept.contains(id) {
                    continue;
                }
                let old = &before.entities[id];
                if old.geometry != info.geometry || old.structure != info.structure {
                    modified.insert(*id);
                }
            }
        }
        for id in &naming.kept {
            relations.push(if modified.contains(id) {
                Relation::Modified { from: *id, to: *id }
            } else {
                Relation::Unchanged { id: *id }
            });
        }
        relations.sort_by_cached_key(Relation::sort_key);
        relations.dedup();
        let outputs: Vec<EntityId> = out.iter().map(|(_, b)| b.topology.body_id()).collect();
        let history = History::new(
            operation,
            OperationKind::PlaneSplit,
            vec![body],
            outputs,
            relations,
            Vec::new(),
        )
        .at_level(context.level);
        if cfg!(debug_assertions) {
            let inputs = [input.entity_set(self.resolution())];
            let outputs: Vec<_> = out
                .iter()
                .map(|(_, b)| b.topology.entity_set(b.resolution()))
                .collect();
            let issues = history::check(&inputs, &outputs, &history);
            assert!(issues.is_empty(), "split history is invalid: {issues:?}");
        }
        Ok((out, history))
    }
}
