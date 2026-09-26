//! Height split and stacked fuse (M3 of `IDENTITY_AND_HISTORY.md`): the first
//! operations whose histories split, merge and delete entities.
//!
//! Both rebuild prisms of the same profile and frame, so every output has the
//! slot layout of its inputs; ids come from the inputs' ids by rule, never from
//! geometry.
use super::{replayable, Solid};
use crate::history::{self, History, Relation};
use crate::identity::{
    AlgorithmLevel, Derivation, EntityId, EntityKind, OperationId, OperationKind, Parent, Role,
};
use crate::topology::{End, Place, Slot, Topology};
use crate::{Error, Result};
use std::collections::BTreeMap;

fn kind_of(slot: Slot) -> EntityKind {
    match slot {
        Slot::Vertex(_) => EntityKind::Vertex,
        Slot::Edge(_) => EntityKind::Edge,
        Slot::Face(_) => EntityKind::Face,
        Slot::Region(_) => EntityKind::Region,
    }
}

fn id_of(t: &Topology, slot: Slot) -> EntityId {
    t.id_of(slot).expect("every prism slot has an id")
}

fn derivation_of(t: &Topology, slot: Slot) -> Derivation {
    t.derivation(id_of(t, slot))
        .expect("every id has a derivation")
        .clone()
}

fn same_layout(a: &Topology, b: &Topology) -> Result<()> {
    if a.layout() != b.layout() || a.layout().is_empty() {
        return Err(Error::InvalidTopology("prism layouts differ"));
    }
    Ok(())
}

impl Solid {
    /// Split by the plane at offset `height` along the frame normal, strictly
    /// between the prism's offsets. The pieces come in axial order (increasing
    /// offset), each extruded in the input's direction.
    ///
    /// History (`OperationKind::HeightSplit`): each piece keeps its own end's
    /// cap, cap edges and cap vertices `Unchanged`; every wall, vertical edge
    /// and the solid region is `Split` into one child per piece (ordinal 0
    /// lower, 1 upper, same role); each piece's cut face is `Generated` from
    /// every wall in profile order, each cut edge from the wall it bounds and
    /// each cut vertex from the vertical edge through it (`CutFace`,
    /// `CutEdge`, `CutVertex`, with the piece ordinal). The input body is
    /// replaced by the two piece bodies.
    pub fn split_at_height(
        &self,
        operation: OperationId,
        height: f64,
    ) -> Result<([Solid; 2], History)> {
        self.split_at_height_at(AlgorithmLevel::CURRENT, operation, height)
    }

    /// [`Solid::split_at_height`] at a recorded algorithm level (H8).
    pub fn split_at_height_at(
        &self,
        level: AlgorithmLevel,
        operation: OperationId,
        height: f64,
    ) -> Result<([Solid; 2], History)> {
        replayable(level)?;
        let (low, high) = (self.start.min(self.end), self.start.max(self.end));
        // Exact comparisons; NaN fails both.
        if !(low < height && height < high) {
            return Err(Error::OutOfDomain("split height outside the prism"));
        }
        let offsets = if self.start < self.end {
            [(self.start, height), (height, self.end)]
        } else {
            [(height, self.end), (self.start, height)]
        };
        let parent = &self.topology;
        let walls: Vec<Parent> = parent
            .layout()
            .iter()
            .filter(|(slot, place)| matches!(slot, Slot::Face(_)) && *place == Place::Swept)
            .map(|(slot, _)| Parent::Entity(id_of(parent, *slot)))
            .collect();
        let body = parent.body_id();
        let mut relations = Vec::new();
        let mut children: BTreeMap<EntityId, Vec<EntityId>> = BTreeMap::new();
        let mut pieces = Vec::new();
        for (k, (start, end)) in offsets.into_iter().enumerate() {
            let mut piece = Self::build(operation, self.profile.clone(), self.frame, start, end)?;
            same_layout(parent, &piece.topology)?;
            let derive = |entity, role, parents| Derivation {
                operation,
                kind: OperationKind::HeightSplit,
                entity,
                role,
                ordinal: k as u32,
                parents,
            };
            let mut named = Vec::new();
            for &(slot, place) in piece.topology.layout() {
                let entity = kind_of(slot);
                let d = match place {
                    Place::Swept => {
                        let whole = id_of(parent, slot);
                        let role = derivation_of(parent, slot).role;
                        let d = derive(entity, role, vec![Parent::Entity(whole)]);
                        children.entry(whole).or_default().push(d.id());
                        d
                    }
                    Place::Low(_) if k == 0 => derivation_of(parent, slot),
                    Place::High(_) if k == 1 => derivation_of(parent, slot),
                    Place::Low(at) | Place::High(at) => {
                        let (role, parents) = match at {
                            End::Cap => (Role::CutFace, walls.clone()),
                            End::Edge(wall) => (
                                Role::CutEdge,
                                vec![Parent::Entity(id_of(parent, Slot::Face(wall)))],
                            ),
                            End::Vertex(vertical) => (
                                Role::CutVertex,
                                vec![Parent::Entity(id_of(parent, Slot::Edge(vertical)))],
                            ),
                        };
                        let d = derive(entity, role, parents);
                        relations.push(Relation::Generated {
                            from: d.parents.clone(),
                            to: d.id(),
                            role,
                        });
                        d
                    }
                };
                named.push((slot, d));
            }
            let body_derivation = derive(EntityKind::Body, Role::Body, vec![Parent::Entity(body)]);
            piece.topology = piece.topology.renamed(body_derivation, named)?;
            piece.operation = operation;
            pieces.push(piece);
        }
        for &(slot, place) in parent.layout() {
            let id = id_of(parent, slot);
            relations.push(match place {
                Place::Swept => Relation::Split {
                    from: id,
                    into: children.remove(&id).unwrap_or_default(),
                },
                _ => Relation::Unchanged { id },
            });
        }
        let [lower, upper]: [Solid; 2] = pieces.try_into().expect("two pieces");
        let mut history = History::new(
            operation,
            OperationKind::HeightSplit,
            vec![body],
            vec![lower.topology.body_id(), upper.topology.body_id()],
            relations,
            Vec::new(),
        );
        history.level = level;
        debug_check(&[self], &[&lower, &upper], &history);
        Ok(([lower, upper], history))
    }

    /// Merge two prisms stacked along their common frame: the same profile
    /// (points and labels) and frame, the same direction, and adjacent
    /// offsets (`self.end == other.start` or `other.end == self.start`,
    /// exactly). The result spans both.
    ///
    /// History (`OperationKind::StackedFuse`): walls, vertical edges and the
    /// solid regions pair up and are `Merged` (parents in body order, ordinal
    /// 0, same role); the shared caps with their edges and vertices are
    /// `Deleted`; the outer ends stay `Unchanged`. Bodies that share an id
    /// cannot be combined (I4).
    pub fn fuse_stacked(&self, other: &Solid, operation: OperationId) -> Result<(Solid, History)> {
        self.fuse_stacked_at(AlgorithmLevel::CURRENT, other, operation)
    }

    /// [`Solid::fuse_stacked`] at a recorded algorithm level (H8).
    pub fn fuse_stacked_at(
        &self,
        level: AlgorithmLevel,
        other: &Solid,
        operation: OperationId,
    ) -> Result<(Solid, History)> {
        replayable(level)?;
        if self.profile != other.profile || self.frame != other.frame {
            return Err(Error::OutOfDomain(
                "stacked fuse needs one profile and one frame",
            ));
        }
        if (self.start < self.end) != (other.start < other.end) {
            return Err(Error::OutOfDomain("stacked fuse needs one direction"));
        }
        let (start, end) = if self.end == other.start {
            (self.start, other.end)
        } else if other.end == self.start {
            (other.start, self.end)
        } else {
            return Err(Error::OutOfDomain("stacked fuse needs adjacent offsets"));
        };
        let (a, b) = (&self.topology, &other.topology);
        if a.body_id() == b.body_id() || a.ids().any(|(id, _)| b.slot_of(id).is_some()) {
            return Err(Error::InvalidTopology("id collision"));
        }
        same_layout(a, b)?;
        let mut fused = Self::build(operation, self.profile.clone(), self.frame, start, end)?;
        same_layout(a, &fused.topology)?;
        let self_is_lower = self.start.min(self.end) < other.start.min(other.end);
        let (lower, upper) = if self_is_lower { (a, b) } else { (b, a) };
        let mut relations = Vec::new();
        let mut named = Vec::new();
        for &(slot, place) in fused.topology.layout() {
            let d = match place {
                Place::Swept => {
                    let from = vec![id_of(a, slot), id_of(b, slot)];
                    let d = Derivation {
                        operation,
                        kind: OperationKind::StackedFuse,
                        entity: kind_of(slot),
                        role: derivation_of(a, slot).role,
                        ordinal: 0,
                        parents: from.iter().copied().map(Parent::Entity).collect(),
                    };
                    relations.push(Relation::Merged { from, into: d.id() });
                    d
                }
                Place::Low(_) => derivation_of(lower, slot),
                Place::High(_) => derivation_of(upper, slot),
            };
            named.push((slot, d));
        }
        for (body, keep_low) in [(lower, true), (upper, false)] {
            for &(slot, place) in body.layout() {
                let id = id_of(body, slot);
                match place {
                    Place::Swept => {}
                    Place::Low(_) if keep_low => relations.push(Relation::Unchanged { id }),
                    Place::High(_) if !keep_low => relations.push(Relation::Unchanged { id }),
                    _ => relations.push(Relation::Deleted { id }),
                }
            }
        }
        let body_derivation = Derivation {
            operation,
            kind: OperationKind::StackedFuse,
            entity: EntityKind::Body,
            role: Role::Body,
            ordinal: 0,
            parents: vec![Parent::Entity(a.body_id()), Parent::Entity(b.body_id())],
        };
        fused.topology = fused.topology.renamed(body_derivation, named)?;
        fused.operation = operation;
        let mut history = History::new(
            operation,
            OperationKind::StackedFuse,
            vec![a.body_id(), b.body_id()],
            vec![fused.topology.body_id()],
            relations,
            Vec::new(),
        );
        history.level = level;
        debug_check(&[self, other], &[&fused], &history);
        Ok((fused, history))
    }
}

/// Every operation checks its history independently in debug builds.
fn debug_check(inputs: &[&Solid], outputs: &[&Solid], history: &History) {
    if cfg!(debug_assertions) {
        let sets = |bodies: &[&Solid]| {
            bodies
                .iter()
                .map(|s| s.topology.entity_set(s.profile.tolerance()))
                .collect::<Vec<_>>()
        };
        let issues = history::check(&sets(inputs), &sets(outputs), history);
        assert!(
            issues.is_empty(),
            "operation history is invalid: {issues:?}"
        );
    }
}
