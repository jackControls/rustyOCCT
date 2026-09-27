//! Bodies without a solid region (S6 of REVIEW_NOTES.md): a planar face
//! from a profile and a wire from a boundary. Their topology is the cell
//! complex's (`TOPOLOGY_MODEL.md`): a face with both sides in the infinite
//! void, or wire edges listed by the void's shell. Their class (D9) is
//! computed from it, never stored.
use crate::history::{self, History, Relation};
use crate::identity::{EntityId, OperationId, OperationKind};
use crate::solid::{replayable, Context};
use crate::topology::{BodyClass, MeasureEnclosure, Topology};
use crate::{Boundary, Frame3, Profile, Result, RigidTransform, Tolerance};

/// How a body was made, kept so rigid motions rebuild it exactly.
#[derive(Debug, Clone, PartialEq)]
enum Construction {
    Face(Box<Profile>),
    Wire {
        boundary: Boundary,
        tolerance: Tolerance,
    },
}

/// An immutable, validated sheet or wire body.
#[derive(Debug, Clone, PartialEq)]
pub struct Body {
    construction: Construction,
    frame: Frame3,
    topology: Topology,
    operation: OperationId,
}

impl Body {
    /// One planar face on the frame's plane bounded by the profile's
    /// boundaries, its normal the frame's (`BRepBuilderAPI_MakeFace` of the
    /// profile's wires on a plane). The face derives from every boundary,
    /// each edge from its segment and each vertex from its point (labels or
    /// stored indices), with the operation kind `MakeFace`; the history is
    /// the construction.
    pub fn face_from_profile_with(
        operation: OperationId,
        profile: Profile,
        frame: Frame3,
    ) -> Result<(Self, History)> {
        Self::face_from_profile_in(&Context::new(operation), profile, frame)
    }

    /// [`Body::face_from_profile_with`] in an operation context (its level;
    /// a construction carries no attributes).
    pub fn face_from_profile_in(
        context: &Context,
        profile: Profile,
        frame: Frame3,
    ) -> Result<(Self, History)> {
        Self::construct(
            context,
            Construction::Face(Box::new(profile)),
            frame,
            OperationKind::MakeFace,
        )
    }

    /// The boundary's edges and vertices on the frame's plane as one closed
    /// wire (`BRepBuilderAPI_MakeWire`), in its stored counter-clockwise
    /// order; a circle is one ring edge. Each edge derives from its segment
    /// and each vertex from its point, with the operation kind `MakeWire`.
    pub fn wire_from_boundary_with(
        operation: OperationId,
        boundary: Boundary,
        frame: Frame3,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        Self::wire_from_boundary_in(&Context::new(operation), boundary, frame, tolerance)
    }

    /// [`Body::wire_from_boundary_with`] in an operation context.
    pub fn wire_from_boundary_in(
        context: &Context,
        boundary: Boundary,
        frame: Frame3,
        tolerance: Tolerance,
    ) -> Result<(Self, History)> {
        Self::construct(
            context,
            Construction::Wire {
                boundary,
                tolerance,
            },
            frame,
            OperationKind::MakeWire,
        )
    }

    fn construct(
        context: &Context,
        construction: Construction,
        frame: Frame3,
        kind: OperationKind,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        replayable(level)?;
        let body = Self::build(operation, construction, frame)?;
        let t = &body.topology;
        let relations = t
            .ids()
            .map(|(id, _)| {
                let d = t.derivation(id).expect("every id has a derivation");
                Relation::Generated {
                    from: d.parents.clone(),
                    to: id,
                    role: d.role,
                }
            })
            .collect();
        let history = History::new(
            operation,
            kind,
            Vec::new(),
            vec![t.body_id()],
            relations,
            Vec::new(),
        )
        .at_level(level);
        body.debug_check(&[], &history);
        Ok((body, history))
    }

    fn build(operation: OperationId, construction: Construction, frame: Frame3) -> Result<Self> {
        let topology = match &construction {
            Construction::Face(profile) => {
                let boundaries: Vec<&Boundary> = profile.boundaries().collect();
                Topology::planar_sheet(&boundaries, frame, profile.tolerance(), true, operation)?
            }
            Construction::Wire {
                boundary,
                tolerance,
            } => Topology::planar_sheet(&[boundary], frame, *tolerance, false, operation)?,
        };
        Ok(Self {
            construction,
            frame,
            topology,
            operation,
        })
    }

    fn debug_check(&self, inputs: &[history::EntitySet], history: &History) {
        if cfg!(debug_assertions) {
            let outputs = [self.topology.entity_set(self.resolution())];
            let issues = history::check(inputs, &outputs, history);
            assert!(
                issues.is_empty(),
                "operation history is invalid: {issues:?}"
            );
        }
    }

    /// Rebuilds the body in the moved frame, keeping the body id and every
    /// entity id; the history reports every entity `Modified` with its id,
    /// and no enclosure falls below the input's (T5).
    pub fn transform_with(
        &self,
        operation: OperationId,
        transform: RigidTransform,
    ) -> Result<(Self, History)> {
        self.transform_in(&Context::new(operation), transform)
    }

    /// [`Body::transform_with`] in an operation context. A face or wire
    /// body carries no attributes, so the context's policies have nothing to
    /// act on.
    pub fn transform_in(
        &self,
        context: &Context,
        transform: RigidTransform,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        replayable(level)?;
        let mut body = Self::build(
            self.operation,
            self.construction.clone(),
            self.frame.transformed(transform, self.resolution())?,
        )?;
        body.topology = body.topology.with_identity_of(&self.topology);
        body.topology
            .raise_enclosures(|id| self.topology.enclosure_bound(id));
        let ids: Vec<EntityId> = self.topology.ids().map(|(id, _)| id).collect();
        let history = History::new(
            operation,
            OperationKind::Transform,
            vec![self.topology.body_id()],
            vec![body.topology.body_id()],
            ids.into_iter()
                .map(|id| Relation::Modified { from: id, to: id })
                .collect(),
            Vec::new(),
        )
        .at_level(level);
        body.debug_check(&[self.topology.entity_set(self.resolution())], &history);
        Ok((body, history))
    }

    /// The profile of a face body; `None` for a wire.
    pub fn profile(&self) -> Option<&Profile> {
        match &self.construction {
            Construction::Face(profile) => Some(profile),
            Construction::Wire { .. } => None,
        }
    }
    /// The boundary of a wire body; `None` for a face.
    pub fn boundary(&self) -> Option<&Boundary> {
        match &self.construction {
            Construction::Face(_) => None,
            Construction::Wire { boundary, .. } => Some(boundary),
        }
    }
    pub fn resolution(&self) -> Tolerance {
        match &self.construction {
            Construction::Face(profile) => profile.tolerance(),
            Construction::Wire { tolerance, .. } => *tolerance,
        }
    }
    pub fn frame(&self) -> Frame3 {
        self.frame
    }
    pub fn topology(&self) -> &Topology {
        &self.topology
    }
    pub fn operation(&self) -> OperationId {
        self.operation
    }
    /// `Sheet` for a face body, `Wire` for a wire body (D9).
    pub fn class(&self) -> BodyClass {
        self.topology.class()
    }
    /// The certified area (face) or length (wire) and centre.
    pub fn measure(&self) -> Option<MeasureEnclosure> {
        self.topology.measure_enclosure()
    }
}
