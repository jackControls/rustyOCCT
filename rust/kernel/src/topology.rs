//! Analytic boundary representation for the supported extrusion builders.
//!
//! Edges are shared by oriented face uses. Each use retains its own parameter
//! curve, including the two distinct parameter curves of a cylindrical seam.
//! Slots (`VertexId`, `EdgeId`, `FaceId`) are dense body-local indices. Every
//! entity also has a value [`EntityId`] derived from how it was made (see
//! `identity.rs`), with maps between slots and ids.
use crate::history::{EntityInfo, EntitySet, Geometry};
use crate::identity::{
    Derivation, EntityId, EntityKind, InputLabel, OperationId, OperationKind, Parent,
    ProfileElement, Role,
};
use crate::profile::BoundaryKind;
use crate::{Error, Frame3, Point2, Point3, Profile, Result, Tolerance, Vec3};
use std::collections::BTreeMap;
use std::f64::consts::TAU;

mod validate;
pub use validate::{EdgeEnd, Entity, Issue, IssueKind};

macro_rules! index_type {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(pub(crate) usize);
        impl $name {
            pub fn new(index: usize) -> Self {
                Self(index)
            }
            pub fn index(self) -> usize {
                self.0
            }
        }
    };
}
index_type!(VertexId);
index_type!(EdgeId);
index_type!(FaceId);

/// An entity's body-local position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Slot {
    Vertex(VertexId),
    Edge(EdgeId),
    Face(FaceId),
}

/// Ids, retained derivations and the slot maps of one body.
#[derive(Debug, Clone, PartialEq)]
struct Identity {
    body: EntityId,
    body_derivation: Derivation,
    vertices: Vec<EntityId>,
    edges: Vec<EntityId>,
    faces: Vec<EntityId>,
    derivations: BTreeMap<EntityId, Derivation>,
    slots: BTreeMap<EntityId, Slot>,
    /// Profile labels, for deriving builder provenance.
    labels: BTreeMap<InputLabel, (u32, ProfileElement)>,
}

impl Identity {
    /// Ids from per-slot derivations; a repeated id is an operation error.
    fn new(
        body: Derivation,
        derivations: Vec<(Slot, Derivation)>,
        labels: BTreeMap<InputLabel, (u32, ProfileElement)>,
    ) -> Result<Self> {
        let mut identity = Self {
            body: body.id(),
            body_derivation: body,
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            derivations: BTreeMap::new(),
            slots: BTreeMap::new(),
            labels,
        };
        let mut ordered = derivations;
        ordered.sort_by_key(|(slot, _)| *slot);
        for (slot, derivation) in ordered {
            let id = derivation.id();
            if id == identity.body || identity.slots.insert(id, slot).is_some() {
                return Err(Error::InvalidTopology("id collision"));
            }
            identity.derivations.insert(id, derivation);
            let list = match slot {
                Slot::Vertex(v) => (&mut identity.vertices, v.0),
                Slot::Edge(e) => (&mut identity.edges, e.0),
                Slot::Face(f) => (&mut identity.faces, f.0),
            };
            if list.0.len() != list.1 {
                return Err(Error::InvalidTopology("slot without a derivation"));
            }
            list.0.push(id);
        }
        Ok(identity)
    }

    /// `External` derivations for caller-supplied parts: one per slot.
    fn external(vertices: usize, edges: usize, faces: usize) -> Result<Self> {
        let d = |entity, ordinal: usize| Derivation {
            operation: OperationId::UNSPECIFIED,
            kind: OperationKind::External,
            entity,
            role: Role::External,
            ordinal: ordinal as u32,
            parents: Vec::new(),
        };
        let slots = (0..vertices)
            .map(|i| (Slot::Vertex(VertexId(i)), d(EntityKind::Vertex, i)))
            .chain((0..edges).map(|i| (Slot::Edge(EdgeId(i)), d(EntityKind::Edge, i))))
            .chain((0..faces).map(|i| (Slot::Face(FaceId(i)), d(EntityKind::Face, i))))
            .collect();
        Self::new(d(EntityKind::Body, 0), slots, BTreeMap::new())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    Reversed,
}

impl Orientation {
    pub fn sign(self) -> f64 {
        if self == Self::Forward {
            1.0
        } else {
            -1.0
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Curve3 {
    LineSegment {
        start: Point3,
        end: Point3,
    },
    /// A full turn from the frame's x axis: the arc with start 0 and sweep TAU.
    Circle {
        frame: Frame3,
        radius: f64,
    },
    /// frame.origin + radius (cos a x + sin a y), a = start + sweep * fraction.
    CircularArc {
        frame: Frame3,
        radius: f64,
        start_angle: f64,
        sweep_angle: f64,
    },
}

impl Curve3 {
    /// Evaluate a normalized edge parameter (0..=1).
    pub fn point(&self, fraction: f64) -> Point3 {
        match self {
            Self::LineSegment { start, end } => *start + (*end - *start) * fraction,
            Self::Circle { frame, radius } => {
                let (sine, cosine) = (fraction * TAU).sin_cos();
                frame.point(Point2::new(radius * cosine, radius * sine), 0.0)
            }
            Self::CircularArc {
                frame,
                radius,
                start_angle,
                sweep_angle,
            } => {
                let (sine, cosine) = (start_angle + sweep_angle * fraction).sin_cos();
                frame.point(Point2::new(radius * cosine, radius * sine), 0.0)
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Curve2 {
    LineSegment {
        start: Point2,
        end: Point2,
    },
    CircularArc {
        center: Point2,
        radius: f64,
        start_angle: f64,
        sweep_angle: f64,
    },
}

impl Curve2 {
    pub fn point(&self, fraction: f64) -> Point2 {
        match self {
            Self::LineSegment { start, end } => Point2::new(
                start.x + (end.x - start.x) * fraction,
                start.y + (end.y - start.y) * fraction,
            ),
            Self::CircularArc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => {
                let (sine, cosine) = (start_angle + sweep_angle * fraction).sin_cos();
                Point2::new(center.x + radius * cosine, center.y + radius * sine)
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Surface {
    Plane(Frame3),
    /// u is an angle in radians, v is axial distance in millimetres.
    Cylinder {
        frame: Frame3,
        radius: f64,
    },
}

impl Surface {
    pub fn point(&self, uv: Point2) -> Point3 {
        match self {
            Self::Plane(frame) => frame.point(uv, 0.0),
            Self::Cylinder { frame, radius } => {
                frame.point(Point2::new(radius * uv.x.cos(), radius * uv.x.sin()), uv.y)
            }
        }
    }
    pub fn normal(&self, uv: Point2) -> Vec3 {
        match self {
            Self::Plane(frame) => frame.normal(),
            Self::Cylinder { frame, .. } => frame.x() * uv.x.cos() + frame.y() * uv.x.sin(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Vertex {
    pub position: Point3,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    pub start: VertexId,
    pub end: VertexId,
    pub curve: Curve3,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Coedge {
    pub edge: EdgeId,
    pub orientation: Orientation,
    /// This edge's curve in the owning face's parameter space, in traversal order.
    pub pcurve: Curve2,
}

/// Builder provenance: useful when translating an application's feature and
/// sketch-entity identities. Boundary 0 is the outer wire; holes start at 1.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceOrigin {
    StartCap,
    EndCap,
    Wall {
        boundary: usize,
        segment: usize,
    },
    /// Supplied through [`Topology::from_parts`] rather than a builder.
    External,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Face {
    pub surface: Surface,
    pub orientation: Orientation,
    /// Outer loop first, then inner loops; traversal follows the oriented face.
    pub loops: Vec<Vec<Coedge>>,
}

impl Face {
    pub fn normal(&self, uv: Point2) -> Vec3 {
        self.surface.normal(uv) * self.orientation.sign()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Topology {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    /// The first shell bounds the solid; later shells bound cavities.
    shells: Vec<Vec<FaceId>>,
    identity: Identity,
}

/// Unvalidated boundary data for [`Topology::from_parts`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TopologyParts {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub faces: Vec<Face>,
    pub shells: Vec<Vec<FaceId>>,
}

impl TopologyParts {
    /// Every issue of the validation contract; empty means valid.
    pub fn check(&self, tolerance: Tolerance) -> Vec<Issue> {
        validate::check(
            &self.vertices,
            &self.edges,
            &self.faces,
            &self.shells,
            tolerance,
        )
    }
}

impl Topology {
    /// Build a topology only if the complete validation contract holds;
    /// otherwise return every issue found.
    pub fn from_parts(
        parts: TopologyParts,
        tolerance: Tolerance,
    ) -> std::result::Result<Self, Vec<Issue>> {
        let issues = parts.check(tolerance);
        if !issues.is_empty() {
            return Err(issues);
        }
        let identity =
            Identity::external(parts.vertices.len(), parts.edges.len(), parts.faces.len())
                .expect("distinct external ordinals give distinct ids");
        Ok(Self {
            vertices: parts.vertices,
            edges: parts.edges,
            faces: parts.faces,
            shells: parts.shells,
            identity,
        })
    }
    /// The body's id; rigid transforms keep it.
    pub fn body_id(&self) -> EntityId {
        self.identity.body
    }
    pub fn body_derivation(&self) -> &Derivation {
        &self.identity.body_derivation
    }
    pub fn id_of(&self, slot: Slot) -> Option<EntityId> {
        match slot {
            Slot::Vertex(v) => self.identity.vertices.get(v.0),
            Slot::Edge(e) => self.identity.edges.get(e.0),
            Slot::Face(f) => self.identity.faces.get(f.0),
        }
        .copied()
    }
    pub fn slot_of(&self, id: EntityId) -> Option<Slot> {
        self.identity.slots.get(&id).copied()
    }
    /// Why the entity has its id.
    pub fn derivation(&self, id: EntityId) -> Option<&Derivation> {
        self.identity.derivations.get(&id)
    }
    /// Every entity id with its slot, in id order.
    pub fn ids(&self) -> impl Iterator<Item = (EntityId, Slot)> + '_ {
        self.identity.slots.iter().map(|(id, slot)| (*id, *slot))
    }
    /// Builder provenance, derived from the face's derivation and the
    /// profile's labels.
    pub fn face_origin(&self, face: FaceId) -> Option<FaceOrigin> {
        let derivation = self.derivation(self.id_of(Slot::Face(face))?)?;
        Some(match derivation.role {
            Role::StartCap => FaceOrigin::StartCap,
            Role::EndCap => FaceOrigin::EndCap,
            Role::Wall => {
                let (boundary, element) = match derivation.parents.first()? {
                    Parent::Label(label) => *self.identity.labels.get(label)?,
                    Parent::Profile { boundary, element } => (*boundary, *element),
                    Parent::Entity(_) => return None,
                };
                let ProfileElement::Segment(segment) = element else {
                    return None;
                };
                FaceOrigin::Wall {
                    boundary: boundary as usize,
                    segment: segment as usize,
                }
            }
            _ => FaceOrigin::External,
        })
    }
    /// Every entity as the history checker sees it.
    pub fn entity_set(&self, tolerance: Tolerance) -> EntitySet {
        let mut entities = BTreeMap::new();
        for (id, slot) in self.ids() {
            let ordinal = self.identity.derivations[&id].ordinal;
            let vid = |v: VertexId| self.identity.vertices[v.0];
            let (kind, geometry, structure) = match slot {
                Slot::Vertex(v) => (
                    EntityKind::Vertex,
                    Geometry::Point(self.vertices[v.0].position),
                    Vec::new(),
                ),
                Slot::Edge(e) => {
                    let edge = &self.edges[e.0];
                    (
                        EntityKind::Edge,
                        Geometry::Curve(edge.curve.clone()),
                        vec![vec![
                            (vid(edge.start), Orientation::Forward),
                            (vid(edge.end), Orientation::Forward),
                        ]],
                    )
                }
                Slot::Face(f) => {
                    let face = &self.faces[f.0];
                    let loops = face
                        .loops
                        .iter()
                        .map(|lp| {
                            lp.iter()
                                .map(|u| (self.identity.edges[u.edge.0], u.orientation))
                                .collect()
                        })
                        .collect();
                    (
                        EntityKind::Face,
                        Geometry::Surface {
                            surface: face.surface.clone(),
                            orientation: face.orientation,
                        },
                        loops,
                    )
                }
            };
            entities.insert(
                id,
                EntityInfo {
                    kind,
                    ordinal,
                    geometry,
                    structure,
                },
            );
        }
        EntitySet {
            body: self.body_id(),
            tolerance,
            entities,
        }
    }
    pub fn shells(&self) -> &[Vec<FaceId>] {
        &self.shells
    }
    /// Every issue of the validation contract; empty means valid.
    pub fn check(&self, tolerance: Tolerance) -> Vec<Issue> {
        validate::check(
            &self.vertices,
            &self.edges,
            &self.faces,
            &self.shells,
            tolerance,
        )
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn vertex(&self, id: VertexId) -> Option<&Vertex> {
        self.vertices.get(id.0)
    }
    pub fn edge(&self, id: EdgeId) -> Option<&Edge> {
        self.edges.get(id.0)
    }
    pub fn face(&self, id: FaceId) -> Option<&Face> {
        self.faces.get(id.0)
    }
    pub fn face_ids(&self) -> impl Iterator<Item = FaceId> {
        (0..self.faces.len()).map(FaceId)
    }
    pub fn edge_ids(&self) -> impl Iterator<Item = EdgeId> {
        (0..self.edges.len()).map(EdgeId)
    }
    pub fn incident_faces(&self, edge: EdgeId) -> Option<Vec<FaceId>> {
        self.edge(edge)?;
        Some(
            self.faces
                .iter()
                .enumerate()
                .filter_map(|(i, face)| {
                    face.loops
                        .iter()
                        .flatten()
                        .any(|coedge| coedge.edge == edge)
                        .then_some(FaceId(i))
                })
                .collect(),
        )
    }
    /// Periodic seams have two uses on the same face, rather than two faces.
    pub fn is_seam(&self, edge: EdgeId) -> Option<bool> {
        self.incident_faces(edge).map(|faces| faces.len() == 1)
    }

    /// Face interiors with holes are not disks: each inner loop subtracts one
    /// from the face's contribution to V - E + F.
    pub fn euler_characteristic(&self) -> i64 {
        self.vertices.len() as i64 - self.edges.len() as i64
            + self
                .faces
                .iter()
                .map(|face| 2 - face.loops.len() as i64)
                .sum::<i64>()
    }

    /// Run the complete validation contract; the error names the first issue.
    /// See [`Topology::check`] for the complete typed report.
    pub fn validate(&self, tolerance: Tolerance) -> Result<()> {
        match self.check(tolerance).first() {
            None => Ok(()),
            Some(issue) => Err(Error::InvalidTopology(issue.kind.name())),
        }
    }

    pub(crate) fn prism(
        profile: &Profile,
        frame: Frame3,
        low: f64,
        high: f64,
        start_is_low: bool,
        operation: OperationId,
    ) -> Result<Self> {
        let tolerance = profile.tolerance();
        // Derivations per slot. "Bottom"/"top" roles name the start/end side.
        let derive = |entity, role, ordinal, parents| Derivation {
            operation,
            kind: OperationKind::Extrude,
            entity,
            role,
            ordinal,
            parents,
        };
        let (low_vertex, high_vertex, low_edge, high_edge, low_cap, high_cap, low_seam) =
            if start_is_low {
                (
                    Role::BottomVertex,
                    Role::TopVertex,
                    Role::BottomEdge,
                    Role::TopEdge,
                    Role::StartCap,
                    Role::EndCap,
                    0,
                )
            } else {
                (
                    Role::TopVertex,
                    Role::BottomVertex,
                    Role::TopEdge,
                    Role::BottomEdge,
                    Role::EndCap,
                    Role::StartCap,
                    1,
                )
            };
        let mut derivations: Vec<(Slot, Derivation)> = Vec::new();
        let mut labels = BTreeMap::new();
        let mut cap_parents = Vec::new();
        for (b, wire) in profile.boundaries().enumerate() {
            let b = b as u32;
            match wire.labels() {
                Some(l) => {
                    cap_parents.push(Parent::Label(l.boundary));
                    labels.insert(l.boundary, (b, ProfileElement::Boundary));
                    for (j, label) in l.segments.iter().enumerate() {
                        labels.insert(*label, (b, ProfileElement::Segment(j as u32)));
                    }
                    for (j, label) in l.vertices.iter().enumerate() {
                        labels.insert(*label, (b, ProfileElement::Vertex(j as u32)));
                    }
                }
                None => cap_parents.push(Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Boundary,
                }),
            }
        }
        derivations.push((
            Slot::Face(FaceId(0)),
            derive(EntityKind::Face, low_cap, 0, cap_parents.clone()),
        ));
        derivations.push((
            Slot::Face(FaceId(1)),
            derive(EntityKind::Face, high_cap, 0, cap_parents),
        ));
        let bottom_frame = Frame3::new(
            frame.point(Point2::default(), low),
            -frame.normal(),
            frame.x(),
            tolerance,
        )?;
        let top_frame = Frame3::new(
            frame.point(Point2::default(), high),
            frame.normal(),
            frame.x(),
            tolerance,
        )?;
        let mut topology = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: vec![
                Face {
                    surface: Surface::Plane(bottom_frame),
                    orientation: Orientation::Forward,
                    loops: Vec::new(),
                },
                Face {
                    surface: Surface::Plane(top_frame),
                    orientation: Orientation::Forward,
                    loops: Vec::new(),
                },
            ],
            shells: Vec::new(),
            identity: Identity::new(
                derive(EntityKind::Body, Role::Body, 0, Vec::new()),
                Vec::new(),
                BTreeMap::new(),
            )?,
        };
        for (boundary, wire) in profile.boundaries().enumerate() {
            let inner = boundary > 0;
            let b = boundary as u32;
            let seg = |j: usize| match wire.labels() {
                Some(l) => Parent::Label(l.segments[j]),
                None => Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Segment(j as u32),
                },
            };
            let vert = |j: usize| match wire.labels() {
                Some(l) => Parent::Label(l.vertices[j]),
                None => Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Vertex(j as u32),
                },
            };
            match &wire.kind {
                BoundaryKind::Polygon(points) => {
                    let bottom = points
                        .iter()
                        .map(|p| topology.add_vertex(frame.point(*p, low)))
                        .collect::<Vec<_>>();
                    let top = points
                        .iter()
                        .map(|p| topology.add_vertex(frame.point(*p, high)))
                        .collect::<Vec<_>>();
                    let count = points.len();
                    let bottom_edges = (0..count)
                        .map(|i| topology.add_line(bottom[i], bottom[(i + 1) % count]))
                        .collect::<Vec<_>>();
                    let top_edges = (0..count)
                        .map(|i| topology.add_line(top[i], top[(i + 1) % count]))
                        .collect::<Vec<_>>();
                    let vertical = (0..count)
                        .map(|i| topology.add_line(bottom[i], top[i]))
                        .collect::<Vec<_>>();
                    for j in 0..count {
                        let (v, e) = (EntityKind::Vertex, EntityKind::Edge);
                        derivations.push((
                            Slot::Vertex(bottom[j]),
                            derive(v, low_vertex, 0, vec![vert(j)]),
                        ));
                        derivations.push((
                            Slot::Vertex(top[j]),
                            derive(v, high_vertex, 0, vec![vert(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(bottom_edges[j]),
                            derive(e, low_edge, 0, vec![seg(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(top_edges[j]),
                            derive(e, high_edge, 0, vec![seg(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(vertical[j]),
                            derive(e, Role::Vertical, 0, vec![vert(j)]),
                        ));
                    }
                    topology.add_cap_loop(0, &bottom_edges, !inner, bottom_frame);
                    topology.add_cap_loop(1, &top_edges, inner, top_frame);
                    for i in 0..count {
                        let j = (i + 1) % count;
                        let (start, end) = if inner { (j, i) } else { (i, j) };
                        let origin = topology.vertices[bottom[start].0].position;
                        let tangent = topology.vertices[bottom[end].0].position - origin;
                        let side_frame =
                            Frame3::new(origin, tangent.cross(frame.normal()), tangent, tolerance)?;
                        let orientation = if inner {
                            Orientation::Reversed
                        } else {
                            Orientation::Forward
                        };
                        let opposite = if inner {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        };
                        let boundary_edges = [
                            (bottom_edges[i], orientation),
                            (vertical[end], Orientation::Forward),
                            (top_edges[i], opposite),
                            (vertical[start], Orientation::Reversed),
                        ];
                        let coedges = boundary_edges
                            .iter()
                            .map(|(edge, orientation)| {
                                topology.plane_use(*edge, *orientation, side_frame)
                            })
                            .collect();
                        derivations.push((
                            Slot::Face(FaceId(topology.faces.len())),
                            derive(EntityKind::Face, Role::Wall, 0, vec![seg(i)]),
                        ));
                        topology.faces.push(Face {
                            surface: Surface::Plane(side_frame),
                            orientation: Orientation::Forward,
                            loops: vec![coedges],
                        });
                    }
                }
                BoundaryKind::Circle { center, radius } => {
                    let cylinder_frame = Frame3::new(
                        frame.point(*center, low),
                        frame.normal(),
                        frame.x(),
                        tolerance,
                    )?;
                    let end_frame = Frame3::new(
                        frame.point(*center, high),
                        frame.normal(),
                        frame.x(),
                        tolerance,
                    )?;
                    let bottom_vertex =
                        topology.add_vertex(cylinder_frame.point(Point2::new(*radius, 0.0), 0.0));
                    let top_vertex =
                        topology.add_vertex(end_frame.point(Point2::new(*radius, 0.0), 0.0));
                    let bottom = topology.add_edge(
                        bottom_vertex,
                        bottom_vertex,
                        Curve3::Circle {
                            frame: cylinder_frame,
                            radius: *radius,
                        },
                    );
                    let top = topology.add_edge(
                        top_vertex,
                        top_vertex,
                        Curve3::Circle {
                            frame: end_frame,
                            radius: *radius,
                        },
                    );
                    let seam = topology.add_line(bottom_vertex, top_vertex);
                    let (v, e) = (EntityKind::Vertex, EntityKind::Edge);
                    derivations.push((
                        Slot::Vertex(bottom_vertex),
                        derive(v, Role::SeamVertex, low_seam, vec![vert(0)]),
                    ));
                    derivations.push((
                        Slot::Vertex(top_vertex),
                        derive(v, Role::SeamVertex, 1 - low_seam, vec![vert(0)]),
                    ));
                    derivations.push((Slot::Edge(bottom), derive(e, low_edge, 0, vec![seg(0)])));
                    derivations.push((Slot::Edge(top), derive(e, high_edge, 0, vec![seg(0)])));
                    derivations.push((Slot::Edge(seam), derive(e, Role::Seam, 0, vec![vert(0)])));
                    derivations.push((
                        Slot::Face(FaceId(topology.faces.len())),
                        derive(EntityKind::Face, Role::Wall, 0, vec![seg(0)]),
                    ));
                    topology.add_cap_loop(0, &[bottom], !inner, bottom_frame);
                    topology.add_cap_loop(1, &[top], inner, top_frame);
                    let (u0, u1, orientation, opposite) = if inner {
                        (TAU, 0.0, Orientation::Reversed, Orientation::Forward)
                    } else {
                        (0.0, TAU, Orientation::Forward, Orientation::Reversed)
                    };
                    let height = high - low;
                    let coedge = |edge, orientation, start, end| Coedge {
                        edge,
                        orientation,
                        pcurve: Curve2::LineSegment { start, end },
                    };
                    let coedges = vec![
                        coedge(
                            bottom,
                            orientation,
                            Point2::new(u0, 0.0),
                            Point2::new(u1, 0.0),
                        ),
                        coedge(
                            seam,
                            Orientation::Forward,
                            Point2::new(u1, 0.0),
                            Point2::new(u1, height),
                        ),
                        coedge(
                            top,
                            opposite,
                            Point2::new(u1, height),
                            Point2::new(u0, height),
                        ),
                        coedge(
                            seam,
                            Orientation::Reversed,
                            Point2::new(u0, height),
                            Point2::new(u0, 0.0),
                        ),
                    ];
                    topology.faces.push(Face {
                        surface: Surface::Cylinder {
                            frame: cylinder_frame,
                            radius: *radius,
                        },
                        orientation: if inner {
                            Orientation::Reversed
                        } else {
                            Orientation::Forward
                        },
                        loops: vec![coedges],
                    });
                }
            }
        }
        topology.shells = vec![topology.face_ids().collect()];
        topology.identity = Identity::new(
            derive(EntityKind::Body, Role::Body, 0, Vec::new()),
            derivations,
            labels,
        )?;
        let identity = &topology.identity;
        if (
            identity.vertices.len(),
            identity.edges.len(),
            identity.faces.len(),
        ) != (
            topology.vertices.len(),
            topology.edges.len(),
            topology.faces.len(),
        ) {
            return Err(Error::InvalidTopology("slot without a derivation"));
        }
        topology.validate(tolerance)?;
        if topology.euler_characteristic() != 2 - 2 * profile.holes().len() as i64 {
            return Err(Error::InvalidTopology("unexpected shell genus"));
        }
        Ok(topology)
    }

    fn add_vertex(&mut self, position: Point3) -> VertexId {
        let id = VertexId(self.vertices.len());
        self.vertices.push(Vertex { position });
        id
    }
    fn add_edge(&mut self, start: VertexId, end: VertexId, curve: Curve3) -> EdgeId {
        let id = EdgeId(self.edges.len());
        self.edges.push(Edge { start, end, curve });
        id
    }
    fn add_line(&mut self, start: VertexId, end: VertexId) -> EdgeId {
        self.add_edge(
            start,
            end,
            Curve3::LineSegment {
                start: self.vertices[start.0].position,
                end: self.vertices[end.0].position,
            },
        )
    }
    fn add_cap_loop(&mut self, face: usize, edges: &[EdgeId], reverse: bool, frame: Frame3) {
        let mut coedges = edges
            .iter()
            .map(|edge| {
                self.plane_use(
                    *edge,
                    if reverse {
                        Orientation::Reversed
                    } else {
                        Orientation::Forward
                    },
                    frame,
                )
            })
            .collect::<Vec<_>>();
        if reverse {
            coedges.reverse();
        }
        self.faces[face].loops.push(coedges);
    }
    fn plane_use(&self, id: EdgeId, orientation: Orientation, frame: Frame3) -> Coedge {
        let edge = &self.edges[id.0];
        let (start, end) = oriented_vertices(edge, orientation);
        let local = |vertex: VertexId| {
            let [x, y, _] = frame.coordinates(self.vertices[vertex.0].position);
            Point2::new(x, y)
        };
        let pcurve = match &edge.curve {
            Curve3::LineSegment { .. } => Curve2::LineSegment {
                start: local(start),
                end: local(end),
            },
            Curve3::Circle {
                frame: circle,
                radius,
            } => {
                let [x, y, _] = frame.coordinates(circle.origin());
                let start = local(start);
                Curve2::CircularArc {
                    center: Point2::new(x, y),
                    radius: *radius,
                    start_angle: (start.y - y).atan2(start.x - x),
                    sweep_angle: TAU
                        * orientation.sign()
                        * circle.normal().dot(frame.normal()).signum(),
                }
            }
            Curve3::CircularArc { .. } => {
                unreachable!("the extrusion builder creates only lines and full circles")
            }
        };
        Coedge {
            edge: id,
            orientation,
            pcurve,
        }
    }
}

fn oriented_vertices(edge: &Edge, orientation: Orientation) -> (VertexId, VertexId) {
    if orientation == Orientation::Forward {
        (edge.start, edge.end)
    } else {
        (edge.end, edge.start)
    }
}
