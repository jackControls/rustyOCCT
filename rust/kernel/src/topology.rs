//! Analytic boundary representation: a cellular partition of space
//! (`TOPOLOGY_MODEL.md`, decisions D1–D12).
//!
//! A body is a set of regions, shells, faces, loops, fins, edges and
//! vertices. Region 0 is the infinite void; bounded regions are solid or void.
//! A shell is one connected boundary component of one region and lists face
//! sides. A face has two sides and stores the shell on each. A loop is an
//! ordered cycle of fins, with a winding number per periodic direction, or a
//! single vertex. A fin is one loop's oriented use of an edge and carries its
//! own pcurve. Each edge stores its fins in radial order. There are no seams:
//! closed curves without vertices are ring edges, and loops on a cylinder
//! close modulo the period.
//!
//! Slots are dense body-local indices. Vertices, edges, faces and bounded
//! regions also have value ids ([`EntityId`]) derived from how they were made
//! (see `identity.rs`); shells, loops, fins and the infinite void are
//! structure.
use crate::attributes::{Attribute, AttributeMap};
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
index_type!(FinId);
index_type!(LoopId);
index_type!(FaceId);
index_type!(ShellId);
index_type!(RegionId);

/// An entity with identity, by its body-local position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Slot {
    Vertex(VertexId),
    Edge(EdgeId),
    Face(FaceId),
    /// A bounded region; region 0, the infinite void, has no id.
    Region(RegionId),
}

/// Ids, retained derivations and the slot maps of one body.
#[derive(Debug, Clone, PartialEq)]
struct Identity {
    body: EntityId,
    body_derivation: Derivation,
    vertices: Vec<EntityId>,
    edges: Vec<EntityId>,
    faces: Vec<EntityId>,
    /// Region r at index r - 1.
    regions: Vec<EntityId>,
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
            regions: Vec::new(),
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
                Slot::Region(r) => (&mut identity.regions, r.0.wrapping_sub(1)),
            };
            if list.0.len() != list.1 {
                return Err(Error::InvalidTopology("slot without a derivation"));
            }
            list.0.push(id);
        }
        Ok(identity)
    }

    /// `External` derivations for caller-supplied parts: one per slot.
    fn external(vertices: usize, edges: usize, faces: usize, regions: usize) -> Result<Self> {
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
            .chain((1..regions).map(|i| (Slot::Region(RegionId(i)), d(EntityKind::Region, i))))
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
    /// `O + (radius + v sin a)(cos u x + sin u y) + v cos a n` with `a` the
    /// half angle, `0 < |a| < pi/2`: v runs along the generatrix, as in
    /// OCCT's `Geom_ConicalSurface`, and the apex is at `v = -radius / sin a`.
    Cone {
        frame: Frame3,
        radius: f64,
        half_angle: f64,
    },
}

impl Surface {
    pub fn point(&self, uv: Point2) -> Point3 {
        match self {
            Self::Plane(frame) => frame.point(uv, 0.0),
            Self::Cylinder { frame, radius } => {
                frame.point(Point2::new(radius * uv.x.cos(), radius * uv.x.sin()), uv.y)
            }
            Self::Cone {
                frame,
                radius,
                half_angle,
            } => {
                let rho = radius + uv.y * half_angle.sin();
                frame.point(
                    Point2::new(rho * uv.x.cos(), rho * uv.x.sin()),
                    uv.y * half_angle.cos(),
                )
            }
        }
    }
    /// The unit normal of the parametrization, `S_u x S_v` normalized; on a
    /// cone, for the nappe where `radius + v sin a` is positive.
    pub fn normal(&self, uv: Point2) -> Vec3 {
        match self {
            Self::Plane(frame) => frame.normal(),
            Self::Cylinder { frame, .. } => frame.x() * uv.x.cos() + frame.y() * uv.x.sin(),
            Self::Cone {
                frame, half_angle, ..
            } => {
                let radial = frame.x() * uv.x.cos() + frame.y() * uv.x.sin();
                radial * half_angle.cos() - frame.normal() * half_angle.sin()
            }
        }
    }
    /// Periodic in u (an angle): cylinders and cones.
    pub fn is_periodic(&self) -> bool {
        !matches!(self, Self::Plane(_))
    }
}

/// Where an enclosure came from (Contract 5 of `IDENTITY_AND_HISTORY.md`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Provenance {
    /// A certified upper bound computed from the stored geometry, or carried
    /// from an operation's inputs.
    Computed,
    /// A tolerance another system stored (an OCCT `.brep` file); the checker
    /// still verifies it.
    Imported,
}

/// An upper bound on a representation gap, in the body's length unit: for a
/// vertex, its distance to the ends of its edges' curves (and to the surface
/// of a vertex loop); for a fin, the distance between the edge curve and the
/// pcurve's image; for a face, the gaps between consecutive fins in its
/// parameter space (angles scaled by the radius). It never exceeds the
/// body's resolution, and `Topology::check` verifies it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Enclosure {
    pub bound: f64,
    pub provenance: Provenance,
}

impl Enclosure {
    pub fn computed(bound: f64) -> Self {
        Self {
            bound,
            provenance: Provenance::Computed,
        }
    }
    pub fn imported(bound: f64) -> Self {
        Self {
            bound,
            provenance: Provenance::Imported,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Vertex {
    pub position: Point3,
    pub enclosure: Option<Enclosure>,
}

/// A curve bounded by vertices, or a ring edge (a closed curve with neither).
#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    pub start: Option<VertexId>,
    pub end: Option<VertexId>,
    pub curve: Curve3,
    /// Every fin of this edge, in radial order about the curve tangent.
    pub fins: Vec<FinId>,
}

impl Edge {
    pub fn is_ring(&self) -> bool {
        self.start.is_none() && self.end.is_none()
    }
}

/// One loop's oriented use of an edge.
#[derive(Debug, Clone, PartialEq)]
pub struct Fin {
    pub edge: EdgeId,
    /// Against the edge curve.
    pub sense: Orientation,
    /// This edge's curve in the owning face's parameter space, in traversal
    /// order. On a periodic surface it lives in the universal cover.
    pub pcurve: Curve2,
    pub enclosure: Option<Enclosure>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Loop {
    /// An ordered cycle of fins. It closes modulo the surface period with
    /// `winding[d]` turns in periodic direction d (u, v); both are zero on a
    /// plane, and a cylinder is periodic in u only.
    Edges { fins: Vec<FinId>, winding: [i32; 2] },
    /// A pole or an immersed vertex.
    Vertex(VertexId),
}

/// Which side of a face: the oriented normal points from the front side's
/// region into the back side's region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Side {
    Front,
    Back,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RegionKind {
    Solid,
    Void,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Shell {
    pub region: RegionId,
    pub sides: Vec<(FaceId, Side)>,
    /// Edges with no fins that belong to this shell.
    pub wire_edges: Vec<EdgeId>,
    pub acorn_vertices: Vec<VertexId>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    pub kind: RegionKind,
    /// The first shell of a bounded region is its outer boundary.
    pub shells: Vec<ShellId>,
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
    /// Against the surface normal.
    pub sense: Orientation,
    /// On a plane the outer loop comes first; traversal follows the oriented face.
    pub loops: Vec<LoopId>,
    pub front: ShellId,
    pub back: ShellId,
    pub enclosure: Option<Enclosure>,
}

impl Face {
    pub fn normal(&self, uv: Point2) -> Vec3 {
        self.surface.normal(uv) * self.sense.sign()
    }
}

/// Certified enclosures `[lo, hi]` of a body's mass properties at unit
/// density (REVIEW_NOTES.md U2): volume, surface area, centre of gravity and
/// the inertia tensor about it, in world axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassEnclosure {
    pub volume: [f64; 2],
    pub surface_area: [f64; 2],
    pub centroid: [[f64; 2]; 3],
    pub inertia: [[[f64; 2]; 3]; 3],
}

impl MassEnclosure {
    /// The midpoints of every enclosure.
    pub fn midpoints(&self) -> crate::MassProperties {
        let mid = |x: [f64; 2]| 0.5 * x[0] + 0.5 * x[1];
        crate::MassProperties {
            volume: mid(self.volume),
            surface_area: mid(self.surface_area),
            centroid: Point3::new(
                mid(self.centroid[0]),
                mid(self.centroid[1]),
                mid(self.centroid[2]),
            ),
            inertia: std::array::from_fn(|a| std::array::from_fn(|b| mid(self.inertia[a][b]))),
        }
    }
    /// The largest half width, relative to each value's magnitude (at least
    /// 1 for the centroid, in the body's length unit).
    pub fn relative_width(&self) -> f64 {
        let rel = |x: [f64; 2], scale: f64| 0.5 * (x[1] - x[0]) / scale.max(f64::MIN_POSITIVE);
        let mut worst = rel(self.volume, self.volume[1].abs())
            .max(rel(self.surface_area, self.surface_area[1].abs()));
        for c in &self.centroid {
            worst = worst.max(rel(*c, c[1].abs().max(1.0)));
        }
        worst
    }
}

/// The counts OCCT's `nbshapes` reports for the same body (a synthesis, see
/// [`Topology::occt_counts`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OcctCounts {
    pub vertices: usize,
    pub edges: usize,
    pub wires: usize,
    pub faces: usize,
    pub shells: usize,
    pub solids: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Topology {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    fins: Vec<Fin>,
    loops: Vec<Loop>,
    faces: Vec<Face>,
    shells: Vec<Shell>,
    /// regions[0] is the infinite void.
    regions: Vec<Region>,
    identity: Identity,
    /// Where each slot of a builder-made prism sits; empty for `from_parts`.
    layout: Vec<(Slot, Place)>,
    /// Opaque attributes by entity id (contract 4), each list sorted by key
    /// with one value per key.
    attributes: AttributeMap,
}

/// Where a slot of a prism sits along its axis (M3): on the lower or the
/// upper end, or spanning the prism.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Place {
    Low(End),
    High(End),
    /// A wall, a vertical edge or the solid region.
    Swept,
}

/// An entity on one end of a prism, with the swept entity it bounds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum End {
    Cap,
    /// A cap edge and the wall it bounds.
    Edge(FaceId),
    /// A cap vertex and the vertical edge through it.
    Vertex(EdgeId),
}

/// Unvalidated boundary data for [`Topology::from_parts`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TopologyParts {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub fins: Vec<Fin>,
    pub loops: Vec<Loop>,
    pub faces: Vec<Face>,
    pub shells: Vec<Shell>,
    pub regions: Vec<Region>,
}

impl TopologyParts {
    fn view(&self) -> validate::View<'_> {
        validate::View {
            vertices: &self.vertices,
            edges: &self.edges,
            fins: &self.fins,
            loops: &self.loops,
            faces: &self.faces,
            shells: &self.shells,
            regions: &self.regions,
        }
    }
    /// Every issue of the validation contract; empty means valid.
    pub fn check(&self, tolerance: Tolerance) -> Vec<Issue> {
        validate::check(&self.view(), tolerance)
    }
    /// The same parts with a measured, computed enclosure on every vertex,
    /// fin and face that has none: a certified upper bound of its gaps from
    /// the stored geometry. Entities whose gaps have no finite certified bound
    /// keep none, and [`TopologyParts::check`] reports them.
    pub fn with_measured_enclosures(mut self) -> Self {
        let m = validate::measure(&self.view());
        fill(&mut self.vertices, &m.vertices, |v| &mut v.enclosure);
        fill(&mut self.fins, &m.fins, |f| &mut f.enclosure);
        fill(&mut self.faces, &m.faces, |f| &mut f.enclosure);
        self
    }
}

fn fill<T>(
    items: &mut [T],
    bounds: &[Option<f64>],
    slot: impl Fn(&mut T) -> &mut Option<Enclosure>,
) {
    for (item, bound) in items.iter_mut().zip(bounds) {
        let enclosure = slot(item);
        if enclosure.is_none() {
            *enclosure = bound.map(Enclosure::computed);
        }
    }
}

impl Topology {
    fn view(&self) -> validate::View<'_> {
        validate::View {
            vertices: &self.vertices,
            edges: &self.edges,
            fins: &self.fins,
            loops: &self.loops,
            faces: &self.faces,
            shells: &self.shells,
            regions: &self.regions,
        }
    }
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
        let identity = Identity::external(
            parts.vertices.len(),
            parts.edges.len(),
            parts.faces.len(),
            parts.regions.len(),
        )
        .expect("distinct external ordinals give distinct ids");
        Ok(Self {
            vertices: parts.vertices,
            edges: parts.edges,
            fins: parts.fins,
            loops: parts.loops,
            faces: parts.faces,
            shells: parts.shells,
            regions: parts.regions,
            identity,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
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
            Slot::Region(r) => self.identity.regions.get(r.0.wrapping_sub(1)),
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
            Role::External => FaceOrigin::External,
            _ => return None,
        })
    }
    /// Every entity with identity as the history checker sees it.
    pub fn entity_set(&self, tolerance: Tolerance) -> EntitySet {
        let mut entities = BTreeMap::new();
        let vid = |v: Option<VertexId>| v.map(|v| self.identity.vertices[v.0]);
        let fid = |f: FaceId| self.identity.faces[f.0];
        for (id, slot) in self.ids() {
            let ordinal = self.identity.derivations[&id].ordinal;
            let (kind, geometry, structure) = match slot {
                Slot::Vertex(v) => (
                    EntityKind::Vertex,
                    Geometry::Point(self.vertices[v.0].position),
                    Vec::new(),
                ),
                Slot::Edge(e) => {
                    let edge = &self.edges[e.0];
                    let ends = [edge.start, edge.end]
                        .into_iter()
                        .filter_map(vid)
                        .map(|v| (v, Orientation::Forward))
                        .collect();
                    (
                        EntityKind::Edge,
                        Geometry::Curve(edge.curve.clone()),
                        vec![ends],
                    )
                }
                Slot::Face(f) => {
                    let face = &self.faces[f.0];
                    let loops = face
                        .loops
                        .iter()
                        .map(|l| match &self.loops[l.0] {
                            Loop::Edges { fins, .. } => fins
                                .iter()
                                .map(|k| {
                                    let fin = &self.fins[k.0];
                                    (self.identity.edges[fin.edge.0], fin.sense)
                                })
                                .collect(),
                            Loop::Vertex(v) => {
                                vec![(self.identity.vertices[v.0], Orientation::Forward)]
                            }
                        })
                        .collect();
                    (
                        EntityKind::Face,
                        Geometry::Surface {
                            surface: face.surface.clone(),
                            orientation: face.sense,
                        },
                        loops,
                    )
                }
                Slot::Region(r) => {
                    let region = &self.regions[r.0];
                    let shells = region
                        .shells
                        .iter()
                        .map(|s| {
                            self.shells[s.0]
                                .sides
                                .iter()
                                .map(|(f, side)| {
                                    let o = if *side == Side::Front {
                                        Orientation::Forward
                                    } else {
                                        Orientation::Reversed
                                    };
                                    (fid(*f), o)
                                })
                                .collect()
                        })
                        .collect();
                    (EntityKind::Region, Geometry::Region(region.kind), shells)
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
    /// Every issue of the validation contract; empty means valid.
    pub fn check(&self, tolerance: Tolerance) -> Vec<Issue> {
        validate::check(&self.view(), tolerance)
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }
    pub fn fins(&self) -> &[Fin] {
        &self.fins
    }
    pub fn loops(&self) -> &[Loop] {
        &self.loops
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn shells(&self) -> &[Shell] {
        &self.shells
    }
    pub fn regions(&self) -> &[Region] {
        &self.regions
    }
    pub fn vertex(&self, id: VertexId) -> Option<&Vertex> {
        self.vertices.get(id.0)
    }
    pub fn edge(&self, id: EdgeId) -> Option<&Edge> {
        self.edges.get(id.0)
    }
    pub fn fin(&self, id: FinId) -> Option<&Fin> {
        self.fins.get(id.0)
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
    /// The fins of one face, loop by loop.
    pub fn face_fins(&self, face: FaceId) -> Vec<Vec<&Fin>> {
        self.faces[face.0]
            .loops
            .iter()
            .map(|l| match &self.loops[l.0] {
                Loop::Edges { fins, .. } => fins.iter().map(|k| &self.fins[k.0]).collect(),
                Loop::Vertex(_) => Vec::new(),
            })
            .collect()
    }
    /// The faces using an edge, one per fin, in the edge's radial order.
    pub fn incident_faces(&self, edge: EdgeId) -> Option<Vec<FaceId>> {
        let edge = self.edge(edge)?;
        let mut owner = BTreeMap::new();
        for (f, face) in self.faces.iter().enumerate() {
            for l in &face.loops {
                if let Loop::Edges { fins, .. } = &self.loops[l.0] {
                    for k in fins {
                        owner.insert(*k, FaceId(f));
                    }
                }
            }
        }
        Some(
            edge.fins
                .iter()
                .filter_map(|k| owner.get(k).copied())
                .collect(),
        )
    }

    /// V - E + sum over faces (2 - loops). A ring edge is a closed cell with
    /// no vertex, contributing nothing, and a vertex loop adds its vertex and
    /// a loop.
    pub fn euler_characteristic(&self) -> i64 {
        let edges = self.edges.iter().filter(|e| !e.is_ring()).count() as i64;
        self.vertices.len() as i64 - edges
            + self
                .faces
                .iter()
                .map(|face| 2 - face.loops.len() as i64)
                .sum::<i64>()
    }

    /// The counts OCCT's `nbshapes` would report for this body, synthesized
    /// by rule (TOPOLOGY_MODEL.md, T2): each face with loops winding a
    /// periodic direction gets one seam edge per such direction, and each
    /// ring edge those loops use gets one seam vertex; a cone face's pole
    /// gets one degenerated edge (its vertex is the pole's); a face's wound
    /// loops form one wire and every other edge loop its own wire; shells
    /// and solids are those of solid regions.
    pub fn occt_counts(&self) -> OcctCounts {
        let mut seams = 0;
        let mut degenerate = 0;
        let mut seam_vertices = std::collections::BTreeSet::new();
        let mut wires = 0;
        for face in &self.faces {
            let mut wound = [false, false];
            for l in &face.loops {
                if let Loop::Edges { fins, winding } = &self.loops[l.0] {
                    if winding == &[0, 0] {
                        wires += 1;
                        continue;
                    }
                    for (d, w) in winding.iter().enumerate() {
                        wound[d] |= *w != 0;
                    }
                    for k in fins {
                        let e = self.fins[k.0].edge;
                        if self.edges[e.0].is_ring() {
                            seam_vertices.insert(e);
                        }
                    }
                }
            }
            seams += wound.iter().filter(|w| **w).count();
            wires += usize::from(wound.iter().any(|w| *w));
            // OCCT closes a cone's band at the apex with a degenerated edge.
            if validate::pole_position(face, &self.loops).is_some() {
                degenerate += 1;
            }
        }
        let solid: Vec<&Region> = self
            .regions
            .iter()
            .filter(|r| r.kind == RegionKind::Solid)
            .collect();
        OcctCounts {
            vertices: self.vertices.len() + seam_vertices.len(),
            edges: self.edges.len() + seams + degenerate,
            wires,
            faces: self.faces.len(),
            shells: solid.iter().map(|r| r.shells.len()).sum(),
            solids: solid.len(),
        }
    }

    /// Run the complete validation contract; the error names the first issue.
    /// See [`Topology::check`] for the complete typed report.
    pub fn validate(&self, tolerance: Tolerance) -> Result<()> {
        match self.check(tolerance).first() {
            None => Ok(()),
            Some(issue) => Err(Error::InvalidTopology(issue.kind.name())),
        }
    }

    /// A reference point for integration: the first vertex, or the first
    /// face's frame origin.
    fn reference_point(&self) -> [f64; 3] {
        if let Some(v) = self.vertices.first() {
            return v.position.to_array();
        }
        match self.faces.first().map(|f| &f.surface) {
            Some(
                Surface::Plane(f)
                | Surface::Cylinder { frame: f, .. }
                | Surface::Cone { frame: f, .. },
            ) => f.origin().to_array(),
            None => [0.0; 3],
        }
    }
    /// Certified mass properties of the body's solid regions (general:
    /// planes, cylinders and cones; REVIEW_NOTES.md U2). `None` when a face
    /// cannot be integrated (an arc pcurve on a cylinder or cone) or a
    /// certified division fails.
    pub fn mass_enclosure(&self) -> Option<MassEnclosure> {
        let e = validate::mass(&self.view(), self.reference_point())?;
        Some(MassEnclosure {
            volume: e.volume,
            surface_area: e.area,
            centroid: e.centroid,
            inertia: e.inertia,
        })
    }
    /// A face's area and centre of gravity, from certified enclosures.
    pub fn face_area_and_centre(&self, face: FaceId) -> Option<(f64, Point3)> {
        let (area, centre) = validate::face_mass(&self.view(), face.0, self.reference_point())?;
        let mid = |x: [f64; 2]| 0.5 * x[0] + 0.5 * x[1];
        Some((
            mid(area),
            Point3::new(mid(centre[0]), mid(centre[1]), mid(centre[2])),
        ))
    }

    /// Every entity's attributes, by id; each list is sorted by key.
    pub fn attributes(&self) -> &AttributeMap {
        &self.attributes
    }
    /// The same topology with `attribute` on entity `id`, replacing any
    /// value under the same key.
    pub fn with_attribute(mut self, id: EntityId, attribute: Attribute) -> Result<Self> {
        if self.slot_of(id).is_none() {
            return Err(Error::OutOfDomain(
                "attribute on an id the body does not have",
            ));
        }
        let list = self.attributes.entry(id).or_default();
        list.retain(|a| a.key != attribute.key);
        list.push(attribute);
        list.sort();
        Ok(self)
    }
    /// An entity's enclosure bound: a vertex's or face's own, an edge's
    /// largest over its fins; `None` for regions or entities without one.
    pub(crate) fn enclosure_bound(&self, id: EntityId) -> Option<f64> {
        match self.slot_of(id)? {
            Slot::Vertex(v) => self.vertices[v.0].enclosure.map(|e| e.bound),
            Slot::Face(f) => self.faces[f.0].enclosure.map(|e| e.bound),
            Slot::Edge(e) => self.edges[e.0]
                .fins
                .iter()
                .filter_map(|k| self.fins[k.0].enclosure.map(|e| e.bound))
                .reduce(f64::max),
            Slot::Region(_) => None,
        }
    }
    /// Raise every vertex, face and fin enclosure to `parents(id)` of its
    /// entity (an edge's for a fin) where that is larger: an output's bound
    /// never falls below its inputs' (Contract 5, T5).
    pub(crate) fn raise_enclosures(&mut self, parents: impl Fn(EntityId) -> Option<f64>) {
        let raise = |e: &mut Option<Enclosure>, b: Option<f64>| {
            if let (Some(e), Some(b)) = (e.as_mut(), b) {
                e.bound = e.bound.max(b);
            }
        };
        let slots: Vec<(EntityId, Slot)> = self.ids().collect();
        for (id, slot) in slots {
            let b = parents(id);
            match slot {
                Slot::Vertex(v) => raise(&mut self.vertices[v.0].enclosure, b),
                Slot::Face(f) => raise(&mut self.faces[f.0].enclosure, b),
                Slot::Edge(e) => {
                    for k in self.edges[e.0].fins.clone() {
                        raise(&mut self.fins[k.0].enclosure, b);
                    }
                }
                Slot::Region(_) => {}
            }
        }
    }
    pub(crate) fn set_attributes(&mut self, attributes: AttributeMap) {
        self.attributes = attributes;
    }
    /// Where each slot of a builder-made prism sits, in slot order.
    pub(crate) fn layout(&self) -> &[(Slot, Place)] {
        &self.layout
    }
    /// The same structure under new ids: every slot gets the given
    /// derivation; a repeated id or an unnamed slot is an error.
    pub(crate) fn renamed(
        mut self,
        body: Derivation,
        derivations: Vec<(Slot, Derivation)>,
    ) -> Result<Self> {
        let named = derivations.len();
        let labels = std::mem::take(&mut self.identity.labels);
        self.identity = Identity::new(body, derivations, labels)?;
        self.attributes.clear();
        let slots =
            self.vertices.len() + self.edges.len() + self.faces.len() + self.regions.len() - 1;
        if named != slots || self.identity.slots.len() != slots {
            return Err(Error::InvalidTopology("slot without a derivation"));
        }
        Ok(self)
    }
    /// The same structure with the ids of `other`, a rigid copy of it.
    pub(crate) fn with_identity_of(mut self, other: &Topology) -> Self {
        debug_assert_eq!(self.layout, other.layout);
        self.identity = other.identity.clone();
        self.attributes.clear();
        self
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
        let (low_vertex, high_vertex, low_edge, high_edge, low_cap, high_cap) = if start_is_low {
            (
                Role::BottomVertex,
                Role::TopVertex,
                Role::BottomEdge,
                Role::TopEdge,
                Role::StartCap,
                Role::EndCap,
            )
        } else {
            (
                Role::TopVertex,
                Role::BottomVertex,
                Role::TopEdge,
                Role::BottomEdge,
                Role::EndCap,
                Role::StartCap,
            )
        };
        let mut derivations: Vec<(Slot, Derivation)> = Vec::new();
        let mut layout: Vec<(Slot, Place)> = vec![
            (Slot::Face(FaceId(0)), Place::Low(End::Cap)),
            (Slot::Face(FaceId(1)), Place::High(End::Cap)),
            (Slot::Region(RegionId(1)), Place::Swept),
        ];
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
            derive(EntityKind::Face, high_cap, 0, cap_parents.clone()),
        ));
        derivations.push((
            Slot::Region(RegionId(1)),
            derive(EntityKind::Region, Role::Region, 0, cap_parents),
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
        // One solid region bounded by shell 0 (every front side) inside the
        // infinite void bounded by shell 1 (every back side).
        let cap = |surface| Face {
            surface,
            sense: Orientation::Forward,
            loops: Vec::new(),
            front: ShellId(0),
            back: ShellId(1),
            enclosure: None,
        };
        let mut topology = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            fins: Vec::new(),
            loops: Vec::new(),
            faces: vec![
                cap(Surface::Plane(bottom_frame)),
                cap(Surface::Plane(top_frame)),
            ],
            shells: Vec::new(),
            regions: Vec::new(),
            identity: Identity::new(
                derive(EntityKind::Body, Role::Body, 0, Vec::new()),
                Vec::new(),
                BTreeMap::new(),
            )?,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
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
            let (forward, reversed) = if inner {
                (Orientation::Reversed, Orientation::Forward)
            } else {
                (Orientation::Forward, Orientation::Reversed)
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
                    // The walls of this boundary follow the faces made so far.
                    let wall = |j: usize| FaceId(topology.faces.len() + j);
                    for j in 0..count {
                        layout.extend([
                            (
                                Slot::Vertex(bottom[j]),
                                Place::Low(End::Vertex(vertical[j])),
                            ),
                            (Slot::Vertex(top[j]), Place::High(End::Vertex(vertical[j]))),
                            (Slot::Edge(bottom_edges[j]), Place::Low(End::Edge(wall(j)))),
                            (Slot::Edge(top_edges[j]), Place::High(End::Edge(wall(j)))),
                            (Slot::Edge(vertical[j]), Place::Swept),
                            (Slot::Face(wall(j)), Place::Swept),
                        ]);
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
                        let fins = [
                            (bottom_edges[i], forward),
                            (vertical[end], Orientation::Forward),
                            (top_edges[i], reversed),
                            (vertical[start], Orientation::Reversed),
                        ]
                        .iter()
                        .map(|(edge, sense)| topology.plane_fin(*edge, *sense, side_frame))
                        .collect();
                        derivations.push((
                            Slot::Face(FaceId(topology.faces.len())),
                            derive(EntityKind::Face, Role::Wall, 0, vec![seg(i)]),
                        ));
                        let l = topology.add_loop(fins, [0, 0]);
                        topology.faces.push(Face {
                            surface: Surface::Plane(side_frame),
                            sense: Orientation::Forward,
                            loops: vec![l],
                            front: ShellId(0),
                            back: ShellId(1),
                            enclosure: None,
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
                    // Two ring edges and the wall; no seam, no vertex.
                    let bottom = topology.add_ring(Curve3::Circle {
                        frame: cylinder_frame,
                        radius: *radius,
                    });
                    let top = topology.add_ring(Curve3::Circle {
                        frame: end_frame,
                        radius: *radius,
                    });
                    let e = EntityKind::Edge;
                    let wall = FaceId(topology.faces.len());
                    layout.extend([
                        (Slot::Edge(bottom), Place::Low(End::Edge(wall))),
                        (Slot::Edge(top), Place::High(End::Edge(wall))),
                        (Slot::Face(wall), Place::Swept),
                    ]);
                    derivations.push((Slot::Edge(bottom), derive(e, low_edge, 0, vec![seg(0)])));
                    derivations.push((Slot::Edge(top), derive(e, high_edge, 0, vec![seg(0)])));
                    derivations.push((
                        Slot::Face(FaceId(topology.faces.len())),
                        derive(EntityKind::Face, Role::Wall, 0, vec![seg(0)]),
                    ));
                    topology.add_cap_loop(0, &[bottom], !inner, bottom_frame);
                    topology.add_cap_loop(1, &[top], inner, top_frame);
                    // The wall's loops wind once around the axis, in opposite
                    // directions, on the universal cover of the cylinder.
                    let (u0, u1, turns) = if inner { (TAU, 0.0, -1) } else { (0.0, TAU, 1) };
                    let height = high - low;
                    let fin = |edge, sense, start, end| Fin {
                        edge,
                        sense,
                        pcurve: Curve2::LineSegment { start, end },
                        enclosure: None,
                    };
                    let lower = topology.add_loop(
                        vec![fin(
                            bottom,
                            forward,
                            Point2::new(u0, 0.0),
                            Point2::new(u1, 0.0),
                        )],
                        [turns, 0],
                    );
                    let upper = topology.add_loop(
                        vec![fin(
                            top,
                            reversed,
                            Point2::new(u1, height),
                            Point2::new(u0, height),
                        )],
                        [-turns, 0],
                    );
                    topology.faces.push(Face {
                        surface: Surface::Cylinder {
                            frame: cylinder_frame,
                            radius: *radius,
                        },
                        sense: forward,
                        loops: vec![lower, upper],
                        front: ShellId(0),
                        back: ShellId(1),
                        enclosure: None,
                    });
                }
            }
        }
        for (k, fin) in topology.fins.iter().enumerate() {
            topology.edges[fin.edge.0].fins.push(FinId(k));
        }
        let fronts = topology.face_ids().map(|f| (f, Side::Front)).collect();
        let backs = topology.face_ids().map(|f| (f, Side::Back)).collect();
        topology.shells = vec![
            Shell {
                region: RegionId(1),
                sides: fronts,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
            Shell {
                region: RegionId(0),
                sides: backs,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
        ];
        topology.regions = vec![
            Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(1)],
            },
            Region {
                kind: RegionKind::Solid,
                shells: vec![ShellId(0)],
            },
        ];
        topology.identity = Identity::new(
            derive(EntityKind::Body, Role::Body, 0, Vec::new()),
            derivations,
            labels,
        )?;
        layout.sort_by_key(|(slot, _)| *slot);
        topology.layout = layout;
        let identity = &topology.identity;
        if (
            identity.vertices.len(),
            identity.edges.len(),
            identity.faces.len(),
            identity.regions.len(),
        ) != (
            topology.vertices.len(),
            topology.edges.len(),
            topology.faces.len(),
            topology.regions.len() - 1,
        ) {
            return Err(Error::InvalidTopology("slot without a derivation"));
        }
        topology.measure_enclosures(tolerance)?;
        topology.validate(tolerance)?;
        if topology.euler_characteristic() != 2 - 2 * profile.holes().len() as i64 {
            return Err(Error::InvalidTopology("unexpected shell genus"));
        }
        Ok(topology)
    }

    /// A right circular cone or frustum (S3 of REVIEW_NOTES.md), as
    /// `BRepPrimAPI_MakeCone(gp_Ax2, bottom, top, height)` makes it: the
    /// bottom at the frame origin, the top at `height` along its normal, a
    /// zero radius an apex. The lateral face is a `Cone` surface with v along
    /// the generatrix from the bottom; a nonzero end is a ring edge bounding
    /// a disc, a zero end a pole (a vertex loop at the apex). Entities derive
    /// from the meridian profile `(0, 0), (bottom, 0), (top, height), (0,
    /// height)` as a revolution: segment 0 (bottom radius) the start cap,
    /// segment 1 (generatrix) the wall, segment 2 (top radius) the end cap,
    /// vertex 1 the bottom ring or apex, vertex 2 the top ring or apex.
    pub(crate) fn cone(
        frame: Frame3,
        bottom: f64,
        top: f64,
        height: f64,
        tolerance: Tolerance,
        operation: OperationId,
    ) -> Result<Self> {
        let tol = tolerance.linear();
        for (value, what) in [
            (bottom, "cone radius"),
            (top, "cone radius"),
            (height, "cone height"),
        ] {
            crate::math::finite(value, what)?;
        }
        if bottom < 0.0 || top < 0.0 || (bottom > 0.0 && bottom <= tol) || (top > 0.0 && top <= tol)
        {
            return Err(Error::Degenerate("cone radius"));
        }
        if bottom == top {
            return Err(Error::OutOfDomain("equal cone radii make a cylinder"));
        }
        if height <= tol {
            return Err(Error::Degenerate("cone height"));
        }
        let derive = |entity, role, parents| Derivation {
            operation,
            kind: OperationKind::Revolve,
            entity,
            role,
            ordinal: 0,
            parents,
        };
        let meridian = |element| Parent::Profile {
            boundary: 0,
            element,
        };
        let half_angle = (top - bottom).atan2(height);
        let slant = height.hypot(top - bottom);
        let mut topology = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            fins: Vec::new(),
            loops: Vec::new(),
            faces: Vec::new(),
            shells: Vec::new(),
            regions: Vec::new(),
            identity: Identity::new(
                derive(EntityKind::Body, Role::Body, Vec::new()),
                Vec::new(),
                BTreeMap::new(),
            )?,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
        };
        let mut derivations: Vec<(Slot, Derivation)> = vec![(
            Slot::Region(RegionId(1)),
            derive(
                EntityKind::Region,
                Role::Region,
                vec![meridian(ProfileElement::Boundary)],
            ),
        )];
        let mut wall_loops = Vec::new();
        // (radius, height, cap role, cap segment, edge role, rim vertex)
        let ends = [
            (bottom, 0.0, Role::StartCap, 0, Role::BottomEdge, 1),
            (top, height, Role::EndCap, 2, Role::TopEdge, 2),
        ];
        for (radius, z, cap_role, segment, edge_role, rim) in ends {
            let upper = z > 0.0;
            if radius == 0.0 {
                let apex = topology.add_vertex(frame.point(Point2::default(), z));
                derivations.push((
                    Slot::Vertex(apex),
                    derive(
                        EntityKind::Vertex,
                        Role::Apex,
                        vec![meridian(ProfileElement::Vertex(rim))],
                    ),
                ));
                topology.loops.push(Loop::Vertex(apex));
                wall_loops.push(LoopId(topology.loops.len() - 1));
                continue;
            }
            let centre = frame.point(Point2::default(), z);
            let normal = if upper {
                frame.normal()
            } else {
                -frame.normal()
            };
            let disc_frame = Frame3::new(centre, normal, frame.x(), tolerance)?;
            let ring_frame = Frame3::new(centre, frame.normal(), frame.x(), tolerance)?;
            let ring = topology.add_ring(Curve3::Circle {
                frame: ring_frame,
                radius,
            });
            derivations.push((
                Slot::Edge(ring),
                derive(
                    EntityKind::Edge,
                    edge_role,
                    vec![meridian(ProfileElement::Vertex(rim))],
                ),
            ));
            let disc = FaceId(topology.faces.len());
            derivations.push((
                Slot::Face(disc),
                derive(
                    EntityKind::Face,
                    cap_role,
                    vec![meridian(ProfileElement::Segment(segment))],
                ),
            ));
            topology.faces.push(Face {
                surface: Surface::Plane(disc_frame),
                sense: Orientation::Forward,
                loops: Vec::new(),
                front: ShellId(0),
                back: ShellId(1),
                enclosure: None,
            });
            // The bottom disc runs against the ring, the top disc with it.
            topology.add_cap_loop(disc.0, &[ring], !upper, disc_frame);
            // The wall: the bottom ring once in +u at v = 0, the top ring once
            // in -u at v = slant, on the universal cover.
            let v = if upper { slant } else { 0.0 };
            let (sense, u0, u1, turns) = if upper {
                (Orientation::Reversed, TAU, 0.0, -1)
            } else {
                (Orientation::Forward, 0.0, TAU, 1)
            };
            let fin = Fin {
                edge: ring,
                sense,
                pcurve: Curve2::LineSegment {
                    start: Point2::new(u0, v),
                    end: Point2::new(u1, v),
                },
                enclosure: None,
            };
            wall_loops.push(topology.add_loop(vec![fin], [turns, 0]));
        }
        let wall = FaceId(topology.faces.len());
        derivations.push((
            Slot::Face(wall),
            derive(
                EntityKind::Face,
                Role::Wall,
                vec![meridian(ProfileElement::Segment(1))],
            ),
        ));
        topology.faces.push(Face {
            surface: Surface::Cone {
                frame,
                radius: bottom,
                half_angle,
            },
            sense: Orientation::Forward,
            loops: wall_loops,
            front: ShellId(0),
            back: ShellId(1),
            enclosure: None,
        });
        for (k, fin) in topology.fins.iter().enumerate() {
            topology.edges[fin.edge.0].fins.push(FinId(k));
        }
        let fronts = topology.face_ids().map(|f| (f, Side::Front)).collect();
        let backs = topology.face_ids().map(|f| (f, Side::Back)).collect();
        topology.shells = vec![
            Shell {
                region: RegionId(1),
                sides: fronts,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
            Shell {
                region: RegionId(0),
                sides: backs,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
        ];
        topology.regions = vec![
            Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(1)],
            },
            Region {
                kind: RegionKind::Solid,
                shells: vec![ShellId(0)],
            },
        ];
        topology.identity = Identity::new(
            derive(EntityKind::Body, Role::Body, Vec::new()),
            derivations,
            BTreeMap::new(),
        )?;
        topology.measure_enclosures(tolerance)?;
        topology.validate(tolerance)?;
        Ok(topology)
    }

    /// A builder's own enclosures: every gap of the constructed geometry,
    /// measured. A construction that cannot be enclosed within the resolution
    /// fails; nothing widens to succeed (Contract 5, T4).
    fn measure_enclosures(&mut self, tolerance: Tolerance) -> Result<()> {
        let m = validate::measure(&self.view());
        let bounds = m.vertices.iter().chain(&m.fins).chain(&m.faces);
        for bound in bounds {
            match bound {
                None => return Err(Error::Unrepresentable("an entity's gap enclosure")),
                Some(b) if *b > tolerance.linear() => return Err(Error::PrecisionLoss),
                Some(_) => {}
            }
        }
        fill(&mut self.vertices, &m.vertices, |v| &mut v.enclosure);
        fill(&mut self.fins, &m.fins, |f| &mut f.enclosure);
        fill(&mut self.faces, &m.faces, |f| &mut f.enclosure);
        Ok(())
    }

    fn add_vertex(&mut self, position: Point3) -> VertexId {
        let id = VertexId(self.vertices.len());
        self.vertices.push(Vertex {
            position,
            enclosure: None,
        });
        id
    }
    fn add_edge(
        &mut self,
        start: Option<VertexId>,
        end: Option<VertexId>,
        curve: Curve3,
    ) -> EdgeId {
        let id = EdgeId(self.edges.len());
        self.edges.push(Edge {
            start,
            end,
            curve,
            fins: Vec::new(),
        });
        id
    }
    fn add_ring(&mut self, curve: Curve3) -> EdgeId {
        self.add_edge(None, None, curve)
    }
    fn add_line(&mut self, start: VertexId, end: VertexId) -> EdgeId {
        self.add_edge(
            Some(start),
            Some(end),
            Curve3::LineSegment {
                start: self.vertices[start.0].position,
                end: self.vertices[end.0].position,
            },
        )
    }
    fn add_loop(&mut self, fins: Vec<Fin>, winding: [i32; 2]) -> LoopId {
        let ids = fins
            .into_iter()
            .map(|fin| {
                self.fins.push(fin);
                FinId(self.fins.len() - 1)
            })
            .collect();
        self.loops.push(Loop::Edges { fins: ids, winding });
        LoopId(self.loops.len() - 1)
    }
    fn add_cap_loop(&mut self, face: usize, edges: &[EdgeId], reverse: bool, frame: Frame3) {
        let sense = if reverse {
            Orientation::Reversed
        } else {
            Orientation::Forward
        };
        let mut fins = edges
            .iter()
            .map(|edge| self.plane_fin(*edge, sense, frame))
            .collect::<Vec<_>>();
        if reverse {
            fins.reverse();
        }
        let l = self.add_loop(fins, [0, 0]);
        self.faces[face].loops.push(l);
    }
    fn plane_fin(&self, id: EdgeId, sense: Orientation, frame: Frame3) -> Fin {
        let edge = &self.edges[id.0];
        let local = |point: Point3| {
            let [x, y, _] = frame.coordinates(point);
            Point2::new(x, y)
        };
        let pcurve = match &edge.curve {
            Curve3::LineSegment { start, end } => {
                let (a, b) = if sense == Orientation::Forward {
                    (*start, *end)
                } else {
                    (*end, *start)
                };
                Curve2::LineSegment {
                    start: local(a),
                    end: local(b),
                }
            }
            Curve3::Circle {
                frame: circle,
                radius,
            } => {
                let [x, y, _] = frame.coordinates(circle.origin());
                let start = local(circle.point(Point2::new(*radius, 0.0), 0.0));
                Curve2::CircularArc {
                    center: Point2::new(x, y),
                    radius: *radius,
                    start_angle: (start.y - y).atan2(start.x - x),
                    sweep_angle: TAU * sense.sign() * circle.normal().dot(frame.normal()).signum(),
                }
            }
            Curve3::CircularArc { .. } => {
                unreachable!("the extrusion builder creates only lines and full circles")
            }
        };
        Fin {
            edge: id,
            sense,
            pcurve,
            enclosure: None,
        }
    }
}
