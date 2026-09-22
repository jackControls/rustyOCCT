//! Analytic boundary representation for the supported extrusion builders.
//!
//! Edges are shared by oriented face uses. Each use retains its own parameter
//! curve, including the two distinct parameter curves of a cylindrical seam.
//! IDs are body-local indices, not persistent names across arbitrary edits.
use crate::profile::BoundaryKind;
use crate::{Error, Frame3, Point2, Point3, Profile, Result, Tolerance, Vec3};
use std::f64::consts::TAU;

macro_rules! index_type {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(pub(crate) usize);
        impl $name {
            pub fn index(self) -> usize {
                self.0
            }
        }
    };
}
index_type!(VertexId);
index_type!(EdgeId);
index_type!(FaceId);

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
    LineSegment { start: Point3, end: Point3 },
    Circle { frame: Frame3, radius: f64 },
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
    Wall { boundary: usize, segment: usize },
}

#[derive(Debug, Clone, PartialEq)]
pub struct Face {
    pub surface: Surface,
    pub orientation: Orientation,
    /// Outer loop first, then inner loops; traversal follows the oriented face.
    pub loops: Vec<Vec<Coedge>>,
    pub origin: FaceOrigin,
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
}

impl Topology {
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

    /// Validate connectivity, opposite edge uses, geometric endpoints,
    /// connectedness and sampled agreement of 3D curves with face pcurves.
    /// This is not a general self-intersection or sewing/healing algorithm.
    pub fn validate(&self, tolerance: Tolerance) -> Result<()> {
        let invalid = Error::InvalidTopology;
        if self.faces.is_empty() || self.edges.is_empty() || self.vertices.is_empty() {
            return Err(invalid("empty shell"));
        }
        let mut vertex_used = vec![false; self.vertices.len()];
        for vertex in &self.vertices {
            vertex.position.checked(tolerance)?;
        }
        for edge in &self.edges {
            for (id, fraction) in [(edge.start, 0.0), (edge.end, 1.0)] {
                let vertex = self.vertex(id).ok_or(invalid("invalid vertex id"))?;
                let evaluated = edge.curve.point(fraction).checked(tolerance)?;
                if vertex.position.distance(evaluated) > tolerance.linear() {
                    return Err(invalid("edge endpoint disagrees with its vertex"));
                }
                vertex_used[id.0] = true;
            }
        }
        if vertex_used.contains(&false) {
            return Err(invalid("unreferenced vertex"));
        }
        let mut uses = vec![Vec::new(); self.edges.len()];
        for (face_index, face) in self.faces.iter().enumerate() {
            if face.loops.is_empty() {
                return Err(invalid("face has no boundary"));
            }
            for wire in &face.loops {
                if wire.is_empty() {
                    return Err(invalid("empty wire"));
                }
                for (i, coedge) in wire.iter().enumerate() {
                    let edge = self.edge(coedge.edge).ok_or(invalid("invalid edge id"))?;
                    let (_, end) = oriented_vertices(edge, coedge.orientation);
                    let next = &wire[(i + 1) % wire.len()];
                    let next_edge = self.edge(next.edge).ok_or(invalid("invalid edge id"))?;
                    if end != oriented_vertices(next_edge, next.orientation).0 {
                        return Err(invalid("open wire"));
                    }
                    uses[coedge.edge.0].push((face_index, coedge.orientation));
                    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                        let parameter = if coedge.orientation == Orientation::Forward {
                            t
                        } else {
                            1.0 - t
                        };
                        let point = edge.curve.point(parameter).checked(tolerance)?;
                        let on_face = face
                            .surface
                            .point(coedge.pcurve.point(t))
                            .checked(tolerance)?;
                        if point.distance(on_face) > tolerance.linear() {
                            return Err(invalid("pcurve disagrees with 3D edge"));
                        }
                    }
                }
            }
        }
        let mut neighbors = vec![Vec::new(); self.faces.len()];
        for pair in uses {
            if pair.len() != 2 || pair[0].1 == pair[1].1 {
                return Err(invalid("edge must have exactly two opposite uses"));
            }
            neighbors[pair[0].0].push(pair[1].0);
            neighbors[pair[1].0].push(pair[0].0);
        }
        let mut seen = vec![false; self.faces.len()];
        let mut pending = vec![0];
        while let Some(face) = pending.pop() {
            if !seen[face] {
                seen[face] = true;
                pending.extend(&neighbors[face]);
            }
        }
        if seen.contains(&false) {
            return Err(invalid("disconnected shell"));
        }
        Ok(())
    }

    pub(crate) fn prism(
        profile: &Profile,
        frame: Frame3,
        low: f64,
        high: f64,
        start_is_low: bool,
    ) -> Result<Self> {
        let tolerance = profile.tolerance();
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
                    origin: if start_is_low {
                        FaceOrigin::StartCap
                    } else {
                        FaceOrigin::EndCap
                    },
                },
                Face {
                    surface: Surface::Plane(top_frame),
                    orientation: Orientation::Forward,
                    loops: Vec::new(),
                    origin: if start_is_low {
                        FaceOrigin::EndCap
                    } else {
                        FaceOrigin::StartCap
                    },
                },
            ],
        };
        for (boundary, wire) in profile.boundaries().enumerate() {
            let inner = boundary > 0;
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
                        topology.faces.push(Face {
                            surface: Surface::Plane(side_frame),
                            orientation: Orientation::Forward,
                            loops: vec![coedges],
                            origin: FaceOrigin::Wall {
                                boundary,
                                segment: i,
                            },
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
                        origin: FaceOrigin::Wall {
                            boundary,
                            segment: 0,
                        },
                    });
                }
            }
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
