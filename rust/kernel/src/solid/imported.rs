//! S9e.4a: imported solids (REVIEW_NOTES.md, "S9e.4 refined").
//!
//! A body read from a `.brep` file or a STEP file has no construction: the
//! converter gives a cell topology and the body's resolution (the largest
//! OCCT tolerance). What it stores is rounded (OCCT writes surfaces and
//! curves with 17 significant digits, vertices with 15; the converter
//! normalizes every frame again), so no stored vertex lies on its faces'
//! stored surfaces exactly. An imported solid whose stored topology is one
//! of the kernel's constructions (a prism of lines, arcs and circles; a
//! sphere, cap or zone; a cone or frustum; a whole torus) is that
//! construction, read off its stored surfaces once (a cap's plane frame and
//! its loops' vertices and arcs rounded once into it; a sphere's, cone's or
//! torus's stored frame and radii exactly, its ends' heights from their
//! rings), matched to the stored topology by S9e.2's geometric match
//! (`curved::matched`): every stored vertex within the resolution of the
//! construction's, edges through their points, faces by their edges and
//! surface kinds. The body is decided on the construction's exact model: a
//! Boolean takes the construction in its place (`polyhedra::build`) and
//! names the result over the stored ids through the match; it classifies
//! and moves with it. Every other body is S9e.4b's.
use super::{replayable, Construction, Context, Solid};
use crate::history::{History, Relation};
use crate::identity::{
    AlgorithmLevel, Derivation, EntityId, EntityKind, InputLabel, OperationId, OperationKind,
    Parent, Role,
};
use crate::solid::split::{q, rational_f64, zero};
use crate::topology::{
    Curve3, EdgeId, FaceId, Loop, Orientation, RegionId, RegionKind, Slot, Surface, Topology,
    VertexId,
};
use crate::{
    Boundary, Error, Frame3, Point2, Point3, Profile, Result, RigidTransform, Segment, Tolerance,
    Vec3,
};
use num_rational::BigRational as R;
use std::collections::BTreeMap;

/// An imported solid's construction and names (S9e.4a).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Imported {
    /// The construction its stored surfaces give, with ids of its own.
    pub(crate) recognized: Box<Solid>,
    /// Each construction entity's stored entity (the match).
    pub(crate) names: BTreeMap<EntityId, EntityId>,
    pub(crate) tolerance: Tolerance,
}

fn general() -> Error {
    Error::OutOfDomain("an imported solid other than a prism, a sphere, a cone or a torus (S9e.4b)")
}

/// Candidates' parallel and perpendicular tests (the match decides).
const ALIGN: f64 = 1e-9;

impl Solid {
    /// A solid from an imported topology (S9e.4a): a `.brep` or STEP
    /// converter's solid and its resolution. Its entities are renamed under
    /// the operation (each slot its ordinal, as `Topology::from_parts` names
    /// them, with the operation's id), and every one is `Generated` in the
    /// history from a label of its slot in the file's cells (its kind above
    /// its ordinal: the import's input). The solid must be one of the kernel's
    /// constructions read off its stored surfaces (a prism of lines, arcs
    /// and circles, a sphere, cap or zone, a cone or frustum, a whole
    /// torus), within the resolution of every stored vertex and edge;
    /// otherwise `OutOfDomain` (S9e.4b; spline faces or edges S9f).
    pub fn imported_with(
        operation: OperationId,
        topology: Topology,
        resolution: Tolerance,
    ) -> Result<(Self, History)> {
        Self::imported_in(&Context::new(operation), topology, resolution)
    }

    /// [`Solid::imported_with`] at a recorded algorithm level (H8).
    pub fn imported_at(
        level: AlgorithmLevel,
        operation: OperationId,
        topology: Topology,
        resolution: Tolerance,
    ) -> Result<(Self, History)> {
        Self::imported_in(&Context::new(operation).at(level), topology, resolution)
    }

    /// [`Solid::imported_with`] in an operation context (its level; an
    /// import carries no attributes).
    pub fn imported_in(
        context: &Context,
        topology: Topology,
        resolution: Tolerance,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        replayable(level)?;
        let topology = renamed(topology, operation)?;
        let recognized = recognize(&topology, resolution, construction_operation(operation))?;
        let names = names(&recognized, &topology, resolution)?;
        let mass = topology
            .mass_enclosure()
            .ok_or(Error::Unrepresentable(
                "an imported solid's mass properties",
            ))?
            .midpoints();
        let bounds = widened(recognized.bounds, &topology);
        let solid = Self {
            frame: recognized.frame,
            start: recognized.start,
            end: recognized.end,
            construction: Construction::Imported(Box::new(Imported {
                recognized: Box::new(recognized),
                names,
                tolerance: resolution,
            })),
            topology,
            mass,
            bounds,
            operation,
        };
        // Each entity generated from its slot in the file's cells (the
        // import's input), a label of its kind and ordinal.
        let t = &solid.topology;
        let relations = t
            .ids()
            .map(|(id, slot)| Relation::Generated {
                from: vec![Parent::Label(slot_label(slot))],
                to: id,
                role: Role::External,
            })
            .collect();
        let history = History::new(
            operation,
            OperationKind::External,
            Vec::new(),
            vec![t.body_id()],
            relations,
            Vec::new(),
        )
        .at_level(level);
        solid.debug_check(&[], &history);
        Ok((solid, history))
    }
}

impl Imported {
    /// The solid moved rigidly: its stored topology moved (as S9b's
    /// results' stored geometry) and its construction rebuilt in the moved
    /// frame, the match kept (both keep their ids).
    pub(crate) fn moved(
        &self,
        topology: &Topology,
        operation: OperationId,
        motion: RigidTransform,
        mass: super::MassProperties,
    ) -> Result<Solid> {
        let tolerance = self.tolerance;
        let parts = topology.moved_parts(motion, tolerance)?;
        let moved =
            Topology::from_parts(parts.with_measured_enclosures(), tolerance).map_err(|i| {
                Error::InvalidTopology(
                    i.first()
                        .map_or("a moved imported solid", |i| i.kind.name()),
                )
            })?;
        let moved = moved.with_identity_of(topology);
        let recognized = self
            .recognized
            .transform_with(self.recognized.operation, motion)?
            .0;
        let bounds = widened(recognized.bounds, &moved);
        Ok(Solid {
            frame: recognized.frame,
            start: recognized.start,
            end: recognized.end,
            construction: Construction::Imported(Box::new(Imported {
                recognized: Box::new(recognized),
                names: self.names.clone(),
                tolerance,
            })),
            topology: moved,
            mass,
            bounds,
            operation,
        })
    }
}

/// A stored slot's label in an import's history: its kind (vertex 0, edge
/// 1, face 2, region 3) above its ordinal.
fn slot_label(slot: Slot) -> InputLabel {
    let (kind, k) = match slot {
        Slot::Vertex(v) => (0u64, v.0),
        Slot::Edge(e) => (1, e.0),
        Slot::Face(f) => (2, f.0),
        Slot::Region(r) => (3, r.0),
    };
    InputLabel((kind << 48) | k as u64)
}

/// The operation whose ids the construction takes: the import's, its bits
/// turned (apart from the import's own ids and from another import's).
fn construction_operation(operation: OperationId) -> OperationId {
    OperationId(!operation.0)
}

/// The stored topology under the operation's ids.
fn renamed(topology: Topology, operation: OperationId) -> Result<Topology> {
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

/// The construction's bounds widened by the stored edges'.
fn widened(bounds: crate::Bounds3, t: &Topology) -> crate::Bounds3 {
    let edges = super::split::edge_bounds(t);
    let (a, b, c, d) = (
        bounds.min.to_array(),
        bounds.max.to_array(),
        edges.min.to_array(),
        edges.max.to_array(),
    );
    crate::Bounds3 {
        min: Point3::new(a[0].min(c[0]), a[1].min(c[1]), a[2].min(c[2])),
        max: Point3::new(b[0].max(d[0]), b[1].max(d[1]), b[2].max(d[2])),
    }
}

/// The construction matched to the stored topology: each construction
/// entity's stored entity, or `OutOfDomain` (S9e.4b) where they differ.
fn names(
    recognized: &Solid,
    stored: &Topology,
    resolution: Tolerance,
) -> Result<BTreeMap<EntityId, EntityId>> {
    let r = &recognized.topology;
    let m = super::boolean::curved::matched::matched(&r.to_parts(), stored, resolution.linear())
        .ok_or(Error::OutOfDomain(
            "an imported solid unlike the construction its stored surfaces give (S9e.4b)",
        ))?;
    let mut out = BTreeMap::new();
    out.insert(r.body_id(), stored.body_id());
    for (id, slot) in r.ids() {
        let to = match slot {
            Slot::Vertex(v) => Slot::Vertex(VertexId(m.vertices[v.0])),
            Slot::Edge(e) => Slot::Edge(EdgeId(m.edges[e.0])),
            Slot::Face(f) => Slot::Face(FaceId(m.faces[f.0])),
            Slot::Region(g) => Slot::Region(g),
        };
        let stored_id = stored
            .id_of(to)
            .ok_or(Error::InvalidTopology("an imported solid's unnamed slot"))?;
        out.insert(id, stored_id);
    }
    Ok(out)
}

// ------------------------------------------------------------------ recognition

/// Exact rationals of a binary64 vector.
fn qv(v: [f64; 3]) -> [R; 3] {
    v.map(q)
}

fn det3(a: &[R; 3], b: &[R; 3], c: &[R; 3]) -> R {
    &a[0] * (&b[1] * &c[2] - &b[2] * &c[1]) - &a[1] * (&b[0] * &c[2] - &b[2] * &c[0])
        + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0])
}

/// A frame's exact affine map on its stored axes: a point's local
/// coordinates `(u, v, w)`, `p = o + u x + v y + w n` (Cramer's rule).
struct Local {
    o: [R; 3],
    x: [R; 3],
    y: [R; 3],
    n: [R; 3],
    det: R,
}

impl Local {
    fn new(f: &Frame3) -> Result<Self> {
        let (o, x, y, n) = (
            qv(f.origin().to_array()),
            qv(f.x().to_array()),
            qv(f.y().to_array()),
            qv(f.normal().to_array()),
        );
        let det = det3(&x, &y, &n);
        if det == zero() {
            return Err(Error::Degenerate("a frame's axes"));
        }
        Ok(Self { o, x, y, n, det })
    }

    fn of(&self, p: Point3) -> [R; 3] {
        let p = qv(p.to_array());
        let d: [R; 3] = std::array::from_fn(|i| &p[i] - &self.o[i]);
        [
            det3(&d, &self.y, &self.n) / &self.det,
            det3(&self.x, &d, &self.n) / &self.det,
            det3(&self.x, &self.y, &d) / &self.det,
        ]
    }

    /// The point at local height `w` on the frame's axis, rounded once.
    fn axis_point(&self, w: &R) -> Point3 {
        let c: [f64; 3] = std::array::from_fn(|i| rational_f64(&(&self.o[i] + w * &self.n[i])));
        Point3::new(c[0], c[1], c[2])
    }
}

/// A local point's first two coordinates rounded once.
fn point2(l: &[R; 3]) -> Point2 {
    Point2::new(rational_f64(&l[0]), rational_f64(&l[1]))
}

/// The face's outward normal at a point of its surface's frame (a plane's
/// normal; the solid region behind the oriented normal's front).
fn outward(t: &Topology, fi: usize) -> Option<Vec3> {
    let f = &t.faces()[fi];
    let Surface::Plane(frame) = &f.surface else {
        return None;
    };
    let oriented = frame.normal() * f.sense.sign();
    let front = &t.regions()[t.shells()[f.front.index()].region.index()];
    Some(if front.kind == RegionKind::Solid {
        oriented
    } else {
        -oriented
    })
}

fn parallel(a: Vec3, b: Vec3) -> bool {
    a.cross(b).length() <= ALIGN * a.length() * b.length()
}

fn perpendicular(a: Vec3, b: Vec3) -> bool {
    a.dot(b).abs() <= ALIGN * a.length() * b.length()
}

/// The construction a stored topology is (S9e.4 refined, (3)).
fn recognize(t: &Topology, tolerance: Tolerance, op: OperationId) -> Result<Solid> {
    if t.faces()
        .iter()
        .any(|f| matches!(f.surface, Surface::BSpline(_)))
        || t.edges()
            .iter()
            .any(|e| matches!(e.curve, Curve3::BSpline(_)))
    {
        return Err(Error::OutOfDomain(
            "an imported solid with spline faces or edges (S9f)",
        ));
    }
    // One solid region with one shell (no cavity), its void twin outside.
    let solids: Vec<_> = t
        .regions()
        .iter()
        .filter(|r| r.kind == RegionKind::Solid)
        .collect();
    if solids.len() != 1 || solids[0].shells.len() != 1 {
        return Err(general());
    }
    let count = |pick: fn(&Surface) -> bool| t.faces().iter().filter(|f| pick(&f.surface)).count();
    let spheres = count(|s| matches!(s, Surface::Sphere { .. }));
    let cones = count(|s| matches!(s, Surface::Cone { .. }));
    let tori = count(|s| matches!(s, Surface::Torus { .. }));
    let planes = count(|s| matches!(s, Surface::Plane(_)));
    let n = t.faces().len();
    let built = match (spheres, cones, tori) {
        (1, 0, 0) if planes == n - 1 => sphere(t, tolerance, op),
        (0, 1, 0) if planes == n - 1 => cone(t, tolerance, op),
        (0, 0, 1) if n == 1 => torus(t, tolerance, op),
        (0, 0, 0) => prism(t, tolerance, op),
        _ => Err(general()),
    };
    // A construction its stored data cannot build (a profile rounded onto
    // itself, a degenerate frame or height) is none of the kernel's: the
    // body is S9e.4b's.
    built.map_err(|e| match e {
        Error::OutOfDomain(_) => e,
        _ => general(),
    })
}

/// A prism: its caps, frame, profile and heights.
fn prism(t: &Topology, tolerance: Tolerance, op: OperationId) -> Result<Solid> {
    let faces = t.faces();
    let normals: Vec<Option<Vec3>> = (0..faces.len()).map(|i| outward(t, i)).collect();
    // The caps: the first pair of planes facing apart with every other face
    // along their normal.
    let mut caps = None;
    'pairs: for i in 0..faces.len() {
        let Some(mi) = normals[i] else { continue };
        for j in i + 1..faces.len() {
            let Some(mj) = normals[j] else { continue };
            if mi.dot(mj) >= 0.0 || !parallel(mi, mj) {
                continue;
            }
            let up = -mi;
            let walls = (0..faces.len()).filter(|&k| k != i && k != j).all(|k| {
                match (&faces[k].surface, normals[k]) {
                    (Surface::Plane(_), Some(m)) => perpendicular(m, up),
                    (Surface::Cylinder { frame, .. }, _) => parallel(frame.normal(), up),
                    _ => false,
                }
            });
            if walls {
                caps = Some((i, j, up));
                break 'pairs;
            }
        }
    }
    let (bottom, top, up) = caps.ok_or_else(general)?;
    let Surface::Plane(stored) = &faces[bottom].surface else {
        return Err(general());
    };
    let frame = if stored.normal().dot(up) > 0.0 {
        *stored
    } else {
        stored.flipped()
    };
    let local = Local::new(&frame)?;
    let mut boundaries = Vec::new();
    for l in &faces[bottom].loops {
        let Loop::Edges { fins, .. } = &t.loops()[l.index()] else {
            return Err(general());
        };
        let mut points = Vec::new();
        let mut segments = Vec::new();
        let mut ring = None;
        for f in fins {
            let fin = &t.fins()[f.index()];
            let edge = &t.edges()[fin.edge.index()];
            let forward = fin.sense == Orientation::Forward;
            if edge.is_ring() {
                if fins.len() != 1 {
                    return Err(general());
                }
                let (Curve3::Circle { frame: c, radius }
                | Curve3::CircularArc {
                    frame: c, radius, ..
                }) = &edge.curve
                else {
                    return Err(general());
                };
                ring = Some((point2(&local.of(c.origin())), *radius));
                continue;
            }
            let start = if forward { edge.start } else { edge.end }.ok_or_else(general)?;
            points.push(point2(&local.of(t.vertices()[start.index()].position)));
            segments.push(match &edge.curve {
                Curve3::LineSegment { .. } => Segment::Line,
                Curve3::CircularArc {
                    frame: c,
                    radius,
                    sweep_angle,
                    ..
                } => {
                    let turns = (*sweep_angle > 0.0) == forward;
                    let about = c.normal().dot(up) > 0.0;
                    Segment::Arc {
                        center: point2(&local.of(c.origin())),
                        radius: *radius,
                        ccw: turns == about,
                    }
                }
                _ => return Err(general()),
            });
        }
        boundaries.push(match ring {
            Some((center, radius)) => Boundary::circle(center, radius, tolerance)?,
            None if segments.iter().all(|s| matches!(s, Segment::Line)) => {
                Boundary::polygon(points, tolerance)?
            }
            None => Boundary::path(points, segments, tolerance)?,
        });
    }
    if boundaries.is_empty() {
        return Err(general());
    }
    // The outer loop first on a plane (`Face::loops`).
    let holes = boundaries.split_off(1);
    let profile = Profile::new(boundaries.remove(0), holes, tolerance)?;
    let Surface::Plane(top_frame) = &faces[top].surface else {
        return Err(general());
    };
    let height = rational_f64(&local.of(top_frame.origin())[2]);
    if height <= 0.0 {
        return Err(general());
    }
    Solid::build(op, profile, frame, 0.0, height)
}

/// An end of a sphere or cone: its local height, its radius (zero at a
/// pole or apex) and whether it faces along the axis.
type End = (R, f64, bool);

/// The discs of a sphere or a cone (each plane face bounded by one ring),
/// in the curved face's frame, and its vertex loops' heights.
fn ends(t: &Topology, curved: usize, local: &Local, axis: Vec3) -> Result<(Vec<End>, Vec<R>)> {
    let mut out = Vec::new();
    for (k, f) in t.faces().iter().enumerate() {
        if k == curved {
            continue;
        }
        let m = outward(t, k).ok_or_else(general)?;
        if !parallel(m, axis) || f.loops.len() != 1 {
            return Err(general());
        }
        let Loop::Edges { fins, .. } = &t.loops()[f.loops[0].index()] else {
            return Err(general());
        };
        let [fin] = fins.as_slice() else {
            return Err(general());
        };
        let edge = &t.edges()[t.fins()[fin.index()].edge.index()];
        let (Curve3::Circle { frame: c, radius }
        | Curve3::CircularArc {
            frame: c, radius, ..
        }) = &edge.curve
        else {
            return Err(general());
        };
        if !edge.is_ring() {
            return Err(general());
        }
        let h = local.of(c.origin())[2].clone();
        out.push((h, *radius, m.dot(axis) > 0.0));
    }
    let mut points = Vec::new();
    for l in &t.faces()[curved].loops {
        if let Loop::Vertex(v) = &t.loops()[l.index()] {
            points.push(local.of(t.vertices()[v.index()].position)[2].clone());
        }
    }
    Ok((out, points))
}

/// A sphere, cap or zone: the stored frame and radius, each disc's
/// latitude from its ring's height.
fn sphere(t: &Topology, tolerance: Tolerance, op: OperationId) -> Result<Solid> {
    let k = t
        .faces()
        .iter()
        .position(|f| matches!(f.surface, Surface::Sphere { .. }))
        .ok_or_else(general)?;
    let Surface::Sphere { frame, radius } = &t.faces()[k].surface else {
        return Err(general());
    };
    let local = Local::new(frame)?;
    let (discs, _) = ends(t, k, &local, frame.normal())?;
    let (mut low, mut high) = (None, None);
    for (h, _, up) in discs {
        let latitude = (rational_f64(&h) / radius).clamp(-1.0, 1.0).asin();
        let slot = if up { &mut high } else { &mut low };
        if slot.replace(latitude).is_some() {
            return Err(general());
        }
    }
    let half = std::f64::consts::FRAC_PI_2;
    Solid::build_sphere(
        op,
        *frame,
        *radius,
        low.unwrap_or(-half),
        high.unwrap_or(half),
        tolerance,
    )
}

/// A cone or frustum: the stored frame moved along its axis to the lower
/// end, the ends' rings' radii (zero at the apex) and their heights.
fn cone(t: &Topology, tolerance: Tolerance, op: OperationId) -> Result<Solid> {
    let k = t
        .faces()
        .iter()
        .position(|f| matches!(f.surface, Surface::Cone { .. }))
        .ok_or_else(general)?;
    let Surface::Cone { frame, .. } = &t.faces()[k].surface else {
        return Err(general());
    };
    let local = Local::new(frame)?;
    let (discs, apexes) = ends(t, k, &local, frame.normal())?;
    let mut all: Vec<(R, f64)> = discs.into_iter().map(|(h, r, _)| (h, r)).collect();
    all.extend(apexes.into_iter().map(|h| (h, 0.0)));
    let [a, b] = <[(R, f64); 2]>::try_from(all).map_err(|_| general())?;
    let (lo, hi) = if a.0 <= b.0 { (a, b) } else { (b, a) };
    let base = if lo.0 == zero() {
        *frame
    } else {
        frame.at(local.axis_point(&lo.0))
    };
    let height = rational_f64(&(&hi.0 - &lo.0));
    Solid::build_cone(op, base, lo.1, hi.1, height, tolerance)
}

/// A whole torus: the stored frame and radii.
fn torus(t: &Topology, tolerance: Tolerance, op: OperationId) -> Result<Solid> {
    let f = &t.faces()[0];
    let Surface::Torus {
        frame,
        major,
        minor,
    } = &f.surface
    else {
        return Err(general());
    };
    if !f.loops.is_empty() {
        return Err(general());
    }
    let tau = std::f64::consts::TAU;
    Solid::build_torus(op, *frame, *major, *minor, 0.0, tau, tau, tolerance)
}

/// Whether every arc of a prism's profile ends on its circle exactly
/// (S9c's model: an imported profile rounded into its cap's frame may not).
pub(crate) fn arcs_on_circles(s: &Solid) -> bool {
    let Construction::Prism(p) = &s.construction else {
        return true;
    };
    p.boundaries().all(|b| match b.path_geometry() {
        Some((points, segments)) => segments.iter().enumerate().all(|(j, seg)| {
            let Segment::Arc { center, radius, .. } = seg else {
                return true;
            };
            let on = |p: Point2| {
                let (dx, dy) = (q(p.x) - q(center.x), q(p.y) - q(center.y));
                &dx * &dx + &dy * &dy == q(*radius) * q(*radius)
            };
            on(points[j]) && on(points[(j + 1) % points.len()])
        }),
        None => true,
    })
}
