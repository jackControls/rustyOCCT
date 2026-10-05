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
//! and moves with it. S9e.4b.2: a body of plane faces and line edges that is
//! no such prism is a polyhedron decided on its stored vertices
//! (`boolean::polyhedra::imported`: S9b.2's stored model, each face the
//! polygon of its stored vertices, cut into exactly planar triangles where
//! they are not coplanar), its entities its stored ones. S9e.4b.3a: a body
//! of one sphere, cylinder or cone face and plane faces is a plane piece,
//! its primitive common its planes' half-spaces (`boolean::curved::pieces`),
//! its model built and matched to its stored topology on import, its
//! entities its stored ones; S9e.4b.3c.3: or another Boolean of the two (the
//! hull less the primitive, the primitive less the hull, their fuse), each
//! form tried in turn, the match the arbiter. Every other body is S9e.4b's.
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
    /// How it is decided: the construction its stored surfaces give, or its
    /// stored vertices (S9e.4b.2).
    pub(crate) recognized: Recognized,
    /// Each construction entity's stored entity (the match; empty for a
    /// polyhedron, whose entities are its stored ones).
    pub(crate) names: BTreeMap<EntityId, EntityId>,
    pub(crate) tolerance: Tolerance,
}

/// What an imported solid is decided on.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Recognized {
    /// S9e.4a: the construction its stored surfaces give, with ids of its
    /// own.
    Construction(Box<Solid>),
    /// S9e.4b.2: a polyhedron other than a prism, on its stored vertices.
    Polyhedron,
    /// S9e.4b.3a: a plane piece of a sphere, a cylinder or a cone, its
    /// primitive common its planes' half-spaces.
    Piece(Box<Piece>),
}

/// S9e.4b.3a: a plane piece's primitive (a whole sphere, or a cylinder or a
/// cone on the curved face's stored frame reaching past the body's ends,
/// with ids of its own) and its planes (each plane face's stored frame and
/// whether its normal `x * y` points into the material).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Piece {
    pub(crate) primitive: Box<Solid>,
    pub(crate) planes: Vec<(Frame3, bool)>,
    /// S9e.4b.3c.3: the Boolean of the primitive and the hull of the planes
    /// the body is.
    pub(crate) form: Form,
}

/// S9e.4b.3c.3: which Boolean of its primitive and the convex hull of its
/// planes a plane piece is (REVIEW_NOTES.md, "S9e.4b.3c.3 refined").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Form {
    /// The primitive common the hull (S9e.4b.3a): a convex piece.
    Common,
    /// The hull less the primitive: a hole, groove, notch or dimple, the
    /// curved face's material outside its quadric.
    Groove,
    /// The primitive less the hull of its planes turned over: a bite.
    Bite,
    /// The hull fused with the primitive: a boss.
    Boss,
}

impl Form {
    /// The forms a body whose curved face's material lies outside its
    /// quadric (`outside`) or inside it may be, in the order tried.
    fn candidates(outside: bool) -> &'static [Form] {
        if outside {
            &[Form::Groove]
        } else {
            &[Form::Common, Form::Bite, Form::Boss]
        }
    }
}

impl Imported {
    /// Whether it is a polyhedron decided on its stored vertices (S9e.4b.2).
    pub(crate) fn polyhedron(&self) -> bool {
        self.recognized == Recognized::Polyhedron
    }

    /// Whether it is a plane piece of a sphere, a cylinder or a cone
    /// (S9e.4b.3a).
    pub(crate) fn piece(&self) -> bool {
        matches!(self.recognized, Recognized::Piece(_))
    }

    /// A point's location: its construction's (within the resolution), or
    /// the polyhedron's stored model's (S9e.4b.2).
    pub(crate) fn classify(&self, solid: &Solid, point: Point3) -> Result<crate::Location> {
        match &self.recognized {
            Recognized::Construction(s) => s.classify(point),
            Recognized::Polyhedron => super::boolean::polyhedra::imported::classify(solid, point),
            Recognized::Piece(p) => p.classify(point, self.tolerance),
        }
    }
}

impl Piece {
    /// A point's location: its primitive's and its planes' sides (within
    /// the resolution), combined by its form.
    fn classify(&self, point: Point3, tolerance: Tolerance) -> Result<crate::Location> {
        use crate::Location;
        let tol = tolerance.linear();
        let mut hull = Location::Inside;
        for (frame, into) in &self.planes {
            let d = frame.coordinates(point)[2];
            let d = if *into { -d } else { d };
            if d > tol {
                hull = Location::Outside;
                break;
            }
            if d >= -tol {
                hull = Location::Boundary;
            }
        }
        let primitive = self.primitive.classify(point)?;
        // Inside 2, on the boundary 1, outside 0: a common the least, a
        // fuse the most, a complement the reverse.
        let rank = |l: Location| match l {
            Location::Inside => 2,
            Location::Boundary => 1,
            Location::Outside => 0,
        };
        let (p, h) = (rank(primitive), rank(hull));
        let r = match self.form {
            Form::Common => p.min(h),
            Form::Groove => h.min(2 - p),
            Form::Bite => p.min(2 - h),
            Form::Boss => p.max(h),
        };
        Ok([Location::Outside, Location::Boundary, Location::Inside][r])
    }
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
    /// torus), within the resolution of every stored vertex and edge, or
    /// (S9e.4b.2) a polyhedron of one shell, decided on its stored vertices;
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
        // S9e.4a's construction where its recognition and match succeed,
        // else (S9e.4b.2) a polyhedron of one shell on its stored vertices.
        let construction = recognize(&topology, resolution, construction_operation(operation))
            .and_then(|r| names(&r, &topology, resolution).map(|n| (r, n)));
        let (recognized, names) = match construction {
            Ok((r, n)) => (Recognized::Construction(Box::new(r)), n),
            Err(e) if planar(&topology) => {
                if !one_shell(&topology) {
                    return Err(match e {
                        Error::OutOfDomain(_) => Error::OutOfDomain(
                            "an imported polyhedron with a cavity or several shells (S9e.4b.4)",
                        ),
                        e => e,
                    });
                }
                (Recognized::Polyhedron, BTreeMap::new())
            }
            // S9e.4b.3a: a plane piece of a sphere, a cylinder or a cone, its
            // primitive common its planes' half-spaces (checked below).
            // S9e.4b.3c.3: or another Boolean of the two, each form tried in
            // turn (its model's match to the stored topology the arbiter).
            Err(e) if piece_shape(&topology) => {
                return Self::imported_piece(context, topology, resolution, e);
            }
            Err(e) => return Err(e),
        };
        Self::imported_as(context, topology, resolution, recognized, names)
    }

    /// S9e.4b.3a: an imported plane piece, its forms tried in turn
    /// (S9e.4b.3c.3): the first whose model matches the stored topology,
    /// else the first refusal other than a mismatch (a primitive the stored
    /// data cannot build S9e.4a's refusal `e`).
    fn imported_piece(
        context: &Context,
        topology: Topology,
        resolution: Tolerance,
        e: Error,
    ) -> Result<(Self, History)> {
        let op = construction_operation(context.operation);
        let k = curved_face(&topology).ok_or_else(general)?;
        let mut refusal = None;
        for &form in Form::candidates(material_outside(&topology, k)) {
            let built = piece(&topology, resolution, op, form).map_err(|x| {
                if x == not_its_primitive() {
                    x
                } else {
                    e.clone()
                }
            });
            let attempt = built.and_then(|p| {
                Self::imported_as(
                    context,
                    topology.clone(),
                    resolution,
                    Recognized::Piece(Box::new(p)),
                    BTreeMap::new(),
                )
            });
            match attempt {
                Ok(done) => return Ok(done),
                Err(e) if e == not_its_primitive() => {}
                Err(e) => {
                    refusal.get_or_insert(e);
                }
            }
        }
        Err(refusal.unwrap_or_else(not_its_primitive))
    }

    /// The imported solid of what its stored topology was recognized as
    /// (S9e.4a), checked.
    fn imported_as(
        context: &Context,
        topology: Topology,
        resolution: Tolerance,
        recognized: Recognized,
        names: BTreeMap<EntityId, EntityId>,
    ) -> Result<(Self, History)> {
        let (level, operation) = (context.level, context.operation);
        let mass = topology
            .mass_enclosure()
            .ok_or(Error::Unrepresentable(
                "an imported solid's mass properties",
            ))?
            .midpoints();
        let (frame, start, end, bounds) = placed(&recognized, &topology);
        let solid = Self {
            frame,
            start,
            end,
            construction: Construction::Imported(Box::new(Imported {
                recognized,
                names,
                tolerance: resolution,
            })),
            topology,
            mass,
            bounds,
            operation,
        };
        // S9e.4b.2: a polyhedron's stored model builds (no face folded).
        // S9e.4b.3a: a piece's model builds and matches its stored topology.
        if let Construction::Imported(i) = &solid.construction {
            if i.piece() {
                super::boolean::curved::pieces::check(&solid)?;
            }
            if i.polyhedron() {
                super::boolean::polyhedra::imported::check(&solid)?;
            }
        }
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
        let recognized = match &self.recognized {
            Recognized::Construction(r) => {
                Recognized::Construction(Box::new(r.transform_with(r.operation, motion)?.0))
            }
            Recognized::Polyhedron => Recognized::Polyhedron,
            // S9e.4b.3a: read off the moved stored topology again.
            Recognized::Piece(p) => Recognized::Piece(Box::new(piece(
                &moved,
                tolerance,
                p.primitive.operation,
                p.form,
            )?)),
        };
        let (frame, start, end, bounds) = placed(&recognized, &moved);
        Ok(Solid {
            frame,
            start,
            end,
            construction: Construction::Imported(Box::new(Imported {
                recognized,
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

/// The solid's frame, heights and bounds: its construction's, the bounds
/// widened by the stored edges'; a polyhedron's the world's frame, its
/// stored edges' bounds and their heights.
fn placed(recognized: &Recognized, t: &Topology) -> (Frame3, f64, f64, crate::Bounds3) {
    match recognized {
        Recognized::Construction(r) => (r.frame, r.start, r.end, widened(r.bounds, t)),
        Recognized::Polyhedron => {
            let b = super::split::edge_bounds(t);
            (Frame3::xy(), b.min.z, b.max.z, b)
        }
        // S9e.4b.3a: the primitive's frame, the bounds its stored edges' and
        // its primitive's between the body's axial ends.
        Recognized::Piece(p) => {
            let r = &p.primitive;
            (r.frame, r.start, r.end, widened(piece_bounds(p, t), t))
        }
    }
}

/// Whether every face is a plane and every edge a line (S9e.4b.2).
fn planar(t: &Topology) -> bool {
    t.faces()
        .iter()
        .all(|f| matches!(f.surface, Surface::Plane(_)))
        && t.edges()
            .iter()
            .all(|e| matches!(e.curve, Curve3::LineSegment { .. }))
}

/// One solid region with one shell (no cavity), as `recognize` requires.
fn one_shell(t: &Topology) -> bool {
    let solids: Vec<_> = t
        .regions()
        .iter()
        .filter(|r| r.kind == RegionKind::Solid)
        .collect();
    solids.len() == 1 && solids[0].shells.len() == 1
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
pub(crate) fn outward(t: &Topology, fi: usize) -> Option<Vec3> {
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
    // S9e.4b.1: its arcs' ends are roundings, taken onto their circles by
    // the exact model.
    let profile = Profile::new(boundaries.remove(0), holes, tolerance)?.with_rounded_arcs();
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

// ------------------------------------------------------------------ plane pieces

/// Whether a stored topology is a plane piece's shape (S9e.4b.3a): one
/// solid region of one shell, one sphere, cylinder or cone face and plane
/// faces, no spline edge.
fn piece_shape(t: &Topology) -> bool {
    let curved: Vec<_> = t
        .faces()
        .iter()
        .filter(|f| !matches!(f.surface, Surface::Plane(_)))
        .collect();
    one_shell(t)
        && t.faces().len() >= 2
        && matches!(
            curved.as_slice(),
            [f] if matches!(
                f.surface,
                Surface::Sphere { .. } | Surface::Cylinder { .. } | Surface::Cone { .. }
            )
        )
        && !t
            .edges()
            .iter()
            .any(|e| matches!(e.curve, Curve3::BSpline(_)))
}

/// The axial range of the stored topology's edges about a frame (the
/// corners of their bounds).
fn axial_range(t: &Topology, frame: &Frame3) -> (f64, f64) {
    let b = super::split::edge_bounds(t);
    let (lo, hi) = (b.min.to_array(), b.max.to_array());
    let mut range = (f64::INFINITY, f64::NEG_INFINITY);
    for k in 0..8 {
        let c = Point3::new(
            if k & 1 == 0 { lo[0] } else { hi[0] },
            if k & 2 == 0 { lo[1] } else { hi[1] },
            if k & 4 == 0 { lo[2] } else { hi[2] },
        );
        let w = frame.coordinates(c)[2];
        range = (range.0.min(w), range.1.max(w));
    }
    range
}

/// A body of one curved face and planes that is no Boolean of its
/// primitive and its planes' hull (S9e.4b.3c.3b).
fn not_its_primitive() -> Error {
    Error::OutOfDomain(
        "an imported plane piece other than one Boolean of its primitive and its planes' hull \
         (S9e.4b.3c.3b)",
    )
}

/// The curved face of a plane piece's shape.
fn curved_face(t: &Topology) -> Option<usize> {
    t.faces()
        .iter()
        .position(|f| !matches!(f.surface, Surface::Plane(_)))
}

/// Whether the curved face's material lies outside its quadric (its outward
/// normal toward the axis or centre: a hole, groove or notch).
fn material_outside(t: &Topology, k: usize) -> bool {
    let curved = &t.faces()[k];
    let front = &t.regions()[t.shells()[curved.front.index()].region.index()];
    let solid = front.kind == RegionKind::Solid;
    (curved.sense.sign() < 0.0) == solid
}

/// The axial range of a face's edges about a frame: their ends' and points'
/// heights (16 a curve).
fn face_range(t: &Topology, k: usize, frame: &Frame3) -> (f64, f64) {
    let mut range = (f64::INFINITY, f64::NEG_INFINITY);
    for l in &t.faces()[k].loops {
        let fins = match &t.loops()[l.index()] {
            Loop::Edges { fins, .. } => fins,
            Loop::Vertex(v) => {
                let w = frame.coordinates(t.vertices()[v.index()].position)[2];
                range = (range.0.min(w), range.1.max(w));
                continue;
            }
        };
        for f in fins {
            let edge = &t.edges()[t.fins()[f.index()].edge.index()];
            for s in 0..=16 {
                let w = frame.coordinates(edge.curve.point(f64::from(s) / 16.0))[2];
                range = (range.0.min(w), range.1.max(w));
            }
        }
    }
    range
}

/// S9e.4b.3c.3: a cylinder's or a cone's caps: the plane faces normal to its
/// axis at an end of its curved face's axial range `(w0, w1)`, their outward
/// normal away from the range (the primitive's own end, the material inside
/// it) or, for a groove, into it (the primitive's end seen from outside: a
/// blind hole's floor), each its exact height along the axis rounded once.
fn caps(
    t: &Topology,
    k: usize,
    frame: &Frame3,
    (w0, w1): (f64, f64),
    form: Form,
    tolerance: Tolerance,
) -> Result<Vec<(usize, f64)>> {
    let local = Local::new(frame)?;
    let tol = tolerance.linear();
    let mut out = Vec::new();
    for (i, f) in t.faces().iter().enumerate() {
        let Surface::Plane(plane) = &f.surface else {
            continue;
        };
        if i == k || !parallel(plane.normal(), frame.normal()) {
            continue;
        }
        let h = rational_f64(&local.of(plane.origin())[2]);
        let up = outward(t, i).ok_or_else(general)?.dot(frame.normal()) > 0.0;
        let away = if (h - w0).abs() <= tol {
            !up
        } else if (h - w1).abs() <= tol {
            up
        } else {
            continue;
        };
        if away != (form == Form::Groove) {
            out.push((i, h));
        }
    }
    Ok(out)
}

/// A plane piece of a form (S9e.4b.3a; S9e.4b.3c.3): its primitive on the
/// curved face's stored frame and its plane faces' stored frames. A common
/// takes the primitive past the body's ends and every plane; the other forms
/// the primitive past its curved face's ends but at its caps (a cylinder's
/// or a cone's), and the planes but the caps', turned over for a bite.
fn piece(t: &Topology, tolerance: Tolerance, op: OperationId, form: Form) -> Result<Piece> {
    let k = curved_face(t).ok_or_else(general)?;
    // The curved face's material inside its quadric (its outward normal the
    // surface's own, away from the axis or centre) but for a groove.
    if material_outside(t, k) != (form == Form::Groove) {
        return Err(not_its_primitive());
    }
    let half = std::f64::consts::FRAC_PI_2;
    // The axial range: the body's (a common's) or the curved face's, past
    // either end by a quarter of it (at least of the radius) but at a cap.
    let range = |frame: &Frame3, radius: f64| -> Result<((f64, f64), Vec<usize>)> {
        let (w0, w1) = if form == Form::Common {
            axial_range(t, frame)
        } else {
            face_range(t, k, frame)
        };
        let margin = 0.25 * (w1 - w0).max(radius);
        let (mut lo, mut hi) = (w0 - margin, w1 + margin);
        let mut capped = Vec::new();
        if form != Form::Common {
            for (i, h) in caps(t, k, frame, (w0, w1), form, tolerance)? {
                if (h - w0).abs() <= (h - w1).abs() {
                    lo = h;
                } else {
                    hi = h;
                }
                capped.push(i);
            }
        }
        Ok(((lo, hi), capped))
    };
    let (primitive, capped) = match &t.faces()[k].surface {
        Surface::Sphere { frame, radius } => (
            Solid::build_sphere(op, *frame, *radius, -half, half, tolerance)?,
            Vec::new(),
        ),
        Surface::Cylinder { frame, radius } => {
            let ((lo, hi), capped) = range(frame, *radius)?;
            let profile = Profile::new(
                Boundary::circle(Point2::new(0.0, 0.0), *radius, tolerance)?,
                Vec::new(),
                tolerance,
            )?;
            (Solid::build(op, profile, *frame, lo, hi)?, capped)
        }
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => {
            // The radius `radius + w tan a`, not past the apex.
            let slope = half_angle.tan();
            let ((mut lo, mut hi), capped) = range(frame, radius.abs())?;
            let apex = -radius / slope;
            let (mut r_lo, mut r_hi) = (radius + lo * slope, radius + hi * slope);
            if slope > 0.0 && lo <= apex {
                lo = apex;
                r_lo = 0.0;
            }
            if slope < 0.0 && hi >= apex {
                hi = apex;
                r_hi = 0.0;
            }
            if r_lo < 0.0 || r_hi < 0.0 {
                return Err(general());
            }
            let base = frame.at(frame.point(Point2::new(0.0, 0.0), lo));
            (
                Solid::build_cone(op, base, r_lo, r_hi, hi - lo, tolerance)?,
                capped,
            )
        }
        _ => return Err(general()),
    };
    let mut planes = Vec::new();
    for (i, f) in t.faces().iter().enumerate() {
        if i == k || capped.contains(&i) {
            continue;
        }
        let Surface::Plane(frame) = &f.surface else {
            return Err(general());
        };
        let out = outward(t, i).ok_or_else(general)?;
        let into = frame.normal().dot(out) < 0.0;
        planes.push((*frame, into != (form == Form::Bite)));
    }
    Ok(Piece {
        primitive: Box::new(primitive),
        planes,
        form,
    })
}

/// A plane piece's bounds: its primitive's between the body's axial ends
/// (a cylinder's or a cone's), else its primitive's.
fn piece_bounds(p: &Piece, t: &Topology) -> crate::Bounds3 {
    let r = &p.primitive;
    let (frame, n) = (r.frame, r.frame.normal());
    let reach = |w: f64, radius: f64| {
        let c = frame.point(Point2::new(0.0, 0.0), w);
        let e = n.to_array().map(|x| radius * (1.0 - x * x).max(0.0).sqrt());
        (c, e)
    };
    let (w0, w1) = axial_range(t, &frame);
    let ends = match &r.construction {
        Construction::Prism(profile) => {
            let Some((_, radius)) = profile.outer().circle_geometry() else {
                return r.bounds;
            };
            [reach(w0, radius), reach(w1, radius)]
        }
        Construction::Cone { bottom, top, .. } => {
            let h = r.end - r.start;
            let at = |w: f64| bottom + (top - bottom) * ((w - r.start) / h).clamp(0.0, 1.0);
            [reach(w0, at(w0)), reach(w1, at(w1))]
        }
        _ => return r.bounds,
    };
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for (c, e) in ends {
        let c = c.to_array();
        for i in 0..3 {
            let pad = 4.0 * f64::EPSILON * (c[i].abs() + e[i]);
            lo[i] = lo[i].min(c[i] - e[i] - pad);
            hi[i] = hi[i].max(c[i] + e[i] + pad);
        }
    }
    crate::Bounds3 {
        min: Point3::new(lo[0], lo[1], lo[2]),
        max: Point3::new(hi[0], hi[1], hi[2]),
    }
}

#[cfg(test)]
mod tests {
    use crate::identity::OperationId;
    use crate::topology::Surface;
    use crate::{Boundary, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3};

    /// A stadium prism written, read back and imported, on frames whose
    /// normal is the first normalization of `(0, 3, 4)` or of `(0, 2, 3)`
    /// (the same bits on every host), given bit for bit. Normalized again,
    /// the first stays under macOS's `hypot` and turns by an ulp under
    /// glibc's (correctly rounded), the second the other way round. A cap
    /// built on the prism's normal normalized again was written an ulp off
    /// the walls' axis on one platform or the other, and the import,
    /// normalizing each stored frame once more, put the construction's
    /// cylinders (on its bottom cap's frame) an ulp off the stored ones: a
    /// Boolean's split walls failed the history's exactly parallel axes
    /// (the boolean target's `crash-26c72abf`, the second frame its
    /// counterpart on macOS). The caps take the prism's axes bit for bit,
    /// so the imported caps and cylinders share their axis on any host.
    #[test]
    fn imported_caps_and_walls_share_their_axis_on_either_platforms_frames() {
        let tol = Tolerance::default();
        let stadium = || {
            let (s, t) = (4.75, 2.0);
            let arc = |x: f64| Segment::Arc {
                center: Point2::new(x, 0.0),
                radius: t,
                ccw: true,
            };
            let points = [(0.0, -t), (s, -t), (s, t), (0.0, t)]
                .map(|(x, y)| Point2::new(x, y))
                .to_vec();
            let segments = vec![Segment::Line, arc(s), Segment::Line, arc(0.0)];
            let outer = Boundary::path(points, segments, tol).unwrap();
            Profile::new(outer, Vec::new(), tol).unwrap()
        };
        let axis = |y: u64, z: u64| Vec3::new(0.0, f64::from_bits(y), f64::from_bits(z));
        // Turned again by glibc's `hypot`, by macOS's.
        for n in [
            axis(0x3fe3_3333_3333_3333, 0x3fe9_9999_9999_999a),
            axis(0x3fe1_c01a_a03b_e895, 0x3fea_a027_f059_dce1),
        ] {
            let x = Vec3::new(1.0, 0.0, 0.0);
            let frame = Frame3::from_axes(Point3::new(1.0, -2.0, 0.5), x, n.cross(x), n);
            let (a, _) = Solid::extrude_with(OperationId(1), stadium(), frame, 0.0, 1.75).unwrap();
            let text = crate::occt_brep::write(a.topology(), a.resolution().linear()).unwrap();
            let doc = crate::occt_brep::read(&text).unwrap();
            let [solid] = <[_; 1]>::try_from(crate::occt_brep::import(&doc).solids).unwrap();
            let (stored, resolution) = (solid.result.unwrap(), solid.tolerance);
            let (imported, _) = Solid::imported_with(OperationId(11), stored, resolution).unwrap();
            // The crash's tool: a cone from its apex up the axis, crossing
            // both walls.
            let base = Frame3::new(frame.point(Point2::new(0.0, 0.0), 0.4375), n, x, tol).unwrap();
            let (cone, _) =
                Solid::cone_with(OperationId(2), base, 0.0, 2.4375, 0.875, tol).unwrap();
            for (out, history) in [
                imported.fuse(OperationId(3), &cone).unwrap(),
                imported.cut(OperationId(4), &cone).unwrap(),
                imported.common(OperationId(5), &cone).unwrap(),
            ] {
                let ins = [
                    imported.topology().entity_set(imported.resolution()),
                    cone.topology().entity_set(cone.resolution()),
                ];
                let outs: Vec<_> = out
                    .iter()
                    .map(|s| s.topology().entity_set(s.resolution()))
                    .collect();
                let issues = crate::history::check(&ins, &outs, &history);
                assert!(issues.is_empty(), "{n:?}: {issues:?}");
            }
            // The imported caps' normals and walls' axes one direction, up
            // to sign, bit for bit.
            let axes: Vec<Vec3> = imported
                .topology()
                .faces()
                .iter()
                .filter_map(|f| match &f.surface {
                    Surface::Cylinder { frame, .. } => Some(frame.normal()),
                    Surface::Plane(p) if p.normal().cross(n).length() < 1e-12 => Some(p.normal()),
                    _ => None,
                })
                .map(|m| if m.dot(n) > 0.0 { m } else { -m })
                .collect();
            assert_eq!(axes.len(), 4, "{n:?}");
            assert!(axes.iter().all(|m| *m == axes[0]), "{n:?}: {axes:?}");
        }
    }
}
