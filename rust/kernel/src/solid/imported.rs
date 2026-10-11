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
//! form tried in turn, the match the arbiter; S9e.4b.3c.3b: else a Boolean
//! tree of its primitive and several convex hulls of its planes (`Tree`), its
//! faces grouped by how its edges bend, each grouping tried in turn.
//! S9e.4b.4b.1: a body of several sphere, cylinder and cone faces whose plane
//! faces are all ends of their primitives is a Boolean chain of those
//! primitives (`primitives_piece`: widest first, each next one cut, in
//! common or fused by its material's side and how its edges with the
//! earlier ones bend; coaxial ones on one frame). S9e.4b.4b.2a: where some
//! plane faces are no primitive's ends, or its faces are tangent along an
//! edge, such a chain led by a prism leaf (`leaves`: S9e.4a's prism read off
//! a cap's loops whose edges' other faces are walls along its normal, its
//! walls' tangent joints its own), each leaf tried in turn. S9e.4b.4b.2b.1:
//! then each such chain with the other plane faces (`chain_piece`'s
//! `extra`): prisms of one cap (`bosses`: a boss read off its one cap, its
//! other end hidden in the part it stands on), flats (a part common the hull
//! of planes it meets along convex edges) and pockets (the hull of planes
//! joined along concave edges turned over, less the primitives standing in
//! it, cut from the chain). S9e.4b.4c.1: a polyhedron may hold cavities (one
//! solid region of several shells), and against curved faces its stored
//! triangles are a leaf of the curved engine (`boolean::curved::meshes`).
//! S9e.4b.4c.2a: a tree's pocket may hold pockets of its own (a tooth, a
//! post or an island standing in it), its hull the planes of its faces that
//! leave its boundary on one side, its other faces its own pockets, nested
//! at most three deep (`pocket`). S9e.4b.4b.2b.2a: or two prism leaves (a
//! box with a stadium through it, a stadium with a box across its top), a
//! leaf's walls another part splits on one surface (coaxial cylinders of one
//! radius), and a lone leaf with the other plane faces (`leaf_pairs`).
//! S9e.4b.4b.2b.2b.1: then the same chains over the leaves whose cap another
//! part or a pocket cuts at its rim too (`rim_leaves`: walls ending below it,
//! a hole a well on another part), a tooth standing in a pocket to a leaf's
//! cap running past it. S9e.4b.4b.2b.2b.2: then each leaf as the lead of a
//! chain with the parts flush with its walls (`flush_parts`: a prism within
//! its heights standing out of its walls, fused; a pocket open sideways,
//! cut), each part's loops completed beyond their runs along the lead's walls
//! (`complete_loop`). S9e.4b.4b.2b.2b.3a: then the stacks, two prisms on one
//! frame whose upper one stands on the lower one's cap (`stack_leaves`: two
//! plane faces facing apart never on one surface): each pair of leaves, the
//! second on the first's frame with its caps taken onto the first's heights
//! (`Levels`: the caps on one plane exactly), and each leaf with the prisms of
//! one cap read with their holes (a lower prism the upper one hides).
//! S9e.4b.4c.2b.1: in any of those chains, a component of the other plane
//! faces of both kinds is a pocket with pockets of its own (`nested_pocket`:
//! a post or an island standing in it, S9e.4b.4c.2a's `pocket` on the
//! chain's axes, nested at most three deep), any refusal of such a chain the
//! one its chain gave before. Every other body is S9e.4b's.
use super::{replayable, Construction, Context, Solid};
use crate::history::{History, Relation};
use crate::identity::{
    AlgorithmLevel, Derivation, EntityId, EntityKind, InputLabel, OperationId, OperationKind,
    Parent, Role,
};
use crate::profile::boolean::Op2;
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
use std::cell::RefCell;
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
    /// S9e.4b.4b.1: a body of several primitives' curved faces, its other
    /// primitives (the tree's leaves `Primitive(1)` on, in its order).
    pub(crate) others: Vec<Solid>,
    pub(crate) planes: Vec<(Frame3, bool)>,
    /// S9e.4b.3c.3: the Boolean of the primitive and the hull of the planes
    /// the body is.
    pub(crate) form: Form,
    /// S9e.4b.3c.3b: the Boolean tree of the primitive and its planes'
    /// hulls the body is (`Form::Tree`; `planes` then empty), and the
    /// choice of its faces' groups it was decomposed by.
    pub(crate) tree: Option<Tree>,
    pub(crate) choice: u32,
    /// S9e.4b.4b.1: a body of several primitives' chain (its other
    /// primitives, or S9e.4b.4b.2b.2a: a lone prism leaf with the other plane
    /// faces), read off its stored topology by `primitives_piece`.
    pub(crate) chain: bool,
    /// S9e.4b.4b.2b.2b.3a: a stack chain's (its leaves' caps on one plane),
    /// whose fused leaves' boundaries a point may lie on with the body's
    /// material all round it (`Tree::rank`'s probes).
    pub(crate) stack: bool,
    /// S9e.4b.4c.2b.1: a chain reading a component of the other plane faces
    /// of both kinds as a pocket with pockets of its own (any refusal of
    /// such a chain is one the body is not).
    pub(crate) nested: bool,
}

/// S9e.4b.3c.3b: a Boolean tree of a plane piece's primitive and convex
/// hulls of its planes (REVIEW_NOTES.md, "S9e.4b.3c.3b refined").
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Tree {
    /// A primitive: the piece's (0), or (S9e.4b.4b.1) another of a body of
    /// several primitives.
    Primitive(usize),
    /// The convex hull of planes: each plane face's stored frame and whether
    /// its normal `x * y` points into the hull.
    Hull(Vec<(Frame3, bool)>),
    /// A Boolean of two trees.
    Op(Op2, Box<Tree>, Box<Tree>),
}

impl Tree {
    /// Its Booleans.
    pub(crate) fn booleans(&self) -> usize {
        match self {
            Tree::Op(_, a, b) => 1 + a.booleans() + b.booleans(),
            _ => 0,
        }
    }

    /// A point's location by rank (inside 2, on the boundary 1, outside 0):
    /// each primitive's `primitives` at a point, a hull's sides within `tol`,
    /// a common the lesser, a fuse the greater (S9e.4b.4b.2b.2b.3a: in a
    /// stack's, `probe`, on both boundaries inside where no probe lies outside
    /// both), a cut the first's common the second's reversed.
    fn rank(
        &self,
        primitives: &dyn Fn(usize, Point3) -> u8,
        point: Point3,
        tol: f64,
        probe: bool,
    ) -> u8 {
        match self {
            Tree::Primitive(i) => primitives(*i, point),
            Tree::Hull(planes) => hull_rank(planes, point, tol),
            Tree::Op(op, a, b) => {
                let (x, y) = (
                    a.rank(primitives, point, tol, probe),
                    b.rank(primitives, point, tol, probe),
                );
                match op {
                    Op2::Common => x.min(y),
                    // S9e.4b.4b.2b.2b.3a: in a stack, on the boundary of
                    // both (its touching caps, or where one's face meets the
                    // other's), inside where no point `16 tol` off it along
                    // the axes and the cube's diagonals lies outside both:
                    // their directions lie within 35.3 degrees of every
                    // direction, so a wedge outside both of 90 degrees (a
                    // stack's walls on its caps) holds one.
                    Op2::Fuse if probe && x == 1 && y == 1 => {
                        let d = 16.0 * tol;
                        let s = d / 3f64.sqrt();
                        let mut off: Vec<Vec3> = Vec::new();
                        for k in 0..3 {
                            let mut v = [0.0; 3];
                            v[k] = d;
                            off.push(Vec3::new(v[0], v[1], v[2]));
                            off.push(Vec3::new(-v[0], -v[1], -v[2]));
                        }
                        for m in 0..8 {
                            let sign = |b: i32| if m >> b & 1 == 1 { -s } else { s };
                            off.push(Vec3::new(sign(0), sign(1), sign(2)));
                        }
                        let outside = off.iter().any(|v| {
                            let q = point + *v;
                            a.rank(primitives, q, tol, probe) == 0
                                && b.rank(primitives, q, tol, probe) == 0
                        });
                        if outside {
                            1
                        } else {
                            2
                        }
                    }
                    Op2::Fuse => x.max(y),
                    Op2::Cut => x.min(2 - y),
                }
            }
        }
    }
}

/// A point's rank in the hull of planes (each frame's plane, `into` its
/// normal pointing into the hull), within `tol`.
fn hull_rank(planes: &[(Frame3, bool)], point: Point3, tol: f64) -> u8 {
    let mut rank = 2;
    for (frame, into) in planes {
        let d = frame.coordinates(point)[2];
        let d = if *into { -d } else { d };
        if d > tol {
            return 0;
        }
        if d >= -tol {
            rank = 1;
        }
    }
    rank
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
    /// S9e.4b.3c.3b: a Boolean tree of the primitive and its planes' hulls.
    Tree,
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
        let h = hull_rank(&self.planes, point, tol);
        // Inside 2, on the boundary 1, outside 0: a common the least, a
        // fuse the most, a complement the reverse.
        let rank = |s: &Solid| -> Result<u8> {
            Ok(match s.classify(point)? {
                Location::Inside => 2,
                Location::Boundary => 1,
                Location::Outside => 0,
            })
        };
        let p = rank(&self.primitive)?;
        let r = match (self.form, &self.tree) {
            (Form::Common, _) => p.min(h),
            (Form::Groove, _) => h.min(2 - p),
            (Form::Bite, _) => p.min(2 - h),
            (Form::Boss, _) => p.max(h),
            (Form::Tree, Some(tree)) => {
                let mut ranks = vec![p];
                for s in &self.others {
                    ranks.push(rank(s)?);
                }
                // A probe's primitives classified there (an undecided one on
                // its boundary).
                let solids: Vec<&Solid> = std::iter::once(&*self.primitive)
                    .chain(&self.others)
                    .collect();
                let at = |i: usize, q: Point3| -> u8 {
                    if q == point {
                        return ranks[i];
                    }
                    match solids[i].classify(q) {
                        Ok(Location::Inside) => 2,
                        Ok(Location::Outside) => 0,
                        _ => 1,
                    }
                };
                tree.rank(&at, point, tol, self.stack)
            }
            (Form::Tree, None) => unreachable!("a tree's piece"),
        };
        Ok([Location::Outside, Location::Boundary, Location::Inside][usize::from(r)])
    }
}

fn general() -> Error {
    Error::OutOfDomain(
        "an imported solid other than a prism, a sphere, a cone, a torus or a Boolean of its \
         primitives (S9e.4b)",
    )
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
            // S9e.4b.4c.1: or with cavities (one solid region, its shells
            // the outer one and its voids').
            Err(e) if planar(&topology) => {
                if !one_region(&topology) {
                    return Err(match e {
                        Error::OutOfDomain(_) => {
                            Error::OutOfDomain("an imported polyhedron of several solids (S9e.4b)")
                        }
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
            // S9e.4b.4b.1: a body of several primitives' curved faces, a
            // Boolean chain of its primitives (S9e.4b.4b.2a: or led by a
            // prism leaf).
            Err(_) if primitives_shape(&topology) => {
                return Self::imported_primitives(context, topology, resolution);
            }
            Err(e) => return Err(e),
        };
        Self::imported_as(context, topology, resolution, recognized, names)
    }

    /// S9e.4b.4b.1: an imported body of several primitives as a Boolean
    /// chain of them; S9e.4b.4b.2a: where some plane faces are no
    /// primitive's ends or its faces are tangent along an edge, each chain
    /// led by a prism leaf in turn; S9e.4b.4b.2b.1: then each such chain
    /// with the other plane faces (prisms of one cap, flats, pockets);
    /// S9e.4b.4b.2b.2a: then each pair of leaves', alone and with them;
    /// S9e.4b.4b.2b.2b.1: then each leaf's and each pair's again over the
    /// leaves cut at their rim too; S9e.4b.4b.2b.2b.2: then each leaf's with
    /// the parts flush with its walls; S9e.4b.4b.2b.2b.3a: then each stack's
    /// (S9e.4b.4c.2b.1: in each, a component of the other plane faces of both
    /// kinds a pocket with pockets of its own, any refusal of such a chain the
    /// one it gave before); the first whose model matches the stored topology,
    /// else the first
    /// refusal an earlier chain's Booleans give, else the body's own tangency,
    /// else S9e.4b.4b.2b.2b.3b's.
    fn imported_primitives(
        context: &Context,
        topology: Topology,
        resolution: Tolerance,
    ) -> Result<(Self, History)> {
        let op = construction_operation(context.operation);
        // The chains whose models were built: a later choice reading the
        // same chain off the body (a leaf read through either of two faces
        // on its cap's plane, a stack's caps) is not built again. Its model
        // would be the same and refused the same, and a later choice's
        // refusal never counts for more than an earlier one's (the rules
        // below ignore as much or more as the choices rise), so the outcome
        // is the same without it.
        let tried: RefCell<Vec<Piece>> = RefCell::new(Vec::new());
        // S9e.4b.4c.2b.1: the chains reading a pocket with pockets of its own
        // whose models were refused so far (at `NESTED_CHAINS`, the later
        // chains refuse such a component as before).
        let refused = std::cell::Cell::new(0);
        let attempt = |choice: u32| -> Option<Result<(Self, History)>> {
            let nesting = refused.get() < NESTED_CHAINS;
            let p = match primitives_piece_with(&topology, resolution, op, choice, nesting) {
                Ok(p) => p,
                Err(e) => return Some(Err(e)),
            };
            if tried.borrow().iter().any(|q| same_chain(q, &p)) {
                return None;
            }
            tried.borrow_mut().push(p.clone());
            #[cfg(test)]
            CHAIN_MODELS.with(|c| c.set(c.get() + 1));
            // S9e.4b.4c.2b.1: any refusal of a chain reading a pocket with
            // pockets of its own is one the body is not: the refusal its
            // chain gave before (that component's), as S9e.4b.4b.2b.2b.1's
            // rim chains keep the body's.
            let nested = p.nested;
            Some(
                Self::imported_as(
                    context,
                    topology.clone(),
                    resolution,
                    Recognized::Piece(Box::new(p)),
                    BTreeMap::new(),
                )
                .map_err(|e| {
                    if nested {
                        refused.set(refused.get() + 1);
                        not_their_faces()
                    } else {
                        e
                    }
                }),
            )
        };
        let alone = match attempt(0).unwrap_or_else(|| Err(not_their_faces())) {
            Err(e) if e == not_their_faces() || e == tangent_body() => e,
            r => return r,
        };
        let (mut refusal, mut tangent) = (None, alone == tangent_body());
        let (choices, plain, rim) = primitives_choices(&topology, resolution);
        for choice in choices.into_iter().skip(1) {
            let Some(attempted) = attempt(choice) else {
                continue;
            };
            match attempted {
                Ok(done) => return Ok(done),
                // S9e.4b.4b.2b.1: a decomposition with the other plane faces
                // whose prisms part a joint tangent along its walls is one the
                // body is not.
                Err(e) if e == tangent_body() && choice >= plain => {}
                // S9e.4b.4b.2b.2b.1: nor any refusal of a chain over the
                // leaves cut at their rim (its Booleans', its tangency, its
                // mismatch), S9e.4b.4b.2b.2b.2: nor of a flush chain: the
                // body's refusal as before where none matches.
                Err(_) if choice >= rim => {}
                Err(e) if e == tangent_body() => tangent = true,
                Err(e) if e == not_their_faces() || e == general() => {}
                Err(e) => {
                    refusal.get_or_insert(e);
                }
            }
        }
        Err(refusal.unwrap_or_else(|| {
            if tangent {
                tangent_body()
            } else {
                not_their_faces()
            }
        }))
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
        // S9e.4b.3c.3b: then the Boolean trees its faces decompose into,
        // each choice of its groups in turn.
        let forms = Form::candidates(material_outside(&topology, k));
        let choices = tree_choices(&topology, resolution, op, k).unwrap_or(1);
        let tried = forms
            .iter()
            .map(|&f| (f, 0))
            .chain((0..choices).map(|c| (Form::Tree, c)));
        for (form, choice) in tried {
            let built = piece(&topology, resolution, op, form, choice).map_err(|x| {
                let tree = form == Form::Tree && matches!(x, Error::Degenerate(_));
                if x == not_its_primitive() || tree {
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
                // A later choice's refusal is a tree the body is not (its
                // first choice's is the decomposition's own).
                Err(_) if form == Form::Tree && choice > 0 => {}
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
            Recognized::Piece(p) if p.chain => Recognized::Piece(Box::new(primitives_piece(
                &moved,
                tolerance,
                p.primitive.operation,
                p.choice,
            )?)),
            Recognized::Piece(p) => Recognized::Piece(Box::new(piece(
                &moved,
                tolerance,
                p.primitive.operation,
                p.form,
                p.choice,
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
            let mut bounds = piece_bounds(r, t);
            // S9e.4b.4b.1: every primitive's.
            for o in &p.others {
                let b = piece_bounds(o, t);
                bounds = crate::Bounds3 {
                    min: Point3::new(
                        bounds.min.x.min(b.min.x),
                        bounds.min.y.min(b.min.y),
                        bounds.min.z.min(b.min.z),
                    ),
                    max: Point3::new(
                        bounds.max.x.max(b.max.x),
                        bounds.max.y.max(b.max.y),
                        bounds.max.z.max(b.max.z),
                    ),
                };
            }
            (r.frame, r.start, r.end, widened(bounds, t))
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

/// One solid region, of one shell or several (S9e.4b.4c.1: a polyhedron's
/// cavities).
fn one_region(t: &Topology) -> bool {
    t.regions()
        .iter()
        .filter(|r| r.kind == RegionKind::Solid)
        .count()
        == 1
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
    let loops: Vec<usize> = (0..faces[bottom].loops.len()).collect();
    prism_on(t, (bottom, top, up), &loops, (None, None), tolerance, op)
}

/// S9e.4a's prism between a bottom and a top cap (`up` the bottom's inward
/// normal): on the bottom's stored frame (turned over to `up`), its profile
/// the bottom's `loops` (indices in its loops, the outer first) rounded once
/// into it, its height the top's; S9e.4b.4b.2b: a prism of one cap (`open`
/// its walls) past its walls' far end by a quarter of their range.
fn prism_on(
    t: &Topology,
    (bottom, top, up): (usize, usize, Vec3),
    loops: &[usize],
    (open, along): (Option<&[usize]>, Option<&Frame3>),
    tolerance: Tolerance,
    op: OperationId,
) -> Result<Solid> {
    completed_on(
        t,
        (bottom, top, up),
        loops,
        (open, along),
        None,
        None,
        tolerance,
        op,
    )
}

/// S9e.4b.4b.2b.2b.3a: a stack chain's leaf: none for its first (on its own
/// cap's frame), else the first leaf's frame and its caps' heights on it,
/// which a later leaf parallel to it stands on. Every leaf of a stack chain
/// stands along its cylinder walls' stored axes (as a prism of one cap does).
struct Levels {
    first: Option<(Frame3, [f64; 2])>,
}

/// S9e.4b.4b.2b.2b.2: what a leaf's runs along the lead part are completed
/// by: the faces across them (`run`), the lead prism, and whether the
/// completion lies in its material (a fuse) or outside it (a cut).
struct Completion<'a> {
    run: &'a [bool],
    lead: &'a Solid,
    into: bool,
}

/// [`prism_on`], its loops' runs of edges across the lead part's walls
/// completed where given (`complete_loop`).
#[allow(clippy::too_many_arguments)]
fn completed_on(
    t: &Topology,
    (bottom, top, up): (usize, usize, Vec3),
    loops: &[usize],
    (open, along): (Option<&[usize]>, Option<&Frame3>),
    complete: Option<&Completion>,
    levels: Option<&Levels>,
    tolerance: Tolerance,
    op: OperationId,
) -> Result<Solid> {
    let faces = t.faces();
    let fins_of = fin_faces(t);
    let Surface::Plane(stored) = &faces[bottom].surface else {
        return Err(general());
    };
    // S9e.4b.4b.2b.1: a prism of one cap on the side its walls' stored axes
    // point to (its cylinders' frames and the model's one way), and a boss
    // parallel to the leaf on its axes, bit for bit.
    let against = open.is_some_and(|walls| {
        walls.iter().any(|&k| {
            matches!(&faces[k].surface, Surface::Cylinder { frame, .. } if frame.normal().dot(up) < 0.0)
        })
    });
    // S9e.4b.4b.2b.2b.3a: a stack chain's leaf likewise, its walls the faces
    // across its loops' edges (its model's cylinders then turning as the
    // stored ones its walls' meetings are stored over).
    let against = against
        || (levels.is_some()
            && open.is_none()
            && loops.iter().any(|&li| {
                let Some(l) = faces[bottom].loops.get(li) else {
                    return false;
                };
                let Loop::Edges { fins, .. } = &t.loops()[l.index()] else {
                    return false;
                };
                fins.iter().any(|f| {
                    let edge = &t.edges()[t.fins()[f.index()].edge.index()];
                    edge.fins.iter().filter(|x| *x != f).any(|x| {
                        matches!(&faces[fins_of[x.index()]].surface,
                            Surface::Cylinder { frame, .. } if frame.normal().dot(up) < 0.0)
                    })
                })
            }));
    let dir = if against { -up } else { up };
    let turned = |f: &Frame3| {
        if f.normal().dot(dir) > 0.0 {
            *f
        } else {
            f.flipped()
        }
    };
    let frame = match along {
        Some(a) if parallel(a.normal(), up) => turned(a).at(stored.origin()),
        _ => turned(stored),
    };
    // S9e.4b.4b.2b.2b.3a: a later leaf of a stack chain on the first leaf's
    // frame itself (its axes bit for bit, turned over where it faces the
    // other way, its origin too), its caps' heights there rounded once and
    // taken onto the first's within the resolution (caps on one stored plane,
    // or rounded apart, one plane exactly); none taken, no stack.
    let mut snapped = None;
    if let (
        Some(Levels {
            first: Some((lf, heights)),
        }),
        None,
        None,
    ) = (levels, open, complete)
    {
        if !parallel(lf.normal(), up) {
            return Err(not_their_chain());
        }
        let Surface::Plane(top_frame) = &faces[top].surface else {
            return Err(general());
        };
        let g = turned(lf);
        let sign = if g.normal().dot(lf.normal()) > 0.0 {
            1.0
        } else {
            -1.0
        };
        let lg = Local::new(&g)?;
        let tol = tolerance.linear();
        let snap = |h: f64| {
            heights
                .iter()
                .map(|x| sign * x)
                .find(|x| (x - h).abs() <= tol)
        };
        let (hb, ht) = (
            rational_f64(&lg.of(stored.origin())[2]),
            rational_f64(&lg.of(top_frame.origin())[2]),
        );
        match (snap(hb), snap(ht)) {
            (None, None) => return Err(not_their_chain()),
            (b, c) => snapped = Some((g, b.unwrap_or(hb), c.unwrap_or(ht))),
        }
    }
    let frame = snapped.map_or(frame, |s| s.0);
    let local = Local::new(&frame)?;
    // S9e.4b.4b.2b.2b.2: whether a point of the profile lies on the
    // completion's side of the lead prism's profile (in its material for a
    // fuse, outside it for a cut), within its tolerance.
    let beyond = |p: Point2| -> bool {
        let Some(c) = complete else { return false };
        let Construction::Prism(profile) = &c.lead.construction else {
            return false;
        };
        let w = c.lead.frame.coordinates(frame.point(p, 0.0));
        matches!(
            (profile.classify(Point2::new(w[0], w[1])), c.into),
            (Ok(crate::Location::Inside), true) | (Ok(crate::Location::Outside), false)
        )
    };
    let mut boundaries = Vec::new();
    for &li in loops {
        let l = &faces[bottom].loops[li];
        let Loop::Edges { fins, .. } = &t.loops()[l.index()] else {
            return Err(general());
        };
        let mut points = Vec::new();
        let mut segments = Vec::new();
        let mut runs = Vec::new();
        let mut walls_of: Vec<Option<usize>> = Vec::new();
        let mut ring = None;
        for f in fins {
            let fin = &t.fins()[f.index()];
            let edge = &t.edges()[fin.edge.index()];
            let forward = fin.sense == Orientation::Forward;
            let across = edge
                .fins
                .iter()
                .find(|x| *x != f)
                .map(|x| fins_of[x.index()]);
            runs.push(
                complete.is_some_and(|c| across.is_some_and(|k| k < c.run.len() && c.run[k])),
            );
            walls_of.push(across);
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
                    let about = c.normal().dot(dir) > 0.0;
                    Segment::Arc {
                        center: point2(&local.of(c.origin())),
                        radius: *radius,
                        ccw: turns == about,
                    }
                }
                _ => return Err(general()),
            });
        }
        // S9e.4b.4b.2b.2b.2: a ring along the lead part completed by its
        // circle offset to the completion's side (`ring_offset`), a path's
        // runs along it by `complete_loop`.
        if let (Some(_), Some((center, radius)), Some(true)) = (complete, ring, runs.first()) {
            points = ring_offset(center, radius, li > 0, &beyond)?;
            segments = vec![Segment::Line; points.len()];
            ring = None;
        } else if complete.is_some() && runs.iter().any(|&r| r) {
            let tol = tolerance.linear();
            let same = |i: usize, j: usize| match (walls_of[i], walls_of[j]) {
                (Some(a), Some(b)) => a == b || one_surface(t, a, b, tol),
                _ => false,
            };
            (points, segments) = complete_loop(points, segments, &runs, li > 0, &beyond, &same)?;
            // A loop bridged into one arc about its whole circle: a ring.
            if let ([_], [Segment::Arc { center, radius, .. }]) =
                (points.as_slice(), segments.as_slice())
            {
                ring = Some((*center, *radius));
            }
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
    let (start, end) = match open {
        Some(walls) => {
            let (lo, hi) = walls.iter().fold((0.0f64, 0.0f64), |(lo, hi), &k| {
                let (a, b) = face_range(t, k, &frame);
                (lo.min(a), hi.max(b))
            });
            if against {
                (1.25 * lo, 0.0)
            } else {
                (0.0, 1.25 * hi)
            }
        }
        // S9e.4b.4b.2b.2b.3a: the bottom's and the top's heights in either
        // order (a stack chain's leaf against its walls' axes).
        None => match snapped {
            Some((_, b, c)) => (b.min(c), b.max(c)),
            None => {
                let Surface::Plane(top_frame) = &faces[top].surface else {
                    return Err(general());
                };
                let h = rational_f64(&local.of(top_frame.origin())[2]);
                if against {
                    (h.min(0.0), h.max(0.0))
                } else {
                    (0.0, h)
                }
            }
        },
    };
    if end - start <= 0.0 {
        return Err(general());
    }
    Solid::build(op, profile, frame, start, end)
}

/// S9e.4b.4b.2b.2b.2: a ring of a leaf's loop along the lead part, completed
/// as `complete_loop` completes a loop whose every edge runs along it: its
/// circle offset by `delta` (the ring's diameter times `sqrt 2` over 64) away
/// from the leaf's cap (outward for an outer loop, inward for a hole), a
/// polygon within `delta / 4` of it, every vertex and side's middle on its
/// side of the lead part.
fn ring_offset(
    center: Point2,
    radius: f64,
    hole: bool,
    inside: &dyn Fn(Point2) -> bool,
) -> Result<Vec<Point2>> {
    let delta = 2.0 * radius * std::f64::consts::SQRT_2 / 64.0;
    let r = if hole { radius - delta } else { radius + delta };
    let step = 2.0 * (1.0 - delta / (4.0 * r)).acos();
    let count = (std::f64::consts::TAU / step).ceil().max(8.0) as usize;
    if count > COMPLETION_ARC {
        return Err(not_their_chain());
    }
    let at = |k: f64| {
        let th = std::f64::consts::TAU * k / count as f64;
        Point2::new(center.x + r * th.cos(), center.y + r * th.sin())
    };
    let points: Vec<Point2> = (0..count).map(|k| at(k as f64)).collect();
    let middles = (0..count).map(|k| {
        let (a, b) = (at(k as f64), at(k as f64 + 1.0));
        Point2::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0)
    });
    if !points.iter().copied().chain(middles).all(inside) {
        return Err(not_their_chain());
    }
    Ok(points)
}

/// S9e.4b.4b.2b.2b.2: the most points a completion takes along one arc.
const COMPLETION_ARC: usize = 64;

/// S9e.4b.4b.2b.2b.2: a leaf's loop (its points and segments, each segment's
/// edge across the lead part's wall or not, `runs`) with each maximal run of
/// edges across the lead part completed beyond it, on the side away from the
/// leaf's cap (into the lead part's material for a fuse, out of it for a
/// cut, as `inside` tells). A run between two walls of the leaf on one
/// surface (`same`: a line or circle the lead part parts), turning alike,
/// whose carrier between them lies on that side is bridged by it, the two one
/// segment (an arc bridged round to its own start, a stored vertex between
/// two faces of one circle, OCCT's seam, split there and at its antipode);
/// else the leaf's walls at the run's ends are extended along their lines or
/// circles to a depth `delta` past the run (`delta` the loop's extent over
/// 64; refused where a wall meets the run at a sine below 1/8), then lines
/// run through the run's joints and arcs offset by `delta` to that side (an
/// arc's polyline within `delta / 4` of its offset, at most `COMPLETION_ARC`
/// points; a joint folding back refused). A loop whose every edge runs along
/// the lead part is that offset alone. Every new point and segment's middle
/// must lie on its side of the lead part (`inside`), else the leaf is no
/// completion of it.
fn complete_loop(
    points: Vec<Point2>,
    segments: Vec<Segment>,
    runs: &[bool],
    hole: bool,
    inside: &dyn Fn(Point2) -> bool,
    same: &dyn Fn(usize, usize) -> bool,
) -> Result<(Vec<Point2>, Vec<Segment>)> {
    let n = points.len();
    if n == 0 || runs.len() != n {
        return Err(general());
    }
    type V2 = (f64, f64);
    let unit = |v: V2| -> V2 {
        // `sqrt`, correctly rounded on every host (`hypot` is not).
        let l = (v.0 * v.0 + v.1 * v.1).sqrt();
        (v.0 / l, v.1 / l)
    };
    let at = |j: usize| -> V2 {
        let p = points[j % n];
        (p.x, p.y)
    };
    let tau = std::f64::consts::TAU;
    // An arc's sweep from `a` to `b` about `c` (signed by its turn).
    let sweep = |c: V2, a: V2, b: V2, ccw: bool| -> f64 {
        let (ta, tb) = ((a.1 - c.1).atan2(a.0 - c.0), (b.1 - c.1).atan2(b.0 - c.0));
        let mut s = if ccw {
            (tb - ta).rem_euclid(tau)
        } else {
            -(ta - tb).rem_euclid(tau)
        };
        if s == 0.0 {
            s = if ccw { tau } else { -tau };
        }
        s
    };
    let arc_of = |j: usize| -> Option<(V2, f64, bool)> {
        match &segments[j % n] {
            Segment::Arc {
                center,
                radius,
                ccw,
            } => Some(((center.x, center.y), *radius, *ccw)),
            _ => None,
        }
    };
    if segments
        .iter()
        .any(|s| !matches!(s, Segment::Line | Segment::Arc { .. }))
    {
        return Err(general());
    }
    // The tangent of segment `j` at its start or end.
    let tangent = |j: usize, end: bool| -> V2 {
        let (a, b) = (at(j), at(j + 1));
        match arc_of(j) {
            None => unit((b.0 - a.0, b.1 - a.1)),
            Some((c, _, ccw)) => {
                let x = if end { b } else { a };
                let u = unit((x.0 - c.0, x.1 - c.1));
                if ccw {
                    (-u.1, u.0)
                } else {
                    (u.1, -u.0)
                }
            }
        }
    };
    // The loop's samples: its orientation and extent.
    let mut samples: Vec<V2> = Vec::new();
    for j in 0..n {
        let (a, b) = (at(j), at(j + 1));
        samples.push(a);
        if let Some((c, r, ccw)) = arc_of(j) {
            let (s, ta) = (sweep(c, a, b, ccw), (a.1 - c.1).atan2(a.0 - c.0));
            for k in 1..8 {
                let th = ta + s * f64::from(k) / 8.0;
                samples.push((c.0 + r * th.cos(), c.1 + r * th.sin()));
            }
        }
    }
    let m = samples.len();
    let area: f64 = (0..m)
        .map(|i| {
            let (a, b) = (samples[i], samples[(i + 1) % m]);
            a.0 * b.1 - a.1 * b.0
        })
        .sum();
    let (lo, hi) = samples.iter().fold(
        (
            (f64::INFINITY, f64::INFINITY),
            (f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
        |(lo, hi), p| {
            (
                (lo.0.min(p.0), lo.1.min(p.1)),
                (hi.0.max(p.0), hi.1.max(p.1)),
            )
        },
    );
    let delta = ((hi.0 - lo.0).powi(2) + (hi.1 - lo.1).powi(2)).sqrt() / 64.0;
    // The completion's side: right of the run where the leaf's cap lies on
    // its left.
    let sigma = if (area > 0.0) != hole { 1.0 } else { -1.0 };
    let side = |t: V2| -> V2 { (sigma * t.1, -sigma * t.0) };
    let refuse = not_their_chain;
    // A run joint's offset point.
    let joint = |k: usize| -> Result<V2> {
        let (n1, n2) = (side(tangent(k + n - 1, true)), side(tangent(k, false)));
        let den = 1.0 + n1.0 * n2.0 + n1.1 * n2.1;
        if den < 0.25 {
            return Err(refuse());
        }
        let x = at(k);
        Ok((
            x.0 + delta * (n1.0 + n2.0) / den,
            x.1 + delta * (n1.1 + n2.1) / den,
        ))
    };
    // An arc run segment's offset points strictly between its ends.
    let arc_points = |k: usize| -> Result<Vec<V2>> {
        let Some((c, r, ccw)) = arc_of(k) else {
            return Ok(Vec::new());
        };
        let (a, b) = (at(k), at(k + 1));
        let u = unit((a.0 - c.0, a.1 - c.1));
        let nn = side(tangent(k, false));
        let r2 = r + delta * (nn.0 * u.0 + nn.1 * u.1).signum();
        if r2 <= delta {
            return Err(refuse());
        }
        let s = sweep(c, a, b, ccw);
        let step = 2.0 * (1.0 - delta / (4.0 * r2)).acos();
        let count = (s.abs() / step).ceil().max(1.0);
        if count > COMPLETION_ARC as f64 {
            return Err(refuse());
        }
        let ta = (a.1 - c.1).atan2(a.0 - c.0);
        let count = count as usize;
        Ok((1..count)
            .map(|i| {
                let th = ta + s * i as f64 / count as f64;
                (c.0 + r2 * th.cos(), c.1 + r2 * th.sin())
            })
            .collect())
    };
    // A wall extended from a run's end `x` (the start of segment `j`, `back`,
    // or its end) by `delta` past the run along it.
    let extend = |j: usize, back: bool, normal: V2| -> Result<V2> {
        let x = if back { at(j) } else { at(j + 1) };
        let t = tangent(j, !back);
        let d = if back { (-t.0, -t.1) } else { t };
        let sine = d.0 * normal.0 + d.1 * normal.1;
        if sine < 0.125 {
            return Err(refuse());
        }
        let e = delta / sine;
        Ok(match arc_of(j) {
            None => (x.0 + e * d.0, x.1 + e * d.1),
            Some((c, r, ccw)) => {
                let phi = e / r;
                if phi > 1.0 {
                    return Err(refuse());
                }
                let turn = if ccw != back { phi } else { -phi };
                let th = (x.1 - c.1).atan2(x.0 - c.0) + turn;
                (c.0 + r * th.cos(), c.1 + r * th.sin())
            }
        })
    };
    let p2 = |v: V2| Point2::new(v.0, v.1);
    let mut checked: Vec<V2> = Vec::new();
    let (mut out_points, mut out_segments) = (Vec::new(), Vec::new());
    if runs.iter().all(|&r| r) {
        // The whole loop along the lead part: its offset alone.
        for k in 0..n {
            out_points.push(joint(k)?);
            out_points.extend(arc_points(k)?);
        }
        checked.extend(out_points.iter().copied());
        let c = out_points.len();
        for i in 0..c {
            let (a, b) = (out_points[i], out_points[(i + 1) % c]);
            checked.push(((a.0 + b.0) / 2.0, (a.1 + b.1) / 2.0));
        }
        if !checked.iter().all(|&v| inside(p2(v))) {
            return Err(refuse());
        }
        let segments = vec![Segment::Line; c];
        return Ok((out_points.into_iter().map(p2).collect(), segments));
    }
    // Each run's completion (from its in-segment's extended end to its
    // out-segment's extended start), by its first segment; or a run whose
    // in-segment and out-segment lie on one wall's line or circle, turning
    // alike, with that carrier between them on the completion's side,
    // bridged by it (the in-segment then ends at the out-segment's end, the
    // two one segment).
    let mut path_after: BTreeMap<usize, Vec<V2>> = BTreeMap::new();
    let mut start_of: BTreeMap<usize, V2> = BTreeMap::new();
    let mut bridged: Vec<usize> = Vec::new();
    for a in 0..n {
        if !runs[a] || runs[(a + n - 1) % n] {
            continue;
        }
        let mut b = a;
        while runs[(b + 1) % n] {
            b = (b + 1) % n;
        }
        let first = (a + n - 1) % n;
        let last = (b + 1) % n;
        let alike = first != last
            && same(first, last)
            && match (arc_of(first), arc_of(last)) {
                (None, None) => true,
                (Some((_, _, x)), Some((_, _, y))) => x == y,
                _ => false,
            };
        if alike {
            let (x_in, x_out) = (at(a), at(last));
            let between: Vec<V2> = match arc_of(first) {
                None => (1..8)
                    .map(|k| {
                        let f = f64::from(k) / 8.0;
                        (
                            x_in.0 + f * (x_out.0 - x_in.0),
                            x_in.1 + f * (x_out.1 - x_in.1),
                        )
                    })
                    .collect(),
                Some((c, r, ccw)) => {
                    let s = sweep(c, x_in, x_out, ccw);
                    let ta = (x_in.1 - c.1).atan2(x_in.0 - c.0);
                    (1..16)
                        .map(|k| {
                            let th = ta + s * f64::from(k) / 16.0;
                            (c.0 + r * th.cos(), c.1 + r * th.sin())
                        })
                        .collect()
                }
            };
            if between.iter().all(|&v| inside(p2(v))) {
                bridged.push(first);
                bridged.push(last);
                continue;
            }
        }
        let y_in = extend(first, false, side(tangent(a, false)))?;
        let y_out = extend(last, true, side(tangent(b, true)))?;
        let mut path = vec![y_in];
        let mut k = a;
        loop {
            path.extend(arc_points(k)?);
            if k == b {
                break;
            }
            k = (k + 1) % n;
            path.push(joint(k)?);
        }
        // The walls' extensions and the new lines (to the out-segment's new
        // start).
        let (x_in, x_out) = (at(a), at(last));
        checked.push(((x_in.0 + y_in.0) / 2.0, (x_in.1 + y_in.1) / 2.0));
        checked.push(((x_out.0 + y_out.0) / 2.0, (x_out.1 + y_out.1) / 2.0));
        let mut all = path.clone();
        all.push(y_out);
        checked.extend(all.iter().copied());
        for w in all.windows(2) {
            checked.push(((w[0].0 + w[1].0) / 2.0, (w[0].1 + w[1].1) / 2.0));
        }
        path_after.insert(first, path);
        start_of.insert(last, y_out);
    }
    if !checked.iter().all(|&v| inside(p2(v))) {
        return Err(refuse());
    }
    // A bridged in-segment runs on to its out-segment's end (through the
    // out-segments it continues in turn); a loop bridged all round is one
    // ring. The output starts at a segment no bridge continues.
    let next_of: BTreeMap<usize, usize> = bridged.chunks(2).map(|w| (w[0], w[1])).collect();
    let Some(j0) = (0..n).find(|&j| !runs[j] && !next_of.values().any(|&l| l == j)) else {
        let j = (0..n).find(|&j| !runs[j]).ok_or_else(general)?;
        return Ok((vec![p2(at(j))], vec![segments[j].clone()]));
    };
    let mut covered = vec![false; n];
    for i in 0..n {
        let j = (j0 + i) % n;
        if runs[j] || covered[j] {
            continue;
        }
        out_points.push(start_of.get(&j).copied().unwrap_or(at(j)));
        out_segments.push(segments[j].clone());
        let mut k = j;
        while let Some(&l) = next_of.get(&k) {
            covered[l] = true;
            k = l;
        }
        if let Some(path) = path_after.get(&k) {
            for &v in path {
                out_points.push(v);
                out_segments.push(Segment::Line);
            }
        }
    }
    // One arc bridged round to its own start (an arc OCCT split at its
    // circle's seam, the run between the halves' far ends): the circle split
    // at that start, which the stored topology keeps, and its antipode on the
    // completion's side.
    if let ([p], [Segment::Arc { center, .. }]) = (out_points.as_slice(), out_segments.as_slice()) {
        let opposite = (2.0 * center.x - p.0, 2.0 * center.y - p.1);
        if !inside(p2(opposite)) {
            return Err(refuse());
        }
        out_points.push(opposite);
        out_segments.push(out_segments[0].clone());
    }
    Ok((out_points.into_iter().map(p2).collect(), out_segments))
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

/// A body of one curved face and planes that is no Boolean tree of its
/// primitive and its planes' hulls (S9e.4b.3c.3b), its pockets nested at most
/// `DEPTH` deep in at most `BOOLEANS` Booleans (S9e.4b.4c.2a).
pub(crate) fn not_its_primitive() -> Error {
    Error::OutOfDomain(
        "an imported plane piece other than a Boolean tree of its primitive and its planes' hulls, \
         its pockets nested at most three deep (S9e.4b.4)",
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
fn piece(
    t: &Topology,
    tolerance: Tolerance,
    op: OperationId,
    form: Form,
    choice: u32,
) -> Result<Piece> {
    let k = curved_face(t).ok_or_else(general)?;
    if form == Form::Tree {
        return Ok(tree_piece(t, tolerance, op, k, choice)?.0);
    }
    // The curved face's material inside its quadric (its outward normal the
    // surface's own, away from the axis or centre) but for a groove.
    if material_outside(t, k) != (form == Form::Groove) {
        return Err(not_its_primitive());
    }
    let (primitive, capped) = primitive_of(t, tolerance, op, k, form)?;
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
        others: Vec::new(),
        planes,
        form,
        tree: None,
        choice: 0,
        chain: false,
        stack: false,
        nested: false,
    })
}

/// A plane piece's primitive (S9e.4b.3a; S9e.4b.3c.3) for a form other than
/// a tree (a tree's: a boss's, or a groove's outside the quadric), and its
/// caps' faces.
fn primitive_of(
    t: &Topology,
    tolerance: Tolerance,
    op: OperationId,
    k: usize,
    form: Form,
) -> Result<(Solid, Vec<usize>)> {
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
    Ok((primitive, capped))
}

// ------------------------------------------------------------------ Boolean trees

/// How an edge bends between its two faces (S9e.4b.3c.3b): convex where
/// the second face leans away from the first's material, concave where
/// toward it, flat within `1e-6` of neither (the faces tangent along it).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Bend {
    Convex,
    Concave,
    Flat,
}

impl Bend {
    fn reversed(self) -> Self {
        match self {
            Bend::Convex => Bend::Concave,
            Bend::Concave => Bend::Convex,
            Bend::Flat => Bend::Flat,
        }
    }
}

/// A face's oriented normal at a point of its surface (a plane's, a
/// sphere's, a cylinder's or a cone's), against the surface where the face
/// is reversed.
fn oriented_normal(t: &Topology, fi: usize, p: Point3) -> Option<Vec3> {
    let f = &t.faces()[fi];
    let n = match &f.surface {
        Surface::Plane(frame) => frame.normal(),
        s @ (Surface::Sphere { frame, .. }
        | Surface::Cylinder { frame, .. }
        | Surface::Cone { frame, .. }) => {
            let c = frame.coordinates(p);
            let u = c[1].atan2(c[0]);
            let v = c[2].atan2(c[0].hypot(c[1]));
            s.normal(Point2::new(u, v))
        }
        _ => return None,
    };
    Some(n * f.sense.sign())
}

/// A face's outward normal at a point (away from its solid region).
fn outward_at(t: &Topology, fi: usize, p: Point3) -> Option<Vec3> {
    let f = &t.faces()[fi];
    let n = oriented_normal(t, fi, p)?;
    let front = &t.regions()[t.shells()[f.front.index()].region.index()];
    Some(if front.kind == RegionKind::Solid {
        n
    } else {
        -n
    })
}

/// Each fin's face.
fn fin_faces(t: &Topology) -> Vec<usize> {
    let mut out = vec![usize::MAX; t.fins().len()];
    for (fi, f) in t.faces().iter().enumerate() {
        for l in &f.loops {
            if let Loop::Edges { fins, .. } = &t.loops()[l.index()] {
                for x in fins {
                    out[x.index()] = fi;
                }
            }
        }
    }
    out
}

/// Each edge of two faces with how it bends from the first to the second:
/// at its curve's middle point, the first face's outward normal against the
/// direction into the second face from the edge (the second's oriented
/// normal across its loop's tangent there, its loops counter-clockwise
/// about it), both ways alike.
fn bends(t: &Topology) -> Result<Vec<(usize, usize, Bend)>> {
    let faces = fin_faces(t);
    let mut out = Vec::new();
    for e in t.edges() {
        let [fa, fb] = e.fins.as_slice() else {
            return Err(not_its_primitive());
        };
        let (a, b) = (faces[fa.index()], faces[fb.index()]);
        if a == usize::MAX || b == usize::MAX || a == b {
            return Err(not_its_primitive());
        }
        let m = e.curve.point(0.5);
        let h = 1e-4;
        let along = e.curve.point(0.5 + h) - e.curve.point(0.5 - h);
        if along.length() == 0.0 {
            return Err(not_its_primitive());
        }
        let tangent = along * (1.0 / along.length());
        let side = |x: usize, y: usize, fin: crate::topology::FinId| -> Option<f64> {
            let n = oriented_normal(t, y, m)?;
            let d = n.cross(tangent * t.fins()[fin.index()].sense.sign());
            Some(outward_at(t, x, m)?.dot(d) / d.length().max(f64::MIN_POSITIVE))
        };
        let (Some(ab), Some(ba)) = (side(a, b, *fb), side(b, a, *fa)) else {
            return Err(not_its_primitive());
        };
        let bend = |v: f64| {
            if v.abs() <= 1e-6 {
                Bend::Flat
            } else if v < 0.0 {
                Bend::Convex
            } else {
                Bend::Concave
            }
        };
        if bend(ab) != bend(ba) {
            return Err(not_its_primitive());
        }
        out.push((a, b, bend(ab)));
    }
    Ok(out)
}

/// A plane face's stored frame and whether its normal `x * y` points into
/// the material (turned over where `over`).
fn plane_of(t: &Topology, i: usize, over: bool) -> Result<(Frame3, bool)> {
    let Surface::Plane(frame) = &t.faces()[i].surface else {
        return Err(general());
    };
    let out = outward(t, i).ok_or_else(general)?;
    let into = frame.normal().dot(out) < 0.0;
    Ok((*frame, into != over))
}

/// The root of a union-find set.
fn find(root: &mut [usize], x: usize) -> usize {
    let mut x = x;
    while root[x] != x {
        root[x] = root[root[x]];
        x = root[x];
    }
    x
}

/// A group's region: its hull's planes and its pockets (faces joined by
/// edges concave in the region's material, components of two faces or more
/// in the order of their first faces, turned over); `turned` for the
/// primitive's group outside the quadric (its faces turned over, its edges'
/// bends reversed). S9e.4b.4c.2a: each pocket its own region (`pocket`),
/// its hull less its own pockets in turn.
type Region = (Vec<(Frame3, bool)>, Vec<Tree>);

/// S9e.4b.4c.2a: how deep pockets nest (a pocket, its pocket, and that
/// one's).
const DEPTH: usize = 3;

/// S9e.4b.4c.2a: the most Booleans a body's tree takes.
const BOOLEANS: usize = 8;

fn group_region(
    t: &Topology,
    bends: &[(usize, usize, Bend)],
    members: &[usize],
    turned: bool,
    tolerance: Tolerance,
) -> Result<Region> {
    let n = t.faces().len();
    let mut root: Vec<usize> = (0..n).collect();
    let mut pocketed = vec![false; n];
    for &(a, b, bend) in bends {
        if !members.contains(&a) || !members.contains(&b) {
            continue;
        }
        let bend = if turned { bend.reversed() } else { bend };
        if bend == Bend::Flat {
            return Err(not_its_primitive());
        }
        if bend == Bend::Concave {
            let (x, y) = (find(&mut root, a), find(&mut root, b));
            root[x.max(y)] = x.min(y);
            pocketed[a] = true;
            pocketed[b] = true;
        }
    }
    let mut hull = Vec::new();
    let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &i in members {
        if pocketed[i] {
            let r = find(&mut root, i);
            components.entry(r).or_default().push(i);
        } else {
            hull.push(plane_of(t, i, turned)?);
        }
    }
    let pockets = components
        .into_values()
        .map(|faces| pocket(t, bends, &faces, !turned, tolerance, 1))
        .collect::<Result<Vec<_>>>()?;
    Ok((one_plane(hull, tolerance.linear()), pockets))
}

/// S9e.4b.4c.2a: a pocket's region (its material the region's turned over,
/// `turned`) at a depth of pockets: its hull the planes of its faces that
/// leave every boundary point of its faces (`boundary_points`) on its
/// material's side within the resolution, less each component of its other
/// faces (joined by their edges, in the order of their first faces) as a
/// pocket of its own, turned over again; at most `DEPTH` deep. A convex
/// pocket is its hull alone (S9e.4b.3c.3b's); a pocket none of whose faces
/// is on its hull, or holding pockets past the depth, is refused. Not by the
/// bends alone: a U island's notch's back wall meets only the notch's sides
/// and the pocket's floor, along edges convex in the pocket's material, yet
/// its plane cuts the pocket in two; by the planes' sides it is the
/// island's, whose region then has the notch as its pocket.
fn pocket(
    t: &Topology,
    bends: &[(usize, usize, Bend)],
    faces: &[usize],
    turned: bool,
    tolerance: Tolerance,
    depth: usize,
) -> Result<Tree> {
    let tol = tolerance.linear();
    let n = t.faces().len();
    let corners: Vec<Point3> = faces.iter().flat_map(|&i| boundary_points(t, i)).collect();
    let mut hull = Vec::new();
    let mut rest = Vec::new();
    for &i in faces {
        let (frame, into) = plane_of(t, i, turned)?;
        let apart = corners.iter().any(|&p| {
            let d = frame.coordinates(p)[2];
            if into {
                d < -tol
            } else {
                d > tol
            }
        });
        if apart {
            rest.push(i);
        } else {
            hull.push((frame, into));
        }
    }
    // A pocket of its own pockets alone bounds nothing.
    if hull.is_empty() {
        return Err(not_its_primitive());
    }
    let mut root: Vec<usize> = (0..n).collect();
    for &(a, b, _) in bends {
        if rest.contains(&a) && rest.contains(&b) {
            let (x, y) = (find(&mut root, a), find(&mut root, b));
            root[x.max(y)] = x.min(y);
        }
    }
    let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &i in &rest {
        let r = find(&mut root, i);
        components.entry(r).or_default().push(i);
    }
    // A pocket's own pockets past the depth: a tree too deep.
    if !components.is_empty() && depth >= DEPTH {
        return Err(not_its_primitive());
    }
    let mut tree = Tree::Hull(one_plane(hull, tol));
    for own in components.into_values() {
        let x = pocket(t, bends, &own, !turned, tolerance, depth + 1)?;
        tree = Tree::Op(Op2::Cut, Box::new(tree), Box::new(x));
    }
    Ok(tree)
}

/// A face's boundary points: its edges' ends and 15 points along each
/// edge's stored curve (`face_range`'s), a vertex loop's vertex.
fn boundary_points(t: &Topology, fi: usize) -> Vec<Point3> {
    let mut out = Vec::new();
    for l in &t.faces()[fi].loops {
        match &t.loops()[l.index()] {
            Loop::Edges { fins, .. } => {
                for f in fins {
                    let e = &t.edges()[t.fins()[f.index()].edge.index()];
                    for v in [e.start, e.end].into_iter().flatten() {
                        out.push(t.vertices()[v.index()].position);
                    }
                    for k in 1..16 {
                        out.push(e.curve.point(f64::from(k) / 16.0));
                    }
                }
            }
            Loop::Vertex(v) => out.push(t.vertices()[v.index()].position),
        }
    }
    out
}

/// Planes of faces within the resolution of one plane (two walls of a U
/// prism on one plane, their stored frames rounded apart) on the first's
/// frame, so the hull takes them as one plane.
fn one_plane(planes: Vec<(Frame3, bool)>, tol: f64) -> Vec<(Frame3, bool)> {
    let mut out: Vec<(Frame3, bool)> = Vec::new();
    for (frame, into) in planes {
        let alike = out.iter().find(|(f, _)| {
            parallel(f.normal(), frame.normal()) && f.coordinates(frame.origin())[2].abs() <= tol
        });
        out.push(match alike {
            Some((f, _)) => (*f, into != (f.normal().dot(frame.normal()) < 0.0)),
            None => (frame, into),
        });
    }
    out
}

/// Each plane face's group (S9e.4b.3c.3b): the primitive's (0: sharing an
/// edge with the curved face convex inside the quadric, concave outside; the
/// caps) or the other (1). An edge of the kind no seam between the groups
/// takes (convex inside the quadric, concave outside: within a hull, or
/// between a hull and its pocket) keeps its faces in one group, every other
/// face joining through such edges; a face reached only across the other
/// kind (a pocket's own edges, or the seam where the groups are fused or
/// cut) lies with the faces it reaches through the first kind, those
/// components each in the other group than the faces they meet (bit clear in
/// `choice`) or the same (bit set), in the order of their first faces.
/// Returns the groups and the number of choices.
fn groups(
    t: &Topology,
    bends: &[(usize, usize, Bend)],
    k: usize,
    capped: &[usize],
    outside: bool,
    choice: u32,
) -> Result<(Vec<Option<usize>>, u32)> {
    let n = t.faces().len();
    let mut group: Vec<Option<usize>> = vec![None; n];
    for &c in capped {
        group[c] = Some(0);
    }
    let mut seeds: Vec<Option<usize>> = vec![None; n];
    for &(a, b, bend) in bends {
        let plane = match (a == k, b == k) {
            (true, false) => b,
            (false, true) => a,
            _ => continue,
        };
        if capped.contains(&plane) {
            continue;
        }
        let g = match bend {
            Bend::Flat => {
                return Err(Error::Degenerate(
                    "an imported plane piece whose curved face is tangent to its plane faces",
                ))
            }
            Bend::Convex => usize::from(outside),
            Bend::Concave => usize::from(!outside),
        };
        match seeds[plane] {
            Some(h) if h != g => return Err(not_its_primitive()),
            _ => seeds[plane] = Some(g),
        }
    }
    for (i, s) in seeds.iter().enumerate() {
        if s.is_some() {
            group[i] = *s;
        }
    }
    // The kind of edge that keeps a face of a group in it: within a hull,
    // or between a hull and its pocket (convex, but concave for the
    // primitive's group outside the quadric, its faces turned over).
    let keeping = |g: usize| {
        if g == 0 && outside {
            Bend::Concave
        } else {
            Bend::Convex
        }
    };
    let plane_edges = || {
        bends
            .iter()
            .filter(move |(a, b, _)| *a != k && *b != k)
            .copied()
    };
    // Through keeping edges, in rounds; both groups at once refused.
    let spread = |group: &mut Vec<Option<usize>>| -> Result<()> {
        loop {
            let before = group.clone();
            let mut changed = false;
            for (i, slot) in group.iter_mut().enumerate() {
                if i == k || slot.is_some() {
                    continue;
                }
                let mut reached = [false; 2];
                for (a, b, bend) in plane_edges() {
                    let other = if a == i {
                        b
                    } else if b == i {
                        a
                    } else {
                        continue;
                    };
                    if let Some(g) = before[other] {
                        if bend == keeping(g) {
                            reached[g] = true;
                        }
                    }
                }
                match reached {
                    [true, true] => return Err(not_its_primitive()),
                    [false, false] => {}
                    _ => {
                        *slot = Some(usize::from(reached[1]));
                        changed = true;
                    }
                }
            }
            if !changed {
                return Ok(());
            }
        }
    };
    spread(&mut group)?;
    // The faces left: components through convex edges, each meeting
    // assigned faces across the other kind.
    let mut root: Vec<usize> = (0..n).collect();
    for (a, b, bend) in plane_edges() {
        if bend == Bend::Convex && group[a].is_none() && group[b].is_none() {
            let (x, y) = (find(&mut root, a), find(&mut root, b));
            root[x.max(y)] = x.min(y);
        }
    }
    let mut components: Vec<usize> = Vec::new();
    for (i, g) in group.iter().enumerate() {
        if i != k && g.is_none() {
            let r = find(&mut root, i);
            if !components.contains(&r) {
                components.push(r);
            }
        }
    }
    if components.len() > 3 {
        return Err(not_its_primitive());
    }
    for (bit, &c) in components.iter().enumerate() {
        let members: Vec<usize> = (0..n)
            .filter(|&i| i != k && group[i].is_none() && find(&mut root, i) == c)
            .collect();
        // The first assigned face it meets, and its group.
        let met = plane_edges().find_map(|(a, b, _)| {
            match (members.contains(&a), members.contains(&b)) {
                (true, false) => group[b].map(|g| (b, g)),
                (false, true) => group[a].map(|g| (a, g)),
                _ => None,
            }
        });
        let Some((m, g)) = met else {
            return Err(not_its_primitive());
        };
        // Met from the other group or from a face of a pocket (one meeting a
        // face of its own group across a pocket's edge), the pocket first;
        // else the seam first.
        let pocket = if g == 0 && outside {
            Bend::Convex
        } else {
            Bend::Concave
        };
        let in_pocket = plane_edges().any(|(a, b, bend)| {
            bend == pocket && ((a == m && group[b] == Some(g)) || (b == m && group[a] == Some(g)))
        });
        let same = (g == 1 || in_pocket) != (choice >> bit & 1 == 1);
        for i in members {
            group[i] = Some(if same { g } else { 1 - g });
        }
    }
    Ok((group, 1 << components.len()))
}

/// S9e.4b.3c.3b: how many trees a body's faces decompose into (its
/// groupings, inside the quadric each two ways).
fn tree_choices(t: &Topology, tolerance: Tolerance, op: OperationId, k: usize) -> Result<u32> {
    let outside = material_outside(t, k);
    let form = if outside { Form::Groove } else { Form::Boss };
    let (_, capped) = primitive_of(t, tolerance, op, k, form)?;
    let (_, grouped) = groups(t, &bends(t)?, k, &capped, outside, 0)?;
    Ok(grouped * if outside { 1 } else { 2 })
}

/// S9e.4b.3c.3b: an imported plane piece as a Boolean tree of its primitive
/// and its planes' convex hulls (REVIEW_NOTES.md, "S9e.4b.3c.3b refined"):
/// its plane faces in the primitive's group (trimming it) and the other
/// (joined to it, or cutting it outside the quadric), each group's region
/// its hull less its pockets; the primitive common its group's region fused
/// with the other's, or the other's less it.
fn tree_piece(
    t: &Topology,
    tolerance: Tolerance,
    op: OperationId,
    k: usize,
    choice: u32,
) -> Result<(Piece, u32)> {
    let outside = material_outside(t, k);
    let form = if outside { Form::Groove } else { Form::Boss };
    let (primitive, capped) = primitive_of(t, tolerance, op, k, form)?;
    let n = t.faces().len();
    let bends = bends(t)?;
    let (group, grouped) = groups(t, &bends, k, &capped, outside, choice)?;
    // Inside the quadric, the union two ways (the second choice of each
    // grouping): the primitive common its group's region fused with the
    // other's, or the primitive common both regions' union (the other's
    // hidden faces within the first, none of them meeting the curved face:
    // a box within a ball's half fused with it), the latter first where no
    // face of the other group meets the curved face.
    let variants = if outside { 1 } else { 2 };
    let touches = bends
        .iter()
        .any(|&(a, b, _)| (a == k && group[b] == Some(1)) || (b == k && group[a] == Some(1)));
    let within = (choice / grouped == 0) != touches;
    let members = |g: usize| -> Vec<usize> {
        (0..n)
            .filter(|&i| i != k && group[i] == Some(g) && !capped.contains(&i))
            .collect()
    };
    let (p_hull, p_pockets) = group_region(t, &bends, &members(0), outside, tolerance)?;
    let (o_hull, o_pockets) = group_region(t, &bends, &members(1), false, tolerance)?;
    let cut = |x: Tree, pockets: Vec<Tree>| {
        pockets
            .into_iter()
            .fold(x, |x, k| Tree::Op(Op2::Cut, Box::new(x), Box::new(k)))
    };
    let others = !o_hull.is_empty() || !o_pockets.is_empty();
    // A group of pockets alone outside the primitive bounds nothing.
    if others && o_hull.is_empty() {
        return Err(not_its_primitive());
    }
    let tree = if !outside && others && within {
        if p_hull.is_empty() {
            return Err(not_its_primitive());
        }
        let x = cut(Tree::Hull(p_hull), p_pockets);
        let y = cut(Tree::Hull(o_hull), o_pockets);
        let union = Tree::Op(Op2::Fuse, Box::new(x), Box::new(y));
        Tree::Op(Op2::Common, Box::new(Tree::Primitive(0)), Box::new(union))
    } else {
        // The primitive common its group's hull, less its pockets.
        let mut x = Tree::Primitive(0);
        if !p_hull.is_empty() {
            x = Tree::Op(Op2::Common, Box::new(x), Box::new(Tree::Hull(p_hull)));
        }
        let x = cut(x, p_pockets);
        if !others {
            // Outside the quadric the primitive's alone bounds nothing;
            // inside, the union's second way is its first.
            if outside || choice / grouped == 1 {
                return Err(not_its_primitive());
            }
            x
        } else {
            let y = cut(Tree::Hull(o_hull), o_pockets);
            if outside {
                Tree::Op(Op2::Cut, Box::new(y), Box::new(x))
            } else {
                Tree::Op(Op2::Fuse, Box::new(x), Box::new(y))
            }
        }
    };
    let choices = grouped * variants;
    if !(1..=BOOLEANS).contains(&tree.booleans()) {
        return Err(not_its_primitive());
    }
    let piece = Piece {
        primitive: Box::new(primitive),
        others: Vec::new(),
        planes: Vec::new(),
        form: Form::Tree,
        tree: Some(tree),
        choice,
        chain: false,
        stack: false,
        nested: false,
    };
    Ok((piece, choices))
}

// ------------------------------------------------------------------ several primitives

/// S9e.4b.4b.1: the curved faces grouped by their stored surfaces (faces on
/// one surface one primitive), in face order.
fn surfaces(t: &Topology) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for (i, f) in t.faces().iter().enumerate() {
        if matches!(f.surface, Surface::Plane(_)) {
            continue;
        }
        match out
            .iter_mut()
            .find(|g| t.faces()[g[0]].surface == f.surface)
        {
            Some(g) => g.push(i),
            None => out.push(vec![i]),
        }
    }
    out
}

/// Whether a stored topology is a body of several primitives' shape
/// (S9e.4b.4b.1): one solid region of one shell, sphere, cylinder, cone and
/// plane faces, the curved ones on two surfaces or more, no spline edge.
fn primitives_shape(t: &Topology) -> bool {
    one_shell(t)
        && t.faces().iter().all(|f| {
            matches!(
                f.surface,
                Surface::Plane(_)
                    | Surface::Sphere { .. }
                    | Surface::Cylinder { .. }
                    | Surface::Cone { .. }
            )
        })
        && surfaces(t).len() >= 2
        && !t
            .edges()
            .iter()
            .any(|e| matches!(e.curve, Curve3::BSpline(_)))
}

/// A body of several primitives no chain of them matches.
pub(crate) fn not_their_chain() -> Error {
    Error::OutOfDomain(
        "an imported body of several primitives other than a Boolean chain of them (S9e.4b.4c)",
    )
}

/// S9e.4b.4b.1: a primitive's stored axis (a cylinder's or a cone's, a
/// sphere's frame's normal through its centre).
fn axis_of(t: &Topology, k: usize) -> Option<Frame3> {
    match &t.faces()[k].surface {
        Surface::Sphere { frame, .. }
        | Surface::Cylinder { frame, .. }
        | Surface::Cone { frame, .. } => Some(*frame),
        _ => None,
    }
}

/// Whether an axis lies on a frame's within the resolution: parallel
/// within `ALIGN`, its origin within `tol` of the frame's axis line.
fn coaxial(reference: &Frame3, axis: &Frame3, tol: f64) -> bool {
    let c = reference.coordinates(axis.origin());
    parallel(reference.normal(), axis.normal()) && c[0].hypot(c[1]) <= tol
}

/// A primitive's heights `(lo, hi)`, its faces' axial range and its caps'
/// faces.
type Span = ((f64, f64), (f64, f64), Vec<usize>);

/// S9e.4b.4b.1: a primitive of a body of several (its faces on one surface):
/// a whole sphere, or a cylinder or a cone over its faces' axial range,
/// past each end by a quarter of it (at least of the radius) but at a cap
/// (`caps`: a plane face normal to its axis at that end facing away from
/// the range, or into it where its material lies outside its quadric),
/// clamped at a cone's apex; on its stored frame, or on `reference` (an
/// earlier primitive's frame whose axis its own lies on within the
/// resolution: coaxial primitives on one axis, their caps' heights rounded
/// alike); its caps' faces and its reach from its axis (a sphere's radius,
/// a cylinder's, a cone's widest end's over its faces' range).
/// S9e.4b.4b.2b.2b.1: no face of `skip` is a cap (a prism's: a tooth standing
/// in a pocket to a leaf's cap runs past it, the cap the leaf's).
#[allow(clippy::too_many_arguments)]
fn primitive_over(
    t: &Topology,
    tolerance: Tolerance,
    op: OperationId,
    faces: &[usize],
    outside: bool,
    reference: Option<&Frame3>,
    past: &[crate::Bounds3],
    skip: &[usize],
) -> Result<(Solid, Vec<usize>, f64)> {
    let half = std::f64::consts::FRAC_PI_2;
    let k = faces[0];
    let form = if outside { Form::Groove } else { Form::Boss };
    let range = |frame: &Frame3, radius: f64| -> Result<Span> {
        let (w0, w1) = faces
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |r, &f| {
                let (a, b) = face_range(t, f, frame);
                (r.0.min(a), r.1.max(b))
            });
        let margin = 0.25 * (w1 - w0).max(radius);
        // Past the other primitives' bounds along the axis, where given.
        let (p0, p1) = past.iter().fold((w0, w1), |r, b| {
            let (a, c) = (b.min.to_array(), b.max.to_array());
            (0..8).fold(r, |r, k| {
                let corner = Point3::new(
                    if k & 1 == 0 { a[0] } else { c[0] },
                    if k & 2 == 0 { a[1] } else { c[1] },
                    if k & 4 == 0 { a[2] } else { c[2] },
                );
                let w = frame.coordinates(corner)[2];
                (r.0.min(w), r.1.max(w))
            })
        });
        let (mut lo, mut hi) = (p0 - margin, p1 + margin);
        let mut capped = Vec::new();
        for (i, h) in caps(t, k, frame, (w0, w1), form, tolerance)? {
            if skip.contains(&i) {
                continue;
            }
            if (h - w0).abs() <= (h - w1).abs() {
                lo = h;
            } else {
                hi = h;
            }
            capped.push(i);
        }
        Ok(((lo, hi), (w0, w1), capped))
    };
    match &t.faces()[k].surface {
        Surface::Sphere { frame, radius } => Ok((
            Solid::build_sphere(op, *frame, *radius, -half, half, tolerance)?,
            Vec::new(),
            *radius,
        )),
        Surface::Cylinder { frame, radius } => {
            let frame = reference.unwrap_or(frame);
            let ((lo, hi), _, capped) = range(frame, *radius)?;
            let profile = Profile::new(
                Boundary::circle(Point2::new(0.0, 0.0), *radius, tolerance)?,
                Vec::new(),
                tolerance,
            )?;
            Ok((Solid::build(op, profile, *frame, lo, hi)?, capped, *radius))
        }
        Surface::Cone {
            frame: stored,
            radius,
            half_angle,
        } => {
            // The radius `radius + (w - c) slope` at height `w` on the
            // frame, `c` the stored origin's height (zero on its own).
            let frame = reference.unwrap_or(stored);
            let along = if frame.normal().dot(stored.normal()) < 0.0 {
                -1.0
            } else {
                1.0
            };
            let c = match reference {
                Some(r) => rational_f64(&Local::new(r)?.of(stored.origin())[2]),
                None => 0.0,
            };
            let slope = along * half_angle.tan();
            let at = |w: f64| radius + (w - c) * slope;
            let ((mut lo, mut hi), (w0, w1), capped) = range(frame, radius.abs())?;
            let reach = at(w0).abs().max(at(w1).abs());
            let apex = c - radius / slope;
            let (mut r_lo, mut r_hi) = (at(lo), at(hi));
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
            Ok((
                Solid::build_cone(op, base, r_lo, r_hi, hi - lo, tolerance)?,
                capped,
                reach,
            ))
        }
        _ => Err(general()),
    }
}

/// S9e.4b.4b.2a: a body of several primitives whose plane faces are not all
/// their ends or a prism leaf's faces; S9e.4b.4b.2b.1: nor a prism's of one
/// cap, a flat's or a pocket's; S9e.4b.4b.2b.2a: nor two leaves';
/// S9e.4b.4b.2b.2b.1: nor leaves' cut at their rim; S9e.4b.4b.2b.2b.2: nor
/// parts flush with a leaf's walls; S9e.4b.4b.2b.2b.3a: nor a stack's
/// (S9e.4b.4b.2b.2b.3b); S9e.4b.4c.2b.1: nor a pocket's with its pockets
/// nested at most three deep (S9e.4b.4c.2b.2: a flat with pockets of its own,
/// a nested pocket holding a primitive).
fn not_their_faces() -> Error {
    Error::OutOfDomain(
        "an imported body of several primitives with plane faces other than their ends, one or \
         two prisms', a flat or a pocket, its pockets nested at most three deep \
         (S9e.4b.4b.2b.2b.3b)",
    )
}

/// A body of several primitives whose own faces are tangent along an edge.
fn tangent_body() -> Error {
    Error::Degenerate(
        "an imported body of several primitives whose faces are tangent along an edge",
    )
}

/// S9e.4b.4b.2a: a prism leaf of a body of several primitives: its caps
/// (bottom, top and the bottom's inward normal), the bottom's loops its
/// profile is read off (indices, the outer first) and its faces (the caps,
/// the walls across those loops' edges and every face on one of their
/// surfaces).
#[derive(Debug, Clone)]
struct Leaf {
    caps: (usize, usize, Vec3),
    loops: Vec<usize>,
    faces: Vec<bool>,
    /// S9e.4b.4b.2b.1: a prism of one cap, its walls (its height a quarter
    /// past their far end).
    open: Option<Vec<usize>>,
    /// S9e.4b.4b.2b.2b.2: the faces across its loops' runs along the lead
    /// part (completed beyond them; empty where none), and whether it is a
    /// void cut from the chain (a pocket open sideways) rather than fused.
    run: Vec<bool>,
    void: bool,
}

/// Whether a face is a wall along `up`: a plane whose outward normal is
/// perpendicular to it, or a cylinder whose axis is parallel to it.
fn wall(t: &Topology, k: usize, up: Vec3) -> bool {
    match &t.faces()[k].surface {
        Surface::Plane(_) => outward(t, k).is_some_and(|m| perpendicular(m, up)),
        Surface::Cylinder { frame, .. } => parallel(frame.normal(), up),
        _ => false,
    }
}

/// S9e.4b.4b.2a: whether an edge of two faces is a prism's walls' joint: a
/// line along a cylinder's stored axis between two walls along it (a fillet
/// tangent to a wall or another fillet).
fn joint(t: &Topology, edge: usize, a: usize, b: usize) -> bool {
    let e = &t.edges()[edge];
    if !matches!(e.curve, Curve3::LineSegment { .. }) {
        return false;
    }
    let along = e.curve.point(1.0) - e.curve.point(0.0);
    [a, b].iter().any(|&k| match &t.faces()[k].surface {
        Surface::Cylinder { frame, .. } => {
            let d = frame.normal();
            parallel(along, d) && wall(t, a, d) && wall(t, b, d)
        }
        _ => false,
    })
}

/// Whether two faces lie on one surface: equal stored surfaces, or planes
/// within the resolution of one plane facing one way (a wall or a cap split
/// by another primitive, its stored frames rounded apart); S9e.4b.4b.2b.2a:
/// or cylinders whose stored axes lie on one line within the resolution
/// with radii within it, their material on one side of their quadric (a
/// wall another part splits, the kernel's writer storing each piece on a
/// cylinder of its own, its origin moved along the axis).
fn one_surface(t: &Topology, a: usize, b: usize, tol: f64) -> bool {
    let (fa, fb) = (&t.faces()[a], &t.faces()[b]);
    if fa.surface == fb.surface {
        return true;
    }
    if let (
        Surface::Cylinder {
            frame: ca,
            radius: ra,
        },
        Surface::Cylinder {
            frame: cb,
            radius: rb,
        },
    ) = (&fa.surface, &fb.surface)
    {
        return coaxial(ca, cb, tol)
            && coaxial(cb, ca, tol)
            && (ra - rb).abs() <= tol
            && material_outside(t, a) == material_outside(t, b);
    }
    let (Surface::Plane(pa), Surface::Plane(pb)) = (&fa.surface, &fb.surface) else {
        return false;
    };
    let (Some(ma), Some(mb)) = (outward(t, a), outward(t, b)) else {
        return false;
    };
    parallel(ma, mb) && ma.dot(mb) > 0.0 && pa.coordinates(pb.origin())[2].abs() <= tol
}

/// S9e.4b.4b.2b.2b.3a: [`one_surface`] but two plane faces facing apart,
/// even on one stored plane (a stack's touching caps, which OCCT writes on
/// one surface): the two prisms' faces apart (`leaves_with`'s `strict`).
fn one_side(t: &Topology, a: usize, b: usize, tol: f64) -> bool {
    let apart = match (outward(t, a), outward(t, b)) {
        (Some(ma), Some(mb)) => ma.dot(mb) < 0.0,
        _ => false,
    };
    !apart && one_surface(t, a, b, tol)
}

/// S9e.4b.4b.2b.2a: whether a face rises from a plane face along `up` (its
/// inward normal): its range along it, on the plane's stored frame turned
/// to it, past the resolution (a wall of a part the plane is cut from lies
/// below it).
fn rises(t: &Topology, plane: usize, k: usize, up: Vec3, tol: f64) -> bool {
    let Surface::Plane(frame) = &t.faces()[plane].surface else {
        return false;
    };
    let frame = if frame.normal().dot(up) > 0.0 {
        *frame
    } else {
        frame.flipped()
    };
    face_range(t, k, &frame).1 > tol
}

/// S9e.4b.4b.2a: the prism leaves a body of several primitives may hold, in
/// the order tried: for each pair of plane faces facing apart (in face
/// order), each as the bottom in turn, the bottom's loops every edge of
/// which another face's wall along its normal holds (its outer loop among
/// them), the walls across them and the faces on their and the caps'
/// surfaces.
fn leaves(t: &Topology, tolerance: Tolerance) -> Vec<Leaf> {
    leaves_with(t, tolerance, false, false)
}

/// S9e.4b.4b.2b.2b.1: the leaves whose top cap another part or a pocket cuts
/// at its rim (`leaves_with`'s `rim`), in the order found, but those `leaves`
/// finds (`known`).
fn rim_leaves(t: &Topology, tolerance: Tolerance, known: &[Leaf]) -> Vec<Leaf> {
    leaves_with(t, tolerance, true, false)
        .into_iter()
        .filter(|l| {
            !known
                .iter()
                .any(|k| k.caps.0 == l.caps.0 && k.caps.1 == l.caps.1 && k.faces == l.faces)
        })
        .collect()
}

/// [`leaves`]; S9e.4b.4b.2b.2b.1: with `rim`, a wall across the bottom's
/// loops may end at or below the top cap's plane instead of reaching it (its
/// faces' range along the bottom's inward normal, on the bottom's stored
/// frame turned to it, within the resolution of the top's height or below:
/// the top cap cut at its rim by another part or a pocket there), a wall of
/// the outer loop reaching the top (a hole's walls may all end on another
/// part: a well). S9e.4b.4b.2b.2b.3a: with `strict`, two plane faces facing
/// apart are never on one surface (`one_side`: a stack's touching caps each
/// its own prism's).
fn leaves_with(t: &Topology, tolerance: Tolerance, rim: bool, strict: bool) -> Vec<Leaf> {
    let tol = tolerance.linear();
    let one_surface = |a: usize, b: usize, tol: f64| {
        if strict {
            one_side(t, a, b, tol)
        } else {
            one_surface(t, a, b, tol)
        }
    };
    let n = t.faces().len();
    let fins = fin_faces(t);
    let normals: Vec<Option<Vec3>> = (0..n)
        .map(|i| match t.faces()[i].surface {
            Surface::Plane(_) => outward(t, i),
            _ => None,
        })
        .collect();
    // The pairs of faces sharing an edge.
    let mut pairs = std::collections::BTreeSet::new();
    for e in t.edges() {
        if let [a, b] = e.fins.as_slice() {
            let (a, b) = (fins[a.index()], fins[b.index()]);
            pairs.insert((a.min(b), a.max(b)));
        }
    }
    let meet = |a: usize, b: usize| pairs.contains(&(a.min(b), a.max(b)));
    let mut out = Vec::new();
    for (i, mi) in normals.iter().enumerate() {
        let Some(mi) = *mi else { continue };
        for (j, mj) in normals.iter().enumerate().skip(i + 1) {
            let Some(mj) = *mj else { continue };
            if mi.dot(mj) >= 0.0 || !parallel(mi, mj) {
                continue;
            }
            for (bottom, top, m) in [(i, j, mi), (j, i, mj)] {
                let up = -m;
                // S9e.4b.4b.2b.2b.1: the bottom's stored frame turned to `up`
                // and the top's height on it.
                let Surface::Plane(stored) = &t.faces()[bottom].surface else {
                    continue;
                };
                let bottom_frame = if stored.normal().dot(up) > 0.0 {
                    *stored
                } else {
                    stored.flipped()
                };
                let height = face_range(t, top, &bottom_frame).1;
                let mut faces = vec![false; n];
                faces[bottom] = true;
                faces[top] = true;
                // The top cap's faces (one split by another primitive).
                let tops: Vec<usize> = (0..n)
                    .filter(|&k| k == top || one_surface(top, k, tol))
                    .collect();
                let (mut loops, mut reaching) = (Vec::new(), false);
                for (li, l) in t.faces()[bottom].loops.iter().enumerate() {
                    let Loop::Edges { fins: lf, .. } = &t.loops()[l.index()] else {
                        continue;
                    };
                    let across: Vec<usize> = lf
                        .iter()
                        .filter_map(|f| {
                            let e = &t.edges()[t.fins()[f.index()].edge.index()];
                            e.fins.iter().find(|x| *x != f).map(|x| fins[x.index()])
                        })
                        .collect();
                    // Walls reaching the top cap (else a prism on another
                    // primitive, S9e.4b.4b.2b's); S9e.4b.4b.2b.2a: through any
                    // face on their surface (a wall another part splits),
                    // each rising from the bottom (a piece below it, of a
                    // part the cap is cut from, is not this prism's).
                    let reach = |k: usize| {
                        (0..n)
                            .filter(|&m| m == k || one_surface(k, m, tol))
                            .any(|m| tops.iter().any(|&c| meet(m, c)))
                    };
                    // S9e.4b.4b.2b.2b.1: or ending at or below the top cap's
                    // plane (the top cap cut at its rim there).
                    let below = |k: usize| rim && face_range(t, k, &bottom_frame).1 <= height + tol;
                    if across.len() == lf.len()
                        && across.iter().all(|&k| {
                            k < n
                                && wall(t, k, up)
                                && (reach(k) || below(k))
                                && rises(t, bottom, k, up, tol)
                        })
                    {
                        reaching |= li == 0 && across.iter().any(|&k| reach(k));
                        loops.push(li);
                        for k in across {
                            faces[k] = true;
                        }
                    }
                }
                // S9e.4b.4b.2b.2b.1: a wall of the outer loop reaching the top
                // (a hole's walls may all end on another part: a well).
                if loops.first() != Some(&0) || (rim && !reaching) {
                    continue;
                }
                let members: Vec<usize> = (0..n).filter(|&k| faces[k]).collect();
                for (k, face) in faces.iter_mut().enumerate() {
                    if !*face && members.iter().any(|&m| one_surface(m, k, tol)) {
                        *face = true;
                    }
                }
                out.push(Leaf {
                    caps: (bottom, top, up),
                    loops,
                    faces,
                    open: None,
                    run: Vec::new(),
                    void: false,
                });
            }
        }
    }
    out
}

/// S9e.4b.4b.2b.2b.2: the parts whose walls run along a lead leaf's walls
/// (REVIEW_NOTES.md, "S9e.4b.4b.2b.2b.2 refined"), in the order found: two
/// plane faces of no lead face, parallel to the lead's caps and strictly
/// between them (on its bottom's stored frame turned up), the second lying
/// from the first along the first's inward normal (a part fused, standing
/// out of the lead's walls) or its outward normal (a void cut, a pocket open
/// sideways), the first's every loop edge across a lead wall (its run) or a
/// wall rising from it to the second (through any face on its surface), at
/// least one across a lead wall and one not; or a void open through the
/// lead's cap: a floor alone whose walls not the lead's rise from it to that
/// cap's plane (within the resolution), a prism of one cap past them
/// (`open`). Parts alike in their own walls (a part read off either cap) are
/// one.
fn flush_parts(t: &Topology, tolerance: Tolerance, lead: &Leaf) -> Vec<Leaf> {
    let tol = tolerance.linear();
    let n = t.faces().len();
    let fins = fin_faces(t);
    let (qb, qt, up) = lead.caps;
    let Surface::Plane(stored) = &t.faces()[qb].surface else {
        return Vec::new();
    };
    let base = if stored.normal().dot(up) > 0.0 {
        *stored
    } else {
        stored.flipped()
    };
    let height = face_range(t, qt, &base).1;
    let within = |k: usize| {
        let (a, b) = face_range(t, k, &base);
        a > tol && b < height - tol
    };
    // The faces across a face's loop's edges (none where a fin has no
    // other).
    let across = |face: usize, li: usize| -> Option<Vec<usize>> {
        let l = t.faces()[face].loops.get(li)?;
        let Loop::Edges { fins: lf, .. } = &t.loops()[l.index()] else {
            return None;
        };
        lf.iter()
            .map(|f| {
                let e = &t.edges()[t.fins()[f.index()].edge.index()];
                e.fins.iter().find(|x| *x != f).map(|x| fins[x.index()])
            })
            .collect()
    };
    let mut pairs = std::collections::BTreeSet::new();
    for e in t.edges() {
        if let [a, b] = e.fins.as_slice() {
            let (a, b) = (fins[a.index()], fins[b.index()]);
            pairs.insert((a.min(b), a.max(b)));
        }
    }
    let meet = |a: usize, b: usize| pairs.contains(&(a.min(b), a.max(b)));
    let mut out: Vec<Leaf> = Vec::new();
    let mut seen: Vec<Vec<usize>> = Vec::new();
    for void in [false, true] {
        for f in 0..n {
            if lead.faces[f] || !matches!(t.faces()[f].surface, Surface::Plane(_)) || !within(f) {
                continue;
            }
            let Some(m) = outward(t, f) else { continue };
            if !parallel(m, up) {
                continue;
            }
            // The part's direction from its first cap: into the material, or
            // (a void) out of it.
            let dir = if void { m } else { -m };
            let Surface::Plane(fs) = &t.faces()[f].surface else {
                continue;
            };
            let floor = if fs.normal().dot(dir) > 0.0 {
                *fs
            } else {
                fs.flipped()
            };
            // A void open across the lead's cap (a pocket open sideways and
            // through its top): every wall across the floor's loops a lead
            // wall or one rising from it to the lead's cap plane.
            if void {
                let far = if dir.dot(up) > 0.0 { qt } else { qb };
                let cap = face_range(t, far, &floor).0;
                let mut walls = Vec::new();
                let mut fit = true;
                for li in 0..t.faces()[f].loops.len() {
                    let Some(ws) = across(f, li) else {
                        fit = false;
                        break;
                    };
                    for &k in &ws {
                        fit &= wall(t, k, dir)
                            && (lead.faces[k]
                                || (rises(t, f, k, dir, tol)
                                    && face_range(t, k, &floor).1 >= cap - tol));
                    }
                    walls.extend(ws);
                }
                let run: Vec<bool> = (0..n)
                    .map(|k| lead.faces[k] && walls.contains(&k))
                    .collect();
                walls.sort_unstable();
                walls.dedup();
                let own: Vec<usize> = walls.iter().copied().filter(|&k| !run[k]).collect();
                if fit && run.iter().any(|&r| r) && !own.is_empty() && !seen.contains(&own) {
                    seen.push(own.clone());
                    out.push(Leaf {
                        caps: (f, f, dir),
                        loops: (0..t.faces()[f].loops.len()).collect(),
                        faces: part_faces(t, &[f], &own, lead, tol),
                        open: Some(own),
                        run,
                        void,
                    });
                    continue;
                }
            }
            for c in 0..n {
                if c == f || lead.faces[c] || !within(c) {
                    continue;
                }
                let Some(mc) = outward(t, c) else { continue };
                if mc.dot(m) >= 0.0 || !parallel(mc, m) || face_range(t, c, &floor).0 <= tol {
                    continue;
                }
                let ceiling: Vec<usize> = (0..n)
                    .filter(|&k| k == c || one_surface(t, c, k, tol))
                    .collect();
                let mut walls = Vec::new();
                let mut fit = true;
                for li in 0..t.faces()[f].loops.len() {
                    let Some(ws) = across(f, li) else {
                        fit = false;
                        break;
                    };
                    for &k in &ws {
                        fit &= wall(t, k, dir)
                            && (lead.faces[k]
                                || (rises(t, f, k, dir, tol)
                                    && (0..n)
                                        .filter(|&x| x == k || one_surface(t, k, x, tol))
                                        .any(|x| ceiling.iter().any(|&y| meet(x, y)))));
                    }
                    walls.extend(ws);
                }
                let run: Vec<bool> = (0..n)
                    .map(|k| lead.faces[k] && walls.contains(&k))
                    .collect();
                if !fit || !run.iter().any(|&r| r) || walls.iter().all(|&k| run[k]) && !void {
                    continue;
                }
                walls.sort_unstable();
                walls.dedup();
                // Read off either cap, a part's own walls are one set.
                let own: Vec<usize> = walls.iter().copied().filter(|&k| !run[k]).collect();
                if seen.contains(&own) {
                    continue;
                }
                seen.push(own);
                let own: Vec<usize> = walls.iter().copied().filter(|&k| !run[k]).collect();
                out.push(Leaf {
                    caps: (f, c, dir),
                    loops: (0..t.faces()[f].loops.len()).collect(),
                    faces: part_faces(t, &[f, c], &own, lead, tol),
                    open: None,
                    run,
                    void,
                });
            }
        }
    }
    out
}

/// S9e.4b.4b.2b.2b.2: a part's faces: its caps, its own walls and every face
/// of no lead face on one of their surfaces.
fn part_faces(t: &Topology, caps: &[usize], own: &[usize], lead: &Leaf, tol: f64) -> Vec<bool> {
    let n = t.faces().len();
    let members: Vec<usize> = caps.iter().chain(own).copied().collect();
    (0..n)
        .map(|k| {
            members.contains(&k)
                || (!lead.faces[k] && members.iter().any(|&x| one_surface(t, x, k, tol)))
        })
        .collect()
}

/// S9e.4b.4b.2b.2b.3a: the leaves a body's stack chains are tried with: the
/// leaves and then those cut at their rim, found as `leaves` and
/// `rim_leaves` find them but with two plane faces facing apart never on one
/// surface (`leaves_with`'s `strict`: a stack's touching caps, which OCCT
/// writes on one surface, each its own prism's).
fn stack_leaves(t: &Topology, tolerance: Tolerance) -> Vec<Leaf> {
    let mut out = leaves_with(t, tolerance, false, true);
    for l in leaves_with(t, tolerance, true, true) {
        if !out
            .iter()
            .any(|k| k.caps.0 == l.caps.0 && k.caps.1 == l.caps.1 && k.faces == l.faces)
        {
            out.push(l);
        }
    }
    out
}

/// S9e.4b.4b.2b.2b.2: the lead leaves a body's flush chains are tried
/// with, in order: each leaf but one alike an earlier in its faces and axis,
/// at most `FLUSH_LEADS`.
fn flush_leads(every: &[Leaf]) -> Vec<usize> {
    let mut out: Vec<usize> = Vec::new();
    for (i, l) in every.iter().enumerate() {
        if out.len() < FLUSH_LEADS
            && !out
                .iter()
                .any(|&j| every[j].faces == l.faces && parallel(every[j].caps.2, l.caps.2))
        {
            out.push(i);
        }
    }
    out
}

/// S9e.4b.4b.2b.2b.2: the most lead leaves a body's flush chains try.
const FLUSH_LEADS: usize = 8;

/// S9e.4b.4b.2b.2b.2: the most parts a flush chain completes.
const FLUSH_PARTS: usize = 4;

/// S9e.4b.4b.2b.2a: how many distinct leaves' face sets two leaves are
/// chosen from (the first, in the order found).
const LEAF_SETS: usize = 6;

/// S9e.4b.4b.2b.2a: the pairs of prism leaves a body may be (indices in
/// `leaves`): for each two of the leaves' first `LEAF_SETS` distinct face
/// sets with no face in common, in their order, the leaves of the two sets
/// whose axes are parallel within `ALIGN` where any are (the first such two),
/// else the first of each; each pair in both orders.
fn leaf_pairs(all: &[Leaf]) -> Vec<[usize; 2]> {
    let mut sets: Vec<Vec<usize>> = Vec::new();
    for (i, l) in all.iter().enumerate() {
        match sets.iter_mut().find(|g| all[g[0]].faces == l.faces) {
            Some(g) => g.push(i),
            None => sets.push(vec![i]),
        }
    }
    sets.truncate(LEAF_SETS);
    let mut out = Vec::new();
    for (x, gx) in sets.iter().enumerate() {
        for gy in &sets[x + 1..] {
            let apart = all[gx[0]]
                .faces
                .iter()
                .zip(&all[gy[0]].faces)
                .all(|(a, b)| !(*a && *b));
            if !apart {
                continue;
            }
            let [i, j] = gx
                .iter()
                .flat_map(|&i| gy.iter().map(move |&j| [i, j]))
                .find(|&[i, j]| parallel(all[i].caps.2, all[j].caps.2))
                .unwrap_or([gx[0], gy[0]]);
            out.push([i, j]);
            out.push([j, i]);
        }
    }
    out
}

/// S9e.4b.4b.1: a body of several primitives' curved faces as a Boolean
/// chain of its primitives (REVIEW_NOTES.md, "S9e.4b.4b refined"): each
/// surface's primitive (`primitive_over`), every plane face an end of one;
/// in the order of their reach from their axes or centres, widest first, the
/// first's material inside its quadric; each next one cut where its
/// material lies outside its quadric, else in common where its faces meet
/// the earlier ones' curved faces along convex edges, fused where along
/// concave edges or none. S9e.4b.4b.2a: or such a chain led by a prism leaf
/// (`leaves`): `choice` 0 the primitives alone, `k` the `k`-th leaf.
/// S9e.4b.4b.2b.1: from `L + 1` on (`L` the leaves) the same chains with
/// the other plane faces (`chain_piece`'s `extra`): `L + 1` the primitives'
/// alone, `L + 1 + k` the `k`-th leaf's. S9e.4b.4b.2b.2a: from `2 (L + 1)`
/// on (`P` the pairs of leaves, `leaf_pairs`) each pair's chain, then from
/// `2 (L + 1) + P` each pair's with the other plane faces. S9e.4b.4b.2b.2b.1:
/// from `N = 2 (L + 1) + 2 P` on, the chains again over every leaf and the
/// leaves cut at their rim (`rim_leaves`; `L'` in all, `P'` their pairs), a
/// tooth to a leaf's cap running past it (`chain_piece`'s `rim`): `N + k` the
/// `k`-th leaf's alone, `N + L' + k` with the other plane faces, `N + 2 L' +
/// k` the `k`-th pair's alone, `N + 2 L' + P' + k` with them.
/// S9e.4b.4b.2b.2b.2: from `M = N + 2 L' + 2 P'` on (`F` the leads,
/// `flush_leads`), `M + k` the `k`-th lead's flush chain (`flush_parts`)
/// alone, `M + F + k` with the other plane faces. S9e.4b.4b.2b.2b.3a: from
/// `S = M + 2 F` on (`L''` the stack leaves, `stack_leaves`, `P''` their pairs),
/// `S + k` the `k`-th pair's stack chain alone, `S + P'' + k` with the other
/// plane faces, `S + 2 P'' + k` the `k`-th stack leaf with them (`chain_piece`'s
/// `stack`).
fn primitives_piece(
    t: &Topology,
    tolerance: Tolerance,
    op: OperationId,
    choice: u32,
) -> Result<Piece> {
    primitives_piece_with(t, tolerance, op, choice, true)
}

/// S9e.4b.4c.2b.1: the most chains of one body reading a pocket with pockets
/// of its own whose models are refused before the later chains refuse such a
/// component as before (each such refusal arranges the chain's Booleans: a
/// body on a frame turned off the world's `z`, its pocket walls within
/// rounding of its cylinders' direction, tried every chain so).
const NESTED_CHAINS: usize = 4;

/// [`primitives_piece`], a component of the other plane faces of both kinds
/// read as a pocket with pockets of its own only where `nesting`.
fn primitives_piece_with(
    t: &Topology,
    tolerance: Tolerance,
    op: OperationId,
    choice: u32,
    nesting: bool,
) -> Result<Piece> {
    let all = leaves(t, tolerance);
    let count = all.len() + 1;
    let c = choice as usize;
    let first = 2 * count + 2 * leaf_pairs(&all).len();
    if c >= first {
        let mut every = all.clone();
        every.extend(rim_leaves(t, tolerance, &all));
        let (lone, pairs) = (every.len(), leaf_pairs(&every));
        let k = c - first;
        // S9e.4b.4b.2b.2b.2: then each lead's flush chain, alone and with
        // the other plane faces.
        if k >= 2 * lone + 2 * pairs.len() {
            let leads = flush_leads(&every);
            let j = k - 2 * lone - 2 * pairs.len();
            // S9e.4b.4b.2b.2b.3a: then each pair's stack chain (its second
            // leaf on the first's frame and levels), alone and with the
            // other plane faces.
            if j >= 2 * leads.len() {
                let s = j - 2 * leads.len();
                let strict = stack_leaves(t, tolerance);
                let pairs = leaf_pairs(&strict);
                let (chosen, extra) = if s < 2 * pairs.len() {
                    let [i, l] = *pairs
                        .get(s % pairs.len().max(1))
                        .ok_or_else(not_their_chain)?;
                    (vec![strict[i].clone(), strict[l].clone()], s >= pairs.len())
                } else {
                    // Each leaf with the other plane faces, a prism of one
                    // cap read with its holes.
                    let leaf = strict
                        .get((s - 2 * pairs.len()) % strict.len().max(1))
                        .cloned()
                        .ok_or_else(not_their_chain)?;
                    (vec![leaf], true)
                };
                let mut piece = chain_piece(
                    t,
                    tolerance,
                    op,
                    &chosen,
                    extra.then_some(&strict[..]),
                    true,
                    false,
                    true,
                    nesting,
                )?;
                piece.choice = choice;
                return Ok(piece);
            }
            let lead = every
                .get(
                    *leads
                        .get(j % leads.len().max(1))
                        .ok_or_else(not_their_chain)?,
                )
                .cloned()
                .ok_or_else(not_their_chain)?;
            let mut chosen = vec![lead.clone()];
            chosen.extend(flush_parts(t, tolerance, &lead));
            if chosen.len() < 2 || chosen.len() > FLUSH_PARTS + 1 {
                return Err(not_their_chain());
            }
            let extra = j >= leads.len();
            let mut piece = chain_piece(
                t,
                tolerance,
                op,
                &chosen,
                extra.then_some(&every[..]),
                true,
                true,
                false,
                nesting,
            )?;
            piece.choice = choice;
            return Ok(piece);
        }
        let (chosen, extra) = if k < 2 * lone {
            let leaf = every
                .get(k % lone.max(1))
                .cloned()
                .ok_or_else(not_their_chain)?;
            (vec![leaf], k >= lone)
        } else {
            let k = k - 2 * lone;
            let [i, j] = *pairs
                .get(k % pairs.len().max(1))
                .ok_or_else(not_their_chain)?;
            (vec![every[i].clone(), every[j].clone()], k >= pairs.len())
        };
        let mut piece = chain_piece(
            t,
            tolerance,
            op,
            &chosen,
            extra.then_some(&every[..]),
            true,
            false,
            false,
            nesting,
        )?;
        piece.choice = choice;
        return Ok(piece);
    }
    let (chosen, extra): (Vec<Leaf>, bool) = if c < 2 * count {
        let k = c % count;
        let leaf = match k {
            0 => Vec::new(),
            k => vec![all.get(k - 1).cloned().ok_or_else(not_their_chain)?],
        };
        (leaf, c >= count)
    } else {
        let pairs = leaf_pairs(&all);
        let k = c - 2 * count;
        let [i, j] = *pairs
            .get(k % pairs.len().max(1))
            .ok_or_else(not_their_chain)?;
        (vec![all[i].clone(), all[j].clone()], k >= pairs.len())
    };
    let mut piece = chain_piece(
        t,
        tolerance,
        op,
        &chosen,
        extra.then_some(&all[..]),
        false,
        false,
        false,
        nesting,
    )?;
    piece.choice = choice;
    Ok(piece)
}

/// S9e.4b.4b.2a: the chains a body of several primitives is tried as, in
/// order: its primitives alone, then each prism leaf; S9e.4b.4b.2b.1: then
/// each leaf's with the other plane faces, then the primitives' alone with
/// them; S9e.4b.4b.2b.2a: then each pair of leaves, then each pair's with
/// the other plane faces (every body S9e.4b.4b.2b.1 imports keeps its
/// chain); S9e.4b.4b.2b.2b.1: then each leaf's and each pair's again over
/// the leaves cut at their rim too, alone and with the other plane faces
/// (every body S9e.4b.4b.2b.2a imports keeps its chain). Returns them, the
/// first of those with the other plane faces or of two leaves (a tangency in
/// such a chain's is one the body is not) and the first over the leaves cut
/// at their rim (any refusal of such a chain's is one too, and of a flush
/// chain's, S9e.4b.4b.2b.2b.2, and of a stack chain's, S9e.4b.4b.2b.2b.3a,
/// tried last).
fn primitives_choices(t: &Topology, tolerance: Tolerance) -> (Vec<u32>, u32, u32) {
    let all = leaves(t, tolerance);
    let count = u32::try_from(all.len() + 1).unwrap_or(u32::MAX / 4);
    let pairs = u32::try_from(leaf_pairs(&all).len()).unwrap_or(0);
    let mut every = all.clone();
    every.extend(rim_leaves(t, tolerance, &all));
    let lone = u32::try_from(every.len()).unwrap_or(0);
    let rim_pairs = u32::try_from(leaf_pairs(&every).len()).unwrap_or(0);
    let leads = u32::try_from(flush_leads(&every).len()).unwrap_or(0);
    // S9e.4b.4b.2b.2b.3a: each pair of the stack leaves alone and with the
    // other plane faces, then each stack leaf with them.
    let strict = stack_leaves(t, tolerance);
    let stacks = u32::try_from(2 * leaf_pairs(&strict).len() + strict.len()).unwrap_or(0);
    let first = 2 * count + 2 * pairs;
    let flush = first + 2 * lone + 2 * rim_pairs;
    let order = (0..count)
        .chain((count + 1..2 * count).chain([count]))
        .chain(2 * count..first)
        .chain(first..flush)
        .chain(flush..flush + 2 * leads)
        .chain(flush + 2 * leads..flush + 2 * leads + stacks)
        .collect();
    (order, count, first)
}

#[cfg(test)]
thread_local! {
    /// The chains' models `imported_primitives` built on this thread.
    static CHAIN_MODELS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Whether two chains read off one body are the same but for their choice
/// (`imported_primitives` builds each model once).
fn same_chain(a: &Piece, b: &Piece) -> bool {
    let Piece {
        primitive,
        others,
        planes,
        form,
        tree,
        choice: _,
        chain,
        stack,
        nested,
    } = a;
    *primitive == b.primitive
        && *others == b.others
        && *planes == b.planes
        && *form == b.form
        && *tree == b.tree
        && *chain == b.chain
        && *stack == b.stack
        && *nested == b.nested
}

/// S9e.4b.4b.2b.1: the prisms of one cap a body holds besides its leaves
/// `lead`, in face order: a plane face no leaf's cap (of `all`) and no face
/// of `lead`'s or of an earlier such prism's, whose outer loop's every edge's
/// other face is a wall along its inward normal (a prism boss: its walls
/// reaching another part's faces, its other cap hidden), its faces it, those
/// walls and every face on their surfaces. S9e.4b.4b.2b.2a: its walls rise
/// from it along that normal (else a pocket's floor, its walls up to its
/// rim).
fn bosses(
    t: &Topology,
    tolerance: Tolerance,
    lead: &[Leaf],
    all: &[Leaf],
    holed: bool,
) -> Vec<Leaf> {
    let tol = tolerance.linear();
    let n = t.faces().len();
    let fins = fin_faces(t);
    let mut taken = vec![false; n];
    for l in lead {
        for (k, x) in l.faces.iter().enumerate() {
            taken[k] |= *x;
        }
    }
    let caps: Vec<usize> = all.iter().flat_map(|l| [l.caps.0, l.caps.1]).collect();
    let mut out = Vec::new();
    for b in 0..n {
        if taken[b] || caps.contains(&b) || !matches!(t.faces()[b].surface, Surface::Plane(_)) {
            continue;
        }
        let Some(m) = outward(t, b) else { continue };
        let up = -m;
        let Some(l) = t.faces()[b].loops.first() else {
            continue;
        };
        let Loop::Edges { fins: lf, .. } = &t.loops()[l.index()] else {
            continue;
        };
        let across: Vec<usize> = lf
            .iter()
            .filter_map(|f| {
                let e = &t.edges()[t.fins()[f.index()].edge.index()];
                e.fins.iter().find(|x| *x != f).map(|x| fins[x.index()])
            })
            .collect();
        let walls = across.len() == lf.len()
            && across.iter().all(|&k| {
                k < n && k != b && !taken[k] && wall(t, k, up) && rises(t, b, k, up, tol)
            });
        if !walls {
            continue;
        }
        // S9e.4b.4b.2b.2b.3a: where `holed`, its holes whose every wall rises
        // from it too (a lower prism whose other cap the upper one hides,
        // its holes closed by it).
        let mut loops = vec![0];
        let mut across = across;
        if holed {
            for (li, l) in t.faces()[b].loops.iter().enumerate().skip(1) {
                let Loop::Edges { fins: lf, .. } = &t.loops()[l.index()] else {
                    continue;
                };
                let hole: Vec<usize> = lf
                    .iter()
                    .filter_map(|f| {
                        let e = &t.edges()[t.fins()[f.index()].edge.index()];
                        e.fins.iter().find(|x| *x != f).map(|x| fins[x.index()])
                    })
                    .collect();
                if hole.len() == lf.len()
                    && hole.iter().all(|&k| {
                        k < n && k != b && !taken[k] && wall(t, k, up) && rises(t, b, k, up, tol)
                    })
                {
                    loops.push(li);
                    across.extend(hole);
                }
            }
        }
        let mut faces = vec![false; n];
        faces[b] = true;
        for &k in &across {
            faces[k] = true;
        }
        let members: Vec<usize> = (0..n).filter(|&k| faces[k]).collect();
        for k in 0..n {
            let alike = if holed {
                members.iter().any(|&m| one_side(t, m, k, tol))
            } else {
                members.iter().any(|&m| one_surface(t, m, k, tol))
            };
            if !faces[k] && !taken[k] && alike {
                faces[k] = true;
            }
        }
        for k in 0..n {
            taken[k] |= faces[k];
        }
        let mut walls = across;
        walls.sort_unstable();
        walls.dedup();
        out.push(Leaf {
            caps: (b, b, up),
            loops,
            faces,
            open: Some(walls),
            run: Vec::new(),
            void: false,
        });
    }
    out
}

/// S9e.4b.4b.2b.1: a plane face's stored frame, taken onto an axis of the
/// chain's (`axes`: a prism's or a cylinder's or a cone's frame) whose
/// direction its normal is parallel to within `ALIGN` (that axis's axes bit
/// for bit at its stored origin), and whether its normal `x * y` points into
/// the material (turned over where `over`).
fn plane_along(t: &Topology, i: usize, over: bool, axes: &[Frame3]) -> Result<(Frame3, bool)> {
    let Surface::Plane(stored) = &t.faces()[i].surface else {
        return Err(general());
    };
    let frame = axes
        .iter()
        .find(|a| parallel(a.normal(), stored.normal()))
        .map_or(*stored, |a| a.at(stored.origin()));
    let out = outward(t, i).ok_or_else(general)?;
    let into = frame.normal().dot(out) < 0.0;
    Ok((frame, into != over))
}

/// [`primitives_piece`] of a chain of primitives, led by a prism leaf where
/// given (S9e.4b.4b.2b.2a: or two, the second next); S9e.4b.4b.2b.1: with
/// the other plane faces where `extra` (the body's leaves): prisms of one cap
/// (`bosses`), each fused after the leaf; each component of the other plane
/// faces (through their edges with each other) convex within itself a flat
/// of the one part it meets along convex edges (that part common the hull of
/// its planes), or concave within itself a pocket (the hull of its planes
/// turned over, less the primitives it meets along concave edges: its teeth)
/// cut from the chain. S9e.4b.4b.2b.2b.1: where `rim`, no prism's face is a
/// primitive's cap (`primitive_over`'s `skip`: a tooth to a leaf's cap).
/// S9e.4b.4b.2b.2b.3a: where `stack`, a later leaf on the first's frame with
/// its caps taken onto the first's heights (`Levels`), the prisms of one cap
/// read with their holes (`bosses`' `holed`), and the body classified with
/// probes off both leaves' boundaries (`Piece::stack`).
/// S9e.4b.4c.2b.1: a component of the other plane faces of both kinds whose
/// edges with the parts are convex alone (its rims) or concave alone with
/// primitives (its teeth) a pocket with pockets of its own (`nested_pocket`);
/// any refusal of a chain reading one so is the one its chain gave before,
/// that component's (`not_their_faces`).
#[allow(clippy::too_many_arguments)]
fn chain_piece(
    t: &Topology,
    tolerance: Tolerance,
    op: OperationId,
    leaf: &[Leaf],
    extra: Option<&[Leaf]>,
    rim: bool,
    flush: bool,
    stack: bool,
    nesting: bool,
) -> Result<Piece> {
    let mut nested = false;
    chain_piece_with(
        t,
        tolerance,
        op,
        leaf,
        extra,
        rim,
        flush,
        stack,
        nesting,
        &mut nested,
    )
    .map_err(|e| if nested { not_their_faces() } else { e })
}

/// [`chain_piece`], `nested` set where it reads a pocket with pockets of its
/// own (only where `nesting`).
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn chain_piece_with(
    t: &Topology,
    tolerance: Tolerance,
    op: OperationId,
    leaf: &[Leaf],
    extra: Option<&[Leaf]>,
    rim: bool,
    flush: bool,
    stack: bool,
    nesting: bool,
    nested: &mut bool,
) -> Result<Piece> {
    let tol = tolerance.linear();
    let n = t.faces().len();
    let bossed = extra.map_or_else(Vec::new, |all| bosses(t, tolerance, leaf, all, stack));
    let prisms: Vec<&Leaf> = leaf.iter().chain(&bossed).collect();
    let prism_of = |k: usize| prisms.iter().position(|l| l.faces[k]);
    let skip: Vec<usize> = if rim {
        (0..n).filter(|&k| prism_of(k).is_some()).collect()
    } else {
        Vec::new()
    };
    let groups: Vec<Vec<usize>> = surfaces(t)
        .into_iter()
        .filter(|g| !g.iter().any(|&k| prism_of(k).is_some()))
        .collect();
    if prisms.is_empty() && groups.len() < 2 {
        return Err(general());
    }
    // Each face's part: a prism's (in order), else a primitive's (its faces
    // and its caps).
    let p = prisms.len();
    let mut part: Vec<Option<usize>> = (0..n).map(prism_of).collect();
    let mut order = Vec::new();
    for (g, faces) in groups.iter().enumerate() {
        let outside = material_outside(t, faces[0]);
        if faces.iter().any(|&k| material_outside(t, k) != outside) {
            return Err(not_their_chain());
        }
        let (_, capped, reach) =
            primitive_over(t, tolerance, op, faces, outside, None, &[], &skip)?;
        for &k in faces.iter().chain(&capped) {
            part[k] = Some(p + g);
        }
        order.push((g, outside, reach));
    }
    // The body's own faces tangent along an edge (a fillet on a primitive's
    // faces, a primitive on another's surface) but a prism's walls' joints
    // (S9e.4b.4b.2a: two walls along one axis tangent along a line along it,
    // a leaf's both); then plane faces other than the primitives' ends and a
    // leaf's are S9e.4b.4b.2b.2b's but the flats and pockets taken (below).
    let bends = bends(t)?;
    for (i, &(a, b, bend)) in bends.iter().enumerate() {
        // S9e.4b.4b.2b.2b.2: in a flush chain, a joint of any one prism
        // (parts of one prism read off several caps share its walls).
        let one = prisms.is_empty()
            || (prism_of(a).is_some() && prism_of(a) == prism_of(b))
            || (flush && prisms.iter().any(|l| l.faces[a] && l.faces[b]));
        // S9e.4b.4b.2b.1: two faces of one surface (a face OCCT split at its
        // surface's seam) are none.
        let split = t.faces()[a].surface == t.faces()[b].surface;
        if bend == Bend::Flat && !split && !(joint(t, i, a, b) && one) {
            return Err(tangent_body());
        }
    }
    let others: Vec<usize> = (0..n)
        .filter(|&i| matches!(t.faces()[i].surface, Surface::Plane(_)) && part[i].is_none())
        .collect();
    // One part alone is a construction of its own, but (S9e.4b.4b.2b.2a) a
    // lone leaf with the other plane faces (its pockets and flats).
    if prisms.len() + groups.len() < 2 && (others.is_empty() || extra.is_none()) {
        return Err(general());
    }
    if !others.is_empty() && extra.is_none() {
        return Err(not_their_faces());
    }
    // S9e.4b.4b.2b.1: with nothing more than the chain alone, its own.
    if extra.is_some() && others.is_empty() && bossed.is_empty() {
        return Err(general());
    }
    // The other plane faces' components, each a flat or a pocket.
    let mut root: Vec<usize> = (0..n).collect();
    for &(a, b, _) in &bends {
        if others.contains(&a) && others.contains(&b) {
            let (x, y) = (find(&mut root, a), find(&mut root, b));
            root[x.max(y)] = x.min(y);
        }
    }
    let mut flats: Vec<Vec<usize>> = vec![Vec::new(); p + groups.len()];
    // Each pocket's faces, its teeth and (S9e.4b.4c.2b.1) whether it holds
    // pockets of its own.
    let mut pockets: Vec<(Vec<usize>, Vec<usize>, bool)> = Vec::new();
    let mut teeth = vec![false; groups.len()];
    let mut seen = Vec::new();
    for &i in &others {
        let r = find(&mut root, i);
        if seen.contains(&r) {
            continue;
        }
        seen.push(r);
        let members: Vec<usize> = others
            .iter()
            .copied()
            .filter(|&k| find(&mut root, k) == r)
            .collect();
        let (mut convex, mut concave) = (false, false);
        let mut met: BTreeMap<usize, (bool, bool)> = BTreeMap::new();
        for &(a, b, bend) in &bends {
            let (x, y) = (members.contains(&a), members.contains(&b));
            if x && y {
                convex |= bend == Bend::Convex;
                concave |= bend == Bend::Concave;
            } else if x || y {
                let other = if x { b } else { a };
                let Some(q) = part[other] else {
                    return Err(not_their_faces());
                };
                let e = met.entry(q).or_default();
                e.0 |= bend == Bend::Convex;
                e.1 |= bend == Bend::Concave;
            }
        }
        if !concave {
            // A flat: the one part it meets along convex edges (the others
            // fused to it along concave ones).
            let trimmed: Vec<usize> = met.iter().filter(|(_, m)| m.0).map(|(q, _)| *q).collect();
            let [x] = trimmed.as_slice() else {
                return Err(not_their_faces());
            };
            flats[*x].extend(members);
        } else if !convex {
            // A pocket: its teeth the primitives it meets along concave
            // edges alone.
            let mut mine = Vec::new();
            for (&q, &(cx, cc)) in &met {
                match (cx, cc) {
                    (false, true) if q >= p => {
                        teeth[q - p] = true;
                        mine.push(q - p);
                    }
                    (true, false) => {}
                    _ => return Err(not_their_faces()),
                }
            }
            pockets.push((members, mine, false));
        } else if nesting
            && met
                .iter()
                .all(|(&q, &(cx, cc))| (cx && !cc) || (!cx && cc && q >= p))
        {
            // S9e.4b.4c.2b.1: of both kinds, its rims and teeth as a pocket's:
            // a pocket with pockets of its own (a post or an island of planes
            // standing in it).
            let mut mine = Vec::new();
            for (&q, &(cx, _)) in &met {
                if !cx {
                    teeth[q - p] = true;
                    mine.push(q - p);
                }
            }
            *nested = true;
            pockets.push((members, mine, true));
        } else {
            return Err(not_their_faces());
        }
    }
    if teeth
        .iter()
        .zip(&flats[p..])
        .any(|(t, f)| *t && !f.is_empty())
    {
        return Err(not_their_faces());
    }
    // Widest first (a stable order: ties in face order), the teeth apart.
    order.sort_by(|a, b| b.2.total_cmp(&a.2));
    let (tooth, order): (Vec<_>, Vec<_>) = order.into_iter().partition(|o| teeth[o.0]);
    if prisms.is_empty() && order.first().is_none_or(|o| o.1) {
        return Err(not_their_chain());
    }
    // The prisms first, on their bottom caps' frames (a boss parallel to the
    // leaf on its axes, bit for bit).
    let mut built: Vec<Solid> = Vec::new();
    for (i, l) in prisms.iter().enumerate() {
        let along = built.first().map(|s| s.frame);
        let open = (l.open.as_deref(), along.as_ref());
        // S9e.4b.4b.2b.2b.2: a part flush with the lead completed through
        // it (or around it, a void).
        let lead = built.first().cloned();
        let complete = lead
            .as_ref()
            .filter(|_| !l.run.is_empty())
            .map(|lead| Completion {
                run: &l.run,
                lead,
                into: !l.void,
            });
        // S9e.4b.4b.2b.2b.3a: a later leaf of a stack chain on the first's
        // frame, its caps taken onto the first's levels.
        let levels = stack.then(|| Levels {
            first: built.first().map(|s| (s.frame, [s.start, s.end])),
        });
        built.push(
            completed_on(
                t,
                l.caps,
                &l.loops,
                open,
                complete.as_ref(),
                levels.as_ref(),
                tolerance,
                mine(op, i),
            )
            .map_err(|e| match e {
                Error::OutOfDomain(_) => e,
                _ => not_their_chain(),
            })?,
        );
    }
    let lead = built.len();
    let lead_frame = built.first().map(|s| s.frame);
    // Each part's faces: a prism's, or a primitive's curved faces (not its
    // caps, as S9e.4b.4b.1's chain), with its flats.
    let faces_of = |q: usize| -> Vec<usize> {
        let own = |k: usize| {
            if q < p {
                part[k] == Some(q)
            } else {
                groups[q - p].contains(&k)
            }
        };
        (0..n)
            .filter(|&k| own(k) || flats[q].contains(&k))
            .collect()
    };
    // Each one's Boolean with the earlier ones and its frame.
    let mut ops: Vec<Op2> = Vec::new();
    let mut axes: Vec<Frame3> = built.iter().map(|s| s.frame).collect();
    let mut frames: Vec<Option<Frame3>> = Vec::new();
    let mut earlier: Vec<usize> = Vec::new();
    let meets = |mine: &[usize], earlier: &[usize]| {
        let (mut convex, mut concave) = (false, false);
        for &(a, b, bend) in &bends {
            if (mine.contains(&a) && earlier.contains(&b))
                || (mine.contains(&b) && earlier.contains(&a))
            {
                convex |= bend == Bend::Convex;
                concave |= bend == Bend::Concave;
            }
        }
        (convex, concave)
    };
    for (q, l) in prisms.iter().enumerate().take(lead) {
        let mine = faces_of(q);
        let (convex, concave) = meets(&mine, &earlier);
        // S9e.4b.4b.2b.2b.2: a flush part fused, a void cut.
        ops.push(if q == 0 {
            Op2::Common
        } else if l.void {
            Op2::Cut
        } else if !l.run.is_empty() {
            Op2::Fuse
        } else if convex && concave {
            return Err(not_their_chain());
        } else if convex {
            Op2::Common
        } else {
            Op2::Fuse
        });
        earlier.extend(mine);
    }
    let chained: Vec<(usize, bool, f64)> = order.iter().chain(&tooth).copied().collect();
    for (i, &(g, outside, _)) in chained.iter().enumerate() {
        let faces = &groups[g];
        let own = axis_of(t, faces[0]).ok_or_else(general)?;
        // A cylinder or a cone on an earlier primitive's axis takes its
        // frame; S9e.4b.4b.2a: one along the leaf's axis its axes, bit for
        // bit, at its stored origin.
        let sphere = matches!(t.faces()[faces[0]].surface, Surface::Sphere { .. });
        let along = lead_frame
            .filter(|f| parallel(f.normal(), own.normal()))
            .map(|f| f.at(own.origin()));
        let reference = (!sphere)
            .then(|| {
                axes[lead..]
                    .iter()
                    .find(|a| coaxial(a, &own, tol))
                    .copied()
                    .or(along)
            })
            .flatten();
        if i < order.len() {
            let mine = faces_of(p + g);
            let (convex, concave) = meets(&mine, &earlier);
            ops.push(if i + lead == 0 {
                Op2::Common
            } else if outside {
                Op2::Cut
            } else if convex && concave {
                return Err(not_their_chain());
            } else if convex {
                Op2::Common
            } else {
                Op2::Fuse
            });
            earlier.extend(mine);
        }
        axes.push(reference.unwrap_or(own));
        frames.push(reference);
    }
    // Each primitive's entities under an operation of its own (the first
    // the body's construction's), past each open end by a quarter of its
    // faces' range; the first, a common's and a cut's past the others' too
    // (what bounds the body there, or what the cut leaves). A tooth as a
    // fused one.
    let mut solids = Vec::new();
    for (i, &(g, outside, _)) in chained.iter().enumerate() {
        let (solid, _, _) = primitive_over(
            t,
            tolerance,
            mine(op, i + lead),
            &groups[g],
            outside,
            frames[i].as_ref(),
            &[],
            &skip,
        )?;
        solids.push(solid);
    }
    let mut bounds: Vec<crate::Bounds3> = solids.iter().map(|s| s.bounds).collect();
    bounds.extend(built.iter().map(|s| s.bounds));
    for (i, &(g, outside, _)) in chained.iter().enumerate() {
        if i >= order.len() || ops[i + lead] == Op2::Fuse {
            continue;
        }
        let others: Vec<crate::Bounds3> = (0..bounds.len())
            .filter(|&j| j != i)
            .map(|j| bounds[j])
            .collect();
        solids[i] = primitive_over(
            t,
            tolerance,
            mine(op, i + lead),
            &groups[g],
            outside,
            frames[i].as_ref(),
            &others,
            &skip,
        )?
        .0;
    }
    // The tree: each part (common its flats' hull), the chain of them, each
    // pocket (less its teeth) cut from it.
    let hull_axes = axes.clone();
    let hull = |members: &[usize], over: bool| -> Result<Vec<(Frame3, bool)>> {
        let planes = members
            .iter()
            .map(|&k| plane_along(t, k, over, &hull_axes))
            .collect::<Result<Vec<_>>>()?;
        Ok(one_plane(planes, tol))
    };
    let index_of = |q: usize| -> usize {
        if q < lead {
            q
        } else {
            lead + chained.iter().position(|o| o.0 == q - p).expect("a part")
        }
    };
    let leaf_of = |q: usize| -> Result<Tree> {
        let x = Tree::Primitive(index_of(q));
        Ok(if flats[q].is_empty() {
            x
        } else {
            Tree::Op(
                Op2::Common,
                Box::new(x),
                Box::new(Tree::Hull(hull(&flats[q], false)?)),
            )
        })
    };
    let parts: Vec<usize> = (0..lead).chain(order.iter().map(|o| p + o.0)).collect();
    let mut tree = leaf_of(parts[0])?;
    for (k, &q) in parts.iter().enumerate().skip(1) {
        tree = Tree::Op(ops[k], Box::new(tree), Box::new(leaf_of(q)?));
    }
    // S9e.4b.4c.2b.1: the nested pockets' Booleans of the chain, together.
    let mut inner = 0;
    for (members, mine, holding) in &pockets {
        let mut pocket = if *holding {
            let (tree, outer) =
                nested_pocket(t, &bends, members, true, &hull_axes, tol, 1, &mut inner)?;
            // A tooth meets the pocket's own hull alone (a primitive standing
            // in a nested pocket is S9e.4b.4c.2b.2's).
            let deeper = |k: usize| members.contains(&k) && !outer.contains(&k);
            for &g in mine {
                if bends.iter().any(|&(a, b, _)| {
                    (deeper(a) && groups[g].contains(&b)) || (deeper(b) && groups[g].contains(&a))
                }) {
                    return Err(not_their_faces());
                }
            }
            tree
        } else {
            Tree::Hull(hull(members, true)?)
        };
        for &g in mine {
            pocket = Tree::Op(Op2::Cut, Box::new(pocket), Box::new(leaf_of(p + g)?));
        }
        tree = Tree::Op(Op2::Cut, Box::new(tree), Box::new(pocket));
    }
    let mut solids = built.into_iter().chain(solids);
    let primitive = solids.next().expect("two primitives");
    Ok(Piece {
        primitive: Box::new(primitive),
        others: solids.collect(),
        planes: Vec::new(),
        form: Form::Tree,
        tree: Some(tree),
        choice: 0,
        chain: true,
        stack,
        nested: *nested,
    })
}

/// S9e.4b.4c.2b.1: a chain's pocket holding pockets of its own, its region
/// at a depth of pockets (S9e.4b.4c.2a's `pocket` on the chain's axes): its
/// hull the planes of its faces (`plane_along`, turned over where `over`)
/// that leave every boundary point of its faces on its material's side
/// within the resolution, less each component of its other faces (joined by
/// their edges, in the order of their first faces) as a pocket of its own,
/// turned over again; at most `DEPTH` deep, its Booleans counted in `inner`
/// (at most `BOOLEANS` a chain). Returns the tree and its hull's faces.
#[allow(clippy::too_many_arguments)]
fn nested_pocket(
    t: &Topology,
    bends: &[(usize, usize, Bend)],
    faces: &[usize],
    over: bool,
    axes: &[Frame3],
    tol: f64,
    depth: usize,
    inner: &mut usize,
) -> Result<(Tree, Vec<usize>)> {
    let n = t.faces().len();
    let corners: Vec<Point3> = faces.iter().flat_map(|&i| boundary_points(t, i)).collect();
    let (mut hull, mut own, mut rest) = (Vec::new(), Vec::new(), Vec::new());
    for &i in faces {
        let (frame, into) = plane_along(t, i, over, axes)?;
        let apart = corners.iter().any(|&p| {
            let d = frame.coordinates(p)[2];
            if into {
                d < -tol
            } else {
                d > tol
            }
        });
        if apart {
            rest.push(i);
        } else {
            hull.push((frame, into));
            own.push(i);
        }
    }
    // A pocket none of whose faces bounds it, or holding pockets past the
    // depth.
    if hull.is_empty() || (!rest.is_empty() && depth >= DEPTH) {
        return Err(not_their_faces());
    }
    let mut root: Vec<usize> = (0..n).collect();
    for &(a, b, _) in bends {
        if rest.contains(&a) && rest.contains(&b) {
            let (x, y) = (find(&mut root, a), find(&mut root, b));
            root[x.max(y)] = x.min(y);
        }
    }
    let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &i in &rest {
        let r = find(&mut root, i);
        components.entry(r).or_default().push(i);
    }
    let mut tree = Tree::Hull(one_plane(hull, tol));
    for faces in components.into_values() {
        *inner += 1;
        if *inner > BOOLEANS {
            return Err(not_their_faces());
        }
        let (x, _) = nested_pocket(t, bends, &faces, !over, axes, tol, depth + 1, inner)?;
        tree = Tree::Op(Op2::Cut, Box::new(tree), Box::new(x));
    }
    Ok((tree, own))
}

/// A chain's `i`-th part's operation: the body's construction's for the
/// first.
fn mine(op: OperationId, i: usize) -> OperationId {
    OperationId(op.0 ^ ((i as u64) << 24))
}

/// A plane piece's bounds: its primitive's between the body's axial ends
/// (a cylinder's or a cone's), else its primitive's.
fn piece_bounds(r: &Solid, t: &Topology) -> crate::Bounds3 {
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

    /// The boolean target's `timeout-d737f94d`'s `STACKS` body (`k = 11 /
    /// 32`, `h = 1`: a square `[-2k, 2k]^2` with a hole of radius `k` under
    /// a stadium of half discs of radius `3k` about `(0, 0)` and `(3k, 0)`
    /// with a hole of radius `k / 2` about `(k / 4, k / 4)`, on the tilted
    /// frame turned about its axis), written, read back and imported: of the
    /// chains tried before its stack's, two pairs read the same chain off it
    /// (the stadium's prism through either face on the touching caps'
    /// plane), and each chain's model is built once (5 of 7 here), the
    /// stack's matching. Built again, they took two of the import's seven
    /// models.
    #[test]
    fn a_stacks_chains_are_each_built_once() {
        use super::{
            primitives_choices, primitives_piece, renamed, same_chain, Recognized, CHAIN_MODELS,
        };
        use crate::solid::Construction;
        let tol = Tolerance::default();
        let tilted = Frame3::new(
            Point3::new(1.0, -2.0, 0.5),
            Vec3::new(0.0, 3.0, 4.0),
            Vec3::new(1.0, 0.0, 0.0),
            tol,
        )
        .unwrap();
        let on = Frame3::new(
            tilted.origin(),
            tilted.normal(),
            tilted.x() * 8.0 + tilted.y() * 15.0,
            tol,
        )
        .unwrap();
        let (k, h) = (2.75 / 8.0, 1.0);
        let prism = |outer: Boundary, hole: Boundary, h0: f64, h1: f64, op: u64| {
            let profile = Profile::new(outer, vec![hole], tol).unwrap();
            Solid::extrude_with(OperationId(op), profile, on, h0, h1)
                .unwrap()
                .0
        };
        let square = [(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0)]
            .map(|(x, y)| Point2::new(x * k, y * k))
            .to_vec();
        let lower = prism(
            Boundary::polygon(square, tol).unwrap(),
            Boundary::circle(Point2::new(0.0, 0.0), k, tol).unwrap(),
            0.0,
            h,
            41,
        );
        let r = 3.0 * k;
        let arc = |x: f64| Segment::Arc {
            center: Point2::new(x, 0.0),
            radius: r,
            ccw: true,
        };
        let stadium = [(0.0, -r), (3.0 * k, -r), (3.0 * k, r), (0.0, r)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec();
        let segments = vec![Segment::Line, arc(3.0 * k), Segment::Line, arc(0.0)];
        let upper = prism(
            Boundary::path(stadium, segments, tol).unwrap(),
            Boundary::circle(Point2::new(0.25 * k, 0.25 * k), 0.5 * k, tol).unwrap(),
            h,
            2.0 * h,
            42,
        );
        let (out, _) = lower.fuse(OperationId(43), &upper).unwrap();
        let [body] = <[Solid; 1]>::try_from(out).unwrap();
        let text = crate::occt_brep::write(body.topology(), body.resolution().linear()).unwrap();
        let doc = crate::occt_brep::read(&text).unwrap();
        let [solid] = <[_; 1]>::try_from(crate::occt_brep::import(&doc).solids).unwrap();
        let (stored, resolution) = (solid.result.unwrap(), solid.tolerance);
        CHAIN_MODELS.with(|c| c.set(0));
        let (imported, _) =
            Solid::imported_with(OperationId(11), stored.clone(), resolution).unwrap();
        let built = CHAIN_MODELS.with(std::cell::Cell::get);
        let Construction::Imported(i) = &imported.construction else {
            panic!("an import")
        };
        let Recognized::Piece(p) = &i.recognized else {
            panic!("a chain")
        };
        assert!(p.stack, "the stack's chain");
        let v = body.mass_properties().volume;
        assert!((imported.mass_properties().volume - v).abs() <= 1e-9 * v);
        // The chains read off the body up to the stack's, as the import
        // reads them.
        let topology = renamed(stored, OperationId(11)).unwrap();
        let op = super::construction_operation(OperationId(11));
        let (choices, _, _) = primitives_choices(&topology, resolution);
        let mut read: Vec<super::Piece> = Vec::new();
        for c in choices {
            if let Ok(q) = primitives_piece(&topology, resolution, op, c) {
                read.push(q);
            }
            if c == p.choice {
                break;
            }
        }
        let distinct = read
            .iter()
            .enumerate()
            .filter(|(j, q)| !read[..*j].iter().any(|x| same_chain(x, q)))
            .count();
        assert!(distinct < read.len(), "{distinct} of {}", read.len());
        assert_eq!(built, distinct);
    }
}
