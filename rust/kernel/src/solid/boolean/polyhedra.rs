//! S9b: Booleans of polyhedral prisms in any relative position (REVIEW_NOTES
//! S9b's decisions, `BOOLEAN.md`).
//!
//! Each input is decided on its construction's exact model: a point `o + u
//! x + v y + w n` in rationals from its frame's stored binary64 origin and
//! axes and its profile's and heights' values, so every model vertex lies
//! exactly on its faces' planes. Each boundary face, cut into convex pieces
//! (a cap's trapezoids, a wall's rectangle), is split by every plane of the
//! other input's faces into fragments no such face crosses; each fragment is
//! classified on both of its sides by an infinitesimal push against the
//! other input's convex cells (closed half-spaces), and kept, oriented to
//! leave the result's material, where the operation's set function differs
//! across it (a fragment on both inputs' boundaries once, as the object's).
//! The kept fragments are made conforming (each edge split at every kept
//! vertex on it), joined across shared edges into maximal faces of one
//! plane and orientation, their edges joined where they run straight on
//! between the same two faces; each connected set of faces is a shell, an
//! outer one (positive exact volume) or a cavity (negative) of the outer one
//! holding it. Vertices and planes are rounded to binary64 once, each solid
//! validated as it is built.
//!
//! Names follow provenance as S9a's: a face continues the input faces its
//! fragments come from facing their way (a cut's tool's, or one facing the
//! other way, it touches); an edge along an input edge continues it, one
//! elsewhere lies on the faces meeting there; a vertex at an input vertex
//! continues it, one elsewhere lies on what it is on.
use super::{At, SlotPlan, What};
use crate::identity::{EntityId, EntityKind, OperationId, Role};
use crate::profile::boolean::{Op2, Operand};
use crate::solid::split::{q, rational_f64, zero};
use crate::solid::{Construction, MassProperties, Solid};
use crate::topology::{
    plane_pcurve, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId, Orientation,
    Region, RegionId, RegionKind, Shell, ShellId, Side, Slot, Surface, Topology, TopologyParts,
    Vertex, VertexId,
};
use crate::{Error, Frame3, Location, Point2, Point3, Profile, Result, Segment, Tolerance, Vec3};
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_rational::BigRational as R;
use std::collections::{BTreeMap, BTreeSet};

type V = [R; 3];

fn sub(a: &V, b: &V) -> V {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}

fn add(a: &V, b: &V) -> V {
    [&a[0] + &b[0], &a[1] + &b[1], &a[2] + &b[2]]
}

fn scale(a: &V, s: &R) -> V {
    [&a[0] * s, &a[1] * s, &a[2] * s]
}

fn dot(a: &V, b: &V) -> R {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

fn cross(a: &V, b: &V) -> V {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

fn neg(a: &V) -> V {
    [-&a[0], -&a[1], -&a[2]]
}

fn is_zero(a: &V) -> bool {
    a.iter().all(|x| *x == zero())
}

fn vq(v: Vec3) -> V {
    [q(v.x), q(v.y), q(v.z)]
}

fn abs(x: &R) -> R {
    if *x < zero() {
        -x
    } else {
        x.clone()
    }
}

fn rounded(p: &V) -> Point3 {
    Point3::new(
        rational_f64(&p[0]),
        rational_f64(&p[1]),
        rational_f64(&p[2]),
    )
}

/// A plane `n . X + d = 0`, outside where positive.
#[derive(Debug, Clone, PartialEq)]
struct Plane {
    n: V,
    d: R,
}

impl Plane {
    fn through(n: V, p: &V) -> Self {
        let d = -dot(&n, p);
        Self { n, d }
    }
    fn eval(&self, p: &V) -> R {
        dot(&self.n, p) + &self.d
    }
    /// The plane scaled so its first nonzero normal coordinate is `+-1`
    /// (`oriented`: keeping the side; else also its sign made positive).
    fn canonical(&self, oriented: bool) -> (V, R) {
        let first = self
            .n
            .iter()
            .find(|x| **x != zero())
            .expect("a plane's normal")
            .clone();
        let k = if oriented {
            if first > zero() {
                first
            } else {
                -first
            }
        } else {
            first
        };
        let inv = R::from_integer(1.into()) / k;
        (scale(&self.n, &inv), &self.d * &inv)
    }
}

/// A polyhedral Boolean's result solid: both inputs, the operation and the
/// solid's index.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Polyhedron {
    pub(super) a: Box<Solid>,
    pub(super) b: Box<Solid>,
    pub(super) op: Op2,
    pub(super) index: usize,
}

/// One result solid: its parts, each slot's provenance and its bounds'
/// height range in the object's frame.
pub(super) struct Component {
    pub(super) parts: TopologyParts,
    pub(super) plans: Vec<SlotPlan>,
}

impl Polyhedron {
    pub(crate) fn tolerance(&self) -> Tolerance {
        self.a.resolution()
    }

    /// The same solid built again from its inputs.
    pub(crate) fn rebuilt(&self, operation: OperationId) -> Result<Solid> {
        let moved = self.clone();
        let mut components = build(&moved)?;
        if moved.index >= components.len() {
            return Err(Error::InvalidTopology(
                "a polyhedral Boolean rebuilt differently",
            ));
        }
        let component = components.swap_remove(moved.index);
        moved.solid(component, operation, None)
    }

    /// The solid moved rigidly: its stored geometry moved (vertices, lines
    /// and planes; the planes' pcurves unchanged), its enclosures measured
    /// again, its inputs' frames moved for classification.
    pub(crate) fn moved(
        &self,
        topology: &Topology,
        operation: OperationId,
        motion: crate::RigidTransform,
        mass: MassProperties,
    ) -> Result<Solid> {
        let parts = topology.moved_parts(motion, self.tolerance())?;
        let mut moved = self.clone();
        *moved.a = self.a.transform_with(self.a.operation, motion)?.0;
        *moved.b = self.b.transform_with(self.b.operation, motion)?.0;
        moved.solid(
            Component {
                parts,
                plans: Vec::new(),
            },
            operation,
            Some(mass),
        )
    }

    /// A component as a solid (external ids).
    pub(super) fn solid(
        &self,
        component: Component,
        operation: OperationId,
        known: Option<MassProperties>,
    ) -> Result<Solid> {
        let topology =
            Topology::from_parts(component.parts.with_measured_enclosures(), self.tolerance())
                .map_err(|issues| {
                    // A containment the validator's rays cannot decide (a
                    // cavity in a solid bounded by a sphere or torus face
                    // with loops, S9d.4b.1) is undecided, not invalid.
                    if !issues.is_empty()
                        && issues
                            .iter()
                            .all(|i| i.kind == crate::topology::IssueKind::UncertifiedContainment)
                    {
                        return Error::ComputationLimit(
                            "a cavity's containment the validator's rays leave undecided",
                        );
                    }
                    Error::InvalidTopology(
                        issues
                            .first()
                            .map_or("a polyhedral Boolean", |i| i.kind.name()),
                    )
                })?;
        let mass = match known {
            Some(mass) => mass,
            None => topology
                .mass_enclosure()
                .ok_or(Error::Unrepresentable(
                    "a polyhedral Boolean's mass properties",
                ))?
                .midpoints(),
        };
        let mut bounds = crate::solid::split::edge_bounds(&topology);
        // A sphere's face bulges past its edges (S9d.1): its whole sphere's
        // box, and every vertex (a pole's vertex loop among them).
        for f in topology.faces() {
            if let crate::topology::Surface::Sphere { frame, radius } = &f.surface {
                let c = frame.origin().to_array();
                let r = radius * (1.0 + 4.0 * f64::EPSILON);
                let (mut lo, mut hi) = (bounds.min.to_array(), bounds.max.to_array());
                for i in 0..3 {
                    lo[i] = lo[i].min(c[i] - r);
                    hi[i] = hi[i].max(c[i] + r);
                }
                bounds.min = crate::Point3::new(lo[0], lo[1], lo[2]);
                bounds.max = crate::Point3::new(hi[0], hi[1], hi[2]);
            }
        }
        let (lo, hi) = (self.a.start.min(self.a.end), self.a.start.max(self.a.end));
        Ok(Solid {
            construction: Construction::Polyhedron(Box::new(self.clone())),
            frame: self.a.frame,
            start: lo,
            end: hi,
            topology,
            mass,
            bounds,
            operation,
        })
    }

    /// Inside where the operation's set function holds of both inputs'
    /// memberships (a point on an input's boundary taking both) and the
    /// point lies within this solid's bounds.
    pub(crate) fn classify(
        &self,
        point: Point3,
        bounds: crate::Bounds3,
        tolerance: Tolerance,
    ) -> Result<Location> {
        let tol = tolerance.linear();
        let [x, y, z] = point.to_array();
        let (lo, hi) = (bounds.min.to_array(), bounds.max.to_array());
        for i in 0..3 {
            let c = [x, y, z][i];
            if c < lo[i] - 2.0 * tol || c > hi[i] + 2.0 * tol {
                return Ok(Location::Outside);
            }
        }
        let options = |l: Location| match l {
            Location::Inside => vec![true],
            Location::Outside => vec![false],
            Location::Boundary => vec![false, true],
        };
        let (la, lb) = (
            options(self.a.classify(point)?),
            options(self.b.classify(point)?),
        );
        let mut seen = BTreeSet::new();
        for a in &la {
            for b in &lb {
                seen.insert(holds(self.op, *a, *b));
            }
        }
        Ok(match (seen.contains(&true), seen.contains(&false)) {
            (true, false) => Location::Inside,
            (true, true) => Location::Boundary,
            _ => Location::Outside,
        })
    }
}

fn holds(op: Op2, a: bool, b: bool) -> bool {
    match op {
        Op2::Fuse => a || b,
        Op2::Cut => a && !b,
        Op2::Common => a && b,
    }
}

// ------------------------------------------------------------------ models

/// Where a fragment's face lies, for joining fragments into faces: an exact
/// oriented plane (a construction's face), or an input face (a stored
/// face, whose triangles' planes differ within rounding) with whether the
/// fragment keeps its orientation.
#[derive(Debug, Clone, PartialEq)]
enum Surf {
    Plane(Box<(V, R)>),
    Face(EntityId, bool),
}

/// A boundary face piece of a model: its plane (outward), convex pieces
/// (outward cycles), the input face it lies on and how fragments of it join.
struct MFace {
    plane: Plane,
    pieces: Vec<Vec<V>>,
    id: EntityId,
    stored: bool,
}

/// How a model decides a point's side.
enum Inside {
    /// Convex cells with disjoint interiors, as planes (inside `<= 0`).
    Cells(Vec<Vec<Plane>>),
    /// A closed surface of triangles, by exact ray parity.
    Mesh(Vec<MeshTri>),
}

/// A model's triangle with its normal (twice its vector area) and its
/// bounding box in floats, widened past rounding (for filtering only).
struct MeshTri {
    lo: [f64; 3],
    hi: [f64; 3],
    /// The corners times their common denominator `den`, and twice the
    /// vector area of those (integers: the tests' signs without
    /// reductions).
    ints: [I3; 3],
    den: BigInt,
    n: I3,
}

type I3 = [BigInt; 3];

impl MeshTri {
    fn new(pts: Vec<V>) -> Self {
        let (lo, hi) = float_box(&pts);
        let den = pts
            .iter()
            .flatten()
            .fold(BigInt::from(1), |l, x| l.lcm(x.denom()));
        let ints: [I3; 3] = std::array::from_fn(|m| {
            std::array::from_fn(|i| (&pts[m][i] * R::from_integer(den.clone())).to_integer())
        });
        let n = icross(&isub(&ints[1], &ints[0]), &isub(&ints[2], &ints[0]));
        Self {
            lo,
            hi,
            ints,
            den,
            n,
        }
    }

    /// Each corner minus `q`, times `den` and `q`'s denominator `w`.
    fn about(&self, (q, w): &(I3, BigInt)) -> [I3; 3] {
        std::array::from_fn(|m| std::array::from_fn(|i| &self.ints[m][i] * w - &self.den * &q[i]))
    }
}

/// A point as integers over a common positive denominator.
fn homogeneous(p: &V) -> (I3, BigInt) {
    let w = p.iter().fold(BigInt::from(1), |l, x| l.lcm(x.denom()));
    let q = std::array::from_fn(|i| (&p[i] * R::from_integer(w.clone())).to_integer());
    (q, w)
}

fn isub(a: &I3, b: &I3) -> I3 {
    std::array::from_fn(|i| &a[i] - &b[i])
}

fn icross(a: &I3, b: &I3) -> I3 {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

fn idot(a: &I3, b: &I3) -> BigInt {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

/// A rational's value within a relative `2^-60` or so, from its leading
/// bits (for filters only: far cheaper than a correctly rounded value).
fn approx(x: &R) -> f64 {
    let top = |v: &num_bigint::BigInt| {
        let shift = v.bits().saturating_sub(64);
        let bits = (v.magnitude() >> shift).to_u64_digits();
        (bits.first().copied().unwrap_or(0) as f64, shift as i32)
    };
    let ((n, sn), (d, sd)) = (top(x.numer()), top(x.denom()));
    let v = n / d * f64::powi(2.0, sn - sd);
    if x.numer().sign() == num_bigint::Sign::Minus {
        -v
    } else {
        v
    }
}

fn approx3(p: &V) -> [f64; 3] {
    [approx(&p[0]), approx(&p[1]), approx(&p[2])]
}

/// A bounding box of exact points in floats, widened by a relative margin
/// far past the error of each coordinate.
fn float_box(pts: &[V]) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in pts {
        let r = approx3(p);
        for k in 0..3 {
            lo[k] = lo[k].min(r[k]);
            hi[k] = hi[k].max(r[k]);
        }
    }
    for k in 0..3 {
        let m = 1e-9 * (1.0 + lo[k].abs().max(hi[k].abs()));
        lo[k] -= m;
        hi[k] += m;
    }
    (lo, hi)
}

fn boxes_meet(a: &([f64; 3], [f64; 3]), b: &([f64; 3], [f64; 3])) -> bool {
    (0..3).all(|k| a.0[k] <= b.1[k] && b.0[k] <= a.1[k])
}

/// An input's exact model, its entities by id.
struct Model {
    faces: Vec<MFace>,
    inside: Inside,
    edges: Vec<(V, V, EntityId)>,
    vertices: Vec<(V, EntityId)>,
    region: EntityId,
    /// Each entity's operand and role, for naming.
    info: BTreeMap<EntityId, (Operand, Role)>,
}

type P2 = (R, R);
/// The region of CCW polygons (outer and holes) as convex CCW polygons:
/// strips between consecutive vertex `u`s, the spanning edges ordered at the
/// strip's middle and paired by parity.
fn trapezoids(polygons: &[Vec<P2>]) -> Vec<Vec<P2>> {
    let mut edges: Vec<(P2, P2)> = Vec::new();
    for poly in polygons {
        for i in 0..poly.len() {
            let (p, q) = (poly[i].clone(), poly[(i + 1) % poly.len()].clone());
            if p.0 != q.0 {
                edges.push(if p.0 < q.0 { (p, q) } else { (q, p) });
            }
        }
    }
    let us: Vec<R> = polygons
        .iter()
        .flatten()
        .map(|p| p.0.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let at =
        |e: &(P2, P2), u: &R| &e.0 .1 + (&e.1 .1 - &e.0 .1) * (u - &e.0 .0) / (&e.1 .0 - &e.0 .0);
    let two = R::from_integer(2.into());
    let mut out = Vec::new();
    for w in us.windows(2) {
        let (u0, u1) = (&w[0], &w[1]);
        let mid = (u0 + u1) / &two;
        let mut span: Vec<&(P2, P2)> = edges
            .iter()
            .filter(|e| e.0 .0 <= *u0 && e.1 .0 >= *u1)
            .collect();
        span.sort_by_key(|e| at(e, &mid));
        for pair in span.chunks(2) {
            let [lo, hi] = [pair[0], pair[1]];
            let pts = [
                (u0.clone(), at(lo, u0)),
                (u1.clone(), at(lo, u1)),
                (u1.clone(), at(hi, u1)),
                (u0.clone(), at(hi, u0)),
            ];
            let mut ring: Vec<P2> = Vec::new();
            for p in pts {
                if ring.last() != Some(&p) {
                    ring.push(p);
                }
            }
            if ring.len() > 1 && ring.first() == ring.last() {
                ring.pop();
            }
            if ring.len() >= 3 {
                out.push(ring);
            }
        }
    }
    out
}

fn model(solid: &Solid, op: Operand) -> Result<Model> {
    match &solid.construction {
        Construction::Prism(profile)
            if profile.boundaries().all(|b| {
                b.polygon_vertices().is_some()
                    || b.path_geometry()
                        .is_some_and(|(_, s)| s.iter().all(|s| matches!(s, Segment::Line)))
            }) =>
        {
            prism_model(solid, profile, op)
        }
        _ => stored_model(solid, op),
    }
}

/// A prism's construction model (S9b.1).
fn prism_model(solid: &Solid, profile: &Profile, op: Operand) -> Result<Model> {
    let keys = super::index(solid)?;
    let id = |at: At, what: What, b: usize, j: usize| -> Result<EntityId> {
        keys.get(&(at, what, b, j))
            .copied()
            .ok_or(Error::InvalidTopology("a prism's entity"))
    };
    let f = solid.frame;
    let (o, x, y, n) = (
        vq(f.origin() - Point3::ORIGIN),
        vq(f.x()),
        vq(f.y()),
        vq(f.normal()),
    );
    let (lo, hi) = (q(solid.start.min(solid.end)), q(solid.start.max(solid.end)));
    let point =
        |u: &R, v: &R, w: &R| add(&add(&o, &scale(&x, u)), &add(&scale(&y, v), &scale(&n, w)));
    let lift = |p: &P2, w: &R| point(&p.0, &p.1, w);
    let mut polys: Vec<Vec<P2>> = Vec::new();
    for boundary in profile.boundaries() {
        let pts: Vec<Point2> = match (boundary.polygon_vertices(), boundary.path_geometry()) {
            (Some(p), _) => p.to_vec(),
            (None, Some((p, segs))) if segs.iter().all(|s| matches!(s, Segment::Line)) => {
                p.to_vec()
            }
            _ => {
                return Err(Error::OutOfDomain(
                    "a Boolean of a prism with arcs and a solid other than a prism (S9c)",
                ))
            }
        };
        polys.push(pts.iter().map(|p| (q(p.x), q(p.y))).collect());
    }
    let traps = trapezoids(&polys);
    let (low_at, high_at) = (point(&zero(), &zero(), &lo), point(&zero(), &zero(), &hi));
    // A cap's points differ by `u x + v y`: its plane's normal is `x * y`
    // (the stored axes are not exactly orthogonal to `n`).
    let up = cross(&x, &y);
    let mut faces = vec![
        MFace {
            plane: Plane::through(neg(&up), &low_at),
            pieces: traps
                .iter()
                .map(|t| t.iter().rev().map(|p| lift(p, &lo)).collect())
                .collect(),
            id: id(At::Low, What::Cap, 0, 0)?,
            stored: false,
        },
        MFace {
            plane: Plane::through(up.clone(), &high_at),
            pieces: traps
                .iter()
                .map(|t| t.iter().map(|p| lift(p, &hi)).collect())
                .collect(),
            id: id(At::High, What::Cap, 0, 0)?,
            stored: false,
        },
    ];
    let (mut edges, mut vertices) = (Vec::new(), Vec::new());
    for (b, poly) in polys.iter().enumerate() {
        let hole = b > 0;
        let count = poly.len();
        for j in 0..count {
            let (p, e) = (&poly[j], &poly[(j + 1) % count]);
            let t = add(&scale(&x, &(&e.0 - &p.0)), &scale(&y, &(&e.1 - &p.1)));
            let outward = cross(&t, &n);
            let (a, c) = (lift(p, &lo), lift(e, &lo));
            let (a2, c2) = (lift(p, &hi), lift(e, &hi));
            let (normal, rect) = if hole {
                (
                    neg(&outward),
                    vec![c.clone(), a.clone(), a2.clone(), c2.clone()],
                )
            } else {
                (outward, vec![a.clone(), c.clone(), c2.clone(), a2.clone()])
            };
            faces.push(MFace {
                plane: Plane::through(normal, &a),
                pieces: vec![rect],
                id: id(At::Swept, What::Wall, b, j)?,
                stored: false,
            });
            edges.push((a.clone(), c.clone(), id(At::Low, What::Edge, b, j)?));
            edges.push((a2.clone(), c2.clone(), id(At::High, What::Edge, b, j)?));
            edges.push((a.clone(), a2.clone(), id(At::Swept, What::Vertical, b, j)?));
            vertices.push((a, id(At::Low, What::Vertex, b, j)?));
            vertices.push((a2, id(At::High, What::Vertex, b, j)?));
        }
    }
    let mut cells = Vec::new();
    for t in &traps {
        let mut planes = vec![
            Plane::through(neg(&up), &low_at),
            Plane::through(up.clone(), &high_at),
        ];
        for i in 0..t.len() {
            let (p, e) = (&t[i], &t[(i + 1) % t.len()]);
            let d = add(&scale(&x, &(&e.0 - &p.0)), &scale(&y, &(&e.1 - &p.1)));
            planes.push(Plane::through(cross(&d, &n), &lift(p, &lo)));
        }
        cells.push(planes);
    }
    Ok(Model {
        faces,
        inside: Inside::Cells(cells),
        edges,
        vertices,
        region: id(At::Swept, What::Region, 0, 0)?,
        info: roles(solid, op),
    })
}

/// Each entity's operand and role.
fn roles(solid: &Solid, op: Operand) -> BTreeMap<EntityId, (Operand, Role)> {
    let t = &solid.topology;
    t.ids()
        .filter_map(|(id, _)| t.derivation(id).map(|d| (id, (op, d.role))))
        .collect()
}

/// The sign of `(b - a) x (c - a)`: in binary64 when its error bound
/// decides it, else exactly.
fn orient(a: &(P2, [f64; 2]), b: &(P2, [f64; 2]), c: &(P2, [f64; 2])) -> Sign {
    let (fa, fb, fc) = (a.1, b.1, c.1);
    let (l, r) = (
        (fb[0] - fa[0]) * (fc[1] - fa[1]),
        (fb[1] - fa[1]) * (fc[0] - fa[0]),
    );
    let det = l - r;
    let bound = 1e-14 * (l.abs() + r.abs());
    if det > bound {
        return Sign::Plus;
    }
    if det < -bound {
        return Sign::Minus;
    }
    let (a, b, c) = (&a.0, &b.0, &c.0);
    let x = (&b.0 - &a.0) * (&c.1 - &a.1) - (&b.1 - &a.1) * (&c.0 - &a.0);
    if x > zero() {
        Sign::Plus
    } else if x < zero() {
        Sign::Minus
    } else {
        Sign::NoSign
    }
}

/// Twice a polygon's signed area, exactly.
fn area_2d(poly: &[P2]) -> R {
    (0..poly.len()).fold(zero(), |acc, i| {
        let (a, b) = (&poly[i], &poly[(i + 1) % poly.len()]);
        acc + &a.0 * &b.1 - &a.1 * &b.0
    })
}

/// A region of polygons (an outer boundary and holes, either orientation)
/// as triangles of its own vertices, counter-clockwise: holes bridged to
/// the boundary by segments meeting nothing else, then ears clipped, the
/// fattest first (no slivers where a fatter ear exists). None when the
/// triangles' area is not the region's exactly (a degenerate polygon).
fn ear_clip(polys: &[Vec<P2>]) -> Option<Vec<[P2; 3]>> {
    type Pt = (P2, [f64; 2]);
    let pt = |p: &P2| -> Pt { (p.clone(), [approx(&p.0), approx(&p.1)]) };
    let areas: Vec<R> = polys.iter().map(|p| area_2d(p)).collect();
    let outer = (0..polys.len()).max_by(|a, b| abs(&areas[*a]).cmp(&abs(&areas[*b])))?;
    let oriented = |i: usize, ccw: bool| -> Vec<Pt> {
        let mut v: Vec<Pt> = polys[i].iter().map(pt).collect();
        if (areas[i] > zero()) != ccw {
            v.reverse();
        }
        v
    };
    let mut poly = oriented(outer, true);
    let mut holes: Vec<Vec<Pt>> = (0..polys.len())
        .filter(|&i| i != outer)
        .map(|i| oriented(i, false))
        .collect();
    // Twice the region's area: the outer boundary's less the holes'.
    let total = abs(&areas[outer]) * R::from_integer(2.into())
        - areas.iter().fold(zero(), |acc, a| acc + abs(a));
    // Whether the segment `pq` meets the edge `ab` anywhere but where it
    // shares an end with it.
    let meets = |p: &Pt, q: &Pt, a: &Pt, b: &Pt| -> bool {
        let shared = |x: &Pt| x.0 == p.0 || x.0 == q.0;
        if shared(a) && shared(b) {
            return a.0 != b.0;
        }
        let (d1, d2) = (orient(p, q, a), orient(p, q, b));
        let (d3, d4) = (orient(a, b, p), orient(a, b, q));
        if d1 != Sign::NoSign && d1 == d2 || d3 != Sign::NoSign && d3 == d4 {
            return false;
        }
        if d1 == Sign::NoSign && d2 == Sign::NoSign {
            // Collinear: overlapping beyond a shared end.
            let key = |x: &Pt| {
                if p.0 .0 != q.0 .0 {
                    x.0 .0.clone()
                } else {
                    x.0 .1.clone()
                }
            };
            let (lo, hi) = {
                let (x, y) = (key(p), key(q));
                if x <= y {
                    (x, y)
                } else {
                    (y, x)
                }
            };
            let (ka, kb) = (key(a), key(b));
            let (alo, ahi) = if ka <= kb { (ka, kb) } else { (kb, ka) };
            return alo < hi && lo < ahi;
        }
        // Touching at an end of one: only a shared end is allowed.
        let touch_a = d1 == Sign::NoSign && !shared(a);
        let touch_b = d2 == Sign::NoSign && !shared(b);
        let touch_p = d3 == Sign::NoSign && !(p.0 == a.0 || p.0 == b.0);
        let touch_q = d4 == Sign::NoSign && !(q.0 == a.0 || q.0 == b.0);
        touch_a
            || touch_b
            || touch_p
            || touch_q
            || (d1 != d2
                && d3 != d4
                && d1 != Sign::NoSign
                && d2 != Sign::NoSign
                && d3 != Sign::NoSign
                && d4 != Sign::NoSign)
    };
    // Whether the direction from a vertex (between `prev` and `next`, the
    // region on the left) to `x` enters the region there.
    let enters = |prev: &Pt, v: &Pt, next: &Pt, x: &Pt| -> bool {
        let (left_in, left_out) = (
            orient(prev, v, x) == Sign::Plus,
            orient(v, next, x) == Sign::Plus,
        );
        if orient(prev, v, next) == Sign::Plus {
            left_in && left_out
        } else {
            left_in || left_out
        }
    };
    holes.sort_by(|a, b| {
        let top = |h: &Vec<Pt>| h.iter().map(|p| p.0 .0.clone()).max();
        top(b).cmp(&top(a))
    });
    for (hi, hole) in holes.iter().enumerate() {
        let mut pairs: Vec<(f64, usize, usize)> = Vec::new();
        for (j, m) in hole.iter().enumerate() {
            for (i, p) in poly.iter().enumerate() {
                let d = (m.1[0] - p.1[0]).powi(2) + (m.1[1] - p.1[1]).powi(2);
                pairs.push((d, i, j));
            }
        }
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        let bridge = pairs.into_iter().find(|&(_, i, j)| {
            let (p, m) = (&poly[i], &hole[j]);
            let (pn, hn) = (poly.len(), hole.len());
            if !enters(&poly[(i + pn - 1) % pn], p, &poly[(i + 1) % pn], m)
                || !enters(&hole[(j + hn - 1) % hn], m, &hole[(j + 1) % hn], p)
            {
                return false;
            }
            let edges = |c: &[Pt]| -> bool {
                (0..c.len()).any(|k| meets(p, m, &c[k], &c[(k + 1) % c.len()]))
            };
            !edges(&poly) && !holes[hi..].iter().any(|h| edges(h))
        })?;
        let (_, i, j) = bridge;
        let mut merged: Vec<Pt> = poly[..=i].to_vec();
        merged.extend(hole[j..].iter().cloned());
        merged.extend(hole[..=j].iter().cloned());
        merged.push(poly[i].clone());
        merged.extend(poly[i + 1..].iter().cloned());
        poly = merged;
    }
    // Ears, fattest first.
    let mut out: Vec<[P2; 3]> = Vec::new();
    while poly.len() > 3 {
        let n = poly.len();
        let mut order: Vec<(f64, usize)> = (0..n)
            .map(|v| {
                let (a, b, c) = (&poly[(v + n - 1) % n].1, &poly[v].1, &poly[(v + 1) % n].1);
                let area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
                let l = |x: &[f64; 2], y: &[f64; 2]| (x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2);
                (
                    area / (l(a, b) + l(b, c) + l(c, a)).max(f64::MIN_POSITIVE),
                    v,
                )
            })
            .collect();
        order.sort_by(|a, b| b.0.total_cmp(&a.0));
        let ear = order.into_iter().map(|(_, v)| v).find(|&v| {
            let (a, b, c) = (&poly[(v + n - 1) % n], &poly[v], &poly[(v + 1) % n]);
            if orient(a, b, c) != Sign::Plus {
                return false;
            }
            poly.iter().all(|x| {
                x.0 == a.0
                    || x.0 == b.0
                    || x.0 == c.0
                    || orient(a, b, x) == Sign::Minus
                    || orient(b, c, x) == Sign::Minus
                    || orient(c, a, x) == Sign::Minus
            })
        })?;
        out.push([
            poly[(ear + n - 1) % n].0.clone(),
            poly[ear].0.clone(),
            poly[(ear + 1) % n].0.clone(),
        ]);
        poly.remove(ear);
    }
    if orient(&poly[0], &poly[1], &poly[2]) == Sign::Plus {
        out.push([poly[0].0.clone(), poly[1].0.clone(), poly[2].0.clone()]);
    }
    let sum = out.iter().fold(zero(), |acc, t| acc + area_2d(t));
    (sum == total).then_some(out)
}

/// A polygon's trapezoids, each side at its strip's `u`s through every
/// corner there (lifted corners of a side are not collinear exactly, so the
/// neighbours' triangles meet at all of them), zipped into triangles
/// between the two sides.
fn zipped(polys: &[Vec<P2>]) -> Vec<[P2; 3]> {
    let mut tris: Vec<[P2; 3]> = Vec::new();
    let traps = trapezoids(polys);
    let mut at_u: BTreeMap<R, BTreeSet<R>> = BTreeMap::new();
    for p in traps.iter().flatten().chain(polys.iter().flatten()) {
        at_u.entry(p.0.clone()).or_default().insert(p.1.clone());
    }
    for trap in &traps {
        let u0 = trap.iter().map(|p| &p.0).min().expect("a corner");
        let u1 = trap.iter().map(|p| &p.0).max().expect("a corner");
        let side = |u: &R| -> Vec<P2> {
            let vs: Vec<&R> = trap.iter().filter(|p| &p.0 == u).map(|p| &p.1).collect();
            let (lo, hi) = (vs.iter().min().unwrap(), vs.iter().max().unwrap());
            at_u[u]
                .range((*lo).clone()..=(*hi).clone())
                .map(|v| (u.clone(), v.clone()))
                .collect()
        };
        let (left, right) = (side(u0), side(u1));
        let (mut i, mut j) = (0, 0);
        while i + 1 < left.len() || j + 1 < right.len() {
            if j + 1 < right.len() && (i + 1 == left.len() || right[j + 1].1 <= left[i + 1].1) {
                tris.push([left[i].clone(), right[j].clone(), right[j + 1].clone()]);
                j += 1;
            } else {
                tris.push([left[i].clone(), right[j].clone(), left[i + 1].clone()]);
                i += 1;
            }
        }
    }
    tris
}

/// A planar body's stored model (S9b.2): its stored vertices as rationals,
/// each face as the triangles of its trapezoids in its projection on the
/// normal's largest coordinate plane (their corners exact points of the
/// face's straight edges), so each triangle is planar exactly; a face's
/// fragments join back by the face. Its side by exact ray parity.
fn stored_model(solid: &Solid, op: Operand) -> Result<Model> {
    let other = || {
        Error::OutOfDomain("a Boolean of a solid with curved faces or edges in any position (S9c)")
    };
    let t = &solid.topology;
    let point = |p: Point3| [q(p.x), q(p.y), q(p.z)];
    let id_of = |slot: Slot| {
        t.id_of(slot)
            .ok_or(Error::InvalidTopology("an unnamed slot"))
    };
    let mut faces = Vec::new();
    let mut triangles: Vec<Vec<V>> = Vec::new();
    for (fi, face) in t.faces().iter().enumerate() {
        let Surface::Plane(frame) = &face.surface else {
            return Err(other());
        };
        let normal = frame.normal() * face.sense.sign();
        let nv = vq(normal);
        // Loops as cycles of exact points.
        let mut cycles: Vec<Vec<V>> = Vec::new();
        for l in &face.loops {
            let Loop::Edges { fins, .. } = &t.loops()[l.index()] else {
                return Err(other());
            };
            let mut cycle = Vec::new();
            for f in fins {
                let fin = &t.fins()[f.index()];
                let edge = &t.edges()[fin.edge.index()];
                if !matches!(edge.curve, Curve3::LineSegment { .. }) {
                    return Err(other());
                }
                let v = if fin.sense == Orientation::Forward {
                    edge.start
                } else {
                    edge.end
                }
                .ok_or_else(other)?;
                cycle.push(point(t.vertices()[v.index()].position));
            }
            cycles.push(cycle);
        }
        // Project on the normal's largest coordinate plane.
        let k = (0..3)
            .max_by(|a, b| abs(&nv[*a]).cmp(&abs(&nv[*b])))
            .expect("a coordinate");
        let (i, j) = [(1, 2), (2, 0), (0, 1)][k];
        let flip = nv[k] < zero();
        let mut polys: Vec<Vec<P2>> = Vec::new();
        let mut lift: BTreeMap<(P2, P2), (V, V)> = BTreeMap::new();
        for cycle in &cycles {
            let poly: Vec<P2> = cycle.iter().map(|p| (p[i].clone(), p[j].clone())).collect();
            for m in 0..cycle.len() {
                let n2 = (m + 1) % cycle.len();
                let (a, b) = (poly[m].clone(), poly[n2].clone());
                let key = if a <= b { (a, b) } else { (b, a) };
                let ends = if poly[m] <= poly[n2] {
                    (cycle[m].clone(), cycle[n2].clone())
                } else {
                    (cycle[n2].clone(), cycle[m].clone())
                };
                lift.insert(key, ends);
            }
            polys.push(poly);
        }
        // A projected point on a projected edge, lifted along its 3D edge.
        let lifted = |p: &P2| -> V {
            for ((a, b), (a3, b3)) in &lift {
                let (du, dv) = (&b.0 - &a.0, &b.1 - &a.1);
                let (pu, pv) = (&p.0 - &a.0, &p.1 - &a.1);
                if &du * &pv - &dv * &pu != zero() {
                    continue;
                }
                let len = &du * &du + &dv * &dv;
                let t = (&pu * &du + &pv * &dv) / &len;
                if t >= zero() && t <= R::from_integer(1.into()) {
                    return add(a3, &scale(&sub(b3, a3), &t));
                }
            }
            unreachable!("a trapezoid's corner on the face's boundary")
        };
        let id = id_of(Slot::Face(FaceId(fi)))?;
        // The face's polygon cut into triangles of its own vertices (ears
        // clipped, holes bridged), so neighbouring faces' triangles share
        // their stored edges exactly; else (the clipping's area check
        // failing) its trapezoids zipped into triangles.
        let tris = match ear_clip(&polys) {
            Some(tris) => tris,
            None => zipped(&polys),
        };
        for t2 in tris {
            let mut tri: Vec<V> = t2.iter().map(lifted).collect();
            if flip {
                tri.reverse();
            }
            let n = area2(&tri);
            faces.push(MFace {
                plane: Plane::through(n, &tri[0]),
                pieces: vec![tri.clone()],
                id,
                stored: true,
            });
            triangles.push(tri);
        }
    }
    let mut edges = Vec::new();
    for (ei, e) in t.edges().iter().enumerate() {
        let Curve3::LineSegment { start, end } = e.curve else {
            return Err(other());
        };
        edges.push((point(start), point(end), id_of(Slot::Edge(EdgeId(ei)))?));
    }
    let vertices = t
        .vertices()
        .iter()
        .enumerate()
        .map(|(vi, v)| Ok((point(v.position), id_of(Slot::Vertex(VertexId(vi)))?)))
        .collect::<Result<Vec<_>>>()?;
    Ok(Model {
        faces,
        inside: Inside::Mesh(triangles.into_iter().map(MeshTri::new).collect()),
        edges,
        vertices,
        region: id_of(Slot::Region(RegionId(1)))?,
        info: roles(solid, op),
    })
}

/// A point strictly inside a convex planar cycle on `plane`: one of short
/// dyadic coordinates near its centroid (its last coordinate solved on the
/// plane) when one lies inside, else the centroid. Short coordinates keep
/// the classification's exact arithmetic small.
fn interior(poly: &[V], plane: &Plane) -> V {
    let n = &plane.n;
    let k = (0..3)
        .max_by(|a, b| abs(&n[*a]).cmp(&abs(&n[*b])))
        .expect("a coordinate");
    let count = poly.len() as f64;
    let centre: Vec<f64> = (0..3)
        .map(|i| poly.iter().map(|p| approx(&p[i])).sum::<f64>() / count)
        .collect();
    let normal = area2(poly);
    let strictly_inside = |x: &V| {
        (0..poly.len()).all(|m| {
            let (a, b) = (&poly[m], &poly[(m + 1) % poly.len()]);
            dot(&cross(&sub(b, a), &sub(x, a)), &normal) > zero()
        })
    };
    // Grids of 2^-bits of the size (a power of two).
    let size = centre
        .iter()
        .fold(1.0f64, |m, x| m.max(x.abs()))
        .log2()
        .ceil() as i32;
    for bits in [20i32, 36, 52] {
        let steps = R::from_integer(2.into()).pow(bits - size);
        let mut x: V = [zero(), zero(), zero()];
        for i in 0..3 {
            if i != k {
                let units = (centre[i] * f64::powi(2.0, bits - size)).round() as i64;
                x[i] = R::from_integer(units.into()) / &steps;
            }
        }
        let rest =
            n[(k + 1) % 3].clone() * &x[(k + 1) % 3] + n[(k + 2) % 3].clone() * &x[(k + 2) % 3];
        x[k] = -(rest + &plane.d) / &n[k];
        if strictly_inside(&x) {
            return x;
        }
    }
    // Else (a sliver) the centroid of the fan's triangle of the shortest
    // corners: inside the triangle, so inside the cycle.
    let bits = |p: &V| -> u64 { p.iter().map(|x| x.numer().bits() + x.denom().bits()).sum() };
    let third = R::new(1.into(), 3.into());
    (1..poly.len() - 1)
        .filter(|&m| {
            !is_zero(&area2(&[
                poly[0].clone(),
                poly[m].clone(),
                poly[m + 1].clone(),
            ]))
        })
        .min_by_key(|&m| bits(&poly[m]) + bits(&poly[m + 1]))
        .map(|m| scale(&add(&add(&poly[0], &poly[m]), &poly[m + 1]), &third))
        .expect("a cycle of positive area")
}

/// Whether `p + e d` lies in the model for every small `e > 0` (never on
/// the model's boundary here: every plane of the other's faces splits the
/// fragments).
fn inside(model: &Model, p: &V, d: &V) -> bool {
    match &model.inside {
        Inside::Cells(cells) => cells.iter().any(|cell| {
            cell.iter().all(|plane| {
                let s = plane.eval(p);
                s < zero() || (s == zero() && dot(&plane.n, d) <= zero())
            })
        }),
        Inside::Mesh(triangles) => {
            // The push along a short dyadic direction on `d`'s side (only
            // the planes holding the fragment hold `p`: they split it),
            // shortened until it crosses no plane of a triangle near it.
            let df = approx3(d);
            let len = df.iter().map(|x| x * x).sum::<f64>().sqrt();
            let grid = R::from_integer((1i64 << 20).into());
            let short: V =
                df.map(|x| R::from_integer(((x / len * 1048576.0).round() as i64).into()) / &grid);
            let dir = if len.is_finite() && len > 0.0 && dot(&short, d) > zero() {
                short
            } else {
                d.clone()
            };
            let pf = approx3(p);
            let size = pf.iter().fold(1.0f64, |m, x| m.max(x.abs()));
            let mut step = R::new(1.into(), (1i64 << 20).into())
                * R::from_integer(((size as i64).max(1) + 1).into());
            let hp = homogeneous(p);
            let q = loop {
                let q = add(p, &scale(&dir, &step));
                let hq = homogeneous(&q);
                let reach = float_box(&[p.clone(), q.clone()]);
                let crosses = triangles.iter().any(|t| {
                    if !boxes_meet(&reach, &(t.lo, t.hi)) {
                        return false;
                    }
                    // The corner's side of the plane from each point.
                    let (s0, s1) = (
                        idot(&t.n, &t.about(&hp)[0]).sign(),
                        idot(&t.n, &t.about(&hq)[0]).sign(),
                    );
                    s1 == Sign::NoSign || (s0 != Sign::NoSign && s0 != s1)
                });
                if !crosses {
                    break q;
                }
                step *= R::new(1.into(), 1024.into());
            };
            parity(triangles, &q).unwrap_or(false)
        }
    }
}

/// Whether a point off a closed surface of triangles lies inside it: the
/// parity of a ray's crossings (the ray's sides of the triangle's edges
/// agreeing and the plane ahead), retried in other directions when the ray
/// meets an edge or a triangle's plane holds it. Triangles whose boxes the
/// ray misses are skipped in floats.
fn parity(triangles: &[MeshTri], q: &V) -> Option<bool> {
    let dirs: [[i64; 3]; 5] = [[7, 3, 5], [2, 11, 13], [17, 5, 3], [3, 19, 7], [23, 29, 31]];
    let qf = approx3(q);
    let hq = homogeneous(q);
    'dirs: for d in dirs {
        let df = d.map(|x| x as f64);
        let d: I3 = d.map(BigInt::from);
        let mut count = 0usize;
        for tri in triangles {
            // The ray's parameters within the box (every direction positive).
            let (mut t0, mut t1) = (0.0f64, f64::INFINITY);
            for k in 0..3 {
                t0 = t0.max((tri.lo[k] - qf[k]) / df[k]);
                t1 = t1.min((tri.hi[k] - qf[k]) / df[k]);
            }
            if t0 > t1 + 1e-9 * (1.0 + t1.abs()) {
                continue;
            }
            // Signs of the integer forms: the plane's value at the first
            // corner from `q` and the ray's sides of the edges.
            let rel = tri.about(&hq);
            let denom = idot(&tri.n, &d).sign();
            let s = idot(&tri.n, &rel[0]).sign();
            let sides: Vec<Sign> = (0..3)
                .map(|m| idot(&d, &icross(&rel[m], &rel[(m + 1) % 3])).sign())
                .collect();
            let pos = sides.contains(&Sign::Plus);
            let neg = sides.contains(&Sign::Minus);
            if pos && neg {
                continue;
            }
            let none = Sign::NoSign;
            // The ray's line meets the closed triangle (or lies in its plane).
            if denom != none && s != none && s != denom {
                continue;
            }
            if denom == none || s == none || sides.contains(&none) {
                if denom == none && s != none {
                    continue;
                }
                continue 'dirs;
            }
            count += 1;
        }
        return Some(count % 2 == 1);
    }
    None
}

/// A convex cycle's two sides of a plane (the part where the plane's value
/// has `keep`'s sign or vanishes).
fn clip(poly: &[V], plane: &Plane, keep: i8) -> Vec<V> {
    let s: Vec<R> = poly
        .iter()
        .map(|p| {
            let v = plane.eval(p);
            if keep < 0 {
                v
            } else {
                -v
            }
        })
        .collect();
    let mut out = Vec::new();
    for i in 0..poly.len() {
        let j = (i + 1) % poly.len();
        let (sp, sq) = (&s[i], &s[j]);
        if *sp <= zero() {
            out.push(poly[i].clone());
        }
        if (*sp < zero() && *sq > zero()) || (*sq < zero() && *sp > zero()) {
            let t = sp / (sp - sq);
            out.push(add(&poly[i], &scale(&sub(&poly[j], &poly[i]), &t)));
        }
    }
    let mut dedup: Vec<V> = Vec::new();
    for p in out {
        if dedup.last() != Some(&p) {
            dedup.push(p);
        }
    }
    while dedup.len() > 1 && dedup.first() == dedup.last() {
        dedup.pop();
    }
    dedup
}

/// Twice a planar cycle's vector area.
fn area2(poly: &[V]) -> V {
    let mut total = [zero(), zero(), zero()];
    for i in 1..poly.len().saturating_sub(1) {
        total = add(
            &total,
            &cross(&sub(&poly[i], &poly[0]), &sub(&poly[i + 1], &poly[0])),
        );
    }
    total
}

/// A kept fragment: its cycle (outward for the result), its source face and
/// whether it keeps that face's orientation.
struct Frag {
    pts: Vec<V>,
    src: (Operand, usize),
    kept_way: bool,
}

// ------------------------------------------------------------------ build

/// Builds a polyhedral Boolean's components in its inputs' frames.
#[allow(clippy::too_many_lines)]
pub(super) fn build(poly: &Polyhedron) -> Result<Vec<Component>> {
    // S9c.1: prisms with arcs in any position.
    if super::curved::applies(poly) {
        return super::curved::build(poly);
    }
    let tolerance = poly.tolerance();
    let models = [model(&poly.a, Operand::A)?, model(&poly.b, Operand::B)?];
    // Each input entity's operand and role.
    let info: BTreeMap<EntityId, (Operand, Role)> = models
        .iter()
        .flat_map(|m| m.info.iter().map(|(k, v)| (*k, *v)))
        .collect();
    let op = poly.op;
    // Each input's faces' planes, unoriented and distinct, with their faces'
    // bounding box (a part away from it meets none of them).
    type Boxed = (Plane, ([f64; 3], [f64; 3]));
    let planes: Vec<Vec<Boxed>> = models
        .iter()
        .map(|m| {
            let mut seen: BTreeMap<(V, R), usize> = BTreeMap::new();
            let mut out: Vec<Boxed> = Vec::new();
            for f in &m.faces {
                let (n, d) = f.plane.canonical(false);
                let b = float_box(&f.pieces.concat());
                match seen.get(&(n.clone(), d.clone())) {
                    Some(&i) => {
                        let e = &mut out[i].1;
                        for k in 0..3 {
                            e.0[k] = e.0[k].min(b.0[k]);
                            e.1[k] = e.1[k].max(b.1[k]);
                        }
                    }
                    None => {
                        seen.insert((n.clone(), d.clone()), out.len());
                        out.push((Plane { n, d }, b));
                    }
                }
            }
            out
        })
        .collect();
    let on_planes: Vec<BTreeSet<(V, R)>> = planes
        .iter()
        .map(|ps| ps.iter().map(|(p, _)| (p.n.clone(), p.d.clone())).collect())
        .collect();
    // Fragments, classified and kept.
    let mut frags: Vec<Frag> = Vec::new();
    for (k, m) in models.iter().enumerate() {
        let other = &models[1 - k];
        let operand = if k == 0 { Operand::A } else { Operand::B };
        for (fi, face) in m.faces.iter().enumerate() {
            let n = &face.plane.n;
            let mut parts: Vec<(Vec<V>, _)> = face
                .pieces
                .iter()
                .map(|p| (p.clone(), float_box(p)))
                .collect();
            for (plane, reach) in &planes[1 - k] {
                let mut next = Vec::new();
                for (part, bounds) in parts {
                    if !boxes_meet(&bounds, reach) {
                        next.push((part, bounds));
                        continue;
                    }
                    let s: Vec<R> = part.iter().map(|p| plane.eval(p)).collect();
                    let crosses = s.iter().any(|x| *x < zero()) && s.iter().any(|x| *x > zero());
                    if crosses {
                        for keep in [-1i8, 1] {
                            let g = clip(&part, plane, keep);
                            // Positive area: the plane crosses the part.
                            if g.len() >= 3 {
                                let b = float_box(&g);
                                next.push((g, b));
                            }
                        }
                    } else {
                        next.push((part, bounds));
                    }
                }
                parts = next;
            }
            let flat = face.plane.canonical(false);
            let coplanar = on_planes[1 - k].contains(&flat);
            for (part, _) in parts {
                let c = interior(&part, &face.plane);
                // Off the other's planes a fragment's point decides both
                // sides at once (its near planes split it).
                let (front, back) = match &other.inside {
                    Inside::Mesh(triangles) if !coplanar => {
                        let v = parity(triangles, &c).unwrap_or(false);
                        (v, v)
                    }
                    _ => (inside(other, &c, n), inside(other, &c, &neg(n))),
                };
                let (keep, behind) = if k == 0 {
                    (
                        holds(op, false, front) != holds(op, true, back),
                        holds(op, true, back),
                    )
                } else {
                    (
                        front == back && holds(op, front, false) != holds(op, back, true),
                        holds(op, back, true),
                    )
                };
                if !keep {
                    continue;
                }
                let mut pts = part;
                if !behind {
                    pts.reverse();
                }
                frags.push(Frag {
                    pts,
                    src: (operand, fi),
                    kept_way: behind,
                });
            }
        }
    }
    if frags.is_empty() {
        return Ok(Vec::new());
    }
    // Conforming: every kept vertex on an edge's interior splits it.
    let mut ids: BTreeMap<V, usize> = BTreeMap::new();
    let mut points: Vec<V> = Vec::new();
    for f in &frags {
        for p in &f.pts {
            if !ids.contains_key(p) {
                ids.insert(p.clone(), points.len());
                points.push(p.clone());
            }
        }
    }
    let approx: Vec<[f64; 3]> = points.iter().map(approx3).collect();
    let homs: Vec<(I3, BigInt)> = points.iter().map(homogeneous).collect();
    let tol = tolerance.linear();
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for f in &frags {
        let mut cycle = Vec::new();
        for i in 0..f.pts.len() {
            let (a, b) = (&f.pts[i], &f.pts[(i + 1) % f.pts.len()]);
            let (ia, ib) = (ids[a], ids[b]);
            cycle.push(ia);
            let ab = sub(b, a);
            let (fa, fb) = (approx[ia], approx[ib]);
            let mut on: Vec<(R, usize)> = Vec::new();
            for (iv, fv) in approx.iter().enumerate() {
                if iv == ia || iv == ib {
                    continue;
                }
                let outside = (0..3)
                    .any(|d| fv[d] < fa[d].min(fb[d]) - tol || fv[d] > fa[d].max(fb[d]) + tol);
                if outside {
                    continue;
                }
                // Far off the line in floats (past the error of the
                // approximations, relative 2^-60) needs no exact test.
                let (u, w) = (
                    [fb[0] - fa[0], fb[1] - fa[1], fb[2] - fa[2]],
                    [fv[0] - fa[0], fv[1] - fa[1], fv[2] - fa[2]],
                );
                let size = (0..3).fold(1.0f64, |m, d| m.max(fa[d].abs()).max(fv[d].abs()));
                let reach =
                    1e-12 * size * (u.iter().chain(&w).map(|x| x.abs()).sum::<f64>() + size);
                let c = [
                    u[1] * w[2] - u[2] * w[1],
                    u[2] * w[0] - u[0] * w[2],
                    u[0] * w[1] - u[1] * w[0],
                ];
                if c.iter().any(|x| x.abs() > reach) {
                    continue;
                }
                // Exactly, in integers over the points' denominators.
                let ((ia3, wa), (ib3, wb), (iv3, wv)) = (&homs[ia], &homs[ib], &homs[iv]);
                let scaled = |x: &I3, wx: &BigInt, y: &I3, wy: &BigInt| -> I3 {
                    std::array::from_fn(|d| wy * &x[d] - wx * &y[d])
                };
                let (iab, iav) = (scaled(ib3, wb, ia3, wa), scaled(iv3, wv, ia3, wa));
                if icross(&iab, &iav).iter().any(|x| x.sign() != Sign::NoSign) {
                    continue;
                }
                let ibv = scaled(iv3, wv, ib3, wb);
                if idot(&iav, &iab).sign() == Sign::Plus && idot(&ibv, &iab).sign() == Sign::Minus {
                    on.push((dot(&sub(&points[iv], a), &ab), iv));
                }
            }
            on.sort();
            cycle.extend(on.into_iter().map(|(_, iv)| iv));
        }
        cycles.push(cycle);
    }
    // Faces: fragments of one oriented plane joined across shared edges.
    let face_of = |(o, fi): (Operand, usize)| &models[usize::from(o == Operand::B)].faces[fi];
    let frag_plane: Vec<Surf> = frags
        .iter()
        .map(|f| {
            let face = face_of(f.src);
            if face.stored {
                Surf::Face(face.id, f.kept_way)
            } else {
                let (n, d) = Plane::through(area2(&f.pts), &f.pts[0]).canonical(true);
                Surf::Plane(Box::new((n, d)))
            }
        })
        .collect();
    let mut uses: BTreeMap<(usize, usize), Vec<(usize, bool)>> = BTreeMap::new();
    for (fi, cycle) in cycles.iter().enumerate() {
        for i in 0..cycle.len() {
            let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
            uses.entry((a.min(b), a.max(b)))
                .or_default()
                .push((fi, a < b));
        }
    }
    let mut parent: Vec<usize> = (0..frags.len()).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut y = x;
        while p[y] != r {
            let n = p[y];
            p[y] = r;
            y = n;
        }
        r
    }
    for u in uses.values() {
        for i in 0..u.len() {
            for j in (i + 1)..u.len() {
                let ((f1, d1), (f2, d2)) = (u[i], u[j]);
                if d1 != d2 && frag_plane[f1] == frag_plane[f2] {
                    let (r1, r2) = (find(&mut parent, f1), find(&mut parent, f2));
                    if r1 != r2 {
                        parent[r1.max(r2)] = r1.min(r2);
                    }
                }
            }
        }
    }
    let mut group_of = vec![0usize; frags.len()];
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut root_group: BTreeMap<usize, usize> = BTreeMap::new();
    for (fi, g) in group_of.iter_mut().enumerate() {
        let r = find(&mut parent, fi);
        let gi = *root_group.entry(r).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        groups[gi].push(fi);
        *g = gi;
    }
    // Remaining uses: a face's own shared edges cancel.
    let mut kept_uses: BTreeMap<(usize, usize), Vec<(usize, bool)>> = BTreeMap::new();
    for (e, u) in &uses {
        let mut rest: Vec<(usize, bool)> = Vec::new();
        for &(f, d) in u {
            if let Some(k) = rest
                .iter()
                .position(|&(g, dd)| group_of[g] == group_of[f] && dd != d)
            {
                rest.remove(k);
            } else {
                rest.push((f, d));
            }
        }
        match rest.len() {
            0 => {}
            2 if rest[0].1 != rest[1].1 && group_of[rest[0].0] != group_of[rest[1].0] => {
                kept_uses.insert(*e, rest);
            }
            2 => return Err(Error::Degenerate("a face meeting itself along an edge")),
            4 => return Err(Error::Degenerate("a result touching itself along an edge")),
            _ => return Err(Error::InvalidTopology("an open polyhedral Boolean")),
        }
    }
    // Vertices where the remaining edges do not run straight on between the
    // same two faces.
    let mut incident: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
    for e in kept_uses.keys() {
        incident.entry(e.0).or_default().push(*e);
        incident.entry(e.1).or_default().push(*e);
    }
    let faces_of = |e: &(usize, usize)| -> BTreeSet<usize> {
        kept_uses[e].iter().map(|u| group_of[u.0]).collect()
    };
    let removable = |v: usize| -> bool {
        let es = &incident[&v];
        if es.len() != 2 || faces_of(&es[0]) != faces_of(&es[1]) {
            return false;
        }
        let other = |e: &(usize, usize)| if e.0 == v { e.1 } else { e.0 };
        let (a, b) = (
            sub(&points[other(&es[0])], &points[v]),
            sub(&points[other(&es[1])], &points[v]),
        );
        is_zero(&cross(&a, &b)) && dot(&a, &b) < zero()
    };
    let mut parts = TopologyParts::default();
    let mut vertex_of: BTreeMap<usize, VertexId> = BTreeMap::new();
    for &v in incident.keys() {
        if !removable(v) {
            vertex_of.insert(v, VertexId(parts.vertices.len()));
            parts.vertices.push(Vertex {
                position: rounded(&points[v]),
                enclosure: None,
            });
        }
    }
    let mut edge_of: BTreeMap<(usize, usize), (EdgeId, bool)> = BTreeMap::new();
    let mut chains: Vec<Vec<(usize, usize)>> = Vec::new();
    // Each edge's end points.
    let mut edge_ends: Vec<(usize, usize)> = Vec::new();
    for &start in vertex_of.keys() {
        for first in incident[&start].clone() {
            if edge_of.contains_key(&first) {
                continue;
            }
            let id = EdgeId(chains.len());
            let mut chain = Vec::new();
            let (mut at, mut e) = (start, first);
            loop {
                let next = if e.0 == at { e.1 } else { e.0 };
                edge_of.insert(e, (id, e.0 == at));
                chain.push(e);
                at = next;
                if vertex_of.contains_key(&at) {
                    break;
                }
                e = *incident[&at]
                    .iter()
                    .find(|x| **x != e)
                    .ok_or(Error::InvalidTopology("an open polyhedral edge"))?;
            }
            edge_ends.push((start, at));
            parts.edges.push(Edge {
                start: Some(vertex_of[&start]),
                end: Some(vertex_of[&at]),
                curve: Curve3::LineSegment {
                    start: rounded(&points[start]),
                    end: rounded(&points[at]),
                },
                fins: Vec::new(),
            });
            chains.push(chain);
        }
    }
    if edge_of.len() != kept_uses.len() {
        return Err(Error::Degenerate("a closed edge without a vertex"));
    }
    // Faces' loops.
    struct Built {
        surface: Surface,
        loops: Vec<Vec<Fin>>,
    }
    let mut built: Vec<Built> = Vec::new();
    for group in &groups {
        // The face's vector area (a stored face's fragments lie on its
        // triangles' planes, some of them slivers) and its point farthest
        // from the first for the axis.
        let f0 = &frags[group[0]];
        let normal = group.iter().fold([zero(), zero(), zero()], |acc, &fi| {
            add(&acc, &area2(&frags[fi].pts))
        });
        let nf = Vec3::new(
            rational_f64(&normal[0]),
            rational_f64(&normal[1]),
            rational_f64(&normal[2]),
        );
        let origin = rounded(&f0.pts[0]);
        let far = group
            .iter()
            .flat_map(|&fi| frags[fi].pts.iter())
            .map(|p| rounded(p) - origin)
            .max_by(|a, b| a.length().total_cmp(&b.length()))
            .expect("a point");
        let hint = far;
        let frame = Frame3::new(origin, nf, hint, tolerance)?;
        let mut used: BTreeMap<EdgeId, Orientation> = BTreeMap::new();
        for &fi in group {
            let cycle = &cycles[fi];
            for i in 0..cycle.len() {
                let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
                let key = (a.min(b), a.max(b));
                if !kept_uses.contains_key(&key) {
                    continue;
                }
                let (e, along) = edge_of[&key];
                let o = if (key.0 == a) == along {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                if used.insert(e, o).is_some_and(|x| x != o) {
                    return Err(Error::Degenerate("a face using an edge both ways"));
                }
            }
        }
        let ends = |e: EdgeId, o: Orientation| {
            let edge = &parts.edges[e.0];
            let (s, t) = (edge.start.expect("a vertex"), edge.end.expect("a vertex"));
            if o == Orientation::Forward {
                (s, t)
            } else {
                (t, s)
            }
        };
        let mut starts: BTreeMap<VertexId, EdgeId> = BTreeMap::new();
        for (&e, &o) in &used {
            if starts.insert(ends(e, o).0, e).is_some() {
                return Err(Error::Degenerate("a face touching itself at a vertex"));
            }
        }
        // A face thinner than the resolution: one of its vertices within it
        // of another vertex or of an edge not ending there (exactly, on the
        // model).
        let tol2 = q(tolerance.linear()) * q(tolerance.linear());
        let corners: BTreeSet<usize> = used
            .keys()
            .flat_map(|e| [edge_ends[e.0].0, edge_ends[e.0].1])
            .collect();
        // Pairs apart by twice the resolution in binary64 (past the
        // approximations' error) need no exact test.
        let tol_f = tolerance.linear();
        let far = |p: [f64; 3], a: [f64; 3], b: [f64; 3]| -> bool {
            let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
            let len2 = ab.iter().map(|x| x * x).sum::<f64>();
            let t = if len2 > 0.0 {
                (ap.iter().zip(&ab).map(|(x, y)| x * y).sum::<f64>() / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let d2 = (0..3).map(|k| (ap[k] - t * ab[k]).powi(2)).sum::<f64>();
            let size = (0..3).fold(1.0f64, |m, k| {
                m.max(p[k].abs()).max(a[k].abs()).max(b[k].abs())
            });
            d2.sqrt() > 2.0 * tol_f + 1e-12 * size
        };
        for &v in &corners {
            let p = &points[v];
            for &w in &corners {
                if w > v && !far(approx[v], approx[w], approx[w]) {
                    let d = sub(&points[w], p);
                    if dot(&d, &d) <= tol2 {
                        return Err(Error::Degenerate("a face thinner than the resolution"));
                    }
                }
            }
            for e in used.keys() {
                let (a, b) = edge_ends[e.0];
                if a == v || b == v || far(approx[v], approx[a], approx[b]) {
                    continue;
                }
                let (pa, ab) = (sub(p, &points[a]), sub(&points[b], &points[a]));
                let t = dot(&pa, &ab);
                let len2 = dot(&ab, &ab);
                let d2 = if t <= zero() {
                    dot(&pa, &pa)
                } else if t >= len2 {
                    let pb = sub(p, &points[b]);
                    dot(&pb, &pb)
                } else {
                    let c = cross(&pa, &ab);
                    dot(&c, &c) / &len2
                };
                if d2 <= tol2 {
                    return Err(Error::Degenerate("a face thinner than the resolution"));
                }
            }
        }
        let mut done: BTreeSet<EdgeId> = BTreeSet::new();
        let mut loops: Vec<(Vec<Fin>, f64)> = Vec::new();
        for (&e, &o) in &used {
            if !done.insert(e) {
                continue;
            }
            let (start, mut to) = ends(e, o);
            let mut cycle = vec![(e, o)];
            while to != start {
                let next = *starts
                    .get(&to)
                    .ok_or(Error::InvalidTopology("an open face loop"))?;
                if !done.insert(next) {
                    return Err(Error::InvalidTopology("a face loop does not close"));
                }
                cycle.push((next, used[&next]));
                to = ends(next, used[&next]).1;
            }
            let mut twice = 0.0;
            let fins: Vec<Fin> = cycle
                .iter()
                .map(|&(e, o)| {
                    let pcurve = plane_pcurve(&parts.edges[e.0].curve, o, frame);
                    let (a, b) = (pcurve.point(0.0), pcurve.point(1.0));
                    twice += a.x * b.y - b.x * a.y;
                    Fin {
                        edge: e,
                        sense: o,
                        pcurve,
                        enclosure: None,
                    }
                })
                .collect();
            loops.push((fins, twice));
        }
        loops.sort_by(|a, b| b.1.total_cmp(&a.1));
        built.push(Built {
            surface: Surface::Plane(frame),
            loops: loops.into_iter().map(|(f, _)| f).collect(),
        });
    }
    // Shells: faces joined by edges; outer (positive volume) or cavities.
    let mut shell_parent: Vec<usize> = (0..groups.len()).collect();
    for u in kept_uses.values() {
        let (g1, g2) = (group_of[u[0].0], group_of[u[1].0]);
        let (r1, r2) = (find(&mut shell_parent, g1), find(&mut shell_parent, g2));
        if r1 != r2 {
            shell_parent[r1.max(r2)] = r1.min(r2);
        }
    }
    let mut shells: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for gi in 0..groups.len() {
        let r = find(&mut shell_parent, gi);
        shells.entry(r).or_default().push(gi);
    }
    let six = R::from_integer(6.into());
    let volume = |faces: &[usize]| -> R {
        let mut v = zero();
        for &gi in faces {
            for &fi in &groups[gi] {
                let p = &frags[fi].pts;
                for i in 1..p.len() - 1 {
                    v += dot(&p[0], &cross(&p[i], &p[i + 1]));
                }
            }
        }
        v / &six
    };
    let (mut outers, mut cavities) = (Vec::new(), Vec::new());
    for faces in shells.into_values() {
        if volume(&faces) > zero() {
            outers.push(faces);
        } else {
            cavities.push(faces);
        }
    }
    // Solids touching at a vertex or along an edge (the edge's four uses
    // are refused above) until the kernel holds non-manifold bodies.
    let mut owner_of: BTreeMap<usize, usize> = BTreeMap::new();
    for (k, faces) in outers.iter().enumerate() {
        for &gi in faces {
            for &fi in &groups[gi] {
                for &v in &cycles[fi] {
                    if owner_of.insert(v, k).is_some_and(|o| o != k) {
                        return Err(Error::Degenerate("two solids touching at a point"));
                    }
                }
            }
        }
    }
    // Each cavity in the outer shell holding one of its points.
    let mut held: Vec<Vec<Vec<usize>>> = vec![Vec::new(); outers.len()];
    for cavity in cavities {
        let p = &frags[groups[cavity[0]][0]].pts[0];
        let owner = outers
            .iter()
            .position(|faces| encloses(&frags, &groups, faces, p).unwrap_or(false))
            .ok_or(Error::InvalidTopology("a cavity outside every shell"))?;
        held[owner].push(cavity);
    }
    // Plans' inputs.
    let tool = |o: Operand| op == Op2::Cut && o == Operand::B;
    let on_segment = |p: &V, a: &V, b: &V| -> bool {
        let ab = sub(b, a);
        let ap = sub(p, a);
        is_zero(&cross(&ab, &ap)) && dot(&ap, &ab) >= zero() && dot(&ap, &ab) <= dot(&ab, &ab)
    };
    let mut out = Vec::new();
    for (index, faces) in outers.iter().enumerate() {
        let mut all: Vec<usize> = faces.clone();
        for c in &held[index] {
            all.extend(c);
        }
        all.sort_unstable();
        let cavity_faces: BTreeSet<usize> = held[index].iter().flatten().copied().collect();
        let face_id: BTreeMap<usize, FaceId> = all
            .iter()
            .enumerate()
            .map(|(i, g)| (*g, FaceId(i)))
            .collect();
        let mut edges: BTreeSet<EdgeId> = BTreeSet::new();
        for &gi in &all {
            for fins in &built[gi].loops {
                edges.extend(fins.iter().map(|f| f.edge));
            }
        }
        let edge_id: BTreeMap<EdgeId, EdgeId> = edges
            .iter()
            .enumerate()
            .map(|(i, e)| (*e, EdgeId(i)))
            .collect();
        let mut verts: BTreeSet<VertexId> = BTreeSet::new();
        for e in &edges {
            verts.extend(parts.edges[e.0].start);
            verts.extend(parts.edges[e.0].end);
        }
        let vertex_id: BTreeMap<VertexId, VertexId> = verts
            .iter()
            .enumerate()
            .map(|(i, v)| (*v, VertexId(i)))
            .collect();
        let mut p = TopologyParts::default();
        for v in &verts {
            p.vertices.push(parts.vertices[v.0].clone());
        }
        for e in &edges {
            let edge = &parts.edges[e.0];
            p.edges.push(Edge {
                start: edge.start.map(|v| vertex_id[&v]),
                end: edge.end.map(|v| vertex_id[&v]),
                curve: edge.curve.clone(),
                fins: Vec::new(),
            });
        }
        let has_cavity = !cavity_faces.is_empty();
        for &gi in &all {
            let inner = cavity_faces.contains(&gi);
            let mut loop_ids = Vec::new();
            for fins in &built[gi].loops {
                let mut fids = Vec::new();
                for fin in fins {
                    let mut fin = fin.clone();
                    fin.edge = edge_id[&fin.edge];
                    let id = FinId(p.fins.len());
                    p.edges[fin.edge.0].fins.push(id);
                    p.fins.push(fin);
                    fids.push(id);
                }
                loop_ids.push(LoopId(p.loops.len()));
                p.loops.push(Loop::Edges {
                    fins: fids,
                    winding: [0, 0],
                });
            }
            p.faces.push(Face {
                surface: built[gi].surface.clone(),
                sense: Orientation::Forward,
                loops: loop_ids,
                front: if inner { ShellId(2) } else { ShellId(0) },
                back: if inner { ShellId(3) } else { ShellId(1) },
                enclosure: None,
            });
        }
        let sides = |inner: bool, side: Side| -> Vec<(FaceId, Side)> {
            all.iter()
                .filter(|g| cavity_faces.contains(g) == inner)
                .map(|g| (face_id[g], side))
                .collect()
        };
        let shell = |region: usize, sides: Vec<(FaceId, Side)>| Shell {
            region: RegionId(region),
            sides,
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        };
        p.shells = vec![
            shell(1, sides(false, Side::Front)),
            shell(0, sides(false, Side::Back)),
        ];
        p.regions = vec![
            Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(1)],
            },
            Region {
                kind: RegionKind::Solid,
                shells: vec![ShellId(0)],
            },
        ];
        if has_cavity {
            p.shells.push(shell(1, sides(true, Side::Front)));
            p.shells.push(shell(2, sides(true, Side::Back)));
            p.regions[1].shells.push(ShellId(2));
            p.regions.push(Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(3)],
            });
        }

        // Plans.
        let tidy = |mut c: Vec<EntityId>, mut t: Vec<EntityId>| {
            c.sort();
            c.dedup();
            t.sort();
            t.dedup();
            t.retain(|x| !c.contains(x));
            (c, t)
        };
        let operand = |id: &EntityId| info.get(id).map(|i| i.0);
        let role_of = |ids: &[EntityId], new: Role| {
            ids.first().and_then(|id| info.get(id)).map_or(new, |i| i.1)
        };
        let mut plans: Vec<SlotPlan> = Vec::new();
        let mut members: BTreeSet<Operand> = BTreeSet::new();
        for &gi in &all {
            let (mut c, mut t) = (Vec::new(), Vec::new());
            for &fi in &groups[gi] {
                let f = &frags[fi];
                let id = face_of(f.src).id;
                if f.kept_way && !tool(f.src.0) {
                    c.push(id);
                } else {
                    t.push(id);
                }
            }
            members.extend(c.iter().filter_map(operand));
            let (c, t) = tidy(c, t);
            let role = role_of(&c, Role::CutFace);
            plans.push((Slot::Face(face_id[&gi]), c, t, EntityKind::Face, role));
        }
        for e in &edges {
            let (mut c, mut t) = (Vec::new(), Vec::new());
            for fe in &chains[e.0] {
                let (a, b) = (&points[fe.0], &points[fe.1]);
                let mut on_input = false;
                for m in &models {
                    for (s, u, id) in &m.edges {
                        if on_segment(a, s, u) && on_segment(b, s, u) {
                            on_input = true;
                            if operand(id).is_some_and(tool) {
                                t.push(*id);
                            } else {
                                c.push(*id);
                            }
                        }
                    }
                }
                if !on_input {
                    for &(fi, _) in &kept_uses[fe] {
                        t.push(face_of(frags[fi].src).id);
                    }
                }
            }
            let (c, t) = tidy(c, t);
            let role = role_of(&c, Role::CutEdge);
            plans.push((Slot::Edge(edge_id[e]), c, t, EntityKind::Edge, role));
        }
        for v in &verts {
            let node = *vertex_of
                .iter()
                .find(|(_, x)| *x == v)
                .map(|(n, _)| n)
                .expect("a vertex's point");
            let point = &points[node];
            let (mut c, mut t) = (Vec::new(), Vec::new());
            for m in &models {
                for (pv, id) in &m.vertices {
                    if pv == point {
                        if operand(id).is_some_and(tool) {
                            t.push(*id);
                        } else {
                            c.push(*id);
                        }
                    }
                }
            }
            if c.is_empty() {
                for m in &models {
                    for (s, u, id) in &m.edges {
                        if on_segment(point, s, u) {
                            t.push(*id);
                        }
                    }
                }
                for e in &incident[&node] {
                    for &(fi, _) in &kept_uses[e] {
                        t.push(face_of(frags[fi].src).id);
                    }
                }
            }
            let (c, t) = tidy(c, t);
            let role = role_of(&c, Role::CutVertex);
            plans.push((Slot::Vertex(vertex_id[v]), c, t, EntityKind::Vertex, role));
        }
        if op == Op2::Cut {
            members = BTreeSet::from([Operand::A]);
        }
        plans.push((
            Slot::Region(RegionId(1)),
            members
                .iter()
                .map(|o| models[usize::from(*o == Operand::B)].region)
                .collect(),
            Vec::new(),
            EntityKind::Region,
            Role::Region,
        ));
        if has_cavity {
            plans.push((
                Slot::Region(RegionId(2)),
                Vec::new(),
                vec![models[1].region],
                EntityKind::Region,
                Role::Region,
            ));
        }
        out.push(Component { parts: p, plans });
    }
    Ok(out)
}

/// Whether a point lies inside a closed shell of fragments: the parity of a
/// ray's crossings, the ray retried in other directions when it meets an
/// edge or lies in a fragment's plane.
fn encloses(frags: &[Frag], groups: &[Vec<usize>], faces: &[usize], p: &V) -> Option<bool> {
    let dirs: [[i64; 3]; 4] = [[7, 3, 5], [2, 11, 13], [17, 5, 3], [3, 19, 7]];
    'dirs: for d in dirs {
        let d: V = d.map(|x| R::from_integer(x.into()));
        let mut count = 0usize;
        for &gi in faces {
            for &fi in &groups[gi] {
                let poly = &frags[fi].pts;
                let n = area2(poly);
                let denom = dot(&n, &d);
                let s = dot(&n, &sub(&poly[0], p));
                if denom == zero() {
                    if s == zero() {
                        continue 'dirs;
                    }
                    continue;
                }
                let t = &s / &denom;
                if t < zero() {
                    continue;
                }
                if t == zero() {
                    continue 'dirs;
                }
                let x = add(p, &scale(&d, &t));
                // Inside the convex cycle: every edge's side the same.
                let mut sign = 0i8;
                for i in 0..poly.len() {
                    let (a, b) = (&poly[i], &poly[(i + 1) % poly.len()]);
                    let side = dot(&cross(&sub(b, a), &sub(&x, a)), &n);
                    let sg = if side > zero() {
                        1
                    } else if side < zero() {
                        -1
                    } else {
                        0
                    };
                    if sg == 0 {
                        continue 'dirs;
                    }
                    if sign == 0 {
                        sign = sg;
                    } else if sg != sign {
                        sign = 2;
                        break;
                    }
                }
                if sign != 2 {
                    count += 1;
                }
            }
        }
        return Some(count % 2 == 1);
    }
    None
}
