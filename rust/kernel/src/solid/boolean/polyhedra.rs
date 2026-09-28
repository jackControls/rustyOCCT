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
use super::stack::{Key, KeyPlan};
use super::{At, What};
use crate::identity::{EntityKind, OperationId, Role};
use crate::profile::boolean::{Op2, Operand};
use crate::solid::split::{q, rational_f64, zero};
use crate::solid::{Construction, MassProperties, Solid};
use crate::topology::{
    plane_pcurve, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId, Orientation,
    Region, RegionId, RegionKind, Shell, ShellId, Side, Slot, Surface, Topology, TopologyParts,
    Vertex, VertexId,
};
use crate::{Error, Frame3, Location, Point2, Point3, Profile, Result, Segment, Tolerance, Vec3};
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

/// One input prism's data, kept to rebuild its result in moved frames.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PrismData {
    pub(super) profile: Profile,
    pub(super) frame: Frame3,
    pub(super) start: f64,
    pub(super) end: f64,
}

impl PrismData {
    pub(super) fn of(solid: &Solid) -> Result<Self> {
        let Construction::Prism(profile) = &solid.construction else {
            return Err(Error::OutOfDomain(
                "a Boolean of solids other than prisms (S9b on)",
            ));
        };
        Ok(Self {
            profile: (**profile).clone(),
            frame: solid.frame,
            start: solid.start,
            end: solid.end,
        })
    }

    /// Where a point lies against the prism, within the resolution.
    fn locate(&self, point: Point3, tolerance: Tolerance) -> Result<Location> {
        let [u, v, w] = self.frame.coordinates(point);
        let tol = tolerance.linear();
        let (lo, hi) = (self.start.min(self.end), self.start.max(self.end));
        use crate::decide::sum_le;
        if !(sum_le(&[lo, -tol], &[w]) && sum_le(&[w], &[hi, tol])) {
            return Ok(Location::Outside);
        }
        let at = self.profile.classify(Point2::new(u, v))?;
        let near = |h: f64| sum_le(&[w, -h], &[tol]) && sum_le(&[h, -w], &[tol]);
        Ok(match at {
            Location::Outside => Location::Outside,
            Location::Boundary => Location::Boundary,
            Location::Inside if near(lo) || near(hi) => Location::Boundary,
            Location::Inside => Location::Inside,
        })
    }
}

/// A polyhedral Boolean's result solid: both inputs, the operation and the
/// solid's index.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Polyhedron {
    pub(super) a: PrismData,
    pub(super) b: PrismData,
    pub(super) op: Op2,
    pub(super) index: usize,
}

/// One result solid: its parts, each slot's provenance and its bounds'
/// height range in the object's frame.
pub(super) struct Component {
    pub(super) parts: TopologyParts,
    pub(super) plans: Vec<KeyPlan>,
}

impl Polyhedron {
    pub(crate) fn tolerance(&self) -> Tolerance {
        self.a.profile.tolerance()
    }

    /// The same solid with both inputs moved.
    pub(crate) fn rebuilt_with(
        &self,
        operation: OperationId,
        motion: Option<crate::RigidTransform>,
        known: Option<MassProperties>,
    ) -> Result<Solid> {
        let mut moved = self.clone();
        if let Some(t) = motion {
            let tol = self.tolerance();
            moved.a.frame = self.a.frame.transformed(t, tol)?;
            moved.b.frame = self.b.frame.transformed(t, tol)?;
        }
        let mut components = build(&moved)?;
        if moved.index >= components.len() {
            return Err(Error::InvalidTopology(
                "a polyhedral Boolean rebuilt differently",
            ));
        }
        let component = components.swap_remove(moved.index);
        moved.solid(component, operation, known)
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
        let tol = self.tolerance();
        let mut parts = topology.to_parts();
        for v in &mut parts.vertices {
            v.position = motion.point(v.position);
        }
        for e in &mut parts.edges {
            let Curve3::LineSegment { start, end } = &mut e.curve else {
                return Err(Error::InvalidTopology(
                    "a polyhedral edge that is not a line",
                ));
            };
            *start = motion.point(*start);
            *end = motion.point(*end);
        }
        for f in &mut parts.faces {
            let Surface::Plane(frame) = &mut f.surface else {
                return Err(Error::InvalidTopology(
                    "a polyhedral face that is not a plane",
                ));
            };
            *frame = frame.transformed(motion, tol)?;
        }
        let mut moved = self.clone();
        moved.a.frame = self.a.frame.transformed(motion, tol)?;
        moved.b.frame = self.b.frame.transformed(motion, tol)?;
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
        let bounds = crate::solid::split::edge_bounds(&topology);
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
            options(self.a.locate(point, tolerance)?),
            options(self.b.locate(point, tolerance)?),
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

/// A boundary face of a model: its plane (outward), convex pieces (outward
/// cycles) and key.
struct MFace {
    plane: Plane,
    pieces: Vec<Vec<V>>,
    key: Key,
}

/// An input's exact model.
struct Model {
    faces: Vec<MFace>,
    /// Convex cells with disjoint interiors, as planes (inside `<= 0`).
    cells: Vec<Vec<Plane>>,
    edges: Vec<(V, V, Key)>,
    vertices: Vec<(V, Key)>,
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

fn model(data: &PrismData, op: Operand) -> Result<Model> {
    let f = data.frame;
    let (o, x, y, n) = (
        vq(f.origin() - Point3::ORIGIN),
        vq(f.x()),
        vq(f.y()),
        vq(f.normal()),
    );
    let (lo, hi) = (q(data.start.min(data.end)), q(data.start.max(data.end)));
    let point =
        |u: &R, v: &R, w: &R| add(&add(&o, &scale(&x, u)), &add(&scale(&y, v), &scale(&n, w)));
    let lift = |p: &P2, w: &R| point(&p.0, &p.1, w);
    let mut polys: Vec<Vec<P2>> = Vec::new();
    for boundary in data.profile.boundaries() {
        let pts: Vec<Point2> = match (boundary.polygon_vertices(), boundary.path_geometry()) {
            (Some(p), _) => p.to_vec(),
            (None, Some((p, segs))) if segs.iter().all(|s| matches!(s, Segment::Line)) => {
                p.to_vec()
            }
            _ => {
                return Err(Error::OutOfDomain(
                    "a Boolean of prisms with arcs in frames with different axes (S9c)",
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
            key: (op, At::Low, What::Cap, 0, 0),
        },
        MFace {
            plane: Plane::through(up.clone(), &high_at),
            pieces: traps
                .iter()
                .map(|t| t.iter().map(|p| lift(p, &hi)).collect())
                .collect(),
            key: (op, At::High, What::Cap, 0, 0),
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
                key: (op, At::Swept, What::Wall, b, j),
            });
            edges.push((a.clone(), c.clone(), (op, At::Low, What::Edge, b, j)));
            edges.push((a2.clone(), c2.clone(), (op, At::High, What::Edge, b, j)));
            edges.push((a.clone(), a2.clone(), (op, At::Swept, What::Vertical, b, j)));
            vertices.push((a, (op, At::Low, What::Vertex, b, j)));
            vertices.push((a2, (op, At::High, What::Vertex, b, j)));
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
        cells,
        edges,
        vertices,
    })
}

/// Whether `p + e d` lies in the model for every small `e > 0` (closed
/// cells: `p + e d` is never on the model's boundary here).
fn inside(model: &Model, p: &V, d: &V) -> bool {
    model.cells.iter().any(|cell| {
        cell.iter().all(|plane| {
            let s = plane.eval(p);
            s < zero() || (s == zero() && dot(&plane.n, d) <= zero())
        })
    })
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
    let tolerance = poly.tolerance();
    let models = [model(&poly.a, Operand::A)?, model(&poly.b, Operand::B)?];
    let op = poly.op;
    // Each input's faces' planes, unoriented and distinct.
    let planes: Vec<Vec<Plane>> = models
        .iter()
        .map(|m| {
            let mut seen = BTreeSet::new();
            let mut out = Vec::new();
            for f in &m.faces {
                let (n, d) = f.plane.canonical(false);
                if seen.insert((n.clone(), d.clone())) {
                    out.push(Plane { n, d });
                }
            }
            out
        })
        .collect();
    // Fragments, classified and kept.
    let mut frags: Vec<Frag> = Vec::new();
    for (k, m) in models.iter().enumerate() {
        let other = &models[1 - k];
        let operand = if k == 0 { Operand::A } else { Operand::B };
        for (fi, face) in m.faces.iter().enumerate() {
            let n = &face.plane.n;
            let mut parts: Vec<Vec<V>> = face.pieces.clone();
            for plane in &planes[1 - k] {
                let mut next = Vec::new();
                for part in parts {
                    let s: Vec<R> = part.iter().map(|p| plane.eval(p)).collect();
                    let crosses = s.iter().any(|x| *x < zero()) && s.iter().any(|x| *x > zero());
                    if crosses {
                        for keep in [-1i8, 1] {
                            let g = clip(&part, plane, keep);
                            if g.len() >= 3 && !is_zero(&area2(&g)) {
                                next.push(g);
                            }
                        }
                    } else {
                        next.push(part);
                    }
                }
                parts = next;
            }
            for part in parts {
                let count = R::from_integer((part.len() as i64).into());
                let c = scale(
                    &part
                        .iter()
                        .fold([zero(), zero(), zero()], |acc, p| add(&acc, p)),
                    &(R::from_integer(1.into()) / count),
                );
                let (front, back) = (inside(other, &c, n), inside(other, &c, &neg(n)));
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
    let approx: Vec<[f64; 3]> = points.iter().map(|p| rounded(p).to_array()).collect();
    let tol = tolerance.linear();
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for f in &frags {
        let mut cycle = Vec::new();
        for i in 0..f.pts.len() {
            let (a, b) = (&f.pts[i], &f.pts[(i + 1) % f.pts.len()]);
            let (ia, ib) = (ids[a], ids[b]);
            cycle.push(ia);
            let ab = sub(b, a);
            let len2 = dot(&ab, &ab);
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
                let av = sub(&points[iv], a);
                if !is_zero(&cross(&ab, &av)) {
                    continue;
                }
                let t = dot(&av, &ab);
                if t > zero() && t < len2 {
                    on.push((t, iv));
                }
            }
            on.sort();
            cycle.extend(on.into_iter().map(|(_, iv)| iv));
        }
        cycles.push(cycle);
    }
    // Faces: fragments of one oriented plane joined across shared edges.
    let frag_plane: Vec<(V, R)> = frags
        .iter()
        .map(|f| {
            let n = area2(&f.pts);
            Plane::through(n, &f.pts[0]).canonical(true)
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
        let f0 = &frags[group[0]];
        let normal = area2(&f0.pts);
        let nf = Vec3::new(
            rational_f64(&normal[0]),
            rational_f64(&normal[1]),
            rational_f64(&normal[2]),
        );
        let origin = rounded(&f0.pts[0]);
        let hint = rounded(&f0.pts[1]) - origin;
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
        for &v in &corners {
            let p = &points[v];
            for &w in &corners {
                if w > v {
                    let d = sub(&points[w], p);
                    if dot(&d, &d) <= tol2 {
                        return Err(Error::Degenerate("a face thinner than the resolution"));
                    }
                }
            }
            for e in used.keys() {
                let (a, b) = edge_ends[e.0];
                if a == v || b == v {
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
        let tidy = |mut c: Vec<Key>, mut t: Vec<Key>| {
            c.sort();
            c.dedup();
            t.sort();
            t.dedup();
            t.retain(|x| !c.contains(x));
            (c, t)
        };
        let face_key =
            |(o, fi): (Operand, usize)| models[usize::from(o == Operand::B)].faces[fi].key;
        let mut plans: Vec<KeyPlan> = Vec::new();
        let mut members: BTreeSet<Operand> = BTreeSet::new();
        for &gi in &all {
            let (mut c, mut t) = (Vec::new(), Vec::new());
            for &fi in &groups[gi] {
                let f = &frags[fi];
                let key = face_key(f.src);
                if f.kept_way && !tool(f.src.0) {
                    c.push(key);
                } else {
                    t.push(key);
                }
            }
            members.extend(c.iter().map(|k| k.0));
            let role = match c.first() {
                Some((_, At::Low, What::Cap, ..)) => Role::StartCap,
                Some((_, At::High, What::Cap, ..)) => Role::EndCap,
                Some(_) => Role::Wall,
                None => Role::CutFace,
            };
            let (c, t) = tidy(c, t);
            plans.push((Slot::Face(face_id[&gi]), c, t, EntityKind::Face, role));
        }
        for e in &edges {
            let (mut c, mut t) = (Vec::new(), Vec::new());
            for fe in &chains[e.0] {
                let (a, b) = (&points[fe.0], &points[fe.1]);
                let mut on_input = false;
                for m in &models {
                    for (s, u, key) in &m.edges {
                        if on_segment(a, s, u) && on_segment(b, s, u) {
                            on_input = true;
                            if tool(key.0) {
                                t.push(*key);
                            } else {
                                c.push(*key);
                            }
                        }
                    }
                }
                if !on_input {
                    for &(fi, _) in &kept_uses[fe] {
                        t.push(face_key(frags[fi].src));
                    }
                }
            }
            let role = match c.first() {
                Some((_, At::Low, ..)) => Role::BottomEdge,
                Some((_, At::High, ..)) => Role::TopEdge,
                Some(_) => Role::Vertical,
                None => Role::CutEdge,
            };
            let (c, t) = tidy(c, t);
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
                for (pv, key) in &m.vertices {
                    if pv == point {
                        if tool(key.0) {
                            t.push(*key);
                        } else {
                            c.push(*key);
                        }
                    }
                }
            }
            if c.is_empty() {
                for m in &models {
                    for (s, u, key) in &m.edges {
                        if on_segment(point, s, u) {
                            t.push(*key);
                        }
                    }
                }
                for e in &incident[&node] {
                    for &(fi, _) in &kept_uses[e] {
                        t.push(face_key(frags[fi].src));
                    }
                }
            }
            let role = match c.first() {
                Some((_, At::Low, ..)) => Role::BottomVertex,
                Some((_, At::High, ..)) => Role::TopVertex,
                _ => Role::CutVertex,
            };
            let (c, t) = tidy(c, t);
            plans.push((Slot::Vertex(vertex_id[v]), c, t, EntityKind::Vertex, role));
        }
        if op == Op2::Cut {
            members = BTreeSet::from([Operand::A]);
        }
        plans.push((
            Slot::Region(RegionId(1)),
            members
                .iter()
                .map(|o| (*o, At::Swept, What::Region, 0, 0))
                .collect(),
            Vec::new(),
            EntityKind::Region,
            Role::Region,
        ));
        if has_cavity {
            plans.push((
                Slot::Region(RegionId(2)),
                Vec::new(),
                vec![(Operand::B, At::Swept, What::Region, 0, 0)],
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
