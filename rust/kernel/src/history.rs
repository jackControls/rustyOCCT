//! Complete operation histories, their composition and an independent checker
//! (see `IDENTITY_AND_HISTORY.md`, contracts 2 and 3).
//!
//! A [`History`] relates every entity of an operation's input bodies to its
//! outputs. [`check`] validates a history against the input and output
//! entities alone; it knows nothing about the algorithm that produced it.
//! [`History::then`] composes histories and [`History::resolve`] answers
//! "what became of this id" without ever choosing a survivor.
use crate::attributes::AttributeOutcome;
use crate::identity::{
    encode_parents, role_code, AlgorithmLevel, EntityId, EntityKind, OperationId, OperationKind,
    Parent, Role,
};
use crate::topology::{Curve3, Orientation, RegionKind, Surface};
use crate::{Point3, Tolerance};
use num_rational::BigRational as R;
use std::collections::{BTreeMap, BTreeSet};

/// How one operation related input entities to output entities.
///
/// `Modified` normally keeps the id (`from == to`). A composed history may
/// relate an input to a different output one-to-one, e.g. a split whose
/// pieces are merged again; that is still `Modified`, never `Unchanged`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Relation {
    Unchanged {
        id: EntityId,
    },
    Modified {
        from: EntityId,
        to: EntityId,
    },
    Generated {
        from: Vec<Parent>,
        to: EntityId,
        role: Role,
    },
    /// Children in canonical ordinal order; no child is a primary survivor.
    Split {
        from: EntityId,
        into: Vec<EntityId>,
    },
    /// Parents in canonical order; no parent is a primary survivor.
    Merged {
        from: Vec<EntityId>,
        into: EntityId,
    },
    Deleted {
        id: EntityId,
    },
}

impl Relation {
    /// Input entity ids this relation accounts for (not Generated parents).
    pub fn sources(&self) -> Vec<EntityId> {
        match self {
            Self::Unchanged { id } | Self::Deleted { id } => vec![*id],
            Self::Modified { from, .. } | Self::Split { from, .. } => vec![*from],
            Self::Merged { from, .. } => from.clone(),
            Self::Generated { .. } => Vec::new(),
        }
    }
    /// Output entity ids this relation produces.
    pub fn targets(&self) -> Vec<EntityId> {
        match self {
            Self::Unchanged { id } => vec![*id],
            Self::Modified { to, .. } | Self::Generated { to, .. } => vec![*to],
            Self::Split { into, .. } => into.clone(),
            Self::Merged { into, .. } => vec![*into],
            Self::Deleted { .. } => Vec::new(),
        }
    }
    fn code(&self) -> u8 {
        match self {
            Self::Unchanged { .. } => 1,
            Self::Modified { .. } => 2,
            Self::Generated { .. } => 3,
            Self::Split { .. } => 4,
            Self::Merged { .. } => 5,
            Self::Deleted { .. } => 6,
        }
    }
    /// Canonical order: encoded sources, relation kind, then targets (and the
    /// role of a generated target). Mirrored by `identity_reference.py`.
    pub fn sort_key(&self) -> Vec<u8> {
        let sources: Vec<Parent> = match self {
            Self::Generated { from, .. } => from.clone(),
            _ => self.sources().into_iter().map(Parent::Entity).collect(),
        };
        let mut key = Vec::new();
        encode_parents(&mut key, &sources);
        key.push(self.code());
        let targets = self.targets();
        key.extend_from_slice(&(targets.len() as u32).to_le_bytes());
        for t in targets {
            key.extend_from_slice(&t.0);
        }
        if let Self::Generated { role, .. } = self {
            key.push(role_code(*role));
        }
        key
    }
}

/// One operation a history covers, with the behaviour version it ran at
/// (H8).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Step {
    pub operation: OperationId,
    pub kind: OperationKind,
    pub level: AlgorithmLevel,
}

/// The complete relation list of one operation (or of a composed chain).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct History {
    pub operation: OperationId,
    pub kind: OperationKind,
    /// Every operation the history covers, oldest first: one step for a
    /// single operation, and the concatenation of both lists for
    /// [`History::then`]. A composite's level is this list (H8); composites
    /// are never replayed.
    pub steps: Vec<Step>,
    pub input_bodies: Vec<EntityId>,
    pub output_bodies: Vec<EntityId>,
    /// Sorted by [`Relation::sort_key`].
    pub relations: Vec<Relation>,
    /// Attribute outcomes (contract 4); not composed by [`History::then`].
    pub attributes: Vec<AttributeOutcome>,
}

/// What became of one input id. The kernel never picks a survivor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Resolution {
    /// Unchanged or modified.
    Same(EntityId),
    /// In canonical ordinal order.
    Split(Vec<EntityId>),
    Merged {
        into: EntityId,
        with: Vec<EntityId>,
    },
    Deleted,
    /// The id was not an input of this history: a caller error.
    Unknown,
}

impl History {
    /// A history with its relations put in canonical order.
    pub fn new(
        operation: OperationId,
        kind: OperationKind,
        input_bodies: Vec<EntityId>,
        output_bodies: Vec<EntityId>,
        mut relations: Vec<Relation>,
        attributes: Vec<AttributeOutcome>,
    ) -> Self {
        relations.sort_by_cached_key(Relation::sort_key);
        Self {
            operation,
            kind,
            steps: vec![Step {
                operation,
                kind,
                level: AlgorithmLevel::FIRST,
            }],
            input_bodies,
            output_bodies,
            relations,
            attributes,
        }
    }

    /// A single operation's level; `None` for a composite, whose level is
    /// its step list.
    pub fn level(&self) -> Option<AlgorithmLevel> {
        match self.steps.as_slice() {
            [step] if self.kind != OperationKind::Composite => Some(step.level),
            _ => None,
        }
    }

    /// The same single-operation history recorded at `level`.
    pub fn at_level(mut self, level: AlgorithmLevel) -> Self {
        for step in &mut self.steps {
            step.level = level;
        }
        self
    }

    pub fn resolve(&self, id: EntityId) -> Resolution {
        for relation in &self.relations {
            match relation {
                Relation::Unchanged { id: x } if *x == id => return Resolution::Same(id),
                Relation::Modified { from, to } if *from == id => return Resolution::Same(*to),
                Relation::Split { from, into } if *from == id => {
                    return Resolution::Split(into.clone())
                }
                Relation::Merged { from, into } if from.contains(&id) => {
                    return Resolution::Merged {
                        into: *into,
                        with: from.iter().copied().filter(|f| *f != id).collect(),
                    }
                }
                Relation::Deleted { id: x } if *x == id => return Resolution::Deleted,
                _ => {}
            }
        }
        Resolution::Unknown
    }

    /// This history followed by `next`, from this history's inputs to
    /// `next`'s outputs. Fails when the histories do not chain, or when the
    /// chain relates entities many-to-many, which no single relation can say.
    pub fn then(&self, next: &History) -> std::result::Result<History, HistoryIssue> {
        let invalid = |id| HistoryIssue {
            kind: HistoryIssueKind::CompositionInvalid,
            id,
        };
        if self.output_bodies != next.input_bodies {
            let body = self
                .output_bodies
                .iter()
                .chain(&next.input_bodies)
                .find(|b| !self.output_bodies.contains(b) || !next.input_bodies.contains(b))
                .or(self.output_bodies.first())
                .or(next.input_bodies.first())
                .copied()
                .unwrap_or(EntityId([0; 16]));
            return Err(invalid(body));
        }
        let a = Flows::of(self);
        let b = Flows::of(next);
        for t in &a.targets {
            if !b.sources.contains(t) {
                return Err(invalid(*t));
            }
        }
        for s in &b.sources {
            if !a.targets.contains(s) {
                return Err(invalid(*s));
            }
        }
        // Where each original input finally goes, in canonical child order.
        let mut finals: BTreeMap<EntityId, Vec<EntityId>> = BTreeMap::new();
        let mut contributors: BTreeMap<EntityId, BTreeSet<EntityId>> = BTreeMap::new();
        for s in &a.sources {
            let mut out: Vec<EntityId> = Vec::new();
            for t in a.flow.get(s).into_iter().flatten() {
                for x in b.flow.get(t).into_iter().flatten() {
                    if !out.contains(x) {
                        out.push(*x);
                    }
                }
            }
            for x in &out {
                contributors.entry(*x).or_default().insert(*s);
            }
            finals.insert(*s, out);
        }
        let mut relations = Vec::new();
        let mut merged_done = BTreeSet::new();
        for (s, out) in &finals {
            match out.as_slice() {
                [] => relations.push(Relation::Deleted { id: *s }),
                [x] => {
                    let from = &contributors[x];
                    if from.len() == 1 {
                        if a.unchanged.contains(s) && b.unchanged.contains(s) && x == s {
                            relations.push(Relation::Unchanged { id: *s });
                        } else {
                            relations.push(Relation::Modified { from: *s, to: *x });
                        }
                    } else if from.iter().all(|c| finals[c] == [*x]) {
                        if merged_done.insert(*x) {
                            relations.push(Relation::Merged {
                                from: from.iter().copied().collect(),
                                into: *x,
                            });
                        }
                    } else {
                        return Err(invalid(*s));
                    }
                }
                many => {
                    if many.iter().any(|x| contributors[x].len() != 1) {
                        return Err(invalid(*s));
                    }
                    relations.push(Relation::Split {
                        from: *s,
                        into: many.to_vec(),
                    });
                }
            }
        }
        for (parents, t, role) in &a.generated {
            for x in b.flow.get(t).into_iter().flatten() {
                relations.push(Relation::Generated {
                    from: parents.clone(),
                    to: *x,
                    role: *role,
                });
            }
        }
        for (parents, x, role) in &b.generated {
            let mut from: Vec<Parent> = Vec::new();
            for p in parents {
                let mapped = match p {
                    Parent::Entity(t) if a.origin.contains_key(t) => a.origin[t].clone(),
                    other => vec![*other],
                };
                for q in mapped {
                    if !from.contains(&q) {
                        from.push(q);
                    }
                }
            }
            relations.push(Relation::Generated {
                from,
                to: *x,
                role: *role,
            });
        }
        let mut composed = History::new(
            next.operation,
            OperationKind::Composite,
            self.input_bodies.clone(),
            next.output_bodies.clone(),
            relations,
            Vec::new(),
        );
        composed.steps = self.steps.iter().chain(&next.steps).copied().collect();
        Ok(composed)
    }
}

/// One history as flows from sources to targets.
struct Flows {
    sources: BTreeSet<EntityId>,
    targets: BTreeSet<EntityId>,
    flow: BTreeMap<EntityId, Vec<EntityId>>,
    unchanged: BTreeSet<EntityId>,
    generated: Vec<(Vec<Parent>, EntityId, Role)>,
    /// Each target expressed through this history's own inputs.
    origin: BTreeMap<EntityId, Vec<Parent>>,
}

impl Flows {
    fn of(h: &History) -> Self {
        let mut f = Flows {
            sources: BTreeSet::new(),
            targets: BTreeSet::new(),
            flow: BTreeMap::new(),
            unchanged: BTreeSet::new(),
            generated: Vec::new(),
            origin: BTreeMap::new(),
        };
        for r in &h.relations {
            f.sources.extend(r.sources());
            f.targets.extend(r.targets());
            match r {
                Relation::Generated { from, to, role } => {
                    f.generated.push((from.clone(), *to, *role));
                    let o = f.origin.entry(*to).or_default();
                    for p in from {
                        if !o.contains(p) {
                            o.push(*p);
                        }
                    }
                }
                Relation::Deleted { .. } => {}
                _ => {
                    if let Relation::Unchanged { id } = r {
                        f.unchanged.insert(*id);
                    }
                    for s in r.sources() {
                        f.flow.entry(s).or_default().extend(r.targets());
                        for t in r.targets() {
                            let o = f.origin.entry(t).or_default();
                            if !o.contains(&Parent::Entity(s)) {
                                o.push(Parent::Entity(s));
                            }
                        }
                    }
                }
            }
        }
        f
    }
}

// ------------------------------------------------------------------ checker

/// The geometry the checker compares.
#[derive(Debug, Clone, PartialEq)]
pub enum Geometry {
    Point(Point3),
    Curve(Curve3),
    Surface {
        surface: Surface,
        orientation: Orientation,
    },
    /// A bounded region; its structure lists its shells' face sides.
    Region(RegionKind),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Family {
    Point,
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    Plane,
    Cylinder,
    Region,
}

impl Geometry {
    pub fn family(&self) -> Family {
        match self {
            Self::Point(_) => Family::Point,
            Self::Curve(Curve3::LineSegment { .. }) => Family::Line,
            Self::Curve(Curve3::EllipseArc { .. }) => Family::Ellipse,
            Self::Curve(Curve3::HyperbolaArc { .. }) => Family::Hyperbola,
            Self::Curve(Curve3::ParabolaArc { .. }) => Family::Parabola,
            Self::Curve(_) => Family::Circle,
            Self::Surface {
                surface: Surface::Plane(_),
                ..
            } => Family::Plane,
            Self::Surface { .. } => Family::Cylinder,
            Self::Region(_) => Family::Region,
        }
    }
}

/// One entity as the checker sees it.
#[derive(Debug, Clone, PartialEq)]
pub struct EntityInfo {
    pub kind: EntityKind,
    /// The ordinal of the entity's derivation.
    pub ordinal: u32,
    pub geometry: Geometry,
    /// Face loops as (edge id, use orientation); an edge's start and end
    /// vertex ids as one two-element loop; empty for a vertex.
    pub structure: Vec<Vec<(EntityId, Orientation)>>,
}

/// Every entity of one body, keyed by id.
#[derive(Debug, Clone, PartialEq)]
pub struct EntitySet {
    pub body: EntityId,
    pub tolerance: Tolerance,
    pub entities: BTreeMap<EntityId, EntityInfo>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum HistoryIssueKind {
    MissingSource,
    DuplicateSource,
    MissingTarget,
    DanglingId,
    DimensionMismatch,
    UnchangedGeometryDiffers,
    ModifiedFamilyDiffers,
    ModifiedIdChanged,
    SplitSupportDiffers,
    MergedSupportDiffers,
    DeletedStillPresent,
    OrdinalNotCanonical,
    CompositionInvalid,
    InvalidArity,
    BodyMismatch,
    DuplicateId,
    RelationOrder,
    StepsInvalid,
}

impl HistoryIssueKind {
    pub fn name(self) -> &'static str {
        use HistoryIssueKind::*;
        match self {
            MissingSource => "missing_source",
            DuplicateSource => "duplicate_source",
            MissingTarget => "missing_target",
            DanglingId => "dangling_id",
            DimensionMismatch => "dimension_mismatch",
            UnchangedGeometryDiffers => "unchanged_geometry_differs",
            ModifiedFamilyDiffers => "modified_family_differs",
            ModifiedIdChanged => "modified_id_changed",
            SplitSupportDiffers => "split_support_differs",
            MergedSupportDiffers => "merged_support_differs",
            DeletedStillPresent => "deleted_still_present",
            OrdinalNotCanonical => "ordinal_not_canonical",
            CompositionInvalid => "composition_invalid",
            InvalidArity => "invalid_arity",
            BodyMismatch => "body_mismatch",
            DuplicateId => "duplicate_id",
            RelationOrder => "relation_order",
            StepsInvalid => "steps_invalid",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HistoryIssue {
    pub kind: HistoryIssueKind,
    /// The offending entity (or body) id.
    pub id: EntityId,
}

impl std::fmt::Display for HistoryIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.kind.name(), self.id)
    }
}

fn r(x: f64) -> Option<R> {
    R::from_float(x)
}

/// Exact |p - q|^2 <= tol^2.
fn near(p: Point3, q: Point3, tol: f64) -> bool {
    let (Some(t), Some(d)) = (
        r(tol),
        (0..3)
            .map(|i| {
                let (a, b) = (r(p.to_array()[i])?, r(q.to_array()[i])?);
                Some((&a - &b) * (&a - &b))
            })
            .sum::<Option<R>>(),
    ) else {
        return false;
    };
    d <= &t * &t
}

/// Exact: point p lies within tol of the infinite line through a and b.
fn on_line(p: Point3, a: Point3, b: Point3, tol: f64) -> bool {
    let conv = |v: Point3| v.to_array().map(r);
    let (pa, aa, ba) = (conv(p), conv(a), conv(b));
    if [&pa, &aa, &ba]
        .iter()
        .any(|v| v.iter().any(Option::is_none))
    {
        return false;
    }
    let get = |v: &[Option<R>; 3], i: usize| v[i].clone().unwrap();
    let d: [R; 3] = std::array::from_fn(|i| get(&ba, i) - get(&aa, i));
    let w: [R; 3] = std::array::from_fn(|i| get(&pa, i) - get(&aa, i));
    let cross = [
        &d[1] * &w[2] - &d[2] * &w[1],
        &d[2] * &w[0] - &d[0] * &w[2],
        &d[0] * &w[1] - &d[1] * &w[0],
    ];
    let c2: R = cross.iter().map(|c| c * c).sum();
    let d2: R = d.iter().map(|c| c * c).sum();
    let t = r(tol).unwrap();
    d2 > R::from_integer(0.into()) && c2 <= &t * &t * d2
}

fn circle_of(c: &Curve3) -> Option<(crate::Frame3, f64)> {
    match c {
        Curve3::Circle { frame, radius } | Curve3::CircularArc { frame, radius, .. } => {
            Some((*frame, *radius))
        }
        Curve3::LineSegment { .. }
        | Curve3::BSpline(_)
        | Curve3::EllipseArc { .. }
        | Curve3::HyperbolaArc { .. }
        | Curve3::ParabolaArc { .. }
        | Curve3::Section(_)
        | Curve3::Meet(_)
        | Curve3::Rise(_)
        | Curve3::Toric(_) => None,
    }
}

/// An ellipse arc's frame and semi-axes.
fn ellipse_of(c: &Curve3) -> Option<(crate::Frame3, f64, f64)> {
    match c {
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            ..
        } => Some((*frame, *major, *minor)),
        _ => None,
    }
}

/// Whether the centre and axes' ends of ellipse `a` lie on ellipse `b`
/// within `tol` (the centre at `b`'s, the ends at most `tol` off `b`'s
/// plane and off its curve, measured along its radius from the centre).
fn on_ellipse(a: &(crate::Frame3, f64, f64), b: &(crate::Frame3, f64, f64), tol: f64) -> bool {
    let (fa, ma, na) = *a;
    let (fb, mb, nb) = *b;
    if !near(fa.origin(), fb.origin(), tol) || mb <= 0.0 || nb <= 0.0 {
        return false;
    }
    let (o, x, y) = (fa.origin(), fa.x(), fa.y());
    [o + x * ma, o + x * -ma, o + y * na, o + y * -na]
        .iter()
        .all(|p| {
            let [u, v, w] = fb.coordinates(*p);
            let rho = (u * u + v * v).sqrt();
            // The ellipse's radius in the point's direction.
            let k = ((u / mb).powi(2) + (v / nb).powi(2)).sqrt();
            w.abs() <= tol && k > 0.0 && (rho - rho / k).abs() <= tol
        })
}

fn exact(v: [f64; 3]) -> Option<[R; 3]> {
    let [x, y, z] = v.map(r);
    Some([x?, y?, z?])
}

fn dot(a: &[R; 3], b: &[R; 3]) -> R {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

fn cross(a: &[R; 3], b: &[R; 3]) -> [R; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

fn minus(a: &[R; 3], b: &[R; 3]) -> [R; 3] {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}

/// Exact: point p lies within tol of the line through o with direction n.
fn on_axis(p: Point3, o: &[R; 3], n: &[R; 3], tol: &R) -> Option<bool> {
    let c = cross(&minus(&exact(p.to_array())?, o), n);
    Some(dot(&c, &c) <= tol * tol * dot(n, n))
}

/// A face `piece` (boundary edges looked up in `edges`) lies on the surface
/// of `whole`, facing the same way. A plane piece faces the whole's normal,
/// and its boundary projected onto its own plane lies within tol of the
/// whole plane: line ends, and circles with normals exactly parallel to
/// both planes by their centres (the distance to the whole plane is affine
/// on the piece's plane, so the projected boundary bounds the face). A
/// cylinder piece has an exactly parallel axis whose distance from the
/// whole's axis plus the radius difference is within tol, so the two surfaces
/// are within tol everywhere. Anything else is not the same support.
fn same_surface(piece: &EntityInfo, whole: &Geometry, tol: f64, edges: &Entities) -> Option<bool> {
    let (
        Geometry::Surface {
            surface: a,
            orientation: oa,
        },
        Geometry::Surface {
            surface: b,
            orientation: ob,
        },
    ) = (&piece.geometry, whole)
    else {
        return Some(false);
    };
    // Planes compare as oriented planes (a normal and a sense are one
    // representation of an oriented plane among two); other surfaces by
    // their senses.
    let planes = matches!((a, b), (Surface::Plane(_), Surface::Plane(_)));
    let splines = matches!((a, b), (Surface::BSpline(_), Surface::BSpline(_)));
    if oa != ob && !planes && !splines {
        return Some(false);
    }
    let tol_f = tol;
    let tol = r(tol)?;
    let zero = R::from_integer(0.into());
    match (a, b) {
        (Surface::Plane(fa), Surface::Plane(fb)) => {
            let oriented = |f: &crate::Frame3, o: &crate::topology::Orientation| {
                let n = f.normal();
                match o {
                    crate::topology::Orientation::Forward => n,
                    crate::topology::Orientation::Reversed => -n,
                }
            };
            let (na, nb) = (
                exact(oriented(fa, oa).to_array())?,
                exact(oriented(fb, ob).to_array())?,
            );
            let (oa, ob) = (
                exact(fa.origin().to_array())?,
                exact(fb.origin().to_array())?,
            );
            if dot(&na, &nb) <= zero || piece.structure.is_empty() {
                return Some(false);
            }
            // A boundary point projected onto the piece's own plane.
            let on_piece = |p: Point3| -> Option<[R; 3]> {
                let p = exact(p.to_array())?;
                let k = dot(&minus(&p, &oa), &na) / dot(&na, &na);
                Some([
                    &p[0] - &k * &na[0],
                    &p[1] - &k * &na[1],
                    &p[2] - &k * &na[2],
                ])
            };
            let within = |q: [R; 3]| {
                let d = dot(&minus(&q, &ob), &nb);
                &d * &d <= &tol * &tol * dot(&nb, &nb)
            };
            for (edge, _) in piece.structure.iter().flatten() {
                let ok = match &edges.get(edge)?.0.geometry {
                    Geometry::Curve(Curve3::LineSegment { start, end }) => {
                        within(on_piece(*start)?) && within(on_piece(*end)?)
                    }
                    // A nonrational spline lies in its poles' hull, and the
                    // projection and the distance are affine (S8b).
                    Geometry::Curve(Curve3::BSpline(span)) => {
                        let curve = span.curve();
                        !curve.is_rational()
                            && curve
                                .poles()
                                .iter()
                                .map(|p| Some(within(on_piece(*p)?)))
                                .collect::<Option<Vec<bool>>>()?
                                .into_iter()
                                .all(|x| x)
                    }
                    // An ellipse (a section's edge, its frame rounded):
                    // its centre and axes' ends within both planes' reach.
                    Geometry::Curve(Curve3::EllipseArc {
                        frame,
                        major,
                        minor,
                        ..
                    }) => {
                        let (o, x, y) = (frame.origin(), frame.x(), frame.y());
                        [
                            o,
                            o + x * *major,
                            o + x * -*major,
                            o + y * *minor,
                            o + y * -*minor,
                        ]
                        .iter()
                        .map(|p| Some(within(on_piece(*p)?)))
                        .collect::<Option<Vec<bool>>>()?
                        .into_iter()
                        .all(|x| x)
                    }
                    // A torus's spiric section (S9d.4a): points along it
                    // within both planes' reach.
                    Geometry::Curve(c @ Curve3::Section(_)) => [0.0, 0.25, 0.5, 0.75, 1.0]
                        .iter()
                        .map(|&t| Some(within(on_piece(c.point(t))?)))
                        .collect::<Option<Vec<bool>>>()?
                        .into_iter()
                        .all(|x| x),
                    // A hyperbola or parabola (a cone's section, S9d.3a):
                    // its frame's plane, points at its own reach, within
                    // both planes' reach.
                    Geometry::Curve(Curve3::HyperbolaArc {
                        frame,
                        major: reach,
                        ..
                    })
                    | Geometry::Curve(Curve3::ParabolaArc {
                        frame,
                        focal: reach,
                        ..
                    }) => {
                        let (o, x, y) = (frame.origin(), frame.x(), frame.y());
                        [o, o + x * *reach, o + y * *reach]
                            .iter()
                            .map(|p| Some(within(on_piece(*p)?)))
                            .collect::<Option<Vec<bool>>>()?
                            .into_iter()
                            .all(|x| x)
                    }
                    // A circle: exactly parallel to both planes, its
                    // centre's distance; else (a section's circle, its
                    // frame rounded: S9d.1) its centre and axes' ends
                    // within both planes' reach, as an ellipse's.
                    Geometry::Curve(c) => {
                        let (frame, radius) = circle_of(c)?;
                        let m = exact(frame.normal().to_array())?;
                        if cross(&m, &na).iter().all(|x| *x == zero)
                            && cross(&m, &nb).iter().all(|x| *x == zero)
                        {
                            within(on_piece(frame.origin())?)
                        } else {
                            let (o, x, y) = (frame.origin(), frame.x(), frame.y());
                            [
                                o,
                                o + x * radius,
                                o + x * -radius,
                                o + y * radius,
                                o + y * -radius,
                            ]
                            .iter()
                            .map(|p| Some(within(on_piece(*p)?)))
                            .collect::<Option<Vec<bool>>>()?
                            .into_iter()
                            .all(|x| x)
                        }
                    }
                    _ => false,
                };
                if !ok {
                    return Some(false);
                }
            }
            Some(true)
        }
        (
            Surface::Cylinder {
                frame: fa,
                radius: ra,
            },
            Surface::Cylinder {
                frame: fb,
                radius: rb,
            },
        ) => {
            let (na, nb) = (
                exact(fa.normal().to_array())?,
                exact(fb.normal().to_array())?,
            );
            // Axis distance plus radius difference within tol keeps the
            // surfaces within tol of each other everywhere.
            let dr = r(*ra)? - r(*rb)?;
            let slack = &tol - if dr < zero { -dr } else { dr };
            Some(
                slack >= zero
                    && cross(&na, &nb).iter().all(|x| *x == zero)
                    && on_axis(fa.origin(), &exact(fb.origin().to_array())?, &nb, &slack)?,
            )
        }
        // A sphere: centres within the tolerance less the radius difference.
        (
            Surface::Sphere {
                frame: fa,
                radius: ra,
            },
            Surface::Sphere {
                frame: fb,
                radius: rb,
            },
        ) => {
            let dr = r(*ra)? - r(*rb)?;
            let slack = &tol - if dr < zero { -dr } else { dr };
            let d = minus(
                &exact(fa.origin().to_array())?,
                &exact(fb.origin().to_array())?,
            );
            Some(slack >= zero && dot(&d, &d) <= &slack * &slack)
        }
        // A cone (S8c): axes parallel and half-angles equal to 1e-12, the
        // piece's origin on the whole's axis and its radius the whole's at
        // that height, both within the tolerance (the surfaces then agree
        // within it over a piece of the whole's extent).
        (
            Surface::Cone {
                frame: fa,
                radius: ra,
                half_angle: ha,
            },
            Surface::Cone {
                frame: fb,
                radius: rb,
                half_angle: hb,
            },
        ) => {
            let (na, nb) = (fa.normal(), fb.normal());
            let d = fa.origin() - fb.origin();
            let along = d.dot(nb);
            let off = (d - nb * along).length();
            let radius = rb + along * hb.tan();
            Some(
                na.cross(nb).length() <= 1e-12
                    && na.dot(nb) > 0.0
                    && (ha - hb).abs() <= 1e-12
                    && off <= tol_f
                    && (ra - radius).abs() <= tol_f,
            )
        }
        // A torus (S8d): axes parallel to 1e-12, centres, major and minor
        // radii within the tolerance.
        (
            Surface::Torus {
                frame: fa,
                major: ma,
                minor: ra,
            },
            Surface::Torus {
                frame: fb,
                major: mb,
                minor: rb,
            },
        ) => {
            let (na, nb) = (fa.normal(), fb.normal());
            Some(
                na.cross(nb).length() <= 1e-12
                    && na.dot(nb) > 0.0
                    && fa.origin().distance(fb.origin()) <= tol_f
                    && (ma - mb).abs() <= tol_f
                    && (ra - rb).abs() <= tol_f,
            )
        }
        // A spline wall (S8b): both extrusions of a curve (degree 1 across
        // two rows of poles) along parallel directions, the piece's rows on
        // the whole's extruded surface over the piece's `u` range.
        // A wall and its reversal in u with the other sense are one
        // oriented surface (a boundary traced the other way round).
        (Surface::BSpline(sa), Surface::BSpline(sb)) => Some(if oa == ob {
            on_extrusion(sa, sb, tol_f)
        } else {
            reversed_u(sb).is_some_and(|rb| on_extrusion(sa, &rb, tol_f))
        }),
        _ => Some(false),
    }
}

/// A spline curve traced backwards: poles reversed, knots mirrored in their
/// range (a parameter `t` becomes `a + b - t`).
fn reversed_curve(c: &crate::BSplineCurve3) -> Option<crate::BSplineCurve3> {
    let (a, b) = (c.knots()[0], *c.knots().last()?);
    let knots: Vec<f64> = c.knots().iter().rev().map(|k| a + b - k).collect();
    let mults: Vec<usize> = c.multiplicities().iter().rev().copied().collect();
    let poles: Vec<Point3> = c.poles().iter().rev().copied().collect();
    let weights: Vec<f64> = c.weights().iter().rev().copied().collect();
    crate::BSplineCurve3::new(c.degree(), poles, Some(weights), knots, mults).ok()
}

/// A spline surface with its u parameter reversed: rows in the opposite
/// order, knots mirrored in their range.
fn reversed_u(s: &crate::BSplineSurface3) -> Option<crate::BSplineSurface3> {
    let u = s.u_knots();
    let (n, m) = (u.pole_count(), s.v_knots().pole_count());
    let (a, b) = (u.knots()[0], *u.knots().last()?);
    let knots: Vec<f64> = u.knots().iter().rev().map(|k| a + b - k).collect();
    let mults: Vec<usize> = u.multiplicities().iter().rev().copied().collect();
    let mut poles = Vec::with_capacity(n * m);
    let mut weights = Vec::with_capacity(n * m);
    for i in (0..n).rev() {
        for j in 0..m {
            poles.push(s.poles()[i * m + j]);
            weights.push(s.weights()[i * m + j]);
        }
    }
    crate::BSplineSurface3::new(
        crate::KnotVector::new(u.degree(), knots, mults).ok()?,
        s.v_knots().clone(),
        poles,
        Some(weights),
    )
    .ok()
}

/// A nonrational spline curve's exact Bézier arcs over `[a, b]` of its
/// parameter, as control points.
fn arcs_over(curve: &crate::BSplineCurve3, a: f64, b: f64) -> Option<Vec<Vec<[R; 3]>>> {
    if curve.is_rational() {
        return None;
    }
    Some(
        curve
            .bezier_arcs_in(a, b)
            .ok()?
            .iter()
            .map(|arc| {
                arc.homogeneous_poles()
                    .iter()
                    .map(|h| [&h[0] / &h[3], &h[1] / &h[3], &h[2] / &h[3]])
                    .collect()
            })
            .collect(),
    )
}

/// A spline curve piece over `[a, b]` of its parameter lies on a whole
/// spline curve: over the same range their exact Bézier arcs pair up and
/// every control point differs by at most `tol` (a restriction keeps its
/// parameter; the difference of two arcs is the arc of the differences).
/// With `across`, only the parts of the differences across it count.
fn arcs_within(
    piece: &crate::BSplineCurve3,
    whole: &crate::BSplineCurve3,
    range: [f64; 2],
    tol: f64,
    across: Option<[f64; 3]>,
) -> bool {
    let (w0, w1) = whole.domain();
    if !(w0 <= range[0] && range[1] <= w1) {
        return false;
    }
    let (Some(pa), Some(wa)) = (
        arcs_over(piece, range[0], range[1]),
        arcs_over(whole, range[0], range[1]),
    ) else {
        return false;
    };
    if pa.len() != wa.len() {
        return false;
    }
    let Some(tol) = r(tol) else { return false };
    let along = across.and_then(exact);
    pa.iter().zip(&wa).all(|(x, y)| {
        x.len() == y.len()
            && x.iter().zip(y).all(|(p, q)| {
                let d = minus(p, q);
                let off2 = match &along {
                    // |d|^2 |n|^2 - (d.n)^2 <= tol^2 |n|^2.
                    Some(n) => {
                        let nn = dot(n, n);
                        let dn = dot(&d, n);
                        (dot(&d, &d) * &nn - &dn * &dn, nn)
                    }
                    None => (dot(&d, &d), R::from_integer(1.into())),
                };
                off2.0 <= &tol * &tol * off2.1
            })
    })
}

/// A spline surface's rows of poles along `u` as curves, when it is an
/// extrusion (degree 1 over two rows, nonrational): the rows and the
/// extrusion vector (the mean of the row differences).
fn extrusion(s: &crate::BSplineSurface3) -> Option<([crate::BSplineCurve3; 2], [f64; 3])> {
    let v = s.v_knots();
    if s.is_rational() || v.degree() != 1 || v.pole_count() != 2 {
        return None;
    }
    let u = s.u_knots();
    let row = |j: usize| {
        let poles: Vec<Point3> = (0..u.pole_count()).map(|i| s.poles()[i * 2 + j]).collect();
        crate::BSplineCurve3::new(
            u.degree(),
            poles,
            None,
            u.knots().to_vec(),
            u.multiplicities().to_vec(),
        )
        .ok()
    };
    let n = u.pole_count() as f64;
    let mut d = [0.0; 3];
    for i in 0..u.pole_count() {
        let (a, b) = (s.poles()[i * 2], s.poles()[i * 2 + 1]);
        d[0] += (b.x - a.x) / n;
        d[1] += (b.y - a.y) / n;
        d[2] += (b.z - a.z) / n;
    }
    Some(([row(0)?, row(1)?], d))
}

/// A spline wall piece lies on a whole one: both extrusions along
/// directions parallel to 1e-12, and each of the piece's rows within `tol`
/// of the whole's first row across the direction over the piece's range.
fn on_extrusion(piece: &crate::BSplineSurface3, whole: &crate::BSplineSurface3, tol: f64) -> bool {
    let (Some((rows, dp)), Some((base, dw))) = (extrusion(piece), extrusion(whole)) else {
        return false;
    };
    let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let cross = [
        dp[1] * dw[2] - dp[2] * dw[1],
        dp[2] * dw[0] - dp[0] * dw[2],
        dp[0] * dw[1] - dp[1] * dw[0],
    ];
    if norm(cross) > 1e-12 * norm(dp) * norm(dw) {
        return false;
    }
    let (u0, u1) = piece.domain().0;
    rows.iter()
        .all(|row| arcs_within(row, &base[0], [u0, u1], tol, Some(dw)))
}

/// `piece` lies on the support (infinite curve, surface or region kind) of
/// `whole`; `edges` holds the piece's boundary edges.
fn same_support(piece: &EntityInfo, whole: &Geometry, tol: f64, edges: &Entities) -> bool {
    match (&piece.geometry, whole) {
        (Geometry::Point(p), Geometry::Point(q)) => near(*p, *q, tol),
        (
            Geometry::Curve(Curve3::LineSegment { start, end }),
            Geometry::Curve(Curve3::LineSegment { start: a, end: b }),
        ) => on_line(*start, *a, *b, tol) && on_line(*end, *a, *b, tol),
        // A spline edge (S8b): a restriction of the whole's curve.
        // Or of its reversal: the whole's curve traced the other way round
        // (poles reversed, knots mirrored), at the piece's parameters.
        (Geometry::Curve(Curve3::BSpline(a)), Geometry::Curve(Curve3::BSpline(b))) => {
            arcs_within(a.curve(), b.curve(), a.range(), tol, None)
                || reversed_curve(b.curve())
                    .is_some_and(|rb| arcs_within(a.curve(), &rb, a.range(), tol, None))
        }
        // One circle: centres, normals (either way round) and radii within
        // tolerance, whatever the frames' x axes (a sphere's rim split at its
        // seam, S9d.1).
        (Geometry::Curve(c), Geometry::Curve(d)) => match (circle_of(c), circle_of(d)) {
            (Some((f, r)), Some((g, s))) => {
                let (m, n) = (f.normal(), g.normal());
                near(f.origin(), g.origin(), tol)
                    && (r - s).abs() <= tol
                    && m.cross(n).length() * r.max(s) <= tol
            }
            // One ellipse (S9e.1: a given result's section edge split by
            // another Boolean): each one's centre and axes' ends on the
            // other within tolerance, whatever the frames' axes' signs.
            _ => match (ellipse_of(c), ellipse_of(d)) {
                (Some(a), Some(b)) => on_ellipse(&a, &b, tol) && on_ellipse(&b, &a, tol),
                _ => false,
            },
        },
        (Geometry::Surface { .. }, Geometry::Surface { .. }) => {
            same_surface(piece, whole, tol, edges).unwrap_or(false)
        }
        (Geometry::Region(a), Geometry::Region(b)) => a == b,
        _ => false,
    }
}

type Entities<'a> = BTreeMap<EntityId, (&'a EntityInfo, f64)>;

/// All entities of several bodies with their tolerances; repeated ids reported.
fn collect<'a>(
    sets: &'a [EntitySet],
    add: &mut dyn FnMut(HistoryIssueKind, EntityId),
) -> Entities<'a> {
    let mut all = BTreeMap::new();
    for set in sets {
        for (id, info) in &set.entities {
            if all.insert(*id, (info, set.tolerance.linear())).is_some() {
                add(HistoryIssueKind::DuplicateId, *id);
            }
        }
    }
    all
}

/// Every issue of the history contract; empty means the history is valid for
/// these inputs and outputs. Sorted and duplicate-free.
pub fn check(inputs: &[EntitySet], outputs: &[EntitySet], history: &History) -> Vec<HistoryIssue> {
    use HistoryIssueKind as K;
    let mut issues = BTreeSet::new();
    let mut add = |kind, id| {
        issues.insert(HistoryIssue { kind, id });
    };
    let ins = collect(inputs, &mut add);
    let outs = collect(outputs, &mut add);
    let bodies = |sets: &[EntitySet]| sets.iter().map(|s| s.body).collect::<Vec<_>>();
    for (given, recorded) in [
        (bodies(inputs), &history.input_bodies),
        (bodies(outputs), &history.output_bodies),
    ] {
        if given != *recorded {
            for body in given.iter().chain(recorded.iter()) {
                if !given.contains(body)
                    || !recorded.contains(body)
                    || given.len() != recorded.len()
                {
                    add(K::BodyMismatch, *body);
                }
            }
        }
    }
    // H8: a single operation is one step naming it; a composite lists every
    // step of at least two operations, none itself a composite, ending with
    // the operation the composite is named after.
    let steps_ok = match (history.kind, history.steps.as_slice()) {
        (OperationKind::Composite, steps) => {
            steps.len() >= 2
                && steps.iter().all(|s| s.kind != OperationKind::Composite)
                && steps.last().map(|s| s.operation) == Some(history.operation)
        }
        (kind, [step]) => step.kind == kind && step.operation == history.operation,
        _ => false,
    };
    if !steps_ok {
        let body = history
            .output_bodies
            .first()
            .or(history.input_bodies.first())
            .copied()
            .unwrap_or(EntityId([0; 16]));
        add(K::StepsInvalid, body);
    }
    let keys: Vec<Vec<u8>> = history.relations.iter().map(Relation::sort_key).collect();
    for (k, pair) in keys.windows(2).enumerate() {
        if pair[0] >= pair[1] {
            let relation = &history.relations[k + 1];
            let id = relation
                .sources()
                .first()
                .or(relation.targets().first())
                .copied()
                .unwrap_or(EntityId([0; 16]));
            add(K::RelationOrder, id);
        }
    }
    let mut source_count: BTreeMap<EntityId, usize> = BTreeMap::new();
    let mut targeted = BTreeSet::new();
    let composite = history.kind == OperationKind::Composite;
    for relation in &history.relations {
        for s in relation.sources() {
            *source_count.entry(s).or_default() += 1;
            if !ins.contains_key(&s) {
                add(K::DanglingId, s);
            }
        }
        for t in relation.targets() {
            targeted.insert(t);
            if !outs.contains_key(&t) {
                add(K::DanglingId, t);
            }
        }
        let input = |id: &EntityId| ins.get(id).map(|x| x.0);
        let output = |id: &EntityId| outs.get(id).map(|x| x.0);
        match relation {
            Relation::Unchanged { id } => {
                if let (Some(a), Some(b)) = (input(id), output(id)) {
                    if a.kind != b.kind {
                        add(K::DimensionMismatch, *id);
                    } else if a.geometry != b.geometry || a.structure != b.structure {
                        add(K::UnchangedGeometryDiffers, *id);
                    }
                }
            }
            Relation::Modified { from, to } => {
                if from != to && !composite {
                    add(K::ModifiedIdChanged, *from);
                }
                if let (Some(a), Some(b)) = (input(from), output(to)) {
                    if a.kind != b.kind {
                        add(K::DimensionMismatch, *from);
                    } else if a.geometry.family() != b.geometry.family() {
                        add(K::ModifiedFamilyDiffers, *from);
                    }
                }
            }
            Relation::Generated { from, to, .. } => {
                if from.is_empty() {
                    add(K::InvalidArity, *to);
                }
                let target = output(to).map(|b| b.kind.dimension());
                for p in from {
                    let parent_dim = match p {
                        Parent::Entity(e) => {
                            if !ins.contains_key(e) {
                                add(K::DanglingId, *e);
                            }
                            input(e).map(|a| a.kind.dimension())
                        }
                        other => other.profile_dimension(),
                    };
                    // A sweep adds at most one dimension; intersections lower it.
                    if let (Some(t), Some(p)) = (target, parent_dim) {
                        if t > p + 1 {
                            add(K::DimensionMismatch, *to);
                        }
                    }
                }
            }
            Relation::Split { from, into } => {
                if into.len() < 2 {
                    add(K::InvalidArity, *from);
                }
                let tol = ins.get(from).map_or(0.0, |x| x.1);
                if let Some(a) = input(from) {
                    for (k, child) in into.iter().enumerate() {
                        if let Some(b) = output(child) {
                            if a.kind != b.kind {
                                add(K::DimensionMismatch, *child);
                            } else if !same_support(b, &a.geometry, tol, &outs) {
                                add(K::SplitSupportDiffers, *child);
                            }
                            if !composite && b.ordinal != k as u32 {
                                add(K::OrdinalNotCanonical, *child);
                            }
                        }
                    }
                }
                let distinct: BTreeSet<_> = into.iter().collect();
                if distinct.len() != into.len() {
                    add(K::InvalidArity, *from);
                }
            }
            Relation::Merged { from, into } => {
                if from.len() < 2 {
                    add(K::InvalidArity, *into);
                }
                let tol = outs.get(into).map_or(0.0, |x| x.1);
                if let Some(b) = output(into) {
                    for parent in from {
                        if let Some(a) = input(parent) {
                            if a.kind != b.kind {
                                add(K::DimensionMismatch, *parent);
                            } else if !same_support(a, &b.geometry, tol, &ins) {
                                add(K::MergedSupportDiffers, *parent);
                            }
                        }
                    }
                }
            }
            Relation::Deleted { id } => {
                if outs.contains_key(id) {
                    add(K::DeletedStillPresent, *id);
                }
            }
        }
    }
    for id in ins.keys() {
        match source_count.get(id) {
            None => add(K::MissingSource, *id),
            Some(n) if *n > 1 => add(K::DuplicateSource, *id),
            _ => {}
        }
    }
    for id in outs.keys() {
        if !targeted.contains(id) {
            add(K::MissingTarget, *id);
        }
    }
    issues.into_iter().collect()
}
