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
    encode_parents, role_code, EntityId, EntityKind, OperationId, OperationKind, Parent, Role,
};
use crate::topology::{Curve3, Orientation, Surface};
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

/// The complete relation list of one operation (or of a composed chain).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct History {
    pub operation: OperationId,
    pub kind: OperationKind,
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
            input_bodies,
            output_bodies,
            relations,
            attributes,
        }
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
        Ok(History::new(
            next.operation,
            OperationKind::Composite,
            self.input_bodies.clone(),
            next.output_bodies.clone(),
            relations,
            Vec::new(),
        ))
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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Family {
    Point,
    Line,
    Circle,
    Plane,
    Cylinder,
}

impl Geometry {
    pub fn family(&self) -> Family {
        match self {
            Self::Point(_) => Family::Point,
            Self::Curve(Curve3::LineSegment { .. }) => Family::Line,
            Self::Curve(_) => Family::Circle,
            Self::Surface {
                surface: Surface::Plane(_),
                ..
            } => Family::Plane,
            Self::Surface { .. } => Family::Cylinder,
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
        Curve3::LineSegment { .. } => None,
    }
}

/// `piece` lies on the support (infinite curve or whole surface) of `whole`.
fn same_support(piece: &Geometry, whole: &Geometry, tol: f64) -> bool {
    match (piece, whole) {
        (Geometry::Point(p), Geometry::Point(q)) => near(*p, *q, tol),
        (
            Geometry::Curve(Curve3::LineSegment { start, end }),
            Geometry::Curve(Curve3::LineSegment { start: a, end: b }),
        ) => on_line(*start, *a, *b, tol) && on_line(*end, *a, *b, tol),
        (Geometry::Curve(c), Geometry::Curve(d)) => {
            circle_of(c).is_some() && circle_of(c) == circle_of(d)
        }
        (Geometry::Surface { .. }, Geometry::Surface { .. }) => piece == whole,
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
                            } else if !same_support(&b.geometry, &a.geometry, tol) {
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
                            } else if !same_support(&a.geometry, &b.geometry, tol) {
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
