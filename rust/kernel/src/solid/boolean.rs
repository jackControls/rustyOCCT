//! S9a: Booleans of two prisms in one frame (`BOOLEAN.md`, S9 of
//! REVIEW_NOTES.md).
//!
//! Both frames' axes are equal bit for bit and their origins differ by a
//! vector whose coordinates in them are binary64, decided exactly; the
//! tool's profile, translated exactly, and the object's are combined in 2D
//! (`profile::boolean`). S9a.1 builds the results that are one prism over
//! one height range: every common (the ranges' intersection), a cut whose
//! tool spans the object's heights (or misses them: the object), a fuse of
//! equal ranges (or of disjoint ranges or profiles: both inputs). Every
//! other stack of slabs is S9a.2 (`OutOfDomain`).
//!
//! Names follow provenance. Each result entity continues the input
//! entities it is a part of (a cap from the caps at its height of the
//! inputs whose material it holds, a wall or cap edge from the walls or cap
//! edges its segment comes from, a vertical edge or cap vertex from its
//! vertex's); one continuing nothing is `Generated` from the entities it
//! lies on (a new vertical edge from the two walls crossing there, a cap
//! edge inside a wall from that wall and the other input's cap). An input
//! continued by one result entity alone is `Unchanged` or `Modified`, by
//! several each continuing it alone `Split`; several inputs continued by
//! one result alone are `Merged`; where several inputs reach several
//! results (coplanar caps over several results, a wall shared in part), no
//! single relation says it: each such result is `Generated` from its
//! parents and those inputs are `Deleted`. An input continued by nothing is
//! `Deleted`.
use super::{replayable, Construction, Context, Solid};
use crate::history::{self, History, Relation};
use crate::identity::{
    Derivation, EntityId, EntityKind, OperationId, OperationKind, Parent, ProfileElement, Role,
};
use crate::profile::boolean::{boolean, Op2, Operand, Origin2, PId, SegRef};
use crate::profile::BoundaryKind;
use crate::solid::split::{q, rational_f64, zero};
use crate::topology::{End, Place, Slot, Topology};
use crate::{Boundary, Error, Point2, Profile, Result, Segment};
use num_rational::BigRational as R;
use std::collections::{BTreeMap, BTreeSet};

/// A result slot's plan: the input entities it continues and those it
/// lies on, its kind and role.
type SlotPlan = (Slot, Vec<EntityId>, Vec<EntityId>, EntityKind, Role);

/// Where a prism's slot sits: its end (low, high) or swept, and what.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum At {
    Low,
    High,
    Swept,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum What {
    Cap,
    Edge,
    Vertex,
    Wall,
    Vertical,
    Region,
}

fn at_of(place: Place) -> (At, Option<What>) {
    match place {
        Place::Low(End::Cap) => (At::Low, Some(What::Cap)),
        Place::High(End::Cap) => (At::High, Some(What::Cap)),
        Place::Low(End::Edge(_)) => (At::Low, Some(What::Edge)),
        Place::High(End::Edge(_)) => (At::High, Some(What::Edge)),
        Place::Low(End::Vertex(_)) => (At::Low, Some(What::Vertex)),
        Place::High(End::Vertex(_)) => (At::High, Some(What::Vertex)),
        Place::Swept => (At::Swept, None),
    }
}

/// A prism's entities by (end, what, boundary, element index).
fn index(solid: &Solid) -> Result<BTreeMap<(At, What, usize, usize), EntityId>> {
    let t = &solid.topology;
    let places: BTreeMap<Slot, Place> = t.layout().iter().copied().collect();
    let profile = match &solid.construction {
        Construction::Prism(p) => p,
        _ => unreachable!("a prism"),
    };
    // Labels back to stored elements.
    let mut labels = BTreeMap::new();
    for (b, boundary) in profile.boundaries().enumerate() {
        if let Some(l) = boundary.labels() {
            labels.insert(l.boundary, (b, ProfileElement::Boundary));
            for (j, x) in l.segments.iter().enumerate() {
                labels.insert(*x, (b, ProfileElement::Segment(j as u32)));
            }
            for (j, x) in l.vertices.iter().enumerate() {
                labels.insert(*x, (b, ProfileElement::Vertex(j as u32)));
            }
        }
    }
    let mut out = BTreeMap::new();
    for (id, slot) in t.ids() {
        let d = t.derivation(id).expect("a derivation");
        let (at, what) = at_of(*places.get(&slot).unwrap_or(&Place::Swept));
        let element = match d.parents.first() {
            Some(Parent::Profile { boundary, element }) => Some((*boundary as usize, *element)),
            Some(Parent::Label(l)) => labels.get(l).copied(),
            _ => None,
        };
        let key = match (d.entity, at, what, element) {
            (EntityKind::Region, ..) => (At::Swept, What::Region, 0, 0),
            (EntityKind::Face, At::Swept, _, Some((b, ProfileElement::Segment(j)))) => {
                (At::Swept, What::Wall, b, j as usize)
            }
            (EntityKind::Edge, At::Swept, _, Some((b, ProfileElement::Vertex(j)))) => {
                (At::Swept, What::Vertical, b, j as usize)
            }
            (EntityKind::Face, at, Some(What::Cap), _) => (at, What::Cap, 0, 0),
            (EntityKind::Edge, at, Some(_), Some((b, ProfileElement::Segment(j)))) => {
                (at, What::Edge, b, j as usize)
            }
            (EntityKind::Vertex, at, Some(_), Some((b, ProfileElement::Vertex(j)))) => {
                (at, What::Vertex, b, j as usize)
            }
            _ => return Err(Error::InvalidTopology("a prism slot of unknown place")),
        };
        out.insert(key, id);
    }
    Ok(out)
}

/// `a` translated by `(du, dv)`, exactly (`None` when a coordinate would
/// round).
fn translated(profile: &Profile, du: f64, dv: f64) -> Result<Option<Profile>> {
    let exact = |x: f64, d: f64| {
        let y = x + d;
        (q(x) + q(d) == q(y)).then_some(y)
    };
    let point =
        |p: Point2| -> Option<Point2> { Some(Point2::new(exact(p.x, du)?, exact(p.y, dv)?)) };
    let tolerance = profile.tolerance();
    let mut out: Vec<Boundary> = Vec::new();
    for b in profile.boundaries() {
        let moved = match &b.kind {
            BoundaryKind::Circle { center, radius } => {
                let Some(c) = point(*center) else {
                    return Ok(None);
                };
                Boundary::circle(c, *radius, tolerance)?
            }
            BoundaryKind::Polygon(points) => {
                let Some(ps) = points.iter().map(|p| point(*p)).collect::<Option<Vec<_>>>() else {
                    return Ok(None);
                };
                Boundary::polygon(ps, tolerance)?
            }
            BoundaryKind::Path { points, segments } => {
                let Some(ps) = points.iter().map(|p| point(*p)).collect::<Option<Vec<_>>>() else {
                    return Ok(None);
                };
                let mut segs = Vec::new();
                for s in segments {
                    segs.push(match s {
                        Segment::Arc {
                            center,
                            radius,
                            ccw,
                        } => {
                            let Some(c) = point(*center) else {
                                return Ok(None);
                            };
                            Segment::Arc {
                                center: c,
                                radius: *radius,
                                ccw: *ccw,
                            }
                        }
                        s => s.clone(),
                    });
                }
                Boundary::path(ps, segs, tolerance)?
            }
        };
        out.push(moved);
    }
    let holes = out.split_off(1);
    Ok(Some(Profile::new(out.remove(0), holes, tolerance)?))
}

/// Whether an exact rational is a binary64, and which.
fn as_f64(x: &R) -> Option<f64> {
    let f = rational_f64(x);
    (q(f) == *x).then_some(f)
}

impl Solid {
    /// `self ∪ other` (S9a): the result's solids, each a maximal connected
    /// solid region, and the history. `OutOfDomain` until S9b for solids
    /// other than prisms in frames with equal axes, and until S9a.2 for a
    /// result that is a stack of prisms of different profiles.
    pub fn fuse(&self, operation: OperationId, other: &Solid) -> Result<(Vec<Solid>, History)> {
        self.boolean_in(&Context::new(operation), other, Op2::Fuse)
    }
    /// `self − tool` (S9a).
    pub fn cut(&self, operation: OperationId, tool: &Solid) -> Result<(Vec<Solid>, History)> {
        self.boolean_in(&Context::new(operation), tool, Op2::Cut)
    }
    /// `self ∩ other` (S9a).
    pub fn common(&self, operation: OperationId, other: &Solid) -> Result<(Vec<Solid>, History)> {
        self.boolean_in(&Context::new(operation), other, Op2::Common)
    }
    /// [`Solid::fuse`] in an operation context.
    pub fn fuse_in(&self, context: &Context, other: &Solid) -> Result<(Vec<Solid>, History)> {
        self.boolean_in(context, other, Op2::Fuse)
    }
    /// [`Solid::cut`] in an operation context.
    pub fn cut_in(&self, context: &Context, tool: &Solid) -> Result<(Vec<Solid>, History)> {
        self.boolean_in(context, tool, Op2::Cut)
    }
    /// [`Solid::common`] in an operation context.
    pub fn common_in(&self, context: &Context, other: &Solid) -> Result<(Vec<Solid>, History)> {
        self.boolean_in(context, other, Op2::Common)
    }

    fn boolean_in(
        &self,
        context: &Context,
        other: &Solid,
        op: Op2,
    ) -> Result<(Vec<Solid>, History)> {
        replayable(context.level)?;
        let (Construction::Prism(pa), Construction::Prism(pb)) =
            (&self.construction, &other.construction)
        else {
            return Err(Error::OutOfDomain(
                "a Boolean of solids other than prisms (S9b on)",
            ));
        };
        let (fa, fb) = (self.frame, other.frame);
        if fa.x() != fb.x() || fa.y() != fb.y() || fa.normal() != fb.normal() {
            return Err(Error::OutOfDomain(
                "a Boolean of prisms in frames with different axes (S9b)",
            ));
        }
        // The origins' offset in the axes, exactly (Cramer's rule).
        let v = |x: [f64; 3]| x.map(q);
        let (x, y, n) = (
            v(fa.x().to_array()),
            v(fa.y().to_array()),
            v(fa.normal().to_array()),
        );
        let d: [R; 3] = {
            let (ob, oa) = (fb.origin().to_array(), fa.origin().to_array());
            std::array::from_fn(|i| q(ob[i]) - q(oa[i]))
        };
        let det3 = |a: &[R; 3], b: &[R; 3], c: &[R; 3]| {
            &a[0] * (&b[1] * &c[2] - &b[2] * &c[1]) - &a[1] * (&b[0] * &c[2] - &b[2] * &c[0])
                + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0])
        };
        let det = det3(&x, &y, &n);
        if det == zero() {
            return Err(Error::Degenerate("a frame's axes"));
        }
        let (du, dv, dw) = (
            det3(&d, &y, &n) / &det,
            det3(&x, &d, &n) / &det,
            det3(&x, &y, &d) / &det,
        );
        let (Some(du), Some(dv), Some(dw)) = (as_f64(&du), as_f64(&dv), as_f64(&dw)) else {
            return Err(Error::OutOfDomain(
                "prisms whose frames' offset is not binary64 in their axes (S9b)",
            ));
        };
        let Some(pb) = translated(pb, du, dv)? else {
            return Err(Error::OutOfDomain(
                "a tool whose profile does not translate exactly (S9b)",
            ));
        };
        let heights = |s: &Solid, dw: f64| -> Option<[f64; 2]> {
            let (lo, hi) = (s.start.min(s.end), s.start.max(s.end));
            let (a, b) = (lo + dw, hi + dw);
            (q(lo) + q(dw) == q(a) && q(hi) + q(dw) == q(b)).then_some([a, b])
        };
        let ha = heights(self, 0.0).expect("its own heights");
        let Some(hb) = heights(other, dw) else {
            return Err(Error::OutOfDomain(
                "a tool whose heights do not translate exactly (S9b)",
            ));
        };
        // S9a.1: the result's profiles and height range, or a shortcut.
        let whole = |s: &Solid| -> Result<Vec<Solid>> { Ok(vec![s.clone()]) };
        let (range, results) = match op {
            Op2::Common => {
                let (lo, hi) = (ha[0].max(hb[0]), ha[1].min(hb[1]));
                if hi <= lo {
                    return self.finish(context, other, op, Vec::new(), Vec::new());
                }
                ([lo, hi], boolean(pa, &pb, op)?)
            }
            Op2::Cut => {
                if hb[1] <= ha[0] || hb[0] >= ha[1] {
                    return self.finish(context, other, op, whole(self)?, Vec::new());
                }
                if !(hb[0] <= ha[0] && hb[1] >= ha[1]) {
                    return Err(Error::OutOfDomain(
                        "a cut leaving a stack of prisms (S9a.2)",
                    ));
                }
                (ha, boolean(pa, &pb, op)?)
            }
            Op2::Fuse => {
                if ha == hb {
                    (ha, boolean(pa, &pb, op)?)
                } else if hb[1] < ha[0] || hb[0] > ha[1] {
                    return self.finish(
                        context,
                        other,
                        op,
                        vec![self.clone(), other.clone()],
                        Vec::new(),
                    );
                } else {
                    // Different ranges: separate only when the profiles
                    // share no area and do not touch.
                    let common = boolean(pa, &pb, Op2::Common)?;
                    let fused = boolean(pa, &pb, Op2::Fuse)?;
                    if common.pieces.is_empty()
                        && fused.pieces.len() == 2
                        && hb[1] != ha[0]
                        && hb[0] != ha[1]
                    {
                        return self.finish(
                            context,
                            other,
                            op,
                            vec![self.clone(), other.clone()],
                            Vec::new(),
                        );
                    }
                    return Err(Error::OutOfDomain(
                        "a fuse of prisms of different heights (S9a.2)",
                    ));
                }
            }
        };
        // Which inputs' material each result holds.
        let mut built = Vec::new();
        for piece in &results.pieces {
            let solid = Solid::build(
                context.operation,
                piece.profile.clone(),
                fa,
                range[0],
                range[1],
            )?;
            let mut members: BTreeSet<Operand> = BTreeSet::new();
            match op {
                Op2::Common => {
                    members.extend([Operand::A, Operand::B]);
                }
                Op2::Cut => {
                    members.insert(Operand::A);
                }
                Op2::Fuse => {
                    for o in piece.segments.iter().flatten() {
                        members.insert(o.from.0);
                        if let Some(s) = o.shared {
                            members.insert(s.0);
                        }
                    }
                }
            }
            built.push((solid, piece, members));
        }
        // A fuse input with no boundary in any result lies inside the other
        // (connected) input's result.
        if op == Op2::Fuse {
            for (operand, other_op) in [(Operand::A, Operand::B), (Operand::B, Operand::A)] {
                if built.iter().any(|(_, _, m)| m.contains(&operand)) {
                    continue;
                }
                let owner = built
                    .iter()
                    .position(|(_, _, m)| m.contains(&other_op))
                    .ok_or(Error::InvalidTopology("a fuse input outside its result"))?;
                built[owner].2.insert(operand);
            }
        }
        let alias = results_alias(pa, &pb);
        let heights_of = |o: Operand| if o == Operand::A { ha } else { hb };
        let (index_a, index_b) = (index(self)?, index(other)?);
        let input = |o: Operand| if o == Operand::A { &index_a } else { &index_b };
        // Each result slot's continued and touched inputs.
        let mut plans: Vec<Vec<SlotPlan>> = Vec::new();
        for (solid, piece, members) in &built {
            let t = &solid.topology;
            let places: BTreeMap<Slot, Place> = t.layout().iter().copied().collect();
            let mut named = Vec::new();
            for (id, slot) in t.ids() {
                let d = t.derivation(id).expect("a derivation");
                let (at, what) = at_of(*places.get(&slot).unwrap_or(&Place::Swept));
                let h = match at {
                    At::Low => Some(range[0]),
                    At::High => Some(range[1]),
                    At::Swept => None,
                };
                // An input's end at this height, if any.
                let end_of = |o: Operand| -> Option<At> {
                    let [lo, hi] = heights_of(o);
                    match h {
                        Some(h) if h == lo => Some(At::Low),
                        Some(h) if h == hi => Some(At::High),
                        _ => None,
                    }
                };
                let element = match d.parents.first() {
                    Some(Parent::Profile { boundary, element }) => {
                        Some((*boundary as usize, *element))
                    }
                    _ => None,
                };
                let (mut continues, mut touches) = (Vec::new(), Vec::new());
                // A cut's tool faces the other way where it bounds the
                // result: what it gives is generated, not continued.
                let keep = |o: Operand,
                            id: EntityId,
                            continues: &mut Vec<EntityId>,
                            touches: &mut Vec<EntityId>| {
                    if op == Op2::Cut && o == Operand::B {
                        touches.push(id);
                    } else {
                        continues.push(id);
                    }
                };
                let seg_refs = |rb: usize, rj: usize| -> Vec<SegRef> {
                    let o: Origin2 = piece.segments[rb][rj];
                    std::iter::once(o.from).chain(o.shared).collect()
                };
                match (d.entity, what, element) {
                    (EntityKind::Region, ..) => {
                        for m in members {
                            continues.push(input(*m)[&(At::Swept, What::Region, 0, 0)]);
                        }
                    }
                    (EntityKind::Face, Some(What::Cap), _) => {
                        for m in members {
                            if let Some(e) = end_of(*m) {
                                continues.push(input(*m)[&(e, What::Cap, 0, 0)]);
                            }
                        }
                    }
                    (EntityKind::Face, None, Some((rb, ProfileElement::Segment(rj)))) => {
                        for (o, b, j, _) in seg_refs(rb, rj as usize) {
                            keep(
                                o,
                                input(o)[&(At::Swept, What::Wall, b, j)],
                                &mut continues,
                                &mut touches,
                            );
                        }
                    }
                    (EntityKind::Edge, Some(_), Some((rb, ProfileElement::Segment(rj)))) => {
                        for (o, b, j, _) in seg_refs(rb, rj as usize) {
                            match end_of(o) {
                                Some(e) => keep(
                                    o,
                                    input(o)[&(e, What::Edge, b, j)],
                                    &mut continues,
                                    &mut touches,
                                ),
                                None => {
                                    // Inside o's wall: where the result's cap
                                    // (another input's) meets it.
                                    touches.push(input(o)[&(At::Swept, What::Wall, b, j)]);
                                    for m in members {
                                        if let Some(e) = end_of(*m) {
                                            touches.push(input(*m)[&(e, What::Cap, 0, 0)]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    (EntityKind::Edge, None, Some((rb, ProfileElement::Vertex(rj)))) => {
                        match piece.vertices[rb][rj as usize] {
                            PId::Vertex(o, b, j) => {
                                keep(
                                    o,
                                    input(o)[&(At::Swept, What::Vertical, b, j)],
                                    &mut continues,
                                    &mut touches,
                                );
                                if let Some(PId::Vertex(o2, b2, j2)) = twin(o, b, j, &alias) {
                                    keep(
                                        o2,
                                        input(o2)[&(At::Swept, What::Vertical, b2, j2)],
                                        &mut continues,
                                        &mut touches,
                                    );
                                }
                            }
                            PId::Cross(k) => {
                                for (o, b, j) in results.crosses[&k] {
                                    touches.push(input(o)[&(At::Swept, What::Wall, b, j)]);
                                }
                            }
                        }
                    }
                    (EntityKind::Vertex, Some(_), Some((rb, ProfileElement::Vertex(rj)))) => {
                        let vertex = piece.vertices[rb][rj as usize];
                        let owners: Vec<(Operand, usize, usize)> = match vertex {
                            PId::Vertex(o, b, j) => {
                                let mut v = vec![(o, b, j)];
                                if let Some(PId::Vertex(o2, b2, j2)) = twin(o, b, j, &alias) {
                                    v.push((o2, b2, j2));
                                }
                                v
                            }
                            PId::Cross(_) => Vec::new(),
                        };
                        for (o, b, j) in &owners {
                            match end_of(*o) {
                                Some(e) => keep(
                                    *o,
                                    input(*o)[&(e, What::Vertex, *b, *j)],
                                    &mut continues,
                                    &mut touches,
                                ),
                                None => {
                                    touches.push(input(*o)[&(At::Swept, What::Vertical, *b, *j)]);
                                    for m in members {
                                        if let Some(e) = end_of(*m) {
                                            touches.push(input(*m)[&(e, What::Cap, 0, 0)]);
                                        }
                                    }
                                }
                            }
                        }
                        if let PId::Cross(k) = vertex {
                            for (o, b, j) in results.crosses[&k] {
                                match end_of(o) {
                                    Some(e) => touches.push(input(o)[&(e, What::Edge, b, j)]),
                                    None => touches.push(input(o)[&(At::Swept, What::Wall, b, j)]),
                                }
                            }
                        }
                    }
                    _ => {
                        return Err(Error::InvalidTopology(
                            "a result slot of unknown provenance",
                        ))
                    }
                }
                continues.sort();
                continues.dedup();
                touches.sort();
                touches.dedup();
                named.push((slot, continues, touches, d.entity, d.role));
            }
            plans.push(named);
        }
        let solids: Vec<Solid> = built.into_iter().map(|(s, ..)| s).collect();
        self.finish(context, other, op, solids, plans)
    }

    /// Names the results by their plans, and the history. With no plans the
    /// results are inputs returned whole (or none).
    fn finish(
        &self,
        context: &Context,
        other: &Solid,
        op: Op2,
        mut results: Vec<Solid>,
        plans: Vec<Vec<SlotPlan>>,
    ) -> Result<(Vec<Solid>, History)> {
        let operation = context.operation;
        let kind = match op {
            Op2::Fuse => OperationKind::Fuse,
            Op2::Cut => OperationKind::Cut,
            Op2::Common => OperationKind::Common,
        };
        let inputs = [self, other];
        let bodies: Vec<EntityId> = inputs.iter().map(|s| s.topology.body_id()).collect();
        let mut relations: Vec<Relation> = Vec::new();
        if plans.is_empty() {
            // Inputs returned whole keep every id; the others are deleted.
            let kept: BTreeSet<EntityId> = results.iter().map(|s| s.topology.body_id()).collect();
            for s in inputs {
                let whole = kept.contains(&s.topology.body_id());
                for (id, _) in s.topology.ids() {
                    relations.push(if whole {
                        Relation::Unchanged { id }
                    } else {
                        Relation::Deleted { id }
                    });
                }
            }
            relations.sort_by_cached_key(Relation::sort_key);
            let outputs = results.iter().map(|s| s.topology.body_id()).collect();
            let history = History::new(operation, kind, bodies, outputs, relations, Vec::new())
                .at_level(context.level);
            debug_check(&inputs, &results.iter().collect::<Vec<_>>(), &history);
            return Ok((results, history));
        }
        // Continuation groups: inputs and results joined by `continues`.
        let mut children: BTreeMap<EntityId, Vec<(usize, usize)>> = BTreeMap::new();
        for (k, named) in plans.iter().enumerate() {
            for (i, (_, continues, ..)) in named.iter().enumerate() {
                for c in continues {
                    children.entry(*c).or_default().push((k, i));
                }
            }
        }
        let all_inputs: Vec<EntityId> = inputs
            .iter()
            .flat_map(|s| s.topology.ids().map(|(id, _)| id))
            .collect();
        let derive = |entity, role, ordinal, parents: Vec<Parent>| Derivation {
            operation,
            kind,
            entity,
            role,
            ordinal,
            parents,
        };
        let input_derivation = |id: EntityId| -> Derivation {
            inputs
                .iter()
                .find_map(|s| s.topology.derivation(id))
                .expect("an input's derivation")
                .clone()
        };
        let mut names: Vec<Vec<Option<Derivation>>> =
            plans.iter().map(|n| vec![None; n.len()]).collect();
        let mut kept: BTreeSet<EntityId> = BTreeSet::new();
        let mut counters: BTreeMap<(EntityKind, Role, Vec<Parent>), u32> = BTreeMap::new();
        let mut fresh = |entity, role, parents: Vec<Parent>| {
            let n = counters.entry((entity, role, parents.clone())).or_insert(0);
            *n += 1;
            derive(entity, role, *n - 1, parents)
        };
        let mut done: BTreeSet<EntityId> = BTreeSet::new();
        for id in &all_inputs {
            if done.contains(id) {
                continue;
            }
            let Some(kids) = children.get(id) else {
                relations.push(Relation::Deleted { id: *id });
                done.insert(*id);
                continue;
            };
            // The group: inputs sharing any of these results.
            let mut group_in: BTreeSet<EntityId> = BTreeSet::from([*id]);
            let mut group_out: BTreeSet<(usize, usize)> = kids.iter().copied().collect();
            loop {
                let before = (group_in.len(), group_out.len());
                for &(k, i) in &group_out.clone() {
                    group_in.extend(plans[k][i].1.iter().copied());
                }
                for x in &group_in.clone() {
                    group_out.extend(children[x].iter().copied());
                }
                if (group_in.len(), group_out.len()) == before {
                    break;
                }
            }
            done.extend(group_in.iter().copied());
            let ins: Vec<EntityId> = group_in.into_iter().collect();
            let outs: Vec<(usize, usize)> = group_out.into_iter().collect();
            match (ins.len(), outs.len()) {
                (1, 1) => {
                    let (k, i) = outs[0];
                    names[k][i] = Some(input_derivation(ins[0]));
                    kept.insert(ins[0]);
                }
                (1, _) => {
                    for (ordinal, &(k, i)) in outs.iter().enumerate() {
                        let (_, _, _, entity, role) = &plans[k][i];
                        names[k][i] = Some(derive(
                            *entity,
                            *role,
                            ordinal as u32,
                            vec![Parent::Entity(ins[0])],
                        ));
                    }
                    let into = outs
                        .iter()
                        .map(|&(k, i)| names[k][i].as_ref().expect("a name").id())
                        .collect();
                    relations.push(Relation::Split { from: ins[0], into });
                }
                (_, 1) => {
                    let (k, i) = outs[0];
                    let (_, _, _, entity, role) = &plans[k][i];
                    let d = fresh(
                        *entity,
                        *role,
                        ins.iter().map(|x| Parent::Entity(*x)).collect(),
                    );
                    relations.push(Relation::Merged {
                        from: ins.clone(),
                        into: d.id(),
                    });
                    names[k][i] = Some(d);
                }
                _ => {
                    for x in &ins {
                        relations.push(Relation::Deleted { id: *x });
                    }
                    for &(k, i) in &outs {
                        let (_, continues, _, entity, role) = &plans[k][i];
                        let parents: Vec<Parent> =
                            continues.iter().map(|x| Parent::Entity(*x)).collect();
                        let d = fresh(*entity, *role, parents.clone());
                        relations.push(Relation::Generated {
                            from: parents,
                            to: d.id(),
                            role: *role,
                        });
                        names[k][i] = Some(d);
                    }
                }
            }
        }
        // Results continuing nothing: generated from what they lie on.
        for (k, named) in plans.iter().enumerate() {
            for (i, (_, continues, touches, entity, role)) in named.iter().enumerate() {
                if !continues.is_empty() {
                    continue;
                }
                let parents: Vec<Parent> = touches.iter().map(|x| Parent::Entity(*x)).collect();
                if parents.is_empty() {
                    return Err(Error::InvalidTopology(
                        "a Boolean result entity from nothing",
                    ));
                }
                let d = fresh(*entity, *role, parents.clone());
                relations.push(Relation::Generated {
                    from: parents,
                    to: d.id(),
                    role: *role,
                });
                names[k][i] = Some(d);
            }
        }
        // Rename, raise enclosures to the parents', and the relations of
        // inputs kept whole.
        let before: Vec<_> = inputs
            .iter()
            .map(|s| s.topology.entity_set(s.resolution()))
            .collect();
        for (k, solid) in results.iter_mut().enumerate() {
            let derivations: Vec<(Slot, Derivation)> = plans[k]
                .iter()
                .zip(&names[k])
                .map(|((slot, ..), d)| (*slot, d.clone().expect("every slot named")))
                .collect();
            let body = derive(
                EntityKind::Body,
                Role::Body,
                k as u32,
                bodies.iter().map(|b| Parent::Entity(*b)).collect(),
            );
            let parents_of: BTreeMap<EntityId, Vec<EntityId>> = plans[k]
                .iter()
                .zip(&names[k])
                .map(|((_, continues, ..), d)| {
                    (d.as_ref().expect("a name").id(), continues.clone())
                })
                .collect();
            let mut topology: Topology = solid.topology.clone().renamed(body, derivations)?;
            topology.raise_enclosures(|id| {
                parents_of
                    .get(&id)?
                    .iter()
                    .filter_map(|p| inputs.iter().find_map(|s| s.topology.enclosure_bound(*p)))
                    .reduce(f64::max)
            });
            solid.topology = topology;
            solid.operation = operation;
        }
        for id in &kept {
            let old = before
                .iter()
                .find_map(|set| set.entities.get(id))
                .expect("an input entity");
            let new = results
                .iter()
                .find_map(|s| {
                    s.topology
                        .entity_set(s.resolution())
                        .entities
                        .get(id)
                        .cloned()
                })
                .expect("a kept entity");
            relations.push(
                if old.geometry == new.geometry && old.structure == new.structure {
                    Relation::Unchanged { id: *id }
                } else {
                    Relation::Modified { from: *id, to: *id }
                },
            );
        }
        relations.sort_by_cached_key(Relation::sort_key);
        relations.dedup();
        let outputs = results.iter().map(|s| s.topology.body_id()).collect();
        let history = History::new(operation, kind, bodies, outputs, relations, Vec::new())
            .at_level(context.level);
        debug_check(&inputs, &results.iter().collect::<Vec<_>>(), &history);
        Ok((results, history))
    }
}

/// The other operand's vertex equal to this one (the arrangement's alias).
fn twin(
    o: Operand,
    b: usize,
    j: usize,
    alias: &BTreeMap<(Operand, usize, usize), (Operand, usize, usize)>,
) -> Option<PId> {
    alias
        .iter()
        .find_map(|(k, v)| {
            if *v == (o, b, j) {
                Some(*k)
            } else if *k == (o, b, j) {
                Some(*v)
            } else {
                None
            }
        })
        .map(|(o, b, j)| PId::Vertex(o, b, j))
}

/// B's stored vertices equal to A's, by position (as the arrangement's).
fn results_alias(
    a: &Profile,
    b: &Profile,
) -> BTreeMap<(Operand, usize, usize), (Operand, usize, usize)> {
    let vertices = |p: &Profile| -> Vec<(usize, usize, Point2)> {
        let mut out = Vec::new();
        for (bi, boundary) in p.boundaries().enumerate() {
            let points: &[Point2] = match &boundary.kind {
                BoundaryKind::Polygon(points) => points,
                BoundaryKind::Path { points, .. } => points,
                BoundaryKind::Circle { .. } => &[],
            };
            for (j, pt) in points.iter().enumerate() {
                out.push((bi, j, *pt));
            }
        }
        out
    };
    let (va, vb) = (vertices(a), vertices(b));
    let mut out = BTreeMap::new();
    for (bb, jb, pb) in &vb {
        if let Some((ba, ja, _)) = va.iter().find(|(_, _, pa)| pa == pb) {
            out.insert((Operand::B, *bb, *jb), (Operand::A, *ba, *ja));
        }
    }
    out
}

fn debug_check(inputs: &[&Solid], outputs: &[&Solid], history: &History) {
    if cfg!(debug_assertions) {
        let ins: Vec<_> = inputs
            .iter()
            .map(|s| s.topology.entity_set(s.resolution()))
            .collect();
        let outs: Vec<_> = outputs
            .iter()
            .map(|s| s.topology.entity_set(s.resolution()))
            .collect();
        let issues = history::check(&ins, &outs, history);
        assert!(issues.is_empty(), "Boolean history is invalid: {issues:?}");
    }
}
