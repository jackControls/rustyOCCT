//! Attributes through operations (M4 of `IDENTITY_AND_HISTORY.md`, contract
//! 4): an operation context carries the policy table, and every input
//! attribute gets one recorded outcome derived from the relation its entity
//! is in. The independent `attributes::check` verifies every result in debug
//! builds.
use super::Solid;
use crate::attributes::{
    self, Attribute, AttributeKey, AttributeMap, AttributeOutcome, OnMerge, OnModify, OnSplit,
    OnTransform, OutcomeResult, PolicyTable,
};
use crate::history::{History, Relation};
use crate::identity::{AlgorithmLevel, EntityId, OperationId, OperationKind};
use crate::{Error, Result};
use std::collections::BTreeMap;

/// Computes a new value for an attribute under an `OnTransform::Recompute`
/// policy: the key, the entity (whose id the transform keeps) and the old
/// value. The kernel records the result; it never interprets values.
pub type Recompute<'a> = &'a dyn Fn(AttributeKey, EntityId, &[u8]) -> Vec<u8>;

/// The context of one operation: its id, the algorithm level it runs at, the
/// attribute policy of every key its inputs carry, and the application's
/// recompute callback.
#[derive(Clone, Copy)]
pub struct Context<'a> {
    pub operation: OperationId,
    pub level: AlgorithmLevel,
    pub policies: &'a PolicyTable,
    pub recompute: Option<Recompute<'a>>,
}

static NO_POLICIES: PolicyTable = PolicyTable(BTreeMap::new());

impl<'a> Context<'a> {
    /// The current level, no policies (bodies without attributes) and no
    /// recompute callback.
    pub fn new(operation: OperationId) -> Self {
        Self {
            operation,
            level: AlgorithmLevel::CURRENT,
            policies: &NO_POLICIES,
            recompute: None,
        }
    }
    pub fn at(self, level: AlgorithmLevel) -> Self {
        Self { level, ..self }
    }
    pub fn with_policies(self, policies: &'a PolicyTable) -> Self {
        Self { policies, ..self }
    }
    pub fn with_recompute(self, recompute: Recompute<'a>) -> Self {
        Self {
            recompute: Some(recompute),
            ..self
        }
    }
}

impl std::fmt::Debug for Context<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Context")
            .field("operation", &self.operation)
            .field("level", &self.level)
            .field("policies", self.policies)
            .field("recompute", &self.recompute.is_some())
            .finish()
    }
}

fn value_of(map: &AttributeMap, id: EntityId, key: AttributeKey) -> Option<&Vec<u8>> {
    map.get(&id)?
        .iter()
        .find(|a| a.key == key)
        .map(|a| &a.value)
}

/// All attributes of several bodies (their ids are disjoint).
fn union(bodies: &[&Solid]) -> AttributeMap {
    bodies
        .iter()
        .flat_map(|b| b.topology.attributes().clone())
        .collect()
}

/// Record every input attribute's outcome in `history` and return each
/// output's attributes, in output order.
pub(crate) fn propagate(
    context: &Context,
    inputs: &[&Solid],
    history: &mut History,
    outputs: &[&Solid],
) -> Result<Vec<AttributeMap>> {
    use OutcomeResult::*;
    let before = union(inputs);
    let mut after = vec![AttributeMap::new(); outputs.len()];
    let relation_of: BTreeMap<EntityId, &Relation> = history
        .relations
        .iter()
        .flat_map(|r| r.sources().into_iter().map(move |s| (s, r)))
        .collect();
    let transform = history.kind == OperationKind::Transform;
    let mut outcomes = Vec::new();
    for (id, list) in &before {
        for a in list {
            let policy = context
                .policies
                .0
                .get(&a.key)
                .ok_or(Error::MissingAttributePolicy(a.key.0))?;
            let relation = relation_of.get(id).ok_or(Error::InvalidTopology(
                "an attribute on an id with no relation",
            ))?;
            let keep = |to: EntityId| (vec![to], Kept, Some(a.value.clone()));
            let dropped = (Vec::new(), Dropped, None);
            let (to, result, value) = match relation {
                Relation::Unchanged { id } => keep(*id),
                Relation::Modified { to, .. } if transform => match policy.on_transform {
                    OnTransform::Keep => keep(*to),
                    OnTransform::Drop => dropped,
                    OnTransform::Recompute => {
                        let f = context.recompute.ok_or(Error::MissingRecompute(a.key.0))?;
                        (vec![*to], Recomputed, Some(f(a.key, *id, &a.value)))
                    }
                },
                Relation::Modified { to, .. } => match policy.on_modify {
                    OnModify::Keep => keep(*to),
                    OnModify::Drop => dropped,
                },
                Relation::Split { into, .. } => match policy.on_split {
                    OnSplit::CopyToAll => (into.clone(), Copied, Some(a.value.clone())),
                    OnSplit::Drop => dropped,
                },
                Relation::Merged { from, into } => match policy.on_merge {
                    OnMerge::Drop => dropped,
                    // A parent without the key counts as unequal.
                    OnMerge::KeepIfAllEqual => {
                        if from
                            .iter()
                            .all(|p| value_of(&before, *p, a.key) == Some(&a.value))
                        {
                            keep(*into)
                        } else {
                            (Vec::new(), Conflict, None)
                        }
                    }
                },
                Relation::Deleted { .. } => dropped,
                Relation::Generated { .. } => unreachable!("generated relations have no source"),
            };
            if let Some(value) = value {
                for t in &to {
                    let k = outputs
                        .iter()
                        .position(|o| o.topology.slot_of(*t).is_some())
                        .ok_or(Error::InvalidTopology("an attribute target in no output"))?;
                    let list = after[k].entry(*t).or_default();
                    if !list.iter().any(|b| b.key == a.key) {
                        list.push(Attribute {
                            key: a.key,
                            value: value.clone(),
                        });
                    }
                }
            }
            outcomes.push(AttributeOutcome {
                key: a.key,
                from: *id,
                to,
                result,
            });
        }
    }
    for map in &mut after {
        for list in map.values_mut() {
            list.sort();
        }
    }
    outcomes.sort();
    history.attributes = outcomes;
    Ok(after)
}

/// The independent attribute checker, in debug builds.
pub(crate) fn debug_check_attributes(
    context: &Context,
    inputs: &[&Solid],
    outputs: &[&Solid],
    history: &History,
) {
    if cfg!(debug_assertions) {
        let issues = attributes::check(&union(inputs), &union(outputs), history, context.policies);
        assert!(
            issues.is_empty(),
            "attribute outcomes are invalid: {issues:?}"
        );
    }
}

impl Solid {
    /// This solid with `attribute` on entity `id`, replacing any value under
    /// the same key; an id the body does not have is an error.
    pub fn with_attribute(&self, id: EntityId, attribute: Attribute) -> Result<Self> {
        let mut out = self.clone();
        out.topology = out.topology.with_attribute(id, attribute)?;
        Ok(out)
    }
    /// Every entity's attributes, by id.
    pub fn attributes(&self) -> &AttributeMap {
        self.topology.attributes()
    }
}
