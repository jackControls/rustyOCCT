//! Attribute propagation policy and its independent checker (see
//! `IDENTITY_AND_HISTORY.md`, contract 4).
//!
//! Attributes are opaque key/value pairs on entities. Every key has a declared
//! [`Policy`]; an operation records one [`AttributeOutcome`] per input
//! attribute in its history. [`check`] derives the outcome each policy
//! requires from the history's relations alone and compares.
use crate::history::{History, Relation};
use crate::identity::{EntityId, OperationKind};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AttributeKey(pub u64);

/// The kernel never interprets the value.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Attribute {
    pub key: AttributeKey,
    pub value: Vec<u8>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnModify {
    Keep,
    Drop,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnTransform {
    Keep,
    Drop,
    /// The application recomputes the value; the kernel records it.
    Recompute,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnSplit {
    CopyToAll,
    Drop,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnMerge {
    KeepIfAllEqual,
    Drop,
}

/// Deletion always drops; generation never creates attributes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Policy {
    pub on_modify: OnModify,
    pub on_transform: OnTransform,
    pub on_split: OnSplit,
    pub on_merge: OnMerge,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PolicyTable(pub BTreeMap<AttributeKey, Policy>);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum OutcomeResult {
    Kept,
    Copied,
    Recomputed,
    Dropped,
    /// Merged parents carried unequal values; the output has no attribute.
    Conflict,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct AttributeOutcome {
    pub key: AttributeKey,
    pub from: EntityId,
    pub to: Vec<EntityId>,
    pub result: OutcomeResult,
}

/// Attributes per entity id, in any order.
pub type AttributeMap = BTreeMap<EntityId, Vec<Attribute>>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AttributeIssueKind {
    MissingPolicy,
    MissingOutcome,
    DuplicateOutcome,
    DanglingOutcome,
    OutcomeContradictsPolicy,
    UntracedAttribute,
    ValueChanged,
    AttributeOnDropped,
    MissingAttribute,
}

impl AttributeIssueKind {
    pub fn name(self) -> &'static str {
        match self {
            Self::MissingPolicy => "missing_policy",
            Self::MissingOutcome => "missing_outcome",
            Self::DuplicateOutcome => "duplicate_outcome",
            Self::DanglingOutcome => "dangling_outcome",
            Self::OutcomeContradictsPolicy => "outcome_contradicts_policy",
            Self::UntracedAttribute => "untraced_attribute",
            Self::ValueChanged => "value_changed",
            Self::AttributeOnDropped => "attribute_on_dropped",
            Self::MissingAttribute => "missing_attribute",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AttributeIssue {
    pub kind: AttributeIssueKind,
    pub key: AttributeKey,
    pub id: EntityId,
}

fn value(map: &AttributeMap, id: EntityId, key: AttributeKey) -> Option<&Vec<u8>> {
    map.get(&id)?
        .iter()
        .find(|a| a.key == key)
        .map(|a| &a.value)
}

/// The outcome a policy requires for attribute `key` on input `from`.
pub fn expected_outcome(
    inputs: &AttributeMap,
    history: &History,
    policy: &Policy,
    key: AttributeKey,
    from: EntityId,
) -> Option<AttributeOutcome> {
    use OutcomeResult::*;
    let outcome = |to: Vec<EntityId>, result| AttributeOutcome {
        key,
        from,
        to,
        result,
    };
    let relation = history
        .relations
        .iter()
        .find(|r| r.sources().contains(&from))?;
    Some(match relation {
        Relation::Unchanged { id } => outcome(vec![*id], Kept),
        Relation::Modified { to, .. } => {
            if history.kind == OperationKind::Transform {
                match policy.on_transform {
                    OnTransform::Keep => outcome(vec![*to], Kept),
                    OnTransform::Drop => outcome(vec![], Dropped),
                    OnTransform::Recompute => outcome(vec![*to], Recomputed),
                }
            } else {
                match policy.on_modify {
                    OnModify::Keep => outcome(vec![*to], Kept),
                    OnModify::Drop => outcome(vec![], Dropped),
                }
            }
        }
        Relation::Split { into, .. } => match policy.on_split {
            OnSplit::CopyToAll => outcome(into.clone(), Copied),
            OnSplit::Drop => outcome(vec![], Dropped),
        },
        Relation::Merged {
            from: parents,
            into,
        } => match policy.on_merge {
            OnMerge::Drop => outcome(vec![], Dropped),
            OnMerge::KeepIfAllEqual => {
                let mine = value(inputs, from, key);
                if parents.iter().all(|p| value(inputs, *p, key) == mine) {
                    outcome(vec![*into], Kept)
                } else {
                    outcome(vec![], Conflict)
                }
            }
        },
        Relation::Deleted { .. } => outcome(vec![], Dropped),
        Relation::Generated { .. } => return None,
    })
}

/// Every issue of the attribute contract for one history. Sorted and
/// duplicate-free.
pub fn check(
    inputs: &AttributeMap,
    outputs: &AttributeMap,
    history: &History,
    table: &PolicyTable,
) -> Vec<AttributeIssue> {
    use AttributeIssueKind as K;
    let mut issues = BTreeSet::new();
    let mut add = |kind, key, id| {
        issues.insert(AttributeIssue { kind, key, id });
    };
    let mut recorded: BTreeMap<(EntityId, AttributeKey), Vec<&AttributeOutcome>> = BTreeMap::new();
    for o in &history.attributes {
        recorded.entry((o.from, o.key)).or_default().push(o);
        if value(inputs, o.from, o.key).is_none() {
            add(K::DanglingOutcome, o.key, o.from);
        }
    }
    // Where each output attribute may legitimately come from.
    let mut carried: BTreeMap<(EntityId, AttributeKey), Vec<(&AttributeOutcome, bool)>> =
        BTreeMap::new();
    for (id, attributes) in inputs {
        for a in attributes {
            let Some(policy) = table.0.get(&a.key) else {
                add(K::MissingPolicy, a.key, *id);
                continue;
            };
            let got = recorded.get(&(*id, a.key)).cloned().unwrap_or_default();
            match got.as_slice() {
                [] => add(K::MissingOutcome, a.key, *id),
                [o] => {
                    let want = expected_outcome(inputs, history, policy, a.key, *id);
                    if want.as_ref() != Some(*o) {
                        add(K::OutcomeContradictsPolicy, a.key, *id);
                    }
                    let copies = matches!(o.result, OutcomeResult::Kept | OutcomeResult::Copied);
                    let lands = matches!(
                        o.result,
                        OutcomeResult::Kept | OutcomeResult::Copied | OutcomeResult::Recomputed
                    );
                    for t in &o.to {
                        if lands {
                            carried.entry((*t, a.key)).or_default().push((o, copies));
                            if value(outputs, *t, a.key).is_none() {
                                add(K::MissingAttribute, a.key, *t);
                            }
                        }
                    }
                }
                _ => add(K::DuplicateOutcome, a.key, *id),
            }
        }
    }
    for (id, attributes) in outputs {
        for a in attributes {
            match carried.get(&(*id, a.key)) {
                None => {
                    // Present on an output that every outcome dropped or
                    // conflicted, or never carried at all.
                    let dropped = history.attributes.iter().any(|o| {
                        o.key == a.key
                            && matches!(o.result, OutcomeResult::Dropped | OutcomeResult::Conflict)
                            && history
                                .relations
                                .iter()
                                .any(|r| r.sources().contains(&o.from) && r.targets().contains(id))
                    });
                    add(
                        if dropped {
                            K::AttributeOnDropped
                        } else {
                            K::UntracedAttribute
                        },
                        a.key,
                        *id,
                    );
                }
                Some(sources) => {
                    for (o, copies) in sources {
                        if *copies && value(inputs, o.from, a.key) != Some(&a.value) {
                            add(K::ValueChanged, a.key, *id);
                        }
                    }
                }
            }
        }
    }
    issues.into_iter().collect()
}
