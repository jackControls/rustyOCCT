//! Enclosures through operations (Contract 5, T5, of IDENTITY_AND_HISTORY.md).
//!
//! Every operation rebuilds its outputs from the exact profile, so each
//! output entity first carries its own measured bound. An entity that
//! continues input entities (`Unchanged`, `Modified`, `Split`, `Merged`) is
//! then raised to the largest of their bounds, so no bound ever resets below
//! an input's, whatever its provenance. A generated entity (a cut face and
//! its edges) has no input of its kind and keeps its own measured bound. A
//! rigid transform rebuilds in a rounded frame: that rounding moves the body
//! but not the gaps between its curves and surfaces, which are measured on
//! the rebuilt geometry.
use super::Solid;
use crate::history::{History, Relation};
use crate::identity::EntityId;
use std::collections::BTreeMap;

pub(crate) fn carry(inputs: &[&Solid], history: &History, outputs: &mut [&mut Solid]) {
    let mut parents: BTreeMap<EntityId, Vec<EntityId>> = BTreeMap::new();
    for r in &history.relations {
        match r {
            Relation::Unchanged { id } => parents.entry(*id).or_default().push(*id),
            Relation::Modified { from, to } => parents.entry(*to).or_default().push(*from),
            Relation::Split { from, into } => {
                for to in into {
                    parents.entry(*to).or_default().push(*from);
                }
            }
            Relation::Merged { from, into } => parents
                .entry(*into)
                .or_default()
                .extend(from.iter().copied()),
            Relation::Generated { .. } | Relation::Deleted { .. } => {}
        }
    }
    let bound = |id: EntityId| {
        inputs
            .iter()
            .find_map(|s| s.topology.slot_of(id).and(s.topology.enclosure_bound(id)))
    };
    for output in outputs.iter_mut() {
        output.topology.raise_enclosures(|id| {
            parents
                .get(&id)?
                .iter()
                .filter_map(|p| bound(*p))
                .reduce(f64::max)
        });
    }
}
