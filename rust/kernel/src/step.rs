//! STEP (ISO 10303-21 exchange structures of AP203, AP214 and AP242): a
//! Part 21 reader and, through OCCT's shape structure, an importer of the
//! solids and surface models the kernel can represent (the STEP import track
//! of `REVIEW_NOTES.md`).
use std::fmt;

mod import;
pub mod part21;

pub use import::{import, Item, StepBody, StepImport};
pub use part21::{read, Exchange, Instance, Parameter, Record};

/// A malformed exchange structure, or one of a schema the importer does not
/// read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StepError {
    /// The text does not follow Part 21 at `line` (1-based).
    Syntax { line: usize, what: &'static str },
    /// Instance `entity` refers to `target`, which no data section defines.
    Reference { entity: u64, target: u64 },
    /// `FILE_SCHEMA` names no schema of AP203, AP214 or AP242.
    Schema(String),
}

impl fmt::Display for StepError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Syntax { line, what } => write!(f, "line {line}: expected {what}"),
            Self::Reference { entity, target } => {
                write!(f, "#{entity} refers to #{target}, which is not defined")
            }
            Self::Schema(name) => write!(f, "unsupported schema: {name}"),
        }
    }
}

impl std::error::Error for StepError {}
