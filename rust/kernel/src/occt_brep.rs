//! OCCT's `.brep` text format (T2 of `TOPOLOGY_MODEL.md`): a reader of every
//! record of versions 1 to 3, a converter of the solids it can represent into
//! the cell model, and a writer of cell topologies back to version 1 text.
//!
//! OCCT is one encoding of the kernel's model, not its definition. The
//! converter merges the two uses of a seam into periodic loops with winding
//! numbers, drops the vertices only seams and closed curves use, and turns
//! each solid into a solid region with its shells; the writer inserts seams
//! and their vertices back by rule. Anything the kernel cannot represent
//! (other surfaces and curves, degenerate edges, non-rigid locations, free
//! shapes, internal orientations, mesh-only faces) is reported by name and
//! counted, never approximated.
use std::fmt;

mod import;
pub mod read;
mod write;

pub use import::{import, Import, ImportedSolid, Rejected};
pub use read::{read, Document};
pub use write::write;

/// A malformed or inconsistent document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrepError {
    /// The text does not follow the format at `line` (1-based).
    Syntax { line: usize, what: &'static str },
    /// A reference to a table entry that does not exist.
    Reference { what: &'static str },
    /// A topology the writer cannot express in the format's structure.
    Unwritable { what: &'static str },
}

impl fmt::Display for BrepError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Syntax { line, what } => write!(f, "line {line}: expected {what}"),
            Self::Reference { what } => write!(f, "invalid reference: {what}"),
            Self::Unwritable { what } => write!(f, "cannot write {what}"),
        }
    }
}

impl std::error::Error for BrepError {}

/// A 3x4 affine matrix, row-major, as a location record stores it: a point
/// maps to `M[..3] p + M[3]` per row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform(pub [f64; 12]);

impl Transform {
    pub const IDENTITY: Self = Self([1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);

    fn at(&self, r: usize, c: usize) -> f64 {
        self.0[4 * r + c]
    }
    pub fn point(&self, p: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|r| {
            self.at(r, 0) * p[0] + self.at(r, 1) * p[1] + self.at(r, 2) * p[2] + self.at(r, 3)
        })
    }
    pub fn vector(&self, v: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|r| self.at(r, 0) * v[0] + self.at(r, 1) * v[1] + self.at(r, 2) * v[2])
    }
    /// `self * other`: `other` applied first, as `TopLoc_Location` composes.
    pub fn times(&self, other: &Self) -> Self {
        let mut out = [0.0; 12];
        for r in 0..3 {
            for c in 0..4 {
                let mut v = if c == 3 { self.at(r, 3) } else { 0.0 };
                for k in 0..3 {
                    v += self.at(r, k) * other.at(k, c);
                }
                out[4 * r + c] = v;
            }
        }
        Self(out)
    }
    fn inverse(&self) -> Result<Self, &'static str> {
        let m = |r: usize, c: usize| self.at(r, c);
        let det = m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1))
            - m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0))
            + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
        if det == 0.0 || !det.is_finite() {
            return Err("an invertible location");
        }
        // The adjugate's transpose over the determinant.
        let inv: [[f64; 3]; 3] = std::array::from_fn(|r| {
            std::array::from_fn(|c| {
                let (r1, r2) = ((c + 1) % 3, (c + 2) % 3);
                let (c1, c2) = ((r + 1) % 3, (r + 2) % 3);
                (m(r1, c1) * m(r2, c2) - m(r1, c2) * m(r2, c1)) / det
            })
        });
        let t: [f64; 3] = std::array::from_fn(|r| {
            -(inv[r][0] * m(0, 3) + inv[r][1] * m(1, 3) + inv[r][2] * m(2, 3))
        });
        Ok(Self([
            inv[0][0], inv[0][1], inv[0][2], t[0], inv[1][0], inv[1][1], inv[1][2], t[1],
            inv[2][0], inv[2][1], inv[2][2], t[2],
        ]))
    }
    /// `self` raised to an integer power; negative powers invert.
    pub fn powered(&self, power: i64) -> Result<Self, &'static str> {
        if power.unsigned_abs() > 1024 {
            return Err("a location power within 1024");
        }
        let base = if power < 0 { self.inverse()? } else { *self };
        let mut out = Self::IDENTITY;
        for _ in 0..power.unsigned_abs() {
            out = base.times(&out);
        }
        Ok(out)
    }
    /// A rotation and translation within `tol` in each matrix entry of
    /// `R^T R = I` with determinant 1: locations that only move shapes.
    pub fn is_rigid(&self, tol: f64) -> bool {
        let col = |c: usize| [self.at(0, c), self.at(1, c), self.at(2, c)];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let (x, y, z) = (col(0), col(1), col(2));
        let cross = [
            x[1] * y[2] - x[2] * y[1],
            x[2] * y[0] - x[0] * y[2],
            x[0] * y[1] - x[1] * y[0],
        ];
        (dot(x, x) - 1.0).abs() <= tol
            && (dot(y, y) - 1.0).abs() <= tol
            && (dot(z, z) - 1.0).abs() <= tol
            && dot(x, y).abs() <= tol
            && dot(y, z).abs() <= tol
            && dot(x, z).abs() <= tol
            && (dot(cross, z) - 1.0).abs() <= tol
            && self.0.iter().all(|v| v.is_finite())
    }
}
