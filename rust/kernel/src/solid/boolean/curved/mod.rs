//! S9c.1: Booleans of prisms whose profiles hold arcs and circles, in any
//! relative position, where every pair of faces meets in lines, circles or
//! ellipses (REVIEW_NOTES.md, S9c).
//!
//! Each input is its exact model (`model.rs`). Vertices are where an edge
//! of one meets a face of the other (and where equal cylinders' ellipses
//! cross); edges are the inputs' edges split there and the sections of two
//! faces inside both (`meet.rs`, `graph.rs`); each face's pieces are traced
//! from them, classified by the other solid's exact membership and kept by
//! the operation's set function; the kept pieces are joined into the
//! result's faces and solids and named (`assemble.rs`). Full circles are
//! split at a rational seam, tried again at another when a meeting falls
//! on it. Coincident surfaces are not yet taken (`OutOfDomain`).
mod assemble;
mod graph;
mod meet;
mod model;
mod num;

use super::polyhedra::{Component, Polyhedron};
use crate::profile::boolean::Operand;
use crate::profile::{BoundaryKind, Segment};
use crate::solid::Construction;
use crate::{Error, Result};
use num_bigint::BigInt;
use num_rational::BigRational as R;

/// Whether S9c.1 takes the pair: two prisms, one of them with an arc.
pub(super) fn applies(poly: &Polyhedron) -> bool {
    let arcs = |s: &crate::Solid| match &s.construction {
        Construction::Prism(p) => p.boundaries().any(|b| match &b.kind {
            BoundaryKind::Circle { .. } => true,
            BoundaryKind::Path { segments, .. } => {
                segments.iter().any(|s| matches!(s, Segment::Arc { .. }))
            }
            BoundaryKind::Polygon(_) => false,
        }),
        _ => false,
    };
    let prism = |s: &crate::Solid| matches!(s.construction, Construction::Prism(_));
    prism(&poly.a) && prism(&poly.b) && (arcs(&poly.a) || arcs(&poly.b))
}

/// The result's components, a full circle's seam tried at several
/// rational points.
pub(super) fn build(poly: &Polyhedron) -> Result<Vec<Component>> {
    for (n, d) in [(2, 7), (3, 11), (5, 13), (7, 19)] {
        let seam = R::new(BigInt::from(n), BigInt::from(d));
        match attempt(poly, &seam) {
            Err(Error::ComputationLimit(m)) if m == graph::SEAM => continue,
            r => return r,
        }
    }
    Err(Error::Degenerate("a meeting at every seam tried"))
}

fn attempt(poly: &Polyhedron, seam: &R) -> Result<Vec<Component>> {
    let a = model::Prism::new(&poly.a, Operand::A, seam)?;
    let b = model::Prism::new(&poly.b, Operand::B, seam)?;
    let arr = graph::arrange([a, b], poly.op)?;
    assemble::assemble(&arr, poly.op)
}
