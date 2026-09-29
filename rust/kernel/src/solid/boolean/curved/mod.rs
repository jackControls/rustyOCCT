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
//! on it. Faces of both inputs on one surface hold each other's edges within
//! them, and their pieces facing one way join.
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
    // Each input's seam apart from the other's: one cylinder shared by both
    // would put both seams on one line.
    let seams = [(2, 7), (3, 11), (5, 13), (7, 19), (11, 23)];
    for k in 0..seams.len() - 1 {
        let r = |(n, d): (i64, i64)| R::new(BigInt::from(n), BigInt::from(d));
        match attempt(poly, &r(seams[k]), &r(seams[k + 1])) {
            Err(Error::ComputationLimit(m)) if m == graph::SEAM => continue,
            r => return r,
        }
    }
    Err(Error::Degenerate("a meeting at every seam tried"))
}

fn attempt(poly: &Polyhedron, seam_a: &R, seam_b: &R) -> Result<Vec<Component>> {
    let a = model::Prism::new(&poly.a, Operand::A, seam_a)?;
    let b = model::Prism::new(&poly.b, Operand::B, seam_b)?;
    let arr = graph::arrange([a, b], poly.op)?;
    assemble::assemble(&arr, poly.op)
}
