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
mod algebraic;
mod assemble;
mod cone;
mod cones;
mod graph;
mod meet;
mod model;
mod num;
mod procedural;
mod sphere;
mod spheres;
mod torus;
mod torus_curved;
mod torus_segment;
mod turned;

use super::polyhedra::{Component, Polyhedron};
use crate::profile::boolean::Operand;
use crate::profile::{BoundaryKind, Segment};
use crate::solid::Construction;
use crate::{Error, Result};
use num_bigint::BigInt;
use num_rational::BigRational as R;

/// Whether a prism's profile holds an arc or a circle.
fn applies_arcs(s: &crate::Solid) -> bool {
    match &s.construction {
        Construction::Prism(p) => p.boundaries().any(|b| match &b.kind {
            BoundaryKind::Circle { .. } => true,
            BoundaryKind::Path { segments, .. } => {
                segments.iter().any(|s| matches!(s, Segment::Arc { .. }))
            }
            BoundaryKind::Polygon(_) => false,
        }),
        _ => false,
    }
}

/// Whether S9c.1 takes the pair: two prisms, one of them with an arc.
pub(super) fn applies(poly: &Polyhedron) -> bool {
    let arcs = applies_arcs;
    let prism = |s: &crate::Solid| matches!(s.construction, Construction::Prism(_));
    let sphere = |s: &crate::Solid| matches!(s.construction, Construction::Sphere { .. });
    let cone = |s: &crate::Solid| matches!(s.construction, Construction::Cone { .. });
    let torus = |s: &crate::Solid| matches!(s.construction, Construction::Torus { .. });
    let quadric = |s: &crate::Solid| prism(s) || sphere(s) || cone(s);
    // S9d.4: a torus against a prism, a sphere, a cone or a torus (parts
    // against curved faces refused in `build`).
    let any = |s: &crate::Solid| quadric(s) || torus(s);
    if (torus(&poly.a) && any(&poly.b)) || (any(&poly.a) && torus(&poly.b)) {
        return true;
    }
    // S9d.3: a cone against a prism, a sphere or a cone (S9d.3b's refused
    // in `build`).
    ((cone(&poly.a) && quadric(&poly.b)) || (quadric(&poly.a) && cone(&poly.b)))
        // S9d.1: a sphere against a prism.
        || (prism(&poly.a) && prism(&poly.b) && (arcs(&poly.a) || arcs(&poly.b)))
        || (sphere(&poly.a) && prism(&poly.b))
        || (prism(&poly.a) && sphere(&poly.b))
        // S9d.2: two spheres.
        || (sphere(&poly.a) && sphere(&poly.b))
}

/// An input's exact model: a prism's, a sphere's (S9d.1) or a cone's
/// (S9d.3a).
fn model_of(s: &crate::Solid, op: Operand, seam: &R) -> Result<model::Prism> {
    match &s.construction {
        Construction::Sphere { .. } => sphere::model(s, op, seam),
        Construction::Cone { .. } => cone::model(s, op, seam),
        Construction::Torus { .. } => torus::model(s, op, seam),
        _ => model::Prism::new(s, op, seam),
    }
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
    // S9d.4a and S9d.4b.1 take a torus (a part too) against a polyhedral
    // prism; S9d.4b.2a a whole torus against a prism with arcs, a sphere or
    // a cone, S9d.4b.2b against another whole torus.
    let torus = |s: &crate::Solid| matches!(s.construction, Construction::Torus { .. });
    let whole = |s: &crate::Solid| match &s.construction {
        Construction::Torus {
            low, high, angle, ..
        } => high - low == std::f64::consts::TAU && *angle == std::f64::consts::TAU,
        _ => false,
    };
    let polyhedral =
        |s: &crate::Solid| matches!(s.construction, Construction::Prism(_)) && !applies_arcs(s);
    if (torus(&poly.a) && !polyhedral(&poly.b) && !whole(&poly.a))
        || (torus(&poly.b) && !polyhedral(&poly.a) && !whole(&poly.b))
    {
        return Err(Error::OutOfDomain(
            "a torus segment or wedge against a curved face (S9d.4b)",
        ));
    }
    let a = model_of(&poly.a, Operand::A, seam_a)?;
    let b = model_of(&poly.b, Operand::B, seam_b)?;
    let arr = graph::arrange([a, b], poly.op)?;
    assemble::assemble(&arr, poly.op)
}
