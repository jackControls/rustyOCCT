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
//! them, and their pieces facing one way join. S9f.1: a spline prism's walls
//! against a prism of lines (`spline_walls.rs`); S9f.2a: against arc,
//! circle and spline walls on an exactly parallel axis
//! (`spline_parallel.rs`); S9f.2b: against cylinders on crossing axes
//! (`spline_crossing.rs`); S9f.3a: against spheres, caps and zones
//! (`spline_sphere.rs`); S9f.3b: against cones and frustums
//! (`spline_cone.rs`).
mod algebraic;
mod assemble;
mod chain;
mod cone;
mod cones;
mod cones_loops;
mod given;
mod graph;
pub(crate) mod matched;
mod meet;
mod model;
mod num;
mod procedural;
mod snapped;
mod sphere;
mod spheres;
mod spheres_turned;
mod spline_cone;
mod spline_crossing;
mod spline_parallel;
mod spline_sphere;
mod spline_walls;
mod torus;
mod torus_curved;
mod torus_parts;
mod torus_segment;
mod triple;
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

/// Whether a profile holds a spline segment.
fn profile_splines(p: &crate::Profile) -> bool {
    p.boundaries().any(|b| match &b.kind {
        BoundaryKind::Path { segments, .. } => {
            segments.iter().any(|s| matches!(s, Segment::Spline(_)))
        }
        _ => false,
    })
}

/// Whether a prism's profile holds a spline segment (S9f.1).
pub(super) fn applies_splines(s: &crate::Solid) -> bool {
    matches!(&s.construction, Construction::Prism(p) if profile_splines(p))
}

/// Whether a solid is a prism of a line profile (S9b.1's polyhedra).
fn line_prism(s: &crate::Solid) -> bool {
    match &s.construction {
        Construction::Prism(p) => p.boundaries().all(|b| match &b.kind {
            BoundaryKind::Polygon(_) => true,
            BoundaryKind::Path { segments, .. } => {
                segments.iter().all(|s| matches!(s, Segment::Line))
            }
            BoundaryKind::Circle { .. } => false,
        }),
        _ => false,
    }
}

/// Whether a Boolean's result (a stack among them) was made from a spline
/// prism, at any depth.
fn spline_leaves(s: &crate::Solid) -> bool {
    match &s.construction {
        Construction::Polyhedron(p) => [&p.a, &p.b]
            .iter()
            .any(|x| applies_splines(x) || spline_leaves(x)),
        Construction::Stack(st) => profile_splines(&st.a) || profile_splines(&st.b),
        _ => false,
    }
}

/// Whether two solids' stored normals are exactly parallel (as rationals).
fn parallel_axes(x: &crate::Solid, y: &crate::Solid) -> bool {
    let v = |s: &crate::Solid| s.frame.normal().to_array().map(crate::solid::split::q);
    let (a, b) = (v(x), v(y));
    let zero = crate::solid::split::zero();
    &a[1] * &b[2] - &a[2] * &b[1] == zero
        && &a[2] * &b[0] - &a[0] * &b[2] == zero
        && &a[0] * &b[1] - &a[1] * &b[0] == zero
}

/// S9f.1 takes a spline prism against a prism of lines in any position,
/// S9f.2a against a prism with arcs, circles or splines on an exactly
/// parallel axis, S9f.2b against a prism with arcs or circles (no spline) on
/// a crossing axis, S9f.3a against a sphere, a cap or a zone, S9f.3b
/// against a cone or a frustum; the others are refused (REVIEW_NOTES.md,
/// "S9f refined", "S9f.3 refined" and "S9f.3b refined").
fn spline_pairs(poly: &Polyhedron) -> Result<()> {
    for (x, y) in [(&poly.a, &poly.b), (&poly.b, &poly.a)] {
        if spline_leaves(x) {
            return Err(Error::OutOfDomain(
                "a Boolean's result with spline walls given to another Boolean (S9f)",
            ));
        }
        if !applies_splines(x) || line_prism(y) {
            continue;
        }
        if matches!(y.construction, Construction::Prism(_)) && parallel_axes(x, y) {
            continue;
        }
        // S9f.2b: against a prism of lines, arcs and circles on a crossing
        // axis (spline walls against cylinder walls).
        if matches!(&y.construction, Construction::Prism(p) if !profile_splines(p)) {
            continue;
        }
        // S9f.3a: against a sphere, a cap or a zone (spline walls against
        // the sphere); S9f.3b: against a cone or a frustum.
        if matches!(
            y.construction,
            Construction::Sphere { .. } | Construction::Cone { .. }
        ) {
            continue;
        }
        return Err(Error::OutOfDomain(match &y.construction {
            Construction::Prism(_) => {
                "spline walls against spline walls on crossing axes (refused, S9f)"
            }
            Construction::Torus { .. } => "a spline prism against a torus (refused, S9f)",
            _ => "a spline prism against a solid other than a prism in any position (S9f)",
        }));
    }
    Ok(())
}

/// Whether S9c.1 takes the pair: two prisms, one of them with an arc.
pub(super) fn applies(poly: &Polyhedron) -> bool {
    let arcs = applies_arcs;
    let prism = |s: &crate::Solid| matches!(s.construction, Construction::Prism(_));
    let sphere = |s: &crate::Solid| matches!(s.construction, Construction::Sphere { .. });
    let cone = |s: &crate::Solid| matches!(s.construction, Construction::Cone { .. });
    let torus = |s: &crate::Solid| matches!(s.construction, Construction::Torus { .. });
    let quadric = |s: &crate::Solid| prism(s) || sphere(s) || cone(s);
    // S9d.4: a torus (a part too, S9d.4c) against a prism, a sphere, a cone
    // or a torus.
    let any = |s: &crate::Solid| quadric(s) || torus(s);
    if (torus(&poly.a) && any(&poly.b)) || (any(&poly.a) && torus(&poly.b)) {
        return true;
    }
    // S9e.1: a Boolean's result from this arrangement given to another
    // Boolean, with a prism or another such result (with a sphere, cone or
    // torus refused in `build`: S9e.3's); S9e.2: a stack or an S9b.1 result
    // where either input has an arc or a curved face, and one solid of
    // several of these.
    if given::applies(&poly.a, &poly.b) || given::applies(&poly.b, &poly.a) {
        return true;
    }
    // S9f.1: a spline prism against a prism (against a solid other than a
    // prism of lines refused in `build`).
    if prism(&poly.a) && prism(&poly.b) && (applies_splines(&poly.a) || applies_splines(&poly.b)) {
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
        // S9e.1: a Boolean's result, on its construction's exact model
        // (S9e.2: a stack too).
        Construction::Polyhedron(_) | Construction::Stack(_) => given::model(s, op),
        Construction::Sphere { .. } => sphere::model(s, op, seam),
        Construction::Cone { .. } => cone::model(s, op, seam),
        Construction::Torus { .. } => torus::model(s, op, seam),
        _ => model::Prism::new(s, op, seam),
    }
}

/// The result's components, a full circle's seam tried at several
/// rational points.
/// Each input's seam apart from the other's: one cylinder shared by both
/// would put both seams on one line.
const SEAMS: [(i64, i64); 5] = [(2, 7), (3, 11), (5, 13), (7, 19), (11, 23)];

pub(super) fn build(poly: &Polyhedron) -> Result<Vec<Component>> {
    // S9f.1: a spline prism against a prism of lines only.
    spline_pairs(poly)?;
    // S9e.1 and S9e.2 take a given result with a prism or another given
    // result, S9e.3a with a sphere, cone or torus too.
    let piece = |s: &crate::Solid| {
        matches!(
            s.construction,
            Construction::Clipped(_) | Construction::Half(_)
        )
    };
    for (x, y) in [(&poly.a, &poly.b), (&poly.b, &poly.a)] {
        if given::applies(x, y) && piece(y) {
            return Err(Error::OutOfDomain(
                "a Boolean's result given to another Boolean with a plane's piece (S9e.4)",
            ));
        }
    }
    let seams = SEAMS;
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
    // prism; S9d.4b.2 a whole torus against a prism with arcs, a sphere, a
    // cone or another whole torus; S9d.4c a part against those too, and a
    // sphere's cap or zone against a torus.
    let mut arr = shared(poly, seam_a, seam_b)?;
    arr.for_op(poly.op);
    assemble::assemble(&arr, poly.op)
}

/// The arrangements kept for the next operations on the same inputs.
const KEPT: usize = 2;

thread_local! {
    /// The last arrangements of two inputs, by their content: the inputs'
    /// and the seams' `Debug` text (which tells every binary64 number
    /// apart). Fuse, cut and common of one pair share one (its pieces'
    /// sides decided, the operation's keeping not), bit for bit what each
    /// would build.
    static ARRANGED: std::cell::RefCell<Vec<(String, Result<graph::Arr>)>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// The models' arrangement every operation on the pair shares, kept.
fn shared(poly: &Polyhedron, seam_a: &R, seam_b: &R) -> Result<graph::Arr> {
    let key = format!("{:?}\n{:?}\n{seam_a}\n{seam_b}", poly.a, poly.b);
    if let Some(hit) = ARRANGED.with(|k| {
        k.borrow()
            .iter()
            .find(|(known, _)| *known == key)
            .map(|(_, arr)| arr.clone())
    }) {
        return hit;
    }
    let arr = (|| {
        let a = model_of(&poly.a, Operand::A, seam_a)?;
        let b = model_of(&poly.b, Operand::B, seam_b)?;
        let mut arr = graph::arrange_shared([a, b])?;
        // S9e.2: a given solid of several kept alone.
        matched::keep_solid(&mut arr)?;
        Ok(arr)
    })();
    ARRANGED.with(|k| {
        let mut k = k.borrow_mut();
        if k.len() >= KEPT {
            k.drain(..1);
        }
        k.push((key, arr.clone()));
    });
    arr
}
