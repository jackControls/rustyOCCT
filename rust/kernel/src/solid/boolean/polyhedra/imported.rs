//! S9e.4b.2: imported polyhedra other than prisms (REVIEW_NOTES.md, "S9e.4b.2
//! refined").
//!
//! A body read from a file whose faces are planes and edges lines, and which
//! is no S9e.4a prism, is decided on its stored vertices: each its binary64
//! point as a rational, each edge the segment between its two stored
//! vertices, each face the polygon of its stored vertices, cut into exactly
//! planar triangles in its projection (S9b.2's stored model of a Boolean's
//! result, `stored_model`), the solid's side by an exact ray's parity. Its
//! stored planes only orient its faces: their common points are not taken
//! (a vertex of four planes has none once they are rounded, and two bodies
//! sharing a face would meet along slivers where their stored points agree).
//! A face whose stored vertices are coplanar exactly is that plane, joined
//! with any face on it; a face they fold is `Degenerate`; against a body
//! with curved faces or edges (S9e.4b.4c.1) its triangles are a leaf of the
//! curved engine (`curved::meshes`). Its shells may hold cavities.
use super::{approx3, parity, stored_model, vq, Inside, Polyhedron};
use crate::profile::boolean::Operand;
use crate::solid::{Construction, Solid};
use crate::topology::{Curve3, Surface, Topology};
use crate::{Error, Location, Point3, Result};

/// Whether a solid is an imported polyhedron decided on its stored vertices.
pub(crate) fn is_imported(s: &Solid) -> bool {
    matches!(&s.construction, Construction::Imported(i) if i.polyhedron())
}

/// Whether every face is a plane and every edge a line.
pub(crate) fn planar(t: &Topology) -> bool {
    t.faces()
        .iter()
        .all(|f| matches!(f.surface, Surface::Plane(_)))
        && t.edges()
            .iter()
            .all(|e| matches!(e.curve, Curve3::LineSegment { .. }))
}

/// Whether a Boolean has an imported polyhedron among its inputs against a
/// body of plane faces and line edges (then decided on the stored models);
/// against curved faces or edges (S9e.4b.4c.1) the curved engine's, the
/// polyhedron its stored triangles (`curved::meshes`).
pub(super) fn involved(poly: &Polyhedron) -> Result<bool> {
    let (a, b) = (is_imported(&poly.a), is_imported(&poly.b));
    if (a && !planar(&poly.b.topology)) || (b && !planar(&poly.a.topology)) {
        return Ok(false);
    }
    Ok(a || b)
}

pub(super) fn folded() -> Error {
    Error::Degenerate("an imported face folded by its stored vertices")
}

/// S9e.4b.4c.1: an imported polyhedron's stored model as the curved
/// engine's leaf: each stored face's triangles (the face's index, its
/// corners counter-clockwise about its outward normal), every corner a
/// stored vertex exactly; `ComputationLimit` where a face's triangles are
/// not of its own vertices (its trapezoids zipped).
pub(crate) fn triangles(
    solid: &Solid,
) -> Result<Vec<(usize, [[num_rational::BigRational; 3]; 3])>> {
    use crate::topology::{FaceId, Slot};
    use std::collections::{BTreeMap, BTreeSet};
    let model = stored_model(solid, Operand::A)?;
    let t = &solid.topology;
    let index: BTreeMap<crate::identity::EntityId, usize> = (0..t.faces().len())
        .filter_map(|i| t.id_of(Slot::Face(FaceId(i))).map(|id| (id, i)))
        .collect();
    let stored: BTreeSet<super::V> = t
        .vertices()
        .iter()
        .map(|v| vq(v.position - Point3::ORIGIN))
        .collect();
    let mut out = Vec::new();
    for face in &model.faces {
        let fi = *index
            .get(&face.id)
            .ok_or(Error::InvalidTopology("an imported polyhedron's face"))?;
        for piece in &face.pieces {
            let [a, b, c] = <[super::V; 3]>::try_from(piece.clone())
                .map_err(|_| Error::InvalidTopology("an imported polyhedron's triangle"))?;
            if ![&a, &b, &c].iter().all(|p| stored.contains(*p)) {
                return Err(Error::ComputationLimit(
                    "an imported polyhedron's face not cut into triangles of its own vertices",
                ));
            }
            out.push((fi, [a, b, c]));
        }
    }
    Ok(out)
}

/// The stored model of an imported polyhedron builds (its faces' triangles
/// facing their way): checked once on import.
pub(crate) fn check(solid: &Solid) -> Result<()> {
    stored_model(solid, Operand::A).map(|_| ())
}

/// A point's location against an imported polyhedron: outside its bounds
/// widened by twice the resolution, `Boundary` within the resolution of a
/// face's triangle (in binary64), else by the exact parity of a ray.
pub(crate) fn classify(solid: &Solid, point: Point3) -> Result<Location> {
    let tol = solid.resolution().linear();
    let p = point.to_array();
    let (lo, hi) = (solid.bounds.min.to_array(), solid.bounds.max.to_array());
    if (0..3).any(|i| p[i] < lo[i] - 2.0 * tol || p[i] > hi[i] + 2.0 * tol) {
        return Ok(Location::Outside);
    }
    let model = stored_model(solid, Operand::A)?;
    for face in &model.faces {
        for tri in &face.pieces {
            let c: Vec<[f64; 3]> = tri.iter().map(approx3).collect();
            if distance(p, [c[0], c[1], c[2]]) <= tol {
                return Ok(Location::Boundary);
            }
        }
    }
    let Inside::Mesh(triangles) = &model.inside else {
        return Err(Error::InvalidTopology("an imported polyhedron's model"));
    };
    match parity(triangles, &vq(point - Point3::ORIGIN)) {
        Some(true) => Ok(Location::Inside),
        Some(false) => Ok(Location::Outside),
        None => Err(Error::ComputationLimit(
            "a point every ray tried meets an imported polyhedron's edge",
        )),
    }
}

/// A point's distance from a triangle.
fn distance(p: [f64; 3], t: [[f64; 3]; 3]) -> f64 {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let segment = |a: [f64; 3], b: [f64; 3]| {
        let (ab, ap) = (sub(b, a), sub(p, a));
        let len2 = dot(ab, ab);
        let s = if len2 > 0.0 {
            (dot(ap, ab) / len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let d = [ap[0] - s * ab[0], ap[1] - s * ab[1], ap[2] - s * ab[2]];
        dot(d, d).sqrt()
    };
    let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
    let nn = dot(n, n);
    if nn > 0.0 {
        let inside = (0..3).all(|i| dot(cross(sub(t[(i + 1) % 3], t[i]), sub(p, t[i])), n) >= 0.0);
        if inside {
            return dot(sub(p, t[0]), n).abs() / nn.sqrt();
        }
    }
    (0..3)
        .map(|i| segment(t[i], t[(i + 1) % 3]))
        .fold(f64::INFINITY, f64::min)
}
