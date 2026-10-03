//! S9e.3a: given results of spheres, cones and tori, with procedural
//! edges, of deeper chains, and against a sphere, cone or torus
//! (REVIEW_NOTES.md, S9e.3's refined decisions).
//!
//! A given result's construction is re-run whichever arrangement decided it
//! (`given.rs`), its leaves any exact model (a given result's too: the
//! construction tree evaluated level by level, at most `MAX_DEPTH`
//! Booleans deep). Every first-arrangement edge on a kept piece's boundary
//! is a model edge with its exact curve and places; the second arrangement
//! reaches a curve's own surface through the edge's faces' views (a
//! sphere's circle its ball, a cone's or a torus's section its funnel or
//! ring). One meeting is new (a plane's section of a cone against a plane,
//! `cone::ConeSec::meet_plane`); a meeting of two curved faces, a cone's
//! section against a curved face, or a torus's section other than a circle
//! of the torus against a curved face, met by another face, is S9e.3b's
//! (none where one of the edge's faces' surfaces is apart from it). A result edge along a whole given edge of a
//! procedural curve keeps that edge's stored curve (`stored`).
use super::graph::{Arr, CurveRef};
use super::model::{Crv, Prism};
use crate::solid::{Construction, Solid};
use crate::topology::Curve3;
use crate::{Error, Point3, Result};

/// The Booleans a given result's construction tree may hold.
pub(super) const MAX_DEPTH: usize = 3;

/// A given result's meeting of two curved faces (or a cone's or a torus's
/// plane section against a curved face) met by a face of the other input:
/// S9e.3b's.
pub(super) const S9E3B: &str =
    "a given result's meeting of two curved faces met by another face (S9e.3b)";

/// The Booleans of a solid's construction tree (its longest path).
pub(super) fn depth(s: &Solid) -> usize {
    match &s.construction {
        Construction::Polyhedron(p) => 1 + depth(&p.a).max(depth(&p.b)),
        Construction::Stack(_) => 1,
        _ => 0,
    }
}

/// A given result's construction within the depth limit.
pub(super) fn check_depth(s: &Solid) -> Result<()> {
    if depth(s) > MAX_DEPTH {
        return Err(Error::ComputationLimit(
            "a given result's construction deeper than three Booleans",
        ));
    }
    Ok(())
}

/// A model edge's own model of a curve's surface: of its faces' views, the
/// first whose model `holds` (a ball, a funnel, a ring).
pub(super) fn own_model(m: &Prism, ei: usize, holds: impl Fn(&Prism) -> bool) -> Option<&Prism> {
    m.edges[ei]
        .faces
        .iter()
        .map(|&f| m.view(f).0)
        .find(|v| holds(v))
}

/// The model whose data rounds a curve of an arrangement edge: for a
/// section, the view of its carrier operand's face; for a given edge, the
/// view of its face that `holds`; else the operand's model.
pub(super) fn curve_model(
    arr: &Arr,
    curve: CurveRef,
    carrier: usize,
    holds: impl Fn(&Prism) -> bool,
) -> &Prism {
    match curve {
        CurveRef::Section(si, _) => {
            let s = &arr.secs[si];
            arr.models[carrier]
                .view(if carrier == 0 { s.fa } else { s.fb })
                .0
        }
        CurveRef::Edge(o, ei) => {
            let m = &arr.models[o];
            if m.given.is_some() {
                own_model(m, ei, holds).unwrap_or(m)
            } else {
                m
            }
        }
    }
}

/// Whether a curve is one only S9e.3a's given models hold as an edge (a
/// section of the first arrangement).
pub(super) fn procedural(crv: &Crv) -> bool {
    matches!(
        crv,
        Crv::Meet(_)
            | Crv::Rise(_)
            | Crv::Toric(_)
            | Crv::Cone(_)
            | Crv::Torus(_)
            | Crv::WallMeet(_)
    )
}

/// Whether two faces' surfaces are apart where a plane meets a cylinder or a
/// sphere (exactly: a plane parallel to the cylinder's axis farther than
/// its radius, a plane farther from the sphere's centre than its radius);
/// none known otherwise.
pub(super) fn planes_apart(a: &Prism, fa: usize, b: &Prism, fb: usize) -> bool {
    use super::model::Surf;
    use super::num::{dot, sub};
    let ((va, ia), (vb, ib)) = (a.view(fa), b.view(fb));
    let (plane, (vq, iq)) = match (&va.faces[ia].surf, &vb.faces[ib].surf) {
        (Surf::Plane { p, m }, _) => ((p, m), (vb, ib)),
        (_, Surf::Plane { p, m }) => ((p, m), (va, ia)),
        _ => return false,
    };
    match &vq.faces[iq].surf {
        Surf::Cyl { c, r, .. } => {
            let [alpha, beta, mu, kappa] =
                super::meet::plane_on_cylinder(&vq.f, c, r, plane.0, plane.1);
            mu == crate::solid::split::zero() && &alpha * &alpha + &beta * &beta < &kappa * &kappa
        }
        Surf::Sphere { c, r } => {
            let off = dot(plane.1, &sub(c, plane.0));
            &off * &off > r * r * dot(plane.1, plane.1)
        }
        _ => false,
    }
}

/// Whether a result edge runs over the whole of the given result's edge
/// that model edge `ei` of operand `o` lies in: its parts are that result
/// edge's model edges (one first-arrangement edge each), every one of them,
/// none split by the second arrangement.
pub(super) fn whole(arr: &Arr, e: &super::assemble::REdge, o: usize, ei: usize) -> bool {
    let m = &arr.models[o];
    let Some(id) = m.edges[ei].id else {
        return false;
    };
    let mine: std::collections::BTreeSet<usize> = (0..m.edges.len())
        .filter(|&k| m.edges[k].id == Some(id))
        .collect();
    let mut seen = std::collections::BTreeSet::new();
    for &(gi, _) in &e.parts {
        match arr.edges[gi].curve {
            CurveRef::Edge(o2, k) if o2 == o && mine.contains(&k) => {
                if !seen.insert(k) {
                    return false;
                }
            }
            _ => return false,
        }
    }
    // Each of them one arrangement edge (not split by the other input).
    seen == mine
        && mine.iter().all(|&k| {
            arr.edges
                .iter()
                .filter(|x| x.curve == CurveRef::Edge(o, k))
                .count()
                == 1
        })
}

/// A stored curve run the other way.
pub(super) fn reversed(c: &Curve3) -> Option<Curve3> {
    let mut c = c.clone();
    match &mut c {
        Curve3::LineSegment { start, end } => std::mem::swap(start, end),
        Curve3::CircularArc {
            start_angle,
            sweep_angle,
            ..
        }
        | Curve3::EllipseArc {
            start_angle,
            sweep_angle,
            ..
        } => {
            *start_angle += *sweep_angle;
            *sweep_angle = -*sweep_angle;
        }
        Curve3::HyperbolaArc { start, sweep, .. } | Curve3::ParabolaArc { start, sweep, .. } => {
            *start += *sweep;
            *sweep = -*sweep;
        }
        Curve3::Section(s) => {
            s.start += s.sweep;
            s.sweep = -s.sweep;
        }
        Curve3::Meet(s) => {
            s.start += s.sweep;
            s.sweep = -s.sweep;
        }
        Curve3::Rise(s) => {
            s.start += s.sweep;
            s.sweep = -s.sweep;
        }
        Curve3::Toric(s) => {
            s.start += s.sweep;
            s.sweep = -s.sweep;
        }
        Curve3::WallMeet(s) => {
            s.start += s.sweep;
            s.sweep = -s.sweep;
        }
        Curve3::Circle { .. } | Curve3::BSpline(_) => return None,
    }
    Some(c)
}

/// A given edge's stored curve for a result edge running over the whole of
/// it from `start` to `end` (rounded points): as stored, or reversed, the
/// one whose ends lie nearer.
pub(super) fn stored(c: &Curve3, ends: Option<(Point3, Point3)>) -> Option<Curve3> {
    let Some((start, end)) = ends else {
        // A ring: as stored.
        return Some(c.clone());
    };
    let d = |a: Point3, b: Point3| (a - b).length();
    let (p0, p1) = (c.point(0.0), c.point(1.0));
    if d(p0, start) + d(p1, end) <= d(p1, start) + d(p0, end) {
        Some(c.clone())
    } else {
        reversed(c)
    }
}
