//! S9e.4b.3b: the kernel's own plane pieces, a prism's oblique piece
//! (`Clipped`, S8a.2) and a cone's, a zone's or a torus's piece (`Half`,
//! S8c.2, S8d.1, S8d.2), against curved faces (REVIEW_NOTES.md, "S9e.4b.3b
//! refined").
//!
//! A split piece is its primitive common the half-space of its plane, on
//! S9e.4b.3a's model (`pieces.rs`): the primitive the split's own (the
//! prism on its frame between its heights, the cone or whole torus on its
//! frame; a zone's or cap's whole sphere with its ends' parallels' planes)
//! and the hull the plane the piece was built on, in the solid's frame
//! through its exact axes (`F(u, v, w) = a u + b v + c w + d` is `m . (X -
//! o) + d` with `m` the plane's coefficients through the frame's inverse),
//! on the piece's side, bounded by a cube about the primitive.
//! Their arrangement's assembly (one solid, or several where a prism's
//! profile is not convex across the plane) is matched to the piece's stored
//! topology by S9e.2's geometric match, so the given model carries the
//! piece's ids and a Boolean's history is over them directly.
use super::assemble::Made;
use super::graph::Arr;
use super::model::*;
use super::num::*;
use super::pieces::{common, cube, hull_of, hull_operation};
use crate::identity::OperationId;
use crate::profile::boolean::{Op2, Operand};
use crate::solid::boolean::polyhedra::Component;
use crate::solid::split::{q, rational_f64, zero, Primitive};
use crate::solid::{Construction, Solid};
use crate::topology::Surface;
use crate::{Error, Frame3, Point3, Result, Vec3};
use num_rational::BigRational as R;
use std::f64::consts::{FRAC_PI_2, TAU};

/// Whether a solid is a split piece (S8a.2's `Clipped`, S8c.2's, S8d.1's
/// and S8d.2's `Half`).
pub(crate) fn is_split(s: &Solid) -> bool {
    matches!(
        s.construction,
        Construction::Clipped(_) | Construction::Half(_)
    )
}

/// The primitive's operation: the piece's, its bits turned another way
/// than an import's construction's.
fn primitive_operation(s: &Solid) -> OperationId {
    OperationId(s.operation.0.rotate_left(17) ^ 0x3b3b_9e4b)
}

/// A plane in a frame's local coordinates (`a u + b v + c w + d`) and the
/// sign of its function in the material.
type Local = ([R; 4], i8);

/// The split's primitive on the piece's frame and the planes bounding the
/// piece in that frame: the plane the split was built on and, for a zone, its
/// ends' planes (its primitive the whole sphere, as S9e.4b.3a's pieces of a
/// sphere: a plane through the axis then meets the sphere through its
/// stored poles as an imported half's does).
fn primitive(s: &Solid) -> Result<(Solid, Vec<Local>)> {
    let op = primitive_operation(s);
    let tolerance = s.resolution();
    match &s.construction {
        Construction::Clipped(c) => {
            let (profile, plane, sign) = c.parts();
            if super::profile_splines(profile) {
                return Err(Error::OutOfDomain(
                    "a spline prism's split piece against curved faces (S9f)",
                ));
            }
            let prism = Solid::build(op, profile.clone(), s.frame, s.start, s.end)?;
            Ok((prism, vec![(plane.clone(), sign)]))
        }
        Construction::Half(h) => {
            let (primitive, _, index) = h.parts();
            let plane = h.built_plane();
            let side = if index == 0 { -1 } else { 1 };
            let mut planes = vec![(plane.clone(), side)];
            let solid = match *primitive {
                Primitive::Cone {
                    bottom,
                    top,
                    height,
                } => Solid::build_cone(op, s.frame, bottom, top, height, tolerance)?,
                Primitive::Zone { radius, low, high } => {
                    // The ends' parallels' planes `w = start` (the material
                    // above) and `w = end` (below), a pole none.
                    let level = |w: f64| [zero(), zero(), int(1), -q(w)];
                    if low > -FRAC_PI_2 {
                        planes.push((level(s.start), 1));
                    }
                    if high < FRAC_PI_2 {
                        planes.push((level(s.end), -1));
                    }
                    Solid::build_sphere(op, s.frame, radius, -FRAC_PI_2, FRAC_PI_2, tolerance)?
                }
                Primitive::Torus { major, minor } => {
                    // S8d.3's spiric pieces are matched by no rule of
                    // S9e.2's: only S8d.1's bands.
                    if plane[0] != zero() || plane[1] != zero() {
                        return Err(Error::OutOfDomain(
                            "a torus's piece by a plane oblique to its axis against curved faces \
                             (refused, S9e.4b.3b)",
                        ));
                    }
                    Solid::build_torus(op, s.frame, major, minor, 0.0, TAU, TAU, tolerance)?
                }
            };
            Ok((solid, planes))
        }
        _ => unreachable!("a split piece"),
    }
}

/// The plane `a u + b v + c w + d` of a frame's local coordinates in the
/// world, exactly: a point on it and its normal turned out of the material
/// (where the function's sign is `-sign`).
fn world_plane(frame: &Frame3, plane: &[R; 4], sign: i8) -> Result<(V, V)> {
    let f = Affine::new(frame)?;
    let [a, b, c, d] = plane;
    // F(X) = l . inv (X - o) + d = m . (X - o) + d, m = inv^T l.
    let unit = |k: usize| {
        let mut e = [zero(), zero(), zero()];
        e[k] = int(1);
        f.local_dir(&e)
    };
    let (e0, e1, e2) = (unit(0), unit(1), unit(2));
    let m: V = [
        a * &e0[0] + b * &e0[1] + c * &e0[2],
        a * &e1[0] + b * &e1[1] + c * &e1[2],
        a * &e2[0] + b * &e2[1] + c * &e2[2],
    ];
    let m2 = dot(&m, &m);
    if m2 == zero() {
        return Err(Error::Degenerate("a split's plane"));
    }
    // The point o - d m / |m|^2.
    let p = sub(&f.o, &scale(&m, &(d / &m2)));
    // Material where F has the sign `sign`: outward along `-sign F`'s
    // gradient.
    Ok((p, if sign < 0 { m } else { neg(&m) }))
}

/// A plane's frame for its hull face's parameters: the piece's stored
/// face on it (its stored frame, bit for bit), else one at the rounded
/// point along the rounded normal.
fn face_frame(s: &Solid, point: &V, normal: &V) -> Result<Frame3> {
    let n = normal.clone().map(|x| rational_f64(&x));
    let n = Vec3::new(n[0], n[1], n[2]);
    let o = point.clone().map(|x| rational_f64(&x));
    let o = Point3::new(o[0], o[1], o[2]);
    let tol = s.resolution().linear();
    for face in s.topology.faces() {
        let Surface::Plane(f) = &face.surface else {
            continue;
        };
        if f.normal().cross(n).length() <= 1e-9 * n.length()
            && (f.origin() - o).dot(n).abs() <= tol * n.length()
        {
            return Ok(*f);
        }
    }
    let hint = if n.cross(Vec3::new(1.0, 0.0, 0.0)).length() > 0.5 * n.length() {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    Frame3::new(o, n, hint, s.resolution())
}

/// The split piece's primitive common its planes' half-spaces, arranged
/// and assembled.
fn arranged(s: &Solid) -> Result<(Arr, Vec<(Component, Made)>)> {
    let (p, planes) = primitive(s)?;
    let mut exact = Vec::new();
    for (plane, sign) in &planes {
        let (point, normal) = world_plane(&s.frame, plane, *sign)?;
        let frame = face_frame(s, &point, &normal)?;
        exact.push(((point, normal), frame));
    }
    let (centre, half) = cube(&p);
    let hull = hull_of(
        &exact,
        centre,
        half,
        Operand::B,
        hull_operation(&p),
        s.resolution(),
    )?;
    common(&p, &hull)
}

thread_local! {
    /// The last split pieces' models, by their content (the piece's
    /// construction, frame, heights and stored topology's `Debug` text).
    static ARRANGED: std::cell::RefCell<Vec<(String, Result<Prism>)>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// A split piece's model (S9e.4b.3b): the given model of its primitive
/// common its planes' half-spaces, matched to its stored topology.
pub(super) fn model(s: &Solid, op: Operand) -> Result<Prism> {
    let key = format!(
        "{:?}\n{:?}\n{:?} {:?}\n{:?}\n{op:?}",
        s.construction, s.frame, s.start, s.end, s.topology
    );
    if let Some(hit) = ARRANGED.with(|k| {
        k.borrow()
            .iter()
            .find(|(known, _)| *known == key)
            .map(|(_, m)| m.clone())
    }) {
        return hit;
    }
    let built = (|| {
        let (arr, out) = arranged(s)?;
        super::given::built(s, op, arr, out, Op2::Common, None)
    })();
    ARRANGED.with(|k| {
        let mut k = k.borrow_mut();
        if k.len() >= 2 {
            k.drain(..1);
        }
        k.push((key, built.clone()));
    });
    built
}

/// A split zone's sphere, its centre and radius (S9e.4b.3c's pieces of one
/// sphere).
pub(super) fn ball_of(s: &Solid) -> Option<(Point3, f64)> {
    match &s.construction {
        Construction::Half(h) => match h.parts().0 {
            Primitive::Zone { radius, .. } => Some((s.frame.origin(), *radius)),
            _ => None,
        },
        _ => None,
    }
}
