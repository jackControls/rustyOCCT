//! S9e.4b.1: an imported prism's arcs' ends taken onto their circles
//! (REVIEW_NOTES.md, "S9e.4b refined").
//!
//! S9e.4a's construction of an imported prism rounds its stored vertices'
//! local coordinates once in its bottom cap's stored frame, and its arcs'
//! centres likewise: a vertex of a turned or tilted cap, or at an angle
//! whose cosine is irrational, then lies off its arc's circle by the
//! rounding. S9c's exact model takes an arc between its profile points
//! exactly on its circle, so the model of such a profile (flagged by the
//! import, `Profile::rounded_arcs`) takes each arc's end onto the circle:
//! the circle's rational point at the binary64 rounding of the end's
//! half-angle tangent, computed exactly. A joint of a line and an arc takes
//! the arc's point (the line's end with it); a joint of two arcs of one
//! circle that circle's point, the same for both; a joint of two lines its
//! rounded point; an end already on its circle stays. The ends move by at
//! most their distance from the circle plus `r 2^-52`, within the
//! resolution (the profile's validation holds every arc's ends within the
//! resolution of its circle). A joint of arcs of two different circles off
//! either is S9e.4b.4's: their common point is a quadratic surd, which a
//! profile segment's rational ends do not hold.
use super::model::P2;
use super::num::int;
use crate::profile::Segment;
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Point2, Result};
use num_rational::BigRational as R;

/// Whether `p` lies on the circle `|p - c| = r` exactly.
fn on(c: &P2, r: &R, p: &P2) -> bool {
    let (dx, dy) = (&p[0] - &c[0], &p[1] - &c[1]);
    &dx * &dx + &dy * &dy == r * r
}

/// The circle's rational point nearest `p`'s direction from its centre:
/// `c + r ((1 - s^2), 2 s) / (1 + s^2)` where `p`'s `dx` is not negative,
/// `c + r (-(1 - s^2), 2 s) / (1 + s^2)` where it is, `s` the binary64
/// rounding of `dy / (r + |dx|)` (the half-angle tangent from the nearer
/// end of the circle's `x` diameter, `|s| <= 1` within rounding); `p`
/// itself where it is on the circle.
pub(super) fn onto(c: &P2, r: &R, p: &P2) -> P2 {
    if on(c, r, p) {
        return p.clone();
    }
    let (dx, dy) = (&p[0] - &c[0], &p[1] - &c[1]);
    let right = dx >= zero();
    let den = if right { r + &dx } else { r - &dx };
    let s = q(rational_f64(&(&dy / &den)));
    let one = int(1);
    let w = &one + &s * &s;
    let a = r * (&one - &s * &s) / &w;
    let b = r * (int(2) * &s) / &w;
    [if right { &c[0] + a } else { &c[0] - a }, &c[1] + b]
}

/// A path's points with every arc's ends on its circle: each point taken
/// onto the circle of the arc ending or starting there.
pub(super) fn points(points: &[Point2], segments: &[Segment]) -> Result<Vec<P2>> {
    let n = points.len();
    let circle = |s: &Segment| match s {
        Segment::Arc { center, radius, .. } => Some(([q(center.x), q(center.y)], q(*radius))),
        _ => None,
    };
    let mut out = Vec::with_capacity(n);
    for j in 0..n {
        let p = [q(points[j].x), q(points[j].y)];
        let (before, after) = (circle(&segments[(j + n - 1) % n]), circle(&segments[j]));
        out.push(match (before, after) {
            (Some(x), Some(y)) if x != y => {
                if on(&x.0, &x.1, &p) && on(&y.0, &y.1, &p) {
                    p
                } else {
                    return Err(Error::OutOfDomain(
                        "an imported prism's arcs of two circles meeting at a joint (S9e.4b.4)",
                    ));
                }
            }
            (Some(x), _) | (None, Some(x)) => onto(&x.0, &x.1, &p),
            (None, None) => p,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> P2 {
        [q(x), q(y)]
    }

    #[test]
    fn ends_land_on_their_circles_near_where_they_were() {
        let (c, r) = (p(3.0, -1.0), q(2.0));
        for (x, y) in [
            (4.732050807568877, 0.0000000000000002),
            (1.267949192431123, 0.0),
            (3.0, 1.0000000000000004),
            (3.0, -3.0),
            (1.0, -1.0),
            (4.414213562373095, 0.41421356237309503),
        ] {
            let e = p(x, y);
            let s = onto(&c, &r, &e);
            assert!(on(&c, &r, &s), "{x} {y}");
            let d = [&s[0] - &e[0], &s[1] - &e[1]];
            let moved = crate::solid::split::rational_f64(&(&d[0] * &d[0] + &d[1] * &d[1])).sqrt();
            assert!(moved <= 1e-14, "{x} {y} moved {moved}");
        }
        // On the circle already: kept bit for bit.
        assert_eq!(onto(&c, &r, &p(5.0, -1.0)), p(5.0, -1.0));
    }

    #[test]
    fn joints_take_their_arcs_points_once() {
        // A stadium's arc end off its circle shares the moved point with
        // its line; two arcs of one circle share theirs; two circles are
        // refused unless the point lies on both.
        let line = Segment::Line;
        let arc = |cx: f64, cy: f64, r: f64| Segment::Arc {
            center: Point2::new(cx, cy),
            radius: r,
            ccw: true,
        };
        let pts = [
            Point2::new(-3.0, -2.0),
            Point2::new(3.0, -2.0000000000000004),
            Point2::new(3.0, 2.0),
            Point2::new(-3.0, 2.0),
        ];
        let segs = [
            line.clone(),
            arc(3.0, 0.0, 2.0),
            line.clone(),
            arc(-3.0, 0.0, 2.0),
        ];
        let out = points(&pts, &segs).unwrap();
        assert!(on(&p(3.0, 0.0), &q(2.0), &out[1]));
        assert_ne!(out[1], p(3.0, -2.0000000000000004));
        assert_eq!(out[2], p(3.0, 2.0));
        let halves = [Point2::new(5.5, 0.1), Point2::new(0.5, 0.1)];
        let both = [
            arc(3.0, 0.1, 2.5000000000000004),
            arc(3.0, 0.1, 2.5000000000000004),
        ];
        let out = points(&halves, &both).unwrap();
        assert!(out
            .iter()
            .all(|x| on(&p(3.0, 0.1), &q(2.5000000000000004), x)));
        let lens = [Point2::new(0.0, -4.0), Point2::new(0.0, 4.0)];
        assert!(points(&lens, &[arc(-3.0, 0.0, 5.0), arc(3.0, 0.0, 5.0)]).is_ok());
        let off = [Point2::new(0.0, -4.000000000000001), Point2::new(0.0, 4.0)];
        assert!(matches!(
            points(&off, &[arc(-3.0, 0.0, 5.0), arc(3.0, 0.0, 5.0)]),
            Err(Error::OutOfDomain(_))
        ));
    }

    /// The imported slot (`fixtures/imported/slot.brep`, its caps and walls
    /// tilted along `(0, 0.6, 0.8)`) with its caps' and cylinders' stored
    /// frames each platform's normalization of that axis: macOS's `hypot`
    /// keeps `(0, 0.5999999999999999, 0.8)` when normalized again, glibc's
    /// (correctly rounded) turns it to `(0, 0.6, 0.8000000000000002)`, which
    /// macOS's turns back. The construction's walls carry its cap's axes bit
    /// for bit, so a Boolean's split walls lie on the stored cylinders'
    /// exactly parallel axes on either frame: the history holds on any host.
    #[test]
    fn split_walls_keep_the_stored_axes_on_either_platforms_frames() {
        use crate::identity::OperationId;
        use crate::topology::{Surface, Topology};
        use crate::{Frame3, Point3, Solid, Vec3};
        let doc =
            crate::occt_brep::read(include_str!("../../../../../fixtures/imported/slot.brep"))
                .unwrap();
        let [solid] = <[_; 1]>::try_from(crate::occt_brep::import(&doc).solids).unwrap();
        let (stored, resolution) = (solid.result.unwrap(), solid.tolerance);
        let axis = |y: u64, z: u64| Vec3::new(0.0, f64::from_bits(y), f64::from_bits(z));
        // macOS's fixed point, glibc's.
        for n in [
            axis(0x3fe3_3333_3333_3332, 0x3fe9_9999_9999_999a),
            axis(0x3fe3_3333_3333_3333, 0x3fe9_9999_9999_999b),
        ] {
            // A frame of normal `+-n`, the stored `x` and `y = n x x`.
            let turned = |f: Frame3| {
                let n = if f.normal().dot(n) > 0.0 { n } else { -n };
                Frame3::from_axes(f.origin(), f.x(), n.cross(f.x()), n)
            };
            let mut parts = stored.to_parts();
            let mut cylinders = 0;
            for face in &mut parts.faces {
                match &mut face.surface {
                    Surface::Cylinder { frame, .. } => {
                        *frame = turned(*frame);
                        cylinders += 1;
                    }
                    Surface::Plane(frame) if frame.normal().cross(n).length() < 1e-12 => {
                        *frame = turned(*frame);
                    }
                    _ => {}
                }
            }
            assert_eq!(cylinders, 2);
            let t = Topology::from_parts(parts.with_measured_enclosures(), resolution).unwrap();
            let (slot, _) = Solid::imported_with(OperationId(91), t, resolution).unwrap();
            let (block, _) = Solid::box_at_with(
                OperationId(92),
                Point3::new(7.0, 0.0, 0.0),
                Vec3::new(5.0, 10.0, 1.5),
                resolution,
            )
            .unwrap();
            for (out, history) in [
                slot.fuse(OperationId(93), &block).unwrap(),
                slot.cut(OperationId(94), &block).unwrap(),
            ] {
                let ins = [
                    slot.topology().entity_set(slot.resolution()),
                    block.topology().entity_set(block.resolution()),
                ];
                let outs: Vec<_> = out
                    .iter()
                    .map(|s| s.topology().entity_set(s.resolution()))
                    .collect();
                let issues = crate::history::check(&ins, &outs, &history);
                assert!(issues.is_empty(), "{n:?}: {issues:?}");
            }
        }
    }
}
