//! S9e.4b.1: an imported prism's arcs' ends taken onto their circles
//! (REVIEW_NOTES.md, "S9e.4b refined"); S9e.4b.4a: its arcs of two circles
//! meeting at a joint taken through their ends ("S9e.4b.4 refined").
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
//! resolution of its circle).
//!
//! S9e.4b.4a: where arcs of two different circles meet at a joint whose
//! rounded point lies off either (a crossing joint: four discs' common, a
//! fillet chain, an arc tangent inside another), their common point is a
//! quadratic surd, which a profile segment's rational ends do not hold.
//! Instead each arc ending there whose circle no other arc of the path
//! shares is taken through its two ends: its circle in the exact model the
//! one through both, of centre `a + rho e` from its start `a`, `e` the
//! rational unit vector nearest the stored centre's direction (the same
//! half-angle tangent) and `rho = |b - a|^2 / (2 (b - a) . e)` exactly, its
//! centre and radius within the resolution of the stored ones. The joint is
//! then its rounded point, or the point of an arc's circle kept (a circle
//! several arcs share, or one meeting lines only). A joint tangent in the
//! body OCCT was given is tangent in no rounded data; the arrangement
//! decides it by its exact turn, as a line tangent to its arc.
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

/// The rational unit vector nearest `d`'s direction where `|d|` is about
/// `len`: `(+-(1 - s^2), 2 s) / (1 + s^2)`, its `x` of `d`'s sign, `s` the
/// binary64 rounding of `d_y / (len + |d_x|)` (the half-angle tangent from
/// the nearer end of the `x` diameter, `|s| <= 1` within rounding).
fn unit(d: &P2, len: &R) -> P2 {
    let right = d[0] >= zero();
    let den = if right { len + &d[0] } else { len - &d[0] };
    let s = q(rational_f64(&(&d[1] / &den)));
    let one = int(1);
    let w = &one + &s * &s;
    let a = (&one - &s * &s) / &w;
    [if right { a } else { -a }, int(2) * &s / &w]
}

/// The circle's rational point nearest `p`'s direction from its centre
/// (`c + r unit(p - c)`); `p` itself where it is on the circle.
pub(super) fn onto(c: &P2, r: &R, p: &P2) -> P2 {
    if on(c, r, p) {
        return p.clone();
    }
    let e = unit(&[&p[0] - &c[0], &p[1] - &c[1]], r);
    [&c[0] + r * &e[0], &c[1] + r * &e[1]]
}

/// S9e.4b.4a: the circle through `a` and `b` whose centre lies from `a` in
/// the rational direction `e` nearest the stored centre `c`'s (`unit`, the
/// distance `|c - a|` about `r`), of radius `rho = |b - a|^2 / (2 (b - a) .
/// e)`, so that `|(b - a) - rho e|^2 = |b - a|^2 - 2 rho (b - a) . e +
/// rho^2 = rho^2`. A chord seen from its end along the radius there has
/// `(b - a) . e > 0`, so for ends within rounding of the stored circle
/// `rho` is positive.
fn through(a: &P2, b: &P2, c: &P2, r: &R) -> Option<(P2, R)> {
    let e = unit(&[&c[0] - &a[0], &c[1] - &a[1]], r);
    let ab = [&b[0] - &a[0], &b[1] - &a[1]];
    let den = int(2) * (&ab[0] * &e[0] + &ab[1] * &e[1]);
    if den <= zero() {
        return None;
    }
    let rho = (&ab[0] * &ab[0] + &ab[1] * &ab[1]) / den;
    Some(([&a[0] + &rho * &e[0], &a[1] + &rho * &e[1]], rho))
}

/// A path's points and its arcs' circles for the exact model, every arc's
/// ends on its circle (the module's rules): the points, and per segment the
/// circle an arc is taken through its ends on (S9e.4b.4a), `None` where it
/// keeps its stored circle. `tol` bounds how far a circle taken through its
/// ends may lie from the stored one (the body's resolution).
#[allow(clippy::type_complexity)]
pub(super) fn path(
    points: &[Point2],
    segments: &[Segment],
    tol: f64,
) -> Result<(Vec<P2>, Vec<Option<(P2, R)>>)> {
    let n = points.len();
    let stored = |j: usize| match &segments[j % n] {
        Segment::Arc { center, radius, .. } => Some((center.x, center.y, *radius)),
        _ => None,
    };
    let circle = |j: usize| stored(j).map(|(x, y, r)| ([q(x), q(y)], q(r)));
    let point = |j: usize| [q(points[j % n].x), q(points[j % n].y)];
    // A circle another arc of the path shares stays (a seam-split circle).
    let shared = |j: usize| {
        let c = stored(j);
        c.is_some() && (0..n).any(|k| k != j % n && stored(k) == c)
    };
    // Joint `j`, between segments `j - 1` and `j`: arcs of two circles
    // meeting at a rounded point off either.
    let crossing = |j: usize| match (circle(j + n - 1), circle(j)) {
        (Some(x), Some(y)) if x != y => {
            let p = point(j);
            !(on(&x.0, &x.1, &p) && on(&y.0, &y.1, &p))
        }
        _ => false,
    };
    let through_ends: Vec<bool> = (0..n)
        .map(|j| stored(j).is_some() && !shared(j) && (crossing(j) || crossing(j + 1)))
        .collect();
    let kept = |j: usize| circle(j).filter(|_| !through_ends[j % n]);
    let mut out = Vec::with_capacity(n);
    for j in 0..n {
        let p = point(j);
        out.push(match (kept(j + n - 1), kept(j)) {
            (Some(x), Some(y)) if x != y => {
                if on(&x.0, &x.1, &p) && on(&y.0, &y.1, &p) {
                    p
                } else {
                    return Err(Error::OutOfDomain(
                        "an imported prism's joint of two circles each holding several arcs (S9e.4b.4)",
                    ));
                }
            }
            (Some(x), _) | (None, Some(x)) => onto(&x.0, &x.1, &p),
            (None, None) => p,
        });
    }
    let mut circles = Vec::with_capacity(n);
    for (j, taken) in through_ends.iter().enumerate() {
        if !taken {
            circles.push(None);
            continue;
        }
        let (c, r) = circle(j).expect("an arc");
        let fit = through(&out[j], &out[(j + 1) % n], &c, &r);
        let off = |x: &R, y: &R| rational_f64(&(x - y)).abs() > tol;
        match fit {
            Some((c2, r2)) if !(off(&c2[0], &c[0]) || off(&c2[1], &c[1]) || off(&r2, &r)) => {
                circles.push(Some((c2, r2)));
            }
            _ => {
                return Err(Error::Degenerate(
                    "an imported prism's arc too short to take through its ends",
                ))
            }
        }
    }
    Ok((out, circles))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> P2 {
        [q(x), q(y)]
    }

    fn arc(cx: f64, cy: f64, r: f64, ccw: bool) -> Segment {
        Segment::Arc {
            center: Point2::new(cx, cy),
            radius: r,
            ccw,
        }
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
        // its line; two arcs of one circle share theirs.
        let line = Segment::Line;
        let pts = [
            Point2::new(-3.0, -2.0),
            Point2::new(3.0, -2.0000000000000004),
            Point2::new(3.0, 2.0),
            Point2::new(-3.0, 2.0),
        ];
        let segs = [
            line.clone(),
            arc(3.0, 0.0, 2.0, true),
            line.clone(),
            arc(-3.0, 0.0, 2.0, true),
        ];
        let (out, circles) = path(&pts, &segs, 1e-7).unwrap();
        assert!(on(&p(3.0, 0.0), &q(2.0), &out[1]));
        assert_ne!(out[1], p(3.0, -2.0000000000000004));
        assert_eq!(out[2], p(3.0, 2.0));
        assert!(circles.iter().all(Option::is_none));
        let halves = [Point2::new(5.5, 0.1), Point2::new(0.5, 0.1)];
        let both = [
            arc(3.0, 0.1, 2.5000000000000004, true),
            arc(3.0, 0.1, 2.5000000000000004, true),
        ];
        let (out, circles) = path(&halves, &both, 1e-7).unwrap();
        assert!(out
            .iter()
            .all(|x| on(&p(3.0, 0.1), &q(2.5000000000000004), x)));
        assert!(circles.iter().all(Option::is_none));
    }

    /// S9e.4b.4a: a lens whose joints lie on both circles keeps them; one
    /// whose joint rounds off either takes both arcs through their ends,
    /// each new circle through both its ends exactly, within rounding of
    /// the stored one.
    #[test]
    fn arcs_of_two_circles_are_taken_through_their_ends() {
        let segs = [arc(-3.0, 0.0, 5.0, true), arc(3.0, 0.0, 5.0, true)];
        let exact = [Point2::new(0.0, -4.0), Point2::new(0.0, 4.0)];
        let (out, circles) = path(&exact, &segs, 1e-7).unwrap();
        assert_eq!(out, vec![p(0.0, -4.0), p(0.0, 4.0)]);
        assert!(circles.iter().all(Option::is_none));
        let off = [Point2::new(0.0, -4.000000000000001), Point2::new(0.0, 4.0)];
        let (out, circles) = path(&off, &segs, 1e-7).unwrap();
        assert_eq!(out, vec![p(0.0, -4.000000000000001), p(0.0, 4.0)]);
        for (j, (stored, c)) in [(p(-3.0, 0.0), 5.0), (p(3.0, 0.0), 5.0)]
            .into_iter()
            .enumerate()
        {
            let (centre, radius) = circles[j].clone().expect("taken through its ends");
            assert!(on(&centre, &radius, &out[j]) && on(&centre, &radius, &out[(j + 1) % 2]));
            let near = |x: &R, y: &R| rational_f64(&(x - y)).abs() <= 1e-14;
            assert!(near(&centre[0], &stored[0]) && near(&centre[1], &stored[1]));
            assert!(near(&radius, &q(c)));
        }
        // Two arcs tangent from either side (an S curve between lines): the
        // joint off either, both arcs taken through their ends, the lines'
        // joints on the new circles.
        let pts = [
            Point2::new(-4.0, -3.0),
            Point2::new(4.0, -3.0),
            Point2::new(4.0, 1.0),
            Point2::new(2.0000000000000004, 3.0),
            Point2::new(0.0, 5.0),
            Point2::new(-4.0, 5.0),
        ];
        let segs = [
            Segment::Line,
            Segment::Line,
            arc(2.0, 1.0, 2.0, true),
            arc(2.0, 5.0, 2.0, false),
            Segment::Line,
            Segment::Line,
        ];
        let (out, circles) = path(&pts, &segs, 1e-7).unwrap();
        assert_eq!(out[3], p(2.0000000000000004, 3.0));
        for j in [2, 3] {
            let (centre, radius) = circles[j].clone().expect("taken through its ends");
            assert!(on(&centre, &radius, &out[j]) && on(&centre, &radius, &out[j + 1]));
        }
        assert!(circles[0].is_none() && circles[5].is_none());
    }

    /// A joint of two circles each holding other arcs is refused (none can
    /// move); an arc whose circle through its ends lies off the stored one
    /// by more than the tolerance is degenerate.
    #[test]
    fn shared_circles_and_short_arcs_are_refused() {
        let split = [
            Point2::new(0.0, -4.000000000000001),
            Point2::new(2.0, 0.0),
            Point2::new(0.0, 4.0),
            Point2::new(-2.0, 0.0),
        ];
        let segs = [
            arc(-3.0, 0.0, 5.0, true),
            arc(-3.0, 0.0, 5.0, true),
            arc(3.0, 0.0, 5.0, true),
            arc(3.0, 0.0, 5.0, true),
        ];
        assert!(matches!(
            path(&split, &segs, 1e-7),
            Err(Error::OutOfDomain(m)) if m.contains("several arcs")
        ));
        let segs = [arc(-3.0, 0.0, 5.0, true), arc(3.0, 0.0, 5.0, true)];
        let far = [Point2::new(0.0, -4.000001), Point2::new(0.0, 4.0)];
        assert!(path(&far, &segs, 1e-5).is_ok());
        assert!(matches!(
            path(&far, &segs, 1e-9),
            Err(Error::Degenerate(m)) if m.contains("too short")
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
