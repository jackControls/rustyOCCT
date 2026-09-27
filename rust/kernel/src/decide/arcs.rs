//! Certified tolerance decisions for profile paths with circular arcs (S5 of
//! REVIEW_NOTES.md). An arc is its centre `c`, its radius `r` (an input
//! term), its end points `a` and `b` (the path's vertices, within tolerance
//! of the circle) and its direction. Its points are the circle's points in
//! the closed sector from `a - c` to `b - c` turning its way, and its ends are
//! taken as the vertices themselves. Distances that are sums of input terms
//! (a point's distance from the centre against `r` and the tolerance, the
//! distance between centres against radii) are decided exactly; the one
//! irrational construction, a line's or circle's crossing of a circle, is
//! located in rational intervals, and an undecided location counts as
//! touching: these decisions screen profiles for validity, as
//! `area_is_degenerate` does.
use super::{distance_ge, distance_le, q, segment_distance_le, sum_le, zero};
use crate::certified::Interval;
use crate::Point2;
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// A circular arc of a profile path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Arc2 {
    pub center: Point2,
    pub radius: f64,
    pub start: Point2,
    pub end: Point2,
    /// Counter-clockwise about the profile's normal.
    pub ccw: bool,
}

/// A piece of a profile boundary for the proximity screen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Piece {
    Line(Point2, Point2),
    Arc(Arc2),
    Circle(Point2, f64),
}

type V = (R, R);

fn vec_r(a: Point2, b: Point2) -> V {
    (q(b.x) - q(a.x), q(b.y) - q(a.y))
}

fn cross(u: &V, v: &V) -> R {
    &u.0 * &v.1 - &u.1 * &v.0
}

fn dot(u: &V, v: &V) -> R {
    &u.0 * &v.0 + &u.1 * &v.1
}

fn neg(u: &V) -> V {
    (-u.0.clone(), -u.1.clone())
}

impl Arc2 {
    /// The sector's bounding directions turning counter-clockwise, from the
    /// first to the second.
    fn ccw_bounds(&self) -> (V, V) {
        let (u, v) = (vec_r(self.center, self.start), vec_r(self.center, self.end));
        if self.ccw {
            (u, v)
        } else {
            (v, u)
        }
    }

    /// Whether the direction `w` (from the centre) lies in the closed
    /// sector, exactly. The zero direction lies in none.
    fn contains(&self, w: &V) -> bool {
        if w.0 == zero() && w.1 == zero() {
            return false;
        }
        let (u, v) = self.ccw_bounds();
        let (uw, wv) = (cross(&u, w), cross(w, &v));
        match cross(&u, &v).cmp(&zero()) {
            // The minor sector from u to v.
            Ordering::Greater => uw >= zero() && wv >= zero(),
            // The complement of the open minor sector from v to u.
            Ordering::Less => !(uw < zero() && wv < zero()),
            // A half turn: the closed half-plane left of u.
            Ordering::Equal => uw >= zero(),
        }
    }

    /// The same for a direction in intervals: `None` when a side is not
    /// certain.
    fn contains_interval(&self, w: &(Interval, Interval)) -> Option<bool> {
        let (u, v) = self.ccw_bounds();
        let ex = |x: &R| Interval::exact(x.clone());
        // cross(u, w) and cross(w, v).
        let uw = ex(&u.0).mul(&w.1).sub(&ex(&u.1).mul(&w.0)).sign()?;
        let wv = w.0.mul(&ex(&v.1)).sub(&w.1.mul(&ex(&v.0))).sign()?;
        Some(match cross(&u, &v).cmp(&zero()) {
            Ordering::Greater => uw != Ordering::Less && wv != Ordering::Less,
            Ordering::Less => !(uw == Ordering::Less && wv == Ordering::Less),
            Ordering::Equal => uw != Ordering::Less,
        })
    }

    /// The unit-free tangent at an end, in the direction of travel.
    fn tangent(&self, at_end: bool) -> V {
        let p = if at_end { self.end } else { self.start };
        let r = vec_r(self.center, p);
        // Counter-clockwise travel turns the radius a quarter left.
        let left = (-r.1.clone(), r.0.clone());
        if self.ccw {
            left
        } else {
            neg(&left)
        }
    }
}

/// `|d - r| <= t` for `d^2` exact, `r` and `t` input terms.
fn radial_within(d2: &R, r: f64, t: f64) -> bool {
    let (r, t) = (q(r), q(t));
    let outer = &r + &t;
    if d2 > &(&outer * &outer) {
        return false;
    }
    let inner = &r - &t;
    inner <= zero() || d2 >= &(&inner * &inner)
}

/// The distance from `p` to the arc is at most `t`.
pub(crate) fn point_arc_within(p: Point2, arc: &Arc2, t: f64) -> bool {
    if distance_le(p, arc.start, &[t]) || distance_le(p, arc.end, &[t]) {
        return true;
    }
    arc.contains(&vec_r(arc.center, p))
        && distance_le(p, arc.center, &[arc.radius, t])
        && distance_ge(p, arc.center, &[arc.radius, -t])
}

fn point_circle_within(p: Point2, c: Point2, r: f64, t: f64) -> bool {
    distance_le(p, c, &[r, t]) && distance_ge(p, c, &[r, -t])
}

/// The point of the segment `s0 s1` nearest to `c`, when it lies strictly
/// inside the segment.
fn foot(c: Point2, s0: Point2, s1: Point2) -> Option<V> {
    let e = vec_r(s0, s1);
    let w = vec_r(s0, c);
    let len2 = dot(&e, &e);
    let t = dot(&w, &e);
    if t <= zero() || t >= len2 {
        return None;
    }
    let k = t / len2;
    Some((q(s0.x) + &k * &e.0, q(s0.y) + &k * &e.1))
}

/// The segment's line crosses the circle: the crossings' parameters along
/// the segment and their directions from the centre, in intervals.
fn line_circle(c: Point2, r: f64, s0: Point2, s1: Point2) -> Vec<(Interval, (Interval, Interval))> {
    let e = vec_r(s0, s1);
    let w = vec_r(c, s0);
    let a = dot(&e, &e);
    let b = dot(&e, &w);
    let cc = dot(&w, &w) - q(r) * q(r);
    let disc = &b * &b - &a * &cc;
    if disc < zero() {
        return Vec::new();
    }
    let root = Interval::exact(disc).sqrt();
    let mut out = Vec::new();
    for sign in [-1, 1] {
        let tau = Interval::exact(-b.clone())
            .add(&if sign < 0 { root.neg() } else { root.clone() })
            .div(&Interval::exact(a.clone()))
            .expect("a segment has length");
        let dir = (
            Interval::exact(w.0.clone()).add(&tau.mul(&Interval::exact(e.0.clone()))),
            Interval::exact(w.1.clone()).add(&tau.mul(&Interval::exact(e.1.clone()))),
        );
        out.push((tau, dir));
    }
    out
}

/// `0 <= tau <= 1` certainly (`Some(true)`), certainly not, or undecided.
fn in_unit(tau: &Interval) -> Option<bool> {
    let one = R::from_integer(1.into());
    if tau.hi() < &zero() || tau.lo() > &one {
        return Some(false);
    }
    if tau.lo() >= &zero() && tau.hi() <= &one {
        return Some(true);
    }
    None
}

fn segment_arc_within(s0: Point2, s1: Point2, arc: &Arc2, t: f64) -> bool {
    if point_arc_within(s0, arc, t)
        || point_arc_within(s1, arc, t)
        || segment_distance_le(arc.start, s0, s1, &[t])
        || segment_distance_le(arc.end, s0, s1, &[t])
    {
        return true;
    }
    // The pair nearest across the interior: the foot of the centre on the
    // segment and the circle's point in its direction.
    if let Some(f) = foot(arc.center, s0, s1) {
        let w = (&f.0 - q(arc.center.x), &f.1 - q(arc.center.y));
        if arc.contains(&w) && radial_within(&dot(&w, &w), arc.radius, t) {
            return true;
        }
    }
    // A crossing inside both pieces, or one not located.
    line_circle(arc.center, arc.radius, s0, s1)
        .iter()
        .any(|(tau, dir)| match in_unit(tau) {
            Some(false) => false,
            Some(true) => arc.contains_interval(dir) != Some(false),
            None => true,
        })
}

/// Crossings of two circles as interval directions from each centre.
fn circle_circle(
    c1: Point2,
    r1: f64,
    c2: Point2,
    r2: f64,
) -> Vec<((Interval, Interval), (Interval, Interval))> {
    let e = vec_r(c1, c2);
    let d2 = dot(&e, &e);
    let (r1, r2) = (q(r1), q(r2));
    // |r1 - r2| <= d <= r1 + r2, squared.
    let (sum, diff) = (&r1 + &r2, &r1 - &r2);
    if d2 == zero() || d2 > &sum * &sum || d2 < &diff * &diff {
        return Vec::new();
    }
    // Along the axis from c1: x = (d^2 + r1^2 - r2^2) / (2 d); across it:
    // h = sqrt(r1^2 - x^2). With k = x / d and m = h / d, a crossing is
    // c1 + k e +- m e_perp.
    let k = (&d2 + &r1 * &r1 - &r2 * &r2) / (R::from_integer(2.into()) * &d2);
    let m2 = &r1 * &r1 / &d2 - &k * &k;
    if m2 < zero() {
        return Vec::new();
    }
    let m = Interval::exact(m2).sqrt();
    let ex = |x: &R| Interval::exact(x.clone());
    let mut out = Vec::new();
    for sign in [-1, 1] {
        let mm = if sign < 0 { m.neg() } else { m.clone() };
        // From c1: k e + mm (-e.y, e.x).
        let from1 = (
            ex(&(&k * &e.0)).sub(&mm.mul(&ex(&e.1))),
            ex(&(&k * &e.1)).add(&mm.mul(&ex(&e.0))),
        );
        // From c2: the same point minus e.
        let from2 = (from1.0.sub(&ex(&e.0)), from1.1.sub(&ex(&e.1)));
        out.push((from1, from2));
    }
    out
}

fn arc_arc_within(a: &Arc2, b: &Arc2, t: f64) -> bool {
    if point_arc_within(b.start, a, t)
        || point_arc_within(b.end, a, t)
        || point_arc_within(a.start, b, t)
        || point_arc_within(a.end, b, t)
    {
        return true;
    }
    arcs_interior_within(a, b, t)
}

/// The interior candidates of two arcs: the line-of-centres pairs and the
/// crossings.
fn arcs_interior_within(a: &Arc2, b: &Arc2, t: f64) -> bool {
    let e = vec_r(a.center, b.center);
    if e.0 == zero() && e.1 == zero() {
        // Concentric: close radii and overlapping sectors.
        let close = sum_le(&[a.radius, -b.radius], &[t]) && sum_le(&[b.radius, -a.radius], &[t]);
        return close
            && (a.contains(&vec_r(b.center, b.start)) || b.contains(&vec_r(a.center, a.start)));
    }
    for s1 in [1i8, -1] {
        for s2 in [1i8, -1] {
            let w1 = if s1 > 0 { e.clone() } else { neg(&e) };
            let w2 = if s2 > 0 { e.clone() } else { neg(&e) };
            if !(a.contains(&w1) && b.contains(&w2)) {
                continue;
            }
            // |d + s2 r2 - s1 r1| <= t.
            let (x1, x2) = (f64::from(s1) * a.radius, f64::from(s2) * b.radius);
            if distance_le(a.center, b.center, &[x1, -x2, t])
                && distance_ge(a.center, b.center, &[x1, -x2, -t])
            {
                return true;
            }
        }
    }
    circle_circle(a.center, a.radius, b.center, b.radius)
        .iter()
        .any(|(d1, d2)| {
            !matches!(
                (a.contains_interval(d1), b.contains_interval(d2)),
                (Some(false), _) | (_, Some(false))
            )
        })
}

/// Two non-adjacent pieces come within `t` of each other (or their
/// proximity is undecided).
pub(crate) fn pieces_within(a: &Piece, b: &Piece, t: f64) -> bool {
    use Piece::*;
    match (a, b) {
        (Line(a0, a1), Line(b0, b1)) => crate::profile::segments_touch(*a0, *a1, *b0, *b1, t),
        (Line(s0, s1), Arc(arc)) | (Arc(arc), Line(s0, s1)) => segment_arc_within(*s0, *s1, arc, t),
        (Arc(x), Arc(y)) => arc_arc_within(x, y, t),
        (Circle(c, r), Line(s0, s1)) | (Line(s0, s1), Circle(c, r)) => {
            segment_distance_le(*c, *s0, *s1, &[*r, t])
                && (distance_ge(*c, *s0, &[*r, -t]) || distance_ge(*c, *s1, &[*r, -t]))
        }
        (Circle(c, r), Arc(arc)) | (Arc(arc), Circle(c, r)) => {
            point_circle_within(arc.start, *c, *r, t)
                || point_circle_within(arc.end, *c, *r, t)
                || circle_arc_interior(*c, *r, arc, t)
        }
        (Circle(c1, r1), Circle(c2, r2)) => {
            let (big, small) = (r1.max(*r2), r1.min(*r2));
            distance_le(*c1, *c2, &[*r1, *r2, t]) && distance_ge(*c1, *c2, &[big, -small, -t])
        }
    }
}

fn circle_arc_interior(c: Point2, r: f64, arc: &Arc2, t: f64) -> bool {
    let e = vec_r(arc.center, c);
    if e.0 == zero() && e.1 == zero() {
        return sum_le(&[arc.radius, -r], &[t]) && sum_le(&[r, -arc.radius], &[t]);
    }
    for s1 in [1i8, -1] {
        let w1 = if s1 > 0 { e.clone() } else { neg(&e) };
        if !arc.contains(&w1) {
            continue;
        }
        for s2 in [1i8, -1] {
            let (x1, x2) = (f64::from(s1) * arc.radius, f64::from(s2) * r);
            if distance_le(arc.center, c, &[x1, -x2, t])
                && distance_ge(arc.center, c, &[x1, -x2, -t])
            {
                return true;
            }
        }
    }
    circle_circle(arc.center, arc.radius, c, r)
        .iter()
        .any(|(d1, _)| arc.contains_interval(d1) != Some(false))
}

/// A piece's end point and the tangent there, in the direction of travel.
fn end_and_tangent(p: &Piece, at_end: bool) -> (Point2, V) {
    match p {
        Piece::Line(a, b) => (if at_end { *b } else { *a }, vec_r(*a, *b)),
        Piece::Arc(arc) => (
            if at_end { arc.end } else { arc.start },
            arc.tangent(at_end),
        ),
        Piece::Circle(..) => unreachable!("a circle is a whole boundary"),
    }
}

fn point_piece_within(p: Point2, piece: &Piece, t: f64) -> bool {
    match piece {
        Piece::Line(a, b) => segment_distance_le(p, *a, *b, &[t]),
        Piece::Arc(arc) => point_arc_within(p, arc, t),
        Piece::Circle(c, r) => point_circle_within(p, *c, *r, t),
    }
}

/// The second meeting of two pieces through the shared point `v` (a line's
/// or circle's other crossing of a circle through `v`, taking `v` on both):
/// the reflection of `v`, exactly.
fn second_meeting(a: &Piece, b: &Piece, v: Point2) -> Option<(V, bool)> {
    let vv = (q(v.x), q(v.y));
    match (a, b) {
        (Piece::Line(..), Piece::Line(..)) => None,
        (Piece::Line(s0, s1), Piece::Arc(arc)) | (Piece::Arc(arc), Piece::Line(s0, s1)) => {
            // v - 2 ((v - c).e / e.e) e along the line.
            let e = vec_r(*s0, *s1);
            let w = vec_r(arc.center, v);
            let k = R::from_integer(2.into()) * dot(&w, &e) / dot(&e, &e);
            let x = (&vv.0 - &k * &e.0, &vv.1 - &k * &e.1);
            // On the segment: its parameter strictly inside.
            let from0 = (&x.0 - q(s0.x), &x.1 - q(s0.y));
            let tau = dot(&from0, &e) / dot(&e, &e);
            let inside = tau > zero() && tau < R::from_integer(1.into());
            let dir = (&x.0 - q(arc.center.x), &x.1 - q(arc.center.y));
            Some((x, inside && arc.contains(&dir)))
        }
        (Piece::Arc(x1), Piece::Arc(x2)) => {
            // The reflection of v across the line of centres.
            let e = vec_r(x1.center, x2.center);
            if e.0 == zero() && e.1 == zero() {
                return None;
            }
            let w = vec_r(x1.center, v);
            let k = dot(&w, &e) / dot(&e, &e);
            let proj = (q(x1.center.x) + &k * &e.0, q(x1.center.y) + &k * &e.1);
            let two = R::from_integer(2.into());
            let x = (&two * &proj.0 - &vv.0, &two * &proj.1 - &vv.1);
            let d1 = (&x.0 - q(x1.center.x), &x.1 - q(x1.center.y));
            let d2 = (&x.0 - q(x2.center.x), &x.1 - q(x2.center.y));
            Some((x, x1.contains(&d1) && x2.contains(&d2)))
        }
        _ => None,
    }
}

/// Two pieces of one path meeting at `v` (the first ends there, the second
/// starts there) are an invalid join: they double back, meet again away
/// from their shared points, or one's far end lies within `t` of the other.
/// `others` are the pieces' other shared points (a two-piece path meets at
/// both ends).
pub(crate) fn adjacent_invalid(a: &Piece, b: &Piece, t: f64, others: &[Point2]) -> bool {
    let (v, ta) = end_and_tangent(a, true);
    let (_, tb) = end_and_tangent(b, false);
    // Doubling back: the tangents exactly opposite.
    if cross(&ta, &tb) == zero() && dot(&ta, &tb) < zero() {
        return true;
    }
    let (a_start, _) = end_and_tangent(a, false);
    let (b_end, _) = end_and_tangent(b, true);
    let shared = |p: Point2| others.iter().any(|o| distance_le(p, *o, &[t]));
    if (!shared(b_end) && point_piece_within(b_end, a, t))
        || (!shared(a_start) && point_piece_within(a_start, b, t))
    {
        return true;
    }
    if let Some((x, on_both)) = second_meeting(a, b, v) {
        let near = |p: Point2| {
            let d = (&x.0 - q(p.x), &x.1 - q(p.y));
            let lim = q(2.0 * t);
            dot(&d, &d) <= &lim * &lim
        };
        if on_both && !near(v) && !others.iter().any(|o| near(*o)) {
            return true;
        }
    }
    // Interior proximity away from the join.
    match (a, b) {
        (Piece::Line(s0, s1), Piece::Arc(arc)) | (Piece::Arc(arc), Piece::Line(s0, s1)) => {
            if let Some(f) = foot(arc.center, *s0, *s1) {
                let w = (&f.0 - q(arc.center.x), &f.1 - q(arc.center.y));
                let fv = (&f.0 - q(v.x), &f.1 - q(v.y));
                let lim = q(2.0 * t);
                let at_join = dot(&fv, &fv) <= &lim * &lim;
                if !at_join && arc.contains(&w) && radial_within(&dot(&w, &w), arc.radius, t) {
                    return true;
                }
            }
            false
        }
        (Piece::Arc(x), Piece::Arc(y)) => arcs_interior_line_of_centres(x, y, t, v, others),
        _ => false,
    }
}

/// Line-of-centres proximity of two adjacent arcs away from their shared
/// points (a tangent join has its candidate there).
fn arcs_interior_line_of_centres(a: &Arc2, b: &Arc2, t: f64, v: Point2, others: &[Point2]) -> bool {
    let e = vec_r(a.center, b.center);
    if e.0 == zero() && e.1 == zero() {
        // Concentric neighbours overlap only by doubling back.
        return false;
    }
    let lim = q(2.0 * t);
    for s1 in [1i8, -1] {
        for s2 in [1i8, -1] {
            let w1 = if s1 > 0 { e.clone() } else { neg(&e) };
            let w2 = if s2 > 0 { e.clone() } else { neg(&e) };
            if !(a.contains(&w1) && b.contains(&w2)) {
                continue;
            }
            let (x1, x2) = (f64::from(s1) * a.radius, f64::from(s2) * b.radius);
            if !(distance_le(a.center, b.center, &[x1, -x2, t])
                && distance_ge(a.center, b.center, &[x1, -x2, -t]))
            {
                continue;
            }
            // The candidate point on a: c1 + s1 r1 e/|e|, in intervals.
            let d = Interval::exact(dot(&e, &e)).sqrt();
            let scale = Interval::from_f64(x1).div(&d);
            let Some(scale) = scale else { return true };
            let p = (
                Interval::from_f64(a.center.x).add(&scale.mul(&Interval::exact(e.0.clone()))),
                Interval::from_f64(a.center.y).add(&scale.mul(&Interval::exact(e.1.clone()))),
            );
            let near = |o: Point2| {
                let dx = p.0.sub(&Interval::from_f64(o.x));
                let dy = p.1.sub(&Interval::from_f64(o.y));
                dx.square().add(&dy.square()).hi() <= &(&lim * &lim)
            };
            if !(near(v) || others.iter().any(|o| near(*o))) {
                return true;
            }
        }
    }
    false
}

/// Crossings (0 or 1 per monotone part) of the ray `x > p.x`, `y = p.y` with
/// a piece, counted half-open in `y` as polygon edges are: a line or a
/// part of an arc between points on different sides of `p.y` crosses once
/// when its crossing lies right of `p`. Arcs split at the circle's top and
/// bottom into parts monotone in `y`. Exact; `p` is not within tolerance of
/// the piece (checked first by the caller).
pub(crate) fn ray_crossings(p: Point2, piece: &Piece) -> u32 {
    use crate::predicates::{orient2d_finite, Orientation2};
    match piece {
        Piece::Line(a, b) => {
            if (a.y > p.y) == (b.y > p.y) {
                return 0;
            }
            let side = orient2d_finite(*a, *b, p);
            u32::from(
                (b.y > a.y && side == Orientation2::CounterClockwise)
                    || (b.y < a.y && side == Orientation2::Clockwise),
            )
        }
        Piece::Arc(arc) => {
            let (cx, cy, r) = (q(arc.center.x), q(arc.center.y), q(arc.radius));
            let (top, bottom) = (&cy + &r, &cy - &r);
            let (ux, uy) = (q(arc.start.x) - &cx, q(arc.start.y) - &cy);
            // The circle's extreme points inside the sector, in travel order.
            let up: V = (zero(), R::from_integer(1.into()));
            let down: V = (zero(), R::from_integer((-1).into()));
            // An extreme point at an end is that end's vertex, not a break.
            let at_end = |w: &V| {
                [arc.start, arc.end].iter().any(|e| {
                    let d = vec_r(arc.center, *e);
                    cross(&d, w) == zero() && dot(&d, w) > zero()
                })
            };
            let (has_top, has_bottom) = (
                arc.contains(&up) && !at_end(&up),
                arc.contains(&down) && !at_end(&down),
            );
            // Travelling counter-clockwise from the right half the top comes
            // first; clockwise, the bottom.
            let right = ux > zero() || (ux == zero() && uy < zero());
            let top_first = if arc.ccw { right } else { !right };
            let mut ys = vec![q(arc.start.y)];
            let order: [(bool, &R); 2] = if top_first {
                [(has_top, &top), (has_bottom, &bottom)]
            } else {
                [(has_bottom, &bottom), (has_top, &top)]
            };
            for (present, y) in order {
                if present {
                    ys.push(y.clone());
                }
            }
            ys.push(q(arc.end.y));
            let py = q(p.y);
            let d = q(p.x) - &cx;
            let dy = &py - &cy;
            let qq = &r * &r - &dy * &dy;
            let mut count = 0;
            for w in ys.windows(2) {
                let (y0, y1) = (&w[0], &w[1]);
                if (y0 > &py) == (y1 > &py) {
                    continue;
                }
                // Rising counter-clockwise (or falling clockwise) runs up the
                // right half of the circle.
                let rising = y1 > y0;
                let right_half = rising == arc.ccw;
                let crossing_right = if right_half {
                    d < zero() || qq > &d * &d
                } else {
                    d < zero() && qq < &d * &d
                };
                count += u32::from(crossing_right);
            }
            count
        }
        Piece::Circle(..) => unreachable!("a path has no circle"),
    }
}
