//! S8b: the proximity screen for profile paths with spline segments.
//!
//! Every segment is a chain of parts (a line, a circular arc over an angle
//! range, or a Bézier arc by its exact control points). Two segments that
//! share no point are apart when every pair of their parts is, which boxes
//! farther apart than the resolution certify, the larger part halved until
//! they do: lines and Bézier arcs by de Casteljau in outward binary64
//! intervals (each control point's box holds the exact one, so the part lies
//! in their hull), arcs by angle. Adjacent segments must not double back at
//! their shared point, their far ends must be apart from the other segment,
//! and a pair of their parts is exempt when both are seen from the shared
//! point, which they keep as an exact end, in disjoint directions (cones of
//! their control points' box corners, or of their own box's; they meet only
//! there), as the arcs' adjacency rules ignore the meeting itself. A pair
//! still undecided at the depth or work limit counts as touching. A spline
//! segment is screened against itself in simple pieces (control polygons
//! turning through less than a half-turn).
use super::arcs::Arc2;
use crate::certified::{Fast, Real};
use crate::Point2;
use num_rational::BigRational as R;

type P = [R; 2];

const MAX_DEPTH: usize = 48;
const BUDGET: usize = 400_000;

fn zero() -> R {
    R::from_integer(0.into())
}

/// The binary64 nearest an exact rational (within its outward enclosure).
pub(crate) fn to_f64(x: &R) -> f64 {
    let (lo, hi) = crate::certified::Interval::exact(x.clone()).bounds_f64();
    if lo == hi {
        lo
    } else {
        0.5 * lo + 0.5 * hi
    }
}

fn q(x: f64) -> R {
    R::from_float(x).expect("a finite binary64")
}

fn pt(p: Point2) -> P {
    [q(p.x), q(p.y)]
}

fn sub(a: &P, b: &P) -> P {
    [&a[0] - &b[0], &a[1] - &b[1]]
}

fn cross(a: &P, b: &P) -> R {
    &a[0] * &b[1] - &a[1] * &b[0]
}

fn dot(a: &P, b: &P) -> R {
    &a[0] * &b[0] + &a[1] * &b[1]
}

fn is_zero(a: &P) -> bool {
    (a[0] == zero()) && (a[1] == zero())
}

/// Outward binary64 bounds of an exact rational.
fn down(x: &R) -> f64 {
    let (lo, _) = crate::certified::Interval::exact(x.clone()).bounds_f64();
    lo
}

fn up(x: &R) -> f64 {
    let (_, hi) = crate::certified::Interval::exact(x.clone()).bounds_f64();
    hi
}

/// A part of a segment.
#[derive(Debug, Clone)]
pub(crate) enum Part {
    Line(P, P),
    /// An arc over angles `[lo, hi]` (`lo < hi`), with its exact end points
    /// at `lo` and `hi` when they are the segment's, and whether the segment
    /// runs from `lo` to `hi`.
    Arc {
        center: Point2,
        radius: f64,
        lo: f64,
        hi: f64,
        ends: [Option<P>; 2],
        forward: bool,
    },
    Bezier(Vec<P>),
}

/// An outward box `[x_lo, x_hi, y_lo, y_hi]`.
type Bx = [f64; 4];

impl Part {
    fn bbox(&self) -> Bx {
        match self {
            Part::Line(a, b) => [
                down(&a[0].clone().min(b[0].clone())),
                up(&a[0].clone().max(b[0].clone())),
                down(&a[1].clone().min(b[1].clone())),
                up(&a[1].clone().max(b[1].clone())),
            ],
            Part::Bezier(cps) => {
                let (mut xl, mut xh, mut yl, mut yh) = (
                    cps[0][0].clone(),
                    cps[0][0].clone(),
                    cps[0][1].clone(),
                    cps[0][1].clone(),
                );
                for c in &cps[1..] {
                    xl = xl.min(c[0].clone());
                    xh = xh.max(c[0].clone());
                    yl = yl.min(c[1].clone());
                    yh = yh.max(c[1].clone());
                }
                [down(&xl), up(&xh), down(&yl), up(&yh)]
            }
            Part::Arc {
                center,
                radius,
                lo,
                hi,
                ..
            } => {
                let angle = Fast::exact_f64(*lo).union(&Fast::exact_f64(*hi));
                let (c, s) = Fast::cos_sin(&angle);
                let r = Fast::exact_f64(*radius);
                let x = Fast::exact_f64(center.x).add(&r.mul(&c));
                let y = Fast::exact_f64(center.y).add(&r.mul(&s));
                let ((xl, xh), (yl, yh)) = (x.bounds_f64(), y.bounds_f64());
                [xl, xh, yl, yh]
            }
        }
    }

    fn halves(&self) -> [Part; 2] {
        match self {
            Part::Line(a, b) => {
                let two = R::from_integer(2.into());
                let m = [(&a[0] + &b[0]) / &two, (&a[1] + &b[1]) / &two];
                [Part::Line(a.clone(), m.clone()), Part::Line(m, b.clone())]
            }
            Part::Arc {
                center,
                radius,
                lo,
                hi,
                ends,
                forward,
            } => {
                let mid = 0.5 * lo + 0.5 * hi;
                [
                    Part::Arc {
                        center: *center,
                        radius: *radius,
                        lo: *lo,
                        hi: mid,
                        ends: [ends[0].clone(), None],
                        forward: *forward,
                    },
                    Part::Arc {
                        center: *center,
                        radius: *radius,
                        lo: mid,
                        hi: *hi,
                        ends: [None, ends[1].clone()],
                        forward: *forward,
                    },
                ]
            }
            Part::Bezier(cps) => {
                // de Casteljau at one half, exactly.
                let two = R::from_integer(2.into());
                let mut rows = vec![cps.clone()];
                while rows.last().expect("a row").len() > 1 {
                    let last = rows.last().expect("a row");
                    let next: Vec<P> = last
                        .windows(2)
                        .map(|w| [(&w[0][0] + &w[1][0]) / &two, (&w[0][1] + &w[1][1]) / &two])
                        .collect();
                    rows.push(next);
                }
                let left: Vec<P> = rows.iter().map(|r| r[0].clone()).collect();
                let right: Vec<P> = rows.iter().rev().map(|r| r[r.len() - 1].clone()).collect();
                [Part::Bezier(left), Part::Bezier(right)]
            }
        }
    }

    /// The part's end points that are exact (a line's and a Bézier arc's
    /// always, an arc's when they are its segment's), an arc's in angle
    /// order.
    fn ends(&self) -> [Option<P>; 2] {
        match self {
            Part::Line(a, b) => [Some(a.clone()), Some(b.clone())],
            Part::Bezier(cps) => [Some(cps[0].clone()), Some(cps[cps.len() - 1].clone())],
            Part::Arc { ends, .. } => ends.clone(),
        }
    }

    /// The same in the segment's order.
    fn segment_ends(&self) -> [Option<P>; 2] {
        match self {
            Part::Arc {
                ends,
                forward: false,
                ..
            } => [ends[1].clone(), ends[0].clone()],
            _ => self.ends(),
        }
    }

    /// Directions from `p` covering the part (all but `p` itself), when
    /// they span less than a half-turn: `None` to refine.
    fn sector(&self, p: &P) -> Option<Sector> {
        let [s, e] = self.ends();
        let at_start = s.as_ref() == Some(p);
        let at_end = e.as_ref() == Some(p);
        let vectors: Vec<P> = if at_start || at_end {
            match self {
                Part::Line(a, b) => vec![sub(if at_start { b } else { a }, p)],
                Part::Bezier(cps) => cps.iter().filter(|c| *c != p).map(|c| sub(c, p)).collect(),
                Part::Arc {
                    center,
                    radius,
                    lo,
                    hi,
                    ..
                } => {
                    if hi - lo >= std::f64::consts::FRAC_PI_2 {
                        return None;
                    }
                    // The arc leaves p along its tangent towards its far
                    // end and stays in the cone of the tangent and the far
                    // end's enclosure.
                    let r = sub(p, &pt(*center));
                    let tangent = if at_start {
                        [-r[1].clone(), r[0].clone()]
                    } else {
                        [r[1].clone(), -r[0].clone()]
                    };
                    let far = if at_start { *hi } else { *lo };
                    let (c, s) = Fast::cos_sin(&Fast::exact_f64(far));
                    let x = Fast::exact_f64(center.x).add(&Fast::exact_f64(*radius).mul(&c));
                    let y = Fast::exact_f64(center.y).add(&Fast::exact_f64(*radius).mul(&s));
                    let ((xl, xh), (yl, yh)) = (x.bounds_f64(), y.bounds_f64());
                    let mut v = vec![tangent];
                    for (cx, cy) in [(xl, yl), (xl, yh), (xh, yl), (xh, yh)] {
                        v.push(sub(&[q(cx), q(cy)], p));
                    }
                    v
                }
            }
        } else {
            let [xl, xh, yl, yh] = self.bbox();
            let (px, py) = (p[0].clone(), p[1].clone());
            let inside = q(xl) <= px && px <= q(xh) && q(yl) <= py && py <= q(yh);
            if inside {
                return None;
            }
            [(xl, yl), (xl, yh), (xh, yl), (xh, yh)]
                .iter()
                .map(|(x, y)| sub(&[q(*x), q(*y)], p))
                .collect()
        };
        Sector::of(&vectors)
    }
}

/// A cone of directions spanning less than a half-turn, `[lo, hi]`
/// counter-clockwise.
struct Sector {
    lo: P,
    hi: P,
}

impl Sector {
    fn of(vectors: &[P]) -> Option<Self> {
        if vectors.is_empty() || vectors.iter().any(is_zero) {
            return None;
        }
        // lo: every vector counter-clockwise of it within a half-turn.
        let within = |a: &P, b: &P| {
            let c = cross(a, b);
            c > zero() || (c == zero() && dot(a, b) > zero())
        };
        let lo = vectors
            .iter()
            .find(|v| vectors.iter().all(|w| w == *v || within(v, w)))?;
        let hi = vectors
            .iter()
            .find(|v| vectors.iter().all(|w| w == *v || within(w, v)))?;
        if !(lo == hi || cross(lo, hi) > zero()) {
            return None;
        }
        Some(Sector {
            lo: lo.clone(),
            hi: hi.clone(),
        })
    }

    fn contains(&self, v: &P) -> bool {
        let (a, b) = (cross(&self.lo, v), cross(v, &self.hi));
        if cross(&self.lo, &self.hi) == zero() {
            return a == zero() && dot(&self.lo, v) > zero();
        }
        a >= zero() && b >= zero()
    }

    fn disjoint(&self, other: &Sector) -> bool {
        !self.contains(&other.lo) && !other.contains(&self.lo)
    }
}

/// Whether two boxes are certainly farther apart than `t`.
fn boxes_apart(a: &Bx, b: &Bx, t: f64) -> bool {
    let gap = |lo: f64, hi: f64| {
        // max(0, lo - hi), rounded down.
        let d = Fast::exact_f64(lo).sub(&Fast::exact_f64(hi)).bounds_f64().0;
        d.max(0.0)
    };
    let gx = gap(b[0], a[1]).max(gap(a[0], b[1]));
    let gy = gap(b[2], a[3]).max(gap(a[2], b[3]));
    let (gx, gy) = (Fast::exact_f64(gx), Fast::exact_f64(gy));
    let d2 = gx.square().add(&gy.square()).bounds_f64().0;
    let t2 = Fast::exact_f64(t).square().bounds_f64().1;
    d2 > t2
}

fn diagonal(b: &Bx) -> f64 {
    (b[1] - b[0]).hypot(b[3] - b[2])
}

/// A segment's parts: a line, an arc (its start angle and signed sweep from
/// its start to its end point) or a spline's exact Bézier arcs.
pub(crate) fn line(a: Point2, b: Point2) -> Vec<Part> {
    vec![Part::Line(pt(a), pt(b))]
}

/// A whole circle as one arc part.
pub(crate) fn circle(center: Point2, radius: f64) -> Vec<Part> {
    vec![Part::Arc {
        center,
        radius,
        lo: 0.0,
        hi: std::f64::consts::TAU,
        ends: [None, None],
        forward: true,
    }]
}

pub(crate) fn arc(arc: &Arc2) -> Vec<Part> {
    let start = (arc.start.y - arc.center.y).atan2(arc.start.x - arc.center.x);
    let sweep = crate::profile::arc_sweep(arc.center, arc.start, arc.end, arc.ccw);
    let (lo, hi, ends) = if sweep > 0.0 {
        (
            start,
            start + sweep,
            [Some(pt(arc.start)), Some(pt(arc.end))],
        )
    } else {
        (
            start + sweep,
            start,
            [Some(pt(arc.end)), Some(pt(arc.start))],
        )
    };
    vec![Part::Arc {
        center: arc.center,
        radius: arc.radius,
        lo,
        hi,
        ends,
        forward: sweep > 0.0,
    }]
}

/// A nonrational spline's Bézier arcs from its start point to its end
/// point.
pub(crate) fn spline(
    span: &crate::topology::SplineSpan<crate::BSplineCurve2>,
) -> Option<Vec<Part>> {
    let arcs = span.curve().as_curve3().bezier_arcs().ok()?;
    let mut parts: Vec<Part> = arcs
        .iter()
        .map(|a| {
            Part::Bezier(
                a.homogeneous_poles()
                    .iter()
                    .map(|h| [&h[0] / &h[3], &h[1] / &h[3]])
                    .collect(),
            )
        })
        .collect();
    if span.is_reversed() {
        parts.reverse();
        for p in &mut parts {
            if let Part::Bezier(cps) = p {
                cps.reverse();
            }
        }
    }
    Some(parts)
}

/// The exact outgoing direction of a segment's parts at its start (`start`)
/// or end: along the segment away from that point.
fn outgoing(parts: &[Part], start: bool) -> Option<P> {
    let part = if start {
        &parts[0]
    } else {
        &parts[parts.len() - 1]
    };
    match part {
        Part::Line(a, b) => Some(if start { sub(b, a) } else { sub(a, b) }),
        Part::Bezier(cps) => {
            let seq: Vec<&P> = if start {
                cps.iter().collect()
            } else {
                cps.iter().rev().collect()
            };
            seq[1..]
                .iter()
                .map(|c| sub(c, seq[0]))
                .find(|v| !is_zero(v))
        }
        Part::Arc { center, .. } => {
            // Away from the point along the segment: towards higher angles
            // from its `lo` end, lower from its `hi` end.
            let [s, e] = part.segment_ends();
            let p = if start { s? } else { e? };
            let at_lo = part.ends()[0].as_ref() == Some(&p);
            let r = sub(&p, &pt(*center));
            Some(if at_lo {
                [-r[1].clone(), r[0].clone()]
            } else {
                [r[1].clone(), -r[0].clone()]
            })
        }
    }
}

/// A part as the screen subdivides it: a line or a Bézier arc by outward
/// boxes of its control points and its exact ends where they are its
/// original part's; an arc as a part.
#[derive(Clone)]
enum Screened {
    Poly {
        boxes: Vec<[Fast; 2]>,
        ends: [Option<P>; 2],
    },
    Arc(Part),
}

fn point_box(p: &P) -> [Fast; 2] {
    [Fast::from_r(&p[0]), Fast::from_r(&p[1])]
}

/// A box's corners, exactly.
fn corners(b: &[Fast; 2]) -> [P; 4] {
    let ((xl, xh), (yl, yh)) = (b[0].bounds_f64(), b[1].bounds_f64());
    [
        [q(xl), q(yl)],
        [q(xl), q(yh)],
        [q(xh), q(yl)],
        [q(xh), q(yh)],
    ]
}

impl Screened {
    fn of(part: &Part) -> Self {
        match part {
            Part::Line(a, b) => Screened::Poly {
                boxes: vec![point_box(a), point_box(b)],
                ends: [Some(a.clone()), Some(b.clone())],
            },
            Part::Bezier(cps) => Screened::Poly {
                boxes: cps.iter().map(point_box).collect(),
                ends: [Some(cps[0].clone()), Some(cps[cps.len() - 1].clone())],
            },
            Part::Arc { .. } => Screened::Arc(part.clone()),
        }
    }

    fn bbox(&self) -> Bx {
        match self {
            Screened::Arc(part) => part.bbox(),
            Screened::Poly { boxes, .. } => {
                let mut out = [
                    f64::INFINITY,
                    f64::NEG_INFINITY,
                    f64::INFINITY,
                    f64::NEG_INFINITY,
                ];
                for b in boxes {
                    let ((xl, xh), (yl, yh)) = (b[0].bounds_f64(), b[1].bounds_f64());
                    out = [
                        out[0].min(xl),
                        out[1].max(xh),
                        out[2].min(yl),
                        out[3].max(yh),
                    ];
                }
                out
            }
        }
    }

    fn halves(&self) -> [Screened; 2] {
        match self {
            Screened::Arc(part) => part.halves().map(Screened::Arc),
            Screened::Poly { boxes, ends } => {
                let half = Fast::exact_f64(0.5);
                let mut rows = vec![boxes.clone()];
                while rows.last().expect("a row").len() > 1 {
                    let next: Vec<[Fast; 2]> = rows
                        .last()
                        .expect("a row")
                        .windows(2)
                        .map(|w| {
                            [
                                w[0][0].add(&w[1][0]).mul(&half),
                                w[0][1].add(&w[1][1]).mul(&half),
                            ]
                        })
                        .collect();
                    rows.push(next);
                }
                let left = rows.iter().map(|r| r[0]).collect();
                let right = rows.iter().rev().map(|r| r[r.len() - 1]).collect();
                [
                    Screened::Poly {
                        boxes: left,
                        ends: [ends[0].clone(), None],
                    },
                    Screened::Poly {
                        boxes: right,
                        ends: [None, ends[1].clone()],
                    },
                ]
            }
        }
    }

    /// Directions from `p` covering the part (all but `p` itself), when they
    /// span less than a half-turn: `None` to refine. From an exact end, the
    /// cone of the other control points' box corners (the part is their
    /// nonnegative combination); otherwise its box's corners, `p` outside
    /// it.
    fn sector(&self, p: &P) -> Option<Sector> {
        let Screened::Poly { boxes, ends } = self else {
            let Screened::Arc(part) = self else {
                unreachable!()
            };
            return part.sector(p);
        };
        let (at_start, at_end) = (ends[0].as_ref() == Some(p), ends[1].as_ref() == Some(p));
        let vectors: Vec<P> = if at_start || at_end {
            let own = if at_start { 0 } else { boxes.len() - 1 };
            boxes
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != own)
                .flat_map(|(_, b)| corners(b))
                .map(|c| sub(&c, p))
                .filter(|v| !is_zero(v))
                .collect()
        } else {
            let [xl, xh, yl, yh] = self.bbox();
            let (px, py) = (p[0].clone(), p[1].clone());
            if q(xl) <= px && px <= q(xh) && q(yl) <= py && py <= q(yh) {
                return None;
            }
            [(xl, yl), (xl, yh), (xh, yl), (xh, yh)]
                .iter()
                .map(|(x, y)| sub(&[q(*x), q(*y)], p))
                .collect()
        };
        Sector::of(&vectors)
    }
}

/// Screen state: pairs to examine and a work budget.
struct Screen {
    work: usize,
}

impl Screen {
    /// Whether every pair of parts is apart (or exempt); `false` when
    /// within `t` or undecided.
    fn pairs(&mut self, a: &[Part], b: &[Part], t: f64, joints: &[P]) -> bool {
        let (a, b): (Vec<Screened>, Vec<Screened>) = (
            a.iter().map(Screened::of).collect(),
            b.iter().map(Screened::of).collect(),
        );
        let mut stack: Vec<(Screened, Screened, usize)> = Vec::new();
        for x in &a {
            for y in &b {
                stack.push((x.clone(), y.clone(), 0));
            }
        }
        while let Some((x, y, depth)) = stack.pop() {
            self.work += 1;
            if self.work > BUDGET {
                return false;
            }
            let (bx, by) = (x.bbox(), y.bbox());
            if boxes_apart(&bx, &by, t) {
                continue;
            }
            // Parts seen from a shared point in disjoint directions meet
            // only there.
            let exempt = joints.iter().any(|j| match (x.sector(j), y.sector(j)) {
                (Some(sx), Some(sy)) => sx.disjoint(&sy),
                _ => false,
            });
            if exempt {
                continue;
            }
            if depth >= MAX_DEPTH {
                return false;
            }
            if diagonal(&bx) >= diagonal(&by) {
                for h in x.halves() {
                    stack.push((h, y.clone(), depth + 1));
                }
            } else {
                for h in y.halves() {
                    stack.push((x.clone(), h, depth + 1));
                }
            }
        }
        true
    }
}

/// Whether two segments sharing no point stay farther apart than `t`.
pub(crate) fn apart(a: &[Part], b: &[Part], t: f64) -> bool {
    Screen { work: 0 }.pairs(a, b, t, &[])
}

/// Whether two segments meeting at `a`'s end and `b`'s start (and, on a
/// two-segment path, also at `b`'s end and `a`'s start) form a valid join:
/// no doubling back at a shared point, far ends apart from the other
/// segment, and every pair of parts apart unless seen from a shared point
/// in disjoint directions.
pub(crate) fn adjacent_valid(a: &[Part], b: &[Part], t: f64, closed_pair: bool) -> bool {
    let joint = match a[a.len() - 1].segment_ends()[1].clone() {
        Some(p) => p,
        None => return false,
    };
    let mut joints = vec![joint];
    // Doubling back: the two outgoing directions at the joint coincide.
    let (Some(da), Some(db)) = (outgoing(a, false), outgoing(b, true)) else {
        return false;
    };
    if cross(&da, &db) == zero() && dot(&da, &db) > zero() {
        return false;
    }
    if closed_pair {
        let Some(other) = b[b.len() - 1].segment_ends()[1].clone() else {
            return false;
        };
        let (Some(da), Some(db)) = (outgoing(a, true), outgoing(b, false)) else {
            return false;
        };
        if cross(&da, &db) == zero() && dot(&da, &db) > zero() {
            return false;
        }
        joints.push(other);
    } else {
        // Far ends apart from the other segment.
        let (Some(fa), Some(fb)) = (
            a[0].segment_ends()[0].clone(),
            b[b.len() - 1].segment_ends()[1].clone(),
        ) else {
            return false;
        };
        let point = |p: &P| vec![Part::Line(p.clone(), p.clone())];
        let mut s = Screen { work: 0 };
        if !s.pairs(&point(&fa), b, t, &[]) || !s.pairs(&point(&fb), a, t, &[]) {
            return false;
        }
    }
    Screen { work: 0 }.pairs(a, b, t, &joints)
}

/// A spline's Bézier arcs cut into pieces whose control polygons turn
/// through less than a half-turn: each is monotone along a direction, so
/// simple. A cusp's halves end at its stationary point and double back
/// there (the self screen refuses them); `None` when a turn does not narrow
/// within the subdivision limit.
pub(crate) fn simple(parts: Vec<Part>) -> Option<Vec<Part>> {
    let narrow = |part: &Part| {
        let Part::Bezier(cps) = part else {
            return true;
        };
        let legs: Vec<P> = cps
            .windows(2)
            .map(|w| sub(&w[1], &w[0]))
            .filter(|d| !is_zero(d))
            .collect();
        Sector::of(&legs).is_some()
    };
    let mut out = Vec::new();
    let mut stack: Vec<(Part, usize)> = parts.into_iter().rev().map(|p| (p, 0)).collect();
    while let Some((part, depth)) = stack.pop() {
        if narrow(&part) {
            out.push(part);
            continue;
        }
        if depth >= SIMPLE_DEPTH {
            return None;
        }
        let [left, right] = part.halves();
        stack.push((right, depth + 1));
        stack.push((left, depth + 1));
    }
    Some(out)
}

/// Subdivisions of a Bézier arc to narrow its turn.
const SIMPLE_DEPTH: usize = 12;

/// Whether a segment's simple pieces stay farther apart than `t` from each
/// other, consecutive ones except at their shared point.
pub(crate) fn self_apart(parts: &[Part], t: f64) -> bool {
    let mut screen = Screen { work: 0 };
    for i in 0..parts.len() {
        for j in (i + 1)..parts.len() {
            let joints: Vec<P> = if j == i + 1 {
                match parts[i].segment_ends()[1].clone() {
                    Some(p) => vec![p],
                    None => return false,
                }
            } else {
                Vec::new()
            };
            if !screen.pairs(
                std::slice::from_ref(&parts[i]),
                std::slice::from_ref(&parts[j]),
                t,
                &joints,
            ) {
                return false;
            }
        }
    }
    true
}

/// Whether a point lies within `t` of a segment's parts (or undecided).
pub(crate) fn point_within(p: Point2, parts: &[Part], t: f64) -> bool {
    let point = [Part::Line(pt(p), pt(p))];
    !Screen { work: 0 }.pairs(&point, parts, t, &[])
}

/// Crossings of the `+x` ray from `p` with a chain's Bézier arcs, half-open
/// in `y` (an arc counts when exactly one of its ends lies strictly above
/// `p.y`), from exact subdivision into pieces monotone in `y`; `None` when
/// undecided (the point within rounding of the curve).
pub(crate) fn ray_crossings(p: Point2, parts: &[Part]) -> Option<u32> {
    let (px, py) = (q(p.x), q(p.y));
    let mut count = 0;
    let mut stack: Vec<(Vec<P>, usize)> = parts
        .iter()
        .filter_map(|part| match part {
            Part::Bezier(cps) => Some((cps.clone(), 0)),
            _ => None,
        })
        .collect();
    while let Some((cps, depth)) = stack.pop() {
        let (ylo, yhi) = cps
            .iter()
            .fold((cps[0][1].clone(), cps[0][1].clone()), |(l, h), c| {
                (l.min(c[1].clone()), h.max(c[1].clone()))
            });
        let xhi = cps.iter().map(|c| c[0].clone()).max().expect("a pole");
        let xlo = cps.iter().map(|c| c[0].clone()).min().expect("a pole");
        // Clear of the ray: off its line or left of the point.
        if yhi < py || ylo > py || xhi < px {
            continue;
        }
        let rises = cps.windows(2).all(|w| w[1][1] >= w[0][1]);
        let falls = cps.windows(2).all(|w| w[1][1] <= w[0][1]);
        let above = |c: &P| c[1] > py;
        let (first, last) = (&cps[0], &cps[cps.len() - 1]);
        if rises || falls {
            if above(first) == above(last) {
                continue;
            }
            if xlo > px {
                count += 1;
                continue;
            }
        }
        if depth >= 64 {
            return None;
        }
        if let [Part::Bezier(a), Part::Bezier(b)] = Part::Bezier(cps).halves() {
            stack.push((a, depth + 1));
            stack.push((b, depth + 1));
        }
    }
    Some(count)
}

/// A polynomial in the power basis over `[0, 1]`.
type Poly = Vec<R>;

fn binomial(n: usize, k: usize) -> R {
    let mut b = R::from_integer(1.into());
    for i in 0..k {
        b = b * R::from_integer(((n - i) as i64).into()) / R::from_integer(((i + 1) as i64).into());
    }
    b
}

/// A Bézier arc's coordinates in the power basis.
fn power(cps: &[P]) -> [Poly; 2] {
    let n = cps.len() - 1;
    let mut out = [vec![zero(); n + 1], vec![zero(); n + 1]];
    for (i, c) in cps.iter().enumerate() {
        // C(n, i) t^i (1 - t)^(n - i) = sum_j C(n, i) C(n - i, j) (-1)^j t^(i + j).
        for j in 0..=(n - i) {
            let mut k = binomial(n, i) * binomial(n - i, j);
            if j % 2 == 1 {
                k = -k;
            }
            for d in 0..2 {
                out[d][i + j] += &k * &c[d];
            }
        }
    }
    out
}

fn pmul(a: &[R], b: &[R]) -> Poly {
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn pder(a: &[R]) -> Poly {
    if a.len() == 1 {
        return vec![zero()];
    }
    a.iter()
        .enumerate()
        .skip(1)
        .map(|(k, x)| x * R::from_integer((k as i64).into()))
        .collect()
}

fn pint(a: &[R]) -> R {
    a.iter()
        .enumerate()
        .map(|(k, x)| x / R::from_integer(((k + 1) as i64).into()))
        .sum()
}

/// Green's-theorem integrals `[A, ∫x, ∫y, ∫x², ∫xy, ∫y²]` of the region left
/// of a chain of Bézier arcs, relative to `anchor`, exactly (the forms of
/// `profile::path_moments`: `A = ½∮(x dy - y dx)`, `∫x = ½∮x² dy`, `∫y =
/// -½∮y² dx`, `∫x² = ⅓∮x³ dy`, `∫xy = ½∮x² y dy`, `∫y² = -⅓∮y³ dx`).
pub(crate) fn green(parts: &[Part], anchor: &P) -> [R; 6] {
    let half = R::new(1.into(), 2.into());
    let third = R::new(1.into(), 3.into());
    let mut out: [R; 6] = std::array::from_fn(|_| zero());
    for part in parts {
        let Part::Bezier(cps) = part else { continue };
        let local: Vec<P> = cps.iter().map(|c| sub(c, anchor)).collect();
        let [x, y] = power(&local);
        let (dx, dy) = (pder(&x), pder(&y));
        let (xx, yy) = (pmul(&x, &x), pmul(&y, &y));
        out[0] += pint(&pmul(&x, &dy)) * &half - pint(&pmul(&y, &dx)) * &half;
        out[1] += pint(&pmul(&xx, &dy)) * &half;
        out[2] -= pint(&pmul(&yy, &dx)) * &half;
        out[3] += pint(&pmul(&pmul(&xx, &x), &dy)) * &third;
        out[4] += pint(&pmul(&pmul(&xx, &y), &dy)) * &half;
        out[5] -= pint(&pmul(&pmul(&yy, &y), &dx)) * &third;
    }
    out
}

/// Bounds on a chain's length: its chord (a lower bound) and its control
/// polygons' length (an upper bound), as squared lengths summed per leg.
pub(crate) fn length_bounds<T: Real>(parts: &[Part]) -> T {
    let (mut lo, mut hi) = (T::exact_f64(0.0), T::exact_f64(0.0));
    for part in parts {
        let Part::Bezier(cps) = part else { continue };
        let chord = sub(&cps[cps.len() - 1], &cps[0]);
        lo = lo.add(&T::from_r(&dot(&chord, &chord)).sqrt());
        for w in cps.windows(2) {
            let d = sub(&w[1], &w[0]);
            hi = hi.add(&T::from_r(&dot(&d, &d)).sqrt());
        }
    }
    lo.union(&hi)
}

/// A chain's length in binary64 (Gauss-Legendre on sixteen pieces of each
/// arc), for the profile's perimeter.
pub(crate) fn length(parts: &[Part]) -> f64 {
    const NODES: [(f64, f64); 5] = [
        (0.0, 0.568_888_888_888_888_9),
        (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
        (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
        (0.906_179_845_938_664, 0.236_926_885_056_189_08),
        (-0.906_179_845_938_664, 0.236_926_885_056_189_08),
    ];
    let mut total = 0.0;
    for part in parts {
        let Part::Bezier(cps) = part else { continue };
        let [x, y] = power(cps);
        let f = |p: &Poly| -> Vec<f64> { pder(p).iter().map(to_f64).collect() };
        let (dx, dy) = (f(&x), f(&y));
        let eval = |c: &[f64], t: f64| c.iter().rev().fold(0.0, |acc, k| acc * t + k);
        let pieces = 16;
        for k in 0..pieces {
            let (a, b) = (k as f64 / pieces as f64, (k + 1) as f64 / pieces as f64);
            for (node, w) in NODES {
                let t = 0.5 * (a + b) + 0.5 * (b - a) * node;
                total += 0.5 * (b - a) * w * eval(&dx, t).hypot(eval(&dy, t));
            }
        }
    }
    total
}
