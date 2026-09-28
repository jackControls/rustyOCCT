//! S9a: Booleans of two profiles (fuse, cut, common) by an exact
//! arrangement of their boundaries (S9 of REVIEW_NOTES.md).
//!
//! Every crossing of two stored segments is decided on the stored data:
//! two lines by their rational crossing, a line and a circle by the exact
//! sign of their discriminant and a quadratic surd, two circles by their
//! radical line. A crossing identifies with a stored vertex it is within the
//! resolution of when that vertex ends one of the two segments (the vertex
//! is then on the other segment within the resolution); two segments
//! overlapping along one line or one circle are split at each other's ends
//! and their pieces there are shared. Every other piece is classified
//! against the other profile at an off-centre point (the certified
//! `Profile::classify`; a piece within the resolution of the other's
//! boundary without being shared is `Degenerate`), the operation keeps its
//! pieces, and they are traced into cycles: each kept piece's end starts
//! exactly one kept piece (two would make the result touch itself there:
//! `Degenerate`), counter-clockwise cycles are outer boundaries and
//! clockwise ones holes, each result validated as a profile.
use super::{boundaries_touch, Boundary, BoundaryKind, Location, Profile, Segment};
mod splines;

use crate::certified::{Interval as I, Real};
use crate::solid::split::spline::Span;
use crate::solid::split::{circle_points, q, rounded, within_arc, zero, ArcPos};
use crate::{Error, Point2, Result, Tolerance};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

/// One of the Boolean's two profiles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Operand {
    A,
    B,
}

/// The operation: `A ∪ B`, `A − B` or `A ∩ B`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Op2 {
    Fuse,
    Cut,
    Common,
}

/// A point of the arrangement: a stored vertex of an operand, or a
/// crossing of two segments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum PId {
    Vertex(Operand, usize, usize),
    Cross(usize),
}

/// A stored segment `(operand, boundary, segment)` and, when cut, the part
/// of it (numbered along its stored direction).
pub(crate) type SegRef = (Operand, usize, usize, Option<usize>);

/// Where a result segment comes from: one operand's stored segment, whole
/// or a part, and the other's too when both boundaries share it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Origin2 {
    pub(crate) from: SegRef,
    pub(crate) shared: Option<SegRef>,
}

/// The Boolean's results and, per crossing, the two stored segments
/// `(operand, boundary, segment)` meeting there.
pub(crate) struct Boolean2 {
    pub(crate) pieces: Vec<Piece2>,
    pub(crate) crosses: BTreeMap<usize, [(Operand, usize, usize); 2]>,
}

/// A result profile with, per boundary (outer first) and per stored
/// segment, its provenance (the input pieces it joins, in order along it)
/// and its start point's.
pub(crate) struct Piece2 {
    pub(crate) profile: Profile,
    pub(crate) segments: Vec<Vec<Vec<Origin2>>>,
    pub(crate) vertices: Vec<Vec<PId>>,
    /// Per boundary (outer first), its arrangement pieces in the order the
    /// region lies on their left, each with whether it runs along its
    /// segment's stored direction (S9a.2's stacks).
    pub(crate) cycles: Vec<Vec<(usize, bool)>>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Shape {
    Line,
    Arc {
        center: Point2,
        radius: f64,
        ccw: bool,
    },
    /// A whole circle, stored counter-clockwise from angle 0.
    Circle {
        center: Point2,
        radius: f64,
    },
    /// A spline segment (S9a.2): its index among the arrangement's splines.
    Spline(usize),
}

/// A stored segment of an operand's boundary.
#[derive(Debug, Clone)]
struct Seg {
    op: Operand,
    b: usize,
    j: usize,
    /// Its stored ends (a circle's are its angle-0 point, twice, and none).
    from: Option<PId>,
    to: Option<PId>,
    p: Point2,
    e: Point2,
    shape: Shape,
}

/// A place along a segment: a line's parameter, or an arc's or circle's
/// position from its start (or angle 0).
#[derive(Debug, Clone)]
enum Key {
    Line(I),
    Arc(ArcPos),
    /// A spline's curve parameter, negated when it runs against the curve
    /// (so it increases along the stored direction).
    Spline(f64),
}

impl Key {
    fn compare(&self, o: &Self) -> Option<Ordering> {
        match (self, o) {
            (Key::Line(a), Key::Line(b)) => a.cmp(b),
            (Key::Arc(a), Key::Arc(b)) => a.compare(b),
            (Key::Spline(a), Key::Spline(b)) => a.partial_cmp(b),
            _ => None,
        }
    }
}

/// A spline's key at its curve parameter `t`.
fn spline_key(span: &Span, t: f64) -> Key {
    Key::Spline(if span.is_reversed() { -t } else { t })
}

fn exact2(p: Point2) -> [I; 2] {
    [I::exact(q(p.x)), I::exact(q(p.y))]
}

/// A line's `a x + b y + d` through `p` and `e`, exactly.
fn line_through(p: Point2, e: Point2) -> [R; 3] {
    let a = q(e.y) - q(p.y);
    let b = q(p.x) - q(e.x);
    let d = -(&a * q(p.x) + &b * q(p.y));
    [a, b, d]
}

/// A point's parameter along `p -> e`, enclosed.
fn line_param(p: Point2, e: Point2, x: &[I; 2]) -> Option<I> {
    let (dx, dy) = (q(e.x) - q(p.x), q(e.y) - q(p.y));
    let len2 = &dx * &dx + &dy * &dy;
    let t = x[0]
        .sub(&I::exact(q(p.x)))
        .mul(&I::exact(dx))
        .add(&x[1].sub(&I::exact(q(p.y))).mul(&I::exact(dy)));
    t.div(&I::exact(len2))
}

/// The arrangement being built.
struct Arrangement {
    tol: f64,
    segs: Vec<Seg>,
    /// Every point's rounded position.
    positions: BTreeMap<PId, Point2>,
    /// A vertex of `B` equal to one of `A`.
    alias: BTreeMap<PId, PId>,
    /// Places where each segment is cut.
    events: Vec<Vec<(Key, PId)>>,
    crosses: usize,
    /// The two segments meeting at each crossing.
    cross_segs: BTreeMap<usize, [(Operand, usize, usize); 2]>,
    /// The spline segments' spans (S9a.2).
    splines: Vec<Span>,
}

impl Arrangement {
    fn canon(&self, p: PId) -> PId {
        *self.alias.get(&p).unwrap_or(&p)
    }

    fn new_cross(&mut self, x: Point2) -> PId {
        let id = PId::Cross(self.crosses);
        self.crosses += 1;
        self.positions.insert(id, x);
        id
    }

    /// The place of an exact point strictly inside segment `s`, or `None`
    /// outside it; `Err` when undecided.
    fn place(&self, s: usize, x: &[I; 2]) -> Result<Option<Key>> {
        let seg = &self.segs[s];
        match seg.shape {
            Shape::Line => {
                let t = line_param(seg.p, seg.e, x)
                    .ok_or(Error::Degenerate("a zero-length segment"))?;
                let (lo, hi) = (t.lo().clone(), t.hi().clone());
                if lo > zero() && hi < R::from_integer(1.into()) {
                    Ok(Some(Key::Line(t)))
                } else if hi < zero() || lo > R::from_integer(1.into()) {
                    Ok(None)
                } else {
                    Err(Error::ComputationLimit("a crossing at a segment's end"))
                }
            }
            Shape::Arc { center, ccw, .. } => {
                Ok(within_arc(center, seg.p, seg.e, ccw, x)?.map(Key::Arc))
            }
            Shape::Circle { center, .. } => {
                let c = exact2(center);
                let r = [x[0].sub(&c[0]), x[1].sub(&c[1])];
                Ok(Some(Key::Arc(ArcPos::of(
                    &[I::exact(R::from_integer(1.into())), I::exact(zero())],
                    r,
                    1,
                )?)))
            }
            // A point within the resolution of a spline, at its nearest
            // parameter strictly inside the segment.
            Shape::Spline(k) => {
                let span = &self.splines[k];
                let [lo, hi] = span.range();
                Ok(splines::nearest(span, rounded(x), self.tol)?
                    .filter(|t| lo < *t && *t < hi)
                    .map(|t| spline_key(span, t)))
            }
        }
    }

    /// Records a meeting of segments `s` and `t` at the exact point `x`:
    /// identified with a stored end of either within the resolution, a cut
    /// of each segment it lies strictly inside.
    fn meet(&mut self, s: usize, t: usize, x: &[I; 2]) -> Result<()> {
        self.meet_keyed(s, t, x, &[])
    }

    /// [`Arrangement::meet`] with the keys of the segments known where the
    /// meeting is a new crossing (a spline's parameter there).
    fn meet_keyed(&mut self, s: usize, t: usize, x: &[I; 2], known: &[(usize, Key)]) -> Result<()> {
        let rx = rounded(x);
        let ends =
            |seg: &Seg| -> Vec<PId> { seg.from.iter().chain(seg.to.iter()).copied().collect() };
        let (ends_s, ends_t) = (ends(&self.segs[s]), ends(&self.segs[t]));
        // A stored end within the resolution: the meeting is that vertex.
        let near: Vec<PId> = ends_s
            .iter()
            .chain(&ends_t)
            .map(|p| self.canon(*p))
            .filter(|p| self.positions[p].distance(rx) <= self.tol)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if near.len() > 1 {
            return Err(Error::Degenerate(
                "a crossing within the resolution of two vertices",
            ));
        }
        let id = match near.first() {
            Some(v) => *v,
            None => {
                let id = self.new_cross(rx);
                if let PId::Cross(k) = id {
                    let key = |x: &Seg| (x.op, x.b, x.j);
                    let (a, b) = (key(&self.segs[s]), key(&self.segs[t]));
                    self.cross_segs
                        .insert(k, if a <= b { [a, b] } else { [b, a] });
                }
                id
            }
        };
        let exact_of = |id: PId, a: &Arrangement| match id {
            PId::Vertex(..) => exact2(a.positions[&id]),
            PId::Cross(_) => x.clone(),
        };
        for seg in [s, t] {
            let own: Vec<PId> = ends(&self.segs[seg])
                .iter()
                .map(|p| self.canon(*p))
                .collect();
            if own.contains(&id) {
                continue;
            }
            let given = matches!(id, PId::Cross(_))
                .then(|| {
                    known
                        .iter()
                        .find(|(k, _)| *k == seg)
                        .map(|(_, key)| key.clone())
                })
                .flatten();
            let placed = match given {
                Some(key) => Some(key),
                None => self.place(seg, &exact_of(id, self))?,
            };
            match placed {
                Some(key) => self.events[seg].push((key, id)),
                // A vertex snapped onto a segment it lies just past.
                None if matches!(id, PId::Vertex(..)) => {
                    return Err(Error::Degenerate(
                        "a vertex within the resolution of a segment's end",
                    ))
                }
                None => {}
            }
        }
        Ok(())
    }

    /// Two lines: their rational crossing, or their overlap cut at each
    /// other's ends.
    fn lines(&mut self, s: usize, t: usize) -> Result<Vec<(usize, usize)>> {
        let (a, b) = (&self.segs[s], &self.segs[t]);
        let (p, e, r, u) = (a.p, a.e, b.p, b.e);
        let v = |x: Point2| [q(x.x), q(x.y)];
        let (p, e, r, u) = (v(p), v(e), v(r), v(u));
        let sub = |x: &[R; 2], y: &[R; 2]| [&x[0] - &y[0], &x[1] - &y[1]];
        let cross = |x: &[R; 2], y: &[R; 2]| &x[0] * &y[1] - &x[1] * &y[0];
        let dot = |x: &[R; 2], y: &[R; 2]| &x[0] * &y[0] + &x[1] * &y[1];
        let (d1, d2) = (sub(&e, &p), sub(&u, &r));
        let den = cross(&d1, &d2);
        let rp = sub(&r, &p);
        if den != zero() {
            let ta = cross(&rp, &d2) / &den;
            let tb = cross(&rp, &d1) / &den;
            let one = R::from_integer(1.into());
            if ta < zero() || ta > one || tb < zero() || tb > one {
                return Ok(Vec::new());
            }
            let x = [
                I::exact(&p[0] + &ta * &d1[0]),
                I::exact(&p[1] + &ta * &d1[1]),
            ];
            self.meet(s, t, &x)?;
            return Ok(Vec::new());
        }
        if cross(&rp, &d1) != zero() {
            return Ok(Vec::new());
        }
        // Collinear: each one's ends strictly inside the other cut it.
        for (on, other) in [(s, t), (t, s)] {
            let (ends, pts) = {
                let o = &self.segs[other];
                ([o.from, o.to], [o.p, o.e])
            };
            for (end, pt) in ends.iter().zip(pts) {
                let Some(end) = end else { continue };
                let id = self.canon(*end);
                let seg = &self.segs[on];
                let (sp, se) = (v(seg.p), v(seg.e));
                let dd = sub(&se, &sp);
                let tt = dot(&sub(&v(pt), &sp), &dd) / dot(&dd, &dd);
                if tt > zero() && tt < R::from_integer(1.into()) {
                    let own = [seg.from, seg.to].map(|x| x.map(|x| self.canon(x)));
                    if !own.contains(&Some(id)) {
                        self.events[on].push((Key::Line(I::exact(tt)), id));
                    }
                }
            }
        }
        Ok(vec![(s, t)])
    }

    /// A spline segment and any other (S9a.2): its crossings with a line
    /// (S8b.3's roots), a circle or another spline; one curve shared with
    /// another spline segment, each one's ends inside the other cutting it.
    fn with_spline(&mut self, s: usize, t: usize) -> Result<Vec<(usize, usize)>> {
        let (sp, other) = if matches!(self.segs[s].shape, Shape::Spline(_)) {
            (s, t)
        } else {
            (t, s)
        };
        let Shape::Spline(k) = self.segs[sp].shape else {
            unreachable!("a spline")
        };
        let span = self.splines[k].clone();
        let exact = |e: &[R; 2]| [I::exact(e[0].clone()), I::exact(e[1].clone())];
        match self.segs[other].shape {
            Shape::Line => {
                let [a, b, d] = line_through(self.segs[other].p, self.segs[other].e);
                let m = crate::solid::split::spline::meets(&span, [&a, &b, &d])?;
                for c in m.crossings {
                    self.meet_keyed(sp, other, &exact(&c.exact), &[(sp, spline_key(&span, c.t))])?;
                }
            }
            Shape::Arc { center, radius, .. } | Shape::Circle { center, radius } => {
                let m = splines::with_circle(&span, center, radius)?;
                for c in m.crossings {
                    self.meet_keyed(sp, other, &exact(&c.exact), &[(sp, spline_key(&span, c.t))])?;
                }
            }
            Shape::Spline(k2) => {
                let span2 = self.splines[k2].clone();
                if splines::same_curve(&span, &span2) {
                    for (on, from) in [(sp, other), (other, sp)] {
                        let ends = [self.segs[from].from, self.segs[from].to];
                        for end in ends.iter().flatten() {
                            let id = self.canon(*end);
                            let own = [self.segs[on].from, self.segs[on].to]
                                .map(|x| x.map(|x| self.canon(x)));
                            if own.contains(&Some(id)) {
                                continue;
                            }
                            if let Some(key) = self.place(on, &exact2(self.positions[&id]))? {
                                self.events[on].push((key, id));
                            }
                        }
                    }
                    return Ok(vec![(sp.min(other), sp.max(other))]);
                }
                for c in splines::with_spline(&span, &span2)? {
                    self.meet_keyed(
                        sp,
                        other,
                        &exact(&c.exact),
                        &[
                            (sp, spline_key(&span, c.t[0])),
                            (other, spline_key(&span2, c.t[1])),
                        ],
                    )?;
                }
            }
        }
        Ok(Vec::new())
    }

    /// A line `s` and a circle or arc `t`.
    fn line_circle(&mut self, s: usize, t: usize) -> Result<()> {
        let (p, e) = (self.segs[s].p, self.segs[s].e);
        let (center, radius) = match self.segs[t].shape {
            Shape::Arc { center, radius, .. } | Shape::Circle { center, radius } => {
                (center, radius)
            }
            Shape::Line | Shape::Spline(_) => unreachable!("a curve"),
        };
        let line = line_through(p, e);
        self.circle_line(s, t, center, radius, &line)
    }

    /// The meetings of circle `(center, radius)` (segment `t`'s, or `s`'s
    /// for two circles) with `line`, recorded on `s` and `t`.
    fn circle_line(
        &mut self,
        s: usize,
        t: usize,
        center: Point2,
        radius: f64,
        line: &[R; 3],
    ) -> Result<()> {
        let [a, b, d] = line;
        let fc = a * q(center.x) + b * q(center.y) + d;
        let ab2 = a * a + b * b;
        let delta = &fc * &fc - q(radius) * q(radius) * &ab2;
        // A tangency cuts nothing: the pieces on either side of it lie on
        // one side of the other boundary (a result touching itself there is
        // refused when it is traced or validated).
        let points: Vec<[I; 2]> = match delta.cmp(&zero()) {
            Ordering::Greater | Ordering::Equal => return Ok(()),
            Ordering::Less => circle_points(center, radius, &[a.clone(), b.clone(), d.clone()])?,
        };
        for x in points {
            // On both segments (strictly, or at an end within the
            // resolution): recorded.
            let on = |a: &Arrangement, k: usize| -> Result<bool> {
                let seg = &a.segs[k];
                let at_end = [seg.from, seg.to]
                    .iter()
                    .flatten()
                    .any(|v| a.positions[&a.canon(*v)].distance(rounded(&x)) <= a.tol);
                Ok(at_end || a.place(k, &x)?.is_some())
            };
            if on(self, s)? && on(self, t)? {
                self.meet(s, t, &x)?;
            }
        }
        Ok(())
    }

    /// Two circles or arcs: one circle's pieces overlapping, or their
    /// radical line's meetings.
    fn circles(&mut self, s: usize, t: usize) -> Result<Vec<(usize, usize)>> {
        let geo = |seg: &Seg| match seg.shape {
            Shape::Arc { center, radius, .. } | Shape::Circle { center, radius } => {
                (center, radius)
            }
            Shape::Line | Shape::Spline(_) => unreachable!("a curve"),
        };
        let ((c1, r1), (c2, r2)) = (geo(&self.segs[s]), geo(&self.segs[t]));
        if c1 == c2 && r1 == r2 {
            // One circle: each one's ends strictly inside the other cut it.
            for (on, other) in [(s, t), (t, s)] {
                let ends = [self.segs[other].from, self.segs[other].to];
                for end in ends.iter().flatten() {
                    let id = self.canon(*end);
                    let own =
                        [self.segs[on].from, self.segs[on].to].map(|x| x.map(|x| self.canon(x)));
                    if own.contains(&Some(id)) {
                        continue;
                    }
                    if let Some(key) = self.place(on, &exact2(self.positions[&id]))? {
                        self.events[on].push((key, id));
                    }
                }
            }
            return Ok(vec![(s, t)]);
        }
        if c1 == c2 {
            return Ok(Vec::new());
        }
        let (x1, y1, x2, y2) = (q(c1.x), q(c1.y), q(c2.x), q(c2.y));
        let two = R::from_integer(2.into());
        let line = [
            &two * (&x2 - &x1),
            &two * (&y2 - &y1),
            (&x1 * &x1 + &y1 * &y1 - q(r1) * q(r1)) - (&x2 * &x2 + &y2 * &y2 - q(r2) * q(r2)),
        ];
        self.circle_line(s, t, c1, r1, &line)?;
        Ok(Vec::new())
    }
}

/// The stored segments of a set of profiles, boundaries numbered across
/// them (each profile's outer first), and which are holes.
fn segments(profiles: &[&Profile], op: Operand, splines: &mut Vec<Span>) -> (Vec<Seg>, Vec<bool>) {
    let mut all = Vec::new();
    let mut holes = Vec::new();
    for profile in profiles {
        let base = holes.len();
        for (k, _) in profile.boundaries().enumerate() {
            holes.push(k > 0);
        }
        for mut seg in segments_of(profile, op, splines) {
            seg.b += base;
            let shift = |id: Option<PId>| {
                id.map(|x| match x {
                    PId::Vertex(o, b, j) => PId::Vertex(o, b + base, j),
                    c => c,
                })
            };
            seg.from = shift(seg.from);
            seg.to = shift(seg.to);
            all.push(seg);
        }
    }
    (all, holes)
}

/// The stored segments of a profile.
fn segments_of(profile: &Profile, op: Operand, splines: &mut Vec<Span>) -> Vec<Seg> {
    let mut out = Vec::new();
    for (b, boundary) in profile.boundaries().enumerate() {
        match &boundary.kind {
            BoundaryKind::Circle { center, radius } => {
                let p = Point2::new(center.x + radius, center.y);
                out.push(Seg {
                    op,
                    b,
                    j: 0,
                    from: None,
                    to: None,
                    p,
                    e: p,
                    shape: Shape::Circle {
                        center: *center,
                        radius: *radius,
                    },
                });
            }
            BoundaryKind::Polygon(points) => {
                let n = points.len();
                for j in 0..n {
                    out.push(Seg {
                        op,
                        b,
                        j,
                        from: Some(PId::Vertex(op, b, j)),
                        to: Some(PId::Vertex(op, b, (j + 1) % n)),
                        p: points[j],
                        e: points[(j + 1) % n],
                        shape: Shape::Line,
                    });
                }
            }
            BoundaryKind::Path { points, segments } => {
                let n = points.len();
                for (j, s) in segments.iter().enumerate() {
                    let shape = match s {
                        Segment::Line => Shape::Line,
                        Segment::Arc {
                            center,
                            radius,
                            ccw,
                        } => Shape::Arc {
                            center: *center,
                            radius: *radius,
                            ccw: *ccw,
                        },
                        Segment::Spline(span) => {
                            splines.push(span.clone());
                            Shape::Spline(splines.len() - 1)
                        }
                    };
                    out.push(Seg {
                        op,
                        b,
                        j,
                        from: Some(PId::Vertex(op, b, j)),
                        to: Some(PId::Vertex(op, b, (j + 1) % n)),
                        p: points[j],
                        e: points[(j + 1) % n],
                        shape,
                    });
                }
            }
        }
    }
    out
}

/// A piece of a stored segment between two arrangement points, in the
/// segment's stored direction.
#[derive(Debug, Clone)]
struct Piece {
    seg: usize,
    part: usize,
    whole: bool,
    from: Option<PId>,
    to: Option<PId>,
    p: Point2,
    e: Point2,
    /// A spline piece's curve parameters at its start and end (S9a.2).
    range: Option<(f64, f64)>,
}

/// A directed result edge, the region on its left.
#[derive(Debug, Clone)]
struct Edge {
    from: Option<PId>,
    to: Option<PId>,
    p: Point2,
    e: Point2,
    shape: Shape,
    origins: Vec<Origin2>,
    /// The arrangement pieces it joins, in order, each with whether it runs
    /// along its segment's stored direction.
    parts: Vec<(usize, bool)>,
    /// A spline edge's curve parameters at its start and end.
    range: Option<(f64, f64)>,
}

impl Piece {
    /// An off-centre point of the piece (fraction `0.4453125` along it).
    fn sample(&self, seg: &Seg, splines: &[Span]) -> Result<Point2> {
        let f = 0.4453125;
        Ok(match seg.shape {
            Shape::Spline(k) => {
                let (t0, t1) = self.range.expect("a spline piece's range");
                splines[k].curve().point(t0 + (t1 - t0) * f)?
            }
            Shape::Line => Point2::new(
                self.p.x + (self.e.x - self.p.x) * f,
                self.p.y + (self.e.y - self.p.y) * f,
            ),
            Shape::Arc {
                center,
                radius,
                ccw,
            } => {
                let a0 = (self.p.y - center.y).atan2(self.p.x - center.x);
                let sweep = super::arc_sweep(center, self.p, self.e, ccw);
                let a = a0 + sweep * f;
                Point2::new(center.x + radius * a.cos(), center.y + radius * a.sin())
            }
            Shape::Circle { center, radius } => {
                if self.from.is_none() {
                    let a = std::f64::consts::TAU * f;
                    return Ok(Point2::new(
                        center.x + radius * a.cos(),
                        center.y + radius * a.sin(),
                    ));
                }
                let a0 = (self.p.y - center.y).atan2(self.p.x - center.x);
                let sweep = super::arc_sweep(center, self.p, self.e, true);
                let a = a0 + sweep * f;
                Point2::new(center.x + radius * a.cos(), center.y + radius * a.sin())
            }
        })
    }
}

/// Two sets of profiles' boundaries arranged: their pieces, each shared
/// with the other set's or placed inside or outside it.
pub(crate) struct Arranged {
    arr: Arrangement,
    pieces: Vec<Piece>,
    partner: BTreeMap<usize, (usize, bool)>,
    place: Vec<Option<Location>>,
    /// Per operand, which of its boundaries are holes.
    holes: [Vec<bool>; 2],
}

/// The Boolean of two profiles in one frame: its result profiles with their
/// provenance, in a deterministic order (by their first segment's origin).
pub(crate) fn boolean(a: &Profile, b: &Profile, op: Op2) -> Result<Boolean2> {
    let arranged = arrange(&[a], &[b], a.tolerance())?;
    select_trace(&arranged, op, false, true, a.tolerance())
}

/// A piece's class against the other set (S9a.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Class {
    Inside,
    Outside,
    /// On the other set's boundary: the partner piece, and whether their
    /// region-left directions agree.
    Shared {
        partner: usize,
        agree: bool,
    },
}

/// A piece's curve, in its segment's stored direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Curve {
    Line,
    Arc {
        center: Point2,
        radius: f64,
        ccw: bool,
    },
    /// A whole circle (uncut), counter-clockwise from angle 0.
    Circle {
        center: Point2,
        radius: f64,
    },
    /// A spline piece: its spline (`Arranged::spline`) and curve
    /// parameters at its start and end.
    Spline {
        index: usize,
        range: (f64, f64),
    },
}

/// One piece of an arrangement, for S9a.2's stacks.
#[derive(Debug, Clone)]
pub(crate) struct PieceView {
    pub(crate) op: Operand,
    pub(crate) b: usize,
    pub(crate) j: usize,
    pub(crate) from: Option<PId>,
    pub(crate) to: Option<PId>,
    pub(crate) p: Point2,
    pub(crate) e: Point2,
    pub(crate) curve: Curve,
    pub(crate) class: Class,
    /// Whether it stands for its place: not the `B` side of a shared pair.
    pub(crate) representative: bool,
    /// Whether the region left and right of its stored direction lies in
    /// `A` and in `B`.
    pub(crate) left: [bool; 2],
    pub(crate) right: [bool; 2],
}

impl Arranged {
    /// Every piece, in segment order and then along each segment.
    pub(crate) fn pieces(&self) -> Vec<PieceView> {
        let hole = |o: Operand, b: usize| self.holes[if o == Operand::A { 0 } else { 1 }][b];
        self.pieces
            .iter()
            .enumerate()
            .map(|(i, piece)| {
                let seg = &self.arr.segs[piece.seg];
                let curve = match seg.shape {
                    Shape::Line => Curve::Line,
                    Shape::Arc {
                        center,
                        radius,
                        ccw,
                    } => Curve::Arc {
                        center,
                        radius,
                        ccw,
                    },
                    Shape::Circle { center, radius } if piece.from.is_some() => Curve::Arc {
                        center,
                        radius,
                        ccw: true,
                    },
                    Shape::Circle { center, radius } => Curve::Circle { center, radius },
                    Shape::Spline(index) => Curve::Spline {
                        index,
                        range: piece.range.expect("a spline piece's range"),
                    },
                };
                let class = match self.partner.get(&i) {
                    Some(&(k, same)) => {
                        let other = &self.arr.segs[self.pieces[k].seg];
                        Class::Shared {
                            partner: k,
                            agree: same == (hole(seg.op, seg.b) == hole(other.op, other.b)),
                        }
                    }
                    None => match self.place[i] {
                        Some(Location::Inside) => Class::Inside,
                        _ => Class::Outside,
                    },
                };
                // Its own profile's material lies left of a stored outer
                // boundary (counter-clockwise), right of a hole; the
                // other's on both sides, neither, or as its partner's.
                let own = !hole(seg.op, seg.b);
                let (other_left, other_right) = match class {
                    Class::Inside => (true, true),
                    Class::Outside => (false, false),
                    Class::Shared { agree, .. } => (own == agree, own != agree),
                };
                let (left, right) = if seg.op == Operand::A {
                    ([own, other_left], [!own, other_right])
                } else {
                    ([other_left, own], [other_right, !own])
                };
                PieceView {
                    op: seg.op,
                    b: seg.b,
                    j: seg.j,
                    from: piece.from,
                    to: piece.to,
                    p: piece.p,
                    e: piece.e,
                    curve,
                    class,
                    representative: seg.op == Operand::A || !self.partner.contains_key(&i),
                    left,
                    right,
                }
            })
            .collect()
    }
    /// A spline piece's spline.
    pub(crate) fn spline(&self, index: usize) -> &Span {
        &self.arr.splines[index]
    }
    /// A point's canonical id (a vertex of `B` equal to one of `A`'s is
    /// `A`'s).
    pub(crate) fn canon(&self, p: PId) -> PId {
        self.arr.canon(p)
    }
    /// A point's rounded position.
    pub(crate) fn position(&self, p: PId) -> Point2 {
        self.arr.positions[&self.arr.canon(p)]
    }
}

/// Where a point lies against a set of disjoint profiles.
fn classify_set(set: &[&Profile], point: Point2) -> Result<Location> {
    for p in set {
        match p.classify(point)? {
            Location::Outside => {}
            at => return Ok(at),
        }
    }
    Ok(Location::Outside)
}

/// Arranges two sets of profiles' boundaries.
pub(crate) fn arrange(a: &[&Profile], b: &[&Profile], tolerance: Tolerance) -> Result<Arranged> {
    let tol = tolerance.linear();
    let mut splines: Vec<Span> = Vec::new();
    let (mut segs, holes_a) = segments(a, Operand::A, &mut splines);
    let (segs_b, holes_b) = segments(b, Operand::B, &mut splines);
    segs.extend(segs_b);
    let mut positions = BTreeMap::new();
    for seg in &segs {
        for (id, pt) in [(seg.from, seg.p), (seg.to, seg.e)] {
            if let Some(id) = id {
                positions.insert(id, pt);
            }
        }
    }
    // B's vertices equal to A's are A's.
    let mut alias = BTreeMap::new();
    for (id, pt) in &positions {
        if let PId::Vertex(Operand::B, ..) = id {
            if let Some((a_id, _)) = positions
                .iter()
                .find(|(k, v)| matches!(k, PId::Vertex(Operand::A, ..)) && **v == *pt)
            {
                alias.insert(*id, *a_id);
            }
        }
    }
    let n = segs.len();
    let mut arr = Arrangement {
        tol,
        segs,
        positions,
        alias,
        events: vec![Vec::new(); n],
        crosses: 0,
        cross_segs: BTreeMap::new(),
        splines,
    };
    // Stored vertices of one within the resolution of the other's (not
    // equal): sub-resolution.
    let ids: Vec<(PId, Point2)> = arr.positions.iter().map(|(k, v)| (*k, *v)).collect();
    for (i, (ka, pa)) in ids.iter().enumerate() {
        for (kb, pb) in &ids[i + 1..] {
            let ops = |k: &PId| match k {
                PId::Vertex(o, ..) => Some(*o),
                PId::Cross(_) => None,
            };
            if ops(ka) != ops(kb) && pa != pb && pa.distance(*pb) <= tol {
                return Err(Error::Degenerate(
                    "two profiles' vertices within the resolution",
                ));
            }
        }
    }
    // Every pair of segments of different operands.
    let mut overlaps: BTreeSet<(usize, usize)> = BTreeSet::new();
    for s in 0..n {
        for t in 0..n {
            if arr.segs[s].op != Operand::A || arr.segs[t].op != Operand::B {
                continue;
            }
            let spline = |x: &Seg| matches!(x.shape, Shape::Spline(_));
            let lines = |x: &Seg| matches!(x.shape, Shape::Line);
            let found = if spline(&arr.segs[s]) || spline(&arr.segs[t]) {
                arr.with_spline(s, t)?
            } else {
                match (lines(&arr.segs[s]), lines(&arr.segs[t])) {
                    (true, true) => arr.lines(s, t)?,
                    (true, false) => {
                        arr.line_circle(s, t)?;
                        Vec::new()
                    }
                    (false, true) => {
                        arr.line_circle(t, s)?;
                        Vec::new()
                    }
                    (false, false) => arr.circles(s, t)?,
                }
            };
            overlaps.extend(found);
        }
    }
    // Every stored vertex of one on a segment of the other (exactly on a
    // line, or within the resolution of a curve) cuts it.
    for s in 0..n {
        let op = arr.segs[s].op;
        let verts: Vec<PId> = arr
            .positions
            .keys()
            .copied()
            .filter(|k| matches!(k, PId::Vertex(o, ..) if *o != op))
            .map(|k| arr.canon(k))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        for v in verts {
            let seg = &arr.segs[s];
            let own = [seg.from, seg.to].map(|x| x.map(|x| arr.canon(x)));
            if own.contains(&Some(v)) || arr.events[s].iter().any(|(_, id)| *id == v) {
                continue;
            }
            let pt = arr.positions[&v];
            let on = match seg.shape {
                Shape::Line => {
                    let [a, b, d] = line_through(seg.p, seg.e);
                    &a * q(pt.x) + &b * q(pt.y) + &d == zero()
                }
                Shape::Arc { center, radius, .. } | Shape::Circle { center, radius } => {
                    (pt.distance(center) - radius).abs() <= tol
                }
                Shape::Spline(k) => splines::nearest(&arr.splines[k], pt, tol)?.is_some(),
            };
            if on {
                if let Some(key) = arr.place(s, &exact2(pt))? {
                    arr.events[s].push((key, v));
                }
            }
        }
    }
    // Pieces along each segment.
    let mut pieces: Vec<Piece> = Vec::new();
    for s in 0..n {
        let mut ev = std::mem::take(&mut arr.events[s]);
        let mut failed = false;
        ev.sort_by(|x, y| {
            x.0.compare(&y.0).unwrap_or_else(|| {
                failed = true;
                Ordering::Equal
            })
        });
        if failed {
            return Err(Error::ComputationLimit("two cuts of a segment in order"));
        }
        ev.dedup_by(|x, y| x.1 == y.1);
        let seg = &arr.segs[s];
        let mut chain: Vec<Option<PId>> = Vec::new();
        match seg.shape {
            // A circle cut at one point only (another's vertex on it) stays
            // whole.
            Shape::Circle { .. } if ev.len() < 2 => {
                pieces.push(Piece {
                    seg: s,
                    part: 0,
                    whole: true,
                    from: None,
                    to: None,
                    p: seg.p,
                    e: seg.p,
                    range: None,
                });
                continue;
            }
            Shape::Circle { .. } => {
                for (_, id) in &ev {
                    chain.push(Some(*id));
                }
                chain.push(Some(ev[0].1));
            }
            _ => {
                chain.push(seg.from.map(|x| arr.canon(x)));
                for (_, id) in &ev {
                    chain.push(Some(*id));
                }
                chain.push(seg.to.map(|x| arr.canon(x)));
            }
        }
        let parts = chain.len() - 1;
        // A spline's curve parameters at its cuts, along its stored
        // direction.
        let params: Option<Vec<f64>> = match seg.shape {
            Shape::Spline(k) => {
                let span = &arr.splines[k];
                let [lo, hi] = span.range();
                let back = span.is_reversed();
                let mut v = vec![if back { hi } else { lo }];
                for (key, _) in &ev {
                    if let Key::Spline(x) = key {
                        v.push(if back { -x } else { *x });
                    }
                }
                v.push(if back { lo } else { hi });
                Some(v)
            }
            _ => None,
        };
        for k in 0..parts {
            let (f, t) = (chain[k], chain[k + 1]);
            let (pf, pt) = (arr.positions[&f.unwrap()], arr.positions[&t.unwrap()]);
            if f == t || pf.distance(pt) <= tol {
                return Err(Error::Degenerate(
                    "a piece of a profile thinner than the resolution",
                ));
            }
            pieces.push(Piece {
                seg: s,
                part: k,
                whole: parts == 1 && !matches!(seg.shape, Shape::Circle { .. }),
                from: f,
                to: t,
                p: pf,
                e: pt,
                range: params.as_ref().map(|v| (v[k], v[k + 1])),
            });
        }
    }
    // Shared pieces: the same ends on one line or one circle.
    let same_support = |x: usize, y: usize| overlaps.contains(&(x.min(y), x.max(y)));
    let mut partner: BTreeMap<usize, (usize, bool)> = BTreeMap::new();
    for (i, pa) in pieces.iter().enumerate() {
        if arr.segs[pa.seg].op != Operand::A {
            continue;
        }
        for (k, pb) in pieces.iter().enumerate() {
            if arr.segs[pb.seg].op != Operand::B || !same_support(pa.seg, pb.seg) {
                continue;
            }
            let same = pa.from == pb.from && pa.to == pb.to;
            let opposite = pa.from == pb.to && pa.to == pb.from;
            let sense = |x: &Piece| match arr.segs[x.seg].shape {
                Shape::Arc { ccw, .. } => Some(ccw),
                _ => Some(true),
            };
            let shared = match (&arr.segs[pa.seg].shape, &arr.segs[pb.seg].shape) {
                (Shape::Line, Shape::Line) | (Shape::Spline(_), Shape::Spline(_)) => {
                    same || opposite
                }
                (
                    Shape::Arc { .. } | Shape::Circle { .. },
                    Shape::Arc { .. } | Shape::Circle { .. },
                ) => {
                    let (ca, cb) = (sense(pa), sense(pb));
                    (same && ca == cb) || (opposite && ca != cb && pa.from.is_some())
                }
                _ => false,
            };
            if shared {
                partner.insert(i, (k, same));
                partner.insert(k, (i, same));
            }
        }
    }
    // Each piece's place against the other profile.
    let mut place: Vec<Option<Location>> = Vec::new();
    for (i, piece) in pieces.iter().enumerate() {
        if partner.contains_key(&i) {
            place.push(None);
            continue;
        }
        let other = if arr.segs[piece.seg].op == Operand::A {
            b
        } else {
            a
        };
        let at = classify_set(other, piece.sample(&arr.segs[piece.seg], &arr.splines)?)?;
        if at == Location::Boundary {
            return Err(Error::Degenerate(
                "a boundary within the resolution of the other profile's",
            ));
        }
        place.push(Some(at));
    }
    Ok(Arranged {
        arr,
        pieces,
        partner,
        place,
        holes: [holes_a, holes_b],
    })
}

/// The operation's pieces of an arrangement, traced into its result
/// profiles.
///
/// With `swap` the operands' roles are exchanged (a cut keeps `B - A`);
/// with `join` a cycle keeps no vertex where it does not turn.
pub(crate) fn select_trace(
    arranged: &Arranged,
    op: Op2,
    swap: bool,
    join: bool,
    tolerance: Tolerance,
) -> Result<Boolean2> {
    let Arranged {
        arr,
        pieces,
        partner,
        place,
        holes,
    } = arranged;
    let hole = |o: Operand, b: usize| holes[if o == Operand::A { 0 } else { 1 }][b];
    // The kept pieces, each reversed when the region lies on its right: a
    // hole's stored direction reversed, and a kept tool's once more.
    let mut kept: Vec<(usize, bool)> = Vec::new();
    for (i, piece) in pieces.iter().enumerate() {
        let seg = &arr.segs[piece.seg];
        let (inside, outside) = (Some(Location::Inside), Some(Location::Outside));
        // The first operand of the rule: A, or B when swapped.
        let first = if swap { Operand::B } else { Operand::A };
        let role = if seg.op == first {
            Operand::A
        } else {
            Operand::B
        };
        let keep_flip: Option<bool> = match partner.get(&i) {
            Some(&(k, same)) => {
                // In region-left directions: a hole's reversed.
                let other = &arr.segs[pieces[k].seg];
                let agree = same == (hole(seg.op, seg.b) == hole(other.op, other.b));
                match (op, role, agree) {
                    (Op2::Fuse | Op2::Common, Operand::A, true) => Some(false),
                    (Op2::Cut, Operand::A, false) => Some(false),
                    _ => None,
                }
            }
            None => match (op, role, place[i]) {
                (Op2::Fuse, _, p) if p == outside => Some(false),
                (Op2::Common, _, p) if p == inside => Some(false),
                (Op2::Cut, Operand::A, p) if p == outside => Some(false),
                (Op2::Cut, Operand::B, p) if p == inside => Some(true),
                _ => None,
            },
        };
        if let Some(flip) = keep_flip {
            kept.push((i, hole(seg.op, seg.b) != flip));
        }
    }
    traced(arranged, &kept, join, tolerance)
}

/// The region of a set function of the operands (`keep(in A, in B)`),
/// traced into profiles (S9a.2's caps between slabs): the pieces with the
/// region on one side only, a piece shared by both boundaries once (as
/// `A`'s).
pub(crate) fn select_with(
    arranged: &Arranged,
    keep: &dyn Fn(bool, bool) -> bool,
    join: bool,
    tolerance: Tolerance,
) -> Result<Boolean2> {
    let mut kept: Vec<(usize, bool)> = Vec::new();
    for (i, view) in arranged.pieces().iter().enumerate() {
        if !view.representative {
            continue;
        }
        let (l, r) = (
            keep(view.left[0], view.left[1]),
            keep(view.right[0], view.right[1]),
        );
        if l != r {
            kept.push((i, r));
        }
    }
    traced(arranged, &kept, join, tolerance)
}

/// Kept pieces `(piece, reversed)`, the region on the left of each as
/// traversed, traced into profiles.
fn traced(
    arranged: &Arranged,
    kept: &[(usize, bool)],
    join: bool,
    tolerance: Tolerance,
) -> Result<Boolean2> {
    let tol = tolerance.linear();
    let Arranged {
        arr,
        pieces,
        partner,
        ..
    } = arranged;
    let seg_ref = |piece: &Piece| -> SegRef {
        let seg = &arr.segs[piece.seg];
        (seg.op, seg.b, seg.j, (!piece.whole).then_some(piece.part))
    };
    let mut edges: Vec<Edge> = Vec::new();
    // Whole circles keep a direction flag in their origin's order: a
    // reversed whole circle is a clockwise cycle.
    let mut circle_cw: BTreeSet<usize> = BTreeSet::new();
    for &(i, reverse) in kept {
        let piece = &pieces[i];
        let seg = &arr.segs[piece.seg];
        let origin = Origin2 {
            from: seg_ref(piece),
            shared: partner.get(&i).map(|&(k, _)| seg_ref(&pieces[k])),
        };
        let shape = match seg.shape {
            Shape::Arc {
                center,
                radius,
                ccw,
            } => Shape::Arc {
                center,
                radius,
                ccw: ccw != reverse,
            },
            Shape::Circle { center, radius } if piece.from.is_some() => Shape::Arc {
                center,
                radius,
                ccw: !reverse,
            },
            other => other,
        };
        let (from, to, p, e) = if reverse {
            (piece.to, piece.from, piece.e, piece.p)
        } else {
            (piece.from, piece.to, piece.p, piece.e)
        };
        if piece.from.is_none() && reverse {
            circle_cw.insert(edges.len());
        }
        edges.push(Edge {
            from,
            to,
            p,
            e,
            shape,
            origins: vec![origin],
            parts: vec![(i, !reverse)],
            range: piece
                .range
                .map(|(a, b)| if reverse { (b, a) } else { (a, b) }),
        });
    }
    // Cycles: each kept piece's end starts exactly one kept piece.
    let mut starts: BTreeMap<PId, usize> = BTreeMap::new();
    for (i, e) in edges.iter().enumerate() {
        if let Some(f) = e.from {
            if starts.insert(f, i).is_some() {
                return Err(Error::Degenerate("a result touching itself at a point"));
            }
        }
    }
    let mut used = vec![false; edges.len()];
    let mut cycles: Vec<(Vec<usize>, bool)> = Vec::new();
    for i in 0..edges.len() {
        if used[i] {
            continue;
        }
        used[i] = true;
        if edges[i].from.is_none() {
            cycles.push((vec![i], !circle_cw.contains(&i)));
            continue;
        }
        let mut cycle = vec![i];
        let mut at = edges[i].to;
        while at != edges[i].from {
            let next = at
                .and_then(|p| starts.get(&p))
                .copied()
                .ok_or(Error::InvalidTopology("an open Boolean cycle"))?;
            if used[next] {
                return Err(Error::InvalidTopology("a Boolean cycle does not close"));
            }
            used[next] = true;
            cycle.push(next);
            at = edges[next].to;
        }
        let area = signed_area(
            &cycle.iter().map(|&k| &edges[k]).collect::<Vec<_>>(),
            &arr.splines,
        );
        cycles.push((cycle, area > 0.0));
    }
    // Outer boundaries and their holes, with no vertex where a boundary
    // does not turn when joined.
    type Cycle = (Traced, Vec<(usize, bool)>);
    let mut outers: Vec<Cycle> = Vec::new();
    let mut holes: Vec<Cycle> = Vec::new();
    for (cycle, ccw) in &cycles {
        let raw: Vec<Edge> = cycle.iter().map(|&k| edges[k].clone()).collect();
        let parts: Vec<(usize, bool)> = raw.iter().flat_map(|e| e.parts.clone()).collect();
        let merged = if join { merged(raw) } else { raw };
        let built = boundary_of(
            &merged.iter().collect::<Vec<_>>(),
            !ccw,
            &arr.splines,
            tolerance,
        )?;
        if *ccw {
            outers.push((built, parts));
        } else {
            holes.push((built, parts));
        }
    }
    let mut owned: Vec<Vec<Cycle>> = vec![Vec::new(); outers.len()];
    for hole in holes {
        let sample = hole.0 .0.sample();
        let owner = outers
            .iter()
            .position(|((o, ..), _)| o.locate(sample, tolerance) == Location::Inside)
            .ok_or(Error::Degenerate("a result's hole outside its boundary"))?;
        owned[owner].push(hole);
    }
    let mut out = Vec::new();
    for (((outer, segs, verts), parts), hs) in outers.into_iter().zip(owned) {
        let mut segments = vec![segs];
        let mut vertices = vec![verts];
        let mut cycles = vec![parts];
        let mut hole_boundaries = Vec::new();
        for ((h, s, v), c) in hs {
            hole_boundaries.push(h);
            segments.push(s);
            vertices.push(v);
            cycles.push(c);
        }
        let profile = Profile::new(outer, hole_boundaries, tolerance).map_err(|e| match e {
            Error::InvalidHole(_) | Error::SelfIntersection | Error::IntersectingBoundaries(..) => {
                Error::Degenerate("a result thinner than the resolution")
            }
            e => e,
        })?;
        out.push(Piece2 {
            profile,
            segments,
            vertices,
            cycles,
        });
    }
    // Separate results must not touch.
    for i in 0..out.len() {
        for k in (i + 1)..out.len() {
            if boundaries_touch(out[i].profile.outer(), out[k].profile.outer(), tol) {
                return Err(Error::Degenerate("two results touching"));
            }
        }
    }
    out.sort_by_key(|p| p.segments[0].iter().flatten().min().copied());
    Ok(Boolean2 {
        pieces: out,
        crosses: arr.cross_segs.clone(),
    })
}

/// A spline piece with its traversal's ends set to the arrangement's
/// points (a crossing of two splines is one point; each restriction's end
/// is within rounding of it).
pub(crate) fn pinned(span: &Span, p: Point2, e: Point2) -> Result<Span> {
    let c = span.curve().as_curve3();
    let mut poles = span.curve().poles();
    let n = poles.len();
    let (first, last) = if span.is_reversed() { (e, p) } else { (p, e) };
    if poles[0] == first && poles[n - 1] == last {
        return Ok(span.clone());
    }
    poles[0] = first;
    poles[n - 1] = last;
    let curve = crate::BSplineCurve2::new(
        c.degree(),
        poles,
        None,
        c.knots().to_vec(),
        c.multiplicities().to_vec(),
    )?;
    let whole = crate::topology::SplineSpan::whole(curve);
    Ok(if span.is_reversed() {
        whole.reversed()
    } else {
        whole
    })
}

/// A traced boundary with its segments' and vertices' provenance.
type Traced = (Boundary, Vec<Vec<Origin2>>, Vec<PId>);

/// A cycle with consecutive collinear lines, and consecutive arcs of one
/// circle in one sense, joined (their shared vertex a point where the
/// boundary does not turn); a cycle left as one arc back to its start is
/// its whole circle.
fn merged(mut cycle: Vec<Edge>) -> Vec<Edge> {
    let joins = |a: &Edge, b: &Edge| match (a.shape, b.shape) {
        (Shape::Line, Shape::Line) => {
            let d = |e: &Edge| [q(e.e.x) - q(e.p.x), q(e.e.y) - q(e.p.y)];
            let (u, v) = (d(a), d(b));
            &u[0] * &v[1] - &u[1] * &v[0] == zero() && &u[0] * &v[0] + &u[1] * &v[1] > zero()
        }
        (
            Shape::Arc {
                center: c1,
                radius: r1,
                ccw: s1,
            },
            Shape::Arc {
                center: c2,
                radius: r2,
                ccw: s2,
            },
        ) => c1 == c2 && r1 == r2 && s1 == s2,
        // Consecutive pieces of one spline segment, running on (S9a.2).
        (Shape::Spline(k1), Shape::Spline(k2)) => {
            k1 == k2 && matches!((a.range, b.range), (Some(x), Some(y)) if x.1 == y.0)
        }
        _ => false,
    };
    loop {
        let n = cycle.len();
        if n < 2 {
            break;
        }
        let Some(i) = (0..n).find(|&i| joins(&cycle[i], &cycle[(i + 1) % n])) else {
            break;
        };
        let j = (i + 1) % n;
        let (a, b) = (cycle[i].clone(), cycle[j].clone());
        let mut origins = a.origins.clone();
        origins.extend(b.origins.iter().copied());
        let mut parts = a.parts.clone();
        parts.extend(b.parts.iter().copied());
        let joined = Edge {
            from: a.from,
            to: b.to,
            p: a.p,
            e: b.e,
            shape: a.shape,
            origins,
            parts,
            range: a.range.zip(b.range).map(|(x, y)| (x.0, y.1)),
        };
        // Keep the cycle's start where it was when the join wraps round.
        if j == 0 {
            cycle.remove(n - 1);
            cycle[0] = joined;
        } else {
            cycle[i] = joined;
            cycle.remove(j);
        }
    }
    if cycle.len() == 1 && cycle[0].from.is_some() && cycle[0].from == cycle[0].to {
        if let Shape::Arc { center, radius, .. } = cycle[0].shape {
            cycle[0].shape = Shape::Circle { center, radius };
            cycle[0].from = None;
            cycle[0].to = None;
        }
    }
    cycle
}

/// A cycle's signed area (binary64, for its orientation).
fn signed_area(cycle: &[&Edge], splines: &[Span]) -> f64 {
    let mut twice = 0.0;
    for e in cycle {
        match e.shape {
            Shape::Spline(k) => {
                twice += e.p.x * e.e.y - e.e.x * e.p.y;
                if let Some((t0, t1)) = e.range {
                    if let Ok(part) = crate::solid::split::spline::piece(&splines[k], t0, t1) {
                        twice += crate::solid::split::spline::twice_area_beyond_chord(&part);
                    }
                }
            }
            Shape::Circle { radius, .. } => {
                twice += std::f64::consts::TAU * radius * radius;
            }
            Shape::Line => twice += e.p.x * e.e.y - e.e.x * e.p.y,
            Shape::Arc {
                center,
                radius,
                ccw,
            } => {
                twice += e.p.x * e.e.y - e.e.x * e.p.y;
                let phi = super::arc_sweep(center, e.p, e.e, ccw);
                twice += radius * radius * (phi - phi.sin());
            }
        }
    }
    twice / 2.0
}

/// A cycle as a stored counter-clockwise boundary (a hole's traversal
/// reversed), with its segments' and vertices' provenance in stored order.
fn boundary_of(
    cycle: &[&Edge],
    hole: bool,
    splines: &[Span],
    tolerance: Tolerance,
) -> Result<Traced> {
    if cycle.len() == 1 && cycle[0].from.is_none() {
        let Shape::Circle { center, radius } = cycle[0].shape else {
            return Err(Error::InvalidTopology("a ring that is not a circle"));
        };
        return Ok((
            Boundary::circle(center, radius, tolerance)?,
            vec![cycle[0].origins.clone()],
            vec![],
        ));
    }
    let mut items: Vec<(PId, Point2, Segment, Vec<Origin2>)> = cycle
        .iter()
        .map(|e| {
            let seg = match e.shape {
                Shape::Spline(k) => {
                    let (t0, t1) = e.range.expect("a spline edge's range");
                    let part = crate::solid::split::spline::piece(&splines[k], t0, t1)?;
                    Segment::Spline(pinned(&part, e.p, e.e)?)
                }
                Shape::Line => Segment::Line,
                Shape::Arc {
                    center,
                    radius,
                    ccw,
                } => Segment::Arc {
                    center,
                    radius,
                    ccw,
                },
                Shape::Circle { .. } => unreachable!("a cut circle is arcs"),
            };
            Ok((e.from.expect("a vertex"), e.p, seg, e.origins.clone()))
        })
        .collect::<Result<_>>()?;
    if hole {
        let n = items.len();
        let mut rev = Vec::with_capacity(n);
        for k in 0..n {
            let (_, _, seg, origin) = items[n - 1 - k].clone();
            let (from, p) = (items[(n - k) % n].0, items[(n - k) % n].1);
            let seg = match seg {
                Segment::Arc {
                    center,
                    radius,
                    ccw,
                } => Segment::Arc {
                    center,
                    radius,
                    ccw: !ccw,
                },
                Segment::Spline(span) => Segment::Spline(span.reversed()),
                s => s,
            };
            rev.push((from, p, seg, origin));
        }
        items = rev;
    }
    let points: Vec<Point2> = items.iter().map(|i| i.1).collect();
    let segments: Vec<Segment> = items.iter().map(|i| i.2.clone()).collect();
    let all_lines = segments.iter().all(|s| matches!(s, Segment::Line));
    let boundary = if all_lines {
        Boundary::polygon(points.clone(), tolerance)
    } else {
        Boundary::path(points.clone(), segments, tolerance)
    }
    .map_err(|e| match e {
        Error::SelfIntersection => Error::Degenerate("a result thinner than the resolution"),
        e => e,
    })?;
    let stored: Vec<Point2> = match &boundary.kind {
        BoundaryKind::Polygon(p) => p.clone(),
        BoundaryKind::Path { points, .. } => points.clone(),
        BoundaryKind::Circle { .. } => Vec::new(),
    };
    if stored != points {
        return Err(Error::InvalidTopology(
            "a Boolean boundary stored in another order",
        ));
    }
    Ok((
        boundary,
        items.iter().map(|i| i.3.clone()).collect(),
        items.iter().map(|i| i.0).collect(),
    ))
}
