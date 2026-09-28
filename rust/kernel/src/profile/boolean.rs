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
use crate::certified::{Interval as I, Real};
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
/// segment, its provenance and its start point's.
pub(crate) struct Piece2 {
    pub(crate) profile: Profile,
    pub(crate) segments: Vec<Vec<Origin2>>,
    pub(crate) vertices: Vec<Vec<PId>>,
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
}

impl Key {
    fn compare(&self, o: &Self) -> Option<Ordering> {
        match (self, o) {
            (Key::Line(a), Key::Line(b)) => a.cmp(b),
            (Key::Arc(a), Key::Arc(b)) => a.compare(b),
            _ => None,
        }
    }
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
        }
    }

    /// Records a meeting of segments `s` and `t` at the exact point `x`:
    /// identified with a stored end of either within the resolution, a cut
    /// of each segment it lies strictly inside.
    fn meet(&mut self, s: usize, t: usize, x: &[I; 2]) -> Result<()> {
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
            match self.place(seg, &exact_of(id, self))? {
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

    /// A line `s` and a circle or arc `t`.
    fn line_circle(&mut self, s: usize, t: usize) -> Result<()> {
        let (p, e) = (self.segs[s].p, self.segs[s].e);
        let (center, radius) = match self.segs[t].shape {
            Shape::Arc { center, radius, .. } | Shape::Circle { center, radius } => {
                (center, radius)
            }
            Shape::Line => unreachable!("a curve"),
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
        let points: Vec<[I; 2]> = match delta.cmp(&zero()) {
            Ordering::Greater => return Ok(()),
            Ordering::Equal => {
                let foot = [q(center.x) - &fc * a / &ab2, q(center.y) - &fc * b / &ab2];
                vec![[I::exact(foot[0].clone()), I::exact(foot[1].clone())]]
            }
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
            Shape::Line => unreachable!("a curve"),
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

/// The stored segments of a profile.
fn segments(profile: &Profile, op: Operand) -> Vec<Seg> {
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
                        Segment::Spline(_) => unreachable!("refused before"),
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
}

/// A directed result edge, the region on its left.
#[derive(Debug, Clone)]
struct Edge {
    from: Option<PId>,
    to: Option<PId>,
    p: Point2,
    e: Point2,
    shape: Shape,
    origin: Origin2,
}

impl Piece {
    /// An off-centre point of the piece (fraction `0.4453125` along it).
    fn sample(&self, seg: &Seg) -> Point2 {
        let f = 0.4453125;
        match seg.shape {
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
                    return Point2::new(center.x + radius * a.cos(), center.y + radius * a.sin());
                }
                let a0 = (self.p.y - center.y).atan2(self.p.x - center.x);
                let sweep = super::arc_sweep(center, self.p, self.e, true);
                let a = a0 + sweep * f;
                Point2::new(center.x + radius * a.cos(), center.y + radius * a.sin())
            }
        }
    }
}

/// The Boolean of two profiles in one frame: its result profiles with their
/// provenance, in a deterministic order (by their first segment's origin).
pub(crate) fn boolean(a: &Profile, b: &Profile, op: Op2) -> Result<Boolean2> {
    let tolerance = a.tolerance();
    let tol = tolerance.linear();
    for p in [a, b] {
        if p.boundaries().any(|bd| match &bd.kind {
            BoundaryKind::Path { segments, .. } => {
                segments.iter().any(|s| matches!(s, Segment::Spline(_)))
            }
            _ => false,
        }) {
            return Err(Error::OutOfDomain("a Boolean of spline profiles (S9a.2)"));
        }
    }
    let mut segs = segments(a, Operand::A);
    segs.extend(segments(b, Operand::B));
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
            let lines = |x: &Seg| matches!(x.shape, Shape::Line);
            let found = match (lines(&arr.segs[s]), lines(&arr.segs[t])) {
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
            Shape::Circle { .. } if ev.is_empty() => {
                pieces.push(Piece {
                    seg: s,
                    part: 0,
                    whole: true,
                    from: None,
                    to: None,
                    p: seg.p,
                    e: seg.p,
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
            let arcs = |x: &Piece| match arr.segs[x.seg].shape {
                Shape::Arc { ccw, .. } => Some(ccw),
                Shape::Circle { .. } => Some(true),
                Shape::Line => None,
            };
            let same = pa.from == pb.from && pa.to == pb.to;
            let opposite = pa.from == pb.to && pa.to == pb.from;
            let shared = match (arcs(pa), arcs(pb)) {
                (None, None) => same || opposite,
                (Some(ca), Some(cb)) => {
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
        let at = other.classify(piece.sample(&arr.segs[piece.seg]))?;
        if at == Location::Boundary {
            return Err(Error::Degenerate(
                "a boundary within the resolution of the other profile's",
            ));
        }
        place.push(Some(at));
    }
    // The kept pieces, directed with the region on their left: a hole's
    // stored direction reversed.
    let mut edges: Vec<Edge> = Vec::new();
    let seg_ref = |piece: &Piece, arr: &Arrangement| -> SegRef {
        let seg = &arr.segs[piece.seg];
        (seg.op, seg.b, seg.j, (!piece.whole).then_some(piece.part))
    };
    let directed = |piece: &Piece, arr: &Arrangement, flip: bool, origin: Origin2| -> Edge {
        let seg = &arr.segs[piece.seg];
        let reverse = (seg.b > 0) != flip;
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
        let shape = match shape {
            Shape::Circle { center, radius } if reverse => Shape::Circle { center, radius },
            s => s,
        };
        Edge {
            from,
            to,
            p,
            e,
            shape,
            origin,
        }
    };
    // Whole circles keep a direction flag in their origin's order: a
    // reversed whole circle is a clockwise cycle.
    let mut circle_cw: BTreeSet<usize> = BTreeSet::new();
    for (i, piece) in pieces.iter().enumerate() {
        let seg = &arr.segs[piece.seg];
        let (inside, outside) = (Some(Location::Inside), Some(Location::Outside));
        let keep_flip: Option<bool> = match partner.get(&i) {
            Some(&(k, same)) => {
                // In region-left directions: a hole's reversed.
                let other = &arr.segs[pieces[k].seg];
                let agree = same == ((seg.b > 0) == (other.b > 0));
                match (op, seg.op, agree) {
                    (Op2::Fuse | Op2::Common, Operand::A, true) => Some(false),
                    (Op2::Cut, Operand::A, false) => Some(false),
                    _ => None,
                }
            }
            None => match (op, seg.op, place[i]) {
                (Op2::Fuse, _, p) if p == outside => Some(false),
                (Op2::Common, _, p) if p == inside => Some(false),
                (Op2::Cut, Operand::A, p) if p == outside => Some(false),
                (Op2::Cut, Operand::B, p) if p == inside => Some(true),
                _ => None,
            },
        };
        let Some(flip) = keep_flip else { continue };
        let origin = Origin2 {
            from: seg_ref(piece, &arr),
            shared: partner.get(&i).map(|&(k, _)| seg_ref(&pieces[k], &arr)),
        };
        let edge = directed(piece, &arr, flip, origin);
        if piece.from.is_none() && ((seg.b > 0) != flip) {
            circle_cw.insert(edges.len());
        }
        edges.push(edge);
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
        let area = signed_area(&cycle.iter().map(|&k| &edges[k]).collect::<Vec<_>>());
        cycles.push((cycle, area > 0.0));
    }
    // Outer boundaries and their holes.
    let mut outers: Vec<Traced> = Vec::new();
    let mut holes: Vec<Traced> = Vec::new();
    for (cycle, ccw) in &cycles {
        let built = boundary_of(
            &cycle.iter().map(|&k| &edges[k]).collect::<Vec<_>>(),
            !ccw,
            tolerance,
        )?;
        if *ccw {
            outers.push(built);
        } else {
            holes.push(built);
        }
    }
    let mut owned: Vec<Vec<Traced>> = vec![Vec::new(); outers.len()];
    for hole in holes {
        let sample = hole.0.sample();
        let owner = outers
            .iter()
            .position(|(o, ..)| o.locate(sample, tolerance) == Location::Inside)
            .ok_or(Error::Degenerate("a result's hole outside its boundary"))?;
        owned[owner].push(hole);
    }
    let mut out = Vec::new();
    for ((outer, segs, verts), hs) in outers.into_iter().zip(owned) {
        let mut segments = vec![segs];
        let mut vertices = vec![verts];
        let mut hole_boundaries = Vec::new();
        for (h, s, v) in hs {
            hole_boundaries.push(h);
            segments.push(s);
            vertices.push(v);
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
    out.sort_by_key(|p| p.segments[0].iter().min().copied());
    Ok(Boolean2 {
        pieces: out,
        crosses: arr.cross_segs,
    })
}

/// A traced boundary with its segments' and vertices' provenance.
type Traced = (Boundary, Vec<Origin2>, Vec<PId>);

/// A cycle's signed area (binary64, for its orientation).
fn signed_area(cycle: &[&Edge]) -> f64 {
    let mut twice = 0.0;
    for e in cycle {
        match e.shape {
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
fn boundary_of(cycle: &[&Edge], hole: bool, tolerance: Tolerance) -> Result<Traced> {
    if cycle.len() == 1 && cycle[0].from.is_none() {
        let Shape::Circle { center, radius } = cycle[0].shape else {
            return Err(Error::InvalidTopology("a ring that is not a circle"));
        };
        return Ok((
            Boundary::circle(center, radius, tolerance)?,
            vec![cycle[0].origin],
            vec![],
        ));
    }
    let mut items: Vec<(PId, Point2, Segment, Origin2)> = cycle
        .iter()
        .map(|e| {
            let seg = match e.shape {
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
            (e.from.expect("a vertex"), e.p, seg, e.origin)
        })
        .collect();
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
        items.iter().map(|i| i.3).collect(),
        items.iter().map(|i| i.0).collect(),
    ))
}
