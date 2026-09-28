//! S8a.2: a prism split by a plane oblique to its axis.
//!
//! In the prism's frame the plane is `F = a u + b v + c w + d = 0` with
//! `c != 0` and `(a, b) != 0`: over a profile point the material below the
//! plane (in `w`) spans `[low, min(high, g)]`, `g = -(a u + b v + d) / c`,
//! and the material above it `[max(low, g), high]`. A lower piece's
//! footprint is therefore a piece of the profile's exact section by the
//! plane's trace on the bottom cap's plane (`F(u, v, low) = 0`, where its
//! height vanishes), and its top is that footprint's section by the trace on
//! the top cap's plane: the top cap where `g >= high`, the cut face where
//! `g <= high`; an upper piece is the same with the ends exchanged. Walls
//! over the footprint's boundary run from its flat end to its creased one:
//! a planar wall's top a line on the plane, a cylindrical wall's an ellipse
//! arc whose pcurve is the graph `v = a0 + a1 cos u + a2 sin u`.
//!
//! Every side, crossing and tangency is decided exactly on the stored data
//! by the two sections; new vertices, heights, ellipses and pcurves are
//! rounded from exact values, and each piece validates as a general body.
//! Entities are named by provenance: an input entity whole in one piece
//! keeps its id (`Unchanged`, or `Modified` when its stored geometry or
//! bounding ids changed), one in several pieces or in parts is `Split`
//! (vertices and edges lying in the plane into copies), and cut vertices,
//! edges and faces are `Generated` from the input entities they cut.
use super::{piece_ordinal, q, rational_f64, zero, PointId, Section, SectionPiece, SegOrigin};
use crate::certified::Interval as I;
use crate::identity::{
    Derivation, EntityId, EntityKind, OperationId, OperationKind, Parent, ProfileElement, Role,
};
use crate::profile::{arc_sweep, BoundaryKind, Segment};
use crate::topology::{
    plane_pcurve, Curve2, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId,
    Orientation, Region, RegionId, RegionKind, Shell, ShellId, Side as FaceSide, Slot, Surface,
    Topology, TopologyParts, Vertex, VertexId,
};
use crate::{Error, Frame3, Point2, Point3, Profile, Result, Side, Solid, Tolerance};
use num_rational::BigRational as R;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

/// An input prism entity by what it derives from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Key {
    /// The cap at the low (`true`) or high end.
    Cap(bool),
    Region,
    /// The wall of stored segment `j` of boundary `b`.
    Wall(usize, usize),
    /// A cap edge at the low or high end.
    CapEdge(bool, usize, usize),
    /// The vertical edge through stored vertex `j` of boundary `b`.
    Vertical(usize, usize),
    /// A cap vertex at the low or high end.
    CapVertex(bool, usize, usize),
}

/// A slot's provenance before naming.
#[derive(Debug, Clone)]
pub(super) enum Occ {
    /// The input entity, whole.
    Whole(Key),
    /// A part of it, keyed within the piece.
    Part(Key, usize),
    /// New, generated from these input entities.
    New(Vec<Key>, Role),
}

/// One piece before naming: its side, whether it lies under the plane along
/// the prism's axis, its footprint, its parts and every slot's provenance.
pub(super) struct Built {
    pub(super) side: Side,
    pub(super) lower: bool,
    pub(super) footprint: Profile,
    pub(super) parts: TopologyParts,
    pub(super) plans: Vec<(Slot, Occ)>,
}

/// A 3D vertex of a piece by its place in the two sections.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum VKey {
    /// A footprint point at the flat end.
    Flat(PointId),
    /// A footprint point off the flat end's trace, at the creased end.
    Crease(PointId),
    /// The creased end's trace crossing footprint segment `(qb, qj)`.
    Cross(usize, usize, usize),
    /// Where the flat end's trace touches footprint circle `(qb, qj)`.
    Touch(usize, usize),
}

/// A 3D edge of a piece.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum EKey {
    Flat(usize, usize),
    Vertical(PointId),
    /// Part `k` of a footprint segment at the creased end.
    Crease(usize, usize, usize),
    /// A chord of the creased end's trace.
    Chord(usize),
}

/// The plane in the prism's frame and the prism's data.
struct Setup<'a> {
    profile: &'a Profile,
    frame: Frame3,
    low: f64,
    high: f64,
    plane: [R; 4],
    tolerance: Tolerance,
    /// The cut plane's frame: its normal up the prism's axis, its x axis
    /// along the plane's steepest ascent (every section ellipse's major axis).
    cut: Frame3,
    /// `|m| / |c|`: an ellipse's major radius over its circle's radius.
    stretch: f64,
    /// The angle of `(a, b)`: an ellipse's angle is its circle's less this.
    alpha: f64,
}

impl Setup<'_> {
    fn f(&self, p: Point2, w: &R) -> R {
        let [a, b, c, d] = &self.plane;
        a * q(p.x) + b * q(p.y) + c * w + d
    }
    /// The plane's height over a profile point, rounded.
    fn g(&self, p: Point2) -> f64 {
        rational_f64(&self.g_exact(p))
    }
    fn g_exact(&self, p: Point2) -> R {
        let [a, b, c, d] = &self.plane;
        -(a * q(p.x) + b * q(p.y) + d) / c
    }
    fn line(&self, w: f64) -> [R; 3] {
        let [a, b, c, d] = &self.plane;
        [a.clone(), b.clone(), d + c * q(w)]
    }
    fn at(&self, p: Point2, w: f64) -> Point3 {
        self.frame.point(p, w)
    }
    fn level_frame(&self, center: Point2, w: f64) -> Result<Frame3> {
        Frame3::new(
            self.frame.point(center, w),
            self.frame.normal(),
            self.frame.x(),
            self.tolerance,
        )
    }
    fn ellipse_frame(&self, center: Point2) -> Result<Frame3> {
        Frame3::new(
            self.frame.point(center, self.g(center)),
            self.cut.normal(),
            self.cut.x(),
            self.tolerance,
        )
    }
    /// A section curve's sinusoid on the cylinder of `center` and `radius`
    /// (its frame at the low end): `v = g(c) - low - (r / c)(a cos u + b sin u)`.
    fn sinusoid(&self, center: Point2, radius: f64) -> [f64; 3] {
        let [a, b, c, d] = &self.plane;
        let g = -(a * q(center.x) + b * q(center.y) + d) / c;
        [
            rational_f64(&(g - q(self.low))),
            rational_f64(&(-(a * q(radius)) / c)),
            rational_f64(&(-(b * q(radius)) / c)),
        ]
    }
}

/// A stored segment of the input profile.
#[derive(Debug, Clone, Copy)]
enum PSeg {
    Line(Point2, Point2),
    Arc {
        center: Point2,
        radius: f64,
        start: Point2,
        /// The stored start angle and sweep (the prism's).
        a0: f64,
        sweep: f64,
    },
    Circle {
        center: Point2,
        radius: f64,
    },
    /// A spline (S8b.3), stored reversed or not: its geometry is the
    /// footprint's piece of it.
    Spline {
        reversed: bool,
    },
}

fn pseg(profile: &Profile, b: usize, j: usize) -> PSeg {
    let boundary = profile.boundaries().nth(b).expect("a boundary");
    match &boundary.kind {
        BoundaryKind::Circle { center, radius } => PSeg::Circle {
            center: *center,
            radius: *radius,
        },
        BoundaryKind::Polygon(points) => PSeg::Line(points[j], points[(j + 1) % points.len()]),
        BoundaryKind::Path { points, segments } => {
            let (s, e) = (points[j], points[(j + 1) % points.len()]);
            match segments[j].clone() {
                Segment::Spline(span) => PSeg::Spline {
                    reversed: span.is_reversed(),
                },
                Segment::Line => PSeg::Line(s, e),
                Segment::Arc {
                    center,
                    radius,
                    ccw,
                } => PSeg::Arc {
                    center,
                    radius,
                    start: s,
                    a0: (s.y - center.y).atan2(s.x - center.x),
                    sweep: arc_sweep(center, s, e, ccw),
                },
            }
        }
    }
}

/// The turn from `a` to `b` about `center` in the direction of `sign`, in
/// `(0, 2 pi]`.
fn turn(center: Point2, a: Point2, b: Point2, sign: f64) -> f64 {
    let (u, v) = (
        Point2::new(a.x - center.x, a.y - center.y),
        Point2::new(b.x - center.x, b.y - center.y),
    );
    let cross = u.x * v.y - u.y * v.x;
    let dot = u.x * v.x + u.y * v.y;
    let mut t = (sign * cross).atan2(dot);
    if t <= 0.0 {
        t += TAU;
    }
    t
}

fn sign_of(x: &R) -> i8 {
    x.cmp(&zero()) as i8
}

fn side_sign(side: Side) -> i8 {
    match side {
        Side::Below => -1,
        Side::Above => 1,
    }
}

/// Whether the line touches an arc of the profile inside it: the plane then
/// touches a cap's edge at a point, pinching a wall between its ends.
fn tangent_inside(profile: &Profile, line: &[R; 3]) -> Result<bool> {
    let [a, b, d] = line;
    let ab2 = a * a + b * b;
    let touches = |center: Point2, radius: f64| {
        let fc = a * q(center.x) + b * q(center.y) + d;
        &fc * &fc == q(radius) * q(radius) * &ab2
    };
    for boundary in profile.boundaries() {
        match &boundary.kind {
            // A circle touching the trace keeps a zero-height vertex there.
            BoundaryKind::Circle { .. } | BoundaryKind::Polygon(_) => {}
            BoundaryKind::Path { points, segments } => {
                let n = points.len();
                for (j, s) in segments.iter().enumerate() {
                    if let Segment::Spline(span) = s {
                        if !super::spline::meets(span, [a, b, d])?.touches.is_empty() {
                            return Ok(true);
                        }
                        continue;
                    }
                    let Segment::Arc {
                        center,
                        radius,
                        ccw,
                    } = *s
                    else {
                        continue;
                    };
                    if !touches(center, radius) {
                        continue;
                    }
                    let fc = a * q(center.x) + b * q(center.y) + d;
                    let foot = [q(center.x) - &fc * a / &ab2, q(center.y) - &fc * b / &ab2];
                    let t = [I::exact(foot[0].clone()), I::exact(foot[1].clone())];
                    let (p, e) = (points[j], points[(j + 1) % n]);
                    if super::same_point(&t, p)? || super::same_point(&t, e)? {
                        continue;
                    }
                    if super::within_arc(center, p, e, ccw, &t)?.is_some() {
                        return Ok(true);
                    }
                }
            }
        }
    }
    Ok(false)
}

/// The split's pieces in their order (those below the plane first, each
/// side by its footprint's first input segment), or only the one at `only`.
pub(super) fn pieces(
    profile: &Profile,
    frame: Frame3,
    start: f64,
    end: f64,
    plane: &[R; 4],
    only: Option<usize>,
) -> Result<Vec<Built>> {
    let tolerance = profile.tolerance();
    let (low, high) = (start.min(end), start.max(end));
    let [a, b, c, _] = plane;
    let (af, bf, cf) = (rational_f64(a), rational_f64(b), rational_f64(c));
    let up = if *c > zero() { 1.0 } else { -1.0 };
    let (x, y, n) = (frame.x(), frame.y(), frame.normal());
    let m = x * af + y * bf + n * cf;
    let hint = x * af + y * bf + n * rational_f64(&(-(a * a + b * b) / c));
    // The traces on the caps' planes are `|c| (high - low) / |(a, b)|`
    // apart: within the resolution, the plane is parallel to the axis over
    // the prism's height to binary64 (a frame whose stored axes are not
    // exactly orthogonal), and its sections' pieces would be thinner than it.
    let span = q(high) - q(low);
    let ab2 = a * a + b * b;
    if c * c * &span * &span <= q(tolerance.linear()) * q(tolerance.linear()) * &ab2 {
        return Err(Error::Degenerate(
            "a plane within the resolution of parallel to the prism's axis",
        ));
    }
    // The cut face's frame at the foot of the perpendicular from a profile
    // point at mid-height: on the plane, near the solid.
    let anchor = match &profile.outer().kind {
        BoundaryKind::Circle { center, .. } => *center,
        BoundaryKind::Polygon(p) => p[0],
        BoundaryKind::Path { points, .. } => points[0],
    };
    let mid = q(low) / R::from_integer(2.into()) + q(high) / R::from_integer(2.into());
    let [_, _, _, d] = plane;
    let f0 = a * q(anchor.x) + b * q(anchor.y) + c * &mid + d;
    let k = &f0 / (&ab2 + c * c);
    let origin = frame.point(
        Point2::new(
            rational_f64(&(q(anchor.x) - &k * a)),
            rational_f64(&(q(anchor.y) - &k * b)),
        ),
        rational_f64(&(mid - &k * c)),
    );
    let cut = Frame3::new(origin, m * up, hint, tolerance)?;
    let stretch = rational_f64(&((a * a + b * b + c * c) / (c * c))).sqrt();
    let setup = Setup {
        profile,
        frame,
        low,
        high,
        plane: plane.clone(),
        tolerance,
        cut,
        stretch,
        alpha: bf.atan2(af),
    };
    for w in [low, high] {
        if tangent_inside(profile, &setup.line(w))? {
            return Err(Error::Degenerate(
                "a plane tangent to a cap's arc or spline edge",
            ));
        }
    }
    // The lower side's sign of F: material under the plane along the axis.
    let lower_sign: i8 = if *c > zero() { -1 } else { 1 };
    // Each end's section, and the footprints in order.
    let mut sections = Vec::new();
    let mut keys: Vec<(Side, SegOrigin, usize, usize)> = Vec::new();
    for lower in [true, false] {
        let wf = if lower { low } else { high };
        let own = if lower { lower_sign } else { -lower_sign };
        let s1 = Section::new(profile, setup.line(wf))?;
        let side = if own < 0 { Side::Below } else { Side::Above };
        for (i, footprint) in s1.pieces.iter().enumerate() {
            if side_sign(footprint.side) == own {
                keys.push((side, footprint.first_origin(), sections.len(), i));
            }
        }
        sections.push((lower, own, s1));
    }
    keys.sort_by_key(|k| (k.0, k.1));
    let mut out = Vec::new();
    for (k, &(side, _, s, i)) in keys.iter().enumerate() {
        if only.is_some_and(|o| o != k) {
            continue;
        }
        let (lower, own, s1) = &sections[s];
        let footprint = &s1.pieces[i];
        let wk = if *lower { high } else { low };
        let s2 = Section::new(&footprint.profile, setup.line(wk))?;
        let (parts, plans) = build_piece(&setup, s1, footprint, &s2, *lower, *own)?;
        out.push(Built {
            side,
            lower: *lower,
            footprint: footprint.profile.clone(),
            parts,
            plans,
        });
    }
    Ok(out)
}

/// A part of a footprint segment at the creased end: its section part,
/// stored start and end, and relation (1 cap side, 0 on the trace, -1 cut).
type CreasePart = (usize, PointId, PointId, i8);

/// A footprint segment's provenance and geometry.
#[derive(Clone, Copy)]
struct QSeg {
    origin: SegOrigin,
    /// The input segment `(b, j)` it lies on, if not a chord.
    input: Option<(usize, usize)>,
    /// Its stored direction agrees with the input segment's.
    same: bool,
    /// Its stored ends (none for a whole circle).
    ends: Option<(PointId, PointId)>,
    /// No height: a chord of the flat trace, or a line along it.
    flat: bool,
}

#[allow(clippy::too_many_lines)]
fn build_piece(
    setup: &Setup,
    s1: &Section,
    q_piece: &SectionPiece,
    s2: &Section,
    lower: bool,
    own: i8,
) -> Result<(TopologyParts, Vec<(Slot, Occ)>)> {
    let tol = setup.tolerance.linear();
    let (wf, wk) = if lower {
        (setup.low, setup.high)
    } else {
        (setup.high, setup.low)
    };
    let (wfr, wkr) = (q(wf), q(wk));
    let pos1 = |p: PointId| s1.positions[&p];
    let on_flat = |p: PointId| match p {
        PointId::Cross(..) => true,
        PointId::Vertex(..) => setup.f(pos1(p), &wfr) == zero(),
    };
    // 1 on the cap side of the creased end, 0 on its trace, -1 cut.
    let crease_rel = |p: PointId| sign_of(&setup.f(pos1(p), &wkr)) * own;
    let q_profile = &q_piece.profile;
    // A footprint segment's spline, in the footprint's stored direction.
    let q_span = |qb: usize, qj: usize| -> Option<super::spline::Span> {
        match &q_profile.boundaries().nth(qb)?.kind {
            BoundaryKind::Path { segments, .. } => match &segments[qj] {
                Segment::Spline(span) => Some(span.clone()),
                _ => None,
            },
            _ => None,
        }
    };
    // The footprint's segments.
    let mut qsegs: BTreeMap<(usize, usize), QSeg> = BTreeMap::new();
    for (qb, boundary) in q_profile.boundaries().enumerate() {
        let count = match &boundary.kind {
            BoundaryKind::Circle { .. } => 1,
            BoundaryKind::Polygon(p) => p.len(),
            BoundaryKind::Path { points, .. } => points.len(),
        };
        for qj in 0..count {
            let origin = q_piece.segments[qb][qj];
            let ends = (!q_piece.vertices[qb].is_empty()).then(|| {
                (
                    q_piece.vertices[qb][qj],
                    q_piece.vertices[qb][(qj + 1) % count],
                )
            });
            let input = match origin {
                SegOrigin::Whole(b, j) | SegOrigin::Part(b, j, _) => Some((b, j)),
                SegOrigin::Chord(_) => None,
            };
            let flat = match (input, ends) {
                (None, _) => true,
                (Some((b, j)), Some((s, e))) => {
                    matches!(pseg(setup.profile, b, j), PSeg::Line(..)) && on_flat(s) && on_flat(e)
                }
                _ => false,
            };
            qsegs.insert(
                (qb, qj),
                QSeg {
                    origin,
                    input,
                    same: input.is_some_and(|(b, _)| (b > 0) == (qb > 0)),
                    ends,
                    flat,
                },
            );
        }
    }
    // The creased end's parts of each footprint segment, in stored order:
    // (part, stored start, stored end, 1 cap / 0 on the trace / -1 cut).
    let mut crease_parts: BTreeMap<(usize, usize), Vec<CreasePart>> = BTreeMap::new();
    for e in &s2.edges {
        let (qb, qj, k) = match e.origin {
            SegOrigin::Whole(qb, qj) => (qb, qj, 0),
            SegOrigin::Part(qb, qj, k) => (qb, qj, k),
            SegOrigin::Chord(_) => continue,
        };
        let (s, t) = if qb > 0 {
            (e.to, e.from)
        } else {
            (e.from, e.to)
        };
        let rel = if e.side == 0 { 0 } else { e.side * own };
        crease_parts
            .entry((qb, qj))
            .or_default()
            .push((k, s, t, rel));
    }
    for parts in crease_parts.values_mut() {
        parts.sort_by_key(|p| p.0);
    }
    // A section point of the creased end as a 3D vertex.
    let s2_key = |p: PointId| -> VKey {
        match p {
            PointId::Vertex(qb, qj) => {
                let pp = q_piece.vertices[qb][qj];
                if on_flat(pp) {
                    VKey::Flat(pp)
                } else {
                    VKey::Crease(pp)
                }
            }
            PointId::Cross(qb, qj, k) => VKey::Cross(qb, qj, k),
        }
    };
    // Footprint circles the flat end's trace touches: the point and its
    // angle, where the wall's height vanishes.
    let mut touch: BTreeMap<(usize, usize), (Point2, f64)> = BTreeMap::new();
    for (&(qb, qj), qs) in &qsegs {
        let (Some((b, j)), None) = (qs.input, qs.ends) else {
            continue;
        };
        let PSeg::Circle { center, radius } = pseg(setup.profile, b, j) else {
            continue;
        };
        let [a, bb, d] = setup.line(wf);
        let ab2 = &a * &a + &bb * &bb;
        let fc = &a * q(center.x) + &bb * q(center.y) + &d;
        if &fc * &fc != q(radius) * q(radius) * &ab2 {
            continue;
        }
        let foot = Point2::new(
            rational_f64(&(q(center.x) - &fc * &a / &ab2)),
            rational_f64(&(q(center.y) - &fc * &bb / &ab2)),
        );
        let angle = (foot.y - center.y).atan2(foot.x - center.x);
        touch.insert(
            (qb, qj),
            (foot, if angle < 0.0 { angle + TAU } else { angle }),
        );
    }
    let mut parts = TopologyParts::default();
    let mut plans: Vec<(Slot, Occ)> = Vec::new();
    let mut vids: BTreeMap<VKey, VertexId> = BTreeMap::new();
    let mut heights: BTreeMap<VKey, f64> = BTreeMap::new();
    // Vertices: every footprint point at the flat end, off its trace at the
    // creased end too, and the creased trace's crossings.
    let mut add_vertex = |key: VKey,
                          p2: Point2,
                          w: f64,
                          occ: Occ,
                          parts: &mut TopologyParts,
                          plans: &mut Vec<(Slot, Occ)>| {
        if vids.contains_key(&key) {
            return;
        }
        let id = VertexId(parts.vertices.len());
        parts.vertices.push(Vertex {
            position: setup.at(p2, w),
            enclosure: None,
        });
        plans.push((Slot::Vertex(id), occ));
        vids.insert(key, id);
        heights.insert(key, w);
    };
    let at_low = |w: f64| w == setup.low;
    for qs in qsegs.values() {
        let Some((p, _)) = qs.ends else { continue };
        let p2 = pos1(p);
        let occ = match p {
            PointId::Vertex(b, j) => Occ::Whole(Key::CapVertex(at_low(wf), b, j)),
            PointId::Cross(b, j, _) => {
                Occ::New(vec![Key::CapEdge(at_low(wf), b, j)], Role::CutVertex)
            }
        };
        add_vertex(VKey::Flat(p), p2, wf, occ, &mut parts, &mut plans);
        if on_flat(p) {
            continue;
        }
        let PointId::Vertex(b, j) = p else {
            unreachable!("a crossing lies on its trace")
        };
        let (w, occ) = match crease_rel(p) {
            -1 => {
                let w = setup.g(p2);
                if (w - wf).abs() <= tol || (w - wk).abs() <= tol {
                    return Err(Error::Degenerate(
                        "a vertex within the resolution of the plane",
                    ));
                }
                (w, Occ::New(vec![Key::Vertical(b, j)], Role::CutVertex))
            }
            _ => (wk, Occ::Whole(Key::CapVertex(at_low(wk), b, j))),
        };
        add_vertex(VKey::Crease(p), p2, w, occ, &mut parts, &mut plans);
    }
    for (&(qb, qj), &(foot, _)) in &touch {
        let (b, j) = qsegs[&(qb, qj)].input.expect("a circle");
        add_vertex(
            VKey::Touch(qb, qj),
            foot,
            wf,
            Occ::New(vec![Key::CapEdge(at_low(wf), b, j)], Role::CutVertex),
            &mut parts,
            &mut plans,
        );
    }
    for (&(qb, qj), list) in &crease_parts {
        for &(_, s, t, _) in list {
            for p in [s, t] {
                if let PointId::Cross(..) = p {
                    let (b, j) = qsegs[&(qb, qj)].input.expect("a crossing on a segment");
                    add_vertex(
                        s2_key(p),
                        s2.positions[&p],
                        wk,
                        Occ::New(vec![Key::CapEdge(at_low(wk), b, j)], Role::CutVertex),
                        &mut parts,
                        &mut plans,
                    );
                }
            }
        }
    }
    let vid = |k: VKey| vids[&k];
    let v3 = |parts: &TopologyParts, k: VKey| parts.vertices[vids[&k].0].position;
    // Edges.
    let mut eids: BTreeMap<EKey, EdgeId> = BTreeMap::new();
    let add_edge = |key: EKey,
                    ends: Option<(VertexId, VertexId)>,
                    curve: Curve3,
                    occ: Occ,
                    parts: &mut TopologyParts,
                    plans: &mut Vec<(Slot, Occ)>,
                    eids: &mut BTreeMap<EKey, EdgeId>| {
        let id = EdgeId(parts.edges.len());
        parts.edges.push(Edge {
            start: ends.map(|e| e.0),
            end: ends.map(|e| e.1),
            curve,
            fins: Vec::new(),
        });
        plans.push((Slot::Edge(id), occ));
        eids.insert(key, id);
    };
    // The angles of a footprint segment's points along its input arc, in
    // the input's direction: the input's own at its stored ends.
    let angles = |b: usize, j: usize, points: &[(PointId, Point2)]| -> Vec<f64> {
        match pseg(setup.profile, b, j) {
            PSeg::Line(..) => vec![0.0; points.len()],
            PSeg::Spline { .. } => unreachable!("a spline's points are by its parameter"),
            PSeg::Arc {
                center,
                start,
                a0,
                sweep,
                ..
            } => {
                let sign = sweep.signum();
                let n_p = setup
                    .profile
                    .boundaries()
                    .nth(b)
                    .map(|bd| match &bd.kind {
                        BoundaryKind::Path { points, .. } => points.len(),
                        _ => 0,
                    })
                    .unwrap_or(0);
                let mut out: Vec<f64> = Vec::new();
                for (k, (id, p)) in points.iter().enumerate() {
                    let u = if *id == PointId::Vertex(b, j) {
                        a0
                    } else if *id == PointId::Vertex(b, (j + 1) % n_p) && k > 0 {
                        a0 + sweep
                    } else if k == 0 {
                        a0 + sign * turn(center, start, *p, sign)
                    } else {
                        out[k - 1] + sign * turn(center, points[k - 1].1, *p, sign)
                    };
                    out.push(u);
                }
                out
            }
            PSeg::Circle { center, .. } => {
                let mut out: Vec<f64> = Vec::new();
                for (k, (_, p)) in points.iter().enumerate() {
                    let u = if k == 0 {
                        let a = (p.y - center.y).atan2(p.x - center.x);
                        if a < 0.0 {
                            a + TAU
                        } else {
                            a
                        }
                    } else {
                        out[k - 1] + turn(center, points[k - 1].1, *p, 1.0)
                    };
                    out.push(u);
                }
                out
            }
        }
    };
    // Per footprint segment: its points along the input's direction with
    // their angles, and each crease part's range in that list.
    struct Chain {
        points: Vec<(VKey, VKey, Point2)>,
        u: Vec<f64>,
        /// (edge, from index, to index, relation), in the input's direction.
        crease: Vec<(usize, usize, usize, i8)>,
        /// Each creased-end section part's edges, in the input's direction.
        parts: BTreeMap<usize, Vec<usize>>,
    }
    let mut chains: BTreeMap<(usize, usize), Chain> = BTreeMap::new();
    for (&(qb, qj), qs) in &qsegs {
        let Some((b, j)) = qs.input else { continue };
        if qs.flat {
            continue;
        }
        let list = crease_parts.get(&(qb, qj)).cloned().unwrap_or_default();
        if let Some(&(foot, u_t)) = touch.get(&(qb, qj)) {
            // From the touch point round the circle (counter-clockwise) and
            // back: the crossings in between, each interval its section
            // part's relation, the part through the touch point in two.
            let center = match pseg(setup.profile, b, j) {
                PSeg::Circle { center, .. } => center,
                _ => unreachable!("a circle"),
            };
            let offset = |p: Point2| turn(center, foot, p, 1.0) % TAU;
            let mut xs: Vec<(f64, PointId)> = list
                .iter()
                .filter(|p| matches!(p.1, PointId::Cross(..)))
                .map(|p| (offset(s2.positions[&p.1]), p.1))
                .collect();
            xs.sort_by(|x, y| x.0.total_cmp(&y.0));
            let mut points = vec![(VKey::Touch(qb, qj), VKey::Touch(qb, qj), foot)];
            let mut u = vec![u_t];
            for (o, x) in &xs {
                points.push((s2_key(*x), s2_key(*x), s2.positions[x]));
                u.push(u_t + o);
            }
            points.push((VKey::Touch(qb, qj), VKey::Touch(qb, qj), foot));
            u.push(u_t + TAU);
            let n = points.len() - 1;
            let mut crease = Vec::new();
            let mut by_part: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for i in 0..n {
                // The section part starting at this interval's start, or
                // (from the touch point) the one ending at its end.
                let part = if i > 0 {
                    list.iter().find(|p| p.1 == xs[i - 1].1)
                } else if let Some(first) = xs.first() {
                    list.iter().find(|p| p.2 == first.1)
                } else {
                    list.first()
                }
                .copied()
                .ok_or(Error::InvalidTopology(
                    "a section part without its crossing",
                ))?;
                crease.push((i, i, i + 1, part.3));
                by_part.entry(part.0).or_default().push(i);
            }
            // The part through the touch point runs from the last interval
            // to the first.
            for list in by_part.values_mut() {
                if list.len() == 2 {
                    list.reverse();
                }
            }
            chains.insert(
                (qb, qj),
                Chain {
                    points,
                    u,
                    crease,
                    parts: by_part,
                },
            );
            continue;
        }
        // Stored-order points: the footprint's ends and the crossings.
        let mut stored: Vec<(PointId, PointId)> = Vec::new(); // (s1 id or s2 id, s2 id)
        match qs.ends {
            Some((s, e)) => {
                stored.push((s, PointId::Vertex(qb, qj)));
                for &(k, _, t, _) in &list {
                    if k + 1 < list.len() {
                        stored.push((t, t));
                    }
                }
                let count = q_piece.vertices[qb].len();
                stored.push((e, PointId::Vertex(qb, (qj + 1) % count)));
            }
            None => {
                // A whole circle: its crossings, or none.
                for &(_, s, _, _) in &list {
                    if let PointId::Cross(..) = s {
                        stored.push((s, s));
                    }
                }
                if let Some(first) = stored.first().copied() {
                    stored.push(first);
                }
            }
        }
        let mut crease: Vec<(usize, usize, usize, i8)> = list
            .iter()
            .enumerate()
            .map(|(i, &(k, _, _, rel))| (k, i, i + 1, rel))
            .collect();
        let mut keyed: Vec<(VKey, VKey, Point2)> = stored
            .iter()
            .map(|&(p, p2)| match p {
                PointId::Vertex(..) if qs.ends.is_some() => {
                    let flat_key = VKey::Flat(p);
                    let crease_key = if on_flat(p) {
                        VKey::Flat(p)
                    } else {
                        VKey::Crease(p)
                    };
                    (flat_key, crease_key, pos1(p))
                }
                _ => (s2_key(p2), s2_key(p2), s2.positions[&p2]),
            })
            .collect();
        if !qs.same && qs.ends.is_some() {
            keyed.reverse();
            let last = keyed.len() - 1;
            for c in &mut crease {
                let (i, t) = (c.1, c.2);
                *c = (c.0, last - t, last - i, c.3);
            }
            crease.reverse();
        }
        let ids: Vec<(PointId, Point2)> = keyed
            .iter()
            .enumerate()
            .map(|(i, k)| {
                let id = match k.0 {
                    VKey::Flat(p) | VKey::Crease(p) => p,
                    // Not an input vertex: any id other than one.
                    VKey::Cross(..) | VKey::Touch(..) => PointId::Cross(usize::MAX, 0, i),
                };
                (id, k.2)
            })
            .collect();
        let spline =
            q_span(qb, qj).filter(|_| matches!(pseg(setup.profile, b, j), PSeg::Spline { .. }));
        let u = if let Some(span) = spline {
            // A spline's points by its curve's parameter, which its pieces
            // share with the input's curve.
            let (first, last) = span.curve().as_curve3().domain();
            let (start, end) = if span.is_reversed() {
                (last, first)
            } else {
                (first, last)
            };
            let n = stored.len();
            let mut us: Vec<f64> = stored
                .iter()
                .enumerate()
                .map(|(i, &(_, p2))| {
                    if i == 0 {
                        start
                    } else if i + 1 == n {
                        end
                    } else {
                        s2.params[&p2]
                    }
                })
                .collect();
            if !qs.same {
                us.reverse();
            }
            us
        } else if qs.ends.is_none() && ids.is_empty() {
            Vec::new()
        } else {
            angles(b, j, &ids)
        };
        let parts_of = crease.iter().map(|c| (c.0, vec![c.0])).collect();
        chains.insert(
            (qb, qj),
            Chain {
                points: keyed,
                u,
                crease,
                parts: parts_of,
            },
        );
    }
    // Flat edges: every footprint segment at the flat end.
    for (&(qb, qj), qs) in &qsegs {
        let occ = match qs.origin {
            SegOrigin::Whole(b, j) => Occ::Whole(Key::CapEdge(at_low(wf), b, j)),
            SegOrigin::Part(b, j, k) => Occ::Part(Key::CapEdge(at_low(wf), b, j), k << 10),
            SegOrigin::Chord(_) => Occ::New(vec![Key::Cap(at_low(wf))], Role::CutEdge),
        };
        let (ends, curve) = match (qs.input, qs.ends) {
            (None, Some((s, e))) => {
                // A chord of the flat trace, from its first point.
                let chord = match qs.origin {
                    SegOrigin::Chord(c) => s1.chords[c],
                    _ => unreachable!(),
                };
                let (from, to) = if chord.0 == s { (s, e) } else { (e, s) };
                let (a, b) = (VKey::Flat(from), VKey::Flat(to));
                (
                    Some((vid(a), vid(b))),
                    Curve3::LineSegment {
                        start: v3(&parts, a),
                        end: v3(&parts, b),
                    },
                )
            }
            (Some((b, j)), ends) => match (pseg(setup.profile, b, j), ends) {
                (PSeg::Line(..), _) => {
                    let (s, e) = ends.expect("a line has ends");
                    let (s, e) = if qs.same { (s, e) } else { (e, s) };
                    let (a, bb) = (VKey::Flat(s), VKey::Flat(e));
                    (
                        Some((vid(a), vid(bb))),
                        Curve3::LineSegment {
                            start: v3(&parts, a),
                            end: v3(&parts, bb),
                        },
                    )
                }
                (PSeg::Spline { .. }, _) => {
                    let chain = &chains[&(qb, qj)];
                    let last = chain.points.len() - 1;
                    let span = q_span(qb, qj).expect("a spline footprint segment");
                    let span = if qs.same { span } else { span.reversed() };
                    (
                        Some((vid(chain.points[0].0), vid(chain.points[last].0))),
                        crate::topology::lifted_spline(&span, setup.frame, wf)?,
                    )
                }
                // An arc, or a part of a circle the flat trace crosses.
                (PSeg::Arc { center, radius, .. }, _)
                | (PSeg::Circle { center, radius }, Some(_)) => {
                    let chain = &chains[&(qb, qj)];
                    let last = chain.points.len() - 1;
                    let (u0, u1) = (chain.u[0], chain.u[last]);
                    (
                        Some((vid(chain.points[0].0), vid(chain.points[last].0))),
                        Curve3::CircularArc {
                            frame: setup.level_frame(center, wf)?,
                            radius,
                            start_angle: u0,
                            sweep_angle: u1 - u0,
                        },
                    )
                }
                (PSeg::Circle { center, radius }, None) => match touch.get(&(qb, qj)) {
                    None => (
                        None,
                        Curve3::Circle {
                            frame: setup.level_frame(center, wf)?,
                            radius,
                        },
                    ),
                    // From the touch point round to it.
                    Some(&(_, u_t)) => {
                        let v = vid(VKey::Touch(qb, qj));
                        (
                            Some((v, v)),
                            Curve3::CircularArc {
                                frame: setup.level_frame(center, wf)?,
                                radius,
                                start_angle: u_t,
                                sweep_angle: TAU,
                            },
                        )
                    }
                },
            },
            (None, None) => unreachable!("a chord has ends"),
        };
        add_edge(
            EKey::Flat(qb, qj),
            ends,
            curve,
            occ,
            &mut parts,
            &mut plans,
            &mut eids,
        );
    }
    // Vertical edges: at every footprint vertex off the flat trace, from its
    // lower end up.
    let mut verticals: BTreeSet<PointId> = BTreeSet::new();
    for qs in qsegs.values() {
        if let Some((p, _)) = qs.ends {
            if !on_flat(p) {
                verticals.insert(p);
            }
        }
    }
    for p in verticals {
        let PointId::Vertex(b, j) = p else {
            unreachable!("a crossing lies on its trace")
        };
        let (f, k) = (VKey::Flat(p), VKey::Crease(p));
        let (lo, hi) = if lower { (f, k) } else { (k, f) };
        let occ = if crease_rel(p) < 0 {
            Occ::Part(Key::Vertical(b, j), 0)
        } else {
            Occ::Whole(Key::Vertical(b, j))
        };
        add_edge(
            EKey::Vertical(p),
            Some((vid(lo), vid(hi))),
            Curve3::LineSegment {
                start: v3(&parts, lo),
                end: v3(&parts, hi),
            },
            occ,
            &mut parts,
            &mut plans,
            &mut eids,
        );
    }
    // Crease edges: each part of a footprint segment at the creased end, on
    // the cap side (or along the trace) at its level, cut on the plane.
    for (&(qb, qj), chain) in &chains {
        let qs = qsegs[&(qb, qj)];
        let (b, j) = qs.input.expect("a wall");
        let single = chain.crease.len() == 1;
        let ring = qs.ends.is_none() && chain.points.is_empty();
        for &(k, i, t, rel) in &chain.crease {
            let (pa, pb) = if ring {
                (VKey::Cross(qb, qj, 0), VKey::Cross(qb, qj, 0))
            } else {
                (chain.points[i].1, chain.points[t].1)
            };
            let occ = if rel >= 0 {
                match (qs.origin, single) {
                    (SegOrigin::Whole(..), true) => Occ::Whole(Key::CapEdge(at_low(wk), b, j)),
                    (SegOrigin::Part(_, _, k1), _) => {
                        Occ::Part(Key::CapEdge(at_low(wk), b, j), (k1 << 10) + k + 1)
                    }
                    _ => Occ::Part(Key::CapEdge(at_low(wk), b, j), k + 1),
                }
            } else {
                Occ::New(vec![Key::Wall(b, j)], Role::CutEdge)
            };
            let ends = (!ring).then(|| (vid(pa), vid(pb)));
            let curve = match (pseg(setup.profile, b, j), rel >= 0) {
                (PSeg::Line(..), _) => Curve3::LineSegment {
                    start: v3(&parts, pa),
                    end: v3(&parts, pb),
                },
                (PSeg::Spline { .. }, cap) => {
                    let span = q_span(qb, qj).expect("a spline footprint segment");
                    let part = super::spline::piece(&span, chain.u[i], chain.u[t])?;
                    if cap {
                        crate::topology::lifted_spline(&part, setup.frame, wk)?
                    } else {
                        // Its ends at their vertices' heights (a crossing's
                        // is its level's), its other poles on the plane.
                        let ends = [
                            (chain.points[i].2, heights[&pa]),
                            (chain.points[t].2, heights[&pb]),
                        ];
                        let height = |p: Point2| {
                            ends.iter()
                                .find(|e| e.0 == p)
                                .map_or_else(|| setup.g(p), |e| e.1)
                        };
                        super::spline::plane_image(&part, setup.frame, &height)?
                    }
                }
                (PSeg::Arc { center, radius, .. }, true) => Curve3::CircularArc {
                    frame: setup.level_frame(center, wk)?,
                    radius,
                    start_angle: chain.u[i],
                    sweep_angle: chain.u[t] - chain.u[i],
                },
                (PSeg::Circle { center, radius }, true) if ring => Curve3::Circle {
                    frame: setup.level_frame(center, wk)?,
                    radius,
                },
                (PSeg::Circle { center, radius }, true) => Curve3::CircularArc {
                    frame: setup.level_frame(center, wk)?,
                    radius,
                    start_angle: chain.u[i],
                    sweep_angle: chain.u[t] - chain.u[i],
                },
                (PSeg::Arc { center, radius, .. } | PSeg::Circle { center, radius }, false) => {
                    let (u0, sweep) = if ring {
                        (0.0, TAU)
                    } else {
                        (chain.u[i], chain.u[t] - chain.u[i])
                    };
                    Curve3::EllipseArc {
                        frame: setup.ellipse_frame(center)?,
                        major: radius * setup.stretch,
                        minor: radius,
                        start_angle: u0 - setup.alpha,
                        sweep_angle: sweep,
                    }
                }
            };
            add_edge(
                EKey::Crease(qb, qj, k),
                ends,
                curve,
                occ,
                &mut parts,
                &mut plans,
                &mut eids,
            );
        }
    }
    // Chords of the creased trace.
    for (c, &(s, t)) in s2.chords.iter().enumerate() {
        let (a, b) = (s2_key(s), s2_key(t));
        add_edge(
            EKey::Chord(c),
            Some((vid(a), vid(b))),
            Curve3::LineSegment {
                start: v3(&parts, a),
                end: v3(&parts, b),
            },
            Occ::New(vec![Key::Cap(at_low(wk))], Role::CutEdge),
            &mut parts,
            &mut plans,
            &mut eids,
        );
    }
    // A creased-end section part's edges in the input's direction (a
    // zero-height segment's is its flat edge).
    let crease_edges = |qb: usize, qj: usize, k: usize| -> Vec<EdgeId> {
        if qsegs[&(qb, qj)].flat {
            vec![eids[&EKey::Flat(qb, qj)]]
        } else {
            chains[&(qb, qj)].parts[&k]
                .iter()
                .map(|e| eids[&EKey::Crease(qb, qj, *e)])
                .collect()
        }
    };
    let crease_edge = |qb: usize, qj: usize, k: usize| eids[&EKey::Crease(qb, qj, k)];
    // Faces.
    let add_face = |surface: Surface,
                    sense: Orientation,
                    loops: Vec<(Vec<Fin>, [i32; 2])>,
                    occ: Occ,
                    parts: &mut TopologyParts,
                    plans: &mut Vec<(Slot, Occ)>| {
        let mut ids = Vec::new();
        for (fins, winding) in loops {
            let mut list = Vec::new();
            for fin in fins {
                parts.edges[fin.edge.0].fins.push(FinId(parts.fins.len()));
                list.push(FinId(parts.fins.len()));
                parts.fins.push(fin);
            }
            ids.push(LoopId(parts.loops.len()));
            parts.loops.push(Loop::Edges {
                fins: list,
                winding,
            });
        }
        let f = FaceId(parts.faces.len());
        parts.faces.push(Face {
            surface,
            sense,
            loops: ids,
            front: ShellId(0),
            back: ShellId(1),
            enclosure: None,
        });
        plans.push((Slot::Face(f), occ));
    };
    let cap_frame = |w: f64| -> Result<Frame3> {
        let normal = if w == setup.low {
            -setup.frame.normal()
        } else {
            setup.frame.normal()
        };
        Frame3::new(
            setup.frame.point(Point2::default(), w),
            normal,
            setup.frame.x(),
            setup.tolerance,
        )
    };
    // A planar region's loops: outer first, material on the left seen from
    // the face's outward side (`up`: along the prism's axis).
    let region_loops = |region: &Profile,
                        up: bool,
                        frame: Frame3,
                        seg: &dyn Fn(usize, usize) -> Vec<(EdgeId, bool)>,
                        parts: &TopologyParts|
     -> Vec<(Vec<Fin>, [i32; 2])> {
        let mut out = Vec::new();
        for (rb, boundary) in region.boundaries().enumerate() {
            let count = match &boundary.kind {
                BoundaryKind::Circle { .. } => 1,
                BoundaryKind::Polygon(p) => p.len(),
                BoundaryKind::Path { points, .. } => points.len(),
            };
            let stored = (rb == 0) == up;
            let order: Vec<usize> = if stored {
                (0..count).collect()
            } else {
                (0..count).rev().collect()
            };
            let fins = order
                .into_iter()
                .flat_map(|rj| {
                    let mut uses = seg(rb, rj);
                    if !stored {
                        uses.reverse();
                    }
                    uses
                })
                .map(|(edge, agrees)| {
                    let sense = if agrees == stored {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    };
                    Fin {
                        edge,
                        sense,
                        pcurve: plane_pcurve(&parts.edges[edge.0].curve, sense, frame),
                        enclosure: None,
                    }
                })
                .collect();
            out.push((fins, [0, 0]));
        }
        out
    };
    // The flat face: the footprint at the flat end.
    let q_seg_edge = |qb: usize, qj: usize| -> Vec<(EdgeId, bool)> {
        let qs = qsegs[&(qb, qj)];
        let agrees = match qs.origin {
            SegOrigin::Chord(c) => qs.ends.map(|e| e.0) == Some(s1.chords[c].0),
            _ => qs.same,
        };
        vec![(eids[&EKey::Flat(qb, qj)], agrees)]
    };
    let whole = q_piece
        .segments
        .iter()
        .flatten()
        .all(|o| matches!(o, SegOrigin::Whole(..)))
        && q_piece.segments.iter().map(Vec::len).sum::<usize>()
            == setup
                .profile
                .boundaries()
                .map(|b| match &b.kind {
                    BoundaryKind::Circle { .. } => 1,
                    BoundaryKind::Polygon(p) => p.len(),
                    BoundaryKind::Path { points, .. } => points.len(),
                })
                .sum::<usize>();
    let flat_frame = cap_frame(wf)?;
    let loops = region_loops(q_profile, !lower, flat_frame, &q_seg_edge, &parts);
    add_face(
        Surface::Plane(flat_frame),
        Orientation::Forward,
        loops,
        if whole {
            Occ::Whole(Key::Cap(at_low(wf)))
        } else {
            Occ::Part(Key::Cap(at_low(wf)), 0)
        },
        &mut parts,
        &mut plans,
    );
    // The creased end's faces: caps on the cap side, cut faces on the plane.
    let crease_frame = cap_frame(wk)?;
    for (index, sub) in s2.pieces.iter().enumerate() {
        let cap = side_sign(sub.side) == own;
        let seg = |rb: usize, rj: usize| -> Vec<(EdgeId, bool)> {
            match sub.segments[rb][rj] {
                SegOrigin::Whole(qb, qj) | SegOrigin::Part(qb, qj, _) => {
                    let k = match sub.segments[rb][rj] {
                        SegOrigin::Part(_, _, k) => k,
                        _ => 0,
                    };
                    // The region's stored direction against the footprint's,
                    // then the footprint's against the edges'.
                    let with_q = (rb > 0) == (qb > 0);
                    let qs = qsegs[&(qb, qj)];
                    let q_agrees = match qs.origin {
                        SegOrigin::Chord(c) => qs.ends.map(|e| e.0) == Some(s1.chords[c].0),
                        _ => qs.same,
                    };
                    let agrees = with_q == q_agrees;
                    let mut edges = crease_edges(qb, qj, k);
                    if !agrees {
                        edges.reverse();
                    }
                    edges.into_iter().map(|e| (e, agrees)).collect()
                }
                SegOrigin::Chord(c) => {
                    let start = sub.vertices[rb].get(rj).copied();
                    vec![(eids[&EKey::Chord(c)], start == Some(s2.chords[c].0))]
                }
            }
        };
        if cap {
            let loops = region_loops(&sub.profile, lower, crease_frame, &seg, &parts);
            add_face(
                Surface::Plane(crease_frame),
                Orientation::Forward,
                loops,
                Occ::Part(Key::Cap(at_low(wk)), index + 1),
                &mut parts,
                &mut plans,
            );
        } else {
            let loops = region_loops(&sub.profile, lower, setup.cut, &seg, &parts);
            // Generated from the faces its boundary runs along.
            let mut from: BTreeSet<Key> = BTreeSet::new();
            for list in &sub.segments {
                for o in list {
                    match *o {
                        SegOrigin::Whole(qb, qj) | SegOrigin::Part(qb, qj, _) => {
                            let qs = qsegs[&(qb, qj)];
                            match (qs.flat, qs.input) {
                                (false, Some((b, j))) => {
                                    from.insert(Key::Wall(b, j));
                                }
                                _ => {
                                    from.insert(Key::Cap(at_low(wf)));
                                }
                            }
                        }
                        SegOrigin::Chord(_) => {
                            from.insert(Key::Cap(at_low(wk)));
                        }
                    }
                }
            }
            add_face(
                Surface::Plane(setup.cut),
                if lower {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                },
                loops,
                Occ::New(from.into_iter().collect(), Role::CutFace),
                &mut parts,
                &mut plans,
            );
        }
    }
    // Walls: over every footprint segment with a height, from its lower
    // boundary along the material-left direction, up its far end, back
    // along its upper boundary and down its near end.
    let height = setup.high - setup.low;
    let v_of = |w: f64| if w == setup.low { 0.0 } else { height };
    for (&(qb, qj), chain) in &chains {
        let qs = qsegs[&(qb, qj)];
        let (b, j) = qs.input.expect("a wall");
        // Material-left runs along the input's stored direction on its
        // outer boundary, against it on a hole.
        let along = b == 0;
        let last = chain.points.len().saturating_sub(1);
        let flat_edge = eids[&EKey::Flat(qb, qj)];
        let sense = |forward: bool| {
            if forward {
                Orientation::Forward
            } else {
                Orientation::Reversed
            }
        };
        // The crease parts in traversal order (along or against the input).
        let crease_in = |forward: bool| -> Vec<(EdgeId, bool, usize, usize, i8)> {
            let mut list: Vec<(EdgeId, bool, usize, usize, i8)> = chain
                .crease
                .iter()
                .map(|&(k, i, t, rel)| {
                    let e = crease_edge(qb, qj, k);
                    if forward {
                        (e, true, i, t, rel)
                    } else {
                        (e, false, t, i, rel)
                    }
                })
                .collect();
            if !forward {
                list.reverse();
            }
            list
        };
        let occ = match (qs.origin, chain.crease.iter().all(|c| c.3 >= 0)) {
            (SegOrigin::Whole(..), true) => Occ::Whole(Key::Wall(b, j)),
            (SegOrigin::Part(_, _, k), _) => Occ::Part(Key::Wall(b, j), k),
            _ => Occ::Part(Key::Wall(b, j), 0),
        };
        match pseg(setup.profile, b, j) {
            PSeg::Line(p0, p1) => {
                // The input wall's own frame.
                let (sa, sb) = if along { (p0, p1) } else { (p1, p0) };
                let origin = setup.at(sa, setup.low);
                let tangent = setup.at(sb, setup.low) - origin;
                let wall = Frame3::new(
                    origin,
                    tangent.cross(setup.frame.normal()),
                    tangent,
                    setup.tolerance,
                )?;
                let (near, far) = if along {
                    (chain.points[0], chain.points[last])
                } else {
                    (chain.points[last], chain.points[0])
                };
                let mut uses: Vec<(EdgeId, bool)> = Vec::new();
                let vertical_up = |p: (VKey, VKey, Point2), up: bool| -> Option<(EdgeId, bool)> {
                    let VKey::Flat(id) = p.0 else { return None };
                    eids.get(&EKey::Vertical(id)).map(|e| (*e, up))
                };
                if lower {
                    uses.push((flat_edge, along));
                    uses.extend(vertical_up(far, true));
                    uses.extend(crease_in(!along).into_iter().map(|c| (c.0, c.1)));
                    uses.extend(vertical_up(near, false));
                } else {
                    uses.extend(crease_in(along).into_iter().map(|c| (c.0, c.1)));
                    uses.extend(vertical_up(far, true));
                    uses.push((flat_edge, !along));
                    uses.extend(vertical_up(near, false));
                }
                let fins = uses
                    .into_iter()
                    .map(|(edge, forward)| Fin {
                        edge,
                        sense: sense(forward),
                        pcurve: plane_pcurve(&parts.edges[edge.0].curve, sense(forward), wall),
                        enclosure: None,
                    })
                    .collect();
                add_face(
                    Surface::Plane(wall),
                    Orientation::Forward,
                    vec![(fins, [0, 0])],
                    occ,
                    &mut parts,
                    &mut plans,
                );
            }
            kind => {
                let ring = qs.ends.is_none() && !touch.contains_key(&(qb, qj));
                // A spline wall's footprint piece (S8b.3).
                let spline = match kind {
                    PSeg::Spline { .. } => q_span(qb, qj),
                    _ => None,
                };
                // The wall's surface, and whether material-left increases
                // its `u`: counter-clockwise on an outer arc when the input's
                // sweep is positive, along a spline's parameter when the
                // input is stored forward.
                let (surface, ccw, a_coef) = match kind {
                    PSeg::Arc {
                        center,
                        radius,
                        sweep,
                        ..
                    } => (
                        Surface::Cylinder {
                            frame: setup.level_frame(center, setup.low)?,
                            radius,
                        },
                        (sweep > 0.0) == along,
                        setup.sinusoid(center, radius),
                    ),
                    PSeg::Circle { center, radius } => (
                        Surface::Cylinder {
                            frame: setup.level_frame(center, setup.low)?,
                            radius,
                        },
                        along,
                        setup.sinusoid(center, radius),
                    ),
                    PSeg::Spline { reversed } => {
                        let span = spline.as_ref().expect("a spline footprint segment");
                        let whole = crate::topology::SplineSpan::whole(span.curve().clone());
                        (
                            crate::topology::spline_wall(
                                &whole,
                                setup.frame,
                                setup.low,
                                setup.high,
                            )?
                            .0,
                            reversed != along,
                            [0.0; 3],
                        )
                    }
                    PSeg::Line(..) => unreachable!("a planar wall"),
                };
                let face_sense = sense(ccw);
                // A point's height on the plane above the wall's low end.
                let plane_v = |p: Point2| rational_f64(&(setup.g_exact(p) - q(setup.low)));
                let line = |from: (f64, f64), to: (f64, f64)| Curve2::LineSegment {
                    start: Point2::new(from.0, from.1),
                    end: Point2::new(to.0, to.1),
                };
                // A crease part's pcurve between two chain indices.
                let crease_pcurve = |from: usize, to: usize, rel: i8| -> Result<Curve2> {
                    let (u0, u1) = if ring && chain.points.is_empty() {
                        if from < to {
                            (0.0, TAU)
                        } else {
                            (TAU, 0.0)
                        }
                    } else {
                        (chain.u[from], chain.u[to])
                    };
                    Ok(if rel >= 0 {
                        line((u0, v_of(wk)), (u1, v_of(wk)))
                    } else if let Some(span) = &spline {
                        // Its ends at their vertices' heights, as the edge's.
                        let ends = [from, to].map(|k| {
                            let (_, key, p) = chain.points[k];
                            (p, v_of_w(heights[&key], setup.low, height))
                        });
                        let v = |p: Point2| {
                            ends.iter()
                                .find(|e| e.0 == p)
                                .map_or_else(|| plane_v(p), |e| e.1)
                        };
                        super::spline::wall_pcurve(&super::spline::piece(span, u0, u1)?, &v)?
                    } else {
                        Curve2::Sinusoid {
                            start: u0,
                            sweep: u1 - u0,
                            a: a_coef,
                        }
                    })
                };
                let flat_pcurve = |forward: bool| -> Curve2 {
                    let (u0, u1) = if ring {
                        (0.0, TAU)
                    } else {
                        (chain.u[0], chain.u[last])
                    };
                    if forward {
                        line((u0, v_of(wf)), (u1, v_of(wf)))
                    } else {
                        line((u1, v_of(wf)), (u0, v_of(wf)))
                    }
                };
                let flat_fin = |forward: bool| Fin {
                    edge: flat_edge,
                    sense: sense(forward),
                    pcurve: flat_pcurve(forward),
                    enclosure: None,
                };
                let crease_fins = |forward: bool| -> Result<Vec<Fin>> {
                    if chain.crease.is_empty() {
                        return Ok(Vec::new());
                    }
                    crease_in(forward)
                        .into_iter()
                        .map(|(edge, f, from, to, rel)| {
                            let (from, to) = if ring && chain.points.is_empty() {
                                if f {
                                    (0, 1)
                                } else {
                                    (1, 0)
                                }
                            } else {
                                (from, to)
                            };
                            Ok(Fin {
                                edge,
                                sense: sense(f),
                                pcurve: crease_pcurve(from, to, rel)?,
                                enclosure: None,
                            })
                        })
                        .collect()
                };
                if ring {
                    // Two loops winding about the axis in opposite senses.
                    let turns = if along { 1 } else { -1 };
                    let (lower_fins, upper_fins) = if lower {
                        (vec![flat_fin(along)], crease_fins(!along)?)
                    } else {
                        (crease_fins(along)?, vec![flat_fin(!along)])
                    };
                    add_face(
                        surface,
                        face_sense,
                        vec![(lower_fins, [turns, 0]), (upper_fins, [-turns, 0])],
                        occ,
                        &mut parts,
                        &mut plans,
                    );
                    continue;
                }
                let (near, far) = if along { (0, last) } else { (last, 0) };
                let vertical = |index: usize, up: bool| -> Option<Fin> {
                    let (fk, ck, _) = chain.points[index];
                    let VKey::Flat(id) = fk else { return None };
                    let edge = *eids.get(&EKey::Vertical(id))?;
                    let (wa, wb) = (heights[&fk], heights[&ck]);
                    let (lo, hi) = if wa <= wb { (wa, wb) } else { (wb, wa) };
                    let (v_lo, v_hi) =
                        (v_of_w(lo, setup.low, height), v_of_w(hi, setup.low, height));
                    let u = chain.u[index];
                    Some(Fin {
                        edge,
                        sense: sense(up),
                        pcurve: if up {
                            line((u, v_lo), (u, v_hi))
                        } else {
                            line((u, v_hi), (u, v_lo))
                        },
                        enclosure: None,
                    })
                };
                let mut fins = Vec::new();
                if lower {
                    fins.push(flat_fin(along));
                    fins.extend(vertical(far, true));
                    fins.extend(crease_fins(!along)?);
                    fins.extend(vertical(near, false));
                } else {
                    fins.extend(crease_fins(along)?);
                    fins.extend(vertical(far, true));
                    fins.push(flat_fin(!along));
                    fins.extend(vertical(near, false));
                }
                add_face(
                    surface,
                    face_sense,
                    vec![(fins, [0, 0])],
                    occ,
                    &mut parts,
                    &mut plans,
                );
            }
        }
    }
    // One solid region inside the void.
    let faces: Vec<FaceId> = (0..parts.faces.len()).map(FaceId).collect();
    parts.shells = vec![
        Shell {
            region: RegionId(1),
            sides: faces.iter().map(|f| (*f, FaceSide::Front)).collect(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        },
        Shell {
            region: RegionId(0),
            sides: faces.iter().map(|f| (*f, FaceSide::Back)).collect(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        },
    ];
    parts.regions = vec![
        Region {
            kind: RegionKind::Void,
            shells: vec![ShellId(1)],
        },
        Region {
            kind: RegionKind::Solid,
            shells: vec![ShellId(0)],
        },
    ];
    plans.push((Slot::Region(RegionId(1)), Occ::Part(Key::Region, 0)));
    Ok((parts, plans))
}

/// A height's `v` on a wall's cylinder (its frame at the low end).
fn v_of_w(w: f64, low: f64, height: f64) -> f64 {
    if w == low {
        0.0
    } else if w == low + height {
        height
    } else {
        w - low
    }
}

/// The input prism's entities by key.
pub(super) fn resolve(
    solid: &Topology,
    profile: &Profile,
    start_is_low: bool,
) -> Result<BTreeMap<Key, EntityId>> {
    let label_of = |b: usize, e: ProfileElement| -> Parent {
        let boundary = profile.boundaries().nth(b).expect("a boundary");
        match boundary.labels() {
            Some(l) => Parent::Label(match e {
                ProfileElement::Boundary => l.boundary,
                ProfileElement::Segment(j) => l.segments[j as usize],
                ProfileElement::Vertex(j) => l.vertices[j as usize],
            }),
            None => Parent::Profile {
                boundary: b as u32,
                element: e,
            },
        }
    };
    let mut by_origin: BTreeMap<(EntityKind, Role, Vec<Parent>), EntityId> = BTreeMap::new();
    let mut singles: BTreeMap<(EntityKind, Role), EntityId> = BTreeMap::new();
    for (id, _) in solid.ids() {
        let d = solid.derivation(id).expect("a derivation");
        by_origin.insert((d.entity, d.role, d.parents.clone()), id);
        singles.insert((d.entity, d.role), id);
    }
    let level = |at_low: bool, low: Role, high: Role| {
        if at_low == start_is_low {
            low
        } else {
            high
        }
    };
    let mut out = BTreeMap::new();
    let missing = || Error::InvalidTopology("a prism entity without its profile element");
    for at_low in [true, false] {
        let cap = level(at_low, Role::StartCap, Role::EndCap);
        out.insert(
            Key::Cap(at_low),
            *singles.get(&(EntityKind::Face, cap)).ok_or_else(missing)?,
        );
    }
    out.insert(
        Key::Region,
        *singles
            .get(&(EntityKind::Region, Role::Region))
            .ok_or_else(missing)?,
    );
    for (b, boundary) in profile.boundaries().enumerate() {
        let (segments, vertices) = match &boundary.kind {
            BoundaryKind::Circle { .. } => (1, 0),
            BoundaryKind::Polygon(p) => (p.len(), p.len()),
            BoundaryKind::Path { points, .. } => (points.len(), points.len()),
        };
        let find = |entity: EntityKind, role: Role, e: ProfileElement| {
            by_origin
                .get(&(entity, role, vec![label_of(b, e)]))
                .copied()
                .ok_or_else(missing)
        };
        for j in 0..segments {
            let seg = ProfileElement::Segment(j as u32);
            out.insert(Key::Wall(b, j), find(EntityKind::Face, Role::Wall, seg)?);
            for at_low in [true, false] {
                let role = level(at_low, Role::BottomEdge, Role::TopEdge);
                out.insert(
                    Key::CapEdge(at_low, b, j),
                    find(EntityKind::Edge, role, seg)?,
                );
            }
        }
        for j in 0..vertices {
            let vertex = ProfileElement::Vertex(j as u32);
            out.insert(
                Key::Vertical(b, j),
                find(EntityKind::Edge, Role::Vertical, vertex)?,
            );
            for at_low in [true, false] {
                let role = level(at_low, Role::BottomVertex, Role::TopVertex);
                out.insert(
                    Key::CapVertex(at_low, b, j),
                    find(EntityKind::Vertex, role, vertex)?,
                );
            }
        }
    }
    Ok(out)
}

/// The derivations of every piece's slots and the history's relations:
/// an input entity whole in exactly one piece keeps its derivation; one in
/// several pieces, or in parts, is split, its children ordered by piece and
/// part; new entities are generated with ordinals by piece.
#[allow(clippy::type_complexity)]
pub(super) fn name(
    operation: OperationId,
    keys: &BTreeMap<Key, EntityId>,
    input: &Topology,
    plans: &[Vec<(Slot, Occ)>],
) -> (
    Vec<Vec<(Slot, Derivation)>>,
    Vec<crate::history::Relation>,
    BTreeSet<EntityId>,
) {
    use crate::history::Relation;
    let mut wholes: BTreeMap<Key, usize> = BTreeMap::new();
    let mut children: BTreeMap<Key, Vec<(usize, usize)>> = BTreeMap::new();
    for (k, named) in plans.iter().enumerate() {
        for (_, occ) in named {
            match occ {
                Occ::Whole(key) => {
                    *wholes.entry(*key).or_default() += 1;
                    children.entry(*key).or_default().push((k, 0));
                }
                Occ::Part(key, local) => {
                    children.entry(*key).or_default().push((k, local + 1));
                }
                Occ::New(..) => {}
            }
        }
    }
    let kept: BTreeSet<Key> = wholes
        .iter()
        .filter(|(key, n)| **n == 1 && children[*key].len() == 1)
        .map(|(key, _)| *key)
        .collect();
    for list in children.values_mut() {
        list.sort();
        list.dedup();
    }
    let mut relations = Vec::new();
    let mut out = Vec::new();
    let mut split: BTreeMap<EntityId, Vec<(u32, EntityId)>> = BTreeMap::new();
    for (k, named) in plans.iter().enumerate() {
        let mut derivations = Vec::new();
        let mut local = 0usize;
        for (slot, occ) in named {
            let d = match occ {
                Occ::Whole(key) if kept.contains(key) => {
                    input.derivation(keys[key]).expect("a derivation").clone()
                }
                Occ::Whole(key) | Occ::Part(key, _) => {
                    let mark = match occ {
                        Occ::Part(_, l) => l + 1,
                        _ => 0,
                    };
                    let from = keys[key];
                    let parent = input.derivation(from).expect("a derivation");
                    let ordinal = children[key].binary_search(&(k, mark)).expect("a child") as u32;
                    let d = Derivation {
                        operation,
                        kind: OperationKind::PlaneSplit,
                        entity: parent.entity,
                        role: parent.role,
                        ordinal,
                        parents: vec![Parent::Entity(from)],
                    };
                    split.entry(from).or_default().push((ordinal, d.id()));
                    d
                }
                Occ::New(from, role) => {
                    let entity = match slot {
                        Slot::Vertex(_) => EntityKind::Vertex,
                        Slot::Edge(_) => EntityKind::Edge,
                        Slot::Face(_) => EntityKind::Face,
                        Slot::Region(_) => EntityKind::Region,
                    };
                    let d = Derivation {
                        operation,
                        kind: OperationKind::PlaneSplit,
                        entity,
                        role: *role,
                        ordinal: piece_ordinal(k, local),
                        parents: from.iter().map(|key| Parent::Entity(keys[key])).collect(),
                    };
                    local += 1;
                    relations.push(Relation::Generated {
                        from: d.parents.clone(),
                        to: d.id(),
                        role: *role,
                    });
                    d
                }
            };
            derivations.push((*slot, d));
        }
        out.push(derivations);
    }
    for (from, mut into) in split {
        into.sort();
        into.dedup();
        relations.push(Relation::Split {
            from,
            into: into.into_iter().map(|(_, id)| id).collect(),
        });
    }
    let kept_ids: BTreeSet<EntityId> = kept.iter().map(|k| keys[k]).collect();
    (out, relations, kept_ids)
}

/// The parts of a piece as a topology with external ids (a rebuild in
/// another frame, whose ids the caller restores).
pub(super) fn unnamed(parts: TopologyParts, tolerance: Tolerance) -> Result<Topology> {
    Topology::from_parts(parts.with_measured_enclosures(), tolerance)
        .map_err(|issues| Error::InvalidTopology(issue_name(&issues)))
}

pub(super) fn issue_name(issues: &[crate::topology::Issue]) -> &'static str {
    issues
        .first()
        .map_or("a split piece's identity", |i| i.kind.name())
}

/// Axis-aligned bounds of a body from its edges (planar and cylindrical
/// faces reach their extremes on their boundaries): lines by their ends,
/// arcs and ellipse arcs by their ends and the extremes inside them.
pub(crate) fn edge_bounds(t: &Topology) -> crate::Bounds3 {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    let mut add = |p: Point3| {
        for (i, x) in p.to_array().into_iter().enumerate() {
            lo[i] = lo[i].min(x);
            hi[i] = hi[i].max(x);
        }
    };
    for edge in t.edges() {
        add(edge.curve.point(0.0));
        add(edge.curve.point(1.0));
        let (frame, rx, ry, start, sweep) = match &edge.curve {
            Curve3::Circle { frame, radius } => (frame, *radius, *radius, 0.0, TAU),
            Curve3::CircularArc {
                frame,
                radius,
                start_angle,
                sweep_angle,
            } => (frame, *radius, *radius, *start_angle, *sweep_angle),
            Curve3::EllipseArc {
                frame,
                major,
                minor,
                start_angle,
                sweep_angle,
            } => (frame, *major, *minor, *start_angle, *sweep_angle),
            _ => continue,
        };
        let (x, y) = (frame.x().to_array(), frame.y().to_array());
        for i in 0..3 {
            // Extremes of rx cos t x_i + ry sin t y_i at t = atan2(ry y_i, rx x_i) (+ pi).
            let t0 = (ry * y[i]).atan2(rx * x[i]);
            for t in [t0, t0 + std::f64::consts::PI] {
                let mut f = (t - start) / sweep;
                f -= f.floor();
                // Every alias of t inside the arc.
                for shift in [-1.0, 0.0, 1.0] {
                    let g = f + shift * TAU / sweep.abs();
                    if (0.0..=1.0).contains(&g) {
                        add(edge.curve.point(g));
                    }
                }
            }
        }
    }
    let pad = |x: f64| x.abs() * 4.0 * f64::EPSILON;
    crate::Bounds3 {
        min: Point3::new(lo[0] - pad(lo[0]), lo[1] - pad(lo[1]), lo[2] - pad(lo[2])),
        max: Point3::new(hi[0] + pad(hi[0]), hi[1] + pad(hi[1]), hi[2] + pad(hi[2])),
    }
}

/// How a split piece was made, kept so a rigid motion rebuilds it exactly:
/// the prism's profile, the plane in its frame, which piece, and the
/// piece's footprint and side for classification.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Clipped {
    profile: Profile,
    plane: [R; 4],
    /// The piece's index in the split's order.
    index: usize,
    footprint: Profile,
    /// The sign of the plane's function inside the piece.
    sign: i8,
}

impl Clipped {
    pub(crate) fn tolerance(&self) -> Tolerance {
        self.profile.tolerance()
    }

    /// The same piece of the same split in another frame, its ids external
    /// (the caller restores them).
    pub(crate) fn rebuilt(
        &self,
        operation: OperationId,
        frame: Frame3,
        start: f64,
        end: f64,
    ) -> Result<Solid> {
        self.rebuilt_with(operation, frame, start, end, None)
    }

    /// [`Clipped::rebuilt`] with its mass properties when known (a rigid
    /// motion's: its source's, moved), which are costly to enclose again.
    pub(crate) fn rebuilt_with(
        &self,
        operation: OperationId,
        frame: Frame3,
        start: f64,
        end: f64,
        known: Option<crate::MassProperties>,
    ) -> Result<Solid> {
        let mut built = pieces(
            &self.profile,
            frame,
            start,
            end,
            &self.plane,
            Some(self.index),
        )?;
        let Some(piece) = built.pop() else {
            return Err(Error::InvalidTopology("a split piece rebuilt differently"));
        };
        let topology = unnamed(piece.parts, self.tolerance())?;
        let mass = match known {
            Some(mass) => mass,
            None => topology
                .mass_enclosure()
                .ok_or(Error::Unrepresentable("a split piece's mass properties"))?
                .midpoints(),
        };
        let bounds = edge_bounds(&topology);
        Ok(Solid {
            construction: super::Construction::Clipped(Box::new(self.clone())),
            frame,
            start,
            end,
            topology,
            mass,
            bounds,
            operation,
        })
    }

    /// Inside where the footprint, the prism's height range and the plane's
    /// side all hold; on the boundary within the resolution of any of them.
    pub(crate) fn classify(
        &self,
        [u, v, w]: [f64; 3],
        low: f64,
        high: f64,
        tolerance: Tolerance,
    ) -> Result<crate::Location> {
        use crate::Location;
        for x in [u, v, w] {
            crate::math::finite(x, "coordinate")?;
        }
        let tol = tolerance.linear();
        let [a, b, c, d] = &self.plane;
        let f = a * q(u) + b * q(v) + c * q(w) + d;
        let norm2 = a * a + b * b + c * c;
        let near_plane = &f * &f <= q(tol) * q(tol) * &norm2;
        let wrong_side = sign_of(&f) * self.sign < 0;
        let footprint = self.footprint.classify(Point2::new(u, v))?;
        let below = crate::decide::sum_gt(&[low, -tol], &[w]);
        let above = crate::decide::sum_gt(&[w], &[high, tol]);
        if footprint == Location::Outside || below || above || (wrong_side && !near_plane) {
            return Ok(Location::Outside);
        }
        let at_end = [low, high].iter().any(|e| {
            crate::decide::sum_le(&[w, -e], &[tol]) && crate::decide::sum_le(&[*e, -w], &[tol])
        });
        Ok(if footprint == Location::Boundary || at_end || near_plane {
            Location::Boundary
        } else {
            Location::Inside
        })
    }
}

impl Solid {
    /// S8a.2: the pieces of a prism split by a plane oblique to its axis,
    /// named by provenance, with their `PlaneSplit` history.
    pub(super) fn split_oblique(
        &self,
        context: &crate::solid::Context,
        profile: &Profile,
        plane: [R; 4],
    ) -> Result<(Vec<(Side, Solid)>, crate::history::History)> {
        use crate::history::{History, Relation};
        let operation = context.operation;
        let built = pieces(profile, self.frame, self.start, self.end, &plane, None)?;
        let keys = resolve(&self.topology, profile, self.start < self.end)?;
        let plans: Vec<Vec<(Slot, Occ)>> = built.iter().map(|b| b.plans.clone()).collect();
        let (derivations, mut relations, kept) = name(operation, &keys, &self.topology, &plans);
        let body = self.topology.body_id();
        let tolerance = self.resolution();
        let lower_sign: i8 = if plane[2] > zero() { -1 } else { 1 };
        let mut pieces = Vec::new();
        for (k, (piece, named)) in built.into_iter().zip(derivations).enumerate() {
            let body_derivation = Derivation {
                operation,
                kind: OperationKind::PlaneSplit,
                entity: EntityKind::Body,
                role: Role::Body,
                ordinal: k as u32,
                parents: vec![Parent::Entity(body)],
            };
            let topology =
                Topology::from_parts_named(piece.parts, tolerance, body_derivation, named)
                    .map_err(|issues| Error::InvalidTopology(issue_name(&issues)))?;
            let mass = topology
                .mass_enclosure()
                .ok_or(Error::Unrepresentable("a split piece's mass properties"))?
                .midpoints();
            let bounds = edge_bounds(&topology);
            let clipped = Clipped {
                profile: profile.clone(),
                plane: plane.clone(),
                index: k,
                footprint: piece.footprint,
                sign: if piece.lower { lower_sign } else { -lower_sign },
            };
            pieces.push((
                piece.side,
                Solid {
                    construction: super::Construction::Clipped(Box::new(clipped)),
                    frame: self.frame,
                    start: self.start,
                    end: self.end,
                    topology,
                    mass,
                    bounds,
                    operation,
                },
            ));
        }
        // An input entity whole in one piece keeps its id: `Unchanged` when
        // its stored geometry and bounding ids are the input's.
        let before = self.topology.entity_set(tolerance);
        let mut modified: BTreeSet<EntityId> = BTreeSet::new();
        for (_, piece) in &pieces {
            let after = piece.topology.entity_set(tolerance);
            for (id, info) in &after.entities {
                if !kept.contains(id) {
                    continue;
                }
                let old = &before.entities[id];
                if old.geometry != info.geometry || old.structure != info.structure {
                    modified.insert(*id);
                }
            }
        }
        for id in &kept {
            relations.push(if modified.contains(id) {
                Relation::Modified { from: *id, to: *id }
            } else {
                Relation::Unchanged { id: *id }
            });
        }
        relations.sort_by_cached_key(Relation::sort_key);
        relations.dedup();
        let outputs: Vec<EntityId> = pieces.iter().map(|(_, s)| s.topology.body_id()).collect();
        let history = History::new(
            operation,
            OperationKind::PlaneSplit,
            vec![body],
            outputs,
            relations,
            Vec::new(),
        );
        let mut history = history.at_level(context.level);
        {
            let mut outs: Vec<&mut Solid> = pieces.iter_mut().map(|(_, s)| s).collect();
            crate::solid::enclose::carry(&[self], &history, &mut outs);
        }
        let outs: Vec<&Solid> = pieces.iter().map(|(_, s)| s).collect();
        let maps = crate::solid::attrs::propagate(context, &[self], &mut history, &outs)?;
        for ((_, s), map) in pieces.iter_mut().zip(maps) {
            s.topology.set_attributes(map);
        }
        let outs: Vec<&Solid> = pieces.iter().map(|(_, s)| s).collect();
        crate::solid::stack::debug_check(&[self], &outs, &history);
        crate::solid::attrs::debug_check_attributes(context, &[self], &outs, &history);
        Ok((pieces, history))
    }
}
