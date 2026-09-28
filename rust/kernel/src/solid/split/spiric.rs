//! S8d.3: a whole torus split by a plane neither normal to its axis nor
//! containing it.
//!
//! In the torus's frame (major `R`, minor `r`) the plane is `F = a x + b y +
//! c z + d`; on the torus `F = C + W cos(v - psi)` with `alpha = a cos u + b
//! sin u`, `C = R alpha + d` and `(W cos psi, W sin psi) = (r alpha, r c)`, so
//! the section exists over `u` where `D = W^2 - C^2 >= 0`: a quadratic in
//! `alpha` with the roots `alpha_{1,2} = (-R d -+ r sqrt(c^2 (R^2 - r^2) +
//! d^2)) / (R^2 - r^2)`, `alpha` ranging over `[-rho, rho]`, `rho = |(a, b)|`.
//! Exactly (surds on the stored data): both roots beyond `[-rho, rho]` cut the
//! tube in two loops winding about the axis (the pieces tube bands); one root
//! inside it in one contractible loop (a cap, and the torus less that disc);
//! both inside it in two loops, winding about the tube when `C` changes sign
//! between the roots (the pieces C-shaped) and two caps otherwise
//! (`OutOfDomain`). A root at `+-rho` is a tangency (`Degenerate`).
//!
//! The loops about the axis are graphs over `u` for the whole turn, those
//! about the tube over `v`; a cap's loop is two graphs over `u` (its two
//! branches) and two over `v` (round its turning points), joined where the
//! slope is one, so no edge reaches a turning point (`Curve3::Section`). Each
//! piece is a general body on the input's torus surface with `Projection`
//! pcurves and planar cut faces. Names: the wall and the region `Split`,
//! everything else `Generated` from the wall.
use super::torus::Band;
use super::{q, zero, Side};
use crate::certified::{Interval as I, Real};
use crate::identity::Role;
use crate::topology::{
    plane_pcurve, Curve2, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId,
    Orientation, Projection, Region, RegionId, RegionKind, Shell, ShellId, Side as FaceSide, Slot,
    Spiric, Surface, TopologyParts, Vertex, VertexId,
};
use crate::{Error, Frame3, Point2, Point3, Result, Tolerance, Vec3};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::f64::consts::{PI, TAU};

/// Lifts recorded on a wall's projection pcurve.
const ANCHORS: usize = 32;

/// The sign of `k + t sqrt(r2)`, `r2 >= 0`, `t` in `{1, -1}`, exactly.
fn surd_sign(k: &R, r2: &R, t: i8) -> Ordering {
    let kk = k * k;
    match (k.cmp(&zero()), t) {
        (Ordering::Equal, _) => {
            if *r2 == zero() {
                Ordering::Equal
            } else if t > 0 {
                Ordering::Greater
            } else {
                Ordering::Less
            }
        }
        (Ordering::Greater, 1) => Ordering::Greater,
        (Ordering::Less, -1) => Ordering::Less,
        (Ordering::Greater, _) => kk.cmp(r2),
        (_, _) => r2.cmp(&kk),
    }
}

/// The sign of `x + y sqrt(e) - z sqrt(f)`, `e, f >= 0`, exactly.
fn sign_two(x: &R, y: &R, e: &R, z: &R, f: &R) -> Ordering {
    let sy = if *y > zero() { 1 } else { -1 };
    let p = if *y == zero() {
        x.cmp(&zero())
    } else {
        surd_sign(x, &(y * y * e), sy)
    };
    let qs = if *f == zero() {
        Ordering::Equal
    } else {
        z.cmp(&zero())
    };
    if p != qs {
        return match (p, qs) {
            (Ordering::Greater, _) | (Ordering::Equal, Ordering::Less) => Ordering::Greater,
            _ => Ordering::Less,
        };
    }
    if p == Ordering::Equal {
        return Ordering::Equal;
    }
    // Same sign s: P - Q has the sign s times that of P^2 - Q^2.
    let k = x * x + y * y * e - z * z * f;
    let r2 = R::from_integer(4.into()) * x * x * y * y * e;
    let t = if (*x > zero()) == (*y > zero()) {
        1
    } else {
        -1
    };
    let squares = surd_sign(&k, &r2, t);
    if p == Ordering::Greater {
        squares
    } else {
        squares.reverse()
    }
}

/// How the plane cuts the tube.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// Two loops winding about the axis.
    Bands,
    /// One contractible loop about `u = phi` (`true`) or `phi + pi`.
    Cap(bool),
    /// Two loops winding about the tube.
    Tubes,
}

/// The classification, exactly.
fn classify(major: f64, minor: f64, plane: &[R; 4]) -> Result<Kind> {
    let [a, b, c, d] = plane;
    let (big, small) = (q(major), q(minor));
    let delta = &big * &big - &small * &small;
    let e = c * c * &delta + d * d;
    let ab2 = a * a + b * b;
    let p = -&big * d;
    // alpha_k = (p + s_k r sqrt(e)) / delta against sigma rho: the sign of
    // p + s_k r sqrt(e) - sigma delta sqrt(ab2).
    let cmp = |s: i64, sigma: i64| {
        sign_two(
            &p,
            &(&small * R::from_integer(s.into())),
            &e,
            &(&delta * R::from_integer(sigma.into())),
            &ab2,
        )
    };
    let (lo_low, lo_high) = (cmp(-1, -1), cmp(-1, 1));
    let (hi_low, hi_high) = (cmp(1, -1), cmp(1, 1));
    if [lo_low, lo_high, hi_low, hi_high].contains(&Ordering::Equal) {
        return Err(Error::Degenerate(
            "a plane tangent to the torus's tube along a turning point",
        ));
    }
    let below = |o: Ordering| o == Ordering::Less;
    match (below(lo_low), below(hi_high), below(lo_high), below(hi_low)) {
        // alpha_1 < -rho and alpha_2 > rho.
        (true, false, _, _) => Ok(Kind::Bands),
        // alpha_1 < -rho < alpha_2 < rho: about the minimum of alpha.
        (true, true, _, false) => Ok(Kind::Cap(false)),
        // -rho < alpha_1 < rho < alpha_2: about the maximum.
        (false, false, true, _) => Ok(Kind::Cap(true)),
        // Both roots inside: C = R alpha + d at each, R p + delta d + s_k R r sqrt(e).
        (false, true, true, false) => {
            let k = &big * &p + &delta * d;
            let r2 = &big * &big * &small * &small * &e;
            if surd_sign(&k, &r2, -1) != surd_sign(&k, &r2, 1) {
                Ok(Kind::Tubes)
            } else {
                Err(Error::OutOfDomain("a torus cut in two caps"))
            }
        }
        _ => Err(Error::Degenerate("a plane missing or touching the tube")),
    }
}

/// The unit plane in the torus's frame, rounded.
fn unit_plane(plane: &[R; 4]) -> Result<[f64; 4]> {
    let [a, b, c, d] = plane;
    let m = I::exact(a * a + b * b + c * c).sqrt();
    let mut out = [0.0; 4];
    for (o, x) in out.iter_mut().zip([a, b, c, d]) {
        let (lo, hi) = I::exact(x.clone())
            .div(&m)
            .ok_or(Error::Degenerate("a plane's normal"))?
            .bounds_f64();
        *o = 0.5 * lo + 0.5 * hi;
    }
    Ok(out)
}

/// The section's geometry in binary64.
struct Geometry {
    frame: Frame3,
    major: f64,
    minor: f64,
    plane: [f64; 4],
}

impl Geometry {
    fn section(&self, over_v: bool, sign: f64, start: f64, sweep: f64) -> Curve3 {
        Curve3::Section(Box::new(Spiric {
            frame: self.frame,
            major: self.major,
            minor: self.minor,
            plane: self.plane,
            over_v,
            sign,
            start,
            sweep,
        }))
    }

    /// `v` on branch `sign` of the graph over `u`.
    fn v_of(&self, u: f64, sign: f64) -> f64 {
        let [a, b, c, d] = self.plane;
        let alpha = a * u.cos() + b * u.sin();
        let (x, y) = (self.minor * alpha, self.minor * c);
        let ratio = -(self.major * alpha + d) / x.hypot(y);
        y.atan2(x) + sign * ratio.clamp(-1.0, 1.0).acos()
    }

    /// `dv / du` of the section at `(u, v)`.
    fn slope(&self, u: f64, v: f64) -> f64 {
        let [a, b, c, _] = self.plane;
        let alpha = a * u.cos() + b * u.sin();
        let fu = (self.major + self.minor * v.cos()) * (-a * u.sin() + b * u.cos());
        let fv = self.minor * (-alpha * v.sin() + c * v.cos());
        -fu / fv
    }

    fn point(&self, u: f64, v: f64) -> Point3 {
        let rho = self.major + self.minor * v.cos();
        self.frame.point(
            Point2::new(rho * u.cos(), rho * u.sin()),
            self.minor * v.sin(),
        )
    }

    /// The torus's outward normal at `(u, v)`, in world coordinates.
    fn outward(&self, u: f64, v: f64) -> Vec3 {
        let ((su, cu), (sv, cv)) = (u.sin_cos(), v.sin_cos());
        let f = self.frame;
        f.x() * (cv * cu) + f.y() * (cv * su) + f.normal() * sv
    }

    fn n_world(&self) -> Vec3 {
        let [a, b, c, _] = self.plane;
        self.frame.x() * a + self.frame.y() * b + self.frame.normal() * c
    }
}

/// An edge of the section: its curve, its ends (none for a ring), and the
/// side on its left over the torus (about its outward normal).
struct SectionEdge {
    curve: Curve3,
    ends: Option<(Point3, Point3)>,
    left: Side,
}

/// The side on a curve's left at its middle: the plane's normal against the
/// torus's outward normal crossed with the tangent.
fn left_of(g: &Geometry, curve: &Curve3) -> Side {
    let (p0, p1) = (curve.point(0.5 - 1e-6), curve.point(0.5 + 1e-6));
    let [x, y, z] = g.frame.coordinates(curve.point(0.5));
    let rho = x.hypot(y);
    let u = y.atan2(x);
    let v = z.atan2(rho - g.major);
    let left = g.outward(u, v).cross(p1 - p0);
    if left.dot(g.n_world()) > 0.0 {
        Side::Above
    } else {
        Side::Below
    }
}

/// `acos(q / |(a, b)|)` of the graph over `v` at `v`.
fn canonical_acos_v(g: &Geometry, v: f64) -> f64 {
    let [a, b, c, d] = g.plane;
    let q = -(c * g.minor * v.sin() + d) / (g.major + g.minor * v.cos());
    (q / a.hypot(b)).clamp(-1.0, 1.0).acos()
}

/// A cap's four edges: the branches over `u` between their switch points and
/// the graphs over `v` round the turning points.
fn cap_edges(g: &Geometry, about_max: bool, tol: f64) -> Result<Vec<SectionEdge>> {
    let [a, b, c, d] = g.plane;
    let rho = a.hypot(b);
    let phi = b.atan2(a);
    let (big, small) = (g.major, g.minor);
    // The root of D inside (-rho, rho), rounded.
    let delta = big * big - small * small;
    let root = (small * (c * c * delta + d * d).sqrt()) / delta;
    let alpha_e = if about_max {
        -big * d / delta - root
    } else {
        -big * d / delta + root
    };
    let ratio = (alpha_e / rho).clamp(-1.0, 1.0);
    let (centre, half) = if about_max {
        (phi, ratio.acos())
    } else {
        (phi + PI, (-ratio).acos())
    };
    if half * (big + small) <= tol {
        return Err(Error::Degenerate("a cap within the resolution of a point"));
    }
    let (ul, ur) = (centre - half, centre + half);
    // On each branch, the switch points nearest each end: where the slope
    // falls to one, found from the end inward.
    let samples = 4096;
    let mut switches = [[0.0f64; 2]; 2];
    for (bi, sign) in [1.0f64, -1.0].into_iter().enumerate() {
        let steep = |u: f64| g.slope(u, g.v_of(u, sign)).abs() > 1.0;
        let at = |k: usize| ul + (ur - ul) * k as f64 / samples as f64;
        // From the left end inward.
        let first = (1..samples)
            .find(|&k| !steep(at(k)))
            .ok_or(Error::Degenerate("a cap's branch steep throughout"))?;
        let last = (1..samples)
            .rev()
            .find(|&k| !steep(at(k)))
            .ok_or(Error::Degenerate("a cap's branch steep throughout"))?;
        let refine = |mut flat: f64, mut sharp: f64| {
            for _ in 0..60 {
                let m = 0.5 * (flat + sharp);
                if steep(m) {
                    sharp = m;
                } else {
                    flat = m;
                }
            }
            flat
        };
        switches[bi] = [
            refine(at(first), at(first - 1)),
            refine(at(last), at(last + 1)),
        ];
        if switches[bi][1] - switches[bi][0] <= 0.0 {
            return Err(Error::Degenerate("a cap's branch without a gentle stretch"));
        }
    }
    let mut out = Vec::new();
    // The branches over u, increasing.
    for (bi, sign) in [1.0f64, -1.0].into_iter().enumerate() {
        let [s0, s1] = switches[bi];
        let curve = g.section(false, sign, s0, s1 - s0);
        let ends = (g.point(s0, g.v_of(s0, sign)), g.point(s1, g.v_of(s1, sign)));
        out.push(SectionEdge {
            left: left_of(g, &curve),
            curve,
            ends: Some(ends),
        });
    }
    // Round each end over v: from the + branch's switch to the - branch's,
    // through the turning point.
    for (end, uend) in [(0usize, ul), (1usize, ur)] {
        let (up, um) = (switches[0][end], switches[1][end]);
        let (vp, vm) = (g.v_of(up, 1.0), g.v_of(um, -1.0));
        // The turning point's v, and the lifts of the two ends about it.
        let vt = g.v_of(uend, 1.0);
        let near = |x: f64| x + TAU * ((vt - x) / TAU).round();
        let (vp, vm) = (near(vp), near(vm));
        let sign = if (uend - phi + PI).rem_euclid(TAU) - PI >= 0.0 {
            1.0
        } else {
            -1.0
        };
        let curve = g.section(true, sign, vp, vm - vp);
        out.push(SectionEdge {
            left: left_of(g, &curve),
            ends: Some((g.point(up, vp), g.point(um, vm))),
            curve,
        });
    }
    Ok(out)
}

/// A plane frame on the section's plane (its origin the torus centre's
/// foot).
fn cut_frame(g: &Geometry, tolerance: Tolerance) -> Result<Frame3> {
    let [a, b, c, d] = g.plane;
    let foot = g.frame.point(Point2::new(-d * a, -d * b), -d * c);
    let n = g.n_world();
    let hint = if n.cross(g.frame.x()).length() > 0.5 {
        g.frame.x()
    } else {
        g.frame.y()
    };
    Frame3::new(foot, n, hint, tolerance)
}

fn sense(forward: bool) -> Orientation {
    if forward {
        Orientation::Forward
    } else {
        Orientation::Reversed
    }
}

/// Chain edge uses into loops (as the conic split does).
fn chain(uses: &[(EdgeId, bool)], parts: &TopologyParts) -> Result<Vec<Vec<(EdgeId, bool)>>> {
    let ends = |&(e, forward): &(EdgeId, bool)| {
        let edge = &parts.edges[e.0];
        if forward {
            (edge.start, edge.end)
        } else {
            (edge.end, edge.start)
        }
    };
    let mut used = vec![false; uses.len()];
    let mut loops = Vec::new();
    for s in 0..uses.len() {
        if used[s] {
            continue;
        }
        used[s] = true;
        let mut list = vec![uses[s]];
        let (first, mut at) = ends(&uses[s]);
        if first.is_none() {
            loops.push(list);
            continue;
        }
        while let Some(next) = (0..uses.len()).find(|&i| !used[i] && ends(&uses[i]).0 == at) {
            used[next] = true;
            list.push(uses[next]);
            at = ends(&uses[next]).1;
        }
        if at != first {
            return Err(Error::InvalidTopology("an open loop in a torus piece"));
        }
        loops.push(list);
    }
    Ok(loops)
}

/// Both pieces (below first).
#[allow(clippy::too_many_lines)]
pub(super) fn pieces(
    frame: Frame3,
    major: f64,
    minor: f64,
    plane: &[R; 4],
    tolerance: Tolerance,
) -> Result<Vec<Band>> {
    let tol = tolerance.linear();
    let kind = classify(major, minor, plane)?;
    let g = Geometry {
        frame,
        major,
        minor,
        plane: unit_plane(plane)?,
    };
    let [a, b, _, d] = g.plane;
    let rho = a.hypot(b);
    // The pieces' thickness: F's extremes over the solid torus.
    let reach = major * rho + minor;
    if d + reach <= tol || d - reach >= -tol {
        return Err(Error::Degenerate("a piece thinner than the resolution"));
    }
    // The section's edges.
    let edges: Vec<SectionEdge> = match kind {
        Kind::Bands => [1.0f64, -1.0]
            .into_iter()
            .map(|sign| {
                let curve = g.section(false, sign, 0.0, TAU);
                SectionEdge {
                    left: left_of(&g, &curve),
                    curve,
                    ends: None,
                }
            })
            .collect(),
        Kind::Tubes => [1.0f64, -1.0]
            .into_iter()
            .map(|sign| {
                let curve = g.section(true, sign, 0.0, TAU);
                SectionEdge {
                    left: left_of(&g, &curve),
                    curve,
                    ends: None,
                }
            })
            .collect(),
        Kind::Cap(about_max) => cap_edges(&g, about_max, tol)?,
    };
    // Where the loops are separated by less than the resolution the pieces
    // are not certified: a coarse check on the rings' closest approach.
    if !matches!(kind, Kind::Cap(_)) {
        let (p, q) = (&edges[0].curve, &edges[1].curve);
        let close = (0..256).any(|k| {
            let t = k as f64 / 256.0;
            (0..256).any(|j| p.point(t).distance(q.point(j as f64 / 256.0)) <= 16.0 * tol)
        });
        if close {
            return Err(Error::Degenerate(
                "a plane within the resolution of a torus's tangent",
            ));
        }
    }
    let cut = cut_frame(&g, tolerance)?;
    let wall_surface = Surface::Torus {
        frame,
        major,
        minor,
    };
    let mut out = Vec::new();
    for side in [Side::Below, Side::Above] {
        let mut parts = TopologyParts::default();
        let mut plans: Vec<(Slot, Option<Role>)> = Vec::new();
        // Vertices by position (the cap's four switch points).
        let mut vertices: Vec<(Point3, VertexId)> = Vec::new();
        let mut vertex =
            |p: Point3, parts: &mut TopologyParts, plans: &mut Vec<(Slot, Option<Role>)>| {
                if let Some((_, id)) = vertices.iter().find(|(q, _)| q.distance(p) <= tol) {
                    return *id;
                }
                let id = VertexId(parts.vertices.len());
                parts.vertices.push(Vertex {
                    position: p,
                    enclosure: None,
                });
                plans.push((Slot::Vertex(id), Some(Role::CutVertex)));
                vertices.push((p, id));
                id
            };
        let mut section_ids = Vec::new();
        for e in &edges {
            let (start, end) = match e.ends {
                Some((p, q)) => (
                    Some(vertex(p, &mut parts, &mut plans)),
                    Some(vertex(q, &mut parts, &mut plans)),
                ),
                None => (None, None),
            };
            let id = EdgeId(parts.edges.len());
            parts.edges.push(Edge {
                start,
                end,
                curve: e.curve.clone(),
                fins: Vec::new(),
            });
            plans.push((Slot::Edge(id), Some(Role::CutEdge)));
            section_ids.push((id, e.left));
        }
        // The wall's uses: forward where this side lies on the edge's left.
        let wall_uses: Vec<(EdgeId, bool)> = section_ids
            .iter()
            .map(|&(id, left)| (id, left == side))
            .collect();
        let cut_uses: Vec<(EdgeId, bool)> = section_ids
            .iter()
            .map(|&(id, left)| (id, left != side))
            .collect();
        let add_face = |surface: Surface,
                        face_sense: Orientation,
                        loops: Vec<(Vec<Fin>, [i32; 2])>,
                        plan: Option<Role>,
                        parts: &mut TopologyParts,
                        plans: &mut Vec<(Slot, Option<Role>)>| {
            let mut ids = Vec::new();
            for (fins, winding) in loops {
                let mut list = Vec::new();
                for fin in fins {
                    parts.edges[fin.edge.0].fins.push(FinId(parts.fins.len()));
                    list.push(FinId(parts.fins.len()));
                    parts.fins.push(fin);
                }
                parts.loops.push(Loop::Edges {
                    fins: list,
                    winding,
                });
                ids.push(LoopId(parts.loops.len() - 1));
            }
            parts.faces.push(Face {
                surface,
                sense: face_sense,
                loops: ids,
                front: ShellId(0),
                back: ShellId(1),
                enclosure: None,
            });
            plans.push((Slot::Face(FaceId(parts.faces.len() - 1)), plan));
        };
        // The wall: each loop walked on the cover.
        let mut wall_loops = Vec::new();
        for list in chain(&wall_uses, &parts)? {
            let mut fins = Vec::new();
            let mut first: Option<Point2> = None;
            let mut at: Option<Point2> = None;
            for &(edge, forward) in &list {
                let curve = parts.edges[edge.0].curve.clone();
                let p = curve.point(if forward { 0.0 } else { 1.0 });
                let mut uv = Projection::inverse(&wall_surface, p)
                    .ok_or(Error::InvalidTopology("a torus without an inverse"))?;
                if let Some(prev) = at {
                    uv.x += TAU * ((prev.x - uv.x) / TAU).round();
                    uv.y += TAU * ((prev.y - uv.y) / TAU).round();
                } else if let Curve3::Section(sec) = &curve {
                    // A ring's lift: the piece lies between the two rings on
                    // the cover, the + branch above the - one on the upper
                    // side (the plane's side where F > 0 holds `psi`, or
                    // `phi`) and below the - one lifted a turn on the lower.
                    if sec.sweep.abs() == TAU {
                        let lift = if side == Side::Below && sec.sign < 0.0 {
                            TAU
                        } else {
                            0.0
                        };
                        if sec.over_v {
                            let u = b.atan2(a) + sec.sign * canonical_acos_v(&g, uv.y);
                            uv.x = u + lift;
                        } else {
                            uv.y = g.v_of(uv.x, sec.sign) + lift;
                        }
                    }
                }
                let projection =
                    Projection::new(curve, wall_surface.clone(), !forward, uv, ANCHORS)
                        .ok_or(Error::Degenerate("a section turning too fast on the torus"))?;
                first.get_or_insert(uv);
                at = Some(projection.point(1.0));
                fins.push(Fin {
                    edge,
                    sense: sense(forward),
                    pcurve: Curve2::Projection(Box::new(projection)),
                    enclosure: None,
                });
            }
            let (s0, s1) = (first.expect("a fin"), at.expect("a fin"));
            let turns = [
                ((s1.x - s0.x) / TAU).round() as i32,
                ((s1.y - s0.y) / TAU).round() as i32,
            ];
            wall_loops.push((fins, turns));
        }
        add_face(
            wall_surface.clone(),
            Orientation::Forward,
            wall_loops,
            None,
            &mut parts,
            &mut plans,
        );
        // The cut faces: one annulus for bands (the outer loop first), one
        // disc per loop otherwise.
        let mut cut_loops: Vec<(Vec<Fin>, f64)> = Vec::new();
        for list in chain(&cut_uses, &parts)? {
            let fins: Vec<Fin> = list
                .iter()
                .map(|&(edge, forward)| Fin {
                    edge,
                    sense: sense(forward),
                    pcurve: plane_pcurve(&parts.edges[edge.0].curve, sense(forward), cut),
                    enclosure: None,
                })
                .collect();
            // The loop's area in the plane (its size orders the annulus).
            let mut area = 0.0;
            for &(edge, forward) in &list {
                let e = &parts.edges[edge.0].curve;
                let n = 64;
                for k in 0..n {
                    let (t0, t1) = (k as f64 / n as f64, (k + 1) as f64 / n as f64);
                    let (t0, t1) = if forward {
                        (t0, t1)
                    } else {
                        (1.0 - t0, 1.0 - t1)
                    };
                    let [x0, y0, _] = cut.coordinates(e.point(t0));
                    let [x1, y1, _] = cut.coordinates(e.point(t1));
                    area += 0.5 * (x0 * y1 - x1 * y0);
                }
            }
            cut_loops.push((fins, area));
        }
        let face_sense = if side == Side::Below {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        if kind == Kind::Bands {
            cut_loops.sort_by(|x, y| y.1.abs().total_cmp(&x.1.abs()));
            add_face(
                Surface::Plane(cut),
                face_sense,
                cut_loops.into_iter().map(|(f, _)| (f, [0, 0])).collect(),
                Some(Role::CutFace),
                &mut parts,
                &mut plans,
            );
        } else {
            for (fins, _) in cut_loops {
                add_face(
                    Surface::Plane(cut),
                    face_sense,
                    vec![(fins, [0, 0])],
                    Some(Role::CutFace),
                    &mut parts,
                    &mut plans,
                );
            }
        }
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
        plans.push((Slot::Region(RegionId(1)), None));
        out.push(Band { side, parts, plans });
    }
    Ok(out)
}
