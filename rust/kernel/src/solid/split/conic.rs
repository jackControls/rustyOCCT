//! S8d.2: a cone, frustum or sphere zone split by a plane neither normal to
//! its axis nor containing it.
//!
//! In the solid's frame the plane is `F = a u + b v + c w + d`. Whether it
//! misses, touches or crosses each end circle is exact (`(c w + d)^2`
//! against `r^2 (a^2 + b^2)` on the stored radius), as is whether it passes
//! through an apex, a pole or a frustum's virtual apex, and which conic it
//! cuts from a cone (the sign of `H^2 |m|^2 - (a^2 + b^2) (H^2 + (r1 -
//! r0)^2)`). The section is an ellipse, hyperbola or parabola on a cone,
//! the frustum's rulings through its virtual apex, a circle on a sphere: an
//! explicit curve (D13), cut at the rim crossings into the arcs inside the
//! solid, with `Projection` pcurves on the wall and exact ones on the cut
//! plane. Each piece is a general body: its part of the wall (bounded by its
//! rim arcs or rings and the section arcs, round an apex or pole on its
//! side), its end discs or their parts closed by chords, and the cut face.
//!
//! A closed section touching a rim keeps one vertex there (the wall a
//! bigon), as S8a.2's touch; a plane touching one rim while crossing the
//! other would pinch a wall mid-loop and is `Degenerate`, as is one within
//! the resolution of a rim's tangent, an apex or a pole, or cutting a piece
//! thinner than it. A plane through an apex or pole off the axis cuts rulings
//! or a circle through a pole and is `OutOfDomain`.
//!
//! Names follow provenance: the wall, the region, and each crossed rim and
//! disc are `Split` into one child per piece (below first); a rim, disc,
//! apex or pole whole in one piece keeps its id (`Modified` when a touch
//! gives it a vertex); crossing and touch vertices are `Generated` from
//! their rim, chords from their disc, section edges from the wall, and the
//! cut face from the wall and the discs it crosses.
use super::meridian::{ends, Built, EndSpec, Plan, Primitive};
use super::{q, zero, Side};
use crate::certified::{Interval as I, Real};
use crate::identity::Role;
use crate::topology::{
    plane_pcurve, Curve2, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId,
    Orientation, Projection, Region, RegionId, RegionKind, Shell, ShellId, Side as FaceSide, Slot,
    Surface, TopologyParts, Vertex, VertexId,
};
use crate::{Error, Frame3, Point2, Point3, Result, Tolerance, Vec3};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

/// Lifts recorded on a wall's projection pcurve.
const ANCHORS: usize = 32;

/// How the plane meets an end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Rim {
    /// An apex or pole (a zero radius), off the plane.
    Point,
    Miss,
    Touch,
    Cross,
}

/// The section's curve: a conic (or a circle) by its parameter, or the two
/// rulings through a frustum's virtual apex.
enum Section {
    Ellipse {
        frame: Frame3,
        major: f64,
        minor: f64,
    },
    Hyperbola {
        frame: Frame3,
        major: f64,
        minor: f64,
    },
    Parabola {
        frame: Frame3,
        focal: f64,
    },
    Circle {
        frame: Frame3,
        radius: f64,
    },
    Rulings,
}

impl Section {
    fn periodic(&self) -> bool {
        matches!(self, Self::Ellipse { .. } | Self::Circle { .. })
    }

    /// The parameter of a point on the curve.
    fn param(&self, p: Point3) -> f64 {
        match self {
            Self::Ellipse {
                frame,
                major,
                minor,
            } => {
                let [x, y, _] = frame.coordinates(p);
                (y / minor).atan2(x / major)
            }
            Self::Hyperbola { frame, minor, .. } => {
                let [_, y, _] = frame.coordinates(p);
                (y / minor).asinh()
            }
            Self::Parabola { frame, .. } => frame.coordinates(p)[1],
            Self::Circle { frame, .. } => {
                let [x, y, _] = frame.coordinates(p);
                y.atan2(x)
            }
            Self::Rulings => 0.0,
        }
    }

    /// The arc from `t0` over `sweep`.
    fn arc(&self, t0: f64, sweep: f64) -> Curve3 {
        match *self {
            Self::Ellipse {
                frame,
                major,
                minor,
            } => Curve3::EllipseArc {
                frame,
                major,
                minor,
                start_angle: t0,
                sweep_angle: sweep,
            },
            Self::Hyperbola {
                frame,
                major,
                minor,
            } => Curve3::HyperbolaArc {
                frame,
                major,
                minor,
                start: t0,
                sweep,
            },
            Self::Parabola { frame, focal } => Curve3::ParabolaArc {
                frame,
                focal,
                start: t0,
                sweep,
            },
            Self::Circle { frame, radius } => Curve3::CircularArc {
                frame,
                radius,
                start_angle: t0,
                sweep_angle: sweep,
            },
            Self::Rulings => unreachable!("rulings are segments between their crossings"),
        }
    }
}

/// A point where the plane meets a rim: its end, its angle and place.
#[derive(Clone, Copy)]
struct Crossing {
    end: usize,
    angle: f64,
    point: Point3,
    t: f64,
}

/// An arc of the section inside the solid: its crossings (none for a
/// ring), its curve, and the side on its left over the wall.
struct SectionArc {
    from: Option<usize>,
    to: Option<usize>,
    curve: Curve3,
    left: Side,
}

/// A part of a rim: its crossings (none for a whole ring), its angles, its
/// side.
struct RimPart {
    end: usize,
    from: Option<usize>,
    to: Option<usize>,
    start: f64,
    sweep: f64,
    side: Side,
}

/// Everything both pieces are cut from.
struct Setup {
    frame: Frame3,
    tolerance: Tolerance,
    ends: [EndSpec; 2],
    rims: [Rim; 2],
    /// An apex's or pole's side at each end.
    points: [Option<Side>; 2],
    crossings: Vec<Crossing>,
    arcs: Vec<SectionArc>,
    rim_parts: Vec<RimPart>,
    /// Each end's whole side when the plane misses or touches it.
    disc_side: [Option<Side>; 2],
    wall: Surface,
    cut: Frame3,
}

fn mid(x: &I) -> f64 {
    let (lo, hi) = x.bounds_f64();
    0.5 * lo + 0.5 * hi
}

fn side_of(f: f64) -> Side {
    if f < 0.0 {
        Side::Below
    } else {
        Side::Above
    }
}

/// Both pieces, below first; `plane` in the solid's frame, not containing
/// its axis nor normal to it.
pub(super) fn pieces(
    primitive: &Primitive,
    frame: Frame3,
    tolerance: Tolerance,
    plane: &[R; 4],
) -> Result<Vec<Built>> {
    let setup = setup(primitive, frame, tolerance, plane)?;
    [Side::Below, Side::Above]
        .into_iter()
        .map(|side| piece(&setup, primitive, side))
        .collect()
}

#[allow(clippy::too_many_lines)]
fn setup(
    primitive: &Primitive,
    frame: Frame3,
    tolerance: Tolerance,
    plane: &[R; 4],
) -> Result<Setup> {
    let tol = tolerance.linear();
    let [a, b, c, d] = plane;
    let ab2 = a * a + b * b;
    let m2 = &ab2 + c * c;
    let ends = ends(primitive);
    // Exact: how the plane meets each end.
    let mut rims = [Rim::Point; 2];
    for (e, end) in ends.iter().enumerate() {
        let k = c * q(end.w) + d;
        rims[e] = if end.radius == 0.0 {
            if k == zero() {
                return Err(Error::OutOfDomain(
                    "a plane through an apex or pole off the axis",
                ));
            }
            Rim::Point
        } else {
            match (&k * &k).cmp(&(q(end.radius) * q(end.radius) * &ab2)) {
                Ordering::Greater => Rim::Miss,
                Ordering::Equal => Rim::Touch,
                Ordering::Less => Rim::Cross,
            }
        };
    }
    let closed = !rims.contains(&Rim::Cross);
    // A closed section winds round the axis on a cone always, on a sphere
    // when the axis crosses the plane inside it.
    let winding = match *primitive {
        Primitive::Zone { radius, .. } => *c != zero() && d * d < q(radius) * q(radius) * c * c,
        _ => true,
    };
    let pinch = matches!(rims, [Rim::Cross, Rim::Touch] | [Rim::Touch, Rim::Cross])
        || (closed && !winding && rims.contains(&Rim::Touch));
    if pinch {
        return Err(Error::Degenerate(
            "a plane tangent to a rim, pinching the wall",
        ));
    }
    // The unit normal and offset in the frame.
    let m = I::exact(m2.clone()).sqrt();
    let unit = |x: &R| -> Result<f64> {
        Ok(mid(&I::exact(x.clone())
            .div(&m)
            .ok_or(Error::Degenerate("a plane's normal"))?))
    };
    let (na, nb, nc, nd) = (unit(a)?, unit(b)?, unit(c)?, unit(d)?);
    let f_at = |p: [f64; 3]| na * p[0] + nb * p[1] + nc * p[2] + nd;
    let world = |p: [f64; 3]| frame.point(Point2::new(p[0], p[1]), p[2]);
    let wvec = |v: [f64; 3]| frame.x() * v[0] + frame.y() * v[1] + frame.normal() * v[2];
    let n_world = wvec([na, nb, nc]);
    let ab = na.hypot(nb);
    // The pieces' thickness: F's extremes over the end circles and, on a
    // sphere, its extreme points between the ends.
    let (mut fmax, mut fmin) = (f64::NEG_INFINITY, f64::INFINITY);
    for end in &ends {
        let k = nc * end.w + nd;
        fmax = fmax.max(k + end.radius * ab);
        fmin = fmin.min(k - end.radius * ab);
        if end.radius == 0.0 && k.abs() <= tol {
            return Err(Error::Degenerate(
                "a plane within the resolution of an apex or pole",
            ));
        }
    }
    if let Primitive::Zone { radius, .. } = *primitive {
        for t in [1.0, -1.0] {
            let w = t * radius * nc;
            if ends[0].w <= w && w <= ends[1].w {
                fmax = fmax.max(nd + t * radius);
                fmin = fmin.min(nd + t * radius);
            }
        }
    }
    if fmax <= tol || fmin >= -tol {
        return Err(Error::Degenerate("a piece thinner than the resolution"));
    }
    // The rims' crossings, by angle about the axis.
    let phi = nb.atan2(na);
    let mut crossings: Vec<Crossing> = Vec::new();
    let mut per_end: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
    for (e, end) in ends.iter().enumerate() {
        if !matches!(rims[e], Rim::Touch | Rim::Cross) {
            continue;
        }
        let k = c * q(end.w) + d;
        let deltas = if rims[e] == Rim::Touch {
            // F = k + r |ab| cos(angle - phi) vanishes where the cosine is
            // -sign(k).
            vec![if k > zero() { PI } else { 0.0 }]
        } else {
            let ratio = I::exact(-k)
                .div(&I::exact(q(end.radius) * q(end.radius) * &ab2).sqrt())
                .ok_or(Error::Degenerate("a rim's radius"))?;
            let delta = mid(&ratio).clamp(-1.0, 1.0).acos();
            // The smaller segment's sagitta: the rim's reach past the chord.
            let sagitta = end.radius * (1.0 - delta.cos().abs());
            if 2.0 * end.radius * delta.sin() <= tol || sagitta <= tol {
                return Err(Error::Degenerate(
                    "a plane within the resolution of a rim's tangent",
                ));
            }
            vec![-delta, delta]
        };
        for delta in deltas {
            let angle = phi + delta;
            let (s, co) = angle.sin_cos();
            per_end[e].push(crossings.len());
            crossings.push(Crossing {
                end: e,
                angle,
                point: world([end.radius * co, end.radius * s, end.w]),
                t: 0.0,
            });
        }
    }
    // The section's curve.
    let axis_x = {
        // The axis's direction projected on the plane: the conic's axis.
        let v = [-nc * na, -nc * nb, 1.0 - nc * nc];
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / l, v[1] / l, v[2] / l]
    };
    let e1_world = wvec(axis_x);
    let (section, wall, origin) = match *primitive {
        Primitive::Cone {
            bottom,
            top,
            height,
        } => {
            let wall = Surface::Cone {
                frame,
                radius: bottom,
                half_angle: (top - bottom).atan2(height),
            };
            // The apex (virtual for a frustum), exactly and rounded.
            let za = -q(bottom) * q(height) / (q(top) - q(bottom));
            let at_apex = c * &za + d;
            let za = super::rational_f64(&za);
            let h = nc * za + nd;
            let foot = [-h * na, -h * nb, za - h * nc];
            // Through the virtual apex, or within the resolution of it (its
            // height rounded): the rulings between the rims' crossings.
            if at_apex == zero() || h.abs() <= tol {
                (Section::Rulings, wall, foot)
            } else {
                // The nappe holding the solid opens along k.
                let k = if top > bottom { 1.0 } else { -1.0 };
                let (kn, k1) = (k * nc, (1.0 - nc * nc).sqrt());
                let e1 = [k * axis_x[0], k * axis_x[1], k * axis_x[2]];
                let cos2 = height * height / (height * height + (top - bottom).powi(2));
                // In the plane's coordinates (x along e1 from the apex's
                // foot, y across): A x^2 + C y^2 + D x + F0 = 0.
                let (aa, cc) = (cos2 - k1 * k1, cos2);
                let dd = 2.0 * h * k1 * kn;
                let f0 = h * h * (cos2 - kn * kn);
                let qh = q(height);
                let kind = (&qh * &qh * &m2)
                    .cmp(&(&ab2 * (&qh * &qh + (q(top) - q(bottom)) * (q(top) - q(bottom)))));
                let at = |x: f64| {
                    [
                        foot[0] + x * e1[0],
                        foot[1] + x * e1[1],
                        foot[2] + x * e1[2],
                    ]
                };
                let size = bottom.max(top).max(height);
                let section = match kind {
                    Ordering::Greater => {
                        let xc = -dd / (2.0 * aa);
                        let g = dd * dd / (4.0 * aa) - f0;
                        let (major, minor) = ((g / aa).sqrt(), (g / cc).sqrt());
                        if !(major.is_finite() && minor > 16.0 * tol) || major > 1e6 * size {
                            return Err(Error::Degenerate(
                                "a cone section within binary64 of a parabola or a point",
                            ));
                        }
                        Section::Ellipse {
                            frame: Frame3::new(world(at(xc)), n_world, wvec(e1), tolerance)?,
                            major,
                            minor,
                        }
                    }
                    Ordering::Less => {
                        let xc = -dd / (2.0 * aa);
                        let g = dd * dd / (4.0 * aa) - f0;
                        let (major, minor) = ((g / aa).sqrt(), (-g / cc).sqrt());
                        if !(major > 16.0 * tol && minor > 16.0 * tol) || major > 1e6 * size {
                            return Err(Error::Degenerate(
                                "a cone section within binary64 of a parabola or its asymptotes",
                            ));
                        }
                        // The branch on the solid's nappe: (p - apex) . k > 0.
                        let sigma = if (xc + major) * k1 - h * kn > 0.0 {
                            1.0
                        } else {
                            -1.0
                        };
                        Section::Hyperbola {
                            frame: Frame3::new(
                                world(at(xc)),
                                n_world,
                                wvec(e1) * sigma,
                                tolerance,
                            )?,
                            major,
                            minor,
                        }
                    }
                    Ordering::Equal => {
                        // C y^2 + D x + F0 = 0: the vertex at -F0 / D,
                        // opening along -sign(D) e1.
                        let sigma = -dd.signum();
                        if dd.abs() / (4.0 * cc) <= 16.0 * tol {
                            return Err(Error::Degenerate(
                                "a plane within the resolution of a frustum's virtual apex",
                            ));
                        }
                        Section::Parabola {
                            frame: Frame3::new(
                                world(at(-f0 / dd)),
                                n_world,
                                wvec(e1) * sigma,
                                tolerance,
                            )?,
                            focal: dd.abs() / (4.0 * cc),
                        }
                    }
                };
                (section, wall, foot)
            }
        }
        Primitive::Zone { radius, .. } => {
            let h = nd;
            let rho2 = I::exact(q(radius) * q(radius) - d * d / &m2);
            let rho = mid(&rho2.sqrt());
            let foot = [-h * na, -h * nb, -h * nc];
            (
                Section::Circle {
                    frame: Frame3::new(world(foot), n_world, e1_world, tolerance)?,
                    radius: rho,
                },
                Surface::Sphere { frame, radius },
                foot,
            )
        }
        Primitive::Torus { .. } => unreachable!("a torus's sections are S8d.3's"),
    };
    // The cut face's x axis is an ellipse's own (its plane pcurve needs it).
    let cut_x = match &section {
        Section::Ellipse { frame, .. } => frame.x(),
        _ => e1_world,
    };
    let cut = Frame3::new(world(origin), n_world, cut_x, tolerance)?;
    for x in &mut crossings {
        x.t = section.param(x.point);
    }
    // The section's arcs inside the solid.
    let height_of = |p: Point3| frame.coordinates(p)[2];
    // How far inside the ends a point lies (negative outside).
    let margin = |p: Point3| {
        let w = height_of(p);
        (w - ends[0].w).min(ends[1].w - w)
    };
    let mut arcs: Vec<SectionArc> = Vec::new();
    let mut spans: Vec<(Option<usize>, Option<usize>, Curve3)> = Vec::new();
    if let Section::Rulings = section {
        // Each ruling joins the bottom crossing to the top one at the same
        // angle about the axis.
        for &i in &per_end[0] {
            let j = *per_end[1]
                .iter()
                .min_by(|&&x, &&y| {
                    let gap = |z: usize| {
                        let g = (crossings[z].angle - crossings[i].angle).rem_euclid(TAU);
                        g.min(TAU - g)
                    };
                    gap(x).total_cmp(&gap(y))
                })
                .ok_or(Error::Degenerate(
                    "a plane within the resolution of a frustum's virtual apex",
                ))?;
            spans.push((
                Some(i),
                Some(j),
                Curve3::LineSegment {
                    start: crossings[i].point,
                    end: crossings[j].point,
                },
            ));
        }
    } else {
        let mut order: Vec<usize> = (0..crossings.len()).collect();
        let key = |i: usize| {
            if section.periodic() {
                crossings[i].t.rem_euclid(TAU)
            } else {
                crossings[i].t
            }
        };
        order.sort_by(|&x, &y| key(x).total_cmp(&key(y)));
        if order.is_empty() {
            if !section.periodic() {
                // An open conic always leaves the solid: missing both rims
                // happens only within rounding of a tangency.
                return Err(Error::Degenerate(
                    "a plane within the resolution of a rim's tangent",
                ));
            }
            spans.push((None, None, section.arc(0.0, TAU)));
        } else {
            let n = order.len();
            let pairs = if section.periodic() { n } else { n - 1 };
            let mut candidates = Vec::new();
            for s in 0..pairs {
                let (i, j) = (order[s], order[(s + 1) % n]);
                let (t0, mut t1) = (key(i), key(j));
                if section.periodic() && t1 <= t0 {
                    t1 += TAU;
                }
                let curve = section.arc(t0, t1 - t0);
                candidates.push((i, j, margin(curve.point(0.5)), curve));
            }
            // Inside and outside alternate at each crossing (not at a touch):
            // anchored at the arc farthest from the ends (an open conic's
            // rays lie outside), so an arc within rounding of an end never
            // decides by its own midpoint.
            let flips = |i: usize| rims[crossings[i].end] == Rim::Cross;
            let mut inside_at: Vec<bool> = vec![false; candidates.len()];
            let anchor = if section.periodic() {
                (0..candidates.len())
                    .max_by(|&x, &y| candidates[x].2.abs().total_cmp(&candidates[y].2.abs()))
                    .expect("an arc")
            } else {
                usize::MAX
            };
            if anchor == usize::MAX {
                let mut status = false;
                for (k, c) in candidates.iter().enumerate() {
                    if flips(c.0) {
                        status = !status;
                    }
                    inside_at[k] = status;
                }
            } else {
                let m = candidates.len();
                inside_at[anchor] = candidates[anchor].2 > 0.0;
                for step in 1..m {
                    let k = (anchor + step) % m;
                    let before = (anchor + step - 1) % m;
                    inside_at[k] = if flips(candidates[k].0) {
                        !inside_at[before]
                    } else {
                        inside_at[before]
                    };
                }
            }
            for ((i, j, _, curve), inside) in candidates.into_iter().zip(inside_at) {
                if inside {
                    spans.push((Some(i), Some(j), curve));
                }
            }
        }
    }
    // The side on each arc's left over the wall (about its outward normal).
    let outward = |p: Point3| -> Vec3 {
        let [x, y, z] = frame.coordinates(p);
        match *primitive {
            Primitive::Cone {
                bottom,
                top,
                height,
            } => {
                let r = x.hypot(y);
                wvec([x / r, y / r, -(top - bottom) / height])
            }
            _ => wvec([x, y, z]),
        }
    };
    for (from, to, curve) in spans {
        let (p0, p1) = (curve.point(0.5 - 1e-6), curve.point(0.5 + 1e-6));
        let tangent = p1 - p0;
        let left = outward(curve.point(0.5)).cross(tangent);
        arcs.push(SectionArc {
            from,
            to,
            curve,
            left: if left.dot(n_world) > 0.0 {
                Side::Above
            } else {
                Side::Below
            },
        });
    }
    // The exact decisions on the stored rims and the rounded section
    // disagree only within rounding of a tangency.
    if arcs.is_empty() {
        return Err(Error::Degenerate(
            "a plane within the resolution of a rim's tangent",
        ));
    }
    // The rims' parts and the discs' sides.
    let mut rim_parts = Vec::new();
    let mut disc_side = [None, None];
    let mut points = [None, None];
    for (e, end) in ends.iter().enumerate() {
        let k = nc * end.w + nd;
        match rims[e] {
            Rim::Point => points[e] = Some(side_of(k)),
            Rim::Miss => {
                disc_side[e] = Some(side_of(k));
                rim_parts.push(RimPart {
                    end: e,
                    from: None,
                    to: None,
                    start: 0.0,
                    sweep: TAU,
                    side: side_of(k),
                });
            }
            Rim::Touch => {
                let i = per_end[e][0];
                disc_side[e] = Some(side_of(k));
                rim_parts.push(RimPart {
                    end: e,
                    from: Some(i),
                    to: Some(i),
                    start: crossings[i].angle,
                    sweep: TAU,
                    side: side_of(k),
                });
            }
            Rim::Cross => {
                // The two arcs lie on opposite sides: the one whose middle is
                // farther from the plane decides.
                let (i, j) = (per_end[e][0], per_end[e][1]);
                let arcs: Vec<(usize, usize, f64, f64, f64)> = [(i, j), (j, i)]
                    .into_iter()
                    .map(|(x, y)| {
                        let start = crossings[x].angle;
                        let sweep = (crossings[y].angle - start).rem_euclid(TAU);
                        let (s, co) = (start + 0.5 * sweep).sin_cos();
                        (
                            x,
                            y,
                            start,
                            sweep,
                            f_at([end.radius * co, end.radius * s, end.w]),
                        )
                    })
                    .collect();
                let first = side_of(if arcs[0].4.abs() >= arcs[1].4.abs() {
                    arcs[0].4
                } else {
                    -arcs[1].4
                });
                for (k, &(x, y, start, sweep, _)) in arcs.iter().enumerate() {
                    let side = match (k, first) {
                        (0, side) => side,
                        (_, Side::Below) => Side::Above,
                        (_, Side::Above) => Side::Below,
                    };
                    rim_parts.push(RimPart {
                        end: e,
                        from: Some(x),
                        to: Some(y),
                        start,
                        sweep,
                        side,
                    });
                }
            }
        }
    }
    Ok(Setup {
        frame,
        tolerance,
        ends,
        rims,
        points,
        crossings,
        arcs,
        rim_parts,
        disc_side,
        wall,
        cut,
    })
}

fn sense(forward: bool) -> Orientation {
    if forward {
        Orientation::Forward
    } else {
        Orientation::Reversed
    }
}

/// Chain edge uses into loops: from each unused use, on along an unused
/// use leaving where the last ended until none is left there (a closed
/// section touching a rim chains through its vertex into one loop).
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
            return Err(Error::InvalidTopology("an open loop in a split piece"));
        }
        loops.push(list);
    }
    Ok(loops)
}

/// A face's loop as built: fins with a winding, or a vertex.
enum LoopSpec {
    Edges(Vec<Fin>, [i32; 2]),
    Vertex(VertexId),
}

#[allow(clippy::too_many_lines)]
fn piece(setup: &Setup, primitive: &Primitive, side: Side) -> Result<Built> {
    let k = u32::from(side == Side::Above);
    let frame = setup.frame;
    let (n, x) = (frame.normal(), frame.x());
    let point_role = match primitive {
        Primitive::Cone { .. } => Role::Apex,
        _ => Role::Pole,
    };
    let ring_role = |e: usize| {
        if e == 0 {
            Role::BottomEdge
        } else {
            Role::TopEdge
        }
    };
    let cap_role = |e: usize| {
        if e == 0 {
            Role::StartCap
        } else {
            Role::EndCap
        }
    };
    let mut parts = TopologyParts::default();
    let mut plans: Vec<(Slot, Plan)> = Vec::new();
    // Vertices on demand: crossings, and an apex or pole on this side.
    let mut vertex_of: BTreeMap<usize, VertexId> = BTreeMap::new();
    let mut vertex = |i: usize, parts: &mut TopologyParts, plans: &mut Vec<(Slot, Plan)>| {
        *vertex_of.entry(i).or_insert_with(|| {
            let id = VertexId(parts.vertices.len());
            parts.vertices.push(Vertex {
                position: setup.crossings[i].point,
                enclosure: None,
            });
            let e = setup.crossings[i].end;
            plans.push((
                Slot::Vertex(id),
                Plan::New(vec![(ring_role(e), e)], Role::CutVertex, 4 * k + i as u32),
            ));
            id
        })
    };
    let add_edge = |start: Option<VertexId>,
                    end: Option<VertexId>,
                    curve: Curve3,
                    plan: Plan,
                    parts: &mut TopologyParts,
                    plans: &mut Vec<(Slot, Plan)>| {
        let id = EdgeId(parts.edges.len());
        parts.edges.push(Edge {
            start,
            end,
            curve,
            fins: Vec::new(),
        });
        plans.push((Slot::Edge(id), plan));
        id
    };
    // The rim parts on this side, and each crossed rim's chord.
    let mut rim_edges: Vec<(usize, EdgeId)> = Vec::new();
    for part in setup.rim_parts.iter().filter(|p| p.side == side) {
        let end = setup.ends[part.end];
        let ring_frame = Frame3::new(frame.point(Point2::default(), end.w), n, x, setup.tolerance)?;
        let (start, stop, curve) = match (part.from, part.to) {
            (Some(i), Some(j)) => (
                Some(vertex(i, &mut parts, &mut plans)),
                Some(vertex(j, &mut parts, &mut plans)),
                Curve3::CircularArc {
                    frame: ring_frame,
                    radius: end.radius,
                    start_angle: part.start,
                    sweep_angle: part.sweep,
                },
            ),
            _ => (
                None,
                None,
                Curve3::Circle {
                    frame: ring_frame,
                    radius: end.radius,
                },
            ),
        };
        let plan = if setup.rims[part.end] == Rim::Cross {
            Plan::Child(ring_role(part.end), part.end)
        } else {
            Plan::Kept(ring_role(part.end), part.end)
        };
        let id = add_edge(start, stop, curve, plan, &mut parts, &mut plans);
        rim_edges.push((part.end, id));
    }
    let mut chords: [Option<EdgeId>; 2] = [None, None];
    for (e, chord) in chords.iter_mut().enumerate() {
        if setup.rims[e] != Rim::Cross || !rim_edges.iter().any(|r| r.0 == e) {
            continue;
        }
        let list: Vec<usize> = (0..setup.crossings.len())
            .filter(|&i| setup.crossings[i].end == e)
            .collect();
        let (i, j) = (list[0], list[1]);
        let (vi, vj) = (
            vertex(i, &mut parts, &mut plans),
            vertex(j, &mut parts, &mut plans),
        );
        *chord = Some(add_edge(
            Some(vi),
            Some(vj),
            Curve3::LineSegment {
                start: setup.crossings[i].point,
                end: setup.crossings[j].point,
            },
            Plan::New(vec![(cap_role(e), e)], Role::CutEdge, k),
            &mut parts,
            &mut plans,
        ));
    }
    // The section's arcs.
    let mut section_edges: Vec<(EdgeId, Side)> = Vec::new();
    for (j, arc) in setup.arcs.iter().enumerate() {
        let start = arc.from.map(|i| vertex(i, &mut parts, &mut plans));
        let stop = arc.to.map(|i| vertex(i, &mut parts, &mut plans));
        let id = add_edge(
            start,
            stop,
            arc.curve.clone(),
            Plan::New(vec![(Role::Wall, 0)], Role::CutEdge, 4 * k + j as u32 + 2),
            &mut parts,
            &mut plans,
        );
        section_edges.push((id, arc.left));
    }
    // An apex or pole on this side.
    let mut point_vertex = None;
    for e in 0..2 {
        if setup.points[e] == Some(side) {
            let id = VertexId(parts.vertices.len());
            parts.vertices.push(Vertex {
                position: frame.point(Point2::default(), setup.ends[e].w),
                enclosure: None,
            });
            plans.push((Slot::Vertex(id), Plan::Kept(point_role, e)));
            point_vertex = Some(id);
        }
    }
    // Uses: the wall runs a bottom rim in +u and a top one in -u, a section
    // arc forward when this side lies on its left; a disc and the cut face
    // use each edge they share with the wall or a disc the other way.
    let mut wall_uses: Vec<(EdgeId, bool)> =
        rim_edges.iter().map(|&(e, id)| (id, e == 0)).collect();
    wall_uses.extend(section_edges.iter().map(|&(id, left)| (id, left == side)));
    let mut disc_uses: [Vec<(EdgeId, bool)>; 2] = [Vec::new(), Vec::new()];
    for &(e, id) in &rim_edges {
        disc_uses[e].push((id, e == 1));
        if let Some(chord) = chords[e] {
            // The chord closes the arc: from its end back to its start.
            let arc = &parts.edges[id.0];
            let (arc_start, arc_end) = if e == 1 {
                (arc.start, arc.end)
            } else {
                (arc.end, arc.start)
            };
            let forward = parts.edges[chord.0].start == arc_end;
            debug_assert_eq!(
                if forward {
                    parts.edges[chord.0].end
                } else {
                    parts.edges[chord.0].start
                },
                arc_start
            );
            disc_uses[e].push((chord, forward));
        }
    }
    let mut cut_uses: Vec<(EdgeId, bool)> = section_edges
        .iter()
        .map(|&(id, left)| (id, left != side))
        .collect();
    for list in &disc_uses {
        for &(id, forward) in list {
            if chords.contains(&Some(id)) {
                cut_uses.push((id, !forward));
            }
        }
    }
    // Faces.
    let add_face = |surface: Surface,
                    face_sense: Orientation,
                    loops: Vec<LoopSpec>,
                    plan: Plan,
                    parts: &mut TopologyParts,
                    plans: &mut Vec<(Slot, Plan)>| {
        let mut ids = Vec::new();
        for spec in loops {
            match spec {
                LoopSpec::Edges(fins, winding) => {
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
                }
                LoopSpec::Vertex(v) => parts.loops.push(Loop::Vertex(v)),
            }
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
    let planar = |surface: Frame3, uses: &[(EdgeId, bool)], parts: &TopologyParts| {
        Ok::<_, Error>(
            chain(uses, parts)?
                .into_iter()
                .map(|list| {
                    let fins = list
                        .iter()
                        .map(|&(edge, forward)| Fin {
                            edge,
                            sense: sense(forward),
                            pcurve: plane_pcurve(
                                &parts.edges[edge.0].curve,
                                sense(forward),
                                surface,
                            ),
                            enclosure: None,
                        })
                        .collect();
                    LoopSpec::Edges(fins, [0, 0])
                })
                .collect::<Vec<_>>(),
        )
    };
    // The wall: each loop walked on the cover from its first vertex.
    let mut wall_loops = Vec::new();
    for list in chain(&wall_uses, &parts)? {
        let mut fins = Vec::new();
        let mut first: Option<f64> = None;
        let mut at: Option<Point2> = None;
        for &(edge, forward) in &list {
            let curve = parts.edges[edge.0].curve.clone();
            let rim = rim_edges.iter().find(|r| r.1 == edge).map(|r| r.0);
            let pcurve = if let Some(e) = rim {
                let v = setup.ends[e].v;
                let (start, sweep) = match curve {
                    Curve3::CircularArc {
                        start_angle,
                        sweep_angle,
                        ..
                    } => (start_angle, sweep_angle),
                    _ => (0.0, TAU),
                };
                let (u0, u1) = if forward {
                    (start, start + sweep)
                } else {
                    (start + sweep, start)
                };
                let shift = at.map_or(0.0, |p| TAU * ((p.x - u0) / TAU).round());
                first.get_or_insert(u0 + shift);
                at = Some(Point2::new(u1 + shift, v));
                Curve2::LineSegment {
                    start: Point2::new(u0 + shift, v),
                    end: Point2::new(u1 + shift, v),
                }
            } else {
                let p = curve.point(if forward { 0.0 } else { 1.0 });
                let mut uv = Projection::inverse(&setup.wall, p)
                    .ok_or(Error::InvalidTopology("a wall without an inverse"))?;
                if let Some(prev) = at {
                    uv.x += TAU * ((prev.x - uv.x) / TAU).round();
                }
                let projection = Projection::new(curve, setup.wall.clone(), !forward, uv, ANCHORS)
                    .ok_or(Error::Degenerate(
                        "a section turning too fast about the axis",
                    ))?;
                first.get_or_insert(uv.x);
                at = Some(projection.point(1.0));
                Curve2::Projection(Box::new(projection))
            };
            fins.push(Fin {
                edge,
                sense: sense(forward),
                pcurve,
                enclosure: None,
            });
        }
        let turns = match (first, at) {
            (Some(u0), Some(end)) => ((end.x - u0) / TAU).round() as i32,
            _ => 0,
        };
        wall_loops.push(LoopSpec::Edges(fins, [turns, 0]));
    }
    if let Some(v) = point_vertex {
        wall_loops.push(LoopSpec::Vertex(v));
    }
    add_face(
        setup.wall.clone(),
        Orientation::Forward,
        wall_loops,
        Plan::Child(Role::Wall, 0),
        &mut parts,
        &mut plans,
    );
    // The discs (or their parts) on this side: outward along -n at the
    // bottom, n at the top.
    for (e, uses) in disc_uses.iter().enumerate() {
        if uses.is_empty() {
            continue;
        }
        let normal = if e == 1 { n } else { n * -1.0 };
        let disc = Frame3::new(
            frame.point(Point2::default(), setup.ends[e].w),
            normal,
            x,
            setup.tolerance,
        )?;
        let loops = planar(disc, uses, &parts)?;
        let plan = if setup.disc_side[e].is_some() {
            Plan::Kept(cap_role(e), e)
        } else {
            Plan::Child(cap_role(e), e)
        };
        add_face(
            Surface::Plane(disc),
            Orientation::Forward,
            loops,
            plan,
            &mut parts,
            &mut plans,
        );
    }
    // The cut face: its normal the plane's, facing out of the lower piece.
    let loops = planar(setup.cut, &cut_uses, &parts)?;
    let mut from = vec![(Role::Wall, 0)];
    for e in 0..2 {
        if setup.rims[e] == Rim::Cross {
            from.push((cap_role(e), e));
        }
    }
    add_face(
        Surface::Plane(setup.cut),
        if side == Side::Below {
            Orientation::Forward
        } else {
            Orientation::Reversed
        },
        loops,
        Plan::New(from, Role::CutFace, k),
        &mut parts,
        &mut plans,
    );
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
    plans.push((Slot::Region(RegionId(1)), Plan::Child(Role::Region, 0)));
    Ok(Built { side, parts, plans })
}
