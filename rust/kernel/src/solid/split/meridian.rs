//! S8c.2: a cone, frustum or sphere zone split by a plane containing its
//! axis.
//!
//! The plane `a u + b v + c w + d` contains the axis when `c = d = 0` (or
//! lies within a quarter of the resolution of it over the solid's height: a
//! frame whose stored axes are not exactly orthogonal). Each piece is a half
//! of the solid over the angles on its side, `u` from where the plane's
//! trace `t = (-b, a) / |(a, b)|` leaves the axis to its opposite: the half
//! wall on the input's own surface, bounded by the end circles' halves (arcs
//! of the input's rings) and two meridians (a cone's rulings, a sphere's
//! great-circle arcs, both `u`-constant pcurves); the halves of the end
//! discs, closed by chords along the trace; and the cut face in the plane,
//! bounded by the meridians and the chords. An apex or pole ends both
//! meridians. Each piece validates as a general body.
//!
//! Names follow provenance: the wall, the region, the discs and the rings
//! are `Split` into one child per piece (below first), an apex or pole into
//! a copy per piece; the chords' ends are `Generated` from the ring they
//! cut, the chords from their disc, the meridians from the wall, the cut
//! face from the wall and the discs.
use super::{q, rational_f64, zero, Side};
use crate::certified::{Interval as I, Real};
use crate::history::{History, Relation};
use crate::identity::{
    Derivation, EntityId, EntityKind, OperationId, OperationKind, Parent, ProfileElement, Role,
};
use crate::solid::{Construction, Context, Solid};
use crate::topology::{
    plane_pcurve, Curve2, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId,
    Orientation, Region, RegionId, RegionKind, Shell, ShellId, Side as FaceSide, Slot, Surface,
    Topology, TopologyParts, Vertex, VertexId,
};
use crate::{Error, Frame3, Point2, Point3, Result, Tolerance, Vec3};
use num_rational::BigRational as R;
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::{FRAC_PI_2, PI};

/// The pieces with their sides, and the history.
type Pieces = (Vec<(Side, Solid)>, History);

/// A revolved solid's ends and its wall, as the S3 builders make them.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Primitive {
    Cone {
        bottom: f64,
        top: f64,
        height: f64,
    },
    Zone {
        radius: f64,
        low: f64,
        high: f64,
    },
    /// A whole torus (S8d.1's bands).
    Torus {
        major: f64,
        minor: f64,
    },
}

/// How a half was made, kept so a rigid motion rebuilds it exactly and it
/// classifies points: the primitive, the plane in its frame, the side.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Half {
    primitive: Primitive,
    tolerance: Tolerance,
    plane: [R; 4],
    /// The piece's index in the split's order.
    index: usize,
}

/// An end of the primitive: its height, its circle's radius (zero at an
/// apex or pole) and, on a sphere, its latitude.
#[derive(Clone, Copy)]
struct EndSpec {
    w: f64,
    radius: f64,
    /// The wall's `v` there: a cone's distance along its generatrix, a
    /// sphere's latitude.
    v: f64,
}

/// A piece before naming: its side, parts and every slot's provenance.
struct Built {
    side: Side,
    parts: TopologyParts,
    plans: Vec<(Slot, Plan)>,
}

/// A slot's provenance: an input entity by role (and end: 0 bottom, 1 top)
/// split into this piece's child, or new from input entities.
#[derive(Clone)]
enum Plan {
    Child(Role, usize),
    New(Vec<(Role, usize)>, Role, u32),
}

fn ends(primitive: &Primitive) -> [EndSpec; 2] {
    match *primitive {
        Primitive::Cone {
            bottom,
            top,
            height,
        } => [
            EndSpec {
                w: 0.0,
                radius: bottom,
                v: 0.0,
            },
            EndSpec {
                w: height,
                radius: top,
                v: height.hypot(top - bottom),
            },
        ],
        Primitive::Zone { radius, low, high } => {
            let at = |latitude: f64, pole: f64| EndSpec {
                w: if latitude == pole {
                    pole.signum() * radius
                } else {
                    radius * latitude.sin()
                },
                radius: if latitude == pole {
                    0.0
                } else {
                    radius * latitude.cos()
                },
                v: latitude,
            };
            [at(low, -FRAC_PI_2), at(high, FRAC_PI_2)]
        }
        Primitive::Torus { major, minor } => [
            EndSpec {
                w: -minor,
                radius: major,
                v: 0.0,
            },
            EndSpec {
                w: minor,
                radius: major,
                v: 0.0,
            },
        ],
    }
}

/// Whether the plane contains the axis, to within a quarter of the
/// resolution over the solid's height.
pub(super) fn contains_axis(plane: &[R; 4], w0: f64, w1: f64, tolerance: Tolerance) -> bool {
    let [a, b, c, d] = plane;
    let m2 = a * a + b * b + c * c;
    let tol = q(tolerance.linear());
    [w0, w1].iter().all(|w| {
        let f = c * q(*w) + d;
        &f * &f * R::from_integer(16.into()) <= &tol * &tol * &m2
    })
}

/// Both halves of the primitive on `frame`, below first.
fn halves(
    primitive: &Primitive,
    frame: Frame3,
    tolerance: Tolerance,
    plane: &[R; 4],
    world_normal: Vec3,
) -> Result<Vec<Built>> {
    let [a, b, _, _] = plane;
    let ab = I::exact(a * a + b * b).sqrt();
    // The trace's direction in the frame: t = (-b, a) / |(a, b)|.
    let unit = |x: &R| {
        let (lo, hi) = I::exact(x.clone())
            .div(&ab)
            .ok_or(Error::Degenerate("a plane normal to the axis"))?
            .bounds_f64();
        Ok::<f64, Error>(0.5 * lo + 0.5 * hi)
    };
    let t = [unit(&-b.clone())?, unit(a)?];
    let angle = t[1].atan2(t[0]);
    // The lower side (F < 0) turns counter-clockwise from t to -t: its
    // angles are `angle` to `angle + pi`; the upper's from -t back to t.
    let mut out = Vec::new();
    for (side, u0) in [(Side::Below, angle), (Side::Above, angle - PI)] {
        out.push(half(primitive, frame, tolerance, side, u0, world_normal)?);
    }
    Ok(out)
}

#[allow(clippy::too_many_lines)]
fn half(
    primitive: &Primitive,
    frame: Frame3,
    tolerance: Tolerance,
    side: Side,
    u0: f64,
    world_normal: Vec3,
) -> Result<Built> {
    let u1 = u0 + PI;
    let [bottom, top] = ends(primitive);
    let (n, x) = (frame.normal(), frame.x());
    let dir = |u: f64| {
        let (s, c) = u.sin_cos();
        Point2::new(c, s)
    };
    let at = |u: f64, e: &EndSpec| {
        let d = dir(u);
        frame.point(Point2::new(e.radius * d.x, e.radius * d.y), e.w)
    };
    let mut parts = TopologyParts::default();
    let mut plans: Vec<(Slot, Plan)> = Vec::new();
    let mut add_vertex = |p: Point3, plan: Plan, parts: &mut TopologyParts| {
        let id = VertexId(parts.vertices.len());
        parts.vertices.push(Vertex {
            position: p,
            enclosure: None,
        });
        plans.push((Slot::Vertex(id), plan));
        id
    };
    let k = u32::from(side == Side::Above);
    // Per end (0 bottom, 1 top): the vertices at u0 and u1 (one at an apex
    // or pole).
    let mut corners: Vec<(VertexId, VertexId)> = Vec::new();
    for (e, end) in [bottom, top].iter().enumerate() {
        if end.radius == 0.0 {
            let role = match primitive {
                Primitive::Cone { .. } => Role::Apex,
                Primitive::Zone { .. } => Role::Pole,
                Primitive::Torus { .. } => unreachable!("a torus's pieces are bands or wedges"),
            };
            let p = frame.point(Point2::default(), end.w);
            let v = add_vertex(p, Plan::Child(role, e), &mut parts);
            corners.push((v, v));
        } else {
            let ring = if e == 0 {
                Role::BottomEdge
            } else {
                Role::TopEdge
            };
            let a = add_vertex(
                at(u0, end),
                Plan::New(vec![(ring, e)], Role::CutVertex, 2 * k),
                &mut parts,
            );
            let b = add_vertex(
                at(u1, end),
                Plan::New(vec![(ring, e)], Role::CutVertex, 2 * k + 1),
                &mut parts,
            );
            corners.push((a, b));
        }
    }
    let pos = |parts: &TopologyParts, v: VertexId| parts.vertices[v.0].position;
    let add_edge = |start: VertexId,
                    end: VertexId,
                    curve: Curve3,
                    plan: Plan,
                    parts: &mut TopologyParts,
                    plans: &mut Vec<(Slot, Plan)>| {
        let id = EdgeId(parts.edges.len());
        parts.edges.push(Edge {
            start: Some(start),
            end: Some(end),
            curve,
            fins: Vec::new(),
        });
        plans.push((Slot::Edge(id), plan));
        id
    };
    // The rings' halves (u0 to u1) and the chords (u1 back to u0).
    let mut arcs: [Option<EdgeId>; 2] = [None, None];
    let mut chords: [Option<EdgeId>; 2] = [None, None];
    for (e, end) in [bottom, top].iter().enumerate() {
        if end.radius == 0.0 {
            continue;
        }
        let (a, b) = corners[e];
        let ring_frame = Frame3::new(frame.point(Point2::default(), end.w), n, x, tolerance)?;
        let ring = if e == 0 {
            Role::BottomEdge
        } else {
            Role::TopEdge
        };
        let cap = if e == 0 { Role::StartCap } else { Role::EndCap };
        arcs[e] = Some(add_edge(
            a,
            b,
            Curve3::CircularArc {
                frame: ring_frame,
                radius: end.radius,
                start_angle: u0,
                sweep_angle: PI,
            },
            Plan::Child(ring, e),
            &mut parts,
            &mut plans,
        ));
        let (pa, pb) = (pos(&parts, a), pos(&parts, b));
        chords[e] = Some(add_edge(
            b,
            a,
            Curve3::LineSegment { start: pb, end: pa },
            Plan::New(vec![(cap, e)], Role::CutEdge, k),
            &mut parts,
            &mut plans,
        ));
    }
    // The meridians at u0 and u1, from the bottom up.
    let mut meridians = [EdgeId(0); 2];
    for (j, u) in [u0, u1].into_iter().enumerate() {
        let (lo, hi) = (
            if j == 0 { corners[0].0 } else { corners[0].1 },
            if j == 0 { corners[1].0 } else { corners[1].1 },
        );
        let curve = match primitive {
            Primitive::Torus { .. } => unreachable!("a torus's pieces are bands or wedges"),
            Primitive::Cone { .. } => Curve3::LineSegment {
                start: pos(&parts, lo),
                end: pos(&parts, hi),
            },
            Primitive::Zone { radius, .. } => {
                // In the plane of the meridian: x its radial direction, y the
                // axis, so the angle is the latitude.
                let d = dir(u);
                let e = x * d.x + frame.y() * d.y;
                Curve3::CircularArc {
                    frame: Frame3::new(frame.origin(), e.cross(n), e, tolerance)?,
                    radius: *radius,
                    start_angle: bottom.v,
                    sweep_angle: top.v - bottom.v,
                }
            }
        };
        meridians[j] = add_edge(
            lo,
            hi,
            curve,
            Plan::New(vec![(Role::Wall, 0)], Role::CutEdge, 2 * k + j as u32),
            &mut parts,
            &mut plans,
        );
    }
    // Faces.
    let add_face = |surface: Surface,
                    sense: Orientation,
                    fins: Vec<Fin>,
                    plan: Plan,
                    parts: &mut TopologyParts,
                    plans: &mut Vec<(Slot, Plan)>| {
        let mut list = Vec::new();
        for fin in fins {
            parts.edges[fin.edge.0].fins.push(FinId(parts.fins.len()));
            list.push(FinId(parts.fins.len()));
            parts.fins.push(fin);
        }
        parts.loops.push(Loop::Edges {
            fins: list,
            winding: [0, 0],
        });
        parts.faces.push(Face {
            surface,
            sense,
            loops: vec![LoopId(parts.loops.len() - 1)],
            front: ShellId(0),
            back: ShellId(1),
            enclosure: None,
        });
        plans.push((Slot::Face(FaceId(parts.faces.len() - 1)), plan));
    };
    let sense = |forward: bool| {
        if forward {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    };
    let line = |a: (f64, f64), b: (f64, f64)| Curve2::LineSegment {
        start: Point2::new(a.0, a.1),
        end: Point2::new(b.0, b.1),
    };
    // The half wall: the bottom arc (u0 to u1), the meridian at u1 up, the
    // top arc back, the meridian at u0 down; an apex or pole ends the
    // meridians there.
    let wall_surface = match *primitive {
        Primitive::Cone {
            bottom: r1,
            top: r2,
            height,
        } => Surface::Cone {
            frame,
            radius: r1,
            half_angle: (r2 - r1).atan2(height),
        },
        Primitive::Zone { radius, .. } => Surface::Sphere { frame, radius },
        Primitive::Torus { .. } => unreachable!("a torus's pieces are bands or wedges"),
    };
    let (vb, vt) = (bottom.v, top.v);
    let mut fins = Vec::new();
    if let Some(e) = arcs[0] {
        fins.push(Fin {
            edge: e,
            sense: Orientation::Forward,
            pcurve: line((u0, vb), (u1, vb)),
            enclosure: None,
        });
    }
    fins.push(Fin {
        edge: meridians[1],
        sense: Orientation::Forward,
        pcurve: line((u1, vb), (u1, vt)),
        enclosure: None,
    });
    if let Some(e) = arcs[1] {
        fins.push(Fin {
            edge: e,
            sense: Orientation::Reversed,
            pcurve: line((u1, vt), (u0, vt)),
            enclosure: None,
        });
    }
    fins.push(Fin {
        edge: meridians[0],
        sense: Orientation::Reversed,
        pcurve: line((u0, vt), (u0, vb)),
        enclosure: None,
    });
    add_face(
        wall_surface,
        Orientation::Forward,
        fins,
        Plan::Child(Role::Wall, 0),
        &mut parts,
        &mut plans,
    );
    // The half discs: outward along -n at the bottom (the chord, then the
    // arc backwards), along n at the top (the arc, then the chord).
    for (e, end) in [bottom, top].iter().enumerate() {
        let (Some(arc), Some(chord)) = (arcs[e], chords[e]) else {
            continue;
        };
        let upper = e == 1;
        let normal = if upper { n } else { n * -1.0 };
        let disc = Frame3::new(frame.point(Point2::default(), end.w), normal, x, tolerance)?;
        let uses = if upper {
            [(arc, true), (chord, true)]
        } else {
            [(chord, false), (arc, false)]
        };
        let fins = uses
            .iter()
            .map(|&(edge, forward)| Fin {
                edge,
                sense: sense(forward),
                pcurve: plane_pcurve(&parts.edges[edge.0].curve, sense(forward), disc),
                enclosure: None,
            })
            .collect();
        let cap = if upper { Role::EndCap } else { Role::StartCap };
        add_face(
            Surface::Plane(disc),
            Orientation::Forward,
            fins,
            Plan::Child(cap, e),
            &mut parts,
            &mut plans,
        );
    }
    // The cut face on the plane (its normal the plane's, facing out of the
    // lower piece): the meridians and the chords, turned to run
    // counter-clockwise about the outward normal.
    let cut = Frame3::new(frame.origin(), world_normal, n, tolerance)?;
    let mut uses: Vec<(EdgeId, bool)> = vec![(meridians[0], true)];
    if let Some(chord) = chords[1] {
        // The top chord runs from u1 to u0: back towards u1 here.
        uses.push((chord, false));
    }
    uses.push((meridians[1], false));
    if let Some(chord) = chords[0] {
        uses.push((chord, true));
    }
    // Their orientation about the frame's normal, by the corners' polygon.
    let corner_list: Vec<Point3> = uses
        .iter()
        .map(|&(edge, forward)| {
            let e = &parts.edges[edge.0];
            pos(&parts, if forward { e.start } else { e.end }.expect("ends"))
        })
        .collect();
    let mut twice = 0.0;
    for i in 0..corner_list.len() {
        let [x0, y0, _] = cut.coordinates(corner_list[i]);
        let [x1, y1, _] = cut.coordinates(corner_list[(i + 1) % corner_list.len()]);
        twice += x0 * y1 - x1 * y0;
    }
    // Outward along the plane's normal for the lower piece.
    let (face_sense, want_ccw) = match side {
        Side::Below => (Orientation::Forward, true),
        Side::Above => (Orientation::Reversed, false),
    };
    if (twice > 0.0) != want_ccw {
        uses.reverse();
        for u in &mut uses {
            u.1 = !u.1;
        }
    }
    let fins = uses
        .iter()
        .map(|&(edge, forward)| Fin {
            edge,
            sense: sense(forward),
            pcurve: plane_pcurve(&parts.edges[edge.0].curve, sense(forward), cut),
            enclosure: None,
        })
        .collect();
    let mut from = vec![(Role::Wall, 0)];
    for (e, chord) in chords.iter().enumerate() {
        if chord.is_some() {
            from.push((if e == 0 { Role::StartCap } else { Role::EndCap }, e));
        }
    }
    add_face(
        Surface::Plane(cut),
        face_sense,
        fins,
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

/// The input's entity for a role at an end (an apex or pole by its place).
fn input_entity(input: &Topology, role: Role, end: usize) -> Option<EntityId> {
    input.ids().map(|(id, _)| id).find(|id| {
        let d = input.derivation(*id).expect("a derivation");
        let top = matches!(
            d.parents.first(),
            Some(Parent::Profile {
                element: ProfileElement::Vertex(2),
                ..
            })
        );
        d.role == role && (!matches!(role, Role::Apex | Role::Pole) || top == (end == 1))
    })
}

impl Half {
    pub(super) fn new(
        primitive: Primitive,
        tolerance: Tolerance,
        plane: [R; 4],
        index: usize,
    ) -> Self {
        Self {
            primitive,
            tolerance,
            plane,
            index,
        }
    }

    /// This piece's solid on `frame` with its topology.
    pub(super) fn solid(
        self,
        frame: Frame3,
        topology: Topology,
        operation: OperationId,
    ) -> Result<Solid> {
        piece_solid(self, frame, topology, operation)
    }

    pub(crate) fn tolerance(&self) -> Tolerance {
        self.tolerance
    }

    /// The same half in another frame, its ids external (the caller
    /// restores them).
    pub(crate) fn rebuilt(&self, operation: OperationId, frame: Frame3) -> Result<Solid> {
        if let Primitive::Torus { major, minor } = self.primitive {
            let (topology, _) = super::torus::rebuilt_band(
                frame,
                major,
                minor,
                &self.plane,
                self.index,
                self.tolerance,
                operation,
            )?;
            return piece_solid(self.clone(), frame, topology, operation);
        }
        let [a, b, c, _] = &self.plane;
        let world = frame.x() * rational_f64(a)
            + frame.y() * rational_f64(b)
            + frame.normal() * rational_f64(c);
        let mut built = halves(&self.primitive, frame, self.tolerance, &self.plane, world)?;
        if self.index >= built.len() {
            return Err(Error::InvalidTopology("a split half rebuilt differently"));
        }
        let piece = built.swap_remove(self.index);
        let topology = Topology::from_parts(piece.parts.with_measured_enclosures(), self.tolerance)
            .map_err(|issues| Error::InvalidTopology(super::oblique::issue_name(&issues)))?;
        piece_solid(self.clone(), frame, topology, operation)
    }

    /// Inside where the primitive and the plane's side both hold; on the
    /// boundary within the resolution of either.
    pub(crate) fn classify(
        &self,
        local: [f64; 3],
        tolerance: Tolerance,
    ) -> Result<crate::Location> {
        use crate::Location;
        for x in local {
            crate::math::finite(x, "coordinate")?;
        }
        let tol = tolerance.linear();
        let inside = match self.primitive {
            Primitive::Cone {
                bottom,
                top,
                height,
            } => crate::decide::cone_location(local, bottom, top, height, tol),
            Primitive::Zone { radius, low, high } => {
                let [w0, w1] = [ends(&self.primitive)[0].w, ends(&self.primitive)[1].w];
                let _ = (low, high);
                crate::decide::sphere_location(local, radius, w0, w1, tol)
            }
            Primitive::Torus { major, minor } => {
                crate::decide::torus_location(local, major, minor, tol)
            }
        };
        let [a, b, c, d] = &self.plane;
        let f = a * q(local[0]) + b * q(local[1]) + c * q(local[2]) + d;
        let m2 = a * a + b * b + c * c;
        let near = &f * &f <= q(tol) * q(tol) * &m2;
        let below = f < zero();
        let own = self.index == 0;
        if inside == 2 || (below != own && !near) {
            return Ok(Location::Outside);
        }
        Ok(if inside == 1 || near {
            Location::Boundary
        } else {
            Location::Inside
        })
    }
}

fn piece_solid(
    half: Half,
    frame: Frame3,
    topology: Topology,
    operation: OperationId,
) -> Result<Solid> {
    let mass = topology
        .mass_enclosure()
        .ok_or(Error::Unrepresentable("a split half's mass properties"))?
        .midpoints();
    let bounds = super::oblique::edge_bounds(&topology);
    let [bottom, top] = ends(&half.primitive);
    Ok(Solid {
        construction: Construction::Half(Box::new(half)),
        frame,
        start: bottom.w,
        end: top.w,
        topology,
        mass,
        bounds,
        operation,
    })
}

impl Solid {
    /// S8c.2's halves; `None` when the plane does not contain the axis.
    pub(super) fn split_meridian(
        &self,
        context: &Context,
        plane: &Frame3,
        coefficients: &[R; 4],
    ) -> Result<Option<Pieces>> {
        let (primitive, tolerance) = match &self.construction {
            Construction::Cone {
                bottom,
                top,
                tolerance,
            } => (
                Primitive::Cone {
                    bottom: *bottom,
                    top: *top,
                    height: self.end,
                },
                *tolerance,
            ),
            Construction::Sphere {
                radius,
                low,
                high,
                tolerance,
            } => (
                Primitive::Zone {
                    radius: *radius,
                    low: *low,
                    high: *high,
                },
                *tolerance,
            ),
            _ => return Ok(None),
        };
        if !contains_axis(coefficients, self.start, self.end, tolerance) {
            return Ok(None);
        }
        let operation = context.operation;
        let built = halves(
            &primitive,
            self.frame,
            tolerance,
            coefficients,
            plane.normal(),
        )?;
        // The input's entities by role and end.
        let input = &self.topology;
        let body = input.body_id();
        let resolve = |role: Role, end: usize| {
            input_entity(input, role, end)
                .ok_or(Error::InvalidTopology("a revolved entity without its role"))
        };
        let mut relations = Vec::new();
        let mut split: BTreeMap<EntityId, Vec<EntityId>> = BTreeMap::new();
        let mut pieces = Vec::new();
        for (k, piece) in built.into_iter().enumerate() {
            let mut derivations = Vec::new();
            for (slot, plan) in &piece.plans {
                let entity = match slot {
                    Slot::Vertex(_) => EntityKind::Vertex,
                    Slot::Edge(_) => EntityKind::Edge,
                    Slot::Face(_) => EntityKind::Face,
                    Slot::Region(_) => EntityKind::Region,
                };
                let d = match plan {
                    Plan::Child(role, end) => {
                        let from = resolve(*role, *end)?;
                        let d = Derivation {
                            operation,
                            kind: OperationKind::PlaneSplit,
                            entity,
                            role: *role,
                            ordinal: k as u32,
                            parents: vec![Parent::Entity(from)],
                        };
                        split.entry(from).or_default().push(d.id());
                        d
                    }
                    Plan::New(from, role, ordinal) => {
                        let parents = from
                            .iter()
                            .map(|(r, e)| resolve(*r, *e).map(Parent::Entity))
                            .collect::<Result<Vec<_>>>()?;
                        let d = Derivation {
                            operation,
                            kind: OperationKind::PlaneSplit,
                            entity,
                            role: *role,
                            ordinal: *ordinal,
                            parents,
                        };
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
            let body_derivation = Derivation {
                operation,
                kind: OperationKind::PlaneSplit,
                entity: EntityKind::Body,
                role: Role::Body,
                ordinal: k as u32,
                parents: vec![Parent::Entity(body)],
            };
            let topology =
                Topology::from_parts_named(piece.parts, tolerance, body_derivation, derivations)
                    .map_err(|issues| {
                        Error::InvalidTopology(super::oblique::issue_name(&issues))
                    })?;
            let half = Half {
                primitive: primitive.clone(),
                tolerance,
                plane: coefficients.clone(),
                index: k,
            };
            pieces.push((
                piece.side,
                piece_solid(half, self.frame, topology, operation)?,
            ));
        }
        for (from, into) in split {
            relations.push(Relation::Split { from, into });
        }
        // Every input entity is split (a half keeps none whole).
        let covered: BTreeSet<EntityId> = relations.iter().flat_map(Relation::sources).collect();
        if input.ids().any(|(id, _)| !covered.contains(&id)) {
            return Err(Error::InvalidTopology(
                "a revolved entity no half accounts for",
            ));
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
        Ok(Some((pieces, history)))
    }
}
