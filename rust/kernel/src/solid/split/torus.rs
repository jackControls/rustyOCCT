//! S8d.1: a whole torus split by a plane normal to its axis or containing
//! it.
//!
//! In the torus's frame the plane is `F = a u + b v + c w + d`. Normal to
//! the axis (`a = b = 0`, or within a quarter of the resolution over the
//! outer equator), it cuts the tube at the height `h = -d / c` in two
//! parallels, the circles at the tube's angles `v1 = asin(h / r)` and `pi -
//! v1` (rounded): the piece above is the band of the tube between them
//! through its top, the piece below the rest, each closed by the planar
//! annulus between the parallels. Both are general bodies on the input's
//! own torus surface, the wall wound once in `u` with the parallels' `v`-
//! constant pcurves. Containing the axis, the plane cuts the tube in its
//! meridian circles on either side: the pieces are the two half-turn
//! wedges of the S3 construction on frames whose x axes point along the
//! plane's trace and against it. Other planes cut the tube in spiric curves
//! (S8d.3); tori other than whole ones wait too.
//!
//! Names follow provenance: the wall and the region are `Split` into one
//! child per piece (below first); everything else is new and `Generated`
//! from the wall (the cut faces, their circles).
use super::meridian::{Half, Primitive};
use super::{q, rational_f64, zero, Side};
use crate::certified::{Interval as I, Real};
use crate::history::{History, Relation};
use crate::identity::{Derivation, EntityId, EntityKind, OperationId, OperationKind, Parent, Role};
use crate::solid::{Construction, Context, Solid};
use crate::topology::{
    plane_pcurve, Curve2, Curve3, Edge, EdgeId, Face, FaceId, Fin, FinId, Loop, LoopId,
    Orientation, Region, RegionId, RegionKind, Shell, ShellId, Side as FaceSide, Slot, Surface,
    Topology, TopologyParts,
};
use crate::{Error, Frame3, Point2, Result, Tolerance};
use num_rational::BigRational as R;
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

/// The pieces with their sides, and the history.
type Pieces = (Vec<(Side, Solid)>, History);

/// A piece's parts and each slot's provenance: split from the input's wall
/// or region, or new from the wall.
pub(super) struct Band {
    pub(super) side: Side,
    pub(super) parts: TopologyParts,
    pub(super) plans: Vec<(Slot, Option<Role>)>,
}

/// The two bands of a whole torus cut at height `h` (below first):
/// `below_is_lower` when the side below the plane is the part under `h`.
pub(super) fn bands(
    frame: Frame3,
    major: f64,
    minor: f64,
    h: f64,
    below_is_lower: bool,
    tolerance: Tolerance,
) -> Result<Vec<Band>> {
    let v1 = (h / minor).asin();
    let (outer, inner) = (major + minor * v1.cos(), major - minor * v1.cos());
    let (n, x) = (frame.normal(), frame.x());
    let centre = frame.point(Point2::default(), h);
    let ring_frame = Frame3::new(centre, n, x, tolerance)?;
    let mut out = Vec::new();
    // The part above the plane along the axis: the tube from v1 through its
    // top to pi - v1; below: from pi - v1 through its bottom to 2 pi + v1.
    for upper in [false, true] {
        let (va, vb) = if upper {
            (v1, PI - v1)
        } else {
            (PI - v1, TAU + v1)
        };
        let mut parts = TopologyParts::default();
        let mut plans: Vec<(Slot, Option<Role>)> = Vec::new();
        // The rings: the one at va, the one at vb (outer at v1, inner at
        // pi - v1).
        let radius_at = |v: f64| {
            if v == v1 || v == TAU + v1 {
                outer
            } else {
                inner
            }
        };
        let mut rings = [EdgeId(0); 2];
        for (j, v) in [va, vb].into_iter().enumerate() {
            rings[j] = EdgeId(parts.edges.len());
            parts.edges.push(Edge {
                start: None,
                end: None,
                curve: Curve3::Circle {
                    frame: ring_frame,
                    radius: radius_at(v),
                },
                fins: Vec::new(),
            });
            plans.push((Slot::Edge(rings[j]), Some(Role::CutEdge)));
        }
        let mut add_face = |surface: Surface,
                            sense: Orientation,
                            loops: Vec<(Vec<Fin>, [i32; 2])>,
                            plan: Option<Role>,
                            parts: &mut TopologyParts| {
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
                sense,
                loops: ids,
                front: ShellId(0),
                back: ShellId(1),
                enclosure: None,
            });
            plans.push((Slot::Face(FaceId(parts.faces.len() - 1)), plan));
        };
        // The band: +u along the ring at va, -u along the one at vb (the
        // region between them on the left, the torus's normal outward).
        let line = |a: (f64, f64), b: (f64, f64)| Curve2::LineSegment {
            start: Point2::new(a.0, a.1),
            end: Point2::new(b.0, b.1),
        };
        let lower_loop = vec![Fin {
            edge: rings[0],
            sense: Orientation::Forward,
            pcurve: line((0.0, va), (TAU, va)),
            enclosure: None,
        }];
        let upper_loop = vec![Fin {
            edge: rings[1],
            sense: Orientation::Reversed,
            pcurve: line((TAU, vb), (0.0, vb)),
            enclosure: None,
        }];
        add_face(
            Surface::Torus {
                frame,
                major,
                minor,
            },
            Orientation::Forward,
            vec![(lower_loop, [1, 0]), (upper_loop, [-1, 0])],
            None,
            &mut parts,
        );
        // The annulus, facing out of the piece: down under the upper part,
        // up over the lower part; its outer loop first, counter-clockwise
        // about that normal, the inner clockwise.
        let normal = if upper { n * -1.0 } else { n };
        let disc = Frame3::new(centre, normal, x, tolerance)?;
        let (outer_ring, inner_ring) = if upper {
            (rings[0], rings[1])
        } else {
            (rings[1], rings[0])
        };
        let facing_up = !upper;
        let fin = |edge: EdgeId, ccw: bool, parts: &TopologyParts| {
            // A ring runs counter-clockwise about n.
            let sense = if ccw == facing_up {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            Fin {
                edge,
                sense,
                pcurve: plane_pcurve(&parts.edges[edge.0].curve, sense, disc),
                enclosure: None,
            }
        };
        let outer_fin = fin(outer_ring, true, &parts);
        let inner_fin = fin(inner_ring, false, &parts);
        add_face(
            Surface::Plane(disc),
            Orientation::Forward,
            vec![(vec![outer_fin], [0, 0]), (vec![inner_fin], [0, 0])],
            Some(Role::CutFace),
            &mut parts,
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
        plans.push((Slot::Region(RegionId(1)), None));
        let side = if upper == below_is_lower {
            Side::Above
        } else {
            Side::Below
        };
        out.push(Band { side, parts, plans });
    }
    out.sort_by_key(|b| b.side);
    Ok(out)
}

impl Solid {
    /// S8d.1's split of a whole torus; `None` when the solid is not one.
    pub(super) fn split_torus(
        &self,
        context: &Context,
        _plane: &Frame3,
        [a, b, c, d]: [R; 4],
    ) -> Result<Option<Pieces>> {
        let Construction::Torus {
            major,
            minor,
            low,
            high,
            angle,
            tolerance,
        } = self.construction
        else {
            return Ok(None);
        };
        if !(high - low == TAU && angle == TAU) {
            return Err(Error::OutOfDomain(
                "a torus band or wedge split by a plane (S8d)",
            ));
        }
        let tol = tolerance.linear();
        let ab2 = &a * &a + &b * &b;
        let m2 = &ab2 + &c * &c;
        // The plane misses or touches the tube when its distance from the
        // core circle `R (cos t, sin t, 0)` is at least `r` everywhere: the
        // core's `F` ranges over `d +- R |(a, b)|`, so when `|d| - R |(a,
        // b)| >= r |m|`, i.e. `d^2 - R^2 |ab|^2 - r^2 |m|^2 >= 2 R r |ab|
        // |m|`, decided by squares.
        let (rr, big) = (q(minor), q(major));
        let lhs = &d * &d - &big * &big * &ab2 - &rr * &rr * &m2;
        let four = R::from_integer(4.into());
        let whole = lhs >= zero() && &lhs * &lhs >= four * &big * &big * &rr * &rr * &ab2 * &m2;
        if whole {
            return Ok(Some(self.unchanged(context)));
        }
        let operation = context.operation;
        let widest = major + minor;
        let normal_to_axis = ab2 == zero()
            || &ab2 * q(widest) * q(widest) * R::from_integer(16.into())
                <= q(tol) * q(tol) * &c * &c;
        let contains_axis = [-minor, minor].iter().all(|w| {
            let f = &c * q(*w) + &d;
            &f * &f * R::from_integer(16.into()) <= q(tol) * q(tol) * &m2
        });
        if normal_to_axis {
            let h = rational_f64(&(-&d / &c));
            if !(-minor < h && h < minor) || minor - h.abs() <= tol {
                return Err(Error::Degenerate(
                    "a cut within the resolution of the tube's top",
                ));
            }
            let below_is_lower = c > zero();
            let built = bands(self.frame, major, minor, h, below_is_lower, tolerance)?;
            return self
                .name_torus_pieces(context, built, tolerance, &[a, b, c, d])
                .map(Some);
        }
        if contains_axis {
            // Two half-turn wedges: the lower side (F < 0) turns counter-
            // clockwise from the trace t = (-b, a) / |(a, b)| to -t.
            let ab = I::exact(ab2).sqrt();
            let unit = |x: &R| {
                let (lo, hi) = I::exact(x.clone())
                    .div(&ab)
                    .ok_or(Error::Degenerate("a plane normal to the axis"))?
                    .bounds_f64();
                Ok::<f64, Error>(0.5 * lo + 0.5 * hi)
            };
            let t = [unit(&-b.clone())?, unit(&a)?];
            let along = self.frame.x() * t[0] + self.frame.y() * t[1];
            let mut pieces = Vec::new();
            for (side, x) in [(Side::Below, along), (Side::Above, along * -1.0)] {
                let frame = Frame3::new(self.frame.origin(), self.frame.normal(), x, tolerance)?;
                let wedge =
                    Solid::build_torus(operation, frame, major, minor, 0.0, TAU, PI, tolerance)?;
                pieces.push((side, wedge));
            }
            return self.name_wedges(context, pieces).map(Some);
        }
        Err(Error::OutOfDomain(
            "a torus by a plane neither normal to nor containing its axis (S8d.3)",
        ))
    }

    /// Name general torus pieces: the wall and region split, the rest new
    /// from the wall.
    fn name_torus_pieces(
        &self,
        context: &Context,
        built: Vec<Band>,
        tolerance: Tolerance,
        plane: &[R; 4],
    ) -> Result<Pieces> {
        let operation = context.operation;
        let input = &self.topology;
        let find = |role: Role| {
            input
                .ids()
                .map(|(id, _)| id)
                .find(|id| input.derivation(*id).map(|d| d.role) == Some(role))
                .ok_or(Error::InvalidTopology("a torus without its wall or region"))
        };
        let (wall, region) = (find(Role::Wall)?, find(Role::Region)?);
        let Construction::Torus { major, minor, .. } = self.construction else {
            unreachable!("a torus")
        };
        let mut relations = Vec::new();
        let mut split: BTreeMap<EntityId, Vec<EntityId>> = BTreeMap::new();
        let mut pieces = Vec::new();
        for (k, band) in built.into_iter().enumerate() {
            let mut derivations = Vec::new();
            let mut local = 0u32;
            for (slot, plan) in &band.plans {
                let entity = match slot {
                    Slot::Vertex(_) => EntityKind::Vertex,
                    Slot::Edge(_) => EntityKind::Edge,
                    Slot::Face(_) => EntityKind::Face,
                    Slot::Region(_) => EntityKind::Region,
                };
                let d = match plan {
                    None => {
                        let from = if entity == EntityKind::Region {
                            region
                        } else {
                            wall
                        };
                        let role = if entity == EntityKind::Region {
                            Role::Region
                        } else {
                            Role::Wall
                        };
                        let d = Derivation {
                            operation,
                            kind: OperationKind::PlaneSplit,
                            entity,
                            role,
                            ordinal: k as u32,
                            parents: vec![Parent::Entity(from)],
                        };
                        split.entry(from).or_default().push(d.id());
                        d
                    }
                    Some(role) => {
                        let d = Derivation {
                            operation,
                            kind: OperationKind::PlaneSplit,
                            entity,
                            role: *role,
                            ordinal: (k as u32) * 16 + local,
                            parents: vec![Parent::Entity(wall)],
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
            let body_derivation = Derivation {
                operation,
                kind: OperationKind::PlaneSplit,
                entity: EntityKind::Body,
                role: Role::Body,
                ordinal: k as u32,
                parents: vec![Parent::Entity(input.body_id())],
            };
            let topology =
                Topology::from_parts_named(band.parts, tolerance, body_derivation, derivations)
                    .map_err(|issues| {
                        Error::InvalidTopology(super::oblique::issue_name(&issues))
                    })?;
            let half = Half::new(
                Primitive::Torus { major, minor },
                tolerance,
                plane.clone(),
                k,
            );
            pieces.push((band.side, half.solid(self.frame, topology, operation)?));
        }
        for (from, into) in split {
            relations.push(Relation::Split { from, into });
        }
        self.finish_split(context, pieces, relations)
    }

    /// Name the wedges built by the S3 construction: the wall and region
    /// split, the discs and circles new from the wall.
    fn name_wedges(&self, context: &Context, mut pieces: Vec<(Side, Solid)>) -> Result<Pieces> {
        let operation = context.operation;
        let input = &self.topology;
        let find = |role: Role| {
            input
                .ids()
                .map(|(id, _)| id)
                .find(|id| input.derivation(*id).map(|d| d.role) == Some(role))
                .ok_or(Error::InvalidTopology("a torus without its wall or region"))
        };
        let (wall, region) = (find(Role::Wall)?, find(Role::Region)?);
        let mut relations = Vec::new();
        let mut split: BTreeMap<EntityId, Vec<EntityId>> = BTreeMap::new();
        for (k, (_, solid)) in pieces.iter_mut().enumerate() {
            let t = &solid.topology;
            let mut derivations = Vec::new();
            let mut local = 0u32;
            for (id, slot) in t.ids().collect::<Vec<_>>() {
                let d = t.derivation(id).expect("a derivation");
                let named = match d.role {
                    Role::Wall | Role::Region => {
                        let from = if d.role == Role::Wall { wall } else { region };
                        let child = Derivation {
                            operation,
                            kind: OperationKind::PlaneSplit,
                            entity: d.entity,
                            role: d.role,
                            ordinal: k as u32,
                            parents: vec![Parent::Entity(from)],
                        };
                        split.entry(from).or_default().push(child.id());
                        child
                    }
                    _ => {
                        let role = match d.entity {
                            EntityKind::Face => Role::CutFace,
                            EntityKind::Edge => Role::CutEdge,
                            _ => Role::CutVertex,
                        };
                        let new = Derivation {
                            operation,
                            kind: OperationKind::PlaneSplit,
                            entity: d.entity,
                            role,
                            ordinal: (k as u32) * 16 + local,
                            parents: vec![Parent::Entity(wall)],
                        };
                        local += 1;
                        relations.push(Relation::Generated {
                            from: new.parents.clone(),
                            to: new.id(),
                            role,
                        });
                        new
                    }
                };
                derivations.push((slot, named));
            }
            let body_derivation = Derivation {
                operation,
                kind: OperationKind::PlaneSplit,
                entity: EntityKind::Body,
                role: Role::Body,
                ordinal: k as u32,
                parents: vec![Parent::Entity(input.body_id())],
            };
            solid.topology = solid
                .topology
                .clone()
                .renamed(body_derivation, derivations)?;
            solid.operation = operation;
        }
        for (from, into) in split {
            relations.push(Relation::Split { from, into });
        }
        self.finish_split(context, pieces, relations)
    }

    /// The history of a split whose relations are complete, the pieces'
    /// enclosures carried, attributes propagated and both checked.
    pub(super) fn finish_split(
        &self,
        context: &Context,
        mut pieces: Vec<(Side, Solid)>,
        mut relations: Vec<Relation>,
    ) -> Result<Pieces> {
        let operation = context.operation;
        relations.sort_by_cached_key(Relation::sort_key);
        relations.dedup();
        let outputs: Vec<EntityId> = pieces.iter().map(|(_, s)| s.topology.body_id()).collect();
        let history = History::new(
            operation,
            OperationKind::PlaneSplit,
            vec![self.topology.body_id()],
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

/// The bands of a whole torus rebuilt in another frame (a rigid motion):
/// the piece at `index`, its ids external.
pub(super) fn rebuilt_band(
    frame: Frame3,
    major: f64,
    minor: f64,
    plane: &[R; 4],
    index: usize,
    tolerance: Tolerance,
    operation: OperationId,
) -> Result<(Topology, Side)> {
    let [_, _, c, d] = plane;
    let h = rational_f64(&(-d / c));
    let mut built = bands(frame, major, minor, h, *c > zero(), tolerance)?;
    if index >= built.len() {
        return Err(Error::InvalidTopology("a torus band rebuilt differently"));
    }
    let band = built.swap_remove(index);
    let topology = Topology::from_parts(band.parts.with_measured_enclosures(), tolerance)
        .map_err(|issues| Error::InvalidTopology(super::oblique::issue_name(&issues)))?;
    let _ = operation;
    Ok((topology, band.side))
}
