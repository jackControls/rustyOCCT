//! `.brep` solids into the cell model, by the rule of
//! `rust/tools/cell_reference.py::to_cell`: a seam pair on a cylinder whose
//! uses run opposite ways, lie one period apart and continue their
//! neighbours in UV merges into periodic loops with winding numbers; the
//! vertices only seams and closed curves use disappear, and closed curves
//! that lose their vertex become ring edges. Each solid becomes a solid
//! region whose outer shell comes first, each shell's opposite sides a void
//! twin. The result must pass `Topology::from_parts`.
use super::read::{self, Data, Document, EdgeRep, Kind, Orient, Sub};
use super::Transform;
use crate::topology::{
    Curve2, Curve3, Edge, EdgeId, Enclosure, Face, FaceId, Fin, FinId, Issue, Loop, LoopId,
    Orientation, Region, RegionId, RegionKind, Shell, ShellId, Side, Surface, Topology,
    TopologyParts, Vertex, VertexId,
};
use crate::{Frame3, Point2, Point3, Tolerance, Vec3};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

/// One solid of a document.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportedSolid {
    /// The solid's shape record (0-based, file order).
    pub record: usize,
    /// The body resolution: the largest vertex, edge or face tolerance of
    /// the solid (per-entity tolerances become enclosures in M5).
    pub tolerance: Tolerance,
    /// The cell topology, or why there is none: unsupported constructs by
    /// name, or the validation issues of the converted cells.
    pub result: Result<Topology, Rejected>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Rejected {
    Unsupported(Vec<&'static str>),
    /// The converted cells, kept for diagnosis, and every validation issue.
    Invalid {
        issues: Vec<Issue>,
        parts: Box<TopologyParts>,
    },
}

/// Everything a document holds, as far as the kernel represents it.
#[derive(Debug, Clone, PartialEq)]
pub struct Import {
    pub solids: Vec<ImportedSolid>,
    /// Every construct the kernel cannot represent, by name, with the number
    /// of records (geometry) or shape uses (topology) that carry it.
    pub unsupported: BTreeMap<&'static str, usize>,
}

fn compose(a: Orient, b: Orient) -> Option<Orient> {
    match (a, b) {
        (Orient::Forward, x) | (x, Orient::Forward) => Some(x),
        (Orient::Reversed, Orient::Reversed) => Some(Orient::Forward),
        _ => None,
    }
}

fn p3(v: [f64; 3]) -> Point3 {
    Point3::new(v[0], v[1], v[2])
}

fn v3(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn same_transform(a: &Transform, b: &Transform) -> bool {
    a.0.iter()
        .zip(&b.0)
        .all(|(x, y)| (x - y).abs() <= 1e-12 * (1.0 + x.abs().max(y.abs())))
}

/// Shape instances by record and placement: the same record reached through
/// different location chains is one instance when the composed placements
/// agree to rounding.
#[derive(Default)]
struct Instances(BTreeMap<usize, Vec<(Transform, usize)>>);

impl Instances {
    fn get(&self, record: usize, t: &Transform) -> Option<usize> {
        self.0
            .get(&record)?
            .iter()
            .find(|(u, _)| same_transform(u, t))
            .map(|(_, i)| *i)
    }
    fn insert(&mut self, record: usize, t: &Transform, index: usize) {
        self.0.entry(record).or_default().push((*t, index));
    }
}

/// A seamed model: OCCT's structure before the seam merge.
struct SUse {
    edge: usize,
    forward: bool,
    pcurve: Curve2,
}
struct SEdge {
    start: usize,
    end: usize,
    curve: Curve3,
    /// OCCT's stored edge tolerance.
    tolerance: f64,
}
struct SFace {
    surface: Surface,
    forward: bool,
    loops: Vec<Vec<SUse>>,
}

struct Walk<'a> {
    doc: &'a Document,
    placement: Tolerance,
    unsupported: Vec<&'static str>,
    tolerance: f64,
    vertices: Vec<Point3>,
    /// OCCT's stored tolerance of each vertex.
    vertex_tolerances: Vec<f64>,
    vertex_keys: Instances,
    edges: Vec<SEdge>,
    edge_keys: Instances,
    faces: Vec<SFace>,
    shells: Vec<Vec<usize>>,
}

fn location(doc: &Document, index: usize) -> Transform {
    if index == 0 {
        Transform::IDENTITY
    } else {
        doc.locations[index - 1]
    }
}

impl Walk<'_> {
    fn no<T>(&mut self, what: &'static str) -> Option<T> {
        self.unsupported.push(what);
        None
    }

    fn vertex(&mut self, record: usize, t: &Transform) -> Option<usize> {
        if let Some(v) = self.vertex_keys.get(record, t) {
            return Some(v);
        }
        let Data::Vertex { tolerance, point } = self.doc.shapes[record].data else {
            return self.no("MalformedVertex");
        };
        self.tolerance = self.tolerance.max(tolerance);
        let v = self.vertices.len();
        self.vertices.push(p3(t.point(point)));
        self.vertex_tolerances.push(tolerance);
        self.vertex_keys.insert(record, t, v);
        Some(v)
    }

    fn edge(&mut self, record: usize, t: &Transform) -> Option<usize> {
        if let Some(e) = self.edge_keys.get(record, t) {
            return Some(e);
        }
        let doc = self.doc;
        let shape = &doc.shapes[record];
        let Data::Edge {
            tolerance,
            degenerated,
            reps,
        } = &shape.data
        else {
            return self.no("MalformedEdge");
        };
        if *degenerated {
            return self.no("DegeneratedEdge");
        }
        self.tolerance = self.tolerance.max(*tolerance);
        let (mut start, mut end) = (None, None);
        for sub in &shape.subs {
            let lt = t.times(&location(doc, sub.location));
            let v = self.vertex(sub.shape, &lt)?;
            match sub.orient {
                Orient::Forward => start = Some(v),
                Orient::Reversed => end = Some(v),
                _ => return self.no("InternalOrExternalVertex"),
            }
        }
        let (Some(start), Some(end)) = (start, end) else {
            return self.no("EdgeWithoutBothVertices");
        };
        let Some((curve, location_index, range)) = reps.iter().find_map(|r| match r {
            EdgeRep::Curve {
                curve,
                location,
                range,
            } => Some((*curve, *location, *range)),
            _ => None,
        }) else {
            return self.no("EdgeWithout3DCurve");
        };
        let ct = t.times(&location(doc, location_index));
        let Some(record_curve) = curve.checked_sub(1).and_then(|c| doc.curves.get(c)) else {
            return self.no("MissingCurve");
        };
        let [f, l] = range;
        let curve = match record_curve {
            read::Curve3::Line { p, d } => Curve3::LineSegment {
                start: p3(ct.point(std::array::from_fn(|i| p[i] + f * d[i]))),
                end: p3(ct.point(std::array::from_fn(|i| p[i] + l * d[i]))),
            },
            read::Curve3::Circle { p, x, y, r, .. } => {
                let (x, y) = (ct.vector(*x), ct.vector(*y));
                let frame =
                    Frame3::new(p3(ct.point(*p)), v3(cross(x, y)), v3(x), self.placement).ok()?;
                // A closed edge is a whole circle: its sweep is exactly one turn.
                let sweep = if start == end && ((l - f) - TAU).abs() <= 1e-12 * TAU {
                    TAU
                } else {
                    l - f
                };
                Curve3::CircularArc {
                    frame,
                    radius: *r,
                    start_angle: f,
                    sweep_angle: sweep,
                }
            }
            read::Curve3::Other(name) => return self.no(name),
        };
        let e = self.edges.len();
        self.edges.push(SEdge {
            start,
            end,
            curve,
            tolerance: *tolerance,
        });
        self.edge_keys.insert(record, t, e);
        Some(e)
    }

    /// The pcurve of edge record `record` (at `t`, our edge `edge`) on the
    /// face's surface, in the edge's direction. A plane may store none: its
    /// pcurve is then the edge's curve in the plane's coordinates, as
    /// `BRep_Tool::CurveOnPlane` derives it.
    #[allow(clippy::too_many_arguments)]
    fn pcurve(
        &mut self,
        record: usize,
        t: &Transform,
        edge: usize,
        surface: &Surface,
        face_surface: usize,
        face_location: &Transform,
        stored: Orient,
    ) -> Option<Curve2> {
        let doc = self.doc;
        let Data::Edge { reps, .. } = &doc.shapes[record].data else {
            return None;
        };
        let rep = reps.iter().find_map(|r| match r {
            EdgeRep::OnSurface {
                pcurves,
                surface,
                location,
                range,
            } if *surface == face_surface
                && same_transform(&t.times(&location_of(doc, *location)), face_location) =>
            {
                Some((pcurves.clone(), *range))
            }
            _ => None,
        });
        let Some((pcurves, [f, mut l])) = rep else {
            return match surface {
                Surface::Plane(plane) => Some(on_plane(plane, &self.edges[edge].curve)),
                _ => self.no("EdgeWithoutPCurve"),
            };
        };
        // A closed circle's sweep was snapped to one turn; its pcurves share
        // the edge's range (SameRange), so theirs is snapped alike.
        if matches!(self.edges[edge].curve, Curve3::CircularArc { sweep_angle, .. } if sweep_angle == TAU)
            && ((l - f) - TAU).abs() <= 1e-12 * TAU
        {
            l = f + TAU;
        }
        // BRep_Tool::CurveOnSurface: the second pcurve of a closed surface
        // serves the edge's reversed use in the unoriented face.
        let index = if pcurves.len() == 2 && stored == Orient::Reversed {
            pcurves[1]
        } else {
            pcurves[0]
        };
        let Some(record_curve) = index.checked_sub(1).and_then(|i| doc.curves2d.get(i)) else {
            return self.no("MissingPCurve");
        };
        Some(match record_curve {
            read::Curve2::Line { p, d } => Curve2::LineSegment {
                start: Point2::new(p[0] + f * d[0], p[1] + f * d[1]),
                end: Point2::new(p[0] + l * d[0], p[1] + l * d[1]),
            },
            read::Curve2::Circle { c, x, y, r } => {
                let base = x[1].atan2(x[0]);
                let turn = if x[0] * y[1] - x[1] * y[0] > 0.0 {
                    1.0
                } else {
                    -1.0
                };
                Curve2::CircularArc {
                    center: Point2::new(c[0], c[1]),
                    radius: *r,
                    start_angle: base + turn * f,
                    sweep_angle: turn * (l - f),
                }
            }
            read::Curve2::Other(name) => return self.no(name),
        })
    }

    fn face(&mut self, record: usize, t: &Transform, oriented: Orient) -> Option<usize> {
        let doc = self.doc;
        let shape = &doc.shapes[record];
        let (tolerance, surface_index, location_index) = match shape.data {
            Data::Face {
                tolerance,
                surface,
                location,
            } => (tolerance, surface, location),
            Data::MeshFace => return self.no("TriangulationOnlyFace"),
            _ => return self.no("MalformedFace"),
        };
        self.tolerance = self.tolerance.max(tolerance);
        let st = t.times(&location(doc, location_index));
        let Some(record_surface) = surface_index
            .checked_sub(1)
            .and_then(|s| doc.surfaces.get(s))
        else {
            return self.no("FaceWithoutSurface");
        };
        let mut indirect = false;
        let surface = match record_surface {
            read::Surface::Plane { p, x, y, .. } => {
                let (x, y) = (st.vector(*x), st.vector(*y));
                Surface::Plane(
                    Frame3::new(p3(st.point(*p)), v3(cross(x, y)), v3(x), self.placement).ok()?,
                )
            }
            read::Surface::Cylinder { p, n, x, y, r } => {
                // An indirect axis (Y = X x N) is the kernel's cylinder about
                // -N with v negated; its surface normal points inward, so the
                // face's sense flips and the oriented normal is unchanged.
                indirect = dot(cross(*x, *y), *n) <= 0.0;
                let axis = if indirect { n.map(|c| -c) } else { *n };
                Surface::Cylinder {
                    frame: Frame3::new(
                        p3(st.point(*p)),
                        v3(st.vector(axis)),
                        v3(st.vector(*x)),
                        self.placement,
                    )
                    .ok()?,
                    radius: *r,
                }
            }
            read::Surface::Other(name) => return self.no(name),
        };
        let mut loops = Vec::new();
        for wire in &shape.subs {
            if doc.shapes[wire.shape].kind != Kind::Wire {
                return self.no("FaceWithNonWireChild");
            }
            let wt = t.times(&location(doc, wire.location));
            let mut uses = Vec::new();
            for e in &doc.shapes[wire.shape].subs {
                let et = wt.times(&location(doc, e.location));
                let stored =
                    compose(e.orient, wire.orient).or_else(|| self.no("InternalOrExternalEdge"))?;
                let traversal = compose(stored, oriented).expect("forward or reversed");
                let edge = self.edge(e.shape, &et)?;
                let pcurve =
                    self.pcurve(e.shape, &et, edge, &surface, surface_index, &st, stored)?;
                let pcurve = if indirect { negate_v(&pcurve) } else { pcurve };
                let forward = traversal == Orient::Forward;
                uses.push(SUse {
                    edge,
                    forward,
                    pcurve: if forward { pcurve } else { reversed(&pcurve) },
                });
            }
            loops.push(self.in_traversal_order(uses));
        }
        // On a plane the outer loop comes first; OCCT stores wires in any
        // order. The outer loop encloses the largest area.
        if matches!(surface, Surface::Plane(_)) && loops.len() > 1 {
            let area = |lp: &Vec<SUse>| {
                lp.iter()
                    .map(|u| {
                        (0..16)
                            .map(|k| {
                                let (a, b) = (
                                    u.pcurve.point(k as f64 / 16.0),
                                    u.pcurve.point((k + 1) as f64 / 16.0),
                                );
                                a.x * b.y - b.x * a.y
                            })
                            .sum::<f64>()
                    })
                    .sum::<f64>()
                    .abs()
            };
            let outer = (0..loops.len())
                .max_by(|a, b| area(&loops[*a]).total_cmp(&area(&loops[*b])))
                .unwrap_or(0);
            loops.swap(0, outer);
        }
        let f = self.faces.len();
        self.faces.push(SFace {
            surface,
            forward: (oriented == Orient::Forward) != indirect,
            loops,
        });
        Some(f)
    }
}

impl Walk<'_> {
    /// A wire's uses in traversal order: OCCT does not require a wire to
    /// store its edges in order. Each next use starts where the last ended;
    /// among several (a seam vertex), the one continuing in UV.
    fn in_traversal_order(&self, mut uses: Vec<SUse>) -> Vec<SUse> {
        let ends = |u: &SUse| {
            let e = &self.edges[u.edge];
            if u.forward {
                (e.start, e.end)
            } else {
                (e.end, e.start)
            }
        };
        let mut out: Vec<SUse> = Vec::with_capacity(uses.len());
        if uses.is_empty() {
            return out;
        }
        out.push(uses.remove(0));
        while !uses.is_empty() {
            let last = out.last().expect("not empty");
            let (_, at) = ends(last);
            let here = last.pcurve.point(1.0);
            let next = (0..uses.len())
                .filter(|k| ends(&uses[*k]).0 == at)
                .min_by(|a, b| {
                    let d = |k: usize| {
                        let p = uses[k].pcurve.point(0.0);
                        (p.x - here.x).powi(2) + (p.y - here.y).powi(2)
                    };
                    d(*a).total_cmp(&d(*b))
                });
            match next {
                Some(k) => out.push(uses.remove(k)),
                // Not a cycle: keep the stored order for the validator.
                None => out.append(&mut uses),
            }
        }
        out
    }
}

/// The pcurve of a curve lying in a plane, in the curve's direction.
fn on_plane(plane: &Frame3, curve: &Curve3) -> Curve2 {
    let uv = |p: Point3| {
        let [x, y, _] = plane.coordinates(p);
        Point2::new(x, y)
    };
    match curve {
        Curve3::LineSegment { start, end } => Curve2::LineSegment {
            start: uv(*start),
            end: uv(*end),
        },
        Curve3::Circle { frame, radius } => on_plane(
            plane,
            &Curve3::CircularArc {
                frame: *frame,
                radius: *radius,
                start_angle: 0.0,
                sweep_angle: TAU,
            },
        ),
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let (x, n) = (frame.x(), plane.normal());
            let a = x.dot(plane.x());
            let b = x.dot(plane.normal().cross(plane.x()));
            let turn = if frame.normal().dot(n) > 0.0 {
                1.0
            } else {
                -1.0
            };
            let base = b.atan2(a);
            Curve2::CircularArc {
                center: uv(frame.origin()),
                radius: *radius,
                start_angle: base + turn * start_angle,
                sweep_angle: turn * sweep_angle,
            }
        }
    }
}

fn negate_v(p: &Curve2) -> Curve2 {
    match p {
        Curve2::LineSegment { start, end } => Curve2::LineSegment {
            start: Point2::new(start.x, -start.y),
            end: Point2::new(end.x, -end.y),
        },
        Curve2::CircularArc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => Curve2::CircularArc {
            center: Point2::new(center.x, -center.y),
            radius: *radius,
            start_angle: -start_angle,
            sweep_angle: -sweep_angle,
        },
    }
}

fn location_of(doc: &Document, index: usize) -> Transform {
    location(doc, index)
}

fn reversed(p: &Curve2) -> Curve2 {
    match p {
        Curve2::LineSegment { start, end } => Curve2::LineSegment {
            start: *end,
            end: *start,
        },
        Curve2::CircularArc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => Curve2::CircularArc {
            center: *center,
            radius: *radius,
            start_angle: start_angle + sweep_angle,
            sweep_angle: -sweep_angle,
        },
    }
}

fn start_of(p: &Curve2) -> Point2 {
    p.point(0.0)
}

fn end_of(p: &Curve2) -> Point2 {
    p.point(1.0)
}

/// A seamed loop split at its seams: the seam edges, and the runs of use
/// indices between them with their windings.
type SeamRuns = (Vec<usize>, Vec<(Vec<usize>, i32)>);

/// Split one seamed loop into runs between seam uses, each with its winding;
/// `None` when the loop has no consistent seam pair.
fn seam_merge(face: &SFace, lp: &[SUse], tol: f64) -> Option<SeamRuns> {
    let Surface::Cylinder { radius, .. } = face.surface else {
        return None;
    };
    let close =
        |a: Point2, b: Point2| ((a.x - b.x) * radius).abs() <= tol && (a.y - b.y).abs() <= tol;
    let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
    for u in lp {
        *counts.entry(u.edge).or_default() += 1;
    }
    let seams: Vec<usize> = counts
        .iter()
        .filter(|(_, n)| **n == 2)
        .map(|(e, _)| *e)
        .collect();
    if seams.is_empty() {
        return None;
    }
    let n = lp.len();
    for &e in &seams {
        let at: Vec<usize> = (0..n).filter(|k| lp[*k].edge == e).collect();
        let (ua, ub) = (&lp[at[0]], &lp[at[1]]);
        if ua.forward == ub.forward {
            return None;
        }
        let (
            Curve2::LineSegment { start: sa, end: ea },
            Curve2::LineSegment { start: sb, end: eb },
        ) = (&ua.pcurve, &ub.pcurve)
        else {
            return None;
        };
        let shift = sa.x - eb.x;
        if (shift.abs() - TAU).abs() * radius > tol
            || ((ea.x - sb.x) - shift).abs() * radius > tol
            || (sa.y - eb.y).abs() > tol
            || (ea.y - sb.y).abs() > tol
        {
            return None;
        }
        for k in at {
            let (prev, next) = (&lp[(k + n - 1) % n], &lp[(k + 1) % n]);
            if !close(end_of(&prev.pcurve), start_of(&lp[k].pcurve))
                || !close(end_of(&lp[k].pcurve), start_of(&next.pcurve))
            {
                return None;
            }
        }
    }
    let first = (0..n).find(|k| seams.contains(&lp[*k].edge))?;
    let mut runs = Vec::new();
    let mut run: Vec<usize> = Vec::new();
    for step in 1..=n {
        let k = (first + step) % n;
        if seams.contains(&lp[k].edge) {
            if !run.is_empty() {
                runs.push(std::mem::take(&mut run));
            }
        } else {
            run.push(k);
        }
    }
    if !run.is_empty() {
        runs.push(run);
    }
    if runs.is_empty() {
        return None;
    }
    let mut out = Vec::new();
    for run in runs {
        let du = end_of(&lp[run[run.len() - 1]].pcurve).x - start_of(&lp[run[0]].pcurve).x;
        let w = (du / TAU).round();
        if w == 0.0 || (du - w * TAU).abs() * radius > tol {
            return None;
        }
        out.push((run, w as i32));
    }
    Some((seams, out))
}

/// The cell complex of a seamed solid, by the rule of `to_cell`.
fn to_cell(walk: Walk, tol: f64) -> TopologyParts {
    let mut face_loops: Vec<Vec<(Vec<&SUse>, i32)>> = Vec::new();
    let mut removed = BTreeSet::new();
    for f in &walk.faces {
        let mut loops = Vec::new();
        for lp in &f.loops {
            match (!lp.is_empty()).then(|| seam_merge(f, lp, tol)).flatten() {
                Some((seams, runs)) => {
                    removed.extend(seams);
                    loops.extend(
                        runs.into_iter()
                            .map(|(run, w)| (run.into_iter().map(|k| &lp[k]).collect(), w)),
                    );
                }
                None => loops.push((lp.iter().collect(), 0)),
            }
        }
        face_loops.push(loops);
    }
    let still_used: BTreeSet<usize> = face_loops
        .iter()
        .flatten()
        .flat_map(|(lp, _)| lp.iter().map(|u| u.edge))
        .collect();
    removed.retain(|e| !still_used.contains(e));
    let closed = |c: &Curve3| matches!(c, Curve3::CircularArc { sweep_angle, .. } if sweep_angle.abs() == TAU);
    let mut other_use = BTreeSet::new();
    for (i, e) in walk.edges.iter().enumerate() {
        if !(removed.contains(&i) || closed(&e.curve) && e.start == e.end) {
            other_use.extend([e.start, e.end]);
        }
    }
    let drop: BTreeSet<usize> = removed
        .iter()
        .flat_map(|&i| [walk.edges[i].start, walk.edges[i].end])
        .filter(|v| !other_use.contains(v))
        .collect();
    let mut parts = TopologyParts::default();
    let mut vertex_map = BTreeMap::new();
    // OCCT accepts a vertex within the larger of its own and each edge's
    // tolerance (BRepCheck_Vertex); that is the imported claim.
    let mut claim = walk.vertex_tolerances.clone();
    for (i, e) in walk.edges.iter().enumerate() {
        if !removed.contains(&i) {
            for v in [e.start, e.end] {
                claim[v] = claim[v].max(e.tolerance);
            }
        }
    }
    for (v, p) in walk.vertices.iter().enumerate() {
        if !drop.contains(&v) {
            vertex_map.insert(v, VertexId::new(parts.vertices.len()));
            parts.vertices.push(Vertex {
                position: *p,
                enclosure: Some(Enclosure::imported(claim[v])),
            });
        }
    }
    let mut edge_map = BTreeMap::new();
    for (i, e) in walk.edges.iter().enumerate() {
        if removed.contains(&i) {
            continue;
        }
        edge_map.insert(i, EdgeId::new(parts.edges.len()));
        let ring = closed(&e.curve) && drop.contains(&e.start);
        parts.edges.push(Edge {
            start: if ring {
                None
            } else {
                vertex_map.get(&e.start).copied()
            },
            end: if ring {
                None
            } else {
                vertex_map.get(&e.end).copied()
            },
            curve: e.curve.clone(),
            fins: Vec::new(),
        });
    }
    // Old shell k is solid shell k (front sides); its twin follows the
    // solid shells: the outer shell's twin bounds the infinite void, a
    // cavity's a bounded void region.
    let n = walk.shells.len();
    let mut owner = BTreeMap::new();
    for (k, s) in walk.shells.iter().enumerate() {
        for f in s {
            owner.entry(*f).or_insert(k);
        }
    }
    for (fi, f) in walk.faces.iter().enumerate() {
        let k = owner.get(&fi).copied().unwrap_or(0);
        let mut loops = Vec::new();
        for (lp, w) in &face_loops[fi] {
            let fins = lp
                .iter()
                .map(|u| {
                    parts.fins.push(Fin {
                        edge: edge_map[&u.edge],
                        sense: if u.forward {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        },
                        pcurve: u.pcurve.clone(),
                        // BRepCheck_Edge validates every use against the
                        // edge tolerance.
                        enclosure: Some(Enclosure::imported(walk.edges[u.edge].tolerance)),
                    });
                    FinId::new(parts.fins.len() - 1)
                })
                .collect();
            parts.loops.push(Loop::Edges {
                fins,
                winding: [*w, 0],
            });
            loops.push(LoopId::new(parts.loops.len() - 1));
        }
        parts.faces.push(Face {
            surface: f.surface.clone(),
            sense: if f.forward {
                Orientation::Forward
            } else {
                Orientation::Reversed
            },
            loops,
            front: ShellId::new(k),
            back: ShellId::new(n + k),
            // Computed after conversion: OCCT stores no bound on UV closure.
            enclosure: None,
        });
    }
    for (k, fin) in parts.fins.iter().enumerate() {
        parts.edges[fin.edge.index()].fins.push(FinId::new(k));
    }
    parts.regions = vec![
        Region {
            kind: RegionKind::Void,
            shells: vec![ShellId::new(n)],
        },
        Region {
            kind: RegionKind::Solid,
            shells: (0..n).map(ShellId::new).collect(),
        },
    ];
    for s in &walk.shells {
        parts.shells.push(Shell {
            region: RegionId::new(1),
            sides: s.iter().map(|f| (FaceId::new(*f), Side::Front)).collect(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        });
    }
    for (k, s) in walk.shells.iter().enumerate() {
        let region = if k == 0 {
            0
        } else {
            parts.regions.push(Region {
                kind: RegionKind::Void,
                shells: vec![ShellId::new(n + k)],
            });
            parts.regions.len() - 1
        };
        parts.shells.push(Shell {
            region: RegionId::new(region),
            sides: s.iter().map(|f| (FaceId::new(*f), Side::Back)).collect(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        });
    }
    parts
}

fn import_solid(doc: &Document, record: usize, t: &Transform, oriented: Orient) -> ImportedSolid {
    let mut walk = Walk {
        doc,
        placement: Tolerance::default(),
        unsupported: Vec::new(),
        tolerance: 0.0,
        vertices: Vec::new(),
        vertex_tolerances: Vec::new(),
        vertex_keys: Instances::default(),
        edges: Vec::new(),
        edge_keys: Instances::default(),
        faces: Vec::new(),
        shells: Vec::new(),
    };
    for shell in &doc.shapes[record].subs {
        if doc.shapes[shell.shape].kind != Kind::Shell {
            walk.unsupported.push("SolidWithNonShellChild");
            continue;
        }
        let Some(so) = compose(shell.orient, oriented) else {
            walk.unsupported.push("InternalOrExternalShell");
            continue;
        };
        let st = t.times(&location(doc, shell.location));
        let mut faces = Vec::new();
        for face in &doc.shapes[shell.shape].subs {
            let Some(fo) = compose(face.orient, so) else {
                walk.unsupported.push("InternalOrExternalFace");
                continue;
            };
            if doc.shapes[face.shape].kind != Kind::Face {
                walk.unsupported.push("ShellWithNonFaceChild");
                continue;
            }
            let ft = st.times(&location(doc, face.location));
            if let Some(f) = walk.face(face.shape, &ft, fo) {
                faces.push(f);
            }
        }
        walk.shells.push(faces);
    }
    if walk.shells.is_empty() {
        walk.unsupported.push("SolidWithoutShell");
    }
    let resolution = Tolerance::new(walk.tolerance.max(Tolerance::default().linear()), 1e-12)
        .unwrap_or_default();
    if !walk.unsupported.is_empty() {
        let mut names = std::mem::take(&mut walk.unsupported);
        names.sort_unstable();
        names.dedup();
        return ImportedSolid {
            record,
            tolerance: resolution,
            result: Err(Rejected::Unsupported(names)),
        };
    }
    // The outer shell comes first: the one whose box holds every other's.
    if walk.shells.len() > 1 {
        let boxes: Vec<[f64; 6]> = walk
            .shells
            .iter()
            .map(|s| {
                let mut b = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
                for f in s {
                    for u in walk.faces[*f].loops.iter().flatten() {
                        for v in [walk.edges[u.edge].start, walk.edges[u.edge].end] {
                            let p = walk.vertices[v].to_array();
                            for i in 0..3 {
                                b[i] = b[i].min(p[i]);
                                b[i + 3] = b[i + 3].max(p[i]);
                            }
                        }
                    }
                }
                b
            })
            .collect();
        let volume =
            |b: &[f64; 6]| (b[3] - b[0]).max(0.0) * (b[4] - b[1]).max(0.0) * (b[5] - b[2]).max(0.0);
        let outer = (0..boxes.len())
            .max_by(|a, b| volume(&boxes[*a]).total_cmp(&volume(&boxes[*b])))
            .unwrap_or(0);
        walk.shells.swap(0, outer);
    }
    let parts = to_cell(walk, resolution.linear()).with_measured_enclosures();
    ImportedSolid {
        record,
        tolerance: resolution,
        result: Topology::from_parts(parts.clone(), resolution).map_err(|issues| {
            Rejected::Invalid {
                issues,
                parts: Box::new(parts),
            }
        }),
    }
}

/// Every solid of a document in the cell model, and every construct the
/// kernel cannot represent, by name.
pub fn import(doc: &Document) -> Import {
    let mut unsupported: BTreeMap<&'static str, usize> = BTreeMap::new();
    for c in &doc.curves {
        if let read::Curve3::Other(name) = c {
            *unsupported.entry(name).or_default() += 1;
        }
    }
    for c in &doc.curves2d {
        if let read::Curve2::Other(name) = c {
            *unsupported.entry(name).or_default() += 1;
        }
    }
    for s in &doc.surfaces {
        if let read::Surface::Other(name) = s {
            *unsupported.entry(name).or_default() += 1;
        }
    }
    let geometry: Vec<&'static str> = unsupported.keys().copied().collect();
    let mut solids = Vec::new();
    let mut stack: Vec<(Sub, Transform, Orient)> =
        vec![(doc.root, Transform::IDENTITY, Orient::Forward)];
    while let Some((sub, parent, o)) = stack.pop() {
        let t = parent.times(&location(doc, sub.location));
        let Some(o) = compose(sub.orient, o) else {
            *unsupported.entry("InternalOrExternalShape").or_default() += 1;
            continue;
        };
        if !t.is_rigid(1e-12) {
            *unsupported.entry("NonRigidLocation").or_default() += 1;
            continue;
        }
        let shape = &doc.shapes[sub.shape];
        match shape.kind {
            Kind::Compound => {
                for s in shape.subs.iter().rev() {
                    stack.push((*s, t, o));
                }
            }
            Kind::Solid => {
                let solid = import_solid(doc, sub.shape, &t, o);
                if let Err(Rejected::Unsupported(names)) = &solid.result {
                    // Geometry is counted by record above; a solid adds its
                    // structural constructs.
                    for name in names.iter().filter(|n| !geometry.contains(*n)) {
                        *unsupported.entry(name).or_default() += 1;
                    }
                }
                solids.push(solid);
            }
            Kind::CompSolid => *unsupported.entry("CompSolid").or_default() += 1,
            Kind::Shell => *unsupported.entry("FreeShell").or_default() += 1,
            Kind::Face => *unsupported.entry("FreeFace").or_default() += 1,
            Kind::Wire => *unsupported.entry("FreeWire").or_default() += 1,
            Kind::Edge => *unsupported.entry("FreeEdge").or_default() += 1,
            Kind::Vertex => *unsupported.entry("FreeVertex").or_default() += 1,
        }
    }
    Import {
        solids,
        unsupported,
    }
}
