//! `.brep` solids into the cell model, by the rule of
//! `rust/tools/cell_reference.py::to_cell`: a seam pair on a cylinder or cone
//! whose uses run opposite ways, lie one period apart and continue their
//! neighbours in UV merges into periodic loops with winding numbers; the
//! vertices only seams and closed curves use disappear, and closed curves
//! that lose their vertex become ring edges. On a cone or sphere, a run that
//! is one degenerated edge (OCCT's pole) becomes a vertex loop of its vertex
//! (S3 of REVIEW_NOTES.md), a sphere bounded only by poles is the whole
//! sphere without loops, and a degenerated edge among other uses is dropped:
//! the loop passes through the pole. A degenerated edge on another surface
//! is unsupported. Each solid becomes a solid
//! region whose outer shell comes first, each shell's opposite sides a void
//! twin. The result must pass `Topology::from_parts`.
use super::read::{self, Data, Document, EdgeRep, Kind, Orient, Sub};
use super::Transform;
use crate::topology::{
    Curve2, Curve3, Edge, EdgeId, Enclosure, Face, FaceId, Fin, FinId, Issue, Loop, LoopId,
    Orientation, Region, RegionId, RegionKind, Shell, ShellId, Side, Surface, Topology,
    TopologyParts, Vertex, VertexId,
};
use crate::topology::{SplineDomain, SplineSpan};
use crate::{
    BSplineCurve2, BSplineCurve3, BSplineSurface3, Error, Frame3, KnotVector, Point2, Point3,
    Tolerance, Vec3,
};
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

/// A shell, face, wire, edge or vertex reached from the root outside any
/// solid (S6 of REVIEW_NOTES.md): a closed shell becomes a shell bounding a
/// void region, an open shell or a face a sheet (both sides of each face in
/// the infinite void), a wire or an edge a wire body and a vertex an acorn.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportedFree {
    /// The shape record (0-based, file order).
    pub record: usize,
    /// Shell, Face, Wire, Edge or Vertex.
    pub kind: Kind,
    /// The largest vertex, edge or face tolerance of the shape.
    pub tolerance: Tolerance,
    pub result: Result<Topology, Rejected>,
}

/// Everything a document holds, as far as the kernel represents it.
#[derive(Debug, Clone, PartialEq)]
pub struct Import {
    pub solids: Vec<ImportedSolid>,
    /// Free shapes, in the root's order (S6).
    pub free: Vec<ImportedFree>,
    /// Every construct the kernel cannot represent, by name, with the number
    /// of records (geometry) or shape uses (topology) that carry it.
    pub unsupported: BTreeMap<&'static str, usize>,
}

/// Two orientations composed; `None` when either is internal or external,
/// which the cell model cannot represent.
fn compose(a: Orient, b: Orient) -> Option<Orient> {
    match (a, b) {
        (Orient::Forward, Orient::Forward) | (Orient::Reversed, Orient::Reversed) => {
            Some(Orient::Forward)
        }
        (Orient::Forward, Orient::Reversed) | (Orient::Reversed, Orient::Forward) => {
            Some(Orient::Reversed)
        }
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
    /// `None` for a degenerated edge.
    curve: Option<Curve3>,
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
        if *degenerated {
            // A point in space: its pcurves carry it. Only a cone's pole may
            // use one (checked with its face).
            if start != end {
                return self.no("DegeneratedEdge");
            }
            let e = self.edges.len();
            self.edges.push(SEdge {
                start,
                end,
                curve: None,
                tolerance: *tolerance,
            });
            self.edge_keys.insert(record, t, e);
            return Some(e);
        }
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
            read::Curve3::BSpline(b) => {
                let poles = b.poles.iter().map(|p| p3(ct.point(*p))).collect();
                let build = if b.periodic {
                    BSplineCurve3::new_periodic
                } else {
                    BSplineCurve3::new
                };
                let curve = match build(
                    b.degree,
                    poles,
                    b.weights.clone(),
                    b.knots.clone(),
                    b.multiplicities.clone(),
                ) {
                    Ok(curve) => curve,
                    Err(Error::LimitExceeded(_)) => return self.no("BSplineControlDataLimit"),
                    Err(_) => return self.no("InvalidBSplineCurve"),
                };
                match spline_span(curve, f, l) {
                    Some(span) => Curve3::BSpline(span),
                    None => return self.no("BSplineRangeOutsideDomain"),
                }
            }
            read::Curve3::Other(name) => return self.no(name),
        };
        let e = self.edges.len();
        self.edges.push(SEdge {
            start,
            end,
            curve: Some(curve),
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
            return match (surface, &self.edges[edge].curve) {
                (Surface::Plane(plane), Some(curve)) => Some(on_plane(plane, curve)),
                _ => self.no("EdgeWithoutPCurve"),
            };
        };
        // A closed circle's sweep was snapped to one turn; its pcurves share
        // the edge's range (SameRange), so theirs is snapped alike. So is a
        // degenerated edge's turn around the apex.
        let turn = match &self.edges[edge].curve {
            Some(Curve3::CircularArc { sweep_angle, .. }) => *sweep_angle == TAU,
            None => true,
            _ => false,
        };
        if turn && ((l - f) - TAU).abs() <= 1e-12 * TAU {
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
            read::Curve2::BSpline(b) => {
                let poles = b.poles.iter().map(|p| Point2::new(p[0], p[1])).collect();
                let build = if b.periodic {
                    BSplineCurve2::new_periodic
                } else {
                    BSplineCurve2::new
                };
                let curve = match build(
                    b.degree,
                    poles,
                    b.weights.clone(),
                    b.knots.clone(),
                    b.multiplicities.clone(),
                ) {
                    Ok(curve) => curve,
                    Err(Error::LimitExceeded(_)) => return self.no("BSplineControlDataLimit"),
                    Err(_) => return self.no("InvalidBSplineCurve2d"),
                };
                match spline_span(curve, f, l) {
                    Some(span) => Curve2::BSpline(span),
                    None => return self.no("BSplineRangeOutsideDomain"),
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
            read::Surface::Cone { p, n, x, y, r, a } => {
                // An indirect axis is the kernel's cone about -N with v and
                // the semi-angle negated: the same points, the normal inward.
                indirect = dot(cross(*x, *y), *n) <= 0.0;
                let axis = if indirect { n.map(|c| -c) } else { *n };
                Surface::Cone {
                    frame: Frame3::new(
                        p3(st.point(*p)),
                        v3(st.vector(axis)),
                        v3(st.vector(*x)),
                        self.placement,
                    )
                    .ok()?,
                    radius: *r,
                    half_angle: if indirect { -*a } else { *a },
                }
            }
            read::Surface::Sphere { p, n, x, y, r } => {
                // An indirect axis is the kernel's sphere about -N with v
                // negated, as for the cylinder.
                indirect = dot(cross(*x, *y), *n) <= 0.0;
                let axis = if indirect { n.map(|c| -c) } else { *n };
                Surface::Sphere {
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
            read::Surface::Torus {
                p,
                n,
                x,
                y,
                major,
                minor,
            } => {
                // An indirect axis is the kernel's torus about -N with v
                // negated, as for the cylinder.
                indirect = dot(cross(*x, *y), *n) <= 0.0;
                let axis = if indirect { n.map(|c| -c) } else { *n };
                Surface::Torus {
                    frame: Frame3::new(
                        p3(st.point(*p)),
                        v3(st.vector(axis)),
                        v3(st.vector(*x)),
                        self.placement,
                    )
                    .ok()?,
                    major: *major,
                    minor: *minor,
                }
            }
            read::Surface::BSpline(b) => {
                // A periodic spline surface's loops would wind with its
                // domain's period, which the model does not have yet.
                if b.periodic[0] || b.periodic[1] {
                    return self.no("PeriodicBSplineSurface");
                }
                let axis = |k: usize| {
                    KnotVector::new(
                        b.degrees[k],
                        b.knots[k].clone(),
                        b.multiplicities[k].clone(),
                    )
                };
                // The kernel's resource limit (MAX_POLES) is its own name.
                let limit = |e: &Error| matches!(e, Error::LimitExceeded(_));
                let (u, v) = match (axis(0), axis(1)) {
                    (Ok(u), Ok(v)) => (u, v),
                    (Err(e), _) | (_, Err(e)) if limit(&e) => {
                        return self.no("BSplineControlDataLimit")
                    }
                    _ => return self.no("InvalidBSplineSurface"),
                };
                let poles = b.poles.iter().map(|p| p3(st.point(*p))).collect();
                match BSplineSurface3::new(u, v, poles, b.weights.clone()) {
                    Ok(s) => Surface::BSpline(s),
                    Err(e) if limit(&e) => return self.no("BSplineControlDataLimit"),
                    Err(_) => return self.no("InvalidBSplineSurface"),
                }
            }
            read::Surface::Other(name) => return self.no(name),
        };
        let cone = matches!(surface, Surface::Cone { .. } | Surface::Sphere { .. });
        let mut loops = Vec::new();
        for wire in &shape.subs {
            if doc.shapes[wire.shape].kind != Kind::Wire {
                return self.no("FaceWithNonWireChild");
            }
            if compose(wire.orient, Orient::Forward).is_none() {
                return self.no("InternalOrExternalWire");
            }
            let wt = t.times(&location(doc, wire.location));
            let mut uses = Vec::new();
            for e in &doc.shapes[wire.shape].subs {
                let et = wt.times(&location(doc, e.location));
                let stored =
                    compose(e.orient, wire.orient).or_else(|| self.no("InternalOrExternalEdge"))?;
                let traversal = compose(stored, oriented).expect("a face is forward or reversed");
                let edge = self.edge(e.shape, &et)?;
                if self.edges[edge].curve.is_none() && !cone {
                    return self.no("DegeneratedEdge");
                }
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
        // On a plane or a (non-periodic) spline surface the outer loop comes
        // first; OCCT stores wires in any order. The outer loop encloses the
        // largest area in UV.
        if matches!(surface, Surface::Plane(_) | Surface::BSpline(_)) && loops.len() > 1 {
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
        // The poles in the plane's coordinates: the curve lies in the plane.
        Curve3::BSpline(span) => {
            let c = span.curve();
            let poles = c.poles().iter().map(|p| uv(*p)).collect();
            let build = if c.is_periodic() {
                BSplineCurve2::new_periodic
            } else {
                BSplineCurve2::new
            };
            let curve = build(
                c.degree(),
                poles,
                Some(c.weights().to_vec()),
                c.knots().to_vec(),
                c.multiplicities().to_vec(),
            )
            .expect("the same basis and weights");
            let [first, last] = span.range();
            Curve2::BSpline(SplineSpan::new(curve, first, last).expect("the same range"))
        }
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
        // A hyperbola or parabola on a plane: its exact projection (S8d.2).
        Curve3::HyperbolaArc { .. } | Curve3::ParabolaArc { .. } | Curve3::Section(_) => {
            crate::topology::plane_pcurve(curve, crate::topology::Orientation::Forward, *plane)
        }
        // Axis-aligned with the plane (the kernel's own ellipses, whose
        // frames share the plane's axes); validation rejects any other.
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => {
            let turn = if frame.normal().dot(plane.normal()) > 0.0 {
                1.0
            } else {
                -1.0
            };
            Curve2::EllipseArc {
                center: uv(frame.origin()),
                major: *major,
                minor: *minor,
                start_angle: turn * start_angle,
                sweep_angle: turn * sweep_angle,
            }
        }
    }
}

fn negate_v(p: &Curve2) -> Curve2 {
    match p {
        Curve2::BSpline(span) => {
            let c = span.curve();
            let c3 = c.as_curve3();
            let poles = c.poles().iter().map(|p| Point2::new(p.x, -p.y)).collect();
            let build = if c3.is_periodic() {
                BSplineCurve2::new_periodic
            } else {
                BSplineCurve2::new
            };
            let curve = build(
                c3.degree(),
                poles,
                Some(c3.weights().to_vec()),
                c3.knots().to_vec(),
                c3.multiplicities().to_vec(),
            )
            .expect("the same basis and weights");
            let [first, last] = span.range();
            let out = SplineSpan::new(curve, first, last).expect("the same range");
            Curve2::BSpline(if span.is_reversed() {
                out.reversed()
            } else {
                out
            })
        }
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
        Curve2::EllipseArc {
            center,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => Curve2::EllipseArc {
            center: Point2::new(center.x, -center.y),
            major: *major,
            minor: *minor,
            start_angle: -start_angle,
            sweep_angle: -sweep_angle,
        },
        Curve2::Sinusoid { start, sweep, a } => Curve2::Sinusoid {
            start: *start,
            sweep: *sweep,
            a: a.map(|x| -x),
        },
        Curve2::Projection(_) => unreachable!("the reader makes no projection pcurves"),
    }
}

fn location_of(doc: &Document, index: usize) -> Transform {
    location(doc, index)
}

fn reversed(p: &Curve2) -> Curve2 {
    match p {
        Curve2::BSpline(span) => Curve2::BSpline(span.reversed()),
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
        Curve2::EllipseArc {
            center,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => Curve2::EllipseArc {
            center: *center,
            major: *major,
            minor: *minor,
            start_angle: start_angle + sweep_angle,
            sweep_angle: -sweep_angle,
        },
        Curve2::Sinusoid { start, sweep, a } => Curve2::Sinusoid {
            start: start + sweep,
            sweep: -sweep,
            a: *a,
        },
        Curve2::Projection(p) => {
            let mut p = (**p).clone();
            p.reversed = !p.reversed;
            p.lifts.reverse();
            Curve2::Projection(Box::new(p))
        }
    }
}

/// A spline over an edge's range: the whole domain when the range is it, a
/// range end within printing precision of a domain end snapped to it, and a
/// periodic range printed a little longer than one period shortened to one,
/// as a circle's sweep is snapped to a turn. `None` when the range leaves
/// the domain.
fn spline_span<C: SplineDomain + Clone>(
    curve: C,
    mut first: f64,
    mut last: f64,
) -> Option<SplineSpan<C>> {
    let (a, b) = curve.domain();
    // OCCT prints ranges with fifteen significant digits and knots with
    // seventeen: a range end that rounds a domain end is that end.
    let near = |x: f64, y: f64| (x - y).abs() <= 1e-12 * (1.0 + (b - a).abs() + y.abs());
    if near(first, a) {
        first = a;
    }
    if near(last, b) {
        last = b;
    }
    if (first, last) == (a, b) {
        return Some(SplineSpan::whole(curve));
    }
    if curve.is_periodic() {
        let period = b - a;
        if last - first > period && last - first <= period * (1.0 + 1e-12) {
            last = first + period;
            // Down to at most one period exactly.
            while SplineSpan::new(curve.clone(), first, last).is_err() && last > first {
                last = next_below(last);
            }
        }
    }
    SplineSpan::new(curve, first, last).ok()
}

/// The largest binary64 below a finite `x`.
fn next_below(x: f64) -> f64 {
    if x > 0.0 {
        f64::from_bits(x.to_bits() - 1)
    } else if x < 0.0 {
        f64::from_bits(x.to_bits() + 1)
    } else {
        -f64::from_bits(1)
    }
}

fn start_of(p: &Curve2) -> Point2 {
    p.point(0.0)
}

fn end_of(p: &Curve2) -> Point2 {
    p.point(1.0)
}

/// A seamed loop split at its seams: the seam edges, and the runs of use
/// indices between them with their windings in u and v.
type SeamRuns = (Vec<usize>, Vec<(Vec<usize>, [i32; 2])>);

/// The length of a unit step in `u` at height `v`: the cylinder's radius,
/// the cone's `|R + v sin a|`, the sphere's `|R cos v|`, the torus's
/// `R + r cos v`.
fn u_scale(surface: &Surface, v: f64) -> Option<f64> {
    match surface {
        Surface::BSpline(_) => None,
        Surface::Cylinder { radius, .. } => Some(*radius),
        Surface::Cone {
            radius, half_angle, ..
        } => Some((radius + v * half_angle.sin()).abs()),
        Surface::Sphere { radius, .. } => Some((radius * v.cos()).abs()),
        Surface::Torus { major, minor, .. } => Some(major + minor * v.cos()),
        Surface::Plane(_) => None,
    }
}

/// The length of a unit step in `v`: a sphere's radius, a torus's minor
/// radius, 1 along rulings.
fn v_scale(surface: &Surface) -> f64 {
    match surface {
        Surface::Sphere { radius, .. } => *radius,
        Surface::Torus { minor, .. } => *minor,
        _ => 1.0,
    }
}

/// Split one seamed loop into runs between seam uses, each with its
/// windings; `None` when the loop has no consistent seam pair. A seam pair's
/// uses lie one period apart in u, or on a torus in v. A torus loop of seams
/// only, in both directions, is the whole torus (no runs).
fn seam_merge(face: &SFace, lp: &[SUse], tol: f64) -> Option<SeamRuns> {
    let surface = &face.surface;
    u_scale(surface, 0.0)?;
    let torus = matches!(surface, Surface::Torus { .. });
    let vs = v_scale(surface);
    let scale = |a: Point2, b: Point2| {
        u_scale(surface, a.y)
            .unwrap_or(0.0)
            .max(u_scale(surface, b.y).unwrap_or(0.0))
    };
    let close = |a: Point2, b: Point2| {
        ((a.x - b.x) * scale(a, b)).abs() <= tol && ((a.y - b.y) * vs).abs() <= tol
    };
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
    let (mut in_u, mut in_v) = (false, false);
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
        let shift = (sa.x - eb.x, sa.y - eb.y);
        let radius = scale(*sa, *ea);
        let same = ((ea.x - sb.x) - shift.0).abs() * radius <= tol
            && ((ea.y - sb.y) - shift.1).abs() * vs <= tol;
        let period_u = (shift.0.abs() - TAU).abs() * radius <= tol && (shift.1 * vs).abs() <= tol;
        let period_v =
            torus && (shift.0 * radius).abs() <= tol && (shift.1.abs() - TAU).abs() * vs <= tol;
        if !same || !(period_u || period_v) {
            return None;
        }
        in_u |= period_u;
        in_v |= period_v;
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
        return (torus && in_u && in_v).then_some((seams, Vec::new()));
    }
    let mut out = Vec::new();
    for run in runs {
        let (a, b) = (
            start_of(&lp[run[0]].pcurve),
            end_of(&lp[run[run.len() - 1]].pcurve),
        );
        let (du, dv) = (b.x - a.x, b.y - a.y);
        let (wu, wv) = ((du / TAU).round(), (dv / TAU).round());
        if (wu == 0.0 && wv == 0.0)
            || (du - wu * TAU).abs() * scale(a, b) > tol
            || (dv - wv * TAU).abs() * vs > tol
        {
            return None;
        }
        out.push((run, [wu as i32, wv as i32]));
    }
    Some((seams, out))
}

/// The uses of a run, a degenerated edge among others dropped: the loop
/// passes through the pole at its vertex, where the UV gap in `u` has no
/// length (`rho = 0`) and the closing chord takes the degenerated edge's
/// place. The vertex stays, used by the neighbours.
fn through_poles<'a>(
    walk: &Walk,
    uses: impl Iterator<Item = &'a SUse>,
    removed: &mut BTreeSet<usize>,
) -> Vec<&'a SUse> {
    let uses: Vec<&SUse> = uses.collect();
    if uses.len() < 2 {
        return uses;
    }
    uses.into_iter()
        .filter(|u| {
            let degenerated = walk.edges[u.edge].curve.is_none();
            if degenerated {
                removed.insert(u.edge);
            }
            !degenerated
        })
        .collect()
}

/// A loop of the cell complex before numbering: runs of seamed uses with
/// their winding, or a pole at a seamed vertex.
enum CellLoop<'a> {
    Edges(Vec<&'a SUse>, [i32; 2]),
    Pole(usize),
}

/// The cell complex of a seamed solid, by the rule of `to_cell`; the name of
/// the unsupported construct when a degenerated edge is not a cone's pole.
/// What the walked shapes become (S6): a solid's shells, a closed shell
/// without a solid, an open sheet, a wire of the listed edges, an acorn.
enum Mode {
    Solid,
    ClosedShell,
    Sheet,
    Wire(Vec<usize>),
    Acorn(usize),
}

fn to_cell(walk: Walk, tol: f64, mode: Mode) -> Result<TopologyParts, &'static str> {
    let mut face_loops: Vec<Vec<CellLoop>> = Vec::new();
    let mut removed = BTreeSet::new();
    let mut poles = BTreeSet::new();
    for f in &walk.faces {
        let mut loops = Vec::new();
        // A closed spline surface's seam (an edge used twice by the face)
        // would need windings with its domain's period (S4e of
        // REVIEW_NOTES.md).
        if matches!(f.surface, Surface::BSpline(_)) {
            let mut used = BTreeSet::new();
            if f.loops.iter().flatten().any(|u| !used.insert(u.edge)) {
                return Err("SeamOnBSplineSurface");
            }
        }
        for lp in &f.loops {
            match (!lp.is_empty()).then(|| seam_merge(f, lp, tol)).flatten() {
                Some((seams, runs)) => {
                    removed.extend(seams);
                    for (run, w) in runs {
                        let e = lp[run[0]].edge;
                        if run.len() == 1 && walk.edges[e].curve.is_none() {
                            removed.insert(e);
                            loops.push(CellLoop::Pole(walk.edges[e].start));
                        } else {
                            loops.push(CellLoop::Edges(
                                through_poles(&walk, run.into_iter().map(|k| &lp[k]), &mut removed),
                                w,
                            ));
                        }
                    }
                }
                None => loops.push(CellLoop::Edges(
                    through_poles(&walk, lp.iter(), &mut removed),
                    [0, 0],
                )),
            }
        }
        // A sphere bounded only by its poles is the whole sphere: no loops,
        // and its pole vertices go with their degenerated edges. Any other
        // pole is a vertex loop.
        let whole = matches!(f.surface, Surface::Sphere { .. })
            && loops.iter().all(|l| matches!(l, CellLoop::Pole(_)));
        if whole {
            loops.clear();
        }
        for l in &loops {
            if let CellLoop::Pole(v) = l {
                poles.insert(*v);
            }
        }
        face_loops.push(loops);
    }
    let still_used: BTreeSet<usize> = face_loops
        .iter()
        .flatten()
        .flat_map(|l| match l {
            CellLoop::Edges(lp, _) => lp.iter().map(|u| u.edge).collect(),
            CellLoop::Pole(_) => Vec::new(),
        })
        .collect();
    if still_used.iter().any(|e| walk.edges[*e].curve.is_none()) {
        return Err("DegeneratedEdge");
    }
    if let Mode::Wire(edges) = &mode {
        if edges.iter().any(|e| walk.edges[*e].curve.is_none()) {
            return Err("DegeneratedEdge");
        }
    }
    removed.retain(|e| !still_used.contains(e));
    let closed = |c: &Option<Curve3>| matches!(c, Some(Curve3::CircularArc { sweep_angle, .. }) if sweep_angle.abs() == TAU);
    let mut other_use = poles;
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
            curve: e.curve.clone().expect("degenerated edges are poles"),
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
        for l in &face_loops[fi] {
            let (lp, w) = match l {
                CellLoop::Edges(lp, w) => (lp, w),
                CellLoop::Pole(v) => {
                    parts.loops.push(Loop::Vertex(vertex_map[v]));
                    loops.push(LoopId::new(parts.loops.len() - 1));
                    continue;
                }
            };
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
            parts.loops.push(Loop::Edges { fins, winding: *w });
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
            front: ShellId::new(if matches!(mode, Mode::Sheet) { 0 } else { k }),
            back: ShellId::new(if matches!(mode, Mode::Sheet) {
                0
            } else {
                n + k
            }),
            // Computed after conversion: OCCT stores no bound on UV closure.
            enclosure: None,
        });
    }
    for (k, fin) in parts.fins.iter().enumerate() {
        parts.edges[fin.edge.index()].fins.push(FinId::new(k));
    }
    let void = |shell: Shell| {
        (
            vec![shell],
            vec![Region {
                kind: RegionKind::Void,
                shells: vec![ShellId::new(0)],
            }],
        )
    };
    let (shells, regions) = match &mode {
        Mode::Sheet => void(Shell {
            region: RegionId::new(0),
            sides: (0..parts.faces.len())
                .flat_map(|f| [(FaceId::new(f), Side::Front), (FaceId::new(f), Side::Back)])
                .collect(),
            wire_edges: Vec::new(),
            acorn_vertices: Vec::new(),
        }),
        Mode::Wire(edges) => void(Shell {
            region: RegionId::new(0),
            sides: Vec::new(),
            wire_edges: edges.iter().map(|e| edge_map[e]).collect(),
            acorn_vertices: Vec::new(),
        }),
        Mode::Acorn(v) => void(Shell {
            region: RegionId::new(0),
            sides: Vec::new(),
            wire_edges: Vec::new(),
            acorn_vertices: vec![vertex_map[v]],
        }),
        Mode::Solid | Mode::ClosedShell => (Vec::new(), Vec::new()),
    };
    if !regions.is_empty() {
        parts.shells = shells;
        parts.regions = regions;
        return Ok(parts);
    }
    parts.regions = vec![
        Region {
            kind: RegionKind::Void,
            shells: vec![ShellId::new(n)],
        },
        Region {
            // A closed shell without a solid bounds a void region (S6).
            kind: if matches!(mode, Mode::ClosedShell) {
                RegionKind::Void
            } else {
                RegionKind::Solid
            },
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
    Ok(parts)
}

fn import_solid(doc: &Document, record: usize, t: &Transform, oriented: Orient) -> ImportedSolid {
    let mut walk = new_walk(doc);
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
    let parts = match to_cell(walk, resolution.linear(), Mode::Solid) {
        Ok(parts) => parts.with_measured_enclosures(),
        Err(name) => {
            return ImportedSolid {
                record,
                tolerance: resolution,
                result: Err(Rejected::Unsupported(vec![name])),
            }
        }
    };
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

fn new_walk(doc: &Document) -> Walk<'_> {
    Walk {
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
    }
}

/// A free shell, face, wire, edge or vertex as a body (S6).
fn import_free(doc: &Document, record: usize, t: &Transform, oriented: Orient) -> ImportedFree {
    let mut walk = new_walk(doc);
    let kind = doc.shapes[record].kind;
    let mode = match kind {
        Kind::Shell => {
            let mut faces = Vec::new();
            for face in &doc.shapes[record].subs {
                let Some(fo) = compose(face.orient, oriented) else {
                    walk.unsupported.push("InternalOrExternalFace");
                    continue;
                };
                if doc.shapes[face.shape].kind != Kind::Face {
                    walk.unsupported.push("ShellWithNonFaceChild");
                    continue;
                }
                let ft = t.times(&location(doc, face.location));
                if let Some(f) = walk.face(face.shape, &ft, fo) {
                    faces.push(f);
                }
            }
            // Closed when every edge with a curve is used exactly twice
            // (a seam's two uses in one face included).
            let mut uses: BTreeMap<usize, usize> = BTreeMap::new();
            for u in faces
                .iter()
                .flat_map(|f| walk.faces[*f].loops.iter().flatten())
            {
                if walk.edges[u.edge].curve.is_some() {
                    *uses.entry(u.edge).or_default() += 1;
                }
            }
            let closed = !uses.is_empty() && uses.values().all(|n| *n == 2);
            walk.shells.push(faces);
            if closed {
                Mode::ClosedShell
            } else {
                Mode::Sheet
            }
        }
        Kind::Face => {
            if let Some(f) = walk.face(record, t, oriented) {
                walk.shells.push(vec![f]);
            }
            Mode::Sheet
        }
        Kind::Wire => {
            let mut edges = Vec::new();
            for e in &doc.shapes[record].subs {
                if compose(e.orient, Orient::Forward).is_none() {
                    walk.unsupported.push("InternalOrExternalEdge");
                    continue;
                }
                if doc.shapes[e.shape].kind != Kind::Edge {
                    walk.unsupported.push("WireWithNonEdgeChild");
                    continue;
                }
                let et = t.times(&location(doc, e.location));
                if let Some(edge) = walk.edge(e.shape, &et) {
                    if !edges.contains(&edge) {
                        edges.push(edge);
                    }
                }
            }
            Mode::Wire(edges)
        }
        Kind::Edge => match walk.edge(record, t) {
            Some(edge) => Mode::Wire(vec![edge]),
            None => Mode::Wire(Vec::new()),
        },
        _ => match walk.vertex(record, t) {
            Some(v) => Mode::Acorn(v),
            None => Mode::Acorn(0),
        },
    };
    let resolution = Tolerance::new(walk.tolerance.max(Tolerance::default().linear()), 1e-12)
        .unwrap_or_default();
    let rejected = |names: Vec<&'static str>| ImportedFree {
        record,
        kind,
        tolerance: resolution,
        result: Err(Rejected::Unsupported(names)),
    };
    if !walk.unsupported.is_empty() {
        let mut names = std::mem::take(&mut walk.unsupported);
        names.sort_unstable();
        names.dedup();
        return rejected(names);
    }
    let parts = match to_cell(walk, resolution.linear(), mode) {
        Ok(parts) => parts.with_measured_enclosures(),
        Err(name) => return rejected(vec![name]),
    };
    ImportedFree {
        record,
        kind,
        tolerance: resolution,
        result: Topology::from_parts(parts.clone(), resolution).map_err(|issues| {
            Rejected::Invalid {
                issues,
                parts: Box::new(parts),
            }
        }),
    }
}

/// Every solid and free shape of a document in the cell model, and every
/// construct the kernel cannot represent, by name.
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
    let mut free = Vec::new();
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
            Kind::Shell | Kind::Face | Kind::Wire | Kind::Edge | Kind::Vertex => {
                let body = import_free(doc, sub.shape, &t, o);
                if let Err(Rejected::Unsupported(names)) = &body.result {
                    for name in names.iter().filter(|n| !geometry.contains(*n)) {
                        *unsupported.entry(name).or_default() += 1;
                    }
                }
                free.push(body);
            }
        }
    }
    Import {
        solids,
        free,
        unsupported,
    }
}
