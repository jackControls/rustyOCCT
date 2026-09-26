//! The writer: a cell topology's solid regions as `.brep` version 1 text,
//! with OCCT's structure inserted by rule. A face that winds around a
//! cylinder gets one seam at the `u` where its ring loops start: its forward
//! use (in the unoriented face) at `u0 + 2 pi` going from the lower ring to
//! the upper one, its reversed use at `u0`; every ring edge gets a vertex
//! where the seam meets it and becomes a closed edge. Every pcurve shares its
//! edge's parameter, as OCCT's SameParameter edges do.
use super::BrepError;
use crate::topology::{
    Curve2, Curve3, EdgeId, FaceId, Loop, Orientation, RegionKind, Side, Surface, Topology,
};
use crate::{Point2, Point3};
use std::collections::BTreeMap;
use std::f64::consts::TAU;
use std::fmt::Write as _;

fn num(x: f64) -> String {
    // Shortest round-trip text; OCCT reads it with operator>>.
    let s = format!("{x:?}");
    s.strip_suffix(".0").map(str::to_string).unwrap_or(s)
}

fn nums(v: &[f64]) -> String {
    v.iter().map(|x| num(*x)).collect::<Vec<_>>().join(" ")
}

fn unwritable(what: &'static str) -> BrepError {
    BrepError::Unwritable { what }
}

/// A 3D curve record with the edge's parameter range.
struct EdgeGeometry {
    record: String,
    range: [f64; 2],
}

/// The 3D record of a curve and its range, the parameter increasing along
/// the edge.
fn curve_record(c: &Curve3, ring_start: Option<f64>) -> EdgeGeometry {
    match c {
        Curve3::LineSegment { start, end } => {
            let d = *end - *start;
            let l = d.length();
            let u = d * (1.0 / l);
            EdgeGeometry {
                record: format!("1 {} {}", nums(&start.to_array()), nums(&u.to_array())),
                range: [0.0, l],
            }
        }
        Curve3::Circle { frame, radius } => curve_record(
            &Curve3::CircularArc {
                frame: *frame,
                radius: *radius,
                start_angle: 0.0,
                sweep_angle: TAU,
            },
            ring_start,
        ),
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => {
            // A negative sweep is the positive sweep about the flipped axis:
            // t = -angle.
            let flip = *sweep_angle < 0.0;
            let (n, y) = if flip {
                (-frame.normal(), -frame.normal().cross(frame.x()))
            } else {
                (frame.normal(), frame.normal().cross(frame.x()))
            };
            let start = match ring_start {
                Some(s) => s,
                None if flip => -start_angle,
                None => *start_angle,
            };
            EdgeGeometry {
                record: format!(
                    "2 {} {} {} {} {}",
                    nums(&frame.origin().to_array()),
                    nums(&n.to_array()),
                    nums(&frame.x().to_array()),
                    nums(&y.to_array()),
                    num(*radius)
                ),
                range: [start, start + sweep_angle.abs()],
            }
        }
    }
}

/// The point of a curve at its own parameter `t` (as `curve_record` writes it).
fn curve_at(c: &Curve3, t: f64) -> Point3 {
    match c {
        Curve3::LineSegment { start, end } => {
            let d = *end - *start;
            *start + d * (t / d.length())
        }
        Curve3::Circle { frame, radius } => {
            frame.point(Point2::new(radius * t.cos(), radius * t.sin()), 0.0)
        }
        Curve3::CircularArc {
            frame,
            radius,
            sweep_angle,
            ..
        } => {
            let a = if *sweep_angle < 0.0 { -t } else { t };
            frame.point(Point2::new(radius * a.cos(), radius * a.sin()), 0.0)
        }
    }
}

/// The 2D record of a fin's pcurve in the edge's parameter `[t0, t1]`.
fn pcurve_record(
    surface: &Surface,
    pcurve: &Curve2,
    forward: bool,
    range: [f64; 2],
) -> Result<String, BrepError> {
    let [t0, t1] = range;
    Ok(match (surface, pcurve) {
        (_, Curve2::LineSegment { start, end }) => {
            // In the edge's direction.
            let (a, b) = if forward {
                (*start, *end)
            } else {
                (*end, *start)
            };
            let (du, dv) = (b.x - a.x, b.y - a.y);
            let span = t1 - t0;
            let (dx, dy) = match surface {
                // On a cylinder a horizontal pcurve runs in u at unit speed
                // against an angle parameter.
                Surface::Cylinder { .. } if dv == 0.0 => (du.signum(), 0.0),
                _ => {
                    let l = du.hypot(dv);
                    (du / l, dv / l)
                }
            };
            let _ = span;
            format!(
                "1 {} {}",
                nums(&[a.x - dx * t0, a.y - dy * t0]),
                nums(&[dx, dy])
            )
        }
        (
            Surface::Plane(_),
            Curve2::CircularArc {
                center,
                radius,
                start_angle,
                sweep_angle,
            },
        ) => {
            let (start, sweep) = if forward {
                (*start_angle, *sweep_angle)
            } else {
                (start_angle + sweep_angle, -sweep_angle)
            };
            let sigma = sweep.signum();
            let alpha = start - sigma * t0;
            let (x, y) = (
                (alpha.cos(), alpha.sin()),
                (-sigma * alpha.sin(), sigma * alpha.cos()),
            );
            format!(
                "2 {} {} {} {}",
                nums(&[center.x, center.y]),
                nums(&[x.0, x.1]),
                nums(&[y.0, y.1]),
                num(*radius)
            )
        }
        (Surface::Cylinder { .. } | Surface::Cone { .. }, Curve2::CircularArc { .. }) => {
            return Err(unwritable("an arc pcurve on a cylinder or cone"))
        }
    })
}

/// Where a seam meets a wound loop: an existing vertex, or the seam vertex
/// of a ring edge.
#[derive(Clone, Copy)]
enum Vertex {
    Shared(usize),
    Ring(usize),
}

#[derive(Clone, Copy)]
struct End {
    vertex: Vertex,
    /// `v` of the seam end.
    height: f64,
}

/// A wound loop: its winding in the unoriented face, its ring edge if it is
/// one, and otherwise each vertex with its position in the face.
type WoundLoop = (i32, Option<usize>, Vec<(usize, Point2)>);

/// A wound face's seam at `u0`, from the loop winding `+u` (in the
/// unoriented face) to the loop winding `-u`.
struct Seam {
    u0: f64,
    bottom: End,
    top: End,
}

struct Records {
    /// Shape records in file order.
    shapes: Vec<String>,
}

impl Records {
    fn push(&mut self, record: String) -> usize {
        self.shapes.push(record);
        self.shapes.len() - 1
    }
}

/// Version 1 `.brep` text of every solid region of `topology`; one solid is
/// the root, several a compound.
pub fn write(topology: &Topology, tolerance: f64) -> Result<String, BrepError> {
    let t = topology;
    let tol = num(tolerance);
    let mut curves2d: Vec<String> = Vec::new();
    let mut curves: Vec<String> = Vec::new();
    let mut surfaces: Vec<String> = Vec::new();
    let mut records = Records { shapes: Vec::new() };
    // Seams: one per wound face, at a u0 where both wound loops have a
    // vertex (a ring loop gets its seam vertex there).
    let mut ring_start: BTreeMap<usize, f64> = BTreeMap::new();
    let mut wound: BTreeMap<usize, Seam> = BTreeMap::new();
    for (fi, face) in t.faces().iter().enumerate() {
        let Surface::Cylinder { frame, radius } = &face.surface else {
            continue;
        };
        let flip = face.sense == Orientation::Reversed;
        let mut loops: Vec<WoundLoop> = Vec::new();
        for l in &face.loops {
            let Loop::Edges { fins, winding } = &t.loops()[l.index()] else {
                continue;
            };
            if winding[1] != 0 {
                return Err(unwritable("a winding in v"));
            }
            if winding[0] == 0 {
                continue;
            }
            if winding[0].abs() != 1 {
                return Err(unwritable("a loop winding more than once"));
            }
            let first = &t.fins()[fins[0].index()];
            let ring = (fins.len() == 1 && t.edges()[first.edge.index()].is_ring())
                .then_some(first.edge.index());
            let mut at = Vec::new();
            for k in fins {
                let f = &t.fins()[k.index()];
                let edge = &t.edges()[f.edge.index()];
                let v = if f.sense == Orientation::Forward {
                    edge.start
                } else {
                    edge.end
                };
                if let Some(v) = v {
                    at.push((v.index(), f.pcurve.point(0.0)));
                }
            }
            if ring.is_none() && at.len() != fins.len() {
                return Err(unwritable("a wound loop mixing rings and vertices"));
            }
            loops.push((if flip { -winding[0] } else { winding[0] }, ring, at));
        }
        if loops.is_empty() {
            continue;
        }
        let (Some(bottom), Some(top), 2) = (
            loops.iter().find(|l| l.0 > 0),
            loops.iter().find(|l| l.0 < 0),
            loops.len(),
        ) else {
            return Err(unwritable("a wound face without one loop each way"));
        };
        let same = |a: f64, b: f64| {
            let d = (a - b).rem_euclid(TAU);
            d.min(TAU - d) * radius <= tolerance
        };
        // The seam position: a vertex both loops have, or any vertex of a
        // loop facing a ring, or where the rings' pcurves start.
        let ring_u = |e: usize| {
            t.face_fins(FaceId::new(fi))
                .into_iter()
                .flatten()
                .find(|f| f.edge.index() == e)
                .map(|f| f.pcurve.point(0.0).x)
        };
        let choice = match (bottom.1, top.1) {
            (Some(e), Some(_)) => ring_u(e).map(|u| (u, None, None)),
            (Some(_), None) => top.2.first().map(|(w, q)| (q.x, None, Some((*w, q.y)))),
            (None, Some(_)) => bottom.2.first().map(|(v, p)| (p.x, Some((*v, p.y)), None)),
            (None, None) => bottom.2.iter().find_map(|(v, p)| {
                top.2
                    .iter()
                    .find(|(_, q)| same(p.x, q.x))
                    .map(|(w, q)| (p.x, Some((*v, p.y)), Some((*w, q.y))))
            }),
        };
        let Some((u, bottom_vertex, top_vertex)) = choice else {
            return Err(unwritable(
                "a wound face whose loops share no seam position",
            ));
        };
        let u0 = u.rem_euclid(TAU);
        let mut end =
            |ring: Option<usize>, vertex: Option<(usize, f64)>| -> Result<End, BrepError> {
                if let Some((v, height)) = vertex {
                    return Ok(End {
                        vertex: Vertex::Shared(v),
                        height,
                    });
                }
                let e = ring.expect("a ring loop");
                let fin = t
                    .face_fins(FaceId::new(fi))
                    .into_iter()
                    .flatten()
                    .find(|f| f.edge.index() == e)
                    .expect("the ring's fin on this face");
                let height = fin.pcurve.point(0.0).y;
                // The ring's own parameter at the seam point.
                let (Curve3::Circle {
                    frame: cf,
                    radius: cr,
                }
                | Curve3::CircularArc {
                    frame: cf,
                    radius: cr,
                    ..
                }) = &t.edges()[e].curve
                else {
                    return Err(unwritable("a ring edge that is not a circle"));
                };
                if (cr - radius).abs() > tolerance {
                    return Err(unwritable("a ring edge off its wall's radius"));
                }
                let dir = frame.x() * u0.cos() + frame.normal().cross(frame.x()) * u0.sin();
                let ny = cf.normal().cross(cf.x());
                let mut angle = dir.dot(ny).atan2(dir.dot(cf.x()));
                if let Curve3::CircularArc { sweep_angle, .. } = &t.edges()[e].curve {
                    if *sweep_angle < 0.0 {
                        angle = -angle;
                    }
                }
                ring_start.entry(e).or_insert(angle);
                Ok(End {
                    vertex: Vertex::Ring(e),
                    height,
                })
            };
        let bottom_end = end(bottom.1, bottom_vertex)?;
        let top_end = end(top.1, top_vertex)?;
        wound.insert(
            fi,
            Seam {
                u0,
                bottom: bottom_end,
                top: top_end,
            },
        );
    }
    // Vertices, then seam vertices of ring edges.
    let mut vertex_record: BTreeMap<usize, usize> = BTreeMap::new();
    let vertex = |p: Point3, records: &mut Records| {
        records.push(format!(
            "Ve\n{tol}\n{}\n0 0\n\n0101101\n*",
            nums(&p.to_array())
        ))
    };
    for (v, x) in t.vertices().iter().enumerate() {
        vertex_record.insert(v, vertex(x.position, &mut records));
    }
    let mut ring_vertex: BTreeMap<usize, usize> = BTreeMap::new();
    for (e, start) in &ring_start {
        let p = curve_at(&t.edges()[*e].curve, *start);
        ring_vertex.insert(*e, vertex(p, &mut records));
    }
    // Surfaces.
    for face in t.faces() {
        if matches!(face.surface, Surface::Cone { .. }) {
            return Err(unwritable("a cone"));
        }
        surfaces.push(match &face.surface {
            Surface::Plane(f) => format!(
                "1 {} {} {} {}",
                nums(&f.origin().to_array()),
                nums(&f.normal().to_array()),
                nums(&f.x().to_array()),
                nums(&f.normal().cross(f.x()).to_array())
            ),
            Surface::Cylinder { frame: f, radius } => format!(
                "2 {} {} {} {} {}",
                nums(&f.origin().to_array()),
                nums(&f.normal().to_array()),
                nums(&f.x().to_array()),
                nums(&f.normal().cross(f.x()).to_array()),
                num(*radius)
            ),
            Surface::Cone { .. } => unreachable!("refused above"),
        });
    }
    // Edge geometry and pcurves per face.
    let geometry: Vec<EdgeGeometry> = t
        .edges()
        .iter()
        .enumerate()
        .map(|(e, edge)| curve_record(&edge.curve, ring_start.get(&e).copied()))
        .collect();
    let mut reps: Vec<Vec<String>> = vec![Vec::new(); t.edges().len()];
    for (fi, face) in t.faces().iter().enumerate() {
        for fins in t.face_fins(FaceId::new(fi)) {
            for fin in fins {
                let e = fin.edge.index();
                let range = geometry[e].range;
                let mut pcurve = fin.pcurve.clone();
                // On a wound face every pcurve lies in [u0, u0 + 2 pi].
                if let (Some(seam), Curve2::LineSegment { start, end }) = (wound.get(&fi), &pcurve)
                {
                    let low = start.x.min(end.x);
                    let k = -((low - seam.u0) / TAU + 1e-9).floor() * TAU;
                    pcurve = Curve2::LineSegment {
                        start: Point2::new(start.x + k, start.y),
                        end: Point2::new(end.x + k, end.y),
                    };
                }
                let forward = fin.sense == Orientation::Forward;
                curves2d.push(pcurve_record(&face.surface, &pcurve, forward, range)?);
                reps[e].push(format!(
                    "2 {} {} 0 {}",
                    curves2d.len(),
                    fi + 1,
                    nums(&range)
                ));
            }
        }
    }
    // Seam edges: from the bottom loop's seam vertex to the top one's.
    let seam_vertex = |end: &End| match end.vertex {
        Vertex::Shared(v) => vertex_record[&v],
        Vertex::Ring(e) => ring_vertex[&e],
    };
    let mut seam_of: BTreeMap<usize, usize> = BTreeMap::new();
    for (fi, seam) in &wound {
        let Surface::Cylinder { frame, radius } = &t.faces()[*fi].surface else {
            unreachable!("wound faces are cylinders");
        };
        let (u0, low, high) = (seam.u0, seam.bottom.height, seam.top.height);
        let at = |v: f64| frame.point(Point2::new(radius * u0.cos(), radius * u0.sin()), v);
        let (a, b) = (at(low), at(high));
        let l = b.distance(a);
        if l <= tolerance {
            return Err(unwritable("a wound face of no height at its seam"));
        }
        let d = (b - a) * (1.0 / l);
        let dv = (high - low).signum();
        curves.push(format!("1 {} {}", nums(&a.to_array()), nums(&d.to_array())));
        let c3 = curves.len();
        curves2d.push(format!("1 {} 0 {}", nums(&[u0 + TAU, low]), num(dv)));
        curves2d.push(format!("1 {} 0 {}", nums(&[u0, low]), num(dv)));
        let (p1, p2) = (curves2d.len() - 1, curves2d.len());
        let edge = records.push(format!(
            "Ed\n {tol} 1 1 0\n1 {c3} 0 0 {}\n3 {p1} {p2} CN {} 0 0 {}\n0\n\n0101000\n{{v+{}}} 0 {{v-{}}} 0 *",
            num(l),
            fi + 1,
            num(l),
            seam_vertex(&seam.bottom),
            seam_vertex(&seam.top),
        ));
        seam_of.insert(*fi, edge);
    }
    // Edges.
    let mut edge_record: BTreeMap<usize, usize> = BTreeMap::new();
    for (e, edge) in t.edges().iter().enumerate() {
        if edge.fins.is_empty() {
            return Err(unwritable("a wire edge"));
        }
        curves.push(geometry[e].record.clone());
        let mut text = format!(
            "Ed\n {tol} 1 1 0\n1 {} 0 {}\n",
            curves.len(),
            nums(&geometry[e].range)
        );
        for r in &reps[e] {
            text.push_str(r);
            text.push('\n');
        }
        text.push_str("0\n\n0101000\n");
        let (s, en) = match (edge.start, edge.end) {
            (Some(s), Some(en)) => (vertex_record[&s.index()], vertex_record[&en.index()]),
            _ => match ring_vertex.get(&e) {
                Some(v) => (*v, *v),
                None => return Err(unwritable("a ring edge on no wound face")),
            },
        };
        let _ = write!(text, "{{v+{s}}} 0 {{v-{en}}} 0 *");
        edge_record.insert(e, records.push(text));
    }
    // Wires and faces, in the unoriented face: a reversed face's loops are
    // traversed backwards.
    let mut face_record: BTreeMap<usize, usize> = BTreeMap::new();
    for (fi, face) in t.faces().iter().enumerate() {
        let flip = face.sense == Orientation::Reversed;
        let mut loops: Vec<(Vec<(usize, bool)>, i32)> = Vec::new();
        for l in &face.loops {
            let Loop::Edges { fins, winding } = &t.loops()[l.index()] else {
                return Err(unwritable("a vertex loop"));
            };
            let mut uses: Vec<(usize, bool)> = fins
                .iter()
                .map(|k| {
                    let f = &t.fins()[k.index()];
                    (f.edge.index(), (f.sense == Orientation::Forward) != flip)
                })
                .collect();
            if flip {
                uses.reverse();
            }
            loops.push((uses, if flip { -winding[0] } else { winding[0] }));
        }
        let mut wires = Vec::new();
        let use_text = |(e, fwd): &(usize, bool)| {
            format!("{{{}{}}} 0", if *fwd { "+" } else { "-" }, edge_record[e])
        };
        if let Some(seam_edge) = seam_of.get(&fi) {
            // One wire: the bottom loop (+u) from the seam vertex, the seam
            // up, the top loop (-u) from its seam vertex, the seam down.
            let seam = &wound[&fi];
            let starting = |uses: &[(usize, bool)], end: &End| -> Vec<(usize, bool)> {
                let at = match end.vertex {
                    Vertex::Ring(_) => 0,
                    Vertex::Shared(v) => uses
                        .iter()
                        .position(|(e, fwd)| {
                            let edge = &t.edges()[*e];
                            let s = if *fwd { edge.start } else { edge.end };
                            s.map(|s| s.index()) == Some(v)
                        })
                        .unwrap_or(0),
                };
                uses[at..].iter().chain(&uses[..at]).copied().collect()
            };
            let (Some((b, _)), Some((u, _))) = (
                loops.iter().find(|(_, w)| *w > 0),
                loops.iter().find(|(_, w)| *w < 0),
            ) else {
                return Err(unwritable("a wound face without both wound loops"));
            };
            let mut text: Vec<String> = starting(b, &seam.bottom).iter().map(use_text).collect();
            text.push(format!("{{+{seam_edge}}} 0"));
            text.extend(starting(u, &seam.top).iter().map(use_text));
            text.push(format!("{{-{seam_edge}}} 0"));
            wires.push(records.push(format!("Wi\n\n0101100\n{} *", text.join(" "))));
            for (uses, _) in loops.iter().filter(|(_, w)| *w == 0) {
                let text: Vec<String> = uses.iter().map(use_text).collect();
                wires.push(records.push(format!("Wi\n\n0101100\n{} *", text.join(" "))));
            }
        } else {
            for (uses, w) in loops {
                if w != 0 {
                    return Err(unwritable("a wound loop on a face without a seam"));
                }
                let text: Vec<String> = uses
                    .iter()
                    .map(|(e, fwd)| {
                        format!("{{{}{}}} 0", if *fwd { "+" } else { "-" }, edge_record[e])
                    })
                    .collect();
                wires.push(records.push(format!("Wi\n\n0101100\n{} *", text.join(" "))));
            }
        }
        let text: Vec<String> = wires.iter().map(|w| format!("{{+{w}}} 0")).collect();
        face_record.insert(
            fi,
            records.push(format!(
                "Fa\n0 {tol} {} 0\n\n0101000\n{} *",
                fi + 1,
                text.join(" ")
            )),
        );
    }
    // Shells and solids: one solid per solid region, its shells' sides with
    // their orientation against the surface.
    let mut solids = Vec::new();
    for region in t.regions().iter().filter(|r| r.kind == RegionKind::Solid) {
        let mut shells = Vec::new();
        for s in &region.shells {
            let shell = &t.shells()[s.index()];
            let text: Vec<String> = shell
                .sides
                .iter()
                .map(|(f, side)| {
                    let face = &t.faces()[f.index()];
                    let forward = (face.sense == Orientation::Forward) == (*side == Side::Front);
                    format!(
                        "{{{}{}}} 0",
                        if forward { "+" } else { "-" },
                        face_record[&f.index()]
                    )
                })
                .collect();
            shells.push(records.push(format!("Sh\n\n0101100\n{} *", text.join(" "))));
        }
        let text: Vec<String> = shells.iter().map(|s| format!("{{+{s}}} 0")).collect();
        solids.push(records.push(format!("So\n\n0100100\n{} *", text.join(" "))));
    }
    let root = match solids.as_slice() {
        [] => return Err(unwritable("a topology without a solid region")),
        [one] => *one,
        many => {
            let text: Vec<String> = many.iter().map(|s| format!("{{+{s}}} 0")).collect();
            records.push(format!("Co\n\n1100000\n{} *", text.join(" ")))
        }
    };
    let _ = EdgeId::new(0);
    // Records are numbered backwards: the first written is n, the last 1.
    let n = records.shapes.len();
    let number = |i: usize| n - i;
    let mut out = String::from(
        "DBRep_DrawableShape\n\nCASCADE Topology V1, (c) Matra-Datavision\nLocations 0\n",
    );
    let _ = writeln!(out, "Curve2ds {}", curves2d.len());
    for c in &curves2d {
        let _ = writeln!(out, "{c}");
    }
    let _ = writeln!(out, "Curves {}", curves.len());
    for c in &curves {
        let _ = writeln!(out, "{c}");
    }
    let _ = writeln!(out, "Polygon3D 0\nPolygonOnTriangulations 0");
    let _ = writeln!(out, "Surfaces {}", surfaces.len());
    for s in &surfaces {
        let _ = writeln!(out, "{s}");
    }
    let _ = writeln!(out, "Triangulations 0\n\nTShapes {n}");
    for record in &records.shapes {
        // Resolve {v+i}, {v-i}, {+i} and {-i} placeholders into numbers.
        let mut text = String::new();
        let mut rest = record.as_str();
        while let Some(open) = rest.find('{') {
            text.push_str(&rest[..open]);
            let close = open + rest[open..].find('}').expect("closed placeholder");
            let inner = &rest[open + 1..close];
            let (sign, index) = match inner.strip_prefix('v') {
                Some(v) => (&v[..1], &v[1..]),
                None => (&inner[..1], &inner[1..]),
            };
            let i: usize = index.parse().expect("record index");
            let _ = write!(text, "{sign}{}", number(i));
            rest = &rest[close + 1..];
        }
        text.push_str(rest);
        let _ = writeln!(out, "{text}");
    }
    let _ = writeln!(out, "\n+{} 0 \n0", number(root));
    Ok(out)
}
