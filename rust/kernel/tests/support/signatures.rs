//! Geometric signatures in the format of the native history probes
//! (rust/tools/occt_history_oracle.cpp, occt_split_merge_oracle.cpp): a vertex
//! point; an edge's curve kind and points at fractions 0, 1/2 and 1; a face's
//! surface kind, area and centroid; a body's volume and centroid. A region is
//! signed as its body. `counts` gives the synthesized OCCT counts.
use rusty_occt::topology::{Curve2, Curve3, Slot, Surface};
use rusty_occt::{Point2, Point3, Solid};

pub fn p3(p: Point3) -> String {
    format!("{:?} {:?} {:?}", p.x, p.y, p.z)
}

/// Twice the signed area and twice the first moments of one pcurve,
/// from the closed forms of 1/2 ∮ (u dv - v du), 1/2 ∮ u² dv, -1/2 ∮ v² du.
fn moments(p: &Curve2) -> [f64; 3] {
    match p {
        Curve2::LineSegment { start: a, end: b } => {
            let (du, dv) = (b.x - a.x, b.y - a.y);
            [
                a.x * b.y - b.x * a.y,
                dv * (a.x * a.x + a.x * du + du * du / 3.0),
                -du * (a.y * a.y + a.y * dv + dv * dv / 3.0),
            ]
        }
        Curve2::CircularArc {
            center: c,
            radius: r,
            start_angle: t0,
            sweep_angle: sweep,
        } => {
            let t1 = t0 + sweep;
            let (s0, c0, s1, c1) = (t0.sin(), t0.cos(), t1.sin(), t1.cos());
            let cube_s = |s: f64| s - s * s * s / 3.0;
            let cube_c = |c: f64| -c + c * c * c / 3.0;
            let sq_c = |t: f64| t / 2.0 + (2.0 * t).sin() / 4.0;
            let sq_s = |t: f64| t / 2.0 - (2.0 * t).sin() / 4.0;
            [
                r * (c.x * (s1 - s0) - c.y * (c1 - c0)) + r * r * sweep,
                r * (c.x * c.x * (s1 - s0)
                    + 2.0 * c.x * r * (sq_c(t1) - sq_c(*t0))
                    + r * r * (cube_s(s1) - cube_s(s0))),
                r * (-c.y * c.y * (c1 - c0)
                    + 2.0 * c.y * r * (sq_s(t1) - sq_s(*t0))
                    + r * r * (cube_c(c1) - cube_c(c0))),
            ]
        }
    }
}

pub fn signature(s: &Solid, slot: Slot) -> String {
    let t = s.topology();
    match slot {
        Slot::Vertex(v) => format!("V {}", p3(t.vertices()[v.index()].position)),
        Slot::Edge(e) => {
            let c = &t.edges()[e.index()].curve;
            let kind = if matches!(c, Curve3::LineSegment { .. }) {
                "line"
            } else {
                "circle"
            };
            format!(
                "E {kind} {} {} {}",
                p3(c.point(0.0)),
                p3(c.point(0.5)),
                p3(c.point(1.0))
            )
        }
        Slot::Face(f) => {
            let face = &t.faces()[f.index()];
            let fins = t.face_fins(f);
            let fins = fins.iter().flatten();
            match face.surface {
                Surface::Plane(frame) => {
                    let mut m = [0.0; 3];
                    for u in fins {
                        let x = moments(&u.pcurve);
                        (0..3).for_each(|i| m[i] += x[i]);
                    }
                    let (area, cu, cv) = (0.5 * m[0], m[1] / m[0], m[2] / m[0]);
                    let c = frame.point(Point2::new(cu, cv), 0.0);
                    format!("F plane {:?} {}", area.abs(), p3(c))
                }
                Surface::Cylinder { frame, radius } => {
                    // A full-turn wall of two ring loops: its area is the
                    // radius times the periodic area -∮ v du, its centroid on
                    // the axis at mid height.
                    let (mut periodic, mut v) = (0.0, (f64::MAX, f64::MIN));
                    for u in fins {
                        let Curve2::LineSegment { start: a, end: b } = u.pcurve else {
                            unreachable!("prism walls have line pcurves");
                        };
                        periodic -= 0.5 * (a.y + b.y) * (b.x - a.x);
                        v = (v.0.min(a.y), v.1.max(a.y));
                    }
                    let c = frame.point(Point2::default(), 0.5 * (v.0 + v.1));
                    format!("F cylinder {:?} {}", radius * periodic.abs(), p3(c))
                }
                Surface::Cone { .. } => {
                    let (area, centre) = t
                        .face_area_and_centre(f)
                        .expect("a cone face has certified mass properties");
                    format!("F cone {area:?} {}", p3(centre))
                }
            }
        }
        // A prism or cone has one solid region: the body.
        Slot::Region(_) => mass(s),
    }
}

pub fn mass(s: &Solid) -> String {
    let m = s.mass_properties();
    format!("S {:?} {}", m.volume, p3(m.centroid))
}

pub fn body(s: &Solid) -> String {
    format!("B {}", mass(s))
}

/// The counts OCCT reports for the same body, seams synthesized.
pub fn counts(s: &Solid) -> String {
    let c = s.topology().occt_counts();
    format!(
        "C {} {} {} {} {} {}",
        c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids
    )
}
