//! Line protocol written by rust/tools/brep_reference.py::encode, shared by
//! the fixture test and the native-comparison probe.
use rusty_occt::topology::{
    Coedge, Curve2, Curve3, Edge, EdgeId, Face, FaceId, Orientation, Surface, TopologyParts,
    Vertex, VertexId,
};
use rusty_occt::{Frame3, Point2, Point3, Tolerance, Vec3};

fn frame(v: &[f64]) -> Frame3 {
    Frame3::new(
        Point3::new(v[0], v[1], v[2]),
        Vec3::new(v[3], v[4], v[5]),
        Vec3::new(v[6], v[7], v[8]),
        Tolerance::default(),
    )
    .expect("fixture frames are valid")
}

fn orientation(word: &str) -> Orientation {
    if word == "F" {
        Orientation::Forward
    } else {
        Orientation::Reversed
    }
}

/// One case block: its name, linear tolerance and unvalidated parts.
pub fn parse(block: &str) -> (String, f64, TopologyParts) {
    let mut name = String::new();
    let mut tolerance = 0.0;
    let mut parts = TopologyParts::default();
    for line in block.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        let num = |i: usize| -> f64 { w[i].parse().unwrap() };
        let nums = |a: usize, n: usize| -> Vec<f64> { (a..a + n).map(num).collect() };
        match w[0] {
            "case" => name = w[1].to_string(),
            "tolerance" => tolerance = num(1),
            "v" => parts.vertices.push(Vertex {
                position: Point3::new(num(1), num(2), num(3)),
            }),
            "e" => {
                let (start, end) = (w[1].parse().unwrap(), w[2].parse().unwrap());
                let curve = if w[3] == "line" {
                    let v = nums(4, 6);
                    Curve3::LineSegment {
                        start: Point3::new(v[0], v[1], v[2]),
                        end: Point3::new(v[3], v[4], v[5]),
                    }
                } else {
                    let v = nums(4, 12);
                    Curve3::CircularArc {
                        frame: frame(&v),
                        radius: v[9],
                        start_angle: v[10],
                        sweep_angle: v[11],
                    }
                };
                parts.edges.push(Edge {
                    start: VertexId::new(start),
                    end: VertexId::new(end),
                    curve,
                });
            }
            "f" => {
                let (surface, o) = if w[1] == "plane" {
                    (Surface::Plane(frame(&nums(2, 9))), w[11])
                } else {
                    let v = nums(2, 10);
                    (
                        Surface::Cylinder {
                            frame: frame(&v),
                            radius: v[9],
                        },
                        w[12],
                    )
                };
                parts.faces.push(Face {
                    surface,
                    orientation: orientation(o),
                    loops: Vec::new(),
                });
            }
            "l" => parts.faces.last_mut().unwrap().loops.push(Vec::new()),
            "u" => {
                let pcurve = if w[3] == "line" {
                    let v = nums(4, 4);
                    Curve2::LineSegment {
                        start: Point2::new(v[0], v[1]),
                        end: Point2::new(v[2], v[3]),
                    }
                } else {
                    let v = nums(4, 5);
                    Curve2::CircularArc {
                        center: Point2::new(v[0], v[1]),
                        radius: v[2],
                        start_angle: v[3],
                        sweep_angle: v[4],
                    }
                };
                let face = parts.faces.last_mut().unwrap();
                face.loops.last_mut().unwrap().push(Coedge {
                    edge: EdgeId::new(w[1].parse().unwrap()),
                    orientation: orientation(w[2]),
                    pcurve,
                });
            }
            "s" => parts.shells.push(
                w[1..]
                    .iter()
                    .map(|f| FaceId::new(f.parse().unwrap()))
                    .collect(),
            ),
            "end" => {}
            other => panic!("unknown row {other}"),
        }
    }
    (name, tolerance, parts)
}
