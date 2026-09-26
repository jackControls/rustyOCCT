//! Line protocol written by rust/tools/cell_reference.py::encode (the cell
//! model): arenas in order (vertices, edges, loops with their fins, faces,
//! shells, regions), shared by the fixture test and the native-comparison probe.
use rusty_occt::topology::{
    Curve2, Curve3, Edge, EdgeId, Enclosure, Face, FaceId, Fin, FinId, Loop, LoopId, Orientation,
    Region, RegionId, RegionKind, Shell, ShellId, Side, Surface, TopologyParts, Vertex, VertexId,
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

fn vertex(word: &str) -> Option<VertexId> {
    (word != "-").then(|| VertexId::new(word.parse().unwrap()))
}

/// The words after `key` up to the next keyword (or the end).
fn section<'a>(w: &[&'a str], key: &str, keys: &[&str]) -> Vec<&'a str> {
    let Some(at) = w.iter().position(|x| *x == key) else {
        return Vec::new();
    };
    w[at + 1..]
        .iter()
        .take_while(|x| !keys.contains(x))
        .copied()
        .collect()
}

/// A declared enclosure: `enc BOUND` (or `enc -`, or no `enc`, for none).
fn enclosure(w: &[&str]) -> Option<Enclosure> {
    let at = w.iter().position(|x| *x == "enc")?;
    (w[at + 1] != "-").then(|| Enclosure::computed(w[at + 1].parse().unwrap()))
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
        let ids = |key: &str| -> Vec<usize> {
            section(&w, key, &["fins", "sides", "wire", "acorn", "loops", "enc"])
                .iter()
                .map(|x| x.parse().unwrap())
                .collect()
        };
        match w[0] {
            "case" => name = w[1].to_string(),
            "tolerance" => tolerance = num(1),
            "v" => parts.vertices.push(Vertex {
                position: Point3::new(num(1), num(2), num(3)),
                enclosure: enclosure(&w),
            }),
            "e" => {
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
                    start: vertex(w[1]),
                    end: vertex(w[2]),
                    curve,
                    fins: ids("fins").into_iter().map(FinId::new).collect(),
                });
            }
            "f" => {
                let (surface, rest) = if w[1] == "plane" {
                    (Surface::Plane(frame(&nums(2, 9))), 11)
                } else if w[1] == "cone" {
                    let v = nums(2, 11);
                    (
                        Surface::Cone {
                            frame: frame(&v),
                            radius: v[9],
                            half_angle: v[10],
                        },
                        13,
                    )
                } else {
                    let v = nums(2, 10);
                    (
                        Surface::Cylinder {
                            frame: frame(&v),
                            radius: v[9],
                        },
                        12,
                    )
                };
                parts.faces.push(Face {
                    surface,
                    sense: orientation(w[rest]),
                    loops: ids("loops").into_iter().map(LoopId::new).collect(),
                    front: ShellId::new(w[rest + 1].parse().unwrap()),
                    back: ShellId::new(w[rest + 2].parse().unwrap()),
                    enclosure: enclosure(&w),
                });
            }
            "l" => {
                parts.loops.push(Loop::Edges {
                    fins: Vec::new(),
                    winding: [w[1].parse().unwrap(), 0],
                });
            }
            "lv" => {
                parts
                    .loops
                    .push(Loop::Vertex(VertexId::new(w[1].parse().unwrap())));
            }
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
                let id = FinId::new(parts.fins.len());
                parts.fins.push(Fin {
                    edge: EdgeId::new(w[1].parse().unwrap()),
                    sense: orientation(w[2]),
                    pcurve,
                    enclosure: enclosure(&w),
                });
                if let Some(Loop::Edges { fins, .. }) = parts.loops.last_mut() {
                    fins.push(id);
                }
            }
            "s" => {
                let sides = section(&w, "sides", &["wire", "acorn"])
                    .iter()
                    .map(|x| {
                        let (f, side) = x.split_once(':').unwrap();
                        let side = if side == "F" { Side::Front } else { Side::Back };
                        (FaceId::new(f.parse().unwrap()), side)
                    })
                    .collect();
                parts.shells.push(Shell {
                    region: RegionId::new(w[1].parse().unwrap()),
                    sides,
                    wire_edges: ids("wire").into_iter().map(EdgeId::new).collect(),
                    acorn_vertices: ids("acorn").into_iter().map(VertexId::new).collect(),
                });
            }
            "r" => parts.regions.push(Region {
                kind: if w[1] == "solid" {
                    RegionKind::Solid
                } else {
                    RegionKind::Void
                },
                shells: w[2..]
                    .iter()
                    .map(|x| ShellId::new(x.parse().unwrap()))
                    .collect(),
            }),
            "end" => {}
            other => panic!("unknown row {other}"),
        }
    }
    (name, tolerance, parts)
}
