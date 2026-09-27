//! Test-only input for compare_procedural_intersections.py: reads the case
//! protocol of `procedural-intersection-cases.txt` on stdin and prints per
//! case `NAME empty`, `NAME point` with its enclosure, or the curve's rows in
//! the reference's order, every number as `lo hi`: `loop` (the range's two
//! enclosures, the points at each end, the points on both branches at the
//! range's middle), `rings` (points at 0 and pi on each branch) or
//! `figure_eight` (the node, and the points opposite it on each branch).
use rusty_occt::intersection::{
    surface_surface, AnalyticItem, Branch, Component, SurfaceIntersection,
};
use rusty_occt::topology::Surface;
use rusty_occt::{Frame3, Point3, Tolerance, Vec3};
use std::io::Read;

fn surface(words: &[&str]) -> Surface {
    let v: Vec<f64> = words[1..].iter().map(|w| w.parse().unwrap()).collect();
    let frame = Frame3::new(
        Point3::new(v[0], v[1], v[2]),
        Vec3::new(v[3], v[4], v[5]),
        Vec3::new(v[6], v[7], v[8]),
        Tolerance::default(),
    )
    .unwrap();
    match words[0] {
        "cylinder" => Surface::Cylinder {
            frame,
            radius: v[9],
        },
        _ => Surface::Sphere {
            frame,
            radius: v[9],
        },
    }
}

const PI: [f64; 2] = [std::f64::consts::PI, 3.1415926535897936];

fn text(values: &[[f64; 2]]) -> String {
    values
        .iter()
        .map(|[lo, hi]| format!("{lo:?} {hi:?}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for block in input.split("\nend").filter(|b| !b.trim().is_empty()) {
        let lines: Vec<Vec<&str>> = block
            .trim()
            .lines()
            .map(|l| l.split_whitespace().collect())
            .collect();
        let name = lines[0][1];
        let (a, b) = (surface(&lines[1][1..]), surface(&lines[2][1..]));
        match surface_surface(&a, &b) {
            Ok(SurfaceIntersection::Empty) => println!("{name} empty"),
            Ok(SurfaceIntersection::Items(items)) => {
                let [AnalyticItem::Point(p)] = items.as_slice() else {
                    panic!("{name}: {items:?}");
                };
                println!("{name} point {}", text(p));
            }
            Ok(SurfaceIntersection::Procedural(c)) => {
                let at = |u: [f64; 2], b| c.point_at(u, b).unwrap().to_vec();
                let mut values: Vec<[f64; 2]> = Vec::new();
                let kind = match c.components() {
                    [Component::Loop { u }] => {
                        values.extend(u);
                        values.extend(at(u[0], Branch::Plus));
                        values.extend(at(u[1], Branch::Plus));
                        let mid = 0.25 * (u[0][0] + u[0][1] + u[1][0] + u[1][1]);
                        values.extend(at([mid, mid], Branch::Plus));
                        values.extend(at([mid, mid], Branch::Minus));
                        "loop"
                    }
                    [Component::Ring { .. }, Component::Ring { .. }] => {
                        for u in [[0.0, 0.0], PI] {
                            values.extend(at(u, Branch::Plus));
                            values.extend(at(u, Branch::Minus));
                        }
                        "rings"
                    }
                    [Component::FigureEight { node }] => {
                        values.extend(at(*node, Branch::Plus));
                        values.extend(at([0.0, 0.0], Branch::Plus));
                        values.extend(at([0.0, 0.0], Branch::Minus));
                        "figure_eight"
                    }
                    other => panic!("{name}: {other:?}"),
                };
                println!("{name} {kind} {}", text(&values));
            }
            other => println!("{name} error {other:?}"),
        }
    }
}
