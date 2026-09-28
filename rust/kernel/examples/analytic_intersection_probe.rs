//! Test-only input for compare_analytic_intersections.py: reads the case
//! protocol of `analytic-intersection-cases.txt` on stdin and prints, per
//! case, `NAME empty|same|not_conic` or one `NAME ITEM` row per item: its
//! kind and every enclosed number as `lo hi`.
use rusty_occt::intersection::{surface_surface, SurfaceIntersection};
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
        "plane" => Surface::Plane(frame),
        "cylinder" => Surface::Cylinder {
            frame,
            radius: v[9],
        },
        "sphere" => Surface::Sphere {
            frame,
            radius: v[9],
        },
        _ => Surface::Cone {
            frame,
            radius: v[9],
            half_angle: v[10],
        },
    }
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
            Ok(SurfaceIntersection::Same) => println!("{name} same"),
            Ok(
                SurfaceIntersection::NotConic
                | SurfaceIntersection::Procedural(_)
                | SurfaceIntersection::Traced(_),
            ) => {
                println!("{name} not_conic")
            }
            Ok(SurfaceIntersection::Items(items)) => {
                for item in items {
                    let values: Vec<String> = item
                        .values()
                        .iter()
                        .map(|[lo, hi]| format!("{lo:?} {hi:?}"))
                        .collect();
                    println!("{name} {} {}", item.kind(), values.join(" "));
                }
            }
            Err(e) => println!("{name} error {e}"),
        }
    }
}
