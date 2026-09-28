//! Test-only input for compare_curve_surface.py: reads the case protocol of
//! `curve-surface-cases.txt` on stdin and prints per case `NAME empty`,
//! `NAME contained`, `NAME limit`, or one row per point sorted by parameter:
//! `NAME point s_lo s_hi x_lo x_hi y_lo y_hi z_lo z_hi crossing|tangent`.
use rusty_occt::intersection::{curve_surface, CurveSurfaceIntersection};
use rusty_occt::topology::{Curve3, Surface};
use rusty_occt::{Error, Frame3, Point3, Tolerance, Vec3};
use std::io::Read;

fn frame(v: &[f64]) -> Frame3 {
    Frame3::new(
        Point3::new(v[0], v[1], v[2]),
        Vec3::new(v[3], v[4], v[5]),
        Vec3::new(v[6], v[7], v[8]),
        Tolerance::default(),
    )
    .unwrap()
}

fn numbers(words: &[&str]) -> Vec<f64> {
    words.iter().map(|w| w.parse().unwrap()).collect()
}

fn curve(words: &[&str]) -> Curve3 {
    let v = numbers(&words[1..]);
    match words[0] {
        "line" => Curve3::LineSegment {
            start: Point3::new(v[0], v[1], v[2]),
            end: Point3::new(v[3], v[4], v[5]),
        },
        _ => Curve3::Circle {
            frame: frame(&v),
            radius: v[9],
        },
    }
}

fn surface(words: &[&str]) -> Surface {
    let v = numbers(&words[1..]);
    let frame = frame(&v);
    match words[0] {
        "plane" => Surface::Plane(frame),
        "sphere" => Surface::Sphere {
            frame,
            radius: v[9],
        },
        "cylinder" => Surface::Cylinder {
            frame,
            radius: v[9],
        },
        "torus" => Surface::Torus {
            frame,
            major: v[9],
            minor: v[10],
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
        let (c, s) = (curve(&lines[1][1..]), surface(&lines[2][1..]));
        match curve_surface(&c, &s) {
            Ok(CurveSurfaceIntersection::Empty) => println!("{name} empty"),
            Ok(CurveSurfaceIntersection::Contained) => println!("{name} contained"),
            Ok(CurveSurfaceIntersection::Points(points)) => {
                for p in points {
                    let mut words = vec![format!("{:?} {:?}", p.parameter[0], p.parameter[1])];
                    words.extend(p.point.iter().map(|[lo, hi]| format!("{lo:?} {hi:?}")));
                    let contact = if p.tangent { "tangent" } else { "crossing" };
                    println!("{name} point {} {contact}", words.join(" "));
                }
            }
            Err(Error::ComputationLimit(_)) => println!("{name} limit"),
            Err(e) => panic!("{name}: {e}"),
        }
    }
}
