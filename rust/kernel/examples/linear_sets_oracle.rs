//! Text bridge for independent complete-set validation; not a runtime API.
use rusty_occt::intersection::{
    linear_intersection, LinearIntersection as I, LinearPrimitive3 as P,
};
use rusty_occt::{Plane3, Point3, Triangle3};
use std::io::{self, BufRead};
fn primitive(w: &mut std::str::SplitWhitespace<'_>) -> P {
    let kind = w.next().expect("kind");
    let mut point = || {
        let mut x = || w.next().expect("coordinate").parse().expect("number");
        Point3::new(x(), x(), x())
    };
    match kind {
        "P" => P::Point(point()),
        "L" => P::Line([point(), point()]),
        "S" => P::Segment([point(), point()]),
        "F" => P::Plane(Plane3::through_points(point(), point(), point()).unwrap()),
        "T" => P::Triangle(Triangle3::new(point(), point(), point()).unwrap()),
        _ => panic!("kind"),
    }
}
fn main() {
    for row in io::stdin().lock().lines() {
        let row = row.unwrap();
        let mut w = row.split_whitespace();
        let name = w.next().unwrap();
        let (a, b) = (primitive(&mut w), primitive(&mut w));
        assert!(w.next().is_none());
        let result = linear_intersection(&a, &b).unwrap();
        let (kind, values) = match &result {
            I::Empty => ("E", Vec::new()),
            I::Point(_) => (
                "P",
                result
                    .vertices()
                    .iter()
                    .map(|p| p.coordinates().to_vec())
                    .collect(),
            ),
            I::Segment(_) => (
                "S",
                result
                    .vertices()
                    .iter()
                    .map(|p| p.coordinates().to_vec())
                    .collect(),
            ),
            I::Polygon(_) => (
                "G",
                result
                    .vertices()
                    .iter()
                    .map(|p| p.coordinates().to_vec())
                    .collect(),
            ),
            I::Line { origin, direction } => {
                ("L", vec![origin.coordinates().to_vec(), direction.to_vec()])
            }
            I::Plane { normal, offset } => ("F", vec![normal.to_vec(), vec![offset.clone()]]),
        };
        print!("{name} R {kind} {}", values.len());
        for v in values {
            for x in v {
                print!(" {}/{}", x.numer(), x.denom());
            }
        }
        println!();
    }
}
