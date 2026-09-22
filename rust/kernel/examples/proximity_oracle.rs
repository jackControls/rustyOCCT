//! Text bridge used only by the independent proximity checks.
use num_rational::BigRational as R;
use rusty_occt::proximity::{closest_points, LinearPrimitive3 as P};
use rusty_occt::{Error, Plane3, Point3, Result, ScalarInterval, Triangle3};
use std::io::{self, BufRead};

fn primitive(w: &mut std::str::SplitWhitespace<'_>) -> Result<P> {
    let kind = w.next().expect("kind");
    let mut point = || {
        let mut number = || w.next().expect("coordinate").parse().expect("number");
        Point3::new(number(), number(), number())
    };
    Ok(match kind {
        "P" => P::Point(point()),
        "L" => P::Line([point(), point()]),
        "S" => P::Segment([point(), point()]),
        "F" => P::Plane(Plane3::through_points(point(), point(), point())?),
        "T" => P::Triangle(Triangle3::new(point(), point(), point())?),
        _ => panic!("invalid kind"),
    })
}
fn rational(r: &R) {
    print!(" {}/{}", r.numer(), r.denom());
}
fn bounds(b: Result<ScalarInterval>) {
    match b {
        Ok(b) => print!(
            " E {:016x} {:016x}",
            b.lower().to_bits(),
            b.upper().to_bits()
        ),
        Err(Error::Unrepresentable(_)) => print!(" U"),
        Err(e) => panic!("unexpected bound failure: {e}"),
    }
}
fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let mut w = line.split_whitespace();
        let name = w.next().expect("case name");
        let a = primitive(&mut w).expect("valid reference primitive");
        let b = primitive(&mut w).expect("valid reference primitive");
        assert!(w.next().is_none());
        let pair = closest_points(&a, &b).expect("valid proximity");
        print!("{name} R");
        rational(pair.exact_squared_distance());
        bounds(pair.squared_distance_bounds());
        bounds(pair.distance_bounds());
        for i in 0..2 {
            print!(" {}", pair.exact_parameters()[i].len());
            for t in &pair.exact_parameters()[i] {
                rational(t);
            }
            for p in &pair.exact_points()[i] {
                rational(p);
            }
        }
        println!();
    }
}
