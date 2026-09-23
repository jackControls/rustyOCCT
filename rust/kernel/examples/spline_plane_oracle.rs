//! Test-only text bridge. Each hit exports all certified bounds and contact data.
use num_rational::BigRational as R;
use rusty_occt::intersection::{
    exact_spline_cylinder_in, exact_spline_plane_in, exact_spline_sphere_in, spline_cylinder_in,
    spline_plane, spline_plane_in, spline_sphere_in, Cylinder3, Plane3, Sphere3,
    SplineSurfaceContact,
};
use rusty_occt::{BSplineCurve3, Point3, Vec3};
use std::io::{self, BufRead};
#[path = "../tests/support/exact_spline_edits.rs"]
mod exact_edits;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let mode = args.first().map_or("binary64", String::as_str);
    if args.len() > 1 || !["binary64", "exact", "refined", "roundtrip"].contains(&mode) {
        return Err("expected binary64, exact, refined or roundtrip".into());
    }
    for line in io::stdin().lock().lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let name = words.next().ok_or("missing label")?;
        let kind = words.next().ok_or("missing kind")?;
        let (surface_kind, kind) = kind.split_once(':').unwrap_or(("plane", kind));
        let degree = words.next().ok_or("missing degree")?.parse()?;
        let count: usize = words.next().ok_or("missing pole count")?.parse()?;
        let knot_count: usize = words.next().ok_or("missing knot count")?.parse()?;
        let mut number = || -> Result<f64, Box<dyn std::error::Error>> {
            Ok(words.next().ok_or("incomplete curve or plane")?.parse()?)
        };
        let mut points = Vec::new();
        for _ in 0..3 {
            points.push(Point3::new(number()?, number()?, number()?));
        }
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for _ in 0..count {
            poles.push(Point3::new(number()?, number()?, number()?));
            weights.push(number()?);
        }
        let mut knots = Vec::new();
        let mut mults = Vec::new();
        for _ in 0..knot_count {
            knots.push(words.next().ok_or("missing knot")?.parse()?);
            mults.push(words.next().ok_or("missing multiplicity")?.parse()?);
        }
        let range = if kind.ends_with('T') {
            Some((
                words.next().ok_or("missing first parameter")?.parse()?,
                words.next().ok_or("missing last parameter")?.parse()?,
            ))
        } else {
            None
        };
        if words.next().is_some() {
            return Err("extra input data".into());
        }
        let curve = match kind {
            "B" | "S" | "BT" | "ST" => {
                BSplineCurve3::new(degree, poles, Some(weights), knots, mults)?
            }
            "P" | "PT" => BSplineCurve3::new_periodic(degree, poles, Some(weights), knots, mults)?,
            _ => return Err("invalid curve kind".into()),
        };
        let (a, b) = range.unwrap_or_else(|| curve.domain());
        let result = if mode != "binary64" {
            let original = curve.to_exact();
            let exact = if mode == "exact" {
                original
            } else {
                let (fine, cuts) = exact_edits::refined(&original);
                if mode == "refined" {
                    fine
                } else {
                    let restored = exact_edits::restored(&fine, &cuts);
                    assert_eq!(restored, original);
                    restored
                }
            };
            let (a, b) = (R::from_float(a).unwrap(), R::from_float(b).unwrap());
            match surface_kind {
                "plane" => exact_spline_plane_in(
                    &exact,
                    &Plane3::through_points(points[0], points[1], points[2])?,
                    &a,
                    &b,
                )?,
                "sphere" => {
                    exact_spline_sphere_in(&exact, &Sphere3::new(points[0], points[2].x)?, &a, &b)?
                }
                "cylinder" => exact_spline_cylinder_in(
                    &exact,
                    &Cylinder3::new(
                        points[0],
                        Vec3::new(points[1].x, points[1].y, points[1].z),
                        points[2].x,
                    )?,
                    &a,
                    &b,
                )?,
                _ => return Err("invalid surface kind".into()),
            }
            .enclosed()?
        } else {
            match surface_kind {
                "plane" => {
                    let plane = Plane3::through_points(points[0], points[1], points[2])?;
                    if range.is_some() {
                        spline_plane_in(&curve, &plane, a, b)?
                    } else {
                        spline_plane(&curve, &plane)?
                    }
                }
                "sphere" => spline_sphere_in(&curve, &Sphere3::new(points[0], points[2].x)?, a, b)?,
                "cylinder" => spline_cylinder_in(
                    &curve,
                    &Cylinder3::new(
                        points[0],
                        Vec3::new(points[1].x, points[1].y, points[1].z),
                        points[2].x,
                    )?,
                    a,
                    b,
                )?,
                _ => return Err("invalid surface kind".into()),
            }
        };
        print!(
            "{name} {} {}",
            result.points().len(),
            result.overlaps().len()
        );
        for point in result.points() {
            for bounds in std::iter::once(point.parameter()).chain(point.coordinate_bounds()) {
                print!(" {:.17e} {:.17e}", bounds.lower(), bounds.upper());
            }
            let contact = match point.contact() {
                SplineSurfaceContact::Crossing => "C",
                SplineSurfaceContact::Tangent => "T",
                SplineSurfaceContact::Boundary => "B",
            };
            let [left, right] = point.multiplicities().map(|m| m.unwrap_or(0));
            print!(" {contact} {left} {right}");
        }
        for overlap in result.overlaps() {
            if range.is_some() {
                for b in overlap.parameter_bounds() {
                    print!(" {:.17e} {:.17e}", b.lower(), b.upper());
                }
            } else {
                let (a, b) = overlap.parameters();
                print!(" {a:.17e} {b:.17e}");
            }
        }
        println!();
    }
    Ok(())
}
