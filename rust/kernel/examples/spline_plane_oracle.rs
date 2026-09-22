//! Test-only text bridge. Each hit exports all certified bounds and contact data.
use rusty_occt::intersection::{spline_plane, spline_plane_in, Plane3, SplinePlaneContact};
use rusty_occt::{BSplineCurve3, Point3};
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let name = words.next().ok_or("missing label")?;
        let kind = words.next().ok_or("missing kind")?;
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
        let plane = Plane3::through_points(points[0], points[1], points[2])?;
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
        let result = if let Some((a, b)) = range {
            spline_plane_in(&curve, &plane, a, b)?
        } else {
            spline_plane(&curve, &plane)?
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
                SplinePlaneContact::Crossing => "C",
                SplinePlaneContact::Tangent => "T",
                SplinePlaneContact::Boundary => "B",
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
