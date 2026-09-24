//! Test-only line protocol shared with the source-pinned OCCT edge/edge probe.
//! Each result parameter is printed as tight curve-parameter, line-parameter
//! and coordinate bounds; exact identities are checked separately in Python.
use rusty_occt::intersection::{
    spline_line_in, spline_segment_in, SplineLinearIntersection, SplineLinearPoint,
};
use rusty_occt::{BSplineCurve3, Point3};
use std::io::{self, BufRead};

fn bounds(point: &SplineLinearPoint) -> rusty_occt::Result<String> {
    let u = point.parameter_bounds()?;
    let s = point.linear_parameter_bounds()?;
    let mut words = vec![u.lower(), u.upper(), s.lower(), s.upper()];
    for c in point.coordinate_bounds()? {
        words.extend([c.lower(), c.upper()]);
    }
    Ok(words
        .iter()
        .map(|x| format!("{x:?}"))
        .collect::<Vec<_>>()
        .join(" "))
}

fn encode(result: &SplineLinearIntersection) -> rusty_occt::Result<String> {
    let mut out = format!("R {} {}", result.points().len(), result.overlaps().len());
    for point in result.points() {
        out += &format!(" P {}", bounds(point)?);
    }
    for overlap in result.overlaps() {
        let [a, b] = overlap.endpoints();
        out += &format!(" I {} {}", bounds(a)?, bounds(b)?);
    }
    Ok(out)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let Some(name) = words.next() else {
            continue;
        };
        let kind = words.next().ok_or("missing kind")?;
        let mut head = [0usize; 4];
        for h in &mut head {
            *h = words.next().ok_or("missing counts")?.parse()?;
        }
        let [degree, periodic, np, nk] = head;
        if np > 4096 || nk > 4122 || periodic > 1 || !matches!(kind, "L" | "S") {
            return Err("input counts or kind".into());
        }
        let values = words.map(str::parse).collect::<Result<Vec<f64>, _>>()?;
        if values.len() != 8 + 4 * np + 2 * nk {
            return Err("spline/linear input length".into());
        }
        let (first, last) = (values[0], values[1]);
        let at = |i: usize| Point3::new(values[i], values[i + 1], values[i + 2]);
        let (a, b) = (at(2), at(5));
        let poles = (0..np).map(|i| at(8 + 4 * i)).collect();
        let weights = (0..np).map(|i| values[11 + 4 * i]).collect();
        let tail = &values[8 + 4 * np..];
        let knots = tail.iter().step_by(2).copied().collect();
        let multiplicities = tail
            .iter()
            .skip(1)
            .step_by(2)
            .map(|&m| m as usize)
            .collect();
        let result = (if periodic == 0 {
            BSplineCurve3::new(degree, poles, Some(weights), knots, multiplicities)
        } else {
            BSplineCurve3::new_periodic(degree, poles, Some(weights), knots, multiplicities)
        })
        .and_then(|curve| {
            if kind == "L" {
                spline_line_in(&curve, a, b, first, last)
            } else {
                spline_segment_in(&curve, a, b, first, last)
            }
        })
        .and_then(|result| encode(&result));
        match result {
            Ok(encoded) => println!("{name} {encoded}"),
            Err(error) => println!("{name} E {error}"),
        }
    }
    Ok(())
}
