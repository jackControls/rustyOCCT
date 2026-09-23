//! Test-only line protocol shared with the source-pinned OCCT projection probe.
use rusty_occt::{proximity::closest_points_on_spline_in, BSplineCurve3, Point3};
use std::io::{self, BufRead};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let timings = std::env::args().any(|a| a == "--timings");
    for line in io::stdin().lock().lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let Some(name) = words.next() else {
            continue;
        };
        let mut next = || words.next().ok_or("missing projection input");
        let degree: usize = next()?.parse()?;
        let periodic: u8 = next()?.parse()?;
        let np: usize = next()?.parse()?;
        let nk: usize = next()?.parse()?;
        if np > 4096 || nk > 4122 || periodic > 1 {
            return Err("input counts".into());
        }
        let first: f64 = next()?.parse()?;
        let last: f64 = next()?.parse()?;
        let point = Point3::new(next()?.parse()?, next()?.parse()?, next()?.parse()?);
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for _ in 0..np {
            poles.push(Point3::new(
                next()?.parse()?,
                next()?.parse()?,
                next()?.parse()?,
            ));
            weights.push(next()?.parse()?);
        }
        let mut knots = Vec::new();
        let mut multiplicities = Vec::new();
        for _ in 0..nk {
            knots.push(next()?.parse()?);
            multiplicities.push(next()?.parse()?);
        }
        if words.next().is_some() {
            return Err("trailing projection input".into());
        }
        let started = Instant::now();
        let result = (if periodic == 0 {
            BSplineCurve3::new(degree, poles, Some(weights), knots, multiplicities)
        } else {
            BSplineCurve3::new_periodic(degree, poles, Some(weights), knots, multiplicities)
        })
        .and_then(|curve| {
            closest_points_on_spline_in(&curve, point, first, last, Default::default())
        });
        if timings {
            eprintln!("{name} solve {:?}", started.elapsed());
        }
        match result {
            Err(error) => println!("{name} E {error}"),
            Ok(result) => {
                let distance = result.squared_distance_bounds()?;
                if timings {
                    eprintln!("{name} distance_bounds {:?}", started.elapsed());
                }
                print!(
                    "{name} R {} {} {:.17e} {:.17e}",
                    result.points().len(),
                    result.intervals().len(),
                    distance.lower(),
                    distance.upper()
                );
                for (i, point) in result.points().iter().enumerate() {
                    let p = point.parameter_bounds()?;
                    print!(" P {:.17e} {:.17e}", p.lower(), p.upper());
                    for coordinate in point.coordinate_bounds()? {
                        print!(" {:.17e} {:.17e}", coordinate.lower(), coordinate.upper());
                    }
                    if timings {
                        eprintln!("{name} point_{i}_bounds {:?}", started.elapsed());
                    }
                }
                for interval in result.intervals() {
                    print!(
                        " I {} {}",
                        interval.parameters()[0],
                        interval.parameters()[1]
                    );
                }
                println!();
            }
        }
    }
    Ok(())
}
