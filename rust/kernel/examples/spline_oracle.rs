use rusty_occt::curve::{BSplineCurve3, BezierCurve3, DerivativeOrder, KnotSide};
use rusty_occt::Point3;
use std::io::{self, BufRead};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let line = line?;
        let mut words = line.split_whitespace();
        let name = words.next().ok_or("missing label")?;
        let kind = words.next().ok_or("missing kind")?;
        let degree = words.next().ok_or("missing degree")?.parse()?;
        let count: usize = words.next().ok_or("missing count")?.parse()?;
        let knot_count: usize = words.next().ok_or("missing knot count")?.parse()?;
        let u = words.next().ok_or("missing parameter")?.parse()?;
        let side = match words.next() {
            Some("L") => KnotSide::Left,
            Some("R") => KnotSide::Right,
            Some("A") => KnotSide::Automatic,
            _ => return Err("invalid side".into()),
        };
        let order = match words.next() {
            Some("0") => DerivativeOrder::Position,
            Some("1") => DerivativeOrder::First,
            Some("2") => DerivativeOrder::Second,
            _ => return Err("invalid order".into()),
        };
        let mut number = || -> Result<f64, Box<dyn std::error::Error>> {
            Ok(words.next().ok_or("incomplete curve")?.parse()?)
        };
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
        if words.next().is_some() {
            return Err("extra curve data".into());
        }
        let value = match kind {
            "B" => BezierCurve3::new(poles, Some(weights))?.evaluate(u, order)?,
            "S" => BSplineCurve3::new(degree, poles, Some(weights), knots, mults)?
                .evaluate(u, order, side)?,
            _ => return Err("invalid curve type".into()),
        };
        print!("{name}");
        for bounds in [
            Some(value.position_bounds()),
            value.derivative_bounds(1),
            value.derivative_bounds(2),
        ]
        .into_iter()
        .flatten()
        {
            for x in bounds {
                print!(" {:.17e}", x.representative());
            }
        }
        println!();
    }
    Ok(())
}
