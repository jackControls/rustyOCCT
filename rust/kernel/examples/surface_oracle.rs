use rusty_occt::curve::{DerivativeOrder as D, KnotSide as S};
use rusty_occt::{BSplineSurface3, BezierSurface3, KnotVector, Point3};
use std::io::{self, BufRead};

fn side(word: &str) -> S {
    match word {
        "A" => S::Automatic,
        "L" => S::Left,
        "R" => S::Right,
        _ => panic!("invalid side"),
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let line = line?;
        let words: Vec<_> = line.split_whitespace().collect();
        let [du, dv, nu, nv, ku, kv, pu, pv]: [usize; 8] =
            std::array::from_fn(|i| words[i + 2].parse().unwrap());
        let u = words[10].parse()?;
        let v = words[11].parse()?;
        let sides = [side(words[12]), side(words[13])];
        let order = match words[14] {
            "0" => D::Position,
            "1" => D::First,
            "2" => D::Second,
            _ => panic!("invalid order"),
        };
        let mut data = words[15..].iter();
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for _ in 0..nu * nv {
            let [x, y, z, w] =
                std::array::from_fn(|_| data.next().unwrap().parse::<f64>().unwrap());
            poles.push(Point3::new(x, y, z));
            weights.push(w);
        }
        let mut axis = |degree, count, periodic| {
            let mut knots = Vec::new();
            let mut mults = Vec::new();
            for _ in 0..count {
                knots.push(data.next().unwrap().parse().unwrap());
                mults.push(data.next().unwrap().parse().unwrap());
            }
            if periodic {
                KnotVector::new_periodic(degree, knots, mults)
            } else {
                KnotVector::new(degree, knots, mults)
            }
        };
        let u_axis = axis(du, ku, pu != 0)?;
        let v_axis = axis(dv, kv, pv != 0)?;
        assert!(data.next().is_none());
        let value = match words[1] {
            "B" => BezierSurface3::new(du, dv, poles, Some(weights))?.evaluate(u, v, order)?,
            "S" => BSplineSurface3::new(u_axis, v_axis, poles, Some(weights))?
                .evaluate(u, v, order, sides)?,
            _ => panic!("invalid kind"),
        };
        print!("{}", words[0]);
        for row in [
            Some(value.position_bounds()),
            value.derivative_bounds(1, 0),
            value.derivative_bounds(0, 1),
            value.derivative_bounds(2, 0),
            value.derivative_bounds(0, 2),
            value.derivative_bounds(1, 1),
        ]
        .into_iter()
        .flatten()
        {
            for x in row {
                print!(" {:.17e}", x.representative());
            }
        }
        println!();
    }
    Ok(())
}
