use num_rational::BigRational as R;
use rusty_occt::{BSplineCurve3, ExactBSplineCurve3, Point3};

pub struct Input {
    pub name: String,
    pub curve: ExactBSplineCurve3,
    pub operations: Vec<(char, R, usize)>,
}
pub fn parse(row: &str) -> Input {
    let mut w = row.split_whitespace();
    let name = w.next().unwrap().to_owned();
    let periodic = w.next().unwrap() == "P";
    let degree: usize = w.next().unwrap().parse().unwrap();
    let np: usize = w.next().unwrap().parse().unwrap();
    let nk: usize = w.next().unwrap().parse().unwrap();
    let no: usize = w.next().unwrap().parse().unwrap();
    let mut poles = Vec::new();
    let mut weights = Vec::new();
    for _ in 0..np {
        let p: [f64; 4] = std::array::from_fn(|_| w.next().unwrap().parse().unwrap());
        poles.push(Point3::new(p[0], p[1], p[2]));
        weights.push(p[3]);
    }
    let mut knots = Vec::new();
    let mut mults = Vec::new();
    for _ in 0..nk {
        knots.push(w.next().unwrap().parse().unwrap());
        mults.push(w.next().unwrap().parse().unwrap());
    }
    let mut operations = Vec::new();
    for _ in 0..no {
        let op = w.next().unwrap().chars().next().unwrap();
        let u: f64 = w.next().unwrap().parse().unwrap();
        let m = w.next().unwrap().parse().unwrap();
        operations.push((op, R::from_float(u).unwrap(), m));
    }
    assert!(w.next().is_none());
    let curve = if periodic {
        BSplineCurve3::new_periodic(degree, poles, Some(weights), knots, mults)
    } else {
        BSplineCurve3::new(degree, poles, Some(weights), knots, mults)
    }
    .unwrap()
    .to_exact();
    Input {
        name,
        curve,
        operations,
    }
}
pub fn execute(input: &Input) -> (Vec<bool>, ExactBSplineCurve3) {
    let mut curve = input.curve.clone();
    let mut flags = Vec::new();
    for (op, u, m) in &input.operations {
        if *op == 'I' {
            curve = curve.insert_knot(u, *m).unwrap();
            flags.push(true);
        } else {
            assert_eq!(*op, 'R');
            if let Some(next) = curve.remove_knot(u, *m).unwrap() {
                curve = next;
                flags.push(true);
            } else {
                flags.push(false);
            }
        }
    }
    (flags, curve)
}
pub fn encode(name: &str, flags: &[bool], curve: &ExactBSplineCurve3) -> String {
    let mut w = vec![name.to_owned(), "R".into(), flags.len().to_string()];
    w.extend(flags.iter().map(|&v| usize::from(v).to_string()));
    w.extend([
        curve.degree().to_string(),
        usize::from(curve.is_periodic()).to_string(),
        curve.homogeneous_poles().len().to_string(),
        curve.knots().len().to_string(),
        curve.domain()[0].to_string(),
        curve.domain()[1].to_string(),
    ]);
    w.extend(
        curve
            .homogeneous_poles()
            .iter()
            .flatten()
            .map(ToString::to_string),
    );
    for (k, m) in curve.knots().iter().zip(curve.multiplicities()) {
        w.extend([k.to_string(), m.to_string()]);
    }
    w.join(" ")
}
