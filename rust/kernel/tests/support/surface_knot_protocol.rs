//! Test-only complete grid protocol. No expected values are computed here.
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::{BSplineSurface3, ExactBSplineSurface3, KnotVector, Point3};
use std::str::{FromStr, SplitWhitespace};
pub struct Input {
    pub name: String,
    pub surface: ExactBSplineSurface3,
    pub operations: Vec<(char, usize, R, usize)>,
}
fn token<T: FromStr>(w: &mut SplitWhitespace<'_>) -> T {
    w.next()
        .expect("missing token")
        .parse()
        .ok()
        .expect("invalid token")
}
pub fn parse(row: &str) -> Input {
    let mut w = row.split_whitespace();
    let name = token(&mut w);
    let degrees: [usize; 2] = std::array::from_fn(|_| token(&mut w));
    let counts: [usize; 2] = std::array::from_fn(|_| token(&mut w));
    let knots: [usize; 2] = std::array::from_fn(|_| token(&mut w));
    let periodic: [bool; 2] = std::array::from_fn(|_| token::<u8>(&mut w) != 0);
    let nop: usize = token(&mut w);
    let mut poles = Vec::new();
    let mut weights = Vec::new();
    for _ in 0..counts[0] * counts[1] {
        poles.push(Point3::new(token(&mut w), token(&mut w), token(&mut w)));
        weights.push(token(&mut w));
    }
    let axes: [KnotVector; 2] = std::array::from_fn(|i| {
        let mut k = Vec::new();
        let mut m = Vec::new();
        for _ in 0..knots[i] {
            k.push(token(&mut w));
            m.push(token(&mut w));
        }
        if periodic[i] {
            KnotVector::new_periodic(degrees[i], k, m)
        } else {
            KnotVector::new(degrees[i], k, m)
        }
        .unwrap()
    });
    let mut operations = Vec::new();
    for _ in 0..nop {
        let op = token(&mut w);
        let axis: char = token(&mut w);
        let at: f64 = token(&mut w);
        let target = token(&mut w);
        operations.push((
            op,
            match axis {
                'U' => 0,
                'V' => 1,
                _ => panic!("axis"),
            },
            R::from_float(at).unwrap(),
            target,
        ));
    }
    assert!(w.next().is_none());
    let [u, v] = axes;
    Input {
        name,
        surface: BSplineSurface3::new(u, v, poles, Some(weights))
            .unwrap()
            .to_exact(),
        operations,
    }
}
pub fn execute(input: &Input) -> (Vec<bool>, ExactBSplineSurface3) {
    let mut surface = input.surface.clone();
    let mut flags = Vec::new();
    for (op, axis, u, m) in &input.operations {
        let next = match (op, axis) {
            ('I', 0) => Some(surface.insert_u_knot(u, *m).unwrap()),
            ('I', 1) => Some(surface.insert_v_knot(u, *m).unwrap()),
            ('R', 0) => surface.remove_u_knot(u, *m).unwrap(),
            ('R', 1) => surface.remove_v_knot(u, *m).unwrap(),
            _ => panic!("operation"),
        };
        flags.push(next.is_some());
        if let Some(next) = next {
            surface = next;
        }
    }
    (flags, surface)
}
pub fn encode(name: &str, flags: &[bool], surface: &ExactBSplineSurface3, compact: bool) -> String {
    let axes = [surface.u_knots(), surface.v_knots()];
    let mut w = vec![
        name.to_owned(),
        if compact { "D" } else { "R" }.into(),
        flags.len().to_string(),
    ];
    w.extend(flags.iter().map(|&x| usize::from(x).to_string()));
    w.extend(surface.degrees().map(|d| d.to_string()));
    w.extend(axes.map(|a| usize::from(a.is_periodic()).to_string()));
    w.extend(surface.pole_counts().map(|n| n.to_string()));
    w.extend(axes.map(|a| a.knots().len().to_string()));
    w.extend(surface.domain().iter().flatten().map(ToString::to_string));
    if compact {
        let mut denominator = BigInt::from(1);
        for x in surface.homogeneous_poles().iter().flatten() {
            let (mut a, mut b) = (denominator.clone(), x.denom().clone());
            while b != BigInt::from(0) {
                (a, b) = (b.clone(), a % b);
            }
            denominator = denominator / a * x.denom();
        }
        w.push(denominator.to_string());
        w.extend(
            surface
                .homogeneous_poles()
                .iter()
                .flatten()
                .map(|x| (x.numer() * (&denominator / x.denom())).to_string()),
        );
    } else {
        w.extend(
            surface
                .homogeneous_poles()
                .iter()
                .flatten()
                .map(ToString::to_string),
        );
    }
    for a in axes {
        for (k, m) in a.knots().iter().zip(a.multiplicities()) {
            w.extend([k.to_string(), m.to_string()]);
        }
    }
    w.join(" ")
}
