//! Test-only text protocol. Expected values are computed outside this module.
use num_rational::BigRational as R;
use rusty_occt::{BSplineCurve3, ExactBezierCurve3, Point3};
use std::str::{FromStr, SplitWhitespace};
pub struct Input {
    pub name: String,
    pub curve: BSplineCurve3,
    pub first: f64,
    pub last: f64,
    pub op: u8,
    pub elevation: usize,
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
    let kind: String = token(&mut w);
    let (degree, np, nk) = (
        token(&mut w),
        token::<usize>(&mut w),
        token::<usize>(&mut w),
    );
    let (first, last, op, elevation) = (token(&mut w), token(&mut w), token(&mut w), token(&mut w));
    let mut poles = Vec::new();
    let mut weights = Vec::new();
    for _ in 0..np {
        poles.push(Point3::new(token(&mut w), token(&mut w), token(&mut w)));
        weights.push(token(&mut w));
    }
    let mut knots = Vec::new();
    let mut mults = Vec::new();
    for _ in 0..nk {
        knots.push(token(&mut w));
        mults.push(token(&mut w));
    }
    assert!(w.next().is_none());
    assert!(matches!(kind.as_str(), "P" | "S" | "B"));
    let build = if kind == "P" {
        BSplineCurve3::new_periodic
    } else {
        BSplineCurve3::new
    };
    Input {
        name,
        curve: build(degree, poles, Some(weights), knots, mults).unwrap(),
        first,
        last,
        op,
        elevation,
    }
}
pub fn apply(curve: ExactBezierCurve3, op: u8, elevation: usize) -> Vec<ExactBezierCurve3> {
    let mut curve = curve;
    assert!(op <= 5);
    if op == 1 || op == 5 {
        let [a, b] = curve.domain();
        let length = b - a;
        curve = curve
            .trim(
                &(a + &length / R::from_integer(4.into())),
                &(a + length * R::new(3.into(), 4.into())),
            )
            .unwrap();
    }
    if op == 4 || op == 5 {
        curve = curve.elevated(elevation).unwrap();
    }
    if op == 3 || op == 5 {
        curve = curve.reversed();
    }
    if op == 2 || op == 5 {
        let [a, b] = curve.domain();
        let cut = a + (b - a) * R::new(3.into(), 8.into());
        curve.split_at(&cut).unwrap().to_vec()
    } else {
        vec![curve]
    }
}
pub fn execute(input: &Input) -> Vec<ExactBezierCurve3> {
    input
        .curve
        .bezier_arcs_in(input.first, input.last)
        .unwrap()
        .into_iter()
        .flat_map(|c| apply(c, input.op, input.elevation))
        .collect()
}
pub fn encode(name: &str, curves: &[ExactBezierCurve3]) -> String {
    let mut w = vec![name.to_owned(), "R".to_owned(), curves.len().to_string()];
    for c in curves {
        w.extend([
            c.degree().to_string(),
            c.domain()[0].to_string(),
            c.domain()[1].to_string(),
        ]);
        w.extend(
            c.homogeneous_poles()
                .iter()
                .flatten()
                .map(ToString::to_string),
        );
    }
    w.join(" ")
}
