//! Test-only text protocol. No expected values are computed here.
use num_rational::BigRational as R;
use rusty_occt::{BSplineSurface3, ExactBezierCurve3, ExactBezierSurface3, KnotVector, Point3};
use std::str::{FromStr, SplitWhitespace};

pub struct Input {
    pub name: String,
    pub surface: BSplineSurface3,
    pub rectangle: [f64; 4],
    pub op: u8,
    pub elevation: [usize; 2],
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
    assert!(matches!(kind.as_str(), "S" | "B"));
    let degrees: [usize; 2] = std::array::from_fn(|_| token(&mut w));
    let counts: [usize; 2] = std::array::from_fn(|_| token(&mut w));
    let knots: [usize; 2] = std::array::from_fn(|_| token(&mut w));
    let periodic: [bool; 2] = std::array::from_fn(|_| token::<u8>(&mut w) != 0);
    let rectangle = std::array::from_fn(|_| token(&mut w));
    let op = token(&mut w);
    let elevation = std::array::from_fn(|_| token(&mut w));
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
    assert!(w.next().is_none());
    let [u, v] = axes;
    Input {
        name,
        surface: BSplineSurface3::new(u, v, poles, Some(weights)).unwrap(),
        rectangle,
        op,
        elevation,
    }
}

pub enum Item {
    Patch(ExactBezierSurface3),
    Curve(ExactBezierCurve3),
}
fn at(domain: &[R; 2], n: i32, d: i32) -> R {
    &domain[0] + (&domain[1] - &domain[0]) * R::new(n.into(), d.into())
}
pub fn apply(mut p: ExactBezierSurface3, op: u8, elevation: [usize; 2]) -> Vec<Item> {
    assert!(op <= 10);
    if op == 1 || op == 10 {
        let [u, v] = p.domain();
        p = p
            .trim(&at(u, 1, 4), &at(u, 3, 4), &at(v, 1, 4), &at(v, 3, 4))
            .unwrap();
    }
    if op == 7 || op == 10 {
        p = p.elevated(elevation[0], elevation[1]).unwrap();
    }
    if op == 4 || op == 10 {
        p = p.u_reversed();
    }
    if op == 5 {
        p = p.v_reversed();
    }
    if op == 6 || op == 10 {
        p = p.exchanged_uv();
    }
    if op == 8 {
        return vec![Item::Curve(p.u_iso(&at(&p.domain()[0], 3, 8)).unwrap())];
    }
    if op == 9 {
        return vec![Item::Curve(p.v_iso(&at(&p.domain()[1], 5, 8)).unwrap())];
    }
    let us = if op == 2 || op == 10 {
        p.split_u_at(&at(&p.domain()[0], 3, 8)).unwrap().to_vec()
    } else {
        vec![p]
    };
    us.into_iter()
        .flat_map(|p| {
            if op == 3 || op == 10 {
                p.split_v_at(&at(&p.domain()[1], 5, 8))
                    .unwrap()
                    .into_iter()
                    .map(Item::Patch)
                    .collect()
            } else {
                vec![Item::Patch(p)]
            }
        })
        .collect()
}
pub fn execute(input: &Input) -> Vec<Item> {
    let [ua, ub, va, vb] = input.rectangle;
    input
        .surface
        .bezier_patches_in(ua, ub, va, vb)
        .unwrap()
        .into_iter()
        .flat_map(|p| apply(p, input.op, input.elevation))
        .collect()
}
#[allow(dead_code)]
pub fn encode(name: &str, items: &[Item]) -> String {
    let mut w = vec![name.to_owned(), "R".to_owned(), items.len().to_string()];
    for item in items {
        match item {
            Item::Patch(p) => {
                w.push("P".to_owned());
                w.extend(p.degrees().map(|d| d.to_string()));
                w.extend(p.domain().iter().flatten().map(ToString::to_string));
                w.extend(
                    p.homogeneous_poles()
                        .iter()
                        .flatten()
                        .map(ToString::to_string),
                );
            }
            Item::Curve(c) => {
                w.extend(["C".to_owned(), c.degree().to_string(), "0".to_owned()]);
                w.extend(c.domain().iter().map(ToString::to_string));
                w.extend(
                    c.homogeneous_poles()
                        .iter()
                        .flatten()
                        .map(ToString::to_string),
                );
            }
        }
    }
    w.join(" ")
}

/// Lossless fixture encoding: one least common control denominator per item.
/// This keeps full independently computed controls without repeating large
/// denominator strings. The live bridge continues to use reduced rationals.
#[allow(dead_code)]
pub fn encode_fixture(name: &str, items: &[Item]) -> String {
    use num_bigint::BigInt;
    let mut w = vec![name.to_owned(), "D".to_owned(), items.len().to_string()];
    for item in items {
        let controls = match item {
            Item::Patch(p) => {
                w.push("P".to_owned());
                w.extend(p.degrees().map(|d| d.to_string()));
                w.extend(p.domain().iter().flatten().map(ToString::to_string));
                p.homogeneous_poles()
            }
            Item::Curve(c) => {
                w.extend(["C".to_owned(), c.degree().to_string(), "0".to_owned()]);
                w.extend(c.domain().iter().map(ToString::to_string));
                c.homogeneous_poles()
            }
        };
        let mut denominator = BigInt::from(1);
        for x in controls.iter().flatten() {
            let (mut a, mut b) = (denominator.clone(), x.denom().clone());
            while b != BigInt::from(0) {
                (a, b) = (b.clone(), a % b);
            }
            denominator = denominator / a * x.denom();
        }
        w.push(denominator.to_string());
        w.extend(
            controls
                .iter()
                .flatten()
                .map(|x| (x.numer() * (&denominator / x.denom())).to_string()),
        );
    }
    w.join(" ")
}
