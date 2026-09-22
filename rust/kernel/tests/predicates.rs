use rusty_occt::predicates::{orient2d, Orientation2};
use rusty_occt::{Error, Point2};

fn sign(orientation: Orientation2) -> i8 {
    match orientation {
        Orientation2::Clockwise => -1,
        Orientation2::Collinear => 0,
        Orientation2::CounterClockwise => 1,
    }
}

#[test]
fn matches_independent_rational_oracle_and_all_permutations() {
    let mut cases = 0;
    for line in include_str!("../../fixtures/orient2d.tsv").lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        assert_eq!(fields.len(), 8);
        let values: Vec<_> = fields[1..7]
            .iter()
            .map(|s| f64::from_bits(u64::from_str_radix(s, 16).unwrap()))
            .collect();
        let points = [
            Point2::new(values[0], values[1]),
            Point2::new(values[2], values[3]),
            Point2::new(values[4], values[5]),
        ];
        let expected: i8 = fields[7].parse().unwrap();
        for (a, b, c, parity) in [
            (0, 1, 2, 1),
            (1, 2, 0, 1),
            (2, 0, 1, 1),
            (1, 0, 2, -1),
            (0, 2, 1, -1),
            (2, 1, 0, -1),
        ] {
            assert_eq!(
                sign(orient2d(points[a], points[b], points[c]).unwrap()),
                parity * expected,
                "{}: permutation {a},{b},{c}",
                fields[0]
            );
        }
        cases += 1;
    }
    assert_eq!(cases, 2417);
}

#[test]
fn cancellation_does_not_turn_a_positive_determinant_into_collinearity() {
    let u = 134_217_728.0_f64;
    // Exact determinant is u*u - (u-1)*(u+1) = 1. Naive f64 gives zero.
    assert_eq!(u * u - (u - 1.0) * (u + 1.0), 0.0);
    assert_eq!(
        orient2d(
            Point2::default(),
            Point2::new(u, u - 1.0),
            Point2::new(u + 1.0, u)
        )
        .unwrap(),
        Orientation2::CounterClockwise
    );
}

#[test]
fn every_coordinate_rejects_nonfinite_input() {
    for index in 0..6 {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut v = [0., 0., 1., 0., 0., 1.];
            v[index] = invalid;
            assert!(matches!(
                orient2d(
                    Point2::new(v[0], v[1]),
                    Point2::new(v[2], v[3]),
                    Point2::new(v[4], v[5])
                ),
                Err(Error::NonFinite(_))
            ));
        }
    }
}

#[test]
fn seeded_integer_cases_match_i128_and_exact_dyadic_transformations() {
    let mut state = 0xa73b_61c9_d421_502f_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % (1 << 40)) as i128 - (1 << 39)
    };
    for case in 0..10_000 {
        let v = [next(), next(), next(), next(), next(), next()];
        let determinant = (v[2] - v[0]) * (v[5] - v[1]) - (v[3] - v[1]) * (v[4] - v[0]);
        let expected = determinant.signum() as i8;
        // These translations and power-of-two scales are exactly representable;
        // arbitrary rounded transforms need not preserve the represented sign.
        let scale = 2.0_f64.powi(case % 81 - 40);
        let p: Vec<_> = v
            .chunks_exact(2)
            .map(|xy| Point2::new((xy[0] + 123) as f64 * scale, (xy[1] - 321) as f64 * scale))
            .collect();
        assert_eq!(
            sign(orient2d(p[0], p[1], p[2]).unwrap()),
            expected,
            "seeded case {case}"
        );
        let reflected: Vec<_> = p.iter().map(|p| Point2::new(-p.x, p.y)).collect();
        assert_eq!(
            sign(orient2d(reflected[0], reflected[1], reflected[2]).unwrap()),
            -expected,
            "reflection {case}"
        );
    }
}
