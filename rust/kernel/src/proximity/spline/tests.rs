use super::*;
use std::{fmt::Debug, str::FromStr};

struct Words<'a>(std::str::SplitWhitespace<'a>);
impl<'a> Words<'a> {
    fn new(line: &'a str) -> Self {
        Self(line.split_whitespace())
    }
    fn take<T: FromStr>(&mut self) -> T
    where
        T::Err: Debug,
    {
        self.0.next().unwrap().parse().unwrap()
    }
    fn polynomial(&mut self) -> Vec<R> {
        let n = self.take();
        (0..n).map(|_| self.take()).collect()
    }
    fn end(&mut self) {
        assert!(self.0.next().is_none());
    }
}

// Test-only affine substitution by polynomial Horner expansion. Expected
// coefficients themselves come from the independent Cox/VAS/resultant oracle.
fn substitute(p: &[R], offset: &R, scale: &R) -> Vec<R> {
    p.iter().rev().fold(vec![], |v, c| {
        subtract(
            &product(&v, &[offset.clone(), scale.clone()]),
            std::slice::from_ref(c),
            &-one(),
        )
    })
}

#[test]
fn complete_minimum_sets_and_curve_equations_match_independent_oracle() {
    let inputs = include_str!("../../../../fixtures/spline-proximity-inputs.txt");
    let expected = include_str!("../../../../fixtures/spline-proximity.tsv");
    let expected: Vec<_> = expected.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(inputs.lines().count(), expected.len());
    assert_eq!(expected.len(), 30);
    let mut total_points = 0;
    let mut total_intervals = 0;
    for (input, expected) in inputs.lines().zip(expected) {
        let mut data = Words::new(input);
        let name: String = data.take();
        let degree = data.take();
        let periodic: usize = data.take();
        let np = data.take();
        let nk = data.take();
        let first: f64 = data.take();
        let last: f64 = data.take();
        let point = Point3::new(data.take(), data.take(), data.take());
        let mut poles = vec![];
        let mut weights = vec![];
        for _ in 0..np {
            poles.push(Point3::new(data.take(), data.take(), data.take()));
            weights.push(data.take());
        }
        let mut knots = vec![];
        let mut multiplicities = vec![];
        for _ in 0..nk {
            knots.push(data.take());
            multiplicities.push(data.take());
        }
        data.end();
        let curve = if periodic == 0 {
            BSplineCurve3::new(degree, poles, Some(weights), knots, multiplicities)
        } else {
            BSplineCurve3::new_periodic(degree, poles, Some(weights), knots, multiplicities)
        }
        .unwrap();
        let exact = curve.to_exact();
        let mut oracle = Words::new(expected);
        assert_eq!(oracle.take::<String>(), name);
        let cells = oracle.take();
        for _ in 0..cells {
            let lo: R = oracle.take();
            let hi: R = oracle.take();
            let h: [Vec<R>; 4] = std::array::from_fn(|_| oracle.polynomial());
            if lo == hi {
                let actual = exact
                    .exact_evaluate(&lo, DerivativeOrder::Position, KnotSide::Automatic)
                    .unwrap();
                for (i, c) in actual.position().coordinates().iter().enumerate() {
                    assert_eq!(
                        *c,
                        evaluate(&h[i], &zero()) / evaluate(&h[3], &zero()),
                        "{name}"
                    );
                }
            } else {
                let spans = exact.knot_vector().spans_in(&lo, &hi, 4096).unwrap();
                assert_eq!(spans.len(), 1, "{name}");
                let span = &spans[0];
                let length = &span.end - &span.start;
                let offset = (&lo - &span.start) / &length;
                let scale = (&hi - &lo) / &length;
                for (actual, expected) in exact.span_polynomial(span.index).iter().zip(h) {
                    assert_eq!(substitute(actual, &offset, &scale), expected, "{name}");
                }
            }
        }
        let result = closest_points_on_spline_in(&curve, point, first, last, Default::default())
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let points: usize = oracle.take();
        let intervals: usize = oracle.take();
        assert_eq!(result.points().len(), points, "{name}");
        assert_eq!(result.intervals().len(), intervals, "{name}");
        for actual in result.points() {
            let a: R = oracle.take();
            let b: R = oracle.take();
            let lo: R = oracle.take();
            let hi: R = oracle.take();
            let p = oracle.polynomial();
            let offset = (&actual.at.span.start - &a) / (&b - &a);
            let scale = &actual.at.span.length / (&b - &a);
            // Exact membership in the independent irreducible equation AND
            // its unique-root interval identifies the parameter, even when
            // two different parameters have identical binary64 enclosures.
            assert!(
                actual
                    .at
                    .root
                    .vanishes_polynomial(&IntPolynomial::from_rationals(&substitute(
                        &p, &offset, &scale
                    ))),
                "{name}"
            );
            assert_ne!(
                actual.compare_parameter(&(&a + (&b - &a) * lo)).unwrap(),
                Less,
                "{name}"
            );
            assert_ne!(
                actual.compare_parameter(&(&a + (&b - &a) * hi)).unwrap(),
                Greater,
                "{name}"
            );
        }
        for actual in result.intervals() {
            assert_eq!(
                actual.parameters(),
                &[oracle.take(), oracle.take()],
                "{name}"
            );
        }
        oracle.end();
        total_points += points;
        total_intervals += intervals;
    }
    assert_eq!((total_points, total_intervals), (66, 3));
}
