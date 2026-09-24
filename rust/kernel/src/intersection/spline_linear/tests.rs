use super::*;
use std::{fmt::Debug, str::FromStr};

struct Words<'a>(std::str::SplitWhitespace<'a>);
impl<'a> Words<'a> {
    fn new(row: &'a str) -> Self {
        Self(row.split_whitespace())
    }
    fn take<T: FromStr>(&mut self) -> T
    where
        T::Err: Debug,
    {
        self.0.next().unwrap().parse().unwrap()
    }
    fn polynomial(&mut self) -> IntPolynomial {
        let n = self.take();
        IntPolynomial::from_rationals(&(0..n).map(|_| self.take()).collect::<Vec<R>>())
    }
    fn end(&mut self) {
        assert!(self.0.next().is_none());
    }
}

#[test]
fn affine_parameter_ordering_matches_global_irreducible_root_oracle() {
    let mut count = 0;
    for row in include_str!("../../../../fixtures/affine-parameters.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let mut words = Words::new(row);
        let name: String = words.take();
        let mut equation = || {
            let n = words.take();
            let lo = words.take();
            let hi = words.take();
            let offset: R = words.take();
            let width: R = words.take();
            let p: Vec<R> = (0..n).map(|_| words.take()).collect();
            let expected: usize = words.take();
            let roots = real::isolate(
                &IntPolynomial::from_rationals(&p),
                lo,
                hi,
                &mut Budget::new(Default::default()),
            )
            .unwrap();
            assert_eq!(roots.len(), expected, "{name}");
            (roots, offset, width)
        };
        let (a, ao, aw) = equation();
        let (b, bo, bw) = equation();
        for x in &a {
            for y in &b {
                let expected: i32 = words.take();
                assert_eq!(
                    x.compare_affine(&ao, &aw, y, &bo, &bw),
                    expected.cmp(&0),
                    "{name}"
                );
                assert_eq!(
                    y.compare_affine(&bo, &bw, x, &ao, &aw),
                    expected.cmp(&0).reverse(),
                    "{name} reversed"
                );
                assert_eq!(x.compare_affine(&ao, &aw, x, &ao, &aw), Equal);
                count += 1;
            }
        }
        words.end();
    }
    assert_eq!(count, 851);
}

#[test]
fn closed_clipping_matches_independent_inequality_sets() {
    let mut totals = (0, 0, 0);
    for row in include_str!("../../../../fixtures/closed-clip.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let mut words = Words::new(row);
        let name: String = words.take();
        let p = words.polynomial();
        let q = words.polynomial();
        let lo = words.take();
        let hi = words.take();
        let np: usize = words.take();
        let ni: usize = words.take();
        let actual = closed_clip(&p, &q, lo, hi, &mut Context::new(Default::default())).unwrap();
        assert_eq!(actual.points.len(), np, "{name}: point count");
        assert_eq!(actual.intervals.len(), ni, "{name}: interval count");
        let mut endpoint = |root: &AlgebraicRoot| {
            let equation = words.polynomial();
            let lo: R = words.take();
            let hi: R = words.take();
            assert_eq!(root.sign_polynomial(&equation), Equal, "{name}: equation");
            assert_ne!(root.compare_rational(&lo), Less, "{name}: lower isolator");
            assert_ne!(
                root.compare_rational(&hi),
                Greater,
                "{name}: upper isolator"
            );
        };
        for point in &actual.points {
            endpoint(point);
        }
        for [lo, hi] in &actual.intervals {
            endpoint(lo);
            endpoint(hi);
        }
        words.end();
        totals.0 += 1;
        totals.1 += np;
        totals.2 += ni;
    }
    assert_eq!(totals, (79, 79, 69));
}

// p(offset + scale*x) by Horner expansion. The fixture's irreducible equation
// and isolator come from the independent Cox/VAS oracle.
fn substitute(p: &[R], offset: &R, scale: &R) -> Vec<R> {
    p.iter().rev().fold(vec![], |v, c| {
        let mut next = add_scaled(std::slice::from_ref(c), &scaled(&v, offset), &one());
        next.push(zero());
        for (i, x) in v.iter().enumerate() {
            next[i + 1] += scale * x;
        }
        next
    })
}

fn binary(value: &R) -> f64 {
    let bound = interval::enclose(|x| value.cmp(&real::rat(x)), "fixture value").unwrap();
    assert_eq!(
        bound.lower(),
        bound.upper(),
        "fixture value is not binary64"
    );
    bound.lower()
}

fn same(actual: &SplineLinearIntersection, expected: &SplineLinearIntersection, name: &str) {
    assert_eq!(actual.points().len(), expected.points().len(), "{name}");
    assert_eq!(actual.overlaps().len(), expected.overlaps().len(), "{name}");
    for (x, y) in actual.points().iter().zip(expected.points()) {
        assert_eq!(x.parameter_cmp(y), Equal, "{name}");
    }
    for (x, y) in actual.overlaps().iter().zip(expected.overlaps()) {
        for (a, b) in x.endpoints().iter().zip(y.endpoints()) {
            assert_eq!(a.parameter_cmp(b), Equal, "{name}");
        }
    }
}

#[test]
fn complete_preimages_match_independent_oracle_and_survive_exact_edits() {
    let inputs = include_str!("../../../../fixtures/spline-linear-inputs.txt");
    let expected = include_str!("../../../../fixtures/spline-linear.tsv");
    let data = |text: &'static str| text.lines().filter(|r| !r.starts_with('#'));
    assert_eq!(data(inputs).count(), data(expected).count());
    let mut totals = (0, 0, 0, 0);
    for (input, expected) in data(inputs).zip(data(expected)) {
        let mut words = Words::new(input);
        let name: String = words.take();
        let kind: String = words.take();
        let is_binary = words.take::<u8>() == 1;
        let degree: usize = words.take();
        let periodic = words.take::<u8>() == 1;
        let np: usize = words.take();
        let nk: usize = words.take();
        let first: R = words.take();
        let last: R = words.take();
        let a: [R; 3] = std::array::from_fn(|_| words.take());
        let b: [R; 3] = std::array::from_fn(|_| words.take());
        let controls: Vec<[R; 4]> = (0..np)
            .map(|_| std::array::from_fn(|_| words.take()))
            .collect();
        let (knots, multiplicities): (Vec<R>, Vec<usize>) = (0..nk)
            .map(|_| (words.take::<R>(), words.take::<usize>()))
            .unzip();
        words.end();
        let homogeneous = controls
            .iter()
            .map(|[x, y, z, w]| [x * w, y * w, z * w, w.clone()])
            .collect();
        let basis = if periodic {
            crate::ExactKnotVector::new_periodic(degree, knots.clone(), multiplicities.clone())
        } else {
            crate::ExactKnotVector::new(degree, knots.clone(), multiplicities.clone())
        }
        .unwrap();
        let curve = ExactBSplineCurve3::from_homogeneous(basis, homogeneous).unwrap();
        let query = |curve: &ExactBSplineCurve3| {
            if kind == "L" {
                exact_spline_line_in(curve, &a, &b, &first, &last)
            } else {
                exact_spline_segment_in(curve, &a, &b, &first, &last)
            }
            .unwrap_or_else(|e| panic!("{name}: {e}"))
        };
        let result = query(&curve);

        let mut oracle = Words::new(expected);
        assert_eq!(oracle.take::<String>(), name);
        let points: usize = oracle.take();
        let intervals: usize = oracle.take();
        assert_eq!(result.points().len(), points, "{name}: points");
        assert_eq!(result.overlaps().len(), intervals, "{name}: intervals");
        let mut check = |actual: &SplineLinearPoint| {
            let first: R = oracle.take();
            let last: R = oracle.take();
            let p = oracle.polynomial();
            let lo: R = oracle.take();
            let hi: R = oracle.take();
            let s_lo: R = oracle.take();
            let s_hi: R = oracle.take();
            let width = &last - &first;
            let offset = (&actual.span.start - &first) / &width;
            let scale = &actual.span.length / &width;
            let p: Vec<R> = p.0.iter().cloned().map(R::from_integer).collect();
            let local = substitute(&p, &offset, &scale);
            // Membership in the independent irreducible equation AND its
            // unique-root isolator identifies the parameter exactly.
            assert!(
                actual
                    .root
                    .vanishes_polynomial(&IntPolynomial::from_rationals(&local)),
                "{name}: equation"
            );
            for (bound, forbidden) in [(&lo, Less), (&hi, Greater)] {
                let u = &first + &width * bound;
                assert_ne!(actual.compare_parameter(&u).unwrap(), forbidden, "{name}");
            }
            assert_ne!(
                actual.compare_linear_parameter(&s_lo).unwrap(),
                Less,
                "{name}"
            );
            assert_ne!(
                actual.compare_linear_parameter(&s_hi).unwrap(),
                Greater,
                "{name}"
            );
        };
        for point in result.points() {
            check(point);
        }
        for overlap in result.overlaps() {
            overlap.endpoints().iter().for_each(&mut check);
        }
        oracle.end();

        if is_binary {
            let poles = controls
                .iter()
                .map(|c| Point3::new(binary(&c[0]), binary(&c[1]), binary(&c[2])))
                .collect();
            let weights = Some(controls.iter().map(|c| binary(&c[3])).collect());
            let knots = knots.iter().map(binary).collect();
            let curve = if periodic {
                BSplineCurve3::new_periodic(degree, poles, weights, knots, multiplicities)
            } else {
                BSplineCurve3::new(degree, poles, weights, knots, multiplicities)
            }
            .unwrap();
            let point = |p: &[R; 3]| Point3::new(binary(&p[0]), binary(&p[1]), binary(&p[2]));
            let (first, last) = (binary(&first), binary(&last));
            let actual = if kind == "L" {
                spline_line_in(&curve, point(&a), point(&b), first, last)
            } else {
                spline_segment_in(&curve, point(&a), point(&b), first, last)
            }
            .unwrap();
            same(&actual, &result, &name);
            totals.2 += 1;
        }

        // Exact knot insertion and elevation preserve the parameter function.
        let [lo, hi] = curve.domain().clone();
        let u = &lo + (&hi - &lo) * R::new(3.into(), 7.into());
        let mut edited = curve.insert_knot(&u, degree).unwrap();
        if degree < 25 {
            edited = edited.elevated(degree + 1).unwrap();
        }
        same(&query(&edited), &result, &name);
        totals.0 += points;
        totals.1 += intervals;
        totals.3 += 1;
    }
    assert_eq!(totals, (92, 20, 43, 60));
}
