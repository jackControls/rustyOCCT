use rusty_occt::intersection::{
    spline_plane, spline_plane_with_options, Plane3, SplinePlaneContact,
};
use rusty_occt::polynomial::RootIsolationOptions;
use rusty_occt::{BSplineCurve3, Error, Point3};
use std::cmp::Ordering::{Equal, Greater, Less};

fn bits(word: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(word, 16).unwrap())
}

#[test]
fn complete_results_match_independent_basis_and_continued_fraction_oracle() {
    let mut count = 0;
    for row in include_str!("../../fixtures/spline-plane.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let (input, expected) = row.split_once(" | ").unwrap();
        let mut w = input.split_whitespace();
        let name = w.next().unwrap();
        let kind = w.next().unwrap();
        let degree = w.next().unwrap().parse().unwrap();
        let np: usize = w.next().unwrap().parse().unwrap();
        let nk: usize = w.next().unwrap().parse().unwrap();
        let mut number = || w.next().unwrap().parse().unwrap();
        let plane: [Point3; 3] = std::array::from_fn(|_| Point3::new(number(), number(), number()));
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for _ in 0..np {
            poles.push(Point3::new(number(), number(), number()));
            weights.push(number());
        }
        let mut knots = Vec::new();
        let mut mults = Vec::new();
        for _ in 0..nk {
            knots.push(w.next().unwrap().parse().unwrap());
            mults.push(w.next().unwrap().parse().unwrap());
        }
        assert!(w.next().is_none());
        let curve = if kind == "P" {
            BSplineCurve3::new_periodic(degree, poles, Some(weights), knots, mults)
        } else {
            BSplineCurve3::new(degree, poles, Some(weights), knots, mults)
        }
        .unwrap();
        let plane = Plane3::through_points(plane[0], plane[1], plane[2]).unwrap();
        let result = spline_plane(&curve, &plane).unwrap_or_else(|e| panic!("{name}: {e}"));
        let mut e = expected.split_whitespace();
        assert_eq!(
            result.points().len(),
            e.next().unwrap().parse::<usize>().unwrap(),
            "{name}"
        );
        assert_eq!(
            result.overlaps().len(),
            e.next().unwrap().parse::<usize>().unwrap(),
            "{name}"
        );
        assert_eq!(
            result.is_disjoint(),
            result.points().is_empty() && result.overlaps().is_empty()
        );
        for point in result.points() {
            for (c, bounds) in std::iter::once(point.parameter())
                .chain(point.coordinate_bounds())
                .enumerate()
            {
                let lo = bits(e.next().unwrap());
                let hi = bits(e.next().unwrap());
                assert_eq!(
                    (bounds.lower(), bounds.upper()),
                    (lo, hi),
                    "{name}, component {c}"
                );
                let compare = |x| {
                    if c == 0 {
                        point.compare_parameter(x)
                    } else {
                        point.compare_coordinate(c - 1, x)
                    }
                    .unwrap()
                };
                assert_eq!(
                    compare(lo),
                    if lo == hi { Equal } else { Greater },
                    "{name}"
                );
                assert_eq!(compare(hi), if lo == hi { Equal } else { Less }, "{name}");
            }
            let contact = match e.next().unwrap() {
                "B" => SplinePlaneContact::Boundary,
                "C" => SplinePlaneContact::Crossing,
                "T" => SplinePlaneContact::Tangent,
                _ => panic!(),
            };
            assert_eq!(point.contact(), contact, "{name}");
            let orders = std::array::from_fn(|_| {
                let n = e.next().unwrap().parse::<usize>().unwrap();
                (n != 0).then_some(n)
            });
            assert_eq!(point.multiplicities(), orders, "{name}");
            assert!(matches!(
                point.compare_coordinate(3, 0.),
                Err(Error::OutOfDomain(_))
            ));
            assert!(matches!(
                point.compare_parameter(f64::NAN),
                Err(Error::NonFinite(_))
            ));
        }
        for overlap in result.overlaps() {
            assert_eq!(
                overlap.parameters(),
                (bits(e.next().unwrap()), bits(e.next().unwrap())),
                "{name}"
            );
        }
        assert!(e.next().is_none(), "{name}");
        count += 1;
    }
    assert_eq!(count, 185);
}

#[test]
fn subdivision_exhaustion_never_publishes_a_partial_intersection() {
    // Two distinct roots: z = 16(t-1/4)(t-3/4).
    let curve = BSplineCurve3::new(
        2,
        vec![
            Point3::new(0., 0., 3.),
            Point3::new(1., 0., -5.),
            Point3::new(2., 0., 3.),
        ],
        None,
        vec![0., 1.],
        vec![3, 3],
    )
    .unwrap();
    let plane = Plane3::through_points(
        Point3::new(0., 0., 0.),
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
    )
    .unwrap();
    assert!(matches!(
        spline_plane_with_options(
            &curve,
            &plane,
            RootIsolationOptions {
                max_subdivisions: 0
            }
        ),
        Err(Error::ComputationLimit(_))
    ));
    assert_eq!(spline_plane(&curve, &plane).unwrap().points().len(), 2);
}
