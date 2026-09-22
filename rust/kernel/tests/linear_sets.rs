#[path = "support/linear_sets.rs"]
mod reference;
use num_rational::BigRational as R;
use rusty_occt::intersection::{
    linear_intersection, ExactPoint3, LinearIntersection as I, LinearPrimitive3 as P,
};
use rusty_occt::{Error, Plane3, Point3, Result, Triangle3};

fn primitive(w: &mut std::str::SplitWhitespace<'_>) -> Result<P> {
    let k = w.next().unwrap();
    let mut p = || {
        let mut x = || f64::from_bits(u64::from_str_radix(w.next().unwrap(), 16).unwrap());
        Point3::new(x(), x(), x())
    };
    Ok(match k {
        "P" => P::Point(p()),
        "L" => P::Line([p(), p()]),
        "S" => P::Segment([p(), p()]),
        "F" => P::Plane(Plane3::through_points(p(), p(), p())?),
        "T" => P::Triangle(Triangle3::new(p(), p(), p())?),
        _ => panic!("kind"),
    })
}
fn exact(w: &mut std::str::SplitWhitespace<'_>) -> I {
    let (kind, n) = (
        w.next().unwrap(),
        w.next().unwrap().parse::<usize>().unwrap(),
    );
    let mut vector = || std::array::from_fn(|_| w.next().unwrap().parse::<R>().unwrap());
    match kind {
        "E" => {
            assert_eq!(n, 0);
            I::Empty
        }
        "P" => {
            assert_eq!(n, 1);
            I::Point(ExactPoint3::from_coordinates(vector()))
        }
        "S" => {
            assert_eq!(n, 2);
            I::Segment([
                ExactPoint3::from_coordinates(vector()),
                ExactPoint3::from_coordinates(vector()),
            ])
        }
        "G" => {
            assert!((3..=6).contains(&n));
            I::Polygon(
                (0..n)
                    .map(|_| ExactPoint3::from_coordinates(vector()))
                    .collect(),
            )
        }
        "L" => {
            assert_eq!(n, 2);
            I::Line {
                origin: ExactPoint3::from_coordinates(vector()),
                direction: vector(),
            }
        }
        "F" => {
            assert_eq!(n, 2);
            I::Plane {
                normal: vector(),
                offset: w.next().unwrap().parse().unwrap(),
            }
        }
        _ => panic!("result kind"),
    }
}
fn permutations(a: &P) -> Vec<P> {
    match a {
        P::Point(_) => vec![a.clone()],
        P::Line(p) => vec![a.clone(), P::Line([p[1], p[0]])],
        P::Segment(p) => vec![a.clone(), P::Segment([p[1], p[0]])],
        P::Plane(_) | P::Triangle(_) => {
            let (plane, p) = match a {
                P::Plane(p) => (true, p.defining_points()),
                P::Triangle(p) => (false, p.vertices()),
                _ => unreachable!(),
            };
            [
                [0, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ]
            .map(|i| {
                if plane {
                    P::Plane(Plane3::through_points(p[i[0]], p[i[1]], p[i[2]]).unwrap())
                } else {
                    P::Triangle(Triangle3::new(p[i[0]], p[i[1]], p[i[2]]).unwrap())
                }
            })
            .to_vec()
        }
    }
}

#[test]
fn complete_sets_match_independent_boundary_oracle_and_all_defining_point_orders() {
    let mut count = 0;
    let mut dimensions = [0; 4];
    for row in include_str!("../../fixtures/linear_sets.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let mut w = row.split_whitespace();
        let name = w.next().unwrap();
        let (a, b) = (primitive(&mut w), primitive(&mut w));
        let result = match (&a, &b) {
            (Ok(a), Ok(b)) => linear_intersection(a, b),
            (Err(e), _) | (_, Err(e)) => Err(e.clone()),
        };
        match w.next().unwrap() {
            "N" => assert!(matches!(result, Err(Error::NonFinite(_))), "{name}"),
            "D" => assert!(matches!(result, Err(Error::Degenerate(_))), "{name}"),
            "R" => {
                let expected = exact(&mut w);
                let actual = result.unwrap_or_else(|e| panic!("{name}: {e}"));
                assert_eq!(actual, expected, "{name}");
                dimensions[actual.dimension().map_or(0, |x| x + 1)] += 1;
                let (a, b) = (a.unwrap(), b.unwrap());
                assert_eq!(
                    reference::reference(&a, &b),
                    actual,
                    "{name} independent Rust boundary oracle"
                );
                reference::check_bounds(&actual);
                assert_eq!(
                    linear_intersection(&b, &a).unwrap(),
                    actual,
                    "{name} swapped"
                );
                // Full permutations on both sides; separate loops avoid redundant
                // 36 identical solves while checking every input ordering.
                for p in permutations(&a) {
                    assert_eq!(
                        linear_intersection(&p, &b).unwrap(),
                        actual,
                        "{name} first permutation"
                    );
                }
                for p in permutations(&b) {
                    assert_eq!(
                        linear_intersection(&a, &p).unwrap(),
                        actual,
                        "{name} second permutation"
                    );
                }
            }
            _ => panic!("status"),
        }
        assert!(w.next().is_none(), "{name}");
        count += 1;
    }
    assert_eq!(count, 684);
    assert!(dimensions.iter().all(|n| *n > 0));
}

#[test]
fn exact_unrepresentable_intersections_survive_finite_conversion_errors() {
    let a = P::Line([Point3::ORIGIN, Point3::new(1., f64::from_bits(1), 0.)]);
    let b = P::Line([Point3::new(0., 1., 0.), Point3::new(1., 1., 0.)]);
    let I::Point(p) = linear_intersection(&a, &b).unwrap() else {
        panic!("point")
    };
    assert!(p.coordinates()[0] > R::from_float(f64::MAX).unwrap());
    assert!(matches!(p.bounds(), Err(Error::Unrepresentable(_))));
    assert!(matches!(p.position(), Err(Error::Unrepresentable(_))));
    assert_eq!(p.coordinate_bounds(1).unwrap().lower(), 1.);
    for axis in [3, usize::MAX] {
        assert!(matches!(
            p.coordinate_bounds(axis),
            Err(Error::OutOfDomain(_))
        ));
    }
    assert_eq!(linear_intersection(&b, &a).unwrap(), I::Point(p));
}
