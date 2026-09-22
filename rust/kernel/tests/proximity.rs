#[path = "support/proximity.rs"]
mod certificate;
use num_rational::BigRational as R;
use rusty_occt::proximity::{closest_points, LinearPrimitive3 as P};
use rusty_occt::{Error, Plane3, Point3, Result, ScalarInterval, Triangle3};
use std::cmp::Ordering;

fn number(word: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(word, 16).unwrap())
}
fn primitive(w: &mut std::str::SplitWhitespace<'_>) -> Result<P> {
    let k = w.next().unwrap();
    let mut p = || {
        Point3::new(
            number(w.next().unwrap()),
            number(w.next().unwrap()),
            number(w.next().unwrap()),
        )
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
fn check_bounds(w: &mut std::str::SplitWhitespace<'_>, bounds: Result<ScalarInterval>) {
    if w.next().unwrap() == "U" {
        assert!(matches!(bounds, Err(Error::Unrepresentable(_))));
    } else {
        let b = bounds.unwrap();
        // Numeric equality deliberately treats both encodings of zero alike.
        assert_eq!(b.lower(), number(w.next().unwrap()));
        assert_eq!(b.upper(), number(w.next().unwrap()));
    }
}

#[test]
fn all_pairings_match_independent_analytic_formulas_and_global_certificates() {
    let mut count = 0;
    for row in include_str!("../../fixtures/proximity.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let mut w = row.split_whitespace();
        let name = w.next().unwrap();
        let a = primitive(&mut w);
        let b = primitive(&mut w);
        let result = match (&a, &b) {
            (Ok(a), Ok(b)) => closest_points(a, b),
            (Err(e), _) | (_, Err(e)) => Err(e.clone()),
        };
        match w.next().unwrap() {
            "D" => assert!(matches!(result, Err(Error::Degenerate(_))), "{name}"),
            "N" => assert!(matches!(result, Err(Error::NonFinite(_))), "{name}"),
            "R" => {
                let pair = result.unwrap_or_else(|e| panic!("{name}: {e}"));
                let expected: R = w.next().unwrap().parse().unwrap();
                assert_eq!(pair.exact_squared_distance(), &expected, "{name}");
                check_bounds(&mut w, pair.squared_distance_bounds());
                check_bounds(&mut w, pair.distance_bounds());
                let (a, b) = (a.unwrap(), b.unwrap());
                certificate::certify(&a, &b, &pair);
                let swapped = closest_points(&b, &a).unwrap();
                assert_eq!(pair.distance_cmp(&swapped), Ordering::Equal, "{name}");
                certificate::certify(&b, &a, &swapped);
                assert_eq!(pair, closest_points(&a, &b).unwrap(), "{name}");
            }
            _ => panic!("fixture status"),
        }
        assert!(w.next().is_none(), "{name}");
        count += 1;
    }
    assert_eq!(count, 554);
}

#[test]
fn invalid_queries_are_typed_and_unrepresentable_outputs_do_not_discard_exact_results() {
    let a = P::Point(Point3::new(1., 0., 0.));
    let b = P::Line([Point3::ORIGIN, Point3::new(f64::from_bits(1), 0., 0.)]);
    let pair = closest_points(&a, &b).unwrap();
    assert_eq!(pair.compare_distance(0.).unwrap(), Ordering::Equal);
    assert_eq!(pair.point(1).unwrap(), Point3::new(1., 0., 0.));
    assert!(matches!(
        pair.parameter_bounds(1),
        Err(Error::Unrepresentable(_))
    ));
    assert_eq!(
        pair.compare_parameter(1, 0, f64::MAX).unwrap(),
        Ordering::Greater
    );
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(pair.compare_distance(x), Err(Error::NonFinite(_))));
        assert!(matches!(
            pair.compare_squared_distance(x),
            Err(Error::NonFinite(_))
        ));
        assert!(matches!(
            pair.compare_coordinate(0, 0, x),
            Err(Error::NonFinite(_))
        ));
        assert!(matches!(
            pair.compare_parameter(1, 0, x),
            Err(Error::NonFinite(_))
        ));
    }
    for i in [2, usize::MAX] {
        assert!(matches!(pair.point(i), Err(Error::OutOfDomain(_))));
        assert!(matches!(pair.point_bounds(i), Err(Error::OutOfDomain(_))));
        assert!(matches!(
            pair.parameter_bounds(i),
            Err(Error::OutOfDomain(_))
        ));
        assert!(matches!(
            pair.compare_coordinate(i, 0, 0.),
            Err(Error::OutOfDomain(_))
        ));
        assert!(matches!(
            pair.compare_parameter(i, 0, 0.),
            Err(Error::OutOfDomain(_))
        ));
    }
    assert!(matches!(
        pair.compare_coordinate(0, 3, 0.),
        Err(Error::OutOfDomain(_))
    ));
    assert!(matches!(
        pair.compare_parameter(0, 0, 0.),
        Err(Error::OutOfDomain(_))
    ));
    assert!(matches!(
        pair.compare_parameter(1, 1, 0.),
        Err(Error::OutOfDomain(_))
    ));
}
