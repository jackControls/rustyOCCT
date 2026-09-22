use rusty_occt::intersection::*;
use rusty_occt::polynomial::{quadratic_roots, QuadraticRoot, QuadraticRoots};
use rusty_occt::{Error, Point3, Result, Vec3};
use std::cmp::Ordering;

fn number(word: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(word, 16).unwrap())
}
fn point(v: &[f64]) -> Point3 {
    Point3::new(v[0], v[1], v[2])
}
fn check_root(root: &QuadraticRoot, words: &mut std::slice::Iter<'_, &str>, label: &str) {
    let multiplicity: u8 = words.next().unwrap().parse().unwrap();
    assert_eq!(root.multiplicity(), multiplicity, "{label}");
    let lower = words.next().unwrap();
    if *lower == "U" {
        assert!(
            matches!(root.bounds(), Err(Error::Unrepresentable(_))),
            "{label}"
        );
        assert!(
            root.compare(-f64::MAX).unwrap() == Ordering::Less
                || root.compare(f64::MAX).unwrap() == Ordering::Greater
        );
    } else {
        let upper = words.next().unwrap();
        let bounds = root.bounds().unwrap();
        assert_eq!(bounds.lower(), number(lower), "{label}");
        assert_eq!(bounds.upper(), number(upper), "{label}");
        assert_ne!(root.compare(bounds.lower()).unwrap(), Ordering::Less);
        assert_ne!(root.compare(bounds.upper()).unwrap(), Ordering::Greater);
        assert!(bounds.representative().is_finite());
    }
}

#[test]
fn quadratic_roots_match_exact_polynomial_sign_oracle() {
    let mut cases = 0;
    for row in include_str!("../../fixtures/quadratic.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let w: Vec<_> = row.split_whitespace().collect();
        let result = quadratic_roots(number(w[1]), number(w[2]), number(w[3])).unwrap();
        let mut values = w[5..].iter();
        match result {
            QuadraticRoots::All => assert_eq!(w[4], "A", "{}", w[0]),
            QuadraticRoots::None => assert_eq!(w[4], "0", "{}", w[0]),
            QuadraticRoots::One(root) => {
                assert_eq!(w[4], "1", "{}", w[0]);
                check_root(&root, &mut values, w[0]);
            }
            QuadraticRoots::Two { lower, upper } => {
                assert_eq!(w[4], "2", "{}", w[0]);
                check_root(&lower, &mut values, w[0]);
                check_root(&upper, &mut values, w[0]);
            }
        }
        assert!(values.next().is_none());
        cases += 1;
    }
    assert_eq!(cases, 448);
}

fn intersection(kind: &str, v: &[f64]) -> Result<CurvedIntersection> {
    let (center, axis, radius, p, q) = (
        point(v),
        Vec3::new(v[3], v[4], v[5]),
        v[6],
        point(&v[7..]),
        point(&v[10..]),
    );
    match kind {
        "s" => line_sphere(p, q, &Sphere3::new(center, radius)?),
        "S" => segment_sphere(p, q, &Sphere3::new(center, radius)?),
        "y" => line_cylinder(p, q, &Cylinder3::new(center, axis, radius)?),
        "Y" => segment_cylinder(p, q, &Cylinder3::new(center, axis, radius)?),
        "c" => line_circle(p, q, &Circle3::new(center, axis, radius)?),
        "C" => segment_circle(p, q, &Circle3::new(center, axis, radius)?),
        _ => panic!("fixture kind"),
    }
}

fn check_hit(hit: CurveHit, expected: &[&str], label: &str) {
    assert_eq!(expected.len(), 9);
    assert_eq!(
        hit.kind,
        match expected[0] {
            "K" => ContactKind::Crossing,
            "T" => ContactKind::Tangent,
            "P" => ContactKind::DegeneratePoint,
            _ => panic!("hit kind"),
        },
        "{label}"
    );
    let bounds = hit.point.bounds();
    let low = bounds.min.to_array();
    let high = bounds.max.to_array();
    let values = [
        hit.point.parameter().lower(),
        hit.point.parameter().upper(),
        low[0],
        high[0],
        low[1],
        high[1],
        low[2],
        high[2],
    ];
    for (a, e) in values.into_iter().zip(&expected[1..]) {
        assert_eq!(a, number(e), "{label}");
    }
    for i in 0..3 {
        let value = hit.point.position().to_array()[i];
        assert!(
            value.is_finite() && low[i] <= value && value <= high[i],
            "{label}"
        );
    }
}

#[test]
fn curved_intersections_match_independent_projection_and_root_oracles() {
    let mut cases = 0;
    for row in include_str!("../../fixtures/curved.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let w: Vec<_> = row.split_whitespace().collect();
        let v: Vec<_> = w[2..15].iter().map(|x| number(x)).collect();
        match intersection(w[0], &v) {
            Err(Error::Degenerate(_)) => assert_eq!(w[15], "D", "{} {}", w[0], w[1]),
            Err(Error::Unrepresentable(_)) => assert_eq!(w[15], "U", "{} {}", w[0], w[1]),
            Err(e) => panic!("{} {}: {e}", w[0], w[1]),
            Ok(CurvedIntersection::Disjoint) => assert_eq!(w[15], "N", "{} {}", w[0], w[1]),
            Ok(CurvedIntersection::Contained) => assert_eq!(w[15], "A", "{} {}", w[0], w[1]),
            Ok(CurvedIntersection::One(hit)) => {
                assert_eq!(w[15], "1", "{} {}", w[0], w[1]);
                check_hit(hit, &w[16..], w[1]);
            }
            Ok(CurvedIntersection::Two { first, second }) => {
                assert_eq!(w[15], "2", "{} {}", w[0], w[1]);
                assert_eq!(w.len(), 34);
                check_hit(first, &w[16..25], w[1]);
                check_hit(second, &w[25..], w[1]);
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 1092);
}

#[test]
fn exact_tangency_is_not_an_epsilon_band() {
    let circle = Circle3::new(Point3::ORIGIN, Vec3::Z, 1.).unwrap();
    let hit = |y| line_circle(Point3::new(-2., y, 0.), Point3::new(2., y, 0.), &circle).unwrap();
    assert!(matches!(
        hit(1.),
        CurvedIntersection::One(CurveHit {
            kind: ContactKind::Tangent,
            ..
        })
    ));
    assert!(matches!(
        hit(f64::from_bits(1.0_f64.to_bits() - 1)),
        CurvedIntersection::Two { .. }
    ));
    assert!(matches!(
        hit(f64::from_bits(1.0_f64.to_bits() + 1)),
        CurvedIntersection::Disjoint
    ));
}

#[test]
fn large_parameter_cancellation_preserves_two_distinct_exact_points() {
    let sphere = Sphere3::new(Point3::ORIGIN, 1.).unwrap();
    let x = 1e100_f64;
    let CurvedIntersection::Two { first, second } = line_sphere(
        Point3::new(x, 0., 0.),
        // Three ulps make the center parameter nonrepresentable; both nearby
        // algebraic roots then share one binary64 enclosure.
        Point3::new(f64::from_bits(x.to_bits() + 3), 0., 0.),
        &sphere,
    )
    .unwrap() else {
        panic!("two roots")
    };
    assert_eq!(first.point.parameter(), second.point.parameter());
    assert_eq!(first.point.position(), Point3::new(-1., 0., 0.));
    assert_eq!(second.point.position(), Point3::new(1., 0., 0.));
    assert_eq!(first.point.bounds().min, first.point.bounds().max);
    assert_eq!(second.point.bounds().min, second.point.bounds().max);
}

#[test]
fn clipped_away_unrepresentable_roots_do_not_discard_a_valid_endpoint() {
    let sphere = Sphere3::new(Point3::new(f64::MAX, 0., 0.), f64::MAX).unwrap();
    let p = Point3::ORIGIN;
    let q = Point3::new(f64::MAX, 0., 0.);
    assert!(matches!(
        line_sphere(p, q, &sphere),
        Err(Error::Unrepresentable(_))
    ));
    let CurvedIntersection::One(hit) = segment_sphere(p, q, &sphere).unwrap() else {
        panic!("endpoint")
    };
    assert_eq!(hit.point.position(), p);
    assert_eq!(hit.point.parameter().lower(), 0.);
    assert_eq!(hit.kind, ContactKind::Crossing);
}

#[test]
fn reversal_and_exact_axis_rescaling_preserve_intersections() {
    let p = Point3::new(-2., 0.5, 0.25);
    let q = Point3::new(2., 0.5, 1.);
    let mut expected = None;
    for scale in [f64::from_bits(1), 2.0_f64.powi(-500), 1., -1., f64::MAX] {
        let cylinder = Cylinder3::new(Point3::ORIGIN, Vec3::new(0., 0., scale), 1.).unwrap();
        let forward = segment_cylinder(p, q, &cylinder).unwrap();
        if let Some(expected) = expected {
            assert_eq!(forward, expected);
        } else {
            expected = Some(forward);
        }
        let CurvedIntersection::Two { first, second } = forward else {
            panic!("two crossings")
        };
        let CurvedIntersection::Two {
            first: back_first,
            second: back_second,
        } = segment_cylinder(q, p, &cylinder).unwrap()
        else {
            panic!("reversed crossings")
        };
        assert_eq!(first.point.bounds(), back_second.point.bounds());
        assert_eq!(second.point.bounds(), back_first.point.bounds());
    }
}

#[test]
fn all_nonfinite_inputs_are_rejected_at_the_public_boundary() {
    for index in 0..3 {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut coefficients = [1., 0., -2.];
            coefficients[index] = invalid;
            assert!(matches!(
                quadratic_roots(coefficients[0], coefficients[1], coefficients[2]),
                Err(Error::NonFinite(_))
            ));
        }
    }
    let QuadraticRoots::Two { lower, .. } = quadratic_roots(1., 0., -2.).unwrap() else {
        panic!()
    };
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(lower.compare(invalid), Err(Error::NonFinite(_))));
    }
    for kind in ["s", "S", "y", "Y", "c", "C"] {
        for index in 0..13 {
            if kind.eq_ignore_ascii_case("s") && (3..6).contains(&index) {
                continue;
            }
            for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut v = [0., 0., 0., 0., 0., 1., 1., -2., 0., 0., 2., 0., 0.];
                v[index] = invalid;
                assert!(
                    matches!(intersection(kind, &v), Err(Error::NonFinite(_))),
                    "{kind} {index}"
                );
            }
        }
    }
}
