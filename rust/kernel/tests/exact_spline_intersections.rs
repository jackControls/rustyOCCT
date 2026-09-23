use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::intersection::*;
use rusty_occt::{Error, ExactBSplineCurve3, ExactKnotVector, Point3, Vec3};
use std::cmp::Ordering::{Equal, Greater, Less};
#[allow(dead_code)]
#[path = "support/knot_reference.rs"]
mod identity;

fn r(n: i64) -> R {
    R::from_integer(n.into())
}
fn q(n: i64, d: i64) -> R {
    R::new(n.into(), d.into())
}
fn power(n: usize) -> R {
    R::from_integer(BigInt::from(1) << n)
}
fn plane() -> Plane3 {
    Plane3::through_points(
        Point3::new(0., 0., 0.),
        Point3::new(1., 0., 0.),
        Point3::new(0., 1., 0.),
    )
    .unwrap()
}
fn sphere() -> Sphere3 {
    Sphere3::new(Point3::new(0., 0., 0.), 1.).unwrap()
}
fn cylinder() -> Cylinder3 {
    Cylinder3::new(Point3::new(0., 0., 0.), Vec3::new(0., 0., 1.), 1.).unwrap()
}
fn clamped(controls: Vec<[R; 4]>, domain: [R; 2]) -> ExactBSplineCurve3 {
    let n = controls.len();
    ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new(n - 1, domain.to_vec(), vec![n, n]).unwrap(),
        controls,
    )
    .unwrap()
}
fn choose(n: usize, k: usize) -> R {
    R::from_integer((0..k).fold(BigInt::from(1), |a, i| a * (n - i) / (i + 1)))
}
// Independently manufacture a complete known-factor numerator in Bernstein
// form. W(t)=1+t is positive; no production spline or root routine is used.
fn factored(roots: &[R], family: usize, domain: [R; 2]) -> ExactBSplineCurve3 {
    let mut coefficients = vec![r(1)];
    for root in roots {
        let mut next = vec![r(0); coefficients.len() + 1];
        for (i, c) in coefficients.iter().enumerate() {
            next[i] -= c * root;
            next[i + 1] += c;
        }
        coefficients = next;
    }
    let p = roots.len();
    let controls = (0..=p)
        .map(|i| {
            let value: R = (0..=i)
                .map(|j| &coefficients[j] * choose(i, j) / choose(p, j))
                .sum();
            let t = R::new(i.into(), p.into());
            let w = r(1) + &t;
            match family {
                0 => [t, r(0), value, w],
                1 => [value, w.clone(), r(0), w],
                2 => [value, w.clone(), t, w],
                _ => unreachable!(),
            }
        })
        .collect();
    clamped(controls, domain)
}
fn intersect(
    curve: &ExactBSplineCurve3,
    family: usize,
    a: &R,
    b: &R,
) -> ExactSplineSurfaceIntersection {
    match family {
        0 => exact_spline_plane_in(curve, &plane(), a, b),
        1 => exact_spline_sphere_in(curve, &sphere(), a, b),
        2 => exact_spline_cylinder_in(curve, &cylinder(), a, b),
        _ => unreachable!(),
    }
    .unwrap()
}

fn check_bounds(actual: rusty_occt::Result<rusty_occt::ScalarInterval>, expected: &str) {
    if expected == "X" {
        assert!(matches!(actual, Err(Error::Unrepresentable(_))));
    } else {
        let (lo, hi) = expected.split_once(':').unwrap();
        let actual = actual.unwrap();
        assert_eq!(
            actual.lower().to_bits(),
            u64::from_str_radix(lo, 16).unwrap()
        );
        assert_eq!(
            actual.upper().to_bits(),
            u64::from_str_radix(hi, 16).unwrap()
        );
    }
}

#[test]
fn arbitrary_rational_inputs_match_independent_fraction_certificates() {
    let mut count = 0;
    for row in include_str!("../../fixtures/exact-spline-intersections.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let (input, output) = row.split_once(" | ").unwrap();
        let mut w = input.split_whitespace();
        let name = w.next().unwrap();
        let family: usize = w.next().unwrap().parse().unwrap();
        let degree = w.next().unwrap().parse().unwrap();
        let periodic = w.next().unwrap() == "1";
        let np: usize = w.next().unwrap().parse().unwrap();
        let nk: usize = w.next().unwrap().parse().unwrap();
        let first: R = w.next().unwrap().parse().unwrap();
        let last: R = w.next().unwrap().parse().unwrap();
        let controls = (0..np)
            .map(|_| std::array::from_fn(|_| w.next().unwrap().parse().unwrap()))
            .collect();
        let mut knots = Vec::new();
        let mut mults = Vec::new();
        for _ in 0..nk {
            knots.push(w.next().unwrap().parse().unwrap());
            mults.push(w.next().unwrap().parse().unwrap());
        }
        assert!(w.next().is_none());
        let basis = if periodic {
            ExactKnotVector::new_periodic(degree, knots, mults)
        } else {
            ExactKnotVector::new(degree, knots, mults)
        }
        .unwrap();
        let curve = ExactBSplineCurve3::from_homogeneous(basis, controls).unwrap();
        let hits = intersect(&curve, family, &first, &last);
        let mut words = output.split_whitespace();
        assert_eq!(
            hits.points().len(),
            words.next().unwrap().parse::<usize>().unwrap(),
            "{name}"
        );
        assert_eq!(
            hits.overlaps().len(),
            words.next().unwrap().parse::<usize>().unwrap(),
            "{name}"
        );
        let mut representable = true;
        for point in hits.points() {
            let values: [R; 4] = std::array::from_fn(|_| words.next().unwrap().parse().unwrap());
            assert_eq!(
                point.compare_parameter(&values[0]).unwrap(),
                Equal,
                "{name}"
            );
            for c in 0..3 {
                assert_eq!(
                    point.compare_coordinate(c, &values[c + 1]).unwrap(),
                    Equal,
                    "{name}"
                );
            }
            let mut point_representable = true;
            for c in 0..4 {
                let expected = words.next().unwrap();
                point_representable &= expected != "X";
                check_bounds(
                    if c == 0 {
                        point.parameter_bounds()
                    } else {
                        point.coordinate_bound(c - 1)
                    },
                    expected,
                );
            }
            assert_eq!(point.enclosed().is_ok(), point_representable, "{name}");
            representable &= point_representable;
            let contact = match words.next().unwrap() {
                "C" => SplineSurfaceContact::Crossing,
                "T" => SplineSurfaceContact::Tangent,
                "B" => SplineSurfaceContact::Boundary,
                _ => panic!("contact"),
            };
            assert_eq!(point.contact(), contact, "{name}");
            let orders = std::array::from_fn(|_| {
                let n: usize = words.next().unwrap().parse().unwrap();
                (n != 0).then_some(n)
            });
            assert_eq!(point.multiplicities(), orders, "{name}");
        }
        for overlap in hits.overlaps() {
            let values: [R; 2] = std::array::from_fn(|_| words.next().unwrap().parse().unwrap());
            assert_eq!(overlap.parameters(), &values, "{name}");
            let bounds: [&str; 2] = std::array::from_fn(|_| words.next().unwrap());
            let finite = bounds.iter().all(|b| *b != "X");
            representable &= finite;
            if finite {
                let actual = overlap.parameter_bounds().unwrap();
                for i in 0..2 {
                    check_bounds(Ok(actual[i]), bounds[i]);
                }
            } else {
                assert!(matches!(
                    overlap.parameter_bounds(),
                    Err(Error::Unrepresentable(_))
                ));
            }
        }
        assert_eq!(hits.enclosed().is_ok(), representable, "{name}");
        assert!(words.next().is_none(), "{name}");
        count += 1;
    }
    assert_eq!(count, 274);
}

#[test]
fn huge_exact_parameters_and_coordinates_survive_failed_enclosure() {
    let m = power(2048);
    let curve = clamped(
        vec![[r(0), r(0), -r(1), r(1)], [r(0), r(0), r(1), r(1)]],
        [m.clone(), &m + r(2)],
    );
    let result = exact_spline_plane(&curve, &plane()).unwrap();
    assert_eq!(result.points().len(), 1);
    let p = &result.points()[0];
    assert_eq!(p.compare_parameter(&(&m + r(1))).unwrap(), Equal);
    assert!(matches!(
        p.parameter_bounds(),
        Err(Error::Unrepresentable(_))
    ));
    assert!(p
        .coordinate_bounds()
        .unwrap()
        .iter()
        .all(|b| b.lower() == 0. && b.upper() == 0.));
    assert!(matches!(result.enclosed(), Err(Error::Unrepresentable(_))));
    assert_eq!(p.compare_parameter(&m).unwrap(), Greater);

    let curve = clamped(
        vec![
            [m.clone(), r(0), -r(1), r(1)],
            [m.clone(), r(0), r(1), r(1)],
        ],
        [r(0), r(1)],
    );
    let result = exact_spline_plane(&curve, &plane()).unwrap();
    let p = &result.points()[0];
    assert_eq!(p.compare_coordinate(0, &m).unwrap(), Equal);
    assert_eq!(p.parameter_bounds().unwrap().lower(), 0.5);
    assert!(matches!(
        p.coordinate_bounds(),
        Err(Error::Unrepresentable(_))
    ));
    assert!(matches!(p.enclosed(), Err(Error::Unrepresentable(_))));
    // Huge controls do not force a conversion error when the actual hit is finite.
    let curve = clamped(
        vec![[-&m, r(0), -r(1), r(1)], [m, r(0), r(1), r(1)]],
        [r(0), r(1)],
    );
    let result = exact_spline_plane(&curve, &plane())
        .unwrap()
        .enclosed()
        .unwrap();
    assert_eq!(result.points()[0].position(), Point3::new(0., 0., 0.));
}

#[test]
fn rational_trim_keeps_distinct_contacts_inside_one_float_interval() {
    let a = q(1, 3);
    let h = r(1) / power(2048);
    let b = &a + r(2) * &h;
    let curve = ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new(1, vec![a.clone(), &a + &h, b.clone()], vec![2, 1, 2]).unwrap(),
        vec![
            [r(0), r(0), r(1), r(1)],
            [r(0), r(0), -r(1), r(1)],
            [r(0), r(0), r(1), r(1)],
        ],
    )
    .unwrap();
    let hits = exact_spline_plane(&curve, &plane()).unwrap();
    let roots = [&a + &h / r(2), &a + q(3, 2) * &h];
    assert_eq!(hits.points().len(), 2);
    assert_eq!(
        hits.points()[0].parameter_bounds().unwrap(),
        hits.points()[1].parameter_bounds().unwrap()
    );
    for (p, t) in hits.points().iter().zip(&roots) {
        assert_eq!(p.compare_parameter(t).unwrap(), Equal);
        assert_eq!(p.contact(), SplineSurfaceContact::Crossing);
        assert_eq!(p.multiplicities(), [Some(1), Some(1)]);
    }
    assert_eq!(hits.points()[0].compare_parameter(&roots[1]).unwrap(), Less);
    let trimmed =
        exact_spline_plane_in(&curve, &plane(), &(&roots[0] + &h / r(8)), &roots[1]).unwrap();
    assert_eq!(trimmed.points().len(), 1);
    assert_eq!(
        trimmed.points()[0].compare_parameter(&roots[1]).unwrap(),
        Equal
    );
    assert_eq!(
        trimmed.points()[0].contact(),
        SplineSurfaceContact::Boundary
    );
    assert_eq!(trimmed.points()[0].multiplicities(), [Some(1), None]);
}

#[test]
fn exact_rational_circle_produces_whole_overlaps_at_huge_parameters() {
    let a = power(2048);
    let b = &a + r(1);
    // (1-t^2, 2t, 0)/(1+t^2) lies identically on both unit surfaces.
    let curve = clamped(
        vec![
            [r(1), r(0), r(0), r(1)],
            [r(1), r(1), r(0), r(1)],
            [r(0), r(2), r(0), r(2)],
        ],
        [a.clone(), b.clone()],
    );
    let elevated = curve.bezier_arcs().unwrap()[0]
        .elevated(25)
        .unwrap()
        .to_bspline();
    for family in [1, 2] {
        for c in [&curve, &elevated] {
            let lo = &a + q(1, 7);
            let hi = &b - q(1, 11);
            let hits = intersect(c, family, &lo, &hi);
            assert!(hits.points().is_empty());
            assert_eq!(hits.overlaps().len(), 1);
            let interval = &hits.overlaps()[0];
            assert_eq!(interval.parameters(), &[lo.clone(), hi.clone()]);
            assert_eq!(interval.compare_parameter(0, &lo).unwrap(), Equal);
            assert_eq!(interval.compare_parameter(1, &hi).unwrap(), Equal);
            assert!(matches!(
                interval.parameter_bounds(),
                Err(Error::Unrepresentable(_))
            ));
            assert!(matches!(hits.enclosed(), Err(Error::Unrepresentable(_))));
        }
    }
}

#[test]
fn degree25_rational_contacts_and_degree50_quadric_orders_survive_edits() {
    let roots = [
        (r(0), 1),
        (q(1, 7), 11),
        (q(1, 3), 7),
        (q(2, 3), 5),
        (r(1), 1),
    ];
    let factors: Vec<_> = roots
        .iter()
        .flat_map(|(t, m)| std::iter::repeat_n(t.clone(), *m))
        .collect();
    let a = q(-2, 3);
    let width = q(5, 7);
    let b = &a + &width;
    for family in 0..3 {
        let original = factored(&factors, family, [a.clone(), b.clone()]);
        let cut = &a + &width * q(1, 7);
        let fine = original.insert_knot(&cut, 25).unwrap();
        assert!(identity::equal(&original, &fine));
        assert_eq!(fine.remove_knot(&cut, 0).unwrap(), Some(original.clone()));
        for c in [&original, &fine] {
            let hits = intersect(c, family, &a, &b);
            assert_eq!(hits.points().len(), roots.len());
            assert!(hits.overlaps().is_empty());
            for (i, (p, (t, m))) in hits.points().iter().zip(&roots).enumerate() {
                assert_eq!(p.compare_parameter(&(&a + &width * t)).unwrap(), Equal);
                let m = m * if family == 0 { 1 } else { 2 };
                assert_eq!(
                    p.multiplicities(),
                    [(i > 0).then_some(m), (i + 1 < roots.len()).then_some(m)]
                );
                assert_eq!(
                    p.contact(),
                    if i == 0 || i + 1 == roots.len() {
                        SplineSurfaceContact::Boundary
                    } else if family == 0 {
                        SplineSurfaceContact::Crossing
                    } else {
                        SplineSurfaceContact::Tangent
                    }
                );
                let value = t / (r(1) + t);
                let expected = match family {
                    0 => [value, r(0), r(0)],
                    1 => [r(0), r(1), r(0)],
                    _ => [r(0), r(1), value],
                };
                for (c, x) in expected.iter().enumerate() {
                    assert_eq!(p.compare_coordinate(c, x).unwrap(), Equal);
                }
            }
            assert!(hits.enclosed().is_ok());
        }
    }
}

#[test]
fn periodic_origin_removal_preserves_contacts_on_the_same_query() {
    let original = ExactBSplineCurve3::from_homogeneous(
        ExactKnotVector::new_periodic(1, vec![r(0), r(1), r(2), r(3)], vec![1; 4]).unwrap(),
        vec![
            [r(0), r(0), r(0), r(1)],
            [r(1), r(0), r(1), r(1)],
            [-r(1), r(0), -r(1), r(1)],
        ],
    )
    .unwrap();
    let shifted = original.remove_knot(&r(0), 0).unwrap().unwrap();
    assert_eq!(shifted.domain(), &[r(1), r(4)]);
    assert!(identity::equal(&original, &shifted));
    for family in 0..3 {
        let before = intersect(&original, family, &-r(3), &r(6))
            .enclosed()
            .unwrap();
        let after = intersect(&shifted, family, &-r(3), &r(6))
            .enclosed()
            .unwrap();
        assert_eq!(before.points().len(), after.points().len());
        for (a, b) in before.points().iter().zip(after.points()) {
            assert_eq!(a.parameter(), b.parameter());
            assert_eq!(a.coordinate_bounds(), b.coordinate_bounds());
            assert_eq!(a.contact(), b.contact());
            assert_eq!(a.multiplicities(), b.multiplicities());
        }
        assert_eq!(before.overlaps(), after.overlaps());
    }
    let hits = exact_spline_plane_in(&shifted, &plane(), &-r(3), &r(6)).unwrap();
    assert_eq!(hits.points().len(), 7);
    for (i, p) in hits.points().iter().enumerate() {
        assert_eq!(
            p.compare_parameter(&(-r(3) + q(3, 2) * r(i as i64)))
                .unwrap(),
            Equal
        );
    }
    let far = power(2048);
    let a = &far - (&far / r(3)).floor() * r(3);
    let hits = exact_spline_plane_in(&shifted, &plane(), &far, &(&far + r(3))).unwrap();
    assert_eq!(hits.points().len(), 2);
    let first = &far + (q(3, 2) - a);
    for (i, p) in hits.points().iter().enumerate() {
        assert_eq!(
            p.compare_parameter(&(&first + q(3, 2) * r(i as i64)))
                .unwrap(),
            Equal
        );
    }
}

#[test]
fn retained_degree25_weighted_contacts_and_parameter_ranges() {
    use std::collections::BTreeMap;
    let inputs = [
        include_bytes!(
            "../../fuzz/regressions/exact_spline_intersections/degree25-negative-huge-domain.bin"
        ),
        include_bytes!(
            "../../fuzz/regressions/exact_spline_intersections/degree25-positive-huge-domain.bin"
        ),
        include_bytes!(
            "../../fuzz/regressions/exact_spline_intersections/degree25-unit-domain.bin"
        ),
        include_bytes!(
            "../../fuzz/regressions/exact_spline_intersections/mutated-rational-contacts.bin"
        ),
        include_bytes!(
            "../../fuzz/regressions/exact_spline_intersections/degree25-sphere-unit-domain.bin"
        ),
    ];
    let minimal = |interval: rusty_occt::ScalarInterval, value: &R| {
        let lo = R::from_float(interval.lower()).unwrap();
        let hi = R::from_float(interval.upper()).unwrap();
        assert!(lo <= *value && *value <= hi);
        if lo != hi {
            assert!(lo < *value && *value < hi);
            // All finite values in these recipes are nonnegative.
            assert_eq!(interval.upper().to_bits(), interval.lower().to_bits() + 1);
        }
    };
    for bytes in inputs {
        let family = usize::from(bytes[0]);
        let domain_mode = bytes[3];
        assert!(family <= 2);
        assert_eq!(&bytes[1..3], &[24, 0]);
        assert_eq!(&bytes[4..16], &[1, 2, 0, 0, 7, 2, 5, 255, 0, 84, 0, 0]);
        let (a, width) = match domain_mode {
            0 => (q(1, 3), q(5, 7)),
            2 => (-power(2048), q(7, 3)),
            5 => (r(0), r(1)),
            6 => (power(1024), power(1024)),
            _ => panic!("unexpected retained domain"),
        };
        let b = &a + &width;
        let roots: Vec<_> = bytes[16..41]
            .iter()
            .map(|v| q(i64::from(v % 10) - 1, 7))
            .collect();
        let weights: Vec<_> = bytes[64..90]
            .iter()
            .map(|v| q(i64::from(1 + v % 31), 11))
            .collect();
        // Known numerator factors; replace the positive Bernstein weights and
        // rotate coordinates. The plane equation is N; the cylinder equation
        // is N^2. Positive Bernstein weights do not introduce other contacts.
        let base = factored(&roots, family, [a.clone(), b.clone()]);
        let controls = base
            .homogeneous_poles()
            .iter()
            .zip(&weights)
            .map(|(h, w)| {
                if family == 0 {
                    [r(0), h[2].clone(), h[0].clone(), w.clone()]
                } else {
                    [w.clone(), h[2].clone(), h[0].clone(), w.clone()]
                }
            })
            .collect();
        let original =
            ExactBSplineCurve3::from_homogeneous(base.knot_vector().clone(), controls).unwrap();
        let cut = &a + &width * q(85, 257);
        let next = &a + &width * q(171, 257);
        let curve = original
            .refined(&[(cut.clone(), 25), (next, 1)])
            .unwrap()
            .remove_knot(&cut, 0)
            .unwrap()
            .unwrap();
        assert!(identity::equal(&original, &curve));
        let cylinder = Cylinder3::new(Point3::new(0., 0., 0.), Vec3::new(0., 1., 0.), 1.).unwrap();
        let hits = if family == 0 {
            let plane = Plane3::through_points(
                Point3::new(0., 0., 0.),
                Point3::new(0., 0., 1.),
                Point3::new(1., 0., 0.),
            )
            .unwrap();
            exact_spline_plane(&curve, &plane).unwrap()
        } else if family == 1 {
            exact_spline_sphere(&curve, &sphere()).unwrap()
        } else {
            exact_spline_cylinder(&curve, &cylinder).unwrap()
        };
        let mut orders = BTreeMap::new();
        for t in roots.into_iter().filter(|t| t >= &r(0) && t <= &r(1)) {
            *orders.entry(t).or_insert(0) += if family == 0 { 1 } else { 2 };
        }
        assert_eq!(hits.points().len(), orders.len());
        assert!(hits.overlaps().is_empty());
        for (point, (t, order)) in hits.points().iter().zip(orders) {
            let parameter = &a + &width * &t;
            assert_eq!(point.compare_parameter(&parameter).unwrap(), Equal);
            if domain_mode == 2 || domain_mode == 6 {
                assert!(matches!(
                    point.parameter_bounds(),
                    Err(Error::Unrepresentable(_))
                ));
            } else {
                minimal(point.parameter_bounds().unwrap(), &parameter);
            }
            let w: R = weights
                .iter()
                .enumerate()
                .map(|(i, w)| {
                    w * choose(25, i) * t.pow(i as i32) * (r(1) - &t).pow((25 - i) as i32)
                })
                .sum();
            let expected = if family == 0 {
                [r(0), r(0), &t / w]
            } else if family == 1 {
                [r(1), r(0), r(0)]
            } else {
                [r(1), &t / w, r(0)]
            };
            for (i, value) in expected.iter().enumerate() {
                assert_eq!(point.compare_coordinate(i, value).unwrap(), Equal);
                minimal(point.coordinate_bound(i).unwrap(), value);
            }
            assert_eq!(
                point.multiplicities(),
                [(t > r(0)).then_some(order), (t < r(1)).then_some(order)]
            );
            assert_eq!(
                point.contact(),
                if t == r(0) || t == r(1) {
                    SplineSurfaceContact::Boundary
                } else if order % 2 == 1 {
                    SplineSurfaceContact::Crossing
                } else {
                    SplineSurfaceContact::Tangent
                }
            );
        }
    }
}

#[test]
fn an_unrepresentable_contact_does_not_discard_other_exact_contacts() {
    let m = power(2048);
    let roots = [r(1) / &m, q(1, 2)];
    let c = factored(&roots, 0, [r(0), m.clone()]);
    let hits = exact_spline_plane(&c, &plane()).unwrap();
    assert_eq!(hits.points().len(), 2);
    assert_eq!(hits.points()[0].compare_parameter(&r(1)).unwrap(), Equal);
    assert_eq!(
        hits.points()[1].compare_parameter(&(&m / r(2))).unwrap(),
        Equal
    );
    assert!(hits.points()[0].enclosed().is_ok());
    assert!(matches!(
        hits.points()[1].enclosed(),
        Err(Error::Unrepresentable(_))
    ));
    assert!(matches!(hits.enclosed(), Err(Error::Unrepresentable(_))));
    assert_eq!(
        hits.points()[1].compare_coordinate(0, &q(1, 3)).unwrap(),
        Equal
    );
}

#[test]
fn exact_domain_validation_normalization_and_work_limits() {
    let curve = clamped(
        vec![[-r(2), r(0), -r(1), r(1)], [r(2), r(0), r(1), r(1)]],
        [r(0), r(1)],
    );
    let bad = R::new_raw(1.into(), 0.into());
    for (a, b) in [
        (-r(1), r(1)),
        (r(0), r(2)),
        (r(1), r(0)),
        (r(1), r(1)),
        (bad.clone(), r(1)),
        (r(0), bad.clone()),
    ] {
        assert!(exact_spline_plane_in(&curve, &plane(), &a, &b).is_err());
        assert!(exact_spline_sphere_in(&curve, &sphere(), &a, &b).is_err());
        assert!(exact_spline_cylinder_in(&curve, &cylinder(), &a, &b).is_err());
    }
    let options = SplineSurfaceOptions {
        max_spans: 0,
        ..Default::default()
    };
    assert!(matches!(
        exact_spline_plane_with_options(&curve, &plane(), options),
        Err(Error::ComputationLimit(_))
    ));
    assert!(matches!(
        exact_spline_sphere_with_options(&curve, &sphere(), options),
        Err(Error::ComputationLimit(_))
    ));
    assert!(matches!(
        exact_spline_cylinder_with_options(&curve, &cylinder(), options),
        Err(Error::ComputationLimit(_))
    ));
    let options = SplineSurfaceOptions {
        root_isolation: rusty_occt::polynomial::RootIsolationOptions {
            max_subdivisions: 0,
        },
        ..Default::default()
    };
    assert!(matches!(
        exact_spline_sphere_with_options(&curve, &sphere(), options),
        Err(Error::ComputationLimit(_))
    ));
    let hits = exact_spline_plane(&curve, &plane()).unwrap();
    let p = &hits.points()[0];
    assert_eq!(
        p.compare_parameter(&R::new_raw((-2).into(), (-4).into()))
            .unwrap(),
        Equal
    );
    assert!(p.compare_parameter(&bad).is_err());
    assert!(p.compare_coordinate(0, &bad).is_err());
    assert!(p.compare_coordinate(3, &r(0)).is_err());
    let constant = clamped(vec![[r(0), r(0), r(0), r(1)]; 2], [r(0), r(1)]);
    let hits = exact_spline_plane(&constant, &plane()).unwrap();
    let overlap = &hits.overlaps()[0];
    assert!(overlap.compare_parameter(0, &bad).is_err());
    assert!(overlap.compare_parameter(2, &r(0)).is_err());
}
