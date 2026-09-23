use rusty_occt::curve::{DerivativeOrder as D, KnotSide as S};
use rusty_occt::{BSplineSurface3, BezierSurface3, Error, KnotVector, Point3};

const PARTIALS: [(usize, usize); 5] = [(1, 0), (0, 1), (2, 0), (0, 2), (1, 1)];
fn number(w: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(w, 16).unwrap())
}
fn side(w: &str) -> S {
    match w {
        "L" => S::Left,
        "R" => S::Right,
        "A" => S::Automatic,
        _ => panic!(),
    }
}

#[test]
fn surface_jets_match_independent_exact_tensor_basis_oracle() {
    let mut count = 0;
    for row in include_str!("../../fixtures/surfaces.tsv")
        .lines()
        .filter(|r| !r.starts_with('#'))
    {
        let words: Vec<_> = row.split_whitespace().collect();
        let name = words[0];
        let [du, dv, nu, nv, ku, kv, pu, pv]: [usize; 8] =
            std::array::from_fn(|i| words[i + 2].parse().unwrap());
        let (u, v) = (number(words[10]), number(words[11]));
        let sides = [side(words[12]), side(words[13])];
        let order = match words[14] {
            "0" => D::Position,
            "1" => D::First,
            "2" => D::Second,
            _ => panic!(),
        };
        let mut data = words[15..].iter();
        let mut poles = Vec::new();
        let mut weights = Vec::new();
        for _ in 0..nu * nv {
            let [x, y, z, w] = std::array::from_fn(|_| number(data.next().unwrap()));
            poles.push(Point3::new(x, y, z));
            weights.push(w);
        }
        let mut axis = |degree, count, periodic| {
            let mut knots = Vec::new();
            let mut mults = Vec::new();
            for _ in 0..count {
                knots.push(number(data.next().unwrap()));
                mults.push(data.next().unwrap().parse().unwrap());
            }
            if periodic {
                KnotVector::new_periodic(degree, knots, mults)
            } else {
                KnotVector::new(degree, knots, mults)
            }
            .unwrap()
        };
        let u_axis = axis(du, ku, pu != 0);
        let v_axis = axis(dv, kv, pv != 0);
        let exact = BSplineSurface3::new(
            u_axis.clone(),
            v_axis.clone(),
            poles.clone(),
            Some(weights.clone()),
        )
        .unwrap()
        .to_exact()
        .evaluate(u, v, order, sides);
        let result = if words[1] == "B" {
            BezierSurface3::new(du, dv, poles, Some(weights))
                .unwrap()
                .evaluate(u, v, order)
        } else {
            BSplineSurface3::new(u_axis, v_axis, poles, Some(weights))
                .unwrap()
                .evaluate(u, v, order, sides)
        };
        assert_eq!(
            result, exact,
            "{name}: exact tensor jet disagrees with original homogeneous de Boor path"
        );
        let status = *data.next().unwrap();
        match result {
            Err(Error::Unrepresentable(_)) => assert_eq!(status, "U", "{name}"),
            Err(Error::OutOfDomain(_)) => assert_eq!(status, "O", "{name}"),
            Err(Error::DiscontinuousDerivative) => assert_eq!(status, "D", "{name}"),
            Err(e) => panic!("{name}: {e}"),
            Ok(value) => {
                assert_eq!(status, "P", "{name}");
                for bounds in std::iter::once(value.position_bounds()).chain(
                    PARTIALS
                        .into_iter()
                        .filter_map(|(a, b)| value.derivative_bounds(a, b)),
                ) {
                    for b in bounds {
                        assert_eq!(b.lower(), number(data.next().unwrap()), "{name}");
                        assert_eq!(b.upper(), number(data.next().unwrap()), "{name}");
                    }
                }
            }
        }
        assert!(data.next().is_none(), "{name}");
        count += 1;
    }
    assert_eq!(count, 791);
}

#[test]
fn bilinear_saddle_has_exact_mixed_derivative_and_original_parameter_units() {
    let poles = vec![
        Point3::ORIGIN,
        Point3::new(0., 1., 0.),
        Point3::new(1., 0., 0.),
        Point3::new(1., 1., 1.),
    ];
    let patch = BezierSurface3::new(1, 1, poles.clone(), None).unwrap();
    let result = patch.evaluate(0.25, 0.75, D::Second).unwrap();
    assert_eq!(result.position(), Point3::new(0.25, 0.75, 0.1875));
    for (partial, expected) in [
        ((1, 0), [1., 0., 0.75]),
        ((0, 1), [0., 1., 0.25]),
        ((2, 0), [0.; 3]),
        ((0, 2), [0.; 3]),
        ((1, 1), [0., 0., 1.]),
    ] {
        assert_eq!(
            result
                .derivative_bounds(partial.0, partial.1)
                .unwrap()
                .map(|b| (b.lower(), b.upper())),
            expected.map(|x| (x, x))
        );
    }
    assert!(result.derivative_bounds(0, 0).is_none());
    assert!(result.derivative_bounds(usize::MAX, 0).is_none());
    let scaled = BSplineSurface3::new(
        KnotVector::new(1, vec![2., 4.], vec![2, 2]).unwrap(),
        KnotVector::new(1, vec![-2., 2.], vec![2, 2]).unwrap(),
        poles,
        None,
    )
    .unwrap()
    .evaluate(2.5, 1., D::Second, [S::Automatic; 2])
    .unwrap();
    assert_eq!(scaled.position_bounds(), result.position_bounds());
    assert_eq!(
        scaled
            .derivative_bounds(1, 1)
            .unwrap()
            .map(|b| b.representative()),
        [0., 0., 0.125]
    );
    assert!(patch
        .evaluate(0.25, 0.75, D::Position)
        .unwrap()
        .derivative_bounds(1, 0)
        .is_none());
    assert!(patch
        .evaluate(0.25, 0.75, D::First)
        .unwrap()
        .derivative_bounds(1, 1)
        .is_none());
}

#[test]
fn a_mixed_derivative_discontinuity_is_not_hidden_by_matching_first_partials() {
    let patch = BSplineSurface3::new(
        KnotVector::new(1, vec![0., 1., 2.], vec![2, 1, 2]).unwrap(),
        KnotVector::new(1, vec![0., 1.], vec![2, 2]).unwrap(),
        [1., 0., 2.]
            .into_iter()
            .flat_map(|a| [Point3::ORIGIN, Point3::new(0., 0., a)])
            .collect(),
        None,
    )
    .unwrap();
    assert!(patch
        .evaluate(1., 0., D::First, [S::Automatic, S::Right])
        .is_ok());
    assert_eq!(
        patch.evaluate(1., 0., D::Second, [S::Automatic, S::Right]),
        Err(Error::DiscontinuousDerivative)
    );
    for (side, value) in [(S::Left, -1.), (S::Right, 2.)] {
        assert_eq!(
            patch
                .evaluate(1., 0., D::Second, [side, S::Right])
                .unwrap()
                .derivative_bounds(1, 1)
                .unwrap()
                .map(|b| b.representative()),
            [0., 0., value]
        );
    }
}

#[test]
fn swapping_parameter_directions_and_rescaling_weights_preserves_the_exact_surface() {
    let u = KnotVector::new_periodic(2, vec![-1., 0., 1., 3.], vec![1, 2, 1, 1]).unwrap();
    let v = KnotVector::new(3, vec![0., 1., 2.], vec![4, 2, 4]).unwrap();
    let (nu, nv) = (u.pole_count(), v.pole_count());
    let poles: Vec<_> = (0..nu * nv)
        .map(|i| Point3::new((i % 7) as f64, (i % 5) as f64, (i % 3) as f64))
        .collect();
    let weights: Vec<_> = (0..nu * nv).map(|i| 1. + (i % 3) as f64).collect();
    let a =
        BSplineSurface3::new(u.clone(), v.clone(), poles.clone(), Some(weights.clone())).unwrap();
    let transpose: Vec<_> = (0..nv)
        .flat_map(|j| (0..nu).map(move |i| i * nv + j))
        .collect();
    let b = BSplineSurface3::new(
        v,
        u,
        transpose.iter().map(|&i| poles[i]).collect(),
        Some(transpose.iter().map(|&i| 8. * weights[i]).collect()),
    )
    .unwrap();
    for (u, v) in [(-1., 0.), (0.375, 0.25), (1., 1.), (7., 2.)] {
        for side in [S::Left, S::Right] {
            let vs = if v == 0. {
                S::Right
            } else if v == 2. {
                S::Left
            } else {
                side
            };
            let x = a.evaluate(u, v, D::Second, [side, vs]).unwrap();
            let y = b.evaluate(v, u, D::Second, [vs, side]).unwrap();
            assert_eq!(x.position_bounds(), y.position_bounds());
            for (i, j) in PARTIALS {
                assert_eq!(x.derivative_bounds(i, j), y.derivative_bounds(j, i));
            }
        }
    }
}

#[test]
fn input_validation_bounds_resources_and_rejects_nonfinite_values() {
    let axis = || KnotVector::new(1, vec![0., 1.], vec![2, 2]).unwrap();
    for degrees in [(0, 1), (1, 26), (usize::MAX, 1)] {
        assert!(BezierSurface3::new(degrees.0, degrees.1, vec![], None).is_err());
    }
    for n in [0, 1, 3, 5] {
        assert!(BSplineSurface3::new(axis(), axis(), vec![Point3::ORIGIN; n], None).is_err());
    }
    for weights in [vec![], vec![1.; 3], vec![1.; 5], vec![0.; 4], vec![-1.; 4]] {
        assert!(
            BSplineSurface3::new(axis(), axis(), vec![Point3::ORIGIN; 4], Some(weights)).is_err()
        );
    }
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            BSplineSurface3::new(axis(), axis(), vec![Point3::ORIGIN; 4], Some(vec![x; 4])),
            Err(Error::NonFinite(_))
        ));
        for c in 0..3 {
            let mut coords = [0.; 3];
            coords[c] = x;
            assert!(matches!(
                BSplineSurface3::new(
                    axis(),
                    axis(),
                    vec![Point3::new(coords[0], coords[1], coords[2]); 4],
                    None
                ),
                Err(Error::NonFinite(_))
            ));
        }
        let s = BSplineSurface3::new(axis(), axis(), vec![Point3::ORIGIN; 4], None).unwrap();
        assert!(matches!(
            s.evaluate(x, 0., D::Position, [S::Automatic; 2]),
            Err(Error::NonFinite(_))
        ));
        assert!(matches!(
            s.evaluate(0., x, D::Position, [S::Automatic; 2]),
            Err(Error::NonFinite(_))
        ));
    }
    let mut mults = vec![1; 65];
    mults[0] = 2;
    mults[64] = 2;
    let large = KnotVector::new(1, (0..65).map(|i| i as f64).collect(), mults).unwrap();
    assert!(matches!(
        BSplineSurface3::new(large.clone(), large, vec![], None),
        Err(Error::LimitExceeded(_))
    ));
    let tiny = BSplineSurface3::new(
        KnotVector::new(1, vec![0., f64::from_bits(1)], vec![2, 2]).unwrap(),
        axis(),
        vec![
            Point3::ORIGIN,
            Point3::ORIGIN,
            Point3::new(1., 0., 0.),
            Point3::new(1., 0., 0.),
        ],
        None,
    )
    .unwrap();
    assert!(tiny
        .evaluate(0., 0., D::Position, [S::Automatic; 2])
        .is_ok());
    assert!(matches!(
        tiny.evaluate(0., 0., D::First, [S::Automatic; 2]),
        Err(Error::Unrepresentable(_))
    ));
}
