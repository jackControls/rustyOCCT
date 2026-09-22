//! Complete coefficient identities and exact jets for extraction/edit sequences.
use crate::byte;
use num_rational::BigRational as R;
use rusty_occt::curve::{BezierExtractionOptions, DerivativeOrder};
use rusty_occt::{BSplineCurve3, Error, Point3};
#[path = "../../kernel/tests/support/bezier_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use crate::bezier_reference as reference;

pub fn check_bezier_editing(data: &[u8]) {
    check(data, |_| {});
}

pub fn profile_bezier_editing(
    data: &[u8],
) -> std::collections::BTreeMap<&'static str, std::time::Duration> {
    let mut times = std::collections::BTreeMap::new();
    let mut at = std::time::Instant::now();
    check(data, |label| {
        *times.entry(label).or_default() += at.elapsed();
        at = std::time::Instant::now();
    });
    times
}

fn check(data: &[u8], mut mark: impl FnMut(&'static str)) {
    let mode = byte(data, 0) % 3;
    let degree = 1 + byte(data, 1) as usize
        % if mode == 0 {
            3
        } else if mode == 1 {
            8
        } else {
            25
        };
    let mult = 1 + byte(data, 2) as usize % degree;
    let kind = byte(data, 3) % 3;
    let op = byte(data, 4) % 6;
    let elevation = (degree + byte(data, 5) as usize % 5).min(25);
    let scale = 2_f64.powi((byte(data, 8) as i32 - 128) * 7);
    let parameter_scale = match byte(data, 9) % 8 {
        0 => f64::from_bits(1),
        1 => 2_f64.powi(-500),
        2 => 2_f64.powi(500),
        _ => 1.,
    };
    let (knots, mults) = match kind {
        0 => (vec![0., 1., 4.], vec![degree + 1, mult, degree + 1]),
        1 => (vec![0., 1., 3., 4.], vec![mult, degree, 1, mult]),
        _ => (
            (0..2 * degree + 4).map(|i| i as f64).collect(),
            vec![1; 2 * degree + 4],
        ),
    };
    let knots = knots.into_iter().map(|k| k * parameter_scale).collect();
    let np = if kind == 1 {
        mults.iter().sum::<usize>() - mults[0]
    } else {
        mults.iter().sum::<usize>() - degree - 1
    };
    let scalar = |i| {
        if mode == 0 {
            f64::from_bits(u64::from_le_bytes(std::array::from_fn(|j| {
                byte(data, 12 + 8 * i + j)
            })))
        } else {
            (byte(data, 12 + i) as i8 as f64) * if mode == 1 { scale } else { 1. }
        }
    };
    let poles: Vec<_> = (0..np)
        .map(|i| Point3::new(scalar(4 * i), scalar(4 * i + 1), scalar(4 * i + 2)))
        .collect();
    let weights: Vec<_> = (0..np)
        .map(|i| {
            if mode == 0 {
                scalar(4 * i + 3)
            } else {
                1. + f64::from(byte(data, 12 + 4 * i + 3) % 8)
            }
        })
        .collect();
    let invalid = poles
        .iter()
        .any(|p| p.to_array().iter().any(|x| !x.is_finite()))
        || weights.iter().any(|x| !x.is_finite() || *x <= 0.);
    let build = if kind == 1 {
        BSplineCurve3::new_periodic
    } else {
        BSplineCurve3::new
    };
    let curve = build(degree, poles, Some(weights), knots, mults);
    if invalid {
        assert!(matches!(
            curve,
            Err(Error::NonFinite(_) | Error::InvalidCurve(_))
        ));
        return;
    }
    let curve = curve.unwrap();
    let (a, b) = curve.domain();
    let (first, last) = if kind == 1 && byte(data, 9) & 8 != 0 {
        (a - (b - a) / 2., b + (b - a) / 2.)
    } else {
        (a, b)
    };
    mark("input construction");
    let exact = curve.bezier_arcs_in(first, last).unwrap();
    mark("kernel extraction");
    let expected = reference::extract(&curve, first, last);
    mark("oracle extraction");
    assert_eq!(exact.len(), expected.len());
    assert!(matches!(
        curve.bezier_arcs_with_options(
            first,
            last,
            BezierExtractionOptions {
                max_arcs: exact.len() - 1
            }
        ),
        Err(Error::ComputationLimit(_))
    ));
    let probe = R::new(byte(data, 11).into(), 255.into());
    mark("limits");
    for (arc, reference) in exact.iter().zip(&expected) {
        reference::check(arc, reference, &probe);
        mark("base coefficient/jet/bounds checks");
        let actual = protocol::apply(arc.clone(), op, elevation);
        mark("kernel edits");
        let expected = reference::apply(reference.clone(), op, elevation);
        mark("oracle edits");
        assert_eq!(actual.len(), expected.len());
        for (a, b) in actual.iter().zip(&expected) {
            reference::check(a, b, &probe);
        }
        mark("edited coefficient/jet/bounds checks");
    }
    // Additional rational cuts vary independently from the native fixed probes.
    let selected = byte(data, 6) as usize % exact.len();
    let arc = &exact[selected];
    let fraction = R::new((1 + u16::from(byte(data, 10))).into(), 257.into());
    let [a, b] = arc.domain();
    let cut = a + (b - a) * fraction;
    let [left, right] = arc.split_at(&cut).unwrap();
    mark("kernel rational split");
    assert_eq!(
        left.homogeneous_poles().last(),
        right.homogeneous_poles().first()
    );
    assert_eq!(
        left.exact_evaluate(&cut, DerivativeOrder::Second).unwrap(),
        right.exact_evaluate(&cut, DerivativeOrder::Second).unwrap()
    );
    reference::check(&left, &expected[selected].trim(a, &cut), &probe);
    reference::check(&right, &expected[selected].trim(&cut, b), &probe);
    assert_eq!(left.reversed().reversed(), left);
    mark("rational split checks");
    if byte(data, 7) & 1 != 0 {
        assert_eq!(
            arc.elevated(elevation).unwrap().split_at(&cut).unwrap(),
            [
                left.elevated(elevation).unwrap(),
                right.elevated(elevation).unwrap()
            ]
        );
    }
    mark("kernel commutation check");
}
