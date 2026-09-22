//! Tensor coefficient identities, exact quotient jets and operation sequences.
use crate::byte;
use num_rational::BigRational as R;
use rusty_occt::surface::BezierPatchExtractionOptions as Options;
use rusty_occt::{BSplineSurface3, Error, KnotVector, Point3};
#[path = "../../kernel/tests/support/surface_editing_protocol.rs"]
#[allow(dead_code)]
mod protocol;
#[path = "../../kernel/tests/support/surface_editing_reference.rs"]
mod reference;

pub fn check_surface_editing(data: &[u8]) {
    check(data, |_| {});
}
pub fn profile_surface_editing(
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
    let bound = match mode {
        0 => 2,
        1 => 5,
        _ => 25,
    };
    let degrees: [usize; 2] = std::array::from_fn(|i| 1 + byte(data, 1 + i) as usize % bound);
    let kinds: [u8; 2] = std::array::from_fn(|i| byte(data, 3 + i) % 3);
    let scales: [f64; 2] = std::array::from_fn(|i| match byte(data, 10 + i) % 8 {
        0 => f64::from_bits(1),
        1 => 2_f64.powi(-500),
        2 => 2_f64.powi(500),
        _ => 1.,
    });
    let axes: [KnotVector; 2] = std::array::from_fn(|i| {
        let d = degrees[i];
        let (knots, mults) = match kinds[i] {
            0 => (vec![0., 1., 3.], vec![d + 1, 1, d + 1]),
            1 => (vec![0., 1., 3.], vec![1, d, 1]),
            _ => (
                (0..2 * d + 4).map(|k| k as f64).collect(),
                vec![1; 2 * d + 4],
            ),
        };
        let knots = knots.into_iter().map(|k| k * scales[i]).collect();
        if kinds[i] == 1 {
            KnotVector::new_periodic(d, knots, mults)
        } else {
            KnotVector::new(d, knots, mults)
        }
        .unwrap()
    });
    let count = axes[0].pole_count() * axes[1].pole_count();
    let coordinate_scale = if mode == 1 {
        2_f64.powi((byte(data, 9) as i32 - 128) * 7)
    } else {
        1.
    };
    let scalar = |i: usize| {
        if mode == 0 {
            f64::from_bits(u64::from_le_bytes(std::array::from_fn(|j| {
                byte(data, 16 + 8 * i + j)
            })))
        } else {
            f64::from(byte(data, 16 + i) as i8) * coordinate_scale
        }
    };
    let poles: Vec<_> = (0..count)
        .map(|i| Point3::new(scalar(4 * i), scalar(4 * i + 1), scalar(4 * i + 2)))
        .collect();
    let weights: Vec<_> = (0..count)
        .map(|i| {
            if mode == 0 {
                scalar(4 * i + 3)
            } else {
                1. + f64::from(byte(data, 16 + 4 * i + 3) % 8)
            }
        })
        .collect();
    let invalid = poles
        .iter()
        .flat_map(|p| p.to_array())
        .any(|x| !x.is_finite())
        || weights.iter().any(|x| !x.is_finite() || *x <= 0.);
    let ranges: [[f64; 2]; 2] = std::array::from_fn(|i| {
        let axis = &axes[i];
        let (a, b) = axis.domain();
        // High degree tensor inputs query one selected span; every operation
        // still checks every coefficient. Lower degrees also exercise full
        // domains and multiple turns. Full high-degree decompositions have
        // complete fixtures and native comparisons in the main test suite.
        let mut range = if degrees[0] * degrees[1] <= 25 && byte(data, 8) & 128 != 0 {
            [a, b]
        } else {
            let spans: Vec<_> = axis
                .knots()
                .windows(2)
                .filter(|k| k[0] >= a && k[1] <= b)
                .collect();
            let k = spans[byte(data, 12 + i) as usize % spans.len()];
            [k[0], k[1]]
        };
        if axis.is_periodic() {
            if degrees[0] * degrees[1] <= 25 && byte(data, 8) & 64 != 0 {
                range = [a - (b - a), b + (b - a)];
            } else if byte(data, 8) & 32 != 0 {
                range = range.map(|x| x - 2. * (b - a));
            }
        }
        range
    });
    let [u, v] = axes;
    let surface = BSplineSurface3::new(u, v, poles, Some(weights));
    if invalid {
        assert!(matches!(
            surface,
            Err(Error::NonFinite(_) | Error::InvalidSurface(_))
        ));
        return;
    }
    let surface = surface.unwrap();
    let [ua, ub] = ranges[0];
    let [va, vb] = ranges[1];
    mark("input construction");
    let patches = surface.bezier_patches_in(ua, ub, va, vb).unwrap();
    mark("kernel extraction");
    let expected = reference::extract(&surface, [ua, ub, va, vb]);
    mark("oracle extraction");
    assert_eq!(patches.len(), expected.len());
    let controls = patches.len() * (degrees[0] + 1) * (degrees[1] + 1);
    for options in [
        Options {
            max_patches: patches.len() - 1,
            max_controls: controls,
        },
        Options {
            max_patches: patches.len(),
            max_controls: controls - 1,
        },
    ] {
        assert!(matches!(
            surface.bezier_patches_with_options(ua, ub, va, vb, options),
            Err(Error::ComputationLimit(_))
        ));
    }
    let probe: [R; 2] = std::array::from_fn(|i| R::new(byte(data, 14 + i).into(), 255.into()));
    for (a, b) in patches.iter().zip(&expected) {
        reference::check(a, b, [&probe[0], &probe[1]]);
    }
    mark("base coefficient jet bounds");
    let selected = byte(data, 8) as usize % patches.len();
    let p = &patches[selected];
    let op = byte(data, 5) % 11;
    let elevation = std::array::from_fn(|i| (degrees[i] + byte(data, 6 + i) as usize % 4).min(25));
    let actual = protocol::apply(p.clone(), op, elevation);
    mark("kernel edits");
    let references = reference::apply(expected[selected].clone(), op, elevation);
    mark("oracle edits");
    assert_eq!(actual.len(), references.len());
    for (a, b) in actual.iter().zip(&references) {
        match (a, b) {
            (protocol::Item::Patch(a), reference::Item::Patch(b)) => {
                reference::check(a, b, [&probe[0], &probe[1]])
            }
            (protocol::Item::Curve(a), reference::Item::Curve(b)) => {
                reference::curves::check(a, b, &probe[0])
            }
            _ => panic!("inconsistent result types"),
        }
    }
    mark("edited coefficient jet bounds");
    let axis = (byte(data, 9) & 1) as usize;
    let fraction = R::new((1 + u16::from(byte(data, 12 + axis))).into(), 257.into());
    let [a, b] = &p.domain()[axis];
    let cut = a + (b - a) * fraction;
    let [left, right] = if axis == 0 {
        p.split_u_at(&cut)
    } else {
        p.split_v_at(&cut)
    }
    .unwrap();
    let iso = |s: &rusty_occt::ExactBezierSurface3| {
        if axis == 0 {
            s.u_iso(&cut)
        } else {
            s.v_iso(&cut)
        }
        .unwrap()
    };
    assert_eq!(iso(&left), iso(&right));
    for split in [&left, &right] {
        reference::check(
            split,
            &expected[selected].trim(split.domain().clone()),
            [&probe[0], &probe[1]],
        );
    }
    assert_eq!(p.exchanged_uv().exchanged_uv(), *p);
    assert_eq!(p.u_reversed().u_reversed(), *p);
    assert_eq!(p.v_reversed().v_reversed(), *p);
    mark("rational split checks");
    if byte(data, 9) & 2 != 0 {
        let raised = p.elevated(elevation[0], elevation[1]).unwrap();
        let split = if axis == 0 {
            raised.split_u_at(&cut)
        } else {
            raised.split_v_at(&cut)
        }
        .unwrap();
        assert_eq!(
            split,
            [
                left.elevated(elevation[0], elevation[1]).unwrap(),
                right.elevated(elevation[0], elevation[1]).unwrap()
            ]
        );
    }
    mark("commutation checks");
}
