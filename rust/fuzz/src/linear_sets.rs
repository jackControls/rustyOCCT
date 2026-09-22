use super::byte;
use num_rational::BigRational as R;
use rusty_occt::intersection::{
    linear_intersection, LinearIntersection as I, LinearPrimitive3 as P,
};
use rusty_occt::proximity::closest_points;
use rusty_occt::{Error, Plane3, Point3, Triangle3};

#[path = "../../kernel/tests/support/linear_sets.rs"]
mod reference;

fn make(kind: u8, points: &[Point3]) -> rusty_occt::Result<P> {
    Ok(match kind {
        0 => P::Point(points[0]),
        1 => P::Line([points[0], points[1]]),
        2 => P::Segment([points[0], points[1]]),
        3 => P::Plane(Plane3::through_points(points[0], points[1], points[2])?),
        4 => P::Triangle(Triangle3::new(points[0], points[1], points[2])?),
        _ => unreachable!(),
    })
}

fn invalid(kind: u8, points: &[Point3]) -> Option<bool> {
    if points
        .iter()
        .flat_map(|p| p.to_array())
        .any(|x| !x.is_finite())
    {
        return Some(true);
    }
    if kind == 1 && points[0] == points[1] {
        return Some(false);
    }
    if kind >= 3 {
        let p: Vec<_> = points
            .iter()
            .map(|p| p.to_array().map(|x| R::from_float(x).unwrap()))
            .collect();
        let a: [R; 3] = std::array::from_fn(|i| &p[1][i] - &p[0][i]);
        let b: [R; 3] = std::array::from_fn(|i| &p[2][i] - &p[0][i]);
        if (0..3).all(|i| &a[(i + 1) % 3] * &b[(i + 2) % 3] == &a[(i + 2) % 3] * &b[(i + 1) % 3]) {
            return Some(false);
        }
    }
    None
}
fn check_error(e: Error, nonfinite: bool) {
    assert!(if nonfinite {
        matches!(e, Error::NonFinite(_))
    } else {
        matches!(e, Error::Degenerate(_))
    });
}

pub fn check_linear_sets(data: &[u8]) {
    // Mode, two kinds, covariance selector, scale, then up to six xyz points.
    let mode = byte(data, 0) % 4;
    let kinds = [byte(data, 1) % 5, byte(data, 2) % 5];
    let counts = kinds.map(|k| {
        if k == 0 {
            1
        } else if k <= 2 {
            2
        } else {
            3
        }
    });
    let mut p: Vec<_> = (0..counts.iter().sum())
        .map(|i| {
            let a: [f64; 3] = std::array::from_fn(|j| {
                let offset = 5 + (3 * i + j) * 8;
                let word = u64::from_le_bytes(std::array::from_fn(|k| byte(data, offset + k)));
                if mode == 0 {
                    f64::from_bits(word)
                } else {
                    (word as i16 as f64) * 2_f64.powi(byte(data, 4) as i32 - 128)
                }
            });
            Point3::new(a[0], a[1], a[2])
        })
        .collect();
    if mode == 2 {
        for p in &mut p {
            p.z = 0.;
        }
    }
    if mode == 3 {
        // Shared vertices and collapsed segments, including non-unique minima.
        p[counts[0]] = p[0];
        if kinds[0] == 2 {
            p[1] = p[0];
        }
    }
    let vertices = [&p[..counts[0]], &p[counts[0]..]];
    let errors: Vec<_> = (0..2).map(|i| invalid(kinds[i], vertices[i])).collect();
    let mut shapes = Vec::new();
    for i in 0..2 {
        match make(kinds[i], vertices[i]) {
            Ok(p) => shapes.push(p),
            Err(e) => {
                check_error(e, errors[i].expect("constructor rejected valid shape"));
                return;
            }
        }
    }
    let result = linear_intersection(&shapes[0], &shapes[1]);
    if let Some(nonfinite) = errors.into_iter().flatten().next() {
        check_error(result.unwrap_err(), nonfinite);
        return;
    }
    let result = result.unwrap();
    assert_eq!(result, reference::reference(&shapes[0], &shapes[1]));
    reference::check_bounds(&result);
    assert_eq!(linear_intersection(&shapes[1], &shapes[0]).unwrap(), result);
    // Reversal changes line parameters and triangle winding, not the set.
    let reversed: Vec<_> = (0..2)
        .map(|i| {
            let p: Vec<_> = vertices[i].iter().rev().copied().collect();
            make(kinds[i], &p).unwrap()
        })
        .collect();
    assert_eq!(
        linear_intersection(&reversed[0], &reversed[1]).unwrap(),
        result
    );
    // A second kernel capability checks incidence independently of both
    // intersection algorithms. Scaled-integer modes bound the extra cost.
    if mode != 0 && byte(data, 3) % 4 == 0 {
        let distance = closest_points(&shapes[0], &shapes[1]).unwrap();
        assert_eq!(
            distance.exact_squared_distance() == &R::from_integer(0.into()),
            !matches!(result, I::Empty)
        );
    }
}
