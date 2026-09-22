//! Independent Cox basis polynomial oracle. Edits use binomial substitutions;
//! returned Bernstein controls are expanded and compared coefficient by coefficient.
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::curve::DerivativeOrder;
use rusty_occt::{BSplineCurve3, Error, ExactBezierCurve3, ScalarInterval};
pub fn integer(n: usize) -> R {
    R::from_integer(n.into())
}
fn rat(x: f64) -> R {
    R::from_float(x).unwrap()
}
pub fn binomial(n: usize, k: usize) -> u64 {
    let mut x = 1_u64;
    for i in 0..k {
        x = x * (n - i) as u64 / (i + 1) as u64;
    }
    x
}

// Clear coefficient denominators before a linear polynomial transform. This
// preserves the independent power-basis algorithm while avoiding thousands of
// repeated rational GCD reductions under the sanitizer.
pub fn integers(p: &[R]) -> (Vec<BigInt>, BigInt) {
    let mut denominator = BigInt::from(1);
    for x in p {
        let (mut a, mut b) = (denominator.clone(), x.denom().clone());
        while b != BigInt::from(0) {
            (a, b) = (b.clone(), a % b);
        }
        denominator = denominator / a * x.denom();
    }
    (
        p.iter()
            .map(|x| x.numer() * (&denominator / x.denom()))
            .collect(),
        denominator,
    )
}
fn powers(x: &BigInt, count: usize) -> Vec<BigInt> {
    let mut result = vec![BigInt::from(1); count];
    for i in 1..count {
        result[i] = &result[i - 1] * x;
    }
    result
}
pub fn add(a: &[R], b: &[R]) -> Vec<R> {
    let mut result = vec![integer(0); a.len().max(b.len())];
    for p in [a, b] {
        for (i, x) in p.iter().enumerate() {
            result[i] += x;
        }
    }
    while result.last() == Some(&integer(0)) {
        result.pop();
    }
    result
}
pub fn linear(p: &[R], a: &R, b: &R) -> Vec<R> {
    if p.is_empty() {
        return Vec::new();
    }
    let mut result = vec![integer(0); p.len() + 1];
    for (i, x) in p.iter().enumerate() {
        result[i] += a * x;
        result[i + 1] += b * x;
    }
    result
}
pub fn substitute(p: &[R], a: &R, b: &R) -> Vec<R> {
    let (p, den) = integers(p);
    let (an, ad, bn, bd) = (
        powers(a.numer(), p.len()),
        powers(a.denom(), p.len()),
        powers(b.numer(), p.len()),
        powers(b.denom(), p.len()),
    );
    (0..p.len())
        .map(|i| {
            let sum: BigInt = (i..p.len())
                .map(|j| &p[j] * binomial(j, i) * &an[j - i] * &ad[p.len() - 1 - j])
                .sum();
            R::new(sum * &bn[i], &den * &ad[p.len() - 1 - i] * &bd[i])
        })
        .collect()
}
pub fn value(p: &[R], t: &R, order: usize) -> R {
    if order >= p.len() {
        return integer(0);
    }
    let (p, den) = integers(p);
    let (tn, td) = (powers(t.numer(), p.len()), powers(t.denom(), p.len()));
    let sum: BigInt = (order..p.len())
        .map(|i| {
            let factor: usize = (i + 1 - order..=i).product();
            &p[i] * factor * &tn[i - order] * &td[p.len() - 1 - i]
        })
        .sum();
    R::new(sum, den * &td[p.len() - 1 - order])
}
pub fn coefficients(curve: &ExactBezierCurve3) -> [Vec<R>; 4] {
    let p = curve.degree();
    std::array::from_fn(|c| {
        let (controls, den) = integers(
            &curve
                .homogeneous_poles()
                .iter()
                .map(|p| p[c].clone())
                .collect::<Vec<_>>(),
        );
        (0..=p)
            .map(|k| {
                let sum: BigInt = (0..=k)
                    .map(|i| {
                        let v = &controls[i] * (binomial(p, i) * binomial(p - i, k - i));
                        if (k - i) % 2 == 0 {
                            v
                        } else {
                            -v
                        }
                    })
                    .sum();
                R::new(sum, den.clone())
            })
            .collect()
    })
}

#[derive(Clone)]
pub struct Arc {
    pub degree: usize,
    pub domain: [R; 2],
    pub coefficients: [Vec<R>; 4],
}
impl Arc {
    pub fn trim(&self, a: &R, b: &R) -> Self {
        let length = &self.domain[1] - &self.domain[0];
        let lo = (a - &self.domain[0]) / &length;
        let scale = (b - a) / length;
        Self {
            degree: self.degree,
            domain: [a.clone(), b.clone()],
            coefficients: std::array::from_fn(|c| substitute(&self.coefficients[c], &lo, &scale)),
        }
    }
    pub fn reversed(&self) -> Self {
        Self {
            coefficients: std::array::from_fn(|c| {
                substitute(&self.coefficients[c], &integer(1), &-integer(1))
            }),
            ..self.clone()
        }
    }
    pub fn jet(&self, u: &R) -> [[R; 3]; 3] {
        let length = &self.domain[1] - &self.domain[0];
        let t = (u - &self.domain[0]) / &length;
        let h: [[R; 4]; 3] = std::array::from_fn(|r| {
            std::array::from_fn(|c| value(&self.coefficients[c], &t, r) / length.pow(r as i32))
        });
        let w = &h[0][3];
        std::array::from_fn(|r| {
            std::array::from_fn(|c| match r {
                0 => &h[0][c] / w,
                1 => (&h[1][c] * w - &h[0][c] * &h[1][3]) / (w * w),
                _ => {
                    (&h[2][c] * w * w
                        - &h[0][c] * &h[2][3] * w
                        - integer(2) * &h[1][c] * w * &h[1][3]
                        + integer(2) * &h[0][c] * &h[1][3] * &h[1][3])
                        / (w * w * w)
                }
            })
        })
    }
}

pub fn extract(curve: &BSplineCurve3, first: f64, last: f64) -> Vec<Arc> {
    let p = curve.degree();
    let base: Vec<_> = curve
        .knots()
        .iter()
        .zip(curve.multiplicities())
        .flat_map(|(&k, &m)| std::iter::repeat_n(rat(k), m))
        .collect();
    let (start, end) = curve.domain();
    let (start, end, first, last) = (rat(start), rat(end), rat(first), rat(last));
    let period = &end - &start;
    let flat = if curve.is_periodic() {
        let n = (base.len() - curve.multiplicities()[0]) as isize;
        let pad = (p + 1 - curve.multiplicities()[0]) as isize;
        (-pad..n + curve.multiplicities()[0] as isize + pad)
            .map(|i| {
                &base[i.rem_euclid(n) as usize]
                    + R::from_integer(BigInt::from(i.div_euclid(n))) * &period
            })
            .collect()
    } else {
        base
    };
    let mut offset = if curve.is_periodic() {
        ((&first - &start) / &period).floor() * &period
    } else {
        integer(0)
    };
    let mut result = Vec::new();
    loop {
        for span in p..flat.len() - p - 1 {
            let (a, b) = (&flat[span], &flat[span + 1]);
            if a >= b || a < &start || b > &end {
                continue;
            }
            let lo = (a + &offset).max(first.clone());
            let hi = (b + &offset).min(last.clone());
            if lo >= hi {
                continue;
            }
            let mut basis: Vec<Vec<R>> = (0..flat.len() - 1)
                .map(|i| {
                    if i == span {
                        vec![integer(1)]
                    } else {
                        Vec::new()
                    }
                })
                .collect();
            for degree in 1..=p {
                basis = (0..basis.len() - 1)
                    .map(|i| {
                        let left = &flat[i + degree] - &flat[i];
                        let right = &flat[i + degree + 1] - &flat[i + 1];
                        let x = if left == integer(0) || basis[i].is_empty() {
                            Vec::new()
                        } else {
                            linear(&basis[i], &((a - &flat[i]) / &left), &((b - a) / &left))
                        };
                        let y = if right == integer(0) || basis[i + 1].is_empty() {
                            Vec::new()
                        } else {
                            linear(
                                &basis[i + 1],
                                &((&flat[i + degree + 1] - a) / &right),
                                &((a - b) / &right),
                            )
                        };
                        add(&x, &y)
                    })
                    .collect();
            }
            let coefficients = std::array::from_fn(|c| {
                let mut sum = vec![integer(0); p + 1];
                for (i, coefficients) in basis.iter().enumerate() {
                    if coefficients.is_empty() {
                        continue;
                    }
                    let j = i % curve.poles().len();
                    let h = rat(curve.weights()[j])
                        * if c == 3 {
                            integer(1)
                        } else {
                            rat(curve.poles()[j].to_array()[c])
                        };
                    for (k, x) in coefficients.iter().enumerate() {
                        sum[k] += &h * x;
                    }
                }
                sum
            });
            let arc = Arc {
                degree: p,
                domain: [a + &offset, b + &offset],
                coefficients,
            };
            result.push(arc.trim(&lo, &hi));
        }
        if !curve.is_periodic() || &end + &offset >= last {
            break;
        }
        offset += &period;
    }
    result
}

pub fn apply(mut arc: Arc, op: u8, elevation: usize) -> Vec<Arc> {
    if op == 1 || op == 5 {
        let [a, b] = &arc.domain;
        arc = arc.trim(
            &(a + (b - a) / integer(4)),
            &(a + (b - a) * R::new(3.into(), 4.into())),
        );
    }
    if op == 4 || op == 5 {
        arc.degree = elevation;
    }
    if op == 3 || op == 5 {
        arc = arc.reversed();
    }
    if op == 2 || op == 5 {
        let [a, b] = &arc.domain;
        let cut = a + (b - a) * R::new(3.into(), 8.into());
        vec![arc.trim(a, &cut), arc.trim(&cut, b)]
    } else {
        vec![arc]
    }
}

pub fn bounds(actual: ScalarInterval, expected: &R) {
    let (lo, hi) = (actual.lower(), actual.upper());
    assert!(lo.is_finite() && hi.is_finite());
    if lo == hi {
        assert_eq!(&rat(lo), expected);
    } else {
        let next = if lo == 0. {
            f64::from_bits(1)
        } else if lo > 0. {
            f64::from_bits(lo.to_bits() + 1)
        } else {
            f64::from_bits(lo.to_bits() - 1)
        };
        assert_eq!(next, hi);
        assert!(rat(lo) < *expected && *expected < rat(hi));
    }
}

pub fn check(actual: &ExactBezierCurve3, expected: &Arc, probe: &R) {
    assert_eq!(actual.degree(), expected.degree);
    assert_eq!(actual.domain(), &expected.domain);
    assert!(actual.homogeneous_poles().iter().all(|p| p[3] > integer(0)));
    let mut coefficients = expected.coefficients.clone();
    for c in &mut coefficients {
        c.resize(actual.degree() + 1, integer(0));
    }
    assert_eq!(
        self::coefficients(actual),
        coefficients,
        "complete homogeneous polynomial identity"
    );
    let u = &expected.domain[0] + (&expected.domain[1] - &expected.domain[0]) * probe;
    let jet = expected.jet(&u);
    let exact = actual.exact_evaluate(&u, DerivativeOrder::Second).unwrap();
    assert_eq!(exact.position().coordinates(), &jet[0]);
    assert_eq!(exact.derivative(1), Some(&jet[1]));
    assert_eq!(exact.derivative(2), Some(&jet[2]));
    let max = rat(f64::MAX);
    let representable = jet.iter().flatten().all(|x| x >= &-&max && x <= &max);
    match exact.enclosed() {
        Ok(enclosed) => {
            assert!(representable);
            for (b, v) in [
                enclosed.position_bounds(),
                enclosed.derivative_bounds(1).unwrap(),
                enclosed.derivative_bounds(2).unwrap(),
            ]
            .iter()
            .zip(&jet)
            {
                for (b, v) in b.iter().zip(v) {
                    bounds(*b, v);
                }
            }
        }
        Err(Error::Unrepresentable(_)) => assert!(!representable),
        other => panic!("unexpected conversion {other:?}"),
    }
}
