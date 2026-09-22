//! Independent Cox tensor/power identities. No production extraction/editing
//! or evaluation routine supplies any expected coefficient or derivative.
pub(crate) use crate::bezier_reference as curves;
use curves::{add, binomial, integer, integers, linear, substitute, value};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::curve::DerivativeOrder;
use rusty_occt::{BSplineSurface3, Error, ExactBezierSurface3, KnotVector};
const PARTIALS: [(usize, usize); 6] = [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1)];
fn rat(x: f64) -> R {
    R::from_float(x).unwrap()
}

struct Span {
    domain: [R; 2],
    indices: Vec<usize>,
    coefficients: Vec<Vec<R>>,
}
fn spans(axis: &KnotVector, first: f64, last: f64) -> Vec<Span> {
    let p = axis.degree();
    let base: Vec<_> = axis
        .knots()
        .iter()
        .zip(axis.multiplicities())
        .flat_map(|(&k, &m)| std::iter::repeat_n(rat(k), m))
        .collect();
    let (start, end) = axis.domain();
    let (start, end, first, last) = (rat(start), rat(end), rat(first), rat(last));
    let period = &end - &start;
    let flat = if axis.is_periodic() {
        let n = (base.len() - axis.multiplicities()[0]) as isize;
        let pad = (p + 1 - axis.multiplicities()[0]) as isize;
        (-pad..n + axis.multiplicities()[0] as isize + pad)
            .map(|i| {
                &base[i.rem_euclid(n) as usize]
                    + R::from_integer(BigInt::from(i.div_euclid(n))) * &period
            })
            .collect()
    } else {
        base
    };
    let mut offset = if axis.is_periodic() {
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
                        let l = &flat[i + degree] - &flat[i];
                        let r = &flat[i + degree + 1] - &flat[i + 1];
                        let x = if l == integer(0) || basis[i].is_empty() {
                            Vec::new()
                        } else {
                            linear(&basis[i], &((a - &flat[i]) / &l), &((b - a) / &l))
                        };
                        let y = if r == integer(0) || basis[i + 1].is_empty() {
                            Vec::new()
                        } else {
                            linear(
                                &basis[i + 1],
                                &((&flat[i + degree + 1] - a) / &r),
                                &((a - b) / &r),
                            )
                        };
                        add(&x, &y)
                    })
                    .collect();
            }
            let coefficients = (span - p..=span)
                .map(|i| {
                    let mut c = substitute(
                        &basis[i],
                        &((&lo - a - &offset) / (b - a)),
                        &((&hi - &lo) / (b - a)),
                    );
                    c.resize(p + 1, integer(0));
                    c
                })
                .collect();
            result.push(Span {
                domain: [lo, hi],
                indices: (span - p..=span).map(|i| i % axis.pole_count()).collect(),
                coefficients,
            });
        }
        if !axis.is_periodic() || &end + &offset >= last {
            break;
        }
        offset += &period;
    }
    result
}

fn map_axis(
    grid: &[[R; 4]],
    degrees: [usize; 2],
    axis: usize,
    mut transform: impl FnMut(&[R]) -> Vec<R>,
) -> (Vec<[R; 4]>, [usize; 2]) {
    let mut out = Vec::new();
    let mut new_degrees = degrees;
    for fixed in 0..=degrees[1 - axis] {
        let indices: Vec<_> = (0..=degrees[axis])
            .map(|i| {
                if axis == 0 {
                    i * (degrees[1] + 1) + fixed
                } else {
                    fixed * (degrees[1] + 1) + i
                }
            })
            .collect();
        let rows: [Vec<R>; 4] = std::array::from_fn(|c| {
            transform(
                &indices
                    .iter()
                    .map(|&i| grid[i][c].clone())
                    .collect::<Vec<_>>(),
            )
        });
        new_degrees[axis] = rows[0].len() - 1;
        if out.is_empty() {
            out.resize_with((new_degrees[0] + 1) * (new_degrees[1] + 1), || {
                std::array::from_fn(|_| integer(0))
            });
        }
        for (i, _) in rows[0].iter().enumerate() {
            let index = if axis == 0 {
                i * (new_degrees[1] + 1) + fixed
            } else {
                fixed * (new_degrees[1] + 1) + i
            };
            out[index] = std::array::from_fn(|c| rows[c][i].clone());
        }
    }
    (out, new_degrees)
}

/// Matrix multiplication after clearing denominators. This is the Cox
/// coefficient sum, not a de Boor/de Casteljau control-point recurrence.
fn multiply(matrix: &[Vec<R>], poles: &[R]) -> Vec<R> {
    let (p, den) = integers(poles);
    matrix
        .iter()
        .map(|row| {
            let (m, d) = integers(row);
            R::new(p.iter().zip(&m).map(|(a, b)| a * b).sum(), &den * d)
        })
        .collect()
}

#[derive(Clone)]
pub struct Patch {
    pub degrees: [usize; 2],
    pub domain: [[R; 2]; 2],
    pub coefficients: Vec<[R; 4]>,
}
impl Patch {
    pub fn trim(&self, domain: [[R; 2]; 2]) -> Self {
        let mut out = self.clone();
        for (axis, range) in domain.iter().enumerate() {
            let length = &self.domain[axis][1] - &self.domain[axis][0];
            let a = (&range[0] - &self.domain[axis][0]) / &length;
            let b = (&range[1] - &range[0]) / length;
            (out.coefficients, out.degrees) = map_axis(&out.coefficients, out.degrees, axis, |p| {
                substitute(p, &a, &b)
            });
        }
        out.domain = domain;
        out
    }
    pub fn reversed(&self, axis: usize) -> Self {
        let (coefficients, degrees) = map_axis(&self.coefficients, self.degrees, axis, |p| {
            substitute(p, &integer(1), &-integer(1))
        });
        Self {
            coefficients,
            degrees,
            domain: self.domain.clone(),
        }
    }
    pub fn exchanged(&self) -> Self {
        Self {
            degrees: [self.degrees[1], self.degrees[0]],
            domain: [self.domain[1].clone(), self.domain[0].clone()],
            coefficients: (0..=self.degrees[1])
                .flat_map(|j| {
                    (0..=self.degrees[0])
                        .map(move |i| self.coefficients[i * (self.degrees[1] + 1) + j].clone())
                })
                .collect(),
        }
    }
    pub fn elevated(&self, target: [usize; 2]) -> Self {
        let mut out = self.clone();
        for (axis, &d) in target.iter().enumerate() {
            (out.coefficients, out.degrees) = map_axis(&out.coefficients, out.degrees, axis, |p| {
                let mut p = p.to_vec();
                p.resize(d + 1, integer(0));
                p
            });
        }
        out
    }
    pub fn iso(&self, axis: usize, t: &R) -> curves::Arc {
        let (coefficients, _) = map_axis(&self.coefficients, self.degrees, axis, |p| {
            vec![value(p, t, 0)]
        });
        curves::Arc {
            degree: self.degrees[1 - axis],
            domain: self.domain[1 - axis].clone(),
            coefficients: std::array::from_fn(|c| {
                coefficients.iter().map(|p| p[c].clone()).collect()
            }),
        }
    }
    pub fn jet(&self, parameters: [&R; 2]) -> [[R; 3]; 6] {
        let length: [R; 2] = std::array::from_fn(|i| &self.domain[i][1] - &self.domain[i][0]);
        let t: [R; 2] = std::array::from_fn(|i| (parameters[i] - &self.domain[i][0]) / &length[i]);
        let h: [[R; 4]; 6] = std::array::from_fn(|r| {
            let (uo, vo) = PARTIALS[r];
            std::array::from_fn(|c| {
                let rows: Vec<_> = self
                    .coefficients
                    .chunks(self.degrees[1] + 1)
                    .map(|row| {
                        value(
                            &row.iter().map(|p| p[c].clone()).collect::<Vec<_>>(),
                            &t[1],
                            vo,
                        )
                    })
                    .collect();
                value(&rows, &t[0], uo) / (length[0].pow(uo as i32) * length[1].pow(vo as i32))
            })
        });
        let w = &h[0][3];
        std::array::from_fn(|r| {
            std::array::from_fn(|c| match r {
                0 => &h[0][c] / w,
                1 | 2 => (&h[r][c] * w - &h[0][c] * &h[r][3]) / (w * w),
                3 | 4 => {
                    let i = r - 2;
                    (&h[r][c] * w * w
                        - &h[0][c] * &h[r][3] * w
                        - integer(2) * &h[i][c] * w * &h[i][3]
                        + integer(2) * &h[0][c] * &h[i][3] * &h[i][3])
                        / (w * w * w)
                }
                _ => {
                    (&h[5][c] * w * w
                        - &h[1][c] * w * &h[2][3]
                        - &h[2][c] * w * &h[1][3]
                        - &h[0][c] * w * &h[5][3]
                        + integer(2) * &h[0][c] * &h[1][3] * &h[2][3])
                        / (w * w * w)
                }
            })
        })
    }
}

pub fn extract(surface: &BSplineSurface3, rectangle: [f64; 4]) -> Vec<Patch> {
    let [ua, ub, va, vb] = rectangle;
    let us = spans(surface.u_knots(), ua, ub);
    let vs = spans(surface.v_knots(), va, vb);
    let degrees = [surface.u_knots().degree(), surface.v_knots().degree()];
    let mut result = Vec::new();
    for u in &us {
        let um: Vec<Vec<R>> = (0..=degrees[0])
            .map(|k| u.coefficients.iter().map(|c| c[k].clone()).collect())
            .collect();
        for v in &vs {
            let mut controls = Vec::new();
            for &i in &u.indices {
                for &j in &v.indices {
                    let index = i * surface.v_knots().pole_count() + j;
                    let p = surface.poles()[index].to_array();
                    let w = rat(surface.weights()[index]);
                    controls.push(std::array::from_fn(|c| {
                        if c == 3 {
                            w.clone()
                        } else {
                            rat(p[c]) * &w
                        }
                    }));
                }
            }
            let vm: Vec<Vec<R>> = (0..=degrees[1])
                .map(|k| v.coefficients.iter().map(|c| c[k].clone()).collect())
                .collect();
            let (controls, d) = map_axis(&controls, degrees, 0, |p| multiply(&um, p));
            let (coefficients, degrees) = map_axis(&controls, d, 1, |p| multiply(&vm, p));
            result.push(Patch {
                degrees,
                coefficients,
                domain: [u.domain.clone(), v.domain.clone()],
            });
        }
    }
    result
}

pub enum Item {
    Patch(Patch),
    Curve(curves::Arc),
}
pub fn apply(mut p: Patch, op: u8, elevation: [usize; 2]) -> Vec<Item> {
    if op == 1 || op == 10 {
        let domain = std::array::from_fn(|i| {
            let [a, b] = &p.domain[i];
            [
                a + (b - a) / integer(4),
                a + (b - a) * R::new(3.into(), 4.into()),
            ]
        });
        p = p.trim(domain);
    }
    if op == 7 || op == 10 {
        p = p.elevated(elevation);
    }
    if op == 4 || op == 10 {
        p = p.reversed(0);
    }
    if op == 5 {
        p = p.reversed(1);
    }
    if op == 6 || op == 10 {
        p = p.exchanged();
    }
    if op == 8 || op == 9 {
        return vec![Item::Curve(p.iso(
            if op == 8 { 0 } else { 1 },
            &R::new(if op == 8 { 3.into() } else { 5.into() }, 8.into()),
        ))];
    }
    let mut ranges: [Vec<[R; 2]>; 2] = std::array::from_fn(|i| vec![p.domain[i].clone()]);
    for (axis, range) in ranges.iter_mut().enumerate() {
        if op == 10 || op == axis as u8 + 2 {
            let [a, b] = &p.domain[axis];
            let cut = a + (b - a) * R::new((3 + 2 * axis).into(), 8.into());
            *range = vec![[a.clone(), cut.clone()], [cut, b.clone()]];
        }
    }
    ranges[0]
        .iter()
        .flat_map(|u| {
            ranges[1]
                .iter()
                .map(|v| Item::Patch(p.trim([u.clone(), v.clone()])))
        })
        .collect()
}

pub fn check(actual: &ExactBezierSurface3, expected: &Patch, probe: [&R; 2]) {
    assert_eq!(actual.degrees(), expected.degrees);
    assert_eq!(actual.domain(), &expected.domain);
    assert!(actual.homogeneous_poles().iter().all(|p| p[3] > integer(0)));
    let mut coefficients = actual.homogeneous_poles().to_vec();
    let degrees = actual.degrees();
    for axis in 0..2 {
        (coefficients, _) = map_axis(&coefficients, degrees, axis, |p| {
            let (p, den) = integers(p);
            let degree = p.len() - 1;
            (0..=degree)
                .map(|k| {
                    let sum: BigInt = (0..=k)
                        .map(|i| {
                            let x = &p[i] * (binomial(degree, i) * binomial(degree - i, k - i));
                            if (k - i) % 2 == 0 {
                                x
                            } else {
                                -x
                            }
                        })
                        .sum();
                    R::new(sum, den.clone())
                })
                .collect()
        });
    }
    assert_eq!(
        coefficients, expected.coefficients,
        "complete homogeneous tensor polynomial identity"
    );
    let at: [R; 2] = std::array::from_fn(|i| {
        &expected.domain[i][0] + (&expected.domain[i][1] - &expected.domain[i][0]) * probe[i]
    });
    let jet = expected.jet([&at[0], &at[1]]);
    let exact = actual
        .exact_evaluate(&at[0], &at[1], DerivativeOrder::Second)
        .unwrap();
    assert_eq!(exact.position().coordinates(), &jet[0]);
    for (i, &(u, v)) in PARTIALS.iter().enumerate().skip(1) {
        assert_eq!(exact.derivative(u, v), Some(&jet[i]));
    }
    let max = rat(f64::MAX);
    let representable = jet.iter().flatten().all(|x| x >= &-&max && x <= &max);
    match exact.enclosed() {
        Ok(enclosed) => {
            assert!(representable);
            for (i, &(u, v)) in PARTIALS.iter().enumerate() {
                let bounds = if i == 0 {
                    enclosed.position_bounds()
                } else {
                    enclosed.derivative_bounds(u, v).unwrap()
                };
                for (b, e) in bounds.iter().zip(&jet[i]) {
                    curves::bounds(*b, e);
                }
            }
        }
        Err(Error::Unrepresentable(_)) => assert!(!representable),
        other => panic!("unexpected enclosure {other:?}"),
    }
}
