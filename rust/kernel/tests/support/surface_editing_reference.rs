//! Independent Cox tensor/power identities. No production extraction/editing
//! or evaluation routine supplies any expected coefficient or derivative.
pub(crate) use crate::bezier_reference as curves;
use curves::{binomial, integer, integers, powers, substitute};
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
            let basis = curves::cox_basis(&flat, span, p);
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

/// Independent power coefficients share one positive denominator throughout
/// tensor transforms. Fractions are reduced only when a scalar is returned;
/// coefficient equality cross-multiplies complete integer grids exactly.
#[derive(Clone)]
struct Grid {
    values: Vec<[BigInt; 4]>,
    denominator: BigInt,
}
impl Grid {
    fn from_rationals(values: &[[R; 4]]) -> Self {
        let (values, denominator) = integers(&values.iter().flatten().cloned().collect::<Vec<_>>());
        Self {
            values: values
                .chunks(4)
                .map(|p| std::array::from_fn(|c| p[c].clone()))
                .collect(),
            denominator,
        }
    }
    fn rational(&self, i: usize, c: usize) -> R {
        R::new(self.values[i][c].clone(), self.denominator.clone())
    }
    fn map(&self, degrees: [usize; 2], axis: usize, matrix: &Matrix) -> (Self, [usize; 2]) {
        let mut new_degrees = degrees;
        new_degrees[axis] = matrix.rows.len() - 1;
        let values = (0..=new_degrees[0])
            .flat_map(|u| {
                (0..=new_degrees[1]).map(move |v| {
                    let (row, fixed) = if axis == 0 { (u, v) } else { (v, u) };
                    std::array::from_fn(|c| {
                        matrix.rows[row]
                            .iter()
                            .enumerate()
                            .map(|(i, m)| {
                                let index = if axis == 0 {
                                    i * (degrees[1] + 1) + fixed
                                } else {
                                    fixed * (degrees[1] + 1) + i
                                };
                                m * &self.values[index][c]
                            })
                            .sum()
                    })
                })
            })
            .collect();
        (
            Self {
                values,
                denominator: &self.denominator * &matrix.denominator,
            },
            new_degrees,
        )
    }
    fn assert_equal(&self, other: &Self) {
        assert_eq!(self.values.len(), other.values.len());
        let scale = R::new(other.denominator.clone(), self.denominator.clone());
        for (a, b) in self
            .values
            .iter()
            .flatten()
            .zip(other.values.iter().flatten())
        {
            assert_eq!(
                a * scale.numer(),
                b * scale.denom(),
                "complete homogeneous tensor polynomial identity"
            );
        }
    }
}

struct Matrix {
    rows: Vec<Vec<BigInt>>,
    denominator: BigInt,
}
impl Matrix {
    fn new(rows: Vec<(Vec<BigInt>, BigInt)>) -> Self {
        let (_, denominator) = integers(
            &rows
                .iter()
                .map(|(_, d)| R::new(1.into(), d.clone()))
                .collect::<Vec<_>>(),
        );
        Self {
            rows: rows
                .into_iter()
                .map(|(r, d)| {
                    let factor = &denominator / d;
                    r.into_iter().map(|x| x * &factor).collect()
                })
                .collect(),
            denominator,
        }
    }
}

/// Binomial affine substitution on power coefficients. Build each exact
/// matrix once per axis/operation, independently of the control data.
fn affine_matrix(degree: usize, a: &R, b: &R) -> Matrix {
    let n = degree + 1;
    let (an, ad, bn, bd) = (
        powers(a.numer(), n),
        powers(a.denom(), n),
        powers(b.numer(), n),
        powers(b.denom(), n),
    );
    Matrix::new(
        (0..n)
            .map(|i| {
                let row = (0..n)
                    .map(|j| {
                        if j < i {
                            BigInt::from(0)
                        } else {
                            binomial(j, i) * &an[j - i] * &ad[degree - j] * &bn[i]
                        }
                    })
                    .collect();
                (row, &ad[degree - i] * &bd[i])
            })
            .collect(),
    )
}

/// Direct monomial derivative sums, sharing powers and input denominators
/// between derivative orders. This does not use a Bernstein recurrence.
fn value_matrix(degree: usize, t: &R, order: usize) -> Matrix {
    let (tn, td) = (powers(t.numer(), degree + 1), powers(t.denom(), degree + 1));
    Matrix::new(
        (0..=order)
            .map(|r| {
                let row = (0..=degree)
                    .map(|i| {
                        if i < r {
                            BigInt::from(0)
                        } else {
                            let factor: usize = (i + 1 - r..=i).product();
                            factor * &tn[i - r] * &td[degree - i]
                        }
                    })
                    .collect();
                (row, td[degree - r.min(degree)].clone())
            })
            .collect(),
    )
}

#[derive(Clone)]
pub struct Patch {
    pub degrees: [usize; 2],
    pub domain: [[R; 2]; 2],
    coefficients: Grid,
}
impl Patch {
    pub fn trim(&self, domain: [[R; 2]; 2]) -> Self {
        let mut out = self.clone();
        for (axis, range) in domain.iter().enumerate() {
            if range == &self.domain[axis] {
                continue;
            }
            let length = &self.domain[axis][1] - &self.domain[axis][0];
            let a = (&range[0] - &self.domain[axis][0]) / &length;
            let b = (&range[1] - &range[0]) / length;
            let matrix = affine_matrix(out.degrees[axis], &a, &b);
            (out.coefficients, out.degrees) = out.coefficients.map(out.degrees, axis, &matrix);
        }
        out.domain = domain;
        out
    }
    pub fn reversed(&self, axis: usize) -> Self {
        let matrix = affine_matrix(self.degrees[axis], &integer(1), &-integer(1));
        let (coefficients, degrees) = self.coefficients.map(self.degrees, axis, &matrix);
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
            coefficients: Grid {
                values: (0..=self.degrees[1])
                    .flat_map(|j| {
                        (0..=self.degrees[0]).map(move |i| {
                            self.coefficients.values[i * (self.degrees[1] + 1) + j].clone()
                        })
                    })
                    .collect(),
                denominator: self.coefficients.denominator.clone(),
            },
        }
    }
    pub fn elevated(&self, target: [usize; 2]) -> Self {
        let values = (0..=target[0])
            .flat_map(|i| {
                (0..=target[1]).map(move |j| {
                    if i <= self.degrees[0] && j <= self.degrees[1] {
                        self.coefficients.values[i * (self.degrees[1] + 1) + j].clone()
                    } else {
                        std::array::from_fn(|_| BigInt::from(0))
                    }
                })
            })
            .collect();
        Self {
            degrees: target,
            domain: self.domain.clone(),
            coefficients: Grid {
                values,
                denominator: self.coefficients.denominator.clone(),
            },
        }
    }
    pub fn iso(&self, axis: usize, t: &R) -> curves::Arc {
        let (coefficients, _) =
            self.coefficients
                .map(self.degrees, axis, &value_matrix(self.degrees[axis], t, 0));
        curves::Arc {
            degree: self.degrees[1 - axis],
            domain: self.domain[1 - axis].clone(),
            coefficients: std::array::from_fn(|c| {
                (0..coefficients.values.len())
                    .map(|i| coefficients.rational(i, c))
                    .collect()
            }),
        }
    }
    pub fn jet(&self, parameters: [&R; 2]) -> [[R; 3]; 6] {
        let length: [R; 2] = std::array::from_fn(|i| &self.domain[i][1] - &self.domain[i][0]);
        let t: [R; 2] = std::array::from_fn(|i| (parameters[i] - &self.domain[i][0]) / &length[i]);
        let vm = value_matrix(self.degrees[1], &t[1], 2);
        let um = value_matrix(self.degrees[0], &t[0], 2);
        let (v_jets, d) = self.coefficients.map(self.degrees, 1, &vm);
        let (jets, _) = v_jets.map(d, 0, &um);
        let h: [[R; 4]; 6] = std::array::from_fn(|r| {
            let (u, v) = PARTIALS[r];
            let units = length[0].pow(u as i32) * length[1].pow(v as i32);
            std::array::from_fn(|c| jets.rational(u * 3 + v, c) / &units)
        });
        // The common scale cancels from every closed quotient formula. Keep
        // integer numerators until the final scalar so extreme weights do not
        // trigger a rational GCD at each intermediate multiply and subtract.
        let (values, _) = integers(&h.iter().flatten().cloned().collect::<Vec<_>>());
        let h: [[BigInt; 4]; 6] =
            std::array::from_fn(|i| std::array::from_fn(|c| values[4 * i + c].clone()));
        let w = &h[0][3];
        let w2 = w * w;
        let w3 = &w2 * w;
        std::array::from_fn(|r| {
            std::array::from_fn(|c| match r {
                0 => R::new(h[0][c].clone(), w.clone()),
                1 | 2 => R::new(&h[r][c] * w - &h[0][c] * &h[r][3], w2.clone()),
                3 | 4 => {
                    let i = r - 2;
                    R::new(
                        &h[r][c] * &w2 - &h[0][c] * &h[r][3] * w - 2 * &h[i][c] * w * &h[i][3]
                            + 2 * &h[0][c] * &h[i][3] * &h[i][3],
                        w3.clone(),
                    )
                }
                _ => R::new(
                    &h[5][c] * &w2
                        - &h[1][c] * w * &h[2][3]
                        - &h[2][c] * w * &h[1][3]
                        - &h[0][c] * w * &h[5][3]
                        + 2 * &h[0][c] * &h[1][3] * &h[2][3],
                    w3.clone(),
                ),
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
        let um = Matrix::new(
            (0..=degrees[0])
                .map(|k| {
                    integers(
                        &u.coefficients
                            .iter()
                            .map(|c| c[k].clone())
                            .collect::<Vec<_>>(),
                    )
                })
                .collect(),
        );
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
            let vm = Matrix::new(
                (0..=degrees[1])
                    .map(|k| {
                        integers(
                            &v.coefficients
                                .iter()
                                .map(|c| c[k].clone())
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect(),
            );
            let (controls, d) = Grid::from_rationals(&controls).map(degrees, 0, &um);
            let (coefficients, degrees) = controls.map(d, 1, &vm);
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

#[allow(dead_code)] // Cargo fixtures use this entry point; fuzz profiling uses the one below.
pub fn check(actual: &ExactBezierSurface3, expected: &Patch, probe: [&R; 2]) {
    check_profiled(actual, expected, probe, |_| {});
}

pub fn check_profiled(
    actual: &ExactBezierSurface3,
    expected: &Patch,
    probe: [&R; 2],
    mut mark: impl FnMut(&'static str),
) {
    assert_eq!(actual.degrees(), expected.degrees);
    assert_eq!(actual.domain(), &expected.domain);
    assert!(actual.homogeneous_poles().iter().all(|p| p[3] > integer(0)));
    let mut coefficients = Grid::from_rationals(actual.homogeneous_poles());
    let degrees = actual.degrees();
    for axis in 0..2 {
        let degree = degrees[axis];
        let matrix = Matrix::new(
            (0..=degree)
                .map(|k| {
                    let row = (0..=k)
                        .map(|i| {
                            let x = BigInt::from(binomial(degree, i) * binomial(degree - i, k - i));
                            if (k - i) % 2 == 0 {
                                x
                            } else {
                                -x
                            }
                        })
                        .collect();
                    (row, BigInt::from(1))
                })
                .collect(),
        );
        (coefficients, _) = coefficients.map(degrees, axis, &matrix);
    }
    coefficients.assert_equal(&expected.coefficients);
    mark("oracle tensor controls");
    let at: [R; 2] = std::array::from_fn(|i| {
        &expected.domain[i][0] + (&expected.domain[i][1] - &expected.domain[i][0]) * probe[i]
    });
    let jet = expected.jet([&at[0], &at[1]]);
    mark("oracle tensor jet");
    let exact = actual
        .exact_evaluate(&at[0], &at[1], DerivativeOrder::Second)
        .unwrap();
    assert_eq!(exact.position().coordinates(), &jet[0]);
    for (i, &(u, v)) in PARTIALS.iter().enumerate().skip(1) {
        assert_eq!(exact.derivative(u, v), Some(&jet[i]));
    }
    mark("kernel tensor jet");
    let max = rat(f64::MAX);
    let representable = jet.iter().flatten().all(|x| x >= &-&max && x <= &max);
    let enclosed = exact.enclosed();
    mark("kernel tensor enclosure");
    match enclosed {
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
    mark("oracle tensor bounds");
}

/// Complete tensor Cox power polynomial on a cell without an interior knot.
/// Exact rational axes/control data remain exact throughout this independent
/// coefficient transform. No kernel de Boor or extraction routine is used.
#[allow(dead_code)]
pub fn exact_patch(
    surface: &rusty_occt::ExactBSplineSurface3,
    domain: [[R; 2]; 2],
    basis: impl Fn(&rusty_occt::ExactKnotVector, &R, &R) -> Vec<Vec<R>>,
) -> Patch {
    let mut coefficients = Grid::from_rationals(surface.homogeneous_poles());
    let mut degrees = surface.pole_counts().map(|n| n - 1);
    for (i, axis) in [surface.u_knots(), surface.v_knots()]
        .into_iter()
        .enumerate()
    {
        let polynomials = basis(axis, &domain[i][0], &domain[i][1]);
        let matrix = Matrix::new(
            (0..=axis.degree())
                .map(|k| integers(&polynomials.iter().map(|c| c[k].clone()).collect::<Vec<_>>()))
                .collect(),
        );
        (coefficients, degrees) = coefficients.map(degrees, i, &matrix);
    }
    Patch {
        degrees,
        domain,
        coefficients,
    }
}
