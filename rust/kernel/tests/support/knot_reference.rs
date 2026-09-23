//! Independent Cox polynomial equations over the entire raw support. The
//! coefficient system also proves whether an exact removal exists. No local
//! insertion, inverse insertion, de Boor or production extraction is used.
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::{ExactBSplineCurve3, ExactKnotVector};
use std::collections::{BTreeMap, BTreeSet};

pub fn integer(n: usize) -> R {
    R::from_integer(n.into())
}

fn flat(basis: &ExactKnotVector) -> Vec<R> {
    let p = basis.degree();
    let base: Vec<_> = basis
        .knots()
        .iter()
        .zip(basis.multiplicities())
        .flat_map(|(k, &m)| std::iter::repeat_n(k.clone(), m))
        .collect();
    if !basis.is_periodic() {
        return base;
    }
    let n = basis.pole_count() as isize;
    let pad = (p + 1 - basis.multiplicities()[0]) as isize;
    let period = &basis.domain()[1] - &basis.domain()[0];
    (-pad..n + basis.multiplicities()[0] as isize + pad)
        .map(|i| {
            &base[i.rem_euclid(n) as usize] + R::from_integer(i.div_euclid(n).into()) * &period
        })
        .collect()
}

fn gcd(mut a: BigInt, mut b: BigInt) -> BigInt {
    while b != BigInt::from(0) {
        (a, b) = (b.clone(), a % b);
    }
    if a < BigInt::from(0) {
        -a
    } else {
        a
    }
}
fn lcm(a: &BigInt, b: &BigInt) -> BigInt {
    a / gcd(a.clone(), b.clone()) * b
}

// Keep each Cox polynomial over one denominator. This is still the same
// power-basis recurrence; reductions no longer occur at every scalar product.
struct Polynomial {
    coefficients: Vec<BigInt>,
    denominator: BigInt,
}
impl Polynomial {
    fn zero() -> Self {
        Self {
            coefficients: Vec::new(),
            denominator: 1.into(),
        }
    }
    fn add_linear(&mut self, row: &Self, a: R, b: R) {
        let linear_den = lcm(a.denom(), b.denom());
        let an = a.numer() * (&linear_den / a.denom());
        let bn = b.numer() * (&linear_den / b.denom());
        let term_den = &row.denominator * linear_den;
        let denominator = lcm(&self.denominator, &term_den);
        let old = &denominator / &self.denominator;
        let scale = &denominator / term_den;
        for x in &mut self.coefficients {
            *x *= &old;
        }
        self.coefficients.resize(
            self.coefficients.len().max(row.coefficients.len() + 1),
            0.into(),
        );
        for (i, x) in row.coefficients.iter().enumerate() {
            self.coefficients[i] += &scale * &an * x;
            self.coefficients[i + 1] += &scale * &bn * x;
        }
        self.denominator = denominator;
    }
}
struct Basis {
    rows: Vec<Vec<BigInt>>,
    denominator: BigInt,
}
fn basis_polynomials(basis: &ExactKnotVector, lo: &R, hi: &R) -> Basis {
    let knots = flat(basis);
    let (lo, hi) = if basis.is_periodic() {
        let period = &basis.domain()[1] - &basis.domain()[0];
        let shift = ((lo - &basis.domain()[0]) / &period).floor() * period;
        (lo - &shift, hi - shift)
    } else {
        (lo.clone(), hi.clone())
    };
    let midpoint = (&lo + &hi) / integer(2);
    let mut rows: Vec<_> = knots
        .windows(2)
        .map(|k| {
            if k[0] <= midpoint && midpoint < k[1] {
                Polynomial {
                    coefficients: vec![1.into()],
                    denominator: 1.into(),
                }
            } else {
                Polynomial::zero()
            }
        })
        .collect();
    for d in 1..=basis.degree() {
        rows = (0..rows.len() - 1)
            .map(|i| {
                let mut row = Polynomial::zero();
                let left = &knots[i + d] - &knots[i];
                let right = &knots[i + d + 1] - &knots[i + 1];
                if left != integer(0) && !rows[i].coefficients.is_empty() {
                    row.add_linear(&rows[i], (&lo - &knots[i]) / &left, (&hi - &lo) / left);
                }
                if right != integer(0) && !rows[i + 1].coefficients.is_empty() {
                    row.add_linear(
                        &rows[i + 1],
                        (&knots[i + d + 1] - &lo) / &right,
                        (&lo - &hi) / right,
                    );
                }
                while row.coefficients.last() == Some(&BigInt::from(0)) {
                    row.coefficients.pop();
                }
                row
            })
            .collect();
    }
    let denominator = rows
        .iter()
        .fold(BigInt::from(1), |d, row| lcm(&d, &row.denominator));
    let mut result = vec![vec![BigInt::from(0); basis.degree() + 1]; basis.pole_count()];
    for (i, row) in rows.iter().enumerate() {
        let factor = &denominator / &row.denominator;
        for (k, x) in row.coefficients.iter().enumerate() {
            result[i % basis.pole_count()][k] += x * &factor;
        }
    }
    Basis {
        rows: result,
        denominator,
    }
}

struct Homogeneous {
    coefficients: Vec<[BigInt; 4]>,
    denominator: BigInt,
}
fn homogeneous(curve: &ExactBSplineCurve3, lo: &R, hi: &R) -> Homogeneous {
    let rows = basis_polynomials(curve.knot_vector(), lo, hi);
    let denominator = curve
        .homogeneous_poles()
        .iter()
        .flatten()
        .fold(BigInt::from(1), |d, x| lcm(&d, x.denom()));
    let poles: Vec<[BigInt; 4]> = curve
        .homogeneous_poles()
        .iter()
        .map(|p| std::array::from_fn(|c| p[c].numer() * (&denominator / p[c].denom())))
        .collect();
    Homogeneous {
        coefficients: (0..=curve.degree())
            .map(|k| {
                std::array::from_fn(|c| {
                    rows.rows
                        .iter()
                        .zip(&poles)
                        .map(|(row, p)| &row[k] * &p[c])
                        .sum()
                })
            })
            .collect(),
        denominator: denominator * rows.denominator,
    }
}
pub fn polynomial(curve: &ExactBSplineCurve3, lo: &R, hi: &R) -> [Vec<R>; 4] {
    let h = homogeneous(curve, lo, hi);
    std::array::from_fn(|c| {
        h.coefficients
            .iter()
            .map(|p| R::new(p[c].clone(), h.denominator.clone()))
            .collect()
    })
}

pub fn partition(original: &ExactKnotVector, candidate: &ExactKnotVector) -> Vec<[R; 2]> {
    let mut cuts = BTreeSet::new();
    if original.is_periodic() {
        let a = &original.domain()[0];
        let period = &original.domain()[1] - a;
        for basis in [original, candidate] {
            for k in basis.knots() {
                cuts.insert(k - ((k - a) / &period).floor() * &period);
            }
        }
        cuts.insert(original.domain()[1].clone());
    } else {
        cuts.extend(original.knots().iter().cloned());
        cuts.extend(candidate.knots().iter().cloned());
    }
    let cuts: Vec<_> = cuts.into_iter().collect();
    cuts.windows(2)
        .map(|k| [k[0].clone(), k[1].clone()])
        .collect()
}

pub fn equal(original: &ExactBSplineCurve3, candidate: &ExactBSplineCurve3) -> bool {
    if original.degree() != candidate.degree()
        || original.is_periodic() != candidate.is_periodic()
        || (original.is_periodic()
            && original.domain()[1].clone() - &original.domain()[0]
                != candidate.domain()[1].clone() - &candidate.domain()[0])
    {
        return false;
    }
    partition(original.knot_vector(), candidate.knot_vector())
        .iter()
        .all(|[a, b]| {
            let (a, b) = (homogeneous(original, a, b), homogeneous(candidate, a, b));
            a.coefficients
                .iter()
                .flatten()
                .zip(b.coefficients.iter().flatten())
                .all(|(x, y)| x * &b.denominator == y * &a.denominator)
        })
}

type Equation = (usize, BTreeMap<usize, BigInt>, [BigInt; 4]);
struct Solution {
    controls: Vec<[R; 4]>,
    integers: Vec<[BigInt; 4]>,
    denominator: BigInt,
}
impl Solution {
    fn new(pivots: &[Equation], count: usize) -> Self {
        let mut controls = vec![std::array::from_fn(|_| integer(0)); count];
        for (j, row, rhs) in pivots.iter().rev() {
            controls[*j] = std::array::from_fn(|c| {
                (R::from_integer(rhs[c].clone())
                    - row
                        .iter()
                        .filter(|(i, _)| **i != *j)
                        .map(|(i, x)| R::from_integer(x.clone()) * &controls[*i][c])
                        .sum::<R>())
                    / R::from_integer(row[j].clone())
            });
        }
        let denominator = controls
            .iter()
            .flatten()
            .fold(BigInt::from(1), |d, x| lcm(&d, x.denom()));
        let integers = controls
            .iter()
            .map(|p| std::array::from_fn(|c| p[c].numer() * (&denominator / p[c].denom())))
            .collect();
        Self {
            controls,
            integers,
            denominator,
        }
    }
    fn satisfies(&self, basis: &Basis, rhs: &Homogeneous, k: usize) -> bool {
        (0..4).all(|c| {
            let sum: BigInt = basis
                .rows
                .iter()
                .zip(&self.integers)
                .map(|(row, p)| &row[k] * &p[c])
                .sum();
            sum * &rhs.denominator
                == &rhs.coefficients[k][c] * &basis.denominator * &self.denominator
        })
    }
}

pub fn recover(original: &ExactBSplineCurve3, basis: &ExactKnotVector) -> Option<Vec<[R; 4]>> {
    let mut pivots: Vec<Equation> = Vec::new();
    let mut solution: Option<Solution> = None;
    for [lo, hi] in partition(original.knot_vector(), basis) {
        let rows = basis_polynomials(basis, &lo, &hi);
        let rhs = homogeneous(original, &lo, &hi);
        for k in 0..=basis.degree() {
            // Once the coefficient equations determine a unique control row,
            // check every remaining equation by exact substitution. Eliminating
            // hundreds of dependent equations again is unnecessary; none of
            // their residual checks is omitted or replaced by point samples.
            if let Some(solved) = &solution {
                if !solved.satisfies(&rows, &rhs, k) {
                    return None;
                }
                continue;
            }
            let mut row: BTreeMap<_, _> = rows
                .rows
                .iter()
                .enumerate()
                .filter(|(_, r)| r[k] != BigInt::from(0))
                .map(|(i, r)| (i, &r[k] * &rhs.denominator))
                .collect();
            let mut answer: [BigInt; 4] =
                std::array::from_fn(|c| &rhs.coefficients[k][c] * &rows.denominator);
            // Normalize each source equation once, then use fraction-free
            // elimination. A previous pivot divides every subsequent minor;
            // verify each division exactly, including under release fuzzing.
            let mut content = BigInt::from(0);
            for x in row.values().chain(&answer) {
                content = gcd(content, x.clone());
                if content == BigInt::from(1) {
                    break;
                }
            }
            if content > BigInt::from(1) {
                for x in row.values_mut().chain(&mut answer) {
                    *x /= &content;
                }
            }
            let mut previous = BigInt::from(1);
            for (j, pivot, values) in &pivots {
                if row.is_empty() {
                    break;
                }
                let a = &pivot[j];
                let b = row.get(j).cloned().unwrap_or_else(|| 0.into());
                for x in row.values_mut() {
                    *x *= a;
                }
                if b != BigInt::from(0) {
                    for (&i, x) in pivot {
                        let y = row.get(&i).cloned().unwrap_or_else(|| 0.into()) - &b * x;
                        if y == BigInt::from(0) {
                            row.remove(&i);
                        } else {
                            row.insert(i, y);
                        }
                    }
                }
                answer = std::array::from_fn(|c| a * &answer[c] - &b * &values[c]);
                if previous != BigInt::from(1) {
                    for x in row.values_mut().chain(&mut answer) {
                        let q = &*x / &previous;
                        assert_eq!(&q * &previous, *x, "exact fraction-free division");
                        *x = q;
                    }
                }
                previous = a.clone();
            }
            if let Some((&j, _)) = row.first_key_value() {
                pivots.push((j, row, answer));
                if pivots.len() == basis.pole_count() {
                    let solved = Solution::new(&pivots, basis.pole_count());
                    assert!(solved.satisfies(&rows, &rhs, k));
                    solution = Some(solved);
                }
            } else if answer.iter().any(|x| x != &BigInt::from(0)) {
                return None;
            }
        }
    }
    assert_eq!(pivots.len(), basis.pole_count(), "independent basis rank");
    Some(solution.expect("independent basis has full rank").controls)
}

/// For valid removal requests which leave more poles than degree.
pub fn removed(curve: &ExactBSplineCurve3, u: &R, target: usize) -> Option<ExactBSplineCurve3> {
    let seam = curve.is_periodic() && (u == &curve.domain()[0] || u == &curve.domain()[1]);
    let u = if seam { &curve.domain()[0] } else { u };
    let mut map: BTreeMap<_, _> = curve
        .knots()
        .iter()
        .cloned()
        .zip(curve.multiplicities().iter().copied())
        .collect();
    if target >= map[u] {
        return Some(curve.clone());
    }
    if target == 0 {
        map.remove(u);
        if seam {
            map.remove(&curve.domain()[1]);
            let (a, m) = map.first_key_value().unwrap();
            map.insert(a + (&curve.domain()[1] - &curve.domain()[0]), *m);
        }
    } else {
        map.insert(u.clone(), target);
        if seam {
            map.insert(curve.domain()[1].clone(), target);
        }
    }
    let (knots, mults) = map.into_iter().unzip();
    let basis = if curve.is_periodic() {
        ExactKnotVector::new_periodic(curve.degree(), knots, mults)
    } else {
        ExactKnotVector::new(curve.degree(), knots, mults)
    }
    .unwrap();
    let controls = recover(curve, &basis)?;
    if controls.iter().any(|c| c[3] <= integer(0)) {
        return None;
    }
    Some(ExactBSplineCurve3::from_homogeneous(basis, controls).unwrap())
}
