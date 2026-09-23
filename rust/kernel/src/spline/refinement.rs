//! Reusable exact Boehm maps for tensor refinement. The same local insertion
//! and five-period support/cropping conventions as ExactBSplineCurve3 are
//! applied to sparse unit controls, then shared across all transverse fields.
//! OCCT references: BSplSLib::InsertKnots and BSplCLib::InsertKnots.
use super::{common_denominator, integer, ExactKnotVector};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::collections::BTreeMap;

type Sparse = Vec<(usize, R)>;
fn blend(left: &Sparse, right: &Sparse, alpha: &R) -> Sparse {
    if alpha == &integer(0) {
        return left.clone();
    }
    if alpha == &integer(1) {
        return right.clone();
    }
    let complement = integer(1) - alpha;
    let mut result: BTreeMap<usize, R> = left.iter().map(|(i, x)| (*i, &complement * x)).collect();
    for (i, x) in right {
        *result.entry(*i).or_insert_with(|| integer(0)) += alpha * x;
    }
    result
        .into_iter()
        .filter(|(_, x)| x != &integer(0))
        .collect()
}
struct Work {
    knots: Vec<R>,
    rows: Vec<Sparse>,
}
impl Work {
    fn new(basis: &ExactKnotVector) -> Self {
        let n = basis.pole_count();
        if !basis.is_periodic() {
            return Self {
                knots: basis.flat.clone(),
                rows: (0..n).map(|i| vec![(i, integer(1))]).collect(),
            };
        }
        let extension = basis.degree() + 1 - basis.multiplicities()[0];
        let base: Vec<_> = basis
            .knots()
            .iter()
            .zip(basis.multiplicities())
            .take(basis.knots().len() - 1)
            .flat_map(|(k, &m)| std::iter::repeat_n(k, m))
            .collect();
        let period = &basis.domain()[1] - &basis.domain()[0];
        let start = -2 * (n as isize);
        let count = 5 * n + extension;
        let knots = (start..start + (count + basis.degree() + 1) as isize)
            .map(|i| {
                let j = i - extension as isize;
                base[j.rem_euclid(n as isize) as usize]
                    + R::from_integer(j.div_euclid(n as isize).into()) * &period
            })
            .collect();
        Self {
            knots,
            rows: (0..count).map(|i| vec![(i % n, integer(1))]).collect(),
        }
    }
    fn insert(&mut self, p: usize, u: &R) {
        let k = self.knots.partition_point(|x| x <= u) - 1;
        let mult = k + 1 - self.knots.partition_point(|x| x < u);
        let (first, last) = (k - p + 1, k - mult);
        let changed: Vec<_> = (first..=last)
            .map(|i| {
                let alpha = (u - &self.knots[i]) / (&self.knots[i + p] - &self.knots[i]);
                debug_assert!(alpha >= integer(0) && alpha <= integer(1));
                blend(&self.rows[i - 1], &self.rows[i], &alpha)
            })
            .collect();
        self.rows.insert(last + 1, self.rows[last].clone());
        self.rows[first..=last].clone_from_slice(&changed);
        self.knots.insert(k + 1, u.clone());
    }
}

pub(crate) struct KnotRefinementTransform {
    rows: Vec<(Vec<(usize, BigInt)>, BigInt)>,
}
impl KnotRefinementTransform {
    /// The caller has validated both bases, requests and final grid limits.
    pub(crate) fn new(
        old: &ExactKnotVector,
        new: &ExactKnotVector,
        changes: &[(R, usize)],
    ) -> Self {
        let mut work = Work::new(old);
        for (u, count) in changes {
            let copies = if old.is_periodic() {
                let period = &old.domain()[1] - &old.domain()[0];
                let last = if u == &old.domain()[0] { 3 } else { 2 };
                (-2..=last)
                    .map(|k| u + R::from_integer(k.into()) * &period)
                    .collect()
            } else {
                vec![u.clone()]
            };
            for u in copies {
                for _ in 0..*count {
                    work.insert(old.degree(), &u);
                }
            }
        }
        let offset = if old.is_periodic() {
            work.knots
                .windows(new.flat.len())
                .position(|w| w == new.flat)
                .expect("refined periodic support contains canonical extension")
        } else {
            debug_assert_eq!(work.knots, new.flat);
            0
        };
        let rows = work.rows[offset..offset + new.pole_count()]
            .iter()
            .map(|row| {
                let denominator = common_denominator(row.iter().map(|(_, x)| x));
                let numerators = row
                    .iter()
                    .map(|(i, x)| (*i, x.numer() * (&denominator / x.denom())))
                    .collect();
                (numerators, denominator)
            })
            .collect();
        Self { rows }
    }
    pub(crate) fn apply(&self, controls: &[[R; 4]]) -> Vec<[R; 4]> {
        let denominator = common_denominator(controls.iter().flatten());
        let integers: Vec<[BigInt; 4]> = controls
            .iter()
            .map(|p| std::array::from_fn(|c| p[c].numer() * (&denominator / p[c].denom())))
            .collect();
        self.rows
            .iter()
            .map(|(row, d)| {
                if row.len() == 1 && row[0].1 == *d {
                    return controls[row[0].0].clone();
                }
                let den = &denominator * d;
                std::array::from_fn(|c| {
                    R::new(
                        row.iter().map(|(i, x)| x * &integers[*i][c]).sum(),
                        den.clone(),
                    )
                })
            })
            .collect()
    }
}
