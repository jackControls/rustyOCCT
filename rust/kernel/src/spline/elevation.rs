//! Exact degree elevation, referencing BSplCLib::IncreaseDegree and its
//! Prautzsch rank-curve averaging. Sparse control maps are shared across all
//! transverse fields. An independent Cox coefficient solver checks the maps.
//!
//! Unclamped cropping is derived from the original active knot indices, rather
//! than OCCT's output-copy loop, which can change domains or emit bad counts.
//! Periodic data is extended over three periods and cropped by the complete
//! canonical extended knot sequence, preserving the original cyclic origin.
use super::refinement::{blend, Sparse};
use super::{integer, ExactKnotVector, KnotRefinementTransform};
use num_rational::BigRational as R;
use std::collections::BTreeMap;

fn expand(knots: &[R], mults: &[usize]) -> Vec<R> {
    knots
        .iter()
        .zip(mults)
        .flat_map(|(k, &m)| std::iter::repeat_n(k.clone(), m))
        .collect()
}

fn insert(p: usize, flat: &mut Vec<R>, rows: &mut Vec<Sparse>, u: &R) {
    let k = flat.partition_point(|x| x <= u) - 1;
    let mult = k + 1 - flat.partition_point(|x| x < u);
    let (first, last) = (k - p + 1, k - mult);
    let changed: Vec<_> = (first..=last)
        .map(|i| {
            let alpha = (u - &flat[i]) / (&flat[i + p] - &flat[i]);
            debug_assert!(alpha >= integer(0) && alpha <= integer(1));
            blend(&rows[i - 1], &rows[i], &alpha)
        })
        .collect();
    rows.insert(last + 1, rows[last].clone());
    rows[first..=last].clone_from_slice(&changed);
    flat.insert(k + 1, u.clone());
}

fn increment(p: usize, knots: &[R], mults: &mut [usize], rows: &[Sparse]) -> Vec<Sparse> {
    let mut sums = vec![BTreeMap::<usize, R>::new(); rows.len() + knots.len() - 1];
    for rank in 0..=p {
        let mut duplicated: Vec<_> = rows
            .iter()
            .enumerate()
            .flat_map(|(i, r)| std::iter::repeat_n(r.clone(), 1 + usize::from(i % (p + 1) == rank)))
            .collect();
        let mut cumulative = 0;
        let mut next = rank + 1;
        let mut missing = Vec::new();
        let selected: Vec<_> = knots
            .iter()
            .zip(mults.iter())
            .map(|(k, &m)| {
                cumulative += m;
                if cumulative >= next {
                    next += p + 1;
                    m + 1
                } else {
                    missing.push(k);
                    m
                }
            })
            .collect();
        let mut flat = expand(knots, &selected);
        debug_assert_eq!(flat.len(), duplicated.len() + p + 2);
        for k in missing {
            insert(p + 1, &mut flat, &mut duplicated, k);
        }
        debug_assert_eq!(duplicated.len(), sums.len());
        for (sum, row) in sums.iter_mut().zip(duplicated) {
            for (i, x) in row {
                *sum.entry(i).or_insert_with(|| integer(0)) += x;
            }
        }
    }
    for m in mults {
        *m += 1;
    }
    let divisor = integer(p + 1);
    sums.into_iter()
        .map(|sum| sum.into_iter().map(|(i, x)| (i, x / &divisor)).collect())
        .collect()
}

pub(crate) struct DegreeElevationTransform(KnotRefinementTransform);

impl DegreeElevationTransform {
    /// Both bases and the complete output grid have already been preflighted.
    pub(crate) fn new(old: &ExactKnotVector, new: &ExactKnotVector) -> Self {
        let p = old.degree();
        let n = old.pole_count();
        let (knots, mut mults, mut rows, prefix, suffix) = if old.is_periodic() {
            let period = &old.domain()[1] - &old.domain()[0];
            let e = p + 1 - old.multiplicities()[0];
            let mut knots = Vec::new();
            let mut mults = Vec::new();
            // n > p puts the active interval of these three periods strictly
            // around the entire fundamental period. Zero clamping at the two
            // outer ends therefore cannot affect it. After elevation, n' > q
            // and e = q + 1 - m' is unchanged: the canonical knot extension
            // lies strictly inside the neighboring periods, away from either
            // clamped end. Matching that full window preserves the cyclic
            // origin. Extra outer periods add work but no supporting basis.
            for turn in -1..=1 {
                let shift = R::from_integer(turn.into()) * &period;
                for (k, &m) in old
                    .knots()
                    .iter()
                    .zip(old.multiplicities())
                    .take(old.knots().len() - 1)
                {
                    knots.push(k + &shift);
                    mults.push(m);
                }
            }
            knots.push(&old.knots()[0] + integer(2) * &period);
            mults.push(old.multiplicities()[0]);
            let rows: Vec<_> = (0..3 * n - e)
                .map(|i| vec![((i + e) % n, integer(1))])
                .collect();
            (knots, mults, rows, e, e)
        } else {
            (
                old.knots().to_vec(),
                old.multiplicities().to_vec(),
                (0..n).map(|i| vec![(i, integer(1))]).collect(),
                p + 1 - old.multiplicities()[0],
                p + 1 - old.multiplicities().last().unwrap(),
            )
        };
        // Empty sparse rows are temporary zero homogeneous controls, never
        // exposed as positive-weight curve or surface data.
        rows.splice(0..0, std::iter::repeat_n(Vec::new(), prefix));
        rows.extend(std::iter::repeat_n(Vec::new(), suffix));
        mults[0] = p + 1;
        *mults.last_mut().unwrap() = p + 1;
        for degree in p..new.degree() {
            rows = increment(degree, &knots, &mut mults, &rows);
        }
        let offset = if old.is_periodic() {
            expand(&knots, &mults)
                .windows(new.flat.len())
                .position(|w| w == new.flat)
                .expect("elevated periodic support contains canonical extension")
        } else {
            let first = old.knots().binary_search(&old.domain()[0]).unwrap();
            let last = old.knots().binary_search(&old.domain()[1]).unwrap();
            let delta = new.degree() - p;
            let offset = prefix + delta * first;
            let tail = suffix + delta * (old.knots().len() - 1 - last);
            let flat = expand(&knots, &mults);
            debug_assert_eq!(&flat[offset..flat.len() - tail], new.flat);
            offset
        };
        let rows = &rows[offset..offset + new.pole_count()];
        debug_assert!(rows.iter().all(|r| !r.is_empty()));
        Self(KnotRefinementTransform::from_rows(rows))
    }

    pub(crate) fn apply(&self, controls: &[[R; 4]]) -> Vec<[R; 4]> {
        self.0.apply(controls)
    }
}
