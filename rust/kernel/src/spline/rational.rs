//! Rational knot data for exact editing, including sub-binary64 knot spacing.
use super::{integer, zero, KnotSpan, KnotVector, Parameter, MAX_DEGREE, MAX_POLES};
use crate::{curve::KnotSide, Error, Result};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::collections::BTreeMap;

/// OCCT knot and periodic pole-order conventions, with exact rational atoms.
/// Degree 1..=25, degree < pole count <=4096. There is no knot snapping.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactKnotVector {
    pub(crate) degree: usize,
    pub(crate) pole_count: usize,
    pub(crate) periodic: bool,
    pub(crate) knots: Vec<R>,
    pub(crate) multiplicities: Vec<usize>,
    pub(crate) flat: Vec<R>,
    pub(crate) domain: [R; 2],
}

pub(crate) fn normalize(x: &R) -> Result<R> {
    if x.denom() == &BigInt::from(0) {
        return Err(Error::InvalidSpline("rational value has zero denominator"));
    }
    Ok(R::new(x.numer().clone(), x.denom().clone()))
}

impl KnotVector {
    /// Preserve all binary64 knot atoms and the periodic extension exactly.
    pub fn to_exact(&self) -> ExactKnotVector {
        ExactKnotVector {
            degree: self.degree,
            pole_count: self.pole_count,
            periodic: self.periodic,
            knots: self.knots.iter().map(|&k| super::rational(k)).collect(),
            multiplicities: self.multiplicities.clone(),
            flat: self.flat.clone(),
            domain: [
                super::rational(self.domain.0),
                super::rational(self.domain.1),
            ],
        }
    }
}

impl ExactKnotVector {
    /// Normalize rationals; reject zero denominators, unordered/duplicate knots,
    /// invalid multiplicities, empty domains and counts outside the family.
    pub fn new(degree: usize, knots: Vec<R>, multiplicities: Vec<usize>) -> Result<Self> {
        Self::build(degree, knots, multiplicities, false)
    }
    /// Equal endpoint multiplicities in 1..=degree; cyclic poles. Evaluation
    /// reduces parameters modulo the period with exact rational arithmetic.
    pub fn new_periodic(degree: usize, knots: Vec<R>, multiplicities: Vec<usize>) -> Result<Self> {
        Self::build(degree, knots, multiplicities, true)
    }
    pub(crate) fn build(
        degree: usize,
        knots: Vec<R>,
        multiplicities: Vec<usize>,
        periodic: bool,
    ) -> Result<Self> {
        if degree == 0 || degree > MAX_DEGREE {
            return Err(Error::InvalidSpline("degree must be in 1..=25"));
        }
        if knots.len() > MAX_POLES + MAX_DEGREE + 1 {
            return Err(Error::LimitExceeded("spline knot data"));
        }
        if knots.len() < 2 || knots.len() != multiplicities.len() {
            return Err(Error::InvalidSpline("knot/multiplicity counts"));
        }
        let knots: Vec<_> = knots.iter().map(normalize).collect::<Result<_>>()?;
        if knots.windows(2).any(|k| k[0] >= k[1]) {
            return Err(Error::InvalidSpline("knots must be strictly increasing"));
        }
        let mut total: usize = 0;
        for (i, &m) in multiplicities.iter().enumerate() {
            let end = i == 0 || i + 1 == knots.len();
            if m == 0 || m > degree + usize::from(end && !periodic) {
                return Err(Error::InvalidSpline("knot multiplicity"));
            }
            total += m;
        }
        if periodic && multiplicities[0] != multiplicities[knots.len() - 1] {
            return Err(Error::InvalidSpline(
                "periodic end multiplicities must agree",
            ));
        }
        let pole_count = total
            .checked_sub(if periodic {
                multiplicities[0]
            } else {
                degree + 1
            })
            .ok_or(Error::InvalidSpline("sum of knot multiplicities"))?;
        if pole_count > MAX_POLES {
            return Err(Error::LimitExceeded("spline control data"));
        }
        if pole_count <= degree {
            return Err(Error::InvalidSpline("more poles than degree required"));
        }
        let mut flat: Vec<_> = knots
            .iter()
            .zip(&multiplicities)
            .flat_map(|(k, &m)| std::iter::repeat_n(k.clone(), m))
            .collect();
        let domain = if periodic {
            let domain = [knots[0].clone(), knots.last().unwrap().clone()];
            let period = &domain[1] - &domain[0];
            let extension = degree + 1 - multiplicities[0];
            let before: Vec<_> = flat[pole_count - extension..pole_count]
                .iter()
                .map(|k| k - &period)
                .collect();
            let after: Vec<_> = flat[multiplicities[0]..multiplicities[0] + extension]
                .iter()
                .map(|k| k + &period)
                .collect();
            flat = before.into_iter().chain(flat).chain(after).collect();
            domain
        } else {
            [flat[degree].clone(), flat[pole_count].clone()]
        };
        if domain[0] >= domain[1] {
            return Err(Error::InvalidSpline("empty parameter domain"));
        }
        Ok(Self {
            degree,
            pole_count,
            periodic,
            knots,
            multiplicities,
            flat,
            domain,
        })
    }
    pub fn degree(&self) -> usize {
        self.degree
    }
    pub fn pole_count(&self) -> usize {
        self.pole_count
    }
    pub fn knots(&self) -> &[R] {
        &self.knots
    }
    pub fn multiplicities(&self) -> &[usize] {
        &self.multiplicities
    }
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }
    pub fn domain(&self) -> &[R; 2] {
        &self.domain
    }
    /// Preflight degree and pole limits without evaluating any control data.
    /// An unclamped axis loses only the exterior flat knots prescribed by its
    /// original active knot indices; the active parameter interval is unchanged.
    pub(crate) fn elevation_plan(&self, degree: usize) -> Result<Self> {
        if degree < self.degree || degree > MAX_DEGREE {
            return Err(Error::InvalidSpline("degree elevation target"));
        }
        let delta = degree - self.degree;
        if delta == 0 {
            return Ok(self.clone());
        }
        let first = self.knots.binary_search(&self.domain[0]).unwrap();
        let last = self.knots.binary_search(&self.domain[1]).unwrap();
        if self.pole_count + delta * (last - first) > MAX_POLES {
            return Err(Error::LimitExceeded("spline control data"));
        }
        let mut left = if self.periodic { 0 } else { delta * first };
        let mut right = if self.periodic {
            0
        } else {
            delta * (self.knots.len() - 1 - last)
        };
        let mut pairs: Vec<_> = self
            .knots
            .iter()
            .zip(&self.multiplicities)
            .filter_map(|(k, &m)| {
                let m = m + delta;
                let take = left.min(m);
                left -= take;
                (m > take).then(|| (k.clone(), m - take))
            })
            .collect();
        for (_, m) in pairs.iter_mut().rev() {
            let take = right.min(*m);
            right -= take;
            *m -= take;
        }
        let (knots, mults) = pairs.into_iter().filter(|(_, m)| *m != 0).unzip();
        let result = Self::build(degree, knots, mults, self.periodic)?;
        debug_assert_eq!(result.domain, self.domain);
        debug_assert_eq!(result.pole_count, self.pole_count + delta * (last - first));
        Ok(result)
    }
    /// Plan both axis and grid limits without touching any control values.
    pub(crate) fn refinement_plan(
        &self,
        requests: &[(R, usize)],
    ) -> Result<(Self, Vec<(R, usize)>)> {
        if requests.len() > MAX_POLES {
            return Err(Error::LimitExceeded("knot refinement requests"));
        }
        let mut requested = BTreeMap::<R, usize>::new();
        for (u, target) in requests {
            let mut u = normalize(u)?;
            if u < self.domain()[0] || u > self.domain()[1] {
                return Err(Error::OutOfDomain("knot refinement parameter"));
            }
            let physical_end =
                !self.is_periodic() && (u == self.knots()[0] || &u == self.knots().last().unwrap());
            if *target > self.degree() + usize::from(physical_end) {
                return Err(Error::InvalidSpline("knot refinement multiplicity"));
            }
            if self.is_periodic() && u == self.domain()[1] {
                u = self.domain()[0].clone();
            }
            requested
                .entry(u)
                .and_modify(|m| *m = (*m).max(*target))
                .or_insert(*target);
        }
        let mut map: BTreeMap<_, _> = self
            .knots()
            .iter()
            .cloned()
            .zip(self.multiplicities().iter().copied())
            .collect();
        let mut changes = Vec::new();
        for (u, target) in requested {
            let current = map.get(&u).copied().unwrap_or(0);
            if target > current {
                map.insert(u.clone(), target);
                if self.is_periodic() && u == self.domain()[0] {
                    map.insert(self.domain()[1].clone(), target);
                }
                changes.push((u, target - current));
            }
        }
        if changes.is_empty() {
            return Ok((self.clone(), changes));
        }
        let (knots, multiplicities) = map.into_iter().unzip();
        let basis =
            ExactKnotVector::build(self.degree(), knots, multiplicities, self.is_periodic())?;
        Ok((basis, changes))
    }
    pub(crate) fn pole_index(&self, index: usize) -> usize {
        if self.periodic {
            index % self.pole_count
        } else {
            index
        }
    }
    pub(crate) fn locate(&self, u: &R, side: KnotSide, order: usize) -> Result<Vec<Parameter>> {
        let mut u = normalize(u)?;
        let [start, end] = &self.domain;
        if !self.periodic
            && (u < *start
                || u > *end
                || (u == *start && side == KnotSide::Left)
                || (u == *end && side == KnotSide::Right))
        {
            return Err(Error::OutOfDomain("spline parameter or requested side"));
        }
        if self.periodic {
            let period = end - start;
            u -= ((&u - start) / &period).floor() * period;
        }
        let seam = self.periodic && u == *start;
        let parameter = |value: R, left: bool| {
            let span = self
                .flat
                .partition_point(|k| if left { k < &value } else { k <= &value })
                - 1;
            Parameter { value, span }
        };
        let mut sides = vec![if seam && side == KnotSide::Left {
            parameter(end.clone(), true)
        } else {
            parameter(
                u.clone(),
                side == KnotSide::Left || (!self.periodic && u == *end),
            )
        }];
        if side == KnotSide::Automatic && order > 0 {
            let multiplicity = if seam {
                self.multiplicities[0]
            } else if u > *start && u < *end {
                self.flat.partition_point(|k| k <= &u) - self.flat.partition_point(|k| k < &u)
            } else {
                0
            };
            if multiplicity > 0 && order > self.degree - multiplicity {
                sides.push(parameter(if seam { end.clone() } else { u }, true));
            }
        }
        Ok(sides)
    }
    pub(crate) fn spans_in(&self, first: &R, last: &R, max_spans: usize) -> Result<Vec<KnotSpan>> {
        let (lower, upper) = (normalize(first)?, normalize(last)?);
        if lower >= upper {
            return Err(Error::OutOfDomain(
                "curve interval must have positive length",
            ));
        }
        if !self.periodic && (lower < self.domain[0] || upper > self.domain[1]) {
            return Err(Error::OutOfDomain("curve interval"));
        }
        let base: Vec<_> = self
            .knots
            .windows(2)
            .filter(|k| k[0] >= self.domain[0] && k[1] <= self.domain[1])
            .map(|k| (&k[0], &k[1], self.flat.partition_point(|x| x <= &k[0]) - 1))
            .collect();
        let period = &self.domain[1] - &self.domain[0];
        let mut count = BigInt::from(0);
        for &(a, b, _) in &base {
            count += if self.periodic {
                (((&upper - a) / &period).ceil() - ((&lower - b) / &period).floor() - integer(1))
                    .to_integer()
            } else {
                BigInt::from(usize::from(a < &upper && b > &lower))
            };
            if count > BigInt::from(max_spans) {
                return Err(Error::ComputationLimit("spline span traversal"));
            }
        }
        let mut offset = if self.periodic {
            ((&lower - &self.domain[0]) / &period).floor() * &period
        } else {
            zero()
        };
        let mut result = Vec::new();
        loop {
            for &(a, b, index) in &base {
                let (start, end) = (a + &offset, b + &offset);
                if start >= upper {
                    break;
                }
                if end <= lower {
                    continue;
                }
                result.push(KnotSpan {
                    index,
                    lower: start.clone().max(lower.clone()),
                    upper: end.clone().min(upper.clone()),
                    start,
                    end,
                });
            }
            if !self.periodic {
                break;
            }
            offset += &period;
            if &self.domain[0] + &offset >= upper {
                break;
            }
        }
        debug_assert_eq!(BigInt::from(result.len()), count);
        Ok(result)
    }
}
