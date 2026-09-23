//! Exact Boehm refinement and inverse insertion. OCCT references:
//! BSplCLib::InsertKnots, RemoveKnot, BoorScheme and AntiBoorScheme.
//! Periodic support is unrolled into five periods before local editing; the
//! desired extended knot sequence selects the canonical cyclic control block.
use super::{
    BSplineCurve3, BezierExtractionOptions, CurveEvaluation, DerivativeOrder, ExactBezierCurve3,
    ExactCurveEvaluation, KnotSide,
};
use crate::{
    spline::{self, ExactKnotVector},
    Error, Result,
};
use num_rational::BigRational as R;

fn integer(n: usize) -> R {
    R::from_integer(n.into())
}

/// An exact positive-weight rational B-spline. Editing preserves all four
/// homogeneous polynomial components, including inactive unclamped controls.
/// No operation rounds knots or control data to binary64.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExactBSplineCurve3 {
    basis: ExactKnotVector,
    controls: Vec<[R; 4]>,
}

impl BSplineCurve3 {
    /// Preserve the input binary64 atoms exactly, including weight scaling.
    pub fn to_exact(&self) -> ExactBSplineCurve3 {
        ExactBSplineCurve3 {
            basis: self.basis.to_exact(),
            controls: self
                .poles
                .iter()
                .zip(&self.weights)
                .map(|(&p, &w)| spline::homogeneous::<1>(p, w).map(|[x]| x))
                .collect(),
        }
    }
}
impl ExactBezierCurve3 {
    /// Exact clamped B-spline on the same rational parameter interval.
    pub fn to_bspline(&self) -> ExactBSplineCurve3 {
        let p = self.degree();
        ExactBSplineCurve3 {
            basis: ExactKnotVector::new(p, self.domain().to_vec(), vec![p + 1, p + 1])
                .expect("validated Bezier data"),
            controls: self.homogeneous_poles().to_vec(),
        }
    }
}

impl ExactBSplineCurve3 {
    /// Controls are `(w*x,w*y,w*z,w)`; every weight must be positive. All
    /// rationals are normalized and zero denominators are rejected.
    pub fn from_homogeneous(basis: ExactKnotVector, controls: Vec<[R; 4]>) -> Result<Self> {
        if controls.len() != basis.pole_count() {
            return Err(Error::InvalidCurve("pole/knot/multiplicity counts"));
        }
        let controls: Vec<_> = controls
            .into_iter()
            .map(|p| {
                let [x, y, z, w] = p;
                Ok([
                    spline::normalize(&x)?,
                    spline::normalize(&y)?,
                    spline::normalize(&z)?,
                    spline::normalize(&w)?,
                ])
            })
            .collect::<Result<_>>()?;
        if controls.iter().any(|p| p[3] <= integer(0)) {
            return Err(Error::InvalidCurve("weights must be positive"));
        }
        Ok(Self { basis, controls })
    }
    pub fn degree(&self) -> usize {
        self.basis.degree()
    }
    pub fn knot_vector(&self) -> &ExactKnotVector {
        &self.basis
    }
    pub fn knots(&self) -> &[R] {
        self.basis.knots()
    }
    pub fn multiplicities(&self) -> &[usize] {
        self.basis.multiplicities()
    }
    pub fn domain(&self) -> &[R; 2] {
        self.basis.domain()
    }
    pub fn is_periodic(&self) -> bool {
        self.basis.is_periodic()
    }
    pub fn homogeneous_poles(&self) -> &[[R; 4]] {
        &self.controls
    }

    pub(crate) fn span_polynomial(&self, span: usize) -> [Vec<R>; 4] {
        let controls = (span - self.degree()..=span)
            .map(|i| self.controls[self.basis.pole_index(i)].clone())
            .collect();
        spline::span_polynomial_from_knots(self.degree(), &self.basis.flat, span, controls)
    }

    /// Raise the total multiplicity. Zero or a smaller valid target is a no-op.
    /// The closed fundamental domain is required; periodic endpoint aliases
    /// refer to one seam. Invalid/oversized requests are errors, never snapped.
    pub fn insert_knot(&self, u: &R, target: usize) -> Result<Self> {
        self.refined(&[(u.clone(), target)])
    }

    /// Immutable, atomic batch refinement. Requests may be unordered; duplicate
    /// knots and periodic seam aliases use the greatest target. At most 4096
    /// requests and 4096 resulting poles, preflighted before local arithmetic.
    pub fn refined(&self, requests: &[(R, usize)]) -> Result<Self> {
        let (basis, changes) = self.basis.refinement_plan(requests)?;
        if changes.is_empty() {
            return Ok(self.clone());
        }
        let mut work = self.unrolled();
        for (u, count) in changes {
            let copies = if self.is_periodic() {
                let period = &self.domain()[1] - &self.domain()[0];
                let last = if u == self.domain()[0] { 3 } else { 2 };
                (-2..=last)
                    .map(|k| &u + R::from_integer(k.into()) * &period)
                    .collect()
            } else {
                vec![u]
            };
            for u in copies {
                for _ in 0..count {
                    work.insert(self.degree(), &u);
                }
            }
        }
        Ok(Self {
            controls: work.crop(&basis),
            basis,
        })
    }

    /// Exact homogeneous removal, with no tolerance. A valid request returns
    /// `None` when inverse insertion disagrees or produces a nonpositive weight.
    /// Missing knots, nonperiodic domain ends/exterior knots, and a resulting
    /// basis outside the supported family are errors. A target >= the existing
    /// multiplicity is a no-op. Removing a periodic seam to zero shifts the
    /// fundamental origin to the next distinct knot, retaining the period.
    pub fn remove_knot(&self, u: &R, target: usize) -> Result<Option<Self>> {
        let mut u = spline::normalize(u)?;
        if self.is_periodic() && u == self.domain()[1] {
            u = self.domain()[0].clone();
        }
        let index = self
            .knots()
            .binary_search(&u)
            .map_err(|_| Error::OutOfDomain("knot removal parameter"))?;
        if !self.is_periodic() && (u <= self.domain()[0] || u >= self.domain()[1]) {
            return Err(Error::OutOfDomain(
                "knot removal requires a strict interior knot",
            ));
        }
        let current = self.multiplicities()[index];
        if target >= current {
            return Ok(Some(self.clone()));
        }
        let mut knots = self.knots().to_vec();
        let mut mults = self.multiplicities().to_vec();
        if self.is_periodic() && index == 0 {
            if target == 0 {
                knots.pop();
                mults.pop();
                knots.remove(0);
                mults.remove(0);
                knots.push(&knots[0] + (&self.domain()[1] - &self.domain()[0]));
                mults.push(mults[0]);
            } else {
                mults[0] = target;
                *mults.last_mut().unwrap() = target;
            }
        } else if target == 0 {
            knots.remove(index);
            mults.remove(index);
        } else {
            mults[index] = target;
        }
        let basis = ExactKnotVector::build(self.degree(), knots, mults, self.is_periodic())?;
        let mut work = self.unrolled();
        let copies: Vec<_> = if self.is_periodic() {
            let period = &self.domain()[1] - &self.domain()[0];
            (-1..=2)
                .map(|k| &u + R::from_integer(k.into()) * &period)
                .collect()
        } else {
            vec![u]
        };
        for u in copies {
            for _ in target..current {
                if !work.remove(self.degree(), &u) {
                    return Ok(None);
                }
            }
        }
        let controls = work.crop(&basis);
        if controls.iter().any(|p| p[3] <= integer(0)) {
            return Ok(None);
        }
        Ok(Some(Self { basis, controls }))
    }

    pub fn exact_evaluate(
        &self,
        u: &R,
        order: DerivativeOrder,
        side: KnotSide,
    ) -> Result<ExactCurveEvaluation> {
        let count = order.count() + 1;
        let evaluate = |at: &spline::Parameter| {
            let poles = (at.span - self.degree()..=at.span)
                .map(|i| {
                    let p = &self.controls[self.basis.pole_index(i)];
                    std::array::from_fn(|c| {
                        std::array::from_fn(|r| if r == 0 { p[c].clone() } else { integer(0) })
                    })
                })
                .collect();
            let h = spline::de_boor_exact(
                self.degree(),
                &self.basis.flat,
                at,
                poles,
                [None, Some((0, 1)), Some((1, 2))],
                count,
            );
            spline::rationalize(h, [(0, 0), (1, 0), (2, 0)], count)
        };
        let sides = self.basis.locate(u, side, order.count())?;
        let values = evaluate(&sides[0]);
        if sides.get(1).is_some_and(|at| evaluate(at) != values) {
            return Err(Error::DiscontinuousDerivative);
        }
        Ok(ExactCurveEvaluation::from_values(values))
    }

    /// Minimal finite binary64 enclosures; exact values remain available when
    /// conversion would overflow. No extrapolation at nonperiodic ends.
    pub fn evaluate(
        &self,
        u: f64,
        order: DerivativeOrder,
        side: KnotSide,
    ) -> Result<CurveEvaluation> {
        crate::math::finite(u, "spline parameter")?;
        self.exact_evaluate(&spline::rational(u), order, side)?
            .enclosed()
    }
    pub fn bezier_arcs(&self) -> Result<Vec<ExactBezierCurve3>> {
        self.bezier_arcs_in(&self.domain()[0], &self.domain()[1])
    }
    pub fn bezier_arcs_in(&self, first: &R, last: &R) -> Result<Vec<ExactBezierCurve3>> {
        self.bezier_arcs_with_options(first, last, BezierExtractionOptions::default())
    }
    /// Rational intervals can span multiple periodic turns. Count the required
    /// arcs before traversal; zero max_arcs rejects every positive interval.
    pub fn bezier_arcs_with_options(
        &self,
        first: &R,
        last: &R,
        options: BezierExtractionOptions,
    ) -> Result<Vec<ExactBezierCurve3>> {
        let spans = self.basis.spans_in(first, last, options.max_arcs)?;
        Ok(spans
            .into_iter()
            .map(|span| {
                let controls: Vec<_> = (span.index - self.degree()..=span.index)
                    .map(|i| self.controls[self.basis.pole_index(i)].clone())
                    .collect();
                let length = &span.end - &span.start;
                let map = spline::BezierSpanTransform::from_knots(
                    self.degree(),
                    &self.basis.flat,
                    span.index,
                    &((&span.lower - &span.start) / &length),
                    &((&span.upper - &span.start) / length),
                );
                ExactBezierCurve3::from_homogeneous(map.apply(&controls), [span.lower, span.upper])
            })
            .collect())
    }

    fn unrolled(&self) -> Work {
        if !self.is_periodic() {
            return Work {
                knots: self.basis.flat.clone(),
                controls: self.controls.clone(),
            };
        }
        let n = self.controls.len();
        let extension = self.degree() + 1 - self.multiplicities()[0];
        let base: Vec<_> = self
            .knots()
            .iter()
            .zip(self.multiplicities())
            .take(self.knots().len() - 1)
            .flat_map(|(k, &m)| std::iter::repeat_n(k, m))
            .collect();
        let period = &self.domain()[1] - &self.domain()[0];
        let start = -2 * (n as isize);
        let count = 5 * n + extension;
        let knots = (start..start + (count + self.degree() + 1) as isize)
            .map(|i| {
                let j = i - extension as isize;
                base[j.rem_euclid(n as isize) as usize]
                    + R::from_integer(j.div_euclid(n as isize).into()) * &period
            })
            .collect();
        let controls = (0..count).map(|i| self.controls[i % n].clone()).collect();
        Work { knots, controls }
    }
}

struct Work {
    knots: Vec<R>,
    controls: Vec<[R; 4]>,
}
impl Work {
    fn insert(&mut self, degree: usize, u: &R) {
        let k = self.knots.partition_point(|x| x <= u) - 1;
        let mult = k + 1 - self.knots.partition_point(|x| x < u);
        let (first, last) = (k - degree + 1, k - mult);
        let row: Vec<_> = (first..=last)
            .map(|i| {
                let alpha = (u - &self.knots[i]) / (&self.knots[i + degree] - &self.knots[i]);
                debug_assert!(alpha >= integer(0) && alpha <= integer(1));
                spline::blend_homogeneous(&self.controls[i - 1], &self.controls[i], &alpha)
            })
            .collect();
        self.controls.insert(last + 1, self.controls[last].clone());
        self.controls[first..=last].clone_from_slice(&row);
        self.knots.insert(k + 1, u.clone());
    }
    fn remove(&mut self, degree: usize, u: &R) -> bool {
        let last_equal = self.knots.partition_point(|x| x <= u) - 1;
        self.knots.remove(last_equal);
        let k = self.knots.partition_point(|x| x <= u) - 1;
        let mult = k + 1 - self.knots.partition_point(|x| x < u);
        let (first, last) = (k - degree + 1, k - mult);
        let mut previous = self.controls[first - 1].clone();
        let mut row = Vec::with_capacity(last - first + 1);
        for i in first..=last {
            let alpha = (u - &self.knots[i]) / (&self.knots[i + degree] - &self.knots[i]);
            debug_assert!(alpha > integer(0) && alpha < integer(1));
            let complement = integer(1) - &alpha;
            previous = std::array::from_fn(|c| {
                (&self.controls[i][c] - &complement * &previous[c]) / &alpha
            });
            row.push(previous.clone());
        }
        if previous != self.controls[last + 1] {
            return false;
        }
        self.controls[first..last].clone_from_slice(&row[..row.len() - 1]);
        self.controls.remove(last);
        true
    }
    fn crop(self, basis: &ExactKnotVector) -> Vec<[R; 4]> {
        if !basis.is_periodic() {
            debug_assert_eq!(self.knots, basis.flat);
            return self.controls;
        }
        let offset = self
            .knots
            .windows(basis.flat.len())
            .position(|w| w == basis.flat)
            .expect("edited periodic support contains the canonical knot extension");
        self.controls[offset..offset + basis.pole_count()].to_vec()
    }
}
