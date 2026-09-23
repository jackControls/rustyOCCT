//! Complete minimum-distance parameter sets on closed rational B-spline ranges.
//! Reference: OCCT Extrema/ExtremaPC stationary equations, span boundaries and
//! endpoint handling. Exact polynomial isolation replaces sampled root searches.
use crate::polynomial::{
    real::{
        self,
        image::{image_roots, ImageBudget},
        Budget, IntPolynomial,
    },
    AlgebraicRoot, RootIsolationOptions,
};
use crate::{
    curve::{DerivativeOrder, KnotSide},
    interval,
    math::finite,
    spline, BSplineCurve3, Bounds3, Error, ExactBSplineCurve3, Point3, Result, ScalarInterval,
};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::{
    cmp::Ordering::{self, Equal, Greater, Less},
    sync::Arc,
};

#[cfg(test)]
mod tests;

/// Independent deterministic work limits. These are not geometric tolerances
/// or hard CPU/memory deadlines. Image limits cover coefficient updates and
/// intermediate integer sizes in the exact tie-resolution construction; root
/// subdivision has its own shared limit. Exhaustion returns no partial result.
#[derive(Debug, Clone, Copy)]
pub struct SplineProximityOptions {
    pub root_isolation: RootIsolationOptions,
    pub max_spans: usize,
    pub max_candidates: usize,
    pub max_image_coefficient_updates: usize,
    pub max_image_coefficient_bits: u64,
}
impl Default for SplineProximityOptions {
    fn default() -> Self {
        Self {
            root_isolation: RootIsolationOptions::default(),
            max_spans: spline::MAX_POLES,
            max_candidates: 65_536,
            max_image_coefficient_updates: 2_000_000,
            max_image_coefficient_bits: 65_536,
        }
    }
}

/// Every globally closest parameter on the queried closed range. Distinct
/// parameter aliases remain distinct even when their geometric points coincide.
/// Whole minimizing intervals replace their endpoints in the isolated list.
#[derive(Debug, Clone)]
pub struct SplineClosestSet {
    points: Vec<SplineClosestPoint>,
    intervals: Vec<SplineClosestInterval>,
    distance: Distance,
}
impl SplineClosestSet {
    pub fn points(&self) -> &[SplineClosestPoint] {
        &self.points
    }
    pub fn intervals(&self) -> &[SplineClosestInterval] {
        &self.intervals
    }
    pub fn compare_squared_distance(&self, value: f64) -> Result<Ordering> {
        finite(value, "squared distance comparison")?;
        Ok(self.distance.compare_rational(&real::rat(value)))
    }
    pub fn compare_squared_distance_exact(&self, value: &R) -> Result<Ordering> {
        Ok(self.distance.compare_rational(&spline::normalize(value)?))
    }
    pub fn compare_distance(&self, value: f64) -> Result<Ordering> {
        finite(value, "distance comparison")?;
        Ok(if value < 0. {
            Greater
        } else {
            let value = real::rat(value);
            self.distance.compare_rational(&(&value * &value))
        })
    }
    pub fn squared_distance_bounds(&self) -> Result<ScalarInterval> {
        let distance = self.distance.for_view();
        interval::enclose(
            |x| distance.compare_rational(&real::rat(x)),
            "squared curve distance",
        )
    }
    pub fn distance_bounds(&self) -> Result<ScalarInterval> {
        let distance = self.distance.for_view();
        interval::enclose(
            |x| {
                if x < 0. {
                    Greater
                } else {
                    let x = real::rat(x);
                    distance.compare_rational(&(&x * &x))
                }
            },
            "curve distance",
        )
    }
    /// Exact comparison between the two globally minimal distances. The
    /// supplied algebraic limits apply to new work in this comparison.
    pub fn distance_cmp(&self, other: &Self, options: SplineProximityOptions) -> Result<Ordering> {
        Context::new(options).compare(&mut self.distance.clone(), &mut other.distance.clone())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplineClosestInterval {
    parameters: [R; 2],
}
impl SplineClosestInterval {
    pub fn parameters(&self) -> &[R; 2] {
        &self.parameters
    }
    pub fn parameter_bounds(&self) -> Result<[ScalarInterval; 2]> {
        let bound = |i: usize| {
            interval::enclose(
                |x| self.parameters[i].cmp(&real::rat(x)),
                "minimum interval endpoint",
            )
        };
        Ok([bound(0)?, bound(1)?])
    }
}

/// Exact closest point, retaining its algebraic parameter and rational curve
/// equation. Floating bounds are optional views, not the point's identity.
#[derive(Debug, Clone)]
pub struct SplineClosestPoint {
    at: Distance,
}
impl SplineClosestPoint {
    pub fn compare_parameter(&self, value: &R) -> Result<Ordering> {
        let value = spline::normalize(value)?;
        Ok(self
            .at
            .root
            .compare_rational(&((value - &self.at.span.start) / &self.at.span.length)))
    }
    pub fn rational_parameter(&self) -> Option<R> {
        self.at
            .root
            .rational_value()
            .map(|x| &self.at.span.start + x * &self.at.span.length)
    }
    pub fn parameter_bounds(&self) -> Result<ScalarInterval> {
        interval::enclose(
            |x| {
                self.at.root.compare_rational(
                    &((real::rat(x) - &self.at.span.start) / &self.at.span.length),
                )
            },
            "closest curve parameter",
        )
    }
    pub fn compare_coordinate(&self, component: usize, value: &R) -> Result<Ordering> {
        let coordinate = self
            .at
            .span
            .homogeneous
            .get(component)
            .filter(|_| component < 3)
            .ok_or(Error::OutOfDomain("closest curve coordinate"))?;
        let value = spline::normalize(value)?;
        Ok(self
            .at
            .root
            .sign_polynomial(&IntPolynomial::from_rationals(&subtract(
                coordinate,
                &self.at.span.homogeneous[3],
                &value,
            ))))
    }
    pub fn coordinate_bounds(&self) -> Result<[ScalarInterval; 3]> {
        let point = Self {
            at: self.at.for_view(),
        };
        let bound = |i: usize| {
            interval::enclose(
                |x| {
                    point
                        .compare_coordinate(i, &real::rat(x))
                        .expect("validated component and finite value")
                },
                "closest curve point",
            )
        };
        Ok([bound(0)?, bound(1)?, bound(2)?])
    }
    pub fn point_bounds(&self) -> Result<Bounds3> {
        let [x, y, z] = self.coordinate_bounds()?;
        Ok(Bounds3 {
            min: Point3::new(x.lower(), y.lower(), z.lower()),
            max: Point3::new(x.upper(), y.upper(), z.upper()),
        })
    }
}

pub fn closest_points_on_spline(curve: &BSplineCurve3, point: Point3) -> Result<SplineClosestSet> {
    let (a, b) = curve.domain();
    closest_points_on_spline_in(curve, point, a, b, SplineProximityOptions::default())
}
/// Include both endpoints; periodic curves may cover multiple turns. The result
/// remains exact even when some coordinates cannot be enclosed by finite f64.
pub fn closest_points_on_spline_in(
    curve: &BSplineCurve3,
    point: Point3,
    first: f64,
    last: f64,
    options: SplineProximityOptions,
) -> Result<SplineClosestSet> {
    for x in point.to_array() {
        finite(x, "projection point")?;
    }
    finite(first, "projection range")?;
    finite(last, "projection range")?;
    closest_points_on_exact_spline_in(
        &curve.to_exact(),
        &point.to_array().map(real::rat),
        &real::rat(first),
        &real::rat(last),
        options,
    )
}
pub fn closest_points_on_exact_spline(
    curve: &ExactBSplineCurve3,
    point: &[R; 3],
) -> Result<SplineClosestSet> {
    closest_points_on_exact_spline_in(
        curve,
        point,
        &curve.domain()[0],
        &curve.domain()[1],
        SplineProximityOptions::default(),
    )
}
/// Minimize squared Euclidean distance over the complete closed range. An exact
/// affine-hull bound can certify all span minima directly; otherwise every
/// stationary root and knot/range endpoint is considered. No speed
/// normalization discards singular stationary points.
/// Degree 1..25 and positive weights follow the curve's existing contract.
pub fn closest_points_on_exact_spline_in(
    curve: &ExactBSplineCurve3,
    point: &[R; 3],
    first: &R,
    last: &R,
    options: SplineProximityOptions,
) -> Result<SplineClosestSet> {
    let first = spline::normalize(first)?;
    let last = spline::normalize(last)?;
    let point = [
        spline::normalize(&point[0])?,
        spline::normalize(&point[1])?,
        spline::normalize(&point[2])?,
    ];
    if first > last {
        return Err(Error::OutOfDomain("reversed projection range"));
    }
    let mut context = Context::new(options);
    let mut winners: Vec<Candidate> = Vec::new();
    if first == last {
        if options.max_spans == 0 {
            return Err(Error::ComputationLimit("spline span traversal"));
        }
        let position = curve
            .exact_evaluate(&first, DerivativeOrder::Position, KnotSide::Automatic)?
            .position()
            .coordinates()
            .clone();
        let h = [
            vec![position[0].clone()],
            vec![position[1].clone()],
            vec![position[2].clone()],
            vec![one()],
        ];
        let span = Arc::new(Span::new(first, one(), h, &point, one(), one()));
        context.consider(
            &mut winners,
            Candidate {
                at: Distance {
                    span,
                    root: AlgebraicRoot::rational(zero()),
                    known_squared_distance: None,
                },
                interval: None,
            },
        )?;
    } else {
        let spans = curve
            .knot_vector()
            .spans_in(&first, &last, options.max_spans)?;
        let mut zero_found = false;
        for (i, at) in spans.into_iter().enumerate() {
            let length = &at.end - &at.start;
            let lower = (&at.lower - &at.start) / &length;
            let upper = (&at.upper - &at.start) / &length;
            let h = curve.span_polynomial(at.index);
            let weights = (at.index - curve.degree()..=at.index)
                .map(|j| &curve.homogeneous_poles()[curve.knot_vector().pole_index(j)][3]);
            let min_weight = weights.clone().min().unwrap().clone();
            let max_weight = weights.max().unwrap().clone();
            let span = Arc::new(Span::new(
                at.start, length, h, &point, min_weight, max_weight,
            ));
            // Orthogonal projection onto the exact affine hull proves a
            // squared-distance lower bound. Every preimage of that projection
            // attains the bound; no other point on the span can improve it.
            let (projection, residual) = span.affine_hull_projection();
            let lower_bound = dot(&projection, &projection);
            if zero_found && lower_bound > zero() {
                continue;
            }
            let mut common = IntPolynomial::from_rationals(&residual[0]);
            for delta in &residual[1..] {
                if common.is_constant() && !common.is_zero() {
                    break;
                }
                common = common.gcd(&IntPolynomial::from_rationals(delta));
            }
            if !common.is_zero() {
                let zeros =
                    real::isolate(&common, lower.clone(), upper.clone(), &mut context.roots)?;
                if !zeros.is_empty() {
                    zero_found |= lower_bound == zero();
                    for mut root in zeros {
                        root.refine_for_signs(16);
                        context.consider(
                            &mut winners,
                            Candidate {
                                at: Distance {
                                    span: span.clone(),
                                    root,
                                    known_squared_distance: Some(lower_bound.clone()),
                                },
                                interval: None,
                            },
                        )?;
                    }
                    continue;
                }
                if zero_found {
                    continue;
                }
            } else {
                zero_found |= lower_bound == zero();
            }
            let f = subtract(
                &product(&derivative(&span.numerator), &span.homogeneous[3]),
                &product(&span.numerator, &derivative(&span.homogeneous[3])),
                &integer(2),
            );
            let f = IntPolynomial::from_rationals(&f);
            if f.is_zero() {
                context.consider(
                    &mut winners,
                    Candidate {
                        at: Distance {
                            span,
                            root: AlgebraicRoot::rational(lower),
                            known_squared_distance: None,
                        },
                        interval: Some([at.lower, at.upper]),
                    },
                )?;
                continue;
            }
            if i == 0 {
                context.consider(
                    &mut winners,
                    Candidate {
                        at: Distance {
                            span: span.clone(),
                            root: AlgebraicRoot::rational(lower.clone()),
                            known_squared_distance: None,
                        },
                        interval: None,
                    },
                )?;
            }
            for mut root in real::isolate(&f, lower.clone(), upper.clone(), &mut context.roots)? {
                if root.compare_rational(&lower) == Equal || root.compare_rational(&upper) == Equal
                {
                    continue;
                }
                root.refine_for_signs(16);
                context.consider(
                    &mut winners,
                    Candidate {
                        at: Distance {
                            span: span.clone(),
                            root,
                            known_squared_distance: None,
                        },
                        interval: None,
                    },
                )?;
            }
            context.consider(
                &mut winners,
                Candidate {
                    at: Distance {
                        span,
                        root: AlgebraicRoot::rational(upper),
                        known_squared_distance: None,
                    },
                    interval: None,
                },
            )?;
        }
    }
    let distance = winners
        .first()
        .expect("continuous curve on nonempty compact range has a minimum")
        .at
        .clone();
    let mut intervals: Vec<SplineClosestInterval> = Vec::new();
    for range in winners.iter().filter_map(|c| c.interval.as_ref()) {
        if let Some(last) = intervals.last_mut().filter(|p| p.parameters[1] == range[0]) {
            last.parameters[1] = range[1].clone();
        } else {
            intervals.push(SplineClosestInterval {
                parameters: range.clone(),
            });
        }
    }
    let mut points: Vec<SplineClosestPoint> = Vec::new();
    for candidate in winners.into_iter().filter(|c| c.interval.is_none()) {
        let point = SplineClosestPoint { at: candidate.at };
        let covered = intervals.iter().any(|range| {
            point.compare_parameter(&range.parameters[0]).unwrap() != Less
                && point.compare_parameter(&range.parameters[1]).unwrap() != Greater
        });
        // A local minimum can be emitted by both adjacent closed spans.
        // Only their exact shared rational boundary can be duplicated; aliases
        // from different turns remain distinct parameters.
        let duplicate = point.rational_parameter().is_some_and(|u| {
            points
                .last()
                .is_some_and(|last| last.compare_parameter(&u).unwrap() == Equal)
        });
        if !covered && !duplicate {
            points.push(point);
        }
    }
    Ok(SplineClosestSet {
        points,
        intervals,
        distance,
    })
}

struct Candidate {
    at: Distance,
    interval: Option<[R; 2]>,
}
#[derive(Debug, Clone)]
struct Distance {
    span: Arc<Span>,
    root: AlgebraicRoot,
    // Set only after common residual roots certify the affine lower bound.
    known_squared_distance: Option<R>,
}
#[derive(Debug)]
struct Span {
    start: R,
    length: R,
    homogeneous: [Vec<R>; 4],
    delta: [Vec<R>; 3],
    numerator: Vec<R>,
    denominator: Vec<R>,
    min_weight: R,
    max_weight: R,
}
impl Span {
    /// The query is the origin in `delta` coordinates. Return the orthogonal
    /// projection of that origin onto the curve's affine hull, and the
    /// homogeneous displacement from that projection. All arithmetic is exact.
    fn affine_hull_projection(&self) -> ([R; 3], [Vec<R>; 3]) {
        let w = &self.homogeneous[3];
        let anchor: [R; 3] =
            std::array::from_fn(|i| self.delta[i].first().cloned().unwrap_or_else(zero) / &w[0]);
        let degree = self
            .delta
            .iter()
            .map(Vec::len)
            .chain([w.len()])
            .max()
            .unwrap();
        let mut direction: Option<[R; 3]> = None;
        let mut normal: Option<[R; 3]> = None;
        for j in 1..degree {
            let v: [R; 3] = std::array::from_fn(|i| {
                self.delta[i].get(j).cloned().unwrap_or_else(zero)
                    - &anchor[i] * w.get(j).cloned().unwrap_or_else(zero)
            });
            if let Some(n) = &normal {
                if dot(n, &v) != zero() {
                    // Full three-dimensional hull: its lower bound is zero.
                    return (std::array::from_fn(|_| zero()), self.delta.clone());
                }
            } else if let Some(u) = &direction {
                let n: [R; 3] = std::array::from_fn(|i| {
                    &u[(i + 1) % 3] * &v[(i + 2) % 3] - &u[(i + 2) % 3] * &v[(i + 1) % 3]
                });
                if n.iter().any(|x| x != &zero()) {
                    normal = Some(n);
                }
            } else if v.iter().any(|x| x != &zero()) {
                direction = Some(v);
            }
        }
        let projection = if let Some(n) = normal {
            let scale = dot(&anchor, &n) / dot(&n, &n);
            n.map(|x| x * &scale)
        } else if let Some(u) = direction {
            let scale = dot(&anchor, &u) / dot(&u, &u);
            std::array::from_fn(|i| &anchor[i] - &u[i] * &scale)
        } else {
            anchor
        };
        let residual = std::array::from_fn(|i| subtract(&self.delta[i], w, &projection[i]));
        (projection, residual)
    }

    fn new(
        start: R,
        length: R,
        homogeneous: [Vec<R>; 4],
        point: &[R; 3],
        min_weight: R,
        max_weight: R,
    ) -> Self {
        let delta: [Vec<R>; 3] =
            std::array::from_fn(|i| subtract(&homogeneous[i], &homogeneous[3], &point[i]));
        let numerator = delta
            .iter()
            .fold(vec![], |a, y| subtract(&a, &product(y, y), &integer(-1)));
        let denominator = product(&homogeneous[3], &homogeneous[3]);
        Self {
            start,
            length,
            homogeneous,
            delta,
            numerator,
            denominator,
            min_weight,
            max_weight,
        }
    }
}
impl Distance {
    fn for_view(&self) -> Self {
        let mut refined = self.clone();
        // Share this bounded refinement among the view's many exact threshold
        // comparisons. A loose isolator otherwise repeats expensive signed
        // algebraic queries for well-separated values. No approximate result
        // is accepted: undecided signs still take the complete exact path.
        refined.root.refine_for_signs(128);
        refined
    }
    fn rational_value(&self) -> Option<R> {
        self.known_squared_distance.clone().or_else(|| {
            self.root
                .rational_value()
                .map(|x| evaluate(&self.span.numerator, x) / evaluate(&self.span.denominator, x))
        })
    }
    fn compare_rational(&self, value: &R) -> Ordering {
        if value < &zero() {
            return Greater;
        }
        if let Some(x) = self.rational_value() {
            return x.cmp(value);
        }
        // A sum of real squares vanishes exactly when every component does.
        // Test the degree-p components instead of the degree-2p squared sum.
        if value == &zero() {
            return if self.span.delta.iter().all(|p| {
                self.root
                    .vanishes_polynomial(&IntPolynomial::from_rationals(p))
            }) {
                Equal
            } else {
                Greater
            };
        }
        self.root
            .sign_polynomial(&IntPolynomial::from_rationals(&subtract(
                &self.span.numerator,
                &self.span.denominator,
                value,
            )))
    }
    fn range(&self) -> (R, R) {
        let (a, b) = self.root.isolator();
        let (mut lo, mut hi) = (zero(), zero());
        for delta in &self.span.delta {
            let (x, y) = range(delta, a, b);
            let xx = &x * &x;
            let yy = &y * &y;
            lo += if x <= zero() && y >= zero() {
                zero()
            } else {
                xx.clone().min(yy.clone())
            };
            hi += xx.max(yy);
        }
        let (wlo, whi) = range(&self.span.homogeneous[3], a, b);
        let wlo = wlo.max(self.span.min_weight.clone());
        let whi = whi.min(self.span.max_weight.clone());
        debug_assert!(zero() < wlo && wlo <= whi);
        (lo / (&whi * &whi), hi / (&wlo * &wlo))
    }
}
struct Context {
    roots: Budget,
    images: ImageBudget,
    candidates_left: usize,
    image_cache: Vec<(Arc<Span>, Vec<AlgebraicRoot>)>,
}
impl Context {
    fn new(options: SplineProximityOptions) -> Self {
        Self {
            roots: Budget::new(options.root_isolation),
            images: ImageBudget::new(
                options.max_image_coefficient_updates,
                options.max_image_coefficient_bits,
            ),
            candidates_left: options.max_candidates,
            image_cache: vec![],
        }
    }
    fn consider(&mut self, winners: &mut Vec<Candidate>, mut candidate: Candidate) -> Result<()> {
        self.candidates_left = self
            .candidates_left
            .checked_sub(1)
            .ok_or(Error::ComputationLimit("spline distance candidates"))?;
        let order = match winners.first_mut() {
            Some(best) => self.compare(&mut candidate.at, &mut best.at)?,
            None => Less,
        };
        if order == Less {
            winners.clear();
        }
        if order != Greater {
            winners.push(candidate);
        }
        Ok(())
    }
    fn compare(&mut self, a: &mut Distance, b: &mut Distance) -> Result<Ordering> {
        if let Some(value) = a.rational_value() {
            return Ok(b.compare_rational(&value).reverse());
        }
        if let Some(value) = b.rational_value() {
            return Ok(a.compare_rational(&value));
        }
        for steps in [0, 8, 16, 32, 64] {
            a.root.refine_for_signs(steps);
            b.root.refine_for_signs(steps);
            let (al, ah) = a.range();
            let (bl, bh) = b.range();
            if ah < bl {
                return Ok(Less);
            }
            if bh < al {
                return Ok(Greater);
            }
            if al == zero() && a.compare_rational(&zero()) == Equal {
                return Ok(b.compare_rational(&zero()).reverse());
            }
            if bl == zero() && b.compare_rational(&zero()) == Equal {
                return Ok(Greater);
            }
            if let Some(value) = a.rational_value() {
                return Ok(b.compare_rational(&value).reverse());
            }
            if let Some(value) = b.rational_value() {
                return Ok(a.compare_rational(&value));
            }
        }
        // Squared distance is nonnegative, so exact zero needs no image.
        if a.compare_rational(&zero()) == Equal {
            return Ok(b.compare_rational(&zero()).reverse());
        }
        if b.compare_rational(&zero()) == Equal {
            return Ok(Greater);
        }
        if a.root.compare_root(&b.root) == Equal {
            let difference = subtract(
                &product(&a.span.numerator, &b.span.denominator),
                &product(&b.span.numerator, &a.span.denominator),
                &one(),
            );
            return Ok(a
                .root
                .sign_polynomial(&IntPolynomial::from_rationals(&difference)));
        }
        // A rational squared distance can occur at an irrational parameter.
        // The continued-fraction candidate is only a proposal: exact polynomial
        // vanishing certifies it before it can decide a comparison.
        let (al, ah) = a.range();
        let (bl, bh) = b.range();
        let value = real::rational_in_interval(&al.max(bl), &ah.min(bh));
        for (selected, other, reverse) in [(&*a, &*b, true), (&*b, &*a, false)] {
            let difference = subtract(&selected.span.numerator, &selected.span.denominator, &value);
            if selected
                .root
                .vanishes_polynomial(&IntPolynomial::from_rationals(&difference))
            {
                let order = other.compare_rational(&value);
                return Ok(if reverse { order.reverse() } else { order });
            }
        }
        Ok(self.image(a)?.compare_root(&self.image(b)?))
    }
    fn image(&mut self, value: &Distance) -> Result<AlgebraicRoot> {
        let index = if let Some(i) = self
            .image_cache
            .iter()
            .position(|(s, _)| Arc::ptr_eq(s, &value.span))
        {
            i
        } else {
            let roots = image_roots(
                &value.root,
                &value.span.numerator,
                &value.span.denominator,
                &mut self.roots,
                &mut self.images,
            )?;
            self.image_cache.push((value.span.clone(), roots));
            self.image_cache.len() - 1
        };
        for candidate in &self.image_cache[index].1 {
            let (a, b) = candidate.isolator();
            if value.compare_rational(a) != Less && value.compare_rational(b) != Greater {
                return Ok(candidate.clone());
            }
        }
        unreachable!("annihilating polynomial contains the selected finite image")
    }
}

fn dot(a: &[R; 3], b: &[R; 3]) -> R {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn zero() -> R {
    integer(0)
}
fn one() -> R {
    integer(1)
}
fn integer(x: i32) -> R {
    R::from_integer(BigInt::from(x))
}
fn subtract(a: &[R], b: &[R], scale: &R) -> Vec<R> {
    let mut r = a.to_vec();
    r.resize(r.len().max(b.len()), zero());
    for (x, y) in r.iter_mut().zip(b) {
        *x -= y * scale;
    }
    trim(&mut r);
    r
}
fn trim(r: &mut Vec<R>) {
    while r.last().is_some_and(|c| c == &zero()) {
        r.pop();
    }
}
fn product(a: &[R], b: &[R]) -> Vec<R> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            r[i + j] += x * y;
        }
    }
    trim(&mut r);
    r
}
fn derivative(a: &[R]) -> Vec<R> {
    a.iter()
        .enumerate()
        .skip(1)
        .map(|(i, x)| x * R::from_integer(BigInt::from(i)))
        .collect()
}
fn evaluate(a: &[R], x: &R) -> R {
    a.iter().rev().fold(zero(), |v, c| v * x + c)
}
fn range(p: &[R], a: &R, b: &R) -> (R, R) {
    p.iter().rev().fold((zero(), zero()), |(lo, hi), c| {
        let products = [&lo * a, &lo * b, &hi * a, &hi * b];
        (
            products.iter().min().unwrap() + c,
            products.iter().max().unwrap() + c,
        )
    })
}
