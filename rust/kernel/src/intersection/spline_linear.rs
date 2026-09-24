//! Complete parameter preimages of lines and closed segments on rational curves.
//! Reference: OCCT IntTools_EdgeEdge and IntTools_CommonPrt. Exact common roots
//! and closed polynomial inequalities replace tolerance-based range searches.
use crate::polynomial::{
    real::{self, Budget, IntPolynomial},
    AlgebraicRoot, RootIsolationOptions,
};
use crate::{
    curve::{DerivativeOrder, KnotSide},
    interval,
    math::finite,
    spline, BSplineCurve3, Bounds3, Error, ExactBSplineCurve3, Point3, Result, ScalarInterval,
};
use num_rational::BigRational as R;
use std::{
    cmp::Ordering::{self, Equal, Greater, Less},
    sync::Arc,
};

#[cfg(test)]
mod tests;

/// Deterministic limits on span traversal, root subdivision and candidate
/// boundaries before consolidation. They are not geometric tolerances or hard
/// CPU/memory deadlines. Exhaustion returns no partial intersection set.
#[derive(Debug, Clone, Copy)]
pub struct SplineLinearOptions {
    pub root_isolation: RootIsolationOptions,
    pub max_spans: usize,
    pub max_candidates: usize,
}
impl Default for SplineLinearOptions {
    fn default() -> Self {
        Self {
            root_isolation: RootIsolationOptions::default(),
            max_spans: spline::MAX_POLES,
            max_candidates: 65_536,
        }
    }
}

/// Every isolated curve parameter and every maximal closed overlap interval.
/// Different parameters survive when their points coincide. Interval endpoints
/// are omitted from the isolated list. Floating bounds are optional views.
#[derive(Debug, Clone)]
pub struct SplineLinearIntersection {
    points: Vec<SplineLinearPoint>,
    overlaps: Vec<SplineLinearOverlap>,
}
impl SplineLinearIntersection {
    pub fn points(&self) -> &[SplineLinearPoint] {
        &self.points
    }
    pub fn overlaps(&self) -> &[SplineLinearOverlap] {
        &self.overlaps
    }
    pub fn is_disjoint(&self) -> bool {
        self.points.is_empty() && self.overlaps.is_empty()
    }
}

/// Maximal closed parameter interval, including constant-point curves and
/// backtracking. Endpoint positions alone need not bound its geometric image.
#[derive(Debug, Clone)]
pub struct SplineLinearOverlap {
    endpoints: [SplineLinearPoint; 2],
}
impl SplineLinearOverlap {
    pub fn endpoints(&self) -> &[SplineLinearPoint; 2] {
        &self.endpoints
    }
    pub fn parameter_bounds(&self) -> Result<[ScalarInterval; 2]> {
        Ok([
            self.endpoints[0].parameter_bounds()?,
            self.endpoints[1].parameter_bounds()?,
        ])
    }
}

/// Algebraic curve parameter with exact coordinate and linear-parameter views.
/// The line/segment parameter is t in A+t(B-A), in original endpoint order.
/// A collapsed segment uses t=0 as its deterministic representative.
#[derive(Debug, Clone)]
pub struct SplineLinearPoint {
    span: Arc<Span>,
    root: AlgebraicRoot,
}
impl SplineLinearPoint {
    pub fn compare_parameter(&self, value: &R) -> Result<Ordering> {
        Ok(self.root.compare_rational(
            &((spline::normalize(value)? - &self.span.start) / &self.span.length),
        ))
    }
    /// Exact ordering in original parameter units, including different spans,
    /// curves, queries and periodic turns. No enclosure decides equality.
    pub fn parameter_cmp(&self, other: &Self) -> Ordering {
        self.root.compare_affine(
            &self.span.start,
            &self.span.length,
            &other.root,
            &other.span.start,
            &other.span.length,
        )
    }
    /// A rational witness when recognized; None does not imply irrationality.
    pub fn rational_parameter(&self) -> Option<R> {
        self.root
            .rational_value()
            .map(|x| &self.span.start + &self.span.length * x)
    }
    pub fn parameter_bounds(&self) -> Result<ScalarInterval> {
        let point = self.for_view();
        interval::enclose(
            |x| {
                point
                    .compare_parameter(&real::rat(x))
                    .expect("finite enclosure probe")
            },
            "spline linear intersection parameter",
        )
    }
    pub fn compare_coordinate(&self, component: usize, value: &R) -> Result<Ordering> {
        if component >= 3 {
            return Err(Error::OutOfDomain("spline linear coordinate"));
        }
        let value = spline::normalize(value)?;
        Ok(self
            .root
            .sign_polynomial(&IntPolynomial::from_rationals(&add_scaled(
                &self.span.homogeneous[component],
                &self.span.homogeneous[3],
                &-value,
            ))))
    }
    pub fn coordinate_bounds(&self) -> Result<[ScalarInterval; 3]> {
        let point = self.for_view();
        let bound = |i| {
            interval::enclose(
                |x| {
                    point
                        .compare_coordinate(i, &real::rat(x))
                        .expect("validated component and finite enclosure probe")
                },
                "spline linear intersection coordinate",
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
    pub fn compare_linear_parameter(&self, value: &R) -> Result<Ordering> {
        let value = spline::normalize(value)?;
        Ok(self
            .root
            .sign_polynomial(&IntPolynomial::from_rationals(&add_scaled(
                &self.span.along,
                &self.span.linear_denominator,
                &-value,
            ))))
    }
    pub fn linear_parameter_bounds(&self) -> Result<ScalarInterval> {
        let point = self.for_view();
        interval::enclose(
            |x| {
                point
                    .compare_linear_parameter(&real::rat(x))
                    .expect("finite enclosure probe")
            },
            "spline intersection line parameter",
        )
    }
    fn for_view(&self) -> Self {
        let mut point = self.clone();
        point.root.refine_for_signs(128);
        point
    }
}

#[derive(Debug)]
struct Span {
    start: R,
    length: R,
    homogeneous: [Vec<R>; 4],
    along: Vec<R>,
    linear_denominator: Vec<R>,
}
impl Span {
    fn new(
        start: R,
        length: R,
        homogeneous: [Vec<R>; 4],
        anchor: &[R; 3],
        direction: &[R; 3],
        norm: &R,
    ) -> (Self, IntPolynomial) {
        let w = &homogeneous[3];
        let delta: [Vec<R>; 3] =
            std::array::from_fn(|i| add_scaled(&homogeneous[i], w, &-&anchor[i]));
        let (equations, along, linear_denominator) = if norm == &zero() {
            (delta, vec![], vec![one()])
        } else {
            let cross: [Vec<R>; 3] = std::array::from_fn(|i| {
                add_scaled(
                    &scaled(&delta[(i + 1) % 3], &direction[(i + 2) % 3]),
                    &delta[(i + 2) % 3],
                    &-&direction[(i + 1) % 3],
                )
            });
            let along = (0..3).fold(vec![], |p, i| add_scaled(&p, &delta[i], &direction[i]));
            (cross, along, scaled(w, norm))
        };
        let common = equations.iter().fold(IntPolynomial::new(vec![]), |g, p| {
            g.gcd(&IntPolynomial::from_rationals(p))
        });
        (
            Self {
                start,
                length,
                homogeneous,
                along,
                linear_denominator,
            },
            common,
        )
    }
}

struct Context {
    roots: Budget,
    candidates_left: usize,
}
impl Context {
    fn new(options: SplineLinearOptions) -> Self {
        Self {
            roots: Budget::new(options.root_isolation),
            candidates_left: options.max_candidates,
        }
    }
    fn charge(&mut self, n: usize) -> Result<()> {
        self.candidates_left = self
            .candidates_left
            .checked_sub(n)
            .ok_or(Error::ComputationLimit("spline linear candidates"))?;
        Ok(())
    }
}

struct ClosedClip {
    points: Vec<AlgebraicRoot>,
    intervals: Vec<[AlgebraicRoot; 2]>,
}
fn right_sign(root: &AlgebraicRoot, polynomial: &IntPolynomial) -> Ordering {
    let mut p = polynomial.clone();
    while !p.is_zero() {
        let sign = root.sign_polynomial(&p);
        if sign != Equal {
            return sign;
        }
        // The first nonzero Taylor derivative fixes the sign immediately to
        // the right. Factorials and positive displacement powers preserve it.
        p = p.derivative();
    }
    Equal
}
fn closed_clip(
    p: &IntPolynomial,
    q: &IntPolynomial,
    lo: R,
    hi: R,
    context: &mut Context,
) -> Result<ClosedClip> {
    let mut boundaries = vec![AlgebraicRoot::rational(lo.clone())];
    if lo != hi {
        boundaries.push(AlgebraicRoot::rational(hi.clone()));
    }
    context.charge(boundaries.len())?;
    for equation in [p, q] {
        if !equation.is_zero() {
            let roots = real::isolate(equation, lo.clone(), hi.clone(), &mut context.roots)?;
            context.charge(roots.len())?;
            boundaries.extend(roots);
        }
    }
    boundaries.sort_by(AlgebraicRoot::compare_root);
    boundaries.dedup_by(|a, b| a.compare_root(b) == Equal);
    let accepted: Vec<_> = boundaries
        .iter()
        .map(|r| r.sign_polynomial(p) != Less && r.sign_polynomial(q) != Less)
        .collect();
    let mut ranges: Vec<(usize, usize)> = vec![];
    for i in 0..boundaries.len() - 1 {
        if right_sign(&boundaries[i], p) != Less && right_sign(&boundaries[i], q) != Less {
            // Non-strict polynomial inequalities define closed sets.
            debug_assert!(accepted[i] && accepted[i + 1]);
            if let Some(last) = ranges.last_mut().filter(|v| v.1 == i) {
                last.1 = i + 1;
            } else {
                ranges.push((i, i + 1));
            }
        }
    }
    let points = (0..boundaries.len())
        .filter(|&i| accepted[i] && !ranges.iter().any(|&(a, b)| a <= i && i <= b))
        .map(|i| boundaries[i].clone())
        .collect();
    let intervals = ranges
        .into_iter()
        .map(|(a, b)| [boundaries[a].clone(), boundaries[b].clone()])
        .collect();
    Ok(ClosedClip { points, intervals })
}

fn intersect(
    curve: &ExactBSplineCurve3,
    endpoints: [&[R; 3]; 2],
    segment: bool,
    first: &R,
    last: &R,
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    let normalize_point = |p: &[R; 3]| -> Result<[R; 3]> {
        Ok([
            spline::normalize(&p[0])?,
            spline::normalize(&p[1])?,
            spline::normalize(&p[2])?,
        ])
    };
    let anchor = normalize_point(endpoints[0])?;
    let end = normalize_point(endpoints[1])?;
    let direction: [R; 3] = std::array::from_fn(|i| &end[i] - &anchor[i]);
    let norm: R = direction.iter().map(|x| x * x).sum();
    if !segment && norm == zero() {
        return Err(Error::Degenerate("line defining points"));
    }
    let (first, last) = (spline::normalize(first)?, spline::normalize(last)?);
    if first > last {
        return Err(Error::OutOfDomain("reversed spline linear range"));
    }
    let mut context = Context::new(options);
    let mut result = SplineLinearIntersection {
        points: vec![],
        overlaps: vec![],
    };
    let mut span = |start, length, h, lo: R, hi: R| -> Result<()> {
        let (span, common) = Span::new(start, length, h, &anchor, &direction, &norm);
        let span = Arc::new(span);
        let at = |mut root: AlgebraicRoot| {
            root.refine_for_signs(16);
            SplineLinearPoint {
                span: span.clone(),
                root,
            }
        };
        if common.is_zero() {
            let clip = if segment && norm != zero() {
                closed_clip(
                    &IntPolynomial::from_rationals(&span.along),
                    &IntPolynomial::from_rationals(&add_scaled(
                        &span.linear_denominator,
                        &span.along,
                        &-one(),
                    )),
                    lo,
                    hi,
                    &mut context,
                )?
            } else if lo == hi {
                context.charge(1)?;
                ClosedClip {
                    points: vec![AlgebraicRoot::rational(lo)],
                    intervals: vec![],
                }
            } else {
                context.charge(2)?;
                ClosedClip {
                    points: vec![],
                    intervals: vec![[AlgebraicRoot::rational(lo), AlgebraicRoot::rational(hi)]],
                }
            };
            result.points.extend(clip.points.into_iter().map(at));
            result
                .overlaps
                .extend(clip.intervals.into_iter().map(|v| SplineLinearOverlap {
                    endpoints: v.map(at),
                }));
        } else {
            let roots = real::isolate(&common, lo, hi, &mut context.roots)?;
            context.charge(roots.len())?;
            let lower = IntPolynomial::from_rationals(&span.along);
            let upper = IntPolynomial::from_rationals(&add_scaled(
                &span.linear_denominator,
                &span.along,
                &-one(),
            ));
            result.points.extend(
                roots
                    .into_iter()
                    .filter(|r| {
                        !segment
                            || (r.sign_polynomial(&lower) != Less
                                && r.sign_polynomial(&upper) != Less)
                    })
                    .map(at),
            );
        }
        Ok(())
    };
    if first == last {
        if options.max_spans == 0 {
            return Err(Error::ComputationLimit("spline span traversal"));
        }
        let value = curve
            .exact_evaluate(&first, DerivativeOrder::Position, KnotSide::Automatic)?
            .position()
            .coordinates()
            .clone();
        span(
            first,
            one(),
            [
                vec![value[0].clone()],
                vec![value[1].clone()],
                vec![value[2].clone()],
                vec![one()],
            ],
            zero(),
            zero(),
        )?;
    } else {
        for piece in curve
            .knot_vector()
            .spans_in(&first, &last, options.max_spans)?
        {
            let length = &piece.end - &piece.start;
            let lo = (&piece.lower - &piece.start) / &length;
            let hi = (&piece.upper - &piece.start) / &length;
            span(
                piece.start,
                length,
                curve.span_polynomial(piece.index),
                lo,
                hi,
            )?;
        }
    }
    // Spans and each local result are already ordered in original parameter
    // units. Consolidate shared closed endpoints without merging point aliases.
    let mut overlaps: Vec<SplineLinearOverlap> = vec![];
    for interval in result.overlaps {
        if let Some(last) = overlaps
            .last_mut()
            .filter(|v| v.endpoints[1].parameter_cmp(&interval.endpoints[0]) != Less)
        {
            if last.endpoints[1].parameter_cmp(&interval.endpoints[1]) == Less {
                last.endpoints[1] = interval.endpoints[1].clone();
            }
        } else {
            overlaps.push(interval);
        }
    }
    result.points.dedup_by(|a, b| a.parameter_cmp(b) == Equal);
    let mut index = 0;
    result.points.retain(|point| {
        while index < overlaps.len() && overlaps[index].endpoints[1].parameter_cmp(point) == Less {
            index += 1;
        }
        index == overlaps.len() || overlaps[index].endpoints[0].parameter_cmp(point) == Greater
    });
    result.overlaps = overlaps;
    Ok(result)
}

fn binary_inputs(a: Point3, b: Point3, first: f64, last: f64) -> Result<[[R; 3]; 2]> {
    for p in [a, b] {
        for x in p.to_array() {
            finite(x, "spline intersection line point")?;
        }
    }
    finite(first, "spline intersection range")?;
    finite(last, "spline intersection range")?;
    Ok([a.to_array().map(real::rat), b.to_array().map(real::rat)])
}
fn add_scaled(a: &[R], b: &[R], scale: &R) -> Vec<R> {
    (0..a.len().max(b.len()))
        .map(|i| {
            a.get(i).cloned().unwrap_or_else(zero) + scale * b.get(i).cloned().unwrap_or_else(zero)
        })
        .collect()
}
fn scaled(p: &[R], scale: &R) -> Vec<R> {
    p.iter().map(|x| x * scale).collect()
}
fn zero() -> R {
    R::from_integer(0.into())
}
fn one() -> R {
    R::from_integer(1.into())
}

/// Complete exact preimage of the infinite line through A and B on the fundamental domain.
/// Equal endpoints return a degenerate-line error.
pub fn spline_line(
    curve: &BSplineCurve3,
    a: Point3,
    b: Point3,
) -> Result<SplineLinearIntersection> {
    spline_line_with_options(curve, a, b, SplineLinearOptions::default())
}
/// Fundamental domain with explicit span, root and candidate work limits.
pub fn spline_line_with_options(
    curve: &BSplineCurve3,
    a: Point3,
    b: Point3,
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    let (first, last) = curve.domain();
    spline_line_in_with_options(curve, a, b, first, last, options)
}
/// Complete closed range, including singletons and multiple periodic turns.
/// Nonperiodic queries must remain inside the active curve domain.
pub fn spline_line_in(
    curve: &BSplineCurve3,
    a: Point3,
    b: Point3,
    first: f64,
    last: f64,
) -> Result<SplineLinearIntersection> {
    spline_line_in_with_options(curve, a, b, first, last, SplineLinearOptions::default())
}
/// Explicit work limits; failure returns no partial result. Exact identities
/// remain usable when their floating coordinate or parameter views fail.
pub fn spline_line_in_with_options(
    curve: &BSplineCurve3,
    a: Point3,
    b: Point3,
    first: f64,
    last: f64,
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    let [a, b] = binary_inputs(a, b, first, last)?;
    intersect(
        &curve.to_exact(),
        [&a, &b],
        false,
        &real::rat(first),
        &real::rat(last),
        options,
    )
}

/// Complete exact preimage of the infinite line through A and B on the fundamental domain.
/// Equal endpoints return a degenerate-line error.
pub fn exact_spline_line(
    curve: &ExactBSplineCurve3,
    a: &[R; 3],
    b: &[R; 3],
) -> Result<SplineLinearIntersection> {
    exact_spline_line_with_options(curve, a, b, SplineLinearOptions::default())
}
/// Fundamental domain with explicit span, root and candidate work limits.
pub fn exact_spline_line_with_options(
    curve: &ExactBSplineCurve3,
    a: &[R; 3],
    b: &[R; 3],
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    let [first, last] = curve.domain();
    exact_spline_line_in_with_options(curve, a, b, first, last, options)
}
/// Complete closed range, including singletons and multiple periodic turns.
/// Nonperiodic queries must remain inside the active curve domain.
pub fn exact_spline_line_in(
    curve: &ExactBSplineCurve3,
    a: &[R; 3],
    b: &[R; 3],
    first: &R,
    last: &R,
) -> Result<SplineLinearIntersection> {
    exact_spline_line_in_with_options(curve, a, b, first, last, SplineLinearOptions::default())
}
/// Explicit work limits; failure returns no partial result. Exact identities
/// remain usable when their floating coordinate or parameter views fail.
pub fn exact_spline_line_in_with_options(
    curve: &ExactBSplineCurve3,
    a: &[R; 3],
    b: &[R; 3],
    first: &R,
    last: &R,
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    intersect(curve, [a, b], false, first, last, options)
}

/// Complete exact preimage of the closed segment from A to B on the fundamental domain.
/// Equal endpoints are treated as a point with linear parameter zero.
pub fn spline_segment(
    curve: &BSplineCurve3,
    a: Point3,
    b: Point3,
) -> Result<SplineLinearIntersection> {
    spline_segment_with_options(curve, a, b, SplineLinearOptions::default())
}
/// Fundamental domain with explicit span, root and candidate work limits.
pub fn spline_segment_with_options(
    curve: &BSplineCurve3,
    a: Point3,
    b: Point3,
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    let (first, last) = curve.domain();
    spline_segment_in_with_options(curve, a, b, first, last, options)
}
/// Complete closed range, including singletons and multiple periodic turns.
/// Nonperiodic queries must remain inside the active curve domain.
pub fn spline_segment_in(
    curve: &BSplineCurve3,
    a: Point3,
    b: Point3,
    first: f64,
    last: f64,
) -> Result<SplineLinearIntersection> {
    spline_segment_in_with_options(curve, a, b, first, last, SplineLinearOptions::default())
}
/// Explicit work limits; failure returns no partial result. Exact identities
/// remain usable when their floating coordinate or parameter views fail.
pub fn spline_segment_in_with_options(
    curve: &BSplineCurve3,
    a: Point3,
    b: Point3,
    first: f64,
    last: f64,
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    let [a, b] = binary_inputs(a, b, first, last)?;
    intersect(
        &curve.to_exact(),
        [&a, &b],
        true,
        &real::rat(first),
        &real::rat(last),
        options,
    )
}

/// Complete exact preimage of the closed segment from A to B on the fundamental domain.
/// Equal endpoints are treated as a point with linear parameter zero.
pub fn exact_spline_segment(
    curve: &ExactBSplineCurve3,
    a: &[R; 3],
    b: &[R; 3],
) -> Result<SplineLinearIntersection> {
    exact_spline_segment_with_options(curve, a, b, SplineLinearOptions::default())
}
/// Fundamental domain with explicit span, root and candidate work limits.
pub fn exact_spline_segment_with_options(
    curve: &ExactBSplineCurve3,
    a: &[R; 3],
    b: &[R; 3],
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    let [first, last] = curve.domain();
    exact_spline_segment_in_with_options(curve, a, b, first, last, options)
}
/// Complete closed range, including singletons and multiple periodic turns.
/// Nonperiodic queries must remain inside the active curve domain.
pub fn exact_spline_segment_in(
    curve: &ExactBSplineCurve3,
    a: &[R; 3],
    b: &[R; 3],
    first: &R,
    last: &R,
) -> Result<SplineLinearIntersection> {
    exact_spline_segment_in_with_options(curve, a, b, first, last, SplineLinearOptions::default())
}
/// Explicit work limits; failure returns no partial result. Exact identities
/// remain usable when their floating coordinate or parameter views fail.
pub fn exact_spline_segment_in_with_options(
    curve: &ExactBSplineCurve3,
    a: &[R; 3],
    b: &[R; 3],
    first: &R,
    last: &R,
    options: SplineLinearOptions,
) -> Result<SplineLinearIntersection> {
    intersect(curve, [a, b], true, first, last, options)
}
