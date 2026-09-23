//! Certified minimum distances for linear sets and point-to-spline queries.
//!
//! All 25 ordered pairings of points, infinite lines, closed segments, infinite
//! planes and closed triangles are supported. Finite binary64 inputs define
//! exact rational geometry; no snapping, angular threshold or modeling tolerance
//! is applied. A collapsed segment is a point; a collapsed line is invalid.
//!
//! A query returns one deterministic minimum witness, including when infinitely
//! many pairs minimize distance. It does not enumerate minima or assert uniqueness.
//! Exact rational points and original affine parameters remain available even
//! when their finite binary64 enclosures cannot be represented. Rounded points
//! need not lie exactly on either operand: use the rational values or bounds.
//!
//! The separate `closest_points_on_spline*` APIs enumerate every globally
//! closest parameter on a closed rational B-spline range, including whole
//! minimizing intervals and repeated periodic aliases. They preserve algebraic
//! identities, expose exact comparisons, and report exhausted work limits.
//!
//! OCCT references: Extrema_ExtElC (line/line normal equations), Extrema_ExtPElC
//! and Extrema_ExtPElS (projection), and BRepExtrema_DistShapeShape (boundary
//! strata). Exact rank and complete face enumeration replace native tolerance
//! decisions. See rust/SOURCE_MAP.md and rust/MATHEMATICS.md.
//!
//! ```
//! use rusty_occt::proximity::{closest_points, LinearPrimitive3 as Shape};
//! use rusty_occt::Point3;
//! use std::cmp::Ordering;
//!
//! let point = Shape::Point(Point3::new(1., 2., 0.));
//! let segment = Shape::Segment([Point3::ORIGIN, Point3::new(4., 0., 0.)]);
//! let pair = closest_points(&point, &segment)?;
//! assert_eq!(pair.compare_distance(2.)?, Ordering::Equal);
//! assert_eq!(pair.compare_parameter(1, 0, 0.25)?, Ordering::Equal);
//! assert_eq!(pair.point(1)?, Point3::new(1., 0., 0.));
//! // Retain pair.exact_points() for further exact incidence calculations.
//! # Ok::<(), rusty_occt::Error>(())
//! ```

use crate::{interval, Bounds3, Error, Plane3, Point3, Result, ScalarInterval, Triangle3};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::cmp::Ordering;

mod spline;
pub use spline::{
    closest_points_on_exact_spline, closest_points_on_exact_spline_in, closest_points_on_spline,
    closest_points_on_spline_in, SplineClosestInterval, SplineClosestPoint, SplineClosestSet,
    SplineProximityOptions,
};

type Vector = [R; 3];
type IntegerVector = [BigInt; 3];

/// The geometric sets accepted by [`closest_points`] and
/// [`crate::intersection::linear_intersection`]. Line and segment
/// parameters use `a + t(b-a)`; plane and triangle parameters use
/// `a + u(b-a) + v(c-a)` in the original defining-point order.
#[derive(Debug, Clone)]
pub enum LinearPrimitive3 {
    Point(Point3),
    Line([Point3; 2]),
    Segment([Point3; 2]),
    Plane(Plane3),
    Triangle(Triangle3),
}

/// One exact minimum pair. Rational coordinates are constructions, not rounded
/// samples; distance is the nonnegative square root of `exact_squared_distance`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClosestPair {
    points: [Vector; 2],
    parameters: [Vec<R>; 2],
    squared_distance: R,
}

impl ClosestPair {
    pub fn exact_points(&self) -> &[Vector; 2] {
        &self.points
    }

    /// Empty for points; one parameter for lines/segments; two for planes/triangles.
    /// A collapsed segment still has one valid parameter, selected deterministically.
    pub fn exact_parameters(&self) -> &[Vec<R>; 2] {
        &self.parameters
    }

    pub fn exact_squared_distance(&self) -> &R {
        &self.squared_distance
    }

    /// Compare two distances exactly, without square roots or conversion to f64.
    pub fn distance_cmp(&self, other: &Self) -> Ordering {
        self.squared_distance.cmp(&other.squared_distance)
    }

    /// Exact ordering of the nonnegative distance and a finite threshold.
    pub fn compare_distance(&self, threshold: f64) -> Result<Ordering> {
        let r = rational(threshold)?;
        Ok(if r < zero() {
            Ordering::Greater
        } else {
            self.squared_distance.cmp(&(&r * &r))
        })
    }

    pub fn compare_squared_distance(&self, threshold: f64) -> Result<Ordering> {
        Ok(self.squared_distance.cmp(&rational(threshold)?))
    }

    /// Minimal finite binary64 enclosure of the distance, including subnormals.
    pub fn distance_bounds(&self) -> Result<ScalarInterval> {
        interval::enclose(
            |x| self.compare_distance(x).expect("finite enclosure probe"),
            "minimum distance",
        )
    }

    /// May be unrepresentable even when the unsquared distance is representable.
    pub fn squared_distance_bounds(&self) -> Result<ScalarInterval> {
        bounds(&self.squared_distance, "squared minimum distance")
    }

    pub fn point_bounds(&self, operand: usize) -> Result<Bounds3> {
        let p = self
            .points
            .get(operand)
            .ok_or(Error::OutOfDomain("closest-pair operand"))?;
        let [x, y, z] = [
            bounds(&p[0], "closest point")?,
            bounds(&p[1], "closest point")?,
            bounds(&p[2], "closest point")?,
        ];
        Ok(Bounds3 {
            min: Point3::new(x.lower(), y.lower(), z.lower()),
            max: Point3::new(x.upper(), y.upper(), z.upper()),
        })
    }

    /// A rounded representative; use `exact_points` for exact incidence.
    pub fn point(&self, operand: usize) -> Result<Point3> {
        let b = self.point_bounds(operand)?;
        let p: [f64; 3] = std::array::from_fn(|i| {
            let lo = b.min.to_array()[i];
            lo + (b.max.to_array()[i] - lo) * 0.5
        });
        Ok(Point3::new(p[0], p[1], p[2]))
    }

    pub fn parameter_bounds(&self, operand: usize) -> Result<Vec<ScalarInterval>> {
        self.parameters
            .get(operand)
            .ok_or(Error::OutOfDomain("closest-pair operand"))?
            .iter()
            .map(|r| bounds(r, "closest-point parameter"))
            .collect()
    }

    pub fn compare_coordinate(&self, operand: usize, component: usize, x: f64) -> Result<Ordering> {
        let r = self
            .points
            .get(operand)
            .and_then(|p| p.get(component))
            .ok_or(Error::OutOfDomain("closest-point coordinate"))?;
        Ok(r.cmp(&rational(x)?))
    }

    pub fn compare_parameter(&self, operand: usize, parameter: usize, x: f64) -> Result<Ordering> {
        let r = self
            .parameters
            .get(operand)
            .and_then(|p| p.get(parameter))
            .ok_or(Error::OutOfDomain("closest-point parameter"))?;
        Ok(r.cmp(&rational(x)?))
    }
}

/// Minimize Euclidean distance over the two represented geometric sets.
///
/// There are at most 49 face pairs and four variables per exact linear system.
/// Singular systems retain a particular solution; lower-dimensional faces cover
/// constrained minima missed by that particular solution. No iterative numerical
/// convergence criterion is used. Integer arithmetic cost depends on input bits.
pub fn closest_points(a: &LinearPrimitive3, b: &LinearPrimitive3) -> Result<ClosestPair> {
    let [a, b] = [Shape::new(a)?, Shape::new(b)?];
    // Input rationals are dyadic, so their largest denominator is a common
    // denominator. Clear it once: the positive scale cancels from the normal
    // equations. Keep integer coefficients through fraction-free elimination.
    let denominator = a
        .points
        .iter()
        .chain(&b.points)
        .flat_map(|p| p.iter())
        .map(R::denom)
        .max()
        .expect("nonempty shapes");
    let [af, bf] = [a.faces(denominator), b.faces(denominator)];
    let mut best: Option<Candidate> = None;
    for fa in &af {
        for fb in &bf {
            let (solution, scale) = stationary(fa, fb);
            let (u, v) = solution.split_at(fa.directions.len());
            if !fa.feasible(u, &scale) || !fb.feasible(v, &scale) {
                continue;
            }
            let parameters = [
                fa.parameters(u, &scale, a.points.len()),
                fb.parameters(v, &scale, b.points.len()),
            ];
            let points = [fa.evaluate(u, &scale), fb.evaluate(v, &scale)];
            let delta = integer_sub(&points[0], &points[1]);
            let candidate = Candidate {
                points,
                parameters,
                squared_distance: integer_dot(&delta, &delta),
                denominator: scale,
            };
            // A stable tie choice amongst these candidates, not a global
            // lexicographic minimum of an unbounded family of minimizers.
            if best
                .as_ref()
                .is_none_or(|prior| candidate.cmp(prior) == Ordering::Less)
            {
                best = Some(candidate);
            }
        }
    }
    // Every nonempty closed polyhedral pair in this domain has a minimum, and
    // the face enumeration contains a feasible stationary representative.
    Ok(best
        .expect("complete face enumeration contains a closest pair")
        .into_result(denominator))
}

/// Homogeneous candidate: points are divided by denominator * input scale,
/// parameters by denominator, and squared distance by their squared product.
/// Postpone fraction reduction until after selecting the minimum witness.
struct Candidate {
    points: [IntegerVector; 2],
    parameters: [Vec<BigInt>; 2],
    squared_distance: BigInt,
    denominator: BigInt,
}
impl Candidate {
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.squared_distance * &other.denominator * &other.denominator)
            .cmp(&(&other.squared_distance * &self.denominator * &self.denominator))
            .then_with(|| {
                // Input scale is common and positive. Parameter counts are
                // fixed by the two original primitives, not the active faces.
                self.points
                    .iter()
                    .flatten()
                    .chain(self.parameters.iter().flatten())
                    .zip(
                        other
                            .points
                            .iter()
                            .flatten()
                            .chain(other.parameters.iter().flatten()),
                    )
                    .map(|(a, b)| (a * &other.denominator).cmp(&(b * &self.denominator)))
                    .find(|c| *c != Ordering::Equal)
                    .unwrap_or(Ordering::Equal)
            })
    }

    fn into_result(self, input_scale: &BigInt) -> ClosestPair {
        let coordinate_denominator = &self.denominator * input_scale;
        ClosestPair {
            points: self
                .points
                .map(|p| p.map(|x| R::new(x, coordinate_denominator.clone()))),
            parameters: self.parameters.map(|p| {
                p.into_iter()
                    .map(|x| R::new(x, self.denominator.clone()))
                    .collect()
            }),
            squared_distance: R::new(
                self.squared_distance,
                &coordinate_denominator * &coordinate_denominator,
            ),
        }
    }
}

fn zero() -> R {
    R::from_integer(BigInt::from(0))
}
fn rational(x: f64) -> Result<R> {
    R::from_float(x).ok_or(Error::NonFinite("proximity coordinate or threshold"))
}
fn bounds(r: &R, what: &'static str) -> Result<ScalarInterval> {
    interval::enclose(
        |x| r.cmp(&rational(x).expect("finite enclosure probe")),
        what,
    )
}
fn integer_sub(a: &IntegerVector, b: &IntegerVector) -> IntegerVector {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn integer_dot(a: &IntegerVector, b: &IntegerVector) -> BigInt {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

struct Shape {
    points: Vec<Vector>,
    bounded: bool,
}
impl Shape {
    fn new(primitive: &LinearPrimitive3) -> Result<Self> {
        use LinearPrimitive3::*;
        let (points, bounded) = match primitive {
            Point(p) => (vec![*p], true),
            Line(p) => (p.to_vec(), false),
            Segment(p) => (p.to_vec(), true),
            Plane(p) => (p.defining_points().to_vec(), false),
            Triangle(t) => (t.vertices().to_vec(), true),
        };
        let points = points
            .into_iter()
            .map(|p| Ok([rational(p.x)?, rational(p.y)?, rational(p.z)?]))
            .collect::<Result<Vec<_>>>()?;
        if matches!(primitive, Line(_)) && points[0] == points[1] {
            return Err(Error::Degenerate("line defining points"));
        }
        Ok(Self { points, bounded })
    }

    fn faces(&self, denominator: &BigInt) -> Vec<Face> {
        let points: Vec<IntegerVector> = self
            .points
            .iter()
            .map(|p| std::array::from_fn(|i| p[i].numer() * (denominator / p[i].denom())))
            .collect();
        let full = (1 << self.points.len()) - 1;
        let first = if self.bounded { 1 } else { full };
        (first..=full)
            .map(|mask| {
                let indices: Vec<_> = (0..self.points.len())
                    .filter(|&i| mask & (1 << i) != 0)
                    .collect();
                let anchor = points[indices[0]].clone();
                let directions = indices[1..]
                    .iter()
                    .map(|&i| integer_sub(&points[i], &anchor))
                    .collect();
                Face {
                    anchor,
                    directions,
                    indices,
                    bounded: self.bounded,
                }
            })
            .collect()
    }
}

struct Face {
    anchor: IntegerVector,
    directions: Vec<IntegerVector>,
    indices: Vec<usize>,
    bounded: bool,
}
impl Face {
    fn feasible(&self, x: &[BigInt], denominator: &BigInt) -> bool {
        !self.bounded
            || (x.iter().all(|v| *v >= BigInt::from(0)) && x.iter().sum::<BigInt>() <= *denominator)
    }

    fn parameters(&self, x: &[BigInt], denominator: &BigInt, vertices: usize) -> Vec<BigInt> {
        let mut weights = vec![BigInt::from(0); vertices];
        weights[self.indices[0]] = denominator - x.iter().sum::<BigInt>();
        for (&i, t) in self.indices[1..].iter().zip(x) {
            weights[i] = t.clone();
        }
        weights[1..].to_vec()
    }

    fn evaluate(&self, x: &[BigInt], denominator: &BigInt) -> IntegerVector {
        std::array::from_fn(|k| {
            &self.anchor[k] * denominator
                + self
                    .directions
                    .iter()
                    .zip(x)
                    .map(|(d, t)| &d[k] * t)
                    .sum::<BigInt>()
        })
    }
}

/// Solve DᵀD x = -Dᵀ(a-b) exactly. This system is always consistent, even
/// with dependent columns: ker(DᵀD)=ker(D), and the right side is in its range.
fn stationary(a: &Face, b: &Face) -> (Vec<BigInt>, BigInt) {
    let columns: Vec<_> = a
        .directions
        .iter()
        .cloned()
        .chain(b.directions.iter().map(|v| v.clone().map(|x| -x)))
        .collect();
    let n = columns.len();
    let delta = integer_sub(&a.anchor, &b.anchor);
    let mut rows: Vec<Vec<BigInt>> = columns
        .iter()
        .map(|c| {
            columns
                .iter()
                .map(|d| integer_dot(c, d))
                .chain([-integer_dot(c, &delta)])
                .collect()
        })
        .collect();
    let mut pivots = Vec::new();
    let mut previous = BigInt::from(1);
    let integer_zero = BigInt::from(0);
    for col in 0..n {
        let row = pivots.len();
        let Some(pivot) = (row..n).find(|&i| rows[i][col] != integer_zero) else {
            continue;
        };
        rows.swap(row, pivot);
        let pivot = rows[row].clone();
        // Bareiss elimination: each quotient is exact by Sylvester's identity.
        // Skipped zero columns do not change the previous nonzero pivot. Even
        // rows with a zero multiplier must receive the pivot/previous scaling.
        for cells in rows.iter_mut().skip(row + 1) {
            let scale = cells[col].clone();
            for (cell, value) in cells[col + 1..].iter_mut().zip(&pivot[col + 1..]) {
                let numerator = &*cell * &pivot[col] - &scale * value;
                debug_assert_eq!(&numerator % &previous, integer_zero);
                *cell = numerator / &previous;
            }
            cells[col] = integer_zero.clone();
        }
        previous = pivot[col].clone();
        pivots.push(col);
    }
    debug_assert!(rows[pivots.len()..]
        .iter()
        .all(|row| row[n] == integer_zero));
    let mut x = vec![integer_zero.clone(); n];
    let Some(&last) = pivots.last() else {
        return (x, BigInt::from(1));
    };
    // The final pivot is the determinant of a nonsingular subsystem. Cramer's
    // rule gives a common denominator for all non-free variables. Backsolve
    // its integer numerators; no per-operation rational reductions are needed.
    let mut denominator = rows[pivots.len() - 1][last].clone();
    for (row, col) in pivots.into_iter().enumerate().rev() {
        let rest: BigInt = rows[row][col + 1..n]
            .iter()
            .zip(&x[col + 1..])
            .map(|(a, x)| a * x)
            .sum();
        let numerator = &rows[row][n] * &denominator - rest;
        debug_assert_eq!(&numerator % &rows[row][col], integer_zero);
        x[col] = numerator / &rows[row][col];
    }
    if denominator < integer_zero {
        denominator = -denominator;
        for v in &mut x {
            *v = -std::mem::take(v);
        }
    }
    (x, denominator)
}
