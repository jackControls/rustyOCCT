//! Certified minimum distances between closed linear geometric sets.
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
//! OCCT references: Extrema_ExtElC (line/line normal equations), Extrema_ExtPElC
//! and Extrema_ExtPElS (projection), and BRepExtrema_DistShapeShape (boundary
//! strata). Exact rank and complete face enumeration replace native tolerance
//! decisions. See rust/SOURCE_MAP.md and rust/MATHEMATICS.md.

use crate::{interval, Bounds3, Error, Plane3, Point3, Result, ScalarInterval, Triangle3};
use num_bigint::BigInt;
use num_rational::BigRational as R;
use std::cmp::Ordering;

type Vector = [R; 3];

/// The geometric sets accepted by [`closest_points`]. Line and segment
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
    let [af, bf] = [a.faces(), b.faces()];
    let mut best: Option<ClosestPair> = None;
    for fa in &af {
        for fb in &bf {
            let solution = stationary(fa, fb);
            let (u, v) = solution.split_at(fa.directions.len());
            if !fa.feasible(u) || !fb.feasible(v) {
                continue;
            }
            let parameters = [
                fa.parameters(u, a.points.len()),
                fb.parameters(v, b.points.len()),
            ];
            let points = [a.evaluate(&parameters[0]), b.evaluate(&parameters[1])];
            let delta = sub(&points[0], &points[1]);
            let candidate = ClosestPair {
                points,
                parameters,
                squared_distance: dot(&delta, &delta),
            };
            // A stable tie choice amongst these candidates, not a global
            // lexicographic minimum of an unbounded family of minimizers.
            if best.as_ref().is_none_or(|prior| {
                (
                    &candidate.squared_distance,
                    &candidate.points,
                    &candidate.parameters,
                ) < (&prior.squared_distance, &prior.points, &prior.parameters)
            }) {
                best = Some(candidate);
            }
        }
    }
    // Every nonempty closed polyhedral pair in this domain has a minimum, and
    // the face enumeration contains a feasible stationary representative.
    Ok(best.expect("complete face enumeration contains a closest pair"))
}

fn zero() -> R {
    R::from_integer(BigInt::from(0))
}
fn one() -> R {
    R::from_integer(BigInt::from(1))
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
fn sub(a: &Vector, b: &Vector) -> Vector {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn dot(a: &Vector, b: &Vector) -> R {
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

    fn faces(&self) -> Vec<Face> {
        let full = (1 << self.points.len()) - 1;
        let first = if self.bounded { 1 } else { full };
        (first..=full)
            .map(|mask| {
                let indices: Vec<_> = (0..self.points.len())
                    .filter(|&i| mask & (1 << i) != 0)
                    .collect();
                let anchor = self.points[indices[0]].clone();
                let directions = indices[1..]
                    .iter()
                    .map(|&i| sub(&self.points[i], &anchor))
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

    fn evaluate(&self, parameters: &[R]) -> Vector {
        std::array::from_fn(|k| {
            &self.points[0][k]
                + parameters
                    .iter()
                    .enumerate()
                    .map(|(i, t)| t * (&self.points[i + 1][k] - &self.points[0][k]))
                    .sum::<R>()
        })
    }
}

struct Face {
    anchor: Vector,
    directions: Vec<Vector>,
    indices: Vec<usize>,
    bounded: bool,
}
impl Face {
    fn feasible(&self, x: &[R]) -> bool {
        !self.bounded || (x.iter().all(|v| *v >= zero()) && x.iter().sum::<R>() <= one())
    }

    fn parameters(&self, x: &[R], vertices: usize) -> Vec<R> {
        let mut weights = vec![zero(); vertices];
        weights[self.indices[0]] = one() - x.iter().sum::<R>();
        for (&i, t) in self.indices[1..].iter().zip(x) {
            weights[i] = t.clone();
        }
        weights[1..].to_vec()
    }
}

/// Solve DᵀD x = -Dᵀ(a-b) exactly. This system is always consistent, even
/// with dependent columns: ker(DᵀD)=ker(D), and the right side is in its range.
fn stationary(a: &Face, b: &Face) -> Vec<R> {
    let columns: Vec<_> = a
        .directions
        .iter()
        .cloned()
        .chain(b.directions.iter().map(|v| v.clone().map(|x| -x)))
        .collect();
    let n = columns.len();
    let delta = sub(&a.anchor, &b.anchor);
    let mut rows: Vec<Vec<R>> = columns
        .iter()
        .map(|c| {
            columns
                .iter()
                .map(|d| dot(c, d))
                .chain([-dot(c, &delta)])
                .collect()
        })
        .collect();
    let mut pivots = Vec::new();
    for col in 0..n {
        let row = pivots.len();
        let Some(pivot) = (row..n).find(|&i| rows[i][col] != zero()) else {
            continue;
        };
        rows.swap(row, pivot);
        let scale = rows[row][col].clone();
        for v in &mut rows[row][col..] {
            *v /= &scale;
        }
        let pivot = rows[row].clone();
        for (i, cells) in rows.iter_mut().enumerate() {
            if i == row {
                continue;
            }
            let scale = cells[col].clone();
            if scale == zero() {
                continue;
            }
            for (cell, value) in cells[col..].iter_mut().zip(&pivot[col..]) {
                *cell -= &scale * value;
            }
        }
        pivots.push(col);
    }
    debug_assert!(rows[pivots.len()..].iter().all(|row| row[n] == zero()));
    let mut x = vec![zero(); n];
    for (row, col) in pivots.into_iter().enumerate() {
        x[col] = rows[row][n].clone();
    }
    x
}
