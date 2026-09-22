//! Complete exact intersections of linear convex sets. See
//! rust/LINEAR_INTERSECTIONS.md for the API and native comparison contracts.
use crate::proximity::LinearPrimitive3 as Primitive;
use crate::{interval, Bounds3, Error, Point3, Result, ScalarInterval};
use num_bigint::BigInt;
use num_rational::BigRational as R;

type Vector = [R; 3];

/// An exact rational construction, including coordinates outside binary64's
/// finite range. Rounded representatives need not preserve incidence.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ExactPoint3(Vector);

impl ExactPoint3 {
    pub fn from_coordinates(coordinates: [R; 3]) -> Self {
        Self(coordinates)
    }
    pub fn coordinates(&self) -> &[R; 3] {
        &self.0
    }
    /// Minimal finite enclosure of one coordinate. Does not modify the exact point.
    pub fn coordinate_bounds(&self, axis: usize) -> Result<ScalarInterval> {
        let r = self.0.get(axis).ok_or(Error::OutOfDomain("point axis"))?;
        interval::enclose(
            |x| r.cmp(&R::from_float(x).expect("finite enclosure probe")),
            "intersection coordinate",
        )
    }
    pub fn bounds(&self) -> Result<Bounds3> {
        let [x, y, z] = [
            self.coordinate_bounds(0)?,
            self.coordinate_bounds(1)?,
            self.coordinate_bounds(2)?,
        ];
        Ok(Bounds3 {
            min: Point3::new(x.lower(), y.lower(), z.lower()),
            max: Point3::new(x.upper(), y.upper(), z.upper()),
        })
    }
    pub fn position(&self) -> Result<Point3> {
        Ok(Point3::new(
            self.coordinate_bounds(0)?.representative(),
            self.coordinate_bounds(1)?.representative(),
            self.coordinate_bounds(2)?.representative(),
        ))
    }
}

/// The entire closed-set intersection, in a deterministic canonical form.
/// Segment endpoints are lexicographically ordered. Polygon vertices are
/// cyclic, without duplicate/collinear vertices; the lexicographically smallest
/// vertex comes first and the smaller of the two cycle directions is selected.
/// These orientation-independent cycles do not carry a B-rep face orientation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LinearIntersection {
    Empty,
    Point(ExactPoint3),
    Segment([ExactPoint3; 2]),
    /// Three to six vertices, including the filled convex interior and boundary.
    Polygon(Vec<ExactPoint3>),
    /// origin + t*direction for every rational/real t. The first nonzero
    /// direction coordinate is one, and the corresponding origin coordinate zero.
    Line {
        origin: ExactPoint3,
        direction: [R; 3],
    },
    /// normal dot x = offset. The first nonzero normal coordinate is one.
    Plane {
        normal: [R; 3],
        offset: R,
    },
}

impl LinearIntersection {
    /// Affine dimension; the empty set has no dimension.
    pub fn dimension(&self) -> Option<usize> {
        match self {
            Self::Empty => None,
            Self::Point(_) => Some(0),
            Self::Segment(_) | Self::Line { .. } => Some(1),
            Self::Polygon(_) | Self::Plane { .. } => Some(2),
        }
    }
    /// Extreme points of a bounded result. Empty for empty/unbounded results.
    pub fn vertices(&self) -> &[ExactPoint3] {
        match self {
            Self::Point(p) => std::slice::from_ref(p),
            Self::Segment(p) => p,
            Self::Polygon(p) => p,
            _ => &[],
        }
    }
}

/// Complete intersection for all point/line/segment/plane/triangle pairings.
/// Finite input floats are exact values: no snapping or modeling tolerance is
/// applied. A collapsed segment is a point; a collapsed line is an error.
/// The exact result survives even when finite output bounds are unrepresentable.
/// Operand order and defining-point permutations do not change the result.
///
/// Source reference: OCCT IntAna_QuadQuadGeo (plane incidence), IntTools_EdgeEdge
/// (line incidence/overlap), and BOPAlgo_BOP/Section (dimension semantics).
/// The exact affine/halfspace algorithm below is independently implemented.
///
/// ```
/// use rusty_occt::intersection::{linear_intersection, LinearIntersection,
///     LinearPrimitive3 as Shape};
/// use rusty_occt::Point3;
/// let a = Shape::Segment([Point3::ORIGIN, Point3::new(4., 0., 0.)]);
/// let b = Shape::Segment([Point3::new(3., 0., 0.), Point3::new(1., 0., 0.)]);
/// let LinearIntersection::Segment(points) = linear_intersection(&a, &b)? else {
///     panic!("expected an overlap");
/// };
/// assert_eq!(points[0].position()?, Point3::new(1., 0., 0.));
/// assert_eq!(points[1].position()?, Point3::new(3., 0., 0.));
/// # Ok::<(), rusty_occt::Error>(())
/// ```
pub fn linear_intersection(a: &Primitive, b: &Primitive) -> Result<LinearIntersection> {
    let mut equations = Vec::new();
    let mut inequalities = Vec::new();
    constraints(a, &mut equations, &mut inequalities)?;
    constraints(b, &mut equations, &mut inequalities)?;
    let Some((origin, basis)) = affine_space(equations) else {
        return Ok(LinearIntersection::Empty);
    };
    // Every input has affine dimension <= 2. Each inequality is n.x >= c.
    let reduced: Vec<_> = inequalities
        .into_iter()
        .map(|(n, c)| {
            let coefficients: Vec<_> = basis.iter().map(|v| dot(&n, v)).collect();
            (coefficients, c - dot(&n, &origin))
        })
        .collect();
    Ok(match basis.len() {
        0 => {
            if reduced.iter().any(|(_, c)| *c > zero()) {
                LinearIntersection::Empty
            } else {
                LinearIntersection::Point(ExactPoint3(origin))
            }
        }
        1 => clip_line(origin, &basis[0], &reduced),
        2 if reduced.is_empty() => {
            let n = cross(&basis[0], &basis[1]);
            let pivot = n.iter().find(|x| **x != zero()).expect("plane rank");
            LinearIntersection::Plane {
                normal: n.clone().map(|x| x / pivot),
                offset: dot(&n, &origin) / pivot,
            }
        }
        2 => {
            // A bounded nonempty 2D polyhedron, including a lower-dimensional
            // contact, is the convex hull of its feasible boundary intersections.
            // At most six halfspaces and fifteen pairs for two triangles.
            let mut points = Vec::new();
            for i in 0..reduced.len() {
                let (u, r) = &reduced[i];
                for (v, s) in &reduced[..i] {
                    let det = &u[0] * &v[1] - &u[1] * &v[0];
                    if det == zero() {
                        continue;
                    }
                    let x = (r * &v[1] - &u[1] * s) / &det;
                    let y = (&u[0] * s - r * &v[0]) / det;
                    if reduced.iter().all(|(w, c)| &w[0] * &x + &w[1] * &y >= *c) {
                        points.push(ExactPoint3(std::array::from_fn(|j| {
                            &origin[j] + &basis[0][j] * &x + &basis[1][j] * &y
                        })));
                    }
                }
            }
            hull(points)
        }
        _ => unreachable!("input affine dimension <= 2"),
    })
}

type Constraint = (Vector, R);
fn constraints(
    primitive: &Primitive,
    equations: &mut Vec<Constraint>,
    inequalities: &mut Vec<Constraint>,
) -> Result<()> {
    let input = match primitive {
        Primitive::Point(p) => vec![*p],
        Primitive::Line(p) | Primitive::Segment(p) => p.to_vec(),
        Primitive::Plane(p) => p.defining_points().to_vec(),
        Primitive::Triangle(p) => p.vertices().to_vec(),
    };
    let p: Vec<Vector> = input
        .into_iter()
        .map(|p| {
            let [x, y, z] = p.to_array().map(|x| {
                R::from_float(x).ok_or(Error::NonFinite("linear intersection coordinate"))
            });
            Ok([x?, y?, z?])
        })
        .collect::<Result<_>>()?;
    if p.len() == 1 || (matches!(primitive, Primitive::Segment(_)) && p[0] == p[1]) {
        for (j, coordinate) in p[0].iter().enumerate() {
            equations.push((unit(j), coordinate.clone()));
        }
        return Ok(());
    }
    let d = sub(&p[1], &p[0]);
    if p.len() == 2 {
        let k = d
            .iter()
            .position(|x| *x != zero())
            .ok_or(Error::Degenerate("line defining points"))?;
        for j in (0..3).filter(|j| *j != k) {
            let mut n = [zero(), zero(), zero()];
            n[j] = d[k].clone();
            n[k] = -d[j].clone();
            let offset = dot(&n, &p[0]);
            equations.push((n, offset));
        }
        if matches!(primitive, Primitive::Segment(_)) {
            let (min, max) = if p[0][k] < p[1][k] {
                (&p[0], &p[1])
            } else {
                (&p[1], &p[0])
            };
            inequalities.push((unit(k), min[k].clone()));
            inequalities.push((unit(k).map(|x| -x), -max[k].clone()));
        }
    } else {
        let normal = cross(&d, &sub(&p[2], &p[0]));
        equations.push((normal.clone(), dot(&normal, &p[0])));
        if matches!(primitive, Primitive::Triangle(_)) {
            for j in 0..3 {
                let inward = cross(&normal, &sub(&p[(j + 1) % 3], &p[j]));
                let offset = dot(&inward, &p[j]);
                inequalities.push((inward, offset));
            }
        }
    }
    Ok(())
}

fn affine_space(equations: Vec<Constraint>) -> Option<(Vector, Vec<Vector>)> {
    let mut rows: Vec<Vec<R>> = equations
        .into_iter()
        .map(|(n, c)| n.into_iter().chain([c]).collect())
        .collect();
    let mut pivots = Vec::new();
    for col in 0..3 {
        let rank = pivots.len();
        let Some(row) = (rank..rows.len()).find(|i| rows[*i][col] != zero()) else {
            continue;
        };
        rows.swap(rank, row);
        let divisor = rows[rank][col].clone();
        for x in &mut rows[rank] {
            *x /= &divisor;
        }
        let pivot = rows[rank].clone();
        for (i, row) in rows.iter_mut().enumerate() {
            if i != rank {
                let factor = row[col].clone();
                for j in col..4 {
                    row[j] -= &factor * &pivot[j];
                }
            }
        }
        pivots.push(col);
    }
    if rows[pivots.len()..].iter().any(|row| row[3] != zero()) {
        return None;
    }
    let mut origin = [zero(), zero(), zero()];
    for (i, &col) in pivots.iter().enumerate() {
        origin[col] = rows[i][3].clone();
    }
    let mut basis = Vec::new();
    for free in (0..3).filter(|col| !pivots.contains(col)) {
        let mut v = unit(free);
        for (i, &col) in pivots.iter().enumerate() {
            v[col] = -rows[i][free].clone();
        }
        basis.push(v);
    }
    Some((origin, basis))
}

fn clip_line(origin: Vector, d: &Vector, constraints: &[(Vec<R>, R)]) -> LinearIntersection {
    let (mut lower, mut upper): (Option<R>, Option<R>) = (None, None);
    for (n, c) in constraints {
        if n[0] == zero() {
            if *c > zero() {
                return LinearIntersection::Empty;
            }
        } else {
            let t = c / &n[0];
            if n[0] > zero() {
                lower = Some(lower.map_or(t.clone(), |x| x.max(t)));
            } else {
                upper = Some(upper.map_or(t.clone(), |x| x.min(t)));
            }
        }
    }
    match (lower, upper) {
        (Some(low), Some(high)) => {
            if low > high {
                LinearIntersection::Empty
            } else {
                hull(vec![
                    ExactPoint3(std::array::from_fn(|j| &origin[j] + &d[j] * &low)),
                    ExactPoint3(std::array::from_fn(|j| &origin[j] + &d[j] * &high)),
                ])
            }
        }
        (None, None) => {
            let k = d.iter().position(|x| *x != zero()).expect("line rank");
            let direction = d.clone().map(|x| x / &d[k]);
            LinearIntersection::Line {
                origin: ExactPoint3(std::array::from_fn(|j| {
                    &origin[j] - &direction[j] * &origin[k]
                })),
                direction,
            }
        }
        _ => unreachable!("a bounded operand bounds both ends; affine operands add no halfspaces"),
    }
}

fn hull(mut points: Vec<ExactPoint3>) -> LinearIntersection {
    points.sort();
    points.dedup();
    match points.len() {
        0 => return LinearIntersection::Empty,
        1 => return LinearIntersection::Point(points.remove(0)),
        _ => {}
    }
    let d = sub(&points[1].0, &points[0].0);
    let normal = points[2..]
        .iter()
        .map(|p| cross(&d, &sub(&p.0, &points[0].0)))
        .find(|n| n.iter().any(|x| *x != zero()));
    let Some(n) = normal else {
        return LinearIntersection::Segment([points[0].clone(), points.last().unwrap().clone()]);
    };
    let drop = n.iter().position(|x| *x != zero()).unwrap();
    let axes: Vec<_> = (0..3).filter(|i| *i != drop).collect();
    let (u, v) = (axes[0], axes[1]);
    points.sort_by(|a, b| a.0[u].cmp(&b.0[u]).then(a.0[v].cmp(&b.0[v])));
    let turn = |a: &ExactPoint3, b: &ExactPoint3, c: &ExactPoint3| {
        (&b.0[u] - &a.0[u]) * (&c.0[v] - &a.0[v]) - (&b.0[v] - &a.0[v]) * (&c.0[u] - &a.0[u])
    };
    let chain = |input: Vec<&ExactPoint3>| {
        let mut out: Vec<ExactPoint3> = Vec::new();
        for p in input {
            while out.len() >= 2 && turn(&out[out.len() - 2], &out[out.len() - 1], p) <= zero() {
                out.pop();
            }
            out.push(p.clone());
        }
        out.pop();
        out
    };
    let mut result = chain(points.iter().collect());
    result.extend(chain(points.iter().rev().collect()));
    let start = result
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.cmp(b))
        .unwrap()
        .0;
    result.rotate_left(start);
    if result[1] > result[result.len() - 1] {
        result[1..].reverse();
    }
    LinearIntersection::Polygon(result)
}
fn zero() -> R {
    R::from_integer(BigInt::from(0))
}
fn unit(j: usize) -> Vector {
    std::array::from_fn(|i| R::from_integer(BigInt::from(u8::from(i == j))))
}
fn sub(a: &Vector, b: &Vector) -> Vector {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn dot(a: &Vector, b: &Vector) -> R {
    (0..3).map(|i| &a[i] * &b[i]).sum()
}
fn cross(a: &Vector, b: &Vector) -> Vector {
    std::array::from_fn(|i| &a[(i + 1) % 3] * &b[(i + 2) % 3] - &a[(i + 2) % 3] * &b[(i + 1) % 3])
}
