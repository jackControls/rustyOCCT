//! Analytic curved primitives. Exact substitution into their implicit equations
//! preserves tangency, root identity, and closed segment membership before any
//! rounding. OCCT references: IntAna_IntConicQuad line/quadric substitution,
//! IntAna_Quadric sphere/cylinder equations, IntAna2d line/circle cases.
use super::IntersectionPoint;
use crate::polynomial::{solve_integer, QuadraticRoot, QuadraticRoots};
use crate::{exact, math::finite, Error, Point3, Result, Vec3};
use num_bigint::BigInt;
use std::cmp::Ordering;

/// Sphere surface with an exact represented center and strictly positive radius.
#[derive(Debug, Clone, Copy)]
pub struct Sphere3 {
    center: Point3,
    radius: f64,
}
impl Sphere3 {
    pub fn new(center: Point3, radius: f64) -> Result<Self> {
        exact::point(center)?;
        finite(radius, "sphere radius")?;
        if radius <= 0. {
            return Err(Error::Degenerate("sphere radius"));
        }
        Ok(Self { center, radius })
    }
    pub fn center(self) -> Point3 {
        self.center
    }
    pub fn radius(self) -> f64 {
        self.radius
    }
}

/// Infinite cylinder surface. Its finite nonzero axis vector is interpreted
/// exactly and need not be normalized. End caps are not part of this primitive.
#[derive(Debug, Clone, Copy)]
pub struct Cylinder3 {
    origin: Point3,
    axis: Vec3,
    radius: f64,
}
impl Cylinder3 {
    pub fn new(origin: Point3, axis: Vec3, radius: f64) -> Result<Self> {
        Sphere3::new(origin, radius)?;
        nonzero_vector(axis, "cylinder axis")?;
        Ok(Self {
            origin,
            axis,
            radius,
        })
    }
    pub fn origin(self) -> Point3 {
        self.origin
    }
    pub fn axis(self) -> Vec3 {
        self.axis
    }
    pub fn radius(self) -> f64 {
        self.radius
    }
}

/// Circle in the exact plane `(x-center) dot normal = 0`, at distance radius
/// from its center. The finite nonzero normal is not rounded by normalization.
#[derive(Debug, Clone, Copy)]
pub struct Circle3 {
    sphere: Sphere3,
    normal: Vec3,
}
impl Circle3 {
    pub fn new(center: Point3, normal: Vec3, radius: f64) -> Result<Self> {
        let sphere = Sphere3::new(center, radius)?;
        nonzero_vector(normal, "circle normal")?;
        Ok(Self { sphere, normal })
    }
    pub fn center(self) -> Point3 {
        self.sphere.center()
    }
    pub fn radius(self) -> f64 {
        self.sphere.radius()
    }
    pub fn normal(self) -> Vec3 {
        self.normal
    }
}

fn nonzero_vector(vector: Vec3, what: &'static str) -> Result<()> {
    for value in vector.to_array() {
        finite(value, what)?;
    }
    if vector.to_array().iter().all(|&x| x == 0.) {
        return Err(Error::Degenerate(what));
    }
    Ok(())
}
fn integer_vector(v: Vec3) -> exact::Vector {
    v.to_array().map(exact::integer)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactKind {
    /// A simple isolated hit. For a circle this also includes a line crossing
    /// its plane at a point on the circle.
    Crossing,
    /// A double quadratic root: the line is tangent to the primitive.
    Tangent,
    /// A zero-length segment exactly on the primitive; no tangent direction.
    DegeneratePoint,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveHit {
    pub point: IntersectionPoint,
    pub kind: ContactKind,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurvedIntersection {
    Disjoint,
    /// A complete line/segment on an infinite cylinder generator.
    Contained,
    One(CurveHit),
    /// Ordered by exact input-line parameter. Distinct points stay distinct
    /// even if their rounded positions or coordinate enclosures coincide.
    Two {
        first: CurveHit,
        second: CurveHit,
    },
}

/// Intersect the infinite line through p,q with a sphere surface.
/// Equal line endpoints are an error. Unrepresentable output coordinates or
/// affine parameters return an error, never a partial or infinite-valued result.
pub fn line_sphere(p: Point3, q: Point3, sphere: &Sphere3) -> Result<CurvedIntersection> {
    intersect(p, q, Primitive::Sphere(sphere), false)
}
/// Closed segment/sphere-surface intersection (not containment in the ball).
/// A zero-length segment on the surface returns `DegeneratePoint` at t=0.
pub fn segment_sphere(p: Point3, q: Point3, sphere: &Sphere3) -> Result<CurvedIntersection> {
    intersect(p, q, Primitive::Sphere(sphere), true)
}
/// Infinite line against the complete infinite cylinder surface, without caps.
pub fn line_cylinder(p: Point3, q: Point3, cylinder: &Cylinder3) -> Result<CurvedIntersection> {
    intersect(p, q, Primitive::Cylinder(cylinder), false)
}
/// Closed segment against the infinite cylinder surface, including contained
/// generators. Roots outside [0,1] are excluded before constructing any bounds.
pub fn segment_cylinder(p: Point3, q: Point3, cylinder: &Cylinder3) -> Result<CurvedIntersection> {
    intersect(p, q, Primitive::Cylinder(cylinder), true)
}
/// Infinite line/circle intersection, in arbitrary 3D orientation. A line
/// crossing the circle's plane is checked against its radius exactly.
pub fn line_circle(p: Point3, q: Point3, circle: &Circle3) -> Result<CurvedIntersection> {
    intersect(p, q, Primitive::Circle(circle), false)
}
/// Closed segment/circle intersection, including tangent and endpoint hits.
pub fn segment_circle(p: Point3, q: Point3, circle: &Circle3) -> Result<CurvedIntersection> {
    intersect(p, q, Primitive::Circle(circle), true)
}

enum Primitive<'a> {
    Sphere(&'a Sphere3),
    Cylinder(&'a Cylinder3),
    Circle(&'a Circle3),
}

fn intersect(
    p: Point3,
    q: Point3,
    primitive: Primitive<'_>,
    segment: bool,
) -> Result<CurvedIntersection> {
    let (p, q) = (exact::point(p)?, exact::point(q)?);
    let point_segment = p == q;
    if point_segment && !segment {
        return Err(Error::Degenerate("line defining points"));
    }
    let (center, radius) = match primitive {
        Primitive::Sphere(s) => (s.center, s.radius),
        Primitive::Cylinder(c) => (c.origin, c.radius),
        Primitive::Circle(c) => (c.center(), c.radius()),
    };
    let w = exact::sub(&p, &exact::point(center)?);
    let direction = exact::sub(&q, &p);
    let radius = exact::integer(radius);
    let (a, b, c) = if let Primitive::Cylinder(cylinder) = primitive {
        let axis = integer_vector(cylinder.axis);
        let wd = exact::cross(&w, &axis);
        let dd = exact::cross(&direction, &axis);
        // |(p-origin + t*d) x axis|² - r²|axis|² = 0, without normalizing.
        (
            exact::dot(&dd, &dd),
            exact::dot(&wd, &dd) * 2_u8,
            exact::dot(&wd, &wd) - &radius * &radius * exact::dot(&axis, &axis),
        )
    } else {
        // |p-center + t*d|²-r² = 0. Coordinate differences never round.
        (
            exact::dot(&direction, &direction),
            exact::dot(&w, &direction) * 2_u8,
            exact::dot(&w, &w) - &radius * &radius,
        )
    };
    let roots = if let Primitive::Circle(circle) = primitive {
        let normal = integer_vector(circle.normal);
        let origin_distance = exact::dot(&normal, &w);
        let plane_direction = exact::dot(&normal, &direction);
        if exact::zero(&plane_direction) {
            if !exact::zero(&origin_distance) {
                return Ok(CurvedIntersection::Disjoint);
            }
            solve_integer(a, b, c)
        } else {
            let num = -origin_distance;
            let den = plane_direction;
            // Substitute the exact rational plane crossing into the sphere
            // equation; testing a rounded point here would corrupt topology.
            if !exact::zero(&(a * &num * &num + b * &num * &den + c * &den * &den)) {
                return Ok(CurvedIntersection::Disjoint);
            }
            QuadraticRoots::One(QuadraticRoot::rational(num, den, 1))
        }
    } else {
        solve_integer(a, b, c)
    };

    let roots = match roots {
        QuadraticRoots::None => return Ok(CurvedIntersection::Disjoint),
        QuadraticRoots::All if point_segment => {
            let zero = QuadraticRoot::rational(BigInt::from(0), BigInt::from(1), 1);
            return Ok(CurvedIntersection::One(CurveHit {
                point: construct(&p, &direction, &zero)?,
                kind: ContactKind::DegeneratePoint,
            }));
        }
        QuadraticRoots::All => return Ok(CurvedIntersection::Contained),
        QuadraticRoots::One(root) => vec![root],
        QuadraticRoots::Two { lower, upper } => vec![lower, upper],
    };
    let mut hits = Vec::with_capacity(2);
    for root in roots {
        if segment
            && (root.value.compare(0., 0) == Ordering::Less
                || root.value.compare(1., 0) == Ordering::Greater)
        {
            continue;
        }
        hits.push(CurveHit {
            point: construct(&p, &direction, &root)?,
            kind: if root.multiplicity() == 2 {
                ContactKind::Tangent
            } else {
                ContactKind::Crossing
            },
        });
    }
    Ok(match hits.as_slice() {
        [] => CurvedIntersection::Disjoint,
        [hit] => CurvedIntersection::One(*hit),
        [first, second] => CurvedIntersection::Two {
            first: *first,
            second: *second,
        },
        _ => unreachable!("quadratic has at most two distinct roots"),
    })
}

fn construct(
    p: &exact::Vector,
    direction: &exact::Vector,
    root: &QuadraticRoot,
) -> Result<IntersectionPoint> {
    let coordinate = |i: usize| {
        root.value
            .affine(&p[i], &direction[i])
            .bounds(-1074, "intersection coordinate")
    };
    Ok(IntersectionPoint {
        coordinates: [coordinate(0)?, coordinate(1)?, coordinate(2)?],
        parameter: root.value.bounds(0, "intersection parameter")?,
    })
}
