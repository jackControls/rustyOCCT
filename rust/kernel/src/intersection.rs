//! Certified intersections of represented lines, segments, planes and triangles.
//!
//! Decisions use exact integer arithmetic. Constructed coordinates and line
//! parameters are enclosed between adjacent finite binary64 numbers (equal
//! when exactly representable). No modeling tolerance or snapping is applied.
//! A rounded representative is not itself an exact point on both operands;
//! use its bounds when propagating construction uncertainty.
//!
//! OCCT reference: IntAna_IntConicQuad::Perform(gp_Lin, gp_Pln). The affine
//! substitution is shared; OCCT's angular-tolerance parallelism is deliberately
//! replaced by exact sidedness here. See rust/MATHEMATICS.md and SOURCE_MAP.md.
pub use crate::interval::ScalarInterval;
use crate::{exact, Bounds3, Error, Point3, Result};
use num_bigint::{BigInt, BigUint, Sign};
use std::cmp::Ordering;

mod curved;
mod exact_spline;
mod linear_sets;
mod spline_linear;
mod spline_plane;
mod spline_quadric;
mod spline_surface;
pub use crate::proximity::LinearPrimitive3;
pub use curved::{
    line_circle, line_cylinder, line_sphere, segment_circle, segment_cylinder, segment_sphere,
    Circle3, ContactKind, CurveHit, CurvedIntersection, Cylinder3, Sphere3,
};
pub use exact_spline::{
    exact_spline_cylinder, exact_spline_cylinder_in, exact_spline_cylinder_in_with_options,
    exact_spline_cylinder_with_options, exact_spline_plane, exact_spline_plane_in,
    exact_spline_plane_in_with_options, exact_spline_plane_with_options, exact_spline_sphere,
    exact_spline_sphere_in, exact_spline_sphere_in_with_options, exact_spline_sphere_with_options,
};
pub use linear_sets::{linear_intersection, ExactPoint3, LinearIntersection};
pub use spline_linear::{
    exact_spline_line, exact_spline_line_in, exact_spline_line_in_with_options,
    exact_spline_line_with_options, exact_spline_segment, exact_spline_segment_in,
    exact_spline_segment_in_with_options, exact_spline_segment_with_options, spline_line,
    spline_line_in, spline_line_in_with_options, spline_line_with_options, spline_segment,
    spline_segment_in, spline_segment_in_with_options, spline_segment_with_options,
    SplineLinearIntersection, SplineLinearOptions, SplineLinearOverlap, SplineLinearPoint,
};
pub use spline_plane::{
    spline_plane, spline_plane_in, spline_plane_in_with_options, spline_plane_with_options,
    SplinePlaneContact, SplinePlaneIntersection, SplinePlaneOptions, SplinePlaneOverlap,
    SplinePlanePoint,
};
pub use spline_quadric::{
    spline_cylinder, spline_cylinder_in, spline_cylinder_in_with_options,
    spline_cylinder_with_options, spline_sphere, spline_sphere_in, spline_sphere_in_with_options,
    spline_sphere_with_options,
};
pub use spline_surface::{
    ExactSplineSurfaceIntersection as ExactSplinePlaneIntersection,
    ExactSplineSurfaceOverlap as ExactSplinePlaneOverlap,
    ExactSplineSurfacePoint as ExactSplinePlanePoint,
};
pub use spline_surface::{
    ExactSplineSurfaceIntersection, ExactSplineSurfaceOverlap, ExactSplineSurfacePoint,
    SplineSurfaceContact, SplineSurfaceIntersection, SplineSurfaceOptions, SplineSurfaceOverlap,
    SplineSurfacePoint,
};

/// Plane defined by three exactly noncollinear finite points.
/// No rounded unit normal is used to define its geometry.
#[derive(Debug, Clone)]
pub struct Plane3 {
    vertices: [Point3; 3],
    anchor: exact::Vector,
    normal: exact::Vector,
}

impl Plane3 {
    pub fn through_points(a: Point3, b: Point3, c: Point3) -> Result<Self> {
        let vertices = [a, b, c];
        let [anchor, b, c] = [exact::point(a)?, exact::point(b)?, exact::point(c)?];
        let normal = exact::cross(&exact::sub(&b, &anchor), &exact::sub(&c, &anchor));
        if normal.iter().all(exact::zero) {
            return Err(Error::Degenerate("plane defining points"));
        }
        Ok(Self {
            vertices,
            anchor,
            normal,
        })
    }

    pub fn defining_points(&self) -> [Point3; 3] {
        self.vertices
    }

    fn distance(&self, p: &exact::Vector) -> BigInt {
        exact::dot(&self.normal, &exact::sub(p, &self.anchor))
    }
}

/// A closed, exactly nondegenerate triangle, including all three edges.
#[derive(Debug, Clone)]
pub struct Triangle3 {
    plane: Plane3,
}

impl Triangle3 {
    pub fn new(a: Point3, b: Point3, c: Point3) -> Result<Self> {
        Ok(Self {
            plane: Plane3::through_points(a, b, c)?,
        })
    }
    pub fn vertices(&self) -> [Point3; 3] {
        self.plane.defining_points()
    }
    pub fn plane(&self) -> &Plane3 {
        &self.plane
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntersectionPoint {
    coordinates: [ScalarInterval; 3],
    parameter: ScalarInterval,
}

impl IntersectionPoint {
    pub fn position(self) -> Point3 {
        let [x, y, z] = self.coordinates.map(ScalarInterval::representative);
        Point3::new(x, y, z)
    }
    pub fn bounds(self) -> Bounds3 {
        let [lx, ly, lz] = self.coordinates.map(ScalarInterval::lower);
        let [ux, uy, uz] = self.coordinates.map(ScalarInterval::upper);
        Bounds3 {
            min: Point3::new(lx, ly, lz),
            max: Point3::new(ux, uy, uz),
        }
    }
    /// Encloses t in p + t(q-p), using the original input endpoint order.
    pub fn parameter(self) -> ScalarInterval {
        self.parameter
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LinePlaneIntersection {
    Disjoint,
    Contained,
    Point(IntersectionPoint),
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SegmentPlaneIntersection {
    Disjoint,
    Contained,
    Point(IntersectionPoint),
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SegmentTriangleIntersection {
    Disjoint,
    Point(IntersectionPoint),
    /// Coplanar overlap endpoints, in the direction of the input segment.
    Overlap {
        start: IntersectionPoint,
        end: IntersectionPoint,
    },
}

/// Infinite line through p and q. Identical points are a degenerate line.
/// Returns `Unrepresentable` if the unique point or its parameter cannot be
/// enclosed with finite binary64 endpoints, even when its exact value exists.
pub fn line_plane(p: Point3, q: Point3, plane: &Plane3) -> Result<LinePlaneIntersection> {
    let (p, q, dp, dq) = distances(p, q, plane)?;
    if p == q {
        return Err(Error::Degenerate("line defining points"));
    }
    let denominator = &dp - &dq;
    if exact::zero(&denominator) {
        return Ok(if exact::zero(&dp) {
            LinePlaneIntersection::Contained
        } else {
            LinePlaneIntersection::Disjoint
        });
    }
    Ok(LinePlaneIntersection::Point(construct(
        &p,
        &q,
        &Ratio::new(dp, denominator),
    )?))
}

/// Closed segment. Equal endpoints are treated as a point at parameter zero.
pub fn segment_plane(p: Point3, q: Point3, plane: &Plane3) -> Result<SegmentPlaneIntersection> {
    let (p, q, dp, dq) = distances(p, q, plane)?;
    if p == q {
        return Ok(if exact::zero(&dp) {
            SegmentPlaneIntersection::Point(construct(&p, &q, &Ratio::endpoint(0))?)
        } else {
            SegmentPlaneIntersection::Disjoint
        });
    }
    if exact::zero(&dp) && exact::zero(&dq) {
        return Ok(SegmentPlaneIntersection::Contained);
    }
    if dp.sign() == dq.sign() {
        return Ok(SegmentPlaneIntersection::Disjoint);
    }
    let denominator = &dp - &dq;
    Ok(SegmentPlaneIntersection::Point(construct(
        &p,
        &q,
        &Ratio::new(dp, denominator),
    )?))
}

/// Closed segment against a closed triangle, including coplanar clipping,
/// tangencies, vertex hits and zero-length segments. Exact decisions precede
/// coordinate rounding. Reversing triangle winding does not change the region.
pub fn segment_triangle(
    p: Point3,
    q: Point3,
    triangle: &Triangle3,
) -> Result<SegmentTriangleIntersection> {
    let (p, q, dp, dq) = distances(p, q, &triangle.plane)?;
    if dp.sign() == dq.sign() && !exact::zero(&dp) {
        return Ok(SegmentTriangleIntersection::Disjoint);
    }
    // Dropping an axis with nonzero exact normal leaves a nondegenerate 2D
    // triangle. Neither a rounded normal nor a floating-point projection enters.
    let drop = triangle
        .plane
        .normal
        .iter()
        .position(|x| !exact::zero(x))
        .unwrap();
    let axes = [[1, 2], [0, 2], [0, 1]][drop];
    let [a, b, c] = triangle
        .vertices()
        .map(|v| exact::point(v).expect("validated triangle"));
    let edge = |a: &exact::Vector, b: &exact::Vector, p: &exact::Vector| {
        (&b[axes[0]] - &a[axes[0]]) * (&p[axes[1]] - &a[axes[1]])
            - (&b[axes[1]] - &a[axes[1]]) * (&p[axes[0]] - &a[axes[0]])
    };
    let positive = edge(&a, &b, &c).sign() == Sign::Plus;
    let values = [(&a, &b), (&b, &c), (&c, &a)].map(|(a, b)| {
        let (ep, eq) = (edge(a, b, &p), edge(a, b, &q));
        if positive {
            (ep, eq)
        } else {
            (-ep, -eq)
        }
    });
    if !exact::zero(&dp) || !exact::zero(&dq) {
        let denominator = &dp - &dq;
        let t = Ratio::new(dp, denominator);
        for (ep, eq) in values {
            if (&ep * &t.den + (eq - &ep) * &t.num).sign() == Sign::Minus {
                return Ok(SegmentTriangleIntersection::Disjoint);
            }
        }
        return Ok(SegmentTriangleIntersection::Point(construct(&p, &q, &t)?));
    }
    // Exact half-plane clipping of t in [0,1]. No division until construction.
    let (mut low, mut high) = (Ratio::endpoint(0), Ratio::endpoint(1));
    for (ep, eq) in values {
        let (outside_p, outside_q) = (ep.sign() == Sign::Minus, eq.sign() == Sign::Minus);
        if outside_p && outside_q {
            return Ok(SegmentTriangleIntersection::Disjoint);
        }
        if outside_p {
            let t = Ratio::new(-&ep, eq - &ep);
            if t.cmp(&low) == Ordering::Greater {
                low = t;
            }
        } else if outside_q {
            let t = Ratio::new(ep.clone(), ep - eq);
            if t.cmp(&high) == Ordering::Less {
                high = t;
            }
        }
    }
    match low.cmp(&high) {
        Ordering::Greater => Ok(SegmentTriangleIntersection::Disjoint),
        Ordering::Equal => Ok(SegmentTriangleIntersection::Point(construct(&p, &q, &low)?)),
        Ordering::Less if p == q => Ok(SegmentTriangleIntersection::Point(construct(
            &p,
            &q,
            &Ratio::endpoint(0),
        )?)),
        Ordering::Less => Ok(SegmentTriangleIntersection::Overlap {
            start: construct(&p, &q, &low)?,
            end: construct(&p, &q, &high)?,
        }),
    }
}

fn distances(
    p: Point3,
    q: Point3,
    plane: &Plane3,
) -> Result<(exact::Vector, exact::Vector, BigInt, BigInt)> {
    let (p, q) = (exact::point(p)?, exact::point(q)?);
    let (dp, dq) = (plane.distance(&p), plane.distance(&q));
    Ok((p, q, dp, dq))
}

struct Ratio {
    num: BigInt,
    den: BigInt,
}
impl Ratio {
    fn new(num: BigInt, den: BigInt) -> Self {
        debug_assert!(!exact::zero(&den));
        if den.sign() == Sign::Minus {
            Self {
                num: -num,
                den: -den,
            }
        } else {
            Self { num, den }
        }
    }
    fn endpoint(value: u8) -> Self {
        Self {
            num: BigInt::from(value),
            den: BigInt::from(1),
        }
    }
    fn cmp(&self, other: &Self) -> Ordering {
        (&self.num * &other.den).cmp(&(&other.num * &self.den))
    }
}

fn construct(p: &exact::Vector, q: &exact::Vector, t: &Ratio) -> Result<IntersectionPoint> {
    let coordinate = |i: usize| {
        enclose(
            &(&p[i] * &t.den + (&q[i] - &p[i]) * &t.num),
            &t.den,
            -1074,
            "intersection coordinate",
        )
    };
    Ok(IntersectionPoint {
        coordinates: [coordinate(0)?, coordinate(1)?, coordinate(2)?],
        parameter: enclose(&t.num, &t.den, 0, "intersection parameter")?,
    })
}

/// Compare n/d * 2^scale to a nonnegative finite binary64 candidate, exactly.
fn compare(n: &BigUint, d: &BigUint, scale: i32, bits: u64) -> Ordering {
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1_u64 << 52) - 1);
    let (mantissa, power) = if exponent == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1_u64 << 52), exponent - 1075)
    };
    let right = d * mantissa;
    if scale >= power {
        (n << (scale - power) as usize).cmp(&right)
    } else {
        n.cmp(&(right << (power - scale) as usize))
    }
}

fn enclose(num: &BigInt, den: &BigInt, scale: i32, what: &'static str) -> Result<ScalarInterval> {
    debug_assert_eq!(den.sign(), Sign::Plus);
    let (n, d) = (num.magnitude(), den.magnitude());
    let negative = num.sign() == Sign::Minus;
    crate::interval::enclose(
        |x| {
            if x == 0. {
                return num.cmp(&BigInt::from(0));
            }
            if x.is_sign_negative() != negative {
                return if negative {
                    Ordering::Less
                } else {
                    Ordering::Greater
                };
            }
            let order = compare(n, d, scale, x.abs().to_bits());
            if negative {
                order.reverse()
            } else {
                order
            }
        },
        what,
    )
}
