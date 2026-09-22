//! Small geometric values; precision is explicit at modeling boundaries.
use crate::{Error, Result};
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerance {
    linear: f64,
    angular: f64,
}

impl Default for Tolerance {
    fn default() -> Self {
        // OCCT Precision::Confusion / Precision::Angular conventions.
        Self {
            linear: 1.0e-7,
            angular: 1.0e-12,
        }
    }
}

impl Tolerance {
    pub fn new(linear: f64, angular: f64) -> Result<Self> {
        if !linear.is_finite()
            || linear <= 0.0
            || !angular.is_finite()
            || angular <= 0.0
            || angular >= 1.0
        {
            return Err(Error::InvalidTolerance);
        }
        Ok(Self { linear, angular })
    }
    pub fn linear(self) -> f64 {
        self.linear
    }
    pub fn angular(self) -> f64 {
        self.angular
    }

    pub(crate) fn resolve(self, coordinates: &[f64]) -> Result<()> {
        for value in coordinates {
            finite(*value, "coordinate")?;
            if value.abs() * (16.0 * f64::EPSILON) > self.linear {
                return Err(Error::PrecisionLoss);
            }
        }
        Ok(())
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Point2 {
    pub x: f64,
    pub y: f64,
}

impl Point2 {
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    pub fn distance(self, other: Self) -> f64 {
        (self.x - other.x).hypot(self.y - other.y)
    }
    pub(crate) fn checked(self, tolerance: Tolerance) -> Result<Self> {
        tolerance.resolve(&[self.x, self.y])?;
        Ok(self)
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub const ORIGIN: Self = Self::new(0.0, 0.0, 0.0);
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn distance(self, other: Self) -> f64 {
        (self - other).length()
    }
    pub fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
    pub(crate) fn checked(self, tolerance: Tolerance) -> Result<Self> {
        tolerance.resolve(&self.to_array())?;
        Ok(self)
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    pub fn cross(self, other: Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
    pub fn length(self) -> f64 {
        self.x.hypot(self.y).hypot(self.z)
    }
    pub fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
    pub fn normalized(self) -> Result<Self> {
        for component in self.to_array() {
            finite(component, "direction")?;
        }
        // Scale first: neither very large nor subnormal nonzero directions
        // should overflow/underflow while being normalized.
        let scale = self.x.abs().max(self.y.abs()).max(self.z.abs());
        if scale == 0.0 {
            return Err(Error::Degenerate("direction"));
        }
        let scaled = self / scale;
        Ok(scaled / scaled.length())
    }
}

impl Add<Vec3> for Point3 {
    type Output = Self;
    fn add(self, rhs: Vec3) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}
impl Sub<Point3> for Point3 {
    type Output = Vec3;
    fn sub(self, rhs: Self) -> Vec3 {
        Vec3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}
impl Add for Vec3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}
impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}
impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.0
    }
}
impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}
impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self::new(self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

/// A right-handed, orthonormal plane frame. Fields cannot be changed independently.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame3 {
    origin: Point3,
    x: Vec3,
    y: Vec3,
    normal: Vec3,
}

impl Frame3 {
    pub const fn xy() -> Self {
        Self {
            origin: Point3::ORIGIN,
            x: Vec3::X,
            y: Vec3::Y,
            normal: Vec3::Z,
        }
    }
    pub fn new(origin: Point3, normal: Vec3, x_hint: Vec3, tolerance: Tolerance) -> Result<Self> {
        origin.checked(tolerance)?;
        let normal = normal.normalized()?;
        let hint = x_hint.normalized()?;
        let y = normal.cross(hint);
        if y.length() <= tolerance.angular {
            return Err(Error::Degenerate("plane axes"));
        }
        let y = y.normalized()?;
        let x = y.cross(normal).normalized()?;
        Ok(Self {
            origin,
            x,
            y,
            normal,
        })
    }
    pub fn origin(self) -> Point3 {
        self.origin
    }
    pub fn x(self) -> Vec3 {
        self.x
    }
    pub fn y(self) -> Vec3 {
        self.y
    }
    pub fn normal(self) -> Vec3 {
        self.normal
    }
    pub fn point(self, point: Point2, height: f64) -> Point3 {
        self.origin + self.x * point.x + self.y * point.y + self.normal * height
    }
    pub fn coordinates(self, point: Point3) -> [f64; 3] {
        let vector = point - self.origin;
        [
            vector.dot(self.x),
            vector.dot(self.y),
            vector.dot(self.normal),
        ]
    }
    pub fn transformed(self, transform: RigidTransform, tolerance: Tolerance) -> Result<Self> {
        Self::new(
            transform.point(self.origin),
            transform.vector(self.normal),
            transform.vector(self.x),
            tolerance,
        )
    }
}

/// Rotation and translation only. This cannot encode scale or reflection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RigidTransform {
    columns: [Vec3; 3],
    translation: Vec3,
}

impl RigidTransform {
    pub const fn identity() -> Self {
        Self {
            columns: [Vec3::X, Vec3::Y, Vec3::Z],
            translation: Vec3::new(0.0, 0.0, 0.0),
        }
    }
    pub fn translation(vector: Vec3) -> Result<Self> {
        for value in vector.to_array() {
            finite(value, "translation")?;
        }
        Ok(Self {
            translation: vector,
            ..Self::identity()
        })
    }
    pub fn rotation(origin: Point3, axis: Vec3, angle: f64) -> Result<Self> {
        for value in origin.to_array() {
            finite(value, "rotation origin")?;
        }
        finite(angle, "rotation angle")?;
        let axis = axis.normalized()?;
        let (sine, cosine) = angle.sin_cos();
        let rotate =
            |v: Vec3| v * cosine + axis.cross(v) * sine + axis * (axis.dot(v) * (1.0 - cosine));
        let mut result = Self {
            columns: [rotate(Vec3::X), rotate(Vec3::Y), rotate(Vec3::Z)],
            ..Self::identity()
        };
        result.translation = origin - result.point(origin);
        for value in result.translation.to_array() {
            finite(value, "rotation translation")?;
        }
        Ok(result)
    }
    pub fn vector(self, vector: Vec3) -> Vec3 {
        self.columns[0] * vector.x + self.columns[1] * vector.y + self.columns[2] * vector.z
    }
    pub fn point(self, point: Point3) -> Point3 {
        Point3::ORIGIN + self.vector(point - Point3::ORIGIN) + self.translation
    }
    pub fn inverse(self) -> Result<Self> {
        let [x, y, z] = self.columns;
        let mut result = Self {
            columns: [
                Vec3::new(x.x, y.x, z.x),
                Vec3::new(x.y, y.y, z.y),
                Vec3::new(x.z, y.z, z.z),
            ],
            ..Self::identity()
        };
        result.translation = -result.vector(self.translation);
        for value in result.translation.to_array() {
            finite(value, "inverse translation")?;
        }
        Ok(result)
    }
    /// Apply `self`, then `next`.
    pub fn then(self, next: Self) -> Result<Self> {
        let result = Self {
            columns: self.columns.map(|column| next.vector(column)),
            translation: next.vector(self.translation) + next.translation,
        };
        for value in result.translation.to_array() {
            finite(value, "composed translation")?;
        }
        Ok(result)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds3 {
    pub min: Point3,
    pub max: Point3,
}

impl Bounds3 {
    /// Conservative overlap of finite, closed axis-aligned boxes.
    /// Each box gets the supplied linear gap, so a separation of at most
    /// twice that gap still overlaps. This is broad phase, not solid collision.
    pub fn intersects(self, other: Self, tolerance: Tolerance) -> Result<bool> {
        for bounds in [self, other] {
            bounds.min.checked(tolerance)?;
            bounds.max.checked(tolerance)?;
            if bounds
                .min
                .to_array()
                .into_iter()
                .zip(bounds.max.to_array())
                .any(|(lo, hi)| lo > hi)
            {
                return Err(Error::Degenerate("inverted bounds"));
            }
        }
        // OCCT Bnd_Box.cxx, Bnd_Box::IsOut: finite/non-open fast path.
        // Both gaps contribute, and touching intervals are not separated.
        let delta = 2.0 * tolerance.linear();
        let (a_min, a_max) = (self.min.to_array(), self.max.to_array());
        let (b_min, b_max) = (other.min.to_array(), other.max.to_array());
        Ok(!(0..3).any(|i| a_min[i] - b_max[i] > delta || b_min[i] - a_max[i] > delta))
    }

    pub fn contains(self, point: Point3, tolerance: Tolerance) -> bool {
        let p = point.to_array();
        let lo = self.min.to_array();
        let hi = self.max.to_array();
        (0..3).all(|i| p[i] >= lo[i] - tolerance.linear && p[i] <= hi[i] + tolerance.linear)
    }
}

pub(crate) fn finite(value: f64, name: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::NonFinite(name))
    }
}

/// Neumaier summation makes profile moments less sensitive to vertex order.
pub(crate) fn sum(values: impl Iterator<Item = f64>) -> f64 {
    let mut total: f64 = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = total + value;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    total + correction
}
