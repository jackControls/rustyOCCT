use std::fmt;

pub type Result<T> = std::result::Result<T, Error>;

/// Invalid geometry is reported before a solid is made available to callers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    NonFinite(&'static str),
    InvalidTolerance,
    Degenerate(&'static str),
    PrecisionLoss,
    Unrepresentable(&'static str),
    SelfIntersection,
    InvalidHole(usize),
    IntersectingBoundaries(usize, usize),
    LimitExceeded(&'static str),
    InvalidTopology(&'static str),
    InvalidCurve(&'static str),
    InvalidSpline(&'static str),
    InvalidSurface(&'static str),
    OutOfDomain(&'static str),
    DiscontinuousDerivative,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite(what) => write!(f, "{what} must be finite"),
            Self::InvalidTolerance => write!(f, "tolerances must be finite and positive; angular tolerance must be less than one radian"),
            Self::Degenerate(what) => write!(f, "degenerate {what}"),
            Self::PrecisionLoss => write!(f, "coordinates cannot resolve the requested linear tolerance"),
            Self::Unrepresentable(what) => write!(f, "{what} cannot be enclosed by finite binary64 values"),
            Self::SelfIntersection => write!(f, "polygon crosses or touches itself within tolerance"),
            Self::InvalidHole(i) => write!(f, "hole {i} is not strictly inside the outer boundary"),
            Self::IntersectingBoundaries(a, b) => write!(f, "boundaries {a} and {b} intersect, touch, or overlap (outer boundary is 0)"),
            Self::LimitExceeded(what) => write!(f, "{what} exceeds the supported input limit"),
            Self::InvalidTopology(what) => write!(f, "invalid topology: {what}"),
            Self::InvalidCurve(what) => write!(f, "invalid curve: {what}"),
            Self::InvalidSpline(what) => write!(f, "invalid spline: {what}"),
            Self::InvalidSurface(what) => write!(f, "invalid surface: {what}"),
            Self::OutOfDomain(what) => write!(f, "outside the supported domain: {what}"),
            Self::DiscontinuousDerivative => write!(f, "requested derivative is discontinuous at the knot or seam; select a side"),
        }
    }
}

impl std::error::Error for Error {}
