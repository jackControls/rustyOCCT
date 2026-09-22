use std::fmt;

pub type Result<T> = std::result::Result<T, Error>;

/// Invalid geometry is reported before a solid is made available to callers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    NonFinite(&'static str),
    InvalidTolerance,
    Degenerate(&'static str),
    PrecisionLoss,
    SelfIntersection,
    InvalidHole(usize),
    IntersectingBoundaries(usize, usize),
    LimitExceeded(&'static str),
    InvalidTopology(&'static str),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite(what) => write!(f, "{what} must be finite"),
            Self::InvalidTolerance => write!(f, "tolerances must be finite and positive; angular tolerance must be less than one radian"),
            Self::Degenerate(what) => write!(f, "degenerate {what}"),
            Self::PrecisionLoss => write!(f, "coordinates cannot resolve the requested linear tolerance"),
            Self::SelfIntersection => write!(f, "polygon crosses or touches itself within tolerance"),
            Self::InvalidHole(i) => write!(f, "hole {i} is not strictly inside the outer boundary"),
            Self::IntersectingBoundaries(a, b) => write!(f, "boundaries {a} and {b} intersect, touch, or overlap (outer boundary is 0)"),
            Self::LimitExceeded(what) => write!(f, "{what} exceeds the supported input limit"),
            Self::InvalidTopology(what) => write!(f, "invalid topology: {what}"),
        }
    }
}

impl std::error::Error for Error {}
