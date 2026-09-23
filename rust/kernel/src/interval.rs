use crate::{Error, Result};
use std::cmp::Ordering;

/// Smallest closed finite binary64 interval containing a certified real value.
/// Distinct algebraic roots can have overlapping intervals: an enclosure does
/// not by itself establish root identity or isolate it from another root.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScalarInterval {
    lower: f64,
    upper: f64,
}

impl ScalarInterval {
    pub fn lower(self) -> f64 {
        self.lower
    }
    pub fn upper(self) -> f64 {
        self.upper
    }
    /// A finite representative inside the interval, not an exact construction.
    pub fn representative(self) -> f64 {
        self.lower + (self.upper - self.lower) * 0.5
    }
}

/// The callback returns the exact ordering of the represented real and x.
/// Search the ordered finite magnitude bit patterns; at most 63 iterations.
pub(crate) fn enclose(
    compare: impl Fn(f64) -> Ordering,
    what: &'static str,
) -> Result<ScalarInterval> {
    let negative = match compare(0.) {
        Ordering::Equal => {
            return Ok(ScalarInterval {
                lower: 0.,
                upper: 0.,
            })
        }
        Ordering::Less => true,
        Ordering::Greater => false,
    };
    let magnitude_cmp = |x: f64| {
        if negative {
            compare(-x).reverse()
        } else {
            compare(x)
        }
    };
    const MAX: u64 = 0x7fef_ffff_ffff_ffff;
    if magnitude_cmp(f64::MAX) == Ordering::Greater {
        return Err(Error::Unrepresentable(what));
    }
    let (mut low, mut high) = (0, MAX);
    while low < high {
        let mid = low + (high - low).div_ceil(2);
        if magnitude_cmp(f64::from_bits(mid)) == Ordering::Less {
            high = mid - 1;
        } else {
            low = mid;
        }
    }
    let upper = if magnitude_cmp(f64::from_bits(low)) == Ordering::Equal {
        low
    } else {
        low + 1
    };
    let (low, high) = (f64::from_bits(low), f64::from_bits(upper));
    Ok(if negative {
        ScalarInterval {
            lower: -high,
            upper: -low,
        }
    } else {
        ScalarInterval {
            lower: low,
            upper: high,
        }
    })
}
