//! The case protocol of `curve-curve-cases.txt` (S7d.1) and the kernel's
//! rows for it, shared by `curve_curve_probe` and `tests/curve_curve.rs`:
//! `empty`, `coincident`, `limit`, or `point sa_lo sa_hi sb_lo sb_hi x_lo x_hi
//! y_lo y_hi z_lo z_hi crossing|tangent`, sorted by the first parameter.
#[path = "curve_surface_protocol.rs"]
#[allow(dead_code)]
mod curve_surface_protocol;
use curve_surface_protocol::{curve, Curve};
use rusty_occt::intersection::{curve_curve, AnalyticCurve, CurveCurveIntersection};
use rusty_occt::Error;

pub struct Case {
    pub name: String,
    pub curves: [AnalyticCurve; 2],
}

fn analytic(c: Curve) -> AnalyticCurve {
    match c {
        Curve::Edge(e) => AnalyticCurve::Edge(e),
        Curve::Conic(c) => AnalyticCurve::Conic(c),
        Curve::Spline(_) => unreachable!("S7d.1 has no splines"),
    }
}

pub fn cases(text: &str) -> Vec<Case> {
    text.split("\nend")
        .filter(|b| !b.trim().is_empty())
        .map(|block| {
            let lines: Vec<Vec<&str>> = block
                .trim()
                .lines()
                .map(|l| l.split_whitespace().collect())
                .collect();
            Case {
                name: lines[0][1].to_string(),
                curves: [
                    analytic(curve(&lines[1][1..])),
                    analytic(curve(&lines[2][1..])),
                ],
            }
        })
        .collect()
}

/// The kernel's rows for a case; `Err` for an unexpected error.
pub fn rows(case: &Case) -> Result<Vec<String>, Error> {
    Ok(match curve_curve(&case.curves[0], &case.curves[1]) {
        Ok(CurveCurveIntersection::Empty) => vec!["empty".into()],
        Ok(CurveCurveIntersection::Coincident) => vec!["coincident".into()],
        Ok(CurveCurveIntersection::Points(points)) => points
            .iter()
            .map(|p| {
                let mut words: Vec<String> = p
                    .parameters
                    .iter()
                    .chain(p.point.iter())
                    .map(|[lo, hi]| format!("{lo:?} {hi:?}"))
                    .collect();
                words.push(if p.tangent { "tangent" } else { "crossing" }.into());
                format!("point {}", words.join(" "))
            })
            .collect(),
        Err(Error::ComputationLimit(_)) => vec!["limit".into()],
        Err(e) => return Err(e),
    })
}
