//! The case protocol of `curve-curve-cases.txt` (S7d) and the kernel's rows
//! for it, shared by `curve_curve_probe` and `tests/curve_curve.rs`:
//! `empty`, `coincident`, `limit`, `overlap s0 s1` (a spline's spans on the
//! conic), or `point sa_lo sa_hi sb_lo sb_hi x_lo x_hi y_lo y_hi z_lo z_hi
//! crossing|tangent`, sorted by the first parameter.
#[path = "curve_surface_protocol.rs"]
#[allow(dead_code)]
mod curve_surface_protocol;
use curve_surface_protocol::{curve, Curve};
use rusty_occt::intersection::{
    curve_curve, spline_curve, AnalyticCurve, CurveCurveIntersection, CurveCurvePoint,
};
use rusty_occt::{BSplineCurve3, Error};

pub struct Case {
    pub name: String,
    /// A spline (S7d.2) is always the first curve.
    pub spline: Option<BSplineCurve3>,
    /// The first curve (unless a spline) and the second.
    pub curves: [Option<AnalyticCurve>; 2],
}

/// A spline, or an analytic curve.
fn split(c: Curve) -> (Option<BSplineCurve3>, Option<AnalyticCurve>) {
    match c {
        Curve::Edge(e) => (None, Some(AnalyticCurve::Edge(e))),
        Curve::Conic(c) => (None, Some(AnalyticCurve::Conic(c))),
        Curve::Spline(s) => (Some(s), None),
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
            let (spline, first) = split(curve(&lines[1][1..]));
            let (none, second) = split(curve(&lines[2][1..]));
            assert!(none.is_none(), "a spline second");
            Case {
                name: lines[0][1].to_string(),
                spline,
                curves: [first, second],
            }
        })
        .collect()
}

fn point_row(p: &CurveCurvePoint) -> String {
    let mut words: Vec<String> = p
        .parameters
        .iter()
        .chain(p.point.iter())
        .map(|[lo, hi]| format!("{lo:?} {hi:?}"))
        .collect();
    words.push(if p.tangent { "tangent" } else { "crossing" }.into());
    format!("point {}", words.join(" "))
}

/// The kernel's rows for a case; `Err` for an unexpected error.
pub fn rows(case: &Case) -> Result<Vec<String>, Error> {
    let second = case.curves[1].as_ref().unwrap();
    if let Some(spline) = &case.spline {
        let found = match spline_curve(spline, second) {
            Err(Error::ComputationLimit(_)) => return Ok(vec!["limit".into()]),
            other => other?,
        };
        let mut rows: Vec<(f64, String)> = found
            .overlaps
            .iter()
            .map(|[a, b]| (*a, format!("overlap {a:?} {b:?}")))
            .chain(
                found
                    .points
                    .iter()
                    .map(|p| (p.parameters[0][0], point_row(p))),
            )
            .collect();
        rows.sort_by(|a, b| a.0.total_cmp(&b.0));
        return Ok(if rows.is_empty() {
            vec!["empty".into()]
        } else {
            rows.into_iter().map(|r| r.1).collect()
        });
    }
    Ok(
        match curve_curve(case.curves[0].as_ref().unwrap(), second) {
            Ok(CurveCurveIntersection::Empty) => vec!["empty".into()],
            Ok(CurveCurveIntersection::Coincident) => vec!["coincident".into()],
            Ok(CurveCurveIntersection::Points(points)) => points.iter().map(point_row).collect(),
            Err(Error::ComputationLimit(_)) => vec!["limit".into()],
            Err(e) => return Err(e),
        },
    )
}
