//! STEP-b of the STEP import track in `REVIEW_NOTES.md`: B-spline curves and
//! surfaces (simple instances, and rational ones as the complex instances
//! OCCT and other writers produce), the file's pcurves on spline surfaces,
//! and where an edge's vertices lie on a spline or on a pcurve's image.
//!
//! Nothing here approximates: a located parameter is only a starting claim,
//! which the validator certifies (the vertices on their edges, each pcurve
//! against its edge at matching fractions), so a wrong location is an import
//! failure with its issues, never a moved edge.
use super::import::{list, real, reference, Build, Named};
use super::part21::{Instance, Parameter};
use crate::occt_brep::read::{BSplineRecord, BSplineSurfaceRecord, Curve2 as Record2};
use crate::{BSplineCurve2, BSplineCurve3, BSplineSurface3, Error, KnotVector, Point2, Point3};

/// The knotless forms (OCCT converts them to knotted B-splines), refused by
/// name for now.
const KNOTLESS: [&str; 6] = [
    "BEZIER_CURVE",
    "UNIFORM_CURVE",
    "QUASI_UNIFORM_CURVE",
    "BEZIER_SURFACE",
    "UNIFORM_SURFACE",
    "QUASI_UNIFORM_SURFACE",
];

/// Whether an instance is a knotted B-spline curve, simple or complex.
pub(super) fn is_curve(i: &Instance) -> bool {
    i.record("B_SPLINE_CURVE_WITH_KNOTS").is_some()
}

/// Whether an instance is a knotted B-spline surface, simple or complex.
pub(super) fn is_surface(i: &Instance) -> bool {
    i.record("B_SPLINE_SURFACE_WITH_KNOTS").is_some()
}

/// The name of a knotless B-spline form the instance has.
pub(super) fn knotless(i: &Instance) -> Option<&'static str> {
    KNOTLESS.iter().find(|n| i.record(n).is_some()).copied()
}

/// Degrees, multiplicities and counts are bounded so malformed data fails
/// rather than allocates; the kernel's constructors apply its own limits.
const LIMIT: i64 = 1 << 20;

fn count(p: &Parameter) -> Named<usize> {
    match p {
        Parameter::Integer(i) if (1..=LIMIT).contains(i) => Ok(*i as usize),
        _ => Err("MalformedEntity"),
    }
}

fn counts(p: &Parameter) -> Named<Vec<usize>> {
    list(p)?.iter().map(count).collect()
}

fn reals(p: &Parameter) -> Named<Vec<f64>> {
    list(p)?.iter().map(real).collect()
}

/// A curve's degree, control points, multiplicities, knots and weights.
struct CurveParts<'a> {
    degree: usize,
    poles: &'a [Parameter],
    multiplicities: Vec<usize>,
    knots: Vec<f64>,
    weights: Option<Vec<f64>>,
}

fn curve_parts(i: &Instance) -> Named<CurveParts<'_>> {
    let with_knots = i
        .record("B_SPLINE_CURVE_WITH_KNOTS")
        .ok_or("MalformedEntity")?;
    // A complex instance: the supertype's record holds the degree and the
    // control points, the subtypes' the knots and the weights.
    if let Some(head) = i.record("B_SPLINE_CURVE") {
        let [degree, poles, _, _, _] = head.parameters.as_slice() else {
            return Err("MalformedEntity");
        };
        let [multiplicities, knots, _] = with_knots.parameters.as_slice() else {
            return Err("MalformedEntity");
        };
        let weights = match i.record("RATIONAL_B_SPLINE_CURVE") {
            Some(r) => match r.parameters.as_slice() {
                [w] => Some(reals(w)?),
                _ => return Err("MalformedEntity"),
            },
            None => None,
        };
        return Ok(CurveParts {
            degree: count(degree)?,
            poles: list(poles)?,
            multiplicities: counts(multiplicities)?,
            knots: reals(knots)?,
            weights,
        });
    }
    let [_, degree, poles, _, _, _, multiplicities, knots, _] = with_knots.parameters.as_slice()
    else {
        return Err("MalformedEntity");
    };
    Ok(CurveParts {
        degree: count(degree)?,
        poles: list(poles)?,
        multiplicities: counts(multiplicities)?,
        knots: reals(knots)?,
        weights: None,
    })
}

/// OCCT's reading of a knot sequence (`StepToGeom`): clamped or unclamped
/// when the multiplicities sum to poles plus degree plus one, periodic when
/// they sum to the poles plus one end's and the ends agree.
fn periodic(
    poles: usize,
    degree: usize,
    multiplicities: &[usize],
    invalid: &'static str,
) -> Named<bool> {
    let sum: usize = multiplicities.iter().sum();
    match (multiplicities.first(), multiplicities.last()) {
        _ if sum == poles + degree + 1 => Ok(false),
        (Some(a), Some(b)) if a == b && sum - a == poles => Ok(true),
        _ => Err(invalid),
    }
}

fn kernel_error(e: Error, invalid: &'static str) -> &'static str {
    match e {
        Error::LimitExceeded(_) => "BSplineControlDataLimit",
        _ => invalid,
    }
}

/// A 3D B-spline curve: OCCT's record and the kernel's curve.
pub(super) fn curve3(b: &Build, n: u64) -> Named<(BSplineRecord<3>, BSplineCurve3)> {
    let i = b.instance(n)?;
    let parts = curve_parts(i)?;
    let poles = parts
        .poles
        .iter()
        .map(|p| b.point(reference(p)?))
        .collect::<Named<Vec<[f64; 3]>>>()?;
    let periodic = periodic(
        poles.len(),
        parts.degree,
        &parts.multiplicities,
        "InvalidBSplineCurve",
    )?;
    let build = if periodic {
        BSplineCurve3::new_periodic
    } else {
        BSplineCurve3::new
    };
    let curve = build(
        parts.degree,
        poles
            .iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect(),
        parts.weights.clone(),
        parts.knots.clone(),
        parts.multiplicities.clone(),
    )
    .map_err(|e| kernel_error(e, "InvalidBSplineCurve"))?;
    Ok((
        BSplineRecord {
            degree: parts.degree,
            periodic,
            poles,
            weights: parts.weights,
            knots: parts.knots,
            multiplicities: parts.multiplicities,
        },
        curve,
    ))
}

/// A 2D B-spline curve in a surface's parameters (no length unit).
fn curve2(b: &Build, n: u64) -> Named<(BSplineRecord<2>, BSplineCurve2)> {
    let i = b.instance(n)?;
    let parts = curve_parts(i)?;
    let poles = parts
        .poles
        .iter()
        .map(|p| b.point2(reference(p)?))
        .collect::<Named<Vec<[f64; 2]>>>()?;
    let periodic = periodic(
        poles.len(),
        parts.degree,
        &parts.multiplicities,
        "InvalidBSplineCurve2d",
    )?;
    let build = if periodic {
        BSplineCurve2::new_periodic
    } else {
        BSplineCurve2::new
    };
    let curve = build(
        parts.degree,
        poles.iter().map(|p| Point2::new(p[0], p[1])).collect(),
        parts.weights.clone(),
        parts.knots.clone(),
        parts.multiplicities.clone(),
    )
    .map_err(|e| kernel_error(e, "InvalidBSplineCurve2d"))?;
    Ok((
        BSplineRecord {
            degree: parts.degree,
            periodic,
            poles,
            weights: parts.weights,
            knots: parts.knots,
            multiplicities: parts.multiplicities,
        },
        curve,
    ))
}

/// A B-spline surface: OCCT's record and the kernel's surface, the control
/// points rows in `u` (the kernel's `u`-major order).
pub(super) fn surface(b: &Build, n: u64) -> Named<(BSplineSurfaceRecord, BSplineSurface3)> {
    let i = b.instance(n)?;
    let with_knots = i
        .record("B_SPLINE_SURFACE_WITH_KNOTS")
        .ok_or("MalformedEntity")?;
    let (du, dv, grid, knots, weights) = if let Some(head) = i.record("B_SPLINE_SURFACE") {
        let [du, dv, grid, _, _, _, _] = head.parameters.as_slice() else {
            return Err("MalformedEntity");
        };
        let [um, vm, uk, vk, _] = with_knots.parameters.as_slice() else {
            return Err("MalformedEntity");
        };
        let weights = match i.record("RATIONAL_B_SPLINE_SURFACE") {
            Some(r) => match r.parameters.as_slice() {
                [w] => Some(list(w)?.iter().map(reals).collect::<Named<Vec<_>>>()?),
                _ => return Err("MalformedEntity"),
            },
            None => None,
        };
        (du, dv, grid, [um, vm, uk, vk], weights)
    } else {
        let [_, du, dv, grid, _, _, _, _, um, vm, uk, vk, _] = with_knots.parameters.as_slice()
        else {
            return Err("MalformedEntity");
        };
        (du, dv, grid, [um, vm, uk, vk], None)
    };
    let degrees = [count(du)?, count(dv)?];
    let rows = list(grid)?;
    let mut poles = Vec::new();
    let mut columns = None;
    for row in rows {
        let row = list(row)?;
        if *columns.get_or_insert(row.len()) != row.len() {
            return Err("MalformedEntity");
        }
        for p in row {
            poles.push(b.point(reference(p)?)?);
        }
    }
    let counts2 = [rows.len(), columns.unwrap_or(0)];
    let weights = match weights {
        Some(w) => {
            if w.len() != counts2[0] || w.iter().any(|row| row.len() != counts2[1]) {
                return Err("MalformedEntity");
            }
            Some(w.concat())
        }
        None => None,
    };
    let multiplicities = [counts(knots[0])?, counts(knots[1])?];
    let knots = [reals(knots[2])?, reals(knots[3])?];
    let periodic = [
        periodic(
            counts2[0],
            degrees[0],
            &multiplicities[0],
            "InvalidBSplineSurface",
        )?,
        periodic(
            counts2[1],
            degrees[1],
            &multiplicities[1],
            "InvalidBSplineSurface",
        )?,
    ];
    let axis = |k: usize| {
        let build = if periodic[k] {
            KnotVector::new_periodic
        } else {
            KnotVector::new
        };
        build(degrees[k], knots[k].clone(), multiplicities[k].clone())
            .map_err(|e| kernel_error(e, "InvalidBSplineSurface"))
    };
    let (u, v) = (axis(0)?, axis(1)?);
    let surface = BSplineSurface3::new(
        u,
        v,
        poles
            .iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect(),
        weights.clone(),
    )
    .map_err(|e| kernel_error(e, "InvalidBSplineSurface"))?;
    Ok((
        BSplineSurfaceRecord {
            degrees,
            periodic,
            counts: counts2,
            poles,
            weights,
            knots,
            multiplicities,
        },
        surface,
    ))
}

/// A pcurve the kernel takes from a file: a line or a B-spline in `(u, v)`.
pub(super) enum Pcurve {
    Line { p: [f64; 2], d: [f64; 2] },
    BSpline(Box<(BSplineRecord<2>, BSplineCurve2)>),
}

impl Pcurve {
    fn uv(&self, t: f64) -> Option<[f64; 2]> {
        match self {
            Self::Line { p, d } => Some([p[0] + t * d[0], p[1] + t * d[1]]),
            Self::BSpline(c) => c.1.point(t).ok().map(|p| [p.x, p.y]),
        }
    }
}

/// The 2D curve of a `PCURVE`, with its basis surface's entity number.
pub(super) fn pcurve(b: &Build, n: u64) -> Named<(u64, Pcurve)> {
    let p = b.params(n, "PCURVE", 3)?;
    let surface = reference(&p[1])?;
    let rep = b.params(reference(&p[2])?, "DEFINITIONAL_REPRESENTATION", 3)?;
    let item = reference(list(&rep[1])?.first().ok_or("MalformedEntity")?)?;
    let i = b.instance(item)?;
    if i.record("LINE").is_some() {
        let l = b.params(item, "LINE", 3)?;
        let origin = b.point2(reference(&l[1])?)?;
        let v = b.params(reference(&l[2])?, "VECTOR", 3)?;
        if real(&v[2])? <= 0.0 {
            return Err("DegenerateDirection");
        }
        let d = b.params(reference(&v[1])?, "DIRECTION", 2)?;
        let [x, y] = list(&d[1])? else {
            return Err("MalformedEntity");
        };
        let (x, y) = (real(x)?, real(y)?);
        let norm = x.hypot(y);
        if !(norm > 0.0 && norm.is_finite()) {
            return Err("DegenerateDirection");
        }
        return Ok((
            surface,
            Pcurve::Line {
                p: origin,
                d: [x / norm, y / norm],
            },
        ));
    }
    if is_curve(i) {
        return Ok((surface, Pcurve::BSpline(Box::new(curve2(b, item)?))));
    }
    Err("PCurveNotDerived")
}

fn distance2(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]) * (a[k] - b[k])).sum()
}

/// The parameter in `[lo, hi]` where `f` comes nearest `target`: an end
/// within `tol` of it exactly, else the best of samples spread over the
/// spans between `breaks`, refined by golden-section search between its
/// neighbours (OCCT's `ShapeAnalysis_Curve::Project` samples and refines
/// alike). Deterministic; a point that does not evaluate is infinitely far.
fn locate(
    f: &dyn Fn(f64) -> Option<[f64; 3]>,
    [lo, hi]: [f64; 2],
    breaks: &[f64],
    target: [f64; 3],
    tol: f64,
) -> f64 {
    let d = |t: f64| f(t).map_or(f64::INFINITY, |p| distance2(p, target));
    let (dlo, dhi) = (d(lo), d(hi));
    if dlo <= tol * tol && dlo <= dhi {
        return lo;
    }
    if dhi <= tol * tol {
        return hi;
    }
    let mut ends = vec![lo];
    ends.extend(breaks.iter().copied().filter(|k| *k > lo && *k < hi));
    ends.push(hi);
    let per = (2048 / ends.len()).clamp(4, 64);
    let mut ts = vec![lo];
    for w in ends.windows(2) {
        for k in 1..=per {
            ts.push(w[0] + (w[1] - w[0]) * k as f64 / per as f64);
        }
    }
    let values: Vec<f64> = ts.iter().map(|t| d(*t)).collect();
    let best = (0..ts.len()).fold(0, |m, k| if values[k] < values[m] { k } else { m });
    let (mut a, mut b) = (ts[best.saturating_sub(1)], ts[(best + 1).min(ts.len() - 1)]);
    const G: f64 = 0.618_033_988_749_894_9;
    let (mut c, mut e) = (b - G * (b - a), a + G * (b - a));
    let (mut dc, mut de) = (d(c), d(e));
    for _ in 0..100 {
        if dc < de {
            b = e;
            (e, de) = (c, dc);
            c = b - G * (b - a);
            dc = d(c);
        } else {
            a = c;
            (c, dc) = (e, de);
            e = a + G * (b - a);
            de = d(e);
        }
    }
    let mid = 0.5 * (a + b);
    if d(mid) <= values[best] {
        mid
    } else {
        ts[best]
    }
}

fn p3(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

/// An edge's range on its B-spline curve, running from `from` to `to`; a
/// closed edge (one vertex) runs the whole domain, or one period from its
/// vertex.
pub(super) fn curve_range(
    curve: &BSplineCurve3,
    from: [f64; 3],
    to: [f64; 3],
    closed: bool,
    tol: f64,
) -> Named<[f64; 2]> {
    let (a, b) = curve.domain();
    let f = |t: f64| curve.point(t).ok().map(p3);
    let breaks = curve.knots();
    if closed && !curve.is_periodic() {
        return Ok([a, b]);
    }
    let first = locate(&f, [a, b], breaks, from, tol);
    if closed {
        return Ok([first, first + (b - a)]);
    }
    let last = locate(&f, [a, b], breaks, to, tol);
    if curve.is_periodic() {
        return Ok([first, if last <= first { last + (b - a) } else { last }]);
    }
    if first >= last {
        return Err("EdgeRangeInverted");
    }
    Ok([first, last])
}

/// A pcurve's record and its range on a spline surface: where the image
/// `S(c(t))` meets the edge's vertices, running from `from` to `to`. A line
/// is taken over its part inside the surface's domain; a B-spline pcurve
/// must run along its parameter (`PCurveAgainstEdge` otherwise).
pub(super) fn pcurve_range(
    surface: &BSplineSurface3,
    pcurve: Pcurve,
    from: [f64; 3],
    to: [f64; 3],
    closed: bool,
    tol: f64,
) -> Named<(Record2, [f64; 2])> {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let image = |t: f64| {
        let [u, v] = pcurve.uv(t)?;
        // Rounding may put a boundary point a hair outside the domain.
        let (u, v) = (u.clamp(u0, u1), v.clamp(v0, v1));
        surface.point(u, v).ok().map(p3)
    };
    let (domain, breaks) = match &pcurve {
        Pcurve::Line { p, d } => {
            let mut range = [f64::NEG_INFINITY, f64::INFINITY];
            for (k, (lo, hi)) in [(u0, u1), (v0, v1)].into_iter().enumerate() {
                if d[k] == 0.0 {
                    if p[k] < lo || p[k] > hi {
                        return Err("PCurveOutsideSurface");
                    }
                    continue;
                }
                let (s, e) = ((lo - p[k]) / d[k], (hi - p[k]) / d[k]);
                range = [range[0].max(s.min(e)), range[1].min(s.max(e))];
            }
            if range[0] >= range[1] || !range.iter().all(|x| x.is_finite()) {
                return Err("PCurveOutsideSurface");
            }
            if closed {
                return Err("PCurveNotDerived");
            }
            (range, Vec::new())
        }
        Pcurve::BSpline(c) => {
            let (a, b) = c.1.as_curve3().domain();
            ([a, b], c.1.as_curve3().knots().to_vec())
        }
    };
    let range = if closed {
        domain
    } else {
        [
            locate(&image, domain, &breaks, from, tol),
            locate(&image, domain, &breaks, to, tol),
        ]
    };
    let record = match pcurve {
        Pcurve::Line { p, d } => {
            if range[0] == range[1] {
                return Err("ZeroLengthEdge");
            }
            Record2::Line { p, d }
        }
        Pcurve::BSpline(c) => {
            if range[0] >= range[1] {
                return Err("PCurveAgainstEdge");
            }
            Record2::BSpline(c.0)
        }
    };
    Ok((record, range))
}
