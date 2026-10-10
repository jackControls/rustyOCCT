//! STEP-b of the STEP import track in `REVIEW_NOTES.md`: B-spline curves and
//! surfaces (simple instances, and rational ones as the complex instances
//! OCCT and other writers produce), the file's pcurves on spline surfaces,
//! and where an edge's vertices lie on a spline or on a pcurve's image.
//!
//! Nothing here approximates the geometry: a located parameter, searched
//! for in binary64, is only a starting claim, which the validator certifies
//! (the vertices on their edges, each pcurve against its edge at matching
//! fractions), so a wrong location is an import failure with its issues,
//! never a moved edge.
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

/// A knot vector for binary64 de Boor: its flat knots (a periodic vector's
/// extended by a period at each end, as the kernel's) and its domain.
struct Knots {
    degree: usize,
    poles: usize,
    periodic: bool,
    domain: (f64, f64),
    flat: Vec<f64>,
}

impl Knots {
    fn new(k: &KnotVector) -> Self {
        let (degree, poles) = (k.degree(), k.pole_count());
        let mut flat: Vec<f64> = k
            .knots()
            .iter()
            .zip(k.multiplicities())
            .flat_map(|(&x, &m)| std::iter::repeat_n(x, m))
            .collect();
        let domain = k.domain();
        if k.is_periodic() {
            let period = domain.1 - domain.0;
            let first = k.multiplicities()[0];
            let extension = degree + 1 - first;
            let before = flat[poles - extension..poles].iter().map(|x| x - period);
            let after = flat[first..first + extension].iter().map(|x| x + period);
            flat = before.chain(flat.iter().copied()).chain(after).collect();
        }
        Self {
            degree,
            poles,
            periodic: k.is_periodic(),
            domain,
            flat,
        }
    }

    /// The parameter (a periodic one taken into the domain) and its span,
    /// the domain's end in the last span; none outside a nonperiodic domain.
    fn span(&self, u: f64) -> Option<(f64, usize)> {
        let (a, b) = self.domain;
        let u = if self.periodic {
            let u = a + (u - a).rem_euclid(b - a);
            if u < b {
                u
            } else {
                a
            }
        } else if (a..=b).contains(&u) {
            u
        } else {
            return None;
        };
        let at_end = !self.periodic && u == b;
        let span = self
            .flat
            .partition_point(|k| if at_end { *k < u } else { *k <= u })
            .checked_sub(1)?;
        (span >= self.degree && span + self.degree < self.flat.len()).then_some((u, span))
    }

    /// The index of the `j`th of the span's `degree + 1` poles.
    fn pole(&self, span: usize, j: usize) -> usize {
        (span - self.degree + j) % self.poles
    }

    /// De Boor's blend of the span's homogeneous poles `d` at `u`.
    fn blend(&self, u: f64, span: usize, d: &mut [[f64; 4]]) -> [f64; 4] {
        let p = self.degree;
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = span - p + j;
                let width = self.flat[i + p + 1 - r] - self.flat[i];
                let alpha = if width > 0.0 {
                    (u - self.flat[i]) / width
                } else {
                    0.0
                };
                let below = d[j - 1];
                for (x, b) in d[j].iter_mut().zip(below) {
                    *x = (1.0 - alpha) * b + alpha * *x;
                }
            }
        }
        d[p]
    }
}

fn homogeneous(p: Point3, w: f64) -> [f64; 4] {
    [p.x * w, p.y * w, p.z * w, w]
}

fn dehomogenized(h: [f64; 4]) -> Option<[f64; 3]> {
    let p = [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
    p.iter().all(|x| x.is_finite()).then_some(p)
}

/// A B-spline curve evaluated in binary64 (the location's search only).
struct FastCurve {
    knots: Knots,
    poles: Vec<[f64; 4]>,
}

impl FastCurve {
    fn new(c: &BSplineCurve3) -> Self {
        Self {
            knots: Knots::new(c.knot_vector()),
            poles: c
                .poles()
                .iter()
                .zip(c.weights())
                .map(|(p, w)| homogeneous(*p, *w))
                .collect(),
        }
    }

    fn point(&self, u: f64) -> Option<[f64; 3]> {
        let (u, span) = self.knots.span(u)?;
        let mut d: Vec<[f64; 4]> = (0..=self.knots.degree)
            .map(|j| self.poles[self.knots.pole(span, j)])
            .collect();
        dehomogenized(self.knots.blend(u, span, &mut d))
    }
}

/// A B-spline surface evaluated in binary64 (the location's search only).
struct FastSurface {
    u: Knots,
    v: Knots,
    poles: Vec<[f64; 4]>,
}

impl FastSurface {
    fn new(s: &BSplineSurface3) -> Self {
        Self {
            u: Knots::new(s.u_knots()),
            v: Knots::new(s.v_knots()),
            poles: s
                .poles()
                .iter()
                .zip(s.weights())
                .map(|(p, w)| homogeneous(*p, *w))
                .collect(),
        }
    }

    fn point(&self, u: f64, v: f64) -> Option<[f64; 3]> {
        let (u, su) = self.u.span(u)?;
        let (v, sv) = self.v.span(v)?;
        let mut row = vec![[0.0; 4]; self.v.degree + 1];
        let mut column: Vec<[f64; 4]> = (0..=self.u.degree)
            .map(|i| {
                let i = self.u.pole(su, i);
                for (j, x) in row.iter_mut().enumerate() {
                    *x = self.poles[i * self.v.poles + self.v.pole(sv, j)];
                }
                self.v.blend(v, sv, &mut row)
            })
            .collect();
        dehomogenized(self.u.blend(u, su, &mut column))
    }
}

/// The parameter in `[lo, hi]` where the curve comes nearest `target`: an
/// end within `tol` of it, its point correctly rounded (`exact`; on a
/// closed curve, where both are, the end for an edge's `last` vertex), else
/// the best of samples spread over the spans between `breaks`, refined by
/// golden-section search between its neighbours (OCCT's
/// `ShapeAnalysis_Curve::Project` samples and refines alike), its points
/// from binary64 de Boor (`fast`): the search's some two hundred points
/// taken exactly cost seconds on a trimmed spline face (the `step` fuzz
/// target's timeout `f692f018`), and its result is only the claim the
/// validator certifies. Deterministic (no libm); a point that does not
/// evaluate is infinitely far.
fn locate(
    exact: &dyn Fn(f64) -> Option<[f64; 3]>,
    fast: &dyn Fn(f64) -> Option<[f64; 3]>,
    [lo, hi]: [f64; 2],
    breaks: &[f64],
    target: [f64; 3],
    tol: f64,
    last: bool,
) -> f64 {
    let at_end = |t: f64| exact(t).is_some_and(|p| distance2(p, target) <= tol * tol);
    let d = |t: f64| fast(t).map_or(f64::INFINITY, |p| distance2(p, target));
    match (at_end(lo), at_end(hi)) {
        (true, true) => return if last { hi } else { lo },
        (true, false) => return lo,
        (false, true) => return hi,
        (false, false) => {}
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
    let search = FastCurve::new(curve);
    let fast = |t: f64| search.point(t);
    let breaks = curve.knots();
    if closed && !curve.is_periodic() {
        return Ok([a, b]);
    }
    let first = locate(&f, &fast, [a, b], breaks, from, tol, false);
    if closed {
        return Ok([first, first + (b - a)]);
    }
    let last = locate(&f, &fast, [a, b], breaks, to, tol, true);
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
    let search = FastSurface::new(surface);
    let trace = match &pcurve {
        Pcurve::BSpline(c) => Some(FastCurve::new(c.1.as_curve3())),
        Pcurve::Line { .. } => None,
    };
    let fast = |t: f64| {
        let [u, v] = match (&pcurve, &trace) {
            (Pcurve::Line { p, d }, _) => [p[0] + t * d[0], p[1] + t * d[1]],
            (_, Some(c)) => c.point(t).map(|p| [p[0], p[1]])?,
            _ => return None,
        };
        search.point(u.clamp(u0, u1), v.clamp(v0, v1))
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
            locate(&image, &fast, domain, &breaks, from, tol, false),
            locate(&image, &fast, domain, &breaks, to, tol, true),
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A closed clamped curve: its start and end are one point.
    fn loop_curve() -> BSplineCurve3 {
        let poles = [
            [1.0, 0.0],
            [1.0, 1.0],
            [-1.0, 1.0],
            [-1.0, -1.0],
            [1.0, -1.0],
            [1.0, 0.0],
        ];
        BSplineCurve3::new(
            2,
            poles.iter().map(|p| Point3::new(p[0], p[1], 0.0)).collect(),
            None,
            vec![0.0, 1.0, 2.0, 3.0, 4.0],
            vec![3, 1, 1, 1, 3],
        )
        .unwrap()
    }

    fn near(a: [f64; 3], b: [f64; 3]) -> bool {
        distance2(a, b).sqrt() <= 1e-13
    }

    /// The search's binary64 de Boor agrees with the kernel's correctly
    /// rounded points: clamped, unclamped, periodic (inside and outside
    /// its period) and rational, on curves and on a surface.
    #[test]
    fn binary64_search_points_agree_with_the_kernels() {
        let periodic = BSplineCurve3::new_periodic(
            3,
            (0..6)
                .map(|k| {
                    let a = k as f64;
                    Point3::new(a.cos(), a.sin(), 0.25 * a)
                })
                .collect(),
            Some(vec![1.0, 2.0, 0.5, 1.0, 3.0, 1.5]),
            vec![0.0, 0.5, 1.5, 2.0, 3.0, 3.25, 4.0],
            vec![1; 7],
        )
        .unwrap();
        let unclamped = BSplineCurve3::new(
            2,
            (0..5)
                .map(|k| Point3::new(k as f64, (k * k) as f64 * 0.5, 1.0))
                .collect(),
            None,
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            vec![1; 8],
        )
        .unwrap();
        for c in [loop_curve(), periodic, unclamped] {
            let fast = FastCurve::new(&c);
            let (a, b) = c.domain();
            for k in 0..=200 {
                let t = a + (b - a) * k as f64 / 200.0;
                assert!(near(fast.point(t).unwrap(), p3(c.point(t).unwrap())), "{t}");
                if c.is_periodic() {
                    let t = t + 2.0 * (b - a);
                    assert!(near(fast.point(t).unwrap(), p3(c.point(t).unwrap())), "{t}");
                }
            }
            if !c.is_periodic() {
                assert_eq!(fast.point(b + 1e-9), None);
            }
        }
        let poles: Vec<Point3> = (0..20)
            .map(|k| {
                let (i, j) = ((k / 5) as f64, (k % 5) as f64);
                Point3::new(3.0 * i, 2.0 * j, (i * j).sin())
            })
            .collect();
        let weights = (0..20).map(|k| 1.0 + 0.125 * (k % 3) as f64).collect();
        let s = BSplineSurface3::new(
            KnotVector::new(3, vec![0.0, 1.0], vec![4, 4]).unwrap(),
            KnotVector::new(2, vec![0.0, 0.25, 0.5, 1.0], vec![3, 1, 1, 3]).unwrap(),
            poles,
            Some(weights),
        )
        .unwrap();
        let fast = FastSurface::new(&s);
        for i in 0..=20 {
            for j in 0..=20 {
                let (u, v) = (i as f64 / 20.0, j as f64 / 20.0);
                assert!(
                    near(fast.point(u, v).unwrap(), p3(s.point(u, v).unwrap())),
                    "{u} {v}"
                );
            }
        }
    }

    /// Locating a vertex takes exact points only to test the ends: the
    /// search (some two hundred points) is binary64 (the `step` target's
    /// timeout `f692f018`, a trimmed face's vertex moved off its curves).
    #[test]
    fn a_located_vertex_is_searched_in_binary64() {
        let c = loop_curve();
        let fast = FastCurve::new(&c);
        let calls = std::cell::Cell::new(0);
        let exact = |t: f64| {
            calls.set(calls.get() + 1);
            c.point(t).ok().map(p3)
        };
        let (a, b) = c.domain();
        let target = p3(c.point(2.5).unwrap());
        let t = locate(
            &exact,
            &|t| fast.point(t),
            [a, b],
            c.knots(),
            target,
            1e-7,
            false,
        );
        assert!((t - 2.5).abs() < 1e-9, "{t}");
        assert_eq!(calls.get(), 2);
        // Off the curve: the nearest point, still two exact points.
        let t = locate(
            &exact,
            &|t| fast.point(t),
            [a, b],
            c.knots(),
            [0.0, 2.0, 0.0],
            1e-7,
            false,
        );
        assert!((t - 1.0).abs() < 1e-6, "{t}");
        assert_eq!(calls.get(), 4);
    }

    #[test]
    fn ranges_run_from_the_first_vertex_to_the_last() {
        let c = loop_curve();
        let start = [1.0, 0.0, 0.0];
        // Two vertices at the closed curve's joint: the whole domain.
        assert_eq!(curve_range(&c, start, start, false, 1e-7), Ok([0.0, 4.0]));
        assert_eq!(curve_range(&c, start, start, true, 1e-7), Ok([0.0, 4.0]));
        // An interior vertex, found to rounding.
        let p = c.point(1.5).unwrap();
        let [first, last] = curve_range(&c, start, p3(p), false, 1e-7).unwrap();
        assert_eq!(first, 0.0);
        assert!((last - 1.5).abs() < 1e-9, "{last}");
        // Against the curve's parameter.
        assert_eq!(
            curve_range(&c, p3(p), [-1.0, -0.5, 0.0], false, 1e-7).map(|_| ()),
            Ok(())
        );
        assert_eq!(
            curve_range(&c, p3(c.point(3.0).unwrap()), p3(p), false, 1e-7),
            Err("EdgeRangeInverted")
        );
    }
}
