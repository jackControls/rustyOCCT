use crate::decide;
use crate::identity::InputLabel;
use crate::math::{finite, sum};
use crate::predicates::{orient2d_finite, Orientation2};
use crate::{Error, Point2, Result, Tolerance};
use std::f64::consts::PI;

const MAX_EDGES: usize = 4096;
const MAX_HOLES: usize = 128;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Location {
    Inside,
    Boundary,
    Outside,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct AreaMoments {
    pub area: f64,
    pub centroid: Point2,
    /// Central integrals of x*x, x*y and y*y over the material region.
    pub second: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum BoundaryKind {
    Polygon(Vec<Point2>),
    Circle {
        center: Point2,
        radius: f64,
    },
    /// Points with a segment from each to the next (S5), at least one an
    /// arc.
    Path {
        points: Vec<Point2>,
        segments: Vec<Segment>,
    },
}

/// A segment of a path boundary from one point to the next (S5).
#[derive(Debug, Clone, PartialEq)]
pub enum Segment {
    Line,
    /// A circular arc about `center` of the given radius, turning
    /// counter-clockwise (`ccw`) or clockwise about the profile's normal.
    /// Both of its points lie within tolerance of the circle.
    Arc {
        center: Point2,
        radius: f64,
        ccw: bool,
    },
    /// A planar nonrational B-spline (S8b) over its whole domain, from its
    /// first pole to its last (or the reverse when flagged): its ends are the
    /// two points exactly.
    Spline(crate::topology::SplineSpan<crate::BSplineCurve2>),
}

impl Segment {
    /// The same segment traversed backwards.
    fn reversed(&self) -> Self {
        match self.clone() {
            Segment::Line => Segment::Line,
            Segment::Arc {
                center,
                radius,
                ccw,
            } => Segment::Arc {
                center,
                radius,
                ccw: !ccw,
            },
            Segment::Spline(span) => Segment::Spline(span.reversed()),
        }
    }
}

/// A spline segment from `a` to `b` checked (S8b): nonrational, degree 1 to
/// 7, over its whole domain, its ends the two points exactly, C1 inside and
/// turning through less than a quarter-turn; its screen parts.
fn spline_parts(
    span: &crate::topology::SplineSpan<crate::BSplineCurve2>,
    a: Point2,
    b: Point2,
) -> Result<Vec<decide::splines::Part>> {
    let curve = span.curve().as_curve3();
    if curve.is_rational() {
        return Err(Error::OutOfDomain("a rational spline profile segment"));
    }
    if !(1..=7).contains(&curve.degree()) {
        return Err(Error::OutOfDomain(
            "a spline profile segment of degree above 7",
        ));
    }
    let (first, last) = curve.domain();
    if span.range() != [first, last] || curve.is_periodic() {
        return Err(Error::InvalidCurve(
            "a spline profile segment over part of its domain",
        ));
    }
    let poles = span.curve().poles();
    let (s, e) = (poles[0], poles[poles.len() - 1]);
    let (s, e) = if span.is_reversed() { (e, s) } else { (s, e) };
    if s != a || e != b {
        return Err(Error::InvalidCurve(
            "a spline segment's ends off its path points",
        ));
    }
    if !crate::topology::validate::continuity::curve_c1(curve, [first, last], false) {
        return Err(Error::InvalidCurve(
            "a spline profile segment not C1 inside",
        ));
    }
    let parts = decide::splines::spline(span).ok_or(Error::InvalidCurve(
        "a spline profile segment's Bézier arcs",
    ))?;
    if !decide::splines::quarter_turn(&parts) {
        return Err(Error::OutOfDomain(
            "a spline profile segment turning through a quarter-turn or more (split it)",
        ));
    }
    Ok(parts)
}

/// Every segment of a path as screen parts.
fn path_parts(points: &[Point2], segments: &[Segment]) -> Result<Vec<Vec<decide::splines::Part>>> {
    let n = points.len();
    segments
        .iter()
        .enumerate()
        .map(|(i, segment)| {
            let (a, b) = (points[i], points[(i + 1) % n]);
            Ok(match segment {
                Segment::Line => decide::splines::line(a, b),
                Segment::Arc {
                    center,
                    radius,
                    ccw,
                } => decide::splines::arc(&decide::arcs::Arc2 {
                    center: *center,
                    radius: *radius,
                    start: a,
                    end: b,
                    ccw: *ccw,
                }),
                Segment::Spline(span) => spline_parts(span, a, b)?,
            })
        })
        .collect()
}

/// The signed sweep of an arc from `a` to `b` about `center`: in `(0, 2π)`
/// turning counter-clockwise, in `(-2π, 0)` clockwise.
pub(crate) fn arc_sweep(center: Point2, a: Point2, b: Point2, ccw: bool) -> f64 {
    let (u, v) = (
        Point2::new(a.x - center.x, a.y - center.y),
        Point2::new(b.x - center.x, b.y - center.y),
    );
    let mut turn = cross(u, v).atan2(u.x * v.x + u.y * v.y);
    if turn <= 0.0 {
        turn += 2.0 * PI;
    }
    if ccw {
        turn
    } else {
        turn - 2.0 * PI
    }
}

/// `∫ cos^m t sin^n t dt` over `[a, b]` for `m + n <= 4`, by the reduction
/// formulas.
fn trig_integral(m: u32, n: u32, a: f64, b: f64) -> f64 {
    let at = |t: f64, m: i32, n: i32| t.cos().powi(m) * t.sin().powi(n);
    match (m, n) {
        (0, 0) => b - a,
        (1, 0) => b.sin() - a.sin(),
        (0, 1) => a.cos() - b.cos(),
        (1, 1) => 0.5 * (b.sin().powi(2) - a.sin().powi(2)),
        (m, n) if n >= 2 => {
            let k = f64::from(m + n);
            -(at(b, m as i32 + 1, n as i32 - 1) - at(a, m as i32 + 1, n as i32 - 1)) / k
                + f64::from(n - 1) / k * trig_integral(m, n - 2, a, b)
        }
        (m, n) => {
            let k = f64::from(m + n);
            (at(b, m as i32 - 1, n as i32 + 1) - at(a, m as i32 - 1, n as i32 + 1)) / k
                + f64::from(m - 1) / k * trig_integral(m - 2, n, a, b)
        }
    }
}

/// Green's-theorem integrals `[A, ∫x, ∫y, ∫x², ∫xy, ∫y²]` of the region
/// left of an arc (centre `(cx, cy)` relative to the anchor) from angle
/// `start` over `sweep`.
fn arc_moments(cx: f64, cy: f64, r: f64, start: f64, sweep: f64) -> [f64; 6] {
    let (a, b) = (start, start + sweep);
    let i = |m: u32, n: u32| trig_integral(m, n, a, b);
    // x = cx + r c, y = cy + r s, dx = -r s dt, dy = r c dt.
    let area = 0.5 * (r * cx * i(1, 0) + r * cy * i(0, 1) + r * r * i(0, 0));
    // ∫x dA = ∮ x²/2 dy.
    let mx = 0.5 * r * (cx * cx * i(1, 0) + 2.0 * cx * r * i(2, 0) + r * r * i(3, 0));
    // ∫y dA = -∮ y²/2 dx.
    let my = 0.5 * r * (cy * cy * i(0, 1) + 2.0 * cy * r * i(0, 2) + r * r * i(0, 3));
    // ∫x² dA = ∮ x³/3 dy.
    let xx = r / 3.0
        * (cx.powi(3) * i(1, 0)
            + 3.0 * cx * cx * r * i(2, 0)
            + 3.0 * cx * r * r * i(3, 0)
            + r.powi(3) * i(4, 0));
    // ∫y² dA = -∮ y³/3 dx.
    let yy = r / 3.0
        * (cy.powi(3) * i(0, 1)
            + 3.0 * cy * cy * r * i(0, 2)
            + 3.0 * cy * r * r * i(0, 3)
            + r.powi(3) * i(0, 4));
    // ∫xy dA = ∮ x² y/2 dy.
    let xy = 0.5
        * r
        * (cx * cx * cy * i(1, 0)
            + cx * cx * r * i(1, 1)
            + 2.0 * cx * r * cy * i(2, 0)
            + 2.0 * cx * r * r * i(2, 1)
            + r * r * cy * i(3, 0)
            + r.powi(3) * i(3, 1));
    [area, mx, my, xx, xy, yy]
}

/// Caller labels for one boundary, the roots of its entities' ids. A polygon
/// has one segment and one vertex label per distinct point (a repeated closing
/// point has none); a circle has one of each (its seam vertex).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryLabels {
    pub boundary: InputLabel,
    pub segments: Vec<InputLabel>,
    pub vertices: Vec<InputLabel>,
}

/// A validated simple closed boundary. Polygon points are stored CCW with the
/// original first vertex retained. A repeated closing point is optional.
/// Equality compares stored geometry and stored-order labels, not the
/// caller's input orientation.
#[derive(Debug, Clone)]
pub struct Boundary {
    pub(crate) kind: BoundaryKind,
    pub(crate) moments: AreaMoments,
    perimeter: f64,
    /// The caller gave the polygon clockwise; storage reversed points[1..].
    reversed: bool,
    /// In stored order.
    labels: Option<BoundaryLabels>,
}

impl PartialEq for Boundary {
    fn eq(&self, other: &Self) -> bool {
        (&self.kind, &self.moments, self.perimeter, &self.labels)
            == (&other.kind, &other.moments, other.perimeter, &other.labels)
    }
}

impl Boundary {
    pub fn polygon(mut points: Vec<Point2>, tolerance: Tolerance) -> Result<Self> {
        if points.len() > MAX_EDGES + 1 {
            return Err(Error::LimitExceeded("polygon vertices"));
        }
        for point in &points {
            point.checked(tolerance)?;
        }
        let tol = tolerance.linear();
        if points.len() > 1 && decide::distance_le(points[0], *points.last().unwrap(), &[tol]) {
            points.pop();
        }
        if points.len() < 3 {
            return Err(Error::Degenerate("polygon"));
        }
        if points.len() > MAX_EDGES {
            return Err(Error::LimitExceeded("polygon vertices"));
        }
        let count = points.len();
        for i in 0..count {
            let (a, b, c) = (points[i], points[(i + 1) % count], points[(i + 2) % count]);
            if decide::distance_le(a, b, &[tol]) {
                return Err(Error::Degenerate("polygon edge"));
            }
            // Adjacent edges may continue along a line, but may not double back.
            if decide::segment_distance_le(a, b, c, &[tol])
                || decide::segment_distance_le(c, a, b, &[tol])
            {
                return Err(Error::SelfIntersection);
            }
            for j in (i + 1)..count {
                if j == i + 1 || (i == 0 && j == count - 1) {
                    continue;
                }
                if segments_touch(a, b, points[j], points[(j + 1) % count], tol) {
                    return Err(Error::SelfIntersection);
                }
            }
        }
        let perimeter = finite(sum(edges(&points).map(|(a, b)| a.distance(b))), "perimeter")?;
        let anchor = points[0];
        let signed_area = sum(edges(&points).map(|(a, b)| orient(anchor, a, b))) * 0.5;
        if !signed_area.is_finite()
            || decide::area_is_degenerate(&[decide::Outline::Polygon(&points)], tol)
        {
            return Err(Error::Degenerate("polygon area"));
        }
        let reversed = signed_area < 0.0;
        if reversed {
            points[1..].reverse();
        }
        let area = signed_area.abs();
        // Integrate in a local frame to avoid subtracting large global moments.
        let local: Vec<_> = points
            .iter()
            .map(|p| Point2::new(p.x - anchor.x, p.y - anchor.y))
            .collect();
        let cx = sum(edges(&local).map(|(a, b)| (a.x + b.x) * cross(a, b))) / (6.0 * area);
        let cy = sum(edges(&local).map(|(a, b)| (a.y + b.y) * cross(a, b))) / (6.0 * area);
        let centered: Vec<_> = local
            .iter()
            .map(|p| Point2::new(p.x - cx, p.y - cy))
            .collect();
        let xx =
            sum(edges(&centered).map(|(a, b)| (a.x * a.x + a.x * b.x + b.x * b.x) * cross(a, b)))
                / 12.0;
        let yy =
            sum(edges(&centered).map(|(a, b)| (a.y * a.y + a.y * b.y + b.y * b.y) * cross(a, b)))
                / 12.0;
        let xy = sum(edges(&centered).map(|(a, b)| {
            (2.0 * a.x * a.y + a.x * b.y + b.x * a.y + 2.0 * b.x * b.y) * cross(a, b)
        })) / 24.0;
        for value in [xx, xy, yy] {
            finite(value, "area moment")?;
        }
        Ok(Self {
            kind: BoundaryKind::Polygon(points),
            moments: AreaMoments {
                area,
                centroid: Point2::new(anchor.x + cx, anchor.y + cy),
                second: [xx, xy, yy],
            },
            perimeter,
            reversed,
            labels: None,
        })
    }

    /// A path of line and circular-arc segments (S5): segment `i` runs from
    /// point `i` to point `i + 1`, the last back to the first. Each arc's
    /// points lie within tolerance of its circle; its sweep is the turn from
    /// the first point's direction to the second's, in its direction. The
    /// path is stored counter-clockwise like a polygon. A path of lines only
    /// is exactly [`Boundary::polygon`] of its points.
    pub fn path(points: Vec<Point2>, segments: Vec<Segment>, tolerance: Tolerance) -> Result<Self> {
        if points.len() != segments.len() {
            return Err(Error::InvalidCurve("one segment per path point"));
        }
        if segments.iter().all(|s| *s == Segment::Line) {
            return Self::polygon(points, tolerance);
        }
        if points.len() > MAX_EDGES {
            return Err(Error::LimitExceeded("path segments"));
        }
        if points.len() < 2 {
            return Err(Error::Degenerate("path"));
        }
        for point in &points {
            point.checked(tolerance)?;
        }
        let tol = tolerance.linear();
        let count = points.len();
        let mut perimeter = 0.0;
        for (i, segment) in segments.iter().enumerate() {
            let (a, b) = (points[i], points[(i + 1) % count]);
            if decide::distance_le(a, b, &[tol]) {
                return Err(Error::Degenerate("path segment"));
            }
            perimeter += match segment {
                Segment::Line => a.distance(b),
                Segment::Spline(span) => decide::splines::length(&spline_parts(span, a, b)?),
                Segment::Arc {
                    center,
                    radius,
                    ccw,
                } => {
                    center.checked(tolerance)?;
                    tolerance.resolve(&[*radius, center.x + radius, center.y + radius])?;
                    if !radius.is_finite() || *radius <= tol {
                        return Err(Error::Degenerate("arc radius"));
                    }
                    for end in [a, b] {
                        if !(decide::distance_le(end, *center, &[*radius, tol])
                            && decide::distance_ge(end, *center, &[*radius, -tol]))
                        {
                            return Err(Error::InvalidCurve("arc point off its circle"));
                        }
                    }
                    radius * arc_sweep(*center, a, b, *ccw).abs()
                }
            };
        }
        let perimeter = finite(perimeter, "perimeter")?;
        // S8b: a path with a spline segment by the spline screen.
        if segments.iter().any(|s| matches!(s, Segment::Spline(_))) {
            let parts = path_parts(&points, &segments)?;
            for i in 0..count {
                let j = (i + 1) % count;
                if !decide::splines::adjacent_valid(&parts[i], &parts[j], tol, count == 2) {
                    return Err(Error::SelfIntersection);
                }
                for k in (i + 2)..count {
                    if (k + 1) % count == i {
                        continue;
                    }
                    if !decide::splines::apart(&parts[i], &parts[k], tol) {
                        return Err(Error::SelfIntersection);
                    }
                }
                if count == 2 {
                    break;
                }
            }
        }
        let pieces = path_pieces(&points, &segments);
        let splined = segments.iter().any(|s| matches!(s, Segment::Spline(_)));
        for i in (0..count).filter(|_| !splined) {
            let j = (i + 1) % count;
            // Two pieces meet at both ends of a two-segment path.
            let others: Vec<Point2> = if count == 2 {
                vec![points[i]]
            } else {
                Vec::new()
            };
            if decide::arcs::adjacent_invalid(&pieces[i], &pieces[j], tol, &others) {
                return Err(Error::SelfIntersection);
            }
            for k in (i + 2)..count {
                if (k + 1) % count == i {
                    continue;
                }
                if decide::arcs::pieces_within(&pieces[i], &pieces[k], tol) {
                    return Err(Error::SelfIntersection);
                }
            }
        }
        if decide::area_is_degenerate(&[decide::Outline::Path(&points, &segments)], tol) {
            return Err(Error::Degenerate("path area"));
        }
        let (moments, signed) = path_moments(&points, &segments)?;
        let reversed = signed < 0.0;
        let (points, segments) = if reversed {
            let mut p = points.clone();
            p[1..].reverse();
            let s = (0..count)
                .map(|j| segments[count - 1 - j].reversed())
                .collect();
            (p, s)
        } else {
            (points, segments)
        };
        Ok(Self {
            kind: BoundaryKind::Path { points, segments },
            moments,
            perimeter,
            reversed,
            labels: None,
        })
    }

    /// Rectangle from (0, 0) to (width, height).
    pub fn rectangle(width: f64, height: f64, tolerance: Tolerance) -> Result<Self> {
        finite(width, "rectangle width")?;
        finite(height, "rectangle height")?;
        if width <= tolerance.linear() || height <= tolerance.linear() {
            return Err(Error::Degenerate("rectangle"));
        }
        Self::polygon(
            vec![
                Point2::new(0.0, 0.0),
                Point2::new(width, 0.0),
                Point2::new(width, height),
                Point2::new(0.0, height),
            ],
            tolerance,
        )
    }

    pub fn circle(center: Point2, radius: f64, tolerance: Tolerance) -> Result<Self> {
        center.checked(tolerance)?;
        tolerance.resolve(&[
            radius,
            center.x + radius,
            center.y + radius,
            center.x - radius,
            center.y - radius,
        ])?;
        if radius <= tolerance.linear() {
            return Err(Error::Degenerate("circle radius"));
        }
        let area = finite(PI * radius * radius, "circle area")?;
        let moment = finite(area * radius * radius / 4.0, "circle moment")?;
        if area <= 0.0 || moment <= 0.0 {
            return Err(Error::Degenerate("circle moments"));
        }
        Ok(Self {
            kind: BoundaryKind::Circle { center, radius },
            moments: AreaMoments {
                area,
                centroid: center,
                second: [moment, 0.0, moment],
            },
            perimeter: 2.0 * PI * radius,
            reversed: false,
            labels: None,
        })
    }

    /// Attach caller labels, given in the caller's input order: segment `i`
    /// runs from input vertex `i` to `i + 1`. They are stored in the
    /// boundary's counter-clockwise order. Labels must be distinct.
    pub fn with_labels(mut self, labels: BoundaryLabels) -> Result<Self> {
        let count = self.segment_count();
        if labels.segments.len() != count || labels.vertices.len() != count {
            return Err(Error::InvalidLabel(
                "one segment and one vertex label per point",
            ));
        }
        let mut seen = std::collections::BTreeSet::new();
        let all = std::iter::once(&labels.boundary)
            .chain(&labels.segments)
            .chain(&labels.vertices);
        if !all.into_iter().all(|l| seen.insert(*l)) {
            return Err(Error::InvalidLabel("duplicate label"));
        }
        let (segments, vertices) = if self.reversed {
            // Stored segment j is input segment n-1-j traversed backwards;
            // stored vertex j is input vertex (n-j) mod n.
            (
                (0..count).map(|j| labels.segments[count - 1 - j]).collect(),
                (0..count)
                    .map(|j| labels.vertices[(count - j) % count])
                    .collect(),
            )
        } else {
            (labels.segments, labels.vertices)
        };
        self.labels = Some(BoundaryLabels {
            boundary: labels.boundary,
            segments,
            vertices,
        });
        Ok(self)
    }
    /// Labels in stored (counter-clockwise) order, if any.
    pub fn labels(&self) -> Option<&BoundaryLabels> {
        self.labels.as_ref()
    }
    pub fn area(&self) -> f64 {
        self.moments.area
    }
    pub fn centroid(&self) -> Point2 {
        self.moments.centroid
    }
    pub fn perimeter(&self) -> f64 {
        self.perimeter
    }
    pub fn polygon_vertices(&self) -> Option<&[Point2]> {
        match &self.kind {
            BoundaryKind::Polygon(points) => Some(points),
            _ => None,
        }
    }
    /// The stored (counter-clockwise) points and segments of a path.
    pub fn path_geometry(&self) -> Option<(&[Point2], &[Segment])> {
        match &self.kind {
            BoundaryKind::Path { points, segments } => Some((points, segments)),
            _ => None,
        }
    }
    /// The number of segments (and of points): one for a circle.
    pub(crate) fn segment_count(&self) -> usize {
        match &self.kind {
            BoundaryKind::Polygon(points) => points.len(),
            BoundaryKind::Circle { .. } => 1,
            BoundaryKind::Path { points, .. } => points.len(),
        }
    }
    /// The boundary as pieces for the proximity screen.
    pub(crate) fn pieces(&self) -> Vec<decide::arcs::Piece> {
        match &self.kind {
            BoundaryKind::Polygon(points) => edges(points)
                .map(|(a, b)| decide::arcs::Piece::Line(a, b))
                .collect(),
            BoundaryKind::Circle { center, radius } => {
                vec![decide::arcs::Piece::Circle(*center, *radius)]
            }
            BoundaryKind::Path { points, segments } => path_pieces(points, segments),
        }
    }
    /// Whether the boundary has a spline segment (S8b).
    fn has_spline(&self) -> bool {
        matches!(&self.kind, BoundaryKind::Path { segments, .. }
            if segments.iter().any(|s| matches!(s, Segment::Spline(_))))
    }

    /// Each segment as screen parts (a polygon's edges, a circle's whole
    /// turn); `None` when a spline segment is not valid.
    fn spline_chains(&self) -> Option<Vec<Vec<decide::splines::Part>>> {
        match &self.kind {
            BoundaryKind::Polygon(points) => Some(
                edges(points)
                    .map(|(a, b)| decide::splines::line(a, b))
                    .collect(),
            ),
            BoundaryKind::Circle { center, radius } => {
                Some(vec![decide::splines::circle(*center, *radius)])
            }
            BoundaryKind::Path { points, segments } => path_parts(points, segments).ok(),
        }
    }

    pub fn circle_geometry(&self) -> Option<(Point2, f64)> {
        match self.kind {
            BoundaryKind::Circle { center, radius } => Some((center, radius)),
            _ => None,
        }
    }

    pub(crate) fn validated(&self, tolerance: Tolerance) -> Result<Self> {
        let mut rebuilt = match &self.kind {
            BoundaryKind::Polygon(points) => Self::polygon(points.clone(), tolerance)?,
            BoundaryKind::Circle { center, radius } => Self::circle(*center, *radius, tolerance)?,
            BoundaryKind::Path { points, segments } => {
                Self::path(points.clone(), segments.clone(), tolerance)?
            }
        };
        // Stored points are already counter-clockwise: keep the caller's
        // orientation record and the stored-order labels.
        if rebuilt.kind != self.kind {
            return Err(Error::Degenerate("polygon"));
        }
        rebuilt.reversed = self.reversed;
        rebuilt.labels = self.labels.clone();
        Ok(rebuilt)
    }
    fn sample(&self) -> Point2 {
        match &self.kind {
            BoundaryKind::Polygon(points) => points[0],
            BoundaryKind::Circle { center, radius } => Point2::new(center.x + radius, center.y),
            BoundaryKind::Path { points, .. } => points[0],
        }
    }
    pub(crate) fn locate(&self, point: Point2, tolerance: Tolerance) -> Location {
        match &self.kind {
            BoundaryKind::Circle { center, radius } => {
                let tol = tolerance.linear();
                if decide::distance_le(point, *center, &[*radius, tol])
                    && decide::distance_ge(point, *center, &[*radius, -tol])
                {
                    Location::Boundary
                } else if !decide::distance_ge(point, *center, &[*radius]) {
                    Location::Inside
                } else {
                    Location::Outside
                }
            }
            // S8b: a path with a spline segment by its screen parts.
            BoundaryKind::Path { points, segments }
                if segments.iter().any(|s| matches!(s, Segment::Spline(_))) =>
            {
                let tol = tolerance.linear();
                let Ok(chains) = path_parts(points, segments) else {
                    return Location::Boundary;
                };
                if chains
                    .iter()
                    .any(|parts| decide::splines::point_within(point, parts, tol))
                {
                    return Location::Boundary;
                }
                let mut crossings = 0;
                for (segment, (i, parts)) in segments.iter().zip(chains.iter().enumerate()) {
                    crossings += match segment {
                        Segment::Spline(_) => match decide::splines::ray_crossings(point, parts) {
                            Some(n) => n,
                            None => return Location::Boundary,
                        },
                        _ => decide::arcs::ray_crossings(point, &path_pieces(points, segments)[i]),
                    };
                }
                if crossings % 2 == 1 {
                    Location::Inside
                } else {
                    Location::Outside
                }
            }
            BoundaryKind::Path { points, segments } => {
                let tol = tolerance.linear();
                let pieces = path_pieces(points, segments);
                if pieces.iter().any(|piece| match piece {
                    decide::arcs::Piece::Line(a, b) => {
                        decide::segment_distance_le(point, *a, *b, &[tol])
                    }
                    decide::arcs::Piece::Arc(arc) => {
                        decide::arcs::point_arc_within(point, arc, tol)
                    }
                    decide::arcs::Piece::Circle(..) => unreachable!("a path has no circle"),
                }) {
                    return Location::Boundary;
                }
                let crossings: u32 = pieces
                    .iter()
                    .map(|piece| decide::arcs::ray_crossings(point, piece))
                    .sum();
                if crossings % 2 == 1 {
                    Location::Inside
                } else {
                    Location::Outside
                }
            }
            BoundaryKind::Polygon(points) => {
                let mut inside = false;
                for (a, b) in edges(points) {
                    if decide::segment_distance_le(point, a, b, &[tolerance.linear()]) {
                        return Location::Boundary;
                    }
                    if (a.y > point.y) != (b.y > point.y) {
                        // Decide which side of the ray crossing the point occupies
                        // without constructing a rounded intersection coordinate.
                        let side = orient2d_finite(a, b, point);
                        if side == Orientation2::Collinear {
                            return Location::Boundary;
                        }
                        if (b.y > a.y && side == Orientation2::CounterClockwise)
                            || (b.y < a.y && side == Orientation2::Clockwise)
                        {
                            inside = !inside;
                        }
                    }
                }
                if inside {
                    Location::Inside
                } else {
                    Location::Outside
                }
            }
        }
    }
}

/// One connected material region. Holes are strictly contained, mutually
/// disjoint and non-nested. Model nested material islands as separate profiles.
#[derive(Debug, Clone, PartialEq)]
pub struct Profile {
    outer: Boundary,
    holes: Vec<Boundary>,
    tolerance: Tolerance,
    pub(crate) moments: AreaMoments,
    perimeter: f64,
}

impl Profile {
    /// Whether any boundary has a spline segment (S8b).
    pub(crate) fn has_spline(&self) -> bool {
        std::iter::once(&self.outer)
            .chain(self.holes.iter())
            .any(|b| match &b.kind {
                BoundaryKind::Path { segments, .. } => {
                    segments.iter().any(|s| matches!(s, Segment::Spline(_)))
                }
                _ => false,
            })
    }
    pub fn new(outer: Boundary, holes: Vec<Boundary>, tolerance: Tolerance) -> Result<Self> {
        if holes.len() > MAX_HOLES {
            return Err(Error::LimitExceeded("profile holes"));
        }
        // Enforce the aggregate limit before repeating quadratic validation.
        let total_edges: usize = std::iter::once(&outer)
            .chain(&holes)
            .map(Boundary::segment_count)
            .sum();
        if total_edges > MAX_EDGES {
            return Err(Error::LimitExceeded("profile edges"));
        }
        let outer = outer.validated(tolerance)?;
        let holes = holes
            .iter()
            .map(|hole| hole.validated(tolerance))
            .collect::<Result<Vec<_>>>()?;
        let boundaries = std::iter::once(&outer).chain(&holes).collect::<Vec<_>>();
        let mut labels = std::collections::BTreeSet::new();
        for b in boundaries.iter().filter_map(|b| b.labels()) {
            let all = std::iter::once(&b.boundary)
                .chain(&b.segments)
                .chain(&b.vertices);
            if !all.into_iter().all(|l| labels.insert(*l)) {
                return Err(Error::InvalidLabel("duplicate label in profile"));
            }
        }
        for (i, hole) in holes.iter().enumerate() {
            if boundaries_touch(&outer, hole, tolerance.linear())
                || outer.locate(hole.sample(), tolerance) != Location::Inside
            {
                return Err(Error::InvalidHole(i));
            }
            for (j, other) in holes.iter().enumerate().take(i) {
                if boundaries_touch(hole, other, tolerance.linear())
                    || hole.locate(other.sample(), tolerance) != Location::Outside
                    || other.locate(hole.sample(), tolerance) != Location::Outside
                {
                    return Err(Error::IntersectingBoundaries(j + 1, i + 1));
                }
            }
        }
        let signed = |i| if i == 0 { 1.0 } else { -1.0 };
        let area = sum(boundaries
            .iter()
            .enumerate()
            .map(|(i, b)| signed(i) * b.area()));
        let perimeter = sum(boundaries.iter().map(|b| b.perimeter()));
        let outlines: Vec<decide::Outline> = boundaries
            .iter()
            .map(|b| match &b.kind {
                BoundaryKind::Polygon(points) => decide::Outline::Polygon(points),
                BoundaryKind::Circle { radius, .. } => decide::Outline::Circle(*radius),
                BoundaryKind::Path { points, segments } => decide::Outline::Path(points, segments),
            })
            .collect();
        if decide::area_is_degenerate(&outlines, tolerance.linear()) {
            return Err(Error::Degenerate("profile material area"));
        }
        let anchor = outer.centroid();
        let cx = sum(boundaries
            .iter()
            .enumerate()
            .map(|(i, b)| signed(i) * b.area() * (b.centroid().x - anchor.x)))
            / area;
        let cy = sum(boundaries
            .iter()
            .enumerate()
            .map(|(i, b)| signed(i) * b.area() * (b.centroid().y - anchor.y)))
            / area;
        let centroid = Point2::new(anchor.x + cx, anchor.y + cy).checked(tolerance)?;
        let mut second = [0.0; 3];
        for (axis, result) in second.iter_mut().enumerate() {
            *result = sum(boundaries.iter().enumerate().map(|(i, b)| {
                let dx = (b.centroid().x - anchor.x) - cx;
                let dy = (b.centroid().y - anchor.y) - cy;
                let shift = [dx * dx, dx * dy, dy * dy][axis];
                signed(i) * (b.moments.second[axis] + b.area() * shift)
            }));
            finite(*result, "profile moment")?;
        }
        if second[0] <= 0.0 || second[2] <= 0.0 {
            return Err(Error::Degenerate("profile moments"));
        }
        Ok(Self {
            outer,
            holes,
            tolerance,
            moments: AreaMoments {
                area,
                centroid,
                second,
            },
            perimeter,
        })
    }
    pub fn outer(&self) -> &Boundary {
        &self.outer
    }
    pub fn holes(&self) -> &[Boundary] {
        &self.holes
    }
    pub fn tolerance(&self) -> Tolerance {
        self.tolerance
    }
    pub fn area(&self) -> f64 {
        self.moments.area
    }
    pub fn centroid(&self) -> Point2 {
        self.moments.centroid
    }
    pub fn perimeter(&self) -> f64 {
        self.perimeter
    }
    pub fn classify(&self, point: Point2) -> Result<Location> {
        point.checked(self.tolerance)?;
        let outer = self.outer.locate(point, self.tolerance);
        if outer != Location::Inside {
            return Ok(outer);
        }
        for hole in &self.holes {
            match hole.locate(point, self.tolerance) {
                Location::Inside => return Ok(Location::Outside),
                Location::Boundary => return Ok(Location::Boundary),
                Location::Outside => {}
            }
        }
        Ok(Location::Inside)
    }
    pub(crate) fn boundaries(&self) -> impl Iterator<Item = &Boundary> {
        std::iter::once(&self.outer).chain(&self.holes)
    }
}

fn edges(points: &[Point2]) -> impl Iterator<Item = (Point2, Point2)> + '_ {
    points
        .iter()
        .copied()
        .zip(points.iter().copied().cycle().skip(1))
        .take(points.len())
}
fn cross(a: Point2, b: Point2) -> f64 {
    a.x * b.y - a.y * b.x
}
fn orient(a: Point2, b: Point2, c: Point2) -> f64 {
    cross(
        Point2::new(b.x - a.x, b.y - a.y),
        Point2::new(c.x - a.x, c.y - a.y),
    )
}
pub(crate) fn segments_touch(a: Point2, b: Point2, c: Point2, d: Point2, tolerance: f64) -> bool {
    let (ab_c, ab_d, cd_a, cd_b) = (
        orient2d_finite(a, b, c),
        orient2d_finite(a, b, d),
        orient2d_finite(c, d, a),
        orient2d_finite(c, d, b),
    );
    let opposite = |left, right| {
        left != Orientation2::Collinear && right != Orientation2::Collinear && left != right
    };
    let crosses = opposite(ab_c, ab_d) && opposite(cd_a, cd_b);
    crosses
        || decide::segment_distance_le(a, c, d, &[tolerance])
        || decide::segment_distance_le(b, c, d, &[tolerance])
        || decide::segment_distance_le(c, a, b, &[tolerance])
        || decide::segment_distance_le(d, a, b, &[tolerance])
}
fn boundaries_touch(a: &Boundary, b: &Boundary, tolerance: f64) -> bool {
    // S8b: with a spline on either, every pair of segments by the screen.
    if let (Some(ca), Some(cb)) = (a.spline_chains(), b.spline_chains()) {
        if a.has_spline() || b.has_spline() {
            return ca
                .iter()
                .any(|x| cb.iter().any(|y| !decide::splines::apart(x, y, tolerance)));
        }
    }
    if matches!(a.kind, BoundaryKind::Path { .. }) || matches!(b.kind, BoundaryKind::Path { .. }) {
        let (pa, pb) = (a.pieces(), b.pieces());
        return pa.iter().any(|x| {
            pb.iter()
                .any(|y| decide::arcs::pieces_within(x, y, tolerance))
        });
    }
    match (&a.kind, &b.kind) {
        (BoundaryKind::Polygon(a), BoundaryKind::Polygon(b)) => {
            edges(a).any(|(a, b_)| edges(b).any(|(c, d)| segments_touch(a, b_, c, d, tolerance)))
        }
        (
            BoundaryKind::Circle {
                center: a,
                radius: ar,
            },
            BoundaryKind::Circle {
                center: b,
                radius: br,
            },
        ) => {
            // |ar - br| exactly: the larger radius minus the smaller.
            let (big, small) = (ar.max(*br), ar.min(*br));
            decide::distance_le(*a, *b, &[*ar, *br, tolerance])
                && decide::distance_ge(*a, *b, &[big, -small, -tolerance])
        }
        (BoundaryKind::Circle { center, radius }, BoundaryKind::Polygon(points))
        | (BoundaryKind::Polygon(points), BoundaryKind::Circle { center, radius }) => edges(points)
            .any(|(a, b)| {
                decide::segment_distance_le(*center, a, b, &[*radius, tolerance])
                    && (decide::distance_ge(*center, a, &[*radius, -tolerance])
                        || decide::distance_ge(*center, b, &[*radius, -tolerance]))
            }),
        _ => unreachable!("paths are handled above"),
    }
}

/// A path's segments as pieces.
fn path_pieces(points: &[Point2], segments: &[Segment]) -> Vec<decide::arcs::Piece> {
    let n = points.len();
    segments
        .iter()
        .enumerate()
        .map(|(i, segment)| {
            let (a, b) = (points[i], points[(i + 1) % n]);
            match segment {
                // A spline's chord (S8b: spline paths are screened by
                // their parts, never by these pieces).
                Segment::Line | Segment::Spline(_) => decide::arcs::Piece::Line(a, b),
                Segment::Arc {
                    center,
                    radius,
                    ccw,
                } => decide::arcs::Piece::Arc(decide::arcs::Arc2 {
                    center: *center,
                    radius: *radius,
                    start: a,
                    end: b,
                    ccw: *ccw,
                }),
            }
        })
        .collect()
}

/// A path's area moments (of its material region, whichever way it runs)
/// and its signed area. Integrals are taken about the first point.
fn path_moments(points: &[Point2], segments: &[Segment]) -> Result<(AreaMoments, f64)> {
    let anchor = points[0];
    let n = points.len();
    let local = |p: Point2| Point2::new(p.x - anchor.x, p.y - anchor.y);
    let mut total = [0.0; 6];
    for (i, segment) in segments.iter().enumerate() {
        let (a, b) = (local(points[i]), local(points[(i + 1) % n]));
        let part = match segment {
            // S8b: exact Bernstein integrals about the anchor, rounded.
            Segment::Spline(span) => {
                let parts = spline_parts(span, points[i], points[(i + 1) % n])?;
                let r = |x: f64| num_rational::BigRational::from_float(x).expect("finite");
                decide::splines::green(&parts, &[r(anchor.x), r(anchor.y)])
                    .map(|x| decide::splines::to_f64(&x))
            }
            Segment::Line => {
                // The arcs' Green forms along the line (the polygon's
                // symmetric per-edge forms agree only over a closed loop).
                let (dx, dy) = (b.x - a.x, b.y - a.y);
                let cube = |p: f64, q: f64| p * p * p + p * p * q + p * q * q + q * q * q;
                let xxy = a.x * a.x * a.y
                    + (a.x * a.x * dy + 2.0 * a.x * dx * a.y) / 2.0
                    + (2.0 * a.x * dx * dy + dx * dx * a.y) / 3.0
                    + dx * dx * dy / 4.0;
                [
                    0.5 * cross(a, b),
                    dy * (a.x * a.x + a.x * b.x + b.x * b.x) / 6.0,
                    -dx * (a.y * a.y + a.y * b.y + b.y * b.y) / 6.0,
                    dy * cube(a.x, b.x) / 12.0,
                    dy * xxy / 2.0,
                    -dx * cube(a.y, b.y) / 12.0,
                ]
            }
            Segment::Arc {
                center,
                radius,
                ccw,
            } => {
                let c = local(*center);
                let start = (a.y - c.y).atan2(a.x - c.x);
                let sweep = arc_sweep(*center, points[i], points[(i + 1) % n], *ccw);
                arc_moments(c.x, c.y, *radius, start, sweep)
            }
        };
        for (t, p) in total.iter_mut().zip(part) {
            *t += p;
        }
    }
    let signed = total[0];
    // A clockwise path integrates to the negated moments.
    let k = if signed < 0.0 { -1.0 } else { 1.0 };
    let [area, mx, my, xx, xy, yy] = total.map(|x| k * x);
    let (cx, cy) = (mx / area, my / area);
    let second = [
        xx - area * cx * cx,
        xy - area * cx * cy,
        yy - area * cy * cy,
    ];
    for value in [area, cx, cy, second[0], second[1], second[2]] {
        finite(value, "area moment")?;
    }
    Ok((
        AreaMoments {
            area,
            centroid: Point2::new(anchor.x + cx, anchor.y + cy),
            second,
        },
        signed,
    ))
}
