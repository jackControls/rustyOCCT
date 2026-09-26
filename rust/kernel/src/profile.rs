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
    Circle { center: Point2, radius: f64 },
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
        let count = self.polygon_vertices().map_or(1, <[Point2]>::len);
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
        };
        // Stored points are already counter-clockwise: keep the caller's
        // orientation record and the stored-order labels.
        if rebuilt.polygon_vertices().map(<[Point2]>::len)
            != self.polygon_vertices().map(<[Point2]>::len)
        {
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
    pub fn new(outer: Boundary, holes: Vec<Boundary>, tolerance: Tolerance) -> Result<Self> {
        if holes.len() > MAX_HOLES {
            return Err(Error::LimitExceeded("profile holes"));
        }
        // Enforce the aggregate limit before repeating quadratic validation.
        let total_edges: usize = std::iter::once(&outer)
            .chain(&holes)
            .map(|boundary| boundary.polygon_vertices().map_or(1, <[Point2]>::len))
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
fn segments_touch(a: Point2, b: Point2, c: Point2, d: Point2, tolerance: f64) -> bool {
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
    }
}
