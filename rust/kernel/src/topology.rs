//! Analytic boundary representation: a cellular partition of space
//! (`TOPOLOGY_MODEL.md`, decisions D1–D12).
//!
//! A body is a set of regions, shells, faces, loops, fins, edges and
//! vertices. Region 0 is the infinite void; bounded regions are solid or void.
//! A shell is one connected boundary component of one region and lists face
//! sides. A face has two sides and stores the shell on each. A loop is an
//! ordered cycle of fins, with a winding number per periodic direction, or a
//! single vertex. A fin is one loop's oriented use of an edge and carries its
//! own pcurve. Each edge stores its fins in radial order. There are no seams:
//! closed curves without vertices are ring edges, and loops on a cylinder
//! close modulo the period.
//!
//! Slots are dense body-local indices. Vertices, edges, faces and bounded
//! regions also have value ids ([`EntityId`]) derived from how they were made
//! (see `identity.rs`); shells, loops, fins and the infinite void are
//! structure.
use crate::attributes::{Attribute, AttributeMap};
use crate::history::{EntityInfo, EntitySet, Geometry};
use crate::identity::{
    Derivation, EntityId, EntityKind, InputLabel, OperationId, OperationKind, Parent,
    ProfileElement, Role,
};
use crate::profile::{BoundaryKind, Segment};
use crate::{
    BSplineCurve2, BSplineCurve3, BSplineSurface3, Error, Frame3, Point2, Point3, Profile, Result,
    Tolerance, Vec3,
};
use std::collections::BTreeMap;
use std::f64::consts::TAU;

pub(crate) mod validate;
pub(crate) use validate::{conic_point_fast, projection_range, section_rates};
pub use validate::{EdgeEnd, Entity, Issue, IssueKind};

macro_rules! index_type {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(pub(crate) usize);
        impl $name {
            pub fn new(index: usize) -> Self {
                Self(index)
            }
            pub fn index(self) -> usize {
                self.0
            }
        }
    };
}
index_type!(VertexId);
index_type!(EdgeId);
index_type!(FinId);
index_type!(LoopId);
index_type!(FaceId);
index_type!(ShellId);
index_type!(RegionId);

/// An entity with identity, by its body-local position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Slot {
    Vertex(VertexId),
    Edge(EdgeId),
    Face(FaceId),
    /// A bounded region; region 0, the infinite void, has no id.
    Region(RegionId),
}

/// Ids, retained derivations and the slot maps of one body.
#[derive(Debug, Clone, PartialEq)]
struct Identity {
    body: EntityId,
    body_derivation: Derivation,
    vertices: Vec<EntityId>,
    edges: Vec<EntityId>,
    faces: Vec<EntityId>,
    /// Region r at index r - 1.
    regions: Vec<EntityId>,
    derivations: BTreeMap<EntityId, Derivation>,
    slots: BTreeMap<EntityId, Slot>,
    /// Profile labels, for deriving builder provenance.
    labels: BTreeMap<InputLabel, (u32, ProfileElement)>,
}

impl Identity {
    /// Ids from per-slot derivations; a repeated id is an operation error.
    fn new(
        body: Derivation,
        derivations: Vec<(Slot, Derivation)>,
        labels: BTreeMap<InputLabel, (u32, ProfileElement)>,
    ) -> Result<Self> {
        let mut identity = Self {
            body: body.id(),
            body_derivation: body,
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            regions: Vec::new(),
            derivations: BTreeMap::new(),
            slots: BTreeMap::new(),
            labels,
        };
        let mut ordered = derivations;
        ordered.sort_by_key(|(slot, _)| *slot);
        for (slot, derivation) in ordered {
            let id = derivation.id();
            if id == identity.body || identity.slots.insert(id, slot).is_some() {
                return Err(Error::InvalidTopology("id collision"));
            }
            identity.derivations.insert(id, derivation);
            let list = match slot {
                Slot::Vertex(v) => (&mut identity.vertices, v.0),
                Slot::Edge(e) => (&mut identity.edges, e.0),
                Slot::Face(f) => (&mut identity.faces, f.0),
                Slot::Region(r) => (&mut identity.regions, r.0.wrapping_sub(1)),
            };
            if list.0.len() != list.1 {
                return Err(Error::InvalidTopology("slot without a derivation"));
            }
            list.0.push(id);
        }
        Ok(identity)
    }

    /// `External` derivations for caller-supplied parts: one per slot.
    fn external(vertices: usize, edges: usize, faces: usize, regions: usize) -> Result<Self> {
        let d = |entity, ordinal: usize| Derivation {
            operation: OperationId::UNSPECIFIED,
            kind: OperationKind::External,
            entity,
            role: Role::External,
            ordinal: ordinal as u32,
            parents: Vec::new(),
        };
        let slots = (0..vertices)
            .map(|i| (Slot::Vertex(VertexId(i)), d(EntityKind::Vertex, i)))
            .chain((0..edges).map(|i| (Slot::Edge(EdgeId(i)), d(EntityKind::Edge, i))))
            .chain((0..faces).map(|i| (Slot::Face(FaceId(i)), d(EntityKind::Face, i))))
            .chain((1..regions).map(|i| (Slot::Region(RegionId(i)), d(EntityKind::Region, i))))
            .collect();
        Self::new(d(EntityKind::Body, 0), slots, BTreeMap::new())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    Reversed,
}

impl Orientation {
    pub fn sign(self) -> f64 {
        if self == Self::Forward {
            1.0
        } else {
            -1.0
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Curve3 {
    LineSegment {
        start: Point3,
        end: Point3,
    },
    /// A full turn from the frame's x axis: the arc with start 0 and sweep TAU.
    Circle {
        frame: Frame3,
        radius: f64,
    },
    /// frame.origin + radius (cos a x + sin a y), a = start + sweep * fraction.
    CircularArc {
        frame: Frame3,
        radius: f64,
        start_angle: f64,
        sweep_angle: f64,
    },
    /// A rational B-spline over a range, the edge fraction mapped affinely
    /// onto it (S4 of REVIEW_NOTES.md). A ring edge's is a full period.
    BSpline(SplineSpan<BSplineCurve3>),
    /// frame.origin + major cos a x + minor sin a y, a = start + sweep *
    /// fraction (OCCT's `Geom_Ellipse`; S8a.2): a plane's section of a
    /// cylinder. A ring edge's sweep is a full turn.
    EllipseArc {
        frame: Frame3,
        major: f64,
        minor: f64,
        start_angle: f64,
        sweep_angle: f64,
    },
    /// frame.origin + major cosh t x + minor sinh t y, t = start + sweep *
    /// fraction (OCCT's `Geom_Hyperbola`; S8d.2): a branch of a plane's
    /// section of a cone.
    HyperbolaArc {
        frame: Frame3,
        major: f64,
        minor: f64,
        start: f64,
        sweep: f64,
    },
    /// frame.origin + t^2 / (4 focal) x + t y, t = start + sweep * fraction
    /// (OCCT's `Geom_Parabola`; S8d.2): a plane's section of a cone parallel
    /// to one of its rulings.
    ParabolaArc {
        frame: Frame3,
        focal: f64,
        start: f64,
        sweep: f64,
    },
    /// A plane's section of a torus (S8d.3), a graph over one of its angles.
    Section(Box<Spiric>),
}

/// A plane's section of a torus as a graph over one of its angles (S8d.3).
/// In `frame` (the torus's, `major` and `minor` its radii) the plane is
/// `a x + b y + c z + d = 0` (`plane`, a unit normal). Over `u` (`over_v`
/// false), `u = start + sweep f` and `v = psi + sign acos(-(major alpha +
/// d) / W)`, with `alpha = a cos u + b sin u` and `(W cos psi, W sin psi) =
/// (minor alpha, minor c)`; over `v`, `v = start + sweep f` and `u = atan2(b,
/// a) + sign acos(q / |(a, b)|)`, with `q = -(c minor sin v + d) / (major +
/// minor cos v)`. The point is `frame.point(((major + minor cos v) cos u,
/// (major + minor cos v) sin u), minor sin v)`. An edge's range keeps each
/// `acos` argument strictly inside `(-1, 1)` (no turning point), so it is
/// analytic.
#[derive(Debug, Clone, PartialEq)]
pub struct Spiric {
    pub frame: Frame3,
    pub major: f64,
    pub minor: f64,
    pub plane: [f64; 4],
    pub over_v: bool,
    pub sign: f64,
    pub start: f64,
    pub sweep: f64,
}

impl Spiric {
    /// The torus's angles `(u, v)` at a fraction (principal values).
    pub fn angles(&self, fraction: f64) -> (f64, f64) {
        let [a, b, c, d] = self.plane;
        let (big, small) = (self.major, self.minor);
        let t = self.start + self.sweep * fraction;
        if self.over_v {
            let q = -(c * small * t.sin() + d) / (big + small * t.cos());
            let u = b.atan2(a) + self.sign * (q / a.hypot(b)).clamp(-1.0, 1.0).acos();
            (u, t)
        } else {
            let alpha = a * t.cos() + b * t.sin();
            let (x, y) = (small * alpha, small * c);
            let ratio = -(big * alpha + d) / x.hypot(y);
            (t, y.atan2(x) + self.sign * ratio.clamp(-1.0, 1.0).acos())
        }
    }

    pub fn point(&self, fraction: f64) -> Point3 {
        let (u, v) = self.angles(fraction);
        let rho = self.major + self.minor * v.cos();
        self.frame.point(
            Point2::new(rho * u.cos(), rho * u.sin()),
            self.minor * v.sin(),
        )
    }
}

/// A pcurve defined as the exact inverse of its face's surface map applied
/// to its fin's edge (D13; S8d.2): at a fraction `f` the surface's
/// parameters of the edge's point at `f` (or `1 - f` for a reversed use),
/// continuous on the surface's cover, `u` (and a torus's `v`) lifted to the
/// representative nearest the interpolation of `lifts` (the parameters at
/// fractions `k / (lifts.len() - 1)`, recorded when the pcurve is made).
#[derive(Debug, Clone, PartialEq)]
pub struct Projection {
    pub curve: Curve3,
    pub surface: Surface,
    pub reversed: bool,
    pub lifts: Vec<Point2>,
}

impl Projection {
    /// The surface's principal parameters of a point, before lifting.
    pub(crate) fn inverse(surface: &Surface, p: Point3) -> Option<Point2> {
        let frame = match surface {
            Surface::Plane(f)
            | Surface::Cylinder { frame: f, .. }
            | Surface::Cone { frame: f, .. }
            | Surface::Sphere { frame: f, .. }
            | Surface::Torus { frame: f, .. } => f,
            Surface::BSpline(_) => return None,
        };
        let [x, y, z] = frame.coordinates(p);
        Some(match surface {
            Surface::Plane(_) => Point2::new(x, y),
            Surface::Cylinder { .. } => Point2::new(y.atan2(x), z),
            Surface::Cone { half_angle, .. } => Point2::new(y.atan2(x), z / half_angle.cos()),
            Surface::Sphere { .. } => Point2::new(y.atan2(x), z.atan2(x.hypot(y))),
            Surface::Torus { major, .. } => Point2::new(y.atan2(x), z.atan2(x.hypot(y) - major)),
            Surface::BSpline(_) => unreachable!("returned above"),
        })
    }

    /// The lift at a fraction: the linear interpolation of the recorded ones.
    pub(crate) fn lift(&self, fraction: f64) -> Point2 {
        let n = self.lifts.len() - 1;
        let x = (fraction.clamp(0.0, 1.0) * n as f64).min(n as f64);
        let k = (x.floor() as usize).min(n.saturating_sub(1));
        let w = x - k as f64;
        let (a, b) = (self.lifts[k], self.lifts[(k + 1).min(n)]);
        Point2::new(a.x + (b.x - a.x) * w, a.y + (b.y - a.y) * w)
    }

    /// Whether the surface is periodic in u (every surface of revolution)
    /// and in v (a torus).
    pub(crate) fn periodic(&self) -> [bool; 2] {
        match self.surface {
            Surface::Plane(_) | Surface::BSpline(_) => [false, false],
            Surface::Torus { .. } => [true, true],
            _ => [true, false],
        }
    }

    pub fn point(&self, fraction: f64) -> Point2 {
        let f = if self.reversed {
            1.0 - fraction
        } else {
            fraction
        };
        let p = self.curve.point(f);
        let Some(mut uv) = Self::inverse(&self.surface, p) else {
            return Point2::new(f64::NAN, f64::NAN);
        };
        let lift = self.lift(fraction);
        let [pu, pv] = self.periodic();
        let near = |x: f64, target: f64| x + TAU * ((target - x) / TAU).round();
        if pu {
            uv.x = near(uv.x, lift.x);
        }
        if pv {
            uv.y = near(uv.y, lift.y);
        }
        uv
    }

    /// A projection of `curve` onto `surface` with `anchors` recorded lifts,
    /// each the representative nearest the previous one (from `start` at
    /// fraction 0); `None` when consecutive anchors turn by more than a
    /// quarter period, where the interpolation would not pin the lift.
    pub(crate) fn new(
        curve: Curve3,
        surface: Surface,
        reversed: bool,
        start: Point2,
        anchors: usize,
    ) -> Option<Self> {
        let mut out = Self {
            curve,
            surface,
            reversed,
            lifts: vec![start],
        };
        let [pu, pv] = out.periodic();
        let steps = anchors.max(1) * 8;
        let mut last = start;
        let mut lifts = vec![start];
        for k in 1..=steps {
            let fraction = k as f64 / steps as f64;
            let f = if reversed { 1.0 - fraction } else { fraction };
            let mut uv = Self::inverse(&out.surface, out.curve.point(f))?;
            let near = |x: f64, target: f64| x + TAU * ((target - x) / TAU).round();
            if pu {
                uv.x = near(uv.x, last.x);
                if (uv.x - last.x).abs() > std::f64::consts::FRAC_PI_2 {
                    return None;
                }
            }
            if pv {
                uv.y = near(uv.y, last.y);
                if (uv.y - last.y).abs() > std::f64::consts::FRAC_PI_2 {
                    return None;
                }
            }
            last = uv;
            if k % 8 == 0 {
                lifts.push(uv);
            }
        }
        out.lifts = lifts;
        Some(out)
    }
}

/// A spline over a closed range of its parameter: its whole domain when the
/// kernel builds it, a sub-range (or, periodic, any range of at most one
/// period) when a file trims it. An edge or pcurve fraction maps affinely
/// onto the range (S4 of REVIEW_NOTES.md).
#[derive(Debug, Clone, PartialEq)]
pub struct SplineSpan<C> {
    curve: C,
    range: [f64; 2],
    /// Traversed from the range's end to its start (a pcurve reversed for a
    /// reversed use, without mirroring its knots, which would round).
    reversed: bool,
}

/// The parameter domain and periodicity of a spline curve.
pub trait SplineDomain {
    fn domain(&self) -> (f64, f64);
    fn is_periodic(&self) -> bool;
}
impl SplineDomain for BSplineCurve3 {
    fn domain(&self) -> (f64, f64) {
        BSplineCurve3::domain(self)
    }
    fn is_periodic(&self) -> bool {
        BSplineCurve3::is_periodic(self)
    }
}
impl SplineDomain for BSplineCurve2 {
    fn domain(&self) -> (f64, f64) {
        self.as_curve3().domain()
    }
    fn is_periodic(&self) -> bool {
        self.as_curve3().is_periodic()
    }
}

impl<C: SplineDomain> SplineSpan<C> {
    /// The whole domain.
    pub fn whole(curve: C) -> Self {
        let (a, b) = curve.domain();
        Self {
            curve,
            range: [a, b],
            reversed: false,
        }
    }
    /// `first < last`, finite; inside the domain, or of at most one period
    /// on a periodic curve.
    pub fn new(curve: C, first: f64, last: f64) -> Result<Self> {
        let (a, b) = curve.domain();
        let fits = if curve.is_periodic() {
            crate::spline::rational(last) - crate::spline::rational(first)
                <= crate::spline::rational(b) - crate::spline::rational(a)
        } else {
            a <= first && last <= b
        };
        if !(first.is_finite() && last.is_finite() && first < last && fits) {
            return Err(Error::OutOfDomain("spline range"));
        }
        Ok(Self {
            curve,
            range: [first, last],
            reversed: false,
        })
    }
    /// The same span traversed the other way.
    pub fn reversed(&self) -> Self
    where
        C: Clone,
    {
        Self {
            curve: self.curve.clone(),
            range: self.range,
            reversed: !self.reversed,
        }
    }
    pub fn is_reversed(&self) -> bool {
        self.reversed
    }
    /// The curve's parameter at a fraction of the span, its ends exact.
    pub fn parameter(&self, fraction: f64) -> f64 {
        let t = if !self.reversed {
            fraction
        } else if fraction == 0.0 {
            1.0
        } else if fraction == 1.0 {
            0.0
        } else {
            1.0 - fraction
        };
        spline_parameter(self.range.into(), t)
    }
    pub fn curve(&self) -> &C {
        &self.curve
    }
    pub fn range(&self) -> [f64; 2] {
        self.range
    }
    /// A full period of a periodic curve: closed, its ends joined inside.
    pub fn is_closed_period(&self) -> bool {
        let (a, b) = self.curve.domain();
        self.curve.is_periodic()
            && crate::spline::rational(self.range[1]) - crate::spline::rational(self.range[0])
                == crate::spline::rational(b) - crate::spline::rational(a)
    }
}

impl SplineSpan<BSplineCurve2> {
    /// The same traversal without the direction flag: a flagged span is its
    /// mirrored curve (knots `k -> a + b - k`, poles and weights reversed,
    /// the range mirrored) run forward. `None` when a mirrored value is not
    /// a binary64 (the writer then cannot express it).
    pub fn unflagged(&self) -> Option<Self> {
        if !self.reversed {
            return Some(self.clone());
        }
        use crate::spline::rational;
        let c = self.curve.as_curve3();
        let (a, b) = c.domain();
        let mirror = |x: f64| -> Option<f64> {
            let m = rational(a) + rational(b) - rational(x);
            let f = a + b - x;
            (rational(f) == m).then_some(f)
        };
        let knots: Vec<f64> = c
            .knots()
            .iter()
            .rev()
            .map(|k| mirror(*k))
            .collect::<Option<_>>()?;
        let mults: Vec<usize> = c.multiplicities().iter().rev().copied().collect();
        let poles: Vec<Point2> = self.curve.poles().into_iter().rev().collect();
        let weights: Vec<f64> = c.weights().iter().rev().copied().collect();
        let curve = if c.is_periodic() {
            BSplineCurve2::new_periodic(c.degree(), poles, Some(weights), knots, mults)
        } else {
            BSplineCurve2::new(c.degree(), poles, Some(weights), knots, mults)
        }
        .ok()?;
        let [first, last] = self.range;
        SplineSpan::new(curve, mirror(last)?, mirror(first)?).ok()
    }
}

/// The spline parameter of an edge or pcurve fraction: affine onto the
/// range, its ends exact.
pub(crate) fn spline_parameter((a, b): (f64, f64), fraction: f64) -> f64 {
    if fraction >= 1.0 {
        b
    } else if fraction <= 0.0 {
        a
    } else {
        (a + (b - a) * fraction).clamp(a, b)
    }
}

impl Curve3 {
    /// Evaluate a normalized edge parameter (0..=1).
    pub fn point(&self, fraction: f64) -> Point3 {
        match self {
            Self::LineSegment { start, end } => *start + (*end - *start) * fraction,
            Self::Circle { frame, radius } => {
                let (sine, cosine) = (fraction * TAU).sin_cos();
                frame.point(Point2::new(radius * cosine, radius * sine), 0.0)
            }
            Self::CircularArc {
                frame,
                radius,
                start_angle,
                sweep_angle,
            } => {
                let (sine, cosine) = (start_angle + sweep_angle * fraction).sin_cos();
                frame.point(Point2::new(radius * cosine, radius * sine), 0.0)
            }
            Self::BSpline(span) => span
                .curve()
                .point(span.parameter(fraction))
                .expect("a finite spline evaluates in its domain"),
            Self::EllipseArc {
                frame,
                major,
                minor,
                start_angle,
                sweep_angle,
            } => {
                let (sine, cosine) = (start_angle + sweep_angle * fraction).sin_cos();
                frame.point(Point2::new(major * cosine, minor * sine), 0.0)
            }
            Self::HyperbolaArc {
                frame,
                major,
                minor,
                start,
                sweep,
            } => {
                let t = start + sweep * fraction;
                frame.point(Point2::new(major * t.cosh(), minor * t.sinh()), 0.0)
            }
            Self::ParabolaArc {
                frame,
                focal,
                start,
                sweep,
            } => {
                let t = start + sweep * fraction;
                frame.point(Point2::new(t * t / (4.0 * focal), t), 0.0)
            }
            Self::Section(s) => s.point(fraction),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Curve2 {
    LineSegment {
        start: Point2,
        end: Point2,
    },
    CircularArc {
        center: Point2,
        radius: f64,
        start_angle: f64,
        sweep_angle: f64,
    },
    /// A planar rational B-spline in the face's (u, v) over a range, the
    /// fraction mapped affinely onto it (S4).
    BSpline(SplineSpan<BSplineCurve2>),
    /// center + (major cos a, minor sin a), a = start + sweep * fraction:
    /// an ellipse with its axes along the plane's u and v (S8a.2).
    EllipseArc {
        center: Point2,
        major: f64,
        minor: f64,
        start_angle: f64,
        sweep_angle: f64,
    },
    /// u = start + sweep * fraction, v = a0 + a1 cos u + a2 sin u: a plane's
    /// section on a cylinder, the graph of its height over the angle (S8a.2).
    Sinusoid {
        start: f64,
        sweep: f64,
        a: [f64; 3],
    },
    /// The exact projection of the fin's own edge onto its face's surface
    /// (D13; S8d.2).
    Projection(Box<Projection>),
}

impl Curve2 {
    pub fn point(&self, fraction: f64) -> Point2 {
        match self {
            Self::LineSegment { start, end } => Point2::new(
                start.x + (end.x - start.x) * fraction,
                start.y + (end.y - start.y) * fraction,
            ),
            Self::CircularArc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => {
                let (sine, cosine) = (start_angle + sweep_angle * fraction).sin_cos();
                Point2::new(center.x + radius * cosine, center.y + radius * sine)
            }
            Self::BSpline(span) => span
                .curve()
                .point(span.parameter(fraction))
                .expect("a finite spline evaluates in its domain"),
            Self::EllipseArc {
                center,
                major,
                minor,
                start_angle,
                sweep_angle,
            } => {
                let (sine, cosine) = (start_angle + sweep_angle * fraction).sin_cos();
                Point2::new(center.x + major * cosine, center.y + minor * sine)
            }
            Self::Sinusoid { start, sweep, a } => {
                let u = start + sweep * fraction;
                let (sine, cosine) = u.sin_cos();
                Point2::new(u, a[0] + a[1] * cosine + a[2] * sine)
            }
            Self::Projection(p) => p.point(fraction),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Surface {
    Plane(Frame3),
    /// u is an angle in radians, v is axial distance in millimetres.
    Cylinder {
        frame: Frame3,
        radius: f64,
    },
    /// `O + (radius + v sin a)(cos u x + sin u y) + v cos a n` with `a` the
    /// half angle, `0 < |a| < pi/2`: v runs along the generatrix, as in
    /// OCCT's `Geom_ConicalSurface`, and the apex is at `v = -radius / sin a`.
    Cone {
        frame: Frame3,
        radius: f64,
        half_angle: f64,
    },
    /// `O + radius (cos v (cos u x + sin u y) + sin v n)`: u the longitude,
    /// v the latitude in `[-pi/2, pi/2]`, as OCCT's `Geom_SphericalSurface`;
    /// the poles are at `v = +-pi/2`.
    Sphere {
        frame: Frame3,
        radius: f64,
    },
    /// `O + (major + minor cos v)(cos u x + sin u y) + minor sin v n`: a ring
    /// torus (`major > minor`), u the longitude, v the angle round the tube
    /// from its outer equator, as OCCT's `Geom_ToroidalSurface`. Periodic in
    /// both; loops may wind in u or in v.
    Torus {
        frame: Frame3,
        major: f64,
        minor: f64,
    },
    /// A rational B-spline surface in its own (u, v) (S4). Its loops do not
    /// wind until the rest of S4 gives windings the domain's period.
    BSpline(BSplineSurface3),
}

impl Surface {
    pub fn point(&self, uv: Point2) -> Point3 {
        match self {
            Self::Plane(frame) => frame.point(uv, 0.0),
            Self::Cylinder { frame, radius } => {
                frame.point(Point2::new(radius * uv.x.cos(), radius * uv.x.sin()), uv.y)
            }
            Self::Cone {
                frame,
                radius,
                half_angle,
            } => {
                let rho = radius + uv.y * half_angle.sin();
                frame.point(
                    Point2::new(rho * uv.x.cos(), rho * uv.x.sin()),
                    uv.y * half_angle.cos(),
                )
            }
            Self::Sphere { frame, radius } => {
                let rho = radius * uv.y.cos();
                frame.point(
                    Point2::new(rho * uv.x.cos(), rho * uv.x.sin()),
                    radius * uv.y.sin(),
                )
            }
            Self::Torus {
                frame,
                major,
                minor,
            } => {
                let rho = major + minor * uv.y.cos();
                frame.point(
                    Point2::new(rho * uv.x.cos(), rho * uv.x.sin()),
                    minor * uv.y.sin(),
                )
            }
            Self::BSpline(surface) => surface
                .point(uv.x, uv.y)
                .expect("a finite spline surface evaluates in its domain"),
        }
    }
    /// The unit normal of the parametrization, `S_u x S_v` normalized; on a
    /// cone, for the nappe where `radius + v sin a` is positive.
    pub fn normal(&self, uv: Point2) -> Vec3 {
        match self {
            Self::Plane(frame) => frame.normal(),
            Self::Cylinder { frame, .. } => frame.x() * uv.x.cos() + frame.y() * uv.x.sin(),
            Self::Cone {
                frame, half_angle, ..
            } => {
                let radial = frame.x() * uv.x.cos() + frame.y() * uv.x.sin();
                radial * half_angle.cos() - frame.normal() * half_angle.sin()
            }
            // Outward, away from the poles, or from the tube's core circle.
            Self::Sphere { frame, .. } | Self::Torus { frame, .. } => {
                let radial = frame.x() * uv.x.cos() + frame.y() * uv.x.sin();
                radial * uv.y.cos() + frame.normal() * uv.y.sin()
            }
            Self::BSpline(surface) => {
                use crate::curve::KnotSide;
                // The right-hand jet, the left-hand one at a closed domain's
                // upper end.
                let side = |x: f64, k: &crate::KnotVector| {
                    if !k.is_periodic() && x >= k.domain().1 {
                        KnotSide::Left
                    } else {
                        KnotSide::Right
                    }
                };
                let sides = [side(uv.x, surface.u_knots()), side(uv.y, surface.v_knots())];
                let e = surface
                    .evaluate(uv.x, uv.y, crate::curve::DerivativeOrder::First, sides)
                    .expect("a finite spline surface evaluates in its domain");
                let d = |u, v| {
                    let b = e.derivative_bounds(u, v).expect("requested");
                    Vec3::new(
                        b[0].representative(),
                        b[1].representative(),
                        b[2].representative(),
                    )
                };
                let n = d(1, 0).cross(d(0, 1));
                n * (1.0 / n.length())
            }
        }
    }
    /// Periodic in v too: a torus.
    pub fn is_periodic_v(&self) -> bool {
        matches!(self, Self::Torus { .. })
    }
    /// Periodic in u (an angle): cylinders, cones, spheres and tori.
    pub fn is_periodic(&self) -> bool {
        !matches!(self, Self::Plane(_) | Self::BSpline(_))
    }
}

/// Where an enclosure came from (Contract 5 of `IDENTITY_AND_HISTORY.md`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Provenance {
    /// A certified upper bound computed from the stored geometry, or carried
    /// from an operation's inputs.
    Computed,
    /// A tolerance another system stored (an OCCT `.brep` file); the checker
    /// still verifies it.
    Imported,
}

/// An upper bound on a representation gap, in the body's length unit: for a
/// vertex, its distance to the ends of its edges' curves (and to the surface
/// of a vertex loop); for a fin, the distance between the edge curve and the
/// pcurve's image; for a face, the gaps between consecutive fins in its
/// parameter space (angles scaled by the radius). It never exceeds the
/// body's resolution, and `Topology::check` verifies it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Enclosure {
    pub bound: f64,
    pub provenance: Provenance,
}

impl Enclosure {
    pub fn computed(bound: f64) -> Self {
        Self {
            bound,
            provenance: Provenance::Computed,
        }
    }
    pub fn imported(bound: f64) -> Self {
        Self {
            bound,
            provenance: Provenance::Imported,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Vertex {
    pub position: Point3,
    pub enclosure: Option<Enclosure>,
}

/// A curve bounded by vertices, or a ring edge (a closed curve with neither).
#[derive(Debug, Clone, PartialEq)]
pub struct Edge {
    pub start: Option<VertexId>,
    pub end: Option<VertexId>,
    pub curve: Curve3,
    /// Every fin of this edge, in radial order about the curve tangent.
    pub fins: Vec<FinId>,
}

impl Edge {
    pub fn is_ring(&self) -> bool {
        self.start.is_none() && self.end.is_none()
    }
}

/// One loop's oriented use of an edge.
#[derive(Debug, Clone, PartialEq)]
pub struct Fin {
    pub edge: EdgeId,
    /// Against the edge curve.
    pub sense: Orientation,
    /// This edge's curve in the owning face's parameter space, in traversal
    /// order. On a periodic surface it lives in the universal cover.
    pub pcurve: Curve2,
    pub enclosure: Option<Enclosure>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Loop {
    /// An ordered cycle of fins. It closes modulo the surface period with
    /// `winding[d]` turns in periodic direction d (u, v); both are zero on a
    /// plane, and a cylinder is periodic in u only.
    Edges { fins: Vec<FinId>, winding: [i32; 2] },
    /// A pole or an immersed vertex.
    Vertex(VertexId),
}

/// Which side of a face: the oriented normal points from the front side's
/// region into the back side's region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Side {
    Front,
    Back,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RegionKind {
    Solid,
    Void,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Shell {
    pub region: RegionId,
    pub sides: Vec<(FaceId, Side)>,
    /// Edges with no fins that belong to this shell.
    pub wire_edges: Vec<EdgeId>,
    pub acorn_vertices: Vec<VertexId>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    pub kind: RegionKind,
    /// The first shell of a bounded region is its outer boundary.
    pub shells: Vec<ShellId>,
}

/// Builder provenance: useful when translating an application's feature and
/// sketch-entity identities. Boundary 0 is the outer wire; holes start at 1.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceOrigin {
    StartCap,
    EndCap,
    Wall {
        boundary: usize,
        segment: usize,
    },
    /// Supplied through [`Topology::from_parts`] rather than a builder.
    External,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Face {
    pub surface: Surface,
    /// Against the surface normal.
    pub sense: Orientation,
    /// On a plane the outer loop comes first; traversal follows the oriented face.
    pub loops: Vec<LoopId>,
    pub front: ShellId,
    pub back: ShellId,
    pub enclosure: Option<Enclosure>,
}

impl Face {
    pub fn normal(&self, uv: Point2) -> Vec3 {
        self.surface.normal(uv) * self.sense.sign()
    }
}

/// Certified enclosures `[lo, hi]` of a sheet's area, a wire's length (zero
/// for an acorn) and its centre (S6).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeasureEnclosure {
    pub measure: [f64; 2],
    pub centre: [[f64; 2]; 3],
}

/// A body's class (D9): computed, never stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodyClass {
    /// Every face between a solid and a void region.
    Solid,
    /// Faces only, no solid region.
    Sheet,
    /// Wire edges only.
    Wire,
    /// One vertex.
    Acorn,
    /// Anything else.
    General,
}

impl BodyClass {
    pub fn name(self) -> &'static str {
        match self {
            BodyClass::Solid => "solid",
            BodyClass::Sheet => "sheet",
            BodyClass::Wire => "wire",
            BodyClass::Acorn => "acorn",
            BodyClass::General => "general",
        }
    }
}

/// Certified enclosures `[lo, hi]` of a body's mass properties at unit
/// density (REVIEW_NOTES.md U2): volume, surface area, centre of gravity and
/// the inertia tensor about it, in world axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassEnclosure {
    pub volume: [f64; 2],
    pub surface_area: [f64; 2],
    pub centroid: [[f64; 2]; 3],
    pub inertia: [[[f64; 2]; 3]; 3],
}

impl MassEnclosure {
    /// The midpoints of every enclosure.
    pub fn midpoints(&self) -> crate::MassProperties {
        let mid = |x: [f64; 2]| 0.5 * x[0] + 0.5 * x[1];
        crate::MassProperties {
            volume: mid(self.volume),
            surface_area: mid(self.surface_area),
            centroid: Point3::new(
                mid(self.centroid[0]),
                mid(self.centroid[1]),
                mid(self.centroid[2]),
            ),
            inertia: std::array::from_fn(|a| std::array::from_fn(|b| mid(self.inertia[a][b]))),
        }
    }
    /// The largest half width, relative to each value's magnitude (at least
    /// 1 for the centroid, in the body's length unit).
    pub fn relative_width(&self) -> f64 {
        let rel = |x: [f64; 2], scale: f64| 0.5 * (x[1] - x[0]) / scale.max(f64::MIN_POSITIVE);
        let mut worst = rel(self.volume, self.volume[1].abs())
            .max(rel(self.surface_area, self.surface_area[1].abs()));
        for c in &self.centroid {
            worst = worst.max(rel(*c, c[1].abs().max(1.0)));
        }
        worst
    }
}

/// The counts OCCT's `nbshapes` reports for the same body (a synthesis, see
/// [`Topology::occt_counts`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OcctCounts {
    pub vertices: usize,
    pub edges: usize,
    pub wires: usize,
    pub faces: usize,
    pub shells: usize,
    pub solids: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Topology {
    vertices: Vec<Vertex>,
    edges: Vec<Edge>,
    fins: Vec<Fin>,
    loops: Vec<Loop>,
    faces: Vec<Face>,
    shells: Vec<Shell>,
    /// regions[0] is the infinite void.
    regions: Vec<Region>,
    identity: Identity,
    /// Where each slot of a builder-made prism sits; empty for `from_parts`.
    layout: Vec<(Slot, Place)>,
    /// Opaque attributes by entity id (contract 4), each list sorted by key
    /// with one value per key.
    attributes: AttributeMap,
}

/// Where a slot of a prism sits along its axis (M3): on the lower or the
/// upper end, or spanning the prism.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Place {
    Low(End),
    High(End),
    /// A wall, a vertical edge or the solid region.
    Swept,
}

/// An entity on one end of a prism, with the swept entity it bounds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum End {
    Cap,
    /// A cap edge and the wall it bounds.
    Edge(FaceId),
    /// A cap vertex and the vertical edge through it.
    Vertex(EdgeId),
}

/// Unvalidated boundary data for [`Topology::from_parts`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TopologyParts {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub fins: Vec<Fin>,
    pub loops: Vec<Loop>,
    pub faces: Vec<Face>,
    pub shells: Vec<Shell>,
    pub regions: Vec<Region>,
}

impl TopologyParts {
    fn view(&self) -> validate::View<'_> {
        validate::View {
            vertices: &self.vertices,
            edges: &self.edges,
            fins: &self.fins,
            loops: &self.loops,
            faces: &self.faces,
            shells: &self.shells,
            regions: &self.regions,
        }
    }
    /// Every issue of the validation contract; empty means valid.
    pub fn check(&self, tolerance: Tolerance) -> Vec<Issue> {
        validate::check(&self.view(), tolerance)
    }
    /// The same parts with a measured, computed enclosure on every vertex,
    /// fin and face that has none: a certified upper bound of its gaps from
    /// the stored geometry. Entities whose gaps have no finite certified bound
    /// keep none, and [`TopologyParts::check`] reports them.
    pub fn with_measured_enclosures(mut self) -> Self {
        let m = validate::measure(&self.view());
        fill(&mut self.vertices, &m.vertices, |v| &mut v.enclosure);
        fill(&mut self.fins, &m.fins, |f| &mut f.enclosure);
        fill(&mut self.faces, &m.faces, |f| &mut f.enclosure);
        self
    }
}

fn fill<T>(
    items: &mut [T],
    bounds: &[Option<f64>],
    slot: impl Fn(&mut T) -> &mut Option<Enclosure>,
) {
    for (item, bound) in items.iter_mut().zip(bounds) {
        let enclosure = slot(item);
        if enclosure.is_none() {
            *enclosure = bound.map(Enclosure::computed);
        }
    }
}

impl Topology {
    /// The same topology with every spline edge traversed along its curve:
    /// an edge whose span is flagged reversed runs from its end to its start
    /// over the unflagged span, each of its fins used the other way (a fin's
    /// pcurve is in its use's direction, so it stays). For writers whose
    /// records follow the curve (S8b).
    pub(crate) fn spline_edges_along_curves(&self) -> Topology {
        let mut t = self.clone();
        for e in 0..t.edges.len() {
            let Curve3::BSpline(span) = &t.edges[e].curve else {
                continue;
            };
            if !span.is_reversed() {
                continue;
            }
            let span = span.reversed();
            let edge = &mut t.edges[e];
            edge.curve = Curve3::BSpline(span);
            std::mem::swap(&mut edge.start, &mut edge.end);
            for f in edge.fins.clone() {
                let fin = &mut t.fins[f.0];
                fin.sense = match fin.sense {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                };
            }
        }
        t
    }
}

impl Topology {
    fn view(&self) -> validate::View<'_> {
        validate::View {
            vertices: &self.vertices,
            edges: &self.edges,
            fins: &self.fins,
            loops: &self.loops,
            faces: &self.faces,
            shells: &self.shells,
            regions: &self.regions,
        }
    }
    /// Build a topology only if the complete validation contract holds;
    /// otherwise return every issue found.
    pub fn from_parts(
        parts: TopologyParts,
        tolerance: Tolerance,
    ) -> std::result::Result<Self, Vec<Issue>> {
        let issues = parts.check(tolerance);
        if !issues.is_empty() {
            return Err(issues);
        }
        let identity = Identity::external(
            parts.vertices.len(),
            parts.edges.len(),
            parts.faces.len(),
            parts.regions.len(),
        )
        .expect("distinct external ordinals give distinct ids");
        Ok(Self {
            vertices: parts.vertices,
            edges: parts.edges,
            fins: parts.fins,
            loops: parts.loops,
            faces: parts.faces,
            shells: parts.shells,
            regions: parts.regions,
            identity,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
        })
    }
    /// Validated parts under the given derivations (an operation's own
    /// body, S8a.2): every vertex, fin and face measured, then the whole
    /// contract checked; the issues otherwise.
    pub(crate) fn from_parts_named(
        parts: TopologyParts,
        tolerance: Tolerance,
        body: Derivation,
        derivations: Vec<(Slot, Derivation)>,
    ) -> std::result::Result<Self, Vec<Issue>> {
        let parts = parts.with_measured_enclosures();
        let issues = parts.check(tolerance);
        if !issues.is_empty() {
            return Err(issues);
        }
        let slots =
            parts.vertices.len() + parts.edges.len() + parts.faces.len() + parts.regions.len() - 1;
        let identity = match Identity::new(body, derivations, BTreeMap::new()) {
            Ok(identity) if identity.slots.len() == slots => identity,
            _ => return Err(Vec::new()),
        };
        Ok(Self {
            vertices: parts.vertices,
            edges: parts.edges,
            fins: parts.fins,
            loops: parts.loops,
            faces: parts.faces,
            shells: parts.shells,
            regions: parts.regions,
            identity,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
        })
    }
    /// The body's id; rigid transforms keep it.
    pub fn body_id(&self) -> EntityId {
        self.identity.body
    }
    pub fn body_derivation(&self) -> &Derivation {
        &self.identity.body_derivation
    }
    pub fn id_of(&self, slot: Slot) -> Option<EntityId> {
        match slot {
            Slot::Vertex(v) => self.identity.vertices.get(v.0),
            Slot::Edge(e) => self.identity.edges.get(e.0),
            Slot::Face(f) => self.identity.faces.get(f.0),
            Slot::Region(r) => self.identity.regions.get(r.0.wrapping_sub(1)),
        }
        .copied()
    }
    pub fn slot_of(&self, id: EntityId) -> Option<Slot> {
        self.identity.slots.get(&id).copied()
    }
    /// Why the entity has its id.
    pub fn derivation(&self, id: EntityId) -> Option<&Derivation> {
        self.identity.derivations.get(&id)
    }
    /// Every entity id with its slot, in id order.
    pub fn ids(&self) -> impl Iterator<Item = (EntityId, Slot)> + '_ {
        self.identity.slots.iter().map(|(id, slot)| (*id, *slot))
    }
    /// Builder provenance, derived from the face's derivation and the
    /// profile's labels.
    pub fn face_origin(&self, face: FaceId) -> Option<FaceOrigin> {
        let derivation = self.derivation(self.id_of(Slot::Face(face))?)?;
        Some(match derivation.role {
            Role::StartCap => FaceOrigin::StartCap,
            Role::EndCap => FaceOrigin::EndCap,
            Role::Wall => {
                let (boundary, element) = match derivation.parents.first()? {
                    Parent::Label(label) => *self.identity.labels.get(label)?,
                    Parent::Profile { boundary, element } => (*boundary, *element),
                    Parent::Entity(_) => return None,
                };
                let ProfileElement::Segment(segment) = element else {
                    return None;
                };
                FaceOrigin::Wall {
                    boundary: boundary as usize,
                    segment: segment as usize,
                }
            }
            Role::External => FaceOrigin::External,
            _ => return None,
        })
    }
    /// Every entity with identity as the history checker sees it.
    pub fn entity_set(&self, tolerance: Tolerance) -> EntitySet {
        let mut entities = BTreeMap::new();
        let vid = |v: Option<VertexId>| v.map(|v| self.identity.vertices[v.0]);
        let fid = |f: FaceId| self.identity.faces[f.0];
        for (id, slot) in self.ids() {
            let ordinal = self.identity.derivations[&id].ordinal;
            let (kind, geometry, structure) = match slot {
                Slot::Vertex(v) => (
                    EntityKind::Vertex,
                    Geometry::Point(self.vertices[v.0].position),
                    Vec::new(),
                ),
                Slot::Edge(e) => {
                    let edge = &self.edges[e.0];
                    let ends = [edge.start, edge.end]
                        .into_iter()
                        .filter_map(vid)
                        .map(|v| (v, Orientation::Forward))
                        .collect();
                    (
                        EntityKind::Edge,
                        Geometry::Curve(edge.curve.clone()),
                        vec![ends],
                    )
                }
                Slot::Face(f) => {
                    let face = &self.faces[f.0];
                    let loops = face
                        .loops
                        .iter()
                        .map(|l| match &self.loops[l.0] {
                            Loop::Edges { fins, .. } => fins
                                .iter()
                                .map(|k| {
                                    let fin = &self.fins[k.0];
                                    (self.identity.edges[fin.edge.0], fin.sense)
                                })
                                .collect(),
                            Loop::Vertex(v) => {
                                vec![(self.identity.vertices[v.0], Orientation::Forward)]
                            }
                        })
                        .collect();
                    (
                        EntityKind::Face,
                        Geometry::Surface {
                            surface: face.surface.clone(),
                            orientation: face.sense,
                        },
                        loops,
                    )
                }
                Slot::Region(r) => {
                    let region = &self.regions[r.0];
                    let shells = region
                        .shells
                        .iter()
                        .map(|s| {
                            self.shells[s.0]
                                .sides
                                .iter()
                                .map(|(f, side)| {
                                    let o = if *side == Side::Front {
                                        Orientation::Forward
                                    } else {
                                        Orientation::Reversed
                                    };
                                    (fid(*f), o)
                                })
                                .collect()
                        })
                        .collect();
                    (EntityKind::Region, Geometry::Region(region.kind), shells)
                }
            };
            entities.insert(
                id,
                EntityInfo {
                    kind,
                    ordinal,
                    geometry,
                    structure,
                },
            );
        }
        EntitySet {
            body: self.body_id(),
            tolerance,
            entities,
        }
    }
    /// Every issue of the validation contract; empty means valid.
    pub fn check(&self, tolerance: Tolerance) -> Vec<Issue> {
        validate::check(&self.view(), tolerance)
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }
    pub fn fins(&self) -> &[Fin] {
        &self.fins
    }
    pub fn loops(&self) -> &[Loop] {
        &self.loops
    }
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    pub fn shells(&self) -> &[Shell] {
        &self.shells
    }
    pub fn regions(&self) -> &[Region] {
        &self.regions
    }
    pub fn vertex(&self, id: VertexId) -> Option<&Vertex> {
        self.vertices.get(id.0)
    }
    pub fn edge(&self, id: EdgeId) -> Option<&Edge> {
        self.edges.get(id.0)
    }
    pub fn fin(&self, id: FinId) -> Option<&Fin> {
        self.fins.get(id.0)
    }
    pub fn face(&self, id: FaceId) -> Option<&Face> {
        self.faces.get(id.0)
    }
    pub fn face_ids(&self) -> impl Iterator<Item = FaceId> {
        (0..self.faces.len()).map(FaceId)
    }
    pub fn edge_ids(&self) -> impl Iterator<Item = EdgeId> {
        (0..self.edges.len()).map(EdgeId)
    }
    /// The fins of one face, loop by loop.
    pub fn face_fins(&self, face: FaceId) -> Vec<Vec<&Fin>> {
        self.faces[face.0]
            .loops
            .iter()
            .map(|l| match &self.loops[l.0] {
                Loop::Edges { fins, .. } => fins.iter().map(|k| &self.fins[k.0]).collect(),
                Loop::Vertex(_) => Vec::new(),
            })
            .collect()
    }
    /// The faces using an edge, one per fin, in the edge's radial order.
    pub fn incident_faces(&self, edge: EdgeId) -> Option<Vec<FaceId>> {
        let edge = self.edge(edge)?;
        let mut owner = BTreeMap::new();
        for (f, face) in self.faces.iter().enumerate() {
            for l in &face.loops {
                if let Loop::Edges { fins, .. } = &self.loops[l.0] {
                    for k in fins {
                        owner.insert(*k, FaceId(f));
                    }
                }
            }
        }
        Some(
            edge.fins
                .iter()
                .filter_map(|k| owner.get(k).copied())
                .collect(),
        )
    }

    /// V - E + sum over faces (2 - loops). A ring edge is a closed cell with
    /// no vertex, contributing nothing, and a vertex loop adds its vertex and
    /// a loop.
    pub fn euler_characteristic(&self) -> i64 {
        let edges = self.edges.iter().filter(|e| !e.is_ring()).count() as i64;
        self.vertices.len() as i64 - edges
            + self
                .faces
                .iter()
                .map(|face| 2 - face.loops.len() as i64)
                .sum::<i64>()
    }

    /// The counts OCCT's `nbshapes` would report for this body, synthesized
    /// by rule (TOPOLOGY_MODEL.md, T2): each face with loops winding a
    /// periodic direction gets one seam edge per such direction, and each
    /// ring edge those loops use gets one seam vertex; a cone face's pole
    /// gets one degenerated edge (its vertex is the pole's); a face's wound
    /// loops form one wire and every other edge loop its own wire; shells
    /// and solids are those of solid regions.
    /// The body's class (D9 of `TOPOLOGY_MODEL.md`, S6): computed from its
    /// regions and what its shells list, never stored.
    pub fn class(&self) -> BodyClass {
        let wires = self.shells.iter().any(|s| !s.wire_edges.is_empty());
        let acorns: usize = self.shells.iter().map(|s| s.acorn_vertices.len()).sum();
        let solid = self.regions.iter().any(|r| r.kind == RegionKind::Solid);
        let kind = |s: ShellId| self.regions[self.shells[s.0].region.0].kind;
        if !self.faces.is_empty() && !wires && acorns == 0 {
            if solid
                && self.faces.iter().all(|f| {
                    (kind(f.front) == RegionKind::Solid) != (kind(f.back) == RegionKind::Solid)
                })
            {
                return BodyClass::Solid;
            }
            if !solid {
                return BodyClass::Sheet;
            }
        }
        if self.faces.is_empty() && !solid {
            if wires && acorns == 0 {
                return BodyClass::Wire;
            }
            if !wires && acorns == 1 {
                return BodyClass::Acorn;
            }
        }
        BodyClass::General
    }

    pub fn occt_counts(&self) -> OcctCounts {
        let mut seams = 0;
        let mut degenerate = 0;
        let mut wires = 0;
        let mut closed_vertices = 0;
        for face in &self.faces {
            let mut wound = [false, false];
            for l in &face.loops {
                if let Loop::Edges { winding, .. } = &self.loops[l.0] {
                    if winding == &[0, 0] {
                        wires += 1;
                        continue;
                    }
                    for (d, w) in winding.iter().enumerate() {
                        wound[d] |= *w != 0;
                    }
                }
            }
            seams += wound.iter().filter(|w| **w).count();
            wires += usize::from(wound.iter().any(|w| *w));
            // OCCT closes a cone's or sphere's band at the pole with a
            // degenerated edge, and has one wherever a loop passes through
            // a pole.
            if validate::pole_position(face, &self.loops).is_some() {
                degenerate += 1;
            }
            degenerate += self.pole_passes(face);
            // A whole sphere (no edge loops) is OCCT's one wire: a seam from
            // pole to pole and a degenerated edge at each pole, their vertices.
            let whole = face
                .loops
                .iter()
                .all(|l| matches!(self.loops[l.0], Loop::Vertex(_)));
            if whole && matches!(face.surface, Surface::Sphere { .. }) {
                seams += 1;
                degenerate += 2;
                wires += 1;
                closed_vertices += 2;
            }
            // A whole torus is OCCT's one wire of its two seams (the meridian
            // and latitude circles through one vertex).
            if whole && matches!(face.surface, Surface::Torus { .. }) {
                seams += 2;
                wires += 1;
                closed_vertices += 1;
            }
        }
        let solid: Vec<&Region> = self
            .regions
            .iter()
            .filter(|r| r.kind == RegionKind::Solid)
            .collect();
        // S6: without a solid, a sheet of one face is a free face and of
        // several one shell per face-bearing shell (a closed shell's two
        // sides, twins, are one); a wire of one edge is a free edge and of
        // several one wire.
        let sheet_shells = if solid.is_empty() && self.faces.len() > 1 {
            let key = |s: &Shell, flip: bool| {
                let mut k: Vec<(usize, bool)> = s
                    .sides
                    .iter()
                    .map(|(f, side)| (f.0, (*side == Side::Front) != flip))
                    .collect();
                k.sort_unstable();
                k
            };
            let mut seen: Vec<Vec<(usize, bool)>> = Vec::new();
            let mut count = 0;
            for shell in self.shells.iter().filter(|s| !s.sides.is_empty()) {
                if !seen.contains(&key(shell, true)) {
                    count += 1;
                }
                seen.push(key(shell, false));
            }
            count
        } else {
            0
        };
        wires += self
            .shells
            .iter()
            .filter(|s| s.wire_edges.len() > 1)
            .count();
        OcctCounts {
            // OCCT closes every closed edge at a vertex (on a cylinder's or
            // cone's seam, on a disc, or alone in a wire).
            vertices: self.vertices.len()
                + self.edges.iter().filter(|e| e.is_ring()).count()
                + closed_vertices,
            edges: self.edges.len() + seams + degenerate,
            wires,
            faces: self.faces.len(),
            shells: solid.iter().map(|r| r.shells.len()).sum::<usize>() + sheet_shells,
            solids: solid.len(),
        }
    }

    /// How often a face's edge loops pass through a pole of its cone or
    /// sphere: consecutive fins meeting there with different `u` (OCCT
    /// joins them by a degenerated edge; the cell model's UV gap there has
    /// no length).
    pub fn pole_passes(&self, face: &Face) -> usize {
        let poles: Vec<f64> = match face.surface {
            Surface::Sphere { .. } => {
                vec![std::f64::consts::FRAC_PI_2, -std::f64::consts::FRAC_PI_2]
            }
            Surface::Cone {
                radius, half_angle, ..
            } => vec![-radius / half_angle.sin()],
            _ => return 0,
        };
        let mut passes = 0;
        for l in &face.loops {
            let Loop::Edges { fins, .. } = &self.loops[l.0] else {
                continue;
            };
            for (k, f) in fins.iter().enumerate() {
                let next = &fins[(k + 1) % fins.len()];
                let a = self.fins[f.0].pcurve.point(1.0);
                let b = self.fins[next.0].pcurve.point(0.0);
                let at_pole = poles
                    .iter()
                    .any(|p| (a.y - p).abs() <= 1e-12 * p.abs().max(1.0));
                let du = (a.x - b.x).rem_euclid(std::f64::consts::TAU);
                if at_pole && du.min(std::f64::consts::TAU - du) > 1e-9 {
                    passes += 1;
                }
            }
        }
        passes
    }

    /// Run the complete validation contract; the error names the first issue.
    /// See [`Topology::check`] for the complete typed report.
    pub fn validate(&self, tolerance: Tolerance) -> Result<()> {
        match self.check(tolerance).first() {
            None => Ok(()),
            Some(issue) => Err(Error::InvalidTopology(issue.kind.name())),
        }
    }

    /// A reference point for integration: the first vertex, or the first
    /// face's frame origin.
    fn reference_point(&self) -> [f64; 3] {
        if let Some(v) = self.vertices.first() {
            return v.position.to_array();
        }
        match self.faces.first().map(|f| &f.surface) {
            Some(
                Surface::Plane(f)
                | Surface::Cylinder { frame: f, .. }
                | Surface::Cone { frame: f, .. }
                | Surface::Sphere { frame: f, .. }
                | Surface::Torus { frame: f, .. },
            ) => f.origin().to_array(),
            Some(Surface::BSpline(s)) => s.poles()[0].to_array(),
            None => [0.0; 3],
        }
    }
    /// Certified mass properties of the body's solid regions (general, U2
    /// and S4d of REVIEW_NOTES.md). `None` when a face cannot be integrated
    /// (an arc pcurve on a curved surface, a periodic spline surface, a
    /// spline piece whose enclosure fails) or a certified division fails.
    pub fn mass_enclosure(&self) -> Option<MassEnclosure> {
        let e = validate::mass(&self.view(), self.reference_point())?;
        Some(MassEnclosure {
            volume: e.volume,
            surface_area: e.area,
            centroid: e.centroid,
            inertia: e.inertia,
        })
    }
    /// Certified measure of a sheet (its area), wire (its length) or acorn
    /// (zero) and its centre (S6); `None` for other classes, a wire with a
    /// spline edge, or a face that cannot be integrated.
    pub fn measure_enclosure(&self) -> Option<MeasureEnclosure> {
        if !matches!(
            self.class(),
            BodyClass::Sheet | BodyClass::Wire | BodyClass::Acorn
        ) {
            return None;
        }
        let e = validate::sheet_measure(&self.view(), self.reference_point())?;
        Some(MeasureEnclosure {
            measure: e.measure,
            centre: e.centre,
        })
    }
    /// A face's area and centre of gravity, from certified enclosures.
    pub fn face_area_and_centre(&self, face: FaceId) -> Option<(f64, Point3)> {
        let (area, centre) = validate::face_mass(&self.view(), face.0, self.reference_point())?;
        let mid = |x: [f64; 2]| 0.5 * x[0] + 0.5 * x[1];
        Some((
            mid(area),
            Point3::new(mid(centre[0]), mid(centre[1]), mid(centre[2])),
        ))
    }

    /// Every entity's attributes, by id; each list is sorted by key.
    pub fn attributes(&self) -> &AttributeMap {
        &self.attributes
    }
    /// The same topology with `attribute` on entity `id`, replacing any
    /// value under the same key.
    pub fn with_attribute(mut self, id: EntityId, attribute: Attribute) -> Result<Self> {
        if self.slot_of(id).is_none() {
            return Err(Error::OutOfDomain(
                "attribute on an id the body does not have",
            ));
        }
        let list = self.attributes.entry(id).or_default();
        list.retain(|a| a.key != attribute.key);
        list.push(attribute);
        list.sort();
        Ok(self)
    }
    /// An entity's enclosure bound: a vertex's or face's own, an edge's
    /// largest over its fins; `None` for regions or entities without one.
    pub(crate) fn enclosure_bound(&self, id: EntityId) -> Option<f64> {
        match self.slot_of(id)? {
            Slot::Vertex(v) => self.vertices[v.0].enclosure.map(|e| e.bound),
            Slot::Face(f) => self.faces[f.0].enclosure.map(|e| e.bound),
            Slot::Edge(e) => self.edges[e.0]
                .fins
                .iter()
                .filter_map(|k| self.fins[k.0].enclosure.map(|e| e.bound))
                .reduce(f64::max),
            Slot::Region(_) => None,
        }
    }
    /// Raise every vertex, face and fin enclosure to `parents(id)` of its
    /// entity (an edge's for a fin) where that is larger: an output's bound
    /// never falls below its inputs' (Contract 5, T5).
    pub(crate) fn raise_enclosures(&mut self, parents: impl Fn(EntityId) -> Option<f64>) {
        let raise = |e: &mut Option<Enclosure>, b: Option<f64>| {
            if let (Some(e), Some(b)) = (e.as_mut(), b) {
                e.bound = e.bound.max(b);
            }
        };
        let slots: Vec<(EntityId, Slot)> = self.ids().collect();
        for (id, slot) in slots {
            let b = parents(id);
            match slot {
                Slot::Vertex(v) => raise(&mut self.vertices[v.0].enclosure, b),
                Slot::Face(f) => raise(&mut self.faces[f.0].enclosure, b),
                Slot::Edge(e) => {
                    for k in self.edges[e.0].fins.clone() {
                        raise(&mut self.fins[k.0].enclosure, b);
                    }
                }
                Slot::Region(_) => {}
            }
        }
    }
    pub(crate) fn set_attributes(&mut self, attributes: AttributeMap) {
        self.attributes = attributes;
    }
    /// Where each slot of a builder-made prism sits, in slot order.
    pub(crate) fn layout(&self) -> &[(Slot, Place)] {
        &self.layout
    }
    /// The same structure under new ids: every slot gets the given
    /// derivation; a repeated id or an unnamed slot is an error.
    pub(crate) fn renamed(
        mut self,
        body: Derivation,
        derivations: Vec<(Slot, Derivation)>,
    ) -> Result<Self> {
        let named = derivations.len();
        let labels = std::mem::take(&mut self.identity.labels);
        self.identity = Identity::new(body, derivations, labels)?;
        self.attributes.clear();
        let slots =
            self.vertices.len() + self.edges.len() + self.faces.len() + self.regions.len() - 1;
        if named != slots || self.identity.slots.len() != slots {
            return Err(Error::InvalidTopology("slot without a derivation"));
        }
        Ok(self)
    }
    /// The same structure with the ids of `other`, a rigid copy of it.
    pub(crate) fn with_identity_of(mut self, other: &Topology) -> Self {
        debug_assert_eq!(self.layout, other.layout);
        self.identity = other.identity.clone();
        self.attributes.clear();
        self
    }

    pub(crate) fn prism(
        profile: &Profile,
        frame: Frame3,
        low: f64,
        high: f64,
        start_is_low: bool,
        operation: OperationId,
    ) -> Result<Self> {
        let tolerance = profile.tolerance();
        // Derivations per slot. "Bottom"/"top" roles name the start/end side.
        let derive = |entity, role, ordinal, parents| Derivation {
            operation,
            kind: OperationKind::Extrude,
            entity,
            role,
            ordinal,
            parents,
        };
        let (low_vertex, high_vertex, low_edge, high_edge, low_cap, high_cap) = if start_is_low {
            (
                Role::BottomVertex,
                Role::TopVertex,
                Role::BottomEdge,
                Role::TopEdge,
                Role::StartCap,
                Role::EndCap,
            )
        } else {
            (
                Role::TopVertex,
                Role::BottomVertex,
                Role::TopEdge,
                Role::BottomEdge,
                Role::EndCap,
                Role::StartCap,
            )
        };
        let mut derivations: Vec<(Slot, Derivation)> = Vec::new();
        let mut layout: Vec<(Slot, Place)> = vec![
            (Slot::Face(FaceId(0)), Place::Low(End::Cap)),
            (Slot::Face(FaceId(1)), Place::High(End::Cap)),
            (Slot::Region(RegionId(1)), Place::Swept),
        ];
        let mut labels = BTreeMap::new();
        let mut cap_parents = Vec::new();
        for (b, wire) in profile.boundaries().enumerate() {
            let b = b as u32;
            match wire.labels() {
                Some(l) => {
                    cap_parents.push(Parent::Label(l.boundary));
                    labels.insert(l.boundary, (b, ProfileElement::Boundary));
                    for (j, label) in l.segments.iter().enumerate() {
                        labels.insert(*label, (b, ProfileElement::Segment(j as u32)));
                    }
                    for (j, label) in l.vertices.iter().enumerate() {
                        labels.insert(*label, (b, ProfileElement::Vertex(j as u32)));
                    }
                }
                None => cap_parents.push(Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Boundary,
                }),
            }
        }
        derivations.push((
            Slot::Face(FaceId(0)),
            derive(EntityKind::Face, low_cap, 0, cap_parents.clone()),
        ));
        derivations.push((
            Slot::Face(FaceId(1)),
            derive(EntityKind::Face, high_cap, 0, cap_parents.clone()),
        ));
        derivations.push((
            Slot::Region(RegionId(1)),
            derive(EntityKind::Region, Role::Region, 0, cap_parents),
        ));
        let bottom_frame = Frame3::new(
            frame.point(Point2::default(), low),
            -frame.normal(),
            frame.x(),
            tolerance,
        )?;
        let top_frame = Frame3::new(
            frame.point(Point2::default(), high),
            frame.normal(),
            frame.x(),
            tolerance,
        )?;
        // One solid region bounded by shell 0 (every front side) inside the
        // infinite void bounded by shell 1 (every back side).
        let cap = |surface| Face {
            surface,
            sense: Orientation::Forward,
            loops: Vec::new(),
            front: ShellId(0),
            back: ShellId(1),
            enclosure: None,
        };
        let mut topology = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            fins: Vec::new(),
            loops: Vec::new(),
            faces: vec![
                cap(Surface::Plane(bottom_frame)),
                cap(Surface::Plane(top_frame)),
            ],
            shells: Vec::new(),
            regions: Vec::new(),
            identity: Identity::new(
                derive(EntityKind::Body, Role::Body, 0, Vec::new()),
                Vec::new(),
                BTreeMap::new(),
            )?,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
        };
        for (boundary, wire) in profile.boundaries().enumerate() {
            let inner = boundary > 0;
            let b = boundary as u32;
            let seg = |j: usize| match wire.labels() {
                Some(l) => Parent::Label(l.segments[j]),
                None => Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Segment(j as u32),
                },
            };
            let vert = |j: usize| match wire.labels() {
                Some(l) => Parent::Label(l.vertices[j]),
                None => Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Vertex(j as u32),
                },
            };
            let (forward, reversed) = if inner {
                (Orientation::Reversed, Orientation::Forward)
            } else {
                (Orientation::Forward, Orientation::Reversed)
            };
            match &wire.kind {
                BoundaryKind::Polygon(points) => {
                    let bottom = points
                        .iter()
                        .map(|p| topology.add_vertex(frame.point(*p, low)))
                        .collect::<Vec<_>>();
                    let top = points
                        .iter()
                        .map(|p| topology.add_vertex(frame.point(*p, high)))
                        .collect::<Vec<_>>();
                    let count = points.len();
                    let bottom_edges = (0..count)
                        .map(|i| topology.add_line(bottom[i], bottom[(i + 1) % count]))
                        .collect::<Vec<_>>();
                    let top_edges = (0..count)
                        .map(|i| topology.add_line(top[i], top[(i + 1) % count]))
                        .collect::<Vec<_>>();
                    let vertical = (0..count)
                        .map(|i| topology.add_line(bottom[i], top[i]))
                        .collect::<Vec<_>>();
                    // The walls of this boundary follow the faces made so far.
                    let wall = |j: usize| FaceId(topology.faces.len() + j);
                    for j in 0..count {
                        layout.extend([
                            (
                                Slot::Vertex(bottom[j]),
                                Place::Low(End::Vertex(vertical[j])),
                            ),
                            (Slot::Vertex(top[j]), Place::High(End::Vertex(vertical[j]))),
                            (Slot::Edge(bottom_edges[j]), Place::Low(End::Edge(wall(j)))),
                            (Slot::Edge(top_edges[j]), Place::High(End::Edge(wall(j)))),
                            (Slot::Edge(vertical[j]), Place::Swept),
                            (Slot::Face(wall(j)), Place::Swept),
                        ]);
                        let (v, e) = (EntityKind::Vertex, EntityKind::Edge);
                        derivations.push((
                            Slot::Vertex(bottom[j]),
                            derive(v, low_vertex, 0, vec![vert(j)]),
                        ));
                        derivations.push((
                            Slot::Vertex(top[j]),
                            derive(v, high_vertex, 0, vec![vert(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(bottom_edges[j]),
                            derive(e, low_edge, 0, vec![seg(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(top_edges[j]),
                            derive(e, high_edge, 0, vec![seg(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(vertical[j]),
                            derive(e, Role::Vertical, 0, vec![vert(j)]),
                        ));
                    }
                    topology.add_cap_loop(0, &bottom_edges, !inner, bottom_frame);
                    topology.add_cap_loop(1, &top_edges, inner, top_frame);
                    for i in 0..count {
                        let j = (i + 1) % count;
                        let (start, end) = if inner { (j, i) } else { (i, j) };
                        let origin = topology.vertices[bottom[start].0].position;
                        let tangent = topology.vertices[bottom[end].0].position - origin;
                        let side_frame =
                            Frame3::new(origin, tangent.cross(frame.normal()), tangent, tolerance)?;
                        let fins = [
                            (bottom_edges[i], forward),
                            (vertical[end], Orientation::Forward),
                            (top_edges[i], reversed),
                            (vertical[start], Orientation::Reversed),
                        ]
                        .iter()
                        .map(|(edge, sense)| topology.plane_fin(*edge, *sense, side_frame))
                        .collect();
                        derivations.push((
                            Slot::Face(FaceId(topology.faces.len())),
                            derive(EntityKind::Face, Role::Wall, 0, vec![seg(i)]),
                        ));
                        let l = topology.add_loop(fins, [0, 0]);
                        topology.faces.push(Face {
                            surface: Surface::Plane(side_frame),
                            sense: Orientation::Forward,
                            loops: vec![l],
                            front: ShellId(0),
                            back: ShellId(1),
                            enclosure: None,
                        });
                    }
                }
                BoundaryKind::Path { points, segments } => {
                    // S5: a polygon's structure, with circular-arc bottom and
                    // top edges and partial cylinder walls on arc segments.
                    let bottom = points
                        .iter()
                        .map(|p| topology.add_vertex(frame.point(*p, low)))
                        .collect::<Vec<_>>();
                    let top = points
                        .iter()
                        .map(|p| topology.add_vertex(frame.point(*p, high)))
                        .collect::<Vec<_>>();
                    let count = points.len();
                    let arc_of =
                        |i: usize, height: f64| -> Result<Option<(Frame3, f64, f64, f64)>> {
                            let Segment::Arc {
                                center,
                                radius,
                                ccw,
                            } = segments[i].clone()
                            else {
                                return Ok(None);
                            };
                            let arc_frame = Frame3::new(
                                frame.point(center, height),
                                frame.normal(),
                                frame.x(),
                                tolerance,
                            )?;
                            let a = points[i];
                            let start = (a.y - center.y).atan2(a.x - center.x);
                            let sweep =
                                crate::profile::arc_sweep(center, a, points[(i + 1) % count], ccw);
                            Ok(Some((arc_frame, radius, start, sweep)))
                        };
                    let mut bottom_edges = Vec::with_capacity(count);
                    let mut top_edges = Vec::with_capacity(count);
                    for i in 0..count {
                        let j = (i + 1) % count;
                        for (ends, height, edges) in [
                            (&bottom, low, &mut bottom_edges),
                            (&top, high, &mut top_edges),
                        ] {
                            let edge = match arc_of(i, height)? {
                                // S8b: a spline segment lifted to the cap.
                                None if matches!(segments[i], Segment::Spline(_)) => {
                                    let Segment::Spline(span) = &segments[i] else {
                                        unreachable!("matched above")
                                    };
                                    topology.add_edge(
                                        Some(ends[i]),
                                        Some(ends[j]),
                                        lifted_spline(span, frame, height)?,
                                    )
                                }
                                None => topology.add_line(ends[i], ends[j]),
                                Some((arc_frame, radius, start, sweep)) => topology.add_edge(
                                    Some(ends[i]),
                                    Some(ends[j]),
                                    Curve3::CircularArc {
                                        frame: arc_frame,
                                        radius,
                                        start_angle: start,
                                        sweep_angle: sweep,
                                    },
                                ),
                            };
                            edges.push(edge);
                        }
                    }
                    let vertical = (0..count)
                        .map(|i| topology.add_line(bottom[i], top[i]))
                        .collect::<Vec<_>>();
                    let wall = |j: usize| FaceId(topology.faces.len() + j);
                    for j in 0..count {
                        layout.extend([
                            (
                                Slot::Vertex(bottom[j]),
                                Place::Low(End::Vertex(vertical[j])),
                            ),
                            (Slot::Vertex(top[j]), Place::High(End::Vertex(vertical[j]))),
                            (Slot::Edge(bottom_edges[j]), Place::Low(End::Edge(wall(j)))),
                            (Slot::Edge(top_edges[j]), Place::High(End::Edge(wall(j)))),
                            (Slot::Edge(vertical[j]), Place::Swept),
                            (Slot::Face(wall(j)), Place::Swept),
                        ]);
                        let (v, e) = (EntityKind::Vertex, EntityKind::Edge);
                        derivations.push((
                            Slot::Vertex(bottom[j]),
                            derive(v, low_vertex, 0, vec![vert(j)]),
                        ));
                        derivations.push((
                            Slot::Vertex(top[j]),
                            derive(v, high_vertex, 0, vec![vert(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(bottom_edges[j]),
                            derive(e, low_edge, 0, vec![seg(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(top_edges[j]),
                            derive(e, high_edge, 0, vec![seg(j)]),
                        ));
                        derivations.push((
                            Slot::Edge(vertical[j]),
                            derive(e, Role::Vertical, 0, vec![vert(j)]),
                        ));
                    }
                    topology.add_cap_loop(0, &bottom_edges, !inner, bottom_frame);
                    topology.add_cap_loop(1, &top_edges, inner, top_frame);
                    let height = high - low;
                    for i in 0..count {
                        let j = (i + 1) % count;
                        let (start, end) = if inner { (j, i) } else { (i, j) };
                        derivations.push((
                            Slot::Face(FaceId(topology.faces.len())),
                            derive(EntityKind::Face, Role::Wall, 0, vec![seg(i)]),
                        ));
                        let (surface, sense, fins) = match arc_of(i, low)? {
                            // S8b: the spline's degree-(p, 1) wall, its
                            // pcurves lines in its (u, v) as an arc's are.
                            None if matches!(segments[i], Segment::Spline(_)) => {
                                let Segment::Spline(span) = &segments[i] else {
                                    unreachable!("matched above")
                                };
                                let (surface, a, b) = spline_wall(span, frame, low, high)?;
                                let fin = |edge, sense, from: (f64, f64), to: (f64, f64)| Fin {
                                    edge,
                                    sense,
                                    pcurve: Curve2::LineSegment {
                                        start: Point2::new(from.0, from.1),
                                        end: Point2::new(to.0, to.1),
                                    },
                                    enclosure: None,
                                };
                                let fins = if inner {
                                    vec![
                                        fin(
                                            bottom_edges[i],
                                            Orientation::Reversed,
                                            (b, 0.0),
                                            (a, 0.0),
                                        ),
                                        fin(
                                            vertical[i],
                                            Orientation::Forward,
                                            (a, 0.0),
                                            (a, height),
                                        ),
                                        fin(
                                            top_edges[i],
                                            Orientation::Forward,
                                            (a, height),
                                            (b, height),
                                        ),
                                        fin(
                                            vertical[j],
                                            Orientation::Reversed,
                                            (b, height),
                                            (b, 0.0),
                                        ),
                                    ]
                                } else {
                                    vec![
                                        fin(
                                            bottom_edges[i],
                                            Orientation::Forward,
                                            (a, 0.0),
                                            (b, 0.0),
                                        ),
                                        fin(
                                            vertical[j],
                                            Orientation::Forward,
                                            (b, 0.0),
                                            (b, height),
                                        ),
                                        fin(
                                            top_edges[i],
                                            Orientation::Reversed,
                                            (b, height),
                                            (a, height),
                                        ),
                                        fin(
                                            vertical[i],
                                            Orientation::Reversed,
                                            (a, height),
                                            (a, 0.0),
                                        ),
                                    ]
                                };
                                // Its normal is the right of increasing u: away
                                // from the material on an outer boundary run
                                // with u, on a hole run against it.
                                let sense = if span.is_reversed() == inner {
                                    Orientation::Forward
                                } else {
                                    Orientation::Reversed
                                };
                                (surface, sense, fins)
                            }
                            None => {
                                let origin = topology.vertices[bottom[start].0].position;
                                let tangent = topology.vertices[bottom[end].0].position - origin;
                                let side_frame = Frame3::new(
                                    origin,
                                    tangent.cross(frame.normal()),
                                    tangent,
                                    tolerance,
                                )?;
                                let fins = [
                                    (bottom_edges[i], forward),
                                    (vertical[end], Orientation::Forward),
                                    (top_edges[i], reversed),
                                    (vertical[start], Orientation::Reversed),
                                ]
                                .iter()
                                .map(|(edge, sense)| topology.plane_fin(*edge, *sense, side_frame))
                                .collect::<Vec<_>>();
                                (Surface::Plane(side_frame), Orientation::Forward, fins)
                            }
                            Some((arc_frame, radius, a, sweep)) => {
                                // The cover's rectangle [a, a + sweep] x
                                // [0, height], its normal leaving the
                                // material: outward on a counter-clockwise
                                // outer arc.
                                let b = a + sweep;
                                let ccw = sweep > 0.0;
                                let sense = if ccw != inner {
                                    Orientation::Forward
                                } else {
                                    Orientation::Reversed
                                };
                                let fin = |edge, sense, from: (f64, f64), to: (f64, f64)| Fin {
                                    edge,
                                    sense,
                                    pcurve: Curve2::LineSegment {
                                        start: Point2::new(from.0, from.1),
                                        end: Point2::new(to.0, to.1),
                                    },
                                    enclosure: None,
                                };
                                let fins = if inner {
                                    vec![
                                        fin(
                                            bottom_edges[i],
                                            Orientation::Reversed,
                                            (b, 0.0),
                                            (a, 0.0),
                                        ),
                                        fin(
                                            vertical[i],
                                            Orientation::Forward,
                                            (a, 0.0),
                                            (a, height),
                                        ),
                                        fin(
                                            top_edges[i],
                                            Orientation::Forward,
                                            (a, height),
                                            (b, height),
                                        ),
                                        fin(
                                            vertical[j],
                                            Orientation::Reversed,
                                            (b, height),
                                            (b, 0.0),
                                        ),
                                    ]
                                } else {
                                    vec![
                                        fin(
                                            bottom_edges[i],
                                            Orientation::Forward,
                                            (a, 0.0),
                                            (b, 0.0),
                                        ),
                                        fin(
                                            vertical[j],
                                            Orientation::Forward,
                                            (b, 0.0),
                                            (b, height),
                                        ),
                                        fin(
                                            top_edges[i],
                                            Orientation::Reversed,
                                            (b, height),
                                            (a, height),
                                        ),
                                        fin(
                                            vertical[i],
                                            Orientation::Reversed,
                                            (a, height),
                                            (a, 0.0),
                                        ),
                                    ]
                                };
                                (
                                    Surface::Cylinder {
                                        frame: arc_frame,
                                        radius,
                                    },
                                    sense,
                                    fins,
                                )
                            }
                        };
                        let l = topology.add_loop(fins, [0, 0]);
                        topology.faces.push(Face {
                            surface,
                            sense,
                            loops: vec![l],
                            front: ShellId(0),
                            back: ShellId(1),
                            enclosure: None,
                        });
                    }
                }
                BoundaryKind::Circle { center, radius } => {
                    let cylinder_frame = Frame3::new(
                        frame.point(*center, low),
                        frame.normal(),
                        frame.x(),
                        tolerance,
                    )?;
                    let end_frame = Frame3::new(
                        frame.point(*center, high),
                        frame.normal(),
                        frame.x(),
                        tolerance,
                    )?;
                    // Two ring edges and the wall; no seam, no vertex.
                    let bottom = topology.add_ring(Curve3::Circle {
                        frame: cylinder_frame,
                        radius: *radius,
                    });
                    let top = topology.add_ring(Curve3::Circle {
                        frame: end_frame,
                        radius: *radius,
                    });
                    let e = EntityKind::Edge;
                    let wall = FaceId(topology.faces.len());
                    layout.extend([
                        (Slot::Edge(bottom), Place::Low(End::Edge(wall))),
                        (Slot::Edge(top), Place::High(End::Edge(wall))),
                        (Slot::Face(wall), Place::Swept),
                    ]);
                    derivations.push((Slot::Edge(bottom), derive(e, low_edge, 0, vec![seg(0)])));
                    derivations.push((Slot::Edge(top), derive(e, high_edge, 0, vec![seg(0)])));
                    derivations.push((
                        Slot::Face(FaceId(topology.faces.len())),
                        derive(EntityKind::Face, Role::Wall, 0, vec![seg(0)]),
                    ));
                    topology.add_cap_loop(0, &[bottom], !inner, bottom_frame);
                    topology.add_cap_loop(1, &[top], inner, top_frame);
                    // The wall's loops wind once around the axis, in opposite
                    // directions, on the universal cover of the cylinder.
                    let (u0, u1, turns) = if inner { (TAU, 0.0, -1) } else { (0.0, TAU, 1) };
                    let height = high - low;
                    let fin = |edge, sense, start, end| Fin {
                        edge,
                        sense,
                        pcurve: Curve2::LineSegment { start, end },
                        enclosure: None,
                    };
                    let lower = topology.add_loop(
                        vec![fin(
                            bottom,
                            forward,
                            Point2::new(u0, 0.0),
                            Point2::new(u1, 0.0),
                        )],
                        [turns, 0],
                    );
                    let upper = topology.add_loop(
                        vec![fin(
                            top,
                            reversed,
                            Point2::new(u1, height),
                            Point2::new(u0, height),
                        )],
                        [-turns, 0],
                    );
                    topology.faces.push(Face {
                        surface: Surface::Cylinder {
                            frame: cylinder_frame,
                            radius: *radius,
                        },
                        sense: forward,
                        loops: vec![lower, upper],
                        front: ShellId(0),
                        back: ShellId(1),
                        enclosure: None,
                    });
                }
            }
        }
        for (k, fin) in topology.fins.iter().enumerate() {
            topology.edges[fin.edge.0].fins.push(FinId(k));
        }
        let fronts = topology.face_ids().map(|f| (f, Side::Front)).collect();
        let backs = topology.face_ids().map(|f| (f, Side::Back)).collect();
        topology.shells = vec![
            Shell {
                region: RegionId(1),
                sides: fronts,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
            Shell {
                region: RegionId(0),
                sides: backs,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
        ];
        topology.regions = vec![
            Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(1)],
            },
            Region {
                kind: RegionKind::Solid,
                shells: vec![ShellId(0)],
            },
        ];
        topology.identity = Identity::new(
            derive(EntityKind::Body, Role::Body, 0, Vec::new()),
            derivations,
            labels,
        )?;
        layout.sort_by_key(|(slot, _)| *slot);
        topology.layout = layout;
        let identity = &topology.identity;
        if (
            identity.vertices.len(),
            identity.edges.len(),
            identity.faces.len(),
            identity.regions.len(),
        ) != (
            topology.vertices.len(),
            topology.edges.len(),
            topology.faces.len(),
            topology.regions.len() - 1,
        ) {
            return Err(Error::InvalidTopology("slot without a derivation"));
        }
        topology.measure_enclosures(tolerance)?;
        topology.validate(tolerance)?;
        if topology.euler_characteristic() != 2 - 2 * profile.holes().len() as i64 {
            return Err(Error::InvalidTopology("unexpected shell genus"));
        }
        Ok(topology)
    }

    /// A face body (`face`: one planar face on `frame` bounded by every
    /// boundary, both sides in the infinite void) or a wire body (the first
    /// boundary's edges as wire edges of the void), S6 of REVIEW_NOTES.md.
    /// Vertices sit at the boundaries' stored points in the frame's plane,
    /// edges are their segments as stored (a circle one ring edge); the face
    /// is the frame's plane, its normal the frame's, its outer loop counter-
    /// clockwise. Every entity is generated from its profile element: the
    /// face from every boundary, an edge from its segment, a vertex from its
    /// point.
    pub(crate) fn planar_sheet(
        boundaries: &[&crate::Boundary],
        frame: Frame3,
        tolerance: Tolerance,
        face: bool,
        operation: OperationId,
    ) -> Result<Self> {
        let kind = if face {
            OperationKind::MakeFace
        } else {
            OperationKind::MakeWire
        };
        let derive = |entity, role, parents| Derivation {
            operation,
            kind,
            entity,
            role,
            ordinal: 0,
            parents,
        };
        let boundaries = if face { boundaries } else { &boundaries[..1] };
        let mut topology = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            fins: Vec::new(),
            loops: Vec::new(),
            faces: Vec::new(),
            shells: Vec::new(),
            regions: Vec::new(),
            identity: Identity::new(
                derive(EntityKind::Body, Role::Body, Vec::new()),
                Vec::new(),
                BTreeMap::new(),
            )?,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
        };
        if face {
            topology.faces.push(Face {
                surface: Surface::Plane(frame),
                sense: Orientation::Forward,
                loops: Vec::new(),
                front: ShellId(0),
                back: ShellId(0),
                enclosure: None,
            });
        }
        let mut derivations: Vec<(Slot, Derivation)> = Vec::new();
        let mut labels = BTreeMap::new();
        let mut face_parents = Vec::new();
        let mut wire_edges = Vec::new();
        for (boundary, wire) in boundaries.iter().enumerate() {
            let b = boundary as u32;
            match wire.labels() {
                Some(l) => {
                    face_parents.push(Parent::Label(l.boundary));
                    labels.insert(l.boundary, (b, ProfileElement::Boundary));
                    for (j, label) in l.segments.iter().enumerate() {
                        labels.insert(*label, (b, ProfileElement::Segment(j as u32)));
                    }
                    for (j, label) in l.vertices.iter().enumerate() {
                        labels.insert(*label, (b, ProfileElement::Vertex(j as u32)));
                    }
                }
                None => face_parents.push(Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Boundary,
                }),
            }
            let seg = |j: usize| match wire.labels() {
                Some(l) => Parent::Label(l.segments[j]),
                None => Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Segment(j as u32),
                },
            };
            let vert = |j: usize| match wire.labels() {
                Some(l) => Parent::Label(l.vertices[j]),
                None => Parent::Profile {
                    boundary: b,
                    element: ProfileElement::Vertex(j as u32),
                },
            };
            let (points, segments): (&[Point2], Option<&[Segment]>) = match &wire.kind {
                BoundaryKind::Polygon(points) => (points.as_slice(), None),
                BoundaryKind::Path { points, segments } => {
                    (points.as_slice(), Some(segments.as_slice()))
                }
                BoundaryKind::Circle { center, radius } => {
                    let circle = Frame3::new(
                        frame.point(*center, 0.0),
                        frame.normal(),
                        frame.x(),
                        tolerance,
                    )?;
                    let edge = topology.add_ring(Curve3::Circle {
                        frame: circle,
                        radius: *radius,
                    });
                    derivations.push((
                        Slot::Edge(edge),
                        derive(EntityKind::Edge, Role::Edge, vec![seg(0)]),
                    ));
                    wire_edges.push(edge);
                    if face {
                        topology.add_cap_loop(0, &[edge], boundary > 0, frame);
                    }
                    continue;
                }
            };
            let count = points.len();
            let vertices = points
                .iter()
                .map(|p| topology.add_vertex(frame.point(*p, 0.0)))
                .collect::<Vec<_>>();
            let mut edges = Vec::with_capacity(count);
            for i in 0..count {
                let j = (i + 1) % count;
                let edge = match segments.map(|s| s[i].clone()) {
                    None | Some(Segment::Line) => topology.add_line(vertices[i], vertices[j]),
                    Some(Segment::Spline(span)) => topology.add_edge(
                        Some(vertices[i]),
                        Some(vertices[j]),
                        lifted_spline(&span, frame, 0.0)?,
                    ),
                    Some(Segment::Arc {
                        center,
                        radius,
                        ccw,
                    }) => {
                        let arc_frame = Frame3::new(
                            frame.point(center, 0.0),
                            frame.normal(),
                            frame.x(),
                            tolerance,
                        )?;
                        let a = points[i];
                        topology.add_edge(
                            Some(vertices[i]),
                            Some(vertices[j]),
                            Curve3::CircularArc {
                                frame: arc_frame,
                                radius,
                                start_angle: (a.y - center.y).atan2(a.x - center.x),
                                sweep_angle: crate::profile::arc_sweep(center, a, points[j], ccw),
                            },
                        )
                    }
                };
                edges.push(edge);
                derivations.push((
                    Slot::Vertex(vertices[i]),
                    derive(EntityKind::Vertex, Role::Vertex, vec![vert(i)]),
                ));
                derivations.push((
                    Slot::Edge(edge),
                    derive(EntityKind::Edge, Role::Edge, vec![seg(i)]),
                ));
            }
            if face {
                topology.add_cap_loop(0, &edges, boundary > 0, frame);
            }
            wire_edges.extend(edges);
        }
        if face {
            derivations.push((
                Slot::Face(FaceId(0)),
                derive(EntityKind::Face, Role::Face, face_parents),
            ));
        }
        for (k, fin) in topology.fins.iter().enumerate() {
            topology.edges[fin.edge.0].fins.push(FinId(k));
        }
        topology.shells = vec![Shell {
            region: RegionId(0),
            sides: if face {
                vec![(FaceId(0), Side::Front), (FaceId(0), Side::Back)]
            } else {
                Vec::new()
            },
            wire_edges: if face { Vec::new() } else { wire_edges },
            acorn_vertices: Vec::new(),
        }];
        topology.regions = vec![Region {
            kind: RegionKind::Void,
            shells: vec![ShellId(0)],
        }];
        topology.identity = Identity::new(
            derive(EntityKind::Body, Role::Body, Vec::new()),
            derivations,
            labels,
        )?;
        topology.measure_enclosures(tolerance)?;
        topology.validate(tolerance)?;
        Ok(topology)
    }

    /// A right circular cone or frustum (S3 of REVIEW_NOTES.md), as
    /// `BRepPrimAPI_MakeCone(gp_Ax2, bottom, top, height)` makes it: the
    /// bottom at the frame origin, the top at `height` along its normal, a
    /// zero radius an apex. The lateral face is a `Cone` surface with v along
    /// the generatrix from the bottom; a nonzero end is a ring edge bounding
    /// a disc, a zero end a pole (a vertex loop at the apex). Entities derive
    /// from the meridian profile `(0, 0), (bottom, 0), (top, height), (0,
    /// height)` as a revolution: segment 0 (bottom radius) the start cap,
    /// segment 1 (generatrix) the wall, segment 2 (top radius) the end cap,
    /// vertex 1 the bottom ring or apex, vertex 2 the top ring or apex.
    pub(crate) fn cone(
        frame: Frame3,
        bottom: f64,
        top: f64,
        height: f64,
        tolerance: Tolerance,
        operation: OperationId,
    ) -> Result<Self> {
        let tol = tolerance.linear();
        for (value, what) in [
            (bottom, "cone radius"),
            (top, "cone radius"),
            (height, "cone height"),
        ] {
            crate::math::finite(value, what)?;
        }
        if bottom < 0.0 || top < 0.0 || (bottom > 0.0 && bottom <= tol) || (top > 0.0 && top <= tol)
        {
            return Err(Error::Degenerate("cone radius"));
        }
        if bottom == top {
            return Err(Error::OutOfDomain("equal cone radii make a cylinder"));
        }
        if height <= tol {
            return Err(Error::Degenerate("cone height"));
        }
        let derive = |entity, role, parents| Derivation {
            operation,
            kind: OperationKind::Revolve,
            entity,
            role,
            ordinal: 0,
            parents,
        };
        let meridian = |element| Parent::Profile {
            boundary: 0,
            element,
        };
        let half_angle = (top - bottom).atan2(height);
        let slant = height.hypot(top - bottom);
        let mut topology = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            fins: Vec::new(),
            loops: Vec::new(),
            faces: Vec::new(),
            shells: Vec::new(),
            regions: Vec::new(),
            identity: Identity::new(
                derive(EntityKind::Body, Role::Body, Vec::new()),
                Vec::new(),
                BTreeMap::new(),
            )?,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
        };
        let mut derivations: Vec<(Slot, Derivation)> = vec![(
            Slot::Region(RegionId(1)),
            derive(
                EntityKind::Region,
                Role::Region,
                vec![meridian(ProfileElement::Boundary)],
            ),
        )];
        let mut wall_loops = Vec::new();
        // (radius, height, cap role, cap segment, edge role, rim vertex)
        let ends = [
            (bottom, 0.0, Role::StartCap, 0, Role::BottomEdge, 1),
            (top, height, Role::EndCap, 2, Role::TopEdge, 2),
        ];
        for (radius, z, cap_role, segment, edge_role, rim) in ends {
            let upper = z > 0.0;
            if radius == 0.0 {
                let apex = topology.add_vertex(frame.point(Point2::default(), z));
                derivations.push((
                    Slot::Vertex(apex),
                    derive(
                        EntityKind::Vertex,
                        Role::Apex,
                        vec![meridian(ProfileElement::Vertex(rim))],
                    ),
                ));
                topology.loops.push(Loop::Vertex(apex));
                wall_loops.push(LoopId(topology.loops.len() - 1));
                continue;
            }
            let centre = frame.point(Point2::default(), z);
            let normal = if upper {
                frame.normal()
            } else {
                -frame.normal()
            };
            let disc_frame = Frame3::new(centre, normal, frame.x(), tolerance)?;
            let ring_frame = Frame3::new(centre, frame.normal(), frame.x(), tolerance)?;
            let ring = topology.add_ring(Curve3::Circle {
                frame: ring_frame,
                radius,
            });
            derivations.push((
                Slot::Edge(ring),
                derive(
                    EntityKind::Edge,
                    edge_role,
                    vec![meridian(ProfileElement::Vertex(rim))],
                ),
            ));
            let disc = FaceId(topology.faces.len());
            derivations.push((
                Slot::Face(disc),
                derive(
                    EntityKind::Face,
                    cap_role,
                    vec![meridian(ProfileElement::Segment(segment))],
                ),
            ));
            topology.faces.push(Face {
                surface: Surface::Plane(disc_frame),
                sense: Orientation::Forward,
                loops: Vec::new(),
                front: ShellId(0),
                back: ShellId(1),
                enclosure: None,
            });
            // The bottom disc runs against the ring, the top disc with it.
            topology.add_cap_loop(disc.0, &[ring], !upper, disc_frame);
            // The wall: the bottom ring once in +u at v = 0, the top ring once
            // in -u at v = slant, on the universal cover.
            let v = if upper { slant } else { 0.0 };
            let (sense, u0, u1, turns) = if upper {
                (Orientation::Reversed, TAU, 0.0, -1)
            } else {
                (Orientation::Forward, 0.0, TAU, 1)
            };
            let fin = Fin {
                edge: ring,
                sense,
                pcurve: Curve2::LineSegment {
                    start: Point2::new(u0, v),
                    end: Point2::new(u1, v),
                },
                enclosure: None,
            };
            wall_loops.push(topology.add_loop(vec![fin], [turns, 0]));
        }
        let wall = FaceId(topology.faces.len());
        derivations.push((
            Slot::Face(wall),
            derive(
                EntityKind::Face,
                Role::Wall,
                vec![meridian(ProfileElement::Segment(1))],
            ),
        ));
        topology.faces.push(Face {
            surface: Surface::Cone {
                frame,
                radius: bottom,
                half_angle,
            },
            sense: Orientation::Forward,
            loops: wall_loops,
            front: ShellId(0),
            back: ShellId(1),
            enclosure: None,
        });
        for (k, fin) in topology.fins.iter().enumerate() {
            topology.edges[fin.edge.0].fins.push(FinId(k));
        }
        let fronts = topology.face_ids().map(|f| (f, Side::Front)).collect();
        let backs = topology.face_ids().map(|f| (f, Side::Back)).collect();
        topology.shells = vec![
            Shell {
                region: RegionId(1),
                sides: fronts,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
            Shell {
                region: RegionId(0),
                sides: backs,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
        ];
        topology.regions = vec![
            Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(1)],
            },
            Region {
                kind: RegionKind::Solid,
                shells: vec![ShellId(0)],
            },
        ];
        topology.identity = Identity::new(
            derive(EntityKind::Body, Role::Body, Vec::new()),
            derivations,
            BTreeMap::new(),
        )?;
        topology.measure_enclosures(tolerance)?;
        topology.validate(tolerance)?;
        Ok(topology)
    }

    /// A sphere or spherical zone on the frame (S3 of REVIEW_NOTES.md), as
    /// `BRepPrimAPI_MakeSphere(gp_Ax2, radius, low, high)` makes it: the
    /// latitudes `-pi/2 <= low < high <= pi/2` in radians, an end at `+-pi/2`
    /// (the binary64 value) being a pole. Faces are the discs of the ends that
    /// are not poles, then the wall. The wall's loops are the rings (bottom
    /// `+u` at `v = low`, top `-u` at `v = high`) and, when exactly one end
    /// is a pole, that pole as a vertex loop; a whole sphere has no loops.
    /// Entities derive from the meridian as for the cone, the arc being
    /// segment 1 and a pole's vertex having role `Pole`.
    pub(crate) fn sphere(
        frame: Frame3,
        radius: f64,
        low: f64,
        high: f64,
        tolerance: Tolerance,
        operation: OperationId,
    ) -> Result<Self> {
        let tol = tolerance.linear();
        for (value, what) in [
            (radius, "sphere radius"),
            (low, "sphere latitude"),
            (high, "sphere latitude"),
        ] {
            crate::math::finite(value, what)?;
        }
        let half = std::f64::consts::FRAC_PI_2;
        if radius <= tol {
            return Err(Error::Degenerate("sphere radius"));
        }
        if !(-half..=half).contains(&low) || !(-half..=half).contains(&high) || low >= high {
            return Err(Error::OutOfDomain(
                "sphere latitudes must satisfy -pi/2 <= low < high <= pi/2",
            ));
        }
        let ends = [
            (low, low == -half, Role::StartCap, 0, Role::BottomEdge, 1),
            (high, high == half, Role::EndCap, 2, Role::TopEdge, 2),
        ];
        let heights = [radius * low.sin(), radius * high.sin()];
        if heights[1] - heights[0] <= tol {
            return Err(Error::Degenerate("sphere zone height"));
        }
        for (latitude, pole, ..) in ends {
            if !pole && radius * latitude.cos() <= tol {
                return Err(Error::Degenerate("sphere cap radius"));
            }
        }
        let derive = |entity, role, parents| Derivation {
            operation,
            kind: OperationKind::Revolve,
            entity,
            role,
            ordinal: 0,
            parents,
        };
        let meridian = |element| Parent::Profile {
            boundary: 0,
            element,
        };
        let mut topology = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            fins: Vec::new(),
            loops: Vec::new(),
            faces: Vec::new(),
            shells: Vec::new(),
            regions: Vec::new(),
            identity: Identity::new(
                derive(EntityKind::Body, Role::Body, Vec::new()),
                Vec::new(),
                BTreeMap::new(),
            )?,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
        };
        let mut derivations: Vec<(Slot, Derivation)> = vec![(
            Slot::Region(RegionId(1)),
            derive(
                EntityKind::Region,
                Role::Region,
                vec![meridian(ProfileElement::Boundary)],
            ),
        )];
        let poles = ends.iter().filter(|e| e.1).count();
        let mut wall_loops = Vec::new();
        for (k, (latitude, pole, cap_role, segment, edge_role, rim)) in ends.into_iter().enumerate()
        {
            let upper = k == 1;
            let z = heights[k];
            if pole {
                // A pole is a vertex only when it closes a band: a whole
                // sphere has no loops.
                if poles == 1 {
                    let at = topology.add_vertex(
                        frame.point(Point2::default(), if upper { radius } else { -radius }),
                    );
                    derivations.push((
                        Slot::Vertex(at),
                        derive(
                            EntityKind::Vertex,
                            Role::Pole,
                            vec![meridian(ProfileElement::Vertex(rim))],
                        ),
                    ));
                    topology.loops.push(Loop::Vertex(at));
                    wall_loops.push(LoopId(topology.loops.len() - 1));
                }
                continue;
            }
            let ring_radius = radius * latitude.cos();
            let centre = frame.point(Point2::default(), z);
            let normal = if upper {
                frame.normal()
            } else {
                -frame.normal()
            };
            let disc_frame = Frame3::new(centre, normal, frame.x(), tolerance)?;
            let ring_frame = Frame3::new(centre, frame.normal(), frame.x(), tolerance)?;
            let ring = topology.add_ring(Curve3::Circle {
                frame: ring_frame,
                radius: ring_radius,
            });
            derivations.push((
                Slot::Edge(ring),
                derive(
                    EntityKind::Edge,
                    edge_role,
                    vec![meridian(ProfileElement::Vertex(rim))],
                ),
            ));
            let disc = FaceId(topology.faces.len());
            derivations.push((
                Slot::Face(disc),
                derive(
                    EntityKind::Face,
                    cap_role,
                    vec![meridian(ProfileElement::Segment(segment))],
                ),
            ));
            topology.faces.push(Face {
                surface: Surface::Plane(disc_frame),
                sense: Orientation::Forward,
                loops: Vec::new(),
                front: ShellId(0),
                back: ShellId(1),
                enclosure: None,
            });
            topology.add_cap_loop(disc.0, &[ring], !upper, disc_frame);
            let (sense, u0, u1, turns) = if upper {
                (Orientation::Reversed, TAU, 0.0, -1)
            } else {
                (Orientation::Forward, 0.0, TAU, 1)
            };
            let fin = Fin {
                edge: ring,
                sense,
                pcurve: Curve2::LineSegment {
                    start: Point2::new(u0, latitude),
                    end: Point2::new(u1, latitude),
                },
                enclosure: None,
            };
            wall_loops.push(topology.add_loop(vec![fin], [turns, 0]));
        }
        let wall = FaceId(topology.faces.len());
        derivations.push((
            Slot::Face(wall),
            derive(
                EntityKind::Face,
                Role::Wall,
                vec![meridian(ProfileElement::Segment(1))],
            ),
        ));
        topology.faces.push(Face {
            surface: Surface::Sphere { frame, radius },
            sense: Orientation::Forward,
            loops: wall_loops,
            front: ShellId(0),
            back: ShellId(1),
            enclosure: None,
        });
        for (k, fin) in topology.fins.iter().enumerate() {
            topology.edges[fin.edge.index()].fins.push(FinId(k));
        }
        let fronts = topology.face_ids().map(|f| (f, Side::Front)).collect();
        let backs = topology.face_ids().map(|f| (f, Side::Back)).collect();
        topology.shells = vec![
            Shell {
                region: RegionId(1),
                sides: fronts,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
            Shell {
                region: RegionId(0),
                sides: backs,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
        ];
        topology.regions = vec![
            Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(1)],
            },
            Region {
                kind: RegionKind::Solid,
                shells: vec![ShellId(0)],
            },
        ];
        topology.identity = Identity::new(
            derive(EntityKind::Body, Role::Body, Vec::new()),
            derivations,
            BTreeMap::new(),
        )?;
        topology.measure_enclosures(tolerance)?;
        topology.validate(tolerance)?;
        Ok(topology)
    }

    /// A torus, v-segment or wedge on the frame (S3 of REVIEW_NOTES.md), as
    /// `BRepPrimAPI_MakeTorus(gp_Ax2, major, minor, low, high, angle)` makes
    /// it: the tube's latitudes `low..high` revolved by `angle`. A whole torus
    /// (`high - low` and `angle` both the binary64 2 pi) is one face without
    /// loops. A v-segment (a full turn) has rings at `low` and `high`
    /// bounding discs, the wall wound in u and reversed when its meridian arc
    /// bulges toward the axis (the solid is the revolved region between the
    /// arc and the axis). A wedge (the whole tube from `v = 0`, `0 < angle <
    /// 2 pi`) has the tube's circles at `u = 0` and `u = angle` bounding
    /// discs, the wall wound in v. Ids follow the meridian as for the cone;
    /// a wedge's discs come from the boundary and its circles from the arc
    /// (the start and end copies of segment 1, as a prism's cap edges).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn torus(
        frame: Frame3,
        major: f64,
        minor: f64,
        low: f64,
        high: f64,
        angle: f64,
        tolerance: Tolerance,
        operation: OperationId,
    ) -> Result<Self> {
        let tol = tolerance.linear();
        for (value, what) in [
            (major, "torus radius"),
            (minor, "torus radius"),
            (low, "torus latitude"),
            (high, "torus latitude"),
            (angle, "torus angle"),
        ] {
            crate::math::finite(value, what)?;
        }
        if minor <= tol || major - minor <= tol {
            return Err(Error::Degenerate("torus radii (a ring torus)"));
        }
        let closed = high - low == TAU;
        let turn = angle == TAU;
        if !(low < high && high - low <= TAU && angle > 0.0 && angle <= TAU) {
            return Err(Error::OutOfDomain(
                "torus latitudes must satisfy low < high <= low + 2 pi, angle in (0, 2 pi]",
            ));
        }
        if !closed && !turn {
            return Err(Error::OutOfDomain("a torus segment of a partial turn"));
        }
        if !turn && low != 0.0 {
            return Err(Error::OutOfDomain("a torus wedge starts its tube at v = 0"));
        }
        let derive = |entity, role, parents| Derivation {
            operation,
            kind: OperationKind::Revolve,
            entity,
            role,
            ordinal: 0,
            parents,
        };
        let meridian = |element| Parent::Profile {
            boundary: 0,
            element,
        };
        let mut topology = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            fins: Vec::new(),
            loops: Vec::new(),
            faces: Vec::new(),
            shells: Vec::new(),
            regions: Vec::new(),
            identity: Identity::new(
                derive(EntityKind::Body, Role::Body, Vec::new()),
                Vec::new(),
                BTreeMap::new(),
            )?,
            layout: Vec::new(),
            attributes: AttributeMap::new(),
        };
        let mut derivations: Vec<(Slot, Derivation)> = vec![(
            Slot::Region(RegionId(1)),
            derive(
                EntityKind::Region,
                Role::Region,
                vec![meridian(ProfileElement::Boundary)],
            ),
        )];
        let mut wall_loops = Vec::new();
        let mut wall_sense = Orientation::Forward;
        if !closed {
            // A v-segment: the meridian region between the arc and the axis
            // runs counterclockwise (the arc bulges away from the axis, the
            // wall's normal points out) exactly when its signed area is
            // positive; its boundary must be simple, the arc's lower end
            // below its upper one in that direction.
            let area = minor * major * (high.sin() - low.sin())
                + minor
                    * minor
                    * ((high - low) / 2.0 + ((2.0 * high).sin() - (2.0 * low).sin()) / 4.0);
            let (z_low, z_high) = (minor * low.sin(), minor * high.sin());
            if area.abs() <= tol * tol || (z_high - z_low).abs() <= tol {
                return Err(Error::Degenerate("torus segment"));
            }
            if (area > 0.0) != (z_high > z_low) {
                return Err(Error::OutOfDomain(
                    "a torus segment whose meridian boundary is not simple",
                ));
            }
            wall_sense = if area > 0.0 {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            let ends = [
                (low, z_low, Role::StartCap, 0, Role::BottomEdge, 1, false),
                (high, z_high, Role::EndCap, 2, Role::TopEdge, 2, true),
            ];
            for (latitude, z, cap_role, segment, edge_role, rim, upper_end) in ends {
                let ring_radius = major + minor * latitude.cos();
                let centre = frame.point(Point2::default(), z);
                let lower = if upper_end { z < z_low } else { z < z_high };
                let normal = if lower {
                    -frame.normal()
                } else {
                    frame.normal()
                };
                let disc_frame = Frame3::new(centre, normal, frame.x(), tolerance)?;
                let ring_frame = Frame3::new(centre, frame.normal(), frame.x(), tolerance)?;
                let ring = topology.add_ring(Curve3::Circle {
                    frame: ring_frame,
                    radius: ring_radius,
                });
                derivations.push((
                    Slot::Edge(ring),
                    derive(
                        EntityKind::Edge,
                        edge_role,
                        vec![meridian(ProfileElement::Vertex(rim))],
                    ),
                ));
                let disc = FaceId(topology.faces.len());
                derivations.push((
                    Slot::Face(disc),
                    derive(
                        EntityKind::Face,
                        cap_role,
                        vec![meridian(ProfileElement::Segment(segment))],
                    ),
                ));
                topology.faces.push(Face {
                    surface: Surface::Plane(disc_frame),
                    sense: Orientation::Forward,
                    loops: Vec::new(),
                    front: ShellId(0),
                    back: ShellId(1),
                    enclosure: None,
                });
                topology.add_cap_loop(disc.0, &[ring], lower, disc_frame);
                // In the oriented wall: +u at low and -u at high for a
                // forward wall, the other way round for a reversed one.
                let plus = upper_end == (wall_sense == Orientation::Reversed);
                let (sense, u0, u1, turns) = if plus {
                    (Orientation::Forward, 0.0, TAU, 1)
                } else {
                    (Orientation::Reversed, TAU, 0.0, -1)
                };
                let fin = Fin {
                    edge: ring,
                    sense,
                    pcurve: Curve2::LineSegment {
                        start: Point2::new(u0, latitude),
                        end: Point2::new(u1, latitude),
                    },
                    enclosure: None,
                };
                wall_loops.push(topology.add_loop(vec![fin], [turns, 0]));
            }
        } else if !turn {
            // A wedge: the tube's circles at u = 0 and u = angle, each
            // bounding a disc in its meridian half-plane.
            let e = |u: f64| frame.x() * u.cos() + frame.normal().cross(frame.x()) * u.sin();
            let ends = [
                (0.0, Role::StartCap, Role::BottomEdge, false),
                (angle, Role::EndCap, Role::TopEdge, true),
            ];
            for (u, cap_role, edge_role, end) in ends {
                let radial = e(u);
                let centre = frame.origin() + radial * major;
                // The circle's parameter is v: x-axis e(u), y-axis n.
                let circle_normal = radial.cross(frame.normal());
                let ring_frame = Frame3::new(centre, circle_normal, radial, tolerance)?;
                let ring = topology.add_ring(Curve3::Circle {
                    frame: ring_frame,
                    radius: minor,
                });
                derivations.push((
                    Slot::Edge(ring),
                    derive(
                        EntityKind::Edge,
                        edge_role,
                        vec![meridian(ProfileElement::Segment(1))],
                    ),
                ));
                // Outward: back along -u at the start, on along +u at the end.
                let tangent = frame.normal().cross(radial);
                let outward = if end { tangent } else { -tangent };
                let disc_frame = Frame3::new(centre, outward, radial, tolerance)?;
                let disc = FaceId(topology.faces.len());
                derivations.push((
                    Slot::Face(disc),
                    derive(
                        EntityKind::Face,
                        cap_role,
                        vec![meridian(ProfileElement::Boundary)],
                    ),
                ));
                topology.faces.push(Face {
                    surface: Surface::Plane(disc_frame),
                    sense: Orientation::Forward,
                    loops: Vec::new(),
                    front: ShellId(0),
                    back: ShellId(1),
                    enclosure: None,
                });
                topology.add_cap_loop(disc.0, &[ring], end, disc_frame);
                // The region [0, angle] x [0, 2 pi] of the wall: up in v at
                // u = angle, down at u = 0.
                let (sense, v0, v1, turns) = if end {
                    (Orientation::Forward, 0.0, TAU, 1)
                } else {
                    (Orientation::Reversed, TAU, 0.0, -1)
                };
                let fin = Fin {
                    edge: ring,
                    sense,
                    pcurve: Curve2::LineSegment {
                        start: Point2::new(u, v0),
                        end: Point2::new(u, v1),
                    },
                    enclosure: None,
                };
                wall_loops.push(topology.add_loop(vec![fin], [0, turns]));
            }
        }
        let wall = FaceId(topology.faces.len());
        derivations.push((
            Slot::Face(wall),
            derive(
                EntityKind::Face,
                Role::Wall,
                vec![meridian(ProfileElement::Segment(1))],
            ),
        ));
        topology.faces.push(Face {
            surface: Surface::Torus {
                frame,
                major,
                minor,
            },
            sense: wall_sense,
            loops: wall_loops,
            front: ShellId(0),
            back: ShellId(1),
            enclosure: None,
        });
        for (k, fin) in topology.fins.iter().enumerate() {
            topology.edges[fin.edge.index()].fins.push(FinId(k));
        }
        let fronts = topology.face_ids().map(|f| (f, Side::Front)).collect();
        let backs = topology.face_ids().map(|f| (f, Side::Back)).collect();
        topology.shells = vec![
            Shell {
                region: RegionId(1),
                sides: fronts,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
            Shell {
                region: RegionId(0),
                sides: backs,
                wire_edges: Vec::new(),
                acorn_vertices: Vec::new(),
            },
        ];
        topology.regions = vec![
            Region {
                kind: RegionKind::Void,
                shells: vec![ShellId(1)],
            },
            Region {
                kind: RegionKind::Solid,
                shells: vec![ShellId(0)],
            },
        ];
        topology.identity = Identity::new(
            derive(EntityKind::Body, Role::Body, Vec::new()),
            derivations,
            BTreeMap::new(),
        )?;
        topology.measure_enclosures(tolerance)?;
        topology.validate(tolerance)?;
        Ok(topology)
    }

    /// A builder's own enclosures: every gap of the constructed geometry,
    /// measured. A construction that cannot be enclosed within the resolution
    /// fails; nothing widens to succeed (Contract 5, T4).
    fn measure_enclosures(&mut self, tolerance: Tolerance) -> Result<()> {
        let m = validate::measure(&self.view());
        let bounds = m.vertices.iter().chain(&m.fins).chain(&m.faces);
        for bound in bounds {
            match bound {
                None => return Err(Error::Unrepresentable("an entity's gap enclosure")),
                Some(b) if *b > tolerance.linear() => return Err(Error::PrecisionLoss),
                Some(_) => {}
            }
        }
        fill(&mut self.vertices, &m.vertices, |v| &mut v.enclosure);
        fill(&mut self.fins, &m.fins, |f| &mut f.enclosure);
        fill(&mut self.faces, &m.faces, |f| &mut f.enclosure);
        Ok(())
    }

    fn add_vertex(&mut self, position: Point3) -> VertexId {
        let id = VertexId(self.vertices.len());
        self.vertices.push(Vertex {
            position,
            enclosure: None,
        });
        id
    }
    fn add_edge(
        &mut self,
        start: Option<VertexId>,
        end: Option<VertexId>,
        curve: Curve3,
    ) -> EdgeId {
        let id = EdgeId(self.edges.len());
        self.edges.push(Edge {
            start,
            end,
            curve,
            fins: Vec::new(),
        });
        id
    }
    fn add_ring(&mut self, curve: Curve3) -> EdgeId {
        self.add_edge(None, None, curve)
    }
    fn add_line(&mut self, start: VertexId, end: VertexId) -> EdgeId {
        self.add_edge(
            Some(start),
            Some(end),
            Curve3::LineSegment {
                start: self.vertices[start.0].position,
                end: self.vertices[end.0].position,
            },
        )
    }
    fn add_loop(&mut self, fins: Vec<Fin>, winding: [i32; 2]) -> LoopId {
        let ids = fins
            .into_iter()
            .map(|fin| {
                self.fins.push(fin);
                FinId(self.fins.len() - 1)
            })
            .collect();
        self.loops.push(Loop::Edges { fins: ids, winding });
        LoopId(self.loops.len() - 1)
    }
    fn add_cap_loop(&mut self, face: usize, edges: &[EdgeId], reverse: bool, frame: Frame3) {
        let sense = if reverse {
            Orientation::Reversed
        } else {
            Orientation::Forward
        };
        let mut fins = edges
            .iter()
            .map(|edge| self.plane_fin(*edge, sense, frame))
            .collect::<Vec<_>>();
        if reverse {
            fins.reverse();
        }
        let l = self.add_loop(fins, [0, 0]);
        self.faces[face].loops.push(l);
    }
    fn plane_fin(&self, id: EdgeId, sense: Orientation, frame: Frame3) -> Fin {
        Fin {
            edge: id,
            sense,
            pcurve: plane_pcurve(&self.edges[id.0].curve, sense, frame),
            enclosure: None,
        }
    }
}

/// A spline profile segment lifted to `height` on `frame` (S8b): its poles
/// placed by the frame, its knots unchanged, run as the segment runs.
pub(crate) fn lifted_spline(
    span: &SplineSpan<crate::BSplineCurve2>,
    frame: Frame3,
    height: f64,
) -> Result<Curve3> {
    let c = span.curve().as_curve3();
    let poles = c
        .poles()
        .iter()
        .map(|p| frame.point(Point2::new(p.x, p.y), height))
        .collect();
    let lifted = crate::BSplineCurve3::new(
        c.degree(),
        poles,
        None,
        c.knots().to_vec(),
        c.multiplicities().to_vec(),
    )?;
    let whole = SplineSpan::whole(lifted);
    Ok(Curve3::BSpline(if span.is_reversed() {
        whole.reversed()
    } else {
        whole
    }))
}

/// A spline segment's wall (S8b): the degree-(p, 1) surface over its knots
/// and `[0, high - low]`, its poles the profile's lifted to `low` and `high`;
/// with the parameters `u` at the segment's start and end.
pub(crate) fn spline_wall(
    span: &SplineSpan<crate::BSplineCurve2>,
    frame: Frame3,
    low: f64,
    high: f64,
) -> Result<(Surface, f64, f64)> {
    let c = span.curve().as_curve3();
    let mut poles = Vec::with_capacity(2 * c.poles().len());
    for p in c.poles() {
        poles.push(frame.point(Point2::new(p.x, p.y), low));
        poles.push(frame.point(Point2::new(p.x, p.y), high));
    }
    let u = c.knot_vector().clone();
    let v = crate::KnotVector::new(1, vec![0.0, high - low], vec![2, 2])?;
    let surface = crate::BSplineSurface3::new(u, v, poles, None)?;
    let (first, last) = c.domain();
    let (a, b) = if span.is_reversed() {
        (last, first)
    } else {
        (first, last)
    };
    Ok((Surface::BSpline(surface), a, b))
}

/// A planar curve's pcurve on the plane of `frame`, in the use's direction:
/// lines by their ends, arcs and circles turned with the two normals,
/// ellipses sharing the plane's x axis (S8a.2).
pub(crate) fn plane_pcurve(curve: &Curve3, sense: Orientation, frame: Frame3) -> Curve2 {
    let local = |point: Point3| {
        let [x, y, _] = frame.coordinates(point);
        Point2::new(x, y)
    };
    match curve {
        Curve3::LineSegment { start, end } => {
            let (a, b) = if sense == Orientation::Forward {
                (*start, *end)
            } else {
                (*end, *start)
            };
            Curve2::LineSegment {
                start: local(a),
                end: local(b),
            }
        }
        Curve3::Circle {
            frame: circle,
            radius,
        } => {
            let [x, y, _] = frame.coordinates(circle.origin());
            let start = local(circle.point(Point2::new(*radius, 0.0), 0.0));
            Curve2::CircularArc {
                center: Point2::new(x, y),
                radius: *radius,
                start_angle: (start.y - y).atan2(start.x - x),
                sweep_angle: TAU * sense.sign() * circle.normal().dot(frame.normal()).signum(),
            }
        }
        Curve3::CircularArc {
            frame: arc,
            radius,
            sweep_angle,
            ..
        } => {
            // S5: the arc in the cap's frame, its sweep's sign by the two
            // normals and the traversal.
            let [x, y, _] = frame.coordinates(arc.origin());
            let (a, b) = (local(curve.point(0.0)), local(curve.point(1.0)));
            let turn = sweep_angle * arc.normal().dot(frame.normal()).signum();
            let (from, sweep) = if sense == Orientation::Forward {
                (a, turn)
            } else {
                (b, -turn)
            };
            Curve2::CircularArc {
                center: Point2::new(x, y),
                radius: *radius,
                start_angle: (from.y - y).atan2(from.x - x),
                sweep_angle: sweep,
            }
        }
        // S8a.2: an ellipse whose frame shares the plane's x axis (the
        // cut face's), its angle turned with the normals.
        Curve3::EllipseArc {
            frame: ellipse,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => {
            let [x, y, _] = frame.coordinates(ellipse.origin());
            let turn = ellipse.normal().dot(frame.normal()).signum();
            let (start, sweep) = if sense == Orientation::Forward {
                (*start_angle, *sweep_angle)
            } else {
                (start_angle + sweep_angle, -sweep_angle)
            };
            Curve2::EllipseArc {
                center: Point2::new(x, y),
                major: *major,
                minor: *minor,
                start_angle: turn * start,
                sweep_angle: turn * sweep,
            }
        }
        // A hyperbola, a parabola or a torus section in the plane: its exact
        // projection (S8d.2, S8d.3), lifted from its start.
        Curve3::HyperbolaArc { .. } | Curve3::ParabolaArc { .. } | Curve3::Section(_) => {
            let reversed = sense == Orientation::Reversed;
            let start = local(curve.point(if reversed { 1.0 } else { 0.0 }));
            Curve2::Projection(Box::new(
                Projection::new(curve.clone(), Surface::Plane(frame), reversed, start, 8)
                    .expect("a conic or section projects onto a plane"),
            ))
        }
        // S8b: a spline in the plane, its poles in the frame's coordinates
        // (the curve's image by the affine map), run as the use runs.
        Curve3::BSpline(span) => {
            let c = span.curve();
            let poles = c
                .poles()
                .iter()
                .map(|p| {
                    let [x, y, _] = frame.coordinates(*p);
                    Point2::new(x, y)
                })
                .collect();
            let weights = c.is_rational().then(|| c.weights().to_vec());
            let planar = crate::BSplineCurve2::new(
                c.degree(),
                poles,
                weights,
                c.knots().to_vec(),
                c.multiplicities().to_vec(),
            )
            .expect("a spline's image is a spline");
            let [first, last] = span.range();
            let within = SplineSpan::new(planar, first, last).expect("the same range");
            let flip = span.is_reversed() != (sense == Orientation::Reversed);
            Curve2::BSpline(if flip { within.reversed() } else { within })
        }
    }
}
