//! Value ids on structure-aware extrusions (M1 of IDENTITY_AND_HISTORY.md).
//! An independent byte encoder and FNV-1a-128 recompute every id from its
//! retained derivation. Rebuilding, rigid motion, reversing the direction and
//! moving labelled points keep every id; permuting label values permutes ids
//! bijectively, parent by parent; counts and roles follow the profile. The
//! same bytes also make a cone (S3 of REVIEW_NOTES.md): its entities follow
//! the meridian (identity_reference.cone_entities), and rebuilding,
//! stretching and rigid motion keep its ids. So do a sphere's or zone's, and
//! a torus's, v-segment's or wedge's (identity_reference.torus_entities).
//! And a filleted rectangle with an optional notch (S5): arcs exactly on
//! their circles, valid, a polygon's counts, ids kept through reversal and
//! rigid motion, the closed-form area and mass inside the certified
//! enclosure, and points classified by the arcs. The profile's face body
//! and its outer boundary's wire body (S6) derive from the same elements as
//! the prism's start cap, bottom edges and bottom vertices, validate as a
//! sheet and a wire with OCCT's counts, keep their ids through rigid motion,
//! and enclose the profile's area.

use libfuzzer_sys::arbitrary::{Result, Unstructured};
use rusty_occt::history::History;
use rusty_occt::identity::{
    Derivation, EntityKind, InputLabel, OperationId, OperationKind, Parent, ProfileElement, Role,
};
use rusty_occt::topology::Slot;
use rusty_occt::{
    Body, Boundary, BoundaryLabels, Frame3, Location, Point2, Point3, Profile, RigidTransform,
    Segment, Solid, Tolerance, Vec3,
};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::TAU;

/// Independent version-1 encoder (see identity.rs for the specification).
fn encode(d: &Derivation) -> Vec<u8> {
    let mut out = b"RSID\x01".to_vec();
    out.extend(d.operation.0.to_le_bytes());
    out.push(match d.kind {
        OperationKind::Extrude => 1,
        OperationKind::Transform => 2,
        OperationKind::External => 3,
        OperationKind::Composite => 4,
        OperationKind::HeightSplit => 5,
        OperationKind::StackedFuse => 6,
        OperationKind::Revolve => 7,
        OperationKind::MakeFace => 8,
        OperationKind::MakeWire => 9,
        OperationKind::PlaneSplit => 10,
    });
    out.push(match d.entity {
        EntityKind::Vertex => 1,
        EntityKind::Edge => 2,
        EntityKind::Face => 3,
        EntityKind::Body => 4,
        EntityKind::Region => 5,
    });
    let roles = [
        Role::StartCap,
        Role::EndCap,
        Role::Wall,
        Role::BottomEdge,
        Role::TopEdge,
        Role::Vertical,
        Role::Seam,
        Role::BottomVertex,
        Role::TopVertex,
        Role::SeamVertex,
        Role::Body,
        Role::External,
        Role::Region,
        Role::CutFace,
        Role::CutEdge,
        Role::CutVertex,
        Role::Apex,
        Role::Pole,
        Role::Face,
        Role::Edge,
        Role::Vertex,
    ];
    out.push(roles.iter().position(|r| *r == d.role).unwrap() as u8 + 1);
    out.extend(d.ordinal.to_le_bytes());
    out.extend((d.parents.len() as u32).to_le_bytes());
    for p in &d.parents {
        match p {
            Parent::Label(l) => {
                out.push(1);
                out.extend(l.0.to_le_bytes());
            }
            Parent::Profile { boundary, element } => {
                out.push(2);
                out.extend(boundary.to_le_bytes());
                let (tag, i) = match element {
                    ProfileElement::Boundary => (0, 0u32),
                    ProfileElement::Segment(i) => (1, *i),
                    ProfileElement::Vertex(i) => (2, *i),
                };
                out.push(tag);
                out.extend(i.to_le_bytes());
            }
            Parent::Entity(id) => {
                out.push(3);
                out.extend(id.0);
            }
        }
    }
    out
}

fn fnv(bytes: &[u8]) -> [u8; 16] {
    let mut h: u128 = 0x6c62272e07bb014262b821756295c58d;
    for b in bytes {
        h ^= u128::from(*b);
        h = h.wrapping_mul(0x0000000001000000000000000000013B);
    }
    h.to_be_bytes()
}

pub(crate) struct Spec {
    pub(crate) tolerance: Tolerance,
    pub(crate) operation: OperationId,
    frame: Frame3,
    start: f64,
    end: f64,
    /// Outer polygon radii (None: a circle) and holes (square or circle).
    outer: Option<Vec<f64>>,
    holes: Vec<(bool, Point2)>,
    clockwise: bool,
    pub(crate) labels: Option<Vec<u64>>,
    pub(crate) transforms: Vec<RigidTransform>,
}

fn unit(u: &mut Unstructured) -> Result<f64> {
    Ok(f64::from(u.arbitrary::<u16>()?) / 65535.0)
}

pub(crate) fn spec(u: &mut Unstructured) -> Result<Option<Spec>> {
    let scale = 2f64.powi(u.int_in_range(-8..=8)?);
    let tolerance = Tolerance::new(1e-9 * scale, 1e-12).unwrap();
    let outer = if u.ratio(1, 6)? {
        None
    } else {
        let n = u.int_in_range(3..=12)?;
        Some(
            (0..n)
                .map(|_| unit(u).map(|r| 0.8 + 0.4 * r))
                .collect::<Result<Vec<_>>>()?,
        )
    };
    let mut holes = Vec::new();
    for k in 0..u.int_in_range(0..=3)? {
        let angle = TAU * (k as f64 + 0.5 * unit(u)?) / 3.0;
        holes.push((
            u.arbitrary()?,
            Point2::new(0.4 * angle.cos(), 0.4 * angle.sin()),
        ));
    }
    let normal = Vec3::new(2.0 * unit(u)? - 1.0, 2.0 * unit(u)? - 1.0, 0.3 + unit(u)?);
    let origin = Point3::new(
        scale * (10.0 * unit(u)? - 5.0),
        scale * 10.0 * unit(u)?,
        scale * -3.0 * unit(u)?,
    );
    let Ok(frame) = Frame3::new(origin, normal, Vec3::X, tolerance) else {
        return Ok(None);
    };
    let (a, b) = (-scale * (0.2 + unit(u)?), scale * (0.2 + unit(u)?));
    let (start, end) = if u.arbitrary()? { (a, b) } else { (b, a) };
    let labels = if u.arbitrary()? {
        // Distinct values: an odd multiplier is invertible modulo 2^64.
        let base: u64 = u.arbitrary()?;
        let step: u64 = u.arbitrary::<u64>()? | 1;
        Some(
            (0..64u64)
                .map(|k| base.wrapping_add(k.wrapping_mul(step)))
                .collect(),
        )
    } else {
        None
    };
    let mut transforms = Vec::new();
    for _ in 0..u.int_in_range(0..=2)? {
        let axis = Vec3::new(unit(u)? + 0.1, unit(u)? - 0.5, unit(u)? - 0.5);
        let r = RigidTransform::rotation(Point3::ORIGIN, axis, TAU * unit(u)?).unwrap();
        let t =
            RigidTransform::translation(Vec3::new(unit(u)?, unit(u)?, unit(u)?) * (7.0 * scale))
                .unwrap();
        transforms.push(r.then(t).unwrap());
    }
    Ok(Some(Spec {
        tolerance,
        operation: OperationId(u.arbitrary()?),
        frame,
        start,
        end,
        outer,
        holes,
        clockwise: u.arbitrary()?,
        labels,
        transforms,
    }))
}

/// Build with labels drawn from `pool` in order; `stretch` moves every point.
pub(crate) fn build(s: &Spec, pool: Option<&[u64]>, stretch: f64, reverse: bool) -> Option<Solid> {
    build_tracked(s, pool, stretch, reverse).map(|(solid, _)| solid)
}

/// The extrusion and its construction history.
pub(crate) fn build_tracked(
    s: &Spec,
    pool: Option<&[u64]>,
    stretch: f64,
    reverse: bool,
) -> Option<(Solid, History)> {
    let scale = s.tolerance.linear() / 1e-9;
    let mut next = 0;
    let mut take = |n: usize| -> Vec<InputLabel> {
        let out = (next..next + n)
            .map(|k| InputLabel(pool.unwrap()[k]))
            .collect();
        next += n;
        out
    };
    let mut label = |b: Boundary, n: usize| -> Boundary {
        if pool.is_none() {
            return b;
        }
        let boundary = take(1)[0];
        let segments = take(n);
        let vertices = take(n);
        b.with_labels(BoundaryLabels {
            boundary,
            segments,
            vertices,
        })
        .unwrap()
    };
    let outer = match &s.outer {
        None => label(
            Boundary::circle(Point2::default(), scale * stretch, s.tolerance).ok()?,
            1,
        ),
        Some(radii) => {
            let n = radii.len();
            let mut points: Vec<Point2> = radii
                .iter()
                .enumerate()
                .map(|(i, r)| {
                    let a = TAU * i as f64 / n as f64;
                    Point2::new(scale * stretch * r * a.cos(), scale * stretch * r * a.sin())
                })
                .collect();
            if s.clockwise {
                points.reverse();
            }
            label(Boundary::polygon(points, s.tolerance).ok()?, n)
        }
    };
    let mut holes = Vec::new();
    for (square, c) in &s.holes {
        let (c, h) = (Point2::new(c.x * scale, c.y * scale), 0.12 * scale);
        let hole = if *square {
            let pts = vec![
                Point2::new(c.x - h, c.y - h),
                Point2::new(c.x + h, c.y - h),
                Point2::new(c.x + h, c.y + h),
                Point2::new(c.x - h, c.y + h),
            ];
            label(Boundary::polygon(pts, s.tolerance).ok()?, 4)
        } else {
            label(Boundary::circle(c, h, s.tolerance).ok()?, 1)
        };
        holes.push(hole);
    }
    let profile = Profile::new(outer, holes, s.tolerance).ok()?;
    let (start, end) = if reverse {
        (s.end, s.start)
    } else {
        (s.start, s.end)
    };
    Solid::extrude_with(s.operation, profile, s.frame, start, end).ok()
}

type Ids = BTreeMap<Slot, (String, Derivation)>;

fn ids(solid: &Solid) -> Ids {
    let t = solid.topology();
    let mut out = BTreeMap::new();
    for (id, slot) in t.ids() {
        let d = t.derivation(id).unwrap().clone();
        // The independent encoder and digest reproduce the id.
        assert_eq!(fnv(&encode(&d)), id.0, "{d:?}");
        assert_eq!(d.encode(), encode(&d));
        assert_eq!(t.id_of(slot), Some(id));
        assert_eq!(t.slot_of(id), Some(slot));
        assert_eq!(d.operation, solid.operation());
        out.insert(slot, (id.to_string(), d));
    }
    // Every bounded region has an id; the infinite void is no entity.
    let total = t.vertices().len() + t.edges().len() + t.faces().len() + t.regions().len() - 1;
    assert_eq!(out.len(), total, "every slot has an id");
    let distinct: BTreeSet<&String> = out.values().map(|v| &v.0).collect();
    assert_eq!(distinct.len(), total, "ids are unique");
    assert!(!distinct.contains(&t.body_id().to_string()));
    out
}

fn id_set(ids: &Ids) -> BTreeSet<String> {
    ids.values().map(|v| v.0.clone()).collect()
}

/// A cone of `Solid::cone_with`, with rigid motions.
pub(crate) struct ConeSpec {
    pub(crate) tolerance: Tolerance,
    pub(crate) operation: OperationId,
    pub(crate) frame: Frame3,
    pub(crate) bottom: f64,
    pub(crate) top: f64,
    pub(crate) height: f64,
    pub(crate) transforms: Vec<RigidTransform>,
}

impl ConeSpec {
    pub(crate) fn build(&self, stretch: f64) -> Option<(Solid, History)> {
        Solid::cone_with(
            self.operation,
            self.frame,
            self.bottom * stretch,
            self.top * stretch,
            self.height * stretch,
            self.tolerance,
        )
        .ok()
    }
}

pub(crate) fn cone_spec(u: &mut Unstructured) -> Result<Option<ConeSpec>> {
    let scale = 2f64.powi(u.int_in_range(-8..=8)?);
    let tolerance = Tolerance::new(1e-9 * scale, 1e-12).unwrap();
    let radius = |u: &mut Unstructured| -> Result<f64> {
        Ok(if u.ratio(1, 3)? {
            0.0
        } else {
            scale * (0.1 + unit(u)?)
        })
    };
    let (bottom, top) = (radius(u)?, radius(u)?);
    let height = scale * (0.1 + unit(u)?);
    let normal = Vec3::new(2.0 * unit(u)? - 1.0, 2.0 * unit(u)? - 1.0, 0.3 + unit(u)?);
    let origin = Point3::new(
        scale * (10.0 * unit(u)? - 5.0),
        scale * 10.0 * unit(u)?,
        scale * -3.0 * unit(u)?,
    );
    let Ok(frame) = Frame3::new(origin, normal, Vec3::X, tolerance) else {
        return Ok(None);
    };
    let mut transforms = Vec::new();
    for _ in 0..u.int_in_range(0..=2)? {
        let axis = Vec3::new(unit(u)? + 0.1, unit(u)? - 0.5, unit(u)? - 0.5);
        let r = RigidTransform::rotation(Point3::ORIGIN, axis, TAU * unit(u)?).unwrap();
        let t =
            RigidTransform::translation(Vec3::new(unit(u)?, unit(u)?, unit(u)?) * (7.0 * scale))
                .unwrap();
        transforms.push(r.then(t).unwrap());
    }
    Ok(Some(ConeSpec {
        tolerance,
        operation: OperationId(u.arbitrary()?),
        frame,
        bottom,
        top,
        height,
        transforms,
    }))
}

/// A sphere or zone of `Solid::sphere_with`, with rigid motions.
pub(crate) struct SphereSpec {
    pub(crate) tolerance: Tolerance,
    pub(crate) operation: OperationId,
    pub(crate) frame: Frame3,
    pub(crate) radius: f64,
    pub(crate) low: f64,
    pub(crate) high: f64,
    pub(crate) transforms: Vec<RigidTransform>,
}

impl SphereSpec {
    pub(crate) fn build(&self, stretch: f64) -> Option<(Solid, History)> {
        Solid::sphere_with(
            self.operation,
            self.frame,
            self.radius * stretch,
            self.low,
            self.high,
            self.tolerance,
        )
        .ok()
    }
}

pub(crate) fn sphere_spec(u: &mut Unstructured) -> Result<Option<SphereSpec>> {
    let half = std::f64::consts::FRAC_PI_2;
    let scale = 2f64.powi(u.int_in_range(-8..=8)?);
    let tolerance = Tolerance::new(1e-9 * scale, 1e-12).unwrap();
    let radius = scale * (0.2 + unit(u)?);
    let latitude = |u: &mut Unstructured| -> Result<f64> { Ok(1.4 * unit(u)? - 0.7) };
    let (low, high) = match u.int_in_range(0..=3)? {
        0 => (-half, half),
        1 => (-half, latitude(u)?),
        2 => (latitude(u)?, half),
        _ => {
            let a = latitude(u)?;
            (a, a + 0.1 + 0.7 * unit(u)?)
        }
    };
    let normal = Vec3::new(2.0 * unit(u)? - 1.0, 2.0 * unit(u)? - 1.0, 0.3 + unit(u)?);
    let origin = Point3::new(
        scale * (10.0 * unit(u)? - 5.0),
        scale * 10.0 * unit(u)?,
        scale * -3.0 * unit(u)?,
    );
    let Ok(frame) = Frame3::new(origin, normal, Vec3::X, tolerance) else {
        return Ok(None);
    };
    let mut transforms = Vec::new();
    for _ in 0..u.int_in_range(0..=2)? {
        let axis = Vec3::new(unit(u)? + 0.1, unit(u)? - 0.5, unit(u)? - 0.5);
        let r = RigidTransform::rotation(Point3::ORIGIN, axis, TAU * unit(u)?).unwrap();
        let t =
            RigidTransform::translation(Vec3::new(unit(u)?, unit(u)?, unit(u)?) * (7.0 * scale))
                .unwrap();
        transforms.push(r.then(t).unwrap());
    }
    Ok(Some(SphereSpec {
        tolerance,
        operation: OperationId(u.arbitrary()?),
        frame,
        radius,
        low,
        high,
        transforms,
    }))
}

fn check_sphere(data: &[u8]) {
    let half = std::f64::consts::FRAC_PI_2;
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = sphere_spec(&mut u) else {
        return;
    };
    let (solid, _) = s.build(1.0).expect("a sphere within its domain builds");
    let base = ids(&solid);
    let t = solid.topology();
    let poles = usize::from(s.low == -half) + usize::from(s.high == half);
    let rings = 2 - poles;
    // A pole is a vertex only when it closes a band; the whole sphere has
    // none.
    assert_eq!(t.vertices().len(), usize::from(poles == 1));
    assert_eq!(t.edges().len(), rings);
    assert_eq!(t.faces().len(), 1 + rings);
    let c = t.occt_counts();
    assert_eq!(
        (c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids),
        (2, 3, 1 + rings, 1 + rings, 1, 1)
    );
    for (_, d) in base.values() {
        assert_eq!(d.kind, OperationKind::Revolve);
        let want = match d.role {
            Role::Pole if s.low == -half => ProfileElement::Vertex(1),
            Role::Pole => ProfileElement::Vertex(2),
            role => meridian_parent(role, 1.0),
        };
        assert_eq!(
            d.parents,
            vec![Parent::Profile {
                boundary: 0,
                element: want,
            }]
        );
    }
    assert_eq!(ids(&s.build(1.0).unwrap().0), base);
    if let Some((stretched, _)) = s.build(1.03) {
        assert_eq!(ids(&stretched), base, "stretching keeps ids");
    }
    let mut moved = solid.clone();
    for transform in &s.transforms {
        if let Ok((next, _)) = moved.transform_with(OperationId::UNSPECIFIED, *transform) {
            assert_eq!(ids(&next), base);
            moved = next;
        }
    }
}

/// A torus, v-segment or wedge of `Solid::torus_with`, with rigid motions.
pub(crate) struct TorusSpec {
    pub(crate) tolerance: Tolerance,
    pub(crate) operation: OperationId,
    pub(crate) frame: Frame3,
    pub(crate) major: f64,
    pub(crate) minor: f64,
    pub(crate) low: f64,
    pub(crate) high: f64,
    pub(crate) angle: f64,
    pub(crate) transforms: Vec<RigidTransform>,
}

impl TorusSpec {
    pub(crate) fn build(&self, stretch: f64) -> Option<(Solid, History)> {
        Solid::torus_with(
            self.operation,
            self.frame,
            self.major * stretch,
            self.minor * stretch,
            self.low,
            self.high,
            self.angle,
            self.tolerance,
        )
        .ok()
    }

    pub(crate) fn closed(&self) -> bool {
        self.high - self.low == TAU
    }
}

pub(crate) fn torus_spec(u: &mut Unstructured) -> Result<Option<TorusSpec>> {
    let scale = 2f64.powi(u.int_in_range(-8..=8)?);
    let tolerance = Tolerance::new(1e-9 * scale, 1e-12).unwrap();
    let minor = scale * (0.2 + unit(u)?);
    let major = minor + scale * (0.1 + unit(u)?);
    let (low, high, angle) = match u.int_in_range(0..=2)? {
        0 => (0.0, TAU, TAU),
        1 => (0.0, TAU, 0.2 + 5.9 * unit(u)?),
        _ => {
            let a = TAU * unit(u)? - std::f64::consts::PI;
            (a, a + 0.1 + 5.9 * unit(u)?, TAU)
        }
    };
    let normal = Vec3::new(2.0 * unit(u)? - 1.0, 2.0 * unit(u)? - 1.0, 0.3 + unit(u)?);
    let origin = Point3::new(
        scale * (10.0 * unit(u)? - 5.0),
        scale * 10.0 * unit(u)?,
        scale * -3.0 * unit(u)?,
    );
    let Ok(frame) = Frame3::new(origin, normal, Vec3::X, tolerance) else {
        return Ok(None);
    };
    let mut transforms = Vec::new();
    for _ in 0..u.int_in_range(0..=2)? {
        let axis = Vec3::new(unit(u)? + 0.1, unit(u)? - 0.5, unit(u)? - 0.5);
        let r = RigidTransform::rotation(Point3::ORIGIN, axis, TAU * unit(u)?).unwrap();
        let t =
            RigidTransform::translation(Vec3::new(unit(u)?, unit(u)?, unit(u)?) * (7.0 * scale))
                .unwrap();
        transforms.push(r.then(t).unwrap());
    }
    Ok(Some(TorusSpec {
        tolerance,
        operation: OperationId(u.arbitrary()?),
        frame,
        major,
        minor,
        low,
        high,
        angle,
        transforms,
    }))
}

fn check_torus(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = torus_spec(&mut u) else {
        return;
    };
    let Some((solid, _)) = s.build(1.0) else {
        // Only a segment whose meridian boundary is not simple is refused.
        assert!(!s.closed());
        return;
    };
    let base = ids(&solid);
    let t = solid.topology();
    let whole = s.closed() && s.angle == TAU;
    let sides = if whole { 0 } else { 2 };
    assert_eq!(t.vertices().len(), 0);
    assert_eq!(t.edges().len(), sides);
    assert_eq!(t.faces().len(), 1 + sides);
    let c = t.occt_counts();
    assert_eq!(
        (c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids),
        if whole {
            (1, 2, 1, 1, 1, 1)
        } else {
            (2, 3, 3, 3, 1, 1)
        }
    );
    for (_, d) in base.values() {
        assert_eq!(d.kind, OperationKind::Revolve);
        assert_eq!(d.ordinal, 0);
        // A segment's sides come from its rim points and radial segments,
        // a wedge's from the arc and the whole meridian.
        let want = match (d.role, s.closed()) {
            (Role::Region, _) | (Role::StartCap | Role::EndCap, true) => ProfileElement::Boundary,
            (Role::Wall, _) | (Role::BottomEdge | Role::TopEdge, true) => {
                ProfileElement::Segment(1)
            }
            (role, false) => meridian_parent(role, 1.0),
            (other, true) => panic!("a torus wedge has no {other:?}"),
        };
        assert_eq!(
            d.parents,
            vec![Parent::Profile {
                boundary: 0,
                element: want,
            }]
        );
    }
    assert_eq!(ids(&s.build(1.0).unwrap().0), base);
    if let Some((stretched, _)) = s.build(1.03) {
        assert_eq!(ids(&stretched), base, "stretching keeps ids");
    }
    let mut moved = solid.clone();
    for transform in &s.transforms {
        if let Ok((next, _)) = moved.transform_with(OperationId::UNSPECIFIED, *transform) {
            assert_eq!(ids(&next), base);
            moved = next;
        }
    }
}

/// The meridian parent every cone entity must have, by role: the rim points
/// 1 and 2, the radial segments 0 and 2, the slant 1, the boundary.
fn meridian_parent(role: Role, bottom: f64) -> ProfileElement {
    match role {
        Role::Region => ProfileElement::Boundary,
        Role::Wall => ProfileElement::Segment(1),
        Role::StartCap => ProfileElement::Segment(0),
        Role::EndCap => ProfileElement::Segment(2),
        Role::BottomEdge => ProfileElement::Vertex(1),
        Role::TopEdge => ProfileElement::Vertex(2),
        Role::Apex if bottom == 0.0 => ProfileElement::Vertex(1),
        Role::Apex => ProfileElement::Vertex(2),
        other => panic!("a cone has no {other:?}"),
    }
}

fn check_cone(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = cone_spec(&mut u) else {
        return;
    };
    let Some((solid, _)) = s.build(1.0) else {
        // Equal radii are a cylinder; two apices bound nothing.
        assert_eq!(s.bottom, s.top);
        return;
    };
    let base = ids(&solid);
    let t = solid.topology();
    let apices = usize::from(s.bottom == 0.0) + usize::from(s.top == 0.0);
    let rings = 2 - apices;
    assert_eq!(t.vertices().len(), apices);
    assert_eq!(t.edges().len(), rings);
    assert_eq!(t.faces().len(), 1 + rings);
    assert_eq!(t.regions().len(), 2);
    let c = t.occt_counts();
    assert_eq!(
        (c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids),
        (2, 3, 1 + rings, 1 + rings, 1, 1)
    );
    for (_, d) in base.values() {
        assert_eq!(d.kind, OperationKind::Revolve);
        assert_eq!(d.ordinal, 0);
        assert_eq!(
            d.parents,
            vec![Parent::Profile {
                boundary: 0,
                element: meridian_parent(d.role, s.bottom),
            }]
        );
    }
    // Rebuilding is deterministic; ids ignore the dimensions.
    assert_eq!(ids(&s.build(1.0).unwrap().0), base);
    if let Some((stretched, _)) = s.build(1.03) {
        assert_eq!(ids(&stretched), base, "stretching keeps ids");
    }
    let mut moved = solid.clone();
    for transform in &s.transforms {
        if let Ok((next, _)) = moved.transform_with(OperationId::UNSPECIFIED, *transform) {
            assert_eq!(ids(&next), base);
            assert_eq!(next.topology().body_id(), solid.topology().body_id());
            moved = next;
        }
    }
}

/// S5: a `w x h` rectangle whose corners are filleted (radius 0 for a sharp
/// corner) and whose top side may carry a concave half-circle notch; every
/// value dyadic, so each arc's points lie exactly on its circle.
pub(crate) fn arc_path(u: &mut Unstructured) -> Result<(Vec<Point2>, Vec<Segment>, f64)> {
    let (w, h) = (
        1.0 + f64::from(u.int_in_range(0u8..=63)?) / 16.0,
        1.0 + f64::from(u.int_in_range(0u8..=63)?) / 16.0,
    );
    let small = w.min(h);
    let mut radius = [0.0; 4];
    for r in &mut radius {
        // Up to a quarter of the shorter side, in 64ths.
        *r = small / 4.0 * f64::from(u.int_in_range(0u8..=63)?) / 64.0;
    }
    let notch = if u.arbitrary()? { small / 8.0 } else { 0.0 };
    let (mut points, mut segments) = (Vec::new(), Vec::new());
    let mut area = w * h;
    let corner = |points: &mut Vec<Point2>,
                  segments: &mut Vec<Segment>,
                  at: Point2,
                  into: Point2,
                  out: Point2,
                  r: f64| {
        // The corner `at`, entered along `into` and left along `out` (unit
        // axis directions).
        if r == 0.0 {
            points.push(at);
            segments.push(Segment::Line);
        } else {
            points.push(Point2::new(at.x - into.x * r, at.y - into.y * r));
            segments.push(Segment::Arc {
                center: Point2::new(at.x - into.x * r + out.x * r, at.y - into.y * r + out.y * r),
                radius: r,
                ccw: true,
            });
            points.push(Point2::new(at.x + out.x * r, at.y + out.y * r));
            segments.push(Segment::Line);
        }
    };
    let (e, n, west, s) = (
        Point2::new(1.0, 0.0),
        Point2::new(0.0, 1.0),
        Point2::new(-1.0, 0.0),
        Point2::new(0.0, -1.0),
    );
    corner(
        &mut points,
        &mut segments,
        Point2::new(0.0, 0.0),
        s,
        e,
        radius[0],
    );
    corner(
        &mut points,
        &mut segments,
        Point2::new(w, 0.0),
        e,
        n,
        radius[1],
    );
    corner(
        &mut points,
        &mut segments,
        Point2::new(w, h),
        n,
        west,
        radius[2],
    );
    if notch > 0.0 {
        points.push(Point2::new(w / 2.0 + notch, h));
        segments.push(Segment::Arc {
            center: Point2::new(w / 2.0, h),
            radius: notch,
            ccw: false,
        });
        points.push(Point2::new(w / 2.0 - notch, h));
        segments.push(Segment::Line);
        area -= std::f64::consts::PI * notch * notch / 2.0;
    }
    corner(
        &mut points,
        &mut segments,
        Point2::new(0.0, h),
        west,
        s,
        radius[3],
    );
    for r in radius {
        area -= r * r * (1.0 - std::f64::consts::PI / 4.0);
    }
    Ok((points, segments, area))
}

fn check_arcs(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok((points, segments, area)) = arc_path(&mut u) else {
        return;
    };
    // A path of lines only is the polygon, covered above.
    if segments.iter().all(|s| *s == Segment::Line) {
        return;
    }
    let Ok(Some(s)) = spec(&mut u) else {
        return;
    };
    let scale = s.tolerance.linear() / 1e-9;
    let scaled: Vec<Point2> = points
        .iter()
        .map(|p| Point2::new(p.x * scale, p.y * scale))
        .collect();
    let arcs: Vec<Segment> = segments
        .iter()
        .map(|g| match *g {
            Segment::Arc {
                center,
                radius,
                ccw,
            } => Segment::Arc {
                center: Point2::new(center.x * scale, center.y * scale),
                radius: radius * scale,
                ccw,
            },
            Segment::Line => Segment::Line,
        })
        .collect();
    let n = scaled.len();
    let boundary = Boundary::path(scaled.clone(), arcs.clone(), s.tolerance)
        .expect("a filleted rectangle is a valid path");
    assert!((boundary.area() - area * scale * scale).abs() <= 1e-12 * area * scale * scale);
    // The same path entered clockwise is stored the same way.
    let mut cw_points = vec![scaled[0]];
    cw_points.extend(scaled[1..].iter().rev());
    let cw_segments: Vec<Segment> = (0..n)
        .map(|j| match arcs[n - 1 - j] {
            Segment::Arc {
                center,
                radius,
                ccw,
            } => Segment::Arc {
                center,
                radius,
                ccw: !ccw,
            },
            g => g,
        })
        .collect();
    let cw = Boundary::path(cw_points, cw_segments, s.tolerance).unwrap();
    assert_eq!(cw.path_geometry(), boundary.path_geometry());
    let build = |b: Boundary, reverse: bool| {
        let profile = Profile::new(b, vec![], s.tolerance).unwrap();
        let (start, end) = if reverse {
            (s.end, s.start)
        } else {
            (s.start, s.end)
        };
        Solid::extrude_with(s.operation, profile, s.frame, start, end).map(|(x, _)| x)
    };
    let Ok(solid) = build(boundary.clone(), false) else {
        return;
    };
    let t = solid.topology();
    assert_eq!(t.check(solid.resolution()), Vec::new());
    assert_eq!(t.vertices().len(), 2 * n);
    assert_eq!(t.edges().len(), 3 * n);
    assert_eq!(t.faces().len(), n + 2);
    let c = t.occt_counts();
    assert_eq!(
        (c.vertices, c.edges, c.wires, c.faces),
        (2 * n, 3 * n, n + 2, n + 2)
    );
    let base = ids(&solid);
    assert_eq!(
        ids(&build(cw, false).unwrap()),
        base,
        "clockwise input keeps ids"
    );
    assert_eq!(id_set(&ids(&build(boundary, true).unwrap())), id_set(&base));
    let mut moved = solid.clone();
    for transform in &s.transforms {
        if let Ok((next, _)) = moved.transform_with(OperationId::UNSPECIFIED, *transform) {
            assert_eq!(ids(&next), base);
            moved = next;
        }
    }
    let m = solid.mass_properties();
    let e = t.mass_enclosure().expect("certified");
    let slack = 1e-12 * m.volume.abs().max(1e-300);
    assert!(e.volume[0] - slack <= m.volume && m.volume <= e.volume[1] + slack);
    // The rectangle's filleted corner lies outside; its centre inside.
    let mid = (s.start + s.end) / 2.0;
    let at = |p: Point2| {
        solid
            .frame()
            .point(Point2::new(p.x * scale, p.y * scale), mid)
    };
    let (w, h) = (
        points[1..].iter().fold(0.0f64, |a, p| a.max(p.x)),
        points.iter().fold(0.0f64, |a, p| a.max(p.y)),
    );
    assert_eq!(
        solid.classify(at(Point2::new(w / 2.0, h / 4.0))).unwrap(),
        Location::Inside
    );
    if let Segment::Arc { radius, .. } = segments[0] {
        if radius * scale > 1e-3 * scale {
            assert_eq!(
                solid.classify(at(Point2::new(0.0, 0.0))).unwrap(),
                Location::Outside
            );
        }
    }
}

pub fn check_identity(data: &[u8]) {
    check_arcs(data);
    check_cone(data);
    check_sphere(data);
    check_torus(data);
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = spec(&mut u) else {
        return;
    };
    let pool = s.labels.as_deref();
    let Some(solid) = build(&s, pool, 1.0, false) else {
        return;
    };
    let base = ids(&solid);
    check_bodies(&solid, &base, &s.transforms);
    // Counts and roles follow the profile structure.
    let boundaries = 1 + s.holes.len();
    let polygon_sides: usize =
        s.outer.as_ref().map_or(0, Vec::len) + 4 * s.holes.iter().filter(|h| h.0).count();
    let circles = usize::from(s.outer.is_none()) + s.holes.iter().filter(|h| !h.0).count();
    let t = solid.topology();
    // A circle sweeps two ring edges and one seamless wall, no vertices.
    assert_eq!(t.vertices().len(), 2 * polygon_sides);
    assert_eq!(t.edges().len(), 3 * polygon_sides + 2 * circles);
    assert_eq!(t.faces().len(), 2 + polygon_sides + circles);
    assert_eq!(t.regions().len(), 2);
    // OCCT's encoding adds a seam and a seam vertex per ring edge per circle.
    let c = t.occt_counts();
    assert_eq!(
        (c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids),
        (
            2 * polygon_sides + 2 * circles,
            3 * polygon_sides + 3 * circles,
            2 * boundaries + polygon_sides + circles,
            2 + polygon_sides + circles,
            1,
            1
        )
    );
    let regions: Vec<&Derivation> = base
        .values()
        .map(|v| &v.1)
        .filter(|d| d.role == Role::Region)
        .collect();
    assert_eq!(regions.len(), 1);
    assert_eq!(regions[0].entity, EntityKind::Region);
    assert_eq!(regions[0].parents.len(), boundaries);
    let caps: Vec<&Derivation> = base
        .values()
        .map(|v| &v.1)
        .filter(|d| matches!(d.role, Role::StartCap | Role::EndCap))
        .collect();
    assert_eq!(caps.len(), 2);
    assert!(caps.iter().all(|d| d.parents.len() == boundaries));
    // Rebuilding is deterministic, and ids ignore geometry and direction.
    assert_eq!(ids(&build(&s, pool, 1.0, false).unwrap()), base);
    if let Some(stretched) = build(&s, pool, 1.03, false) {
        assert_eq!(ids(&stretched), base, "moving points keeps ids");
    }
    let reversed = ids(&build(&s, pool, 1.0, true).unwrap());
    assert_eq!(id_set(&reversed), id_set(&base), "reversing keeps ids");
    // Rigid motion keeps every id, slot and derivation.
    let mut moved = solid.clone();
    for transform in &s.transforms {
        if let Ok(next) = moved
            .transform_with(OperationId::UNSPECIFIED, *transform)
            .map(|(s, _)| s)
        {
            assert_eq!(ids(&next), base);
            assert_eq!(next.topology().body_id(), solid.topology().body_id());
            moved = next;
        }
    }
    // Permuting label values permutes ids bijectively, parent by parent.
    if let Some(pool) = pool {
        let mut shuffled = pool.to_vec();
        shuffled.rotate_left(1 + (data.len() % 7));
        let map: BTreeMap<u64, u64> = pool.iter().copied().zip(shuffled.iter().copied()).collect();
        let permuted = ids(&build(&s, Some(&shuffled), 1.0, false).unwrap());
        for (slot, (_, d)) in &base {
            let (_, e) = &permuted[slot];
            let expect: Vec<Parent> = d
                .parents
                .iter()
                .map(|p| match p {
                    Parent::Label(l) => Parent::Label(InputLabel(map[&l.0])),
                    other => *other,
                })
                .collect();
            assert_eq!(e.parents, expect);
            assert_eq!((e.role, e.ordinal, e.entity), (d.role, d.ordinal, d.entity));
        }
        let unlabelled = ids(&build(&s, None, 1.0, false).unwrap());
        assert!(id_set(&unlabelled).is_disjoint(&id_set(&base)));
    }
}

fn body_ids(body: &Body) -> Ids {
    let t = body.topology();
    let mut out = BTreeMap::new();
    for (id, slot) in t.ids() {
        let d = t.derivation(id).unwrap().clone();
        assert_eq!(fnv(&encode(&d)), id.0, "{d:?}");
        assert_eq!(t.slot_of(id), Some(slot));
        assert_eq!(d.operation, body.operation());
        out.insert(slot, (id.to_string(), d));
    }
    let total = t.vertices().len() + t.edges().len() + t.faces().len();
    assert_eq!(out.len(), total, "every slot has an id");
    out
}

/// The parents of every derivation of a role, sorted.
fn parents_of(ids: &Ids, role: Role) -> Vec<Vec<Parent>> {
    let mut out: Vec<Vec<Parent>> = ids
        .values()
        .filter(|v| v.1.role == role)
        .map(|v| v.1.parents.clone())
        .collect();
    out.sort();
    out
}

/// S6: the prism's profile as a face body and its outer boundary as a wire.
fn check_bodies(solid: &Solid, prism: &Ids, transforms: &[RigidTransform]) {
    let profile = solid.profile().unwrap().clone();
    let tolerance = profile.tolerance();
    let op = solid.operation();
    // The prism validated, so its start cap's face does.
    let (face, history) = Body::face_from_profile_with(op, profile.clone(), solid.frame())
        .expect("the prism's profile makes a face");
    let t = face.topology();
    assert_eq!(t.check(tolerance), Vec::new());
    assert_eq!(face.class().name(), "sheet");
    assert_eq!(history.relations.len(), t.ids().count());
    let ids = body_ids(&face);
    assert!(ids.values().all(|v| v.1.kind == OperationKind::MakeFace));
    // The same elements as the prism's start cap, bottom (start-side) edges
    // and vertices, whichever way it was extruded.
    assert_eq!(
        parents_of(&ids, Role::Face),
        parents_of(prism, Role::StartCap)
    );
    assert_eq!(
        parents_of(&ids, Role::Edge),
        parents_of(prism, Role::BottomEdge)
    );
    assert_eq!(
        parents_of(&ids, Role::Vertex),
        parents_of(prism, Role::BottomVertex)
    );
    let boundaries = 1 + profile.holes().len();
    let rings = t.edges().iter().filter(|e| e.is_ring()).count();
    let c = t.occt_counts();
    assert_eq!(
        (c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids),
        (
            t.vertices().len() + rings,
            t.edges().len(),
            boundaries,
            1,
            0,
            0
        )
    );
    let area = profile.area();
    let m = face.measure().expect("a planar face integrates");
    let slack = 1e-12 * area.abs().max(1e-300);
    assert!(m.measure[0] - slack <= area && area <= m.measure[1] + slack);
    let mut moved = face.clone();
    for transform in transforms {
        if let Ok((next, _)) = moved.transform_with(OperationId::UNSPECIFIED, *transform) {
            assert_eq!(body_ids(&next), ids);
            moved = next;
        }
    }
    let (wire, _) =
        Body::wire_from_boundary_with(op, profile.outer().clone(), solid.frame(), tolerance)
            .expect("the outer boundary makes a wire");
    let w = wire.topology();
    assert_eq!(w.check(tolerance), Vec::new());
    assert_eq!(wire.class().name(), "wire");
    let wire_ids = body_ids(&wire);
    assert!(wire_ids
        .values()
        .all(|v| v.1.kind == OperationKind::MakeWire && v.1.role != Role::Face));
    let c = w.occt_counts();
    let n = w.edges().len();
    assert_eq!(
        (c.edges, c.wires, c.faces, c.shells),
        (n, usize::from(n > 1), 0, 0)
    );
    assert!(wire.measure().is_some());
}
