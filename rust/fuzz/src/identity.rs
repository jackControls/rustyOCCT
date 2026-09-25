//! Value ids on structure-aware extrusions (M1 of IDENTITY_AND_HISTORY.md).
//! An independent byte encoder and FNV-1a-128 recompute every id from its
//! retained derivation. Rebuilding, rigid motion, reversing the direction and
//! moving labelled points keep every id; permuting label values permutes ids
//! bijectively, parent by parent; counts and roles follow the profile.
use libfuzzer_sys::arbitrary::{Result, Unstructured};
use rusty_occt::history::History;
use rusty_occt::identity::{
    Derivation, EntityKind, InputLabel, OperationId, OperationKind, Parent, ProfileElement, Role,
};
use rusty_occt::topology::Slot;
use rusty_occt::{
    Boundary, BoundaryLabels, Frame3, Point2, Point3, Profile, RigidTransform, Solid, Tolerance,
    Vec3,
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
    });
    out.push(match d.entity {
        EntityKind::Vertex => 1,
        EntityKind::Edge => 2,
        EntityKind::Face => 3,
        EntityKind::Body => 4,
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
    let total = t.vertices().len() + t.edges().len() + t.faces().len();
    assert_eq!(out.len(), total, "every slot has an id");
    let distinct: BTreeSet<&String> = out.values().map(|v| &v.0).collect();
    assert_eq!(distinct.len(), total, "ids are unique");
    assert!(!distinct.contains(&t.body_id().to_string()));
    out
}

fn id_set(ids: &Ids) -> BTreeSet<String> {
    ids.values().map(|v| v.0.clone()).collect()
}

pub fn check_identity(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(s)) = spec(&mut u) else {
        return;
    };
    let pool = s.labels.as_deref();
    let Some(solid) = build(&s, pool, 1.0, false) else {
        return;
    };
    let base = ids(&solid);
    // Counts and roles follow the profile structure.
    let boundaries = 1 + s.holes.len();
    let polygon_sides: usize =
        s.outer.as_ref().map_or(0, Vec::len) + 4 * s.holes.iter().filter(|h| h.0).count();
    let circles = usize::from(s.outer.is_none()) + s.holes.iter().filter(|h| !h.0).count();
    let t = solid.topology();
    assert_eq!(t.vertices().len(), 2 * polygon_sides + 2 * circles);
    assert_eq!(t.edges().len(), 3 * polygon_sides + 3 * circles);
    assert_eq!(t.faces().len(), 2 + polygon_sides + circles);
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
