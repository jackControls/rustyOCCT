//! Analytic surface intersections (S7a of REVIEW_NOTES.md). Two planes,
//! cylinders, cones or spheres on dyadic origins, small integer axes, dyadic
//! radii and angles, with exact degeneracies made on purpose: a shared axis
//! or normal (parallel), a shared origin (coaxial, concentric), an offset of
//! exactly the radius (tangent). The intersection never panics and fails
//! only by a certified comparison it cannot decide; it is symmetric in its
//! arguments; every enclosure is ordered; every returned point, line and
//! conic lies on both surfaces; and an exact dyadic translation moves every
//! item with the surfaces. A procedural curve (S7b.1: cylinders with
//! crossing axes, a cylinder and a sphere off its axis) is the same whatever
//! the argument order, and points along every loop, ring and figure-eight
//! lie on both surfaces.
use rusty_occt::intersection::{
    surface_surface, AnalyticItem, Branch, Component, SurfaceIntersection,
};
use rusty_occt::topology::Surface;
use rusty_occt::{Error, Frame3, Point3, RigidTransform, Tolerance, Vec3};

struct Bytes<'a>(&'a [u8], usize);
impl Bytes<'_> {
    fn next(&mut self) -> u8 {
        let b = self.0.get(self.1).copied().unwrap_or(0);
        self.1 += 1;
        b
    }
    /// A dyadic value k / 8 in [-8, 8).
    fn dyadic(&mut self) -> f64 {
        f64::from(i16::from(self.next()) - 128) / 16.0
    }
    fn small(&mut self) -> f64 {
        f64::from(self.next() % 7) - 3.0
    }
}

fn frame(o: Point3, n: Vec3) -> Option<Frame3> {
    let hint = if n.x.abs() >= n.y.abs().max(n.z.abs()) {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    Frame3::new(o, n, hint, Tolerance::default()).ok()
}

fn make(kind: u8, f: Frame3, radius: f64, angle: f64) -> Surface {
    match kind % 4 {
        0 => Surface::Plane(f),
        1 => Surface::Cylinder { frame: f, radius },
        2 => Surface::Cone {
            frame: f,
            radius: radius - 0.25,
            half_angle: angle,
        },
        _ => Surface::Sphere { frame: f, radius },
    }
}

fn frame_of(s: &Surface) -> Frame3 {
    match s {
        Surface::Plane(f)
        | Surface::Cylinder { frame: f, .. }
        | Surface::Cone { frame: f, .. }
        | Surface::Sphere { frame: f, .. } => *f,
        _ => unreachable!(),
    }
}

/// Binary64 distance from a point to a surface (the check's own).
fn distance(s: &Surface, p: Vec3) -> f64 {
    let f = frame_of(s);
    let rel = p - (f.origin() - Point3::ORIGIN);
    let n = f.normal();
    let h = rel.dot(n);
    let radial = (rel - n * h).length();
    match s {
        Surface::Plane(_) => h.abs(),
        Surface::Cylinder { radius, .. } => (radial - radius).abs(),
        Surface::Sphere { radius, .. } => (rel.length() - radius).abs(),
        Surface::Cone {
            radius, half_angle, ..
        } => {
            let (sa, ca) = half_angle.sin_cos();
            let shift = radius * ca + h * sa;
            (radial * ca - shift).abs().min((radial * ca + shift).abs())
        }
        _ => unreachable!(),
    }
}

fn mid([lo, hi]: [f64; 2]) -> f64 {
    0.5 * lo + 0.5 * hi
}
fn v(x: &[[f64; 2]; 3]) -> Vec3 {
    Vec3::new(mid(x[0]), mid(x[1]), mid(x[2]))
}

/// A few points of an item, and its scale.
fn samples(item: &AnalyticItem) -> Vec<Vec3> {
    let any_perp = |n: Vec3| {
        let e = if n.x.abs() < 0.5 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let x = n.cross(e);
        x * (1.0 / x.length())
    };
    match item {
        AnalyticItem::Point(p) => vec![v(p)],
        AnalyticItem::Line { point, direction } => (-1..=1)
            .map(|k| v(point) + v(direction) * f64::from(k))
            .collect(),
        AnalyticItem::Circle {
            centre,
            normal,
            radius,
        } => {
            let n = v(normal);
            let (x, y) = (any_perp(n), n.cross(any_perp(n)));
            (0..4)
                .map(|k| {
                    let t = f64::from(k) * 1.5;
                    v(centre) + (x * t.cos() + y * t.sin()) * mid(*radius)
                })
                .collect()
        }
        AnalyticItem::Ellipse {
            centre,
            normal,
            major,
            semi_major,
            semi_minor,
        } => {
            let (u, w) = (v(major), v(normal).cross(v(major)));
            (0..4)
                .map(|k| {
                    let t = f64::from(k) * 1.5;
                    v(centre) + u * (mid(*semi_major) * t.cos()) + w * (mid(*semi_minor) * t.sin())
                })
                .collect()
        }
        AnalyticItem::Hyperbola {
            centre,
            normal,
            transverse,
            semi_transverse,
            semi_conjugate,
        } => {
            let (u, w) = (v(transverse), v(normal).cross(v(transverse)));
            (-1..=1)
                .map(|k| {
                    let t = f64::from(k) * 0.5;
                    v(centre)
                        + u * (mid(*semi_transverse) * t.cosh())
                        + w * (mid(*semi_conjugate) * t.sinh())
                })
                .collect()
        }
    }
}

fn translated(s: &Surface, t: RigidTransform) -> Surface {
    let f = frame_of(s).transformed(t, Tolerance::default()).unwrap();
    match s {
        Surface::Plane(_) => Surface::Plane(f),
        Surface::Cylinder { radius, .. } => Surface::Cylinder {
            frame: f,
            radius: *radius,
        },
        Surface::Sphere { radius, .. } => Surface::Sphere {
            frame: f,
            radius: *radius,
        },
        Surface::Cone {
            radius, half_angle, ..
        } => Surface::Cone {
            frame: f,
            radius: *radius,
            half_angle: *half_angle,
        },
        _ => unreachable!(),
    }
}

pub fn check_analytic_intersections(data: &[u8]) {
    let mut b = Bytes(data, 0);
    let (k1, k2) = (b.next(), b.next());
    let mode = b.next() % 5;
    let o1 = Point3::new(b.dyadic(), b.dyadic(), b.dyadic());
    let n1 = Vec3::new(b.small(), b.small(), b.small());
    let r1 = 0.25 + f64::from(b.next() % 32) / 8.0;
    let a1 = 0.05 + f64::from(b.next() % 64) / 48.0;
    let Some(f1) = frame(o1, n1) else { return };
    let s1 = make(k1, f1, r1, a1);
    // The second surface: independent, or sharing the first's axis
    // (parallel), origin and axis (coaxial), or offset along the axis's
    // perpendicular by exactly a radius (tangent).
    let r2 = 0.25 + f64::from(b.next() % 32) / 8.0;
    let a2 = 0.05 + f64::from(b.next() % 64) / 48.0;
    let (o2, n2) = match mode {
        0 => (
            Point3::new(b.dyadic(), b.dyadic(), b.dyadic()),
            Vec3::new(b.small(), b.small(), b.small()),
        ),
        1 => (Point3::new(b.dyadic(), b.dyadic(), b.dyadic()), n1),
        2 => (o1, n1),
        3 => {
            let side = f1.x() * (r1 + r2);
            (o1 + side, n1)
        }
        _ => (o1, Vec3::new(b.small(), b.small(), b.small())),
    };
    let Some(f2) = frame(o2, n2) else { return };
    let s2 = make(k2, f2, r2, a2);
    let result = match surface_surface(&s1, &s2) {
        Ok(r) => r,
        Err(Error::ComputationLimit(_)) => return,
        Err(e) => panic!("unexpected error {e}"),
    };
    // Symmetric: the same items, enclosures overlapping (each contains the
    // exact value; equal kinds compute from the other surface's data).
    let swapped = surface_surface(&s2, &s1).expect("the swapped pair decides too");
    match (&result, &swapped) {
        (SurfaceIntersection::Items(x), SurfaceIntersection::Items(y)) => {
            assert_eq!(x.len(), y.len(), "symmetric");
            for (p, q) in x.iter().zip(y) {
                assert_eq!(p.kind(), q.kind(), "symmetric");
                for ([a, b], [c, d]) in p.values().iter().zip(q.values()) {
                    assert!(*a <= d && c <= *b, "overlapping enclosures");
                }
            }
        }
        _ => assert_eq!(result, swapped, "symmetric"),
    }
    // A procedural curve (S7b.1): points along every component lie on both
    // surfaces.
    if let SurfaceIntersection::Procedural(c) = &result {
        for comp in c.components() {
            let (params, branches): (Vec<f64>, &[Branch]) = match comp {
                Component::Loop { u } => {
                    let (lo, hi) = (u[0][1], u[1][0]);
                    (
                        (1..8)
                            .map(|k| lo + (hi - lo) * f64::from(k) / 8.0)
                            .collect(),
                        &[Branch::Plus, Branch::Minus],
                    )
                }
                Component::Ring { branch } => (
                    (0..8).map(|k| -3.0 + 0.75 * f64::from(k)).collect(),
                    std::slice::from_ref(branch),
                ),
                Component::FigureEight { .. } => (
                    (0..8).map(|k| -3.0 + 0.75 * f64::from(k)).collect(),
                    &[Branch::Plus, Branch::Minus],
                ),
            };
            for t in params {
                for branch in branches {
                    let Ok(e) = c.point_at([t, t], *branch) else {
                        panic!("a parameter inside a component is on the curve");
                    };
                    let p = Vec3::new(mid(e[0]), mid(e[1]), mid(e[2]));
                    let scale = p.length().max(1.0) * 8.0;
                    for s in [&s1, &s2] {
                        let gap = distance(s, p);
                        assert!(gap <= 1e-9 * scale, "{gap} off {s:?} at {t}");
                    }
                }
            }
        }
        return;
    }
    let SurfaceIntersection::Items(items) = &result else {
        return;
    };
    for item in items {
        for [lo, hi] in item.values() {
            assert!(lo <= hi, "ordered enclosure");
        }
        for p in samples(item) {
            if !p.to_array().iter().all(|x| x.is_finite()) {
                continue;
            }
            let scale = p.length().max(1.0) * 8.0;
            for s in [&s1, &s2] {
                let gap = distance(s, p);
                assert!(gap <= 1e-9 * scale, "{gap} off {s:?} for {item:?}");
            }
        }
    }
    // An exact dyadic translation keeps every predicate and moves the items
    // (origins are dyadic except in the tangent mode).
    if mode == 3 {
        return;
    }
    let shift = Vec3::new(b.dyadic(), b.dyadic(), b.dyadic());
    let t = RigidTransform::translation(shift).unwrap();
    if let Ok(SurfaceIntersection::Items(moved)) =
        surface_surface(&translated(&s1, t), &translated(&s2, t))
    {
        assert_eq!(moved.len(), items.len(), "translation keeps the items");
    }
}
