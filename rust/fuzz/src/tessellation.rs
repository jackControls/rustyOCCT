//! Tessellation (T-a and T-b of REVIEW_NOTES.md) of structure-aware bodies:
//! prisms of polygons with square and circular holes, filleted and notched
//! arc profiles, their face bodies, cones, spheres and zones, tori,
//! v-segments and wedges, each possibly moved rigidly, and (T-b) prisms
//! with spline sides and face bodies on spline surfaces (`spline.rs`), at a
//! deflection from an eighth to a 2048th of the body's size (a 64th for
//! splines) and an angle from 0.2 rad (0.35 for splines) to π/2. Every
//! result must be a mesh (no error on a buildable body): distinct nodes per
//! triangle, every mesh edge in two triangles traversed once each way (a
//! face body's boundary edges once), the Euler characteristic of the body's
//! boundary, every reported bound within the request, 12 barycentric
//! samples of every triangle within its reported bound of its face's whole
//! surface (closed-form distances, independent of the kernel's evaluation;
//! on a spline surface the centroids of up to 64 triangles per face, by
//! projection), a solid's enclosed volume within the deflection times the
//! areas of its mass properties (a spline prism's by quadrature of its
//! profile), and the same mesh again, bit for bit (a spline body's up to
//! 500 triangles).
use crate::identity::{arc_path, build, cone_spec, spec, sphere_spec, torus_spec};
use libfuzzer_sys::arbitrary::{Result, Unstructured};
use rusty_occt::identity::OperationId;
use rusty_occt::tessellation::{tessellate, Mesh, Parameters};
use rusty_occt::topology::{Surface, Topology};
use rusty_occt::{Body, Boundary, Frame3, Point3, Profile, Solid, Tolerance, Vec3};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::FRAC_PI_2;

mod spline;

enum Built {
    Solid(Solid),
    Face(Body, usize),
    /// T-b: a spline prism (a solid, with its volume and area by
    /// quadrature) or a spline face body.
    Spline(Topology, Option<(f64, f64)>),
}

fn unit(u: &mut Unstructured) -> Result<f64> {
    Ok(f64::from(u.arbitrary::<u16>()?) / 65535.0)
}

fn moved(solid: Solid, transforms: &[rusty_occt::RigidTransform]) -> Solid {
    match transforms.first() {
        Some(t) => match solid.transform_with(OperationId(2), *t) {
            Ok((moved, _)) => moved,
            Err(_) => solid,
        },
        None => solid,
    }
}

fn body(u: &mut Unstructured) -> Result<Option<Built>> {
    Ok(match u.int_in_range(0u8..=8)? {
        7 => spline::prism(u)?.map(|(t, v, a)| Built::Spline(t, Some((v, a)))),
        8 => spline::sheet(u)?.map(|t| Built::Spline(t, None)),
        0 => {
            let Some(s) = spec(u)? else { return Ok(None) };
            build(&s, None, 1.0, false).map(|solid| Built::Solid(moved(solid, &s.transforms)))
        }
        kind @ (1 | 2) => {
            let (points, segments, _) = arc_path(u)?;
            let scale = 2f64.powi(u.int_in_range(-6..=6)?);
            let tolerance = Tolerance::new(1e-9 * scale, 1e-12).unwrap();
            let points = points
                .iter()
                .map(|p| rusty_occt::Point2::new(p.x * scale, p.y * scale))
                .collect();
            let segments = segments
                .into_iter()
                .map(|g| match g {
                    rusty_occt::Segment::Arc {
                        center,
                        radius,
                        ccw,
                    } => rusty_occt::Segment::Arc {
                        center: rusty_occt::Point2::new(center.x * scale, center.y * scale),
                        radius: radius * scale,
                        ccw,
                    },
                    line => line,
                })
                .collect();
            let Ok(outer) = Boundary::path(points, segments, tolerance) else {
                return Ok(None);
            };
            let Ok(profile) = Profile::new(outer, Vec::new(), tolerance) else {
                return Ok(None);
            };
            let normal = Vec3::new(2.0 * unit(u)? - 1.0, 2.0 * unit(u)? - 1.0, 0.3 + unit(u)?);
            let Ok(frame) =
                Frame3::new(Point3::new(scale, -scale, 0.0), normal, Vec3::X, tolerance)
            else {
                return Ok(None);
            };
            if kind == 1 {
                let height = scale * (0.1 + 2.0 * unit(u)?);
                Solid::extrude_with(OperationId(1), profile, frame, 0.0, height)
                    .ok()
                    .map(|(s, _)| Built::Solid(s))
            } else {
                Body::face_from_profile_with(OperationId(1), profile, frame)
                    .ok()
                    .map(|(b, _)| Built::Face(b, 0))
            }
        }
        3 => {
            let Some(s) = cone_spec(u)? else {
                return Ok(None);
            };
            s.build(1.0)
                .map(|(solid, _)| Built::Solid(moved(solid, &s.transforms)))
        }
        4 => {
            let Some(s) = sphere_spec(u)? else {
                return Ok(None);
            };
            s.build(1.0)
                .map(|(solid, _)| Built::Solid(moved(solid, &s.transforms)))
        }
        5 => {
            let Some(s) = torus_spec(u)? else {
                return Ok(None);
            };
            s.build(1.0)
                .map(|(solid, _)| Built::Solid(moved(solid, &s.transforms)))
        }
        6 => {
            // A prism's profile as a face body.
            let Some(s) = spec(u)? else { return Ok(None) };
            let Some(solid) = build(&s, None, 1.0, false) else {
                return Ok(None);
            };
            let profile = solid.profile().unwrap().clone();
            let holes = profile.holes().len();
            Body::face_from_profile_with(OperationId(1), profile, solid.frame())
                .ok()
                .map(|(b, _)| Built::Face(b, holes))
        }
        _ => unreachable!("kinds 0 to 8"),
    })
}

/// The distance from a point to a face's whole surface, in closed form.
fn surface_distance(surface: &Surface, p: Point3) -> f64 {
    let local = |frame: &Frame3| {
        let w = p - frame.origin();
        let z = w.dot(frame.normal());
        ((w.dot(w) - z * z).max(0.0).sqrt(), z)
    };
    match surface {
        Surface::Plane(frame) => (p - frame.origin()).dot(frame.normal()).abs(),
        Surface::Cylinder { frame, radius } => (local(frame).0 - radius).abs(),
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => {
            let (rho, z) = local(frame);
            let (s, c) = half_angle.sin_cos();
            [rho, -rho]
                .iter()
                .map(|r| ((r - radius) * c - z * s).abs())
                .fold(f64::INFINITY, f64::min)
        }
        Surface::Sphere { frame, radius } => ((p - frame.origin()).length() - radius).abs(),
        Surface::Torus {
            frame,
            major,
            minor,
        } => {
            let (rho, z) = local(frame);
            ((rho - major).hypot(z) - minor).abs()
        }
        Surface::BSpline(_) => f64::INFINITY,
    }
}

/// The contract; returns the Euler characteristic, the number of boundary
/// edges, the area and the enclosed volume.
fn contract(t: &Topology, mesh: &Mesh, p: Parameters) -> (i64, usize, f64, f64) {
    assert!(mesh.deflection <= p.deflection() && mesh.angle <= p.angle());
    let mut directed: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for tri in &mesh.triangles {
        assert!(tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]);
        for k in 0..3 {
            *directed.entry((tri[k], tri[(k + 1) % 3])).or_default() += 1;
        }
    }
    let mut boundary = 0;
    let mut undirected = BTreeSet::new();
    for (&(a, b), &k) in &directed {
        assert_eq!(k, 1, "a mesh edge used twice in one direction");
        boundary += usize::from(!directed.contains_key(&(b, a)));
        undirected.insert((a.min(b), a.max(b)));
    }
    let used: BTreeSet<usize> = mesh.triangles.iter().flatten().copied().collect();
    let extent = mesh
        .nodes
        .iter()
        .flat_map(|p| p.to_array())
        .fold(1.0f64, |e, x| e.max(x.abs()));
    let eps = 1e-12 * extent;
    let (mut area, mut volume) = (0.0, 0.0);
    let origin = Point3::new(0.0, 0.0, 0.0);
    for f in &mesh.faces {
        let surface = &t.faces()[f.face.index()].surface;
        let projection = match surface {
            Surface::BSpline(s) => Some(spline::Projection::new(s)),
            _ => None,
        };
        let stride = f.triangles.len().div_ceil(64).max(1);
        for k in f.triangles.clone() {
            let x = mesh.triangles[k].map(|n| mesh.nodes[n]);
            let bound = mesh.triangle_bounds[k];
            assert!(bound.deflection <= p.deflection() && bound.angle <= p.angle());
            if let Some(projection) = &projection {
                if (k - f.triangles.start) % stride == 0 {
                    let centroid = Point3::new(
                        (x[0].x + x[1].x + x[2].x) / 3.0,
                        (x[0].y + x[1].y + x[2].y) / 3.0,
                        (x[0].z + x[1].z + x[2].z) / 3.0,
                    );
                    let d = projection.distance(centroid);
                    assert!(d <= bound.deflection + eps, "centroid {d} beyond {bound:?}");
                }
            }
            for i in 0..=4 {
                for j in 0..=4 - i {
                    let l = [i as f64 / 4.0, j as f64 / 4.0, (4 - i - j) as f64 / 4.0];
                    if l.contains(&1.0) {
                        continue;
                    }
                    let s = Point3::new(
                        l[0] * x[0].x + l[1] * x[1].x + l[2] * x[2].x,
                        l[0] * x[0].y + l[1] * x[1].y + l[2] * x[2].y,
                        l[0] * x[0].z + l[1] * x[1].z + l[2] * x[2].z,
                    );
                    if projection.is_some() {
                        continue;
                    }
                    let d = surface_distance(surface, s);
                    assert!(d <= bound.deflection + eps, "sample {d} beyond {bound:?}");
                }
            }
            area += (x[1] - x[0]).cross(x[2] - x[0]).length() / 2.0;
            volume += (x[0] - origin).dot((x[1] - origin).cross(x[2] - origin)) / 6.0;
        }
    }
    let euler = used.len() as i64 - undirected.len() as i64 + mesh.triangles.len() as i64;
    (euler, boundary, area, volume)
}

pub fn check_tessellation(data: &[u8]) {
    let mut u = Unstructured::new(data);
    let Ok(Some(built)) = body(&mut u) else {
        return;
    };
    // Splines: a deflection down to a 64th of the size and an angle from
    // 0.35 rad, which keep an input within seconds under the sanitizer.
    let (finest, sharpest) = if matches!(built, Built::Spline(..)) {
        (6, 1)
    } else {
        (11, 0)
    };
    let (Ok(k), Ok(a)) = (u.int_in_range(3u8..=finest), u.int_in_range(sharpest..=5u8)) else {
        return;
    };
    let angle = [0.2, 0.35, 0.5, 0.8, 1.2, FRAC_PI_2][a as usize];
    let topology = match &built {
        Built::Solid(s) => s.topology(),
        Built::Face(b, _) => b.topology(),
        Built::Spline(t, _) => t,
    };
    let points: Vec<Point3> = topology
        .vertices()
        .iter()
        .map(|v| v.position)
        .chain(
            topology
                .edges()
                .iter()
                .flat_map(|e| [0.0, 0.25, 0.5, 0.75].map(|f| e.curve.point(f))),
        )
        .collect();
    let mut size = 0.0f64;
    for p in &points {
        for q in &points {
            size = size.max(p.distance(*q));
        }
    }
    if let Built::Solid(s) = &built {
        size = size.max(s.mass_properties().volume.abs().cbrt());
    }
    let p = Parameters::new(size * 2f64.powi(-i32::from(k)), angle).unwrap();
    let mesh = tessellate(topology, p).unwrap_or_else(|e| panic!("{e}"));
    let (euler, boundary, area, volume) = contract(topology, &mesh, p);
    match &built {
        Built::Solid(s) => {
            let whole_torus = topology.faces().len() == 1
                && matches!(topology.faces()[0].surface, Surface::Torus { .. });
            let holes = s.profile().map_or(0, |pr| pr.holes().len()) as i64;
            assert_eq!(euler, if whole_torus { 0 } else { 2 - 2 * holes });
            assert_eq!(boundary, 0, "not closed");
            let m = s.mass_properties();
            let slack =
                p.deflection() * (m.surface_area + area) * (1.0 + 1e-9) + 1e-9 * m.volume.abs();
            assert!(
                (volume - m.volume).abs() <= slack,
                "volume {volume} against {}",
                m.volume
            );
        }
        Built::Face(_, holes) => {
            assert_eq!(euler, 1 - *holes as i64);
            assert!(boundary >= 3 * (1 + holes));
        }
        Built::Spline(_, Some((exact, surface))) => {
            assert_eq!(euler, 2);
            assert_eq!(boundary, 0, "not closed");
            // The certified mass enclosure of a spline wall costs seconds
            // under the sanitizer: the profile's quadrature instead.
            let slack = p.deflection() * (surface + area) * (1.0 + 1e-6) + 1e-9 * exact;
            assert!(
                (volume - exact).abs() <= slack,
                "volume {volume} against {exact}"
            );
        }
        Built::Spline(_, None) => {
            assert_eq!(euler, 1);
            assert!(boundary >= 3);
        }
    }
    // A spline body's second tessellation is a sanitizer's seconds: the
    // smaller meshes only.
    if !matches!(built, Built::Spline(..)) || mesh.triangles.len() <= 500 {
        assert_eq!(tessellate(topology, p).unwrap(), mesh, "deterministic");
    }
}
