//! Tessellation (T-a and T-b of REVIEW_NOTES.md) against the independent
//! reference's expectations (fixtures/tessellation-*.{txt,tsv} from
//! tools/generate_tessellation_fixtures.py), and the contract on the valid
//! spline cases of the B-rep fixtures and on imported corpus solids.
#[path = "support/brep_protocol.rs"]
mod brep_protocol;
#[path = "support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use rusty_occt::curve::{DerivativeOrder, KnotSide};
use rusty_occt::occt_brep::{import, read};
use rusty_occt::tessellation::{tessellate, Mesh, Parameters};
use rusty_occt::topology::{Surface, Topology};
use rusty_occt::{BSplineSurface3, Error, Point3, Tolerance, Vec3};
use std::collections::BTreeMap;

/// Barycentric samples: the order-4 lattice without the vertices.
fn samples() -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for i in 0..=4 {
        for j in 0..=4 - i {
            let k = 4 - i - j;
            if i < 4 && j < 4 && k < 4 {
                out.push([i as f64 / 4.0, j as f64 / 4.0, k as f64 / 4.0]);
            }
        }
    }
    out
}

/// The distance from a point to a face's whole surface, in closed form
/// (independent of the kernel's evaluation): a plane, the revolved line of
/// a cylinder or cone (either nappe), a sphere, a torus.
fn surface_distance(surface: &Surface, p: Point3) -> f64 {
    let local = |frame: &rusty_occt::Frame3| {
        let w = p - frame.origin();
        let z = w.dot(frame.normal());
        let rho = (w.dot(w) - z * z).max(0.0).sqrt();
        (rho, z)
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
            // The generatrix through (radius, 0) along (sin a, cos a), and
            // its mirror image.
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

/// Distances to a spline surface over its domain by projection: the
/// nearest points of a grid with four samples per knot span in each
/// direction (at most 257), each followed by Gauss-Newton iterations clamped
/// to the domain, with the kernel's exact evaluation (independent of the
/// tessellation's cells). A projection can only overstate a distance.
struct SplineDistance<'a> {
    surface: &'a BSplineSurface3,
    grid: Vec<(Point3, f64, f64)>,
}

impl<'a> SplineDistance<'a> {
    fn new(surface: &'a BSplineSurface3) -> Self {
        let ((u0, u1), (v0, v1)) = surface.domain();
        let count = |k: &rusty_occt::KnotVector| (4 * (k.knots().len() - 1)).clamp(16, 256);
        let (nu, nv) = (count(surface.u_knots()), count(surface.v_knots()));
        let mut grid = Vec::new();
        for i in 0..=nu {
            for j in 0..=nv {
                let (u, v) = (
                    (u0 + (u1 - u0) * i as f64 / nu as f64).min(u1),
                    (v0 + (v1 - v0) * j as f64 / nv as f64).min(v1),
                );
                grid.push((surface.point(u, v).unwrap(), u, v));
            }
        }
        Self { surface, grid }
    }

    fn distance(&self, p: Point3) -> f64 {
        let mut order: Vec<(f64, f64, f64)> = self
            .grid
            .iter()
            .map(|(q, u, v)| (q.distance(p), *u, *v))
            .collect();
        order.sort_by(|a, b| a.0.total_cmp(&b.0));
        order
            .iter()
            .take(4)
            .map(|&(d, u, v)| d.min(self.descend(p, u, v)))
            .fold(f64::INFINITY, f64::min)
    }

    fn descend(&self, p: Point3, mut u: f64, mut v: f64) -> f64 {
        let ((u0, u1), (v0, v1)) = self.surface.domain();
        let vector = |b: [rusty_occt::ScalarInterval; 3]| {
            Vec3::new(
                b[0].representative(),
                b[1].representative(),
                b[2].representative(),
            )
        };
        let mut best = f64::INFINITY;
        for _ in 0..20 {
            let e = self
                .surface
                .evaluate(u, v, DerivativeOrder::First, [KnotSide::Automatic; 2])
                .unwrap();
            let r = e.position() - p;
            best = best.min(r.length());
            let (su, sv) = (
                vector(e.derivative_bounds(1, 0).unwrap()),
                vector(e.derivative_bounds(0, 1).unwrap()),
            );
            let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
            let (g, h) = (-r.dot(su), -r.dot(sv));
            let det = a * c - b * b;
            if det.is_nan() || det <= 0.0 {
                break;
            }
            let (du, dv) = ((g * c - h * b) / det, (a * h - b * g) / det);
            let (nu, nv) = ((u + du).clamp(u0, u1), (v + dv).clamp(v0, v1));
            if (nu - u).abs() <= 1e-15 * (u1 - u0) && (nv - v).abs() <= 1e-15 * (v1 - v0) {
                break;
            }
            (u, v) = (nu, nv);
        }
        best.min(self.surface.point(u, v).unwrap().distance(p))
    }
}

fn kind(surface: &Surface) -> &'static str {
    match surface {
        Surface::Plane(_) => "plane",
        Surface::Cylinder { .. } => "cylinder",
        Surface::Cone { .. } => "cone",
        Surface::Sphere { .. } => "sphere",
        Surface::Torus { .. } => "torus",
        Surface::BSpline(_) => "spline",
    }
}

fn extent(t: &Topology, mesh: &Mesh) -> f64 {
    let mut e = 1.0f64;
    for p in t
        .vertices()
        .iter()
        .map(|v| v.position)
        .chain(mesh.nodes.iter().copied())
    {
        e = e.max(p.x.abs()).max(p.y.abs()).max(p.z.abs());
    }
    e
}

struct Measured {
    euler: i64,
    boundary_edges: usize,
    area: f64,
    volume: f64,
}

/// The contract, measured without the kernel's own check: distinct nodes
/// per triangle, every mesh edge in two triangles in opposite directions
/// (a face body's boundary edges once), each sample within its triangle's
/// bound of the face's whole surface, every bound within the request.
fn contract(name: &str, t: &Topology, mesh: &Mesh, p: Parameters, solid: bool) -> Measured {
    assert!(
        mesh.deflection <= p.deflection(),
        "{name}: {}",
        mesh.deflection
    );
    assert!(mesh.angle <= p.angle(), "{name}: {}", mesh.angle);
    let mut directed: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for tri in &mesh.triangles {
        assert!(
            tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2],
            "{name}"
        );
        for k in 0..3 {
            *directed.entry((tri[k], tri[(k + 1) % 3])).or_default() += 1;
        }
    }
    let mut boundary = 0;
    let mut undirected = BTreeMap::new();
    for (&(a, b), &k) in &directed {
        assert_eq!(k, 1, "{name}: a mesh edge used twice in one direction");
        if !directed.contains_key(&(b, a)) {
            boundary += 1;
        }
        undirected.insert((a.min(b), a.max(b)), ());
    }
    if solid {
        assert_eq!(boundary, 0, "{name}: not closed");
    }
    let used: std::collections::BTreeSet<usize> =
        mesh.triangles.iter().flatten().copied().collect();
    let eps = 1e-12 * extent(t, mesh);
    let (mut area, mut volume) = (0.0, 0.0);
    for f in &mesh.faces {
        let surface = &t.faces()[f.face.index()].surface;
        // A spline face: the centroids of at most 100 of its triangles,
        // evenly spread, by projection.
        let spline = match surface {
            Surface::BSpline(s) => Some(SplineDistance::new(s)),
            _ => None,
        };
        let stride = if spline.is_some() {
            f.triangles.len().div_ceil(100).max(1)
        } else {
            1
        };
        for k in f.triangles.clone() {
            let x = mesh.triangles[k].map(|n| mesh.nodes[n]);
            let bound = mesh.triangle_bounds[k];
            assert!(
                bound.deflection <= p.deflection() && bound.angle <= p.angle(),
                "{name}"
            );
            let lattice = if spline.is_none() {
                samples()
            } else if (k - f.triangles.start) % stride == 0 {
                vec![[1.0 / 3.0; 3]]
            } else {
                Vec::new()
            };
            for l in lattice {
                let s = Point3::new(
                    l[0] * x[0].x + l[1] * x[1].x + l[2] * x[2].x,
                    l[0] * x[0].y + l[1] * x[1].y + l[2] * x[2].y,
                    l[0] * x[0].z + l[1] * x[1].z + l[2] * x[2].z,
                );
                let d = match &spline {
                    Some(spline) => spline.distance(s),
                    None => surface_distance(surface, s),
                };
                assert!(
                    d <= bound.deflection + eps,
                    "{name}: face {} ({}) sample {d} beyond bound {bound:?}",
                    f.face.index(),
                    kind(surface)
                );
            }
            let n: Vec3 = (x[1] - x[0]).cross(x[2] - x[0]);
            area += n.length() / 2.0;
            let o = Point3::new(0.0, 0.0, 0.0);
            volume += (x[0] - o).dot((x[1] - o).cross(x[2] - o)) / 6.0;
        }
    }
    for e in &mesh.edges {
        assert_eq!(e.segments.len() + 1, e.nodes.len(), "{name}");
        assert!(e
            .segments
            .iter()
            .all(|s| s.deflection <= p.deflection() && s.angle <= p.angle()));
    }
    Measured {
        euler: used.len() as i64 - undirected.len() as i64 + mesh.triangles.len() as i64,
        boundary_edges: boundary,
        area,
        volume,
    }
}

fn settings(block: &str) -> Vec<(String, f64, f64)> {
    block
        .lines()
        .filter(|l| l.starts_with("mesh "))
        .map(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            (
                w[1].to_string(),
                w[2].parse().unwrap(),
                w[3].parse().unwrap(),
            )
        })
        .collect()
}

#[test]
fn every_fixture_mesh_meets_the_reference_expectations() {
    let expected: BTreeMap<&str, Vec<&str>> =
        include_str!("../../fixtures/tessellation-expected.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let row: Vec<&str> = l.split('\t').collect();
                (row[0], row[1..].to_vec())
            })
            .collect();
    let text = include_str!("../../fixtures/tessellation-cases.txt");
    let mut meshes = 0;
    for block in text.split("\nend").filter(|b| !b.trim().is_empty()) {
        let body: String = block
            .trim()
            .lines()
            .filter(|l| !l.starts_with("mesh "))
            .map(|l| format!("{l}\n"))
            .collect();
        let spec = identity_protocol::parse(&body);
        let row = &expected[spec.name.as_str()];
        let solid = row[0] == "solid";
        let (euler, loops): (i64, usize) = (row[1].parse().unwrap(), row[2].parse().unwrap());
        let (area, volume): (f64, f64) = (row[3].parse().unwrap(), row[4].parse().unwrap());
        for (setting, deflection, angle) in settings(block) {
            let p = Parameters::new(deflection, angle).unwrap();
            let name = format!("{} {setting}", spec.name);
            let (topology, mesh) = match spec.make {
                Some(_) => {
                    let body = identity_protocol::build_body(&spec).0;
                    let mesh = body.tessellate(p).unwrap_or_else(|e| panic!("{name}: {e}"));
                    (body.topology().clone(), mesh)
                }
                None => {
                    let s = identity_protocol::build(&spec);
                    let mesh = s.tessellate(p).unwrap_or_else(|e| panic!("{name}: {e}"));
                    (s.topology().clone(), mesh)
                }
            };
            let m = contract(&name, &topology, &mesh, p, solid);
            assert_eq!(m.euler, euler, "{name}");
            if !solid {
                assert!(m.boundary_edges >= 3 * loops, "{name}");
            }
            // The mesh points lie within the deflection of the surface, so
            // the enclosed volume is within it times the areas.
            if solid {
                assert!(
                    (m.volume - volume).abs() <= deflection * (area + m.area) * (1.0 + 1e-9),
                    "{name}: volume {} against {volume}",
                    m.volume
                );
            }
            assert!(
                (m.area - area).abs() <= 0.05 * area,
                "{name}: area {}",
                m.area
            );
            // Determinism: the same mesh again, bit for bit.
            let again = tessellate(&topology, p).unwrap();
            assert_eq!(again, mesh, "{name}");
            meshes += 1;
        }
    }
    assert_eq!(meshes, 56);
}

#[test]
fn finer_settings_converge() {
    // A sphere: the enclosed volume and area approach the exact ones.
    let text = include_str!("../../fixtures/tessellation-cases.txt");
    let block = text
        .split("\nend")
        .find(|b| b.contains("case sphere "))
        .unwrap();
    let body: String = block
        .trim()
        .lines()
        .filter(|l| !l.starts_with("mesh "))
        .map(|l| format!("{l}\n"))
        .collect();
    let solid = identity_protocol::build(&identity_protocol::parse(&body));
    let exact = 4.0 / 3.0 * std::f64::consts::PI * 125.0;
    let mut last = f64::INFINITY;
    for deflection in [0.2, 0.05, 0.0125] {
        let mesh = solid
            .tessellate(Parameters::new(deflection, 0.5).unwrap())
            .unwrap();
        let m = contract(
            "sphere",
            solid.topology(),
            &mesh,
            Parameters::new(deflection, 0.5).unwrap(),
            true,
        );
        let error = (m.volume - exact).abs();
        assert!(error < last, "{deflection}: {error} after {last}");
        last = error;
    }
}

#[test]
fn parameters_are_checked() {
    for (d, a) in [(0.0, 0.5), (-1.0, 0.5), (0.1, 0.0), (0.1, 1.6), (0.1, -0.2)] {
        assert!(
            matches!(Parameters::new(d, a), Err(Error::OutOfDomain(_))),
            "{d} {a}"
        );
    }
    for (d, a) in [(f64::NAN, 0.5), (0.1, f64::INFINITY)] {
        assert!(
            matches!(Parameters::new(d, a), Err(Error::NonFinite(_))),
            "{d} {a}"
        );
    }
    assert!(Parameters::new(1e-3, std::f64::consts::FRAC_PI_2).is_ok());
}

#[test]
fn wire_bodies_give_polylines_and_budgets_are_typed() {
    let specs = identity_protocol::cases(include_str!("../../fixtures/identity-sheet-cases.txt"));
    let wire = specs
        .iter()
        .find(|s| s.make.as_deref() == Some("wire"))
        .unwrap();
    let body = identity_protocol::build_body(wire).0;
    let mesh = body
        .tessellate(Parameters::new(0.01, 0.5).unwrap())
        .unwrap();
    assert!(mesh.triangles.is_empty());
    assert_eq!(mesh.edges.len(), body.topology().edges().len());
    assert!(mesh.edges.iter().all(|e| e.nodes.len() >= 2));
    // A deflection far below the body's size runs out of budget, typed.
    let spec = specs.iter().find(|s| s.name == "face_disc").unwrap();
    let disc = identity_protocol::build_body(spec).0;
    assert!(matches!(
        disc.tessellate(Parameters::new(1e-15, 0.5).unwrap()),
        Err(Error::ComputationLimit(_))
    ));
}

/// T-b: the twelve spline fixture bodies at both settings against the
/// reference's expectations (tessellation-spline-*.{txt,tsv}).
#[test]
fn every_spline_fixture_mesh_meets_the_reference_expectations() {
    let expected: BTreeMap<&str, Vec<&str>> =
        include_str!("../../fixtures/tessellation-spline-expected.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let row: Vec<&str> = l.split('\t').collect();
                (row[0], row[1..].to_vec())
            })
            .collect();
    let text = include_str!("../../fixtures/tessellation-spline-cases.txt");
    let mut meshes = 0;
    for block in text.split("\nend").filter(|b| !b.trim().is_empty()) {
        let body: String = block
            .trim()
            .lines()
            .filter(|l| !l.starts_with("mesh "))
            .map(|l| format!("{l}\n"))
            .collect();
        let (name, tolerance, parts) = brep_protocol::parse(body.trim());
        let topology = Topology::from_parts(parts, Tolerance::new(tolerance, 1e-12).unwrap())
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let row = &expected[name.as_str()];
        let solid = row[0] == "solid";
        let (euler, loops): (i64, usize) = (row[1].parse().unwrap(), row[2].parse().unwrap());
        let (area, volume): (f64, f64) = (row[3].parse().unwrap(), row[4].parse().unwrap());
        for (setting, deflection, angle) in settings(block) {
            let p = Parameters::new(deflection, angle).unwrap();
            let label = format!("{name} {setting}");
            let mesh = tessellate(&topology, p).unwrap_or_else(|e| panic!("{label}: {e}"));
            let m = contract(&label, &topology, &mesh, p, solid);
            assert_eq!(m.euler, euler, "{label}");
            if !solid {
                assert!(m.boundary_edges >= 3 * loops, "{label}");
            } else {
                assert!(
                    (m.volume - volume).abs() <= deflection * (area + m.area) * (1.0 + 1e-9),
                    "{label}: volume {} against {volume}",
                    m.volume
                );
            }
            assert!(
                (m.area - area).abs() <= 0.05 * area,
                "{label}: area {}",
                m.area
            );
            assert_eq!(tessellate(&topology, p).unwrap(), mesh, "{label}");
            meshes += 1;
        }
    }
    assert_eq!(meshes, 24);
}

/// T-b: every valid case of the B-rep fixtures with spline geometry
/// (ranges, reversed spans, rational and unclamped edges, a spline cap,
/// spline walls and holes, a tiny far corner, a spline sheet) meets the
/// contract, a solid's volume within the deflection times the areas of its
/// certified mass enclosure.
#[test]
fn valid_spline_cases_of_the_brep_fixtures_meet_the_contract() {
    let mut meshed = 0;
    for block in include_str!("../../fixtures/brep-cases.txt")
        .split("\nend")
        .filter(|b| b.contains("bspline"))
    {
        let (name, tolerance, parts) = brep_protocol::parse(block.trim());
        let Ok(t) = Topology::from_parts(parts, Tolerance::new(tolerance, 1e-12).unwrap()) else {
            continue;
        };
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in t.vertices().iter().map(|v| v.position) {
            let a = p.to_array();
            for k in 0..3 {
                lo[k] = lo[k].min(a[k]);
                hi[k] = hi[k].max(a[k]);
            }
        }
        let size = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
        let p = Parameters::new(size * 2e-3, 0.5).unwrap();
        let mesh = tessellate(&t, p).unwrap_or_else(|e| panic!("{name}: {e}"));
        let solid = t
            .regions()
            .iter()
            .any(|r| r.kind == rusty_occt::topology::RegionKind::Solid);
        let m = contract(&name, &t, &mesh, p, solid);
        if let (true, Some(mass)) = (solid, t.mass_enclosure()) {
            let [v0, v1] = mass.volume;
            let slack = p.deflection() * (mass.surface_area[1] + m.area) * (1.0 + 1e-9);
            assert!(
                m.volume >= v0 - slack && m.volume <= v1 + slack,
                "{name}: {} not in {:?}",
                m.volume,
                mass.volume
            );
        }
        meshed += 1;
    }
    assert_eq!(meshed, 17);
}

/// The document without the free shapes directly under its root compound
/// (as `occt_brep.rs` does): only the solids are needed here.
fn solids_only(mut doc: rusty_occt::occt_brep::Document) -> rusty_occt::occt_brep::Document {
    use rusty_occt::occt_brep::read::Kind;
    let root = doc.root.shape;
    if doc.shapes[root].kind == Kind::Compound {
        let kept: Vec<_> = doc.shapes[root]
            .subs
            .iter()
            .filter(|s| {
                matches!(
                    doc.shapes[s.shape].kind,
                    Kind::Solid | Kind::CompSolid | Kind::Compound
                )
            })
            .copied()
            .collect();
        doc.shapes[root].subs = kept;
    }
    doc
}

/// Every corpus solid the importer certifies (planes, cylinders, cones,
/// spheres, tori and, since T-b, splines, with general loops, seams merged
/// into windings, holes): closed, oriented, within its bounds, its volume
/// within the deflection times the areas of the certified mass enclosure.
#[test]
fn imported_corpus_solids_meet_the_contract() {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/occ/");
    let mut files: Vec<&str> = include_str!("../../fixtures/brep-io-expected.tsv")
        .lines()
        .filter_map(|l| l.strip_prefix("file\t"))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    files.dedup();
    let mut meshed = 0;
    for name in files {
        let text = std::fs::read_to_string(format!("{root}{name}")).unwrap();
        let im = import(&solids_only(read(&text).unwrap()));
        for solid in &im.solids {
            let Ok(t) = &solid.result else { continue };
            // A spline solid's certified mass enclosure takes minutes: its
            // mesh is checked without the volume band.
            let spline = t
                .faces()
                .iter()
                .any(|f| matches!(f.surface, Surface::BSpline(_)))
                || t.edges()
                    .iter()
                    .any(|e| matches!(e.curve, rusty_occt::topology::Curve3::BSpline(_)));
            let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
            // The box of the edges' points (ring edges have no vertex).
            let points = t
                .edges()
                .iter()
                .flat_map(|e| [0.0, 0.25, 0.5, 0.75].map(|f| e.curve.point(f)));
            for p in points {
                let a = p.to_array();
                for k in 0..3 {
                    lo[k] = lo[k].min(a[k]);
                    hi[k] = hi[k].max(a[k]);
                }
            }
            // A whole sphere or torus has no edge: its volume's cube root.
            let cube = if spline {
                0.0
            } else {
                t.mass_enclosure().map_or(0.0, |m| m.volume[1].abs().cbrt())
            };
            let size = (0..3).map(|k| hi[k] - lo[k]).fold(cube.max(1e-3), f64::max);
            let p = Parameters::new(size * 2e-3, 0.5).unwrap();
            let label = format!("{name} {}", solid.record);
            let mesh = tessellate(t, p).unwrap_or_else(|e| panic!("{label}: {e}"));
            let m = contract(&label, t, &mesh, p, true);
            if let Some(mass) = (!spline).then(|| t.mass_enclosure()).flatten() {
                let [v0, v1] = mass.volume;
                let [_, a1] = mass.surface_area;
                let slack = p.deflection() * (a1 + m.area) * (1.0 + 1e-9);
                assert!(
                    m.volume >= v0 - slack && m.volume <= v1 + slack,
                    "{label}: {} not in {:?}",
                    m.volume,
                    mass.volume
                );
            }
            meshed += 1;
        }
    }
    // The 58 certified corpus solids, the four with spline geometry among
    // them.
    assert_eq!(meshed, 58);
}
