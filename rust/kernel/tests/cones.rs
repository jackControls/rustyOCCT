//! The cone builder (S3 of REVIEW_NOTES.md) against the independent
//! primitive reference: OCCT's synthesized counts, certified mass
//! properties, faces (type, area, centre), ring edges (length, point at the
//! middle parameter) and apex vertices, for every cone of
//! `primitive-cases.txt`. The kernel stores rounded frame axes, angle and
//! slant, so numbers compare within 1e-12 of the case's scale; the mass
//! enclosures themselves are certified for the stored geometry.
use rusty_occt::identity::OperationId;
use rusty_occt::topology::{Curve3, FaceId, Surface};
use rusty_occt::{Frame3, Point3, Solid, Tolerance, Vec3};
use std::collections::BTreeMap;

const BOUND: f64 = 1e-12;

struct Expected {
    counts: Vec<usize>,
    props: Vec<f64>,
    faces: Vec<(String, f64, [f64; 3])>,
    edges: Vec<(String, String, f64, [f64; 3])>,
    vertices: Vec<[f64; 3]>,
}

fn expected() -> BTreeMap<String, Expected> {
    let mut out: BTreeMap<String, Expected> = BTreeMap::new();
    for line in include_str!("../../fixtures/primitive-expected.tsv").lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut w = line.split('\t');
        let (name, kind, values) = (w.next().unwrap(), w.next().unwrap(), w.next().unwrap());
        let v: Vec<&str> = values.split_whitespace().collect();
        let e = out.entry(name.to_string()).or_insert_with(|| Expected {
            counts: vec![],
            props: vec![],
            faces: vec![],
            edges: vec![],
            vertices: vec![],
        });
        let f = |s: &str| s.parse::<f64>().unwrap();
        let p = |a: &[&str]| [f(a[0]), f(a[1]), f(a[2])];
        match kind {
            "counts" => e.counts = v.iter().map(|x| x.parse().unwrap()).collect(),
            "props" => e.props = v.iter().map(|x| f(x)).collect(),
            "face" => e.faces.push((v[0].into(), f(v[1]), p(&v[2..5]))),
            "edge" => e
                .edges
                .push((v[0].into(), v[1].into(), f(v[2]), p(&v[3..6]))),
            _ => e.vertices.push(p(&v[..3])),
        }
    }
    out
}

fn near(a: [f64; 3], b: [f64; 3], scale: f64) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() <= BOUND * scale)
}

#[test]
fn cones_match_the_independent_reference() {
    let want = expected();
    let mut checked = 0;
    for line in include_str!("../../fixtures/primitive-cases.txt").lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        let n: Vec<f64> = w[2..].iter().map(|x| x.parse().unwrap()).collect();
        let name = w[1];
        let e = &want[name];
        let tol = Tolerance::default();
        let frame = Frame3::new(
            Point3::new(n[0], n[1], n[2]),
            Vec3::new(n[3], n[4], n[5]),
            Vec3::new(n[6], n[7], n[8]),
            tol,
        )
        .unwrap();
        let (solid, _) = Solid::cone_with(OperationId(1), frame, n[9], n[10], n[11], tol)
            .unwrap_or_else(|err| panic!("{name}: {err}"));
        let t = solid.topology();
        assert!(t.check(tol).is_empty(), "{name}: {:?}", t.check(tol));
        let size = [n[0], n[1], n[2], n[9], n[10], n[11], 1.0]
            .iter()
            .fold(0.0f64, |m, x| m.max(x.abs()));
        let c = t.occt_counts();
        assert_eq!(
            vec![c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids],
            e.counts,
            "{name}"
        );
        // Certified enclosure within the bound of the exact reference.
        let m = t.mass_enclosure().expect("certified");
        let (volume, area) = (e.props[0], e.props[1]);
        let within = |x: [f64; 2], v: f64, scale: f64| {
            x[0] - BOUND * scale <= v && v <= x[1] + BOUND * scale
        };
        assert!(
            within(m.volume, volume, volume.abs().max(size.powi(3))),
            "{name} volume {m:?}"
        );
        assert!(
            within(m.surface_area, area, area.max(size * size)),
            "{name} area"
        );
        for i in 0..3 {
            assert!(
                within(m.centroid[i], e.props[2 + i], size),
                "{name} centroid"
            );
        }
        let entries = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
        for (k, (a, b)) in entries.iter().enumerate() {
            assert!(
                within(m.inertia[*a][*b], e.props[5 + k], size.powi(5)),
                "{name} inertia {a}{b}"
            );
        }
        assert_eq!(solid.mass_properties(), m.midpoints());
        // Faces by type, area and centre.
        let mut faces: Vec<_> = e.faces.clone();
        for f in 0..t.faces().len() {
            let kind = match t.faces()[f].surface {
                Surface::Plane(_) => "plane",
                Surface::Cone { .. } => "cone",
                Surface::Cylinder { .. } => "cylinder",
            };
            let (a, centre) = t.face_area_and_centre(FaceId::new(f)).expect("face");
            let k = faces
                .iter()
                .position(|(k, area, c)| {
                    k == kind
                        && (a - area).abs() <= BOUND * area.max(size * size)
                        && near(centre.to_array(), *c, size)
                })
                .unwrap_or_else(|| panic!("{name}: face {f} {kind} {a} {centre:?}"));
            faces.remove(k);
        }
        assert!(faces.is_empty(), "{name}: faces left {faces:?}");
        // Ring edges (OCCT's seam and degenerated edge are its structure).
        for edge in t.edges() {
            let Curve3::Circle { radius, .. } = edge.curve else {
                panic!("{name}: a cone's edges are rings");
            };
            let mid = edge.curve.point(0.5).to_array();
            let length = std::f64::consts::TAU * radius;
            assert!(
                e.edges.iter().any(|(d, closed, l, m)| d == "regular"
                    && closed == "closed"
                    && (length - l).abs() <= BOUND * l.max(size)
                    && near(mid, *m, size)),
                "{name}: ring {length} {mid:?}"
            );
        }
        for v in t.vertices() {
            assert!(
                e.vertices
                    .iter()
                    .any(|p| near(v.position.to_array(), *p, size)),
                "{name}: apex {:?}",
                v.position
            );
        }
        checked += 1;
    }
    assert_eq!(checked, want.len());
}
