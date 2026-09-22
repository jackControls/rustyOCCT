//! Reproducible generated geometry tests. These are mathematical/metamorphic
//! checks, independent of OCCT and any application document or DTO format.
use rusty_occt::*;
use std::f64::consts::TAU;

struct Generator(u64);
impl Generator {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
    fn signed(&mut self) -> f64 {
        2.0 * self.unit() - 1.0
    }
}

fn close(actual: f64, expected: f64, unit: f64, case: usize, property: &str) {
    assert!(actual.is_finite() && expected.is_finite());
    let budget = 2e-10 * (unit + expected.abs());
    assert!(
        (actual - expected).abs() <= budget,
        "case {case}, {property}: {actual} != {expected}, budget {budget}"
    );
}

#[test]
fn generated_prisms_preserve_geometric_invariants() {
    let mut rng = Generator(0xc0de_5eed_0cc7_2026);
    for case in 0..256 {
        let scale = 2.0_f64.powi((case % 25) as i32 - 12);
        let tolerance = Tolerance::new(1e-9 * scale, 1e-12).unwrap();
        let count = 8 + case % 13;
        let vertices: Vec<_> = (0..count)
            .map(|i| {
                let angle = TAU * i as f64 / count as f64;
                let radius = scale * (0.8 + 0.6 * rng.unit());
                Point2::new(radius * angle.cos(), radius * angle.sin())
            })
            .collect();
        let holes = if case % 2 == 0 {
            vec![]
        } else {
            vec![Boundary::circle(Point2::default(), 0.2 * scale, tolerance).unwrap()]
        };
        let profile = Profile::new(
            Boundary::polygon(vertices.clone(), tolerance).unwrap(),
            holes.clone(),
            tolerance,
        )
        .unwrap();
        let normal = Vec3::new(rng.signed(), rng.signed(), 0.5 + rng.unit());
        let origin = Point3::new(
            10.0 * scale * rng.signed(),
            10.0 * scale * rng.signed(),
            10.0 * scale * rng.signed(),
        );
        let frame = Frame3::new(origin, normal, Vec3::X, tolerance).unwrap();
        let low = -scale * (0.5 + rng.unit());
        let high = scale * (0.5 + rng.unit());
        let solid = Solid::extrude(profile.clone(), frame, low, high).unwrap();
        let mass = solid.mass_properties();
        solid.topology().validate(tolerance).unwrap();
        assert_eq!(
            solid.topology().euler_characteristic(),
            2 - 2 * holes.len() as i64,
            "case {case}"
        );
        assert!(mass.volume > 0.0 && mass.surface_area > 0.0);
        assert!(solid.bounds().contains(mass.centroid, tolerance));

        // A planar cut conserves volume and first moments; exactly two copies
        // of the cut's material area are introduced into total surface area.
        let cut = low + (high - low) * (0.2 + 0.6 * rng.unit());
        let first = Solid::extrude(profile.clone(), frame, low, cut)
            .unwrap()
            .mass_properties();
        let second = Solid::extrude(profile.clone(), frame, cut, high)
            .unwrap()
            .mass_properties();
        close(
            first.volume + second.volume,
            mass.volume,
            scale.powi(3),
            case,
            "volume conservation",
        );
        close(
            first.surface_area + second.surface_area,
            mass.surface_area + 2.0 * profile.area(),
            scale.powi(2),
            case,
            "cut surfaces",
        );
        for axis in 0..3 {
            close(
                first.volume * first.centroid.to_array()[axis]
                    + second.volume * second.centroid.to_array()[axis],
                mass.volume * mass.centroid.to_array()[axis],
                scale.powi(4),
                case,
                "first moment conservation",
            );
        }

        // Permuting the start vertex and reversing winding changes storage and
        // summation order, but must preserve the material region and its mass.
        let mut reordered = vertices;
        reordered.rotate_left(1 + case % (count - 1));
        reordered.reverse();
        let alternate = Solid::extrude(
            Profile::new(
                Boundary::polygon(reordered, tolerance).unwrap(),
                holes,
                tolerance,
            )
            .unwrap(),
            frame,
            high,
            low,
        )
        .unwrap();
        let alternate_mass = alternate.mass_properties();
        close(
            alternate_mass.volume,
            mass.volume,
            scale.powi(3),
            case,
            "winding volume",
        );
        close(
            alternate_mass.surface_area,
            mass.surface_area,
            scale.powi(2),
            case,
            "winding area",
        );

        let rotation = RigidTransform::rotation(
            Point3::ORIGIN,
            Vec3::new(0.3 + rng.unit(), rng.signed(), rng.signed()),
            TAU * rng.signed(),
        )
        .unwrap();
        let transform = rotation
            .then(
                RigidTransform::translation(Vec3::new(20.0 * scale, -40.0 * scale, 30.0 * scale))
                    .unwrap(),
            )
            .unwrap();
        let moved = solid.transformed(transform).unwrap();
        let moved_mass = moved.mass_properties();
        close(
            moved_mass.volume,
            mass.volume,
            scale.powi(3),
            case,
            "rigid volume",
        );
        close(
            moved_mass.surface_area,
            mass.surface_area,
            scale.powi(2),
            case,
            "rigid area",
        );
        let expected_centroid = transform.point(mass.centroid).to_array();
        for (actual, expected) in moved_mass
            .centroid
            .to_array()
            .into_iter()
            .zip(expected_centroid)
        {
            close(actual, expected, scale, case, "rigid centroid");
        }

        // Central inertia transforms as R I R^T; translations do not enter.
        let columns = [
            rotation.vector(Vec3::X).to_array(),
            rotation.vector(Vec3::Y).to_array(),
            rotation.vector(Vec3::Z).to_array(),
        ];
        for i in 0..3 {
            for j in 0..3 {
                let mut expected = 0.0;
                for a in 0..3 {
                    for b in 0..3 {
                        expected += columns[a][i] * mass.inertia[a][b] * columns[b][j];
                    }
                }
                close(
                    moved_mass.inertia[i][j],
                    expected,
                    scale.powi(5),
                    case,
                    "rigid inertia",
                );
                close(
                    alternate_mass.inertia[i][j],
                    mass.inertia[i][j],
                    scale.powi(5),
                    case,
                    "winding inertia",
                );
            }
        }
        for _ in 0..32 {
            let local = Point2::new(1.6 * scale * rng.signed(), 1.6 * scale * rng.signed());
            let height = (low + high) * 0.5 + (high - low) * rng.signed();
            let point = frame.point(local, height);
            let classification = solid.classify(point).unwrap();
            assert_eq!(
                alternate.classify(point).unwrap(),
                classification,
                "case {case}: winding classification"
            );
            assert_eq!(
                moved.classify(transform.point(point)).unwrap(),
                classification,
                "case {case}: rigid classification"
            );
            let direction = [rng.signed(), rng.signed(), rng.signed()];
            let mut quadratic_form = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    quadratic_form += direction[i] * mass.inertia[i][j] * direction[j];
                }
            }
            assert!(
                quadratic_form >= 0.0,
                "case {case}: inertia must be positive semidefinite"
            );
        }
        // Analytic edge samples must lie in both the original and moved bounds.
        for shape in [&solid, &moved] {
            for edge in shape.topology().edges() {
                for sample in 0..=16 {
                    let p = edge.curve.point(sample as f64 / 16.0);
                    assert!(
                        shape.bounds().contains(p, tolerance),
                        "case {case}: bounds missed edge"
                    );
                    assert_eq!(
                        shape.classify(p).unwrap(),
                        Location::Boundary,
                        "case {case}: edge classification"
                    );
                }
            }
        }
    }
}
