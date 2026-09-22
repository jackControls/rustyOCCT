#![no_main]
use rusty_occt::*;
use rusty_occt_fuzz::{byte, decode};
use std::f64::consts::TAU;

fn close(a: f64, b: f64, unit: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!((a - b).abs() <= 2e-9 * (unit + b.abs()), "{a} != {b}");
}

fn validate(solid: &Solid, tolerance: Tolerance, scale: f64) {
    solid.topology().validate(tolerance).unwrap();
    let mass = solid.mass_properties();
    assert!(mass.volume.is_finite() && mass.volume > 0.);
    assert!(mass.surface_area.is_finite() && mass.surface_area > 0.);
    for i in 0..3 {
        assert!(mass.centroid.to_array()[i].is_finite());
        assert!(mass.inertia[i][i].is_finite() && mass.inertia[i][i] >= 0.);
        for j in 0..3 {
            close(mass.inertia[i][j], mass.inertia[j][i], scale.powi(5));
        }
    }
    assert_eq!(
        solid.topology().euler_characteristic(),
        2 - 2 * solid.profile().holes().len() as i64
    );
    for edge in solid.topology().edges() {
        for t in [0., 0.5, 1.] {
            let p = edge.curve.point(t);
            assert!(solid.bounds().contains(p, tolerance));
            assert_eq!(solid.classify(p).unwrap(), Location::Boundary);
        }
    }
}

libfuzzer_sys::fuzz_target!(|data: &[u8]| {
    // Invalid/unbounded numeric input is a separate mode; accepted solids must
    // still satisfy their public invariants. No errors are converted to passes
    // in the valid structured mode below.
    if byte(data, 0) == 255 {
        let mut raw = data.to_vec();
        raw[0] = 0;
        let p = decode(&raw);
        let tol = Tolerance::default();
        if let Ok(boundary) =
            Boundary::polygon(p[..3].iter().map(|p| Point2::new(p.x, p.y)).collect(), tol)
        {
            if let Ok(profile) = Profile::new(boundary, vec![], tol) {
                if let Ok(solid) = Solid::extrude(profile, Frame3::xy(), p[3].z, p[4].z) {
                    solid.topology().validate(tol).unwrap();
                    let mass = solid.mass_properties();
                    assert!(mass.volume.is_finite() && mass.volume > 0.);
                    assert!(mass.surface_area.is_finite() && mass.surface_area > 0.);
                }
            }
        }
        return;
    }
    let scale = 2.0_f64.powi((byte(data, 1) % 49) as i32 - 24);
    let tolerance = Tolerance::new(scale * 1e-9, 1e-12).unwrap();
    let count = 3 + (byte(data, 2) % 22) as usize;
    let vertices: Vec<_> = (0..count)
        .map(|i| {
            let angle = TAU * i as f64 / count as f64;
            let radius = scale * (1. + byte(data, 4 + i) as f64 / 255.);
            Point2::new(radius * angle.cos(), radius * angle.sin())
        })
        .collect();
    let outer = if byte(data, 3).is_multiple_of(3) {
        Boundary::circle(Point2::default(), 2. * scale, tolerance).unwrap()
    } else {
        Boundary::polygon(vertices, tolerance).unwrap()
    };
    let holes = if byte(data, 3).is_multiple_of(2) {
        vec![]
    } else {
        vec![Boundary::circle(Point2::default(), 0.1 * scale, tolerance).unwrap()]
    };
    let profile = Profile::new(outer, holes, tolerance).unwrap();
    let start = -scale;
    let end = scale * (1. + byte(data, 28) as f64 / 255.);
    let mut solid = Solid::extrude(profile, Frame3::xy(), start, end).unwrap();
    let original = solid.mass_properties();
    for step in 0..1 + (byte(data, 29) % 8) as usize {
        let offset = 30 + 8 * step;
        let unit = |i| byte(data, offset + i) as f64 / 255.;
        match byte(data, offset) % 4 {
            0 | 1 => {
                let transform = if byte(data, offset).is_multiple_of(4) {
                    RigidTransform::rotation(
                        Point3::ORIGIN,
                        Vec3::new(0.25 + unit(1), unit(2) - 0.5, unit(3) - 0.5),
                        TAU * unit(4),
                    )
                    .unwrap()
                } else {
                    RigidTransform::translation(Vec3::new(
                        scale * (unit(1) - 0.5),
                        scale * (unit(2) - 0.5),
                        scale * (unit(3) - 0.5),
                    ))
                    .unwrap()
                };
                let centroid = transform.point(solid.mass_properties().centroid);
                // Probe chosen well away from all boundary tolerance bands.
                let probe = solid.frame().point(Point2::new(0.25 * scale, 0.), 0.);
                let classification = solid.classify(probe).unwrap();
                solid = solid.transformed(transform).unwrap();
                assert_eq!(
                    solid.classify(transform.point(probe)).unwrap(),
                    classification
                );
                for (a, b) in solid
                    .mass_properties()
                    .centroid
                    .to_array()
                    .into_iter()
                    .zip(centroid.to_array())
                {
                    close(a, b, scale);
                }
            }
            2 => {
                let cut = start + (end - start) * (0.2 + 0.6 * unit(1));
                let left = Solid::extrude(solid.profile().clone(), solid.frame(), start, cut)
                    .unwrap()
                    .mass_properties();
                let right = Solid::extrude(solid.profile().clone(), solid.frame(), cut, end)
                    .unwrap()
                    .mass_properties();
                close(left.volume + right.volume, original.volume, scale.powi(3));
                close(
                    left.surface_area + right.surface_area,
                    original.surface_area + 2. * solid.profile().area(),
                    scale.powi(2),
                );
                for i in 0..3 {
                    close(
                        left.volume * left.centroid.to_array()[i]
                            + right.volume * right.centroid.to_array()[i],
                        original.volume * solid.mass_properties().centroid.to_array()[i],
                        scale.powi(4),
                    );
                }
            }
            _ => {
                let mut profile = solid.profile().clone();
                if let Some(vertices) = profile.outer().polygon_vertices() {
                    let mut reordered = vertices.to_vec();
                    reordered.reverse();
                    reordered.rotate_left(1);
                    profile = Profile::new(
                        Boundary::polygon(reordered, tolerance).unwrap(),
                        profile.holes().to_vec(),
                        tolerance,
                    )
                    .unwrap();
                }
                solid = Solid::extrude(profile, solid.frame(), end, start).unwrap();
            }
        }
        validate(&solid, tolerance, scale);
        close(
            solid.mass_properties().volume,
            original.volume,
            scale.powi(3),
        );
        close(
            solid.mass_properties().surface_area,
            original.surface_area,
            scale.powi(2),
        );
    }
});
