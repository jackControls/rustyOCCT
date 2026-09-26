use rusty_occt::identity::OperationId;
use rusty_occt::topology::{FaceId, FaceOrigin, Loop, Surface};
use rusty_occt::*;
use std::f64::consts::{FRAC_PI_2, PI};

fn close(actual: f64, expected: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn point(actual: Point3, expected: Point3) {
    for (a, e) in actual.to_array().into_iter().zip(expected.to_array()) {
        close(a, e);
    }
}
fn polygon(points: &[[f64; 2]]) -> Boundary {
    Boundary::polygon(
        points.iter().map(|p| Point2::new(p[0], p[1])).collect(),
        Tolerance::default(),
    )
    .unwrap()
}
fn profile(outer: Boundary, holes: Vec<Boundary>) -> Profile {
    Profile::new(outer, holes, Tolerance::default()).unwrap()
}

#[test]
fn box_has_exact_properties_and_closed_oriented_topology() {
    let t = Tolerance::default();
    let solid = Solid::cuboid(20.0, 30.0, 40.0, t).unwrap();
    let mass = solid.mass_properties();
    close(mass.volume, 24000.0);
    close(mass.surface_area, 5200.0);
    point(mass.centroid, Point3::new(10.0, 15.0, 20.0));
    close(
        mass.inertia[0][0],
        24000.0 * (30.0 * 30.0 + 40.0 * 40.0) / 12.0,
    );
    close(
        mass.inertia[1][1],
        24000.0 * (20.0 * 20.0 + 40.0 * 40.0) / 12.0,
    );
    close(
        mass.inertia[2][2],
        24000.0 * (20.0 * 20.0 + 30.0 * 30.0) / 12.0,
    );
    close(mass.inertia[0][1], 0.0);
    let topology = solid.topology();
    assert_eq!(
        (
            topology.vertices().len(),
            topology.edges().len(),
            topology.faces().len()
        ),
        (8, 12, 6)
    );
    assert_eq!(topology.euler_characteristic(), 2);
    topology.validate(t).unwrap();
    for edge in topology.edge_ids() {
        assert_eq!(topology.incident_faces(edge).unwrap().len(), 2);
    }
    assert_eq!(
        solid.classify(Point3::new(1.0, 2.0, 3.0)).unwrap(),
        Location::Inside
    );
    assert_eq!(
        solid.classify(Point3::new(0.0, 2.0, 3.0)).unwrap(),
        Location::Boundary
    );
    assert_eq!(
        solid.classify(Point3::new(0.0, 2.0, 41.0)).unwrap(),
        Location::Outside
    );
    point(solid.bounds().min, Point3::ORIGIN);
    point(solid.bounds().max, Point3::new(20.0, 30.0, 40.0));
}

#[test]
fn signed_box_dimensions_follow_occt_minimum_corner_convention() {
    let t = Tolerance::default();
    for x in [-1.0, 1.0] {
        for y in [-1.0, 1.0] {
            for z in [-1.0, 1.0] {
                let solid = Solid::box_at(
                    Point3::new(5.0, 7.0, 9.0),
                    Vec3::new(x * 2.0, y * 3.0, z * 4.0),
                    t,
                )
                .unwrap();
                point(
                    solid.bounds().min,
                    Point3::new(
                        5.0 + (x * 2.0).min(0.0),
                        7.0 + (y * 3.0).min(0.0),
                        9.0 + (z * 4.0).min(0.0),
                    ),
                );
                close(solid.mass_properties().volume, 24.0);
                solid.topology().validate(t).unwrap();
            }
        }
    }
    assert!(Solid::box_at(Point3::ORIGIN, Vec3::new(0.0, 1.0, 1.0), t).is_err());
    assert!(Solid::box_at(Point3::ORIGIN, Vec3::new(f64::NAN, 1.0, 1.0), t).is_err());
}

#[test]
fn aabb_overlap_is_symmetric_and_respects_both_gaps() {
    let t = Tolerance::default();
    let a = Bounds3 {
        min: Point3::ORIGIN,
        max: Point3::new(1.0, 1.0, 1.0),
    };
    for axis in 0..3 {
        for (gap, expected) in [(0.0, true), (1.9e-7, true), (2.1e-7, false), (1.0, false)] {
            let mut lo = [0.0; 3];
            let mut hi = [1.0; 3];
            lo[axis] = 1.0 + gap;
            hi[axis] = 2.0 + gap;
            let b = Bounds3 {
                min: Point3::new(lo[0], lo[1], lo[2]),
                max: Point3::new(hi[0], hi[1], hi[2]),
            };
            assert_eq!(a.intersects(b, t).unwrap(), expected);
            assert_eq!(b.intersects(a, t).unwrap(), expected);
        }
    }
    assert!(a.intersects(a, t).unwrap());
    assert!(a
        .intersects(
            Bounds3 {
                min: a.max,
                max: a.min
            },
            t
        )
        .is_err());
    assert!(a
        .intersects(
            Bounds3 {
                min: Point3::new(f64::NAN, 0.0, 0.0),
                max: a.max
            },
            t
        )
        .is_err());
}

#[test]
fn concave_profile_preserves_the_missing_corner() {
    let outer = polygon(&[[0., 0.], [4., 0.], [4., 1.], [1., 1.], [1., 4.], [0., 4.]]);
    let solid = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile(outer, vec![]),
        Frame3::xy(),
        -2.,
        3.,
    )
    .map(|(s, _)| s)
    .unwrap();
    let mass = solid.mass_properties();
    close(mass.volume, 35.0);
    close(mass.surface_area, 94.0);
    point(mass.centroid, Point3::new(9.5 / 7.0, 9.5 / 7.0, 0.5));
    assert_eq!(
        solid.classify(Point3::new(2., 2., 0.)).unwrap(),
        Location::Outside
    );
    assert_eq!(
        solid.classify(Point3::new(0.5, 2., 0.)).unwrap(),
        Location::Inside
    );
    assert_eq!(solid.topology().faces().len(), 8);
}

#[test]
fn polygon_hole_is_real_topology_and_subtracts_mass() {
    let p = profile(
        polygon(&[[0., 0.], [10., 0.], [10., 8.], [0., 8.]]),
        vec![polygon(&[[2., 2.], [4., 2.], [4., 4.], [2., 4.]])],
    );
    let solid = Solid::extrude_with(OperationId::UNSPECIFIED, p, Frame3::xy(), 0., 5.)
        .map(|(s, _)| s)
        .unwrap();
    let mass = solid.mass_properties();
    close(mass.volume, 380.0);
    close(mass.surface_area, 372.0);
    point(mass.centroid, Point3::new(388.0 / 76.0, 308.0 / 76.0, 2.5));
    assert_eq!(
        solid.classify(Point3::new(3., 3., 2.)).unwrap(),
        Location::Outside
    );
    assert_eq!(
        solid.classify(Point3::new(2., 3., 2.)).unwrap(),
        Location::Boundary
    );
    assert_eq!(
        solid.classify(Point3::new(3., 3., 0.)).unwrap(),
        Location::Outside
    );
    assert_eq!(solid.topology().euler_characteristic(), 0);
    assert_eq!(solid.topology().faces()[0].loops.len(), 2);
}

#[test]
fn cylinder_is_seamless_with_ring_edges_and_wound_loops() {
    let solid = Solid::cylinder(Frame3::xy(), 3., 0., 8., Tolerance::default()).unwrap();
    let mass = solid.mass_properties();
    close(mass.volume, 72. * PI);
    close(mass.surface_area, 66. * PI);
    close(mass.inertia[2][2], 0.5 * mass.volume * 9.);
    close(mass.inertia[0][0], mass.volume * (3. * 9. + 64.) / 12.);
    let topology = solid.topology();
    // No seam and no vertex: two ring edges bound the three faces.
    assert_eq!(
        (
            topology.vertices().len(),
            topology.edges().len(),
            topology.faces().len()
        ),
        (0, 2, 3)
    );
    assert!(topology
        .edges()
        .iter()
        .all(|e| e.is_ring() && e.fins.len() == 2));
    let face = &topology.faces()[2];
    assert!(matches!(face.surface, Surface::Cylinder { .. }));
    let windings: Vec<[i32; 2]> = face
        .loops
        .iter()
        .map(|l| match &topology.loops()[l.index()] {
            Loop::Edges { winding, .. } => *winding,
            Loop::Vertex(_) => panic!("no vertex loops"),
        })
        .collect();
    assert_eq!(windings, vec![[1, 0], [-1, 0]]);
    // OCCT reports the same body with a seam edge and two seam vertices.
    let counts = topology.occt_counts();
    assert_eq!(
        (
            counts.vertices,
            counts.edges,
            counts.wires,
            counts.faces,
            counts.shells,
            counts.solids
        ),
        (2, 3, 3, 3, 1, 1)
    );
    close(face.normal(Point2::new(0., 4.)).x, 1.);
    assert_eq!(
        solid.classify(Point3::new(3., 0., 4.)).unwrap(),
        Location::Boundary
    );
    assert_eq!(
        solid.classify(Point3::new(2.9, 0., 4.)).unwrap(),
        Location::Inside
    );
}

#[test]
fn tube_has_inward_hole_wall_and_annular_caps() {
    let t = Tolerance::default();
    let outer = Boundary::circle(Point2::default(), 5., t).unwrap();
    let hole = Boundary::circle(Point2::default(), 3., t).unwrap();
    let solid = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile(outer, vec![hole]),
        Frame3::xy(),
        0.,
        10.,
    )
    .map(|(s, _)| s)
    .unwrap();
    let mass = solid.mass_properties();
    close(mass.volume, 160. * PI);
    close(mass.surface_area, 192. * PI);
    close(mass.inertia[2][2], mass.volume * (25. + 9.) / 2.);
    close(
        mass.inertia[0][0],
        mass.volume * (3. * (25. + 9.) + 100.) / 12.,
    );
    assert_eq!(solid.topology().euler_characteristic(), 0);
    close(
        solid.topology().faces()[3].normal(Point2::new(0., 5.)).x,
        -1.,
    );
    assert_eq!(
        solid.classify(Point3::new(0., 0., 5.)).unwrap(),
        Location::Outside
    );
    assert_eq!(
        solid.classify(Point3::new(4., 0., 5.)).unwrap(),
        Location::Inside
    );
}

#[test]
fn plate_with_mixed_holes_has_correct_genus_and_centroid() {
    let t = Tolerance::default();
    let outer = Boundary::rectangle(40., 20., t).unwrap();
    let holes = vec![
        Boundary::circle(Point2::new(10., 10.), 2., t).unwrap(),
        polygon(&[[28., 8.], [32., 8.], [32., 12.], [28., 12.]]),
    ];
    let solid = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile(outer, holes),
        Frame3::xy(),
        -2.,
        2.,
    )
    .map(|(s, _)| s)
    .unwrap();
    let area = 800. - 4. * PI - 16.;
    close(solid.mass_properties().volume, 4. * area);
    point(
        solid.mass_properties().centroid,
        Point3::new((16000. - 40. * PI - 480.) / area, 10., 0.),
    );
    assert_eq!(solid.topology().euler_characteristic(), -2);
    assert_eq!(solid.topology().faces().len(), 11);
}

#[test]
fn winding_and_duplicate_closure_do_not_change_the_solid() {
    let a = polygon(&[[0., 0.], [4., 0.], [4., 3.], [0., 3.]]);
    let b = polygon(&[[0., 0.], [0., 3.], [4., 3.], [4., 0.], [0., 0.]]);
    assert_eq!(a, b);
    let a = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile(a, vec![]),
        Frame3::xy(),
        -3.,
        7.,
    )
    .map(|(s, _)| s)
    .unwrap();
    let b = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile(b, vec![]),
        Frame3::xy(),
        7.,
        -3.,
    )
    .map(|(s, _)| s)
    .unwrap();
    assert_eq!(a.mass_properties(), b.mass_properties());
    assert_eq!(a.bounds(), b.bounds());
    assert_eq!(
        a.topology().face_origin(FaceId::new(0)),
        Some(FaceOrigin::StartCap)
    );
    assert_eq!(
        b.topology().face_origin(FaceId::new(0)),
        Some(FaceOrigin::EndCap)
    );
    close(b.topology().faces()[0].normal(Point2::default()).z, -1.);
}

#[test]
fn arbitrary_plane_and_rigid_motion_preserve_geometry() {
    let t = Tolerance::default();
    let solid = Solid::cuboid(2., 4., 6., t).unwrap();
    let rotation = RigidTransform::rotation(Point3::ORIGIN, Vec3::Y, FRAC_PI_2).unwrap();
    let transform = rotation
        .then(RigidTransform::translation(Vec3::new(10., 20., 30.)).unwrap())
        .unwrap();
    let moved = solid
        .transform_with(OperationId::UNSPECIFIED, transform)
        .map(|(s, _)| s)
        .unwrap();
    point(moved.mass_properties().centroid, Point3::new(13., 22., 29.));
    point(moved.bounds().min, Point3::new(10., 20., 28.));
    point(moved.bounds().max, Point3::new(16., 24., 30.));
    close(
        moved.mass_properties().inertia[0][0],
        solid.mass_properties().inertia[2][2],
    );
    close(
        moved.mass_properties().inertia[2][2],
        solid.mass_properties().inertia[0][0],
    );
    assert_eq!(
        moved
            .classify(transform.point(Point3::new(1., 1., 1.)))
            .unwrap(),
        Location::Inside
    );
    let restored = moved
        .transform_with(OperationId::UNSPECIFIED, transform.inverse().unwrap())
        .map(|(s, _)| s)
        .unwrap();
    point(
        restored.mass_properties().centroid,
        solid.mass_properties().centroid,
    );
    let origins = |s: &Solid| {
        s.topology()
            .face_ids()
            .map(|f| s.topology().face_origin(f))
            .collect::<Vec<_>>()
    };
    assert_eq!(origins(&restored), origins(&solid));
    // Rigid motion keeps every id and derivation.
    for moved in [&moved, &restored] {
        assert_eq!(moved.topology().body_id(), solid.topology().body_id());
        assert!(moved.topology().ids().eq(solid.topology().ids()));
    }
}

#[test]
fn tilted_circle_bounds_are_analytic() {
    let t = Tolerance::default();
    let frame = Frame3::new(
        Point3::new(10., 20., 30.),
        Vec3::new(1., 1., 1.),
        Vec3::X,
        t,
    )
    .unwrap();
    let solid = Solid::cylinder(frame, 4., -2., 7., t).unwrap();
    let extent = 4. * (2_f64 / 3.).sqrt();
    for (i, origin) in [10., 20., 30.].into_iter().enumerate() {
        close(
            solid.bounds().min.to_array()[i],
            origin - 2. / 3_f64.sqrt() - extent,
        );
        close(
            solid.bounds().max.to_array()[i],
            origin + 7. / 3_f64.sqrt() + extent,
        );
    }
}

#[test]
fn large_translation_does_not_corrupt_local_moments() {
    let p = polygon(&[
        [1e6, 1e6],
        [1e6 + 2., 1e6],
        [1e6 + 2., 1e6 + 4.],
        [1e6, 1e6 + 4.],
    ]);
    let solid = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile(p, vec![]),
        Frame3::xy(),
        0.,
        6.,
    )
    .map(|(s, _)| s)
    .unwrap();
    let reference = Solid::cuboid(2., 4., 6., Tolerance::default()).unwrap();
    close(solid.mass_properties().volume, 48.);
    for i in 0..3 {
        for j in 0..3 {
            close(
                solid.mass_properties().inertia[i][j],
                reference.mass_properties().inertia[i][j],
            );
        }
    }
}

#[test]
fn rejects_crossing_touching_and_backtracking_polygons() {
    let t = Tolerance::default();
    for points in [
        vec![[0., 0.], [3., 3.], [0., 3.], [3., 0.]],
        vec![[0., 0.], [3., 0.], [1., 0.], [1., 3.], [0., 3.]],
        vec![[0., 0.], [3., 0.], [3., 3.], [1., 0.], [0., 3.]],
        vec![[0., 0.], [0., 0.], [3., 0.], [0., 3.]],
    ] {
        assert!(Boundary::polygon(
            points
                .into_iter()
                .map(|p| Point2::new(p[0], p[1]))
                .collect(),
            t
        )
        .is_err());
    }
    assert!(Boundary::polygon(
        vec![
            Point2::new(0., 0.),
            Point2::new(1., 1.),
            Point2::new(2., 2.)
        ],
        t
    )
    .is_err());
}

#[test]
fn accepts_collinear_forward_edges_without_losing_provenance() {
    let p = polygon(&[[0., 0.], [2., 0.], [4., 0.], [4., 3.], [0., 3.]]);
    let solid = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile(p, vec![]),
        Frame3::xy(),
        0.,
        5.,
    )
    .map(|(s, _)| s)
    .unwrap();
    close(solid.mass_properties().volume, 60.);
    assert_eq!(solid.topology().faces().len(), 7);
}

#[test]
fn rejects_tangent_crossing_external_and_nested_holes() {
    let t = Tolerance::default();
    let rectangle = || Boundary::rectangle(10., 10., t).unwrap();
    let circle = |x, y, r| Boundary::circle(Point2::new(x, y), r, t).unwrap();
    for hole in [
        circle(1., 5., 1.),
        circle(0., 5., 2.),
        circle(20., 20., 1.),
        circle(5., 5., 20.),
    ] {
        assert!(matches!(
            Profile::new(rectangle(), vec![hole], t),
            Err(Error::InvalidHole(0))
        ));
    }
    assert!(Profile::new(rectangle(), vec![circle(3., 5., 2.), circle(7., 5., 2.)], t).is_err());
    assert!(Profile::new(rectangle(), vec![circle(5., 5., 3.), circle(5., 5., 1.)], t).is_err());
    let hole = polygon(&[[2., 2.], [8., 2.], [8., 8.], [2., 8.]]);
    assert!(Profile::new(rectangle(), vec![hole.clone(), circle(5., 5., 1.)], t).is_err());
    assert!(Profile::new(rectangle(), vec![hole, circle(8., 5., 1.)], t).is_err());
}

#[test]
fn polygon_hole_in_circle_and_circle_hole_in_concave_polygon() {
    let t = Tolerance::default();
    let p = profile(
        Boundary::circle(Point2::default(), 5., t).unwrap(),
        vec![polygon(&[[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])],
    );
    close(
        Solid::extrude_with(OperationId::UNSPECIFIED, p, Frame3::xy(), 0., 1.)
            .map(|(s, _)| s)
            .unwrap()
            .mass_properties()
            .volume,
        25. * PI - 4.,
    );
    let concave = polygon(&[[0., 0.], [4., 0.], [4., 1.], [1., 1.], [1., 4.], [0., 4.]]);
    assert!(Profile::new(
        concave,
        vec![Boundary::circle(Point2::new(0.7, 0.7), 0.5, t).unwrap()],
        t
    )
    .is_err());
}

#[test]
fn tolerances_are_explicit_and_revalidated_across_boundaries() {
    let small = Tolerance::new(1e-10, 1e-12).unwrap();
    let large = Tolerance::new(1e-3, 1e-12).unwrap();
    let outer = Boundary::rectangle(1., 1., small).unwrap();
    let near = Boundary::circle(Point2::new(0.10001, 0.5), 0.1, small).unwrap();
    assert!(Profile::new(outer.clone(), vec![near.clone()], small).is_ok());
    assert!(Profile::new(outer, vec![near], large).is_err());
    let narrow = Boundary::rectangle(0.0001, 1., small).unwrap();
    assert!(Profile::new(narrow, vec![], large).is_err());
    let solid = Solid::cuboid(1., 1., 1., large).unwrap();
    assert_eq!(
        solid.classify(Point3::new(-0.0005, 0.5, 0.5)).unwrap(),
        Location::Boundary
    );
    assert_eq!(
        solid.classify(Point3::new(-0.002, 0.5, 0.5)).unwrap(),
        Location::Outside
    );
}

#[test]
fn rejects_nonfinite_degenerate_and_unresolvable_inputs() {
    let t = Tolerance::default();
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(Boundary::circle(Point2::default(), value, t).is_err());
        assert!(Boundary::rectangle(value, 2., t).is_err());
        assert!(Solid::cylinder(Frame3::xy(), 2., 0., value, t).is_err());
        assert!(Tolerance::new(value, 1e-12).is_err());
        assert!(RigidTransform::rotation(Point3::ORIGIN, Vec3::X, value).is_err());
        assert!(Solid::cuboid(1., 1., 1., t)
            .unwrap()
            .classify(Point3::new(value, 0., 0.))
            .is_err());
    }
    assert!(Tolerance::new(0., 1e-12).is_err());
    assert!(Tolerance::new(1e-7, 1.).is_err());
    assert!(Boundary::circle(Point2::default(), -1., t).is_err());
    assert!(Solid::cuboid(1., 1., 0., t).is_err());
    assert!(Solid::cylinder(Frame3::xy(), 1., 2., 2., t).is_err());
    assert!(Frame3::new(Point3::ORIGIN, Vec3::Z, Vec3::Z, t).is_err());
    assert!(Frame3::new(Point3::ORIGIN, Vec3::default(), Vec3::X, t).is_err());
    assert!(matches!(
        Boundary::circle(Point2::new(1e15, 0.), 1., t),
        Err(Error::PrecisionLoss)
    ));
}

#[test]
fn directions_normalize_without_overflow_or_underflow() {
    close(
        Vec3::new(1e308, 1e308, 0.).normalized().unwrap().length(),
        1.,
    );
    close(Vec3::new(1e-310, 0., 0.).normalized().unwrap().x, 1.);
}

#[test]
fn kernel_values_are_send_and_sync() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Solid>();
    assert_send_sync::<Profile>();
}

#[test]
fn every_wall_normal_points_out_of_material_including_holes() {
    let t = Tolerance::default();
    let outer = polygon(&[
        [0., 0.],
        [20., 0.],
        [20., 6.],
        [6., 6.],
        [6., 20.],
        [0., 20.],
    ]);
    let holes = vec![
        Boundary::circle(Point2::new(3., 12.), 1., t).unwrap(),
        polygon(&[[10., 2.], [14., 2.], [14., 4.], [10., 4.]]),
    ];
    let frame = Frame3::new(
        Point3::new(10., 20., 30.),
        Vec3::new(1., 2., 3.),
        Vec3::X,
        t,
    )
    .unwrap();
    let solid = Solid::extrude_with(
        OperationId::UNSPECIFIED,
        profile(outer, holes),
        frame,
        3.,
        -2.,
    )
    .map(|(s, _)| s)
    .unwrap();
    for (index, face) in solid.topology().faces().iter().enumerate().skip(2) {
        let origin = solid.topology().face_origin(FaceId::new(index));
        // A planar wall's four corners, or a cylinder wall's two ring-loop starts.
        let corners = solid
            .topology()
            .face_fins(FaceId::new(index))
            .iter()
            .flatten()
            .map(|c| c.pcurve.point(0.0))
            .collect::<Vec<_>>();
        let n = corners.len() as f64;
        let uv = Point2::new(
            corners.iter().map(|p| p.x).sum::<f64>() / n,
            corners.iter().map(|p| p.y).sum::<f64>() / n,
        );
        let on_face = face.surface.point(uv);
        let normal = face.normal(uv);
        assert_eq!(
            solid.classify(on_face + normal * 1e-4).unwrap(),
            Location::Outside,
            "{origin:?}"
        );
        assert_eq!(
            solid.classify(on_face + normal * -1e-4).unwrap(),
            Location::Inside,
            "{origin:?}"
        );
    }
}

#[test]
fn microscopic_moment_underflow_is_an_error() {
    let t = Tolerance::new(1e-160, 1e-12).unwrap();
    assert!(Boundary::circle(Point2::default(), 1e-150, t).is_err());
}
