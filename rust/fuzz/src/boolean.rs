//! Booleans of two prisms in one frame (S9a of REVIEW_NOTES.md): the split
//! target's line and arc profiles (rectangles, regular polygons, stadiums,
//! U shapes, squares with round or square holes) on dyadic sizes, the tool
//! offset by a dyadic vector in the axis-aligned frame (its origin moved) or
//! sharing the tilted frame's origin, its heights equal to the object's,
//! spanning them, overlapping them, disjoint from them, inside them or on
//! them, their profiles lines and arcs or the split target's splines
//! (S9a.2). Fuse, cut and common never panic and fail only as documented (a
//! cavity among several solids, S9a.2's, or arcs in frames with different
//! axes, S9c's; S9b turns, leans or tilts the tool's frame, S9c.2a stands it
//! on its side (perpendicular cylinders), S9d.1 makes it a sphere, a cap or
//! a zone, S9d.3a a cone or frustum, S9d.4a a whole torus; a result
//! thinner than the
//! resolution or touching itself; an undecided comparison); each result
//! validates as it is built and its history passes the independent check
//! (debug builds); when all three succeed their volumes agree,
//! `V(A ∪ B) = V(A) + V(B) - V(A ∩ B)` and `V(A - B) = V(A) - V(A ∩ B)`;
//! every result moves rigidly with its ids, and each of its vertices
//! classifies on its boundary; each operation's first result is an input
//! again (S9b.2) against a turned box, with the same identities.
use crate::analytic_intersections::Bytes;
use crate::split::{profile, spline_profile};
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Frame3, Point3, Solid, Tolerance, Vec3};

pub fn check_boolean(data: &[u8]) {
    let mut b = Bytes(data, 0);
    let (ka, kb, flags) = (b.next(), b.next(), b.next());
    let (s1, t1) = (
        1.0 + f64::from(b.next() % 16) / 4.0,
        0.5 + f64::from(b.next() % 16) / 8.0,
    );
    let (s2, t2) = (
        1.0 + f64::from(b.next() % 16) / 4.0,
        0.5 + f64::from(b.next() % 16) / 8.0,
    );
    let (dx, dy) = (
        f64::from(b.next() % 33) / 8.0 - 2.0,
        f64::from(b.next() % 33) / 8.0 - 2.0,
    );
    let h = 0.5 + f64::from(b.next() % 8) / 4.0;
    let (Some(mut pa), Some(mut pb)) = (profile(ka, s1, t1), profile(kb, s2, t2)) else {
        return;
    };
    let tolerance = Tolerance::default();
    let tilted = flags & 1 == 1;
    let (fa, fb) = if tilted {
        let Ok(f) = Frame3::new(
            Point3::new(1.0, -2.0, 0.5),
            Vec3::new(0.0, 3.0, 4.0),
            Vec3::new(1.0, 0.0, 0.0),
            tolerance,
        ) else {
            return;
        };
        (f, f)
    } else {
        let Ok(f) = Frame3::new(
            Point3::new(dx, dy, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            tolerance,
        ) else {
            return;
        };
        (Frame3::xy(), f)
    };
    // S9a.2 added heights inside the object's (pockets, cavities) and on
    // its top (touching stacks), chosen by a byte after the others.
    let (lo, hi) = match (b.next() % 3, (flags >> 1) % 4) {
        (1, _) => (h / 4.0, h * 0.75),
        (2, _) => (h, h + 1.0),
        (_, 0) => (0.0, h),
        (_, 1) => (-1.0, h + 1.0),
        (_, 2) => (h / 2.0, h * 1.5),
        _ => (h + 1.0, h + 2.0),
    };
    // S9b: the tool's frame turned about the axis, leaning or tilted
    // (frames with different axes), chosen by a byte after the others.
    // S9c.2a: or stood on its side (its axis along x, an exact frame),
    // by the byte's top bit.
    let pick = b.next();
    let fb = match (tilted, pick % 4) {
        (false, 0) if pick >= 128 => {
            let Ok(f) = Frame3::new(
                Point3::new(dx - h / 2.0, dy, h / 2.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                tolerance,
            ) else {
                return;
            };
            f
        }
        (false, turn @ 1..=3) => {
            let (normal, x) = match turn {
                1 => (Vec3::new(0.0, 0.0, 1.0), Vec3::new(3.0, 4.0, 0.0)),
                2 => (Vec3::new(3.0, 0.0, 4.0), Vec3::new(0.0, 1.0, 0.0)),
                _ => (Vec3::new(0.0, 3.0, 4.0), Vec3::new(1.0, 0.0, 0.0)),
            };
            let Ok(f) = Frame3::new(Point3::new(dx, dy, h / 2.0), normal, x, tolerance) else {
                return;
            };
            f
        }
        _ => fb,
    };
    // S9a.2: the split target's spline profiles for the object, the tool
    // or both, chosen by a byte after the others.
    let spline_byte = b.next();
    let splines = spline_byte % 4;
    if splines & 1 == 1 {
        let Some(p) = spline_profile(ka, s1, t1) else {
            return;
        };
        pa = p;
    }
    if splines & 2 == 2 {
        let Some(p) = spline_profile(kb, s2, t2) else {
            return;
        };
        pb = p;
    }
    let Ok((a, _)) = Solid::extrude_with(OperationId(1), pa, fa, 0.0, h) else {
        return;
    };
    // S9d.1: a sphere, a cap or a zone for the tool in its frame, by the
    // spline byte's two top bits.
    let tool = if spline_byte >= 192 {
        let half = std::f64::consts::FRAC_PI_2;
        let (low, high) = [(-half, half), (-half, 0.0), (-0.5, 0.75), (0.25, half)]
            [usize::from((spline_byte >> 2) % 4)];
        let Ok((tool, _)) = Solid::sphere_with(OperationId(2), fb, 0.75 * s2, low, high, tolerance)
        else {
            return;
        };
        tool
    } else if (144..160).contains(&spline_byte) {
        // S9d.4a: a whole torus about the tool's frame, its radii by the
        // byte's low bits.
        let big = 0.75 * s2;
        let small = big * [0.25, 0.375, 0.5, 0.625][usize::from(spline_byte % 4)];
        let turn = std::f64::consts::TAU;
        let Ok((tool, _)) =
            Solid::torus_with(OperationId(2), fb, big, small, 0.0, turn, turn, tolerance)
        else {
            return;
        };
        tool
    } else if spline_byte >= 160 {
        // S9d.3a: a cone or frustum over the tool's heights, its radii by
        // the byte's next bits.
        let r = 0.75 * s2;
        let (bottom, top) =
            [(r, 0.0), (0.0, r), (r, r / 2.0), (r / 2.0, r)][usize::from((spline_byte >> 2) % 4)];
        let Ok(base) = Frame3::new(
            fb.point(rusty_occt::Point2::new(0.0, 0.0), lo),
            fb.normal(),
            fb.x(),
            tolerance,
        ) else {
            return;
        };
        let Ok((tool, _)) = Solid::cone_with(OperationId(2), base, bottom, top, hi - lo, tolerance)
        else {
            return;
        };
        tool
    } else {
        let Ok((tool, _)) = Solid::extrude_with(OperationId(2), pb, fb, lo, hi) else {
            return;
        };
        tool
    };
    let run = |r: Result<(Vec<Solid>, rusty_occt::history::History), Error>| -> Option<Vec<Solid>> {
        match r {
            Ok((out, _)) => Some(out),
            Err(Error::Degenerate(_) | Error::ComputationLimit(_)) => None,
            // A cavity among several solids (S9a.2's), or arcs in frames
            // with different axes (S9c's).
            // Splines along one curve of different forms, or a spline span
            // along a line (S9a.2's).
            Err(Error::OutOfDomain(m))
                if m.contains("cavity")
                    || m.contains("S9c")
                    || m.contains("S9d")
                    || m.contains("different forms")
                    || m.contains("along the plane") =>
            {
                None
            }
            Err(e) => panic!("unexpected error {e}"),
        }
    };
    let fused = run(a.fuse(OperationId(3), &tool));
    let cut = run(a.cut(OperationId(4), &tool));
    let common = run(a.common(OperationId(5), &tool));
    let volume = |out: &[Solid]| -> f64 { out.iter().map(|s| s.mass_properties().volume).sum() };
    let (va, vb) = (a.mass_properties().volume, tool.mass_properties().volume);
    let near = |x: f64, y: f64| (x - y).abs() <= 1e-9 * x.abs().max(y.abs()).max(1.0);
    if let (Some(f), Some(m)) = (&fused, &common) {
        assert!(
            near(volume(f), va + vb - volume(m)),
            "fuse {} for {va} + {vb} - {}",
            volume(f),
            volume(m)
        );
    }
    if let (Some(c), Some(m)) = (&cut, &common) {
        assert!(
            near(volume(c), va - volume(m)),
            "cut {} for {va} - {}",
            volume(c),
            volume(m)
        );
    }
    // S9b.2: an operation's first result is an input again, against a
    // turned box about the object's origin (a fresh operation's ids).
    let turned = Frame3::new(
        fa.point(rusty_occt::Point2::new(0.5, 0.25), h / 3.0),
        fa.normal(),
        fa.x() * 3.0 + fa.y() * 4.0,
        tolerance,
    )
    .ok()
    .and_then(|f| {
        let square = rusty_occt::Boundary::polygon(
            vec![
                rusty_occt::Point2::new(-1.0, -1.0),
                rusty_occt::Point2::new(1.0, -1.0),
                rusty_occt::Point2::new(1.0, 1.0),
                rusty_occt::Point2::new(-1.0, 1.0),
            ],
            tolerance,
        )
        .ok()?;
        let profile = rusty_occt::Profile::new(square, vec![], tolerance).ok()?;
        Solid::extrude_with(OperationId(7), profile, f, 0.0, h).ok()
    });
    if let Some((box_, _)) = turned {
        // One operation's first result, chosen by a byte after the others,
        // of at most 12 faces, cut by the box and in common with it: exact
        // fragments of larger stored models took up to 165 s an input under
        // ASan (fuzz/regressions/README.md); the kernel's tests take them.
        let chosen = [&fused, &cut, &common][usize::from(b.next() % 3)];
        let small = |s: &&Solid| s.topology().faces().len() <= 12;
        if let Some(first) = chosen.as_ref().and_then(|out| out.first()).filter(small) {
            let (c, m) = (
                run(first.cut(OperationId(9), &box_)),
                run(first.common(OperationId(10), &box_)),
            );
            if let (Some(c), Some(m)) = (&c, &m) {
                let v1 = first.mass_properties().volume;
                assert!(near(volume(c), v1 - volume(m)), "chained cut");
            }
        }
    }
    let motion = rusty_occt::RigidTransform::rotation(
        Point3::new(0.5, -1.0, 2.0),
        Vec3::new(1.0, 2.0, 2.0),
        0.5,
    )
    .expect("a rotation");
    for out in [&fused, &cut, &common].into_iter().flatten() {
        for piece in out {
            let (moved, _) = piece
                .transform_with(OperationId(6), motion)
                .expect("a result moves rigidly");
            let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
            assert_eq!(ids(piece), ids(&moved), "a moved result keeps its ids");
            for v in piece.topology().vertices() {
                assert_eq!(
                    piece.classify(v.position).expect("a vertex classifies"),
                    rusty_occt::Location::Boundary,
                    "a result's vertex on its boundary"
                );
            }
        }
    }
}
