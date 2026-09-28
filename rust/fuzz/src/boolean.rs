//! Booleans of two prisms in one frame (S9a of REVIEW_NOTES.md): the split
//! target's line and arc profiles (rectangles, regular polygons, stadiums,
//! U shapes, squares with round or square holes) on dyadic sizes, the tool
//! offset by a dyadic vector in the axis-aligned frame (its origin moved) or
//! sharing the tilted frame's origin, its heights equal to the object's,
//! spanning them, overlapping them or disjoint from them. Fuse, cut and
//! common never panic and fail only as documented (a stack of slabs,
//! S9a.2's; a result thinner than the resolution or touching itself; an
//! undecided comparison); each result validates as it is built and its
//! history passes the independent check (debug builds); when all three
//! succeed their volumes agree, `V(A ∪ B) = V(A) + V(B) - V(A ∩ B)` and
//! `V(A - B) = V(A) - V(A ∩ B)`; every result moves rigidly with its ids.
use crate::analytic_intersections::Bytes;
use crate::split::profile;
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
    let (Some(pa), Some(pb)) = (profile(ka, s1, t1), profile(kb, s2, t2)) else {
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
    let (lo, hi) = match (flags >> 1) % 4 {
        0 => (0.0, h),
        1 => (-1.0, h + 1.0),
        2 => (h / 2.0, h * 1.5),
        _ => (h + 1.0, h + 2.0),
    };
    let Ok((a, _)) = Solid::extrude_with(OperationId(1), pa, fa, 0.0, h) else {
        return;
    };
    let Ok((tool, _)) = Solid::extrude_with(OperationId(2), pb, fb, lo, hi) else {
        return;
    };
    let run = |r: Result<(Vec<Solid>, rusty_occt::history::History), Error>| -> Option<Vec<Solid>> {
        match r {
            Ok((out, _)) => Some(out),
            Err(Error::OutOfDomain(_) | Error::Degenerate(_) | Error::ComputationLimit(_)) => None,
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
        }
    }
}
