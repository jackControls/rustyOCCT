//! S9d.2's near miss at a vertex (REVIEW_NOTES.md, "A vertex within the
//! resolution of a face"; `curved/near.rs`'s `vertex_faces` and
//! `polyhedra.rs`'s `vertex_near_misses`): an input vertex whose boundary
//! leaves it strictly on one side of the other input's plane, cylinder or
//! cone surface, within the resolution of it, crossing it or missing it
//! with the gap outside either input, is `Degenerate`. A vertex on the face
//! exactly, one far through it, and one whose input runs along the surface
//! there (U1's box on a rod's tangent plane) evaluate as before.
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Point2, Point3, Profile, RigidTransform, Solid, Tolerance, Vec3,
};

const NEAR_VERTEX: &str = "a vertex within the resolution of a face (S9d.2)";

fn tol() -> Tolerance {
    Tolerance::default()
}

fn frame(o: [f64; 3], n: [f64; 3], x: [f64; 3]) -> Frame3 {
    Frame3::new(
        Point3::new(o[0], o[1], o[2]),
        Vec3::new(n[0], n[1], n[2]),
        Vec3::new(x[0], x[1], x[2]),
        tol(),
    )
    .unwrap()
}

/// A level frame, one turned off the world's z (the boolean target's) and
/// one turned about it by 30 degrees.
fn frames() -> [(&'static str, Frame3); 3] {
    [
        ("level", frame([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])),
        (
            "turned",
            frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
        ),
        (
            "about z",
            frame([0.0; 3], [0.0, 0.0, 1.0], [3f64.sqrt() / 2.0, 0.5, 0.0]),
        ),
    ]
}

fn add(a: [f64; 3], b: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] + k * b[0], a[1] + k * b[1], a[2] + k * b[2]]
}

fn neg(a: [f64; 3]) -> [f64; 3] {
    a.map(|x| -x)
}

fn square(lo: f64, hi: f64) -> Boundary {
    Boundary::polygon(
        [(lo, lo), (hi, lo), (hi, hi), (lo, hi)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec(),
        tol(),
    )
    .unwrap()
}

fn prism(f: Frame3, outer: Boundary, lo: f64, hi: f64, id: u64) -> Solid {
    Solid::extrude_with(
        OperationId(id),
        Profile::new(outer, vec![], tol()).unwrap(),
        f,
        lo,
        hi,
    )
    .unwrap()
    .0
}

/// A cube of side `s` on its corner `o`, each of its three edges from it
/// leaving along `t` (at `1 / sqrt(2)`, `1 / sqrt(3)` and `1 / sqrt(6)` of
/// its length), `(t, p1, p2)` a right-handed orthonormal basis: its frame's
/// normal `t + p2`, its `x` `t - p1 - p2` and so its `y` `t + 2 p1 - p2`.
fn cube(o: [f64; 3], t: [f64; 3], p1: [f64; 3], p2: [f64; 3], s: f64) -> Solid {
    let n = add(t, p2, 1.0);
    let x = add(add(t, p1, -1.0), p2, -1.0);
    prism(frame(o, n, x), square(0.0, s), 0.0, s, 9)
}

/// Every operation, either way round, refused with `why`.
fn refused(what: &str, a: &Solid, b: &Solid, why: &str) {
    for (op, out) in [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ] {
        assert!(
            matches!(&out, Err(Error::Degenerate(m)) if *m == why),
            "{what} {op}: {:?}",
            out.map(|x| x.0.len())
        );
    }
}

fn volume(s: &[Solid]) -> f64 {
    s.iter().map(|x| x.mass_properties().volume).sum()
}

/// Fuse, cut and common evaluate, the fuse `fused` solids, in the pair
/// identities; the common's volume.
fn evaluates(what: &str, a: &Solid, b: &Solid, fused: usize) -> f64 {
    let f = a.fuse(OperationId(3), b).unwrap().0;
    let c = a.cut(OperationId(4), b).unwrap().0;
    let back = b.cut(OperationId(4), a).unwrap().0;
    let m = a.common(OperationId(5), b).unwrap().0;
    assert_eq!(f.len(), fused, "{what}");
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let close = |x: f64, y: f64| (x - y).abs() <= 1e-9 * (va + vb);
    let common = volume(&m);
    assert!(close(volume(&f), va + vb - common), "{what} fuse");
    assert!(close(volume(&c), va - common), "{what} cut");
    assert!(close(volume(&back), vb - common), "{what} cut back");
    common
}

/// The gaps within the resolution (`1e-7`), either side.
const WITHIN: [f64; 6] = [1e-12, -1e-12, 1e-11, -1e-11, 1e-9, -1e-9];

/// A cube's corner within the resolution of a rod's wall (radius 1) is
/// `Degenerate`: outside it, at `1 + d` from the axis, the fuse was two
/// solids and the common empty; inside it, its corner at `1 - d` toward the
/// wall, the cut a cavity behind a wall thinner than the resolution; a
/// corner through the wall that deep was refused as a piece or a result
/// thinner than the resolution. So is a corner by a frustum's wall (radius
/// 2 to 1 over a height of 3, the corner on the normal through its point at
/// height 1.5). Each in a level frame, one turned off the world's z and one
/// turned about it, the corner placed in binary64 (its gap's rounding far
/// below `1e-12`); at `1e-6`, beyond the resolution, every operation
/// evaluates.
#[test]
fn a_corner_within_the_resolution_of_a_cylinder_or_cone_is_degenerate() {
    for (name, f) in frames() {
        let (u, v, w) = (f.x().to_array(), f.y().to_array(), f.normal().to_array());
        let rod = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 3.0, tol())
            .unwrap()
            .0;
        let at = |r: f64| f.point(Point2::new(r, 0.0), 0.0).to_array();
        let outside = |d: f64| cube(at(1.0 + d), u, v, w, 1.0);
        let inside = |d: f64| cube(at(1.0 - d), neg(u), neg(v), w, 0.5);
        for d in WITHIN {
            refused(&format!("{name} rod {d:e}"), &rod, &outside(d), NEAR_VERTEX);
            refused(
                &format!("{name} in the rod {d:e}"),
                &rod,
                &inside(d),
                NEAR_VERTEX,
            );
        }
        assert_eq!(
            evaluates(&format!("{name} rod"), &rod, &outside(1e-6), 2),
            0.0
        );
        let common = evaluates(&format!("{name} in the rod"), &rod, &inside(1e-6), 1);
        assert!((common - 0.125).abs() <= 1e-12, "{name} in the rod");
        // A frustum's wall: its outward normal at local `(1.5, 0, 1.5)` is
        // `(3, 0, 1) / sqrt(10)`, its generatrix `(-1, 0, 3) / sqrt(10)`.
        let cone = Solid::cone_with(OperationId(1), f, 2.0, 1.0, 3.0, tol())
            .unwrap()
            .0;
        let s10 = 10f64.sqrt();
        let normal = add(u.map(|x| 3.0 * x / s10), w, 1.0 / s10);
        let along = add(u.map(|x| -x / s10), w, 3.0 / s10);
        let foot = f.point(Point2::new(1.5, 0.0), 1.5).to_array();
        let by = |d: f64| cube(add(foot, normal, d), normal, v, along, 1.0);
        for d in WITHIN {
            refused(&format!("{name} frustum {d:e}"), &cone, &by(d), NEAR_VERTEX);
        }
        assert_eq!(
            evaluates(&format!("{name} frustum"), &cone, &by(1e-6), 2),
            0.0
        );
    }
}

/// A box's corner within the resolution of another box's face is
/// `Degenerate` (the polyhedral engine, S9b): a slab's top and a cube on
/// its corner `d` above it were fused into two solids, the corner `d` into
/// it refused as a face thinner than the resolution, and in the turned
/// frame a corner placed on the top but for its rounding fused into two
/// solids; in the curved engine a cube's corner by a rod's cap alike. Each
/// in the three frames; at `1e-6` every operation evaluates.
#[test]
fn a_corner_within_the_resolution_of_a_plane_face_is_degenerate() {
    for (name, f) in frames() {
        let (u, v, w) = (f.x().to_array(), f.y().to_array(), f.normal().to_array());
        let slab = prism(f, square(-2.0, 2.0), -2.0, 0.0, 1);
        let above = |d: f64| cube(f.point(Point2::new(0.1, 0.2), d).to_array(), w, u, v, 1.0);
        for d in WITHIN {
            refused(&format!("{name} slab {d:e}"), &slab, &above(d), NEAR_VERTEX);
        }
        if name == "turned" {
            // On the top but for the turned frame's rounding.
            refused("turned slab 0", &slab, &above(0.0), NEAR_VERTEX);
        }
        assert_eq!(
            evaluates(&format!("{name} slab"), &slab, &above(1e-6), 2),
            0.0
        );
        let rod = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 0.0, tol())
            .unwrap()
            .0;
        for d in WITHIN {
            refused(
                &format!("{name} rod's cap {d:e}"),
                &rod,
                &above(d),
                NEAR_VERTEX,
            );
        }
        assert_eq!(
            evaluates(&format!("{name} rod's cap"), &rod, &above(1e-6), 2),
            0.0
        );
    }
}

/// What stays no contact of a vertex's own evaluates as before the rule:
/// a corner on a face exactly (an incidence: the slab's cut and common, its
/// fuse two solids touching at a point refused as before; on the rod's
/// wall, refused as before by the incidences' own rule, the vertex on the
/// other's face retried at every seam), a corner `0.1` through the face, and
/// a box whose face lies on a rod's tangent plane with its corner on the
/// ruling within rounding (U1's configuration, turned by 30 degrees: its
/// edges and face running along the surface, the faces' own tangency).
#[test]
fn a_corner_on_through_or_along_a_face_evaluates_as_before() {
    for (name, f) in frames() {
        let (u, v, w) = (f.x().to_array(), f.y().to_array(), f.normal().to_array());
        let slab = prism(f, square(-2.0, 2.0), -2.0, 0.0, 1);
        let above = |d: f64| cube(f.point(Point2::new(0.1, 0.2), d).to_array(), w, u, v, 1.0);
        // On the top exactly where the frame keeps its height (the turned
        // frame's corner lies off it by rounding: a near miss).
        if name != "turned" {
            let touching = above(0.0);
            let cut = volume(&slab.cut(OperationId(4), &touching).unwrap().0);
            assert!((cut - 32.0).abs() <= 1e-12, "{name} touching cut {cut}");
            assert!(slab.common(OperationId(5), &touching).unwrap().0.is_empty());
            assert!(matches!(
                slab.fuse(OperationId(3), &touching),
                Err(Error::Degenerate("two solids touching at a point"))
            ));
        }
        evaluates(&format!("{name} slab through"), &slab, &above(-0.1), 1);
        let rod = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 3.0, tol())
            .unwrap()
            .0;
        let at = |r: f64| f.point(Point2::new(r, 0.0), 0.0).to_array();
        evaluates(
            &format!("{name} rod through"),
            &rod,
            &cube(at(0.9), u, v, w, 1.0),
            1,
        );
    }
    // The level rod's wall through the corner exactly.
    let level = frames()[0].1;
    let rod = Solid::cylinder_with(OperationId(1), level, 1.0, -3.0, 3.0, tol())
        .unwrap()
        .0;
    let on = cube(
        [0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        1.0,
    );
    refused(
        "on the rod's wall",
        &rod,
        &on,
        "a meeting at every seam tried",
    );
    // U1's configuration: a rod of radius 1 and height 2, and boxes on its
    // tangent plane `y = -1` turned by 30 degrees about its axis.
    let rod = Solid::cylinder_with(OperationId(1), level, 1.0, 0.0, 2.0, tol())
        .unwrap()
        .0;
    let block = |o: [f64; 3], s: [f64; 3]| {
        let b = Solid::box_at(
            Point3::new(o[0], o[1], o[2]),
            Vec3::new(s[0], s[1], s[2]),
            tol(),
        )
        .unwrap();
        let t = RigidTransform::rotation(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            30f64.to_radians(),
        )
        .unwrap();
        b.transform_with(OperationId(8), t).unwrap().0
    };
    let tau = std::f64::consts::TAU;
    // Outside, its corners `(0, -1, 0.5)` and `(0, -1, 1.5)` on the ruling.
    let ruled = block([0.0, -2.0, 0.5], [1.0, 1.0, 1.0]);
    assert_eq!(evaluates("on the ruling", &rod, &ruled, 2), 0.0);
    let cut = volume(&rod.cut(OperationId(4), &ruled).unwrap().0);
    assert!((cut - tau).abs() <= 1e-12, "on the ruling cut {cut}");
    // U1's box, `[-r, r] x [-1, r] x [0, 2]` with `r = sqrt(3) / 2`.
    let r = 3f64.sqrt() / 2.0;
    evaluates(
        "U1",
        &rod,
        &block([-r, -1.0, 0.0], [2.0 * r, 1.0 + r, 2.0]),
        1,
    );
}
