//! S9d.2's near misses between the inputs' edges and vertices
//! (REVIEW_NOTES.md, "Edges and vertices within the resolution of each
//! other"; `curved/near.rs` and `polyhedra.rs`): a vertex within the
//! resolution of the other input's vertex or edge, two edges passing within
//! it of each other, and a cone's apex within it of a face, each input's
//! boundary leaving its point strictly away from the other's (missing, the
//! gap outside either input) or heading strictly toward it (crossing), are
//! `Degenerate`; so is a ball off a rod's rim at the rim's seam, where the
//! edge rules took no point. Exact crossings and incidences, configurations
//! well apart, U1's and Y5's box on a rod's tangent plane and the vertex
//! rule's controls evaluate as before.
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Point2, Point3, Profile, RigidTransform, Solid, Tolerance, Vec3,
};

const NEAR_VV: &str = "a vertex within the resolution of a vertex of the other input (S9d.2)";
const NEAR_VE: &str = "a vertex within the resolution of an edge of the other input (S9d.2)";
const NEAR_EE: &str = "two edges within the resolution of each other (S9d.2)";
const NEAR_APEX: &str = "a cone's apex within the resolution of a face (S9d.3)";
const NEAR_SPHERE_EDGE: &str = "a sphere within the resolution of tangency to an edge (S9d.1)";

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

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    a.map(|x| x / l)
}

/// A frame's point and direction of local coordinates.
fn at(f: &Frame3, l: [f64; 3]) -> [f64; 3] {
    f.point(Point2::new(l[0], l[1]), l[2]).to_array()
}

fn dir(f: &Frame3, l: [f64; 3]) -> [f64; 3] {
    let (x, y, n) = (f.x().to_array(), f.y().to_array(), f.normal().to_array());
    add(add(x.map(|c| c * l[0]), y, l[1]), n, l[2])
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
/// its length), `(t, p1, p2)` a right-handed orthonormal basis.
fn cube(o: [f64; 3], t: [f64; 3], p1: [f64; 3], p2: [f64; 3], s: f64) -> Solid {
    let n = add(t, p2, 1.0);
    let x = add(add(t, p1, -1.0), p2, -1.0);
    prism(frame(o, n, x), square(0.0, s), 0.0, s, 9)
}

/// A cube on its corner `o` leaving along the unit `t`, `p1` a unit normal
/// to it, turned about `t` by 0.4.
fn corner_cube(o: [f64; 3], t: [f64; 3], p1: [f64; 3]) -> Solid {
    // Turned about `t` by 0.4 (no face along the other's planes or axis).
    let q = cross(t, p1);
    let (c, s) = (0.4f64.cos(), 0.4f64.sin());
    let p1 = add(p1.map(|x| c * x), q, s);
    cube(o, t, p1, cross(t, p1), 1.0)
}

/// A unit block standing on its edge through `o` along the unit `e` (from
/// `-1` to `1` along it), its faces leaving the edge along `(t + p) /
/// sqrt(2)` and `(t - p) / sqrt(2)`, `(e, t, p)` with `e x t = -p`.
fn edge_block(o: [f64; 3], e: [f64; 3], t: [f64; 3], p: [f64; 3]) -> Solid {
    let x = unit(add(t, p, 1.0));
    prism(frame(o, e, x), square(0.0, 1.0), -1.0, 1.0, 9)
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

/// Every operation refused as degenerate (a crossing within the resolution:
/// the new rules' reason or an existing rule's), the reasons returned.
fn degenerate(what: &str, a: &Solid, b: &Solid) -> Vec<&'static str> {
    [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ]
    .into_iter()
    .map(|(op, out)| match out {
        Err(Error::Degenerate(m)) => m,
        other => panic!("{what} {op}: {:?}", other.map(|x| x.0.len())),
    })
    .collect()
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

/// The misses within the resolution (`1e-7`).
const MISSES: [f64; 3] = [1e-12, 1e-11, 1e-9];

/// A cube's corner within the resolution of a slab's top edge or corner,
/// its foot on the top's plane just past the face (diagonally off the
/// edge, and off the corner), and of a rod's rim, is `Degenerate`: the fuse
/// was two solids and the common empty where the inputs lie no farther
/// apart than the resolution. Each in the three frames, the corner placed
/// in binary64 (its gap's rounding far below `1e-12`); a corner through the
/// edge that deep is refused too; at `1e-6` every operation evaluates.
#[test]
fn a_corner_within_the_resolution_of_an_edge_or_corner_is_degenerate() {
    for (name, f) in frames() {
        let slab = prism(f, square(-2.0, 2.0), -2.0, 0.0, 1);
        // Off the top edge along `v` at `u = 2`.
        let e = at(&f, [2.0, 0.3, 0.0]);
        let t = unit(dir(&f, [1.0, 0.0, 1.0]));
        let p1 = dir(&f, [0.0, 1.0, 0.0]);
        let by_edge = |d: f64| corner_cube(add(e, t, d), t, p1);
        for d in MISSES {
            refused(&format!("{name} edge {d:e}"), &slab, &by_edge(d), NEAR_VE);
            degenerate(&format!("{name} edge {:e}", -d), &slab, &by_edge(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} edge"), &slab, &by_edge(1e-6), 2),
            0.0
        );
        // Off the top corner `(2, 2, 0)`.
        let c = at(&f, [2.0, 2.0, 0.0]);
        let t = unit(dir(&f, [1.0, 1.0, 1.0]));
        let p1 = unit(dir(&f, [1.0, -1.0, 0.0]));
        let by_corner = |d: f64| corner_cube(add(c, t, d), t, p1);
        for d in MISSES {
            refused(
                &format!("{name} corner {d:e}"),
                &slab,
                &by_corner(d),
                NEAR_VV,
            );
            degenerate(&format!("{name} corner {:e}", -d), &slab, &by_corner(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} corner"), &slab, &by_corner(1e-6), 2),
            0.0
        );
        // Off a rod's top rim at the angle 0.7, diagonally.
        let rod = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 0.0, tol())
            .unwrap()
            .0;
        let (cs, sn) = (0.7f64.cos(), 0.7f64.sin());
        let r = at(&f, [cs, sn, 0.0]);
        let t = unit(dir(&f, [cs, sn, 1.0]));
        let p1 = dir(&f, [-sn, cs, 0.0]);
        let by_rim = |d: f64| corner_cube(add(r, t, d), t, p1);
        for d in MISSES {
            refused(&format!("{name} rim {d:e}"), &rod, &by_rim(d), NEAR_VE);
            degenerate(&format!("{name} rim {:e}", -d), &rod, &by_rim(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} rim"), &rod, &by_rim(1e-6), 2),
            0.0
        );
    }
}

/// Two edges passing within the resolution of each other (skew) are
/// `Degenerate`: a block standing on its edge across a slab's top edge, its
/// edge `d` off the slab's along their common normal (the diagonal between
/// the top and the wall), its faces leaving it away from the slab, was
/// fused with it into two solids; so was one across a rod's rim (a line and
/// a circle, a conic in the turned frames). Crossing that deep is refused
/// too; at `1e-6` every operation evaluates.
#[test]
fn edges_within_the_resolution_of_each_other_are_degenerate() {
    for (name, f) in frames() {
        let slab = prism(f, square(-2.0, 2.0), -2.0, 0.0, 1);
        let e = at(&f, [2.0, 0.3, 0.0]);
        let t = unit(dir(&f, [1.0, 0.0, 1.0]));
        let along = unit(dir(&f, [1.0, 0.0, -1.0]));
        let j = dir(&f, [0.0, 1.0, 0.0]);
        let across = |d: f64| edge_block(add(e, t, d), along, t, j);
        for d in MISSES {
            refused(&format!("{name} edges {d:e}"), &slab, &across(d), NEAR_EE);
            degenerate(&format!("{name} edges {:e}", -d), &slab, &across(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} edges"), &slab, &across(1e-6), 2),
            0.0
        );
        let rod = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 0.0, tol())
            .unwrap()
            .0;
        let (cs, sn) = (0.7f64.cos(), 0.7f64.sin());
        let r = at(&f, [cs, sn, 0.0]);
        let t = unit(dir(&f, [cs, sn, 1.0]));
        let along = unit(dir(&f, [cs, sn, -1.0]));
        let tangent = dir(&f, [-sn, cs, 0.0]);
        let over = |d: f64| edge_block(add(r, t, d), along, t, tangent);
        for d in MISSES {
            refused(&format!("{name} rim {d:e}"), &rod, &over(d), NEAR_EE);
            degenerate(&format!("{name} rim {:e}", -d), &rod, &over(-d));
        }
        assert_eq!(evaluates(&format!("{name} rim"), &rod, &over(1e-6), 2), 0.0);
    }
}

/// A cone's apex within the resolution of a face is `Degenerate`: a cone
/// standing `d` under a slab's bottom was fused with it into two solids;
/// its apex that far through the slab is refused too; one lying on its side
/// with its apex `d` off a rod's wall or through it is refused as before;
/// at `1e-6` every operation evaluates.
#[test]
fn an_apex_within_the_resolution_of_a_face_is_degenerate() {
    for (name, f) in frames() {
        let cone = Solid::cone_with(OperationId(1), f, 1.0, 0.0, 1.0, tol())
            .unwrap()
            .0;
        let slab = |d: f64| prism(f, square(-2.0, 2.0), 1.0 + d, 2.0, 2);
        for d in MISSES {
            refused(&format!("{name} slab {d:e}"), &cone, &slab(d), NEAR_APEX);
            refused(
                &format!("{name} slab {:e}", -d),
                &cone,
                &slab(-d),
                NEAR_APEX,
            );
        }
        assert_eq!(
            evaluates(&format!("{name} slab"), &cone, &slab(1e-6), 2),
            0.0
        );
        let rod = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 3.0, tol())
            .unwrap()
            .0;
        let lying = |d: f64| {
            let g = frame(
                at(&f, [2.0 + d, 0.0, 0.5]),
                dir(&f, [-1.0, 0.0, 0.0]),
                dir(&f, [0.0, 0.0, 1.0]),
            );
            Solid::cone_with(OperationId(2), g, 0.5, 0.0, 1.0, tol())
                .unwrap()
                .0
        };
        // Refused as before by S9c's node rule (the cone's meeting with the
        // rod's cylinder within the resolution of a node), which comes first.
        for d in MISSES {
            degenerate(&format!("{name} rod {d:e}"), &rod, &lying(d));
            degenerate(&format!("{name} rod {:e}", -d), &rod, &lying(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} rod"), &rod, &lying(1e-6), 2),
            0.0
        );
    }
}

/// A rod's rim is two arcs meeting at its seam (`(45, 28) / 53` for a
/// Boolean's first input), where the edge rules took no point and the
/// vertex rules none either (a seam's vertex): a ball diagonally off the
/// rim whose nearest point is the seam, `d` from it, was fused with the rod
/// into two solids, where the cut of the ball by the rod (its seam
/// elsewhere) was refused; and a cube's corner there likewise. Both are
/// `Degenerate` now in every operation; at `1e-6` every operation
/// evaluates.
#[test]
fn a_near_miss_at_a_rims_seam_is_degenerate() {
    let level = frames()[0].1;
    let rod = Solid::cylinder_with(OperationId(1), level, 1.0, -3.0, 0.0, tol())
        .unwrap()
        .0;
    // The ball's centre `(45, 28, 21) / 32`, `(21 / 32) sqrt(2)` from the
    // seam `(45, 28, 0) / 53`, its point nearest the rim.
    let reach = 21.0 / 32.0 * 2f64.sqrt();
    let ball = |r: f64| {
        Solid::sphere_with(
            OperationId(9),
            frame(
                [45.0 / 32.0, 28.0 / 32.0, 21.0 / 32.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
            ),
            r,
            -r,
            r,
            tol(),
        )
        .unwrap()
        .0
    };
    for d in MISSES {
        refused(
            &format!("ball {d:e}"),
            &rod,
            &ball(reach - d),
            NEAR_SPHERE_EDGE,
        );
        degenerate(&format!("ball {:e}", -d), &rod, &ball(reach + d));
    }
    assert_eq!(evaluates("ball", &rod, &ball(reach - 1e-6), 2), 0.0);
    // A cube's corner at `(45, 28, 53) u - (0, 0, 1)`, `u` a short dyadic
    // just past `1 / 53`: its point nearest the rim the seam exactly, at
    // `sqrt(2) (53 u - 1)`.
    let t = unit([45.0 / 53.0, 28.0 / 53.0, 1.0]);
    let p1 = [-28.0 / 53.0, 45.0 / 53.0, 0.0];
    let corner = |d: f64| {
        let k = ((1.0 + d / 2f64.sqrt()) * 2f64.powi(50) / 53.0).ceil();
        let u = k / 2f64.powi(50);
        let o = [45.0 * u, 28.0 * u, 53.0 * u - 1.0];
        corner_cube(o, t, p1)
    };
    for d in MISSES {
        refused(&format!("seam corner {d:e}"), &rod, &corner(d), NEAR_VE);
    }
    assert_eq!(evaluates("seam corner", &rod, &corner(1e-6), 2), 0.0);
}

/// What is no near miss of the edges' and vertices' own evaluates as
/// before: edges and corners `0.1` apart, U1's and Y5's boxes on a rod's
/// tangent plane turned by 30 degrees (the bottom edge tangent to the rod's
/// rim but for rounding, the box's face along it), a box on a rod's
/// tangent plane with its corners on the ruling, a corner through a slab's
/// top, and an exact crossing of two edges and a corner on an edge exactly
/// (incidences, refused by their own rules as before).
#[test]
fn edges_and_vertices_apart_on_or_along_evaluate_as_before() {
    for (name, f) in frames() {
        let slab = prism(f, square(-2.0, 2.0), -2.0, 0.0, 1);
        let e = at(&f, [2.0, 0.3, 0.0]);
        let t = unit(dir(&f, [1.0, 0.0, 1.0]));
        let along = unit(dir(&f, [1.0, 0.0, -1.0]));
        let j = dir(&f, [0.0, 1.0, 0.0]);
        assert_eq!(
            evaluates(
                &format!("{name} edges apart"),
                &slab,
                &edge_block(add(e, t, 0.1), along, t, j),
                2
            ),
            0.0
        );
        assert_eq!(
            evaluates(
                &format!("{name} corner apart"),
                &slab,
                &corner_cube(add(e, t, 0.1), t, j),
                2
            ),
            0.0
        );
        evaluates(
            &format!("{name} corner through"),
            &slab,
            &corner_cube(add(e, t, -0.1), t, j),
            1,
        );
    }
    // Exactly on: the level frame's edge point `(2, 0.3, 0)`.
    let level = frames()[0].1;
    let slab = prism(level, square(-2.0, 2.0), -2.0, 0.0, 1);
    let t = unit([1.0, 0.0, 1.0]);
    for (what, other) in [
        (
            "edges crossing",
            edge_block([2.0, 0.3, 0.0], unit([1.0, 0.0, -1.0]), t, [0.0, 1.0, 0.0]),
        ),
        (
            "corner on the edge",
            corner_cube([2.0, 0.3, 0.0], t, [0.0, 1.0, 0.0]),
        ),
    ] {
        for (op, out) in [
            ("fuse", slab.fuse(OperationId(3), &other)),
            ("cut", slab.cut(OperationId(4), &other)),
            ("common", slab.common(OperationId(5), &other)),
        ] {
            if let Err(Error::Degenerate(m)) = &out {
                assert!(![NEAR_VV, NEAR_VE, NEAR_EE].contains(m), "{what} {op}: {m}");
            }
        }
    }
    // U1 and Y5: a rod of radius 1 and height 2, and boxes on its tangent
    // plane `y = -1` turned by 30 degrees about its axis.
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
    let r = 3f64.sqrt() / 2.0;
    for (what, h) in [("U1", 2.0), ("Y5", 1.0)] {
        evaluates(
            what,
            &rod,
            &block([-r, -1.0, 0.0], [2.0 * r, 1.0 + r, h]),
            1,
        );
    }
    // Outside, its corners `(0, -1, 0.5)` and `(0, -1, 1.5)` on the ruling.
    let ruled = block([0.0, -2.0, 0.5], [1.0, 1.0, 1.0]);
    assert_eq!(evaluates("on the ruling", &rod, &ruled, 2), 0.0);
}
