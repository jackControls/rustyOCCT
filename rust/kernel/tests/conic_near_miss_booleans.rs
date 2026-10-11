//! S9d.2's near misses between two conic or circle edges (REVIEW_NOTES.md,
//! "Two conic edges within the resolution of each other"; `curved/near.rs`):
//! two rims passing within the resolution of each other at a point, each
//! input's boundary leaving its point strictly away from the other's
//! (missing, the gap outside either input) or heading strictly toward it
//! (crossing), are `Degenerate` as "two edges within the resolution of each
//! other (S9d.2)": a rod's rim by another rod's, by a cone's and by a
//! block's arc edge (a rounded corner), in a level frame, one turned off the
//! world's z and one turned about it. Exact rim crossings and incidences,
//! rims well apart, coaxial rims on one plane, U1's box on a rod's tangent
//! plane and the residue's controls evaluate as before.
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Point2, Point3, Profile, RigidTransform, Segment, Solid, Tolerance,
    Vec3,
};

const NEAR_EE: &str = "two edges within the resolution of each other (S9d.2)";
const NEAR_LIMIT: &str = "a near miss's least distance within rounding of the resolution";

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

fn scale(a: [f64; 3], k: f64) -> [f64; 3] {
    a.map(|x| k * x)
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

/// The frame of a rim through `q` whose tangent there is normal to both the
/// unit `t` (the gap's direction) and the unit `along` (the other rim's
/// tangent), its face directions from `q` leaving along `t`: the wall
/// along `(along + t) / sqrt(2)` (the frame's normal, its body above its
/// rim) and the cap toward its centre along `(t - along) / sqrt(2)`, the
/// centre `r` from `q` that way; the frame's `x` axis turned so that `q`
/// lies at the angle `angle` (off the rim's seam).
fn rim_frame(q: [f64; 3], along: [f64; 3], t: [f64; 3], r: f64, angle: f64) -> Frame3 {
    let n = unit(add(along, t, 1.0));
    let inward = unit(add(t, along, -1.0));
    let e1 = scale(inward, -1.0);
    let e2 = cross(n, e1);
    let x = add(scale(e1, angle.cos()), e2, -angle.sin());
    frame(add(q, inward, r), n, x)
}

/// A rod of radius `r` and height 1 whose bottom rim passes through `q`
/// (`rim_frame`).
fn rim_rod(q: [f64; 3], along: [f64; 3], t: [f64; 3], r: f64) -> Solid {
    Solid::cylinder_with(
        OperationId(2),
        rim_frame(q, along, t, r, 1.0),
        r,
        0.0,
        1.0,
        tol(),
    )
    .unwrap()
    .0
}

/// A block with one rounded vertical edge (radius `r`), its bottom cap's
/// arc edge through `q` at the arc's middle (`rim_frame`): the profile
/// `(0, -r)`, an arc about the origin to `(r, 0)`, then `(r, 2)`, `(-2, 2)`
/// and `(-2, -r)`, its arc's middle `(1, -1) / sqrt(2)` toward `q`.
fn rounded_block(q: [f64; 3], along: [f64; 3], t: [f64; 3], r: f64) -> Solid {
    let f = rim_frame(q, along, t, r, -std::f64::consts::FRAC_PI_4);
    let points = [(0.0, -r), (r, 0.0), (r, 2.0), (-2.0, 2.0), (-2.0, -r)]
        .map(|(x, y)| Point2::new(x, y))
        .to_vec();
    let mut segments = vec![Segment::Line; 5];
    segments[0] = Segment::Arc {
        center: Point2::new(0.0, 0.0),
        radius: r,
        ccw: true,
    };
    let outer = Boundary::path(points, segments, tol()).unwrap();
    Solid::extrude_with(
        OperationId(2),
        Profile::new(outer, vec![], tol()).unwrap(),
        f,
        0.0,
        1.0,
    )
    .unwrap()
    .0
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
/// the new rule's reason or an existing rule's), the reasons returned.
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

/// No operation refused with the new rule's reasons (other outcomes as
/// they are).
fn not_near(what: &str, a: &Solid, b: &Solid) {
    for (op, out) in [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ] {
        match &out {
            Err(Error::Degenerate(m)) => assert!(*m != NEAR_EE, "{what} {op}: {m}"),
            Err(Error::ComputationLimit(m)) => assert!(*m != NEAR_LIMIT, "{what} {op}: {m}"),
            _ => {}
        }
    }
}

/// The misses within the resolution (`1e-7`).
const MISSES: [f64; 3] = [1e-12, 1e-11, 1e-9];

/// A rod's top rim (radius 1, the frame's `z` from -3 to 0) at the angle
/// 0.7: its point, its tangent and the direction leaving both its faces
/// (between its cap and its wall).
fn rod_rim(f: &Frame3) -> (Solid, [f64; 3], [f64; 3], [f64; 3]) {
    let rod = Solid::cylinder_with(OperationId(1), *f, 1.0, -3.0, 0.0, tol())
        .unwrap()
        .0;
    let (cs, sn) = (0.7f64.cos(), 0.7f64.sin());
    (
        rod,
        at(f, [cs, sn, 0.0]),
        unit(dir(f, [-sn, cs, 0.0])),
        unit(dir(f, [cs, sn, 1.0])),
    )
}

/// Two rods' rims crossing each other `d` apart (their tangents at right
/// angles, the gap normal to both, each rim's faces leaving it away from
/// the other) are `Degenerate`: the fuse was two solids and the common
/// empty where the inputs lie no farther apart than the resolution; so is a
/// rod's rim crossing the other rod's that deep. At `1e-6` every operation
/// evaluates.
#[test]
fn two_rims_within_the_resolution_of_each_other_are_degenerate() {
    for (name, f) in frames() {
        let (rod, p, along, t) = rod_rim(&f);
        let other = |d: f64| rim_rod(add(p, t, d), along, t, 0.5);
        for d in MISSES {
            refused(&format!("{name} rims {d:e}"), &rod, &other(d), NEAR_EE);
            degenerate(&format!("{name} rims {:e}", -d), &rod, &other(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} rims"), &rod, &other(1e-6), 2),
            0.0
        );
    }
}

/// A rod's rim by a cone's base rim (a frustum of radii 1 and 0.5 over a
/// height of 1), and by a block's arc edge (a rounded vertical edge of
/// radius 0.5, on its bottom cap), `d` apart and crossing that deep, are
/// `Degenerate`; at `1e-6` every operation evaluates.
#[test]
fn a_rim_within_the_resolution_of_a_cone_rim_or_an_arc_edge_is_degenerate() {
    for (name, f) in frames() {
        let cone = Solid::cone_with(OperationId(1), f, 1.0, 0.5, 1.0, tol())
            .unwrap()
            .0;
        let (cs, sn) = (0.7f64.cos(), 0.7f64.sin());
        let p = at(&f, [cs, sn, 0.0]);
        let along = unit(dir(&f, [-sn, cs, 0.0]));
        // Leaving the base (toward the axis) and the wall (toward the top)
        // away.
        let t = unit(dir(&f, [cs, sn, -0.5]));
        let by_cone = |d: f64| rim_rod(add(p, t, d), along, t, 0.5);
        for d in MISSES {
            refused(&format!("{name} cone {d:e}"), &cone, &by_cone(d), NEAR_EE);
            degenerate(&format!("{name} cone {:e}", -d), &cone, &by_cone(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} cone"), &cone, &by_cone(1e-6), 2),
            0.0
        );
        let (rod, p, along, t) = rod_rim(&f);
        let block = |d: f64| rounded_block(add(p, t, d), along, t, 0.5);
        for d in MISSES {
            refused(&format!("{name} arc {d:e}"), &rod, &block(d), NEAR_EE);
            degenerate(&format!("{name} arc {:e}", -d), &rod, &block(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} arc"), &rod, &block(1e-6), 2),
            0.0
        );
    }
}

/// What is no near miss of two rims' own evaluates as before: rims `0.1`
/// apart and `0.1` through, `1.5e-7` apart (beyond the resolution), two rims crossing exactly (but for rounding),
/// coaxial rims on one plane (two rods stacked, of one radius and of two),
/// U1's box on a rod's tangent plane turned by 30 degrees.
#[test]
fn rims_apart_crossing_or_on_one_plane_evaluate_as_before() {
    for (name, f) in frames() {
        let (rod, p, along, t) = rod_rim(&f);
        let other = |d: f64| rim_rod(add(p, t, d), along, t, 0.5);
        assert_eq!(
            evaluates(&format!("{name} apart"), &rod, &other(0.1), 2),
            0.0
        );
        evaluates(&format!("{name} through"), &rod, &other(-0.1), 1);
        // Beyond the resolution but within twice it: apart, bracketed.
        assert_eq!(
            evaluates(&format!("{name} beyond"), &rod, &other(1.5e-7), 2),
            0.0
        );
        not_near(&format!("{name} crossing"), &rod, &other(0.0));
        // Stacked on the rod's top, coaxial.
        for r in [1.0, 0.5, 2.0] {
            let top = Solid::cylinder_with(OperationId(2), f, r, 0.0, 1.0, tol())
                .unwrap()
                .0;
            not_near(&format!("{name} stacked {r}"), &rod, &top);
        }
    }
    let level = frames()[0].1;
    let rod = Solid::cylinder_with(OperationId(1), level, 1.0, 0.0, 2.0, tol())
        .unwrap()
        .0;
    let r = 3f64.sqrt() / 2.0;
    let b = Solid::box_at(
        Point3::new(-r, -1.0, 0.0),
        Vec3::new(2.0 * r, 1.0 + r, 2.0),
        tol(),
    )
    .unwrap();
    let turn = RigidTransform::rotation(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        30f64.to_radians(),
    )
    .unwrap();
    let u1 = b.transform_with(OperationId(8), turn).unwrap().0;
    evaluates("U1", &rod, &u1, 1);
}
