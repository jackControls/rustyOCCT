//! S9's near misses of faces with no rule of their own (REVIEW_NOTES.md,
//! "Faces with no near-miss rule of their own"; `curved/near.rs`): two
//! cylinders' or cones' walls, a torus or a spline wall against another
//! input's face, and edges and vertices against a torus or a spline wall,
//! within the resolution of touching at a point, along a line or along a
//! circle, missing with the gap outside either input or crossing, are
//! `Degenerate`, in a level frame, one turned about the world's z and one
//! tilted off it. Exact tangencies and incidences (the faces' own rules),
//! pairs well apart and just beyond the resolution, and U1's and Y5's box
//! on a rod's tangent plane evaluate or are refused as before.
use rusty_occt::identity::OperationId;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, RigidTransform, Segment,
    Solid, Tolerance, Vec3,
};

const NEAR_WALLS: &str =
    "a cylinder's or cone's wall within the resolution of tangency to the other input's (S9c)";
const NEAR_TORUS: &str =
    "a torus within the resolution of tangency to the other input's surface (S9d.4)";
const NEAR_SPLINE: &str =
    "a spline wall within the resolution of tangency to the other input's surface (S9f)";
const NEAR_FACE: &str = "an edge within the resolution of tangency to a face (S9d.2)";
const NEAR_VERTEX_FACE: &str = "a vertex within the resolution of a face (S9d.2)";
const NEAR_LIMIT: &str = "a near miss's least distance within rounding of the resolution";
const NODE: &str = "two cylinders' section within the resolution of a node";
const TANGENCY: &str = "a tangency between the inputs (S9c)";
const TURNING: &str = "a torus meeting's turning point not covered";
const AXIS: &str = "a plane within rounding of a spline wall's axis";

/// The misses within the resolution (`1e-7`).
const MISSES: [f64; 3] = [1e-12, 1e-11, 1e-9];

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

/// A level frame, one turned about the world's z by 30 degrees and one
/// tilted off it (the boolean target's).
fn frames() -> [(&'static str, Frame3); 3] {
    [
        ("level", frame([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])),
        (
            "turned",
            frame([0.0; 3], [0.0, 0.0, 1.0], [3f64.sqrt() / 2.0, 0.5, 0.0]),
        ),
        (
            "tilted",
            frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
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

/// A frame at a frame's local point with its local normal and `x`.
fn local_frame(f: &Frame3, o: [f64; 3], n: [f64; 3], x: [f64; 3]) -> Frame3 {
    frame(at(f, o), dir(f, n), dir(f, x))
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Boundary {
    Boundary::polygon(
        [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
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

fn rod(f: Frame3, r: f64, lo: f64, hi: f64, id: u64) -> Solid {
    Solid::cylinder_with(OperationId(id), f, r, lo, hi, tol())
        .unwrap()
        .0
}

/// A ball of radius `r` about the local point `c` of `f`.
fn ball(f: &Frame3, c: [f64; 3], r: f64) -> Solid {
    let h = std::f64::consts::FRAC_PI_2;
    Solid::sphere_with(
        OperationId(2),
        local_frame(f, c, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        r,
        -h,
        h,
        tol(),
    )
    .unwrap()
    .0
}

/// A whole torus of radii 3 and 1 on a frame.
fn torus_on(f: Frame3, id: u64) -> Solid {
    let pi = std::f64::consts::PI;
    Solid::torus_with(OperationId(id), f, 3.0, 1.0, -pi, pi, 2.0 * pi, tol())
        .unwrap()
        .0
}

/// The cone of radius 1 over a height of 1 on a frame (its apex at `(0, 0,
/// 1)`), the point `(0.5, 0, 0.5)` of its wall and its outward normal there.
fn cone(f: &Frame3) -> Solid {
    Solid::cone_with(OperationId(1), *f, 1.0, 0.0, 1.0, tol())
        .unwrap()
        .0
}

const CONE_P: [f64; 3] = [0.5, 0.0, 0.5];

fn cone_n() -> [f64; 3] {
    unit([1.0, 0.0, 1.0])
}

/// A frustum (radii 0.05 to 0.35 over a length of 1) whose wall's
/// generatrix on its side toward `-m` runs along the unit `e` through the
/// local point `p` of `f` (from half its length before `p`), its axis
/// tilted from `e` toward the unit `m` (normal to `e`): its wall lies on
/// `m`'s side of the plane through that generatrix normal to `m`.
fn frustum(f: &Frame3, p: [f64; 3], e: [f64; 3], m: [f64; 3]) -> Solid {
    let (r0, r1, h) = (0.05f64, 0.35f64, 1.0f64);
    let ta = (r1 - r0) / h;
    let ca = 1.0 / (1.0 + ta * ta).sqrt();
    let sa = ta * ca;
    let a = add(e.map(|x| x * ca), m, sa);
    let mm = add(e.map(|x| -x * sa), m, ca);
    let b0 = add(add(p, e, -0.5 * ca), mm, r0);
    Solid::cone_with(OperationId(2), local_frame(f, b0, a, mm), r0, r1, h, tol())
        .unwrap()
        .0
}

/// A spline dome `v = u (4 - u) / 2` over `[0, 4]`, heights 0 to 2, on a
/// frame: its wall's crest `(2, 2)` with its outward normal `(0, 1)`.
fn dome(f: &Frame3, id: u64) -> Solid {
    let curve = BSplineCurve2::new(
        2,
        vec![
            Point2::new(4.0, 0.0),
            Point2::new(2.0, 4.0),
            Point2::new(0.0, 0.0),
        ],
        None,
        vec![0.0, 1.0],
        vec![3, 3],
    )
    .unwrap();
    let outer = Boundary::path(
        vec![Point2::new(0.0, 0.0), Point2::new(4.0, 0.0)],
        vec![
            Segment::Line,
            Segment::Spline(rusty_occt::topology::SplineSpan::whole(curve)),
        ],
        tol(),
    )
    .unwrap();
    prism(*f, outer, 0.0, 2.0, id)
}

/// A cube of side `s` on its corner `o`, each of its three edges from it
/// leaving along the unit `t` (at `1 / sqrt(3)` of its length), turned
/// about `t` by 0.4 from a basis of `t` and the world's axes.
fn corner_cube(o: [f64; 3], t: [f64; 3], s: f64) -> Solid {
    let h = if t[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let p1 = unit(cross(t, h));
    let q = cross(t, p1);
    let (c, sn) = (0.4f64.cos(), 0.4f64.sin());
    let p1 = add(p1.map(|x| c * x), q, sn);
    let p2 = cross(t, p1);
    let n = add(t, p2, 1.0);
    let x = add(add(t, p1, -1.0), p2, -1.0);
    prism(frame(o, n, x), rect(0.0, 0.0, s, s), 0.0, s, 9)
}

/// A block standing on its edge through `o` along the unit `e` (from `-s`
/// to `s` along it), its faces leaving the edge along `(t + p) / sqrt(2)`
/// and `(t - p) / sqrt(2)`, `s` wide, `(e, t, p)` with `p = t x e`.
fn edge_block(o: [f64; 3], e: [f64; 3], t: [f64; 3], s: f64) -> Solid {
    let p = cross(t, e);
    let x = unit(add(t, p, 1.0));
    prism(frame(o, e, x), rect(0.0, 0.0, s, s), -s, s, 9)
}

/// Every operation, either way round, refused with `why`, or with one of
/// `earlier` (the faces' own meetings' refusals or limits, first where
/// they hold).
fn refused(what: &str, a: &Solid, b: &Solid, why: &str, earlier: &[&str]) {
    for (op, out) in [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ] {
        assert!(
            match &out {
                Err(Error::Degenerate(m)) => *m == why || earlier.contains(m),
                Err(Error::ComputationLimit(m)) => earlier.contains(m),
                _ => false,
            },
            "{what} {op}: {:?}",
            out.map(|x| x.0.len())
        );
    }
}

/// Every operation refused (a crossing within the resolution: the new
/// rules' reason or an existing rule's, `Degenerate` or a
/// `ComputationLimit` of the meetings).
fn crossing_refused(what: &str, a: &Solid, b: &Solid) {
    for (op, out) in [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ] {
        assert!(
            matches!(&out, Err(Error::Degenerate(_) | Error::ComputationLimit(_))),
            "{what} {op}: {:?}",
            out.map(|x| x.0.len())
        );
    }
}

fn volume(s: &[Solid]) -> f64 {
    s.iter().map(|x| x.mass_properties().volume).sum()
}

/// Fuse, both cuts and common evaluate, the fuse `fused` solids, in the
/// pair identities; the common's volume.
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

/// No operation refused with the new rules' reasons or their limit (other
/// outcomes as they are).
fn not_near(what: &str, a: &Solid, b: &Solid) {
    for (op, out) in [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ] {
        match &out {
            Err(Error::Degenerate(m)) => assert!(
                ![
                    NEAR_WALLS,
                    NEAR_TORUS,
                    NEAR_SPLINE,
                    NEAR_FACE,
                    NEAR_VERTEX_FACE
                ]
                .contains(m),
                "{what} {op}: {m}"
            ),
            Err(Error::ComputationLimit(m)) => assert!(*m != NEAR_LIMIT, "{what} {op}: {m}"),
            _ => {}
        }
    }
}

/// A configuration placed at a gap `d` (negative: crossing that deep) of
/// a fixed body: what it is, the fixed body, the placed body and the
/// rule's reason.
type Placed<'a> = dyn Fn(f64) -> Solid + 'a;

/// Every operation a `ComputationLimit` with the rule's limit.
fn limited(what: &str, a: &Solid, b: &Solid) {
    for (op, out) in [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ] {
        assert!(
            matches!(&out, Err(Error::ComputationLimit(m)) if *m == NEAR_LIMIT),
            "{what} {op}: {:?}",
            out.map(|x| x.0.len())
        );
    }
}

/// How a configuration `1.5e-7` apart (beyond the resolution, within twice
/// it) is taken: certified apart (no new refusal), or a contact along a
/// cone's ruling or a torus's circle the interval bound does not reach
/// (`ComputationLimit`).
#[derive(Clone, Copy)]
enum Beyond {
    Apart,
    Limit,
}

/// The rule's evidence for one configuration: every miss within the
/// resolution refused with `why`, the crossings `1e-12` and `1e-9` deep
/// refused, the pair evaluating at `1e-6` (`fused` solids), and at `1.5e-7`
/// as `beyond`.
fn check(what: &str, body: &Solid, placed: &Placed, why: &str, fused: usize, beyond: Beyond) {
    check_or(what, body, placed, why, fused, beyond, &[], true);
}

/// `check`, the faces' own meetings' refusals or limits `earlier` first
/// where they hold.
fn check_or(
    what: &str,
    body: &Solid,
    placed: &Placed,
    why: &str,
    fused: usize,
    beyond: Beyond,
    earlier: &[&str],
    deep: bool,
) {
    for d in MISSES {
        refused(&format!("{what} {d:e}"), body, &placed(d), why, earlier);
        if d != 1e-11 && (deep || d < 1e-9) {
            crossing_refused(&format!("{what} {:e}", -d), body, &placed(-d));
        }
    }
    match beyond {
        Beyond::Apart => not_near(&format!("{what} beyond"), body, &placed(1.5e-7)),
        Beyond::Limit => limited(&format!("{what} beyond"), body, &placed(1.5e-7)),
    }
    evaluates(&format!("{what} 1e-6"), body, &placed(1e-6), fused);
}

/// Two rods crossing at right angles whose walls miss by `d` across their
/// common perpendicular (S9c.2a's extents decided them apart before the
/// node rule) and a rod across a cone's wall: `Degenerate`, the fuse two
/// solids and the common empty before; their crossings refused; at `1e-6`
/// every operation evaluates. Rods touching exactly are the faces'
/// tangency, as before.
#[test]
fn cylinders_and_cones_within_the_resolution_at_a_point_are_degenerate() {
    for (name, f) in frames() {
        let a = rod(f, 1.0, -2.0, 2.0, 1);
        let across = |d: f64| {
            rod(
                local_frame(&f, [-3.0, 1.5 + d, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                0.5,
                0.0,
                6.0,
                2,
            )
        };
        check(
            &format!("{name} rods"),
            &a,
            &across,
            NEAR_WALLS,
            2,
            Beyond::Apart,
        );
        not_near(&format!("{name} rods touching"), &a, &across(0.0));
        let c = cone(&f);
        let n = cone_n();
        let by_cone = |d: f64| {
            rod(
                local_frame(
                    &f,
                    add(add(CONE_P, n, 0.25 + d), [0.0, 1.0, 0.0], -1.0),
                    [0.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0],
                ),
                0.25,
                0.0,
                2.0,
                2,
            )
        };
        check(
            &format!("{name} rod by cone"),
            &c,
            &by_cone,
            NEAR_WALLS,
            2,
            Beyond::Apart,
        );
    }
}

/// Frusta by a cone's and a rod's walls at a point: `Degenerate`, the fuse
/// two solids and the common empty before.
#[test]
fn frusta_within_the_resolution_at_a_point_are_degenerate() {
    for (name, f) in frames() {
        let a = rod(f, 1.0, -2.0, 2.0, 1);
        let c = cone(&f);
        let n = cone_n();
        let frustum_by_cone = |d: f64| frustum(&f, add(CONE_P, n, d), [0.0, 1.0, 0.0], n);
        check(
            &format!("{name} frustum by cone"),
            &c,
            &frustum_by_cone,
            NEAR_WALLS,
            2,
            Beyond::Apart,
        );
        let frustum_by_rod =
            |d: f64| frustum(&f, [1.0 + d, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        check(
            &format!("{name} frustum by rod"),
            &a,
            &frustum_by_rod,
            NEAR_WALLS,
            2,
            Beyond::Apart,
        );
    }
}

/// A rod along another in a frame turned about its axis (parallel axes: a
/// contact along a ruling) and a box's face along a rod's ruling:
/// `Degenerate`; the faces' own tangency where they touch exactly.
#[test]
fn cylinders_along_a_ruling_within_the_resolution_are_degenerate() {
    for (name, f) in frames() {
        let a = rod(f, 1.0, -2.0, 2.0, 1);
        let g = local_frame(&f, [0.0; 3], [0.0, 0.0, 1.0], [0.6, 0.8, 0.0]);
        let along = |d: f64| {
            rod(
                local_frame(&g, [1.5 + d, 0.0, -1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                0.5,
                0.0,
                2.0,
                2,
            )
        };
        check(
            &format!("{name} rods along"),
            &a,
            &along,
            NEAR_WALLS,
            2,
            Beyond::Apart,
        );
        let face = |d: f64| {
            prism(
                local_frame(&f, [1.0 + d, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                rect(-0.5, -0.5, 0.5, 0.5),
                0.0,
                1.0,
                2,
            )
        };
        check(
            &format!("{name} box along rod"),
            &a,
            &face,
            NEAR_WALLS,
            2,
            Beyond::Apart,
        );
    }
}

/// A ball by a torus's outer equator and on its top and a second torus by
/// it (their outer equators): `Degenerate`, the fuse two solids and the
/// common empty before.
#[test]
fn a_torus_within_the_resolution_is_degenerate() {
    for (name, f) in frames() {
        let t = torus_on(f, 1);
        let by = |d: f64| ball(&f, [4.75 + d, 0.0, 0.0], 0.75);
        check(
            &format!("{name} ball by torus"),
            &t,
            &by,
            NEAR_TORUS,
            2,
            Beyond::Apart,
        );
        let on = |d: f64| ball(&f, [3.0, 0.0, 1.75 + d], 0.75);
        check(
            &format!("{name} ball on torus"),
            &t,
            &on,
            NEAR_TORUS,
            2,
            Beyond::Apart,
        );
        let second = |d: f64| {
            torus_on(
                local_frame(&f, [8.0 + d, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                2,
            )
        };
        check_or(
            &format!("{name} torus by torus"),
            &t,
            &second,
            NEAR_TORUS,
            2,
            Beyond::Apart,
            &[TURNING],
            true,
        );
    }
}

/// A frustum by a torus's equator, and a ball and a rod in its hole
/// (coaxial: a contact along a circle): `Degenerate`, the fuse two solids
/// and the common empty before.
#[test]
fn a_frustum_by_a_torus_or_a_body_in_its_hole_is_degenerate() {
    for (name, f) in frames() {
        let t = torus_on(f, 1);
        let frustum_by =
            |d: f64| frustum(&f, [4.0 + d, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        check(
            &format!("{name} frustum by torus"),
            &t,
            &frustum_by,
            NEAR_TORUS,
            2,
            Beyond::Apart,
        );
        let hole_ball = |d: f64| ball(&f, [0.0; 3], 2.0 - d);
        check(
            &format!("{name} ball in torus hole"),
            &t,
            &hole_ball,
            NEAR_TORUS,
            2,
            Beyond::Limit,
        );
        let hole_rod = |d: f64| {
            rod(
                local_frame(&f, [0.0, 0.0, -1.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                2.0 - d,
                0.0,
                3.0,
                2,
            )
        };
        // The rod crossing that deep is left out in the tilted frame.
        check_or(
            &format!("{name} rod in torus hole"),
            &t,
            &hole_rod,
            NEAR_TORUS,
            2,
            Beyond::Limit,
            &[],
            // A rod crossing the hole's equator by `1e-9` in the tilted
            // frame runs for minutes in the torus meeting's verification
            // (`ToricCrv::verified`), before any near-miss rule, at the base
            // as now (REVIEW_NOTES.md, "Faces with no near-miss rule of
            // their own, done").
            name != "tilted",
        );
    }
}

/// A rod by a torus's equator and a plane by it at a point: the torus's
/// near-node rule refuses them near tangency, with a band other than the
/// resolution: at `5e-8` and `9e-8` they evaluated (two solids); now
/// `Degenerate` there too, the node rule's reason first where it holds.
#[test]
fn a_torus_beyond_its_node_band_within_the_resolution_is_degenerate() {
    for (name, f) in frames() {
        let t = torus_on(f, 1);
        let by_rod = |d: f64| {
            rod(
                local_frame(&f, [4.5 + d, 0.0, -1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
                0.5,
                0.0,
                2.0,
                2,
            )
        };
        let n = [0.6, 0.0, 0.8];
        let by_plane = |d: f64| {
            prism(
                local_frame(&f, add([3.0, 0.0, 0.0], n, 1.0 + d), n, [0.0, 1.0, 0.0]),
                rect(-1.0, -1.0, 1.0, 1.0),
                0.0,
                1.0,
                2,
            )
        };
        let cases: [(&str, &Placed); 2] = [("rod", &by_rod), ("plane", &by_plane)];
        for (what, placed) in cases {
            let what = format!("{name} {what} by torus");
            // The node rule's tangency first where it holds (a band of its
            // own, the frame's), else the rule's.
            for d in [5e-8, 9e-8] {
                refused(
                    &format!("{what} {d:e}"),
                    &t,
                    &placed(d),
                    NEAR_TORUS,
                    &[NODE, TANGENCY],
                );
            }
            for d in MISSES {
                let out = t.fuse(OperationId(3), &placed(d));
                assert!(matches!(out, Err(Error::Degenerate(_))), "{what} {d:e}");
            }
            // The node rule's band reaches `1e-6` in some frames.
            not_near(&format!("{what} beyond"), &t, &placed(1.5e-7));
            not_near(&format!("{what} 1e-6"), &t, &placed(1e-6));
        }
    }
}

/// A rod across and along a spline dome's wall at its crest (in a frame
/// turned about its axis: an exactly parallel axis, S9f.2a), a ball by it,
/// a frustum across and along it, a slab on its crest's tangent plane and
/// a second dome along it in a frame turned by a half turn (parallel
/// axes): `Degenerate`, the fuse two solids and the common empty before.
#[test]
fn a_spline_wall_within_the_resolution_is_degenerate() {
    for (name, f) in frames() {
        let w = dome(&f, 1);
        let across = |d: f64| {
            rod(
                local_frame(&f, [1.0, 2.25 + d, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                0.25,
                0.0,
                2.0,
                2,
            )
        };
        check(
            &format!("{name} rod across spline"),
            &w,
            &across,
            NEAR_SPLINE,
            2,
            Beyond::Apart,
        );
        let along = |d: f64| {
            rod(
                local_frame(&f, [2.0, 2.25 + d, 0.5], [0.0, 0.0, 1.0], [0.6, 0.8, 0.0]),
                0.25,
                0.0,
                1.0,
                2,
            )
        };
        check(
            &format!("{name} rod along spline"),
            &w,
            &along,
            NEAR_SPLINE,
            2,
            Beyond::Apart,
        );
        let by = |d: f64| ball(&f, [2.0, 2.5 + d, 1.0], 0.5);
        check(
            &format!("{name} ball by spline"),
            &w,
            &by,
            NEAR_SPLINE,
            2,
            Beyond::Apart,
        );
        let tilt = unit([1.0, 0.0, 0.4]);
        let frustum_across = |d: f64| frustum(&f, [2.0, 2.0 + d, 1.0], tilt, [0.0, 1.0, 0.0]);
        check(
            &format!("{name} frustum across spline"),
            &w,
            &frustum_across,
            NEAR_SPLINE,
            2,
            Beyond::Apart,
        );
        let frustum_along =
            |d: f64| frustum(&f, [2.0, 2.0 + d, 1.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]);
        check(
            &format!("{name} frustum along spline"),
            &w,
            &frustum_along,
            NEAR_SPLINE,
            2,
            Beyond::Limit,
        );
        let slab = |d: f64| {
            prism(
                local_frame(&f, [2.0, 2.0 + d, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]),
                rect(-0.5, -0.5, 0.5, 0.5),
                0.0,
                0.5,
                2,
            )
        };
        // Off the level frame the slab's side faces hold the wall's axis
        // but for rounding: refused by that rule from some gap on, `1e-6`
        // among them, as before.
        if name == "level" {
            check(
                &format!("{name} slab on spline"),
                &w,
                &slab,
                NEAR_SPLINE,
                2,
                Beyond::Apart,
            );
        } else {
            for d in MISSES {
                refused(
                    &format!("{name} slab on spline {d:e}"),
                    &w,
                    &slab(d),
                    NEAR_SPLINE,
                    &[AXIS],
                );
            }
            not_near(&format!("{name} slab on spline 1e-6"), &w, &slab(1e-6));
        }
        let second = |d: f64| {
            dome(
                &local_frame(&f, [4.0, 4.0 + d, 0.5], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]),
                2,
            )
        };
        // The tilted frame's second dome's axis is its normal normalized
        // again, not exactly parallel: two spline walls on crossing axes,
        // refused by design.
        if name != "tilted" {
            check(
                &format!("{name} spline along spline"),
                &w,
                &second,
                NEAR_SPLINE,
                2,
                Beyond::Limit,
            );
        }
    }
}

/// A box's edge and a cube's corner by a torus's outer equator and by a
/// spline dome's crest, leaving the gap away: `Degenerate` with the edge
/// and vertex rules' reasons, the fuse two solids and the common empty
/// before.
#[test]
fn edges_and_vertices_within_the_resolution_of_a_torus_or_spline_are_degenerate() {
    for (name, f) in frames() {
        let t = torus_on(f, 1);
        let (p, out) = (at(&f, [4.0, 0.0, 0.0]), dir(&f, [1.0, 0.0, 0.0]));
        let e = dir(&f, [0.0, 0.0, 1.0]);
        let edge = |d: f64| edge_block(add(p, out, d), e, out, 0.5);
        check(
            &format!("{name} edge by torus"),
            &t,
            &edge,
            NEAR_FACE,
            2,
            Beyond::Apart,
        );
        let corner = |d: f64| corner_cube(add(p, out, d), out, 0.5);
        check(
            &format!("{name} corner by torus"),
            &t,
            &corner,
            NEAR_VERTEX_FACE,
            2,
            Beyond::Apart,
        );
        let w = dome(&f, 1);
        let (p, out) = (at(&f, [2.0, 2.0, 1.0]), dir(&f, [0.0, 1.0, 0.0]));
        let e = unit(dir(&f, [1.0, 0.0, 1.0]));
        let edge = |d: f64| edge_block(add(p, out, d), e, out, 0.3);
        check(
            &format!("{name} edge by spline"),
            &w,
            &edge,
            NEAR_FACE,
            2,
            Beyond::Apart,
        );
        let corner = |d: f64| corner_cube(add(p, out, d), out, 0.3);
        check(
            &format!("{name} corner by spline"),
            &w,
            &corner,
            NEAR_VERTEX_FACE,
            2,
            Beyond::Apart,
        );
    }
}

/// What is no near miss of the rules' own evaluates or is refused as
/// before: a ball, a rod and a second torus touching a torus exactly, rods
/// touching exactly, pairs `0.1` apart, U1's box overlapping a rod and Y5's
/// kind (a box on a rod's tangent plane turned by 30 degrees, its gap of
/// rounding's size), each in the pair identities where it evaluates.
#[test]
fn tangencies_incidences_and_pairs_apart_evaluate_as_before() {
    for (name, f) in frames() {
        let t = torus_on(f, 1);
        not_near(
            &format!("{name} ball on torus"),
            &t,
            &ball(&f, [4.75, 0.0, 0.0], 0.75),
        );
        assert_eq!(
            evaluates(
                &format!("{name} ball apart from torus"),
                &t,
                &ball(&f, [4.85, 0.0, 0.0], 0.75),
                2
            ),
            0.0
        );
        let w = dome(&f, 1);
        assert_eq!(
            evaluates(
                &format!("{name} rod apart from spline"),
                &w,
                &rod(
                    local_frame(&f, [1.0, 2.35, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                    0.25,
                    0.0,
                    2.0,
                    2,
                ),
                2
            ),
            0.0
        );
        let a = rod(f, 1.0, -2.0, 2.0, 1);
        assert_eq!(
            evaluates(
                &format!("{name} rods apart"),
                &a,
                &rod(
                    local_frame(&f, [-3.0, 1.6, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                    0.5,
                    0.0,
                    6.0,
                    2,
                ),
                2
            ),
            0.0
        );
    }
    let level = frames()[0].1;
    let rod1 = rod(level, 1.0, 0.0, 2.0, 1);
    let turn = RigidTransform::rotation(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        30f64.to_radians(),
    )
    .unwrap();
    let r = 3f64.sqrt() / 2.0;
    let u1 = Solid::box_at(
        Point3::new(-r, -1.0, 0.0),
        Vec3::new(2.0 * r, 1.0 + r, 2.0),
        tol(),
    )
    .unwrap()
    .transform_with(OperationId(8), turn)
    .unwrap()
    .0;
    evaluates("U1", &rod1, &u1, 1);
    let outside = Solid::box_at(
        Point3::new(-0.5, -2.0, 0.0),
        Vec3::new(1.0, 1.0, 2.0),
        tol(),
    )
    .unwrap()
    .transform_with(OperationId(8), turn)
    .unwrap()
    .0;
    assert_eq!(evaluates("Y5's kind", &rod1, &outside, 2), 0.0);
}
