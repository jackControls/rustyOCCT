//! S9d.2's near misses at edges of other curves, between points of two
//! fields and at a cone's apex against a spline wall (REVIEW_NOTES.md,
//! "Other curves, points of two fields and an apex against a spline or
//! torus face"; `curved/near.rs`): a cube's corner or a block's edge within
//! the resolution of a given result's meeting of two curved faces, a cone's
//! or a torus's section or a spline prism's cap edge, two balls' zones
//! whose rims pass within it of each other, a given result's vertex of one
//! quadratic field by another's of another field or by a zone's rim, and a
//! cone's apex by a spline prism's wall, each input's boundary leaving its
//! point strictly away from the other's (missing, the gap outside either
//! input) or heading strictly toward it (crossing), are `Degenerate`, in a
//! level frame, one turned about the world's z and one tilted off it.
//! Incidences but for rounding, pairs well apart, U1's box on a rod's
//! tangent plane and an apex by a ball or a torus evaluate or are refused
//! as before.
use rusty_occt::identity::OperationId;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, RigidTransform, Segment,
    Solid, Tolerance, Vec3,
};

const NEAR_VV: &str = "a vertex within the resolution of a vertex of the other input (S9d.2)";
const NEAR_VE: &str = "a vertex within the resolution of an edge of the other input (S9d.2)";
const NEAR_EE: &str = "two edges within the resolution of each other (S9d.2)";
const NEAR_APEX: &str = "a cone's apex within the resolution of a face (S9d.3)";
const NEAR_SPHERE_VERTEX: &str = "a sphere within the resolution of a vertex (S9d.1)";
const NEAR_LIMIT: &str = "a near miss's least distance within rounding of the resolution";
const NEAR_FIELDS: &str = "a near miss's gap along a face between points of two fields";

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

/// A ball's zone of radius `r` between the latitudes `lo` and `hi` whose
/// bottom rim passes through `q` with its tangent there normal to the unit
/// `t` (the gap's direction) and the unit `along`, its cap and its wall
/// leaving `q` along `t` (the frame's normal `(along + t) / sqrt(2)`).
fn rim_zone(q: [f64; 3], along: [f64; 3], t: [f64; 3], r: f64, lo: f64, hi: f64) -> Solid {
    let n = unit(add(along, t, 1.0));
    let inward = unit(add(t, along, -1.0));
    let c = add(add(q, inward, r * lo.cos()), n, -r * lo.sin());
    Solid::sphere_with(OperationId(2), frame(c, n, inward), r, lo, hi, tol())
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

/// Every operation refused (a crossing within the resolution: the new
/// rules' reason or an existing rule's), the reasons returned.
fn degenerate(what: &str, a: &Solid, b: &Solid) -> Vec<String> {
    [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ]
    .into_iter()
    .map(|(op, out)| match out {
        Err(Error::Degenerate(m)) => m.to_string(),
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

/// No operation refused with the new rules' reasons (other outcomes as
/// they are).
fn not_near(what: &str, a: &Solid, b: &Solid) {
    for (op, out) in [
        ("fuse", a.fuse(OperationId(3), b)),
        ("cut", a.cut(OperationId(4), b)),
        ("cut back", b.cut(OperationId(4), a)),
        ("common", a.common(OperationId(5), b)),
    ] {
        match &out {
            Err(Error::Degenerate(m)) => {
                assert!(![NEAR_VV, NEAR_VE, NEAR_EE].contains(m), "{what} {op}: {m}")
            }
            Err(Error::ComputationLimit(m)) => {
                assert!(![NEAR_LIMIT, NEAR_FIELDS].contains(m), "{what} {op}: {m}")
            }
            _ => {}
        }
    }
}

/// The misses within the resolution (`1e-7`).
const MISSES: [f64; 3] = [1e-12, 1e-11, 1e-9];

/// A given result's edge of another curve at a point: the body, the
/// point, the curve's tangent there and the direction leaving both its
/// faces (their outward normals' bisector), all in the world.
struct OnCurve {
    what: &'static str,
    body: Solid,
    p: [f64; 3],
    tangent: [f64; 3],
    t: [f64; 3],
}

/// The edges of other curves in a frame `f` (local coordinates): a rod's
/// hole's rim (a rod less a crossing rod: a meeting of two cylinders) and
/// a dimple's rim (a rod less a ball: a cylinder's meeting with a sphere)
/// at `(1, 0, 0.5)`, their faces' normals `(1, 0, 0)` and `(0, 0, -1)`; a
/// cone's section by a plane through `(0, 0, 0.4)` normal to `(0, 0.6,
/// 0.8)` at `(0.6, 0, 0.4)`; a torus's section by such a plane through `(0,
/// 0, 0.5)` at `(3 + sqrt(0.75), 0, 0.5)`; a spline dome prism's top cap
/// edge at its apex `(2, 2, 2)`. Built once in each frame for the file's
/// tests.
fn curves(f: &Frame3) -> &'static [OnCurve] {
    static BUILT: [std::sync::OnceLock<Vec<OnCurve>>; 3] = [
        std::sync::OnceLock::new(),
        std::sync::OnceLock::new(),
        std::sync::OnceLock::new(),
    ];
    let k = frames()
        .iter()
        .position(|(_, g)| g == f)
        .expect("one of the frames");
    BUILT[k].get_or_init(|| build_curves(f))
}

fn build_curves(f: &Frame3) -> Vec<OnCurve> {
    let pi = std::f64::consts::PI;
    let mut out = Vec::new();
    let body = rod(*f, 1.0, -3.0, 3.0, 1);
    let bar = rod(
        local_frame(f, [-3.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        0.5,
        0.0,
        6.0,
        2,
    );
    let bite = Solid::sphere_with(
        OperationId(2),
        local_frame(f, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        0.5,
        -pi / 2.0,
        pi / 2.0,
        tol(),
    )
    .unwrap()
    .0;
    for (what, tool) in [("hole's rim", &bar), ("dimple's rim", &bite)] {
        out.push(OnCurve {
            what,
            body: body.cut(OperationId(6), tool).unwrap().0.remove(0),
            p: at(f, [1.0, 0.0, 0.5]),
            tangent: unit(dir(f, [0.0, 1.0, 0.0])),
            t: unit(dir(f, [1.0, 0.0, -1.0])),
        });
    }
    let plane = |o: f64, half: f64| {
        prism(
            local_frame(f, [0.0, 0.0, o], [0.0, 0.6, 0.8], [1.0, 0.0, 0.0]),
            rect(-half, -half, half, half),
            0.0,
            3.0,
            2,
        )
    };
    let cone = Solid::cone_with(OperationId(1), *f, 1.0, 0.0, 1.0, tol())
        .unwrap()
        .0;
    let h = std::f64::consts::FRAC_1_SQRT_2;
    out.push(OnCurve {
        what: "cone's section",
        body: cone
            .cut(OperationId(6), &plane(0.4, 3.0))
            .unwrap()
            .0
            .remove(0),
        p: at(f, [0.6, 0.0, 0.4]),
        tangent: unit(dir(f, [0.0, 0.8, -0.6])),
        t: unit(dir(f, [h, 0.6, h + 0.8])),
    });
    let torus = Solid::torus_with(OperationId(1), *f, 3.0, 1.0, -pi, pi, 2.0 * pi, tol())
        .unwrap()
        .0;
    let s = 0.75f64.sqrt();
    let n1 = [s, 0.0, 0.5];
    let n2 = [0.0, 0.6, 0.8];
    out.push(OnCurve {
        what: "torus's section",
        body: torus
            .cut(OperationId(6), &plane(0.5, 6.0))
            .unwrap()
            .0
            .remove(0),
        p: at(f, [3.0 + s, 0.0, 0.5]),
        tangent: unit(dir(f, cross(n1, n2))),
        t: unit(dir(f, add(n1, n2, 1.0))),
    });
    out.push(OnCurve {
        what: "spline's cap edge",
        body: dome(f),
        p: at(f, [2.0, 2.0, 2.0]),
        tangent: unit(dir(f, [1.0, 0.0, 0.0])),
        t: unit(dir(f, [0.0, 1.0, 1.0])),
    });
    out
}

/// A spline dome `v = u (4 - u) / 2` over `[0, 4]`, heights 0 to 2.
fn dome(f: &Frame3) -> Solid {
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
    prism(*f, outer, 0.0, 2.0, 1)
}

/// A cube's corner within the resolution of a given result's meeting of
/// two curved faces (a hole's rim, a dimple's rim), a cone's or a torus's
/// oblique section or a spline prism's cap edge, diagonally off it (its
/// faces' and the corner's edges leaving the gap away), is `Degenerate`:
/// the fuse was two solids and the common empty where the inputs lie no
/// farther apart than the resolution. A corner through the edge that deep
/// is refused too; at `1e-6` every operation evaluates.
#[test]
fn a_corner_off_another_curve_is_degenerate() {
    for (name, f) in frames() {
        for c in curves(&f) {
            let what = format!("{name} {}", c.what);
            let corner = |d: f64| corner_cube(add(c.p, c.t, d), c.t, 0.25);
            for d in MISSES {
                refused(&format!("{what} {d:e}"), &c.body, &corner(d), NEAR_VE);
                degenerate(&format!("{what} {:e}", -d), &c.body, &corner(-d));
            }
            assert_eq!(evaluates(&what, &c.body, &corner(1e-6), 2), 0.0);
        }
    }
}

/// A small rod (radius and height 0.1) whose bottom rim passes through `q`
/// with its tangent there normal to the unit `t` (the gap's direction) and
/// the unit `along` (the other curve's tangent), its cap and its wall
/// leaving `q` along `t`; `q` at the angle 1 of its frame (off its seam).
fn rim_rod(q: [f64; 3], along: [f64; 3], t: [f64; 3]) -> Solid {
    let r = 0.1;
    let n = unit(add(along, t, 1.0));
    let inward = unit(add(t, along, -1.0));
    let e1 = inward.map(|x| -x);
    let e2 = cross(n, e1);
    let x = add(e1.map(|c| c * 1f64.cos()), e2, -1f64.sin());
    rod(frame(add(q, inward, r), n, x), r, 0.0, r, 2)
}

/// A block's edge across a hole's rim or a spline prism's cap edge (normal
/// to the curve's tangent and to the gap), and a small rod's rim across a
/// hole's or a dimple's rim or a spline prism's cap edge (a conic against
/// another curve), `d` off it with their faces leaving it away, are
/// `Degenerate`; so are they crossing it that deep; at `1e-6` every
/// operation evaluates.
#[test]
fn an_edge_across_another_curve_is_degenerate() {
    for (name, f) in frames() {
        for c in curves(&f) {
            if !["hole's rim", "dimple's rim", "spline's cap edge"].contains(&c.what) {
                continue;
            }
            let e = unit(cross(c.t, c.tangent));
            let across = |d: f64| edge_block(add(c.p, c.t, d), e, c.t, 0.1);
            let rim = |d: f64| rim_rod(add(c.p, c.t, d), c.tangent, c.t);
            let tools: [(&str, &dyn Fn(f64) -> Solid); 2] = [("edge", &across), ("rim", &rim)];
            for (kind, tool) in tools {
                if kind == "edge" && c.what == "dimple's rim" {
                    continue;
                }
                let what = format!("{name} {} {kind}", c.what);
                for d in MISSES {
                    refused(&format!("{what} {d:e}"), &c.body, &tool(d), NEAR_EE);
                    degenerate(&format!("{what} {:e}", -d), &c.body, &tool(-d));
                }
                assert_eq!(evaluates(&what, &c.body, &tool(1e-6), 2), 0.0);
            }
        }
    }
}

/// A ball's zone (radius 1, latitudes -0.3 to 0.4) and its top rim's point
/// at the angle 0.7, the rim's tangent there and the direction leaving its
/// cap and its wall (`(cos 0.7, sin 0.7, 2) / sqrt(5)`, local).
fn zone_rim(f: &Frame3) -> (Solid, [f64; 3], [f64; 3], [f64; 3]) {
    let zone = Solid::sphere_with(OperationId(1), *f, 1.0, -0.3, 0.4, tol())
        .unwrap()
        .0;
    let (th, al) = (0.4f64, 0.7f64);
    let (cs, sn) = (al.cos(), al.sin());
    (
        zone,
        at(f, [th.cos() * cs, th.cos() * sn, th.sin()]),
        unit(dir(f, [-sn, cs, 0.0])),
        unit(dir(f, [cs, sn, 2.0])),
    )
}

/// Two balls' zones whose rims cross at right angles `d` apart (sphere
/// circles: points of two quadratic fields), each rim's faces leaving it
/// away from the other, are `Degenerate`: the fuse was two solids and the
/// common empty; so are rims crossing that deep. At `1e-6` every operation
/// evaluates.
#[test]
fn two_sphere_circles_within_the_resolution_are_degenerate() {
    for (name, f) in frames() {
        let (zone, p, along, t) = zone_rim(&f);
        let other = |d: f64| rim_zone(add(p, t, d), along, t, 0.6, -0.2, 0.3);
        for d in MISSES {
            refused(&format!("{name} zones {d:e}"), &zone, &other(d), NEAR_EE);
            degenerate(&format!("{name} zones {:e}", -d), &zone, &other(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} zones"), &zone, &other(1e-6), 2),
            0.0
        );
    }
}

/// A block bitten by a ball at its corner, on a frame `g` whose origin is
/// the vertex where the ball's sphere meets the block's top front edge: the
/// block `[-x, 2 - x] x [0, 2] x [-2, 0]` (`x = 0.3 + sqrt(0.87)`) less a
/// ball of radius 1 about `(0.3 - x, 0.2, 0.3)`, the vertex a point of
/// `Q(sqrt(0.87))` (but for the placement's rounding), its edges leaving it
/// along `(1, 0, 0)`, the top's arc `(0.2, sqrt(0.87), 0)` and the front's
/// arc `(-0.3, 0, -sqrt(0.87))`, all against `(-1, -1, 1)`. `back`: the
/// same mirrored through the vertex with a ball of radius 0.9 (a point of
/// `Q(sqrt(0.68))`), its edges leaving along `(-1, -1, 1)`.
fn bitten(g: &Frame3, back: bool) -> Solid {
    let x = 0.3 + 0.87f64.sqrt();
    let (block, c, r) = if back {
        (
            prism(*g, rect(x - 2.0, -2.0, x, 0.0), 0.0, 2.0, 3),
            [0.68f64.sqrt(), -0.2, -0.3],
            0.9,
        )
    } else {
        (
            prism(*g, rect(-x, 0.0, 2.0 - x, 2.0), -2.0, 0.0, 1),
            [0.3 - x, 0.2, 0.3],
            1.0,
        )
    };
    let ball = Solid::sphere_with(
        OperationId(if back { 4 } else { 2 }),
        local_frame(g, c, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        r,
        -std::f64::consts::FRAC_PI_2,
        std::f64::consts::FRAC_PI_2,
        tol(),
    )
    .unwrap()
    .0;
    block
        .cut(OperationId(if back { 7 } else { 6 }), &ball)
        .unwrap()
        .0
        .remove(0)
}

/// A vector turned about the unit `k` by `a` (Rodrigues).
fn turned(v: [f64; 3], k: [f64; 3], a: f64) -> [f64; 3] {
    let (c, s) = (a.cos(), a.sin());
    let kv = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
    add(add(v.map(|x| x * c), cross(k, v), s), k, kv * (1.0 - c))
}

/// `bitten`'s frame in `f` (its vertex at the local `(1, 0.5, 0.25)`), the
/// direction leaving its vertex, and the frame of its mirrored partner
/// whose vertex lies `d` along it, turned about it by 0.5 (no face of one
/// parallel to the other's).
fn bitten_frames(f: &Frame3) -> (Frame3, [f64; 3], impl Fn(f64) -> Frame3) {
    let g = local_frame(f, [1.0, 0.5, 0.25], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let t = unit(dir(&g, [-1.0, -1.0, 1.0]));
    let (v, n, x) = (
        g.origin().to_array(),
        g.normal().to_array(),
        g.x().to_array(),
    );
    let partner = move |d: f64| frame(add(v, t, d), turned(n, t, 0.5), turned(x, t, 0.5));
    (g, t, partner)
}

/// A given result's vertex of one quadratic field within the resolution of
/// another's of another field (two blocks bitten by balls, their vertices
/// `d` apart along the direction leaving both) or of a ball's zone's rim (a
/// sphere circle) is `Degenerate`: the fuse was two solids and the common
/// empty. Crossing that deep is refused too; at `1e-6` every operation
/// evaluates.
#[test]
fn vertices_of_two_fields_within_the_resolution_are_degenerate() {
    for (name, f) in frames() {
        let (g, t, partner) = bitten_frames(&f);
        let body = bitten(&g, false);
        let back = |d: f64| bitten(&partner(d), true);
        for d in MISSES {
            refused(&format!("{name} vertices {d:e}"), &body, &back(d), NEAR_VV);
            degenerate(&format!("{name} vertices {:e}", -d), &body, &back(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} vertices"), &body, &back(1e-6), 2),
            0.0
        );
        let v = g.origin().to_array();
        let along = unit(cross(t, g.normal().to_array()));
        let zone = |d: f64| rim_zone(add(v, t, d), along, t, 0.5, -0.2, 0.3);
        for d in MISSES {
            refused(&format!("{name} zone {d:e}"), &body, &zone(d), NEAR_VE);
            degenerate(&format!("{name} zone {:e}", -d), &body, &zone(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} zone"), &body, &zone(1e-6), 2),
            0.0
        );
    }
}

/// A cone (radius 0.5 over a height of 1) whose apex stands `d` off a spline
/// prism's wall, opening away from it, is `Degenerate` with the apex rule's
/// reason: the fuse was two solids and the common empty; its apex that far
/// through the wall is refused too; at `1e-6` every operation evaluates.
/// Against a ball the apex is held by S9d.1's vertex rule, and against a
/// torus the torus's meeting with the cone refuses it, as before.
#[test]
fn an_apex_within_the_resolution_of_a_spline_wall_is_degenerate() {
    let pi = std::f64::consts::PI;
    for (name, f) in frames() {
        let apex_cone = |a: [f64; 3], open: [f64; 3]| {
            let h = if open[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            let g = local_frame(&f, add(a, open, 1.0), open.map(|x| -x), cross(open, h));
            Solid::cone_with(OperationId(8), g, 0.5, 0.0, 1.0, tol())
                .unwrap()
                .0
        };
        let wall = dome(&f);
        // Its axis leaning off the wall's normal (its base no plane along
        // the wall's axis), every ruling leaving the wall.
        let lean = [0.0, 1.0, 0.3].map(|x: f64| x / 1.09f64.sqrt());
        let off = |d: f64| apex_cone([2.0, 2.0 + d, 1.0], lean);
        for d in MISSES {
            refused(&format!("{name} wall {d:e}"), &wall, &off(d), NEAR_APEX);
            degenerate(&format!("{name} wall {:e}", -d), &wall, &off(-d));
        }
        assert_eq!(
            evaluates(&format!("{name} wall"), &wall, &off(1e-6), 2),
            0.0
        );
        let ball = Solid::sphere_with(OperationId(7), f, 1.0, -pi / 2.0, pi / 2.0, tol())
            .unwrap()
            .0;
        let torus = Solid::torus_with(OperationId(7), f, 3.0, 1.0, -pi, pi, 2.0 * pi, tol())
            .unwrap()
            .0;
        for d in [1e-12, -1e-12, 1e-9] {
            refused(
                &format!("{name} ball {d:e}"),
                &ball,
                &apex_cone([1.0 + d, 0.0, 0.0], [1.0, 0.0, 0.0]),
                NEAR_SPHERE_VERTEX,
            );
            let by_torus = apex_cone([4.0 + d, 0.0, 0.0], [1.0, 0.0, 0.0]);
            for out in [
                torus.fuse(OperationId(3), &by_torus),
                torus.common(OperationId(5), &by_torus),
            ] {
                assert!(
                    matches!(out, Err(Error::Degenerate(_) | Error::ComputationLimit(_))),
                    "{name} torus {d:e}: {:?}",
                    out.map(|x| x.0.len())
                );
            }
        }
    }
}

/// What is no near miss of the rules' own evaluates as before: a corner
/// on another curve but for rounding (an incidence), a corner `0.1` off
/// it and `1.5e-7` off it (beyond the resolution), zones' rims and vertices
/// of two fields `0.1` apart, and U1's box on a rod's tangent plane.
#[test]
fn incidences_and_pairs_apart_evaluate_as_before() {
    for (name, f) in frames() {
        for c in curves(&f) {
            let what = format!("{name} {}", c.what);
            let corner = |d: f64| corner_cube(add(c.p, c.t, d), c.t, 0.25);
            not_near(&format!("{what} on"), &c.body, &corner(0.0));
            assert_eq!(
                evaluates(&format!("{what} apart"), &c.body, &corner(0.1), 2),
                0.0
            );
            assert_eq!(
                evaluates(&format!("{what} beyond"), &c.body, &corner(1.5e-7), 2),
                0.0
            );
        }
        let (zone, p, along, t) = zone_rim(&f);
        let other = rim_zone(add(p, t, 0.1), along, t, 0.6, -0.2, 0.3);
        assert_eq!(
            evaluates(&format!("{name} zones apart"), &zone, &other, 2),
            0.0
        );
        let other = rim_zone(add(p, t, 1.5e-7), along, t, 0.6, -0.2, 0.3);
        assert_eq!(
            evaluates(&format!("{name} zones beyond"), &zone, &other, 2),
            0.0
        );
        not_near(
            &format!("{name} zones on"),
            &zone,
            &rim_zone(p, along, t, 0.6, -0.2, 0.3),
        );
        let (g, _, partner) = bitten_frames(&f);
        assert_eq!(
            evaluates(
                &format!("{name} vertices apart"),
                &bitten(&g, false),
                &bitten(&partner(0.1), true),
                2
            ),
            0.0
        );
    }
    let level = frames()[0].1;
    let rod = rod(level, 1.0, 0.0, 2.0, 1);
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
