//! S9d's near misses at cylinders, cones and faces' edges (REVIEW_NOTES.md,
//! "A ball within the resolution of a cylinder or cone, and edges within it
//! of a face"; `curved/near.rs`): within the resolution of tangency,
//! crossing is `Degenerate`, and missing where the gap lies outside either
//! input; a gap inside both inputs evaluates, and so does every case beyond
//! the resolution. `tests/sphere_booleans.rs` holds the earlier rules (a
//! sphere against a plane face, a sphere, an edge or a vertex).
use num_rational::BigRational as R;
use rusty_occt::identity::OperationId;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance,
    Vec3,
};

type V = [R; 3];

/// A ball placed at a gap `d` of its radius.
type Placed<'a> = dyn Fn(f64) -> Solid + 'a;

const NEAR_QUADRIC: &str =
    "a sphere within the resolution of tangency to a cylinder or cone (S9d.2)";
const NEAR_FACE: &str = "an edge within the resolution of tangency to a face (S9d.2)";
const NEAR_EDGE: &str = "a sphere within the resolution of tangency to an edge (S9d.1)";

fn q(x: f64) -> R {
    R::from_float(x).unwrap()
}

fn sub(a: &V, b: &V) -> V {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}

fn dot(a: &V, b: &V) -> R {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

fn cross(a: &V, b: &V) -> V {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

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

/// The level frame and the boolean target's tilted one.
fn frames() -> [(&'static str, Frame3); 2] {
    [
        ("level", frame([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])),
        (
            "turned",
            frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
        ),
    ]
}

/// A point's exact coordinates on a frame's stored axes.
fn local(f: &Frame3, p: [f64; 3]) -> V {
    let v = |a: [f64; 3]| a.map(q);
    let (o, x, y, n) = (
        v(f.origin().to_array()),
        v(f.x().to_array()),
        v(f.y().to_array()),
        v(f.normal().to_array()),
    );
    let det = dot(&x, &cross(&y, &n));
    let inv = [cross(&y, &n), cross(&n, &x), cross(&x, &y)];
    let d = sub(&v(p), &o);
    [0, 1, 2].map(|k| dot(&inv[k], &d) / &det)
}

fn world(f: &Frame3, u: f64, v: f64, w: f64) -> [f64; 3] {
    f.point(Point2::new(u, v), w).to_array()
}

fn ball_in(f: Frame3, r: f64) -> Solid {
    let half = std::f64::consts::FRAC_PI_2;
    Solid::sphere_with(OperationId(9), f, r, -half, half, tol())
        .unwrap()
        .0
}

fn ball(c: [f64; 3], r: f64) -> Solid {
    ball_in(frame(c, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]), r)
}

fn square(s: f64) -> Boundary {
    Boundary::polygon(
        [(-s, -s), (s, -s), (s, s), (-s, s)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec(),
        tol(),
    )
    .unwrap()
}

fn prism(f: Frame3, outer: Boundary, holes: Vec<Boundary>, lo: f64, hi: f64, id: u64) -> Solid {
    Solid::extrude_with(
        OperationId(id),
        Profile::new(outer, holes, tol()).unwrap(),
        f,
        lo,
        hi,
    )
    .unwrap()
    .0
}

/// A centre near `c0` (its coordinates `held` kept) whose `dist2` lies
/// within `[(want + lo)^2, (want + hi)^2]` for the band about `gap` (a
/// factor of two either way), over its coordinates' last bits.
fn search(
    c0: [f64; 3],
    held: [bool; 3],
    dist2: &dyn Fn([f64; 3]) -> R,
    want: f64,
    gap: f64,
) -> [f64; 3] {
    let step = |x: f64, k: i64| f64::from_bits((x.to_bits() as i64 + k) as u64);
    let (lo, hi) = if gap > 0.0 {
        (gap / 2.0, gap * 2.0)
    } else {
        (gap * 2.0, gap / 2.0)
    };
    let (a, b) = (q(want) + q(lo), q(want) + q(hi));
    let (a, b) = (&a * &a, &b * &b);
    let span = |k: i64, i: usize| if held[i] { 0..=0 } else { -k..=k };
    for k in 0..=24i64 {
        for i in span(k, 0) {
            for j in span(k, 1) {
                for l in span(k, 2) {
                    if i.abs().max(j.abs()).max(l.abs()) != k {
                        continue;
                    }
                    let c = [step(c0[0], i), step(c0[1], j), step(c0[2], l)];
                    let d2 = dist2(c);
                    if a <= d2 && d2 <= b {
                        return c;
                    }
                }
            }
        }
    }
    panic!("no centre near {c0:?} for a gap of {gap:e}");
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
/// identities.
fn evaluates(what: &str, a: &Solid, b: &Solid, fused: usize) {
    let f = a.fuse(OperationId(3), b).unwrap().0;
    let c = a.cut(OperationId(4), b).unwrap().0;
    let m = a.common(OperationId(5), b).unwrap().0;
    assert_eq!(f.len(), fused, "{what}");
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let close = |x: f64, y: f64| (x - y).abs() <= 1e-9 * (va + vb);
    assert!(close(volume(&f), va + vb - volume(&m)), "{what} fuse");
    assert!(close(volume(&c), va - volume(&m)), "{what} cut");
}

/// Fuse and both cuts evaluate (a shallow crossing, whose common lens the
/// validator may refuse), the fuse one solid: `fuse = cut + |b|` and `cut
/// back = |b| - (|a| - cut)`.
fn crossing_evaluates(what: &str, a: &Solid, b: &Solid) {
    let f = a.fuse(OperationId(3), b).unwrap().0;
    let c = volume(&a.cut(OperationId(4), b).unwrap().0);
    let back = volume(&b.cut(OperationId(4), a).unwrap().0);
    assert_eq!(f.len(), 1, "{what}");
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let close = |x: f64, y: f64| (x - y).abs() <= 1e-9 * (va + vb);
    assert!(close(volume(&f), c + vb), "{what} fuse");
    assert!(close(back, vb - (va - c)), "{what} cut back");
}

/// A sphere within the resolution of tangency to a cylinder or a cone face
/// is `Degenerate`, as at a plane: a ball missing a rod's wall by `1e-12`
/// to `1e-9` of its radius was fused with it into two solids, a ball inside
/// a rod or a block's hole was cut from it behind a wall thinner than the
/// resolution, a ball by a frustum's wall alike, and one crossing a rod's
/// wall that deep evaluated, was refused by the validator, or ran the
/// validator's integrals for hours. Each ball's distance from the axis (or
/// the frustum's generatrix) is placed by an exact search in the frame's
/// local coordinates at `R +- r (1 + d)`, `d` from `1e-15` to `1e-9` either
/// side, in a level and a turned frame: every operation is refused. At
/// `1e-6` (beyond the resolution) every one evaluates, and so do crossings
/// that deep (fuse and cuts; a minute each before); a gap inside both
/// inputs (a block's cavity inside a rod by its wall) evaluates.
#[test]
fn a_ball_within_the_resolution_of_a_cylinder_or_cone_is_degenerate() {
    let within = [1e-15, -1e-15, 1e-12, -1e-12, 1e-9, -1e-9];
    let r = 0.75;
    for (name, f) in frames() {
        let rod = |radius: f64| {
            Solid::cylinder_with(OperationId(1), f, radius, -3.0, 3.0, tol())
                .unwrap()
                .0
        };
        let radial = |c: [f64; 3]| {
            let l = local(&f, c);
            &l[0] * &l[0] + &l[1] * &l[1]
        };
        // A ball by the axis at `R + side r (1 + d)`.
        let by = |big: f64, side: f64, d: f64| {
            let c0 = world(&f, big + side * r * (1.0 + d), 0.0, 0.25);
            let c = search(c0, [false; 3], &radial, big + side * r, side * r * d);
            ball(c, r)
        };
        let (thin, thick) = (rod(1.0), rod(2.0));
        let hole = Boundary::circle(Point2::new(0.0, 0.0), 2.0, tol()).unwrap();
        let holed = prism(f, square(4.0), vec![hole], -3.0, 3.0, 1);
        // A frustum, radius 2 at 0 and 1 at 3: a ball of radius `rb` by
        // its wall's point at height 1.5, along `(3, 1) / sqrt(10)`.
        let frustum = Solid::cone_with(OperationId(1), f, 2.0, 1.0, 3.0, tol())
            .unwrap()
            .0;
        let cone = |side: f64, rb: f64, d: f64| {
            let k = 10f64.sqrt();
            let s = side * rb * (1.0 + d);
            let c0 = world(&f, 1.5 + s * 3.0 / k, 0.0, 1.5 + s / k);
            let dist2 = |c: [f64; 3]| {
                let l = local(&f, c);
                let e = (&l[0] - q(2.0)) * q(3.0) + &l[2];
                &e * &e / q(10.0)
            };
            ball(search(c0, [false, true, false], &dist2, rb, rb * d), rb)
        };
        for d in within {
            let at = |what: &str| format!("{name} {what} {d:e}");
            refused(&at("rod"), &thin, &by(1.0, 1.0, d), NEAR_QUADRIC);
            refused(&at("in a rod"), &thick, &by(2.0, -1.0, d), NEAR_QUADRIC);
            refused(&at("hole"), &holed, &by(2.0, -1.0, d), NEAR_QUADRIC);
            refused(&at("frustum"), &frustum, &cone(1.0, r, d), NEAR_QUADRIC);
            refused(
                &at("in a frustum"),
                &frustum,
                &cone(-1.0, 0.3, d),
                NEAR_QUADRIC,
            );
        }
        evaluates(&format!("{name} rod"), &thin, &by(1.0, 1.0, 1e-6), 2);
        evaluates(&format!("{name} in a rod"), &thick, &by(2.0, -1.0, 1e-6), 1);
        evaluates(&format!("{name} hole"), &holed, &by(2.0, -1.0, 1e-6), 2);
        evaluates(&format!("{name} frustum"), &frustum, &cone(1.0, r, 1e-6), 2);
        let inner = cone(-1.0, 0.3, 1e-6);
        evaluates(&format!("{name} in a frustum"), &frustum, &inner, 1);
        crossing_evaluates(&format!("{name} rod crossed"), &thin, &by(1.0, 1.0, -1e-6));
        crossing_evaluates(
            &format!("{name} hole crossed"),
            &holed,
            &by(2.0, -1.0, -1e-6),
        );
        // A block's cavity inside the thick rod by its wall: the gap inside
        // both inputs, no contact.
        let block = prism(f, square(6.0), vec![], -5.0, 5.0, 7);
        for d in [1e-12, 1e-9] {
            let hollow = block.cut(OperationId(2), &by(2.0, -1.0, d)).unwrap().0;
            assert_eq!(hollow.len(), 1);
            let fuse = thick.fuse(OperationId(3), &hollow[0]).unwrap().0;
            let common = thick.common(OperationId(5), &hollow[0]).unwrap().0;
            let back = hollow[0].cut(OperationId(4), &thick).unwrap().0;
            let counts = (fuse.len(), common.len(), back.len());
            assert_eq!(counts, (1, 1, 1), "{name} cavity {d:e}");
            let vt = thick.mass_properties().volume;
            let vh = hollow[0].mass_properties().volume;
            let close = |x: f64, y: f64| (x - y).abs() <= 1e-9 * (vt + vh);
            assert!(
                close(volume(&fuse), vt + vh - volume(&common)),
                "{name} {d:e}"
            );
            assert!(close(volume(&back), vh - volume(&common)), "{name} {d:e}");
        }
    }
}

/// An input edge's conic or circle within the resolution of tangency to a
/// plane or a cylinder face of the other input is `Degenerate`, as at a
/// sphere: a dome (a hemisphere turned by 0.3 about `x`) whose rim circle
/// misses a box's top by `1e-15` to `1e-9` of its radius, a rod turned
/// alike whose bottom cap circle misses it, and a rod's cap circle missing
/// another rod's wall (level and turned) were fused with it into two
/// solids, crossings that deep refused by the validator, as a piece thinner
/// than the resolution, or evaluated. Each is placed by an exact search
/// over the dome's or rod's height (the circle's least distance from the
/// plane `|c_z| / sqrt(x_z^2 + y_z^2)` on its exact axes) or in binary64
/// (the cap circle on the rod's wall, gaps far wider than its rounding):
/// every operation is refused; at `1e-6` every one evaluates. A gap inside
/// both inputs evaluates: a block's blind tilted hole (a given result)
/// whose floor's rim misses a slab's floor inside the block by `1e-12`;
/// the same rim over a block's top there is refused, a wall thinner than
/// the resolution.
#[test]
fn an_edge_within_the_resolution_of_a_face_is_degenerate() {
    let half = std::f64::consts::FRAC_PI_2;
    let within = [1e-15, -1e-15, 1e-12, -1e-12, 1e-9, -1e-9];
    let tilt = [0.0, -(0.3f64).sin(), (0.3f64).cos()];
    let lid = Solid::box_at(
        Point3::new(-4.0, -4.0, -4.0),
        Vec3::new(8.0, 8.0, 4.0),
        tol(),
    )
    .unwrap();
    let v = |a: [f64; 3]| a.map(q);
    // A frame turned as `tilt` about `c` whose circle of radius `rr` about
    // its origin has its lowest point `rr (1 + d)` above `z = z0`: its
    // height `|c_z - z0| / sqrt(x_z^2 + y_z^2)` placed by the search.
    let lowest = |rr: f64, d: f64, z0: f64| {
        let g0 = frame([0.25, 0.5, 0.0], tilt, [1.0, 0.0, 0.0]);
        let (x, y) = (v(g0.x().to_array()), v(g0.y().to_array()));
        let k2 = &x[2] * &x[2] + &y[2] * &y[2];
        let dist2 = |c: [f64; 3]| {
            let z = q(c[2]) - q(z0);
            &z * &z / &k2
        };
        let s0 = (1.0 - tilt[2] * tilt[2]).sqrt();
        let c0 = [0.25, 0.5, z0 + s0 * rr * (1.0 + d)];
        frame(
            search(c0, [true, true, false], &dist2, rr, rr * d),
            tilt,
            [1.0, 0.0, 0.0],
        )
    };
    let dome = |d: f64| {
        Solid::sphere_with(OperationId(9), lowest(1.5, d, 0.0), 1.5, 0.0, half, tol())
            .unwrap()
            .0
    };
    let rod = |d: f64| {
        Solid::cylinder_with(OperationId(9), lowest(1.0, d, 0.0), 1.0, 0.0, 3.0, tol())
            .unwrap()
            .0
    };
    for d in within {
        refused(&format!("dome {d:e}"), &lid, &dome(d), NEAR_FACE);
        refused(&format!("rod's cap {d:e}"), &lid, &rod(d), NEAR_FACE);
    }
    evaluates("dome", &lid, &dome(1e-6), 2);
    evaluates("rod's cap", &lid, &rod(1e-6), 2);
    // A rod along the local `x` whose end cap circle (at `u = 0.5`, radius
    // 0.5) misses a rod of radius 1 along the local `w`.
    for (name, f) in frames() {
        let wall = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 3.0, tol())
            .unwrap()
            .0;
        let capped = |d: f64| {
            let y0 = (0.75f64).sqrt() + 0.5 * (1.0 + d);
            let g = Frame3::new(f.point(Point2::new(0.5, y0), 0.0), f.x(), f.y(), tol()).unwrap();
            Solid::cylinder_with(OperationId(9), g, 0.5, 0.0, 3.0, tol())
                .unwrap()
                .0
        };
        for d in [1e-12, -1e-12, 1e-9, -1e-9] {
            refused(
                &format!("{name} cap on a wall {d:e}"),
                &wall,
                &capped(d),
                NEAR_FACE,
            );
        }
        evaluates(&format!("{name} cap on a wall"), &wall, &capped(1e-6), 2);
    }
    // A block with a blind hole: the turned rod's lower part, its floor's
    // rim's lowest point `1e-12` above `z = -2`; a slab over `z = -2`
    // inside the block (the gap inside both) and a block below it.
    let xy = frame([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let block = prism(xy, square(4.0), vec![], -4.0, 0.0, 1);
    let pocket = Solid::cylinder_with(
        OperationId(2),
        lowest(1.0, 1e-12, -2.0),
        1.0,
        0.0,
        6.0,
        tol(),
    )
    .unwrap()
    .0;
    let blind = block.cut(OperationId(6), &pocket).unwrap().0;
    assert_eq!(blind.len(), 1);
    let slab = |lo: f64, hi: f64| prism(xy, square(3.0), vec![], lo, hi, 11);
    let over = slab(-2.0, -1.0);
    let fuse = blind[0].fuse(OperationId(3), &over).unwrap().0;
    let common = blind[0].common(OperationId(5), &over).unwrap().0;
    let cut = blind[0].cut(OperationId(4), &over).unwrap().0;
    assert_eq!((fuse.len(), common.len()), (1, 1));
    let (va, vb) = (
        blind[0].mass_properties().volume,
        over.mass_properties().volume,
    );
    assert!((volume(&fuse) - (va + vb - volume(&common))).abs() <= 1e-9 * (va + vb));
    assert!((volume(&cut) - (va - volume(&common))).abs() <= 1e-9 * (va + vb));
    refused(
        "a blind hole's floor over a block",
        &blind[0],
        &slab(-3.0, -2.0),
        NEAR_FACE,
    );
}

/// A ball within the resolution of an edge on a meeting of two curved
/// faces, a torus's section or a spline curve is `Degenerate` with the
/// edge rule's reason (S9d.1's `edge_near_misses`, its nearest points on
/// such curves found in binary64 at rational parameters, `curved/near.rs`):
/// a ball on the rim of a rod's hole (a rod less a crossing rod), of a
/// rod's dimple (a rod less a ball), of a torus less a block above its
/// equator and of a spline dome prism's top cap was fused with it into two
/// solids. Each ball is placed in binary64 on the edge's normals' bisector
/// at `r (1 + d)` from the edge's point, gaps far wider than the placement's
/// rounding: `+-1e-12` and `+-1e-9` refused, `1e-6` evaluating.
#[test]
fn a_ball_within_the_resolution_of_a_meeting_or_spline_edge_is_degenerate() {
    let pi = std::f64::consts::PI;
    for (name, f) in frames() {
        let rod = Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 3.0, tol())
            .unwrap()
            .0;
        // The rim's point `(1, 0, 0.5)` with its faces' normals `(1, 0, 0)`
        // and `(0, 0, -1)`; a ball of radius 0.2 on their bisector.
        let rim_ball = |d: f64| {
            let k = 0.2 * (1.0 + d) * std::f64::consts::FRAC_1_SQRT_2;
            ball(world(&f, 1.0 + k, 0.0, 0.5 - k), 0.2)
        };
        let across =
            Frame3::new(f.point(Point2::new(-3.0, 0.0), 0.0), f.x(), f.y(), tol()).unwrap();
        let bar = Solid::cylinder_with(OperationId(2), across, 0.5, 0.0, 6.0, tol())
            .unwrap()
            .0;
        // The biting ball's axis off the rod's: its poles off the wall.
        let bite = ball_in(
            Frame3::new(f.point(Point2::new(1.0, 0.0), 0.0), f.x(), f.y(), tol()).unwrap(),
            0.5,
        );
        // A spline dome `v = u (4 - u) / 2` over `[0, 4]`, heights 0 to 2:
        // its top cap's apex `(2, 2, 2)`, normals `(0, 1, 0)` and `(0, 0, 1)`.
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
        let spline = prism(f, outer, vec![], 0.0, 2.0, 1);
        let apex_ball = |d: f64| {
            let k = 0.25 * (1.0 + d) * std::f64::consts::FRAC_1_SQRT_2;
            ball(world(&f, 2.0, 2.0 + k, 2.0 + k), 0.25)
        };
        // A torus (3, 1) less a block above `w = 0.5`: its outer section
        // circle's point `(3 + sqrt(0.75), 0, 0.5)`, normals `(sqrt(0.75),
        // 0, 0.5)` and `(0, 0, 1)`.
        let torus = Solid::torus_with(OperationId(1), f, 3.0, 1.0, -pi, pi, 2.0 * pi, tol())
            .unwrap()
            .0;
        let above = prism(f, square(5.0), vec![], 0.5, 2.0, 2);
        let ring_ball = |d: f64| {
            let s = (0.75f64).sqrt();
            let len = (0.75f64 + 2.25).sqrt();
            let rr = 0.25 * (1.0 + d);
            ball(
                world(&f, 3.0 + s + rr * s / len, 0.0, 0.5 + rr * 1.5 / len),
                0.25,
            )
        };
        let cases: [(&str, Solid, &Placed); 4] = [
            (
                "hole's rim",
                rod.cut(OperationId(6), &bar).unwrap().0.remove(0),
                &rim_ball,
            ),
            (
                "dimple's rim",
                rod.cut(OperationId(6), &bite).unwrap().0.remove(0),
                &rim_ball,
            ),
            ("spline apex", spline, &apex_ball),
            (
                "torus section",
                torus.cut(OperationId(6), &above).unwrap().0.remove(0),
                &ring_ball,
            ),
        ];
        for (what, body, at) in &cases {
            for d in [1e-12, -1e-12, 1e-9, -1e-9] {
                refused(&format!("{name} {what} {d:e}"), body, &at(d), NEAR_EDGE);
            }
            evaluates(&format!("{name} {what}"), body, &at(1e-6), 2);
        }
    }
}
