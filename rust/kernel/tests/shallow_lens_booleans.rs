//! The shallow lens (HANDOFF item 5 (d), REVIEW_NOTES.md's "The shallow
//! lens refined"): a ball of radius 0.75 crossing a rod's wall (radius 1,
//! heights -3 to 3) by `1e-5` to `1.4e-7` of its radius (the last `1.05e-7`
//! deep, just beyond the resolution), and the controls `0.05` and `0.25`,
//! in the level frame and the boolean target's turned one, against
//! `tools/shallow_lens_reference.py`'s 30-digit references
//! (`fixtures/shallow-lens-cases.tsv`, from
//! `tools/generate_shallow_lens_fixtures.py`). Every operation evaluates,
//! one solid each; every volume and area enclosure holds its reference and
//! is at most `1e-9` wide (of the value, or absolutely below one); and the
//! results are consistent: fuse and common enclose `|A| + |B|` together,
//! each cut and the common the cut body's own. The fixture's frames are
//! the kernel's stored ones, bit for bit, and so is the rod's wall's origin
//! (`o - 3 n`), which the reference integrates on. The other procedural
//! curves' shallow loops (a ball by a frustum's wall, two rods crossing,
//! balls, a rod and a slab by a torus) evaluate consistently, and a lens on
//! a spline wall, which the narrow tier does not reach, is a
//! `ComputationLimit`.
use num_bigint::BigInt;
use num_rational::BigRational as R;
use rusty_occt::identity::OperationId;
use rusty_occt::topology::{SplineSpan, Surface};
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance,
    Vec3,
};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};

const FIXTURE: &str = include_str!("../../fixtures/shallow-lens-cases.tsv");

fn tol() -> Tolerance {
    Tolerance::default()
}

fn hex(word: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(word, 16).unwrap())
}

fn hex3(text: &str) -> [f64; 3] {
    let v: Vec<f64> = text.split(' ').map(hex).collect();
    [v[0], v[1], v[2]]
}

/// A decimal's exact rational.
fn decimal(text: &str) -> R {
    let (mantissa, exponent) = match text.split_once(['e', 'E']) {
        Some((m, e)) => (m, e.parse::<i32>().unwrap()),
        None => (text, 0),
    };
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits: BigInt = format!("{whole}{fraction}").parse().unwrap();
    let shift = exponent - fraction.len() as i32;
    let ten = BigInt::from(10);
    if shift >= 0 {
        R::from_integer(digits * ten.pow(shift as u32))
    } else {
        R::new(digits, ten.pow((-shift) as u32))
    }
}

fn q(x: f64) -> R {
    R::from_float(x).unwrap()
}

/// A frame of the fixture: its construction's origin, normal and hint, and
/// its stored axes and the rod's wall's origin as the reference computes
/// them.
struct FrameRow {
    origin: [f64; 3],
    axes: [[f64; 3]; 3],
    bottom: [f64; 3],
}

struct Case {
    name: String,
    frame: String,
    centre: [f64; 3],
    refs: BTreeMap<String, R>,
}

const KEYS: [&str; 8] = [
    "lens_volume",
    "lens_area",
    "lens_wall_area",
    "lens_sphere_area",
    "rod_volume",
    "rod_area",
    "ball_volume",
    "ball_area",
];

fn fixture() -> (BTreeMap<String, FrameRow>, Vec<Case>) {
    let mut frames = BTreeMap::new();
    let mut cases = Vec::new();
    for line in FIXTURE.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "frame" => {
                frames.insert(
                    f[1].to_string(),
                    FrameRow {
                        origin: hex3(f[2]),
                        axes: [hex3(f[3]), hex3(f[4]), hex3(f[5])],
                        bottom: hex3(f[6]),
                    },
                );
            }
            "case" => cases.push(Case {
                name: f[1].to_string(),
                frame: f[2].to_string(),
                centre: hex3(f[4]),
                refs: KEYS
                    .iter()
                    .zip(&f[5..])
                    .map(|(k, v)| (k.to_string(), decimal(v)))
                    .collect(),
            }),
            other => panic!("unknown row {other}"),
        }
    }
    (frames, cases)
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

/// The fixture's frames' constructions (as the generator's `FRAMES`).
fn construction(name: &str) -> Frame3 {
    match name {
        "level" => frame([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        "turned" => frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
        other => panic!("unknown frame {other}"),
    }
}

fn rod(f: Frame3) -> Solid {
    Solid::cylinder_with(OperationId(1), f, 1.0, -3.0, 3.0, tol())
        .unwrap()
        .0
}

fn ball(c: [f64; 3]) -> Solid {
    let half = std::f64::consts::FRAC_PI_2;
    let f = frame(c, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    Solid::sphere_with(OperationId(9), f, 0.75, -half, half, tol())
        .unwrap()
        .0
}

/// A body's volume and area enclosures.
fn measures(s: &Solid) -> [[R; 2]; 2] {
    let e = s.topology().mass_enclosure().expect("mass properties");
    [
        [q(e.volume[0]), q(e.volume[1])],
        [q(e.surface_area[0]), q(e.surface_area[1])],
    ]
}

/// At most `1e-9` wide, of the value or absolutely below one.
fn narrow(e: &[R; 2]) -> bool {
    let one = R::from_integer(BigInt::from(1));
    let abs = |x: &R| {
        if *x < R::from_integer(BigInt::from(0)) {
            -x.clone()
        } else {
            x.clone()
        }
    };
    let big = abs(&e[0]).max(abs(&e[1])).max(one);
    &e[1] - &e[0] <= big * R::new(BigInt::from(1), BigInt::from(1_000_000_000))
}

/// The four operations' results, each one solid's measures, by name.
type Results = BTreeMap<&'static str, Measures>;

/// A body's volume and area enclosures.
type Measures = [[R; 2]; 2];

fn run(rod: &Solid, ball: &Solid) -> Result<Results, String> {
    let mut out = BTreeMap::new();
    for (op, result) in [
        ("fuse", rod.fuse(OperationId(3), ball)),
        ("cut", rod.cut(OperationId(4), ball)),
        ("cut back", ball.cut(OperationId(4), rod)),
        ("common", rod.common(OperationId(5), ball)),
    ] {
        let solids = result.map_err(|e| format!("{op}: {e:?}"))?.0;
        if solids.len() != 1 {
            return Err(format!("{op}: {} solids", solids.len()));
        }
        out.insert(op, measures(&solids[0]));
    }
    Ok(out)
}

/// `f` over the cases on a few threads, in their order.
fn each<C: Sync, T: Send>(cases: &[C], f: impl Fn(&C) -> T + Sync) -> Vec<T> {
    let next = AtomicUsize::new(0);
    let workers = std::thread::available_parallelism()
        .map_or(2, |n| n.get())
        .clamp(1, 6);
    let mut out: Vec<(usize, T)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    let mut mine = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::SeqCst);
                        let Some(case) = cases.get(i) else { break };
                        mine.push((i, f(case)));
                    }
                    mine
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("a worker"))
            .collect()
    });
    out.sort_by_key(|x| x.0);
    out.into_iter().map(|x| x.1).collect()
}

/// The fixture's frames are the kernel's stored ones, and the rod's wall's
/// frame origin is the reference's `o - 3 n`.
#[test]
fn the_fixtures_frames_are_the_kernels() {
    let (frames, _) = fixture();
    for (name, row) in &frames {
        let f = construction(name);
        let stored = [f.x(), f.y(), f.normal()].map(|v| v.to_array().map(f64::to_bits));
        assert_eq!(f.origin().to_array(), row.origin, "{name}");
        assert_eq!(stored, row.axes.map(|v| v.map(f64::to_bits)), "{name}");
        let wall = rod(f)
            .topology()
            .faces()
            .iter()
            .find_map(|face| match &face.surface {
                Surface::Cylinder { frame, .. } => Some(*frame),
                _ => None,
            })
            .expect("a wall");
        assert_eq!(
            wall.origin().to_array().map(f64::to_bits),
            row.bottom.map(f64::to_bits),
            "{name}"
        );
        let axes = [wall.x(), wall.y(), wall.normal()].map(|v| v.to_array().map(f64::to_bits));
        assert_eq!(axes, row.axes.map(|v| v.map(f64::to_bits)), "{name}");
    }
}

/// Every operation of every case evaluates to one solid whose volume and
/// area enclosures hold the reference's and are at most `1e-9` wide; the
/// bodies' own too; and fuse with common, each cut with common, enclose
/// the bodies' sums together. Before the narrow tier the commons of the
/// lenses `1e-6` and shallower were refused (`uncertified_shell_orientation`
/// or their mass properties unrepresentable), and the deeper ones' volumes
/// enclosed a hundredth wide.
#[test]
fn shallow_lenses_evaluate_within_their_references() {
    let (_, cases) = fixture();
    let outcomes = each(&cases, |case| {
        let r = rod(construction(&case.frame));
        let b = ball(case.centre);
        (measures(&r), measures(&b), run(&r, &b))
    });
    let mut failures = Vec::new();
    for (case, (a, b, results)) in cases.iter().zip(outcomes) {
        let refs = &case.refs;
        let results = match results {
            Ok(results) => results,
            Err(e) => {
                failures.push(format!("{}: {e}", case.name));
                continue;
            }
        };
        let (lens, wall, sphere) = (
            &refs["lens_volume"],
            &refs["lens_wall_area"],
            &refs["lens_sphere_area"],
        );
        let (rod_v, rod_a) = (&refs["rod_volume"], &refs["rod_area"]);
        let (ball_v, ball_a) = (&refs["ball_volume"], &refs["ball_area"]);
        let mut want: Vec<(&str, &Measures, [R; 2])> = vec![
            ("rod", &a, [rod_v.clone(), rod_a.clone()]),
            ("ball", &b, [ball_v.clone(), ball_a.clone()]),
        ];
        let expected: [(&str, [R; 2]); 4] = [
            (
                "fuse",
                [rod_v + ball_v - lens, rod_a - wall + (ball_a - sphere)],
            ),
            ("cut", [rod_v - lens, rod_a - wall + sphere]),
            ("cut back", [ball_v - lens, ball_a - sphere + wall]),
            ("common", [lens.clone(), wall + sphere]),
        ];
        for (op, values) in expected {
            want.push((op, &results[op], values));
        }
        for (what, got, values) in want {
            for (k, (e, v)) in got.iter().zip(values.iter()).enumerate() {
                let measure = ["volume", "area"][k];
                if !(e[0] <= *v && *v <= e[1]) {
                    failures.push(format!(
                        "{} {what} {measure}: [{}, {}] misses {}",
                        case.name,
                        rat(&e[0]),
                        rat(&e[1]),
                        rat(v)
                    ));
                }
                if !narrow(e) {
                    failures.push(format!(
                        "{} {what} {measure}: [{}, {}] wide",
                        case.name,
                        rat(&e[0]),
                        rat(&e[1])
                    ));
                }
            }
        }
        // Consistency: the enclosures' sums meet.
        let sum = |x: &[R; 2], y: &[R; 2]| [&x[0] + &y[0], &x[1] + &y[1]];
        let meet = |x: &[R; 2], y: &[R; 2]| x[0] <= y[1] && y[0] <= x[1];
        let common = &results["common"];
        for k in 0..2 {
            let bodies = sum(&a[k], &b[k]);
            let pairs = [
                (
                    "fuse + common = A + B",
                    sum(&results["fuse"][k], &common[k]),
                    bodies,
                ),
                (
                    "cut + common = A",
                    sum(&results["cut"][k], &common[k]),
                    a[k].clone(),
                ),
                (
                    "cut back + common = B",
                    sum(&results["cut back"][k], &common[k]),
                    b[k].clone(),
                ),
            ];
            // Areas only for the fuse: a cut's faces and the common's share
            // the lens's surface, counted on both.
            for (what, x, y) in pairs.iter().take(if k == 0 { 3 } else { 1 }) {
                if !meet(x, y) {
                    failures.push(format!("{} {what}: apart", case.name));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

fn turned() -> Frame3 {
    construction("turned")
}

fn level() -> Frame3 {
    construction("level")
}

fn world(f: &Frame3, u: f64, v: f64, w: f64) -> [f64; 3] {
    f.point(Point2::new(u, v), w).to_array()
}

fn direction(f: &Frame3, u: f64, v: f64, w: f64) -> [f64; 3] {
    let (a, b) = (world(f, u, v, w), world(f, 0.0, 0.0, 0.0));
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn ball_of(c: [f64; 3], r: f64) -> Solid {
    let half = std::f64::consts::FRAC_PI_2;
    let f = frame(c, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    Solid::sphere_with(OperationId(9), f, r, -half, half, tol())
        .unwrap()
        .0
}

fn torus(f: Frame3) -> Solid {
    let turn = std::f64::consts::TAU;
    Solid::torus_with(OperationId(1), f, 2.5, 1.0, 0.0, turn, turn, tol())
        .unwrap()
        .0
}

/// A shallow crossing of another procedural curve: its name and inputs.
type Crossing = (&'static str, fn(f64) -> (Solid, Solid), f64);

/// The other procedural curves' shallow loops (REVIEW_NOTES.md's "The
/// shallow lens refined" measures them): a ball of radius 0.75 by a
/// frustum's wall (a meeting over the cone's angle), a rod of radius 0.5
/// crossing a rod's wall at right angles (two cylinders' meeting), balls
/// at a torus's outer equator and its top, a rod of radius 0.5 by its
/// outer wall (torus meetings) and a slab's face across its equator (a
/// spiric section), by `1.5e-7` or `1e-6` of the smaller radius: every
/// operation evaluates, the commons (refused before, as an undecided
/// orientation or an unrepresentable mass) within `1e-9` and positive, and
/// fuse and common, each cut and common, enclose the bodies' sums
/// together.
#[test]
fn other_procedural_curves_shallow_loops_evaluate_consistently() {
    let crossings: Vec<Crossing> = vec![
        (
            "ball by a frustum, level",
            |d| frustum_ball(level(), d),
            1.5e-7,
        ),
        (
            "ball by a frustum, turned",
            |d| frustum_ball(turned(), d),
            1.5e-7,
        ),
        (
            "rod across a rod, level",
            |d| {
                let f = level();
                let off = 1.0 + 0.5 * (1.0 - d);
                let g = frame(
                    world(&f, off, -2.0, 0.25),
                    direction(&f, 0.0, 1.0, 0.0),
                    direction(&f, 1.0, 0.0, 0.0),
                );
                let small = Solid::cylinder_with(OperationId(2), g, 0.5, 0.0, 4.0, tol())
                    .unwrap()
                    .0;
                (rod(f), small)
            },
            1.5e-7,
        ),
        (
            "ball at a torus's equator, turned",
            |d| {
                let f = turned();
                (
                    torus(f),
                    ball_of(world(&f, 3.5 + 0.75 * (1.0 - d), 0.0, 0.0), 0.75),
                )
            },
            1.5e-7,
        ),
        (
            "ball on a torus's top, turned",
            |d| {
                let f = turned();
                let c = world(&f, 1.5, 2.0, 1.0 + 0.75 * (1.0 - d));
                (torus(f), ball_of(c, 0.75))
            },
            1.5e-7,
        ),
        (
            "rod by a torus, level",
            |d| {
                let f = level();
                let g = frame(
                    world(&f, 3.5 + 0.5 * (1.0 - d), 0.0, -2.0),
                    direction(&f, 0.0, 0.0, 1.0),
                    direction(&f, 1.0, 0.0, 0.0),
                );
                let rod = Solid::cylinder_with(OperationId(2), g, 0.5, 0.0, 4.0, tol())
                    .unwrap()
                    .0;
                (torus(f), rod)
            },
            1e-6,
        ),
        (
            "slab across a torus's equator, level",
            |d| {
                let f = level();
                let g = frame(
                    world(&f, 3.5 - d, 0.0, 0.0),
                    direction(&f, 1.0, 0.0, 0.0),
                    direction(&f, 0.0, 1.0, 0.0),
                );
                (torus(f), slab(g, 2.0, 0.0, 1.0))
            },
            1.5e-7,
        ),
    ];
    let outcomes = each(&crossings, |(name, make, d)| {
        let (a, b) = make(*d);
        let volume = |s: &Solid| measures(s)[0].clone();
        let one = |op: &str, r: Result<(Vec<Solid>, rusty_occt::history::History), Error>| {
            let solids = r.map_err(|e| format!("{name} {op}: {e:?}"))?.0;
            let mut total = [
                R::from_integer(BigInt::from(0)),
                R::from_integer(BigInt::from(0)),
            ];
            for s in &solids {
                let v = volume(s);
                total = [&total[0] + &v[0], &total[1] + &v[1]];
            }
            Ok::<_, String>((solids.len(), total))
        };
        let fuse = one("fuse", a.fuse(OperationId(3), &b))?;
        let cut = one("cut", a.cut(OperationId(4), &b))?;
        let back = one("cut back", b.cut(OperationId(4), &a))?;
        let common = one("common", a.common(OperationId(5), &b))?;
        Ok::<_, String>((volume(&a), volume(&b), fuse, cut, back, common))
    });
    let mut failures = Vec::new();
    let zero = R::from_integer(BigInt::from(0));
    for ((name, _, d), outcome) in crossings.iter().zip(outcomes) {
        let (a, b, fuse, cut, back, common) = match outcome {
            Ok(x) => x,
            Err(e) => {
                failures.push(e);
                continue;
            }
        };
        let sum = |x: &[R; 2], y: &[R; 2]| [&x[0] + &y[0], &x[1] + &y[1]];
        let meet = |x: &[R; 2], y: &[R; 2]| x[0] <= y[1] && y[0] <= x[1];
        if common.0 != 1 || common.1[0] <= zero || !narrow(&common.1) {
            failures.push(format!(
                "{name} {d:e}: common {} solids [{}, {}]",
                common.0,
                rat(&common.1[0]),
                rat(&common.1[1])
            ));
        }
        for (what, x, y) in [
            (
                "fuse + common = A + B",
                sum(&fuse.1, &common.1),
                sum(&a, &b),
            ),
            ("cut + common = A", sum(&cut.1, &common.1), a.clone()),
            ("cut back + common = B", sum(&back.1, &common.1), b.clone()),
        ] {
            if !meet(&x, &y) {
                failures.push(format!("{name} {d:e}: {what} apart"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

fn frustum_ball(f: Frame3, d: f64) -> (Solid, Solid) {
    // Radius 2 at 0 and 1 at 3: the ball by the wall's point at height 1.5,
    // along `(3, 1) / sqrt(10)`.
    let frustum = Solid::cone_with(OperationId(1), f, 2.0, 1.0, 3.0, tol())
        .unwrap()
        .0;
    let k = 10f64.sqrt();
    let s = 0.75 * (1.0 - d);
    let c = world(&f, 1.5 + s * 3.0 / k, 0.0, 1.5 + s / k);
    (frustum, ball_of(c, 0.75))
}

fn slab(f: Frame3, s: f64, lo: f64, hi: f64) -> Solid {
    let square = Boundary::polygon(
        [(-s, -s), (s, -s), (s, s), (-s, s)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec(),
        tol(),
    )
    .unwrap();
    Solid::extrude_with(
        OperationId(2),
        Profile::new(square, vec![], tol()).unwrap(),
        f,
        lo,
        hi,
    )
    .unwrap()
    .0
}

/// A lens on a spline wall stays refused, now a `ComputationLimit`: a rod
/// of radius 0.5 along `x` crossing the ridge of a dome prism (the profile
/// `y = x (4 - x) / 2` over `[0, 4]`, a quadratic spline, extruded along
/// `z`) by `1.4e-7` of its radius. Its meeting runs on the spline wall
/// (S9f.2b), whose fluxes the narrow tier does not narrow: the common's
/// orientation stays undecided (`uncertified_shell_orientation` before),
/// the fuse and both cuts evaluate.
#[test]
fn a_shallow_lens_on_a_spline_wall_is_a_computation_limit() {
    let poles = [(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)].map(|(x, y)| Point2::new(x, y));
    let curve =
        BSplineCurve2::new(2, poles.to_vec(), None, vec![0.0, 1.0], vec![3, 3]).expect("a spline");
    let outer = Boundary::path(
        vec![Point2::new(0.0, 0.0), Point2::new(4.0, 0.0)],
        vec![Segment::Line, Segment::Spline(SplineSpan::whole(curve))],
        tol(),
    )
    .expect("a dome");
    let dome = Solid::extrude_with(
        OperationId(1),
        Profile::new(outer, vec![], tol()).expect("a profile"),
        Frame3::xy(),
        -2.0,
        2.0,
    )
    .expect("a prism")
    .0;
    let r = 0.5;
    let f = frame(
        [0.5, 2.0 + r * (1.0 - 1.4e-7), 0.25],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    );
    let rod = Solid::cylinder_with(OperationId(2), f, r, 0.0, 3.0, tol())
        .unwrap()
        .0;
    assert!(matches!(
        dome.common(OperationId(5), &rod),
        Err(Error::ComputationLimit(m))
            if m == "a shell's orientation the validator's integrals leave undecided"
    ));
    for (op, out) in [
        ("fuse", dome.fuse(OperationId(3), &rod)),
        ("cut", dome.cut(OperationId(4), &rod)),
        ("cut back", rod.cut(OperationId(4), &dome)),
    ] {
        assert_eq!(out.map(|x| x.0.len()).ok(), Some(1), "{op}");
    }
}

/// A rational's twenty significant digits, for messages.
fn rat(x: &R) -> String {
    let zero = BigInt::from(0);
    if *x.numer() == zero {
        return "0".into();
    }
    let ten = BigInt::from(10);
    let (sign, mut v) = if *x.numer() < zero {
        ("-", -x.clone())
    } else {
        ("", x.clone())
    };
    let mut exponent = 0i32;
    let (lo, hi) = (R::from_integer(ten.pow(19)), R::from_integer(ten.pow(20)));
    while v < lo {
        v *= R::from_integer(ten.clone());
        exponent -= 1;
    }
    while v >= hi {
        v /= R::from_integer(ten.clone());
        exponent += 1;
    }
    let digits = v.to_integer().to_string();
    format!("{sign}{}.{}e{}", &digits[..1], &digits[1..], exponent + 19)
}
