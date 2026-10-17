//! A whole torus of radii 3 and 1 and a rod about its axis whose wall
//! crosses the torus's inner equator by `g` (radius `2 + g`, coaxial: a
//! meeting near a double root along a whole circle), in a level frame, one
//! turned about its axis and one tilted off the world's axes (the rod's
//! frame built in the torus's, so its axis is the torus's within rounding).
//! In the tilted frame a crossing of `1e-9` ran for over half an hour in
//! the torus meeting's verification (`ToricCrv::verified`), its runs split
//! and their pieces verified again up to 48 times, each verification capped
//! but not their total; the meeting's verifications now draw on one budget
//! (REVIEW_NOTES.md, "A torus meeting's verification budget"). Its
//! crossings beyond the resolution took minutes an operation there in
//! their results' validation (the meeting's jets in the validator's
//! integrals, widening a hundredfold an order by the double root); they
//! take about a second now (REVIEW_NOTES.md, "A tilted torus meeting's
//! validation").
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Frame3, Point2, Point3, Solid, Tolerance, Vec3};

const TANGENCY: &str = "a tangency between the inputs (S9c)";
const COAXIAL: &str = "a torus's coaxial circle left unverified";
const STALLS: &str = "a torus meeting's trace stalls";
const BUDGET: &str = "a torus meeting's verification over its budget";

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
/// tilted off it (the boolean fuzz target's).
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

/// The frame's point at local coordinates `l`.
fn at(f: &Frame3, l: [f64; 3]) -> [f64; 3] {
    f.point(Point2::new(l[0], l[1]), l[2]).to_array()
}

/// The frame's direction of local coordinates `l`.
fn dir(f: &Frame3, l: [f64; 3]) -> [f64; 3] {
    let (x, y, n) = (f.x().to_array(), f.y().to_array(), f.normal().to_array());
    std::array::from_fn(|i| x[i] * l[0] + y[i] * l[1] + n[i] * l[2])
}

/// The whole torus of radii 3 and 1 on `f`.
fn torus(f: Frame3) -> Solid {
    let pi = std::f64::consts::PI;
    Solid::torus_with(OperationId(1), f, 3.0, 1.0, -pi, pi, 2.0 * pi, tol())
        .unwrap()
        .0
}

/// The rod of radius `2 + g` about `f`'s axis from `-1.5` to `1.5` (its
/// frame at the local point `(0, 0, -1.5)`, with `f`'s normal and `x`).
fn rod(f: &Frame3, g: f64) -> Solid {
    let rf = frame(
        at(f, [0.0, 0.0, -1.5]),
        dir(f, [0.0, 0.0, 1.0]),
        dir(f, [1.0, 0.0, 0.0]),
    );
    Solid::cylinder_with(OperationId(2), rf, 2.0 + g, 0.0, 3.0, tol())
        .unwrap()
        .0
}

fn volume(s: &[Solid]) -> f64 {
    s.iter().map(|x| x.mass_properties().volume).sum()
}

/// The four operations' outcomes: fuse, cut, the cut back and common.
fn outcomes(t: &Solid, r: &Solid) -> [Result<Vec<Solid>, Error>; 4] {
    [
        t.fuse(OperationId(3), r),
        t.cut(OperationId(4), r),
        r.cut(OperationId(4), t),
        t.common(OperationId(5), r),
    ]
    .map(|out| out.map(|x| x.0))
}

/// Every operation refused, with one of `why`.
fn refused(what: &str, t: &Solid, r: &Solid, why: &[Error]) {
    for (op, out) in ["fuse", "cut", "cut back", "common"]
        .iter()
        .zip(outcomes(t, r))
    {
        match out {
            Err(e) => assert!(why.contains(&e), "{what} {op}: {e:?}"),
            Ok(s) => panic!("{what} {op}: {} solids, not {why:?}", s.len()),
        }
    }
}

/// The crossings beyond the resolution, `1e-6`, `8e-7` and `5e-7`, and
/// their reference volumes (mpmath, 25 digits, rounded once): the rod's
/// `3 pi rho^2` at `rho = 2 + g` in binary64, and its common with the torus
/// by Pappus over the tube's segment, `2 pi` times the integral of `2 x
/// sqrt(1 - (x - 3)^2)` over `[2, rho]` (also in closed form).
const BEYOND: [(f64, f64, f64); 3] = [
    (1e-6, 37.699_149_542_198_79, 2.369_537_922_945_086e-8),
    (8e-7, 37.699_142_002_373_02, 1.695_503_267_562_629_4e-8),
    (5e-7, 37.699_130_692_635_8, 8.377_581_039_647_674e-9),
];

/// The torus's volume, `6 pi^2`.
const TORUS: f64 = 59.217_626_406_536_15;

/// One solid's certified volume, at most `1e-9` wide (relative to its size,
/// at least one), holding the reference within `slack` (its rounding).
fn certified(what: &str, s: &Solid, reference: f64, slack: f64) {
    let [lo, hi] = s
        .topology()
        .mass_enclosure()
        .unwrap_or_else(|| panic!("{what}: its mass undecided"))
        .volume;
    let size = lo.abs().max(hi.abs()).max(1.0);
    assert!(hi - lo <= 1e-9 * size, "{what}: wide [{lo}, {hi}]");
    assert!(
        lo - slack <= reference && reference <= hi + slack,
        "{what}: [{lo}, {hi}] against {reference}"
    );
}

/// A crossing beyond the resolution (`BEYOND`'s `(g, rod, common)`): each
/// operation one solid, its certified volume holding the reference's (the
/// sums' within their rounding), the fuse and the common together the
/// inputs, the cut and the common the torus, the cut back and the common
/// the rod.
fn evaluates(what: &str, t: &Solid, r: &Solid, (_, rod, common): (f64, f64, f64)) {
    let (vt, vr) = (t.mass_properties().volume, r.mass_properties().volume);
    let close = |x: f64, y: f64| (x - y).abs() <= 1e-12 * (TORUS + rod);
    assert!(close(vt, TORUS) && close(vr, rod), "{what}: {vt} {vr}");
    let rounding = |x: f64| x.abs().max(1.0) * f64::EPSILON;
    certified(&format!("{what} torus"), t, TORUS, rounding(TORUS));
    certified(&format!("{what} rod"), r, rod, rounding(rod));
    let [fuse, cut, back, both] = outcomes(t, r).map(|out| out.expect(what));
    let sums = 4.0 * rounding(TORUS + rod);
    for (op, s, reference, slack) in [
        ("fuse", &fuse, TORUS + rod - common, sums),
        ("cut", &cut, TORUS - common, sums),
        ("cut back", &back, rod - common, sums),
        ("common", &both, common, rounding(common)),
    ] {
        assert_eq!(s.len(), 1, "{what} {op}");
        assert!(close(volume(s), reference), "{what} {op} {}", volume(s));
        certified(&format!("{what} {op}"), &s[0], reference, slack);
    }
    let [vf, vc, vb, vm] = [&fuse, &cut, &back, &both].map(|s| volume(s));
    assert!(close(vf + vm, vt + vr), "{what} fuse and common");
    assert!(close(vc + vm, vt), "{what} cut and common");
    assert!(close(vb + vm, vr), "{what} cut back and common");
}

/// The rod tangent along the inner equator (`g = 0`) and crossing it by
/// `1e-12`, `1e-9` and `5e-8`. Coaxial exactly (level, turned): the
/// tangency `Degenerate`, the crossings within the resolution circles the
/// box verification leaves unverified (`ComputationLimit`). Coaxial within
/// rounding (tilted), each a `ComputationLimit`: the tangency and `1e-12` a
/// stalled trace; `1e-9` the meeting over its verification budget where
/// the platform's `hypot` is macOS's (about a second a meeting in release;
/// over half an hour before), a stalled trace under glibc's correctly
/// rounded one (the rod's frame normalized again); `5e-8` over the budget
/// under both. The crossings beyond the resolution are
/// `rod_crossings_beyond_the_resolution_evaluate_*`.
#[test]
fn a_rod_crossing_a_torus_inner_equator_is_decided_quickly() {
    let limit = |why: &'static str| Error::ComputationLimit(why);
    for (name, f) in frames() {
        let t = torus(f);
        let coaxial = name != "tilted";
        let at = |g: f64| format!("{name} {g:e}");
        let within: [(f64, Vec<Error>); 4] = if coaxial {
            [
                (0.0, vec![Error::Degenerate(TANGENCY)]),
                (1e-12, vec![limit(COAXIAL)]),
                (1e-9, vec![limit(COAXIAL)]),
                (5e-8, vec![limit(COAXIAL)]),
            ]
        } else {
            [
                (0.0, vec![limit(STALLS)]),
                (1e-12, vec![limit(STALLS)]),
                (1e-9, vec![limit(BUDGET), limit(STALLS)]),
                (5e-8, vec![limit(BUDGET)]),
            ]
        };
        for (g, why) in within {
            refused(&at(g), &t, &rod(&f, g), &why);
        }
    }
}

/// The crossings beyond the resolution in one frame: each operation one
/// solid of the reference's volumes.
fn beyond(frame: usize) {
    let (name, f) = frames()[frame];
    let t = torus(f);
    for case in BEYOND {
        evaluates(&format!("{name} {:e}", case.0), &t, &rod(&f, case.0), case);
    }
}

#[test]
fn rod_crossings_beyond_the_resolution_evaluate_level() {
    beyond(0);
}

#[test]
fn rod_crossings_beyond_the_resolution_evaluate_turned() {
    beyond(1);
}

/// The tilted frame's: coaxial within rounding, its meeting traced and
/// verified in pieces, two closed curves near a double root. Validating its
/// results' faces took minutes an operation (100 to 274 s at `1e-6`, some
/// 200 s a fuse at `5e-7` and `8e-7`): the meeting's jets in the
/// validator's integrals, from the other surface's functionals' squares
/// whose terms in the torus's angle cancel only exactly, widened a
/// hundredfold an order by the division by the small slope there, so each
/// integral took some 65,000 pieces; from `G`'s harmonic forms, about a
/// second (REVIEW_NOTES.md, "A tilted torus meeting's validation").
#[test]
fn rod_crossings_beyond_the_resolution_evaluate_tilted() {
    beyond(2);
}
