//! A whole torus of radii 3 and 1 and a rod about its axis whose wall
//! crosses the torus's inner equator by `g` (radius `2 + g`, coaxial: a
//! meeting near a double root along a whole circle), in a level frame, one
//! turned about its axis and one tilted off the world's axes (the rod's
//! frame built in the torus's, so its axis is the torus's within rounding).
//! In the tilted frame a crossing of `1e-9` ran for over half an hour in
//! the torus meeting's verification (`ToricCrv::verified`), its runs split
//! and their pieces verified again up to 48 times, each verification capped
//! but not their total; the meeting's verifications now draw on one budget
//! (REVIEW_NOTES.md, "A torus meeting's verification budget").
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

/// The crossing of `1e-6` (beyond the resolution): each operation one
/// solid, of the reference's volumes (mpmath, 25 digits: the torus
/// `6 pi^2`, the rod `3 pi rho^2` at `rho = 2 + 1e-6` in binary64, their
/// common by Pappus over the tube's segment, `2 pi` times the integral of
/// `2 x sqrt(1 - (x - 3)^2)` over `[2, rho]`), the fuse and the common
/// together the inputs, the cut and the common the torus.
fn evaluates(what: &str, t: &Solid, r: &Solid) {
    const TORUS: f64 = 59.217_626_406_536_15;
    const ROD: f64 = 37.699_149_542_198_79;
    const COMMON: f64 = 2.369_537_922_945_086e-8;
    let (vt, vr) = (t.mass_properties().volume, r.mass_properties().volume);
    let close = |x: f64, y: f64| (x - y).abs() <= 1e-12 * (TORUS + ROD);
    assert!(close(vt, TORUS) && close(vr, ROD), "{what}: {vt} {vr}");
    let [fuse, cut, back, common] = outcomes(t, r).map(|out| out.expect(what));
    for (op, s) in [
        ("fuse", &fuse),
        ("cut", &cut),
        ("cut back", &back),
        ("common", &common),
    ] {
        assert_eq!(s.len(), 1, "{what} {op}");
    }
    let [vf, vc, vb, vm] = [&fuse, &cut, &back, &common].map(|s| volume(s));
    assert!(close(vm, COMMON), "{what} common {vm}");
    assert!(close(vf, TORUS + ROD - COMMON), "{what} fuse {vf}");
    assert!(close(vc, TORUS - COMMON), "{what} cut {vc}");
    assert!(close(vb, ROD - COMMON), "{what} cut back {vb}");
    assert!(close(vf + vm, vt + vr), "{what} fuse and common");
    assert!(close(vc + vm, vt), "{what} cut and common");
    assert!(close(vb + vm, vr), "{what} cut back and common");
}

/// The rod tangent along the inner equator (`g = 0`) and crossing it by
/// `1e-12`, `1e-9`, `5e-8` and `1e-6`. Coaxial exactly (level, turned): the
/// tangency `Degenerate`, the crossings within the resolution circles the
/// box verification leaves unverified (`ComputationLimit`), and `1e-6`
/// evaluating. Coaxial within rounding (tilted), each a `ComputationLimit`:
/// the tangency and `1e-12` a stalled trace; `1e-9` the meeting over its
/// verification budget where the platform's `hypot` is macOS's (about a
/// second a meeting in release; over half an hour before), a stalled trace
/// under glibc's correctly rounded one (the rod's frame normalized
/// again); `5e-8` over the budget under both; its `1e-6` is
/// `a_tilted_rod_crossing_beyond_the_resolution_evaluates`.
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
        if coaxial {
            evaluates(&at(1e-6), &t, &rod(&f, 1e-6));
        }
    }
}

/// The tilted frame's crossing of `1e-6`: it evaluates as the others do,
/// its meeting's pieces verified at once, but validating its results'
/// faces (their projections' integrals along the meeting near a double
/// root) takes minutes in release (HANDOFF.md).
#[test]
#[ignore]
fn a_tilted_rod_crossing_beyond_the_resolution_evaluates() {
    let (_, f) = frames()[2];
    evaluates("tilted 1e-6", &torus(f), &rod(&f, 1e-6));
}
