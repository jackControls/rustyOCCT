//! The fix of S9e.4b.4d.2c.1a's consistency (REVIEW_NOTES.md, "S9e.4b.4d.2c.1a's
//! consistency fixed"): an imported turned body of `OCC485a`'s form (a
//! bottom and a step within the resolution of planes, a fillet and a round
//! off their snapped tangencies, an axis tilted by 2^-27 rad, resolution
//! 3e-3; `fixtures/imported/consistency_ribbed.brep`, written by OCCT from
//! `tools/generate_turned_consistency_fixtures.py`) is stored as its model,
//! so that its Booleans with a rib through its step, fillet, shaft and round
//! and with boxes across its cylinders, its fillet and its round evaluate,
//! agree with one another within their enclosures (the fuse and the common
//! in both orders, each cut and the common summing to its input, the fuse
//! to both inputs less the common), keep `y = 0` where both inputs are
//! symmetric in `y`, and give the model's measures
//! (`fixtures/boolean-turned-consistency-expected.tsv`). The other turned
//! fixtures are checked the same way against a box across each.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::identity::OperationId;

use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};

/// A partner of the body: its solid, its volume and its common's with the
/// model (the reference's), and whether it is symmetric in `y`.
struct Case {
    name: String,
    partner: Solid,
    volume: f64,
    common: f64,
    symmetric: bool,
}

fn number(w: &str) -> f64 {
    w.parse().unwrap_or_else(|_| panic!("a number: {w}"))
}

/// The body (its path and the model's volume) and every partner.
fn fixture() -> ((String, f64), Vec<Case>) {
    let tol = Tolerance::default();
    let mut body = None;
    let mut cases = Vec::new();
    for line in include_str!("../../fixtures/boolean-turned-consistency-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (tag, row) = line.split_once('\t').expect("a tagged row");
        let w: Vec<&str> = row.split(' ').collect();
        if tag == "body" {
            body = Some((w[1].to_string(), number(w[3])));
            continue;
        }
        let v: Vec<f64> = w[2..8].iter().map(|x| number(x)).collect();
        let partner = match w[1] {
            "rib" => {
                // The triangle `(x0, w0)`, `(x1, w0)`, `(x0, w1)` in `(x, z)`
                // from `y = y1` down to `y0`.
                let points = vec![
                    Point2::new(v[0], v[1]),
                    Point2::new(v[2], v[1]),
                    Point2::new(v[0], v[3]),
                ];
                let profile =
                    Profile::new(Boundary::polygon(points, tol).unwrap(), Vec::new(), tol).unwrap();
                let frame = Frame3::new(
                    Point3::new(0.0, v[5], 0.0),
                    Vec3::new(0.0, -1.0, 0.0),
                    Vec3::new(1.0, 0.0, 0.0),
                    tol,
                )
                .unwrap();
                Solid::extrude_with(OperationId(92), profile, frame, 0.0, v[5] - v[4])
                    .unwrap()
                    .0
            }
            "box" => {
                Solid::box_at_with(
                    OperationId(92),
                    Point3::new(v[0], v[1], v[2]),
                    Vec3::new(v[3], v[4], v[5]),
                    tol,
                )
                .unwrap()
                .0
            }
            other => panic!("unknown partner {other}"),
        };
        let symmetric = match w[1] {
            "rib" => v[4] == -v[5],
            _ => v[1] == -(v[1] + v[4]),
        };
        cases.push(Case {
            name: w[0].to_string(),
            partner,
            volume: number(w[9]),
            common: number(w[11]),
            symmetric,
        });
    }
    (body.expect("the body row"), cases)
}

fn imported(path: &str) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(path);
    Solid::imported_with(OperationId(91), topology, resolution).map(|(s, _)| s)
}

/// Certified enclosures of a set of solids: its volume and its moment along
/// `y` (each solid's volume times its centre's `y`, the volume's and the
/// centre's ends multiplied as intervals).
#[derive(Debug, Clone, Copy)]
struct Measure {
    volume: [f64; 2],
    moment: [f64; 2],
}

fn measure(out: &[Solid]) -> Result<Measure, String> {
    let mut m = Measure {
        volume: [0.0; 2],
        moment: [0.0; 2],
    };
    for s in out {
        let issues = s.topology().check(s.resolution());
        if !issues.is_empty() {
            return Err(format!("an invalid result: {issues:?}"));
        }
        let e = s
            .topology()
            .mass_enclosure()
            .ok_or("a result's mass undecided")?;
        let products = [
            e.volume[0] * e.centroid[1][0],
            e.volume[0] * e.centroid[1][1],
            e.volume[1] * e.centroid[1][0],
            e.volume[1] * e.centroid[1][1],
        ];
        m.volume[0] += e.volume[0];
        m.volume[1] += e.volume[1];
        m.moment[0] += products.iter().copied().fold(f64::INFINITY, f64::min);
        m.moment[1] += products.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    }
    Ok(m)
}

/// Whether two enclosures of one value meet, within `slack`.
fn meet(a: [f64; 2], b: [f64; 2], slack: f64) -> bool {
    a[0] - slack <= b[1] && b[0] - slack <= a[1]
}

fn add(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[1], a[1] - b[0]]
}

/// The six Booleans of `a` and `b` checked against one another: the fuse and
/// the common in both orders, `cut(a, b) + common = a`, `cut(b, a) + common
/// = b`, `fuse = a + b - common`; each result valid, its enclosures at most
/// 1e-9 wide (relative), and where both inputs are symmetric in `y` every
/// result's moment along `y` zero within its enclosure. The common's
/// enclosure.
fn consistent(a: &Solid, b: &Solid, symmetric: bool) -> Result<[f64; 2], String> {
    let run = |r: Result<(Vec<Solid>, rusty_occt::history::History), Error>, what: &str| {
        r.map_err(|e| format!("{what}: {e}"))
            .and_then(|(out, _)| measure(&out).map_err(|e| format!("{what}: {e}")))
    };
    let op = OperationId(93);
    let ab = run(a.common(op, b), "common")?;
    let ba = run(b.common(op, a), "the common the other way")?;
    let cut_ab = run(a.cut(op, b), "the cut")?;
    let cut_ba = run(b.cut(op, a), "the other cut")?;
    let fuse_ab = run(a.fuse(op, b), "the fuse")?;
    let fuse_ba = run(b.fuse(op, a), "the fuse the other way")?;
    let inputs = [
        measure(std::slice::from_ref(a))?,
        measure(std::slice::from_ref(b))?,
    ];
    let scale = inputs[0].volume[1].max(inputs[1].volume[1]);
    let slack = 1e-12 * scale;
    let mut bad = Vec::new();
    let all = [
        ("common", ab),
        ("common'", ba),
        ("cut", cut_ab),
        ("cut'", cut_ba),
        ("fuse", fuse_ab),
        ("fuse'", fuse_ba),
    ];
    for (what, m) in &all {
        if m.volume[1] - m.volume[0] > 1e-9 * m.volume[1].abs().max(1.0) {
            bad.push(format!("{what}: a wide enclosure {:?}", m.volume));
        }
        if symmetric && !meet(m.moment, [0.0, 0.0], slack) {
            bad.push(format!("{what}: its moment along y {:?}, not 0", m.moment));
        }
    }
    let checks = [
        ("the common's two orders", ab.volume, ba.volume),
        ("the fuse's two orders", fuse_ab.volume, fuse_ba.volume),
        (
            "cut + common = a",
            add(cut_ab.volume, ab.volume),
            inputs[0].volume,
        ),
        (
            "cut' + common = b",
            add(cut_ba.volume, ab.volume),
            inputs[1].volume,
        ),
        (
            "fuse = a + b - common",
            fuse_ab.volume,
            sub(add(inputs[0].volume, inputs[1].volume), ab.volume),
        ),
    ];
    for (what, x, y) in checks {
        if !meet(x, y, slack) {
            bad.push(format!("{what}: {x:?} against {y:?}"));
        }
    }
    if bad.is_empty() {
        Ok(ab.volume)
    } else {
        Err(bad.join("; "))
    }
}

/// The body is its model (the reference's volume by Pappus within 1e-9),
/// and every partner's six Booleans with it evaluate, agree with one another
/// and keep `y = 0` where the partner is symmetric in `y`, the common the
/// model's.
#[test]
fn the_ribbed_body_s_booleans_are_consistent() {
    let ((path, volume), cases) = fixture();
    let body = imported(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let own = measure(std::slice::from_ref(&body)).unwrap();
    let near = |x: f64, e: [f64; 2]| meet(e, [x, x], 1e-9 * x.abs().max(1.0));
    assert!(near(volume, own.volume), "{} for {volume}", own.volume[0]);
    assert!(
        meet(own.moment, [0.0, 0.0], 1e-12 * volume),
        "{:?}",
        own.moment
    );
    let mut failures = Vec::new();
    for case in &cases {
        let partner = measure(std::slice::from_ref(&case.partner)).unwrap();
        if !near(case.volume, partner.volume) {
            failures.push(format!("{}: the partner {:?}", case.name, partner.volume));
        }
        match consistent(&body, &case.partner, case.symmetric) {
            Ok(common) if near(case.common, common) => {}
            Ok(common) => failures.push(format!(
                "{}: the common {common:?} for {}",
                case.name, case.common
            )),
            Err(e) => failures.push(format!("{}: {e}", case.name)),
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The other imported turned fixtures, each against a box across it (from
/// inside past its side, across its middle heights): their six Booleans
/// evaluate and agree with one another.
#[test]
fn other_turned_bodies_booleans_are_consistent() {
    let tol = Tolerance::default();
    let mut failures = Vec::new();
    for name in [
        "turned_fillet",
        "turned_pin",
        "turned_step",
        "fillet_step",
        "fillet_collar",
        "fillet_refit",
        "flat_pin",
        "flat_collar",
        "profile_crown",
        "profile_point",
        "profile_steps",
        "profile_waist",
        "tilted_waist",
        "tilted_pin",
    ] {
        let body = match imported(&format!("imported/{name}.brep")) {
            Ok(b) => b,
            Err(e) => {
                failures.push(format!("{name}: {e}"));
                continue;
            }
        };
        let b = body.bounds();
        let (lo, hi) = (b.min.to_array(), b.max.to_array());
        let at = |k: usize, f: f64| lo[k] + f * (hi[k] - lo[k]);
        let start = Point3::new(at(0, 0.5625), at(1, 0.3125), at(2, 0.28125));
        let size = Vec3::new(
            hi[0] - lo[0],
            (hi[1] - lo[1]) * 0.40625,
            (hi[2] - lo[2]) * 0.46875,
        );
        let partner = Solid::box_at_with(OperationId(92), start, size, tol)
            .unwrap()
            .0;
        if let Err(e) = consistent(&body, &partner, false) {
            failures.push(format!("{name}: {e}"));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
