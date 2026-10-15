//! The consistency identities of a Boolean of two solids (REVIEW_NOTES.md,
//! "S9e.4b.4d.2c.1a's consistency fixed"): the six Booleans of `a` and `b`
//! (the common and the fuse both ways, both cuts) each valid at its own
//! resolution with enclosures at most 1e-9 wide, `cut(a, b) + common = a`,
//! `cut(b, a) + common = b`, `fuse = a + b - common`, and where both inputs
//! are symmetric in `y` every result's moment along `y` zero within its
//! enclosure; and the partners of a turned body of the consistency
//! generator's rows (`tools/generate_turned_consistency_fixtures.py`).
//! Shared by `turned_consistency_booleans.rs` and
//! `turned_crossing_booleans.rs`, each declaring `mod protocol` (the
//! Boolean protocol's support) at its root.
use rusty_occt::identity::OperationId;
use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};

/// A partner of a body: its solid, its volume and its common's with the
/// model (the reference's), and whether it is symmetric in `y`.
pub struct Case {
    pub name: String,
    pub partner: Solid,
    pub volume: f64,
    pub common: f64,
    pub symmetric: bool,
}

fn number(w: &str) -> f64 {
    w.parse().unwrap_or_else(|_| panic!("a number: {w}"))
}

/// A body's row (its path and the model's volume) and its partners, from
/// the consistency generator's rows: `body NAME PATH volume V`, then per
/// partner `case NAME KIND ... volume V common C`, the partner a `rib X0
/// W0 X1 W1 Y0 Y1` (the prism over the triangle `(X0, W0)`, `(X1, W0)`,
/// `(X0, W1)` in `(x, z)` from `y = Y0` to `Y1`), a `box X Y Z DX DY DZ`
/// (`Solid::box_at`'s) or a `post CX CY R Z0 Z1` (the cylinder of radius
/// `R` about `(CX, CY)` parallel to `z` from `Z0` to `Z1`).
pub fn partners(text: &str) -> ((String, f64), Vec<Case>) {
    let tol = Tolerance::default();
    let mut body = None;
    let mut cases = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let (tag, row) = line.split_once('\t').expect("a tagged row");
        let w: Vec<&str> = row.split(' ').collect();
        if tag == "body" {
            body = Some((w[1].to_string(), number(w[3])));
            continue;
        }
        let at = w.iter().position(|x| *x == "volume").expect("a volume");
        let v: Vec<f64> = w[2..at].iter().map(|x| number(x)).collect();
        let (partner, symmetric) = match w[1] {
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
                let s = Solid::extrude_with(OperationId(92), profile, frame, 0.0, v[5] - v[4])
                    .unwrap()
                    .0;
                (s, v[4] == -v[5])
            }
            "box" => {
                let s = Solid::box_at_with(
                    OperationId(92),
                    Point3::new(v[0], v[1], v[2]),
                    Vec3::new(v[3], v[4], v[5]),
                    tol,
                )
                .unwrap()
                .0;
                (s, v[1] == -(v[1] + v[4]))
            }
            "post" => {
                let profile = Profile::new(
                    Boundary::circle(Point2::new(v[0], v[1]), v[2], tol).unwrap(),
                    Vec::new(),
                    tol,
                )
                .unwrap();
                let frame = Frame3::new(
                    Point3::new(0.0, 0.0, v[3]),
                    Vec3::new(0.0, 0.0, 1.0),
                    Vec3::new(1.0, 0.0, 0.0),
                    tol,
                )
                .unwrap();
                let s = Solid::extrude_with(OperationId(92), profile, frame, 0.0, v[4] - v[3])
                    .unwrap()
                    .0;
                (s, v[1] == 0.0)
            }
            other => panic!("unknown partner {other}"),
        };
        cases.push(Case {
            name: w[0].to_string(),
            partner,
            volume: number(w[at + 1]),
            common: number(w[at + 3]),
            symmetric,
        });
    }
    (body.expect("the body row"), cases)
}

/// An imported solid of `fixtures/`.
pub fn imported(path: &str) -> Result<Solid, Error> {
    let (topology, resolution) = crate::protocol::imported_topology(path);
    Solid::imported_with(OperationId(91), topology, resolution).map(|(s, _)| s)
}

/// Certified enclosures of a set of solids: its volume and its moment along
/// `y` (each solid's volume times its centre's `y`, the volume's and the
/// centre's ends multiplied as intervals).
#[derive(Debug, Clone, Copy)]
pub struct Measure {
    pub volume: [f64; 2],
    pub moment: [f64; 2],
}

pub fn measure(out: &[Solid]) -> Result<Measure, String> {
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
pub fn meet(a: [f64; 2], b: [f64; 2], slack: f64) -> bool {
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
pub fn consistent(a: &Solid, b: &Solid, symmetric: bool) -> Result<[f64; 2], String> {
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

/// A body of the consistency generator's rows is its model (the reference's
/// volume by Pappus within 1e-9, its moment along `y` zero), and every
/// partner's six Booleans with it evaluate, agree with one another and keep
/// `y = 0` where the partner is symmetric in `y`, the common the model's
/// within 1e-9. The failures, one line each.
pub fn check_partners(text: &str) -> Vec<String> {
    let ((path, volume), cases) = partners(text);
    let body = imported(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let own = measure(std::slice::from_ref(&body)).unwrap();
    let near = |x: f64, e: [f64; 2]| meet(e, [x, x], 1e-9 * x.abs().max(1.0));
    let mut failures = Vec::new();
    if !near(volume, own.volume) {
        failures.push(format!("{path}: {} for {volume}", own.volume[0]));
    }
    if !meet(own.moment, [0.0, 0.0], 1e-12 * volume) {
        failures.push(format!("{path}: its moment along y {:?}", own.moment));
    }
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
    failures
}
