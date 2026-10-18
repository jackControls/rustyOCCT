//! A Boolean's resolution (REVIEW_NOTES.md, "A Boolean's resolution"): a
//! Boolean works and returns its result at its inputs' coarser resolution,
//! whatever their order. Imported bodies (their resolution their file's
//! tolerance, `1e-7` to `3e-3`) against a kernel box at `1e-7` in both
//! orders: each result valid at its own resolution, the coarser one, the
//! two orders' commons and fuses alike, every history complete, results
//! deterministic, and each Boolean's history and result ids
//! those of the kernel before the rule (`fixtures/boolean-resolution-
//! histories.tsv`, written by the kernel at `37949a87`, whose results in the
//! tool's order were invalid at their own resolution: the imported body's
//! kept entities' enclosures, its file's tolerance, past the object's);
//! kernel prisms of two resolutions on each engine (one frame, a stack, any
//! position, arcs) at the coarser in both orders, an input returned whole
//! keeping its own, and a near miss decided at the coarser in both orders.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, History};
use rusty_occt::identity::OperationId;
use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
use std::collections::BTreeMap;

/// Imported bodies of several primitives (planes, cylinders, spheres, a
/// turned body's), their resolutions their files' tolerances.
const BODIES: &[&str] = &[
    "bitten",
    "consistency_ribbed",
    "corner_shaft",
    "crossed_joint",
    "deep_post",
    "flush_arms",
    "form_dee_boss",
    "incidence_half",
    "lune",
    "nested_post",
    "octant_tilt",
    "part_stack",
    "rim_well",
    "sphere_wedge",
    "stack_step",
    "tilted_skewbox",
    "two_through",
    "zone_wedge",
];

const OPS: [&str; 3] = ["common", "cut", "fuse"];

fn tolerance(linear: f64) -> Tolerance {
    Tolerance::new(linear, 1e-12).unwrap()
}

fn body(name: &str) -> Solid {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(91), topology, resolution)
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .0
}

/// A box at `1e-7` over most of a body's bounds, off its corners by
/// fractions no face of the body lies at.
fn partner(body: &Solid) -> Solid {
    let b = body.bounds();
    let (lo, hi) = (b.min.to_array(), b.max.to_array());
    let o: [f64; 3] = std::array::from_fn(|i| lo[i] + (hi[i] - lo[i]) * 0.3137);
    let s: [f64; 3] = std::array::from_fn(|i| (hi[i] - lo[i]) * 0.9071);
    Solid::box_at_with(
        OperationId(92),
        Point3::new(o[0], o[1], o[2]),
        Vec3::new(s[0], s[1], s[2]),
        tolerance(1e-7),
    )
    .unwrap()
    .0
}

fn boolean(a: &Solid, op: &str, b: &Solid) -> Result<(Vec<Solid>, History), Error> {
    let id = OperationId(93);
    match op {
        "common" => a.common(id, b),
        "cut" => a.cut(id, b),
        _ => a.fuse(id, b),
    }
}

/// FNV-1a of a text: stable across Rust releases, unlike `DefaultHasher`.
fn fnv(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    })
}

/// A Boolean's outcome as a digest: its history and each result's ids by
/// slot, with `points` its vertices' bits too (the pinned digests, written
/// before R13 when macOS's `hypot` and glibc's differed by an ulp, leave
/// them out; `construction_bits.rs` pins the bodies' bits); or its error.
fn digest(r: &Result<(Vec<Solid>, History), Error>, points: bool) -> String {
    match r {
        Err(e) => format!("error {e}"),
        Ok((out, h)) => {
            let mut text = format!("{h:?}");
            for s in out {
                let ids: Vec<_> = s.topology().ids().collect();
                text += &format!("\n{:?} {ids:?}", s.topology().body_id());
                if points {
                    let bits: Vec<[u64; 3]> = s
                        .topology()
                        .vertices()
                        .iter()
                        .map(|v| v.position.to_array().map(f64::to_bits))
                        .collect();
                    text += &format!(" {bits:?}");
                }
            }
            format!("{:016x}", fnv(&text))
        }
    }
}

/// The results' volume enclosure, each valid at its own resolution, that
/// resolution `expected` (an input returned whole its own).
fn checked(
    out: &[Solid],
    inputs: [&Solid; 2],
    expected: Tolerance,
    history: &History,
) -> Result<[f64; 2], String> {
    let mut v = [0.0, 0.0];
    for s in out {
        let whole = inputs
            .iter()
            .find(|x| x.topology().body_id() == s.topology().body_id());
        let want = whole.map_or(expected, |x| x.resolution());
        if s.resolution() != want {
            return Err(format!("the resolution {:?} for {want:?}", s.resolution()));
        }
        let issues = s.topology().check(s.resolution());
        if !issues.is_empty() {
            return Err(format!("invalid at its own resolution: {issues:?}"));
        }
        let e = s.topology().mass_enclosure().ok_or("a mass undecided")?;
        v[0] += e.volume[0];
        v[1] += e.volume[1];
    }
    let ins: Vec<_> = inputs
        .iter()
        .map(|s| s.topology().entity_set(s.resolution()))
        .collect();
    let outs: Vec<_> = out
        .iter()
        .map(|s| s.topology().entity_set(s.resolution()))
        .collect();
    let issues = history::check(&ins, &outs, history);
    if !issues.is_empty() {
        return Err(format!("the history: {issues:?}"));
    }
    Ok(v)
}

fn meet(a: [f64; 2], b: [f64; 2], slack: f64) -> bool {
    a[0] - slack <= b[1] && b[0] - slack <= a[1]
}

/// A Boolean's row: the body, the order, the operation, the digest and a
/// failure.
type Row = (String, String, String, String, Option<String>);

/// Each body's and its box's Booleans: `(body, order, op)` to the digest.
fn outcomes(failures: &mut Vec<String>) -> BTreeMap<(String, String, String), String> {
    let rows: Vec<Vec<Row>> = std::thread::scope(|scope| {
        let handles: Vec<_> = BODIES
            .iter()
            .map(|name| {
                scope.spawn(move || {
                    let b = body(name);
                    let x = partner(&b);
                    let coarser = b.resolution().coarser(x.resolution());
                    let mut rows = Vec::new();
                    let mut volumes = BTreeMap::new();
                    for (order, p, q) in [("object", &b, &x), ("tool", &x, &b)] {
                        for op in OPS {
                            let r = boolean(p, op, q);
                            let again = boolean(p, op, q);
                            let d = digest(&r, false);
                            let mut failure = None;
                            if digest(&again, true) != digest(&r, true) {
                                failure = Some("not deterministic".to_string());
                            }
                            if let Ok((out, h)) = &r {
                                match checked(out, [p, q], coarser, h) {
                                    Ok(v) => {
                                        volumes.insert((order, op), v);
                                    }
                                    Err(e) => failure = Some(e),
                                }
                            }
                            rows.push((
                                name.to_string(),
                                order.to_string(),
                                op.to_string(),
                                d,
                                failure,
                            ));
                        }
                    }
                    // The common and the fuse alike in both orders.
                    for op in ["common", "fuse"] {
                        if let (Some(u), Some(v)) =
                            (volumes.get(&("object", op)), volumes.get(&("tool", op)))
                        {
                            if !meet(*u, *v, 1e-12 * u[1].abs().max(1.0)) {
                                rows.push((
                                    name.to_string(),
                                    "both".to_string(),
                                    op.to_string(),
                                    String::new(),
                                    Some(format!("the orders' volumes {u:?} and {v:?}")),
                                ));
                            }
                        }
                    }
                    rows
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut out = BTreeMap::new();
    for (name, order, op, d, failure) in rows.into_iter().flatten() {
        if let Some(f) = failure {
            failures.push(format!("{name} {order} {op}: {f}"));
        }
        if !d.is_empty() {
            out.insert((name, order, op), d);
        }
    }
    out
}

/// Every imported body against its box, both orders: each result valid at
/// the coarser resolution, the orders alike, histories complete, results
/// deterministic; each outcome (history, result ids) the kernel's before
/// the rule. `WRITE_DIGESTS=1` prints the rows instead.
#[test]
fn imported_bodies_and_boxes_at_the_coarser_resolution() {
    let mut failures = Vec::new();
    let got = outcomes(&mut failures);
    if std::env::var("WRITE_DIGESTS").is_ok() {
        for ((name, order, op), d) in &got {
            println!("{name}\t{order} {op} {d}");
        }
        return;
    }
    let mut expected = BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-resolution-histories.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        let (order, rest) = row.split_once(' ').unwrap();
        let (op, d) = rest.split_once(' ').unwrap();
        expected.insert(
            (name.to_string(), order.to_string(), op.to_string()),
            d.to_string(),
        );
    }
    assert_eq!(expected.len(), BODIES.len() * 6, "the fixture's rows");
    for (k, d) in &got {
        if expected.get(k) != Some(d) {
            failures.push(format!("{k:?}: {d} for {:?}", expected.get(k)));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

fn prism(op: u64, points: &[(f64, f64)], lo: f64, hi: f64, frame: Frame3, res: f64) -> Solid {
    let pts = points.iter().map(|(x, y)| Point2::new(*x, *y)).collect();
    let profile = Profile::new(
        Boundary::polygon(pts, tolerance(res)).unwrap(),
        Vec::new(),
        tolerance(res),
    )
    .unwrap();
    Solid::extrude_with(OperationId(op), profile, frame, lo, hi)
        .unwrap()
        .0
}

fn rod(op: u64, c: (f64, f64), r: f64, lo: f64, hi: f64, frame: Frame3, res: f64) -> Solid {
    let profile = Profile::new(
        Boundary::circle(Point2::new(c.0, c.1), r, tolerance(res)).unwrap(),
        Vec::new(),
        tolerance(res),
    )
    .unwrap();
    Solid::extrude_with(OperationId(op), profile, frame, lo, hi)
        .unwrap()
        .0
}

fn world() -> Frame3 {
    Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tolerance(1e-7),
    )
    .unwrap()
}

/// A frame turned off every axis (a Pythagorean normal), so a pair is the
/// polyhedral or the curved engine's.
fn turned() -> Frame3 {
    Frame3::new(
        Point3::new(0.25, 0.125, -0.5),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tolerance(1e-7),
    )
    .unwrap()
}

/// Kernel prisms of resolutions `1e-7` and `1e-5` on each engine: in one
/// frame over one height range (S9a), as a stack (S9a.2), turned (S9b) and
/// with an arc (S9c); each result at `1e-5` in both orders, valid there,
/// the orders' commons and fuses alike; an input a cut misses returned
/// whole at its own.
#[test]
fn kernel_prisms_of_two_resolutions_on_each_engine() {
    let square = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
    let shifted = [(1.0, 0.5), (3.0, 0.5), (3.0, 2.5), (1.0, 2.5)];
    // Across the turned frame, holding the box's point `(1.5, 1, 0.5)`.
    let across = [(1.0, -0.5), (3.0, -0.5), (3.0, 1.5), (1.0, 1.5)];
    let pairs: Vec<(&str, Solid, Solid)> = vec![
        (
            "one frame",
            prism(1, &square, 0.0, 1.0, world(), 1e-7),
            prism(2, &shifted, 0.0, 1.0, world(), 1e-5),
        ),
        (
            "a stack",
            prism(1, &square, 0.0, 1.0, world(), 1e-7),
            prism(2, &shifted, 0.5, 1.5, world(), 1e-5),
        ),
        (
            "turned",
            prism(1, &square, 0.0, 1.0, world(), 1e-7),
            prism(2, &across, 0.75, 1.75, turned(), 1e-5),
        ),
        (
            "an arc",
            prism(1, &square, 0.0, 1.0, world(), 1e-7),
            rod(2, (1.0, 0.5), 0.625, 0.75, 1.75, turned(), 1e-5),
        ),
    ];
    let coarser = tolerance(1e-5);
    let mut failures = Vec::new();
    for (what, a, b) in &pairs {
        let mut volumes = BTreeMap::new();
        for (order, p, q) in [("object", a, b), ("tool", b, a)] {
            for op in OPS {
                match boolean(p, op, q) {
                    Ok((out, h)) => match checked(&out, [p, q], coarser, &h) {
                        Ok(v) => {
                            // Every result built (none an input whole).
                            if out.is_empty() || out.iter().any(|s| s.resolution() != coarser) {
                                failures.push(format!("{what} {order} {op}: not built at 1e-5"));
                            }
                            volumes.insert((order, op), v);
                        }
                        Err(e) => failures.push(format!("{what} {order} {op}: {e}")),
                    },
                    Err(e) => failures.push(format!("{what} {order} {op}: {e}")),
                }
            }
        }
        for op in ["common", "fuse"] {
            let (u, v) = (volumes.get(&("object", op)), volumes.get(&("tool", op)));
            if let (Some(u), Some(v)) = (u, v) {
                if !meet(*u, *v, 1e-12) {
                    failures.push(format!("{what} {op}: {u:?} and {v:?}"));
                }
            }
        }
    }
    // A cut missing its object returns it whole, at its own resolution.
    let far = prism(3, &square, 5.0, 6.0, world(), 1e-5);
    let (out, _) = pairs[0].1.cut(OperationId(93), &far).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].resolution(), tolerance(1e-7));
    assert!(failures.is_empty(), "{failures:#?}");
}

/// A near miss decided at the coarser resolution whatever the order: a
/// prism at `1e-5` whose top lies `4e-6` above a prism's at `1e-7` in one
/// frame leaves a slab thinner than the coarser resolution in the fuse,
/// refused in both orders (before the rule, taken at `1e-7` with the finer
/// object first and refused the other way).
#[test]
fn a_near_miss_is_decided_at_the_coarser_resolution_in_both_orders() {
    let square = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
    let shifted = [(1.0, 0.5), (3.0, 0.5), (3.0, 2.5), (1.0, 2.5)];
    let fine = prism(1, &square, 0.0, 1.0, world(), 1e-7);
    let coarse = prism(2, &shifted, 0.0, 1.000004, world(), 1e-5);
    for (p, q) in [(&fine, &coarse), (&coarse, &fine)] {
        match p.fuse(OperationId(93), q) {
            Err(Error::Degenerate(_)) => {}
            other => panic!("{:?}", other.map(|(out, _)| out.len())),
        }
    }
}
