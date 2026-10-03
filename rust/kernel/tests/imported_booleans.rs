//! S9e.4a: imported solids (bodies OCCT wrote to `.brep` files, read and
//! converted, without a construction) given to Booleans, against the
//! independent reference (`fixtures/boolean-imported-*` from
//! `tools/generate_imported_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`): each case's Boolean (or chain) with its imported
//! inputs made by `Solid::imported_with`, decided on the construction their
//! stored surfaces give. Each case runs once (on a few threads) for the
//! checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Location, Point3, Solid};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!("../../fixtures/boolean-imported-cases.txt"))
}

/// `f` over the cases on at most six threads, in the cases' order.
fn each<T: Send>(f: impl Fn(&protocol::Case) -> T + Sync) -> Vec<(String, T)> {
    let cases = cases();
    let next = AtomicUsize::new(0);
    let workers = std::thread::available_parallelism()
        .map_or(2, |n| n.get())
        .clamp(1, 6);
    let mut out: Vec<(usize, String, T)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    let mut mine = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::SeqCst);
                        let Some(case) = cases.get(i) else { break };
                        mine.push((i, case.name.clone(), f(case)));
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
    out.into_iter().map(|(_, n, t)| (n, t)).collect()
}

/// Every case's run (the last Boolean's inputs, results and history), once.
fn runs() -> &'static [(String, Result<protocol::Run, Error>)] {
    static RUNS: OnceLock<Vec<(String, Result<protocol::Run, Error>)>> = OnceLock::new();
    RUNS.get_or_init(|| each(protocol::run))
}

/// The reference's kind and `(solids, [volume, area])` per case.
type Expected = std::collections::BTreeMap<String, (String, Option<(usize, [f64; 2])>)>;

fn expected() -> Expected {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-imported-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let entry = expect
            .entry(name.to_string())
            .or_insert((String::new(), None));
        match w[0] {
            "expect" => entry.0 = w[1].to_string(),
            "result" => {
                let v: Vec<f64> = w[2..4].iter().map(|x| x.parse().unwrap()).collect();
                entry.1 = Some((w[1].parse::<usize>().unwrap(), [v[0], v[1]]));
            }
            _ => {}
        }
    }
    expect
}

/// A result's volume and area enclosures within the reference's, each at
/// most `1e-9` wide (relative).
fn within(out: &[Solid], count: usize, v: [f64; 2]) -> Result<(), String> {
    let mut sums = [0.0f64; 4];
    for s in out {
        let m = s
            .topology()
            .mass_enclosure()
            .ok_or("a result's mass undecided")?;
        for (k, x) in [
            m.volume[0],
            m.volume[1],
            m.surface_area[0],
            m.surface_area[1],
        ]
        .into_iter()
        .enumerate()
        {
            sums[k] += x;
        }
    }
    if out.len() != count {
        return Err(format!("{} solids for {count}", out.len()));
    }
    let near = |x: f64, lo: f64, hi: f64| {
        let slack = 1e-9 * x.abs().max(1.0);
        lo - slack <= x && x <= hi + slack
    };
    let narrow = |lo: f64, hi: f64| hi - lo <= 1e-9 * lo.abs().max(hi.abs()).max(1.0);
    if !narrow(sums[0], sums[1]) || !narrow(sums[2], sums[3]) {
        return Err(format!("wide enclosures {sums:?}"));
    }
    if !near(v[0], sums[0], sums[1]) || !near(v[1], sums[2], sums[3]) {
        return Err(format!("volume and area {v:?} against {sums:?}"));
    }
    Ok(())
}

#[test]
fn every_case_matches_the_reference() {
    let expect = expected();
    let mut failures = Vec::new();
    for (name, run) in runs() {
        let (kind, want) = &expect[name];
        match (kind.as_str(), run) {
            ("degenerate", Err(Error::Degenerate(_))) => continue,
            ("unsupported", Err(Error::OutOfDomain(_))) => continue,
            ("degenerate" | "unsupported", other) => {
                failures.push(format!(
                    "{name}: {:?} not refused as {kind}",
                    other.as_ref().map(|_| ())
                ));
                continue;
            }
            (_, Err(e)) => {
                failures.push(format!("{name}: {e}"));
                continue;
            }
            ("empty", Ok((_, _, out, _))) => {
                if !out.is_empty() {
                    failures.push(format!("{name}: {} solids, not empty", out.len()));
                }
                continue;
            }
            _ => {}
        }
        let Ok((_, _, out, _)) = run else {
            unreachable!("refused above")
        };
        let (count, v) = want.expect("a result");
        if let Err(e) = within(out, count, v) {
            failures.push(format!("{name}: {e}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures: {failures:#?}",
        failures.len()
    );
}

/// The turned profile's arc, rounded off its circle in its cap's frame, is
/// refused as S9e.4b's; the declared tangencies as `Degenerate`.
#[test]
fn refusals_name_their_reasons() {
    for (name, run) in runs() {
        if name.starts_with("dee_turn") {
            assert!(
                matches!(run, Err(Error::OutOfDomain(m)) if m.contains("round off its circle")),
                "{name}: {run:?}"
            );
        }
        if ["cyl_tangent", "ball_touch", "box_kiss"]
            .iter()
            .any(|p| name.starts_with(p))
        {
            assert!(matches!(run, Err(Error::Degenerate(_))), "{name}: {run:?}");
        }
    }
}

/// Each Boolean's history is complete over its inputs, the imported ones'
/// stored ids among them: every entity of every input resolved.
#[test]
fn histories_are_complete_over_the_imported_ids() {
    let failures: Vec<String> = runs()
        .iter()
        .filter_map(|(name, run)| run.as_ref().ok().map(|r| (name, r)))
        .flat_map(|(name, (a, b, out, h))| {
            let mut bad = Vec::new();
            let ins = [
                a.topology().entity_set(a.resolution()),
                b.topology().entity_set(b.resolution()),
            ];
            let outs: Vec<_> = out
                .iter()
                .map(|s| s.topology().entity_set(s.resolution()))
                .collect();
            let issues = history::check(&ins, &outs, h);
            if !issues.is_empty() {
                bad.push(format!("{name}: {issues:?}"));
            }
            let inputs: std::collections::BTreeSet<_> = [a, b]
                .iter()
                .flat_map(|s| s.topology().ids().map(|(id, _)| id))
                .collect();
            for s in [a, b] {
                for (id, _) in s.topology().ids() {
                    if matches!(h.resolve(id), Resolution::Unknown)
                        && !h.relations.iter().any(|r| r.sources().contains(&id))
                    {
                        bad.push(format!("{name}: an input entity unnamed"));
                    }
                }
            }
            // No relation names an entity of neither input (a
            // construction's own id).
            for r in &h.relations {
                if r.sources().iter().any(|id| !inputs.contains(id)) {
                    bad.push(format!("{name}: a relation from outside the inputs"));
                }
            }
            bad
        })
        .collect();
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
fn results_are_deterministic_and_move_rigidly() {
    use rusty_occt::{RigidTransform, Vec3};
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    let first: std::collections::BTreeMap<&str, _> =
        runs().iter().map(|(n, r)| (n.as_str(), r)).collect();
    let failures: Vec<String> = each(|case| {
        let Ok((_, _, out, h)) = first[case.name.as_str()] else {
            return Vec::new();
        };
        let mut bad = Vec::new();
        let (_, _, again, h2) = protocol::run(case).unwrap();
        if *h != h2 {
            bad.push(format!("{}: histories differ", case.name));
        }
        for (x, y) in out.iter().zip(&again) {
            if ids(x) != ids(y) || x.topology().edges() != y.topology().edges() {
                bad.push(format!("{}: results differ", case.name));
            }
        }
        for s in out {
            let (moved, _) = s
                .transform_with(OperationId(900), motion)
                .unwrap_or_else(|e| panic!("{}: {e}", case.name));
            if ids(&moved) != ids(s) {
                bad.push(format!("{}: moved ids", case.name));
            }
            let (v0, v1) = (s.mass_properties().volume, moved.mass_properties().volume);
            if (v0 - v1).abs() > 1e-9 * v0.abs().max(1.0) {
                bad.push(format!("{}: moved volume {v0} {v1}", case.name));
            }
            for v in moved.topology().vertices() {
                if moved.classify(v.position).unwrap() != Location::Boundary {
                    bad.push(format!("{}: a moved vertex off the boundary", case.name));
                }
            }
        }
        bad
    })
    .into_iter()
    .flat_map(|(_, b)| b)
    .collect();
    assert!(failures.is_empty(), "{failures:#?}");
}

/// Both inputs moved rigidly keep the reference's volumes: an imported
/// solid moves with its stored geometry and its construction (by binary64
/// steps, so the turned profiles' arcs keep their ends on their circles).
#[test]
fn moved_inputs_keep_their_volumes() {
    use rusty_occt::{RigidTransform, Vec3};
    let motion = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let expect = expected();
    let all = cases();
    for name in [
        "box_slab_cut",
        "tilt_box_fuse",
        "box_cyl_common",
        "dee_slab_cut",
        "halves_box_fuse",
        "plate_bore_cut",
        "ball_box_common",
        "dome_box_cut",
        "frustum_box_fuse",
        "spike_ball_cut",
        "ring_pin_cut",
        "both_fuse",
    ] {
        let case = all.iter().find(|c| c.name == name).unwrap();
        let (a, b, _, _) = protocol::run_first(case).unwrap();
        let a = a.transform_with(OperationId(801), motion).unwrap().0;
        let b = b.transform_with(OperationId(802), motion).unwrap().0;
        let out = match case.op.as_str() {
            "fuse" => a.fuse(case.operation, &b),
            "cut" => a.cut(case.operation, &b),
            _ => a.common(case.operation, &b),
        }
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .0;
        let (count, v) = expect[name].1.unwrap();
        within(&out, count, v).unwrap_or_else(|e| panic!("{name}: {e}"));
    }
}

/// A body file's imported solid under an operation.
fn body(name: &str, op: u64) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(op), topology, resolution).map(|(s, _)| s)
}

/// Every body is its construction: its mass the construction's in closed
/// form, its stored vertices on its boundary and points off it classified,
/// its ids the stored topology's under the import's operation (two imports
/// apart), no profile.
#[test]
fn imported_solids_are_their_constructions() {
    use std::f64::consts::PI;
    let dee = 6.0 * (100.0 + 12.5 * PI);
    let bodies = [
        ("box", 400.0, [5.0, 5.0, 2.0]),
        ("box_tilt", 120.0, [5.0, 4.2, 0.1]),
        ("cyl", 54.0 * PI, [5.0, 5.0, 2.0]),
        ("cyl_side", 27.0 * PI, [5.0, 5.0, 2.0]),
        ("dee", dee, [5.0, 5.0, 3.0]),
        ("halves", 31.25 * PI, [5.0, 5.0, 2.5]),
        ("plate", 3.0 * (144.0 - 4.0 * PI), [4.0, 0.0, 1.5]),
        ("ball", 36.0 * PI, [5.0, 5.0, 4.0]),
        ("dome", 128.0 * PI / 3.0, [0.0, 0.0, 1.0]),
        ("frustum", 52.0 * PI / 3.0, [5.0, 5.0, 2.0]),
        ("spike", 20.0 * PI / 3.0, [0.0, 0.0, 2.0]),
        ("ring", 22.5 * PI * PI, [5.0, 0.0, 0.0]),
        ("dee_turn", dee, [3.0, 4.0, 3.0]),
    ];
    for (name, volume, inside) in bodies {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let m = s.mass_properties();
        assert!(
            (m.volume - volume).abs() <= 1e-9 * volume,
            "{name}: {} for {volume}",
            m.volume
        );
        assert!(s.profile().is_none(), "{name}");
        let [x, y, z] = inside;
        assert_eq!(
            s.classify(Point3::new(x, y, z)).unwrap(),
            Location::Inside,
            "{name}"
        );
        assert_eq!(
            s.classify(Point3::new(x + 50.0, y, z)).unwrap(),
            Location::Outside,
            "{name}"
        );
        for v in s.topology().vertices() {
            assert_eq!(
                s.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
        let other = body(name, 92).unwrap();
        let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
        assert!(
            ids(&s).iter().all(|id| !ids(&other).contains(id)),
            "{name}: two imports share ids"
        );
    }
}

/// A STEP body is imported the same way: OCCT's STEP solids of the import
/// track's fixtures, decided on their constructions (a box, cylinders, a
/// sphere, a hemisphere, cones and a torus, a prism of lines, a plate with a
/// hole), each cut by a box through it, the cut and the common its volume;
/// a body with a cavity and a spline prism refused.
#[test]
fn step_solids_are_imported_too() {
    use rusty_occt::step;
    let solids = |name: &str| -> Vec<(rusty_occt::topology::Topology, rusty_occt::Tolerance)> {
        let path = format!("{}/../fixtures/step/{name}.stp", env!("CARGO_MANIFEST_DIR"));
        let text = std::fs::read_to_string(path).unwrap();
        let imported = step::import(&step::read(text.as_bytes()).unwrap()).unwrap();
        imported
            .bodies
            .into_iter()
            .filter(|b| b.item == step::Item::Solid)
            .map(|b| (b.result.unwrap(), b.tolerance))
            .collect()
    };
    for name in [
        "box",
        "cylinder",
        "cylinder_tilted",
        "sphere",
        "hemisphere",
        "cone",
        "frustum",
        "torus",
        "l_prism",
        "plate_hole",
    ] {
        for (k, (t, tol)) in solids(name).into_iter().enumerate() {
            let s = Solid::imported_with(OperationId(91 + k as u64), t, tol)
                .unwrap_or_else(|e| panic!("{name}: {e}"))
                .0;
            let b = s.bounds();
            let mid = |i: usize| (b.min.to_array()[i] + b.max.to_array()[i]) / 2.0;
            // A box over the solid's lower half, a little wider.
            let (lo, hi) = (b.min.to_array(), b.max.to_array());
            let pad = 0.25 * (0..3).map(|i| hi[i] - lo[i]).fold(0.0, f64::max);
            let cutter = Solid::box_at_with(
                OperationId(200),
                Point3::new(lo[0] - pad, lo[1] - pad, lo[2] - pad),
                rusty_occt::Vec3::new(
                    hi[0] - lo[0] + 2.0 * pad,
                    hi[1] - lo[1] + 2.0 * pad,
                    mid(2) - lo[2] + pad + 0.123,
                ),
                tol,
            )
            .unwrap()
            .0;
            let volume =
                |out: Vec<Solid>| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
            let c = volume(
                s.cut(OperationId(201), &cutter)
                    .unwrap_or_else(|e| panic!("{name}: {e}"))
                    .0,
            );
            let m = volume(
                s.common(OperationId(202), &cutter)
                    .unwrap_or_else(|e| panic!("{name}: {e}"))
                    .0,
            );
            let v = s.mass_properties().volume;
            assert!((c + m - v).abs() <= 1e-9 * v, "{name}: {c} + {m} for {v}");
            assert!(c > 0.0 && m > 0.0, "{name}: {c} {m}");
        }
    }
    for (name, why) in [("box_void", "S9e.4b"), ("bspline_prism", "S9f")] {
        for (t, tol) in solids(name) {
            match Solid::imported_with(OperationId(91), t, tol) {
                Err(Error::OutOfDomain(m)) => assert!(m.contains(why), "{name}: {m}"),
                other => panic!("{name}: {:?}", other.map(|_| ())),
            }
        }
    }
}

/// The kernel stores the frames the reference swept, bit for bit (the
/// solids built from rows; the imported ones' are the converter's).
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-imported-frames.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let w: Vec<&str> = row.split('\t').collect();
        let case = cases.iter().find(|c| c.name == w[0]).unwrap();
        let (which, axis) = w[1].split_once(' ').unwrap();
        let k: usize = which.trim_start_matches("solid").parse().unwrap();
        let spec = match k {
            0 => &case.object,
            1 => &case.tool,
            _ => &case.then.as_ref().unwrap().third,
        };
        let frame = protocol::build(spec).frame();
        let v = match axis {
            "n" => frame.normal(),
            "x" => frame.x(),
            _ => frame.y(),
        };
        let got = [v.x, v.y, v.z].map(hex).join(" ");
        assert_eq!(got, w[2], "{} {}", w[0], w[1]);
    }
}
