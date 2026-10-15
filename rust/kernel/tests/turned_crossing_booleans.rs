//! S9e.4b.4d.2c.2: imported turned bodies given to Booleans whose partners
//! cross their smooth edges (the tangent circles where a fillet or a round
//! meets its plane, cylinder or cone), against the independent reference
//! (`fixtures/boolean-turned-crossings-*` from
//! `tools/generate_turned_crossings_boolean_fixtures.py`; the bodies
//! `fixtures/imported/crossing_ribbed.brep`, `OCC485a`'s form on exact data
//! at its resolution, and the earlier steps' `fillet_step`, `fillet_collar`,
//! `cone_fillet_round` and `cone_fillet_fillet`), among them the known
//! failure of "S9e.4b.4d.2c.1a's consistency fixed" (a rib whose side planes
//! cross a round's tangent circle); a box whose side touches a circle
//! refused (a tangency between the inputs). The arrangement's meetings of a
//! plane or a quadric with a torus are found once per torus (two tori of
//! one input, a fillet and a round, each its own). Each group's six
//! Booleans agree with one another (the consistency identities), and the
//! known failure's partners on `consistency_ribbed.brep` itself give the
//! model's commons (`fixtures/boolean-turned-consistency-crossings.tsv`).
#[path = "support/turned_consistency.rs"]
#[allow(dead_code)]
mod consistency;
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Location, Point3, RigidTransform, Solid, Vec3};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-turned-crossings-cases.txt"
    ))
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

/// The reference's kind, reason and `(solids, [volume, area])` per case.
type Expected = BTreeMap<String, (String, String, Option<(usize, [f64; 2])>)>;

fn expected() -> Expected {
    let mut expect = BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-turned-crossings-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let entry = expect
            .entry(name.to_string())
            .or_insert((String::new(), String::new(), None));
        match w[0] {
            "expect" => entry.0 = w[1].to_string(),
            "reason" => entry.1 = row["reason ".len()..].to_string(),
            "result" => {
                let v: Vec<f64> = w[2..4].iter().map(|x| x.parse().unwrap()).collect();
                entry.2 = Some((w[1].parse::<usize>().unwrap(), [v[0], v[1]]));
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

/// Every case as declared: solid within the reference's enclosures, the
/// degenerate refused for their declared reasons.
#[test]
fn every_case_matches_the_reference() {
    let expect = expected();
    let mut failures = Vec::new();
    for (name, run) in runs() {
        let (kind, reason, want) = &expect[name];
        match (kind.as_str(), run) {
            ("degenerate", Err(Error::Degenerate(m))) => {
                if m != reason {
                    failures.push(format!("{name}: refused as {m}, not {reason}"));
                }
                continue;
            }
            ("degenerate", other) => {
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

/// The solid groups: each partner symmetric in `y` (the body about the
/// world's `z`).
const GROUPS: [&str; 8] = [
    "ribbed_post",
    "ribbed_rib",
    "ribbed_box",
    "ribbed_skew",
    "collar_box",
    "step_ring",
    "fillet_rod",
    "round_rib",
];

/// Each group's six Booleans (the common and the fuse both ways, both cuts)
/// valid and consistent within their enclosures: each cut and the common
/// summing to its input, the fuse to both inputs less the common, both
/// orders equal, the moment along `y` zero where the partner is symmetric in
/// `y` (all but the skew cube); the common the reference's.
#[test]
fn each_group_is_consistent() {
    let expect = expected();
    let all = cases();
    let failures: Vec<String> = std::thread::scope(|scope| {
        let handles: Vec<_> = GROUPS
            .iter()
            .map(|group| {
                let case = all
                    .iter()
                    .find(|c| c.name == format!("{group}_common"))
                    .expect("a common case");
                let expect = &expect;
                scope.spawn(move || {
                    let (a, b, _, _) = protocol::run_first(case).unwrap();
                    let (_, [v, _]) = expect[&case.name].2.unwrap();
                    let symmetric = *group != "ribbed_skew";
                    match consistency::consistent(&a, &b, symmetric) {
                        Ok(c) if consistency::meet(c, [v, v], 1e-9 * v) => None,
                        Ok(c) => Some(format!("{group}: the common {c:?} for {v}")),
                        Err(e) => Some(format!("{group}: {e}")),
                    }
                })
            })
            .collect();
        handles
            .into_iter()
            .filter_map(|h| h.join().expect("a worker"))
            .collect()
    });
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The known failure on its own body: `consistency_ribbed.brep` (`OCC485a`'s
/// form, its near-flat bottom and step, its refitted round and its tilted
/// axis) and the rib whose side planes cross its round's tangent circle,
/// the box and the rod across both its fillet and its round, each failing
/// an open Boolean of arcs before the meetings were found per torus: every
/// partner's six Booleans consistent, its common the model's
/// (`tools/generate_turned_consistency_fixtures.py`'s integral).
#[test]
fn the_known_failure_is_consistent_on_its_own_body() {
    let failures = consistency::check_partners(include_str!(
        "../../fixtures/boolean-turned-consistency-crossings.tsv"
    ));
    assert!(failures.is_empty(), "{failures:#?}");
}

/// Each Boolean's history is complete over its inputs, the imported ones'
/// stored ids among them: every entity of every input resolved, and no
/// relation from an id outside the inputs.
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
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    let first: BTreeMap<&str, _> = runs().iter().map(|(n, r)| (n.as_str(), r)).collect();
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

/// Both inputs moved by a dyadic translation keep the reference's volumes
/// (every solid group's common and the known failure's fuse and cut); moved
/// by a rotation that rounds their frames, a partner's face parallel or
/// normal to the axis lies within rounding of it (S9's rule, refused), and
/// the skew cube's three cases evaluate turned with the reference's volumes.
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    let mut names: Vec<String> = GROUPS.iter().map(|g| format!("{g}_common")).collect();
    names.extend(
        [
            "ribbed_rib_fuse",
            "ribbed_rib_cut",
            "ribbed_skew_fuse",
            "ribbed_skew_cut",
        ]
        .map(String::from),
    );
    let failures: Vec<String> = std::thread::scope(|scope| {
        let handles: Vec<_> = names
            .iter()
            .map(|name| {
                let case = all.iter().find(|c| &c.name == name).unwrap();
                let expect = &expect;
                scope.spawn(move || {
                    let mut bad = Vec::new();
                    for (k, motion) in [shift, turn].into_iter().enumerate() {
                        let (a, b, _, _) = protocol::run_first(case).unwrap();
                        let a = a.transform_with(OperationId(801), motion).unwrap().0;
                        let b = b.transform_with(OperationId(802), motion).unwrap().0;
                        let run = match case.op.as_str() {
                            "fuse" => a.fuse(case.operation, &b),
                            "cut" => a.cut(case.operation, &b),
                            _ => a.common(case.operation, &b),
                        };
                        let (count, [v, _]) = expect[name].2.unwrap();
                        let skew = name.starts_with("ribbed_skew");
                        match run {
                            Err(Error::Degenerate(_)) if k == 1 && !skew => {}
                            Err(e) => bad.push(format!("{name} {k}: {e}")),
                            Ok((out, _)) => {
                                let got: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
                                if out.len() != count || (got - v).abs() > 1e-9 * v {
                                    bad.push(format!(
                                        "{name} {k}: {} solids, {got} for {v}",
                                        out.len()
                                    ));
                                }
                            }
                        }
                    }
                    bad
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("a worker"))
            .collect()
    });
    assert!(failures.is_empty(), "{failures:#?}");
}
