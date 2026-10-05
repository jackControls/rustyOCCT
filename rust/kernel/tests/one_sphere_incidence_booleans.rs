//! S9e.4b.3c.2: exact incidences of two pieces of one sphere given to a
//! Boolean (a vertex, a circle or a line of both inputs, plane faces on one
//! plane with overlapping edges: the DRAW survey's `so1` and `so2`, `so2`
//! and `so3`, `so5` and `so2`), against the independent reference
//! (`fixtures/boolean-one-sphere-incidence-*` from
//! `tools/generate_one_sphere_incidence_boolean_fixtures.py`, the bodies
//! under `fixtures/imported/`). Each case runs once (on a few threads) for
//! the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, History, Resolution};
use rusty_occt::identity::{EntityId, OperationId};
use rusty_occt::topology::Slot;
use rusty_occt::{Error, Location, Point3, RigidTransform, Solid, Vec3};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-one-sphere-incidence-cases.txt"
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

/// The reference's kind and `(solids, [volume, area])` per case.
type Expected = BTreeMap<String, (String, Option<(usize, [f64; 2])>)>;

fn expected() -> Expected {
    let mut expect = BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-one-sphere-incidence-expected.tsv")
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
            ("degenerate", other) => {
                failures.push(format!(
                    "{name}: {:?} not refused as degenerate",
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

/// The declared refusals name S9's reasons: the near copies and the two
/// octants touching along their axis two faces within the resolution of one
/// plane, the wedge's axis edge inside the half's face an edge of one input
/// on a face of the other; none refuses as an incidence of S9e.4b.3c.2.
#[test]
fn refusals_name_their_reasons() {
    for (name, run) in runs() {
        let reason = match run {
            Err(Error::Degenerate(m) | Error::OutOfDomain(m)) => *m,
            _ => "",
        };
        if name.starts_with("near_") || name.starts_with("quadrants") {
            assert_eq!(
                reason, "two faces within the resolution of one plane",
                "{name}"
            );
        }
        if name.starts_with("half_wedge") {
            assert_eq!(
                reason, "an edge of one input on a face of the other",
                "{name}"
            );
        }
        assert!(!reason.contains("S9e.4b.3c.2"), "{name}: {reason}");
    }
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

/// What an input entity became in the result.
fn targets(h: &History, id: EntityId) -> Vec<EntityId> {
    match h.resolve(id) {
        Resolution::Same(x) => vec![x],
        Resolution::Split(v) => v,
        Resolution::Merged { into, .. } => vec![into],
        Resolution::Deleted | Resolution::Unknown => Vec::new(),
    }
}

/// The ids of a solid's edges with both ends on the vertical line through
/// `(5, 5)` (a wedge's axis edge), and of its vertices at a point.
fn axis_edges(s: &Solid) -> Vec<EntityId> {
    let t = s.topology();
    let on_axis = |p: Point3| (p.x - 5.0).abs() < 1e-9 && (p.y - 5.0).abs() < 1e-9;
    t.edges()
        .iter()
        .enumerate()
        .filter(|(_, e)| {
            [e.start, e.end]
                .iter()
                .all(|v| v.is_some_and(|v| on_axis(t.vertices()[v.index()].position)))
        })
        .filter_map(|(k, _)| t.id_of(Slot::Edge(rusty_occt::topology::EdgeId::new(k))))
        .collect()
}

fn vertex_at(s: &Solid, at: Point3) -> Option<EntityId> {
    let t = s.topology();
    let k = t
        .vertices()
        .iter()
        .position(|v| (v.position - at).length() < 1e-9)?;
    t.id_of(Slot::Vertex(rusty_occt::topology::VertexId::new(k)))
}

/// An edge or a vertex of both inputs continues both: the two wedges'
/// common (`so2` and `so3`) has one axis edge, which both axis edges
/// become, and one corner and one pole, which both inputs' become; the
/// higher wedge's (`so5`'s) axis edge inside the other's becomes the
/// common's.
#[test]
fn entities_of_both_continue_both() {
    let all: BTreeMap<&str, _> = runs().iter().map(|(n, r)| (n.as_str(), r)).collect();
    let Ok((a, b, out, h)) = all["wedge_wedge_common"] else {
        panic!("the wedges' common")
    };
    let [result] = &out[..] else {
        panic!("one solid")
    };
    let axis = axis_edges(result);
    assert_eq!(axis.len(), 1);
    for s in [a, b] {
        let mine = axis_edges(s);
        assert_eq!(mine.len(), 1);
        assert_eq!(targets(h, mine[0]), axis);
        for at in [Point3::new(5.0, 5.0, 4.0), Point3::new(5.0, 5.0, 9.0)] {
            let v = vertex_at(s, at).expect("a corner and a pole");
            assert_eq!(targets(h, v), vec![vertex_at(result, at).unwrap()]);
        }
    }
    let Ok((a, _, out, h)) = all["high_wedge_common"] else {
        panic!("the higher wedge's common")
    };
    let [result] = &out[..] else {
        panic!("one solid")
    };
    let mine = axis_edges(a);
    assert_eq!(mine.len(), 1);
    assert_eq!(targets(h, mine[0]), axis_edges(result));
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

/// Both inputs moved by exact motions (a dyadic translation, a quarter turn
/// about a world axis: the stored zeros kept) keep the reference's volumes;
/// turned by a rotation that rounds their frames, the incidences are within
/// the resolution and no longer exact: refused as `Degenerate`, or within
/// the reference.
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let quarter = RigidTransform::quarter_turn(Point3::new(1.0, 0.0, 0.0), 0, 1).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    for name in [
        "hemi_wedge_cut",
        "wedge_wedge_fuse",
        "high_wedge_cut",
        "half_octant_common",
        "x_wedges_cut",
        "cap_wedge_fuse",
    ] {
        let case = all.iter().find(|c| c.name == name).unwrap();
        for (k, motion) in [shift, quarter, turn].into_iter().enumerate() {
            let (a, b, _, _) = protocol::run_first(case).unwrap();
            let a = a.transform_with(OperationId(801), motion).unwrap().0;
            let b = b.transform_with(OperationId(802), motion).unwrap().0;
            let run = match case.op.as_str() {
                "fuse" => a.fuse(case.operation, &b),
                "cut" => a.cut(case.operation, &b),
                _ => a.common(case.operation, &b),
            };
            let (count, [v, _]) = expect[name].1.unwrap();
            match run {
                Err(Error::Degenerate(_)) if k == 2 => continue,
                Err(e) => panic!("{name} {k}: {e}"),
                Ok((out, _)) => {
                    let got: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
                    assert_eq!(out.len(), count, "{name} {k}");
                    assert!((got - v).abs() <= 1e-9 * v, "{name} {k}: {got} for {v}");
                }
            }
        }
    }
}
