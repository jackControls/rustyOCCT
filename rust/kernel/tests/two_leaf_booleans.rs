//! S9e.4b.4b.2b.2a: imported bodies of several cylinder faces and plane faces
//! that are two prism leaves or a lone leaf with a pocket (a box with a
//! stadium through it, a stadium with a box through it across its axis, a stadium plate
//! less a turned box from its top) given to Booleans, against the
//! independent reference (`fixtures/boolean-two-leaves-*` from
//! `tools/generate_two_leaves_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); the kernel's own through and the `boolean` fuzz
//! target's `LEAVES` bodies built, written by its writer, read back and
//! imported. Each case runs once (on a few threads) for the checks that read
//! its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, History, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Location, Point2, Point3, Profile, RigidTransform, Segment, Solid,
    Tolerance, Vec3,
};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!("../../fixtures/boolean-two-leaves-cases.txt"))
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
    for line in include_str!("../../fixtures/boolean-two-leaves-expected.tsv")
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
/// degenerate and unsupported refused for their declared reasons.
#[test]
fn every_case_matches_the_reference() {
    let expect = expected();
    let mut failures = Vec::new();
    for (name, run) in runs() {
        let (kind, reason, want) = &expect[name];
        match (kind.as_str(), run) {
            ("degenerate", Err(Error::Degenerate(m)))
            | ("unsupported", Err(Error::OutOfDomain(m))) => {
                if m != reason {
                    failures.push(format!("{name}: refused as {m}, not {reason}"));
                }
                continue;
            }
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

/// Both inputs moved by a dyadic translation and by a rotation that rounds
/// their frames keep the reference's volumes (a turn moving an exact
/// incidence within rounding is refused as `Degenerate`).
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    for name in [
        "through_rod_cut",
        "cross_rod_common",
        "pocketed_rod_fuse",
        "rod_cross_cut",
        "cross_pocketed_fuse",
    ] {
        let case = all.iter().find(|c| c.name == name).unwrap();
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
            match run {
                Err(Error::Degenerate(_)) if k == 1 => continue,
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

fn body(name: &str, op: u64) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(op), topology, resolution).map(|(s, _)| s)
}

fn tolerance() -> Tolerance {
    Tolerance::new(1e-7, 1e-12).unwrap()
}

fn frame(o: [f64; 3], n: [f64; 3], x: [f64; 3]) -> Frame3 {
    let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
    Frame3::new(Point3::new(o[0], o[1], o[2]), v(n), v(x), tolerance()).unwrap()
}

/// A stadium of half discs of radius `r` about `(x0, y)` and `(x1, y)`,
/// counter-clockwise.
fn stadium(x0: f64, x1: f64, y: f64, r: f64) -> Boundary {
    let arc = |x: f64| Segment::Arc {
        center: Point2::new(x, y),
        radius: r,
        ccw: true,
    };
    let points = [(x0, y - r), (x1, y - r), (x1, y + r), (x0, y + r)]
        .map(|(x, y)| Point2::new(x, y))
        .to_vec();
    Boundary::path(
        points,
        vec![Segment::Line, arc(x1), Segment::Line, arc(x0)],
        tolerance(),
    )
    .unwrap()
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Boundary {
    let points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        .map(|(x, y)| Point2::new(x, y))
        .to_vec();
    Boundary::polygon(points, tolerance()).unwrap()
}

fn prism(b: Boundary, f: Frame3, h0: f64, h1: f64, op: u64) -> Solid {
    Solid::extrude_with(
        OperationId(op),
        Profile::new(b, vec![], tolerance()).unwrap(),
        f,
        h0,
        h1,
    )
    .unwrap()
    .0
}

/// Every body imports as its two leaves or its lone leaf with its pocket
/// (its model matched on import), its volume its closed form, its stored
/// vertices on its boundary, points in its material inside and points
/// outside it outside; S9e.4b.4b.2a's notch (both caps cut at its rim) is
/// refused as S9e.4b.4b.2b.2b.2's.
#[test]
fn imported_bodies_are_their_leaves() {
    use std::f64::consts::PI;
    let through = frame([2.0, 1.0, 0.5], [0.0, 0.0, 1.0], [7.0, 24.0, 0.0]);
    let world = frame([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let pocketed = frame([7.0, -1.5, 0.75], [0.0, 0.0, 1.0], [8.0, 15.0, 0.0]);
    let at = |f: Frame3, u: f64, v: f64, w: f64| f.point(Point2::new(u, v), w).to_array();
    type Points = Vec<[f64; 3]>;
    let bodies: [(&str, f64, Points, Points); 3] = [
        (
            "two_through",
            120.0 + 4.5 * PI,
            vec![at(through, 1.0, 1.0, 1.0), at(through, 4.0, 3.0, -0.5)],
            vec![at(through, 1.0, 1.0, -0.5), at(through, 4.0, 0.5, 2.5)],
        ),
        (
            "two_cross",
            64.0 + 32.0 * PI / 3.0 - 2.0 * 3.0f64.sqrt(),
            vec![at(world, 5.0, -1.5, 2.5), at(world, -3.0, 0.0, 1.5)],
            vec![at(world, -3.0, 0.0, 2.5), at(world, 7.0, 1.5, 1.5)],
        ),
        (
            "two_pocketed",
            61.5 + 15.625 * PI,
            vec![at(pocketed, 0.5, 0.0, 0.5), at(pocketed, 0.0, 0.0, 2.0)],
            vec![at(pocketed, 3.0, 0.0, 2.0), at(pocketed, 3.0, 0.0, 2.75)],
        ),
    ];
    for (name, volume, inside, outside) in bodies {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let m = s.mass_properties();
        assert!(
            (m.volume - volume).abs() <= 1e-9 * volume,
            "{name}: {} for {volume}",
            m.volume
        );
        for (points, want) in [(inside, Location::Inside), (outside, Location::Outside)] {
            for p in points {
                assert_eq!(
                    s.classify(Point3::new(p[0], p[1], p[2])).unwrap(),
                    want,
                    "{name} {p:?}"
                );
            }
        }
        for v in s.topology().vertices() {
            assert_eq!(
                s.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
    }
    match body("leaf_notch", 91) {
        Err(Error::OutOfDomain(m)) => assert_eq!(
            m,
            "an imported body of several primitives with plane faces other than their ends, one or \
             two prisms', a flat or a pocket (S9e.4b.4b.2b.2b.2)"
        ),
        other => panic!("leaf_notch: {:?}", other.map(|_| ())),
    }
}

/// A solid written by the kernel's writer, read back and imported.
fn reimported(s: &Solid, op: u64) -> Result<Solid, Error> {
    use rusty_occt::occt_brep::{import, read, write};
    let text = write(s.topology(), s.resolution().linear()).unwrap();
    let doc = read(&text).unwrap();
    let [solid] = <[_; 1]>::try_from(import(&doc).solids).unwrap();
    Solid::imported_with(OperationId(op), solid.result.unwrap(), solid.tolerance).map(|(s, _)| s)
}

/// The `boolean` fuzz target's `LEAVES` bodies (`leaves_body`, `k = s /
/// 8`): the box of half side `5k` over `[0, h]` fused with a stadium of half
/// discs of radius `2k` about `(-k, k)` and `(k, k)` over `[-h / 2, 3h /
/// 2]`, or a stadium of half discs of radius `3k` about `(-2k, 0)` and `(2k,
/// 0)` over `[0, h]` less a square of half side `3k / 2` turned about its
/// axis over `[h / 2, 2h]`.
fn leaves_body(f: Frame3, s: f64, h: f64, pocket: bool) -> Solid {
    let k = s / 8.0;
    let one = |r: rusty_occt::Result<(Vec<Solid>, History)>| {
        let (out, _) = r.unwrap();
        assert_eq!(out.len(), 1);
        out.into_iter().next().unwrap()
    };
    if pocket {
        let plate = prism(stadium(-2.0 * k, 2.0 * k, 0.0, 3.0 * k), f, 0.0, h, 28);
        let turned = Frame3::new(
            f.origin(),
            f.normal(),
            f.x() * 3.0 + f.y() * 4.0,
            tolerance(),
        )
        .unwrap();
        let a = 1.5 * k;
        let hole = prism(square(-a, -a, a, a), turned, h / 2.0, 2.0 * h, 29);
        one(plate.cut(OperationId(30), &hole))
    } else {
        let a = 5.0 * k;
        let block = prism(square(-a, -a, a, a), f, 0.0, h, 28);
        let rod = prism(stadium(-k, k, k, 2.0 * k), f, -h / 2.0, 1.5 * h, 29);
        one(block.fuse(OperationId(30), &rod))
    }
}

/// The kernel's own bodies, written, read back and imported, each its two
/// leaves or its lone leaf with its pocket: its volume, its stored vertices
/// on its boundary, and its Booleans with a turned box the kernel's own
/// results'. The through on a turned frame (as the fixture's: the
/// kernel's writer stores each piece of the stadium's walls on a cylinder of
/// its own, coaxial), and the `boolean` fuzz target's `LEAVES` bodies on the
/// world's axes and on its tilted frame (the measurement's first results of
/// a box fused with a stadium through it and of a stadium plate less the
/// chained turned box; the kernel refuses its own pocketed plate on the
/// tilted frame, the turned square's walls within rounding of its axis, so
/// the stage takes it on the world's axes alone).
#[test]
fn kernel_bodies_imported_as_leaves() {
    let one = |r: rusty_occt::Result<(Vec<Solid>, History)>| {
        let (out, _) = r.unwrap();
        assert_eq!(out.len(), 1);
        out.into_iter().next().unwrap()
    };
    let skew2 = frame([3.0, 4.5, 5.0], [4.0, 4.0, 7.0], [1.0, -8.0, 4.0]);
    let tilted = frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let through = one(prism(square(0.0, 0.0, 8.0, 6.0), skew2, 0.0, 2.0, 1).fuse(
        OperationId(2),
        &prism(stadium(2.0, 6.0, 3.0, 1.5), skew2, -1.0, 3.0, 3),
    ));
    let probe_at = |at: Point3| {
        // Turned off every body's axis (a probe along the tilted frame's
        // axis lies within rounding of the import's leaves, on their stored
        // cap's frame, where the kernel's own body's are on that frame).
        let f = Frame3::new(
            at,
            Vec3::new(2.0, 3.0, 6.0),
            Vec3::new(3.0, -6.0, 2.0),
            tolerance(),
        )
        .unwrap();
        prism(square(-1.0, -1.0, 1.0, 1.0), f, -3.0, 3.0, 16)
    };
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    let near = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(1.0);
    let bodies = [
        ("through", through, skew2.point(Point2::new(6.5, 3.5), 1.5)),
        (
            "leaves",
            leaves_body(Frame3::xy(), 2.0, 1.0, false),
            Point3::new(0.75, 0.5, 1.25),
        ),
        (
            "tilted leaves",
            leaves_body(tilted, 3.0, 1.25, false),
            tilted.point(Point2::new(0.75, 0.5), 1.25),
        ),
        (
            "pocket",
            leaves_body(Frame3::xy(), 2.0, 1.0, true),
            Point3::new(0.5, 0.25, 0.75),
        ),
    ];
    for (name, s, at) in bodies {
        let i = reimported(&s, 50).unwrap_or_else(|e| panic!("{name}: {e}"));
        let (v0, v1) = (s.mass_properties().volume, i.mass_properties().volume);
        assert!(near(v0, v1), "{name}: {v1} for {v0}");
        for v in i.topology().vertices() {
            assert_eq!(
                i.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
        let probe = probe_at(at);
        for (direct, again) in [
            (
                s.fuse(OperationId(60), &probe),
                i.fuse(OperationId(61), &probe),
            ),
            (
                s.cut(OperationId(60), &probe),
                i.cut(OperationId(61), &probe),
            ),
            (
                s.common(OperationId(60), &probe),
                i.common(OperationId(61), &probe),
            ),
        ] {
            match (direct, again) {
                (Ok(d), Ok(a)) => {
                    let (d, a) = (volume(Ok(d)), volume(Ok(a)));
                    assert!(near(d, a), "{name}: {a} for {d}");
                }
                (Err(d), Err(a)) => assert_eq!(d, a, "{name}"),
                (d, a) => panic!("{name}: {:?} and {:?}", d.map(|_| ()), a.map(|_| ())),
            }
        }
    }
}
