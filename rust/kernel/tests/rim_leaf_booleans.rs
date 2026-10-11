//! S9e.4b.4b.2b.2b.1: imported bodies of several cylinder faces and plane
//! faces that are a prism leaf whose cap a pocket or another part cuts at its
//! rim (its walls ending below it, a hole a well on the other part), or a
//! tooth standing in such a pocket to the leaf's cap (a stadium plate less a
//! square with a round hole across its rim, a stadium plate fused with a
//! square with a square hole over its arc end) given to Booleans, against the
//! independent reference
//! (`fixtures/boolean-rim-leaves-*` from
//! `tools/generate_rim_leaves_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); the `boolean` fuzz target's `RIMS` bodies built,
//! written by its writer, read back and imported. Each case runs once (on a
//! few threads) for the checks that read its result.
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
    protocol::cases(include_str!("../../fixtures/boolean-rim-leaves-cases.txt"))
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
    for line in include_str!("../../fixtures/boolean-rim-leaves-expected.tsv")
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
        "post_rod_cut",
        "well_rod_common",
        "well_rod_fuse",
        "rod_post_cut",
        "post_well_fuse",
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

/// Every body imports as its leaf cut at its rim with its pocket and post or
/// its plate and its ring with a well (its model matched on import), its
/// volume its closed form, its stored vertices on its boundary, points in its
/// material inside and points outside it outside; S9e.4b.4b.2a's notch (both
/// caps cut at its rim) is refused as S9e.4b.4b.2b.2b.3's.
#[test]
fn imported_bodies_are_their_leaves() {
    use std::f64::consts::PI;
    let post = frame([2.0, 1.0, 0.5], [0.0, 0.0, 1.0], [8.0, 15.0, 0.0]);
    let well = frame([1.0, 0.5, 0.75], [0.0, 0.0, 1.0], [7.0, 24.0, 0.0]);
    let at = |f: Frame3, u: f64, v: f64, w: f64| f.point(Point2::new(u, v), w).to_array();
    type Points = Vec<[f64; 3]>;
    let bodies: [(&str, f64, Points, Points); 2] = [
        (
            "rim_post",
            60.0 + 12.4375 * PI,
            vec![
                at(post, 4.0, 0.0, 1.5),
                at(post, -0.5, 0.25, 2.0),
                at(post, -1.5, -1.5, 0.5),
            ],
            vec![at(post, 0.0, 1.5, 2.0), at(post, 8.0, 0.0, 2.75)],
        ),
        (
            "rim_well",
            93.0 + 9.375 * PI + 6.25 * 0.6f64.asin(),
            vec![at(well, 8.0, 0.0, 2.5), at(well, 6.0, 0.0, 1.5)],
            vec![at(well, 6.0, 0.0, 2.5), at(well, 9.5, 0.0, 0.5)],
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
             two prisms', a flat or a pocket (S9e.4b.4b.2b.2b.3b)"
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

/// The `boolean` fuzz target's `RIMS` bodies (`rims_body`, `k = s / 8`): a
/// stadium plate of half discs of radius `3k` about `(0, 0)` and `(4k, 0)`
/// over `[0, h]` less the square `[-6k, 2k] x [-6k, 6k]` with a hole of
/// radius `3k / 2` about `(-k, k / 2)` over `[h / 2, 2h]`; or a stadium plate
/// of half discs of radius `5k / 2` about `(0, 0)` and `(4k, 0)` over `[0, h]`
/// on the frame turned about its axis (its `x` along `8 x + 15 y`) fused with
/// the square `[5k / 4, 9k] x [-3k / 2, 3k / 2]` with the hole `[7k /
/// 2, 11k / 2] x [-k / 2, k / 2]` over `[h / 2, 3h / 2]` (its long walls
/// across the plate's arc at `x = 6k`).
fn rims_body(f: Frame3, s: f64, h: f64, well: bool) -> Solid {
    let k = s / 8.0;
    let tol = tolerance();
    let stadium = |cx: f64, cy: f64, length: f64, r: f64| {
        let arc = |x: f64| Segment::Arc {
            center: Point2::new(x, cy),
            radius: r,
            ccw: true,
        };
        let points = [
            (cx, cy - r),
            (cx + length, cy - r),
            (cx + length, cy + r),
            (cx, cy + r),
        ]
        .map(|(x, y)| Point2::new(x, y))
        .to_vec();
        Boundary::path(
            points,
            vec![Segment::Line, arc(cx + length), Segment::Line, arc(cx)],
            tol,
        )
        .unwrap()
    };
    let square = |x0: f64, y0: f64, x1: f64, y1: f64| {
        let points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec();
        Boundary::polygon(points, tol).unwrap()
    };
    let prism = |outer: Boundary, holes: Vec<Boundary>, h0: f64, h1: f64, op: u64| {
        Solid::extrude_with(
            OperationId(op),
            Profile::new(outer, holes, tol).unwrap(),
            f,
            h0,
            h1,
        )
        .unwrap()
        .0
    };
    let one = |r: rusty_occt::Result<(Vec<Solid>, History)>| {
        let (out, _) = r.unwrap();
        assert_eq!(out.len(), 1);
        out.into_iter().next().unwrap()
    };
    let plate = |r: f64| prism(stadium(0.0, 0.0, 4.0 * k, r), vec![], 0.0, h, 31);
    if well {
        let turned = Frame3::new(f.origin(), f.normal(), f.x() * 8.0 + f.y() * 15.0, tol).unwrap();
        let prism = |outer: Boundary, holes: Vec<Boundary>, h0: f64, h1: f64, op: u64| {
            let profile = Profile::new(outer, holes, tol).unwrap();
            Solid::extrude_with(OperationId(op), profile, turned, h0, h1)
                .unwrap()
                .0
        };
        let plate = prism(stadium(0.0, 0.0, 4.0 * k, 2.5 * k), vec![], 0.0, h, 31);
        let hole = square(3.5 * k, -0.5 * k, 5.5 * k, 0.5 * k);
        let outer = square(1.25 * k, -1.5 * k, 9.0 * k, 1.5 * k);
        let ring = prism(outer, vec![hole], h / 2.0, 1.5 * h, 32);
        one(plate.fuse(OperationId(33), &ring))
    } else {
        let hole = Boundary::circle(Point2::new(-k, 0.5 * k), 1.5 * k, tol).unwrap();
        let outer = square(-6.0 * k, -6.0 * k, 2.0 * k, 6.0 * k);
        let tool = prism(outer, vec![hole], h / 2.0, 2.0 * h, 32);
        one(plate(3.0 * k).cut(OperationId(33), &tool))
    }
}

/// The kernel's own bodies, written, read back and imported, each its leaf
/// cut at its rim with its post or its plate and ring with a well: its
/// volume, its stored vertices on its boundary, and its Booleans with a
/// turned box the kernel's own results'. The `boolean` fuzz target's `RIMS`
/// bodies on the world's axes and on its tilted frame.
#[test]
fn kernel_bodies_imported_as_leaves() {
    let tilted = frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let probe_at = |at: Point3| {
        // Turned off every body's axis (a probe along the tilted frame's
        // axis lies within rounding of the import's leaves, on their stored
        // caps' frames, where the kernel's own body's are on that frame).
        let f = Frame3::new(
            at,
            Vec3::new(2.0, 3.0, 6.0),
            Vec3::new(3.0, -6.0, 2.0),
            tolerance(),
        )
        .unwrap();
        let points = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
            .map(|(x, y)| Point2::new(x, y))
            .to_vec();
        let square = Boundary::polygon(points, tolerance()).unwrap();
        Solid::extrude_with(
            OperationId(16),
            Profile::new(square, vec![], tolerance()).unwrap(),
            f,
            -3.0,
            3.0,
        )
        .unwrap()
        .0
    };
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    let near = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(1.0);
    let world = Frame3::xy();
    let bodies = [
        (
            "post",
            rims_body(world, 2.0, 1.0, false),
            Point3::new(0.25, 0.5, 0.75),
        ),
        (
            "well",
            rims_body(world, 2.0, 1.0, true),
            Point3::new(1.25, 0.25, 1.0),
        ),
        (
            "tilted post",
            rims_body(tilted, 3.0, 1.25, false),
            tilted.point(Point2::new(0.25, 0.5), 0.75),
        ),
        (
            "tilted well",
            rims_body(tilted, 3.0, 1.25, true),
            tilted.point(Point2::new(1.5, 0.25), 1.0),
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
