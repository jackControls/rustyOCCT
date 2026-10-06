//! S9e.4b.4b.2a: imported bodies of several sphere, cylinder and cone faces
//! led by a prism leaf (a prism of lines and arcs with both its caps, its
//! walls' joints tangent where its corners are rounded: `bfuse_complex/K1`'s
//! rounded box less a bore, a stadium plate with a boss, a dimple, a conical
//! pocket, a dome) given to Booleans, against the independent reference
//! (`fixtures/boolean-prism-leaves-*` from
//! `tools/generate_prism_leaves_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); the kernel's own rounded plates with a bore and a
//! boss built, written by its writer, read back and imported. Each case runs
//! once (on a few threads) for the checks that read its result.
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
    protocol::cases(include_str!(
        "../../fixtures/boolean-prism-leaves-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-prism-leaves-expected.tsv")
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
        "bore_rod_cut",
        "boss_box_common",
        "dimple_ball_fuse",
        "pocket_rod_cut",
        "dome_box_cut",
        "boss_pocket_common",
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

/// A ball's cap beyond a plane at distance `d` from its centre.
fn cap(r: f64, d: f64) -> f64 {
    let t = r - d;
    std::f64::consts::PI * t * t * (3.0 * r - t) / 3.0
}

/// Every body imports as its chain led by its prism (its model matched on
/// import), its volume its closed form, its stored vertices on its boundary,
/// points in its material inside and points its chain removes outside; the
/// post (on a fillet's circle) is refused for its own tangency and the notch
/// (both caps cut at its rim: no prism leaf) as S9e.4b.4b.2b's.
#[test]
fn imported_bodies_are_their_chains() {
    use std::f64::consts::PI;
    let skew = |o: [f64; 3]| frame(o, [8.0, 4.0, 1.0], [-1.0, 4.0, -8.0]);
    let skew2 = |o: [f64; 3]| frame(o, [4.0, 4.0, 7.0], [1.0, -8.0, 4.0]);
    let bore = frame([1.0, 2.0, 0.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let boss = skew2([3.0, 4.5, 5.0]);
    let pocket = skew(boss.point(Point2::new(1.5, -0.5), -1.25).to_array());
    let dome = skew2([1.0, -1.0, 0.5]);
    let at = |f: Frame3, u: f64, v: f64, w: f64| f.point(Point2::new(u, v), w).to_array();
    // The rounded square of side 9 (its corners' arcs of radius 2) and the
    // stadium's areas.
    let square = 65.0 + 4.0 * PI;
    let stadium = 30.0 + 6.25 * PI;
    // Each body's name, closed-form volume, and points inside and outside.
    type Points = Vec<[f64; 3]>;
    let bodies: [(&str, f64, Points, Points); 5] = [
        (
            "leaf_bore",
            6.0 * square - 20.25 * PI,
            vec![at(bore, 1.0, 4.5, 3.0), at(bore, 4.5, 4.5, 0.5)],
            vec![at(bore, 4.5, 4.5, 3.0), at(bore, 0.25, 0.25, 3.0)],
        ),
        (
            "leaf_boss",
            2.0 * stadium + 3.0 * PI * 1.5625,
            vec![at(boss, 3.0, 0.0, 4.5), at(boss, -2.0, 0.0, 1.0)],
            vec![at(boss, 5.0, 0.0, 3.0), at(boss, 3.0, 0.0, 5.5)],
        ),
        (
            "leaf_dimple",
            3.0 * square - cap(2.5, 1.5),
            vec![[4.5, 4.5, 1.5], [1.0, 4.5, 2.5]],
            vec![[4.5, 4.5, 2.5], [0.25, 0.25, 1.0]],
        ),
        (
            "leaf_pocket",
            3.0 * stadium - PI * 2.0 / 3.0 * 3.25,
            vec![at(pocket, 3.0, 0.0, 0.5), at(pocket, 3.0, 1.75, 2.5)],
            vec![at(pocket, 3.0, 0.0, 2.0), at(pocket, 3.0, 2.75, 1.5)],
        ),
        (
            "leaf_dome",
            3.0 * square + cap(2.0, 0.5),
            vec![at(dome, 4.5, 4.5, 4.0), at(dome, 1.0, 4.5, 1.5)],
            vec![at(dome, 4.5, 4.5, 4.75), at(dome, 6.5, 4.5, 3.5)],
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
    match body("leaf_post", 91) {
        Err(Error::Degenerate(m)) => assert_eq!(
            m,
            "an imported body of several primitives whose faces are tangent along an edge"
        ),
        other => panic!("leaf_post: {:?}", other.map(|_| ())),
    }
    match body("leaf_notch", 91) {
        Err(Error::OutOfDomain(m)) => assert_eq!(
            m,
            "an imported body of several primitives with plane faces other than their ends or a \
             prism's (S9e.4b.4b.2b)"
        ),
        other => panic!("leaf_notch: {:?}", other.map(|_| ())),
    }
}

/// A rounded square `[0, side]^2`, its corners rounded by arcs of radius `r`.
fn rounded(side: f64, r: f64) -> Profile {
    let arc = |x: f64, y: f64| Segment::Arc {
        center: Point2::new(x, y),
        radius: r,
        ccw: true,
    };
    let s = side;
    let points = [
        (r, 0.0),
        (s - r, 0.0),
        (s, r),
        (s, s - r),
        (s - r, s),
        (r, s),
        (0.0, s - r),
        (0.0, r),
    ]
    .map(|(x, y)| Point2::new(x, y))
    .to_vec();
    let segments = vec![
        Segment::Line,
        arc(s - r, r),
        Segment::Line,
        arc(s - r, s - r),
        Segment::Line,
        arc(r, s - r),
        Segment::Line,
        arc(r, r),
    ];
    let outer = Boundary::path(points, segments, tolerance()).unwrap();
    Profile::new(outer, Vec::new(), tolerance()).unwrap()
}

/// A cylinder of radius `r` about `(cx, cy)` on a frame over `[h0, h1]`.
fn rod(f: Frame3, cx: f64, cy: f64, r: f64, h0: f64, h1: f64, op: u64) -> Solid {
    let profile = Profile::new(
        Boundary::circle(Point2::new(cx, cy), r, tolerance()).unwrap(),
        vec![],
        tolerance(),
    )
    .unwrap();
    Solid::extrude_with(OperationId(op), profile, f, h0, h1)
        .unwrap()
        .0
}

/// A solid written by the kernel's writer, read back and imported.
fn reimported(s: &Solid, op: u64) -> Solid {
    use rusty_occt::occt_brep::{import, read, write};
    let text = write(s.topology(), s.resolution().linear()).unwrap();
    let doc = read(&text).unwrap();
    let [solid] = <[_; 1]>::try_from(import(&doc).solids).unwrap();
    Solid::imported_with(OperationId(op), solid.result.unwrap(), solid.tolerance)
        .unwrap()
        .0
}

/// The kernel's own rounded plates on the world's axes, written, read back
/// and imported: a boss (a rounded plate fused with a rod along its axis from
/// inside it to a cap), a dimple and a dome (less and fused with a ball
/// meeting its top cap). Each its chain led by its prism: its volume, and its
/// Booleans with a turned rod the kernel's own results'. K1's shape (a rounded
/// plate less a rod along `y` through two walls) the writer refuses: it
/// writes circles alone, and the kernel stores the bore's rims as ellipses of
/// equal axes (S9e.4b.4b.1's amendment (a)).
#[test]
fn kernel_bodies_imported_as_chains() {
    let one = |r: rusty_occt::Result<(Vec<Solid>, History)>| {
        let (out, _) = r.unwrap();
        assert_eq!(out.len(), 1);
        out.into_iter().next().unwrap()
    };
    let world = frame([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let along_y = frame([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let plate = |h: f64, op: u64| {
        Solid::extrude_with(OperationId(op), rounded(9.0, 2.0), world, 0.0, h)
            .unwrap()
            .0
    };
    let bored = one(plate(6.0, 1).cut(OperationId(2), &rod(along_y, 3.0, 4.5, 1.5, -1.0, 10.0, 3)));
    let ball = |c: [f64; 3], r: f64, op: u64| {
        let half = std::f64::consts::FRAC_PI_2;
        let f = frame(c, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        Solid::sphere_with(OperationId(op), f, r, -half, half, tolerance())
            .unwrap()
            .0
    };
    let dimpled = one(plate(3.0, 7).cut(OperationId(8), &ball([4.5, 4.5, 4.5], 2.5, 9)));
    let domed = one(plate(3.0, 10).fuse(OperationId(11), &ball([4.5, 4.5, 2.5], 2.0, 12)));
    let bossed = one(plate(2.0, 4).fuse(OperationId(5), &rod(world, 3.25, 4.5, 1.25, 1.0, 5.0, 6)));
    let probe = rod(
        frame([3.5, 2.25, 1.75], [0.0, 3.0, 4.0], [3.0, 0.0, 4.0]),
        0.0,
        0.0,
        1.25,
        -5.0,
        5.0,
        14,
    );
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    // K1's shape, its bore's rims written as circles of equal axes, is
    // refused by the writer.
    assert!(rusty_occt::occt_brep::write(bored.topology(), bored.resolution().linear()).is_err());
    for (name, s) in [("bossed", bossed), ("dimpled", dimpled), ("domed", domed)] {
        let i = reimported(&s, 50);
        let (v0, v1) = (s.mass_properties().volume, i.mass_properties().volume);
        assert!((v0 - v1).abs() <= 1e-9 * v0, "{name}: {v1} for {v0}");
        for v in i.topology().vertices() {
            assert_eq!(
                i.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
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
                    assert!((d - a).abs() <= 1e-9 * d.max(1.0), "{name}: {a} for {d}");
                }
                (Err(d), Err(a)) => assert_eq!(d, a, "{name}"),
                (d, a) => panic!("{name}: {:?} and {:?}", d.map(|_| ()), a.map(|_| ())),
            }
        }
    }
}
