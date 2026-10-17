//! S9e.4b.4c.2b.2: imported bodies of several primitives that are a Boolean
//! chain only in the inside-first order (a rod less a wider ball, its
//! dimple through the rod's rim or inside it; a turned flange whose spigot
//! stands in its recess, with a lug) given to Booleans, against the
//! independent reference (`fixtures/boolean-inside-first-*` from
//! `tools/generate_inside_first_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); each group's six Booleans consistent with one
//! another; the kernel's own bodies of the dataset's two forms written by
//! its writer, read back and imported. Each case runs once (on a few
//! threads) for the checks that read its result.
#[path = "support/turned_consistency.rs"]
#[allow(dead_code)]
mod consistency;
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, History, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Location, Point2, Point3, Profile, RigidTransform, Solid, Tolerance,
    Vec3,
};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-inside-first-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-inside-first-expected.tsv")
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
        "dimple_rod_cut",
        "cup_box_common",
        "flange_rod_fuse",
        "rod_dimple_cut",
        "flange_cyl_fuse",
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

/// The groups of one Boolean whose six Booleans are checked.
const GROUPS: [&str; 6] = [
    "dimple_rod",
    "cup_box",
    "flange_rod",
    "flange_box",
    "rod_dimple",
    "flange_cyl",
];

/// Each group's six Booleans (the common and the fuse both ways, both cuts)
/// valid at the larger of the inputs' resolutions and consistent within
/// their enclosures: each cut and the common summing to its input, the fuse
/// to both inputs less the common, both orders equal; the common the
/// reference's. (A Boolean whose object is a construction and whose tool the
/// imported body takes the construction's smaller resolution and keeps the
/// body's stored entities, valid at the body's: every imported body of
/// several primitives alike, the earlier steps' leaves among them.)
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
                    let at = if a.resolution().linear() >= b.resolution().linear() {
                        a.resolution()
                    } else {
                        b.resolution()
                    };
                    match consistency::consistent_at(&a, &b, false, Some(at)) {
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

/// Every body imports as its chain in the inside-first order (its model
/// matched on import): its volume the closed form, points in its material
/// inside (the rod above its dimple, the spigot, the lug), points in its
/// dimple, its recess and beside its lug outside, its stored vertices on its
/// boundary.
#[test]
fn imported_bodies_are_their_inside_first_chains() {
    use std::f64::consts::PI;
    let at = |o: [f64; 3], x: f64, y: f64, z: f64| [o[0] + x, o[1] + y, o[2] + z];
    let (dimple, cup, flange) = ([1.0, 2.0, 0.5], [-1.0, 0.5, 1.0], [3.0, 3.0, 0.5]);
    let lug = 6.25 * 0.75f64.atan() - 0.75;
    type Points = Vec<[f64; 3]>;
    let bodies: [(&str, f64, Points, Points); 3] = [
        (
            "inside_dimple",
            54.0 * PI - 14.0 * PI / 3.0,
            vec![at(dimple, 0.0, 0.0, 1.5), at(dimple, 2.5, 0.0, 1.0)],
            vec![at(dimple, 0.0, 0.0, 0.5), at(dimple, 0.0, 0.0, -0.25)],
        ),
        (
            "inside_cup",
            54.0 * PI - 29.0 * PI / 24.0,
            vec![at(cup, 0.0, 0.0, 1.0), at(cup, 2.5, 0.0, 0.125)],
            vec![at(cup, 0.0, 0.0, 0.25), at(cup, 0.0, 0.0, -0.25)],
        ),
        (
            "inside_flange",
            47.0 * PI / 2.0 + lug,
            vec![
                at(flange, 0.0, 0.0, -0.25),
                at(flange, 3.5, 0.0, 0.25),
                at(flange, 1.75, 0.0, 2.0),
            ],
            vec![
                at(flange, 2.5, 0.0, 0.25),
                at(flange, 0.0, 2.5, 2.0),
                at(flange, 3.0, 0.0, -0.25),
            ],
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
}

/// A solid written by the kernel's writer, read back and imported.
fn reimported(s: &Solid, op: u64) -> Result<Solid, Error> {
    use rusty_occt::occt_brep::{import, read, write};
    let text = write(s.topology(), s.resolution().linear()).unwrap();
    let doc = read(&text).unwrap();
    let [solid] = <[_; 1]>::try_from(import(&doc).solids).unwrap();
    Solid::imported_with(OperationId(op), solid.result.unwrap(), solid.tolerance).map(|(s, _)| s)
}

fn one(r: rusty_occt::Result<(Vec<Solid>, History)>) -> Solid {
    let (out, _) = r.unwrap();
    assert_eq!(out.len(), 1);
    out.into_iter().next().unwrap()
}

/// The kernel's own body of the dimple's form, written, read back and
/// imported (its chain only in the inside-first order): a rod of radius 3
/// over `[0, 6]` less a ball of radius 5 about its axis's point at height
/// `-9/2` (the cup's form: the kernel's writer writes no meeting of a ball
/// and a rod's wall, `OCC14536-sh1`'s rim, nor the flange's spigot's ring on
/// its recess's ceiling, a circle the kernel stores as an ellipse); its own
/// volume, its stored vertices on its boundary, its Booleans with a box
/// turned off every axis the kernel's own results'.
#[test]
fn a_kernel_body_imported_inside_first() {
    let beside = frame([1.0, -2.0, 0.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let f = Frame3::new(
        beside.point(Point2::new(2.5, 0.5), 1.5),
        Vec3::new(2.0, 3.0, 6.0),
        Vec3::new(3.0, -6.0, 2.0),
        tolerance(),
    )
    .unwrap();
    let square = Boundary::polygon(
        vec![
            Point2::new(-0.5, -0.5),
            Point2::new(0.5, -0.5),
            Point2::new(0.5, 0.5),
            Point2::new(-0.5, 0.5),
        ],
        tolerance(),
    )
    .unwrap();
    let probe = Solid::extrude_with(
        OperationId(16),
        Profile::new(square, vec![], tolerance()).unwrap(),
        f,
        -1.5,
        1.5,
    )
    .unwrap()
    .0;
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    let near = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(1.0);
    let s = {
        let profile = Profile::new(
            Boundary::circle(Point2::new(0.0, 0.0), 3.0, tolerance()).unwrap(),
            vec![],
            tolerance(),
        )
        .unwrap();
        let rod = Solid::extrude_with(OperationId(52), profile, beside, 0.0, 6.0)
            .unwrap()
            .0;
        let c = beside.point(Point2::new(0.0, 0.0), -4.5).to_array();
        let centre = frame(c, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let half = std::f64::consts::FRAC_PI_2;
        let ball = Solid::sphere_with(OperationId(53), centre, 5.0, -half, half, tolerance())
            .unwrap()
            .0;
        one(rod.cut(OperationId(54), &ball))
    };
    let i = reimported(&s, 50).unwrap();
    let (v0, v1) = (s.mass_properties().volume, i.mass_properties().volume);
    assert!(near(v0, v1), "{v1} for {v0}");
    for v in i.topology().vertices() {
        assert_eq!(i.classify(v.position).unwrap(), Location::Boundary);
    }
    for (direct, again) in [
        (
            s.fuse(OperationId(70), &probe),
            i.fuse(OperationId(71), &probe),
        ),
        (
            s.cut(OperationId(70), &probe),
            i.cut(OperationId(71), &probe),
        ),
        (
            s.common(OperationId(70), &probe),
            i.common(OperationId(71), &probe),
        ),
    ] {
        let (d, a) = (volume(direct), volume(again));
        assert!(near(d, a), "{a} for {d}");
    }
}
