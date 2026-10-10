//! S9e.4b.4c.2a: imported bodies of one cylinder or cone face and plane faces
//! that are a Boolean tree of their primitive and convex hulls of their
//! planes whose pockets hold pockets of their own (S9e.4b.3c.3b's tooth, a
//! post standing in a pocket, a U island whose notch is a pocket again, on
//! cylinders),
//! given to Booleans, against the independent reference
//! (`fixtures/boolean-deep-trees-*` from
//! `tools/generate_deep_trees_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); the kernel's own such bodies (the `boolean` fuzz
//! target's `NESTED` ones among them) written, read back and imported, and a
//! tree four deep refused. Each case runs once (on a few threads) for the
//! checks that read its result.
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
    protocol::cases(include_str!("../../fixtures/boolean-deep-trees-cases.txt"))
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
    for line in include_str!("../../fixtures/boolean-deep-trees-expected.tsv")
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
        "tooth_pin_cut",
        "post_rod_common",
        "well_rod_fuse",
        "rod_post_cut",
        "post_cyl_cut",
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

/// Every body imports as its tree (its model matched on import), its volume
/// its closed form where it has one (the tooth's the reference's), its stored
/// vertices on its boundary, points in its material inside and points in its
/// pockets (each depth's) outside.
#[test]
fn imported_bodies_are_their_trees() {
    use std::f64::consts::PI;
    let post = frame([4.875, 3.125, -1.25], [0.0, 0.0, 1.0], [8.0, 15.0, 0.0]);
    let well = frame([1.0, -1.0, 0.0], [0.0, 0.0, 1.0], [12.0, 5.0, 0.0]);
    let at = |f: Frame3, u: f64, v: f64, w: f64| f.point(Point2::new(u, v), w).to_array();
    // The tooth's volume: its cases' cut and common with the pin.
    let expect = expected();
    let tooth =
        expect["tooth_pin_cut"].2.unwrap().1[0] + expect["tooth_pin_common"].2.unwrap().1[0];
    let bodies = [
        (
            "form_tooth",
            tooth,
            vec![[2.5, 3.75, 3.25], [3.5, 4.0, 3.25], [4.0, 4.0, 1.0]],
            vec![[1.5, 3.75, 3.25], [2.5, 2.75, 3.25], [2.5, 5.0, 4.0]],
        ),
        (
            "deep_post",
            96.0 * PI - 63.0,
            vec![
                at(post, 0.0, 0.0, 4.5),
                at(post, 3.0, 0.0, 4.5),
                at(post, 1.5, 0.0, 2.0),
            ],
            vec![at(post, 1.75, 0.0, 4.5), at(post, -0.5, 1.5, 3.5)],
        ),
        (
            "deep_well",
            104.0 * PI - 72.0,
            vec![
                at(well, 0.0, 1.0, 4.0),
                at(well, -1.0, 0.0, 5.0),
                at(well, 0.0, 0.0, 1.5),
                at(well, 3.5, 0.0, 5.0),
            ],
            vec![
                at(well, 0.5, 0.0, 4.0),
                at(well, 2.0, 0.0, 4.0),
                at(well, 0.0, -2.0, 3.0),
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

/// A prism of a profile with holes on a frame over `[h0, h1]`.
fn prism(f: Frame3, outline: Boundary, holes: Vec<Boundary>, h0: f64, h1: f64, op: u64) -> Solid {
    let profile = Profile::new(outline, holes, tolerance()).unwrap();
    Solid::extrude_with(OperationId(op), profile, f, h0, h1)
        .unwrap()
        .0
}

fn polygon(points: &[(f64, f64)]) -> Boundary {
    let pts: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    Boundary::polygon(pts, tolerance()).unwrap()
}

fn square(a: f64) -> Boundary {
    polygon(&[(-a, -a), (a, -a), (a, a), (-a, a)])
}

/// A U of half side `a` less its notch `[-b, a] x [-c, c]`, open toward `+x`.
fn u_shape(a: f64, b: f64, c: f64) -> Boundary {
    polygon(&[
        (-a, -a),
        (a, -a),
        (a, -c),
        (-b, -c),
        (-b, c),
        (a, c),
        (a, a),
        (-a, a),
    ])
}

/// A solid written by the kernel's writer, read back and imported.
fn reimported(s: &Solid, op: u64) -> Result<Solid, Error> {
    use rusty_occt::occt_brep::{import, read, write};
    let text = write(s.topology(), s.resolution().linear()).unwrap();
    let doc = read(&text).unwrap();
    let [solid] = <[_; 1]>::try_from(import(&doc).solids).unwrap();
    Solid::imported_with(OperationId(op), solid.result.unwrap(), solid.tolerance).map(|x| x.0)
}

fn one(r: rusty_occt::Result<(Vec<Solid>, History)>) -> Solid {
    let (out, _) = r.unwrap();
    assert_eq!(out.len(), 1);
    out.into_iter().next().unwrap()
}

/// The kernel's own bodies of nested pockets, written by its writer, read
/// back and imported: a post and a U island in a square pocket of a
/// cylinder (the `boolean` fuzz target's `NESTED` bodies, on the world's
/// axes and turned about its `z`) and of a frustum; each its tree, its
/// volume, its stored vertices on its boundary, and its Booleans with a
/// turned box the kernel's own results', all 15 evaluating. A U island whose notch holds a
/// tongue (pockets four deep) is refused with the limit's text.
#[test]
fn kernel_bodies_imported_as_trees() {
    let tol = tolerance();
    let world = frame([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    // Turned about the world's `z`: walls along a cylinder's axis turned off
    // it (the fuzz target's tilted frame) lie within rounding of its
    // direction under a correctly rounded `hypot` (S9's refusal of a plane
    // within rounding of a cylinder's direction), on macOS's exactly.
    let turned = frame([1.0, -2.0, 0.5], [0.0, 0.0, 1.0], [3.0, 4.0, 0.0]);
    let beside = frame([1.5, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    // The fuzz target's: a cylinder of radius `s` over `[0, h + 1]` less a
    // square of half side `5k` with a hole over `[(h + 1) / 2, h + 2]`, `k =
    // s / 8`.
    let nested = |f: Frame3, s: f64, h: f64, u: bool, op: u64| {
        let k = s / 8.0;
        let hole = if u {
            u_shape(3.0 * k, k, k)
        } else {
            square(2.0 * k)
        };
        let top = h + 1.0;
        let cylinder = Solid::cylinder_with(OperationId(op), f, s, 0.0, top, tol)
            .unwrap()
            .0;
        let ring = prism(f, square(5.0 * k), vec![hole], top / 2.0, top + 1.0, op + 1);
        one(cylinder.cut(OperationId(op + 2), &ring))
    };
    let frustum = Solid::cone_with(OperationId(30), world, 5.0, 4.0, 6.0, tol)
        .unwrap()
        .0;
    let well = prism(
        world,
        square(2.5),
        vec![u_shape(1.5, 0.5, 0.5)],
        2.5,
        7.0,
        31,
    );
    let bodies = [
        ("post", nested(world, 4.0, 1.75, false, 10)),
        ("island", nested(beside, 3.0, 1.25, true, 13)),
        ("turned post", nested(turned, 3.0, 1.5, false, 16)),
        ("turned island", nested(turned, 4.0, 0.75, true, 19)),
        ("well", one(frustum.cut(OperationId(32), &well))),
    ];
    let probe = prism(
        frame([0.25, -0.5, 1.25], [3.0, 0.0, 4.0], [0.0, 1.0, 0.0]),
        square(1.25),
        vec![],
        0.0,
        2.5,
        40,
    );
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    let mut evaluated = 0;
    for (name, s) in &bodies {
        let i = reimported(s, 50).unwrap_or_else(|e| panic!("{name}: {e}"));
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
                    evaluated += 1;
                }
                (Err(d), Err(a)) => assert_eq!(d, a, "{name}"),
                (d, a) => panic!("{name}: {:?} and {:?}", d.map(|_| ()), a.map(|_| ())),
            }
        }
    }
    assert_eq!(evaluated, 15);
    // Four deep: the U island's notch holds a tongue from its back wall.
    let tongue = polygon(&[
        (-1.5, -1.5),
        (1.5, -1.5),
        (1.5, -0.75),
        (-0.75, -0.75),
        (-0.75, -0.25),
        (0.75, -0.25),
        (0.75, 0.25),
        (-0.75, 0.25),
        (-0.75, 0.75),
        (1.5, 0.75),
        (1.5, 1.5),
        (-1.5, 1.5),
    ]);
    let cylinder = Solid::cylinder_with(OperationId(70), world, 4.0, 0.0, 6.0, tol)
        .unwrap()
        .0;
    let deep = one(cylinder.cut(
        OperationId(72),
        &prism(world, square(2.5), vec![tongue], 3.0, 7.0, 71),
    ));
    match reimported(&deep, 50) {
        Err(Error::OutOfDomain(m)) => assert_eq!(
            m,
            "an imported plane piece other than a Boolean tree of its primitive and its planes' \
             hulls, its pockets nested at most three deep (S9e.4b.4)"
        ),
        other => panic!("four deep: {:?}", other.map(|_| ())),
    }
}
