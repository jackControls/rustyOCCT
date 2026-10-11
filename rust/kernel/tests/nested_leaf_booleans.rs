//! S9e.4b.4c.2b.1: imported bodies of several cylinder faces and plane
//! faces that are a prism leaf with a pocket holding pockets of its own (a
//! stadium plate less a square with a square hole, a post standing in the
//! pocket; less a square with a U-shaped hole, a U island whose notch is a
//! pocket of its own) given to Booleans, against the independent reference
//! (`fixtures/boolean-nested-leaves-*` from
//! `tools/generate_nested_leaves_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); the `boolean` fuzz target's `NESTED_LEAVES` bodies
//! and the kernel's own two leaves with such a pocket built, written by its
//! writer, read back and imported. Each case runs once (on a few threads) for
//! the checks that read its result.
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
        "../../fixtures/boolean-nested-leaves-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-nested-leaves-expected.tsv")
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
        "island_rod_common",
        "island_wall_fuse",
        "rod_post_cut",
        "island_cyl_fuse",
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

/// The refusal of a body of several primitives with other plane faces than
/// this step and the earlier ones take.
const NOT_THEIR_FACES: &str = "an imported body of several primitives with plane faces other than \
                               their ends, one or two prisms', a flat or a pocket, its pockets \
                               nested at most three deep (S9e.4b.4b.2b.2b.3b)";

/// Every body imports as its plate with its nested pockets (its model matched
/// on import), its volume its closed form, its stored vertices on its
/// boundary, points in its material inside (the post, the island's arms) and
/// points in each depth's pockets outside (the pocket round the post or the
/// island, the island's notch); S9e.4b.4b.2a's notch (both caps cut at its
/// rim) is refused with the new text.
#[test]
fn imported_bodies_are_their_nested_pockets() {
    use std::f64::consts::PI;
    let post = frame([1.0, 0.5, 0.25], [0.0, 0.0, 1.0], [8.0, 15.0, 0.0]);
    let island = frame([3.3, 3.885, 0.625], [0.0, 0.0, 1.0], [12.0, 5.0, 0.0]);
    let at = |f: Frame3, u: f64, v: f64, w: f64| f.point(Point2::new(u, v), w).to_array();
    type Points = Vec<[f64; 3]>;
    let bodies: [(&str, f64, Points, Points); 2] = [
        (
            "nested_post",
            90.0 + 27.0 * PI,
            vec![
                at(post, 3.0, 0.0, 2.5),
                at(post, 0.0, 0.0, 1.0),
                at(post, 7.0, 0.0, 2.5),
            ],
            vec![
                at(post, 1.5, 0.0, 2.5),
                at(post, 3.0, 1.5, 2.0),
                at(post, 3.0, 0.0, 3.5),
            ],
        ),
        (
            "nested_island",
            72.0 + 27.0 * PI,
            vec![
                at(island, 2.0, 1.0, 2.0),
                at(island, 4.0, 1.0, 2.5),
                at(island, 3.5, 0.0, 0.5),
            ],
            vec![
                at(island, 3.5, 0.0, 2.0),
                at(island, 5.0, 2.0, 2.0),
                at(island, 0.8, 0.0, 2.0),
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
    match body("leaf_notch", 91) {
        Err(Error::OutOfDomain(m)) => assert_eq!(m, NOT_THEIR_FACES),
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

fn polygon(points: &[(f64, f64)]) -> Boundary {
    let points = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    Boundary::polygon(points, tolerance()).unwrap()
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Boundary {
    polygon(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
}

fn stadium(cx: f64, cy: f64, length: f64, r: f64) -> Boundary {
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
        tolerance(),
    )
    .unwrap()
}

fn prism(f: Frame3, outer: Boundary, holes: Vec<Boundary>, h0: f64, h1: f64, op: u64) -> Solid {
    Solid::extrude_with(
        OperationId(op),
        Profile::new(outer, holes, tolerance()).unwrap(),
        f,
        h0,
        h1,
    )
    .unwrap()
    .0
}

fn one(r: rusty_occt::Result<(Vec<Solid>, History)>) -> Solid {
    let (out, _) = r.unwrap();
    assert_eq!(out.len(), 1);
    out.into_iter().next().unwrap()
}

/// The frame turned about `f`'s axis, its `x` along `8 x + 15 y`.
fn turned(f: Frame3) -> Frame3 {
    Frame3::new(
        f.origin(),
        f.normal(),
        f.x() * 8.0 + f.y() * 15.0,
        tolerance(),
    )
    .unwrap()
}

/// The U-shaped hole `[3/2, 9/2] x [-3/2, 3/2]` less its notch `[5/2, 9/2] x
/// [-1/2, 1/2]`, scaled by `k`.
fn u_hole(k: f64) -> Boundary {
    let p = [
        (1.5, -1.5),
        (4.5, -1.5),
        (4.5, -0.5),
        (2.5, -0.5),
        (2.5, 0.5),
        (4.5, 0.5),
        (4.5, 1.5),
        (1.5, 1.5),
    ];
    polygon(&p.map(|(x, y)| (x * k, y * k)))
}

/// The `boolean` fuzz target's `NESTED_LEAVES` bodies (`nested_leaves_body`,
/// `k = s / 8`, on the frame turned about `f`'s axis): the plate of half
/// discs of radius `3k` about `(0, 0)` and `(6k, 0)` over `[0, h]` less, over
/// `[h / 2, h + 1]`, the square `[k, 5k] x [-2k, 2k]` with the hole `[2k, 4k]
/// x [-k, k]` (a post in a pocket) or, where `u`, the square `[k / 2, 11k /
/// 2] x [-5k / 2, 5k / 2]` with the U-shaped hole (a U island, its notch a
/// pocket of its own).
fn nested_leaves_body(f: Frame3, s: f64, h: f64, u: bool) -> Solid {
    let k = s / 8.0;
    let f = turned(f);
    let plate = prism(f, stadium(0.0, 0.0, 6.0 * k, 3.0 * k), vec![], 0.0, h, 52);
    let (outer, hole) = if u {
        (square(0.5 * k, -2.5 * k, 5.5 * k, 2.5 * k), u_hole(k))
    } else {
        (
            square(k, -2.0 * k, 5.0 * k, 2.0 * k),
            square(2.0 * k, -k, 4.0 * k, k),
        )
    };
    let ring = prism(f, outer, vec![hole], h / 2.0, h + 1.0, 53);
    one(plate.cut(OperationId(54), &ring))
}

/// The kernel's own bodies of a pocket with pockets of its own, written,
/// read back and imported, each its plate (or its two leaves) less its
/// nested pockets: its volume, its stored vertices on its boundary, and its
/// Booleans with a box turned off every axis the kernel's own results'. The
/// `boolean` fuzz target's `NESTED_LEAVES` bodies on the world's axes and
/// on a frame turned about the world's `z` (on the target's tilted frame
/// refused as before the step), and two leaves (a box with a stadium through it) less a square with a square
/// hole: a pocket within a pocket of a body of two leaves, which no one
/// Boolean of OCCT's `write` blocks holds. A plate whose U island's notch
/// holds a tongue (pockets four deep) is refused as before the step.
#[test]
fn kernel_bodies_imported_with_their_nested_pockets() {
    let world = Frame3::xy();
    let beside = frame([1.5, 1.0, 0.0], [0.0, 0.0, 1.0], [3.0, 4.0, 0.0]);
    let tilted = frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let probe_at = |at: Point3| {
        // Turned off every body's axis.
        let f = Frame3::new(
            at,
            Vec3::new(2.0, 3.0, 6.0),
            Vec3::new(3.0, -6.0, 2.0),
            tolerance(),
        )
        .unwrap();
        prism(f, square(-0.5, -0.5, 0.5, 0.5), vec![], -1.5, 1.5, 16)
    };
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    let near = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(1.0);
    // Two leaves less a square with a square hole within the box: the box
    // `[-5, 5]^2` over `[0, 2]` with the stadium of half discs of radius 2
    // about `(-1, 1)` and `(1, 1)` over `[-1, 3]` through it, less the square
    // `[-4, 4] x [-9/2, -2]` with the hole `[-2, 2] x [-7/2, -3]` over `[1,
    // 3]`, on the frame turned about the world's `z`.
    let two = {
        let f = turned(world);
        let boxed = prism(f, square(-5.0, -5.0, 5.0, 5.0), vec![], 0.0, 2.0, 60);
        let through = prism(f, stadium(-1.0, 1.0, 2.0, 2.0), vec![], -1.0, 3.0, 61);
        let leaves = one(boxed.fuse(OperationId(62), &through));
        let ring = prism(
            f,
            square(-4.0, -4.5, 4.0, -2.0),
            vec![square(-2.0, -3.5, 2.0, -3.0)],
            1.0,
            3.0,
            63,
        );
        one(leaves.cut(OperationId(64), &ring))
    };
    let bodies = [
        (
            "post",
            nested_leaves_body(world, 8.0, 1.5, false),
            turned(world).point(Point2::new(4.0, 1.0), 1.0),
        ),
        (
            "island",
            nested_leaves_body(world, 8.0, 1.5, true),
            turned(world).point(Point2::new(3.5, 0.5), 1.25),
        ),
        (
            "turned post",
            nested_leaves_body(beside, 6.0, 1.25, false),
            turned(beside).point(Point2::new(2.25, 0.75), 1.0),
        ),
        (
            "turned island",
            nested_leaves_body(beside, 6.0, 1.25, true),
            turned(beside).point(Point2::new(2.5, 0.5), 1.0),
        ),
        (
            "two leaves",
            two,
            turned(world).point(Point2::new(0.5, -3.25), 1.75),
        ),
    ];
    let mut evaluated = 0;
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
            match (direct, again) {
                (Ok(d), Ok(a)) => {
                    let (d, a) = (volume(Ok(d)), volume(Ok(a)));
                    assert!(near(d, a), "{name}: {a} for {d}");
                    evaluated += 1;
                }
                (Err(d), Err(a)) => assert_eq!(d, a, "{name}"),
                (d, a) => panic!("{name}: {:?} and {:?}", d.map(|_| ()), a.map(|_| ())),
            }
        }
    }
    assert_eq!(evaluated, 15);
    // On the target's tilted frame the pocket's walls, stored along the
    // plate's rounded axis, lie within rounding of its cylinders' direction
    // (S9's `Degenerate("a plane within rounding of a cylinder's
    // direction")` in the nested chain's Booleans): the body's refusal as
    // before the step, its leaf's tangent joints.
    for u in [false, true] {
        match reimported(&nested_leaves_body(tilted, 6.0, 1.25, u), 50) {
            Err(e) => assert_eq!(
                e,
                Error::Degenerate(
                    "an imported body of several primitives whose faces are tangent along an edge"
                ),
                "tilted {u}"
            ),
            other => panic!("tilted {u}: {:?}", other.map(|_| ())),
        }
    }
    // Four deep: the U island's notch holds a tongue from its back wall.
    let tongue = polygon(&[
        (1.5, -1.5),
        (4.5, -1.5),
        (4.5, -0.75),
        (2.25, -0.75),
        (2.25, -0.25),
        (3.75, -0.25),
        (3.75, 0.25),
        (2.25, 0.25),
        (2.25, 0.75),
        (4.5, 0.75),
        (4.5, 1.5),
        (1.5, 1.5),
    ]);
    let f = turned(world);
    let plate = prism(f, stadium(0.0, 0.0, 6.0, 3.0), vec![], 0.0, 3.0, 80);
    let deep = one(plate.cut(
        OperationId(82),
        &prism(f, square(0.5, -2.5, 5.5, 2.5), vec![tongue], 1.0, 4.0, 81),
    ));
    match reimported(&deep, 50) {
        Err(e) => assert_eq!(
            e,
            Error::Degenerate(
                "an imported body of several primitives whose faces are tangent along an edge"
            ),
            "four deep"
        ),
        other => panic!("four deep: {:?}", other.map(|_| ())),
    }
}
