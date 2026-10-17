//! S9e.4b.4b.2b.2b.3b: imported bodies of several cylinder faces and plane
//! faces that are a prism leaf whose caps a leaf crossing it splits (a
//! stadium plate with a stadium slab through it on crossing axes) given to
//! Booleans, against the independent reference
//! (`fixtures/boolean-crossed-leaves-*` from
//! `tools/generate_crossed_leaves_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); each group's six Booleans consistent with one
//! another; the `boolean` fuzz target's measured body (the kernel's own
//! fuse of a plate and a slab on the world's axes) written by its writer,
//! read back and imported; a plate two slabs split refused as before. Each
//! case runs once (on a few threads) for the checks that read its result.
#[path = "support/turned_consistency.rs"]
#[allow(dead_code)]
mod consistency;
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
        "../../fixtures/boolean-crossed-leaves-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-crossed-leaves-expected.tsv")
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
        "joint_rod_cut",
        "joint_arc_common",
        "waist_rod_fuse",
        "rod_joint_cut",
        "joint_cyl_fuse",
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
const GROUPS: [&str; 5] = [
    "joint_rod",
    "joint_arc",
    "waist_rod",
    "rod_joint",
    "joint_cyl",
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

/// Every body imports as its plate fused with its slab (its model matched on
/// import): its volume the closed form (the plate's and the slab's sections
/// less their common part), points inside the hidden joints (the plate's
/// joint inside the slab), in each piece and in the slab above the plate
/// inside, points above a piece and below the slab's lower lobe beside the
/// plate outside, its stored vertices on its boundary.
#[test]
fn imported_bodies_are_their_plates_and_slabs() {
    use std::f64::consts::PI;
    let joint = frame([5.25, 1.625, 0.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let waist = frame([-2.0, 1.0, 0.25], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let at = |f: Frame3, u: f64, v: f64, w: f64| f.point(Point2::new(u, v), w).to_array();
    // The joint's common part: the plate's section over `[-1/4, 3/4]`, its
    // half disc's `x sqrt(4 - x^2) + 4 asin(x / 2)` from `-1/4` to 0 and the
    // strip `3/4 x 4`, over `3/2`.
    let chord = 0.25 * (4.0f64 - 0.0625).sqrt() + 4.0 * (0.125f64).asin();
    let slab = 16.0 + 4.0 * PI;
    type Points = Vec<[f64; 3]>;
    let bodies: [(&str, f64, Points, Points); 2] = [
        (
            "crossed_joint",
            1.5 * (16.0 + 4.0 * PI) - 1.5 * (chord + 3.0) + slab,
            vec![
                at(joint, 0.0, -1.9, 0.75),
                at(joint, 0.0, 1.9, 1.4),
                at(joint, 3.0, 0.0, 0.75),
                at(joint, -1.5, 0.0, 0.75),
                at(joint, 0.25, 1.0, 2.5),
            ],
            vec![
                at(joint, 3.0, 0.0, 2.0),
                at(joint, -1.5, 0.0, -0.25),
                at(joint, 0.25, 0.0, 5.0),
            ],
        ),
        (
            "crossed_waist",
            1.5 * (32.0 + 4.0 * PI) - 6.0 + slab,
            vec![
                at(waist, 3.5, -1.99, 0.75),
                at(waist, 1.0, 0.0, 0.75),
                at(waist, 7.0, 1.0, 0.75),
                at(waist, 3.5, 4.0, 0.75),
            ],
            vec![
                at(waist, 1.0, 0.0, 2.0),
                at(waist, 7.0, 3.0, 0.75),
                at(waist, 3.5, 0.0, 4.0),
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

fn prism(f: Frame3, outer: Boundary, h0: f64, h1: f64, op: u64) -> Solid {
    Solid::extrude_with(
        OperationId(op),
        Profile::new(outer, vec![], tolerance()).unwrap(),
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

/// The `boolean` fuzz target's measured body (the fuse of corpus input
/// `8400c1e8`): the stadium plate of half discs of radius 2 about `(0, 0)`
/// and `(4, 0)` over `[0, 3/2]` on the world's axes fused with the stadium
/// slab normal to the world's `x` over `[-1/4, 3/4]`, half discs of radius 2
/// about `(y, z) = (-1, 3/4)` and `(3, 3/4)`; at `x0` along the plate the
/// slab crosses (`0`: its line-arc joints, `2`: its straight middle).
fn plate_and_slab(x0: f64) -> Solid {
    let plate = prism(Frame3::xy(), stadium(0.0, 0.0, 4.0, 2.0), 0.0, 1.5, 52);
    let side = frame([x0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let slab = prism(side, stadium(-1.0, 0.75, 4.0, 2.0), -0.25, 0.75, 53);
    one(plate.fuse(OperationId(54), &slab))
}

/// The kernel's own bodies written, read back and imported: the fuzz
/// target's measured body and the slab across the plate's middle, each its
/// plate and slab (its volume the kernel's own, its stored vertices on its
/// boundary), their Booleans with a box turned off every axis the kernel's
/// own results'; the plate split by two slabs (three pieces, two leads)
/// refused as before the step.
#[test]
fn kernel_bodies_imported_as_crossed_leaves() {
    let probe_at = |at: Point3| {
        let f = Frame3::new(
            at,
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
        prism(f, square, -1.5, 1.5, 16)
    };
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    let near = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(1.0);
    let mut evaluated = 0;
    for (name, x0, at) in [
        ("joint", 0.0, Point3::new(1.3, 1.1, 1.3)),
        ("middle", 2.0, Point3::new(3.4, 0.4, 1.2)),
    ] {
        let s = plate_and_slab(x0);
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
    assert_eq!(evaluated, 6);
    // Two slabs across the plate: its caps three pieces each, the middle
    // one's loop running along both slabs; neither lead's crossed leaf
    // closes, so no crossed chain is its model and the body's refusal is the
    // one it had before the step.
    let plate = prism(Frame3::xy(), stadium(0.0, 0.0, 4.0, 2.0), 0.0, 1.5, 52);
    let mut two = plate;
    // (The second slab's half discs of radius 5/2 about `(-5/4, 3/4)` and
    // `(13/4, 3/4)`, its walls on no surface of the first's.)
    for (k, x0, sections) in [
        (0, 0.0, stadium(-1.0, 0.75, 4.0, 2.0)),
        (1, 3.0, stadium(-1.25, 0.75, 4.5, 2.5)),
    ] {
        let side = frame([x0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let slab = prism(side, sections, -0.25, 0.75, 53 + k);
        two = one(two.fuse(OperationId(60 + k), &slab));
    }
    match reimported(&two, 50) {
        Err(Error::Degenerate(m)) => assert_eq!(
            m,
            "an imported body of several primitives whose faces are tangent along an edge"
        ),
        other => panic!("two slabs: {:?}", other.map(|_| ())),
    }
}
