//! S9e.4b.3c.3b: imported bodies of one sphere or cylinder face and plane
//! faces that are a Boolean tree of their primitive and several convex hulls
//! of their planes (a U prism's groove and boss, a cylinder bitten twice by
//! a U, a hemisphere's disc fused with or cut from a box, a flattened
//! cylindrical boss) given to Booleans, against the independent reference
//! (`fixtures/boolean-piece-trees-*` from
//! `tools/generate_piece_trees_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); the `boolean` fuzz target's two bodies, a
//! cylinder's oblique end and a primitive bitten twice built by the kernel
//! and imported. Each case runs once (on a few threads) for the checks that
//! read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, History, Resolution};
use rusty_occt::identity::{EntityId, OperationId};
use rusty_occt::{
    Boundary, Error, Frame3, Location, Point2, Point3, Profile, RigidTransform, Solid, Tolerance,
    Vec3,
};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!("../../fixtures/boolean-piece-trees-cases.txt"))
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
    for line in include_str!("../../fixtures/boolean-piece-trees-expected.tsv")
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

/// What an input entity became in the result.
fn targets(h: &History, id: EntityId) -> Vec<EntityId> {
    match h.resolve(id) {
        Resolution::Same(x) => vec![x],
        Resolution::Split(v) => v,
        Resolution::Merged { into, .. } => vec![into],
        Resolution::Deleted | Resolution::Unknown => Vec::new(),
    }
}

/// The U's two arm ends on one plane (one plane of its hull, their stored
/// frames rounded apart) each continue into the rod's fuse with the groove,
/// as themselves or split, never merged into one another.
#[test]
fn faces_on_one_plane_keep_their_ids() {
    let all: BTreeMap<&str, _> = runs().iter().map(|(n, r)| (n.as_str(), r)).collect();
    let Ok((a, _, out, h)) = all["u_scoop_rod_fuse"] else {
        panic!("the groove and the rod's fuse")
    };
    let t = a.topology();
    // The stored faces normal to the U's `x` at its open end: the arms'.
    let x = Vec3::new(-7.0, -4.0, 4.0) * (1.0 / 9.0);
    let end = Point3::new(1.0, 2.0, 0.5) + x * 12.0;
    let mut arms = Vec::new();
    for (k, f) in t.faces().iter().enumerate() {
        if let rusty_occt::topology::Surface::Plane(frame) = &f.surface {
            let along = frame.normal().cross(x).length() < 1e-9;
            if along && (frame.origin() - end).dot(x).abs() < 1e-9 {
                arms.push(
                    t.id_of(rusty_occt::topology::Slot::Face(
                        rusty_occt::topology::FaceId::new(k),
                    ))
                    .unwrap(),
                );
            }
        }
    }
    assert_eq!(arms.len(), 2, "the U's arm ends");
    let outs: Vec<Vec<EntityId>> = arms.iter().map(|&id| targets(h, id)).collect();
    assert!(outs.iter().all(|o| !o.is_empty()), "{outs:?}");
    assert!(outs[0].iter().all(|x| !outs[1].contains(x)), "{outs:?}");
    assert!(!out.is_empty());
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
        "u_scoop_rod_cut",
        "u_boss_slab_common",
        "bites_box_cut",
        "cap_boss_rod_fuse",
        "rod_cap_pocket_cut",
        "dee_boss_ball_common",
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

/// Every body imports as its tree (its model matched on import), its
/// volume its closed form where it has one, its stored vertices on its
/// boundary, a point in its material inside and one in a part its tree
/// removes outside; the U's notch is refused for its own tangency and the
/// tooth (a pocket within a pocket) as S9e.4b.4's.
#[test]
fn imported_bodies_are_their_trees() {
    use std::f64::consts::PI;
    let segment = |r: f64, d: f64| r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let skew4 = frame([1.0, 2.0, 0.5], [4.0, 1.0, 8.0], [-7.0, -4.0, 4.0]);
    let skew = frame([1.0, 2.0, 0.5], [8.0, 4.0, 1.0], [-1.0, 4.0, -8.0]);
    let skew2 = frame([0.5, 1.0, 1.5], [4.0, 4.0, 7.0], [1.0, -8.0, 4.0]);
    let at = |f: Frame3, u: f64, v: f64, w: f64| f.point(Point2::new(u, v), w).to_array();
    let bodies = [
        (
            "form_u_scoop",
            None,
            at(skew4, 2.0, 2.0, 3.0),
            at(skew4, 8.0, 1.5, 5.5),
        ),
        (
            "form_u_boss",
            Some(304.0 + 6.75 * PI),
            [2.5, 4.0, 6.0],
            [8.0, 4.0, 2.0],
        ),
        ("form_bites", None, [4.0, 4.0, 3.0], [2.5, 6.0, 3.0]),
        (
            "form_cap_boss",
            Some(384.0 + 26.0 * PI / 3.0),
            at(skew, 4.0, 4.0, -0.5),
            at(skew, 0.5, 0.5, -0.5),
        ),
        (
            "form_cap_pocket",
            Some(256.0 - 46.0 * PI / 3.0),
            at(skew2, 1.0, 1.0, 1.0),
            at(skew2, 4.0, 4.0, 3.0),
        ),
        (
            "form_dee_boss",
            Some(300.0 + 6.0 * (3.125 * PI - segment(2.5, 1.5) / 2.0)),
            [8.0, 1.0, 4.5],
            [8.0, 3.0, 4.5],
        ),
    ];
    for (name, volume, inside, outside) in bodies {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let m = s.mass_properties();
        if let Some(volume) = volume {
            assert!(
                (m.volume - volume).abs() <= 1e-9 * volume,
                "{name}: {} for {volume}",
                m.volume
            );
        }
        for (p, want) in [(inside, Location::Inside), (outside, Location::Outside)] {
            assert_eq!(
                s.classify(Point3::new(p[0], p[1], p[2])).unwrap(),
                want,
                "{name} {p:?}"
            );
        }
        for v in s.topology().vertices() {
            assert_eq!(
                s.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
    }
    match body("form_u_notch", 91) {
        Err(Error::Degenerate(m)) => assert_eq!(
            m,
            "an imported plane piece whose curved face is tangent to its plane faces"
        ),
        other => panic!("form_u_notch: {:?}", other.map(|_| ())),
    }
    match body("form_tooth", 91) {
        Err(Error::OutOfDomain(m)) => assert_eq!(
            m,
            "an imported plane piece other than a Boolean tree of its primitive and its planes' \
             hulls (S9e.4b.4)"
        ),
        other => panic!("form_tooth: {:?}", other.map(|_| ())),
    }
}

/// A prism of a profile on a frame over `[h0, h1]`.
fn prism(f: Frame3, outline: Boundary, h0: f64, h1: f64, op: u64) -> Solid {
    let profile = Profile::new(outline, vec![], tolerance()).unwrap();
    Solid::extrude_with(OperationId(op), profile, f, h0, h1)
        .unwrap()
        .0
}

fn polygon(points: &[(f64, f64)]) -> Boundary {
    let pts: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    Boundary::polygon(pts, tolerance()).unwrap()
}

fn rectangle(x0: f64, y0: f64, x1: f64, y1: f64) -> Boundary {
    polygon(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
}

/// A U: the rectangle `[x0, x1] x [y0, y1]` less the slot `[s0, x1] x [t0,
/// t1]`.
fn u_shape(x0: f64, y0: f64, x1: f64, y1: f64, s0: f64, t0: f64, t1: f64) -> Boundary {
    polygon(&[
        (x0, y0),
        (x1, y0),
        (x1, t0),
        (s0, t0),
        (s0, t1),
        (x1, t1),
        (x1, y1),
        (x0, y1),
    ])
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

/// The kernel's own bodies of trees, imported: the `boolean` fuzz target's
/// two (a box fused with a ball's half whose disc plane cuts it; a U prism
/// less a ball about its end's corner), a rod bitten twice by a U, written
/// by its writer and read back; a cylinder's piece of an oblique plane
/// fused with a box (its end an ellipse, which the reader does not take),
/// its topology imported directly. Each its tree: its volume, and its
/// Booleans with a turned box, the kernel's own results'.
#[test]
fn kernel_bodies_imported_as_trees() {
    let tol = tolerance();
    let half = std::f64::consts::FRAC_PI_2;
    let world = frame([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let one = |r: rusty_occt::Result<(Vec<Solid>, History)>| {
        let (out, _) = r.unwrap();
        assert_eq!(out.len(), 1);
        out.into_iter().next().unwrap()
    };
    // The fuzz target's ball's half: the ball's frame turned about `z`.
    let ball_half = Solid::sphere_with(
        OperationId(2),
        frame([0.375, -0.125, 1.125], [0.0, 0.0, 1.0], [0.6, 0.8, 0.0]),
        2.8125,
        -half,
        0.0,
        tol,
    )
    .unwrap()
    .0;
    let block = prism(world, rectangle(0.0, 0.0, 1.25, 2.25), 0.0, 2.25, 1);
    // The fuzz target's U: on a turned frame, a ball about its end's corner
    // (its frame's origin) through its whole thickness.
    let g = frame([1.0, -2.0, 0.5], [0.0, 0.6, 0.8], [1.0, 0.0, 0.0]);
    let u = prism(
        g,
        u_shape(0.0, 0.0, 11.25, 3.0, 3.75, 1.0, 2.0),
        0.0,
        0.75,
        3,
    );
    let corner = Solid::sphere_with(OperationId(4), g, 1.6875, -half, half, tol)
        .unwrap()
        .0;
    // A rod bitten twice by a U along its axis between its caps.
    let rod = Solid::cylinder_with(
        OperationId(5),
        frame([4.0, 4.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        3.0,
        0.0,
        6.0,
        tol,
    )
    .unwrap()
    .0;
    let fork = prism(
        world,
        u_shape(-3.0, 1.5, 3.5, 6.5, -1.0, 2.5, 5.5),
        2.0,
        4.5,
        6,
    );
    // A cylinder's piece of an oblique plane fused with a box holding its
    // lower end: its end an ellipse.
    let cylinder = Solid::cylinder_with(
        OperationId(10),
        frame([5.0, 5.0, -1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        2.0,
        0.0,
        8.0,
        tol,
    )
    .unwrap()
    .0;
    let (pieces, _) = cylinder
        .split_by_plane(
            OperationId(11),
            frame([5.0, 5.0, 4.0], [0.6, 0.0, 0.8], [0.0, 1.0, 0.0]),
        )
        .unwrap();
    let low = pieces
        .into_iter()
        .map(|(_, s)| s)
        .min_by(|a, b| a.bounds().min.z.total_cmp(&b.bounds().min.z))
        .unwrap();
    let base = prism(world, rectangle(1.0, 1.0, 9.0, 9.0), -2.0, 1.0, 12);
    let oblique = one(base.fuse(OperationId(13), &low));
    let written = [
        ("ball's half", one(block.fuse(OperationId(7), &ball_half))),
        ("U's corner", one(u.cut(OperationId(8), &corner))),
        ("rod bitten twice", one(rod.cut(OperationId(9), &fork))),
    ];
    let mut bodies: Vec<(&str, Solid, Solid)> = written
        .into_iter()
        .map(|(name, s)| {
            let i = reimported(&s, 50);
            (name, s, i)
        })
        .collect();
    let i = Solid::imported_with(
        OperationId(50),
        oblique.topology().clone(),
        oblique.resolution(),
    )
    .unwrap()
    .0;
    bodies.push(("oblique end", oblique, i));
    let probe = prism(
        frame([2.5, 2.0, 0.25], [0.0, 3.0, 4.0], [3.0, 0.0, 4.0]),
        rectangle(0.0, 0.0, 3.0, 3.0),
        0.0,
        2.5,
        14,
    );
    let volume = |r: rusty_occt::Result<(Vec<Solid>, History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    for (name, s, i) in &bodies {
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
