//! S9e.4b.3a: imported plane pieces of a sphere (bodies OCCT wrote to
//! `.brep` files: a ball's wedges, lunes and halves between planes through
//! its axis and normal to it), given to Booleans, against the independent
//! reference (`fixtures/boolean-imported-pieces-*` from
//! `tools/generate_imported_pieces_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`): each case's Boolean (or chain) with its imported
//! inputs made by `Solid::imported_with`, decided as the given model of the
//! piece's primitive common the half-spaces of its planes. Cylinders' and
//! cones' pieces, which no `.brep` file the reader takes can hold (OCCT's
//! ellipse records), are the kernel's own split pieces' topologies imported.
//! Each case runs once (on a few threads) for the checks that read its
//! result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Location, Point2, Point3, Profile, RigidTransform, Solid, Tolerance,
    Vec3,
};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-imported-pieces-cases.txt"
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
type Expected = std::collections::BTreeMap<String, (String, Option<(usize, [f64; 2])>)>;

fn expected() -> Expected {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-imported-pieces-expected.tsv")
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
            ("unsupported", Err(Error::OutOfDomain(_))) => continue,
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

/// The declared refusals name their reasons: two pieces of one sphere and
/// a piece not convex in its planes as S9e.4b.3c's, the flush box and the
/// touching box as S9's.
#[test]
fn refusals_name_their_reasons() {
    for (name, run) in runs() {
        let reason = match run {
            Err(Error::Degenerate(m) | Error::OutOfDomain(m)) => *m,
            _ => "",
        };
        if name.starts_with("one_sphere") {
            assert_eq!(
                reason, "faces of both inputs on one sphere (S9e.4b.3c)",
                "{name}"
            );
        }
        if name.starts_with("bitten_box") {
            assert!(
                matches!(run, Err(Error::OutOfDomain(_))) && reason.contains("S9e.4b.3c"),
                "{name}: {reason}"
            );
        }
        if name.starts_with("tilt_flush") {
            assert_eq!(
                reason, "two faces within the resolution of one plane",
                "{name}"
            );
        }
        if name.starts_with("half_touch") {
            assert!(reason.contains("tangency"), "{name}: {reason}");
        }
    }
}

/// Each Boolean's history is complete over its inputs, the imported ones'
/// stored ids among them: every entity of every input resolved, and no
/// relation from an id outside the inputs (a primitive's or a hull's).
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
    let first: std::collections::BTreeMap<&str, _> =
        runs().iter().map(|(n, r)| (n.as_str(), r)).collect();
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

/// Both inputs moved rigidly keep the reference's volumes, translated and
/// turned: the piece read again off its moved stored topology.
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    for name in [
        "octant_box_cut",
        "box_octant_common",
        "tilt_rod_fuse",
        "lune_ball_cut",
        "half_rod_cut",
        "zone_box_common",
        "pieces_fuse",
    ] {
        let case = all.iter().find(|c| c.name == name).unwrap();
        for (k, motion) in [shift, turn].into_iter().enumerate() {
            let (a, b, _, _) = protocol::run_first(case).unwrap();
            let a = a.transform_with(OperationId(801), motion).unwrap().0;
            let b = b.transform_with(OperationId(802), motion).unwrap().0;
            let out = match case.op.as_str() {
                "fuse" => a.fuse(case.operation, &b),
                "cut" => a.cut(case.operation, &b),
                _ => a.common(case.operation, &b),
            }
            .unwrap_or_else(|e| panic!("{name} {k}: {e}"))
            .0;
            let (count, [v, _]) = expect[name].1.unwrap();
            let got: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
            assert_eq!(out.len(), count, "{name} {k}");
            assert!((got - v).abs() <= 1e-9 * v, "{name} {k}: {got} for {v}");
        }
    }
}

/// A body file's imported solid under an operation.
fn body(name: &str, op: u64) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(op), topology, resolution).map(|(s, _)| s)
}

/// Every piece imports (its model matched on import), its stored mass its
/// closed form where it has one, its stored vertices on its boundary and
/// points in and off it classified; the ball below one plane is S9e.4a's
/// cap and the ball less a box's corner S9e.4b.3c's.
#[test]
fn imported_bodies_are_pieces() {
    use std::f64::consts::PI;
    // Each with a point inside (along its corner's or wedge's middle) and
    // that point's mirror in its centre outside.
    let bodies = [
        (
            "octant",
            Some(PI * 125.0 / 6.0),
            [5.0, 5.0, 4.0],
            [5.48, 7.41, 3.52],
        ),
        (
            "octant_tilt",
            Some(PI * 125.0 / 6.0),
            [5.25, 8.5, 1.75],
            [5.41, 6.74, 3.51],
        ),
        ("upper", None, [3.0, 4.0, 1.0], [5.7, 3.84, 3.22]),
        (
            "lune",
            Some(PI * 64.0 / 3.0),
            [5.0, 5.0, 5.0],
            [4.21, 6.73, 4.37],
        ),
        (
            "half",
            Some(PI * 128.0 / 3.0),
            [5.0, 5.0, 5.0],
            [5.22, 3.22, 5.89],
        ),
        ("zone_wedge", None, [5.0, 5.0, 4.0], [3.23, 6.41, 5.06]),
    ];
    for (name, volume, centre, inside) in bodies {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let m = s.mass_properties();
        if let Some(volume) = volume {
            assert!(
                (m.volume - volume).abs() <= 1e-9 * volume,
                "{name}: {} for {volume}",
                m.volume
            );
        }
        let [x, y, z] = inside;
        assert_eq!(
            s.classify(Point3::new(x, y, z)).unwrap(),
            Location::Inside,
            "{name}"
        );
        let mirror = [0, 1, 2].map(|i| 2.0 * centre[i] - inside[i]);
        for out in [[x + 50.0, y, z], mirror] {
            assert_eq!(
                s.classify(Point3::new(out[0], out[1], out[2])).unwrap(),
                Location::Outside,
                "{name} {out:?}"
            );
        }
        for v in s.topology().vertices() {
            assert_eq!(
                s.classify(v.position).unwrap(),
                Location::Boundary,
                "{name}"
            );
        }
        // Two imports keep their stored ids apart.
        let t = body(name, 92).unwrap();
        let a: std::collections::BTreeSet<_> = s.topology().ids().map(|(i, _)| i).collect();
        assert!(t.topology().ids().all(|(i, _)| !a.contains(&i)), "{name}");
    }
    body("octant_low", 91).expect("S9e.4a's cap");
    match body("bitten", 91) {
        Err(Error::OutOfDomain(m)) => assert!(m.contains("S9e.4b.3c"), "{m}"),
        other => panic!("bitten: {:?}", other.map(|_| ())),
    }
}

fn tolerance() -> Tolerance {
    Tolerance::new(1e-7, 1e-12).unwrap()
}

fn frame(o: [f64; 3], n: [f64; 3], x: [f64; 3]) -> Frame3 {
    let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
    Frame3::new(Point3::new(o[0], o[1], o[2]), v(n), v(x), tolerance()).unwrap()
}

/// A box `[0, w] x [0, d]` on a frame over `[0, h]`.
fn cuboid(f: Frame3, w: f64, d: f64, h: f64) -> Solid {
    let pts = [(0.0, 0.0), (w, 0.0), (w, d), (0.0, d)].map(|(x, y)| Point2::new(x, y));
    let profile = Profile::new(
        Boundary::polygon(pts.to_vec(), tolerance()).unwrap(),
        vec![],
        tolerance(),
    )
    .unwrap();
    Solid::extrude_with(OperationId(40), profile, f, 0.0, h)
        .unwrap()
        .0
}

/// A cylinder's and a frustum's pieces of an oblique plane (S8a.2, S8d.2)
/// and a zone's halves (S8c.2), the kernel's own split pieces, their
/// topologies imported: each a plane piece (its model matched on import),
/// its Booleans with a box and a ball obeying the pair identities
/// `V(A u B) + V(A n B) = V(A) + V(B)` and `V(A - B) = V(A) - V(A n B)`,
/// its common with a box holding it its own volume.
#[test]
fn kernel_split_pieces_imported() {
    let up = [0.0, 0.0, 1.0];
    let ex = [1.0, 0.0, 0.0];
    let tol = tolerance();
    let cylinder = Solid::cylinder_with(
        OperationId(1),
        frame([5.0, 5.0, -1.0], up, ex),
        3.0,
        0.0,
        8.0,
        tol,
    )
    .unwrap()
    .0;
    let frustum = Solid::cone_with(
        OperationId(3),
        frame([5.0, 5.0, 0.0], up, ex),
        3.0,
        1.0,
        5.0,
        tol,
    )
    .unwrap()
    .0;
    let zone = Solid::sphere_with(
        OperationId(5),
        frame([5.0, 5.0, 3.0], up, ex),
        4.0,
        -0.5,
        0.75,
        tol,
    )
    .unwrap()
    .0;
    let splits = [
        (
            &cylinder,
            frame([5.0, 5.0, 3.0], [0.6, 0.0, 0.8], [0.0, 1.0, 0.0]),
        ),
        (&frustum, frame([5.0, 5.0, 2.5], [0.0, 0.6, 0.8], ex)),
        (&zone, frame([5.0, 5.0, 3.0], ex, [0.0, 1.0, 0.0])),
    ];
    let partners = [
        cuboid(frame([3.0, 3.0, 1.0], up, ex), 4.0, 3.0, 2.5),
        Solid::sphere_with(
            OperationId(41),
            frame([6.0, 6.5, 1.75], up, ex),
            1.5,
            -std::f64::consts::FRAC_PI_2,
            std::f64::consts::FRAC_PI_2,
            tol,
        )
        .unwrap()
        .0,
    ];
    let holder = cuboid(frame([-5.0, -5.0, -5.0], up, ex), 20.0, 20.0, 20.0);
    let volume = |r: rusty_occt::Result<(Vec<Solid>, history::History)>| -> f64 {
        r.unwrap()
            .0
            .iter()
            .map(|s| s.mass_properties().volume)
            .sum()
    };
    for (k, (solid, plane)) in splits.into_iter().enumerate() {
        let (pieces, _) = solid.split_by_plane(OperationId(2), plane).unwrap();
        assert_eq!(pieces.len(), 2, "split {k}");
        for (j, (_, piece)) in pieces.iter().enumerate() {
            let (s, _) = Solid::imported_with(
                OperationId(31),
                piece.topology().clone(),
                piece.resolution(),
            )
            .unwrap_or_else(|e| panic!("split {k} piece {j}: {e}"));
            let va = s.mass_properties().volume;
            let all = volume(s.common(OperationId(50), &holder));
            assert!(
                (all - va).abs() <= 1e-9 * va,
                "split {k} piece {j}: {all} for {va}"
            );
            for (i, b) in partners.iter().enumerate() {
                let vb = b.mass_properties().volume;
                let f = volume(s.fuse(OperationId(50), b));
                let c = volume(s.cut(OperationId(50), b));
                let m = volume(s.common(OperationId(50), b));
                let size = va + vb;
                assert!(
                    (f + m - va - vb).abs() <= 1e-9 * size,
                    "split {k} piece {j} partner {i}"
                );
                assert!(
                    (c - (va - m)).abs() <= 1e-9 * size,
                    "split {k} piece {j} partner {i}"
                );
                assert!(
                    m > 0.0 && m < va.min(vb),
                    "split {k} piece {j} partner {i}: {m}"
                );
            }
        }
    }
}

/// The octant against a whole torus across its sphere face (S9d.4b.2's
/// meeting on the piece's sphere): the pair identities, and the common
/// within both.
#[test]
fn a_piece_against_a_torus() {
    let s = body("octant", 91).unwrap();
    let tau = std::f64::consts::TAU;
    let torus = Solid::torus_with(
        OperationId(92),
        frame([6.0, 7.25, 3.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        2.0,
        0.75,
        0.0,
        tau,
        tau,
        tolerance(),
    )
    .unwrap()
    .0;
    let volume = |r: rusty_occt::Result<(Vec<Solid>, history::History)>| -> f64 {
        let out = r.unwrap().0;
        assert_eq!(out.len(), 1);
        out.iter().map(|s| s.mass_properties().volume).sum()
    };
    let (va, vb) = (s.mass_properties().volume, torus.mass_properties().volume);
    let f = volume(s.fuse(OperationId(93), &torus));
    let c = volume(s.cut(OperationId(93), &torus));
    let m = volume(s.common(OperationId(93), &torus));
    assert!(
        (f + m - va - vb).abs() <= 1e-9 * (va + vb),
        "{f} {m} {va} {vb}"
    );
    assert!((c - (va - m)).abs() <= 1e-9 * va, "{c} {m} {va}");
    assert!(m > 0.0 && m < vb, "{m}");
}

/// The kernel stores the frames the reference swept, bit for bit (the
/// solids built from rows; the imported ones' are the converter's).
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-imported-pieces-frames.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let w: Vec<&str> = row.split('\t').collect();
        let case = cases.iter().find(|c| c.name == w[0]).unwrap();
        let (which, axis) = w[1].split_once(' ').unwrap();
        let k: usize = which.trim_start_matches("solid").parse().unwrap();
        let spec = match k {
            0 => &case.object,
            1 => &case.tool,
            _ => &case.then.as_ref().unwrap().third,
        };
        let frame = protocol::build(spec).frame();
        let v = match axis {
            "n" => frame.normal(),
            "x" => frame.x(),
            _ => frame.y(),
        };
        let got = [v.x, v.y, v.z].map(hex).join(" ");
        assert_eq!(got, w[2], "{} {}", w[0], w[1]);
    }
}
