//! S9e.4b.3b: the kernel's own plane pieces (S8's split pieces: a prism's
//! oblique piece, `Clipped`, and a cone's, a zone's or a torus's piece,
//! `Half`) given to Booleans with curved faces, against the independent
//! reference (`fixtures/boolean-split-pieces-*` from
//! `tools/generate_split_pieces_boolean_fixtures.py`): each case's Boolean
//! (or chain) with its split inputs made by `Solid::split_by_plane` (the
//! protocol's `split` row), decided as the given model of the split's
//! primitive common the half-spaces of its planes. Each case runs once (on
//! a few threads) for the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{
    Boundary, Error, Frame3, Location, Point2, Point3, Profile, RigidTransform, Side, Solid,
    Tolerance, Vec3,
};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-split-pieces-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-split-pieces-expected.tsv")
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

/// The declared refusals name their reasons: the frustum's half through its
/// axis its apex, the box on the cut plane one plane (the zone's half and
/// its own ball, S9e.4b.3c.1's, evaluate).
#[test]
fn refusals_name_their_reasons() {
    for (name, run) in runs() {
        let reason = match run {
            Err(Error::Degenerate(m) | Error::OutOfDomain(m)) => *m,
            _ => "",
        };
        if name.starts_with("cone_axis") {
            assert_eq!(reason, "a plane through a cone's apex", "{name}");
        }
        if name.starts_with("cyl_flush") {
            assert_eq!(
                reason, "two faces within the resolution of one plane",
                "{name}"
            );
        }
    }
}

/// Each Boolean's history is complete over its inputs, the split pieces'
/// ids among them: every entity of every input resolved, and no relation
/// from an id outside the inputs (a primitive's or a hull's).
#[test]
fn histories_are_complete_over_the_pieces_ids() {
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
/// turned: the piece's model built again from the moved piece. (A plane
/// exactly along a rod's axis, `cyl_rod`'s and `block_rod`'s, is within
/// rounding of it once turned: S9's `Degenerate`.)
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    for name in [
        "cyl_box_cut",
        "box_cyl_common",
        "cyl_ball_fuse",
        "dee_rod_cut",
        "frustum_ball_cut",
        "zone_box_common",
        "zone_cone_fuse",
        "cut_rod_fuse",
        "band_ball_common",
        "pair_cut",
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

/// The kernel stores the frames the reference swept, bit for bit: the
/// solids built from rows and the splits' planes.
#[test]
fn stored_frames_are_the_reference_inputs() {
    let hex = |x: f64| format!("{:016x}", x.to_bits());
    let cases = cases();
    for row in include_str!("../../fixtures/boolean-split-pieces-frames.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let w: Vec<&str> = row.split('\t').collect();
        let case = cases.iter().find(|c| c.name == w[0]).unwrap();
        let (which, axis) = w[1].split_once(' ').unwrap();
        let (kind, k) = which.split_at(5);
        let k: usize = k.parse().unwrap();
        let spec = match k {
            0 => &case.object,
            1 => &case.tool,
            _ => &case.then.as_ref().unwrap().third,
        };
        let frame = if kind == "split" {
            let (f, _) = spec.split.unwrap();
            Frame3::new(
                Point3::new(f[0], f[1], f[2]),
                Vec3::new(f[3], f[4], f[5]),
                Vec3::new(f[6], f[7], f[8]),
                spec.tolerance,
            )
            .unwrap()
        } else {
            protocol::build(spec).frame()
        };
        let v = match axis {
            "n" => frame.normal(),
            "x" => frame.x(),
            _ => frame.y(),
        };
        let got = [v.x, v.y, v.z].map(hex).join(" ");
        assert_eq!(got, w[2], "{} {}", w[0], w[1]);
    }
}

fn tolerance() -> Tolerance {
    Tolerance::new(1e-7, 1e-12).unwrap()
}

fn frame(o: [f64; 3], n: [f64; 3], x: [f64; 3]) -> Frame3 {
    let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
    Frame3::new(Point3::new(o[0], o[1], o[2]), v(n), v(x), tolerance()).unwrap()
}

/// A prism of a polygon on a frame over `[0, h]`.
fn prism(op: u64, f: Frame3, pts: &[(f64, f64)], h: f64) -> Solid {
    let pts: Vec<Point2> = pts.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    let profile = Profile::new(
        Boundary::polygon(pts, tolerance()).unwrap(),
        vec![],
        tolerance(),
    )
    .unwrap();
    Solid::extrude_with(OperationId(op), profile, f, 0.0, h)
        .unwrap()
        .0
}

/// A box `[0, w] x [0, d]` on a frame over `[0, h]`.
fn cuboid(op: u64, f: Frame3, w: f64, d: f64, h: f64) -> Solid {
    prism(op, f, &[(0.0, 0.0), (w, 0.0), (w, d), (0.0, d)], h)
}

type Outcome = rusty_occt::Result<(Vec<Solid>, history::History)>;

fn volume(r: Outcome) -> f64 {
    r.unwrap()
        .0
        .iter()
        .map(|s| s.mass_properties().volume)
        .sum()
}

/// The pair identities `V(A u B) + V(A n B) = V(A) + V(B)` and `V(A - B) =
/// V(A) - V(A n B)`, the common's volume returned.
fn identities(what: &str, a: &Solid, b: &Solid) -> f64 {
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let f = volume(a.fuse(OperationId(50), b));
    let c = volume(a.cut(OperationId(50), b));
    let m = volume(a.common(OperationId(50), b));
    let size = va + vb;
    assert!((f + m - va - vb).abs() <= 1e-9 * size, "{what}: {f} {m}");
    assert!((c - (va - m)).abs() <= 1e-9 * size, "{what}: {c} {m}");
    m
}

/// A box holding a solid: its common with any piece of it is the piece.
fn holder() -> Solid {
    cuboid(
        70,
        frame([-20.0, -20.0, -20.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        60.0,
        60.0,
        60.0,
    )
}

/// A U profile split by a plane across both arms' tops: two pieces above
/// it, each a solid of several of its primitive common the plane's
/// half-space (S9e.2's sorting by solid), against a rod through one arm and
/// a box holding both.
#[test]
fn several_pieces_on_one_side() {
    let up = [0.0, 0.0, 1.0];
    let ex = [1.0, 0.0, 0.0];
    let u = [
        (0.0, 0.0),
        (9.0, 0.0),
        (9.0, 6.0),
        (6.0, 6.0),
        (6.0, 2.0),
        (3.0, 2.0),
        (3.0, 6.0),
        (0.0, 6.0),
    ];
    // Planar pieces, the curved engine's against the rod's curved wall.
    let solid = prism(9, frame([0.0, 0.0, 0.0], up, ex), &u, 4.0);
    let (pieces, _) = solid
        .split_by_plane(OperationId(2), frame([4.5, 4.0, 3.0], [0.0, 3.0, 4.0], ex))
        .unwrap();
    let above: Vec<&Solid> = pieces
        .iter()
        .filter(|(s, _)| *s == Side::Above)
        .map(|(_, p)| p)
        .collect();
    assert_eq!(above.len(), 2);
    let rod = Solid::cylinder_with(
        OperationId(42),
        frame([0.0, 5.25, 3.5], ex, [0.0, 1.0, 0.0]),
        0.375,
        -1.0,
        10.0,
        tolerance(),
    )
    .unwrap()
    .0;
    for (k, piece) in above.iter().enumerate() {
        let m = identities(&format!("U {k}"), piece, &rod);
        assert!(m > 0.0, "U {k}: {m}");
        let all = volume(piece.common(OperationId(50), &holder()));
        let own = piece.mass_properties().volume;
        assert!((all - own).abs() <= 1e-9 * own, "U {k}: {all} for {own}");
    }
}

/// A hemisphere's halves by a plane through its axis: the section through
/// its pole, on the whole sphere's model with its base's plane (S9e.4b.3a's
/// model of a sphere's piece), against a box across both.
#[test]
fn a_caps_halves_through_its_pole() {
    let up = [0.0, 0.0, 1.0];
    let ex = [1.0, 0.0, 0.0];
    let hp = std::f64::consts::FRAC_PI_2;
    let cap = Solid::sphere_with(
        OperationId(5),
        frame([5.0, 5.0, 3.0], up, ex),
        4.0,
        0.0,
        hp,
        tolerance(),
    )
    .unwrap()
    .0;
    let (pieces, _) = cap
        .split_by_plane(OperationId(2), frame([5.0, 5.0, 3.0], ex, [0.0, 1.0, 0.0]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    let b = cuboid(43, frame([3.25, 3.5, 4.0], up, ex), 4.5, 3.0, 4.5);
    let mut commons = Vec::new();
    for (k, (_, piece)) in pieces.iter().enumerate() {
        let own = piece.mass_properties().volume;
        assert!((own - std::f64::consts::PI * 64.0 / 3.0).abs() <= 1e-9 * own);
        commons.push(identities(&format!("cap {k}"), piece, &b));
    }
    // The halves' commons add up to the cap's.
    let whole = volume(cap.common(OperationId(50), &b));
    assert!(
        (commons[0] + commons[1] - whole).abs() <= 1e-9 * whole,
        "{commons:?} for {whole}"
    );
}

/// A split piece against an imported piece (S9e.4b.3a's octant) and a given
/// result (a box less a rod), as object and tool: the pair identities.
#[test]
fn pieces_against_other_bodies() {
    let up = [0.0, 0.0, 1.0];
    let ex = [1.0, 0.0, 0.0];
    let cylinder = Solid::cylinder_with(
        OperationId(1),
        frame([6.25, 6.5, -0.5], up, ex),
        2.5,
        0.0,
        8.0,
        tolerance(),
    )
    .unwrap()
    .0;
    let (pieces, _) = cylinder
        .split_by_plane(
            OperationId(2),
            frame([6.25, 6.5, 3.0], [3.0, 0.0, 4.0], [0.0, 1.0, 0.0]),
        )
        .unwrap();
    let low = &pieces[0].1;
    let (topology, resolution) = protocol::imported_topology("imported/octant.brep");
    let octant = Solid::imported_with(OperationId(91), topology, resolution)
        .unwrap()
        .0;
    let m = identities("octant", low, &octant);
    assert!(m > 0.0);
    let block = cuboid(44, frame([2.5, 3.25, 0.5], up, ex), 3.5, 3.5, 3.0);
    let rod = Solid::cylinder_with(
        OperationId(45),
        frame([4.25, 0.0, 2.25], [0.0, 1.0, 0.0], ex),
        0.625,
        0.0,
        10.0,
        tolerance(),
    )
    .unwrap()
    .0;
    let given = block.cut(OperationId(46), &rod).unwrap().0.remove(0);
    let m = identities("given", low, &given);
    assert!(m > 0.0);
    let m = identities("given swapped", &given, low);
    assert!(m > 0.0);
}

/// A torus's spiric piece and a spline prism's piece stay refused, named.
#[test]
fn refused_pieces() {
    let up = [0.0, 0.0, 1.0];
    let ex = [1.0, 0.0, 0.0];
    let tau = std::f64::consts::TAU;
    let hp = std::f64::consts::FRAC_PI_2;
    let torus = Solid::torus_with(
        OperationId(7),
        frame([5.0, 5.0, 3.0], up, ex),
        3.0,
        1.0,
        0.0,
        tau,
        tau,
        tolerance(),
    )
    .unwrap()
    .0;
    let (pieces, _) = torus
        .split_by_plane(OperationId(2), frame([5.0, 5.0, 3.25], [0.0, 0.6, 0.8], ex))
        .unwrap();
    let ball = Solid::sphere_with(
        OperationId(41),
        frame([7.0, 5.5, 3.5], up, ex),
        1.0,
        -hp,
        hp,
        tolerance(),
    )
    .unwrap()
    .0;
    for (_, p) in &pieces {
        match p.common(OperationId(50), &ball) {
            Err(Error::OutOfDomain(m)) => assert!(m.contains("(refused, S9e.4b.3b)"), "{m}"),
            other => panic!("a spiric piece: {:?}", other.map(|_| ())),
        }
    }
}

/// A frustum's pieces by a plane through its virtual apex's rounded height,
/// across its axis (the split target's replay): the split takes the plane
/// within the resolution of the apex and cuts rulings; the model's exact
/// plane misses the apex by less than the resolution, its section's
/// branches within rounding of those rulings, and S9d.3a refuses it,
/// `Degenerate`, for both pieces (the piece below it once failed its
/// model's assembly as an open Boolean, `InvalidTopology`).
#[test]
fn pieces_by_a_plane_within_rounding_of_a_virtual_apex() {
    let tol = Tolerance::default();
    let (bottom, top, height) = (3.5, 0.25, 1.75);
    let f = Frame3::new(
        Point3::new(0.1875, 6.0, -8.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let frustum = Solid::cone_with(OperationId(1), f, bottom, top, height, tol)
        .unwrap()
        .0;
    let apex = f.point(Point2::new(0.0, 0.0), -bottom * height / (top - bottom));
    let plane = Frame3::new(
        apex,
        Vec3::new(1.0, -0.75, -0.5),
        Vec3::new(0.0, 1.0, 0.0),
        tol,
    )
    .unwrap();
    let (pieces, _) = frustum.split_by_plane(OperationId(2), plane).unwrap();
    assert_eq!(pieces.len(), 2);
    let all: f64 = pieces.iter().map(|(_, p)| p.mass_properties().volume).sum();
    let whole = frustum.mass_properties().volume;
    assert!((all - whole).abs() <= 1e-9 * whole, "{all} for {whole}");
    let refused = |r: Outcome| match r {
        Err(Error::Degenerate(m)) => {
            assert_eq!(m, "a plane within the resolution of a cone's apex")
        }
        other => panic!("{:?}", other.map(|_| ())),
    };
    for (_, piece) in &pieces {
        refused(piece.common(OperationId(4), &holder()));
        let c = piece.mass_properties().centroid;
        let turned =
            Frame3::new(c, Vec3::new(1.0, 2.0, 2.0), Vec3::new(2.0, 1.0, -2.0), tol).unwrap();
        let other = cuboid(5, turned, 3.5, 3.5, 3.5);
        refused(piece.fuse(OperationId(6), &other));
        refused(piece.cut(OperationId(6), &other));
        refused(piece.common(OperationId(6), &other));
    }
}

/// Each fixture's split piece classifies as its construction: its stored
/// vertices on its boundary, a point far off outside; its common with a
/// box holding it is itself.
#[test]
fn split_pieces_are_their_own() {
    let mut seen = std::collections::BTreeSet::new();
    for case in cases() {
        for spec in [&case.object, &case.tool] {
            let Some((f, _)) = spec.split else { continue };
            if !seen.insert(format!("{f:?}")) {
                continue;
            }
            let piece = protocol::build(spec);
            for v in piece.topology().vertices() {
                assert_eq!(
                    piece.classify(v.position).unwrap(),
                    Location::Boundary,
                    "{}",
                    case.name
                );
            }
            assert_eq!(
                piece.classify(Point3::new(90.0, 0.0, 0.0)).unwrap(),
                Location::Outside
            );
            let own = piece.mass_properties().volume;
            match piece.common(OperationId(50), &holder()) {
                Ok((out, _)) => {
                    let all: f64 = out.iter().map(|s| s.mass_properties().volume).sum();
                    assert!((all - own).abs() <= 1e-9 * own, "{}: {all}", case.name);
                }
                // The frustum's half through its axis: its apex.
                Err(Error::Degenerate(m)) => {
                    assert_eq!(m, "a plane through a cone's apex", "{}", case.name)
                }
                Err(e) => panic!("{}: {e}", case.name),
            }
        }
    }
    assert_eq!(seen.len(), 8);
}

/// The engine's rules the split target's pieces reached (each also an
/// imported piece's, S9e.4b.3a's model): a zone's rim as a plane's section
/// of its whole sphere stored as a ring the other way round; a cap's half
/// through its pole against a box across the seam, its hole lifted with
/// the loop through the pole; a frustum's oblique piece, its section stored
/// as a ring the other way round about the axis. Each piece's common with
/// a box holding it is itself, and its fuse, cut and common with a turned
/// box obey the pair identities.
#[test]
fn engine_rules_the_split_target_reached() {
    let tol = Tolerance::default();
    let f = |o: [f64; 3], n: [f64; 3], x: [f64; 3]| {
        Frame3::new(
            Point3::new(o[0], o[1], o[2]),
            Vec3::new(n[0], n[1], n[2]),
            Vec3::new(x[0], x[1], x[2]),
            tol,
        )
        .unwrap()
    };
    let zone = Solid::sphere_with(
        OperationId(1),
        f([-1.375, -8.0, 1.25], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        4.25,
        -1.375,
        0.75,
        tol,
    )
    .unwrap()
    .0;
    let zone_plane = f(
        [-0.3125, -9.46875, -0.015625],
        [-2.0, 1.0, -2.5],
        [0.0, 1.0, 0.0],
    );
    let tilted = f([-0.6875, 1.875, 4.4375], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let cap = Solid::sphere_with(
        OperationId(1),
        tilted,
        3.0,
        0.25,
        std::f64::consts::FRAC_PI_2,
        tol,
    )
    .unwrap()
    .0;
    let cap_plane = f(
        [-0.6875, 1.875, 4.4375],
        [1.0, 0.4, -0.29999999999999993],
        [0.0, 1.0, 0.0],
    );
    let frustum = Solid::cone_with(
        OperationId(1),
        f([-7.0, 7.9375, 7.9375], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
        3.75,
        1.5,
        4.25,
        tol,
    )
    .unwrap()
    .0;
    let frustum_plane = f(
        [-5.015625, 9.921875, 9.921875],
        [0.0, -3.0, -2.5],
        [1.0, 0.0, 0.0],
    );
    let holder = Solid::box_at(
        Point3::new(-40.0, -40.0, -40.0),
        Vec3::new(80.0, 80.0, 80.0),
        tol,
    )
    .unwrap();
    for (name, solid, plane) in [
        ("zone", &zone, zone_plane),
        ("cap", &cap, cap_plane),
        ("frustum", &frustum, frustum_plane),
    ] {
        let (pieces, _) = solid.split_by_plane(OperationId(2), plane).unwrap();
        assert_eq!(pieces.len(), 2, "{name}");
        for (k, (_, piece)) in pieces.iter().enumerate() {
            let (imported, _) = Solid::imported_with(
                OperationId(31),
                piece.topology().clone(),
                piece.resolution(),
            )
            .unwrap();
            for (how, p) in [("split", piece), ("imported", &imported)] {
                let what = format!("{name} {k} {how}");
                let own = p.mass_properties().volume;
                let all = volume(p.common(OperationId(4), &holder));
                assert!((all - own).abs() <= 1e-9 * own, "{what}: {all} for {own}");
                // A box turned about the piece's centroid, a corner there.
                let c = p.mass_properties().centroid;
                let turned =
                    Frame3::new(c, Vec3::new(1.0, 2.0, 2.0), Vec3::new(2.0, 1.0, -2.0), tol)
                        .unwrap();
                identities(&what, p, &cuboid(5, turned, 3.5, 3.5, 3.5));
            }
        }
    }
}
