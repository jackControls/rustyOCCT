//! S9e.4b.4c.1: imported polyhedra (bodies of plane faces and line edges OCCT
//! wrote to `.brep` files, decided on their stored vertices) against curved
//! faces, their results given to further Booleans, and imported polyhedra
//! with a cavity, against the independent reference
//! (`fixtures/boolean-polyhedra-curved-*` from
//! `tools/generate_polyhedra_curved_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/`); the kernel's own polyhedra written by its writer,
//! read back and imported. Each case runs once (on a few threads) for the
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
    protocol::cases(include_str!(
        "../../fixtures/boolean-polyhedra-curved-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-polyhedra-curved-expected.tsv")
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
        "pyramid_rod_cut",
        "wedge_ball_common",
        "octa_cone_fuse",
        "notched_ball_cut",
        "tetra_ball_fuse",
        "hollow_ball_fuse",
        "cavity_slab_cut",
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

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Boundary {
    let points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        .map(|(x, y)| Point2::new(x, y))
        .to_vec();
    Boundary::polygon(points, tolerance()).unwrap()
}

fn prism(outer: Boundary, f: Frame3, h0: f64, h1: f64, op: u64) -> Solid {
    let p = Profile::new(outer, vec![], tolerance()).unwrap();
    Solid::extrude_with(OperationId(op), p, f, h0, h1)
        .unwrap()
        .0
}

fn rod(f: Frame3, r: f64, h0: f64, h1: f64, op: u64) -> Solid {
    prism(
        Boundary::circle(Point2::new(0.0, 0.0), r, tolerance()).unwrap(),
        f,
        h0,
        h1,
        op,
    )
}

fn ball(o: [f64; 3], r: f64, op: u64) -> Solid {
    let half = std::f64::consts::FRAC_PI_2;
    Solid::sphere_with(
        OperationId(op),
        frame(o, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        r,
        -half,
        half,
        tolerance(),
    )
    .unwrap()
    .0
}

fn volume(out: &[Solid]) -> f64 {
    out.iter().map(|s| s.mass_properties().volume).sum()
}

/// Every body imports as a polyhedron on its stored vertices (the hollow box
/// and this step's cavity with their cavities), its volume its
/// construction's, its stored vertices on its boundary, points in its
/// material inside and points in its cavity or off it outside.
#[test]
fn imported_bodies_are_their_meshes() {
    let turn30 = frame(
        [1.0, 2.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.8660254037844386, 0.5, 0.0],
    );
    let tilt = frame([1.75, 6.375, 3.0], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let at = |f: Frame3, u: f64, v: f64, w: f64| f.point(Point2::new(u, v), w).to_array();
    type Points = Vec<[f64; 3]>;
    let bodies: [(&str, f64, Points, Points); 8] = [
        (
            "pyramid",
            128.0,
            vec![[5.0, 4.6, 2.3]],
            vec![[5.0, 4.6, 30.0]],
        ),
        (
            "truncated",
            112.0,
            vec![[3.8748, 4.0207, 4.0]],
            vec![[3.8748, 40.0, 4.0]],
        ),
        (
            "wedge",
            368.0 / 3.0,
            vec![[4.6756, 4.9047, 3.587]],
            vec![[4.6756, 4.9047, 30.0]],
        ),
        (
            "tetra",
            56.0,
            vec![[3.6896, 4.5736, 1.0]],
            vec![[3.6896, 4.5736, 6.0]],
        ),
        (
            "octa",
            256.0 / 3.0,
            vec![[5.0, 5.0, 4.0]],
            vec![[8.5, 8.5, 4.0]],
        ),
        (
            "notched",
            394.48066666666665,
            vec![[4.9369, 4.9531, 1.9886]],
            vec![[9.5, 9.5, 3.5]],
        ),
        (
            "hollow",
            936.0,
            vec![[1.5, 5.0, 5.0]],
            vec![[5.0, 5.0, 5.0], [5.0, 5.0, 12.0]],
        ),
        (
            "poly_cavity",
            465.0,
            vec![at(turn30, 1.0, 1.0, 1.0), at(turn30, 9.0, 7.0, 5.0)],
            vec![at(tilt, 1.5, 1.25, 1.0), at(turn30, 5.0, 4.0, 7.0)],
        ),
    ];
    for (name, v, inside, outside) in bodies {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let m = s.mass_properties();
        assert!(
            (m.volume - v).abs() <= 1e-9 * v,
            "{name}: {} for {v}",
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
fn reimported(s: &Solid, op: u64) -> Solid {
    use rusty_occt::occt_brep::{import, read, write};
    let text = write(s.topology(), s.resolution().linear()).unwrap();
    let doc = read(&text).unwrap();
    let [solid] = <[_; 1]>::try_from(import(&doc).solids).unwrap();
    Solid::imported_with(OperationId(op), solid.result.unwrap(), solid.tolerance)
        .unwrap()
        .0
}

/// The kernel's own polyhedra written, read back and imported (no prism: a
/// box less a box on a turned frame across its corner, and a box less a box
/// on another turned frame inside it, a cavity), each given to Booleans with
/// a leaning rod and a ball: the kernel's own results' volumes (its results
/// given to the curved engine on their constructions), or the same refusal
/// (the hollow box fused with the rod through its cavity alone).
#[test]
fn kernel_polyhedra_imported_against_curved_faces() {
    let one = |what: &str, r: rusty_occt::Result<(Vec<Solid>, History)>| {
        let (out, _) = r.unwrap_or_else(|e| panic!("{what}: {e}"));
        assert_eq!(out.len(), 1, "{what}");
        out.into_iter().next().unwrap()
    };
    let world = frame([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let turned = frame([8.2, 7.9, 2.4], [1.0, 2.0, 2.0], [2.0, -2.0, 1.0]);
    let notched = one(
        "notched",
        prism(square(0.0, 0.0, 10.0, 10.0), world, 0.0, 4.0, 1).cut(
            OperationId(2),
            &prism(square(0.0, 0.0, 6.0, 6.0), turned, 0.0, 6.0, 3),
        ),
    );
    let inner = frame([2.5, 2.0, 2.0], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let hollow = one(
        "hollow",
        prism(square(0.0, 0.0, 9.0, 7.0), world, 0.0, 5.0, 4).cut(
            OperationId(5),
            &prism(square(0.0, 0.0, 3.0, 2.5), inner, 0.0, 1.5, 6),
        ),
    );
    assert_eq!(hollow.topology().shells().len(), 4);
    let lean = frame([3.25, 3.5, -1.0], [1.0, 2.0, 8.0], [8.0, 0.0, -1.0]);
    let partners = [
        ("rod", rod(lean, 1.25, 0.0, 9.0, 7)),
        ("ball", ball([8.125, 3.75, 3.5], 1.75, 8)),
    ];
    let near = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(1.0);
    let mut evaluated = Vec::new();
    for (name, s) in [("notched", notched), ("hollow", hollow)] {
        let i = reimported(&s, 50);
        let (v0, v1) = (s.mass_properties().volume, i.mass_properties().volume);
        assert!(near(v0, v1), "{name}: {v1} for {v0}");
        for (partner, p) in &partners {
            for (direct, again) in [
                (s.fuse(OperationId(60), p), i.fuse(OperationId(61), p)),
                (s.cut(OperationId(60), p), i.cut(OperationId(61), p)),
                (s.common(OperationId(60), p), i.common(OperationId(61), p)),
                (p.cut(OperationId(60), &s), p.cut(OperationId(61), &i)),
            ] {
                match (direct, again) {
                    (Ok(d), Ok(a)) => {
                        let (d, a) = (volume(&d.0), volume(&a.0));
                        assert!(near(d, a), "{name} {partner}: {a} for {d}");
                        evaluated.push(format!("{name} {partner}"));
                    }
                    (Err(d), Err(a)) => assert_eq!(d, a, "{name} {partner}"),
                    (d, a) => panic!(
                        "{name} {partner}: {:?} and {:?}",
                        d.map(|_| ()),
                        a.map(|_| ())
                    ),
                }
            }
        }
    }
    // Every pair but the hollow box's fuse with the rod through its
    // cavity, its void a ring the validator's rays leave undecided.
    assert_eq!(evaluated.len(), 15, "{evaluated:?}");
}

/// A cavity split in two by a slab (each cavity a shell and a void region
/// of its own), the result given to a ball crossing the outer wall: the
/// identities hold and the solids validate with both cavities.
#[test]
fn cavities_split_and_given_again() {
    let hollow = body("hollow", 91).unwrap();
    let slab = prism(
        square(-1.0, -1.0, 11.0, 11.0),
        frame([0.0, 0.0, 4.5], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        0.0,
        1.0,
        92,
    );
    let (fused, _) = hollow.fuse(OperationId(93), &slab).unwrap();
    let [fused] = <[Solid; 1]>::try_from(fused).unwrap();
    // The outer shell, two cavities and their voids' shells.
    assert_eq!(fused.topology().shells().len(), 6);
    let v = fused.mass_properties().volume;
    assert!((v - 996.0).abs() <= 1e-9 * v, "{v}");
    let b = ball([10.0, 5.0, 2.0], 1.5, 94);
    let (f, _) = fused.fuse(OperationId(95), &b).unwrap();
    let (c, _) = fused.cut(OperationId(96), &b).unwrap();
    let (m, _) = fused.common(OperationId(97), &b).unwrap();
    let (vf, vc, vm) = (volume(&f), volume(&c), volume(&m));
    let (va, vb) = (fused.mass_properties().volume, b.mass_properties().volume);
    assert!((vf + vm - va - vb).abs() <= 1e-9 * va, "{vf} {vm}");
    assert!((vc + vm - va).abs() <= 1e-9 * va, "{vc} {vm}");
    for s in f.iter().chain(&c) {
        assert_eq!(s.topology().shells().len(), 6);
    }
}
