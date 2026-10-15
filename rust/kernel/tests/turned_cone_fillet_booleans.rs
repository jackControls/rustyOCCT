//! S9e.4b.4d.2c.1b: imported turned bodies whose fillets and rounds (tori)
//! join a cone smoothly (`OCC496b`'s third torus's kind: a round between a
//! cylinder and a cone, a fillet under a conical boss on a base, a widening
//! cone rounded into a top disc, the first's round refitted by 2^-20 at a
//! resolution of 1e-3 with a cone of no Pythagorean direction), each torus
//! snapped tangent to its lines, given to Booleans whose partners cross no
//! smooth edge, against the independent reference
//! (`fixtures/boolean-turned-cone-fillets-*` from
//! `tools/generate_turned_cone_fillets_boolean_fixtures.py`, the bodies under
//! `fixtures/imported/cone_fillet_*.brep`); a ball touching a round refused
//! (a tangency between the inputs). Each case runs once (on a few threads)
//! for the checks that read its result.
#[path = "support/boolean_protocol.rs"]
#[allow(dead_code)]
mod protocol;
use rusty_occt::history::{self, Resolution};
use rusty_occt::identity::OperationId;
use rusty_occt::{Error, Location, Point3, RigidTransform, Solid, Vec3};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

fn cases() -> Vec<protocol::Case> {
    protocol::cases(include_str!(
        "../../fixtures/boolean-turned-cone-fillets-cases.txt"
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
    for line in include_str!("../../fixtures/boolean-turned-cone-fillets-expected.tsv")
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

/// A list of solids' summed volume enclosure.
fn volume_of(out: &[Solid]) -> [f64; 2] {
    out.iter().fold([0.0, 0.0], |acc, s| {
        let v = s.topology().mass_enclosure().expect("a mass").volume;
        [acc[0] + v[0], acc[1] + v[1]]
    })
}

/// Whether two sums of enclosures meet (their widths and a rounding of the
/// sums apart at most).
fn agree(x: [f64; 2], y: [f64; 2]) -> bool {
    let slack = 1e-12 * x[1].abs().max(y[1].abs()).max(1.0);
    x[0] <= y[1] + slack && y[0] <= x[1] + slack
}

fn add(x: [f64; 2], y: [f64; 2]) -> [f64; 2] {
    [x[0] + y[0], x[1] + y[1]]
}

/// Each group's results consistent with one another and with their inputs
/// within their certified enclosures: `A + B` the fuse's and the common's,
/// `A` the cut's and the common's, `B` the reversed cut's and the common's,
/// the fuse in either order the same; inputs symmetric under a mirror give
/// results whose centroids lie on its plane (the rod and the ball on the
/// planes `x = 0` and `y = 0`, each body about the world's `z`).
#[test]
fn results_are_consistent() {
    let first: BTreeMap<&str, _> = runs().iter().map(|(n, r)| (n.as_str(), r)).collect();
    let mut failures = Vec::new();
    for group in ["round_rod", "top_box", "refit_box", "fillet_ball"] {
        let get = |op: &str| match first[format!("{group}_{op}").as_str()] {
            Ok(run) => run,
            Err(e) => panic!("{group}_{op}: {e}"),
        };
        let (a, b, fuse, _) = get("fuse");
        let (_, _, cut, _) = get("cut");
        let (_, _, common, _) = get("common");
        let va = a.topology().mass_enclosure().unwrap().volume;
        let vb = b.topology().mass_enclosure().unwrap().volume;
        let (vf, vc, vn) = (volume_of(fuse), volume_of(cut), volume_of(common));
        if !agree(add(vf, vn), add(va, vb)) {
            failures.push(format!(
                "{group}: fuse and common {vf:?} {vn:?} for {va:?} {vb:?}"
            ));
        }
        if !agree(add(vc, vn), va) {
            failures.push(format!("{group}: cut and common {vc:?} {vn:?} for {va:?}"));
        }
        let (back, _) = b.cut(OperationId(97), a).unwrap();
        if !agree(add(volume_of(&back), vn), vb) {
            failures.push(format!("{group}: reversed cut and common for {vb:?}"));
        }
        let (swapped, _) = b.fuse(OperationId(97), a).unwrap();
        if swapped.len() != fuse.len() || !agree(volume_of(&swapped), vf) {
            failures.push(format!("{group}: the fuse in either order"));
        }
        let mirror = match group {
            "round_rod" => Some(0),
            "fillet_ball" => Some(1),
            _ => None,
        };
        if let Some(k) = mirror {
            for (op, out) in [("fuse", fuse), ("cut", cut), ("common", common)] {
                for s in out.iter() {
                    let c = s.topology().mass_enclosure().unwrap().centroid[k];
                    if !(c[0] <= 1e-12 && -1e-12 <= c[1]) {
                        failures.push(format!("{group}_{op}: centroid {k} {c:?} off the mirror"));
                    }
                }
            }
        }
    }
    // Each body alone: its centroid on its axis.
    for name in [
        "cone_fillet_round",
        "cone_fillet_fillet",
        "cone_fillet_top",
        "cone_fillet_refit",
    ] {
        let s = body(name, 91).unwrap();
        let c = s.topology().mass_enclosure().unwrap().centroid;
        for (k, ck) in c.iter().enumerate().take(2) {
            if !(ck[0] <= 1e-12 && -1e-12 <= ck[1]) {
                failures.push(format!("{name}: centroid {k} {ck:?} off the axis"));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
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
/// their frames keep the reference's volumes, the imported bodies read again
/// off their moved stored topologies (their tori snapped tangent to their
/// lines again, their cones' lines exact, every primitive on the turned
/// axis frame, the touching circles exact on it); a partner's face parallel
/// or normal to the axis in the construction lies within rounding of it once
/// turned (the rod's and the `TILT` box's), refused for the turn: the ball's
/// and the skew box's cases evaluate turned.
#[test]
fn moved_inputs_keep_their_volumes() {
    let shift = RigidTransform::translation(Vec3::new(0.5, -0.25, 1.0)).unwrap();
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    let expect = expected();
    let all = cases();
    let mut turned = 0;
    for name in [
        "round_rod_cut",
        "fillet_ball_common",
        "top_box_fuse",
        "refit_box_common",
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
                    turned += k;
                }
            }
        }
    }
    assert!(turned >= 2, "{turned} turned cases evaluated");
}

fn body(name: &str, op: u64) -> Result<Solid, Error> {
    let (topology, resolution) = protocol::imported_topology(&format!("imported/{name}.brep"));
    Solid::imported_with(OperationId(op), topology, resolution).map(|(s, _)| s)
}

/// A model profile's segment to the next point: a line, or an arc about
/// `(cr, cw)` of radius `t`, counter-clockwise in `(r, w)` where `ccw`.
#[derive(Clone, Copy)]
enum Seg {
    Line,
    Arc(f64, f64, f64, bool),
}

/// A profile's volume by Pappus: `2 pi` times its `int r^2 / 2 dw` along its
/// boundary (counter-clockwise), each arc's in closed form.
fn pappus(points: &[(f64, f64, Seg)]) -> f64 {
    let mut sum = 0.0;
    for p in points.windows(2) {
        let ((a, w0, seg), (b, w1, _)) = (p[0], p[1]);
        sum += match seg {
            Seg::Line => (w1 - w0) * (a * a + a * b + b * b) / 6.0,
            Seg::Arc(cr, cw, t, ccw) => {
                let th0 = (w0 - cw).atan2(a - cr);
                let mut th1 = (w1 - cw).atan2(b - cr);
                let tau = std::f64::consts::TAU;
                if ccw && th1 < th0 {
                    th1 += tau;
                }
                if !ccw && th1 > th0 {
                    th1 -= tau;
                }
                let prim = |th: f64| {
                    let (s, c) = th.sin_cos();
                    t * (cr * cr * s / 2.0
                        + cr * t * (th / 2.0 + s * c / 2.0)
                        + t * t * (s - s * s * s / 3.0) / 2.0)
                };
                prim(th1) - prim(th0)
            }
        };
    }
    2.0 * std::f64::consts::PI * sum
}

const D: f64 = 1.0 / 1_048_576.0;

/// A model profile: its points and their segments to the next.
type Profile = Vec<(f64, f64, Seg)>;

/// The bodies' model profiles `(r, w, segment)` from the axis out and back
/// to it (about the world's `z`), as the generator's: each round's or
/// fillet's point on its cone in the direction `((1 - s^2), 2 s) / (1 +
/// s^2)` from its centre, the cone that direction's tangent there.
fn profiles() -> Vec<(&'static str, Profile)> {
    use Seg::{Arc, Line};
    // The refit's direction's half-angle tangent on the multiples of 2^-40
    // (the generator's and the kernel's snap of OCCT's), its point on the
    // cone and the cone's radius at the top.
    let s = 32_184_485_779.0 / 137_438_953_472.0;
    let e = [(1.0 - s * s) / (1.0 + s * s), 2.0 * s / (1.0 + s * s)];
    let t1 = [4.0 + 2.0 * e[0], 4.0 + D + 2.0 * e[1]];
    let top = t1[0] - e[1] / e[0] * (8.0 - t1[1]);
    vec![
        (
            "cone_fillet_round",
            vec![
                (0.0, 0.0, Line),
                (6.0, 0.0, Line),
                (6.0, 4.0, Arc(4.0, 4.0, 2.0, true)),
                (98.0 / 17.0, 84.0 / 17.0, Line),
                (1054.0 / 255.0, 8.0, Line),
                (0.0, 8.0, Line),
            ],
        ),
        (
            "cone_fillet_fillet",
            vec![
                (0.0, 0.0, Line),
                (8.0, 0.0, Line),
                (8.0, 2.0, Line),
                (5.0, 2.0, Arc(5.0, 3.0, 1.0, false)),
                (70.0 / 17.0, 43.0 / 17.0, Line),
                (26.0 / 15.0, 7.0, Line),
                (0.0, 7.0, Line),
            ],
        ),
        (
            "cone_fillet_top",
            vec![
                (0.0, 0.0, Line),
                (3.0, 0.0, Line),
                (8.6, 4.2, Arc(8.0, 5.0, 1.0, true)),
                (8.0, 6.0, Line),
                (0.0, 6.0, Line),
            ],
        ),
        (
            "cone_fillet_refit",
            vec![
                (0.0, 0.0, Line),
                (6.0, 0.0, Line),
                (6.0, 4.0 + D, Arc(4.0, 4.0 + D, 2.0, true)),
                (t1[0], t1[1], Line),
                (top, 8.0, Line),
                (0.0, 8.0, Line),
            ],
        ),
    ]
}

/// The profile's radius at a height inside its range (its outer point).
fn radius_at(points: &[(f64, f64, Seg)], w: f64) -> f64 {
    let mut best = 0.0f64;
    for p in points.windows(2) {
        let ((a, w0, seg), (b, w1, _)) = (p[0], p[1]);
        if w1 == w0 || w < w0.min(w1) || w > w0.max(w1) {
            continue;
        }
        best = best.max(match seg {
            Seg::Line => a + (b - a) * (w - w0) / (w1 - w0),
            Seg::Arc(cr, cw, t, _) => {
                let d = (t * t - (w - cw) * (w - cw)).max(0.0).sqrt();
                if a.max(b) > cr {
                    cr + d
                } else {
                    cr - d
                }
            }
        });
    }
    best
}

/// A box about the world's axes holding every body.
fn holder(op: u64) -> Solid {
    use rusty_occt::{Boundary, Frame3, Point2, Profile, Tolerance};
    let tol = Tolerance::default();
    let square = [(-16.0, -16.0), (16.0, -16.0), (16.0, 16.0), (-16.0, 16.0)]
        .map(|(x, y)| Point2::new(x, y))
        .to_vec();
    let frame = Frame3::new(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    let profile = Profile::new(Boundary::polygon(square, tol).unwrap(), Vec::new(), tol).unwrap();
    Solid::extrude_with(OperationId(op), profile, frame, -16.0, 16.0)
        .unwrap()
        .0
}

/// Every body imports as its sharp profile's lattice with its fillets'
/// pieces and less its rounds' (each a slant box, frusta of cones normal to
/// the cone, less its torus snapped tangent to its cone and its plane or
/// cylinder): its volume and its common with a box holding it its model
/// profile's by Pappus (the refit's stored faces' within the cone's turn),
/// points within the profile inside and beyond it outside at heights across
/// its torus, its stored vertices on its boundary; moved by a turn, alike
/// (its axis frame turned, every primitive on it).
#[test]
fn imported_bodies_are_their_model_profiles() {
    let turn = RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
        .unwrap();
    for (name, model) in profiles() {
        let s = body(name, 91).unwrap_or_else(|e| panic!("{name}: {e}"));
        let want = pappus(&model);
        let v = s.mass_properties().volume;
        assert!((v - want).abs() <= 1e-9 * want, "{name}: {v} for {want}");
        let (out, _) = s.common(OperationId(95), &holder(94)).unwrap();
        let got: f64 = out.iter().map(|x| x.mass_properties().volume).sum();
        assert_eq!(out.len(), 1, "{name}");
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "{name}: {got} for {want}"
        );
        let top = model.iter().map(|p| p.1).fold(0.0f64, f64::max);
        for k in 1..24 {
            let w = top * (f64::from(k) + 0.37) / 24.0;
            let r = radius_at(&model, w);
            for (j, angle) in [0.3f64, 2.1, 4.4].into_iter().enumerate() {
                let (c, sn) = (angle.cos(), angle.sin());
                let inside = Point3::new(0.95 * r * c, 0.95 * r * sn, w);
                let outside = Point3::new((r + 0.05) * c, (r + 0.05) * sn, w);
                assert_eq!(
                    s.classify(inside).unwrap(),
                    Location::Inside,
                    "{name} {k} {j}"
                );
                assert_eq!(
                    s.classify(outside).unwrap(),
                    Location::Outside,
                    "{name} {k} {j}"
                );
            }
        }
        let (moved, _) = s.transform_with(OperationId(900), turn).unwrap();
        let v = moved.mass_properties().volume;
        assert!((v - want).abs() <= 1e-9 * want, "{name}: moved {v}");
        let (out, _) = moved.common(OperationId(95), &holder(94)).unwrap();
        let got: f64 = out.iter().map(|x| x.mass_properties().volume).sum();
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "{name}: moved common {got} for {want}"
        );
        for x in [&s, &moved] {
            for v in x.topology().vertices() {
                assert_eq!(
                    x.classify(v.position).unwrap(),
                    Location::Boundary,
                    "{name}"
                );
            }
        }
    }
}

/// A torus meeting cones at corners (`turned_corner`) keeps the general
/// refusal: a torus at a corner is no turned body's (`OCC496b`'s first two,
/// S9e.4b.4d.2c.1c's); a round off its tangency (`turned_round_off`) keeps
/// its refusal by name.
#[test]
fn a_torus_at_a_corner_stays_refused() {
    const SMOOTH: &str =
        "an imported body whose torus faces join their neighbours smoothly (S9e.4b.4d)";
    const GENERAL: &str = "an imported solid other than a prism, a sphere, a cone, a torus or a \
                           Boolean of its primitives (S9e.4b)";
    match body("turned_corner", 91) {
        Err(Error::OutOfDomain(m)) => assert_eq!(m, GENERAL),
        other => panic!("turned_corner: {:?}", other.map(|_| ())),
    }
    match body("turned_round_off", 91) {
        Err(Error::OutOfDomain(m)) => assert_eq!(m, SMOOTH),
        other => panic!("turned_round_off: {:?}", other.map(|_| ())),
    }
}
