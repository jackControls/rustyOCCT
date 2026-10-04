//! S9f.3b: Booleans of spline prisms against cones and frustums in any
//! relative position, against the independent reference
//! (`fixtures/boolean-spline-cone-*` from
//! `tools/generate_spline_cone_boolean_fixtures.py`): rulings meeting both
//! nappes (`A < 0`, the apex inside the prism and outside it), loops turning
//! back inside the faces on graphs over the height and branches over the
//! spline's parameter (`A > 0`), a rim's points on a wall's generatrices in
//! a tower field, as the reference; the apex on a wall, a wall along the
//! cone's ruling and a tangency `Degenerate` with the decisions' reasons;
//! cones in turned frames by inclusion and exclusion.
use rusty_occt::identity::OperationId;
use rusty_occt::topology::SplineSpan;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance,
    Vec3,
};

#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

const CASES: &str = include_str!("../../fixtures/boolean-spline-cone-cases.txt");
const EXPECTED: &str = include_str!("../../fixtures/boolean-spline-cone-expected.tsv");

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Each case's kind, its refusal's reason and, for a result, its solids,
/// volume, area and centre.
type Expect = (String, String, Option<(usize, [f64; 5])>);

fn expected(text: &str) -> std::collections::BTreeMap<String, Expect> {
    let mut expect = std::collections::BTreeMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let entry = expect
            .entry(name.to_string())
            .or_insert((String::new(), String::new(), None));
        match w[0] {
            "expect" => entry.0 = w[1].to_string(),
            "reason" => entry.1 = row["reason ".len()..].to_string(),
            "result" => {
                let v: Vec<f64> = w[2..7].iter().map(|x| x.parse().unwrap()).collect();
                entry.2 = Some((
                    w[1].parse::<usize>().unwrap(),
                    [v[0], v[1], v[2], v[3], v[4]],
                ));
            }
            _ => {}
        }
    }
    expect
}

/// Every case as the reference declares it: degenerate ones refused with
/// the decisions' reason, the others with the reference's solids and their
/// volume, area and centre within enclosures narrower than 1e-9 relative
/// (a wide one would hold the reference yet report another).
#[test]
fn every_case_matches_the_reference() {
    let expect = expected(EXPECTED);
    let mut failures = Vec::new();
    for case in protocol::cases(CASES) {
        let (kind, reason, want) = &expect[&case.name];
        if kind == "degenerate" {
            match protocol::run(&case) {
                Err(Error::Degenerate(m)) if m == reason => {}
                other => failures.push(format!(
                    "{}: {:?} not refused as {reason}",
                    case.name,
                    other.map(|r| r.2.len())
                )),
            }
            continue;
        }
        let rows = match protocol::rows(&case) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{}: {e}", case.name));
                continue;
            }
        };
        if kind == "empty" {
            if rows != ["empty"] {
                failures.push(format!("{}: {rows:?} not empty", case.name));
            }
            continue;
        }
        let (count, v) = want.expect("a result");
        let solids: Vec<Vec<f64>> = rows
            .iter()
            .map(|r| {
                r.split(' ')
                    .skip(1)
                    .take(10)
                    .map(|x| x.parse().unwrap_or(f64::NAN))
                    .collect()
            })
            .collect();
        if solids.len() != count || solids.iter().any(|s| s.len() < 10) {
            failures.push(format!("{}: {rows:?} for {count} solids", case.name));
            continue;
        }
        let sum = |i: usize| solids.iter().map(|s| s[i]).sum::<f64>();
        let near = |x: f64, lo: f64, hi: f64| {
            let slack = 1e-9 * x.abs().max(1.0);
            lo - slack <= x && x <= hi + slack
        };
        let narrow = |lo: f64, hi: f64| hi - lo <= 1e-9 * lo.abs().max(hi.abs()).max(1.0);
        if !narrow(sum(0), sum(1)) || !narrow(sum(2), sum(3)) {
            failures.push(format!(
                "{}: wide enclosures [{}, {}], [{}, {}]",
                case.name,
                sum(0),
                sum(1),
                sum(2),
                sum(3)
            ));
        }
        if !near(v[0], sum(0), sum(1)) || !near(v[1], sum(2), sum(3)) {
            failures.push(format!(
                "{}: volume {} area {} against [{}, {}], [{}, {}]",
                case.name,
                v[0],
                v[1],
                sum(0),
                sum(1),
                sum(2),
                sum(3)
            ));
        }
        for i in 0..3 {
            let moment = |pick: fn(f64, f64) -> f64| {
                solids
                    .iter()
                    .map(|s| {
                        let (vl, vh, cl, ch) = (s[0], s[1], s[4 + 2 * i], s[5 + 2 * i]);
                        pick(pick(vl * cl, vl * ch), pick(vh * cl, vh * ch))
                    })
                    .sum::<f64>()
            };
            let (lo, hi) = (moment(f64::min), moment(f64::max));
            let m = v[0] * v[2 + i];
            let slack = 1e-9 * m.abs().max(v[0]).max(1.0);
            if !(lo - slack <= m && m <= hi + slack) || hi - lo > slack {
                failures.push(format!(
                    "{}: centre {i} {} against [{}, {}]",
                    case.name,
                    v[2 + i],
                    lo / v[0],
                    hi / v[0]
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures: {failures:#?}",
        failures.len()
    );
}

fn run(name: &str) -> protocol::Run {
    let case = protocol::cases(CASES)
        .into_iter()
        .find(|c| c.name == name)
        .unwrap();
    protocol::run(&case).unwrap_or_else(|e| panic!("{name}: {e}"))
}

fn meetings(out: &[Solid]) -> Vec<rusty_occt::topology::WallMeet> {
    use rusty_occt::topology::Curve3;
    out.iter()
        .flat_map(|s| s.topology().edges().to_vec())
        .filter_map(|e| match e.curve {
            Curve3::WallMeet(m) => Some(*m),
            _ => None,
        })
        .collect()
}

/// The other surface's (a cone's) function at a point, `|(w . x2, w . y2)|
/// - (R + (w . n2) tan a2)`.
fn cone_gap(m: &rusty_occt::topology::WallMeet, p: Point3) -> f64 {
    let w = p - m.other.origin();
    let r = m.other_radius + w.dot(m.other.normal()) * m.other_half_angle.tan();
    w.dot(m.other.x()).hypot(w.dot(m.other.y())) - r
}

/// A thin frustum across the dome's axis piercing its arch (`A > 0`) meets
/// it in a loop: about each of its two turning points a graph over the
/// wall's `v` (a window inside one knot span), graphs over `u` between
/// them; a frustum on the bulge's axis (`A < 0`) in one branch over `u` on
/// its own nappe, no window. Every piece's points lie on the wall and on the
/// cone (its stored half angle's), on the cone's face.
#[test]
fn loops_and_branches_on_the_cone() {
    for (name, windows) in [("dome_pierce_common", 2), ("bulge_frustum_common", 0)] {
        let (_, _, out, _) = run(name);
        let meets = meetings(&out);
        assert!(!meets.is_empty(), "{name}");
        assert!(
            meets
                .iter()
                .all(|m| m.other_half_angle != 0.0 && !m.other_sphere),
            "{name}"
        );
        let over_v: Vec<_> = meets.iter().filter_map(|m| m.window).collect();
        assert_eq!(over_v.len(), windows, "{name}: {meets:?}");
        for m in &meets {
            for k in 0..=8 {
                let f = f64::from(k) / 8.0;
                let p = m.point(f);
                let (u, v) = m.parameters(f);
                assert!(cone_gap(m, p).abs() < 1e-12, "{name}: {p:?}");
                // (At a cap the height may round past the wall's domain.)
                let ((u0, u1), (v0, v1)) = m.wall.domain();
                let (uc, vc) = (u.clamp(u0, u1), v.clamp(v0, v1));
                assert!((v - vc).abs() < 1e-12, "{name}: {v}");
                assert!(
                    (m.wall.point(uc, vc).unwrap() - p).length() < 1e-12,
                    "{name}"
                );
                // On the cone's own nappe: its radius term positive.
                let w = p - m.other.origin();
                let r = m.other_radius + w.dot(m.other.normal()) * m.other_half_angle.tan();
                assert!(r > 0.0, "{name}: {p:?}");
                if let Some([a, b]) = m.window {
                    assert!(a < u && u < b, "{name}: {u} in {a} {b}");
                    let knots = m.wall.u_knots().knots();
                    assert!(!knots.iter().any(|k| a < *k && *k < b), "{name}");
                }
            }
        }
    }
}

/// A frustum across the dome's axis whose bottom rim's plane `y = 5/4`
/// holds the dome's axis: the rim meets the arch's generatrices at `x = 2
/// +- sqrt(3/2)` in a tower field: the common's vertices there on the rim
/// and on the dome, inside the prism's heights.
#[test]
fn a_rims_tower_points_are_vertices() {
    let (_, b, out, _) = run("dome_side_cone_common");
    let c = b.frame().origin();
    let on_rim: Vec<Point3> = out
        .iter()
        .flat_map(|s| s.topology().vertices().to_vec())
        .map(|v| v.position)
        .filter(|p| {
            (p.y - 1.25).abs() < 1e-12 && ((p.x - c.x).hypot(p.z - c.z) - 1.5).abs() < 1e-12
        })
        .collect();
    assert!(on_rim.len() >= 2, "{on_rim:?}");
    for p in &on_rim {
        assert!(((p.x - 2.0).abs() - 1.5f64.sqrt()).abs() < 1e-12, "{p:?}");
        assert!(p.z > 0.0 && p.z < 3.0, "{p:?}");
    }
}

#[test]
#[ignore]
fn debug_errors() {
    let want = std::env::var("CASE").unwrap_or_default();
    for case in protocol::cases(CASES) {
        if !case.name.contains(&want) {
            continue;
        }
        let t = std::time::Instant::now();
        match protocol::run(&case) {
            Ok((_, _, out, _)) => println!(
                "{} ({:?}): {} solids {:?}",
                case.name,
                t.elapsed(),
                out.len(),
                out.iter()
                    .map(|s| s.mass_properties().volume)
                    .collect::<Vec<_>>()
            ),
            Err(e) => println!("{} ({:?}): {e:?}", case.name, t.elapsed()),
        }
    }
}

#[test]
fn fixture_histories_are_complete() {
    for case in protocol::cases(CASES) {
        let Ok((a, b, out, h)) = protocol::run(&case) else {
            continue;
        };
        let ins = [
            a.topology().entity_set(a.resolution()),
            b.topology().entity_set(b.resolution()),
        ];
        let outs: Vec<_> = out
            .iter()
            .map(|s| s.topology().entity_set(s.resolution()))
            .collect();
        let issues = history::check(&ins, &outs, &h);
        assert!(issues.is_empty(), "{}: {issues:?}", case.name);
    }
}

#[test]
fn results_are_deterministic_and_move_rigidly() {
    use rusty_occt::{Location, RigidTransform};
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    for case in protocol::cases(CASES) {
        let Ok((_, _, out, h)) = protocol::run(&case) else {
            continue;
        };
        let (_, _, again, h2) = protocol::run(&case).unwrap();
        assert_eq!(h, h2, "{}", case.name);
        for (x, y) in out.iter().zip(&again) {
            assert_eq!(ids(x), ids(y), "{}", case.name);
            assert_eq!(
                format!("{:?}", x.topology().vertices()),
                format!("{:?}", y.topology().vertices()),
                "{}",
                case.name
            );
        }
        for s in &out {
            let (moved, _) = s
                .transform_with(OperationId(900), motion)
                .unwrap_or_else(|e| panic!("{}: {e}", case.name));
            assert_eq!(ids(&moved), ids(s), "{}", case.name);
            let (v0, v1) = (s.mass_properties().volume, moved.mass_properties().volume);
            assert!(
                (v0 - v1).abs() <= 1e-9 * v0.abs().max(1.0),
                "{}: {v0} {v1}",
                case.name
            );
            for v in moved.topology().vertices() {
                assert_eq!(
                    moved.classify(v.position).unwrap(),
                    Location::Boundary,
                    "{}",
                    case.name
                );
            }
        }
    }
}

// ------------------------------------------------------------ beyond the fixtures

fn spline_segment(poles: &[(f64, f64)], knots: Vec<f64>, mults: Vec<usize>) -> Segment {
    let degree = mults[0] - 1;
    let poles = poles.iter().map(|(x, y)| Point2::new(*x, *y)).collect();
    let curve = BSplineCurve2::new(degree, poles, None, knots, mults).expect("a spline");
    Segment::Spline(SplineSpan::whole(curve))
}

/// The dome `y = x (4 - x) / 2` over `[0, 4]`, closed by its base.
fn dome() -> Profile {
    let outer = Boundary::path(
        vec![Point2::new(0.0, 0.0), Point2::new(4.0, 0.0)],
        vec![
            Segment::Line,
            spline_segment(
                &[(4.0, 0.0), (2.0, 4.0), (0.0, 0.0)],
                vec![0.0, 1.0],
                vec![3, 3],
            ),
        ],
        tol(),
    )
    .expect("the dome");
    Profile::new(outer, vec![], tol()).expect("a profile")
}

fn frame(origin: (f64, f64, f64), normal: (f64, f64, f64), x: (f64, f64, f64)) -> Frame3 {
    Frame3::new(
        Point3::new(origin.0, origin.1, origin.2),
        Vec3::new(normal.0, normal.1, normal.2),
        Vec3::new(x.0, x.1, x.2),
        tol(),
    )
    .expect("a frame")
}

fn prism(id: u64, profile: Profile, f: Frame3, lo: f64, hi: f64) -> Solid {
    Solid::extrude_with(OperationId(id), profile, f, lo, hi)
        .expect("a prism")
        .0
}

fn agree(a: &Solid, b: &Solid) -> f64 {
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let volume = |r: Vec<Solid>| r.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let fused = volume(a.fuse(OperationId(5), b).expect("a fuse").0);
    let common = volume(a.common(OperationId(6), b).expect("a common").0);
    let cut = volume(a.cut(OperationId(7), b).expect("a cut").0);
    assert!(
        (fused - (va + vb - common)).abs() <= 1e-9 * fused,
        "{fused} {va} {vb} {common}"
    );
    assert!(
        (cut - (va - common)).abs() <= 1e-9 * va,
        "{cut} {va} {common}"
    );
    common
}

/// Cones in turned frames (stored axes not exactly orthonormal, the half
/// angle's tangent not the model's slope): a frustum leaning across the
/// dome's arch (`A > 0`) and one in the tilted frame along the dome's axis
/// within its half angle (`A < 0`), each result validated as it is built,
/// the three operations by inclusion and exclusion.
#[test]
fn cones_in_turned_frames_agree() {
    let a = prism(1, dome(), Frame3::xy(), 0.0, 3.0);
    let lean = frame((2.5, 0.5, 1.5), (0.0, 4.0, 3.0), (1.0, 0.0, 0.0));
    let (k, _) = Solid::cone_with(OperationId(2), lean, 0.75, 0.5, 3.0, tol()).expect("a frustum");
    assert!(agree(&a, &k) > 0.0);
    let tilt = frame((2.0, 1.0, -1.0), (0.0, 0.25, 1.0), (1.0, 0.0, 0.0));
    let (k, _) = Solid::cone_with(OperationId(3), tilt, 2.0, 0.5, 5.0, tol()).expect("a frustum");
    assert!(agree(&a, &k) > 0.0);
}
