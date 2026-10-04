//! S9f.3a: Booleans of spline prisms against spheres and hemispheres in any
//! relative position, against the independent reference
//! (`fixtures/boolean-spline-sphere-*` from
//! `tools/generate_spline_sphere_boolean_fixtures.py`): loops turning back
//! inside the faces on graphs over the height, branches over the spline's
//! parameter, hemispheres' split great circles and rims on spline walls in
//! tower fields, as the reference; tangencies, turning points at knots and
//! on faces' boundaries `Degenerate` with the decisions' reasons; zones in
//! turned frames by inclusion and exclusion; and the decisions' other
//! refusal: a spline prism against a cone (S9f.3b).
use rusty_occt::identity::OperationId;
use rusty_occt::topology::SplineSpan;
use rusty_occt::{
    BSplineCurve2, Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance,
    Vec3,
};

#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

const CASES: &str = include_str!("../../fixtures/boolean-spline-sphere-cases.txt");
const EXPECTED: &str = include_str!("../../fixtures/boolean-spline-sphere-expected.tsv");

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

/// A sphere straddling a spline wall meets it in a loop: about each of its
/// two turning points a graph over the wall's `v` on the sphere (a
/// `Curve3::WallMeet` with `other_sphere` and a window inside one knot
/// span), graphs over `u` between them; a large sphere over the dome in one
/// branch over `u` from edge to edge, no window. Every piece's points lie
/// on the wall and on the sphere.
#[test]
fn loops_and_branches_on_the_sphere() {
    for (name, windows) in [("bulge_ball_common", 2), ("dome_ball_common", 0)] {
        let (_, _, out, _) = run(name);
        let meets = meetings(&out);
        assert!(!meets.is_empty(), "{name}");
        assert!(meets.iter().all(|m| m.other_sphere), "{name}");
        let over_v: Vec<_> = meets.iter().filter_map(|m| m.window).collect();
        assert_eq!(over_v.len(), windows, "{name}: {meets:?}");
        for m in &meets {
            let c = m.other.origin();
            for k in 0..=8 {
                let f = f64::from(k) / 8.0;
                let p = m.point(f);
                let (u, v) = m.parameters(f);
                assert!(
                    ((p - c).length() - m.other_radius).abs() < 1e-12,
                    "{name}: {p:?}"
                );
                assert!((m.wall.point(u, v).unwrap() - p).length() < 1e-12, "{name}");
                if let Some([a, b]) = m.window {
                    assert!(a < u && u < b, "{name}: {u} in {a} {b}");
                    let knots = m.wall.u_knots().knots();
                    assert!(!knots.iter().any(|k| a < *k && *k < b), "{name}");
                }
            }
        }
    }
}

/// A hemisphere on its side whose rim's plane `x = 10.6` holds the bulge's
/// axis: the rim meets the wall's generatrices there in a tower field, found
/// in one field (a primitive element): the common's vertices on that plane
/// and on the sphere, inside the wall's heights.
#[test]
fn a_rims_tower_points_are_vertices() {
    let (_, b, out, _) = run("side_hemi_bulge_common");
    let centre = b.frame().origin();
    let on_rim: Vec<Point3> = out
        .iter()
        .flat_map(|s| s.topology().vertices().to_vec())
        .map(|v| v.position)
        .filter(|p| (p.x - 10.6).abs() < 1e-12 && ((*p - centre).length() - 0.9).abs() < 1e-12)
        .collect();
    assert!(on_rim.len() >= 2, "{on_rim:?}");
    for p in &on_rim {
        // On the bulge: x = 10 + (2/3) y (1 - y/6).
        let x = 10.0 + 2.0 / 3.0 * p.y * (1.0 - p.y / 6.0);
        assert!((x - p.x).abs() < 1e-12, "{p:?}");
        assert!(p.z > 0.0 && p.z < 3.0);
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

/// Zones in turned frames (latitudes whose heights round): their rims in
/// planes oblique to the dome's wall meet it along their creases, their end
/// discs cut it in creases, the sphere meets it in a loop; and a cap whose
/// rim's plane holds the wall's axis (on its side). Each result validated as
/// it is built, the three operations by inclusion and exclusion.
#[test]
fn zones_and_caps_in_turned_frames_agree() {
    let a = prism(1, dome(), Frame3::xy(), 0.0, 3.0);
    let tilt = frame((2.0, 1.75, 1.5), (0.0, 3.0, 4.0), (1.0, 0.0, 0.0));
    let (zone, _) =
        Solid::sphere_with(OperationId(2), tilt, 1.0, -0.5, 0.75, tol()).expect("a zone");
    assert!(agree(&a, &zone) > 0.0);
    let side = frame((2.25, 1.75, 1.5), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0));
    let half = std::f64::consts::FRAC_PI_2;
    let (cap, _) =
        Solid::sphere_with(OperationId(3), side, 1.0, -half, 0.25, tol()).expect("a cap");
    assert!(agree(&a, &cap) > 0.0);
}

/// A spline prism against a cone stays refused (S9f.3b's).
#[test]
fn a_cone_against_a_spline_prism_stays_refused() {
    let a = prism(1, dome(), Frame3::xy(), 0.0, 3.0);
    let (cone, _) = Solid::cone_with(
        OperationId(2),
        frame((2.0, 1.5, -1.0), (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
        1.0,
        0.5,
        5.0,
        tol(),
    )
    .expect("a cone");
    let r = a.fuse(OperationId(5), &cone);
    assert!(
        matches!(&r, Err(Error::OutOfDomain(m)) if m.contains("S9f.3b")),
        "{:?}",
        r.map(|x| x.0.len())
    );
}
