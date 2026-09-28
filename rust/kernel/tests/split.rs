//! S8: prisms split by a plane, against the independent reference
//! (`fixtures/split-*.txt|tsv` from `tools/generate_split_fixtures.py`).
#[path = "support/split_protocol.rs"]
mod protocol;
use protocol::{cases, plane_frame, rows, split, Case};
use rusty_occt::history::{self, Relation};
use rusty_occt::identity::OperationId;
use rusty_occt::{Frame3, Point3, Side, Tolerance, Vec3};
use std::collections::BTreeMap;

fn all() -> Vec<Case> {
    cases(include_str!("../../fixtures/split-cases.txt"))
}

fn expected() -> BTreeMap<String, Vec<(String, [f64; 5])>> {
    let mut out: BTreeMap<String, Vec<(String, [f64; 5])>> = BTreeMap::new();
    for line in include_str!("../../fixtures/split-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let v: Vec<f64> = w[2..7].iter().map(|x| x.parse().unwrap()).collect();
        out.entry(name.to_string())
            .or_default()
            .push((w[1].to_string(), [v[0], v[1], v[2], v[3], v[4]]));
    }
    out
}

fn frame_of(case: &Case) -> Frame3 {
    let f = case.spec.frame;
    Frame3::new(
        Point3::new(f[0], f[1], f[2]),
        Vec3::new(f[3], f[4], f[5]),
        Vec3::new(f[6], f[7], f[8]),
        case.spec.tolerance,
    )
    .unwrap()
}

/// Whether the plane is oblique to the prism's axis (S8a.2), in binary64:
/// exact components along the stored axes are either zero or not far from.
fn oblique(case: &Case) -> bool {
    let f = frame_of(case);
    let m = Vec3::new(case.plane[3], case.plane[4], case.plane[5]);
    let (a, b, c) = (m.dot(f.x()), m.dot(f.y()), m.dot(f.normal()));
    c != 0.0 && (a != 0.0 || b != 0.0)
}

/// The kernel stores the frames the reference used, bit for bit.
#[test]
fn stored_frames_are_the_reference_inputs() {
    let all: BTreeMap<String, Case> = all()
        .into_iter()
        .map(|c| (c.spec.name.clone(), c))
        .collect();
    for row in include_str!("../../fixtures/split-frames.tsv")
        .lines()
        .skip(1)
    {
        let w: Vec<&str> = row.split('\t').collect();
        let f = frame_of(&all[w[0]]);
        let got = match w[1] {
            "n" => f.normal(),
            "x" => f.x(),
            _ => f.y(),
        };
        let want: Vec<u64> = w[2]
            .split(' ')
            .map(|h| u64::from_str_radix(h, 16).unwrap())
            .collect();
        assert_eq!(got.to_array().map(f64::to_bits).to_vec(), want, "{}", w[0]);
    }
}

/// Each side's volumes, areas and centre (summed over its pieces) contain
/// the reference's; only planes oblique to the axis wait for S8a.2.
#[test]
fn every_case_matches_the_reference() {
    let want = expected();
    let mut failures = Vec::new();
    for case in all() {
        let name = case.spec.name.clone();
        let got = rows(&case).unwrap_or_else(|e| panic!("{name}: {e}"));
        if got == ["unsupported"] {
            if !oblique(&case) {
                failures.push(format!("{name}: unsupported"));
            }
            continue;
        }
        let mut sides: BTreeMap<String, [[f64; 2]; 5]> = BTreeMap::new();
        for row in &got {
            let w: Vec<&str> = row.split(' ').collect();
            let v: Vec<f64> = w[2..12].iter().map(|x| x.parse().unwrap()).collect();
            let s = sides.entry(w[1].to_string()).or_insert([[0.0; 2]; 5]);
            s[0] = [s[0][0] + v[0], s[0][1] + v[1]];
            s[1] = [s[1][0] + v[2], s[1][1] + v[3]];
            for i in 0..3 {
                let p = [
                    v[0] * v[4 + 2 * i],
                    v[0] * v[5 + 2 * i],
                    v[1] * v[4 + 2 * i],
                    v[1] * v[5 + 2 * i],
                ];
                let (lo, hi) = (
                    p.iter().copied().fold(f64::MAX, f64::min),
                    p.iter().copied().fold(f64::MIN, f64::max),
                );
                s[2 + i] = [s[2 + i][0] + lo, s[2 + i][1] + hi];
            }
        }
        let rows = &want[&name];
        if rows[0].0 == "whole" {
            if got.len() != 1 {
                failures.push(format!("{name}: split where whole"));
            }
            let s: Vec<[[f64; 2]; 5]> = sides.into_values().collect();
            sides = BTreeMap::from([("whole".to_string(), s[0])]);
        }
        for (side, v) in rows {
            let Some(s) = sides.get(side) else {
                failures.push(format!("{name}: no {side}"));
                continue;
            };
            let inside =
                |x: f64, [lo, hi]: [f64; 2]| lo - 1e-20 * x.abs() <= x && x <= hi + 1e-20 * x.abs();
            let moment = |i: usize| v[0] * v[2 + i];
            let near = |x: f64, [lo, hi]: [f64; 2]| {
                lo - 1e-12 * x.abs().max(1.0) <= x && x <= hi + 1e-12 * x.abs().max(1.0)
            };
            if !inside(v[0], s[0])
                || !inside(v[1], s[1])
                || (0..3).any(|i| !near(moment(i), s[2 + i]))
            {
                failures.push(format!("{name} {side}: {s:?} misses {v:?}"));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// Every split's history passes the independent history check, covers
/// every input entity, and is deterministic.
#[test]
fn histories_are_complete_and_deterministic() {
    for case in all() {
        let name = case.spec.name.clone();
        let Ok((solid, pieces, h)) = split(&case) else {
            continue;
        };
        let sets_in = [solid.topology().entity_set(solid.resolution())];
        let sets_out: Vec<_> = pieces
            .iter()
            .map(|(_, p)| p.topology().entity_set(p.resolution()))
            .collect();
        let issues = history::check(&sets_in, &sets_out, &h);
        assert!(issues.is_empty(), "{name}: {issues:?}");
        let covered: std::collections::BTreeSet<_> =
            h.relations.iter().flat_map(Relation::sources).collect();
        for (id, _) in solid.topology().ids() {
            assert!(covered.contains(&id), "{name}: {id:?} has no relation");
        }
        let (_, again, h2) = split(&case).unwrap();
        assert_eq!(h, h2, "{name}");
        for ((s1, p1), (s2, p2)) in pieces.iter().zip(&again) {
            assert_eq!(s1, s2);
            let ids =
                |p: &rusty_occt::Solid| p.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
            assert_eq!(ids(p1), ids(p2), "{name}");
        }
    }
}

/// Vertical planes through a square prism at many offsets: two prisms whose
/// volumes add up, the smaller side first by side.
#[test]
fn vertical_cuts_conserve_volume() {
    let square = rusty_occt::Boundary::rectangle(10.0, 10.0, Tolerance::default()).unwrap();
    let profile = rusty_occt::Profile::new(square, vec![], Tolerance::default()).unwrap();
    let (solid, _) =
        rusty_occt::Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 4.0).unwrap();
    for k in 1..20 {
        let x = -5.0 + 10.0 * k as f64 / 20.0;
        let plane = plane_frame([x, 0.0, 0.0, 1.0, 0.25 * (k % 3) as f64, 0.0]);
        let (pieces, _) = solid.split_by_plane(OperationId(2), plane).unwrap();
        let total: f64 = pieces.iter().map(|(_, p)| p.mass_properties().volume).sum();
        assert!((total - 400.0).abs() < 1e-9, "{k}: {total}");
        assert!(pieces.iter().any(|(s, _)| *s == Side::Below));
    }
}

/// A round hole the plane misses stays whole in its piece (a fuzz find: the
/// hole's whole circle has no stored vertex).
#[test]
fn a_round_hole_off_the_plane_stays_whole() {
    let t = Tolerance::default();
    let outer = rusty_occt::Boundary::rectangle(10.0, 10.0, t).unwrap();
    let hole = rusty_occt::Boundary::circle(rusty_occt::Point2::new(5.0, 5.0), 2.0, t).unwrap();
    let profile = rusty_occt::Profile::new(outer, vec![hole], t).unwrap();
    let (solid, _) =
        rusty_occt::Solid::extrude_with(OperationId(1), profile, Frame3::xy(), 0.0, 1.0).unwrap();
    let (pieces, _) = solid
        .split_by_plane(OperationId(2), plane_frame([8.0, 0.0, 0.0, 1.0, 0.0, 0.0]))
        .unwrap();
    assert_eq!(pieces.len(), 2);
    let volumes: Vec<f64> = pieces
        .iter()
        .map(|(_, p)| p.mass_properties().volume)
        .collect();
    let hole_area = std::f64::consts::PI * 4.0;
    assert!(
        (volumes[0] - (80.0 - hole_area)).abs() < 1e-9,
        "{volumes:?}"
    );
    assert!((volumes[1] - 20.0).abs() < 1e-9, "{volumes:?}");
}
