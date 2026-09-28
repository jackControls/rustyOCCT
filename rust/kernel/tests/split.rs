//! S8: prisms split by a plane, against the independent reference
//! (`fixtures/split-*.txt|tsv` from `tools/generate_split_fixtures.py`).
#[path = "support/split_protocol.rs"]
mod protocol;
use protocol::{
    build_primitive, cases, plane_frame, primitive_cases, primitive_split, rows, split, Case,
};
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
/// the reference's, for planes normal, parallel and oblique to the axis.
#[test]
fn every_case_matches_the_reference() {
    let want = expected();
    let mut failures = Vec::new();
    for case in all() {
        let name = case.spec.name.clone();
        let got = rows(&case).unwrap_or_else(|e| panic!("{name}: {e}"));
        if got == ["unsupported"] || got == ["limit"] {
            failures.push(format!("{name}: {got:?}"));
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

fn case(name: &str) -> Case {
    all()
        .into_iter()
        .find(|c| c.spec.name == name)
        .unwrap_or_else(|| panic!("no case {name}"))
}

/// S8a.2: a rigid motion of an oblique piece rebuilds it in the new frame
/// with every id, and its volume.
#[test]
fn oblique_pieces_move_rigidly() {
    for name in [
        "sq_tilted",
        "stadium_tilted",
        "disc_wall_ellipse",
        "slot_hole_tilted",
    ] {
        let (_, pieces, _) = split(&case(name)).unwrap();
        let motion = rusty_occt::RigidTransform::rotation(
            Point3::new(1.0, -2.0, 0.5),
            Vec3::new(0.3, -0.4, 0.8),
            0.7,
        )
        .unwrap();
        for (_, piece) in &pieces {
            let (moved, h) = piece.transform_with(OperationId(901), motion).unwrap();
            let ids =
                |s: &rusty_occt::Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
            assert_eq!(ids(piece), ids(&moved), "{name}");
            assert_eq!(h.relations.len(), ids(piece).len());
            let (a, b) = (
                piece.mass_properties().volume,
                moved.mass_properties().volume,
            );
            assert!((a - b).abs() <= 1e-9 * a, "{name}: {a} {b}");
        }
    }
}

/// S8a.2: a piece classifies points by its footprint, the prism's height
/// and the plane's side: its centre of mass inside a convex piece, points
/// across the plane outside, points on the cut face on the boundary.
#[test]
fn oblique_pieces_classify() {
    use rusty_occt::Location;
    let c = case("sq_tilted");
    let (_, pieces, _) = split(&c).unwrap();
    let plane = plane_frame(c.plane);
    for (side, piece) in &pieces {
        let centre = piece.mass_properties().centroid;
        assert_eq!(piece.classify(centre).unwrap(), Location::Inside);
        let n = plane.normal() * if *side == Side::Below { 1.0 } else { -1.0 };
        // Mirrored through the plane: the other side.
        let d = (centre - plane.origin()).dot(n);
        let across = centre + n * (-2.0 * d);
        assert_eq!(piece.classify(across).unwrap(), Location::Outside);
        let on = centre + n * (-d);
        assert_eq!(piece.classify(on).unwrap(), Location::Boundary);
    }
}

/// S8a.2: oblique pieces tessellate within the request, ellipse edges and
/// their sinusoid pcurves included.
#[test]
fn oblique_pieces_tessellate() {
    let parameters = rusty_occt::tessellation::Parameters::new(0.05, 0.5).unwrap();
    for name in [
        "sq_tilted",
        "stadium_tilted",
        "disc_wall_ellipse",
        "disc_through_caps",
    ] {
        let (_, pieces, _) = split(&case(name)).unwrap();
        for (_, piece) in &pieces {
            let mesh = piece
                .tessellate(parameters)
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(mesh.deflection <= 0.05, "{name}: {}", mesh.deflection);
            assert!(!mesh.triangles.is_empty());
        }
    }
}

/// S8a.2: a piece with ellipse edges writes only where OCCT has records:
/// planar pieces write, a cylinder's section (a sinusoid pcurve) does not.
#[test]
fn oblique_pieces_write_where_occt_has_records() {
    let (_, pieces, _) = split(&case("sq_tilted")).unwrap();
    for (_, piece) in &pieces {
        rusty_occt::occt_brep::write(piece.topology(), 1e-7).unwrap();
    }
    let (_, pieces, _) = split(&case("disc_wall_ellipse")).unwrap();
    for (_, piece) in &pieces {
        assert!(rusty_occt::occt_brep::write(piece.topology(), 1e-7).is_err());
    }
}

/// S8c's fixtures: every cone, frustum, sphere and zone builds, and its
/// certified volume contains the reference's total of both sides.
#[test]
fn primitive_fixtures_build_with_the_reference_volume() {
    let mut totals: BTreeMap<String, f64> = BTreeMap::new();
    for line in include_str!("../../fixtures/split-primitive-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        let v: f64 = row.split(' ').nth(2).unwrap().parse().unwrap();
        *totals.entry(name.to_string()).or_default() += v;
    }
    let cases = primitive_cases(include_str!("../../fixtures/split-primitive-cases.txt"));
    assert_eq!(cases.len(), totals.len());
    for case in cases {
        let solid = build_primitive(&case);
        let m = solid.topology().mass_enclosure().unwrap();
        let v = totals[&case.name];
        let slack = 1e-12 * v;
        assert!(
            m.volume[0] - slack <= v && v <= m.volume[1] + slack,
            "{}: {:?} {v}",
            case.name,
            m.volume
        );
        // Its split answers (pieces, or a later sub-step's `unsupported`).
        let rows = protocol::primitive_rows(&case).unwrap_or_else(|e| panic!("{}: {e}", case.name));
        assert!(!rows.is_empty(), "{}", case.name);
    }
}

/// A side's summed enclosures: volume, area and the three first moments.
type SideSums = ([f64; 2], [f64; 2], [[f64; 2]; 3]);

/// The primitive cases' expected sides: (side, [volume, area, cx, cy, cz]).
fn primitive_expected(text: &str) -> BTreeMap<String, Vec<(String, [f64; 5])>> {
    let mut out: BTreeMap<String, Vec<(String, [f64; 5])>> = BTreeMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let (name, row) = line.split_once('\t').unwrap();
        let w: Vec<&str> = row.split(' ').collect();
        let v: Vec<f64> = w[2..7].iter().map(|x| x.parse().unwrap()).collect();
        out.entry(name.to_string())
            .or_default()
            .push((w[1].to_string(), [v[0], v[1], v[2], v[3], v[4]]));
    }
    out
}

/// S8c: every cone, frustum, sphere and zone the kernel splits has each
/// side's volume and area inside the sums of its pieces' enclosures and
/// its centre inside theirs (the conic and circle cases since S8d.2).
#[test]
fn primitive_splits_match_the_reference() {
    check_primitive_sides(
        include_str!("../../fixtures/split-primitive-cases.txt"),
        include_str!("../../fixtures/split-primitive-expected.tsv"),
        &[],
    );
}

/// S8d.2: cones, frusta, zones and caps by planes neither normal to their
/// axes nor containing them: every conic configuration (closed, one rim,
/// both rims, touching, the rulings through a virtual apex) likewise.
#[test]
fn conic_splits_match_the_reference() {
    check_primitive_sides(
        include_str!("../../fixtures/split-conic-cases.txt"),
        include_str!("../../fixtures/split-conic-expected.tsv"),
        &[],
    );
}

/// S8d.1: every whole torus the kernel splits (normal to its axis or
/// containing it) likewise; the spiric sections wait for S8d.3.
#[test]
fn torus_splits_match_the_reference() {
    check_primitive_sides(
        include_str!("../../fixtures/split-torus-cases.txt"),
        include_str!("../../fixtures/split-torus-expected.tsv"),
        &[
            "torus_parallel_outer",
            "torus_parallel_inner",
            "torus_oblique",
            "torus_oblique_tube",
        ],
    );
}

/// S8d.3: tori cut in spiric sections (bands, caps, C-shaped pieces),
/// refused until their code lands.
#[test]
fn spiric_splits_match_the_reference() {
    check_primitive_sides(
        include_str!("../../fixtures/split-spiric-cases.txt"),
        include_str!("../../fixtures/split-spiric-expected.tsv"),
        &[
            "torus_gentle",
            "torus_band_offset",
            "torus_top_band",
            "torus_top_cap",
            "torus_cap_outer",
            "torus_peanut",
            "torus_hole_slice",
            "torus_two_ovals",
            "torus_steep",
            "torus_skew_ovals",
            "torus_small_tube",
            "torus_tilted_spiric",
            "torus_tilted_band",
        ],
    );
}

fn check_primitive_sides(cases: &str, expected: &str, later: &[&str]) {
    let want = primitive_expected(expected);
    let mut failures = Vec::new();
    for case in primitive_cases(cases) {
        let name = case.name.clone();
        let pieces = match primitive_split(&case) {
            Ok((_, pieces, _)) => pieces,
            Err(rusty_occt::Error::OutOfDomain(_)) if later.contains(&name.as_str()) => continue,
            Err(e) => {
                failures.push(format!("{name}: {e}"));
                continue;
            }
        };
        let rows = &want[&name];
        let mut sides: BTreeMap<String, SideSums> = BTreeMap::new();
        for (side, piece) in &pieces {
            let key = if rows[0].0 == "whole" {
                "whole".to_string()
            } else {
                format!("{side:?}").to_lowercase()
            };
            let m = piece.topology().mass_enclosure().unwrap();
            let e = sides
                .entry(key)
                .or_insert(([0.0; 2], [0.0; 2], [[0.0; 2]; 3]));
            e.0 = [e.0[0] + m.volume[0], e.0[1] + m.volume[1]];
            e.1 = [e.1[0] + m.surface_area[0], e.1[1] + m.surface_area[1]];
            for i in 0..3 {
                let p = [
                    m.volume[0] * m.centroid[i][0],
                    m.volume[0] * m.centroid[i][1],
                    m.volume[1] * m.centroid[i][0],
                    m.volume[1] * m.centroid[i][1],
                ];
                e.2[i] = [
                    e.2[i][0] + p.iter().copied().fold(f64::MAX, f64::min),
                    e.2[i][1] + p.iter().copied().fold(f64::MIN, f64::max),
                ];
            }
        }
        for (side, v) in rows {
            let Some((vol, area, moments)) = sides.get(side) else {
                failures.push(format!("{name}: no {side}"));
                continue;
            };
            let inside = |x: f64, [lo, hi]: [f64; 2], rel: f64| {
                lo - rel * x.abs().max(1.0) <= x && x <= hi + rel * x.abs().max(1.0)
            };
            if !inside(v[0], *vol, 1e-20)
                || !inside(v[1], *area, 1e-20)
                || (0..3).any(|i| !inside(v[0] * v[2 + i], moments[i], 1e-12))
            {
                failures.push(format!("{name} {side}: {vol:?} {area:?} miss {v:?}"));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// S8c and S8d: the primitives' and tori's split histories pass the
/// independent check, cover every input entity and repeat exactly.
#[test]
fn primitive_histories_are_complete_and_deterministic() {
    let text = [
        include_str!("../../fixtures/split-primitive-cases.txt"),
        include_str!("../../fixtures/split-conic-cases.txt"),
        include_str!("../../fixtures/split-torus-cases.txt"),
    ]
    .concat();
    for case in primitive_cases(&text) {
        let Ok((solid, pieces, h)) = primitive_split(&case) else {
            continue;
        };
        let name = &case.name;
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
        let (_, again, h2) = primitive_split(&case).unwrap();
        assert_eq!(h, h2, "{name}");
        assert_eq!(pieces.len(), again.len(), "{name}");
    }
}

/// S8d's fixtures: every whole torus builds, and its certified volume
/// contains the reference's total of both sides.
#[test]
fn torus_fixtures_build_with_the_reference_volume() {
    let mut totals: BTreeMap<String, f64> = BTreeMap::new();
    for line in include_str!("../../fixtures/split-torus-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let (name, row) = line.split_once('\t').unwrap();
        let v: f64 = row.split(' ').nth(2).unwrap().parse().unwrap();
        *totals.entry(name.to_string()).or_default() += v;
    }
    let cases = primitive_cases(include_str!("../../fixtures/split-torus-cases.txt"));
    assert_eq!(cases.len(), totals.len());
    for case in cases {
        let solid = build_primitive(&case);
        let m = solid.topology().mass_enclosure().unwrap();
        let v = totals[&case.name];
        assert!(
            m.volume[0] - 1e-12 * v <= v && v <= m.volume[1] + 1e-12 * v,
            "{}: {:?} {v}",
            case.name,
            m.volume
        );
        let rows = protocol::primitive_rows(&case).unwrap_or_else(|e| panic!("{}: {e}", case.name));
        assert!(!rows.is_empty(), "{}", case.name);
    }
}
