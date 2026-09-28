//! S9a: Booleans of prisms in one frame, against the independent reference
//! (`fixtures/boolean-*.txt|tsv` from `tools/generate_boolean_fixtures.py`)
//! and on hand-made cases.
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history::{self, History, Relation};
use rusty_occt::identity::OperationId;
use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Boundary {
    Boundary::polygon(
        vec![
            Point2::new(x0, y0),
            Point2::new(x1, y0),
            Point2::new(x1, y1),
            Point2::new(x0, y1),
        ],
        tol(),
    )
    .unwrap()
}

fn prism(outer: Boundary, holes: Vec<Boundary>, frame: Frame3, lo: f64, hi: f64, op: u64) -> Solid {
    let profile = Profile::new(outer, holes, tol()).unwrap();
    Solid::extrude_with(OperationId(op), profile, frame, lo, hi)
        .unwrap()
        .0
}

fn check(a: &Solid, b: &Solid, out: &[Solid], h: &History) {
    let ins = [
        a.topology().entity_set(a.resolution()),
        b.topology().entity_set(b.resolution()),
    ];
    let outs: Vec<_> = out
        .iter()
        .map(|s| s.topology().entity_set(s.resolution()))
        .collect();
    let issues = history::check(&ins, &outs, h);
    assert!(issues.is_empty(), "{issues:?}");
}

fn volume(out: &[Solid]) -> f64 {
    out.iter().map(|s| s.mass_properties().volume).sum()
}

#[test]
fn overlapping_boxes() {
    let a = prism(rect(0.0, 0.0, 4.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(rect(2.0, 1.0, 6.0, 3.0), vec![], Frame3::xy(), 0.0, 1.0, 2);
    let (f, h) = a.fuse(OperationId(3), &b).unwrap();
    assert_eq!(f.len(), 1);
    assert!((volume(&f) - 14.0).abs() < 1e-12, "{}", volume(&f));
    check(&a, &b, &f, &h);
    let (c, h) = a.cut(OperationId(4), &b).unwrap();
    assert_eq!(c.len(), 1);
    assert!((volume(&c) - 6.0).abs() < 1e-12, "{}", volume(&c));
    check(&a, &b, &c, &h);
    let (m, h) = a.common(OperationId(5), &b).unwrap();
    assert_eq!(m.len(), 1);
    assert!((volume(&m) - 2.0).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
}

#[test]
fn a_box_and_a_cylinder() {
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 2.0, 1);
    let circle = Boundary::circle(Point2::new(4.0, 2.0), 1.0, tol()).unwrap();
    let b = prism(circle, vec![], Frame3::xy(), -1.0, 3.0, 2);
    let half = std::f64::consts::PI / 2.0;
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!((volume(&m) - 2.0 * half).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
    let (c, h) = a.cut(OperationId(4), &b).unwrap();
    assert!(
        (volume(&c) - (32.0 - 2.0 * half)).abs() < 1e-12,
        "{}",
        volume(&c)
    );
    check(&a, &b, &c, &h);
}

#[test]
fn a_lens_of_two_cylinders() {
    let c1 = Boundary::circle(Point2::new(0.0, 0.0), 1.0, tol()).unwrap();
    let c2 = Boundary::circle(Point2::new(1.0, 0.0), 1.0, tol()).unwrap();
    let a = prism(c1, vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(c2, vec![], Frame3::xy(), 0.0, 1.0, 2);
    // Two unit circles one apart: the lens is 2π/3 - √3/2.
    let lens = 2.0 * std::f64::consts::PI / 3.0 - 3f64.sqrt() / 2.0;
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!((volume(&m) - lens).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
    let (f, h) = a.fuse(OperationId(4), &b).unwrap();
    assert!((volume(&f) - (2.0 * std::f64::consts::PI - lens)).abs() < 1e-12);
    check(&a, &b, &f, &h);
}

#[test]
fn a_hole_through_a_box() {
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(
        Boundary::circle(Point2::new(2.0, 2.0), 1.0, tol()).unwrap(),
        vec![],
        Frame3::xy(),
        -1.0,
        2.0,
        2,
    );
    let (c, h) = a.cut(OperationId(3), &b).unwrap();
    assert_eq!(c.len(), 1);
    assert_eq!(c[0].profile().unwrap().holes().len(), 1);
    assert!((volume(&c) - (16.0 - std::f64::consts::PI)).abs() < 1e-12);
    check(&a, &b, &c, &h);
}

#[test]
fn identical_and_disjoint_and_inside() {
    let a = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let same = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 2);
    let (f, h) = a.fuse(OperationId(3), &same).unwrap();
    assert!((volume(&f) - 4.0).abs() < 1e-12);
    check(&a, &same, &f, &h);
    let (c, h) = a.cut(OperationId(4), &same).unwrap();
    assert!(c.is_empty());
    check(&a, &same, &c, &h);
    let far = prism(rect(5.0, 0.0, 6.0, 1.0), vec![], Frame3::xy(), 0.0, 1.0, 5);
    let (f, h) = a.fuse(OperationId(6), &far).unwrap();
    assert_eq!(f.len(), 2);
    check(&a, &far, &f, &h);
    let (m, h) = a.common(OperationId(7), &far).unwrap();
    assert!(m.is_empty());
    check(&a, &far, &m, &h);
    let inner = prism(rect(0.5, 0.5, 1.5, 1.5), vec![], Frame3::xy(), 0.0, 1.0, 8);
    let (f, h) = a.fuse(OperationId(9), &inner).unwrap();
    assert!((volume(&f) - 4.0).abs() < 1e-12);
    check(&a, &inner, &f, &h);
    let (m, h) = a.common(OperationId(10), &inner).unwrap();
    assert!((volume(&m) - 1.0).abs() < 1e-12);
    check(&a, &inner, &m, &h);
}

#[test]
fn boxes_sharing_part_of_a_wall() {
    let a = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(rect(2.0, 1.0, 4.0, 3.0), vec![], Frame3::xy(), 0.0, 1.0, 2);
    let (f, h) = a.fuse(OperationId(3), &b).unwrap();
    assert_eq!(f.len(), 1);
    assert!((volume(&f) - 8.0).abs() < 1e-12);
    check(&a, &b, &f, &h);
}

#[test]
fn frames_with_equal_axes() {
    // An offset origin whose coordinates are binary64 in the axes.
    let other = Frame3::new(
        Point3::new(1.0, 0.5, 0.25),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    let a = prism(rect(0.0, 0.0, 4.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], other, -0.25, 0.75, 2);
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!((volume(&m) - 2.0 * 1.5).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
    // The tilted frame, one origin.
    let tilted = Frame3::new(
        Point3::new(1.0, -2.0, 0.5),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    let a = prism(rect(0.0, 0.0, 4.0, 2.0), vec![], tilted, 0.0, 1.0, 4);
    let b = prism(rect(1.0, 0.5, 3.0, 2.5), vec![], tilted, 0.0, 1.0, 5);
    let (f, h) = a.fuse(OperationId(6), &b).unwrap();
    assert!(
        (volume(&f) - (8.0 + 4.0 - 3.0)).abs() < 1e-9,
        "{}",
        volume(&f)
    );
    check(&a, &b, &f, &h);
    // An offset that rounds in the tilted axes: S9b's polyhedra.
    let shifted = Frame3::new(
        tilted.point(Point2::new(1.0, 0.5), 0.25),
        Vec3::new(0.0, 3.0, 4.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol(),
    )
    .unwrap();
    let c = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], shifted, 0.0, 1.0, 7);
    let (m, h) = a.common(OperationId(8), &c).unwrap();
    assert!((volume(&m) - 2.25).abs() < 1e-9, "{}", volume(&m));
    check(&a, &c, &m, &h);
}

#[test]
fn an_arc_profile_cut_by_a_box() {
    let stadium = Boundary::path(
        vec![
            Point2::new(0.0, -1.0),
            Point2::new(3.0, -1.0),
            Point2::new(3.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        vec![
            Segment::Line,
            Segment::Arc {
                center: Point2::new(3.0, 0.0),
                radius: 1.0,
                ccw: true,
            },
            Segment::Line,
            Segment::Arc {
                center: Point2::new(0.0, 0.0),
                radius: 1.0,
                ccw: true,
            },
        ],
        tol(),
    )
    .unwrap();
    let a = prism(stadium, vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(
        rect(1.0, -2.0, 2.0, 2.0),
        vec![],
        Frame3::xy(),
        -1.0,
        2.0,
        2,
    );
    let (c, h) = a.cut(OperationId(3), &b).unwrap();
    assert_eq!(c.len(), 2);
    let area = 6.0 + std::f64::consts::PI;
    assert!((volume(&c) - (area - 2.0)).abs() < 1e-12, "{}", volume(&c));
    check(&a, &b, &c, &h);
}

fn circle(x: f64, y: f64, r: f64) -> Boundary {
    Boundary::circle(Point2::new(x, y), r, tol()).unwrap()
}

#[test]
fn tangent_cylinders() {
    use rusty_occt::Error;
    // Outside each other, touching at (1, 0): a fuse would touch itself
    // along a line, the common is empty.
    let a = prism(circle(0.0, 0.0, 1.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(circle(2.0, 0.0, 1.0), vec![], Frame3::xy(), 0.0, 1.0, 2);
    assert!(matches!(
        a.fuse(OperationId(3), &b),
        Err(Error::Degenerate(_))
    ));
    let (m, h) = a.common(OperationId(4), &b).unwrap();
    assert!(m.is_empty());
    check(&a, &b, &m, &h);
    let (c, h) = a.cut(OperationId(5), &b).unwrap();
    assert!((volume(&c) - std::f64::consts::PI).abs() < 1e-12);
    check(&a, &b, &c, &h);
    // Inside, touching at (2, 0): the fuse is the big one, the common the
    // small one, the cut would touch itself.
    let big = prism(circle(0.0, 0.0, 2.0), vec![], Frame3::xy(), 0.0, 1.0, 6);
    let small = prism(circle(1.0, 0.0, 1.0), vec![], Frame3::xy(), 0.0, 1.0, 7);
    let (f, h) = big.fuse(OperationId(8), &small).unwrap();
    assert!((volume(&f) - 4.0 * std::f64::consts::PI).abs() < 1e-12);
    check(&big, &small, &f, &h);
    let (m, h) = big.common(OperationId(9), &small).unwrap();
    assert!((volume(&m) - std::f64::consts::PI).abs() < 1e-12);
    check(&big, &small, &m, &h);
    assert!(matches!(
        big.cut(OperationId(10), &small),
        Err(Error::Degenerate(_))
    ));
}

#[test]
fn arcs_of_one_circle() {
    // Two half-discs of the unit circle overlapping in a quarter.
    let half = |a0: f64, op: u64| {
        let (c, s) = (a0.cos(), a0.sin());
        let p = Point2::new(c, s);
        let q = Point2::new(-c, -s);
        let boundary = Boundary::path(
            vec![p, q],
            vec![
                Segment::Arc {
                    center: Point2::new(0.0, 0.0),
                    radius: 1.0,
                    ccw: true,
                },
                Segment::Line,
            ],
            tol(),
        )
        .unwrap();
        prism(boundary, vec![], Frame3::xy(), 0.0, 1.0, op)
    };
    let a = half(0.0, 1);
    let b = half(std::f64::consts::FRAC_PI_2, 2);
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!(
        (volume(&m) - std::f64::consts::FRAC_PI_4).abs() < 1e-9,
        "{}",
        volume(&m)
    );
    check(&a, &b, &m, &h);
    let (f, h) = a.fuse(OperationId(4), &b).unwrap();
    assert!(
        (volume(&f) - 3.0 * std::f64::consts::FRAC_PI_4).abs() < 1e-9,
        "{}",
        volume(&f)
    );
    check(&a, &b, &f, &h);
}

#[test]
fn a_vertex_on_the_others_edge() {
    // A triangle whose apex lies on the square's top edge.
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let tri = Boundary::polygon(
        vec![
            Point2::new(1.0, 6.0),
            Point2::new(3.0, 6.0),
            Point2::new(2.0, 4.0),
        ],
        tol(),
    )
    .unwrap();
    let b = prism(tri, vec![], Frame3::xy(), 0.0, 1.0, 2);
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert!(m.is_empty());
    check(&a, &b, &m, &h);
    let tri = Boundary::polygon(
        vec![
            Point2::new(1.0, 2.0),
            Point2::new(3.0, 2.0),
            Point2::new(2.0, 4.0),
        ],
        tol(),
    )
    .unwrap();
    let inner = prism(tri, vec![], Frame3::xy(), 0.0, 1.0, 4);
    let (m, h) = a.common(OperationId(5), &inner).unwrap();
    assert!((volume(&m) - 2.0).abs() < 1e-12);
    check(&a, &inner, &m, &h);
}

#[test]
fn a_holed_box_and_a_tool() {
    let a = prism(
        rect(0.0, 0.0, 6.0, 6.0),
        vec![circle(3.0, 3.0, 1.0)],
        Frame3::xy(),
        0.0,
        1.0,
        1,
    );
    // Through the hole: the common is the box's part around it, in two.
    let b = prism(
        rect(2.5, -1.0, 3.5, 7.0),
        vec![],
        Frame3::xy(),
        -1.0,
        2.0,
        2,
    );
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    assert_eq!(m.len(), 2);
    let band = 0.75f64.sqrt() + std::f64::consts::PI / 3.0;
    let expected = 6.0 - band;
    assert!((volume(&m) - expected).abs() < 1e-12, "{}", volume(&m));
    check(&a, &b, &m, &h);
    let (c, h) = a.cut(OperationId(4), &b).unwrap();
    assert_eq!(c.len(), 2);
    assert!((volume(&c) - (36.0 - std::f64::consts::PI - expected)).abs() < 1e-12);
    check(&a, &b, &c, &h);
}

/// Faces, edges and vertices of each solid.
fn counts(out: &[Solid]) -> Vec<[usize; 3]> {
    out.iter()
        .map(|s| {
            let t = s.topology();
            [t.faces().len(), t.edges().len(), t.vertices().len()]
        })
        .collect()
}

#[test]
fn stacks_of_boxes() {
    // S9a.2: a tower on a box, and a pocket in it.
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 1.0, 1);
    let b = prism(rect(1.0, 1.0, 3.0, 3.0), vec![], Frame3::xy(), 0.5, 2.0, 2);
    let (f, h) = a.fuse(OperationId(3), &b).unwrap();
    assert!((volume(&f) - 20.0).abs() < 1e-12, "{}", volume(&f));
    assert_eq!(counts(&f), [[11, 24, 16]]);
    check(&a, &b, &f, &h);
    let (c, h) = a.cut(OperationId(4), &b).unwrap();
    assert!((volume(&c) - 14.0).abs() < 1e-12, "{}", volume(&c));
    assert_eq!(counts(&c), [[11, 24, 16]]);
    check(&a, &b, &c, &h);
    let (m, h) = a.common(OperationId(5), &b).unwrap();
    assert!((volume(&m) - 2.0).abs() < 1e-12);
    check(&a, &b, &m, &h);
    // The stack classifies, moves rigidly with its ids and volume.
    use rusty_occt::Location;
    let at = |x, y, z| c[0].classify(Point3::new(x, y, z)).unwrap();
    assert_eq!(at(2.0, 2.0, 0.25), Location::Inside);
    assert_eq!(at(2.0, 2.0, 0.75), Location::Outside);
    assert_eq!(at(2.0, 2.0, 0.5), Location::Boundary);
    assert_eq!(at(0.5, 0.5, 0.75), Location::Inside);
    assert_eq!(at(0.5, 0.5, 1.0), Location::Boundary);
    let at = |x, y, z| f[0].classify(Point3::new(x, y, z)).unwrap();
    assert_eq!(at(2.0, 2.0, 1.5), Location::Inside);
    assert_eq!(at(0.5, 0.5, 1.5), Location::Outside);
    assert_eq!(at(3.0, 2.0, 1.5), Location::Boundary);
    let turn = rusty_occt::RigidTransform::rotation(
        Point3::new(1.0, 0.0, 0.0),
        Vec3::new(1.0, 2.0, 2.0),
        0.5,
    )
    .unwrap();
    let (moved, _) = f[0].transform_with(OperationId(6), turn).unwrap();
    assert!((moved.mass_properties().volume - 20.0).abs() < 1e-9);
    let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    assert_eq!(ids(&moved), ids(&f[0]));
}

/// Per case: its declared kind and the reference's solid count and totals.
type Expected = std::collections::BTreeMap<String, (String, Option<(usize, [f64; 5])>)>;

/// Every fixture: a result's solids' sums of enclosures contain the
/// reference's volume, area and moments; its solid count is the
/// reference's; a degenerate case is refused and an S9a.2 stack is
/// `OutOfDomain`.
#[test]
fn every_case_matches_the_reference() {
    let mut expect: Expected = Default::default();
    for line in include_str!("../../fixtures/boolean-expected.tsv")
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
                let v: Vec<f64> = w[2..7].iter().map(|x| x.parse().unwrap()).collect();
                entry.1 = Some((w[1].parse().unwrap(), [v[0], v[1], v[2], v[3], v[4]]));
            }
            _ => {}
        }
    }
    let mut failures = Vec::new();
    for case in protocol::cases(include_str!("../../fixtures/boolean-cases.txt")) {
        let (kind, want) = &expect[&case.name];
        let rows = protocol::rows(&case).unwrap_or_else(|e| panic!("{}: {e}", case.name));
        match kind.as_str() {
            "degenerate" => {
                if rows != ["refused"] {
                    failures.push(format!("{}: {rows:?} not refused", case.name));
                }
                continue;
            }
            "empty" => {
                if rows != ["empty"] {
                    failures.push(format!("{}: {rows:?} not empty", case.name));
                }
                continue;
            }
            _ => {}
        }
        let (count, v) = want.expect("a result");
        let solids: Vec<Vec<f64>> = rows
            .iter()
            .map(|r| {
                r.split(' ')
                    .skip(1)
                    .take(10)
                    .map(|x| x.parse().unwrap())
                    .collect()
            })
            .collect();
        if solids.len() != count {
            failures.push(format!(
                "{}: {} solids for {count}",
                case.name,
                solids.len()
            ));
            continue;
        }
        let sum = |i: usize| solids.iter().map(|s| s[i]).sum::<f64>();
        let inside =
            |x: f64, lo: f64, hi: f64| lo - 1e-20 * x.abs() <= x && x <= hi + 1e-20 * x.abs();
        if !inside(v[0], sum(0), sum(1)) || !inside(v[1], sum(2), sum(3)) {
            failures.push(format!("{}: measures miss {v:?}", case.name));
        }
        for i in 0..3 {
            let (mut lo, mut hi) = (0.0, 0.0);
            for s in &solids {
                let p = [
                    s[0] * s[4 + 2 * i],
                    s[0] * s[5 + 2 * i],
                    s[1] * s[4 + 2 * i],
                    s[1] * s[5 + 2 * i],
                ];
                lo += p.iter().copied().fold(f64::MAX, f64::min);
                hi += p.iter().copied().fold(f64::MIN, f64::max);
            }
            let m = v[0] * v[2 + i];
            let allow = 1e-12 * m.abs().max(1.0);
            if !(lo - allow <= m && m <= hi + allow) {
                failures.push(format!("{}: moment {i} misses", case.name));
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// Every fixture's history passes the independent check, covers every input
/// entity and repeats exactly.
#[test]
fn fixture_histories_are_complete_and_deterministic() {
    for case in protocol::cases(include_str!("../../fixtures/boolean-cases.txt")) {
        let Ok((a, b, out, h)) = protocol::run(&case) else {
            continue;
        };
        check(&a, &b, &out, &h);
        let covered: std::collections::BTreeSet<_> =
            h.relations.iter().flat_map(Relation::sources).collect();
        for s in [&a, &b] {
            for (id, _) in s.topology().ids() {
                assert!(
                    covered.contains(&id),
                    "{}: {id:?} has no relation",
                    case.name
                );
            }
        }
        let (_, _, again, h2) = protocol::run(&case).unwrap();
        assert_eq!(h, h2, "{}", case.name);
        let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
        for (x, y) in out.iter().zip(&again) {
            assert_eq!(ids(x), ids(y), "{}", case.name);
        }
    }
}

/// A fuse of meeting ranges is one input only when the other lies inside
/// it: a tool filling the object's hole over part of the height, or a
/// profile with the same outer and a hole, is a stack (a fuzz find: the
/// fused boundary came from one input alone, its hole's filled).
#[test]
fn a_fuse_filling_a_hole_is_a_stack() {
    let holed = prism(
        rect(-4.0, -4.0, 4.0, 4.0),
        vec![circle(0.0, 0.0, 1.5)],
        Frame3::xy(),
        -1.0,
        3.0,
        1,
    );
    let plug = prism(circle(0.0, 0.0, 2.0), vec![], Frame3::xy(), 0.0, 2.0, 2);
    // S9a.2: the plug fills the hole over part of its height.
    let pi = std::f64::consts::PI;
    for (f, h) in [
        plug.fuse(OperationId(3), &holed).unwrap(),
        holed.fuse(OperationId(4), &plug).unwrap(),
    ] {
        assert_eq!(f.len(), 1);
        assert!(
            (volume(&f) - (256.0 - 4.5 * pi)).abs() < 1e-9,
            "{}",
            volume(&f)
        );
        check(&plug, &holed, &f, &h);
    }
    let solid = prism(
        rect(-4.0, -4.0, 4.0, 4.0),
        vec![],
        Frame3::xy(),
        0.0,
        1.0,
        5,
    );
    let (f, h) = solid.fuse(OperationId(6), &holed).unwrap();
    assert!(
        (volume(&f) - (256.0 - 6.75 * pi)).abs() < 1e-9,
        "{}",
        volume(&f)
    );
    check(&solid, &holed, &f, &h);
    // Inside in 2D and in height: the holder, unchanged.
    let small = prism(rect(2.0, 2.0, 3.0, 3.0), vec![], Frame3::xy(), 0.0, 1.0, 7);
    let (f, h) = holed.fuse(OperationId(8), &small).unwrap();
    assert_eq!(f.len(), 1);
    assert_eq!(f[0].topology().body_id(), holed.topology().body_id());
    check(&holed, &small, &f, &h);
}

#[test]
fn a_result_is_an_input_again() {
    // A Boolean's result, renamed, is a prism another Boolean takes.
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 2.0, 1);
    let b = prism(rect(2.0, 2.0, 6.0, 6.0), vec![], Frame3::xy(), 0.0, 2.0, 2);
    let (f, _) = a.fuse(OperationId(3), &b).unwrap();
    let c = prism(rect(1.0, 1.0, 5.0, 5.0), vec![], Frame3::xy(), -1.0, 3.0, 4);
    let (cut, h) = f[0].cut(OperationId(5), &c).unwrap();
    assert!((volume(&cut) - 28.0).abs() < 1e-12, "{}", volume(&cut));
    check(&f[0], &c, &cut, &h);
    let up = rusty_occt::RigidTransform::translation(Vec3::new(0.0, 0.0, 0.5)).unwrap();
    let moved = cut[0].transform_with(OperationId(7), up).unwrap().0;
    let e = prism(rect(0.0, 0.0, 6.0, 6.0), vec![], Frame3::xy(), 0.5, 2.5, 8);
    let (again, h) = moved.fuse(OperationId(9), &e).unwrap();
    assert!((volume(&again) - 72.0).abs() < 1e-12, "{}", volume(&again));
    check(&moved, &e, &again, &h);
    // A result and an input it keeps entities of share those ids.
    assert!(matches!(
        moved.fuse(OperationId(10), &f[0]),
        Err(rusty_occt::Error::InvalidLabel(_))
    ));
}

#[test]
fn inputs_sharing_ids_are_refused() {
    // Built by one operation, two prisms share every id: the history could
    // not tell them apart.
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 2.0, 1);
    let b = prism(rect(2.0, 2.0, 6.0, 6.0), vec![], Frame3::xy(), 0.0, 2.0, 1);
    for r in [a.fuse(OperationId(3), &b), a.cut(OperationId(3), &a)] {
        assert!(
            matches!(r, Err(rusty_occt::Error::InvalidLabel(_))),
            "{:?}",
            r.map(|x| x.0.len())
        );
    }
}

#[test]
fn a_cut_leaving_a_cavity() {
    // S9a.2: a tool inside the object in 2D and in height leaves a closed
    // void: one solid of two shells, the void a bounded region.
    use rusty_occt::topology::RegionKind;
    use rusty_occt::Location;
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 4.0, 1);
    let b = prism(circle(2.0, 2.0, 1.0), vec![], Frame3::xy(), 1.0, 3.0, 2);
    let (c, h) = a.cut(OperationId(3), &b).unwrap();
    assert_eq!(c.len(), 1);
    let v = 64.0 - 2.0 * std::f64::consts::PI;
    assert!((volume(&c) - v).abs() < 1e-9, "{}", volume(&c));
    check(&a, &b, &c, &h);
    let t = c[0].topology();
    let kinds: Vec<RegionKind> = t.regions().iter().map(|r| r.kind).collect();
    assert_eq!(
        kinds,
        [RegionKind::Void, RegionKind::Solid, RegionKind::Void]
    );
    assert_eq!(t.regions()[1].shells.len(), 2);
    let at = |x, y, z| c[0].classify(Point3::new(x, y, z)).unwrap();
    assert_eq!(at(2.0, 2.0, 2.0), Location::Outside);
    assert_eq!(at(2.0, 2.0, 0.5), Location::Inside);
    assert_eq!(at(2.0, 2.0, 1.0), Location::Boundary);
    assert_eq!(at(3.0, 2.0, 2.0), Location::Boundary);
    // The fuse is the object.
    let (f, h) = a.fuse(OperationId(4), &b).unwrap();
    assert_eq!(f[0].topology().body_id(), a.topology().body_id());
    check(&a, &b, &f, &h);
}

#[test]
fn a_tool_through_a_round_wall() {
    // Fuzzing (boolean, S9a.2): a holed square over the middle of a
    // stadium's height leaves a window in its round wall, a loop on the
    // cylinder's sheet of the wall's outer loop.
    let (s, t) = (1.5, 2.375);
    let stadium = Boundary::path(
        vec![
            Point2::new(0.0, -t),
            Point2::new(s, -t),
            Point2::new(s, t),
            Point2::new(0.0, t),
        ],
        vec![
            Segment::Line,
            Segment::Arc {
                center: Point2::new(s, 0.0),
                radius: t,
                ccw: true,
            },
            Segment::Line,
            Segment::Arc {
                center: Point2::new(0.0, 0.0),
                radius: t,
                ccw: true,
            },
        ],
        tol(),
    )
    .unwrap();
    let a = prism(stadium, vec![], Frame3::xy(), 0.0, 0.5, 1);
    let b = prism(
        rect(-5.5, -6.0, 2.5, 2.0),
        vec![circle(-1.5, -2.0, 1.625)],
        Frame3::xy(),
        0.125,
        0.375,
        2,
    );
    let (c, h) = a.cut(OperationId(3), &b).unwrap();
    check(&a, &b, &c, &h);
    let (f, h) = a.fuse(OperationId(4), &b).unwrap();
    check(&a, &b, &f, &h);
    let (m, h) = a.common(OperationId(5), &b).unwrap();
    check(&a, &b, &m, &h);
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    assert!((volume(&c) - (va - volume(&m))).abs() < 1e-9);
    assert!((volume(&f) - (va + vb - volume(&m))).abs() < 1e-9);
}

#[test]
fn a_box_inscribed_in_a_cylinder_over_part_of_its_height() {
    // Upstream `bopfuse_simple/Z5`: the box's corners on the cylinder within
    // the resolution, so the cylinder less the box touches itself there;
    // that refused cut still says the cylinder is not inside the box, and
    // the box lies inside the cylinder: the fuse is the cylinder.
    let r = std::f64::consts::SQRT_2 / 2.0;
    let a = prism(circle(0.0, 0.0, 1.0), vec![], Frame3::xy(), 0.0, 2.0, 1);
    let b = prism(rect(-r, -r, r, r), vec![], Frame3::xy(), 0.0, 1.0, 2);
    for (x, y) in [(&a, &b), (&b, &a)] {
        let (f, h) = x.fuse(OperationId(3), y).unwrap();
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].topology().body_id(), a.topology().body_id());
        check(x, y, &f, &h);
    }
}
