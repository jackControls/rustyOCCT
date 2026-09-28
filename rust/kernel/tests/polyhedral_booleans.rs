//! S9b: Booleans of polyhedral prisms in any relative position, against the
//! independent reference (`fixtures/boolean-polyhedra-*.txt|tsv` from
//! `tools/generate_polyhedral_fixtures.py`) and on hand-made cases.
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history::{self, History, Relation};
use rusty_occt::identity::OperationId;
use rusty_occt::topology::RegionKind;
use rusty_occt::{Boundary, Frame3, Location, Point2, Point3, Profile, Solid, Tolerance, Vec3};

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

fn frame(origin: [f64; 3], normal: [f64; 3], x: [f64; 3]) -> Frame3 {
    let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
    Frame3::new(
        Point3::new(origin[0], origin[1], origin[2]),
        v(normal),
        v(x),
        tol(),
    )
    .unwrap()
}

fn prism(outer: Boundary, holes: Vec<Boundary>, f: Frame3, lo: f64, hi: f64, op: u64) -> Solid {
    let profile = Profile::new(outer, holes, tol()).unwrap();
    Solid::extrude_with(OperationId(op), profile, f, lo, hi)
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
fn every_case_matches_the_reference() {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-polyhedra-expected.tsv")
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
    let mut failures = Vec::new();
    for case in protocol::cases(include_str!("../../fixtures/boolean-polyhedra-cases.txt")) {
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
                    .take(4)
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
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
fn fixture_histories_are_complete_and_deterministic() {
    for case in protocol::cases(include_str!("../../fixtures/boolean-polyhedra-cases.txt")) {
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

#[test]
fn a_turned_box_through_a_box() {
    // A 2 x 2 x 2 box turned about z by (0.6, 0.8) and centred at the
    // world box's corner: a quarter of it inside.
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], Frame3::xy(), 0.0, 4.0, 1);
    let turn = frame([4.0, 4.0, 1.0], [0.0, 0.0, 1.0], [3.0, 4.0, 0.0]);
    let b = prism(rect(-1.0, -1.0, 1.0, 1.0), vec![], turn, 0.0, 2.0, 2);
    let (m, h) = a.common(OperationId(3), &b).unwrap();
    check(&a, &b, &m, &h);
    let (f, h) = a.fuse(OperationId(4), &b).unwrap();
    check(&a, &b, &f, &h);
    let (c, h) = a.cut(OperationId(5), &b).unwrap();
    check(&a, &b, &c, &h);
    let (va, vb, vm) = (
        a.mass_properties().volume,
        b.mass_properties().volume,
        volume(&m),
    );
    assert!((volume(&f) - (va + vb - vm)).abs() < 1e-9);
    assert!((volume(&c) - (va - vm)).abs() < 1e-9);
    // The tool's quarter inside the box: its centre's quadrant.
    assert!((vm - 2.0).abs() < 1e-9, "{vm}");
    // A result classifies, and moves rigidly with its ids and volume.
    let at = |s: &Solid, x, y, z| s.classify(Point3::new(x, y, z)).unwrap();
    assert_eq!(at(&c[0], 3.9, 3.9, 2.0), Location::Outside);
    assert_eq!(at(&c[0], 1.0, 1.0, 2.0), Location::Inside);
    assert_eq!(at(&f[0], 4.5, 4.5, 2.0), Location::Inside);
    assert_eq!(at(&f[0], 4.0, 1.0, 2.0), Location::Boundary);
    let motion = rusty_occt::RigidTransform::rotation(
        Point3::new(1.0, 0.0, 0.0),
        Vec3::new(1.0, 2.0, 2.0),
        0.5,
    )
    .unwrap();
    let (moved, _) = f[0].transform_with(OperationId(6), motion).unwrap();
    assert!((moved.mass_properties().volume - volume(&f)).abs() < 1e-9);
    let ids = |s: &Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    assert_eq!(ids(&moved), ids(&f[0]));
    for v in moved.topology().vertices() {
        assert_eq!(moved.classify(v.position).unwrap(), Location::Boundary);
    }
}

#[test]
fn a_tilted_box_inside_leaves_a_cavity() {
    let a = prism(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![],
        Frame3::xy(),
        0.0,
        5.0,
        1,
    );
    let tilt = frame([4.0, 4.0, 1.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let b = prism(rect(0.0, 0.0, 2.0, 2.0), vec![], tilt, 0.0, 2.0, 2);
    let (c, h) = a.cut(OperationId(3), &b).unwrap();
    assert_eq!(c.len(), 1);
    check(&a, &b, &c, &h);
    let kinds: Vec<RegionKind> = c[0].topology().regions().iter().map(|r| r.kind).collect();
    assert_eq!(
        kinds,
        [RegionKind::Void, RegionKind::Solid, RegionKind::Void]
    );
    assert!((volume(&c) - 492.0).abs() < 1e-9, "{}", volume(&c));
}

#[test]
fn a_tilted_bar_cuts_a_box_in_two() {
    let a = prism(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![],
        Frame3::xy(),
        0.0,
        5.0,
        1,
    );
    let tilt = frame([-1.0, 5.0, 2.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let b = prism(rect(0.0, -1.0, 12.0, 1.0), vec![], tilt, -5.0, 5.0, 2);
    let (c, h) = a.cut(OperationId(3), &b).unwrap();
    assert_eq!(c.len(), 2);
    check(&a, &b, &c, &h);
    assert!((volume(&c) - 375.0).abs() < 1e-9, "{}", volume(&c));
}

#[test]
fn touching_solids_are_refused() {
    // Along an edge (the tool's lowest edge on the box's top face) and at a
    // vertex: two solids sharing an edge or a point.
    let a = prism(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![],
        Frame3::xy(),
        0.0,
        5.0,
        1,
    );
    let tilt = frame([3.0, 5.0, 5.0], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let edge = prism(rect(0.0, -3.0, 4.0, 0.0), vec![], tilt, 0.0, 3.0, 2);
    assert!(matches!(
        a.fuse(OperationId(3), &edge),
        Err(rusty_occt::Error::Degenerate(_))
    ));
    // The cut and the common are the box and nothing.
    let (c, h) = a.cut(OperationId(4), &edge).unwrap();
    check(&a, &edge, &c, &h);
    assert!((volume(&c) - 500.0).abs() < 1e-9);
    let (m, _) = a.common(OperationId(5), &edge).unwrap();
    assert!(m.is_empty());
}

/// The three operations of two solids, each result's history checked, and
/// the volume identities `V(A u B) = V(A) + V(B) - V(A n B)`, `V(A - B) =
/// V(A) - V(A n B)`.
fn identities(a: &Solid, b: &Solid, op: u64) -> [f64; 3] {
    let (f, h) = a.fuse(OperationId(op), b).unwrap();
    check(a, b, &f, &h);
    let (c, h) = a.cut(OperationId(op + 1), b).unwrap();
    check(a, b, &c, &h);
    let (m, h) = a.common(OperationId(op + 2), b).unwrap();
    check(a, b, &m, &h);
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let (vf, vc, vm) = (volume(&f), volume(&c), volume(&m));
    assert!((vf - (va + vb - vm)).abs() < 1e-9, "{vf} {va} {vb} {vm}");
    assert!((vc - (va - vm)).abs() < 1e-9, "{vc} {va} {vm}");
    [vf, vc, vm]
}

#[test]
fn results_are_inputs_again() {
    // S9b.2: a stack, an S9b result and a plane's piece as inputs.
    let a = prism(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![],
        Frame3::xy(),
        0.0,
        5.0,
        1,
    );
    // Two square pockets cut in turn (the first leaves a stack).
    let p1 = prism(rect(1.0, 1.0, 3.0, 3.0), vec![], Frame3::xy(), 3.0, 6.0, 2);
    let (s, _) = a.cut(OperationId(3), &p1).unwrap();
    let p2 = prism(rect(6.0, 6.0, 8.0, 8.0), vec![], Frame3::xy(), 2.0, 6.0, 4);
    let (twice, h) = s[0].cut(OperationId(5), &p2).unwrap();
    check(&s[0], &p2, &twice, &h);
    assert!(
        (volume(&twice) - (500.0 - 8.0 - 12.0)).abs() < 1e-9,
        "{}",
        volume(&twice)
    );
    // A turned box against the stack.
    let turn = frame([5.0, 5.0, 4.0], [0.0, 0.0, 1.0], [3.0, 4.0, 0.0]);
    let t = prism(rect(-2.0, -2.0, 2.0, 2.0), vec![], turn, 0.0, 3.0, 6);
    let [_, _, vm] = identities(&s[0], &t, 10);
    assert!(vm > 0.0);
    // An S9b result against a prism.
    let (poly, _) = a.fuse(OperationId(20), &t).unwrap();
    let c = prism(
        rect(-1.0, -1.0, 4.0, 4.0),
        vec![],
        Frame3::xy(),
        -1.0,
        9.0,
        21,
    );
    identities(&poly[0], &c, 30);
    // A plane's piece of a box against a box.
    let plane = Frame3::new(
        Point3::new(5.0, 5.0, 2.5),
        Vec3::new(1.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        tol(),
    )
    .unwrap();
    let (pieces, _) = a.split_by_plane(OperationId(40), plane).unwrap();
    let piece = &pieces[0].1;
    let d = prism(
        rect(4.0, 4.0, 12.0, 12.0),
        vec![],
        Frame3::xy(),
        1.0,
        4.0,
        41,
    );
    identities(piece, &d, 50);
}

#[test]
fn a_tilted_stack_with_collinear_cap_edges_is_an_input() {
    // A block on a box's side, rising above it, in a tilted frame: the box's
    // cap meets the block's wall along one frame line in three edges, a
    // vertical run in the cap's projection whose stored vertices are not
    // collinear exactly. The stored model's triangles meet at all of them.
    let tilt = frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
    let a = prism(rect(0.0, 0.0, 4.0, 4.0), vec![], tilt, 0.0, 2.0, 1);
    let block = prism(rect(4.0, 1.0, 6.0, 3.0), vec![], tilt, 1.0, 3.0, 2);
    let (stack, _) = a.fuse(OperationId(3), &block).unwrap();
    assert_eq!(stack.len(), 1);
    let turned = Frame3::new(
        tilt.point(Point2::new(4.0, 2.0), 1.5),
        tilt.normal(),
        tilt.x() * 3.0 + tilt.y() * 4.0,
        tol(),
    )
    .unwrap();
    let b = prism(rect(-1.0, -1.0, 1.0, 1.0), vec![], turned, 0.0, 2.0, 4);
    let [_, _, vm] = identities(&stack[0], &b, 10);
    assert!(vm > 0.0);
}
