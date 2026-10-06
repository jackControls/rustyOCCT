//! S9c.1: Booleans of prisms with arcs in any relative position, against
//! the independent reference (`fixtures/boolean-curved-*` from
//! `tools/generate_curved_boolean_fixtures.py`).
#[path = "support/boolean_protocol.rs"]
mod protocol;
use rusty_occt::history;

/// None left to a later sub-step.
const LATER: &[&str] = &[];

/// Equal cylinders with meeting axes in stored turned frames: the models'
/// extents across the common perpendicular are equal exactly, a double
/// tangency (S9c.2b's decisions: `Degenerate`).
const NODES: &[&str] = &[
    "steinmetz_oblique_fuse",
    "steinmetz_oblique_common",
    "steinmetz_tilted_common",
];

#[test]
fn every_case_matches_the_reference() {
    let mut expect = std::collections::BTreeMap::new();
    for line in include_str!("../../fixtures/boolean-curved-expected.tsv")
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
    for case in protocol::cases(include_str!("../../fixtures/boolean-curved-cases.txt")) {
        let (kind, want) = &expect[&case.name];
        let rows = match protocol::rows(&case) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{}: {e}", case.name));
                continue;
            }
        };
        if NODES.contains(&case.name.as_str()) {
            if rows != ["refused"] {
                failures.push(format!("{}: {rows:?} not refused", case.name));
            }
            continue;
        }
        if LATER.contains(&case.name.as_str()) {
            if rows != ["unsupported"] {
                failures.push(format!("{}: {rows:?} not unsupported", case.name));
            }
            continue;
        }
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
                    .map(|x| x.parse().unwrap_or(f64::NAN))
                    .collect()
            })
            .collect();
        if solids.len() != count || solids.iter().any(|s| s.len() < 4) {
            failures.push(format!("{}: {rows:?} for {count} solids", case.name));
            continue;
        }
        let sum = |i: usize| solids.iter().map(|s| s[i]).sum::<f64>();
        let near = |x: f64, lo: f64, hi: f64| {
            let slack = 1e-9 * x.abs().max(1.0);
            lo - slack <= x && x <= hi + slack
        };
        // The enclosures narrow: a wide one (its midpoint the reported
        // value) would hold the reference yet report another.
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
                "{}: volume {v:?} against [{}, {}], area [{}, {}]",
                case.name,
                sum(0),
                sum(1),
                sum(2),
                sum(3)
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} failures: {failures:#?}",
        failures.len()
    );
}

#[test]
fn fixture_histories_are_complete() {
    for case in protocol::cases(include_str!("../../fixtures/boolean-curved-cases.txt")) {
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
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Location, Point3, RigidTransform, Vec3};
    let motion =
        RigidTransform::rotation(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 2.0), 0.5)
            .unwrap();
    let ids = |s: &rusty_occt::Solid| s.topology().ids().map(|(id, _)| id).collect::<Vec<_>>();
    for case in protocol::cases(include_str!("../../fixtures/boolean-curved-cases.txt")) {
        let Ok((_, _, out, h)) = protocol::run(&case) else {
            continue;
        };
        let (_, _, again, h2) = protocol::run(&case).unwrap();
        assert_eq!(h, h2, "{}", case.name);
        for (x, y) in out.iter().zip(&again) {
            assert_eq!(ids(x), ids(y), "{}", case.name);
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

/// Parallel circular cylinders (exact frames) where the second input's cap
/// circle crosses the first's wall (found by S9e.1): the second's arc meets
/// the first's cylinder where its own circle meets the first's, both in its
/// own frame (read in the first's frame, its crossings were missed and the
/// result left open). The common is the circles' lens over the heights both
/// hold, whichever input comes first.
#[test]
fn a_second_inputs_arc_across_a_parallel_cylinder() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let disc = |cx: f64, cy: f64, r: f64| {
        Profile::new(
            Boundary::circle(Point2::new(cx, cy), r, tol).unwrap(),
            vec![],
            tol,
        )
        .unwrap()
    };
    let (a, _) =
        Solid::extrude_with(OperationId(1), disc(0.0, 0.0, 3.0), Frame3::xy(), 0.0, 5.0).unwrap();
    let down = Frame3::new(
        Point3::new(0.0, 0.0, 6.0),
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::new(1.0, 0.0, 0.0),
        tol,
    )
    .unwrap();
    // Over z in [2, 6], its circle about (4, 0.5) of radius 2.
    let (b, _) = Solid::extrude_with(OperationId(2), disc(4.0, -0.5, 2.0), down, 0.0, 4.0).unwrap();
    let (r1, r2, d) = (3.0f64, 2.0f64, 16.25f64.sqrt());
    let xc = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let (a1, a2) = ((xc / r1).acos(), ((d - xc) / r2).acos());
    let lens = r1 * r1 * (a1 - a1.sin() * a1.cos()) + r2 * r2 * (a2 - a2.sin() * a2.cos());
    let volume = |out: Vec<Solid>| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    for common in [
        volume(a.common(OperationId(3), &b).unwrap().0),
        volume(b.common(OperationId(4), &a).unwrap().0),
    ] {
        assert!((common - 3.0 * lens).abs() <= 1e-9 * lens, "{common}");
    }
    let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
    let fuse = volume(a.fuse(OperationId(5), &b).unwrap().0);
    let cut = volume(a.cut(OperationId(6), &b).unwrap().0);
    assert!((fuse - (va + vb - 3.0 * lens)).abs() <= 1e-9 * fuse);
    assert!((cut - (va - 3.0 * lens)).abs() <= 1e-9 * cut);
}

/// A circle as two arcs (two faces on one cylinder, a periodic face split
/// at a seam) against a box crossing the arcs' joint, in frames with equal
/// axes and an offset that rounds (the curved engine): the same solids as
/// the circle's. A point on the arcs' shared chord, run either way, was
/// taken inside both circular segments and the box's cap piece inside the
/// circle kept, the fuse left open (S9e.4a found it; latent since S9c.1).
#[test]
fn two_arcs_of_one_circle_are_its_circle() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3};
    let tol = Tolerance::default();
    let at = |x: f64, z: f64| {
        Frame3::new(
            Point3::new(x, 0.0, z),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            tol,
        )
        .unwrap()
    };
    let arc = Segment::Arc {
        center: Point2::new(5.0, 5.0),
        radius: 2.5,
        ccw: true,
    };
    let halves = Boundary::path(
        vec![Point2::new(7.5, 5.0), Point2::new(2.5, 5.0)],
        vec![arc.clone(), arc],
        tol,
    )
    .unwrap();
    let circle = Boundary::circle(Point2::new(5.0, 5.0), 2.5, tol).unwrap();
    let square = Boundary::polygon(
        vec![
            Point2::new(5.5, 2.0),
            Point2::new(9.0, 2.0),
            Point2::new(9.0, 8.0),
            Point2::new(5.5, 8.0),
        ],
        tol,
    )
    .unwrap();
    let prism = |b: Boundary, op: u64, f: Frame3, h: [f64; 2]| {
        let profile = Profile::new(b, vec![], tol).unwrap();
        Solid::extrude_with(OperationId(op), profile, f, h[0], h[1])
            .unwrap()
            .0
    };
    let tool = prism(square, 2, at(0.0, 1.0), [0.0, 3.0]);
    let volume = |out: Vec<Solid>| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    for (k, profile) in [halves, circle].into_iter().enumerate() {
        let object = prism(profile, 1, at(0.1, 0.0), [0.0, 5.0]);
        let f = volume(object.fuse(OperationId(3), &tool).unwrap().0);
        let c = volume(object.cut(OperationId(4), &tool).unwrap().0);
        let m = volume(object.common(OperationId(5), &tool).unwrap().0);
        let (va, vb) = (
            object.mass_properties().volume,
            tool.mass_properties().volume,
        );
        assert!((f - (va + vb - m)).abs() <= 1e-9 * f, "{k}: fuse {f}");
        assert!((c - (va - m)).abs() <= 1e-9 * va, "{k}: cut {c}");
        assert!(
            (m - 23.478130341535252).abs() <= 1e-9 * m,
            "{k}: common {m}"
        );
    }
}

/// A slab's round hole and a rod of its radius whose axis crosses the
/// hole's at right angles (exact frames, S9c.1's two ellipses): the walls
/// are tangent where the ellipses cross. The cut keeps the rod's wall's
/// outside there and evaluates (the slab less the hole, less the rod over
/// the slab's width, plus the Steinmetz common of the hole and the rod,
/// `16 r^3 / 3`); the fuse's void is the hole's two halves either side of
/// the rod, touching at those points, and the common the rod's two halves
/// either side of the hole: `Degenerate`, as a Steinmetz cut is. The fuse
/// failed validation instead (`InvalidTopology("non_manifold_vertex")`,
/// the boolean fuzz target's crash) until the assembly checked each shell
/// for touching itself at a vertex. Crossing at another angle (turned
/// frames) the pair is a tangency, in the tilted frame a near node, all
/// three operations refused; radii apart by more than the resolution
/// evaluate with the pair identities.
#[test]
fn a_hole_and_a_rod_of_its_radius_crossing_it_touch_at_two_points() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let t = Tolerance::default();
    let fr = |o: [f64; 3], n: [f64; 3], x: [f64; 3]| {
        Frame3::new(
            Point3::new(o[0], o[1], o[2]),
            Vec3::new(n[0], n[1], n[2]),
            Vec3::new(x[0], x[1], x[2]),
            t,
        )
        .unwrap()
    };
    let holed = |f: Frame3, r: f64| {
        let square = Boundary::polygon(
            [(-3.0, -3.0), (3.0, -3.0), (3.0, 3.0), (-3.0, 3.0)]
                .map(|(x, y)| Point2::new(x, y))
                .to_vec(),
            t,
        )
        .unwrap();
        let hole = Boundary::circle(Point2::new(0.0, 0.0), r, t).unwrap();
        let profile = Profile::new(square, vec![hole], t).unwrap();
        Solid::extrude_with(OperationId(1), profile, f, -2.0, 2.0)
            .unwrap()
            .0
    };
    let rod = |f: Frame3, r: f64| {
        let disc = Boundary::circle(Point2::new(0.0, 0.0), r, t).unwrap();
        let profile = Profile::new(disc, vec![], t).unwrap();
        Solid::extrude_with(OperationId(2), profile, f, -4.0, 4.0)
            .unwrap()
            .0
    };
    let volume = |out: &[Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    let refused = |r: Result<(Vec<Solid>, history::History), Error>, why: &str| match r {
        Err(Error::Degenerate(m)) => assert_eq!(m, why),
        Err(e) => panic!("{e:?}, not {why}"),
        Ok(_) => panic!("evaluated, not {why}"),
    };
    let xy = Frame3::xy();
    let side = fr([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let (a, b) = (holed(xy, 1.0), rod(side, 1.0));
    let itself = "a result touching itself at a vertex";
    let two = "solids touching at a vertex";
    refused(a.fuse(OperationId(3), &b), itself);
    refused(a.common(OperationId(5), &b), two);
    let (cut, _) = a.cut(OperationId(4), &b).unwrap();
    let want = 144.0 - 10.0 * std::f64::consts::PI + 16.0 / 3.0;
    assert_eq!(cut.len(), 1);
    assert!(
        (volume(&cut) - want).abs() <= 1e-9 * want,
        "{}",
        volume(&cut)
    );
    // The other way round: the same fuse and common; the rod less the slab
    // is its two ends past the slab and the Steinmetz solid in the hole.
    refused(b.fuse(OperationId(3), &a), itself);
    refused(b.common(OperationId(5), &a), two);
    let (rest, _) = b.cut(OperationId(4), &a).unwrap();
    let want = 2.0 * std::f64::consts::PI + 16.0 / 3.0;
    assert_eq!(rest.len(), 3);
    assert!(
        (volume(&rest) - want).abs() <= 1e-9 * want,
        "{}",
        volume(&rest)
    );
    // Turned frames: the tilted frame's near node, an oblique crossing's
    // tangency, in every operation.
    let tilted = (
        fr([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
        fr([1.0, -2.0, 0.5], [1.0, 0.0, 0.0], [0.0, 3.0, 4.0]),
    );
    let oblique = (xy, fr([0.0; 3], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]));
    for ((fa, fb), why) in [
        (
            tilted,
            "two cylinders' section within the resolution of a node",
        ),
        (oblique, "a tangency between the inputs (S9c)"),
    ] {
        let (a, b) = (holed(fa, 1.0), rod(fb, 1.0));
        refused(a.fuse(OperationId(3), &b), why);
        refused(a.cut(OperationId(4), &b), why);
        refused(a.common(OperationId(5), &b), why);
    }
    // A rod of radius `1 + 2^-10` (two rings about the hole, about `0.09`
    // apart where the walls nearly touch): every operation evaluates.
    for (fa, fb) in [(xy, side), tilted, oblique] {
        let (a, b) = (holed(fa, 1.0), rod(fb, 1.0 + 1.0 / 1024.0));
        let f = volume(&a.fuse(OperationId(3), &b).unwrap().0);
        let c = volume(&a.cut(OperationId(4), &b).unwrap().0);
        let m = volume(&a.common(OperationId(5), &b).unwrap().0);
        let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
        assert!((f - (va + vb - m)).abs() <= 1e-9 * f, "fuse {f}");
        assert!((c - (va - m)).abs() <= 1e-9 * va, "cut {c}");
    }
}

/// The boolean fuzz target's input
/// (`fuzz/regressions/boolean/replay-47992d6f….bin`): a square with a round
/// hole of radius 1 and a stadium of that radius on its side, its first
/// arc's axis crossing the hole's at its height's middle, so the walls are
/// tangent at the arc's middle. The fuse and the common (each one solid,
/// the stadium's flat end joining the common's two halves) touch
/// themselves there; the cut evaluates.
#[test]
fn a_stadium_on_its_side_through_an_equal_hole() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{
        Boundary, Error, Frame3, Point2, Point3, Profile, Segment, Solid, Tolerance, Vec3,
    };
    let t = Tolerance::default();
    let p = |x: f64, y: f64| Point2::new(x, y);
    let square = Boundary::polygon(
        vec![p(-2.5, -2.5), p(2.5, -2.5), p(2.5, 2.5), p(-2.5, 2.5)],
        t,
    )
    .unwrap();
    let hole = Boundary::circle(p(0.0, 0.0), 1.0, t).unwrap();
    let arc = |x: f64| Segment::Arc {
        center: p(x, 0.0),
        radius: 1.0,
        ccw: true,
    };
    let stadium = Boundary::path(
        vec![p(0.0, -1.0), p(1.25, -1.0), p(1.25, 1.0), p(0.0, 1.0)],
        vec![Segment::Line, arc(1.25), Segment::Line, arc(0.0)],
        t,
    )
    .unwrap();
    let h = 2.25;
    let side = Frame3::new(
        Point3::new(-h / 2.0, 0.0, h / 2.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        t,
    )
    .unwrap();
    let a = Profile::new(square, vec![hole], t).unwrap();
    let (a, _) = Solid::extrude_with(OperationId(1), a, Frame3::xy(), 0.0, h).unwrap();
    let b = Profile::new(stadium, vec![], t).unwrap();
    let (b, _) = Solid::extrude_with(OperationId(2), b, side, -1.0, h + 1.0).unwrap();
    for r in [a.fuse(OperationId(3), &b), a.common(OperationId(5), &b)] {
        match r {
            Err(Error::Degenerate(m)) => assert_eq!(m, "a result touching itself at a vertex"),
            other => panic!("{:?}", other.map(|(o, _)| o.len())),
        }
    }
    // The stadium lies within the square's slab; within the hole it holds
    // the hole's half towards its flat side over its height 2 and half the
    // Steinmetz solid of the hole and its arc, `8 / 3`.
    let pi = std::f64::consts::PI;
    let common = (2.5 + pi) * 4.25 - (pi + 8.0 / 3.0);
    let want = (25.0 - pi) * h - common;
    let (cut, _) = a.cut(OperationId(4), &b).unwrap();
    assert_eq!(cut.len(), 1);
    let v = cut[0].mass_properties().volume;
    assert!((v - want).abs() <= 1e-9 * want, "{v} for {want}");
}

/// Radii `1` and `1 + 2^-k` crossing (a rod and a rod, a slab's hole and a
/// rod, either radius the larger; the slab also standing on the plane
/// through the near nodes, its vertices there) in the exact, tilted and
/// oblique frames: the two branches of the walls' meeting about
/// `2 sqrt(2^(1 - k)) / A` apart where they nearly touch. A stored
/// meeting's height there is uncertain by about `eps L^2 / sqrt(D)`
/// (`turned::conditioned_node`); from `k` = 43 the validator's enclosures
/// of a ring's closing point or of a vertex there exceeded the resolution,
/// and results failed validation (`InvalidTopology`,
/// `enclosure_exceeds_resolution`, nearer with `uncertified_loop_winding`)
/// in each frame, up to `k` = 52. The band is a near node now in every
/// operation, and `k` = 40 evaluates with the pair identities where it
/// failed.
#[test]
fn crossing_cylinders_near_a_node_are_refused_or_valid() {
    use rusty_occt::identity::OperationId;
    use rusty_occt::{Boundary, Error, Frame3, Point2, Point3, Profile, Solid, Tolerance, Vec3};
    let t = Tolerance::default();
    let fr = |o: [f64; 3], n: [f64; 3], x: [f64; 3]| {
        Frame3::new(
            Point3::new(o[0], o[1], o[2]),
            Vec3::new(n[0], n[1], n[2]),
            Vec3::new(x[0], x[1], x[2]),
            t,
        )
        .unwrap()
    };
    let holed = |f: Frame3, r: f64, z0: f64| {
        let square = Boundary::polygon(
            [(-3.0, -3.0), (3.0, -3.0), (3.0, 3.0), (-3.0, 3.0)]
                .map(|(x, y)| Point2::new(x, y))
                .to_vec(),
            t,
        )
        .unwrap();
        let hole = Boundary::circle(Point2::new(0.0, 0.0), r, t).unwrap();
        let profile = Profile::new(square, vec![hole], t).unwrap();
        Solid::extrude_with(OperationId(1), profile, f, z0, 2.0)
            .unwrap()
            .0
    };
    let rod = |f: Frame3, r: f64, id: u64| {
        let disc = Boundary::circle(Point2::new(0.0, 0.0), r, t).unwrap();
        let profile = Profile::new(disc, vec![], t).unwrap();
        Solid::extrude_with(OperationId(id), profile, f, -4.0, 4.0)
            .unwrap()
            .0
    };
    let xy = Frame3::xy();
    let exact = (xy, fr([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]));
    let tilted = (
        fr([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
        fr([1.0, -2.0, 0.5], [1.0, 0.0, 0.0], [0.0, 3.0, 4.0]),
    );
    let oblique = (xy, fr([0.0; 3], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]));
    // The first operand (a rod, a slab from -2, or one from the nodes'
    // plane) and the rod: radii 1 and 1 + 2^-k, or swapped.
    #[derive(Clone, Copy)]
    enum First {
        Rod,
        Slab,
        Standing,
    }
    let pair = |(fa, fb): (Frame3, Frame3), first: First, swapped: bool, k: i32| {
        let e = 1.0 + 2f64.powi(-k);
        let (ra, rb) = if swapped { (e, 1.0) } else { (1.0, e) };
        let a = match first {
            First::Rod => rod(fa, ra, 1),
            First::Slab => holed(fa, ra, -2.0),
            First::Standing => holed(fa, ra, 0.0),
        };
        (a, rod(fb, rb, 2))
    };
    let node = "two cylinders' section within the resolution of a node";
    for k in 42..=52 {
        for frames in [exact, tilted, oblique] {
            for first in [First::Rod, First::Slab, First::Standing] {
                for swapped in [false, true] {
                    let (a, b) = pair(frames, first, swapped, k);
                    for r in [
                        a.fuse(OperationId(3), &b),
                        a.cut(OperationId(4), &b),
                        a.common(OperationId(5), &b),
                    ] {
                        match r {
                            Err(Error::Degenerate(m)) => assert_eq!(m, node, "k = {k}"),
                            other => panic!("k = {k}: {:?}", other.map(|(o, _)| o.len())),
                        }
                    }
                }
            }
        }
    }
    // Where results failed validation nearer: the swapped rods in exact
    // frames (rings over the thinner rod closing at its near node), the
    // oblique rods, and the standing slab's vertices at the near nodes.
    let volume = |out: &[Solid]| out.iter().map(|s| s.mass_properties().volume).sum::<f64>();
    for (frames, first, swapped) in [
        (exact, First::Rod, true),
        (oblique, First::Rod, false),
        (exact, First::Standing, true),
    ] {
        let (a, b) = pair(frames, first, swapped, 40);
        let f = volume(&a.fuse(OperationId(3), &b).unwrap().0);
        let c = volume(&a.cut(OperationId(4), &b).unwrap().0);
        let m = volume(&a.common(OperationId(5), &b).unwrap().0);
        let (va, vb) = (a.mass_properties().volume, b.mass_properties().volume);
        assert!((f - (va + vb - m)).abs() <= 1e-9 * f, "fuse {f}");
        assert!((c - (va - m)).abs() <= 1e-9 * va, "cut {c}");
    }
}
