//! STEP import (STEP-a and STEP-b of the STEP import track in
//! `REVIEW_NOTES.md`): every authored fixture of `fixtures/step` imports as
//! the independent reference (`step-expected.tsv`, from `step_reference.py`)
//! describes it, except the quarter-arc rational cylinder, which R4's C1 rule
//! refuses, and every construct outside the sub-steps is reported by name.
use rusty_occt::identity::OperationKind;
use rusty_occt::occt_brep::Rejected;
use rusty_occt::step::{self, Item, StepError, StepImport};
use rusty_occt::topology::{BodyClass, IssueKind};
use std::path::PathBuf;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../fixtures")
}

fn fixture(name: &str) -> String {
    std::fs::read_to_string(fixtures().join("step").join(format!("{name}.stp"))).unwrap()
}

fn import(text: &str) -> Result<StepImport, StepError> {
    step::import(&step::read(text.as_bytes())?)
}

fn edit(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "{from}");
    text.replacen(from, to, 1)
}

/// The kernel against the closed forms: the files' decimal data are
/// binary64, so a surface and its edges agree only to rounding.
const BOUND: f64 = 1e-12;

fn inside([lo, hi]: [f64; 2], value: f64, size: f64) {
    let slack = BOUND * size.abs().max(1.0);
    assert!(
        lo - slack <= value && value <= hi + slack,
        "{value} outside [{lo}, {hi}]"
    );
}

/// The fixture the kernel's validator refuses: its circles are rational
/// quarter arcs whose double knots are not C1 in homogeneous form (R4).
const NOT_C1: &str = "rational_cylinder";

#[test]
fn fixtures_import_as_the_reference() {
    let table = std::fs::read_to_string(fixtures().join("step-expected.tsv")).unwrap();
    let mut cases = 0;
    for line in table.lines().skip(1) {
        let w: Vec<&str> = line.split('\t').collect();
        let (name, entity, class) = (w[0], w[1].parse::<u64>().unwrap(), w[2]);
        if name == NOT_C1 {
            continue;
        }
        let counts: Vec<usize> = w[3].split(' ').map(|c| c.parse().unwrap()).collect();
        let centre: Vec<f64> = w[6].split(' ').map(|c| c.parse().unwrap()).collect();
        let area: f64 = w[5].parse().unwrap();
        let imported = import(&fixture(name)).unwrap();
        assert!(imported.unsupported.is_empty(), "{name}");
        let body = imported
            .bodies
            .iter()
            .find(|b| b.entity == entity)
            .unwrap_or_else(|| panic!("{name} #{entity}"));
        let t = body
            .result
            .as_ref()
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        // The files' 1e-7 mm, in metres and inches too (one rounding).
        assert!((body.tolerance.linear() - 1e-7).abs() <= 1e-22, "{name}");
        let c = t.occt_counts();
        // A sheet of one face is synthesized as a free face (S6), where
        // OCCT's reader keeps the surface model's shell.
        let mut counts = counts;
        if class == "sheet" && counts[3] == 1 {
            counts[4] = 0;
        }
        assert_eq!(
            [c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids],
            counts.as_slice(),
            "{name}"
        );
        // Every slot is an External derivation, as `.brep` import gives.
        assert!(t
            .ids()
            .all(|(id, _)| t.derivation(id).unwrap().kind == OperationKind::External));
        let size = centre.iter().fold(area.sqrt(), |m, x| m.max(x.abs()));
        if class == "solid" {
            assert_eq!(
                (body.item, t.class()),
                (Item::Solid, BodyClass::Solid),
                "{name}"
            );
            let m = t.mass_enclosure().unwrap();
            inside(m.volume, w[4].parse().unwrap(), w[4].parse().unwrap());
            inside(m.surface_area, area, area);
            for (i, x) in centre.iter().enumerate() {
                inside(m.centroid[i], *x, size);
            }
        } else {
            assert_eq!(
                (body.item, t.class()),
                (Item::Shell, BodyClass::Sheet),
                "{name}"
            );
            let m = t.measure_enclosure().unwrap();
            inside(m.measure, area, area);
            for (i, x) in centre.iter().enumerate() {
                inside(m.centre[i], *x, size);
            }
        }
        cases += 1;
    }
    assert_eq!(cases, 29);
}

#[test]
fn import_is_deterministic_and_cavities_are_regions() {
    for name in [
        "box_void",
        "sphere",
        "elbow",
        "syntax",
        "two_solids",
        "cylinder_oblique",
        "bspline_prism",
        "bspline_trimmed",
        "rational_cylinder",
    ] {
        let text = fixture(name);
        assert_eq!(import(&text).unwrap(), import(&text).unwrap(), "{name}");
    }
    let imported = import(&fixture("box_void")).unwrap();
    let t = imported.bodies[0].result.as_ref().unwrap();
    // The infinite void, the solid and the cavity.
    assert_eq!(t.regions().len(), 3);
    assert_eq!(t.shells().len(), 4);
}

/// A torus band whose `+v` ring starts a rounding below the latitude seam
/// the `.brep` writer puts at `v0` (the `step` fuzz target's first finding,
/// `fuzz/regressions/step`): written and read back with its counts.
#[test]
fn a_torus_band_round_trips_through_brep() {
    for meridian in ["(0.,1.,0.)", "(0.,-1000.0000010000001,0.)"] {
        let text = edit(
            &fixture("elbow"),
            "#31=DIRECTION('',(0.,1.,0.));",
            &format!("#31=DIRECTION('',{meridian});"),
        );
        let imported = import(&text).unwrap();
        let body = &imported.bodies[0];
        let t = body.result.as_ref().unwrap();
        let written = rusty_occt::occt_brep::write(t, body.tolerance.linear()).unwrap();
        let back = rusty_occt::occt_brep::import(&rusty_occt::occt_brep::read(&written).unwrap());
        let again = back.solids[0].result.as_ref().unwrap();
        assert_eq!(again.occt_counts(), t.occt_counts(), "{meridian}");
    }
}

#[test]
fn schemas_of_ap203_ap214_and_ap242_only() {
    let text = fixture("box");
    let schema = "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));";
    assert!(import(&edit(
        &text,
        schema,
        "FILE_SCHEMA(('AUTOMOTIVE_DESIGN_CC2 { 1 2 10303 214 -1 1 5 4 }'));"
    ))
    .is_ok());
    assert_eq!(
        import(&edit(
            &text,
            schema,
            "FILE_SCHEMA(('STRUCTURAL_ANALYSIS_DESIGN'));"
        )),
        Err(StepError::Schema("STRUCTURAL_ANALYSIS_DESIGN".into()))
    );
}

fn rejected(text: &str) -> (Vec<&'static str>, StepImport) {
    let imported = import(text).unwrap();
    let names = match &imported.bodies[0].result {
        Err(Rejected::Unsupported(names)) => names.clone(),
        other => panic!("{other:?}"),
    };
    (names, imported)
}

#[test]
fn unsupported_constructs_are_named_and_counted() {
    let cylinder = fixture("cylinder");
    // An ellipse in the cap's plane (5 by 4 about the axis) is no section
    // of the wall: no exact pcurve there (STEP-b).
    let (names, imported) = rejected(&edit(
        &cylinder,
        "#24=CIRCLE('',#23,5.);",
        "#24=ELLIPSE('',#23,5.,4.);",
    ));
    assert_eq!(names, ["PCurveNotDerived"]);
    assert_eq!(imported.unsupported.get("PCurveNotDerived"), Some(&1));
    // An oblique line on a cylinder: no straight pcurve.
    let (names, _) = rejected(&edit(
        &cylinder,
        "#53=DIRECTION('',(0.,0.,12.));",
        "#53=DIRECTION('',(0.,1.,12.));",
    ));
    assert_eq!(names, ["PCurveNotDerived"]);
    // Two ring loops on the wall and no seam.
    let two = edit(
        &edit(
            &cylinder,
            "#61=EDGE_LOOP('',(#52,#58,#59,#60));",
            "#61=EDGE_LOOP('',(#52));\n#162=EDGE_LOOP('',(#59));\n#163=FACE_BOUND('',#162,.T.);",
        ),
        "#68=ADVANCED_FACE('',(#62),#67,.T.);",
        "#68=ADVANCED_FACE('',(#62,#163),#67,.T.);",
    );
    assert_eq!(rejected(&two).0, ["PeriodicFaceWithoutSeam"]);
    let vertex_loop = edit(
        &cylinder,
        "#34=ADVANCED_FACE('',(#28),#33,.T.);",
        "#34=ADVANCED_FACE('',(#28,#164),#33,.T.);\n#165=VERTEX_LOOP('',#19);\n#164=FACE_BOUND('',#165,.T.);",
    );
    assert_eq!(rejected(&vertex_loop).0, ["VERTEX_LOOP"]);
    let units = edit(
        &cylinder,
        "SI_UNIT(.MILLI.,.METRE.)",
        "SI_UNIT(.EXA.,.METRE.)",
    );
    assert_eq!(rejected(&units).0, ["UnitsUndetermined"]);
}

#[test]
fn placements_are_reported_and_bodies_kept() {
    let text = edit(
        &fixture("cylinder"),
        "ENDSEC;\nEND-ISO-10303-21;",
        "#200=MAPPED_ITEM('',#201,#17);\n#201=REPRESENTATION_MAP(#17,#71);\nENDSEC;\nEND-ISO-10303-21;",
    );
    let imported = import(&text).unwrap();
    assert_eq!(imported.unsupported.get("AssemblyPlacement"), Some(&1));
    assert!(imported.bodies[0].result.is_ok());
}

#[test]
fn the_uncertainty_is_each_entitys_claim() {
    let text = fixture("cylinder");
    let with = |value: &str| {
        let imported = import(&edit(
            &text,
            "LENGTH_MEASURE(1.E-07)",
            &format!("LENGTH_MEASURE({value})"),
        ))
        .unwrap();
        imported.bodies[0].tolerance.linear()
    };
    assert_eq!(with("1.E-03"), 1e-3);
    // Floored at OCCT's confusion and capped at read.maxprecision.val.
    assert_eq!(with("1.E-12"), 1e-7);
    assert_eq!(with("5."), 1.0);
    // In metres, scaled to millimetres.
    let box_metre = import(&fixture("box_metre")).unwrap();
    assert_eq!(box_metre.bodies[0].tolerance.linear(), 1e-10 * 1000.0);
}

#[test]
fn inaccurate_data_is_an_import_failure_with_its_issues() {
    let text = edit(
        &fixture("cylinder"),
        "#35=CARTESIAN_POINT('',(5.,0.,12.));",
        "#35=CARTESIAN_POINT('',(5.,0.,12.001));",
    );
    let imported = import(&text).unwrap();
    let Err(Rejected::Invalid { issues, .. }) = &imported.bodies[0].result else {
        panic!("{:?}", imported.bodies[0].result);
    };
    assert!(issues.iter().any(|i| i.kind == IssueKind::VertexOffCurve));
}

/// The quarter-arc form of a rational circle (knots of multiplicity 2 at
/// the quarter points, weights `sqrt(2) / 2`) is C1 as a rational curve
/// but not in homogeneous form, which R4 certifies: the rational cylinder
/// imports as cells the validator refuses, on its four arcs, the cap
/// pcurves derived from them and its two walls, and nothing else.
#[test]
fn quarter_arc_rational_circles_are_not_c1() {
    let imported = import(&fixture(NOT_C1)).unwrap();
    assert!(imported.unsupported.is_empty());
    let Err(Rejected::Invalid { issues, .. }) = &imported.bodies[0].result else {
        panic!("{:?}", imported.bodies[0].result);
    };
    let kinds: std::collections::BTreeSet<_> = issues.iter().map(|i| i.kind).collect();
    assert_eq!(
        kinds.into_iter().collect::<Vec<_>>(),
        [
            IssueKind::EdgeNotC1,
            IssueKind::PcurveNotC1,
            IssueKind::FaceNotC1
        ]
    );
    assert_eq!(
        issues
            .iter()
            .filter(|i| i.kind == IssueKind::EdgeNotC1)
            .count(),
        4
    );
}

/// Every simple B-spline curve or surface of a file rewritten as the complex
/// rational instance OCCT writes, all its weights `weight`.
fn rational(text: &str, weight: &str) -> String {
    let mut out = String::new();
    for line in text.lines() {
        let curve = line.split_once("=B_SPLINE_CURVE_WITH_KNOTS('',");
        let surface = line.split_once("=B_SPLINE_SURFACE_WITH_KNOTS('',");
        let changed = if let Some((head, rest)) = curve {
            let body = rest.strip_suffix(",.UNSPECIFIED.);").unwrap();
            let (points, knots) = body.split_once(",.UNSPECIFIED.,.F.,.F.,").unwrap();
            let weights = vec![weight; points.matches('#').count()].join(",");
            format!(
                "{head}=( BOUNDED_CURVE() B_SPLINE_CURVE({points},.UNSPECIFIED.,.F.,.F.) \
                 B_SPLINE_CURVE_WITH_KNOTS({knots},.UNSPECIFIED.) CURVE() \
                 GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE(({weights})) \
                 REPRESENTATION_ITEM('') );"
            )
        } else if let Some((head, rest)) = surface {
            let body = rest.strip_suffix(",.UNSPECIFIED.);").unwrap();
            let (grid, knots) = body.split_once(",.UNSPECIFIED.,.F.,.F.,.F.,").unwrap();
            let rows = grid.matches("(#").count();
            let row = vec![weight; grid.matches('#').count() / rows].join(",");
            let weights = vec![format!("({row})"); rows].join(",");
            format!(
                "{head}=( BOUNDED_SURFACE() B_SPLINE_SURFACE({grid},.UNSPECIFIED.,.F.,.F.,.F.) \
                 B_SPLINE_SURFACE_WITH_KNOTS({knots},.UNSPECIFIED.) \
                 GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({weights})) \
                 REPRESENTATION_ITEM('') SURFACE() );"
            )
        } else {
            line.to_string()
        };
        out.push_str(&changed);
        out.push('\n');
    }
    out
}

/// The plate's and the prism's splines as complex rational instances with
/// equal weights: the same curves and surfaces, so the same bodies.
#[test]
fn rational_complex_instances_are_their_splines() {
    for name in ["bspline_plate", "bspline_prism"] {
        let plain = import(&fixture(name)).unwrap();
        let text = rational(&fixture(name), "2.");
        assert!(text.contains("RATIONAL_B_SPLINE_CURVE((2.,2.,2.,2.,2.))"));
        let weighted = import(&text).unwrap();
        let (a, b) = (
            plain.bodies[0].result.as_ref().unwrap(),
            weighted.bodies[0].result.as_ref().unwrap(),
        );
        assert_eq!(a.occt_counts(), b.occt_counts(), "{name}");
        if name == "bspline_prism" {
            assert!(text.contains("RATIONAL_B_SPLINE_SURFACE(((2.,2.),"));
            let (m, n) = (a.mass_enclosure().unwrap(), b.mass_enclosure().unwrap());
            inside(n.volume, 264.0, 264.0);
            assert!((m.surface_area[0] - n.surface_area[0]).abs() <= 1e-9);
        } else {
            inside(b.measure_enclosure().unwrap().measure, 66.0, 66.0);
        }
    }
}

/// A rational quadratic arc of a third of a circle (weights 1, 1/2, 1)
/// bounding the plate over its chord: the circular segment's area and
/// centroid in closed form.
#[test]
fn a_rational_arc_bounds_a_plate() {
    let text = edit(
        &edit(
            &fixture("bspline_plate"),
            "#33=B_SPLINE_CURVE_WITH_KNOTS('',3,(#28,#29,#30,#31,#32),.UNSPECIFIED.,.F.,.F.,(4,1,4),(0.,1.,2.),.UNSPECIFIED.);",
            "#33=( BOUNDED_CURVE() B_SPLINE_CURVE(2,(#28,#30,#32),.UNSPECIFIED.,.F.,.F.) \
             B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.) CURVE() \
             GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((1.,0.5,1.)) \
             REPRESENTATION_ITEM('') );",
        ),
        "#30=CARTESIAN_POINT('',(5.,10.,0.));",
        "#30=CARTESIAN_POINT('',(5.,8.660254037844386,0.));",
    );
    let imported = import(&text).unwrap();
    let t = imported.bodies[0].result.as_ref().unwrap();
    // Chord 10, half-angle pi / 3: radius 10 / sqrt 3, the centre r / 2
    // below the chord.
    let (theta, r) = (std::f64::consts::FRAC_PI_3, 10.0 / 3f64.sqrt());
    let area = r * r * (theta - theta.sin() * theta.cos());
    let y = -r * theta.cos()
        + 4.0 * r * theta.sin().powi(3) / (3.0 * (2.0 * theta - (2.0 * theta).sin()));
    let m = t.measure_enclosure().unwrap();
    // The middle pole's decimal is 5 sqrt 3 to 16 digits.
    let near = |[lo, hi]: [f64; 2], x: f64| lo - 1e-13 <= x && x <= hi + 1e-13;
    assert!(near(m.measure, area), "{:?} {area}", m.measure);
    assert!(
        near(m.centre[0], 5.0) && near(m.centre[1], y),
        "{:?} {y}",
        m.centre
    );
}

/// STEP-b's refusals, each by name: a knotless B-spline form, an edge whose
/// vertices run against its curve, an edge on a spline surface without the
/// face's pcurve, a B-spline pcurve running against its edge.
#[test]
fn step_b_refusals_are_named() {
    let plate = fixture("bspline_plate");
    let (names, _) = rejected(&edit(
        &plate,
        "#33=B_SPLINE_CURVE_WITH_KNOTS('',3,(#28,#29,#30,#31,#32),.UNSPECIFIED.,.F.,.F.,(4,1,4),(0.,1.,2.),.UNSPECIFIED.);",
        "#33=BEZIER_CURVE('',4,(#28,#29,#30,#31,#32),.UNSPECIFIED.,.F.,.F.);",
    ));
    assert_eq!(names, ["BEZIER_CURVE"]);
    let (names, _) = rejected(&edit(
        &plate,
        "#34=EDGE_CURVE('',#21,#19,#33,.T.);",
        "#34=EDGE_CURVE('',#19,#21,#33,.T.);",
    ));
    assert_eq!(names, ["EdgeRangeInverted"]);
    let patch = fixture("bspline_patch");
    // The first edge's pcurve names another surface.
    let (names, _) = rejected(&edit(
        &patch,
        "#50=PCURVE('',#33,#49);",
        "#50=PCURVE('',#44,#49);",
    ));
    assert_eq!(names, ["PCurveNotDerived"]);
    let (names, _) = rejected(&edit(
        &patch,
        "#62=B_SPLINE_CURVE_WITH_KNOTS('',1,(#60,#61),",
        "#62=B_SPLINE_CURVE_WITH_KNOTS('',1,(#61,#60),",
    ));
    assert_eq!(names, ["PCurveAgainstEdge"]);
}

/// A pcurve that is not its edge's image at matching fractions is an import
/// failure with its issues, never reparameterised: the trim's middle pole
/// moved along `u` keeps the pcurve's ends and bends its parameter.
#[test]
fn a_pcurve_off_its_edge_is_an_import_failure() {
    let text = edit(
        &fixture("bspline_trimmed"),
        "#84=CARTESIAN_POINT('',(0.5,0.75));",
        "#84=CARTESIAN_POINT('',(0.25,0.75));",
    );
    let imported = import(&text).unwrap();
    let Err(Rejected::Invalid { issues, .. }) = &imported.bodies[0].result else {
        panic!("{:?}", imported.bodies[0].result);
    };
    assert!(issues.iter().any(|i| i.kind == IssueKind::PcurveOffEdge));
}

/// The spline bodies write as `.brep` and read back with their counts: a
/// line pcurve on a spline surface, whose length is not its edge's span,
/// runs over the span (the `step` fuzz target's first STEP-b finding,
/// `fuzz/regressions/step`).
#[test]
fn spline_bodies_round_trip_through_brep() {
    for name in [
        "bspline_plate",
        "bspline_prism",
        "bspline_patch",
        "bspline_trimmed",
    ] {
        let imported = import(&fixture(name)).unwrap();
        let body = &imported.bodies[0];
        let t = body.result.as_ref().unwrap();
        let written = rusty_occt::occt_brep::write(t, body.tolerance.linear()).unwrap();
        let back = rusty_occt::occt_brep::import(&rusty_occt::occt_brep::read(&written).unwrap());
        let again = back
            .solids
            .iter()
            .map(|s| &s.result)
            .chain(back.free.iter().map(|f| &f.result))
            .next()
            .unwrap()
            .as_ref()
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));
        assert_eq!(again.occt_counts(), t.occt_counts(), "{name}");
    }
}
