//! STEP import (STEP-a of the STEP import track in `REVIEW_NOTES.md`): every
//! authored fixture of `fixtures/step` imports as the independent reference
//! (`step-expected.tsv`, from `step_reference.py`) describes it, and every
//! construct outside the sub-step is reported by name.
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

#[test]
fn fixtures_import_as_the_reference() {
    let table = std::fs::read_to_string(fixtures().join("step-expected.tsv")).unwrap();
    let mut cases = 0;
    for line in table.lines().skip(1) {
        let w: Vec<&str> = line.split('\t').collect();
        let (name, entity, class) = (w[0], w[1].parse::<u64>().unwrap(), w[2]);
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
    assert_eq!(cases, 23);
}

#[test]
fn import_is_deterministic_and_cavities_are_regions() {
    for name in ["box_void", "sphere", "elbow", "syntax", "two_solids"] {
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
    let (names, imported) = rejected(&edit(
        &cylinder,
        "#24=CIRCLE('',#23,5.);",
        "#24=ELLIPSE('',#23,5.,4.);",
    ));
    assert_eq!(names, ["ELLIPSE"]);
    assert_eq!(imported.unsupported.get("ELLIPSE"), Some(&1));
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
