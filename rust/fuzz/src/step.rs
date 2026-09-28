//! STEP import (the STEP import track of REVIEW_NOTES.md). The first byte
//! picks the input: the rest as a whole file, the rest as the data section
//! of a file with a valid header, or one of the authored fixtures
//! (`fixtures/step`) with up to six mutations of its instances (numbers
//! moved, enumerations flipped, references redirected, entity names and
//! tokens replaced, instances deleted, duplicated or swapped, the text
//! truncated). Reading never panics: malformed text is a typed `StepError`.
//! Importing never panics, gives the same result twice, and every body is
//! a topology that validates, and writes and reads back as `.brep` with its
//! counts when the writer expresses it (the `.brep` reader does not read
//! ellipse records yet, so a body written with one is not read back), or is
//! rejected with its unsupported names or its validation issues. An
//! unmutated fixture imports every body, except the quarter-arc rational
//! cylinder, which the validator refuses (R4). STEP-b's fixtures bring
//! ellipses, B-spline curves and surfaces, rational complex instances and
//! the file's pcurves into the mutations.
use libfuzzer_sys::arbitrary::{Result, Unstructured};
use rusty_occt::occt_brep::{self, Rejected};
use rusty_occt::step::{self, StepImport};

const FIXTURES: [&str; 29] = [
    include_str!("../../fixtures/step/box.stp"),
    include_str!("../../fixtures/step/box_ap203.stp"),
    include_str!("../../fixtures/step/box_ap242.stp"),
    include_str!("../../fixtures/step/box_flipped.stp"),
    include_str!("../../fixtures/step/box_inch.stp"),
    include_str!("../../fixtures/step/box_metre.stp"),
    include_str!("../../fixtures/step/box_void.stp"),
    include_str!("../../fixtures/step/cone.stp"),
    include_str!("../../fixtures/step/cylinder.stp"),
    include_str!("../../fixtures/step/cylinder_tilted.stp"),
    include_str!("../../fixtures/step/elbow.stp"),
    include_str!("../../fixtures/step/frustum.stp"),
    include_str!("../../fixtures/step/frustum_degree.stp"),
    include_str!("../../fixtures/step/half_cylinder.stp"),
    include_str!("../../fixtures/step/hemisphere.stp"),
    include_str!("../../fixtures/step/l_prism.stp"),
    include_str!("../../fixtures/step/open_box.stp"),
    include_str!("../../fixtures/step/plate_hole.stp"),
    include_str!("../../fixtures/step/sphere.stp"),
    include_str!("../../fixtures/step/syntax.stp"),
    include_str!("../../fixtures/step/torus.stp"),
    include_str!("../../fixtures/step/two_solids.stp"),
    include_str!("../../fixtures/step/ellipse_sheet.stp"),
    include_str!("../../fixtures/step/cylinder_oblique.stp"),
    include_str!("../../fixtures/step/bspline_plate.stp"),
    include_str!("../../fixtures/step/bspline_prism.stp"),
    include_str!("../../fixtures/step/bspline_patch.stp"),
    include_str!("../../fixtures/step/bspline_trimmed.stp"),
    include_str!("../../fixtures/step/rational_cylinder.stp"),
];

/// The fixture whose unmutated body the validator refuses: rational
/// quarter arcs are not C1 in homogeneous form (R4).
const NOT_C1: usize = 28;

/// Tokens a mutation may put in place of another: other entities, loop and
/// shell kinds, flags, references, special numbers and punctuation.
const WORDS: [&str; 41] = [
    "ELLIPSE",
    "B_SPLINE_CURVE_WITH_KNOTS",
    "B_SPLINE_SURFACE_WITH_KNOTS",
    "RATIONAL_B_SPLINE_CURVE",
    "RATIONAL_B_SPLINE_SURFACE",
    "BEZIER_CURVE",
    "SURFACE_CURVE",
    "PCURVE",
    "DEFINITIONAL_REPRESENTATION",
    "VERTEX_LOOP",
    "POLY_LOOP",
    "FACE_SURFACE",
    "ORIENTED_FACE",
    "ORIENTED_CLOSED_SHELL",
    "OPEN_SHELL",
    "SEAM_CURVE",
    "TRIMMED_CURVE",
    "CYLINDRICAL_SURFACE",
    "CONICAL_SURFACE",
    "SPHERICAL_SURFACE",
    "TOROIDAL_SURFACE",
    "PLANE",
    "CIRCLE",
    "LINE",
    "MAPPED_ITEM",
    ".T.",
    ".F.",
    ".U.",
    "$",
    "*",
    "#1",
    "#99999",
    "0.",
    "-1.",
    "1.E308",
    "1.E-300",
    "'",
    "(",
    ")",
    ",",
    ".KILO.",
];

const HEADER: &str =
    "ISO-10303-21;HEADER;FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));ENDSEC;DATA;\n#9000=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n";
const FOOTER: &str = "\nENDSEC;END-ISO-10303-21;\n";

/// Punctuation are tokens of their own; everything between them is one.
fn tokens(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut word = String::new();
    for c in line.chars() {
        if "(),;=".contains(c) {
            if !word.is_empty() {
                out.push(std::mem::take(&mut word));
            }
            out.push(c.to_string());
        } else {
            word.push(c);
        }
    }
    if !word.is_empty() {
        out.push(word);
    }
    out
}

fn mutate(u: &mut Unstructured, text: &str) -> Result<String> {
    let mut lines: Vec<Vec<String>> = text.lines().map(tokens).collect();
    for _ in 0..u.int_in_range(1..=6)? {
        if lines.is_empty() {
            break;
        }
        let at = u.choose_index(lines.len())?;
        let pick =
            |u: &mut Unstructured, line: &[String], f: fn(&str) -> bool| -> Result<Option<usize>> {
                let found: Vec<usize> = (0..line.len()).filter(|k| f(&line[*k])).collect();
                Ok(if found.is_empty() {
                    None
                } else {
                    Some(found[u.choose_index(found.len())?])
                })
            };
        match u.int_in_range(0..=8)? {
            // A number moves, flips sign, grows or vanishes.
            0 | 1 => {
                if let Some(k) = pick(u, &lines[at], |w| w.trim().parse::<f64>().is_ok())? {
                    let x: f64 = lines[at][k].trim().parse().unwrap_or(0.0);
                    let y = match u.int_in_range(0..=3)? {
                        0 => x * (1.0 + 1e-9),
                        1 => -x,
                        2 => x * 1000.0,
                        _ => 0.0,
                    };
                    lines[at][k] = format!("{y:?}").replace("e", "E");
                }
            }
            // An enumeration flips.
            2 => {
                if let Some(k) = pick(u, &lines[at], |w| w == ".T." || w == ".F.")? {
                    lines[at][k] = if lines[at][k] == ".T." { ".F." } else { ".T." }.into();
                }
            }
            // A reference points elsewhere, maybe nowhere.
            3 => {
                if let Some(k) = pick(u, &lines[at], |w| w.starts_with('#'))? {
                    lines[at][k] = format!("#{}", u.int_in_range(1..=lines.len() as u32 + 2)?);
                }
            }
            4 => {
                lines.remove(at);
            }
            5 => {
                let copy = lines[at].clone();
                lines.insert(at, copy);
            }
            6 => {
                let other = u.choose_index(lines.len())?;
                lines.swap(at, other);
            }
            7 => {
                if !lines[at].is_empty() {
                    let k = u.choose_index(lines[at].len())?;
                    lines[at][k] = u.choose(&WORDS)?.to_string();
                }
            }
            _ => lines.truncate(at),
        }
    }
    Ok(lines
        .iter()
        .map(|l| l.concat())
        .collect::<Vec<_>>()
        .join("\n"))
}

/// Whatever reads and imports keeps the contract of the module doc.
fn check_text(bytes: &[u8]) -> Option<StepImport> {
    let x = step::read(bytes).ok()?;
    let imported = step::import(&x).ok()?;
    assert_eq!(
        step::import(&x).as_ref().ok(),
        Some(&imported),
        "import is deterministic"
    );
    for body in &imported.bodies {
        match &body.result {
            Ok(t) => {
                assert!(
                    t.check(body.tolerance).is_empty(),
                    "an imported body validates"
                );
                if let Ok(text) = occt_brep::write(t, body.tolerance.linear()) {
                    let doc = occt_brep::read(&text).expect("the writer's text reads");
                    let back = occt_brep::import(&doc);
                    // The reader names ellipse records unsupported (STEP-b).
                    if back.unsupported.contains_key("Ellipse")
                        || back.unsupported.contains_key("Ellipse2d")
                    {
                        continue;
                    }
                    let again = back
                        .solids
                        .iter()
                        .map(|s| &s.result)
                        .chain(back.free.iter().map(|f| &f.result))
                        .next()
                        .expect("a written body")
                        .as_ref()
                        .expect("a written body imports");
                    assert_eq!(again.occt_counts(), t.occt_counts());
                }
            }
            Err(Rejected::Unsupported(names)) => assert!(!names.is_empty()),
            Err(Rejected::Invalid { issues, .. }) => assert!(!issues.is_empty()),
        }
    }
    Some(imported)
}

pub fn check_step(data: &[u8]) {
    match data.first() {
        None => {}
        Some(0) => {
            check_text(&data[1..]);
        }
        Some(1) => {
            let mut text = HEADER.as_bytes().to_vec();
            text.extend_from_slice(&data[1..]);
            text.extend_from_slice(FOOTER.as_bytes());
            check_text(&text);
        }
        Some(_) => {
            let mut u = Unstructured::new(&data[1..]);
            let Ok(pick) = u.choose_index(FIXTURES.len()) else {
                return;
            };
            let base = FIXTURES[pick];
            if u.ratio(1, 8).unwrap_or(true) {
                let imported = check_text(base.as_bytes()).expect("a fixture reads and imports");
                assert!(
                    imported.bodies.iter().all(|b| if pick == NOT_C1 {
                        matches!(b.result, Err(Rejected::Invalid { .. }))
                    } else {
                        b.result.is_ok()
                    }),
                    "an unmutated fixture imports every body"
                );
                return;
            }
            if let Ok(text) = mutate(&mut u, base) {
                check_text(text.as_bytes());
            }
        }
    }
}
