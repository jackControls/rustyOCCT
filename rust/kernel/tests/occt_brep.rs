//! OCCT `.brep` interop (T2): the reader's typed errors, the writer and the
//! round trip of every prism fixture through the cell model.
#[path = "support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use rusty_occt::occt_brep::{import, read, write, BrepError};

#[test]
fn every_prism_round_trips_to_the_same_cells_and_text() {
    let specs = identity_protocol::cases(include_str!("../../fixtures/identity-cases.txt"));
    let mut checked = 0;
    for spec in &specs {
        let solid = identity_protocol::build(spec);
        let tol = solid.resolution().linear();
        let text = write(solid.topology(), tol).unwrap();
        let doc = read(&text).unwrap();
        let im = import(&doc);
        assert!(
            im.unsupported.is_empty(),
            "{}: {:?}",
            spec.name,
            im.unsupported
        );
        assert_eq!(im.solids.len(), 1);
        let back = im.solids[0]
            .result
            .as_ref()
            .unwrap_or_else(|e| panic!("{}: {e:?}", spec.name));
        assert_eq!(
            back.occt_counts(),
            solid.topology().occt_counts(),
            "{}",
            spec.name
        );
        // Vertices come back bit for bit (shortest round-trip text).
        let mut a: Vec<[u64; 3]> = solid
            .topology()
            .vertices()
            .iter()
            .map(|v| v.position.to_array().map(f64::to_bits))
            .collect();
        let mut b: Vec<[u64; 3]> = back
            .vertices()
            .iter()
            .map(|v| v.position.to_array().map(f64::to_bits))
            .collect();
        a.sort();
        b.sort();
        assert_eq!(a, b, "{}", spec.name);
        // Import numbers entities in discovery order and re-derives plane
        // normals and pcurve directions, so rewriting is not a text fixed
        // point; a second round trip keeps the cells, counts and vertices.
        let second = write(back, tol).unwrap();
        let again = import(&read(&second).unwrap());
        let again = again.solids[0].result.as_ref().unwrap();
        assert_eq!(again.occt_counts(), back.occt_counts(), "{}", spec.name);
        assert!(again.vertices().iter().zip(back.vertices()).all(|(x, y)| x
            .position
            .to_array()
            .map(f64::to_bits)
            == y.position.to_array().map(f64::to_bits)));
        checked += 1;
    }
    assert_eq!(checked, specs.len());
}

#[test]
fn malformed_documents_are_typed_errors_with_their_line() {
    let good = write(
        identity_protocol::build(
            &identity_protocol::cases(include_str!("../../fixtures/identity-cases.txt"))[0],
        )
        .topology(),
        1e-7,
    )
    .unwrap();
    assert!(read(&good).is_ok());
    for (bad, line) in [
        (
            good.replacen("CASCADE Topology V1", "CASCADE Topology V9", 1),
            3,
        ),
        (good.replacen("Curves", "Curvez", 1), 0),
        (good.replacen("TShapes", "TShapes 99999999999", 1), 0),
    ] {
        match read(&bad) {
            Err(BrepError::Syntax { line: l, .. }) => assert!(line == 0 || l == line, "{l}"),
            other => panic!("{other:?}"),
        }
    }
    // Every truncation fails cleanly.
    for cut in (0..good.len()).step_by(97) {
        let _ = read(&good[..cut]);
    }
    // An edge's 3D curve representation `1 curve location first last`
    // referring beyond the tables (the first `brep_io` fuzz finding).
    let edge = good.find("Ed\n").unwrap();
    let rep = edge + good[edge..].find("\n1 ").unwrap() + 1;
    let words: Vec<&str> = good[rep..].split(' ').take(3).collect();
    for (curve, location) in [("999", words[2]), (words[1], "1")] {
        let bad = format!(
            "{}1 {curve} {location}{}",
            &good[..rep],
            &good[rep + words.join(" ").len()..]
        );
        match read(&bad) {
            Err(BrepError::Reference { .. }) => {}
            other => panic!("{other:?}"),
        }
    }
}

const GEOMETRY: [&str; 21] = [
    "BSplineCurve",
    "BSplineCurve2d",
    "BSplineSurface",
    "BezierCurve",
    "BezierCurve2d",
    "BezierSurface",
    "Ellipse",
    "Ellipse2d",
    "Hyperbola",
    "Hyperbola2d",
    "NonRingToroidalSurface",
    "OffsetCurve",
    "OffsetCurve2d",
    "OffsetSurface",
    "Parabola",
    "Parabola2d",
    "RectangularTrimmedSurface",
    "SurfaceOfLinearExtrusion",
    "SurfaceOfRevolution",
    "TrimmedCurve",
    "TrimmedCurve2d",
];

/// The data/occ corpus against brep_io_reference.py: the same unsupported
/// geometry records, the same solids, a cell topology exactly for the
/// solids whose structure is representable, and synthesized counts equal to
/// OCCT's distinct subshapes of the original seamed solid.
/// A sub-shape reference made internal or external (the cell model has no
/// such orientation) imports without a panic, and the solid holding it is
/// never certified: OCCT's `bug21246.brep` has internal edges and faces.
#[test]
fn internal_and_external_references_are_unsupported() {
    let text = include_str!("../../../data/occ/wedge_ok.brep");
    let lines: Vec<&str> = text.lines().collect();
    let shapes = lines.iter().position(|l| l.starts_with("TShapes")).unwrap();
    let mut mutated = 0;
    for (k, line) in lines.iter().enumerate().skip(shapes + 1) {
        let words: Vec<&str> = line.split_whitespace().collect();
        if words.last() != Some(&"*") || words.len() < 3 {
            continue;
        }
        for (w, word) in words.iter().enumerate().step_by(2) {
            if !word.starts_with(['+', '-']) {
                continue;
            }
            for mark in ["i", "e"] {
                let mut changed = words.clone();
                let replaced = format!("{mark}{}", &word[1..]);
                changed[w] = &replaced;
                let mut all = lines.clone();
                let joined = changed.join(" ");
                all[k] = &joined;
                let im = import(&read(&(all.join("\n") + "\n")).unwrap());
                assert!(
                    im.solids.iter().all(|s| s.result.is_err()),
                    "line {} word {w} marked {mark}",
                    k + 1
                );
                mutated += 1;
            }
        }
    }
    assert!(mutated >= 40, "{mutated}");
}

#[test]
fn corpus_matches_the_independent_reader() {
    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/occ/");
    let expected = include_str!("../../fixtures/brep-io-expected.tsv");
    let mut files: Vec<(&str, &str, Vec<Vec<&str>>)> = Vec::new();
    for line in expected.lines() {
        let row: Vec<&str> = line.split('\t').collect();
        match row[0] {
            "file" => files.push((row[1], row[2], Vec::new())),
            "solid" => files.last_mut().unwrap().2.push(row[2..].to_vec()),
            other => panic!("row kind {other}"),
        }
    }
    let mut imported = 0;
    let mut not_certified: Vec<(String, String, Vec<&'static str>)> = Vec::new();
    for (name, geometry, solids) in &files {
        let text = std::fs::read_to_string(format!("{root}{name}")).unwrap();
        let im = import(&read(&text).unwrap_or_else(|e| panic!("{name}: {e}")));
        let names: Vec<String> = im
            .unsupported
            .iter()
            .filter(|(k, _)| GEOMETRY.contains(k))
            .map(|(k, v)| format!("{k}={v}"))
            .collect();
        let names = if names.is_empty() {
            "-".to_string()
        } else {
            names.join(",")
        };
        assert_eq!(&names, geometry, "{name}");
        assert_eq!(im.solids.len(), solids.len(), "{name}");
        for (solid, row) in im.solids.iter().zip(solids) {
            assert_eq!(solid.record.to_string(), row[0], "{name}");
            match (&solid.result, row[1]) {
                (Ok(t), "representable") => {
                    let c = t.occt_counts();
                    let got = [c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids]
                        .map(|n| n.to_string());
                    assert_eq!(got.as_slice(), &row[2..], "{name} {}", row[0]);
                    // Written back, it reads as one certified solid with the
                    // same counts.
                    let text = write(t, solid.tolerance.linear())
                        .unwrap_or_else(|e| panic!("{name} {}: {e}", row[0]));
                    let again = import(&read(&text).unwrap());
                    assert_eq!(again.solids.len(), 1, "{name} {}", row[0]);
                    let back = again.solids[0]
                        .result
                        .as_ref()
                        .unwrap_or_else(|e| panic!("{name} {} written back: {e:?}", row[0]));
                    assert_eq!(back.occt_counts(), c, "{name} {}", row[0]);
                    imported += 1;
                }
                (Err(rusty_occt::occt_brep::Rejected::Unsupported(_)), "unsupported") => {}
                // Representable structure the validator does not certify
                // (S4): pinned below with its issue kinds.
                (Err(rusty_occt::occt_brep::Rejected::Invalid { issues, .. }), "representable") => {
                    let mut kinds: Vec<&'static str> =
                        issues.iter().map(|i| i.kind.name()).collect();
                    kinds.sort_unstable();
                    kinds.dedup();
                    not_certified.push((name.to_string(), row[0].to_string(), kinds));
                }
                (result, verdict) => {
                    panic!("{name} {}: {verdict} but {result:?}", row[0])
                }
            }
        }
    }
    assert_eq!(imported, 58);
    // What the validator cannot certify yet: a sphere face whose loop passes
    // both poles between two seam pairs (the seam merge keeps them), and
    // containment in a cylinder face's loops with spline pcurves.
    let pinned = [
        ("Ball.brep", "108", vec!["seam_edge"]),
        ("Motor-c.brep", "378", vec!["uncertified_containment"]),
    ];
    let pinned: Vec<(String, String, Vec<&str>)> = pinned
        .into_iter()
        .map(|(f, r, k)| (f.to_string(), r.to_string(), k))
        .collect();
    assert_eq!(not_certified, pinned);
}
