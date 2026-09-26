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
        let tol = solid.profile().tolerance().linear();
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

const GEOMETRY: [&str; 23] = [
    "BSplineCurve",
    "BSplineCurve2d",
    "BSplineSurface",
    "BezierCurve",
    "BezierCurve2d",
    "BezierSurface",
    "ConicalSurface",
    "Ellipse",
    "Ellipse2d",
    "Hyperbola",
    "Hyperbola2d",
    "OffsetCurve",
    "OffsetCurve2d",
    "OffsetSurface",
    "Parabola",
    "Parabola2d",
    "RectangularTrimmedSurface",
    "SphericalSurface",
    "SurfaceOfLinearExtrusion",
    "SurfaceOfRevolution",
    "ToroidalSurface",
    "TrimmedCurve",
    "TrimmedCurve2d",
];

/// The data/occ corpus against brep_io_reference.py: the same unsupported
/// geometry records, the same solids, a cell topology exactly for the
/// solids whose structure is representable, and synthesized counts equal to
/// OCCT's distinct subshapes of the original seamed solid.
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
                    imported += 1;
                }
                (Err(rusty_occt::occt_brep::Rejected::Unsupported(_)), "unsupported") => {}
                (result, verdict) => {
                    panic!("{name} {}: {verdict} but {result:?}", row[0])
                }
            }
        }
    }
    assert_eq!(imported, 29);
}
