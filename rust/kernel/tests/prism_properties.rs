//! Mass properties, bounds and classification of every prism fixture are
//! bitwise those of `fixtures/prism-properties.tsv`, on every host (R14 of
//! the review of 18 October 2026). The kernel's constructions are
//! correctly rounded (R13), so one committed fixture serves every host; it
//! replaces the pre-migration baseline that each host regenerated at
//! `26fc457f` (T1's requirement, met at `e4adb869`).
#[path = "support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
#[path = "support/properties.rs"]
mod properties;
#[allow(dead_code)]
mod support;
use identity_protocol::{build, cases};
use rusty_occt::identity::OperationId;

pub fn rows() -> Vec<String> {
    let mut out = Vec::new();
    for spec in cases(include_str!("../../fixtures/identity-cases.txt")) {
        // The prisms (the fixture predates cones, S3).
        if spec.cone.is_some() || spec.sphere.is_some() || spec.torus.is_some() {
            continue;
        }
        let mut solid = build(&spec);
        out.push(properties::row(&spec.name, &solid));
        for (k, transform) in spec.transforms.iter().enumerate() {
            solid = solid.transform_with(OperationId(1), *transform).unwrap().0;
            out.push(properties::row(&format!("{}@{k}", spec.name), &solid));
        }
    }
    for (name, solid) in support::solids(include_str!("../../fixtures/prisms.txt")).unwrap() {
        out.push(properties::row(&name, &solid));
    }
    out
}

/// `RUSTY_WRITE_PRISM_PROPERTIES=1` writes the rows instead.
#[test]
fn every_prism_property_is_bitwise_the_fixtures() {
    let got = rows();
    if std::env::var_os("RUSTY_WRITE_PRISM_PROPERTIES").is_some() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../fixtures/prism-properties.tsv"
        );
        std::fs::write(path, got.join("\n") + "\n").unwrap();
    }
    let want: Vec<&str> = include_str!("../../fixtures/prism-properties.tsv")
        .lines()
        .collect();
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(&want) {
        assert_eq!(g, w);
    }
}
