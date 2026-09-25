//! Mass properties, bounds and classification of every prism fixture are
//! bitwise equal to the baseline recorded before the cell-complex migration.
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

#[test]
fn every_prism_property_is_bitwise_unchanged() {
    let got = rows();
    if std::env::var_os("RUSTY_WRITE_PROPERTY_BASELINE").is_some() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../fixtures/prism-properties-baseline.tsv"
        );
        std::fs::write(path, got.join("\n") + "\n").unwrap();
    }
    let want: Vec<&str> = include_str!("../../fixtures/prism-properties-baseline.tsv")
        .lines()
        .collect();
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(&want) {
        assert_eq!(g, w);
    }
}
