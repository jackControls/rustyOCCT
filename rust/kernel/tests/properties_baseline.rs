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
        // The baseline predates cones (S3); it covers the prisms.
        if spec.cone.is_some() {
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

/// The committed baseline was recorded on this host before the migration
/// (at 26fc457f). Frames and rotations use the platform's trigonometry, so
/// other hosts differ in the last bits; CI regenerates the baseline on each
/// such host from that revision and passes it in `RUSTY_PROPERTY_BASELINE`.
const RECORDED_ON: (&str, &str) = ("macos", "aarch64");

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
    let host = (std::env::consts::OS, std::env::consts::ARCH);
    let text = match std::env::var_os("RUSTY_PROPERTY_BASELINE") {
        Some(path) => std::fs::read_to_string(path).unwrap(),
        None if host == RECORDED_ON => {
            include_str!("../../fixtures/prism-properties-baseline.tsv").to_string()
        }
        None => {
            eprintln!("no property baseline recorded on {host:?}; set RUSTY_PROPERTY_BASELINE");
            return;
        }
    };
    let want: Vec<&str> = text.lines().collect();
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(&want) {
        assert_eq!(g, w);
    }
}
