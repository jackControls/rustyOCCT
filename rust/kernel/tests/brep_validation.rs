//! The generic B-rep validator against independently generated expectations.
#[path = "support/brep_protocol.rs"]
mod brep_protocol;
use brep_protocol::parse;
use rusty_occt::Tolerance;
use std::collections::BTreeMap;

#[test]
fn complete_issue_sets_match_the_independent_oracle() {
    let expected: BTreeMap<&str, Vec<&str>> = include_str!("../../fixtures/brep-expected.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let (name, issues) = l.split_once('\t').unwrap();
            (name, issues.split(';').filter(|s| !s.is_empty()).collect())
        })
        .collect();
    let cases = include_str!("../../fixtures/brep-cases.txt");
    let mut checked = 0;
    let mut valid = 0;
    let mut failures = Vec::new();
    for block in cases.split("\nend").filter(|b| !b.trim().is_empty()) {
        let (name, tolerance, parts) = parse(block.trim());
        let tolerance = Tolerance::new(tolerance, 1e-12).unwrap();
        let mut actual: Vec<String> = parts
            .check(tolerance)
            .iter()
            .map(|i| i.to_string())
            .collect();
        actual.sort();
        let mut want: Vec<String> = expected[name.as_str()]
            .iter()
            .map(|s| s.to_string())
            .collect();
        want.sort();
        if actual != want {
            failures.push(format!(
                "{name}\n  expected {want:?}\n  actual   {actual:?}"
            ));
        }
        if want.is_empty() {
            valid += 1;
        }
        checked += 1;
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!((checked, valid), (expected.len(), 21));
}
