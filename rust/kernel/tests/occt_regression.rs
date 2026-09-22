mod support;

/// Expected values come from the separately compiled OCCT oracle, not Rust.
/// This keeps the compatibility corpus runnable without any C++ dependencies.
#[test]
fn matches_recorded_occt_geometry_and_classification() {
    let fixtures = include_str!("../../fixtures/prisms.txt");
    let baseline = include_str!("../../fixtures/occt-baseline.tsv");
    let actual = support::evaluate(fixtures).unwrap();
    let rows = baseline.lines().collect::<Vec<_>>();
    assert_eq!(actual.len(), rows.len());
    let mut classifications = 0;
    for (actual, expected) in actual.iter().zip(rows) {
        let mut tokens = expected.split_whitespace();
        assert_eq!(tokens.next().unwrap(), actual.name);
        let values = tokens
            .map(|value| value.parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(values.len(), actual.values.len());
        assert!(values.len() >= 23);
        for (index, (a, e)) in actual.values.iter().zip(values).enumerate() {
            assert!(a.is_finite() && e.is_finite());
            if index >= 20 {
                assert_eq!(
                    *a, e,
                    "{}: topology/classification column {index}",
                    actual.name
                );
            } else {
                let absolute = if index >= 11 { 1e-6 } else { 1e-7 };
                let budget = absolute + 2e-9 * e.abs();
                assert!(
                    (a - e).abs() <= budget,
                    "{}: property column {index}: {a} != {e} (+/- {budget})",
                    actual.name
                );
            }
        }
        classifications += actual.values.len() - 23;
    }
    assert_eq!(actual.len(), 66);
    assert_eq!(classifications, 2292);
}
