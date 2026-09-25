//! Bitwise mass properties, bounds and point classification of every prism
//! fixture, for fixtures/prism-properties-baseline.tsv. The baseline was
//! recorded before the cell-complex migration (T1 of TOPOLOGY_MODEL.md), which
//! must leave every value bitwise unchanged.
use rusty_occt::{Location, Point3, Solid};

fn bits(values: &[f64]) -> String {
    values
        .iter()
        .map(|v| format!("{:016x}", v.to_bits()))
        .collect::<Vec<_>>()
        .join(" ")
}

/// One row: volume, area, centroid, inertia and bounds as binary64 bits, then
/// the classification of a 5x5x5 grid over the bounds grown by a quarter.
pub fn row(name: &str, solid: &Solid) -> String {
    let m = solid.mass_properties();
    let b = solid.bounds();
    let mut values = vec![m.volume, m.surface_area];
    values.extend(m.centroid.to_array());
    values.extend(m.inertia.into_iter().flatten());
    values.extend(b.min.to_array());
    values.extend(b.max.to_array());
    let (lo, hi) = (b.min.to_array(), b.max.to_array());
    let mut classes = String::new();
    for i in 0..5 {
        for j in 0..5 {
            for k in 0..5 {
                let at = |axis: usize, n: usize| {
                    let span = hi[axis] - lo[axis];
                    lo[axis] - 0.25 * span + 1.5 * span * n as f64 / 4.0
                };
                let p = Point3::new(at(0, i), at(1, j), at(2, k));
                classes.push(match solid.classify(p).unwrap() {
                    Location::Outside => 'o',
                    Location::Boundary => 'b',
                    Location::Inside => 'i',
                });
            }
        }
    }
    format!("{name}\t{}\t{classes}", bits(&values))
}
