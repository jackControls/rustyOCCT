//! Test-only protocol for the source-pinned BRepCheck comparison: reads case
//! blocks and prints each case's complete, sorted issue list.
#[path = "../tests/support/brep_protocol.rs"]
mod brep_protocol;
use rusty_occt::Tolerance;
use std::io::Read;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for block in input.split("\nend").filter(|b| !b.trim().is_empty()) {
        let (name, tolerance, parts) = brep_protocol::parse(block.trim());
        let tolerance = Tolerance::new(tolerance, 1e-12).unwrap();
        let mut issues: Vec<String> = parts
            .check(tolerance)
            .iter()
            .map(|i| i.to_string())
            .collect();
        issues.sort();
        println!("{name}\t{}", issues.join(";"));
    }
}
