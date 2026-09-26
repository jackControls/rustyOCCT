//! Test-only protocol for the source-pinned BRepCheck comparison: reads case
//! blocks and prints each case's complete, sorted issue list, or with the
//! `counts` argument each valid case's synthesized OCCT counts (vertices,
//! edges, wires, faces, shells, solids) and `-` for an invalid one.
#[path = "../tests/support/brep_protocol.rs"]
mod brep_protocol;
use rusty_occt::topology::Topology;
use rusty_occt::Tolerance;
use std::io::Read;

fn main() {
    let counts = std::env::args().nth(1).as_deref() == Some("counts");
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for block in input.split("\nend").filter(|b| !b.trim().is_empty()) {
        let (name, tolerance, parts) = brep_protocol::parse(block.trim());
        let tolerance = Tolerance::new(tolerance, 1e-12).unwrap();
        if counts {
            let row = match Topology::from_parts(parts, tolerance) {
                Ok(t) => {
                    let c = t.occt_counts();
                    format!(
                        "{} {} {} {} {} {}",
                        c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids
                    )
                }
                Err(_) => "-".into(),
            };
            println!("{name}\t{row}");
            continue;
        }
        let mut issues: Vec<String> = parts
            .check(tolerance)
            .iter()
            .map(|i| i.to_string())
            .collect();
        issues.sort();
        println!("{name}\t{}", issues.join(";"));
    }
}
