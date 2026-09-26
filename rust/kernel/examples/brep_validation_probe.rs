//! Test-only protocol for the source-pinned BRepCheck comparison: reads case
//! blocks and prints each case's complete, sorted issue list, or with the
//! `counts` argument each valid case's synthesized OCCT counts (vertices,
//! edges, wires, faces, shells, solids) and `-` for an invalid one, or with
//! `enclosures` each valid case's largest measured vertex, fin and face
//! enclosure (M5) of its parts with every declared bound removed.
#[path = "../tests/support/brep_protocol.rs"]
mod brep_protocol;
use rusty_occt::topology::Topology;
use rusty_occt::Tolerance;
use std::io::Read;

fn main() {
    let mode = std::env::args().nth(1);
    let counts = mode.as_deref() == Some("counts");
    let enclosures = mode.as_deref() == Some("enclosures");
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for block in input.split("\nend").filter(|b| !b.trim().is_empty()) {
        let (name, tolerance, parts) = brep_protocol::parse(block.trim());
        let tolerance = Tolerance::new(tolerance, 1e-12).unwrap();
        if enclosures {
            let row = if parts.check(tolerance).is_empty() {
                let mut bare = parts;
                bare.vertices.iter_mut().for_each(|v| v.enclosure = None);
                bare.fins.iter_mut().for_each(|f| f.enclosure = None);
                bare.faces.iter_mut().for_each(|f| f.enclosure = None);
                let m = bare.with_measured_enclosures();
                let largest = |bounds: Vec<Option<f64>>| {
                    bounds
                        .into_iter()
                        .map(|b| b.expect("valid parts are measurable"))
                        .fold(0.0, f64::max)
                };
                format!(
                    "{:?} {:?} {:?}",
                    largest(
                        m.vertices
                            .iter()
                            .map(|v| v.enclosure.map(|e| e.bound))
                            .collect()
                    ),
                    largest(
                        m.fins
                            .iter()
                            .map(|f| f.enclosure.map(|e| e.bound))
                            .collect()
                    ),
                    largest(
                        m.faces
                            .iter()
                            .map(|f| f.enclosure.map(|e| e.bound))
                            .collect()
                    )
                )
            } else {
                "-".into()
            };
            println!("{name}\t{row}");
            continue;
        }
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
