//! Test-only protocol for the source-pinned BRepCheck comparison: reads case
//! blocks and prints each case's complete, sorted issue list, or with the
//! `counts` argument each valid case's synthesized OCCT counts (vertices,
//! edges, wires, faces, shells, solids) and `-` for an invalid one, or with
//! `enclosures` each valid case's largest measured vertex, fin and face
//! enclosure (M5) of its parts with every declared bound removed, or with
//! `mass` each valid case's certified mass enclosure (volume, area,
//! centroid, inertia about it row by row, each as `lo hi`) and `-` when it
//! is invalid or not integrated, or with `measure` each valid case's class
//! and, for a sheet, wire or acorn, its certified area or length and centre
//! (S6, each as `lo hi`).
#[path = "../tests/support/brep_protocol.rs"]
mod brep_protocol;
use rusty_occt::topology::Topology;
use rusty_occt::Tolerance;
use std::io::Read;

fn main() {
    let mode = std::env::args().nth(1);
    let counts = mode.as_deref() == Some("counts");
    let enclosures = mode.as_deref() == Some("enclosures");
    let mass = mode.as_deref() == Some("mass");
    let measure = mode.as_deref() == Some("measure");
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
        if mass {
            let row = match Topology::from_parts(parts, tolerance)
                .ok()
                .and_then(|t| t.mass_enclosure())
            {
                Some(m) => {
                    let mut values = vec![m.volume, m.surface_area];
                    values.extend(m.centroid);
                    values.extend(m.inertia.iter().flatten().copied());
                    values
                        .iter()
                        .map(|[lo, hi]| format!("{lo:?} {hi:?}"))
                        .collect::<Vec<_>>()
                        .join(" ")
                }
                None => "-".into(),
            };
            println!("{name}\t{row}");
            continue;
        }
        if measure {
            let row = match Topology::from_parts(parts, tolerance) {
                Ok(t) => {
                    let mut row = t.class().name().to_string();
                    if let Some(m) = t.measure_enclosure() {
                        for [lo, hi] in std::iter::once(m.measure).chain(m.centre) {
                            row += &format!(" {lo:?} {hi:?}");
                        }
                    }
                    row
                }
                Err(_) => "-".into(),
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
