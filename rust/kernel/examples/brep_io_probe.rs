//! Test-only input for compare_brep_io.py.
//!
//! `prisms OUT` reads the identity case protocol on stdin, writes each prism
//! as `OUT/NAME.brep` and prints `NAME V E W F SH SO volume area cx cy cz`:
//! synthesized counts and the kernel's exact mass properties.
//! `corpus OUT` reads `.brep` paths on stdin, writes every solid it imports
//! as `OUT/FILE-RECORD.brep` and prints `FILE RECORD V E W F SH SO`.
#[path = "../tests/support/identity_protocol.rs"]
#[allow(dead_code)]
mod identity_protocol;
use rusty_occt::occt_brep::{import, read, write};
use rusty_occt::topology::OcctCounts;
use std::io::Read;

fn counts(c: OcctCounts) -> String {
    format!(
        "{} {} {} {} {} {}",
        c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (mode, out) = (args[1].as_str(), std::path::Path::new(&args[2]));
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    match mode {
        "prisms" => {
            for spec in identity_protocol::cases(&input) {
                let solid = identity_protocol::build(&spec);
                let tol = solid.profile().tolerance().linear();
                let text = write(solid.topology(), tol).unwrap();
                std::fs::write(out.join(format!("{}.brep", spec.name)), text).unwrap();
                let m = solid.mass_properties();
                println!(
                    "{} {} {:?} {:?} {:?} {:?} {:?}",
                    spec.name,
                    counts(solid.topology().occt_counts()),
                    m.volume,
                    m.surface_area,
                    m.centroid.x,
                    m.centroid.y,
                    m.centroid.z
                );
            }
        }
        "corpus" => {
            for path in input.lines().filter(|l| !l.trim().is_empty()) {
                let path = std::path::Path::new(path);
                let name = path.file_name().unwrap().to_string_lossy();
                let doc = read(&std::fs::read_to_string(path).unwrap()).unwrap();
                for solid in import(&doc).solids {
                    let Ok(topology) = &solid.result else {
                        continue;
                    };
                    let text = write(topology, solid.tolerance.linear()).unwrap();
                    let file = format!("{name}-{}.brep", solid.record);
                    std::fs::write(out.join(file), text).unwrap();
                    println!("{name} {} {}", solid.record, counts(topology.occt_counts()));
                }
            }
        }
        other => panic!("mode {other}"),
    }
}
