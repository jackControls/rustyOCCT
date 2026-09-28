//! Test-only input for compare_step.py.
//!
//! Reads STEP file paths on stdin and imports each (`step::read`,
//! `step::import`). Per body it prints `NAME ENTITY KIND ok V E W F SH SO`
//! (KIND `solid` or `sheet`, the synthesized OCCT counts) and the certified
//! measures as `lo hi` pairs: a solid's volume, area and centre, a sheet's
//! area and centre; or `NAME ENTITY KIND unsupported NAMES` or `NAME ENTITY
//! KIND invalid ISSUES`. A file that does not read prints `NAME error
//! MESSAGE`.
use rusty_occt::occt_brep::Rejected;
use rusty_occt::step::{self, Item};
use std::io::Read;

fn pairs(values: impl IntoIterator<Item = [f64; 2]>) -> String {
    values
        .into_iter()
        .map(|[lo, hi]| format!("{lo:?} {hi:?}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for path in input.lines().filter(|l| !l.trim().is_empty()) {
        let path = std::path::Path::new(path);
        let name = path.file_stem().unwrap().to_string_lossy();
        let imported = step::read(&std::fs::read(path).unwrap()).and_then(|x| step::import(&x));
        let imported = match imported {
            Ok(imported) => imported,
            Err(error) => {
                println!("{name} error {error}");
                continue;
            }
        };
        for body in imported.bodies {
            let kind = match body.item {
                Item::Solid => "solid",
                Item::Shell => "sheet",
            };
            let head = format!("{name} {} {kind}", body.entity);
            match &body.result {
                Ok(t) => {
                    let c = t.occt_counts();
                    let counts = format!(
                        "{} {} {} {} {} {}",
                        c.vertices, c.edges, c.wires, c.faces, c.shells, c.solids
                    );
                    let measures = match body.item {
                        Item::Solid => t.mass_enclosure().map(|m| {
                            pairs([m.volume, m.surface_area].into_iter().chain(m.centroid))
                        }),
                        Item::Shell => t
                            .measure_enclosure()
                            .map(|m| pairs(std::iter::once(m.measure).chain(m.centre))),
                    };
                    println!("{head} ok {counts} {}", measures.unwrap_or("-".into()));
                }
                Err(Rejected::Unsupported(names)) => {
                    println!("{head} unsupported {}", names.join(" "))
                }
                Err(Rejected::Invalid { issues, .. }) => println!(
                    "{head} invalid {}",
                    issues
                        .iter()
                        .map(|i| i.to_string())
                        .collect::<Vec<_>>()
                        .join(" ")
                ),
            }
        }
    }
}
