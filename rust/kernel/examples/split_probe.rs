//! Test-only input for compare_split.py: reads the case protocol of
//! `split-cases.txt` (and the lines of `split-primitive-cases.txt`) on stdin and prints per case `NAME limit`, `NAME
//! unsupported` (a plane or solid of a later sub-step), or per piece `NAME
//! piece below|above vol_lo vol_hi area_lo area_hi cx_lo cx_hi cy_lo cy_hi
//! cz_lo cz_hi faces edges vertices`.
#[path = "../tests/support/split_protocol.rs"]
mod protocol;
use std::io::Read;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    // Prism blocks, then (S8c) one cone or sphere per line.
    let prisms: String = input
        .lines()
        .filter(|l| !(l.starts_with("cone ") || l.starts_with("sphere ")))
        .map(|l| format!("{l}\n"))
        .collect();
    for case in protocol::cases(&prisms) {
        let rows = protocol::rows(&case).unwrap_or_else(|e| panic!("{}: {e}", case.spec.name));
        for row in rows {
            println!("{} {row}", case.spec.name);
        }
    }
    for case in protocol::primitive_cases(&input) {
        let rows = protocol::primitive_rows(&case).unwrap_or_else(|e| panic!("{}: {e}", case.name));
        for row in rows {
            println!("{} {row}", case.name);
        }
    }
}
