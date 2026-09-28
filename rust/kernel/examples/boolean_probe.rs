//! Test-only input for compare_boolean.py: reads the case protocol of
//! `boolean-cases.txt` on stdin and prints per case `NAME limit`, `NAME
//! unsupported`, `NAME refused`, `NAME empty`, or per result solid `NAME
//! solid vol_lo vol_hi area_lo area_hi cx_lo cx_hi cy_lo cy_hi cz_lo cz_hi
//! faces edges vertices`.
#[path = "../tests/support/boolean_protocol.rs"]
mod protocol;
use std::io::Read;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    for case in protocol::cases(&input) {
        let rows = protocol::rows(&case).unwrap_or_else(|e| panic!("{}: {e}", case.name));
        for row in rows {
            println!("{} {row}", case.name);
        }
    }
}
