//! Test-only input for compare_curve_curve.py: reads the case protocol of
//! `curve-curve-cases.txt` on stdin and prints per case `NAME empty`,
//! `NAME coincident`, `NAME limit`, or one row per point sorted by the first
//! curve's parameter: `NAME point sa_lo sa_hi sb_lo sb_hi x_lo x_hi y_lo
//! y_hi z_lo z_hi crossing|tangent`.
#[path = "../tests/support/curve_curve_protocol.rs"]
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
