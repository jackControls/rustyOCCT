//! Test-only input for compare_curve_surface.py: reads the case protocol of
//! `curve-surface-cases.txt` on stdin and prints per case `NAME empty`,
//! `NAME contained`, `NAME limit`, or one row per point or overlap sorted by
//! parameter: `NAME point s_lo s_hi x_lo x_hi y_lo y_hi z_lo z_hi
//! crossing|tangent`, `NAME overlap s0 s1` (a spline's spans on the surface).
#[path = "../tests/support/curve_surface_protocol.rs"]
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
