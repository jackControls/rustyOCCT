//! Full surface knot observations, including exact controls and removal flags.
#[path = "../tests/support/surface_knot_protocol.rs"]
mod protocol;
use std::io::{self, BufRead};
fn main() {
    for row in io::stdin().lock().lines() {
        let input = protocol::parse(&row.unwrap());
        let (flags, surface) = protocol::execute(&input);
        println!("{}", protocol::encode(&input.name, &flags, &surface, false));
    }
}
