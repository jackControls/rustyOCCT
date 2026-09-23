//! Complete exact knot editing observations, including operation success flags.
#[path = "../tests/support/knot_protocol.rs"]
mod protocol;
use std::io::{self, BufRead};
fn main() {
    for row in io::stdin().lock().lines() {
        let input = protocol::parse(&row.unwrap());
        let (flags, curve) = protocol::execute(&input);
        println!("{}", protocol::encode(&input.name, &flags, &curve));
    }
}
