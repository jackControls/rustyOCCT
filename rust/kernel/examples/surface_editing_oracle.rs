//! Optional test driver; native OCCT is a separate process.
#[path = "../tests/support/surface_editing_protocol.rs"]
mod protocol;
use std::io::{self, BufRead};
fn main() {
    for row in io::stdin().lock().lines() {
        let input = protocol::parse(&row.unwrap());
        println!(
            "{}",
            protocol::encode(&input.name, &protocol::execute(&input))
        );
    }
}
