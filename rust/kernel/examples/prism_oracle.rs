//! Test-only differential oracle driver; not an application interchange format.
#[path = "../tests/support/mod.rs"]
mod support;
use std::io::{self, Read};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    for observation in support::evaluate(&input)? {
        println!("{observation}");
    }
    Ok(())
}
