//! Replay a proximity input with its independent global optimality certificate.
//! Times include the kernel and oracle; they are diagnostic, not a CI gate.
use std::{fs, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut paths = std::env::args().skip(1).peekable();
    if paths.peek().is_none() {
        return Err("provide at least one proximity input path".into());
    }
    for path in paths {
        let input = fs::read(&path)?;
        let start = Instant::now();
        rusty_occt_fuzz::check_proximity(&input);
        println!(
            "{path}: {:.6}s, {} bytes, oracle passed",
            start.elapsed().as_secs_f64(),
            input.len()
        );
    }
    Ok(())
}
