//! Replay saved .brep interop inputs without mutation.
use std::{fs, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut paths = std::env::args().skip(1).peekable();
    if paths.peek().is_none() {
        return Err("provide at least one brep_io input path".into());
    }
    for path in paths {
        let input = fs::read(&path)?;
        let start = Instant::now();
        rusty_occt_fuzz::check_brep_io(&input);
        println!(
            "{path}: {:.6}s, {} bytes, oracle passed",
            start.elapsed().as_secs_f64(),
            input.len()
        );
    }
    Ok(())
}
