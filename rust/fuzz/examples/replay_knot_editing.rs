//! Diagnostic replay includes kernel work and all independent oracle checks.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut paths = std::env::args().skip(1).peekable();
    if paths.peek().is_none() {
        return Err("provide a saved knot editing input path".into());
    }
    for path in paths {
        let data = std::fs::read(&path)?;
        let start = std::time::Instant::now();
        let times = rusty_occt_fuzz::profile_knot_editing(&data);
        println!(
            "{path}: {:.6}s, {} bytes, oracle passed",
            start.elapsed().as_secs_f64(),
            data.len()
        );
        for (label, time) in times {
            println!("  {label}: {:.6}s", time.as_secs_f64());
        }
    }
    Ok(())
}
