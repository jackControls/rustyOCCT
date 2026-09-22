//! Diagnostic replay includes both kernel and complete independent oracle work.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut paths = std::env::args().skip(1).peekable();
    if paths.peek().is_none() {
        return Err("provide a saved Bezier editing input path".into());
    }
    for path in paths {
        let data = std::fs::read(&path)?;
        let start = std::time::Instant::now();
        rusty_occt_fuzz::check_bezier_editing(&data);
        println!(
            "{path}: {:.6}s, {} bytes, oracle passed",
            start.elapsed().as_secs_f64(),
            data.len()
        );
    }
    Ok(())
}
