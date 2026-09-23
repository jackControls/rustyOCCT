fn main() -> Result<(), Box<dyn std::error::Error>> {
    let paths: Vec<_> = std::env::args().skip(1).collect();
    if paths.is_empty() {
        return Err("provide saved exact intersection input paths".into());
    }
    for path in paths {
        let data = std::fs::read(&path)?;
        let start = std::time::Instant::now();
        let times = rusty_occt_fuzz::profile_exact_spline_intersections(&data);
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
