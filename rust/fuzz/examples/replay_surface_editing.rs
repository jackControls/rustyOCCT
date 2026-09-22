fn main() {
    for path in std::env::args().skip(1) {
        let data = std::fs::read(&path).unwrap();
        let start = std::time::Instant::now();
        let stages = rusty_occt_fuzz::profile_surface_editing(&data);
        println!("{path}: {:?} {stages:?}", start.elapsed());
    }
}
