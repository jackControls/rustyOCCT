#![no_main]
libfuzzer_sys::fuzz_target!(|data: &[u8]| {
    rusty_occt_fuzz::check_exact_spline_intersections(data);
});
