#![no_main]

libfuzzer_sys::fuzz_target!(|data: &[u8]| {
    rusty_occt_fuzz::check_surface_knots(data);
    #[cfg(feature = "asan-allocator")]
    purge_allocator_between_inputs();
});

// The pinned libFuzzer purges unused ASan allocator memory during mutation,
// but not during seed replay. Large exact tensor checks can exhaust the RSS
// limit during replay even after their allocations have all been dropped.
// Apply the same sanitizer operation, at most once per second, after a whole
// input has completed. Live allocations and the per-input RSS limit remain.
// This feature is enabled only by an explicitly address-sanitized fuzz build;
// ordinary replay examples and stable checks have no sanitizer dependency.
#[cfg(feature = "asan-allocator")]
fn purge_allocator_between_inputs() {
    use std::sync::Mutex;
    use std::time::{Duration, Instant};

    unsafe extern "C" {
        fn __sanitizer_purge_allocator();
    }
    static LAST_PURGE: Mutex<Option<Instant>> = Mutex::new(None);
    let mut last = LAST_PURGE.lock().unwrap();
    let now = Instant::now();
    if last.is_none_or(|time| now.duration_since(time) >= Duration::from_secs(1)) {
        // SAFETY: the campaign runner enables this feature together with
        // AddressSanitizer, which exports this no-argument allocator API.
        unsafe { __sanitizer_purge_allocator() };
        *last = Some(now);
    }
}
