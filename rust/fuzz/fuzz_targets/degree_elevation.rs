#![no_main]

#[cfg(feature = "asan-allocator")]
#[path = "support/allocator.rs"]
mod allocator;

libfuzzer_sys::fuzz_target!(|data: &[u8]| {
    rusty_occt_fuzz::check_degree_elevation(data);
    #[cfg(feature = "asan-allocator")]
    allocator::purge_between_inputs();
});
