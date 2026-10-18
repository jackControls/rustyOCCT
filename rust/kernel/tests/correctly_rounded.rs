//! R13 (the review of 18 October 2026): the kernel's elementary functions
//! are correctly rounded, so every host builds the same bits. Every row of
//! `fixtures/correctly-rounded.tsv` (mpmath at 2400 bits, rounded once;
//! `rust/tools/generate_correctly_rounded_fixtures.py`) is reproduced bit
//! for bit, the two quantities X1 found among them.
use rusty_occt::cr;

fn f64_of(word: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(word, 16).unwrap())
}

#[test]
fn every_function_matches_the_independent_oracle_bitwise() {
    let mut failures = Vec::new();
    let mut rows = 0;
    for line in include_str!("../../fixtures/correctly-rounded.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let words: Vec<&str> = line.split('\t').collect();
        let (name, label) = (words[0], words[1]);
        let args: Vec<f64> = words[2..words.len() - 1]
            .iter()
            .map(|w| f64_of(w))
            .collect();
        let want = f64_of(words[words.len() - 1]);
        let x = args[0];
        let got = match name {
            "sin" => cr::sin(x),
            "cos" => cr::cos(x),
            "tan" => cr::tan(x),
            "atan" => cr::atan(x),
            "atan2" => cr::atan2(x, args[1]),
            "asin" => cr::asin(x),
            "acos" => cr::acos(x),
            "hypot" => cr::hypot(x, args[1]),
            "exp" => cr::exp(x),
            "ln" => cr::ln(x),
            "log2" => cr::log2(x),
            "sinh" => cr::sinh(x),
            "cosh" => cr::cosh(x),
            "asinh" => cr::asinh(x),
            "atanh" => cr::atanh(x),
            other => panic!("function {other}"),
        };
        rows += 1;
        if got.to_bits() != want.to_bits() {
            failures.push(format!("{name} {label} {args:?}: {got:e} for {want:e}"));
        }
        // sin_cos agrees with sin and cos.
        if name == "sin" || name == "cos" {
            let (s, c) = cr::sin_cos(x);
            let one = if name == "sin" { s } else { c };
            assert_eq!(one.to_bits(), got.to_bits(), "sin_cos {x:e}");
        }
    }
    assert!(rows > 1900, "{rows} rows");
    assert!(failures.is_empty(), "{failures:#?}");
}
