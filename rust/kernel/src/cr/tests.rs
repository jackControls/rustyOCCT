//! The tables against the rational tier, the `Fast` tier's every decision
//! against the rational tier's on deterministic samples, and the special
//! values. The bits themselves are pinned against an independent mpmath
//! oracle in `tests/correctly_rounded.rs`.
use super::*;

fn dd_r(x: Dd) -> R {
    exact(x.hi) + exact(x.lo)
}

/// `|a - b| <= 2^k |b|`.
fn within(a: &R, b: &R, k: i64) -> bool {
    let d = a - b;
    let d = if d < int(0) { -d } else { d };
    let m = if b < &int(0) { -b.clone() } else { b.clone() };
    d <= m * pow2_r(k)
}

fn inside(x: &R, i: &Interval, slack: &R) -> bool {
    &(i.lo() - slack) <= x && x <= &(i.hi() + slack)
}

fn fact(n: i64) -> R {
    (1..=n).fold(int(1), |p, m| p * int(m))
}

#[test]
fn tables_agree_with_the_rational_tier() {
    for (n, c) in INV_FACT.iter().enumerate() {
        assert!(
            within(&dd_r(*c), &(int(1) / fact(n as i64)), -105),
            "1/{n}!"
        );
    }
    for (k, c) in INV_ODD.iter().enumerate() {
        assert!(
            within(
                &dd_r(*c),
                &R::new(1.into(), (2 * k as i64 + 1).into()),
                -105
            ),
            "1/(2{k}+1)"
        );
    }
    type Exact = Box<dyn Fn(i64) -> R>;
    let tails: [(&[f64], Exact); 7] = [
        (
            &SIN_TAIL,
            Box::new(|i| int(if (i + 7) % 2 == 0 { 1 } else { -1 }) / fact(2 * (i + 7) + 1)),
        ),
        (
            &COS_TAIL,
            Box::new(|i| int(if (i + 7) % 2 == 0 { 1 } else { -1 }) / fact(2 * (i + 7))),
        ),
        (
            &ATAN_TAIL,
            Box::new(|i| {
                R::new(
                    if (i + 4) % 2 == 0 {
                        1.into()
                    } else {
                        (-1).into()
                    },
                    (2 * (i + 4) + 1).into(),
                )
            }),
        ),
        (
            &LN_TAIL,
            Box::new(|i| R::new(1.into(), (2 * (i + 7) + 1).into())),
        ),
        (&EXP_TAIL, Box::new(|i| int(1) / fact(i + 11))),
        (&SINH_TAIL, Box::new(|i| int(1) / fact(2 * (i + 6) + 1))),
        (&COSH_TAIL, Box::new(|i| int(1) / fact(2 * (i + 6)))),
    ];
    for (t, f) in &tails {
        for (i, c) in t.iter().enumerate() {
            assert!(within(&exact(*c), &f(i as i64), -48), "tail {i}: {c}");
        }
    }
    let pi = certified::pi();
    let half = R::new(1.into(), 2.into());
    let half_pi = pi.scale(&half);
    let parts = HALF_PI_PARTS.iter().fold(int(0), |s, p| s + exact(*p));
    assert!(inside(&parts, &half_pi, &pow2_r(-158)), "pi/2's parts");
    for p in &HALF_PI_PARTS[..3] {
        // At most 33 significant bits: k p exact for |k| < 2^20.
        assert_eq!(p.to_bits() & ((1 << 20) - 1), 0, "{p}");
    }
    assert!(inside(&dd_r(PI), pi, &pow2_r(-104)), "pi");
    assert!(inside(&dd_r(HALF_PI), &half_pi, &pow2_r(-105)), "pi/2");
    let fixed = R::new(half_pi_fixed().clone(), BigInt::from(1) << BIG);
    assert!(inside(&fixed, &half_pi, &pow2_r(-180)), "pi/2's 1344 bits");
    for (j, c) in ATAN_TABLE.iter().enumerate().skip(1) {
        let a = certified::atan(&Interval::exact(R::new((j as i64).into(), 16.into())));
        assert!(inside(&dd_r(*c), &a, &pow2_r(-106)), "atan({j}/16)");
    }
    let ln2 = certified::ln2();
    let parts = LN2_PARTS.iter().fold(int(0), |s, p| s + exact(*p));
    assert!(inside(&parts, ln2, &pow2_r(-138)), "ln 2's parts");
    for p in &LN2_PARTS[..2] {
        assert_eq!(p.to_bits() & ((1 << 11) - 1), 0, "{p}");
    }
    assert!(inside(&dd_r(LN2), ln2, &pow2_r(-106)), "ln 2");
}

/// Deterministic samples: a linear congruential generator's bits mapped
/// into each function's domain, with every binade weighted alike.
struct Samples(u64);

impl Samples {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }
    /// A number with exponent in `lo..hi`, random significand and sign.
    fn binade(&mut self, lo: i32, hi: i32, signed: bool) -> f64 {
        let r = self.next();
        let e = lo + (r % (hi - lo) as u64) as i32;
        let m = 1.0 + (self.next() >> 12) as f64 * pow2(-52);
        let v = m * pow2(e);
        if signed && r & (1 << 40) != 0 {
            -v
        } else {
            v
        }
    }
}

fn count() -> usize {
    std::env::var("RUSTY_CR_SAMPLES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(60)
}

/// The function's binary64 is the rational tier's wherever that decides.
fn check(name: &str, got: f64, want: Option<f64>, arg: &str) {
    if let Some(w) = want {
        assert_eq!(
            got.to_bits(),
            w.to_bits(),
            "{name}({arg}): {got:e} against {w:e}"
        );
    }
}

#[test]
fn every_fast_decision_is_the_rational_tiers() {
    let n = count();
    let mut s = Samples(0x5eed_0fc0_ffee);
    for _ in 0..n {
        let x = s.binade(-26, 22, true);
        let red = reduce_precise(x);
        let (si, ci) = sin_cos_rational(&red);
        let (gs, gc) = sin_cos(x);
        check("sin", gs, round_enclosure(&si), &format!("{x:e}"));
        check("cos", gc, round_enclosure(&ci), &format!("{x:e}"));
        check(
            "tan",
            tan(x),
            si.div(&ci).as_ref().and_then(round_enclosure),
            &format!("{x:e}"),
        );
    }
    for _ in 0..n {
        let x = s.binade(-27, 30, true);
        let i = certified::atan(&Interval::from_f64(x));
        check("atan", atan(x), round_enclosure(&i), &format!("{x:e}"));
        let y = s.binade(-30, 30, true);
        let i = certified::atan2(&Interval::from_f64(y), &Interval::from_f64(x));
        check(
            "atan2",
            atan2(y, x),
            i.as_ref().and_then(round_enclosure),
            &format!("{y:e}, {x:e}"),
        );
    }
    for _ in 0..n {
        let x = s.binade(-26, 0, true);
        if x.abs() >= 1.0 {
            continue;
        }
        let c = cosine_interval(x);
        let i = certified::atan2(&Interval::from_f64(x), &c);
        check(
            "asin",
            asin(x),
            i.as_ref().and_then(round_enclosure),
            &format!("{x:e}"),
        );
        let i = certified::atan2(&c, &Interval::from_f64(x));
        check(
            "acos",
            acos(x),
            i.as_ref().and_then(round_enclosure),
            &format!("{x:e}"),
        );
    }
    for _ in 0..n {
        let x = s.binade(-40, 40, false);
        let y = s.binade(-40, 40, true);
        let h = hypot(x, y);
        // The exact squares' answer.
        assert_eq!(
            h.to_bits(),
            hypot_exact(x.abs().max(y.abs()), x.abs().min(y.abs())).to_bits(),
            "hypot({x:e}, {y:e})"
        );
    }
    for _ in 0..n {
        let x = s.binade(-50, 9, true);
        if x.abs() > 700.0 {
            continue;
        }
        let (k, i) = exp_rational(x);
        check(
            "exp",
            exp(x),
            round_enclosure(&i).map(|v| ldexp(v, k)),
            &format!("{x:e}"),
        );
        if x.abs() < 40.0 && x.abs() >= pow2(-27) {
            let a = x.abs();
            let (k, i) = exp_rational(a);
            let e = i.scale(&pow2_r(i64::from(k)));
            let half = R::new(1.into(), 2.into());
            let inv = Interval::exact(int(1)).div(&e).unwrap();
            let ch = e.add(&inv).scale(&half);
            let sh = e.sub(&inv).scale(&half);
            check("cosh", cosh(x), round_enclosure(&ch), &format!("{x:e}"));
            let sign = if x < 0.0 { -1.0 } else { 1.0 };
            check(
                "sinh",
                sinh(x),
                round_enclosure(&sh).map(|v| v * sign),
                &format!("{a:e}"),
            );
        }
    }
    for _ in 0..n {
        let x = s.binade(-1000, 1000, false);
        let i = certified::ln(&Interval::from_f64(x));
        check(
            "ln",
            ln(x),
            i.as_ref().and_then(round_enclosure),
            &format!("{x:e}"),
        );
        let i = i.and_then(|l| l.div(certified::ln2()));
        check(
            "log2",
            log2(x),
            i.as_ref().and_then(round_enclosure),
            &format!("{x:e}"),
        );
    }
}

/// `|v - i| / |i|` against an enclosure's midpoint, as a binary64.
fn rel_err(v: Dd, i: &Interval) -> f64 {
    let m = (i.lo() + i.hi()) / int(2);
    let d = dd_r(v) - &m;
    if d == int(0) {
        return 0.0;
    }
    approx(&(d / m)).map_or(0.0, f64::abs)
}

/// The analyses' bounds hold on samples with room to spare: the quick
/// sine's and cosine's errors stay under `2^-65` (their budget `2^-62`),
/// the tabled evaluations' under `2^-68` (budget `2^-64`), the series'
/// under `2^-92` (budget `2^-88`). The reduced arguments carry low parts,
/// as the reductions' do.
#[test]
fn evaluation_errors_stay_inside_their_budgets() {
    let n = 4 * count();
    let mut s = Samples(0xb0d9e7);
    let (mut worst_q, mut worst_t, mut worst_s) = (0.0f64, 0.0f64, 0.0f64);
    for _ in 0..n {
        // A reduced argument and its sine and cosine.
        let hi = s.binade(-12, 0, true) * 0.78;
        let (h, l) = fast_two_sum(hi, s.binade(-60, -54, true) * hi.abs());
        let r = Dd { hi: h, lo: l };
        let (ci, si) = certified::cos_sin(&dd_r(r));
        let (sq, cq) = sin_cos_quick(r);
        let (st, ct) = sin_cos_tabled(r);
        let (ss, cs) = sin_cos_reduced(r);
        worst_q = worst_q.max(rel_err(sq, &si)).max(rel_err(cq, &ci));
        worst_t = worst_t.max(rel_err(st, &si)).max(rel_err(ct, &ci));
        worst_s = worst_s.max(rel_err(ss, &si)).max(rel_err(cs, &ci));
        // A ratio's arctangent.
        let (y, x) = (s.binade(-12, 0, false), s.binade(0, 1, false));
        let (y, x) = (y.min(x), x.max(y));
        let i = certified::atan2(&Interval::from_f64(y), &Interval::from_f64(x)).unwrap();
        worst_t = worst_t.max(rel_err(atan_tabled(Dd::f(y), Dd::f(x)), &i));
        worst_s = worst_s.max(rel_err(atan_unit(Dd::f(y).div(Dd::f(x))), &i));
    }
    assert!(worst_q < pow2(-65), "quick {worst_q:e}");
    assert!(worst_t < pow2(-68), "tabled {worst_t:e}");
    assert!(worst_s < pow2(-92), "series {worst_s:e}");
}

/// The 128-bit integer `hypot` is the exact one, from any nearby guess,
/// and Pythagorean triples (the ties' and exact results' source) come out
/// exact at any scale.
#[test]
fn hypots_in_integers_are_the_exact_ones() {
    let mut s = Samples(0x0004_b907);
    for _ in 0..20 * count() {
        let a = s.binade(0, 1, false);
        let b = s.binade(-9, 1, false).min(a);
        let guess = next_up((a * a + b * b).sqrt());
        assert_eq!(
            hypot_scaled(a, b, guess).to_bits(),
            hypot_exact(a, b).to_bits(),
            "{a:e} {b:e}"
        );
    }
    for u in 2..40u32 {
        for v in 1..u {
            let (p, q, r) = (u * u - v * v, 2 * u * v, u * u + v * v);
            for k in [-1060, -300, -7, 0, 33, 900] {
                let f = |n: u32| f64::from(n) * 2f64.powi(k);
                assert_eq!(hypot(f(p), f(q)), f(r), "{p} {q} {r} 2^{k}");
            }
        }
    }
}

#[test]
fn inverse_hyperbolic_decisions_are_the_rational_tiers() {
    let n = count();
    let mut s = Samples(0x00a5_17a7);
    for _ in 0..n {
        let x = s.binade(-27, 40, true);
        let i = certified::asinh(&Interval::from_f64(x));
        check("asinh", asinh(x), round_enclosure(&i), &format!("{x:e}"));
        let r = s.binade(-27, 0, true);
        if r.abs() < 1.0 {
            let q = exact(r);
            let i = certified::ln(&Interval::exact((int(1) + &q) / (int(1) - &q)))
                .map(|l| l.scale(&R::new(1.into(), 2.into())));
            check(
                "atanh",
                atanh(r),
                i.as_ref().and_then(round_enclosure),
                &format!("{r:e}"),
            );
        }
    }
}

/// Arguments within an ulp or two of multiples of `pi/2` (the sine's or
/// cosine's tiny results), from the 1344-bit reduction, and huge ones.
#[test]
fn sines_near_multiples_of_half_pi_and_of_huge_arguments() {
    let mut xs = Vec::new();
    for k in [1i64, 2, 3, 4, 7, 100, 355, 1_000_003, 1 << 30] {
        let v = std::f64::consts::FRAC_PI_2 * k as f64;
        let mut x = next_down(next_down(v));
        for _ in 0..5 {
            xs.push(x);
            x = next_up(x);
        }
    }
    xs.extend([
        1e22,
        -1e22,
        1e300,
        2.0f64.powi(1023),
        5.319372648326541e255,
        6381956970095103.0 * 2.0f64.powi(797),
    ]);
    for x in xs {
        let red = reduce_precise(x);
        let (si, ci) = sin_cos_rational(&red);
        let (s, c) = sin_cos(x);
        check("sin", s, round_enclosure(&si), &format!("{x:e}"));
        check("cos", c, round_enclosure(&ci), &format!("{x:e}"));
        assert!(s.abs() <= 1.0 && c.abs() <= 1.0);
    }
}

#[test]
fn special_values() {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI as P};
    let nan = f64::NAN;
    let inf = f64::INFINITY;
    let same = |a: f64, b: f64| a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan());
    let cases: &[(&str, f64, f64)] = &[
        ("sin -0", sin(-0.0), -0.0),
        ("cos -0", cos(-0.0), 1.0),
        ("tan -0", tan(-0.0), -0.0),
        ("sin inf", sin(inf), nan),
        ("atan inf", atan(inf), FRAC_PI_2),
        ("atan -inf", atan(-inf), -FRAC_PI_2),
        ("atan2 +0 +0", atan2(0.0, 0.0), 0.0),
        ("atan2 -0 +0", atan2(-0.0, 0.0), -0.0),
        ("atan2 +0 -0", atan2(0.0, -0.0), P),
        ("atan2 -0 -0", atan2(-0.0, -0.0), -P),
        ("atan2 +0 -1", atan2(0.0, -1.0), P),
        ("atan2 1 0", atan2(1.0, 0.0), FRAC_PI_2),
        ("atan2 -1 -0", atan2(-1.0, -0.0), -FRAC_PI_2),
        ("atan2 1 -inf", atan2(1.0, -inf), P),
        ("atan2 -1 inf", atan2(-1.0, inf), -0.0),
        ("atan2 inf inf", atan2(inf, inf), FRAC_PI_4),
        ("atan2 -inf -inf", atan2(-inf, -inf), -3.0 * FRAC_PI_4),
        ("atan2 inf 1", atan2(inf, 1.0), FRAC_PI_2),
        ("atan2 1 1", atan2(1.0, 1.0), FRAC_PI_4),
        ("atan2 -1 -1", atan2(-1.0, -1.0), -3.0 * FRAC_PI_4),
        ("asin 1", asin(1.0), FRAC_PI_2),
        ("asin -1", asin(-1.0), -FRAC_PI_2),
        ("asin 2", asin(2.0), nan),
        ("acos 1", acos(1.0), 0.0),
        ("acos -1", acos(-1.0), P),
        ("acos 0", acos(0.0), FRAC_PI_2),
        ("hypot inf nan", hypot(inf, nan), inf),
        ("hypot nan 1", hypot(nan, 1.0), nan),
        ("hypot -0 0", hypot(-0.0, 0.0), 0.0),
        ("hypot 3 4", hypot(3.0, 4.0), 5.0),
        ("hypot -5 12", hypot(-5.0, 12.0), 13.0),
        ("exp 0", exp(0.0), 1.0),
        ("exp -inf", exp(-inf), 0.0),
        ("exp 1000", exp(1000.0), inf),
        ("ln 1", ln(1.0), 0.0),
        ("ln 0", ln(0.0), -inf),
        ("ln -1", ln(-1.0), nan),
        ("log2 1024", log2(1024.0), 10.0),
        ("log2 min subnormal", log2(f64::from_bits(1)), -1074.0),
        ("log2 0.5", log2(0.5), -1.0),
        ("sinh -0", sinh(-0.0), -0.0),
        ("cosh 0", cosh(0.0), 1.0),
        ("cosh -inf", cosh(-inf), inf),
        ("sinh -inf", sinh(-inf), -inf),
    ];
    for (name, got, want) in cases {
        assert!(same(*got, *want), "{name}: {got:e} against {want:e}");
    }
    // The hypot of a Pythagorean triple's subnormal and huge multiples.
    let t = f64::from_bits(3);
    assert_eq!(hypot(t, f64::from_bits(4)), f64::from_bits(5));
    assert_eq!(hypot(3.0 * pow2(1000), 4.0 * pow2(1000)), 5.0 * pow2(1000));
    assert_eq!(hypot(f64::MAX, f64::MAX), inf);
}
