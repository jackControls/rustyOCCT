//! General mass properties with certified enclosures (REVIEW_NOTES.md U2).
//!
//! Volume, first and second moments come from the divergence theorem as
//! surface integrals over the boundary of every solid region, and surface
//! area and face centres from `|S_u x S_v|`. On a face, every integrand is a
//! polynomial in `(u, v)` on a plane, or in `v` times a trigonometric
//! polynomial in `u` on a cylinder or cone; Green's theorem turns the face
//! integral into `-loop integral of F du`, `F` the integrand's antiderivative
//! in `v`. Along a line pcurve that is a polynomial in the line's parameter
//! times `cos` and `sin` of integer multiples of `u`, integrated exactly by
//! recursion on the moments; along a plane arc a trigonometric polynomial,
//! integrated through its exact Fourier expansion. Everything is evaluated
//! in the validator's certified tiers, so the results are enclosures, not
//! estimates. Integration runs in coordinates relative to a reference point
//! of the body, so far-away bodies do not cancel.
//!
//! On a cone face with a pole the antiderivative starts at the apex, where it
//! vanishes, so the pole needs no term of its own; on any other face it
//! starts at `v = 0`, since a nearly cylindrical cone's far apex would make
//! the terms of a band's two loops cancel.
use super::*;

/// Integrands, per face: the volume term, the three first moments, the six
/// second moments, `|N|` and the three `p |N|` of the face centre.
const TERMS: usize = 14;

/// `sum of c u^a v^b` on a plane, keyed `(a, b)`.
type Planar<T> = BTreeMap<(u8, u8), T>;
/// `sum of c v^k cos^a u sin^b u` on a surface of revolution, keyed `(k, a, b)`.
type Rev<T> = BTreeMap<(u8, u8, u8), T>;
/// `sum of c cos^a u sin^b u cos^c v sin^d v` on a sphere, keyed
/// `(a, b, c, d)`.
type Sph<T> = BTreeMap<(u8, u8, u8, u8), T>;

fn add_to<K: Ord + Copy, T: Real>(map: &mut BTreeMap<K, T>, key: K, value: T) {
    let entry = map.entry(key).or_insert_with(|| c(0.0));
    *entry = entry.add(&value);
}

fn planar_mul<T: Real>(a: &Planar<T>, b: &Planar<T>) -> Planar<T> {
    let mut out = BTreeMap::new();
    for ((i, j), x) in a {
        for ((k, l), y) in b {
            add_to(&mut out, (i + k, j + l), x.mul(y));
        }
    }
    out
}

fn rev_mul<T: Real>(a: &Rev<T>, b: &Rev<T>) -> Rev<T> {
    let mut out = BTreeMap::new();
    for ((k, i, j), x) in a {
        for ((l, m, n), y) in b {
            add_to(&mut out, (k + l, i + m, j + n), x.mul(y));
        }
    }
    out
}

fn sph_mul<T: Real>(a: &Sph<T>, b: &Sph<T>) -> Sph<T> {
    let mut out = BTreeMap::new();
    for ((i, j, k, l), x) in a {
        for ((m, n, o, p), y) in b {
            add_to(&mut out, (i + m, j + n, k + o, l + p), x.mul(y));
        }
    }
    out
}

fn scaled<K: Ord + Copy, T: Real>(a: &BTreeMap<K, T>, s: &T) -> BTreeMap<K, T> {
    a.iter().map(|(k, x)| (*k, x.mul(s))).collect()
}

fn summed<K: Ord + Copy, T: Real>(parts: &[BTreeMap<K, T>]) -> BTreeMap<K, T> {
    let mut out = BTreeMap::new();
    for p in parts {
        for (k, x) in p {
            add_to(&mut out, *k, x.clone());
        }
    }
    out
}

/// The fourteen integrands at a point (or over a box): position `p`
/// relative to the reference, parametric normal `n` and `|n|`, in the order
/// of `integrands`.
pub(super) fn point_integrands<T: Real>(p: &[T; 3], n: &[T; 3], norm: &T) -> [T; TERMS] {
    let third = T::from_r(&R::new(1.into(), 3.into()));
    let half = c::<T>(0.5);
    let sq = |i: usize| p[i].square();
    let volume = p[0]
        .mul(&n[0])
        .add(&p[1].mul(&n[1]))
        .add(&p[2].mul(&n[2]))
        .mul(&third);
    let first = |i: usize| sq(i).mul(&n[i]).mul(&half);
    let second = |i: usize| sq(i).mul(&p[i]).mul(&n[i]).mul(&third);
    let mixed = |i: usize, j: usize| sq(i).mul(&p[j]).mul(&n[i]).mul(&half);
    [
        volume,
        first(0),
        first(1),
        first(2),
        second(0),
        second(1),
        second(2),
        mixed(0, 1),
        mixed(1, 2),
        mixed(2, 0),
        norm.clone(),
        p[0].mul(norm),
        p[1].mul(norm),
        p[2].mul(norm),
    ]
}

/// `sum of c v^k cos^a u sin^b u` over boxes.
fn rev_eval<T: Real>(f: &Rev<T>, u: &T, v: &T) -> T {
    let (co, si) = T::cos_sin(u);
    let power = |x: &T, n: u8| (0..n).fold(c::<T>(1.0), |acc, _| acc.mul(x));
    let mut total = c::<T>(0.0);
    for ((k, a, b), x) in f {
        total = total.add(
            &x.mul(&power(v, *k))
                .mul(&power(&co, *a))
                .mul(&power(&si, *b)),
        );
    }
    total
}

/// The fourteen integrands from the position `p` relative to the reference,
/// the parametric normal `N = S_u x S_v` and `|N|`, all polynomials of one
/// kind.
fn integrands<K: Ord + Copy, T: Real>(
    p: &[BTreeMap<K, T>; 3],
    n: &[BTreeMap<K, T>; 3],
    norm: &BTreeMap<K, T>,
    mul: impl Fn(&BTreeMap<K, T>, &BTreeMap<K, T>) -> BTreeMap<K, T>,
) -> [BTreeMap<K, T>; TERMS] {
    let third = T::from_r(&R::new(1.into(), 3.into()));
    let half = c::<T>(0.5);
    let sq = |i: usize| mul(&p[i], &p[i]);
    let volume = scaled(
        &summed(&[mul(&p[0], &n[0]), mul(&p[1], &n[1]), mul(&p[2], &n[2])]),
        &third,
    );
    let first = |i: usize| scaled(&mul(&sq(i), &n[i]), &half);
    let second = |i: usize| scaled(&mul(&mul(&sq(i), &p[i]), &n[i]), &third);
    // Mixed moments: x^2 y n_x / 2, y^2 z n_y / 2, z^2 x n_z / 2.
    let mixed = |i: usize, j: usize| scaled(&mul(&mul(&sq(i), &p[j]), &n[i]), &half);
    [
        volume,
        first(0),
        first(1),
        first(2),
        second(0),
        second(1),
        second(2),
        mixed(0, 1),
        mixed(1, 2),
        mixed(2, 0),
        norm.clone(),
        mul(&p[0], norm),
        mul(&p[1], norm),
        mul(&p[2], norm),
    ]
}

/// The exact Fourier expansion of `cos^a x sin^b x`: `(f, alpha, beta)` for
/// `alpha cos(f x) + beta sin(f x)`, `f >= 0`. Every coefficient is a sum of
/// at most `2^(a+b)` terms `+-1 / 2^(a+b)`, a dyadic rational with a
/// numerator below `2^(a+b)`, so binary64 holds it exactly for `a + b <= 52`
/// (the integrands here have `a + b <= 8`). Memoized per thread.
/// `(frequency, cos coefficient, sin coefficient)` terms of one expansion.
type Expansion = std::rc::Rc<Vec<(i64, f64, f64)>>;

fn fourier(a: u8, b: u8) -> Expansion {
    assert!(
        u32::from(a) + u32::from(b) <= 52,
        "exact dyadic coefficients"
    );
    thread_local! {
        static MEMO: std::cell::RefCell<BTreeMap<(u8, u8), Expansion>> =
            const { std::cell::RefCell::new(BTreeMap::new()) };
    }
    if let Some(hit) = MEMO.with(|m| m.borrow().get(&(a, b)).cloned()) {
        return hit;
    }
    // Frequency to (cos, sin) coefficients, negative frequencies allowed.
    let mut terms: BTreeMap<i64, (f64, f64)> = BTreeMap::new();
    terms.insert(0, (1.0, 0.0));
    let step = |terms: &BTreeMap<i64, (f64, f64)>, by_sin: bool| {
        let mut out: BTreeMap<i64, (f64, f64)> = BTreeMap::new();
        let mut put = |f: i64, cs: f64, sn: f64| {
            let e = out.entry(f).or_insert((0.0, 0.0));
            e.0 += cs;
            e.1 += sn;
        };
        for (f, (cs, sn)) in terms {
            if by_sin {
                // cos(fx) sin x = [sin((f+1)x) - sin((f-1)x)] / 2
                // sin(fx) sin x = [cos((f-1)x) - cos((f+1)x)] / 2
                put(f + 1, -sn * 0.5, cs * 0.5);
                put(f - 1, sn * 0.5, -(cs * 0.5));
            } else {
                // cos(fx) cos x = [cos((f+1)x) + cos((f-1)x)] / 2
                // sin(fx) cos x = [sin((f+1)x) + sin((f-1)x)] / 2
                put(f + 1, cs * 0.5, sn * 0.5);
                put(f - 1, cs * 0.5, sn * 0.5);
            }
        }
        out
    };
    for _ in 0..a {
        terms = step(&terms, false);
    }
    for _ in 0..b {
        terms = step(&terms, true);
    }
    // Fold negative frequencies: cos(-gx) = cos(gx), sin(-gx) = -sin(gx).
    let mut folded: BTreeMap<i64, (f64, f64)> = BTreeMap::new();
    for (f, (cs, sn)) in terms {
        let (g, sn) = if f < 0 { (-f, -sn) } else { (f, sn) };
        let e = folded.entry(g).or_insert((0.0, 0.0));
        e.0 += cs;
        e.1 += if g == 0 { 0.0 } else { sn };
    }
    let out: Expansion = std::rc::Rc::new(
        folded
            .into_iter()
            .filter(|(_, (cs, sn))| *cs != 0.0 || *sn != 0.0)
            .map(|(f, (cs, sn))| (f, cs, sn))
            .collect(),
    );
    MEMO.with(|m| m.borrow_mut().insert((a, b), out.clone()));
    out
}

/// `C(n, k)`, exact in binary64 for the small `n` here.
fn binomial(n: u32, k: u32) -> f64 {
    let mut out = 1u64;
    for i in 0..u64::from(k) {
        out = out * (u64::from(n) - i) / (i + 1);
    }
    out as f64
}

/// `integral over [0, 1] of the polynomial with these coefficients`.
fn integrate_t<T: Real>(coeffs: &[T]) -> Option<T> {
    let mut total = c::<T>(0.0);
    for (k, x) in coeffs.iter().enumerate() {
        total = total.add(&x.div(&c(k as f64 + 1.0))?);
    }
    Some(total)
}

fn t_mul<T: Real>(a: &[T], b: &[T]) -> Vec<T> {
    let mut out = vec![c::<T>(0.0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = out[i + j].add(&x.mul(y));
        }
    }
    out
}

fn t_pow<T: Real>(base: &[T], n: u8) -> Vec<T> {
    let mut out = vec![c::<T>(1.0)];
    for _ in 0..n {
        out = t_mul(&out, base);
    }
    out
}

/// `-integral of F du` along the line from `a` to `b` on a plane.
fn planar_line<T: Real>(f: &Planar<T>, a: &V2<T>, b: &V2<T>) -> Option<T> {
    let du = b[0].sub(&a[0]);
    let u = [a[0].clone(), du.clone()];
    let v = [a[1].clone(), b[1].sub(&a[1])];
    let mut total = c::<T>(0.0);
    for ((i, j), x) in f {
        let poly = t_mul(&t_pow(&u, *i), &t_pow(&v, *j));
        total = total.add(&x.mul(&integrate_t(&poly)?));
    }
    Some(total.mul(&du).neg())
}

/// `-integral of f du` along a spline pcurve: on each Bézier piece with
/// homogeneous `(U, V, W)`, `f(U/W, V/W) du = F̂ (U' W - U W') / W^(d+2)`
/// with `F̂ = Σ f_ij U^i V^j W^(d-i-j)`, `d` the total degree, enclosed by
/// `bernstein::quotient_integral`.
fn planar_spline<T: Real>(
    f: &Planar<T>,
    curve: &crate::topology::SplineSpan<crate::BSplineCurve2>,
    about: &[R; 2],
) -> Option<T> {
    use super::bernstein::{
        derivative, difference, power, product, quotient_integral, scaled, span_arcs, sum,
    };
    let d = f.keys().map(|(i, j)| usize::from(i + j)).max().unwrap_or(0);
    let mut total = c::<T>(0.0);
    for (_, _, arc) in &span_arcs(curve)? {
        let poles = arc.homogeneous_poles();
        let uniform = poles.iter().all(|p| p[3] == poles[0][3]);
        // Translated by `-about`, exactly: `U - a W`, `V - b W`.
        let [u, v, w]: [Vec<T>; 3] = [0, 1, 3].map(|k| {
            poles
                .iter()
                .map(|p| match k {
                    0 | 1 => T::from_r(&(&p[k] - &about[k] * &p[3])),
                    _ => T::from_r(&p[k]),
                })
                .collect()
        });
        let mut hat = vec![c::<T>(0.0)];
        for ((i, j), x) in f {
            let (i, j) = (usize::from(*i), usize::from(*j));
            let term = product(
                &product(&power(&u, i), &power(&v, j)),
                &power(&w, d - i - j),
            );
            hat = sum(&hat, &scaled(&term, x));
        }
        let du = difference(&product(&derivative(&u), &w), &product(&u, &derivative(&w)));
        total = total.add(&quotient_integral(
            &product(&hat, &du),
            &w,
            (d + 2) as i32,
            uniform,
        )?);
    }
    Some(total.neg())
}

/// `integral over [t0, t0 + sweep] of cos^p t sin^q t dt`, exactly.
fn trig_integral<T: Real>(cos_power: u8, sin_power: u8, t0: &T, t1: &T) -> Option<T> {
    let mut total = c::<T>(0.0);
    for &(f, alpha, beta) in fourier(cos_power, sin_power).iter() {
        if f == 0 {
            total = total.add(&c::<T>(alpha).mul(&t1.sub(t0)));
            continue;
        }
        let g = c::<T>(f as f64);
        let (c0, s0) = T::cos_sin(&t0.mul(&g));
        let (c1, s1) = T::cos_sin(&t1.mul(&g));
        let cos_part = c::<T>(alpha).mul(&s1.sub(&s0));
        let sin_part = c::<T>(beta).mul(&c0.sub(&c1));
        total = total.add(&cos_part.add(&sin_part).div(&g)?);
    }
    Some(total)
}

/// `-integral of F du` along a plane arc or axis-aligned ellipse arc:
/// `u = cx + rx cos t`, `v = cy + ry sin t`, so `-F du = F rx sin t dt`.
fn planar_arc<T: Real>(
    f: &Planar<T>,
    center: [T; 2],
    [rx, ry]: [f64; 2],
    start: f64,
    sweep: f64,
) -> Option<T> {
    let rad = c::<T>(rx);
    // Polynomials in (cos t, sin t).
    let [cu, cv] = center;
    let u: Planar<T> = [((0, 0), cu), ((1, 0), rad.clone())].into_iter().collect();
    let v: Planar<T> = [((0, 0), cv), ((0, 1), c(ry))].into_iter().collect();
    let pow = |base: &Planar<T>, n: u8| {
        let mut out: Planar<T> = [((0, 0), c(1.0))].into_iter().collect();
        for _ in 0..n {
            out = planar_mul(&out, base);
        }
        out
    };
    let t0 = c::<T>(start);
    let t1 = t0.add(&c(sweep));
    let mut total = c::<T>(0.0);
    for ((i, j), x) in f {
        for ((p, qq), y) in planar_mul(&pow(&u, *i), &pow(&v, *j)) {
            let term = trig_integral::<T>(p, qq + 1, &t0, &t1)?;
            total = total.add(&x.mul(&y).mul(&rad).mul(&term));
        }
    }
    Some(total)
}

/// `-integral of F du` along `v = a0 + a1 cos u + a2 sin u` from `start`
/// over `sweep` on a cylinder, for every `F` of `fs`: each monomial
/// `v^k cos^p u sin^q u` is a polynomial in `(cos u, sin u)` integrated
/// exactly.
fn rev_sinusoid<T: Real>(fs: &[Rev<T>], start: f64, sweep: f64, a: &[f64; 3]) -> Option<Vec<T>> {
    let u0 = c::<T>(start);
    let u1 = u0.add(&c(sweep));
    let v: Planar<T> = [((0, 0), c(a[0])), ((1, 0), c(a[1])), ((0, 1), c(a[2]))]
        .into_iter()
        .collect();
    let mut powers: Vec<Planar<T>> = vec![[((0, 0), c(1.0))].into_iter().collect()];
    let mut integrals: BTreeMap<(u8, u8), T> = BTreeMap::new();
    let mut out = Vec::with_capacity(fs.len());
    for f in fs {
        let mut total = c::<T>(0.0);
        for ((k, p, qq), x) in f {
            while powers.len() <= usize::from(*k) {
                let next = planar_mul(powers.last().expect("one"), &v);
                powers.push(next);
            }
            for ((i, j), y) in &powers[usize::from(*k)] {
                let key = (i + p, j + qq);
                let value = match integrals.get(&key) {
                    Some(value) => value.clone(),
                    None => {
                        let value = trig_integral::<T>(key.0, key.1, &u0, &u1)?;
                        integrals.insert(key, value.clone());
                        value
                    }
                };
                total = total.add(&x.mul(y).mul(&value));
            }
        }
        out.push(total.neg());
    }
    Some(out)
}

/// `-integral of F du` along the line from `a` to `b` on a surface of
/// revolution, for every `F` of `fs` at once (they share the line's
/// trigonometric values and moments); `None` when `du` may be zero without
/// being zero.
fn rev_lines<T: Real>(fs: &[Rev<T>], a: &V2<T>, b: &V2<T>) -> Option<Vec<T>> {
    let d = b[0].sub(&a[0]);
    if d.sign()? == Ordering::Equal {
        return Some(vec![c(0.0); fs.len()]);
    }
    let m = b[1].sub(&a[1]).div(&d)?;
    let mut ends: BTreeMap<i64, (T, T, T, T)> = BTreeMap::new();
    let mut cache: BTreeMap<(i64, u32), (T, T)> = BTreeMap::new();
    // J(f, j) = (integral over [0, d] of w^j cos(f(u0 + w)), ... sin ...),
    // by the recursion on j from the endpoint values.
    let mut moments = |f: i64, j: u32| -> Option<(T, T)> {
        if let Some(x) = cache.get(&(f, j)) {
            return Some(x.clone());
        }
        let out = if f == 0 {
            let power = (0..=j).fold(c::<T>(1.0), |p, _| p.mul(&d));
            (power.div(&c(f64::from(j) + 1.0))?, c(0.0))
        } else {
            let g = c::<T>(f as f64);
            let (c0, s0, c1, s1) = match ends.get(&f) {
                Some(x) => x.clone(),
                None => {
                    let (c0, s0) = T::cos_sin(&a[0].mul(&g));
                    let (c1, s1) = T::cos_sin(&b[0].mul(&g));
                    ends.insert(f, (c0.clone(), s0.clone(), c1.clone(), s1.clone()));
                    (c0, s0, c1, s1)
                }
            };
            let mut cs = s1.sub(&s0).div(&g)?;
            let mut sn = c0.sub(&c1).div(&g)?;
            let mut power = c::<T>(1.0);
            for k in 1..=j {
                power = power.mul(&d);
                let k_over = c::<T>(f64::from(k)).div(&g)?;
                let next_cs = power.mul(&s1).div(&g)?.sub(&k_over.mul(&sn));
                let next_sn = power.mul(&c1).div(&g)?.neg().add(&k_over.mul(&cs));
                cs = next_cs;
                sn = next_sn;
            }
            (cs, sn)
        };
        cache.insert((f, j), out.clone());
        Some(out)
    };
    // v^k = (v0 + m w)^k = sum over j of C(k, j) v0^(k-j) m^j w^j.
    let mut v0_powers = vec![c::<T>(1.0)];
    let mut m_powers = vec![c::<T>(1.0)];
    let mut out = Vec::with_capacity(fs.len());
    for f in fs {
        let mut total = c::<T>(0.0);
        for ((k, p, qq), x) in f {
            let k = usize::from(*k);
            while v0_powers.len() <= k {
                let next = v0_powers.last().expect("one").mul(&a[1]);
                v0_powers.push(next);
                let next = m_powers.last().expect("one").mul(&m);
                m_powers.push(next);
            }
            let expansion = fourier(*p, *qq);
            for j in 0..=k {
                let coef = c::<T>(binomial(k as u32, j as u32))
                    .mul(&v0_powers[k - j])
                    .mul(&m_powers[j]);
                let mut term = c::<T>(0.0);
                for &(freq, alpha, beta) in expansion.iter() {
                    let (jc, js) = moments(freq, j as u32)?;
                    term = term.add(&c::<T>(alpha).mul(&jc).add(&c::<T>(beta).mul(&js)));
                }
                total = total.add(&x.mul(&coef).mul(&term));
            }
        }
        out.push(total.neg());
    }
    Some(out)
}

/// An enclosure of every point of the segment between two enclosed values.
fn hull<T: Real>(a: &T, b: &T) -> T {
    let mid = a.add(b).mul(&c(0.5));
    let half = b.sub(a).mul(&c(0.5));
    let mid_half = half.midpoint();
    let magnitude = if mid_half < int(0) {
        -mid_half
    } else {
        mid_half
    };
    let reach = magnitude + half.radius();
    mid.widen(&reach)
}

/// `-integral of F du` for `u` from `u0` to `u1` at a fixed `v0`, for every
/// `F(u, v) = integral from lower to v of f dv` of `fs`: `F(u, v0)` is a
/// trigonometric polynomial in `u`, integrated exactly.
fn sph_parallel<T: Real>(fs: &[Sph<T>], u0: &T, u1: &T, v0: &T, lower: &T) -> Option<Vec<T>> {
    let mut in_v: BTreeMap<(u8, u8), T> = BTreeMap::new();
    let mut in_u: BTreeMap<(u8, u8), T> = BTreeMap::new();
    let mut out = Vec::with_capacity(fs.len());
    for f in fs {
        let mut collapsed: BTreeMap<(u8, u8), T> = BTreeMap::new();
        for ((a, b, cc, d), x) in f {
            let tv = match in_v.get(&(*cc, *d)) {
                Some(t) => t.clone(),
                None => {
                    let t = trig_integral::<T>(*cc, *d, lower, v0)?;
                    in_v.insert((*cc, *d), t.clone());
                    t
                }
            };
            add_to(&mut collapsed, (*a, *b), x.mul(&tv));
        }
        let mut total = c::<T>(0.0);
        for ((a, b), y) in collapsed {
            let tu = match in_u.get(&(a, b)) {
                Some(t) => t.clone(),
                None => {
                    let t = trig_integral::<T>(a, b, u0, u1)?;
                    in_u.insert((a, b), t.clone());
                    t
                }
            };
            total = total.add(&y.mul(&tu));
        }
        out.push(total.neg());
    }
    Some(out)
}

/// `-integral of F du` along the line from `a` to `b` on a sphere, for every
/// `F(u, v) = integral from lower to v of f dv` of `fs`. A meridian
/// (`du = 0`) contributes nothing and a parallel (`dv = 0`) is exact; any
/// other line (a chord closing a gap) is enclosed by `-du` times `F` over the
/// segment's bounding box, which contains the mean value of `F`.
fn sph_lines<T: Real>(fs: &[Sph<T>], a: &V2<T>, b: &V2<T>, lower: &T) -> Option<Vec<T>> {
    let du = b[0].sub(&a[0]);
    if du.sign()? == Ordering::Equal {
        return Some(vec![c(0.0); fs.len()]);
    }
    let dv = b[1].sub(&a[1]);
    if dv.sign() == Some(Ordering::Equal) {
        return sph_parallel(fs, &a[0], &b[0], &a[1], lower);
    }
    let (uu, vv) = (hull(&a[0], &b[0]), hull(&a[1], &b[1]));
    let (cu, su) = T::cos_sin(&uu);
    let power = |x: &T, n: u8| (0..n).fold(c::<T>(1.0), |acc, _| acc.mul(x));
    let mut out = Vec::with_capacity(fs.len());
    for f in fs {
        let mut total = c::<T>(0.0);
        for ((i, j, k, l), x) in f {
            let tv = trig_integral::<T>(*k, *l, lower, &vv)?;
            total = total.add(&x.mul(&power(&cu, *i)).mul(&power(&su, *j)).mul(&tv));
        }
        out.push(total.mul(&du).neg());
    }
    Some(out)
}

/// The sphere's integrands from `p = d + R q` and `N = R^2 cos v q`, `q` the
/// unit radial direction, and `|N| = R^2 cos v`.
fn sphere_terms<T: Real>(fr: &FrameV<T>, radius: f64, d: &V3<T>) -> [Sph<T>; TERMS] {
    let rad = c::<T>(radius);
    let r2 = rad.square();
    let p: [Sph<T>; 3] = std::array::from_fn(|i| {
        [
            ((0, 0, 0, 0), d[i].clone()),
            ((1, 0, 1, 0), rad.mul(&fr.x[i])),
            ((0, 1, 1, 0), rad.mul(&fr.y[i])),
            ((0, 0, 0, 1), rad.mul(&fr.n[i])),
        ]
        .into_iter()
        .collect()
    });
    let n: [Sph<T>; 3] = std::array::from_fn(|i| {
        [
            ((1, 0, 2, 0), r2.mul(&fr.x[i])),
            ((0, 1, 2, 0), r2.mul(&fr.y[i])),
            ((0, 0, 1, 1), r2.mul(&fr.n[i])),
        ]
        .into_iter()
        .collect()
    });
    let norm: Sph<T> = [((0, 0, 1, 0), r2)].into_iter().collect();
    integrands(&p, &n, &norm, sph_mul)
}

/// The torus's integrands from `p = d + (R + r cos v) e(u) + r sin v n`,
/// `N = r (R + r cos v) q` (`q = cos v e(u) + sin v n`, outward from the
/// tube) and `|N| = r (R + r cos v)`.
fn torus_terms<T: Real>(fr: &FrameV<T>, major: f64, minor: f64, d: &V3<T>) -> [Sph<T>; TERMS] {
    let (big, r) = (c::<T>(major), c::<T>(minor));
    let (rr, rbig) = (r.square(), r.mul(&big));
    let p: [Sph<T>; 3] = std::array::from_fn(|i| {
        [
            ((0, 0, 0, 0), d[i].clone()),
            ((1, 0, 0, 0), big.mul(&fr.x[i])),
            ((1, 0, 1, 0), r.mul(&fr.x[i])),
            ((0, 1, 0, 0), big.mul(&fr.y[i])),
            ((0, 1, 1, 0), r.mul(&fr.y[i])),
            ((0, 0, 0, 1), r.mul(&fr.n[i])),
        ]
        .into_iter()
        .collect()
    });
    let n: [Sph<T>; 3] = std::array::from_fn(|i| {
        [
            ((1, 0, 1, 0), rbig.mul(&fr.x[i])),
            ((1, 0, 2, 0), rr.mul(&fr.x[i])),
            ((0, 1, 1, 0), rbig.mul(&fr.y[i])),
            ((0, 1, 2, 0), rr.mul(&fr.y[i])),
            ((0, 0, 0, 1), rbig.mul(&fr.n[i])),
            ((0, 0, 1, 1), rr.mul(&fr.n[i])),
        ]
        .into_iter()
        .collect()
    });
    let norm: Sph<T> = [((0, 0, 0, 0), rbig), ((0, 0, 1, 0), rr)]
        .into_iter()
        .collect();
    integrands(&p, &n, &norm, sph_mul)
}

/// The integral of every integrand of `fs` over a sphere or torus face (its
/// loops, closed by chords, carry the orientation). Loops winding in `u` or
/// none: `-loop integral of F du`, `F` from the pole on a face with one, the
/// south pole on another sphere face, `v = 0` on a torus. Loops winding in
/// `v` (a torus): `loop integral of G dv`, `G` from `u = 0`: the same
/// routine with `u` and `v` exchanged (a meridian then integrates exactly, a
/// parallel contributes nothing). No edge loops: the whole surface, bounded
/// on the cover by the north pole's line (a sphere) or by both periods (a
/// torus).
fn trig_face<T: Real>(face: &Face, loops: &[Lp], fs: &[Sph<T>]) -> Option<Vec<T>> {
    let torus = matches!(face.surface, Surface::Torus { .. });
    let turns: i32 = loops.iter().map(|lp| lp.winding).sum();
    let wound_u = loops.iter().any(|lp| lp.winding != 0);
    let wound_v = loops.iter().any(|lp| lp.winding_v != 0);
    if wound_u && wound_v {
        return None;
    }
    let two_pi = T::from_r(&(pi().midpoint() * int(2))).widen(&(pi().radius() * int(2)));
    let sign = if face.sense == Orientation::Forward {
        1.0
    } else {
        -1.0
    };
    let mut totals = vec![c::<T>(0.0); fs.len()];
    let mut accumulate = |values: Vec<T>| {
        for (t, v) in totals.iter_mut().zip(values) {
            *t = t.add(&v);
        }
    };
    if loops.is_empty() {
        let whole = if torus {
            sph_parallel(fs, &two_pi, &c(0.0), &two_pi, &c(0.0))?
        } else {
            sph_parallel(fs, &two_pi, &c(0.0), &half_pi(), &half_pi::<T>().neg())?
        };
        accumulate(whole.iter().map(|x| x.mul(&c(sign))).collect());
        return Some(totals);
    }
    if wound_v {
        // Exchange u and v: keys (a, b, c, d) -> (c, d, a, b), points
        // (u, v) -> (v, u); the routine then gives -loop integral of G dv.
        let swapped: Vec<Sph<T>> = fs
            .iter()
            .map(|f| {
                f.iter()
                    .map(|((a, b, cc, d), x)| ((*cc, *d, *a, *b), x.clone()))
                    .collect()
            })
            .collect();
        let flip = |p: &V2<T>| -> V2<T> { [p[1].clone(), p[0].clone()] };
        for lp in loops {
            for u in &lp.fins {
                let Curve2::LineSegment { start, end } = &u.pcurve else {
                    return None;
                };
                let values = sph_lines(
                    &swapped,
                    &[c(start.y), c(start.x)],
                    &[c(end.y), c(end.x)],
                    &c(0.0),
                )?;
                accumulate(values.iter().map(|x| x.neg()).collect());
            }
            for (a, b) in chords::<T>(lp) {
                let values = sph_lines(&swapped, &flip(&a), &flip(&b), &c(0.0))?;
                accumulate(values.iter().map(|x| x.neg()).collect());
            }
        }
        return Some(totals);
    }
    let lower = if torus {
        c(0.0)
    } else if turns.abs() == 1 {
        pole_v::<T>(&face.surface, pole_north(turns, face.sense))?
    } else {
        half_pi::<T>().neg()
    };
    for lp in loops {
        for u in &lp.fins {
            let (start, end) = match &u.pcurve {
                Curve2::LineSegment { start, end } => (start, end),
                // -∫ F(u, v) du along a spline, F each integrand's
                // antiderivative in v from `lower`, enclosed (S4d).
                Curve2::BSpline(spline) => {
                    let mut values = Vec::with_capacity(fs.len());
                    for f in fs {
                        let g = |uu: &T, v: &T| sph_antiderivative(f, uu, v, &lower);
                        values.push(
                            super::bernstein::green_integral(spline, super::SPLINE_DEPTH, &g)?
                                .neg(),
                        );
                    }
                    accumulate(values);
                    continue;
                }
                Curve2::CircularArc { .. }
                | Curve2::EllipseArc { .. }
                | Curve2::Sinusoid { .. } => return None,
            };
            accumulate(sph_lines(
                fs,
                &[c(start.x), c(start.y)],
                &[c(end.x), c(end.y)],
                &lower,
            )?);
        }
        for (a, b) in chords::<T>(lp) {
            let values = match sph_lines(fs, &a, &b, &lower) {
                Some(values) => values,
                None => fs
                    .iter()
                    .map(|f| {
                        let g = |uu: &T, v: &T| sph_antiderivative(f, uu, v, &lower);
                        super::chord_enclosure(&a, &b, &g)
                    })
                    .collect::<Option<_>>()?,
            };
            accumulate(values);
        }
    }
    Some(totals)
}

/// `F(u, v) = ∫_lower^v f(u, s) ds` of a trigonometric polynomial, over
/// boxes: each monomial's `cos^a u sin^b u` times the exact integral of
/// `cos^c sin^d` from `lower` to `v`.
fn sph_antiderivative<T: Real>(f: &Sph<T>, u: &T, v: &T, lower: &T) -> Option<T> {
    let (co, si) = T::cos_sin(u);
    let power = |x: &T, n: u8| (0..n).fold(c::<T>(1.0), |acc, _| acc.mul(x));
    let mut total = c::<T>(0.0);
    for ((a, b, cc, d), x) in f {
        let along_u = power(&co, *a).mul(&power(&si, *b));
        total = total.add(&x.mul(&along_u).mul(&trig_integral(*cc, *d, lower, v)?));
    }
    Some(total)
}

/// The orientation flux of a sphere or torus face: the integral of
/// `(S - origin).(S_u x S_v)` over it (as `face_flux`).
pub(super) fn sphere_flux<T: Real>(face: &Face, loops: &[Lp], origin: &V3<T>) -> Option<T> {
    let terms = match &face.surface {
        Surface::Sphere { frame: f, radius } => {
            let fr = frame::<T>(f);
            sphere_terms(&fr, *radius, &vsub(&fr.o, origin))
        }
        Surface::Torus {
            frame: f,
            major,
            minor,
        } => {
            let fr = frame::<T>(f);
            torus_terms(&fr, *major, *minor, &vsub(&fr.o, origin))
        }
        _ => return None,
    };
    // S.N = 3 times the volume integrand.
    let flux = scaled(&terms[0], &c(3.0));
    trig_face(face, loops, &[flux])?.pop()
}

/// The fourteen face integrals over the face region (loops carry its
/// orientation, as in `face_flux`), relative to `reference`.
fn face_integrals<T: Real>(face: &Face, loops: &[Lp], reference: &V3<T>) -> Option<[T; TERMS]> {
    let mut totals: [T; TERMS] = std::array::from_fn(|_| c(0.0));
    let mut accumulate = |values: [T; TERMS]| {
        for (t, v) in totals.iter_mut().zip(values) {
            *t = t.add(&v);
        }
    };
    match &face.surface {
        Surface::Plane(f) => {
            let fr = frame::<T>(f);
            // In UV relative to a point of the face (exact, `(u0, v0)`),
            // so a face far from its frame's origin keeps its enclosures at
            // its own scale: `p = O + u0 X + v0 Y - reference + u' X + v' Y`.
            let o: [R; 2] = match loops.first().and_then(|lp| lp.fins.first()) {
                Some(u) => pcurve_at::<T>(&u.pcurve, 0.0).map(|x| x.midpoint()),
                None => [R::from_integer(0.into()), R::from_integer(0.into())],
            };
            let (ou, ov) = (T::from_r(&o[0]), T::from_r(&o[1]));
            let d: V3<T> = std::array::from_fn(|i| {
                fr.o[i]
                    .sub(&reference[i])
                    .add(&ou.mul(&fr.x[i]))
                    .add(&ov.mul(&fr.y[i]))
            });
            let at = |x: f64, y: f64| [c::<T>(x).sub(&ou), c::<T>(y).sub(&ov)];
            let p: [Planar<T>; 3] = std::array::from_fn(|i| {
                [
                    ((0, 0), d[i].clone()),
                    ((1, 0), fr.x[i].clone()),
                    ((0, 1), fr.y[i].clone()),
                ]
                .into_iter()
                .collect()
            });
            let n: [Planar<T>; 3] =
                std::array::from_fn(|i| [((0, 0), fr.n[i].clone())].into_iter().collect());
            let one: Planar<T> = [((0, 0), c(1.0))].into_iter().collect();
            // The v-antiderivative from 0.
            let anti: Vec<Planar<T>> = integrands(&p, &n, &one, planar_mul)
                .iter()
                .map(|f| {
                    f.iter()
                        .map(|((a, b), x)| Some(((*a, b + 1), x.div(&c(f64::from(*b) + 1.0))?)))
                        .collect::<Option<Planar<T>>>()
                })
                .collect::<Option<_>>()?;
            for lp in loops {
                for u in &lp.fins {
                    let values: [Option<T>; TERMS] = std::array::from_fn(|k| match &u.pcurve {
                        Curve2::LineSegment { start, end } => {
                            planar_line(&anti[k], &at(start.x, start.y), &at(end.x, end.y))
                        }
                        Curve2::CircularArc {
                            center,
                            radius,
                            start_angle,
                            sweep_angle,
                        } => planar_arc(
                            &anti[k],
                            at(center.x, center.y),
                            [*radius; 2],
                            *start_angle,
                            *sweep_angle,
                        ),
                        Curve2::EllipseArc {
                            center,
                            major,
                            minor,
                            start_angle,
                            sweep_angle,
                        } => planar_arc(
                            &anti[k],
                            at(center.x, center.y),
                            [*major, *minor],
                            *start_angle,
                            *sweep_angle,
                        ),
                        Curve2::BSpline(spline) => planar_spline(&anti[k], spline, &o),
                        Curve2::Sinusoid { .. } => None,
                    });
                    accumulate(values.try_map_all()?);
                }
                for (a, b) in chords::<T>(lp) {
                    let shift = |x: &V2<T>| [x[0].sub(&ou), x[1].sub(&ov)];
                    let (a, b) = (shift(&a), shift(&b));
                    let values: [Option<T>; TERMS] =
                        std::array::from_fn(|k| planar_line(&anti[k], &a, &b));
                    accumulate(values.try_map_all()?);
                }
            }
        }
        Surface::Sphere { frame: f, radius } => {
            let fr = frame::<T>(f);
            let d = vsub(&fr.o, reference);
            let terms = sphere_terms(&fr, *radius, &d);
            accumulate(trig_face(face, loops, &terms)?.try_into().ok()?);
        }
        Surface::Torus {
            frame: f,
            major,
            minor,
        } => {
            let fr = frame::<T>(f);
            let d = vsub(&fr.o, reference);
            let terms = torus_terms(&fr, *major, *minor, &d);
            accumulate(trig_face(face, loops, &terms)?.try_into().ok()?);
        }
        Surface::BSpline(surface) => {
            let values = super::spline_flux::spline_face_integrals(surface, loops, reference)?;
            accumulate(values.try_into().ok()?);
        }
        Surface::Cylinder { frame: f, radius }
        | Surface::Cone {
            frame: f, radius, ..
        } => {
            let fr = frame::<T>(f);
            let (cth, sth) = match &face.surface {
                Surface::Cone { half_angle, .. } => T::cos_sin(&c(*half_angle)),
                _ => (c(1.0), c(0.0)),
            };
            let d = vsub(&fr.o, reference);
            let rad = c::<T>(*radius);
            // rho(v) = R + v sin a, e(u) = x cos u + y sin u, h(v) = v cos a.
            let rho: Rev<T> = [((0, 0, 0), rad.clone()), ((1, 0, 0), sth.clone())]
                .into_iter()
                .collect();
            let e: [Rev<T>; 3] = std::array::from_fn(|i| {
                [((0, 1, 0), fr.x[i].clone()), ((0, 0, 1), fr.y[i].clone())]
                    .into_iter()
                    .collect()
            });
            let p: [Rev<T>; 3] = std::array::from_fn(|i| {
                let mut out = rev_mul(&rho, &e[i]);
                add_to(&mut out, (0, 0, 0), d[i].clone());
                add_to(&mut out, (1, 0, 0), cth.mul(&fr.n[i]));
                out
            });
            // N = rho (e cos a - n sin a).
            let n: [Rev<T>; 3] = std::array::from_fn(|i| {
                let mut dir = scaled(&e[i], &cth);
                add_to(&mut dir, (0, 0, 0), sth.mul(&fr.n[i]).neg());
                rev_mul(&rho, &dir)
            });
            // The v-antiderivative from the apex on a face with a pole (zero
            // there), from 0 otherwise: any constant cancels over a band's
            // loops, and a far apex would cancel badly.
            let poled = matches!(face.surface, Surface::Cone { .. })
                && loops.iter().map(|lp| lp.winding).sum::<i32>().abs() == 1;
            let lower = if poled {
                apex_v::<T>(&face.surface)?
            } else {
                c(0.0)
            };
            let anti: Vec<Rev<T>> = integrands(&p, &n, &rho, rev_mul)
                .iter()
                .map(|f| {
                    let mut out: Rev<T> = BTreeMap::new();
                    for ((k, a, b), x) in f {
                        let scale = c::<T>(f64::from(*k) + 1.0);
                        add_to(&mut out, (k + 1, *a, *b), x.div(&scale)?);
                        let base = (0..=*k).fold(c::<T>(1.0), |acc, _| acc.mul(&lower));
                        add_to(&mut out, (0, *a, *b), x.mul(&base).div(&scale)?.neg());
                    }
                    Some(out)
                })
                .collect::<Option<_>>()?;
            for lp in loops {
                for u in &lp.fins {
                    let values = match &u.pcurve {
                        Curve2::LineSegment { start, end } => {
                            rev_lines(&anti, &[c(start.x), c(start.y)], &[c(end.x), c(end.y)])?
                        }
                        // -∫ F(u, v) du along a spline, enclosed (S4d).
                        Curve2::BSpline(spline) => {
                            let g = |uu: &T, v: &T| {
                                Some(anti.iter().map(|f| rev_eval(f, uu, v)).collect())
                            };
                            super::bernstein::green_integrals(spline, super::SPLINE_DEPTH, &g)?
                                .iter()
                                .map(|x| x.neg())
                                .collect()
                        }
                        Curve2::Sinusoid { start, sweep, a } => {
                            rev_sinusoid(&anti, *start, *sweep, a)?
                        }
                        Curve2::CircularArc { .. } | Curve2::EllipseArc { .. } => return None,
                    };
                    accumulate(values.try_into().ok()?);
                }
                for (a, b) in chords::<T>(lp) {
                    // A closing chord whose `du` straddles zero (a loop that
                    // ends within rounding of its start a turn on) is
                    // enclosed over its box.
                    let values = match rev_lines(&anti, &a, &b) {
                        Some(values) => values,
                        None => anti
                            .iter()
                            .map(|f| {
                                super::chord_enclosure(&a, &b, &|u: &T, v: &T| {
                                    Some(rev_eval(f, u, v))
                                })
                            })
                            .collect::<Option<Vec<T>>>()?,
                    };
                    accumulate(values.try_into().ok()?);
                }
            }
        }
    }
    Some(totals)
}

trait TryMapAll<T> {
    fn try_map_all(self) -> Option<[T; TERMS]>;
}
impl<T> TryMapAll<T> for [Option<T>; TERMS] {
    fn try_map_all(self) -> Option<[T; TERMS]> {
        let mut out = Vec::with_capacity(TERMS);
        for v in self {
            out.push(v?);
        }
        out.try_into().ok()
    }
}

fn resolved<'a>(view: &View<'a>) -> Vec<Vec<Lp<'a>>> {
    view.faces
        .iter()
        .map(|face| {
            face.loops
                .iter()
                .filter_map(|l| match &view.loops[l.0] {
                    Loop::Edges {
                        fins: list,
                        winding,
                    } => Some(Lp {
                        fins: list.iter().map(|k| &view.fins[k.0]).collect(),
                        winding: winding[0],
                        winding_v: winding[1],
                    }),
                    Loop::Vertex(_) => None,
                })
                .collect()
        })
        .collect()
}

/// A certified enclosure `[lo, hi]` of each property of the solid regions.
pub(crate) struct Enclosed {
    pub volume: [f64; 2],
    pub area: [f64; 2],
    pub centroid: [[f64; 2]; 3],
    pub inertia: [[[f64; 2]; 3]; 3],
}

fn solve<T: Real>(view: &View, reference: [f64; 3]) -> Option<Enclosed> {
    let origin = v3::<T>(reference);
    let loops = resolved(view);
    let mut region: [T; 10] = std::array::from_fn(|_| c(0.0));
    let mut area = c::<T>(0.0);
    let mut counted = BTreeSet::new();
    for r in view.regions.iter().filter(|r| r.kind == RegionKind::Solid) {
        for s in &r.shells {
            for (f, side) in &view.shells[s.0].sides {
                let values = face_integrals::<T>(&view.faces[f.0], &loops[f.0], &origin)?;
                for (t, v) in region.iter_mut().zip(values.iter()) {
                    *t = match side {
                        Side::Front => t.add(v),
                        Side::Back => t.sub(v),
                    };
                }
                if counted.insert(f.0) {
                    // The face's area regardless of how its loops run.
                    area = match values[10].sign()? {
                        Ordering::Less => area.sub(&values[10]),
                        _ => area.add(&values[10]),
                    };
                }
            }
        }
    }
    let [volume, mx, my, mz, sxx, syy, szz, sxy, syz, sxz] = region;
    let centre = [mx.div(&volume)?, my.div(&volume)?, mz.div(&volume)?];
    // Inertia about the reference, then about the centroid: J = integral of
    // p p^T, I = tr(J) 1 - J, shifted by V (|c|^2 1 - c c^T).
    let j = [
        [sxx.clone(), sxy.clone(), sxz.clone()],
        [sxy, syy.clone(), syz.clone()],
        [sxz, syz, szz],
    ];
    let trace = j[0][0].add(&j[1][1]).add(&j[2][2]);
    let c2 = vdot(&centre, &centre);
    let inertia: [[T; 3]; 3] = std::array::from_fn(|a| {
        std::array::from_fn(|b| {
            let delta = if a == b { c::<T>(1.0) } else { c::<T>(0.0) };
            let about_reference = delta.mul(&trace).sub(&j[a][b]);
            let shift = delta.mul(&c2).sub(&centre[a].mul(&centre[b])).mul(&volume);
            about_reference.sub(&shift)
        })
    });
    let bounds = |x: &T| {
        let (lo, hi) = x.bounds_f64();
        (lo.is_finite() && hi.is_finite()).then_some([lo, hi])
    };
    let world: [T; 3] = std::array::from_fn(|i| centre[i].add(&origin[i]));
    Some(Enclosed {
        volume: bounds(&volume)?,
        area: bounds(&area)?,
        centroid: [bounds(&world[0])?, bounds(&world[1])?, bounds(&world[2])?],
        inertia: {
            let mut out = [[[0.0; 2]; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    out[a][b] = bounds(&inertia[a][b])?;
                }
            }
            out
        },
    })
}

/// Mass properties of every solid region: binary64 intervals first, rational
/// intervals when those cannot bound a result.
pub(crate) fn mass(view: &View, reference: [f64; 3]) -> Option<Enclosed> {
    super::bernstein::clear_memo();
    solve::<Fast>(view, reference).or_else(|| solve::<I>(view, reference))
}

/// A face's area and centre of gravity, certified.
pub(crate) fn face_mass(
    view: &View,
    face: usize,
    reference: [f64; 3],
) -> Option<([f64; 2], [[f64; 2]; 3])> {
    fn one<T: Real>(
        view: &View,
        face: usize,
        reference: [f64; 3],
    ) -> Option<([f64; 2], [[f64; 2]; 3])> {
        let origin = v3::<T>(reference);
        let loops = resolved(view);
        let values = face_integrals::<T>(&view.faces[face], &loops[face], &origin)?;
        let area = values[10].clone();
        let bounds = |x: &T| {
            let (lo, hi) = x.bounds_f64();
            (lo.is_finite() && hi.is_finite()).then_some([lo, hi])
        };
        let centre: [T; 3] = [
            values[11].div(&area)?.add(&origin[0]),
            values[12].div(&area)?.add(&origin[1]),
            values[13].div(&area)?.add(&origin[2]),
        ];
        let magnitude = match area.sign()? {
            Ordering::Less => area.neg(),
            _ => area,
        };
        Some((
            bounds(&magnitude)?,
            [
                bounds(&centre[0])?,
                bounds(&centre[1])?,
                bounds(&centre[2])?,
            ],
        ))
    }
    one::<Fast>(view, face, reference).or_else(|| one::<I>(view, face, reference))
}

/// A certified enclosure of a sheet's area or a wire's length and its
/// centre (S6).
pub(crate) struct SheetMeasure {
    pub measure: [f64; 2],
    pub centre: [[f64; 2]; 3],
}

/// The length of an edge's curve and its first moment relative to
/// `reference`, in closed form: a segment's midpoint, an arc's
/// `r s (o sweep + r ((sin a1 - sin a0) x + (cos a0 - cos a1) y))`, `s`
/// the sweep's sign. `None` for a spline curve.
fn curve_moments<T: Real>(curve: &Curve3, reference: &V3<T>) -> Option<(T, V3<T>)> {
    match curve {
        Curve3::LineSegment { start, end } => {
            let (a, b) = (v3::<T>(start.to_array()), v3::<T>(end.to_array()));
            let length = vdot(&vsub(&b, &a), &vsub(&b, &a)).sqrt();
            let half = c::<T>(0.5);
            let mid: V3<T> = std::array::from_fn(|i| a[i].add(&b[i]).mul(&half).sub(&reference[i]));
            Some((length.clone(), mid.map(|m| m.mul(&length))))
        }
        Curve3::Circle { frame: f, radius } => {
            let fr = frame::<T>(f);
            let pi = T::from_r(crate::certified::pi().lo())
                .union(&T::from_r(crate::certified::pi().hi()));
            let length = c::<T>(2.0).mul(&pi).mul(&c(*radius));
            let d = vsub(&fr.o, reference);
            Some((length.clone(), d.map(|x| x.mul(&length))))
        }
        Curve3::CircularArc {
            frame: f,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let fr = frame::<T>(f);
            let r = c::<T>(*radius);
            let a0 = c::<T>(*start_angle);
            let sweep = c::<T>(*sweep_angle);
            let (c0, s0) = T::cos_sin(&a0);
            let (c1, s1) = T::cos_sin(&a0.add(&sweep));
            let sign = c::<T>(sweep_angle.signum());
            let length = r.mul(&sweep).mul(&sign);
            let (ds, dc) = (s1.sub(&s0), c0.sub(&c1));
            let d = vsub(&fr.o, reference);
            let moment: V3<T> = std::array::from_fn(|i| {
                d[i].mul(&sweep)
                    .add(&r.mul(&ds.mul(&fr.x[i]).add(&dc.mul(&fr.y[i]))))
                    .mul(&r)
                    .mul(&sign)
            });
            Some((length, moment))
        }
        // An ellipse's length is an elliptic integral.
        Curve3::BSpline(_) | Curve3::EllipseArc { .. } => None,
    }
}

fn solve_measure<T: Real>(view: &View, reference: [f64; 3]) -> Option<SheetMeasure> {
    let origin = v3::<T>(reference);
    let mut total = c::<T>(0.0);
    let mut moment: V3<T> = std::array::from_fn(|_| c(0.0));
    let mut add = |measure: T, m: V3<T>| {
        total = total.add(&measure);
        moment = std::array::from_fn(|i| moment[i].add(&m[i]));
    };
    if !view.faces.is_empty() {
        let loops = resolved(view);
        for (f, face) in view.faces.iter().enumerate() {
            let values = face_integrals::<T>(face, &loops[f], &origin)?;
            let m = [values[11].clone(), values[12].clone(), values[13].clone()];
            match values[10].sign()? {
                Ordering::Less => add(values[10].neg(), m.map(|x| x.neg())),
                _ => add(values[10].clone(), m),
            }
        }
    } else {
        let edges: BTreeSet<usize> = view
            .shells
            .iter()
            .flat_map(|s| s.wire_edges.iter().map(|e| e.0))
            .collect();
        for e in edges {
            let (length, m) = curve_moments::<T>(&view.edges[e].curve, &origin)?;
            add(length, m);
        }
    }
    let bounds = |x: &T| {
        let (lo, hi) = x.bounds_f64();
        (lo.is_finite() && hi.is_finite()).then_some([lo, hi])
    };
    let centre: V3<T> = if total.sign()? == Ordering::Equal {
        // An acorn: its vertex.
        origin
    } else {
        let mut out: V3<T> = std::array::from_fn(|_| c(0.0));
        for i in 0..3 {
            out[i] = moment[i].div(&total)?.add(&origin[i]);
        }
        out
    };
    Some(SheetMeasure {
        measure: bounds(&total)?,
        centre: [
            bounds(&centre[0])?,
            bounds(&centre[1])?,
            bounds(&centre[2])?,
        ],
    })
}

/// A sheet's area or a wire's length and its centre, certified.
pub(crate) fn sheet_measure(view: &View, reference: [f64; 3]) -> Option<SheetMeasure> {
    super::bernstein::clear_memo();
    solve_measure::<Fast>(view, reference).or_else(|| solve_measure::<I>(view, reference))
}
