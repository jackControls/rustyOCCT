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
/// `alpha cos(f x) + beta sin(f x)`, `f >= 0`.
fn fourier(a: u8, b: u8) -> Vec<(i64, R, R)> {
    // Frequency to (cos, sin) coefficients, negative frequencies allowed.
    let mut terms: BTreeMap<i64, (R, R)> = BTreeMap::new();
    terms.insert(0, (int(1), int(0)));
    let half = R::new(1.into(), 2.into());
    let step = |terms: &BTreeMap<i64, (R, R)>, by_sin: bool| {
        let mut out: BTreeMap<i64, (R, R)> = BTreeMap::new();
        let mut put = |f: i64, cs: R, sn: R| {
            let e = out.entry(f).or_insert((int(0), int(0)));
            e.0 += cs;
            e.1 += sn;
        };
        for (f, (cs, sn)) in terms {
            if by_sin {
                // cos(fx) sin x = [sin((f+1)x) - sin((f-1)x)] / 2
                // sin(fx) sin x = [cos((f-1)x) - cos((f+1)x)] / 2
                put(f + 1, -sn * &half, cs * &half);
                put(f - 1, sn * &half, -(cs * &half));
            } else {
                // cos(fx) cos x = [cos((f+1)x) + cos((f-1)x)] / 2
                // sin(fx) cos x = [sin((f+1)x) + sin((f-1)x)] / 2
                put(f + 1, cs * &half, sn * &half);
                put(f - 1, cs * &half, sn * &half);
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
    let mut folded: BTreeMap<i64, (R, R)> = BTreeMap::new();
    for (f, (cs, sn)) in terms {
        let (g, sn) = if f < 0 { (-f, -sn) } else { (f, sn) };
        let e = folded.entry(g).or_insert((int(0), int(0)));
        e.0 += cs;
        e.1 += if g == 0 { int(0) } else { sn };
    }
    folded
        .into_iter()
        .filter(|(_, (cs, sn))| *cs != int(0) || *sn != int(0))
        .map(|(f, (cs, sn))| (f, cs, sn))
        .collect()
}

fn binomial(n: u32, k: u32) -> R {
    let mut out = int(1);
    for i in 0..k {
        out = out * int(i64::from(n - i)) / int(i64::from(i + 1));
    }
    out
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

/// `integral over [t0, t0 + sweep] of cos^p t sin^q t dt`, exactly.
fn trig_integral<T: Real>(cos_power: u8, sin_power: u8, t0: &T, t1: &T) -> Option<T> {
    let mut total = c::<T>(0.0);
    for (f, alpha, beta) in fourier(cos_power, sin_power) {
        if f == 0 {
            total = total.add(&q::<T>(&alpha).mul(&t1.sub(t0)));
            continue;
        }
        let g = c::<T>(f as f64);
        let (c0, s0) = T::cos_sin(&t0.mul(&g));
        let (c1, s1) = T::cos_sin(&t1.mul(&g));
        let cos_part = q::<T>(&alpha).mul(&s1.sub(&s0));
        let sin_part = q::<T>(&beta).mul(&c0.sub(&c1));
        total = total.add(&cos_part.add(&sin_part).div(&g)?);
    }
    Some(total)
}

/// `-integral of F du` along a plane arc: `u = cx + r cos t`,
/// `v = cy + r sin t`, so `-F du = F r sin t dt`.
fn planar_arc<T: Real>(
    f: &Planar<T>,
    center: [f64; 2],
    radius: f64,
    start: f64,
    sweep: f64,
) -> Option<T> {
    let rad = c::<T>(radius);
    // Polynomials in (cos t, sin t).
    let u: Planar<T> = [((0, 0), c(center[0])), ((1, 0), rad.clone())]
        .into_iter()
        .collect();
    let v: Planar<T> = [((0, 0), c(center[1])), ((0, 1), rad.clone())]
        .into_iter()
        .collect();
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

/// `-integral of F du` along the line from `a` to `b` on a surface of
/// revolution; `None` when `du` may be zero without being zero.
fn rev_line<T: Real>(f: &Rev<T>, a: &V2<T>, b: &V2<T>) -> Option<T> {
    let d = b[0].sub(&a[0]);
    if d.sign()? == Ordering::Equal {
        return Some(c(0.0));
    }
    let m = b[1].sub(&a[1]).div(&d)?;
    let mut total = c::<T>(0.0);
    let mut cache: BTreeMap<(i64, u32), (T, T)> = BTreeMap::new();
    // J(f, j) = (integral over [0, d] of w^j cos(f(u0 + w)), ... sin ...).
    let mut moments = |f: i64, j: u32| -> Option<(T, T)> {
        if let Some(x) = cache.get(&(f, j)) {
            return Some(x.clone());
        }
        let out = if f == 0 {
            let power = (0..=j).fold(c::<T>(1.0), |p, _| p.mul(&d));
            (power.div(&c(f64::from(j) + 1.0))?, c(0.0))
        } else {
            let g = c::<T>(f as f64);
            let (c0, s0) = T::cos_sin(&a[0].mul(&g));
            let (c1, s1) = T::cos_sin(&b[0].mul(&g));
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
    for ((k, p, qq), x) in f {
        // v^k = (v0 + m w)^k = sum over j of C(k, j) v0^(k-j) m^j w^j.
        for j in 0..=u32::from(*k) {
            let coef = q::<T>(&binomial(u32::from(*k), j))
                .mul(&(0..u32::from(*k) - j).fold(c::<T>(1.0), |acc, _| acc.mul(&a[1])))
                .mul(&(0..j).fold(c::<T>(1.0), |acc, _| acc.mul(&m)));
            for (freq, alpha, beta) in fourier(*p, *qq) {
                let (jc, js) = moments(freq, j)?;
                let term = q::<T>(&alpha).mul(&jc).add(&q::<T>(&beta).mul(&js));
                total = total.add(&x.mul(&coef).mul(&term));
            }
        }
    }
    Some(total.neg())
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
            let d = vsub(&fr.o, reference);
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
                            planar_line(&anti[k], &[c(start.x), c(start.y)], &[c(end.x), c(end.y)])
                        }
                        Curve2::CircularArc {
                            center,
                            radius,
                            start_angle,
                            sweep_angle,
                        } => planar_arc(
                            &anti[k],
                            [center.x, center.y],
                            *radius,
                            *start_angle,
                            *sweep_angle,
                        ),
                    });
                    accumulate(values.try_map_all()?);
                }
                for (a, b) in chords::<T>(lp) {
                    let values: [Option<T>; TERMS] =
                        std::array::from_fn(|k| planar_line(&anti[k], &a, &b));
                    accumulate(values.try_map_all()?);
                }
            }
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
                    let Curve2::LineSegment { start, end } = &u.pcurve else {
                        return None;
                    };
                    let values: [Option<T>; TERMS] = std::array::from_fn(|k| {
                        rev_line(&anti[k], &[c(start.x), c(start.y)], &[c(end.x), c(end.y)])
                    });
                    accumulate(values.try_map_all()?);
                }
                for (a, b) in chords::<T>(lp) {
                    let values: [Option<T>; TERMS] =
                        std::array::from_fn(|k| rev_line(&anti[k], &a, &b));
                    accumulate(values.try_map_all()?);
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
