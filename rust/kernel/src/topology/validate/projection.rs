//! Certified evaluation of the conic edges and the `Projection` pcurves of
//! S8d.2 (D13), in either tier, pointwise and as Taylor jets.
//!
//! A `Projection` is its surface's inverse applied to its edge: the frame
//! coordinates of the edge's point, then `u = atan2(y, x)` (turned to the
//! lift: `ref + atan2` of the vector rotated back by the recorded lift, so no
//! branch cut is ever near) and `v` as the surface needs (`z`, `z / cos a`,
//! `atan2(z, rho)`, `atan2(z, rho - R)`). Integrals along it use `jet::
//! integrate` over the fraction.
use super::{c, V2, V3};
use crate::certified::{Fast, Real};
use crate::jet::{integrate_many, Jet};
use crate::topology::{Curve3, Meet, Projection, Rise, Spiric, Surface, Toric};
use crate::Frame3;

/// Integration widths and depths for the integrals along projections.
pub(super) const ORDER: usize = 12;
pub(super) const WIDTH: f64 = 1e-12;
pub(super) const DEPTH: usize = 40;
/// The width of a sign decision's integrals when tried loosely first
/// (`tiered_integral`): the validator's areas and fluxes are signs, not
/// measures, and a projection's pieces near a turning point of its meeting
/// halve many more times for `WIDTH` than for this.
pub(super) const LOOSE: f64 = 1e-6;

thread_local! {
    /// The absolute width of the sign integrals along projections: `WIDTH`,
    /// or `LOOSE` within `with_sign_width`.
    static SIGN_WIDTH: std::cell::Cell<f64> = const { std::cell::Cell::new(WIDTH) };
}

/// `f` with the sign integrals' width `width`, restored after it.
pub(super) fn with_sign_width<X>(width: f64, f: impl FnOnce() -> X) -> X {
    struct Restore(f64);
    impl Drop for Restore {
        fn drop(&mut self) {
            SIGN_WIDTH.with(|w| w.set(self.0));
        }
    }
    let _restore = Restore(SIGN_WIDTH.with(|w| w.replace(width)));
    f()
}

/// Whether the pcurve projects the fin's own edge onto its face's surface,
/// in the fin's direction: its deviation is zero by definition.
pub(super) fn own(p: &Projection, curve: &Curve3, surface: &Surface, forward: bool) -> bool {
    p.curve == *curve && p.surface == *surface && p.reversed != forward
}

fn frame_of(s: &Surface) -> Option<&Frame3> {
    match s {
        Surface::Plane(f)
        | Surface::Cylinder { frame: f, .. }
        | Surface::Cone { frame: f, .. }
        | Surface::Sphere { frame: f, .. }
        | Surface::Torus { frame: f, .. } => Some(f),
        Surface::BSpline(_) => None,
    }
}

/// A point of a conic edge (hyperbola or parabola) at a fraction, enclosed.
pub(super) fn conic_point<T: Real>(curve: &Curve3, t: f64) -> Option<V3<T>> {
    let jet = curve_jet(curve, &Jet::variable(c::<T>(t), 0))?;
    Some(jet.map(|j| j.c[0].clone()))
}

/// A world coordinate's jet from a frame's origin and axes and the local
/// coordinates' jets.
fn world<T: Real>(f: &Frame3, x: &Jet<T>, y: &Jet<T>) -> [Jet<T>; 3] {
    let (o, ex, ey) = (f.origin().to_array(), f.x().to_array(), f.y().to_array());
    std::array::from_fn(|i| {
        x.scale(&c(ex[i]))
            .add(&y.scale(&c(ey[i])))
            .add_constant(&c(o[i]))
    })
}

/// The jets of a hyperbola's or parabola's world point in the fraction.
fn conic_jet<T: Real>(curve: &Curve3, fraction: &Jet<T>) -> Option<[Jet<T>; 3]> {
    Some(match curve {
        Curve3::HyperbolaArc {
            frame,
            major,
            minor,
            start,
            sweep,
        } => {
            let t = fraction.scale(&c(*sweep)).add_constant(&c(*start));
            let (ch, sh) = t.cosh_sinh();
            world(frame, &ch.scale(&c(*major)), &sh.scale(&c(*minor)))
        }
        Curve3::ParabolaArc {
            frame,
            focal,
            start,
            sweep,
        } => {
            let t = fraction.scale(&c(*sweep)).add_constant(&c(*start));
            let x = t.square().scale(&c::<T>(0.25).div(&c(*focal))?);
            world(frame, &x, &t)
        }
        _ => return None,
    })
}

/// The jets of any line, circle, arc, ellipse, hyperbola or parabola edge's
/// world point in the fraction.
pub(super) fn curve_jet<T: Real>(curve: &Curve3, fraction: &Jet<T>) -> Option<[Jet<T>; 3]> {
    Some(match curve {
        Curve3::LineSegment { start, end } => {
            let (s, e) = (start.to_array(), end.to_array());
            std::array::from_fn(|i| {
                fraction
                    .scale(&c::<T>(e[i]).sub(&c(s[i])))
                    .add_constant(&c(s[i]))
            })
        }
        Curve3::Circle { frame, radius } => {
            let a = fraction.scale(&c(std::f64::consts::TAU));
            let (co, si) = a.cos_sin();
            world(frame, &co.scale(&c(*radius)), &si.scale(&c(*radius)))
        }
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let a = fraction
                .scale(&c(*sweep_angle))
                .add_constant(&c(*start_angle));
            let (co, si) = a.cos_sin();
            world(frame, &co.scale(&c(*radius)), &si.scale(&c(*radius)))
        }
        Curve3::EllipseArc {
            frame,
            major,
            minor,
            start_angle,
            sweep_angle,
        } => {
            let a = fraction
                .scale(&c(*sweep_angle))
                .add_constant(&c(*start_angle));
            let (co, si) = a.cos_sin();
            world(frame, &co.scale(&c(*major)), &si.scale(&c(*minor)))
        }
        Curve3::HyperbolaArc { .. } | Curve3::ParabolaArc { .. } => conic_jet(curve, fraction)?,
        Curve3::Section(s) => section_jet(s, fraction)?,
        Curve3::Meet(m) => meet_jet(m, fraction)?.1,
        Curve3::Rise(m) => rise_jet(m, fraction)?,
        Curve3::Toric(m) => toric_jet(m, fraction)?.1,
        // S9f.2b: on the knot spans its base meets.
        Curve3::WallMeet(m) => super::wall_meet::jet(m, fraction, None)?.1,
        Curve3::BSpline(_) => return None,
    })
}

/// A jet vector's dot product with a constant axis.
fn dot_axis<T: Real>(v: &[Jet<T>; 3], axis: [f64; 3]) -> Jet<T> {
    v[0].scale(&c(axis[0]))
        .add(&v[1].scale(&c(axis[1])))
        .add(&v[2].scale(&c(axis[2])))
}

/// A torus meeting's point along its parameter's jet `t` as `P0 + C P1 + S
/// P2` (`C` and `S` the other angle's cosine and sine), and the quadric's
/// affine functionals there, each `L0 + C L1 + S L2` with its sign in `G =
/// sum sign L^2`; for another torus (S9d.4b.2b) its local coordinates' and
/// its constants `R2^2 - r2^2` and `4 R2^2`: `G = (L1^2 + L2^2 + L3^2 + k)^2
/// - 4 R2^2 (L1^2 + L2^2)`.
type ToricParts<T> = ([[Jet<T>; 3]; 3], Vec<([Jet<T>; 3], bool)>, Option<(T, T)>);

fn toric_parts<T: Real>(m: &Toric, t: &Jet<T>) -> Option<ToricParts<T>> {
    let n = t.order();
    let (big, small) = (c::<T>(m.major), c::<T>(m.minor));
    let (ct, st) = t.cos_sin();
    let (o, x, y, nn) = (
        m.frame.origin().to_array(),
        m.frame.x().to_array(),
        m.frame.y().to_array(),
        m.frame.normal().to_array(),
    );
    let konst = |v: f64| Jet::constant(c::<T>(v), n);
    let p: [[Jet<T>; 3]; 3] = if m.over_v {
        // t = v: P = o + r sv n + cu (R + r cv) x + su (R + r cv) y.
        let rho = ct.scale(&small).add_constant(&big);
        let z = st.scale(&small);
        [
            std::array::from_fn(|k| z.scale(&c(nn[k])).add_constant(&c(o[k]))),
            std::array::from_fn(|k| rho.scale(&c(x[k]))),
            std::array::from_fn(|k| rho.scale(&c(y[k]))),
        ]
    } else {
        // t = u: P = o + R e + cv r e + sv r n, e = cu x + su y.
        let e: [Jet<T>; 3] = std::array::from_fn(|k| ct.scale(&c(x[k])).add(&st.scale(&c(y[k]))));
        [
            std::array::from_fn(|k| e[k].scale(&big).add_constant(&c(o[k]))),
            std::array::from_fn(|k| e[k].scale(&small)),
            std::array::from_fn(|k| konst(nn[k]).scale(&small)),
        ]
    };
    // The origin's offset first (its cancellation exact for the constant
    // part): `f . (P0 - o2)`.
    let o2 = m.other.origin().to_array();
    let rel: [Jet<T>; 3] = std::array::from_fn(|k| p[0][k].add_constant(&c::<T>(-o2[k])));
    let lin = |f: [f64; 3]| -> [Jet<T>; 3] {
        [dot_axis(&rel, f), dot_axis(&p[1], f), dot_axis(&p[2], f)]
    };
    let zero = || Jet::constant(c::<T>(0.0), n);
    let mut out = Vec::new();
    if m.other_minor > 0.0 {
        for f in [m.other.x(), m.other.y(), m.other.normal()] {
            out.push((lin(f.to_array()), true));
        }
        let (big, small) = (c::<T>(m.other_radius), c::<T>(m.other_minor));
        let big2 = big.square();
        let torus = (big2.sub(&small.square()), big2.mul(&T::exact_f64(4.0)));
        return Some((p, out, Some(torus)));
    }
    if m.other_sphere {
        for f in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            out.push((lin(f), true));
        }
        out.push(([konst(m.other_radius), zero(), zero()], false));
    } else {
        out.push((lin(m.other.x().to_array()), true));
        out.push((lin(m.other.y().to_array()), true));
        let r = c::<T>(m.other_radius);
        let rad = if m.other_half_angle == 0.0 {
            [Jet::constant(r, n), zero(), zero()]
        } else {
            let (ca, sa) = T::cos_sin(&c(m.other_half_angle));
            let tan = sa.div(&ca)?;
            let l = lin(m.other.normal().to_array());
            [
                l[0].scale(&tan).add_constant(&r),
                l[1].scale(&tan),
                l[2].scale(&tan),
            ]
        };
        out.push((rad, false));
    }
    Some((p, out, None))
}

/// Another torus's `G = S^2 - 4 R2^2 P` and its derivative from the
/// functionals' values `v` and derivatives `d` (`S = |v|^2 + k`, `P = v1^2
/// + v2^2`).
fn torus_g<T: Real>(v: &[T], d: &[T], (k, four): &(T, T)) -> (T, T) {
    let two = T::exact_f64(2.0);
    let p = v[0].square().add(&v[1].square());
    let dp = v[0].mul(&d[0]).add(&v[1].mul(&d[1])).mul(&two);
    let s = p.add(&v[2].square()).add(k);
    let ds = dp.add(&v[2].mul(&d[2]).mul(&two));
    (
        s.square().sub(&four.mul(&p)),
        s.mul(&ds).mul(&two).sub(&four.mul(&dp)),
    )
}

/// `G` and its derivative in `s` at scalars: the functionals' constant
/// terms and the other angle's enclosed cosine and sine.
fn toric_g<T: Real>(fs: &[([Jet<T>; 3], bool)], torus: Option<&(T, T)>, cs: &(T, T)) -> (T, T) {
    let (co, si) = cs;
    if let Some(t) = torus {
        let v: Vec<T> = fs
            .iter()
            .map(|(l, _)| l[0].c[0].add(&l[1].c[0].mul(co)).add(&l[2].c[0].mul(si)))
            .collect();
        let d: Vec<T> = fs
            .iter()
            .map(|(l, _)| l[2].c[0].mul(co).sub(&l[1].c[0].mul(si)))
            .collect();
        return torus_g(&v, &d, t);
    }
    let mut g = T::exact_f64(0.0);
    let mut d = T::exact_f64(0.0);
    for (l, plus) in fs {
        let v = l[0].c[0].add(&l[1].c[0].mul(co)).add(&l[2].c[0].mul(si));
        let dv = l[2].c[0].mul(co).sub(&l[1].c[0].mul(si));
        let (gv, dd) = (v.square(), v.mul(&dv).mul(&T::exact_f64(2.0)));
        if *plus {
            g = g.add(&gv);
            d = d.add(&dd);
        } else {
            g = g.sub(&gv);
            d = d.sub(&dd);
        }
    }
    (g, d)
}

/// `G_ss` and `G_sf` over a box: the functionals' values and first
/// derivatives in `f` (`fs`, jets of order one or more over the base) at
/// the other angle's enclosed cosine and sine, each functional `v = L0 + C
/// L1 + S L2` with `v_s = C L2 - S L1`, `v_ss = -(C L1 + S L2)`, `v_f` and
/// `v_sf` from the jets' first coefficients.
fn toric_second<T: Real>(
    fs: &[([Jet<T>; 3], bool)],
    torus: Option<&(T, T)>,
    cs: &(T, T),
) -> (T, T) {
    let (co, si) = cs;
    let two = T::exact_f64(2.0);
    let parts: Vec<[T; 5]> = fs
        .iter()
        .map(|(l, _)| {
            let v = l[0].c[0].add(&l[1].c[0].mul(co)).add(&l[2].c[0].mul(si));
            let vs = l[2].c[0].mul(co).sub(&l[1].c[0].mul(si));
            let vss = l[1].c[0].mul(co).add(&l[2].c[0].mul(si)).neg();
            let vf = l[0].c[1].add(&l[1].c[1].mul(co)).add(&l[2].c[1].mul(si));
            let vsf = l[2].c[1].mul(co).sub(&l[1].c[1].mul(si));
            [v, vs, vss, vf, vsf]
        })
        .collect();
    // Per functional: (v^2)_s, (v^2)_ss, (v^2)_f, (v^2)_sf, halved.
    let square = |[v, vs, vss, vf, vsf]: &[T; 5]| {
        [
            v.mul(vs),
            vs.square().add(&v.mul(vss)),
            v.mul(vf),
            vf.mul(vs).add(&v.mul(vsf)),
        ]
    };
    if let Some((k, four)) = torus {
        let [a, b, e] = [square(&parts[0]), square(&parts[1]), square(&parts[2])];
        let sum = |i: usize, with_normal: bool| {
            let p = a[i].add(&b[i]).mul(&two);
            if with_normal {
                p.add(&e[i].mul(&two))
            } else {
                p
            }
        };
        let (pss, psf) = (sum(1, false), sum(3, false));
        let (qs, qss, qf, qsf) = (sum(0, true), sum(1, true), sum(2, true), sum(3, true));
        let q = parts[0][0]
            .square()
            .add(&parts[1][0].square())
            .add(&parts[2][0].square())
            .add(k);
        let gss = qs.square().add(&q.mul(&qss)).mul(&two).sub(&four.mul(&pss));
        let gsf = qf.mul(&qs).add(&q.mul(&qsf)).mul(&two).sub(&four.mul(&psf));
        return (gss, gsf);
    }
    let (mut gss, mut gsf) = (T::exact_f64(0.0), T::exact_f64(0.0));
    for (part, (_, plus)) in parts.iter().zip(fs) {
        let [_, ss, _, sf] = square(part);
        let (ss, sf) = (ss.mul(&two), sf.mul(&two));
        if *plus {
            gss = gss.add(&ss);
            gsf = gsf.add(&sf);
        } else {
            gss = gss.sub(&ss);
            gsf = gsf.sub(&sf);
        }
    }
    (gss, gsf)
}

/// `G_f` and `G_ff` at the functionals' jets (`G_ff` from their second
/// coefficients when they have them, else zero is not claimed: `None`).
fn toric_ff<T: Real>(
    fs: &[([Jet<T>; 3], bool)],
    torus: Option<&(T, T)>,
    cs: &(T, T),
) -> (T, Option<T>) {
    let (co, si) = cs;
    let two = T::exact_f64(2.0);
    let second = fs.iter().all(|(l, _)| l.iter().all(|j| j.order() >= 2));
    let at = |l: &[Jet<T>; 3], k: usize| l[0].c[k].add(&l[1].c[k].mul(co)).add(&l[2].c[k].mul(si));
    // Per functional: v, v_f, v_ff (twice the second coefficient).
    let parts: Vec<[T; 3]> = fs
        .iter()
        .map(|(l, _)| {
            let vff = if second {
                at(l, 2).mul(&two)
            } else {
                T::exact_f64(0.0)
            };
            [at(l, 0), at(l, 1), vff]
        })
        .collect();
    // (v^2)_f and (v^2)_ff, halved.
    let square = |[v, vf, vff]: &[T; 3]| [v.mul(vf), vf.square().add(&v.mul(vff))];
    if let Some((k, four)) = torus {
        let [a, b, e] = [square(&parts[0]), square(&parts[1]), square(&parts[2])];
        let (pf, pff) = (a[0].add(&b[0]).mul(&two), a[1].add(&b[1]).mul(&two));
        let (qf, qff) = (pf.add(&e[0].mul(&two)), pff.add(&e[1].mul(&two)));
        let q = parts[0][0]
            .square()
            .add(&parts[1][0].square())
            .add(&parts[2][0].square())
            .add(k);
        let gf = q.mul(&qf).mul(&two).sub(&four.mul(&pf));
        let gff = qf.square().add(&q.mul(&qff)).mul(&two).sub(&four.mul(&pff));
        return (gf, second.then_some(gff));
    }
    let (mut gf, mut gff) = (T::exact_f64(0.0), T::exact_f64(0.0));
    for (part, (_, plus)) in parts.iter().zip(fs) {
        let [f, ff] = square(part);
        let (f, ff) = (f.mul(&two), ff.mul(&two));
        if *plus {
            gf = gf.add(&f);
            gff = gff.add(&ff);
        } else {
            gf = gf.sub(&f);
            gff = gff.sub(&ff);
        }
    }
    (gf, second.then_some(gff))
}

/// The narrower of two enclosures of one value.
fn narrower<T: Real>(a: T, b: T) -> T {
    let width = |x: &T| {
        let (lo, hi) = x.bounds_f64();
        hi - lo
    };
    if width(&b) < width(&a) {
        b
    } else {
        a
    }
}

/// A torus's meeting's angles and world point, as jets.
type ToricJet<T> = ([Jet<T>; 2], [Jet<T>; 3]);

/// The jets of a torus's meeting with a quadric (S9d.4b.2) in the fraction:
/// the other angle `s` enclosed over the base by interval Newton inside the
/// window (a unique root of `G` there), then its coefficients by the implicit
/// function theorem, term by term: the `k`-th coefficient of `G` is linear in
/// `s_k` with the slope `G_s` at the base, so `s_k` is minus the rest of that
/// coefficient (with `s_k` zero, `cos s`, `sin s` and the functionals
/// continued by their recurrences) over the slope, each an enclosure at every
/// point of the base (inclusion isotone).
fn toric_jet<T: Real>(m: &Toric, fraction: &Jet<T>) -> Option<ToricJet<T>> {
    let n = fraction.order();
    let t = fraction.scale(&c(m.sweep)).add_constant(&c(m.start));
    let (p, fs, torus) = toric_parts(m, &t)?;
    let torus = torus.as_ref();
    let (lo, hi) = t.c[0].bounds_f64();
    let star = m.root_at(0.5 * lo + 0.5 * hi);
    if !star.is_finite() {
        return None;
    }
    let s0 = c::<T>(star);
    // Interval Newton in its mean-value form over the base: `G` at the
    // base's middle and at `s*`, plus `G_f` over the base times its
    // half-width, over `G_s` (the natural extension of `G` over a wide
    // base cancels badly).
    let (flo, fhi) = fraction.c[0].bounds_f64();
    let fmid = 0.5 * flo + 0.5 * fhi;
    let point = flo == fhi;
    let at_mid = |f: f64| {
        Jet::variable(c::<T>(f), 0)
            .scale(&c(m.sweep))
            .add_constant(&c(m.start))
    };
    let (_, fs_mid, _) = toric_parts(m, &at_mid(fmid))?;
    // The functionals over the base to order two (their first and second
    // derivatives for the mean-value forms), and to order one at its middle.
    let (_, fs_one, _) = if point {
        (p.clone(), fs.clone(), None)
    } else {
        toric_parts(
            m,
            &Jet::variable(fraction.c[0].clone(), 2)
                .scale(&c(m.sweep))
                .add_constant(&c(m.start)),
        )?
    };
    let spread = fraction.c[0].sub(&c(fmid));
    let (g0, d0) = toric_g(&fs_mid, torus, &T::cos_sin(&s0));
    let size = |x: &T| {
        let (a, b) = x.bounds_f64();
        a.abs().max(b.abs())
    };
    let least = |x: &T| {
        let (a, b) = x.bounds_f64();
        if a > 0.0 {
            a
        } else if b < 0.0 {
            -b
        } else {
            0.0
        }
    };
    let slope0 = least(&d0);
    if slope0 <= 0.0 || slope0.is_nan() {
        return None;
    }
    // G's change over the base at s*: G_f times the spread.
    let drift = |cs: &(T, T)| -> T {
        let (co, si) = cs;
        if let Some(t) = torus {
            // G_f of another torus: its values and their derivatives in f.
            let (v, vf): (Vec<T>, Vec<T>) = fs_one
                .iter()
                .map(|(l, _)| {
                    let v = l[0].c[0].add(&l[1].c[0].mul(co)).add(&l[2].c[0].mul(si));
                    let vf = if point {
                        T::exact_f64(0.0)
                    } else {
                        l[0].c[1].add(&l[1].c[1].mul(co)).add(&l[2].c[1].mul(si))
                    };
                    (v, vf)
                })
                .unzip();
            return torus_g(&v, &vf, t).1.mul(&spread);
        }
        let mut out = T::exact_f64(0.0);
        for (l, plus) in &fs_one {
            let v = l[0].c[0].add(&l[1].c[0].mul(co)).add(&l[2].c[0].mul(si));
            let vf = if point {
                T::exact_f64(0.0)
            } else {
                l[0].c[1].add(&l[1].c[1].mul(co)).add(&l[2].c[1].mul(si))
            };
            let x = v.mul(&vf).mul(&T::exact_f64(2.0));
            out = if *plus { out.add(&x) } else { out.sub(&x) };
        }
        out.mul(&spread)
    };
    let mut delta =
        4.0 * (size(&g0) + size(&drift(&T::cos_sin(&s0)))) / slope0 + 1e-15 * (1.0 + star.abs());
    let mut root = None;
    for _ in 0..12 {
        if !delta.is_finite() || delta > 0.5 {
            return None;
        }
        let v = c::<T>(star - delta).union(&c(star + delta));
        let cs = T::cos_sin(&v);
        let (_, d) = toric_g(&fs_one, torus, &cs);
        // G_s over the box also in its mean-value form about `(fmid, s*)`
        // (the natural extension cancels badly); the narrower holds it.
        let d = if point {
            d
        } else {
            let (gss, gsf) = toric_second(&fs_one, torus, &cs);
            let spread_s = v.sub(&s0);
            narrower(d, d0.add(&gss.mul(&spread_s)).add(&gsf.mul(&spread)))
        };
        if let Some(q) = g0.add(&drift(&cs)).div(&d) {
            let next = s0.sub(&q);
            let (nlo, nhi) = next.bounds_f64();
            if nlo > star - delta && nhi < star + delta {
                root = Some(next);
                break;
            }
        }
        delta *= 4.0;
    }
    let s_base = root?;
    // Inside the window: the edge's own root.
    let (slo, shi) = s_base.bounds_f64();
    if !(slo > m.window[0] && shi < m.window[1]) {
        return None;
    }
    let (c0, si0) = T::cos_sin(&s_base);
    let (_, slope) = toric_g(&fs, torus, &(c0.clone(), si0.clone()));
    // Over a wide base, `G_s` and `G_f` along the meeting (the slope and
    // the first coefficient's rest) also in their mean-value forms about
    // `(fmid, s*)`, over the box from there to every point of the meeting
    // over the base (`s*` is binary64's root, not quite the meeting's);
    // the narrower of each holds it.
    let (slope, first) = if point {
        (slope, None)
    } else {
        let hull = s_base.union(&s0);
        let cs_hull = T::cos_sin(&hull);
        let (gss, gsf) = toric_second(&fs_one, torus, &cs_hull);
        let spread_s = s_base.sub(&s0);
        let slope = narrower(slope, d0.add(&gss.mul(&spread_s)).add(&gsf.mul(&spread)));
        let (_, fs_mid1, _) = toric_parts(
            m,
            &Jet::variable(c::<T>(fmid), 1)
                .scale(&c(m.sweep))
                .add_constant(&c(m.start)),
        )?;
        let (gf0, _) = toric_ff(&fs_mid1, torus, &T::cos_sin(&s0));
        let first = toric_ff(&fs_one, torus, &cs_hull)
            .1
            .map(|gff| gf0.add(&gff.mul(&spread)).add(&gsf.mul(&spread_s)));
        (slope, first)
    };
    let zero = || T::exact_f64(0.0);
    let cauchy = |x: &[T], y: &[T], k: usize| T::convolve(x, y, k);
    let mut s = vec![s_base];
    let (mut co, mut sn) = (vec![c0.clone()], vec![si0.clone()]);
    // Each functional's series along the meeting: `L0 + C L1 + S L2`.
    let mut ls: Vec<Vec<T>> = fs
        .iter()
        .map(|(l, _)| vec![l[0].c[0].add(&l[1].c[0].mul(&c0)).add(&l[2].c[0].mul(&si0))])
        .collect();
    // Another torus's `S = |L|^2 + k` (and `P = L1^2 + L2^2`) along it.
    let sums = |ls: &[Vec<T>], k: usize| -> (T, T) {
        let p = cauchy(&ls[0], &ls[0], k).add(&cauchy(&ls[1], &ls[1], k));
        (p.add(&cauchy(&ls[2], &ls[2], k)), p)
    };
    let mut ss = Vec::new();
    if let Some((kc, _)) = torus {
        ss.push(sums(&ls, 0).0.add(kc));
    }
    for k in 1..=n {
        let kk = T::exact_f64(k as f64);
        // cos and sin continued with s_k = 0: k S_k = sum C_i d_{k-1-i},
        // k C_k = -sum S_i d_{k-1-i}, d_j = (j + 1) s_{j+1} (i = 0 holds
        // s_k).
        let (mut a_s, mut a_c) = (zero(), zero());
        for i in 1..k {
            let d = s[k - i].mul(&T::exact_f64((k - i) as f64));
            a_s = a_s.add(&co[i].mul(&d));
            a_c = a_c.add(&sn[i].mul(&d));
        }
        co.push(a_c.neg().div(&kk)?);
        sn.push(a_s.div(&kk)?);
        let mut g = zero();
        for ((l, plus), series) in fs.iter().zip(ls.iter_mut()) {
            let x = l[0].c[k]
                .add(&cauchy(&l[1].c, &co, k))
                .add(&cauchy(&l[2].c, &sn, k));
            series.push(x);
            if torus.is_none() {
                let sq = cauchy(series, series, k);
                g = if *plus { g.add(&sq) } else { g.sub(&sq) };
            }
        }
        if let Some((_, four)) = torus {
            let (sk, pk) = sums(&ls, k);
            ss.push(sk);
            g = cauchy(&ss, &ss, k).sub(&four.mul(&pk));
        }
        // The first coefficient's rest is `G_f` along the meeting times
        // the base variable's own first coefficient (`-1` reversed).
        let g = match (&first, k) {
            (Some(mv), 1) => narrower(g, mv.mul(&fraction.c[1])),
            _ => g,
        };
        let sk = g.div(&slope)?.neg();
        co[k] = co[k].sub(&si0.mul(&sk));
        sn[k] = sn[k].add(&c0.mul(&sk));
        for ((l, _), series) in fs.iter().zip(ls.iter_mut()) {
            let fix = l[2].c[0].mul(&c0).sub(&l[1].c[0].mul(&si0)).mul(&sk);
            series[k] = series[k].add(&fix);
        }
        if torus.is_some() {
            ss[k] = sums(&ls, k).0;
        }
        s.push(sk);
    }
    let s = Jet { c: s };
    let (cj, sj) = (Jet { c: co }, Jet { c: sn });
    let point: [Jet<T>; 3] =
        std::array::from_fn(|k| p[0][k].add(&cj.mul(&p[1][k])).add(&sj.mul(&p[2][k])));
    let (u, v) = if m.over_v { (s, t) } else { (t, s) };
    Some(([u, v], point))
}

/// `acos(x)` on jets: the angle of `(x, sqrt(1 - x^2))`, in `[0, pi]`.
fn acos_jet<T: Real>(x: &Jet<T>) -> Option<Jet<T>> {
    let one = Jet::constant(c::<T>(1.0), x.order());
    let y = one.sub(&x.square()).sqrt()?;
    angle_near(&y, x, std::f64::consts::FRAC_PI_2)
}

/// The jets of a torus section's angles `(u, v)` in the fraction (S8d.3),
/// each near its own principal value at the base.
fn section_angles<T: Real>(s: &Spiric, fraction: &Jet<T>) -> Option<[Jet<T>; 2]> {
    let [a, b, cc, d] = s.plane.map(c::<T>);
    let (big, small) = (c::<T>(s.major), c::<T>(s.minor));
    let t = fraction.scale(&c(s.sweep)).add_constant(&c(s.start));
    let (lo, hi) = t.c[0].bounds_f64();
    let base = 0.5 * lo + 0.5 * hi;
    let sign = c::<T>(s.sign);
    Some(if s.over_v {
        let (cv, sv) = t.cos_sin();
        let q = sv
            .scale(&cc.mul(&small))
            .add_constant(&d)
            .neg()
            .div(&cv.scale(&small).add_constant(&big))?;
        let ab = a.square().add(&b.square()).sqrt();
        let x = q.scale(&c::<T>(1.0).div(&ab)?);
        // The angle of (a, b) turned back by its binary64 value, so the
        // branch cut stays opposite (b = 0 with a < 0 lies on it).
        let reference = s.plane[1].atan2(s.plane[0]);
        let (co, si) = T::cos_sin(&c(reference));
        let (ar, br) = (a.mul(&co).add(&b.mul(&si)), b.mul(&co).sub(&a.mul(&si)));
        let phi = T::atan2(&br, &ar)?.add(&c(reference));
        [acos_jet(&x)?.scale(&sign).add_constant(&phi), t]
    } else {
        let (cu, su) = t.cos_sin();
        let alpha = cu.scale(&a).add(&su.scale(&b));
        let x = alpha.scale(&small);
        let y = Jet::constant(small.mul(&cc), t.order());
        let w = x.square().add(&y.square()).sqrt()?;
        let ratio = alpha.scale(&big).add_constant(&d).neg().div(&w)?;
        let (ab, cf) = (
            s.plane[0] * base.cos() + s.plane[1] * base.sin(),
            s.plane[2],
        );
        let psi = angle_near(&y, &x, (s.minor * cf).atan2(s.minor * ab))?;
        [t, psi.add(&acos_jet(&ratio)?.scale(&sign))]
    })
}

/// The jets of a torus section's world point in the fraction (S8d.3).
fn section_jet<T: Real>(s: &Spiric, fraction: &Jet<T>) -> Option<[Jet<T>; 3]> {
    let (big, small) = (c::<T>(s.major), c::<T>(s.minor));
    let [u, v] = section_angles(s, fraction)?;
    let ((cu, su), (cv, sv)) = (u.cos_sin(), v.cos_sin());
    let rho = cv.scale(&small).add_constant(&big);
    let [x, y] = [rho.mul(&cu), rho.mul(&su)];
    let mut out = world(&s.frame, &x, &y);
    let n = s.frame.normal().to_array();
    let z = sv.scale(&small);
    for (k, o) in out.iter_mut().enumerate() {
        *o = o.add(&z.scale(&c(n[k])));
    }
    Some(out)
}

/// A meeting's `(u, v)` on its carrier and its world point, as jets.
type MeetJet<T> = ([Jet<T>; 2], [Jet<T>; 3]);

/// The jets of two cylinders' meeting (S9c.2) in the fraction: the
/// carrier's angle and ruling height, and the world point.
fn meet_jet<T: Real>(m: &Meet, fraction: &Jet<T>) -> Option<MeetJet<T>> {
    let u = fraction.scale(&c(m.sweep)).add_constant(&c(m.start));
    let (co, si) = u.cos_sin();
    let r = c::<T>(m.radius);
    let foot = world(&m.frame, &co.scale(&r), &si.scale(&r));
    let o2 = m.other.origin().to_array();
    let (x, y, n) = (
        m.frame.x().to_array(),
        m.frame.y().to_array(),
        m.frame.normal().to_array(),
    );
    let tan = |a: f64| -> Option<T> {
        let (ca, sa) = T::cos_sin(&c(a));
        sa.div(&ca)
    };
    // The ruling's direction: along the axis, leaning out on a cone
    // (S9d.3b).
    let t = tan(m.half_angle)?;
    let dir: [Jet<T>; 3] = std::array::from_fn(|k| {
        co.scale(&c(x[k]))
            .add(&si.scale(&c(y[k])))
            .scale(&t)
            .add_constant(&c(n[k]))
    });
    // w = foot - o2 and the direction along an axis.
    let along = |v: &[Jet<T>; 3], axis: &[f64; 3], shift: bool| {
        let mut out = Jet::constant(c::<T>(0.0), u.order());
        for k in 0..3 {
            let vk = if shift {
                v[k].add_constant(&c::<T>(o2[k]).neg())
            } else {
                v[k].clone()
            };
            out = out.add(&vk.scale(&c(axis[k])));
        }
        out
    };
    let axes: Vec<[f64; 3]> = if m.other_sphere {
        vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    } else {
        vec![m.other.x().to_array(), m.other.y().to_array()]
    };
    let ws: Vec<Jet<T>> = axes.iter().map(|a| along(&foot, a, true)).collect();
    let ns: Vec<Jet<T>> = axes.iter().map(|a| along(&dir, a, false)).collect();
    let zero = || Jet::constant(c::<T>(0.0), u.order());
    let mut a = ns.iter().fold(zero(), |acc, x| acc.add(&x.square()));
    let mut b = ws
        .iter()
        .zip(&ns)
        .fold(zero(), |acc, (w, x)| acc.add(&w.mul(x)));
    let mut cc = ws.iter().fold(zero(), |acc, w| acc.add(&w.square()));
    if m.other_sphere {
        cc = cc.add_constant(&c::<T>(m.other_radius).mul(&c(m.other_radius)).neg());
    } else {
        // The other's radius along its axis: `r2 + t2 (w . n2)` (a cone).
        let n2 = m.other.normal().to_array();
        let t2 = tan(m.other_half_angle)?;
        let r0 = along(&foot, &n2, true)
            .scale(&t2)
            .add_constant(&c(m.other_radius));
        let rd = along(&dir, &n2, false).scale(&t2);
        a = a.sub(&rd.square());
        b = b.sub(&r0.mul(&rd));
        cc = cc.sub(&r0.square());
    }
    let d = b.square().sub(&cc.mul(&a));
    // `(-b + s sqrt(d)) / a`, or `c / (-b - s sqrt(d))` where that cancels
    // less (the binary64 curve's own choice, by its midpoints).
    let sq = d.sqrt()?.scale(&c(m.sign));
    let (p, q) = (sq.sub(&b), sq.neg().sub(&b));
    let mid = |j: &Jet<T>| {
        let (lo, hi) = j.c[0].bounds_f64();
        (0.5 * lo + 0.5 * hi).abs()
    };
    let v = if mid(&p) >= mid(&q) {
        p.div(&a)?
    } else {
        cc.div(&q)?
    };
    let mut point = foot;
    for (k, p) in point.iter_mut().enumerate() {
        *p = p.add(&v.mul(&dir[k]));
    }
    Some(([u, v], point))
}

/// The jets of a cylinder's and a sphere's meeting over the height (S9d.2b)
/// in the fraction: `u = phi + sign acos(g(w) / rho)`, the world point.
fn rise_jet<T: Real>(m: &Rise, fraction: &Jet<T>) -> Option<[Jet<T>; 3]> {
    let w = fraction.scale(&c(m.sweep)).add_constant(&c(m.start));
    let ([a, b], g) = m.coefficients();
    let rho = c::<T>(a).square().add(&c::<T>(b).square()).sqrt();
    // The carrier's radius at the height (a cone's varies, S9d.3b.2).
    let (ca, sa) = T::cos_sin(&c(m.half_angle));
    let radius = w.scale(&sa.div(&ca)?).add_constant(&c(m.radius));
    let q = w
        .square()
        .scale(&c(g[2]))
        .add(&w.scale(&c(g[1])))
        .add_constant(&c(g[0]))
        .div(&radius.scale(&rho))?;
    // `phi` is the curve's binary64 constant (its definition's `atan2`).
    let phi = c::<T>(b.atan2(a));
    let u = acos_jet(&q)?.scale(&c(m.sign)).add_constant(&phi);
    let (co, si) = u.cos_sin();
    let mut out = world(&m.frame, &co.mul(&radius), &si.mul(&radius));
    let n = m.frame.normal().to_array();
    for (k, o) in out.iter_mut().enumerate() {
        *o = o.add(&w.scale(&c(n[k])));
    }
    Some(out)
}

/// `atan2(y, x)` near `reference`: the reference plus the angle of the
/// vector turned back by it, so the principal branch's cut stays opposite.
fn angle_near<T: Real>(y: &Jet<T>, x: &Jet<T>, reference: f64) -> Option<Jet<T>> {
    let (co, si) = T::cos_sin(&c(reference));
    let xr = x.scale(&co).add(&y.scale(&si));
    let yr = y.scale(&co).sub(&x.scale(&si));
    let theta0 = T::atan2(&yr.c[0], &xr.c[0])?.add(&c(reference));
    Jet::atan2(&yr, &xr, theta0)
}

/// The jets of a projection's `(u, v)` in the fraction: `fraction` is the
/// pcurve's own (the edge's is `1 - fraction` for a reversed use).
pub(super) fn projection_jet<T: Real>(p: &Projection, fraction: &Jet<T>) -> Option<[Jet<T>; 2]> {
    projection_jet_pinned(p, fraction, None)
}

/// `projection_jet`, a spline wall's meeting (S9f.2b) evaluated on the
/// wall's knot span `pin` when given (its integrals' pieces, split at the
/// knots exactly).
fn projection_jet_pinned<T: Real>(
    p: &Projection,
    fraction: &Jet<T>,
    pin: Option<usize>,
) -> Option<[Jet<T>; 2]> {
    let f = if p.reversed {
        fraction.neg().add_constant(&c(1.0))
    } else {
        fraction.clone()
    };
    // A spline wall's meeting on its own wall: its own parameters (S9f.2b).
    if let (Curve3::WallMeet(m), Surface::BSpline(s)) = (&p.curve, &p.surface) {
        if m.wall == *s {
            let ([u, v], _) = super::wall_meet::jet(m, &f, pin)?;
            return Some([u, v]);
        }
    }
    // A torus section on its own torus: its angles, lifted (S8d.3).
    if let (
        Curve3::Section(sec),
        Surface::Torus {
            frame,
            major,
            minor,
        },
    ) = (&p.curve, &p.surface)
    {
        if sec.frame == *frame && sec.major == *major && sec.minor == *minor {
            let [u, v] = section_angles(sec, &f)?;
            let (lo, hi) = fraction.c[0].bounds_f64();
            let lift = p.lift(0.5 * lo + 0.5 * hi);
            let near = |j: Jet<T>, target: f64| {
                let (a, b) = j.c[0].bounds_f64();
                let k = ((target - (0.5 * a + 0.5 * b)) / std::f64::consts::TAU).round();
                j.add_constant(&c(k * std::f64::consts::TAU))
            };
            return Some([near(u, lift.x), near(v, lift.y)]);
        }
    }
    // A torus's meeting on its own torus: its angles, lifted (S9d.4b.2).
    if let (
        Curve3::Toric(m),
        Surface::Torus {
            frame,
            major,
            minor,
        },
    ) = (&p.curve, &p.surface)
    {
        if m.frame == *frame && m.major == *major && m.minor == *minor {
            let ([u, v], _) = toric_jet(m, &f)?;
            let (lo, hi) = fraction.c[0].bounds_f64();
            let lift = p.lift(0.5 * lo + 0.5 * hi);
            let near = |j: Jet<T>, target: f64| {
                let (a, b) = j.c[0].bounds_f64();
                let k = ((target - (0.5 * a + 0.5 * b)) / std::f64::consts::TAU).round();
                j.add_constant(&c(k * std::f64::consts::TAU))
            };
            return Some([near(u, lift.x), near(v, lift.y)]);
        }
    }
    // Two cylinders' meeting on its carrier: its angle, lifted (S9c.2); on
    // a cone carrier the height over the cosine of its half angle
    // (S9d.3b).
    let carrier = match &p.surface {
        Surface::Cylinder { frame, radius } => Some((frame, *radius, 0.0)),
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => Some((frame, *radius, *half_angle)),
        _ => None,
    };
    if let (Curve3::Meet(m), Some((frame, radius, half))) = (&p.curve, carrier) {
        if m.frame == *frame && m.radius == radius && m.half_angle == half {
            let ([u, v], _) = meet_jet(m, &f)?;
            let v = if half == 0.0 {
                v
            } else {
                v.scale(&c::<T>(1.0).div(&T::cos_sin(&c(half)).0)?)
            };
            let (lo, hi) = fraction.c[0].bounds_f64();
            let lift = p.lift(0.5 * lo + 0.5 * hi);
            let (a, b) = u.c[0].bounds_f64();
            let k = ((lift.x - (0.5 * a + 0.5 * b)) / std::f64::consts::TAU).round();
            return Some([u.add_constant(&c(k * std::f64::consts::TAU)), v]);
        }
    }
    let point = match (&p.curve, pin) {
        (Curve3::WallMeet(m), Some(_)) => super::wall_meet::jet(m, &f, pin)?.1,
        _ => curve_jet(&p.curve, &f)?,
    };
    let frame = frame_of(&p.surface)?;
    let (o, x, y, n) = (
        frame.origin().to_array(),
        frame.x().to_array(),
        frame.y().to_array(),
        frame.normal().to_array(),
    );
    let rel: [Jet<T>; 3] = std::array::from_fn(|i| point[i].add_constant(&c::<T>(-o[i])));
    let dot = |axis: [f64; 3]| {
        rel[0]
            .scale(&c(axis[0]))
            .add(&rel[1].scale(&c(axis[1])))
            .add(&rel[2].scale(&c(axis[2])))
    };
    let (lx, ly, lz) = (dot(x), dot(y), dot(n));
    let (lo, hi) = fraction.c[0].bounds_f64();
    let lift = p.lift(0.5 * lo + 0.5 * hi);
    Some(match &p.surface {
        Surface::Plane(_) => [lx, ly],
        Surface::Cylinder { .. } => [angle_near(&ly, &lx, lift.x)?, lz],
        Surface::Cone { half_angle, .. } => {
            let (ca, _) = T::cos_sin(&c(*half_angle));
            [
                angle_near(&ly, &lx, lift.x)?,
                lz.div(&Jet::constant(ca, lz.order()))?,
            ]
        }
        Surface::Sphere { .. } => {
            let rho = lx.square().add(&ly.square()).sqrt()?;
            [
                angle_near(&ly, &lx, lift.x)?,
                angle_near(&lz, &rho, lift.y)?,
            ]
        }
        Surface::Torus { major, .. } => {
            let rho = lx.square().add(&ly.square()).sqrt()?;
            let r = rho.add_constant(&c::<T>(-*major));
            [angle_near(&ly, &lx, lift.x)?, angle_near(&lz, &r, lift.y)?]
        }
        Surface::BSpline(_) => return None,
    })
}

/// A projection's `(u, v)` at a fraction, enclosed.
pub(super) fn projection_at<T: Real>(p: &Projection, t: f64) -> Option<V2<T>> {
    let [u, v] = projection_jet(p, &Jet::variable(c::<T>(t), 0))?;
    Some([u.c[0].clone(), v.c[0].clone()])
}

/// An integrand along a projection: from the jets of `u`, `v`, `u'`, `v'`.
pub(super) type AlongIntegrand<'a, T> =
    &'a dyn Fn(&Jet<T>, &Jet<T>, &Jet<T>, &Jet<T>) -> Option<Jet<T>>;

/// Integrands along a projection: from the jets of `u`, `v`, `u'`, `v'`.
pub(super) type AlongIntegrands<'a, T> =
    &'a dyn Fn(&Jet<T>, &Jet<T>, &Jet<T>, &Jet<T>) -> Option<Vec<Jet<T>>>;

/// An enclosure of `integral_0^1 g(u, v, u', v') df` along a projection,
/// `g` given the jets of `u`, `v` and their derivatives in the fraction.
pub(super) fn integrate_along<T: Real>(p: &Projection, g: AlongIntegrand<'_, T>) -> Option<T> {
    along(
        p,
        1,
        &|u, v, du, dv| g(u, v, du, dv).map(|j| vec![j]),
        false,
    )?
    .pop()
}

/// The integrals of `n` integrands along a projection at once, its jets
/// computed once per piece, each to `WIDTH`, relative to its scale for mass
/// moments and absolute for sign decisions.
pub(super) fn integrate_along_many<T: Real>(
    p: &Projection,
    n: usize,
    relative: bool,
    g: AlongIntegrands<'_, T>,
) -> Option<Vec<T>> {
    along(p, n, g, relative)
}

/// A memo of projections' jets in the binary64 tier (the integrals along
/// one pcurve, for the validator's signs and the mass's moments, visit the
/// same dyadic pieces of its fraction): keyed by the projection's content
/// (its `Debug` text, which tells every binary64 number apart, `-0.0` from
/// `0.0` among them), a piece's bounds and the order. A hit returns the
/// value the evaluation would, bit for bit.
struct JetMemo {
    ids: std::collections::HashMap<String, u64>,
    next: u64,
    jets: std::collections::HashMap<PieceKey, Option<[Jet<Fast>; 2]>>,
}

/// A memo entry's key: the projection's id, a piece's bounds' bits and
/// the jets' order.
type PieceKey = (u64, u64, u64, usize);

/// Projections (their texts, up to some kilobytes each) and pieces kept
/// before the memo starts again (ids are never reused).
const MEMO_PROJECTIONS: usize = 1 << 10;
const MEMO_LIMIT: usize = 1 << 14;

thread_local! {
    static MEMO: std::cell::RefCell<JetMemo> = std::cell::RefCell::new(JetMemo {
        ids: std::collections::HashMap::new(),
        next: 0,
        jets: std::collections::HashMap::new(),
    });
}

/// The memo's id of a projection's content, `None` for content the text
/// cannot tell apart (a NaN).
pub(super) fn memo_id(p: &Projection) -> Option<u64> {
    let key = format!("{p:?}");
    if key.contains("NaN") {
        return None;
    }
    MEMO.with(|m| {
        let mut m = m.borrow_mut();
        if let Some(id) = m.ids.get(&key) {
            return Some(*id);
        }
        if m.ids.len() >= MEMO_PROJECTIONS {
            m.ids.clear();
            m.jets.clear();
        }
        let id = m.next;
        m.next += 1;
        m.ids.insert(key, id);
        Some(id)
    })
}

/// `projection_jet` of the variable about `base` to `order`, through the
/// memo in the binary64 tier (`id` from `memo_id`).
pub(super) fn projection_jet_about<T: Real>(
    p: &Projection,
    id: Option<u64>,
    base: &T,
    order: usize,
) -> Option<[Jet<T>; 2]> {
    let fresh = || projection_jet(p, &Jet::variable(base.clone(), order));
    let (Some(id), Some(fast)) = (id, base.as_fast()) else {
        return fresh();
    };
    let (lo, hi) = fast.bounds_f64();
    let key = (id, lo.to_bits(), hi.to_bits(), order);
    let to_t = |j: &Jet<Fast>| -> Option<Jet<T>> {
        Some(Jet {
            c: j.c.iter().map(|x| T::of_fast(*x)).collect::<Option<_>>()?,
        })
    };
    let hit = MEMO.with(|m| m.borrow().jets.get(&key).cloned());
    let value = match hit {
        Some(value) => value,
        None => {
            let value = projection_jet(p, &Jet::variable(fast, order));
            MEMO.with(|m| {
                let mut m = m.borrow_mut();
                if m.jets.len() >= MEMO_LIMIT {
                    m.jets.clear();
                }
                m.jets.insert(key, value.clone());
            });
            value
        }
    };
    let [u, v] = value?;
    Some([to_t(&u)?, to_t(&v)?])
}

fn along<T: Real>(
    p: &Projection,
    n: usize,
    g: AlongIntegrands<'_, T>,
    relative: bool,
) -> Option<Vec<T>> {
    // Rational jets of this order grow without bound: binary64 intervals
    // only (their failure is reported, never retried exactly).
    if T::EXACT {
        return None;
    }
    // A spline wall's meeting (S9f.2b): piece by piece between its wall's
    // knots.
    if matches!(p.curve, Curve3::WallMeet(_)) {
        let mut totals = vec![T::exact_f64(0.0); n];
        wall_pieces(
            p,
            n,
            relative,
            &|_, u, v, du, dv| g(u, v, du, dv),
            &mut totals,
        )?;
        return Some(totals);
    }
    let id = memo_id(p);
    let integrand = |f: &Jet<T>| {
        // One order more, so the derivatives keep the order asked for.
        let [u, v] = projection_jet_about(p, id, &f.c[0], f.order() + 1)?;
        let (du, dv) = (u.derivative(), v.derivative());
        let cut = |j: &Jet<T>| Jet {
            c: j.c[..=f.order()].to_vec(),
        };
        g(&cut(&u), &cut(&v), &cut(&du), &cut(&dv))
    };
    let width = if relative {
        WIDTH
    } else {
        SIGN_WIDTH.with(|w| w.get())
    };
    integrate_many(&integrand, n, 0.0, 1.0, ORDER, width, DEPTH, relative)
}

/// An integrand along a spline wall's meeting's piece on knot span `k`:
/// from `k` and the jets of `u`, `v`, `u'`, `v'` (in the pcurve's
/// fraction).
pub(super) type WallIntegrands<'a, T> =
    &'a dyn Fn(usize, &Jet<T>, &Jet<T>, &Jet<T>, &Jet<T>) -> Option<Vec<Jet<T>>>;

/// The integrals along a projection of a spline wall's meeting (S9f.2b),
/// added to `totals`: over each piece of its fraction between the wall's
/// knots, exactly (`wall_meet::pieces`), its jets on that knot span's
/// polynomial, the piece's fraction `g = ga + tau (gb - ga)` integrated
/// over `tau` in `[0, 1]` (the integrands, forms in `u'` and `v'`, given the
/// fraction's derivatives).
pub(super) fn wall_pieces<T: Real>(
    p: &Projection,
    n: usize,
    relative: bool,
    g: WallIntegrands<'_, T>,
    totals: &mut [T],
) -> Option<()> {
    let Curve3::WallMeet(m) = &p.curve else {
        return None;
    };
    // Binary64 intervals only, as `along`.
    if T::EXACT {
        return None;
    }
    let one = num_rational::BigRational::from_integer(1.into());
    for (fa, fb, k) in super::wall_meet::pieces(m)? {
        let (ga, gb) = if p.reversed {
            (&one - &fb, &one - &fa)
        } else {
            (fa, fb)
        };
        let len = &gb - &ga;
        let (start, scale) = (T::from_r(&ga), T::from_r(&len));
        let per = T::from_r(&(&one / &len));
        let integrand = |tau: &Jet<T>| {
            let t1 = Jet::variable(tau.c[0].clone(), tau.order() + 1);
            let at = t1.scale(&scale).add_constant(&start);
            let [u, v] = projection_jet_pinned(p, &at, Some(k))?;
            let (du, dv) = (u.derivative().scale(&per), v.derivative().scale(&per));
            let cut = |j: &Jet<T>| Jet {
                c: j.c[..=tau.order()].to_vec(),
            };
            Some(
                g(k, &cut(&u), &cut(&v), &cut(&du), &cut(&dv))?
                    .iter()
                    .map(|x| x.scale(&scale))
                    .collect(),
            )
        };
        let width = if relative {
            WIDTH
        } else {
            SIGN_WIDTH.with(|w| w.get())
        };
        let values = integrate_many(&integrand, n, 0.0, 1.0, ORDER, width, DEPTH, relative)?;
        for (t, x) in totals.iter_mut().zip(values) {
            *t = t.add(&x);
        }
    }
    Some(())
}

/// Crossings of the `+u` ray from `p` with a projection pcurve (half-open
/// in `v`: a piece counts when exactly one of its ends lies strictly above
/// `p.v`), from certified pieces: none where `v` or `u` keeps clear of the
/// ray, one where `v` is monotone (its derivative's enclosure excludes zero)
/// with its ends on either side and `u` right of `p` all along; others are
/// bisected, at most 40 times. `None` when undecided.
pub(super) fn crossings<T: Real>(pr: &Projection, p: &V2<T>) -> Option<u32> {
    use std::cmp::Ordering;
    if T::EXACT {
        return None;
    }
    let at = |f: f64| projection_at::<T>(pr, f);
    let mut count = 0;
    let mut stack = vec![(0.0f64, 1.0f64, 0usize)];
    while let Some((lo, hi, depth)) = stack.pop() {
        let base = T::exact_f64(lo).union(&T::exact_f64(hi));
        let mid = 0.5 * lo + 0.5 * hi;
        // Not enclosed over the whole piece: halves.
        let Some([u, v]) = projection_jet(pr, &Jet::variable(base, 1)) else {
            if depth >= 40 || !(lo < mid && mid < hi) {
                return None;
            }
            stack.push((lo, mid, depth + 1));
            stack.push((mid, hi, depth + 1));
            continue;
        };
        let dv = v.c[0].sub(&p[1]).sign();
        let du = u.c[0].sub(&p[0]).sign();
        // Clear of the ray over the whole piece.
        if matches!(dv, Some(Ordering::Less) | Some(Ordering::Greater))
            || du == Some(Ordering::Less)
        {
            continue;
        }
        let monotone = matches!(
            v.c[1].sign(),
            Some(Ordering::Less) | Some(Ordering::Greater)
        );
        if monotone {
            let (a, b) = (at(lo)?, at(hi)?);
            let above = |x: &T| Some(x.sub(&p[1]).sign()? == Ordering::Greater);
            if above(&a[1])? == above(&b[1])? {
                continue;
            }
            if du == Some(Ordering::Greater) {
                count += 1;
                continue;
            }
        }
        if depth >= 40 || !(lo < mid && mid < hi) {
            return None;
        }
        stack.push((lo, mid, depth + 1));
        stack.push((mid, hi, depth + 1));
    }
    Some(count)
}

/// Crossings of the `+v` ray from `p` with a projection pcurve over every
/// `u` alias `k TAU` (a cylinder's cover; S9c.2), each `+1` where it runs
/// in `-u` and `-1` in `+u` (half-open in `u`, as a segment's), from
/// certified pieces: a piece counts at an alias when `v` lies above `p`
/// all along it, `u` is monotone and its ends lie on either side; one
/// clear of every alias in `u`, or below `p` in `v`, counts nothing;
/// others are bisected, at most 40 times. `None` when undecided.
pub(super) fn cover_crossings<T: Real>(pr: &Projection, p: &V2<T>) -> Option<Vec<i64>> {
    use std::cmp::Ordering;
    use std::f64::consts::TAU;
    if T::EXACT {
        return None;
    }
    let mut out = Vec::new();
    let mut stack = vec![(0.0f64, 1.0f64, 0usize)];
    while let Some((lo, hi, depth)) = stack.pop() {
        let mid = 0.5 * lo + 0.5 * hi;
        let base = T::exact_f64(lo).union(&T::exact_f64(hi));
        let decided = (|| -> Option<Vec<i64>> {
            let [u, v] = projection_jet(pr, &Jet::variable(base, 1))?;
            let ((ul, uh), (pl, ph)) = (u.c[0].bounds_f64(), p[0].bounds_f64());
            let kmin = ((pl - uh) / TAU).floor() as i64 - 1;
            let kmax = ((ph - ul) / TAU).ceil() as i64 + 1;
            if kmax - kmin > 64 {
                return None;
            }
            let dv = v.c[0].sub(&p[1]).sign();
            let monotone = matches!(
                u.c[1].sign(),
                Some(Ordering::Less) | Some(Ordering::Greater)
            );
            let mut hits = Vec::new();
            for k in kmin..=kmax {
                let uk = p[0].sub(&T::exact_f64(TAU * k as f64));
                if matches!(
                    u.c[0].sub(&uk).sign(),
                    Some(Ordering::Less) | Some(Ordering::Greater)
                ) || dv == Some(Ordering::Less)
                {
                    continue;
                }
                if dv != Some(Ordering::Greater) || !monotone {
                    return None;
                }
                let (a, b) = (projection_at::<T>(pr, lo)?, projection_at::<T>(pr, hi)?);
                let above = |x: &T| Some(x.sub(&uk).sign()? == Ordering::Greater);
                let (ra, rb) = (above(&a[0])?, above(&b[0])?);
                if ra != rb {
                    hits.push(if rb { -1 } else { 1 });
                }
            }
            Some(hits)
        })();
        match decided {
            Some(hits) => out.extend(hits),
            None => {
                if depth >= 40 || !(lo < mid && mid < hi) {
                    return None;
                }
                stack.push((lo, mid, depth + 1));
                stack.push((mid, hi, depth + 1));
            }
        }
    }
    Some(out)
}

/// The certified ranges `[u_lo, u_hi, v_lo, v_hi]` a projection reaches:
/// its enclosures over `pieces` equal parts of the fraction.
pub(crate) fn projection_range(p: &Projection, pieces: usize) -> Option<[f64; 4]> {
    use crate::certified::Fast;
    let mut out = [
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ];
    for k in 0..pieces {
        let (a, b) = (k as f64 / pieces as f64, (k + 1) as f64 / pieces as f64);
        let base = Fast::exact_f64(a).union(&Fast::exact_f64(b));
        let [u, v] = projection_jet(p, &Jet::variable(base, 0))?;
        let ((ul, uh), (vl, vh)) = (u.c[0].bounds_f64(), v.c[0].bounds_f64());
        out = [
            out[0].min(ul),
            out[1].max(uh),
            out[2].min(vl),
            out[3].max(vh),
        ];
    }
    out.iter().all(|x| x.is_finite()).then_some(out)
}

/// A conic edge's or a torus section's point at a fraction in binary64
/// intervals (tessellation).
pub(crate) fn conic_point_fast(curve: &Curve3, t: f64) -> Option<[crate::certified::Fast; 3]> {
    conic_point(curve, t)
}

/// A curve's rates over its fraction (S8d.3's sections): the largest
/// `|C''|` and the largest `|C''| / |C'|` (the tangent's turn per unit
/// fraction), from interval jets on `pieces` equal parts.
pub(crate) fn section_rates(curve: &Curve3, pieces: usize) -> Option<(f64, f64)> {
    use crate::certified::Fast;
    let (mut second, mut turn) = (0.0f64, 0.0f64);
    // A spline wall's meeting's pieces split at its wall's knots too
    // (S9f.2b: its second derivative jumps there).
    let mut cuts: Vec<f64> = (0..=pieces).map(|k| k as f64 / pieces as f64).collect();
    if let Curve3::WallMeet(m) = curve {
        cuts.extend(super::wall_meet::knot_fractions(m));
        cuts.sort_by(f64::total_cmp);
        cuts.dedup();
    }
    for w in cuts.windows(2) {
        let (a, b) = (w[0], w[1]);
        let base = Fast::exact_f64(a).union(&Fast::exact_f64(b));
        let jet = curve_jet(curve, &Jet::variable(base, 2))?;
        let norm = |k: usize| {
            let parts = jet.iter().map(|j| {
                let (lo, hi) = j.c[k].bounds_f64();
                (
                    lo.abs().max(hi.abs()),
                    if lo > 0.0 {
                        lo
                    } else if hi < 0.0 {
                        -hi
                    } else {
                        0.0
                    },
                )
            });
            parts.fold((0.0f64, 0.0f64), |(big, small), (b, s)| {
                (big + b * b, small + s * s)
            })
        };
        let (d2, _) = norm(2);
        let (_, d1_low) = norm(1);
        let c2 = 2.0 * d2.sqrt();
        let c1 = d1_low.sqrt();
        if !(c2.is_finite() && c1 > 0.0) {
            return None;
        }
        second = second.max(c2);
        turn = turn.max(c2 / c1);
    }
    Some((second * (1.0 + 1e-12), turn * (1.0 + 1e-12)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::certified::Fast;
    use crate::{Point3, Tolerance, Vec3};

    fn frame(o: [f64; 3], n: [f64; 3], x: [f64; 3]) -> Frame3 {
        Frame3::new(
            Point3::new(o[0], o[1], o[2]),
            Vec3::new(n[0], n[1], n[2]),
            Vec3::new(x[0], x[1], x[2]),
            Tolerance::default(),
        )
        .unwrap()
    }

    /// A thin pipe (radius 1 about x, offset 0.5 along y) through a thick
    /// one (radius 2 about z): a ring over the thin one's angle lies on both
    /// cylinders, and its jets enclose the binary64 points and differences.
    #[test]
    fn meetings_lie_on_both_cylinders() {
        let thin = frame([0.0, 0.5, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let thick = frame([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        for sign in [1.0, -1.0] {
            let m = Meet {
                frame: thin,
                radius: 1.0,
                other: thick,
                half_angle: 0.0,
                other_radius: 2.0,
                other_sphere: false,
                other_half_angle: 0.0,
                sign,
                start: 0.0,
                sweep: std::f64::consts::TAU,
            };
            for k in 0..=16 {
                let f = k as f64 / 16.0;
                let p = m.point(f);
                let [x, y, _] = thick.coordinates(p);
                assert!((x.hypot(y) - 2.0).abs() < 1e-14, "{p:?}");
                let [_, b, c] = [p.x, p.y - 0.5, p.z];
                assert!((b.hypot(c) - 1.0).abs() < 1e-14);
                assert_eq!(p.x.signum(), sign);
                let (_, jet) = meet_jet(&m, &Jet::variable(Fast::exact_f64(f), 2)).unwrap();
                let h = 1e-6;
                let (q0, q1) = (m.point(f - h), m.point(f + h));
                for (i, j) in jet.iter().enumerate() {
                    let at = [p.x, p.y, p.z][i];
                    let (lo, hi) = j.c[0].bounds_f64();
                    assert!(lo - 1e-15 <= at && at <= hi + 1e-15);
                    let slope = ([q1.x, q1.y, q1.z][i] - [q0.x, q0.y, q0.z][i]) / (2.0 * h);
                    let (lo, hi) = j.c[1].bounds_f64();
                    assert!((slope - 0.5 * (lo + hi)).abs() < 1e-6, "{slope} {lo} {hi}");
                }
            }
        }
    }

    /// A cylinder of radius 1 about z and a sphere of radius 1.5 centred at
    /// (1, 0, 0): its loop over the height, on both surfaces, with jets
    /// enclosing the binary64 points and differences.
    #[test]
    fn rises_lie_on_both_surfaces() {
        let cyl = frame([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        for sign in [1.0, -1.0] {
            let m = Rise {
                frame: cyl,
                radius: 1.0,
                half_angle: 0.0,
                centre: Point3::new(1.0, 0.0, 0.0),
                sphere_radius: 1.5,
                sign,
                start: -0.5,
                sweep: 1.0,
            };
            for k in 0..=8 {
                let f = k as f64 / 8.0;
                let p = m.point(f);
                assert!((p.x.hypot(p.y) - 1.0).abs() < 1e-14);
                assert!(((p - Point3::new(1.0, 0.0, 0.0)).length() - 1.5).abs() < 1e-14);
                let jet = rise_jet(&m, &Jet::variable(Fast::exact_f64(f), 1)).unwrap();
                let h = 1e-6;
                let (q0, q1) = (m.point(f - h), m.point(f + h));
                for (i, j) in jet.iter().enumerate() {
                    let at = [p.x, p.y, p.z][i];
                    let (lo, hi) = j.c[0].bounds_f64();
                    assert!(lo - 1e-14 <= at && at <= hi + 1e-14);
                    let slope = ([q1.x, q1.y, q1.z][i] - [q0.x, q0.y, q0.z][i]) / (2.0 * h);
                    let (lo, hi) = j.c[1].bounds_f64();
                    assert!((slope - 0.5 * (lo + hi)).abs() < 1e-6, "{slope} {lo} {hi}");
                }
            }
        }
    }

    /// A rod of radius 0.5 along x through a torus of radii 2.5 and 1
    /// (S9d.4b.2): over `u` near 0 its meeting's top root lies near `v =
    /// pi / 6`, on both surfaces; its jets enclose the binary64 points and
    /// differences, and over an interval base every point in it.
    #[test]
    fn toric_meetings_lie_on_both_surfaces() {
        let torus = frame([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let rod = frame([-4.5, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let m = Toric {
            frame: torus,
            major: 2.5,
            minor: 1.0,
            other: rod,
            other_radius: 0.5,
            other_sphere: false,
            other_half_angle: 0.0,
            other_minor: 0.0,
            over_v: false,
            window: [0.2, 0.9],
            start: -0.05,
            sweep: 0.1,
        };
        for k in 0..=8 {
            let f = k as f64 / 8.0;
            let p = m.point(f);
            assert!((p.y.hypot(p.z) - 0.5).abs() < 1e-14, "{p:?}");
            let rho = p.x.hypot(p.y);
            assert!(((rho - 2.5).hypot(p.z) - 1.0).abs() < 1e-14, "{p:?}");
            let (_, jet) = toric_jet(&m, &Jet::variable(Fast::exact_f64(f), 2)).unwrap();
            let h = 1e-6;
            let (q0, q1) = (m.point(f - h), m.point(f + h));
            for (i, j) in jet.iter().enumerate() {
                let at = [p.x, p.y, p.z][i];
                let (lo, hi) = j.c[0].bounds_f64();
                assert!(lo - 1e-14 <= at && at <= hi + 1e-14);
                let slope = ([q1.x, q1.y, q1.z][i] - [q0.x, q0.y, q0.z][i]) / (2.0 * h);
                let (lo, hi) = j.c[1].bounds_f64();
                assert!((slope - 0.5 * (lo + hi)).abs() < 1e-6, "{slope} {lo} {hi}");
            }
        }
        let base = Fast::exact_f64(0.25).union(&Fast::exact_f64(0.5));
        let (_, jet) = toric_jet(&m, &Jet::variable(base, 3)).unwrap();
        for f in [0.25, 0.3, 0.4, 0.5] {
            let p = m.point(f).to_array();
            for (i, j) in jet.iter().enumerate() {
                let (lo, hi) = j.c[0].bounds_f64();
                assert!(lo - 1e-14 <= p[i] && p[i] <= hi + 1e-14, "{f} {i}");
            }
        }
    }

    /// A torus of radii 1.25 and 0.5 about x ringing the tube of one of
    /// radii 2.5 and 1 about z at `u = pi / 2` (S9d.4b.2b): over `v` near 0
    /// the ring's side at `u` near `pi / 2 + 0.124` lies on both tori; its
    /// jets enclose the binary64 points and differences, and over an
    /// interval base every point in it.
    #[test]
    fn toric_meetings_with_a_torus_lie_on_both() {
        let torus = frame([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let ring = frame([0.0, 2.5, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let half = std::f64::consts::FRAC_PI_2;
        let m = Toric {
            frame: torus,
            major: 2.5,
            minor: 1.0,
            other: ring,
            other_radius: 1.25,
            other_sphere: false,
            other_half_angle: 0.0,
            other_minor: 0.5,
            over_v: true,
            window: [half + 0.05, half + 0.2],
            start: -0.05,
            sweep: 0.1,
        };
        for k in 0..=8 {
            let f = k as f64 / 8.0;
            let p = m.point(f);
            let rho = p.x.hypot(p.y);
            assert!(((rho - 2.5).hypot(p.z) - 1.0).abs() < 1e-14, "{p:?}");
            let core = (p.y - 2.5).hypot(p.z);
            assert!(((core - 1.25).hypot(p.x) - 0.5).abs() < 1e-14, "{p:?}");
            let (_, jet) = toric_jet(&m, &Jet::variable(Fast::exact_f64(f), 2)).unwrap();
            let h = 1e-6;
            let (q0, q1) = (m.point(f - h), m.point(f + h));
            for (i, j) in jet.iter().enumerate() {
                let at = [p.x, p.y, p.z][i];
                let (lo, hi) = j.c[0].bounds_f64();
                assert!(lo - 1e-14 <= at && at <= hi + 1e-14);
                let slope = ([q1.x, q1.y, q1.z][i] - [q0.x, q0.y, q0.z][i]) / (2.0 * h);
                let (lo, hi) = j.c[1].bounds_f64();
                assert!((slope - 0.5 * (lo + hi)).abs() < 1e-6, "{slope} {lo} {hi}");
            }
        }
        let base = Fast::exact_f64(0.25).union(&Fast::exact_f64(0.5));
        let (_, jet) = toric_jet(&m, &Jet::variable(base, 3)).unwrap();
        for f in [0.25, 0.3, 0.4, 0.5] {
            let p = m.point(f).to_array();
            for (i, j) in jet.iter().enumerate() {
                let (lo, hi) = j.c[0].bounds_f64();
                assert!(lo - 1e-14 <= p[i] && p[i] <= hi + 1e-14, "{f} {i}");
            }
        }
    }
}
