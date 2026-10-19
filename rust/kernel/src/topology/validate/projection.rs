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
use crate::cr::Rn as _;
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

/// The narrow tier's width relative to an integrand's scale (`2^-30`) and
/// the floor of that scale (`2^-80`, `with_narrow`): its integrals, signs
/// and mass terms alike, are enclosed relative to their integrands' own
/// sizes, where a shell thinner than its enclosures at the tight tier's
/// absolute widths (a shallow crossing's lens, HANDOFF item 5 (d)) is
/// decided; the floor keeps an integrand vanishing but for rounding from
/// halving its pieces to the depth limit.
pub(super) const NARROW_WIDTH: f64 = 9.313_225_746_154_785e-10;
pub(super) const NARROW_FLOOR: f64 = 8.271_806_125_530_277e-25;

thread_local! {
    /// The absolute width of the sign integrals along projections: `WIDTH`,
    /// or `LOOSE` within `with_sign_width`.
    static SIGN_WIDTH: std::cell::Cell<f64> = const { std::cell::Cell::new(WIDTH) };
    /// Whether the narrow tier runs (`with_narrow`).
    static NARROW: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// `f` in the narrow tier, restored after it: the integrals along
/// projections relative to their integrands' scales down to
/// `NARROW_FLOOR`, and the faces' integrands along them in centred forms
/// about a loop's own height (`narrow`'s callers).
pub(super) fn with_narrow<X>(f: impl FnOnce() -> X) -> X {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            NARROW.with(|w| w.set(self.0));
        }
    }
    let _restore = Restore(NARROW.with(|w| w.replace(true)));
    f()
}

/// Whether the narrow tier runs.
pub(super) fn narrow() -> bool {
    NARROW.with(|w| w.get())
}

/// Two-variable integrands' jets: `(u, v)` to each integrand's.
pub(super) type Pair<'a, T> = &'a dyn Fn(&Jet<T>, &Jet<T>) -> Option<Vec<Jet<T>>>;

/// `f` at the jets `(u, v)` in centred form (the narrow tier): `f` at point
/// jets inside them plus its partials `fu` and `fv` over them times the
/// deviations from the points, as the mean value theorem gives it in the
/// algebra of truncated series (each partial over the segment between the
/// point and any jet inside the enclosures lies in the partial over the
/// enclosures). Where a pcurve's point is enclosed widely along its face
/// (a shallow meeting's height from a small discriminant) and the
/// integrand is a sum of large terms that cancel, the natural form's
/// widths are the terms' times the point's; the centred form's are the
/// integrand's rate of change times it. Natural where no point splits
/// off (the rational tier); a partial is not taken against a deviation
/// that is exactly zero (a carrier's own angle).
pub(super) fn centred<T: Real>(
    u: &Jet<T>,
    v: &Jet<T>,
    f: Pair<'_, T>,
    fu: Pair<'_, T>,
    fv: Pair<'_, T>,
) -> Option<Vec<Jet<T>>> {
    let (Some((um, du)), Some((vm, dv))) = (u.centre(), v.centre()) else {
        return f(u, v);
    };
    let mut out = f(&um, &vm)?;
    for (deviation, partial) in [(&du, fu), (&dv, fv)] {
        if deviation.is_zero() {
            continue;
        }
        let rates = partial(u, v)?;
        if rates.len() != out.len() {
            return None;
        }
        for (o, r) in out.iter_mut().zip(&rates) {
            *o = o.add(&r.mul(deviation));
        }
    }
    Some(out)
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

/// A torus meeting's point along its parameter's jet `t` as `P0 + C P1 + S
/// P2` (`C` and `S` the other angle's cosine and sine).
fn toric_point_parts<T: Real>(m: &Toric, t: &Jet<T>) -> [[Jet<T>; 3]; 3] {
    let n = t.order();
    let (big, small) = (c::<T>(m.major), c::<T>(m.minor));
    let (ct, st) = t.cos_sin();
    let (o, x, y, nn) = (
        m.frame.origin().to_array(),
        m.frame.x().to_array(),
        m.frame.y().to_array(),
        m.frame.normal().to_array(),
    );
    if m.over_v {
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
            std::array::from_fn(|k| Jet::constant(c::<T>(nn[k]), n).scale(&small)),
        ]
    }
}

/// An angle's harmonics `[1, cos a, sin a, cos 2a, sin 2a]`, enclosed.
type Harmonics<T> = [T; 5];

fn harmonics<T: Real>(a: &T) -> Harmonics<T> {
    let (co, si) = T::cos_sin(a);
    let (co2, si2) = T::cos_sin(&a.mul(&T::exact_f64(2.0)));
    [T::exact_f64(1.0), co, si, co2, si2]
}

/// The harmonics' derivatives in the angle: `[0, -sin a, cos a, -2 sin 2a,
/// 2 cos 2a]`.
fn harmonics_rate<T: Real>(h: &Harmonics<T>) -> Harmonics<T> {
    let two = T::exact_f64(2.0);
    [
        T::exact_f64(0.0),
        h[2].neg(),
        h[1].clone(),
        h[4].mul(&two).neg(),
        h[3].mul(&two),
    ]
}

/// A torus meeting's `G` (S9d.4b.2) as trigonometric polynomials of degree
/// two in each of the torus's angles, `sum w[a][b] D_a(s) D_b(t)` over the
/// harmonics `D` of the other angle `s` (`a`) and of the parameter `t`
/// (`b`): a quadric's `G` itself, or another torus's (S9d.4b.2b) `Q = |w|^2`
/// and `P = (w . x2)^2 + (w . y2)^2` with its constants `k = R2^2 - r2^2`
/// and `4 R2^2`, `G = (Q + k)^2 - 4 R2^2 P`. Each of the other surface's
/// affine functionals of the torus's point is `alpha + beta z + rho (xi cos
/// u + eta sin u)` (`rho = R + r cos v`, `z = r sin v`), and its square's
/// constants are summed over the functionals before any jet: a pair coaxial
/// within rounding has its terms in `u` of rounding's size, so its jets
/// along a meeting near a double root do not carry the functionals' own
/// terms in `u` (of the torus's size, cancelling in their squares' sum only
/// exactly) through every order's division by the small slope there.
struct Harmonic<T> {
    /// Each form's `w[a][b]`.
    forms: Vec<[[T; 5]; 5]>,
    /// Another torus's `R2^2 - r2^2` and `4 R2^2`.
    torus: Option<(T, T)>,
}

impl<T: Real> Harmonic<T> {
    fn new(m: &Toric) -> Option<Self> {
        let zero = || T::exact_f64(0.0);
        let two = T::exact_f64(2.0);
        let (o, x, y, nn) = (
            m.frame.origin().to_array(),
            m.frame.x().to_array(),
            m.frame.y().to_array(),
            m.frame.normal().to_array(),
        );
        let o2 = m.other.origin().to_array();
        // The origins' offset first (exact where their coordinates are
        // close).
        let d: [T; 3] = std::array::from_fn(|k| c::<T>(o[k]).sub(&c(o2[k])));
        let dot =
            |a: &[T; 3], f: &[f64; 3]| (0..3).fold(zero(), |acc, k| acc.add(&a[k].mul(&c(f[k]))));
        let (xs, ys, ns) = (x.map(c::<T>), y.map(c::<T>), nn.map(c::<T>));
        // The functional `a + k (w . f)` as `[alpha, beta, xi, eta]`.
        let functional = |f: [f64; 3], k: &T, a: &T| -> [T; 4] {
            [
                a.add(&dot(&d, &f).mul(k)),
                dot(&ns, &f).mul(k),
                dot(&xs, &f).mul(k),
                dot(&ys, &f).mul(k),
            ]
        };
        let small = c::<T>(m.minor);
        let rho = [c::<T>(m.major), small.clone(), zero()];
        let rho2 = trig_product(&rho, &rho);
        // `sum sign L^2` as `g[i][j] D_i(u) D_j(v)`: `L = A(v) + rho(v)
        // E(u)`, `A = alpha + beta r sin v`, `E = xi cos u + eta sin u`, so
        // `L^2 = A^2 + 2 A rho E + rho^2 E^2`, the `E^2` summed first.
        let form = |fs: &[([T; 4], bool)]| -> [[T; 5]; 5] {
            let mut g: [[T; 5]; 5] = std::array::from_fn(|_| std::array::from_fn(|_| zero()));
            let mut e2: [T; 5] = std::array::from_fn(|_| zero());
            for ([alpha, beta, xi, eta], plus) in fs {
                let acc = |x: &T, y: &T| if *plus { x.add(y) } else { x.sub(y) };
                let a = [alpha.clone(), zero(), beta.mul(&small)];
                let e = [zero(), xi.clone(), eta.clone()];
                let (a2, ar, ee) = (
                    trig_product(&a, &a),
                    trig_product(&a, &rho),
                    trig_product(&e, &e),
                );
                for j in 0..5 {
                    g[0][j] = acc(&g[0][j], &a2[j]);
                    e2[j] = acc(&e2[j], &ee[j]);
                    g[1][j] = acc(&g[1][j], &xi.mul(&ar[j]).mul(&two));
                    g[2][j] = acc(&g[2][j], &eta.mul(&ar[j]).mul(&two));
                }
            }
            for (row, e) in g.iter_mut().zip(&e2) {
                for (x, r) in row.iter_mut().zip(&rho2) {
                    *x = x.add(&e.mul(r));
                }
            }
            // `[a][b]`: the other angle's harmonics first.
            if m.over_v {
                g
            } else {
                std::array::from_fn(|a| std::array::from_fn(|b| g[b][a].clone()))
            }
        };
        let one = T::exact_f64(1.0);
        let axes = [m.other.x(), m.other.y(), m.other.normal()].map(|f| f.to_array());
        if m.other_minor > 0.0 {
            let fs = axes.map(|f| (functional(f, &one, &zero()), true));
            let big2 = c::<T>(m.other_radius).square();
            let small2 = c::<T>(m.other_minor).square();
            return Some(Self {
                forms: vec![form(&fs), form(&fs[..2])],
                torus: Some((big2.sub(&small2), big2.mul(&T::exact_f64(4.0)))),
            });
        }
        let r = c::<T>(m.other_radius);
        let constant = |a: T| ([a, zero(), zero(), zero()], false);
        let fs = if m.other_sphere {
            let mut fs: Vec<([T; 4], bool)> = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
                .iter()
                .map(|f| (functional(*f, &one, &zero()), true))
                .collect();
            fs.push(constant(r));
            fs
        } else {
            let mut fs = vec![
                (functional(axes[0], &one, &zero()), true),
                (functional(axes[1], &one, &zero()), true),
            ];
            if m.other_half_angle == 0.0 {
                fs.push(constant(r));
            } else {
                // The other's radius along its axis: `r2 + t2 (w . n2)`.
                let (ca, sa) = T::cos_sin(&c(m.other_half_angle));
                fs.push((functional(axes[2], &sa.div(&ca)?, &r), false));
            }
            fs
        };
        Some(Self {
            forms: vec![form(&fs)],
            torus: None,
        })
    }

    /// Each form's coefficients of the other angle's harmonics at the
    /// parameter's: `kappa[a] = sum_b w[a][b] D_b(t)`.
    fn kappa(&self, dt: &Harmonics<T>) -> Vec<Harmonics<T>> {
        self.forms
            .iter()
            .map(|w| {
                std::array::from_fn(|a| {
                    (1..5).fold(w[a][0].clone(), |acc, b| acc.add(&w[a][b].mul(&dt[b])))
                })
            })
            .collect()
    }

    /// Each form's value at the other angle's harmonics (or their rates).
    fn values(kappa: &[Harmonics<T>], ds: &Harmonics<T>) -> Vec<T> {
        kappa
            .iter()
            .map(|k| (1..5).fold(k[0].mul(&ds[0]), |acc, a| acc.add(&k[a].mul(&ds[a]))))
            .collect()
    }

    /// `G` from its forms' values.
    fn g(&self, v: &[T]) -> T {
        match &self.torus {
            Some((k, four)) => v[0].add(k).square().sub(&four.mul(&v[1])),
            None => v[0].clone(),
        }
    }

    /// `G`'s derivative from its forms' values and derivatives.
    fn rate(&self, v: &[T], d: &[T]) -> T {
        match &self.torus {
            Some((k, four)) => v[0]
                .add(k)
                .mul(&d[0])
                .mul(&T::exact_f64(2.0))
                .sub(&four.mul(&d[1])),
            None => d[0].clone(),
        }
    }

    /// `G` and its partials at enclosed angles `s` and `t`: `[G, G_s, G_t,
    /// G_ss, G_st, G_tt]` (the mean-value forms over a wide base).
    fn partials(&self, s: &T, t: &T) -> [T; 6] {
        let (ds, dt) = (harmonics(s), harmonics(t));
        let (ds1, dt1) = (harmonics_rate(&ds), harmonics_rate(&dt));
        let (ds2, dt2) = (harmonics_rate(&ds1), harmonics_rate(&dt1));
        let at = |w: &[[T; 5]; 5], a: &Harmonics<T>, b: &Harmonics<T>| {
            (0..5).fold(T::exact_f64(0.0), |acc, i| {
                let row = (1..5).fold(w[i][0].mul(&b[0]), |r, j| r.add(&w[i][j].mul(&b[j])));
                acc.add(&row.mul(&a[i]))
            })
        };
        let f: Vec<[T; 6]> = self
            .forms
            .iter()
            .map(|w| {
                [
                    at(w, &ds, &dt),
                    at(w, &ds1, &dt),
                    at(w, &ds, &dt1),
                    at(w, &ds2, &dt),
                    at(w, &ds1, &dt1),
                    at(w, &ds, &dt2),
                ]
            })
            .collect();
        let Some((k, four)) = &self.torus else {
            return f[0].clone();
        };
        // `(Q + k)^2 - 4 R2^2 P` and its partials.
        let two = T::exact_f64(2.0);
        let (q, p) = (&f[0], &f[1]);
        let qk = q[0].add(k);
        let first = |i: usize| qk.mul(&q[i]).mul(&two).sub(&four.mul(&p[i]));
        let second = |i: usize, a: usize, b: usize| {
            q[a].mul(&q[b])
                .add(&qk.mul(&q[i]))
                .mul(&two)
                .sub(&four.mul(&p[i]))
        };
        [
            qk.square().sub(&four.mul(&p[0])),
            first(1),
            first(2),
            second(3, 1, 1),
            second(4, 1, 2),
            second(5, 2, 2),
        ]
    }
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

/// A torus meeting's jets kept in the binary64 tier: the meeting's bits, and
/// the variable's base's bits and order.
type ToricKey = (Vec<u64>, u64, u64, usize);

/// Torus meetings' jets kept before the memo starts again.
const TORIC_LIMIT: usize = 1 << 15;

thread_local! {
    /// The jets of torus meetings (`toric_jet_of`) in the binary64 tier, of
    /// the variable `B + s` about a base `B`, by the meeting's every binary64
    /// number and flag and the base's bits and the order: an edge's two
    /// faces' pcurves (on its own torus and on the other surface), the
    /// validator's signs, areas and fluxes, its measures and the mass's
    /// moments integrate along the same edge and visit the same pieces of
    /// its fraction. A hit is the value the evaluation gives.
    static TORIC: std::cell::RefCell<std::collections::HashMap<ToricKey, Option<ToricJet<Fast>>>> =
        std::cell::RefCell::new(std::collections::HashMap::new());
}

/// A torus meeting's binary64 numbers' and flags' bits.
fn toric_bits(m: &Toric) -> Vec<u64> {
    let mut out = Vec::with_capacity(34);
    for f in [&m.frame, &m.other] {
        for v in [
            f.origin().to_array(),
            f.x().to_array(),
            f.y().to_array(),
            f.normal().to_array(),
        ] {
            out.extend(v.iter().map(|x| x.to_bits()));
        }
    }
    out.extend(
        [
            m.major,
            m.minor,
            m.other_radius,
            m.other_half_angle,
            m.other_minor,
            m.window[0],
            m.window[1],
            m.start,
            m.sweep,
        ]
        .iter()
        .map(|x| x.to_bits()),
    );
    out.push(u64::from(m.other_sphere) | (u64::from(m.over_v) << 1));
    out
}

/// `toric_jet_of` in the binary64 tier through its memo (`TORIC`), of a
/// variable `B + s` or `B - s` about a base `B` (a pcurve's fraction on a
/// reversed use: `1 - f`): the latter's jets are the former's in `-s`, each
/// odd coefficient negated. Binary64 intervals round outward to nearest on
/// both sides alike, so the recurrences give the same values to the same
/// bits with every odd term negated: the evaluation's own result (the two
/// fins of an edge on its two faces run opposite ways). Any other jet, or
/// another tier, is evaluated as it is.
fn toric_jet<T: Real>(m: &Toric, fraction: &Jet<T>) -> Option<ToricJet<T>> {
    let exact = |x: &T, v: f64| {
        x.as_fast().is_some_and(|x| {
            let (lo, hi) = x.bounds_f64();
            lo == v && hi == v
        })
    };
    let (Some(base), true) = (
        fraction.c[0].as_fast(),
        fraction.order() >= 1 && fraction.c[2..].iter().all(|x| exact(x, 0.0)),
    ) else {
        return toric_jet_of(m, fraction);
    };
    let reflected = if exact(&fraction.c[1], 1.0) {
        false
    } else if exact(&fraction.c[1], -1.0) {
        true
    } else {
        return toric_jet_of(m, fraction);
    };
    let order = fraction.order();
    let (lo, hi) = base.bounds_f64();
    let key = (toric_bits(m), lo.to_bits(), hi.to_bits(), order);
    let hit = TORIC.with(|t| t.borrow().get(&key).cloned());
    let value = match hit {
        Some(value) => value,
        None => {
            let value = toric_jet_of(m, &Jet::variable(base, order));
            TORIC.with(|t| {
                let mut t = t.borrow_mut();
                if t.len() >= TORIC_LIMIT {
                    t.clear();
                }
                t.insert(key, value.clone());
            });
            value
        }
    };
    let to_t = |j: &Jet<Fast>| -> Option<Jet<T>> {
        Some(Jet {
            c: j.c
                .iter()
                .enumerate()
                .map(|(k, x)| T::of_fast(if reflected && k % 2 == 1 { x.neg() } else { *x }))
                .collect::<Option<_>>()?,
        })
    };
    let ([u, v], p) = value.as_ref()?;
    let out = (
        [to_t(u)?, to_t(v)?],
        [to_t(&p[0])?, to_t(&p[1])?, to_t(&p[2])?],
    );
    Some(out)
}

/// The jets of a torus's meeting with a quadric or another torus (S9d.4b.2)
/// in the fraction, from `G`'s harmonic forms (`Harmonic`): the other angle
/// `s` enclosed over the base by interval Newton inside the window (a unique
/// root of `G` there), then its coefficients by the implicit function
/// theorem, term by term: the `k`-th coefficient of `G` along the meeting is
/// its terms without the other angle's `k`-th harmonics' coefficients
/// (`z`), plus those, which are linear in `C_k` and `S_k` (the cosine's and
/// sine's) with `G`'s partials at the base; and `C_k = C P - S q`, `S_k = S
/// P + C q` with `P = C C_k + S S_k` from the lower orders (`C^2 + S^2 =
/// 1`) and `q = C S_k - S C_k`, so `q = -(z + K P) / G_s` (`K = C G_C + S
/// G_S`, `G_s = C G_S - S G_C`) and `s_k = q - (C S_k' - S C_k')` (`C_k'`,
/// `S_k'` continued with `s_k` zero): each term taken once, an enclosure at
/// every point of the base (inclusion isotone). Over a wide base, `G` at
/// `s*`, `G_s` and the first coefficient are also taken in mean-value forms
/// about the base's middle, the narrower of two enclosures of one value.
fn toric_jet_of<T: Real>(m: &Toric, fraction: &Jet<T>) -> Option<ToricJet<T>> {
    let n = fraction.order();
    let t = fraction.scale(&c(m.sweep)).add_constant(&c(m.start));
    let h = Harmonic::<T>::new(m)?;
    let two = T::exact_f64(2.0);
    let zero = || T::exact_f64(0.0);
    // The parameter's harmonics' jets (the first, `1`, constant).
    let (ct, st) = t.cos_sin();
    let (ct2, st2) = t.scale(&two).cos_sin();
    let tj = [ct, st, ct2, st2];
    let dt: Harmonics<T> = [
        T::exact_f64(1.0),
        tj[0].c[0].clone(),
        tj[1].c[0].clone(),
        tj[2].c[0].clone(),
        tj[3].c[0].clone(),
    ];
    let kappa = h.kappa(&dt);
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
    let (lo, hi) = t.c[0].bounds_f64();
    let star = m.root_at(0.5 * lo + 0.5 * hi);
    if !star.is_finite() {
        return None;
    }
    let s0 = c::<T>(star);
    // Interval Newton over the base: `G` at `s*` over all of it, and `G_s`
    // over the base and the box; over a wide base each also in its
    // mean-value form about `(t_mid, s*)` (another torus's `(Q + k)^2 - 4
    // R2^2 P` cancels over it), the narrower holding it.
    let point = lo == hi;
    let mid = c::<T>(0.5 * lo + 0.5 * hi);
    let spread = t.c[0].sub(&mid);
    let at_mid = h.partials(&s0, &mid);
    let g_star = {
        let d = harmonics(&s0);
        let natural = h.g(&Harmonic::values(&kappa, &d));
        if point {
            natural
        } else {
            let over = h.partials(&s0, &t.c[0]);
            narrower(natural, at_mid[0].add(&over[2].mul(&spread)))
        }
    };
    let slope0 = least(&at_mid[1]);
    if slope0 <= 0.0 || slope0.is_nan() {
        return None;
    }
    // `G_s` over the base at angles `v` about `s*`: natural, or with its
    // mean-value form.
    let rate_over = |v: &T| {
        let d = harmonics(v);
        let values = Harmonic::values(&kappa, &d);
        let natural = h.rate(&values, &Harmonic::values(&kappa, &harmonics_rate(&d)));
        if point {
            return (natural, None);
        }
        let b = h.partials(v, &t.c[0]);
        let form = at_mid[1]
            .add(&b[3].mul(&v.sub(&s0)))
            .add(&b[4].mul(&spread));
        (narrower(natural, form), Some(b))
    };
    let mut delta = 4.0 * size(&g_star) / slope0 + 1e-15 * (1.0 + star.abs());
    let mut root = None;
    for _ in 0..12 {
        if !delta.is_finite() || delta > 0.5 {
            return None;
        }
        let (rate, _) = rate_over(&c::<T>(star - delta).union(&c(star + delta)));
        if let Some(q) = g_star.div(&rate) {
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
    let ds = harmonics(&s_base);
    let [_, c0, si0, c20, s20] = ds.clone();
    // Each form's `G_C` and `G_S` combined: `K` and `G_s` per form, then
    // `G`'s through another torus's `(Q + k)^2 - 4 R2^2 P`.
    let per_form = |k: &Harmonics<T>| {
        let kk = k[1]
            .mul(&c0)
            .add(&k[2].mul(&si0))
            .add(&k[3].mul(&c20).mul(&two))
            .add(&k[4].mul(&s20).mul(&two));
        let gs = k[2]
            .mul(&c0)
            .sub(&k[1].mul(&si0))
            .sub(&k[3].mul(&s20).mul(&two))
            .add(&k[4].mul(&c20).mul(&two));
        (kk, gs)
    };
    let parts: Vec<(T, T)> = kappa.iter().map(per_form).collect();
    // The forms' values along the meeting (another torus's `Q` series).
    let mut values: Vec<Vec<T>> = Harmonic::values(&kappa, &ds)
        .into_iter()
        .map(|v| vec![v])
        .collect();
    // `G`'s rate in each form's value (`2 (Q + k)` and `-4 R2^2`).
    let weights: Vec<T> = match &h.torus {
        Some((k, four)) => vec![values[0][0].add(k).mul(&two), four.neg()],
        None => vec![T::exact_f64(1.0)],
    };
    let combine = |xs: &[T]| {
        xs.iter()
            .zip(&weights)
            .fold(zero(), |acc, (x, w)| acc.add(&w.mul(x)))
    };
    let kk = combine(&parts.iter().map(|p| p.0.clone()).collect::<Vec<_>>());
    let slope = combine(&parts.iter().map(|p| p.1.clone()).collect::<Vec<_>>());
    // Over a wide base, `G_s` along the meeting and its first coefficient's
    // rest (`G_t t'`) also in their mean-value forms about `(t_mid, s*)`,
    // over the box from there to every point of the meeting over the base
    // (`s*` is binary64's root, not quite the meeting's).
    let (slope, first) = match rate_over(&s_base.union(&s0)) {
        (_, Some(b)) => {
            let spread_s = s_base.sub(&s0);
            let slope = narrower(
                slope,
                at_mid[1].add(&b[3].mul(&spread_s)).add(&b[4].mul(&spread)),
            );
            let first = t.c.get(1).map(|rate| {
                at_mid[2]
                    .add(&b[5].mul(&spread))
                    .add(&b[4].mul(&spread_s))
                    .mul(rate)
            });
            (slope, first)
        }
        _ => (slope, None),
    };
    let ks = kk.div(&slope)?;
    // Each form's `K_b = sum_a w[a][b] D_a(s)` along the meeting, for the
    // parameter's harmonics `b` past the first.
    let mut kb: Vec<[Vec<T>; 4]> = h
        .forms
        .iter()
        .map(|w| {
            std::array::from_fn(|b| {
                vec![(1..5).fold(w[0][b + 1].clone(), |acc, a| {
                    acc.add(&w[a][b + 1].mul(&ds[a]))
                })]
            })
        })
        .collect();
    let mut s = vec![s_base];
    let (mut co, mut sn) = (vec![c0.clone()], vec![si0.clone()]);
    for k in 1..=n {
        let order = T::exact_f64(k as f64);
        // cos and sin continued with s_k = 0: k S_k = sum C_i d_{k-1-i},
        // k C_k = -sum S_i d_{k-1-i}, d_j = (j + 1) s_{j+1} (i = 0 holds
        // s_k).
        let (mut a_s, mut a_c) = (zero(), zero());
        for i in 1..k {
            let d = s[k - i].mul(&T::exact_f64((k - i) as f64));
            a_s = a_s.add(&co[i].mul(&d));
            a_c = a_c.add(&sn[i].mul(&d));
        }
        let (ca, sb) = (a_c.neg().div(&order)?, a_s.div(&order)?);
        let qq = c0.mul(&sb).sub(&si0.mul(&ca));
        // The lower orders' parts of `C^2 + S^2`, `C^2 - S^2` and `2 C S`.
        let (mut pp, mut m2c, mut m2s) = (zero(), zero(), zero());
        for i in 1..k {
            let (cc, ss) = (co[i].mul(&co[k - i]), sn[i].mul(&sn[k - i]));
            pp = pp.add(&cc).add(&ss);
            m2c = m2c.add(&cc).sub(&ss);
            m2s = m2s.add(&co[i].mul(&sn[k - i]));
        }
        let pp = pp.mul(&T::exact_f64(-0.5));
        let m2s = m2s.mul(&two);
        // Each form's `k`-th coefficient without the other angle's `k`-th
        // harmonics' (`z`), and with their lower orders' parts (`z2`).
        let zs: Vec<T> = kb
            .iter()
            .map(|ks| {
                let mut z = zero();
                for (b, series) in ks.iter().enumerate() {
                    for (i, x) in series.iter().enumerate() {
                        z = z.add(&tj[b].c[k - i].mul(x));
                    }
                }
                z
            })
            .collect();
        let z2s: Vec<T> = zs
            .iter()
            .zip(&kappa)
            .map(|(z, kp)| z.add(&kp[3].mul(&m2c)).add(&kp[4].mul(&m2s)))
            .collect();
        let mut z = combine(&z2s);
        if let (Some(form), 1) = (&first, k) {
            z = narrower(z, form.clone());
        }
        if h.torus.is_some() {
            // `Q`'s own lower orders in `(Q + k)^2`.
            let q = &values[0];
            for i in 1..k {
                z = z.add(&q[i].mul(&q[k - i]));
            }
        }
        let q = z.div(&slope)?.add(&ks.mul(&pp)).neg();
        let ck = c0.mul(&pp).sub(&si0.mul(&q));
        let snk = si0.mul(&pp).add(&c0.mul(&q));
        let c2k = c20.mul(&pp).sub(&s20.mul(&q)).mul(&two).add(&m2c);
        let s2k = s20.mul(&pp).add(&c20.mul(&q)).mul(&two).add(&m2s);
        let dk = [zero(), ck.clone(), snk.clone(), c2k, s2k];
        for (f, w) in h.forms.iter().enumerate() {
            for (b, series) in kb[f].iter_mut().enumerate() {
                series.push((1..5).fold(zero(), |acc, a| acc.add(&w[a][b + 1].mul(&dk[a]))));
            }
            let own = (1..5).fold(zs[f].clone(), |acc, a| acc.add(&kappa[f][a].mul(&dk[a])));
            values[f].push(own);
        }
        s.push(q.sub(&qq));
        co.push(ck);
        sn.push(snk);
    }
    let s = Jet { c: s };
    let (cj, sj) = (Jet { c: co }, Jet { c: sn });
    let p = toric_point_parts(m, &t);
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
        let reference = s.plane[1].atan2_rn(s.plane[0]);
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
            s.plane[0] * base.cos_rn() + s.plane[1] * base.sin_rn(),
            s.plane[2],
        );
        let psi = angle_near(&y, &x, (s.minor * cf).atan2_rn(s.minor * ab))?;
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

/// The product of two affine forms `k0 + kc cos u + ks sin u` as a form of
/// degree two: its constants of `1`, `cos u`, `sin u`, `cos 2u`, `sin 2u`.
fn trig_product<T: Real>(p: &[T; 3], q: &[T; 3]) -> [T; 5] {
    let half = c::<T>(0.5);
    let (cc, ss) = (p[1].mul(&q[1]), p[2].mul(&q[2]));
    let (cs, sc) = (p[1].mul(&q[2]), p[2].mul(&q[1]));
    [
        p[0].mul(&q[0]).add(&cc.add(&ss).mul(&half)),
        p[0].mul(&q[1]).add(&p[1].mul(&q[0])),
        p[0].mul(&q[2]).add(&p[2].mul(&q[0])),
        cc.sub(&ss).mul(&half),
        cs.add(&sc).mul(&half),
    ]
}

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
    // `foot - o2` and the direction along an axis, as affine forms in the
    // angle, `k0 + kc cos u + ks sin u` (interval constants); `a`, `b` and
    // `c` are sums of their products, forms of degree two in `u` whose
    // constants are combined before any jet: the same functions, but a
    // jet of `cos^2 + sin^2` over a piece is wider than one, and a coaxial
    // pair's forms (constant in `u`) took thousands of pieces a turn.
    let o1 = m.frame.origin().to_array();
    let dot = |p: &[f64; 3], e: &[f64; 3]| {
        (0..3).fold(c::<T>(0.0), |acc, k| acc.add(&c::<T>(p[k]).mul(&c(e[k]))))
    };
    let foot_on = |e: &[f64; 3]| -> [T; 3] {
        let k0 = (0..3).fold(c::<T>(0.0), |acc, k| {
            acc.add(&c::<T>(o1[k]).sub(&c(o2[k])).mul(&c(e[k])))
        });
        [k0, r.mul(&dot(&x, e)), r.mul(&dot(&y, e))]
    };
    let dir_on = |e: &[f64; 3]| -> [T; 3] { [dot(&n, e), t.mul(&dot(&x, e)), t.mul(&dot(&y, e))] };
    let axes: Vec<[f64; 3]> = if m.other_sphere {
        vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    } else {
        vec![m.other.x().to_array(), m.other.y().to_array()]
    };
    let zero = || [0; 5].map(|_| c::<T>(0.0));
    let (mut a, mut b, mut cc) = (zero(), zero(), zero());
    let add = |acc: &mut [T; 5], p: [T; 5], plus: bool| {
        for (x, y) in acc.iter_mut().zip(p) {
            *x = if plus { x.add(&y) } else { x.sub(&y) };
        }
    };
    for e in &axes {
        let (w, d) = (foot_on(e), dir_on(e));
        add(&mut a, trig_product(&d, &d), true);
        add(&mut b, trig_product(&w, &d), true);
        add(&mut cc, trig_product(&w, &w), true);
    }
    if m.other_sphere {
        cc[0] = cc[0].sub(&c::<T>(m.other_radius).mul(&c(m.other_radius)));
    } else {
        // The other's radius along its axis: `r2 + t2 (w . n2)` (a cone).
        let n2 = m.other.normal().to_array();
        let t2 = tan(m.other_half_angle)?;
        let [k0, kc, ks] = foot_on(&n2).map(|k| k.mul(&t2));
        let r0 = [k0.add(&c(m.other_radius)), kc, ks];
        let rd = dir_on(&n2).map(|k| k.mul(&t2));
        add(&mut a, trig_product(&rd, &rd), false);
        add(&mut b, trig_product(&r0, &rd), false);
        add(&mut cc, trig_product(&r0, &r0), false);
    }
    let (co2, si2) = u.scale(&c(2.0)).cos_sin();
    let jet = |k: &[T; 5]| {
        co.scale(&k[1])
            .add(&si.scale(&k[2]))
            .add(&co2.scale(&k[3]))
            .add(&si2.scale(&k[4]))
            .add_constant(&k[0])
    };
    // `d = b^2 - a c` over a piece also in its mean-value form about the
    // angle's middle, from the forms' values there and their derivatives
    // over the piece (the natural extension's terms cancel: wide by their
    // own size, which the square root's recurrence amplified to millions of
    // times the coefficients' size).
    let (ulo, uhi) = u.c[0].bounds_f64();
    let d_mid = if ulo < uhi {
        let um = 0.5 * ulo + 0.5 * uhi;
        let value = |k: &[T; 5], x: &T| {
            let ((c1, s1), (c2, s2)) = (T::cos_sin(x), T::cos_sin(&x.mul(&c(2.0))));
            k[0].add(&k[1].mul(&c1))
                .add(&k[2].mul(&s1))
                .add(&k[3].mul(&c2))
                .add(&k[4].mul(&s2))
        };
        let slope = |k: &[T; 5], x: &T| {
            let ((c1, s1), (c2, s2)) = (T::cos_sin(x), T::cos_sin(&x.mul(&c(2.0))));
            k[2].mul(&c1)
                .sub(&k[1].mul(&s1))
                .add(&k[4].mul(&c2).sub(&k[3].mul(&s2)).mul(&c(2.0)))
        };
        let (m, whole) = (c::<T>(um), u.c[0].clone());
        let d_at = value(&b, &m)
            .square()
            .sub(&value(&a, &m).mul(&value(&cc, &m)));
        let (bw, aw, cw) = (value(&b, &whole), value(&a, &whole), value(&cc, &whole));
        let d_u = bw
            .mul(&slope(&b, &whole))
            .mul(&c(2.0))
            .sub(&slope(&a, &whole).mul(&cw))
            .sub(&aw.mul(&slope(&cc, &whole)));
        Some(d_at.add(&d_u.mul(&whole.sub(&m))))
    } else {
        None
    };
    // The same functions as the plain products of the dot products' jets
    // too (narrower where the forms' constants are wide against their
    // values: near a loop's turning point), each coefficient the narrower.
    let (a, b, cc) = {
        let along = |v: &[Jet<T>; 3], e: &[f64; 3], shift: bool| {
            (0..3).fold(Jet::constant(c::<T>(0.0), u.order()), |acc, k| {
                let vk = if shift {
                    v[k].add_constant(&c::<T>(o2[k]).neg())
                } else {
                    v[k].clone()
                };
                acc.add(&vk.scale(&c(e[k])))
            })
        };
        let zero = || Jet::constant(c::<T>(0.0), u.order());
        let (mut pa, mut pb, mut pc) = (zero(), zero(), zero());
        for e in &axes {
            let (w, d) = (along(&foot, e, true), along(&dir, e, false));
            pa = pa.add(&d.square());
            pb = pb.add(&w.mul(&d));
            pc = pc.add(&w.square());
        }
        if m.other_sphere {
            pc = pc.add_constant(&c::<T>(m.other_radius).mul(&c(m.other_radius)).neg());
        } else {
            let n2 = m.other.normal().to_array();
            let t2 = tan(m.other_half_angle)?;
            let r0 = along(&foot, &n2, true)
                .scale(&t2)
                .add_constant(&c(m.other_radius));
            let rd = along(&dir, &n2, false).scale(&t2);
            pa = pa.sub(&rd.square());
            pb = pb.sub(&r0.mul(&rd));
            pc = pc.sub(&r0.square());
        }
        let pick = |x: Jet<T>, y: Jet<T>| Jet {
            c: x.c
                .into_iter()
                .zip(y.c)
                .map(|(p, q)| narrower(p, q))
                .collect(),
        };
        (pick(jet(&a), pa), pick(jet(&b), pb), pick(jet(&cc), pc))
    };
    let mut d = b.square().sub(&cc.mul(&a));
    if let Some(at) = d_mid {
        d.c[0] = narrower(d.c[0].clone(), at);
    }
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
///
/// `acos(q)`, `q = g / s` with `s = rho radius(w)`, is the angle of `(g,
/// sqrt((s - g) (s + g)))`, each of `g`, `s - g` and `s + g` a quadratic in
/// `w` taken about the middle of the piece's heights: on a loop shallower
/// than its size (a sphere crossing the cylinder by `1e-6` of its radius)
/// `|q|` stays within about that depth of one, and `g` taken in `w` far
/// from the frame's origin spreads over a piece by its terms' slopes, not
/// its own, so `|q|` passed one on all but pieces of `2^-26` of the fraction
/// and the integrals along it ran for an hour (the same function, its
/// enclosures narrower).
fn rise_jet<T: Real>(m: &Rise, fraction: &Jet<T>) -> Option<[Jet<T>; 3]> {
    let w = fraction.scale(&c(m.sweep)).add_constant(&c(m.start));
    let ([a, b], g) = m.coefficients();
    let rho = c::<T>(a).square().add(&c::<T>(b).square()).sqrt();
    // The carrier's radius at the height (a cone's varies, S9d.3b.2).
    let (ca, sa) = T::cos_sin(&c(m.half_angle));
    let slope = sa.div(&ca)?;
    let radius = w.scale(&slope).add_constant(&c(m.radius));
    // Each quadratic `k0 + k1 w + k2 w^2` about the heights' middle `wm`.
    let (lo, hi) = w.c[0].bounds_f64();
    let wm = c::<T>(0.5 * lo + 0.5 * hi);
    let dw = w.add_constant(&wm.neg());
    let about = |k: [T; 3]| {
        let k0 = k[0].add(&k[1].mul(&wm)).add(&k[2].mul(&wm).mul(&wm));
        let k1 = k[1].add(&k[2].mul(&wm).mul(&c(2.0)));
        dw.square()
            .scale(&k[2])
            .add(&dw.scale(&k1))
            .add_constant(&k0)
    };
    let g = g.map(c::<T>);
    let s = [rho.mul(&c(m.radius)), rho.mul(&slope), c(0.0)];
    let gq = about(g.clone());
    let minus = about([0, 1, 2].map(|k| s[k].sub(&g[k])));
    let plus = about([0, 1, 2].map(|k| s[k].add(&g[k])));
    let y = minus.mul(&plus).sqrt()?;
    // `q = g / s` is defined where `s` certainly is not zero (the carrier's
    // apex, a sphere centred on its axis: none), its sign taken into `g`.
    let x = match radius.scale(&rho).c[0].sign() {
        Some(std::cmp::Ordering::Greater) => gq,
        Some(std::cmp::Ordering::Less) => gq.neg(),
        _ => return None,
    };
    // `phi` is the curve's binary64 constant (its definition's `atan2`).
    let phi = c::<T>(b.atan2_rn(a));
    let u = angle_near(&y, &x, std::f64::consts::FRAC_PI_2)?
        .scale(&c(m.sign))
        .add_constant(&phi);
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
    /// A spline wall's meeting's pieces' jets (`wall_piece_jet`), keyed by
    /// the knot span too.
    walls: std::collections::HashMap<WallKey, Option<[Jet<Fast>; 2]>>,
}

/// A wall piece's key: the projection's id, the knot span, the variable's
/// base's bits and the order.
type WallKey = (u64, usize, u64, u64, usize);

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
        walls: std::collections::HashMap::new(),
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
            m.walls.clear();
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
    let (width, floor) = if narrow() {
        (NARROW_WIDTH, Some(NARROW_FLOOR))
    } else if relative {
        (WIDTH, Some(1.0))
    } else {
        (SIGN_WIDTH.with(|w| w.get()), None)
    };
    integrate_many(&integrand, n, 0.0, 1.0, ORDER, width, DEPTH, floor)
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
    let id = memo_id(p);
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
            let [u, v] = wall_piece_jet(p, id, k, &tau.c[0], tau.order(), &at)?;
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
        let values = integrate_many(
            &integrand,
            n,
            0.0,
            1.0,
            ORDER,
            width,
            DEPTH,
            relative.then_some(1.0),
        )?;
        for (t, x) in totals.iter_mut().zip(values) {
            *t = t.add(&x);
        }
    }
    Some(())
}

/// A spline wall's meeting's projection's jets on knot span `k`'s piece at
/// `at` (the piece's fraction, the variable about `base` to `order`),
/// through the memo in the binary64 tier: the validator's signs and areas
/// and the mass's moments and fluxes integrate along the same pieces, each
/// to its own width, and visit the same dyadic parts of them.
fn wall_piece_jet<T: Real>(
    p: &Projection,
    id: Option<u64>,
    k: usize,
    base: &T,
    order: usize,
    at: &Jet<T>,
) -> Option<[Jet<T>; 2]> {
    let fresh = || projection_jet_pinned(p, at, Some(k));
    let (Some(id), Some(fast)) = (id, base.as_fast()) else {
        return fresh();
    };
    let (lo, hi) = fast.bounds_f64();
    let key = (id, k, lo.to_bits(), hi.to_bits(), order);
    let to_t = |j: &Jet<Fast>| -> Option<Jet<T>> {
        Some(Jet {
            c: j.c.iter().map(|x| T::of_fast(*x)).collect::<Option<_>>()?,
        })
    };
    if let Some(hit) = MEMO.with(|m| m.borrow().walls.get(&key).cloned()) {
        let [u, v] = hit?;
        return Some([to_t(&u)?, to_t(&v)?]);
    }
    let value = fresh();
    let to_fast = |j: &Jet<T>| -> Option<Jet<Fast>> {
        Some(Jet {
            c: j.c.iter().map(|x| x.as_fast()).collect::<Option<_>>()?,
        })
    };
    let stored = match &value {
        None => Some(None),
        Some([u, v]) => match (to_fast(u), to_fast(v)) {
            (Some(u), Some(v)) => Some(Some([u, v])),
            _ => None,
        },
    };
    if let Some(stored) = stored {
        MEMO.with(|m| {
            let mut m = m.borrow_mut();
            if m.walls.len() >= MEMO_LIMIT {
                m.walls.clear();
            }
            m.walls.insert(key, stored);
        });
    }
    value
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
        if !crate::budget::unit_ok("a ray crossing's piece") {
            return None;
        }
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
    cover_crossings_on(pr, p, false)
}

/// `cover_crossings` with the parameters' roles exchanged (`swap`): the
/// `+u` ray's crossings over every `v` alias, `p` given as `(v, u)` (a torus
/// face wound in `v`, S9e.3b).
pub(super) fn cover_crossings_on<T: Real>(
    pr: &Projection,
    p: &V2<T>,
    swap: bool,
) -> Option<Vec<i64>> {
    use std::cmp::Ordering;
    use std::f64::consts::TAU;
    if T::EXACT {
        return None;
    }
    let at = |f: f64| -> Option<V2<T>> {
        let [a, b] = projection_at::<T>(pr, f)?;
        Some(if swap { [b, a] } else { [a, b] })
    };
    let mut out = Vec::new();
    let mut stack = vec![(0.0f64, 1.0f64, 0usize)];
    while let Some((lo, hi, depth)) = stack.pop() {
        if !crate::budget::unit_ok("a ray crossing's piece") {
            return None;
        }
        let mid = 0.5 * lo + 0.5 * hi;
        let base = T::exact_f64(lo).union(&T::exact_f64(hi));
        let decided = (|| -> Option<Vec<i64>> {
            let [u, v] = projection_jet(pr, &Jet::variable(base, 1))?;
            let [u, v] = if swap { [v, u] } else { [u, v] };
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
                let (a, b) = (at(lo)?, at(hi)?);
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
                assert!((x.hypot_rn(y) - 2.0).abs() < 1e-14, "{p:?}");
                let [_, b, c] = [p.x, p.y - 0.5, p.z];
                assert!((b.hypot_rn(c) - 1.0).abs() < 1e-14);
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
                assert!((p.x.hypot_rn(p.y) - 1.0).abs() < 1e-14);
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

    /// A ball of radius 0.75 crossing a rod of radius 1 (its frame's origin
    /// 3.25 below the ball) by `1e-12` of its radius meets it in a loop
    /// `1.8e-6` high, `|q|` within `3e-13` of one: its rise's jets over the
    /// whole piece are defined and enclose its points (`acos(q)` taken with
    /// `g`, `s - g` and `s + g` about the heights' middle: as written in `w`
    /// they spread by `2e-5` over the piece and left the jets undefined on
    /// all but pieces of `2^-26` of it, the integrals along it an hour).
    #[test]
    fn a_shallow_loops_rise_has_jets_over_the_whole_piece() {
        let cyl = frame([0.0, 0.0, -3.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let m = Rise {
            frame: cyl,
            radius: 1.0,
            half_angle: 0.0,
            centre: Point3::new(1.74999999999925, 0.0, 0.25),
            sphere_radius: 0.75,
            sign: -1.0,
            start: 3.2499990814005173,
            sweep: 1.83719896584833e-6,
        };
        let whole = Fast::exact_f64(0.0).union(&Fast::exact_f64(1.0));
        let jet = rise_jet(&m, &Jet::variable(whole, ORDER + 1)).expect("jets over the piece");
        for k in 0..=8 {
            let p = m.point(k as f64 / 8.0);
            for (i, j) in jet.iter().enumerate() {
                let (lo, hi) = j.c[0].bounds_f64();
                assert!(
                    lo <= [p.x, p.y, p.z][i] && [p.x, p.y, p.z][i] <= hi,
                    "{k} {i}"
                );
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
            assert!((p.y.hypot_rn(p.z) - 0.5).abs() < 1e-14, "{p:?}");
            let rho = p.x.hypot_rn(p.y);
            assert!(((rho - 2.5).hypot_rn(p.z) - 1.0).abs() < 1e-14, "{p:?}");
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

    /// Two tori on parallel axes, the fuzz target's in its tilted frame
    /// (S9d.4b.2b): to order 13 at a point the meeting's jets stay narrow
    /// (each step's linear forms in its provisional terms, each taken
    /// once; the plain recurrence's corrections counted them twice, a
    /// width growing twentyfold an order, `+-7e4` at the top), and they
    /// enclose the binary64 points over an interval base.
    #[test]
    fn toric_jets_stay_narrow_to_high_orders() {
        let m = Toric {
            frame: frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
            major: 1.3125,
            minor: 0.8203125,
            other: frame(
                [1.0416666666666667, -2.125, 0.9],
                [0.0, 3.0, 4.0],
                [1.0, 0.0, 0.0],
            ),
            other_radius: 1.5,
            other_sphere: false,
            other_half_angle: 0.0,
            other_minor: 0.5625,
            over_v: false,
            window: [-1.007081089716003, 1.9214080444695982],
            start: 0.9550889022507751,
            sweep: 1.4751893912374456,
        };
        let width = |x: &Fast| {
            let (lo, hi) = x.bounds_f64();
            hi - lo
        };
        for f in [0.25, 0.5, 0.75] {
            let ([_, v], point) = toric_jet(&m, &Jet::variable(Fast::exact_f64(f), 13)).unwrap();
            let p = m.point(f).to_array();
            for (i, j) in point.iter().enumerate() {
                let (lo, hi) = j.c[0].bounds_f64();
                assert!(lo - 1e-14 <= p[i] && p[i] <= hi + 1e-14, "{f} {i}");
            }
            if f == 0.5 {
                assert!(width(&v.c[13]) < 100.0 && width(&point[0].c[13]) < 100.0);
            }
        }
        let base = Fast::exact_f64(0.5).union(&Fast::exact_f64(0.5 + 1.0 / 128.0));
        let (_, jet) = toric_jet(&m, &Jet::variable(base, 13)).unwrap();
        for f in [0.5, 0.503, 0.5078125] {
            let p = m.point(f).to_array();
            for (i, j) in jet.iter().enumerate() {
                let (lo, hi) = j.c[0].bounds_f64();
                assert!(lo - 1e-14 <= p[i] && p[i] <= hi + 1e-14, "{f} {i}");
            }
        }
    }

    /// A rod of radius `2 + 1e-6` about a tilted torus's axis (radii 3 and
    /// 1; its frame built in the torus's, coaxial within rounding): its
    /// meeting's ring above `v = pi`, `1.4e-3` from the double root. The
    /// functionals' squares' jets, whose terms in `u` cancel only exactly,
    /// widened a hundredfold an order there (`16` at order 12 at a point,
    /// none over a piece of `1e-3`, some 65,000 pieces an integral); the
    /// harmonic forms' stay within `1e-7` to order 13 at points and over a
    /// sixty-fourth of the turn (`2e-8` at order 12: the frames' rounding in
    /// `cos 2u`'s coefficient, `1e-16`, over the slope `6e-3`, times `(4
    /// pi)^k / k!`), and enclose the binary64 points.
    #[test]
    fn a_meeting_near_a_double_root_coaxial_within_rounding_stays_narrow() {
        let torus = frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
        let base = torus.point(crate::Point2::new(0.0, 0.0), -1.5).to_array();
        let rod = frame(base, torus.normal().to_array(), torus.x().to_array());
        let m = Toric {
            frame: torus,
            major: 3.0,
            minor: 1.0,
            other: rod,
            other_radius: 2.0 + 1e-6,
            other_sphere: false,
            other_half_angle: 0.0,
            other_minor: 0.0,
            over_v: false,
            window: [std::f64::consts::PI + 1e-9, 6.0],
            start: 0.0,
            sweep: std::f64::consts::TAU,
        };
        let width = |x: &Fast| {
            let (lo, hi) = x.bounds_f64();
            hi - lo
        };
        let near = |jet: &[Jet<Fast>; 3], f: f64| {
            let p = m.point(f).to_array();
            for (i, j) in jet.iter().enumerate() {
                let (lo, hi) = j.c[0].bounds_f64();
                assert!(lo - 1e-14 <= p[i] && p[i] <= hi + 1e-14, "{f} {i}");
            }
        };
        for f in [0.25, 0.5, 0.75] {
            let ([_, v], point) = toric_jet(&m, &Jet::variable(Fast::exact_f64(f), 13)).unwrap();
            let (_, at) = m.angles(f);
            let (lo, hi) = v.c[0].bounds_f64();
            assert!(lo <= at && at <= hi, "{f} {lo} {hi} {at}");
            assert!((at - std::f64::consts::PI - 1.414e-3).abs() < 1e-5, "{at}");
            assert!(v.c.iter().all(|x| width(x) < 1e-7), "{f}");
            near(&point, f);
        }
        let base = Fast::exact_f64(0.5).union(&Fast::exact_f64(0.5 + 1.0 / 64.0));
        let ([_, v], point) = toric_jet(&m, &Jet::variable(base, 13)).unwrap();
        assert!(v.c.iter().all(|x| width(x) < 1e-7));
        for f in [0.5, 0.505, 0.5 + 1.0 / 64.0] {
            near(&point, f);
        }
    }

    /// A reversed use's variable `B - s` (a pcurve's `1 - f`): the memo's
    /// jets, the ones in `B + s` with every odd coefficient negated, are
    /// the evaluation's own bit for bit (zeros aside, whose signs no bound
    /// takes), over points and pieces, to order 13.
    #[test]
    fn toric_jets_of_a_reversed_variable_are_the_reflected_ones() {
        let m = Toric {
            frame: frame([1.0, -2.0, 0.5], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]),
            major: 1.3125,
            minor: 0.8203125,
            other: frame(
                [1.0416666666666667, -2.125, 0.9],
                [0.0, 3.0, 4.0],
                [1.0, 0.0, 0.0],
            ),
            other_radius: 1.5,
            other_sphere: false,
            other_half_angle: 0.0,
            other_minor: 0.5625,
            over_v: false,
            window: [-1.007081089716003, 1.9214080444695982],
            start: 0.9550889022507751,
            sweep: 1.4751893912374456,
        };
        let bits = |j: &Jet<Fast>| -> Vec<u64> {
            j.c.iter()
                .flat_map(|x| {
                    let (lo, hi) = x.bounds_f64();
                    [(lo + 0.0).to_bits(), (hi + 0.0).to_bits()]
                })
                .collect()
        };
        let bases = [
            Fast::exact_f64(0.375),
            Fast::exact_f64(0.5).union(&Fast::exact_f64(0.5 + 1.0 / 128.0)),
            Fast::exact_f64(0.75).union(&Fast::exact_f64(0.75 + 1.0 / 1024.0)),
        ];
        for base in bases {
            for order in [0, 1, 12, 13] {
                let reversed = Jet::variable(base, order)
                    .neg()
                    .add_constant(&Fast::exact_f64(1.0));
                let ([u, v], point) = toric_jet(&m, &reversed).unwrap();
                let ([du, dv], dpoint) = toric_jet_of(&m, &reversed).unwrap();
                assert_eq!(bits(&u), bits(&du), "{order}");
                assert_eq!(bits(&v), bits(&dv), "{order}");
                for i in 0..3 {
                    assert_eq!(bits(&point[i]), bits(&dpoint[i]), "{order} {i}");
                }
            }
        }
    }

    /// Two coaxial cones' meeting, a circle (S9d.3c's pair in one frame):
    /// over a piece its quadratic's coefficients are forms of degree two in
    /// the angle, constant here, and its discriminant's natural extension
    /// over the piece is replaced by its mean-value form, so the jets stay
    /// narrow to order 13 over a sixty-fourth of a turn (the plain
    /// products of the dot products' jets took thousands of pieces a turn);
    /// they enclose the binary64 points.
    #[test]
    fn coaxial_meetings_stay_narrow_over_a_piece() {
        let axis = frame([0.0, 0.0, 0.0], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
        let below = frame([0.0, -0.6, -0.8], [0.0, 3.0, 4.0], [1.0, 0.0, 0.0]);
        let width = |x: &Fast| {
            let (lo, hi) = x.bounds_f64();
            hi - lo
        };
        // The ring `1 - w / 2 = 0.8 (w + 1)` at `w = 2 / 13` (the other
        // root's branch on the other nappe).
        let m = Meet {
            frame: axis,
            radius: 1.0,
            half_angle: (-0.5f64).atan_rn(),
            other: below,
            other_radius: 0.0,
            other_sphere: false,
            other_half_angle: 0.8f64.atan_rn(),
            sign: -1.0,
            start: 0.0,
            sweep: std::f64::consts::TAU,
        };
        let base = Fast::exact_f64(0.3).union(&Fast::exact_f64(0.3 + 1.0 / 64.0));
        let (_, jet) = meet_jet(&m, &Jet::variable(base, 13)).unwrap();
        for f in [0.3, 0.31, 0.3 + 1.0 / 64.0] {
            let p = m.point(f);
            let height = p.y * 0.6 + p.z * 0.8;
            assert!((height - 2.0 / 13.0).abs() < 1e-14, "{p:?}");
            for (i, j) in jet.iter().enumerate() {
                let (lo, hi) = j.c[0].bounds_f64();
                let at = p.to_array()[i];
                assert!(lo - 1e-14 <= at && at <= hi + 1e-14, "{f} {i}");
            }
        }
        assert!(jet.iter().all(|j| width(&j.c[13]) < 1.0));
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
            let rho = p.x.hypot_rn(p.y);
            assert!(((rho - 2.5).hypot_rn(p.z) - 1.0).abs() < 1e-14, "{p:?}");
            let core = (p.y - 2.5).hypot_rn(p.z);
            assert!(((core - 1.25).hypot_rn(p.x) - 0.5).abs() < 1e-14, "{p:?}");
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

    /// The narrow tier's centred form (the shallow lens, HANDOFF item 5
    /// (d)): `F(u, v) = (u + 100)^2 - 200 u - 10^4 + v^2 cos u`, a sum of
    /// terms of `10^4` cancelling to `u^2 + v^2 cos u`, at jets whose
    /// coefficients are enclosed `2e-10` wide (a shallow meeting's point
    /// along its face): the natural form is some `4e-8` wide, the centred
    /// one below `1e-10` (the point's rounding), and every point jet inside
    /// the enclosures (their coefficients' ends and middles) has its value
    /// meeting the centred form's coefficient by coefficient, its exact
    /// value inside it.
    #[test]
    fn centred_forms_hold_every_jet_inside_and_narrow_a_cancelling_sum() {
        let f = |u: &Jet<Fast>, v: &Jet<Fast>| {
            let (co, _) = u.cos_sin();
            let big = u.add_constant(&Fast::exact_f64(100.0)).square();
            let value = big
                .sub(&u.scale(&Fast::exact_f64(200.0)))
                .add_constant(&Fast::exact_f64(-1e4))
                .add(&v.square().mul(&co));
            Some(vec![value])
        };
        let fu = |u: &Jet<Fast>, v: &Jet<Fast>| {
            let (_, si) = u.cos_sin();
            let d = u
                .add_constant(&Fast::exact_f64(100.0))
                .scale(&Fast::exact_f64(2.0))
                .add_constant(&Fast::exact_f64(-200.0))
                .sub(&v.square().mul(&si));
            Some(vec![d])
        };
        let fv = |u: &Jet<Fast>, v: &Jet<Fast>| {
            let (co, _) = u.cos_sin();
            Some(vec![v.scale(&Fast::exact_f64(2.0)).mul(&co)])
        };
        let wide = |centre: f64, slope: f64| {
            let w = 1e-10;
            let mut j = Jet::variable(Fast::exact_f64(centre), 4);
            j.c[0] = Fast::exact_f64(centre - w).union(&Fast::exact_f64(centre + w));
            j.c[1] = Fast::exact_f64(slope - w).union(&Fast::exact_f64(slope + w));
            j.c[2] = Fast::exact_f64(-w).union(&Fast::exact_f64(w));
            j
        };
        let (u, v) = (wide(1e-3, 0.5), wide(-2e-3, 0.25));
        let natural = f(&u, &v).unwrap().remove(0);
        let narrow = centred(&u, &v, &f, &fu, &fv).unwrap().remove(0);
        let width = |j: &Jet<Fast>| {
            let (lo, hi) = j.c[0].bounds_f64();
            hi - lo
        };
        assert!(width(&natural) > 1e-8, "{}", width(&natural));
        assert!(width(&narrow) < 1e-10, "{}", width(&narrow));
        for pick in 0..27 {
            let point = |j: &Jet<Fast>, k: usize| Jet {
                c: j.c
                    .iter()
                    .enumerate()
                    .map(|(i, x)| {
                        let (lo, hi) = x.bounds_f64();
                        let choice = if i < 3 {
                            (pick / [1, 3, 9][i] + k) % 3
                        } else {
                            1
                        };
                        Fast::exact_f64([lo, 0.5 * lo + 0.5 * hi, hi][choice])
                    })
                    .collect(),
            };
            let (pu, pv) = (point(&u, 0), point(&v, 1));
            let value = f(&pu, &pv).unwrap().remove(0);
            for (x, y) in value.c.iter().zip(&narrow.c) {
                let ((xl, xh), (yl, yh)) = (x.bounds_f64(), y.bounds_f64());
                assert!(
                    xl <= yh && yl <= xh,
                    "{pick}: [{xl}, {xh}] apart from [{yl}, {yh}]"
                );
            }
            // The value itself, `u^2 + v^2 cos u` without the cancelling
            // terms, at the point (to within its binary64 rounding).
            let (u0, v0) = (pu.c[0].bounds_f64().0, pv.c[0].bounds_f64().0);
            let exact = u0 * u0 + v0 * v0 * u0.cos_rn();
            let (yl, yh) = narrow.c[0].bounds_f64();
            assert!(
                yl - 1e-18 <= exact && exact <= yh + 1e-18,
                "{pick}: {exact}"
            );
        }
    }
}
