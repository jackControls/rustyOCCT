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
use crate::certified::Real;
use crate::jet::{integrate_many, Jet};
use crate::topology::{Curve3, Meet, Projection, Rise, Spiric, Surface};
use crate::Frame3;

/// Integration widths and depths for the integrals along projections.
pub(super) const ORDER: usize = 12;
pub(super) const WIDTH: f64 = 1e-12;
pub(super) const DEPTH: usize = 40;

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
        Curve3::BSpline(_) => return None,
    })
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
    let f = if p.reversed {
        fraction.neg().add_constant(&c(1.0))
    } else {
        fraction.clone()
    };
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
    let point = curve_jet(&p.curve, &f)?;
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
    let integrand = |f: &Jet<T>| {
        // One order more, so the derivatives keep the order asked for.
        let longer = Jet::variable(f.c[0].clone(), f.order() + 1);
        let [u, v] = projection_jet(p, &longer)?;
        let (du, dv) = (u.derivative(), v.derivative());
        let cut = |j: &Jet<T>| Jet {
            c: j.c[..=f.order()].to_vec(),
        };
        g(&cut(&u), &cut(&v), &cut(&du), &cut(&dv))
    };
    integrate_many(&integrand, n, 0.0, 1.0, ORDER, WIDTH, DEPTH, relative)
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
    for k in 0..pieces {
        let (a, b) = (k as f64 / pieces as f64, (k + 1) as f64 / pieces as f64);
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
}
