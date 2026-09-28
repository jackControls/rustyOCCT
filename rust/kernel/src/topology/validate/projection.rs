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
use crate::topology::{Curve3, Projection, Surface};
use crate::Frame3;

/// Integration widths and depths for the integrals along projections.
pub(super) const ORDER: usize = 14;
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
    let jet = conic_jet(curve, &Jet::variable(c::<T>(t), 0))?;
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
        Curve3::BSpline(_) => return None,
    })
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

/// A conic edge's point at a fraction in binary64 intervals (tessellation).
pub(crate) fn conic_point_fast(curve: &Curve3, t: f64) -> Option<[crate::certified::Fast; 3]> {
    conic_point(curve, t)
}
