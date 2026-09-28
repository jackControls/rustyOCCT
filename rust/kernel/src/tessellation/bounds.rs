//! Certified quantities of a tessellation, in the outward-rounded binary64
//! tier of `certified.rs` (trigonometry included): points of surfaces and
//! curves with their rounding, gaps between a node and the surface at its
//! parameter point, and the deflection and normal turn of each triangle and
//! segment (`MATHEMATICS.md`, tessellation).
//!
//! For a parameter triangle with extents `U`, `V` and `a >= |S_uu|`,
//! `b >= |S_uv|`, `c >= |S_vv|` over its box, linear interpolation deviates
//! from the surface by at most `(a U^2 + 2 b U V + c V^2) / 8`; the normal
//! turns by at most `n_u U + n_v V`. The coefficients are closed forms of
//! each analytic surface times `sigma`, a bound on the stored frame's
//! spectral norm (Gershgorin on its Gram matrix).
use crate::certified::{Fast, Real};
use crate::topology::{Curve3, Surface};
use crate::{Frame3, Point2, Point3, Vec3};

type P = [Fast; 3];

fn c(x: f64) -> Fast {
    Fast::exact_f64(x)
}

/// An upper bound of the enclosed value.
pub(super) fn upper(x: &Fast) -> f64 {
    x.bounds_f64().1
}

/// An upper bound of the enclosed value's magnitude.
fn magnitude(x: &Fast) -> f64 {
    let (lo, hi) = x.bounds_f64();
    (-lo).max(hi)
}

fn vector(v: Vec3) -> P {
    [c(v.x), c(v.y), c(v.z)]
}

fn add(a: &P, b: &P) -> P {
    [a[0].add(&b[0]), a[1].add(&b[1]), a[2].add(&b[2])]
}

fn scale(a: &P, k: &Fast) -> P {
    [a[0].mul(k), a[1].mul(k), a[2].mul(k)]
}

fn dot(a: &P, b: &P) -> Fast {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

/// A bound on the spectral norm of the frame's axes as columns: any
/// combination `p x + q y + r n` is at most this times `|(p, q, r)|`.
pub(super) fn frame_norm(frame: &Frame3) -> f64 {
    let axes = [vector(frame.x()), vector(frame.y()), vector(frame.normal())];
    let mut worst = 0.0f64;
    for i in 0..3 {
        let mut row = dot(&axes[i], &axes[i]);
        for (j, axis) in axes.iter().enumerate() {
            if j != i {
                row = row.add(&c(magnitude(&dot(&axes[i], axis))));
            }
        }
        worst = worst.max(upper(&row));
    }
    upper(&c(worst).sqrt())
}

/// `origin + p x + q y + r n` in the tier.
fn combine(frame: &Frame3, p: &Fast, q: &Fast, r: &Fast) -> P {
    let o = frame.origin();
    let origin = [c(o.x), c(o.y), c(o.z)];
    let out = add(&origin, &scale(&vector(frame.x()), p));
    let out = add(&out, &scale(&vector(frame.y()), q));
    add(&out, &scale(&vector(frame.normal()), r))
}

/// The exact point of a plane, cylinder, cone, sphere or torus at a
/// parameter point, enclosed.
pub(super) fn surface_point(s: &Surface, uv: Point2) -> Option<P> {
    let zero = c(0.0);
    Some(match s {
        Surface::Plane(frame) => combine(frame, &c(uv.x), &c(uv.y), &zero),
        Surface::Cylinder { frame, radius } => {
            let (cu, su) = Fast::cos_sin(&c(uv.x));
            let r = c(*radius);
            combine(frame, &r.mul(&cu), &r.mul(&su), &c(uv.y))
        }
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => {
            let (cu, su) = Fast::cos_sin(&c(uv.x));
            let (ca, sa) = Fast::cos_sin(&c(*half_angle));
            let rho = c(*radius).add(&c(uv.y).mul(&sa));
            combine(frame, &rho.mul(&cu), &rho.mul(&su), &c(uv.y).mul(&ca))
        }
        Surface::Sphere { frame, radius } => {
            let (cu, su) = Fast::cos_sin(&c(uv.x));
            let (cv, sv) = Fast::cos_sin(&c(uv.y));
            let rho = c(*radius).mul(&cv);
            combine(frame, &rho.mul(&cu), &rho.mul(&su), &c(*radius).mul(&sv))
        }
        Surface::Torus {
            frame,
            major,
            minor,
        } => {
            let (cu, su) = Fast::cos_sin(&c(uv.x));
            let (cv, sv) = Fast::cos_sin(&c(uv.y));
            let rho = c(*major).add(&c(*minor).mul(&cv));
            combine(frame, &rho.mul(&cu), &rho.mul(&su), &c(*minor).mul(&sv))
        }
        Surface::BSpline(_) => return None,
    })
}

/// The exact point of a line, circle or arc at an edge fraction, enclosed.
pub(super) fn curve_point(curve: &Curve3, t: f64) -> Option<P> {
    Some(match curve {
        Curve3::LineSegment { start, end } => {
            let t = c(t);
            let s = [c(start.x), c(start.y), c(start.z)];
            let e = [c(end.x), c(end.y), c(end.z)];
            let d = [e[0].sub(&s[0]), e[1].sub(&s[1]), e[2].sub(&s[2])];
            add(&s, &scale(&d, &t))
        }
        Curve3::Circle { frame, radius } => {
            let (ca, sa) = Fast::cos_sin(&Fast::two_pi().mul(&c(t)));
            let r = c(*radius);
            combine(frame, &r.mul(&ca), &r.mul(&sa), &c(0.0))
        }
        Curve3::CircularArc {
            frame,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let a = c(*start_angle).add(&c(*sweep_angle).mul(&c(t)));
            let (ca, sa) = Fast::cos_sin(&a);
            let r = c(*radius);
            combine(frame, &r.mul(&ca), &r.mul(&sa), &c(0.0))
        }
        Curve3::BSpline(_) => return None,
    })
}

/// The binary64 point nearest the enclosure's middle, and an upper bound of
/// its distance from the enclosed point.
pub(super) fn settle(p: &P) -> (Point3, f64) {
    let mid = |x: &Fast| {
        let (lo, hi) = x.bounds_f64();
        0.5 * lo + 0.5 * hi
    };
    let point = Point3::new(mid(&p[0]), mid(&p[1]), mid(&p[2]));
    (point, gap(point, p))
}

/// An upper bound of the distance from `x` to the enclosed point.
pub(super) fn gap(x: Point3, p: &P) -> f64 {
    let d = [c(x.x).sub(&p[0]), c(x.y).sub(&p[1]), c(x.z).sub(&p[2])];
    let squared = d[0].square().add(&d[1].square()).add(&d[2].square());
    if upper(&squared) == 0.0 {
        return 0.0;
    }
    let g = upper(&squared.sqrt());
    if g.is_finite() {
        g
    } else {
        f64::INFINITY
    }
}

/// Upper bounds `[a, b, c, n_u, n_v]` over a parameter box with `v` in
/// `[v0, v1]`: `a >= |S_uu|`, `b >= |S_uv|`, `c >= |S_vv|`, `n_u >= |N_u|`,
/// `n_v >= |N_v|` (the normals before the frame allowance).
pub(super) fn coefficients(s: &Surface, v0: f64, v1: f64) -> [f64; 5] {
    let vm = 0.5 * v0 + 0.5 * v1;
    let half = upper(&c(v1).sub(&c(vm))).max(upper(&c(vm).sub(&c(v0))));
    // |cos| and |sin| are 1-Lipschitz: bounds over [v0, v1] from its middle.
    let (cv, sv) = Fast::cos_sin(&c(vm));
    let sup_cos = upper(&c(magnitude(&cv)).add(&c(half))).min(1.0);
    let sup_sin = upper(&c(magnitude(&sv)).add(&c(half))).min(1.0);
    let times = |k: f64, x: f64| upper(&c(k).mul(&c(x)));
    match s {
        Surface::Plane(_) | Surface::BSpline(_) => [0.0; 5],
        Surface::Cylinder { frame, radius } => {
            [times(frame_norm(frame), radius.abs()), 0.0, 0.0, 1.0, 0.0]
        }
        Surface::Cone {
            frame,
            radius,
            half_angle,
        } => {
            let (ca, sa) = Fast::cos_sin(&c(*half_angle));
            let rho = |v: f64| magnitude(&c(*radius).add(&c(v).mul(&sa)));
            let sigma = frame_norm(frame);
            [
                times(sigma, rho(v0).max(rho(v1))),
                times(sigma, magnitude(&sa)),
                0.0,
                magnitude(&ca),
                0.0,
            ]
        }
        Surface::Sphere { frame, radius } => {
            let sigma = frame_norm(frame);
            let r = times(sigma, radius.abs());
            [times(r, sup_cos), times(r, sup_sin), r, sup_cos, 1.0]
        }
        Surface::Torus {
            frame,
            major,
            minor,
        } => {
            let sigma = frame_norm(frame);
            let rho = upper(&c(major.abs()).add(&c(minor.abs()).mul(&c(sup_cos))));
            let r = times(sigma, minor.abs());
            [times(sigma, rho), times(r, sup_sin), r, sup_cos, 1.0]
        }
    }
}

/// Relative allowance on normal turns for the stored frames' departure from
/// orthonormality (at most a few units in the last place).
const FRAME_ALLOWANCE: f64 = 1e-12;

/// Lower bounds `(|S_u|, |S_v|)` over `v` in `[v0, v1]`, before the frame
/// allowance; zero where `S_u` may vanish (a pole in the box).
fn speeds(s: &Surface, v0: f64, v1: f64) -> (f64, f64) {
    let lower = |x: &Fast| x.bounds_f64().0;
    match s {
        Surface::Cylinder { radius, .. } => (radius.abs(), 1.0),
        Surface::Cone {
            radius, half_angle, ..
        } => {
            let (_, sa) = Fast::cos_sin(&c(*half_angle));
            let rho = |v: f64| c(*radius).add(&c(v).mul(&sa)).bounds_f64();
            let ((a0, a1), (b0, b1)) = (rho(v0), rho(v1));
            // rho is affine in v: no sign change between the ends.
            let g = if a0 > 0.0 && b0 > 0.0 {
                a0.min(b0)
            } else if a1 < 0.0 && b1 < 0.0 {
                (-a1).min(-b1)
            } else {
                0.0
            };
            (g, 1.0)
        }
        Surface::Sphere { radius, .. } => {
            // cos is concave on [-pi/2, pi/2]: its least value at an end.
            if !(v0 > -std::f64::consts::FRAC_PI_2 && v1 < std::f64::consts::FRAC_PI_2) {
                return (0.0, radius.abs());
            }
            let least = lower(&Fast::cos_sin(&c(v0)).0).min(lower(&Fast::cos_sin(&c(v1)).0));
            (
                lower(&c(radius.abs()).mul(&c(least.max(0.0)))),
                radius.abs(),
            )
        }
        Surface::Torus { major, minor, .. } => {
            let vm = 0.5 * v0 + 0.5 * v1;
            let half = upper(&c(v1).sub(&c(vm))).max(upper(&c(vm).sub(&c(v0))));
            let least = lower(&Fast::cos_sin(&c(vm)).0.sub(&c(half))).max(-1.0);
            let g = lower(&c(major.abs()).add(&c(minor.abs()).mul(&c(least))));
            (g.max(0.0), minor.abs())
        }
        Surface::Plane(_) | Surface::BSpline(_) => (0.0, 0.0),
    }
}

/// Relative allowance for the stored frames' departure from orthogonality
/// in the tangential correction (their Gram defect is a few units in the
/// last place).
const CORRECTION_ALLOWANCE: f64 = 1e-9;

/// The certified deflection (without node gaps) and normal turn of a
/// parameter triangle `p` (`pole`: its vertex at a pole, where the normal
/// is undefined and every `u` maps to the same point).
///
/// Two maps are certified and the smaller bound is kept. The linear one,
/// `Σ λ_i X_i ↦ S(Σ λ_i p_i)`, deviates by `(a U^2 + 2 b U V + c V^2) / 8`.
/// On every analytic surface `S_uv` is parallel to `S_u`, so it has no
/// normal component: the normal part of that deviation is at most
/// `(a U^2 + 2 t b U V + c V^2) / 8` with `t` the triangle's normal turn,
/// and the tangential part `e_t` is taken up by moving the parameter point by
/// `(|e_t| / |S_u|, |e_t| / |S_v|)`, which leaves at most the second-order
/// remainder `(a' du^2 + 2 b' du dv + c' dv^2) / 2` over the moved box.
pub(super) fn triangle_bound(s: &Surface, p: [Point2; 3], pole: Option<usize>) -> (f64, f64) {
    let range = |f: fn(&Point2) -> f64| {
        let (lo, hi) = p
            .iter()
            .map(f)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
                (l.min(x), h.max(x))
            });
        (lo, hi, upper(&c(hi).sub(&c(lo))))
    };
    let (_, _, u) = range(|q| q.x);
    let (v0, v1, v) = range(|q| q.y);
    let [a, b, cc, nu, nv] = coefficients(s, v0, v1);
    let (fu, fv) = (c(u), c(v));
    let quadratic = |a: f64, b: f64, cc: f64, x: &Fast, y: &Fast| {
        c(a).mul(&x.square())
            .add(&c(2.0 * b).mul(&x.mul(y)))
            .add(&c(cc).mul(&y.square()))
    };
    let linear = quadratic(a, b, cc, &fu, &fv).mul(&c(0.125));
    let turn = upper(
        &c(nu)
            .mul(&fu)
            .add(&c(nv).mul(&fv))
            .mul(&c(1.0 + FRAME_ALLOWANCE)),
    );
    let linear_bound = upper(&linear);
    if let Some(k) = pole {
        return (linear_bound.min(fan_bound(s, p, k)), turn);
    }
    if b == 0.0 || !linear_bound.is_finite() {
        return (linear_bound, turn);
    }
    // |N(p) - N(xi)| <= min(turn, 2) for points of the triangle.
    let t = turn.min(2.0);
    let normal = quadratic(a, upper(&c(t).mul(&c(b))), cc, &fu, &fv).mul(&c(0.125));
    let (g, sv) = speeds(s, v0, v1);
    if !(g > 0.0 && sv > 0.0) {
        return (linear_bound, turn);
    }
    let slack = c(1.0 + CORRECTION_ALLOWANCE);
    let du = linear.div(&c(g)).map(|x| x.mul(&slack));
    let dv = linear.div(&c(sv)).map(|x| x.mul(&slack));
    let (Some(du), Some(dv)) = (du, dv) else {
        return (linear_bound, turn);
    };
    let reach = upper(&dv);
    let [a2, b2, c2, _, _] = coefficients(s, v0 - reach, v1 + reach);
    let remainder = quadratic(a2, b2, c2, &du, &dv).mul(&c(0.5));
    let corrected = upper(&normal.add(&remainder));
    (linear_bound.min(corrected), turn)
}

/// The fan map of a triangle with vertex `k` at a pole: a point with
/// weights `λ` goes to `S(u_μ, Σ λ_i v_i)`, `u_μ` interpolating the other
/// two vertices' `u` in the ratio of their weights (continuous, since every
/// `u` maps to the pole, and linear on the opposite edge). It deviates by at
/// most the opposite edge's linear interpolation, `(a U^2 + 2 b U V' + c
/// V'^2) / 8`, plus the chord of a meridian, `c V^2 / 8`, with `V'` the
/// opposite edge's `v` extent and `V` the triangle's.
fn fan_bound(s: &Surface, p: [Point2; 3], k: usize) -> f64 {
    let (m, n) = ((k + 1) % 3, (k + 2) % 3);
    let extent = |x: f64, y: f64| upper(&c(x).sub(&c(y))).max(upper(&c(y).sub(&c(x))));
    let u = c(extent(p[m].x, p[n].x));
    let edge_v = c(extent(p[m].y, p[n].y));
    let v = c(extent(p[m].y, p[k].y).max(extent(p[n].y, p[k].y)));
    let (v0, v1) = p
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), q| {
            (a.min(q.y), b.max(q.y))
        });
    let [a, b, cc, _, _] = coefficients(s, v0, v1);
    let edge = c(a)
        .mul(&u.square())
        .add(&c(2.0 * b).mul(&u.mul(&edge_v)))
        .add(&c(cc).mul(&edge_v.square()));
    upper(&edge.add(&c(cc).mul(&v.square())).mul(&c(0.125)))
}

/// The certified deflection (without node gaps) and tangent turn of an edge
/// segment over a fraction step `dt`.
pub(super) fn segment_bound(curve: &Curve3, dt: f64) -> (f64, f64) {
    let (frame, radius, sweep) = match curve {
        Curve3::Circle { frame, radius } => (frame, *radius, Fast::two_pi()),
        Curve3::CircularArc {
            frame,
            radius,
            sweep_angle,
            ..
        } => (frame, *radius, c(sweep_angle.abs())),
        _ => return (0.0, 0.0),
    };
    let angle = sweep.mul(&c(dt));
    let deviation = c(frame_norm(frame))
        .mul(&c(radius.abs()))
        .mul(&angle.square())
        .mul(&c(0.125));
    (upper(&deviation), upper(&angle))
}

/// `a + b` rounded up.
pub(super) fn sum_up(a: f64, b: f64) -> f64 {
    upper(&c(a).add(&c(b)))
}
