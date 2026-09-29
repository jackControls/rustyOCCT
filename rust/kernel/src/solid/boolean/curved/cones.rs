//! S9d.3b.1: a cone against a cylinder, a sphere or another cone
//! (REVIEW_NOTES.md, S9d.3b refined). One input is the carrier: a cylinder
//! before a cone, whose ruling's leading coefficient `A` (a cylinder's a
//! constant, a cone's a quadratic form in its angle) has no real root.
//! Its discriminant `D = B^2 - A C` positive all round gives two rings over
//! its angle, negative all round the surfaces apart; changing sign on
//! every carrier is a loop (S9d.3b.2's). Two cones whose quadrics differ by
//! an affine function (their `A` zero for every ruling: parallel axes,
//! equal slopes) meet on that plane, S9d.3a's section of one.
use super::meet::{tangency, CylPair};
use super::model::*;
use super::num::*;
use super::procedural::{other_cone, other_of, other_sphere, MeetCrv, Other, Quartic, Ruled};
use super::turned::{near_node, negative_chart, roots, ruled_discriminant, Chart, Form};
use crate::solid::split::{q, rational_f64, zero};
use crate::{Error, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// A pair of a cone and a curved face meeting in loops: S9d.3b.2's.
fn loops_later() -> Error {
    Error::OutOfDomain("a cone meeting a curved face in a loop (S9d.3b.2)")
}

/// A face as a carrier: its ruled frame (a cylinder or a cone), or none.
fn ruled_of<'a>(m: &'a Prism, f: usize, axis: &'a P2, zero_k: &'a R) -> Option<Ruled<'a>> {
    match &m.faces[f].surf {
        Surf::Cyl { c, r, .. } => Some(Ruled {
            f: &m.f,
            c,
            r,
            k: zero_k,
        }),
        Surf::Cone { b, k } => Some(Ruled {
            f: &m.f,
            c: axis,
            r: b,
            k,
        }),
        _ => None,
    }
}

/// A face's quadric as the other of a carrier.
fn other_face(m: &Prism, f: usize) -> Other {
    match &m.faces[f].surf {
        Surf::Cyl { c, r, .. } => other_of(&m.f, c, r),
        Surf::Sphere { c, r } => other_sphere(c, r),
        Surf::Cone { b, k } => other_cone(&m.f, b, k),
        Surf::Plane { .. } | Surf::Torus => unreachable!("a quadric face"),
    }
}

/// Whether a form vanishes on the whole circle (`cos^2 + sin^2 = 1`
/// reduced: its chart's polynomial zero and its antipode's value zero).
fn vanishes(a: &Form) -> bool {
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    super::turned::trim(a.poly(&chart)).is_empty() && a.value(&[int(-1), zero()]) == zero()
}

/// Whether a form has no real zero on the circle (its chart's polynomial
/// no real root, its antipode nonzero).
fn nowhere_zero(a: &Form) -> bool {
    if let Some(k) = a.as_constant() {
        return k != zero();
    }
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    a.value(&[int(-1), zero()]) != zero() && matches!(roots(&a.poly(&chart)), Ok(r) if r.is_empty())
}

/// A lower bound of `|A|` over the turn, for the near-node test (sampled
/// at rational directions, halved: a guard, not a decision).
fn a_bound(a: &Form) -> R {
    if let Some(k) = a.as_constant() {
        return if k < zero() { -k } else { k };
    }
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let least = (-16..=16)
        .map(|i| {
            let v = rational_f64(&a.value(&chart.at(&(int(i) / int(4)))));
            v.abs()
        })
        .fold(f64::INFINITY, f64::min);
    q(least / 2.0)
}

/// How a cone's face and a curved face of the other input meet (`fa` of
/// operand 0, `fb` of operand 1).
pub(super) fn cone_pair(ms: [&Prism; 2], faces: [usize; 2], res: f64) -> Result<CylPair> {
    if (0..2).any(|k| matches!(ms[k].faces[faces[k]].surf, Surf::Torus)) {
        return Err(Error::OutOfDomain("a torus against a curved face (S9d.4b)"));
    }
    let axis: P2 = [zero(), zero()];
    let zero_k = zero();
    // Carriers: cylinders first, then cones.
    let mut order: Vec<usize> = (0..2)
        .filter(|&k| ruled_of(ms[k], faces[k], &axis, &zero_k).is_some())
        .collect();
    order.sort_by_key(|&k| matches!(ms[k].faces[faces[k]].surf, Surf::Cone { .. }));
    let mut parallel = None;
    for k in order {
        let ruled = ruled_of(ms[k], faces[k], &axis, &zero_k).expect("a ruled face");
        // The carrier's apex (a cone's, virtual for a frustum) on the other
        // quadric: every ruling meets it there.
        let other = other_face(ms[1 - k], faces[1 - k]);
        if *ruled.k != zero() {
            let w = -ruled.r / ruled.k;
            let apex = qv(&ruled.f.point(&zero(), &zero(), &w));
            if other.value(&apex).sign() == Ordering::Equal {
                return Err(Error::Degenerate(
                    "a cone's apex on the other input's surface",
                ));
            }
        }
        let (a, d) = ruled_discriminant(ruled, &other);
        if vanishes(&a) {
            // Every ruling asymptotic to the other: the quadrics differ by
            // an affine function (two cones with parallel axes and equal
            // slopes).
            parallel = Some(k);
            continue;
        }
        if !nowhere_zero(&a) {
            continue;
        }
        // A ruling within rounding of the other's asymptotic direction (a
        // cylinder along a cone's ruling): one branch past any binary64
        // edge.
        let scale = a.max_coefficient();
        if rational_f64(&a_bound(&a)) <= 1e-12 * rational_f64(&scale).max(1.0) {
            return Err(Error::Degenerate(
                "a ruling within rounding of the other's direction",
            ));
        }
        if vanishes(&d) {
            // Tangent all along a curve (a sphere inscribed in a cone).
            return Err(tangency());
        }
        let chart = negative_chart(&d)?;
        let test = chart.clone().unwrap_or(Chart {
            c0: int(1),
            s0: zero(),
        });
        // Two branches, or the surfaces, within the resolution of touching.
        near_node(&a_bound(&a), &d, &test, res)?;
        match chart {
            None => {
                let piece = |plus| MeetCrv::ruled(k, ruled, &other, plus, None);
                return Ok(CylPair::Quartic(Box::new(Quartic {
                    pieces: vec![piece(true), piece(false)],
                    switches: Vec::new(),
                })));
            }
            Some(chart) => {
                if roots(&d.poly(&chart))?.is_empty() {
                    // No ruling meets the other quadric.
                    return Ok(CylPair::Apart);
                }
                // A cone and a sphere in loops (S9d.3b.2): graphs over the
                // cone's height about its rulings' tangencies, over its
                // angle about its circles', as S9d.2b's.
                if let Surf::Sphere { c, r } = &ms[1 - k].faces[faces[1 - k]].surf {
                    return super::spheres::loops(k, ruled, c, r, &other, &d, &chart);
                }
            }
        }
    }
    if let Some(k) = parallel {
        return plane_of(ms, faces, k);
    }
    // A carrier whose rulings reach the other's asymptotic directions (its
    // `A`'s simple real roots) and meet it twice elsewhere (S9d.3b.2): its
    // turn split at those directions, where one branch runs to infinity
    // (past the solids' ends) and the other through a finite switch.
    for k in 0..2 {
        let Some(ruled) = ruled_of(ms[k], faces[k], &axis, &zero_k) else {
            continue;
        };
        let other = other_face(ms[1 - k], faces[1 - k]);
        let (a, b, cc) = super::turned::ruled_quadratic(ruled, &other);
        let d = b.mul(&b).sub(&a.mul(&cc));
        if vanishes(&a) || vanishes(&d) {
            continue;
        }
        let Ok(dirs) = circle_roots(&a) else {
            continue;
        };
        if dirs.is_empty() || negative_chart(&d)?.is_some() {
            continue;
        }
        near_node(
            &a_bound(&a),
            &d,
            &Chart {
                c0: int(1),
                s0: zero(),
            },
            res,
        )?;
        let mut pieces = Vec::new();
        let mut switches = Vec::new();
        let n = dirs.len();
        for i in 0..n {
            let range = [dirs[i].clone(), dirs[(i + 1) % n].clone()];
            for plus in [true, false] {
                pieces.push(MeetCrv::ruled(k, ruled, &other, plus, Some(range.clone())));
            }
            // The finite branch's point there: `w = -C / 2B`.
            let cs = &dirs[i];
            let bv = b.value_q(cs);
            if bv.sign() == Ordering::Equal {
                return Err(tangency());
            }
            let w = b
                .value_q(cs)
                .scale(&int(2))
                .recip()
                .ok_or(tangency())?
                .mul(&cc.value_q(cs).neg());
            switches.push(ruling_point(ruled, cs, &w));
        }
        return Ok(CylPair::Quartic(Box::new(Quartic { pieces, switches })));
    }
    Err(loops_later())
}

/// A form's zeros on the circle, counter-clockwise from the angle `-pi`:
/// the directions of its chart's real roots (algebraic, `Q(alpha)`), the
/// antipode `(-1, 0)` last where it vanishes. A repeated root is refused.
fn circle_roots(a: &Form) -> Result<Vec<[Qd; 2]>> {
    use super::num::{Gen, K};
    use std::sync::Arc;
    let chart = Chart {
        c0: int(1),
        s0: zero(),
    };
    let p = super::turned::trim(a.poly(&chart));
    let mut out = Vec::new();
    for root in roots(&p)? {
        let g = Arc::new(Gen {
            poly: p.clone(),
            root,
        });
        let t = K::generator(&g);
        let den = t.mul(&t).add(&K::Rat(int(1)));
        let inv = den
            .recip()
            .ok_or(Error::ComputationLimit("a direction at infinity"))?;
        let cos = K::Rat(int(1)).sub(&t.mul(&t)).mul(&inv);
        let sin = t.scale(&int(2)).mul(&inv);
        out.push([Qd::of(cos), Qd::of(sin)]);
    }
    if a.value(&[int(-1), zero()]) == zero() {
        out.push([Qd::rat(int(-1)), Qd::rat(zero())]);
    }
    Ok(out)
}

/// A ruling's point at `w` over a direction of any field.
fn ruling_point(ruled: Ruled, cs: &[Qd; 2], w: &Qd) -> QV {
    let f = ruled.f;
    let base = qadd(
        &qv(&f.point(&ruled.c[0], &ruled.c[1], &zero())),
        &qadd(
            &qscale(&f.x, &cs[0].scale(ruled.r)),
            &qscale(&f.y, &cs[1].scale(ruled.r)),
        ),
    );
    let dir = qadd(
        &qv(&f.n),
        &qadd(
            &qscale(&f.x, &cs[0].scale(ruled.k)),
            &qscale(&f.y, &cs[1].scale(ruled.k)),
        ),
    );
    qadd(&base, &dir.map(|x| x.mul(w)))
}

/// A quadric's quadratic, linear and constant parts, `p^T M p + l . p + c`.
fn parts(o: &Other) -> ([[R; 3]; 3], V, R) {
    let mut m: [[R; 3]; 3] = std::array::from_fn(|_| std::array::from_fn(|_| zero()));
    let mut l: V = [zero(), zero(), zero()];
    let mut c = zero();
    for (g, e) in o.g.iter().zip(&o.e) {
        for (i, row) in m.iter_mut().enumerate() {
            for (j, x) in row.iter_mut().enumerate() {
                *x += &g[i] * &g[j];
            }
            l[i] -= int(2) * e * &g[i];
        }
        c += e * e;
    }
    // Less (t h . p + (r - t e_h))^2.
    let r0 = &o.r - &o.t * &o.eh;
    for (i, row) in m.iter_mut().enumerate() {
        for (j, x) in row.iter_mut().enumerate() {
            *x -= &o.t * &o.t * &o.h[i] * &o.h[j];
        }
        l[i] -= int(2) * &o.t * &r0 * &o.h[i];
    }
    c -= &r0 * &r0;
    (m, l, c)
}

/// Two cones whose quadrics differ by an affine function: their meeting
/// lies on that plane, the section of the carrier `k`'s cone.
fn plane_of(ms: [&Prism; 2], faces: [usize; 2], k: usize) -> Result<CylPair> {
    let (m1, l1, c1) = parts(&other_face(ms[k], faces[k]));
    let (m2, l2, c2) = parts(&other_face(ms[1 - k], faces[1 - k]));
    // M2 = lambda M1.
    let (i, j) = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .find(|&(i, j)| m1[i][j] != zero())
        .ok_or(Error::Degenerate("a cone's quadric"))?;
    let lambda = &m2[i][j] / &m1[i][j];
    if (0..3).any(|i| (0..3).any(|j| m2[i][j] != &lambda * &m1[i][j])) {
        return Err(loops_later());
    }
    let m: V = [0, 1, 2].map(|i| &l1[i] * &lambda - &l2[i]);
    let c = &c1 * &lambda - &c2;
    if is_zero(&m) {
        return if c == zero() {
            Err(Error::Degenerate("two cones on one surface"))
        } else {
            Ok(CylPair::Apart)
        };
    }
    // m . p + c = 0 through p0 = -c m / |m|^2.
    let p0 = scale(&m, &(-&c / dot(&m, &m)));
    Ok(CylPair::Plane { k, p: p0, m })
}
