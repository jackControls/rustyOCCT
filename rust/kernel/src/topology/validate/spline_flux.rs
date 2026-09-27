//! The orientation flux `∫∫ S·(S_u × S_v) du dv` of a face on a nonrational
//! spline surface (S4d of REVIEW_NOTES.md). On each Bézier patch, in local
//! coordinates `(ū, v̄)` on the unit square, the integrand `X·(X_ū × X_v̄)`
//! (poles divided by the common weight) is a tensor Bernstein polynomial;
//! its antiderivative `H(ū, v̄)` from `v̄ = 0`, plus the full columns
//! below the patch, `G = Σ H(ū, 1) + H(ū, v̄)`, is the global
//! antiderivative in `v` in local units. By Green's theorem the flux is
//! `-∮ G dū` along the face's loops, each piece in one patch: a line pcurve
//! is split exactly where it crosses a patch boundary, a spline pcurve's
//! Bézier piece must lie in one patch by its control points, and a chord's
//! ends must certainly share a patch. A rational surface, a periodic one or
//! a piece across a patch boundary is not decided.
use super::bernstein::{
    c, derivative, difference, lift, pcurve_arcs, piece_of, power, product, quotient_integral, r,
    ratio, scaled, sum, Bern,
};
use super::Lp;
use crate::certified::Real;
use crate::surface::ExactBezierSurface3;
use crate::topology::Curve2;
use crate::BSplineSurface3;
use num_rational::BigRational as R;
use std::cmp::Ordering;

/// A tensor Bernstein polynomial, `[i][j]` in `u` then `v`.
type Tensor<T> = Vec<Vec<T>>;

fn tensor_product<T: Real>(a: &Tensor<T>, b: &Tensor<T>) -> Tensor<T> {
    // Multiply along v for every pair of u rows, then combine along u with
    // the univariate product weights.
    let (m, n) = (a.len() - 1, b.len() - 1);
    let degree_v = a[0].len() + b[0].len() - 2;
    let mut out: Tensor<T> = vec![vec![T::exact_f64(0.0); degree_v + 1]; m + n + 1];
    for (i, row_a) in a.iter().enumerate() {
        for (k, row_b) in b.iter().enumerate() {
            let along_v = product(row_a, row_b);
            // The weight of B_i^m B_k^n in B_{i+k}^{m+n}.
            let mut indicator_a = vec![T::exact_f64(0.0); m + 1];
            indicator_a[i] = T::exact_f64(1.0);
            let mut indicator_b = vec![T::exact_f64(0.0); n + 1];
            indicator_b[k] = T::exact_f64(1.0);
            let weight = product(&indicator_a, &indicator_b)[i + k].clone();
            for (j, x) in along_v.iter().enumerate() {
                out[i + k][j] = out[i + k][j].add(&x.mul(&weight));
            }
        }
    }
    out
}

fn tensor_sum<T: Real>(a: &Tensor<T>, b: &Tensor<T>) -> Tensor<T> {
    // Equal shapes only (the three terms of a triple product have them).
    a.iter()
        .zip(b)
        .map(|(x, y)| x.iter().zip(y).map(|(p, q)| p.add(q)).collect())
        .collect()
}

/// The partial derivatives of a tensor along u and along v.
fn partial_u<T: Real>(a: &Tensor<T>) -> Tensor<T> {
    let m = (a.len() - 1) as f64;
    a.windows(2)
        .map(|w| {
            w[1].iter()
                .zip(&w[0])
                .map(|(p, q)| p.sub(q).mul(&T::exact_f64(m)))
                .collect()
        })
        .collect()
}

fn partial_v<T: Real>(a: &Tensor<T>) -> Tensor<T> {
    a.iter().map(derivative).collect()
}

/// The patch's flux antiderivative `H(ū, v̄) = ∫_0^v̄ X·(X_ū × X_v̄) ds`.
fn antiderivative<T: Real>(patch: &ExactBezierSurface3) -> Tensor<T> {
    let [du, dv] = patch.degrees();
    let poles = patch.homogeneous_poles();
    let x: [Tensor<T>; 3] = std::array::from_fn(|k| {
        (0..=du)
            .map(|i| {
                (0..=dv)
                    .map(|j| {
                        let p = &poles[i * (dv + 1) + j];
                        c(&(&p[k] / &p[3]))
                    })
                    .collect()
            })
            .collect()
    });
    let xu = x.clone().map(|t| partial_u(&t));
    let xv = x.clone().map(|t| partial_v(&t));
    let neg = |t: Tensor<T>| -> Tensor<T> {
        t.into_iter()
            .map(|row| row.into_iter().map(|q| q.neg()).collect())
            .collect()
    };
    let cross: [Tensor<T>; 3] = std::array::from_fn(|k| {
        let (a, b) = ((k + 1) % 3, (k + 2) % 3);
        tensor_sum(
            &tensor_product(&xu[a], &xv[b]),
            &neg(tensor_product(&xu[b], &xv[a])),
        )
    });
    let f = (1..3).fold(tensor_product(&x[0], &cross[0]), |acc, k| {
        tensor_sum(&acc, &tensor_product(&x[k], &cross[k]))
    });
    // ∫_0^v̄ Σ c_j B_j^n = Σ_k (Σ_{j<k} c_j)/(n+1) B_k^{n+1}.
    f.iter()
        .map(|row| {
            let scale = c::<T>(&ratio(1, row.len() as i64));
            let mut acc = T::exact_f64(0.0);
            let mut out = vec![T::exact_f64(0.0)];
            for x in row {
                acc = acc.add(x);
                out.push(acc.mul(&scale));
            }
            out
        })
        .collect()
}

/// `H(ū, 1)`: the last coefficient of every row.
fn column<T: Real>(h: &Tensor<T>) -> Bern<T> {
    h.iter().map(|row| row[row.len() - 1].clone()).collect()
}

/// `K(ū(τ), v̄(τ)) P_w^(m + n)` for homogeneous local coordinates
/// `(Ū, V̄, P_w)` with `ū = Ū/P_w`.
fn composed<T: Real>(k: &Tensor<T>, big_u: &Bern<T>, big_v: &Bern<T>, pw: &Bern<T>) -> Bern<T> {
    let (m, n) = (k.len() - 1, k[0].len() - 1);
    let (rest_u, rest_v) = (difference(pw, big_u), difference(pw, big_v));
    let binomial = |n: usize, i: usize| {
        (0..i).fold(ratio(1, 1), |acc, j| {
            acc * ratio((n - j) as i64, (j + 1) as i64)
        })
    };
    let bu: Vec<Bern<T>> = (0..=m)
        .map(|i| {
            scaled(
                &product(&power(big_u, i), &power(&rest_u, m - i)),
                &c(&binomial(m, i)),
            )
        })
        .collect();
    let bv: Vec<Bern<T>> = (0..=n)
        .map(|j| {
            scaled(
                &product(&power(big_v, j), &power(&rest_v, n - j)),
                &c(&binomial(n, j)),
            )
        })
        .collect();
    let mut out = vec![T::exact_f64(0.0)];
    for i in 0..=m {
        let row = (0..=n).fold(vec![T::exact_f64(0.0)], |acc, j| {
            sum(&acc, &scaled(&bv[j], &k[i][j]))
        });
        out = sum(&out, &product(&bu[i], &row));
    }
    out
}

struct Patch<T> {
    domain: [[R; 2]; 2],
    /// `H` plus the full columns below, as `G` in local units.
    g: Tensor<T>,
}

/// One homogeneous `(U, V, W)` piece in global `(u, v)`, as exact
/// Bernstein coordinates, with its exact hull for locating its patch.
fn locate<'a, T: Real>(patches: &'a [Patch<T>], points: &[(R, R)]) -> Option<&'a Patch<T>> {
    patches.iter().find(|p| {
        let [[u0, u1], [v0, v1]] = &p.domain;
        points
            .iter()
            .all(|(u, v)| u0 <= u && u <= u1 && v0 <= v && v <= v1)
    })
}

/// `-∫ G dū` along one piece in one patch.
fn piece_flux<T: Real>(patch: &Patch<T>, uv: &[Bern<T>; 3], uniform: bool) -> Option<T> {
    let [[u0, u1], [v0, v1]] = &patch.domain;
    let [pu, pv, pw] = uv;
    let local = |x: &Bern<T>, lo: &R, hi: &R| {
        scaled(
            &difference(x, &scaled(pw, &c(lo))),
            &c(&(ratio(1, 1) / (hi - lo))),
        )
    };
    let (big_u, big_v) = (local(pu, u0, u1), local(pv, v0, v1));
    let (m, n) = (patch.g.len() - 1, patch.g[0].len() - 1);
    let k = composed(&patch.g, &big_u, &big_v, pw);
    // dū = (Ū' W - Ū W') / W^2.
    let du = difference(
        &product(&derivative(&big_u), pw),
        &product(&big_u, &derivative(pw)),
    );
    let integrand: Bern<T> = product(&k, &du).iter().map(|x| x.neg()).collect();
    quotient_integral(&integrand, pw, (m + n + 2) as i32, uniform)
}

/// The flux of a face on a nonrational, nonperiodic spline surface, or
/// `None` when it is not decided.
pub(super) fn spline_face_flux<T: Real>(surface: &BSplineSurface3, loops: &[Lp]) -> Option<T> {
    if surface.is_rational() || surface.u_knots().is_periodic() || surface.v_knots().is_periodic() {
        return None;
    }
    let bezier = surface.bezier_patches().ok()?;
    let mut patches: Vec<Patch<T>> = Vec::new();
    for q in &bezier {
        let h = antiderivative::<T>(q);
        // The columns below: patches of the same u-span with lower v.
        let [us, vs] = q.domain().clone();
        let mut g = h;
        for below in &bezier {
            let [bu, bv] = below.domain();
            if *bu == us && bv[1] <= vs[0] {
                let col = column(&antiderivative::<T>(below));
                for (row, x) in g.iter_mut().zip(&col) {
                    for y in row.iter_mut() {
                        *y = y.add(x);
                    }
                }
            }
        }
        patches.push(Patch {
            domain: q.domain().clone(),
            g,
        });
    }
    let mut total = T::exact_f64(0.0);
    for lp in loops {
        for fin in &lp.fins {
            let arcs = pcurve_arcs(&fin.pcurve)?;
            // Split lines where they cross a patch boundary.
            let mut cuts = vec![ratio(0, 1), ratio(1, 1)];
            for (from, to, _) in &arcs {
                cuts.extend([from.clone(), to.clone()]);
            }
            if let Curve2::LineSegment { start, end } = &fin.pcurve {
                let (a, b) = ([r(start.x), r(start.y)], [r(end.x), r(end.y)]);
                for p in &patches {
                    for (axis, [lo, hi]) in p.domain.iter().enumerate() {
                        let d = &b[axis] - &a[axis];
                        if d == ratio(0, 1) {
                            continue;
                        }
                        for k in [lo, hi] {
                            let t = (k - &a[axis]) / &d;
                            if t > ratio(0, 1) && t < ratio(1, 1) {
                                cuts.push(t);
                            }
                        }
                    }
                }
            }
            cuts.sort();
            cuts.dedup();
            for w in cuts.windows(2) {
                let piece = piece_of(&arcs, &w[0], &w[1])?;
                let points: Vec<(R, R)> = (0..piece[3].len())
                    .map(|i| (&piece[0][i] / &piece[3][i], &piece[1][i] / &piece[3][i]))
                    .collect();
                let patch = locate(&patches, &points)?;
                let uniform = piece[3].iter().all(|x| *x == piece[3][0]);
                let uv = [lift(&piece[0]), lift(&piece[1]), lift(&piece[3])];
                total = total.add(&piece_flux(patch, &uv, uniform)?);
            }
        }
        // Chords closing the loop's gaps: their ends must certainly share a
        // patch.
        for (a, b) in super::chords::<T>(lp) {
            let inside = |p: &Patch<T>, x: &[T; 2]| {
                p.domain.iter().zip(x).all(|([lo, hi], v)| {
                    matches!(v.cmp(&c(lo)), Some(Ordering::Greater | Ordering::Equal))
                        && matches!(v.cmp(&c(hi)), Some(Ordering::Less | Ordering::Equal))
                })
            };
            let patch = patches.iter().find(|p| inside(p, &a) && inside(p, &b))?;
            let one = vec![T::exact_f64(1.0), T::exact_f64(1.0)];
            let uv = [
                vec![a[0].clone(), b[0].clone()],
                vec![a[1].clone(), b[1].clone()],
                one,
            ];
            total = total.add(&piece_flux(patch, &uv, true)?);
        }
    }
    Some(total)
}

/// Strips of the `v` domain for the enclosed flux.
const STRIPS: usize = 32;
/// Pieces of each line or spline pcurve piece for the enclosed flux.
const PIECES: usize = 5;

/// The flux of a face on any nonperiodic spline surface (rational, or with
/// pcurve pieces across patches), enclosed (S4d): with
/// `f = S·(S_u × S_v)` enclosed over boxes from the surfaces' jets (the
/// quotient rule covers weights), the antiderivative in `v` from the
/// domain's start over a box is the sum over the strips below it of
/// `Δs f(box × strip)` plus the partial strip, and the face's flux is
/// `-∮ G du`, enclosed piece by piece. Widths are `O(Δs)` and `O(h)`.
pub(super) fn enclosed_face_flux<T: Real>(surface: &BSplineSurface3, loops: &[Lp]) -> Option<T> {
    use super::spline_taylor::{lift_patches, spline_jet, Patch};
    if surface.u_knots().is_periodic() || surface.v_knots().is_periodic() {
        return None;
    }
    let patches: Vec<Patch<T>> = lift_patches(&surface.bezier_patches().ok()?);
    let ((_, _), (va, vb)) = surface.domain();
    let (va, vb) = (r(va), r(vb));
    let edges: Vec<R> = (0..=STRIPS)
        .map(|k| &va + (&vb - &va) * ratio(k as i64, STRIPS as i64))
        .collect();
    let f = |u: &T, v: &T| -> Option<T> {
        let [s, su, sv, ..] = spline_jet(&patches, u, v)?;
        let n: [T; 3] = std::array::from_fn(|k| {
            let (a, b) = ((k + 1) % 3, (k + 2) % 3);
            su[a].mul(&sv[b]).sub(&su[b].mul(&sv[a]))
        });
        Some(s[0].mul(&n[0]).add(&s[1].mul(&n[1])).add(&s[2].mul(&n[2])))
    };
    // G over a box: whole strips below the box's lowest v, then the part
    // from that strip's start to its highest v.
    let g = |u: &T, v: &T| -> Option<T> {
        let (lo, hi) = super::bernstein::ends(v);
        let mut total = T::exact_f64(0.0);
        // Below the domain's start, G is negative: minus the integral up
        // to it.
        if lo < edges[0] {
            let top = hi.clone().min(edges[0].clone());
            let below = c::<T>(&lo).union(&c(&edges[0]));
            let length = c::<T>(&(&lo - &edges[0])).union(&c(&(&top - &edges[0])));
            total = total.add(&f(u, &below)?.mul(&length));
        }
        let mut k = 0;
        while k < STRIPS && edges[k + 1] <= lo {
            let strip = c::<T>(&edges[k]).union(&c(&edges[k + 1]));
            total = total.add(&f(u, &strip)?.mul(&c(&(&edges[k + 1] - &edges[k]))));
            k += 1;
        }
        let start = edges[k.min(STRIPS)].clone();
        if hi > start {
            let rest = c::<T>(&start).union(&c(&hi));
            let length = T::exact_f64(0.0).union(&c(&(&hi - &start)));
            total = total.add(&f(u, &rest)?.mul(&length));
        }
        Some(total)
    };
    let mut total = T::exact_f64(0.0);
    for lp in loops {
        for fin in &lp.fins {
            match &fin.pcurve {
                Curve2::LineSegment { start, end } => {
                    // Pieces of the segment: -du times G over each box.
                    let (a, b) = ([r(start.x), r(start.y)], [r(end.x), r(end.y)]);
                    let n = 1 << PIECES;
                    for k in 0..n {
                        let at = |t: R| -> [R; 2] {
                            std::array::from_fn(|i| &a[i] + (&b[i] - &a[i]) * &t)
                        };
                        let (p, q) = (at(ratio(k, n)), at(ratio(k + 1, n)));
                        let du = c::<T>(&(&q[0] - &p[0]));
                        let bu = c::<T>(&p[0]).union(&c(&q[0]));
                        let bv = c::<T>(&p[1]).union(&c(&q[1]));
                        total = total.add(&g(&bu, &bv)?.mul(&du).neg());
                    }
                }
                Curve2::BSpline(spline) => {
                    total = total.add(&super::bernstein::green_integral(spline, PIECES, &g)?.neg());
                }
                Curve2::CircularArc { .. } => return None,
            }
        }
        for (a, b) in super::chords::<T>(lp) {
            total = total.add(&super::chord_enclosure(&a, &b, &g)?);
        }
    }
    Some(total)
}
