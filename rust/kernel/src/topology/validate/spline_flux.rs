//! Green's theorem on spline surfaces (S4d of REVIEW_NOTES.md): the
//! orientation flux `∫∫ S·(S_u × S_v) du dv` of a face, and its mass
//! integrands. On a nonrational surface, on each Bézier patch in local
//! coordinates `(ū, v̄)` on the unit square, an integrand polynomial in the
//! position `X` (poles divided by the common weight) and `X_ū × X_v̄` is a
//! tensor Bernstein polynomial; its antiderivative `H(ū, v̄)` from
//! `v̄ = 0`, plus the full columns below the patch, `G = Σ H(ū, 1) +
//! H(ū, v̄)`, is the global antiderivative in `v` in local units. The
//! integral is `-∮ G dū` along the face's loops, each piece in one patch: a
//! line pcurve is split exactly where it crosses a patch boundary, a spline
//! pcurve's Bézier piece must lie in one patch by its control points, and a
//! chord's ends must certainly share a patch. Otherwise (a rational surface,
//! a piece across a patch boundary, an integrand with `|N|`) the
//! antiderivative is enclosed over strips of the `v` domain.
use super::bernstein::{
    c, derivative, difference, lift, pcurve_arcs, piece_of, power, product, quotient_integral, r,
    ratio, scaled, sum, value, Bern,
};
use super::spline_taylor::{lift_patches_about, spline_jet1, Patch as Jets};
use super::Lp;
use crate::certified::Real;
use crate::jet::Jet;
use crate::surface::ExactBezierSurface3;
use crate::topology::Curve2;
use crate::BSplineSurface3;
use num_rational::BigRational as R;

/// `S` (relative to the integration's origin), `S_u` and `S_v` over a box.
type Jet1<T> = [[T; 3]; 3];

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

/// The position relative to `origin` and the parametric normal
/// `X_ū × X_v̄` of a nonrational patch, as tensors in local coordinates.
fn patch_frame<T: Real>(
    patch: &ExactBezierSurface3,
    origin: &[T; 3],
) -> ([Tensor<T>; 3], [Tensor<T>; 3]) {
    let [du, dv] = patch.degrees();
    let poles = patch.homogeneous_poles();
    let x: [Tensor<T>; 3] = std::array::from_fn(|k| {
        (0..=du)
            .map(|i| {
                (0..=dv)
                    .map(|j| {
                        let p = &poles[i * (dv + 1) + j];
                        c::<T>(&(&p[k] / &p[3])).sub(&origin[k])
                    })
                    .collect()
            })
            .collect()
    });
    let xu = x.clone().map(|t| partial_u(&t));
    let xv = x.clone().map(|t| partial_v(&t));
    let cross: [Tensor<T>; 3] = std::array::from_fn(|k| {
        let (a, b) = ((k + 1) % 3, (k + 2) % 3);
        tensor_sum(
            &tensor_product(&xu[a], &xv[b]),
            &tensor_scaled(&tensor_product(&xu[b], &xv[a]), &T::exact_f64(-1.0)),
        )
    });
    (x, cross)
}

fn tensor_scaled<T: Real>(t: &Tensor<T>, s: &T) -> Tensor<T> {
    t.iter()
        .map(|row| row.iter().map(|q| q.mul(s)).collect())
        .collect()
}

/// The flux integrand `(X - origin)·(X_ū × X_v̄)` of a patch.
fn flux_integrand<T: Real>(patch: &ExactBezierSurface3, origin: &[T; 3]) -> Vec<Tensor<T>> {
    let (x, n) = patch_frame::<T>(patch, origin);
    let f = (1..3).fold(tensor_product(&x[0], &n[0]), |acc, k| {
        tensor_sum(&acc, &tensor_product(&x[k], &n[k]))
    });
    vec![f]
}

/// The ten volume and moment integrands of the mass properties (the
/// polynomial ones, in `mass.rs`'s order) of a patch, relative to `origin`.
fn moment_integrands<T: Real>(patch: &ExactBezierSurface3, origin: &[T; 3]) -> Vec<Tensor<T>> {
    let (p, n) = patch_frame(patch, origin);
    let third = c::<T>(&ratio(1, 3));
    let half = T::exact_f64(0.5);
    let sq = |i: usize| tensor_product(&p[i], &p[i]);
    let volume = tensor_scaled(
        &tensor_sum(
            &tensor_sum(&tensor_product(&p[0], &n[0]), &tensor_product(&p[1], &n[1])),
            &tensor_product(&p[2], &n[2]),
        ),
        &third,
    );
    let first = |i: usize| tensor_scaled(&tensor_product(&sq(i), &n[i]), &half);
    let second = |i: usize| {
        tensor_scaled(
            &tensor_product(&tensor_product(&sq(i), &p[i]), &n[i]),
            &third,
        )
    };
    let mixed = |i: usize, j: usize| {
        tensor_scaled(
            &tensor_product(&tensor_product(&sq(i), &p[j]), &n[i]),
            &half,
        )
    };
    vec![
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
    ]
}

/// `H(ū, v̄) = ∫_0^v̄ f(ū, s) ds` along every row:
/// `∫_0^v̄ Σ c_j B_j^n = Σ_k (Σ_{j<k} c_j)/(n+1) B_k^{n+1}`.
fn v_antiderivative<T: Real>(f: &Tensor<T>) -> Tensor<T> {
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
/// The Bernstein basis of degree `m` in the homogeneous local coordinate
/// `X` (with `P_w - X` its complement), each `C(m, i) X^i (P_w - X)^(m-i)`.
fn basis<T: Real>(m: usize, x: &Bern<T>, rest: &Bern<T>) -> Vec<Bern<T>> {
    let binomial = |n: usize, i: usize| {
        (0..i).fold(ratio(1, 1), |acc, j| {
            acc * ratio((n - j) as i64, (j + 1) as i64)
        })
    };
    (0..=m)
        .map(|i| {
            scaled(
                &product(&power(x, i), &power(rest, m - i)),
                &c(&binomial(m, i)),
            )
        })
        .collect()
}

/// `K(ū(τ), v̄(τ)) P_w^(m + n)` from the bases of `ū` and `v̄` of the
/// tensor's degrees.
fn composed<T: Real>(k: &Tensor<T>, bu: &[Bern<T>], bv: &[Bern<T>]) -> Bern<T> {
    let (m, n) = (k.len() - 1, k[0].len() - 1);
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
    /// `H` plus the full columns below, as `G` in local units, per
    /// integrand.
    g: Vec<Tensor<T>>,
}

/// The patch holding every one of a piece's points, if any; else one
/// holding its ends whose box the piece's curve certainly stays in (its
/// homogeneous coordinates against each bound, exactly nonnegative: a
/// control polygon may leave a box its curve keeps to, as a crease's does
/// where it nearly touches a cap).
fn locate<'a, T: Real>(
    patches: &'a [Patch<T>],
    points: &[(R, R)],
    piece: &[Vec<R>; 4],
) -> Option<&'a Patch<T>> {
    let holds = |p: &Patch<T>, pts: &[(R, R)]| {
        let [[u0, u1], [v0, v1]] = &p.domain;
        pts.iter()
            .all(|(u, v)| u0 <= u && u <= u1 && v0 <= v && v <= v1)
    };
    if let Some(p) = patches.iter().find(|p| holds(p, points)) {
        return Some(p);
    }
    let ends = [points[0].clone(), points[points.len() - 1].clone()];
    patches.iter().find(|p| {
        holds(p, &ends)
            && p.domain.iter().enumerate().all(|(axis, [lo, hi])| {
                let x = &piece[axis];
                let w = &piece[3];
                let above: Vec<R> = x.iter().zip(w).map(|(x, w)| x - lo * w).collect();
                let below: Vec<R> = x.iter().zip(w).map(|(x, w)| hi * w - x).collect();
                super::bernstein::nonnegative(&above) && super::bernstein::nonnegative(&below)
            })
    })
}

/// A tensor's value over a box of its local coordinates (de Casteljau in
/// interval arithmetic, any box, inside the unit square or not).
fn tensor_at<T: Real>(g: &Tensor<T>, u: &T, v: &T) -> T {
    let along_v: Bern<T> = g.iter().map(|row| value(row, v)).collect();
    value(&along_v, u)
}

/// The patch holding a spline pcurve piece's points but for slivers across
/// its `u` boundaries at most `2^-20` of its width (a crease's pcurve over
/// a wall's knot line: the identity in `u` rounds, so a piece ending on a
/// knot line at a rounding step past it, S8b.3's (c) and S9f.1), and per
/// integrand a bound of the error of integrating the slivers on that
/// patch's polynomial: on a sliver `δ` wide the piece's `u` turns at most
/// `deg + 1` times, so `∫|du| ≤ 2 (deg + 1) δ` there, and the true `-∫ G
/// dū` (the neighbour's `G`) and the computed one each lie within that
/// times their `|G|` over the sliver's box. `None` for a piece across a `v`
/// boundary or beyond the surface.
fn sliver_patch<'a, T: Real>(
    patches: &'a [Patch<T>],
    points: &[(R, R)],
    piece: &[Vec<R>; 4],
) -> Option<(&'a Patch<T>, Vec<T>)> {
    let span = |f: &dyn Fn(&(R, R)) -> &R| -> [R; 2] {
        let lo = points.iter().map(f).min().expect("a point").clone();
        let hi = points.iter().map(f).max().expect("a point").clone();
        [lo, hi]
    };
    let (mut us, mut vs) = (span(&|p| &p.0), span(&|p| &p.1));
    // Inside the surface's domain: the hull, or the curve where its hull
    // leaves the domain (a crease nearly touching a cap, as `locate`).
    for (axis, range) in [&mut us, &mut vs].into_iter().enumerate() {
        let lo = patches.iter().map(|q| &q.domain[axis][0]).min()?.clone();
        let hi = patches.iter().map(|q| &q.domain[axis][1]).max()?.clone();
        let (x, w) = (&piece[axis], &piece[3]);
        if range[0] < lo {
            let off: Vec<R> = x.iter().zip(w).map(|(x, w)| x - &lo * w).collect();
            if !super::bernstein::nonnegative(&off) {
                return None;
            }
            range[0] = lo;
        }
        if range[1] > hi {
            let off: Vec<R> = x.iter().zip(w).map(|(x, w)| &hi * w - x).collect();
            if !super::bernstein::nonnegative(&off) {
                return None;
            }
            range[1] = hi;
        }
    }
    let two = ratio(2, 1);
    let mid = (&us[0] + &us[1]) / &two;
    let sliver = ratio(1, 1 << 20);
    let holds_v = |q: &Patch<T>| q.domain[1][0] <= vs[0] && vs[1] <= q.domain[1][1];
    let patch = patches
        .iter()
        .find(|q| holds_v(q) && q.domain[0][0] <= mid && mid <= q.domain[0][1])?;
    let [u0, u1] = &patch.domain[0];
    let width = u1 - u0;
    let limit = &width * &sliver;
    if u0 - &us[0] > limit || &us[1] - u1 > limit {
        return None;
    }
    let degree = piece[0].len().saturating_sub(1) as i64;
    let count = patch.g.len();
    let mut bound = vec![ratio(0, 1); count];
    // |G| over a box of global (u, v), in a patch's local units.
    let sup = |q: &Patch<T>, ua: &R, ub: &R| -> Option<Vec<R>> {
        let [[a0, a1], [b0, b1]] = &q.domain;
        let local = |x: &R, lo: &R, hi: &R| c::<T>(&((x - lo) / (hi - lo)));
        let u = local(ua, a0, a1).union(&local(ub, a0, a1));
        let v = local(&vs[0], b0, b1).union(&local(&vs[1], b0, b1));
        q.g.iter()
            .map(|g| {
                let (lo, hi) = tensor_at(g, &u, &v).bounds_f64();
                let m = lo.abs().max(hi.abs());
                if m.is_finite() {
                    R::from_float(m)
                } else {
                    None
                }
            })
            .collect()
    };
    for (outside, edge, other_edge) in [
        (u0 - &us[0], u0.clone(), 1usize),
        (&us[1] - u1, u1.clone(), 0usize),
    ] {
        if outside <= ratio(0, 1) {
            continue;
        }
        let neighbour = patches
            .iter()
            .find(|q| q.domain[0][other_edge] == edge && holds_v(q))?;
        let (lo, hi) = if other_edge == 1 {
            (us[0].clone(), edge.clone())
        } else {
            (edge.clone(), us[1].clone())
        };
        let mine = sup(patch, &lo, &hi)?;
        let theirs = sup(neighbour, &lo, &hi)?;
        let wn = &neighbour.domain[0][1] - &neighbour.domain[0][0];
        let turns = ratio(2 * (degree + 1), 1) * &outside;
        for k in 0..count {
            bound[k] += &turns * (&mine[k] / &width + &theirs[k] / &wn);
        }
    }
    let extra = bound.iter().map(|e| T::exact_f64(0.0).widen(e)).collect();
    Some((patch, extra))
}

/// `-∫ G dū` along one piece in one patch, for every integrand. The piece
/// is homogeneous `(U, V, W)` in global `(u, v)`, as exact Bernstein
/// coordinates. With `high`, a rational piece is integrated by the certified
/// quadrature (F8) where it applies.
fn piece_integrals<T: Real>(
    patch: &Patch<T>,
    uv: &[Bern<T>; 3],
    uniform: bool,
    high: bool,
) -> Option<Vec<T>> {
    let [[u0, u1], [v0, v1]] = &patch.domain;
    let [pu, pv, pw] = uv;
    let local = |x: &Bern<T>, lo: &R, hi: &R| {
        scaled(
            &difference(x, &scaled(pw, &c(lo))),
            &c(&(ratio(1, 1) / (hi - lo))),
        )
    };
    let (big_u, big_v) = (local(pu, u0, u1), local(pv, v0, v1));
    // dū = (Ū' W - Ū W') / W^2.
    let du = difference(
        &product(&derivative(&big_u), pw),
        &product(&big_u, &derivative(pw)),
    );
    // The integrands share degrees: each basis is formed once per piece.
    let (rest_u, rest_v) = (difference(pw, &big_u), difference(pw, &big_v));
    let mut bases_u: std::collections::BTreeMap<usize, Vec<Bern<T>>> = Default::default();
    let mut bases_v: std::collections::BTreeMap<usize, Vec<Bern<T>>> = Default::default();
    for g in &patch.g {
        let (m, n) = (g.len() - 1, g[0].len() - 1);
        bases_u
            .entry(m)
            .or_insert_with(|| basis(m, &big_u, &rest_u));
        bases_v
            .entry(n)
            .or_insert_with(|| basis(n, &big_v, &rest_v));
    }
    patch
        .g
        .iter()
        .map(|g| {
            let (m, n) = (g.len() - 1, g[0].len() - 1);
            let k = composed(g, &bases_u[&m], &bases_v[&n]);
            let integrand: Bern<T> = product(&k, &du).iter().map(|x| x.neg()).collect();
            let power = (m + n + 2) as i32;
            if high && !uniform {
                if let Some(x) = super::quadrature::quotient_integral(&integrand, pw, power) {
                    return Some(x);
                }
            }
            quotient_integral(&integrand, pw, power, uniform)
        })
        .collect()
}

fn add_all<T: Real>(total: &mut Option<Vec<T>>, values: Vec<T>) {
    *total = Some(match total.take() {
        None => values,
        Some(t) => t.iter().zip(&values).map(|(a, b)| a.add(b)).collect(),
    });
}

/// `-∮ G dū` of the integrands `integrands(patch)` over a face on a
/// nonrational, nonperiodic spline surface, exactly up to the tier's
/// rounding (with `high`, rational pcurve pieces by the certified
/// quadrature), or `None` when a piece does not lie in one patch.
fn exact_green<T: Real>(
    surface: &BSplineSurface3,
    loops: &[Lp],
    integrands: &dyn Fn(&ExactBezierSurface3) -> Vec<Tensor<T>>,
    high: bool,
) -> Option<Vec<T>> {
    if surface.is_rational() || surface.u_knots().is_periodic() || surface.v_knots().is_periodic() {
        return None;
    }
    let bezier = surface.bezier_patches().ok()?;
    let anti: Vec<Vec<Tensor<T>>> = bezier
        .iter()
        .map(|q| integrands(q).iter().map(v_antiderivative).collect())
        .collect();
    let mut patches: Vec<Patch<T>> = Vec::new();
    for (q, h) in bezier.iter().zip(&anti) {
        // The columns below: patches of the same u-span with lower v.
        let [us, vs] = q.domain().clone();
        let mut g = h.clone();
        for (below, hb) in bezier.iter().zip(&anti) {
            let [bu, bv] = below.domain();
            if *bu == us && bv[1] <= vs[0] {
                for (gk, hk) in g.iter_mut().zip(hb) {
                    let col = column(hk);
                    for (row, x) in gk.iter_mut().zip(&col) {
                        for y in row.iter_mut() {
                            *y = y.add(x);
                        }
                    }
                }
            }
        }
        patches.push(Patch {
            domain: q.domain().clone(),
            g,
        });
    }
    let mut total: Option<Vec<T>> = None;
    for lp in loops {
        for fin in &lp.fins {
            // A spline wall's meeting's projection onto this wall (S9f.2b):
            // `-G(ū, v̄) dū` along it, piece by piece between the knots on
            // each knot span's patch, with jets.
            if let Curve2::Projection(pr) = &fin.pcurve {
                add_all(&mut total, wall_green(&patches, pr, surface)?);
                continue;
            }
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
                let (patch, extra) = match locate(&patches, &points, &piece) {
                    Some(patch) => (patch, None),
                    None => {
                        let (patch, extra) = sliver_patch(&patches, &points, &piece)?;
                        (patch, Some(extra))
                    }
                };
                let uniform = piece[3].iter().all(|x| *x == piece[3][0]);
                let uv = [lift(&piece[0]), lift(&piece[1]), lift(&piece[3])];
                add_all(&mut total, piece_integrals(patch, &uv, uniform, high)?);
                if let Some(extra) = extra {
                    add_all(&mut total, extra);
                }
            }
        }
        // Chords closing the loop's gaps: their ends must certainly share a
        // patch (S9f.2b: within a sliver of it, `chord_patch`).
        for (a, b) in super::chords::<T>(lp) {
            let at = super::chord_patch(patches.iter().map(|p| &p.domain), &a, &b)?;
            let patch = &patches[at];
            let one = vec![T::exact_f64(1.0), T::exact_f64(1.0)];
            let uv = [
                vec![a[0].clone(), b[0].clone()],
                vec![a[1].clone(), b[1].clone()],
                one,
            ];
            add_all(&mut total, piece_integrals(patch, &uv, true, high)?);
        }
    }
    let count = patches.first()?.g.len();
    Some(total.unwrap_or_else(|| vec![T::exact_f64(0.0); count]))
}

/// `-∫ G dū` along a spline wall's meeting's projection onto its own wall
/// (S9f.2b), for every integrand: on each piece between the wall's knots
/// (`projection::wall_pieces`, exact splits) the patch of its knot span,
/// `G` there at the piece's local `(ū, v̄)` jets.
fn wall_green<T: Real>(
    patches: &[Patch<T>],
    pr: &crate::topology::Projection,
    surface: &BSplineSurface3,
) -> Option<Vec<T>> {
    let crate::topology::Curve3::WallMeet(m) = &pr.curve else {
        return None;
    };
    if m.wall != *surface {
        return None;
    }
    let spans = super::wall_meet::spans(m)?;
    let count = patches.first()?.g.len();
    let mut totals = vec![T::exact_f64(0.0); count];
    super::projection::wall_pieces(
        pr,
        count,
        true,
        &|k, u, v, du, _| {
            let span = spans.spans.get(k)?;
            let patch = patches.iter().find(|q| q.domain[0] == span.u)?;
            let [[u0, u1], [v0, v1]] = &patch.domain;
            let iu = c::<T>(&(ratio(1, 1) / (u1 - u0)));
            let iv = c::<T>(&(ratio(1, 1) / (v1 - v0)));
            let ub = u.add_constant(&c::<T>(&-u0)).scale(&iu);
            let vb = v.add_constant(&c::<T>(&-v0)).scale(&iv);
            let dub = du.scale(&iu);
            Some(
                patch
                    .g
                    .iter()
                    .map(|g| tensor_jet(g, &ub, &vb).mul(&dub).neg())
                    .collect(),
            )
        },
        &mut totals,
    )?;
    Some(totals)
}

/// A Bernstein polynomial at a jet by de Casteljau: about a point base in
/// the form `(1 - t) a + t b`, whose enclosures do not grow with the levels
/// where `t` lies in `[0, 1]` (`a + t (b - a)` widens them by up to `1 + 2
/// t` a level: a tensor of degree nine `10^4` times); over a range as `a +
/// t (b - a)` (`wall_meet::bernstein`'s choice).
fn casteljau_jet<T: Real>(mut row: Vec<Jet<T>>, t: &Jet<T>) -> Jet<T> {
    let point = super::quadrature::Num::<T>::sharp(t);
    let s = Jet::constant(T::exact_f64(1.0), t.order()).sub(t);
    for last in (1..row.len()).rev() {
        for i in 0..last {
            row[i] = if point {
                row[i].mul(&s).add(&row[i + 1].mul(t))
            } else {
                row[i].add(&t.mul(&row[i + 1].sub(&row[i])))
            };
        }
    }
    row.swap_remove(0)
}

/// A tensor Bernstein polynomial at jets of `(ū, v̄)` (rows in `v̄`, then
/// `ū`), as `casteljau_jet`.
fn tensor_jet<T: Real>(g: &Tensor<T>, u: &Jet<T>, v: &Jet<T>) -> Jet<T> {
    let n = u.order();
    casteljau_jet(
        g.iter()
            .map(|row| casteljau_jet(row.iter().map(|x| Jet::constant(x.clone(), n)).collect(), v))
            .collect(),
        u,
    )
}

/// The flux of a face on a nonrational, nonperiodic spline surface
/// relative to `origin`, or `None` when it is not decided.
pub(super) fn spline_face_flux<T: Real>(
    surface: &BSplineSurface3,
    loops: &[Lp],
    origin: &[T; 3],
) -> Option<T> {
    exact_green(surface, loops, &|q| flux_integrand(q, origin), false)?.pop()
}

/// Strips of the `v` domain for the enclosed flux.
const STRIPS: usize = 32;
/// Pieces of each line or spline pcurve piece for the enclosed flux.
const PIECES: usize = 5;

/// `-∮ G du` of the integrands `f(jet)` (from a surface jet over a box)
/// over a face on any nonperiodic spline surface (rational, or with pcurve
/// pieces across patches), enclosed (S4d): the antiderivative in `v` from
/// the domain's start over a box is the sum over the strips below it of
/// `Δs f(box × strip)` plus the partial strip, and `-∮ G du` is enclosed
/// piece by piece. Widths are `O(Δs)` and `O(h)`.
fn enclosed_green<T: Real>(
    surface: &BSplineSurface3,
    loops: &[Lp],
    origin: &[T; 3],
    count: usize,
    f: &dyn Fn(&Jet1<T>) -> Option<Vec<T>>,
) -> Option<Vec<T>> {
    if surface.u_knots().is_periodic() || surface.v_knots().is_periodic() {
        return None;
    }
    // The jets' position is relative to `origin` (a point: its midpoint,
    // plus the zero remainder).
    let centre: [R; 3] = std::array::from_fn(|k| origin[k].midpoint());
    let rest: [T; 3] = std::array::from_fn(|k| c::<T>(&centre[k]).sub(&origin[k]));
    let patches: Vec<Jets<T>> = lift_patches_about(&surface.bezier_patches().ok()?, Some(&centre));
    let ((_, _), (va, vb)) = surface.domain();
    let (va, vb) = (r(va), r(vb));
    let edges: Vec<R> = (0..=STRIPS)
        .map(|k| &va + (&vb - &va) * ratio(k as i64, STRIPS as i64))
        .collect();
    let at = |u: &T, v: &T| -> Option<Vec<T>> {
        let mut jet = spline_jet1(&patches, u, v)?;
        for (x, r) in jet[0].iter_mut().zip(&rest) {
            *x = x.add(r);
        }
        f(&jet)
    };
    let axpy = |total: &mut Vec<T>, values: Vec<T>, s: &T| {
        for (t, x) in total.iter_mut().zip(values) {
            *t = t.add(&x.mul(s));
        }
    };
    // G over a box: whole strips below the box's lowest v, then the part
    // from that strip's start to its highest v.
    let g = |u: &T, v: &T| -> Option<Vec<T>> {
        let (lo, hi) = super::bernstein::ends(v);
        let mut total = vec![T::exact_f64(0.0); count];
        // Below the domain's start, G is negative: minus the integral up
        // to it.
        if lo < edges[0] {
            let top = hi.clone().min(edges[0].clone());
            let below = c::<T>(&lo).union(&c(&edges[0]));
            let length = c::<T>(&(&lo - &edges[0])).union(&c(&(&top - &edges[0])));
            axpy(&mut total, at(u, &below)?, &length);
        }
        let mut k = 0;
        while k < STRIPS && edges[k + 1] <= lo {
            let strip = c::<T>(&edges[k]).union(&c(&edges[k + 1]));
            axpy(&mut total, at(u, &strip)?, &c(&(&edges[k + 1] - &edges[k])));
            k += 1;
        }
        let start = edges[k.min(STRIPS)].clone();
        if hi > start {
            let rest = c::<T>(&start).union(&c(&hi));
            let length = T::exact_f64(0.0).union(&c(&(&hi - &start)));
            axpy(&mut total, at(u, &rest)?, &length);
        }
        Some(total)
    };
    let mut total = vec![T::exact_f64(0.0); count];
    let minus_one = T::exact_f64(-1.0);
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
                        axpy(&mut total, g(&bu, &bv)?, &du.neg());
                    }
                }
                Curve2::BSpline(spline) => {
                    let values = super::bernstein::green_integrals(spline, PIECES, &g)?;
                    axpy(&mut total, values, &minus_one);
                }
                Curve2::CircularArc { .. }
                | Curve2::EllipseArc { .. }
                | Curve2::Sinusoid { .. }
                | Curve2::Projection(_) => return None,
            }
        }
        for (a, b) in super::chords::<T>(lp) {
            // -du times G over the chord's box, which holds its mean.
            let du = b[0].sub(&a[0]);
            let values = g(&a[0].union(&b[0]), &a[1].union(&b[1]))?;
            axpy(&mut total, values, &du.neg());
        }
    }
    Some(total)
}

/// The flux of a face on any nonperiodic spline surface relative to
/// `origin`, enclosed by strips (S4d).
pub(super) fn enclosed_face_flux<T: Real>(
    surface: &BSplineSurface3,
    loops: &[Lp],
    origin: &[T; 3],
) -> Option<T> {
    let flux = |jet: &Jet1<T>| {
        let [p, su, sv] = jet;
        let n = cross(su, sv);
        Some(vec![p[0]
            .mul(&n[0])
            .add(&p[1].mul(&n[1]))
            .add(&p[2].mul(&n[2]))])
    };
    enclosed_green(surface, loops, origin, 1, &flux)?.pop()
}

fn cross<T: Real>(a: &[T; 3], b: &[T; 3]) -> [T; 3] {
    std::array::from_fn(|k| {
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        a[i].mul(&b[j]).sub(&a[j].mul(&b[i]))
    })
}

/// The fourteen mass integrals of a face on a nonperiodic spline surface
/// relative to `origin` (S4d): on a nonrational surface whose pieces each
/// lie in one patch, the ten volume and moment terms exactly and the four
/// with `|N|` by the certified quadrature (F8); otherwise all fourteen by
/// the quadrature. Where it cannot run, the strips.
pub(super) fn spline_face_integrals<T: Real>(
    surface: &BSplineSurface3,
    loops: &[Lp],
    origin: &[T; 3],
) -> Option<Vec<T>> {
    let terms = |jet: &Jet1<T>, all: bool| {
        let [p, su, sv] = jet;
        let n = cross(su, sv);
        let norm = n[0].square().add(&n[1].square()).add(&n[2].square()).sqrt();
        let values = super::mass::point_integrands(p, &n, &norm);
        Some(if all {
            values.to_vec()
        } else {
            values[10..].to_vec()
        })
    };
    use super::quadrature::spline_face;
    if let Some(mut exact) = exact_green(surface, loops, &|q| moment_integrands(q, origin), true) {
        exact.extend(
            spline_face(surface, loops, origin, false)
                .or_else(|| enclosed_green(surface, loops, origin, 4, &|jet| terms(jet, false)))?,
        );
        return Some(exact);
    }
    spline_face(surface, loops, origin, true)
        .or_else(|| enclosed_green(surface, loops, origin, 14, &|jet| terms(jet, true)))
}
