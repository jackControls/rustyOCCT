//! Constrained Delaunay triangulation of a face's planar chart.
//!
//! Points are binary64 chart coordinates; every combinatorial decision
//! (point location, which side of a segment, convexity, crossings) uses the
//! exact orientation predicate, so the triangulation is valid whatever the
//! rounding of the chart. The in-circle test only improves shape: it is
//! Shewchuk's filtered determinant, and an edge is flipped only when the
//! filter certifies the violation, so Lawson flipping terminates and never
//! cycles on cocircular points. Constraints are recovered by Sloan's
//! flipping; a constraint that crosses another, repeats one or passes
//! through a vertex is a conflict for the caller (its chords are too
//! coarse), never repaired here. Triangles are classified inside the domain
//! by the parity of the constraints crossed from the enclosing triangle.
//! Everything is deterministic: a fixed-seed pseudo-random walk, fixed
//! orders, no hashing.
use crate::predicates::{orient2d_finite, Orientation2};
use crate::Point2;
use std::collections::VecDeque;

pub(super) const NONE: u32 = u32::MAX;
/// The enclosing triangle's corners are the first three vertices.
pub(super) const CORNERS: u32 = 3;

#[derive(Debug, Clone, Copy)]
struct Tri {
    /// Counter-clockwise.
    v: [u32; 3],
    /// Neighbour across the edge opposite `v[i]`.
    n: [u32; 3],
    /// Whether the edge opposite `v[i]` is a constraint.
    fixed: [bool; 3],
    inside: bool,
    alive: bool,
    /// Incremented on every change, so a queued check can be skipped.
    stamp: u32,
}

impl Tri {
    /// The same triangle with index `i` first.
    fn rotated(self, i: usize) -> Self {
        let r = |k: usize| (k + i) % 3;
        Self {
            v: [self.v[r(0)], self.v[r(1)], self.v[r(2)]],
            n: [self.n[r(0)], self.n[r(1)], self.n[r(2)]],
            fixed: [self.fixed[r(0)], self.fixed[r(1)], self.fixed[r(2)]],
            ..self
        }
    }
}

/// A constraint that cannot be inserted as given.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Conflict;

enum Location {
    In(u32),
    Edge(u32, usize),
    Vertex,
    Outside,
}

pub(super) struct Cdt {
    points: Vec<Point2>,
    tris: Vec<Tri>,
    /// One triangle incident to each vertex.
    star: Vec<u32>,
    seed: u64,
    last: u32,
}

fn orient(a: Point2, b: Point2, c: Point2) -> i8 {
    match orient2d_finite(a, b, c) {
        Orientation2::CounterClockwise => 1,
        Orientation2::Clockwise => -1,
        Orientation2::Collinear => 0,
    }
}

/// Whether `d` lies certainly inside the circle through the counter-clockwise
/// `a`, `b`, `c`: Shewchuk's incircle with his first error bound
/// (`iccerrboundA`); `false` when the filter cannot decide.
fn in_circle(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let (adx, ady) = (a.x - d.x, a.y - d.y);
    let (bdx, bdy) = (b.x - d.x, b.y - d.y);
    let (cdx, cdy) = (c.x - d.x, c.y - d.y);
    let (bc, cb) = (bdx * cdy, cdx * bdy);
    let (ca, ac) = (cdx * ady, adx * cdy);
    let (ab, ba) = (adx * bdy, bdx * ady);
    let alift = adx * adx + ady * ady;
    let blift = bdx * bdx + bdy * bdy;
    let clift = cdx * cdx + cdy * cdy;
    let det = alift * (bc - cb) + blift * (ca - ac) + clift * (ab - ba);
    let permanent = (bc.abs() + cb.abs()) * alift
        + (ca.abs() + ac.abs()) * blift
        + (ab.abs() + ba.abs()) * clift;
    const EPS: f64 = f64::EPSILON * 0.5;
    const BOUND: f64 = (10.0 + 96.0 * EPS) * EPS;
    det.is_finite() && permanent.is_finite() && det > BOUND * permanent
}

impl Cdt {
    /// An enclosing triangle around `points`' box, then every point in
    /// order (vertex `CORNERS + k` is `points[k]`); `Err` on a repeated point.
    pub(super) fn new(points: &[Point2]) -> Result<Self, Conflict> {
        let (mut lo, mut hi) = (
            Point2::new(f64::MAX, f64::MAX),
            Point2::new(-f64::MAX, -f64::MAX),
        );
        for p in points {
            lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
            hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
        }
        let size = (hi.x - lo.x).max(hi.y - lo.y).max(f64::MIN_POSITIVE.sqrt()) * 32.0;
        let (cx, cy) = (0.5 * (lo.x + hi.x), 0.5 * (lo.y + hi.y));
        let mut cdt = Self {
            points: vec![
                Point2::new(cx - size, cy - size),
                Point2::new(cx + size, cy - size),
                Point2::new(cx, cy + size),
            ],
            tris: vec![Tri {
                v: [0, 1, 2],
                n: [NONE; 3],
                fixed: [false; 3],
                inside: false,
                alive: true,
                stamp: 0,
            }],
            star: vec![0, 0, 0],
            seed: 0x9e37_79b9_7f4a_7c15,
            last: 0,
        };
        for p in points {
            cdt.insert(*p)?;
        }
        Ok(cdt)
    }

    pub(super) fn point(&self, v: u32) -> Point2 {
        self.points[v as usize]
    }
    pub(super) fn slots(&self) -> u32 {
        self.tris.len() as u32
    }
    pub(super) fn alive_inside(&self, t: u32) -> bool {
        let t = &self.tris[t as usize];
        t.alive && t.inside
    }
    pub(super) fn vertices(&self, t: u32) -> [u32; 3] {
        self.tris[t as usize].v
    }
    pub(super) fn fixed(&self, t: u32) -> [bool; 3] {
        self.tris[t as usize].fixed
    }
    pub(super) fn stamp(&self, t: u32) -> u32 {
        self.tris[t as usize].stamp
    }

    fn random(&mut self) -> u64 {
        self.seed ^= self.seed << 13;
        self.seed ^= self.seed >> 7;
        self.seed ^= self.seed << 17;
        self.seed
    }

    fn touch(&mut self, t: u32, touched: &mut Vec<u32>) {
        let tri = &mut self.tris[t as usize];
        tri.stamp = tri.stamp.wrapping_add(1);
        touched.push(t);
    }

    /// Visibility walk from `start` (a random edge order per step, so it
    /// cannot cycle), with a scan of every triangle as the fallback.
    fn locate(&mut self, p: Point2, start: u32) -> Location {
        let mut t = if self.tris[start as usize].alive {
            start
        } else {
            self.last
        };
        let limit = 4 * self.tris.len() + 64;
        'walk: for _ in 0..limit {
            let tri = self.tris[t as usize];
            let offset = (self.random() % 3) as usize;
            let mut zero = [false; 3];
            for k in 0..3 {
                let i = (offset + k) % 3;
                let a = self.points[tri.v[(i + 1) % 3] as usize];
                let b = self.points[tri.v[(i + 2) % 3] as usize];
                match orient(a, b, p) {
                    -1 if tri.n[i] == NONE => return Location::Outside,
                    -1 => {
                        t = tri.n[i];
                        continue 'walk;
                    }
                    0 => zero[i] = true,
                    _ => {}
                }
            }
            return Self::classify_location(t, zero);
        }
        for (k, tri) in self.tris.iter().enumerate() {
            if !tri.alive {
                continue;
            }
            let mut zero = [false; 3];
            let mut outside = false;
            for (i, z) in zero.iter_mut().enumerate() {
                let a = self.points[tri.v[(i + 1) % 3] as usize];
                let b = self.points[tri.v[(i + 2) % 3] as usize];
                match orient(a, b, p) {
                    -1 => outside = true,
                    0 => *z = true,
                    _ => {}
                }
            }
            if !outside {
                return Self::classify_location(k as u32, zero);
            }
        }
        Location::Outside
    }

    fn classify_location(t: u32, zero: [bool; 3]) -> Location {
        match zero.iter().filter(|z| **z).count() {
            0 => Location::In(t),
            1 => Location::Edge(t, zero.iter().position(|z| *z).unwrap()),
            _ => Location::Vertex,
        }
    }

    fn insert(&mut self, p: Point2) -> Result<(u32, Vec<u32>), Conflict> {
        let mut touched = Vec::new();
        let v = match self.locate(p, self.last) {
            Location::In(t) => self.split_triangle(t, p, &mut touched),
            Location::Edge(t, i) => self.split_edge(t, i, p, &mut touched),
            Location::Vertex | Location::Outside => return Err(Conflict),
        };
        Ok((v, touched))
    }

    fn add_point(&mut self, p: Point2) -> u32 {
        self.points.push(p);
        self.star.push(NONE);
        (self.points.len() - 1) as u32
    }

    fn new_tri(&mut self, tri: Tri) -> u32 {
        self.tris.push(tri);
        (self.tris.len() - 1) as u32
    }

    fn set(&mut self, slot: u32, tri: Tri) {
        let stamp = self.tris[slot as usize].stamp;
        self.tris[slot as usize] = Tri { stamp, ..tri };
        for v in tri.v {
            self.star[v as usize] = slot;
        }
    }

    /// Point the neighbour `t`'s link that referred to `old` at `new`.
    fn relink(&mut self, t: u32, old: u32, new: u32) {
        if t == NONE {
            return;
        }
        let tri = &mut self.tris[t as usize];
        for k in 0..3 {
            if tri.n[k] == old {
                tri.n[k] = new;
                return;
            }
        }
    }

    fn split_triangle(&mut self, t: u32, p: Point2, touched: &mut Vec<u32>) -> u32 {
        let v = self.add_point(p);
        let old = self.tris[t as usize];
        let [a, b, c] = old.v;
        let [na, nb, nc] = old.n;
        let [fa, fb, fc] = old.fixed;
        let t1 = self.new_tri(old);
        let t2 = self.new_tri(old);
        let make = |v: [u32; 3], n: [u32; 3], f: bool| Tri {
            v,
            n,
            fixed: [false, false, f],
            inside: old.inside,
            alive: true,
            stamp: 0,
        };
        self.set(t, make([a, b, v], [t1, t2, nc], fc));
        self.set(t1, make([b, c, v], [t2, t, na], fa));
        self.set(t2, make([c, a, v], [t, t1, nb], fb));
        self.relink(na, t, t1);
        self.relink(nb, t, t2);
        for s in [t, t1, t2] {
            self.touch(s, touched);
        }
        self.last = t;
        self.legalize(vec![(t, 2), (t1, 2), (t2, 2)], touched);
        v
    }

    fn split_edge(&mut self, t: u32, i: usize, p: Point2, touched: &mut Vec<u32>) -> u32 {
        let v = self.add_point(p);
        let tt = self.tris[t as usize].rotated(i);
        let u = tt.n[0];
        let [a, b, c] = tt.v;
        let uu = {
            let raw = self.tris[u as usize];
            let j = (0..3).find(|k| raw.n[*k] == t).unwrap();
            raw.rotated(j)
        };
        // uu = (d, c, b): its edge opposite d is (c, b).
        let d = uu.v[0];
        let t2 = self.new_tri(tt);
        let u2 = self.new_tri(uu);
        let make = |v: [u32; 3], n: [u32; 3], fixed: [bool; 3], inside: bool| Tri {
            v,
            n,
            fixed,
            inside,
            alive: true,
            stamp: 0,
        };
        // t: (a, b, v); t2: (a, v, c); u: (d, c, v); u2: (d, v, b).
        self.set(
            t,
            make(
                [a, b, v],
                [u2, t2, tt.n[2]],
                [false, false, tt.fixed[2]],
                tt.inside,
            ),
        );
        self.set(
            t2,
            make(
                [a, v, c],
                [u, tt.n[1], t],
                [false, tt.fixed[1], false],
                tt.inside,
            ),
        );
        self.set(
            u,
            make(
                [d, c, v],
                [t2, u2, uu.n[2]],
                [false, false, uu.fixed[2]],
                uu.inside,
            ),
        );
        self.set(
            u2,
            make(
                [d, v, b],
                [t, uu.n[1], u],
                [false, uu.fixed[1], false],
                uu.inside,
            ),
        );
        self.relink(tt.n[1], t, t2);
        self.relink(uu.n[1], u, u2);
        for s in [t, t2, u, u2] {
            self.touch(s, touched);
        }
        self.last = t;
        self.legalize(vec![(t, 2), (t2, 1), (u, 2), (u2, 1)], touched);
        v
    }

    /// Flip the edge opposite `v[i]` of `t`: the quad must be strictly
    /// convex. Returns the two slots, the old `t` now holding the new
    /// triangle on `v[i]`'s side of the new edge first.
    fn flip(&mut self, t: u32, i: usize, touched: &mut Vec<u32>) -> (u32, u32) {
        let tt = self.tris[t as usize].rotated(i);
        let u = tt.n[0];
        let uu = {
            let raw = self.tris[u as usize];
            let j = (0..3).find(|k| raw.n[*k] == t).unwrap();
            raw.rotated(j)
        };
        let [a, b, c] = tt.v;
        let d = uu.v[0];
        let make = |v: [u32; 3], n: [u32; 3], fixed: [bool; 3]| Tri {
            v,
            n,
            fixed,
            inside: tt.inside,
            alive: true,
            stamp: 0,
        };
        // uu = (d, c, b): (b, d) is opposite c (1), (d, c) opposite b (2).
        // t: (a, b, d); u: (a, d, c).
        self.set(
            t,
            make(
                [a, b, d],
                [uu.n[1], u, tt.n[2]],
                [uu.fixed[1], false, tt.fixed[2]],
            ),
        );
        self.set(
            u,
            make(
                [a, d, c],
                [uu.n[2], tt.n[1], t],
                [uu.fixed[2], tt.fixed[1], false],
            ),
        );
        self.relink(tt.n[1], t, u);
        self.relink(uu.n[1], u, t);
        self.touch(t, touched);
        self.touch(u, touched);
        (t, u)
    }

    /// Lawson flipping of the given edges and those it exposes.
    fn legalize(&mut self, mut stack: Vec<(u32, usize)>, touched: &mut Vec<u32>) {
        let mut budget = 64 * self.tris.len() + 1024;
        while let Some((t, i)) = stack.pop() {
            budget -= 1;
            if budget == 0 {
                return;
            }
            let tri = self.tris[t as usize];
            let u = tri.n[i];
            if !tri.alive || tri.fixed[i] || u == NONE {
                continue;
            }
            let uu = self.tris[u as usize];
            let Some(j) = (0..3).find(|k| uu.n[*k] == t) else {
                continue;
            };
            let d = uu.v[j];
            let [a, b, c] = tri.rotated(i).v;
            let p = |v: u32| self.points[v as usize];
            if !in_circle(p(a), p(b), p(c), p(d)) {
                continue;
            }
            // A certified violation: the quad is strictly convex.
            let (t, u) = self.flip(t, i, touched);
            // t = (a, b, d), u = (a, d, c): check the edges away from a.
            stack.push((t, 0));
            stack.push((u, 0));
        }
    }

    /// The triangle and index of the edge from `x` to `y` in counter-clockwise
    /// order (so the triangle lies left of `x -> y`).
    pub(super) fn directed_edge(&self, x: u32, y: u32) -> Option<(u32, usize)> {
        self.around(x).into_iter().find_map(|t| {
            let v = self.tris[t as usize].v;
            let k = v.iter().position(|w| *w == x)?;
            (v[(k + 1) % 3] == y).then_some((t, (k + 2) % 3))
        })
    }

    /// Every triangle incident to `x`.
    fn around(&self, x: u32) -> Vec<u32> {
        let start = self.star[x as usize];
        let mut out = Vec::new();
        if start == NONE {
            return out;
        }
        // Clockwise through the edge (x, next), then counter-clockwise from
        // the start when a hull edge stops the turn.
        let limit = self.tris.len();
        let mut t = start;
        for _ in 0..limit {
            out.push(t);
            let v = self.tris[t as usize].v;
            let k = v.iter().position(|w| *w == x).unwrap();
            let next = self.tris[t as usize].n[(k + 2) % 3];
            if next == NONE {
                break;
            }
            if next == start {
                return out;
            }
            t = next;
        }
        let mut t = start;
        for _ in 0..limit {
            let v = self.tris[t as usize].v;
            let k = v.iter().position(|w| *w == x).unwrap();
            let next = self.tris[t as usize].n[(k + 1) % 3];
            if next == NONE || next == start {
                return out;
            }
            out.push(next);
            t = next;
        }
        out
    }

    fn set_fixed(&mut self, t: u32, i: usize) {
        self.tris[t as usize].fixed[i] = true;
        let u = self.tris[t as usize].n[i];
        if u != NONE {
            let j = (0..3).find(|k| self.tris[u as usize].n[*k] == t).unwrap();
            self.tris[u as usize].fixed[j] = true;
        }
    }

    /// Insert the constraint `a`–`b`.
    pub(super) fn constrain(&mut self, a: u32, b: u32) -> Result<(), Conflict> {
        if a == b {
            return Err(Conflict);
        }
        let existing = self
            .directed_edge(a, b)
            .or_else(|| self.directed_edge(b, a));
        if let Some((t, i)) = existing {
            if self.tris[t as usize].fixed[i] {
                return Err(Conflict);
            }
            self.set_fixed(t, i);
            return Ok(());
        }
        let mut queue: VecDeque<(u32, u32)> = self.crossed(a, b)?.into();
        let p = |cdt: &Self, v: u32| cdt.points[v as usize];
        let crosses = |cdt: &Self, x: u32, y: u32| {
            let (pa, pb, px, py) = (p(cdt, a), p(cdt, b), p(cdt, x), p(cdt, y));
            orient(pa, pb, px) * orient(pa, pb, py) < 0
                && orient(px, py, pa) * orient(px, py, pb) < 0
        };
        let mut touched = Vec::new();
        let mut budget = 64 * (queue.len() + 16) * (queue.len() + 16);
        while let Some((x, y)) = queue.pop_front() {
            budget = budget.checked_sub(1).ok_or(Conflict)?;
            let (t, i) = self.directed_edge(x, y).ok_or(Conflict)?;
            let tri = self.tris[t as usize];
            let u = tri.n[i];
            if u == NONE {
                return Err(Conflict);
            }
            let w = tri.v[i];
            let j = (0..3)
                .find(|k| self.tris[u as usize].n[*k] == t)
                .ok_or(Conflict)?;
            let z = self.tris[u as usize].v[j];
            let (px, py, pw, pz) = (p(self, x), p(self, y), p(self, w), p(self, z));
            let convex = orient(pw, pz, px) * orient(pw, pz, py) < 0;
            if !convex {
                queue.push_back((x, y));
                continue;
            }
            self.flip(t, i, &mut touched);
            if crosses(self, w, z) {
                queue.push_back((w, z));
            }
        }
        let (t, i) = self.directed_edge(a, b).ok_or(Conflict)?;
        self.set_fixed(t, i);
        Ok(())
    }

    /// The edges the open segment `a`–`b` crosses, from `a`; a conflict when
    /// it passes through a vertex or crosses a constraint.
    fn crossed(&self, a: u32, b: u32) -> Result<Vec<(u32, u32)>, Conflict> {
        let (pa, pb) = (self.points[a as usize], self.points[b as usize]);
        let side = |v: u32| orient(pa, pb, self.points[v as usize]);
        let ahead = |v: u32| {
            let q = self.points[v as usize];
            (q.x - pa.x) * (pb.x - pa.x) + (q.y - pa.y) * (pb.y - pa.y) > 0.0
        };
        let mut current = None;
        for t in self.around(a) {
            let tri = self.tris[t as usize];
            let k = tri.v.iter().position(|w| *w == a).unwrap();
            let (x, y) = (tri.v[(k + 1) % 3], tri.v[(k + 2) % 3]);
            for w in [x, y] {
                if side(w) == 0 && ahead(w) {
                    return Err(Conflict);
                }
            }
            if side(x) < 0 && side(y) > 0 {
                current = Some((t, k));
                break;
            }
        }
        let (mut t, k) = current.ok_or(Conflict)?;
        // The crossed edge is opposite `k`: (right, left) of a -> b.
        let mut edge = k;
        let mut out = Vec::new();
        for _ in 0..self.tris.len() {
            let tri = self.tris[t as usize];
            let (x, y) = (tri.v[(edge + 1) % 3], tri.v[(edge + 2) % 3]);
            if tri.fixed[edge] {
                return Err(Conflict);
            }
            out.push((x, y));
            let u = tri.n[edge];
            if u == NONE {
                return Err(Conflict);
            }
            let uu = self.tris[u as usize];
            let j = (0..3).find(|m| uu.n[*m] == t).ok_or(Conflict)?;
            let w = uu.v[j];
            if w == b {
                return Ok(out);
            }
            // uu = (w, y, x) counter-clockwise from j.
            match side(w) {
                0 => return Err(Conflict),
                s if s < 0 => edge = (j + 2) % 3, // (w, y) next: opposite x
                _ => edge = (j + 1) % 3,          // (x, w) next: opposite y
            }
            t = u;
        }
        Err(Conflict)
    }

    /// Lawson flipping over every unconstrained edge.
    pub(super) fn make_delaunay(&mut self) {
        let mut stack = Vec::new();
        for t in 0..self.tris.len() as u32 {
            if self.tris[t as usize].alive {
                stack.extend([(t, 0), (t, 1), (t, 2)]);
            }
        }
        let mut touched = Vec::new();
        self.legalize(stack, &mut touched);
    }

    /// Mark the triangles inside the constraints: an odd number of
    /// constraints crossed from the enclosing triangle's corners. `Err` when
    /// a corner's triangle would be inside.
    pub(super) fn classify(&mut self) -> Result<(), Conflict> {
        let mut depth = vec![u32::MAX; self.tris.len()];
        let mut queue = VecDeque::new();
        for (t, tri) in self.tris.iter().enumerate() {
            if tri.alive && tri.v.iter().any(|v| *v < CORNERS) {
                depth[t] = 0;
                queue.push_back(t as u32);
            }
        }
        while let Some(t) = queue.pop_front() {
            let tri = self.tris[t as usize];
            for i in 0..3 {
                let u = tri.n[i];
                if u == NONE {
                    continue;
                }
                let d = depth[t as usize] + u32::from(tri.fixed[i]);
                if d < depth[u as usize] {
                    depth[u as usize] = d;
                    if tri.fixed[i] {
                        queue.push_back(u);
                    } else {
                        queue.push_front(u);
                    }
                }
            }
        }
        for (t, tri) in self.tris.iter_mut().enumerate() {
            tri.inside = tri.alive && depth[t] % 2 == 1;
            if tri.inside && tri.v.iter().any(|v| *v < CORNERS) {
                return Err(Conflict);
            }
        }
        Ok(())
    }

    /// Insert a point strictly inside the domain, starting the search at
    /// `hint`: `None` (nothing changed) unless it falls inside an inside
    /// triangle or on an unconstrained edge between two. Returns the vertex
    /// and every triangle created or changed.
    pub(super) fn insert_inside(&mut self, p: Point2, hint: u32) -> Option<(u32, Vec<u32>)> {
        if !(p.x.is_finite() && p.y.is_finite()) {
            return None;
        }
        let mut touched = Vec::new();
        let v = match self.locate(p, hint) {
            Location::In(t) if self.tris[t as usize].inside => {
                self.split_triangle(t, p, &mut touched)
            }
            Location::Edge(t, i) => {
                let tri = self.tris[t as usize];
                let u = tri.n[i];
                if tri.fixed[i] || !tri.inside || u == NONE || !self.tris[u as usize].inside {
                    return None;
                }
                self.split_edge(t, i, p, &mut touched)
            }
            _ => return None,
        };
        Some((v, touched))
    }
}
