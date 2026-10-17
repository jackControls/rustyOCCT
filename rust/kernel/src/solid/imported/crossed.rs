//! S9e.4b.4b.2b.2b.3b: a prism leaf whose caps another leaf crossing it
//! splits (REVIEW_NOTES.md, "S9e.4b.4b.2b.2b.3b refined").
//!
//! A stadium plate with a slab through it on crossing axes stores each of
//! the plate's caps as several faces on one plane, each bounded by the
//! plate's walls and by the slab's faces (its runs along the slab), and the
//! plate's walls' joints inside the slab not stored at all. Such a leaf is
//! read off all the faces of its bottom cap's plane at once
//! (`crossed_leaves`): every face of the group one loop whose every edge is a
//! wall along the cap's inward normal or a face of the lead (a run along
//! it), at least one run each, its top the group of faces facing up that its
//! walls reach. Its profile (`crossed_prism`) is the group's loops with each
//! run dropped and the walls at its ends extended through the lead until
//! they meet the next piece's walls: two lines at their common point, a line
//! and a circle (or two circles) at their meeting nearest the run within
//! the lead, tangent ones at the foot of the circle's centre (the joints of
//! a stadium's line and arc, hidden in the lead), two walls on one line or
//! circle one segment; each end's continuation unique, its join and its
//! continuations' middles inside the lead at the leaf's bottom, middle and
//! top, and the pieces one closed loop, else no such leaf. Rounded once into
//! the bottom face's stored frame (turned to the cap's inward normal, or
//! against it where its cylinder walls' stored axes point so, as a stack's
//! leaf), the joints in binary64, the arcs' ends taken onto their circles by
//! the exact model (S9e.4b.1); the leaf is fused after the lead.
use super::*;

/// A loop's piece: its segments' start points, the segments and its end.
#[derive(Debug, Clone)]
struct Chain {
    points: Vec<Point2>,
    segments: Vec<Segment>,
    end: Point2,
}

/// A segment's carrier: a line through a point along a direction, or a
/// circle turning one way.
#[derive(Debug, Clone, Copy)]
enum Carrier {
    Line { p: Point2, d: (f64, f64) },
    Circle { c: Point2, r: f64, ccw: bool },
}

/// How a piece's end continues into the next piece's start: through a
/// joint, or along one carrier (the two segments one).
#[derive(Debug, Clone, Copy, PartialEq)]
enum Join {
    At(Point2),
    Along,
}

fn carrier(segment: &Segment, a: Point2, b: Point2) -> Carrier {
    match segment {
        Segment::Arc {
            center,
            radius,
            ccw,
        } => Carrier::Circle {
            c: *center,
            r: *radius,
            ccw: *ccw,
        },
        _ => Carrier::Line {
            p: a,
            d: (b.x - a.x, b.y - a.y),
        },
    }
}

fn sub(a: Point2, b: Point2) -> (f64, f64) {
    (a.x - b.x, a.y - b.y)
}

fn dot(a: (f64, f64), b: (f64, f64)) -> f64 {
    a.0 * b.0 + a.1 * b.1
}

fn cross(a: (f64, f64), b: (f64, f64)) -> f64 {
    a.0 * b.1 - a.1 * b.0
}

/// `sqrt` of a sum of squares, correctly rounded on every host (`hypot` is
/// not).
fn norm(a: (f64, f64)) -> f64 {
    (a.0 * a.0 + a.1 * a.1).sqrt()
}

/// The turn from `a` to `b` about a circle's centre in its sense, in
/// `(-pi, pi]`.
fn turn(c: Point2, ccw: bool, a: Point2, b: Point2) -> f64 {
    let (u, v) = (sub(a, c), sub(b, c));
    let t = cross(u, v).atan2(dot(u, v));
    if ccw {
        t
    } else {
        -t
    }
}

/// The meetings of two carriers: a line's and a circle's, or two circles',
/// tangent within `tol` at one point.
fn meetings(a: Carrier, b: Carrier, tol: f64) -> Vec<Point2> {
    match (a, b) {
        (Carrier::Line { p, d }, Carrier::Line { p: q, d: e }) => {
            let den = cross(d, e);
            if den.abs() <= 1e-12 * norm(d) * norm(e) {
                return Vec::new();
            }
            let s = cross(sub(q, p), e) / den;
            vec![Point2::new(p.x + s * d.0, p.y + s * d.1)]
        }
        (Carrier::Line { p, d }, Carrier::Circle { c, r, .. })
        | (Carrier::Circle { c, r, .. }, Carrier::Line { p, d }) => {
            let l = norm(d);
            let u = (d.0 / l, d.1 / l);
            let s = dot(sub(c, p), u);
            let f = Point2::new(p.x + s * u.0, p.y + s * u.1);
            let h = norm(sub(c, f));
            if (h - r).abs() <= tol {
                vec![f]
            } else if h < r {
                let k = (r * r - h * h).sqrt();
                vec![
                    Point2::new(f.x - k * u.0, f.y - k * u.1),
                    Point2::new(f.x + k * u.0, f.y + k * u.1),
                ]
            } else {
                Vec::new()
            }
        }
        (Carrier::Circle { c, r, .. }, Carrier::Circle { c: e, r: s, .. }) => {
            let w = sub(e, c);
            let d = norm(w);
            if d <= tol {
                return Vec::new();
            }
            let u = (w.0 / d, w.1 / d);
            // The radical line's foot along the centres' line.
            let a = (d * d + r * r - s * s) / (2.0 * d);
            let f = Point2::new(c.x + a * u.0, c.y + a * u.1);
            let h2 = r * r - a * a;
            if (d - (r + s)).abs() <= tol || (d - (r - s).abs()).abs() <= tol {
                vec![f]
            } else if h2 > 0.0 {
                let k = h2.sqrt();
                vec![
                    Point2::new(f.x - k * u.1, f.y + k * u.0),
                    Point2::new(f.x + k * u.1, f.y - k * u.0),
                ]
            } else {
                Vec::new()
            }
        }
    }
}

/// Whether two carriers are one: lines through one another's points along
/// one direction, circles of one centre, radius and sense, within `tol`.
fn one_carrier(a: Carrier, b: Carrier, tol: f64) -> bool {
    match (a, b) {
        (Carrier::Line { p, d }, Carrier::Line { p: q, d: e }) => {
            let l = norm(d);
            dot(d, e) > 0.0
                && cross(d, e).abs() <= 1e-12 * l * norm(e)
                && (cross(sub(q, p), d) / l).abs() <= tol
        }
        (Carrier::Circle { c, r, ccw }, Carrier::Circle { c: e, r: s, ccw: w }) => {
            ccw == w && norm(sub(c, e)) <= tol && (r - s).abs() <= tol
        }
        _ => false,
    }
}

/// Whether `m` lies on from `from` along a carrier (a line's direction, a
/// circle's sense, at most half a turn), within `tol`.
fn ahead(k: Carrier, from: Point2, m: Point2, tol: f64) -> bool {
    match k {
        Carrier::Line { d, .. } => dot(sub(m, from), d) / norm(d) >= -tol,
        Carrier::Circle { c, r, ccw } => {
            let t = turn(c, ccw, from, m);
            t >= -tol / r && t < std::f64::consts::PI
        }
    }
}

/// The point midway along a carrier from `a` to `b` (a line's middle, an
/// arc's middle turning in its sense).
fn middle(k: Carrier, a: Point2, b: Point2) -> Point2 {
    match k {
        Carrier::Line { .. } => Point2::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y)),
        Carrier::Circle { c, r, ccw } => {
            let t0 = (a.y - c.y).atan2(a.x - c.x);
            let t = t0 + 0.5 * if ccw { 1.0 } else { -1.0 } * turn(c, ccw, a, b);
            Point2::new(c.x + r * t.cos(), c.y + r * t.sin())
        }
    }
}

/// The crossed chains a body is tried with: each leaf of `leaves` (each
/// face set once) leading the leaves it crosses, where any.
pub(super) fn chains(t: &Topology, tolerance: Tolerance) -> Vec<Vec<Leaf>> {
    let mut seen: Vec<Vec<bool>> = Vec::new();
    let mut out = Vec::new();
    for lead in leaves(t, tolerance) {
        if seen.contains(&lead.faces) {
            continue;
        }
        seen.push(lead.faces.clone());
        let crossed = crossed_leaves(t, tolerance, &lead);
        if !crossed.is_empty() {
            let mut chain = vec![lead];
            chain.extend(crossed);
            out.push(chain);
        }
    }
    out
}

/// The leaves of several faces on one plane `lead` crosses
/// (S9e.4b.4b.2b.2b.3b), in face order: for each group of two plane faces
/// or more on one plane facing one way (`one_side`), none the lead's and its
/// normal not along the lead's, every face one loop each of whose edges'
/// other face is a lead face (a run, at least one a face) or a wall along the
/// group's inward normal (no lead face) rising from it, the walls all
/// meeting one group of faces facing the other way (the top). Its faces: the
/// groups, the walls and every face on their surfaces.
pub(super) fn crossed_leaves(t: &Topology, tolerance: Tolerance, lead: &Leaf) -> Vec<Leaf> {
    let tol = tolerance.linear();
    let n = t.faces().len();
    let fins = fin_faces(t);
    let plane = |k: usize| matches!(t.faces()[k].surface, Surface::Plane(_));
    let across = |face: usize| -> Option<Vec<usize>> {
        let [l] = t.faces()[face].loops.as_slice() else {
            return None;
        };
        let Loop::Edges { fins: lf, .. } = &t.loops()[l.index()] else {
            return None;
        };
        lf.iter()
            .map(|f| {
                let e = &t.edges()[t.fins()[f.index()].edge.index()];
                e.fins.iter().find(|x| *x != f).map(|x| fins[x.index()])
            })
            .collect()
    };
    let mut out = Vec::new();
    let mut seen = vec![false; n];
    for b in 0..n {
        if seen[b] || lead.faces[b] || !plane(b) {
            continue;
        }
        let Some(m) = outward(t, b) else { continue };
        if parallel(m, lead.caps.2) {
            continue;
        }
        let up = -m;
        let group: Vec<usize> = (0..n)
            .filter(|&k| !lead.faces[k] && plane(k) && (k == b || one_side(t, b, k, tol)))
            .collect();
        for &k in &group {
            seen[k] = true;
        }
        if group.len() < 2 {
            continue;
        }
        let mut walls = Vec::new();
        let mut fit = true;
        for &k in &group {
            let Some(ws) = across(k) else {
                fit = false;
                break;
            };
            let runs = ws.iter().filter(|&&w| w < n && lead.faces[w]).count();
            fit &= runs > 0
                && ws.iter().all(|&w| {
                    w < n
                        && (lead.faces[w]
                            || (wall(t, w, up) && rises(t, k, w, up, tol) && !group.contains(&w)))
                });
            walls.extend(ws.into_iter().filter(|&w| w < n && !lead.faces[w]));
        }
        walls.sort_unstable();
        walls.dedup();
        if !fit || walls.is_empty() {
            continue;
        }
        // The top: the faces facing up the walls meet, one group.
        let mut tops = Vec::new();
        for e in t.edges() {
            let [x, y] = e.fins.as_slice() else { continue };
            let (x, y) = (fins[x.index()], fins[y.index()]);
            for (w, c) in [(x, y), (y, x)] {
                if walls.contains(&w)
                    && !lead.faces[c]
                    && plane(c)
                    && outward(t, c).is_some_and(|mc| parallel(mc, up) && mc.dot(up) > 0.0)
                {
                    tops.push(c);
                }
            }
        }
        tops.sort_unstable();
        tops.dedup();
        let Some(&top) = tops.first() else { continue };
        if !tops.iter().all(|&c| c == top || one_side(t, top, c, tol)) {
            continue;
        }
        let mut faces = vec![false; n];
        for &k in group.iter().chain(&tops).chain(&walls) {
            faces[k] = true;
        }
        let members: Vec<usize> = (0..n).filter(|&k| faces[k]).collect();
        for (k, face) in faces.iter_mut().enumerate() {
            if !*face && !lead.faces[k] && members.iter().any(|&x| one_side(t, x, k, tol)) {
                *face = true;
            }
        }
        // Read off either cap, one leaf.
        if out.iter().any(|o: &Leaf| o.faces == faces) {
            continue;
        }
        out.push(Leaf {
            caps: (b, top, up),
            loops: vec![0],
            faces,
            open: None,
            run: Vec::new(),
            void: false,
            pieces: group,
        });
    }
    out
}

/// A crossed leaf's prism (`crossed_leaves`), its pieces' runs along `lead`
/// bridged through it; refused `not_their_chain` where they do not close.
pub(super) fn crossed_prism(
    t: &Topology,
    leaf: &Leaf,
    lead: &Solid,
    tolerance: Tolerance,
    op: OperationId,
) -> Result<Solid> {
    let tol = tolerance.linear();
    let faces = t.faces();
    let fins_of = fin_faces(t);
    let (bottom, top, up) = leaf.caps;
    let Surface::Plane(stored) = &faces[bottom].surface else {
        return Err(general());
    };
    let Surface::Plane(top_frame) = &faces[top].surface else {
        return Err(general());
    };
    // Its pieces' loops: points, segments and runs (as `completed_on`).
    let mut pieces: Vec<(Vec<Point2>, Vec<Segment>, Vec<bool>)> = Vec::new();
    // On the bottom's stored frame turned to `up` (or against it where its
    // cylinder walls' stored axes point so).
    let mut against = false;
    for &b in &leaf.pieces {
        let [l] = faces[b].loops.as_slice() else {
            return Err(not_their_chain());
        };
        let Loop::Edges { fins, .. } = &t.loops()[l.index()] else {
            return Err(not_their_chain());
        };
        for f in fins {
            let edge = &t.edges()[t.fins()[f.index()].edge.index()];
            against |= edge.fins.iter().filter(|x| *x != f).any(|x| {
                matches!(&faces[fins_of[x.index()]].surface,
                    Surface::Cylinder { frame, .. } if frame.normal().dot(up) < 0.0)
            });
        }
    }
    let dir = if against { -up } else { up };
    let frame = if stored.normal().dot(dir) > 0.0 {
        *stored
    } else {
        stored.flipped()
    };
    let local = Local::new(&frame)?;
    let h = rational_f64(&local.of(top_frame.origin())[2]);
    let (lo, hi) = (h.min(0.0), h.max(0.0));
    if hi - lo <= 0.0 {
        return Err(general());
    }
    for &b in &leaf.pieces {
        let [l] = faces[b].loops.as_slice() else {
            return Err(not_their_chain());
        };
        let Loop::Edges { fins, .. } = &t.loops()[l.index()] else {
            return Err(not_their_chain());
        };
        let (mut points, mut segments, mut runs) = (Vec::new(), Vec::new(), Vec::new());
        for f in fins {
            let fin = &t.fins()[f.index()];
            let edge = &t.edges()[fin.edge.index()];
            if edge.is_ring() {
                return Err(not_their_chain());
            }
            let forward = fin.sense == Orientation::Forward;
            let across = edge
                .fins
                .iter()
                .find(|x| *x != f)
                .map(|x| fins_of[x.index()]);
            runs.push(across.is_some_and(|k| k < faces.len() && !leaf.faces[k]));
            let start = if forward { edge.start } else { edge.end }.ok_or_else(general)?;
            points.push(point2(&local.of(t.vertices()[start.index()].position)));
            segments.push(match &edge.curve {
                Curve3::LineSegment { .. } => Segment::Line,
                Curve3::CircularArc {
                    frame: c,
                    radius,
                    sweep_angle,
                    ..
                } => {
                    let turns = (*sweep_angle > 0.0) == forward;
                    let about = c.normal().dot(dir) > 0.0;
                    Segment::Arc {
                        center: point2(&local.of(c.origin())),
                        radius: *radius,
                        ccw: turns == about,
                    }
                }
                _ => return Err(general()),
            });
        }
        pieces.push((points, segments, runs));
    }
    // Each piece's chains between its runs.
    let mut chains: Vec<Chain> = Vec::new();
    for (points, segments, runs) in &pieces {
        let n = points.len();
        let Some(first) = (0..n).find(|&j| runs[j] && !runs[(j + 1) % n]) else {
            return Err(not_their_chain());
        };
        let mut j = (first + 1) % n;
        let mut current: Option<Chain> = None;
        for _ in 0..n {
            if runs[j] {
                if let Some(mut c) = current.take() {
                    c.end = points[j];
                    chains.push(c);
                }
            } else {
                let c = current.get_or_insert_with(|| Chain {
                    points: Vec::new(),
                    segments: Vec::new(),
                    end: points[j],
                });
                c.points.push(points[j]);
                c.segments.push(segments[j].clone());
            }
            j = (j + 1) % n;
        }
        if current.is_some() {
            return Err(not_their_chain());
        }
    }
    if chains.len() < 2 || chains.len() > 16 {
        return Err(not_their_chain());
    }
    // Inside the lead at the leaf's bottom, middle and top.
    let inside = |p: Point2| {
        [lo, 0.5 * (lo + hi), hi].iter().all(|&w| {
            matches!(
                lead.classify(frame.point(p, w)),
                Ok(crate::Location::Inside)
            )
        })
    };
    let last = |c: &Chain| -> Carrier {
        let k = c.points.len() - 1;
        carrier(&c.segments[k], c.points[k], c.end)
    };
    let first = |c: &Chain| -> Carrier {
        let b = c.points.get(1).copied().unwrap_or(c.end);
        carrier(&c.segments[0], c.points[0], b)
    };
    // Each chain's end continued into one chain's start.
    let mut next: Vec<(usize, Join)> = Vec::new();
    for a in &chains {
        let (ka, pa) = (last(a), a.end);
        let mut found: Vec<(usize, Join)> = Vec::new();
        for (bi, b) in chains.iter().enumerate() {
            let (kb, pb) = (first(b), b.points[0]);
            if one_carrier(ka, kb, tol) {
                if ahead(ka, pa, pb, tol) && inside(middle(ka, pa, pb)) {
                    found.push((bi, Join::Along));
                }
                continue;
            }
            for m in meetings(ka, kb, tol) {
                if ahead(ka, pa, m, tol)
                    && ahead(kb, m, pb, tol)
                    && inside(m)
                    && inside(middle(ka, pa, m))
                    && inside(middle(kb, m, pb))
                {
                    found.push((bi, Join::At(m)));
                }
            }
        }
        let [only] = found.as_slice() else {
            return Err(not_their_chain());
        };
        next.push(*only);
    }
    // One closed loop through every chain.
    let mut order = vec![0usize];
    while order.len() < chains.len() {
        let k = next[*order.last().expect("a chain")].0;
        if order.contains(&k) {
            return Err(not_their_chain());
        }
        order.push(k);
    }
    if next[*order.last().expect("a chain")].0 != 0 {
        return Err(not_their_chain());
    }
    let incoming = |c: usize| next.iter().find(|x| x.0 == c).map(|x| x.1);
    let mut points = Vec::new();
    let mut segments = Vec::new();
    for &c in &order {
        let chain = &chains[c];
        for (i, (p, s)) in chain.points.iter().zip(&chain.segments).enumerate() {
            match (i, incoming(c)) {
                (0, Some(Join::Along)) => continue,
                (0, Some(Join::At(m))) => points.push(m),
                _ => points.push(*p),
            }
            segments.push(s.clone());
        }
    }
    if points.len() < 2 {
        return Err(not_their_chain());
    }
    let outer = if segments.iter().all(|s| matches!(s, Segment::Line)) {
        Boundary::polygon(points, tolerance)?
    } else {
        Boundary::path(points, segments, tolerance)?
    };
    let profile = Profile::new(outer, Vec::new(), tolerance)?.with_rounded_arcs();
    Solid::build(op, profile, frame, lo, hi)
}
