//! S8c.1: a cone, frustum or sphere zone split by a plane normal to its
//! axis, and a whole sphere by any plane.
//!
//! In the solid's frame the plane is `F = a u + b v + c w + d`. Which side
//! the solid lies on is exact: `F`'s extremes over a cone are over its end
//! circles (`c w + d +- r |(a, b)|`, an apex at `c w + d`), over a zone at
//! the sphere's extreme points in the directions `+-(a, b, c)` when they lie
//! between its ends and otherwise over its end circles, each a quadratic
//! surd's sign. A plane normal to the axis (`a = b = 0`) cuts at `w = -d/c`,
//! rounded: the pieces are the same primitive between exact heights (a
//! cone's radius there rounded, a zone's latitude `asin(w / R)` rounded).
//! A whole sphere has every axis: its pieces are two caps on a frame whose
//! axis is the plane's normal, cut at the plane's signed distance from the
//! centre. Other planes cut a cone in a conic and a zone in a circle whose
//! pcurves are transcendental graphs over the angle, which wait for D13's
//! procedural edges (S8d).
//!
//! Names follow provenance: the wall and the region are `Split` into one
//! child per piece (below first), an end the piece keeps is the input's
//! (`Unchanged`, or `Modified` when its stored geometry moved by a
//! rounding), the cut disc and its ring are `Generated` from the wall, and a
//! whole sphere's new poles are `Generated` from it with the role `Pole`.
use super::{q, rational_f64, zero, Side};
use crate::certified::{Interval as I, Real};
use crate::history::{History, Relation};
use crate::identity::{
    Derivation, EntityId, EntityKind, OperationKind, Parent, ProfileElement, Role,
};
use crate::solid::{Construction, Context, Solid};
use crate::{Error, Frame3, Point2, Result};
use num_rational::BigRational as R;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::f64::consts::FRAC_PI_2;

/// The sign of `k + t sqrt(r2)`, `r2 >= 0`, `t` in `{1, -1}`, exactly.
fn surd_sign(k: &R, r2: &R, t: i8) -> Ordering {
    let kk = k * k;
    match (k.cmp(&zero()), t) {
        (Ordering::Equal, _) => {
            if *r2 == zero() {
                Ordering::Equal
            } else if t > 0 {
                Ordering::Greater
            } else {
                Ordering::Less
            }
        }
        (Ordering::Greater, 1) => Ordering::Greater,
        (Ordering::Less, -1) => Ordering::Less,
        // k - sqrt(r2), k > 0.
        (Ordering::Greater, _) => kk.cmp(r2),
        // sqrt(r2) - |k|, k < 0.
        (_, _) => r2.cmp(&kk),
    }
}

/// The pieces with their sides, and the history.
type Pieces = (Vec<(Side, Solid)>, History);

/// How an end of a piece came about.
#[derive(Clone, Copy, PartialEq, Eq)]
enum End {
    /// The input's end at this place (its roles match the piece's).
    Kept,
    /// The cut: a disc and its ring.
    Cut,
    /// A new pole (a whole sphere's cap).
    Pole,
}

/// Whether a derivation's apex or pole is at the top (meridian vertex 2).
fn at_top_vertex(d: &Derivation) -> bool {
    matches!(
        d.parents.first(),
        Some(Parent::Profile {
            element: ProfileElement::Vertex(2),
            ..
        })
    )
}

impl Solid {
    /// S8c.1's split of a cone, frustum or sphere zone; `None` when the
    /// solid is not one of them.
    pub(super) fn split_revolved(
        &self,
        context: &Context,
        plane: &Frame3,
        [a, b, c, d]: [R; 4],
    ) -> Result<Option<Pieces>> {
        let (tolerance, sphere) = match &self.construction {
            Construction::Cone { tolerance, .. } => (*tolerance, None),
            Construction::Sphere {
                radius,
                low,
                high,
                tolerance,
            } => (*tolerance, Some((*radius, *low, *high))),
            _ => return Ok(None),
        };
        let tol = tolerance.linear();
        let ab2 = &a * &a + &b * &b;
        let (w0, w1) = (q(self.start), q(self.end));
        // Signs of F at its extremes over the solid.
        let mut signs: Vec<Ordering> = Vec::new();
        let circle = |w: &R, r2: &R, signs: &mut Vec<Ordering>| {
            let k = &c * w + &d;
            for t in [1i8, -1] {
                signs.push(surd_sign(&k, &(r2 * &ab2), t));
            }
        };
        match (&self.construction, sphere) {
            (Construction::Cone { bottom, top, .. }, _) => {
                circle(&w0, &(q(*bottom) * q(*bottom)), &mut signs);
                circle(&w1, &(q(*top) * q(*top)), &mut signs);
            }
            (_, Some((radius, low, high))) => {
                let r2 = q(radius) * q(radius);
                // The end circles (a pole's radius zero).
                for (w, pole) in [(&w0, low == -FRAC_PI_2), (&w1, high == FRAC_PI_2)] {
                    let rho2 = if pole { zero() } else { &r2 - w * w };
                    circle(w, &rho2, &mut signs);
                }
                // The extreme points `t R (a, b, c) / |m|` (F = d + t R |m|)
                // when their heights `t R c / |m|` lie between the ends.
                let m2 = &ab2 + &c * &c;
                for t in [1i8, -1] {
                    let k = q(radius) * &c * R::from_integer(t.into());
                    // The sign of `t R c - end |m|`.
                    let against = |end: &R| {
                        surd_sign(&k, &(end * end * &m2), if *end >= zero() { -1 } else { 1 })
                    };
                    if against(&w0) != Ordering::Less && against(&w1) != Ordering::Greater {
                        signs.push(surd_sign(&d, &(&r2 * &m2), t));
                    }
                }
            }
            _ => unreachable!("a cone or a sphere"),
        }
        let above = signs.contains(&Ordering::Greater);
        let below = signs.contains(&Ordering::Less);
        if !(above && below) {
            return Ok(Some(self.unchanged(context)));
        }
        let operation = context.operation;
        let n = self.frame.normal();
        // Normal to the axis to within a quarter of the resolution over the
        // solid's widest circle (a frame whose stored axes are not exactly
        // orthogonal): the flat cut at the axis's crossing lies within it
        // of the plane everywhere, so the pieces are those of the normal
        // plane.
        let widest = match (&self.construction, sphere) {
            (Construction::Cone { bottom, top, .. }, _) => bottom.max(*top),
            (_, Some((radius, ..))) => radius,
            _ => unreachable!("a cone or a sphere"),
        };
        let normal = ab2 == zero()
            || &ab2 * q(widest) * q(widest) * R::from_integer(16.into())
                <= q(tol) * q(tol) * &c * &c;
        // The pieces with their sides and their ends' provenance (bottom,
        // top).
        let mut pieces: Vec<(Side, Solid, [End; 2])> = Vec::new();
        if normal {
            // Normal to the axis: the cut at w = -d / c.
            let h = -&d / &c;
            let height = rational_f64(&h);
            let (lo, hi) = (self.start.min(self.end), self.start.max(self.end));
            if !(lo < height && height < hi) {
                return Err(Error::Degenerate("a split within binary64 of a cap"));
            }
            if height - lo <= tol || hi - height <= tol {
                return Err(Error::Degenerate("a piece thinner than the resolution"));
            }
            let (lower, upper) = match (&self.construction, sphere) {
                (Construction::Cone { bottom, top, .. }, _) => {
                    let r = q(*bottom) + (q(*top) - q(*bottom)) * q(height) / q(self.end);
                    let r = rational_f64(&r);
                    if r <= tol {
                        return Err(Error::Degenerate("a cut within the resolution of the apex"));
                    }
                    let lifted = Frame3::new(
                        self.frame.point(Point2::default(), height),
                        n,
                        self.frame.x(),
                        tolerance,
                    )?;
                    (
                        Solid::build_cone(operation, self.frame, *bottom, r, height, tolerance)?,
                        Solid::build_cone(
                            operation,
                            lifted,
                            r,
                            *top,
                            self.end - height,
                            tolerance,
                        )?,
                    )
                }
                (_, Some((radius, low, high))) => {
                    let latitude = (height / radius).asin();
                    if !(low < latitude && latitude < high) {
                        return Err(Error::Degenerate("a split within binary64 of a cap"));
                    }
                    (
                        Solid::build_sphere(
                            operation, self.frame, radius, low, latitude, tolerance,
                        )?,
                        Solid::build_sphere(
                            operation, self.frame, radius, latitude, high, tolerance,
                        )?,
                    )
                }
                _ => unreachable!("a cone or a sphere"),
            };
            // F = c (w - h): the lower piece lies below the plane when c > 0.
            let (sl, su) = if c > zero() {
                (Side::Below, Side::Above)
            } else {
                (Side::Above, Side::Below)
            };
            // A whole sphere has no pole entities: its caps' poles are new.
            let end = match sphere {
                Some((_, low, high)) if low == -FRAC_PI_2 && high == FRAC_PI_2 => End::Pole,
                _ => End::Kept,
            };
            pieces.push((sl, lower, [end, End::Cut]));
            pieces.push((su, upper, [End::Cut, end]));
        } else if let Some((radius, low, high)) = sphere {
            if !(low == -FRAC_PI_2 && high == FRAC_PI_2) {
                return Err(Error::OutOfDomain(
                    "a sphere zone by a plane not normal to its axis (S8d)",
                ));
            }
            // Caps on a frame along the plane's normal, cut at the plane's
            // signed distance from the centre: F(centre) = d, so the plane
            // lies at -d / |m| along it.
            let m = plane.normal();
            let hint = if m.cross(self.frame.x()).length() > 0.5 {
                self.frame.x()
            } else {
                self.frame.y()
            };
            let frame = Frame3::new(self.frame.origin(), m, hint, tolerance)?;
            let m2 = &ab2 + &c * &c;
            let distance = I::exact(-d.clone())
                .div(&I::exact(m2).sqrt())
                .ok_or(Error::Degenerate("a plane's normal"))?;
            let (lo, hi) = distance.bounds_f64();
            let height = 0.5 * lo + 0.5 * hi;
            if radius - height.abs() <= tol {
                return Err(Error::Degenerate("a cap thinner than the resolution"));
            }
            let latitude = (height / radius).asin();
            // The frame's normal is the plane's: the cap under the latitude
            // lies below it.
            pieces.push((
                Side::Below,
                Solid::build_sphere(operation, frame, radius, -FRAC_PI_2, latitude, tolerance)?,
                [End::Pole, End::Cut],
            ));
            pieces.push((
                Side::Above,
                Solid::build_sphere(operation, frame, radius, latitude, FRAC_PI_2, tolerance)?,
                [End::Cut, End::Pole],
            ));
        } else {
            return Err(Error::OutOfDomain(
                "a cone by a plane not normal to its axis (S8d)",
            ));
        }
        pieces.sort_by_key(|p| p.0);
        self.rename_revolved(context, pieces).map(Some)
    }

    /// Name the pieces by provenance and report the history.
    fn rename_revolved(
        &self,
        context: &Context,
        mut pieces: Vec<(Side, Solid, [End; 2])>,
    ) -> Result<Pieces> {
        let operation = context.operation;
        let parent = &self.topology;
        let body = parent.body_id();
        let role_of = |id: EntityId| parent.derivation(id).map(|d| d.role);
        let find = |role: Role| {
            parent
                .ids()
                .map(|(id, _)| id)
                .find(|id| role_of(*id) == Some(role))
        };
        let missing = || Error::InvalidTopology("a revolved solid without its wall or region");
        let (wall, region) = (
            find(Role::Wall).ok_or_else(missing)?,
            find(Role::Region).ok_or_else(missing)?,
        );
        // An input end entity by role (an apex or pole by its place).
        let kept = |role: Role, top: bool| -> Option<EntityId> {
            parent.ids().map(|(id, _)| id).find(|id| {
                let d = parent.derivation(*id).expect("a derivation");
                d.role == role
                    && (!matches!(role, Role::Apex | Role::Pole) || at_top_vertex(d) == top)
            })
        };
        let mut relations = Vec::new();
        let mut kept_ids: BTreeSet<EntityId> = BTreeSet::new();
        let mut split: Vec<(EntityId, Vec<EntityId>)> =
            vec![(wall, Vec::new()), (region, Vec::new())];
        for (k, (_, solid, ends)) in pieces.iter_mut().enumerate() {
            let t = &solid.topology;
            let mut derivations = Vec::new();
            for (id, slot) in t.ids().collect::<Vec<_>>() {
                let d = t.derivation(id).expect("a derivation");
                let top = match d.role {
                    Role::EndCap | Role::TopEdge => Some(true),
                    Role::StartCap | Role::BottomEdge => Some(false),
                    Role::Apex | Role::Pole => Some(at_top_vertex(d)),
                    _ => None,
                };
                let named = match (d.role, top) {
                    (Role::Wall | Role::Region, _) => {
                        let from = if d.role == Role::Wall { wall } else { region };
                        let child = Derivation {
                            operation,
                            kind: OperationKind::PlaneSplit,
                            entity: d.entity,
                            role: d.role,
                            ordinal: k as u32,
                            parents: vec![Parent::Entity(from)],
                        };
                        let list = &mut split.iter_mut().find(|s| s.0 == from).expect("a parent").1;
                        list.push(child.id());
                        child
                    }
                    (role, Some(top)) => match ends[usize::from(top)] {
                        End::Kept => {
                            let id = kept(role, top)
                                .ok_or(Error::InvalidTopology("a kept end without its entity"))?;
                            kept_ids.insert(id);
                            parent.derivation(id).expect("a derivation").clone()
                        }
                        End::Cut | End::Pole => {
                            let role = match (ends[usize::from(top)], d.entity) {
                                (End::Pole, _) => Role::Pole,
                                (_, EntityKind::Face) => Role::CutFace,
                                _ => Role::CutEdge,
                            };
                            let new = Derivation {
                                operation,
                                kind: OperationKind::PlaneSplit,
                                entity: d.entity,
                                role,
                                ordinal: k as u32,
                                parents: vec![Parent::Entity(wall)],
                            };
                            relations.push(Relation::Generated {
                                from: new.parents.clone(),
                                to: new.id(),
                                role,
                            });
                            new
                        }
                    },
                    _ => return Err(Error::InvalidTopology("a revolved slot of unknown role")),
                };
                derivations.push((slot, named));
            }
            let body_derivation = Derivation {
                operation,
                kind: OperationKind::PlaneSplit,
                entity: EntityKind::Body,
                role: Role::Body,
                ordinal: k as u32,
                parents: vec![Parent::Entity(body)],
            };
            solid.topology = solid
                .topology
                .clone()
                .renamed(body_derivation, derivations)?;
            solid.operation = operation;
        }
        for (from, into) in split {
            relations.push(Relation::Split { from, into });
        }
        let before = parent.entity_set(self.resolution());
        let mut modified = BTreeSet::new();
        for (_, piece, _) in &pieces {
            let after = piece.topology.entity_set(piece.resolution());
            for (id, info) in &after.entities {
                if kept_ids.contains(id) {
                    let old = &before.entities[id];
                    if old.geometry != info.geometry || old.structure != info.structure {
                        modified.insert(*id);
                    }
                }
            }
        }
        for id in &kept_ids {
            relations.push(if modified.contains(id) {
                Relation::Modified { from: *id, to: *id }
            } else {
                Relation::Unchanged { id: *id }
            });
        }
        relations.sort_by_cached_key(Relation::sort_key);
        relations.dedup();
        let mut pieces: Vec<(Side, Solid)> = pieces.into_iter().map(|(s, p, _)| (s, p)).collect();
        let outputs: Vec<EntityId> = pieces.iter().map(|(_, s)| s.topology.body_id()).collect();
        let history = History::new(
            operation,
            OperationKind::PlaneSplit,
            vec![body],
            outputs,
            relations,
            Vec::new(),
        );
        let mut history = history.at_level(context.level);
        {
            let mut outs: Vec<&mut Solid> = pieces.iter_mut().map(|(_, s)| s).collect();
            crate::solid::enclose::carry(&[self], &history, &mut outs);
        }
        let outs: Vec<&Solid> = pieces.iter().map(|(_, s)| s).collect();
        let maps = crate::solid::attrs::propagate(context, &[self], &mut history, &outs)?;
        for ((_, s), map) in pieces.iter_mut().zip(maps) {
            s.topology.set_attributes(map);
        }
        let outs: Vec<&Solid> = pieces.iter().map(|(_, s)| s).collect();
        crate::solid::stack::debug_check(&[self], &outs, &history);
        crate::solid::attrs::debug_check_attributes(context, &[self], &outs, &history);
        Ok((pieces, history))
    }
}
