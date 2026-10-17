//! S9e.4b.4c.2b.2: the inside-first order of a body of several primitives'
//! chains (REVIEW_NOTES.md, "S9e.4b.4c.2b.2 refined").
//!
//! A chain takes its prisms first (its leaves, then its prisms of one cap)
//! and then its primitives widest first, the first of them inside its
//! quadric. Two of the dataset's bodies are chains only in another order: a
//! rod whose end is a dimple of a ball wider than it (`OCC14536-sh1`: the
//! widest primitive, first, has its material outside its sphere) and a
//! flange whose spigot, a prism of one cap, stands in a recess its outer
//! cylinder cuts (`CTO904_cts21347a`: the spigot built before the recess is
//! cut, so the cut takes it away). Every chain is therefore tried once more
//! after all of them (`late`): its primitives inside their quadrics first,
//! widest first, then those outside, widest first (`order`), and its prisms
//! of one cap after the primitives (`sequence`). A body that imports keeps
//! the chain it imported by; any refusal of a later chain is one the body
//! is not.

/// A chain's primitives `(group, outside, reach)` in the inside-first order:
/// those whose material lies inside their quadric first, then those
/// outside, each widest first (stable: ties in their order).
pub(super) fn order(primitives: &mut [(usize, bool, f64)]) {
    primitives.sort_by(|a, b| a.1.cmp(&b.1).then(b.2.total_cmp(&a.2)));
}

/// The parts' places in the inside-first tree: of `parts` (the prisms, the
/// first `leaves` of them leaves and the rest up to `lead` prisms of one
/// cap, then the primitives), the leaves, then the primitives, then the
/// prisms of one cap.
pub(super) fn sequence(leaves: usize, lead: usize, parts: usize) -> Vec<usize> {
    (0..leaves).chain(lead..parts).chain(leaves..lead).collect()
}

/// A choice's chain: past `total` the chains again in the inside-first
/// order (`late`), the choice among them.
pub(super) fn split(choice: u32, total: u32) -> (bool, u32) {
    if choice >= total {
        (true, choice - total)
    } else {
        (false, choice)
    }
}

#[cfg(test)]
mod tests {
    use super::super::{primitives_choices, renamed, Recognized};
    use super::{order, sequence, split};
    use crate::identity::OperationId;
    use crate::solid::Construction;
    use crate::Solid;

    /// A fixture body's chain: its choice and the total of the choices
    /// before the inside-first ones.
    fn chain_of(text: &str) -> (u32, u32) {
        let doc = crate::occt_brep::read(text).unwrap();
        let [solid] = <[_; 1]>::try_from(crate::occt_brep::import(&doc).solids).unwrap();
        let (stored, resolution) = (solid.result.unwrap(), solid.tolerance);
        let (s, _) = Solid::imported_with(OperationId(11), stored.clone(), resolution).unwrap();
        let Construction::Imported(i) = &s.construction else {
            panic!("an import")
        };
        let Recognized::Piece(p) = &i.recognized else {
            panic!("a piece")
        };
        assert!(p.chain, "a chain");
        let topology = renamed(stored, OperationId(11)).unwrap();
        (p.choice, primitives_choices(&topology, resolution).3)
    }

    /// The earlier steps' bodies of several primitives keep their chains
    /// (before the inside-first ones); this step's take inside-first ones.
    #[test]
    fn earlier_bodies_keep_their_chains() {
        for text in [
            include_str!("../../../../fixtures/imported/chain_cup.brep"),
            include_str!("../../../../fixtures/imported/chain_shaft.brep"),
            include_str!("../../../../fixtures/imported/leaf_boss.brep"),
            include_str!("../../../../fixtures/imported/two_through.brep"),
            include_str!("../../../../fixtures/imported/rim_well.brep"),
            include_str!("../../../../fixtures/imported/flush_arms.brep"),
            include_str!("../../../../fixtures/imported/stack_step.brep"),
            include_str!("../../../../fixtures/imported/nested_post.brep"),
            include_str!("../../../../fixtures/imported/crossed_joint.brep"),
        ] {
            let (choice, total) = chain_of(text);
            assert!(choice < total, "{choice} of {total}");
        }
        for text in [
            include_str!("../../../../fixtures/imported/inside_dimple.brep"),
            include_str!("../../../../fixtures/imported/inside_flange.brep"),
        ] {
            let (choice, total) = chain_of(text);
            assert!(choice >= total, "{choice} of {total}");
        }
    }

    #[test]
    fn inside_first_places_the_cuts_after_and_the_bosses_last() {
        let mut p = vec![
            (0, true, 5.0),
            (1, false, 3.0),
            (2, false, 4.0),
            (3, true, 6.0),
        ];
        order(&mut p);
        assert_eq!(p.iter().map(|x| x.0).collect::<Vec<_>>(), [2, 1, 3, 0]);
        assert_eq!(sequence(1, 3, 5), [0, 3, 4, 1, 2]);
        assert_eq!(sequence(0, 0, 2), [0, 1]);
        assert_eq!(split(3, 7), (false, 3));
        assert_eq!(split(9, 7), (true, 2));
    }
}
