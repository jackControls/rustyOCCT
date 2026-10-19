//! S10's invariant (`BUDGETS.md`, `FUZZING.md`): an operation run under a
//! budget either succeeds, identically to its unlimited run, with its work
//! within the budget, or is refused for work having counted exactly one
//! unit past it; and the budget fixtures' work counts.
use rusty_occt::budget::{Cause, Exhaustion};
use rusty_occt::{Budget, Error};
use std::cell::{Cell, RefCell};
use std::fmt::Debug;

thread_local! {
    /// The work counted by every metered operation of the current input, in
    /// order (`recorded`).
    static WORK: RefCell<Option<Vec<(&'static str, u64)>>> = const { RefCell::new(None) };
    /// The current input's hash and the metered operations so far.
    static INPUT: Cell<(u64, u64)> = const { Cell::new((0, 0)) };
}

/// Runs `check` and returns the work counted by each metered operation it
/// ran, in order (`replay_budgets`, the budget fixtures).
pub fn recorded(check: impl FnOnce()) -> Vec<(&'static str, u64)> {
    WORK.with(|w| *w.borrow_mut() = Some(Vec::new()));
    check();
    WORK.with(|w| w.borrow_mut().take()).unwrap_or_default()
}

/// Starts an input: the budgets its metered operations draw (FNV-1a of its
/// bytes).
pub(crate) fn begin(data: &[u8]) {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for &x in data {
        h = (h ^ u64::from(x)).wrapping_mul(0x0000_0100_0000_01b3);
    }
    INPUT.with(|i| i.set((h ^ (h >> 29), 0)));
}

/// The next draw of the current input.
fn draw() -> u64 {
    INPUT.with(|i| {
        let (h, n) = i.get();
        i.set((h, n + 1));
        let mut x = h ^ (n + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        x = (x ^ (x >> 31)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x ^ (x >> 27)
    })
}

/// Whether this is the input's first metered operation.
pub(crate) fn first() -> bool {
    INPUT.with(|i| i.get().1 == 0)
}

/// Which of `n` operations the current input checks under a budget.
pub(crate) fn pick(n: u64) -> u64 {
    INPUT.with(|i| i.get().0) % n
}

/// A budget below `2^24` drawn log-uniformly (0 included), before the
/// operation's work is known.
fn units(h: u64) -> u64 {
    let k = h % 25;
    if k == 0 {
        0
    } else {
        (1 << (k - 1)) + (h >> 8) % (1 << (k - 1))
    }
}

/// `op` under an unlimited budget, its result returned and its work
/// recorded; when `checked`, first under a budget drawn from the input,
/// checking S10's invariant against the unlimited run:
/// - work `w` within the budget `n`: the same result, `w` units counted;
/// - `w > n`: refused for work at `stage`, exactly `n + 1` units counted.
///
/// The limited run comes first, so a refusal is met inside computations
/// the kernel keeps (and must then not keep), and a success leaves them
/// kept for the unlimited run to charge again.
pub(crate) fn metered<T: PartialEq + Debug>(
    label: &'static str,
    checked: bool,
    op: impl Fn(&Budget) -> Result<T, Error>,
) -> Result<T, Error> {
    let limited = checked.then(|| {
        let budget = Budget::work(units(draw()));
        let out = op(&budget);
        (budget, out)
    });
    let full = Budget::unlimited();
    let out = op(&full);
    let w = full.spent();
    assert!(full.exhausted().is_none(), "an unlimited budget refused");
    if let Some((budget, limited)) = limited {
        let n = budget.units();
        if w <= n {
            assert_eq!(budget.spent(), w, "{label}: the work under a budget of {n}");
            assert!(
                budget.exhausted().is_none(),
                "{label}: refused within its budget"
            );
            assert!(
                limited == out,
                "{label}: {limited:?} under a budget of {n}, {out:?} unlimited ({w} units)"
            );
        } else {
            assert_eq!(
                budget.spent(),
                n + 1,
                "{label}: a refusal's work under {n} ({w} unlimited)"
            );
            match limited {
                Err(Error::BudgetExhausted(Exhaustion {
                    cause: Cause::Work,
                    work,
                    ..
                })) if work == n + 1 => {}
                other => panic!("{label}: {other:?} under a budget of {n} units, {w} unlimited"),
            }
        }
    }
    WORK.with(|r| {
        if let Some(r) = r.borrow_mut().as_mut() {
            r.push((label, w));
        }
    });
    out
}
