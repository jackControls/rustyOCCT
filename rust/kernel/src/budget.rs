//! Work budgets, deadlines and cancellation (S10, `BUDGETS.md`).
//!
//! An operation runs under a [`Budget`]: a number of work units, an
//! optional wall-clock deadline and an optional [`Cancel`] token. The
//! kernel counts one unit at each of its loop points (an arrangement's
//! step, a verification box, a refinement piece, a certified integral's
//! piece, a chain attempt, a Newton step, ...) and checks the budget there.
//! A budget that runs out makes the operation return
//! [`Error::BudgetExhausted`] with an [`Exhaustion`]: what ran out, the
//! point it ran out at and the units counted; never a partial body.
//!
//! Work units are deterministic: the same operation on the same inputs
//! counts the same units on every host, in every run and whatever the
//! process computed before (a value the kernel keeps from an earlier
//! operation is charged the units its computation counted). A refusal for
//! work is therefore part of the contract: with a budget of `n` units an
//! operation that counts `w` succeeds, identically, when `w <= n`, and is
//! refused with exactly `n + 1` units counted when `w > n`. A deadline or a
//! cancellation is not deterministic and says so ([`Cause::WallTime`],
//! [`Cause::Cancelled`]).
//!
//! Without a budget nothing is counted and nothing is refused: every
//! operation keeps its result.
use crate::{Error, Result};
use std::cell::{Cell, RefCell};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, OnceLock};
use std::time::{Duration, Instant};

/// A cancellation token: cloned to another thread, [`Cancel::cancel`]
/// there makes every operation running under a budget holding it refuse at
/// its next work unit.
#[derive(Clone, Debug, Default)]
pub struct Cancel(Arc<AtomicBool>);

impl Cancel {
    pub fn new() -> Self {
        Self::default()
    }
    /// Requests cancellation; operations see it at their next work unit.
    pub fn cancel(&self) {
        self.0.store(true, Ordering::Relaxed);
    }
    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }
}

/// What ran out.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Cause {
    /// The work units: deterministic, the same on every host and run.
    Work,
    /// The deadline passed: depends on the host and its load.
    WallTime,
    /// The budget's [`Cancel`] token was cancelled.
    Cancelled,
}

/// The partial diagnosis of an operation refused by its budget.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Exhaustion {
    pub cause: Cause,
    /// The work-unit point that refused (`BUDGETS.md`'s table).
    pub stage: &'static str,
    /// The units the budget had counted, the refused one included: for
    /// [`Cause::Work`] always the budget's units plus one.
    pub work: u64,
}

impl Exhaustion {
    /// Whether the refusal is deterministic (work units), the same on every
    /// host and run.
    pub fn deterministic(&self) -> bool {
        self.cause == Cause::Work
    }
}

impl std::fmt::Display for Exhaustion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let cause = match self.cause {
            Cause::Work => "work units",
            Cause::WallTime => "wall time",
            Cause::Cancelled => "cancelled",
        };
        write!(
            f,
            "{cause} at {} after {} work units",
            self.stage, self.work
        )
    }
}

#[derive(Debug, Default)]
struct State {
    spent: AtomicU64,
    exhausted: OnceLock<Exhaustion>,
}

/// A work budget: units, an optional deadline and an optional cancellation
/// token. Clones share one count: every operation run under a budget or a
/// clone of it spends from the same units, and once exhausted the budget
/// refuses every later operation. It is `Send` and `Sync`; operations on
/// several threads may share one, though only one operation at a time
/// counts deterministically.
#[derive(Clone, Debug)]
pub struct Budget {
    units: u64,
    deadline: Option<Instant>,
    cancel: Option<Cancel>,
    state: Arc<State>,
}

impl Default for Budget {
    fn default() -> Self {
        Self::unlimited()
    }
}

impl Budget {
    /// Counts work and never refuses for it (a deadline or a token may
    /// still be added).
    pub fn unlimited() -> Self {
        Self::work(u64::MAX)
    }
    /// At most `units` work units (`u64::MAX`: unlimited).
    pub fn work(units: u64) -> Self {
        Self {
            units,
            deadline: None,
            cancel: None,
            state: Arc::default(),
        }
    }
    /// Refuses at the first work unit at or after `deadline`.
    pub fn with_deadline(self, deadline: Instant) -> Self {
        Self {
            deadline: Some(deadline),
            ..self
        }
    }
    /// A deadline `wall` from now.
    pub fn with_wall_time(self, wall: Duration) -> Self {
        let deadline = Instant::now().checked_add(wall);
        Self { deadline, ..self }
    }
    /// Refuses at the first work unit after `cancel` is cancelled.
    pub fn with_cancel(self, cancel: &Cancel) -> Self {
        Self {
            cancel: Some(cancel.clone()),
            ..self
        }
    }
    /// The units this budget allows.
    pub fn units(&self) -> u64 {
        self.units
    }
    /// The units counted so far by every operation run under this budget.
    pub fn spent(&self) -> u64 {
        self.state.spent.load(Ordering::Relaxed)
    }
    /// Why this budget refused, once it has.
    pub fn exhausted(&self) -> Option<Exhaustion> {
        self.state.exhausted.get().copied()
    }

    /// Runs `f` with this budget for every kernel operation it calls on
    /// this thread (operations without a [`crate::Context`], such as
    /// validation, tessellation, file import and intersections, take their
    /// budget this way). When the budget refuses, `f`'s result, whatever it
    /// is, becomes [`Error::BudgetExhausted`]: never a partial result.
    pub fn run<T, E: From<Error>>(
        &self,
        f: impl FnOnce() -> std::result::Result<T, E>,
    ) -> std::result::Result<T, E> {
        if let Err(e) = self.check("an operation's start") {
            return Err(e.into());
        }
        let out = {
            let _installed = Installed::new(self.clone());
            f()
        };
        match self.exhausted() {
            Some(e) => Err(Error::BudgetExhausted(e).into()),
            None => out,
        }
    }

    /// Refuses when exhausted, cancelled or past the deadline, counting
    /// nothing.
    fn check(&self, stage: &'static str) -> Result<()> {
        self.charge(stage, 0)
    }

    fn refuse(&self, cause: Cause, stage: &'static str, work: u64) -> Error {
        let e = self
            .state
            .exhausted
            .get_or_init(|| Exhaustion { cause, stage, work });
        Error::BudgetExhausted(*e)
    }

    fn charge(&self, stage: &'static str, units: u64) -> Result<()> {
        if let Some(e) = self.state.exhausted.get() {
            return Err(Error::BudgetExhausted(*e));
        }
        if self.cancel.as_ref().is_some_and(Cancel::is_cancelled) {
            return Err(self.refuse(Cause::Cancelled, stage, self.spent()));
        }
        let after = self
            .state
            .spent
            .fetch_add(units, Ordering::Relaxed)
            .saturating_add(units);
        if after > self.units {
            // The refused unit is counted, not done: `units + 1` in all.
            let work = self.units + 1;
            self.state.spent.store(work, Ordering::Relaxed);
            return Err(self.refuse(Cause::Work, stage, work));
        }
        if self.deadline.is_some_and(|d| Instant::now() >= d) {
            return Err(self.refuse(Cause::WallTime, stage, after));
        }
        Ok(())
    }
}

/// A run of one point's units in a recorded computation.
type Run = (&'static str, u64);

thread_local! {
    /// Whether a unit does anything: a budget installed or a computation
    /// being recorded.
    static ARMED: Cell<bool> = const { Cell::new(false) };
    /// The budgets installed on this thread, the innermost last.
    static INSTALLED: RefCell<Vec<Budget>> = const { RefCell::new(Vec::new()) };
    /// The recordings of kept computations, the innermost last.
    static RECORDING: RefCell<Vec<Vec<Run>>> = const { RefCell::new(Vec::new()) };
}

fn rearm() {
    let installed = INSTALLED.with(|s| !s.borrow().is_empty());
    let recording = RECORDING.with(|r| !r.borrow().is_empty());
    ARMED.with(|a| a.set(installed || recording));
}

struct Installed;

impl Installed {
    fn new(budget: Budget) -> Self {
        INSTALLED.with(|s| s.borrow_mut().push(budget));
        rearm();
        Installed
    }
}

impl Drop for Installed {
    fn drop(&mut self) {
        INSTALLED.with(|s| s.borrow_mut().pop());
        rearm();
    }
}

/// Runs `f` under `budget` when there is one (an operation with a
/// [`crate::Context`]): see [`Budget::run`]. Without one, an operation
/// called inside another keeps the outer operation's budget.
pub(crate) fn within<T>(budget: Option<&Budget>, f: impl FnOnce() -> Result<T>) -> Result<T> {
    match budget {
        Some(b) => b.run(f),
        None => f(),
    }
}

/// One work unit at the point `stage`: refused when the budget is spent,
/// past its deadline or cancelled. Without a budget it costs a
/// thread-local read.
#[inline]
pub(crate) fn unit(stage: &'static str) -> Result<()> {
    if ARMED.with(Cell::get) {
        charge(stage, 1)
    } else {
        Ok(())
    }
}

/// [`unit`] where a refusal cannot be returned at once: `false` when the
/// budget refused (the operation's result becomes the refusal at its end,
/// whatever the caller does next).
#[inline]
pub(crate) fn unit_ok(stage: &'static str) -> bool {
    unit(stage).is_ok()
}

#[inline(never)]
fn charge(stage: &'static str, units: u64) -> Result<()> {
    RECORDING.with(|r| {
        if let Some(top) = r.borrow_mut().last_mut() {
            match top.last_mut() {
                Some((s, n)) if std::ptr::eq(*s, stage) => *n += units,
                _ => top.push((stage, units)),
            }
        }
    });
    INSTALLED.with(|s| match s.borrow().last() {
        Some(b) => b.charge(stage, units),
        None => Ok(()),
    })
}

/// Whether the budget installed on this thread, if any, has refused.
pub(crate) fn refused() -> bool {
    ARMED.with(Cell::get)
        && INSTALLED.with(|s| s.borrow().last().is_some_and(|b| b.exhausted().is_some()))
}

/// The units a kept computation counted, charged again on each later use
/// so that counts do not depend on what was kept.
#[derive(Clone, Debug, Default)]
pub(crate) struct Recorded(Arc<[Run]>);

impl Recorded {
    /// Charges the recorded units again (a kept value's use), refusing
    /// exactly where the computation would have.
    pub(crate) fn replay(&self) -> Result<()> {
        if !ARMED.with(Cell::get) {
            return Ok(());
        }
        for &(stage, units) in self.0.iter() {
            charge(stage, units)?;
        }
        Ok(())
    }
}

/// Runs a computation whose value is kept for later operations, recording
/// its units. The value may be kept only when the record is `Some`: `None`
/// when the budget refused during it (its value may be that refusal or a
/// consequence of it).
pub(crate) fn recorded<T>(f: impl FnOnce() -> T) -> (T, Option<Recorded>) {
    struct Recording;
    impl Drop for Recording {
        fn drop(&mut self) {
            let runs = RECORDING.with(|r| r.borrow_mut().pop()).unwrap_or_default();
            // An enclosing recording counts these units too.
            RECORDING.with(|r| {
                if let Some(top) = r.borrow_mut().last_mut() {
                    for (stage, units) in runs {
                        match top.last_mut() {
                            Some((s, n)) if std::ptr::eq(*s, stage) => *n += units,
                            _ => top.push((stage, units)),
                        }
                    }
                }
            });
            rearm();
        }
    }
    RECORDING.with(|r| r.borrow_mut().push(Vec::new()));
    rearm();
    let guard = Recording;
    let value = f();
    let runs = RECORDING.with(|r| r.borrow().last().cloned().unwrap_or_default());
    drop(guard);
    let clean = !refused();
    (value, clean.then(|| Recorded(runs.into())))
}
