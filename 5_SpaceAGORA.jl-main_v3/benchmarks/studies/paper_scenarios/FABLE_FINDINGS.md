# Fable Findings — Performance Review of the Parallelization Hot Paths

Findings from a full read of the parallelization-relevant source (July 2026,
by Claude Fable during the paper_scenarios suite build). Part 1 lists the fixes
that were implemented and verified alongside this suite; Part 2 documents the
deferred opportunities — each with location, mechanism, and a suggested approach —
for future reference. Line numbers are as of the commit introducing this file.

---

## Part 1 — Implemented fixes

### F1. Flat-queue partials reduction was satellite-major (cache-hostile)
`src/simulation/engine/dynamics_rhs.jl`, reduction after the flat effector queue.
`partials` is column-major `6×N×W`; the old loop nest iterated satellites outer,
workers inner, so each inner iteration jumped `6N` doubles. Swapped to
worker-major/satellite-inner: unit-stride over both `partials` and `totals`.
At 1024 sat × 8 workers this is ~400 KB of traffic per RHS evaluation made
sequential instead of scattered.

### F2. Flat-queue scratch zeroing now uses contiguous `fill!` when exact-size
Same file, `_ensure_rhs_flat_effector_scratch!`. The steady-state case (buffer
dimensions exactly match) now zeroes with `fill!` (memset) instead of a strided
view broadcast.

### F3. Work-item sort comparator re-derived costs O(n log n) times
Same file, `_prepare_rhs_flat_work_items!`. Item cost depends only on `eff_idx`
(≤ a handful of distinct values), but the `by=` comparator re-ran the cost-model
lookup per comparison. Costs are now resolved once into an `ntuple` and the
comparator indexes it.

### F4. Robot-arm effector scan ran per satellite per RHS call
`_apply_coupled_robot_arm_rhs!` linearly scanned the control + dynamic effector
tuples with `hasproperty` reflection for every satellite on every RHS evaluation —
in simulations with no robot arm anywhere (i.e., every constellation/MC case in
this suite). Added `SharedBuffers.robot_arm_present`, a lazily-computed
run-constant flag (`src/core/types/runtime_types.jl`), checked first.

### F5. Process-route Monte Carlo re-serialized the sample closure per sample
`src/simulation/campaigns/monte_carlo.jl`, `_run_monte_carlo_process`. Each
`remotecall_fetch` shipped `f` — which can capture a large configuration — anew.
Now dispatches a wrapping closure through a `Distributed.CachingPool`, so the
closure crosses the wire once per worker and is referenced by identity afterwards
(`clear!` drops it after the campaign).

### F6. Redundant `sort!` in both MC runners
`samples[index] = sample` plus the in-order filter already produce index-sorted
output; the trailing `sort!(completed; by=s->s.index)` was a no-op cost. Removed
in the threaded and process runners.

Verification: package precompiles clean; scenario_worker smoke runs pass for the
constellation (vacuum + GRAM look-ahead), MC-threads, and MC-process paths;
before/after timings on the LEO constellation cases are recorded in the suite
results. (See also the repo test suites run in the same session.)

---

## Part 2 — Deferred findings (documented, not implemented)

### D1. Flat work-item list is rebuilt every RHS call
`dynamics_rhs.jl`, `_prepare_rhs_flat_work_items!` — O(N·E) reconstruction (plus
conditional sort) on every RHS evaluation, though the active-satellite set and
effector tuple change rarely (deactivation events only). **Approach:** cache the
packed item list in `SharedBuffers` behind an active-set generation counter
(bump on any `is_active` mutation) and a cost-model revision; rebuild only on
invalidation. Interacts with F3: once cached, the sort cost disappears entirely
between invalidations.

### D2. `threaded_reduce` has no persistent-pool variant and allocates per call
`src/parallel/policy/thread_execution.jl`. The per-satellite effector reduction
(`:per_satellite_effector_reduce` route) allocates a `Vector` of partial
accumulators and pays raw `Threads.@spawn` latency on every invocation, inside
the innermost RHS. `threaded_foreach` got persistent pools and a spin-barrier
variant; the reduce path never did. **Approach:** mirror
`threaded_foreach_worker_persistent` with per-worker pre-allocated accumulator
slots keyed by (scope, source), combining serially after the barrier.

### D3. Persistent-pool dispatch allocates a NamedTuple per worker per call
Same file, `_threaded_foreach_persistent!` — each call `put!`s a fresh 6-field
NamedTuple into each worker's request channel under `run_lock` (~1–5 µs per
dispatch). The spin-barrier path (~10–50 ns) exists but only for the harmonics
batch. **Approach:** generalize the spin-barrier opt-in to `:rhs_flat_queue` and
`:density_callback` sources; the pool machinery is already source-keyed.

### D4. Density callback does full kinematics for out-of-atmosphere satellites
`src/simulation/callbacks/density_callbacks/runtime.jl`. The `DiscreteCallback`
fires every accepted step (`condition = true`) and `_stage_environment_kinematics`
(planet-frame transform included) runs for every satellite before any
in-atmosphere gate. For eccentric/aerobraking-style orbits with long coast arcs
most of that work computes densities nothing reads. Also, `affect!` re-derives
the thread decision and batch-model selection every step; both are run-constant.
**Approach:** early-skip satellites flagged out-of-atmosphere (needs care: heat
-rate paths read the buffered density, so the skip must leave a defined stale
value and the `density_sample_t` stamp unset), and cache the decision/batch
selection in `SharedBuffers` at run start. Not implemented now because the
correctness surface (buffered-density consumers) needs a dedicated parity pass.

### D5. GRAM track-cache refreshes convoy on the global GRAM lock
`src/simulation/callbacks/gram_track_cache/refresh.jl` — each satellite's refresh
holds `GRAM_LOCK` for a full native `generate_trajectory` call. Constellation
members launched with similar epochs expire their caches near-simultaneously, so
integration threads serialize behind the lock in bursts. This is the "GRAM's own
global lock is the next ceiling" item confirmed by the constellation
thread-imbalance study. **Approaches** (in increasing ambition): stagger cache
horizons per satellite (e.g., ±10% jitter) so expiries de-phase; batch several
satellites' track requests into one locked native call; prefetch the next
segment asynchronously before expiry so hot-path queries never block.

### D6. Ensemble members each pay a full configuration `deepcopy`
`src/simulation/campaigns/constellation_ensemble.jl`. Root cause: solves mutate
shared model state (e.g., `planet.L_PI` in the planet-frame callback), so
concurrent members must isolate the whole configuration. **Approach
(architectural):** move mutable frame state out of model objects into per-run
buffers (`SharedBuffers` already exists for exactly this); then members could
share one immutable configuration and the O(members × config) copy cost — memory
and time — disappears. Would also remove `run_simulation`'s own `isolate_state`
copy for single runs.

### D7. Clenshaw harmonics evaluator exists but is undecided
Root-level `test_clenshaw_evaluator.jl` benchmarks a `calcForceTorque_clenshaw`
variant against the AssociatedLegendrePolynomials-based implementation. Clenshaw
summation avoids materializing the full Legendre matrix (better cache behavior at
L50) and would remove the mutable per-model workspace state that currently forces
per-worker pool copies in the SIMD batch kernel. **Action:** decide land-or-drop
before the paper freezes harmonics numbers; if landed, re-run S1/S5.

### D8. `rhs_kind` symbol comparisons inside the flat per-satellite closure
`dynamics_rhs.jl`, `_spacecraft_dynamics_flat_constellation_effector_queue!` —
`rhs_kind == :implicit` etc. branch at runtime inside the per-satellite loop.
Minor. **Approach:** `Val{rhs_kind}` specialization or hoisting the branch into
three loop bodies, if profiling ever shows it matters.

### D9. Per-step policy re-decision in `_prefill_environment_samples!`
`dynamics_rhs.jl` calls `_density_callback_thread_decision` on every flat RHS
call with atmosphere prefill. Run-constant inputs; cacheable alongside D4's
decision caching.
