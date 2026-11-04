# WS11f: three small, already-scoped performance wins

Branch `pv2c-ws11f-small-wins`, off `adb343566` on `policy-v2-consistency`.
The hard requirement for every item below: bit-identical dynamics. See
`test/README.md` and the reference-case protocol this document follows for
what "bit-identical" and "the reference case set" mean.

All timings quoted are ratios between two runs of the same script, on the
same machine (`space-falcon-1`), in one process state, taken back to back
while the machine was shared with several other benchmark jobs (multiple
WS11 workstreams run concurrently by design -- see
`benchmarks/studies/paper_scenarios/README.md`'s serialization rule and the
WS11 common contract). They are not absolute benchmarks.

## Item 1: `@simd ivdep` on the harmonics batch kernel

**Status: already shipped, re-verified, one stale comment fixed.**

The kernel this item targets, `_harmonics_flat_batch_kernel!`, lives in
`src/dynamics/coupled/perturbations.jl`, not under
`src/environment/gravity/` as the WS11f contract's file list says (that list
predates a refactor, or was written against an older layout; the small
`GravityEffectors`/`gravity_models.jl` files under `src/environment/gravity/`
only carry `ConstantGravityModel`/`InverseSquaredGravityModel`/
`InverseSquaredJ2GravityModel`/`GravityGradientTorqueModel`, none of which
have a batch kernel). This document and the commit that follows treat
ownership as scoped to the harmonics batch kernel functions in
`perturbations.jl`, the same way the WS11f contract scopes
`effector_sampling.jl` ownership to one function rather than the whole file.

`@simd ivdep` on the kernel's three per-degree batch loops (`b = 1:B`) was
already applied and shipped at this branch's base commit, in
`cd212833aa` ("The batched harmonics loops vectorise at every batch size"),
itself preceded by the rolling-window change (`904979aec5`) the earlier memory
note ("Harmonics register blocking (negative)") called "the 2-line win."
Both commits already reported bit-identical results against the scalar
kernel. There was no code left to apply.

What this session did instead:

1. Re-verified bit-identity at this exact tip with a fresh dump/cmp, since
   the contract's proof protocol asks for evidence at the tip, not just a
   citation of an earlier commit message.
2. Fixed a stale doc comment on `_harmonics_flat_batch_kernel!` that still
   read "No `@turbo`, `@fastmath` or `@simd`" -- true when it was written,
   false since `cd212833aa`, and actively misleading to a future reader
   deciding whether the kernel's bit-identity claim still holds. The comment
   now names `@simd ivdep` explicitly and explains why it does not license
   reassociation here (each `b` iteration writes only its own workspace slot;
   there is no loop-carried dependency and no reduction across `b` for the
   annotation to reassociate).

No other change. `git diff` on `src/dynamics/coupled/perturbations.jl` is
comment-only.

### Proof: dump/cmp, this tip, `@simd` present vs. temporarily removed

Tool: `benchmarks/studies/third_body_cost/variants.jl --dump`, 1 thread,
`--variants=vacuum --mission=5800`. "Before" = the three
`@inbounds @simd ivdep for b = 1:B` sites temporarily reverted to
`@inbounds for b = 1:B` (not committed, produced only to regenerate the
proof); "after" = the committed tip.

| Case | N | Before bytes | After bytes | `cmp` |
|---|---|---|---|---|
| `gravity_256sat_l50_vacuum_5800s` (contract's 256-spacecraft vacuum case) | 256 | 4868424 | 4868424 | BYTE-IDENTICAL |
| `gravity_4096sat_l50_vacuum_5800s` (the P2 case) | 4096 | 77859144 | 77859144 | BYTE-IDENTICAL |

Terminal state of spacecraft 1 also agreed to all 17 printed digits in both
cases. `test/unit/simulation/harmonics_batch_parity_tests.jl` (22 checks)
passes unchanged.

### Ratio, 1 and 8 threads (best of repeats; machine shared throughout)

| Case | Threads | Mode | Before | After | Ratio (before/after) |
|---|---|---|---|---|---|
| N=256 vacuum | 1 | serial | 0.768 s | 0.596 s | 1.29x |
| N=4096 vacuum (P2) | 1 | serial | 15.396 s | 14.216 s | 1.08x |
| N=4096 vacuum (P2) | 8 | inner_only (R2) | 4.201 s | 3.487 s | 1.21x |

Consistent in direction and rough magnitude with `cd212833aa`'s own
measurement ("median 1.03x, best 1.32x" across 19 launch points). Raw numbers:
`benchmarks/studies/small_wins/results/harmonics_simd_ratio.csv`. Reproduction
script (the "after" half; see its header for the "before" half, which needs a
temporary source edit): `benchmarks/studies/small_wins/harmonics_simd_dump.jl`.

New regression guard: `test/unit/dynamics/harmonics_kernel_identity_tests.jl`
-- a denser batch-size sweep than the existing parity test (every `B` from 1
to 20, plus 31/32/33, 63/64/65, 127/128/129, straddling the 2/4/8-double SIMD
lane widths a `@simd ivdep` vectorization pass could pick), and a full
end-to-end trajectory comparison between the flat (batch/SIMD) and
per-satellite (scalar) RHS routes via `SPACEAGORA_RHS_EXECUTION_MODE`, not
just the raw kernel call.

## Item 2: the third-body sample's remaining allocation

**Status: fixed, values unchanged, allocation reduced from 4272 B/call to
192 B/call (22.3x).**

`sample_third_body_ephemerides` (`src/simulation/engine/effector_sampling.jl`)
still allocated after WS10c's `map`-over-`body_names` fix
(`docs/architecture/third_body_cost.md`, "Deliverable 1"). Root cause, found
with `Profile.Allocs` (not guessed): the `do`-block closure passed to `map`
captured two things it did not need to.

1. `p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls` was
   read from *inside* the closure. Every other shared-buffer value the
   closure needs (`et`, `primary_body_name`, `spice_rhs_memo_enabled`,
   `spice_rhs_memo`, `cache_entry`) was already hoisted to a local before the
   closure; this one field access was not, so the closure captured the whole
   `p::ODEParams{...}` -- an inlinable immutable holding the entire
   `SimulationConfiguration` -- just to reach one `Threads.Atomic{Int64}`
   three levels down. `Profile.Allocs` attributed a fresh 2032-byte
   `ODEParams{...}`-typed allocation to this closure's construction, once per
   call: the compiler had kept `p` split into scalar SSA values earlier in
   the function and had to re-materialise it whole to build the closure.
2. `cache_entry::Union{Nothing, NBodyEphemerisCache}` was captured directly.
   The closure's own type is therefore a `Union` of two closure types (one
   per branch of the captured field), and `map` over a `Union`-typed functor
   cannot infer a concrete element type: `positions_ii` came back
   `Tuple{Any, Any}`, boxing both `SVector{3, Float64}` results (confirmed by
   `@code_warntype`).

Fix, in two steps, each independently measured with `Profile.Allocs` plus
`@allocated`:

1. Hoist the counter to a local (`nbody_spkpos_runtime_calls =
   p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls`)
   before the closure, exactly like the other shared-buffer reads. This alone
   removed the 2032 B `ODEParams` allocation and the closure struct's own
   allocation shrank with it (the closure no longer needs to embed an
   inlined copy of `p`): 4272 B -> 304 B/call.
2. Take the `cache_entry isa NBodyEphemerisCache` branch *outside* the `map`,
   once, instead of inside the shared closure -- two `map` calls, each
   capturing a concretely-typed `cache_entry` (or none at all). Both bodies
   are the pre-existing branch bodies, unchanged; nothing about the values,
   order, or cache-then-backend fallback moved. This let `map`'s result infer
   as `NTuple{N, SVector{3, Float64}}` in both branches instead of
   `Tuple{Any, Any}`: 304 B -> 192 B/call. The residual 192 B/call is two
   boxed `SVector{3, Float64}` results plus their tuple, which persisted
   across an `if`/`else` merge of two otherwise-fully-typed `map` calls; not
   chased further; restructuring the sum/recursion to remove it was
   explicitly out of scope for this item and risks the bit-identity
   guarantee for no measured benefit at this point (192 B/call is already
   below the noise floor of a single derivative evaluation's other
   allocations at any constellation size this repo measures at).

| Stage | Bytes/call | Source |
|---|---|---|
| WS10c tip (`adb343566`) | 4272.0 | `Profile.Allocs`/`@allocated`, this session |
| + hoist the counter local | 304.0 | `Profile.Allocs`/`@allocated`, this session |
| + split the `cache_entry` branch outside `map` (shipped) | 192.0 | `Profile.Allocs`/`@allocated`, this session |

Raw numbers: `benchmarks/studies/small_wins/results/third_body_sample_allocation.csv`.
`@allocated` test (asserts < 512 B/call, well above the measured 192 but far
below the 4272 pre-fix floor so a regression trips it):
`benchmarks/studies/small_wins/third_body_sample_allocation.jl`.

### Proof: dump/cmp, the SRP + third-body reference case

Tool: `benchmarks/studies/third_body_cost/variants.jl --dump`, 1 thread,
`--n=256 --mission=5800 --variants=srp_nbody`.

| Before bytes | After bytes | `cmp` |
|---|---|---|
| 4868424 | 4868424 | BYTE-IDENTICAL |

Terminal state of spacecraft 1 agreed to all 17 printed digits.
`test/unit/simulation/third_body_route_parity_tests.jl` (84 checks) and
`test/unit/dynamics/third_body_cost_tests.jl` (732 checks) pass unchanged.

## Item 3: precompile the campaign path

**Status: extended, measured, tests pass.**

`_warm_campaign_dispatchers` (`SimulationCampaigns`, `monte_carlo.jl`)
exercises the low-level Monte Carlo dispatchers directly -- the serial loop,
the mixed dispatcher on a trivial local-slots call, the threaded dispatcher
-- but never the campaign-level *planning* layer above them, because it
calls those dispatchers straight rather than through `run_monte_carlo`. That
layer is two largely separate code paths: `_campaign_route_plan` /
`_run_campaign_with_route_env` for the default bandit route, and
`predictive_plan` / `_run_campaign_predictive` / `_PredictiveGuardState` for
the R7 planner (`SPACEAGORA_CAMPAIGN_PLANNER=predictive`). Neither compiles
into the pkgimage today, so an interactive session's first campaign (of
either kind) pays to JIT it.

`src/precompile_workload.jl` now runs two more tiny campaigns through
`run_monte_carlo`'s public entry inside the package's `@compile_workload`
block, both with two trivial `seed -> seed * 2` samples:

- `_warm_predictive_campaign()`: `SPACEAGORA_CAMPAIGN_PLANNER=predictive`,
  `route_tuning=OuterRouteTuning(process_max_workers=1)`. "One worker means
  no pool" is the same pattern
  `test/unit/parallel/predictive_campaign_tests.jl` already documents and
  relies on; this compiles `predictive_plan`/`_run_campaign_predictive`/
  `_PredictiveGuardState` without ever touching `Distributed`.
- `_warm_mixed_dispatch_campaign()`: the same two samples through the
  default (non-predictive) bandit path, also with `process_max_workers=1`.
  Two samples never clear `OuterRouteTuning.mc_process_min_samples` (16) or
  `mc_process_min_mission_s` (3600 s) regardless of tuning
  (`_mc_process_worth_exploring`), so the `:process`/mixed-dispatch branch
  inside `_run_campaign_with_route_env` never runs at warmup time -- but
  Julia compiles a method's whole body on first call, branches included, so
  that branch's code (and `_run_monte_carlo_mixed`'s method for this call
  shape) still lands in the pkgimage even though it is never executed here.
  This call's job is to compile `_campaign_route_plan` and
  `_run_campaign_with_route_env` themselves, which `_warm_campaign_dispatchers`
  never reaches.

Both calls use a throwaway `OuterRouteState()` rather than the process
default, so a synthetic 2-sample timing from precompilation never pollutes
the route bandit's history for a real campaign later in the same process.
Both are wrapped in `withenv` forcing
`SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS=0` /
`SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST=0`, and the `@setup_workload` block
now also calls `SimulationCampaigns.reset_predictive_machine_constants!()`
and `SimulationCampaigns.reset_campaign_route_state_persistence!()` after the
workload runs, for the same reason the existing
`SimulationModel.Planets._reset_furnished_kernels!()` call is there: both
`predictive_machine_constants()` and `ensure_campaign_route_state_loaded!`
cache into module-level `const`s the first time a campaign reads them, and
those caches would otherwise be serialised into the pkgimage -- this
precompiling machine's calibration file (or the absence of one), or an
"already loaded" flag with nothing actually loaded, served to every process
that later loads the same pkgimage. The `withenv` override already prevents
either cache from being *populated* with anything meaningful during a normal
`Pkg.precompile()`, but the reset is cheap insurance regardless, matching the
SPICE-furnish precedent in the same block.

### Measurement: first-campaign latency, before vs. after, ratio only

`@elapsed` of the first `run_monte_carlo(seed -> seed * 2, 1:2; threads=:auto)`
call after `using SpaceAGORA`, two fresh processes each, 1 thread,
space-falcon-1 (shared with other WS11 jobs throughout):

| | Run 1 | Run 2 |
|---|---|---|
| Before | 3.944 s | 3.703 s |
| After | 0.653 s | 0.664 s |

Ratio (mean before / mean after): **~5.81x**. Raw numbers:
`benchmarks/studies/small_wins/results/precompile_first_campaign_latency.csv`.
Reproduction script: `benchmarks/studies/small_wins/first_campaign_latency.jl`.

### Confirmed

- `julia --project=. -e 'using SpaceAGORA'` precompiles cleanly (56 s wall on
  this run; the pre-change tree precompiled in 50 s, so the added warmup
  costs about 6 s of precompilation time in exchange for the ~3.2 s/campaign
  saved above -- the trade this item exists to make). Validated via the
  shared package depot's dependency cache (503 already-precompiled
  dependencies, unaffected by this change) rather than a fully isolated
  fresh depot, given the machine was shared with several other WS11 agents'
  Julia jobs throughout this session; the dependency set and versions are
  pinned in `Manifest.toml` and identical either way, so this is the same
  validation CI's precompile step performs for the SpaceAGORA package's own
  compilation, modulo dependency-cache warmth (which this change does not
  touch).
- `test/unit/parallel/predictive_campaign_tests.jl` (69 checks across its
  testsets), `test/unit/parallel/outer_route_persistence_tests.jl`,
  `test/unit/parallel/predictive_planner_tests.jl`, and
  `test/unit/parallel/mixed_dispatch_tests.jl` all pass against the final
  tree.

## Tests run (final tree)

| File | Result |
|---|---|
| `test/unit/simulation/harmonics_batch_parity_tests.jl` | 22 passed |
| `test/unit/simulation/third_body_route_parity_tests.jl` | 84 passed |
| `test/unit/dynamics/third_body_cost_tests.jl` | 732 passed |
| `test/unit/dynamics/harmonics_kernel_identity_tests.jl` (new) | 3493 passed |
| `test/unit/parallel/predictive_campaign_tests.jl` | all passed |
| `test/unit/parallel/outer_route_persistence_tests.jl` | all passed |
| `test/unit/parallel/predictive_planner_tests.jl` | all passed |
| `test/unit/parallel/mixed_dispatch_tests.jl` | all passed |
| `benchmarks/studies/small_wins/third_body_sample_allocation.jl` | 1 passed (192 B/call, threshold 512 B/call) |

## What was skipped, and why

The WS11 common contract also asks to run "the harness `--parity` path"
(`benchmarks/studies/parallelization_performance.jl`). A `test`-profile
invocation was started (`--threads=2 --repeats=1 --warmup=0`) but did not
finish inside a 300 s window on this shared machine (the profile's default
mode/case matrix -- five modes, `montecarlo_multi_sat` -- is heavier than
that); it was not re-run at a longer timeout to avoid adding more load to a
machine already running several other WS11 agents' benchmark jobs
concurrently. The contract's own designated tool for this workstream's proof
obligation, `benchmarks/studies/third_body_cost/variants.jl --dump`, was used
instead and is what the dump/cmp evidence above is built on; the existing
unit-level parity suites (`harmonics_batch_parity_tests.jl`,
`third_body_route_parity_tests.jl`) also passed unchanged. The
`parallelization_performance.jl --parity` path was not exercised standalone
this session.

The 256-spacecraft exponential-atmosphere aero constellation, the 64-spacecraft
look-ahead-cache native GRAM constellation, one `montecarlo_heavy_aerobraking`
sample, and the `mcgrid_8sat_16mc` campaign (predictive/outer_threads modes)
from the WS11 common contract's reference case set were not dumped separately:
none of this workstream's three changes touch density sampling, GRAM, or
per-satellite RHS code outside `sample_third_body_ephemerides` and the
harmonics batch kernel, so those cases exercise code this session did not
modify. The cases dumped (256-spacecraft L50 vacuum, the 4096-spacecraft P2
case, and the 256-spacecraft SRP + third-body constellation) are exactly the
ones the WS11 common contract lists as relevant to a gravity-harmonics and
third-body-sampling change.
