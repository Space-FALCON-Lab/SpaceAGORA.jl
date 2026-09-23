# Aerodynamics and density at constellation scale

Design record for batching the aerodynamic effector and the RHS-side
atmosphere pre-sample over a constellation, under the hard requirement
that the dynamics stay **bit-identical** to `adb343566` on
`policy-v2-consistency`. A change that moves any saved state of any trajectory
by one bit is not shipped, however large the speedup.

Status: implemented and measured. Every number below comes from a run of the
scripts in `benchmarks/studies/aero_batch/` on `space-falcon-1`, with the raw
rows in `benchmarks/studies/aero_batch/results/`; timings are ratios between
back-to-back runs on a machine shared with other Julia jobs, never absolute
benchmarks. Section 7 has the results.

## 1. The measured problem, and what is and is not sourced

SOURCED, from the P6 calibration on TRX50, job `20260922-054218-1797061`
(single thread, 4096 spacecraft, 100 s mission): the exponential-atmosphere
aero case `aero_4096sat_l50_expatm_100s` cost 16.27 s against 11.71 s for the
vacuum case `gravity_4096sat_l50_vacuum_5800s` at the same constellation size.

DERIVED from those two numbers: the aero case is 1.39x the vacuum case, and
the 4.56 s difference is 28 percent of the aero case's serial wall time. That
difference is everything the atmosphere adds — the `AerodynamicCoefficientfM`
effector, the RHS-side atmosphere pre-sample, the density `DiscreteCallback`,
and the extra solver steps the smaller `dt_max_orbit` of that case implies.
Splitting the 4.56 s among those four is exactly what deliverable 1 (the
attribution run) exists to do, and nothing in this document assumes a split
before that run has happened.

ASSUMED, and named so it can be checked: that the split at 1024 spacecraft
resembles the split at 4096. The attribution runs at 1024 because that is the
size the reference case set uses and the size this workstation can hold; the
4096 number above is quoted only as the origin of the task.

## 2. Current per-spacecraft evaluation and its order of operations

The constellation cases that carry aerodynamics run the flat constellation
effector queue (`_spacecraft_dynamics_flat_constellation_effector_queue!`,
`src/simulation/engine/dynamics_rhs.jl`). For an effector tuple of
`(GravitationalHarmonicsModel, AerodynamicCoefficientfM)`, one derivative
evaluation does this, in this order:

1. `_prefill_shared_body_samples!` — serial. Warms the SPICE memo for solar
   and third-body positions and computes the harmonics LPI once. With no SRP
   and no third body it only takes the harmonics branch.

2. `_prefill_environment_samples!(p, t, sc_state; atmosphere=true)` — the
   atmosphere pre-sample. `l_pi` is computed once, then
   `threaded_foreach_worker_persistent(:rhs_atmosphere, num_sats, ...)` runs,
   per satellite:
   - `sample_planet_frame_with_lpi` → `_planet_relative_state`, `rtolatlong`;
     writes `pos_pp`, `vel_pp`, `alt_m`, `lat_rad`, `lon_rad` into the
     per-satellite component buffers;
   - `_sample_atmosphere_from_planet_frame(...; write_buffers=true)`
     (`src/simulation/engine/effector_sampling.jl`), which per satellite
     re-reads `_callback_env_config(p)`, re-resolves
     `_density_model_for_sat(p, sat_idx)`, re-derives
     `target_include_j2` (a short-circuited scan of the effector tuple),
     re-reads `cache_cfg`/`stats_enabled`/`caches`, extracts `pos_ii`,
     `vel_ii` and the mass from the state view again, calls
     `_density_state_from_kinematics!`, and writes
     `densities`/`temperatures`/`winds`/`density_sample_t` through
     `_write_density_buffers!`.

   The GRAM process density service (`_rhs_density_service_candidate` /
   `_rhs_density_service_fill!`) already has a batch point here, decided once
   per evaluation; it is not taken for an analytic model.

3. `_accumulate_dynamic_effectors_flat_slots!`:
   - `_prefill_rhs_flat_state_samples!` fills the contiguous
     `rhs_flat_state_pos_ii` / `_vel_ii` / `_mass_kg` (and, under
     `orientation_sim`, `_q_ib` / `_omega_body`) buffers;
   - the harmonics pre-pass (`_accumulate_harmonics_flat_batch!`) writes
     `slots[1:3, harmonics_eff_idx, sat_idx]` for every active satellite;
   - aero is **not** pre-passed. It stays in the per-`(satellite, effector)`
     flat queue: `_prepare_rhs_flat_work_items!` emits one work item per
     active satellite, and `threaded_foreach_worker_persistent(:rhs_flat_queue,
     count_items, ...)` dispatches them one at a time. Each item then does:
     `dynamic_effectors[eff_idx]` (a runtime index into a *heterogeneous*
     tuple), `_wrench_method_available(effector)` (a `hasmethod` call on a
     value whose type is not known at that point),
     `_rhs_flat_state_sample_from_buffers`, and the `@noinline`
     `_evaluate_dynamic_effector`, which is a dynamic dispatch.
   - `_evaluate_dynamic_effector` calls
     `sample_environment_with_reusable_buffers`, which for
     `environment_requirements(::AerodynamicCoefficientfM) =
     (planet_frame=true, atmosphere=true)` reads the prefilled planet frame
     (`sample_buffered_planet_frame`) and the prefilled atmosphere
     (`sample_buffered_atmosphere`) and builds an `EnvironmentSample`, then
     calls `wrench_caching!(model, state_sample, env, t, p, sat_idx)`.
   - `wrench_caching!` calls `_aero_pure_wrench(:fm, x, env, nothing,
     model.fixed_attitude_incidence)` and then `_store_aero_caches!`.

4. `_aero_pure_wrench` (`src/dynamics/coupled/aerodynamic_wrench_models.jl`)
   is the arithmetic. Per satellite: the vacuum/zero-airspeed early returns,
   `cross`/`norm` for `h_pp`, `latlongtoNED`, the wind-relative velocity,
   `lift_pp_hat`/`drag_pp_hat`/`cross_pp_hat`, `l_pi'`, `θ_body`, then a loop
   over the spacecraft's links accumulating `force_ii`, `torque_body`,
   `drag_ii`, `lift_ii`, `cross_ii` into `MVector`s with `.+=`, each link
   costing an `aerodynamic_coefficient_fM` evaluation (the Hart et al. closed
   forms: several `erf`, `exp`, `sqrt` and trigonometric calls).

5. `_reduce_flat_effector_slots!` sums the per-effector slots into the
   per-satellite totals in effector order.

Two things are per-satellite in step 3 that are properties of the *run*, not
of the satellite: the effector's identity (hence its concrete type, its
`environment_requirements`, and whether it has a `wrench` method) and the
density model. Every one of the 1024 items re-derives both.

## 3. Proposed batch layout

Two independent changes. Each is measured and proved on its own, and either
can be dropped without the other.

### A. An aero pre-pass slot in the flat route

`src/simulation/engine/dynamics_rhs.jl` gains `_aero_prepass_effector`, a
trait in the shape of the existing `_harmonics_prepass_effector`, true for
`AerodynamicCoefficientfM` and false for everything else, plus a pre-pass
`_accumulate_aero_flat_batch!(slots, eff_idx, model, sc_state, p, t, plan)`
run in `_accumulate_dynamic_effectors_flat_slots!` after the harmonics
pre-pass. Aero is then removed from the flat queue's selection mask exactly
the way the harmonics and batchable effectors already are, so no work item is
emitted for it and nothing is evaluated twice. The harmonics pre-pass is not
touched.

The pre-pass takes contiguous slices of the active-satellite list, one per
worker, and each worker runs, for each satellite in its slice, exactly:

```
state_sample = _rhs_flat_state_sample_from_buffers(shared_buffers, spacecraft, sat_idx, orientation_sim)
planet_frame = _sample_reusable_planet_frame(req, sc_view, p, sat_idx, t)
atmosphere   = _sample_reusable_atmosphere(req, sc_view, planet_frame, p, sat_idx, t)
env          = EnvironmentSample(planet, planet_frame, atmosphere, nothing, nothing)
force, torque, drag, lift, cross = _aero_pure_wrench(:fm, state_sample, env, nothing, incidence)
_store_vector_cache!(drag_cache, sat_idx, drag)   # then lift, then cross
slots[1:6, eff_idx, sat_idx] = (force..., torque...)
```

This is what the queue item's `sample_environment_with_reusable_buffers` →
`wrench_caching!` chain does for this effector, written out: the environment
sample for its requirements (planet frame and atmosphere, nothing else), the
same `_aero_pure_wrench` method with the per-link hook off (the trait
guarantees it), and the same three cache stores in the same order. Two
differences, neither arithmetic. `model` is a concrete `AerodynamicCoefficientfM`
in the signature instead of an element pulled at runtime out of a
heterogeneous tuple, so nothing dispatches dynamically per spacecraft and the
worker pool is entered once per slice, not once per spacecraft. And the body
lives in a separate, concretely typed `@noinline` function
(`_aero_prepass_satellite!`) with the environment sample built with concrete
field types; as a closure body calling the generic builder, it boxed the
Union-typed intermediates and, once, a copy of `ODEParams` per spacecraft.

`_aero_pure_wrench` is called, not inlined. Inlining it into the pre-pass was
tried and **rejected**: the state histories stayed byte-identical, but the
parity suite found the drag cache (a recorded output, not state) a last bit
different on a few spacecraft. StaticArrays' matrix-vector product is written
with `muladd`, and whether LLVM fuses a `muladd` into an FMA depends on the
code around it, so the same source can round differently once inlined into a
new caller. Called, it is the same compiled method instance `wrench_caching!`
reaches.

Per-satellite independence is what makes the slicing safe: each satellite
writes only `slots[:, eff_idx, sat_idx]` and its own entries of the drag/lift/
cross save caches, so the number of workers cannot change any value, and the
reduction that follows is unchanged.

Not batched, and left on the queue:
- `per_link_atmosphere` models (`_per_link_enabled(model)`), which call a
  density query per link and carry a `maxlog=1` warning whose firing order
  would otherwise become worker-dependent;
- `AerodynamicCoefficientConstant`, `AerodynamicCoefficientNoBallisticFlight`
  and `AerodynamicCoefficientMeshSurrogate`, which are left for
  later;
- the IMEX-partitioned call (`partition !== nothing`), like every other
  pre-pass.

Nothing in `_aero_pure_wrench` changes in change A. The arithmetic is the same
function, called the same number of times, with the same inputs.

#### Which route a run is actually on, and the serial admission

Change A only takes effect where the flat constellation queue is taken, and a
serial aero constellation is not on it. `_rhs_execution_plan`
(`src/simulation/engine/setup.jl`) admits the flat route at a thread budget of
one only when `_rhs_all_prepass_effectors(dynamic_effectors)` holds, and that
predicate lists the batchable kernels and the harmonics pre-pass, so a
`(GravitationalHarmonicsModel, AerodynamicCoefficientfM)` stack fails it today
and routes to `:satellite_batch`. A forced `SPACEAGORA_RHS_EXECUTION_MODE=flat`
does not override it either: that branch has the same `budget <= 1` fallback.

The batch was therefore first shipped on its own, with 1-thread ratios of
1.00 as predicted (section 7), and the route admission followed as a separate
change with its own proof: `_rhs_all_prepass_effectors` now also accepts
`_aero_prepass_effector`, so at budget 1 the auto route sends a
`(harmonics, aero)` constellation to the flat queue, where the queue is never
built, the harmonics coefficient sweep is batched and the atmosphere is sampled
once per evaluation. Kept separate because it changes the *route* a serial
aero constellation takes, which is a different claim from "the batch
reproduces the queue". Three things about it:

- The budget-1 branch is gated by `flat_min_sats` (24 by default) only. The
  64-spacecraft-per-worker floor of `_rhs_flat_default_admits` is applied in
  the single-harmonics branch when `viable_workers >= 2` and in the generic
  multi-thread flat branch, both of which a budget-1 plan never reaches; at an
  allotment of one the flat route spawns no workers, so there is nothing for
  the floor to amortize.
- An explicit `SPACEAGORA_RHS_EXECUTION_MODE=flat` request at budget 1 still
  falls back to the per-satellite batch; only the auto route admits it.
- A per-link-atmosphere fM model is not a pre-pass effector, so a stack
  carrying one stays on the per-satellite batch at budget 1.

The same is true of the RHS-side atmosphere pre-sample that change B batches:
`_prefill_environment_samples!` is called only from
`_spacecraft_dynamics_flat_constellation_effector_queue!`. On the per-satellite
route the atmosphere is sampled inside each effector's own
`sample_environment_with_reusable_buffers` call and there is no batch point to
take. Change B therefore has the same route scope as change A.

### B. One density query per evaluation for a shared analytic atmosphere

`src/simulation/engine/effector_sampling.jl` gains
`_uniform_light_density_model` and `_fill_uniform_light_atmosphere!`, used by
`_prefill_environment_samples!`. When the constellation shares one light,
cache-free density model, the pre-sample's loop fills the planet frames only,
and one `getDensityBatch!` over the just-filled `alt`/`lat`/`lon` buffers
writes `shared_buffers.densities` / `temperatures` / `winds`, followed by
`_write_density_time_buffers!`. Otherwise the per-spacecraft sample runs
exactly as before. (A second, "hoisted per-spacecraft" route was designed but
not built: the uniform route covers the analytic cases, and every other case
keeps today's code untouched.)

The uniform route is taken only when all of the following hold:

- every spacecraft is active (the per-spacecraft path skips inactive ones and
  leaves their buffers alone; a batch over 1:N would not);
- GRAM runtime statistics are off (the per-spacecraft path counts calls);
- the model is not a GRAM grid snapshot;

- every active satellite resolves to the same density model object
  (`_density_batch_model_for_callback` returns non-`nothing`);
- `!density_model_work_is_heavy(model)` — i.e. not native GRAM. GRAM's
  locking is not touched here, and serializing native
  GRAM into a one-thread `getDensityBatch!` would be a regression rather than
  a win;
- `!cb_env.density_freeze_per_step` and `!cb_env.vacuum_gram_cache_enabled` —
  both make the per-satellite value depend on per-satellite cache state
  rather than on the current position;
- `!_gram_track_cache_enabled(cb_env.gram_track_cache, model)` — same reason,
  and the same guard the density callback's own batch route already applies;
- the density service route (`batch_atmosphere`) was not selected;
- `densities`, `temperatures`, `winds` and `density_sample_t` are all at
  least `num_sats` long, which is the same length guard
  `_write_density_buffers!` applies per satellite.

For `ExponentialAtmosphereModel`, `getDensityBatch!`'s specialized method
evaluates `_exponential_density(model.ρ_ref, model.h_ref, model.H, h)`,
`model.temperature_k` and a zero wind — the same three expressions, in the
same order, that the scalar `getDensity` method evaluates. The `Float64(hs[i])`
in the batch method is the identity on the `Vector{Float64}` altitude buffer.
The two are therefore the same arithmetic, which is what the parity test
asserts on a sweep of altitudes rather than on the author's reading of the
source.

This change is in the *density sampling path*, not in the density models:
`src/environment/atmosphere/density_models.jl` is not edited.

One more allocation in the same file, found by the attribution: every
atmosphere sample evaluated `_uses_j2_gravity_effector(dynamic_effectors)`,
which iterates the heterogeneous effector tuple and boxed an element per
spacecraft per sample. `_effectors_use_j2_gravity` peels the tuple instead and
returns the same Bool. This one is on both routes, which is why the serial
allocation halves (section 7).

## 4. Why every operation stays in the same order, with no reassociation and no FMA

The bit-identity argument is deliberately structural rather than numerical.
Neither change rewrites an expression.

- **No expression is rewritten.** Change A calls `wrench_caching!`, which
  calls `_aero_pure_wrench`, unmodified; change B either calls the same
  `_density_state_from_kinematics!` with hoisted arguments or calls
  `getDensityBatch!`, whose specialized methods are the scalar methods'
  expressions verbatim. There is no new sum, no new product and no new
  temporary, so there is nothing for the compiler to reassociate.
- **No `@fastmath`, no `@turbo`, no `@simd`.** The pre-pass loop is a plain
  `for` over a satellite range. `@simd ivdep` was appropriate on the harmonics
  batch kernel because its inner loops are independent per batch slot with no
  loop-carried reduction; here the per-satellite body contains branches, calls
  to `erf`/`exp`/`atan` and writes through `p`, so an annotation would buy
  nothing and would be a license this document would then have to justify.
  None is used.
- **No FMA is introduced, and no arithmetic is inlined into a new caller.**
  Julia does not contract `a*b+c` into `fma` unless `@fastmath` or an explicit
  `muladd`/`fma` is written, and neither change writes one. But library code
  the wrench already calls does: StaticArrays' matrix-vector product uses
  `muladd`, which LLVM *may* fuse, and whether it does depends on the
  surrounding code. So "the same source" is not enough; the arithmetic must
  stay in the same compiled method. That is why `_aero_pure_wrench` is called
  from the pre-pass rather than inlined into it (inlining it measurably moved
  a last bit of the drag cache; section 3A).
- **Accumulation order within a satellite is untouched.** `_aero_pure_wrench`
  still accumulates `force_ii`, `drag_ii`, `lift_ii`, `cross_ii` and
  `torque_body` over `spacecraft.links` in link order with `.+=` into
  `MVector`s. The link-level collect-then-sum in the `calcForceTorque` path is
  not on the route this change takes and is not modified.
- **Accumulation order across satellites does not exist.** Each satellite's
  six numbers go into its own `slots[:, eff_idx, sat_idx]` and are summed
  later by the unchanged `_reduce_flat_effector_slots_range!`, in effector
  order, from zero. Changing how satellites are grouped into workers therefore
  cannot change a result — which is the same property the harmonics pre-pass
  and the SRP/n-body/inverse-square pre-passes already rest on.
- **Worker count independence.** The pre-pass reads only per-satellite inputs
  and writes only per-satellite outputs. The one shared mutable structure it
  touches is the drag/lift/cross save cache, and `_store_vector_cache!`
  already gives every satellite its own index.
- **The queue's evaluation of aero is removed, not duplicated.** Aero is
  added to the same selection-mask exclusion and the same
  `_count_flat_queue_only_effectors` accounting the other pre-passes use, so
  the effector is evaluated exactly once per satellite per evaluation, as
  before.

Where bit-identity cannot be argued structurally the case is excluded rather
than approximated: per-link atmosphere, native GRAM, the frozen and look-ahead
density modes, and the IMEX partition all keep today's path.

## 5. What is measured, and how

Measurement rules: explicit
`--threads`, at most 8, one Julia process at a time, `uptime` and a `julia`
process check before any timing, and ratios between two runs in one process
state — never absolute benchmarks.

**Deliverable 1, attribution** (`benchmarks/studies/aero_batch/attribution.jl`,
1 thread, 1024 spacecraft, exponential atmosphere): a `Profile` sample of the
`aero_1024sat_l50_expatm_*` case, reduced to the share of samples under
(a) `_aero_pure_wrench` and below, (b) the RHS-side atmosphere pre-sample
(`_prefill_environment_samples!` and below), (c) the density `DiscreteCallback`
body, (d) the flat-queue plumbing (`_evaluate_dynamic_effector`,
`_prepare_rhs_flat_work_items!`, the worker dispatch), and (e) the harmonics
pre-pass, plus the same case's vacuum counterpart as the baseline. Allocation
per spacecraft per derivative evaluation comes from the run's total
`@timed` bytes divided by `sol.stats.nf * num_sats`, reported alongside a
`Profile.Allocs` attribution of the largest sites.

**Deliverables 3 and 4, identity and ratios**
(`benchmarks/studies/aero_batch/variants.jl`): the same script dumps the full
state history (`--dump`) and times solves (`--repeats`). Identity cases, all
at 1 thread so the two dumps are comparable by construction:

| Case | N | Why |
|---|---|---|
| `aero_64sat_l50_expatm_100s` | 64 | the 64-spacecraft exponential aero reference constellation |
| `atmo256_exponential_10min` | 256 | the catalog's registered 256-spacecraft exponential rung |
| `aero_1024sat_l50_expatm_100s` | 1024 | the 1024-spacecraft exponential case |
| `montecarlo_heavy_aerobraking` | 1 | Mars, exponential atmosphere, inverse-square gravity — the single-spacecraft path, which must be untouched |
| `multi_64_high_fidelity` | 64 | harmonics + SRP + aero, i.e. aero next to a *batchable* effector in the same tuple |

"Before" is this branch at its base commit, `adb343566`, run before any source
change; "after" is the tip. Both dumps are byte-compared with `cmp`. Any
difference at all stops the change.

Ratios: 256 and 1024 spacecraft, 1 thread and 8 threads, before and after,
back to back in one process state, reported as before/after ratios with the
raw seconds in `benchmarks/studies/aero_batch/results/*.csv` (force-added,
since the repository gitignores `*.csv`). The 1-thread rows are reported even
though the route caveat above predicts 1.00 for them: a 1-thread ratio that is
*not* 1.00 would mean the change reached a route it was not supposed to reach,
which is worth knowing.

The attribution profile is taken twice for the same reason — once serial
(`--mode=serial --threads=1`, the per-satellite route the P6 serial number
was measured on) and once at 8 threads (`--mode=inner_only --threads=8`, the
flat route the pre-passes live on). They are different code and neither
substitutes for the other.

The GRAM look-ahead constellation and the `mcgrid_8sat_16mc` campaign from the
reference set are covered by the identity argument rather than by
a dump only if neither change can reach them; if either can, they are dumped
too, and the final report says which were skipped and why.

## 6. What the parity test checks (implemented)

`test/unit/dynamics/aero_batch_parity_tests.jl`, standalone (`using
SpaceAGORA`), in the shape of
`test/unit/simulation/harmonics_batch_parity_tests.jl` — `===` on `Float64`,
not `≈`:

1. **The aero pre-pass reproduces the queue's per-satellite call.** For a
   constellation of several spacecraft with varied altitudes and velocities,
   the slot values written by `_accumulate_aero_flat_batch!` are compared
   component by component with `wrench_caching!` called directly per
   satellite, for `fixed_attitude_incidence` in `:max_drag`, `:attitude` and
   `:tumbling_average`, and for single-link and multi-link spacecraft.
2. **Worker count does not change a value.** The same comparison across
   several worker allotments and several satellite counts, including counts
   that do not divide evenly into slices, and including a satellite that is
   inactive.
3. **The vacuum and zero-airspeed early returns still return zero** through
   the batch, so a satellite above the atmosphere contributes exactly zero and
   not a denormal.
4. **`getDensityBatch!` equals `getDensity`, bit for bit,** for
   `ExponentialAtmosphereModel` and `PiecewiseExponentialAtmosphereModel` over
   an altitude sweep that crosses the piecewise model's layer boundaries —
   the assumption change B's uniform-light route rests on, asserted rather
   than read off the source.
5. **The two RHS routes agree end to end.** The same short constellation
   solve under `SPACEAGORA_RHS_EXECUTION_MODE=satellite` and `=flat` produces
   the same terminal state bit for bit, which is the whole-trajectory version
   of checks 1-3.

The dump-and-`cmp` evidence of section 5 remains the primary proof; the unit
test is the regression guard that keeps it true.

## 7. Results

### Attribution (before the change)

`attribution.jl --n=1024 --mission=100`, base commit, two routes. Shares are of
one profile's samples, attributed to the innermost phase marker.

Serial (`--mode=serial --threads=1`, i.e. the per-satellite route):

| Phase | Share |
|---|---|
| harmonics, per-satellite scalar kernel | 66.6% |
| density `DiscreteCallback` | 18.2% |
| stage heat rates | 11.0% |
| aero effector (dispatch + wrench) | 4.1% |
| other | 0.2% |

The serial exponential-aero case costs 7.3x the vacuum rung per derivative
evaluation, and almost none of that is aerodynamics: at one thread the
`(harmonics, aero)` stack is routed to `:satellite_batch`, so harmonics runs
the per-satellite scalar kernel instead of the batched pre-pass the vacuum
rung gets. That is the route issue of section 3A. Allocation: 2471 B per
spacecraft per evaluation against 255 B for vacuum, 72% of it at
`getindex(::Tuple)` — runtime indexing of the heterogeneous effector tuple.

At 8 threads (`--mode=inner_only`, the flat route) the time profile is
dominated by worker-idle samples and is not informative; the allocation
profile is: 81% at `getindex(::Tuple)`, which is the flat queue's per-item
`dynamic_effectors[eff_idx]`, exactly what the aero pre-pass removes.
5094 B per spacecraft per evaluation.

Raw: `results/attribution_1024_*_before.csv`, flat profiles
`profile_aero_1024sat_l50_expatm_100s_*_before.txt`.

### Bit-identity

"Before" is `adb343566`; "after" is the tip. All pairs compared with `cmp`.

State histories (`variants.jl --dump`), every saved time, every component of
every spacecraft:

| Case | N | Threads, route | Bytes | `cmp` |
|---|---|---|---|---|
| `aero_64sat_l50_expatm_100s` | 64 | 1, satellite | 106704 | identical |
| `atmo256_exponential_10min` | 256 | 1, satellite | 2065392 | identical |
| `aero_1024sat_l50_expatm_100s` | 1024 | 1, satellite | 1638600 | identical |
| `multi_64_high_fidelity` | 64 | 1, satellite | 751032 | identical |
| `montecarlo_heavy_aerobraking` (Mars) | 1 | 1, satellite | 1619136 | identical |
| `aero_64sat_l50_expatm_100s` | 64 | 8, flat (forced) | 106704 | identical |
| `atmo256_exponential_10min` | 256 | 8, flat (forced) | 2065392 | identical |
| `aero_1024sat_l50_expatm_100s` | 1024 | 8, flat (forced) | 1638600 | identical |
| `multi_64_high_fidelity` | 64 | 8, flat (forced) | 751032 | identical |

The 8-thread base dump was taken twice and the two copies were identical, so
the flat route is itself deterministic at 8 threads. At the base commit the
1-thread satellite-route dump and the 8-thread flat-route dump of each case
were also byte-identical: the two routes produce the same trajectory, which is
the evidence for admitting aero to the budget-1 flat route in `setup.jl`.

Derivative, drag/lift/cross caches and density/temperature/wind/sample-time
buffers (`rhs_dump.jl`, six evaluations per case along a fixed-step path), for
the four constellation cases above: identical before and after on the flat
route at 8 threads and on the satellite route at 1 thread.

Not dumped: the look-ahead-cache native GRAM constellation, the harmonics-only
vacuum and SRP + third-body constellations, and the `mcgrid_8sat_16mc`
campaign. Neither change reaches them except through the J2 check, whose Bool
is unchanged: the aero pre-pass requires `AerodynamicCoefficientfM` in the
stack, and the uniform atmosphere route excludes native GRAM, GRAM grids, the
look-ahead and track caches, and freeze-per-step. `mcgrid` runs
`(harmonics, aero)` on an exponential atmosphere, i.e. the same code path as
the dumped cases.

### Ratios

`variants.jl`, `atmo256_exponential_10min` and
`aero_1024sat_l50_expatm_600s`, five repeats per process, best of five, in the
order after, before, after, before. Ratio is before over after.

| Case | N | Threads | Round 1 | Round 2 | Allocation, before / after |
|---|---|---|---|---|---|
| `atmo256_exponential_10min` | 256 | 1 | 1.05 | 0.97 | 2.01 |
| `aero_1024sat_l50_expatm_600s` | 1024 | 1 | 1.01 | 0.92 | 2.02 |
| `atmo256_exponential_10min` | 256 | 8 | 1.18 | 1.22 | 1.59 |
| `aero_1024sat_l50_expatm_600s` | 1024 | 8 | 1.36 | 1.33 | 2.06 |

At one thread wall time is unchanged within noise, as the route analysis
predicts; allocation halves anyway because of the J2 check. At eight threads,
on the flat route, the 1024-spacecraft case is 1.33-1.36x faster and allocates
half as much. Raw: `results/ratio_{1,8}t_{before,after}_r{1,2}.csv`.

### Follow-up changes

Three changes followed the batch, each proved against the commit before it.
Note that the catalog's P6 density cases (`aero_<N>sat_l50_expatm_<S>s`) were
moved onto a shared constellation geometry between the batch and these
changes, so their dumps below are not comparable with the table above; each
before/after pair uses the same catalog.

**Serial route admission** (`setup.jl`, section 3A). State histories at 1
thread:

| Case | N | Bytes | `cmp` |
|---|---|---|---|
| `aero_64sat_l50_expatm_100s` | 64 | 106704 | identical |
| `aero_256sat_l50_expatm_100s` | 256 | 458976 | identical |
| `aero_1024sat_l50_expatm_100s` | 1024 | 1704144 | identical |
| `atmo256_exponential_10min` | 256 | 2065392 | identical |
| `montecarlo_heavy_aerobraking` (Mars) | 1 | 1619136 | identical |

The recorded results tables (every save field, including drag/lift/cross,
heat rate, winds and the visualization density column), dumped by
`outputs_dump.jl` for the 64 and 256 aero cases and `atmo256`, are
byte-identical too. `rhs_dump.jl` on the 256 and 1024 cases shows the
derivative and the drag/lift/cross caches identical and the shared density
buffers different: the per-satellite route never writes them inside an RHS
call (they keep their initial values between density-callback firings), the
flat route writes the current stage's sample. No output reads them between
those points — the density callback rewrites them before every save, which is
why the recorded density column does not move — and the integrated state does
not depend on them.

Serial ratio, `variants.jl --mode=serial --threads=1`, best of five, rounds
ordered after, before, after, before; ratio is before over after:

| Case | N | Round 1 | Round 2 | Allocation, before / after |
|---|---|---|---|---|
| `aero_256sat_l50_expatm_600s` | 256 | 2.26 | 2.43 | 0.62 |
| `atmo256_exponential_10min` | 256 | 2.30 | 2.43 | 0.62 |
| `aero_1024sat_l50_expatm_600s` | 1024 | 2.68 | 2.79 | 0.63 |

Step counts are unchanged (same `nf`). The serial run allocates about 1.6x
more after the admission: the flat route's slot, state and planet-frame
buffers and its per-call planning replace the per-satellite loop's stack
temporaries. It is still 2.3-2.8x faster.
Raw: `results/serial_route_1t_{before,after}_r{1,2}.csv`.

**Cost-model mirror** (`src/parallel/cost/work_counts.jl`):
`flat_queue_node_effector` now knows the aero pre-pass (and its per-link
exclusion), so `constellation_work_counts` stops charging N aero queue nodes
the queue no longer builds. Routing input only; no dynamics path reads it.

**Slot-reduction mask** (`dynamics_rhs.jl`): `_reduce_flat_effector_slots_range!`
evaluates `_flat_slot_selected` once per effector into an `NTuple{N, Bool}`
instead of indexing the heterogeneous effector tuple once per spacecraft per
effector; loops, statements and summation order unchanged. State histories of
the 1024 and 256 aero cases and `multi_64_high_fidelity` at 8 threads (forced
flat) and of the 1024 and 256 cases at 1 thread, and `rhs_dump.jl` outputs of
the 1024 case and `multi_64_high_fidelity` at 8 threads, are byte-identical to
the commit before. With the flat route forced at a fixed budget of 8, per
derivative evaluation at 1024 spacecraft, allocation falls from 2322 to 1121 B
per spacecraft and the evaluation from 5241 to 4790-4926 us (one process each,
back to back); on `atmo256` from 3570 to 2370 B. A whole 600 s solve of the
1024 case at 8 threads in `inner_only` mode: best of three 2.342 s against
2.480 s (1.06x), allocation 1.79 GiB against 5.27 GiB.
Raw: `results/mask_8t_{before,after}.csv`.

Whole-solve allocation at 8 threads under `inner_only` depends on the plan
the run's calibration chose, and varies between processes: one pair of 100 s
solves in the mask dumps showed the opposite order (0.297 GiB with the old
reduction, 0.805 GiB with the mask). The fixed-budget per-evaluation numbers
above are the controlled comparison; the whole-solve 8-thread allocation
ratios in the table earlier in this section should be read with this in mind.

