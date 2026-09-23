# Aerodynamics and density at constellation scale

Design record for batching the aerodynamic effector and the RHS-side
atmosphere pre-sample over a constellation, under the hard requirement
that the dynamics stay **bit-identical** to `adb343566` on
`policy-v2-consistency`. A change that moves any saved state of any trajectory
by one bit is not shipped, however large the speedup.

Status when this section was written: design only. Every row marked TO BE
MEASURED is filled from a run of the scripts in `benchmarks/studies/aero_batch/`,
never from an estimate.

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
env          = sample_environment_with_reusable_buffers(req, model, sc_view, p, sat_idx, t)
force, torque = SimulationModel.wrench_caching!(model, state_sample, env, t, p, sat_idx)
slots[1:6, eff_idx, sat_idx] = (force..., torque...)
```

— the same four calls the queue item makes, in the same order, with the same
arguments. The difference is that `model` is a concrete
`AerodynamicCoefficientfM` in the pre-pass's signature rather than an element
pulled at runtime out of a heterogeneous tuple, so `_wrench_method_available`,
`environment_requirements`, `solver_partition` and the whole
`sample_environment_with_reusable_buffers` → `wrench_caching!` →
`_aero_pure_wrench` chain specialize and inline instead of dispatching
dynamically once per satellite; and that the worker pool is entered once per
slice instead of once per satellite.

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

#### Which route a run is actually on, and the one line outside these files

Change A only takes effect where the flat constellation queue is taken, and a
serial aero constellation is not on it. `_rhs_execution_plan`
(`src/simulation/engine/setup.jl`) admits the flat route at a thread budget of
one only when `_rhs_all_prepass_effectors(dynamic_effectors)` holds, and that
predicate lists the batchable kernels and the harmonics pre-pass, so a
`(GravitationalHarmonicsModel, AerodynamicCoefficientfM)` stack fails it today
and routes to `:satellite_batch`. A forced `SPACEAGORA_RHS_EXECUTION_MODE=flat`
does not override it either: that branch has the same `budget <= 1` fallback.

Consequently:

- At a multi-thread budget (the 8-thread measurements, and every constellation
  run that the router sends to the flat queue) change A applies.
- At one thread it does not, and the 1-thread ratio is expected to be 1.00
  unless `_rhs_all_prepass_effectors` also learns the aero trait.

That is a one-line change in `setup.jl`, and it is deliberately kept out of
this change, for two reasons. It changes the *route* a serial aero
constellation takes, which is a different claim from "the batch reproduces the
queue" and needs its own before/after dump; and conflating the two would make
a single ratio unreadable. The parity test's route-equivalence check (section
6, item 5) is the evidence that would support it.

The same is true of the RHS-side atmosphere pre-sample that change B batches:
`_prefill_environment_samples!` is called only from
`_spacecraft_dynamics_flat_constellation_effector_queue!`. On the per-satellite
route the atmosphere is sampled inside each effector's own
`sample_environment_with_reusable_buffers` call and there is no batch point to
take. Change B therefore has the same route scope as change A.

### B. One density-model resolution per evaluation instead of per satellite

`src/simulation/engine/effector_sampling.jl` gains a batched atmosphere
pre-sample used by `_prefill_environment_samples!`. Per evaluation, once:
resolve `cb_env`, the uniform density model, `cache_cfg`, `stats_enabled`,
`target_include_j2` and `caches`; then either

- the **uniform-light route**: a single `getDensityBatch!` over the
  already-filled `alt`/`lat`/`lon` component buffers into
  `shared_buffers.densities` / `temperatures` / `winds`, followed by
  `_write_density_time_buffers!`; or
- the **hoisted per-satellite route**: the existing per-satellite call, with
  the six run-scoped values passed in rather than re-derived inside the loop.

The uniform-light route is taken only when all of the following hold, and the
hoisted route (which is a pure hoist, no branch change at all) otherwise:

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
  nothing and would be a licence this document would then have to justify.
  None is used.
- **No FMA is introduced.** Julia does not contract `a*b+c` into `fma` unless
  `@fastmath` or an explicit `muladd`/`fma` is written. Neither appears in
  either change, and neither change edits an arithmetic expression that could
  be contracted.
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
`Profile.Allocs` attribution of the largest sites. TO BE MEASURED.

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

## 6. What the parity test checks

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
