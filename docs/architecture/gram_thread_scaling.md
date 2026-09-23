# Native GRAM under threads: the global lock and the isolated pool

Native GRAM is single-threaded. Every density query in a process is serialized,
by default, on one lock (`RuntimeServices.GRAM_LOCK`, which is the same object as
`SPICE_LOCK` — GRAM's statically linked CSPICE exports the same internal symbols
as SpaceAGORA's own SPICE bindings, so the two cannot be separated). A threaded
constellation therefore does not scale on a GRAM atmosphere: the threads queue.

`SPACEAGORA_GRAM_ISOLATED_POOL` is the alternative that has been in the tree,
default off, since its SPICE-corruption crash was root-caused and fixed. It
replaces the one shared `GRAMAtmosphereModel` with `workers` independent
`deepcopy`ed models, each behind its own `ReentrantLock`, and evaluates a batch
of satellites across them. Nothing had measured it, and nothing had checked that
two GRAM instances return the same number.

This page records what the measurement found. The scripts are in
`benchmarks/studies/gram_thread_scaling/`; the committed CSVs under its
`results/` are the rows these tables are read from.

## Is the pool bit-identical to the locked path?

Yes, and it has to be, because the alternative is that native GRAM's returned
density depends on the sequence of calls an instance has seen — which would make
*any* concurrency change to this path a dynamics change.

It does not. `GRAMSuite._gram_density_state_native` reads `atmos.density`, the
mean field, not `perturbedDensity`, and `_gram_wind_mode()` resolves to
`:nominal` by default, so the returned state is a function of position, time and
epoch alone. GRAM's per-call perturbation random walk advances, but nothing in
the returned tuple reads it.

Three pieces of evidence, all exact comparisons on the `reinterpret`ed bits,
never a tolerance:

| Check | What it compares | Result |
|---|---|---|
| `density_grid_parity.jl`, replay control | the locked batch call run twice on the same grid | identical |
| `density_grid_parity.jl`, per instance | each pool instance's serial evaluation of the whole grid against the locked reference | identical: 4 instances over 256 states, then 2 and 8 instances over 2048 states |
| `density_grid_parity.jl`, batch | the shipped `getDensityBatch!` against the shipped `_gram_isolated_pool_batch_eval!`, threaded | identical for density, temperature and all three wind components |

and on whole trajectories, via `dump_states.jl` (every saved time and every
component of every spacecraft, raw `Float64`, `cmp`ed byte for byte):

| Dump | Bytes | `cmp` |
|---|---|---|
| 64 spacecraft, look-ahead, locked vs pool width 4, 4 threads | 98 496 | identical |
| 64 spacecraft, freeze-per-step, locked vs pool width 4, 4 threads | 98 496 | identical |
| 64 spacecraft, look-ahead, locked vs pool width 2, 8 threads | 98 496 | identical |
| 64 spacecraft, look-ahead, locked vs pool width 8, 8 threads | 98 496 | identical |
| the WS11 reference case (`ppc_constellation`, EI 120 km), locked vs pool width 4 | 98 496 | identical |

The locked dump is also identical between the 4-thread and the 8-thread run, so
the thread count is not moving the trajectory either.

One real difference was found and fixed, and it was not a rounding difference.
The locked scalar path floors altitude at `-30.0` m before calling GRAM
(`EM.getDensity(::GRAMAtmosphereModel, ...)` and `EM._gram_point_density`);
the pool's `_gram_isolated_pool_density_state` did not. Native GRAM rejects a
height below −31 m with *"Height below -31 meters. This is an unrecoverable
error."*, so below that altitude the locked path returned a density and the
pooled path aborted the solve. A run that reaches the surface — an entry, a
landing — completed with the pool off and died with it on. The pool now applies
the same floor. `src/simulation/callbacks/density_callbacks/gram_process_batch.jl`
sends unclamped altitudes to the process-backed density service the same way and
has not been changed here.

The fix cannot move any shipped trajectory, for a reason stronger than the
dumps: the pool is off by default, and with it off `_gram_isolated_pool_batch_eval!`
returns before it reaches the clamped line at all. The dumps above were taken
before and after the change and compared byte for byte anyway, with the pool
both off and on.

## Is the pool faster?

RESULTS_PLACEHOLDER

## The interlock that makes the pool unreachable below 16 threads

This is the single most load-bearing thing to know before touching this path.

`_gram_isolated_pool_batch_eval!` takes its width from

```
_density_callback_thread_decision(p, args, num_sats; heavy_work=true).allotment
```

(`src/simulation/callbacks/density_callbacks/runtime.jl`), and that decision goes
through `ParallelPolicy.auto_thread_min_budget(:density_callback)`
(`src/parallel/policy/env_config.jl`), which is **16**. On any process with fewer
than 16 threads the decision returns `allotment = 1`, the pooled call fails its
own `workers > 1` guard, returns `false`, and the run takes the locked path —
with no warning, whatever `SPACEAGORA_GRAM_ISOLATED_POOL` is set to. Measured
directly on a 24-thread workstation at 4 threads and 64 spacecraft: the decision
returns `(use_threads = false, allotment = 1)`, and the locked and pooled arms
record the same `gram_density` lock-acquisition count to the unit.

The floor's own comment says it exists because "native/point GRAM is serialized
behind a process-wide lock, so oversubscribing it below a reasonably high thread
count wastes cycles fighting for that lock". That reasoning is sound for the
locked path and circular for the pool: the lock is why the gate exists, and the
pool is what removes the lock. The study runs both arms with
`SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET=1` so there is something to
measure; in the locked arm that setting can only widen the kinematics pre-fill,
never the GRAM evaluation, which `getDensityBatch!` performs serially at any
width.

`test/unit/environment/gram_isolated_pool_tests.jl` pins the floor so a change to
it is a deliberate act rather than a silent one.

## Which paths the pool can actually reach

Worth stating plainly, because three of the four plausible answers are wrong.

* **The density callback's batch route** (`density_callbacks/runtime.jl`) is the
  one place the pool evaluates GRAM. It is taken when the batch route is
  selected and the *track* cache is off.
* **The vacuum-predicted look-ahead cache** (`density_callbacks/vacuum_predicted_gram.jl`,
  `SPACEAGORA_VACUUM_GRAM_CACHE`) never reaches the pool *for its own knots*:
  `_build_vacuum_gram_cache!` queries them through scalar `getDensity` on the
  shared model, under the global lock, one satellite at a time. Turning the
  look-ahead cache on does not turn the pool off, though, because the staged
  density callback still takes the batch route — its route selection consults
  the *track* cache, not this one. So in look-ahead mode the pool evaluates the
  callback's staged densities while the RHS reads the spline, and the two are
  produced by different mechanisms. Measured on a 64-spacecraft,
  300-480 km band, 10 s run at 4 threads: enabling the pool at width 4 moved 384 of
  1664 `gram_density` lock acquisitions off the shared lock, and 384 is exactly
  what the same run's freeze-per-step arm takes in total — that is, the
  callback's whole share and nothing else.
* **The GRAM track cache** (`callbacks/gram_track_cache/refresh.jl`,
  `SPACEAGORA_GRAM_TRACK_CACHE`) reaches `_gram_isolated_pool_batch_eval!` only
  on its fallback branch, the one taken when the native driver does not expose
  `generate_trajectory`. The vendored GRAM Suite 2.0 does expose it, so on this
  repository's own GRAM that branch is dead and every refresh is one locked
  native `generate_trajectory` call.
* **The constellation RHS aero path** samples through
  `_density_model_for_sat(p, sat_idx)`, i.e. the per-satellite instance vector
  (`SPACEAGORA_GRAM_PER_SAT_INSTANCES`), not the pool.

## The pool is not the only instance-isolation mechanism in the tree

Three switches share the same premise — that independent native GRAM instances
may be called concurrently as long as each single instance is serialized — and
they apply it in three different places. They are easy to confuse and they do
not compose the way the names suggest.

* `SPACEAGORA_GRAM_ISOLATED_POOL` builds per-*worker* instances inside the
  density callback's batch call, and hands each one its own lock explicitly. It
  therefore ignores `SPACEAGORA_GRAM_LOCK_SCOPE` entirely: the pool's calls
  never reach `_gram_call_lock`, so they are off the shared lock whatever that
  variable says.
* `SPACEAGORA_GRAM_PER_SAT_INSTANCES` builds per-*satellite* instances
  (`_initialize_density_model_instances!` in `simulation/engine/setup.jl`), and
  those are what the constellation RHS aero path actually samples through. They
  are the pool's natural counterpart for the RHS, and they were the path whose
  fresh clones hit CSPICE concurrently; `setup.jl` now gives them the same
  single-threaded warm-up the pool build uses.
* `SPACEAGORA_GRAM_LOCK_SCOPE=model` changes which lock the *scalar* call sites
  take, from the shared one to the wrapper's own `instance_lock`. On its own it
  buys nothing, because one shared model still serializes on its own lock; it is
  only useful in combination with per-satellite or per-sample instances. Its
  occupancy is also deliberately not recorded in the native-lock counters, so a
  run using it reads as having almost no GRAM lock time — which is correct but
  easy to misread as a speedup.

## Two configuration traps in the existing benchmark cases

Both were found while building this study and both change what the P6/S2 GRAM
traces are understood to measure.

1. `ppc_constellation(planet, n)` places member *i* at `540 + 2(i-1)` km
   apoapsis altitude. At `n = 1024` its upper members are above 2000 km, where
   `getDensity(::GRAMAtmosphereModel, ...)` returns zero *without calling GRAM*.
   The fraction of the constellation that touches GRAM therefore falls as the
   size axis is swept, which is the axis a scaling claim is made along.
2. Every member of that constellation starts far above the 120 km entry
   interface `ppc_build_config` sets, so `in_atmosphere` is false for all of them
   and the vacuum-predicted look-ahead cache never builds. The
   `aero_<N>sat_l50_gram_lookahead_<S>s` case therefore does not exercise the
   look-ahead cache; it differs from the `gram_process` case only in
   `SPACEAGORA_DENSITY_FREEZE_PER_STEP`.

`run_scaling_config.jl` builds its own 300–480 km constellation and puts the
entry interface above it for exactly these two reasons.

## Refresh de-phasing: not available

`benchmarks/studies/paper_scenarios/FABLE_FINDINGS.md` D5 proposes staggering
cache horizons per satellite (±10 % jitter) so that expiries de-phase and the
refreshes stop convoying on the global lock.

That cannot be done without changing the dynamics, so it was not done. The
horizon sets the knot spacing directly — `h = horizon_s / (n_pts - 1)` in
`_build_vacuum_gram_cache!`, `Δt = dt_segment / (n - 1)` in the track cache's
refresh — so jittering it changes which `(altitude, latitude, longitude, time)`
points GRAM is asked about, which changes the spline, which changes every
interpolated density the RHS reads. The refresh schedule is part of the
dynamics.

The other two approaches in the same finding are not ruled out by this argument:
batching several satellites' track requests into one locked native call, and
prefetching the next segment before expiry, both leave the sampled points
untouched. Neither is implemented here.
