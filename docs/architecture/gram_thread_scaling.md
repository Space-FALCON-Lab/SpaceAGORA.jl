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

At 1024 spacecraft and four threads or more, yes, by up to 1.90x. At 256 it
loses at every thread count and every width, by as much as a factor of two. The
axis that decides is how many native GRAM calls a single callback makes, and the
threshold sits between those two sizes.

All of it on the workstation (12 physical cores, 24 threads), one thread count
per process, each group solved back to back in one process state with the arms
alternating across three repeats and the minimum taken. Ratios only:
locked ÷ pool, so above 1 means the pool is faster. Lock hold and wait are the
`gram_density` site's, over the locked run. Raw rows in
`benchmarks/studies/gram_thread_scaling/results/`.

| threads | spacecraft | density path | locked (s) | pool 2 | pool 4 | pool 8 | locked lock hold (s) | locked lock wait (s) | pool 8 acquisitions |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 256 | freeze | 0.42 | 1.00 | 1.00 | 0.99 | 0.16 | 0.00 | 6144 |
| 1 | 256 | lookahead | 0.66 | 1.01 | 0.99 | 0.99 | 0.33 | 0.00 | 11520 |
| 1 | 1024 | freeze | 1.81 | 0.99 | 1.01 | 0.99 | 0.55 | 0.00 | 24576 |
| 1 | 1024 | lookahead | 3.23 | 1.15 | 1.01 | 1.15 | 1.24 | 0.00 | 46080 |
| 2 | 256 | freeze | 0.31 | 0.90 | 0.91 | 0.84 | 0.16 | 0.00 | 0 |
| 2 | 256 | lookahead | 0.54 | 0.86 | 0.89 | 0.89 | 0.36 | 0.10 | 5120 |
| 2 | 1024 | freeze | 1.20 | 0.96 | 0.93 | 1.02 | 0.54 | 0.00 | 0 |
| 2 | 1024 | lookahead | 2.05 | 1.06 | 0.96 | 0.85 | 1.24 | 0.37 | 20480 |
| 4 | 256 | freeze | 0.25 | 0.88 | 0.67 | 0.71 | 0.17 | 0.00 | 0 |
| 4 | 256 | lookahead | 0.51 | 0.92 | 0.81 | 0.80 | 0.36 | 0.33 | 5120 |
| 4 | 1024 | freeze | 0.83 | 1.37 | 1.65 | 1.41 | 0.55 | 0.00 | 0 |
| 4 | 1024 | lookahead | 1.76 | 1.17 | 1.29 | 1.10 | 1.35 | 1.34 | 20480 |
| 8 | 256 | freeze | 0.21 | 0.83 | 0.69 | 0.49 | 0.16 | 0.00 | 0 |
| 8 | 256 | lookahead | 0.52 | 0.90 | 0.79 | 0.68 | 0.42 | 1.78 | 5120 |
| 8 | 1024 | freeze | 0.88 | 1.64 | 1.90 | 1.76 | 0.55 | 0.00 | 0 |
| 8 | 1024 | lookahead | 1.76 | 1.13 | 1.27 | 1.20 | 1.32 | 3.76 | 20480 |

The one-thread block is the null control and reads as one. The pool needs at
least two workers, so at one thread it declines and both arms are the same code:
the acquisition counts are identical to the unit, and the ratios come back within
1 % of unity — except the 1024-spacecraft look-ahead row, whose wall time is
bimodal between about 2.80 s and 3.23 s in *both* arms. That 15 % spread is the
resolution floor for that one cell, and nothing smaller should be read out of it.

Above one thread the pool engages, and the acquisition counts say exactly how
far it reaches. In freeze-per-step the shared lock's `gram_density` count goes to
zero: every native GRAM call has moved to a per-worker instance behind a
per-worker lock. In look-ahead mode it falls from 11 520 to 5 120 and from 46 080
to 20 480 — the callback's share moves, and the look-ahead cache's own knot
queries stay on the shared lock, because `_build_vacuum_gram_cache!` runs them
scalar on the shared model.

### The lock was not the problem where the pool helps

The most useful column is the one that is zero. In freeze-per-step the locked
path's lock *wait* is 0.00 s at every thread count and both sizes, against a hold
of 0.16–0.55 s. There is no queueing: `getDensityBatch!` evaluates the whole
batch from the one thread that entered the callback, so the shared lock is taken
6 144 or 24 576 times and never contended. What the pool buys there is not lock
relief, it is parallelism over work that was simply serial.

And where the lock *is* contended, the pool does not reach it. In look-ahead mode
the locked wait grows with the thread count — 0.10, 0.33, 1.78 s at 256
spacecraft, 0.37, 1.34, 3.76 s at 1024 — and the pooled arms carry the same wait
or more (up to 4.39 s), because the calls doing the waiting are the look-ahead
cache's, which the pool never touches. That is the honest reading of the
look-ahead column: its 1.13–1.29x is the callback's share being parallelized
while the contended path underneath is unchanged.

### Why 256 loses, and it is not the build cost

Every pooled solve constructs its own instances, so the obvious explanation for
the 256-spacecraft loss is a fixed build that a 100 s mission cannot amortize.
The mission-length control says otherwise:

| mission (s) | spacecraft | density path | threads | locked (s) | pool 4 |
|---:|---:|---|---:|---:|---:|
| 100 | 256 | freeze | 8 | 0.24 | 0.68 |
| 100 | 256 | lookahead | 8 | 0.50 | 0.75 |
| 1000 | 256 | freeze | 8 | 1.64 | 0.80 |
| 1000 | 256 | lookahead | 8 | 2.08 | 0.88 |

Ten times the mission moves the ratio from 0.68 to 0.80 and from 0.75 to 0.88 —
better, and still a loss. Meanwhile four times the constellation moves it from
0.69 to 1.90 at the same thread count. So the governing quantity is the native
GRAM work inside one callback invocation, not the number of invocations: the
threaded fan-out has a fixed per-invocation cost, and at 256 spacecraft the
0.17 s of native GRAM spread across the whole run is not enough to cover it.

### What ships

`SPACEAGORA_GRAM_ISOLATED_POOL` now defaults to `auto` rather than `off`, at a
threshold of 1024 and a width capped at 4. Each of those three numbers is
SOURCED from the table above and argued where it is defined, in
`density_callbacks/config.jl`: 1024 because 256 loses in every cell measured and
1024 wins at 4 and 8 threads in both density paths; 4 because it is the fastest
width in all four winning cells and each further instance is another native GRAM
image resident in the process.

Two changes were needed to make that default mean anything.

The first is the width. The pool asks for its width with `lock_free=true`, which
routes it past the `:density_callback` source's 16-thread minimum budget — the
floor described in the next section, which exists because native GRAM is
serialized on the shared lock and therefore does not apply to workers that each
hold their own instance. Without this the default would be inert on any machine
with fewer than 16 threads, including the one it was measured on.

The second is a guard, and it came out of measuring a case the sweep does not
cover. A constellation that is *configured* with a GRAM atmosphere but never
reaches it — 1024 spacecraft above the entry interface on a non-keplerian run —
was **1.83x slower** with the pool on (0.19 s locked against 0.35 s pooled at 8
threads), because `_ensure_gram_isolated_pool!` builds its four native GRAM
models before the per-item gate ever runs and then nothing calls them.
`_gram_isolated_pool_native_count` now counts, in one pass over the staged
altitudes, how many items would really reach GRAM, and the threshold is applied
to that rather than to the spacecraft count. With the guard that case is 0.96x,
within the run-to-run spread, and the 1024-spacecraft in-atmosphere win is
unchanged at 1.77x and 1.19x.

The defaults are bit-identical to the old ones where they change behavior, which
is the only claim that matters here: a 1024-spacecraft freeze-per-step run and a
1024-spacecraft look-ahead run, each dumped with the pool explicitly off and then
with nothing set at all, are byte for byte the same (1 573 056 bytes each), and
so is the 64-spacecraft reference case, which the new threshold leaves on the
locked path.

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
