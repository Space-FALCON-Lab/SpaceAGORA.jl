# Where the third-body and SRP cost goes at constellation scale

## The observation

The P6 calibration of the parallelization figure runs timed one serial solve of
each 4096-spacecraft thread-scaling trace on TRX50 (job `20260922-054218-1797061`,
single thread):

| Case | Mission | Serial wall |
|---|---|---|
| `gravity_4096sat_l50_vacuum_5800s` | 5800 s | 11.71 s |
| `gravity_4096sat_l50_srp_nbody_vacuum_5800s` | 5800 s | 365.04 s |
| `gravity_4096sat_l20_vacuum_19700s` | 19700 s | 11.82 s |

Adding solar radiation pressure and Sun/Moon third-body gravity multiplies the
serial solve by 31 at the same constellation size and mission length, while the
other traces sit within a factor of three of each other. The standing suspicion
was per-spacecraft SPICE ephemeris lookups inside the RHS.

That suspicion is wrong. **The run makes zero runtime SPICE calls.** The cost is
three separate effects stacked on top of each other, and only the smallest of
them has anything to do with ephemerides.

## How it was measured

`benchmarks/studies/third_body_cost/variants.jl` builds the P6 trace's
constellation, tolerances, step cap and mission through the
`parallelization_performance` study's own builders and changes exactly one thing
per rung:

| Rung | Dynamic effectors |
|---|---|
| `vacuum` | degree-50 harmonics |
| `nbody` | + `NBodyGravityModel(("Sun", "Moon"))` |
| `srp` | + `SolarRadiationPressureModel(1.2, 12.0)` |
| `srp_nbody` | + both (the P6 trace) |

`--solver=` forces the solver family, `--rhs=` forces the RHS execution route,
and `--mode=` picks the parallel profile, which is how the three effects below
are separated from one another. All rungs of one invocation are solved back to
back in one process; everything quoted here is a **ratio between rungs of the
same invocation**, never an absolute benchmark — the workstation was shared with
another measurement session throughout.

Unless stated otherwise: `space-falcon-1`, one thread, R0 (`serial`) mode,
N = 256, 5800 s mission, best of three repeats.

The effect reproduces at both sizes this workstation can hold, without needing
4096 spacecraft:

| N | `srp_nbody` / `vacuum` | `nf` ratio | per-evaluation ratio |
|---|---|---|---|
| 256 | 49.7x | 5.58x | 8.9x |
| 1024 | 36.4x | 3.79x | 9.6x |

(TRX50's 4096-spacecraft P6 pair is 31x.) The `nf` ratio moves between sizes
because the point at which `AutoTsit5` commits to Rodas5P depends on the
trajectory mix; the per-evaluation ratio, which is what the RHS controls, does
not.

## What the 31x is made of

### 1. The solver family changes — about 8x

`SolverConfig.solver_mode = :auto_stiff` (the harness default) has a fast path:
`_auto_stiff_smooth_gravity_eligible` in `src/simulation/engine/solver_policy.jl`
routes a configuration whose effectors are all *smooth gravity*
(inverse-square, J2, harmonics, n-body) with no atmosphere and no solar
requirement straight to plain `Tsit5()`. Anything else gets
`AutoTsit5(Rodas5P(autodiff=AutoFiniteDiff(), …))`.

`SolarRadiationPressureModel` is not in that effector list, and it declares
`environment_requirements(...).solar = true`, so adding SRP disqualifies the run
twice over. The autoswitching algorithm then detects stiffness and commits to
Rodas5P:

| Rung | Solver | `nf` | `njacs` | `nw` | accepted | rejected |
|---|---|---|---|---|---|---|
| `vacuum` | Tsit5 | 2074 | 0 | 0 | 296 | 0 |
| `nbody` | Tsit5 | 2074 | 0 | 0 | 296 | 0 |
| `srp` | AutoTsit5(Rodas5P), switched | 11582 | 639 | 639 | 643 | 7 |
| `srp_nbody` | AutoTsit5(Rodas5P), switched | 11582 | 639 | 639 | 645 | 5 |

That is 5.58x the derivative evaluations, plus the W-matrix machinery itself: a
CPU profile of the 5800 s `srp_nbody` solve puts 21116 of 30846 samples (68%)
under the Rosenbrock `perform_step!` and 18777 (61%) under `calc_W!`, with
sparse-matrix coloring and KLU factorizations appearing alongside.

The implicit solver buys nothing here. Forcing `--solver=tsit5` on the same
`srp_nbody` case gives **the same step sequence as the vacuum rung** — `nf` 2074,
296 accepted steps, zero rejections — and a terminal state that agrees with the
Rodas5P answer to 2e-11 relative in position and 6e-11 in velocity:

```
auto_stiff  pos = [-6853023.5257259775, 114073.98673278769, 877877.92669214134]
tsit5       pos = [-6853023.5258762725, 114073.98697414460, 877877.92688973306]
```

The dynamics are not stiff; the eligibility list simply does not know that SRP
is smooth enough for the explicit solver. Cost of the misrouting on this case:
`srp_nbody` is 53.7x the vacuum rung under `auto_stiff` and 6.7x under `tsit5`.

### 2. The RHS route falls off the batched harmonics kernel — about 4x

`_resolve_rhs_execution_plan` in `src/simulation/engine/setup.jl` routes a
single-harmonics constellation through `:flat_constellation_effector_queue`
*even at a thread budget of one* — the branch explicitly allows `budget <= 1`,
because the flat route's harmonics pre-pass loads each (degree, order)
coefficient once and sweeps it across a slice of spacecraft.

Every other multi-effector configuration hits `budget <= 1 → :satellite_batch`
first, and `satellite_batch` evaluates each spacecraft's whole effector chain
independently, so the coefficient table is re-traversed per spacecraft. Adding a
second effector of any kind is therefore enough to lose the batched kernel at
one thread. Measured on the vacuum rung alone, so that the route is the only
variable (`--rhs=satellite` against the default route, best of three each):
1.952 s versus 0.562 s, i.e. 3.5x, in line with the ~4x quoted in that
function's own comment.

The flat route is also where the *shared* environment values are hoisted:
`_prefill_shared_body_samples!` computes the Sun and third-body positions once
per derivative evaluation, and `_accumulate_nbody_flat_batch!` /
`_accumulate_srp_flat_batch!` read them once and sweep the constellation. None
of that runs under `satellite_batch`.

This is why the gap closes when threads are available. Same case, same machine,
8 threads, R2 (`inner_only`) instead of R0:

| Rung | 1 thread, R0 | 8 threads, R2 |
|---|---|---|
| `vacuum` | 1.00x | 1.00x |
| `nbody` | 5.76x | 2.10x |
| `srp_nbody` | 49.7x | 17.5x |

with byte-identical trajectories between the two. The P6 serial baseline is
therefore measuring the route's worst case, not the force model's cost.

### 3. The third-body sample is rebuilt per spacecraft — about 1.5x

Under `satellite_batch`, `sample_environment_with_reusable_buffers`
(`src/simulation/engine/effector_sampling.jl`) reuses a prefilled planet frame,
atmosphere and *solar* sample when one is available, but there is no such reuse
for third bodies: `sample_third_body_ephemerides` runs per spacecraft, per
evaluation, unconditionally.

With the degree-50 harmonics rung as the control, adding the Sun/Moon third body
costs 1.5x per derivative evaluation at the same route and the same step
sequence, and raises the allocation of the whole solve from 0.09 GiB to
3.71 GiB — about 7 KiB per spacecraft per evaluation for two vectors that depend
only on time. A CPU profile of that rung (7606 samples inside the RHS) breaks
down as:

| Path | Samples | Share of RHS |
|---|---|---|
| harmonics effector | 4213 | 55% |
| third-body effector | 2403 | 32% |
| — of which `sample_third_body_ephemerides` | 2191 | 29% |
| — — of which SPICE name canonicalization | 741 | 10% |
| — — of which the ephemeris cache lookup itself | 25 | 0.3% |

The n-body arithmetic is a rounding error; the interpolated cache lookup is a
rounding error; the sample construction around them is not.

**Ephemeris calls per RHS evaluation: zero.** `spice_runtime_counters` reports
`nbody_spkpos_runtime_calls = 0` and `srp_spkpos_runtime_calls = 0` for every
rung above, at every size measured. The pre-solve `NBodyEphemerisCache` and
`SRPSunEphemerisCache` already absorb the SPICE traffic, and the per-spacecraft
work that remains is name resolution and sample construction, not I/O and not a
lock.

## What was changed

`_spice_query_name` (`src/dynamics/coupled/perturbations.jl`) now interns its
result. It is called once per third body per spacecraft per derivative
evaluation from the per-spacecraft sampling path — and once per spacecraft for
the SRP primary body — and `strip`/`lowercase`/`replace` allocated three short
strings every time. The table is published copy-on-write behind a lock, so
readers only ever walk a dictionary that is never mutated afterwards, which is
what makes it safe to call from every RHS worker. The interned value is exactly
the string the uncached expression produces.

Measured before/after, one thread, back to back (N = 256 best of three,
N = 1024 best of two):

| N | Rung | Wall before | Wall after | Allocation before | Allocation after |
|---|---|---|---|---|---|
| 256 | `vacuum` (control) | 0.580 s | 0.562 s | 0.07 GiB | 0.07 GiB |
| 256 | `nbody` | 3.597 s | 3.242 s | 3.71 GiB | 3.36 GiB |
| 256 | `srp_nbody` | 31.168 s | 27.966 s | 30.67 GiB | 26.85 GiB |
| 1024 | `vacuum` (control) | 2.878 s | 2.788 s | 0.30 GiB | 0.30 GiB |
| 1024 | `nbody` | 19.204 s | 19.378 s | 14.83 GiB | 13.41 GiB |
| 1024 | `srp_nbody` | 106.452 s | 101.546 s | 83.09 GiB | 72.72 GiB |

The allocation reduction is deterministic and consistent at both sizes (−9.6%
on `nbody`, −12.5% on `srp_nbody`). The wall-clock gain is clear at N = 256
(−9.9% and −10.3%) and at or below this workstation's noise floor at N = 1024,
where the vacuum control itself moves 3.1% between the two runs and the `nbody`
rung comes back 0.9% *slower*; at that size the solve is dominated by memory
traffic the change does not touch. The vacuum rung never takes the path and is
carried only as a control.
In the profile, the name-resolution frames drop out entirely (741 samples to
below the 20-sample floor, leaving 31 samples of dictionary lookup) and the
third-body effector's share of the RHS falls from 2403 to 1791 samples.

**Bit-identity.** The full state history — every saved step time and every
component of every spacecraft's state — was dumped to a raw `Float64` file
before and after the change for all three rungs and compared byte for byte:

```
vacuum:    BYTE-IDENTICAL (4868424 bytes)
nbody:     BYTE-IDENTICAL (4868424 bytes)
srp_nbody: BYTE-IDENTICAL (10589232 bytes)
```

The N = 1024 runs above report the same terminal states before and after to all
17 digits as well.

`test/unit/dynamics/third_body_cost_tests.jl` pins the two properties this
relies on: the interned name equals the uncached one for a spread of inputs and
costs zero allocations across 4096 resolutions (that assertion fails against the
pre-change implementation, which is what makes it a guard rather than a
description), and the planet-frame lookup asks a counting stand-in ephemerides
backend exactly once per evaluation for a 64-spacecraft constellation, returning
the uncached matrix bit for bit.

## What is still on the table

Three changes outside this file's ownership would take the rest, in descending
order of size:

1. **Let SRP take the explicit solver.** Adding `SolarRadiationPressureModel` to
   `_auto_stiff_smooth_gravity_effector` and allowing `req.solar` for it in
   `_auto_stiff_smooth_gravity_reject_reason` (`solver_policy.jl`) removes the
   5.58x in derivative evaluations and the whole Rosenbrock W path. The
   trajectory evidence above says the two solvers agree to 1e-10 relative on this
   case; a configuration that genuinely needs the implicit solver still has
   `solver_mode=:auto_stiff` available through the eligibility switch
   `SolverConfig.auto_stiff_gravity_tsit5`. The one thing to re-check before
   applying it: `eclipse_area_calc` is a piecewise conical-shadow model, so the
   SRP acceleration has a kink at umbra and penumbra boundaries. On this case
   (500-550 km LEO, `reltol = abstol = 1e-9`, 20 s step cap) Tsit5 crossed them
   with zero rejected steps, but a much larger step cap deserves its own look.
2. **Reach the batched route at one thread.** The `budget <= 1 → satellite_batch`
   guard in `setup.jl` predates the batchable-effector pre-passes; the flat route
   at `allotment = 1` spawns no tasks, which is exactly why the single-harmonics
   branch is already allowed to take it.
3. **Prefill the third-body sample like the solar one.** `_prefill_shared_body_samples!`
   already computes it once per evaluation and throws the result away; giving it
   an `rhs_third_body_*` buffer pair and an `_sample_reusable_third_bodies` reader
   would make the per-spacecraft path a buffer read, the same shape as
   `_sample_reusable_solar`.

## The three follow-ups, applied and measured

All three items from the list above were implemented on
`pv2c-ws10c-thirdbody-engine`. Each was measured on its own: the ladder was run
back to back before and after each change separately, one thread, same process
state, N = 256 (best of three) and N = 1024 (best of two), on `space-falcon-1`
while it was shared with another measurement session -- so every number below is
a ratio between two runs of the same ladder, never an absolute benchmark.

### What changed

1. `sample_third_body_ephemerides` (`src/simulation/engine/effector_sampling.jl`)
   builds its position tuple with `map` over `model.body_names` instead of
   `ntuple` over `length(model.body_names)`. The body-name tuple's length is part
   of `NBodyGravityModel`'s type, so the mapped construction infers; the
   `ntuple`-over-a-runtime-`Int` form did not, and allocated the tuple and its
   boxed closure once per spacecraft per derivative evaluation. Same bodies, same
   order, same values.
2. `_rhs_execution_plan_uncached` (`src/simulation/engine/setup.jl`) admits a
   whole class of stacks to the flat route at `budget <= 1`, not just the
   single-harmonics one: if every effector is served by a flat-route pre-pass --
   `_batchable_effector` (n-body, SRP, non-gradient inverse-square) or
   `_harmonics_prepass_effector` -- and the constellation clears
   `env.flat_min_sats`, the route is `:flat_constellation_effector_queue` at
   allotment 1 with the serial effector decision. The new predicate is
   `_rhs_all_prepass_effectors`, and it dispatches on the same two traits the
   flat driver itself uses, so it cannot drift from what that driver pre-passes.
3. `_auto_stiff_smooth_gravity_effector` / `_auto_stiff_smooth_gravity_reject_reason`
   (`src/simulation/engine/solver_policy.jl`) accept `SolarRadiationPressureModel`
   as a smooth effector and no longer reject it for declaring `req.solar` -- the
   solar sample is the Sun position that SRP itself consumes. Any other
   smooth-gravity effector asking for solar samples still rejects, and
   `SolverConfig.auto_stiff_gravity_tsit5=false`
   (`SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5=0`) still disables the whole fast path
   for a configuration that wants the implicit solver.

### Why the route change is safe at one thread

The flat route at allotment 1 spawns no tasks, and for a pre-pass-only stack it
does not even build its work queue:

- `_accumulate_dynamic_effectors_flat_slots!` runs the batchable pre-pass and
  the harmonics pre-pass and returns at `_count_flat_queue_only_effectors(...) == 0`,
  before `_build_constellation_execution_plan!`.
- `_accumulate_harmonics_flat_batch!` caps its worker count through
  `thread_worker_count(active_sats, 1) == 1` and runs the slice inline.
- `_reduce_flat_effector_slots!` and the final per-satellite assembly
  (`threaded_foreach(..., plan.allotment)`) both take their serial branch at one
  worker.

`_rhs_flat_supported(env, ...)` is still required, so the effector set must be
one the flat route supports at all; `outer_active` needs no separate guard here,
because the clamp the other branches apply to protect an enclosing outer split
-- allotment 1 -- is this branch's only allotment.

### Deliverable 1: the type-stable third-body sample

| N | Rung | Wall before | Wall after | Alloc before | Alloc after |
|---|---|---|---|---|---|
| 256 | `vacuum` (control) | 0.566 s | 0.583 s | 0.07 GiB | 0.07 GiB |
| 256 | `nbody` | 3.748 s | 3.367 s | 3.36 GiB | 3.09 GiB |
| 256 | `srp_nbody` | 29.254 s | 27.323 s | 26.85 GiB | 24.75 GiB |
| 1024 | `vacuum` (control) | 2.890 s | 2.788 s | 0.30 GiB | 0.30 GiB |
| 1024 | `nbody` | 18.580 s | 18.497 s | 13.41 GiB | 12.37 GiB |
| 1024 | `srp_nbody` | 98.493 s | 98.301 s | 72.72 GiB | 67.02 GiB |

The allocation reduction is the measurement that means something: -7.8% of the
whole solve on every rung that takes the path, at both sizes (3.09/3.36,
24.75/26.85, 12.37/13.41, 67.02/72.72 -- all 0.921 to 0.922). Divided by the
work, that is about 0.5 KiB per spacecraft per derivative evaluation at both
sizes (546 B at N = 256, 526 B at N = 1024, derived from the printed GiB totals
and `nf` x N). The wall-clock effect is 0.898x on `nbody` and 0.934x on
`srp_nbody` at N = 256 and inside the noise floor at N = 1024, where the vacuum
control itself moves 0.965x between the two runs.

**Bit-identity.** The full state history of all three rungs at N = 256 was
dumped before and after (`variants.jl --dump`) and compared byte for byte:

```
vacuum:    BYTE-IDENTICAL (4868424 bytes)
nbody:     BYTE-IDENTICAL (4868424 bytes)
srp_nbody: BYTE-IDENTICAL (10589232 bytes)
```

The N = 1024 terminal states agree to all 17 printed digits as well.

### Deliverable 2: the RHS route at a thread budget of one

| N | Rung | Wall before | Wall after | Ratio | Alloc before | Alloc after |
|---|---|---|---|---|---|---|
| 256 | `vacuum` (control) | 0.583 s | 0.562 s | 0.96x | 0.07 GiB | 0.07 GiB |
| 256 | `nbody` | 3.367 s | 1.242 s | 0.37x | 3.09 GiB | 1.26 GiB |
| 256 | `srp_nbody` | 27.323 s | 6.779 s | 0.25x | 24.75 GiB | 13.89 GiB |
| 1024 | `vacuum` (control) | 2.788 s | 2.871 s | 1.03x | 0.30 GiB | 0.30 GiB |
| 1024 | `nbody` | 18.497 s | 4.100 s | 0.22x | 12.37 GiB | 4.95 GiB |
| 1024 | `srp_nbody` | 98.301 s | 21.741 s | 0.22x | 67.02 GiB | 36.95 GiB |

The vacuum rung is single-harmonics and already took the flat route, which is
why it is carried here as a control and does not move. Everything else is the
batched kernel and the shared body samples coming back: 2.7x on `nbody` and
4.0x on `srp_nbody` at N = 256, 4.5x on both at N = 1024, at an unchanged step
sequence (`nf`, accepted and rejected steps identical before and after).

**Parity.** The flat route and the per-satellite route reach the same sums by
different paths, so this is the change where parity had to be shown rather than
argued. Three independent checks:

```
dump/cmp, N = 256, before vs after:
  vacuum:    BYTE-IDENTICAL (4868424 bytes)
  nbody:     BYTE-IDENTICAL (4868424 bytes)
  srp_nbody: BYTE-IDENTICAL (10589232 bytes)
```

- `test/unit/simulation/harmonics_batch_parity_tests.jl` passes unchanged (the
  batch kernel still reproduces the scalar kernel bit for bit).
- The `parallelization_performance` harness's trajectory parity, run before and
  after on `stack32_e2_srp` (harmonics + SRP, 32 spacecraft: a stack the new
  branch admits) and `stack32_e4_nbody` (which carries aerodynamics and
  therefore still takes the per-satellite route), serial reference against
  `full_smart` at two threads, 128 sampled states: `pass=true` with
  `pos_rel_max = vel_rel_max = 0.0` in both runs.

The N = 1024 terminal states are also identical to all 17 printed digits before
and after.

### Deliverable 3 (measurement) and deliverable 4: SRP on the explicit solver

| N | Rung | Wall before | Wall after | Ratio |
|---|---|---|---|---|
| 256 | `vacuum` (control) | 0.562 s | 0.567 s | 1.01x |
| 256 | `nbody` (control) | 1.242 s | 1.249 s | 1.01x |
| 256 | `srp_nbody` | 6.779 s | 0.829 s | 0.12x |
| 1024 | `vacuum` (control) | 2.871 s | 2.953 s | 1.03x |
| 1024 | `nbody` (control) | 4.100 s | 4.175 s | 1.02x |
| 1024 | `srp_nbody` | 21.741 s | 4.451 s | 0.20x |

Only the SRP-carrying rung changes, and it changes because the solver does:
`nf` falls from 11582 to 2074 at N = 256 (7856 to 2074 at N = 1024) and the
step sequence becomes the vacuum rung's -- 296 accepted, 0 rejected.

**This one is not bit-identical, by construction:** it is a different
integrator. The two answers were compared directly on the 256-spacecraft
`srp_nbody` rung, in one process, with the fast path disabled
(`SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5=0`, i.e. AutoTsit5 switching to Rodas5P)
and enabled (Tsit5), at the default 20 s step cap and again at 60 s, because
`eclipse_area_calc` is a piecewise conical-shadow model and the SRP
acceleration has a kink at the umbra and penumbra boundaries:

| Step cap | Solver | Wall | `nf` | `njacs` | accepted | rejected |
|---|---|---|---|---|---|---|
| 20 s | AutoTsit5(Rodas5P), switched | 7.312 s | 11582 | 639 | 645 | 5 |
| 20 s | Tsit5 | 0.844 s | 2074 | 0 | 296 | 0 |
| 60 s | AutoTsit5(Rodas5P), switched | 7.211 s | 11474 | 633 | 638 | 6 |
| 60 s | Tsit5 | 0.594 s | 1500 | 0 | 214 | 0 |

Maximum difference over all 256 spacecraft's terminal states, explicit against
implicit:

| Step cap | Position | Velocity |
|---|---|---|
| 20 s | 5.33e-10 (relative to the state's norm), 7.61e-08 worst single component | 5.02e-10, 9.76e-08 worst component |
| 60 s | 1.54e-09 (relative to the norm), 2.90e-07 worst component | 1.46e-09, 2.93e-07 worst component |

The component-wise figure is the larger of the two only because individual
position and velocity components pass through zero; the norm-relative figure is
the meaningful one. Tripling the step cap does not trip the shadow kink on this
case: Tsit5 still takes zero rejected steps at 60 s, while the implicit solver
takes six. Both comparisons are at `reltol = abstol = 1e-9` in 500-550 km LEO.

### Where the 4096-spacecraft case stands

Cumulatively, on the ladder this file has used throughout (one thread, back to
back, same process):

| N | Rung | Before all three | After all three | Ratio |
|---|---|---|---|---|
| 256 | `srp_nbody` | 29.254 s | 0.829 s | 35.3x faster |
| 256 | `nbody` | 3.748 s | 1.249 s | 3.0x faster |
| 1024 | `srp_nbody` | 98.493 s | 4.451 s | 22.1x faster |
| 1024 | `nbody` | 18.580 s | 4.175 s | 4.4x faster |

The quantity the original observation was about -- `srp_nbody` against `vacuum`
at the same size and mission -- falls from 51.7x to 1.46x at N = 256 and from
34.1x to 1.51x at N = 1024. Adding SRP and Sun/Moon third-body gravity to a
degree-50 constellation now costs about half as much again as the harmonics
alone, which is what the force model is worth.

With threads the remaining gap is smaller still. Same case, 8 threads, R2
(`inner_only`), N = 256, after all three changes: `vacuum` 0.294 s, `nbody`
0.530 s (1.80x), `srp_nbody` 0.556 s (1.89x), against 2.10x and 17.5x for the
same rungs before this work.

### Still open

The third item on the list above -- giving `_prefill_shared_body_samples!` an
`rhs_third_body_*` buffer pair so the per-spacecraft path is a buffer read --
was *not* implemented, and the route change makes it matter less than it did: a
pre-pass-only stack never calls `sample_third_body_ephemerides` at all now,
because `_accumulate_nbody_flat_batch!` gathers the body positions once per
evaluation. It still applies to stacks that keep the per-satellite route, such
as anything carrying aerodynamics, where the mapped sample of deliverable 1 is
the only improvement that path has had.

## Reproducing

```bash
# the attribution ladder (ratios only; one Julia process at a time)
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl \
    --n=256 --mission=5800 --repeats=3

# the same ladder with the solver or the route pinned
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl \
    --n=256 --mission=5800 --variants=vacuum,srp_nbody --solver=tsit5
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl \
    --n=256 --mission=5800 --variants=vacuum,nbody --rhs=satellite

# with threads, through the inner-parallel profile
julia --project=. --threads=8 benchmarks/studies/third_body_cost/variants.jl \
    --n=256 --mission=5800 --mode=inner_only --repeats=2

# the P6 pair as the calibration runs it, plus a CPU profile
julia --project=. --threads=1 benchmarks/studies/third_body_cost/attribution.jl \
    --n=256 --mission=5800 --profile
```

```bash
# the same rung on the implicit solver, for the explicit-vs-implicit comparison:
# the env var disables the auto-stiff fast path that now admits SRP
SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5=0 julia --project=. --threads=1 \
    benchmarks/studies/third_body_cost/variants.jl \
    --n=256 --mission=5800 --variants=vacuum,srp_nbody
```

`--dump=<prefix>` writes each rung's full state history as raw `Float64` for a
byte-for-byte trajectory comparison across two invocations.

The profiles quoted above are committed next to the scripts:
`profile_tree_srp_nbody_256sat_5800s_before.txt` and
`profile_flat_srp_nbody_256sat_5800s_before.txt` are the P6 case as the
calibration runs it, and the `*_auto_stiff_satellite_{before,after}.txt` pairs
are the `vacuum` and `nbody` rungs on the per-satellite route, which is where
the name-resolution share is read off.
