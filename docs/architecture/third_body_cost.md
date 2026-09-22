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

`--dump=<prefix>` writes each rung's full state history as raw `Float64` for a
byte-for-byte trajectory comparison across two invocations.

The profiles quoted above are committed next to the scripts:
`profile_tree_srp_nbody_256sat_5800s_before.txt` and
`profile_flat_srp_nbody_256sat_5800s_before.txt` are the P6 case as the
calibration runs it, and the `*_auto_stiff_satellite_{before,after}.txt` pairs
are the `vacuum` and `nbody` rungs on the per-satellite route, which is where
the name-resolution share is read off.
