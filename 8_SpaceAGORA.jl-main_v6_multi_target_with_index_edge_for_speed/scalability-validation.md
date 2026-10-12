# V6 interlink optimization and validation

Date: 2026-10-05. Julia 1.12.1, Linux, CPU reported by Julia as `znver5`.
Comparisons ran sequentially with one Julia thread. See [plan.md](plan.md) for
requirements and implementation tasks. Legacy, v5 and v5 copy source files were
not modified; no dependencies were added.

## What changed

1. [Interlink storage](I_src/1.4_satellite_interlink/inter_link_models.jl) retains
   dictionary lookup and the same mutable connection states, while adding indexed
   keys/connections, score and availability vectors, signed availability transitions,
   and per-spacecraft incident-edge lists. A bulk availability overload updates in
   place; the single-key overload remains. Position differences use `SVector`.
2. [Scoring policies](I_src/1.4_satellite_interlink/scheduling_policies.jl) iterate
   indexed connections. The three existing objectives, including first-index SMA,
   are unchanged. All candidates are rescored; unchanged availability is not used
   to skip scoring or selection.
3. [Force evaluation](I_src/1.3_dynamics/coupled/force_torque_models/laser_force_effectors.jl)
   visits only incident links. Across all N spacecraft this changes edge inspections
   from O(N*E) to O(N+E), retaining continuous stage-dependent forces and the same law.
4. [Matching](I_src/1.4_satellite_interlink/scheduling_models.jl) caches sorted
   terminal numbering, a 2-by-E integer endpoint matrix, adjacency and work buffers.
   Iterative DFS identifies components of the available, positive-score graph.
   Trees use bottom-up dynamic programming and iterative reconstruction; cyclic
   components are passed to NetworkX Blossom. No exponential subset fallback exists.
   Both solvers compare exact BigInt weights scaled from Float64 scores.
5. Existing tests were extended in [test/runtests.jl](test/runtests.jl); one
   [phase benchmark](test/benchmark_interlinks.jl) was added. The existing
   [Case 1 example](II_examples/compare_pdf_case1.jl) now supports all three version
   labels and first-run/five-warmed-run timings. The new
  [Walker example](II_examples/walker_interlinks.jl) defaults to snapshots; the
  subsequent force-coupled experiment below uses its separate `engine` mode.

The tree fast path was added only after the baseline showed 6.835 seconds for
matching a 4,000-terminal chain. Tree matching is linear in component size in
arithmetic operations; BigInt work has an additional bit-width cost. Dense cyclic
matching remains polynomial Blossom work, not a new dense-graph algorithm.

### API and numerical constraints

- Change topology through `register_candidate!`, which invalidates the matching
  cache. Do not directly insert, delete or replace dictionary entries or edit
  indexed/cache arrays. Dictionary lookups and edits to existing connection state
  still work; direct score/availability state changes are reread by selection.
- `available[e]` is the last availability refresh; `availability_transition[e]`
  is new minus old: +1 entered, -1 left, 0 unchanged. Values are overwritten, not
  accumulated. Single-key refresh affects only that edge. Score mirrors refresh
  during scoring/selection. These are not propagation matrices.
- Cached endpoints identify spacecraft-terminal pairs, so different terminals on
  one spacecraft remain distinct matching vertices.
- Equal-optimum tree selections may differ from Blossom-only tie choices. Sorted
  traversal makes tree ties deterministic, including reversed registration order.
- No changes were made to accepted-step callbacks, timestep limits, laser physics,
  policy selection defaults, or output cadence. No one-minute simulation callback
  was introduced. Existing v6 range, battery, temperature and activity checks remain.
- The copy's optional eligibility hooks and active-link penalties were not ported:
  they are not part of the current v6 API and are not needed for these optimizations.

## Phase benchmark: before and after

[Raw baseline](III_output/scalability/interlinks_baseline.csv),
[optimized measurements](III_output/scalability/interlinks_optimized.csv), and
[complete comparison](III_output/scalability/interlinks_comparison.csv).

The baseline was captured from v6 before these production edits, not from v5 copy.
Both measurements use five warmed calls per stage, pre-call GC, the same synthetic
positions, velocities, masses and `:gve_sma_first_index` scores. Fixture creation,
loading, JIT warm-up, validation and CSV writing are excluded. Force timing covers
all spacecraft. These graph cases are synthetic and separate from the Walker case.

| Graph | N / E | Stage | Before (ms) | After (ms) |
|---|---|---|---:|---:|
| Star | 4000 / 3999 | Availability | 0.51255 | 0.06701 |
| Star | 4000 / 3999 | Scoring | 0.39003 | 0.08639 |
| Star | 4000 / 3999 | Matching | 12.5537 | 1.55547 |
| Star | 4000 / 3999 | All-spacecraft forces | 36.0703 | 0.03879 |
| Chain | 4000 / 3999 | Matching | 6835.056 | 1.80371 |
| Chain | 4000 / 3999 | All-spacecraft forces | 40.8082 | 0.05728 |
| Two-terminal chain | 400 / 399 | Matching | 178.641 | 0.25789 |
| Dense | 100 / 4950 | Matching | 199.582 | 198.395 |
| Dense | 200 / 19900 | Matching | 1600.461 | 1611.380 |
| Dense | 200 / 19900 | All-spacecraft forces | 25.4056 | 0.16977 |

Availability allocated zero Julia bytes in every optimized case, versus 319,920
bytes for the original 4,000-node star and 1,592,000 bytes for the dense 200-node
case. Scoring and forces allocated zero in both versions. Matching still allocates:
exact Rational/BigInt arithmetic and returned selections are not allocation-free.
Star matching fell from 5,533,064 to 3,316,056 Julia bytes; dense 200-node matching
fell from 26,136,912 to 19,946,656 bytes. Julia allocation counts do not include all
Python-side allocations and are not peak resident memory.

Dense matching showed no meaningful runtime improvement (the 200-node measurement
was about 0.7% slower). Sparse-stage speedups are not end-to-end propagation speedups.

## PDF Case 1: legacy, v5 and v6

The shared v6 example was executed with three different `--project` arguments;
each summary records the actual loaded package path. No environment activation
inside that example overrides the selected project.

Configuration: one target at 1050 km, 20 helpers at 1000 km, 227 kg each, J2 gravity,
10,000 W / B=100 / 200 km links, target-only `:gve_sma`, Tsit5, 1e-12 tolerances,
10-second maximum steps, and ten initial target periods (63,713.403604987354 s).
Each version ran laser-on and laser-off once cold and five more times warm in the
same fresh process. Every run starts from a newly constructed configuration.

| Version | Laser-on cold (s) | Laser-on warm median (s) | Laser-off warm median (s) | Final target delta SMA (m) | Gain over laser-off (m) |
|---|---:|---:|---:|---:|---:|
| Legacy | 19.5499 | 0.506021 | 0.849407 | 126.233850247 | 126.173751447 |
| v5 | 32.1523 | 0.908066 | 0.834085 | 126.349023105 | 126.288924257 |
| v6 | 31.6095 | 0.832691 | 0.842530 | 126.349023105 | 126.288924257 |

V6's median laser-on time was **8.30% lower than v5** (1.091x speedup).
The warmed ranges were 0.8900-0.9724 s for v5 and 0.8223-0.8842 s for v6.
Median laser-on Julia allocations decreased from 476,387,472 to 456,311,296 bytes.
Laser-off timing was essentially unchanged and had outliers around 1.25 s.
Legacy was faster here, but it is a different implementation, not an exact
numerical reference that v6 must reproduce. First-use compilation dominates cold
times. Timings cover the measured simulation call, including its setup/callbacks;
they exclude outer configuration creation, initial sampling setup and CSV export.

Numerical checks after reading saved outputs:

- **v5 and v6 match exactly** at all 6,373 sampled times for all 21 spacecraft,
  for both laser-on and laser-off SMA histories.
- Their complete laser schedules are identical: 7,642 rows, including initialization,
  with 7,641 accepted steps and no rejected steps.
- Legacy has 7,968 accepted laser-on steps. Its maximum laser-on SMA difference
  from v6 over all sampled spacecraft/times is 0.247767854 m; the final target
  difference is 0.115172857 m. These are measured differences, not certified errors.
- Laser-off final target change is 0.060098800 m for legacy and 0.060098847 m for
  v5/v6; the maximum all-spacecraft difference is 2.44e-7 m.
- Every run reached its requested endpoint, produced finite state and consistent
  saved samples, and passed 11 assertions: 18 comparisons, 198 assertions total.
  Cold and warm repeats produced identical final target SMA values.

Saved summaries: [version_comparison.csv](III_output/pdf_case1_comparison/version_comparison.csv)
and [numerical_differences.csv](III_output/pdf_case1_comparison/numerical_differences.csv).
The same directory contains per-version timing rows, SMA histories, target
laser-on/off comparisons, schedules and detailed run summaries. Numerical histories
are from the last warmed pass; the timing files retain all six passes.

## Walker snapshot configuration

| Setting | Value |
|---|---|
| Walker-Delta T/P/F | 4000 / 80 / 25 |
| Satellites per plane | 50 |
| Seed | 20261003 |
| Circular altitude / inclination | 550 km / 53 degrees |
| RAAN / anomaly offsets | 292.43964213039555 / 129.40018559361184 degrees |
| Target / terminals per satellite | Vector index 1 / 1 |
| Candidate links | 3,999 target-to-helper links |
| Experiment range / default model range | 2,000 km / 200 km |
| Policy | Target-only `:gve_sma` |
| Scheduling times | 0, 60, 120, ..., 3600 s (61 snapshots) |

The generator preserves Walker spacing and randomizes the phasing and common
offsets with a seeded MersenneTwister. Each snapshot is evaluated analytically
from the original circular two-body state. Laser selection does not modify later
snapshots. There is **no ODE integration or one-minute callback** in this snapshot
experiment. The subsequently requested 20-second full integration is reported
separately below.

### Walker results

| Measurement | Observed |
|---|---:|
| Spacecraft / saved initial-element rows | 4,000 / 4,000 |
| Available links per snapshot | 64-138 |
| Selected links per snapshot | 1 at all 61 snapshots |
| Distinct selected partners | 13 |
| Median availability update | 0.033222 ms |
| Median scoring | 0.007564 ms |
| Median matching | 0.063008 ms |
| Median total scheduling | 0.103805 ms |
| Maximum total scheduling | 0.215102 ms |
| Median scheduling allocation | 40,352 Julia bytes |
| Entered / left availability after warm-up | 525 / 587 |

One schedule call warms up the implementation before measurement. The time-zero
transition vector is therefore zero. Measurements exclude analytic propagation,
validation, generation and exports. Each snapshot checks finite states, terminal
exclusivity, active/available consistency and equality with the maximum positive
score of the target-centered star. Saved CSVs were read back and passed 16 checks.

Outputs in [III_output/walker_interlinks/seed_20261003](III_output/walker_interlinks/seed_20261003):

- `constellation.csv`: generator metadata.
- `spacecraft.csv`: all 4,000 initial elements, plane/slot indices, masses and terminal counts.
- `snapshot_metrics.csv`: 61 rows of timings, allocations, counts and signed-transition counts.
- `snapshot_schedule.csv`: the selected endpoints at all 61 times.
- `snapshot_summary.csv`: aggregate measurements and counts.

This is a sparse single-target workload, not all-pairs matching. Omitting
non-target links is valid only for this fixed-target objective; the general
first-index policy must retain its relevant candidate pairs. One selected link
is the correct result with one target terminal. No laser-induced orbital gain
or long-duration force-coupled scalability is claimed from these snapshots.

## Follow-up: 4,000-satellite force-coupled integration

On 2026-10-05, the same seeded Walker 4000/80/25 constellation was propagated
for **20 seconds with the full simulation engine**, not analytic snapshots.
All 4,000 spacecraft have integrated positions and velocities. The unchanged
laser RHS applies forces to velocity derivatives at each solver stage and
integrates laser delta-v and SMA diagnostics. Scheduling runs at initialization
and every accepted integration step; selected links are held between those steps.

The scenario retains 550 km altitude, 53-degree inclination, seed 20261003,
227 kg spacecraft, 10,000 W / B=100 / 2,000 km links, one terminal per spacecraft,
3,999 target-to-helper candidates and target-only `:gve_sma`. It uses two-body
gravity, no atmosphere, no attitude simulation, Tsit5, relative tolerance 1e-9,
absolute tolerance 1e-11 and a 5-second maximum step. This is not the J2 Case 1
scenario or a one-hour integration.

Only the existing Walker example and regression tests needed code changes:
`run_walker_engine` constructs and validates the full-engine run;
`run_walker_engine_case` exports it through `engine [duration_seconds]`.
The default snapshot invocation and production dynamics/scheduler code are unchanged.

| Measurement | Observed |
|---|---:|
| Spacecraft / candidate links | 4,000 / 3,999 |
| Integrated time / return code | 20 s / Success |
| Accepted / rejected steps | 7 / 0 |
| Schedule samples, including initialization | 8 |
| Active links | 1 at all 8 scheduling times |
| Selected partner vector indices | 50 and 1183 |
| Physical target SMA change | 0.535945602693 m |
| Integrated target laser SMA change | 0.535945602582 m |
| Norm of integrated target laser delta-v | 0.000293820024327 m/s |
| Maximum SMA/laser-diagnostic residual over all saved spacecraft states | 7.45058e-9 m |
| Measured simulation-call wall time | 27.532231 s |
| Julia allocated bytes | 3,528,234,320 |
| GC time inside measured call | 0.501410 s |
| Exported trajectory rows | 32,000 |

The timing is one fresh-process run with one Julia thread and **includes first-use
compilation**. It excludes package loading, outer constellation/initial-state
construction, validation and CSV writing. No warmed throughput claim is made.
Allocated bytes are cumulative Julia allocations, including compilation, not
peak memory. The delta-v above is the norm of the integrated laser acceleration,
not the total velocity change including gravity.

Saved results:

- [constellation.csv](III_output/walker_interlinks/seed_20261003/force_coupled_20.0s/constellation.csv): seeded generator parameters.
- [engine_summary.csv](III_output/walker_interlinks/seed_20261003/force_coupled_20.0s/engine_summary.csv): runtime, solver settings, counts and physical/diagnostic results.
- [trajectory.csv](III_output/walker_interlinks/seed_20261003/force_coupled_20.0s/trajectory.csv): all 4,000 spacecraft at eight times, with position, velocity, mass, SMA and laser integrals.
- [accepted_step_schedule.csv](III_output/walker_interlinks/seed_20261003/force_coupled_20.0s/accepted_step_schedule.csv): selected endpoint pairs at each scheduling time.
- [accepted_step_counts.csv](III_output/walker_interlinks/seed_20261003/force_coupled_20.0s/accepted_step_counts.csv): active-link counts, including initialization.

Duplicate solver timestamps from before/after callbacks are exported once, using
the last stored state. The eight distinct times are approximately 0, 0.012742,
0.151544, 1.313630, 6.090144, 11.090144, 16.090144 and 20 seconds.

The saved files passed **49 assertions**: correct seed/counts/endpoint, finite
states, complete spacecraft coverage, accepted-step history, positive laser
effect and SMA/diagnostic consistency. At every saved time, availability and
scores were recomputed from the exported positions/velocities; the saved link
was eligible, positive-scoring and the maximum-score candidate. The in-run checks
also verify terminal exclusivity and strictly increasing scheduling times.

The small force-coupled regression uses 40 spacecraft for two seconds and adds
10 tests. Agreement of physical SMA change with its laser integral confirms
force application in this two-body scenario, but is not a timestep-convergence
study. This remains a sparse single-target laser network: all 4,000 trajectories
are propagated, while only one target-helper pair fires at a time. It does not
establish all-pairs or long-duration integration performance.

## Follow-up: bounded legacy/v5/v6 Walker comparison

The same 4,000-satellite, 20-second force-coupled case was compared on 2026-10-05
using [II_examples/compare_walker_integration.jl](II_examples/compare_walker_integration.jl).
This measures propagation and accepted-step scheduling together, not just the
selected tree-matching functions. Only this comparison example was added;
production code and the legacy/v5 sources were not changed.

Each version ran sequentially in a fresh, single-threaded Julia process with a
**120-second total process limit**, including loading, compilation, exports and
repeats (TERM at the limit, KILL after a further five-second grace period).
**All three finished without a timeout**, each completing one first run and
five warmed runs. All 18 integrations reached 20 seconds with 7 accepted steps,
zero rejected steps, and one selected link at all eight scheduling times.

The runner reconstructs the constellation from its saved Walker metadata and
checks every initial Cartesian state against the previous v6 trajectory. All
versions use two-body gravity, Tsit5, 1e-9/1e-11 tolerances, a 5-second maximum
step, the same laser settings and the same target/helper candidates. A shared
SavingCallback records all 4,000 SMA values at 0, 1, ..., 20 seconds without
imposing additional integration stops.

| Version | First measured call (s) | Warm median (s) | Warm min-max (s) | Final target SMA change (m) |
|---|---:|---:|---:|---:|
| Legacy | 18.928473 | 0.221136 | 0.182392-0.741497 | 0.535846548 |
| v5 | 32.666690 | 2.391378 | 2.346296-2.489050 | 0.535945603 |
| v6 | 29.868080 | 0.097402 | 0.096927-0.097880 | 0.535945603 |

For this workload, v6's warmed median was **24.55x faster than v5** and **2.27x
faster than legacy**. Median Julia allocations were approximately 2,153.6 MB for
legacy, 54.33 MB for v5 and 53.42 MB for v6 (decimal MB). Legacy had a GC-related
timing outlier; the full range is shown rather than treating one run as definitive.
Allocations are cumulative, not peak memory, and exclude some Python/native work.

First-call times include compilation inside the measured simulation call; warmed
times do not include that first-use cost. Package loading, outer configuration
construction, assertions and CSV export are outside the measured call. Common
SMA sampling and each version's required callbacks are inside it. This additional
sampling means the first-call timings are not the identical harness used for the
earlier 27.532-second v6-only result. No general all-pairs scaling claim follows
from this sparse single-target case.

### Numerical comparison

- **v5 and v6 are exactly equal** in all 84,000 saved SMA values (21 per spacecraft),
  all 8,000 initial/final state rows, and the complete eight-row schedule.
  V6 also reproduces the original v6-only initial/final states exactly.
- Legacy uses the earlier reference workflow's `laser_impulse_callback` and
  `laser_link_scheduler_callback`. The impulse callback applies discrete velocity
  kicks; it is not the same numerical formulation as v5/v6's continuous laser RHS.
- Relative to v6, legacy's final target SMA is 0.000099055 m lower. The largest
  SMA difference at the shared one-second sample times is 0.120970 m.
- The maximum final position difference over all spacecraft is 0.000652870 m;
  the maximum final velocity difference is 0.00000611892 m/s.
- Legacy selects the same sequence of helpers, but scheduling timestamps differ
  by up to 0.394444 s. These are measured differences, not a convergence certificate.
- All 18 runs passed 13 assertions each (**234 checks**), and reading back and
  comparing the saved outputs passed **37 additional checks**.

Saved aggregate tables:
[version_comparison.csv](III_output/walker_interlinks/seed_20261003/force_coupled_20.0s/version_comparison/version_comparison.csv)
and [numerical_differences.csv](III_output/walker_interlinks/seed_20261003/force_coupled_20.0s/version_comparison/numerical_differences.csv).
The same directory contains each version's six timing rows, 84,000-row SMA table,
8,000-row initial/final state table and accepted-step schedule. Raw trajectory
tables represent the final completed warmed pass; timing rows preserve all passes.
The original v6-only output files were retained unchanged.

Fixture note: the earlier spacecraft element CSV has `_deg` angular headers but
contains the radians stored inside `InitialCondition`. This comparison does not
use those mislabeled columns; it uses the degree-valued Walker metadata and
verifies identical Cartesian states. That unrelated CSV was not changed here.

## Follow-up: 15-minute v6 force-coupled run

On 2026-10-05, the existing v6 runner completed **900 seconds (15 minutes)** of
full force-coupled integration for the same 4,000-satellite Walker constellation.
No source changes were needed. Seed 20261003, the 3,999 target-to-helper
candidates, one terminal per spacecraft, laser settings, two-body gravity,
Tsit5, 1e-9/1e-11 tolerances and the 5-second maximum step were retained.
Scheduling still runs at initialization and after every accepted integration
step; laser acceleration is evaluated during the internal solver stages.

| Measurement | Observed |
|---|---:|
| Integrated time / return code | 900 s / Success |
| Accepted / rejected steps | 197 / 13 |
| Scheduling calls, including initialization | 198 |
| Active links | 1 at all 198 scheduling times |
| Link switches / distinct helpers | 38 / 20 |
| Minimum / median / maximum accepted step | 0.012742 / 5 / 5 s |
| Physical target SMA increase | 24.112013897859 m |
| Integrated target laser SMA change | 24.112013875254 m |
| Norm of integrated target laser delta-v | 0.012696483716 m/s |
| Maximum SMA/diagnostic residual over all saved spacecraft states | 4.93601e-8 m |
| Measured simulation-call wall time | 30.146711 s |
| Julia allocated bytes | 4,041,800,712 |
| GC time inside measured call | 1.007044 s |
| Exported trajectory rows | 792,000 |
| Trajectory CSV size | 160,998,365 bytes (about 161 MB) |

The 13 rejected steps are adaptive-solver retries, not extra scheduling calls;
197 accepted steps plus initialization give 198 calls. The solver's largest
recorded step is 5.000000000000001 s because of floating-point rounding.

This was one fresh-process, single-threaded run. The **30.15-second measurement
includes first-use compilation** inside the simulation call and excludes package
loading, outer configuration/initial-state construction, assertions and CSV
writing. The entire process, including exports, completed within its 120-second
wall-clock limit. It was not timed as a warmed benchmark. Allocated bytes are
cumulative Julia allocations, not peak memory.

Reading back the saved files passed **812 assertions**. These verify the seed,
4,000 spacecraft at every saved time, exact agreement with the earlier initial
state, finite values, final time, step limits, schedule counts and SMA/diagnostic
consistency. At all 198 scheduling times, availability and scores were recomputed
from the exported states; each selected link was eligible, positive-scoring and
the highest-score candidate. In-run checks also verified terminal exclusivity.
This is still a sparse single-target scenario, not an all-pairs network or a
timestep-convergence study. The existing regression suite was not rerun because
no Julia source files changed for this experiment.

The previous 20-second and snapshot outputs remain separate and unchanged.
New output files:

- [constellation.csv](III_output/walker_interlinks/seed_20261003/force_coupled_900.0s/constellation.csv): seeded constellation metadata.
- [engine_summary.csv](III_output/walker_interlinks/seed_20261003/force_coupled_900.0s/engine_summary.csv): solver settings, timings and final diagnostics.
- [trajectory.csv](III_output/walker_interlinks/seed_20261003/force_coupled_900.0s/trajectory.csv): all 4,000 spacecraft at 198 distinct times.
- [accepted_step_schedule.csv](III_output/walker_interlinks/seed_20261003/force_coupled_900.0s/accepted_step_schedule.csv): selected endpoint pairs.
- [accepted_step_counts.csv](III_output/walker_interlinks/seed_20261003/force_coupled_900.0s/accepted_step_counts.csv): active-link counts at each scheduling time.

## Follow-up: five-orbit v6 force-coupled run

On 2026-10-05, v6 completed the same 4,000-satellite force-coupled scenario for
**five initial orbital periods**. At the initial 550 km altitude, the two-body
period is 5,738.992358075132 s, so the requested duration was
28,694.96179037566 s, or about 7.970823 simulated hours. This fixes duration using
the initial semimajor axis; it does not count actual revolutions after laser
forces perturb each orbit.

No source changes were needed. Seed 20261003, 80 planes, 50 spacecraft per plane,
3,999 target-to-helper candidates, one terminal per spacecraft, laser parameters,
two-body gravity, Tsit5, 1e-9/1e-11 tolerances and the 5-second maximum step were
unchanged. Scheduling remained at initialization and each accepted step, with
continuous laser acceleration during solver stages. Full accepted-step trajectory
output was retained after checking available memory and disk capacity.

| Measurement | Observed |
|---|---:|
| Integrated duration / return code | 28,694.96179037566 s / Success |
| Initial orbital periods | 5 |
| Accepted / rejected steps | 5,916 / 70 |
| Scheduling calls, including initialization | 5,917 |
| Active links | 1 at all 5,917 scheduling times |
| Link switches / distinct helpers | 1,207 / 74 |
| Minimum / median / maximum accepted step | 0.012742 / 5 / 5 s |
| Physical target SMA increase | 768.858372814022 m |
| Integrated target laser SMA change | 768.858372125488 m |
| Norm of integrated target laser delta-v | 0.000024602324317 m/s |
| Maximum SMA/diagnostic residual over all saved spacecraft states | 7.66478e-7 m |
| Measured simulation-call wall time | 90.086077 s |
| Julia allocated bytes | 19,030,627,528 |
| GC time inside measured call | 4.279273 s |
| Exported trajectory rows | 23,668,000 |
| Trajectory CSV size | 4,862,393,035 bytes (about 4.86 GB) |

The delta-v entry is the norm of a vector integral in the inertial frame, not the
integral of acceleration magnitude. Its components can cancel over several
orbits even while laser work produces a positive SMA change.

The **90.09-second runtime includes first-use compilation** inside one fresh,
single-threaded simulation call. It excludes package loading, outer initial-state
construction, verification and CSV writing. The full process, including the
large export, finished within its 600-second wall-clock safeguard. No warmed
timing was measured. Allocated bytes are cumulative, not peak memory usage.

Validation read the trajectory in 64 chunks to bound memory use and passed
**23,755 assertions**. Checks covered all saved values being finite, exact initial
state agreement with the earlier 20-second run, spacecraft ordering and mass,
complete 4,000-spacecraft coverage at all 5,917 times, final time, step limits,
and physical SMA versus integrated laser diagnostics. Every saved selected link
was rechecked for eligibility, positive score and optimal score using the
exported positions and velocities. All checks passed. This remains a sparse
single-target experiment, not an all-pairs or timestep-convergence result; no
additional regression-suite run was needed because Julia source was unchanged.

Saved files, separate from all previous experiments:

- [constellation.csv](III_output/walker_interlinks/seed_20261003/force_coupled_28694.96179037566s/constellation.csv): constellation metadata.
- [engine_summary.csv](III_output/walker_interlinks/seed_20261003/force_coupled_28694.96179037566s/engine_summary.csv): solver settings, timing and final diagnostics.
- [trajectory.csv](III_output/walker_interlinks/seed_20261003/force_coupled_28694.96179037566s/trajectory.csv): all spacecraft at each accepted-step time and initialization.
- [accepted_step_schedule.csv](III_output/walker_interlinks/seed_20261003/force_coupled_28694.96179037566s/accepted_step_schedule.csv): selected endpoint pairs.
- [accepted_step_counts.csv](III_output/walker_interlinks/seed_20261003/force_coupled_28694.96179037566s/accepted_step_counts.csv): active-link counts.
- [validation_summary.csv](III_output/walker_interlinks/seed_20261003/force_coupled_28694.96179037566s/validation_summary.csv): period, coverage, switches, file size and validation statistics.

## Follow-up: all-pairs first-index policy at 1,000 km

**The requested 20-orbit integration was attempted but did not complete.** This
follow-up uses a different topology from every target-centered Walker run above:
any two distinct spacecraft can be linked. It is not another one-active-link run.

### Requested configuration

- Exactly 4,000 spacecraft, Walker 4000/80/25, seed 20261003, 53-degree inclination.
- Initial circular altitude 1,000 km, one terminal per spacecraft and 227 kg each.
- All 7,998,000 unordered spacecraft pairs registered; no candidate-pair pruning.
- `:gve_sma_first_index`, dispatched to `score_gve_sma_first_index!`: each edge is
  scored for its lower spacecraft vector index, independently of `target_idx`.
- Existing 2,000 km range and 10,000 W / B=100 laser parameters retained.
- Continuous force-coupled two-body propagation, Tsit5, 1e-9/1e-11 tolerances,
  5-second maximum steps, scheduling at initialization and every accepted step.
- Requested 20 initial periods: 126,142.37875677619 s, about 35.03955 simulated
  hours, based on an initial period of 6,307.11893783881 s.

### Minimal example changes

[II_examples/walker_interlinks.jl](II_examples/walker_interlinks.jl) now accepts
`all_pairs` and `policy_name` generator options, preserving the previous defaults.
The metadata records both. Selection validation permits multiple simultaneous
links, still enforces terminal exclusivity, and uses a set for active-link lookup
instead of scanning the selected list for each of millions of candidates.
The star-only maximum-edge check remains restricted to star experiments.

The engine runner accepts an optional wall-time limit and records requested and
actual duration separately, plus `completed` and `termination_reason`. A new
`all_pairs [orbits] [wall_limit_seconds]` CLI mode selects the requested 1,000 km
first-index configuration. Partial results have a distinct output directory and
cannot be mistaken for a completed 20-orbit trajectory. The guard checks wall
time at accepted steps; it cannot interrupt an in-progress exact matching call.
Production scoring, matching, force and solver implementations were not changed.

### Full-size initial scheduling measurement

| Stage or count | Observed |
|---|---:|
| Registered pairs | 7,998,000 |
| Registration time | 6.035208 s |
| Initially available pairs | 180,660 |
| Availability update | 0.089847 s |
| Positive-score candidate edges | 85,706 |
| First-index scoring | 0.143894 s |
| Exact matching | 114.915893 s |
| Selected simultaneous links | 2,000 |

These are first-use stage timings, not warmed medians. The full graph has cyclic
positive-score components, so v6 invokes exact NetworkX Blossom instead of using
only its tree fast path. All 4,000 terminals were matched, including 1,999 links
that do not involve spacecraft 1. Scoring the other spacecraft was not disabled.

Measurements and the initial schedule were saved before attempting integration:
[preflight_stages.csv](III_output/walker_interlinks/seed_20261003/all_pairs_1000km_gve_sma_first_index_20orbits/preflight_stages.csv),
[preflight_status.csv](III_output/walker_interlinks/seed_20261003/all_pairs_1000km_gve_sma_first_index_20orbits/preflight_status.csv),
[initial_schedule.csv](III_output/walker_interlinks/seed_20261003/all_pairs_1000km_gve_sma_first_index_20orbits/initial_schedule.csv),
and [constellation.csv](III_output/walker_interlinks/seed_20261003/all_pairs_1000km_gve_sma_first_index_20orbits/constellation.csv).

### Bounded integration result: incomplete

| Measurement | Observed |
|---|---:|
| Requested initial orbital periods | 20 |
| Actual simulated duration | 0.006812013264 s |
| Completed / return code | false / Terminated |
| Stop reason | wall_time_limit |
| Soft wall-time budget | 120 s |
| Measured simulation-call time | 266.713468 s |
| Accepted / rejected steps | 1 / 0 |
| Scheduling calls, including initialization | 2 |
| Active links at each scheduling time | 2,000 |
| Physical spacecraft-1 SMA change | 0.000199378468 m |
| Maximum physical/diagnostic SMA residual | 8.29436e-9 m |
| Julia allocated bytes | 7,708,847,800 |
| Exported trajectory rows | 8,000 |

The 120-second guard is cooperative, so the initialization and accepted-step
matching work can cause it to overshoot. An external 480-second hard process
limit was also applied; it was not reached. The measured call includes first-use
compilation and excludes outer constellation construction and CSV export.

At a 5-second maximum step, the requested duration needs at least **25,229
accepted steps and 25,230 scheduling calls**. If each call cost the measured
114.9 seconds, matching alone would take approximately **33.6 days**. This is
only a constant-cost extrapolation, not a benchmark prediction: graph structure,
scores, compilation and matching cost change with time. It is sufficient evidence
that the current all-pairs exact-matching path cannot be treated like the earlier
fast star runs. No approximate matching, fewer spacecraft, pair pruning, larger
steps or less-frequent scheduling were silently substituted.

The full-size initial schedule and partial integration are valid outputs, but
**there is no completed 20-orbit result**. Completing that experiment requires a
substantial improvement to cyclic matching performance or an explicitly approved
change to scheduling requirements. Longer storage should also be revisited before
retaining every full-state output for all 20 orbits.

Saved partial results:

- [engine_summary.csv](III_output/walker_interlinks/seed_20261003/force_coupled_126142.37875677619s_1000.0km_all_pairs_gve_sma_first_index/engine_summary.csv): explicitly incomplete status, requested/actual durations and settings.
- [constellation.csv](III_output/walker_interlinks/seed_20261003/force_coupled_126142.37875677619s_1000.0km_all_pairs_gve_sma_first_index/constellation.csv): altitude, all-pairs flag and scoring policy.
- [trajectory.csv](III_output/walker_interlinks/seed_20261003/force_coupled_126142.37875677619s_1000.0km_all_pairs_gve_sma_first_index/trajectory.csv): initial and first accepted-step states for all 4,000 spacecraft.
- [accepted_step_schedule.csv](III_output/walker_interlinks/seed_20261003/force_coupled_126142.37875677619s_1000.0km_all_pairs_gve_sma_first_index/accepted_step_schedule.csv): 4,000 selected-link rows across the two scheduling times.
- [accepted_step_counts.csv](III_output/walker_interlinks/seed_20261003/force_coupled_126142.37875677619s_1000.0km_all_pairs_gve_sma_first_index/accepted_step_counts.csv): 2,000 active links per scheduling call.
- [validation_summary.csv](III_output/walker_interlinks/seed_20261003/force_coupled_126142.37875677619s_1000.0km_all_pairs_gve_sma_first_index/validation_summary.csv): saved-output validation metrics.

The saved outputs passed **41 checks** covering correct configuration and partial
status, finite physical states, initial altitude, complete spacecraft coverage,
unique terminal use, range and positive independently recomputed lower-index SMA
rates. The selected rates were all at least 0.0274572 m/s and distances at most
1,999,966.63 m. Validation did not repeat the expensive global exact optimization.
The small 40-spacecraft fixture additionally confirms that changing `target_idx`
does not affect first-index scores, and tests both normal completion and early
termination. It uses a larger test-only range to exercise a dense small graph.

## Regression verification

**335/335 assertions passed** in the complete v6 suite after the all-pairs
follow-up (312 passed before that addition), including:

- Existing policy, continuous-force, accepted/rejected-step and one-hour small-system comparisons.
- 27 indexed-storage tests: shared states, signed transitions, range boundaries,
  activity/battery/forbidden pairs, score mirrors, force equivalence and late registration.
- 85 exact-matching tests and 23 cached-hybrid tests: brute-force small graphs,
  near-tied/subnormal weights, odd cycles, mixed components, topology invalidation,
  deterministic ties, multiple terminals and 4,000-node star/chain/disjoint graphs.
- 20 Walker tests: seeded reproducibility, circular geometry, one-period analytic
  propagation, snapshot times, schedule shape and invalid inputs.
- 10 force-coupled Walker tests: full-engine completion, scheduling history,
  actual SMA change versus laser diagnostics and invalid run durations/step limits.
- 23 all-pairs Walker tests: complete pair registration, per-pair first-index
  scoring, multiple non-target links, altitude, full integration and wall-limited
  partial-result handling.

The previous small-system 10 s versus 5 s comparison was unchanged:
maximum delta-SMA difference 0.146279331 m, delta-v difference 0.000138100594 m/s,
position difference 0.713948850 m, and two link switches in each run.

## Reproduce

Run from the outer workspace root. The benchmark label names its CSV; it does
not select an older implementation. `interlinks_baseline.csv` is the preserved
pre-edit measurement and should not be overwritten with an optimized run.

```sh
SPACEAGORA_RHS_CALIBRATE=off julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/test/runtests.jl
SPACEAGORA_RHS_CALIBRATE=off julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/test/benchmark_interlinks.jl optimized
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off julia --project=2_SpaceAGORA.jl --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/compare_pdf_case1.jl legacy 10 10 benchmark
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off julia --project=8_SpaceAGORA.jl-main_v5 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/compare_pdf_case1.jl v5 10 10 benchmark
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/compare_pdf_case1.jl v6 10 10 benchmark
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/walker_interlinks.jl
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/walker_interlinks.jl engine 20
```

For the 15-minute v6 integration with the same two-minute wall-clock safeguard:

```sh
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off timeout --signal=TERM --kill-after=5s 120s julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/walker_interlinks.jl engine 900
```

For five initial orbits, retaining full accepted-step output and allowing up to
ten wall-clock minutes for integration plus export:

```sh
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off timeout --signal=TERM --kill-after=5s 600s julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/walker_interlinks.jl engine 28694.961790375659
```

The five-orbit validation summary was generated after reading back the saved
files; the runner command alone does not regenerate that validation table.

For the bounded all-pairs, first-index 20-orbit attempt at 1,000 km:

```sh
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off timeout --signal=TERM --kill-after=5s 480s julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/walker_interlinks.jl all_pairs 20 120
```

The last argument is the cooperative simulation wall-time budget in seconds.
Always inspect `completed`, `simulated_duration_s` and `termination_reason` in
the run summary: a normal process exit can represent a saved partial trajectory.
The separate preflight and validation summary tables came from the measured
checks described above, not from this runner command.

For the bounded three-version Walker comparison, the preceding snapshot and
v6 engine commands supply the saved metadata and reference trajectory:

```sh
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off timeout --signal=TERM --kill-after=5s 120s julia --project=8_SpaceAGORA.jl-main_v6 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/compare_walker_integration.jl v6
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off timeout --signal=TERM --kill-after=5s 120s julia --project=8_SpaceAGORA.jl-main_v5 --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/compare_walker_integration.jl v5
JULIA_NUM_THREADS=1 SPACEAGORA_RHS_CALIBRATE=off timeout --signal=TERM --kill-after=5s 120s julia --project=2_SpaceAGORA.jl --startup-file=no 8_SpaceAGORA.jl-main_v6/II_examples/compare_walker_integration.jl legacy
```

Exit status 124 indicates the process hit its limit; 137 can indicate the forced
kill after the grace period. Completed passes write their timing rows before the
next pass starts. No such timeout occurred in the measurements above.

The two aggregate Case 1 CSVs and the phase/Walker summary CSVs were generated
from the saved raw results after validation, as were the three-version Walker
aggregate tables. Rerunning the commands refreshes raw files, not those aggregate
tables. This report records this measurement set, not a hardware-independent
performance guarantee.