# Precompile workload for the paper benchmark harness

Every point of the paper harness is a fresh Julia process
(`ppc_worker_cmd` in `../../parallelization_performance/execution.jl`), and so is
every process-pool worker a Monte Carlo point starts. Each one loads the packages
and then compiles the solver and right-hand-side specializations of the case it
is about to run before it can time anything. This directory builds a package
image that already holds those specializations, and the harness loads it into
every worker when asked (`SPACEAGORA_PPB_WORKLOAD=auto` or `1`).
It is off by default: with the image loaded the timed repeats measured 1 to 4
percent slower on four of five points (table below), and the default for
benchmark timing is the configuration the archived runs used.

## What it covers

`SpaceAGORAPaperWorkload/src/points.jl` derives the point list from the harness's
own definitions: every case of every P phase (`PAPER_BENCHMARK_PHASES` in
`../cli.jl`, P1 through P6p) under every mode the phase runs it with, plus the
phase's parity points. Each case is built by `ppc_single_config` at the size its
name says (ComponentArrays puts the satellite count into the state vector's type,
so each size is its own specialization), under the mode's route environment
(`ppc_mode_env_pairs`), through the same worker entry points a harness worker
calls. Campaign cases run two samples, which reaches the Monte Carlo runner
(`run_monte_carlo` with `threads=:auto`, including the policy_v2 and predictive
planners) and the pool worker's sample task. A case added to a P phase is picked
up without an edit; `test/gates/ci_ppb_workload_coverage_gate.jl` fails if a
P-phase (case, mode) is ever missing.

The workload runs each point on the harness's `test` profile, whose mission is
10 s for every P-phase case (ASSUMED: long enough to take several steps at every
case's step cap, short enough that the 4096-spacecraft cases stay cheap; only the
types matter). At the time of writing that is 72 points.

What it leaves out, and why:

- **Native-GRAM cases** (P6 trace 5, `aero_4096sat_l50_gram_lookahead_100s`, and
  P6p, `aero_4096sat_l50_gram_process_100s`; in general every case matched by
  `PPC_GRAM_LIVE_CASES` in `cases.jl`). Warming them would need GRAMSuite as a
  dependency of the workload package, which makes the image unloadable on a
  machine without the GRAM build, and it would run native GRAM inside the
  precompile process: GRAMSuite includes its native wrapper module on first
  construction, an evaluation into a closed module during incremental
  compilation, and any model reachable from the workload module's globals would
  be serialized into the image with its native handles. This is the same class
  of hazard `src/precompile_workload.jl` describes for SPICE furnish state. The
  harness launches these points without the workload.
- **The harness's own functions and closures.** A worker runs the harness code in
  its `Main`; the image holds the specializations of SpaceAGORA and its
  dependencies that the harness calls into (the solve on each configuration
  type, the ODE integrator on each state type, the campaign planner, CSV writing
  of the row type), not `Main`'s functions. Compiling those still costs a few
  seconds per process (measured below).
- **Threaded code paths as such.** Julia precompiles single-threaded, so the
  bodies of the tasks a threaded route spawns are compiled at run time; the
  per-satellite kernels they call are in the image.

What the workload keeps out of the image: it runs with its working directory and
every runtime state path (`SPACEAGORA_OUTER_ROUTE_STATE_PATH`,
`SPACEAGORA_PARALLEL_POLICY_STATE_PATH`, `SPACEAGORA_COST_CONSTANTS_PATH`,
`SPACEAGORA_RHS_CALIBRATION_PATH`, `SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH`) in a
scratch directory, so building never reads or writes this machine's learned
routing state; it caps every mode's process pool at one worker so nothing is
spawned during precompilation, and fails the point if a worker appears; and it
resets its own copies of the harness caches and SpaceAGORA's furnish set before
the image is written.

## Building

```bash
bash benchmarks/studies/paper_parallelization_benchmarks/workload/build_workload.sh
```

Three stages (`setup`, `precompile`, `check`; `all` is the default):

1. `setup_workload_env.jl` creates `output/paper_workload/env`: the repository's
   exact `Manifest.toml` plus path entries for SpaceAGORA and this package, added
   offline with `PRESERVE_ALL`; it fails if any package version moved or anything
   beyond the workload package was added. The repository's own `Project.toml`
   and `Manifest.toml` are never written.
2. Precompiles `SpaceAGORAPaperWorkload` (which runs the workload) in a
   `systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0` scope when
   systemd provides one, under `/usr/bin/time -v`. It rebuilds SpaceAGORA's own
   image first if that is stale, which is the image the harness loads anyway.
3. Asks, the same way the harness does, whether the image is current.

Environment: `SPACEAGORA_PPB_WORKLOAD_ENV` (environment directory),
`SPACEAGORA_PPB_WORKLOAD_MEMORY_MAX` (scope cap), `SPACEAGORA_PPB_WORKLOAD_SCOPE=0`
(no scope), `SPACEAGORA_PPB_WORKLOAD_LOG` (log directory).

The image is invalidated by anything that invalidates SpaceAGORA's image and by
any edit to the harness files it includes (`cli.jl`, `modes.jl`, `cases.jl`,
`trajectory_parity.jl`, `reporting.jl`, `execution.jl`, the paper `cli.jl`), so
rebuild after changing any of them. A stale image is never used; see below.

## How the harness uses it

The environment is stacked behind the repository project rather than replacing
it. A worker keeps `--project=<repo>`, so it loads the same SpaceAGORA image,
resolves every other package the same way and has the same active project
(which the runtime's process pool and state paths key on) as a worker launched
without the workload; the environment only adds the workload package, second on
`JULIA_LOAD_PATH`.

- **Controller** (`ppc_run_controller`): once per run, `ppc_workload_env()`
  starts a probe process exactly the way a worker is started and asks
  `Base.isprecompiled` about the workload package. That fails whenever SpaceAGORA
  or any other dependency changed since the image was built. When the image is
  current, `ppc_worker_cmd` adds `JULIA_LOAD_PATH="@:<env>:@v#.#:@stdlib"`,
  `SPACEAGORA_PPC_WORKLOAD=1` and `SPACEAGORA_PROCESS_WORKER_PRELOAD` to every
  worker it launches, except for the native-GRAM cases.
- **Per-point worker**: `parallelization_performance/cli.jl` loads the package
  when `SPACEAGORA_PPC_WORKLOAD=1`, after checking again that the image is
  current (a stale image would otherwise be rebuilt inside the point).
- **Harness process pool** (`ppc_ensure_process_workers!`): pool workers inherit
  the worker's environment, and load the package through the same `cli.jl`
  include. They load GRAMSuite only when the coordinator has it, i.e. only for a
  native-GRAM case: the GRAM extension's methods invalidate SpaceAGORA's
  compiled solve path. Measured on this workstation with the image loaded, one
  `--threads=1` worker running `independent_1sat_1hr`: first sample 19.9 s with
  GRAMSuite loaded, 1.9 s without.
- **Runtime process pool** (`SpaceAGORA.ParallelProcess.ensure_process_workers!`):
  new workers load every package named in `SPACEAGORA_PROCESS_WORKER_PRELOAD`
  right after SpaceAGORA, when it is precompiled for their environment. Their
  load path already carries the workload environment (`_process_worker_load_path`
  copies the coordinator's).

Turning it on: `SPACEAGORA_PPB_WORKLOAD=auto` uses the image when it is current
and otherwise runs without it; `SPACEAGORA_PPB_WORKLOAD=1` (or `on`) makes a
missing or stale image an error. The default, `0`, never loads it.
The controller prints `precompile_workload=<env>` or `precompile_workload=off`
at the start of every run.

Remote jobs: `scripts/remote/spaceagora-remote push` builds the image once per
job, after `Pkg.instantiate()` and before the command, when asked
(`--workload auto|on|off`, default `off`; `auto` builds it when the command runs the paper harness; the log is `workload_build.log` in the
job directory and the exit code `WORKLOAD_BUILD_EXIT` in `job.meta`). The
package name starts with `SpaceAGORA`, so the runner's per-job purge of
`shared_depot/compiled/v*/SpaceAGORA*` removes it together with SpaceAGORA's
image, the environment lives in the release's own `output/`, and the harness's
freshness probe refuses any image built against a different SpaceAGORA. A failed
build costs startup time, never correctness.

## Validating

```bash
julia --project=. --startup-file=no \
  benchmarks/studies/paper_parallelization_benchmarks/workload/validate_workload.jl <out_dir> --pairs=3
```

launches five representative harness workers (a P1 single-spacecraft point, a
P2 4096-spacecraft point, P3 on the process route and on policy_v2, and a P5
grid point) through `ppc_worker_cmd`, with and without the workload, and writes
`workload_validation.csv`. It uses two opt-in hooks in the harness:
`SPACEAGORA_PPC_STARTUP_TRACE=1` prints a time stamp when a worker has loaded,
when its first solve returns and when its first timed repeat starts, and
`SPACEAGORA_PPC_DUMP_STATE_DIR=<dir>` writes every sample's step times and final
state as raw Float64 for a byte comparison. Run it on a quiet machine, inside a
memory-capped scope, with nothing else timing.

## Measured effect

Every number in this section was measured on one workstation (AMD Ryzen 9
9900X, 12 cores / 24 threads, 60 GB, Julia 1.12.1) unless it says otherwise.
Nothing here was measured on the benchmark box.

### Why it matters

DERIVED from the archived cold run `trx50_ppb_cold_20260922_002429` (TRX50,
P1-P5, 11 repeats, 3 warm-ups, 199 points), from its raw CSV and manifest:
the run spanned 9.11 h, of which the timed repeats sum to 3.88 h and three
warm-ups per point (3 x the point's median timed wall) to 1.05 h. Taking, for
consecutive points, the gap between their `timestamp_utc` stamps minus the
earlier point's timed repeats and warm-ups leaves a median of 53 s per point
(4.18 h summed over 198 gaps), and a median of 128-134 s for `outer_process`,
`policy_v2` and `predictive` points, whose pool workers are fresh processes too;
37-47 s for the in-process routes. That remainder is process startup: package
loading plus compilation.

### Startup (this workstation)

`validate_workload.jl --pairs=3`: five representative points, each launched
three times per variant, alternating which variant goes first; the values are
medians over the three launches (MEASURED; per-launch values in
`results/space-falcon-1_20260924/workload_validation.csv`). Times are seconds
from process start. "First timed" is when the first timed repeat starts, after
the warm-up solve and, for a campaign, the pool's provisioning and warm-up; it
includes one warm-up solve's compute in both variants.

| Point (phase, case, mode, threads x workers) | Import stock / workload | First timed stock / workload | Ratio | Worker wall stock / workload | Ratio | Saved |
|---|---|---|---|---|---|---|
| P1 `gravity_1sat_l50_vacuum_4150000s`, serial, 8 x 1 | 5.5 / 7.3 | 49.7 / 20.9 | 2.4x | 76.5 / 47.1 | 1.6x | 29 s |
| P2 `gravity_4096sat_l50_vacuum_5800s`, inner_only, 8 x 1 | 5.5 / 7.3 | 39.5 / 15.5 | 2.6x | 49.2 / 25.3 | 1.9x | 24 s |
| P3 `independent_1sat_1hr`, outer_process, 256 samples, 4 x 4 | 5.6 / 7.3 | 135.3 / 27.0 | 5.0x | 151.8 / 42.9 | 3.5x | 109 s |
| P3 `independent_1sat_1hr`, policy_v2, 256 samples, 4 x 4 | 5.6 / 7.3 | 139.1 / 30.6 | 4.6x | 150.1 / 41.0 | 3.7x | 109 s |
| P5 `mcgrid_16sat_8mc`, outer_inner_static, 4 x 2 | 5.7 / 7.2 | 62.5 / 13.7 | 4.6x | 78.3 / 28.8 | 2.7x | 50 s |

Loading the image adds about 1.7 s to every process (import column). Most of
what remains between import and the first timed repeat is compute (the P1
warm-up solve alone is ~8.4 s) and, for the harness's own functions, compilation
of `Main` code: 4.5 s of the P1 worker's 7.4 s of remaining compile time was
`main_parallelization_performance` itself (MEASURED with `--trace-compile-timing`).
On the two P3 points, a pool worker's time from loaded to its first sample
returned fell from 60.0-61.5 s to 1.9-3.1 s (MEASURED, 4 workers x 3 launches
each); what is left of those points' startup is the coordinator's own start and
the pool workers' process start and package loading (11-12 s, DERIVED from the
same traces: the last pool worker's first sample minus the coordinator's first
solve, minus the 2-3 s sample).

A controller run of one P5-style point plus seven parity points
(`parallelization_performance.jl full --cases=mcgrid_16sat_8mc --modes=outer_process
--threads=4 --process-workers=2 --repeats=2`) took 4 m 18 s with the workload and
8 m 26 s with `SPACEAGORA_PPB_WORKLOAD=0` (MEASURED, one run each).

### Results are byte-identical

Every launch above ran with `SPACEAGORA_PPC_DUMP_STATE_DIR`, which writes each
sample's step times and final state as raw Float64. In all 15 stock/workload
pairs every file compared equal byte for byte: 3 files per launch for P1 and P2,
40 for P5, 1284 for each P3 campaign (every warm-up and timed sample on the
coordinator and on the pool workers).

### Timed repeats run slightly slower

The timed medians are not unchanged. Median timed wall, stock / workload,
medians over the three launches, and the per-launch shift (MEASURED):

| Point | Timed median stock / workload | Shift per launch pair |
|---|---|---|
| P1 1 sat, serial | 8.36 / 8.42 s | -0.2%, +2.2%, +1.3% |
| P2 4096 sat, inner_only | 2.74 / 2.87 s | +2.0%, -9.5%, +14.1% (launch-to-launch spread dominates) |
| P3 outer_process | 2.46 / 2.54 s | +4.3%, +3.7%, +2.3% |
| P3 policy_v2 | 1.63 / 1.69 s | +2.7%, +4.3%, +3.2% |
| P5 outer_inner_static | 2.69 / 2.74 s | +1.2%, +2.0%, +2.0% |

Allocations are identical and GC time is the same within noise; the extra time
is in the samples' own compute (the rows' `sample_wall_time_sum_s`). The likely
cause, NOT verified here, is that code loaded from a package image reaches other
functions and constants through relocation slots where freshly compiled code
embeds them directly. Whatever the cause, it applies to every timed repeat of a
run that uses the image: do not put rows measured with the workload and rows
measured without it in the same comparison. The route-to-route comparisons
within one run all carry it, at 1-4% on the points measured here.

### Precompile cost (this workstation)

`build_workload.sh all`, three builds: 3 m 40 s wall for the precompile stage
each time, peak RSS 9.18-9.29 GB (`/usr/bin/time -v`, MEASURED), well inside the
16 GB scope; setup takes about 4 s. The resulting image is 250 MB. When
SpaceAGORA's own image is stale as well (always the case on a remote job, whose
depot is purged), the precompile stage rebuilds it first; that took 48 s and
3.8 GB peak on its own here (MEASURED, one build).

### Values that are assumed

- The workload's 10 s mission per point (the harness's `test` profile). ASSUMED
  enough to take several steps and reach every code path the case's types
  need; the result above is the evidence that it does.
- Two samples per campaign case. ASSUMED enough to reach the campaign runner and
  its planners; the sample count does not enter a type.
- The 16 GB default memory cap of the build scope, chosen as the local
  measurement rule's cap; the measured peak is 9.3 GB.
