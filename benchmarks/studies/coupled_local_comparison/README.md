# Coupled-case comparison on one host

Measures, on whatever machine runs it, the comparisons that the October 2026
coupled-case experiment reported from a different machine: the PR #222 speedup
(controller readout specialization and initial-force recording), optimized
serial against static and adaptive parallel routes, and how the result changes
with spacecraft count. Numbers from another host are not comparable with these,
so nothing here is checked against them; every ratio is taken between runs made
on the same host by the same script.

## Workload

`stack<N>_e6_actuated_saved`, defined in
`benchmarks/studies/parallelization_performance/cases.jl`. It is the B11
`stack<N>_e6_actuated` rung with trajectory output turned on:

- degree-20 Earth gravity (EarthGGM05C), Sun and Moon point-mass gravity,
  solar radiation pressure, `AerodynamicCoefficientfM` aerodynamics with the
  analytic `ExponentialAtmosphereModel` (no native GRAM);
- 6-DOF attitude with four panels per spacecraft, an LVLH cascade attitude
  controller, and one `MagneticMomentumManagerModel` control effector per
  spacecraft at 1 Hz;
- `results=true`, so the saving callback and per-step solution storage are on.

The mission is 120 s under the `smoke` profile at every N (10 s under `test`),
so the sweep varies only N.

What is assumed rather than taken from the original experiment: the saving
cadence is the harness default `data_rate = 10 s`, and the controller gains are
the B11 routing-benchmark values, which make no flight-fidelity claim (see the
comment on `multi_8sat_magnetorquer_attitude`). Spacecraft do not interact with
each other.

The control path calls every control effector for every spacecraft, so its cost
grows as N squared. PR #222 made each call cheaper without changing that, so
expect the larger-N rungs to be dominated by it.

## Code under comparison

| Tag | Commit | What it is |
|---|---|---|
| `pre` | `ea42d2946` | main immediately before PR #222 |
| `post` | `402018b98` | the PR #222 merge; its only `src/` change from `pre` is #222 |

The harness is identical at both commits. `run.sh` copies three harness files
from this checkout into both trees, since the case, the per-call state dumps
and `--parity-cases=none` were added after them.

## Running

```bash
R=/path/to/worktrees
PRE_TREE=$R/coupled-ea42d2946 POST_TREE=$R/coupled-402018b98 \
  benchmarks/studies/coupled_local_comparison/run.sh A <outdir>   # 32 spacecraft, both codes
PRE_TREE=... POST_TREE=... \
  benchmarks/studies/coupled_local_comparison/run.sh B <outdir>   # N sweep, post only
julia --project=. benchmarks/studies/coupled_local_comparison/compare.jl <outdir>
```

Both trees need the `data/` links described in `CLAUDE.md` and an instantiated,
precompiled environment before the first timed run.

Stage A runs `serial` (R0) at one thread, then at 2, 4, 8 and 12 threads the
static routes `inner_only` (R2), `rhs_satellite`, `rhs_per_satellite`,
`rhs_flat` (R7 with the RHS plan pinned and calibration off) and the adaptive
`predictive` route (R7, the policy the paper describes). Stage B runs serial and
the same five routes at each `SWEEP_THREADS` count (default 12) for N = 32, 64,
128, 256 on the post-#222 code. N = 32 repeats Stage A's case so that a sweep
run on another host carries its own anchor; compare N only within one host. Override with
`A_CASE` (Stage A's case), `THREADS`, `ROUTES`, `SWEEP_N`, `SWEEP_THREADS`, `REPEATS` (default 6),
`WARMUP` (default 3), `PROFILE` and `MEMORY_MAX` (per-point cap, default
16G).

Each point is a separate controller invocation with its own state-dump
directory, under a 16 GB cap. A finished point is skipped when the script is
rerun. Two guards keep other work off the CPU while a point is timed. The point
holds every lock in `LOCKS` (by default the agent Julia lock and the paper
repository's heavy-job lock). Once it holds them, it measures each Julia
process's CPU use over `QUIET_SAMPLE_S` (3 s) from `/proc` and refuses to start
if any exceeds `QUIET_MAX_CORES` (0.1 core), which catches jobs that take no
lock, or if `MemAvailable` is under `MIN_AVAILABLE_GB` (24 GB), since swapping
would distort the timings. A refused point releases the locks, waits `QUIET_RETRY_S` (120 s) and
retries for up to `QUIET_WAIT_S` (8 h). The guard sees only Julia processes;
the harness's own load-average check (50% headroom) covers the rest coarsely.

## Output

`compare.jl` writes `summary.csv` and prints:

- every point: median, minimum and maximum wall time of the timed repeats,
  speedup over the same code's serial run, allocation and GC time per solve,
  and the RHS plan the route installed;
- adaptive against the fastest static route at each thread count;
- the PR #222 effect, pre over post median for the same route and thread count,
  with allocation before and after;
- bit identity: every timed solve dumps its step times and final state
  (warm-up solves do not), and the summary reports whether all repeats within a
  point, all routes
  and thread counts of one code, and the two codes produce the same bytes.
  Saved force fields are not in the dump; the initial-force correction is
  covered by `test/unit/simulation/initial_force_output_tests.jl`.
