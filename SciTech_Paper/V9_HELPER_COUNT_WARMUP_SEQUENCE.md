# v9: Changing Helper Count in One Julia Process

## Question and method

Does running 10 helpers first reduce the first-use cost of subsequent
20-helper and 30-helper scenarios?

As requested, exactly **three simulations** were run, in this order:

1. 10 helpers, once.
2. 20 helpers, once.
3. 30 helpers, once.

All three ran in **one fresh Julia process**, with no preliminary simulation
and no repeated scenario. The package and runner were loaded once. Each
simulation built fresh initial states and models. The same runner function
accepted helper count as a runtime keyword; its source was not re-included
or redefined between cases.

Physical and numerical settings match the previous
[100-orbit SMA comparison](RUNTIME_COMPARISON_WITH_V5.md):

- Helpers: circular 1000 km altitude, 0 degree inclination, evenly phased.
- Target: circular 1050 km altitude, 1 degree inclination, initial anomaly,
  RAAN, and argument of periapsis zero.
- 100 initial target periods: 637,134.036049874 simulated seconds.
- Target-only `:gve_sma` objective, target-helper links only, one terminal each.
- Mass 227 kg each; laser power 10 kW, amplification 100, range 200 km.
- Central gravity only; no J2 or drag.
- Native v9 `run_simulation`, Tsit5, relative tolerance 1e-9,
  absolute tolerance 1e-11, maximum step 10 s.

The runs were sequential with one Julia thread and one OpenBLAS thread.
A full garbage collection was requested before each call, outside its
timer. Environment: Julia 1.12.1, Linux, AMD Ryzen 7 9700X 8-Core Processor,
2026-10-18. Other activity on this shared machine was not controlled.

## Results

| Execution order | Helpers | Native simulation-call wall time (s) | Complete runner wrapper (s) |
|---|---:|---:|---:|
| First | 10 | **10.086720279** | 10.762283794 |
| Second | 20 | **5.648067218** | 5.687911015 |
| Third | 30 | **5.912026199** | 5.951026474 |

The initial runner/package load took a separate **6.452295805 s**.
It is excluded from both columns. Neither column measures full Julia
process startup-to-exit time.

The primary timer includes all work inside native `run_simulation`,
including any compilation and initialization. The outer wrapper additionally
includes configuration building, post-processing, CSV writing, and printing.

### Interpretation

- The 20-helper call took **44.00% less time** than the first 10-helper call,
  despite containing more spacecraft.
- The 30-helper call took **41.39% less time** than the first 10-helper call.
- The 30-helper call took **4.67% more time** than the 20-helper call.

**The observed sequence is consistent with useful warm-up carrying over
between helper counts.** The first 10-helper scenario does not force the
next two scenarios to pay the entire first-process cost again.

However, the later calls must not be labeled fully warmed runs. v9's
`ComponentVector` state layout changes with constellation size, so changing
helper count can introduce new concrete types and compilation. Shared
compiled code remains reusable. We did not instrument compilation time
in this experiment.

Increasing helper count also changes the computational workload, geometry,
and selected links. Thus, shared compilation, new specialization,
simulation work, and system noise are mixed in these wall times. The
30-helper call need not be faster than the 20-helper call.

There is only one measurement per count, no fresh-process 20/30-helper
control, and no same-count warm repeat, as requested. This experiment
therefore cannot quantify the isolated compilation savings, establish
warm runtime for each count, or provide a general helper-count scaling law.

**Practical conclusion:** keeping a scenario sweep in one Julia process
can reduce later first-use overhead. It does not guarantee that each
successive scenario gets faster or that every new helper count is fully warm.

## Physical results and verification

| Helpers | SMA increase (m) | Final inclination (deg) | Scheduler records |
|---|---:|---:|---:|
| 10 | 400.593135499 | 1.000027344982 | 63724 |
| 20 | 797.918490289 | 0.999994161422 | 63724 |
| 30 | 1196.500667118 | 0.999972534895 | 63723 |

All three returned `Success` and reached exactly the requested end time.
Validation checked saved-summary agreement, finite final metrics,
initial inclination, helper counts, policy, duration, schedule CSV shape,
monotonic timestamps, valid helper indices, normalized orbit times,
active-update/activation counts, and final helper.
These checks read the existing outputs and did not rerun any simulation.

The shared native runner now supports `run_case(; helper_count=...)`.
Omitting the keyword retains the original configured helper count and
output paths. No production v9 source was changed.

## Saved code and results

- [Sequence benchmark](benchmark_v9_helper_count_sequence.jl)
- [Shared native runner](native_interlink_comparison.jl)
- [Timings, execution order, and process metadata](runtime_benchmarks/v9_helper_count_sequence/sequence.toml)
- [Per-case physical summaries and schedules](runtime_benchmarks/v9_helper_count_sequence/results/)

To reproduce from the workspace root:

```sh
GKSwstype=100 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --project="9_SpaceAGORA.jl-main_full integration" \
  SciTech_Paper/benchmark_v9_helper_count_sequence.jl
```

This starts one process and executes each of the three counts once.
It overwrites only this sequence's outputs, not the earlier benchmark data.
