# Heap contention and process-worker memory

The original investigation, with its measurements and limitations, is kept
at this same relative path in the lab's private research archive (see
[Record a result](../../../docs/src/user/studies_benchmarks.md#record-a-result)).
This page keeps the current methods and settings. An old result does not show
that a later revision performs the same or gives the same science.

Allocation rate, garbage collection, native-library memory and retained results
are different costs. Measure each one before blaming a thread-scaling limit or
growing worker memory on it. A result at one thread count, workload or host
says nothing about another.

## Allocation attribution

`benchmarks/studies/heap_contention/alloc_profile.jl` warms a solve up, times
a second one with `@timed`, and profiles a third with `Profile.Allocs`. It
assigns each sampled allocation to the innermost SpaceAGORA source frame; use
that bucket, not the file that started the solve, to decide who owns the
allocation. Solver work and callback storage are counted apart from writing
results.

Allocation sampling keeps stacks and can use far more memory than the
simulation alone. The tool warns when the sample rate is above 0.01 for any
case other than the small `independent_1sat_1hr`; its usage example for the
large `montecarlo_heavy_aerobraking` case uses `--sample-rate=0.005`. Sampled
bucket fractions are estimates, not exact totals, and the
`estimated_bytes_at_1x` column is labelled as an estimate. Put a hard memory
limit on these probes. The Linux wrappers supplied here use:

```bash
systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
    julia --project=. ...
```

That wrapper needs a systemd user scope. Choose a limit that suits the host
and the workload; a Julia heap-size hint does not replace it.

`benchmarks/studies/heap_contention/dump_state.jl` writes the saved times and
full states as raw `Float64` for comparison across revisions. Keep the inputs,
solver settings and seeds identical, and check which route actually ran.
Checkpoint and gravity-backbone paths need cases that use them; an ordinary
solve cannot vouch for a change to a path it never runs. The debug
initial-derivative problem is built only when
`SPACEAGORA_DEBUG_INITIAL_DERIVATIVE` turns on its debug branch.

## Opt-in state thinning

These settings affect written CSV/results output, not the `ODESolution` returned
by `run_simulation(...; return_solution=true)`:

| Setting | Default | Effect |
|---|---|---|
| `SPACEAGORA_RESULTS_THIN_STRIDE=k` | `1` | Retain every k-th saved state and always the final state. |
| `SPACEAGORA_RESULTS_FINAL_ONLY=1` | off | Retain only the final state; takes priority over the stride. |

With the defaults the time and data vectors are unchanged. Thinning picks from
the saved states without interpolating. It cuts the work of writing results;
it does not undo what the solver and callbacks already accumulated during the
solve, so it does not lower the solver's peak allocation.
`test/unit/simulation/solution_storage_tests.jl` checks which states are kept
and what is written, including that the final time is kept.

## Collector settings

`benchmarks/studies/heap_contention/collector_grid.sh` runs a series of new
Julia processes with combinations of `--gcthreads` and `--heap-size-hint`,
using `collector_grid.jl` for repeated warmed-up solves. Its grid and its Linux
memory wrapper are choices for this study, not run-time defaults.

The parallelization harness accepts `SPACEAGORA_PPC_WORKER_GCTHREADS` (unset:
no GC-thread flag) and `SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT` (unset: the
derived worker hint below). For example, on a suitable machine with a memory
limit:

```bash
SPACEAGORA_PPC_WORKER_GCTHREADS=2 SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT=4G \
  julia --project=examples benchmarks/studies/parallelization_performance.jl \
  --cases=independent_1sat_1hr,montecarlo_heavy_aerobraking \
  --modes=outer_threads --threads=8,16,32 --warmup=1
```

Record the CPU affinity, physical and logical cores, Julia version, GC
settings, other load on the machine and the full timing distribution. Compare
arms repeated often enough on the target machine. Neither a run with fewer
threads nor a noisy shared machine tells you how a larger thread count
behaves.

## Process-pool heap hints

The shared formula in `src/orchestration/parallel/process/worker_pool.jl` is:

```julia
clamp(floor(Int, share * Sys.total_memory() / max(1, pool_size + 1)),
      2 * 1024^3, 64 * 1024^3)
```

`SPACEAGORA_POOL_WORKER_HEAP_HINT_SHARE` defaults to `0.5`. The memory is
split between the workers and their coordinator, and the share leaves room for
native images, the operating system and other work. The share and the 2 GiB
floor and 64 GiB cap are chosen policy values, not bounds derived from any
workload. A hint is a soft target for the collector and can be exceeded.
Because of the floor, on a machine with little memory the workers together can
exceed the share, so check the result there.

Three launch sites use this formula:

| Process | Implementation | Override |
|---|---|---|
| General campaign pool | `ParallelProcess._process_worker_exeflags` | `SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT` |
| Harness per-point controller | `ppc_worker_cmd` / `_ppc_worker_gc_flags` | `SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT` |
| Harness Distributed pool | `_ppc_pool_worker_exeflags` | `SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT` |

At each site, `off` removes the hint, any other non-empty value is passed to
Julia as given, and unset uses the derived value. The controller splits the
memory by its configured worker budget, even for a point that runs serially.
`test/unit/parallel/pool_worker_heap_hint_tests.jl` checks the formula and all
three command builders; it does not measure real worker memory.

`pool_rss_probe.sh` compares arms with and without the hint and then calls
`summarize_pool_rss.py`. Its Linux wrapper caps memory at 16 GiB. The sampler
logs every Julia process on the host, and the summary counts only the probe's
coordinator and its direct `--worker` children. Still inspect the process tree
and the run logs, failed arms included. A hint above the observed peak never
came into play, so it cannot show that the hint limits growth. Measure the
worker and sample counts you care about instead of extrapolating from a
smaller run.

## Native memory and validation limits

A native GRAM model holds memory outside Julia's managed heap. Check epoch
rebuilding, model reuse, finalizer timing and allocator retention separately
from the GC settings.
[Pool-worker memory growth](../../../docs/architecture/gram_thread_scaling.md#pool-worker-memory-growth)
describes the model-lifecycle controls and their diagnostic tools.

Collector flags do not change the equations, but that alone is not a
trajectory comparison. Claims about performance, memory or identical bytes need
evidence from the actual revision, routing mode, inputs and machine. Keep the
raw results and run details with the private study record.
