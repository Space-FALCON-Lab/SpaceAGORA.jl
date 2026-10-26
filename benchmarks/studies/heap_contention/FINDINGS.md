# Heap contention and process-worker memory

The original investigation, including its measurements and limitations, is
preserved at this same relative path in the private archive described by the
[research-record policy](../src/user/studies_benchmarks.md#research-records-and-generated-results).
This page keeps current methods and configuration; historical results do not
qualify performance or scientific equivalence at a later revision.

Allocation rate, garbage collection, native-library memory and retained results
are different costs. Measure them separately before attributing a thread-scaling
limit or worker-memory growth to one of them. A result at one thread count,
workload or host does not qualify another regime.

## Allocation attribution

`benchmarks/studies/heap_contention/alloc_profile.jl` warms a solve, times a
second with `@timed`, and profiles a third with `Profile.Allocs`. It assigns
sampled allocations to the innermost SpaceAGORA source frame. Use the observed
bucket, rather than the file that launched the solve, to identify ownership.
Solver work and callback storage are separate from results persistence.

Allocation sampling retains stacks and can consume substantially more memory
than the simulation alone. The tool warns about high sample rates outside the
small P3 case; the large P4 case's usage example uses `--sample-rate=0.005`.
Sampled bucket fractions are estimates, not exact totals; its
`estimated_bytes_at_1x` column is labeled accordingly. Use a hard memory limit
for these probes. The supplied Linux wrappers use:

```bash
systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
    julia --project=. ...
```

That wrapper requires a systemd user scope. Choose a limit appropriate to the
host and workload; a Julia heap hint is not a replacement for it.

`benchmarks/studies/heap_contention/dump_state.jl` writes saved times and full
states as raw `Float64` for comparison across revisions. Preserve identical
inputs, solver settings and seeds and verify the route actually exercised.
Checkpoint/backbone paths need cases that use them; an ordinary solve cannot
qualify a change to an unexercised path. The debug initial-derivative problem is
constructed only when `SPACEAGORA_DEBUG_INITIAL_DERIVATIVE` enables its debug
branch.

## Opt-in state thinning

These settings affect written CSV/results output, not the `ODESolution` returned
by `run_simulation(...; return_solution=true)`:

| Setting | Default | Effect |
|---|---|---|
| `SPACEAGORA_RESULTS_THIN_STRIDE=k` | `1` | Retain every k-th saved state and always the final state. |
| `SPACEAGORA_RESULTS_FINAL_ONLY=1` | off | Retain only the final state; takes priority over the stride. |

The default returns the same time/data vectors. Thinning selects existing saved
states without interpolation. It reduces results-writing work; it does not
remove solver/callback accumulation already performed during the solve and
should not be advertised as a reduction of the solver's peak allocation.
`test/unit/simulation/solution_storage_tests.jl` checks the selection contract
and written output, including preservation of the final time.

## Collector settings

`benchmarks/studies/heap_contention/collector_grid.sh` runs a sequence of fresh
Julia processes with combinations of `--gcthreads` and `--heap-size-hint`, using
`collector_grid.jl` for repeated warmed solves. Its current grid and Linux
memory wrapper are study choices, not universal runtime defaults.

The regular parallelization harness accepts
`SPACEAGORA_PPC_WORKER_GCTHREADS` (unset: no GC-thread flag) and
`SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT` (unset: the derived worker hint below).
For example, with a suitable machine and memory limit:

```bash
SPACEAGORA_PPC_WORKER_GCTHREADS=2 SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT=4G \
  julia --project=examples benchmarks/studies/parallelization_performance.jl \
  --cases=independent_1sat_1hr,montecarlo_heavy_aerobraking \
  --modes=outer_threads --threads=8,16,32 --warmup=1
```

Record CPU affinity, physical and logical cores, Julia version, GC settings,
concurrent load and the full timing distribution. Compare sufficiently repeated
arms on the target machine. Neither a lower-thread experiment nor a noisy
shared-machine run establishes behavior at a larger thread count.

## Process-pool heap hints

The shared formula in `src/orchestration/parallel/process/worker_pool.jl` is:

```julia
clamp(floor(Int, share * Sys.total_memory() / max(1, pool_size + 1)),
      2 * 1024^3, 64 * 1024^3)
```

`SPACEAGORA_POOL_WORKER_HEAP_HINT_SHARE` defaults to `0.5`. The split accounts
for workers plus their coordinator; the share reserves headroom for native
images, the OS and other work. The share and 2/64 GiB floor/cap are assumed
policy values, not workload-derived memory bounds. A hint is a soft collector
target and can be exceeded. The floor also means small-memory deployments must
review the result instead of assuming the aggregate stays within the share.

Three launch sites use this formula:

| Process | Implementation | Override |
|---|---|---|
| General campaign pool | `ParallelProcess._process_worker_exeflags` | `SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT` |
| Harness per-point controller | `ppc_worker_cmd` / `_ppc_worker_gc_flags` | `SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT` |
| Harness Distributed pool | `_ppc_pool_worker_exeflags` | `SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT` |

At each site, `off` removes the hint, a nonempty explicit value is passed to
Julia, and unset uses the derived value. The controller uses its configured
worker budget for the split, even for a point that runs serially.
`test/unit/parallel/pool_worker_heap_hint_tests.jl` checks the formula and all
three command builders; those checks do not measure actual worker memory.

`pool_rss_probe.sh` compares hinted and unhinted arms and calls
`summarize_pool_rss.py`. Its default Linux wrapper caps memory at 16 GiB and its
sampler observes all Julia-named processes, so unrelated Julia jobs contaminate
the totals. Inspect the process tree and run logs, including failed arms. A
hint above the observed peak was not binding and cannot demonstrate that the
hint controls growth. Measure the intended worker and sample counts rather
than extrapolating from a smaller run.

## Native memory and validation limits

A native GRAM model has memory outside Julia's managed heap. Epoch rebuilding,
model reuse, finalizer timing and allocator retention must be checked separately
from GC settings. See [Pool-worker memory growth](gram_thread_scaling.md#pool-worker-memory-growth)
for the model-lifecycle controls and their diagnostic tools.

Collector flags do not change the equations, but that observation alone is not
a fresh trajectory comparison. Performance, memory and byte-identity claims
need evidence from the actual revision, routing mode, inputs and machine.
Keep the raw results and run metadata with the private study record.
