# Parallelization measurement and validation

The original investigation, including its measurements and limitations, is
preserved at this same relative path in the private archive described by the
[research-record policy](../src/user/studies_benchmarks.md#research-records-and-generated-results).
This page keeps current methods and configuration; historical results do not
qualify performance or scientific equivalence at a later revision.

This reference keeps the measurement requirements used by the parallelization
harness. Historical branch-specific fixes, timing tables and unresolved
hypotheses remain in the original investigation; they are not a current defect
list or a performance qualification for new code.

## Timing scope

- Forward `--warmup` into each worker process. A repeated point in a fresh
  subprocess still pays compilation unless that subprocess warms the work.
- For steady-state throughput, provision and warm the Distributed pool before
  timing. Report process startup or cold compilation separately when those
  costs are part of the intended user workflow.
- Record physical cores, SMT threads, CPU affinity and concurrent load. Make
  the resource ladder explicit; `SPACEAGORA_PPC_PHYSICAL_CORES` overrides the
  harness's physical-core detection.
- Record actual outer concurrency. `SPACEAGORA_OUTER_PARALLEL_ACTIVE` means an
  enclosing split owns the resource budget; a single constellation solve does
  not become outer-parallel merely because its mode has a parallel name.
- Static profiles pin their routes. Adaptive profiles must retain the
  calibration mechanism they ship with when measuring shipped behavior.
  Disabling that mechanism answers a different question and must be labeled.

The timing and mode implementation lives in
`benchmarks/studies/parallelization_performance/{cli,modes,execution}.jl`.
Record the revision and relevant environment with every result, including the
calibration/corrections store and whether the run began cold or reused it.

## Compare equivalent work

Validate the workload before interpreting a scaling curve. Check the force and
actuator stack, mission length, output enabled/disabled state, native calls,
cache path, sample count, actual worker width and inner budget. A setting with
no reader in the exercised path is not an experimental axis. In particular,
verify saved output rather than assuming `num_steps_to_save` changes the solver
or saving callback cadence.

The P6 density-case names have survived changes in their definitions; see
[GRAM benchmark configuration](gram_thread_scaling.md#benchmark-configuration-and-interpretation).
Do not combine old and new rows on the strength of a case name. Likewise,
benchmark exclusions recorded in an old investigation must be requalified
before they are lifted or described as current failures.

The paper phase definitions and
case guide (SpaceAGORAParallelizationStudies)
now live in SpaceAGORAParallelizationStudies.
The larger B7–B14 families cover constellation size, Monte Carlo throughput,
atmosphere fidelity, force/actuator stacks, interacting propagation, fixed
resource splits, and mission/output work. Preview excludes these heavy
families; use the dry run to inspect the resolved cases and ladder.

## Allocation and routing investigations

Distinguish allocation attribution from a demonstrated scaling bottleneck.
Removing allocations does not establish that allocation caused a thread
ceiling. Measure RHS work, scheduling, reduction, callback/storage and collection
costs separately. [Heap contention](heap_contention.md) describes the retained
profiling and state-dump tools.

For nested work, limiting the inner width and changing the chosen algorithm
are separate decisions. A true serial baseline must use an inner budget of one
where required, rather than silently launching parallel work from a process
that happens to have multiple threads. Record `SPACEAGORA_INNER_THREAD_BUDGET`
and `SPACEAGORA_RHS_CALIBRATE` alongside the route. The harmonics spin barrier
uses `SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=auto` by default to compare its
arms within a solve; `1` and `0` force them for controlled measurements. Record the router's actual selection
and widths instead of inferring them from a profile label.

## Repeating an evaluation

From the repository root, choose a supported case and resource budget for the
machine. This example requests warm-up explicitly:

```bash
SPACEAGORA_PPC_PROFILE=full SPACEAGORA_PPC_CASES=independent_1sat_1hr \
SPACEAGORA_PPC_MODES=serial,outer_threads SPACEAGORA_PPC_THREADS=1,4 \
  julia --project=examples benchmarks/studies/parallelization_performance.jl full --warmup=1

```

Use the tool usage headers for allocation probes and raw state dumps. The
historical investigation's temporary `steady.jl` and `allocprobe.jl` scripts
were not tracked and are not supported repository entrypoints. Neither this
page nor the archive reconstructs their missing bytes.

## Scientific and transferability limits

A passing route test must actually exercise distinct routes. Fixed seeds and
repeatability within a route do not establish bit identity across reduction
orders, thread counts or native-model histories. Compare saved times and full
states, report the observed difference, and use only the existing scientifically
justified acceptance criterion. If it is insufficient, raise that separately;
this documentation cleanup does not relax a parity gate or introduce a new
tolerance.

Allocation or callback changes also need coverage of the paths they affect.
Changing state-modification notifications or saved-state behavior can change
solver stepping and must not be inferred safe from timing alone. Native GRAM
comparisons require the stated wrapper, data, wind mode and history semantics.

Retest on the target machine and final code revision before making performance
claims. Calibration constants, contention, SMT effects and saturation points
are machine- and workload-dependent. Keep unresolved hypotheses distinct from
measured results and verified source behavior.
