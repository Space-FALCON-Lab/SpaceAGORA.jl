# Parallelization measurement and validation

The original investigation, with its measurements and limitations, is kept
at this same relative path in the lab's private research archive (see
[Record a result](../../../docs/src/user/studies_benchmarks.md#record-a-result)).
This page keeps the measurement rules the parallelization harness follows. The
fixes, timing tables and open hypotheses of that investigation stay in the
archive; they are neither a list of current defects nor evidence about how new
code performs.

## Timing scope

- Pass `--warmup` on to each worker process. A point repeated in a new
  subprocess still pays for compilation unless that subprocess warms up first.
- For steady-state throughput, start and warm up the Distributed pool before
  timing. Report process start-up or cold compilation separately when those
  costs are part of what a user would see.
- Record physical cores, SMT threads, CPU affinity and other load. State the
  resource ladder explicitly; `SPACEAGORA_PPC_PHYSICAL_CORES` overrides the
  harness's physical-core detection.
- Record the outer concurrency that actually ran.
  `SPACEAGORA_OUTER_PARALLEL_ACTIVE` means an enclosing split owns the resource
  budget; a parallel-sounding mode name does not make a single constellation
  solve outer-parallel.
- Static profiles pin their routes. To measure the shipped behavior of an
  adaptive profile, keep the calibration it ships with; turning calibration off
  answers a different question, so label it.

The timing and modes are implemented in
`benchmarks/studies/parallelization_performance/{cli,modes,execution}.jl`.
Record the revision and the relevant environment with every result, including
the calibration and corrections store and whether the run started cold or
reused it.

## Compare equivalent work

Check the workload before you read a scaling curve: the force and actuator
stack, mission length, whether output is on, native calls, cache path, sample
count, the worker width that actually ran and the inner budget. A setting that
nothing on the exercised path reads is not an experimental axis. In
particular, look at the saved output rather than assuming that
`num_steps_to_save` changes the solver's steps or the saving callback.

The P6 density cases kept their names while their definitions changed; see
[GRAM benchmark configuration](../../../docs/architecture/gram_thread_scaling.md#benchmark-configuration-and-interpretation).
Do not combine old and new rows because they share a case name. In the same
way, a benchmark exclusion recorded in an old investigation must be checked
again before it is lifted or reported as a current failure.

The paper's phase definitions and case guide, including the larger B7 to B14
families and the runner's preview mode, are in the
SpaceAGORAParallelizationStudies repository.

## Allocation and routing investigations

Knowing where allocations happen is not the same as showing a scaling
bottleneck: removing allocations does not prove that they caused a thread
ceiling. Measure the right-hand-side work, scheduling, reduction, callbacks and
storage, and garbage collection separately.
[Heap contention](../heap_contention/FINDINGS.md) describes the profiling and
state-dump tools kept here.

For nested work, limiting the inner width and changing the algorithm are
separate decisions. A real serial baseline uses an inner budget of one where
needed, instead of quietly starting parallel work from a process that happens
to have several threads. Record `SPACEAGORA_INNER_THREAD_BUDGET` and
`SPACEAGORA_RHS_CALIBRATE` with the route. The harmonics spin barrier defaults
to `SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=auto`, which compares its two arms
within a solve; `1` and `0` force one arm for controlled measurements. Record
the router's actual choice and widths instead of inferring them from a profile
label.

## Repeating an evaluation

From the repository root, pick a supported case and a resource budget that
fits the machine. This example asks for a warm-up explicitly:

```bash
SPACEAGORA_PPC_PROFILE=full SPACEAGORA_PPC_CASES=independent_1sat_1hr \
SPACEAGORA_PPC_MODES=serial,outer_threads SPACEAGORA_PPC_THREADS=1,4 \
  julia --project=examples benchmarks/studies/parallelization_performance.jl full --warmup=1
```

The usage notes at the top of each tool cover the allocation probes and raw
state dumps. The original investigation also used two scratch scripts,
`steady.jl` and `allocprobe.jl`, that were never tracked; they are not part of
the repository and cannot be recovered.

## Scientific and transferability limits

A route test that passes must actually run different routes. Fixed seeds and
repeatability within one route do not show that results are bit-identical
across reduction orders, thread counts or native-model histories. Compare the
saved times and full states, report the difference you see, and judge it only
by the existing, scientifically justified acceptance criterion. If that
criterion is not enough, raise it separately; do not loosen a parity gate or
add a tolerance on the way.

Changes to allocations or callbacks also need tests on the paths they touch.
Changing state-modification notifications or what is saved can change the
solver's steps, so timing alone cannot show such a change is safe. Native GRAM
comparisons need the stated wrapper, data, wind mode and history handling.

Retest on the target machine at the final revision before claiming anything
about performance. Calibration constants, contention, SMT effects and
saturation points depend on the machine and the workload. Keep open hypotheses
apart from measured results and from behavior confirmed in the source.
