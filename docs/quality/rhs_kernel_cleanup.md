# Shared per-spacecraft dynamics kernels

`src/orchestration/simulation/engine/dynamics_rhs.jl` owns one private
per-spacecraft operation for each existing physics partition: full, slow,
implicit atmosphere, explicit remainder and fast control. Serial and satellite
batch loops call the same operation for that partition. The partitions remain
separate because their physics and derivative assignments differ.

Execution-plan selection, shared-sample prefill, batch sizing, synchronization
and interlink application remain at their existing outer levels. Early returns
skip only the current spacecraft. Bounds context propagates through the private
helper so extraction does not add checks to the established inbounds route.
Public names and module ownership are unchanged.

## Qualification of the extraction

The new parity test exercises all five partitions with orientation enabled and
disabled, including an inactive spacecraft and nonzero atmospheric heat rates.
A four-thread run passed 88 assertions. Focused wheel, split-solver, external
propagation, articulation and compliant-attachment regressions also passed.
The original and extracted implementations produced identical derivatives in
20 before/after combinations of partition, orientation and serial/batch route.

Allocation measurements use the callback environment snapshot captured by normal
simulation setup. A hand-built parameter fixture without that snapshot instead
repeatedly parses the environment; its initial single-call measurement showed
an intermittent difference and was investigated rather than treated as a pass.
Allocation profiling of the original and extracted code recorded identical
allocation counts and sizes. In a separate run without other local Julia jobs,
93 alternating paired batches covering the investigated threaded full, slow and
fast-control routes each measured 100 allocation samples and 20,000 timed calls.
Every pair had matching allocation minima and maxima. Paired timing-ratio medians
were 0.966, 1.017 and 0.991 respectively, with variable individual timings.

These checks support unchanged tested derivatives and no consistent timing
regression in the investigated routes. They do not establish a general speedup,
zero allocations, or full campaign performance acceptance. Independent review
confirmed the extracted bodies match the originals and found no remaining
adoption blocker after the targeted measurements.
