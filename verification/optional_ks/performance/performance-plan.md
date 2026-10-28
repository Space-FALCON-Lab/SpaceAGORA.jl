# Bounded matched-accuracy comparison policy

Prepared before calibration/timing. Do not execute performance work until
current-source correctness, lifecycle and installation gates have passed.
This is a separate benchmark experiment; original acceptance targets and all
eight coarse/medium misses stay unchanged and visible. No performance result
can repair an acceptance failure.

Cases: two_body and J2_plus_drag, identical baseline physical inputs and
0:60:86400 s comparison grid, same independent refined Cartesian references.
Target: both maximum sampled position <=1 m and velocity <=0.001 m/s; mass/heat
and physical endpoint/clock checks still apply. KS starts from the qualified
fine settings. Normal DP8 gets its own calibrated tolerances, not the KS
tolerances. Calibration records all candidate errors, including misses.

DP8 initial candidates: relative tolerance 1e-10,3e-10,1e-9,3e-9,1e-8,3e-8,
1e-7,3e-7,1e-6; absolute tolerance=relative/10; dtmax=600 s in both phases.
Select a passing candidate whose achieved maximum position AND velocity error
are each within a factor of two of KS. Prefer smallest summed squared log error
ratio to KS, not the quickest observed candidate. If none matches, refine the
logarithmic tolerance brackets (geometric midpoint) before timing; preserve all
calibration results. Do not change the common accuracy targets. If matching
remains impossible, report that limitation and do not claim matched accuracy.

Timing: fresh process per case/method, one Julia thread and one BLAS thread,
sequential workers with no other compilation/qualification workers active.
Record actual hardware, runtime, package/source hashes and resolved versions.
Measure package load separately; first full 24h run separately from three warm
runs. Time the normal run_simulation call with return_solution=true, results
files disabled, full physical Cartesian dense solution retained for both
methods. Configuration construction, independent reference, grid evaluation
and evidence-file writing remain outside the measured call. This deliberately
measures the full normal-run route including setup and dense storage, not just
stepping or the controller. Report allocation bytes, GC time, achieved sampled
errors and all repeats.

Use Julia1.12 Base.@timed's wall, compile_time and recompile_time fields. The
installed base/timing.jl688–714 is the API source. Do not subtract aggregate
compilation seconds from wall time. Record first-use wall/compilation separately.
If a nominal warm call still compiles materially (>1% wall), preserve it as a
warm-up and collect at most two additional calls to obtain three stable calls;
otherwise report residual compilation, not a pure execution comparison.

Summarize warm median and range and method ratio only within the declared
accuracy match. A slower KS result is valid evidence. No extrapolation to other
models, outputs, machines, spacecraft counts or near-collision trajectories.

## Pre-execution review corrections

The amended implementation freezes the complete physical source/Projects,
actual included fixture and inputs, all benchmark scripts/policy, active
Project/Manifest with resolved versions, BLAS/settings and source-bound
correctness evidence before execution. It asserts that immutable identity after
each calibration candidate and on every timed worker persistence; it does not
relabel prior rows. The final KS Base.IO annotation-only delta is explicitly
linked to its independent equivalence review and 4,170 whole-package plus316
installed-adapter checks. Calibration recomputes KS fine accuracy on the current
installed source for both benchmark cases; old matrix/source records stay intact.

Candidate starts are written before execution; incomplete endpoints/retcodes are
rejected before interpolation. Exceptions are recorded and abort calibration for
investigation. Missing KS clock values on DP8 are marked not applicable. All
position/velocity error ratios require finite positive errors; no error floor is
introduced. Physics settings precede imports and Julia/BLAS are both single
threaded for calibration and timing. First-use runs use timed eval to include
compilation (and explicitly its eval overhead); warm runs do not.

The final summarizer requires exactly three stable calls per case/method and
checks EVERY actual stable KS/DP8 cross-pair in BOTH position and velocity within
factor2 before computing any median execution ratio. It reports all first-use
and unstable warm-up records separately, warm medians/ranges and raw records.
Failure of this paired gate means no matched-accuracy speed ratio is justified.

## Recorded dormant background process

The initial inventory found an unrelated Julia process (PID55147, parent1)
asleep after more than six days, with CPU total0:29.93 unchanged across
observations and zero reported CPU usage. It belongs to another historical task
and is not terminated here. The first accuracy-only calibration was started
before honoring the stricter inventory gate, then interrupted and preserved as
attempt1; it supplies no accepted measurement.

Independent operational review approved treating that confirmed dormant process
as documented background. All our qualification/install/compiler handles must
be terminal. The process guardian records PID/state/CPU counters and RSS before,
during (2s observation interval) and after each calibration/timed worker, plus
initial memory-pressure/swap context. It allows only that documented sleeping
process with unchanged CPU counter and descendants of the current worker. Any
new unknown Julia process, waking background process or changed CPU counter
invalidates/stops the worker; its data cannot support a ratio. This controls
serial experimental work with documented background; it does not claim an
otherwise empty machine or detect every shorter-than-observation-interval task.
