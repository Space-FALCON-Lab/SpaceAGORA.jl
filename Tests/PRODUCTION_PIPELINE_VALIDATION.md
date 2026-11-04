# Production study validation

Validated on 2026-10-18 using the actual
[production runner](../Simulation/Code/README.md). The production format is
schema 2: **1001 trajectory rows**, a **10 s maximum integration step**, exact
event boundaries, and encounter products sampled every 10 s plus boundaries.
The reduced output cadence does not change the scheduler.

## Results

Three complete 100-orbit cases passed:

| Scheduler | Target altitude / inclination | Helpers | Persistent numerical output | Runner elapsed time |
|---|---|---:|---:|---:|
| SMA | 850 km / 1 deg | 50 | 3.10 MB | 41.63 s |
| Inclination | 1050 km / 1 deg | 50 | 2.87 MB | 41.67 s |
| SMA | 1050 km / 1 deg | 300 | 18.82 MB | 289.49 s |

Helpers start at 1000 km and 0 deg; target initial true anomaly is 1 deg.
Times above include runner setup, propagation, recording, and finalization,
but exclude Julia/package startup and the separate fresh-file verifier.
Validation instrumentation was enabled. These are measurements, not runtime
guarantees for other cases or machines.

The 301-spacecraft process took **5 min 2 s** including startup, with
**2.084 GB (1.94 GiB) peak RSS**. Fresh-file verification took **11.17 s**.

### Storage

The 301-spacecraft output contains:

| File | Size |
|---|---:|
| `trajectory.feather` | 12.68 MB |
| `longest_encounter.feather` | 5.25 MB |
| `encounters.feather` | 0.16 MB |
| `laser_firing.feather` | 0.08 MB |
| `firing_range_audit.feather` | 0.01 MB |
| Optional pilot-only `schedule.feather` | 0.63 MB |

This is about **98% smaller** than the previous 784.40 MB dense-output pilot.
The raw temporary 10 s state spool is approximately 0.92 GB and is removed
after encounter products are written. It does not accumulate across cases.
Measured scaling from these pilots suggests roughly **1.8-2.0 GB** for the grid,
not a guaranteed bound. The disk-space guard remains enabled.

## Checks passed

- **391 production fixture assertions:** the 210-case grid, initial phase and
  overlap guard, force/diagnostic derivatives, zero-force reset, multiple
  crossings, empty outputs, 1001-row cadence, and dense encounter summaries.
- **323 interlink assertions** and **232 inclination-scheduler assertions**;
  the stale example include and the single-thread flat-route expectation were
  repaired. The interlink suite also passed with two Julia threads.
- **17 prior analytical recorder fixtures** still pass after moving the shared
  recorder into production code.
- Short SMA and inclination pilots with/without event recording have exactly
  equal saved values and scheduler histories.
- The uninstrumented production CLI also reproduces the compact SMA pilot's
  saved values and schedule exactly and passes fresh-file verification.
- Zero-power and coplanar inclination cases save zero laser diagnostics and no
  firing. The short laser-on target differs from laser-off by **0.11914 m**.
- Fresh processes verify all numerical files, checksums, times, spacecraft IDs,
  initial/final states, event boundaries, range-audit partitions, signed-T sums,
  pair impulse balance, J2 energy/work balance, and axial angular momentum.
- The longest all-pair and target-helper tracks reproduce their stored extrema
  exactly. Their rows use the 10 s grid plus exact boundaries.
- Compact versus dense outputs for all three full-duration pilots have
  **identical final states, diagnostics, firing schedules, and original event
  columns**. All 27 comparison assertions pass.
- Batch tests cover start/finish/failure log entries and elapsed times, duplicate
  process locks, safe completed-case reuse, stop-after-case, and low-disk pause:
  **28 assertions passed**.

### Independent RTN/work validation

The solver integrates the actual applied laser acceleration in the instantaneous
RTN frame and each spacecraft's optical work. A separate composite-Simpson
calculation uses dense Cartesian states and the held link for each completed
step; it does not read the solver's RTN/work states to construct its reference.
The reference is also sampled at all 1001 output times and checked against the
saved Feather columns in a fresh process.

| Full case | Maximum net RTN discrepancy | Maximum signed-T discrepancy | Maximum per-spacecraft work discrepancy |
|---|---:|---:|---:|
| SMA, 50 helpers | 9.42e-14 m/s | 2.61e-10 m/s | 1.39e-7 J |
| Inclination, 50 helpers | 3.56e-14 m/s | 1.89e-9 m/s | 2.81e-8 J |
| SMA, 300 helpers | 1.49e-13 m/s | 3.55e-8 m/s | 2.16e-7 J |

Net RTN tolerance is 1e-8 m/s. The positive/negative T integrals have a kink at
each sign change; their independent-agreement tolerance is **1e-7 m/s**
(0.1 micrometre/s). An earlier 1e-8 m/s check was too strict for the measured
3.57e-8 m/s signed-T difference in the dense 301-spacecraft case. That result
is retained, not hidden or described as an error bound.

Work tolerance is `0.1 J + 1e-8 * maximum(abs(final spacecraft work))`.
Saved signed-T sum/monotonicity checks use 1e-8 m/s without clamping values.
Mechanical verification allows `max(20, spacecraft_count) J + 1e-4*abs(total work)`;
maximum measured energy/work residuals were **3.26, 3.43, and 16.32 J** respectively.
Mass-weighted ECI impulse tolerance is 1e-8 kg m/s and relative axial angular
momentum tolerance is 1e-9.

An early development test exposed RTN derivatives being left nonzero after
firing stopped. The derivatives are now cleared explicitly, with a regression
test. Those two failed development outputs are marked `validation_failed`.
They are not production results.

## Evidence and reproduction

- [Compact test results and logs](results/production_pipeline_compact/)
- [Earlier dense-format full pilots](results/production_pipeline/dense_v1/)
- [Production output](../Simulation/Output/)
- [Readable scenario progress log](../Simulation/Output/simulation_progress.log)

Run commands from the `Tests` directory:

```bash
julia --startup-file=no --threads=1 production_pipeline_tests.jl fixtures
julia --startup-file=no --threads=1 run_short_production_pilots.jl
julia --startup-file=no --threads=1 check_production_on_off.jl
julia --startup-file=no --threads=1 compare_compact_and_dense.jl
julia --startup-file=no test_production_batch_control.jl
julia --startup-file=no assess_production_resources.jl
```

Pilot runners refuse to overwrite outputs; use a fresh output root when
reproducing a run. For example:

```bash
julia --startup-file=no --threads=1 production_pipeline_tests.jl \
  pilot /absolute/path/to/new-output 1050 1 300 gve_sma 100 recorded
julia --startup-file=no --threads=1 ../Simulation/Code/verify_output.jl \
  /absolute/path/to/case-directory
julia --startup-file=no --threads=1 production_pipeline_tests.jl \
  verify-reference /absolute/path/to/diagnostic-report.toml
```

Geometric recording consistency and diagnostic agreement do not imply exact
physical switching or study-wide orbital error bounds. The accepted 10 s
scheduler remains unchanged.
