# Saved paper case versus the new integration

## What was run

The [test script](compare_saved_paper_case.jl) runs the new integration against
the requested [saved case](../../SpaceAGORA.jl/output/paper_plot_mode/h1000km_t1050km/ih0.0deg_it1.0deg/N50/T509707s/gve_sma_e0.0000_nu1.0000/).
All new outputs are in [the result directory](results/saved_paper_h1000_t1050_i1_N50_80orbits/).
No production simulator code or old outputs were changed.

| Setting | Value |
|---|---|
| Target | 1050 km altitude, initial inclination 1 deg, phase 1 deg |
| Helpers | 50 at 1000 km, inclination 0 deg, evenly phased |
| Initial conditions | Exact saved Cartesian states of all 51 spacecraft |
| Mass | 227 kg each, checked against saved masses |
| Duration | 509707.22883989883 s: **80**, not 50, initial target periods |
| Scheduler | `gve_sma`, target index 1, one terminal, target-helper links only |
| Laser | 10 kW, B=100, 200 km selection range; constant BP/c |
| Gravity | Central gravity plus J2; no atmosphere/drag |
| Epoch | 2026-01-01 00:00:00 |
| Solver | Tsit5; relative/absolute tolerances 1e-12 |
| Step caps | 10, 5, and 2.5 s, run sequentially |
| Saved states | Every 10 s plus exact endpoints: 50,972 rows per full run |

This is a comparison case, not a change to the planned 50-orbit grid with the
target fixed at 1000 km.

The old manifest confirms duration, spacecraft count, row count, and file hash,
but does not archive all settings or the source revision. Laser settings,
epoch, gravity choice, and tolerances were inferred from the current old runner.
The initial states and duration are verified matches, not inferred from the
folder name. Full historical configuration equivalence cannot be certified.

## Orbital-element results

These are final-minus-initial **osculating** elements, including J2 effects.
They are not laser-on minus laser-off differences.

| Quantity | Old saved | New, cap 10 s | New, cap 5 s | New, cap 2.5 s |
|---|---:|---:|---:|---:|
| SMA change (m) | 1626.752582 | 1626.203490 | 1627.380818 | 1628.096966 |
| Eccentricity change | 0.0005923270 | 0.0005951516 | 0.0005950725 | 0.0005949712 |
| Inclination change (deg) | -0.001169669 | -0.001171099 | -0.001172765 | -0.001172535 |
| RAAN change (deg) | -34.244739 | -34.322403 | -34.322596 | -34.322676 |

Inclination can decrease in this test because it uses the **SMA** scheduler,
not the inclination scheduler.

The new solver-integrated direct laser SMA contributions are respectively
1632.658419, 1633.838006, and 1634.555652 m. They differ from total SMA change
because the latter also includes gravity/J2 effects.

### Step sensitivity

- Total SMA: 10-to-5 s changes by 1.177328 m; 5-to-2.5 s by 0.716147 m.
- The full three-run SMA spread is 1.893476 m, about 0.12% of the gain.
- Final target position: 10-to-5 s differs by 457.582 m; 5-to-2.5 s by 349.978 m.
- Final target velocity: those differences are 0.451914 and 0.345928 m/s.
- The nominal 10 s cap actually produced a maximum accepted step of 7.952370 s
  and a median of 7.601109 s because of the tight solver tolerances.

The SMA result is relatively insensitive for this case. Do not infer the same
relative accuracy for position, encounter timing, or near-cancelling RTN totals.
Position/velocity differences include accumulated orbital phase differences;
they are not cumulative laser delta-v.
These are measured sensitivities, not a rigorous error bound or a certification
of the full grid. Changing the step cap also changes scheduler update times.

## Why the old and new trajectories differ

Despite similar SMA gains, the new 5 s run ends 19,534.200 m from the old
target position, with a 0.077857 deg RAAN-change difference.

Historical code in the old repository at commit
`b0cf1b52c759fffe5218f73238e65c13d96d78a0` uses `planet.Rp_m` in the J2 term
in `src/environment/gravity/gravity_models.jl`.
The [new formula](<../9_SpaceAGORA.jl-main_full integration/src/models/environment/gravity/central_acceleration.jl>)
uses `planet.Rp_e`, consistent with the coefficient's equatorial-radius convention.

- Equatorial radius: 6,378,136.6 m.
- Mean radius: 6,371,008.366666667 m.
- Nominal J2: 0.00108263.

For diagnosis only, an additional full 5 s run used
`J2_effective = J2 * (R_mean / R_equatorial)^2`.
This reproduces the old J2 force term without changing library code, orbital
radii, or the main runs. It is **not** the gravity setting recommended for the study.

| Metric | Old saved | New, historical-radius diagnostic |
|---|---:|---:|
| SMA change (m) | 1626.752582 | 1627.448882 |
| RAAN change (deg) | -34.244739 | -34.244910 |
| Final position difference from old (m) | -- | 478.200 |
| Final velocity difference from old (m/s) | -- | 0.472785 |

This strongly supports the J2 radius convention as the main source of the
large old/new trajectory disagreement. It does not prove the exact historical
run settings, and some differences remain. Do not change the new J2 model
merely to reproduce an old result.

## RTN and work diagnostics

| Cumulative component (m/s) | Old recorded | New, cap 10 s | New, cap 5 s | New, cap 2.5 s |
|---|---:|---:|---:|---:|
| R | 0.569419007 | 0.568970896 | 0.568536537 | 0.568212197 |
| T | 0.802360495 | 0.803952685 | 0.804534331 | 0.804888259 |
| N | -0.038831683 | -0.039079649 | -0.038976537 | -0.038794821 |

**The new test diagnostics are quadrature estimates, not the proposed RTN ODE states.**
At each accepted step, the test takes the increment of the solver-integrated
ECI laser delta-v and projects it into an estimated midpoint RTN frame.
It estimates work from mass times that increment dotted with midpoint velocity.
Thus it includes the actual applied impulse, including held-link firing beyond
range, without applying any additional force or velocity kick.

Positive/negative along-track increments and target/aggregate-helper work are
saved in a separate accepted-step Feather table. The positive/negative split
is stepwise, not an exact split at every within-step sign crossing.
The planned continuously integrated RTN states remain future production work.

The current old tracker uses accepted-step endpoint accumulation with a range
check; it does not provide the same quadrature/held-link accounting.
The archived file does not identify the exact tracker version.
Treat its recorded RTN values as comparison data, not ground truth.
In particular, the N differences do not decrease monotonically with refinement.

The old trajectory has only 1001 samples, approximately 509.707 s apart.
Its sampled active-helper column cannot recover every firing switch or short
encounter. This test does **not** claim to compare exact encounter durations,
and does not implement the planned encounter or auxiliary range-audit files.

## Saved outputs and checks

Each full-run subdirectory contains:

- `simulation_results.feather`: all-spacecraft ECI positions/velocities,
  integrated ECI laser delta-v, and integrated direct laser SMA contribution.
- `accepted_step_schedule.feather`: selected helper and accepted-step time.
- `accepted_step_diagnostics.feather`: target RTN estimates, signed along-track
  estimates, and target/aggregate-helper work estimates.
- `summary.toml`: settings, endpoint metrics, diagnostic definitions, and hashes.
- Native manifest and CSV output. The current native writer emits the CSV even
  with `save_csv=false`; this test leaves its output intact.

The top-level [comparison CSV](results/saved_paper_h1000_t1050_i1_N50_80orbits/comparison.csv)
and [comparison TOML](results/saved_paper_h1000_t1050_i1_N50_80orbits/comparison.toml)
contain numerical differences. Target-element time series are also saved as CSV.

The smoke test and each of the four full runs passed 227 checks. Checks include
successful completion, exact initial states, final time, 10 s output cadence,
finite states, step caps with floating-point roundoff allowance, one active link
at most, signed along-track consistency, positive laser SMA contribution, and
mass-weighted ECI laser impulse balance below 1e-8 kg m/s.
The comparison command reopens persisted files in a fresh process and checks
their timestamps, schedules, summaries, and diagnostic totals: **548 checks passed**.

An early test assertion incorrectly required exactly two solver saves; there
are duplicate initial saves. Another used too-tight roundoff allowance when
subtracting large timestamps. These test-only assertions were corrected and
the affected runs passed afterward; neither required a dynamics change.

Existing library-wide stale-include/routing test issues were not modified.
These comparison checks are separate from that suite.
Full solver step histories and animation were disabled. The measured 10 s and
5 s processes peaked at about 2.6 GiB RSS; this does not establish memory use
for the planned 301-spacecraft case. Wall-clock timings are optional diagnostics,
not required study outputs.

## Reproduce

From this `Tests` directory, with Julia 1.12.1:

```sh
julia --startup-file=no --threads=1 compare_saved_paper_case.jl smoke
julia --startup-file=no --threads=1 compare_saved_paper_case.jl run 10
julia --startup-file=no --threads=1 compare_saved_paper_case.jl run 5
julia --startup-file=no --threads=1 compare_saved_paper_case.jl run 2.5
julia --startup-file=no --threads=1 compare_saved_paper_case.jl legacy-radius
julia --startup-file=no --threads=1 compare_saved_paper_case.jl compare
```

The script activates the new integration's environment automatically. Simulation
commands refuse to overwrite existing result subdirectories; archive the existing
results before repeating a run. The comparison command can be repeated directly.
