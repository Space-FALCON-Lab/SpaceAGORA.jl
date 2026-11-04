# Absolute R/N campaign

This is a separate, fresh campaign. The original `Simulation` workspace and its
results are not changed or reused. The interrupted old 700-helper campaign stays
stopped; its four verified cases remain intact.

## Exact production grid

- Helpers: initially circular, 1000 km, inclination 0 degrees, evenly phased.
- Targets: 900, 950, 1000, 1050, 1100 km; inclinations 0, 0.5, 1 degree.
- Helper counts: 1, 100, 200, 300, 400, 500, 600, 700, 800.
- Schedulers: `gve_sma` and `gve_inc`.
- Total: **270 cases**, each lasting **100 initial target periods**.

The initial target anomaly remains 1 degree. J2 gravity, 227 kg spacecraft,
10 kW lasers, magnification 100, 200 km selection range, target-helper links and
one terminal per spacecraft are unchanged. Tsit5 uses relative/absolute
tolerances `1e-12` and a maximum step of **10 seconds**. Scheduling occurs at
initialization and after accepted steps; selected links stay applied between
updates, including recorded out-of-range firing.

## Diagnostics and schema 3

The shared library's existing `interlink_diagnostics=true` option now integrates
positive and negative **R, T and N** laser acceleration contributions. Its
default is still `false`. The four added states start at zero, do not add forces,
and use the instantaneous RTN frame and applied laser acceleration only.
Gravity is never included in these integrals.

Every trajectory row saves the target (spacecraft `N+1`) columns:

```text
target_dv_R_positive_mps    target_dv_R_negative_mps
target_dv_T_positive_mps    target_dv_T_negative_mps
target_dv_N_positive_mps    target_dv_N_negative_mps
```

Negative contributions remain negative. For each component:

```text
signed net           = positive + negative
accumulated absolute = positive - negative
```

No extra absolute-value states are needed. Signed RTN, fixed-ECI laser delta-v,
optical work, all-pair encounters, longest tracks and firing/range audits remain
available. The library still maintains its direct laser SMA diagnostic internally.

There are **1001 uniformly spaced trajectory rows**, including both endpoints.
For a production run, **row 801 is 80 orbits** and **row 1001 is 100 orbits**.
There are no separate 80-orbit runs. Short encounters and extrema still use
accepted-step boundaries and the temporary 10-second Cartesian spool, not the
sparse trajectory. Finalization removes only its own temporary spool.

Schema 2 outputs and older enabled-diagnostic checkpoints lack the new states
and are explicitly rejected. Disabled-diagnostic callers/checkpoints retain
their previous layout. The new verifier checks exact column names/order, units,
target indexing, initial zeros, endpoints, branch signs, monotonicity, net sums,
no-authority cases and the existing physical/encounter audits.

Initial verification limits are `1e-8 m/s` for R/N branch signs and monotonicity,
the existing `1e-7 m/s` allowance for T, and `1e-8 m/s` for each sum-to-net
identity. Small dense-interpolation reversals are reported, not clamped.
Independent pilot validation and measured residuals are stored below.

## Encounter speed, range rate, encounter duration and laser-on distributions

### SMA, 500 helpers: speed at closest approach

[write_sma500_closest_approach.jl](write_sma500_closest_approach.jl) updates only
the encounter-speed section of
[the 500-helper SMA report](../Output/images/for_paper/parameter_pages_sma_500_helpers.md).
It uses one relative-speed magnitude at minimum helper-target separation per
observed encounter, not minimum speed or entry speed. The other three parameters
remain unchanged, including signed range rate **at entry**.

The nine 1000/1050/1100 km, 0/0.5/1 degree cases are replayed with the existing
solver and scheduler into a separate output directory. Original simulation
files are never overwritten. Closest approach is refined by locating
approaching-to-receding range-rate roots in accepted steps (at most 10 seconds)
using dense solver output, comparing them with step and encounter boundaries.
Clipped encounters use the closest point inside their observed interval.
This avoids interpolation of the sparse 1001-row saved trajectory.

Before publishing a case, the script requires bit-identical saved trajectories
and scheduler fingerprints, checks both boundary states of every encounter,
and verifies that refined separations do not exceed the original sampled minima.
Interior minima must have near-zero range rate. Per-event results, full-precision
summaries and verification/checksum records are retained under
`Output/encounter_evaluation/sma500_closest_approach/`.

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/write_sma500_closest_approach.jl
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_sma500_closest_approach.jl
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_sma500_closest_approach.jl --campaign
```

The campaign check validates all four report tables against their event records
and independently compares closest-approach states with each original
10-second longest-encounter track.

For individual extraction use `--case ALTITUDE_KM INCLINATION_DEG`; use
`--report-only` to regenerate the speed table from all nine verified case CSVs.
Replay directories must not already exist; the runner refuses to overwrite them.

### Full campaign distributions

[evaluate_encounters.jl](evaluate_encounters.jl) evaluates all 270 saved, verified
100-orbit cases without rerunning propagation. It checks settings and checksums,
selects **target-helper pairs only**, and exports event CSVs, a case-summary CSV
with source provenance, and distribution statistics grouped by initial target
inclination, initial target altitude, helper count and scheduler. It then reloads
the CSVs to generate [summary.md](../Output/encounter_evaluation/summary.md).

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/evaluate_encounters.jl
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_encounter_evaluation.jl --campaign
```

Speed and signed range rate are **sampled per-encounter minima/maxima**, not
instantaneous time-weighted distributions. Durations distinguish geometric
encounters, uninterrupted actual firing (including held out-of-range firing),
and summed within-range firing per encounter. Main statistics are event-weighted;
boundary-clipped events are flagged and separate unclipped-only statistics are
included. Zero-firing cases remain in case summaries. See the generated report
for definitions, units, counts, percentiles and interpretation limits.

Options: `--input PATH`, `--output PATH`, and `--report-only` (regenerate Markdown
from existing distribution and case-summary CSVs).

### Histograms and empirical cumulative distributions: 1000/1050/1100 km

[plot_encounter_distributions.jl](plot_encounter_distributions.jl) reads the
evaluation CSVs and validates the complete selected **162-case** grid. It writes
four figures, each as PNG and PDF, to
[encounter_distributions_1000_1050_1100km](../Output/images/encounter_distributions_1000_1050_1100km):
encounter speed, signed range rate, encounter duration, and uninterrupted laser-on
duration. Each figure has altitude columns (1000, 1050, 1100 km), histograms above,
and exact empirical cumulative distribution functions (ECDFs) below.

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/plot_encounter_distributions.jl
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_encounter_distribution_plots.jl
```

Events are pooled over initial inclination, helper count and scheduler, using the
full 100-orbit runs and target-helper pairs only. Speed and range rate show
separately normalized sampled minimum/maximum distributions. The 40 histogram
bins are shared across altitudes within each figure; heights are **percent of
events in a bin**, not probability density. The uninterrupted laser-on figure
separates `gve_inc` and `gve_sma` into two color-matched histogram/ECDF series,
each normalized within that scheduler and altitude. Scheduler event counts are
shown in the legends; zero-firing cases do not create firing intervals.
Durations use equal logarithmic-width
bins and logarithmic x axes to retain the full long tails. ECDFs use every event,
without smoothing, downsampling or binning. All-observed distributions include
boundary-clipped intervals; the geometric-duration ECDF also shows an unclipped-only sensitivity
curve, not a censoring-corrected lifetime estimate. Counts and clipping are labeled.
The 27 zero-firing cases remain in case coverage but contribute no fictional
zero-duration firing interval. Laser-on intervals include held out-of-range firing,
not just the within-encounter accumulated laser-on quantity.

The folder also contains `histogram_bins.csv`, `distribution_quantiles.csv`,
`selected_cases.csv`, and `input_checksums.csv` for provenance. Options:
`--input PATH` (evaluation CSV directory), `--output PATH` (figure directory).
These standalone plots do not overwrite the existing orbital/delta-v figures.

#### Exactly one observation per geometric encounter

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/evaluate_encounters.jl
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/plot_encounter_distributions.jl --per-encounter
```

This writes four additional PNG/PDF figures and plot-data CSVs in the
[`per_encounter` subfolder](../Output/images/encounter_distributions_1000_1050_1100km/per_encounter),
preserving the original plots. Each altitude has exactly one histogram and one
ECDF for each of these quantities:

- Relative speed **at encounter entry**, from the saved start boundary velocities.
- Signed range rate **at encounter entry**, from the saved start boundary positions/velocities.
- Observed geometric encounter duration.
- **Total laser-on within that geometric encounter**, summed across firing segments,
  including zero for every encounter that received no firing.

All four use the same encounters, pooled across schedulers, inclinations and
helper counts. This is not a time-average: for encounters already in range at
simulation start, the entry metrics use the initial saved state. New CSV columns
`entry_relative_speed_mps` and `entry_range_rate_mps` retain these definitions;
the original sampled extrema remain available.
The per-encounter laser axis uses `log10(1 + t)` spacing with ticks labeled in
seconds, so **zero is retained**, not dropped or replaced by an arbitrary small
duration. Its ECDF includes the jump at zero. No unclipped-only overlay is added
in this mode, ensuring a single empirical curve per subplot.

## Four-page parameter min/max report

[write_selected_altitude_overview.jl](write_selected_altitude_overview.jl) generates a
compact all-scenario [summary](
../Output/encounter_evaluation/selected_altitude_summary.md) for target altitudes
1000, 1050 and 1100 km. It first provides a combined table pooling every event
across the three altitudes, followed by per-altitude tables. Each table provides
Max/Avg/Min rows for entry speed, signed entry range rate, encounter duration,
and uninterrupted laser-on duration separated by scheduler. All inclinations and helper counts are included;
encounter metrics pool both scheduler-generated trajectories. Means are
event-weighted, so the report explicitly cautions against reading them as
equal-scenario means or scheduler-effect estimates. The report includes event and
boundary-clipping counts, plus CSV provenance. Combined statistics are also
exported in `selected_altitude_combined_statistics.csv`.
The combined section also repeats its table below with a compact presentation
rounding: five significant figures for speed, range rate and laser-on duration,
and six for encounter duration. Five significant figures are used for range rate
as the requested precision was unspecified; the assumption is labeled in the report.
The first combined table instead prints round-trip Float64 values; the selected-altitude
CSV files retain full precision as well.

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/write_selected_altitude_overview.jl
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_selected_altitude_overview.jl
```

[write_encounter_parameter_pages.jl](write_encounter_parameter_pages.jl) generates
[parameter_pages.md](../Output/encounter_evaluation/parameter_pages.md) and its
full-precision `parameter_pages_min_max.csv` (now including `avg`) and source-checksum CSV.
Each of four page sections has a conclusion and one combined table matching the
requested sketch: grouped altitude headers (1000, 1050, 1100 km), Min/Avg/Max
subcolumns, stacked SMA/INC scheduler blocks with vertical labels, and inclination
rows 0, 0.5, 1.0 deg in each block. A diagonal corner labels inclination/statistic.
Tables are generated as PNG/PDF and embedded using standard Markdown image links,
avoiding unsupported merged HTML cells. Keep the images beside the Markdown.
Avg is the arithmetic mean over events in each cell, not an equal-case average
or a time average within an encounter. Empty groups show N/A for all three statistics.
The four numbered page sections are separated by horizontal rules; physical
print pagination depends on the Markdown viewer/exporter.

Page order is entry speed, entry signed range rate, geometric encounter duration,
and uninterrupted laser-on duration. This explicitly interprets the repeated
"encounter duration" in the page-1/page-3 request as speed on page 1. Min/Max pool
all events over the nine helper counts in each cell. Clipped events remain
included; no-firing cells are N/A. These are not within-encounter speed/range-rate
extrema or accumulated per-encounter laser-on values.

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/write_encounter_parameter_pages.jl
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_encounter_parameter_pages.jl
```

Options: `--input PATH` (evaluation CSV directory), `--output PATH` (report
directory). Report tables and conclusions are built by reloading the exported
table CSV. No simulation or existing plot output is changed.

## Plot the completed campaign

From the project root:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/plot_all_results.jl
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_result_plots.jl
```

[plot_all_results.jl](plot_all_results.jl) reads the new campaign once and generates
**30 figures, each in PNG and PDF**, for both 80 and 100 initial target periods:
18 full-altitude figures, 8 selected-altitude figures and 4 paired-branch figures.
It does not rerun simulations or read the old campaign. The 80-orbit result is the
exact saved row 801 of each 100-orbit trajectory, not an interpolated endpoint.
All 270 verified cases must be present, with matching settings and trajectory
checksums. All saved RTN rows are checked for finite values, initial zeros, branch
signs, monotonicity and positive-plus-negative consistency, using the verifier's
tolerances. Small tolerated branch reversals are reported without clamping.

The output folders under [../Output/images](../Output/images) are:

| Folder | Contents |
|---|---|
| `orbital_elements/1_seperate_plots_by_sceduler_type` | Four OE rows, separate SMA and INC scheduler grids |
| `orbital_elements/2_separate_plots_by_N-sceduler_h-i` | Pooled OE scatter views versus helpers and target altitude; no averaging or jitter |
| `orbital_elements/3_combined_plots` | Four OE rows with both schedulers |
| `delta_v_RTN` | Three signed-net rows: R, T, N |
| `delta_v_RTN/1_absolute_RTN` | Three accumulated-absolute rows: R, T, N |
| `delta_v_RTN/2_signed_RT_positive_negative_N` | Four rows: signed R, signed T, positive N, negative N |
| `delta_v_RTN/3_positive_negative_RTN` | Six rows: R+, R-, T+, T-, N+, N- |

OE CSVs are in `orbital_elements`; each RTN family has its own endpoint CSVs.
Every endpoint CSV includes scheduler, altitude, inclination, helper count, orbit
window, exact saved time, 1-based row and source trajectory path. Each family also
includes a 15-row `continuous_access_thresholds.csv`.

### Selected-altitude variants

These additional folders contain only target altitudes **1000, 1050, 1100 km**,
in that left-to-right order, for both 80 and 100 orbits:

| Folder | Same layout as |
|---|---|
| `orbital_elements/4_combined_1000_1050_1100km` | OE folder 3: combined four-row grids |
| `delta_v_RTN/5_absolute_RTN_1000_1050_1100km` | RTN folder 1: accumulated absolute R/T/N |
| `delta_v_RTN/6_signed_RT_positive_negative_N_1000_1050_1100km` | RTN folder 2: signed R/T and N+/N- |
| `delta_v_RTN/7_positive_negative_RTN_1000_1050_1100km` | RTN folder 3: six positive/negative branch rows |

Each new folder includes PNG/PDF figures, two endpoint CSVs with **162 cases per
window**, and a **9-row** continuous-access threshold CSV. Y-axis limits retain
the full-campaign ranges for direct comparison with the existing plots.
The canvas keeps enough width for the shared legend even with three columns.
No RTN folder 4 is created; the requested numbering 5/6/7 is preserved.

To generate only these new variants without overwriting the original figures:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/plot_selected_altitudes.jl
```

[plot_selected_altitudes.jl](plot_selected_altitudes.jl) accepts the same
`--input`, `--output`, and `--orbits` options as the all-results command.

### Folder 8: paired positive/negative branches, separate schedulers

[plot_paired_dv_RTN.jl](plot_paired_dv_RTN.jl) generates four figures in
`delta_v_RTN/8_paired_RTN_by_scheduler`: SMA and INC separately, each at 80 and
100 orbits. Each figure has three rows (R, T, N), with both positive and negative
integrals on each axis, and columns 1000, 1050, 1100 km from left to right.
These are separate contributions, not their sum or absolute value. Each component
uses a common range covering both branches over the full campaign, both schedulers
and both requested windows. A horizontal zero line separates the signs.

The **MANUAL APPEARANCE SETTINGS** section at the beginning of the script defines
`LINE_STYLES[(scheduler, inclination, component, branch)]` for every line.
Scheduler keys are `"gve_sma"` / `"gve_inc"`, inclinations are `0.0`, `0.5`, `1.0`,
components are `:R`, `:T`, `:N`, and branches are `:positive` / `:negative`.
Edit the default settings or use the commented single-line override example to
change color, `:solid` / `:dash` / `:dot` / `:dashdot`, width, marker shape, size,
fill and outline. Legend samples automatically follow each line's settings.
By default positive branches use filled circles, negative branches use white-filled
circles with colored outlines; SMA is solid and INC dotted. Threshold colors/styles
are configured separately so they remain distinguishable from data curves.

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Simulation_absolute_RN/Code/plot_paired_dv_RTN.jl
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_paired_RTN_plots.jl
```

The standalone command only writes folder 8 by default. It supports `--input`,
`--output` (the folder containing these four figures), and `--orbits`.
Each figure has a matching **81-case CSV**, with signed branches and source
provenance; the folder also includes a nine-row continuous-access threshold CSV.

Combined grids retain the original six colors, solid SMA curves, dotted INC
curves and grid lines. Dashed vertical lines mark
`N_cont = ceil(pi / alpha_min)` for each initial inclination, where
`alpha_min = acos((rh^2 + rt^2 - Lmax^2) / (2*rh*rt*cos(delta_i)))`.
This is a fixed-circular, range-only geometric estimate, not a guarantee under
J2, changing geometry, scheduler constraints or held-link firing. Thresholds
outside the simulated helper range are annotated instead of extending the axis.

Each delta-v quantity has one shared y-axis range across every altitude column,
both schedulers and all requested orbit windows (80 and 100 by default). Signed
net, accumulated absolute, positive and negative branches are scaled separately
for each R/T/N direction; the same quantity keeps its range in different plot
families. Bounds are rounded outward to readable values and include zero and all
endpoints. Nonnegative groups start at zero and nonpositive groups end at zero;
values are never clipped or clamped.

Each orbital-element change (`da_m`, `de`, `di_deg`, `draan_deg`) likewise shares
one y-axis range across altitudes, both schedulers and all requested windows.
Separate-scheduler grids, combined grids and both pooled scatter views use the
same limits for the same quantity. OE bounds are rounded outward around all
finite endpoints, without forcing zero into a range such as the negative RAAN
changes. Undefined equatorial RAAN values are excluded from the range calculation
and remain omitted from the plots.

Accumulated absolute delta-v is **positive minus negative**, not `abs(signed net)`.
Negative branches stay negative. All delta-v quantities are laser-only integrals
in the instantaneous target RTN frame, in m/s. OE changes include J2 and laser;
RAAN is unwrapped through every saved row and omitted for equatorial targets.

Optional arguments are `--input PATH`, `--output PATH` and `--orbits 80,100`.
For the all-results command, output defaults to this workspace's `Output/images`.
[plot_orbital_elements.jl](plot_orbital_elements.jl) and
[plot_pooled_orbital_elements.jl](plot_pooled_orbital_elements.jl) also remain usable
as standalone plotting commands.

## Isolation, validation and reproducibility

- Production data: [../Output/paper_plot_mode](../Output/paper_plot_mode).
- Production controls: [../Output/_batch_absolute_RN](../Output/_batch_absolute_RN).
- New tests/scripts: [../../Tests/absolute_RN](../../Tests/absolute_RN).
- Pilot outputs: [../../Tests/absolute_RN/PilotOutput](../../Tests/absolute_RN/PilotOutput).
- Reports/benchmarks: [../../Tests/absolute_RN/results](../../Tests/absolute_RN/results).
- Pre-change source/configuration snapshot:
  [source_snapshot_pre_absolute_RN](../../Tests/absolute_RN/results/source_snapshot_pre_absolute_RN).
  Its archive, SHA256 and per-file manifest include the current local source
  changes; no clean Git repository is assumed. Package versions and asset
  manifests are preserved; external package caches/binary assets are not copied.

Tests cover the ODE RHS and disabled callers, zero power, target indexing,
independent dense-state quadrature of net/positive/negative/absolute R/T/N,
physical agreement with the frozen pre-change library, both schedulers and a
large helper population. Old/new adaptive steps need not be bitwise identical.
Launch is conditional on successful tests, fresh-file pilot verification and
resource checks.

From the project root, rerun the focused suite with:

```bash
julia --startup-file=no --threads=1 --gcthreads=1 Tests/absolute_RN/test_absolute_RN.jl
```

[pilot.jl](../../Tests/absolute_RN/pilot.jl) refuses to overwrite a pilot label.
[run_benchmarks.sh](../../Tests/absolute_RN/run_benchmarks.sh) documents the
serial old/new benchmark protocol; choose new labels for a new measurement.
[validate_pilots.jl](../../Tests/absolute_RN/validate_pilots.jl) performs fresh-file,
reference-integral and physical comparisons. [preflight.jl](../../Tests/absolute_RN/preflight.jl)
checks the complete evidence set, source snapshot, old N700 checksums and
whole-campaign resource capacity before creating the validation report.

The [validation report](../../Tests/absolute_RN/results/VALIDATION_REPORT.md)
contains the measured runtime and resource estimates.
[confirm_launch.jl](../../Tests/absolute_RN/confirm_launch.jl), run while the batch
is active, waits for the first production verification and checks continued
execution in tmux. It writes the
[launch confirmation](../../Tests/absolute_RN/results/launch_confirmation.toml)
without interrupting the campaign.

## Launch, monitor, stop and resume

Run from the project root:

```bash
bash Simulation_absolute_RN/Code/launch_study.sh start
bash Simulation_absolute_RN/Code/launch_study.sh status
tail -f Simulation_absolute_RN/Output/simulation_progress_absolute_RN.log
tmux attach -t oracle-study-absolute-rn
bash Simulation_absolute_RN/Code/launch_study.sh stop-after-case
bash Simulation_absolute_RN/Code/launch_study.sh resume
```

Detach from tmux with `Ctrl-b`, then `d`. SSH may disconnect; the lab computer
must remain powered on and awake.

Campaign ID: `absolute_rn`. The batch runs one single-threaded simulation and
then one fresh verifier at a time; Julia/GC, BLAS and OpenMP thread limits are
explicit. Its own lock prevents duplicate batches. Stop-after-case finishes
simulation and verification before stopping. Resume re-verifies compatible
complete cases; it never treats interrupted output as complete.

The 5 GiB reserve, pre-case spool/trajectory capacity estimate and in-run disk
guard remain enabled. No old data is automatically deleted to free space.
Runtime-source changes stop the batch before the next case. Standalone
`Code/plot_*.jl` scripts remain excluded from that guard; other simulation,
configuration, verification and launch changes remain guarded.
