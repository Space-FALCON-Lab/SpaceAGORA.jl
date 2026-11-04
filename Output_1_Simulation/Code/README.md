# ORACLE study runner

Each study case runs for 100 initial target periods, with a maximum
integration step of 10 s, and 1001 evenly spaced trajectory rows including both endpoints.
The completed campaigns cover 462 cases: seven target altitudes and eleven
helper counts through 800, with helpers initially at 1000 km.
The new 700-helper campaign adds 42 cases without rerunning the completed cases.
Target settings, helper counts, and schedulers follow the
[study plan](../../Plans/Combined_Maneuverability_and_Encounter_Study_Plan.md).

## Additional 700-helper campaign

`helpers700` runs **42 new cases**: 7 target altitudes x 3 target inclinations
x 2 schedulers, with **700 helpers plus one target**. Targets start at 850, 900,
950, 1000, 1050, 1100, or 1150 km; helpers start at **1000 km**. Inclinations,
100-orbit duration, 10 s step cap, and all recording settings are unchanged.
It does not rerun, reverify, or modify the 462 completed cases.

From this directory:

```bash
bash launch_study.sh start ../Output helpers700
bash launch_study.sh status ../Output helpers700
tail -f ../Output/simulation_progress_helpers700.log
bash launch_study.sh stop-after-case ../Output helpers700
bash launch_study.sh resume ../Output helpers700
```

Results use `N700/` under the existing `../Output/paper_plot_mode/h1000km_t.../`
hierarchy. The separate tmux session is `oracle-study-helpers700`; its progress
log is `simulation_progress_helpers700.log`, with status, detailed logs, lock,
and stop flag in `_batch_helpers700/`. SSH can disconnect; the lab computer
must stay powered on and awake.
Direct CLI: `julia run_grid.jl ../Output helpers700`.

The shared [extension tests](../../Tests/extended_helper_tests.jl) accept
`helpers700`; validation evidence is retained under
[`Tests/results/helper_extension700/`](../../Tests/results/helper_extension700/).
Two-orbit 700-helper SMA and inclination pilots passed independent RTN/work
and fresh-file checks. The inclination pilot matched the uninstrumented
production CLI exactly. Both pilots used about 2.30 GB peak RSS; this is not
a full-duration memory bound.
Every full-duration case is verified in a fresh process before the next starts.
Existing plots remain scoped to the completed 462-case grid and report the
exclusion of 700-helper cases.

## Additional 600/800-helper campaign

`helpers600800` runs **84 new cases**: 7 target altitudes x 3 target inclinations
x 2 helper counts x 2 schedulers. Targets start at **850, 900, 950, 1000, 1050,
1100, or 1150 km**; helpers start at **1000 km**. There are 601 or 801 spacecraft
including the target. All other simulation and recording settings are unchanged.
The campaign does not rerun, reverify, or modify the 378 completed cases.

From this directory:

```bash
bash launch_study.sh start ../Output helpers600800
bash launch_study.sh status ../Output helpers600800
tail -f ../Output/simulation_progress_helpers600800.log
bash launch_study.sh stop-after-case ../Output helpers600800
bash launch_study.sh resume ../Output helpers600800
```

Results use `N600/` and `N800/` under the existing
`../Output/paper_plot_mode/h1000km_t.../` hierarchy.
The campaign has its own `oracle-study-helpers600800` tmux session,
`simulation_progress_helpers600800.log`, and `_batch_helpers600800/` status,
detailed logs, lock, and stop flag. It survives SSH disconnection; the lab
computer must stay powered on and awake.
Direct CLI: `julia run_grid.jl ../Output helpers600800`.

The shared [extension tests](../../Tests/extended_helper_tests.jl) accept
`helpers600800`; validation evidence is retained under
[`Tests/results/helper_extension600800/`](../../Tests/results/helper_extension600800/).
Two-orbit 600-helper SMA and 800-helper inclination pilots passed independent
RTN/work and fresh-file verification. The 800-helper pilot also matched the
uninstrumented production batch CLI exactly. The short pilots used at most
2.46 GB peak RSS; this is a measurement, not a full-duration memory bound.
Every full-duration case is verified in a fresh process before the next starts.

## Additional 900/1100 km target-altitude campaign

`targets9001100km` runs **108 new cases**: 2 target altitudes x
3 target inclinations x 9 helper counts x 2 schedulers, with helpers at **1000 km**.
Helper counts are **1, 50, 100, 150, 200, 250, 300, 400, and 500** at each altitude.
Target inclinations stay at 0, 0.5, and 1 degree. All other simulation and recording
settings are unchanged. This campaign does not rerun, reverify, or modify the
270 completed cases at target altitudes 850, 950, 1000, 1050, and 1150 km.

From this directory:

```bash
bash launch_study.sh start ../Output targets9001100km
bash launch_study.sh status ../Output targets9001100km
tail -f ../Output/simulation_progress_targets9001100km.log
bash launch_study.sh stop-after-case ../Output targets9001100km
bash launch_study.sh resume ../Output targets9001100km
```

Results go under `../Output/paper_plot_mode/h1000km_t900km/` and
`../Output/paper_plot_mode/h1000km_t1100km/`, using the same remaining folder/file layout.
The separate tmux session is `oracle-study-targets9001100km`; its progress log is
`simulation_progress_targets9001100km.log`, with status, detailed logs, lock,
and stop flag in `_batch_targets9001100km/`.
The lab computer must remain powered on and awake; the SSH client can disconnect.
Direct CLI: `julia run_grid.jl ../Output targets9001100km`.

Altitude-grid and isolation tests are in
[target_altitude_tests.jl](../../Tests/target_altitude_tests.jl).
Validation evidence is retained under
[`Tests/results/target_altitudes/`](../../Tests/results/target_altitudes/).
Two-orbit pilots (900 km target / 400 helpers / SMA and 1100 km target /
500 helpers / inclination), both with helpers initially at 1000 km, passed
independent RTN/work and fresh-file verification. The uninstrumented batch
CLI reproduced both pilots' saved values and firing schedules exactly.
Each full case is verified in a fresh process before the batch continues.

The earlier `helpers9001100km` campaign used the wrong altitude interpretation
and is stopped. Its 35 verified cases, one interrupted case, and separate log
are retained, but are **not part of the intended study**. Do not resume it for
this study or mix its results with the corrected target-altitude campaign.

## Additional 400/500-helper campaign

The `helpers400500` campaign runs only **60 new cases** (5 altitudes x
3 inclinations x 2 helper counts x 2 schedulers), with the same 100-orbit
duration and output format. It does not rerun, reverify, or modify the original
1-300-helper results. New cases use the existing output hierarchy under
`N400/` and `N500/` (401 and 501 spacecraft including the target).

From this directory:

```bash
bash launch_study.sh start ../Output helpers400500
bash launch_study.sh status ../Output helpers400500
tail -f ../Output/simulation_progress_helpers400500.log
bash launch_study.sh stop-after-case ../Output helpers400500
bash launch_study.sh resume ../Output helpers400500
```

This campaign has its own `oracle-study-helpers400500` tmux session,
`simulation_progress_helpers400500.log`, and `_batch_helpers400500/`
status, detailed logs, lock, and stop flag. As with the original batch, it
survives SSH disconnection but not machine shutdown or sleep.
Direct CLI: `julia run_grid.jl ../Output helpers400500`.
Omitting the campaign selector still selects only the original 210-case grid.
The orbital-element plot command below combines these two campaigns, the
900/1100 km target-altitude campaign, and the 600/800-helper campaign (462 cases).

Extension validation is retained under
[`Tests/results/helper_extension/`](../../Tests/results/helper_extension/):
60-case selection and original-setting checks, separate-log/lock tests, and
two-orbit 400-helper SMA / 500-helper inclination pilots. Both pilots passed
independent RTN/work and fresh-file verification. Every full-duration case is
also verified by the batch before it continues.

## Start and monitor

From this directory:

```bash
bash launch_study.sh start
bash launch_study.sh status
tail -f ../Output/_batch/grid.log
```

Open `Output/simulation_progress.log` for the readable scenario log. It appends
UTC start/finish timestamps, the scenario folder, and wall-clock elapsed seconds
for simulation and verification separately. It also records reuse, failures,
and disk-space pauses. `SCENARIO_VERIFIED` means both stages succeeded.
Entries are flushed immediately and preserved across batch restarts.

`start` uses a detached `tmux` session named `oracle-study`, so an SSH disconnect
does not stop the batch. Alternatively, a process manager can run
`bash launch_study.sh run` in its own detached session. Neither method survives
a machine shutdown by itself.

Each case and its fresh-file verification use separate, single-threaded Julia
processes. Logs are in `Output/_batch/logs/`; the active case's `progress.toml`
records simulation progress and memory use.

```bash
bash launch_study.sh stop-after-case
bash launch_study.sh resume
```

Stopping this way finishes the current case and verification. Resume verifies
and reuses completed cases with matching settings and simulation-source checksums.
Verifier/launcher-only changes do not require resimulating trajectories; existing
files are always checked again by the current verifier.
Standalone `Code/plot_*.jl` scripts are also excluded from the batch source-change
guard and reuse compatibility checks: they read results but are not loaded by
simulation or verification. Their hashes are still recorded in metadata.
Interrupted/failed cases are retained and cause an explicit stop; they are not
silently overwritten or resumed partway through an orbit. Investigate them and
move the affected case directory to an archive before retrying.

The batch pauses before a new case if free space cannot cover its temporary
dense-state spool, final trajectory allowance, and a 5 GiB reserve.
Free space and use `resume` to continue.
An unexpected drop below 5 GiB during a case stops that case explicitly and
retains its partial output. Failed verification or changed source code also
stops the batch. It never deletes old results. Keep the machine awake.

## Output layout

```text
Output/
  simulation_progress.log
  _batch/
    grid.log
    status.toml
    logs/
  paper_plot_mode/
    h1000km_t1050km/
      ih0.0deg_it1.0deg/
        N300/
          T637134s/
            gve_sma_e0.0000_nu1.0000/
              trajectory.feather
              encounters.feather
              laser_firing.feather
              firing_range_audit.feather
              longest_encounter.feather
              run_info.toml
              progress.toml
              verification.toml
```

Folder duration is rounded for readability; metadata stores the exact duration.
Files are written as `.partial` until closed successfully. A run is usable only
after `run_info.toml` reports `status = "verified"`.

The trajectory is a compressed Feather V2/Arrow IPC file with 512-row batches,
not an unbounded in-memory table. Read it with `Arrow.Table` or
`pyarrow.feather.read_table`; use `Arrow.Stream` to process one batch at a time.
No CSVs, plots, or animations are generated during simulation.

## Orbital-element plots (80 and 100 orbits)

The plot command combines the 462 cases through 800 helpers at 1000 km
and targets at 850, 900, 950, 1000, 1050, 1100, or 1150 km. It reports and excludes
the stopped helper-altitude cases, whose helpers are not at 1000 km.
After all included cases are verified, run from this directory:

```bash
julia --startup-file=no --threads=1 plot_orbital_elements.jl
```

[plot_orbital_elements.jl](plot_orbital_elements.jl) activates the existing study
library environment, reads the saved data without changing it, and generates
four 4-by-7 figures: `gve_inc` and `gve_sma`, each at 80 and 100 initial target
periods. The layout follows SpaceAGORA's `plot_oe_from_feather.jl`: rows are
changes in semi-major axis, eccentricity, inclination, and RAAN; columns are
target altitudes 1150, 1100, 1050, 1000, 950, 900, and 850 km. Figure width
scales with the number of columns to preserve panel readability. Blue, red, and green curves
represent initial inclination differences 0, 0.5, and 1 degree. All eleven saved
helper populations (N=1, 50, 100, 150, 200, 250, 300, 400, 500, 600, and 800) are used.
To avoid overlapping labels on the expanded linear x-axis, ticks label N=1
and the simulated multiples of 100; markers still show every saved population.

PNG/PDF figures and numerical CSVs are saved under
[`../Output/images/orbital_elements/`](../Output/images/orbital_elements/).
Names include the scheduler and `80orbits` or `100orbits`.
`orbital_element_deltas_*orbits.csv` records all 462 target endpoints, including
the source trajectory, exact time, and one-based saved row.
`continuous_access_thresholds.csv` records the 21 analytical estimates.
Optional arguments:

```bash
julia --startup-file=no --threads=1 plot_orbital_elements.jl \
  --input ../Output/paper_plot_mode --output ../Output/images/orbital_elements \
  --orbits 80,100
```

The saved 100-orbit history contains 1001 uniformly spaced samples. Thus 80
initial target periods correspond exactly (within floating-point time precision)
to row 801, and 100 periods to row 1001. This is the 80-orbit prefix of the
100-orbit simulation, not a rescaling of final changes or a separately integrated
80-orbit run. No interpolation is used; a requested endpoint not present in the
saved grid is an error. A separately rerun integration can differ numerically.
The reader requires the full verified grid, checks trajectory SHA-256 hashes,
and uses the target index **N+1**, not `sc1` as in the reference data format.
Orbital elements use the gravitational parameter and initial state in metadata.

Changes are osculating final-minus-initial values, including both J2 and laser
effects, not laser-only changes or orbit averages. RAAN is unwrapped through
every saved row to avoid artificial 360-degree jumps in the reference script.
Equatorial RAAN is undefined and omitted from its plot row (`NaN` in CSV,
`raan_defined=false`), rather than assigning a physical meaning to zero.

### Dashed continuous-access estimates

For pooled scatter comparisons instead of the altitude-column figures, see
[Pooled orbital-element scatter plots](#pooled-orbital-element-scatter-plots).

For the initial circular geometry, the script evaluates the supplied expression
`alpha_min = acos((r_h^2 + r_t^2 - L_max^2)/(2*r_h*r_t*cos(delta_i)))`,
with radii equal to the metadata Earth radius plus altitude and `L_max = 200 km`.
Evenly spaced helpers cover every phase when `2*pi/N <= 2*alpha_min`, giving
**`N_cont = ceil(pi/alpha_min)`**. Dashed lines use the same inclination colors
as the simulation curves and mark the minimum integer population, not the
unrounded `pi/alpha_min` value (both are exported).

The simulated x-axis scale now extends through 800 helpers. The 648-helper
estimate for 850 km / 1 degree is drawn as a dashed line. The 972-helper estimate
for 1150 km / 1 degree remains annotated above the plotted range. Invalid restricted geometry
is rejected; equality or failure of the strict worst-case range condition
produces no finite population. These are fixed-geometry, range-only estimates,
not guarantees of uninterrupted laser firing in the evolving J2/laser simulation.

Targeted tests, including saved-data endpoint checks:

```bash
julia --startup-file=no --threads=1 ../../Tests/test_orbital_element_plots.jl
```

### Cumulative RTN delta-v plots

Run [plot_combined_dv_RTN.jl](plot_combined_dv_RTN.jl) for combined-scheduler
figures at 80 and 100 initial target periods:

```bash
julia --startup-file=no --threads=1 plot_combined_dv_RTN.jl
```

Outputs are PNG/PDF figures and numerical endpoint CSVs in
[`../Output/images/delta_v_RTN/`](../Output/images/delta_v_RTN/).
The three rows are cumulative radial R, along-track T, and normal N delta-v
in m/s. Seven altitude columns, all eleven helper populations, the six exact
scheduler/inclination colors, solid SMA curves, dotted INC curves, grid lines,
and continuous-access thresholds match the combined OE figures.
The same `--input`, `--output`, and `--orbits` options are supported.

Values come directly from `target_dv_R_mps`, `target_dv_T_mps`, and
`target_dv_N_mps` at saved rows 801 and 1001, after verification-status and
trajectory-checksum checks. These are signed, solver-integrated laser-only
acceleration integrals in the target's instantaneous RTN frame. They are not
absolute-value totals, Cartesian velocity differences, or a projection of the
final cumulative ECI impulse. Negative values and cancellations are retained.
Equatorial N values are valid and included (unlike undefined equatorial RAAN).
No simulation rerun or interpolation is needed.

```bash
julia --startup-file=no --threads=1 ../../Tests/test_combined_dv_RTN_plots.jl
```

### Combined scheduler altitude-column plots

Run [plot_combined_orbital_elements.jl](plot_combined_orbital_elements.jl) for
one combined 4-by-7 figure per orbit window (80 and 100):

```bash
julia --startup-file=no --threads=1 plot_combined_orbital_elements.jl
```

PNG and PDF files with `combined` in their names are saved in
[`../Output/images/orbital_elements/3_combined_plots/`](../Output/images/orbital_elements/3_combined_plots/).
Existing separate-scheduler figures are retained. The same `--input`, `--output`,
and `--orbits` options apply. Both schedulers share each panel's y-axis.
SMA curves are solid and INC curves are dotted, including their legend samples;
both retain circular data markers. Background grid lines are enabled.

| Scheduler | Initial inclination | Color |
|---|---|---|
| SMA | 0 degrees | `#652884` |
| SMA | 0.5 degrees | `#4E5689` |
| SMA | 1 degree | `#018B38` |
| INC | 0 degrees | `#B46DA9` |
| INC | 0.5 degrees | `#41B9C1` |
| INC | 1 degree | `#57AF37` |

Continuous-access thresholds depend on geometry, not scheduler. They appear
once per inclination using the darker SMA colors and apply to both schedulers.
Equatorial RAAN remains undefined and omitted. The combined renderer reuses the
original grid plotting function and verified trajectory reader.

```bash
julia --startup-file=no --threads=1 ../../Tests/test_combined_orbital_element_plots.jl
```

### Pooled orbital-element scatter plots

Run [plot_pooled_orbital_elements.jl](plot_pooled_orbital_elements.jl) to generate
four figures, each with stacked panels for changes in a, e, i, and RAAN:

```bash
julia --startup-file=no --threads=1 plot_pooled_orbital_elements.jl
```

- **OE vs number of helpers:** all seven target altitudes and three initial
  inclinations are pooled. SMA scheduler points are red; inclination scheduler
  points are blue.
- **OE vs target altitude:** all eleven helper populations and both schedulers
  are pooled. Initial inclinations 0, 0.5, and 1 degree are orange, green, and
  purple, respectively.

Both views are generated at 80 and 100 initial target periods as PNG and PDF in
[`../Output/images/orbital_elements/2_separate_plots_by_N-sceduler_h-i/`](../Output/images/orbital_elements/2_separate_plots_by_N-sceduler_h-i/).
The same `--input`, `--output`, and `--orbits` options as the grid plotter are
supported. The script reuses its verified trajectory reader and orbital-element
calculations without changing the existing figures or simulation data.

Each of the first three panels plots all 462 cases per orbit window.
The RAAN panel plots 308 non-equatorial cases; equatorial RAAN is undefined.
These are individual, partially transparent scatter points: no averaging,
jitter, connecting lines, or additional distinctions within a color group.
Coincident points can overlap. Changes include J2 and laser effects. The pooled
views omit geometry-specific continuous-access lines to avoid conflating
different geometries.

Plot grouping, coordinates, colors, and error handling are tested with:

```bash
julia --startup-file=no --threads=1 ../../Tests/test_pooled_orbital_element_plots.jl
```

The 1001 trajectory rows are about 10-11 simulated minutes apart in a full run.
Do not use them to reconstruct short encounters or accurate speed extrema.
Exact event boundaries are saved separately. During propagation, a temporary
binary spool captures position/velocity every 10 s without changing solver steps.
At finalization it supplies sampled encounter extrema and the longest-encounter
tracks, then is deleted. It is not retained as a second full trajectory.
For 301 spacecraft this temporary spool is about 0.92 GB; memory remains bounded.

### Columns and conventions

- `time`: simulated seconds from the epoch in metadata.
- `scN_pos_1:3`, `scN_vel_1:3`: J2000 ECI position (m) and velocity (m/s).
- `scN_laser_dv_1:3`: cumulative ECI laser acceleration integrals (m/s).
- `scN_laser_work_J`: cumulative optical work on each spacecraft (J).
- `target_dv_R_mps`, `target_dv_T_mps`, `target_dv_N_mps`: cumulative laser
  acceleration integrals in the target's changing RTN frame (m/s).
- `target_dv_T_positive_mps`, `target_dv_T_negative_mps`: signed along-track
  contributions; they sum to the net T value.

Spacecraft IDs equal vector indices: helpers are `1:N`, target is `N+1`.
The runner enables `SimulationConfiguration(interlink_diagnostics=true)`.
This opt-in library feature integrates RTN, signed T, and work as diagnostic
ODE states without applying another force. It defaults to false for existing
configurations; its default save fields are also available to other runners.
Only the target's RTN columns are persisted here. Direct laser SMA is not saved.

Encounters include all spacecraft pairs. Firing records include actual force
application, including the accepted out-of-range delay. `start_encounter_id` is
context at firing start, not a promise that the whole firing is inside that
encounter. The segment-level range audit provides the complete association.
Zero-force selections are not firing. Event recording is passive: it adds no
solver stops or range-cutoff callbacks.

`encounters.feather` also stores sampled minimum/maximum relative speed, range
rate, and separation, each with its time, plus the number of samples used.
Samples comprise the exact boundary states and the interior 10 s grid.
`longest_encounter.feather` retains those dense states for the longest all-pair
encounter and longest target-helper encounter. The two boolean scope columns
identify them; when both are the same encounter, the track is stored only once.
No encounters produces typed empty tables, not invented extrema.

`run_info.toml` includes exact initial Cartesian states, settings, units,
constants, source/environment checksums, completion status, and resource totals.
Selected-link history is optional; validation pilots may additionally retain
`schedule.feather`.

## Run or verify one case

```bash
julia --startup-file=no --threads=1 run_case.jl ../Output 1050 1 50 gve_sma
julia --startup-file=no --threads=1 verify_output.jl /absolute/path/to/case
```

An optional sixth argument changes the number of initial target periods for a
pilot; a seventh changes laser power. The eighth selects the helper altitude;
keep its default **1000 km** for this study. The target altitude is the second
argument. For example, a target at 1100 km with 500 helpers at 1000 km:

```bash
julia --startup-file=no --threads=1 run_case.jl ../Output 1100 1 500 gve_inc 100 10000 1000
```

Use a separate output root for shortened or otherwise nonproduction tests.

The verifier checks checksums, cadence, endpoints, IDs, finite values, event
boundaries, exact range-audit coverage, signed T sums, pair impulse balance,
J2 energy/work balance, and axial angular momentum. The independent RTN/work
quadrature and recorder-control comparisons are in
[production_pipeline_tests.jl](../../Tests/production_pipeline_tests.jl).
Signed T monotonicity checks allow 1e-7 m/s numerical differences, matching the
independent signed-component allowance; the sum-to-net check remains at 1e-8 m/s.
The saved solver values are not clamped. A reversal above the former 1e-8 m/s
monotonicity threshold emits a warning and sets `signed_T_monotonicity_warning`
in `verification.toml`. That report records the measured maximum reversal,
`signed_T_tolerance_mps` (monotonicity), and `signed_T_sum_tolerance_mps` separately.
Independent quadrature checks use 1e-8 m/s for net R/T/N and 1e-7 m/s for the
nonsmooth positive/negative T integrals.
The latter is 0.1 micrometre/s, far below the accepted step-based firing effects.
Tsit5 accepted steps, not only dense interpolation, can slightly reverse these
signed integrals at a zero crossing of along-track acceleration. An exact replay
of helper-extension case 49 reproduced a 1.7610411064e-8 m/s reversal in the
506743.3997--506750.7517 s accepted step. Independent, zero-crossing-split force
quadrature over the affected saved interval found a 2.1768755e-8 m/s signed-part
error while the net T increment agreed within 8.4e-16 m/s. This allowance accepts
bounded diagnostic integration error; it does not correct the signed integrals
or change the dynamics, solver settings, or firing schedule.
Every saved RTN/work row in validation pilots is checked in a fresh process
against the independent dense-state calculation.
The total energy/work tolerance is `max(20, spacecraft_count)` J plus
`1e-4 * abs(total_work)` to scale the absolute allowance with constellation size.
