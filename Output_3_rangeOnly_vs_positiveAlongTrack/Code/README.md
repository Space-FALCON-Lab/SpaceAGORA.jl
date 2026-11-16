# Production range-only / positive-along-track comparison

Additional [last-two-orbit averaged figures](../Output/images/last_two_orbits_average/)
are available for the 80- and 100-orbit windows. See the shared
[averaging method and reproduction guide](../../LAST_TWO_ORBITS_AVERAGES.md).
The original endpoint figures and source simulations are unchanged.

This campaign reuses the [Output 2 production runner](../../Output_2_Simulation_absolute_RN/Code/StudyRunner.jl),
not the older 80-orbit pilot integrator. Existing pilot and Output 2 results are
not overwritten. The original pilot remains under `h1050km_t1000km`; production
uses `h1000km_t1000km`, `h1000km_t1050km`, and `h1000km_t1100km`.

## Confirmed grid and numerical settings

- Helpers: 1000 km, inclination 0 deg, circular and evenly phased.
- Target: 1000/1050/1100 km, inclination 0/0.5/1 deg, circular, anomaly 1 deg.
- Helper counts: 1, 100, 200, 300, 400, 500, 600, 700, 800.
- Policies: `naive_next_entering` (range-only), `positive_along_track`.
- **162 laser-on cases**, plus **81 matched no-laser cases** in
  [Output 4](../../Output_4_no_laser/Code/README.md). One baseline per geometry
  and helper count is shared by the two laser-on policies.
- Each case spans **100 initial target periods**. J2, 227 kg spacecraft,
  10 kW lasers, magnification 100, 200 km selection range, one terminal and
  target-helper-only links match Output 2. Baseline power is zero.
- Tsit5, relative/absolute tolerances `1e-12`, maximum accepted step **10 s**.
  All helpers receive gravity and the active helper receives the equal-and-opposite
  laser force.
- Scheduling at initialization and after accepted steps. Selected links remain
  applied between updates, including out-of-range firing. This deliberately
  matches Output 2, not the old pilot's 1 s/continuous-range-cutoff behavior.
- **1001 uniform trajectory rows**; rows 801 and 1001 give the 80- and 100-orbit
  endpoints. Encounter products use the existing 10 s spool and exact boundaries.

Range-only retains its helper while in range; replacement considers newly
entering helpers, not every already-available helper. It can intentionally be
idle while another helper is already in range. Positive-along-track retains its
helper while in range and its force direction has positive projection onto the
target's instantaneous RTN transverse axis. Replacement maximizes that positive
directional projection, not a GVE rate. Deterministic ties use helper index.
The sign gate is at scheduling updates, not continuously between updates.

Both runners share the [single-target policy implementation](../../9_SpaceAGORA.jl-main_full%20integration/src/models/satellite_interlink/4_single_target_policies.jl)
through [scheduling_policies.jl](scheduling_policies.jl); the original pilot's
default state-machine behavior and range-gated pair-force function are retained.

## SSH-independent full job

Run from the project root:

```bash
bash Output_3_rangeOnly_vs_positiveAlongTrack/Code/launch_study.sh start
bash Output_3_rangeOnly_vs_positiveAlongTrack/Code/launch_study.sh status
```

Manual `start` uses tmux. The `run` action is the foreground entry point for
an external detached-process launcher; it redirects stdin from `/dev/null`,
logs stdout/stderr to disk, and holds an exclusive `flock`. A detached launch
survives an ordinary SSH disconnect, but not a reboot or machine shutdown.

[run_jobs.jl](run_jobs.jl) automatically runs:

1. All Output 3 simulations, each followed by fresh-process output verification.
2. Output 3 paper figures, endpoint CSVs, and encounter evaluation.
3. All Output 4 no-laser simulations and verification, then their figures/evaluation.
4. [Laser-minus-no-laser data and figures](../../Output_5_Enable_Disable_Laser_compare/Code/README.md).

Only one simulation runs at once, with one Julia/BLAS/GC thread. Disk guards
reserve 5 GiB plus the next case's temporary spool and output estimate. The
completed-case spool is removed. Source hashes are checked between cases and
phases: do not edit simulation, launcher, or plotting code while the job runs.
Failures are explicit, stop the job, and preserve completed cases and diagnostics.

Status/logs are under `Output/_comparison_job/`:

- `status.toml`: overall phase, PID, output paths and completed campaign counts.
- `job.log`: controller and case-level completion messages.
- `plotting.log`: figure/evaluation progress and errors.
- `../_batch_range_along/status.toml`: current case and child PID.
- `../_batch_range_along/logs/`: individual simulation and verification logs.
- Each active case has `progress.toml` with simulation time and fraction.

Pause/resume:

```bash
bash Output_3_rangeOnly_vs_positiveAlongTrack/Code/launch_study.sh stop-after-case
bash Output_3_rangeOnly_vs_positiveAlongTrack/Code/launch_study.sh resume
```

Pause waits for the current case and its verifier; it does not kill a simulation.
Completed matching cases are reverified on resume. Incomplete/failed cases are
never overwritten automatically: inspect their logs before moving that specific
case aside and retrying. If the source-safe resume guard rejects older code,
do not bypass it or mix those results into this comparison.

## Saved products

Production uses the same hierarchy and schema 3 as Output 2:

```text
Output/paper_plot_mode/h1000km_t1050km/ih0.0deg_it0.5deg/N100/T...s/
  naive_next_entering_e0.0000_nu1.0000/
  positive_along_track_e0.0000_nu1.0000/
```

Each leaf contains `trajectory.feather`, `encounters.feather`,
`laser_firing.feather`, `firing_range_audit.feather`, `longest_encounter.feather`,
and `run_info.toml`, `verification.toml`, `progress.toml`. The target is spacecraft
`N+1`. Full spacecraft position/velocity, ECI laser impulse per unit mass,
optical work, and signed/positive/negative R/T/N laser integrals are retained.
Metadata includes settings, source hashes, checksums and solver/resource statistics.

Generated summaries:

- `Output/production_case_summary.csv`: case provenance, events and firing totals.
- `Output/production_endpoint_summary.csv`: complete 80/100-orbit endpoint data.
- `Output/encounter_evaluation/`: per-event CSVs, distributions and a generated
  report using the same definitions as Output 2, including no-firing cases.
- `Output/images/for_paper/`: 100-orbit OE and split-N RTN figures, each PNG/PDF,
  full-precision source CSVs, access thresholds and a figure manifest.
- `Output/images/80orbits/`: corresponding 80-orbit figures and tables.

The paper renderer retains the original four-row, three-altitude-column layout,
colors, typography and scales, with correct new scheduler legends and data-fitted
limits (no clipping). Equatorial RAAN stays undefined and is omitted.

## Reproduction and tests

```bash
# Run only Output 3; this command is foreground.
julia --startup-file=no --threads=1 --gcthreads=1 Output_3_rangeOnly_vs_positiveAlongTrack/Code/run_grid.jl
# Regenerate its plots/evaluation without propagation.
julia --startup-file=no --threads=1 --gcthreads=1 Output_3_rangeOnly_vs_positiveAlongTrack/Code/plot_campaign.jl
# Focused fixtures, full-duration small cases, high-helper pilot and renderer checks.
julia --startup-file=no --threads=1 --gcthreads=1 Output_3_rangeOnly_vs_positiveAlongTrack/Code/test_production_campaign.jl
```

The tests check real saved outputs, exact 80/100-orbit sampling, zero-laser
diagnostics, mismatch rejection, common-frame velocity subtraction, RAAN
unwrapping, every plotted curve, and unchanged saved Output 2 GVE behavior.
Temporary test cases/figures are cleaned up; the validation result is retained in
`Tests/Output_3_rangeOnly_vs_positiveAlongTrack/validation/production_campaign/validation.toml`.

## Firing time versus target argument of latitude

[plot_firing_latitude.jl](./plot_firing_latitude.jl) compares the **300- and
800-helper range-only runs**, target altitude **1100 km**, initial inclination
**1 degree**, over all **100 initial target reference periods**. It writes only
to [this diagnostic folder](../Output/images/firing_latitude_h1100_i1_range_only/);
it does not change trajectories, schedules, or the existing endpoint plots.

- [Overlaid firing-time distribution](../Output/images/firing_latitude_h1100_i1_range_only/laser_on_time_vs_argument_of_latitude.png):
  hours of laser-on time in each 10-degree bin of instantaneous target argument
  of latitude, not an event-start count or geographic latitude.
- [Inclination contributions and effect per firing second](../Output/images/firing_latitude_h1100_i1_range_only/inclination_contribution_vs_argument_of_latitude.png).
- [Combined firing time and inclination contribution](../Output/images/firing_latitude_h1100_i1_range_only/laser_on_time_and_inclination_vs_argument_of_latitude.png):
  solid steps use the left axis (laser-on hours); dashed lines with circles use
  the right axis (laser inclination change in `1e-3` degrees), both per 10-degree
  bin. Purple remains 300 helpers and gold 800 helpers. The independent axes
  compare orbital phase, not numerical magnitudes or equality at curve crossings.
  Node labels are retained, and the original separate figures remain available.
  Both combined versions place a black-outlined, one-row, four-column legend below the axes,
  identifying each helper count's solid on-time and dashed inclination curves.
  Legend samples use the same rendered line widths, dash style, and circle size
  as the plotted curves, rather than the backend's enlarged automatic symbols.
- [Clean combined figure](../Output/images/firing_latitude_h1100_i1_range_only/laser_on_time_and_inclination_vs_argument_of_latitude_clean.png):
  the same dual-axis curves without a title, node labels, node reference lines,
  or helper-plane approximation note. Numeric ticks, axis labels, helper legend,
  ordinary grid lines, and the dotted zero-inclination reference remain.
- [Target latitude versus selected-helper relative phase](../Output/images/firing_latitude_h1100_i1_range_only/joint_target_latitude_helper_relative_phase.png):
  hours per 10-degree by 0.1-degree cell; both panels use the same color scale.

Each figure also has a PDF. The accompanying
[latitude bins](../Output/images/firing_latitude_h1100_i1_range_only/latitude_bins.csv),
[joint bins](../Output/images/firing_latitude_h1100_i1_range_only/joint_latitude_relative_phase.csv),
and [case summary](../Output/images/firing_latitude_h1100_i1_range_only/case_summary.csv)
retain unrounded values. Latitude-bin tables also include mean link distance,
helper signed distance from the target plane, normal acceleration, and relative
phase.

### Method and checks

The sparse trajectory samples cannot resolve all firing/latitude boundaries.
Instead, each recorded firing interval is reconstructed from its saved helper
and target boundary states using shared central-plus-J2 gravity and the recorded
equal/opposite **BP/c** force. The selected link stays active throughout that
recorded interval, including the short out-of-range portions held by the
production scheduler. No new scheduling decisions are made.

Dense ODE output locates angular-bin crossings; the code splits intervals there
and accumulates their actual durations. The simultaneous inclination integral is

```text
di_laser/dt = (r / |r x v|) cos(u) a_N
a_N = -(F/m) (r_helper . N_target) / |r_target - r_helper|
```

Here `u` is measured in the instantaneous target orbital plane from its ascending
node. The helper-relative phase is
`atan2(r_helper . T_target, r_helper . R_target)`: the helper direction projected
into that same plane, positive ahead of the target and negative behind. It is
**not** the difference between the spacecrafts' separate arguments of latitude.
The helper's signed cross-track position must also be considered.

Except for the clean combined figure, all panels mark the target's **ascending node relative to the helper
plane near 0 degrees (repeated at 360)** and **descending node near 180 degrees**
with dashed vertical lines and labeled ticks. These are intersections of the
two **orbital planes**, not exact spacecraft encounters. The ascending mutual
node points along `N_helper x N_target`; the target crosses toward the positive
helper-plane normal there. Its opposite is the descending node.

Helpers start equatorial but can tilt slightly under laser recoil, so the
labels are explicitly approximate. The actual node is computed from both
spacecrafts' instantaneous angular momenta at every firing-replay solver state.
The plotter requires its signed offset from 0 degrees to stay within
0.05 degrees, checks sampled-extrema convergence between 10- and 5-second
replays, and saves each case's minimum/maximum in the case summary.
The descending-node offsets are the same, added to 180 degrees.

Vern9 replays with 10-second and 5-second step caps, both at `1e-12` absolute and
relative tolerances, check convergence; plots use the 5-second results.
Maximum latitude-bin on-time disagreement was below `1.1e-8` seconds, and maximum
per-bin laser inclination disagreement was below `4.1e-16` degrees.
Replayed endpoints agreed within `1.02e-5` metres and `1.04e-8` m/s.
Integrated normal impulse is checked against the recorded diagnostic.
Laser plus J2 inclination contributions close the saved endpoint change within
`1e-14` degrees. J2 during idle gaps is obtained from their saved boundary
inclination changes. Time conservation, convergence, and source/output checksums
are recorded in [validation.toml](../Output/images/firing_latitude_h1100_i1_range_only/validation.toml).

### What explains the dip

| Quantity over 100 reference orbits | 300 helpers | 800 helpers |
|---|---:|---:|
| Laser-on time, hours | 122.2012 | 164.5891 |
| Duty fraction | 0.68356 | 0.92067 |
| Positive laser inclination contribution, degrees | +0.003641368 | +0.008611429 |
| Negative laser inclination contribution, degrees | -0.009231693 | -0.009202873 |
| Net laser inclination change, degrees | -0.005590325 | -0.000591444 |
| J2 inclination contribution, degrees | -0.000639590 | -0.000738397 |
| Total inclination change, degrees | -0.006229915 | -0.001329841 |

The 300-helper case is missing much of the **inclination-raising exposure**:
about **40.85 versus 81.35 hours** in `u = 0-90` and `180-270` degrees combined.
In the predominantly lowering quadrants (`90-180` and `270-360` degrees), exposure
is much closer: **81.35 versus 83.24 hours**. Quadrants are approximate sign
regions; contributions in the table are classified by the actual instantaneous
Gauss rate, including small sign changes near the nodes.

The per-firing-second inclination curves are fairly similar, so the
nonuniform exposure explains most of this comparison. Negative contributions
are nearly equal; the 300-helper run has far less positive cancellation.
Exposure alone is not a general predictor: link direction, helper cross-track
position, distance, the `cos(u)` factor, and J2 all matter. The joint plot shows
different selected-helper phase patterns, which affect both access timing and
force geometry. This diagnoses the recorded schedules, not a universal
300-helper resonance or optimality of either schedule.

### Reproduction

Run from the repository root with the existing Julia environment:

```bash
OPENBLAS_NUM_THREADS=1 julia --startup-file=no --threads=1 --gcthreads=1 \
  --project="9_SpaceAGORA.jl-main_full integration" \
  Output_3_rangeOnly_vs_positiveAlongTrack/Code/plot_firing_latitude.jl
OPENBLAS_NUM_THREADS=1 julia --startup-file=no --threads=1 --gcthreads=1 \
  --project="9_SpaceAGORA.jl-main_full integration" \
  Output_3_rangeOnly_vs_positiveAlongTrack/Code/test_firing_latitude.jl --saved
```

The second command tests angle/frame definitions, wraparound, exact synthetic
bin crossings, duration conservation, mutual-node orientation and crossing signs,
and Gauss-force signs, then checks the
generated CSV shapes, bin marginals, numerical identities, and checksums.
Omit `--saved` to run only the synthetic tests without generated files.
