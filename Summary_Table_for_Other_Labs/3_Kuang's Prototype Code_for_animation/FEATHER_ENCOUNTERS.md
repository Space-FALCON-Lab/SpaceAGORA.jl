# Feather Encounter Reports

`extract_feather_encounters.jl` reads a three-file scenario bundle and writes
CSV files without changing the input files. It examines **every
unordered pair** of spacecraft, including helper-helper pairs, once (`sc_a < sc_b`).
Spacecraft IDs refer to the Feather columns: `sc1` is the target and `sc2` onward
are helpers, not the original prototype CSV ordering.

## Running

From the workspace root, run either implementation:

```sh
julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/test16_feather.jl"
julia --startup-file=no --project=2_SpaceAGORA.jl 2_SpaceAGORA.jl/ORACLE/run_case2_laser_links.jl --feather-only
```

Add `--smoke` to either command for a 60-second run. Smoke disables plots/animation
and keeps persistent artifacts under `output/smoke`, separate from production.
SpaceAGORA's `--feather-only` skips CSV generation and plots, but writes all three
Feather files and the manifest. Without it, CSV conversion runs automatically.

Each implementation writes:

```text
output/feather/<scenario>/trajectory.feather
output/feather/<scenario>/geometry_encounters.feather
output/feather/<scenario>/laser_on.feather
```

There is no `full` folder. Smoke uses `output/smoke/feather/<scenario>/` instead.
Scenario names include implementation, helper count, helper/target altitudes and
inclinations, target eccentricity/anomaly, duration, range, power, magnification,
mass, schedule, J2 setting, and the 10-second output interval. SpaceAGORA also
includes beta and eta. Identical scenario names replace the previous run's files;
old `output/test16-feather/` results are not migrated or deleted.

CSV reports stay in the implementation's own output tree:

```text
1_Kuang's Prototype Code/output/CSV/<scenario>/
1_Kuang's Prototype Code/output/smoke/CSV/<scenario>/
2_SpaceAGORA.jl/output/CSV/<scenario>/
2_SpaceAGORA.jl/output/smoke/CSV/<scenario>/
```

For a bundle under `<root>/feather/<scenario>/`, the standalone reader defaults
to `<root>/CSV/<scenario>/`. This also preserves custom output roots and smoke
separation. Legacy files outside this layout retain the prototype CSV default.

Within each scenario CSV directory, direct tables and derived reports are separate:

```text
sim_output/
  trajectory.csv
  geometry_encounters.csv
  laser_on.csv
analysis/
  encounters.csv
  encounter_states.csv
  extrema.csv
```

SpaceAGORA's additional `summary.csv` goes under `analysis/`. The `output_dir`
keyword and CLI override specify the scenario root, not either subfolder.
Returned report paths point into `analysis/`. Existing flat CSVs are not deleted
or moved automatically; use the subfolder files for newly generated results.

## Run Both Models

Edit `TEST_CASE` near the top of the workspace-root
[`run_comparison_v2.jl`](../run_comparison_v2.jl), then use the single-case entrypoint from the workspace root:

```sh
julia --startup-file=no run_comparison.jl
julia --startup-file=no run_comparison.jl --smoke
```

The script activates the existing SpaceAGORA Julia environment; dependencies must
already be installed. Both models receive a snapshot of the same case definition.
They run sequentially in separate Julia processes, with no plots or animation.
The SpaceAGORA orbit count is calculated from the requested duration in seconds,
so actual durations agree to floating-point precision. Shared settings include
helper count, orbit altitudes/inclinations, target eccentricity/anomaly, laser
range/power/magnification, and mass. Schedules remain separate editable fields;
J2 is enabled, drag disabled, and SpaceAGORA beta/eta are fixed at 1.

After both workers exit, the coordinator reads their saved analysis CSVs and writes
[`comparison_summary.md`](../comparison_summary.md) in the workspace root. It lists
conditions, summarizes each model's extrema and encounters, and compares extrema
and pairwise durations with absolute and percentage differences. Encounter matches
use spacecraft pair and chronological occurrence; absent values are unavailable,
not zero. The report explains model differences and does not claim equivalence.

`--smoke` runs at most 60 seconds, saves simulation/CSV outputs under each project's
`output/smoke`, and writes `comparison_summary_smoke.md` at the workspace root.
Root-level Markdown reports are the intentional exception to smoke-folder routing.
Reports and generated outputs for an identical scenario are overwritten on rerun.
If a worker fails, the command fails and does not generate a new combined report;
an older report may still exist. Changing only solver settings can reuse the same
scenario directory, so archive prior runs when comparing solver configurations.

The coordinator does not generate the duplicate legacy `timeseries_*.csv`.
The standalone test16 scripts retain their existing legacy-CSV behavior.

## Standalone Extraction

The extractor accepts either a scenario directory or its trajectory file:

```sh
julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/extract_feather_encounters.jl" "PATH/TO/SCENARIO"
julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/test_feather_encounters.jl"
```

The range is read from bundle metadata. An optional second argument supplies range
in km (required for old files); an optional third argument overrides the CSV folder.
A supplied range must match recorded events. Missing interval files or mismatched
scenario metadata are errors, not a silent fallback to estimates.

The prototype automatically runs this conversion and also retains its legacy wide
CSV in the scenario CSV directory. From Julia, include the extractor and call
`extract_feather_encounters(path; output_dir=...)`. The returned object contains
the summary DataFrames and paths. The optional `maximum_range_m` keyword uses metres.
Both directories must remain together in this workspace: the prototype uses the
shared recorder in `2_SpaceAGORA.jl/ORACLE/functions/FeatherRecording.jl`, and
SpaceAGORA's ORACLE runner uses the prototype's CSV reader.

## Trajectory Feather

`trajectory.feather` preserves SpaceAGORA's native non-attitude ORACLE column
names/order: `time`, `sc*_pos_1/2/3`, `sc*_vel_1/2/3`, native diagnostic fields,
and laser diagnostic snapshots. Positions are ECI metres and velocities are m/s.
The prototype leaves unsupported geodetic/atmospheric/thermal diagnostics missing.

Both runners output `0, 10, 20, ...` seconds and explicitly save the final time,
which can make the last interval shorter than 10 seconds. Solver steps remain
adaptive, with the original prototype step policy and `opts.dt_max_s` for SpaceAGORA.
Event callbacks do not add rows to this grid. SpaceAGORA's `timeseries_points`
controls output buffering, not cadence. Its default ten-orbit duration is about
63071.1894 seconds; the prototype retains 63071 seconds. Matching cadence does not
silently change either duration or make different physics/schedules equivalent.

## Geometry and Laser Feather

`geometry_encounters.feather` has one row per range encounter.
`laser_on.feather` has one row per uninterrupted nonzero laser interaction for a
pair. One encounter may contain zero, one, or several laser intervals. Both use:

| Column | Meaning |
|---|---|
| `interval_id` | Unique row ID within that file |
| `encounter_id` | Geometry ID; equals its interval ID in the geometry file |
| `sc_a`, `sc_b` | Unordered spacecraft IDs (`sc_a < sc_b`) |
| `start_time_s`, `end_time_s`, `duration_s` | Recorded interval bounds and their difference |
| `start_clipped`, `end_clipped` | Already active initially / still active at simulation end |
| `start_sc_a_r_x/y/z`, `start_sc_a_v_x/y/z` | A's state at interval start |
| `start_sc_b_r_x/y/z`, `start_sc_b_v_x/y/z` | B's state at interval start |
| `end_sc_a_r_x/y/z`, `end_sc_a_v_x/y/z` | A's state at interval end |
| `end_sc_b_r_x/y/z`, `end_sc_b_v_x/y/z` | B's state at interval end |

The slash notation abbreviates separate columns. States come from the integrator
at the event, not interpolation between Feather samples. Event IDs are assigned
in detection order. Empty event files retain the same typed schema.
Metadata records source, scenario, maximum range, timing method, units, and frame.
SpaceAGORA's manifest lists all three Feather files with sizes and SHA-256 hashes.

The reader writes `trajectory.csv`, `geometry_encounters.csv`, and `laser_on.csv`
as direct table conversions, plus the three summaries below.

## Definitions

For spacecraft A and B, at the same timestamp:

$$
\Delta\mathbf{r}=\mathbf{r}_B-\mathbf{r}_A,\qquad
\Delta\mathbf{v}=\mathbf{v}_B-\mathbf{v}_A
$$

$$
\rho=\|\Delta\mathbf{r}\|,\qquad
v_{\rm rel}=\|\Delta\mathbf{v}\|,\qquad
\dot\rho=\frac{\Delta\mathbf{r}\cdot\Delta\mathbf{v}}{\rho}
$$

- Relative speed is a nonnegative magnitude, not the difference of speed magnitudes.
- Signed range rate is **negative when approaching** and **positive when separating**.
  Swapping A and B does not change the sign. At zero separation it is undefined
  and saved as an empty CSV cell, not zero.
- A geometric encounter is a connected time interval where range is less than
  or equal to the supplied maximum range. It does not impose Earth occultation,
  minimum range, pointing, scheduling, or laser availability constraints.

## extrema.csv

Eight rows: minimum and maximum for each of `relative_speed_mps`, `range_rate_mps`,
`geometry_duration_s`, and `laser_on_s`.

| Column | Meaning |
|---|---|
| `metric` | Quantity name, including units |
| `extremum` | `min` or `max` |
| `value` | Extremal value; empty if no applicable data |
| `sc_a`, `sc_b` | Pair producing that value |
| `time_s` | Original sample time for speed/range rate; encounter start for durations |
| `encounter_id` | Foreign key for duration rows; empty for speed/range rate |

Speed and range-rate extrema cover **all pairs and all original saved times**,
including times outside laser range. They are sampled extrema, not continuous-time
optimizations. Ties select the first pair/time encountered. Duration extrema cover
all observed encounters, including clipped ones and zero-duration boundary touches.
They are not necessarily the extrema of complete physical encounters. For complete
observed windows only, filter out either clipped flag in `encounters.csv` first.
Laser duration extrema exclude unknown values but include known zero values.

## encounters.csv

One row per observed pairwise geometric encounter.

| Column | Meaning |
|---|---|
| `encounter_id` | Recorded geometry encounter ID |
| `sc_a`, `sc_b` | Feather spacecraft IDs |
| `start_time_s`, `end_time_s` | Recorded contact boundaries, in Feather time coordinates |
| `geometry_duration_s` | End minus start, in seconds |
| `laser_on_s` | Sum of durations of linked laser-on intervals, including zero if none |
| `start_clipped`, `end_clipped` | Window touches the first/last saved time; completeness is unknown |
| `maximum_range_m` | Range threshold used for this extraction |

Laser time sums all on intervals within one encounter, so interruptions do not
create new geometric encounters. A clipped duration covers only the portion in
the file. A sampled touch of the boundary can have zero duration.

## encounter_states.csv

Join to `encounters.csv` using `encounter_id` for the durations. Each encounter
contains its entry and exit times plus every original saved time strictly inside.
A zero-duration encounter has one row. Both spacecraft are recorded at each time.

| Column | Meaning |
|---|---|
| `encounter_id`, `sc_a`, `sc_b` | Encounter key and spacecraft IDs |
| `time_s` | Seconds in the source simulation's time coordinates |
| `interpolated` | False for bundle data; legacy inferred boundary rows can be true |
| `geometry_range_m` | Norm of relative position from this row's states |
| `relative_speed_mps`, `range_rate_mps` | Computed from the states in this row |
| `laser_on` | Membership in a recorded laser-on interval at this time |
| `sc_a_r_x`, `sc_a_r_y`, `sc_a_r_z` | Spacecraft A ECI position, metres |
| `sc_a_v_x`, `sc_a_v_y`, `sc_a_v_z` | Spacecraft A ECI velocity, metres/second |
| `sc_b_r_x`, `sc_b_r_y`, `sc_b_r_z` | Spacecraft B ECI position, metres |
| `sc_b_v_x`, `sc_b_v_y`, `sc_b_v_z` | Spacecraft B ECI velocity, metres/second |

## Accuracy and Laser Interpretation

Geometry uses a vector continuous callback for all pairwise `range - maximum_range`
roots, with 20 interpolation checks per accepted step and callback `abstol=1e-8`.
Coincident pair boundaries are reconciled within a 1e-6 metre range tolerance.
This callback tolerance is a condition-value tolerance, not a promised timestamp
error. The prototype also monitors minimum-range and LOS-clearance boundaries.
On a boundary, a 1-microsecond forward position/time probe selects the right-hand
laser state; timestamps and saved states stay at the detected root. The small probe
is a numerical side-selection convention, not a new output sampling interval.

The prototype recomputes its nonzero force pairs at gate events and accepted steps.
Its default `gve_schedule=:none` changes eligibility at the monitored gates; GVE
score-based switches are checked at accepted steps, not continuously root-solved.
SpaceAGORA retains its accepted-step scheduler, and records nonzero selected links
after scheduler updates as well as range crossings. These times reflect the model's
discrete scheduling resolution, not hypothetical continuous optimal switching.
Its existing impulse/force integration is unchanged; interval duration is not an
independent measurement of delivered energy or impulse. Event-induced accepted
steps can change discrete scheduler decisions compared with older runs.

Interval laser-on time is no longer reconstructed by holding 10-second trajectory
statuses. Intervals use `[start,end)`; a still-on final clipped endpoint is reported
as on in the state CSV. Trajectory `laser_active_helper` remains a convenience
snapshot and is not the authority for interval totals. Root detection remains
numerical: grazing/tangent contacts and arbitrarily brief crossings can require
smaller solver steps and further convergence checks. Speed extrema remain sampled.

## Legacy Single-File Input

Old files named `simulation_results.feather` remain readable when given a range.
Without interval files the reader uses linearly interpolated sampled ranges and
left-held helper snapshots, retaining the names `laser_on_estimate_s` and
`laser_on_held`. Such files cannot recover exact event times or missed encounters.
They are not silently upgraded to the event-recorded format.

No encounters produces header-only encounter/state CSVs and empty duration extrema.
CSV empty cells mean unavailable/undefined data, whereas numeric zero means zero.