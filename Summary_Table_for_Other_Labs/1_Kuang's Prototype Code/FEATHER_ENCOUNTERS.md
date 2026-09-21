# Feather Encounter Reports

`extract_feather_encounters.jl` converts a three-file Feather bundle into CSV
tables and encounter reports without modifying the inputs. It examines **all
unordered spacecraft pairs**, including helper-helper pairs (`sc_a < sc_b`).
New outputs use `sc1` through `scN` for the N helpers and `sc(N+1)` for the
target, matching the animation. Trajectory metadata records `target_id`;
older files without it retain the target-first interpretation.

## Running

From the workspace root, run either implementation:

```sh
julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/test16_feather.jl"
julia --startup-file=no --project=2_SpaceAGORA.jl 2_SpaceAGORA.jl/ORACLE/run_case2_laser_links.jl --feather-only
```

SpaceAGORA's `--feather-only` skips CSV generation and plots. Omit it to generate
CSVs automatically; the prototype converts automatically.

Each implementation writes:

```text
output/feather/<scenario>/trajectory.feather
output/feather/<scenario>/geometry_encounters.feather
output/feather/<scenario>/laser_on.feather
```

CSV reports stay in each implementation's output tree:

```text
1_Kuang's Prototype Code/output/CSV/<scenario>/
2_SpaceAGORA.jl/output/CSV/<scenario>/
```

Each scenario separates direct conversions from derived reports:

```text
sim_output/
  trajectory.csv
  encounter_record.csv
  laser_link_record.csv
analysis/
  encounters.csv
  encounter_states.csv
  extrema.csv
```

SpaceAGORA also writes `analysis/summary.csv`. Identical scenario names overwrite
previous outputs; archive runs before changing solver settings, which may reuse
the same folder name.

## Run Both Models

Edit `TEST_CASE` near the top of the workspace-root
[`run_comparison_v2.jl`](../run_comparison_v2.jl), then use the single-case entrypoint from the workspace root:

```sh
julia --startup-file=no run_comparison.jl
```

The script activates the SpaceAGORA environment and runs the models sequentially
in separate processes. Dependencies must already be installed. Both use shared
case settings and the same calculated duration; schedules are separately editable.
J2 is enabled, drag disabled, and SpaceAGORA beta/eta are fixed at 1.

Results are compared in [`comparison_summary.md`](../comparison_summary.md), using
spacecraft pair and chronological occurrence to match encounters. Missing values
are not zero, and the report does not establish model equivalence. A failed run
does not replace the report, so an older report may remain.

## Standalone Extraction

The extractor accepts either a scenario directory or its trajectory file:

```sh
julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/extract_feather_encounters.jl" "PATH/TO/SCENARIO"
```

Range comes from bundle metadata. Optional arguments are range in km (required
for legacy files), then the destination CSV folder. For bundles, supplied range
and scenario metadata must match the recorded events; missing interval files fail.
The default destination is `<root>/CSV/<scenario>/` for input under
`<root>/feather/<scenario>/`. The override specifies the scenario root, not
`analysis/` or `sim_output/`.

From Julia, call `extract_feather_encounters(path; output_dir=...)`; the optional
`maximum_range_m` keyword uses metres. Keep both implementation directories together:
they share the recorder and extractor.

## Trajectory Feather

`trajectory.feather` contains `time`, `sc*_pos_1/2/3`, `sc*_vel_1/2/3`, and model
diagnostics. Positions are ECI metres; velocities are m/s. Unsupported prototype
diagnostics are missing.

Samples are saved every 10 seconds by default plus the final time. Solver steps remain
adaptive; event callbacks do not add trajectory rows. SpaceAGORA's
`timeseries_points` controls buffering, not cadence.

## Geometry and Laser Feather

The CSV copies are `sim_output/encounter_record.csv` and
`sim_output/laser_link_record.csv`, respectively; Feather filenames are unchanged.

`geometry_encounters.feather` has one row per range encounter.
`laser_on.feather` has one row per uninterrupted nonzero laser interaction for a
pair, with its own unique `laser_id`. Switching off and on again creates a new
laser ID, even within the same geometry window. The columns are:

| Column | Meaning |
|---|---|
| `laser_id` | Unique laser-link activation ID within this run; laser table only |
| `encounter_id` | Unique geometry encounter ID; links laser intervals to their encounter |
| `sc_a`, `sc_b` | Unordered spacecraft IDs (`sc_a < sc_b`) |
| `start_time_s`, `end_time_s`, `duration_s` | Recorded interval bounds and their difference |
| `start_clipped`, `end_clipped` | Already active initially / still active at simulation end |
| `start_sc_a_r_x/y/z`, `start_sc_a_v_x/y/z` | A's state at interval start |
| `start_sc_b_r_x/y/z`, `start_sc_b_v_x/y/z` | B's state at interval start |
| `end_sc_a_r_x/y/z`, `end_sc_a_v_x/y/z` | A's state at interval end |
| `end_sc_b_r_x/y/z`, `end_sc_b_v_x/y/z` | B's state at interval end |

Slash notation abbreviates separate columns. Boundary states come from the
integrator at each event, not interpolation between saved trajectory samples.
IDs follow detection order and restart for each run; empty files retain their
schema. Laser rows are sorted by `laser_id`. Their `encounter_id` is only a parent
geometry reference: separate laser activations never share a `laser_id`.

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
`encounter_duration_s`, and `laser_on_s`.

| Column | Meaning |
|---|---|
| `metric` | Quantity name, including units |
| `extremum` | `min` or `max` |
| `value` | Extremal value; empty if no applicable data |
| `sc_a`, `sc_b` | Pair producing that value |
| `time_s` | Encounter sample or boundary time for speed/range rate; encounter start for durations |
| `encounter_id` | Foreign key to the encounter producing the extremum |

Speed/range-rate extrema use **only states within encounter intervals**, including
laser-link pairs. They are computed from the same states as `encounter_states.csv`:
saved samples inside each interval and recorded entry/exit states. Legacy inputs
use inferred contact windows with interpolated boundary states. Out-of-encounter
times and unrelated pairs are excluded. These are sampled, not continuous-time
extrema. Ties select the first match.
Duration extrema include clipped encounters and zero-duration touches. Filter out
either clipped flag for complete windows. Unknown laser durations are excluded;
known zeros are included.

The state CSV selected by the maximum `encounter_duration_s` row is also copied
from `analysis/encounter_states/encounter_<id>.csv` to `analysis/encounter_<id>.csv`.
The original file and encounter ID are preserved. No copy is made when there are
no encounters. The duration column in `encounters.csv` remains
`geometry_duration_s`; only the extrema metric label is renamed.

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

Laser interruptions do not split geometric encounters. Clipped durations cover
only the observed portion; boundary touches may have zero duration.

## encounter_states/encounter_<id>.csv

Both prototype and SpaceAGORA analysis exports include one file per encounter
under `analysis/encounter_states/`, named using the unchanged encounter ID.
Each file contains these 16 columns, in order:

| Columns | Meaning |
|---|---|
| `spacecraft_a`, `spacecraft_b` | Spacecraft IDs |
| `time_s` | Simulation time in seconds |
| `spacecraft_a_r_x_m`, `spacecraft_a_r_y_m`, `spacecraft_a_r_z_m` | Spacecraft A ECI position, metres |
| `spacecraft_a_v_x_m_s`, `spacecraft_a_v_y_m_s`, `spacecraft_a_v_z_m_s` | Spacecraft A ECI velocity, metres/second |
| `spacecraft_b_r_x_m`, `spacecraft_b_r_y_m`, `spacecraft_b_r_z_m` | Spacecraft B ECI position, metres |
| `spacecraft_b_v_x_m_s`, `spacecraft_b_v_y_m_s`, `spacecraft_b_v_z_m_s` | Spacecraft B ECI velocity, metres/second |
| `laser_on` | Whether any laser link in the simulation is firing at `time_s`: `true` or `false`, regardless of the pair listed in this file |

Rows retain entry/exit states and saved times inside the encounter. Regeneration
removes stale `encounter_<id>.csv` files; no encounters produces an empty folder.
The combined CSV below remains available for existing report consumers.

For bundles, `laser_on` uses all recorded firing intervals: start inclusive, end exclusive,
except a clipped final endpoint remains on. Rows are not added at firing transitions,
so brief firing intervals between saved samples may not appear as `true` rows.
Legacy single-file inputs use the held activity estimate; unavailable activity is blank.
This simulation-wide flag differs from the pair-specific `laser_on` in the combined
CSV below. Encounter laser-on durations also remain pair-specific.

## encounter_states.csv

Join on `encounter_id`. Rows contain both spacecraft at entry, exit, and saved
times strictly inside the encounter. Zero-duration encounters have one row.

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

Geometry boundaries use continuous event detection; laser intervals reflect each
model's gates and scheduler. Score-based scheduling is evaluated at accepted steps,
not continuously optimized. Grazing or brief contacts may require smaller steps
and convergence checks; event tolerance does not guarantee timestamp accuracy.

Laser totals use recorded intervals, not held trajectory snapshots.
Intervals use `[start,end)`, except a still-on final clipped endpoint is marked on
in the state CSV. Duration alone does not measure delivered energy or impulse.

## Legacy Single-File Input

Legacy `simulation_results.feather` files require a range. They estimate boundaries
from sampled ranges and laser activity from held snapshots (`laser_on_estimate_s`,
`laser_on_held`); they cannot recover missed encounters or exact event times.

No encounters yields header-only encounter/state CSVs and empty values for all extrema.
Empty cells mean unavailable or undefined; numeric zero means zero.