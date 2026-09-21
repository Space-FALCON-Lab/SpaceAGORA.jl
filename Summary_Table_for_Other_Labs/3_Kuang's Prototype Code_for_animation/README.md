This version is just a cleaned version of ver3. Ver3 is already a fully functioning version.

## Feather Output

`test16_options.jl` and `test16_feather.jl` save three files per scenario:
`trajectory.feather`, `geometry_encounters.feather`, and `laser_on.feather`.
Production files go to `output/feather/<scenario>/`; all smoke artifacts go under
`output/smoke/`. CSV reports go to `output/CSV/<scenario>/` or
`output/smoke/CSV/<scenario>/`. No new `full` folder is created.
Within each CSV scenario, direct simulation tables go in `sim_output/` and derived
reports go in `analysis/`. To define and run one case in both models, edit `TEST_CASE` in the
workspace-root [run_comparison_v2.jl](../run_comparison_v2.jl) and run
`julia --startup-file=no run_comparison.jl`. See
[FEATHER_ENCOUNTERS.md](FEATHER_ENCOUNTERS.md#run-both-models) for smoke mode,
report contents, and output behavior.

## Saved Animation

### Replay Saved Feather Data

For VS Code's **Julia: Execute Active File in REPL**, open
[replay_feather_animation.jl](replay_feather_animation.jl), edit `OPTIONS` near the
top (`bundle_directory`, `duration_seconds`, and `animation_fps`), then execute
the file. It activates the shared Julia environment and starts replay using those
options, ignoring any unrelated REPL arguments. Including the whole file manually
in an interactive REPL also starts replay; noninteractive includes do not.

GLMakie requires a working graphics display. On a headless remote machine, use the
`xvfb-run` terminal command below unless the Julia REPL already has a display.

Regenerate the requested `N10`, 1000/1050 km, `T63071s`, `none` case without
running the simulation:

```sh
xvfb-run -a julia --startup-file=no --project=2_SpaceAGORA.jl "3_Kuang's Prototype Code_for_animation/replay_feather_animation.jl"
```

The default input is that case's folder under this folder's `output/feather/`.
The result is `output/videos/<scenario>/animation_replay.mp4`, with 100 seconds
of playback at 30 fps. The original animation and Feather files are untouched.
To choose another prototype bundle, playback duration, and frame rate:

```sh
xvfb-run -a julia --startup-file=no --project=2_SpaceAGORA.jl "3_Kuang's Prototype Code_for_animation/replay_feather_animation.jl" "PATH_TO_BUNDLE" 500 30
```

All three files are required: `trajectory.feather`, `geometry_encounters.feather`,
and `laser_on.feather`. The loader supports old target-first prototype bundles
and newer bundles with explicit `target_id` metadata. Helpers are displayed first
and the target last. It uses cubic Hermite interpolation of saved positions and
velocities; this is not the original solver's dense interpolation. Links use
recorded intervals, not a rerun of the scheduler. Intervals are start-inclusive
and end-exclusive, except a clipped final endpoint remains visible.

The current projection, radial scaling, and reactive legend styles are reused.
The original full-width scene and compact header layout are preserved.
The second title line shows the saved scenario's helper count, initial orbital
settings, and scheduler in normal weight; time and laser status appear below it.
Grey dashed encounter lines remain visible during laser firing, with the solid
laser drawn over them. The legend's `Encounter detected` and `Laser active`
entries can both be active; an active laser always implies a detected encounter.
For programmatic use from a noninteractive script, include the script and call
`FeatherAnimationReplay.replay_feather_animation(directory; duration_seconds=100.0, animation_fps=30, render_options...)`.
Earth radius defaults to the prototype's 6378137 m when absent from metadata.

### Five-Case Comparison Batch

From the workspace root, run both models with ten-second trajectory recording
and no animation:

```sh
julia --startup-file=no --project=2_SpaceAGORA.jl run_comparison_v2.jl
```

Each case has one target, helpers at 1000 km with zero inclination,
and a duration of 50 initial target orbital periods:

| Case | Helpers | Target altitude (km) | Target inclination | Scheduler (both models) |
|---|---|---|---|---|
| 1 | 20 | 1050 | 0 deg | gve_sma |
| 2 | 20 | 1000 | 5 deg | gve_sma |
| 4 | 100 | 1050 | 0 deg | gve_sma |
| 5 | 20 | 1000 | 5 deg | gve_inc |
| 6 | 20 | 1150 | 0 deg | gve_sma |

Each entry in `COMPARISON_CASES_V2` specifies its own `target_altitude_km`.

Requested case 3 duplicates case 2 and is run once. Other simulation settings
inherit `TEST_CASE` in `run_comparison_v2.jl`, which owns the shared settings,
model runners, and report writer. The original `run_comparison.jl` is a single-case
compatibility entrypoint that imports v2; v2 does not depend on it. Videos use the same
prototype simulation as the comparison, not a separate animation simulation.
`comparison_animation.jl` supplies the current trails, projections, radial scaling,
and legend renderer; changes to `test16_options.jl` do not configure this batch.

Output includes every ten simulation seconds plus the final endpoint. Default
scenario folders end in `dt10s`; rerunning the same scenario reuses that path. The batch
index and each comparison report record simulated duration and model run wall
time. Wall time includes setup, simulation, recording, analysis, first-call
compilation, and any requested rendering, but excludes worker startup/imports.
It is not a solver-only benchmark. The solver step cap is unchanged.

Videos are optional: append `--video` and prefix the command with `xvfb-run -a`
on a headless server. Full videos last 500 seconds at 30 fps and are saved in
`1_Kuang's Prototype Code/output/videos/<scenario>/animation.mp4`.
The batch index `comparison_summary_v2.md` links each per-case report and any requested video.
The original single-case comparison report is not overwritten.

Append `--smoke` to test all five cases with 60 simulated seconds. With `--video`,
smoke videos last 2 seconds at 2 fps. Smoke output uses each model's
`output/smoke/` directory and reports ending in `_smoke.md`.

### Single-Case Animation

This folder's `test16_options.jl` uses `OracleOptions` for the case settings.
`schedule` controls both the simulation and scenario name; `orbits` determines
duration from the initial target period. `dt_max_s` caps solver steps, independently
of the default 10-second output cadence. With the current 1050 km target, ten orbits
span about 63713.4083 seconds.

`result_plots` controls diagnostic plots. `feather_only=true` skips CSVs, plots,
and video; otherwise `animate` controls video export, with `duration_seconds`,
`animation_fps`, and `show_earth` controlling playback and rendering.

Copied SpaceAGORA compatibility fields are restricted to their supported defaults:
`planet=:earth`, `paper_grid=false`, `beta=eta=1`, and `timeseries_points=1001`.
Other values raise an error instead of being silently ignored. This prototype
does not use `timeseries_points` to set sampling. For short programmatic runs,
include the script and call `run_animation_case(OracleOptions(...))`.

Add `--video` to the comparison command to save the prototype animation as an
MP4 without opening a playback window:

```sh
julia --startup-file=no run_comparison.jl --video
```

On Linux without a graphical display, use the installed virtual display:

```sh
xvfb-run -a julia --startup-file=no run_comparison.jl --video
xvfb-run -a julia --startup-file=no run_comparison.jl --smoke --video
```

Videos are saved to `output/videos/<scenario>/animation.mp4`, or
`output/smoke/videos/<scenario>/animation.mp4` for smoke runs, and linked in the
comparison report. Playback is 100 seconds at 30 fps (2 seconds for smoke runs),
covering the entire simulated interval. The existing animation's link indicators
are geometry-based visualizations, not a replay of recorded scheduler events.
Without `--video`, comparison runs do not load graphics dependencies.

Set `VIDEO_OPTIONS` near the top of [run_comparison_v2.jl](../run_comparison_v2.jl)
to select which satellite trails appear. For only helper 1's trail:

```julia
const VIDEO_OPTIONS = (
	helper_trails = [1],
	target_trail = false,
)
```

- `helper_trails = true`: all helper trails (default).
- `helper_trails = false` or `Int[]`: no helper trails.
- `helper_trails = [1, 3]`: only helpers 1 and 3.
- `target_trail = true` or `false`: independently show or hide the target trail.

Helper IDs run from 1 through `TEST_CASE.helpers` in prototype orbit order;
helper 1 starts at zero anomaly. These are not Feather spacecraft IDs, where
`sc1` is the target and `sc2` is helper 1. Satellite markers, laser links, and
target RTN direction lines are unaffected. Rerun with `--video` after changing
the options; the scenario's existing MP4 is replaced. These settings affect
only rendering, not the saved trajectory data or simulation.

In `test16_options.jl`, setting `animate=true` also saves an MP4 instead of
displaying it; its `OracleOptions` has the same `helper_trails` and `target_trail`
settings. Direct callers of the animation function can pass these keywords plus
`output_file="path/to/animation.mp4"` and optionally `duration_seconds` and
`animation_fps`; omitting `output_file` retains interactive playback.

This folder's animation evaluates the simulation scheduler at each displayed
state: solid green/orange links indicate active lasers, and grey dashed links
indicate configured pairs that meet range and line-of-sight conditions, including
pairs with active lasers. Solid laser lines overlap the encounter dashes, with
both retaining the original line thickness.
Out-of-range or blocked pairs have no link. Scheduling and eligibility
use physical coordinates, before radial exaggeration. This is frame-by-frame
evaluation; the Feather replay script instead uses the recorded event intervals.

## Helper-Referenced Radial Scaling

In this folder's `test16_options.jl`, `radial_exaggeration = 40.0` now applies
a **radial** stretch above the helper's initial altitude, at every orbital
angle. Earth and all radii at or below that reference remain unchanged.
`1.0` disables the mapping and restores true-scale Cartesian axes.

For real altitude `altitude = radius - earth_radius`, the displayed altitude is:

```julia
altitude <= reference_altitude ? altitude :
	reference_altitude + radial_exaggeration * (altitude - reference_altitude)
```

`reference_altitude` is `reference_altitude_km * 1000`, defaulting to helper 1's
initial altitude. It must be nonnegative. The scenario runner supplies the
helper altitude, currently 1000 km. Each position vector is multiplied by
`display_radius / radius`, preserving its direction. Thus circles remain
circles, including inclined orbits, and equal altitudes have equal displayed
radii everywhere. Variable-altitude orbits are still reshaped by the radial
mapping; the reference stays fixed throughout the animation.

With the current settings, a helper at 1000 km stays at its true radius, while
a target at 1050 km appears at 3000 km altitude. The real 50 km radial gap is
drawn as 2000 km at every orbital angle, but the labels still read 1000 km and
1050 km. Unlike the previous logarithmic mapping, this expands their separation
without moving the helper orbit away from Earth. Below-reference altitudes
remain unscaled. Setting the reference to zero instead stretches all altitudes
above Earth's surface.

The outer X/Y/Z background axes carry **signed geocentric reference ticks in km**.
Internal rulers and callouts are removed. Approximately seven reference values
are selected in exact multiples of 100 km within the central 90% of each axis.
Each value is then mapped through the radial stretch to position its tick and
grid line. Labels are not rounded approximations of existing grid positions.
Grid spacing can therefore be uneven, especially across the stretch boundary.
Duplicate reference values are removed. Earth's centre and the plot centre read
**(0, 0, 0)**; opposite sides have opposite signs. Orbit geometry is unchanged.

These are not globally true Cartesian coordinate scales: the magnitude measures
distance from Earth's centre where an orbit crosses the corresponding central
axis. For example, a 7400 km reference tick represents exactly that geocentric
radius, not the nearby 1050 km-altitude orbit's radius of 7428.137 km.
A tick does not measure arbitrary off-axis satellite coordinates.
The background grid is a visual reference, not a set of constant-altitude
surfaces. With `radial_exaggeration = 1`, ordinary Cartesian km ticks return.

`axis_limit_km` is a physical geocentric radius whose mapped value sets the
symmetric scene bounds, not a real X/Y/Z coordinate limit. `nothing` fits the
scene automatically. Use the same limit, stretch, reference altitude, and figure
size to keep Earth's on-screen size consistent across runs. Too-small limits
are rejected.

`radial_ticks_km` selects actual-altitude references to include within the scene
extent, not the individual tick positions. This runner defaults to
`[500.0, 1000.0, 1010.0, 1020.0, 1030.0, 1040.0, 1050.0]`; use `nothing` for
automatic ticks, or supply another nonempty list of nonnegative altitudes.
The numerical axis labels report signed geocentric references, not altitudes or rendering distances.
Every numbered tick uses the inverse of the satellites' radial mapping. Fixed
axis limits must contain both the orbits and the requested altitude references.

Eligible cavity links are orange; single-pass links are green. Both use
linewidth 2, matching the original prototype's laser styling.
Each laser legend key follows all eligible links of its type, not just the
first configured pair. It is orange (cavity) or green (single-pass) while
active, and grey with `OFF` when none is active. Active labels list satellite
IDs, for example `Open-cavity pair: satellites 10-11`; simultaneous pairs are
comma-separated. The laser legend keys are status indicators, not visibility
toggles for individual beams.
The video header shows ON/OFF based on physical LOS/range and the animation's
pair-selection logic, not recorded scheduler events. No beam is drawn for an
ineligible pair. For the 1000/1050 km, 200 km-range, 63071-second case, saved
laser intervals occupy approximately the first 4.1 and last 5.5 seconds of
the 100-second playback, so no beam is expected through most of that video.

Physics, laser eligibility, reported distances, and CSV/Feather data remain unscaled.
Markers, trails, link endpoints, and RTN arrow origins use display positions;
RTN directions use the real states. Existing video files are replaced on rerun.
Straight links are schematic connections between displayed endpoints, and RTN
arrows indicate physical directions; neither is a distance ruler in this view.

Run this copy (not the root comparison runner):

```sh
xvfb-run -a julia --startup-file=no --project=2_SpaceAGORA.jl "3_Kuang's Prototype Code_for_animation/test16_options.jl"
```

The current user scenario remains unchanged; set the altitude and scale values
for your intended case before running. For example, do not apply 40x to a
1000/5000 km case expecting the same framing as the 1000/1050 km example.

Run the focused rendering tests with:

```sh
xvfb-run -a julia --startup-file=no --project=2_SpaceAGORA.jl "3_Kuang's Prototype Code_for_animation/test_animation.jl"
```

These tests also generate `output/smoke/videos/helper_referenced_orbit_test.mp4` and
`.png`, showing a synthetic inclined circular orbit over a full revolution.
This preview is a geometric regression fixture, not the full scenario run.

## Three-Plane Projections

`show_projections = true` in this folder's `test16_options.jl` enables dashed
trail projections and smaller satellite position markers on the three plot
boundary planes: XY on the bottom face, XZ and YZ on the two positive side
faces. They are offset just inside the boundaries to avoid depth flicker.
Set it to `false` to hide all projections. Direct animation calls default to
`false` for compatibility.

`helper_projections` and `target_projections` independently enable each group's
projected trails and markers. The runner is configured for target-only projections:

```julia
show_projections::Bool   = true
helper_projections::Bool = false
target_projections::Bool = true
```

Set both group switches to `true` for both groups; `show_projections=false`
overrides both. Direct animation calls default both group switches to `true`
for compatibility. Original 3D satellite markers are unaffected.
Within enabled groups, projected trails still obey `helper_trails` and
`target_trail`; projected position markers remain visible even when trails
are hidden. Rerun the scenario to update an existing video.
Projections use the same time, trail length, colours, and exaggerated display
coordinates as the 3D scene. They are orthographic copies on offset planes,
not additional orbits or physical positions. Earth, lasers, and the radial
altitude axis are not projected, and physics and exported data are unchanged.
For equatorial orbits, XZ and YZ projections naturally collapse to lines.

## Export Dependencies

Use the neighboring SpaceAGORA project environment for Arrow, CSV, DataFrames,
and solver dependencies. The projects share the event recorder in SpaceAGORA's
ORACLE functions folder and the CSV extractor in this prototype folder.

`functions/13_Feather_Write.jl` provides `save_timeseries_feather(sol, p; feather_dir)`.
Include it after the existing physics and diagnostic helpers to use it in other
prototype runners. It uses the same Arrow IPC file encoding, column names, and
column order as SpaceAGORA's native non-attitude ORACLE output:

- `time` is in seconds; `sc1` is the target, and `sc2` onward are helpers.
	This differs from the prototype CSV, where helpers come first and the target last.
- `sc*_pos_1/2/3` and `sc*_vel_1/2/3` hold ECI states in metres and metres/second.
- Mass, periapsis altitude, target `dv_r/t/n_accumulated`, and
	`laser_active_helper` are populated. Helper IDs use the output spacecraft order;
	zero means no active helper. Delta-v uses the existing prototype's trapezoidal
	RTN integration, not SpaceAGORA's callback accumulator.
- Altitude, latitude, longitude, wind, drag, lift, cross, heat rate, and heat load
	are nullable Float64 columns containing `missing`: the prototype does not record
	SpaceAGORA's geodetic/atmospheric/thermal diagnostics. They are not zero-valued
	measurements. The Arrow schema metadata records these unavailable fields and
	the original satellite IDs.
- Trajectory samples are saved every 10 seconds, plus the final endpoint. The
	adaptive solver chooses its steps independently of output cadence. No attitude columns or
	SpaceAGORA bundle manifest are generated by the prototype. The exporter
	supports one target, matching the ORACLE case. Enable `record_events=true`
	in `run_open_cavity_multi` to generate the interval files.
- Geometry entry/exit and laser state changes are recorded during integration,
	independently of trajectory output. See [FEATHER_ENCOUNTERS.md](FEATHER_ENCOUNTERS.md)
	for interval schemas, root/accepted-step timing semantics, and accuracy limits.

## Test16 Feather Case

Run from the workspace root:

```sh
julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/test16_feather.jl" --smoke
julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/test16_feather.jl"
julia --startup-file=no --project=2_SpaceAGORA.jl "1_Kuang's Prototype Code/test_feather_encounters.jl"
```

The headless case is derived from `test16_options.jl`, retaining ten helpers,
1050/1000 km helper/target altitudes, the original laser parameters, J2 setting,
and `gve_schedule=:none`. The default duration is the original 63,071 seconds;
`--smoke` runs 60 seconds. Neither plots nor animation are loaded, and the runner's
image cleanup is confined to a temporary directory.

Each run uses a scenario folder identifying the implementation, orbit and laser
settings, duration, and output cadence, then verifies column order, types, every spacecraft's
state mapping, mass, periapsis altitude, delta-v, helper IDs, and CSV round-trip.
These tests validate export compatibility, not numerical equivalence between the
prototype and SpaceAGORA's different physics/scheduling implementations.
Old output folders are left untouched. Delete `output/smoke` when its results are
no longer needed; production runs do not write there.