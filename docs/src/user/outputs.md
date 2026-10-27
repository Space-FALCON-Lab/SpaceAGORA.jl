# Simulation Outputs

This page describes the files `run_simulation` writes, how to read them back,
and what checkpointing and resuming do to them. It assumes you have already run
a simulation. The quickstart example is the shortest way to get one:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

Related pages:

- [Output Columns](output_columns.md)
- [The Integrated State](integrated_state.md) (which columns are integrated and which are derived)
- [Simulation Configuration](simulation_configuration.md)
- [Verification Study](verification_study.md)
- [Recipes](recipes.md)

## Output files

With the default settings, `run_simulation` writes three files to
`simulation_settings.results_directory`:

```text
output/
  simulation_results.csv         ← tabular results, one row per saved time
  simulation_results.feather     ← the same table in Apache Arrow (Feather) format
  simulation_results.manifest.toml ← metadata: row count, spacecraft count, file hashes
```

The manifest records `schema_version` (`"1"` for this layout), `created_utc`,
`mission_time_s`, `steps` (the number of rows), `spacecraft_count`,
`orientation_sim`, `first_time_s` (the first row's time) and an `[atmosphere]`
table that identifies the density model. Under `[files]` it lists each data
file with its path, size and SHA-256 hash. A resumed run's manifest also
records `resumed_from_s`, the checkpoint time it resumed from (see
[Checkpointing](#Checkpointing)). `created_utc` is the wall-clock time the file
was written, not the simulation epoch. The manifest does not store the epoch,
which is the `InitialTime` you configured, so keep it with your script or read
it from the returned configuration.

A run with interlinks also writes its link schedule,
`simulation_results_link_schedule.csv` and
`simulation_results_link_schedule.feather` (see
[the link schedule's columns](output_columns.md#Link-schedule)). The manifest
lists both under `[files]` and describes the schedule in an `[interlinks]`
table:

- `policy` and `score_unit`: the scheduling policy's type and the unit of its
  score;
- `rows`: the number of rows in the schedule;
- `record_available` and `record_terms`: whether the history recorded the
  available links (the contact rows) and a weighted policy's term scores;
- `history_start_s` and `history_end_s`: the times of the interlink history's
  first and last samples, the span the schedule covers (absent when the
  history has no samples);
- `terms`, for a `WeightedPolicy`: each top-level term's `policy`,
  `score_unit`, `weight` and `scale`, in the order of the `term<k>_start`
  columns.

A `run_mission` run with operational segments writes its segment history to
`simulation_results_mission.json`. When it also saves its scene, a
`[mission_history]` entry in the manifest ties the history to these files (see
[Define, run and inspect a mission](mission_workflows.md)). Writing new results
under the same prefix removes that entry, so a saved mission history is never
attached to another run's trajectory.

Set `simulation_settings.save_csv = false` to write the Feather file and the
manifest without the CSV. A CSV left by an earlier run stays in the folder, but
it is not listed in the new manifest or in `SimulationResults.files`. Set
`simulation_settings.results = false` to write no result files at all.

The developer switch `SPACEAGORA_SAVE_BUNDLE=0` turns off the Feather files and
the manifest. With it and `save_csv=false` no results file is written, so pass
`return_results=true` to get the table in memory. A scene sidecar or
checkpoints that you enable explicitly are still written.

## Getting the results in memory

Pass `return_results=true` to get a [`SimulationResults`](@ref) back. Its
`table` is a `DataFrame` with the same columns and values as the CSV. The table
is built even when `simulation_settings.results = false`, in which case no
result files are written. `configuration` is the configuration that actually
ran. Under the default `isolate_state=true` that is a deep copy, so controller
state and logs, such as the RPO command log, are on it and not on the object
you passed in. `files` lists the output files this call wrote, and `solution`
holds the ODE solution when you also pass `return_solution=true` (otherwise it
is `nothing`). `return_results=true` cannot be combined with
`return_solver_metadata=true`.

Checkpoints you enable are still written, and `files` lists them only when
this call wrote them. Resuming from a checkpoint at the mission's end solves
nothing: it writes no result files and does not report the old checkpoint
files as new output. Its `table` holds the whole run, rebuilt from the rows
kept beside the checkpoint, or is empty when the checkpoint holds no usable
rows.

For a configuration `args` you have already built:

```julia
res = run_simulation(args; return_results=true)
res.table                      # DataFrame, one row per saved step
res.files                      # e.g. ["output/simulation_results.csv", ...]
res.configuration              # the (isolated) configuration that ran

# Read back a setting that actually ran, for example spacecraft 1's inclination.
# Keplerian InitialCondition fields are stored in SI units and radians.
ic = res.configuration.dynamics_model.spacecraft[1].initial_condition
rad2deg(ic.i)                  # inclination in degrees

res = run_simulation(args; return_results=true, return_solution=true)
res.solution.t
```

## Where the examples write, and how to keep runs apart

Most repository examples build their configuration with `make_example_config`.
Its `results_directory` is `<repository root>/output`, or the directory named
by `SPACEAGORA_CLI_OUTPUT_DIR` when that variable is set; the CLI's
`--output-dir` option sets the variable. The result file names never change,
so two runs that share a results directory overwrite each other's files. To
give a run its own directory, start it through the CLI from the repository
root:

```text
julia --project=. src/interfaces/cli/main.jl run --example=AGORA_Earth_NoGRAM.jl --output-dir=output/earth-nogram
```

`--output-dir` works for the scripts that take their directory from
`make_example_config`: the first-run scripts (`AGORA_Basic_Quickstart.jl`,
`AGORA_Earth_NoGRAM.jl`), the control and torque tests, and the mission
scripts `AGORA_Basic_GRAMEarth.jl`, `AGORA_Odyssey.jl`, `AGORA_Vex.jl`,
`AGORA_Titan.jl`, `AGORA_Magellan.jl`, `AGORA_LOFTID.jl`,
`AGORA_Mars_NoGRAM.jl` and `Earth_Thruster_Test.jl`. It has no effect on the
scripts that choose their own directory. `AGORA_Earth.jl`,
`AGORA_Keplerian.jl` and `CYGNSS_test.jl` write to `output/`,
`Earth_Navigation.jl` to `output/observer_noDA/`, `AGORA_Earth_Aerobraking.jl`
to `output/earth_aerobraking/`, and the RPO, robot-arm, compliant multibody and
interlink demos to folders of their own under `output/`. Running one of these
twice overwrites its previous results. In your own scripts, pass
`results_directory=` to `make_example_config` or set it on
`SimulationSettings`.

`odyssey_surrogate.jl` runs in its own environment and manages its own files.
The comparison writes to `odyssey_surrogate_results/` in the current directory
unless you pass `--output=DIR`, refuses a directory that already exists, and
prints the absolute path when it finishes (see the
[Odyssey walkthrough](../tutorials/odyssey_surrogate.md)). The Mars RAAN
campaign lives in the separate SpaceAGORAMissionReconstructions repository,
which documents its launcher and output paths.

Some examples write more than the result files. The quickstart saves four PNG
plots under `<results_directory>/plots/`, and the RPO examples also write HTML
plots. `Solar_Panel_Cloth_Deployment_Demo.jl` writes four HTML files to its
demo folder instead of the three result files. `AGORA_Earth_MonteCarlo.jl`
turns result saving off and writes no result files, even with `--output-dir`;
it prints the numbers of successful and failed samples and the elapsed time to
the terminal.

Smoke mode (`--smoke` on the CLI, or `SPACEAGORA_EXAMPLE_SMOKE=1` for a script
run directly) shortens the mission, normally to at most 120 s and one orbit.
Most examples apply the shared smoke settings, which also turn off plots, the
scene sidecar and checkpoints, and keep results only when
`SPACEAGORA_EXAMPLE_SMOKE_RESULTS=1`, which the CLI sets for you. Those results
go to the `--output-dir` directory when one is given, and otherwise to
`output/` under the current working directory, where they replace any earlier
results. Choose a separate directory to keep earlier results.

### Visualization scene sidecar (opt-in)

Setting `simulation_settings.save_visualization_scene = true` adds a scene
sidecar and some result columns, which the post-hoc 3D viewer described in
[Interactive Visualization](visualization.md) reads:

```text
output/
  simulation_results_scene.json  ← planet, epoch, spacecraft geometry, frame rotation samples
```

The sidecar holds what the viewer needs besides the trajectory rows: the
epoch; the central body's radii, spin and texture key; a table of
J2000-to-body-fixed quaternions sampled over the mission span, so the viewer
needs no SPICE kernels; the atmosphere's entry interface and, for models that
can be sampled, a density profile; and every spacecraft's links as boxes built
from `Link.dims`, with thruster, facet and joint glyphs. A run with interlinks
also names its link schedule, and a run with ground stations records their
positions. The run writes the sidecar after its results, and only when both the
flag and `results` are on; the other output files are the same as without the
flag.

With the flag on, the run also saves the following fields when it can supply
them. [Output Columns](output_columns.md#Visualization-fields) lists their
columns and units.

- `link_pose`, for spacecraft with non-root links: `sc{N}_link_pose_1` to
  `sc{N}_link_pose_7n`, holding `[rx, ry, rz, qx, qy, qz, qw]` for each
  non-root link relative to the root bus, in the link order the sidecar
  records. Articulated panels then play back as they moved.
- `thruster_direction`, when a spacecraft carries a thruster gimbal:
  `sc{N}_thruster_direction_1` to `sc{N}_thruster_direction_3k`, every
  thruster's exhaust direction in its link's frame, so the viewer turns the
  gimballed thruster's cone and plume as the gimbal moved.
- `thruster_level`: `sc{N}_thruster_level_{k}`, the firing level of thruster
  `k` from 0 to 1, for spacecraft with thrusters whose control effector reports
  per-thruster levels or runs in `:mounted` actuation. Thrusters are numbered
  in the sidecar's order: the links in order, and each link's thrusters in
  order.
- `sun_dir`: `sun_dir_1` to `sun_dir_3`, the unit vector from the planet's
  center to the Sun in the inertial frame of the saved positions. There is one
  direction per row, shared by every spacecraft. It is saved when the run's
  ephemerides can locate the Sun at the start epoch, and the viewer's sun
  lighting reads it.
- `arm_pose`, for a spacecraft with a compliant robot-arm plan:
  `sc{N}_arm_pose_1` to `sc{N}_arm_pose_7m`, holding
  `[rx, ry, rz, qx, qy, qz, qw]` for each arm link: the link's center of mass
  relative to the spacecraft position in inertial meters, and the link's
  inertial attitude quaternion.
- `density`, for a run with an atmosphere model: `sc{N}_density`, the
  atmospheric density along the trajectory in kg/m³.

An explicit `save_fields` list receives the visualization fields it lacks,
appended once after the fields you named. To read the sidecar back:

```julia
using SpaceAGORA
scene = read_visualization_scene("output/simulation_results_scene.json")
scene.planet.name, length(scene.spacecraft), scene.spacecraft[1].bounding_radius_m
```

Two kinds of columns are often viewed with the sidecar but do not depend on
it. Articulated spacecraft and spacecraft with compliant attachments save
their joint states and body poses in default fields, listed under
[Articulated spacecraft and attachments](output_columns.md#Articulated-spacecraft-and-attachments);
the `link_pose` field still holds the configured rigid geometry. The plume
columns described under
[Plume interaction](output_columns.md#Plume-interaction) belong to the default
field list whenever the plume-surface effector is in the run. An explicit
list that omits them does not write them.

### Interactive viewer page

HTML pages come from the optional SpaceAGORAVisualizer package: load it with
`using SpaceAGORAVisualizer` before asking for a page (see
[Install the optional viewer](visualization.md#Install-the-optional-viewer)).
With the sidecar present, `export_visualization` builds a self-contained HTML
page next to the results. `run_simulation(args; visualization=true)` turns the
sidecar on and builds the page in the same call:

```julia
run_simulation(args; visualization=true)
# or, from an existing run that saved its sidecar
export_visualization("output/simulation_results"; max_frames=2000, frame=:inertial)
```

`SPACEAGORA_VISUALIZATION=1` makes `visualization=true` the default for every
`run_simulation` call in the process, so a script can opt in without being
edited. The process still needs SpaceAGORAVisualizer loaded; without it the
run stops with an error before it starts. The CLI's
`run --example=<script> --visualize` sets the variable and loads the viewer
from a prepared visualizer environment, as shown in
[Make your first page](visualization.md#Make-your-first-page). Without the
viewer package, `export_viewer_payload("output/simulation_results")` writes
the same page data as portable JSON, `simulation_results_viewer.json`.

```text
output/
  simulation_results_viewer.html  ← three.js viewer: textured globe, markers, trails, timeline
```

The page embeds three.js, the planet's surface texture from
`data/scene/textures/` (Earth, Mars, Venus, Titan and the Moon are registered
in its `manifest.toml`), the sidecar and a decimated copy of the trajectory, so
it opens from disk in a browser without a server. `max_frames` and
`data_budget_mb` bound the embedded rows, and the page's information panel
shows the resulting cadence. `trail_orbits` sets how much history trails
behind each spacecraft (three orbits by default; the page can change it).
Zooming in on a spacecraft, or selecting it and pressing F, replaces its marker
with the link-box assembly and its thruster and facet glyphs.
`models=Dict(id => "model.glb")` draws a CAD model instead; `stl=` is the older
spelling for STL files.

For a run with interlinks, the sidecar names the run's link schedule, and the
export adds its windows to the page as an `interlinks` block after checking the
schedule against the results manifest. A renderer that draws them shows each
link while it is active. The windows share `data_budget_mb` with the
trajectory, and `interlinks=false` leaves them out. A ground-station scene saved
without station geometry, by an earlier SpaceAGORA, cannot place the station
ends of its links, so exporting them is refused; pass `interlinks=false` to
export that scene's trajectories without the links.
[Interactive Visualization](visualization.md) covers the options, ensembles
and textures, and
[the independent renderer guide](https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl/blob/76f9037edb61e19e7314eb4db9bb7772a42f2ddc/viewer/README.md)
documents the controls and how to develop the viewer itself.

## Loading results in Julia

```julia
using CSV, DataFrames
df = CSV.read("output/simulation_results.csv", DataFrame)
```

The Feather file holds the same table in Arrow format and keeps the column
types:

```julia
using Arrow, DataFrames
df = DataFrame(Arrow.Table("output/simulation_results.feather"))
```

## Column reference

The meaning, units and frame of every column in the saved table are on
[Output Columns](output_columns.md).

## Multi-spacecraft runs

When several spacecraft are passed to `DynamicsModel`, each one contributes
its own set of columns: spacecraft 1 produces `sc1_*` columns, spacecraft 2
produces `sc2_*`, and so on. The number is the spacecraft's position in the
configuration, not its `SpacecraftModel.id`.

## Long layout for many spacecraft

A run of hundreds of spacecraft has thousands of `sc{i}_` columns. Set
`simulation_settings.results_layout = :long` to write the per-spacecraft
results as a long table instead, with one row per spacecraft and saved time and
the same values. For example, these settings go in the run's
`SimulationConfiguration`:

```julia
settings = SimulationSettings(results_layout=:long, results_directory="output/fleet")
```

```text
output/fleet/
  simulation_results_spacecraft.csv      ← one row per spacecraft and saved time
  simulation_results_spacecraft.feather
  simulation_results_run.csv             ← one row per saved time: the run-wide fields
  simulation_results_run.feather
  simulation_results.manifest.toml       ← schema_version "2", layout "long"
```

The spacecraft table starts with `spacecraft` (the spacecraft's position in the
configuration, the `i` of `sc{i}_`), `spacecraft_id` (its
`SpacecraftModel.id`), `plane` and `slot` (its constellation slot, `missing`
without one) and `time`. The per-spacecraft fields follow, named as on
[Output Columns](output_columns.md) without the `sc{i}_` prefix: `pos_1`,
`vel_1`, `mass` and so on. Rows are time-major: every spacecraft at the first
saved time, then every spacecraft at the next. Use `(spacecraft, time)` as the
key rather than the row order.

The run table holds `time` and the run-wide fields, such as `sun_dir_1..3` and
run-wide fields of your own. It is written only when the run has such fields;
a run without them removes the run table files that an earlier long run left
in the folder.

The manifest has `schema_version = "2"` and `layout = "long"` and lists the
files above. As in the wide layout, the CSVs follow `save_csv`, and the Feather
files and manifest follow the bundle switch. No `simulation_results.csv` or
`simulation_results.feather` is written, so a reader that expects the wide
files finds none. Wide files that an earlier wide run left in the folder stay
there; the manifest names this run's files. When two wide runs write the same
folder at once, the one that sees the other's CSV appear or change while it
writes also saves its own table under a unique name,
`simulation_results.<UTC time>.<process id>.<random>.csv`. The long layout
makes no such copy. In both layouts the run that finishes last leaves its
files.

In memory, `SimulationResults.table` is the spacecraft table, and
[`run_table`](@ref) returns the run table (`nothing` without run-wide fields).
For a wide run, [`long_table`](@ref) rearranges `table` into the long
spacecraft table.

A field's `per_satellite` flag decides which table it goes to, and its values
are flattened into columns as in the wide table, spacecraft by spacecraft. A
field whose values differ in shape between spacecraft, such as the thruster
levels of a spacecraft with thrusters next to one without, gets the columns of
every spacecraft, with `missing` where a spacecraft has none. The default
layout is `results_layout = :wide`, whose files are as described in
[Output files](#Output-files).

The long layout reserves `spacecraft`, `spacecraft_id`, `plane`, `slot` and
`time` for row identity in the spacecraft table, and `time` in the run table.
A custom field whose flattened column name collides with one of these is
refused with an `ArgumentError`; choose another `column_prefix` or keep the
wide layout. The same check applies when `long_table` converts wide results.

The long table and its Arrow file record which spacecraft own each column,
separately from the values. A field that is present but whose measurements are
all `missing` remains a channel with gaps, and an absent field does not become
a channel. This information travels in the Arrow metadata, so it stays with a
results file that is moved. Long files written without it, by development
versions before it existed, fall back to treating a column whose values are all
`missing` for one spacecraft as absent for that spacecraft. A CSV alone does
not carry this information.

### Streaming at each checkpoint

A long-layout run with checkpoints (see [Checkpointing](#Checkpointing)) that
writes its results with the bundle streams them instead of keeping every saved
row in memory until the end.

At each checkpoint the run appends the rows saved since the last one, except
the last row, to the spacecraft and run tables. The last row stays in memory
and is kept beside the checkpoint, so the run holds about one segment of rows
at a time. The tables are appended as Arrow IPC streams beside the checkpoint,
`simulation_checkpoint.results_spacecraft.arrows` and
`simulation_checkpoint.results_run.arrows`, and to the CSVs. While the run is
going, the manifest shows `streaming = true`, the number of saved times written
so far (`steps`) and the time reached (`time_reached_s`), and `Arrow.Table`
can read the rows written so far from the streams.

When the run ends, or a solve fails, the last rows are appended and the Feather
files and the manifest are written from the streams, as for a run that does not
stream. The files hold the same rows, values and column types, and the CSVs the
same bytes. Each Feather file has one record batch per append.
`SimulationResults.table` and [`run_table`](@ref) then read the Feather files
in place rather than holding a copy, so their columns are read-only; use
`DataFrame(results.table)` for a copy you can change.

Thinning (`SPACEAGORA_RESULTS_THIN_STRIDE`, `SPACEAGORA_RESULTS_FINAL_ONLY`)
counts rows over the whole run, so a streaming run keeps the same rows as a run
that does not stream.

The spacecraft that own each column stay fixed across checkpoints, including
across a resume. A change of ownership stops the run before the changed batch
is appended. The ownership record travels with the streams, the final Arrow
files and the returned tables; streams written without it fall back as
described above.

Each column keeps the type its first checkpoint's rows gave it, and later rows
must fit that type: numbers fit a column that also allows `missing`. Any other
change, such as a field that is `missing` through the first checkpoint and a
number after it, or a field whose flattened columns change, stops the run with
an error that names the column. A resume of such a run stops the same way.
Make the field return values of one type (NaN rather than `missing`, for
example), or run it without checkpoints or in the wide layout.

A resumed run cuts the streams and CSVs back to its checkpoint, dropping rows
that a crashed run wrote after it, and carries on appending, so its files hold
the whole run. A resume that solves nothing, because its checkpoint is at the
mission's end, leaves the files and streams as they are. The exception is a run
killed while it wrote its Feather files: a resume with that run's save fields
then writes them from the streams, whether or not it returns results. A resume
that solves but cannot include the earlier rows renames the streams beside the
checkpoint (`...before_resume_at_<t>s.arrows`), as it renames the rows file,
immediately before it writes its first replacement checkpoint. If it fails before
that checkpoint, the old checkpoint and its streams are still usable.

A run without checkpoints, a run that returns results without writing them, and
a run with `SPACEAGORA_SAVE_BUNDLE=0` keep their rows in memory, as a wide run
does.

The viewer reads both layouts. For a long run it reads the spacecraft table as
views, without copying it, and builds the same viewer frames as for the same
run in the wide layout. `run_monte_carlo_visualization` reads each sample's
wide results file, so each sample refuses `results_layout = :long`. The
telemetry verification study runs its cases with settings of its own, in the
wide layout.

## Output rate

The `data_rate` field of `MissionConfiguration` sets the output sample cadence
in seconds of simulated time; the default is `10.0`. A smaller `data_rate`
gives more rows and larger files, and a larger one gives smaller files with
coarser time resolution.

The `num_steps_to_save` field is kept for configuration compatibility. It does
not control buffering, flushing or peak memory use; `data_rate` sets the output
cadence.

## Checkpointing

Turn on periodic checkpoints in `SimulationSettings`. This fragment sets the
three checkpoint fields and leaves the others at their defaults:

```julia
SimulationSettings(
    checkpoint_enabled      = true,
    checkpoint_interval_s   = 300.0,         # every 300 s of simulated time
    checkpoint_directory    = "output/checkpoints"
)
```

The checkpoint directory defaults to `results_directory/checkpoints` when
`checkpoint_directory` is empty, and `checkpoint_interval_s` must be positive.
Each checkpoint writes `simulation_checkpoint.bin` (the state and run state)
and `simulation_checkpoint.manifest.toml`; a run that saves results also keeps
`simulation_checkpoint.results.bin`, and a run with interlinks
`simulation_checkpoint.interlinks.bin`, as described below.

A checkpointed run is solved in segments of `checkpoint_interval_s`. Its
navigation, guidance and control effectors tick at the same times as in a run
solved in one piece, including a tick that falls on a segment boundary, which
comes before the checkpoint written there. A tick and a boundary that differ
only by floating-point rounding, within 100 units in the last place, count as
one time, and the tick is taken at the boundary, once. For example, `3 * 0.1`
is 0.30000000000000004, and a 0.3 s interval ends the first segment at 0.3.
Likewise, a segment that would end that close to the run's end ends at the
run's end.

To resume a crashed or interrupted run, set `resume_from_checkpoint = true` in
`SimulationSettings`. The resumed run's effectors keep the same tick times,
counted from the run's start. A run resumed from a checkpoint written at a
segment boundary starts with the first tick after it. A run resumed from the
checkpoint written at the end of a run, to extend that run, first takes the
ticks due at that time, which the ended run did not take, and then saves the
row at that time again, as a run never interrupted saves it after those ticks.
Each checkpoint records which periodic models had not yet taken the tick due at
its time, also when the run stopped early there, and a run resumed from it
takes those models' ticks first, once. A checkpoint written by SpaceAGORA 0.4.0
or earlier holds no such record. It is read as that version left it, so a tick
an ulp or so after the checkpoint's time, which 0.4.0 left to the next segment,
is taken when the resumed run starts.

A resumed run's results hold the whole run, not only the part after the
resume. A checkpointed run that saves results (written, or returned with
`return_results=true`) keeps the rows it has saved beside its checkpoint, in
`simulation_checkpoint.results.bin`, and adds the new rows at each checkpoint.
A resumed run takes those rows back first. Its results files, its manifest
(`steps`, `first_time_s`) and its returned table then match those of a run
checkpointed at the same times and never interrupted, also after a crash, whose
partial results the resumed run completes. Rows that a crashed run wrote after
its last checkpoint are dropped, since the resumed run computes them again. A
resume from a checkpoint at or past the mission's end solves nothing and writes
no results files; its returned table holds the rows up to the checkpoint, or
none when the checkpoint holds no usable rows.

The link schedule holds the whole run too. A checkpointed run with interlinks
keeps its interlink history beside its checkpoint, in
`simulation_checkpoint.interlinks.bin`, and a resumed run takes it back first,
so its schedule keeps the windows recorded before the resume. With the same
configuration, the same checkpoint times and deterministic models whose
complete run state is restored, the resumed run's history and schedule match
those of the uninterrupted checkpointed run. A run without checkpoints can take
different adaptive solver steps, and so sample the schedule at different
times. A checkpoint with no usable history (written before SpaceAGORA kept it,
or whose history file is missing, damaged or another run's) resumes with a
warning: the schedule then covers only the resumed interval, and its
`[interlinks].history_start_s` is the resume time, `resumed_from_s`. If such a
resume solves nothing, its schedule is empty instead. A resume that changes
`record_available` or `record_terms` from the checkpoint's recording settings
is refused before any result stream changes, also when the history file is
missing.

The checkpoint also keeps the active interlink matching at its time, separately
from the history file, so a full receiving store keeps its carried-outflow
credit when scheduling resumes. Availability is checked again at the
checkpoint time. With backlog scheduling, a checkpoint written without that
matching warns and starts with zero receiver-outflow credit.

The rows file stays with the checkpoint files after the run; a later run in the
same checkpoint directory that is not a resume replaces it. It holds every row
the run saves, before any thinning, so it grows with the run. A long-layout run
that streams its results (see
[Streaming at each checkpoint](#Streaming-at-each-checkpoint)) keeps only its
last row there, and the rest in its results streams, which also stay beside the
checkpoint. A resume in the other layout cannot include rows kept that way: it
warns and writes to `resumed_at_<t>s`, as described below. A run that neither
writes nor returns results keeps no rows file, unless it runs the gravity
backbone, which builds its results from those rows. The checkpoint and its rows
file are written with Julia's `Serialization`, whose format can change between
Julia versions, so resume with the Julia version that wrote them. A resumed run
that cannot read the rows warns, as below.

A resumed run cannot include the earlier rows when its checkpoint holds none,
or when they were saved with other save fields or another number of
spacecraft. A checkpoint holds none when SpaceAGORA 0.5.1 or earlier wrote it
(which a run with guidance, navigation or control models refuses anyway, as
described below), or when the run that wrote it did not save results. The
resumed run then warns, and its results start at the resume. They are written
to `resumed_at_<t>s` inside `results_directory` (for example
`output/resumed_at_2400.0s/`), so the earlier results files stay as they were,
and the returned `configuration` names that folder as its `results_directory`.

The rows file beside the checkpoint is kept as well, renamed to
`simulation_checkpoint.results.before_resume_at_<t>s.bin`: after a run that was
killed, it may hold the only copy of the earlier rows, since such a run writes
no results files. A resumed run that saves no results renames the file the same
way rather than removing it. These renames happen when the first replacement
checkpoint is ready, so a resume that fails before then leaves the old
checkpoint and its rows usable, and a resume that solves nothing renames
nothing. If a kept name already exists, a numbered name is used. A normal resume
does not look at kept files, and restoring the original save fields does not
bring them back. Keep a kept rows file together with its checkpoint if you plan
a manual recovery; there is no public reader for the format.

Besides the state, every checkpoint holds the run state that the callbacks keep
on the spacecraft: the batteries' empty and full flags and shedding tier, the
data stores' empty and full flags, the actuator faults' activity and held
levels, the sensors' random streams, latest measurements and gyro biases, the
hinges' and gimbals' angles, and the ramp factor and held thrust of thrusters
with firing dynamics. A run with batteries, faults, sensors, data stores or
thrusters with firing dynamics refuses a checkpoint written without them, or
for a different configuration of them; a run with firing dynamics therefore
refuses a checkpoint written by SpaceAGORA 0.4 or earlier. A checkpoint without
hinge or gimbal angles resumes with a warning, from their current angles. A
guidance, navigation or control model, or a sensor kind, that keeps run state
of its own adds it to the file through `checkpoint_state` (see
[Run state in checkpoint files](../extensibility.md#Run-state-in-checkpoint-files));
without it, a resumed run starts that state from its reset.

A run with the burn scheduler (`BaseThrusterModel`) also records the latest
maneuver commands, the accepted burn plans and which requests have been
consumed. New requests are processed on the configured control schedule, after
any guidance due at the same time. A checkpoint segment or an atmosphere exit
does not re-plan an accepted burn. The resumed run keeps the burn's window and
takes any pending control tick once, and an accepted request is not flown again
after it completes.

Guidance requests are identified by their Δv, direction and source orbit, so a
new source orbit allows a new request even when the Δv is the same. A
configured Δv is consumed once unless it changes. An explicit `BurnPlanFamily`
installation supplies its own window; the engine registers its boundaries after
guidance, also between control ticks, and keeps these plans separate from
ordinary requests.

The burn scheduler refuses checkpoints written before it processed requests on
the control schedule, also when no burn was configured: their pending clock
positions and consumed requests cannot be inferred safely. Resume such a
checkpoint with the SpaceAGORA version that wrote it, or start a new run. The
refusal comes before any result stream changes.

The checkpoint also records the orbit each spacecraft is on, once one is past
its first, so a resumed run continues the count. A `MissionOrbits` run then
ends after the same number of orbits, and the campaign maneuver guidance, the
pass-indexed density models and the periapsis pulses see the orbits they would
have seen. A checkpoint written by an earlier SpaceAGORA holds no count, and a
run resumed from it counts from the first orbit again. If such a checkpoint was
written after the first orbit, the resumed count is behind: a `MissionOrbits`
run ends late, and the campaign maneuver guidance repeats the maneuvers of
orbits already flown. Run it again from its start, or from a checkpoint the
current version wrote.

With the atmosphere-entry limit (`SPACEAGORA_ENTRY_TARGET_COUNT`), the
checkpoint also keeps each spacecraft's count of completed entries, once one
has entered, so a resumed run ends at the same entry as an uninterrupted
checkpointed run. A checkpoint written without this count restarts the count
at zero, so a run resumed after an entry can end late. Rerun from the start to
get a checkpoint that includes the count.

The checkpoint also records which spacecraft are active. A spacecraft lost to
an impact or a touchdown before the checkpoint stays inactive in the resumed
run, as it would in a run that was never interrupted: it is not propagated, its
guidance, navigation and control models treat it as inactive (`skip_inactive`),
and its impact or touchdown does not happen again. A checkpoint written by an
earlier SpaceAGORA records this only for an interlinked run; any other run
resumed from such a checkpoint starts with every spacecraft active.

A run with scoped dynamic effectors (`scoped`) also records each effector's
kind and the run positions it acts on, and resumes only with the same scopes.
It refuses a checkpoint whose scopes differ, and a checkpoint written without
scoped effectors (by a run without them, or by an earlier SpaceAGORA). A run
without scoped effectors adds nothing to its checkpoint, and refuses a
checkpoint written by a run with them.

Every checkpoint also records the spacecraft each guidance, navigation and
control model acts on: its `spacecraft_scope` as the run-start initializers
left it. That is every spacecraft, the run positions of a bound model
(`BoundSpacecraft`, which a model declared on a spacecraft has once it is bound
to that spacecraft), or a coupled model's primary and partners. A resumed run
refuses a checkpoint whose models act on other spacecraft, or that has a
different number of guidance, navigation or control models, before any model's
run state is restored and whether or not the models keep run state. A momentum
manager's saved command therefore never goes to a manager of another
spacecraft. A checkpoint written by SpaceAGORA 0.5.1 or earlier holds no such
record, so a run with any guidance, navigation or control model refuses it
before anything is restored: resume it with SpaceAGORA 0.5.1, or run from the
start. A run without those models resumes from it.

## Campaign and multi-process runs

`run_monte_carlo` calls `run_simulation` once for each sample, with the
`SimulationConfiguration` your per-seed function returns, so each sample's
`results_directory` and `SimulationSettings` are whatever that configuration
sets. `run_constellation_ensemble` runs each spacecraft as its own member; when
result saving or checkpoints are on, each member writes to a
`sat_<index>_id_<spacecraft id>` folder inside `results_directory`. In both
cases the file layout above applies to each sample or member, whether the
outer route is serial, threaded or process-backed. See
[Parallel Execution](parallel_execution.md) for how to choose and configure the
outer route.

`run_constellation_ensemble(args; merge_results=true)` also merges its members
into one table with the long layout's columns. With `results=true` it writes
`constellation_members.feather` (and `.csv` with `save_csv`) and
`constellation_summary.csv` beside the member folders. See
[One table for the whole ensemble](constellations.md#One-table-for-the-whole-ensemble).


## In-memory constellation recorder

For output at a fixed cadence without keeping the solver's full step history,
attach the optional recorder to an existing `SimulationConfiguration` named
`args`:

```julia
using SpaceAGORA
using SpaceAGORA.SimulationModel: TrajectoryRecorder, get_trajectory_recorder_callback,
    trajectory_times, trajectory_positions, trajectory_save_data

recorder = TrajectoryRecorder(args)
run_simulation(args; return_solution=false,
    extra_callbacks=(get_trajectory_recorder_callback(recorder),))
times = trajectory_times(recorder)
positions = trajectory_positions(recorder)  # component x spacecraft x sample
```

The cadence defaults to `args.mission_configuration.data_rate`, in seconds.
The recorder preallocates its arrays and grows them when needed, so its memory
still grows with the number of spacecraft and recorded samples. Read the
returned views after the solve, and take new views after resetting or reusing
the recorder. `trajectory_save_data(recorder)` turns the recorded samples into
ordinary saved snapshots when you need them; that conversion builds one
dictionary per sample.

The recorder does not change the run's file output; use
`SimulationSettings(results=false)` when the in-memory result is all you need.
With an explicit `save_fields` list, every field is recorded through its own
getter, including a custom getter that uses a built-in field's name. With the
default fields, those that have a typed array filler (position, velocity, the
geodetic coordinates, periapsis altitude, mass, wind, the aerodynamic forces,
heat rate and load, and the attitude quaternion) go straight into typed arrays,
and the other default fields are recorded through their getters.

### Backlog scheduling and stored data

`DataBacklogPolicy` scores a candidate link with
`min(link_rate, held/horizon + current_production, free_space/horizon + current_outflow)`.
The receiver's outflow is sampled before the new selection, from links that
were already active and are still available, with the same carried-rate caps
as the runtime. It is an estimate for scoring; the data-flow calculation still
keeps stores from overfilling. The saved transfer-rate and stored and lost bit
columns report the flows actually carried, not this estimated score. A full
store with no outflow still scores zero room.
