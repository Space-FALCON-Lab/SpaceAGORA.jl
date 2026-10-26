# Simulation Outputs

Use this page when you want to understand the files written by
`run_simulation`, the column names in the CSV, or how to load results for
post-processing.

This page is for users who have already run a simulation and now need to
interpret the output.

Shortest successful command:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Output Columns](output_columns.md)
- [The Integrated State](integrated_state.md) (which columns are integrated and which are derived)
- [Simulation Configuration](simulation_configuration.md)
- [Verification Study](verification_study.md)
- [Recipes](recipes.md)

## Output files

With the default output settings, `run_simulation` writes three files to the
`results_directory`:

```text
output/
  simulation_results.csv         ← tabular results, one row per time step
  simulation_results.feather     ← same data in Apache Arrow / Feather format
  simulation_results.manifest.toml ← metadata: row count, spacecraft count, file hashes
```

A run with interlinks also writes its link schedule,
`simulation_results_link_schedule.csv` and
`simulation_results_link_schedule.feather` (see
[the link schedule's columns](output_columns.md#Link-schedule)). The manifest lists
both under `[files]`, and describes the schedule under `[interlinks]`:

- `policy` and `score_unit`: the scheduling policy's type and its score's unit;
- `rows`: the schedule's number of rows;
- `record_available` and `record_terms`: whether the history recorded the
  available links (the contact rows) and a weighted policy's term scores;
- `history_start_s` and `history_end_s`: the times of the interlink history's
  first and last samples, the span the schedule covers (absent when it has no
  samples);
- `terms`, for a `WeightedPolicy`: each top-level term's `policy`, `score_unit`,
  `weight` and `scale`, in the order of the `term<k>_start` columns.

Set `simulation_settings.save_csv = false` to save the Feather file and manifest
without writing a CSV. A CSV left by an earlier run is preserved, but is not
listed in the new manifest or in `SimulationResults.files`.
Set `simulation_settings.results = false` to disable result-file output.
The manifest records
`schema_version`, `created_utc`, `mission_time_s`, `steps`,
`spacecraft_count`, `first_time_s` (the first row's time), and SHA-256 hashes
for each data file. A resumed run's manifest also records `resumed_from_s`, the
checkpoint time it resumed from (see [Checkpointing](#Checkpointing)).
`created_utc` is the wall-clock time the file was written, not the simulation
epoch. The epoch is the `InitialTime` you configured; the manifest does not
store it, so keep it with your script or read it from the returned
configuration.

The developer override `SPACEAGORA_SAVE_BUNDLE=0` disables Feather and manifest
output. With this override and `save_csv=false`, use `return_results=true` to
receive the table in memory. Explicitly enabled scene sidecars and checkpoints
can still be written.

## Getting the results in memory

Pass `return_results=true` to get a [`SimulationResults`](@ref) back instead of
reading files. Its `table` is the same `DataFrame` the CSV would contain
(same columns and values), and it is built even when
`simulation_settings.results = false`, in which case result files are not written.
Explicitly enabled checkpoints can still be written; they are listed in `files`
only when this call writes them. Resuming an already completed checkpoint does
not report the old checkpoint files as new output, and writes no results files:
its `table` holds the whole run from the rows kept beside the checkpoint, or is
empty when the checkpoint holds no usable rows.
`configuration` is the configuration that actually ran: under the default
`isolate_state=true` that is the deep copy, so controller state and logs, such
as the RPO command log, are reachable there and not on the object you passed
in. `files` lists the output files the run wrote, and `solution` holds the ODE
solution when `return_solution=true` is passed as well (otherwise `nothing`).
`return_results=true` cannot be combined with `return_solver_metadata=true`.

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

Most repository examples build their configuration with `make_example_config`,
whose `results_directory` is `<repository root>/output` unless
`SPACEAGORA_CLI_OUTPUT_DIR` is set; the CLI's `--output-dir` sets exactly that
variable. With the default result settings, the three file names do not
change, so:

- two runs using the same results directory overwrite those files;
- `julia --project=. src/interfaces/cli/main.jl run --example=<script> --output-dir=output/<name>`
  gives each run its own directory, for the scripts that take their directory
  from `make_example_config`: the first-run scripts (`AGORA_Basic_Quickstart.jl`,
  `AGORA_Earth_NoGRAM.jl`), the controls and torque
  tests, and the mission scripts `AGORA_Basic_GRAMEarth.jl`, `AGORA_Odyssey.jl`,
  `AGORA_Vex.jl`, `AGORA_Titan.jl`, `AGORA_Magellan.jl`, `AGORA_LOFTID.jl`,
  `AGORA_Mars_NoGRAM.jl` and `Earth_Thruster_Test.jl`;
- `--output-dir` has no effect on the scripts that choose their own directory:
  `AGORA_Earth.jl` and `AGORA_Keplerian.jl` write to `output/` directly,
  `Earth_Navigation.jl` to `output/observer_noDA/`,
  `AGORA_Earth_Aerobraking.jl` to `output/earth_aerobraking/`, and the RPO,
  robot-arm and compliant multibody demos to their own `output/<demo>/` directories; running
  one of them twice overwrites its previous results;
- `odyssey_surrogate.jl` runs in its own environment and writes its own files:
  the comparison writes to `odyssey_surrogate_results/` in the current
  directory unless you pass `--output=DIR`, refuses an existing directory, and
  prints the absolute path when it finishes (see the
  [Odyssey walkthrough](../tutorials/odyssey_surrogate.md));
- for your own scripts, pass `results_directory=` to `make_example_config` or
  set it on `SimulationSettings`.

The Mars RAAN campaign and its output convention belong to
SpaceAGORAMissionReconstructions.
Its launcher and output paths are documented there.

The quickstart example additionally saves four PNG plots under
`<results_directory>/plots/`. The RPO examples also generate HTML plots.
`Solar_Panel_Cloth_Deployment_Demo.jl` writes four HTML files in its demo
directory instead of the three simulation result files.

`AGORA_Earth_MonteCarlo.jl` prints its successful and failed sample counts and
elapsed time to the terminal. It disables result saving, so it writes no
result files, including when `--output-dir` is supplied.

Smoke mode (`--smoke` on the CLI, or `SPACEAGORA_EXAMPLE_SMOKE=1` for a script)
uses a short mission, normally at most 120 s and one orbit. An explicit
`--output-dir` is honored through `SPACEAGORA_CLI_OUTPUT_DIR`. Without that
override, smoke mode writes to `output/` under the current working directory.
Choose a separate destination to keep earlier results: reusing the same
directory replaces its previous result files. Results are kept only when
`SPACEAGORA_EXAMPLE_SMOKE_RESULTS=1`, which the CLI sets for you.

### Visualization scene sidecar (opt-in)

Setting `simulation_settings.save_visualization_scene = true` adds one more
file and one more group of columns, both consumed by the post-hoc 3D viewer
described in [Interactive Visualization](visualization.md):

```text
output/
  simulation_results_scene.json  ← planet, epoch, spacecraft geometry, frame rotation samples
```

The sidecar carries everything the viewer needs besides the trajectory rows:
the central body's radii, spin and texture key, a table of J2000-to-body-fixed
quaternions sampled over the mission span (so the viewer needs no SPICE
kernels), and every spacecraft's links as boxes built from `Link.dims`, with
thruster, facet and joint glyphs. It is written after the results bundle by
the same run, and only when the flag is set; the default output is unchanged.

With the flag on, spacecraft that have non-root links also get the
`link_pose` save field: columns `sc{N}_link_pose_1` .. `sc{N}_link_pose_7n`
holding `[rx, ry, rz, qx, qy, qz, qw]` per non-root link, in the link order
recorded in the sidecar. These are the link poses relative to the root bus
at each saved step, so articulated panels play back as they moved. A
spacecraft that carries a thruster gimbal also gets the `thruster_direction`
field (`sc{N}_thruster_direction_1` .. `_3k`): every thruster's exhaust
direction in its link's frame at each saved step, so the viewer turns the
gimballed thruster's cone and plume as the gimbal moved. Reading the sidecar
back:

```julia
using SpaceAGORA
scene = read_visualization_scene("output/simulation_results_scene.json")
scene.planet.name, length(scene.spacecraft), scene.spacecraft[1].bounding_radius_m
```

The flag adds four more fields when the run can supply them, each listed in
on [Output Columns](output_columns.md#Visualization-fields):

- `thruster_level`: `sc{N}_thruster_level_{k}`, the firing level of thruster
  `k` from 0 to 1, for a spacecraft with thrusters whose control effector
  reports per-thruster levels. Thrusters are numbered in the sidecar's order:
  the spacecraft's links in order, each link's thrusters in order.
- `sun_dir`: `sun_dir_1` .. `sun_dir_3`, the unit vector from the planet's
  centre to the Sun in the inertial frame of the saved positions. One
  direction per row, shared by every spacecraft, written when the run's
  ephemerides can locate the Sun at the start epoch. The viewer's sun lighting
  reads it.
- `arm_pose`: `sc{N}_arm_pose_1` .. `sc{N}_arm_pose_7m` for a spacecraft with
  a compliant robot-arm plan, `[rx, ry, rz, qx, qy, qz, qw]` per arm link: the
  link centre of mass relative to the spacecraft position, in inertial metres,
  and the link's inertial attitude quaternion.
- Articulated spacecraft (any non-fixed `Joint`) add four fields to the default
  set: `joint_q` (`sc{N}_joint_q_1` ..) and `joint_qd` (`sc{N}_joint_qd_1` ..)
  list the non-fixed joints in `spacecraft.joints` order (a ball joint takes four
  and three columns); `articulated_link_pose` (`sc{N}_articulated_link_pose_1` ..
  `7n`) is `[rx, ry, rz, qx, qy, qz, qw]` per link in `spacecraft.links` order:
  the link's COM in inertial metres and its scalar-last body-to-inertial
  attitude, from the joint kinematics (links merged through `:fixed` joints
  included); `system_com` (`sc{N}_system_com_1..3`) is the inertial center of
  mass of the whole spacecraft. The `link_pose` viewer field still holds the
  configured rigid geometry.
- Spacecraft with compliant attachments add `attachment_pose`
  (`sc{N}_attachment_pose_1` .. `7N`): `[rx, ry, rz, qx, qy, qz, qw]` per attachment body in
  attachment order, the body COM in inertial metres (`pos + att_r`) and its scalar-last body-to-inertial
  attitude. `system_com` (`sc{N}_system_com_1..3`) includes the attachment bodies (rigid
  spacecraft get the column from the attachment fields; articulated ones from the articulated
  fields). Spacecraft without attachments save neither.
- `density`: `sc{N}_density`, the atmospheric density along the trajectory in
  kg/m³, for a run with an atmosphere model.

An explicit `save_fields` list receives the visualization fields it lacks,
appended once after the fields you named. The plume columns below are not
visualization fields: they belong to the default field list whenever the
plume-surface effector is in the run, with or without the sidecar, and an
explicit list that omits them does not write them.

### Interactive viewer page

With the sidecar present, `export_visualization` builds a self-contained
HTML page next to the bundle, or pass `visualization=true` to `run_simulation`
to do both in one go:

```julia
run_simulation(args; visualization=true)
# or, from an existing flagged run
export_visualization("output/simulation_results"; max_frames=2000, frame=:inertial)
```

Setting `SPACEAGORA_VISUALIZATION=1` in the environment does the same for a
script you do not want to edit, for example
`SPACEAGORA_VISUALIZATION=1 julia --project=examples examples/AGORA_Earth_NoGRAM.jl`.

```text
output/
  simulation_results_viewer.html  ← three.js viewer: textured globe, markers, trails, timeline
```

The page embeds three.js, the surface texture from `data/scene/textures/`
(Earth, Mars, Venus, Titan, Moon are registered in `manifest.toml`), the
sidecar and a decimated copy of the trajectory, so it opens from disk in any
modern browser without a server. `max_frames` and `data_budget_mb` bound the
embedded rows; the info panel shows the effective cadence. `trail_orbits`
sets how much history trails behind each spacecraft (default three orbits,
selectable on the page). Zooming in on a spacecraft, or selecting it and
pressing F, swaps its marker for the link-box assembly with thruster and
facet glyphs; `stl=Dict(id => "model.stl")` draws a CAD mesh instead.
For a run with interlinks, the sidecar names the run's link schedule, and the
exported payload carries its windows (an `interlinks` block, checked against the
results manifest): a renderer that draws them shows each link while it is active.
They share `data_budget_mb` with the trajectory, and `interlinks=false` leaves
them out. Older saved ground-station scenes may lack station geometry; exporting
their links is refused because their endpoints cannot be placed. Pass
`interlinks=false` to export those scenes' trajectories without links.
The [Interactive Visualization](visualization.md) page covers the options,
ensembles and textures; [the independent renderer guide](https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl/blob/76f9037edb61e19e7314eb4db9bb7772a42f2ddc/viewer/README.md) documents the controls and how to
develop the viewer itself.

## Loading results in Julia

```julia
using CSV, DataFrames
df = CSV.read("output/simulation_results.csv", DataFrame)
```

Or from the Feather file (faster for large runs):

```julia
using Arrow, DataFrames
df = DataFrame(Arrow.Table("output/simulation_results.feather"))
```

## Column reference

The meaning, units and frame of every column in the saved table are on
[Output Columns](output_columns.md).

## Multi-spacecraft runs

When multiple spacecraft are passed to `DynamicsModel`, each spacecraft
contributes its own column set. Spacecraft 1 produces `sc1_*` columns,
spacecraft 2 produces `sc2_*`, and so on. The spacecraft position in the
configuration sets the numeric suffix; `SpacecraftModel.id` does not.

## Long layout for many spacecraft

A run of hundreds of spacecraft has tens of thousands of `sc{i}_` columns. Set
`simulation_settings.results_layout = :long` to write the per-spacecraft
results as a long table instead: one row per spacecraft and saved time, with
the same values.

```julia
settings = SM.SimulationSettings(results_layout=:long, results_directory="output/fleet")
```

```text
output/fleet/
  simulation_results_spacecraft.csv      ← one row per spacecraft and saved time
  simulation_results_spacecraft.feather
  simulation_results_run.csv             ← one row per saved time: the run-wide fields
  simulation_results_run.feather
  simulation_results.manifest.toml       ← schema_version "2", layout "long"
```

- **The spacecraft table** starts with `spacecraft` (the spacecraft's position in
  the configuration, the `i` of `sc{i}_`), `spacecraft_id` (its
  `SpacecraftModel.id`), `plane` and `slot` (its constellation slot, `missing`
  without one) and `time`. The per-spacecraft fields follow, named as on
  [Output Columns](output_columns.md) without the `sc{i}_` prefix: `pos_1`,
  `vel_1`, `mass` and so on. Its rows are time-major: every spacecraft at the
  first saved time, then at the next. Use `(spacecraft, time)` as the key rather
  than the row order.
- **The run table** holds `time` and the run-wide fields, such as `sun_dir_1..3`
  and run-wide fields of your own. It is written only when the run has one, and
  a run without one removes the run table files an earlier long run left there.
- **The manifest** has `schema_version = "2"` and `layout = "long"`, and lists
  the files above. The CSVs follow `save_csv`, the Feather files and manifest
  follow the bundle, as in the wide layout. No `simulation_results.csv` or
  `.feather` is written, so a reader expecting the wide file finds none. Wide
  files left in the folder by an earlier wide run stay there; the manifest names
  this run's files. When two runs write the same folder at once, the wide layout
  keeps the other run's CSV under a second name; the long layout does not, and
  the run that finishes last leaves its files.
- **In memory**, `SimulationResults.table` is the spacecraft table and
  [`run_table`](@ref) returns the run table (`nothing` without run-wide fields).
  For a wide run, [`long_table`](@ref) rearranges `table` into the long
  spacecraft table.

Which table a field goes to follows its `per_satellite` flag, and its values
are flattened into columns as in the wide table, spacecraft by spacecraft. A
field whose values differ in shape between spacecraft, such as the thruster
levels of a spacecraft with thrusters beside one without, has the columns of
every spacecraft, `missing` where a spacecraft has none. The default stays
`results_layout = :wide`, whose files are unchanged.

Long output reserves `spacecraft`, `spacecraft_id`, `plane`, `slot` and `time`
for spacecraft row identities, and `time` in the run table. A custom field whose
flattened name collides with one of these columns is refused with an
`ArgumentError`; choose another `column_prefix` or keep the wide layout.
This also applies when converting wide results with `long_table`.

The long table and its Arrow file retain which spacecraft own each column,
separately from its values. A present channel whose measurements are all `missing`
remains a channel with gaps; an absent field does not become a channel. This
schema information is carried in Arrow metadata, so it stays with a moved results
file. Older draft long files without it retain the earlier fallback: columns
whose values are all `missing` on one spacecraft are treated as absent. CSV alone
does not retain this schema metadata.

### Streaming at each checkpoint

A long-layout run with checkpoints (see [Checkpointing](#Checkpointing)) that
writes results streams them: it does not keep its saved rows in memory until it
ends.

- **At each checkpoint**, the rows saved since the last one, except the last row,
  are appended to the spacecraft and run tables. The last row stays in memory
  and is kept beside the checkpoint, so the run holds about one segment of rows
  at a time.
- **The tables grow as the run goes.** They are appended to as Arrow IPC streams
  beside the checkpoint (`simulation_checkpoint.results_spacecraft.arrows` and
  `simulation_checkpoint.results_run.arrows`) and to the CSVs. The manifest
  shows `streaming = true`, the number of saved times written and the time
  reached. Until the run finishes, `Arrow.Table` reads its rows from the
  streams.
- **When the run ends, or a solve fails,** the last rows are appended and the
  Feather files and the manifest are written from the streams, as for a run that
  does not stream. The files hold the same rows, values and column types, and
  the CSVs the same bytes. A Feather file has one record batch per checkpoint.
- **`SimulationResults.table` and [`run_table`](@ref)** read the Feather files
  in place instead of holding a copy. Their columns are read-only; use
  `DataFrame(results.table)` for a copy you can change.
- **Thinning** (`SPACEAGORA_RESULTS_THIN_STRIDE`, `SPACEAGORA_RESULTS_FINAL_ONLY`)
  counts rows over the whole run, so it keeps the rows a run that does not
  stream keeps.
- **Column ownership.** The spacecraft owning each column stay fixed across
  checkpoints, including a resume. Changing that ownership stops the run before
  appending the changed batch. Ownership metadata stays with the streams, final
  Arrow files and returned tables. Older draft streams without it keep the
  legacy fallback described above.
- **Column types.** Each column keeps the type its first checkpoint's rows gave
  it, and later rows must fit it: numbers fit a column that also allows
  `missing`. A column whose type changes otherwise (a field that is `missing`
  through the first checkpoint and a number after it, or whose flattened
  columns change) stops the run with an error naming the column, and so does a
  resume of it. Make such a field return values of one type (NaN rather than
  `missing`, say), or run it without checkpoints or in the wide layout.
- **A resumed run** cuts the streams and CSVs back to its checkpoint, so rows
  written after it by a run that crashed are dropped, and goes on appending. Its
  files then hold the whole run. A resume that solves nothing (its checkpoint
  is at the mission's end) leaves the files and the streams as they are, unless
  the run was killed while writing its Feather files: then a resume with the
  run's save fields writes them from the streams, whether or not it returns
  results. A resume that solves and cannot include the earlier rows keeps the
  streams beside the checkpoint under other names
  (`...before_resume_at_<t>s.arrows`), as it keeps the rows file, immediately
  before its first replacement checkpoint is written. A failure before that
  checkpoint leaves the old checkpoint and its streams usable.

A run without checkpoints, a run that returns results without writing them, and
a run with `SPACEAGORA_SAVE_BUNDLE=0` keep their rows in memory, as a wide run
does.

The viewer reads both layouts. For a long run it reads the spacecraft table as
views, without copying it, and its page equals that of the same run in the wide
layout. `run_monte_carlo_visualization` reads each sample's wide results file
and refuses `results_layout = :long`. The telemetry verification study runs its
cases with settings of its own, in the wide layout.

## Output rate

The `data_rate` field in `MissionConfiguration` sets the output sample cadence
in seconds of simulated time. The default is `10.0` seconds. Decreasing
`data_rate` produces more rows and larger files; increasing it reduces output
size but loses temporal resolution.

The `num_steps_to_save` field is retained for configuration compatibility. It
does not currently control buffering, flushing or peak memory use. Use
`data_rate` to set the output cadence.

## Checkpointing

Enable periodic checkpointing via `SimulationSettings`:

```julia
SM.SimulationSettings(
    checkpoint_enabled      = true,
    checkpoint_interval_s   = 300.0,         # every 300 s of simulated time
    checkpoint_directory    = "output/checkpoints"
)
```

A checkpointed run is solved in segments of `checkpoint_interval_s`. Its
navigation, guidance and control effectors tick at the same times as in a run
solved in one piece, including a tick that falls on a segment boundary, which
comes before the checkpoint written there. A tick and a boundary that differ only
by floating-point rounding, within 100 units in the last place (`3 * 0.1` is
0.30000000000000004, and a 0.3 s interval ends the first segment at 0.3), count
as one time: the tick is taken at the boundary, once. Likewise, a segment that
would end that close to the run's end ends at the run's end.

A crashed or interrupted run can be resumed by setting
`resume_from_checkpoint = true` in `SimulationSettings`. The resumed run's
effectors keep the same tick times, from the run's start. A run resumed from a
checkpoint written at a segment boundary starts with the first tick after it;
one resumed from the checkpoint written at the end of a run, to extend that
run, first takes the ticks due at that time, which the ended run did not take.
Each checkpoint records which periodic models had not yet taken the tick due at
its time, also when the run stopped early there, and a run resumed from it takes
those models' ticks first, once. A checkpoint written by SpaceAGORA 0.4.0 or
earlier holds no such record; it is read as that version left it, which also
means a tick an ulp or so after its time, which 0.4.0 left to the next segment, is
taken when the resumed run starts.

A resumed run's results hold the whole run, not only the part after the resume.
A checkpointed run that saves results (written or returned with
`return_results=true`) keeps the rows it has saved beside its checkpoint, in
`simulation_checkpoint.results.bin`, and adds the new rows at each checkpoint. A
resumed run takes those rows back first. Its results files, its manifest
(`steps`, `first_time_s`) and its returned table then match those of a run
checkpointed at the same times and never interrupted, including after a crash,
whose partial results the resumed run completes. Rows written after the last
checkpoint by a run that crashed are dropped, since the resumed run computes them
again. The link schedule holds the whole run too: a checkpointed run with
interlinks keeps its interlink history beside its checkpoint, in
`simulation_checkpoint.interlinks.bin`, and a resumed run takes it back first,
so its schedule retains the earlier recorded windows. With the same configuration,
checkpoint schedule and deterministic models whose complete run state is restored,
its history and schedule match the uninterrupted checkpointed run. A run with no checkpoints can have different
adaptive solver steps and therefore different scheduling samples. A checkpoint that
holds no record of that history (written before SpaceAGORA kept it), or whose
history file is missing or another run's, resumes with a warning: its schedule
then covers the resumed interval only, and its `[interlinks].history_start_s` is
the resume's time, `resumed_from_s`. If it solves nothing, its schedule is empty
instead. Changing `record_available` or `record_terms` from the checkpoint's
recording settings is refused before changing result streams, including when the
history sidecar is missing. A run extended from the end of a finished run takes the ticks due there
first (see above), and saves the row at that time again after them, as a run never
interrupted saves it. A resume from a checkpoint at or past the mission's end
solves nothing and writes no results files; its returned table holds the rows up
to the checkpoint, or none when the checkpoint holds no usable rows.

The checkpoint also retains the previous active interlink matching independently
of the history sidecar, so a full receiving store keeps its carried-outflow credit
when scheduling resumes. Availability is checked again at the checkpoint time.
With backlog scheduling, older checkpoints without that matching warn and begin
with zero receiver-outflow credit.

The rows file stays with the checkpoint files after the run; a later run in the
same checkpoint directory that is not a resume replaces it. It holds every row the
run saves, before any thinning, so it grows with the run, to about the size of the
results CSV. A long-layout run that streams its results (see
[Streaming at each checkpoint](#Streaming-at-each-checkpoint)) keeps only its last
row there: the rest is in its results streams, which stay beside the checkpoint
too. A resume in another layout cannot include rows kept that way: it warns and
writes to `resumed_at_<t>s`, as below. A run that neither writes nor returns results keeps no rows file,
except one with the gravity backbone, which builds its results from those rows. The
checkpoint and its rows file are written with Julia's `Serialization`, whose
format can change between Julia versions: resume with the Julia version that wrote
them. A resumed run that cannot read the rows warns, as below.

A resumed run cannot include the earlier rows when its checkpoint holds none (a
checkpoint written by SpaceAGORA 0.5.1 or earlier, which a run with guidance,
navigation or control models refuses, as described below, or by a run that did not
save results), or when they were saved with other save fields or another number of
spacecraft. It then warns, and its results start at the resume: they are written
to `resumed_at_<t>s` inside `results_directory` (for example
`output/resumed_at_2400.0s/`), so the earlier results files are kept as they were.
The rows file beside the checkpoint is kept too, renamed to
`simulation_checkpoint.results.before_resume_at_<t>s.bin`: after a run that was
killed, it may hold the only copy of the earlier rows, since such a run writes no
results files. A resumed run that saves no results renames it the same way rather
than removing it. These renames happen when the first replacement checkpoint
is ready, so a resume that fails before then leaves the old checkpoint and its
rows usable. A resume that solves nothing does not rename them. If a kept name
already exists, a numbered name is used. Normal resume does not search those kept
files: restoring the original save-field settings alone does not recover them.
Keep the rows file and its matching checkpoint together for an explicitly
reviewed recovery procedure; there is no public reader for the kept-file format.
The returned `configuration` names that folder as its `results_directory`.
The checkpoint directory defaults to `results_directory/checkpoints` when left
empty. Besides
the state, the checkpoint holds the run state the callbacks keep on the
spacecraft: the batteries' empty and full flags and shedding tier, the
actuator faults' activity and held levels, the sensors' random streams,
latest measurements and gyro biases, the hinges' and gimbals' angles, and the
ramp factor and held thrust of thrusters with firing dynamics. A run with
batteries, faults, sensors or thrusters with firing dynamics refuses a
checkpoint written without them, or for a different configuration of them; a
run with firing dynamics therefore refuses a checkpoint written by SpaceAGORA
0.4 or earlier. A guidance, navigation or
control model, or a sensor kind, that keeps run state of its own adds it to
the file through `checkpoint_state` (see
[Run state in checkpoint files](../extensibility.md#Run-state-in-checkpoint-files));
without it, a resumed run starts that state from its reset.

A run with the burn scheduler (`BaseThrusterModel`) also records the latest
maneuver commands, accepted burn plans and which requests have been consumed.
New requests are processed on the configured control schedule, after guidance due
at the same time. A checkpoint segment or atmosphere exit does not re-plan an
accepted burn. The resumed run preserves its window and takes any pending control
tick once; the same accepted request is not flown again after completion.

Guidance requests are identified by their Δv, direction and source orbit. A new
source orbit permits a new request even when its Δv is unchanged. A configured
Δv is consumed once unless changed. An explicit `BurnPlanFamily` installation
supplies its own window; the engine registers its boundaries after guidance, even
between control ticks. These plans remain separate from ordinary requests.

Burn-scheduler checkpoints written before this control-clock policy are refused:
their pending clock positions and consumed requests cannot be inferred safely.
Resume those checkpoints with the SpaceAGORA version that wrote them, or start a
new run. The refusal occurs before result streams are changed.

The checkpoint also records the orbit each spacecraft is on, once one is past its
first, so a resumed run continues the count: a `MissionOrbits` run ends after the
same number of orbits, and the campaign maneuver guidance, the pass-indexed
density models and the periapsis pulses see the orbits they would have seen. A
checkpoint written by an earlier SpaceAGORA holds no count, and a run resumed from
it counts from the first orbit again, as before. If such a checkpoint was written
after the first orbit, the resumed run's count is behind: a `MissionOrbits` run
ends late, and the campaign maneuver guidance repeats the maneuvers of the orbits
already flown. Run it again from its start, or from a checkpoint this version
wrote.

The atmosphere-entry limit (`SPACEAGORA_ENTRY_TARGET_COUNT`) also retains each
spacecraft's completed entry count, once one has entered. A resumed run then
ends at the same entry as an uninterrupted checkpointed run. Older checkpoints
without this count retain their earlier behavior: counting restarts at zero,
so a run resumed after an entry can end late. Rerun from the start to obtain a
checkpoint that includes the count.

The checkpoint also records which spacecraft are active. A spacecraft lost to
an impact or a touchdown before the checkpoint stays inactive in the resumed
run, as it would in a run that was never interrupted: it is not propagated, its
guidance, navigation and control models treat it as inactive (`skip_inactive`),
and its impact or touchdown does not happen again. A checkpoint written by an
earlier SpaceAGORA records this only for an interlinked run; any other run
resumed from such a checkpoint starts with every spacecraft active, as before.

A run with scoped dynamic effectors (`scoped`) also records each effector's kind
and the run positions it acts on, and resumes only with the same scopes. It
refuses a checkpoint whose scopes differ, and a checkpoint written without scoped
effectors (by a run without them, or by an earlier SpaceAGORA). A run without
scoped effectors records nothing more than before.

Every checkpoint also records the spacecraft each guidance, navigation and
control model acts on, its `spacecraft_scope` as the run-start initializers left
it: every spacecraft, the run positions of a bound model (`BoundSpacecraft`, which
a model declared on a spacecraft has once bound to it), or a coupled model's
primary and partners. A resumed run refuses a checkpoint whose models act on other
spacecraft, or whose guidance, navigation or control models are another number,
before any model's run state is restored, whether or not the models keep run
state. So a momentum manager's saved command never goes to a manager of another
spacecraft. A checkpoint written by SpaceAGORA 0.5.1 or earlier holds no such
record, so a run with any guidance, navigation or control model refuses it before
anything is restored: resume it with SpaceAGORA 0.5.1, or run from the start. A
run without those models resumes from it as before.

## Campaign and multi-process runs

`run_monte_carlo` and `run_constellation_ensemble` samples each call
`run_simulation` independently (per-sample `results_directory` and
`SimulationSettings` come from whatever `SimulationConfiguration` your
per-seed closure returns), so the file layout above applies per sample
whether the outer route is serial, threaded, or process-backed. See
[Parallel Execution](parallel_execution.md) for how to select and configure
outer routing.

`run_constellation_ensemble(args; merge_results=true)` also merges its members
into one table in the long layout's columns, and with `results=true` writes
`constellation_members.feather` (`.csv` with `save_csv`) and
`constellation_summary.csv` beside the member folders. See
[One table for the whole ensemble](constellations.md#One-table-for-the-whole-ensemble).


## In-memory constellation recorder

For fixed-cadence output without retaining the solver's complete step history,
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
The recorder preallocates arrays and grows them when needed; memory still grows
with the number of spacecraft and recorded samples. Read the returned views
after the solve, and obtain new views after resetting or reusing the recorder.
`trajectory_save_data(recorder)` materializes ordinary saved snapshots when
needed; that conversion allocates dictionaries at the output boundary.

This callback leaves configured file output unchanged. Use
`SimulationSettings(results=false)` when an in-memory result is sufficient.
Explicit `save_fields` use their supplied getters, including a custom getter
with a built-in field name. Additional default fields without a specialized
array filler also retain their getters. Timing improvements depend on the run
and should be measured before adopting this recorder for a campaign.

### Backlog scheduling and stored data

`DataBacklogPolicy` scores a candidate with
`min(link_rate, held/horizon + current_production, free_space/horizon + current_outflow)`.
Receiver outflow is sampled before the new selection, from already-active links
that remain available, with the runtime's carried-rate caps. It is an estimate
for scoring; the data-flow calculation still prevents stores from overfilling.
The saved transfer-rate and stored/lost-bit columns report actual carried flows,
not this estimated score. A full store with no outflow still scores zero room.
