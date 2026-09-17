# Simulation Outputs

Use this page when you want to understand the files written by
`run_simulation`, the column names in the CSV, or how to load results for
post-processing.

This page is for users who have already run a simulation and now need to
interpret the output.

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [The Integrated State](integrated_state.md) (which columns are integrated and which are derived)
- [Simulation Configuration](simulation_configuration.md)
- [Verification Study](verification_study.md)
- [Recipes](recipes.md)

## Output files

When `simulation_settings.results = true` (the default), `run_simulation`
writes three files to the `results_directory`:

```text
output/
  simulation_results.csv         ← tabular results, one row per time step
  simulation_results.feather     ← same data in Apache Arrow / Feather format
  simulation_results.manifest.toml ← metadata: row count, spacecraft count, file hashes
```

The Feather file is always written. The CSV is written when
`simulation_settings.save_csv = true` (also the default). The manifest records
`schema_version`, `created_utc`, `mission_time_s`, `steps`,
`spacecraft_count`, and SHA-256 hashes for each data file.

## Where the examples write, and how to keep runs apart

Most repository examples build their configuration with `make_example_config`,
whose `results_directory` is `<repository root>/output` unless
`SPACEAGORA_CLI_OUTPUT_DIR` is set; the CLI's `--output-dir` sets exactly that
variable. With the default result settings, the three file names do not
change, so:

- two runs using the same results directory overwrite those files;
- `julia --project=. src/cli/main.jl run --example=<script> --output-dir=output/<name>`
  gives each run its own directory, for the scripts that take their directory
  from `make_example_config`: the first-run scripts (`AGORA_Basic_Quickstart.jl`,
  `AGORA_Earth_NoGRAM.jl`), the controls and torque
  tests, and the mission scripts `AGORA_Basic_GRAMEarth.jl`, `AGORA_Odyssey.jl`,
  `AGORA_Vex.jl`, `AGORA_Titan.jl`, `AGORA_Magellan.jl`, `AGORA_LOFTID.jl`,
  `AGORA_Mars_NoGRAM.jl` and `Earth_Thruster_Test.jl`;
- `--output-dir` has no effect on the scripts that choose their own directory:
  `AGORA_Earth.jl`, `AGORA_Keplerian.jl` and `Earth_Navigation.jl` write to
  `output/` directly, `AGORA_Earth_Aerobraking.jl` to `output/earth_aerobraking/`,
  `AGORA_Mars_RAAN_Scenario.jl` to `output/mars_raan_scenario/`, and the RPO,
  robot-arm and cloth demos to their own `output/<demo>/` directories; running
  one of them twice overwrites its previous results;
- for your own scripts, pass `results_directory=` to `make_example_config` or
  set it on `SimulationSettings`.

The quickstart example additionally saves four PNG plots under
`<results_directory>/plots/`. The RPO examples also generate HTML plots.
`Solar_Panel_Cloth_Deployment_Demo.jl` writes four HTML files in its demo
directory instead of the three simulation result files.

`AGORA_Earth_MonteCarlo.jl` prints its successful and failed sample counts and
elapsed time to the terminal. It disables result saving, so it writes no
result files, including when `--output-dir` is supplied.

Smoke mode (`--smoke` on the CLI, or `SPACEAGORA_EXAMPLE_SMOKE=1` for a script)
shortens the mission to at most 120 s and one orbit. It does not honour
`--output-dir` or `SPACEAGORA_CLI_OUTPUT_DIR`: the smoke configuration sets
`results_directory` to `output/` under the current working directory, so a
smoke run replaces the results of a previous full run in that `output/`
(run from the repository root, that is the same `output/` the quickstart
writes to). Results are kept only when `SPACEAGORA_EXAMPLE_SMOKE_RESULTS=1`,
which the CLI sets for you.

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

## CSV column reference

All per-spacecraft columns use the prefix `sc{N}_` where `N` is the
spacecraft `id` field (1-based). A single-spacecraft run produces `sc1_*`
columns.

### Time

| Column | Unit | Description |
|---|---|---|
| `time` | s | Elapsed simulation time from epoch |

### Position and velocity (inertial frame)

| Column | Unit | Description |
|---|---|---|
| `sc1_pos_1` | m | Inertial position X component |
| `sc1_pos_2` | m | Inertial position Y component |
| `sc1_pos_3` | m | Inertial position Z component |
| `sc1_vel_1` | m/s | Inertial velocity X component |
| `sc1_vel_2` | m/s | Inertial velocity Y component |
| `sc1_vel_3` | m/s | Inertial velocity Z component |

### Geodetic state

| Column | Unit | Description |
|---|---|---|
| `sc1_altitude` | m | Altitude above the reference ellipsoid |
| `sc1_latitude_deg` | deg | Geodetic latitude |
| `sc1_longitude_deg` | deg | Longitude |
| `sc1_periapsis_altitude` | m | Current osculating periapsis altitude |

### Mass

| Column | Unit | Description |
|---|---|---|
| `sc1_mass` | kg | Total spacecraft mass (structure + remaining propellant) |

### Atmosphere and aerodynamics

| Column | Unit | Description |
|---|---|---|
| `sc1_wind_1` | m/s | Atmospheric wind vector X |
| `sc1_wind_2` | m/s | Atmospheric wind vector Y |
| `sc1_wind_3` | m/s | Atmospheric wind vector Z |
| `sc1_drag_1` | N | Aerodynamic drag force X (planet frame) |
| `sc1_drag_2` | N | Aerodynamic drag force Y |
| `sc1_drag_3` | N | Aerodynamic drag force Z |
| `sc1_lift_1` | N | Aerodynamic lift force X |
| `sc1_lift_2` | N | Aerodynamic lift force Y |
| `sc1_lift_3` | N | Aerodynamic lift force Z |
| `sc1_cross_1` | N | Aerodynamic cross force X |
| `sc1_cross_2` | N | Aerodynamic cross force Y |
| `sc1_cross_3` | N | Aerodynamic cross force Z |

### Thermal

| Column | Unit | Description |
|---|---|---|
| `sc1_heat_rate` | W/cm² | Largest of the per-link stagnation heat rates (the built-in Maxwellian model returns W/cm²) |
| `sc1_heat_load` | J/cm² | Largest of the per-link accumulated heat loads (each link integrates its own rate, without an area-unit conversion); not the sum over links |

### Attitude (orientation_sim only)

These columns are present only when `mission_configuration.orientation_sim = true`:

| Column | Unit | Description |
|---|---|---|
| `sc1_q_1` | — | Attitude quaternion component `x` (scalar-last convention `[x, y, z, w]`, inertial to body) |
| `sc1_q_2` | — | Attitude quaternion component `y` |
| `sc1_q_3` | — | Attitude quaternion component `z` |
| `sc1_q_4` | — | Attitude quaternion scalar component `w`; `[0, 0, 0, 1]` is the identity attitude |

## Multi-spacecraft runs

When multiple spacecraft are passed to `DynamicsModel`, each spacecraft
contributes its own column set. Spacecraft 1 produces `sc1_*` columns,
spacecraft 2 produces `sc2_*`, and so on. The `id` field on each
`SpacecraftModel` sets the numeric suffix.

## Output rate

The `data_rate` field in `MissionConfiguration` sets the output sample cadence
in seconds of simulated time. The default is `10.0` seconds. Decreasing
`data_rate` produces more rows and larger files; increasing it reduces output
size but loses temporal resolution.

The `num_steps_to_save` field controls how many time steps are buffered in
memory before a flush. For very long runs, tuning this can reduce peak memory
use.

## Checkpointing

Enable periodic checkpointing via `SimulationSettings`:

```julia
SM.SimulationSettings(
    checkpoint_enabled      = true,
    checkpoint_interval_s   = 300.0,         # every 300 s of simulated time
    checkpoint_directory    = "output/checkpoints"
)
```

A crashed or interrupted run can be resumed by setting
`resume_from_checkpoint = true` in `SimulationSettings`. The checkpoint
directory defaults to `results_directory/checkpoints` when left empty.

## Campaign and multi-process runs

`run_monte_carlo` and `run_constellation_ensemble` samples each call
`run_simulation` independently (per-sample `results_directory` and
`SimulationSettings` come from whatever `SimulationConfiguration` your
per-seed closure returns), so the file layout above applies per sample
whether the outer route is serial, threaded, or process-backed. See
[Parallel Execution](parallel_execution.md) for how to select and configure
outer routing.


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
