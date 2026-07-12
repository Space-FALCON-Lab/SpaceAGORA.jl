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
| `sc1_heat_rate` | W/m² | Instantaneous stagnation heat rate |
| `sc1_heat_load` | J/m² | Accumulated heat load (time-integral of heat rate) |

### Attitude (orientation_sim only)

These columns are present only when `mission_configuration.orientation_sim = true`:

| Column | Unit | Description |
|---|---|---|
| `sc1_q_1` | — | Attitude quaternion component 1 (scalar-first convention) |
| `sc1_q_2` | — | Attitude quaternion component 2 |
| `sc1_q_3` | — | Attitude quaternion component 3 |
| `sc1_q_4` | — | Attitude quaternion component 4 |

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
