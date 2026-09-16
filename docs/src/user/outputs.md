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

## Choosing what gets saved

`default_save_fields(args)` builds the column set above, and `run_simulation`
uses it whenever `save_fields` is not passed. To add more columns, pass the same
list with an `extra` argument naming the built-in fields you want:

```julia
run_simulation(args; save_fields=default_save_fields(args; extra=(:orbital_elements, :gravity_accel)))
```

`extra` accepts built-in names, `SaveField`s of your own, or any mix of the two,
and a field the defaults already cover is skipped rather than duplicated. The
example scripts take the same list through `run_and_report(args; save_fields=...)`.

`available_save_fields()` returns every name that can be requested, and
`save_field(name, args)` builds a single one if you would rather assemble the
list yourself. Passing a name that is not built in throws and lists what is.

### Built-in fields

| Name | Columns | Unit | What it holds |
|---|---|---|---|
| `:orbital_elements` | `sc{i}_orbital_elements_1..6` | mixed | Osculating classical elements |
| `:gravity_accel` | `sc{i}_gravity_accel_1..3` | m/s² | Inertial acceleration from the gravity effectors |
| `:aero_accel` | `sc{i}_aero_accel_1..3` | m/s² | Inertial acceleration from the aerodynamic effectors |
| `:total_accel` | `sc{i}_total_accel_1..3` | m/s² | Inertial acceleration from every dynamic effector |
| `:quaternion` | `sc{i}_q_1..4` | — | Attitude; already a default when `orientation_sim = true` |

A name requested from a run that cannot supply it fails at the first saved
sample rather than writing an empty column. `:quaternion` without
`orientation_sim` is the usual case; `:total_accel` is the other, when the run
carries a dynamic effector that defines no `wrench` method and the column would
otherwise under-report without saying so.

### Orbital elements

| Column | Unit | Description |
|---|---|---|
| `sc1_orbital_elements_1` | m | Semimajor axis |
| `sc1_orbital_elements_2` | — | Eccentricity |
| `sc1_orbital_elements_3` | deg | Inclination |
| `sc1_orbital_elements_4` | deg | Right ascension of the ascending node |
| `sc1_orbital_elements_5` | deg | Argument of periapsis |
| `sc1_orbital_elements_6` | deg | True anomaly |

Osculating, computed from the same inertial position and velocity the run
integrates. Circular and equatorial states fall back to argument of latitude
and true longitude in the usual way, so the angles stay defined as eccentricity
or inclination approaches zero. Note that at very small eccentricity the
periapsis direction is nearly undefined and genuinely does spin fast, so `ω`
and `ν` sweep several times per orbit while their sum, the argument of latitude,
advances once — that is the orbit, not a numerical artifact.

### Accelerations

The three acceleration fields share one machine and differ only in which of the
run's dynamic effectors they sum:

| Field | Effectors summed |
|---|---|
| `:gravity_accel` | Constant, inverse-square, J2, spherical harmonics, third-body |
| `:aero_accel` | The aerodynamic coefficient models; zero without an atmosphere |
| `:total_accel` | Every dynamic effector in the run |

All three are inertial, in m/s², written as `sc{i}_<field>_1..3`, and none
include the control model's contribution — thrust is not in `:total_accel`.

They are re-evaluated at each saved sample through the same `wrench` hooks the
right-hand side calls, rather than read from a right-hand-side buffer. This is
deliberate. The right-hand side runs every solver stage of every step, accepted
or rejected, so caching would pay on the order of 10⁵ stores to serve 10³
samples; and the cached value would be the one the last stage happened to
leave, at a different time and state than the row it landed in. Recomputing
costs one extra evaluation per sample and puts the value that belongs with the
row in the row.

### Writing your own

A `SaveField` is a name plus a getter:

```julia
SaveField(name, getter; per_satellite=false, column_prefix=String(name))
```

`getter(u, t, integrator)` is called at each saved sample. Return one entry per
spacecraft when `per_satellite` is true — written as `sc{i}_{column_prefix}`,
with `_1.._n` appended when the entry is a vector — or a single value for the
whole run otherwise. Names must be unique within a run's save set.

```julia
speed = SaveField(
    :speed,
    (u, t, integrator) -> [norm(SpaceAGORA.SimulationEngine._state_velocity_ii(u, i))
                           for i in eachindex(integrator.p.args.dynamics_model.spacecraft)];
    per_satellite=true
)
run_simulation(args; save_fields=default_save_fields(args; extra=(:orbital_elements, speed)))
```

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
