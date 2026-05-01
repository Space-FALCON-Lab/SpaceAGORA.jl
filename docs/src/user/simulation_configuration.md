# Building a Simulation Configuration

Use this page when you want to assemble a `SimulationConfiguration` from
scratch rather than relying on the `make_example_config` shortcut used in the
repository examples.

This page is for users who have already completed the quickstart and now want
to write their own scenario, change the spacecraft, switch atmosphere models, or
understand what each configuration field controls.

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Atmosphere Models](atmosphere_models.md)
- [Solver Configuration](solver_configuration.md)
- [Extensibility](../extensibility.md)

## The top-level struct

`SimulationConfiguration` is the single object passed to `run_simulation`. It
is composed from several nested structs, all accessed through
`SpaceAGORA.SimulationModel` (abbreviated `SM` in the examples):

```julia
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel

config = SM.SimulationConfiguration(
    file_paths             = SM.FilePaths(),
    simulation_settings    = SM.SimulationSettings(...),
    mission_configuration  = SM.MissionConfiguration(...),
    environment_model      = SM.EnvironmentModel(...),
    dynamics_model         = SM.DynamicsModel([spacecraft], effectors),
    guidance_model         = SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model       = SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model          = SM.ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time           = SM.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances = SM.IntegrationTolerances()
)

run_simulation(config)
```

## InitialTime

Specifies the simulation epoch. All fields default to the J2000 epoch
(2000-01-01 00:00:00).

```julia
SM.InitialTime(
    year   = 2024,
    month  = 5,
    day    = 27,
    hour   = 5,
    minute = 0,
    second = 0.0
)
```

This epoch is used as the reference for SPICE-backed ephemerides and for
`NRLMSISE00AtmosphereModel` when `use_space_indices=true`.

## InitialCondition

`InitialCondition` accepts two input modes. Use whichever is more natural for
your scenario. All angular inputs are in degrees; the constructor converts to
radians internally.

**Mode 1: apoapsis and periapsis radii**

```julia
ic = SM.InitialCondition(
    ra  = planet.Rp_e + 1_200e3,  # apoapsis radius from planet center, m
    rp  = planet.Rp_e + 400e3,    # periapsis radius from planet center, m
    i   = 28.5,                   # inclination, degrees
    ω   = 10.0,                   # argument of periapsis, degrees
    Ω   = 20.0,                   # right ascension of ascending node, degrees
    ν   = 0.0                     # true anomaly at epoch, degrees (default 0.0)
)
```

**Mode 2: semi-major axis and eccentricity**

```julia
ic = SM.InitialCondition(
    a   = planet.Rp_e + 800e3,
    e   = 0.001,
    i   = 28.5,
    ω   = 10.0,
    Ω   = 20.0,
    ν   = 0.0
)
```

You cannot mix both modes in the same call; providing both `ra`/`rp` and
`a`/`e` raises an error.

**Cartesian initial condition**

For non-Keplerian starts or when state is available in an inertial Cartesian
frame:

```julia
ic = SM.CartesianInitialCondition(
    pos     = [6_778_137.0, 0.0, 0.0],   # inertial position, m
    vel     = [0.0, 7784.0, 0.0],         # inertial velocity, m/s
    q       = [1.0, 0.0, 0.0, 0.0],       # unit quaternion (scalar-first)
    ang_vel = [0.0, 0.0, 0.0]             # body angular velocity, rad/s
)
```

## MissionConfiguration

Controls the termination condition, time horizon, and output cadence.

```julia
SM.MissionConfiguration(
    mission_type     = SM.MissionTime,    # SM.MissionTime or SM.MissionOrbits
    mission_time     = 3600.0 * 12.0,    # total simulated seconds (MissionTime)
    number_of_orbits = 1,                # orbit count (MissionOrbits)
    keplerian        = true,             # use Keplerian orbit mode
    orientation_sim  = false,            # include attitude dynamics
    num_steps_to_save = 1000,            # output buffer size before flush
    data_rate        = 10.0              # output sample cadence, seconds
)
```

When `mission_type = SM.MissionOrbits`, the simulation terminates after
`number_of_orbits` complete orbits. When `keplerian = true`, the integrator
uses Keplerian two-phase integration (free-flight + drag-pass). Set
`keplerian = false` to keep uniform integration throughout (useful for
continuous atmosphere scenarios or full entry trajectories).

## EnvironmentModel

Wires together the planet, atmosphere, ephemerides, and thermal model. All
four of these must be provided; there are no defaults for `planet` or
`density_model`.

```julia
planet = SM.make_no_gram_planet(:earth)

SM.EnvironmentModel(
    planet            = planet,
    EI                = 120.0,                   # entry interface altitude, km
    density_model     = SM.NoAtmosphereModel(),
    ephemerides_model = SM.SimpleEphemeridesModel(),
    thermal_model     = SM.MaxwellianHeat(
                            thermal_accomodation_factor=1.0,
                            planet=planet
                        ),
    topography        = false,
    wind              = false
)
```

`EI` (entry interface) controls the altitude threshold below which atmospheric
drag is applied. For orbit-only scenarios, set it high enough that the
spacecraft never crosses it (e.g., `120.0` km for Earth).

Set `wind = true` to request wind vectors from the atmosphere model; note that
the open-data models (`NoAtmosphereModel`, `ExponentialAtmosphereModel`,
`PiecewiseExponentialAtmosphereModel`) always return zero wind regardless.

For the supported atmosphere models and their constructors, see
[Atmosphere Models](atmosphere_models.md).

## DynamicsModel

Holds the list of spacecraft and the tuple of force/torque effectors:

```julia
SM.DynamicsModel([spacecraft], (SM.InverseSquaredJ2GravityModel(),))
```

Multiple spacecraft can be passed in the array for parallel propagation. The
effectors tuple accepts any combination of `AbstractForceTorqueModel`
implementations. Common built-in effectors:

- `SM.InverseSquaredGravityModel()` — inverse-square point-mass gravity
- `SM.InverseSquaredJ2GravityModel()` — point-mass gravity with J2 oblateness

## SimulationSettings

Controls output behavior and diagnostics:

```julia
SM.SimulationSettings(
    results            = true,           # write output files
    verbose            = false,          # print solver diagnostics
    results_directory  = "output",       # directory for CSV and bundle output
    generate_plots     = false,          # reserved; disable for CLI runs
    save_csv           = true,           # write CSV alongside the Feather bundle
    checkpoint_enabled = false,          # periodic checkpoint for restart safety
    checkpoint_interval_s = 300.0,      # checkpoint cadence, simulated seconds
    checkpoint_directory  = "",         # defaults to results_directory/checkpoints
    resume_from_checkpoint = false       # resume from latest checkpoint if present
)
```

Set `results = false` to run without writing any output (useful for
performance profiling or validation-only runs).

## FilePaths

`FilePaths` holds paths to licensed external asset directories. For no-GRAM
runs, the defaults are fine and this struct does not need to be set explicitly:

```julia
SM.FilePaths(
    results              = "Results",
    GRAM                 = "data/GRAMSuite.jl/GRAM Suite 2.0",
    SPICE                = "data/GRAMSuite.jl/GRAM Suite 2.0/SPICE",
    topography_harmonics = "data/Topography_harmonics_data",
    gravity_harmonics    = "data/Gravity_harmonics_data"
)
```

## Using `make_example_config`

For quick studies and all repository examples, `make_example_config` from
`SpaceAGORA.TelemetryVerification` assembles the configuration in one call:

```julia
import SpaceAGORA.TelemetryVerification: make_example_config, make_three_body_spacecraft, run_and_report

spacecraft = make_three_body_spacecraft(
    bus_dims        = (2.05, 2.05, 2.8),
    panel_dims      = (0.01, 2.85, 1.0),
    bus_mass        = 620.0,
    panel_mass_each = 10.0,
    panel_offset_y  = 2.05/2.0 + 2.85/2.0,
    ic              = SM.InitialCondition(ra=..., rp=..., i=28.5, ω=10.0, Ω=20.0, ν=0.0),
    prop_mass       = 200.0,
    id              = 1
)

config = make_example_config(
    planet             = SM.make_no_gram_planet(:earth),
    spacecraft         = spacecraft,
    mission_time       = 3600.0 * 12.0,
    initial_time       = SM.InitialTime(year=2024, month=1, day=1),
    dynamic_effectors  = (SM.InverseSquaredJ2GravityModel(),),
    density_model      = SM.NoAtmosphereModel(),
    ephemerides_model  = SM.SimpleEphemeridesModel(),
    orientation_sim    = false,
    keplerian          = true,
    EI_km              = 120.0,
    verbose            = true
)

run_and_report(config)
```

`make_example_config` is not part of the stable root `SpaceAGORA` export — it
lives in `SpaceAGORA.TelemetryVerification` and is imported by all repository
examples through `examples/common.jl`. It is appropriate for examples and
quick studies; for production use, build `SimulationConfiguration` directly.

## `make_three_body_spacecraft`

Constructs a three-body spacecraft: a main bus plus two symmetric solar
panels. This is the geometry used by all repository examples.

```julia
make_three_body_spacecraft(
    bus_dims        = (x, y, z),         # bus bounding box, m
    panel_dims      = (t, span, chord),  # panel thickness, half-span, chord, m
    bus_mass        = 620.0,             # kg
    panel_mass_each = 10.0,              # kg per panel
    panel_offset_y  = offset,            # panel center offset from bus center, m
    ic              = SM.InitialCondition(...),
    prop_mass       = 200.0,             # propellant mass, kg
    id              = 1                  # spacecraft ID, used in output column names
)
```

The `id` field determines the column prefix in the output CSV: spacecraft 1
gets `sc1_*` columns, spacecraft 2 gets `sc2_*`, and so on.
