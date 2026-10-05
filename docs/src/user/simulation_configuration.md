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
- [Adding a Force or Torque of Your Own](custom_effector.md)
- [Stopping a Simulation on a Condition](stop_conditions.md)
- [Extensibility](../extensibility.md)

## The top-level struct

`SimulationConfiguration` is the single object passed to `run_simulation`. It
is composed from several nested structs. The common setup types
(`SimulationConfiguration`, `MissionConfiguration`, `EnvironmentModel`,
`DynamicsModel`, `SpacecraftModel`, `Link`, `InitialTime`,
`IntegrationTolerances`, the inverse-square gravity effectors, and the
`make_example_config` helper) are exported from the root module, so
`using SpaceAGORA` alone is enough. Less common types are reached through
`SpaceAGORA.SimulationModel` (abbreviated `SM` below).

```julia
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel

config = SimulationConfiguration(
    file_paths             = SM.FilePaths(),
    simulation_settings    = SimulationSettings(...),
    mission_configuration  = MissionConfiguration(...),
    environment_model      = EnvironmentModel(...),
    dynamics_model         = DynamicsModel([spacecraft], effectors),
    guidance_model         = GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model       = NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model          = ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time           = InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances = IntegrationTolerances()
)

run_simulation(config)
```

The guidance, navigation and control models are optional. Leave them out and
each defaults to an empty model, the same as the zero-argument constructors
`SM.GuidanceModel()`, `SM.NavigationModel()` and `SM.ControlModel()` (no
effectors, no rates). The minimal form is:

```julia
config = SM.SimulationConfiguration(
    simulation_settings   = SM.SimulationSettings(...),
    mission_configuration = SM.MissionConfiguration(...),
    environment_model     = SM.EnvironmentModel(...),
    dynamics_model        = SM.DynamicsModel([spacecraft], effectors),
    initial_time          = SM.InitialTime(year=2024, month=1, day=1),
)
```

`environment_model`, `dynamics_model` and `initial_time` stay required;
pass explicit models only when you attach guidance, navigation or control
effectors.

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

`InitialCondition` accepts orbital elements, apsis radii, or apsis altitudes.
The keyword-only constructor supports Modes 1 and 2; passing a planet as the
first positional argument selects Mode 3. All angular inputs are in degrees;
the constructor converts to radians internally.

**Mode 1: apoapsis and periapsis radii**

```julia
ic = SM.InitialCondition(
    ra  = planet.Rp_e + 1_200e3,  # apoapsis radius from planet center, m
    rp  = planet.Rp_e + 400e3,    # periapsis radius from planet center, m
    i   = 28.5,                   # inclination, degrees
    ω   = 10.0,                   # argument of periapsis, degrees
    Ω   = 20.0,                   # right ascension of ascending node, degrees
    ν   = 0.0                     # start at periapsis; omitted ν defaults to 180.0 (apoapsis)
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

For the `a`/`e` form, omitting `ν` starts at periapsis (`0.0` degrees).
For the `ra`/`rp` form, omitting `ν` starts at apoapsis (`180.0` degrees), so
set `ν=0.0` explicitly when you want to start at periapsis.

Choose one pair of inputs per call. The current keyword-only constructor
gives `ra`/`rp` precedence if `a`/`e` are also supplied; it does not reject
that combination. Supplying only one of `ra` and `rp` raises an error.

**Mode 3: apoapsis and periapsis altitudes above the ellipsoid**

```julia
planet = SM.make_no_gram_planet(:earth)
initial_time = SM.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0)
ephemerides_model = SM.SimpleEphemeridesModel()

ic = SM.InitialCondition(
    planet;
    ra = 1_200e3,  # apoapsis altitude above the reference ellipsoid, m
    hp = 400e3,    # periapsis altitude above the reference ellipsoid, m
    i = 28.5,
    ω = 10.0,
    Ω = 20.0,
    ν = 0.0,      # start at periapsis; omitted ν defaults to 180.0 (apoapsis)
    initial_time=initial_time,
    ephemerides_model=ephemerides_model,
)
```

In this form, both `ra` and `hp` are altitudes in metres. The positional
`planet` argument selects this meaning of `ra`; in Mode 1, `ra` is a radius
from the planet's centre. The constructor finds the radius at each apsis
whose geodetic altitude matches the requested value, then computes `a` and
`e`. Both altitudes must be nonnegative, and the resulting apoapsis radius
must exceed the periapsis radius.

Use the same `initial_time` and `ephemerides_model` in the simulation
configuration so the constructor uses the intended inertial-to-planet-fixed
frame. Advanced callers can supply that rotation directly as `L_PI`, which
takes precedence. With neither an explicit rotation nor an initial time,
the constructor uses an initialized `planet.L_PI`, or the identity rotation
if that matrix is unavailable or all zeros.

**Cartesian initial condition**

For non-Keplerian starts or when state is available in an inertial Cartesian
frame:

```julia
using StaticArrays: SVector

ic = SM.CartesianInitialCondition(
    [6_778_137.0, 0.0, 0.0],              # inertial position, m (positional)
    [0.0, 7784.0, 0.0];                   # inertial velocity, m/s (positional)
    q       = SVector(0.0, 0.0, 0.0, 1.0),   # unit quaternion, scalar-last [x, y, z, w]; the identity attitude
    ang_vel = SVector(0.0, 0.0, 0.0)         # bus angular velocity, rad/s
)
```

Position and velocity are positional arguments; `q` and `ang_vel` are keyword
arguments that must be `SVector`s and default to the values shown, so
`SM.CartesianInitialCondition(pos, vel)` is also valid.

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

Wires together the planet, atmosphere, ephemerides, and thermal model. Only
`planet`, `EI`, `density_model`, and `thermal_model` are required; the rest
have defaults: `ephemerides_model` defaults to `SM.SpiceEphemeridesModel()`,
`wind` defaults to `true`, and `topo_degree`/`topo_order` default to `90`.

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

`thermal_model` is one of three heating models, each evaluated per thermal
link (panel) and integrated into that link's heat-load state:

- `SM.MaxwellianHeat(thermal_accomodation_factor, planet)`: free-molecular
  Maxwellian heat flux, for rarefied aerobraking corridors.
- `SM.SuttonGravesHeat(planet=planet, nose_radius_m=0.5)`: Sutton–Graves
  stagnation-point convective heating, `q = k √(ρ/r_n) v³`, with `k` defaulting
  to the planet's coefficient (`planet.k`).
- `SM.TabularHeat("aerothermal.csv")`: a vehicle-level flux interpolated from an
  aerothermal database tabulated on a velocity × density grid (CSV columns
  `velocity_m_s`, `density_kg_m3`, `heat_rate_W_cm2`), or
  `SM.TabularHeat(velocities, densities, heat_rates)` from arrays.

`EI` (entry interface) is the altitude at which the integrator switches
between its orbit and atmosphere step-size/tolerance regimes (see
`IntegrationTolerances`). It is not a force gate: whenever a non-vacuum
`density_model` is configured together with an aerodynamic effector, drag is
evaluated at every altitude and the density model itself decides where the
force becomes negligible.

Set `wind = true` to request wind vectors from the atmosphere model; note that
the open-data models (`NoAtmosphereModel`, `ExponentialAtmosphereModel`,
`PiecewiseExponentialAtmosphereModel`) always return zero wind regardless.

With `wind = false` the simulation treats the atmosphere as co-rotating with the
planet for every model: density queries are made with `wind=false`, and the
wind used for the atmosphere-relative velocity in aerodynamics and guidance, and
recorded in the `wind` output, is zero. This matters for native GRAM, whose own
`wind=false` query still returns its nominal (mean) winds, and for GRAM grid
snapshots, which return stored winds either way; the simulation discards those
values. Density and temperature are unaffected. A density-only run also stays
eligible for the automatic native-GRAM pools, which otherwise keep perturbed
wind requests on the locked path (see [Parallel Execution](parallel_execution.md)).

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
    generate_plots     = true,           # generate plots after simulation (default)
    generate_filenames = false,          # embed run parameters in output filenames (default)
    normalize          = false,          # legacy compatibility flag; typed run_simulation propagates SI state directly (default)
    save_csv           = true,           # write CSV alongside the Feather bundle
    save_visualization_scene = false,    # write the viewer scene sidecar and link_pose columns (see Simulation Outputs)
    checkpoint_enabled = false,          # periodic checkpoint for restart safety
    checkpoint_interval_s = 300.0,      # checkpoint cadence, simulated seconds
    checkpoint_directory  = "",         # defaults to results_directory/checkpoints
    resume_from_checkpoint = false       # resume from latest checkpoint if present
)
```

Set `results = false` to run without writing any output (useful for
performance profiling or validation-only runs). Set `generate_plots = false`
to skip plot generation, which is typically what you want for CLI/batch runs
and performance studies.

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

## Panel angles and heating

With a built-in aerodynamic model, heating uses the current panel geometry. If
`orientation_sim = true`, it follows the propagated spacecraft attitude and each
panel's orientation. Otherwise it follows the aerodynamic model's fixed-attitude
incidence policy. Panel-control changes take effect in heating without requiring
a force calculation first. The aerodynamic scale factor is not applied directly
to heat rates; it can still change heating by changing the trajectory.

`AerodynamicCoefficientfM.fixed_attitude_incidence` governs both aerodynamics
and heating when attitude is not propagated. `AerodynamicCoefficientConstant`
and `AerodynamicCoefficientNoBallisticFlight` always use `:max_drag` in that
case. With propagated attitude, incidence follows the wind-relative flow;
without it, the selected fixed-attitude policy supplies the geometric angle.

A configuration without a recognized built-in aerodynamic model keeps its
existing `Link.α` heating input. A custom effector does not become a geometric
owner by declaring `environment_requirements(model).atmosphere = true`; its
heating input remains the angle maintained by that effector or controller.
If a built-in model is also present, it determines geometric incidence.
Multiple built-in models must use the same fixed-attitude policy when attitude
is not propagated; conflicting policies raise an error at the first thermal
sample because there is no single angle for heating to use.

For a controlled panel, heating follows its executed orientation. The resulting
incidence equals the panel's commanded angle for the supported geometry with a
panel offset along the body y axis and a flow-aligned root. Other offsets or
root attitudes can give a different incidence without changing the stored
command. In the Odyssey energy-depletion example, saved maximum-link heat
columns include the uncontrolled bus, while the controller's panel heat limits
apply only to its controlled panels.

## Source ownership for contributors

The existing configuration API is assembled by `src/simulation/config/configuration.jl`.
This is a source-file organization; users still access the same types through
`SpaceAGORA.SimulationModel`. No new configuration wrapper is required.

| Source file under `src/simulation/config/` | Responsibility |
| --- | --- |
| `run_settings.jl` | Epoch, mission duration/orbits, sampling, paths, output and checkpoint settings |
| `solver_settings.jl` | `SolverConfig` and `IntegrationTolerances` |
| `environment_settings.jl` | Select and compose environmental models |
| `constellation_configuration.jl` | Existing `DynamicsModel`: member spacecraft and selected dynamic effectors |
| `simulation_configuration.jl` | Final container and `_with_configuration` helper |

`constellation_configuration.jl` is included inside the existing `SpacecraftModels`
module to preserve the identity of `DynamicsModel`; the other definitions remain
inside `SimConfig`. A one-spacecraft run and a constellation use the same collection.
The constructor retains the supplied spacecraft vector. It does not generate orbital
layouts, schedule activities or introduce additional collection validation.

Execution policy remains in `src/simulation/engine/config/`: `SimulationEngineConfig`
composes parallel, solver, runtime-policy and artifact settings. It uses the same
`SolverConfig` definition, not a second solver type. Output and checkpoint path
derivation stays in `src/io/config/`; solver environment settings are parsed in
`src/simulation/engine/adapters/from_env.jl`. These are distinct responsibilities
from assembling a scenario.

`_with_configuration` makes a shallow update and preserves unspecified references.
Runtime state isolation remains the engine's responsibility. Moving the definitions
changes neither those semantics nor typed-solver precedence over environment settings.
