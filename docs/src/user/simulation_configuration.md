# Configuring a Simulation

This page shows how to assemble a `SimulationConfiguration` yourself, instead
of through the `make_example_config` shortcut that the repository examples
use, and what each part of the configuration controls. Read it after the
quickstart, when you want a scenario of your own, a different spacecraft or
another atmosphere model.

The quickstart example builds its configuration with `make_example_config`;
running it is the quickest way to see a complete configuration work:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Atmosphere Models](atmosphere_models.md)
- [Solver Configuration](solver_configuration.md)
- [Adding a Force or Torque of Your Own](custom_effector.md)
- [Stopping a Simulation on a Condition](stop_conditions.md)
- [Extensibility](../extensibility.md)

## The components at a glance

![How a configuration nests: make_example_config or a MissionPlan builder produces a SimulationConfiguration. It holds InitialTime, MissionConfiguration, SimulationSettings and FilePaths; the EnvironmentModel with planet, atmosphere, ephemerides and thermal model; the guidance, navigation and control models; the DynamicsModel with its force and torque models and one SpacecraftModel per vehicle, each with an initial condition, links joined by joints and hardware mounted on the links; and the solver settings. run_simulation or run_mission turns it into results.](../assets/figures/configuration_anatomy.svg)

The table lists the same components with the section that describes each.

| Component | What it sets | Section |
| --- | --- | --- |
| `InitialTime` | the simulation epoch | [InitialTime](#InitialTime) |
| `InitialCondition` | each spacecraft's starting orbit, from orbital elements, apsis radii or apsis altitudes | [InitialCondition](#InitialCondition) |
| `MissionConfiguration` | when the run ends (a time or an orbit count) and the output cadence | [MissionConfiguration](#MissionConfiguration) |
| `EnvironmentModel` | the planet, atmosphere, winds, ephemerides and thermal model | [EnvironmentModel](#EnvironmentModel) |
| `DynamicsModel` | the spacecraft and their force and torque models, gravity included | [DynamicsModel](#DynamicsModel), [Gravity and Forces](gravity_forces.md) |
| `SimulationSettings` | result files, results directory, checkpoints and verbosity | [SimulationSettings](#SimulationSettings) |
| `FilePaths` | paths to external asset directories, which the run does not read | [FilePaths](#FilePaths) |
| `make_example_config` | a complete configuration from a planet, a spacecraft, a duration and an epoch | [Using `make_example_config`](#Using-make_example_config) |

Panel angles have their own section further down this page. Articulated and
compliant vehicles and interlinks have pages of their own:
[Articulated and Compliant Vehicles](articulated_compliant.md) and
[Interlinks](interlinks.md).

## The top-level struct

`SimulationConfiguration` is the one object passed to `run_simulation`, and it
is built from the nested structs this page describes. The common setup types
(`SimulationConfiguration`, `SimulationSettings`, `MissionConfiguration`,
`EnvironmentModel`, `DynamicsModel`, the guidance, navigation and control
models, `SpacecraftModel`, `Link`, `InitialCondition`, `InitialTime`,
`IntegrationTolerances`, `SolverConfig`, the gravity models and the
`make_example_config` helper) are exported by SpaceAGORA, so
`using SpaceAGORA` is enough for them. Less common types, such as `FilePaths`,
`SpiceEphemeridesModel` and the heating models other than `MaxwellianHeat`,
are reached through `SpaceAGORA.SimulationModel`, abbreviated `SM` on this
page.

In outline, with `...` standing for the arguments the sections below describe
and `spacecraft` and `effectors` as in [DynamicsModel](#DynamicsModel):

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

The guidance, navigation and control models are optional. Left out, each
defaults to the empty model, the same as the zero-argument constructors
`SM.GuidanceModel()`, `SM.NavigationModel()` and `SM.ControlModel()` (no
effectors, no rates). `file_paths`, `simulation_settings`,
`mission_configuration`, `integration_tolerances` and `solver_config` have
defaults too, while `environment_model`, `dynamics_model` and `initial_time`
are required. A shorter form, in the same outline:

```julia
config = SM.SimulationConfiguration(
    simulation_settings   = SM.SimulationSettings(...),
    mission_configuration = SM.MissionConfiguration(...),
    environment_model     = SM.EnvironmentModel(...),
    dynamics_model        = SM.DynamicsModel([spacecraft], effectors),
    initial_time          = SM.InitialTime(year=2024, month=1, day=1),
)
```

Pass explicit GNC models only when you attach guidance, navigation or control
effectors.

## InitialTime

The simulation epoch, in UTC. Every field has a default, and the defaults give
2000-01-01 00:00:00:

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

The run measures time from this epoch. It anchors the ephemerides model, SPICE
or analytic; the date and time at which `NRLMSISE00AtmosphereModel` is
evaluated, including the space-weather indices it looks up when
`use_space_indices=true`; and the construction epoch of a native GRAM model.

## InitialCondition

`InitialCondition` takes orbital elements, apsis radii or apsis altitudes. Its
keyword-only form covers Modes 1 and 2 below; passing a planet as the first
positional argument selects Mode 3. Angles are in degrees and are converted to
radians internally. Every form also accepts the mean anomaly `M` (degrees,
elliptic orbits only) in place of the true anomaly `ν`. The examples use the
no-GRAM Earth as `planet`.

**Mode 1: apoapsis and periapsis radii**

```julia
planet = SM.make_no_gram_planet(:earth)

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

With `a` and `e`, an omitted `ν` starts the orbit at periapsis (`0.0`
degrees). With `ra` and `rp`, it starts at apoapsis (`180.0` degrees), so set
`ν=0.0` explicitly to start at periapsis.

Give one pair of inputs per call. The keyword-only constructor does not reject
`ra` and `rp` given together with `a` and `e`: it uses `ra` and `rp` and
ignores the others. Giving only one of `ra` and `rp` raises an error.

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

In this form `ra` and `hp` are both altitudes in meters. The positional
`planet` argument is what gives `ra` this meaning; in Mode 1, `ra` is a radius
from the planet's center. The constructor finds the radius at each apsis whose
geodetic altitude matches the requested value, then computes `a` and `e`. Both
altitudes must be nonnegative, and the resulting apoapsis radius must exceed
the periapsis radius. This mode therefore cannot express a circular orbit:
equal apoapsis and periapsis altitudes give equal radii, which are refused
(rounding occasionally lets an orbit with eccentricity near zero through). For
a circular orbit, use Mode 1 with `ra` equal to `rp`, or Mode 2 with
`e = 0.0`.

The apsis directions depend on the rotation from the inertial to the
planet-fixed frame, so pass the same `initial_time` and `ephemerides_model` as
the simulation configuration. Given `initial_time` without `ephemerides_model`,
the constructor uses `SM.SpiceEphemeridesModel()`, which needs SPICE kernels.
Advanced callers can supply the rotation directly as `L_PI`, which takes
precedence. With neither `L_PI` nor `initial_time`, the constructor uses
`planet.L_PI` if it has been initialized, and the identity rotation if that
matrix is missing or all zeros.

**Cartesian initial condition**

For a start that is not given as orbital elements, or a state you already have
in inertial (J2000) Cartesian coordinates:

```julia
ic = SM.CartesianInitialCondition(
    [6_778_137.0, 0.0, 0.0],              # inertial position, m (positional)
    [0.0, 7784.0, 0.0];                   # inertial velocity, m/s (positional)
    q       = [0.0, 0.0, 0.0, 1.0],       # unit quaternion, scalar-last [x, y, z, w]; the identity attitude
    ang_vel = [0.0, 0.0, 0.0]             # bus angular velocity, rad/s
)
```

Position and velocity are positional arguments. `q` and `ang_vel` are keywords
that take any vector or tuple of 4 or 3 numbers (an `SVector` too) and default
to the values shown, so `SM.CartesianInitialCondition(pos, vel)` is also
valid. `InitialCondition` takes `q` and `ang_vel` the same way.

## MissionConfiguration

`MissionConfiguration` sets when the run ends, whether attitude is propagated
and how often results are saved:

```julia
SM.MissionConfiguration(
    mission_type     = SM.MissionTime,    # SM.MissionTime or SM.MissionOrbits
    mission_time     = 3600.0 * 12.0,    # simulated seconds; the run never goes past it
    number_of_orbits = 1,                # completed orbits (MissionOrbits)
    keplerian        = true,             # native GRAM only: see below
    orientation_sim  = false,            # propagate attitude
    num_steps_to_save = 1000,            # not read by the run
    data_rate        = 10.0              # interval between saved rows, seconds
)
```

Without arguments the fields default to `MissionTime`, a `mission_time` of
1,080,000 s (12.5 days), one orbit, `keplerian = true`, no attitude,
`num_steps_to_save = 1000` and `data_rate = 10.0`. `mission_time`,
`number_of_orbits`, `num_steps_to_save` and `data_rate` must be positive.

With `mission_type = SM.MissionOrbits`, the run ends once every active
spacecraft has completed `number_of_orbits` orbits; an orbit is counted each
time a spacecraft passes apoapsis. `mission_time` still bounds the run in this
mode, so set it long enough for the orbits you ask for.

`keplerian` matters only with the native GRAM atmosphere
(`GRAMAtmosphereModel`, and the `GRAMAtmosphereModelSurrogate` wrapper around
it). With `keplerian = true`, the default, GRAM is queried at every altitude
up to 2000 km. With `keplerian = false`, densities above the entry interface
`EI` come from a polynomial fit for the planet instead of from GRAM. Above
2000 km these models return vacuum either way. The flag does not change how the
integrator steps: the switch between the orbit and atmosphere step limits and
tolerances follows `EI` (see [EnvironmentModel](#EnvironmentModel) and
[Solver Configuration](solver_configuration.md)).

To size `mission_time` from an orbit count, use the planet's gravitational
parameter `planet.μ` in m³/s², a field of every built-in planet including the
one returned by `make_no_gram_planet(:earth)`. One Keplerian period is
`2π * sqrt(a^3 / planet.μ)` seconds for a semi-major axis `a` in meters, so a
400 km circular Earth orbit (`a = planet.Rp_e + 400e3`) takes about
92.6 minutes. `make_example_config` always builds a `MissionTime` mission, so
pass the period times the number of orbits as its `mission_time`; to stop
after a number of completed orbits instead, construct `MissionConfiguration`
yourself with `mission_type = SM.MissionOrbits`.

## EnvironmentModel

`EnvironmentModel` brings together the planet, atmosphere, ephemerides and
thermal model. Only `planet`, `EI`, `density_model` and `thermal_model` are
required. The others default to `ephemerides_model = SM.SpiceEphemeridesModel()`,
`wind = true`, `topography = false`, and `topo_degree` and `topo_order` of
`90`. `EI`, `topo_degree` and `topo_order` must be nonnegative.

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

`thermal_model` is one of three heating models. Each gives a heat rate that is
integrated into the heat-load state of every thermal link (panel):

- `SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet)`:
  free-molecular Maxwellian heat flux, for rarefied aerobraking corridors,
  evaluated with each link's incidence.
- `SM.SuttonGravesHeat(planet=planet, nose_radius_m=0.5)`: Sutton–Graves
  stagnation-point convective heating, `q = k √(ρ/r_n) v³`. `nose_radius_m`
  is required, and `k` defaults to the planet's coefficient (`planet.k`).
- `SM.TabularHeat("aerothermal.csv")`: a vehicle-level flux interpolated from
  an aerothermal database tabulated on a velocity × density grid (CSV columns
  `velocity_m_s`, `density_kg_m3`, `heat_rate_W_cm2`), or
  `SM.TabularHeat(velocities, densities, heat_rates)` from arrays.

`SuttonGravesHeat` and `TabularHeat` apply the same rate to every link.

`EI`, the entry interface altitude in km, is where the integrator switches
between its orbit and atmosphere step limits and tolerances (see
`IntegrationTolerances` on [Solver Configuration](solver_configuration.md)).
The boundary is a sphere `EI` km above the planet's equatorial radius. It does
not switch forces on or off: whenever a density model other than
`NoAtmosphereModel()` is configured together with an aerodynamic effector,
drag is evaluated at every altitude, and the density model itself decides
where the force becomes negligible.

Set `wind = true` to use the wind vectors the atmosphere model returns. The
open-data models (`NoAtmosphereModel`, `ExponentialAtmosphereModel`,
`PiecewiseExponentialAtmosphereModel`) always return zero wind.

With `wind = false`, the simulation treats the atmosphere as rotating with the
planet, whatever the model: it queries density with `wind=false`, and the wind
it uses for the atmosphere-relative velocity in aerodynamics and guidance, and
saves in the `wind` output, is zero. This matters for native GRAM, whose own
`wind=false` query still returns its nominal (mean) winds, and for GRAM grid
snapshots, which return stored winds either way; the simulation discards those
values. Density and temperature are unaffected. A run without winds also stays
eligible for the automatic native-GRAM pools, which keep runs that request
perturbed winds on the locked path (see
[Parallel Execution](parallel_execution.md)).

For the supported atmosphere models and their constructors, see
[Atmosphere Models](atmosphere_models.md).

### Ephemerides models

`SM.SpiceEphemeridesModel()` takes the planet's orientation and the Sun's
direction from the run's SPICE kernels. `SimpleEphemeridesModel()` needs no
kernels: it turns the planet about its spin axis (Earth by its Greenwich mean
sidereal time) and gives the Sun's direction for Earth only. An ephemerides
model of your own subtypes `SpaceAGORA.AbstractEphemeridesModel` and
implements five hooks, each taking the model first:

| Hook | Returns |
|---|---|
| `ephemerides_epoch_seconds(model, initial_time)` | the run's start, in seconds past J2000 on the model's timeline |
| `ephemerides_frame_lpi(model, planet, et)` | the rotation `L` from the inertial (J2000) frame to the planet-fixed frame, `r_pf = L * r_ii` |
| `ephemerides_sun_ii(model, planet, et)` | the unit vector from the planet's center to the Sun, inertial, or `nothing` |
| `ephemerides_requires_spice(model)` | whether the model reads SPICE kernels |
| `ephemerides_cache_key(model)` | a tuple; models with equal keys must give the same frames and times |

`et` is the run's start plus the elapsed time, in seconds. The two built-in
models, in `src/models/environment/ephemerides/simple_ephemerides.jl`,
implement every hook. The
[Extension Hook Reference](../generated/extension_hook_reference.md) lists
when the engine calls them. Code written for the older names, which took the
model last (`planet_frame_lpi` and two others, removed in SpaceAGORA 0.5.0),
is converted as
[Migrating Extension Hooks](../reference/migrating_extension_hooks.md#Removed-hook-names)
shows.

## DynamicsModel

`DynamicsModel` holds the vector of spacecraft and the tuple of force and
torque models (effectors). For a `spacecraft` built as in
[Using `make_example_config`](#Using-make_example_config):

```julia
SM.DynamicsModel([spacecraft], (SM.InverseSquaredJ2GravityModel(),))
```

Pass several spacecraft in the vector to propagate them together in one run.
The effector tuple takes any combination of `AbstractForceTorqueModel`
implementations, and an entry `scoped(effector, scope)` acts only on the
spacecraft its scope names. The built-in gravity models, spherical-harmonic
gravity among them, and how force models combine are on
[Gravity and Forces](gravity_forces.md).

## SimulationSettings

`SimulationSettings` controls the result files, checkpoints and diagnostics.
All its fields, with their defaults except where a comment says otherwise:

```julia
SM.SimulationSettings(
    results            = true,           # write result files
    verbose            = false,          # print run diagnostics
    results_directory  = "output",       # directory for the CSV and Feather results
    generate_plots     = true,           # not read by run_simulation (see below)
    generate_filenames = false,          # not read by run_simulation
    normalize          = false,          # keep false: the run works in SI units and refuses true
    save_csv           = true,           # write CSV alongside the Feather bundle
    save_visualization_scene = false,    # write the viewer scene sidecar and link_pose columns (see Simulation Outputs)
    checkpoint_enabled = false,          # write periodic checkpoints for restarts
    checkpoint_interval_s = 300.0,       # checkpoint interval, simulated seconds
    checkpoint_directory  = "",          # "" means results_directory/checkpoints
    resume_from_checkpoint = false,      # resume from the latest checkpoint if one exists
    articulated_live_pose_loads = false, # articulated spacecraft: per-link loads at the live link poses
    results_layout = :wide               # :long writes one row per spacecraft and saved time (see Simulation Outputs)
)
```

Set `results = false` to run without writing any files, for example when
profiling or when only a returned solution matters. `normalize = true` is
refused because `run_simulation` propagates the state in SI units
(`SPACEAGORA_ALLOW_TYPED_NORMALIZE=1` turns the refusal into a warning, for
legacy checks only). With checkpoints enabled, `checkpoint_interval_s` must be
positive. With `resume_from_checkpoint = true` and no checkpoint on disk, the
run warns and starts from the initial conditions. `articulated_live_pose_loads`
is described on
[Articulated and Compliant Vehicles](articulated_compliant.md).

`run_simulation` never reads `generate_plots`; the flag is only carried along
when settings are copied. `make_example_config` sets it to `false`, which
changes nothing about what the run writes: the CSV, Feather and manifest
files. The verification study switches its own plots with the command-line
option `--plots=`, and the four PNG files the quickstart saves under
`output/plots/` come from that example script's own plotting code, which reads
the saved CSV. For any other run, make plots from the saved files as
[Simulation Outputs](outputs.md) describes.

## FilePaths

`FilePaths` records where external assets are kept, but `run_simulation` does
not read it: the models that need assets take their paths directly (a planet
constructor such as `Earth("", SPICE_PATH)` takes the SPICE directory, and
`GravitationalHarmonicsModel` its coefficient file). The defaults are fine for
every run, and the field can be left out:

```julia
SM.FilePaths(
    results              = "Results",
    GRAM                 = "data/GRAMSuite.jl/GRAM Suite 2.0",
    SPICE                = "data/GRAMSuite.jl/GRAM Suite 2.0/SPICE",
    topography_harmonics = "data/planets/topography",
    gravity_harmonics    = "data/planets/gravity"
)
```

## Using `make_example_config`

For quick studies, `make_example_config` assembles the whole configuration in
one call. It and `make_three_body_spacecraft` are exported by SpaceAGORA, so
`using SpaceAGORA` is the only import:

```julia
using SpaceAGORA

planet = make_no_gram_planet(:earth)

spacecraft = make_three_body_spacecraft(
    bus_dims        = (2.05, 2.05, 2.8),
    panel_dims      = (0.01, 2.85, 1.0),
    bus_mass        = 620.0,
    panel_mass_each = 10.0,
    panel_offset_y  = 2.05/2.0 + 2.85/2.0,
    ic              = InitialCondition(ra=planet.Rp_e + 500e3, rp=planet.Rp_e + 500e3,
                                     i=28.5, ω=10.0, Ω=20.0, ν=0.0),
    prop_mass       = 200.0,
    id              = 1
)

config = make_example_config(
    planet             = planet,
    spacecraft         = spacecraft,
    mission_time       = 3600.0 * 12.0,
    initial_time       = InitialTime(year=2024, month=1, day=1),
    dynamic_effectors  = (InverseSquaredJ2GravityModel(),),
    density_model      = NoAtmosphereModel(),
    ephemerides_model  = SimpleEphemeridesModel(),
    orientation_sim    = false,
    keplerian          = true,
    EI_km              = 120.0,
    verbose            = true
)

run_simulation(config)
```

Left to its defaults, `make_example_config` uses J2 gravity, no atmosphere,
SPICE ephemerides, `MaxwellianHeat` with an accommodation factor of 1,
`EI_km = 300.0`, no winds or topography, a `MissionTime` mission, no GNC
models, and the tolerances listed on
[Solver Configuration](solver_configuration.md#Practical-recipes). It writes
results (`results = true`) to the `output/` directory at the repository root.
For a first script built on it, see [Your First Script](first_script.md).
Build `SimulationConfiguration` directly, as described above, when you need to
choose its components yourself.

## `make_three_body_spacecraft`

Builds a three-body spacecraft: a bus and two identical solar panels, placed
symmetrically on either side of the bus along the body y axis. Most repository
examples use this geometry. Its keywords, with letters standing for your own
values:

```julia
make_three_body_spacecraft(
    bus_dims        = (x, y, z),         # bus bounding box, m
    panel_dims      = (t, span, chord),  # panel thickness, span of one panel (half the array), chord, m
    bus_mass        = 620.0,             # kg
    panel_mass_each = 10.0,              # kg per panel
    panel_offset_y  = offset,            # distance along body y from bus center to each panel center, m
    ic              = SM.InitialCondition(...),
    prop_mass       = 200.0,             # propellant mass, kg (default 0.0)
    id              = 1                  # spacecraft identifier (default 1); CSV prefixes use list position
)
```

It also takes `reflection_coefficient` (default `1.0`), `bus_ram_face`
(`:legacy` or `:frontal`, which bus face sets the reference area) and optional
attitude quaternions for the bus and each panel; the comments in
`src/orchestration/simulation/config/example_configuration.jl` explain them.

Saved columns use the spacecraft's one-based position in the list passed to
`DynamicsModel`, not its `id`. For `DynamicsModel([spacecraft_a, spacecraft_b],
effectors)`, the first spacecraft gets `sc1_*` columns and the second gets
`sc2_*`, even if their IDs are `42` and `7`. A single spacecraft with `id=42`
still gets `sc1_*` columns. See [Output Columns](output_columns.md).

## Panel angles and heating

With a built-in aerodynamic model (`AerodynamicCoefficientfM`,
`AerodynamicCoefficientConstant` or `AerodynamicCoefficientNoBallisticFlight`),
heating uses the current panel geometry. If `orientation_sim = true`, it
follows the propagated spacecraft attitude and each panel's orientation;
otherwise it follows the aerodynamic model's fixed-attitude incidence policy. A
change of panel command reaches the heating without waiting for a force
evaluation. A drag scale such as the verification study's `cd_scale` multiplies
the aerodynamic force only; it changes heating only by changing the
trajectory.

`AerodynamicCoefficientfM`'s `fixed_attitude_incidence` (`:max_drag` by
default, `:attitude` or `:tumbling_average`) sets the incidence for both
aerodynamics and heating when attitude is not propagated.
`AerodynamicCoefficientConstant` and `AerodynamicCoefficientNoBallisticFlight`
always use `:max_drag` in that case. With propagated attitude, incidence
follows the flow relative to the spacecraft; without it, the fixed-attitude
policy supplies the geometric angle.

A configuration without one of these built-in models heats with the stored
`Link.α` angle. A custom effector does not take over the incidence by
declaring `environment_requirements(model).atmosphere = true`; its heating
input stays the angle that the effector or a controller maintains. If a
built-in model is also present, it sets the incidence. When attitude is not
propagated, several built-in models must share one fixed-attitude policy;
conflicting policies raise an error at the first thermal sample, because there
is then no single angle for heating to use.

For a controlled panel, heating follows the orientation the panel actually
has. Its incidence equals the panel's commanded angle for the supported
geometry, a panel offset along the body y axis with a flow-aligned root; other
offsets or root attitudes can give a different incidence without changing the
stored command. In the Odyssey energy-depletion example, the saved
maximum-link heat columns include the uncontrolled bus, while the controller's
panel heat limits apply only to the panels it controls.
