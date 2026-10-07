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

## Joint types and articulated spacecraft

A spacecraft is **articulated** if and only if at least one of its `Joint`s is not
`:fixed`. Articulated spacecraft integrate their joint coordinates beside the bus state;
every other spacecraft, including one whose joints are all `:fixed`, runs exactly as before
(no tree is built and nothing is validated). Run `examples/Articulated_Panels_Demo.jl` for a
working case.

### Declaring joints

```julia
bus = Link(root=true, m=620.0, dims=MVector(2.05, 2.05, 2.8))
panel = Link(m=10.0, dims=MVector(0.01, 2.85, 1.0), r=MVector(0.0, 2.45, 0.0))   # COM in the bus frame
add_joint!(sc, Joint(bus, SVector(0.0, 1.025, 0.0), panel, SVector(0.0, -1.425, 0.0);
    joint_type=:hinge, axis=[1.0, 0.0, 0.0], stiffness=60.0, damping=4.0, initial_q=0.2))
```

`Joint(link1, p1ᵇ, link2, p2ᵇ; joint_type=...)` connects the parent `link1` to the child `link2`
at a joint point given in the parent frame (`p1ᵇ`) and in the child frame (`p2ᵇ`). Joint
coordinate 0 is the configured geometry: the links' `r` (COM, bus frame) and `q` (attitude
relative to the bus, scalar-last). For non-fixed joints `p1ᵇ` and `p2ᵇ` must map to the same
point in that geometry (1e-9 m); fixed joints merge at the configured geometry as given.

| `joint_type` | Coordinates | `axis` | `stiffness`, `damping` | `rest` | `initial_q`, `initial_qd` |
|---|---|---|---|---|---|
| `:fixed` (default) | none; the child is merged into the parent body | unused | unused | unused | unused |
| `:hinge` | 1 angle (rad) about `axis` | required, parent frame, normalized | scalars (N m/rad, N m s/rad) | angle (default 0) | scalars (default 0) |
| `:slide` | 1 displacement (m) along `axis` | required, parent frame, normalized | scalars (N/m, N s/m) | displacement (default 0) | scalars (default 0) |
| `:ball` | scalar-last quaternion (3 rate components) | unused | scalar or 3x3 PSD matrix | quaternion (default identity) | quaternion (default identity), parent-frame angular velocity 3-vector (default 0) |

Gains must be `>= 0`. Springs and dampers act in joint space: `τ = -k (q - rest) - c q̇`; for a
ball joint the spring uses the axis-angle of `rest⁻¹ ⊗ q` and the damper the relative angular
velocity (the 3x3 stiffness is exact for an isotropic `k I`, first-order for an anisotropic one).
Existing joints keep working and are `:fixed`; the legacy `Kx/Kt/Cx/Ct` fields are unused by this
feature. Links whose joint is `:fixed` move rigidly with their parent body.

### Conventions and state

- The root body is the root link plus every link merged into it through `:fixed` joints plus
  the propellant. The engine's `pos` and `vel` are the **root composite center of mass**; the
  system center of mass is saved as `sc{i}_system_com_*`.
- Quaternions are scalar-last and are the body-to-inertial rotation, like the bus `q`.
- `joint_q` and `joint_qd` are appended to the spacecraft's state, listing the non-fixed
  joints in `spacecraft.joints` order. A ball joint has four coordinates and three rates, and
  its rate is the angular velocity of the child relative to its parent in the parent frame.
- Saved outputs: `joint_q`, `joint_qd`, `articulated_link_pose` (every link's inertial COM and
  attitude from the joint kinematics) and `system_com`; see [Simulation Outputs](outputs.md).
- Checkpoint and resume work: the joint state is part of the checkpointed state.

### Mass and inertia

- `mass` in the state keeps its meaning: total spacecraft mass, dry plus propellant, so mass flow
  and `sc1_mass` are unchanged. The root body mass at run time is the state mass minus the
  moving bodies' masses (which never change).
- The link masses and inertias are used. The spacecraft's `inertia_tensor` is **ignored** in
  articulated mode, and `dry_mass` must equal the sum of the link masses (checked at setup).
- Propellant is a point mass at the root composite COM: it adds mass and no inertia and never
  moves the COM, so the root inertia stays the configured composite.

### Loads in v1

Every non-gravity load (aerodynamics, SRP, thrusters, control torques, third-body gravity) is
evaluated by the existing effectors from the rigid configured geometry and applied to the root
body. Gravity from the position-only gravity models (point mass, J2, harmonics) is evaluated per
body at its own COM, so the gravity-gradient effect of the layout, including joint motion, emerges
from the dynamics; a gravity effector's own `gravity_gradient` flag is superseded by this. A
separate `GravityGradientTorqueModel` still uses the spacecraft `inertia_tensor` on the root: a
known v1 limitation. Per-link aerodynamics and SRP on moving links, joint motors and joint limits
are not implemented.

### Solver advice and supported routes

| Route | Support |
|---|---|
| Single-RHS first-order solver modes `:tsit5`, `:auto_stiff`, `:rodas5p`, `:dp8` | supported |
| `:split_imex`, `:multirate`, `:symplectic`, `:gravity_backbone_split` | refused with an `ArgumentError` |
| Serial, `satellite_batch` and per-satellite RHS routes | supported (one workspace per spacecraft) |
| Flat constellation queue | automatically rerouted to the per-satellite route; forcing `SPACEAGORA_RHS_EXECUTION_MODE=flat` is refused |
| Process routes and `isolate_state` copies | supported: runtime data is rebuilt from the configuration inside each run |
| `orientation_sim=false`, reaction wheels, a robot-arm effector on the same spacecraft | refused |
| Kinematic panel-angle control (`SolarPanelAngleOfAttackControlModel`) | works on links of the root body exactly as before; an error on links of a moving body until joint motors exist |
| Mixing articulated and rigid spacecraft in one run | refused (the state holds equal-sized blocks); articulated spacecraft must share link and joint counts |

Joint dynamics can be stiff: for a stiff hinge keep `dt_max_orbit` below a fraction of the hinge
period, and use tight tolerances; `:dp8` was used for the conservation checks. The joint
coordinates use the quaternion and angular-rate tolerances.

!!! note
    Per-link loads, motors and limits are not part of this release; see "Loads in v1" above.

## Compliant attachments

A `CompliantAttachment` mounts a `CompliantMultibodyModel` (a cloth panel mesh, a
flexible appendage) on a spacecraft link. The attachment bodies move under their own
compliant joint springs, dampers and actuators and exchange forces and torques with the
link in both directions. `examples/Cloth_Panel_Attachment_Demo.jl` runs the four-panel cloth
deployment of `Solar_Panel_Cloth_Deployment_Demo.jl` inside `run_simulation` this way.

```julia
build = build_rectangular_compliant_grid(3, 4; anchor_index=1)       # or any CompliantTopologyBuild / model
att = CompliantAttachment(;
    model=build,                       # a CompliantMultibodyModel, or a build (its state is the initial state)
    link=bus,                          # root, a fixed-merged link, or a link of a moving body
    mount_point=(0.5, 0.0, 0.0),       # mount frame origin in the link frame (m)
    mount_quaternion=(0.0, 0.0, 0.0, 1.0),   # mount frame orientation relative to the link frame
    joint_actuators=actuators,         # CompliantJointActuator list
    rest_schedule=(out, t) -> (out[1] = rest_quaternion(t); nothing),   # optional
)
sc = SpacecraftModel(; links=[bus], root=bus, inertia_tensor=bus.inertia, initial_condition=ic, attachments=[att])
```

### Frames and conventions

- The **mount frame** sits at `mount_point` in the link frame with orientation
  `mount_quaternion` relative to it. The link frame origin is the link center of mass (the
  spacecraft position for the root link; the link's `r` and `q` relative to the bus otherwise)
  and its axes are the link body axes. In an articulated spacecraft the link frame follows the
  link's dynamic body.
- The model's joints with `parent == 0` attach to the mount frame instead of a fixed base: the
  model's `base_position` and `base_quaternion` are **ignored**. At least one such joint is
  required. A topology build's own state, given for its base pose, is converted to mount-frame
  coordinates; an explicit `initial_state` is already in the mount frame (13 numbers per body:
  position, scalar-last quaternion, velocity, body-frame angular velocity, relative to a mount at
  rest). With neither, the model starts at its rest geometry (a spanning tree from the mount joints;
  pass a build's state for meshes with closed loops).
- The run state (`att_r`, `att_q`, `att_v`, `att_ω`) holds each body's position and velocity RELATIVE to the
  spacecraft's `pos` and `vel` (inertial axes) and its absolute attitude and body rate; the initial state is the
  mount-frame state moved by the mount's pose and rigid-body velocity about the bus. The saved `attachment_pose`
  and `system_com` are inertial (`pos + att_r`). (The robot-arm state `arm_*` still holds absolute positions: a known follow-up.)
- `rest_schedule(out::Vector{SVector{4,Float64}}, t)` writes every joint's rest quaternion for time `t`
  (seconds since the start) into `out`; `nothing` keeps each joint's own rest quaternion. It is called
  on every RHS evaluation with the stage time, so a smoothstep deployment is followed exactly (a
  standalone stepper that freezes the rest per step lags by half a step), and the run is
  allocation-free when the schedule is.

### How the loads reach the spacecraft

The joint loads and the gravity difference are computed from relative quantities. The bus (or articulated base)
acceleration, which includes the reactions, is then known, and each body's relative acceleration is
`F_i/m_i + g(r_base + r_rel) - g(r_base) + (g(r_base) - a_base)`: gravity is evaluated at each body's absolute position
but only its difference from the base gravity enters, so the large common term never reaches the relative acceleration.

- The compliant joint springs, dampers and actuators use the same math as `compliant_joint_loads`
  (`compliant_joint_loads_in_place!` is the allocation-free variant the engine calls).
- The reaction on the mount, a force at the mount point plus a torque, goes to the link it is mounted
  on. **Rigid spacecraft:** the force is added to the bus force (inertial) and the torque, including
  the `r_mount x F` lever, to the bus torque (body frame) before the bus equations. **Articulated
  spacecraft:** the body kinematics are computed first, then the attachment loads, and the reaction enters
  `articulated_dynamics!` as a per-body wrench (`body_force_world`, `body_torque_world`, about the body
  COM) that is mapped to the generalized forces through each body's Jacobian columns, so the hinge,
  slide and ball coordinates and the root respond to the attachment.
- **Loads in v1:** gravity acts on every attachment body at its own position (the position-only
  gravity effectors). No aerodynamics, SRP, thermal or other effector acts on attachment bodies.

### Mass bookkeeping

Attachment bodies are not links: they are not in `links`, `dry_mass` or the state `mass`, and carry
their own masses (`attachment_total_mass(sc)`). The system total is the spacecraft plus the attachments;
the saved `system_com` includes them. A body modeled as a `:fixed` link instead adds to `dry_mass`.
Thrust, mass flow and every other effector see the spacecraft mass only.

### Solver advice and supported routes

| Route | Support |
|---|---|
| `:tsit5`, `:auto_stiff`, `:rodas5p`, `:dp8`, rigid or articulated spacecraft | supported |
| `:split_imex`, `:multirate`, `:symplectic`, `:gravity_backbone_split` | refused with an `ArgumentError` |
| Serial and per-satellite RHS routes, `isolate_state` copies, checkpoint and resume | supported |
| Flat constellation queue | automatically rerouted to the per-satellite route; forcing `SPACEAGORA_RHS_EXECUTION_MODE=flat` is refused |
| `orientation_sim=false`, a robot-arm effector on the same spacecraft | refused |
| Attachment link not in the spacecraft, model without a joint to the mount | refused |
| Mixing spacecraft with different attachment body counts (or with and without attachments) in one run | refused (equal-sized state blocks) |

Cloth meshes are stiff in general. Use `:auto_stiff` or `:rodas5p` when the fastest attachment
mode (about `sqrt(k/m)` for the translational springs and `sqrt(k_rot/I)` for the rotational ones,
plus the damping rate `c/m`) is much faster than the orbit and attitude dynamics of interest; the
explicit modes are fine for mild meshes such as the demo (about 30 rad/s). The attitude and
angular-rate tolerances apply to `att_q` and `att_ω`. Because the attachment state is relative to the bus, spring forces carry no roundoff from orbital-radius
positions; a 2e5 N/m attachment integrates with ordinary tolerances.

### Limits

No per-body aerodynamics, SRP or thermal loads, no contact, no attachment-to-attachment joints,
and no attachment on a spacecraft with the robot-arm effector. The `arm_*` robot-arm path is separate
and unchanged.

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
