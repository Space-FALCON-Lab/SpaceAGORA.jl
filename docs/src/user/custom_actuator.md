# Adding an Actuator of Your Own

Use this page when your spacecraft carries hardware SpaceAGORA does not model
yet (a cold-gas nozzle, an electric thruster, a different kind of wheel) and
you want a controller to command it while the engine works out what it does to
the vehicle.

You describe a piece of hardware once, as an *actuator kind*: the force and
torque it puts on the vehicle at a commanded level, the bounds of that level,
and what it consumes. You mount it on a spacecraft link, write a controller
that only decides the level, and run a simulation in which the engine applies
the force and torque, clips the level to the hardware's limits and takes the
propellant out of the vehicle's mass. Then you check the result against the
rocket equation. The page assumes you have run the [Quickstart](quickstart.md)
and can read a short Julia script, and it follows the same pattern as
[Adding a Force or Torque of Your Own](custom_effector.md).

Saved as `my_actuator.jl`, the complete example below runs with:

```text
julia --project=. my_actuator.jl
```

Related pages:

- [Extensibility](../extensibility.md) (the interface, and the built-in kinds)
- [Adding a Force or Torque of Your Own](custom_effector.md)
- [Configuring a Simulation](simulation_configuration.md)

## The actuator in brief

The kind is a type and two methods, plus one optional method. In this outline
the method bodies are elided (`...`); the complete script below has them:

```julia
struct ColdGasNozzle <: SpaceAGORA.AbstractActuator
    axis::SVector{3, Float64}        # thrust direction in the link frame (unit vector)
    location::SVector{3, Float64}    # nozzle position in the link frame, m
    max_thrust_n::Float64
    isp_s::Float64
end

SpaceAGORA.actuator_limits(n::ColdGasNozzle) = (0.0, n.max_thrust_n)
SpaceAGORA.actuator_wrench(n::ColdGasNozzle, level, ctx::SpaceAGORA.ActuatorContext) = ...
SpaceAGORA.actuator_consumption(n::ColdGasNozzle, level, ctx) = ...
```

The controller declares that it commands mounted hardware and holds one level
per mounted actuator:

```julia
SpaceAGORA.actuation_mode(::BurnFor) = :mounted
SpaceAGORA.actuator_levels(c::BurnFor, spacecraft_index::Int) = c.levels
```

The nozzle is mounted with `add_actuator!(spacecraft, root, nozzle)`.

## The complete script

Save it as `my_actuator.jl` in the repository root or in a
[separate user project](installation_environment.md#Your-own-project), and run
it with that project's environment. It uses only names exported by
`SpaceAGORA`. The same script is `examples/walkthrough_env/my_actuator.jl`. CI
runs it from that separate project on every pull request, and a docs gate
checks that the code on this page matches it character for character.

```julia
using SpaceAGORA
using StaticArrays
using LinearAlgebra

# --- 1. The actuator kind ---------------------------------------------------------
struct ColdGasNozzle <: SpaceAGORA.AbstractActuator
    axis::SVector{3, Float64}        # thrust direction in the link frame (unit vector)
    location::SVector{3, Float64}    # nozzle position in the link frame, m
    max_thrust_n::Float64
    isp_s::Float64
end

SpaceAGORA.actuator_limits(n::ColdGasNozzle) = (0.0, n.max_thrust_n)

function SpaceAGORA.actuator_wrench(n::ColdGasNozzle, level, ctx::SpaceAGORA.ActuatorContext)
    force_body = ctx.link_to_body * (level * n.axis)            # newtons, body frame
    arm_body = ctx.link_to_body * n.location + ctx.link_r_body  # metres from the body reference point
    return force_body, cross(arm_body, force_body)
end

SpaceAGORA.actuator_consumption(n::ColdGasNozzle, level, ctx) =
    (mass_rate_kg_s = propellant_mass_rate(level, n.isp_s), power_w = 0.0)

# --- 2. A controller that commands it ----------------------------------------------
mutable struct BurnFor <: SpaceAGORA.AbstractControlEffectorModel
    thrust_n::Float64           # commanded thrust while burning, N
    burn_s::Float64             # burn duration, s
    levels::Vector{Float64}     # one level per mounted actuator, held between ticks
    initial_levels::Vector{Float64}  # the levels each run starts from
end

# Built with one vector of levels: each run starts from a copy of it.
BurnFor(thrust_n, burn_s, levels) = BurnFor(thrust_n, burn_s, collect(Float64, levels), collect(Float64, levels))

SpaceAGORA.actuation_mode(::BurnFor) = :mounted
SpaceAGORA.actuator_levels(c::BurnFor, spacecraft_index::Int) = c.levels

function SpaceAGORA.update_control!(c::BurnFor, u, p, t::Float64, spacecraft_index::Int)
    c.levels[1] = t < c.burn_s ? c.thrust_n : 0.0
    return nothing
end

# At each run's start, the levels the controller was built with.
function SpaceAGORA.initialize_control!(c::BurnFor, u, p, t)
    copy!(c.levels, c.initial_levels)
    return nothing
end

# The held levels through a checkpoint: saved with it, put back on a resumed run.
SpaceAGORA.checkpoint_state(c::BurnFor) = (levels=copy(c.levels),)
function SpaceAGORA.restore_checkpoint_state!(c::BurnFor, data)
    copy!(c.levels, data.levels)
    return nothing
end

# --- 3. A spacecraft carrying the nozzle, far from any planet ---------------------
root = Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
ic = CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
spacecraft = SpacecraftModel(links=[root], root=root, prop_mass=20.0, inertia_tensor=root.inertia,
                             initial_condition=ic, id=1)
nozzle = ColdGasNozzle(SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 0.5, 65.0)
add_actuator!(spacecraft, root, nozzle)

thrust_n, burn_s = 0.4, 300.0
planet = make_no_gram_planet(:earth)
config = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=true, verbose=false, results_directory=joinpath(@__DIR__, "output", "actuator_walkthrough"),
        generate_plots=false, generate_filenames=false, normalize=false, save_csv=true),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=600.0, number_of_orbits=1,
        keplerian=false, orientation_sim=true, num_steps_to_save=600, data_rate=1.0),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], ()),              # no gravity: only the nozzle acts
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(BurnFor(thrust_n, burn_s, [thrust_n]),), control_rates=[1.0]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

results = run_simulation(config; return_results=true)

# --- 4. Check it by hand -----------------------------------------------------------------
last = results.table[end, :]
m0 = spacecraft.dry_mass + spacecraft.prop_mass
used_kg = propellant_mass_rate(thrust_n, nozzle.isp_s) * burn_s
dv_rocket = nozzle.isp_s * 9.80665 * log(m0 / (m0 - used_kg))
println("propellant used: $(round(m0 - last.sc1_mass, digits=6)) kg, expected $(round(used_kg, digits=6)) kg")
println("speed gained:    $(round(last.sc1_vel_1, digits=6)) m/s, rocket equation $(round(dv_rocket, digits=6)) m/s")
```

## Choices in the script

**The kind says what the hardware does, not when.** `actuator_wrench` returns
the force (newtons) and the torque (newton-meters) the nozzle puts on the
vehicle at a given level, both in the spacecraft body frame, with the torque
taken about the body reference point. `ctx` carries where the mounting link
sits in that frame (`ctx.link_to_body`, `ctx.link_r_body`) and the vehicle's
mass, attitude and rate, so the same nozzle works on a panel as well as on the
bus.

**Limits and consumption belong to the hardware.** `actuator_limits` bounds
the level: the engine clips every commanded level to it before it evaluates
the wrench, so asking this nozzle for 9 N gives 0.5 N. `actuator_consumption`
returns the propellant mass rate (kg/s) and the electrical power (W).
`propellant_mass_rate(thrust, isp)` is the rocket formula every thrust model
of the package uses, `thrust / (Isp g0)`.

**The controller only decides levels.** `BurnFor` implements no force, torque
or mass-flow hook. Declaring `actuation_mode(c) = :mounted` tells the engine
to derive those from the actuators mounted on the spacecraft, at the levels
`actuator_levels` returns: one level per mounted actuator, in mounting order.
The built-in kinds come first in that order (thrusters, reaction wheels,
torque rods, laser terminals, hinges, thruster gimbals, then control moment
gyros and their inner gimbals), followed by kinds of your own. The engine
calls `actuator_levels` at every evaluation of the equations of motion, so
`update_control!` computes the levels at each control tick and the method only
returns the stored values.

Each run starts from `initial_levels`, the levels the controller was built
with: `initialize_control!` copies them into `levels`, so a second run on the
same controller starts as the first did. They are held until the first tick,
one control period after the start. `checkpoint_state` and
`restore_checkpoint_state!` carry the held levels through a checkpoint, so a
run resumed from a checkpoint file holds the same levels as the uninterrupted
run until its next tick (see
[Run state in checkpoint files](../extensibility.md#Run-state-in-checkpoint-files)).

**Force in the body frame, applied in inertial space.** The nozzle pushes along
body +x. The engine rotates the force with the current attitude at every
evaluation, so a thrust held between ticks turns with the vehicle, as a
body-fixed nozzle does. Mounted actuation needs that attitude state: with
`orientation_sim=false`, the run stops with an `ArgumentError` the first time
the engine evaluates the mounted actuators. Here the attitude stays at the
identity, so body +x is inertial +x.

**No gravity, on purpose.** The dynamics tuple is empty and the spacecraft is
far from Earth, so the nozzle is the only thing acting and the result has a
closed form.

## Expected result

On a fresh clone with no GRAM or SPICE data, the script prints:

```text
propellant used: 0.188255 kg, expected 0.188255 kg
speed gained:    0.230811 m/s, rocket equation 0.230811 m/s
```

The expected propellant is the mass rate at 0.4 N and a specific impulse of
65 s, over the 300 s burn. The speed is compared with the rocket equation
`Δv = Isp g0 ln(m0 / (m0 - Δm))` for the 520 kg vehicle (the bus and its
propellant). Both values are read from the last row of the results table. CI
checks both agreements to the six decimal places printed.

## Variations

- **Turn the vehicle.** Give the initial condition a different attitude, or
  mount the nozzle on a rotated link: the thrust follows the body.
- **Offset the nozzle.** A nonzero `location` gives a torque. Run with a
  torque-free placement first and compare the attitude columns.
- **Mount several.** Two nozzles mean two levels. `allocation_matrix(spacecraft)`
  and `allocate_actuators` turn a desired force and torque into levels.
- **Use a built-in kind.** `Thruster`, `ReactionWheel`, `Magnetorquer`,
  `LaserTerminal` and `PanelHinge`, among others, follow the same contract; see
  [Extensibility](../extensibility.md#Add-an-actuator).
  `examples/Earth_CubeSat_Mounted_Actuators.jl` flies reaction wheels, torque
  rods and a laser terminal together.
- **Let the nozzle feel its surroundings.** Define
  `SpaceAGORA.actuator_requirements(::ColdGasNozzle) = ActuatorRequirements(atmosphere=true)`
  and the context carries the sampled atmosphere in `ctx.environment`.
- **Make it fail.** `add_fault!(spacecraft, ActuatorFault(nozzle; kind=:failed,
  start_s=100.0))` sets the nozzle's level to zero from 100 s into the burn,
  and the speed gained drops accordingly.
- **Read what it did.** The results table has `sc1_actuator_level_1`, the
  level the nozzle applied at each saved step, and `sc1_actuator_power`, the
  power it reported through `actuator_consumption` (zero here).

## Reference: the interface

| Symbol | Role |
|---|---|
| `SpaceAGORA.AbstractActuator` | Supertype of every actuator kind. |
| `SpaceAGORA.actuator_wrench(a, level, ctx)` | Required: `(force_body, torque_body)` at a saturated level. |
| `SpaceAGORA.actuator_limits(a)` | Required: `(lower, upper)` bounds of the level. |
| `SpaceAGORA.actuator_consumption(a, level, ctx)` | Optional: `(mass_rate_kg_s, power_w)`; the default consumes nothing. |
| `SpaceAGORA.actuator_saturate`, `actuator_fraction`, `actuator_unit_wrench` | Optional, with defaults. |
| `SpaceAGORA.actuator_requirements(a)` | Optional: the environment the kind needs sampled, an `ActuatorRequirements`; the default requests nothing. |
| `SpaceAGORA.ActuatorContext` | The mounting link's pose, the vehicle's mass, attitude, rate and magnetic field, and the sampled environment the kind asked for. |
| `SpaceAGORA.add_actuator!(spacecraft, link, actuator)` | Mounting. |
| `SpaceAGORA.actuation_mode(c)`, `SpaceAGORA.actuator_levels(c, spacecraft_index)` | The controller side. |
| `SpaceAGORA.mounted_actuators`, `allocation_matrix`, `allocate_actuators`, `mounted_wrench` | The mounted set and allocation. |

The controller's own hooks, such as `update_control!(c, u, p, t, spacecraft_index)`,
follow the rules on [Extension Conventions](../reference/extension_conventions.md):
the time comes before the spacecraft index, and the index is named
`spacecraft_index`. The
[Extension Hook Reference](../generated/extension_hook_reference.md) lists
when the engine calls each actuator and control hook.
