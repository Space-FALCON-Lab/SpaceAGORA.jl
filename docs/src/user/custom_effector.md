# Adding a Force or Torque of Your Own

Use this page when you want the simulation to apply a force or a torque that
SpaceAGORA does not already model: a thruster with your own throttle law, a
drag model from a paper, a tether, a magnetic torquer with your own field
model, or a disturbance you want to study in isolation.

This page is for a student who has run the [Quickstart](quickstart.md) and can
read a short Julia script. You do not need to have read SpaceAGORA's source.

Shortest successful command (the complete example below saved as a file):

```text
julia --project=. my_effector.jl
```

What to read next:

- [The Integrated State](integrated_state.md)
- [Stopping a Simulation on a Condition](stop_conditions.md)
- [Simulation Configuration](simulation_configuration.md)
- [Extensibility](../extensibility.md) (the full interface reference)

## What you will accomplish

You will write a *dynamic effector*: a small Julia type that, at every
evaluation of the equations of motion, is handed the spacecraft's current
state and returns the force and the torque it exerts. You will attach it to a
simulation, run the simulation with and without it, and check that its effect
on the orbit matches a hand calculation. The whole example runs on a fresh
clone with no GRAM or SPICE data.

## The lines that teach the new action

An effector is a type with one method. This one is a thrust that always
points along the velocity, which is the simplest effector whose effect can be
predicted by hand: it does work on the spacecraft, so the orbit's energy rises
at a known rate.

```julia
struct AlongTrackThrust <: SpaceAGORA.AbstractForceTorqueModel
    force_n::Float64          # thrust magnitude, newtons, always along the velocity
end

function SpaceAGORA.wrench(model::AlongTrackThrust, x::SpaceAGORA.StateSample,
                           env::SpaceAGORA.EnvironmentSample, t::Float64)
    direction = x.vel_ii / norm(x.vel_ii)          # unit vector along the inertial velocity
    force_ii = model.force_n * direction            # newtons, inertial frame
    torque_body = SVector{3, Float64}(0.0, 0.0, 0.0)  # newton-metres, bus-fixed frame
    return force_ii, torque_body
end
```

It is attached by listing it with the other effectors of the simulation;
gravity is an effector too:

```julia
SM.DynamicsModel([spacecraft], (SM.InverseSquaredGravityModel(), AlongTrackThrust(0.052)))
```

Everything else in the script below is the standard setup of a simulation.

## The complete script

Save it as `my_effector.jl` in the repository root and run it with
`--project=.`. Two names come from SpaceAGORA rather than from you: `SM` is a
short alias for the module that holds every configuration type, and
`make_three_body_spacecraft` builds the bus-plus-two-panels spacecraft that
all repository examples use (it is imported from the module that owns the
example helpers, as the [Simulation Configuration](simulation_configuration.md)
page describes). Everything between `build_config`'s parentheses is the
standard configuration explained on that page.

```julia
using SpaceAGORA
using StaticArrays
using LinearAlgebra
using CSV, DataFrames
const SM = SpaceAGORA.SimulationModel
import SpaceAGORA.TelemetryVerification: make_three_body_spacecraft

# --- 1. The effector -------------------------------------------------------------
struct AlongTrackThrust <: SpaceAGORA.AbstractForceTorqueModel
    force_n::Float64          # thrust magnitude, newtons, always along the velocity
end

function SpaceAGORA.wrench(model::AlongTrackThrust, x::SpaceAGORA.StateSample,
                           env::SpaceAGORA.EnvironmentSample, t::Float64)
    direction = x.vel_ii / norm(x.vel_ii)          # unit vector along the inertial velocity
    force_ii = model.force_n * direction            # newtons, inertial frame
    torque_body = SVector{3, Float64}(0.0, 0.0, 0.0)  # newton-metres, body frame
    return force_ii, torque_body
end

# --- 2. A simulation that uses it ---------------------------------------------------
planet = SM.make_no_gram_planet(:earth)
altitude_m = 400e3
a0 = planet.Rp_e + altitude_m
period_s = 2π * sqrt(a0^3 / planet.μ)

function build_config(effectors; results_directory::String)
    ic = SM.InitialCondition(a=a0, e=0.001, i=51.6, ω=0.0, Ω=0.0, ν=0.0)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=ic, prop_mass=0.0, id=1)
    return SM.SimulationConfiguration(
        file_paths=SM.FilePaths(),
        simulation_settings=SM.SimulationSettings(
            results=true, verbose=false, results_directory=results_directory,
            generate_plots=false, generate_filenames=false, normalize=false, save_csv=true),
        mission_configuration=SM.MissionConfiguration(
            mission_type=SM.MissionTime, mission_time=period_s, number_of_orbits=1,
            keplerian=true, orientation_sim=false, num_steps_to_save=1000, data_rate=10.0),
        environment_model=SM.EnvironmentModel(
            planet=planet, EI=120.0, density_model=SM.NoAtmosphereModel(),
            ephemerides_model=SM.SimpleEphemeridesModel(),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=SM.DynamicsModel([spacecraft], effectors),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances())
end

thrust = AlongTrackThrust(0.052)     # 0.052 N on 520 kg: 1e-4 m/s^2
config_without = build_config((SM.InverseSquaredGravityModel(),); results_directory="output/effector_walkthrough/without")
config_with    = build_config((SM.InverseSquaredGravityModel(), thrust); results_directory="output/effector_walkthrough/with")

run_simulation(config_without)
run_simulation(config_with)

# --- 3. Verify the effect --------------------------------------------------------
function semi_major_axis(csv_path)
    df = CSV.read(csv_path, DataFrame)
    last = df[end, :]
    r = norm((last.sc1_pos_1, last.sc1_pos_2, last.sc1_pos_3))
    v = norm((last.sc1_vel_1, last.sc1_vel_2, last.sc1_vel_3))
    energy = v^2 / 2 - planet.μ / r          # specific orbital energy, J/kg
    return -planet.μ / (2 * energy), last.time, last.sc1_mass
end

a_without, t_end, mass_kg = semi_major_axis("output/effector_walkthrough/without/simulation_results.csv")
a_with, _, _ = semi_major_axis("output/effector_walkthrough/with/simulation_results.csv")
accel = thrust.force_n / mass_kg
v0 = sqrt(planet.μ / a0)
expected_da = 2 * a0^2 / planet.μ * accel * v0 * t_end
observed_da = a_with - a_without
println("mass used by the engine: $(mass_kg) kg, thrust acceleration: $(accel) m/s^2")
println("semi-major axis after one orbit without thrust: $(round(a_without / 1e3, digits=3)) km")
println("semi-major axis after one orbit with thrust:    $(round(a_with / 1e3, digits=3)) km")
println("expected rise: $(round(expected_da, digits=1)) m, observed rise: $(round(observed_da, digits=1)) m, ratio $(round(observed_da / expected_da, digits=3))")
```

## The important choices, in plain language

**A type and one method.** An effector is a Julia `struct` that is a subtype
of `SpaceAGORA.AbstractForceTorqueModel` and carries its own parameters (here
the thrust magnitude). What makes it an effector is one method of the function
`SpaceAGORA.wrench`. The engine calls that method with four arguments:

- `model`, your struct;
- `x`, a *state sample*: the spacecraft's position `x.pos_ii` and velocity
  `x.vel_ii` (metres and metres per second, in the planet-centred inertial
  frame), its current mass `x.mass_kg`, and, when attitude is propagated, its
  attitude quaternion `x.q_ib` and body angular velocity `x.ω_body` (both
  `nothing` otherwise);
- `env`, an *environment sample* with the planet model and, on request, the
  atmosphere, the planet-fixed frame, the Sun and third bodies (see the
  reference below);
- `t`, the elapsed simulation time in seconds.

It returns a pair: the **force in newtons in the inertial frame**, and the
**torque in newton-metres in the bus-fixed frame**. Return zero vectors for
the part you do not model.

**Force, not acceleration.** The engine sums the forces of all effectors and
divides by the spacecraft's current mass itself. If your physics is naturally
an acceleration (a gravity model, a perturbation from a paper), multiply by
`x.mass_kg` before returning, which is what the built-in solar-radiation model
does. Returning an acceleration where a force is expected under-drives the
spacecraft by a factor of its mass.

**Frames.** The subscript `ii` on `pos_ii`, `vel_ii` and `force_ii` means
"inertial frame, inertial components": the non-rotating frame centred on the
planet in which the orbit is integrated. Torque is expressed in the frame
fixed to the spacecraft bus, the same frame in which the spacecraft's inertia
tensor is given, so the two enter the attitude equations together. If your
force is known in the bus frame, rotate it into the inertial frame with the
attitude quaternion first. With `orientation_sim=false` no attitude is
propagated and `x.q_ib` is `nothing`; the spacecraft's configured geometry and
the fixed orientations of its links still exist and are what the built-in
aerodynamic models use in that mode.

**Which method the engine uses.** The engine asks whether a method of
`wrench` accepts the argument types `(YourType, StateSample,
EnvironmentSample, Float64)`. Annotating the model argument with your type is
what makes the method yours. The other three arguments may be left untyped, or
annotated with the types above or with types that contain them (`t::Real`
works); all three forms were checked. A method annotated with a type that does
not accept the argument (say `x::AbstractVector`) is not a match, so for an
ordinary custom effector the engine falls back to the older `calcForceTorque`
hook and the run stops with a `MethodError` that names `calcForceTorque`.

**Attaching it.** Effectors are passed as a tuple to `DynamicsModel` (or to
`make_example_config`'s `dynamic_effectors`). Order does not matter for the
physics: the forces are summed. Gravity is an effector like any other, so
`InverseSquaredGravityModel()` is listed explicitly.

**Two runs, not one.** The cleanest way to see what an effector does is to run
the same configuration with and without it and compare the results files. The
example uses two-body gravity, so without thrust the orbit's semi-major axis
stays constant to the metre over one orbit and every change with thrust is the
effector's.

## Expected result

Tested on `main` at commit `80240c2b` (September 2026), on a fresh clone with
no GRAM or SPICE data. The script runs for about a minute and printed:

```text
mass used by the engine: 520.0 kg, thrust acceleration: 9.999999999999999e-5 m/s^2
semi-major axis after one orbit without thrust: 6778.137 km
semi-major axis after one orbit with thrust:    6779.118 km
expected rise: 981.8 m, observed rise: 981.9 m, ratio 1.0
```

These are the recorded numbers for this configuration, not universal values.
The mass is the bus plus two panels. The expected rise comes from the rate at
which a tangential acceleration raises a near-circular orbit,
`da/dt = 2 a² a_t v / μ`, integrated over one period; the observed rise is
read from the final rows of the two results files, with the semi-major axis
computed from position and velocity through the orbital energy. Success is a
ratio close to one. A ratio near `1/520` means you returned an acceleration
instead of a force; a negative rise means the force points against the
velocity; a `MethodError` naming `calcForceTorque` means the `wrench` method
was not matched (check its annotations).

Both runs write the standard three files under their own directory, see
[Simulation Outputs](outputs.md). The thrust does not appear as a column of
its own: only the built-in aerodynamic model writes per-force diagnostics.

## What you can change next

- **Point the thrust somewhere else.** Cross-track (`cross(pos, vel)`
  normalised) changes the inclination; radial (`x.pos_ii / norm(x.pos_ii)`)
  changes the eccentricity. The two-run comparison still works.
- **Make it depend on time or position.** `t` is elapsed seconds; a burn
  window is one `if`. `norm(x.pos_ii) - env.planet.Rp_e` is the height above
  a sphere of the planet's equatorial radius, which is not the same as the
  ellipsoidal altitude saved in the results column `sc1_altitude` (see
  [The Integrated State](integrated_state.md)); choose the surface your
  condition is meant to refer to.
- **Use the atmosphere.** Declare what your effector needs and the engine
  samples it for you at every evaluation:

  ```julia
  SpaceAGORA.environment_requirements(::MyDrag) =
      SpaceAGORA.EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
  ```

  Then `env.atmosphere.rho_kg_m3` is the density at the current position and
  `env.planet_frame.vel_pp` is the velocity expressed in the planet-fixed
  (rotating) frame. That is the air-relative velocity only when the wind is
  zero; a drag model must subtract the wind, `env.atmosphere.wind_pp`, after
  bringing it into the same frame (the engine's thermal callback shows the
  conversion from the wind's east, north and up components, in that order). With
  `NoAtmosphereModel()` the density is zero; `ExponentialAtmosphereModel(planet)`
  gives a simple open-data profile with zero wind.
- **Add a torque.** Return a non-zero second element and run with
  `orientation_sim=true`; the attitude columns `sc1_q_*` then appear in the
  results, and `x.q_ib` and `x.ω_body` are available inside `wrench`.
- **Give it state of its own.** Keep parameters in the struct fields. If the
  effector must remember something between evaluations (an integrator of its
  own, a mode), prefer a control effector (see [Extensibility](../extensibility.md)):
  `wrench` is meant to be a pure function of its four arguments, because the
  solver evaluates it at trial points it later rejects.

## Reference: the interface as it exists today

This section names the types and functions; the walkthrough above did not
need them by name. Tested on the walkthrough commit: the `wrench` method
matching (untyped arguments, the exact annotations, `t::Real`, and the
fallback on a non-matching annotation) and the force-to-acceleration
behaviour (through the expected rise). Inspected in the source, not run: the
field lists of the sample types, the harmonics exception, `solver_partition`,
and the gravity-backbone hooks.

| Symbol | Role |
|---|---|
| `SpaceAGORA.AbstractForceTorqueModel` | Supertype of every dynamic effector. |
| `SpaceAGORA.wrench(model, x::StateSample, env::EnvironmentSample, t::Float64)` | The hook. Returns `(force_ii::SVector{3}, torque_body::SVector{3})`, SI units, force inertial, torque bus-fixed. Must be pure. |
| `SpaceAGORA.StateSample` | Fields `pos_ii`, `vel_ii` (m, m/s, inertial), `mass_kg`, `q_ib` (inertial-to-body quaternion, scalar-last `[x, y, z, w]`, or `nothing`), `ω_body` (rad/s, or `nothing`), `spacecraft` (the typed spacecraft model, for geometry and inertia). |
| `SpaceAGORA.EnvironmentSample` | Fields `planet` (always), `planet_frame` (`alt_m`, `lat_rad`, `lon_rad`, `pos_pp`, `vel_pp`, `l_pi`: position and velocity in the planet-fixed frame and the inertial-to-planet-fixed rotation matrix), `atmosphere` (`rho_kg_m3`, `temperature_k`, `wind_pp` as east, north, up components, in that order), `solar` (`sun_pos_ii`), `third_bodies` (`names`, `positions_ii`). Each optional field is `nothing` unless requested. |
| `SpaceAGORA.environment_requirements(model)` | Returns `EffectorEnvironmentRequirements(planet_frame=, atmosphere=, solar=, third_body_names=)`; the default requests nothing. |
| `SpaceAGORA.calcForceTorque(model, x, p, i)` | The older hook, kept for existing models; `x` is the spacecraft's raw state block, `p` the engine parameters, `i` the spacecraft index. Same return convention. Used when no `wrench` method matches; the built-in gravitational-harmonics model is routed to it deliberately for its buffer reuse. |
| `SpaceAGORA.solver_partition(model)` | `:explicit` (default) or `:implicit`; only matters under the `split_imex` solver mode. |
| `SM.DynamicsModel([spacecraft], effectors)` | Where the effector tuple is attached. |

How the engine uses the result: the forces of all effectors are summed per
spacecraft and divided by the current mass to give the translational
acceleration; the torques are summed in the bus-fixed frame and drive the
attitude equations when `orientation_sim=true`, and are ignored otherwise.
Solver modes other than the default are described on
[Solver Configuration](solver_configuration.md); the gravity-backbone hooks
listed on the Extensibility page are for effectors that take part in that
special mode and are not needed for an ordinary force.
