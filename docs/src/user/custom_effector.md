# Adding a Force or Torque of Your Own

Use this page when you want the simulation to apply a force or a torque that
SpaceAGORA does not already model: a thruster with your own throttle law, a
drag model from a paper, a tether, a magnetic torquer with your own field
model, or a disturbance you want to study in isolation.

You will write a *dynamic effector*: a small Julia type that the engine hands
the spacecraft's current state at every evaluation of the equations of
motion, and that returns the force and the torque it exerts. You attach it to
a simulation, run the simulation with and without it, and compare its effect
on the orbit with a hand calculation. The page assumes you have run the
[Quickstart](quickstart.md) and can read a short Julia script; you do not need
to have read SpaceAGORA's source. The example needs no GRAM or SPICE data, so
it runs on a fresh clone.

Saved as `my_effector.jl`, the complete example below runs with:

```text
julia --project=. my_effector.jl
```

Related pages:

- [The Integrated State](integrated_state.md)
- [Stopping a Simulation on a Condition](stop_conditions.md)
- [Configuring a Simulation](simulation_configuration.md)
- [Extensibility](../extensibility.md) (the full interface reference)

## The effector in brief

An effector is a type with one method. This one is a thrust that always
points along the velocity. It is the simplest effector whose effect can be
predicted by hand: it does work on the spacecraft, so the orbital energy rises
at a known rate.

```julia
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
```

To attach it, list it with the simulation's other effectors in the
`DynamicsModel`; gravity is an effector too. This line is a fragment of the
script below, where `spacecraft` is defined:

```julia
DynamicsModel([spacecraft], (InverseSquaredGravityModel(), AlongTrackThrust(0.052)))
```

Everything else in the script is the standard setup of a simulation.

## The complete script

Save it as `my_effector.jl` in the repository root or in a
[separate user project](installation_environment.md#Your-own-project), and run
it with that project's environment. It uses only names that `using SpaceAGORA`
exports: `make_three_body_spacecraft` builds the bus with its two panels, and
the configuration types come from the same import, so no internal module or
package source file is loaded. Everything between `build_config`'s parentheses
is the standard configuration described on
[Configuring a Simulation](simulation_configuration.md). The same script is
`examples/walkthrough_env/my_effector.jl`. CI runs it from that separate
project on every pull request, and a docs gate checks that the code on this
page matches it character for character.

```julia
using SpaceAGORA
using StaticArrays
using LinearAlgebra
using CSV, DataFrames

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
planet = make_no_gram_planet(:earth)
altitude_m = 400e3
a0 = planet.Rp_e + altitude_m
period_s = 2π * sqrt(a0^3 / planet.μ)

function build_config(effectors; results_directory::String)
    ic = InitialCondition(a=a0, e=0.001, i=51.6, ω=0.0, Ω=0.0, ν=0.0)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=ic, prop_mass=0.0, id=1)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=true, verbose=false, results_directory=results_directory,
            generate_plots=false, generate_filenames=false, normalize=false, save_csv=true),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, mission_time=period_s, number_of_orbits=1,
            keplerian=true, orientation_sim=false, num_steps_to_save=1000, data_rate=10.0),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false),
        dynamics_model=DynamicsModel([spacecraft], effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances())
end

thrust = AlongTrackThrust(0.052)     # 0.052 N on 520 kg: 1e-4 m/s^2
output_root = joinpath(@__DIR__, "output", "effector_walkthrough")
config_without = build_config((InverseSquaredGravityModel(),); results_directory=joinpath(output_root, "without"))
config_with    = build_config((InverseSquaredGravityModel(), thrust); results_directory=joinpath(output_root, "with"))

results_without = run_simulation(config_without; return_results=true)
results_with = run_simulation(config_with; return_results=true)
results_with.table          # named state columns, also written to the CSV
results_with.files          # files this call wrote
results_with.configuration  # isolated configuration that actually ran

# --- 3. Verify the effect --------------------------------------------------------
function semi_major_axis(csv_path)
    df = CSV.read(csv_path, DataFrame)
    last = df[end, :]
    r = norm((last.sc1_pos_1, last.sc1_pos_2, last.sc1_pos_3))
    v = norm((last.sc1_vel_1, last.sc1_vel_2, last.sc1_vel_3))
    energy = v^2 / 2 - planet.μ / r          # specific orbital energy, J/kg
    return -planet.μ / (2 * energy), last.time, last.sc1_mass
end

a_without, t_end, mass_kg = semi_major_axis(joinpath(output_root, "without", "simulation_results.csv"))
a_with, _, _ = semi_major_axis(joinpath(output_root, "with", "simulation_results.csv"))
accel = thrust.force_n / mass_kg
v0 = sqrt(planet.μ / a0)
expected_da = 2 * a0^2 / planet.μ * accel * v0 * t_end
observed_da = a_with - a_without
println("mass used by the engine: $(mass_kg) kg, thrust acceleration: $(accel) m/s^2")
println("semi-major axis after one orbit without thrust: $(round(a_without / 1e3, digits=3)) km")
println("semi-major axis after one orbit with thrust:    $(round(a_with / 1e3, digits=3)) km")
println("expected rise: $(round(expected_da, digits=1)) m, observed rise: $(round(observed_da, digits=1)) m, ratio $(round(observed_da / expected_da, digits=3))")
```

## Choices in the script

**A type and one method.** An effector is a Julia `struct` that subtypes
`SpaceAGORA.AbstractForceTorqueModel` and carries its own parameters, here the
thrust magnitude. What makes it an effector is a method of the function
`SpaceAGORA.wrench`, which the engine calls with four arguments:

- `model`, your struct;
- `x`, a *state sample*: the spacecraft's position `x.pos_ii` and velocity
  `x.vel_ii` (meters and meters per second, in the planet-centered inertial
  frame), its current mass `x.mass_kg`, and, when attitude is propagated, its
  attitude quaternion `x.q_ib` and body angular velocity `x.ω_body` (both
  `nothing` otherwise);
- `env`, an *environment sample*: the planet model and, on request, the
  atmosphere, the planet-fixed frame, the Sun and third bodies (see the
  reference below);
- `t`, the elapsed simulation time in seconds.

The method returns a pair: the force in newtons in the inertial frame, and the
torque in newton-meters in the bus-fixed frame. Return a zero vector for the
part you do not model.

**Force, not acceleration.** The engine sums the forces of all effectors and
divides by the spacecraft's current mass itself. If your physics is naturally
an acceleration (a gravity model, a perturbation from a paper), multiply by
`x.mass_kg` before returning it, as the built-in gravity and solar-radiation
models do. An acceleration returned where a force is expected under-drives the
spacecraft by a factor equal to its mass.

**Frames.** The subscript `ii` on `pos_ii`, `vel_ii` and `force_ii` means
"inertial frame, inertial components": the non-rotating frame centered on the
planet in which the orbit is integrated. Torque is expressed in the frame
fixed to the spacecraft bus, the frame in which the spacecraft's inertia
tensor is given, so the two enter the attitude equations together. If your
force is known in the bus frame, rotate it into the inertial frame with the
attitude quaternion first. With `orientation_sim=false` no attitude is
propagated and `x.q_ib` is `nothing`. The spacecraft's configured geometry and
the fixed orientations of its links still exist, and the built-in aerodynamic
models use them in that mode.

**Which method the engine uses.** Before integration starts, the engine checks
that a method of `wrench` accepts the argument types `(YourType, StateSample,
EnvironmentSample, Float64)`. Annotating the model argument with your type is
what makes the method yours. The other three arguments can be left untyped,
annotated with those types, or annotated with wider types that contain them
(`t::Real` works). A method whose annotation does not accept the argument, say
`x::AbstractVector`, does not match. The run then stops before integration
with an `ArgumentError` that names the effector and says it has methods of
`wrench`, but none that accepts the arguments the run passes.

**Attaching it.** Effectors are passed as a tuple to `DynamicsModel`, or to
`make_example_config` as `dynamic_effectors`. Their order does not matter for
the physics, because the forces are summed. Gravity is an effector like any
other, so `InverseSquaredGravityModel()` is listed explicitly; a tuple passed
to `make_example_config` replaces its default of J2 gravity alone, so it must
include gravity too. An entry `scoped(effector, scope)` acts on some
spacecraft only: see
[Effectors on part of a constellation](constellations.md#Effectors-on-part-of-a-constellation).

**Two runs, not one.** The clearest way to see what an effector does is to run
the same configuration with and without it and compare the results files. The
example uses two-body gravity, so without thrust the semi-major axis stays
constant to the meter over one orbit, and every change in the run with thrust
comes from the effector.

## Expected result

On a fresh clone with no GRAM or SPICE data, the script prints:

```text
mass used by the engine: 520.0 kg, thrust acceleration: 9.999999999999999e-5 m/s^2
semi-major axis after one orbit without thrust: 6778.137 km
semi-major axis after one orbit with thrust:    6779.118 km
expected rise: 981.8 m, observed rise: 981.9 m, ratio 1.0
```

These numbers belong to this configuration; they are not universal values.
The mass is the bus plus its two panels. The expected rise comes from the rate
at which a tangential acceleration raises a near-circular orbit,
`da/dt = 2 a² a_t v / μ`, integrated over one period. The observed rise is
read from the final rows of the two results files, with the semi-major axis
computed from position and velocity through the orbital energy. The estimate
holds the initial semi-major axis and circular speed fixed while the simulated
orbit evolves, so the ratio is a diagnostic comparison with a first-order
estimate. It is not a pass/fail criterion or a bound on physical accuracy, and
the printed `1.0` is rounded to three decimal places, not a claim of exact
agreement. CI checks that both runs reach the end of the orbit, that the
engine uses the 520 kg vehicle and that the thrust raises the orbit; it does
not test the ratio.

For this configuration, a ratio near `1/520` suggests that the method returned
an acceleration instead of a force, and a negative rise suggests that the
force points against the velocity. Check the force's units and direction when
the comparison looks wrong. An `ArgumentError` before the run saying the
effector has no matching `wrench` method means the method's annotations do not
accept the arguments the run passes.

Each run writes the standard three files to its own directory under
`output/effector_walkthrough` beside the script; see
[Simulation Outputs](outputs.md). The thrust does not get a column of its own:
of the built-in models, only the aerodynamic ones write per-force
diagnostics.

## Variations

- **Point the thrust somewhere else.** Cross-track (`cross(pos, vel)`,
  normalized) changes the inclination; radial (`x.pos_ii / norm(x.pos_ii)`)
  changes the eccentricity. The two-run comparison still works.
- **Make it depend on time or position.** `t` is elapsed seconds, so a burn
  window is one `if`. `norm(x.pos_ii) - env.planet.Rp_e` is the height above a
  sphere of the planet's equatorial radius. That is not the ellipsoidal
  altitude saved in the results column `sc1_altitude` (see
  [The Integrated State](integrated_state.md)); choose the surface your
  condition is meant to refer to.
- **Use the atmosphere.** Declare what your effector needs and the engine
  samples it at every evaluation:

  ```julia
  SpaceAGORA.environment_requirements(::MyDrag) =
      SpaceAGORA.EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
  ```

  Then `env.atmosphere.rho_kg_m3` is the density at the current position and
  `env.planet_frame.vel_pp` is the velocity in the planet-fixed (rotating)
  frame. That is the air-relative velocity only when the wind is zero. A drag
  model must subtract the wind, `env.atmosphere.wind_pp`, after bringing it
  into the same frame: the field holds the east, north and up components, in
  that order, and the engine's thermal callback shows the conversion. With
  `NoAtmosphereModel()` the density is zero; `ExponentialAtmosphereModel(planet)`
  gives a simple analytic profile with zero wind.
- **Add a torque.** Return a nonzero second element and run with
  `orientation_sim=true`. The attitude columns `sc1_q_*` then appear in the
  results, and `x.q_ib` and `x.ω_body` are available inside `wrench`.
- **Give it state of its own.** Keep parameters in the struct's fields. If the
  effector must remember something between evaluations (an integrator of its
  own, a mode), use a control effector instead (see
  [Extensibility](../extensibility.md)). `wrench` is meant to be a pure
  function of its four arguments, because the solver also evaluates it at
  trial points that it later rejects.

## Reference: the interface

This section names the types and functions that the walkthrough used without
naming them. The field lists below are read from the source.

| Symbol | Role |
|---|---|
| `SpaceAGORA.AbstractForceTorqueModel` | Supertype of every dynamic effector. |
| `SpaceAGORA.wrench(model, x::StateSample, env::EnvironmentSample, t::Float64)` | The hook. Returns `(force_ii::SVector{3}, torque_body::SVector{3})`, SI units, force inertial, torque bus-fixed. Must be pure. |
| `SpaceAGORA.StateSample` | Fields `pos_ii`, `vel_ii` (m, m/s, inertial), `mass_kg`, `q_ib` (inertial-to-body quaternion, scalar-last `[x, y, z, w]`, or `nothing`), `ω_body` (rad/s, or `nothing`), `spacecraft` (the typed spacecraft model, for geometry and inertia), `spacecraft_index` (the spacecraft's position in the run, from 1; 0 for a sample built outside a run). |
| `SpaceAGORA.EnvironmentSample` | Fields `planet` (always), `planet_frame` (`alt_m`, `lat_rad`, `lon_rad`, `pos_pp`, `vel_pp`, `l_pi`: position and velocity in the planet-fixed frame and the inertial-to-planet-fixed rotation matrix), `atmosphere` (`rho_kg_m3`, `temperature_k`, `wind_pp` as east, north, up components, in that order), `solar` (`sun_pos_ii`), `third_bodies` (`names` as declared and `positions_ii`, each body's position relative to the run's planet in meters along the J2000 axes, in the declared order; a run with third bodies needs `SpiceEphemeridesModel()`). Each optional field is `nothing` unless requested. |
| `SpaceAGORA.environment_requirements(model)` | Returns `EffectorEnvironmentRequirements(planet_frame=, atmosphere=, solar=, third_body_names=)`; the default requests nothing. |
| `SpaceAGORA.is_gravitational(model)` | `true` for a model of gravity, whose force an accelerometer does not feel; the default is `false`, and the built-in gravity models return `true`. A gravity model of your own adds a method returning `true`. |
| `SpaceAGORA.solver_partition(model)` | `:explicit` (default) or `:implicit`; only matters under the `split_imex` solver mode. |
| `DynamicsModel([spacecraft], effectors)` | Where the effector tuple is attached. |

The rules every hook follows, such as the argument order, the frames, the units
and the signs, are on [Extension Conventions](../reference/extension_conventions.md),
and the [Extension Hook Reference](../generated/extension_hook_reference.md) lists
when the engine calls each hook. A model written for `calcForceTorque`, the
force hook that `wrench` replaced, does not run: SpaceAGORA 0.5.0 removed that
hook, and [Migrating Extension Hooks](../reference/migrating_extension_hooks.md)
shows how to move it.

The engine sums the forces of all effectors per spacecraft and divides by the
current mass to get the translational acceleration. It sums the torques in the
bus-fixed frame; they drive the attitude equations when `orientation_sim=true`
and are ignored otherwise. Solver modes other than the default are described
on [Solver Configuration](solver_configuration.md). The gravity-backbone hooks
listed on the Extensibility page are for effectors that take part in that
special mode; an ordinary force does not need them.

## Advanced examples with internal models

The two examples below, a gravity-gradient torque and replayed reaction-wheel
speeds, use models from the internal `SpaceAGORA.SimulationModel` namespace.
They are not part of the root-level interface the walkthrough above uses, and
the [public API policy](../public_api_policy.md) does not cover their names,
so they can change between releases.

### Add gravity-gradient torque without another gravity force

`SM.GravityGradientTorqueModel()` supplies the central body's gravity-gradient
attitude torque on its own, which is what you need when a separate model
already supplies the orbital gravity, such as `SM.GravitationalHarmonicsModel`.

These effectors apply inverse-squared gravity and its attitude torque as
separate terms:

```julia
using SpaceAGORA
SM = SpaceAGORA.SimulationModel

effectors = (SM.InverseSquaredGravityModel(), SM.GravityGradientTorqueModel())
# Pass effectors to SM.DynamicsModel and set orientation_sim=true.
```

The model returns zero force and the body-frame torque
`3μ/r³ (r̂ × J r̂)`. Here `J` is the spacecraft's body-frame inertia tensor,
`r̂` is the unit vector from the central body to the spacecraft, expressed in
the spacecraft body frame, and `μ` is the body's gravitational parameter. The
torque is in N m. The model uses the central-field calculation, so pairing it
with harmonic gravity does not add nonspherical corrections to the torque.
Without a propagated attitude it returns zero torque, and
`GravityGradientTorqueModel(gravity_gradient=false)` turns it off.

Enable attitude propagation with `orientation_sim=true`. Do not also set
`gravity_gradient=true` on an analytic gravity model in the same run, because
the torque would then be applied twice. When the run uses only an analytic
gravity model, that model's `gravity_gradient=true` option requests the same
central-field torque.

### Replaying reaction-wheel speeds

For recorded or prescribed wheel-speed histories, use
`SM.ReactionWheelMomentumModel`. It supplies the reaction torque the wheels
put on the spacecraft body, including the gyroscopic term, and needs no
atmosphere, GRAM, SPICE or additional package.

```julia
using SpaceAGORA
using LinearAlgebra
const SM = SpaceAGORA.SimulationModel

# Seconds, signed radians per second, body-frame axes, and kg m².
times = [0.0, 5.0, 10.0]
speeds = [0.0 0.0 0.0; 50.0 0.0 0.0; 100.0 0.0 0.0]
wheels = SM.ReactionWheelMomentumModel(times, speeds, Matrix{Float64}(I, 3, 3), 2e-5;
    spacecraft_index=1)
```

Include `wheels` in the configuration's `DynamicsModel` effector tuple and set
`orientation_sim=true`. `spacecraft_index` is the spacecraft's position in the
run's spacecraft vector, not its ID. Keep the same physical wheels out of the
commanded `rw_assembly`, or their momentum is counted twice. The model replays
a history; it does not enforce actuator torque or momentum limits. When the
table's clock differs from the run's, the `time_offset_s` keyword gives the
table time minus the run time.

`SM.wheel_momentum_body(wheels, t)` evaluates the wheel momentum at a run time.
Speeds follow a natural cubic spline through the samples and hold their end
values outside the recorded interval, where their derivative is zero, so the
derivative can jump at either end. Choose a table that covers the run unless
those end transitions are part of the forcing you intend.

With the default state isolation, the `wheels` object you pass keeps its
diagnostic state unchanged. With `isolate_state=false`, `wheels.state` records
the last solver stage evaluated, which may be a rejected trial, so use the
query functions at saved times for analysis. Calling `SpaceAGORA.wrench` on
the model directly is pure and applies to whatever sample you pass; selecting
the spacecraft by `spacecraft_index` and updating `wheels.state` happen only in
the engine's own evaluation.
