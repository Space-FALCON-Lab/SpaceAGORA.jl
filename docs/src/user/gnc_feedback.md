# Feedback Controllers and Estimators

This guide is for models whose tick computes a command or an estimate from the
current state. To keep such a model in a package of its own, follow the
[rate-damper package recipe](feedback_package.md). A planner that proposes a
product for checking belongs in [Planning Algorithms](gnc_planning.md). Both
kinds follow the [shared house rules](gnc_algorithm.md#House-rules-for-both-profiles).

The examples are complete scripts in `examples/gnc_algorithm_env/`, a separate
project that uses only the public API. CI runs them, and a docs check fails when
a code block on this page differs from the matching region of its script. From
the SpaceAGORA checkout:

```bash
julia --startup-file=no --compiled-modules=existing --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/setup.jl
julia --startup-file=no --compiled-modules=existing --depwarn=error --threads=2 --project=examples/gnc_algorithm_env examples/gnc_algorithm_env/feedback_controller.jl
```

Run a whole script, imports included, before you change any of its excerpts. To
try the rate-damper blocks in a fresh Julia session, load these first:

```julia
using SpaceAGORA
using SpaceAGORA.Conformance: conformance_fixture, test_control_model
using StaticArrays
using Test
```

The estimator blocks also need `test_navigation_model` from
`SpaceAGORA.Conformance`.

## The feedback profile

A feedback controller has two parts. The tick,
`update_control!(model, u, p, t, spacecraft_index)`, runs once per control
period: it computes the command from the state `u`, limits it and holds it in
the model. The output hooks, `control_wrench(model, x, p, t)`,
`control_propellant_rate(model, x, p, t)` and
`control_wheel_torque(model, x, p, t)`, run at every evaluation of the equations
of motion. They return the held command and never run the law. The force is in
the inertial frame and the torque in the body frame; the propellant rate is
positive while propellant is used.

A controller that drives mounted hardware declares
`actuation_mode(model) = :mounted` and returns the hardware levels from
`actuator_levels(model, spacecraft_index)`. The engine then does not call
`control_wrench` or `control_propellant_rate` for it; see
[Adding an Actuator of Your Own](custom_actuator.md).

Estimators and guidance models tick the same way, with
`update_navigation!(model, u, p, t, spacecraft_index)` and
`update_guidance!(model, u, p, t, spacecraft_index)`, and keep their result for
the controller to read. When ticks fall at the same time, sensors sample first,
then navigation models tick, then control, then guidance. Two guidance models
tick before control instead: `AlgorithmGuidance`, the generic planning binding,
and the EDG guidance model when its modes include targeting. The control tick at
the same time then uses the product they install.

### Declarations

A model declares what it acts on and what it assumes. Every declaration has a
default that matches the engine's usual behavior, so declare only what differs.
The engine reads the declarations before the run and at each tick.

| Declaration | Default | Declare it when |
|---|---|---|
| `spacecraft_scope(model)` | `AllSpacecraft()`: ticked for every spacecraft | the model acts on some spacecraft only: `BoundSpacecraft(indices)`; `CoupledSpacecraft(primary, partners)` for a model ticked for `primary` that also reads `partners`; or a `SpacecraftSelection` (see [Effectors on part of a constellation](constellations.md#Effectors-on-part-of-a-constellation)) |
| `gnc_period_s(model)` | `nothing` | the law assumes a tick period, such as the step it was discretized with. The engine then refuses a run whose control rate differs |
| `skip_inactive(model)` | `true`: once a spacecraft is no longer active, after an impact or a touchdown, the engine stops ticking the model for it and calls `release_command!(model, t, spacecraft_index)` at each tick time instead | the model must keep ticking for an inactive spacecraft: return `false` |
| `held_command(model, spacecraft_index)`, `held_command_time(model, spacecraft_index)` | `nothing`, `NaN` | the model holds a command: return it as an immutable value, and return the time of the tick that produced it |
| `gnc_links(model)` | none | the model shares an object with a partner; see [Pairing through a link](#Pairing-through-a-link) |

A model that keeps the default `skip_inactive` supplies its off command, such as
zero thrust, by implementing `release_command!`. The fallback does nothing, so a
model without its own method keeps the command it last held.

During the run the engine calls a model only for the spacecraft its scope
names. That covers its ticks and `release_command!`, its output hooks
(`control_wrench`, `control_propellant_rate`, `control_wheel_torque`), its levels
(`actuator_levels`, `control_thruster_levels`), `register_tstops!` and
`atmosphere_transition!`. For a `CoupledSpacecraft` scope these calls go to the
primary only: the partners are read, not acted on. A few calls still reach every
spacecraft. Before the run the engine asks `touchdown_spec` about every
spacecraft, and `control_thruster_levels` too when it sizes the result columns,
and direct callers such as the conformance kit may call any of these hooks for
any spacecraft. So return zeros or `nothing` for a spacecraft outside the scope.

The engine fixes each model's scope when the run starts, after the run-start
initializers, and keeps it for the whole run. A `SpacecraftSelection` is
resolved at that point, with its `where` function called once per spacecraft.
The engine reads the scope again at each tick (an event-driven model's once per
accepted step) and stops the run if it has changed. Set a model's spacecraft in
the model itself, in `bind_spacecraft` or in an initializer, as the RPO planner
and controller do from the spacecraft IDs. A selection's `where` function is
compared by identity, so a model whose scope is a `SpacecraftSelection` should
keep the selection in a field and return that object. A selection built anew at
each call, with a new closure, is a different scope from the one the run took.

A control model that is not ticked periodically declares
`tick_schedule(model) = :event_driven` and supplies its own callbacks through
`gnc_callbacks(model, num_sats, args, period_s)`, as the burn scheduler does.
The engine adds those callbacks to the run unchanged. They run outside the
wrapper that handles periodic ticks, so each callback has to do that wrapper's
work itself:

- act only on the spacecraft that `spacecraft_scope` names;
- follow the model's `skip_inactive` declaration. With `true`, leave a
  spacecraft alone once `p.is_active` says it is inactive, and release its
  command as `release_command!` does; with `false`, keep acting as for an active
  spacecraft;
- time its events on the run's clock, not from the start of a solve, because a
  checkpointed run is solved in segments and can resume from a file;
- keep any state it needs in the model, where `checkpoint_state` can save it.
  The engine does not save a callback's own state, and the
  [checkpoint table](rpo_planner_pilot.md#Checkpoints) leaves custom callbacks
  outside its guarantee.

A periodic model can add solver stops after each tick with
`register_tstops!(integrator, model, spacecraft_index)`, which the engine calls
for each spacecraft in the model's scope, and can record when its clock starts
with `gnc_clock_started!(model, period_s, origin_s)`. Most models need neither.
[Migrating Extension Hooks](../reference/migrating_extension_hooks.md) lists how
these rules changed between releases.

### Worked example A: a rate damper

`RateDamper` damps a spacecraft's body rate. Its law is a torque proportional to
the rate, limited on each axis. The script is
`examples/gnc_algorithm_env/feedback_controller.jl`.

```julia
Base.@kwdef mutable struct RateDamper <: SpaceAGORA.AbstractControlEffectorModel
    spacecraft_index::Int = 1           # the spacecraft it acts on, its position in the run
    gain_n_m_s::Float64 = 2.0           # torque per unit body rate, N·m per rad/s
    max_torque_n_m::Float64 = 1.0       # torque limit on each body axis, N·m
    period_s::Float64 = 1.0             # the tick period the law assumes, s
    # Run state: the held command and the time of the tick that set it.
    torque_body::SVector{3, Float64} = zero(SVector{3, Float64})
    tick_s::Float64 = NaN
end

# Declarations: what the controller acts on and what it assumes.
SpaceAGORA.spacecraft_scope(c::RateDamper) = BoundSpacecraft((c.spacecraft_index,))
SpaceAGORA.gnc_period_s(c::RateDamper) = c.period_s
SpaceAGORA.held_command(c::RateDamper, spacecraft_index) = c.torque_body
SpaceAGORA.held_command_time(c::RateDamper, spacecraft_index) = c.tick_s
# While its spacecraft is inactive (after an impact or a touchdown), the engine calls
# this at each tick in place of the tick (`skip_inactive`): the off command.
function SpaceAGORA.release_command!(c::RateDamper, t, spacecraft_index)
    c.torque_body, c.tick_s = zero(SVector{3, Float64}), t
    return nothing
end

# At each run's start: fresh run state, so two runs of one configuration start equal.
function SpaceAGORA.initialize_control!(c::RateDamper, u, p, t)
    c.torque_body, c.tick_s = zero(SVector{3, Float64}), NaN
    return nothing
end

# The tick (shape T), once per control period: run the law, limit it, hold the command.
function SpaceAGORA.update_control!(c::RateDamper, u, p, t, spacecraft_index)
    ω_body = SVector{3, Float64}(u.sc[spacecraft_index].ω)
    c.torque_body = clamp.(-c.gain_n_m_s .* ω_body, -c.max_torque_n_m, c.max_torque_n_m)
    c.tick_s = t
    return nothing
end

# The output (shape R), at every evaluation of the equations of motion: the held
# command, never the law. The engine calls it for the spacecraft in the scope
# only; direct callers may ask about others, so it checks the index.
function SpaceAGORA.control_wrench(c::RateDamper, x::StateSample, p, t)
    zero3 = zero(SVector{3, Float64})
    x.spacecraft_index == c.spacecraft_index || return zero3, zero3
    return zero3, c.torque_body             # force_ii (N), torque_body (N·m)
end

# A checkpoint must carry the held output between ticks, including its timestamp.
# These values are immutable; restoring them cannot share mutable buffers.
SpaceAGORA.checkpoint_state(c::RateDamper) =
    (torque_body=c.torque_body, tick_s=c.tick_s)
function SpaceAGORA.restore_checkpoint_state!(c::RateDamper, data)
    c.torque_body, c.tick_s = data.torque_body, data.tick_s
    return nothing
end

# Declared on a spacecraft: a copy for that spacecraft, with fresh run state.
SpaceAGORA.bind_spacecraft(c::RateDamper, spacecraft_index::Int) = RateDamper(;
    spacecraft_index, gain_n_m_s=c.gain_n_m_s, max_torque_n_m=c.max_torque_n_m, period_s=c.period_s)
```

The run below uses a cube far from any planet, spinning about one axis, with
nothing else acting on it:

```julia
# A 500 kg cube far from any planet, spinning about body z; only the controller acts.
root = Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
ω0 = 0.05                                    # initial rate, rad/s
ic = CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0);
                               q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, ω0))
spacecraft = SpacecraftModel(links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
                             initial_condition=ic, id=1)
damper = RateDamper()
mission_s = 120.0
planet = make_no_gram_planet(:earth)
config = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=mission_s, number_of_orbits=1,
        keplerian=false, orientation_sim=true, num_steps_to_save=120, data_rate=1.0),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], ()),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
    control_model=ControlModel(control_effectors=(damper,), control_rates=[damper.period_s]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

solution = run_simulation(config; return_solution=true)
```

The result can be checked by hand. On a body whose inertia is the same about
every axis, a torque held for one control period changes the rate by a fixed
fraction:

```julia
# By hand: the first tick is one period after the start, and each held torque
# -k ω(t_n) acts for one period on a cube (no gyroscopic torque), so after the
# ticks at 1, 2, ..., N s the rate is ω0 (1 - k Δt / I)^N.
inertia = root.inertia[3, 3]
ticks = floor(Int, mission_s / damper.period_s) - 1        # ticks before the end
expected = ω0 * (1 - damper.gain_n_m_s * damper.period_s / inertia)^ticks
final = solution.u[end].sc[1].ω[3]
println("body rate after $(mission_s) s: $(final) rad/s, by hand $(expected) rad/s")
@test final ≈ expected rtol=1e-6
```

The law runs in the tick only. `update_control!` reads the rate from
`u.sc[spacecraft_index]` (`ω` is the body rate in rad/s), limits the torque and
stores it. `control_wrench` returns the stored torque at every evaluation, so
the vehicle sees a zero-order hold of the command, as it would from flight
software. The first tick comes one control period after the start; until then
the initial held command acts, which here is zero.

The declarations state the model's assumptions. Because of `gnc_period_s`, the
engine refuses a run whose rate in `ControlModel(control_rates=...)` is not the
period the gain was tuned for. `held_command` and `held_command_time` let code
outside the model read the held command and the time it was set.

`initialize_control!` sets the held torque back to zero when a run starts.
`run_simulation` normally runs a copy of the configuration, but with
`isolate_state=false` the same model object runs again, and it has to start as
it did the first time.

To use the model, add it to `ControlModel(control_effectors=(...),
control_rates=[...])`, or declare it on a spacecraft with
`SpacecraftModel(...; control=ControlModel(...))`. In the second case the engine
calls `bind_spacecraft(model, spacecraft_index)`, which returns a copy for that
spacecraft with fresh run state.

The spacecraft index, gain, torque limit and period are configuration. The held
torque and the tick time are run state that the model owns: initialization
clears them, so configure the next run by changing the gain or the limit, not
the last output. The checkpoint hooks carry that output across a restart
between ticks. A model whose run also writes values a user can configure
follows the [conditional reset rule](gnc_algorithm.md#House-rules-for-both-profiles).

### The lifecycle hooks

Besides its tick, a guidance, navigation or control model can implement hooks
that the engine calls around the ticks. Each default does nothing.

| Hook | When the engine calls it | Implement it to |
|---|---|---|
| `preflight_guidance(model, args; isolate_state)`, `preflight_control(model, args; isolate_state)` | before the run, on the caller's configuration | refuse a configuration the model cannot run, with an `ArgumentError`; the hook must change nothing |
| `initialize_navigation!(model, u, p, t)`, `initialize_guidance!(model, u, p, t)`, `initialize_control!(model, u, p, t)` | at the start of each run, after each link's `link_initialize!`: first for the navigation models, then guidance, then control | reset the model's run state so that two runs of one configuration with `isolate_state=false` start equal, while keeping a change a user made between the runs (see the house rules) |
| `requires_atmosphere_events(model)`, `atmosphere_transition!(model, u, p, t, spacecraft_index, inside)` | `atmosphere_transition!`, for guidance and control models, when a spacecraft in the model's scope enters (`inside` is `true`) or leaves the atmosphere. If any model's `requires_atmosphere_events` returns `true`, the engine tracks these crossings in that run; otherwise only in runs that need them for other reasons, such as different solver tolerances inside the atmosphere | react to the entry or exit |
| `before_reference_control!(producer, consumer, u, p, t)` | before each periodic tick of a control model, for each model that produces a link the control model consumes | bring what the producer writes up to date, as the RPO planner binding checks that its reference is still current |

What a run leaves in a model stays readable after the run, because the reset
happens only when the next run starts. `checkpoint_state` and
`restore_checkpoint_state!`, shown in example A, carry run state through a
checkpoint; see [Run state in checkpoint files](../extensibility.md#Run-state-in-checkpoint-files).
A custom sensor resets its own state in `initialize_sensor!`
([Extensibility](../extensibility.md#Add-a-sensor)).

### Worked example B: an estimator and a link

An estimator ticks and holds a value as a controller does, but its output is an
estimate for other models to read. `GyroRateEstimator` reads a gyro and holds
the measured rate. `EstimatedRateDamper` is example A's controller, reading that
estimate instead of the true state. The script is
`examples/gnc_algorithm_env/estimator.jl`.

#### Pairing through a link

The two models share one object, a link. Each declares its part with
`gnc_links(model)`: `Produces(link)` for the model that writes it and
`Consumes(link)` for a model that reads it. The engine groups the run's models
by link, and the grouping survives the deep copy of the configuration. Before
the run it checks the number of producers and consumers of each link against
`link_cardinality(link)` and calls `check_link(link, producers, consumers, args)`,
where you can refuse settings that disagree. At the start of every run it calls
`link_initialize!(link, producers, consumers, u, p, t)` before any model's own
initializer, so each run starts with the link in the state `link_initialize!`
gives it. Here that is an empty estimate.

```julia
# The estimate the navigation model writes and the controller reads: one shared object.
mutable struct RateEstimate
    ω_body::SVector{3, Float64}     # estimated body rate, rad/s
    time_s::Float64                 # time of the measurement behind it, s into the run
end
RateEstimate() = RateEstimate(zero(SVector{3, Float64}), NaN)

# One estimator writes it; any number of controllers read it. Each run starts empty.
SpaceAGORA.link_cardinality(::RateEstimate) = (producers = 1:1, consumers = 0:typemax(Int))
function SpaceAGORA.link_initialize!(e::RateEstimate, producers, consumers, u, p, t)
    e.ω_body, e.time_s = zero(SVector{3, Float64}), NaN
    return nothing
end
```

The estimator reads the latest measurement of the gyro mounted on its
spacecraft with `sensor_measurement` and `sensor_measurement_time`:

```julia
struct GyroRateEstimator <: SpaceAGORA.AbstractNavigationModel
    spacecraft_index::Int
    gyro::Int                       # the gyro's position among the spacecraft's sensors
    estimate::RateEstimate
end
SpaceAGORA.spacecraft_scope(n::GyroRateEstimator) = BoundSpacecraft((n.spacecraft_index,))
SpaceAGORA.gnc_links(n::GyroRateEstimator) = (Produces(n.estimate),)

# The tick: read the latest measurement and hold it as the estimate. Sensors
# sample before navigation ticks at the same time.
function SpaceAGORA.update_navigation!(n::GyroRateEstimator, u, p, t, spacecraft_index)
    spacecraft = p.args.dynamics_model.spacecraft[spacecraft_index]
    n.estimate.ω_body = SVector{3, Float64}(sensor_measurement(spacecraft, n.gyro))
    n.estimate.time_s = sensor_measurement_time(spacecraft, n.gyro)
    return nothing
end

# The producer saves the estimate. Restore its values into the existing link so
# every consumer still reads the same shared object.
SpaceAGORA.checkpoint_state(n::GyroRateEstimator) =
    (ω_body=n.estimate.ω_body, time_s=n.estimate.time_s)
function SpaceAGORA.restore_checkpoint_state!(n::GyroRateEstimator, data)
    n.estimate.ω_body, n.estimate.time_s = data.ω_body, data.time_s
    return nothing
end
```

The controller declares that it consumes the estimate:

```julia
# The rate damper of example A, reading the estimate instead of the true state.
Base.@kwdef mutable struct EstimatedRateDamper <: SpaceAGORA.AbstractControlEffectorModel
    spacecraft_index::Int = 1
    gain_n_m_s::Float64 = 2.0
    max_torque_n_m::Float64 = 1.0
    period_s::Float64 = 1.0
    estimate::RateEstimate
    torque_body::SVector{3, Float64} = zero(SVector{3, Float64})
    tick_s::Float64 = NaN
end
SpaceAGORA.spacecraft_scope(c::EstimatedRateDamper) = BoundSpacecraft((c.spacecraft_index,))
SpaceAGORA.gnc_period_s(c::EstimatedRateDamper) = c.period_s
SpaceAGORA.gnc_links(c::EstimatedRateDamper) = (Consumes(c.estimate),)
SpaceAGORA.held_command(c::EstimatedRateDamper, spacecraft_index) = c.torque_body
SpaceAGORA.held_command_time(c::EstimatedRateDamper, spacecraft_index) = c.tick_s
function SpaceAGORA.release_command!(c::EstimatedRateDamper, t, spacecraft_index)
    c.torque_body, c.tick_s = zero(SVector{3, Float64}), t
    return nothing
end
function SpaceAGORA.initialize_control!(c::EstimatedRateDamper, u, p, t)
    c.torque_body, c.tick_s = zero(SVector{3, Float64}), NaN
    return nothing
end

function SpaceAGORA.update_control!(c::EstimatedRateDamper, u, p, t, spacecraft_index)
    c.torque_body = clamp.(-c.gain_n_m_s .* c.estimate.ω_body, -c.max_torque_n_m, c.max_torque_n_m)
    c.tick_s = t
    return nothing
end

function SpaceAGORA.control_wrench(c::EstimatedRateDamper, x::StateSample, p, t)
    zero3 = zero(SVector{3, Float64})
    x.spacecraft_index == c.spacecraft_index || return zero3, zero3
    return zero3, c.torque_body
end

# The controller saves its held command; the producer owns the shared estimate.
SpaceAGORA.checkpoint_state(c::EstimatedRateDamper) =
    (torque_body=c.torque_body, tick_s=c.tick_s)
function SpaceAGORA.restore_checkpoint_state!(c::EstimatedRateDamper, data)
    c.torque_body, c.tick_s = data.torque_body, data.tick_s
    return nothing
end
```

With an exact gyro sampled at the control ticks, the pair reproduces example
A's result. At each tick time the gyro samples first, then the estimator ticks,
then the controller.

```julia
root = Link(root=true, m=500.0, dims=MVector(1.0, 1.0, 1.0))
ω0 = 0.05
ic = CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0);
                               q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, ω0))
spacecraft = SpacecraftModel(links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
                             initial_condition=ic, id=1)
add_sensor!(spacecraft, root, RateGyro(period_s=1.0))      # no noise and no bias: exact rates
estimate = RateEstimate()
estimator = GyroRateEstimator(1, 1, estimate)
damper = EstimatedRateDamper(estimate=estimate)
mission_s = 120.0
planet = make_no_gram_planet(:earth)
estimator_configuration(mission_time_s) = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime, mission_time=mission_time_s, number_of_orbits=1,
        keplerian=false, orientation_sim=true, num_steps_to_save=120, data_rate=1.0),
    environment_model=EnvironmentModel(
        planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false, wind=false),
    dynamics_model=DynamicsModel([spacecraft], ()),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(navigation_effectors=(estimator,), navigation_rates=[1.0]),
    control_model=ControlModel(control_effectors=(damper,), control_rates=[damper.period_s]),
    initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
    integration_tolerances=IntegrationTolerances())

config = estimator_configuration(mission_s)
solution = run_simulation(config; return_solution=true)

# With an exact gyro sampled at the control ticks, the controller sees the true
# rate at each tick, and the run matches example A's closed form.
inertia = root.inertia[3, 3]
ticks = floor(Int, mission_s / damper.period_s) - 1
expected = ω0 * (1 - damper.gain_n_m_s * damper.period_s / inertia)^ticks
final = solution.u[end].sc[1].ω[3]
println("body rate after $(mission_s) s with the estimate: $(final) rad/s, by hand $(expected) rad/s")
@test final ≈ expected rtol=1e-6
```

A gyro with noise or a bias gives an estimate that differs from the true rate;
the sensor kinds are described under [Extensibility](../extensibility.md#Add-a-sensor).

## Test your algorithm

`SpaceAGORA.Conformance` runs the engine's own machinery on your model. For a
guidance, navigation or control model it runs the engine's pre-run hook check
and then runs the fixture several ways and compares the results: repeated in
sequence; after a `deepcopy`; with Julia's global random stream seeded
differently; and, with two or more threads, twice at the same time. Two runs of
one configuration with `isolate_state=false` test the reset at a run's start.
It also checks the copies `bind_spacecraft` makes and a checkpoint round trip.
A control model is checked as well for pure output hooks, a propellant rate
that is never negative, and zero outputs for spacecraft outside its scope.
Each `check_*` function returns a report, and with `using Test` loaded each
`test_*` function records the same checks as tests. Each example script ends
with them:

```julia
# The engine's own checks on the controller: hooks, reproducibility, copies, the
# reset at each run's start, output-hook purity and the copies `bind_spacecraft` makes.
@testset "RateDamper conformance" begin
    # Use the spinning case above: a zero-rate fixture would not expose a lost torque.
    report = test_control_model(damper;
        fixture=conformance_fixture(configuration=config))
    # Require the checkpoint check itself to pass on this nonzero-spin fixture.
    @test only(filter(check -> check.id === :checkpoint_round_trip,
                      report.checks)).status === :pass
end
```

The estimator demonstration runs for 120 seconds. Its conformance checks use a
separate 13-second configuration with the same spinning spacecraft, shared
estimate and mounted gyro. The checkpoint round trip writes its checkpoint
halfway through, at 6.5 seconds, between two ticks, so the resumed run has to
carry the estimate and the held command through to the next update:

```julia
@testset "Estimator and controller conformance" begin
    # Keep the student demonstration at 120 s. This separate short configuration
    # checkpoints at 6.5 s, between ticks, while retaining the same nonzero-spin pair.
    conformance_config = estimator_configuration(13.0)
    fixture = conformance_fixture(configuration=conformance_config)
    reports = (test_navigation_model(estimator; fixture=fixture),
               test_control_model(damper; fixture=fixture))
    for report in reports
        @test only(filter(check -> check.id === :checkpoint_round_trip,
                          report.checks)).status === :pass
    end
end
```

Run the checks with two or more threads so that concurrent runs are compared
with sequential ones; with one thread that check is skipped. A pass certifies
interface behavior only. The method, its parameters and any claims made from it
remain your responsibility.

To keep your algorithm in a repository of your own, tested against a recorded
SpaceAGORA commit, start from a contributor starter; see
[Building Your Own Research Repository](research_repository.md).
