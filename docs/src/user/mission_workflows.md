# Define, run and inspect a mission

A mission starts with the same `SimulationConfiguration` used by every
SpaceAGORA run. A `MissionPlan` gives that configuration a name and a reusable
builder, and can collect operational segments, descriptive phases and event
callback factories. All segments share one integration and one continuous
physical state.

The workflow has three steps: define the configuration, run it into a fresh
output directory, then inspect the returned data or open the saved viewer
payload. The viewer can run without the simulation package.

![One mission run on a timeline. One integration carries a continuous state. Operational segments initial_coast, observation and final_coast change commands at transitions at 120 and 180 seconds, where the exit action runs before the entry action. Descriptive phases only name time windows. An event callback, built fresh for the run, can end it before the mission horizon. Afterwards run.history, the saved mission file and the viewer show the segments that actually ran.](../assets/figures/mission_workflow_timeline.svg)

The diagram uses the three-segment example and the two phase windows below;
the stop event is placed for illustration.

## Build the configuration inside a function

Put spacecraft, controllers and other mutable model state inside the builder,
so that each call of the plan starts from a fresh configuration. The engine's
normal state isolation also stays in place, so the returned configuration
holds the objects that actually ran.

This complete example uses analytic Earth gravity and needs no atmosphere or
SPICE data. All numerical choices are made in the configuration builder, and
the mission wrapper keeps them.

```julia
using SpaceAGORA

function build_coast(; duration_s=300.0)
    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=InitialCondition(ra=7.0e6, rp=7.0e6, i=45.0, ω=0.0, Ω=0.0, ν=0.0),
        prop_mass=0.0, id=1)
    return make_example_config(
        planet=planet, spacecraft=spacecraft, mission_time=Float64(duration_s),
        initial_time=InitialTime(year=2024, month=1, day=1),
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        orientation_sim=false, keplerian=true, verbose=false,
        solver_config=SolverConfig(solver_mode=:tsit5, parallel=false))
end

mission = MissionPlan("Earth coast", build_coast)
run = run_mission(mission;
    parameters=(duration_s=300.0,), output_dir="output/earth-coast-001")

run.results.table                # the saved samples, already available in memory
run.results.configuration        # the isolated configuration used by the engine
run.files                        # files written by this run
run.payload_path                 # portable viewer JSON
```

`output_dir` must not exist or must be empty: the workflow refuses a populated
directory rather than replace an earlier result, so give each run of the plan
its own destination. The builder must return a configuration with a finite,
positive `mission_time`. `run_mission` also refuses checkpoint writing and
resume; use `run_simulation` directly for those.

The wrapper changes only the output settings: it turns results on, sets the
output directory and sets the scene-output flag. The builder's dynamics, GNC,
epoch, solver, tolerances and output sampling are kept. The
`MissionConfiguration` still sets the overall horizon or the orbit-count
termination. A final operational segment that completes can end the run
earlier. Descriptive phases do not change the run, and an event shortens it
only when its callback implements a stop.

## Inspect data or open the viewer later

By default the workflow saves the normal result files, the scene sidecar and
`simulation_results_viewer.json`. The JSON holds the data the renderer needs
and can be viewed after the simulation process has ended. To keep it small,
the default payload leaves out surface textures; a payload with textures can
be exported later from the saved bundle without running the simulation again.

Install the standalone viewer as the [visualization guide](visualization.md)
describes, then use:

```julia
using SpaceAGORAVisualizer
SpaceAGORAVisualizer.visualize("output/earth-coast-001")
```

You can move or copy that run directory and pass its new path to the viewer.
The `MissionRun` object keeps the original paths, and it is not needed to
view the moved payload. The viewer does not need the configuration builder,
SpaceAGORA, GRAM or MuJoCo to display the saved data.

With both packages loaded, `SpaceAGORAVisualizer.visualize(run)` opens the
saved payload directly. To build the HTML viewer page with other export
options, call `export_visualization(run)`. It reads the saved scene and result
bundle and does not propagate again. To save only the result tables and their
bundle, set `save_scene=false` and `save_payload=false` together. A payload
requires scene output, and scene output requires the result bundle, which is
written by default (`SPACEAGORA_SAVE_BUNDLE=0` turns it off). The wrapper
rejects an incompatible output request before it creates the destination.

`run.results.solution` is `nothing` unless you pass `return_solution=true`.
The returned table follows the configured output sampling, so not every
accepted solver step or event instant appears as a row.

## Divide a mission into operational segments

A `MissionSegment` defines what happens between two transitions. The first
segment starts with the simulation. When its exit condition is met, SpaceAGORA
calls its exit action, enters the next segment and continues from the same
position, velocity, attitude and mass. It does not build a new spacecraft or
restart its controllers.

This sequence uses the coast configuration above:

```julia
mission = MissionPlan("Three coast segments", build_coast;
    segments=(
        MissionSegment(:initial_coast; until=MissionDeadline(120.0)),
        MissionSegment(:observation; until=MissionDuration(60.0)),
        MissionSegment(:final_coast),
    ))
run = run_mission(mission; output_dir="output/three-segments-001")
run.history.segments
run.history.transitions
run.history_path
```

The transitions occur at 120 and 180 seconds. `MissionDeadline` is an elapsed
time measured from the start of the simulation; `MissionDuration` is measured
from entry into its segment. Both boundaries are registered as solver stops,
so a transition does not wait for the next output sample or controller tick.
This example coasts throughout, because it supplies no command actions.

Only the final segment may omit an exit condition; it then lasts until the
mission horizon or another engine stop. A final segment with an exit condition
ends the simulation when that condition is met. Segment names must be unique,
durations must be positive, and a deadline must lie after the segment's
entry. A deadline or duration longer than `mission_time` is refused before
the run starts, and a boundary that would fall past the horizon from the
segment's actual entry time is refused when that segment is entered.

### Change a command at a transition

Each segment can provide a `build_actions` function. It receives the isolated
configuration that will actually run, after spacecraft GNC has been composed,
and returns a named tuple with `on_enter` and `on_exit` functions (either may
be left out). Capture this run's controller or reference in that function,
not a controller from an earlier run. Both actions receive a context with
`time_s`, `configuration`, a copy of the solver `state`, and immutable
numerical spacecraft samples in `spacecraft`.

For example, a controller with a mutable `commanded_thrust_N` field can be
commanded like this:

```julia
burn = MissionSegment(:maneuver; until=MissionDuration(20.0),
    build_actions=config -> begin
        controller = only(config.control_model.control_effectors)
        (
            on_enter=ctx -> (controller.commanded_thrust_N = 2.0),
            on_exit=ctx -> (controller.commanded_thrust_N = 0.0),
        )
    end)
```

That field belongs to this example controller's interface; SpaceAGORA
controllers in general do not have it. The controller must turn its command
into force, torque, mass flow and any displayed thruster firing. An exit
action must clear every held actuator command that should stop. A change to a
target reference takes effect when its controller next evaluates; a change to
a held actuator command can take effect at the boundary itself. The segment
executor does not disable a controller or reinitialize a planner on its own.

At a segment transition, the exit and entry actions run before the usual
discrete control and saving callbacks at the same time, and at a boundary
between two segments the exit action runs first. The last active segment's
exit action also runs when the engine reaches the horizon or stops early,
after the final normal output snapshot. Use these actions to change commands
and references, not to replace the vehicle, environment, solver or state
layout. Changing the copied `state` does not modify the integrated state, and
a segment boundary never creates an instantaneous velocity change.

### Transition on a physical condition

A `MissionCondition` is a scalar function of the current context. For a
configuration with a chaser and a target, this condition ends the segment when
their separation falls through 10 meters:

```julia
using LinearAlgebra: norm

approach = MissionSegment(:approach;
    until=MissionCondition(
        ctx -> norm(ctx.spacecraft[2].pos_ii - ctx.spacecraft[1].pos_ii) - 10.0;
        direction=:falling))
```

Positions and velocities in these samples are inertial SI values. The samples
also hold `mass_kg`, and the attitude and angular velocity when they are
available. `direction` is `:rising` by default. A falling condition must start
strictly positive and a rising one strictly negative; a condition already
satisfied at entry is rejected rather than silently skipping the segment. The
condition must return a finite real number, and must be deterministic and
free of side effects. The solver locates its zero crossing, so saving more
frames does not make the transition more precise.

The configured models and state dimensions stay fixed through a segmented run.
Segmented runs support the `:tsit5`, `:auto_stiff`, `:rodas5p` and `:dp8`
solver modes; other solver modes and external propagators, including MuJoCo,
are rejected, and the workflow never substitutes a different solver.
Checkpoint and resume are not supported. A servicing sequence can change
commands for models that are already configured, but it cannot transfer
physics ownership between SpaceAGORA and an external simulator.

### Inspect what actually ran

`run.history` holds the actual segment intervals and transition times, with a
reason for each transition. An early stop closes the active segment, and later
segments are absent. The same history is saved as
`simulation_results_mission.json` and included in the portable viewer
payload. The result manifest links the history to that run's scene and
trajectory, with content checks so that an old record cannot be attached to
new results; move the complete output folder, including its manifest, to keep
this link. A history file that is not bound this way is never attached
automatically. The history is also saved when scene and payload output are
disabled, without being associated with any older scene in that folder.

The standalone viewer can select a segment or transition within the saved
playback range and jump to its time, and its current-segment label follows
playback and scrubbing. A boundary outside that range stays listed with its
actual time, but its jump control is disabled because no pose was saved
there. Boundaries keep their recorded times even when trajectory rows are
sampled more sparsely; the viewer interpolates between saved states. Export
keeps the boundary rows that exist and refuses a frame budget too small to
hold them. This history describes the operational segments that ran;
`MissionPhase` windows and declared `MissionEvent` names are separate
metadata.

## Name phases without changing the equations

A `MissionPhase` names an elapsed-time window, for interpreting one run:

```julia
mission = MissionPlan("Earth coast", build_coast;
    phases=(
        MissionPhase(:initial_coast; window_s=(0.0, 120.0)),
        MissionPhase(:later_coast; window_s=(120.0, 300.0)),
    ))
```

These windows are observation metadata. They do not switch force models,
activate a controller, reset state, add solver stops or split the run into
separate solves. The phases are available in `run.phases`, including when an
existing stop condition ends the run before a later window is reached. Windows
must be finite, ordered, nonoverlapping and within the configured mission
time; gaps between them are allowed.

The engine has its own atmosphere entry and exit behavior, solver regimes and
controller lifecycle rules, and a phase name does not override them. Use
operational segments for explicit command transitions. Replacing physical
models during a run is not supported.

## Give events fresh callback state

A `MissionEvent` stores a function that builds a solver callback. SpaceAGORA
calls that factory once for each run, after epoch alignment, spacecraft GNC
composition and state isolation, with the configuration the engine will use.
Mutable callback state belongs inside the factory, so running the same plan
twice does not reuse counters or flags from the first run.

This example stops at the first accepted step at or after a requested elapsed
time, using a discrete callback:

```julia
using OrdinaryDiffEq: DiscreteCallback, terminate!

stop_event = MissionEvent(:coast_limit, config -> begin
    fired = Ref(false)
    limit_s = min(120.0, config.mission_configuration.mission_time)
    condition(u, t, integrator) = !fired[] && t >= limit_s
    function stop!(integrator)
        fired[] = true
        terminate!(integrator)
    end
    DiscreteCallback(condition, stop!; save_positions=(false, false))
end)

mission = MissionPlan("Earth coast with stop", build_coast; events=(stop_event,))
run = run_mission(mission; output_dir="output/earth-coast-stop-001",
    return_solution=true)
run.results.solution.t[end]
```

This is an accepted-step condition, not an exact-time event: the stop can come
later than `limit_s`, and the last saved table row can come before the stop.
See [Stopping on a Condition](stop_conditions.md) for continuous and discrete
conditions and their timing. A callback that registers exact time stops keeps
that behavior inside a `MissionEvent`. A state anchor still replaces the
specified state; wrapping it in a `MissionEvent` does not turn it into a
physical burn.

Factories must return a new callback or callback set for each call, and must
build it without performing I/O or changing the configuration. Their callback
closures may refer to the isolated controller and model objects in `config`;
a callback shared with another run defeats this. The built-in callbacks keep
their order, and mission callbacks go where `run_simulation` puts extra
callbacks.

`run.event_names` records the events that were configured, not which
conditions fired; event diagnostics stay with the callback or controller that
owns them. Invalid factories are rejected before SpaceAGORA creates run
outputs. A failure during propagation follows the engine's usual behavior:
partial outputs can remain, and the exception is rethrown rather than
returned as a successful `MissionRun`.
