# How the ORACLE Laser Callbacks Run

Source references checked on 2026-09-14 for `7_SpaceAGORA.jl-main_v4`. Line numbers may shift after later edits. This is a source-based walkthrough, not a report of a successful simulation run.

## Main Idea

The application creates callback objects, attaches them to an ODE problem, and starts the solver. The solver library then checks their conditions and invokes their actions during integration.

**Registering a callback does not execute its action.** There is no direct call to the laser `affect!(integrator)` inside the application's `run_simulation` function.

## 1. Create the Laser Callbacks

In [oracle_laser_links.jl:211](../7_SpaceAGORA.jl-main_v4/II_examples/oracle_laser_links.jl#L211), the driver constructs the impulse callback; the next line constructs the scheduler callback:

```julia
impulse_cb = laser_impulse_callback(constellation, impulse_tracker, opts.mass_kg)
scheduler_cb = laser_link_scheduler_callback(constellation)
```

- Input: constellation, impulse tracker, and spacecraft mass for the impulse callback; constellation for the scheduler.
- Output: two `DiscreteCallback` objects that contain conditions and actions.
- At this stage: the callback objects are created, but no laser velocity kick is applied.

The builders return these objects at [laser_link_effectors.jl:424](../7_SpaceAGORA.jl-main_v4/I_src/1.3_dynamics/coupled/force_torque_models/laser_link_effectors.jl#L424) and [laser_link_effectors.jl:328](../7_SpaceAGORA.jl-main_v4/I_src/1.3_dynamics/coupled/force_torque_models/laser_link_effectors.jl#L328).

## 2. Pass Them to the Simulation

At [oracle_laser_links.jl:227](../7_SpaceAGORA.jl-main_v4/II_examples/oracle_laser_links.jl#L227), the driver supplies:

```julia
extra_callbacks = (impulse_cb, scheduler_cb)
```

This tuple is an argument to `run_simulation`. Its order places the impulse callback before the scheduler callback.

## 3. Assemble the Complete Callback Set

At [execution.jl:246](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/engine/execution.jl#L246), `run_simulation` calls `SimulationModel.get_callbacks(...)`, forwarding `extra_callbacks` at [execution.jl:252](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/engine/execution.jl#L252).

The implementation begins at [assembly.jl:201](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/callbacks/density_callbacks/assembly.jl#L201). It selects built-in callbacks for the configured run, such as impact detection, planet-frame updates, and result saving.

At [assembly.jl:260](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/callbacks/density_callbacks/assembly.jl#L260), it appends the extra callbacks, then returns the set at [assembly.jl:262](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/callbacks/density_callbacks/assembly.jl#L262):

```julia
callbacks = _append_callbacks(callbacks, extra_callbacks)
return CallbackSet(callbacks...)
```

Input: simulation configuration, effectors, saving inputs, and extra callbacks. Output: one `CallbackSet`. The actions still have not run.

## 4. Attach the Set to the ODE Problem

Inside `_build_typed_solver_problem`, [execution.jl:66](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/engine/execution.jl#L66) attaches callbacks when using a Jacobian prototype:

```julia
return ODEProblem(f, u0, tspan, p; callback=callbacks)
```

The equivalent branch without that prototype is at [execution.jl:68](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/engine/execution.jl#L68).

The resulting problem combines the dynamics function, initial state, time span, parameters, and callback set. Attaching callbacks still does not execute their actions.

## 5. Start the Solver

For the normal, non-checkpointed run, [execution.jl:433](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/engine/execution.jl#L433) calls `_solve_with_solver_policy(...)`. The checkpointed path calls it once per segment at [execution.jl:378](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/engine/execution.jl#L378).

For the uncached explicit solver path used by this example, the policy eventually calls `solve(prob, alg; ...)` at one of these locations:

- [solver_policy.jl:354](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/engine/solver_policy.jl#L354): no explicit maximum-iteration limit.
- [solver_policy.jl:361](../7_SpaceAGORA.jl-main_v4/I_src/5_simulation/engine/solver_policy.jl#L361): maximum-iteration limit supplied.

**This is where control enters the solver library.** Its internal integration loop invokes the registered callbacks; these application lines start that process rather than directly calling each callback action.

## 6. Initialize and Execute the Laser Actions

The scheduler has an initialization hook at [laser_link_effectors.jl:326](../7_SpaceAGORA.jl-main_v4/I_src/1.3_dynamics/coupled/force_torque_models/laser_link_effectors.jl#L326). It calls `choose_active_links!` during solver setup so the initial active-link selection is available.

Both laser callbacks use conditions that return `true`. During normal stepping, their discrete actions therefore run after accepted steps, not at every intermediate ODE derivative evaluation.

### Impulse Action

The nested `affect!(integrator)` starts at [laser_link_effectors.jl:379](../7_SpaceAGORA.jl-main_v4/I_src/1.3_dynamics/coupled/force_torque_models/laser_link_effectors.jl#L379).

- Input: the solver's live `integrator`, plus the constellation, tracker, and mass captured by the callback's closure.
- Action: read active endpoint pairs and their current states, calculate the elapsed-time laser kick, update both endpoints' velocities, and record cumulative RTN delta-V.
- Output: mutations to `integrator.u` and the tracker, rather than a returned trajectory.

### Scheduler Action

At [laser_link_effectors.jl:325](../7_SpaceAGORA.jl-main_v4/I_src/1.3_dynamics/coupled/force_torque_models/laser_link_effectors.jl#L325):

```julia
affect!(integrator) = choose_active_links!(constellation, integrator)
```

It reads the current spacecraft states and updates `constellation.possible_links` and `constellation.active_links` according to range and scheduling policy.

## Execution Order Summary

```text
Driver creates impulse_cb and scheduler_cb
    -> passes them to run_simulation as extra_callbacks
    -> get_callbacks combines them with the built-in callbacks
    -> ODEProblem stores the callback set
    -> solve starts integration and callback initialization
        -> scheduler initializes active links
        -> ODE evaluates derivatives and accepts a step
        -> impulse action applies kicks using the existing active links
        -> scheduler action chooses links for the next step
        -> repeat until termination
    -> solve returns the solution
```

This summary shows the laser-related sequence; built-in callbacks also execute according to their own types and conditions, and terminal events can end the run. The impulse-before-scheduler order comes from the supplied tuple, not from where the functions are defined in the source file.

In short, the application registers the callbacks and starts `solve`; the solver library invokes their actions while integration is running.