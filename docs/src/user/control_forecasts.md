# Candidate control forecasts

A controller can test a candidate command plan against a copy of the current
simulation. SpaceAGORA owns propagation, configured forces, event callbacks and
the simulation clock. The companion owns the candidate plan, its constraints
and the decision to accept it. A forecast does not install a plan in the live run.

This interface is experimental and deliberately limited. It supports one
spacecraft with one periodic controller, Cartesian dynamics, instantaneous
panel geometry, and the native `split_imex` or `auto_stiff` solver. Supported
forces are central/J2 gravity, gravity harmonics, third-body gravity, solar
radiation pressure and free-molecular aerodynamics. Atmosphere choices are
polynomial-fit or no atmosphere; thermal calculations use `MaxwellianHeat`.
Other configurations are rejected.

## Capture before applying the command

Call `capture_control_forecast(u, p, t, controller)` at the beginning of that
controller's `update_control!`, before changing its command, plan, geometry or
history. The controller must implement checkpoint hooks and report its previous
command time through `held_command_time`.

The current control tick remains pending in the copy. Propagation applies the
candidate at that tick exactly once. The original ephemeris epoch, elapsed time
and control-clock origin are preserved; the incoming state is not assigned a
new calendar start time. This is a specific pre-command boundary, not a general
checkpoint mechanism for partially completed callbacks.

## Evaluate in a dedicated process

Use `with_control_forecast_worker(f, snapshot; packages=())` to create a local
Julia worker for the search. Within `f`, call the supplied
`evaluate(candidate; duration_s, sample_period_s)` function. List companion
packages defining controller types as `Base.PkgId` values in `packages`.

Each evaluation starts from the same snapshot and a copy of the candidate.
Disable replanning in that candidate. The result contains saved rows, final
state, start/end times and the solver return code. The saved rows lie on the
run's output grid, counted in sample periods from the snapshot's clock origin,
plus the snapshot time and the forecast end. Inspect the return code and
required physical events; reaching the requested duration does not itself mean
the plan is feasible. The helper removes its worker on return or failure.

The lower-level `run_control_forecast` must also run in a dedicated process.
Do not run unrelated simulations there. Forecast calls are serialized within
the process because SPICE uses global state; nested forecasts are rejected.

## Limits

Linked producers, navigation, guidance, attitude dynamics, joints, mounted
devices, power devices, onboard data models, faults, external propagation and
interlinks are unsupported. Controller state must be ordinary copied data;
closures, remote references, runtime handles and arrays with nonstandard or
shared storage are rejected. This check is not a sandbox for untrusted Julia
code.

SPICE kernel files must be directly loaded and remain unchanged. Their hashes
and load order are retained, and worker caches are cleared before replay.
Meta-kernels and programmatic kernel-pool assignments are unsupported. Scoped
engine overrides are rejected. Resolved solver settings and `SPACEAGORA_`
environment values are frozen for the forecast.

The caller still needs a bounded search, an explicit rejection policy and
qualification against its intended mission. Forecast equivalence is separate
from controller accuracy, constraint satisfaction and mission completion.

API reference: [`ControlForecastSnapshot`](@ref),
[`capture_control_forecast`](@ref), [`run_control_forecast`](@ref),
[`with_control_forecast_worker`](@ref).
