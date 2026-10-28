# Solver Configuration

This page explains the integrator settings: the tolerances and step limits in
`IntegrationTolerances`, which part of the state each one governs, and the
solver modes that `SolverConfig` selects. Turn to it when a run is not accurate
enough, or takes too long, with the defaults.

The quickstart example runs with the tolerances `make_example_config` sets,
listed under [Practical recipes](#Practical-recipes):

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Configuring a Simulation](simulation_configuration.md)
- [Parallel Execution](parallel_execution.md)

## IntegrationTolerances

`IntegrationTolerances` is the `integration_tolerances` field of
`SimulationConfiguration`, and SpaceAGORA exports it. Its fields and their
defaults:

```julia
IntegrationTolerances(
    reltol               = 1e-9,
    abstol               = 1e-11,
    reltol_orbit         = 1e-6,
    abstol_orbit         = 1e-8,
    reltol_atmosphere    = 1e-7,
    abstol_atmosphere    = 1e-9,
    reltol_quaternion    = 1e-9,
    abstol_quaternion    = 1e-11,
    reltol_mass          = 1e-8,
    abstol_mass          = 1e-10,
    reltol_heat_load     = 1e-7,
    abstol_heat_load     = 1e-9,
    reltol_angular_rate  = 1e-8,
    abstol_angular_rate  = 1e-10,
    dt_max               = 1.0,
    dt_max_orbit         = 30.0,
    dt_max_atmosphere    = 1.0
)
```

`reltol`, `abstol` and `dt_max` stay in the struct but the run does not read
them. The sections below describe the others.

## What each tolerance controls

The orbit pair is the baseline: it applies to every state component that has
no pair of its own. Inside the atmosphere phase the atmosphere pair takes its
place. A run has an atmosphere phase only when it has a density model other
than `NoAtmosphereModel()` or a built-in aerodynamic model
(`AerodynamicCoefficientfM`, `AerodynamicCoefficientConstant` or
`AerodynamicCoefficientNoBallisticFlight`), and the phase lasts while any
active spacecraft is inside the entry interface, a sphere `EI` km above the
planet's equatorial radius.

| Fields | Apply to | Defaults |
| --- | --- | --- |
| `reltol_orbit`, `abstol_orbit` | every state component without a pair of its own, outside the atmosphere phase or in a run without one | `1e-6`, `1e-8` |
| `reltol_atmosphere`, `abstol_atmosphere` | the same components, inside the atmosphere phase | `1e-7`, `1e-9` |
| `reltol_quaternion`, `abstol_quaternion` | the attitude quaternion when `orientation_sim = true`, and the joint coordinates and attachment attitudes of articulated and compliant spacecraft | `1e-9`, `1e-11` |
| `reltol_angular_rate`, `abstol_angular_rate` | the body angular rate when `orientation_sim = true`, and joint and attachment rates | `1e-8`, `1e-10` |
| `reltol_mass`, `abstol_mass` | each spacecraft's mass state | `1e-8`, `1e-10` |
| `reltol_heat_load`, `abstol_heat_load` | each link's accumulated heat load | `1e-7`, `1e-9` |

A mass, heat-load or angular-rate tolerance set to `0.0` follows the pair of
the current phase, orbit or atmosphere; a nonzero value keeps its value in
both phases, and a negative one is refused. The quaternion pair has no such
fallback.

The orbit defaults suit most orbit propagation and aerobraking runs. Tighten
them to `1e-8` / `1e-10` for high-accuracy orbit work, or loosen them to
`1e-5` / `1e-7` for long Monte Carlo campaigns where absolute accuracy matters
less. The atmosphere defaults are tighter because aerodynamic forces change
quickly during a drag pass. Loosening them is the first thing to try when long
entry simulations are too slow, and `1e-6` / `1e-8` is often enough for
preliminary aerobraking studies.

With `orientation_sim = true`, the run also projects the attitude quaternion
back to unit norm after every accepted step, so the quaternion tolerances
govern attitude accuracy rather than the norm. If attitude drifts in long runs,
tighten them to `1e-10` / `1e-12`. The mass defaults suit typical thruster
burns, and the heat-load defaults suit entry and aerobraking; on orbit-only
runs the heat-load pair matters little.

## Step size limits

`dt_max_orbit` (default 30 s) is the largest step the integrator takes outside
the atmosphere phase, and in every run without one. Raising it lets long
low-drag arcs take larger steps. It does not set the output cadence: rows are
saved every `data_rate` seconds of `MissionConfiguration` either way. In the
`:gravity_backbone_split` mode it is also the fixed step, unless
`SolverConfig.gravity_backbone_dt_s` sets another.

`dt_max_atmosphere` (default 1 s) is the largest step inside the atmosphere
phase. It is the step limit that matters most for entry and aerobraking
accuracy; for steep entries, or when aerodynamic forces change quickly, reduce
it to `0.5` or `0.1` s.

`dt_max` (default 1 s) is not read by the run.

In adaptive runs, one integrator advances every spacecraft, so it uses the
atmosphere step limit and tolerances whenever any active spacecraft is inside
the entry interface. It picks the phase settings at the start, at a checkpoint
restart and when it reuses an integrator, and updates them when a spacecraft
crosses the entry interface or is deactivated.

In the `:multirate` route, each subsolve starts from the step limit and
tolerances that the route supplies; a crossing callback can still apply the
phase settings within that subsolve.

The fixed-step gravity-backbone driver (`:gravity_backbone_split`) takes no
phase updates. Gravity-backbone runs sample no atmosphere: `run_simulation`
refuses one whose density model is not `NoAtmosphereModel()`, however the
solver mode was chosen.

## Practical recipes

For faster long orbit runs, at reduced accuracy:

```julia
IntegrationTolerances(
    reltol_orbit      = 1e-5,
    abstol_orbit      = 1e-7,
    dt_max_orbit      = 60.0
)
```

For a higher-accuracy aerobraking pass:

```julia
IntegrationTolerances(
    reltol_orbit       = 1e-8,
    abstol_orbit       = 1e-10,
    reltol_atmosphere  = 1e-9,
    abstol_atmosphere  = 1e-11,
    dt_max_atmosphere  = 0.5
)
```

The tolerances `make_example_config` sets:

```julia
IntegrationTolerances(
    reltol_orbit       = 1e-8,
    abstol_orbit       = 1e-8,
    dt_max_orbit       = 20.0,
    reltol_atmosphere  = 1e-8,
    abstol_atmosphere  = 1e-8,
    dt_max_atmosphere  = 0.2
)
```

Fields not listed keep their struct defaults.

## Solver policy and partitioning

The integration mode is a setting; SpaceAGORA does not choose it for you. The
default is `solver_mode = :tsit5`, the explicit route. A solver configuration
given on the simulation configuration takes precedence. In outline, with `...`
standing for the other arguments:

```julia
args = SimulationConfiguration(...; solver_config=SolverConfig(solver_mode=:split_imex))
```

When `solver_config` is left as `nothing`, `run_simulation` reads the engine
settings in effect when it is called, so `SPACEAGORA_SOLVER_MODE` (or a
`SimulationEngineConfig` override) applies. The mode chosen either way decides
both the integrator and the callbacks the run installs. The accepted modes are
`:tsit5`, `:gravity_backbone_split`, `:dp8`, `:auto_stiff`, `:rodas5p`,
`:split_imex`, `:multirate` and `:multirate_mri`. IMEX (`:split_imex`) and
multirate integration run only when selected this way; nothing in the dynamics
model or routing profile switches a run to them.

`SolverConfig` checks its settings when it is built. An unknown
`solver_mode`, `split_imex_solver`, `multirate_slow_solver` or
`multirate_fast_solver` throws an `ArgumentError` that lists the accepted
names, and a non-positive `maxiters`, `multirate_fast_substeps` or step
(`gravity_backbone_dt_s`, `multirate_slow_dt_s`) throws one as well. The
matching environment variables (`SPACEAGORA_SOLVER_MODE`,
`SPACEAGORA_SPLIT_IMEX_SOLVER`, `SPACEAGORA_MULTIRATE_SLOW_SOLVER`,
`SPACEAGORA_MULTIRATE_FAST_SOLVER`) accept the same names and a few
alternative spellings. `:symplectic` is refused, because `run_simulation` has
no path for it; use `:gravity_backbone_split`, which integrates the gravity
core with the fixed-step KahanLi8 symplectic method.

Within `:split_imex`, each force or torque effector's `solver_partition` hook
decides whether it goes on the explicit or the implicit side of the split. The
hook does not select the mode. See [Extensibility](../extensibility.md) if you
are writing a custom effector and need to control its partition.

### Nonlinear solves in split IMEX mode

`:kencarp4` is the default algorithm for `solver_mode = :split_imex`. It
refreshes the Newton Jacobian on every nonlinear iteration
(`NLNewton(always_new=true)`), which avoids an accuracy loss that the
library's default Jacobian reuse showed on the nonlinear atmospheric
trajectories tested. The componentwise integration tolerances and the choice
of a dense or sparse linear solver still apply. `:kencarp47`, `:kencarp58` and
the other solver modes keep the library's own nonlinear-solver settings.

This does not guarantee a bound on the endpoint error: split problems whose
parts strongly cancel can still lose accuracy. Validate the chosen method and
tolerances on the trajectory you intend to run. Refreshing the Jacobian and
its factorization on every iteration adds work, and the cost for large
constellations has not been measured.

## Coupled multirate integration

For rigid spacecraft connected by interlinks, select `:multirate_mri`. The
solver settings below are complete; the configuration is in outline:

```julia
solver = SolverConfig(
    solver_mode=:multirate_mri,
    multirate_slow_dt_s=0.5,
    multirate_fast_substeps=8,
)
args = SimulationConfiguration(...; solver_config=solver)
```

This opt-in route uses OrdinaryDiffEq's explicit MRI-GARK-ERK45a method. The
slow component holds orbit and environmental forces and torques, quaternion
kinematics, gyroscopic terms and custom dynamic effectors. Control and
interlink forces and torques, power, data production and transfer, data stores
and ground-station counters form the fast component. `solver_partition` does
not move a custom effector into the fast component. Both components act on the
same spacecraft states. One integrator owns the callbacks, so a guidance tick,
sensor sample, store boundary or other event is processed on one physical
timeline. Fast attitude motion can still limit the outer step, because its
kinematics are in the slow component.

`multirate_slow_dt_s` caps the outer step. If it is omitted, the cap is
`min(dt_max_orbit, 2.0)` seconds. The step limit of the current phase,
`dt_max_orbit` or `dt_max_atmosphere`, can shorten it further.
`multirate_fast_substeps` (default 8) is the number of RK4 microsteps per MRI
stage, not per outer step.

!!! warning "Fast microstep accuracy needs a separate convergence check"
    The outer method adapts its step using an embedded estimate of the slow
    coupling error. That estimate does not control the error or stability of
    the fast RK4 microsteps. Tightening the tolerances alone is not enough: an
    under-resolved fast controller can return `Success` with badly wrong
    intermediate states.

    Choose `multirate_fast_substeps` and `multirate_slow_dt_s` to resolve the
    fastest time constant in the fast component. Increase the substep count
    and reduce the outer cap separately, and compare the full trajectory and
    the event and output quantities with a suitable single-rate reference. A
    plausible final state alone does not establish convergence. For stiff fast
    dynamics, use a full-state stiff solver such as `:rodas5p` or
    `:auto_stiff`.

For the scalar fast decay `y′ = -λy` with `λ > 0`, each RK4 microstep must
satisfy approximately `λ * δt ≤ 2.785` for linear stability. In
OrdinaryDiffEqMultirate 3.0.3, the version in use, the method's five stages
have equal duration, so `δt = outer_step / (5 * multirate_fast_substeps)`.
This scalar bound is neither an accuracy target nor a stability test for a
nonlinear coupled model; use a smaller step and verify convergence. The
package compat entry holds the dependency to the 3.0 patch series
(`OrdinaryDiffEqMultirate = "~3.0.3"`).

The `multirate_slow_solver` and `multirate_fast_solver` fields configure only
the `:multirate` route, a Strang splitting. They do not change MRI's internal
methods or enable automatic switching to a stiff solver. MRI does not reuse an
integrator through `solver_cache`.

The MRI route supports rigid-spacecraft interlinks, batteries, data stores,
mounted sensors and normal checkpoint continuation. When you check a resumed
run against an uninterrupted one, keep the solver mode and numerical settings
the same: switching between MRI and a single-rate mode is accepted without a
warning and does not promise the same trajectory. A same-mode MRI restart
recovers identical solver states, but saved rows at checkpoint boundaries can
differ by floating-point rounding.

MRI uses Hermite dense output. A custom discrete callback that jumps the state
can affect interpolated output rows earlier in that accepted step, and putting
the callback time in `tstops` does not remove this. Validate saved rows around
such jumps rather than inferring their accuracy from the accepted solver
states.

Articulated spacecraft, compliant attachments, externally propagated
spacecraft and operational mission segments accept only `:tsit5`,
`:auto_stiff`, `:rodas5p` and `:dp8`, so they refuse this mode, and coupled
robot-arm effectors refuse it too. MRI is a numerical option, not a promise
that coupled runs will be faster; measure accuracy and run time on the case
you intend to run.

## Optional KS propagation candidate

A separate optional KSPropagator extension provides
`SolverConfig(solver_mode=:ks)` for the bounded one-spacecraft development
candidate. SpaceAGORA supplies configured physical forces and retains physical
time, scheduling, events, checkpoints and outputs. KSPropagator supplies the
coordinate transformation and coupled numerical stepping. See
[`verification/optional_ks/README.md`](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/feature/optional-ks-propagation/verification/optional_ks/README.md)
for the exact KS installation source, supported models, guards and reproducible
qualification commands. Consult the candidate's source-specific report and capability table for the
verified initial scope, failed settings and installation limitations. Performance
and final delivery review remain pending.
