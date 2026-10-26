# Solver Configuration

Use this page when you need to tune integration tolerances, step size limits,
or understand which tolerance parameter affects which part of the simulation.

This page is for users who are hitting accuracy or performance issues and need
to adjust the solver beyond the defaults.

Shortest successful command:

```text
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

What to read next:

- [Simulation Configuration](simulation_configuration.md)
- [Parallel Execution](parallel_execution.md)

## IntegrationTolerances

`IntegrationTolerances` is passed as the `integration_tolerances` field of
`SimulationConfiguration`. The defaults are set for typical orbit and
aerobraking accuracy:

```julia
SM.IntegrationTolerances(
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

## What each tolerance controls

**`reltol` / `abstol`**
Global fallback tolerances used for state components not covered by a more
specific pair. Tightening these is a blunt instrument; prefer the per-subsystem
pairs below.

**`reltol_orbit` / `abstol_orbit`**
Applied to the translational position and velocity states during the
free-flight (orbital) phase. The defaults (`1e-6` / `1e-8`) are appropriate
for most orbit propagation and aerobraking scenarios. Tighten to `1e-8` /
`1e-10` for high-accuracy orbit determination studies; loosen to `1e-5` /
`1e-7` to speed up long Monte Carlo runs where absolute accuracy is less
critical.

**`reltol_atmosphere` / `abstol_atmosphere`**
Applied during the atmosphere-coupled (drag-pass) phase. The tighter defaults
here (`1e-7` / `1e-9` compared to orbit) reflect that aerodynamic forces
change rapidly during a drag pass. Loosening these is the first lever to reach
for when long entry simulations are too slow. Values of `1e-6` / `1e-8` are
often acceptable for preliminary aerobraking studies.

**`reltol_quaternion` / `abstol_quaternion`**
Applied to the attitude quaternion when `orientation_sim = true`. The defaults
keep the quaternion norm near unity throughout. If attitude drift is visible in
long runs, tighten to `1e-10` / `1e-12`.

**`reltol_mass` / `abstol_mass`**
Applied to the propellant mass state. The defaults are appropriate for typical
thruster burn scenarios.

**`reltol_heat_load` / `abstol_heat_load`**
Applied to the accumulated heat load integral. The defaults are appropriate
for entry and aerobraking scenarios; these have less impact on orbit-only runs.

**`reltol_angular_rate` / `abstol_angular_rate`**
Applied to the body angular rate when `orientation_sim = true`.

## Step size limits

**`dt_max`**
Global maximum step size in simulated seconds. This is the fallback if no
phase-specific limit applies. Rarely needs adjustment.

**`dt_max_orbit`**
Maximum step during the free-flight orbital phase (default 30 s). Increasing
this speeds up long low-drag phases at the cost of temporal resolution in the
output. The output `data_rate` in `MissionConfiguration` is independent —
`dt_max_orbit` controls the integrator, not the output cadence.

**`dt_max_atmosphere`**
Maximum step during the atmosphere-coupled phase (default 1 s). This is the
most important step limit for entry and aerobraking accuracy. For steep entries
or when aerodynamic forces are changing rapidly, reduce to `0.1` or `0.5` s.

For adaptive atmospheric propagation, the integrator uses the atmospheric
cap and tolerances whenever any active spacecraft is inside the configured
atmospheric boundary. It selects these settings at startup, checkpoint restart
and integrator reuse, then updates them on crossings and spacecraft deactivation.
Explicit per-component mass, heat-load and attitude tolerance overrides are preserved.

Explicit multirate subsolve overrides restore their supplied cap and tolerances
at every subsolve. Installed crossing callbacks can apply phase settings within
that subsolve. Previously, a crossing inside a reused subsolve carried its phase
settings into later subsolves, so physical multirate runs with atmospheric
crossings change.

The fixed-step gravity-backbone driver is exempt from adaptive
phase updates. Gravity-backbone runs sample no atmosphere: `run_simulation`
refuses one whose density model is not `NoAtmosphereModel()`, however the
solver mode was chosen.

## Practical recipes

**Faster long orbit runs (reduced accuracy):**

```julia
SM.IntegrationTolerances(
    reltol_orbit      = 1e-5,
    abstol_orbit      = 1e-7,
    dt_max_orbit      = 60.0
)
```

**Higher-accuracy aerobraking pass:**

```julia
SM.IntegrationTolerances(
    reltol_orbit       = 1e-8,
    abstol_orbit       = 1e-10,
    reltol_atmosphere  = 1e-9,
    abstol_atmosphere  = 1e-11,
    dt_max_atmosphere  = 0.5
)
```

**Default used by `make_example_config`:**

```julia
SM.IntegrationTolerances(
    reltol_orbit       = 1e-8,
    abstol_orbit       = 1e-8,
    dt_max_orbit       = 20.0,
    reltol_atmosphere  = 1e-8,
    abstol_atmosphere  = 1e-8,
    dt_max_atmosphere  = 0.2
)
```

All other fields take their struct defaults when not specified.

## Solver policy and partitioning

The integration mode is a setting, not an automatic choice. The default is
`solver_mode = :tsit5`, the explicit route. An explicit solver configuration on
the simulation configuration wins:

```julia
args = SimulationConfiguration(...; solver_config=SolverConfig(solver_mode=:split_imex))
```

When `solver_config` is left as `nothing`, `run_simulation` reads the effective
engine settings at call time, so `SPACEAGORA_SOLVER_MODE` (or a
`SimulationEngineConfig` override) applies. The mode chosen this way decides
both the integrator and the callbacks the run installs. Accepted modes are `:tsit5`,
`:gravity_backbone_split`, `:dp8`, `:auto_stiff`, `:rodas5p`, `:split_imex`
`:multirate` and `:multirate_mri`. IMEX (`:split_imex`) and multirate integration run only when
selected this way; nothing in the dynamics model or routing profile switches a
run to them.

`SolverConfig` checks its settings when it is built: an unknown `solver_mode`,
`split_imex_solver` or `multirate_*_solver` name, or a non-positive
`maxiters`, `multirate_fast_substeps` or step size, throws an `ArgumentError`
that names the accepted values. The matching environment variables
(`SPACEAGORA_SOLVER_MODE`, `SPACEAGORA_SPLIT_IMEX_SOLVER`,
`SPACEAGORA_MULTIRATE_SLOW_SOLVER`, `SPACEAGORA_MULTIRATE_FAST_SOLVER`) accept
the same names plus a few alternative spellings. `:symplectic` is refused:
`run_simulation` has no path for it. Use `:gravity_backbone_split`, which
integrates the gravity core with the fixed-step KahanLi8 symplectic method.

Within `:split_imex`, the solver partition hooks on force/torque effectors
(`solver_partition`) determine whether an effector goes on the explicit or
implicit side of the split. They do not select the mode. See
[Extensibility](../extensibility.md) if you are building a custom effector and
need to control partitioning explicitly.

### Nonlinear solves in split IMEX mode

`solver_mode = :split_imex` retains `:kencarp4` as its default algorithm.
KenCarp4 explicitly refreshes the Newton Jacobian on each nonlinear iteration using
`NLNewton(always_new=true)`. This avoids the accuracy loss observed with
library Jacobian-reuse defaults in the tested nonlinear atmospheric trajectories.
Componentwise integration tolerances and the dense or
sparse linear-solver choice are preserved.

This policy changes split-IMEX KenCarp4 numerical results relative to the earlier
library-default policy. It does not establish a general endpoint-error bound:
strongly cancelling split problems retain a separately documented accuracy
limitation. Validate the chosen method and tolerances for the intended
trajectory. KenCarp47, KenCarp58 and other solver modes retain their existing
nonlinear-solver policy. Refreshing Jacobians and their factorizations adds work;
the runtime cost for large constellations has not been measured.


## Coupled multirate integration

For rigid spacecraft connected by interlinks, select `:multirate_mri`:

```julia
solver = SolverConfig(
    solver_mode=:multirate_mri,
    multirate_slow_dt_s=0.5,
    multirate_fast_substeps=8,
)
args = SimulationConfiguration(...; solver_config=solver)
```

This opt-in route uses OrdinaryDiffEq's explicit MRI-GARK45a method. The slow
component includes orbit and environmental forces and torques, quaternion
kinematics, gyroscopic terms and custom dynamic effectors. Control and interlink
forces and torques, power, data production and transfer, stores and station
counters form the fast component. `solver_partition` does not move a custom
effector into MRI's fast component. Both components act on the same spacecraft
states. One integrator owns the callbacks, so a guidance tick, sensor sample,
store boundary or other event is processed on one physical timeline. Fast
attitude motion can still limit the outer step because its kinematics are slow.

`multirate_slow_dt_s` limits the outer step. If omitted, its cap is
`min(dt_max_orbit, 2.0)` seconds; orbital and atmospheric step limits can shorten
it further. `multirate_fast_substeps` is the number of RK4 microsteps **per MRI
stage**, not per outer step.

!!! warning "Fast microstep accuracy needs a separate convergence check"
    The outer method adapts its step using an embedded estimate of the slow
    coupling error. That estimate does **not** control the error or stability
    of the fast RK4 microsteps. Tightening `reltol` and `abstol` alone is not
    sufficient: an under-resolved fast controller can return `Success` with
    severely wrong intermediate states.

    Choose `multirate_fast_substeps` and `multirate_slow_dt_s` to resolve the
    fastest time constant in the fast component. Increase the substep count
    and reduce the outer cap separately, and compare the full trajectory and
    event/output quantities against a suitable single-rate reference. A
    plausible final state alone does not establish convergence. Use a
    full-state stiff solver such as `:rodas5p` or `:auto_stiff` for stiff fast
    dynamics.

For the scalar fast decay `y′ = -λy` with `λ > 0`, each RK4 microstep must
satisfy approximately `λ * δt ≤ 2.785` for linear stability. In the selected
OrdinaryDiffEqMultirate 3.0.3 method, the five stages have equal duration, so
`δt = outer_step / (5 * multirate_fast_substeps)`. This scalar stability bound
is not an accuracy target or a general stability test for a nonlinear coupled
model; use a smaller step and verify convergence. The dependency is constrained
to the 3.0 patch series while this implementation is qualified.

The `multirate_slow_solver` and `multirate_fast_solver` fields configure only the
original `:multirate` Strang route. They do not change MRI's internal methods or
enable automatic switching to a stiff solver. MRI does not reuse an integrator
through `solver_cache`.

The MRI route supports rigid-spacecraft interlinks, batteries, data stores,
mounted sensors and normal checkpoint continuation. Keep the solver mode and
numerical settings unchanged when checking resumed runs against uninterrupted
runs: switching between MRI and the ordinary single-rate modes is currently
accepted without a warning and does not promise the same trajectory. Same-mode
MRI restart checks recover identical solver states, but saved rows at checkpoint
boundaries can differ by floating-point rounding.

MRI uses Hermite dense output. A custom discrete callback that jumps the state
can affect interpolated output rows earlier in that accepted step; merely
putting the callback time in `tstops` does not remove this limitation. Validate
saved rows around such jumps rather than inferring their accuracy from the
accepted solver states.

Articulated spacecraft, compliant attachments, coupled robot-arm effectors,
external propagation and mission segments retain their setup refusals. The
existing `:multirate` route and default `:tsit5` behavior are unchanged. This is
a numerical option, not a promise that coupled runs will be faster; measure
accuracy and runtime on the intended case.
