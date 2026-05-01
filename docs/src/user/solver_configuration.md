# Solver Configuration

Use this page when you need to tune integration tolerances, step size limits,
or understand which tolerance parameter affects which part of the simulation.

This page is for users who are hitting accuracy or performance issues and need
to adjust the solver beyond the defaults.

Shortest successful command:

```text
julia --project=. examples/AGORA_Basic_Quickstart.jl
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

The solver automatically selects between explicit, IMEX (implicit-explicit),
and multirate strategies based on the dynamics model and profile. The solver
partition hooks on force/torque effectors (`solver_partition`) determine
whether an effector goes on the explicit or implicit side of `split_imex`. For
most users this is transparent; see [Extensibility](../extensibility.md) if you
are building a custom effector and need to control partitioning explicitly.
