# Solver selection and numerical tolerances, shared with engine runtime settings.
# Included inside SimConfig by configuration.jl.

# 2.1. Solver Configuration
"""
    SolverConfig

Typed runtime configuration for solver selection, fixed-step settings, and
multirate/IMEX integration policy. Fields correspond to the `SPACEAGORA_SOLVER_*`
environment variables.

When `solver_config` is `nothing` on a `SimulationConfiguration`, `run_simulation`
reads the effective config from the active environment at call time (respecting any
`SimulationEngineConfig` overrides). Set this field explicitly to pin solver behavior
independent of environment variables.

`split_imex` uses the atmosphere-implicit IMEX partition. `multirate` keeps the
control-focused split path. `gravity_backbone_split` is a fixed-step symplectic
gravity-backbone mode; it is not a fully symplectic whole-system solve.

# Parallel execution

`parallel = true` is the one switch for parallel execution. It lets SpaceAGORA
choose how to use the threads Julia was started with, for this run only:

- a single `run_simulation` call threads its callbacks, effectors and
  right-hand side where the adaptive inner policy predicts a gain;
- `run_constellation_ensemble` (whose members carry this configuration) and
  `run_monte_carlo(...; parallel=true)` let the predictive campaign planner
  choose between serial, threaded and process-worker execution before the
  first sample runs.

Start Julia with threads (`julia --threads=auto`) for the flag to have
anything to use. On the first parallel run on a machine SpaceAGORA measures
the machine's cost constants once (a few seconds) and stores them under
`output/parallel_policy_state/`; later runs reuse them. Settings are scoped
to the call: nothing is left behind in the process environment afterwards.
With `parallel = false` (the default) behavior is unchanged: runs and
campaigns are serial unless a thread count is given explicitly.

```julia
args = SimulationConfiguration(...; solver_config=SolverConfig(parallel=true))
run_simulation(args)
```
"""
Base.@kwdef struct SolverConfig
    solver_mode::Symbol = :tsit5
    maxiters::Union{Nothing, Int} = nothing
    symplectic_dt_s::Union{Nothing, Float64} = nothing
    gravity_backbone_dt_s::Union{Nothing, Float64} = nothing
    split_imex_solver::Symbol = :kencarp4
    multirate_slow_dt_s::Union{Nothing, Float64} = nothing
    multirate_fast_substeps::Int = 8
    multirate_slow_solver::Symbol = :tsit5
    multirate_fast_solver::Symbol = :auto_stiff
    auto_stiff_gravity_tsit5::Bool = true
    auto_stiff_switch_max::Int = 50
    parallel::Bool = false
end

# 2.3. Integration Tolerances
@kwdef struct IntegrationTolerances
    reltol::Float64 = 1e-9
    abstol::Float64 = 1e-11
    reltol_orbit::Float64 = 1e-6
    abstol_orbit::Float64 = 1e-8
    reltol_atmosphere::Float64 = 1e-7
    abstol_atmosphere::Float64 = 1e-9
    reltol_quaternion::Float64 = 1e-9
    abstol_quaternion::Float64 = 1e-11
    reltol_mass::Float64 = 1e-8
    abstol_mass::Float64 = 1e-10
    reltol_heat_load::Float64 = 1e-7
    abstol_heat_load::Float64 = 1e-9
    reltol_angular_rate::Float64 = 1e-8
    abstol_angular_rate::Float64 = 1e-10
    dt_max::Float64 = 1.0
    dt_max_orbit::Float64 = 30.0
    dt_max_atmosphere::Float64 = 1.0
end # struct IntegrationTolerances
