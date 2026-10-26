# Solver selection and numerical tolerances, shared with engine runtime settings.
# Included inside SimConfig by configuration.jl.

# Solver names SolverConfig accepts. The SPACEAGORA_SOLVER_* parsers
# (engine/adapters/from_env.jl) check against these same tuples, so the typed
# and environment paths cannot accept different names.
const _SOLVER_MODES = (:tsit5, :gravity_backbone_split, :dp8, :auto_stiff, :rodas5p, :split_imex, :multirate, :multirate_mri)
const _SPLIT_IMEX_SOLVERS = (:kencarp4, :kencarp47, :kencarp58)
const _MULTIRATE_SOLVERS = (:tsit5, :dp8, :auto_stiff, :rodas5p, :kencarp4)

# run_simulation never builds the partitioned second-order problem a standalone
# symplectic solve needs, so :symplectic is refused where it is configured
# rather than failing inside the run.
const _SYMPLECTIC_MODE_MESSAGE =
    "run_simulation has no partitioned (second-order) path for the :symplectic mode. " *
    "Use :gravity_backbone_split, which integrates the gravity core with the fixed-step KahanLi8 " *
    "symplectic method; gravity_backbone_dt_s (SPACEAGORA_GRAVITY_BACKBONE_DT_S) sets its step."

_solver_name_list(names) = join((repr(name) for name in names), ", ")

function _check_solver_name(field::String, value::Symbol, names)
    value in names && return nothing
    throw(ArgumentError("Unsupported SolverConfig.$field=$(repr(value)). Use one of: $(_solver_name_list(names))."))
end

function _check_positive(field::String, value)
    (value === nothing || value > 0) && return nothing
    throw(ArgumentError("SolverConfig.$field must be positive, got $(value)."))
end

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
control-focused Strang split path. `multirate_mri` uses explicit MRI-GARK45a
with a single callback timeline for rigid spacecraft, including interlinks,
batteries and data stores. Its microsteps use RK4; `multirate_slow_solver`
and `multirate_fast_solver` apply only to the original `multirate` mode.
`gravity_backbone_split` is a fixed-step symplectic
gravity-backbone mode; it is not a fully symplectic whole-system solve.

The constructor checks the settings and throws an `ArgumentError` naming the
accepted values for an unknown name:

- `solver_mode`: `:tsit5`, `:gravity_backbone_split`, `:dp8`, `:auto_stiff`,
  `:rodas5p`, `:split_imex`, `:multirate`, `:multirate_mri`. `:symplectic` is refused; use
  `:gravity_backbone_split`.
- `split_imex_solver`: `:kencarp4`, `:kencarp47`, `:kencarp58`.
- `multirate_slow_solver`, `multirate_fast_solver`: `:tsit5`, `:dp8`,
  `:auto_stiff`, `:rodas5p`, `:kencarp4`.
- `maxiters`, `multirate_fast_substeps` and the `*_dt_s` steps must be positive
  (`nothing` keeps a default).

`symplectic_dt_s` belonged to the `:symplectic` mode and has no effect; it is
kept so existing constructor calls still work.

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

    function SolverConfig(solver_mode, maxiters, symplectic_dt_s, gravity_backbone_dt_s,
                          split_imex_solver, multirate_slow_dt_s, multirate_fast_substeps,
                          multirate_slow_solver, multirate_fast_solver, auto_stiff_gravity_tsit5,
                          auto_stiff_switch_max, parallel)
        cfg = new(solver_mode, maxiters, symplectic_dt_s, gravity_backbone_dt_s,
                  split_imex_solver, multirate_slow_dt_s, multirate_fast_substeps,
                  multirate_slow_solver, multirate_fast_solver, auto_stiff_gravity_tsit5,
                  auto_stiff_switch_max, parallel)
        cfg.solver_mode === :symplectic && throw(ArgumentError(
            "SolverConfig.solver_mode=:symplectic is not supported: " * _SYMPLECTIC_MODE_MESSAGE))
        _check_solver_name("solver_mode", cfg.solver_mode, _SOLVER_MODES)
        _check_solver_name("split_imex_solver", cfg.split_imex_solver, _SPLIT_IMEX_SOLVERS)
        _check_solver_name("multirate_slow_solver", cfg.multirate_slow_solver, _MULTIRATE_SOLVERS)
        _check_solver_name("multirate_fast_solver", cfg.multirate_fast_solver, _MULTIRATE_SOLVERS)
        _check_positive("maxiters", cfg.maxiters)
        _check_positive("multirate_fast_substeps", cfg.multirate_fast_substeps)
        _check_positive("symplectic_dt_s", cfg.symplectic_dt_s)
        _check_positive("gravity_backbone_dt_s", cfg.gravity_backbone_dt_s)
        _check_positive("multirate_slow_dt_s", cfg.multirate_slow_dt_s)
        return cfg
    end
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
