## Canonical aggregator: no behavior ownership.
module SimulationEngine

using ..SimulationModel
import ..RuntimeServices
import DiffEqBase
import ADTypes: AutoFiniteDiff
using SPICE
using SparseArrays
# KLUFactorization: sparse linsolve for the auto_stiff Rodas5P branch
# (solver_policy.jl's _rodas5p_alg) when a block-diagonal jac_prototype is
# available (multi-satellite path). LinearSolve is already resolved
# transitively via OrdinaryDiffEq, but was never imported by name here, so
# any multi-satellite case with an atmosphere-requiring effector (not
# eligible for the plain-Tsit5 smooth-gravity bypass) crashed with
# `UndefVarError: KLUFactorization not defined` the moment OrdinaryDiffEq's
# stiffness autoswitch actually engaged Rodas5P -- independent of solver
# correctness, purely a missing import.
import LinearSolve: KLUFactorization

include(joinpath(@__DIR__, "config", "parallel_config.jl"))
include(joinpath(@__DIR__, "config", "solver_config.jl"))
include(joinpath(@__DIR__, "config", "runtime_policy_config.jl"))
include(joinpath(@__DIR__, "config", "artifact_config.jl"))
include(joinpath(@__DIR__, "config", "simulation_engine_config.jl"))

include(joinpath(@__DIR__, "adapters", "from_env.jl"))
include(joinpath(@__DIR__, "adapters", "from_simulation_configuration.jl"))

include(joinpath(@__DIR__, "effector_sampling.jl"))
include(joinpath(@__DIR__, "state_access.jl"))
include(joinpath(@__DIR__, "setup.jl"))
include(joinpath(@__DIR__, "solver_policy.jl"))
include(joinpath(@__DIR__, "dynamics_rhs.jl"))
include(joinpath(@__DIR__, "rhs_calibration.jl"))
include(joinpath(@__DIR__, "rhs_cost_probe.jl"))
include(joinpath(@__DIR__, "persistence.jl"))
include(joinpath(@__DIR__, "resume_checkpoint.jl"))
include(joinpath(@__DIR__, "reporting.jl"))
include(joinpath(@__DIR__, "execution.jl"))
include(joinpath(@__DIR__, "public_api.jl"))

export ParallelConfig
export SolverConfig
export RuntimePolicyConfig
export ArtifactConfig
export SimulationEngineConfig
export simulation_engine_config_from_env
export run_simulation
export SolverIntegratorCache
export prewarm_nbody_ephemeris_cache
export load_nbody_ephemeris_cache!

end # module SimulationEngine
