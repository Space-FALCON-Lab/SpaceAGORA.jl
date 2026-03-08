## Canonical aggregator: no behavior ownership.
module SimulationEngine

using ..SimulationModel
import ..RuntimeServices
using SPICE

include(joinpath(@__DIR__, "config", "parallel_config.jl"))
include(joinpath(@__DIR__, "config", "solver_config.jl"))
include(joinpath(@__DIR__, "config", "runtime_policy_config.jl"))
include(joinpath(@__DIR__, "config", "artifact_config.jl"))
include(joinpath(@__DIR__, "config", "simulation_engine_config.jl"))

include(joinpath(@__DIR__, "adapters", "from_env.jl"))
include(joinpath(@__DIR__, "adapters", "from_simulation_configuration.jl"))

include(joinpath(@__DIR__, "effector_sampling.jl"))
include(joinpath(@__DIR__, "setup.jl"))
include(joinpath(@__DIR__, "solver_policy.jl"))
include(joinpath(@__DIR__, "dynamics_rhs.jl"))
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
export prewarm_nbody_ephemeris_cache
export load_nbody_ephemeris_cache!

end # module SimulationEngine
