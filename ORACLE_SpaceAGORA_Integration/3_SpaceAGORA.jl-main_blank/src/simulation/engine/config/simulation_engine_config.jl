"""
    SimulationEngineConfig

Top-level typed runtime configuration consumed by `SimulationEngine.run_simulation`.
It composes parallel, solver, runtime policy, and artifact settings into one
configuration object.
"""
Base.@kwdef struct SimulationEngineConfig
    parallel::ParallelConfig = ParallelConfig()
    solver::SolverConfig = SolverConfig()
    runtime_policy::RuntimePolicyConfig = RuntimePolicyConfig()
    artifacts::ArtifactConfig = ArtifactConfig()
    env_overrides::Dict{String, String} = Dict{String, String}()
end
