Base.@kwdef struct SimulationEngineConfig
    parallel::ParallelConfig = ParallelConfig()
    solver::SolverConfig = SolverConfig()
    runtime_policy::RuntimePolicyConfig = RuntimePolicyConfig()
    artifacts::ArtifactConfig = ArtifactConfig()
    env_overrides::Dict{String, String} = Dict{String, String}()
end
