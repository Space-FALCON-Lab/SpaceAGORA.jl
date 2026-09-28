const _RUN_SIMULATION_TYPED_BOUNDARY_DEPRECATION =
    "run_simulation now requires a typed SimulationConfiguration at the execution boundary. " *
    "Passing other inputs is deprecated; construct a SimulationConfiguration first."

@inline function _require_simulation_configuration(args)
    args isa SimulationConfiguration || throw(ArgumentError("run_simulation expects SimulationConfiguration; got $(typeof(args))."))
    return args
end

@noinline function _depwarn_untyped_run_simulation(args)
    Base.depwarn(
        _RUN_SIMULATION_TYPED_BOUNDARY_DEPRECATION * " Got $(typeof(args)).",
        :run_simulation;
        force=true
    )
    return nothing
end

function run_simulation(config::SimulationEngineConfig, args::SimulationConfiguration; kwargs...)
    return _with_engine_env_overrides(config, () -> run_simulation(args; kwargs...))
end

function run_simulation(config::SimulationEngineConfig, args; kwargs...)
    _depwarn_untyped_run_simulation(args)
    typed_args = _require_simulation_configuration(args)
    return _with_engine_env_overrides(config, () -> run_simulation(typed_args; kwargs...))
end

function run_simulation(args; kwargs...)
    _depwarn_untyped_run_simulation(args)
    typed_args = _require_simulation_configuration(args)
    return run_simulation(typed_args; kwargs...)
end

function prewarm_nbody_ephemeris_cache(config::SimulationEngineConfig, args; kwargs...)
    return _with_engine_env_overrides(config, () -> prewarm_nbody_ephemeris_cache(args; kwargs...))
end

"""
    prewarm_nbody_ephemeris_cache(args; dt_s=nothing, mission_end_s=nothing, save_path=nothing) -> cache
    prewarm_nbody_ephemeris_cache(config, args; dt_s=nothing, mission_end_s=nothing, save_path=nothing) -> cache

Precompute and register a process-local N-body SPICE ephemeris cache for later
[`run_simulation`](@ref) calls. This is intended for Monte Carlo campaigns that
reuse the same third-body set, start epoch, mission span, and cache sample
spacing across many runs. The returned cache is keyed by the same deterministic
boundary that the runtime setup already uses.

If `save_path` is provided, the cache is also serialized to disk so other Julia
worker processes can call [`load_nbody_ephemeris_cache!`](@ref) and reuse the
same precomputed ephemeris without rebuilding it from SPICE.
"""
function prewarm_nbody_ephemeris_cache(
    args;
    dt_s::Union{Nothing, Real}=nothing,
    mission_end_s::Union{Nothing, Real}=nothing,
    save_path::Union{Nothing, AbstractString}=nothing
)
    return _prewarm_nbody_ephemeris_cache(
        args;
        dt_s=dt_s,
        mission_end_s=mission_end_s,
        save_path=save_path
    )
end

"""
    load_nbody_ephemeris_cache!(path; replace=true) -> cache

Load a serialized N-body ephemeris cache created by
[`prewarm_nbody_ephemeris_cache`](@ref) and register it in the current Julia
process so later [`run_simulation`](@ref) calls can reuse it. This is intended
for multi-process Monte Carlo campaigns where each worker should load the same
precomputed SPICE cache once before running many trajectories.
"""
function load_nbody_ephemeris_cache!(path::AbstractString; replace::Bool=true)
    return _load_nbody_ephemeris_cache!(String(path); replace=replace)
end
