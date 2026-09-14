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

function load_nbody_ephemeris_cache!(path::AbstractString; replace::Bool=true)
    return _load_nbody_ephemeris_cache!(String(path); replace=replace)
end
