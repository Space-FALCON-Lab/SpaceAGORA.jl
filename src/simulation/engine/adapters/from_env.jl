@inline _env_bool(v::Bool) = v ? "1" : "0"
const _engine_active_config_ref = Ref{Union{Nothing, SimulationEngineConfig}}(nothing)
const _engine_active_overrides_ref = Ref{Union{Nothing, Dict{String, String}}}(nothing)

function _parse_bool(raw, default::Bool)
    raw === nothing && return default
    token = lowercase(strip(String(raw)))
    token in ("1", "true", "yes", "on") && return true
    token in ("0", "false", "no", "off") && return false
    return default
end

function _parse_solver_mode_sym(raw::String)::Symbol
    mode = lowercase(strip(raw))
    mode in ("tsit5", "default", "") && return :tsit5
    mode in ("symplectic", "kahanli8", "verlet") && return :symplectic
    mode in ("gravity_backbone_split", "gravity-backbone-split", "gravity_backbone", "gravity-backbone") && return :gravity_backbone_split
    mode in ("auto_stiff", "auto-stiff", "autostiff", "auto") && return :auto_stiff
    mode in ("rodas5p", "rodas", "stiff") && return :rodas5p
    mode in ("split_imex", "split-imex", "split", "imex") && return :split_imex
    mode in ("multirate", "multirate_split", "split_multirate", "mr") && return :multirate
    mode in ("dp8", "dormandprince8", "dop8") && return :dp8
    throw(ArgumentError(
        "Unsupported SPACEAGORA_SOLVER_MODE='$raw'. Use one of: tsit5, symplectic, gravity_backbone_split, dp8, auto_stiff, rodas5p, split_imex, multirate."
    ))
end

function _parse_multirate_solver_sym(raw::String, env_name::String)::Symbol
    mode = lowercase(strip(raw))
    mode in ("tsit5", "tsit", "default") && return :tsit5
    mode in ("auto_stiff", "auto-stiff", "autostiff", "auto") && return :auto_stiff
    mode in ("rodas5p", "rodas", "stiff") && return :rodas5p
    mode in ("kencarp4", "ken4") && return :kencarp4
    mode in ("dp8", "dormandprince8", "dop8") && return :dp8
    throw(ArgumentError(
        "Unsupported $(env_name)='$raw'. Use one of: tsit5, dp8, auto_stiff, rodas5p, kencarp4."
    ))
end

function _parse_split_imex_solver_sym(raw::String)::Symbol
    mode = lowercase(strip(raw))
    mode in ("kencarp4", "ken4", "default") && return :kencarp4
    mode in ("kencarp47", "ken47") && return :kencarp47
    mode in ("kencarp58", "ken58") && return :kencarp58
    throw(ArgumentError(
        "Unsupported SPACEAGORA_SPLIT_IMEX_SOLVER='$raw'. Use one of: kencarp4, kencarp47, kencarp58."
    ))
end

"""
    _solver_config_from_env([env_get]) -> SolverConfig

Build a typed `SolverConfig` from `SPACEAGORA_SOLVER_*` environment variables.
`env_get(name, default)` defaults to `_engine_env_get`, which respects any active
`SimulationEngineConfig` overrides.
"""
function _parse_positive_float_opt(raw::String, env_name::String)::Union{Nothing, Float64}
    s = strip(raw)
    isempty(s) && return nothing
    v = tryparse(Float64, s)
    v !== nothing || throw(ArgumentError("$(env_name) must be a positive number, got '$raw'."))
    v > 0.0 || throw(ArgumentError("$(env_name) must be > 0.0, got $v."))
    return v
end

function _solver_config_from_env(env_get=_engine_env_get)::SolverConfig
    solver_mode = _parse_solver_mode_sym(env_get("SPACEAGORA_SOLVER_MODE", "tsit5"))

    raw_maxiters = strip(env_get("SPACEAGORA_SOLVER_MAXITERS", ""))
    maxiters = if isempty(raw_maxiters)
        nothing
    else
        v = tryparse(Int, raw_maxiters)
        v !== nothing || throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be a positive integer, got '$raw_maxiters'."))
        v > 0 || throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be > 0, got $v."))
        v
    end

    symplectic_dt_s = _parse_positive_float_opt(env_get("SPACEAGORA_SYMPLECTIC_DT_S", ""), "SPACEAGORA_SYMPLECTIC_DT_S")
    gravity_backbone_dt_s = _parse_positive_float_opt(env_get("SPACEAGORA_GRAVITY_BACKBONE_DT_S", ""), "SPACEAGORA_GRAVITY_BACKBONE_DT_S")

    split_imex_solver = _parse_split_imex_solver_sym(env_get("SPACEAGORA_SPLIT_IMEX_SOLVER", "kencarp4"))

    multirate_slow_dt_s = _parse_positive_float_opt(env_get("SPACEAGORA_MULTIRATE_SLOW_DT_S", ""), "SPACEAGORA_MULTIRATE_SLOW_DT_S")

    multirate_fast_substeps = let v = tryparse(Int, env_get("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS", "8"))
        v !== nothing || throw(ArgumentError("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS must be a positive integer."))
        v > 0 || throw(ArgumentError("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS must be > 0, got $v."))
        v
    end

    multirate_slow_solver = _parse_multirate_solver_sym(env_get("SPACEAGORA_MULTIRATE_SLOW_SOLVER", "tsit5"), "SPACEAGORA_MULTIRATE_SLOW_SOLVER")
    multirate_fast_solver = _parse_multirate_solver_sym(env_get("SPACEAGORA_MULTIRATE_FAST_SOLVER", "auto_stiff"), "SPACEAGORA_MULTIRATE_FAST_SOLVER")

    auto_stiff_gravity_tsit5 = SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5", true)
    auto_stiff_switch_max = SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_AUTO_STIFF_SWITCH_MAX", 50)

    return SolverConfig(
        solver_mode=solver_mode,
        maxiters=maxiters,
        symplectic_dt_s=symplectic_dt_s,
        gravity_backbone_dt_s=gravity_backbone_dt_s,
        split_imex_solver=split_imex_solver,
        multirate_slow_dt_s=multirate_slow_dt_s,
        multirate_fast_substeps=multirate_fast_substeps,
        multirate_slow_solver=multirate_slow_solver,
        multirate_fast_solver=multirate_fast_solver,
        auto_stiff_gravity_tsit5=auto_stiff_gravity_tsit5,
        auto_stiff_switch_max=auto_stiff_switch_max,
    )
end

"""
    simulation_engine_config_from_env([env=ENV]) -> SimulationEngineConfig

Build a typed `SimulationEngineConfig` from the supported `SPACEAGORA_*`
environment variables. This is the adapter boundary for environment-driven
runtime control.

# Examples
```jldoctest
julia> config = simulation_engine_config_from_env(Dict(
           "SPACEAGORA_PARALLEL_PROFILE" => "R2",
           "SPACEAGORA_SAVE_BUNDLE" => "0",
       ));

julia> (config.parallel.profile, config.artifacts.save_bundle)
("R2", false)
```
"""
function simulation_engine_config_from_env(env::AbstractDict{<:Any, <:Any}=ENV)::SimulationEngineConfig
    parallel = ParallelConfig(
        profile=String(get(env, "SPACEAGORA_PARALLEL_PROFILE", "")),
        outer_parallel_active=_parse_bool(get(env, "SPACEAGORA_OUTER_PARALLEL_ACTIVE", nothing), false),
        parallel_policy_adaptive=_parse_bool(get(env, "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE", nothing), false),
        effector_parallel_mode=String(get(env, "SPACEAGORA_EFFECTOR_PARALLEL", "auto")),
        rhs_batch_parallel_mode=String(get(env, "SPACEAGORA_RHS_BATCH_PARALLEL", "auto")),
        density_callback_parallel_mode=String(get(env, "SPACEAGORA_DENSITY_CALLBACK_PARALLEL", "auto")),
        control_callback_parallel_mode=String(get(env, "SPACEAGORA_CONTROL_CALLBACK_PARALLEL", "auto")),
        thermal_callback_parallel_mode=String(get(env, "SPACEAGORA_THERMAL_CALLBACK_PARALLEL", "auto"))
    )

    env_get = (name, default) -> String(get(env, name, default))
    solver = _solver_config_from_env(env_get)

    runtime_policy = RuntimePolicyConfig(
        warn_normalize=_parse_bool(get(env, "SPACEAGORA_WARN_NORMALIZE", nothing), true),
        allow_typed_normalize=_parse_bool(get(env, "SPACEAGORA_ALLOW_TYPED_NORMALIZE", nothing), false),
        gram_per_sat_instances=_parse_bool(get(env, "SPACEAGORA_GRAM_PER_SAT_INSTANCES", nothing), false),
        srp_ephemeris_cache=_parse_bool(get(env, "SPACEAGORA_SRP_EPHEMERIS_CACHE", nothing), true),
        nbody_ephemeris_cache=_parse_bool(get(env, "SPACEAGORA_NBODY_EPHEMERIS_CACHE", nothing), true),
        planet_frame_cache=_parse_bool(get(env, "SPACEAGORA_PLANET_FRAME_CACHE", nothing), true),
        spice_rhs_memo=_parse_bool(get(env, "SPACEAGORA_SPICE_RHS_MEMO", nothing), true)
    )

    artifacts = ArtifactConfig(
        save_bundle=_parse_bool(get(env, "SPACEAGORA_SAVE_BUNDLE", nothing), true),
        warn_deprecated_config=_parse_bool(get(env, "SPACEAGORA_WARN_DEPRECATED_CONFIG", nothing), true)
    )

    return SimulationEngineConfig(
        parallel=parallel,
        solver=solver,
        runtime_policy=runtime_policy,
        artifacts=artifacts
    )
end

@inline function _engine_env_get(name::String, default::String="")::String
    active_overrides = _engine_active_overrides_ref[]
    if active_overrides !== nothing
        return String(get(active_overrides, name, default))
    end
    return String(get(ENV, name, default))
end

@inline function _engine_env_haskey(name::String)::Bool
    active_overrides = _engine_active_overrides_ref[]
    if active_overrides !== nothing
        return haskey(active_overrides, name)
    end
    return haskey(ENV, name)
end

function _engine_env_overrides(config::SimulationEngineConfig)::Dict{String, String}
    overrides = Dict{String, String}(
        "SPACEAGORA_WARN_NORMALIZE" => _env_bool(config.runtime_policy.warn_normalize),
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE" => _env_bool(config.runtime_policy.allow_typed_normalize),
        "SPACEAGORA_GRAM_PER_SAT_INSTANCES" => _env_bool(config.runtime_policy.gram_per_sat_instances),
        "SPACEAGORA_SRP_EPHEMERIS_CACHE" => _env_bool(config.runtime_policy.srp_ephemeris_cache),
        "SPACEAGORA_NBODY_EPHEMERIS_CACHE" => _env_bool(config.runtime_policy.nbody_ephemeris_cache),
        "SPACEAGORA_PLANET_FRAME_CACHE" => _env_bool(config.runtime_policy.planet_frame_cache),
        "SPACEAGORA_SPICE_RHS_MEMO" => _env_bool(config.runtime_policy.spice_rhs_memo),
        "SPACEAGORA_SAVE_BUNDLE" => _env_bool(config.artifacts.save_bundle),
        "SPACEAGORA_WARN_DEPRECATED_CONFIG" => _env_bool(config.artifacts.warn_deprecated_config),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _env_bool(config.parallel.outer_parallel_active),
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => _env_bool(config.parallel.parallel_policy_adaptive),
        "SPACEAGORA_EFFECTOR_PARALLEL" => config.parallel.effector_parallel_mode,
        "SPACEAGORA_RHS_BATCH_PARALLEL" => config.parallel.rhs_batch_parallel_mode,
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => config.parallel.density_callback_parallel_mode,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => config.parallel.control_callback_parallel_mode,
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => config.parallel.thermal_callback_parallel_mode,
        "SPACEAGORA_SOLVER_MODE" => string(config.solver.solver_mode),
        "SPACEAGORA_SPLIT_IMEX_SOLVER" => string(config.solver.split_imex_solver),
        "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => string(config.solver.multirate_fast_substeps),
        "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => string(config.solver.multirate_slow_solver),
        "SPACEAGORA_MULTIRATE_FAST_SOLVER" => string(config.solver.multirate_fast_solver),
        "SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5" => _env_bool(config.solver.auto_stiff_gravity_tsit5),
        "SPACEAGORA_AUTO_STIFF_SWITCH_MAX" => string(config.solver.auto_stiff_switch_max),
    )

    !isempty(config.parallel.profile) && (overrides["SPACEAGORA_PARALLEL_PROFILE"] = config.parallel.profile)
    !(config.solver.maxiters === nothing) && (overrides["SPACEAGORA_SOLVER_MAXITERS"] = string(config.solver.maxiters))
    !(config.solver.symplectic_dt_s === nothing) && (overrides["SPACEAGORA_SYMPLECTIC_DT_S"] = string(config.solver.symplectic_dt_s))
    !(config.solver.gravity_backbone_dt_s === nothing) && (overrides["SPACEAGORA_GRAVITY_BACKBONE_DT_S"] = string(config.solver.gravity_backbone_dt_s))
    !(config.solver.multirate_slow_dt_s === nothing) && (overrides["SPACEAGORA_MULTIRATE_SLOW_DT_S"] = string(config.solver.multirate_slow_dt_s))

    merge!(overrides, config.env_overrides)
    return overrides
end

function _with_engine_env_overrides(config::SimulationEngineConfig, f::Function)
    overrides = _engine_env_overrides(config)
    previous_config = _engine_active_config_ref[]
    previous_overrides = _engine_active_overrides_ref[]
    _engine_active_config_ref[] = config
    _engine_active_overrides_ref[] = overrides
    isempty(overrides) && return try
        f()
    finally
        _engine_active_config_ref[] = previous_config
        _engine_active_overrides_ref[] = previous_overrides
    end

    previous = Dict{String, Union{Nothing, String}}()
    for (k, v) in overrides
        previous[k] = haskey(ENV, k) ? ENV[k] : nothing
        ENV[k] = String(v)
    end

    try
        return f()
    finally
        for (k, old) in previous
            if old === nothing
                delete!(ENV, k)
            else
                ENV[k] = old
            end
        end
        _engine_active_config_ref[] = previous_config
        _engine_active_overrides_ref[] = previous_overrides
    end
end

@inline _with_engine_env_overrides(f::Function, config::SimulationEngineConfig) =
    _with_engine_env_overrides(config, f)
