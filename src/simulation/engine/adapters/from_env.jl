@inline _env_bool(v::Bool) = v ? "1" : "0"

function _parse_bool(raw, default::Bool)
    raw === nothing && return default
    token = lowercase(strip(String(raw)))
    token in ("1", "true", "yes", "on") && return true
    token in ("0", "false", "no", "off") && return false
    return default
end

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

    raw_maxiters = strip(String(get(env, "SPACEAGORA_SOLVER_MAXITERS", "")))
    maxiters = isempty(raw_maxiters) ? nothing : try
        parse(Int, raw_maxiters)
    catch
        nothing
    end

    raw_mr_dt = strip(String(get(env, "SPACEAGORA_MULTIRATE_SLOW_DT_S", "")))
    mr_dt = isempty(raw_mr_dt) ? nothing : try
        parse(Float64, raw_mr_dt)
    catch
        nothing
    end

    solver = SolverConfig(
        mode=String(get(env, "SPACEAGORA_SOLVER_MODE", "")),
        maxiters=maxiters,
        split_imex_solver=String(get(env, "SPACEAGORA_SPLIT_IMEX_SOLVER", "kencarp4")),
        multirate_fast_substeps=try
            parse(Int, String(get(env, "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS", "8")))
        catch
            8
        end,
        multirate_slow_dt_s=mr_dt,
        multirate_slow_solver=String(get(env, "SPACEAGORA_MULTIRATE_SLOW_SOLVER", "tsit5")),
        multirate_fast_solver=String(get(env, "SPACEAGORA_MULTIRATE_FAST_SOLVER", "auto_stiff"))
    )

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
        "SPACEAGORA_SPLIT_IMEX_SOLVER" => config.solver.split_imex_solver,
        "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => string(config.solver.multirate_fast_substeps),
        "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => config.solver.multirate_slow_solver,
        "SPACEAGORA_MULTIRATE_FAST_SOLVER" => config.solver.multirate_fast_solver
    )

    !isempty(config.parallel.profile) && (overrides["SPACEAGORA_PARALLEL_PROFILE"] = config.parallel.profile)
    !isempty(config.solver.mode) && (overrides["SPACEAGORA_SOLVER_MODE"] = config.solver.mode)
    !(config.solver.maxiters === nothing) && (overrides["SPACEAGORA_SOLVER_MAXITERS"] = string(config.solver.maxiters))
    !(config.solver.multirate_slow_dt_s === nothing) && (overrides["SPACEAGORA_MULTIRATE_SLOW_DT_S"] = string(config.solver.multirate_slow_dt_s))

    merge!(overrides, config.env_overrides)
    return overrides
end

function _with_engine_env_overrides(config::SimulationEngineConfig, f::Function)
    overrides = _engine_env_overrides(config)
    isempty(overrides) && return f()

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
    end
end
