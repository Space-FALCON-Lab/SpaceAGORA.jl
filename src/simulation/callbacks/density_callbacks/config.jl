@inline function _parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline _gram_track_cache_ignore_time_window() = _parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW", true)
@inline _gram_track_cache_target_use_j2() = _parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_TARGET_USE_J2", true)

@inline function _gram_entry_target_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_ENTRY_TARGET_MODE", "allen_eggers")))
    if mode in ("off", "none", "0", "false", "no")
        return :off
    elseif mode in ("allen_eggers", "allen-eggers", "allen", "ae", "on", "1", "true", "yes", "auto")
        return :allen_eggers
    end
    throw(ArgumentError("Unsupported SPACEAGORA_GRAM_ENTRY_TARGET_MODE='$mode'. Use one of: off, allen_eggers."))
end

@inline _gram_entry_target_cd() = max(0.05, _parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_CD", 1.5))
@inline _gram_entry_target_dt() = max(0.05, _parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_DT_S", 0.5))
@inline function _gram_entry_target_max_steps()::Int
    raw = strip(get(ENV, "SPACEAGORA_GRAM_ENTRY_TARGET_MAX_STEPS", "400"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_GRAM_ENTRY_TARGET_MAX_STEPS must be an integer value, got '$raw'"))
    end
    return max(8, parsed)
end

@inline function _density_callback_parallel_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL")
end

@inline function _density_callback_thread_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD", 8)
end

@inline function _density_callback_allow_with_outer()::Bool
    return _parse_bool_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER", false)
end

@inline function _density_batch_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_DENSITY_BATCH_PARALLEL")
end

@inline function _density_batch_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_DENSITY_BATCH_THRESHOLD", 2)
end

@inline function _density_batch_enabled(num_sats::Int)::Bool
    mode = _density_batch_mode()
    if mode == :off
        return false
    elseif mode == :on
        return num_sats > 0
    end
    return num_sats >= _density_batch_threshold()
end

@inline function _gram_isolated_pool_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_GRAM_ISOLATED_POOL"; default="off")
end

@inline function _gram_isolated_pool_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD", 4)
end

@inline function _gram_isolated_pool_max_workers()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS", max(1, Threads.nthreads()))
end

@inline function _gram_isolated_pool_enabled(num_items::Int)::Bool
    mode = _gram_isolated_pool_mode()
    if mode == :off
        return false
    elseif mode == :on
        return num_items > 0
    end
    return Threads.nthreads() > 1 && num_items >= _gram_isolated_pool_threshold()
end

@inline function _control_callback_parallel_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL")
end

@inline function _control_callback_thread_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD", 8)
end

@inline function _control_callback_allow_with_outer()::Bool
    return _parse_bool_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER", false)
end

@inline function _thermal_callback_parallel_mode()::Symbol
    if haskey(ENV, "SPACEAGORA_THERMAL_CALLBACK_PARALLEL")
        return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_THERMAL_CALLBACK_PARALLEL")
    end
    return _density_callback_parallel_mode()
end

@inline function _thermal_callback_thread_threshold()::Int
    if haskey(ENV, "SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD")
        return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD", 8)
    end
    return _density_callback_thread_threshold()
end

@inline function _thermal_callback_allow_with_outer()::Bool
    if haskey(ENV, "SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER")
        return _parse_bool_env("SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER", false)
    end
    return _density_callback_allow_with_outer()
end

@inline function _callback_outer_parallel_hint()::Bool
    return ParallelPolicy.outer_parallel_active()
end

# Extend this for custom user density models as needed:
# SimulationModel.SimulationCallbacks.density_model_threadsafe(::MyDensityModel) = true
@inline density_model_threadsafe(::AbstractDensityModel)::Bool = false
@inline density_model_threadsafe(::NoAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.ExponentialAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.PolynomialFitAtmosphereModel)::Bool = true
# GRAM C-wrapper calls are serialized inside getDensity via RuntimeServices.GRAM_LOCK.
@inline density_model_threadsafe(::EnvironmentModels.GRAMAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.GRAMAtmosphereModelSurrogate)::Bool = true

@inline _is_gram_density_model(model)::Bool =
    model isa EnvironmentModels.GRAMAtmosphereModel ||
    model isa EnvironmentModels.GRAMAtmosphereModelSurrogate

@inline function _gram_track_trajectory_supported(density_model)::Bool
    _is_gram_density_model(density_model) || return false
    hasproperty(density_model, :gram) || return false
    hasproperty(density_model, :gram_atmosphere) || return false
    gram_driver = try
        getproperty(density_model, :gram)
    catch
        return false
    end
    return hasproperty(gram_driver, :generate_trajectory)
end

@inline function _density_callback_thread_decision(args::SimulationConfiguration, num_sats::Int)
    mode = _density_callback_parallel_mode()
    outer_active = _callback_outer_parallel_hint()
    allow_with_outer = _density_callback_allow_with_outer()

    model = args.environment_model.density_model
    model_threadsafe = density_model_threadsafe(model)
    if !model_threadsafe && !_parse_bool_env("SPACEAGORA_DENSITY_CALLBACK_ASSUME_THREADSAFE", false)
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=_density_callback_thread_threshold(),
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:density_callback
    )
    return (use_threads=policy.use_threads, allotment=policy.allotment, mode=mode, policy_applied=true)
end

@inline function _density_callback_use_threads(args::SimulationConfiguration, num_sats::Int)::Bool
    return _density_callback_thread_decision(args, num_sats).use_threads
end

# Extend this for custom user control models as needed:
# SimulationModel.SimulationCallbacks.control_model_threadsafe(::MyControlModel) = true
@inline control_model_threadsafe(::Any)::Bool = false
@inline control_model_threadsafe(::BaseThrusterModel)::Bool = true

@inline function _control_callback_thread_decision(control_model, num_sats::Int, use_invokelatest::Bool)
    if use_invokelatest
        return (use_threads=false, allotment=1, mode=:off, policy_applied=false)
    end
    mode = _control_callback_parallel_mode()
    outer_active = _callback_outer_parallel_hint()
    allow_with_outer = _control_callback_allow_with_outer()

    model_threadsafe = control_model_threadsafe(control_model)
    if !model_threadsafe && !_parse_bool_env("SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE", false)
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=_control_callback_thread_threshold(),
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:control_callback
    )
    return (use_threads=policy.use_threads, allotment=policy.allotment, mode=mode, policy_applied=true)
end

@inline function _control_callback_use_threads(control_model, num_sats::Int, use_invokelatest::Bool)::Bool
    return _control_callback_thread_decision(control_model, num_sats, use_invokelatest).use_threads
end

@inline function _thermal_callback_thread_decision(num_sats::Int)
    mode = _thermal_callback_parallel_mode()
    outer_active = _callback_outer_parallel_hint()
    allow_with_outer = _thermal_callback_allow_with_outer()
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=_thermal_callback_thread_threshold(),
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:thermal_callback
    )
    return (use_threads=policy.use_threads, allotment=policy.allotment, mode=mode, policy_applied=true)
end
