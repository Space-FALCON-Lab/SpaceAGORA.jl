@inline function _parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function _gram_track_cache_ignore_time_window()::Bool
    return _parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW", true)
end

@inline function _gram_track_cache_target_use_j2()::Bool
    return _parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_TARGET_USE_J2", true)
end

# See CallbackEnvConfig.density_freeze_per_step docstring for the rationale.
@inline function _density_freeze_per_step_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_DENSITY_FREEZE_PER_STEP", false)
end

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

# Per-satellite cost class for the density callback's thread decision.
#
# `thread_policy_decision` already carries a light-work guard -- `heavy_only`
# with `heavy_work` false pins `use_threads` to false -- and the effector policy
# has used it from the start, for exactly this reason: a dispatch that costs
# more than the work it hands out makes the callback slower, not faster. The
# density callback never passed it. Its decision asked only "are there enough
# satellites", so eight spacecraft on an analytic atmosphere threaded a loop
# whose body is a 3x3 rotation, a lat/lon conversion and one `exp`, once per
# accepted step, through a persistent-pool round trip per worker.
#
# Heavy means the per-satellite body can reach a native GRAM evaluation or a
# GRAM track-cache/surrogate lookup -- tens of microseconds and up, which is
# what the threaded path was built for and where it is measured to win. Every
# closed-form atmosphere is light. So is the batch pre-fill loop whatever the
# model is: that loop only stages altitude/latitude/longitude and the density
# evaluation itself happens afterwards, on one thread, inside
# `getDensityBatch!`.
#
# `SPACEAGORA_DENSITY_CALLBACK_PARALLEL=on` is unaffected -- an explicit `on`
# forces threads ahead of the guard -- so a caller that wants the dispatch
# measured on light work can still ask for it.
@inline density_model_work_is_heavy(::AbstractDensityModel)::Bool = false
@inline density_model_work_is_heavy(::EnvironmentModels.GRAMAtmosphereModel)::Bool = true
@inline density_model_work_is_heavy(::EnvironmentModels.GRAMAtmosphereModelSurrogate)::Bool = true

"""
    _density_callback_work_is_heavy(p, num_sats) -> Bool

True when at least one satellite's density evaluation is expensive enough to be
worth a threaded dispatch.  Reads the per-satellite model vector when the run
installed one and the configured model otherwise, and treats a run with the
vacuum-predicted GRAM cache enabled as heavy: that path rebuilds a spline over
the look-ahead trajectory inside the callback body.
"""
@inline function _density_callback_work_is_heavy(p, num_sats::Int)::Bool
    env = _callback_env_config(p)
    env.vacuum_gram_cache_enabled && return true
    if p !== nothing && hasproperty(p, :shared_buffers)
        models = p.shared_buffers.density_models
        if !isempty(models)
            limit = min(num_sats, length(models))
            @inbounds for i in 1:limit
                density_model_work_is_heavy(models[i]) && return true
            end
            num_sats <= length(models) && return false
        end
    end
    p === nothing && return false
    return density_model_work_is_heavy(p.args.environment_model.density_model)
end

# Extend this for custom user density models as needed:
# SimulationModel.SimulationCallbacks.density_model_threadsafe(::MyDensityModel) = true
@inline density_model_threadsafe(::AbstractDensityModel)::Bool = false
@inline density_model_threadsafe(::NoAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.ExponentialAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.PiecewiseExponentialAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.PolynomialFitAtmosphereModel)::Bool = true
# GRAM C-wrapper calls are serialized inside getDensity via RuntimeServices.GRAM_LOCK.
@inline density_model_threadsafe(::EnvironmentModels.GRAMAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.GRAMAtmosphereModelSurrogate)::Bool = true

@inline _is_gram_density_model(model)::Bool =
    model isa EnvironmentModels.GRAMAtmosphereModel ||
    model isa EnvironmentModels.GRAMAtmosphereModelSurrogate

# The wrapper types forward properties to the GRAMSuite core, whose `gram`
# field holds the native GRAM Julia wrapper module. `hasproperty` on a Module
# only sees exported names, so module drivers are probed with `isdefined` —
# the same capability check GRAMSuite itself uses (e.g. for `get_winds_state`).
@inline function _gram_track_trajectory_supported(density_model)::Bool
    _is_gram_density_model(density_model) || return false
    hasproperty(density_model, :gram) || return false
    hasproperty(density_model, :gram_atmosphere) || return false
    gram_driver = try
        getproperty(density_model, :gram)
    catch
        return false
    end
    if gram_driver isa Module
        return isdefined(gram_driver, :generate_trajectory)
    end
    return hasproperty(gram_driver, :generate_trajectory)
end

"""
    _snapshot_callback_env_config() -> CallbackEnvConfig

Resolve every env-derived knob consulted per callback invocation (and per
RHS-side atmosphere sample) into a typed snapshot.  Built once at
run_simulation setup; hot paths read plain struct fields via
`_callback_env_config(p)` instead of re-parsing ENV.
"""
function _snapshot_callback_env_config()::CallbackEnvConfig
    return CallbackEnvConfig(
        _gram_track_cache_config(),
        _gram_runtime_stats_enabled(),
        _gram_track_cache_ignore_time_window(),
        _gram_track_cache_target_use_j2(),
        _density_freeze_per_step_enabled(),
        _vacuum_gram_cache_enabled(),
        _vacuum_gram_cache_npoints(),
        _vacuum_gram_cache_horizon_s(),
        _vacuum_gram_cache_deviation_m(),
        _density_callback_parallel_mode(),
        _density_callback_thread_threshold(),
        _density_callback_allow_with_outer(),
        _parse_bool_env("SPACEAGORA_DENSITY_CALLBACK_ASSUME_THREADSAFE", false),
        _density_batch_mode(),
        _density_batch_threshold(),
        _gram_isolated_pool_mode(),
        _gram_isolated_pool_threshold(),
        _gram_isolated_pool_max_workers(),
        _control_callback_parallel_mode(),
        _control_callback_thread_threshold(),
        _control_callback_allow_with_outer(),
        _parse_bool_env("SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE", false),
        _thermal_callback_parallel_mode(),
        _thermal_callback_thread_threshold(),
        _thermal_callback_allow_with_outer(),
    )
end

# Run-scoped snapshot accessor.  Falls back to live ENV parsing when the
# snapshot is unset (hand-constructed ODEParams in unit tests / withenv probes).
@inline function _callback_env_config(p)::CallbackEnvConfig
    if p !== nothing && hasproperty(p, :shared_buffers)
        sb = getproperty(p, :shared_buffers)
        if hasproperty(sb, :callback_env_config)
            cfg = sb.callback_env_config[]
            cfg === nothing || return cfg
        end
    end
    return _snapshot_callback_env_config()
end

# Policy snapshot accessor: `nothing` (→ live reads in thread_policy_decision)
# when the run has not installed a snapshot.
@inline function _policy_env_config(p)::Union{Nothing, PolicyDecisionEnvConfig}
    if p !== nothing && hasproperty(p, :shared_buffers)
        sb = getproperty(p, :shared_buffers)
        if hasproperty(sb, :policy_env_config)
            return sb.policy_env_config[]
        end
    end
    return nothing
end

@inline function _density_batch_enabled(env::CallbackEnvConfig, num_sats::Int)::Bool
    mode = env.density_batch_mode
    if mode == :off
        return false
    elseif mode == :on
        return num_sats > 0
    end
    return num_sats >= env.density_batch_threshold
end

@inline function _gram_isolated_pool_enabled(env::CallbackEnvConfig, num_items::Int)::Bool
    mode = env.gram_isolated_pool_mode
    if mode == :off
        return false
    elseif mode == :on
        return num_items > 0
    end
    return Threads.nthreads() > 1 && num_items >= env.gram_isolated_pool_threshold
end

# `heavy_work` defaults to true because the cost of the loop body is a property
# of the call site, not of this function: the callback's batch route stages
# kinematics only, its per-satellite route runs a full density evaluation, and
# the RHS-side atmosphere pre-fill (dynamics_rhs.jl) samples density inline. A
# caller that knows its body is light says so; the default answers the older,
# narrower question -- would the policy thread this if the work were worth
# threading -- and so leaves every existing call site's behavior unchanged.
@inline function _density_callback_thread_decision(
    args::SimulationConfiguration,
    num_sats::Int;
    heavy_work::Bool=true
)
    return _density_callback_thread_decision(nothing, args, num_sats; heavy_work=heavy_work)
end

@inline function _density_callback_thread_decision(
    p,
    args::SimulationConfiguration,
    num_sats::Int;
    heavy_work::Bool=true
)
    env = _callback_env_config(p)
    penv = _policy_env_config(p)
    mode = env.density_parallel_mode
    # A width pinned by the pre-solve sweep (V2) short-circuits the per-call
    # decision. It is only ever set after the static decision below said
    # "thread this", so the thread-safety check has already passed for it.
    if p !== nothing && hasproperty(p, :shared_buffers)
        sb = getproperty(p, :shared_buffers)
        if hasproperty(sb, :density_callback_width)
            w = sb.density_callback_width[]
            if w > 0
                return (use_threads=w > 1, allotment=w, mode=mode, policy_applied=false)
            end
        end
    end
    outer_active = penv === nothing ? _callback_outer_parallel_hint() : penv.outer_parallel_active
    allow_with_outer = env.density_allow_with_outer

    model = args.environment_model.density_model
    model_threadsafe = density_model_threadsafe(model)
    if !model_threadsafe && !env.density_assume_threadsafe
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    # Native/point GRAM is serialized behind a process-wide lock (GRAM_LOCK), so
    # oversubscribing it below a reasonably high thread count wastes cycles
    # fighting for that lock -- the :density_callback source's 16-thread floor
    # exists for that case. A lock-free model (e.g. GRAMAtmosphereModelSurrogate)
    # has no such cost, so it gets the general default floor instead via a
    # separate source category, rather than being held to the same 16-thread gate
    # for no reason (see PARALLELIZATION_CURRENT_STATE.md / Finding 1).
    source = model isa EnvironmentModels.GRAMAtmosphereModel ? :density_callback : :density_callback_lockfree
    # heavy_only is passed unconditionally: the guard only bites when the caller
    # says the per-satellite body is light, and `:on` overrides it either way.
    # See density_model_work_is_heavy for what "light" costs here.
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=env.density_thread_threshold,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        heavy_only=true,
        heavy_work=heavy_work,
        source=source,
        env=penv
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

@inline function _control_callback_thread_decision(control_model, num_sats::Int)
    return _control_callback_thread_decision(nothing, control_model, num_sats)
end

@inline function _control_callback_thread_decision(p, control_model, num_sats::Int)
    env = _callback_env_config(p)
    penv = _policy_env_config(p)
    mode = env.control_parallel_mode
    outer_active = penv === nothing ? _callback_outer_parallel_hint() : penv.outer_parallel_active
    allow_with_outer = env.control_allow_with_outer

    model_threadsafe = control_model_threadsafe(control_model)
    if !model_threadsafe && !env.control_assume_threadsafe
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=env.control_thread_threshold,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:control_callback,
        env=penv
    )
    return (use_threads=policy.use_threads, allotment=policy.allotment, mode=mode, policy_applied=true)
end

@inline function _control_callback_use_threads(control_model, num_sats::Int)::Bool
    return _control_callback_thread_decision(control_model, num_sats).use_threads
end

# `heavy_work` is a keyword with the same contract as the density version, and
# for the same reason: the thermal callback's per-satellite body reads its
# density from `shared_buffers` rather than evaluating a model
# (`_compute_stage_heat_rates!` is called with `use_buffered_density=true`), so
# what it costs is one `sample_planet_frame` plus one `getHeatRate` per link --
# light for a single-link spacecraft.
#
# Measured on the same 8-spacecraft one-hour shape at 24 threads, alternating
# arms in one process: density and thermal both off, 0.61-0.67 s; density auto
# and thermal off, 0.57-0.62 s; thermal auto, 1.14-2.04 s whichever way density
# is set. All twelve runs produced identical final states. So after the density
# guard the thermal callback is what is left of this shape's inner-threading
# cost, and `get_thermal_callback`'s dispatch is the call site that would have
# to pass its own classification here -- one keyword at
# src/simulation/callbacks/thermal_callbacks.jl:82, in a file this change does
# not own. The default stays `true` so that call site's behavior is unchanged
# until someone measures what link count makes the dispatch worth paying for.
@inline function _thermal_callback_thread_decision(num_sats::Int; heavy_work::Bool=true)
    return _thermal_callback_thread_decision(nothing, num_sats; heavy_work=heavy_work)
end

@inline function _thermal_callback_thread_decision(p, num_sats::Int; heavy_work::Bool=true)
    env = _callback_env_config(p)
    penv = _policy_env_config(p)
    mode = env.thermal_parallel_mode
    outer_active = penv === nothing ? _callback_outer_parallel_hint() : penv.outer_parallel_active
    allow_with_outer = env.thermal_allow_with_outer
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=env.thermal_thread_threshold,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        heavy_only=true,
        heavy_work=heavy_work,
        source=:thermal_callback,
        env=penv
    )
    return (use_threads=policy.use_threads, allotment=policy.allotment, mode=mode, policy_applied=true)
end
