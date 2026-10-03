@inline function _coerce_env_bool(v::Bool)::String
    return v ? "1" : "0"
end

@inline function _profile_outer_backend_token(backend::Symbol)::String
    if backend == :none
        return "none"
    elseif backend == :threads
        return "threads"
    elseif backend == :process
        return "process"
    elseif backend == :auto
        return "auto"
    end
    throw(ArgumentError("Unsupported outer backend '$backend'."))
end

@inline function _machine_parallel_class()::Symbol
    override_raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_HARDWARE_CLASS", "auto")))
    if override_raw in ("small", "medium", "large")
        return Symbol(override_raw)
    end
    cpu_threads = Sys.CPU_THREADS
    if cpu_threads >= 24
        return :large
    elseif cpu_threads >= 12
        return :medium
    end
    return :small
end

@inline function _inner_hint_defaults(cfg::ParallelProfileConfig)::NamedTuple{(:exploration, :min_samples), Tuple{Float64, Int}}
    if cfg.profile != R5
        return (exploration=1.5, min_samples=2)
    end
    machine_class = _machine_parallel_class()
    if machine_class == :large
        return (exploration=1.8, min_samples=3)
    elseif machine_class == :medium
        return (exploration=1.5, min_samples=2)
    end
    return (exploration=1.3, min_samples=2)
end

@inline function _env_or_default(name::String, fallback::String; preserve_existing::Bool)::String
    if !preserve_existing
        return fallback
    end
    existing = strip(get(ENV, name, ""))
    return isempty(existing) ? fallback : existing
end

"""
    profile_env_pairs(profile; preserve_existing=true, outer_parallel_active=false)

Return the `SPACEAGORA_*` environment variable pairs implied by a parallel
profile configuration.
"""
function profile_env_pairs(
    profile_in;
    preserve_existing::Bool=true,
    outer_parallel_active::Bool=false
)::Vector{Pair{String, String}}
    cfg = profile_config(profile_in)
    hint_defaults = _inner_hint_defaults(cfg)
    rhs_batch_default = cfg.profile == R0 ? "off" : "auto"
    # R0 is defined as truly serial, so do not preserve a pre-existing RHS batch override.
    rhs_batch_preserve_existing = preserve_existing && cfg.profile != R0
    return Pair{String, String}[
        "SPACEAGORA_PARALLEL_PROFILE" => _env_or_default(
            "SPACEAGORA_PARALLEL_PROFILE",
            parallel_profile_name(cfg.profile);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => _env_or_default(
            "SPACEAGORA_PERF_PARALLEL_BACKEND",
            _profile_outer_backend_token(cfg.outer_backend);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE" => _env_or_default(
            "SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE",
            _coerce_env_bool(cfg.outer_route_adaptive);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _env_or_default(
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE",
            _coerce_env_bool(outer_parallel_active);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE",
            _coerce_env_bool(cfg.inner_adaptive);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => _env_or_default(
            "SPACEAGORA_DENSITY_CALLBACK_PARALLEL",
            cfg.density_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => _env_or_default(
            "SPACEAGORA_CONTROL_CALLBACK_PARALLEL",
            cfg.control_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => _env_or_default(
            "SPACEAGORA_THERMAL_CALLBACK_PARALLEL",
            cfg.thermal_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_MULTIBODY_PARALLEL" => _env_or_default(
            "SPACEAGORA_MULTIBODY_PARALLEL",
            cfg.multibody_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_EFFECTOR_PARALLEL" => _env_or_default(
            "SPACEAGORA_EFFECTOR_PARALLEL",
            cfg.effector_mode;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_RHS_BATCH_PARALLEL" => _env_or_default(
            "SPACEAGORA_RHS_BATCH_PARALLEL",
            rhs_batch_default;
            preserve_existing=rhs_batch_preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER",
            cfg.inner_scheduler;
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_WINDOW" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_WINDOW",
            string(cfg.adaptive_window);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD",
            _coerce_env_bool(cfg.adaptive_control_tail_guard);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD",
            _coerce_env_bool(cfg.adaptive_measured_reward);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS",
            _coerce_env_bool(cfg.persistent_hints);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST",
            _coerce_env_bool(cfg.persistent_state_persist);
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION",
            string(round(hint_defaults.exploration; digits=3));
            preserve_existing=preserve_existing
        ),
        "SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES" => _env_or_default(
            "SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES",
            string(hint_defaults.min_samples);
            preserve_existing=preserve_existing
        )
    ]
end

"""
    with_parallel_profile(f, profile; kwargs...)
    with_parallel_profile(profile, f; kwargs...)

Execute `f` with the environment variables implied by the selected parallel
profile applied via `withenv`.
"""
function with_parallel_profile(
    f::Function,
    profile_in;
    preserve_existing::Bool=true,
    outer_parallel_active::Bool=false
)
    env_pairs = profile_env_pairs(
        profile_in;
        preserve_existing=preserve_existing,
        outer_parallel_active=outer_parallel_active
    )
    return withenv(env_pairs...) do
        f()
    end
end

function with_parallel_profile(
    profile_in,
    f::Function;
    preserve_existing::Bool=true,
    outer_parallel_active::Bool=false
)
    return with_parallel_profile(
        f,
        profile_in;
        preserve_existing=preserve_existing,
        outer_parallel_active=outer_parallel_active
    )
end
