using Base.ScopedValues: ScopedValue, with

# The routing layer's ParallelProfiles module, found in this module's ancestry
# at call time rather than imported: SimulationEngine is also included on its
# own (the include-order contract suite loads it into a sandbox without the
# routing layer), and only the paths that apply a profile need it. Same lookup
# ParallelCost uses for cgroup_cpu_quota.
function _parallel_profiles_module()::Module
    mod = @__MODULE__
    while true
        isdefined(mod, :ParallelProfiles) && return getproperty(mod, :ParallelProfiles)
        parent = parentmodule(mod)
        parent === mod && break
        mod = parent
    end
    error("ParallelProfiles not found in module ancestry for SimulationEngine.")
end

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
function _parse_float_opt(raw::String, env_name::String)::Union{Nothing, Float64}
    s = strip(raw)
    isempty(s) && return nothing
    v = tryparse(Float64, s)
    v === nothing && throw(ArgumentError("$(env_name) must be a positive number, got '$s'."))
    v > 0.0 || throw(ArgumentError("$(env_name) must be a positive number, got '$s'."))
    return v
end

# `strict=true` (used by `_active_solver_config`, the basis for the `_solver_*`
# unit-test-facing accessors) propagates malformed `SPACEAGORA_SOLVER_*` values as
# ArgumentError, matching this repo's other env parsers. `strict=false` (used by
# `simulation_engine_config_from_env`'s general/introspection contract) instead
# swallows a malformed individual knob and falls back to its default, so a typo'd
# solver knob can't crash construction of the whole SimulationEngineConfig.
_parse_or_default(f::Function, strict::Bool, default) = strict ? f() : (try
    f()
catch e
    e isa ArgumentError ? default : rethrow()
end)

function _solver_config_from_env(env_get=_engine_env_get; strict::Bool=true)::SolverConfig
    solver_mode = _parse_or_default(strict, :tsit5) do
        _parse_solver_mode_sym(env_get("SPACEAGORA_SOLVER_MODE", "tsit5"))
    end

    maxiters = _parse_or_default(strict, nothing) do
        raw_maxiters = strip(env_get("SPACEAGORA_SOLVER_MAXITERS", ""))
        if isempty(raw_maxiters)
            nothing
        else
            v = tryparse(Int, raw_maxiters)
            v === nothing && throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be a positive integer, got '$raw_maxiters'."))
            v > 0 || throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be a positive integer, got '$raw_maxiters'."))
            v
        end
    end

    symplectic_dt_s = _parse_or_default(strict, nothing) do
        _parse_float_opt(env_get("SPACEAGORA_SYMPLECTIC_DT_S", ""), "SPACEAGORA_SYMPLECTIC_DT_S")
    end
    gravity_backbone_dt_s = _parse_or_default(strict, nothing) do
        _parse_float_opt(env_get("SPACEAGORA_GRAVITY_BACKBONE_DT_S", ""), "SPACEAGORA_GRAVITY_BACKBONE_DT_S")
    end

    split_imex_solver = _parse_or_default(strict, :kencarp4) do
        _parse_split_imex_solver_sym(env_get("SPACEAGORA_SPLIT_IMEX_SOLVER", "kencarp4"))
    end

    multirate_slow_dt_s = _parse_or_default(strict, nothing) do
        _parse_float_opt(env_get("SPACEAGORA_MULTIRATE_SLOW_DT_S", ""), "SPACEAGORA_MULTIRATE_SLOW_DT_S")
    end

    multirate_fast_substeps = _parse_or_default(strict, 8) do
        raw_fast_substeps = strip(env_get("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS", "8"))
        v = tryparse(Int, raw_fast_substeps)
        v === nothing && throw(ArgumentError("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS must be a positive integer, got '$raw_fast_substeps'."))
        v > 0 || throw(ArgumentError("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS must be a positive integer, got '$raw_fast_substeps'."))
        v
    end

    multirate_slow_solver = _parse_or_default(strict, :tsit5) do
        _parse_multirate_solver_sym(env_get("SPACEAGORA_MULTIRATE_SLOW_SOLVER", "tsit5"), "SPACEAGORA_MULTIRATE_SLOW_SOLVER")
    end
    multirate_fast_solver = _parse_or_default(strict, :auto_stiff) do
        _parse_multirate_solver_sym(env_get("SPACEAGORA_MULTIRATE_FAST_SOLVER", "auto_stiff"), "SPACEAGORA_MULTIRATE_FAST_SOLVER")
    end

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
function simulation_engine_config_from_env(env::AbstractDict{<:Any, <:Any}=ENV; solver_strict::Bool=false)::SimulationEngineConfig
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
    solver = _solver_config_from_env(env_get; strict=solver_strict)

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

# Adapter variants for knobs that are not part of the canonical override set
# (e.g. the solver SAVE_* switches): consult the active override dict first,
# then fall back to the process environment.  This preserves the historical
# behavior where plain ENV settings for such knobs are honored even inside an
# active SimulationEngineConfig override scope.
@inline function _engine_env_get_with_env_fallback(name::String, default::String)::String
    active_overrides = _engine_active_overrides_ref[]
    if active_overrides !== nothing && haskey(active_overrides, name)
        return String(active_overrides[name])
    end
    return String(get(ENV, name, default))
end

@inline function _engine_env_haskey_with_env_fallback(name::String)::Bool
    active_overrides = _engine_active_overrides_ref[]
    if active_overrides !== nothing && haskey(active_overrides, name)
        return true
    end
    return haskey(ENV, name)
end

const _PARALLEL_CONFIG_DEFAULTS = ParallelConfig()

# The SPACEAGORA_* pairs a ParallelConfig contributes.
#
# Without a profile this is exactly the fixed set it always wrote. With one,
# the profile is expanded the way `with_parallel_profile` would with
# `preserve_existing=false` -- the whole bundle, not just
# SPACEAGORA_PARALLEL_PROFILE, which nothing downstream expands on its own --
# and a field then overrides the profile only where it was set to something
# other than its default, so the defaults no longer clobber the profile's
# adaptive policy and callback modes.
function _parallel_config_env_pairs(parallel::ParallelConfig)::Vector{Pair{String, String}}
    fixed = Pair{String, String}[
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _env_bool(parallel.outer_parallel_active),
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => _env_bool(parallel.parallel_policy_adaptive),
        "SPACEAGORA_EFFECTOR_PARALLEL" => parallel.effector_parallel_mode,
        "SPACEAGORA_RHS_BATCH_PARALLEL" => parallel.rhs_batch_parallel_mode,
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => parallel.density_callback_parallel_mode,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => parallel.control_callback_parallel_mode,
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => parallel.thermal_callback_parallel_mode,
    ]
    isempty(strip(parallel.profile)) && return fixed
    pairs = _parallel_profiles_module().profile_env_pairs(
        parallel.profile;
        preserve_existing=false,
        outer_parallel_active=parallel.outer_parallel_active
    )
    defaults = _PARALLEL_CONFIG_DEFAULTS
    explicit = Pair{String, String}[]
    parallel.parallel_policy_adaptive != defaults.parallel_policy_adaptive &&
        push!(explicit, fixed[2])
    for (i, field) in enumerate((:effector_parallel_mode, :rhs_batch_parallel_mode,
                                 :density_callback_parallel_mode, :control_callback_parallel_mode,
                                 :thermal_callback_parallel_mode))
        getfield(parallel, field) != getfield(defaults, field) && push!(explicit, fixed[2 + i])
    end
    return vcat(pairs, explicit)
end

function _engine_env_overrides(
    config::SimulationEngineConfig;
    parallel_flag::Bool=config.solver.parallel
)::Dict{String, String}
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
        "SPACEAGORA_SOLVER_MODE" => string(config.solver.solver_mode),
        "SPACEAGORA_SPLIT_IMEX_SOLVER" => string(config.solver.split_imex_solver),
        "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS" => string(config.solver.multirate_fast_substeps),
        "SPACEAGORA_MULTIRATE_SLOW_SOLVER" => string(config.solver.multirate_slow_solver),
        "SPACEAGORA_MULTIRATE_FAST_SOLVER" => string(config.solver.multirate_fast_solver),
        "SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5" => _env_bool(config.solver.auto_stiff_gravity_tsit5),
        "SPACEAGORA_AUTO_STIFF_SWITCH_MAX" => string(config.solver.auto_stiff_switch_max),
    )

    for (k, v) in _parallel_config_env_pairs(config.parallel)
        overrides[k] = v
    end
    # SolverConfig(parallel=true) on the engine config: the flag's profile
    # replaces whatever routing the ParallelConfig described. A nested call
    # (inside an enclosing outer split or an already-resolved parallel run)
    # contributes nothing; see `_parallel_flag_applies`.
    if _parallel_flag_applies(parallel_flag)
        prof = strip(config.parallel.profile)
        if !isempty(prof) &&
           _parallel_profiles_module().parse_parallel_profile(prof) != _parallel_profiles_module().PARALLEL_FLAG_PROFILE
            throw(ArgumentError(
                "SolverConfig(parallel=true) selects the parallel settings itself; " *
                "it cannot be combined with ParallelConfig(profile=\"$(prof)\"). " *
                "Set one or the other."
            ))
        end
        for (k, v) in _parallel_profiles_module().parallel_flag_env_pairs()
            overrides[k] = v
        end
    end
    !(config.solver.maxiters === nothing) && (overrides["SPACEAGORA_SOLVER_MAXITERS"] = string(config.solver.maxiters))
    !(config.solver.symplectic_dt_s === nothing) && (overrides["SPACEAGORA_SYMPLECTIC_DT_S"] = string(config.solver.symplectic_dt_s))
    !(config.solver.gravity_backbone_dt_s === nothing) && (overrides["SPACEAGORA_GRAVITY_BACKBONE_DT_S"] = string(config.solver.gravity_backbone_dt_s))
    !(config.solver.multirate_slow_dt_s === nothing) && (overrides["SPACEAGORA_MULTIRATE_SLOW_DT_S"] = string(config.solver.multirate_slow_dt_s))

    merge!(overrides, config.env_overrides)
    return overrides
end

function _with_engine_env_overrides(
    config::SimulationEngineConfig,
    f::Function;
    parallel_flag::Bool=config.solver.parallel
)
    overrides = _engine_env_overrides(config; parallel_flag=parallel_flag)
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

@inline _with_engine_env_overrides(f::Function, config::SimulationEngineConfig; kwargs...) =
    _with_engine_env_overrides(config, f; kwargs...)

# ── SolverConfig(parallel=true) ──────────────────────────────────────────────
#
# The flag's whole effect is an environment scope: the profile it names
# (ParallelProfiles.PARALLEL_FLAG_PROFILE) applied around one run or campaign
# and restored afterwards, exception or not. It lives here because this file is
# the engine's only sanctioned reader and writer of ENV.

# True inside a run or campaign whose flag has already been resolved, so a
# member run of a parallel campaign, or the typed run inside an engine-config
# scope, does not apply it (or calibrate) a second time. Scoped, so it follows
# the tasks a campaign spawns and nothing else.
const _PARALLEL_FLAG_RESOLVED = ScopedValue(false)
const _PARALLEL_FLAG_ONE_THREAD_NOTED = Base.Threads.Atomic{Bool}(false)

"""
    _parallel_flag_applies(flag) -> Bool

Whether `SolverConfig(parallel=flag)` should apply its profile here: the flag is
set, no enclosing run or campaign has resolved it already, and this is not a
nested call inside an enclosing outer split (a threaded split's
`SPACEAGORA_OUTER_PARALLEL_ACTIVE`, or a process-pool worker), which keeps the
split's own environment exactly as it did before the flag existed.
"""
@inline function _parallel_flag_applies(flag::Bool)::Bool
    flag || return false
    _PARALLEL_FLAG_RESOLVED[] && return false
    return !_parallel_profiles_module().parallel_flag_nested()
end

"""
    _parallel_flag_prepare!()

One-time work before the first parallel run: create this machine's cost
constants if they are missing (never when a current file exists), and say once
when Julia has a single thread, since the flag then has only the process pool
to use.
"""
function _parallel_flag_prepare!()::Nothing
    SimulationModel.ParallelCost.ensure_machine_constants!()
    if Base.Threads.nthreads() == 1 && !Base.Threads.atomic_xchg!(_PARALLEL_FLAG_ONE_THREAD_NOTED, true)
        @info "SolverConfig(parallel=true) with a single Julia thread: only process-worker " *
              "campaign routes are available. Start Julia with `--threads=auto` to use threads."
    end
    return nothing
end

"""
    _with_parallel_flag(f, flag::Bool)

Run `f()` under the parallel flag's environment when
[`_parallel_flag_applies`](@ref)`(flag)`, restoring every variable (and any
active engine-config override set) afterwards; otherwise just call `f()`.
"""
function _with_parallel_flag(f::Function, flag::Bool)
    _parallel_flag_applies(flag) || return (flag ? with(f, _PARALLEL_FLAG_RESOLVED => true) : f())
    _parallel_flag_prepare!()
    pairs = _parallel_profiles_module().parallel_flag_env_pairs()
    previous_overrides = _engine_active_overrides_ref[]
    if previous_overrides !== nothing
        # Inside a SimulationEngineConfig scope the engine's own reads go to
        # the override set, so it must carry the flag's values too.
        merged = copy(previous_overrides)
        for (k, v) in pairs
            merged[k] = v
        end
        _engine_active_overrides_ref[] = merged
    end
    try
        return withenv(pairs...) do
            with(f, _PARALLEL_FLAG_RESOLVED => true)
        end
    finally
        _engine_active_overrides_ref[] = previous_overrides
    end
end
