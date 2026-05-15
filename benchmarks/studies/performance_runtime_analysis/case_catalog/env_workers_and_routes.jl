function _inner_hint_layer_report_df(spec::ProfileSpec, hw)::DataFrame
    rows = SimulationModel.ParallelPolicy.hint_layer_stats_snapshot(
        profile=_active_parallel_profile_token(),
        machine=_active_parallel_machine_token()
    )
    stamp_utc = string(now(UTC))
    if isempty(rows)
        return DataFrame(
            runtime_profile=String[],
            parallel_profile=String[],
            hint_machine=String[],
            hardware_class=String[],
            machine_label=String[],
            host_name=String[],
            layer=String[],
            signature_count=Int[],
            choice_count=Int[],
            samples_total=Int[],
            successes_total=Int[],
            failures_total=Int[],
            elapsed_mean_ns=Float64[],
            elapsed_std_ns=Float64[],
            confidence_mean=Float64[],
            regret_mean_ns=Float64[],
            exploration_c=Float64[],
            min_samples=Int[],
            state_path=String[],
            timestamp_utc=String[]
        )
    end
    return DataFrame([
        (
            runtime_profile=spec.name,
            parallel_profile=String(row.profile),
            hint_machine=String(row.machine),
            hardware_class=hw.hardware_class,
            machine_label=hw.machine_label,
            host_name=hw.host_name,
            layer=String(row.layer),
            signature_count=Int(row.signature_count),
            choice_count=Int(row.choice_count),
            samples_total=Int(row.samples_total),
            successes_total=Int(row.successes_total),
            failures_total=Int(row.failures_total),
            elapsed_mean_ns=Float64(row.elapsed_mean_ns),
            elapsed_std_ns=Float64(row.elapsed_std_ns),
            confidence_mean=Float64(row.confidence_mean),
            regret_mean_ns=Float64(row.regret_mean_ns),
            exploration_c=Float64(row.exploration_c),
            min_samples=Int(row.min_samples),
            state_path=String(row.state_path),
            timestamp_utc=stamp_utc
        )
        for row in rows
    ])
end

@inline function _total_links(case::BenchmarkCase)::Int
    total = 0
    @inbounds for sc in case.args_template.dynamics_model.spacecraft
        total += length(sc.links)
    end
    return total
end

@inline function _has_atmo_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa AerodynamicCoefficientConstant || effector isa AerodynamicCoefficientfM || effector isa AerodynamicCoefficientNoBallisticFlight
            return true
        end
    end
    return false
end

@inline function _has_nbody_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa NBodyGravityModel
            return true
        end
    end
    return false
end

@inline function _has_srp_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa SolarRadiationPressureModel && effector.A > 0.0
            return true
        end
    end
    return false
end

@inline function _max_harmonics_degree(effectors::Tuple)::Int
    degree = 0
    @inbounds for effector in effectors
        if effector isa GravitationalHarmonicsModel
            degree = max(degree, effector.L)
        end
    end
    return degree
end

@inline function _max_links_per_sat(case::BenchmarkCase)::Int
    max_links = 1
    @inbounds for sc in case.args_template.dynamics_model.spacecraft
        max_links = max(max_links, length(sc.links))
    end
    return max_links
end

@inline function _density_model_family(density_model)::String
    if density_model isa NoAtmosphereModel
        return "none"
    elseif density_model isa GRAMAtmosphereModelSurrogate
        return "gram_surrogate"
    elseif density_model isa GRAMAtmosphereModel
        return "gram_point"
    elseif density_model isa ExponentialAtmosphereModel
        return "exponential"
    elseif density_model isa SimulationModel.EnvironmentModels.PolynomialFitAtmosphereModel
        return "polyfit"
    elseif density_model isa NRLMSISE00AtmosphereModel
        return "nrlmsise00"
    end
    return lowercase(string(nameof(typeof(density_model))))
end

@inline function _read_bool_property_safe(obj, name::Symbol)::Union{Nothing, Bool}
    if !hasproperty(obj, name)
        return nothing
    end
    value = try
        getproperty(obj, name)
    catch
        return nothing
    end
    return value isa Bool ? value : nothing
end

@inline function _gram_surrogate_flag(density_model)::Bool
    if density_model isa GRAMAtmosphereModelSurrogate
        return true
    end
    # Fallback for custom wrappers exposing a boolean surrogate toggle.
    flag = _read_bool_property_safe(density_model, :gram_offline_surrogate)
    return isnothing(flag) ? false : flag
end

@inline function _gram_static_grid_flag(density_model)::Bool
    if !(density_model isa GRAMAtmosphereModel || density_model isa GRAMAtmosphereModelSurrogate)
        return false
    end
    for key in (:use_static_grid, :static_grid, :gram_static_grid)
        flag = _read_bool_property_safe(density_model, key)
        if !(flag === nothing)
            return flag
        end
    end
    return false
end

@inline function _env_bool_token(raw::AbstractString)::Union{Nothing, Bool}
    token = lowercase(strip(String(raw)))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    return nothing
end

@inline function _case_env_value(case::BenchmarkCase, name::String)::Union{Nothing, String}
    for pair in case.env_overrides
        if first(pair) == name
            return String(last(pair))
        end
    end
    return nothing
end

@inline function _env_bool_override_for_case(case::BenchmarkCase, name::String)::Union{Nothing, Bool}
    value = _case_env_value(case, name)
    value === nothing && return nothing
    return _env_bool_token(value)
end

@inline function _gram_static_grid_flag_for_case(case::BenchmarkCase, density_model)::Bool
    if !(density_model isa GRAMAtmosphereModel || density_model isa GRAMAtmosphereModelSurrogate)
        return false
    end
    override = _env_bool_override_for_case(case, "SPACEAGORA_GRAM_STATIC_GRID")
    if !(override === nothing)
        return override
    end
    return _gram_static_grid_flag(density_model)
end

@inline function _gram_track_cache_mode_for_case(case::BenchmarkCase)::String
    override = _case_env_value(case, "SPACEAGORA_GRAM_TRACK_CACHE")
    if !(override === nothing)
        token = lowercase(strip(override))
        if token in ("on", "off", "auto")
            return token
        end
    end
    return lowercase(strip(get(ENV, "SPACEAGORA_GRAM_TRACK_CACHE", "off")))
end

@inline function _density_backend_bucket(
    density_family::AbstractString,
    gram_surrogate_enabled::Bool,
    gram_static_grid_enabled::Bool,
    gram_track_cache_mode::AbstractString
)::String
    family = lowercase(strip(String(density_family)))
    cache_on = lowercase(strip(String(gram_track_cache_mode))) in ("on", "auto")
    if family in ("gram_point", "gram")
        return gram_static_grid_enabled ? "gram_static_grid_or_cached_surrogate" : "gram_point_to_point"
    elseif family in ("gram_surrogate", "gram_offline_surrogate")
        if gram_static_grid_enabled || cache_on
            return "gram_static_grid_or_cached_surrogate"
        end
        return gram_surrogate_enabled ? "gram_surrogate" : "gram_point_to_point"
    elseif family == "none"
        return "non_gram"
    elseif occursin("gram", family)
        return gram_static_grid_enabled ? "gram_static_grid_or_cached_surrogate" : "gram_point_to_point"
    end
    return "non_gram"
end

@inline function _solver_mode_for_outer_route(case::BenchmarkCase, profile_name::String)::String
    if !(case.solver_mode_override === nothing)
        return String(case.solver_mode_override)
    end
    if !(case.split_imex_solver_override === nothing)
        return "split_imex:" * String(case.split_imex_solver_override)
    end
    return _perf_solver_mode_env(profile_name)
end

@inline function _min_positive_rate(rates::AbstractVector{<:Real})::Float64
    best = Inf
    @inbounds for rate in rates
        value = Float64(rate)
        if isfinite(value) && value > 0.0
            best = min(best, value)
        end
    end
    return isfinite(best) ? best : 0.0
end

@inline function _has_aero_dynamic_effector(effectors::Tuple)::Bool
    return _has_atmo_dynamic_effector(effectors)
end

@inline function _dynamic_effector_cost_class(effectors::Tuple, control_effector_count::Int)::String
    effector_count = length(effectors)
    has_nbody = _has_nbody_dynamic_effector(effectors)
    has_srp = _has_srp_dynamic_effector(effectors)
    has_aero = _has_aero_dynamic_effector(effectors)
    harmonics_degree = _max_harmonics_degree(effectors)

    if has_nbody || harmonics_degree >= 20 || effector_count >= 5
        return "heavy"
    end
    if has_srp || has_aero || harmonics_degree > 0 || effector_count >= 2 || control_effector_count >= 2
        return "medium"
    end
    return "light"
end

@inline function _case_outer_threads_safe(case::BenchmarkCase)::Bool
    # Guard against native-library aborts seen under outer-threaded case concurrency.
    density_model = case.args_template.environment_model.density_model
    if density_model isa GRAMAtmosphereModel || density_model isa GRAMAtmosphereModelSurrogate
        return false
    end
    return !_has_nbody_dynamic_effector(case.args_template.dynamics_model.dynamic_effectors)
end

function _split_threadsafe_indices(
    selected::Vector{BenchmarkCase},
    indices::Vector{Int}
)::Tuple{Vector{Int}, Vector{Int}}
    thread_indices = Int[]
    serial_indices = Int[]
    for idx in indices
        if _case_outer_threads_safe(selected[idx])
            push!(thread_indices, idx)
        else
            push!(serial_indices, idx)
        end
    end
    return thread_indices, serial_indices
end

function _reset_outer_route_history!()
    _PERF_OUTER_ROUTE_STATE[] = ParallelProfiles.OuterRouteState()
    return nothing
end

function _load_outer_route_history!(spec::ProfileSpec, outdir::String)
    _reset_outer_route_history!()
    persist = _outer_route_state_persist_enabled() && _outer_route_adaptive_enabled()
    path = _outer_route_state_cache_path(spec, outdir)
    if !persist
        return (persist=false, reset_requested=false, path=path, loaded_rows=0, loaded_signatures=0)
    end
    reset_requested = _outer_route_state_reset_requested()
    if reset_requested
        return (persist=true, reset_requested=true, path=path, loaded_rows=0, loaded_signatures=0)
    end
    loaded = try
        ParallelProfiles.load_outer_route_state!(_PERF_OUTER_ROUTE_STATE[], path; replace=true)
    catch err
        @warn "[perf] Failed to load outer-route adaptive cache; starting with empty history." path=path exception=(err, catch_backtrace())
        (path=path, rows=0, signatures=0)
    end
    return (
        persist=true,
        reset_requested=false,
        path=loaded.path,
        loaded_rows=loaded.rows,
        loaded_signatures=loaded.signatures
    )
end

function _save_outer_route_history!(spec::ProfileSpec, path::String, enabled::Bool)
    enabled || return (path=path, rows=0, signatures=0)
    metadata = Dict{String, Any}(
        "profile" => spec.name,
        "machine_label" => _machine_label(),
        "hardware_class" => _hardware_class_name(),
        "cpu_threads" => Sys.CPU_THREADS,
        "julia_threads" => Threads.nthreads()
    )
    saved = try
        ParallelProfiles.save_outer_route_state(
            _PERF_OUTER_ROUTE_STATE[],
            path;
            metadata=metadata
        )
    catch err
        @warn "[perf] Failed to persist outer-route adaptive cache." path=path exception=(err, catch_backtrace())
        return (path=path, rows=0, signatures=0)
    end
    return saved
end

@inline function _outer_route_tuning()::ParallelProfiles.OuterRouteTuning
    return ParallelProfiles.OuterRouteTuning(
        inner_sat_threshold=_priority_inner_sat_threshold(),
        inner_link_threshold=_priority_inner_link_threshold(),
        outer_light_sat_threshold=_priority_outer_light_sat_threshold(),
        outer_light_link_threshold=_priority_outer_light_link_threshold(),
        outer_light_mission_threshold_s=_priority_outer_light_mission_threshold_s(),
        spice_constellation_process_enabled=_priority_spice_constellation_process_enabled(),
        spice_constellation_min_sats=_priority_spice_constellation_min_sats(),
        adaptive_enabled=_outer_route_adaptive_enabled(),
        adaptive_min_samples=_outer_route_min_samples(),
        failure_penalty_s=_outer_route_failure_penalty_s(),
        mc_process_min_samples=_outer_route_mc_process_min_samples(),
        mc_process_min_mission_s=_outer_route_mc_process_min_mission_s(),
        trace=_outer_route_trace_enabled()
    )
end

@inline function _outer_route_features(
    case::BenchmarkCase;
    spec::Union{Nothing, ProfileSpec}=nothing
)::ParallelProfiles.OuterRouteFeatures
    profile_name = isnothing(spec) ? lowercase(strip(get(ENV, "SPACEAGORA_PERF_PROFILE", "quick"))) : spec.name
    n_sats = length(case.args_template.dynamics_model.spacecraft)
    n_links = _total_links(case)
    max_links_per_sat = _max_links_per_sat(case)
    mission_time_s = case.args_template.mission_configuration.mission_time
    dt_max_orbit_s = case.args_template.integration_tolerances.dt_max_orbit
    dynamic_effectors = case.args_template.dynamics_model.dynamic_effectors
    density_model = case.args_template.environment_model.density_model
    has_nbody = _has_nbody_dynamic_effector(dynamic_effectors)
    has_srp = _has_srp_dynamic_effector(dynamic_effectors)
    harmonics_degree = _max_harmonics_degree(dynamic_effectors)
    has_control = !isempty(case.args_template.control_model.control_effectors)
    control_effector_count = length(case.args_template.control_model.control_effectors)
    orientation_on = case.args_template.mission_configuration.orientation_sim
    density_family = _density_model_family(density_model)
    solver_mode = _solver_mode_for_outer_route(case, profile_name)
    control_rate_s = _min_positive_rate(case.args_template.control_model.control_rates)
    guidance_rate_s = _min_positive_rate(case.args_template.guidance_model.guidance_rates)
    navigation_rate_s = _min_positive_rate(case.args_template.navigation_model.navigation_rates)
    gram_surrogate_enabled = _gram_surrogate_flag(density_model)
    gram_static_grid_enabled = _gram_static_grid_flag_for_case(case, density_model)
    thermal_enabled = _has_atmo_dynamic_effector(dynamic_effectors) || !(density_model isa NoAtmosphereModel)
    dynamic_effector_count = length(dynamic_effectors)
    effector_cost_class = _dynamic_effector_cost_class(dynamic_effectors, control_effector_count)
    mc_samples = isnothing(spec) ? 0 : spec.montecarlo_samples
    return ParallelProfiles.OuterRouteFeatures(
        category=case.category,
        n_sats=n_sats,
        n_links=n_links,
        max_links_per_sat=max_links_per_sat,
        mission_time_s=mission_time_s,
        has_nbody=has_nbody,
        has_srp=has_srp,
        harmonics_degree=harmonics_degree,
        has_control=has_control,
        orientation_on=orientation_on,
        density_family=density_family,
        solver_mode=solver_mode,
        dt_max_orbit_s=dt_max_orbit_s,
        control_rate_s=control_rate_s,
        guidance_rate_s=guidance_rate_s,
        navigation_rate_s=navigation_rate_s,
        gram_surrogate_enabled=gram_surrogate_enabled,
        gram_static_grid_enabled=gram_static_grid_enabled,
        control_effector_count=control_effector_count,
        thermal_enabled=thermal_enabled,
        dynamic_effector_count=dynamic_effector_count,
        effector_cost_class=effector_cost_class,
        montecarlo_samples=mc_samples
    )
end

@inline function _priority_outer_route_montecarlo(case::BenchmarkCase, spec::Union{Nothing, ProfileSpec}=nothing)::Symbol
    features = _outer_route_features(case; spec=spec)
    return ParallelProfiles.default_outer_route(
        features;
        tuning=_outer_route_tuning(),
        machine_class=_machine_parallel_class(),
        threads_available=perf_threads_available(),
        parallel_enabled=perf_parallel_enabled()
    )
end

@inline function _outer_route_candidates(case::BenchmarkCase; spec::Union{Nothing, ProfileSpec}=nothing)::Vector{Symbol}
    features = _outer_route_features(case; spec=spec)
    return ParallelProfiles.outer_route_candidates(
        features;
        tuning=_outer_route_tuning(),
        machine_class=_machine_parallel_class(),
        threads_available=perf_threads_available(),
        parallel_enabled=perf_parallel_enabled()
    )
end

@inline function _outer_route_signature(case::BenchmarkCase)::String
    return ParallelProfiles.outer_route_signature(_outer_route_features(case))
end

@inline function _outer_route_stats_snapshot(signature::String)
    return ParallelProfiles.outer_route_stats_snapshot(_PERF_OUTER_ROUTE_STATE[], signature)
end

function _record_outer_route_feedback!(
    case::BenchmarkCase,
    rows::AbstractVector{<:NamedTuple};
    route::Symbol
)
    features = _outer_route_features(case)
    successes = 0
    failures = 0
    elapsed_success_s = 0.0
    elapsed_success_sq_sum_s = 0.0
    for row in rows
        if !hasproperty(row, :solve_success)
            continue
        end
        if row.solve_success === true
            if hasproperty(row, :total_time_s) && row.total_time_s isa Real && isfinite(Float64(row.total_time_s))
                elapsed_s = Float64(row.total_time_s)
                successes += 1
                elapsed_success_s += elapsed_s
                elapsed_success_sq_sum_s += elapsed_s^2
            end
        elseif row.solve_success === false
            failures += 1
        end
    end
    ParallelProfiles.record_outer_route_feedback!(
        _PERF_OUTER_ROUTE_STATE[],
        features;
        route=route,
        successes=successes,
        failures=failures,
        elapsed_success_s=elapsed_success_s,
        elapsed_success_sq_sum_s=elapsed_success_sq_sum_s,
        tuning=_outer_route_tuning()
    )
    return nothing
end

@inline function _priority_outer_route(case::BenchmarkCase)::Symbol
    features = _outer_route_features(case)
    return ParallelProfiles.default_outer_route(
        features;
        tuning=_outer_route_tuning(),
        machine_class=_machine_parallel_class(),
        threads_available=perf_threads_available(),
        parallel_enabled=perf_parallel_enabled()
    )
end

@inline function parallel_priority_plan(case::BenchmarkCase, outer_route::Symbol)::ParallelPriorityPlan
    resolved_outer = outer_route == :auto ? _priority_outer_route(case) : outer_route
    n_sats = length(case.args_template.dynamics_model.spacecraft)
    n_links = _total_links(case)
    dynamic_effectors = case.args_template.dynamics_model.dynamic_effectors
    has_atmo = _has_atmo_dynamic_effector(dynamic_effectors)
    has_nbody = _has_nbody_dynamic_effector(dynamic_effectors)
    has_density = has_atmo || !(case.args_template.environment_model.density_model isa NoAtmosphereModel)
    has_control = !isempty(case.args_template.control_model.control_effectors)
    sat_threshold = _priority_inner_sat_threshold()
    link_threshold = _priority_inner_link_threshold()

    if resolved_outer != :none
        return ParallelPriorityPlan(
            outer_route=resolved_outer,
            density_mode="off",
            control_mode="off",
            multibody_mode="off"
        )
    end

    density_mode = has_density && n_sats >= sat_threshold ? "auto" : "off"
    control_mode = has_control && n_sats >= sat_threshold ? "auto" : "off"
    multibody_mode = (has_atmo || has_nbody) && n_links >= link_threshold ? "auto" : "off"

    return ParallelPriorityPlan(
        outer_route=resolved_outer,
        density_mode=density_mode,
        control_mode=control_mode,
        multibody_mode=multibody_mode
    )
end

@inline function _existing_or_policy_env(name::String, value::String)::String
    # Only respect user-provided baseline overrides captured at startup.
    # This avoids inheriting transient values left by other benchmark phases.
    baseline = get(_PERF_POLICY_ENV_BASELINE, name, nothing)
    if baseline === nothing
        return value
    end
    token = lowercase(strip(baseline))
    valid = if name == "SPACEAGORA_OUTER_PARALLEL_ACTIVE"
        token in ("1", "0", "true", "false", "yes", "no", "on", "off")
    else
        token in ("off", "none", "serial", "0", "false", "no", "on", "thread", "threads", "1", "true", "yes", "auto")
    end
    return valid ? baseline : value
end

@inline function parallel_priority_env_pairs(plan::ParallelPriorityPlan)::Vector{Pair{String, String}}
    outer_flag = plan.outer_route == :none ? "0" : "1"
    serial_inner_off = plan.density_mode == "off" && plan.control_mode == "off" && plan.multibody_mode == "off"
    rhs_batch_mode = (plan.outer_route == :none && serial_inner_off) ? "off" : "auto"
    return Pair{String, String}[
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _existing_or_policy_env("SPACEAGORA_OUTER_PARALLEL_ACTIVE", outer_flag),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => _existing_or_policy_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL", plan.density_mode),
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => _existing_or_policy_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL", plan.control_mode),
        "SPACEAGORA_MULTIBODY_PARALLEL" => _existing_or_policy_env("SPACEAGORA_MULTIBODY_PARALLEL", plan.multibody_mode),
        "SPACEAGORA_RHS_BATCH_PARALLEL" => _existing_or_policy_env("SPACEAGORA_RHS_BATCH_PARALLEL", rhs_batch_mode),
    ]
end

@inline function _merged_env_pairs(
    base_pairs::Vector{Pair{String, String}},
    override_pairs::Vector{Pair{String, String}}
)::Vector{Pair{String, String}}
    isempty(override_pairs) && return base_pairs
    merged = copy(base_pairs)
    for pair in override_pairs
        key = first(pair)
        idx = findfirst(p -> first(p) == key, merged)
        if idx === nothing
            push!(merged, key => last(pair))
        else
            merged[idx] = key => last(pair)
        end
    end
    return merged
end

@inline function case_env_pairs(
    case::BenchmarkCase,
    plan::ParallelPriorityPlan
)::Vector{Pair{String, String}}
    return _merged_env_pairs(parallel_priority_env_pairs(plan), case.env_overrides)
end

@inline function _env_pairs_key(env_pairs::Vector{Pair{String, String}})::String
    isempty(env_pairs) && return ""
    return join([string(first(p), "=", last(p)) for p in env_pairs], ";")
end

function _thread_plan_groups(
    cases::Vector{BenchmarkCase},
    indices::Vector{Int},
    outer_route::Symbol
)::Vector{Tuple{Vector{Pair{String, String}}, Vector{Tuple{Int, ParallelPriorityPlan}}}}
    grouped_pairs = Dict{String, Vector{Pair{String, String}}}()
    grouped_payload = Dict{String, Vector{Tuple{Int, ParallelPriorityPlan}}}()
    ordered_keys = String[]
    for idx in indices
        plan = parallel_priority_plan(cases[idx], outer_route)
        env_pairs = case_env_pairs(cases[idx], plan)
        key = _env_pairs_key(env_pairs)
        if !haskey(grouped_payload, key)
            grouped_pairs[key] = env_pairs
            grouped_payload[key] = Tuple{Int, ParallelPriorityPlan}[]
            push!(ordered_keys, key)
        end
        push!(grouped_payload[key], (idx, plan))
    end
    return [(grouped_pairs[key], grouped_payload[key]) for key in ordered_keys]
end

function auto_backend_for_case(
    case::BenchmarkCase;
    spec::Union{Nothing, ProfileSpec}=nothing
)::Symbol
    features = _outer_route_features(case; spec=spec)
    return ParallelProfiles.select_outer_route!(
        _PERF_OUTER_ROUTE_STATE[],
        features;
        tuning=_outer_route_tuning(),
        machine_class=_machine_parallel_class(),
        threads_available=perf_threads_available(),
        parallel_enabled=perf_parallel_enabled()
    )
end

@inline function perf_process_workers_target()::Int
    raw = strip(get(ENV, "SPACEAGORA_PERF_PROCS", "0"))
    if isempty(raw) || raw == "0"
        return max(1, Sys.CPU_THREADS - 1)
    end
    n = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_PERF_PROCS must be an integer, got '$raw'"))
    end
    return max(1, n)
end

@inline function _is_valid_worker_project(path::String)::Bool
    return isdir(path) && (
        isfile(joinpath(path, "Project.toml")) ||
        isfile(joinpath(path, "JuliaProject.toml"))
    )
end

@inline function _project_candidate_path(raw::AbstractString)::String
    token = strip(String(raw))
    isempty(token) && return ""
    return normpath(isabspath(token) ? token : joinpath(REPO_ROOT, token))
end

function _resolve_perf_worker_project_path()::String
    checked = String[]
    override_raw = strip(get(ENV, "SPACEAGORA_PERF_WORKER_PROJECT", ""))
    if !isempty(override_raw)
        override_path = _project_candidate_path(override_raw)
        push!(checked, override_path)
        if _is_valid_worker_project(override_path)
            return override_path
        end
    end

    for candidate in (joinpath(REPO_ROOT, ".AGORA"), REPO_ROOT)
        candidate_path = normpath(candidate)
        candidate_path in checked || push!(checked, candidate_path)
        if _is_valid_worker_project(candidate_path)
            return candidate_path
        end
    end

    throw(ArgumentError(
        "Unable to resolve process-worker Julia project. " *
        "Set SPACEAGORA_PERF_WORKER_PROJECT to a valid project directory " *
        "(must contain Project.toml or JuliaProject.toml). " *
        "Checked: $(join(checked, ", "))."
    ))
end

const _perf_workers_initialized = Ref(false)
const _perf_worker_planet_cache = Ref{Any}(nothing)
const _perf_worker_mars_cache = Ref{Any}(nothing)
const _perf_harmonics20_cache = Ref{Any}(nothing)
const _perf_harmonics50_cache = Ref{Any}(nothing)
const _perf_mars_gram_point_density_cache = Ref{Any}(nothing)

function ensure_perf_workers!()
    _perf_workers_initialized[] && return nothing
    target_workers = perf_process_workers_target()
    missing_workers = target_workers - nworkers()
    if missing_workers > 0
        worker_project = _resolve_perf_worker_project_path()
        addprocs(
            missing_workers;
            exeflags=`--startup-file=no --project=$(worker_project)`
        )
    end
    script_path = abspath(joinpath(dirname(@__DIR__), "..", "performance_runtime_analysis.jl"))
    @everywhere workers() include($script_path)
    init_tasks = map(workers()) do w
        @async begin
            remotecall_wait(perf_worker_planet, w)
            remotecall_wait(perf_worker_mars, w)
            remotecall_wait(perf_worker_prime_gram_bindings!, w)
        end
    end
    foreach(wait, init_tasks)
    _perf_workers_initialized[] = true
    return nothing
end

function perf_worker_planet()::Earth
    cached = _perf_worker_planet_cache[]
    if cached === nothing
        cached = Earth("", SPICE_PATH)
        _perf_worker_planet_cache[] = cached
    end
    return cached::Earth
end

function perf_worker_mars()::Mars
    cached = _perf_worker_mars_cache[]
    if cached === nothing
        cached = Mars("", SPICE_PATH)
        _perf_worker_mars_cache[] = cached
    end
    return cached::Mars
end

function perf_worker_prime_gram_bindings!()::Nothing
    # Julia 1.12 is stricter about global bindings during cross-process
    # deserialization; eagerly load both GRAM wrappers on each worker.
    _build_earth_gram_point_density()
    _build_mars_gram_point_density()
    return nothing
end

function _perf_harmonics20_model(planet::Earth)
    cached = _perf_harmonics20_cache[]
    if cached === nothing
        cached = GravitationalHarmonicsModel(20, 20, EARTH_HARMONICS_FILE, planet)
        _perf_harmonics20_cache[] = cached
    end
    return deepcopy(cached)
end

function _perf_harmonics50_model(planet::Earth)
    cached = _perf_harmonics50_cache[]
    if cached === nothing
        cached = GravitationalHarmonicsModel(50, 50, EARTH_HARMONICS_FILE, planet)
        _perf_harmonics50_cache[] = cached
    end
    return deepcopy(cached)
end

function _perf_mars_gram_point_density_model()
    cached = _perf_mars_gram_point_density_cache[]
    if cached === nothing
        cached = _build_mars_gram_point_density()
        _perf_mars_gram_point_density_cache[] = cached
    end
    return deepcopy(cached)
end

@inline function orbital_period_seconds(spacecraft::SpacecraftModel, planet::AbstractPlanet)::Float64
    a = spacecraft.initial_condition.a
    if !isfinite(a) || a <= 0.0
        throw(ArgumentError("Invalid semimajor axis for orbital period calculation: $a"))
    end
    return 2π * sqrt(a^3 / planet.μ)
end

