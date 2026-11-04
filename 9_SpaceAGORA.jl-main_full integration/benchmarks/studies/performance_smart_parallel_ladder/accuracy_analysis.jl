function _zero_crossings(
    signal::Vector{Float64},
    times::Vector{Float64};
    direction::Symbol
)::Vector{Float64}
    n = min(length(signal), length(times))
    n >= 2 || return Float64[]
    crossings = Float64[]
    for i in 2:n
        v0 = signal[i - 1]
        v1 = signal[i]
        t0 = times[i - 1]
        t1 = times[i]
        if !(isfinite(v0) && isfinite(v1) && isfinite(t0) && isfinite(t1))
            continue
        end
        if direction == :up
            if v0 < 0.0 && v1 >= 0.0
                push!(crossings, _linear_zero_crossing_time(t0, t1, v0, v1))
            end
        elseif direction == :down
            if v0 > 0.0 && v1 <= 0.0
                push!(crossings, _linear_zero_crossing_time(t0, t1, v0, v1))
            end
        else
            throw(ArgumentError("Unsupported crossing direction '$direction'. Use :up or :down."))
        end
    end
    return crossings
end

function _primary_event_metrics(sol, case::BenchmarkCase)::NamedTuple
    if isempty(sol.t) || isempty(sol.u)
        return (
            periapsis_first_t=missing,
            periapsis_count=0,
            interface_first_t=missing,
            interface_count=0
        )
    end

    times = Float64.(sol.t)
    radial_signal = Float64[]
    interface_signal = Float64[]
    has_interface = hasproperty(case.args_template.environment_model, :EI)
    entry_interface_m = has_interface ? Float64(case.args_template.environment_model.EI) * 1e3 : NaN
    planet_radius_m = Float64(case.args_template.environment_model.planet.Rp_e)

    for state in sol.u
        sc = state.sc[1]
        push!(radial_signal, dot(sc.pos, sc.vel))
        if has_interface
            altitude_m = norm(sc.pos) - planet_radius_m
            push!(interface_signal, altitude_m - entry_interface_m)
        end
    end

    periapsis_times = _zero_crossings(radial_signal, times; direction=:up)
    interface_times = has_interface ? _zero_crossings(interface_signal, times; direction=:down) : Float64[]
    return (
        periapsis_first_t=isempty(periapsis_times) ? missing : periapsis_times[1],
        periapsis_count=length(periapsis_times),
        interface_first_t=isempty(interface_times) ? missing : interface_times[1],
        interface_count=length(interface_times)
    )
end

function _mean_control_isp_s(control_effectors::Tuple)::Union{Missing, Float64}
    values = Float64[]
    for effector in control_effectors
        if !hasproperty(effector, :Isp)
            continue
        end
        isp = getproperty(effector, :Isp)
        if isp isa Real
            v = Float64(isp)
            if isfinite(v) && v > 0.0
                push!(values, v)
            end
        elseif isp isa AbstractArray
            for item in isp
                item isa Real || continue
                v = Float64(item)
                if isfinite(v) && v > 0.0
                    push!(values, v)
                end
            end
        end
    end
    isempty(values) && return missing
    return mean(values)
end

function _primary_control_metrics(sol, case::BenchmarkCase)::NamedTuple
    if isempty(sol.u)
        return (
            initial_mass_kg=missing,
            final_mass_kg=missing,
            propellant_used_kg=missing,
            control_impulse_equivalent_ns=missing
        )
    end
    first_sc = sol.u[1].sc[1]
    last_sc = sol.u[end].sc[1]
    if !(hasproperty(first_sc, :m) && hasproperty(last_sc, :m))
        return (
            initial_mass_kg=missing,
            final_mass_kg=missing,
            propellant_used_kg=missing,
            control_impulse_equivalent_ns=missing
        )
    end
    mass0 = Float64(getproperty(first_sc, :m))
    massf = Float64(getproperty(last_sc, :m))
    if !(isfinite(mass0) && isfinite(massf))
        return (
            initial_mass_kg=missing,
            final_mass_kg=missing,
            propellant_used_kg=missing,
            control_impulse_equivalent_ns=missing
        )
    end
    propellant_used = max(0.0, mass0 - massf)
    mean_isp = _mean_control_isp_s(case.args_template.control_model.control_effectors)
    impulse_equivalent = if mean_isp isa Missing
        missing
    else
        propellant_used * 9.80665 * mean_isp
    end
    return (
        initial_mass_kg=mass0,
        final_mass_kg=massf,
        propellant_used_kg=propellant_used,
        control_impulse_equivalent_ns=impulse_equivalent
    )
end

function _run_ladder_case_solution(
    case::BenchmarkCase,
    rung::LadderRungSpec,
    config::SmartLadderConfig
)::NamedTuple
    args_run = deepcopy(case.args_template)
    env_pairs = _ladder_env_pairs(rung, config)
    resolved_outer_route = if rung.backend == "none"
        :none
    elseif rung.backend == "threads"
        :threads
    elseif rung.backend == "process"
        :process
    else
        auto_backend_for_case(case; spec=config.profile)
    end
    outer_active = resolved_outer_route == :none ? "0" : "1"
    push!(env_pairs, "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => outer_active)
    started_ns = time_ns()
    try
        run_once = () -> _run_perf_simulation(
            args_run;
            return_solution=true,
            return_solver_metadata=false,
            profile_name=config.profile.name,
            solver_mode_override=case.solver_mode_override,
            split_imex_solver_override=case.split_imex_solver_override,
            entry_target_count_override=case.entry_target_count_override
        )
        result = withenv(env_pairs...) do
            if isempty(case.env_overrides)
                run_once()
            else
                withenv(case.env_overrides...) do
                    run_once()
                end
            end
        end
        elapsed_s = (time_ns() - started_ns) / 1e9
        solve_ok = _solve_success_for_case(result.solution, case)
        if !solve_ok
            return (
                ok=false,
                elapsed_s=elapsed_s,
                result=nothing,
                reason="retcode=$(result.solution.retcode)"
            )
        end
        return (
            ok=true,
            elapsed_s=elapsed_s,
            result=result,
            reason=missing
        )
    catch err
        if err isa InterruptException
            rethrow()
        end
        elapsed_s = (time_ns() - started_ns) / 1e9
        return (
            ok=false,
            elapsed_s=elapsed_s,
            result=nothing,
            reason="$(typeof(err)): $(sprint(showerror, err))"
        )
    end
end

@inline function _default_deep_accuracy_case_names(spec::ProfileSpec)::Vector{String}
    base = String[
        "single_j2",
        "single_entry_earth_nominal",
        "thermal_8sat_panel12_aero",
        "multi_sat_control_8sat_thruster",
        "proximity_2sat_orientation_fullstack_gnc_highrate",
        "effector_6sat_dual_srp_stack",
    ]
    if spec.name == "full"
        append!(base, [
            "articulated_1sat_panel28_fullstack",
            "thermal_aerobrake_mars_panel16",
            "long_constellation_12sat",
        ])
    end
    return base
end

function _deep_accuracy_case_catalog(spec::ProfileSpec)::Dict{String, BenchmarkCase}
    planet = perf_worker_planet()
    selected = selected_cases(spec, build_cases(spec, planet))
    return Dict(case.name => case for case in selected)
end

function _configured_deep_accuracy_case_names(
    spec::ProfileSpec,
    catalog::Dict{String, BenchmarkCase}
)::Vector{String}
    raw = strip(get(ENV, "SPACEAGORA_SMART_LADDER_DEEP_ACCURACY_CASES", ""))
    requested = if isempty(raw)
        _default_deep_accuracy_case_names(spec)
    else
        [strip(token) for token in split(raw, ",") if !isempty(strip(token))]
    end
    filtered = String[]
    for name in requested
        if haskey(catalog, name) && !(name in filtered)
            push!(filtered, name)
        end
    end
    return filtered
end

@inline function _ks_distance(reference::Vector{Float64}, candidate::Vector{Float64})
    n_ref = length(reference)
    n_cmp = length(candidate)
    (n_ref > 0 && n_cmp > 0) || return missing
    ref = sort(reference)
    cmp = sort(candidate)
    i = 1
    j = 1
    d = 0.0
    while i <= n_ref || j <= n_cmp
        x = if j > n_cmp || (i <= n_ref && ref[i] <= cmp[j])
            ref[i]
        else
            cmp[j]
        end
        while i <= n_ref && ref[i] <= x
            i += 1
        end
        while j <= n_cmp && cmp[j] <= x
            j += 1
        end
        cdf_ref = (i - 1) / n_ref
        cdf_cmp = (j - 1) / n_cmp
        d = max(d, abs(cdf_ref - cdf_cmp))
    end
    return d
end

@inline function _finite_metric_values(df::DataFrame, column::Symbol)::Vector{Float64}
    (column in names(df)) || return Float64[]
    values = Float64[]
    for value in df[!, column]
        value isa Real || continue
        v = Float64(value)
        isfinite(v) || continue
        push!(values, v)
    end
    return values
end

function _distribution_parity_metrics(
    reference::Vector{Float64},
    candidate::Vector{Float64}
)::NamedTuple
    if isempty(reference) || isempty(candidate)
        return (
            mean_rel_delta=missing,
            std_rel_delta=missing,
            p90_rel_delta=missing,
            ks_distance=missing
        )
    end
    reference_mean = mean(reference)
    candidate_mean = mean(candidate)
    reference_std = std(reference; corrected=false)
    candidate_std = std(candidate; corrected=false)
    reference_p90 = quantile(reference, 0.9)
    candidate_p90 = quantile(candidate, 0.9)
    return (
        mean_rel_delta=_safe_rel_error(reference_mean, candidate_mean),
        std_rel_delta=_safe_rel_error(reference_std, candidate_std),
        p90_rel_delta=_safe_rel_error(reference_p90, candidate_p90),
        ks_distance=_ks_distance(reference, candidate)
    )
end

function _prepare_accuracy_sample_table(raw_df::DataFrame)::DataFrame
    required = (
        :final_primary_pos_norm_m,
        :final_primary_vel_norm_mps,
        :final_primary_mass_kg
    )
    for column in required
        _has_column(raw_df, column) || return DataFrame()
    end

    rows = NamedTuple[]
    for row in eachrow(raw_df)
        if hasproperty(row, :solve_success) && row.solve_success !== true
            continue
        end
        push!(rows, (
            match_key=_match_key(row),
            final_primary_pos_norm_m=getproperty(row, :final_primary_pos_norm_m),
            final_primary_vel_norm_mps=getproperty(row, :final_primary_vel_norm_mps),
            final_primary_mass_kg=getproperty(row, :final_primary_mass_kg)
        ))
    end

    sample_df = DataFrame(rows)
    nrow(sample_df) == 0 && return sample_df
    return combine(
        groupby(sample_df, :match_key),
        :final_primary_pos_norm_m => _mean_skipmissing => :final_primary_pos_norm_m,
        :final_primary_vel_norm_mps => _mean_skipmissing => :final_primary_vel_norm_mps,
        :final_primary_mass_kg => _mean_skipmissing => :final_primary_mass_kg
    )
end

function _build_accuracy_parity_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _baseline_artifact_or_nothing(artifacts)
    baseline === nothing && return DataFrame()
    baseline_samples = _prepare_accuracy_sample_table(baseline.raw_df)
    if nrow(baseline_samples) == 0
        return DataFrame()
    end
    baseline_samples = select(
        baseline_samples,
        :match_key,
        :final_primary_pos_norm_m => :r0_final_primary_pos_norm_m,
        :final_primary_vel_norm_mps => :r0_final_primary_vel_norm_mps,
        :final_primary_mass_kg => :r0_final_primary_mass_kg
    )

    rows = NamedTuple[]
    for artifact in artifacts
        rung_samples = _prepare_accuracy_sample_table(artifact.raw_df)
        nrow(rung_samples) == 0 && continue
        joined = innerjoin(
            baseline_samples,
            select(
                rung_samples,
                :match_key,
                :final_primary_pos_norm_m => :rung_final_primary_pos_norm_m,
                :final_primary_vel_norm_mps => :rung_final_primary_vel_norm_mps,
                :final_primary_mass_kg => :rung_final_primary_mass_kg
            ),
            on=:match_key
        )
        nrow(joined) == 0 && continue

        pos_err = [_safe_rel_error(r0, rt) for (r0, rt) in zip(joined.r0_final_primary_pos_norm_m, joined.rung_final_primary_pos_norm_m)]
        vel_err = [_safe_rel_error(r0, rt) for (r0, rt) in zip(joined.r0_final_primary_vel_norm_mps, joined.rung_final_primary_vel_norm_mps)]
        mass_err = [_safe_rel_error(r0, rt) for (r0, rt) in zip(joined.r0_final_primary_mass_kg, joined.rung_final_primary_mass_kg)]

        pos_stats = _error_stats(pos_err)
        vel_stats = _error_stats(vel_err)
        mass_stats = _error_stats(mass_err)
        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            samples=nrow(joined),
            pos_rel_err_median_pct=ismissing(pos_stats.median) ? missing : 100.0 * pos_stats.median,
            pos_rel_err_p90_pct=ismissing(pos_stats.p90) ? missing : 100.0 * pos_stats.p90,
            pos_rel_err_max_pct=ismissing(pos_stats.max) ? missing : 100.0 * pos_stats.max,
            vel_rel_err_median_pct=ismissing(vel_stats.median) ? missing : 100.0 * vel_stats.median,
            vel_rel_err_p90_pct=ismissing(vel_stats.p90) ? missing : 100.0 * vel_stats.p90,
            vel_rel_err_max_pct=ismissing(vel_stats.max) ? missing : 100.0 * vel_stats.max,
            mass_rel_err_median_pct=ismissing(mass_stats.median) ? missing : 100.0 * mass_stats.median,
            mass_rel_err_p90_pct=ismissing(mass_stats.p90) ? missing : 100.0 * mass_stats.p90,
            mass_rel_err_max_pct=ismissing(mass_stats.max) ? missing : 100.0 * mass_stats.max
        ))
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, :mode)
    end
    return df
end

function _build_deep_accuracy_parity_table(
    artifacts::Vector{ModeRunArtifacts},
    rungs::Vector{LadderRungSpec},
    rung_label_by_mode::Dict{Symbol, String},
    config::SmartLadderConfig
)::DataFrame
    if !_parse_bool_token(get(ENV, "SPACEAGORA_SMART_LADDER_DEEP_ACCURACY", "1"))
        return DataFrame()
    end

    sample_count = try
        parsed = parse(Int, strip(get(ENV, "SPACEAGORA_SMART_LADDER_DEEP_ACCURACY_SAMPLES", "64")))
        max(8, parsed)
    catch
        64
    end

    catalog = _deep_accuracy_case_catalog(config.profile)
    isempty(catalog) && return DataFrame()
    case_names = _configured_deep_accuracy_case_names(config.profile, catalog)
    isempty(case_names) && return DataFrame()

    rung_by_mode = Dict(r.mode => r for r in rungs)
    baseline_rung = get(rung_by_mode, :serial, nothing)
    baseline_rung === nothing && return DataFrame()

    baseline_runs = Dict{String, NamedTuple}()
    for case_name in case_names
        baseline_runs[case_name] = _run_ladder_case_solution(catalog[case_name], baseline_rung, config)
    end

    rows = NamedTuple[]
    for artifact in artifacts
        rung = get(rung_by_mode, artifact.mode, nothing)
        rung === nothing && continue

        pos_rms_rel = Union{Missing, Float64}[]
        pos_max_rel = Union{Missing, Float64}[]
        vel_rms_rel = Union{Missing, Float64}[]
        vel_max_rel = Union{Missing, Float64}[]
        periapsis_time_abs_err_s = Union{Missing, Float64}[]
        interface_time_abs_err_s = Union{Missing, Float64}[]
        propellant_rel_err = Union{Missing, Float64}[]
        impulse_rel_err = Union{Missing, Float64}[]
        periapsis_count_abs_err = Union{Missing, Float64}[]
        interface_count_abs_err = Union{Missing, Float64}[]
        callback_exact_matches = 0
        callback_comparisons = 0
        compared_cases = 0
        failed_cases = 0

        for case_name in case_names
            baseline_run = get(baseline_runs, case_name, nothing)
            if baseline_run === nothing || baseline_run.ok !== true
                failed_cases += 1
                continue
            end
            case = catalog[case_name]
            rung_run = artifact.mode == :serial ? baseline_run : _run_ladder_case_solution(case, rung, config)
            if rung_run.ok !== true
                failed_cases += 1
                continue
            end

            sol_r0 = baseline_run.result.solution
            sol_rt = rung_run.result.solution
            n_sats = length(case.args_template.dynamics_model.spacecraft)
            orientation = case.args_template.mission_configuration.orientation_sim
            traj = _trajectory_delta_metrics(sol_r0, sol_rt, n_sats, orientation; n_samples=sample_count)
            push!(pos_rms_rel, traj.pos_rel_rms)
            push!(pos_max_rel, traj.pos_rel_max)
            push!(vel_rms_rel, traj.vel_rel_rms)
            push!(vel_max_rel, traj.vel_rel_max)

            events_r0 = _primary_event_metrics(sol_r0, case)
            events_rt = _primary_event_metrics(sol_rt, case)
            peri_count_delta = abs(events_rt.periapsis_count - events_r0.periapsis_count)
            interface_count_delta = abs(events_rt.interface_count - events_r0.interface_count)
            push!(periapsis_time_abs_err_s, _safe_abs_error(events_r0.periapsis_first_t, events_rt.periapsis_first_t))
            push!(interface_time_abs_err_s, _safe_abs_error(events_r0.interface_first_t, events_rt.interface_first_t))
            push!(periapsis_count_abs_err, Float64(peri_count_delta))
            push!(interface_count_abs_err, Float64(interface_count_delta))
            callback_comparisons += 1
            if peri_count_delta == 0 && interface_count_delta == 0
                callback_exact_matches += 1
            end

            control_r0 = _primary_control_metrics(sol_r0, case)
            control_rt = _primary_control_metrics(sol_rt, case)
            push!(propellant_rel_err, _safe_rel_error(control_r0.propellant_used_kg, control_rt.propellant_used_kg))
            push!(impulse_rel_err, _safe_rel_error(control_r0.control_impulse_equivalent_ns, control_rt.control_impulse_equivalent_ns))
            compared_cases += 1
        end

        pos_rms_stats = _error_stats_finite(pos_rms_rel)
        pos_max_stats = _error_stats_finite(pos_max_rel)
        vel_rms_stats = _error_stats_finite(vel_rms_rel)
        vel_max_stats = _error_stats_finite(vel_max_rel)
        peri_time_stats = _error_stats_finite(periapsis_time_abs_err_s)
        interface_time_stats = _error_stats_finite(interface_time_abs_err_s)
        propellant_stats = _error_stats_finite(propellant_rel_err)
        impulse_stats = _error_stats_finite(impulse_rel_err)
        peri_count_stats = _error_stats_finite(periapsis_count_abs_err)
        interface_count_stats = _error_stats_finite(interface_count_abs_err)

        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            probe_cases_requested=length(case_names),
            probe_cases_compared=compared_cases,
            probe_cases_failed=failed_cases,
            trajectory_samples_per_case=sample_count,
            traj_pos_rel_rms_median_pct=ismissing(pos_rms_stats.median) ? missing : 100.0 * pos_rms_stats.median,
            traj_pos_rel_rms_p90_pct=ismissing(pos_rms_stats.p90) ? missing : 100.0 * pos_rms_stats.p90,
            traj_pos_rel_max_pct=ismissing(pos_max_stats.max) ? missing : 100.0 * pos_max_stats.max,
            traj_vel_rel_rms_median_pct=ismissing(vel_rms_stats.median) ? missing : 100.0 * vel_rms_stats.median,
            traj_vel_rel_rms_p90_pct=ismissing(vel_rms_stats.p90) ? missing : 100.0 * vel_rms_stats.p90,
            traj_vel_rel_max_pct=ismissing(vel_max_stats.max) ? missing : 100.0 * vel_max_stats.max,
            periapsis_time_abs_err_median_s=peri_time_stats.median,
            periapsis_time_abs_err_p90_s=peri_time_stats.p90,
            periapsis_time_abs_err_max_s=peri_time_stats.max,
            interface_time_abs_err_median_s=interface_time_stats.median,
            interface_time_abs_err_p90_s=interface_time_stats.p90,
            interface_time_abs_err_max_s=interface_time_stats.max,
            propellant_rel_err_median_pct=ismissing(propellant_stats.median) ? missing : 100.0 * propellant_stats.median,
            propellant_rel_err_p90_pct=ismissing(propellant_stats.p90) ? missing : 100.0 * propellant_stats.p90,
            propellant_rel_err_max_pct=ismissing(propellant_stats.max) ? missing : 100.0 * propellant_stats.max,
            control_impulse_rel_err_median_pct=ismissing(impulse_stats.median) ? missing : 100.0 * impulse_stats.median,
            control_impulse_rel_err_p90_pct=ismissing(impulse_stats.p90) ? missing : 100.0 * impulse_stats.p90,
            control_impulse_rel_err_max_pct=ismissing(impulse_stats.max) ? missing : 100.0 * impulse_stats.max,
            periapsis_count_abs_err_median=peri_count_stats.median,
            periapsis_count_abs_err_max=peri_count_stats.max,
            interface_count_abs_err_median=interface_count_stats.median,
            interface_count_abs_err_max=interface_count_stats.max,
            callback_exact_match_pct=callback_comparisons > 0 ? 100.0 * callback_exact_matches / callback_comparisons : missing
        ))
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, :mode)
    end
    return df
end

function _build_montecarlo_distribution_parity_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _baseline_artifact_or_nothing(artifacts)
    baseline === nothing && return DataFrame()
    baseline_raw = baseline.raw_df
    if !(_has_column(baseline_raw, :category) && _has_column(baseline_raw, :scenario) && _has_column(baseline_raw, :solve_success))
        return DataFrame()
    end
    baseline_success = [value === true for value in baseline_raw.solve_success]
    baseline_mc = baseline_raw[
        (String.(baseline_raw.category) .== "montecarlo") .& baseline_success,
        :
    ]
    nrow(baseline_mc) == 0 && return DataFrame()

    scenarios = sort(unique(String.(baseline_mc.scenario)))
    rows = NamedTuple[]

    for artifact in artifacts
        raw_df = artifact.raw_df
        if !(_has_column(raw_df, :category) && _has_column(raw_df, :scenario) && _has_column(raw_df, :solve_success))
            continue
        end
        rung_success = [value === true for value in raw_df.solve_success]
        rung_mc = raw_df[
            (String.(raw_df.category) .== "montecarlo") .& rung_success,
            :
        ]
        nrow(rung_mc) == 0 && continue

        for scenario_name in scenarios
            baseline_scenario = baseline_mc[String.(baseline_mc.scenario) .== scenario_name, :]
            rung_scenario = rung_mc[String.(rung_mc.scenario) .== scenario_name, :]
            if nrow(baseline_scenario) == 0 || nrow(rung_scenario) == 0
                continue
            end

            pos_metrics = _distribution_parity_metrics(
                _finite_metric_values(baseline_scenario, :final_primary_pos_norm_m),
                _finite_metric_values(rung_scenario, :final_primary_pos_norm_m)
            )
            vel_metrics = _distribution_parity_metrics(
                _finite_metric_values(baseline_scenario, :final_primary_vel_norm_mps),
                _finite_metric_values(rung_scenario, :final_primary_vel_norm_mps)
            )
            mass_metrics = _distribution_parity_metrics(
                _finite_metric_values(baseline_scenario, :final_primary_mass_kg),
                _finite_metric_values(rung_scenario, :final_primary_mass_kg)
            )
            event_metrics = _distribution_parity_metrics(
                _finite_metric_values(baseline_scenario, :terminal_time_s),
                _finite_metric_values(rung_scenario, :terminal_time_s)
            )

            rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
            push!(rows, (
                rung=rung_label,
                mode=String(artifact.mode),
                scenario=scenario_name,
                baseline_samples=nrow(baseline_scenario),
                rung_samples=nrow(rung_scenario),
                pos_mean_rel_delta_pct=ismissing(pos_metrics.mean_rel_delta) ? missing : 100.0 * pos_metrics.mean_rel_delta,
                pos_std_rel_delta_pct=ismissing(pos_metrics.std_rel_delta) ? missing : 100.0 * pos_metrics.std_rel_delta,
                pos_p90_rel_delta_pct=ismissing(pos_metrics.p90_rel_delta) ? missing : 100.0 * pos_metrics.p90_rel_delta,
                pos_ks_distance=pos_metrics.ks_distance,
                vel_mean_rel_delta_pct=ismissing(vel_metrics.mean_rel_delta) ? missing : 100.0 * vel_metrics.mean_rel_delta,
                vel_std_rel_delta_pct=ismissing(vel_metrics.std_rel_delta) ? missing : 100.0 * vel_metrics.std_rel_delta,
                vel_p90_rel_delta_pct=ismissing(vel_metrics.p90_rel_delta) ? missing : 100.0 * vel_metrics.p90_rel_delta,
                vel_ks_distance=vel_metrics.ks_distance,
                mass_mean_rel_delta_pct=ismissing(mass_metrics.mean_rel_delta) ? missing : 100.0 * mass_metrics.mean_rel_delta,
                mass_std_rel_delta_pct=ismissing(mass_metrics.std_rel_delta) ? missing : 100.0 * mass_metrics.std_rel_delta,
                mass_p90_rel_delta_pct=ismissing(mass_metrics.p90_rel_delta) ? missing : 100.0 * mass_metrics.p90_rel_delta,
                mass_ks_distance=mass_metrics.ks_distance,
                event_time_mean_rel_delta_pct=ismissing(event_metrics.mean_rel_delta) ? missing : 100.0 * event_metrics.mean_rel_delta,
                event_time_std_rel_delta_pct=ismissing(event_metrics.std_rel_delta) ? missing : 100.0 * event_metrics.std_rel_delta,
                event_time_p90_rel_delta_pct=ismissing(event_metrics.p90_rel_delta) ? missing : 100.0 * event_metrics.p90_rel_delta,
                event_time_ks_distance=event_metrics.ks_distance
            ))
        end
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, [:mode, :scenario])
    end
    return df
end

function _build_route_mix_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    rows = NamedTuple[]
    for artifact in artifacts
        raw_df = artifact.raw_df
        total_rows = nrow(raw_df)
        successful_rows = _has_column(raw_df, :solve_success) ? count(==(true), raw_df.solve_success) : total_rows
        route_vals = _has_column(raw_df, :outer_route) ? lowercase.(string.(raw_df.outer_route)) : fill("unknown", total_rows)
        none_count = count(==("none"), route_vals)
        threads_count = count(==("threads"), route_vals)
        process_count = count(==("process"), route_vals)
        other_count = max(0, total_rows - (none_count + threads_count + process_count))
        denom = total_rows > 0 ? total_rows : 1
        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            total_rows=total_rows,
            successful_rows=successful_rows,
            none_count=none_count,
            threads_count=threads_count,
            process_count=process_count,
            other_count=other_count,
            none_pct=100.0 * none_count / denom,
            threads_pct=100.0 * threads_count / denom,
            process_pct=100.0 * process_count / denom,
            other_pct=100.0 * other_count / denom
        ))
    end
    df = DataFrame(rows)
    sort!(df, :mode)
    return df
end

function _gate_counts(df::Union{Nothing, DataFrame})::NamedTuple
    if df === nothing
        return (rows=0, pass_rows=0, pass_rate_pct=missing)
    end
    rows = nrow(df)
    if rows == 0 || !_has_column(df, :pass_all)
        return (rows=rows, pass_rows=0, pass_rate_pct=missing)
    end
    pass_rows = count(Bool.(df.pass_all))
    return (rows=rows, pass_rows=pass_rows, pass_rate_pct=100.0 * pass_rows / rows)
end

function _build_fidelity_parity_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String},
    config::SmartLadderConfig
)::DataFrame
    rows = NamedTuple[]
    for artifact in artifacts
        raw_df = artifact.raw_df
        total_rows = nrow(raw_df)
        success_rows = _has_column(raw_df, :solve_success) ? count(==(true), raw_df.solve_success) : total_rows
        success_rate = total_rows > 0 ? 100.0 * success_rows / total_rows : missing
        split_counts = _gate_counts(artifact.split_gate_df)

        agg_dir = dirname(artifact.raw_path)
        multirate_csv_path = _latest_artifact_path_optional(
            agg_dir,
            "multirate_rollout_gate_agg",
            config.profile.name,
            ".csv"
        )
        multirate_df = isnothing(multirate_csv_path) ? nothing : CSV.read(multirate_csv_path, DataFrame)
        multirate_counts = _gate_counts(multirate_df)

        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            sample_rows=total_rows,
            success_rows=success_rows,
            success_rate_pct=success_rate,
            split_gate_rows=split_counts.rows,
            split_gate_pass_rows=split_counts.pass_rows,
            split_gate_pass_rate_pct=split_counts.pass_rate_pct,
            multirate_gate_rows=multirate_counts.rows,
            multirate_gate_pass_rows=multirate_counts.pass_rows,
            multirate_gate_pass_rate_pct=multirate_counts.pass_rate_pct
        ))
    end
    df = DataFrame(rows)
    sort!(df, :mode)
    return df
end
