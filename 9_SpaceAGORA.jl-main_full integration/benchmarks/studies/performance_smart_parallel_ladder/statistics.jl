function _build_protocol_summary_df(
    config::SmartLadderConfig,
    rungs::Vector{LadderRungSpec},
    pass_results::Vector{LadderPassResult};
    primary_solver_label::String=""
)::DataFrame
    tags = _protocol_mode_tags(config)
    rows = NamedTuple[]
    solver_labels = unique([result.solver_label for result in pass_results])
    sort!(solver_labels)
    for solver_label in solver_labels
        solver_runs = [result for result in pass_results if result.solver_label == solver_label]
        solver_mode = isempty(solver_runs) ? missing : (solver_runs[1].solver_mode === nothing ? missing : String(solver_runs[1].solver_mode))
        for rung in rungs
            runs = [result for result in solver_runs if result.rung.label == rung.label]
            isempty(runs) && continue
            total_elapsed = Float64[result.artifact.elapsed_s for result in runs]
            bench_elapsed = Float64[result.artifact.bench_elapsed_s for result in runs]
            orbit_elapsed = Float64[result.artifact.orbit_elapsed_s for result in runs]
            split_elapsed = Float64[result.artifact.split_gate_elapsed_s for result in runs]

            total_stats = _pass_stats(total_elapsed)
            bench_stats = _pass_stats(bench_elapsed)
            orbit_stats = _pass_stats(orbit_elapsed)
            split_stats = _pass_stats(split_elapsed)

            push!(rows, (
                profile=config.profile.name,
                solver_label=solver_label,
                solver_mode=solver_mode,
                solver_axis=String(config.solver_axis),
                solver_primary=(solver_label == primary_solver_label),
                rung=rung.label,
                mode=String(rung.mode),
                matrix=String(rung.matrix),
                backend=rung.backend,
                inner_adaptive=rung.inner_adaptive,
                outer_route_adaptive=rung.outer_route_adaptive,
                pass_count=length(runs),
                passes_requested=config.passes,
                randomize_rung_order=config.randomize_rung_order,
                random_seed=config.random_seed,
                start_mode=tags.start_mode,
                cache_mode=tags.cache_mode,
                cache_cold=tags.cache_cold,
                outer_state_reset=tags.outer_state_reset,
                inner_state_reset=tags.inner_state_reset,
                total_elapsed_mean_s=total_stats.mean_s,
                total_elapsed_std_s=total_stats.std_s,
                total_elapsed_cv_pct=total_stats.cv_pct,
                total_elapsed_ci95_low_s=total_stats.ci95_low_s,
                total_elapsed_ci95_high_s=total_stats.ci95_high_s,
                bench_elapsed_mean_s=bench_stats.mean_s,
                bench_elapsed_std_s=bench_stats.std_s,
                bench_elapsed_cv_pct=bench_stats.cv_pct,
                bench_elapsed_ci95_low_s=bench_stats.ci95_low_s,
                bench_elapsed_ci95_high_s=bench_stats.ci95_high_s,
                orbit_elapsed_mean_s=orbit_stats.mean_s,
                orbit_elapsed_std_s=orbit_stats.std_s,
                orbit_elapsed_cv_pct=orbit_stats.cv_pct,
                orbit_elapsed_ci95_low_s=orbit_stats.ci95_low_s,
                orbit_elapsed_ci95_high_s=orbit_stats.ci95_high_s,
                split_gate_elapsed_mean_s=split_stats.mean_s,
                split_gate_elapsed_std_s=split_stats.std_s,
                split_gate_elapsed_cv_pct=split_stats.cv_pct,
                split_gate_elapsed_ci95_low_s=split_stats.ci95_low_s,
                split_gate_elapsed_ci95_high_s=split_stats.ci95_high_s
            ))
        end
    end
    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, [:solver_label, :rung])
    end
    return df
end

@inline function _key_token(value)::String
    if value isa Missing || value === nothing
        return "_"
    elseif value isa AbstractFloat
        v = Float64(value)
        return isfinite(v) ? string(round(v; digits=6)) : "_"
    end
    return string(value)
end

@inline function _mission_family(row)::String
    category = hasproperty(row, :category) ? lowercase(string(getproperty(row, :category))) : ""
    scenario = hasproperty(row, :scenario) ? lowercase(string(getproperty(row, :scenario))) : ""
    mission_time_s = hasproperty(row, :mission_time_s) && getproperty(row, :mission_time_s) isa Real ? Float64(getproperty(row, :mission_time_s)) : 0.0
    satellites = hasproperty(row, :satellites) && getproperty(row, :satellites) isa Integer ? Int(getproperty(row, :satellites)) : 1
    dynamic_effectors = hasproperty(row, :dynamic_effectors) ? lowercase(string(getproperty(row, :dynamic_effectors))) : ""
    control_effectors = hasproperty(row, :control_effectors) ? lowercase(string(getproperty(row, :control_effectors))) : ""

    if category == "montecarlo" || occursin("montecarlo", scenario)
        return mission_time_s <= 7200.0 ? "mc_short" : "mc_long"
    end
    if category == "thermal_entry" || occursin("thermal_aerobrake", scenario)
        return "thermal_entry"
    end
    if category == "srp_heavy" || occursin("srp_heavy", scenario)
        return "srp_heavy"
    end
    if category == "articulated_multibody" || occursin("articulated", scenario)
        return "articulated_multibody"
    end
    if category == "multi_sat_control" || (satellites >= 2 && control_effectors != "" && control_effectors != "none")
        return "multi_sat_control"
    end
    if category == "long_constellation" || occursin("long_constellation", scenario)
        return "long_constellation"
    end
    if category == "effector_stress" || occursin("effector", scenario)
        return "effector_stress"
    end
    if satellites >= 2
        return "multi_sat"
    end
    if occursin("nbody", dynamic_effectors) || occursin("harmonic", dynamic_effectors)
        return "high_fidelity_nbody"
    end
    if mission_time_s > 7200.0 || category == "mission_length"
        return "long"
    end
    return "short_light"
end

@inline function _is_thermal_scenario(row)::Bool
    category = hasproperty(row, :category) ? lowercase(string(getproperty(row, :category))) : ""
    scenario = hasproperty(row, :scenario) ? lowercase(string(getproperty(row, :scenario))) : ""
    return category in ("thermal_stress", "thermal_entry") || occursin("thermal", scenario)
end

@inline function _match_key(row)::String
    parts = (
        _key_token(hasproperty(row, :pass) ? getproperty(row, :pass) : missing),
        _key_token(hasproperty(row, :scenario) ? getproperty(row, :scenario) : missing),
        _key_token(hasproperty(row, :category) ? getproperty(row, :category) : missing),
        _key_token(hasproperty(row, :repeat) ? getproperty(row, :repeat) : missing),
        _key_token(hasproperty(row, :seed) ? getproperty(row, :seed) : missing),
        _key_token(hasproperty(row, :mission_time_s) ? getproperty(row, :mission_time_s) : missing),
        _key_token(hasproperty(row, :satellites) ? getproperty(row, :satellites) : missing)
    )
    return join(parts, "|")
end

function _baseline_artifact(artifacts::Vector{ModeRunArtifacts})::ModeRunArtifacts
    idx = findfirst(a -> a.mode == :serial, artifacts)
    idx === nothing && error("Smart ladder requires a serial baseline rung with mode=:serial.")
    return artifacts[idx]
end

function _baseline_artifact_or_nothing(
    artifacts::Vector{ModeRunArtifacts}
)::Union{Nothing, ModeRunArtifacts}
    idx = findfirst(a -> a.mode == :serial, artifacts)
    idx === nothing && return nothing
    return artifacts[idx]
end

function _prepare_speed_sample_table(raw_df::DataFrame)::DataFrame
    rows = NamedTuple[]
    for row in eachrow(raw_df)
        if hasproperty(row, :solve_success) && row.solve_success !== true
            continue
        end
        total_time = hasproperty(row, :total_time_s) ? getproperty(row, :total_time_s) : missing
        if !(total_time isa Real) || !isfinite(Float64(total_time)) || Float64(total_time) <= 0.0
            continue
        end
        push!(rows, (
            match_key=_match_key(row),
            total_time_s=Float64(total_time),
            mission_family=_mission_family(row)
        ))
    end
    sample_df = DataFrame(rows)
    nrow(sample_df) == 0 && return sample_df
    return combine(
        groupby(sample_df, :match_key),
        :total_time_s => mean => :total_time_s,
        :mission_family => (v -> first(v)) => :mission_family
    )
end

@inline function _is_layer_attribution_mode(mode::Symbol)::Bool
    return startswith(String(mode), "attr_")
end

@inline function _layer_set_from_mode(mode::Symbol)::String
    token = replace(String(mode), "attr_" => "")
    token == "outer_only" && return "outer_only"
    return replace(token, "_" => "+")
end

@inline function _layer_kind_from_mode(mode::Symbol)::String
    token = replace(String(mode), "attr_" => "")
    token == "outer_only" && return "baseline"
    return occursin("_", token) ? "pairwise" : "single_layer"
end

@inline function _artifact_by_mode(
    artifacts::Vector{ModeRunArtifacts},
    mode::Symbol
)::Union{Nothing, ModeRunArtifacts}
    idx = findfirst(a -> a.mode == mode, artifacts)
    return idx === nothing ? nothing : artifacts[idx]
end

function _build_layer_attribution_speedup_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _artifact_by_mode(artifacts, :attr_outer_only)
    baseline === nothing && return DataFrame()

    rows = NamedTuple[]
    for artifact in artifacts
        if !_is_layer_attribution_mode(artifact.mode)
            continue
        end
        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            layer_kind=_layer_kind_from_mode(artifact.mode),
            layer_set=_layer_set_from_mode(artifact.mode),
            total_elapsed_s=artifact.elapsed_s,
            run_benchmarks_elapsed_s=artifact.bench_elapsed_s,
            run_per_orbit_elapsed_s=artifact.orbit_elapsed_s,
            total_speedup_vs_outer_only=_safe_ratio(baseline.elapsed_s, artifact.elapsed_s),
            run_benchmarks_speedup_vs_outer_only=_safe_ratio(baseline.bench_elapsed_s, artifact.bench_elapsed_s),
            run_per_orbit_speedup_vs_outer_only=_safe_ratio(baseline.orbit_elapsed_s, artifact.orbit_elapsed_s)
        ))
    end
    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, [:layer_kind, :mode])
    end
    return df
end

function _build_layer_attribution_mission_family_speedup_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _artifact_by_mode(artifacts, :attr_outer_only)
    baseline === nothing && return DataFrame()

    baseline_samples = _prepare_speed_sample_table(baseline.raw_df)
    nrow(baseline_samples) == 0 && return DataFrame()
    baseline_samples = select(baseline_samples, :match_key, :total_time_s => :outer_only_total_time_s)

    rows = NamedTuple[]
    for artifact in artifacts
        if !_is_layer_attribution_mode(artifact.mode) || artifact.mode == :attr_outer_only
            continue
        end
        rung_samples = _prepare_speed_sample_table(artifact.raw_df)
        nrow(rung_samples) == 0 && continue
        joined = innerjoin(
            baseline_samples,
            select(rung_samples, :match_key, :mission_family, :total_time_s => :rung_total_time_s),
            on=:match_key
        )
        nrow(joined) == 0 && continue
        joined[!, :speedup_vs_outer_only] = [
            Float64(base) / Float64(rt)
            for (base, rt) in zip(joined.outer_only_total_time_s, joined.rung_total_time_s)
        ]
        joined = joined[isfinite.(joined.speedup_vs_outer_only) .& (joined.speedup_vs_outer_only .> 0.0), :]
        nrow(joined) == 0 && continue

        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        for grp in groupby(joined, :mission_family)
            vals = Float64.(grp.speedup_vs_outer_only)
            sample_count = length(vals)
            slower_count = count(v -> v < 1.0, vals)
            push!(rows, (
                rung=rung_label,
                mode=String(artifact.mode),
                layer_kind=_layer_kind_from_mode(artifact.mode),
                layer_set=_layer_set_from_mode(artifact.mode),
                mission_family=String(first(grp.mission_family)),
                samples=sample_count,
                median_speedup_vs_outer_only=median(vals),
                p90_speedup_vs_outer_only=quantile(vals, 0.9),
                worst_slowdown_x=maximum(1.0 ./ vals),
                slower_share_pct=100.0 * slower_count / sample_count
            ))
        end
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, [:layer_kind, :layer_set, :mission_family])
    end
    return df
end

function _build_layer_pairwise_synergy_table(
    artifacts::Vector{ModeRunArtifacts}
)::DataFrame
    baseline = _artifact_by_mode(artifacts, :attr_outer_only)
    baseline === nothing && return DataFrame()

    baseline_samples = _prepare_speed_sample_table(baseline.raw_df)
    nrow(baseline_samples) == 0 && return DataFrame()
    baseline_samples = select(
        baseline_samples,
        :match_key,
        :mission_family,
        :total_time_s => :outer_total
    )

    pair_defs = (
        (pair=:attr_density_thermal, a=:attr_density, b=:attr_thermal),
        (pair=:attr_density_multibody, a=:attr_density, b=:attr_multibody),
        (pair=:attr_control_effector, a=:attr_control, b=:attr_effector),
        (pair=:attr_multibody_effector, a=:attr_multibody, b=:attr_effector),
    )

    rows = NamedTuple[]
    for def in pair_defs
        pair_art = _artifact_by_mode(artifacts, def.pair)
        a_art = _artifact_by_mode(artifacts, def.a)
        b_art = _artifact_by_mode(artifacts, def.b)
        if pair_art === nothing || a_art === nothing || b_art === nothing
            continue
        end

        a_samples = select(_prepare_speed_sample_table(a_art.raw_df), :match_key, :total_time_s => :a_total)
        b_samples = select(_prepare_speed_sample_table(b_art.raw_df), :match_key, :total_time_s => :b_total)
        pair_samples = select(_prepare_speed_sample_table(pair_art.raw_df), :match_key, :total_time_s => :pair_total)
        if nrow(a_samples) == 0 || nrow(b_samples) == 0 || nrow(pair_samples) == 0
            continue
        end

        joined = innerjoin(baseline_samples, a_samples, on=:match_key)
        joined = innerjoin(joined, b_samples, on=:match_key)
        joined = innerjoin(joined, pair_samples, on=:match_key)
        nrow(joined) == 0 && continue

        speedup_a = Float64[]
        speedup_b = Float64[]
        speedup_pair = Float64[]
        predicted_independent = Float64[]
        synergy_ratio = Float64[]
        for row in eachrow(joined)
            s_a = Float64(row.outer_total) / Float64(row.a_total)
            s_b = Float64(row.outer_total) / Float64(row.b_total)
            s_pair = Float64(row.outer_total) / Float64(row.pair_total)
            pred = s_a * s_b
            if !(isfinite(s_a) && isfinite(s_b) && isfinite(s_pair) && isfinite(pred))
                continue
            end
            if s_a <= 0.0 || s_b <= 0.0 || s_pair <= 0.0 || pred <= 0.0
                continue
            end
            push!(speedup_a, s_a)
            push!(speedup_b, s_b)
            push!(speedup_pair, s_pair)
            push!(predicted_independent, pred)
            push!(synergy_ratio, s_pair / pred)
        end
        isempty(synergy_ratio) && continue

        sample_count = length(synergy_ratio)
        slower_share = 100.0 * count(v -> v < 1.0, speedup_pair) / sample_count
        synergy_positive_share = 100.0 * count(v -> v > 1.0, synergy_ratio) / sample_count
        push!(rows, (
            pair_mode=String(def.pair),
            pair_set=_layer_set_from_mode(def.pair),
            samples=sample_count,
            median_pair_speedup_vs_outer_only=median(speedup_pair),
            p90_pair_speedup_vs_outer_only=quantile(speedup_pair, 0.9),
            median_predicted_independent_speedup=median(predicted_independent),
            p90_predicted_independent_speedup=quantile(predicted_independent, 0.9),
            median_synergy_ratio=median(synergy_ratio),
            p90_synergy_ratio=quantile(synergy_ratio, 0.9),
            min_synergy_ratio=minimum(synergy_ratio),
            max_synergy_ratio=maximum(synergy_ratio),
            synergy_positive_share_pct=synergy_positive_share,
            slower_than_outer_only_share_pct=slower_share
        ))
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, :pair_mode)
    end
    return df
end

function _build_vs_r0_speedup_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _baseline_artifact_or_nothing(artifacts)
    baseline === nothing && return DataFrame()
    rows = NamedTuple[]
    for artifact in artifacts
        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            backend=artifact.backend,
            total_elapsed_s=artifact.elapsed_s,
            run_benchmarks_elapsed_s=artifact.bench_elapsed_s,
            run_per_orbit_elapsed_s=artifact.orbit_elapsed_s,
            total_speedup_vs_r0=_safe_ratio(baseline.elapsed_s, artifact.elapsed_s),
            run_benchmarks_speedup_vs_r0=_safe_ratio(baseline.bench_elapsed_s, artifact.bench_elapsed_s),
            run_per_orbit_speedup_vs_r0=_safe_ratio(baseline.orbit_elapsed_s, artifact.orbit_elapsed_s)
        ))
    end
    df = DataFrame(rows)
    sort!(df, :mode)
    return df
end

function _build_mission_family_speedup_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _baseline_artifact_or_nothing(artifacts)
    baseline === nothing && return DataFrame()
    baseline_samples = _prepare_speed_sample_table(baseline.raw_df)
    if nrow(baseline_samples) == 0
        return DataFrame()
    end
    baseline_samples = select(baseline_samples, :match_key, :total_time_s => :r0_total_time_s)

    rows = NamedTuple[]
    for artifact in artifacts
        artifact.mode == :serial && continue
        rung_samples = _prepare_speed_sample_table(artifact.raw_df)
        nrow(rung_samples) == 0 && continue
        joined = innerjoin(
            baseline_samples,
            select(rung_samples, :match_key, :mission_family, :total_time_s => :rung_total_time_s),
            on=:match_key
        )
        nrow(joined) == 0 && continue
        joined[!, :speedup_vs_r0] = [Float64(r0) / Float64(rt) for (r0, rt) in zip(joined.r0_total_time_s, joined.rung_total_time_s)]
        joined = joined[isfinite.(joined.speedup_vs_r0) .& (joined.speedup_vs_r0 .> 0.0), :]
        nrow(joined) == 0 && continue

        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        for grp in groupby(joined, :mission_family)
            vals = Float64.(grp.speedup_vs_r0)
            sample_count = length(vals)
            slower_count = count(v -> v < 1.0, vals)
            push!(rows, (
                rung=rung_label,
                mode=String(artifact.mode),
                mission_family=String(first(grp.mission_family)),
                samples=sample_count,
                median_speedup_vs_r0=median(vals),
                p90_speedup_vs_r0=quantile(vals, 0.9),
                worst_slowdown_x=maximum(1.0 ./ vals),
                slower_share_pct=100.0 * slower_count / sample_count
            ))
        end
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, [:rung, :mission_family])
    end
    return df
end

function _prepare_thermal_speed_sample_table(raw_df::DataFrame)::DataFrame
    rows = NamedTuple[]
    for row in eachrow(raw_df)
        if hasproperty(row, :solve_success) && row.solve_success !== true
            continue
        end
        total_time = hasproperty(row, :total_time_s) ? getproperty(row, :total_time_s) : missing
        if !(total_time isa Real) || !isfinite(Float64(total_time)) || Float64(total_time) <= 0.0
            continue
        end
        push!(rows, (
            match_key=_match_key(row),
            total_time_s=Float64(total_time),
            is_thermal=_is_thermal_scenario(row)
        ))
    end
    sample_df = DataFrame(rows)
    nrow(sample_df) == 0 && return sample_df
    return combine(
        groupby(sample_df, :match_key),
        :total_time_s => mean => :total_time_s,
        :is_thermal => (v -> any(Bool.(v))) => :is_thermal
    )
end

@inline function _thermal_runtime_share_pct(raw_df::DataFrame)
    total_time_s = 0.0
    thermal_time_s = 0.0
    for row in eachrow(raw_df)
        if hasproperty(row, :solve_success) && row.solve_success !== true
            continue
        end
        total_time = hasproperty(row, :total_time_s) ? getproperty(row, :total_time_s) : missing
        if !(total_time isa Real) || !isfinite(Float64(total_time)) || Float64(total_time) <= 0.0
            continue
        end
        total_time_s += Float64(total_time)
        if _is_thermal_scenario(row)
            thermal_time_s += Float64(total_time)
        end
    end
    total_time_s > 0.0 || return missing
    return 100.0 * thermal_time_s / total_time_s
end

function _build_thermal_contribution_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _baseline_artifact_or_nothing(artifacts)
    baseline === nothing && return DataFrame()
    baseline_samples = _prepare_thermal_speed_sample_table(baseline.raw_df)
    if nrow(baseline_samples) == 0 || !_has_column(baseline_samples, :is_thermal)
        return DataFrame()
    end
    baseline_samples = baseline_samples[baseline_samples.is_thermal .== true, :]
    nrow(baseline_samples) == 0 && return DataFrame()
    baseline_samples = select(
        baseline_samples,
        :match_key,
        :total_time_s => :r0_total_time_s
    )

    rows = NamedTuple[]
    for artifact in artifacts
        rung_samples = _prepare_thermal_speed_sample_table(artifact.raw_df)
        if nrow(rung_samples) == 0 || !_has_column(rung_samples, :is_thermal)
            continue
        end
        rung_samples = rung_samples[rung_samples.is_thermal .== true, :]
        nrow(rung_samples) == 0 && continue

        joined = innerjoin(
            baseline_samples,
            select(rung_samples, :match_key, :total_time_s => :rung_total_time_s),
            on=:match_key
        )
        nrow(joined) == 0 && continue
        joined[!, :speedup_vs_r0] = [Float64(r0) / Float64(rt) for (r0, rt) in zip(joined.r0_total_time_s, joined.rung_total_time_s)]
        joined = joined[isfinite.(joined.speedup_vs_r0) .& (joined.speedup_vs_r0 .> 0.0), :]
        nrow(joined) == 0 && continue

        vals = Float64.(joined.speedup_vs_r0)
        sample_count = length(vals)
        slower_count = count(v -> v < 1.0, vals)
        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            thermal_samples=sample_count,
            median_speedup_vs_r0=median(vals),
            p90_speedup_vs_r0=quantile(vals, 0.9),
            best_speedup_vs_r0=maximum(vals),
            worst_slowdown_x=maximum(1.0 ./ vals),
            slower_share_pct=100.0 * slower_count / sample_count,
            thermal_runtime_share_pct=_thermal_runtime_share_pct(artifact.raw_df)
        ))
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, :mode)
    end
    return df
end

@inline function _mean_skipmissing(values)
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return mean(Float64.(vec))
end

@inline function _safe_rel_error(r0, value)
    if !(r0 isa Real) || !(value isa Real)
        return missing
    end
    r0_f = Float64(r0)
    value_f = Float64(value)
    if !isfinite(r0_f) || !isfinite(value_f)
        return missing
    end
    denom = max(abs(r0_f), eps(Float64))
    return abs(value_f - r0_f) / denom
end

@inline function _error_stats(values)::NamedTuple
    vec = collect(skipmissing(values))
    isempty(vec) && return (median=missing, p90=missing, max=missing)
    return (
        median=median(vec),
        p90=quantile(vec, 0.9),
        max=maximum(vec)
    )
end

@inline function _error_stats_finite(values)::NamedTuple
    vec = Float64[]
    for value in values
        value isa Missing && continue
        value isa Real || continue
        v = Float64(value)
        isfinite(v) || continue
        push!(vec, v)
    end
    isempty(vec) && return (median=missing, p90=missing, max=missing)
    return (
        median=median(vec),
        p90=quantile(vec, 0.9),
        max=maximum(vec)
    )
end

@inline function _safe_abs_error(reference, value)
    if !(reference isa Real) || !(value isa Real)
        return missing
    end
    r = Float64(reference)
    v = Float64(value)
    if !isfinite(r) || !isfinite(v)
        return missing
    end
    return abs(v - r)
end

@inline function _linear_zero_crossing_time(t0::Float64, t1::Float64, v0::Float64, v1::Float64)::Float64
    dv = v1 - v0
    if !isfinite(dv) || abs(dv) <= eps(Float64)
        return t1
    end
    alpha = clamp((-v0) / dv, 0.0, 1.0)
    return t0 + alpha * (t1 - t0)
end
