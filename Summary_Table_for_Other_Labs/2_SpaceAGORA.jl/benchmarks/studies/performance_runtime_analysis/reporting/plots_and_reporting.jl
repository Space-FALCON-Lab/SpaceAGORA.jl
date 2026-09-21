function _plot_wrapped_label(name::AbstractString; width::Int=16, max_lines::Int=3)::String
    words = split(_plot_label(name))
    isempty(words) && return ""
    lines = String[]
    current = ""
    for word in words
        if isempty(current)
            current = word
        elseif ncodeunits(current) + 1 + ncodeunits(word) <= width
            current *= " " * word
        else
            push!(lines, current)
            current = word
        end
    end
    push!(lines, current)
    if length(lines) > max_lines
        head = lines[1:(max_lines - 1)]
        tail = join(lines[max_lines:end], " ")
        push!(head, tail)
        lines = head
    end
    return join(lines, "\n")
end

@inline function _plot_ready()::Bool
    return myid() == 1 && isdefined(Main, :Plots)
end

const _runtime_plot_theme_applied = Ref(false)

function _ensure_runtime_plot_theme!()
    !_plot_ready() && return nothing
    _runtime_plot_theme_applied[] && return nothing

    Plots.theme(:ggplot2)
    Plots.default(
        dpi=220,
        lw=2,
        ms=5,
        markerstrokewidth=1.3,
        markerstrokecolor=:black,
        titlefont=Plots.font(24),
        guidefont=Plots.font(16),
        tickfont=Plots.font(12),
        legend_font=Plots.font(11),
        legend_background_color=:white,
        legend_foreground_color=:black,
        gridalpha=0.24,
        minorgrid=false,
        framestyle=:box
    )
    _runtime_plot_theme_applied[] = true
    return nothing
end

@inline function _plot_margins(;
    size::Tuple{Int, Int}=(2200, 1100),
    left_mm::Int=18,
    right_mm::Int=18,
    top_mm::Int=8,
    bottom_mm::Int=20,
    legend=false,
    xrotation::Real=0
)
    return (
        size=size,
        left_margin=left_mm * Plots.mm,
        right_margin=right_mm * Plots.mm,
        top_margin=top_mm * Plots.mm,
        bottom_margin=bottom_mm * Plots.mm,
        legend=legend,
        xrotation=xrotation,
        framestyle=:box,
        gridalpha=0.24,
        legend_background_color=:white,
        legend_foreground_color=:black
    )
end

function _plot_metric_pairs(df::DataFrame, label_col::Symbol, metric_col::Symbol; axis_labels::Bool=false)
    labels = String[]
    values = Float64[]
    for row in eachrow(df)
        value = row[metric_col]
        if !(value isa Missing)
            f = Float64(value)
            if isfinite(f)
                raw = string(row[label_col])
                push!(labels, axis_labels ? _plot_axis_label(raw) : _plot_label(raw))
                push!(values, f)
            end
        end
    end
    return labels, values
end

@inline function _has_row_fields(row, fields::Vector{Symbol})::Bool
    for field in fields
        if row[field] isa Missing
            return false
        end
    end
    return true
end

function _save_runtime_plot!(
    artifacts::Vector{String},
    plt,
    outdir::String,
    basename::String,
    spec::ProfileSpec,
    stamp::String
)
    path = joinpath(outdir, "$(basename)_$(spec.name)_$(stamp).png")
    try
        Plots.savefig(plt, path)
        push!(artifacts, path)
    catch err
        @warn "[perf] failed to save plot $basename: $(_perf_error_text(err))"
    end
    return nothing
end

@inline function _paper_figure_pack_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_PAPER_FIGURE_PACK", true)
end

@inline function _paper_ladder_outdir()::String
    raw = strip(get(ENV, "SPACEAGORA_PERF_PAPER_LADDER_OUTDIR", joinpath(DEFAULT_OUTPUT_DIR, "smart_parallel_ladder")))
    return normpath(isabspath(raw) ? raw : joinpath(REPO_ROOT, raw))
end

@inline function _paper_cross_machine_outdir()::String
    raw = strip(get(ENV, "SPACEAGORA_PERF_PAPER_CROSS_MACHINE_OUTDIR", joinpath(DEFAULT_OUTPUT_DIR, "smart_parallel_ladder_cross_machine")))
    return normpath(isabspath(raw) ? raw : joinpath(REPO_ROOT, raw))
end

function _latest_profile_artifact_optional(
    outdir::String,
    prefix::String,
    profile::String;
    suffix::String=".csv"
)::Union{Nothing, String}
    isdir(outdir) || return nothing
    pattern = "$(prefix)_$(profile)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

function _read_optional_dataframe(path::Union{Nothing, String}, label::String)::Union{Nothing, DataFrame}
    path === nothing && return nothing
    try
        return CSV.read(path, DataFrame)
    catch err
        @warn "[perf] Failed to read optional paper-figure artifact; skipping." label=label path=path exception=(err, catch_backtrace())
        return nothing
    end
end

function _paper_figure_external_data(spec::ProfileSpec)
    ladder_outdir = _paper_ladder_outdir()
    cross_outdir = _paper_cross_machine_outdir()

    layer_attr_speedup_path = _latest_profile_artifact_optional(
        ladder_outdir,
        "smart_parallel_ladder_layer_attribution_speedup",
        spec.name
    )
    deep_accuracy_path = _latest_profile_artifact_optional(
        ladder_outdir,
        "smart_parallel_ladder_deep_accuracy_parity",
        spec.name
    )
    montecarlo_parity_path = _latest_profile_artifact_optional(
        ladder_outdir,
        "smart_parallel_ladder_montecarlo_distribution_parity",
        spec.name
    )
    route_mix_path = _latest_profile_artifact_optional(
        ladder_outdir,
        "smart_parallel_ladder_route_mix",
        spec.name
    )
    cross_speedup_summary_path = _latest_profile_artifact_optional(
        cross_outdir,
        "smart_parallel_ladder_cross_machine_speedup_summary",
        spec.name
    )
    cross_adaptive_regret_summary_path = _latest_profile_artifact_optional(
        cross_outdir,
        "smart_parallel_ladder_cross_machine_adaptive_regret_summary",
        spec.name
    )
    cross_route_mix_summary_path = _latest_profile_artifact_optional(
        cross_outdir,
        "smart_parallel_ladder_cross_machine_route_mix_summary",
        spec.name
    )

    return (
        layer_attribution_speedup_df=_read_optional_dataframe(layer_attr_speedup_path, "layer_attribution_speedup"),
        deep_accuracy_df=_read_optional_dataframe(deep_accuracy_path, "deep_accuracy_parity"),
        montecarlo_parity_df=_read_optional_dataframe(montecarlo_parity_path, "montecarlo_distribution_parity"),
        route_mix_df=_read_optional_dataframe(route_mix_path, "route_mix"),
        cross_speedup_summary_df=_read_optional_dataframe(cross_speedup_summary_path, "cross_machine_speedup_summary"),
        cross_adaptive_regret_summary_df=_read_optional_dataframe(cross_adaptive_regret_summary_path, "cross_machine_adaptive_regret_summary"),
        cross_route_mix_summary_df=_read_optional_dataframe(cross_route_mix_summary_path, "cross_machine_route_mix_summary"),
        layer_attribution_speedup_path=layer_attr_speedup_path,
        deep_accuracy_path=deep_accuracy_path,
        montecarlo_parity_path=montecarlo_parity_path,
        route_mix_path=route_mix_path,
        cross_speedup_summary_path=cross_speedup_summary_path,
        cross_adaptive_regret_summary_path=cross_adaptive_regret_summary_path,
        cross_route_mix_summary_path=cross_route_mix_summary_path
    )
end

function _sorted_orbit_groups(df::DataFrame, metric::Symbol)
    multiplier_col = :mission_time_multiplier in names(df) ? :mission_time_multiplier : :orbit_count
    groups = collect(groupby(df, :scenario))
    sort!(groups; by=g -> begin
        local_df = DataFrame(g)
        sort!(local_df, multiplier_col)
        value = local_df[1, metric]
        return value isa Missing ? Inf : -Float64(value)
    end)
    return groups
end

function generate_runtime_plots(
    outdir::String,
    spec::ProfileSpec,
    stamp::String,
    raw_df::DataFrame,
    summary_df::DataFrame,
    orbit_summary_df::DataFrame,
    entry_duration_summary_df::DataFrame=DataFrame();
    split_gate_df::Union{Nothing, DataFrame}=nothing,
    multirate_gate_df::Union{Nothing, DataFrame}=nothing,
    inner_hint_layer_df::Union{Nothing, DataFrame}=nothing,
    density_backend_breakdown_df::Union{Nothing, DataFrame}=nothing
)::Vector{String}
    !_plot_ready() && return String[]
    _ensure_runtime_plot_theme!()

    plot_artifacts = String[]
    success_summary = summary_df[summary_df.samples_success .> 0, :]
    paper_external = _paper_figure_pack_enabled() ? _paper_figure_external_data(spec) : nothing

    # 1) Mean and p90 runtime per scenario.
    totals_df = success_summary[.!ismissing.(success_summary.total_time_mean_s), :]
    if nrow(totals_df) > 0
        labels = _plot_axis_label.(String.(totals_df.scenario))
        means = Float64.(totals_df.total_time_mean_s)
        p90_vals = [_plot_number(v) for v in totals_df.total_time_p90_s]
        valid_p90 = findall(isfinite, p90_vals)

        plt = Plots.bar(
            labels,
            means;
            label="Mean total time",
            color="#5b8fb9",
            title="Scenario Runtime: Mean + P90 (bars sorted by mean)",
            xlabel="Scenario",
            ylabel="Wall Time [s]",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        if !isempty(valid_p90)
            Plots.scatter!(
                plt,
                labels[valid_p90],
                p90_vals[valid_p90];
                marker=:diamond,
                color=:black,
                label="P90 total time"
            )
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_totals", spec, stamp)
    end

    # 2) Relative speedup vs selected baseline scenario.
    speedup_df = success_summary[.!ismissing.(success_summary.speedup_vs_baseline), :]
    if nrow(speedup_df) > 0
        labels = _plot_axis_label.(String.(speedup_df.scenario))
        speedups = Float64.(speedup_df.speedup_vs_baseline)
        plt = Plots.bar(
            labels,
            speedups;
            label="Speedup",
            color="#4f9d69",
            title="Relative Speedup (higher is faster)",
            xlabel="Scenario",
            ylabel="Speedup vs $(PERF_BASELINE_SCENARIO) [x]",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        Plots.hline!(plt, [1.0]; color=:black, linestyle=:dash, label="Baseline = 1x")
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_speedup", spec, stamp)
    end

    # 3) Robustness-adjusted performance (all attempts) + failure rate.
    if (:penalized_expected_wall_time_s in Symbol.(names(summary_df))) &&
       (:failure_rate in Symbol.(names(summary_df)))
        robust_df = summary_df[.!ismissing.(summary_df.penalized_expected_wall_time_s), :]
        if nrow(robust_df) > 0
            sort!(robust_df, :penalized_expected_wall_time_s)
            labels = _plot_axis_label.(String.(robust_df.scenario))
            penalized = Float64.(robust_df.penalized_expected_wall_time_s)
            failure_pct = [
                v isa Missing ? 0.0 : (100.0 * Float64(v))
                for v in robust_df.failure_rate
            ]
            p1 = Plots.bar(
                labels,
                penalized;
                color="#3f7fb3",
                label=false,
                title="Robustness-Adjusted Performance (All Attempts)",
                xlabel="Scenario",
                ylabel="Penalized Expected Wall Time [s / requested run]",
                _plot_margins(size=(2500, 760), bottom_mm=92, right_mm=32)...
            )
            p2 = Plots.bar(
                labels,
                failure_pct;
                color="#bc4749",
                label=false,
                title="Failure Rate by Scenario (All Attempts)",
                xlabel="Scenario",
                ylabel="Failure Rate [%]",
                _plot_margins(size=(2500, 760), bottom_mm=92, right_mm=32)...
            )
            plt = Plots.plot(
                p1,
                p2;
                layout=(2, 1),
                size=(2500, 1600),
                left_margin=20 * Plots.mm,
                right_margin=18 * Plots.mm
            )
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_robustness_adjusted", spec, stamp)
        end
    end

    # 4) Runtime variability by scenario.
    variability_df = success_summary[
        .!ismissing.(success_summary.total_time_min_s) .&
        .!ismissing.(success_summary.total_time_mean_s) .&
        .!ismissing.(success_summary.total_time_p90_s) .&
        .!ismissing.(success_summary.total_time_max_s), :
    ]
    if nrow(variability_df) > 0
        labels = _plot_axis_label.(String.(variability_df.scenario))
        x = collect(1:length(labels))
        mins = Float64.(variability_df.total_time_min_s)
        means = Float64.(variability_df.total_time_mean_s)
        p90s = Float64.(variability_df.total_time_p90_s)
        maxs = Float64.(variability_df.total_time_max_s)
        plt = Plots.plot(
            x,
            means;
            label="Mean",
            color="#3d83b8",
            marker=:circle,
            xticks=(x, labels),
            title="Runtime Variability by Scenario (Min/Mean/P90/Max)",
            xlabel="Scenario",
            ylabel="Wall Time [s]",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        for i in eachindex(x)
            Plots.plot!(
                plt,
                [x[i], x[i]],
                [mins[i], maxs[i]];
                color=:gray40,
                linewidth=2,
                label=i == 1 ? "Min-Max range" : ""
            )
        end
        Plots.scatter!(plt, x, p90s; marker=:diamond, color=:black, label="P90")
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_variability", spec, stamp)
    end

    # 5) Configuration copy + solve breakdown.
    labels_breakdown = String[]
    copy_vals = Float64[]
    solve_vals = Float64[]
    for row in eachrow(success_summary)
        if _has_row_fields(row, [:copy_time_mean_s, :solve_time_mean_s])
            push!(labels_breakdown, _plot_axis_label(String(row.scenario)))
            push!(copy_vals, Float64(row.copy_time_mean_s))
            push!(solve_vals, Float64(row.solve_time_mean_s))
        end
    end
    if !isempty(labels_breakdown)
        plt = Plots.bar(
            labels_breakdown,
            hcat(copy_vals, solve_vals);
            label=["Copy time" "Solve time"],
            color=["#9fb3c8" "#2a9d8f"],
            title="Runtime Breakdown (Configuration Copy + Solve)",
            xlabel="Scenario",
            ylabel="Wall Time [s]",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_breakdown_copy_solve", spec, stamp)
    end

    # 6) Allocation footprint by scenario (memory + call count).
    alloc_df = success_summary[
        .!ismissing.(success_summary.total_bytes_mean_mb) .&
        .!ismissing.(success_summary.solve_alloc_mean), :
    ]
    if nrow(alloc_df) > 0
        labels = _plot_axis_label.(String.(alloc_df.scenario))
        x = collect(1:length(labels))
        bytes_mb = Float64.(alloc_df.total_bytes_mean_mb)
        alloc_million = Float64.(alloc_df.solve_alloc_mean) ./ 1e6
        p1 = Plots.bar(
            x,
            bytes_mb;
            color="#d17a4f",
            label=false,
            xticks=(x, fill("", length(x))),
            ylabel="Allocated Memory [MB]",
            title="Allocation Footprint by Scenario",
            _plot_margins(size=(2500, 780), bottom_mm=10)...
        )
        p2 = Plots.bar(
            x,
            alloc_million;
            color="#6999a1",
            label=false,
            xticks=(x, labels),
            xlabel="Scenario",
            ylabel="Allocation Calls [million]",
            _plot_margins(size=(2500, 980), bottom_mm=88)...
        )
        plt = Plots.plot(
            p1,
            p2;
            layout=(2, 1),
            size=(2500, 1760),
            left_margin=20 * Plots.mm,
            right_margin=20 * Plots.mm
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_memory_alloc", spec, stamp)
    end

    # 7) Integrator workload and rejection pressure.
    solver_df = success_summary[
        .!ismissing.(success_summary.accepted_steps_mean) .&
        .!ismissing.(success_summary.saved_points_mean), :
    ]
    if nrow(solver_df) > 0
        labels = _plot_axis_label.(String.(solver_df.scenario))
        x = collect(1:length(labels))
        accepted = Float64.(solver_df.accepted_steps_mean)
        saved = Float64.(solver_df.saved_points_mean)
        rejected = [v isa Missing ? 0.0 : Float64(v) for v in solver_df.rejected_steps_mean]
        rejection_ratio = [acc + rej <= 0.0 ? 0.0 : 100.0 * rej / (acc + rej) for (acc, rej) in zip(accepted, rejected)]
        p1 = Plots.plot(
            x,
            accepted;
            label="Accepted steps",
            color="#2e7d32",
            marker=:circle,
            xticks=(x, fill("", length(x))),
            ylabel="Steps / Saved Points",
            title="Integrator Workload per Scenario",
            _plot_margins(size=(2500, 780), bottom_mm=10, legend=:outertopright)...
        )
        Plots.plot!(p1, x, saved; color="#375fd2", marker=:diamond, label="Saved points")
        p2 = Plots.bar(
            x,
            rejection_ratio;
            color="#cc6666",
            label=false,
            xticks=(x, labels),
            title="Solver Rejection Pressure",
            xlabel="Scenario",
            ylabel="Rejected Step Ratio [%]",
            _plot_margins(size=(2500, 980), bottom_mm=88)...
        )
        plt = Plots.plot(
            p1,
            p2;
            layout=(2, 1),
            size=(2500, 1760),
            left_margin=20 * Plots.mm,
            right_margin=20 * Plots.mm
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_solver_workload", spec, stamp)
    end

    # 8) Throughput ranking with per-satellite markers.
    throughput_df = success_summary[.!ismissing.(success_summary.sim_seconds_per_wall_second_mean), :]
    if nrow(throughput_df) > 0
        sort!(throughput_df, :sim_seconds_per_wall_second_mean, rev=true)
        labels = _plot_axis_label.(String.(throughput_df.scenario))
        x = collect(1:length(labels))
        global_tp = Float64.(throughput_df.sim_seconds_per_wall_second_mean)
        sat_tp = [v isa Missing ? NaN : Float64(v) for v in throughput_df.satellite_sim_seconds_per_wall_second_mean]
        plt = Plots.bar(
            x,
            global_tp;
            color="#4f9d69",
            label="Global throughput",
            xticks=(x, labels),
            title="Simulation Throughput Ranking",
            xlabel="Scenario",
            ylabel="Sim Seconds / Wall Second",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        valid_sat = findall(isfinite, sat_tp)
        if !isempty(valid_sat)
            Plots.scatter!(
                plt,
                x[valid_sat],
                sat_tp[valid_sat];
                marker=:utriangle,
                color=:black,
                label="Per-satellite throughput"
            )
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_throughput", spec, stamp)
    end

    # 9) Satellite-count scaling (measured vs ideal linear).
    sat_df = success_summary[
        (success_summary.category .== "satellite_scaling") .&
        .!ismissing.(success_summary.satellites) .&
        .!ismissing.(success_summary.total_time_mean_s), :
    ]
    if nrow(sat_df) > 0
        sort!(sat_df, :satellites)
        sat_count = Int.(sat_df.satellites)
        measured = Float64.(sat_df.total_time_mean_s)
        base_runtime = measured[1]
        ideal = base_runtime .* (sat_count ./ sat_count[1])
        plt = Plots.plot(
            sat_count,
            measured;
            marker=:circle,
            color="#3f7fb3",
            label="Measured runtime",
            title="Satellite-Count Scaling (Measured vs Ideal Linear)",
            xlabel="Number of Satellites",
            ylabel="Mean Runtime [s]",
            _plot_margins(size=(2200, 1200), bottom_mm=20, right_mm=35, legend=:topleft)...
        )
        Plots.plot!(
            plt,
            sat_count,
            ideal;
            marker=:diamond,
            color=:black,
            linestyle=:dash,
            label="Ideal linear from smallest satellite-count case"
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_satellite_scaling", spec, stamp)
    end

    # 10) Dynamics fidelity ladder (absolute + relative).
    fidelity_order = [
        "single_j2",
        "single_nbody_sun_moon",
        "single_harmonics_l20",
        "single_harmonics_l50",
    ]
    fidelity_labels = String[]
    fidelity_values = Float64[]
    for scenario in fidelity_order
        idx = findfirst(==(scenario), success_summary.scenario)
        if idx !== nothing
            value = success_summary.total_time_mean_s[idx]
            if !(value isa Missing)
                if scenario == "single_j2"
                    push!(fidelity_labels, "J2")
                elseif scenario == "single_nbody_sun_moon"
                    push!(fidelity_labels, "NBody\nSun+Moon")
                elseif scenario == "single_harmonics_l20"
                    push!(fidelity_labels, "Harmonics\nL20")
                elseif scenario == "single_harmonics_l50"
                    push!(fidelity_labels, "Harmonics\nL50")
                else
                    push!(fidelity_labels, _plot_label(scenario))
                end
                push!(fidelity_values, Float64(value))
            end
        end
    end
    if length(fidelity_values) >= 2
        x = collect(1:length(fidelity_values))
        baseline = fidelity_values[1]
        relative = baseline <= 0.0 ? fill(NaN, length(fidelity_values)) : fidelity_values ./ baseline
        p_abs = Plots.bar(
            x,
            fidelity_values;
            color="#7668c7",
            label=false,
            xticks=(x, fill("", length(x))),
            ylabel="Mean Runtime [s]",
            title="Dynamics Fidelity Ladder (Absolute Runtime)",
            _plot_margins(size=(2200, 760), bottom_mm=10)...
        )
        p_rel = Plots.bar(
            x,
            relative;
            color="#d67c1c",
            label=false,
            xticks=(x, fidelity_labels),
            xlabel="Dynamics Model",
            ylabel="Runtime / Baseline [x]",
            title="Dynamics Fidelity Ladder (Relative Runtime)",
            _plot_margins(size=(2200, 940), bottom_mm=56)...
        )
        Plots.hline!(p_rel, [1.0]; color=:black, linestyle=:dash, label=false)
        plt = Plots.plot(
            p_abs,
            p_rel;
            layout=(2, 1),
            size=(2200, 1700),
            left_margin=20 * Plots.mm,
            right_margin=18 * Plots.mm
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_fidelity_ladder", spec, stamp)
    end

    # 11) Monte Carlo runtime distribution with mean and p90 lines.
    mc_df = raw_df[(raw_df.category .== "montecarlo") .& (raw_df.solve_success .== true), :]
    mc_times = [Float64(v) for v in mc_df.total_time_s if !(v isa Missing)]
    if !isempty(mc_times)
        mc_mean = mean(mc_times)
        mc_p90 = quantile(mc_times, 0.9)
        plt = Plots.histogram(
            mc_times;
            bins=min(20, max(5, round(Int, sqrt(length(mc_times))))),
            color="#7b63c6",
            alpha=0.65,
            label="Samples",
            title="Monte Carlo Runtime Distribution",
            xlabel="Total Wall Time [s]",
            ylabel="Count",
            _plot_margins(size=(2200, 1300), right_mm=52, legend=:outertopright)...
        )
        Plots.vline!(plt, [mc_mean]; color=:red, linewidth=2, label="Mean")
        Plots.vline!(plt, [mc_p90]; color=:black, linewidth=2, linestyle=:dash, label="P90")
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_montecarlo_hist", spec, stamp)
    end

    # 12) Monte Carlo seed trace with mean and p90 lines.
    mc_seed_df = mc_df[.!ismissing.(mc_df.seed), :]
    if nrow(mc_seed_df) > 0
        sort!(mc_seed_df, :seed)
        seeds = Int.(mc_seed_df.seed)
        totals = Float64.(mc_seed_df.total_time_s)
        seed_mean = mean(totals)
        seed_p90 = quantile(totals, 0.9)
        plt = Plots.plot(
            seeds,
            totals;
            marker=:circle,
            color="#356a97",
            label="Seed runtime",
            title="Monte Carlo Runtime by Seed",
            xlabel="Monte Carlo Seed",
            ylabel="Total Wall Time [s]",
            _plot_margins(size=(2200, 1300), right_mm=52, legend=:outertopright)...
        )
        Plots.hline!(plt, [seed_mean]; color=:red, linewidth=2, label="Mean")
        Plots.hline!(plt, [seed_p90]; color=:black, linewidth=2, linestyle=:dash, label="P90")
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_montecarlo_seed_trace", spec, stamp)
    end

    orbit_required_cols = (:samples_success, :total_time_mean_s, :orbits_per_wall_second_mean, :time_per_orbit_mean_s)
    orbit_valid = if nrow(orbit_summary_df) > 0 && all(col -> col in names(orbit_summary_df), orbit_required_cols)
        sweep_multiplier_col = :mission_time_multiplier in names(orbit_summary_df) ? :mission_time_multiplier : :orbit_count
        orbit_summary_df[
            (orbit_summary_df.samples_success .> 0) .&
            .!ismissing.(orbit_summary_df[!, sweep_multiplier_col]), :
        ]
    else
        DataFrame()
    end

    # 13) Mission-time sweep runtime scaling.
    orbit_scaling_df = nrow(orbit_valid) > 0 ? orbit_valid[.!ismissing.(orbit_valid.total_time_mean_s), :] : DataFrame()
    if nrow(orbit_scaling_df) > 0
        multiplier_col = :mission_time_multiplier in names(orbit_scaling_df) ? :mission_time_multiplier : :orbit_count
        time_per_unit_col = :time_per_baseline_period_mean_s in names(orbit_scaling_df) ? :time_per_baseline_period_mean_s : :time_per_orbit_mean_s
        plt = Plots.plot(;
            title="Mission-Time Sweep Runtime Scaling",
            xlabel="Mission-Time Multiplier [x baseline period]",
            ylabel="Mean Time per Baseline-Period Unit [s]",
            _plot_margins(size=(2300, 1200), right_mm=72, legend=:outerright)...
        )
        for grp in _sorted_orbit_groups(orbit_scaling_df, :total_time_mean_s)
            local_df = DataFrame(grp)
            sort!(local_df, multiplier_col)
            x = Int.(local_df[!, multiplier_col])
            y = Float64.(local_df[!, time_per_unit_col])
            Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_scaling", spec, stamp)
    end

    # 14) Mission-time sweep efficiency scaling.
    orbit_eff_df = nrow(orbit_valid) > 0 ? orbit_valid[.!ismissing.(orbit_valid.orbits_per_wall_second_mean), :] : DataFrame()
    if nrow(orbit_eff_df) > 0
        multiplier_col = :mission_time_multiplier in names(orbit_eff_df) ? :mission_time_multiplier : :orbit_count
        throughput_col = :baseline_periods_per_wall_second_mean in names(orbit_eff_df) ? :baseline_periods_per_wall_second_mean : :orbits_per_wall_second_mean
        plt = Plots.plot(;
            title="Mission-Time Sweep Efficiency Scaling",
            xlabel="Mission-Time Multiplier [x baseline period]",
            ylabel="Baseline-Period Units / Wall-sec",
            _plot_margins(size=(2300, 1200), right_mm=72, legend=:outerright)...
        )
        for grp in _sorted_orbit_groups(orbit_eff_df, throughput_col)
            local_df = DataFrame(grp)
            sort!(local_df, multiplier_col)
            x = Int.(local_df[!, multiplier_col])
            y = Float64.(local_df[!, throughput_col])
            Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_efficiency", spec, stamp)
    end

    # 15) Mission-time sweep time heatmap.
    heat_df = nrow(orbit_valid) > 0 ? orbit_valid[.!ismissing.(orbit_valid.time_per_orbit_mean_s), :] : DataFrame()
    if nrow(heat_df) > 0
        multiplier_col = :mission_time_multiplier in names(heat_df) ? :mission_time_multiplier : :orbit_count
        heat_value_col = :time_per_baseline_period_mean_s in names(heat_df) ? :time_per_baseline_period_mean_s : :time_per_orbit_mean_s
        scenario_names = unique(String.(heat_df.scenario))
        scenario_order = sort(scenario_names; by=sc -> begin
            vals = [
                Float64(v)
                for v in heat_df[heat_df.scenario .== sc, heat_value_col]
                if !(v isa Missing)
            ]
            return isempty(vals) ? Inf : -mean(vals)
        end)
        multipliers = sort(unique(Int.(heat_df[!, multiplier_col])))
        z = fill(NaN, length(scenario_order), length(multipliers))
        for row in eachrow(heat_df)
            si = findfirst(==(String(row.scenario)), scenario_order)
            oi = findfirst(==(Int(row[ multiplier_col ])), multipliers)
            if si !== nothing && oi !== nothing
                z[si, oi] = Float64(row[heat_value_col])
            end
        end
        plt = Plots.heatmap(
            multipliers,
            1:length(scenario_order),
            z;
            color=Plots.cgrad([:lightsteelblue1, :mediumpurple3, "#3a0a2a"]),
            colorbar_title="s / baseline-period unit",
            colorbar=true,
            yticks=(1:length(scenario_order), _plot_wrapped_label.(scenario_order)),
            xlabel="Mission-Time Multiplier [x baseline period]",
            ylabel="Scenario",
            title="Mission-Time Sweep Heatmap [s / baseline-period unit]",
            _plot_margins(size=(2300, 1400), left_mm=72, right_mm=52, bottom_mm=24)...
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_heatmap", spec, stamp)
    end

    # 16) Entry-duration sweep trends (separate from per-orbit mission-time sweep).
    if nrow(entry_duration_summary_df) > 0 &&
       (:entry_run_role in names(entry_duration_summary_df)) &&
       (:entry_atmospheric_interface_count in names(entry_duration_summary_df))
        entry_measured = entry_duration_summary_df[
            (entry_duration_summary_df.entry_run_role .== "measured") .&
            .!ismissing.(entry_duration_summary_df.entry_atmospheric_interface_count), :
        ]
        if nrow(entry_measured) > 0
            if :passage_duration_mean_s in names(entry_measured)
                passage_df = entry_measured[.!ismissing.(entry_measured.passage_duration_mean_s), :]
                if nrow(passage_df) > 0
                    plt = Plots.plot(;
                        title="Entry-Duration Sweep: Passage Duration",
                        xlabel="Atmospheric-Interface Count",
                        ylabel="Passage Duration [s]",
                        _plot_margins(size=(2200, 1200), right_mm=72, legend=:outerright)...
                    )
                    for grp in groupby(passage_df, :scenario)
                        local_df = DataFrame(grp)
                        sort!(local_df, :entry_atmospheric_interface_count)
                        x = Int.(local_df.entry_atmospheric_interface_count)
                        y = Float64.(local_df.passage_duration_mean_s)
                        Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
                    end
                    _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_entry_duration_passage", spec, stamp)
                end
            end

            if :wall_time_per_passage_mean_s in names(entry_measured)
                wall_df = entry_measured[.!ismissing.(entry_measured.wall_time_per_passage_mean_s), :]
                if nrow(wall_df) > 0
                    plt = Plots.plot(;
                        title="Entry-Duration Sweep: Wall Time per Entry Passage",
                        xlabel="Atmospheric-Interface Count",
                        ylabel="Wall Time per Passage [s]",
                        _plot_margins(size=(2200, 1200), right_mm=72, legend=:outerright)...
                    )
                    for grp in groupby(wall_df, :scenario)
                        local_df = DataFrame(grp)
                        sort!(local_df, :entry_atmospheric_interface_count)
                        x = Int.(local_df.entry_atmospheric_interface_count)
                        y = Float64.(local_df.wall_time_per_passage_mean_s)
                        Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
                    end
                    _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_entry_duration_wall_per_passage", spec, stamp)
                end
            end

            if :event_time_abs_error_mean_s in names(entry_measured)
                error_df = entry_measured[.!ismissing.(entry_measured.event_time_abs_error_mean_s), :]
                if nrow(error_df) > 0
                    plt = Plots.plot(;
                        title="Entry-Duration Sweep: Event-Time Absolute Error vs Serial Reference",
                        xlabel="Atmospheric-Interface Count",
                        ylabel="|Event-Time Error| [s]",
                        _plot_margins(size=(2200, 1200), right_mm=72, legend=:outerright)...
                    )
                    for grp in groupby(error_df, :scenario)
                        local_df = DataFrame(grp)
                        sort!(local_df, :entry_atmospheric_interface_count)
                        x = Int.(local_df.entry_atmospheric_interface_count)
                        y = Float64.(local_df.event_time_abs_error_mean_s)
                        Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
                    end
                    _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_entry_duration_event_time_error", spec, stamp)
                end
            end
        end
    end

    # 16) SPICE call budget by scenario (N-body + SRP + planet frame).
    spice_df = success_summary[
        .!ismissing.(success_summary.spice_calls_total_mean) .&
        .!ismissing.(success_summary.total_time_mean_s), :
    ]
    if nrow(spice_df) > 0
        labels = _plot_axis_label.(String.(spice_df.scenario))
        nbody_calls = [v isa Missing ? 0.0 : Float64(v) for v in spice_df.nbody_spkpos_total_calls_mean]
        srp_calls = [v isa Missing ? 0.0 : Float64(v) for v in spice_df.srp_spkpos_total_calls_mean]
        pxform_calls = [v isa Missing ? 0.0 : Float64(v) for v in spice_df.planet_pxform_total_calls_mean]
        total_calls = nbody_calls .+ srp_calls .+ pxform_calls
        calls_per_wall_s = [
            (tt <= 0.0) ? NaN : tc / tt
            for (tc, tt) in zip(total_calls, Float64.(spice_df.total_time_mean_s))
        ]
        if any(>(0.0), total_calls)
            p1 = Plots.bar(
                labels,
                hcat(nbody_calls, srp_calls, pxform_calls);
                label=["N-body spkpos" "SRP spkpos" "Planet pxform"],
                color=["#2f80ed" "#56b870" "#d17a4f"],
                bar_position=:stack,
                ylabel="Mean SPICE Calls / Successful Run",
                title="SPICE Call Budget by Scenario",
                _plot_margins(size=(2500, 760), bottom_mm=10, right_mm=62, legend=:outertopright)...
            )
            p2 = Plots.bar(
                labels,
                calls_per_wall_s;
                color="#7b63c6",
                label=false,
                xlabel="Scenario",
                ylabel="Calls / Wall-sec",
                title="SPICE Call Rate by Scenario",
                _plot_margins(size=(2500, 940), bottom_mm=88, right_mm=42)...
            )
            plt = Plots.plot(
                p1,
                p2;
                layout=(2, 1),
                size=(2500, 1700),
                left_margin=20 * Plots.mm,
                right_margin=20 * Plots.mm
            )
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_spice_budget", spec, stamp)
        end
    end

    # 17) Paper pack: adaptivity behavior (route mix + confidence/regret + best-fixed regret).
    if _paper_figure_pack_enabled()
        panels = Any[]

        route_mix_df = (paper_external === nothing) ? nothing : paper_external.route_mix_df
        if !(route_mix_df === nothing) &&
           nrow(route_mix_df) > 0 &&
           (:none_pct in names(route_mix_df)) &&
           (:threads_pct in names(route_mix_df)) &&
           (:process_pct in names(route_mix_df))
            local_df = copy(route_mix_df)
            if :mode in names(local_df)
                mask = [_is_adaptive_mode_token(v) for v in local_df.mode]
                if any(mask)
                    local_df = local_df[mask, :]
                end
            end
            if nrow(local_df) > 0
                labels = if :rung in names(local_df)
                    _plot_axis_label.(String.(local_df.rung))
                elseif :mode in names(local_df)
                    _plot_axis_label.(String.(local_df.mode))
                else
                    _plot_axis_label.(string.("adaptive_", collect(1:nrow(local_df))))
                end
                none_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.none_pct]
                threads_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.threads_pct]
                process_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.process_pct]
                p_route = Plots.bar(
                    labels,
                    hcat(none_pct, threads_pct, process_pct);
                    label=["none" "threads" "process"],
                    bar_position=:stack,
                    color=["#9aa4b2" "#3f7fb3" "#2f8f5b"],
                    title="Adaptive Route-Choice Distribution",
                    xlabel="Adaptive Mode / Rung",
                    ylabel="Route Share [%]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_route)
            end
        elseif nrow(raw_df) > 0 && (:outer_route in names(raw_df))
            counts = Dict("none" => 0, "threads" => 0, "process" => 0, "other" => 0)
            for v in raw_df.outer_route
                token = lowercase(strip(String(v)))
                if haskey(counts, token)
                    counts[token] += 1
                else
                    counts["other"] += 1
                end
            end
            denom = max(1, nrow(raw_df))
            labels = ["none", "threads", "process", "other"]
            vals = [100.0 * counts[label] / denom for label in labels]
            p_route = Plots.bar(
                labels,
                vals;
                color=["#9aa4b2", "#3f7fb3", "#2f8f5b", "#c06c84"],
                title="Observed Outer-Route Distribution",
                xlabel="Outer Route",
                ylabel="Share of Runs [%]",
                _plot_margins(size=(2400, 760), bottom_mm=42, right_mm=42, legend=false)...
            )
            push!(panels, p_route)
        end

        if !(inner_hint_layer_df === nothing) &&
           nrow(inner_hint_layer_df) > 0 &&
           (:layer in names(inner_hint_layer_df)) &&
           (:confidence_mean in names(inner_hint_layer_df)) &&
           (:regret_mean_ns in names(inner_hint_layer_df))
            hint_df = inner_hint_layer_df[
                .!ismissing.(inner_hint_layer_df.confidence_mean) .&
                .!ismissing.(inner_hint_layer_df.regret_mean_ns), :
            ]
            if nrow(hint_df) > 0
                labels = _plot_axis_label.(String.(hint_df.layer))
                confidence = Float64.(hint_df.confidence_mean)
                regret_ms = Float64.(hint_df.regret_mean_ns) ./ 1e6
                p_hint = Plots.plot(
                    labels,
                    confidence;
                    marker=:circle,
                    color="#355070",
                    label="Mean confidence width",
                    title="Inner Adaptive Hint Confidence/Regret by Layer",
                    xlabel="Inner Layer",
                    ylabel="Confidence (lower is tighter)",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                Plots.plot!(
                    p_hint,
                    labels,
                    regret_ms;
                    marker=:diamond,
                    color="#bc4749",
                    label="Mean regret [ms]"
                )
                push!(panels, p_hint)
            end
        end

        regret_summary_df = (paper_external === nothing) ? nothing : paper_external.cross_adaptive_regret_summary_df
        if !(regret_summary_df === nothing) &&
           nrow(regret_summary_df) > 0 &&
           (:adaptive_mode in names(regret_summary_df)) &&
           (:mean_time_regret_pct in names(regret_summary_df)) &&
           (:win_rate_pct in names(regret_summary_df))
            local_df = regret_summary_df
            labels = _plot_axis_label.(String.(local_df.adaptive_mode))
            mean_regret = [v isa Missing ? NaN : Float64(v) for v in local_df.mean_time_regret_pct]
            win_rate = [v isa Missing ? NaN : Float64(v) for v in local_df.win_rate_pct]
            p_regret = Plots.plot(
                labels,
                mean_regret;
                marker=:circle,
                color="#d17a4f",
                label="Mean time regret vs best fixed [%]",
                title="Adaptive Regret vs Best Fixed (Cross-Machine)",
                xlabel="Adaptive Mode",
                ylabel="Regret / Win Rate [%]",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
            )
            Plots.plot!(
                p_regret,
                labels,
                win_rate;
                marker=:utriangle,
                color="#2a9d8f",
                label="Win rate vs best fixed [%]"
            )
            push!(panels, p_regret)
        end

        if !isempty(panels)
            plt = if length(panels) == 1
                panels[1]
            else
                Plots.plot(
                    panels...;
                    layout=(length(panels), 1),
                    size=(2500, 660 * length(panels)),
                    left_margin=20 * Plots.mm,
                    right_margin=20 * Plots.mm
                )
            end
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_adaptivity_behavior", spec, stamp)
        end
    end

    # 18) Paper pack: explicit layer-attribution speedup panel.
    if _paper_figure_pack_enabled() && !(paper_external === nothing)
        layer_df = paper_external.layer_attribution_speedup_df
        if !(layer_df === nothing) &&
           nrow(layer_df) > 0 &&
           (:layer_set in names(layer_df)) &&
           (:total_speedup_vs_outer_only in names(layer_df))
            speed_by_layer = Dict{String, Float64}()
            for row in eachrow(layer_df)
                layer = lowercase(strip(String(row.layer_set)))
                speed = row.total_speedup_vs_outer_only
                if speed isa Missing
                    continue
                end
                speed_f = Float64(speed)
                if !isfinite(speed_f)
                    continue
                end
                if !haskey(speed_by_layer, layer) || speed_f > speed_by_layer[layer]
                    speed_by_layer[layer] = speed_f
                end
            end
            ordered_layers = ["outer_only", "density", "thermal", "control", "multibody", "effector"]
            labels = String[]
            values = Float64[]
            for layer in ordered_layers
                if haskey(speed_by_layer, layer)
                    push!(labels, _plot_axis_label(layer))
                    push!(values, speed_by_layer[layer])
                end
            end
            if !isempty(values)
                plt = Plots.bar(
                    labels,
                    values;
                    color="#4f9d69",
                    title="Layer Attribution: Speedup vs Outer-Only Baseline",
                    xlabel="Layer Stack",
                    ylabel="Speedup vs outer-only [x]",
                    _plot_margins(size=(2300, 980), bottom_mm=60, right_mm=42, legend=false)...
                )
                Plots.hline!(plt, [1.0]; color=:black, linestyle=:dash, label=false)
                _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_layer_attribution", spec, stamp)
            end
        end
    end

    # 19) Paper pack: density-backend benchmark panel.
    if _paper_figure_pack_enabled()
        density_df = density_backend_breakdown_df === nothing ? summarize_density_backend_breakdown(raw_df) : density_backend_breakdown_df
        if nrow(density_df) > 0 &&
           (:density_backend_bucket in names(density_df)) &&
           (:total_time_mean_s in names(density_df)) &&
           (:sim_seconds_per_wall_second_mean in names(density_df))
            bucket_order = [
                "gram_point_to_point",
                "gram_surrogate",
                "gram_static_grid_or_cached_surrogate",
                "non_gram"
            ]
            labels = String[]
            mean_time = Float64[]
            throughput = Float64[]
            success_rate = Float64[]
            for bucket in bucket_order
                idx = findfirst(==(bucket), String.(density_df.density_backend_bucket))
                idx === nothing && continue
                row = density_df[idx, :]
                push!(labels, _plot_axis_label(bucket))
                push!(mean_time, row.total_time_mean_s isa Missing ? NaN : Float64(row.total_time_mean_s))
                push!(throughput, row.sim_seconds_per_wall_second_mean isa Missing ? NaN : Float64(row.sim_seconds_per_wall_second_mean))
                push!(success_rate, row.success_rate_pct isa Missing ? NaN : Float64(row.success_rate_pct))
            end
            if !isempty(labels)
                p_time = Plots.bar(
                    labels,
                    mean_time;
                    color="#7f8c8d",
                    ylabel="Mean Runtime [s]",
                    title="Density Backend Benchmark: Runtime",
                    xticks=(1:length(labels), fill("", length(labels))),
                    _plot_margins(size=(2400, 600), bottom_mm=10, right_mm=42)...
                )
                p_tp = Plots.bar(
                    labels,
                    throughput;
                    color="#2f8f5b",
                    ylabel="Sim sec / wall sec",
                    title="Density Backend Benchmark: Throughput",
                    xticks=(1:length(labels), fill("", length(labels))),
                    _plot_margins(size=(2400, 600), bottom_mm=10, right_mm=42)...
                )
                p_sr = Plots.bar(
                    labels,
                    success_rate;
                    color="#3f7fb3",
                    xlabel="Density Backend Bucket",
                    ylabel="Success Rate [%]",
                    title="Density Backend Benchmark: Solve Success",
                    _plot_margins(size=(2400, 760), bottom_mm=64, right_mm=42)...
                )
                plt = Plots.plot(
                    p_time,
                    p_tp,
                    p_sr;
                    layout=(3, 1),
                    size=(2500, 1880),
                    left_margin=20 * Plots.mm,
                    right_margin=20 * Plots.mm
                )
                _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_density_backend", spec, stamp)
            end
        end
    end

    # 20) Paper pack: accuracy-parity suite.
    if _paper_figure_pack_enabled() && !(paper_external === nothing)
        panels = Any[]
        deep_df = paper_external.deep_accuracy_df
        if !(deep_df === nothing) &&
           nrow(deep_df) > 0 &&
           (:rung in names(deep_df)) &&
           (:traj_pos_rel_rms_median_pct in names(deep_df)) &&
           (:traj_vel_rel_rms_median_pct in names(deep_df))
            local_df = copy(deep_df)
            if :mode in names(local_df)
                sort!(local_df, :mode)
            end
            labels = _plot_axis_label.(String.(local_df.rung))
            pos_rms = [v isa Missing ? NaN : Float64(v) for v in local_df.traj_pos_rel_rms_median_pct]
            vel_rms = [v isa Missing ? NaN : Float64(v) for v in local_df.traj_vel_rel_rms_median_pct]
            p_traj = Plots.bar(
                labels,
                hcat(pos_rms, vel_rms);
                label=["Pos RMS rel err (median %)" "Vel RMS rel err (median %)"],
                color=["#2a9d8f" "#3f7fb3"],
                title="Accuracy Parity: Trajectory RMS Relative Error",
                xlabel="Rung",
                ylabel="Relative Error [%]",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
            )
            push!(panels, p_traj)

            if (:periapsis_time_abs_err_p90_s in names(local_df)) && (:interface_time_abs_err_p90_s in names(local_df))
                peri_p90 = [v isa Missing ? NaN : Float64(v) for v in local_df.periapsis_time_abs_err_p90_s]
                interface_p90 = [v isa Missing ? NaN : Float64(v) for v in local_df.interface_time_abs_err_p90_s]
                p_event = Plots.bar(
                    labels,
                    hcat(peri_p90, interface_p90);
                    label=["Periapsis timing abs err (P90)" "Interface timing abs err (P90)"],
                    color=["#d17a4f" "#e9c46a"],
                    title="Accuracy Parity: Event-Time Error",
                    xlabel="Rung",
                    ylabel="Absolute Error [s]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_event)
            end

            if (:propellant_rel_err_p90_pct in names(local_df)) && (:control_impulse_rel_err_p90_pct in names(local_df))
                prop_p90 = [v isa Missing ? NaN : Float64(v) for v in local_df.propellant_rel_err_p90_pct]
                impulse_p90 = [v isa Missing ? NaN : Float64(v) for v in local_df.control_impulse_rel_err_p90_pct]
                p_control = Plots.bar(
                    labels,
                    hcat(prop_p90, impulse_p90);
                    label=["Propellant rel err (P90 %)" "Control impulse rel err (P90 %)"],
                    color=["#8d99ae" "#6d597a"],
                    title="Accuracy Parity: Control/Propellant",
                    xlabel="Rung",
                    ylabel="Relative Error [%]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_control)
            end

            if :callback_exact_match_pct in names(local_df)
                callback_match = [v isa Missing ? NaN : Float64(v) for v in local_df.callback_exact_match_pct]
                p_callback = Plots.bar(
                    labels,
                    callback_match;
                    color="#457b9d",
                    title="Accuracy Parity: Callback-State Exact Match",
                    xlabel="Rung",
                    ylabel="Exact Match [%]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=42, legend=false)...
                )
                push!(panels, p_callback)
            end
        end

        mc_df = paper_external.montecarlo_parity_df
        if !(mc_df === nothing) &&
           nrow(mc_df) > 0 &&
           (:mode in names(mc_df)) &&
           (:rung in names(mc_df)) &&
           (:event_time_ks_distance in names(mc_df))
            agg = combine(
                groupby(mc_df, [:mode, :rung]),
                :event_time_ks_distance => (v -> _safe_stat(v, median)) => :event_ks_median,
                :pos_ks_distance => (v -> _safe_stat(v, median)) => :pos_ks_median,
                :vel_ks_distance => (v -> _safe_stat(v, median)) => :vel_ks_median,
                :mass_ks_distance => (v -> _safe_stat(v, median)) => :mass_ks_median
            )
            sort!(agg, :mode)
            labels = _plot_axis_label.(String.(agg.rung))
            event_ks = [v isa Missing ? NaN : Float64(v) for v in agg.event_ks_median]
            pos_ks = [v isa Missing ? NaN : Float64(v) for v in agg.pos_ks_median]
            vel_ks = [v isa Missing ? NaN : Float64(v) for v in agg.vel_ks_median]
            mass_ks = [v isa Missing ? NaN : Float64(v) for v in agg.mass_ks_median]
            p_mc = Plots.bar(
                labels,
                hcat(event_ks, pos_ks, vel_ks, mass_ks);
                label=["Event-time KS" "Position KS" "Velocity KS" "Mass KS"],
                color=["#bc4749" "#2a9d8f" "#3f7fb3" "#6d597a"],
                title="Accuracy Parity: Monte Carlo Distribution KS Distance",
                xlabel="Rung",
                ylabel="KS Distance",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
            )
            push!(panels, p_mc)
        end

        if isempty(panels)
            gate_parts = DataFrame[]
            if !(split_gate_df === nothing) && nrow(split_gate_df) > 0
                local_df = copy(split_gate_df)
                local_df[!, :gate_label] = fill("split", nrow(local_df))
                push!(gate_parts, local_df)
            end
            if !(multirate_gate_df === nothing) && nrow(multirate_gate_df) > 0
                local_df = copy(multirate_gate_df)
                local_df[!, :gate_label] = fill("multirate", nrow(local_df))
                push!(gate_parts, local_df)
            end
            if !isempty(gate_parts)
                gate_df = vcat(gate_parts...; cols=:union)
                labels = _plot_axis_label.(String.(gate_df.scenario))
                pos_max = [v isa Missing ? NaN : Float64(v) for v in gate_df.pos_rel_max]
                vel_max = [v isa Missing ? NaN : Float64(v) for v in gate_df.vel_rel_max]
                p_gate = Plots.bar(
                    labels,
                    hcat(pos_max, vel_max);
                    label=["Pos rel max" "Vel rel max"],
                    color=["#2a9d8f" "#3f7fb3"],
                    title="Accuracy Gate Fallback: Trajectory Relative Error Maxima",
                    xlabel="Scenario",
                    ylabel="Relative Error",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_gate)
            end
        end

        if !isempty(panels)
            plt = if length(panels) == 1
                panels[1]
            else
                Plots.plot(
                    panels...;
                    layout=(length(panels), 1),
                    size=(2500, 660 * length(panels)),
                    left_margin=20 * Plots.mm,
                    right_margin=20 * Plots.mm
                )
            end
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_accuracy_suite", spec, stamp)
        end
    end

    # 21) Paper pack: cross-machine comparison panel.
    if _paper_figure_pack_enabled() && !(paper_external === nothing)
        panels = Any[]
        speedup_summary_df = paper_external.cross_speedup_summary_df
        if !(speedup_summary_df === nothing) &&
           nrow(speedup_summary_df) > 0 &&
           (:rung in names(speedup_summary_df)) &&
           (:median_speedup_vs_r0 in names(speedup_summary_df))
            local_df = copy(speedup_summary_df)
            sort!(local_df, :median_speedup_vs_r0, rev=true)
            labels = _plot_axis_label.(String.(local_df.rung))
            med_speedup = [v isa Missing ? NaN : Float64(v) for v in local_df.median_speedup_vs_r0]
            p_speedup = Plots.bar(
                labels,
                med_speedup;
                color="#4f9d69",
                title="Cross-Machine Median Speedup vs R0",
                xlabel="Rung",
                ylabel="Median Speedup [x]",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=42, legend=false)...
            )
            Plots.hline!(p_speedup, [1.0]; color=:black, linestyle=:dash, label=false)
            push!(panels, p_speedup)
        end

        regret_summary_df = paper_external.cross_adaptive_regret_summary_df
        if !(regret_summary_df === nothing) &&
           nrow(regret_summary_df) > 0 &&
           (:adaptive_mode in names(regret_summary_df)) &&
           (:mean_time_regret_pct in names(regret_summary_df))
            local_df = regret_summary_df
            labels = _plot_axis_label.(String.(local_df.adaptive_mode))
            regret_pct = [v isa Missing ? NaN : Float64(v) for v in local_df.mean_time_regret_pct]
            win_rate = (:win_rate_pct in names(local_df)) ?
                [v isa Missing ? NaN : Float64(v) for v in local_df.win_rate_pct] :
                fill(NaN, nrow(local_df))
            p_regret = Plots.plot(
                labels,
                regret_pct;
                marker=:circle,
                color="#d17a4f",
                label="Mean time regret [%]",
                title="Cross-Machine Adaptive Regret vs Best Fixed",
                xlabel="Adaptive Mode",
                ylabel="Regret / Win Rate [%]",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
            )
            Plots.plot!(p_regret, labels, win_rate; marker=:diamond, color="#2a9d8f", label="Win rate [%]")
            push!(panels, p_regret)
        end

        route_mix_summary_df = paper_external.cross_route_mix_summary_df
        if !(route_mix_summary_df === nothing) &&
           nrow(route_mix_summary_df) > 0 &&
           (:none_pct_mean in names(route_mix_summary_df)) &&
           (:threads_pct_mean in names(route_mix_summary_df)) &&
           (:process_pct_mean in names(route_mix_summary_df))
            local_df = copy(route_mix_summary_df)
            if :mode in names(local_df)
                mask = [_is_adaptive_mode_token(v) for v in local_df.mode]
                if any(mask)
                    local_df = local_df[mask, :]
                end
            end
            if nrow(local_df) > 0
                labels = if :rung in names(local_df)
                    _plot_axis_label.(String.(local_df.rung))
                elseif :mode in names(local_df)
                    _plot_axis_label.(String.(local_df.mode))
                else
                    _plot_axis_label.(string.("adaptive_", collect(1:nrow(local_df))))
                end
                none_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.none_pct_mean]
                threads_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.threads_pct_mean]
                process_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.process_pct_mean]
                p_mix = Plots.bar(
                    labels,
                    hcat(none_pct, threads_pct, process_pct);
                    label=["none" "threads" "process"],
                    bar_position=:stack,
                    color=["#9aa4b2" "#3f7fb3" "#2f8f5b"],
                    title="Cross-Machine Adaptive Route Mix",
                    xlabel="Adaptive Mode / Rung",
                    ylabel="Mean Route Share [%]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_mix)
            end
        end

        if !isempty(panels)
            plt = if length(panels) == 1
                panels[1]
            else
                Plots.plot(
                    panels...;
                    layout=(length(panels), 1),
                    size=(2500, 660 * length(panels)),
                    left_margin=20 * Plots.mm,
                    right_margin=20 * Plots.mm
                )
            end
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_cross_machine", spec, stamp)
        end
    end

    return plot_artifacts
end

@inline function _fmt(v; digits::Int=3)
    if v isa Missing
        return "n/a"
    elseif v isa AbstractFloat
        return isfinite(v) ? string(round(v; digits=digits)) : "n/a"
    else
        return string(v)
    end
end

function _scenario_metric(summary_df::DataFrame, scenario::String, metric::Symbol)
    idx = findfirst(==(scenario), summary_df.scenario)
    idx === nothing && return nothing
    return summary_df[idx, metric]
end
