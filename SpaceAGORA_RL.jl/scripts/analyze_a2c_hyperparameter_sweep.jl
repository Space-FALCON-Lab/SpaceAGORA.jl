#!/usr/bin/env julia

ENV["GKSwstype"] = "100"

using CSV
using DataFrames
using Plots
using Statistics
using TOML

const PACKAGE_DIR = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_SWEEP_ROOT = joinpath(
    PACKAGE_DIR,
    "outputs",
    "hyperparameter_sweeps",
    "a2c_one_factor_20260814",
)

const PARAMETER_SPECS = Dict(
    "entropy_coefficient" => (
        directory="entropy_coefficient",
        parameter="Entropy coefficient",
        default_nominal_value="0.01",
    ),
    "successful_case_repetitions" => (
        directory="successful_case_repetitions",
        parameter="Successful case repetitions",
        default_nominal_value="9",
    ),
    "learning_rate" => (
        directory="learning_rate",
        parameter="Learning rate",
        default_nominal_value="1e-4",
    ),
    "rollout_length" => (
        directory="rollout_length",
        parameter="Rollout length",
        default_nominal_value="10",
    ),
)

function _usage(io::IO=stdout)
    println(io, """
Usage:
  julia --project=SpaceAGORA_RL.jl \\
    SpaceAGORA_RL.jl/scripts/analyze_a2c_hyperparameter_sweep.jl \\
    [SWEEP_ROOT] [COMMON_EVALUATION_DIR]

The defaults are the a2c_one_factor_20260814 sweep directory and its
common_evaluation subdirectory.
""")
end

function _variant_records(definition)
    return Dict(String(record["name"]) => record for record in definition["variants"])
end

function _comparison_specs(definition, variants)
    baseline_name = String(get(definition, "baseline_variant", "nominal"))
    haskey(variants, baseline_name) || error("unknown baseline variant $baseline_name")
    comparison_names = if haskey(definition, "comparison_variants")
        String.(definition["comparison_variants"])
    else
        [
            String(record["name"])
            for record in definition["variants"]
            if String(record["name"]) != baseline_name
        ]
    end
    return map(comparison_names) do name
        haskey(variants, name) || error("unknown comparison variant $name")
        record = variants[name]
        parameter_name = String(record["parameter"])
        haskey(PARAMETER_SPECS, parameter_name) ||
            error("unsupported sweep parameter $parameter_name")
        parameter = PARAMETER_SPECS[parameter_name]
        return (
            variant=name,
            directory=String(get(record, "analysis_directory", parameter.directory)),
            parameter=parameter.parameter,
            nominal_value=String(get(
                record,
                "nominal_value",
                parameter.default_nominal_value,
            )),
            changed_value=String(record["changed_value"]),
        )
    end
end

function _completed_runs(runs_dir, variants)
    runs = Dict{String,String}()
    for (name, record) in variants
        if haskey(record, "run_dir")
            run_dir = abspath(String(record["run_dir"]))
            isfile(joinpath(run_dir, "manifest.toml")) ||
                error("recorded run for sweep variant $name has no manifest: $run_dir")
            isfile(joinpath(run_dir, "checkpoint_validation", "best_validation_checkpoint.txt")) ||
                error("recorded run for sweep variant $name has no best checkpoint: $run_dir")
            runs[name] = run_dir
            continue
        end
        config_path = abspath(String(record["config_path"]))
        matches = String[]
        for run_dir in readdir(runs_dir; join=true)
            manifest_path = joinpath(run_dir, "manifest.toml")
            isfile(manifest_path) || continue
            manifest = TOML.parsefile(manifest_path)
            abspath(String(get(manifest, "config_path", ""))) == config_path || continue
            isfile(joinpath(run_dir, "checkpoint_validation", "best_validation_checkpoint.txt")) ||
                continue
            push!(matches, run_dir)
        end
        isempty(matches) && error("no completed run for sweep variant $name")
        runs[name] = last(sort(matches))
    end
    return runs
end

function _evaluation_sources(manifest)
    return Dict(
        abspath(String(source["config_path"])) => source
        for source in manifest["source_runs"]
    )
end

function _comparison_row(comparison, source)
    rows = comparison[comparison.policy .== String(source["policy_label"]), :]
    nrow(rows) == 1 || error("expected one comparison row for $(source["policy_label"])")
    return rows[1, :]
end

function _validation_plot(path, comparison, nominal_run, changed_run,
                          nominal_step, changed_step)
    nominal = CSV.read(
        joinpath(nominal_run, "checkpoint_validation", "checkpoint_validation_summary.csv"),
        DataFrame,
    )
    changed = CSV.read(
        joinpath(changed_run, "checkpoint_validation", "checkpoint_validation_summary.csv"),
        DataFrame,
    )
    nominal_label = "Nominal ($(comparison.nominal_value))"
    changed_label = "Changed ($(comparison.changed_value))"

    p = plot(
        xlabel="Training step (thousands)",
        ylabel="Success rate (%)",
        title="$(comparison.parameter): checkpoint validation",
        ylims=(0, 100),
        legend=:outerbottom,
        legend_columns=2,
        size=(1050, 650),
        left_margin=7Plots.mm,
        bottom_margin=10Plots.mm,
    )
    colors = Dict("conservative" => :steelblue, "tolerant" => :darkorange)
    for (frame, role_label, style) in (
        (nominal, nominal_label, :solid),
        (changed, changed_label, :dash),
    )
        for mode in ("conservative", "tolerant")
            rows = frame[frame.mode .== mode, :]
            plot!(
                p,
                rows.global_step ./ 1000,
                rows.success_percent;
                label="$role_label — $mode",
                color=colors[mode],
                linestyle=style,
                linewidth=2,
                marker=:circle,
                markersize=3,
            )
        end
    end
    for (frame, step, color) in (
        (nominal, nominal_step, :steelblue),
        (changed, changed_step, :darkorange),
    )
        rows = frame[frame.mode .== "conservative", :]
        index = argmin(abs.(rows.global_step .- step))
        scatter!(
            p,
            [rows.global_step[index] / 1000],
            [rows.success_percent[index]];
            label=false,
            color=color,
            marker=:star5,
            markersize=10,
        )
    end
    savefig(p, path)
end

function _metric_panel(labels, values, title, ylabel; ylims=:auto)
    return bar(
        labels,
        values;
        title=title,
        ylabel=ylabel,
        legend=false,
        color=[:steelblue, :darkorange, :seagreen],
        xrotation=15,
        ylims=ylims,
        bottom_margin=10Plots.mm,
        left_margin=6Plots.mm,
    )
end

function _performance_plot(path, comparison, nominal, changed, pr_drl)
    labels = [
        "Nominal\n($(comparison.nominal_value))",
        "Changed\n($(comparison.changed_value))",
        "PR-DRL",
    ]
    rows = [nominal, changed, pr_drl]
    panels = [
        _metric_panel(
            labels,
            Float64[row.target_reached_10km_percent for row in rows],
            "Target reached within 10 km",
            "Campaigns (%)";
            ylims=(0, 100),
        ),
        _metric_panel(
            labels,
            Float64[row.mean_abs_final_target_error_km for row in rows],
            "Final target error",
            "Kilometers",
        ),
        _metric_panel(
            labels,
            Float64[row.mean_thermal_violations for row in rows],
            "Thermal violations",
            "Mean per campaign",
        ),
        _metric_panel(
            labels,
            Float64[row.mean_total_delta_v_mps for row in rows],
            "Total aerobraking ΔV",
            "m/s",
        ),
        _metric_panel(
            labels,
            Float64[row.mean_maneuver_count for row in rows],
            "Maneuver count",
            "Mean per campaign",
        ),
        _metric_panel(
            labels,
            Float64[row.mean_duration_days for row in rows],
            "Campaign duration",
            "Days",
        ),
    ]
    figure = plot(
        panels...;
        layout=(2, 3),
        size=(1500, 900),
        plot_title="$(comparison.parameter): nominal vs changed best checkpoint",
    )
    savefig(figure, path)
end

function _metrics_record(role, label, source, row)
    return (
        role=role,
        label=label,
        run_dir=String(source["run_dir"]),
        checkpoint_path=String(source["checkpoint_path"]),
        checkpoint_global_step=Int(source["checkpoint_global_step"]),
        episodes=Int(row.episodes),
        target_reached_10km_percent=Float64(row.target_reached_10km_percent),
        mean_abs_final_target_error_km=Float64(row.mean_abs_final_target_error_km),
        mean_total_delta_v_mps=Float64(row.mean_total_delta_v_mps),
        mean_thermal_violations=Float64(row.mean_thermal_violations),
        mean_maneuver_count=Float64(row.mean_maneuver_count),
        mean_duration_days=Float64(row.mean_duration_days),
    )
end

function _paired_success_stats(episode_metrics, nominal_source, changed_source)
    nominal_key = String(nominal_source["policy_key"])
    changed_key = String(changed_source["policy_key"])
    nominal = sort(episode_metrics[episode_metrics.policy .== nominal_key, :], :episode)
    changed = sort(episode_metrics[episode_metrics.policy .== changed_key, :], :episode)
    nominal.episode == changed.episode || error("evaluation episodes are not paired")
    nominal_success = Bool.(nominal.success)
    changed_success = Bool.(changed.success)
    differences = Float64.(changed_success) .- Float64.(nominal_success)
    standard_error = length(differences) > 1 ? std(differences) / sqrt(length(differences)) : 0.0
    return (
        changed_only_successes=count(changed_success .& .!nominal_success),
        nominal_only_successes=count(nominal_success .& .!changed_success),
        paired_delta_success_pp=100 * mean(differences),
        paired_delta_95ci_low_pp=100 * (mean(differences) - 1.96 * standard_error),
        paired_delta_95ci_high_pp=100 * (mean(differences) + 1.96 * standard_error),
    )
end

function _summary_record(comparison, nominal, changed, paired)
    percent_delta(changed_value, nominal_value) =
        nominal_value == 0 ? NaN : 100 * (changed_value - nominal_value) / nominal_value
    return (
        parameter=comparison.parameter,
        nominal_value=comparison.nominal_value,
        changed_value=comparison.changed_value,
        nominal_checkpoint_step=nominal.checkpoint_global_step,
        changed_checkpoint_step=changed.checkpoint_global_step,
        nominal_success_percent=nominal.target_reached_10km_percent,
        changed_success_percent=changed.target_reached_10km_percent,
        delta_success_pp=changed.target_reached_10km_percent -
                         nominal.target_reached_10km_percent,
        changed_only_successes=paired.changed_only_successes,
        nominal_only_successes=paired.nominal_only_successes,
        paired_delta_95ci_low_pp=paired.paired_delta_95ci_low_pp,
        paired_delta_95ci_high_pp=paired.paired_delta_95ci_high_pp,
        nominal_target_error_km=nominal.mean_abs_final_target_error_km,
        changed_target_error_km=changed.mean_abs_final_target_error_km,
        target_error_change_percent=percent_delta(
            changed.mean_abs_final_target_error_km,
            nominal.mean_abs_final_target_error_km,
        ),
        nominal_thermal_violations=nominal.mean_thermal_violations,
        changed_thermal_violations=changed.mean_thermal_violations,
        thermal_violations_change_percent=percent_delta(
            changed.mean_thermal_violations,
            nominal.mean_thermal_violations,
        ),
        nominal_delta_v_mps=nominal.mean_total_delta_v_mps,
        changed_delta_v_mps=changed.mean_total_delta_v_mps,
        delta_v_change_percent=percent_delta(
            changed.mean_total_delta_v_mps,
            nominal.mean_total_delta_v_mps,
        ),
        nominal_maneuver_count=nominal.mean_maneuver_count,
        changed_maneuver_count=changed.mean_maneuver_count,
        maneuver_count_change_percent=percent_delta(
            changed.mean_maneuver_count,
            nominal.mean_maneuver_count,
        ),
        nominal_duration_days=nominal.mean_duration_days,
        changed_duration_days=changed.mean_duration_days,
        duration_change_percent=percent_delta(
            changed.mean_duration_days,
            nominal.mean_duration_days,
        ),
    )
end

function main(args=ARGS)
    any(==("--help"), args) && return _usage()
    length(args) <= 2 || throw(ArgumentError("too many arguments; use --help for usage"))
    sweep_root = abspath(isempty(args) ? DEFAULT_SWEEP_ROOT : args[1])
    evaluation_dir = abspath(length(args) < 2 ? joinpath(sweep_root, "common_evaluation") : args[2])

    definition = TOML.parsefile(joinpath(sweep_root, "sweep_definition.toml"))
    variants = _variant_records(definition)
    comparisons = _comparison_specs(definition, variants)
    baseline_name = String(get(definition, "baseline_variant", "nominal"))
    runs = _completed_runs(joinpath(sweep_root, "runs"), variants)
    evaluation_manifest = TOML.parsefile(joinpath(evaluation_dir, "evaluation_manifest.toml"))
    sources = _evaluation_sources(evaluation_manifest)
    comparison_metrics = CSV.read(
        joinpath(evaluation_dir, "multi_run_flight_performance_comparison.csv"),
        DataFrame,
    )
    episode_metrics = CSV.read(
        joinpath(evaluation_dir, "trained_policies_vs_aads", "episode_metrics.csv"),
        DataFrame,
    )

    source_for(name) = sources[abspath(String(variants[name]["config_path"]))]
    nominal_source = source_for(baseline_name)
    nominal_row = _comparison_row(comparison_metrics, nominal_source)
    pr_source = only(source for source in values(sources) if source["algorithm"] == "pr_drl")
    pr_row = _comparison_row(comparison_metrics, pr_source)

    analysis_dir = joinpath(sweep_root, "analysis")
    mkpath(analysis_dir)
    summaries = NamedTuple[]
    for comparison in comparisons
        changed_source = source_for(comparison.variant)
        changed_row = _comparison_row(comparison_metrics, changed_source)
        output_dir = joinpath(analysis_dir, comparison.directory)
        mkpath(output_dir)

        nominal_record = _metrics_record(
            "nominal",
            "Nominal ($(comparison.parameter) = $(comparison.nominal_value))",
            nominal_source,
            nominal_row,
        )
        changed_record = _metrics_record(
            "changed",
            "Changed ($(comparison.parameter) = $(comparison.changed_value))",
            changed_source,
            changed_row,
        )
        pr_record = _metrics_record("benchmark", "PR-DRL", pr_source, pr_row)
        CSV.write(
            joinpath(output_dir, "best_checkpoint_evaluation_metrics.csv"),
            DataFrame([nominal_record, changed_record, pr_record]),
        )

        _validation_plot(
            joinpath(output_dir, "checkpoint_validation_success.png"),
            comparison,
            runs[baseline_name],
            runs[comparison.variant],
            nominal_record.checkpoint_global_step,
            changed_record.checkpoint_global_step,
        )
        _performance_plot(
            joinpath(output_dir, "best_checkpoint_flight_performance.png"),
            comparison,
            nominal_row,
            changed_row,
            pr_row,
        )
        paired = _paired_success_stats(episode_metrics, nominal_source, changed_source)
        push!(summaries, _summary_record(comparison, nominal_record, changed_record, paired))
    end
    summary_path = joinpath(analysis_dir, "relative_performance_summary.csv")
    CSV.write(summary_path, DataFrame(summaries))
    println("wrote hyperparameter sweep analysis to ", analysis_dir)
    println("relative performance summary: ", summary_path)
    return analysis_dir
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch error
        showerror(stderr, error, catch_backtrace())
        println(stderr)
        exit(1)
    end
end
