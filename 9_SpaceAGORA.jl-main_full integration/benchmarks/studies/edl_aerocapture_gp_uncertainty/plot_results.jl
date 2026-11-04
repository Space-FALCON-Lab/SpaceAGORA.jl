const STUDY_ROOT = @__DIR__
const DEFAULT_SUMMARY = joinpath(STUDY_ROOT, "output", "summary_metrics.csv")
const DEFAULT_PLOT_DIR = joinpath(STUDY_ROOT, "plots")

using CSV
using DataFrames
using Plots
using Statistics

function _require_summary(path::String)
    isfile(path) || error("Missing summary metrics CSV at $path")
    return CSV.read(path, DataFrame)
end

function _require_pointwise(path::String)
    isfile(path) || error("Missing pointwise predictions CSV at $path. Run the study with --save-pointwise true.")
    return CSV.read(path, DataFrame)
end

function _aggregate(df::DataFrame, metric::Symbol)
    grouped = combine(groupby(df, [:parameterization, :kernel, :mean_basis])) do sub
        DataFrame(value=mean(sub[!, metric]))
    end
    grouped.label = string.(grouped.parameterization, " / ", grouped.kernel, " / ", grouped.mean_basis)
    sort!(grouped, :value)
    return grouped
end

function _model_label(row)
    return row.kernel == "gram_baseline" ? "GRAM prior" : "$(row.parameterization) / $(row.kernel) / $(row.mean_basis)"
end

function _plot_metric_bars(df::DataFrame, metric::Symbol, out_dir::String)
    agg = _aggregate(df, metric)
    agg.label = [_model_label(row) for row in eachrow(agg)]
    plt = bar(
        agg.label,
        agg.value;
        legend=false,
        xlabel="Model",
        ylabel=string(metric),
        xrotation=35,
        size=(1100, 600),
        title="EDL aerocapture GP study: $(metric)",
    )
    savefig(plt, joinpath(out_dir, string(metric, ".png")))
end

function _plot_profiles(summary_df::DataFrame, pointwise_df::DataFrame, out_dir::String, case_id::AbstractString)
    case_metrics = filter(:case_id => ==(case_id), summary_df)
    candidates = filter(:kernel => !=("gram_baseline"), case_metrics)
    best = first(sort(candidates, :weighted_rmse_log))
    gp = filter(
        row -> row.case_id == case_id &&
               row.parameterization == best.parameterization &&
               row.kernel == best.kernel &&
               row.mean_basis == best.mean_basis,
        pointwise_df,
    )
    baseline = filter(
        row -> row.case_id == case_id && row.kernel == "gram_baseline",
        pointwise_df,
    )
    nrow(gp) == nrow(baseline) || error("Pointwise baseline and GP predictions do not align for $case_id.")
    sort!(gp, :alt_m)
    sort!(baseline, :alt_m)
    alt_km = gp.alt_m .* 1.0e-3
    truth = gp.truth_density

    density_plot = plot(
        truth,
        alt_km;
        label="Seeded GRAM truth",
        xscale=:log10,
        linewidth=3,
        xlabel="Density (kg/m³)",
        ylabel="Altitude (km)",
        title="$case_id: density prediction comparison",
        size=(900, 700),
    )
    plot!(density_plot, baseline.pred_mean_density, alt_km; label="GRAM prior", linewidth=2, linestyle=:dash)
    plot!(density_plot, gp.pred_mean_density, alt_km; label="Best GP: $(_model_label(best))", linewidth=2)
    savefig(density_plot, joinpath(out_dir, "density_prediction_comparison.png"))

    baseline_error = 100.0 .* (baseline.pred_mean_density .- truth) ./ truth
    gp_error = 100.0 .* (gp.pred_mean_density .- truth) ./ truth
    error_plot = plot(
        baseline_error,
        alt_km;
        label="GRAM prior",
        linewidth=2,
        linestyle=:dash,
        xlabel="Relative density error (%)",
        ylabel="Altitude (km)",
        title="$case_id: density prediction error",
        size=(900, 700),
    )
    vline!(error_plot, [0.0]; label=false, color=:black, alpha=0.5)
    plot!(error_plot, gp_error, alt_km; label="Best GP: $(_model_label(best))", linewidth=2)
    savefig(error_plot, joinpath(out_dir, "density_relative_error.png"))
    return best
end

function plot_summary(
    ;
    summary_path::String=DEFAULT_SUMMARY,
    pointwise_path::String=joinpath(dirname(summary_path), "pointwise_predictions.csv"),
    out_dir::String=DEFAULT_PLOT_DIR,
    case_id::Union{Nothing, AbstractString}=nothing,
)
    df = _require_summary(summary_path)
    mkpath(out_dir)

    for metric in (:weighted_rmse, :weighted_rmse_log, :weighted_nlpd, :weighted_coverage_1sigma, :weighted_coverage_2sigma)
        _plot_metric_bars(df, metric, out_dir)
    end
    selected_case = isnothing(case_id) ? first(df.case_id) : case_id
    best = _plot_profiles(df, _require_pointwise(pointwise_path), out_dir, selected_case)
    return (out_dir=out_dir, case_id=selected_case, best_model=(best.parameterization, best.kernel, best.mean_basis))
end

function _plot_args(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        startswith(args[i], "--") || error("Unsupported argument '$(args[i])'.")
        body = args[i][3:end]
        if occursin("=", body)
            k, v = split(body, "="; limit=2)
            opts[k] = v
            i += 1
        else
            i + 1 <= length(args) || error("Missing value for --$body.")
            opts[body] = args[i + 1]
            i += 2
        end
    end
    return opts
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = _plot_args(collect(ARGS))
    summary_path = get(opts, "summary-path", DEFAULT_SUMMARY)
    plot_summary(;
        summary_path,
        pointwise_path=get(opts, "pointwise-path", joinpath(dirname(summary_path), "pointwise_predictions.csv")),
        out_dir=get(opts, "out-dir", DEFAULT_PLOT_DIR),
        case_id=get(opts, "case-id", nothing),
    )
end
