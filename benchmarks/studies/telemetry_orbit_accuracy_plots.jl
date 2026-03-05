const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_SUMMARY = joinpath(REPO_ROOT, "output", "telemetry_orbit_accuracy_summary.csv")
const DEFAULT_ERRORS = joinpath(REPO_ROOT, "output", "telemetry_orbit_accuracy_errors.csv")
const DEFAULT_OUTDIR = joinpath(REPO_ROOT, "output", "telemetry_plots")

using CSV
using DataFrames
using Dates
using Printf
using Statistics
using Plots
using Plots.PlotMeasures: mm

function _parse_args(argv::Vector{String})
    summary = DEFAULT_SUMMARY
    errors = DEFAULT_ERRORS
    outdir = DEFAULT_OUTDIR
    for arg in argv
        if startswith(arg, "--summary=")
            summary = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--errors=")
            errors = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        else
            throw(ArgumentError("Unsupported argument '$arg'. Supported: --summary=..., --errors=..., --outdir=..."))
        end
    end
    return (
        summary=abspath(summary),
        errors=abspath(errors),
        outdir=abspath(outdir)
    )
end

@inline function _label_for_row(row)::String
    return string(row.scenario, "/", row.event)
end

function _save_metrics_plot(summary_df::DataFrame, outdir::String)
    ordered = sort(summary_df, [:scenario, :event])
    labels = [_label_for_row(r) for r in eachrow(ordered)]
    x = collect(1:nrow(ordered))
    p = plot(
        title="Telemetry Accuracy Metrics (km)",
        xlabel="scenario/event",
        ylabel="km (log scale)",
        yscale=:log10,
        legend=:topright,
        size=(2200, 1100),
        bottom_margin=24mm,
        left_margin=16mm,
        right_margin=10mm,
        top_margin=8mm,
        dpi=150
    )
    scatter!(p, x, ordered.mae_km; label="MAE", marker=:circle, markersize=6)
    scatter!(p, x, ordered.rmse_km; label="RMSE", marker=:diamond, markersize=6)
    scatter!(p, x, ordered.max_abs_km; label="Max |error|", marker=:star5, markersize=7)
    plot!(p, x, ordered.p95_abs_km; label="P95 |error|", linewidth=2)
    plot!(
        p;
        xticks=(x, labels),
        xrotation=30,
        xtickfont=font(8),
        minorgrid=true,
        gridalpha=0.25
    )
    out = joinpath(outdir, "accuracy_metrics.png")
    savefig(p, out)
    return out
end

function _save_tolerance_margin_plot(summary_df::DataFrame, outdir::String)
    ordered = sort(summary_df, [:scenario, :event])
    labels = [_label_for_row(r) for r in eachrow(ordered)]
    x = collect(1:nrow(ordered))
    abs_ratio = ordered.max_abs_km ./ ordered.limit_max_abs_km
    nmae_ratio = ordered.nmae ./ ordered.limit_nmae
    y_max = max(1.1, maximum(vcat(abs_ratio, nmae_ratio)) * 1.15)

    p = plot(
        title="Tolerance Margin Ratios (<= 1.0 passes)",
        xlabel="scenario/event",
        ylabel="ratio to tolerance",
        legend=:topright,
        size=(2200, 1100),
        bottom_margin=24mm,
        left_margin=16mm,
        right_margin=10mm,
        top_margin=8mm,
        dpi=150
    )
    scatter!(p, x, abs_ratio; label="max_abs / limit_max_abs", marker=:utriangle, markersize=7)
    scatter!(p, x, nmae_ratio; label="nmae / limit_nmae", marker=:diamond, markersize=7)
    hline!(p, [1.0]; linestyle=:dash, color=:red, label="pass threshold")
    plot!(
        p;
        xticks=(x, labels),
        xrotation=30,
        xtickfont=font(8),
        ylims=(0.0, y_max),
        minorgrid=true,
        gridalpha=0.25
    )
    out = joinpath(outdir, "tolerance_margins.png")
    savefig(p, out)
    return out
end

function _save_error_traces_plot(summary_df::DataFrame, errors_df::DataFrame, outdir::String)
    keys = sort(unique(select(summary_df, [:scenario, :event])), [:scenario, :event])
    n = nrow(keys)
    cols = 2
    rows = cld(n, cols)
    p = plot(
        layout=(rows, cols),
        size=(2200, 420 * rows),
        dpi=140,
        bottom_margin=8mm,
        left_margin=10mm,
        right_margin=8mm,
        top_margin=8mm
    )

    for i in 1:n
        scenario = keys.scenario[i]
        event = keys.event[i]
        sub = errors_df[(errors_df.scenario .== scenario) .& (errors_df.event .== event), :]
        if nrow(sub) == 0
            plot!(p[i], title="$(scenario)/$(event): no points", legend=false)
            continue
        end
        sort!(sub, :telemetry_axis)
        axis_units = "index"
        matches = summary_df[(summary_df.scenario .== scenario) .& (summary_df.event .== event), :]
        if nrow(matches) > 0
            axis_units = String(matches.axis_units[1])
        end
        plot!(
            p[i],
            sub.telemetry_axis,
            sub.error_km;
            label="error (km)",
            linewidth=1.3,
            marker=:circle,
            markersize=2.5,
            alpha=0.85
        )
        hline!(p[i], [0.0]; color=:black, linestyle=:dash, linewidth=1, label=false)
        plot!(
            p[i];
            title="$(scenario)/$(event)",
            xlabel="telemetry axis ($(axis_units))",
            ylabel="error (km)",
            legend=false,
            gridalpha=0.25
        )
    end

    out = joinpath(outdir, "error_traces.png")
    savefig(p, out)
    return out
end

function _save_runtime_plot(summary_df::DataFrame, outdir::String)
    by_scenario = combine(
        groupby(summary_df, :scenario),
        :simulation_runtime_s => first => :simulation_runtime_s
    )
    sort!(by_scenario, :simulation_runtime_s; rev=true)

    labels = [String(s) for s in by_scenario.scenario]
    p = bar(
        labels,
        by_scenario.simulation_runtime_s;
        title="Simulation Runtime by Scenario",
        xlabel="scenario",
        ylabel="seconds",
        legend=false,
        xrotation=20,
        size=(1500, 900),
        bottom_margin=18mm,
        left_margin=14mm,
        right_margin=8mm,
        top_margin=8mm,
        dpi=150
    )
    out = joinpath(outdir, "runtime_by_scenario.png")
    savefig(p, out)
    return out
end

function main(argv::Vector{String})
    cfg = _parse_args(argv)
    isfile(cfg.summary) || throw(ArgumentError("Missing summary CSV: $(cfg.summary)"))
    isfile(cfg.errors) || throw(ArgumentError("Missing errors CSV: $(cfg.errors)"))

    summary_df = CSV.read(cfg.summary, DataFrame)
    errors_df = CSV.read(cfg.errors, DataFrame)
    nrow(summary_df) > 0 || throw(ArgumentError("Summary CSV has no rows: $(cfg.summary)"))
    mkpath(cfg.outdir)

    metrics_png = _save_metrics_plot(summary_df, cfg.outdir)
    margins_png = _save_tolerance_margin_plot(summary_df, cfg.outdir)
    traces_png = _save_error_traces_plot(summary_df, errors_df, cfg.outdir)
    runtime_png = _save_runtime_plot(summary_df, cfg.outdir)

    ts = now(UTC)
    total_runtime = hasproperty(summary_df, :total_runtime_s) ? summary_df.total_runtime_s[1] : NaN
    @printf("Plots generated at %s\n", string(ts))
    @printf("Summary rows: %d, error rows: %d\n", nrow(summary_df), nrow(errors_df))
    @printf("Total simulation runtime reported in summary: %.3f s\n", total_runtime)
    println("Files:")
    println("  $metrics_png")
    println("  $margins_png")
    println("  $traces_png")
    println("  $runtime_png")
end

main(ARGS)
