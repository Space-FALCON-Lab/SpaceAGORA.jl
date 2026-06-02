#!/usr/bin/env julia

# Plot ORACLE Case 2 open-cavity laser-link time-history grids.
#
# Usage:
#   julia --project=. scripts/plotting/oracle_case2_laser_grid_plots.jl \
#     --input output/oracle_case2_laser_links/case2_laser_timeseries.csv

using Arrow
using CSV
using DataFrames
get!(ENV, "GKSwstype", "100")
using Plots
using Printf

gr()

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_INPUT = joinpath(REPO_ROOT, "output", "oracle_case2_laser_links", "case2_laser_timeseries.csv")
const DEFAULT_OUTDIR = joinpath(REPO_ROOT, "output", "oracle_case2_laser_links", "plots")
const PAPER_ALTITUDE_ORDER_KM = [1150.0, 1050.0, 1000.0, 950.0, 850.0]
const PAPER_INCLINATION_DELTAS_DEG = [0.0, 0.5, 1.0]

Base.@kwdef struct OracleCase2PlotOptions
    input::String = DEFAULT_INPUT
    outdir::String = DEFAULT_OUTDIR
    format::String = "both"
    plot_set::String = "all"
    schedule::String = "all"
    helpers::Int = 0
    time_axis::String = "orbit"
    element_mode::String = "laser"
end

function _usage()
    return """
    Usage:
      julia --project=. scripts/plotting/oracle_case2_laser_grid_plots.jl [options]

    Options:
      --input PATH              Time-series CSV, Feather, or Arrow file from ORACLE/run_case2_laser_links.jl.
      --outdir PATH             Output directory for figures.
      --format png|pdf|both     Figure format to write. Default: both.
      --plot-set all|velocity|elements
      --schedule NAME|all       Filter schedule column. Default: all.
      --helpers N               Filter helper count. Default: max helper count in the file.
      --time-axis orbit|seconds Plot x-axis in target orbits or seconds. Default: orbit.
      --element-mode laser|total Use laser-only element deltas or total deltas. Default: laser.
    """
end

function _parse_options(argv)::OracleCase2PlotOptions
    opts = Dict{Symbol, Any}()
    i = 1
    while i <= length(argv)
        arg = argv[i]
        if arg == "--help" || arg == "-h"
            println(_usage())
            exit(0)
        end
        startswith(arg, "--") || throw(ArgumentError("Unexpected argument '$arg'.\n$(_usage())"))
        i < length(argv) || throw(ArgumentError("Missing value for $arg."))
        key = Symbol(replace(arg[3:end], '-' => '_'))
        val = argv[i + 1]
        if key in (:input, :outdir)
            opts[key] = abspath(val)
        elseif key == :helpers
            opts[key] = parse(Int, val)
        elseif key in (:format, :plot_set, :schedule, :time_axis, :element_mode)
            opts[key] = val
        else
            throw(ArgumentError("Unknown option $arg.\n$(_usage())"))
        end
        i += 2
    end
    parsed = OracleCase2PlotOptions(; opts...)
    parsed.format in ("png", "pdf", "both") ||
        throw(ArgumentError("--format must be png, pdf, or both."))
    parsed.plot_set in ("all", "velocity", "elements") ||
        throw(ArgumentError("--plot-set must be all, velocity, or elements."))
    parsed.time_axis in ("orbit", "seconds") ||
        throw(ArgumentError("--time-axis must be orbit or seconds."))
    parsed.element_mode in ("laser", "total") ||
        throw(ArgumentError("--element-mode must be laser or total."))
    return parsed
end

function _read_table(path::String)::DataFrame
    isfile(path) || throw(ArgumentError("Missing ORACLE Case 2 time-series file: $path"))
    ext = lowercase(splitext(path)[2])
    if ext == ".csv"
        return CSV.read(path, DataFrame)
    elseif ext in (".feather", ".arrow")
        return DataFrame(Arrow.Table(path))
    end
    throw(ArgumentError("Unsupported input extension '$ext'. Use .csv, .feather, or .arrow."))
end

function _require_columns(df::DataFrame, names_needed::Vector{Symbol})
    missing = [name for name in names_needed if name ∉ propertynames(df)]
    isempty(missing) || throw(ArgumentError("Input file is missing required columns: $(join(string.(missing), ", "))"))
    return nothing
end

function _prepare_timeseries(df::DataFrame, opts::OracleCase2PlotOptions)::DataFrame
    _require_columns(df, [
        :time_s, :orbit, :helpers, :helper_altitude_km, :helper_inclination_delta_deg, :schedule,
        :dv_r_mps, :dv_t_mps, :dv_n_mps,
        :da_m, :de, :di_deg, :draan_deg,
        :laser_da_m, :laser_de, :laser_di_deg, :laser_draan_deg,
    ])
    nrow(df) > 0 || throw(ArgumentError("Input file has no rows: $(opts.input)"))
    if opts.schedule != "all"
        df = df[string.(df.schedule) .== opts.schedule, :]
        nrow(df) > 0 || throw(ArgumentError("No rows match --schedule $(opts.schedule)."))
    end
    helpers = opts.helpers
    if helpers <= 0
        helpers = maximum(Int.(df.helpers))
        @printf("Using helpers=%d from %s\n", helpers, opts.input)
    end
    df = df[Int.(df.helpers) .== helpers, :]
    nrow(df) > 0 || throw(ArgumentError("No rows match helpers=$helpers."))
    sort!(df, [:helper_altitude_km, :helper_inclination_delta_deg, :time_s])
    return df
end

function _ordered_present(values, preferred_order::Vector{Float64}; rev::Bool=false)
    present = sort(unique(Float64.(values)); rev=rev)
    ordered = [value for value in preferred_order if any(isapprox(value, x; atol=1e-8) for x in present)]
    extras = [value for value in present if !any(isapprox(value, x; atol=1e-8) for x in ordered)]
    return vcat(ordered, extras)
end

function _axes(df::DataFrame)
    altitudes = _ordered_present(df.helper_altitude_km, PAPER_ALTITUDE_ORDER_KM; rev=true)
    inclinations = _ordered_present(df.helper_inclination_delta_deg, PAPER_INCLINATION_DELTAS_DEG)
    return altitudes, inclinations
end

@inline _x_column(opts::OracleCase2PlotOptions)::Symbol = opts.time_axis == "orbit" ? :orbit : :time_s
@inline _x_label(opts::OracleCase2PlotOptions)::String = opts.time_axis == "orbit" ? "time [target orbits]" : "time [s]"

function _case_title(df::DataFrame)::String
    schedule = length(unique(string.(df.schedule))) == 1 ? string(df.schedule[1]) : "mixed schedules"
    helpers = length(unique(Int.(df.helpers))) == 1 ? string(Int(df.helpers[1])) : "mixed"
    target = :target_altitude_km in propertynames(df) ? @sprintf("%.0f km", Float64(df.target_altitude_km[1])) : "1000 km"
    return "ORACLE Case 2, target altitude $target, helpers=$helpers, $schedule"
end

function _subdf(df::DataFrame, altitude_km::Float64, dinc_deg::Float64)::DataFrame
    mask = map(eachrow(df)) do row
        isapprox(Float64(row.helper_altitude_km), altitude_km; atol=1e-8) &&
        isapprox(Float64(row.helper_inclination_delta_deg), dinc_deg; atol=1e-8)
    end
    return df[mask, :]
end

function _save_formats(fig, outbase::String, format::String)::Vector{String}
    written = String[]
    if format in ("png", "both")
        path = outbase * ".png"
        savefig(fig, path)
        push!(written, path)
    end
    if format in ("pdf", "both")
        path = outbase * ".pdf"
        savefig(fig, path)
        push!(written, path)
    end
    return written
end

function plot_velocity_grid(df::DataFrame, opts::OracleCase2PlotOptions)::Vector{String}
    altitudes, inclinations = _axes(df)
    xcol = _x_column(opts)
    subplots = Plots.Plot[]
    for altitude in altitudes
        for dinc in inclinations
            sub = _subdf(df, altitude, dinc)
            title = @sprintf("helper altitude %.0f km, dinc %.1f deg", altitude, dinc)
            p = plot(
                title=title,
                xlabel=altitude == last(altitudes) ? _x_label(opts) : "",
                ylabel=dinc == first(inclinations) ? "dV [m/s]" : "",
                legend=(altitude == first(altitudes) && dinc == first(inclinations)) ? :best : false,
                framestyle=:box,
                gridalpha=0.25,
            )
            if nrow(sub) > 0
                sort!(sub, :time_s)
                plot!(p, sub[!, xcol], sub.dv_r_mps; label="R", linewidth=1.7, color=:red)
                plot!(p, sub[!, xcol], sub.dv_t_mps; label="T", linewidth=1.7, color=:blue)
                plot!(p, sub[!, xcol], sub.dv_n_mps; label="N", linewidth=1.7, color=:green)
            end
            push!(subplots, p)
        end
    end
    fig = plot(
        subplots...;
        layout=(length(altitudes), length(inclinations)),
        size=(520 * length(inclinations), 285 * length(altitudes)),
        dpi=160,
        plot_title=_case_title(df) * " velocity RTN",
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=7Plots.mm,
    )
    return _save_formats(fig, joinpath(opts.outdir, "oracle_case2_velocity_rtn_time_grid"), opts.format)
end

function _element_metrics(opts::OracleCase2PlotOptions)
    if opts.element_mode == "laser"
        return (
            (:laser_da_m, "da [m]"),
            (:laser_de, "de"),
            (:laser_di_deg, "di [deg]"),
            (:laser_draan_deg, "dRAAN [deg]"),
        )
    end
    return (
        (:da_m, "da [m]"),
        (:de, "de"),
        (:di_deg, "di [deg]"),
        (:draan_deg, "dRAAN [deg]"),
    )
end

function _line_style_for_inclination(idx::Int)
    colors = (:black, :blue, :red, :green, :purple)
    styles = (:solid, :dash, :dot, :dashdot, :solid)
    return colors[mod1(idx, length(colors))], styles[mod1(idx, length(styles))]
end

function plot_element_grid(df::DataFrame, opts::OracleCase2PlotOptions)::Vector{String}
    altitudes, inclinations = _axes(df)
    metrics = _element_metrics(opts)
    xcol = _x_column(opts)
    subplots = Plots.Plot[]
    for altitude in altitudes
        for (metric, label) in metrics
            p = plot(
                title=@sprintf("helper altitude %.0f km, %s", altitude, label),
                xlabel=altitude == last(altitudes) ? _x_label(opts) : "",
                ylabel=metric == metrics[1][1] ? label : "",
                legend=(altitude == first(altitudes) && metric == metrics[1][1]) ? :best : false,
                framestyle=:box,
                gridalpha=0.25,
            )
            for (idx, dinc) in enumerate(inclinations)
                sub = _subdf(df, altitude, dinc)
                nrow(sub) == 0 && continue
                sort!(sub, :time_s)
                color, linestyle = _line_style_for_inclination(idx)
                plot!(
                    p,
                    sub[!, xcol],
                    sub[!, metric];
                    label=@sprintf("dinc %.1f deg", dinc),
                    linewidth=1.7,
                    color=color,
                    linestyle=linestyle,
                )
            end
            push!(subplots, p)
        end
    end
    suffix = opts.element_mode == "laser" ? "laser_only" : "total"
    fig = plot(
        subplots...;
        layout=(length(altitudes), length(metrics)),
        size=(520 * length(metrics), 285 * length(altitudes)),
        dpi=160,
        plot_title=_case_title(df) * " orbital element changes",
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=7Plots.mm,
    )
    return _save_formats(fig, joinpath(opts.outdir, "oracle_case2_orbital_elements_time_grid_$suffix"), opts.format)
end

function main(argv=ARGS)
    opts = _parse_options(String.(argv))
    mkpath(opts.outdir)
    raw_df = _read_table(opts.input)
    df = _prepare_timeseries(raw_df, opts)

    written = String[]
    if opts.plot_set in ("all", "velocity")
        append!(written, plot_velocity_grid(df, opts))
    end
    if opts.plot_set in ("all", "elements")
        append!(written, plot_element_grid(df, opts))
    end

    println("ORACLE Case 2 plots generated from $(opts.input)")
    println("Rows plotted: $(nrow(df))")
    println("Files:")
    foreach(path -> println("  $path"), written)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
