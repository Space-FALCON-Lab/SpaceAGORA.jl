module PeriapsisCorridorPlots

using CSV
using DataFrames
using Plots

export write_corridor_plots

const _SECONDS_PER_DAY = 86_400.0

function _read_case(case_dir::AbstractString)
    passes_path = joinpath(case_dir, "passes.csv")
    summary_path = joinpath(case_dir, "summary.csv")
    isfile(passes_path) && isfile(summary_path) || return nothing
    burns_path = joinpath(case_dir, "burns.csv")
    passes = CSV.read(passes_path, DataFrame)
    summary = CSV.read(summary_path, DataFrame)
    burns_present = isfile(burns_path)
    burns = burns_present ? CSV.read(burns_path, DataFrame) : DataFrame()
    isempty(summary) && return nothing
    return (name=basename(case_dir), passes=passes, burns=burns,
            burns_present=burns_present, summary=summary)
end

function _case_data(output_root::AbstractString)
    isdir(output_root) || return NamedTuple[]
    cases = NamedTuple[]
    for name in sort(readdir(output_root))
        case_dir = joinpath(output_root, name)
        isdir(case_dir) || continue
        data = _read_case(case_dir)
        data === nothing || push!(cases, data)
    end
    return cases
end

function _column(df::AbstractDataFrame, names::Symbol...)
    for name in names
        name in propertynames(df) && return name
    end
    return nothing
end

function _series_xy(df::AbstractDataFrame, x_name::Symbol, y_name::Symbol)
    x_values, y_values = Float64[], Float64[]
    for (x, y) in zip(df[!, x_name], df[!, y_name])
        (ismissing(x) || ismissing(y)) && continue
        xf, yf = Float64(x), Float64(y)
        isfinite(xf) && isfinite(yf) || continue
        push!(x_values, xf)
        push!(y_values, yf)
    end
    return x_values, y_values
end

function _case_label(case)
    summary = case.summary
    status_col = _column(summary, :mission_result, :run_status)
    status = status_col === nothing || ismissing(summary[1, status_col]) ? "" :
        " ($(summary[1, status_col]))"
    smoke_col = _column(summary, :smoke)
    smoke = smoke_col !== nothing && !ismissing(summary[1, smoke_col]) &&
        Bool(summary[1, smoke_col])
    return case.name * (smoke ? " (smoke)" : status)
end

function _physical_periapsis_plot(output_root, cases, config)
    mission = config.values["mission"]
    minimum_altitude = Float64(mission["corridor_min_altitude_m"])
    maximum_altitude = Float64(mission["corridor_max_altitude_m"])
    path = joinpath(output_root, "physical_periapsis_corridor.png")
    p = plot(; title="Physical periapsis at radial events",
             xlabel="Elapsed time (days)", ylabel="Periapsis geodetic altitude (km)",
             legend=:best, grid=true)
    hline!(p, [minimum_altitude / 1000]; label="Lower corridor", color=:firebrick,
           linestyle=:dash, linewidth=1.5)
    hline!(p, [maximum_altitude / 1000]; label="Upper corridor", color=:firebrick,
           linestyle=:dash, linewidth=1.5)
    plotted = false
    for case in cases
        df = case.passes
        :time_s in propertynames(df) && :physical_periapsis_altitude_m in propertynames(df) || continue
        time_s, altitude_m = _series_xy(df, :time_s, :physical_periapsis_altitude_m)
        isempty(time_s) && continue
        plot!(p, time_s ./ _SECONDS_PER_DAY, altitude_m ./ 1000;
              marker=:circle, markersize=2.5, linewidth=1.5,
              label=_case_label(case))
        plotted = true
    end
    plotted || annotate!(p, 0.5, 0.5, text("No physical periapsis-event rows are available", 10, :gray, :center))
    savefig(p, path)
    return path
end

function _eccentricity_vector_plot(output_root, cases)
    path = joinpath(output_root, "eccentricity_vector_by_pass.png")
    p = plot(; title="Eccentricity-vector coordinates by periapsis pass",
             xlabel="k = e cos(ω)", ylabel="h = e sin(ω)",
             aspect_ratio=:equal, legend=:best, grid=true)
    plotted = false
    using_mean = false
    for case in cases
        df = case.passes
        e_col = _column(df, :mean_e, :e_mean, :e)
        omega_col = _column(df, :mean_omega_rad, :omega_mean_rad, :omega_rad)
        (e_col === nothing || omega_col === nothing) && continue
        e_values, omega_values = _series_xy(df, e_col, omega_col)
        isempty(e_values) && continue
        mean_values = e_col in (:mean_e, :e_mean) && omega_col in (:mean_omega_rad, :omega_mean_rad)
        # The J2 map is a diagnostic coordinate conversion, not validated for the high-e cases.
        source_label = !mean_values ? "osculating" :
            case.name == "vacuum" ? "J2 map (diagnostic)" : "J2 map (diagnostic, not validated at high e)"
        plot!(p, e_values .* cos.(omega_values), e_values .* sin.(omega_values);
              marker=:circle, markersize=2.5, linewidth=1.4,
              label=string(_case_label(case), " (", source_label, ")"))
        plotted = true
    end
    plotted || annotate!(p, 0.5, 0.5, text("No eccentricity and argument-of-periapsis rows are available", 10, :gray, :center))
    savefig(p, path)
    return path
end

function _cumulative_delivered_delta_v_plot(output_root, cases)
    path = joinpath(output_root, "cumulative_delivered_delta_v.png")
    p = plot(; title="Cumulative delivered delta-v",
             xlabel="Elapsed time (days)", ylabel="Cumulative delivered delta-v (m/s)",
             legend=:best, grid=true)
    plotted = false
    for case in cases
        df = case.burns
        if isempty(df)
            # An empty/missing burn table by itself is ambiguous. Draw a zero
            # line only when an actual full-mission summary explicitly records
            # no burns and zero delivered delta-v, and the burn file exists.
            summary = case.summary
            verified_col = _column(summary, :full_mission_verified)
            count_col = _column(summary, :maintenance_burn_count)
            total_dv_col = _column(summary, :delivered_delta_v_mps)
            elapsed_col = _column(summary, :elapsed_horizon_s)
            if case.burns_present && verified_col !== nothing && count_col !== nothing &&
               total_dv_col !== nothing && elapsed_col !== nothing &&
               !any(ismissing, (summary[1, verified_col], summary[1, count_col],
                                summary[1, total_dv_col], summary[1, elapsed_col])) &&
               Bool(summary[1, verified_col]) && Int(summary[1, count_col]) == 0 &&
               Float64(summary[1, total_dv_col]) == 0.0 &&
               isfinite(Float64(summary[1, elapsed_col])) && Float64(summary[1, elapsed_col]) > 0.0
                elapsed_days = Float64(summary[1, elapsed_col]) / _SECONDS_PER_DAY
                plot!(p, [0.0, elapsed_days], [0.0, 0.0]; linewidth=1.6,
                      label=_case_label(case))
                plotted = true
            end
            continue
        end
        dv_col = _column(df, :delivered_delta_v_mps, :actual_delivered_delta_v_mps)
        time_col = _column(df, :execution_stop_time_s, :execution_start_time_s, :execution_time_s)
        (dv_col === nothing || time_col === nothing) && continue
        times, increments = _series_xy(df, time_col, dv_col)
        isempty(times) && continue
        order = sortperm(times)
        times, increments = times[order], abs.(increments[order])
        cumulative = cumsum(increments)
        plot!(p, vcat(0.0, times) ./ _SECONDS_PER_DAY, vcat(0.0, cumulative);
              marker=:circle, markersize=3, linewidth=1.6,
              label=_case_label(case))
        plotted = true
    end
    plotted || annotate!(p, 0.5, 0.5, text("No delivered-burn records are available", 10, :gray, :center))
    savefig(p, path)
    return path
end

function _dimensionless_context_plot(output_root, cases)
    path = joinpath(output_root, "dimensionless_context_by_pass.png")
    p = plot(layout=(3, 1), size=(1000, 900),
             title="Dimensionless context at periapsis by pass", grid=true,
             legend=:outerright)
    panels = ((:Pi_J2, "Π_J2"), (:Pi_J3, "|Π_J3|"), (:Pi_drag, "Π_drag"))
    for (panel_idx, (column, title)) in enumerate(panels)
        plot!(p[panel_idx]; title, xlabel="Elapsed time (days)", ylabel=title, grid=true)
        plotted = false
        for case in cases
            df = case.passes
            :time_s in propertynames(df) && column in propertynames(df) || continue
            times, values = _series_xy(df, :time_s, column)
            isempty(times) && continue
            plot!(p[panel_idx], times ./ _SECONDS_PER_DAY, values;
                  marker=:circle, markersize=2, linewidth=1.3,
                  label=_case_label(case))
            plotted = true
        end
        plotted || annotate!(p[panel_idx], 0.5, 0.5, text("No $title rows available", 9, :gray, :center))
    end
    savefig(p, path)
    return path
end

function _map_vs_physical_plot(output_root, cases)
    path = joinpath(output_root, "map_vs_physical_periapsis_radius.png")
    p = plot(; title="J2 map (diagnostic) minus physical periapsis radius",
             titlefontsize=11, xlabel="Elapsed time (days)", ylabel="Map - physical radius (m)",
             legend=false, grid=true, size=(900, 500), left_margin=6Plots.mm,
             top_margin=3Plots.mm, bottom_margin=4Plots.mm)
    plotted = false
    col = :map_minus_physical_periapsis_radius_m
    for case in cases
        df = case.passes
        :time_s in propertynames(df) && col in propertynames(df) || continue
        times, diffs = _series_xy(df, :time_s, col)
        isempty(times) && continue
        hline!(p, [0.0]; label="", color=:gray, linestyle=:dash)
        plot!(p, times ./ _SECONDS_PER_DAY, diffs; marker=:circle, markersize=2.5,
              linewidth=1.5, label="")
        plotted = true
    end
    plotted || annotate!(p, 0.5, 0.5, text("No map-vs-physical rows are available", 10, :gray, :center))
    savefig(p, path)
    return path
end

"""Generate corridor, eccentricity-vector, delivered-delta-v, and context PNGs from case CSVs."""
function write_corridor_plots(outputroot::AbstractString, config)
    mkpath(outputroot)
    cases = _case_data(outputroot)
    return (
        physical_periapsis=_physical_periapsis_plot(outputroot, cases, config),
        eccentricity_vector=_eccentricity_vector_plot(outputroot, cases),
        cumulative_delivered_delta_v=_cumulative_delivered_delta_v_plot(outputroot, cases),
        dimensionless_context=_dimensionless_context_plot(outputroot, cases),
        map_vs_physical_periapsis=_map_vs_physical_plot(outputroot, cases),
    )
end

end # module PeriapsisCorridorPlots
