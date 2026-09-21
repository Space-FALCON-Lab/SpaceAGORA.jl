const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using CSV
using DataFrames
using LinearAlgebra
using Plots
using Plots.PlotMeasures: mm
using Printf
using Statistics

include(joinpath(REPO_ROOT, "examples", "common.jl"))
setup_gram_example!()

using GRAMSuite

const DEFAULT_INITIAL_TIME = SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
const DEFAULT_PLANETS = ("earth", "mars", "venus", "titan")
const DEFAULT_ALTITUDE_RANGES_KM = Dict(
    "earth" => (0.0, 150.0),
    "mars" => (0.0, 180.0),
    "venus" => (0.0, 180.0),
    "titan" => (0.0, 700.0)
)

function _has_arg(args::Vector{String}, name::String)::Bool
    prefix = name * "="
    return any(arg -> startswith(arg, prefix), args)
end

function _arg_value(args::Vector{String}, name::String, default)
    prefix = name * "="
    for arg in args
        startswith(arg, prefix) && return parse(typeof(default), arg[length(prefix) + 1:end])
    end
    return default
end

function _arg_value(args::Vector{String}, name::String, default::String)
    prefix = name * "="
    for arg in args
        startswith(arg, prefix) && return String(arg[length(prefix) + 1:end])
    end
    return default
end

function _planet_list(args::Vector{String})
    if _has_arg(args, "--planet")
        return (lowercase(_arg_value(args, "--planet", "mars")),)
    end
    raw = _arg_value(args, "--planets", join(DEFAULT_PLANETS, ","))
    return Tuple(lowercase(strip(p)) for p in split(raw, ",") if !isempty(strip(p)))
end

function _altitude_range_km(args::Vector{String}, planet_name::String)
    default_h0, default_h1 = get(DEFAULT_ALTITUDE_RANGES_KM, planet_name, (0.0, 150.0))
    h0_km = _arg_value(args, "--h0-km", default_h0)
    h1_km = _arg_value(args, "--h1-km", default_h1)
    return h0_km, h1_km
end

@inline function _relerr(a::Float64, b::Float64)::Float64
    return abs(a - b) / max(abs(a), eps(Float64))
end

function _linear_samples(;
    n_points::Int,
    h0_km::Float64,
    h1_km::Float64,
    lat0_deg::Float64,
    lat1_deg::Float64,
    lon0_deg::Float64,
    lon1_deg::Float64,
    t0_s::Float64,
    t1_s::Float64
)
    n_points >= 2 || throw(ArgumentError("n_points must be at least 2."))
    denom = n_points - 1
    return [
        (
            elapsed_time_s=t0_s + (k - 1) * (t1_s - t0_s) / denom,
            height_km=h0_km + (k - 1) * (h1_km - h0_km) / denom,
            latitude_deg=lat0_deg + (k - 1) * (lat1_deg - lat0_deg) / denom,
            longitude_deg=lon0_deg + (k - 1) * (lon1_deg - lon0_deg) / denom
        )
        for k in 1:n_points
    ]
end

function _query_samples_with_midpoints(trajectory_samples)
    rows = NamedTuple[]
    for k in eachindex(trajectory_samples)
        s = trajectory_samples[k]
        push!(
            rows,
            merge(
                s,
                (
                    is_trajectory_node=true,
                    trajectory_segment_index=k == length(trajectory_samples) ? max(1, k - 1) : k,
                    trajectory_segment_fraction=k == length(trajectory_samples) ? 1.0 : 0.0
                )
            )
        )
        if k < length(trajectory_samples)
            sn = trajectory_samples[k + 1]
            push!(
                rows,
                (
                    elapsed_time_s=0.5 * (s.elapsed_time_s + sn.elapsed_time_s),
                    height_km=0.5 * (s.height_km + sn.height_km),
                    latitude_deg=0.5 * (s.latitude_deg + sn.latitude_deg),
                    longitude_deg=0.5 * (s.longitude_deg + sn.longitude_deg),
                    is_trajectory_node=false,
                    trajectory_segment_index=k,
                    trajectory_segment_fraction=0.5
                )
            )
        end
    end
    return rows
end

function _point_state_sample(model, sample)
    gram = model.gram
    atmos = model.gram_atmosphere
    set_position! = Base.invokelatest(getproperty, gram, :set_position!)
    update! = Base.invokelatest(getproperty, gram, :update!)
    get_dynamics_state = Base.invokelatest(getproperty, gram, :get_dynamics_state)
    get_winds_state = Base.invokelatest(getproperty, gram, :get_winds_state)

    return lock(RuntimeServices.GRAM_LOCK) do
        Base.invokelatest(
            set_position!,
            atmos;
            height=sample.height_km,
            latitude=sample.latitude_deg,
            longitude=sample.longitude_deg,
            elapsed_time=sample.elapsed_time_s
        )
        err = Base.invokelatest(update!, atmos)
        if err != 0
            get_error_message = Base.invokelatest(getproperty, gram, :get_error_message)
            throw(ErrorException("GRAM update failed (code=$err): $(Base.invokelatest(get_error_message))"))
        end
        dyn = Base.invokelatest(get_dynamics_state, atmos)
        winds = Base.invokelatest(get_winds_state, atmos)
        return (
            rho=Float64(dyn.density),
            temp=Float64(dyn.temperature),
            wind_e=Float64(winds.perturbedEWWind),
            wind_n=Float64(winds.perturbedNSWind),
            wind_u=Float64(winds.perturbedVerticalWind),
            wind_nominal_e=Float64(winds.ewWind),
            wind_nominal_n=Float64(winds.nsWind),
            wind_nominal_u=Float64(winds.verticalWind)
        )
    end
end

function _point_density_rows(model, samples)
    rows = Vector{NamedTuple}(undef, length(samples))
    elapsed_s = @elapsed begin
        for (k, sample) in enumerate(samples)
            rows[k] = _point_state_sample(model, sample)
        end
    end
    return rows, elapsed_s
end

function _trajectory_density_rows(model, samples)
    first_sample = first(samples)
    last_sample = last(samples)
    denom = length(samples) - 1

    trajectory = nothing
    elapsed_s = @elapsed begin
        trajectory_fn = Base.invokelatest(getproperty, model.gram, :generate_trajectory)
        trajectory = Base.invokelatest(
            trajectory_fn,
            model.gram_atmosphere;
            initial_height=first_sample.height_km,
            initial_latitude=first_sample.latitude_deg,
            initial_longitude=first_sample.longitude_deg,
            initial_elapsed_time=first_sample.elapsed_time_s,
            delta_height=(last_sample.height_km - first_sample.height_km) / denom,
            delta_latitude=(last_sample.latitude_deg - first_sample.latitude_deg) / denom,
            delta_longitude=(last_sample.longitude_deg - first_sample.longitude_deg) / denom,
            delta_elapsed_time=(last_sample.elapsed_time_s - first_sample.elapsed_time_s) / denom,
            n_points=length(samples),
            update_initial_perturbations=true
        )
    end

    rows = [
        (
            rho=Float64(pt.dynamics.density),
            temp=Float64(pt.dynamics.temperature),
            wind_e=Float64(pt.winds.perturbedEWWind),
            wind_n=Float64(pt.winds.perturbedNSWind),
            wind_u=Float64(pt.winds.perturbedVerticalWind),
            wind_nominal_e=Float64(pt.winds.ewWind),
            wind_nominal_n=Float64(pt.winds.nsWind),
            wind_nominal_u=Float64(pt.winds.verticalWind)
        )
        for pt in trajectory
    ]
    return rows, elapsed_s
end

const TRAJECTORY_STATE_FIELDS = (
    :rho,
    :temp,
    :wind_e,
    :wind_n,
    :wind_u,
    :wind_nominal_e,
    :wind_nominal_n,
    :wind_nominal_u
)

@inline _lerp_value(a::Float64, b::Float64, x::Float64)::Float64 = a + x * (b - a)

function _interpolate_trajectory_row(trajectory_rows, segment_index::Int, x::Float64)
    idx = clamp(segment_index, 1, length(trajectory_rows) - 1)
    row0 = trajectory_rows[idx]
    row1 = trajectory_rows[idx + 1]
    return NamedTuple{TRAJECTORY_STATE_FIELDS}((
        _lerp_value(Float64(getproperty(row0, field)), Float64(getproperty(row1, field)), x)
        for field in TRAJECTORY_STATE_FIELDS
    ))
end

function _comparison_table(query_samples, point_rows, trajectory_rows)
    rows = NamedTuple[]
    for k in eachindex(query_samples)
        sample = query_samples[k]
        p = point_rows[k]
        tr = _interpolate_trajectory_row(
            trajectory_rows,
            sample.trajectory_segment_index,
            sample.trajectory_segment_fraction
        )
        point_wind = (p.wind_e, p.wind_n, p.wind_u)
        trajectory_wind = (tr.wind_e, tr.wind_n, tr.wind_u)
        point_nominal_wind = (p.wind_nominal_e, p.wind_nominal_n, p.wind_nominal_u)
        trajectory_nominal_wind = (tr.wind_nominal_e, tr.wind_nominal_n, tr.wind_nominal_u)
        wind_abs_err_mps = norm((
            p.wind_e - tr.wind_e,
            p.wind_n - tr.wind_n,
            p.wind_u - tr.wind_u
        ))
        wind_nominal_abs_err_mps = norm((
            p.wind_nominal_e - tr.wind_nominal_e,
            p.wind_nominal_n - tr.wind_nominal_n,
            p.wind_nominal_u - tr.wind_nominal_u
        ))
        push!(
            rows,
            merge(
                (sample_index=k,),
                sample,
                (
                    rho_full=p.rho,
                    rho_trajectory=tr.rho,
                    trajectory_interpolated=!sample.is_trajectory_node,
                    trajectory_segment_index=sample.trajectory_segment_index,
                    trajectory_segment_fraction=sample.trajectory_segment_fraction,
                    rho_abs_err=abs(p.rho - tr.rho),
                    rho_rel_err=_relerr(p.rho, tr.rho),
                    temp_full=p.temp,
                    temp_trajectory=tr.temp,
                    temp_abs_err=abs(p.temp - tr.temp),
                    wind_e_full_mps=p.wind_e,
                    wind_n_full_mps=p.wind_n,
                    wind_u_full_mps=p.wind_u,
                    wind_e_trajectory_mps=tr.wind_e,
                    wind_n_trajectory_mps=tr.wind_n,
                    wind_u_trajectory_mps=tr.wind_u,
                    wind_full_mps=norm(point_wind),
                    wind_trajectory_mps=norm(trajectory_wind),
                    wind_abs_err_mps=wind_abs_err_mps,
                    wind_rel_err=wind_abs_err_mps / max(norm(point_wind), eps(Float64)),
                    wind_nominal_e_full_mps=p.wind_nominal_e,
                    wind_nominal_n_full_mps=p.wind_nominal_n,
                    wind_nominal_u_full_mps=p.wind_nominal_u,
                    wind_nominal_e_trajectory_mps=tr.wind_nominal_e,
                    wind_nominal_n_trajectory_mps=tr.wind_nominal_n,
                    wind_nominal_u_trajectory_mps=tr.wind_nominal_u,
                    wind_nominal_full_mps=norm(point_nominal_wind),
                    wind_nominal_trajectory_mps=norm(trajectory_nominal_wind),
                    wind_nominal_abs_err_mps=wind_nominal_abs_err_mps,
                    wind_nominal_rel_err=wind_nominal_abs_err_mps / max(norm(point_nominal_wind), eps(Float64))
                )
            )
        )
    end
    return DataFrame(rows)
end

function _comparison_table(planet_name::String, query_samples, point_rows, trajectory_rows)
    table = _comparison_table(query_samples, point_rows, trajectory_rows)
    table.planet = fill(planet_name, nrow(table))
    select!(table, :planet, Not(:planet))
    return table
end

function _summary_table(planet_name::String, trajectory_nodes::Int, point_elapsed_s::Float64, trajectory_elapsed_s::Float64, table::DataFrame)
    return DataFrame([(
        planet=planet_name,
        trajectory_nodes=trajectory_nodes,
        comparison_samples=nrow(table),
        interpolated_comparison_samples=count(table.trajectory_interpolated),
        full_point_eval_s=point_elapsed_s,
        trajectory_eval_s=trajectory_elapsed_s,
        speedup=point_elapsed_s / max(trajectory_elapsed_s, 1e-9),
        rho_abs_max=maximum(table.rho_abs_err),
        rho_abs_mean=mean(table.rho_abs_err),
        rho_rel_max=maximum(table.rho_rel_err),
        rho_rel_mean=mean(table.rho_rel_err),
        temp_abs_max=maximum(table.temp_abs_err),
        temp_abs_mean=mean(table.temp_abs_err),
        wind_abs_max_mps=maximum(table.wind_abs_err_mps),
        wind_abs_mean_mps=mean(table.wind_abs_err_mps),
        wind_rel_max=maximum(table.wind_rel_err),
        wind_rel_mean=mean(table.wind_rel_err),
        wind_nominal_abs_max_mps=maximum(table.wind_nominal_abs_err_mps),
        wind_nominal_abs_mean_mps=mean(table.wind_nominal_abs_err_mps),
        wind_nominal_rel_max=maximum(table.wind_nominal_rel_err),
        wind_nominal_rel_mean=mean(table.wind_nominal_rel_err)
    )])
end

@inline function _plot_relerr_floor(x::Float64)::Float64
    return x > 0.0 ? x : eps(Float64)
end

function _error_axis_plot(; xlabel::String, title::String)
    return plot(
        xlabel=xlabel,
        ylabel="Altitude (km)",
        xscale=:log10,
        legend=false,
        grid=true,
        framestyle=:box,
        title=title,
        guidefontsize=11,
        tickfontsize=9,
        titlefontsize=12,
        left_margin=10mm,
        right_margin=6mm,
        bottom_margin=16mm,
        top_margin=7mm
    )
end

function _value_axis_plot(; xlabel::String, title::String)
    return plot(
        xlabel=xlabel,
        ylabel="Altitude (km)",
        legend=false,
        grid=true,
        framestyle=:box,
        title=title,
        guidefontsize=11,
        tickfontsize=9,
        titlefontsize=12,
        left_margin=10mm,
        right_margin=6mm,
        bottom_margin=16mm,
        top_margin=7mm
    )
end

function _save_error_plot(table::DataFrame, plot_path::String)
    rho_plot = _error_axis_plot(
        xlabel="Density relative error",
        title="Density"
    )
    temp_plot = _error_axis_plot(
        xlabel="Temperature error (K)",
        title="Temperature"
    )
    wind_plot = _error_axis_plot(
        xlabel="Wind vector error (m/s)",
        title="Wind Absolute"
    )
    wind_rel_plot = _error_axis_plot(
        xlabel="Perturbed wind relative error",
        title="Perturbed Wind Relative"
    )
    nominal_wind_plot = _error_axis_plot(
        xlabel="Nominal wind vector error (m/s)",
        title="Nominal Wind Absolute"
    )
    nominal_wind_rel_plot = _error_axis_plot(
        xlabel="Nominal wind relative error",
        title="Nominal Wind Relative"
    )

    for planet_name in unique(table.planet)
        rows = sort(table[table.planet .== planet_name, :], :height_km)
        label = uppercasefirst(String(planet_name))
        plot!(
            rho_plot,
            _plot_relerr_floor.(rows.rho_rel_err),
            rows.height_km;
            label=label,
            linewidth=2,
            marker=:circle,
            markersize=3
        )
        plot!(
            temp_plot,
            _plot_relerr_floor.(rows.temp_abs_err),
            rows.height_km;
            label=label,
            linewidth=2,
            marker=:circle,
            markersize=3
        )
        plot!(
            wind_plot,
            _plot_relerr_floor.(rows.wind_abs_err_mps),
            rows.height_km;
            label=label,
            linewidth=2,
            marker=:circle,
            markersize=3
        )
        plot!(
            wind_rel_plot,
            _plot_relerr_floor.(rows.wind_rel_err),
            rows.height_km;
            label=label,
            linewidth=2,
            marker=:circle,
            markersize=3
        )
        plot!(
            nominal_wind_plot,
            _plot_relerr_floor.(rows.wind_nominal_abs_err_mps),
            rows.height_km;
            label=label,
            linewidth=2,
            marker=:circle,
            markersize=3
        )
        plot!(
            nominal_wind_rel_plot,
            _plot_relerr_floor.(rows.wind_nominal_rel_err),
            rows.height_km;
            label=label,
            linewidth=2,
            marker=:circle,
            markersize=3
        )
    end

    plt = plot(
        rho_plot,
        temp_plot,
        wind_plot,
        wind_rel_plot,
        nominal_wind_plot,
        nominal_wind_rel_plot;
        layout=(2, 3),
        size=(1900, 1100),
        dpi=180,
        legend=:outerright,
        plot_title="Full GRAM vs Trajectory GRAM Errors",
        plot_titlefontsize=14,
        margin=8mm
    )

    mkpath(dirname(plot_path))
    savefig(plt, plot_path)
    return plot_path
end

function _plot_wind_pair!(
    plt,
    rows::DataFrame,
    x_full,
    x_trajectory;
    label::String,
    color_idx::Int,
    show_label::Bool
)
    plot!(
        plt,
        x_full,
        rows.height_km;
        label=show_label ? "$label full" : false,
        color=color_idx,
        linewidth=2,
        linestyle=:solid
    )
    plot!(
        plt,
        x_trajectory,
        rows.height_km;
        label=show_label ? "$label trajectory" : false,
        color=color_idx,
        linewidth=2,
        linestyle=:dash
    )
    return plt
end

function _wind_columns(nominal::Bool)
    if nominal
        return (
            speed_full=:wind_nominal_full_mps,
            speed_trajectory=:wind_nominal_trajectory_mps,
            e_full=:wind_nominal_e_full_mps,
            e_trajectory=:wind_nominal_e_trajectory_mps,
            n_full=:wind_nominal_n_full_mps,
            n_trajectory=:wind_nominal_n_trajectory_mps,
            u_full=:wind_nominal_u_full_mps,
            u_trajectory=:wind_nominal_u_trajectory_mps,
            title="True GRAM vs Trajectory GRAM Nominal Wind Values"
        )
    end
    return (
        speed_full=:wind_full_mps,
        speed_trajectory=:wind_trajectory_mps,
        e_full=:wind_e_full_mps,
        e_trajectory=:wind_e_trajectory_mps,
        n_full=:wind_n_full_mps,
        n_trajectory=:wind_n_trajectory_mps,
        u_full=:wind_u_full_mps,
        u_trajectory=:wind_u_trajectory_mps,
        title="True GRAM vs Trajectory GRAM Perturbed Wind Values"
    )
end

function _save_wind_values_plot(table::DataFrame, plot_path::String; nominal::Bool=false)
    cols = _wind_columns(nominal)
    speed_plot = _value_axis_plot(
        xlabel="Wind speed (m/s)",
        title="Speed"
    )
    ew_plot = _value_axis_plot(
        xlabel="EW wind (m/s)",
        title="East-West"
    )
    ns_plot = _value_axis_plot(
        xlabel="NS wind (m/s)",
        title="North-South"
    )
    vertical_plot = _value_axis_plot(
        xlabel="Vertical wind (m/s)",
        title="Vertical"
    )

    for (idx, planet_name) in enumerate(unique(table.planet))
        rows = sort(table[table.planet .== planet_name, :], :height_km)
        label = uppercasefirst(String(planet_name))
        _plot_wind_pair!(
            speed_plot,
            rows,
            rows[!, cols.speed_full],
            rows[!, cols.speed_trajectory];
            label=label,
            color_idx=idx,
            show_label=true
        )
        _plot_wind_pair!(
            ew_plot,
            rows,
            rows[!, cols.e_full],
            rows[!, cols.e_trajectory];
            label=label,
            color_idx=idx,
            show_label=false
        )
        _plot_wind_pair!(
            ns_plot,
            rows,
            rows[!, cols.n_full],
            rows[!, cols.n_trajectory];
            label=label,
            color_idx=idx,
            show_label=false
        )
        _plot_wind_pair!(
            vertical_plot,
            rows,
            rows[!, cols.u_full],
            rows[!, cols.u_trajectory];
            label=label,
            color_idx=idx,
            show_label=false
        )
    end

    plt = plot(
        speed_plot,
        ew_plot,
        ns_plot,
        vertical_plot;
        layout=(2, 2),
        size=(1400, 1100),
        dpi=180,
        legend=:outerright,
        plot_title=cols.title,
        plot_titlefontsize=14,
        margin=8mm
    )

    mkpath(dirname(plot_path))
    savefig(plt, plot_path)
    return plot_path
end

function _run_planet_comparison(
    planet_name::String;
    n_points::Int,
    h0_km::Float64,
    h1_km::Float64,
    lat0_deg::Float64,
    lat1_deg::Float64,
    lon0_deg::Float64,
    lon1_deg::Float64,
    t0_s::Float64,
    t1_s::Float64
)
    density_model = Base.invokelatest(
        SM.GRAMAtmosphereModel;
        planet_name=planet_name,
        initial_time=DEFAULT_INITIAL_TIME
    )
    trajectory_samples = _linear_samples(
        n_points=n_points,
        h0_km=h0_km,
        h1_km=h1_km,
        lat0_deg=lat0_deg,
        lat1_deg=lat1_deg,
        lon0_deg=lon0_deg,
        lon1_deg=lon1_deg,
        t0_s=t0_s,
        t1_s=t1_s
    )
    query_samples = _query_samples_with_midpoints(trajectory_samples)

    point_rows, point_elapsed_s = _point_density_rows(density_model, query_samples)
    trajectory_rows, trajectory_elapsed_s = _trajectory_density_rows(density_model, trajectory_samples)
    table = _comparison_table(planet_name, query_samples, point_rows, trajectory_rows)
    summary = _summary_table(planet_name, length(trajectory_samples), point_elapsed_s, trajectory_elapsed_s, table)
    return summary, table
end

function run_comparison(args::Vector{String}=ARGS)
    planet_names = _planet_list(args)
    n_points = _arg_value(args, "--n", 80)
    lat0_deg = _arg_value(args, "--lat0-deg", -10.0)
    lat1_deg = _arg_value(args, "--lat1-deg", 5.0)
    lon0_deg = _arg_value(args, "--lon0-deg", 20.0)
    lon1_deg = _arg_value(args, "--lon1-deg", 45.0)
    t0_s = _arg_value(args, "--t0-s", 0.0)
    t1_s = _arg_value(args, "--t1-s", 600.0)

    println("Comparing point-by-point GRAM density against GRAM trajectory generation by altitude")
    println(@sprintf(
        "planets=%s trajectory_nodes/planet=%d comparison_points/planet=%d lat=[%.3f, %.3f] deg lon=[%.3f, %.3f] deg t=[%.3f, %.3f] s",
        join(planet_names, ","),
        n_points,
        2 * n_points - 1,
        lat0_deg,
        lat1_deg,
        lon0_deg,
        lon1_deg,
        t0_s,
        t1_s
    ))

    summary_tables = DataFrame[]
    sample_tables = DataFrame[]
    for planet_name in planet_names
        h0_km, h1_km = _altitude_range_km(args, planet_name)
        println(@sprintf("Running %s altitude sweep: %.3f km to %.3f km", planet_name, h0_km, h1_km))
        summary_i, table_i = _run_planet_comparison(
            planet_name;
            n_points=n_points,
            h0_km=h0_km,
            h1_km=h1_km,
            lat0_deg=lat0_deg,
            lat1_deg=lat1_deg,
            lon0_deg=lon0_deg,
            lon1_deg=lon1_deg,
            t0_s=t0_s,
            t1_s=t1_s
        )
        push!(summary_tables, summary_i)
        push!(sample_tables, table_i)
    end

    summary = vcat(summary_tables...)
    table = vcat(sample_tables...)

    mkpath(joinpath(REPO_ROOT, "output"))
    summary_path = joinpath(REPO_ROOT, "output", "gram_full_vs_trajectory_density_summary.csv")
    samples_path = joinpath(REPO_ROOT, "output", "gram_full_vs_trajectory_density_samples.csv")
    plot_path = joinpath(REPO_ROOT, "output", "gram_full_vs_trajectory_density_error_by_altitude.png")
    wind_values_plot_path = joinpath(REPO_ROOT, "output", "gram_full_vs_trajectory_wind_values_by_altitude.png")
    nominal_wind_values_plot_path = joinpath(REPO_ROOT, "output", "gram_full_vs_trajectory_nominal_wind_values_by_altitude.png")
    CSV.write(summary_path, summary)
    CSV.write(samples_path, table)
    _save_error_plot(table, plot_path)
    _save_wind_values_plot(table, wind_values_plot_path)
    _save_wind_values_plot(table, nominal_wind_values_plot_path; nominal=true)

    println("\nSummary:")
    show(summary, allrows=true, allcols=true)
    println("\n\nSaved CSVs:")
    println("  summary: $summary_path")
    println("  samples: $samples_path")
    println("  plot:    $plot_path")
    println("  wind:    $wind_values_plot_path")
    println("  nominal: $nominal_wind_values_plot_path")
    return (
        summary=summary,
        samples=table,
        plot=plot_path,
        wind_plot=wind_values_plot_path,
        nominal_wind_plot=nominal_wind_values_plot_path
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_comparison()
end
