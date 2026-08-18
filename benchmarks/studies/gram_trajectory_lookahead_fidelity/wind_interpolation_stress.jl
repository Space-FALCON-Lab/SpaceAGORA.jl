const WIND_STRESS_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(@__DIR__, "main.jl"))

using CSV
using DataFrames
using LinearAlgebra
using Plots
using Plots.PlotMeasures: mm
using Printf
using Statistics

const WIND_STRESS_POINTS = (4, 5, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128)
const WIND_STRESS_HORIZON_S = 900.0
const WIND_STRESS_DEFAULT_OUTDIR = joinpath(WIND_STRESS_REPO_ROOT, "output", "gram_wind_interpolation_stress")

function _wind_stress_usage()
    return """
    Usage:
      julia --project=. benchmarks/studies/gram_trajectory_lookahead_fidelity/wind_interpolation_stress.jl [options]

    Options:
      --planets mars,earth,venus,titan
      --planet mars
      --points 4,5,6,8,12,16,24,32,48,64,96,128
      --methods linear,pchip,cubic,akima
      --horizon-s 900
      --out-dir output/gram_wind_interpolation_stress
      --scenarios diagonal,shear_lat,shear_lon,reverse_alt,time_offset
    """
end

function _stress_parse_cli(args::Vector{String})
    if any(arg -> arg in ("-h", "--help"), args)
        println(_wind_stress_usage())
        exit(0)
    end
    return _parse_cli(args)
end

function _stress_point_counts(opts::Dict{String, String})
    raw = get(opts, "points", join(WIND_STRESS_POINTS, ","))
    points = Tuple(parse(Int, strip(p)) for p in split(raw, ",") if !isempty(strip(p)))
    isempty(points) && throw(ArgumentError("At least one point count is required."))
    any(<(2), points) && throw(ArgumentError("All point counts must be at least 2."))
    return points
end

function _scenario_list(opts::Dict{String, String})
    raw = get(opts, "scenarios", "diagonal,shear_lat,shear_lon,reverse_alt,time_offset")
    scenarios = Tuple(lowercase(strip(s)) for s in split(raw, ",") if !isempty(strip(s)))
    supported = Set(("diagonal", "shear_lat", "shear_lon", "reverse_alt", "time_offset"))
    unsupported = [s for s in scenarios if !(s in supported)]
    isempty(unsupported) || throw(ArgumentError("Unsupported scenario(s): $(join(unsupported, ", "))."))
    return scenarios
end

function _planet_alt_bounds_km(planet_name::String)
    planet_name == "earth" && return 90.0, 150.0
    planet_name == "mars" && return 85.0, 155.0
    planet_name == "venus" && return 105.0, 180.0
    planet_name == "titan" && return 300.0, 800.0
    return 90.0, 160.0
end

function _wind_stress_case(planet_name::String, scenario::String)
    h_low, h_high = _planet_alt_bounds_km(planet_name)
    if scenario == "diagonal"
        return (
            scenario=scenario,
            h0_km=h_high,
            h1_km=h_low,
            lat0_deg=-55.0,
            lat1_deg=55.0,
            lon0_deg=-130.0,
            lon1_deg=140.0,
            t0_s=0.0,
        )
    elseif scenario == "shear_lat"
        return (
            scenario=scenario,
            h0_km=0.55 * h_high + 0.45 * h_low,
            h1_km=0.45 * h_high + 0.55 * h_low,
            lat0_deg=-75.0,
            lat1_deg=75.0,
            lon0_deg=25.0,
            lon1_deg=35.0,
            t0_s=900.0,
        )
    elseif scenario == "shear_lon"
        return (
            scenario=scenario,
            h0_km=0.70 * h_high + 0.30 * h_low,
            h1_km=0.35 * h_high + 0.65 * h_low,
            lat0_deg=10.0,
            lat1_deg=18.0,
            lon0_deg=-170.0,
            lon1_deg=170.0,
            t0_s=1800.0,
        )
    elseif scenario == "reverse_alt"
        return (
            scenario=scenario,
            h0_km=h_low,
            h1_km=h_high,
            lat0_deg=42.0,
            lat1_deg=-42.0,
            lon0_deg=80.0,
            lon1_deg=-95.0,
            t0_s=2700.0,
        )
    elseif scenario == "time_offset"
        return (
            scenario=scenario,
            h0_km=h_high,
            h1_km=h_low,
            lat0_deg=5.0,
            lat1_deg=-35.0,
            lon0_deg=15.0,
            lon1_deg=210.0,
            t0_s=7200.0,
        )
    end
    throw(ArgumentError("Unknown wind stress scenario '$scenario'."))
end

function _stress_samples(case; horizon_s::Float64)
    stop_s = floor(Int, horizon_s)
    samples = [
        begin
            x = Float64(t_s) / horizon_s
            (
                elapsed_time_s=case.t0_s + Float64(t_s),
                height_km=_lerp(case.h0_km, case.h1_km, x),
                latitude_deg=_lerp(case.lat0_deg, case.lat1_deg, x),
                longitude_deg=_lerp(case.lon0_deg, case.lon1_deg, x),
            )
        end
        for t_s in 0:stop_s
    ]
    if samples[end].elapsed_time_s < case.t0_s + horizon_s
        push!(
            samples,
            (
                elapsed_time_s=case.t0_s + horizon_s,
                height_km=case.h1_km,
                latitude_deg=case.lat1_deg,
                longitude_deg=case.lon1_deg,
            )
        )
    end
    return samples
end

function _trajectory_states_for_case(model, case; n_points::Int, horizon_s::Float64)
    denom = n_points - 1
    trajectory_fn = Base.invokelatest(getproperty, model.gram, :generate_trajectory)
    trajectory = Base.invokelatest(
        trajectory_fn,
        model.gram_atmosphere;
        initial_height=case.h0_km,
        initial_latitude=case.lat0_deg,
        initial_longitude=case.lon0_deg,
        initial_elapsed_time=case.t0_s,
        delta_height=(case.h1_km - case.h0_km) / denom,
        delta_latitude=(case.lat1_deg - case.lat0_deg) / denom,
        delta_longitude=(case.lon1_deg - case.lon0_deg) / denom,
        delta_elapsed_time=horizon_s / denom,
        n_points=n_points,
        update_initial_perturbations=true,
    )
    return [
        (
            rho=Float64(pt.dynamics.density),
            temp=Float64(pt.dynamics.temperature),
            wind_e_mps=Float64(pt.winds.perturbedEWWind),
            wind_n_mps=Float64(pt.winds.perturbedNSWind),
            wind_u_mps=Float64(pt.winds.perturbedVerticalWind),
        )
        for pt in trajectory
    ]
end

function _finite_quantile(values::Vector{Float64}, p::Float64)
    xs = sort!(filter(isfinite, values))
    isempty(xs) && return NaN
    idx = clamp(Int(ceil(p * length(xs))), 1, length(xs))
    return xs[idx]
end

function _wind_error_summary(truth_states, interp, samples, case, horizon_s::Float64)
    abs_err = Float64[]
    rel_err = Float64[]
    ew_abs = Float64[]
    ns_abs = Float64[]
    up_abs = Float64[]
    sizehint!(abs_err, length(samples))
    sizehint!(rel_err, length(samples))

    @inbounds for idx in eachindex(samples)
        truth = truth_states[idx]
        x = (samples[idx].elapsed_time_s - case.t0_s) / horizon_s
        approx, _, _ = _interpolated_trajectory_state(interp, x)
        components = (
            truth.wind_e_mps - approx.wind_e_mps,
            truth.wind_n_mps - approx.wind_n_mps,
            truth.wind_u_mps - approx.wind_u_mps,
        )
        err = norm(components)
        truth_norm = norm((truth.wind_e_mps, truth.wind_n_mps, truth.wind_u_mps))
        if isfinite(err)
            push!(abs_err, err)
            push!(rel_err, err / max(truth_norm, eps(Float64)))
            push!(ew_abs, abs(components[1]))
            push!(ns_abs, abs(components[2]))
            push!(up_abs, abs(components[3]))
        end
    end

    return (
        samples=length(samples),
        finite_wind_samples=length(abs_err),
        wind_abs_mean_mps=isempty(abs_err) ? NaN : mean(abs_err),
        wind_abs_p95_mps=_finite_quantile(abs_err, 0.95),
        wind_abs_max_mps=isempty(abs_err) ? NaN : maximum(abs_err),
        wind_rel_mean=isempty(rel_err) ? NaN : mean(rel_err),
        wind_rel_p95=_finite_quantile(rel_err, 0.95),
        wind_rel_max=isempty(rel_err) ? NaN : maximum(rel_err),
        wind_e_abs_p95_mps=_finite_quantile(ew_abs, 0.95),
        wind_n_abs_p95_mps=_finite_quantile(ns_abs, 0.95),
        wind_u_abs_p95_mps=_finite_quantile(up_abs, 0.95),
    )
end

function _run_wind_stress_planet(planet_name::String, opts::Dict{String, String})
    point_counts = _stress_point_counts(opts)
    methods = _interp_methods(opts)
    scenarios = _scenario_list(opts)
    horizon_s = _get_float(opts, "horizon-s", WIND_STRESS_HORIZON_S)

    model = Base.invokelatest(
        SM.GRAMAtmosphereModel;
        planet_name=planet_name,
        initial_time=DEFAULT_INITIAL_TIME,
    )

    rows = NamedTuple[]
    for scenario in scenarios
        case = _wind_stress_case(planet_name, scenario)
        samples = _stress_samples(case; horizon_s=horizon_s)
        truth_states, direct_eval_s = _direct_gram_states(model, samples)

        for n_points in point_counts
            trajectory_states = nothing
            trajectory_eval_s = @elapsed begin
                trajectory_states = _trajectory_states_for_case(model, case; n_points=n_points, horizon_s=horizon_s)
            end

            for method in methods
                interp = _make_interpolator(trajectory_states, method)
                interpolation_eval_s = @elapsed summary = _wind_error_summary(truth_states, interp, samples, case, horizon_s)
                push!(
                    rows,
                    merge(
                        (
                            planet=planet_name,
                            scenario=scenario,
                            interpolation_method=method,
                            trajectory_points=n_points,
                            horizon_s=horizon_s,
                            h0_km=case.h0_km,
                            h1_km=case.h1_km,
                            lat0_deg=case.lat0_deg,
                            lat1_deg=case.lat1_deg,
                            lon0_deg=case.lon0_deg,
                            lon1_deg=case.lon1_deg,
                            t0_s=case.t0_s,
                            direct_eval_s=direct_eval_s,
                            trajectory_eval_s=trajectory_eval_s,
                            interpolation_eval_s=interpolation_eval_s,
                        ),
                        summary
                    )
                )
            end
        end
    end
    return DataFrame(rows)
end

function _wind_metric_plot(; ylabel::String, title::String)
    return plot(
        xlabel="Trajectory interpolation points",
        ylabel=ylabel,
        xscale=:log10,
        yscale=:log10,
        legend=false,
        grid=true,
        framestyle=:box,
        title=title,
        guidefontsize=11,
        tickfontsize=9,
        titlefontsize=12,
        left_margin=10mm,
        bottom_margin=14mm,
        top_margin=7mm,
    )
end

function _save_wind_stress_plot(df::DataFrame, plot_path::String)
    p95_plot = _wind_metric_plot(ylabel="Wind vector p95 error (m/s)", title="P95")
    max_plot = _wind_metric_plot(ylabel="Wind vector max error (m/s)", title="Max")
    rel_plot = _wind_metric_plot(ylabel="Wind vector p95 relative error", title="P95 Relative")

    grouped = combine(
        groupby(df, [:planet, :interpolation_method, :trajectory_points]),
        :wind_abs_p95_mps => mean => :wind_abs_p95_mps,
        :wind_abs_max_mps => maximum => :wind_abs_max_mps,
        :wind_rel_p95 => mean => :wind_rel_p95,
    )

    for (idx, planet_name) in enumerate([p for p in DEFAULT_PLANETS if p in unique(grouped.planet)])
        for method in DEFAULT_INTERP_METHODS
            rows = sort(grouped[(grouped.planet .== planet_name) .& (grouped.interpolation_method .== method), :], :trajectory_points)
            isempty(rows) && continue
            label = "$(uppercasefirst(String(planet_name))) $(method)"
            linestyle = _method_linestyle(method)
            plot!(p95_plot, rows.trajectory_points, _plot_floor.(rows.wind_abs_p95_mps); label=label, color=idx, linestyle=linestyle, linewidth=2)
            plot!(max_plot, rows.trajectory_points, _plot_floor.(rows.wind_abs_max_mps); label=label, color=idx, linestyle=linestyle, linewidth=2)
            plot!(rel_plot, rows.trajectory_points, _plot_floor.(rows.wind_rel_p95); label=label, color=idx, linestyle=linestyle, linewidth=2)
        end
    end

    plt = plot(
        p95_plot,
        max_plot,
        rel_plot;
        layout=(1, 3),
        size=(1800, 560),
        dpi=180,
        legend=:outerright,
        plot_title="GRAM Wind Interpolation Stress Test",
        plot_titlefontsize=14,
        margin=8mm,
    )
    mkpath(dirname(plot_path))
    savefig(plt, plot_path)
    return plot_path
end

function _ranking_table(df::DataFrame)
    winners = combine(groupby(df, [:planet, :scenario, :trajectory_points])) do g
        finite = filter(row -> isfinite(row.wind_abs_p95_mps), g)
        if nrow(finite) == 0
            return (interpolation_method="none", wind_abs_p95_mps=NaN)
        end
        row = finite[argmin(finite.wind_abs_p95_mps), :]
        return (
            interpolation_method=row.interpolation_method,
            wind_abs_p95_mps=row.wind_abs_p95_mps,
        )
    end
    return combine(groupby(winners, [:planet, :interpolation_method]), nrow => :wins)
end

function run_wind_interpolation_stress(args::Vector{String}=ARGS)
    opts = _stress_parse_cli(args)
    planets = _planet_list(opts)
    out_dir = abspath(_get(opts, "out-dir", WIND_STRESS_DEFAULT_OUTDIR))

    println("GRAM wind interpolation stress test")
    println("planets=$(join(planets, ","))")
    println("points=$(join(_stress_point_counts(opts), ","))")
    println("methods=$(join(_interp_methods(opts), ","))")
    println("scenarios=$(join(_scenario_list(opts), ","))")
    println(@sprintf("sample_dt_s=1.000 horizon_s=%.3f", _get_float(opts, "horizon-s", WIND_STRESS_HORIZON_S)))

    frames = DataFrame[]
    for planet_name in planets
        println("Running $(planet_name)...")
        push!(frames, _run_wind_stress_planet(planet_name, opts))
    end
    df = vcat(frames...)

    mkpath(out_dir)
    csv_path = joinpath(out_dir, "gram_wind_interpolation_stress.csv")
    pdf_path = joinpath(out_dir, "gram_wind_interpolation_stress.pdf")
    png_path = joinpath(out_dir, "gram_wind_interpolation_stress.png")
    ranking_path = joinpath(out_dir, "gram_wind_interpolation_stress_rankings.csv")

    CSV.write(csv_path, df)
    _save_wind_stress_plot(df, pdf_path)
    _save_wind_stress_plot(df, png_path)
    rankings = _ranking_table(df)
    CSV.write(ranking_path, rankings)

    println("\nP95 winner counts:")
    show(rankings, allrows=true, allcols=true)
    println("\n\nSaved:")
    println("  csv: $csv_path")
    println("  rankings: $ranking_path")
    println("  pdf: $pdf_path")
    println("  png: $png_path")
    return (results=df, rankings=rankings)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_wind_interpolation_stress()
end
