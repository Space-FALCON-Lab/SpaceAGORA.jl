const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

using CSV
using DataFrames
using LinearAlgebra
using Plots
using Plots.PlotMeasures: mm
using Printf

include(joinpath(REPO_ROOT, "examples", "common.jl"))
setup_gram_example!()

gr()

const DEFAULT_INITIAL_TIME = SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
const DEFAULT_PLANETS = ("mars", "earth", "venus", "titan")
const DEFAULT_POINTS = (3, 4, 5, 6, 8, 12, 16, 24, 32, 48, 64, 128, 256, 512, 1024)
const DEFAULT_INTERP_METHODS = ("linear", "pchip", "cubic", "akima")
const DEFAULT_HORIZON_S = 900.0
const DEFAULT_ALTITUDE_RANGES_KM = Dict(
    "earth" => (130.0, 90.0),
    "mars" => (145.0, 90.0),
    "venus" => (170.0, 110.0),
    "titan" => (700.0, 350.0),
)

const STATE_FIELDS = (
    :rho,
    :temp,
    :wind_e_mps,
    :wind_n_mps,
    :wind_u_mps,
)

function _usage()
    return """
    Usage:
      julia --project=. benchmarks/studies/gram_trajectory_lookahead_fidelity/main.jl [options]

    Options:
      --planets mars,earth,venus,titan
      --planet mars
      --points 3,4,5,6,8,12,16,24,32,48,64,128,256,512,1024
      --methods linear,pchip,cubic,akima
      --horizon-s 900
      --out-dir output/gram_trajectory_lookahead_fidelity
      --h0-km 130 --h1-km 90
      --lat0-deg -15 --lat1-deg 15
      --lon0-deg 20 --lon1-deg 45
    """
end

function _parse_cli(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            println(_usage())
            exit(0)
        end
        startswith(arg, "--") || throw(ArgumentError("Unsupported argument '$arg'. Use --key=value or --key value."))
        body = arg[3:end]
        if occursin("=", body)
            key, value = split(body, "=", limit=2)
            opts[key] = value
        else
            key = body
            if i == length(args) || startswith(args[i + 1], "--")
                opts[key] = "true"
            else
                i += 1
                opts[key] = args[i]
            end
        end
        i += 1
    end
    return opts
end

function _get(opts::Dict{String, String}, key::String, default::String)::String
    return get(opts, key, default)
end

function _get_float(opts::Dict{String, String}, key::String, default::Float64)::Float64
    return parse(Float64, get(opts, key, string(default)))
end

function _planet_list(opts::Dict{String, String})
    if haskey(opts, "planet")
        return (lowercase(strip(opts["planet"])),)
    end
    raw = get(opts, "planets", join(DEFAULT_PLANETS, ","))
    planets = Tuple(lowercase(strip(p)) for p in split(raw, ",") if !isempty(strip(p)))
    isempty(planets) && throw(ArgumentError("At least one planet is required."))
    return planets
end

function _point_counts(opts::Dict{String, String})
    raw = get(opts, "points", join(DEFAULT_POINTS, ","))
    points = Tuple(parse(Int, strip(p)) for p in split(raw, ",") if !isempty(strip(p)))
    isempty(points) && throw(ArgumentError("At least one interpolation point count is required."))
    any(<(2), points) && throw(ArgumentError("All point counts must be at least 2."))
    return points
end

function _interp_methods(opts::Dict{String, String})
    raw = get(opts, "methods", join(DEFAULT_INTERP_METHODS, ","))
    methods = Tuple(lowercase(strip(m)) for m in split(raw, ",") if !isempty(strip(m)))
    isempty(methods) && throw(ArgumentError("At least one interpolation method is required."))
    supported = Set(DEFAULT_INTERP_METHODS)
    unsupported = [m for m in methods if !(m in supported)]
    isempty(unsupported) || throw(ArgumentError("Unsupported interpolation method(s): $(join(unsupported, ", ")). Supported: $(join(DEFAULT_INTERP_METHODS, ", "))."))
    return methods
end

function _altitude_range_km(opts::Dict{String, String}, planet_name::String)
    default_h0, default_h1 = get(DEFAULT_ALTITUDE_RANGES_KM, planet_name, (150.0, 90.0))
    h0_km = _get_float(opts, "h0-km", default_h0)
    h1_km = _get_float(opts, "h1-km", default_h1)
    return h0_km, h1_km
end

@inline function _lerp(a::Float64, b::Float64, x::Float64)::Float64
    return a + x * (b - a)
end

function _sample_at_fraction(;
    x::Float64,
    h0_km::Float64,
    h1_km::Float64,
    lat0_deg::Float64,
    lat1_deg::Float64,
    lon0_deg::Float64,
    lon1_deg::Float64,
    horizon_s::Float64,
)
    return (
        elapsed_time_s=x * horizon_s,
        height_km=_lerp(h0_km, h1_km, x),
        latitude_deg=_lerp(lat0_deg, lat1_deg, x),
        longitude_deg=_lerp(lon0_deg, lon1_deg, x),
    )
end

function _second_samples(;
    h0_km::Float64,
    h1_km::Float64,
    lat0_deg::Float64,
    lat1_deg::Float64,
    lon0_deg::Float64,
    lon1_deg::Float64,
    horizon_s::Float64,
)
    horizon_s >= 1.0 || throw(ArgumentError("--horizon-s must be at least 1 second."))
    stop_s = floor(Int, horizon_s)
    samples = [
        _sample_at_fraction(
            x=Float64(t_s) / horizon_s,
            h0_km=h0_km,
            h1_km=h1_km,
            lat0_deg=lat0_deg,
            lat1_deg=lat1_deg,
            lon0_deg=lon0_deg,
            lon1_deg=lon1_deg,
            horizon_s=horizon_s,
        )
        for t_s in 0:stop_s
    ]
    if samples[end].elapsed_time_s < horizon_s
        push!(
            samples,
            _sample_at_fraction(
                x=1.0,
                h0_km=h0_km,
                h1_km=h1_km,
                lat0_deg=lat0_deg,
                lat1_deg=lat1_deg,
                lon0_deg=lon0_deg,
                lon1_deg=lon1_deg,
                horizon_s=horizon_s,
            )
        )
    end
    return samples
end

function _direct_gram_state(model, sample)
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
            elapsed_time=sample.elapsed_time_s,
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
            wind_e_mps=Float64(winds.perturbedEWWind),
            wind_n_mps=Float64(winds.perturbedNSWind),
            wind_u_mps=Float64(winds.perturbedVerticalWind),
        )
    end
end

function _direct_gram_states(model, samples)
    state_type = NamedTuple{STATE_FIELDS, Tuple{Float64, Float64, Float64, Float64, Float64}}
    states = Vector{state_type}(undef, length(samples))
    elapsed_s = @elapsed begin
        for (idx, sample) in enumerate(samples)
            states[idx] = _direct_gram_state(model, sample)
        end
    end
    return states, elapsed_s
end

function _trajectory_states(model; n_points::Int, h0_km::Float64, h1_km::Float64, lat0_deg::Float64, lat1_deg::Float64, lon0_deg::Float64, lon1_deg::Float64, horizon_s::Float64)
    denom = n_points - 1
    trajectory_fn = Base.invokelatest(getproperty, model.gram, :generate_trajectory)
    trajectory = Base.invokelatest(
        trajectory_fn,
        model.gram_atmosphere;
        initial_height=h0_km,
        initial_latitude=lat0_deg,
        initial_longitude=lon0_deg,
        initial_elapsed_time=0.0,
        delta_height=(h1_km - h0_km) / denom,
        delta_latitude=(lat1_deg - lat0_deg) / denom,
        delta_longitude=(lon1_deg - lon0_deg) / denom,
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

function _interpolated_trajectory_state(states, query_fraction::Float64)
    scaled = query_fraction * (length(states) - 1)
    idx = clamp(floor(Int, scaled) + 1, 1, length(states) - 1)
    x = clamp(scaled - (idx - 1), 0.0, 1.0)
    row0 = states[idx]
    row1 = states[idx + 1]
    state = NamedTuple{STATE_FIELDS}((
        _lerp(Float64(getproperty(row0, field)), Float64(getproperty(row1, field)), x)
        for field in STATE_FIELDS
    ))
    return state, idx, x
end

@inline function _hermite_value(y0::Float64, y1::Float64, m0::Float64, m1::Float64, x::Float64)::Float64
    x2 = x * x
    x3 = x2 * x
    return (2x3 - 3x2 + 1) * y0 +
           (x3 - 2x2 + x) * m0 +
           (-2x3 + 3x2) * y1 +
           (x3 - x2) * m1
end

function _uniform_secants(values::AbstractVector{<:NamedTuple}, field::Symbol)
    n = length(values)
    return [Float64(getproperty(values[i + 1], field)) - Float64(getproperty(values[i], field)) for i in 1:(n - 1)]
end

@inline function _linear_field(values, field::Symbol, idx::Int, x::Float64)::Float64
    return _lerp(Float64(getproperty(values[idx], field)), Float64(getproperty(values[idx + 1], field)), x)
end

function _cubic_slopes_uniform(d::Vector{Float64})
    n = length(d) + 1
    m = Vector{Float64}(undef, n)
    n == 2 && return [d[1], d[1]]
    m[1] = d[1]
    @inbounds for i in 2:(n - 1)
        m[i] = 0.5 * (d[i - 1] + d[i])
    end
    m[n] = d[end]
    return m
end

@inline function _pchip_endpoint_slope(d1::Float64, d2::Float64)::Float64
    m = 0.5 * (3.0 * d1 - d2)
    if sign(m) != sign(d1)
        return 0.0
    elseif sign(d1) != sign(d2) && abs(m) > abs(3.0 * d1)
        return 3.0 * d1
    end
    return m
end

function _pchip_slopes_uniform(d::Vector{Float64})
    n = length(d) + 1
    m = Vector{Float64}(undef, n)
    if n == 2
        m[1] = d[1]
        m[2] = d[1]
        return m
    end

    m[1] = _pchip_endpoint_slope(d[1], d[2])
    @inbounds for i in 2:(n - 1)
        d0 = d[i - 1]
        d1 = d[i]
        if d0 == 0.0 || d1 == 0.0 || sign(d0) != sign(d1)
            m[i] = 0.0
        else
            m[i] = 2.0 * d0 * d1 / (d0 + d1)
        end
    end
    m[n] = _pchip_endpoint_slope(d[end], d[end - 1])
    return m
end

function _akima_slopes_uniform(d::Vector{Float64})
    n = length(d) + 1
    n <= 4 && return _cubic_slopes_uniform(d)

    ext = Vector{Float64}(undef, n + 4)
    # ext[k + 3] stores the interval slope d[k] for k=1:n-1, plus two
    # linearly extrapolated slopes on each side.
    @inbounds begin
        for k in 1:(n - 1)
            ext[k + 3] = d[k]
        end
        ext[3] = 2.0 * ext[4] - ext[5]
        ext[2] = 2.0 * ext[3] - ext[4]
        ext[1] = 2.0 * ext[2] - ext[3]
        ext[n + 3] = 2.0 * ext[n + 2] - ext[n + 1]
        ext[n + 4] = 2.0 * ext[n + 3] - ext[n + 2]
    end

    m = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        # Akima derivative at node i from surrounding interval slopes:
        # d[i-2], d[i-1], d[i], d[i+1], using extrapolated ends.
        dm2 = ext[i + 1]
        dm1 = ext[i + 2]
        d0 = ext[i + 3]
        d1 = ext[i + 4]
        w1 = abs(d1 - d0)
        w2 = abs(dm1 - dm2)
        denom = w1 + w2
        m[i] = denom == 0.0 ? 0.5 * (dm1 + d0) : (w1 * dm1 + w2 * d0) / denom
    end
    return m
end

function _hermite_slopes(values, field::Symbol, method::String)
    d = _uniform_secants(values, field)
    return if method == "pchip"
        _pchip_slopes_uniform(d)
    elseif method == "akima"
        _akima_slopes_uniform(d)
    elseif method == "cubic"
        _cubic_slopes_uniform(d)
    else
        throw(ArgumentError("Unsupported Hermite interpolation method '$method'."))
    end
end

function _make_interpolator(values, method::String)
    slopes = method == "linear" ? nothing : Dict(field => _hermite_slopes(values, field, method) for field in STATE_FIELDS)
    return (values=values, method=method, slopes=slopes)
end

function _hermite_field(interp, field::Symbol, idx::Int, x::Float64)::Float64
    values = interp.values
    slopes = interp.slopes[field]
    y0 = Float64(getproperty(values[idx], field))
    y1 = Float64(getproperty(values[idx + 1], field))
    return _hermite_value(y0, y1, slopes[idx], slopes[idx + 1], x)
end

function _interpolated_trajectory_state(interp::NamedTuple, query_fraction::Float64)
    values = interp.values
    if interp.method == "linear"
        return _interpolated_trajectory_state(values, query_fraction)
    end
    scaled = query_fraction * (length(values) - 1)
    idx = clamp(floor(Int, scaled) + 1, 1, length(values) - 1)
    x = clamp(scaled - (idx - 1), 0.0, 1.0)
    state = NamedTuple{STATE_FIELDS}((
        _hermite_field(interp, field, idx, x)
        for field in STATE_FIELDS
    ))
    return state, idx, x
end

@inline function _relerr(truth::Float64, approx::Float64)::Float64
    return abs(truth - approx) / max(abs(truth), eps(Float64))
end

function _mean_error_summary(truth_states, interp, samples, horizon_s::Float64)
    n = length(samples)
    rho_abs_sum = 0.0
    rho_rel_sum = 0.0
    temp_abs_sum = 0.0
    temp_rel_sum = 0.0
    wind_abs_sum = 0.0
    wind_rel_sum = 0.0
    rho_count = 0
    temp_count = 0
    wind_count = 0
    rho_abs_max = 0.0
    rho_rel_max = 0.0
    temp_abs_max = 0.0
    temp_rel_max = 0.0
    wind_abs_max = 0.0
    wind_rel_max = 0.0

    @inbounds for idx in 1:n
        truth = truth_states[idx]
        sample = samples[idx]
        approx, _, _ = _interpolated_trajectory_state(
            interp,
            sample.elapsed_time_s / horizon_s
        )
        truth_wind_norm = norm((truth.wind_e_mps, truth.wind_n_mps, truth.wind_u_mps))
        wind_abs_err = norm((
            truth.wind_e_mps - approx.wind_e_mps,
            truth.wind_n_mps - approx.wind_n_mps,
            truth.wind_u_mps - approx.wind_u_mps,
        ))
        rho_abs_err = abs(truth.rho - approx.rho)
        rho_rel_err = _relerr(truth.rho, approx.rho)
        temp_abs_err = abs(truth.temp - approx.temp)
        temp_rel_err = _relerr(truth.temp, approx.temp)
        wind_rel_err = wind_abs_err / max(truth_wind_norm, eps(Float64))

        if isfinite(rho_abs_err) && isfinite(rho_rel_err)
            rho_count += 1
            rho_abs_sum += rho_abs_err
            rho_rel_sum += rho_rel_err
            rho_abs_max = max(rho_abs_max, rho_abs_err)
            rho_rel_max = max(rho_rel_max, rho_rel_err)
        end
        if isfinite(temp_abs_err) && isfinite(temp_rel_err)
            temp_count += 1
            temp_abs_sum += temp_abs_err
            temp_rel_sum += temp_rel_err
            temp_abs_max = max(temp_abs_max, temp_abs_err)
            temp_rel_max = max(temp_rel_max, temp_rel_err)
        end
        if isfinite(wind_abs_err) && isfinite(wind_rel_err)
            wind_count += 1
            wind_abs_sum += wind_abs_err
            wind_rel_sum += wind_rel_err
            wind_abs_max = max(wind_abs_max, wind_abs_err)
            wind_rel_max = max(wind_rel_max, wind_rel_err)
        end
    end

    return (
        samples=n,
        rho_samples=rho_count,
        temp_samples=temp_count,
        wind_samples=wind_count,
        rho_abs_mean=rho_count == 0 ? NaN : rho_abs_sum / rho_count,
        rho_rel_mean=rho_count == 0 ? NaN : rho_rel_sum / rho_count,
        temp_abs_mean_k=temp_count == 0 ? NaN : temp_abs_sum / temp_count,
        temp_rel_mean=temp_count == 0 ? NaN : temp_rel_sum / temp_count,
        wind_abs_mean_mps=wind_count == 0 ? NaN : wind_abs_sum / wind_count,
        wind_rel_mean=wind_count == 0 ? NaN : wind_rel_sum / wind_count,
        rho_abs_max=rho_count == 0 ? NaN : rho_abs_max,
        rho_rel_max=rho_count == 0 ? NaN : rho_rel_max,
        temp_abs_max_k=temp_count == 0 ? NaN : temp_abs_max,
        temp_rel_max=temp_count == 0 ? NaN : temp_rel_max,
        wind_abs_max_mps=wind_count == 0 ? NaN : wind_abs_max,
        wind_rel_max=wind_count == 0 ? NaN : wind_rel_max,
    )
end

function _run_planet(planet_name::String, opts::Dict{String, String})
    point_counts = _point_counts(opts)
    interp_methods = _interp_methods(opts)
    horizon_s = _get_float(opts, "horizon-s", DEFAULT_HORIZON_S)

    h0_km, h1_km = _altitude_range_km(opts, planet_name)
    lat0_deg = _get_float(opts, "lat0-deg", -15.0)
    lat1_deg = _get_float(opts, "lat1-deg", 15.0)
    lon0_deg = _get_float(opts, "lon0-deg", 20.0)
    lon1_deg = _get_float(opts, "lon1-deg", 45.0)
    samples = _second_samples(
        h0_km=h0_km,
        h1_km=h1_km,
        lat0_deg=lat0_deg,
        lat1_deg=lat1_deg,
        lon0_deg=lon0_deg,
        lon1_deg=lon1_deg,
        horizon_s=horizon_s,
    )

    model = Base.invokelatest(
        SM.GRAMAtmosphereModel;
        planet_name=planet_name,
        initial_time=DEFAULT_INITIAL_TIME,
    )
    truth_states, truth_elapsed_s = _direct_gram_states(model, samples)

    rows = NamedTuple[]
    for n_points in point_counts
        states = nothing
        trajectory_elapsed_s = @elapsed begin
            states = _trajectory_states(
                model;
                n_points=n_points,
                h0_km=h0_km,
                h1_km=h1_km,
                lat0_deg=lat0_deg,
                lat1_deg=lat1_deg,
                lon0_deg=lon0_deg,
                lon1_deg=lon1_deg,
                horizon_s=horizon_s,
            )
        end
        for method in interp_methods
            interp = _make_interpolator(states, method)
            interp_elapsed_s = @elapsed summary = _mean_error_summary(truth_states, interp, samples, horizon_s)
            push!(
                rows,
                merge(
                    (
                        planet=planet_name,
                        interpolation_method=method,
                        trajectory_points=n_points,
                        horizon_s=horizon_s,
                        h0_km=h0_km,
                        h1_km=h1_km,
                        lat0_deg=lat0_deg,
                        lat1_deg=lat1_deg,
                        lon0_deg=lon0_deg,
                        lon1_deg=lon1_deg,
                        direct_eval_s=truth_elapsed_s,
                        trajectory_eval_s=trajectory_elapsed_s,
                        interpolation_eval_s=interp_elapsed_s,
                    ),
                    summary
                )
            )
        end
    end
    return DataFrame(rows)
end

@inline function _plot_floor(x::Float64)::Float64
    return x > 0.0 ? x : eps(Float64)
end

function _method_linestyle(method::AbstractString)
    method == "linear" && return :solid
    method == "pchip" && return :dash
    method == "cubic" && return :dot
    method == "akima" && return :dashdot
    return :solid
end

function _metric_plot(; ylabel::String, title::String)
    return plot(
        xlabel="Trajectory interpolation points",
        ylabel=ylabel,
        yscale=:log10,
        xscale=:log10,
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

function _save_absolute_plot(df::DataFrame, plot_path::String)
    rho_plot = _metric_plot(ylabel="Mean density absolute error (kg/m^3)", title="Density")
    temp_plot = _metric_plot(ylabel="Mean temperature absolute error (K)", title="Temperature")
    wind_plot = _metric_plot(ylabel="Mean wind vector error (m/s)", title="Wind")

    planets = [p for p in DEFAULT_PLANETS if p in unique(df.planet)]
    for (idx, planet_name) in enumerate(planets)
        for method in DEFAULT_INTERP_METHODS
            rows = sort(df[(df.planet .== planet_name) .& (df.interpolation_method .== method), :], :trajectory_points)
            isempty(rows) && continue
            label = "$(uppercasefirst(String(planet_name))) $(method)"
            linestyle = _method_linestyle(method)
            plot!(
                rho_plot,
                rows.trajectory_points,
                _plot_floor.(rows.rho_abs_mean);
                label=label,
                linewidth=2,
                linestyle=linestyle,
                marker=:circle,
                markersize=4,
                color=idx,
            )
            plot!(
                temp_plot,
                rows.trajectory_points,
                _plot_floor.(rows.temp_abs_mean_k);
                label=label,
                linewidth=2,
                linestyle=linestyle,
                marker=:circle,
                markersize=4,
                color=idx,
            )
            plot!(
                wind_plot,
                rows.trajectory_points,
                _plot_floor.(rows.wind_abs_mean_mps);
                label=label,
                linewidth=2,
                linestyle=linestyle,
                marker=:circle,
                markersize=4,
                color=idx,
            )
        end
    end

    horizon_s = first(df.horizon_s)
    plt = plot(
        wind_plot,
        rho_plot,
        temp_plot;
        layout=(1, 3),
        size=(1800, 560),
        dpi=180,
        legend=:outerright,
        plot_title=@sprintf("Mean Absolute GRAM Trajectory Fidelity by Interpolation Method, 1 s Samples (%.1f s Horizon)", horizon_s),
        plot_titlefontsize=14,
        margin=8mm,
    )

    mkpath(dirname(plot_path))
    savefig(plt, plot_path)
    return plot_path
end

function _save_relative_plot(df::DataFrame, plot_path::String)
    rho_plot = _metric_plot(ylabel="Mean density relative error", title="Density")
    temp_plot = _metric_plot(ylabel="Mean temperature relative error", title="Temperature")
    wind_plot = _metric_plot(ylabel="Mean wind vector relative error", title="Wind")

    planets = [p for p in DEFAULT_PLANETS if p in unique(df.planet)]
    for (idx, planet_name) in enumerate(planets)
        for method in DEFAULT_INTERP_METHODS
            rows = sort(df[(df.planet .== planet_name) .& (df.interpolation_method .== method), :], :trajectory_points)
            isempty(rows) && continue
            label = "$(uppercasefirst(String(planet_name))) $(method)"
            linestyle = _method_linestyle(method)
            plot!(
                rho_plot,
                rows.trajectory_points,
                _plot_floor.(rows.rho_rel_mean);
                label=label,
                linewidth=2,
                linestyle=linestyle,
                marker=:circle,
                markersize=4,
                color=idx,
            )
            plot!(
                temp_plot,
                rows.trajectory_points,
                _plot_floor.(rows.temp_rel_mean);
                label=label,
                linewidth=2,
                linestyle=linestyle,
                marker=:circle,
                markersize=4,
                color=idx,
            )
            plot!(
                wind_plot,
                rows.trajectory_points,
                _plot_floor.(rows.wind_rel_mean);
                label=label,
                linewidth=2,
                linestyle=linestyle,
                marker=:circle,
                markersize=4,
                color=idx,
            )
        end
    end

    horizon_s = first(df.horizon_s)
    plt = plot(
        wind_plot,
        rho_plot,
        temp_plot;
        layout=(1, 3),
        size=(1800, 560),
        dpi=180,
        legend=:outerright,
        plot_title=@sprintf("Mean Relative GRAM Trajectory Fidelity by Interpolation Method, 1 s Samples (%.1f s Horizon)", horizon_s),
        plot_titlefontsize=14,
        margin=8mm,
    )

    mkpath(dirname(plot_path))
    savefig(plt, plot_path)
    return plot_path
end

function run_study(args::Vector{String}=ARGS)
    opts = _parse_cli(args)
    planets = _planet_list(opts)
    out_dir = abspath(_get(opts, "out-dir", joinpath(REPO_ROOT, "output", "gram_trajectory_lookahead_fidelity")))

    println("GRAM trajectory lookahead fidelity study")
    println("planets=$(join(planets, ",")) points=$(join(_point_counts(opts), ","))")
    println("methods=$(join(_interp_methods(opts), ","))")
    println(@sprintf(
        "sample_dt_s=1.000 horizon_s=%.3f",
        _get_float(opts, "horizon-s", DEFAULT_HORIZON_S)
    ))

    frames = DataFrame[]
    for planet_name in planets
        println("Running $(planet_name)...")
        push!(frames, _run_planet(planet_name, opts))
    end
    df = vcat(frames...)

    mkpath(out_dir)
    csv_path = joinpath(out_dir, "gram_trajectory_lookahead_fidelity.csv")
    pdf_path = joinpath(out_dir, "gram_trajectory_lookahead_fidelity.pdf")
    png_path = joinpath(out_dir, "gram_trajectory_lookahead_fidelity.png")
    rel_pdf_path = joinpath(out_dir, "gram_trajectory_lookahead_fidelity_relative.pdf")
    rel_png_path = joinpath(out_dir, "gram_trajectory_lookahead_fidelity_relative.png")
    CSV.write(csv_path, df)
    _save_absolute_plot(df, pdf_path)
    _save_absolute_plot(df, png_path)
    _save_relative_plot(df, rel_pdf_path)
    _save_relative_plot(df, rel_png_path)

    println("\nSummary:")
    show(select(df, :planet, :interpolation_method, :trajectory_points, :samples, :wind_abs_mean_mps, :wind_rel_mean, :rho_abs_mean, :rho_rel_mean, :temp_abs_mean_k, :temp_rel_mean), allrows=true, allcols=true)
    println("\n\nSaved:")
    println("  csv: $csv_path")
    println("  pdf: $pdf_path")
    println("  png: $png_path")
    println("  relative_pdf: $rel_pdf_path")
    println("  relative_png: $rel_png_path")
    return df
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_study()
end
