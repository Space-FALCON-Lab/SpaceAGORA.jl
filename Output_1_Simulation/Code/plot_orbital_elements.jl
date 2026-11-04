#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module OrbitalElementPlots

include("StudyConfig.jl")
using .StudyConfig
using Arrow, CSV, LinearAlgebra, Printf, TOML
get!(ENV, "GKSwstype", "100")
using Plots

const ALTITUDES = Tuple(sort(collect((campaign_settings("main").target_altitudes_km...,
    campaign_settings("targets9001100km").target_altitudes_km...)); rev=true))
const INCLINATIONS = (0.0, 0.5, 1.0)
const COUNTS = (campaign_settings("main").helper_counts...,
    campaign_settings("helpers400500").helper_counts...,
    campaign_settings("helpers600800").helper_counts...)
const POLICIES = ("gve_inc", "gve_sma")
const CASE_COUNT = length(ALTITUDES) * length(INCLINATIONS) * length(COUNTS) * length(POLICIES)
const COLORS = (:blue, :red, :green)
const COMBINED_COLORS = (
    gve_sma=("#652884", "#4E5689", "#018B38"),
    gve_inc=("#B46DA9", "#41B9C1", "#57AF37"),
)
const INPUT = joinpath(DEFAULT_OUTPUT, "paper_plot_mode")
const OUTPUT = joinpath(DEFAULT_OUTPUT, "images", "orbital_elements")

"""Range-only estimate for equally spaced helpers in one fixed circular plane."""
function continuous_access(rh, rt, range_m, delta_i)
    all(isfinite, (rh, rt, range_m, delta_i)) ||
        throw(ArgumentError("Access geometry must be finite"))
    rh > 0 && rt > 0 && range_m > 0 ||
        throw(ArgumentError("Radii and range must be positive"))
    0 <= delta_i < pi / 2 ||
        throw(ArgumentError("The restricted geometry requires 0 <= delta_i < pi/2"))
    numerator = rh^2 + rt^2 - range_m^2
    numerator > 0 ||
        throw(ArgumentError("The restricted geometry requires rh^2 + rt^2 - range^2 > 0"))
    worst_range = sqrt((rh - rt)^2 + 4rh * rt * sin(delta_i / 2)^2)
    range_m > worst_range ||
        return (alpha_min_rad=0.0, continuous_helpers=Inf, minimum_helpers=missing)
    alpha = acos(clamp(numerator / (2rh * rt * cos(delta_i)), -1.0, 1.0))
    alpha > 0 || error("Accessible arc is below floating-point resolution")
    population = pi / alpha
    (alpha_min_rad=alpha, continuous_helpers=population,
        minimum_helpers=ceil(Int, population))
end

"""Osculating Keplerian elements; equatorial RAAN is undefined (NaN)."""
function orbital_elements(position, velocity, mu)
    all(isfinite, position) && all(isfinite, velocity) && isfinite(mu) && mu > 0 ||
        throw(ArgumentError("Invalid Cartesian state or gravitational parameter"))
    radius = norm(position)
    radius > 0 || throw(ArgumentError("Position must be nonzero"))
    h = cross(position, velocity)
    norm(h) > 0 || throw(ArgumentError("Angular momentum must be nonzero"))
    energy = dot(velocity, velocity) / 2 - mu / radius
    energy < 0 || throw(ArgumentError("Expected a bound orbit"))
    a = -mu / (2energy)
    e = norm(cross(velocity, h) / mu - position / radius)
    inclination = atan(hypot(h[1], h[2]), h[3])
    raan = hypot(h[1], h[2]) <= 1e-12 * norm(h) ? NaN : atan(h[1], -h[2])
    (a=a, e=e, i=inclination, raan=raan)
end

signed_angle(angle) = atan(sin(angle), cos(angle))

function exact_sample_index(times, requested_time)
    index = searchsortedfirst(times, requested_time)
    tolerance = 32eps(max(abs(requested_time), 1.0))
    for candidate in (index, index - 1)
        if 1 <= candidate <= length(times) &&
                abs(times[candidate] - requested_time) <= tolerance
            return candidate
        end
    end
    error("No saved trajectory row at t=$requested_time s; interpolation is intentionally disabled")
end

function read_case(directory, info, windows; include_rtn=false)
    settings = info["settings"]
    info["status"] == "verified" || error("Unverified simulation: $directory")
    duration = settings["duration_s"]
    times = collect(range(0.0, duration; length=info["trajectory_rows"]))
    period = duration / settings["initial_target_periods"]
    indices = [exact_sample_index(times, window * period) for window in windows]
    path = joinpath(directory, "trajectory.feather")
    filehash(path) == info["checksums"]["trajectory.feather"] ||
        error("Trajectory checksum mismatch: $path")
    target = settings["target_index"]
    initial = info["initial_states"][target]
    initial["index"] == target && initial["id"] == settings["target_id"] ||
        error("Target metadata mismatch: $path")
    oe0 = orbital_elements(initial["position_m"], initial["velocity_mps"], settings["mu_m3_s2"])
    pos_columns = ntuple(axis -> Symbol("sc$(target)_pos_$axis"), 3)
    vel_columns = ntuple(axis -> Symbol("sc$(target)_vel_$axis"), 3)
    samples = Dict{Int,NamedTuple}()
    previous_raan = oe0.raan
    delta_raan = isfinite(oe0.raan) ? 0.0 : NaN
    row = 0
    for batch in Arrow.Stream(path)
        if include_rtn
            all(name -> hasproperty(batch, name),
                (:target_dv_R_mps, :target_dv_T_mps, :target_dv_N_mps)) ||
                error("Missing cumulative target RTN diagnostics: $path")
        end
        for local_row in eachindex(batch.time)
            row += 1
            row <= length(times) || error("Extra trajectory rows: $path")
            isapprox(batch.time[local_row], times[row]; atol=32eps(duration), rtol=0) ||
                error("Unexpected trajectory time at row $row: $path")
            position = [getproperty(batch, name)[local_row] for name in pos_columns]
            velocity = [getproperty(batch, name)[local_row] for name in vel_columns]
            if row == 1
                isapprox(position, initial["position_m"]; atol=1e-8, rtol=0) &&
                    isapprox(velocity, initial["velocity_mps"]; atol=1e-8, rtol=0) ||
                    error("Initial target state mismatch: $path")
                if include_rtn
                    all(value -> isfinite(value) && abs(value) <= 1e-12,
                        (batch.target_dv_R_mps[1], batch.target_dv_T_mps[1], batch.target_dv_N_mps[1])) ||
                        error("Cumulative target RTN diagnostics must start at zero: $path")
                end
            end
            oe = orbital_elements(position, velocity, settings["mu_m3_s2"])
            # Unwrap every saved row, not just the two requested endpoints.
            delta_raan += signed_angle(oe.raan - previous_raan)
            previous_raan = oe.raan
            if row in indices
                samples[row] = (time_s=Float64(batch.time[local_row]), row=row,
                    da_m=oe.a - oe0.a, de=oe.e - oe0.e,
                    di_deg=rad2deg(oe.i - oe0.i), draan_deg=rad2deg(delta_raan),
                    raan_defined=isfinite(delta_raan))
                if include_rtn
                    rtn = (dv_R_mps=Float64(batch.target_dv_R_mps[local_row]),
                        dv_T_mps=Float64(batch.target_dv_T_mps[local_row]),
                        dv_N_mps=Float64(batch.target_dv_N_mps[local_row]))
                    all(isfinite, rtn) || error("Nonfinite target RTN endpoint at row $row: $path")
                    samples[row] = merge(samples[row], rtn)
                end
            end
        end
    end
    row == length(times) || error("Incomplete trajectory: $path")
    [merge((orbits=window,), samples[index]) for (window, index) in zip(windows, indices)]
end

function load_grid(input, windows; include_rtn=false)
    isdir(input) || error("Input directory does not exist: $input")
    cases = Dict{Tuple{Int,Float64,Int,String},Tuple{String,Dict{String,Any}}}()
    excluded_cases = 0
    for (directory, _, files) in walkdir(input)
        "run_info.toml" in files || continue
        info = TOML.parsefile(joinpath(directory, "run_info.toml"))
        s = info["settings"]
        if s["helper_altitude_km"] != 1000 || s["target_altitude_km"] ∉ ALTITUDES ||
                s["helper_count"] ∉ COUNTS
            excluded_cases += 1
            continue
        end
        key = (Int(s["target_altitude_km"]), Float64(s["target_inclination_deg"]),
            Int(s["helper_count"]), String(s["scheduler"]))
        haskey(cases, key) && error("Ambiguous duplicate simulation $key in $input")
        cases[key] = (directory, info)
    end
    excluded_cases > 0 && @info "Plotting the $CASE_COUNT-case grid through $(maximum(COUNTS)) helpers at 1000 km; other cases excluded" cases=excluded_cases
    length(cases) == CASE_COUNT ||
        error("Expected the complete $CASE_COUNT-case study, found $(length(cases))")
    rows = NamedTuple[]
    thresholds = NamedTuple[]
    for altitude in ALTITUDES, inclination in INCLINATIONS
        push!(thresholds, merge((target_altitude_km=altitude, delta_i_deg=inclination,
            helper_radius_m=EARTH_RADIUS_M + 1000e3,
            target_radius_m=EARTH_RADIUS_M + altitude * 1000, range_m=RANGE_M),
            continuous_access(EARTH_RADIUS_M + 1000e3,
                EARTH_RADIUS_M + altitude * 1000, RANGE_M, deg2rad(inclination))))
        for policy in POLICIES, helpers in COUNTS
            key = (altitude, inclination, helpers, policy)
            haskey(cases, key) || error("Missing simulation: $key")
            directory, info = cases[key]
            expected = case_settings(StudyCase(altitude, inclination, helpers, Symbol(policy)))
            for field in ("target_altitude_km", "helper_altitude_km", "target_inclination_deg",
                    "helper_inclination_deg", "helper_count", "target_index", "target_id",
                    "eccentricity", "RAAN_deg", "argument_of_periapsis_deg",
                    "target_true_anomaly_deg", "helper_phase_rule", "initial_target_periods",
                    "duration_s", "trajectory_rows", "mu_m3_s2", "earth_equatorial_radius_m",
                    "range_m", "scheduler", "J2", "laser_power_W", "mass_kg")
                info["settings"][field] == expected[field] ||
                    error("Unexpected $field in $directory")
            end
            info["trajectory_rows"] == expected["trajectory_rows"] ||
                error("Unexpected trajectory row count: $directory")
            for result in read_case(directory, info, windows; include_rtn)
                push!(rows, merge((scheduler=policy, target_altitude_km=altitude,
                    delta_i_deg=inclination, helpers=helpers), result,
                    (source_trajectory=joinpath(directory, "trajectory.feather"),)))
            end
        end
        println("Read target altitude $altitude km, inclination $inclination deg")
    end
    rows, thresholds
end

function ylimits(values, field)
    finite = filter(isfinite, values)
    isempty(finite) && error("No finite values for $field")
    lo, hi = extrema(finite)
    minimum_span = field == :da_m ? 1.0 : field == :de ? 1e-7 : 1e-5
    padding = max((hi - lo) * 0.15, minimum_span)
    (lo - padding, hi + padding)
end

function make_plot(rows, thresholds, policy, window, output; quantity=:oe)
    quantity in (:oe, :rtn) || throw(ArgumentError("Unsupported plot quantity: $quantity"))
    policy in (POLICIES..., "combined") || throw(ArgumentError("Unsupported scheduler: $policy"))
    combined = policy == "combined"
    policies = combined ? ("gve_sma", "gve_inc") : (policy,)
    max_helpers = maximum(COUNTS)
    fields = quantity == :oe ? (:da_m, :de, :di_deg, :draan_deg) :
        (:dv_R_mps, :dv_T_mps, :dv_N_mps)
    labels = quantity == :oe ?
        ("\u0394a (m)", "\u0394e", "\u0394i (deg)", "\u0394\u03a9 (deg)") :
        ("\u0394v_R (m/s)", "\u0394v_T (m/s)", "\u0394v_N (m/s)")
    n_rows = length(fields)
    subplots = Plots.Plot[]
    for (r, field) in enumerate(fields), (c, altitude) in enumerate(ALTITUDES)
        selected = filter(row -> row.scheduler in policies && row.orbits == window &&
            row.target_altitude_km == altitude, rows)
        limits = ylimits([getproperty(row, field) for row in selected], field)
        panel = plot(; legend=false, grid=true, gridalpha=0.18, framestyle=:axes,
            title=r == 1 ? "h_t = $altitude km" : "", titlefontsize=12,
            ylabel=c == 1 ? labels[r] : "", xlabel=r == n_rows ? "N_helpers" : "",
            xlims=(-8, max_helpers + 15),
            xticks=filter(n -> n == 1 || n % 100 == 0, collect(COUNTS)),
            xrotation=60, ylims=limits,
            tickfontsize=10, guidefontsize=12, tick_direction=:out,
            left_margin=3Plots.mm, right_margin=2Plots.mm,
            top_margin=2Plots.mm, bottom_margin=3Plots.mm)
        for scheduler in policies, (k, inclination) in enumerate(INCLINATIONS)
            curve = sort(filter(row -> row.scheduler == scheduler &&
                row.delta_i_deg == inclination, selected); by=row -> row.helpers)
            [row.helpers for row in curve] == collect(COUNTS) ||
                error("Incomplete or duplicate curve: $scheduler, $altitude km, $inclination deg, $window orbits")
            values = [getproperty(row, field) for row in curve]
            quantity == :rtn && !all(isfinite, values) &&
                error("Nonfinite cumulative RTN values: $scheduler, $altitude km, $inclination deg")
            color = combined ? getproperty(COMBINED_COLORS, Symbol(scheduler))[k] : COLORS[k]
            any(isfinite, values) && plot!(panel, [row.helpers for row in curve], values;
                color, linestyle=combined && scheduler == "gve_inc" ? :dot : :solid,
                lw=1.6, marker=:circle, ms=3, markerstrokewidth=0)
        end
        for (k, inclination) in enumerate(INCLINATIONS)
            color = combined ? COMBINED_COLORS.gve_sma[k] : COLORS[k]
            threshold = only(filter(t -> t.target_altitude_km == altitude &&
                t.delta_i_deg == inclination, thresholds))
            population = threshold.minimum_helpers
            if !ismissing(population) && population <= max_helpers
                vline!(panel, [population]; color, linestyle=:dash, lw=1.2, alpha=0.75)
            elseif r == 1
                label = ismissing(population) ? "No finite N_cont" : "N_cont = $population > $max_helpers"
                annotate!(panel, max_helpers + 5, limits[2] - 0.045 * (limits[2] - limits[1]),
                    text("$inclination deg: $label", 9, color, :right))
            end
        end
        push!(subplots, panel)
    end
    footer = plot(; framestyle=:none, ticks=nothing, legend=false, grid=false,
        xlims=(0, 1), ylims=(0, 1), margin=0Plots.mm)
    for (j, scheduler) in enumerate(policies), (k, inclination) in enumerate(INCLINATIONS)
        x = 0.03 + (k - 1) * 0.22
        y = combined ? 0.90 - (j - 1) * 0.24 : 0.80
        color = combined ? getproperty(COMBINED_COLORS, Symbol(scheduler))[k] : COLORS[k]
        label = combined ? "$(scheduler == "gve_sma" ? "SMA" : "INC"), i_initial = $inclination deg" :
            "\u0394i_initial = $inclination deg"
        plot!(footer, [x, x + 0.035], [y, y]; color, lw=1.6,
            linestyle=combined && scheduler == "gve_inc" ? :dot : :solid,
            marker=:circle, ms=3, markerstrokewidth=0)
        annotate!(footer, x + 0.045, y, text(label, 11, :left))
    end
    plot!(footer, [0.73, 0.77], [0.80, 0.80]; color=:black, linestyle=:dash, lw=1.2)
    annotate!(footer, 0.78, 0.80, text("N_cont = ceil(\u03c0 / \u03b1_min)", 11, :left))
    combined && annotate!(footer, 0.73, 0.60,
        text("Dashed: inclination colors from SMA; applies to both schedulers", 9, :left))
    annotate!(footer, 0.5, 0.43,
        text("$policy  |  $window initial target periods  |  L_max = 200 km  |  fixed circular geometry estimate", 11))
    note = quantity == :oe ?
        "Osculating changes include J2 + laser; RAAN is unwrapped and omitted for equatorial orbits." :
        "Cumulative signed laser-only delta-v in the instantaneous target RTN frame (m/s); no absolute-value accumulation."
    annotate!(footer, 0.5, 0.09, text(note, 10))
    n_cols = length(ALTITUDES)
    layout = combined ? (@layout [grid(n_rows, n_cols); footer{0.16h}]) :
        (@layout [grid(n_rows, n_cols); footer{0.12h}])
    height = (combined ? 1280 : 1200) - 250 * (4 - n_rows)
    figure = plot(subplots..., footer; layout, size=(320 * n_cols, height), dpi=160,
        fontfamily="Computer Modern", left_margin=10Plots.mm, right_margin=3Plots.mm)
    plot!(figure; subplot=length(subplots) + 1, top_margin=10Plots.mm)
    prefix = quantity == :oe ? "OE" : "dv_RTN"
    stem = "$(prefix)_from_feather_$(policy)_e0.0000_nu1.0000_$(window)orbits_h_cols"
    for extension in ("png", "pdf")
        path = joinpath(output, "$stem.$extension")
        savefig(figure, path)
        println("Saved $path")
    end
    figure
end

function parse_options(args)
    input, output, windows = INPUT, OUTPUT, [80, 100]
    length(args) % 2 == 0 ||
        error("Usage: plot_orbital_elements.jl [--input DIR] [--output DIR] [--orbits 80,100]")
    for index in 1:2:length(args)
        flag, value = args[index:index+1]
        if flag == "--input"
            input = abspath(value)
        elseif flag == "--output"
            output = abspath(value)
        elseif flag == "--orbits"
            windows = parse.(Int, split(value, ","))
            !isempty(windows) && all(w -> 1 <= w <= 100, windows) &&
                length(unique(windows)) == length(windows) ||
                error("Orbit windows must be unique integers between 1 and 100")
        else
            error("Unknown option: $flag")
        end
    end
    (; input, output, windows)
end

function main(args=ARGS)
    options = parse_options(args)
    rows, thresholds = load_grid(options.input, options.windows)
    mkpath(options.output)
    for window in options.windows
        CSV.write(joinpath(options.output, "orbital_element_deltas_$(window)orbits.csv"),
            filter(row -> row.orbits == window, rows))
    end
    CSV.write(joinpath(options.output, "continuous_access_thresholds.csv"), thresholds)
    gr()
    for policy in POLICIES, window in options.windows
        make_plot(rows, thresholds, policy, window, options.output)
    end
    println("Validated $(length(rows)) endpoints from $CASE_COUNT verified, checksum-matched trajectories.")
    nothing
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    OrbitalElementPlots.main()
end
