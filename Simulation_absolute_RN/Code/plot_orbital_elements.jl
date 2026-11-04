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

const ALTITUDES = Tuple(reverse(campaign_settings("absolute_rn").target_altitudes_km))
const INCLINATIONS = (0.0, 0.5, 1.0)
const COUNTS = campaign_settings("absolute_rn").helper_counts
const POLICIES = ("gve_inc", "gve_sma")
const CASE_COUNT = length(ALTITUDES) * length(INCLINATIONS) * length(COUNTS) * length(POLICIES)
const COLORS = (:blue, :red, :green)
const COMBINED_COLORS = (
    gve_sma=("#652884", "#4E5689", "#018B38"),
    gve_inc=("#B46DA9", "#41B9C1", "#57AF37"),
)
const INPUT = joinpath(DEFAULT_OUTPUT, "paper_plot_mode")
const OUTPUT = joinpath(DEFAULT_OUTPUT, "images", "orbital_elements")
const RTN_FIELDS = Tuple(Symbol(replace(string(name), "target_" => ""))
    for name in TARGET_DIAGNOSTIC_COLUMNS)

function rtn_values(batch, index)
    values = NamedTuple{RTN_FIELDS}(Tuple(Float64(getproperty(batch, name)[index])
        for name in TARGET_DIAGNOSTIC_COLUMNS))
    all(isfinite, values) || error("Nonfinite cumulative RTN diagnostic at batch row $index")
    absolute = NamedTuple{(:dv_R_absolute_mps, :dv_T_absolute_mps, :dv_N_absolute_mps)}(
        Tuple(getproperty(values, Symbol("dv_$(component)_positive_mps")) -
            getproperty(values, Symbol("dv_$(component)_negative_mps")) for component in ("R", "T", "N")))
    merge(values, absolute)
end

function check_rtn(values, previous; context="")
    residuals = ntuple(3) do k
        component = ("R", "T", "N")[k]
        net = getproperty(values, Symbol("dv_$(component)_mps"))
        positive = getproperty(values, Symbol("dv_$(component)_positive_mps"))
        negative = getproperty(values, Symbol("dv_$(component)_negative_mps"))
        last_positive = getproperty(previous, Symbol("dv_$(component)_positive_mps"))
        last_negative = getproperty(previous, Symbol("dv_$(component)_negative_mps"))
        tolerance = component == "T" ? 1e-7 : 1e-8
        positive >= -tolerance && negative <= tolerance ||
            error("$component branches have invalid signs $context")
        abs(net - positive - negative) < 1e-8 ||
            error("$component branches do not sum to the signed net $context")
        residual = max(0.0, last_positive - positive, negative - last_negative)
        residual <= tolerance || error("$component branch reversal exceeds tolerance $context")
        residual
    end
    residuals
end

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
    settings["schema_version"] == SCHEMA_VERSION ||
        error("Expected schema $SCHEMA_VERSION in $directory")
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
    previous_rtn = NamedTuple{RTN_FIELDS}(ntuple(_ -> 0.0, length(RTN_FIELDS)))
    max_reversal = zeros(3)
    row = 0
    for batch in Arrow.Stream(path)
        if include_rtn
            all(name -> hasproperty(batch, name), TARGET_DIAGNOSTIC_COLUMNS) ||
                error("Missing cumulative target RTN diagnostics: $path")
        end
        for local_row in eachindex(batch.time)
            row += 1
            row <= length(times) || error("Extra trajectory rows: $path")
            isapprox(batch.time[local_row], times[row]; atol=32eps(duration), rtol=0) ||
                error("Unexpected trajectory time at row $row: $path")
            position = [getproperty(batch, name)[local_row] for name in pos_columns]
            velocity = [getproperty(batch, name)[local_row] for name in vel_columns]
            if include_rtn
                rtn = rtn_values(batch, local_row)
                max_reversal .= max.(max_reversal,
                    check_rtn(rtn, previous_rtn; context="at row $row in $path"))
                previous_rtn = rtn
            end
            if row == 1
                isapprox(position, initial["position_m"]; atol=1e-8, rtol=0) &&
                    isapprox(velocity, initial["velocity_mps"]; atol=1e-8, rtol=0) ||
                    error("Initial target state mismatch: $path")
                if include_rtn
                    all(iszero, rtn) ||
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
                    samples[row] = merge(samples[row], rtn)
                end
            end
        end
    end
    row == length(times) || error("Incomplete trajectory: $path")
    if any(max_reversal .> (1e-9, 1e-8, 1e-9))
        @warn "Small RTN branch reversals within verifier tolerances; no clamping applied" path max_reversal
    end
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
            for field in keys(expected)
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
    all(value -> isfinite(value) || (field == :draan_deg && isnan(value)), values) ||
        error("Invalid orbital-element values for $field")
    finite = filter(isfinite, values)
    isempty(finite) && error("No finite values for $field")
    lo, hi = extrema(finite)
    if lo == hi
        padding = field == :da_m ? 1.0 : field == :de ? 1e-7 : 1e-5
        lo, hi = lo - padding, hi + padding
    end
    _, lower, upper = Plots.PlotUtils.optimize_ticks(lo, hi; strict_span=false, k_ideal=7)
    (lower, upper)
end

function dv_ylimits(values, field)
    !isempty(values) && all(isfinite, values) ||
        error("Expected nonempty, finite delta-v values for $field")
    lo, hi = extrema(values)
    lo, hi = min(0.0, lo), max(0.0, hi)
    if lo == hi
        name = string(field)
        lo = endswith(name, "_positive_mps") || endswith(name, "_absolute_mps") ? 0.0 : -1e-5
        hi = endswith(name, "_negative_mps") ? 0.0 : 1e-5
    end
    _, lower, upper = Plots.PlotUtils.optimize_ticks(lo, hi; strict_span=false, k_ideal=7)
    (lower, upper)
end

function quantity_spec(quantity)
    quantity == :oe && return (
        fields=(:da_m, :de, :di_deg, :draan_deg),
        labels=("\u0394a (m)", "\u0394e", "\u0394i (deg)", "\u0394\u03a9 (deg)"),
        prefix="OE", note="Osculating changes include J2 + laser; RAAN is unwrapped and omitted for equatorial orbits.")
    quantity == :rtn && return (
        fields=(:dv_R_mps, :dv_T_mps, :dv_N_mps),
        labels=("\u0394v_R (m/s)", "\u0394v_T (m/s)", "\u0394v_N (m/s)"),
        prefix="dv_RTN", note="Cumulative signed laser-only delta-v in instantaneous target RTN; cancellation retained.")
    quantity == :absolute && return (
        fields=(:dv_R_absolute_mps, :dv_T_absolute_mps, :dv_N_absolute_mps),
        labels=("\u0394v_R abs (m/s)", "\u0394v_T abs (m/s)", "\u0394v_N abs (m/s)"),
        prefix="dv_RTN_absolute", note="Accumulated absolute laser delta-v = positive minus negative; NOT the absolute value of the signed net.")
    quantity == :split_n && return (
        fields=(:dv_R_mps, :dv_T_mps, :dv_N_positive_mps, :dv_N_negative_mps),
        labels=("\u0394v_R (m/s)", "\u0394v_T (m/s)", "\u0394v_N+ (m/s)", "\u0394v_N- (m/s)"),
        prefix="dv_RTN_split_N", note="R and T: signed net. N+: integral of max(a_N,0). N-: integral of min(a_N,0), retained negative.")
    quantity == :branches && return (
        fields=(:dv_R_positive_mps, :dv_R_negative_mps, :dv_T_positive_mps,
            :dv_T_negative_mps, :dv_N_positive_mps, :dv_N_negative_mps),
        labels=("\u0394v_R+ (m/s)", "\u0394v_R- (m/s)", "\u0394v_T+ (m/s)",
            "\u0394v_T- (m/s)", "\u0394v_N+ (m/s)", "\u0394v_N- (m/s)"),
        prefix="dv_RTN_branches", note="Laser-only branches: positive = integral of max(a,0); negative = integral of min(a,0). Negative signs retained.")
    throw(ArgumentError("Unsupported plot quantity: $quantity"))
end

function make_plot(rows, thresholds, policy, window, output; quantity=:oe, altitudes=ALTITUDES)
    spec = quantity_spec(quantity)
    !isempty(altitudes) && all(in(ALTITUDES), altitudes) &&
        length(unique(altitudes)) == length(altitudes) ||
        throw(ArgumentError("Plot altitudes must be a nonempty, unique subset of $ALTITUDES"))
    policy in (POLICIES..., "combined") || throw(ArgumentError("Unsupported scheduler: $policy"))
    combined = policy == "combined"
    policies = combined ? ("gve_sma", "gve_inc") : (policy,)
    max_helpers = maximum(COUNTS)
    fields, labels = spec.fields, spec.labels
    n_rows = length(fields)
    # Share each field across all supplied windows, altitudes and schedulers.
    limit_function = quantity == :oe ? ylimits : dv_ylimits
    shared_limits = Dict(field => limit_function([getproperty(row, field) for row in rows], field)
        for field in fields)
    subplots = Plots.Plot[]
    for (r, field) in enumerate(fields), (c, altitude) in enumerate(altitudes)
        selected = filter(row -> row.scheduler in policies && row.orbits == window &&
            row.target_altitude_km == altitude, rows)
        limits = shared_limits[field]
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
            quantity != :oe && !all(isfinite, values) &&
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
        text("Dashed: same geometry for SMA / INC", 9, :left))
    annotate!(footer, 0.5, 0.43,
        text("$policy  |  $window initial target periods  |  L_max = 200 km  |  fixed circular geometry estimate", 11))
    annotate!(footer, 0.5, 0.09, text(spec.note, 10))
    n_cols = length(altitudes)
    layout = combined ? (@layout [grid(n_rows, n_cols); footer{0.16h}]) :
        (@layout [grid(n_rows, n_cols); footer{0.12h}])
    height = (combined ? 1280 : 1200) - 250 * (4 - n_rows)
    figure = plot(subplots..., footer; layout, size=(max(1600, 320 * n_cols), height), dpi=160,
        fontfamily="Computer Modern", left_margin=10Plots.mm, right_margin=3Plots.mm)
    plot!(figure; subplot=length(subplots) + 1, top_margin=10Plots.mm)
    prefix = spec.prefix
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
