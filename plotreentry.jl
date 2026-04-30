using CSV
using DataFrames
using Dates
using LinearAlgebra
using Plots
using Printf
using Statistics

const REPO_ROOT = @__DIR__

const MIN_DRAG_SIM_CSV = joinpath(REPO_ROOT, "output", "simulation_results.csv")
const MEDIUM_DRAG_SIM_CSV = joinpath(REPO_ROOT, "output", "medium_drag_config", "simulation_results.csv")
const MAX_DRAG_SIM_CSV = joinpath(REPO_ROOT, "output", "max_drag_config", "simulation_results.csv")
const TLE_FILE = joinpath(REPO_ROOT, "UHF_TM_notebook", "tle_opssat.txt")
const OUT_PNG = joinpath(REPO_ROOT, "output", "ops_sat", "ops_sat_altitude_comparison.png")
const MEDIUM_OUT_PNG = joinpath(REPO_ROOT, "output", "ops_sat", "ops_sat_medium_drag_altitude_comparison.png")

const EARTH_RADIUS_KM = 6378.137
const EARTH_MU_M3_S2 = 3.986004418e14

function tle_epoch(line1::AbstractString)
    year = 2000 + parse(Int, line1[19:20])
    day_of_year = parse(Float64, strip(line1[21:32]))
    whole_days = floor(Int, day_of_year) - 1
    frac_day = day_of_year - floor(day_of_year)
    millis = round(Int, frac_day * 86_400_000)
    return DateTime(year, 1, 1) + Day(whole_days) + Millisecond(millis)
end

function tle_mean_altitude_km(line2::AbstractString)
    mean_motion_rev_day = parse(Float64, strip(line2[53:63]))
    n = mean_motion_rev_day * 2π / 86_400.0
    a_m = cbrt(EARTH_MU_M3_S2 / n^2)
    return a_m / 1000.0 - EARTH_RADIUS_KM
end

function tle_elements(line1::AbstractString, line2::AbstractString)
    mean_motion_rev_day = parse(Float64, strip(line2[53:63]))
    n = mean_motion_rev_day * 2π / 86_400.0
    a_m = cbrt(EARTH_MU_M3_S2 / n^2)
    e = parse(Float64, "0." * strip(line2[27:33]))
    mean_anomaly = deg2rad(parse(Float64, strip(line2[44:51])))

    eccentric_anomaly = mean_anomaly
    for _ in 1:20
        eccentric_anomaly -= (
            eccentric_anomaly - e * sin(eccentric_anomaly) - mean_anomaly
        ) / (1.0 - e * cos(eccentric_anomaly))
    end

    true_anomaly = mod(rad2deg(2.0 * atan(
        sqrt(1.0 + e) * sin(eccentric_anomaly / 2.0),
        sqrt(1.0 - e) * cos(eccentric_anomaly / 2.0),
    )), 360.0)

    return (;
        epoch=tle_epoch(line1),
        a=a_m,
        e=e,
        i=parse(Float64, strip(line2[9:16])),
        Ω=parse(Float64, strip(line2[18:25])),
        ω=parse(Float64, strip(line2[35:42])),
        ν=true_anomaly,
    )
end

function load_tle_history(path::String)
    lines = [strip(x) for x in readlines(path) if !isempty(strip(x))]
    timestamps = DateTime[]
    altitudes_km = Float64[]
    elements = NamedTuple[]

    for i in 1:2:length(lines)
        i + 1 <= length(lines) || break
        line1 = lines[i]
        line2 = lines[i + 1]
        startswith(line1, "1 ") || continue
        startswith(line2, "2 ") || continue

        push!(timestamps, tle_epoch(line1))
        push!(altitudes_km, tle_mean_altitude_km(line2))
        push!(elements, tle_elements(line1, line2))
    end

    perm = sortperm(timestamps)
    timestamps = timestamps[perm]
    altitudes_km = altitudes_km[perm]
    elements = elements[perm]

    t0 = first(timestamps)
    days = Dates.value.(timestamps .- t0) ./ 86_400_000.0

    return (; timestamps, days, altitudes_km, elements)
end

function sim_mean_altitude_km(df::DataFrame, idx::Integer)
    x = df.sc1_pos_1[idx]
    y = df.sc1_pos_2[idx]
    z = df.sc1_pos_3[idx]
    vx = df.sc1_vel_1[idx]
    vy = df.sc1_vel_2[idx]
    vz = df.sc1_vel_3[idx]

    r_m = sqrt(x^2 + y^2 + z^2)
    v2_m2_s2 = vx^2 + vy^2 + vz^2

    specific_energy = 0.5 * v2_m2_s2 - EARTH_MU_M3_S2 / r_m
    a_m = -EARTH_MU_M3_S2 / (2.0 * specific_energy)
    return a_m / 1000.0 - EARTH_RADIUS_KM
end

function load_sim_history(path::String, tle_days::AbstractVector; include_final_point::Bool=false)
    df = CSV.read(path, DataFrame)
    sim_days = Float64.(df.time) ./ 86_400.0

    sampled_days = Float64[]
    sampled_altitude_km = Float64[]

    for t_day in tle_days
        Float64(t_day) > last(sim_days) && continue

        idx = argmin(abs.(sim_days .- Float64(t_day)))

        push!(sampled_days, sim_days[idx])
        push!(sampled_altitude_km, sim_mean_altitude_km(df, idx))
    end

    if include_final_point && (isempty(sampled_days) || last(sampled_days) < last(sim_days))
        push!(sampled_days, last(sim_days))
        push!(sampled_altitude_km, sim_mean_altitude_km(df, nrow(df)))
    end

    return (; days=sampled_days, altitude_km=sampled_altitude_km)
end

function interpolate_value(xs::AbstractVector, ys::AbstractVector, x::Real)
    x < first(xs) && return missing
    x > last(xs) && return missing

    idx = searchsortedlast(xs, x)
    idx == length(xs) && return ys[end]

    x0 = xs[idx]
    x1 = xs[idx + 1]
    y0 = ys[idx]
    y1 = ys[idx + 1]

    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end

function centered_moving_average(values::AbstractVector{<:Real}, window::Integer)
    window > 0 || throw(ArgumentError("moving-average window must be positive"))
    half_window = window ÷ 2
    smoothed = Float64[]
    sizehint!(smoothed, length(values))

    for i in eachindex(values)
        lo = max(firstindex(values), i - half_window)
        hi = min(lastindex(values), i + half_window)
        push!(smoothed, mean(@view values[lo:hi]))
    end

    return smoothed
end


function print_initial_state_check(tle)
    first_tle = tle.elements[1]

    ic = (;
        a=6_706_404.4155,
        e=0.0000465,
        i=97.4411,
        Ω=326.1064,
        ω=216.2541,
        ν=143.8726,
    )

    println("\nInitial-state sanity check against first TLE epoch")
    println("TLE epoch: $(first_tle.epoch)")
    println("Element              IC              TLE       IC - TLE")

    for field in (:a, :e, :i, :Ω, :ω, :ν)
        @printf(
            "%-7s %15.8f %15.8f %14.8g\n",
            String(field),
            getfield(ic, field),
            getfield(first_tle, field),
            getfield(ic, field) - getfield(first_tle, field),
        )
    end
end

function print_discrepancy_summary(label::AbstractString, tle, sim)
    slice_days = [5.0, 10.0, 15.0, 20.0]

    println("\nAltitude discrepancy summary: $label")
    println("Day      TLE km      AGORA km   AGORA - TLE km")

    for day in slice_days
        tle_alt = interpolate_value(tle.days, tle.altitudes_km, day)
        sim_alt = interpolate_value(sim.days, sim.altitude_km, day)

        if ismissing(tle_alt) || ismissing(sim_alt)
            @printf("%4.0f %11s %11s %16s\n", day, "n/a", "n/a", "n/a")
        else
            @printf("%4.0f %11.3f %11.3f %16.3f\n", day, tle_alt, sim_alt, sim_alt - tle_alt)
        end
    end

    comparison_end_day = min(20.0, last(tle.days), last(sim.days))

    tle_decay_rate = (
        interpolate_value(tle.days, tle.altitudes_km, 0.0) -
        interpolate_value(tle.days, tle.altitudes_km, comparison_end_day)
    ) / comparison_end_day

    sim_decay_rate = (
        interpolate_value(sim.days, sim.altitude_km, 0.0) -
        interpolate_value(sim.days, sim.altitude_km, comparison_end_day)
    ) / comparison_end_day

    @printf(
        "\nDecay-rate comparison over days 0-%.1f: TLE %.3f km/day, AGORA %.3f km/day, ratio %.2f\n",
        comparison_end_day,
        tle_decay_rate,
        sim_decay_rate,
        sim_decay_rate / tle_decay_rate,
    )
end

function add_sim_case!(p, sim; raw_label::AbstractString, smooth_label::AbstractString, color, connect_raw::Bool=false)
    sim_smoothed_altitude_km = centered_moving_average(sim.altitude_km, 5)

    if connect_raw
        plot!(
            p,
            sim.days,
            sim.altitude_km;
            label = raw_label,
            color = color,
            linewidth = 1.5,
            alpha = 0.7,
        )
    end

    scatter!(
        p,
        sim.days,
        sim.altitude_km;
        label = connect_raw ? "" : raw_label,
        color = color,
        markersize = 2.5,
        alpha = 0.35,
    )

    plot!(
        p,
        sim.days,
        sim_smoothed_altitude_km;
        label = smooth_label,
        color = color,
        linewidth = 2.0,
        linestyle = :dash,
    )

    return p
end

function main()
    tle = load_tle_history(TLE_FILE)
    min_drag_sim = load_sim_history(MIN_DRAG_SIM_CSV, tle.days)
    medium_drag_sim = load_sim_history(MEDIUM_DRAG_SIM_CSV, tle.days; include_final_point=true)
    max_drag_sim = load_sim_history(MAX_DRAG_SIM_CSV, tle.days; include_final_point=true)

    print_initial_state_check(tle)
    print_discrepancy_summary("minimum-drag case", tle, min_drag_sim)
    print_discrepancy_summary("medium-drag case", tle, medium_drag_sim)
    print_discrepancy_summary("maximum-drag case", tle, max_drag_sim)

    p = plot(
        tle.days,
        tle.altitudes_km;
        label = "OPS-SAT TLE",
        color = :blue,
        linewidth = 2.5,
        xlabel = "Time [day]",
        ylabel = "Altitude [km]",
        title = "OPS-SAT Altitude",
        size = (1100, 650),
        dpi = 160,
        gridalpha = 0.3,
        minorgrid = true,
        legend = :outertopright,
        left_margin = 8Plots.mm,
    )

    add_sim_case!(
        p,
        min_drag_sim;
        raw_label = "AGORA min drag (raw)",
        smooth_label = "AGORA min drag (5-point moving average)",
        color = :red,
        connect_raw = false,
    )

    add_sim_case!(
        p,
        max_drag_sim;
        raw_label = "AGORA max drag (raw)",
        smooth_label = "AGORA max drag (5-point moving average)",
        color = :green,
        connect_raw = true,
    )

    mkpath(dirname(OUT_PNG))
    savefig(p, OUT_PNG)

    display(p)
    println("Saved plot to: $OUT_PNG")

    p_medium = plot(
        tle.days,
        tle.altitudes_km;
        label = "OPS-SAT TLE",
        color = :blue,
        linewidth = 2.5,
        xlabel = "Time [day]",
        ylabel = "Altitude [km]",
        title = "OPS-SAT Medium-Drag Altitude Comparison",
        size = (1100, 650),
        dpi = 160,
        gridalpha = 0.3,
        minorgrid = true,
        legend = :outertopright,
        left_margin = 8Plots.mm,
    )

    add_sim_case!(
        p_medium,
        medium_drag_sim;
        raw_label = "AGORA medium drag (raw)",
        smooth_label = "AGORA medium drag (5-point moving average)",
        color = :orange,
        connect_raw = false,
    )

    savefig(p_medium, MEDIUM_OUT_PNG)

    display(p_medium)
    println("Saved medium-drag plot to: $MEDIUM_OUT_PNG")
end

main()
