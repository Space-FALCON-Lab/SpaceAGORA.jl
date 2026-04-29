using CSV
using DataFrames
using Plots
using PlotlyJS: Plot, scatter3d, surface, Layout, attr
using Roots
using SPICE
using StaticArrays
using LinearAlgebra
using Logging
using Dates

Base.@kwdef struct AerobrakingMissionSpiceConfig
    mission_name::String
    target::String
    central_body::String
    observer::String
    body_fixed_frame::String
    kernel_relpaths::Tuple{Vararg{String}}
    surface_color::String = "rgb(201, 96, 56)"
end

const MARS_ODYSSEY_SPICE_CONFIG = AerobrakingMissionSpiceConfig(
    mission_name="Mars Odyssey",
    target="MARS ODYSSEY",
    central_body="MARS",
    observer="MARS BARYCENTER",
    body_fixed_frame="IAU_MARS",
    kernel_relpaths=("spk/missions/m01_ab_v2.bsp",),
    surface_color="rgb(201, 96, 56)"
)

const VEX_SPICE_CONFIG = AerobrakingMissionSpiceConfig(
    mission_name="Venus Express",
    target="VENUS EXPRESS",
    central_body="VENUS",
    observer="VENUS BARYCENTER",
    body_fixed_frame="IAU_VENUS",
    kernel_relpaths=(
        "spk/missions/ORVV_T19_140501000000_00546.BSP",
        "spk/missions/ORVV_T19_140601000000_00546.BSP",
        "spk/missions/ORVV_T19_140701000000_00551.BSP",
    ),
    surface_color="rgb(218, 170, 105)"
)

const _ACTIVE_MISSION_SPICE_CONFIG = Ref{AerobrakingMissionSpiceConfig}(MARS_ODYSSEY_SPICE_CONFIG)

function _set_aerobraking_mission_spice_config!(config::AerobrakingMissionSpiceConfig)
    _ACTIVE_MISSION_SPICE_CONFIG[] = config
    return config
end

_active_mission_spice_config() = _ACTIVE_MISSION_SPICE_CONFIG[]

function _mission_kernel_paths(spice_path::String, config::AerobrakingMissionSpiceConfig=_active_mission_spice_config())
    paths = map(relpath -> joinpath(spice_path, split(relpath, '/')...), config.kernel_relpaths)
    missing = filter(!isfile, paths)
    isempty(missing) || throw(ArgumentError("$(config.mission_name) mission SPICE kernel(s) not found: $(abspath.(missing))."))
    return collect(paths)
end

function _furnish_mission_kernels(spice_path::String, config::AerobrakingMissionSpiceConfig=_active_mission_spice_config())
    kernels = _mission_kernel_paths(spice_path, config)
    lock(RuntimeServices.SPICE_LOCK) do
        for kernel in kernels
            furnsh(kernel)
        end
    end
    return kernels
end

function _spice_relative_state_km(et::Float64, config::AerobrakingMissionSpiceConfig=_active_mission_spice_config())
    return lock(RuntimeServices.SPICE_LOCK) do
        sc_state_km, _ = spkezr(config.target, et, "J2000", "NONE", config.observer)
        body_state_km, _ = spkezr(config.central_body, et, "J2000", "NONE", config.observer)
        pos_km = SVector{3, Float64}(sc_state_km[1], sc_state_km[2], sc_state_km[3]) -
            SVector{3, Float64}(body_state_km[1], body_state_km[2], body_state_km[3])
        vel_kms = SVector{3, Float64}(sc_state_km[4], sc_state_km[5], sc_state_km[6]) -
            SVector{3, Float64}(body_state_km[4], body_state_km[5], body_state_km[6])
        pos_km, vel_kms
    end
end

function _spice_initial_condition_from_config(
    initial_time::InitialTime,
    spice_path::String,
    config::AerobrakingMissionSpiceConfig=_active_mission_spice_config()
)::CartesianInitialCondition
    et = _initial_time_et(initial_time)
    kernels = _furnish_mission_kernels(spice_path, config)
    pos_km, vel_kms = _spice_relative_state_km(et, config)
    pos_m = pos_km .* 1e3
    vel_mps = vel_kms .* 1e3
    println("Initialized $(config.mission_name) from SPICE kernel(s): $(join(abspath.(kernels), ", "))")
    println("  Epoch (ET seconds past J2000): $et")
    println("  Observer: $(config.observer), frame: J2000")
    println("  Position (m): $pos_m")
    println("  Velocity (m/s): $vel_mps")
    return CartesianInitialCondition(pos_m, vel_mps)
end

function _nearest_apoapsis_initial_condition_from_spice(
    initial_time::InitialTime,
    spice_path::String,
    config::AerobrakingMissionSpiceConfig=_active_mission_spice_config();
    search_window_s::Float64=14.0 * 24.0 * 3600.0,
    coarse_step_s::Float64=300.0
)::CartesianInitialCondition
    initial_et = _initial_time_et(initial_time)
    kernels = _furnish_mission_kernels(spice_path, config)
    offsets = collect(-search_window_s:coarse_step_s:search_window_s)
    rv = [_spice_radial_velocity(Float64(offset), initial_et, config) for offset in offsets]
    candidates = Float64[]
    for i in 1:(length(offsets) - 1)
        f0 = rv[i]
        f1 = rv[i + 1]
        isfinite(f0) && isfinite(f1) || continue
        f0 == 0.0 && push!(candidates, Float64(offsets[i]))
        if f0 >= 0.0 && f1 < 0.0
            push!(
                candidates,
                find_zero(t -> _spice_radial_velocity(Float64(t), initial_et, config), (offsets[i], offsets[i + 1]), Roots.Brent(); xatol=1e-9)
            )
        end
    end
    isempty(candidates) && throw(ArgumentError("No $(config.mission_name) SPICE apoapsis found within +/- $(search_window_s) seconds of the requested initial time."))
    apo_offset = candidates[argmin(abs.(candidates))]
    apo_et = initial_et + apo_offset
    pos_km, vel_kms = _spice_relative_state_km(apo_et, config)
    pos_m = pos_km .* 1e3
    vel_mps = vel_kms .* 1e3
    println("Initialized $(config.mission_name) from nearest SPICE apoapsis using kernel(s): $(join(abspath.(kernels), ", "))")
    println("  Requested epoch ET: $initial_et")
    println("  Apoapsis epoch ET: $apo_et")
    println("  Apoapsis offset from requested epoch (s): $apo_offset")
    println("  Observer: $(config.observer), frame: J2000")
    println("  Position (m): $pos_m")
    println("  Velocity (m/s): $vel_mps")
    return CartesianInitialCondition(pos_m, vel_mps)
end

function _initial_time_from_et(et::Float64)::InitialTime
    utc = lock(RuntimeServices.SPICE_LOCK) do
        et2utc(et, "ISOC", 6)
    end
    date_part, time_part = split(utc, 'T')
    year, month, day = parse.(Int, split(date_part, '-'))
    hour_str, minute_str, second_str = split(time_part, ':')
    return InitialTime(
        year=year,
        month=month,
        day=day,
        hour=parse(Int, hour_str),
        minute=parse(Int, minute_str),
        second=parse(Float64, second_str)
    )
end

function _nearest_apoapsis_initial_time_and_condition_from_spice(
    initial_time::InitialTime,
    spice_path::String,
    config::AerobrakingMissionSpiceConfig=_active_mission_spice_config();
    search_window_s::Float64=14.0 * 24.0 * 3600.0,
    coarse_step_s::Float64=300.0
)
    initial_et = _initial_time_et(initial_time)
    kernels = _furnish_mission_kernels(spice_path, config)
    offsets = collect(-search_window_s:coarse_step_s:search_window_s)
    rv = [_spice_radial_velocity(Float64(offset), initial_et, config) for offset in offsets]
    candidates = Float64[]
    for i in 1:(length(offsets) - 1)
        f0 = rv[i]
        f1 = rv[i + 1]
        isfinite(f0) && isfinite(f1) || continue
        f0 == 0.0 && push!(candidates, Float64(offsets[i]))
        if f0 >= 0.0 && f1 < 0.0
            push!(
                candidates,
                find_zero(t -> _spice_radial_velocity(Float64(t), initial_et, config), (offsets[i], offsets[i + 1]), Roots.Brent(); xatol=1e-9)
            )
        end
    end
    isempty(candidates) && throw(ArgumentError("No $(config.mission_name) SPICE apoapsis found within +/- $(search_window_s) seconds of the requested initial time."))
    apo_offset = candidates[argmin(abs.(candidates))]
    apo_et = initial_et + apo_offset
    pos_km, vel_kms = _spice_relative_state_km(apo_et, config)
    pos_m = pos_km .* 1e3
    vel_mps = vel_kms .* 1e3
    apo_initial_time = _initial_time_from_et(apo_et)
    println("Initialized $(config.mission_name) from nearest SPICE apoapsis using kernel(s): $(join(abspath.(kernels), ", "))")
    println("  Requested epoch ET: $initial_et")
    println("  Apoapsis epoch ET: $apo_et")
    println("  Apoapsis offset from requested epoch (s): $apo_offset")
    println("  Apoapsis InitialTime: $apo_initial_time")
    println("  Observer: $(config.observer), frame: J2000")
    println("  Position (m): $pos_m")
    println("  Velocity (m/s): $vel_mps")
    return (initial_time=apo_initial_time, ic=CartesianInitialCondition(pos_m, vel_mps), apo_offset_s=apo_offset)
end


function _initial_time_et(initial_time::InitialTime)::Float64
    timestamp = string(
        lpad(initial_time.year, 4, '0'), "-",
        lpad(initial_time.month, 2, '0'), "-",
        lpad(initial_time.day, 2, '0'), "T",
        lpad(initial_time.hour, 2, '0'), ":",
        lpad(initial_time.minute, 2, '0'), ":",
        lpad(round(Int, floor(initial_time.second)), 2, '0')
    )
    fractional_seconds = Float64(initial_time.second) - floor(Float64(initial_time.second))
    return lock(RuntimeServices.SPICE_LOCK) do
        str2et(timestamp) + fractional_seconds
    end
end

function _mars_odyssey_initial_condition_from_spice(initial_time::InitialTime, spice_path::String)::CartesianInitialCondition
    _set_aerobraking_mission_spice_config!(MARS_ODYSSEY_SPICE_CONFIG)
    return _spice_initial_condition_from_config(initial_time, spice_path, MARS_ODYSSEY_SPICE_CONFIG)
end

function _require_float_column(df::DataFrame, name::Symbol)::Vector{Float64}
    return Float64.(df[!, name])
end

function _correction_reset_times_for_results(csv_path::String)::Vector{Float64}
    summary_path = joinpath(dirname(dirname(csv_path)), "leg_summary.csv")
    isfile(summary_path) || return Float64[]
    summary = CSV.read(summary_path, DataFrame)
    :reset_time_s in propertynames(summary) || return Float64[]
    correction_applied = :correction_applied in propertynames(summary) ? Bool.(summary.correction_applied) : trues(nrow(summary))
    return Float64[
        Float64(t) for (t, applied) in zip(summary.reset_time_s, correction_applied)
        if applied && isfinite(Float64(t))
    ]
end

function _interval_touches_correction(t0::Float64, t1::Float64, correction_times::Vector{Float64}; atol::Float64=1e-6)::Bool
    isempty(correction_times) && return false
    lo, hi = minmax(t0, t1)
    for correction_t in correction_times
        if lo - atol <= correction_t <= hi + atol
            return true
        end
    end
    return false
end

function _derive_orbit_extrema_from_results(csv_path::String, planet)
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))

    df = CSV.read(csv_path, DataFrame)
    nrow(df) >= 3 || throw(ArgumentError("Need at least 3 saved samples to derive periapsis/apoapsis extrema."))

    time_s = _require_float_column(df, :time)
    correction_times = _correction_reset_times_for_results(csv_path)
    pos_x_m = _require_float_column(df, :sc1_pos_1)
    pos_y_m = _require_float_column(df, :sc1_pos_2)
    pos_z_m = _require_float_column(df, :sc1_pos_3)
    vel_x_mps = _require_float_column(df, :sc1_vel_1)
    vel_y_mps = _require_float_column(df, :sc1_vel_2)
    vel_z_mps = _require_float_column(df, :sc1_vel_3)
    altitude_km = _require_float_column(df, :sc1_altitude) ./ 1e3
    latitude_deg = _require_float_column(df, :sc1_latitude_deg)
    longitude_deg = _require_float_column(df, :sc1_longitude_deg)
    radius_km = [
        norm(SVector{3, Float64}(pos_x_m[i], pos_y_m[i], pos_z_m[i])) / 1e3
        for i in eachindex(pos_x_m)
    ]
    radial_velocity_mps = [
        dot(
            SVector{3, Float64}(pos_x_m[i], pos_y_m[i], pos_z_m[i]),
            SVector{3, Float64}(vel_x_mps[i], vel_y_mps[i], vel_z_mps[i])
        ) / (radius_km[i] * 1e3)
        for i in eachindex(pos_x_m)
    ]

    peri_alt_km = Float64[]
    peri_lat_deg = Float64[]
    peri_lon_deg = Float64[]
    apo_alt_km = Float64[]
    apo_radius_km = Float64[]

    if radial_velocity_mps[1] <= 0.0
        push!(apo_alt_km, altitude_km[1])
        push!(apo_radius_km, radius_km[1])
    end

    for i in 2:(length(altitude_km) - 1)
        if _interval_touches_correction(time_s[i - 1], time_s[i], correction_times) ||
                _interval_touches_correction(time_s[i], time_s[i + 1], correction_times)
            continue
        end
        prev_alt = altitude_km[i - 1]
        curr_alt = altitude_km[i]
        next_alt = altitude_km[i + 1]
        if curr_alt <= prev_alt && curr_alt < next_alt
            push!(peri_alt_km, curr_alt)
            push!(peri_lat_deg, latitude_deg[i])
            push!(peri_lon_deg, longitude_deg[i])
        end
    end
    for i in 1:(length(radial_velocity_mps) - 1)
        _interval_touches_correction(time_s[i], time_s[i + 1], correction_times) && continue
        rv0 = radial_velocity_mps[i]
        rv1 = radial_velocity_mps[i + 1]
        if !(isfinite(rv0) && isfinite(rv1))
            continue
        end
        if !(rv0 > 0.0 && rv1 <= 0.0)
            continue
        end
        frac = rv0 / (rv0 - rv1)
        push!(apo_alt_km, altitude_km[i] + frac * (altitude_km[i + 1] - altitude_km[i]))
        push!(apo_radius_km, radius_km[i] + frac * (radius_km[i + 1] - radius_km[i]))
    end

    isempty(peri_alt_km) && throw(ArgumentError("No periapsis extrema found in saved aerobraking trajectory."))
    isempty(apo_alt_km) && throw(ArgumentError("No apoapsis extrema found in saved aerobraking trajectory."))

    return (
        peri=(orbit=collect(1:length(peri_alt_km)), altitude_km=peri_alt_km, latitude_deg=peri_lat_deg, longitude_deg=peri_lon_deg),
        apo=(orbit=collect(1:length(apo_alt_km)), altitude_km=apo_alt_km, radius_km=apo_radius_km)
    )
end

function _has_orbit_extrema_for_plots(args::SimulationConfiguration, planet)::Bool
    csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
    try
        _derive_orbit_extrema_from_results(csv_path, planet)
        return true
    catch err
        if err isa ArgumentError
            @warn "Skipping orbit-extrema plot subset because the saved trajectory does not contain both periapsis and apoapsis extrema." reason=sprint(showerror, err)
            return false
        end
        rethrow()
    end
end

function _state_pos_vel(u, sat_idx::Int=1)
    return SVector{3, Float64}(u.sc[sat_idx].pos), SVector{3, Float64}(u.sc[sat_idx].vel)
end

function _radial_velocity_from_solution(sol, t::Float64, sat_idx::Int=1)::Float64
    pos, vel = _state_pos_vel(sol(t), sat_idx)
    return dot(pos, vel) / norm(pos)
end

function _empty_periapsis_event_dataframe()
    return DataFrame(
        orbit=Int[],
        time_s=Float64[],
        altitude_km=Float64[],
        latitude_deg=Float64[],
        longitude_deg=Float64[],
        radius_km=Float64[],
        radial_velocity_mps=Float64[]
    )
end

function _periapsis_event_dataframe(sol, args::SimulationConfiguration, planet; sat_idx::Int=1, allow_empty::Bool=false)
    ephemerides_model = args.environment_model.ephemerides_model
    et_start = SimulationModel.ephemerides_time_seconds(args.initial_time, ephemerides_model)
    event_time_s = Float64[]
    altitude_km = Float64[]
    latitude_deg = Float64[]
    longitude_deg = Float64[]
    radius_km = Float64[]
    radial_velocity_mps = Float64[]

    ts = Float64.(sol.t)
    length(ts) >= 2 || throw(ArgumentError("Need at least two solver time points to locate periapsis events."))

    for i in 1:(length(ts) - 1)
        t0 = ts[i]
        t1 = ts[i + 1]
        t1 > t0 || continue
        f0 = _radial_velocity_from_solution(sol, t0, sat_idx)
        f1 = _radial_velocity_from_solution(sol, t1, sat_idx)
        if !(isfinite(f0) && isfinite(f1))
            continue
        end
        if !(f0 <= 0.0 && f1 > 0.0)
            continue
        end

        t_peri = if f0 == 0.0
            t0
        else
            find_zero(t -> _radial_velocity_from_solution(sol, Float64(t), sat_idx), (t0, t1), Roots.Brent(); xatol=1e-9)
        end
        if !isempty(event_time_s) && abs(t_peri - event_time_s[end]) <= 1e-6
            continue
        end

        pos, vel = _state_pos_vel(sol(t_peri), sat_idx)
        et = et_start + t_peri
        rp, _ = SimulationModel.SimulationCallbacks.r_intor_p!(pos, vel, planet, et, ephemerides_model)
        alt, lat, lon = SimulationModel.SimulationCallbacks.rtolatlong(rp, planet, ephemerides_model)

        push!(event_time_s, Float64(t_peri))
        push!(altitude_km, alt / 1e3)
        push!(latitude_deg, rad2deg(lat))
        push!(longitude_deg, rad2deg(lon))
        push!(radius_km, norm(pos) / 1e3)
        push!(radial_velocity_mps, _radial_velocity_from_solution(sol, Float64(t_peri), sat_idx))
    end

    if isempty(event_time_s)
        allow_empty && return _empty_periapsis_event_dataframe()
        throw(ArgumentError("No event-located periapsis crossings found in aerobraking solution."))
    end
    return DataFrame(
        orbit=collect(1:length(event_time_s)),
        time_s=event_time_s,
        altitude_km=altitude_km,
        latitude_deg=latitude_deg,
        longitude_deg=longitude_deg,
        radius_km=radius_km,
        radial_velocity_mps=radial_velocity_mps
    )
end

function _save_periapsis_event_table(sol, args::SimulationConfiguration, planet; allow_empty::Bool=false)
    results_dir = args.simulation_settings.results_directory
    mkpath(results_dir)
    events = _periapsis_event_dataframe(sol, args, planet; allow_empty=allow_empty)
    event_path = joinpath(results_dir, "periapsis_events.csv")
    CSV.write(event_path, events)
    println("Saved event-located periapsis table to $(abspath(event_path))")
    return event_path
end

function _load_periapsis_events(results_dir::String)::DataFrame
    event_path = joinpath(results_dir, "periapsis_events.csv")
    isfile(event_path) || throw(ArgumentError("Periapsis event table not found at $(abspath(event_path)). Run _save_periapsis_event_table before plotting periapsis lat/lon."))
    events = CSV.read(event_path, DataFrame)
    required = (:orbit, :time_s, :altitude_km, :latitude_deg, :longitude_deg)
    missing_cols = setdiff(required, Symbol.(names(events)))
    isempty(missing_cols) || throw(ArgumentError("Periapsis event table missing columns: $(missing_cols)."))
    return events
end

function _spice_radial_velocity(
    time_s::Float64,
    initial_et::Float64,
    config::AerobrakingMissionSpiceConfig=_active_mission_spice_config()
)::Float64
    pos_km, vel_kms = _spice_relative_state_km(initial_et + time_s, config)
    return dot(pos_km, vel_kms) / norm(pos_km)
end

function _spice_periapsis_event_dataframe(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String; allow_empty::Bool=false)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))
    _furnish_mission_kernels(spice_path)

    df = CSV.read(csv_path, DataFrame)
    time_s = _require_float_column(df, :time)
    length(time_s) >= 2 || throw(ArgumentError("Need at least two saved time samples to bracket SPICE periapsis events."))

    ephemerides_model = args.environment_model.ephemerides_model
    initial_et = _initial_time_et(initial_time)
    event_time_s = Float64[]
    altitude_km = Float64[]
    latitude_deg = Float64[]
    longitude_deg = Float64[]
    radius_km = Float64[]
    radial_velocity_kms = Float64[]

    @inbounds for i in 1:(length(time_s) - 1)
        t0 = time_s[i]
        t1 = time_s[i + 1]
        t1 > t0 || continue
        f0 = _spice_radial_velocity(t0, initial_et)
        f1 = _spice_radial_velocity(t1, initial_et)
        if !(isfinite(f0) && isfinite(f1))
            continue
        end
        if !(f0 <= 0.0 && f1 > 0.0)
            continue
        end

        t_peri = if f0 == 0.0
            t0
        else
            find_zero(t -> _spice_radial_velocity(Float64(t), initial_et), (t0, t1), Roots.Brent(); xatol=1e-9)
        end
        if !isempty(event_time_s) && abs(t_peri - event_time_s[end]) <= 1e-6
            continue
        end

        et = initial_et + t_peri
        pos_km, vel_kms = _spice_relative_state_km(et)
        pos_m = pos_km .* 1e3
        vel_mps = vel_kms .* 1e3
        rp, _ = SimulationModel.SimulationCallbacks.r_intor_p!(pos_m, vel_mps, planet, et, ephemerides_model)
        alt, lat, lon = SimulationModel.SimulationCallbacks.rtolatlong(rp, planet, ephemerides_model)

        push!(event_time_s, Float64(t_peri))
        push!(altitude_km, alt / 1e3)
        push!(latitude_deg, rad2deg(lat))
        push!(longitude_deg, rad2deg(lon))
        push!(radius_km, norm(pos_km))
        push!(radial_velocity_kms, _spice_radial_velocity(Float64(t_peri), initial_et))
    end

    if isempty(event_time_s)
        allow_empty && return DataFrame(
            orbit=Int[],
            time_s=Float64[],
            altitude_km=Float64[],
            latitude_deg=Float64[],
            longitude_deg=Float64[],
            radius_km=Float64[],
            radial_velocity_kms=Float64[]
        )
        throw(ArgumentError("No SPICE periapsis crossings found over saved aerobraking time span."))
    end
    return DataFrame(
        orbit=collect(1:length(event_time_s)),
        time_s=event_time_s,
        altitude_km=altitude_km,
        latitude_deg=latitude_deg,
        longitude_deg=longitude_deg,
        radius_km=radius_km,
        radial_velocity_kms=radial_velocity_kms
    )
end

function _spice_apoapsis_event_dataframe(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String; allow_empty::Bool=false)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))
    _furnish_mission_kernels(spice_path)

    df = CSV.read(csv_path, DataFrame)
    time_s = _require_float_column(df, :time)
    length(time_s) >= 2 || throw(ArgumentError("Need at least two saved time samples to bracket SPICE apoapsis events."))

    ephemerides_model = args.environment_model.ephemerides_model
    initial_et = _initial_time_et(initial_time)
    event_time_s = Float64[]
    altitude_km = Float64[]
    latitude_deg = Float64[]
    longitude_deg = Float64[]
    radius_km = Float64[]
    radial_velocity_kms = Float64[]

    if _spice_radial_velocity(0.0, initial_et) <= 0.0
        pos_km, vel_kms = _spice_relative_state_km(initial_et)
        pos_m = pos_km .* 1e3
        vel_mps = vel_kms .* 1e3
        rp, _ = SimulationModel.SimulationCallbacks.r_intor_p!(pos_m, vel_mps, planet, initial_et, ephemerides_model)
        alt, lat, lon = SimulationModel.SimulationCallbacks.rtolatlong(rp, planet, ephemerides_model)

        push!(event_time_s, 0.0)
        push!(altitude_km, alt / 1e3)
        push!(latitude_deg, rad2deg(lat))
        push!(longitude_deg, rad2deg(lon))
        push!(radius_km, norm(pos_km))
        push!(radial_velocity_kms, _spice_radial_velocity(0.0, initial_et))
    end

    @inbounds for i in 1:(length(time_s) - 1)
        t0 = time_s[i]
        t1 = time_s[i + 1]
        t1 > t0 || continue
        f0 = _spice_radial_velocity(t0, initial_et)
        f1 = _spice_radial_velocity(t1, initial_et)
        if !(isfinite(f0) && isfinite(f1))
            continue
        end
        if !(f0 >= 0.0 && f1 < 0.0)
            continue
        end

        t_apo = if f0 == 0.0
            t0
        else
            find_zero(t -> _spice_radial_velocity(Float64(t), initial_et), (t0, t1), Roots.Brent(); xatol=1e-9)
        end
        if !isempty(event_time_s) && abs(t_apo - event_time_s[end]) <= 1e-6
            continue
        end

        et = initial_et + t_apo
        pos_km, vel_kms = _spice_relative_state_km(et)
        pos_m = pos_km .* 1e3
        vel_mps = vel_kms .* 1e3
        rp, _ = SimulationModel.SimulationCallbacks.r_intor_p!(pos_m, vel_mps, planet, et, ephemerides_model)
        alt, lat, lon = SimulationModel.SimulationCallbacks.rtolatlong(rp, planet, ephemerides_model)

        push!(event_time_s, Float64(t_apo))
        push!(altitude_km, alt / 1e3)
        push!(latitude_deg, rad2deg(lat))
        push!(longitude_deg, rad2deg(lon))
        push!(radius_km, norm(pos_km))
        push!(radial_velocity_kms, _spice_radial_velocity(Float64(t_apo), initial_et))
    end

    if isempty(event_time_s)
        allow_empty && return DataFrame(
            orbit=Int[],
            time_s=Float64[],
            altitude_km=Float64[],
            latitude_deg=Float64[],
            longitude_deg=Float64[],
            radius_km=Float64[],
            radial_velocity_kms=Float64[]
        )
        throw(ArgumentError("No SPICE apoapsis crossings found over saved aerobraking time span."))
    end
    return DataFrame(
        orbit=collect(1:length(event_time_s)),
        time_s=event_time_s,
        altitude_km=altitude_km,
        latitude_deg=latitude_deg,
        longitude_deg=longitude_deg,
        radius_km=radius_km,
        radial_velocity_kms=radial_velocity_kms
    )
end

function _save_spice_periapsis_event_table(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String; allow_empty::Bool=false)
    results_dir = args.simulation_settings.results_directory
    mkpath(results_dir)
    events = _spice_periapsis_event_dataframe(args, planet, initial_time, spice_path; allow_empty=allow_empty)
    event_path = joinpath(results_dir, "spice_periapsis_events.csv")
    CSV.write(event_path, events)
    println("Saved SPICE periapsis table to $(abspath(event_path))")
    return event_path
end

function _save_spice_apoapsis_event_table(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String; allow_empty::Bool=false)
    results_dir = args.simulation_settings.results_directory
    mkpath(results_dir)
    events = _spice_apoapsis_event_dataframe(args, planet, initial_time, spice_path; allow_empty=allow_empty)
    event_path = joinpath(results_dir, "spice_apoapsis_events.csv")
    CSV.write(event_path, events)
    println("Saved SPICE apoapsis table to $(abspath(event_path))")
    return event_path
end

function _load_spice_periapsis_events(results_dir::String)::DataFrame
    event_path = joinpath(results_dir, "spice_periapsis_events.csv")
    isfile(event_path) || return DataFrame()
    events = CSV.read(event_path, DataFrame)
    required = (:orbit, :time_s, :altitude_km, :latitude_deg, :longitude_deg)
    missing_cols = setdiff(required, Symbol.(names(events)))
    isempty(missing_cols) || throw(ArgumentError("SPICE periapsis event table missing columns: $(missing_cols)."))
    return events
end

function _load_spice_apoapsis_events(results_dir::String)::DataFrame
    event_path = joinpath(results_dir, "spice_apoapsis_events.csv")
    isfile(event_path) || return DataFrame()
    events = CSV.read(event_path, DataFrame)
    required = (:orbit, :time_s, :altitude_km, :latitude_deg, :longitude_deg)
    missing_cols = setdiff(required, Symbol.(names(events)))
    isempty(missing_cols) || throw(ArgumentError("SPICE apoapsis event table missing columns: $(missing_cols)."))
    return events
end

function _split_maneuver_orbits(odyssey_schedule, max_orbit_plotted::Int)
    lower_orbits = Int[]
    raise_orbits = Int[]
    for (orbit, signed_delta_v) in zip(odyssey_schedule.maneuver_orbit_number, odyssey_schedule.maneuver_Δv)
        orbit > max_orbit_plotted && continue
        # Maneuvers are centered near apoapsis. Positive signed Δv entries are
        # periapsis-raise maneuvers; negative signed Δv entries lower periapsis.
        if signed_delta_v > 0.0
            push!(raise_orbits, orbit)
        elseif signed_delta_v < 0.0
            push!(lower_orbits, orbit)
        end
    end
    return raise_orbits, lower_orbits
end

function _save_apoapsis_periapsis_plot(args::SimulationConfiguration, planet, odyssey_schedule)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    extrema = _derive_orbit_extrema_from_results(csv_path, planet)
    peri_events = _load_periapsis_events(results_dir)
    spice_peri_events = _load_spice_periapsis_events(results_dir)
    spice_apo_events = _load_spice_apoapsis_events(results_dir)
    max_orbit_plotted = max(nrow(peri_events), length(extrema.apo.orbit), nrow(spice_peri_events), nrow(spice_apo_events))
    raise_orbits, lower_orbits = _split_maneuver_orbits(odyssey_schedule, max_orbit_plotted)

    p = plot(
        peri_events.orbit,
        peri_events.altitude_km;
        xlabel="Orbit Number",
        ylabel="Periapsis Altitude (km)",
        label="Simulation Periapsis",
        linewidth=2.5,
        marker=:circle,
        color=:dodgerblue,
        grid=true,
        legend=:topright
    )
    if nrow(spice_peri_events) > 0
        plot!(
            p,
            spice_peri_events.orbit,
            spice_peri_events.altitude_km;
            label="SPICE Periapsis",
            linewidth=2.0,
            linestyle=:dash,
            marker=:xcross,
            color=:navy
        )
    end
    if !isempty(raise_orbits)
        vline!(
            p,
            raise_orbits;
            color=:seagreen3,
            linestyle=:dash,
            linewidth=1.5,
            alpha=0.7,
            label="Periapsis Raise Maneuver"
        )
    end
    if !isempty(lower_orbits)
        vline!(
            p,
            lower_orbits;
            color=:darkorange2,
            linestyle=:dashdot,
            linewidth=1.5,
            alpha=0.7,
            label="Periapsis Lower Maneuver"
        )
    end
    apo_axis = twinx()
    plot!(
        apo_axis,
        extrema.apo.orbit,
        extrema.apo.altitude_km;
        ylabel="Apoapsis Altitude (km)",
        label=false,
        linewidth=2.5,
        marker=:diamond,
        color=:crimson
    )
    plot!(
        p,
        [NaN],
        [NaN];
        label="Simulation Apoapsis",
        linewidth=2.5,
        marker=:diamond,
        color=:crimson
    )
    if nrow(spice_apo_events) > 0
        plot!(
            apo_axis,
            spice_apo_events.orbit,
            spice_apo_events.altitude_km;
            label=false,
            linewidth=2.0,
            linestyle=:dash,
            marker=:xcross,
            color=:darkred
        )
        plot!(
            p,
            [NaN],
            [NaN];
            label="SPICE Apoapsis",
            linewidth=2.0,
            linestyle=:dash,
            marker=:xcross,
            color=:darkred
        )
    end

    plot_path = joinpath(results_dir, "apoapsis_periapsis_vs_orbit.png")
    savefig(p, plot_path)
    println("Saved apoapsis/periapsis plot to $(abspath(plot_path))")
    return plot_path
end

function _save_trajectory_components_plot(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    time_hr = time_s ./ 3600.0

    p = plot(
        time_hr,
        sim_x_km;
        xlabel="Mission Time (hr)",
        ylabel="Position Component (km)",
        label="x",
        linewidth=2.0,
        color=:royalblue,
        grid=true,
        legend=:topright,
        title="Inertial Trajectory Components"
    )
    plot!(p, time_hr, sim_y_km; label="y", linewidth=2.0, color=:firebrick)
    plot!(p, time_hr, sim_z_km; label="z", linewidth=2.0, color=:seagreen)
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "trajectory_components_vs_time.png")
    savefig(p, plot_path)
    println("Saved trajectory-components plot to $(abspath(plot_path))")
    return plot_path
end

function _save_reference_sphere_altitude_plot(args::SimulationConfiguration; reference_radius_km::Float64=6378.0)
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    time_hr = time_s ./ 3600.0
    altitude_km = [
        sqrt(sim_x_km[i]^2 + sim_y_km[i]^2 + sim_z_km[i]^2) - reference_radius_km
        for i in eachindex(time_s)
    ]

    p = plot(
        time_hr,
        altitude_km;
        xlabel="Mission Time (hr)",
        ylabel="Altitude Above 6378 km Sphere (km)",
        label="Altitude",
        linewidth=2.0,
        color=:darkcyan,
        grid=true,
        legend=:topright,
        title="Altitude Relative to 6378 km Reference Sphere"
    )
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "altitude_vs_time_reference_sphere_6378km.png")
    savefig(p, plot_path)
    println("Saved reference-sphere altitude plot to $(abspath(plot_path))")
    return plot_path
end

function _initial_time_datetime(initial_time::InitialTime)::DateTime
    whole_seconds = floor(Int, Float64(initial_time.second))
    millis = round(Int, (Float64(initial_time.second) - whole_seconds) * 1e3)
    base = DateTime(
        Int(initial_time.year),
        Int(initial_time.month),
        Int(initial_time.day),
        Int(initial_time.hour),
        Int(initial_time.minute),
        whole_seconds
    )
    return base + Millisecond(millis)
end

function _parse_tle_timestamp(line1::AbstractString)::DateTime
    epoch_year = parse(Int, strip(line1[19:20]))
    epoch_day = parse(Float64, strip(line1[21:32]))
    year_full = epoch_year >= 57 ? 1900 + epoch_year : 2000 + epoch_year
    base = DateTime(year_full, 1, 1, 0, 0, 0)
    millis = round(Int, (epoch_day - 1.0) * 86400.0 * 1e3)
    return base + Millisecond(millis)
end

function _parse_telemetry_timestamp(value)::DateTime
    s = strip(string(value))
    sep = occursin('T', s) ? 'T' : ' '
    occursin(sep, s) || throw(ArgumentError("Unable to parse telemetry timestamp '$s'."))
    date_part, time_part = split(s, sep; limit=2)
    frac_us = 0
    if occursin('.', time_part)
        whole, frac = split(time_part, '.'; limit=2)
        frac_digits = replace(frac, r"[^0-9]" => "")
        frac_digits = isempty(frac_digits) ? "0" : frac_digits
        frac_trimmed = first(frac_digits, min(length(frac_digits), 6))
        frac_us = parse(Int, rpad(frac_trimmed, 6, '0'))
        time_part = whole
    end
    base = DateTime(string(date_part, " ", time_part), dateformat"yyyy-mm-dd HH:MM:SS")
    return base + Microsecond(frac_us)
end

function _load_tle_history(path::String)
    raw_lines = [strip(line) for line in readlines(path) if !isempty(strip(line))]
    length(raw_lines) >= 2 || throw(ArgumentError("TLE file at $(abspath(path)) must contain at least two lines."))
    first_line = raw_lines[1]
    step = startswith(first_line, "1 ") ? 2 : 3
    entries = NamedTuple[]
    idx = 1
    while idx <= length(raw_lines) - 1
        line1_idx = step == 2 ? idx : idx + 1
        line2_idx = step == 2 ? idx + 1 : idx + 2
        line2_idx <= length(raw_lines) || break
        line1 = raw_lines[line1_idx]
        line2 = raw_lines[line2_idx]
        startswith(line1, "1 ") || (idx += step; continue)
        startswith(line2, "2 ") || (idx += step; continue)
        push!(entries, (timestamp=_parse_tle_timestamp(line1), line1=line1, line2=line2))
        idx += step
    end
    isempty(entries) && throw(ArgumentError("No valid TLE pairs found in $(abspath(path))."))
    sort!(entries, by=x -> x.timestamp)
    return entries
end

@inline function _tle_mean_altitude_km(line2::AbstractString; reference_radius_km::Float64=6378.0)
    mean_motion_rev_day = parse(Float64, strip(line2[53:63]))
    n_rad_s = mean_motion_rev_day * 2π / 86400.0
    μ_earth_km3_s2 = 398600.4418
    semi_major_axis_km = cbrt(μ_earth_km3_s2 / (n_rad_s^2))
    return semi_major_axis_km - reference_radius_km
end

function _closest_tle_entry(entries, target::DateTime)
    timestamps = getfield.(entries, :timestamp)
    idx = searchsortedfirst(timestamps, target)
    if idx <= 1
        return entries[1]
    elseif idx > length(entries)
        return entries[end]
    end
    prev_entry = entries[idx - 1]
    next_entry = entries[idx]
    dt_prev = abs(Dates.value(target - prev_entry.timestamp))
    dt_next = abs(Dates.value(next_entry.timestamp - target))
    return dt_prev <= dt_next ? prev_entry : next_entry
end

@inline function _interp_linear_samples(x::Vector{Float64}, y::Vector{Float64}, xq::Vector{Float64})::Vector{Float64}
    out = Vector{Float64}(undef, length(xq))
    j = 1
    @inbounds for i in eachindex(xq)
        xi = xq[i]
        while j < length(x) - 1 && x[j + 1] < xi
            j += 1
        end
        if xi <= x[1]
            out[i] = y[1]
        elseif xi >= x[end]
            out[i] = y[end]
        else
            x0 = x[j]
            x1 = x[j + 1]
            y0 = y[j]
            y1 = y[j + 1]
            frac = (xi - x0) / (x1 - x0)
            out[i] = y0 + frac * (y1 - y0)
        end
    end
    return out
end

function _save_telemetry_simulation_altitude_plot(
    args::SimulationConfiguration,
    initial_time::InitialTime;
    telemetry_csv_path::String=joinpath(REPO_ROOT, "data", "telemetry", "OPS-SAT1", "UHF_TM_notebook", "uhf_telemetry.csv"),
    tle_path::String=joinpath(REPO_ROOT, "data", "telemetry", "OPS-SAT1", "UHF_TM_notebook", "tle_opssat.txt"),
    reference_radius_km::Float64=6378.0
)
    results_dir = args.simulation_settings.results_directory
    sim_time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    sim_altitude_km = [
        sqrt(sim_x_km[i]^2 + sim_y_km[i]^2 + sim_z_km[i]^2) - reference_radius_km
        for i in eachindex(sim_time_s)
    ]
    sim_start = _initial_time_datetime(initial_time)
    sim_end = sim_start + Millisecond(round(Int, sim_time_s[end] * 1e3))

    telemetry_df = CSV.read(telemetry_csv_path, DataFrame)
    :Timestamp in propertynames(telemetry_df) || throw(ArgumentError("Telemetry CSV missing Timestamp column: $(abspath(telemetry_csv_path))."))
    telemetry_times = _parse_telemetry_timestamp.(telemetry_df.Timestamp)
    in_window = [sim_start <= t <= sim_end for t in telemetry_times]
    any(in_window) || throw(ArgumentError("No telemetry timestamps overlap the simulation window $(sim_start) to $(sim_end)."))
    telemetry_times = telemetry_times[in_window]
    telemetry_order = sortperm(telemetry_times)
    telemetry_times = telemetry_times[telemetry_order]

    tle_entries = _load_tle_history(tle_path)
    telemetry_altitude_km = Vector{Float64}(undef, length(telemetry_times))
    telemetry_elapsed_s = Vector{Float64}(undef, length(telemetry_times))
    @inbounds for i in eachindex(telemetry_times)
        entry = _closest_tle_entry(tle_entries, telemetry_times[i])
        telemetry_altitude_km[i] = _tle_mean_altitude_km(entry.line2; reference_radius_km=reference_radius_km)
        telemetry_elapsed_s[i] = Dates.value(telemetry_times[i] - sim_start) / 1000
    end
    sim_interp_altitude_km = _interp_linear_samples(sim_time_s, sim_altitude_km, telemetry_elapsed_s)

    p = plot(
        telemetry_times,
        telemetry_altitude_km;
        xlabel="Timestamp (UTC)",
        ylabel="Altitude Above $(reference_radius_km) km Sphere (km)",
        label="Telemetry / Closest-TLE",
        linewidth=1.2,
        marker=:circle,
        markersize=2.5,
        color=:black,
        alpha=0.7,
        grid=true,
        legend=:topright,
        title="Telemetry and Simulation Altitude at Corresponding Timestamps"
    )
    plot!(
        p,
        telemetry_times,
        sim_interp_altitude_km;
        label="Simulation",
        linewidth=2.2,
        color=:darkcyan
    )

    plot_path = joinpath(results_dir, "telemetry_vs_simulation_altitude_6378km.png")
    savefig(p, plot_path)
    println("Saved telemetry/simulation altitude comparison plot to $(abspath(plot_path))")
    return plot_path
end

function _save_orbital_elements_plot(args::SimulationConfiguration, planet)
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    _, sim_vx_kms, sim_vy_kms, sim_vz_kms = _simulation_velocity_samples(args)
    sim_oe = _orbital_elements_from_state_samples(
        sim_x_km,
        sim_y_km,
        sim_z_km,
        sim_vx_kms,
        sim_vy_kms,
        sim_vz_kms,
        planet
    )
    time_hr = time_s ./ 3600.0

    specs = (
        (:semi_major_axis_km, "Semi-major Axis (km)", :royalblue),
        (:eccentricity, "Eccentricity", :firebrick),
        (:inclination_deg, "Inclination (deg)", :seagreen),
        (:raan_deg, "RAAN (deg)", :darkorange2),
        (:arg_periapsis_deg, "Argument of Periapsis (deg)", :purple),
        (:true_anomaly_deg, "True Anomaly (deg)", :teal)
    )

    subplots = Plots.Plot[]
    for (field, ylabel, color) in specs
        sp = plot(
            time_hr,
            getproperty(sim_oe, field);
            ylabel=ylabel,
            label=false,
            linewidth=2.0,
            color=color,
            grid=true
        )
        push!(subplots, sp)
    end
    xlabel!(subplots[5], "Mission Time (hr)")
    xlabel!(subplots[6], "Mission Time (hr)")

    p = plot(
        subplots...;
        layout=(3, 2),
        size=(1300, 1000),
        plot_title="Orbital Elements vs Time"
    )
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "orbital_elements_vs_time.png")
    savefig(p, plot_path)
    println("Saved orbital-elements plot to $(abspath(plot_path))")
    return plot_path
end

function _save_apoapsis_periapsis_trace_plot(args::SimulationConfiguration, planet)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    extrema = _derive_orbit_extrema_from_results(csv_path, planet)

    p = plot(
        extrema.peri.orbit,
        extrema.peri.altitude_km;
        xlabel="Orbit Number",
        ylabel="Periapsis Altitude (km)",
        label="Periapsis",
        linewidth=2.5,
        marker=:circle,
        color=:dodgerblue,
        grid=true,
        legend=:topright,
        title="Apoapsis and Periapsis Traces"
    )

    apo_axis = twinx()
    plot!(
        apo_axis,
        extrema.apo.orbit,
        extrema.apo.altitude_km;
        ylabel="Apoapsis Altitude (km)",
        label=false,
        linewidth=2.5,
        marker=:diamond,
        color=:crimson
    )
    plot!(
        p,
        [NaN],
        [NaN];
        label="Apoapsis",
        linewidth=2.5,
        marker=:diamond,
        color=:crimson
    )

    plot_path = joinpath(results_dir, "apoapsis_periapsis_trace_vs_orbit.png")
    savefig(p, plot_path)
    println("Saved apoapsis/periapsis trace plot to $(abspath(plot_path))")
    return plot_path
end

function _save_periapsis_latlon_plot(args::SimulationConfiguration, planet, odyssey_schedule)
    results_dir = args.simulation_settings.results_directory
    peri_events = _load_periapsis_events(results_dir)
    spice_peri_events = _load_spice_periapsis_events(results_dir)
    if nrow(peri_events) == 0
        @warn "Skipping periapsis latitude/longitude plot because no event-located periapsis crossings were found." results_dir
        return nothing
    end
    max_orbit_plotted = max(nrow(peri_events), nrow(spice_peri_events))
    raise_orbits, lower_orbits = _split_maneuver_orbits(odyssey_schedule, max_orbit_plotted)

    p = plot(
        peri_events.orbit,
        peri_events.latitude_deg;
        xlabel="Orbit Number",
        ylabel="Periapsis Latitude (deg)",
        label="Periapsis Latitude",
        linewidth=2.5,
        marker=:circle,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    if nrow(spice_peri_events) > 0
        plot!(
            p,
            spice_peri_events.orbit,
            spice_peri_events.latitude_deg;
            label="SPICE Periapsis Latitude",
            linewidth=2.0,
            linestyle=:dash,
            marker=:xcross,
            color=:navy
        )
    end
    if !isempty(raise_orbits)
        vline!(
            p,
            raise_orbits;
            color=:seagreen3,
            linestyle=:dash,
            linewidth=1.5,
            alpha=0.7,
            label="Periapsis Raise Maneuver"
        )
    end
    if !isempty(lower_orbits)
        vline!(
            p,
            lower_orbits;
            color=:darkorange2,
            linestyle=:dashdot,
            linewidth=1.5,
            alpha=0.7,
            label="Periapsis Lower Maneuver"
        )
    end
    lon_axis = twinx()
    plot!(
        lon_axis,
        peri_events.orbit,
        peri_events.longitude_deg;
        ylabel="Periapsis Longitude (deg)",
        label=false,
        linewidth=2.5,
        marker=:diamond,
        color=:firebrick
    )
    plot!(
        p,
        [NaN],
        [NaN];
        label="Periapsis Longitude",
        linewidth=2.5,
        marker=:diamond,
        color=:firebrick
    )
    if nrow(spice_peri_events) > 0
        plot!(
            lon_axis,
            spice_peri_events.orbit,
            spice_peri_events.longitude_deg;
            label=false,
            linewidth=2.0,
            linestyle=:dash,
            marker=:xcross,
            color=:darkred
        )
        plot!(
            p,
            [NaN],
            [NaN];
            label="SPICE Periapsis Longitude",
            linewidth=2.0,
            linestyle=:dash,
            marker=:xcross,
            color=:darkred
        )
    end

    plot_path = joinpath(results_dir, "periapsis_latitude_longitude_vs_orbit.png")
    savefig(p, plot_path)
    println("Saved periapsis latitude/longitude plot to $(abspath(plot_path))")
    return plot_path
end

function _save_apoapsis_radius_plot(args::SimulationConfiguration, planet, odyssey_schedule)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    extrema = _derive_orbit_extrema_from_results(csv_path, planet)
    max_orbit_plotted = length(extrema.apo.orbit)
    raise_orbits, lower_orbits = _split_maneuver_orbits(odyssey_schedule, max_orbit_plotted)

    p = plot(
        extrema.apo.orbit,
        extrema.apo.radius_km;
        xlabel="Orbit Number",
        ylabel="Apoapsis Radius (km)",
        label="Apoapsis Radius",
        linewidth=2.5,
        marker=:diamond,
        color=:darkmagenta,
        grid=true,
        legend=:topright
    )
    if !isempty(raise_orbits)
        vline!(
            p,
            raise_orbits;
            color=:seagreen3,
            linestyle=:dash,
            linewidth=1.5,
            alpha=0.7,
            label="Periapsis Raise Maneuver"
        )
    end
    if !isempty(lower_orbits)
        vline!(
            p,
            lower_orbits;
            color=:darkorange2,
            linestyle=:dashdot,
            linewidth=1.5,
            alpha=0.7,
            label="Periapsis Lower Maneuver"
        )
    end

    plot_path = joinpath(results_dir, "apoapsis_radius_vs_orbit.png")
    savefig(p, plot_path)
    println("Saved apoapsis radius plot to $(abspath(plot_path))")
    return plot_path
end

function _save_drag_along_velocity_plot(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))

    df = CSV.read(csv_path, DataFrame)
    altitude_m = _require_float_column(df, :sc1_altitude)
    vel_x = _require_float_column(df, :sc1_vel_1)
    vel_y = _require_float_column(df, :sc1_vel_2)
    vel_z = _require_float_column(df, :sc1_vel_3)
    drag_x = _require_float_column(df, :sc1_drag_1)
    drag_y = _require_float_column(df, :sc1_drag_2)
    drag_z = _require_float_column(df, :sc1_drag_3)

    peri_orbit = Int[]
    peri_drag_mag_n = Float64[]
    orbit_number = 0
    @inbounds for i in 2:(length(altitude_m) - 1)
        if altitude_m[i] <= altitude_m[i - 1] && altitude_m[i] < altitude_m[i + 1]
            orbit_number += 1
            vel = SVector{3, Float64}(vel_x[i], vel_y[i], vel_z[i])
            drag = SVector{3, Float64}(drag_x[i], drag_y[i], drag_z[i])
            vel_mag = norm(vel)
            drag_projection = vel_mag > 0.0 ? -dot(drag, vel / vel_mag) : 0.0
            push!(peri_orbit, orbit_number)
            push!(peri_drag_mag_n, max(abs(drag_projection), eps(Float64)))
        end
    end

    isempty(peri_orbit) && throw(ArgumentError("No periapsis samples found for drag plot."))

    p = plot(
        peri_orbit,
        peri_drag_mag_n;
        xlabel="Orbit Number",
        ylabel="|Drag Along Velocity| (N)",
        yscale=:log10,
        label="Periapsis Drag",
        linewidth=2.0,
        marker=:circle,
        color=:teal,
        grid=true,
        legend=:topright
    )

    plot_path = joinpath(results_dir, "drag_along_velocity_at_periapsis_vs_orbit.png")
    savefig(p, plot_path)
    println("Saved periapsis drag log plot to $(abspath(plot_path))")
    return plot_path
end

@inline function _perpendicular_component_magnitude(force::SVector{3, Float64}, velocity::SVector{3, Float64})::Float64
    vel_mag = norm(velocity)
    vel_mag > 0.0 || return 0.0
    vel_hat = velocity / vel_mag
    return norm(force - dot(force, vel_hat) * vel_hat)
end

function _save_aero_sideways_components_plot(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))

    df = CSV.read(csv_path, DataFrame)
    altitude_m = _require_float_column(df, :sc1_altitude)
    vel_x = _require_float_column(df, :sc1_vel_1)
    vel_y = _require_float_column(df, :sc1_vel_2)
    vel_z = _require_float_column(df, :sc1_vel_3)
    drag_x = _require_float_column(df, :sc1_drag_1)
    drag_y = _require_float_column(df, :sc1_drag_2)
    drag_z = _require_float_column(df, :sc1_drag_3)
    lift_x = _require_float_column(df, :sc1_lift_1)
    lift_y = _require_float_column(df, :sc1_lift_2)
    lift_z = _require_float_column(df, :sc1_lift_3)
    cross_x = _require_float_column(df, :sc1_cross_1)
    cross_y = _require_float_column(df, :sc1_cross_2)
    cross_z = _require_float_column(df, :sc1_cross_3)

    peri_orbit = Int[]
    peri_drag_sideways_n = Float64[]
    peri_lift_sideways_n = Float64[]
    peri_cross_sideways_n = Float64[]
    orbit_number = 0
    @inbounds for i in 2:(length(altitude_m) - 1)
        if altitude_m[i] <= altitude_m[i - 1] && altitude_m[i] < altitude_m[i + 1]
            orbit_number += 1
            velocity = SVector{3, Float64}(vel_x[i], vel_y[i], vel_z[i])
            drag = SVector{3, Float64}(drag_x[i], drag_y[i], drag_z[i])
            lift = SVector{3, Float64}(lift_x[i], lift_y[i], lift_z[i])
            cross_force = SVector{3, Float64}(cross_x[i], cross_y[i], cross_z[i])
            push!(peri_orbit, orbit_number)
            push!(peri_drag_sideways_n, max(_perpendicular_component_magnitude(drag, velocity), eps(Float64)))
            push!(peri_lift_sideways_n, max(_perpendicular_component_magnitude(lift, velocity), eps(Float64)))
            push!(peri_cross_sideways_n, max(_perpendicular_component_magnitude(cross_force, velocity), eps(Float64)))
        end
    end

    isempty(peri_orbit) && throw(ArgumentError("No periapsis samples found for aero sideways-component plot."))

    p = plot(
        peri_orbit,
        peri_drag_sideways_n;
        xlabel="Orbit Number",
        ylabel="Perpendicular Force Magnitude (N)",
        yscale=:log10,
        label="Drag Perpendicular",
        linewidth=2.0,
        marker=:circle,
        color=:teal,
        grid=true,
        legend=:topright
    )
    plot!(
        p,
        peri_orbit,
        peri_lift_sideways_n;
        label="Lift Perpendicular",
        linewidth=2.0,
        marker=:diamond,
        color=:firebrick
    )
    plot!(
        p,
        peri_orbit,
        peri_cross_sideways_n;
        label="Cross Perpendicular",
        linewidth=2.0,
        marker=:utriangle,
        color=:darkmagenta
    )

    plot_path = joinpath(results_dir, "aero_sideways_components_at_periapsis_vs_orbit.png")
    savefig(p, plot_path)
    println("Saved periapsis aero sideways-component plot to $(abspath(plot_path))")
    return plot_path
end

function _simulation_position_samples(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))

    df = CSV.read(csv_path, DataFrame)
    time_s = _require_float_column(df, :time)
    sim_x_km = _require_float_column(df, :sc1_pos_1) ./ 1e3
    sim_y_km = _require_float_column(df, :sc1_pos_2) ./ 1e3
    sim_z_km = _require_float_column(df, :sc1_pos_3) ./ 1e3
    return time_s, sim_x_km, sim_y_km, sim_z_km
end

function _simulation_velocity_samples(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))

    df = CSV.read(csv_path, DataFrame)
    time_s = _require_float_column(df, :time)
    sim_vx_kms = _require_float_column(df, :sc1_vel_1) ./ 1e3
    sim_vy_kms = _require_float_column(df, :sc1_vel_2) ./ 1e3
    sim_vz_kms = _require_float_column(df, :sc1_vel_3) ./ 1e3
    return time_s, sim_vx_kms, sim_vy_kms, sim_vz_kms
end

function _orbital_elements_from_state_samples(
    x_km::AbstractVector{<:Real},
    y_km::AbstractVector{<:Real},
    z_km::AbstractVector{<:Real},
    vx_kms::AbstractVector{<:Real},
    vy_kms::AbstractVector{<:Real},
    vz_kms::AbstractVector{<:Real},
    planet
)
    n = length(x_km)
    length(y_km) == n && length(z_km) == n && length(vx_kms) == n && length(vy_kms) == n && length(vz_kms) == n ||
        throw(ArgumentError("Position and velocity sample vectors must have matching lengths."))

    semi_major_axis_km = Vector{Float64}(undef, n)
    eccentricity = Vector{Float64}(undef, n)
    inclination_deg = Vector{Float64}(undef, n)
    raan_deg = Vector{Float64}(undef, n)
    arg_periapsis_deg = Vector{Float64}(undef, n)
    true_anomaly_deg = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        r_m = SVector{3, Float64}(x_km[i], y_km[i], z_km[i]) .* 1e3
        v_mps = SVector{3, Float64}(vx_kms[i], vy_kms[i], vz_kms[i]) .* 1e3
        oe = SpaceAGORA.TelemetryVerification.rvtoorbitalelement(r_m, v_mps, planet)
        semi_major_axis_km[i] = oe[1] / 1e3
        eccentricity[i] = oe[2]
        inclination_deg[i] = rad2deg(oe[3])
        raan_deg[i] = mod(rad2deg(oe[4]), 360.0)
        arg_periapsis_deg[i] = mod(rad2deg(oe[5]), 360.0)
        true_anomaly_deg[i] = mod(rad2deg(oe[6]), 360.0)
    end

    return (
        semi_major_axis_km=semi_major_axis_km,
        eccentricity=eccentricity,
        inclination_deg=inclination_deg,
        raan_deg=raan_deg,
        arg_periapsis_deg=arg_periapsis_deg,
        true_anomaly_deg=true_anomaly_deg
    )
end

function _spice_relative_position_samples(time_s::AbstractVector{<:Real}, initial_time::InitialTime, spice_path::String)
    _furnish_mission_kernels(spice_path)
    initial_et = _initial_time_et(initial_time)
    spice_x_km = Vector{Float64}(undef, length(time_s))
    spice_y_km = Vector{Float64}(undef, length(time_s))
    spice_z_km = Vector{Float64}(undef, length(time_s))
    @inbounds for i in eachindex(time_s)
        rel_pos_km, _ = _spice_relative_state_km(initial_et + time_s[i])
        spice_x_km[i] = rel_pos_km[1]
        spice_y_km[i] = rel_pos_km[2]
        spice_z_km[i] = rel_pos_km[3]
    end
    return spice_x_km, spice_y_km, spice_z_km
end

function _spice_relative_velocity_samples(time_s::AbstractVector{<:Real}, initial_time::InitialTime, spice_path::String)
    _furnish_mission_kernels(spice_path)
    initial_et = _initial_time_et(initial_time)
    spice_vx_kms = Vector{Float64}(undef, length(time_s))
    spice_vy_kms = Vector{Float64}(undef, length(time_s))
    spice_vz_kms = Vector{Float64}(undef, length(time_s))
    @inbounds for i in eachindex(time_s)
        _, rel_vel_kms = _spice_relative_state_km(initial_et + time_s[i])
        spice_vx_kms[i] = rel_vel_kms[1]
        spice_vy_kms[i] = rel_vel_kms[2]
        spice_vz_kms[i] = rel_vel_kms[3]
    end
    return spice_vx_kms, spice_vy_kms, spice_vz_kms
end

function _simulation_planet_fixed_state_samples(args::SimulationConfiguration, planet, initial_time::InitialTime)
    ephemerides_model = args.environment_model.ephemerides_model
    initial_et = SimulationModel.ephemerides_time_seconds(initial_time, ephemerides_model)
    time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    _, sim_vx_kms, sim_vy_kms, sim_vz_kms = _simulation_velocity_samples(args)
    n = length(time_s)
    fixed_x_km = Vector{Float64}(undef, n)
    fixed_y_km = Vector{Float64}(undef, n)
    fixed_z_km = Vector{Float64}(undef, n)
    fixed_vx_kms = Vector{Float64}(undef, n)
    fixed_vy_kms = Vector{Float64}(undef, n)
    fixed_vz_kms = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        et = initial_et + time_s[i]
        pos_m = SVector{3, Float64}(sim_x_km[i], sim_y_km[i], sim_z_km[i]) .* 1e3
        vel_mps = SVector{3, Float64}(sim_vx_kms[i], sim_vy_kms[i], sim_vz_kms[i]) .* 1e3
        pos_p_m, vel_p_mps = SimulationModel.SimulationCallbacks.r_intor_p!(pos_m, vel_mps, planet, et, ephemerides_model)
        fixed_x_km[i] = pos_p_m[1] / 1e3
        fixed_y_km[i] = pos_p_m[2] / 1e3
        fixed_z_km[i] = pos_p_m[3] / 1e3
        fixed_vx_kms[i] = vel_p_mps[1] / 1e3
        fixed_vy_kms[i] = vel_p_mps[2] / 1e3
        fixed_vz_kms[i] = vel_p_mps[3] / 1e3
    end
    return time_s, fixed_x_km, fixed_y_km, fixed_z_km, fixed_vx_kms, fixed_vy_kms, fixed_vz_kms
end

function _spice_planet_fixed_state_samples(
    time_s::AbstractVector{<:Real},
    planet,
    initial_time::InitialTime,
    spice_path::String
)
    _furnish_mission_kernels(spice_path)
    ephemerides_model = SpiceEphemeridesModel()
    initial_et = _initial_time_et(initial_time)
    n = length(time_s)
    fixed_x_km = Vector{Float64}(undef, n)
    fixed_y_km = Vector{Float64}(undef, n)
    fixed_z_km = Vector{Float64}(undef, n)
    fixed_vx_kms = Vector{Float64}(undef, n)
    fixed_vy_kms = Vector{Float64}(undef, n)
    fixed_vz_kms = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        et = initial_et + time_s[i]
        pos_km, vel_kms = _spice_relative_state_km(et)
        pos_m = pos_km .* 1e3
        vel_mps = vel_kms .* 1e3
        pos_p_m, vel_p_mps = SimulationModel.SimulationCallbacks.r_intor_p!(pos_m, vel_mps, planet, et, ephemerides_model)
        fixed_x_km[i] = pos_p_m[1] / 1e3
        fixed_y_km[i] = pos_p_m[2] / 1e3
        fixed_z_km[i] = pos_p_m[3] / 1e3
        fixed_vx_kms[i] = vel_p_mps[1] / 1e3
        fixed_vy_kms[i] = vel_p_mps[2] / 1e3
        fixed_vz_kms[i] = vel_p_mps[3] / 1e3
    end
    return fixed_x_km, fixed_y_km, fixed_z_km, fixed_vx_kms, fixed_vy_kms, fixed_vz_kms
end

function _local_extremum_times(time_s::AbstractVector{<:Real}, values::AbstractVector{<:Real}; kind::Symbol)
    times = Float64[]
    length(values) >= 3 || return times
    @inbounds for i in 2:(length(values) - 1)
        prev_val = values[i - 1]
        curr_val = values[i]
        next_val = values[i + 1]
        is_extremum = if kind === :min
            curr_val <= prev_val && curr_val < next_val
        elseif kind === :max
            curr_val >= prev_val && curr_val > next_val
        else
            throw(ArgumentError("Unsupported extremum kind: $kind"))
        end
        is_extremum && push!(times, Float64(time_s[i]))
    end
    return times
end

function _threshold_crossing_times(
    time_s::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
    threshold::Float64;
    direction::Symbol
)
    times = Float64[]
    length(values) >= 2 || return times
    @inbounds for i in 1:(length(values) - 1)
        t0 = Float64(time_s[i])
        t1 = Float64(time_s[i + 1])
        y0 = Float64(values[i]) - threshold
        y1 = Float64(values[i + 1]) - threshold
        t1 > t0 || continue
        if direction === :down
            crosses = y0 > 0.0 && y1 <= 0.0
        elseif direction === :up
            crosses = y0 <= 0.0 && y1 > 0.0
        else
            throw(ArgumentError("Unsupported crossing direction: $direction"))
        end
        crosses || continue
        denom = y1 - y0
        x = abs(denom) <= eps(Float64) ? 0.0 : -y0 / denom
        push!(times, t0 + clamp(x, 0.0, 1.0) * (t1 - t0))
    end
    return times
end

function _trajectory_marker_times(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    isfile(csv_path) || throw(ArgumentError("Simulation results CSV not found at $(abspath(csv_path))."))
    df = CSV.read(csv_path, DataFrame)
    time_s = _require_float_column(df, :time)
    altitude_m = _require_float_column(df, :sc1_altitude)

    peri_event_path = joinpath(results_dir, "periapsis_events.csv")
    periapsis_times_s = if isfile(peri_event_path)
        peri_events = CSV.read(peri_event_path, DataFrame)
        :time_s in Symbol.(names(peri_events)) ? Float64.(peri_events.time_s) : Float64[]
    else
        _local_extremum_times(time_s, altitude_m; kind=:min)
    end

    EI_m = args.environment_model.EI * 1e3
    return (
        apoapsis_s=_local_extremum_times(time_s, altitude_m; kind=:max),
        periapsis_s=periapsis_times_s,
        atmosphere_entry_s=_threshold_crossing_times(time_s, altitude_m, EI_m; direction=:down),
        atmosphere_exit_s=_threshold_crossing_times(time_s, altitude_m, EI_m; direction=:up)
    )
end

function _add_vertical_markers!(p, marker_times)
    if !isempty(marker_times.apoapsis_s)
        vline!(p, marker_times.apoapsis_s ./ 3600.0; label="Apoapsis", color=:purple, linestyle=:dash, linewidth=1.5, alpha=0.75)
    end
    if !isempty(marker_times.periapsis_s)
        vline!(p, marker_times.periapsis_s ./ 3600.0; label="Periapsis", color=:black, linestyle=:dot, linewidth=1.8, alpha=0.8)
    end
    if !isempty(marker_times.atmosphere_entry_s)
        vline!(p, marker_times.atmosphere_entry_s ./ 3600.0; label="Atmosphere Entry", color=:darkorange2, linestyle=:dashdot, linewidth=1.6, alpha=0.8)
    end
    if !isempty(marker_times.atmosphere_exit_s)
        vline!(p, marker_times.atmosphere_exit_s ./ 3600.0; label="Atmosphere Exit", color=:seagreen4, linestyle=:dashdot, linewidth=1.6, alpha=0.8)
    end
    return p
end

function _rtn_error_components(
    err_x::AbstractVector{<:Real},
    err_y::AbstractVector{<:Real},
    err_z::AbstractVector{<:Real},
    ref_x::AbstractVector{<:Real},
    ref_y::AbstractVector{<:Real},
    ref_z::AbstractVector{<:Real},
    ref_vx::AbstractVector{<:Real},
    ref_vy::AbstractVector{<:Real},
    ref_vz::AbstractVector{<:Real}
)
    n = length(err_x)
    err_r = Vector{Float64}(undef, n)
    err_t = Vector{Float64}(undef, n)
    err_n = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        r_ref = SVector{3, Float64}(ref_x[i], ref_y[i], ref_z[i])
        v_ref = SVector{3, Float64}(ref_vx[i], ref_vy[i], ref_vz[i])
        err = SVector{3, Float64}(err_x[i], err_y[i], err_z[i])
        r_mag = norm(r_ref)
        h_vec = cross(r_ref, v_ref)
        h_mag = norm(h_vec)
        if r_mag <= eps(Float64) || h_mag <= eps(Float64)
            throw(ArgumentError("Cannot construct RTN frame at sample $i: degenerate reference position/angular momentum."))
        end
        r_hat = r_ref / r_mag
        n_hat = h_vec / h_mag
        t_hat = cross(n_hat, r_hat)
        err_r[i] = max(abs(dot(err, r_hat)), eps(Float64))
        err_t[i] = max(abs(dot(err, t_hat)), eps(Float64))
        err_n[i] = max(abs(dot(err, n_hat)), eps(Float64))
    end
    return err_r, err_t, err_n
end

function _save_trajectory_comparison_plot(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String)
    config = _active_mission_spice_config()
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    spice_x_km, spice_y_km, spice_z_km = _spice_relative_position_samples(time_s, initial_time, spice_path)

    θ = range(0.0, 2pi; length=60)
    ϕ = range(0.0, pi; length=30)
    body_radius_km = planet.Rp_e / 1e3
    body_x = [body_radius_km * cos(θj) * sin(ϕi) for ϕi in ϕ, θj in θ]
    body_y = [body_radius_km * sin(θj) * sin(ϕi) for ϕi in ϕ, θj in θ]
    body_z = [body_radius_km * cos(ϕi) for ϕi in ϕ, θj in θ]

    traces = [
        surface(
            x=body_x,
            y=body_y,
            z=body_z,
            name=config.central_body,
            opacity=0.4,
            showscale=false,
            colorscale=[[0.0, config.surface_color], [1.0, config.surface_color]],
            hoverinfo="skip"
        ),
        scatter3d(
            x=sim_x_km,
            y=sim_y_km,
            z=sim_z_km,
            mode="lines",
            name="Simulation",
            line=attr(color="dodgerblue", width=5)
        ),
        scatter3d(
            x=spice_x_km,
            y=spice_y_km,
            z=spice_z_km,
            mode="lines",
            name="SPICE",
            line=attr(color="black", width=4)
        )
    ]
    layout = Layout(
        title="$(config.mission_name) Trajectory vs SPICE",
        scene=attr(
            xaxis=attr(title="x (km)", showgrid=true),
            yaxis=attr(title="y (km)", showgrid=true),
            zaxis=attr(title="z (km)", showgrid=true),
            aspectmode="data"
        ),
        legend=attr(x=0.02, y=0.98),
        margin=attr(l=0, r=0, b=0, t=40)
    )
    p = Plot(traces, layout)

    plot_path = joinpath(results_dir, "trajectory_3d_vs_spice.html")
    open(plot_path, "w") do io
        show(io, MIME("text/html"), p; include_plotlyjs="cdn", full_html=true)
    end
    println("Saved interactive 3D trajectory comparison plot to $(abspath(plot_path))")
    return plot_path
end

function _save_planet_fixed_trajectory_comparison_plot(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String)
    config = _active_mission_spice_config()
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km, _, _, _ = _simulation_planet_fixed_state_samples(args, planet, initial_time)
    spice_x_km, spice_y_km, spice_z_km, _, _, _ = _spice_planet_fixed_state_samples(time_s, planet, initial_time, spice_path)

    θ = range(0.0, 2pi; length=60)
    ϕ = range(0.0, pi; length=30)
    body_radius_km = planet.Rp_e / 1e3
    body_x = [body_radius_km * cos(θj) * sin(ϕi) for ϕi in ϕ, θj in θ]
    body_y = [body_radius_km * sin(θj) * sin(ϕi) for ϕi in ϕ, θj in θ]
    body_z = [body_radius_km * cos(ϕi) for ϕi in ϕ, θj in θ]

    traces = [
        surface(
            x=body_x,
            y=body_y,
            z=body_z,
            name=config.central_body,
            opacity=0.4,
            showscale=false,
            colorscale=[[0.0, config.surface_color], [1.0, config.surface_color]],
            hoverinfo="skip"
        ),
        scatter3d(
            x=sim_x_km,
            y=sim_y_km,
            z=sim_z_km,
            mode="lines",
            name="Simulation $(config.body_fixed_frame)",
            line=attr(color="dodgerblue", width=5)
        ),
        scatter3d(
            x=spice_x_km,
            y=spice_y_km,
            z=spice_z_km,
            mode="lines",
            name="SPICE $(config.body_fixed_frame)",
            line=attr(color="black", width=4)
        )
    ]
    layout = Layout(
        title="$(config.mission_name) Planet-Fixed Trajectory vs SPICE ($(config.body_fixed_frame))",
        scene=attr(
            xaxis=attr(title="$(config.body_fixed_frame) x (km)", showgrid=true),
            yaxis=attr(title="$(config.body_fixed_frame) y (km)", showgrid=true),
            zaxis=attr(title="$(config.body_fixed_frame) z (km)", showgrid=true),
            aspectmode="data"
        ),
        legend=attr(x=0.02, y=0.98),
        margin=attr(l=0, r=0, b=0, t=40)
    )
    p = Plot(traces, layout)

    plot_path = joinpath(results_dir, "planet_fixed_trajectory_3d_vs_spice.html")
    open(plot_path, "w") do io
        show(io, MIME("text/html"), p; include_plotlyjs="cdn", full_html=true)
    end
    println("Saved $(config.body_fixed_frame) 3D trajectory comparison plot to $(abspath(plot_path))")
    return plot_path
end

function _save_inertial_position_error_plot(args::SimulationConfiguration, initial_time::InitialTime, spice_path::String)
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    spice_x_km, spice_y_km, spice_z_km = _spice_relative_position_samples(time_s, initial_time, spice_path)

    err_x_km = max.(abs.(sim_x_km .- spice_x_km), eps(Float64))
    err_y_km = max.(abs.(sim_y_km .- spice_y_km), eps(Float64))
    err_z_km = max.(abs.(sim_z_km .- spice_z_km), eps(Float64))
    time_hr = time_s ./ 3600.0

    p = plot(
        time_hr,
        err_x_km;
        xlabel="Mission Time (hr)",
        ylabel="|Simulation - SPICE Position Error| (km)",
        yscale=:log10,
        label="|x error|",
        linewidth=2.0,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    plot!(
        p,
        time_hr,
        err_y_km;
        label="|y error|",
        linewidth=2.0,
        color=:firebrick
    )
    plot!(
        p,
        time_hr,
        err_z_km;
        label="|z error|",
        linewidth=2.0,
        color=:seagreen
    )
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "inertial_position_error_vs_spice.png")
    savefig(p, plot_path)
    println("Saved inertial position error plot to $(abspath(plot_path))")
    return plot_path
end

function _save_planet_fixed_position_error_plot(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String)
    config = _active_mission_spice_config()
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km, _, _, _ = _simulation_planet_fixed_state_samples(args, planet, initial_time)
    spice_x_km, spice_y_km, spice_z_km, _, _, _ = _spice_planet_fixed_state_samples(time_s, planet, initial_time, spice_path)

    err_x_km = max.(abs.(sim_x_km .- spice_x_km), eps(Float64))
    err_y_km = max.(abs.(sim_y_km .- spice_y_km), eps(Float64))
    err_z_km = max.(abs.(sim_z_km .- spice_z_km), eps(Float64))
    time_hr = time_s ./ 3600.0

    p = plot(
        time_hr,
        err_x_km;
        xlabel="Mission Time (hr)",
        ylabel="|Simulation - SPICE $(config.body_fixed_frame) Position Error| (km)",
        yscale=:log10,
        label="|$(config.body_fixed_frame) x error|",
        linewidth=2.0,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    plot!(p, time_hr, err_y_km; label="|$(config.body_fixed_frame) y error|", linewidth=2.0, color=:firebrick)
    plot!(p, time_hr, err_z_km; label="|$(config.body_fixed_frame) z error|", linewidth=2.0, color=:seagreen)
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "planet_fixed_position_error_vs_spice.png")
    savefig(p, plot_path)
    println("Saved $(config.body_fixed_frame) position error plot to $(abspath(plot_path))")
    return plot_path
end

function _save_planet_fixed_axis_trajectory_plot(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String)
    config = _active_mission_spice_config()
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km, _, _, _ = _simulation_planet_fixed_state_samples(args, planet, initial_time)
    spice_x_km, spice_y_km, spice_z_km, _, _, _ = _spice_planet_fixed_state_samples(time_s, planet, initial_time, spice_path)
    time_hr = time_s ./ 3600.0

    px = plot(
        time_hr,
        sim_x_km;
        ylabel="$(config.body_fixed_frame) x (km)",
        label="Simulation",
        linewidth=2.0,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    plot!(px, time_hr, spice_x_km; label="SPICE", linewidth=2.0, linestyle=:dash, color=:black)

    py = plot(
        time_hr,
        sim_y_km;
        ylabel="$(config.body_fixed_frame) y (km)",
        label="Simulation",
        linewidth=2.0,
        color=:firebrick,
        grid=true,
        legend=:topright
    )
    plot!(py, time_hr, spice_y_km; label="SPICE", linewidth=2.0, linestyle=:dash, color=:black)

    pz = plot(
        time_hr,
        sim_z_km;
        xlabel="Mission Time (hr)",
        ylabel="$(config.body_fixed_frame) z (km)",
        label="Simulation",
        linewidth=2.0,
        color=:seagreen,
        grid=true,
        legend=:topright
    )
    plot!(pz, time_hr, spice_z_km; label="SPICE", linewidth=2.0, linestyle=:dash, color=:black)

    p = plot(
        px,
        py,
        pz;
        layout=(3, 1),
        size=(1100, 900),
        plot_title="$(config.mission_name) $(config.body_fixed_frame) Axis Trajectories vs SPICE"
    )

    plot_path = joinpath(results_dir, "planet_fixed_axis_trajectories_vs_spice.png")
    savefig(p, plot_path)
    println("Saved $(config.body_fixed_frame) axis trajectory plot to $(abspath(plot_path))")
    return plot_path
end

function _save_orbital_elements_comparison_plot(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String)
    config = _active_mission_spice_config()
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    time_v_s, sim_vx_kms, sim_vy_kms, sim_vz_kms = _simulation_velocity_samples(args)
    length(time_s) == length(time_v_s) && all(time_s .== time_v_s) ||
        throw(ArgumentError("Position and velocity samples do not share the same time grid."))

    spice_x_km, spice_y_km, spice_z_km = _spice_relative_position_samples(time_s, initial_time, spice_path)
    spice_vx_kms, spice_vy_kms, spice_vz_kms = _spice_relative_velocity_samples(time_s, initial_time, spice_path)
    sim_oe = _orbital_elements_from_state_samples(
        sim_x_km,
        sim_y_km,
        sim_z_km,
        sim_vx_kms,
        sim_vy_kms,
        sim_vz_kms,
        planet
    )
    spice_oe = _orbital_elements_from_state_samples(
        spice_x_km,
        spice_y_km,
        spice_z_km,
        spice_vx_kms,
        spice_vy_kms,
        spice_vz_kms,
        planet
    )
    time_hr = time_s ./ 3600.0

    specs = (
        (:semi_major_axis_km, "Semi-major Axis (km)", :royalblue),
        (:eccentricity, "Eccentricity", :firebrick),
        (:inclination_deg, "Inclination (deg)", :seagreen),
        (:raan_deg, "RAAN (deg)", :darkorange2),
        (:arg_periapsis_deg, "Argument of Periapsis (deg)", :purple),
        (:true_anomaly_deg, "True Anomaly (deg)", :teal)
    )
    subplots = Plots.Plot[]
    for (field, ylabel, color) in specs
        sp = plot(
            time_hr,
            getproperty(sim_oe, field);
            label="Simulation",
            ylabel=ylabel,
            linewidth=2.0,
            color=color,
            grid=true,
            legend=:topright
        )
        plot!(
            sp,
            time_hr,
            getproperty(spice_oe, field);
            label="SPICE",
            linewidth=2.0,
            linestyle=:dash,
            color=:black
        )
        push!(subplots, sp)
    end
    xlabel!(subplots[5], "Mission Time (hr)")
    xlabel!(subplots[6], "Mission Time (hr)")

    p = plot(
        subplots...;
        layout=(3, 2),
        size=(1300, 1000),
        plot_title="$(config.mission_name) Orbital Elements vs SPICE"
    )
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "orbital_elements_vs_spice.png")
    savefig(p, plot_path)
    println("Saved orbital-elements comparison plot to $(abspath(plot_path))")
    return plot_path
end

function _save_inertial_velocity_error_plot(args::SimulationConfiguration, initial_time::InitialTime, spice_path::String)
    results_dir = args.simulation_settings.results_directory
    time_s, sim_vx_kms, sim_vy_kms, sim_vz_kms = _simulation_velocity_samples(args)
    spice_vx_kms, spice_vy_kms, spice_vz_kms = _spice_relative_velocity_samples(time_s, initial_time, spice_path)

    err_vx_kms = max.(abs.(sim_vx_kms .- spice_vx_kms), eps(Float64))
    err_vy_kms = max.(abs.(sim_vy_kms .- spice_vy_kms), eps(Float64))
    err_vz_kms = max.(abs.(sim_vz_kms .- spice_vz_kms), eps(Float64))
    time_hr = time_s ./ 3600.0

    p = plot(
        time_hr,
        err_vx_kms;
        xlabel="Mission Time (hr)",
        ylabel="|Simulation - SPICE Velocity Error| (km/s)",
        yscale=:log10,
        label="|vx error|",
        linewidth=2.0,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    plot!(
        p,
        time_hr,
        err_vy_kms;
        label="|vy error|",
        linewidth=2.0,
        color=:firebrick
    )
    plot!(
        p,
        time_hr,
        err_vz_kms;
        label="|vz error|",
        linewidth=2.0,
        color=:seagreen
    )
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "inertial_velocity_error_vs_spice.png")
    savefig(p, plot_path)
    println("Saved inertial velocity error plot to $(abspath(plot_path))")
    return plot_path
end

function _save_planet_fixed_velocity_error_plot(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String)
    config = _active_mission_spice_config()
    results_dir = args.simulation_settings.results_directory
    time_s, _, _, _, sim_vx_kms, sim_vy_kms, sim_vz_kms = _simulation_planet_fixed_state_samples(args, planet, initial_time)
    _, _, _, spice_vx_kms, spice_vy_kms, spice_vz_kms = _spice_planet_fixed_state_samples(time_s, planet, initial_time, spice_path)

    err_vx_kms = max.(abs.(sim_vx_kms .- spice_vx_kms), eps(Float64))
    err_vy_kms = max.(abs.(sim_vy_kms .- spice_vy_kms), eps(Float64))
    err_vz_kms = max.(abs.(sim_vz_kms .- spice_vz_kms), eps(Float64))
    time_hr = time_s ./ 3600.0

    p = plot(
        time_hr,
        err_vx_kms;
        xlabel="Mission Time (hr)",
        ylabel="|Simulation - SPICE $(config.body_fixed_frame) Velocity Error| (km/s)",
        yscale=:log10,
        label="|$(config.body_fixed_frame) vx error|",
        linewidth=2.0,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    plot!(p, time_hr, err_vy_kms; label="|$(config.body_fixed_frame) vy error|", linewidth=2.0, color=:firebrick)
    plot!(p, time_hr, err_vz_kms; label="|$(config.body_fixed_frame) vz error|", linewidth=2.0, color=:seagreen)
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "planet_fixed_velocity_error_vs_spice.png")
    savefig(p, plot_path)
    println("Saved $(config.body_fixed_frame) velocity error plot to $(abspath(plot_path))")
    return plot_path
end

function _save_rtn_position_error_plot(args::SimulationConfiguration, initial_time::InitialTime, spice_path::String)
    results_dir = args.simulation_settings.results_directory
    time_s, sim_x_km, sim_y_km, sim_z_km = _simulation_position_samples(args)
    spice_x_km, spice_y_km, spice_z_km = _spice_relative_position_samples(time_s, initial_time, spice_path)
    spice_vx_kms, spice_vy_kms, spice_vz_kms = _spice_relative_velocity_samples(time_s, initial_time, spice_path)

    err_r_km, err_t_km, err_n_km = _rtn_error_components(
        sim_x_km .- spice_x_km,
        sim_y_km .- spice_y_km,
        sim_z_km .- spice_z_km,
        spice_x_km,
        spice_y_km,
        spice_z_km,
        spice_vx_kms,
        spice_vy_kms,
        spice_vz_kms
    )
    time_hr = time_s ./ 3600.0

    p = plot(
        time_hr,
        err_r_km;
        xlabel="Mission Time (hr)",
        ylabel="|Simulation - SPICE Position Error| (km)",
        yscale=:log10,
        label="|R error|",
        linewidth=2.0,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    plot!(p, time_hr, err_t_km; label="|T error|", linewidth=2.0, color=:firebrick)
    plot!(p, time_hr, err_n_km; label="|N error|", linewidth=2.0, color=:seagreen)
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "rtn_position_error_vs_spice.png")
    savefig(p, plot_path)
    println("Saved RTN position error plot to $(abspath(plot_path))")
    return plot_path
end

function _save_rtn_velocity_error_plot(args::SimulationConfiguration, initial_time::InitialTime, spice_path::String)
    results_dir = args.simulation_settings.results_directory
    time_s, _, _, _ = _simulation_position_samples(args)
    _, sim_vx_kms, sim_vy_kms, sim_vz_kms = _simulation_velocity_samples(args)
    spice_x_km, spice_y_km, spice_z_km = _spice_relative_position_samples(time_s, initial_time, spice_path)
    spice_vx_kms, spice_vy_kms, spice_vz_kms = _spice_relative_velocity_samples(time_s, initial_time, spice_path)

    err_r_kms, err_t_kms, err_n_kms = _rtn_error_components(
        sim_vx_kms .- spice_vx_kms,
        sim_vy_kms .- spice_vy_kms,
        sim_vz_kms .- spice_vz_kms,
        spice_x_km,
        spice_y_km,
        spice_z_km,
        spice_vx_kms,
        spice_vy_kms,
        spice_vz_kms
    )
    time_hr = time_s ./ 3600.0

    p = plot(
        time_hr,
        err_r_kms;
        xlabel="Mission Time (hr)",
        ylabel="|Simulation - SPICE Velocity Error| (km/s)",
        yscale=:log10,
        label="|R error|",
        linewidth=2.0,
        color=:royalblue,
        grid=true,
        legend=:topright
    )
    plot!(p, time_hr, err_t_kms; label="|T error|", linewidth=2.0, color=:firebrick)
    plot!(p, time_hr, err_n_kms; label="|N error|", linewidth=2.0, color=:seagreen)
    _add_vertical_markers!(p, _trajectory_marker_times(args))

    plot_path = joinpath(results_dir, "rtn_velocity_error_vs_spice.png")
    savefig(p, plot_path)
    println("Saved RTN velocity error plot to $(abspath(plot_path))")
    return plot_path
end

function _open_plot_in_browser(plot_path::String)
    try
        run(`xdg-open $plot_path`)
    catch err
        @warn "Failed to open interactive trajectory plot automatically." plot_path exception=(err, catch_backtrace())
    end
    return nothing
end
