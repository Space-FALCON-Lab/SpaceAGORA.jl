include(joinpath(@__DIR__, "common.jl"))
include(joinpath(REPO_ROOT, "src", "mission", "operations", "maneuver_plans.jl"))
using CSV
using DataFrames
using Plots
using PlotlyJS: Plot, scatter3d, surface, Layout, attr
using Roots
using SPICE
using StaticArrays
using LinearAlgebra
using Logging

include(joinpath(@__DIR__, "aerobraking_mission_plot_utils.jl"))

const ORIGINAL_MISSION_TIME_S = 3_600.0 * 900.0
const DEFAULT_RESET_DELAY_S = 600.0
const BURN_STOP_EPS_S = 1.0

function _parse_positive_env_float(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    parsed = tryparse(Float64, raw)
    parsed !== nothing && parsed > 0.0 || throw(ArgumentError("$name must be > 0.0, got '$raw'."))
    return parsed
end

function _parse_positive_env_int(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    parsed = tryparse(Int, raw)
    parsed !== nothing && parsed > 0 || throw(ArgumentError("$name must be a positive integer, got '$raw'."))
    return parsed
end

function _mission_time_config(args::SimulationConfiguration, mission_time_s::Float64)::SimulationConfiguration
    mc = args.mission_configuration
    mission_cfg = MissionConfiguration(
        mission_type=MissionTime,
        keplerian=mc.keplerian,
        number_of_orbits=mc.number_of_orbits,
        mission_time=mission_time_s,
        orientation_sim=mc.orientation_sim,
        num_steps_to_save=mc.num_steps_to_save,
        data_rate=mc.data_rate
    )
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=args.simulation_settings,
        mission_configuration=mission_cfg,
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

function _runtime_settings(
    base_settings::SimulationSettings,
    results_dir::String,
    checkpoint_dir::String,
    checkpoint_interval_s::Float64,
    resume::Bool
)::SimulationSettings
    return SimulationSettings(
        results=true,
        verbose=base_settings.verbose,
        results_directory=results_dir,
        generate_plots=false,
        generate_filenames=false,
        normalize=false,
        save_csv=true,
        checkpoint_enabled=true,
        checkpoint_interval_s=max(1.0, checkpoint_interval_s),
        checkpoint_directory=checkpoint_dir,
        resume_from_checkpoint=resume
    )
end

function _with_runtime_settings(args::SimulationConfiguration, settings::SimulationSettings)::SimulationConfiguration
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=settings,
        mission_configuration=args.mission_configuration,
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

function _apoapsis_timing_from_state_s(u, planet)
    pos = SVector{3, Float64}(u.sc[1].pos)
    vel = SVector{3, Float64}(u.sc[1].vel)
    a, e, _, _, _, ν = SimulationModel.ControlHooks.rvtoorbitalelement(pos, vel, planet)
    if !(isfinite(a) && isfinite(e) && isfinite(ν) && a > 0.0 && 0.0 <= e < 1.0)
        throw(ArgumentError("Cannot estimate next apoapsis from non-elliptic state: a=$a, e=$e, ν=$ν."))
    end
    n = sqrt(planet.μ / a^3)
    E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(ν / 2.0))
    M = mod(E - e * sin(E), 2π)
    dt = mod(π - M, 2π) / n
    period = 2π / n
    return (time_to_next_apoapsis_s=(dt <= 1e-6 ? period : dt), period_s=period)
end

function _time_to_correction_apoapsis_s(u, planet, orbits_between_corrections::Int)::Float64
    orbits_between_corrections > 0 || throw(ArgumentError("orbits_between_corrections must be > 0."))
    timing = _apoapsis_timing_from_state_s(u, planet)
    return timing.time_to_next_apoapsis_s + (orbits_between_corrections - 1) * timing.period_s
end

function _checkpoint_active_burn_stop(ckpt)::Float64
    runtime_state = hasproperty(ckpt, :runtime_state) ? ckpt.runtime_state : nothing
    runtime_state === nothing && return NaN
    hasproperty(runtime_state, :control_effectors) || return NaN
    stop_max = NaN
    for control_snapshot in runtime_state.control_effectors
        control_snapshot === nothing && continue
        hasproperty(control_snapshot, :start_burn_time) || continue
        hasproperty(control_snapshot, :stop_burn_time) || continue
        for (start_t, stop_t) in zip(control_snapshot.start_burn_time, control_snapshot.stop_burn_time)
            if isfinite(start_t) && isfinite(stop_t) && stop_t > start_t && ckpt.t <= stop_t + 1e-9
                stop_max = isfinite(stop_max) ? max(stop_max, Float64(stop_t)) : Float64(stop_t)
            end
        end
    end
    return stop_max
end

function _rewrite_checkpoint_state_compat!(args::SimulationConfiguration, t::Float64, u_state)
    if isdefined(SimulationEngine, :_rewrite_checkpoint_state!)
        return Base.invokelatest(SimulationEngine._rewrite_checkpoint_state!, args, t, u_state)
    end
    @warn "Active SimulationEngine does not expose runtime-aware checkpoint rewriting; rewriting state-only checkpoint. Restart Julia or reload SpaceAGORA to enable runtime checkpoint state."
    return SimulationModel.IOSerialization._write_checkpoint!(args, t, u_state, "1")
end

function _reset_checkpoint_position_velocity_from_spice!(
    args::SimulationConfiguration,
    initial_time::InitialTime,
    spice_path::String
)
    ckpt = SimulationEngine._load_checkpoint(args)
    ckpt === nothing && throw(ArgumentError("No checkpoint was written for reset."))
    reset_u = deepcopy(ckpt.u)
    et = _initial_time_et(initial_time) + ckpt.t
    pos_km, vel_kms = _spice_relative_state_km(et, MARS_ODYSSEY_SPICE_CONFIG)
    reset_u.sc[1].pos .= pos_km .* 1e3
    reset_u.sc[1].vel .= vel_kms .* 1e3
    _rewrite_checkpoint_state_compat!(args, ckpt.t, reset_u)
    return reset_u
end

function _replace_last_result_sample_state!(
    csv_path::String,
    args::SimulationConfiguration,
    reset_u,
    reset_t::Float64
)
    isfile(csv_path) || throw(ArgumentError("Cannot update reset sample; result CSV not found at $(abspath(csv_path))."))
    df = CSV.read(csv_path, DataFrame)
    nrow(df) > 0 || return csv_path
    row_idx = nrow(df)
    if abs(Float64(df.time[row_idx]) - reset_t) > 1e-6
        @warn "Skipping reset-sample CSV correction because the last saved sample is not the reset time." csv_path last_time=Float64(df.time[row_idx]) reset_t
        return csv_path
    end

    pos = SVector{3, Float64}(reset_u.sc[1].pos)
    vel = SVector{3, Float64}(reset_u.sc[1].vel)
    df.sc1_pos_1[row_idx] = pos[1]
    df.sc1_pos_2[row_idx] = pos[2]
    df.sc1_pos_3[row_idx] = pos[3]
    df.sc1_vel_1[row_idx] = vel[1]
    df.sc1_vel_2[row_idx] = vel[2]
    df.sc1_vel_3[row_idx] = vel[3]

    planet = args.environment_model.planet
    ephemerides_model = args.environment_model.ephemerides_model
    et = _initial_time_et(args.initial_time) + reset_t
    rp, _ = SimulationModel.SimulationCallbacks.r_intor_p!(pos, vel, planet, et, ephemerides_model)
    altitude_m, lat_rad, lon_rad = SimulationModel.SimulationCallbacks.rtolatlong(rp, planet, ephemerides_model)
    if :sc1_altitude in propertynames(df)
        df.sc1_altitude[row_idx] = altitude_m
    end
    if :sc1_latitude_deg in propertynames(df)
        df.sc1_latitude_deg[row_idx] = rad2deg(lat_rad)
    end
    if :sc1_longitude_deg in propertynames(df)
        df.sc1_longitude_deg[row_idx] = rad2deg(lon_rad)
    end
    if :sc1_periapsis_altitude in propertynames(df)
        a, e, _, _, _, _ = SimulationModel.ControlHooks.rvtoorbitalelement(pos, vel, planet)
        df.sc1_periapsis_altitude[row_idx] = a * (1.0 - e) - planet.Rp_e
    end

    CSV.write(csv_path, df)
    return csv_path
end

function _append_part_csv!(parts::Vector{String}, path::String)
    isfile(path) || throw(ArgumentError("Expected leg result CSV at $(abspath(path))."))
    push!(parts, path)
    return path
end

function _stitch_result_parts(parts::Vector{String}, final_csv::String)
    isempty(parts) && throw(ArgumentError("No result parts were generated."))
    stitched = DataFrame()
    last_t = -Inf
    for path in parts
        df = CSV.read(path, DataFrame)
        if nrow(df) == 0
            continue
        end
        df = df[Float64.(df.time) .> last_t .+ 1e-9, :]
        if nrow(df) == 0
            continue
        end
        append!(stitched, df; cols=:union)
        last_t = Float64(stitched.time[end])
    end
    mkpath(dirname(final_csv))
    CSV.write(final_csv, stitched)
    return stitched
end

function _save_stitched_periapsis_events_from_results(args::SimulationConfiguration)
    results_dir = args.simulation_settings.results_directory
    csv_path = joinpath(results_dir, "simulation_results.csv")
    df = CSV.read(csv_path, DataFrame)
    event_path = joinpath(results_dir, "periapsis_events.csv")
    if nrow(df) < 3
        CSV.write(event_path, DataFrame(
            orbit=Int[],
            time_s=Float64[],
            altitude_km=Float64[],
            latitude_deg=Float64[],
            longitude_deg=Float64[],
            radius_km=Float64[],
            radial_velocity_mps=Float64[]
        ))
        return event_path
    end

    time_s = Float64.(df.time)
    pos_x = Float64.(df.sc1_pos_1)
    pos_y = Float64.(df.sc1_pos_2)
    pos_z = Float64.(df.sc1_pos_3)
    vel_x = Float64.(df.sc1_vel_1)
    vel_y = Float64.(df.sc1_vel_2)
    vel_z = Float64.(df.sc1_vel_3)
    altitude_km = Float64.(df.sc1_altitude) ./ 1e3
    latitude_deg = Float64.(df.sc1_latitude_deg)
    longitude_deg = Float64.(df.sc1_longitude_deg)
    radius_km = sqrt.(pos_x .^ 2 .+ pos_y .^ 2 .+ pos_z .^ 2) ./ 1e3
    correction_times = _correction_reset_times_for_results(csv_path)
    radial_velocity = similar(time_s)
    for idx in eachindex(time_s)
        pos = SVector{3, Float64}(pos_x[idx], pos_y[idx], pos_z[idx])
        vel = SVector{3, Float64}(vel_x[idx], vel_y[idx], vel_z[idx])
        radial_velocity[idx] = dot(pos, vel) / norm(pos)
    end

    event_time_s = Float64[]
    event_altitude_km = Float64[]
    event_latitude_deg = Float64[]
    event_longitude_deg = Float64[]
    event_radius_km = Float64[]
    event_radial_velocity_mps = Float64[]
    for idx in 1:(length(time_s) - 1)
        _interval_touches_correction(time_s[idx], time_s[idx + 1], correction_times) && continue
        rv0 = radial_velocity[idx]
        rv1 = radial_velocity[idx + 1]
        if !(isfinite(rv0) && isfinite(rv1) && rv0 < 0.0 && rv1 >= 0.0)
            continue
        end
        frac = rv0 / (rv0 - rv1)
        push!(event_time_s, time_s[idx] + frac * (time_s[idx + 1] - time_s[idx]))
        push!(event_altitude_km, altitude_km[idx] + frac * (altitude_km[idx + 1] - altitude_km[idx]))
        push!(event_latitude_deg, latitude_deg[idx] + frac * (latitude_deg[idx + 1] - latitude_deg[idx]))
        push!(event_longitude_deg, longitude_deg[idx] + frac * (longitude_deg[idx + 1] - longitude_deg[idx]))
        push!(event_radius_km, radius_km[idx] + frac * (radius_km[idx + 1] - radius_km[idx]))
        push!(event_radial_velocity_mps, 0.0)
    end

    events = DataFrame(
        orbit=collect(1:length(event_time_s)),
        time_s=event_time_s,
        altitude_km=event_altitude_km,
        latitude_deg=event_latitude_deg,
        longitude_deg=event_longitude_deg,
        radius_km=event_radius_km,
        radial_velocity_mps=event_radial_velocity_mps
    )
    CSV.write(event_path, events)
    println("Saved stitched periapsis event table to $(abspath(event_path))")
    return event_path
end

function _save_final_plots(args::SimulationConfiguration, planet, initial_time::InitialTime, spice_path::String, odyssey_schedule, smoke_mode::Bool)
    _save_stitched_periapsis_events_from_results(args)
    spice_periapsis_event_path = _save_spice_periapsis_event_table(args, planet, initial_time, spice_path; allow_empty=smoke_mode)
    spice_apoapsis_event_path = _save_spice_apoapsis_event_table(args, planet, initial_time, spice_path; allow_empty=smoke_mode)
    if _has_orbit_extrema_for_plots(args, planet)
        _save_apoapsis_periapsis_plot(args, planet, odyssey_schedule)
        _save_apoapsis_radius_plot(args, planet, odyssey_schedule)
        try
            _save_periapsis_latlon_plot(args, planet, odyssey_schedule)
        catch err
            if err isa ArgumentError
                @warn "Skipping final periapsis latitude/longitude plot; no event-located simulation periapsis table was generated." reason=sprint(showerror, err)
            else
                rethrow()
            end
        end
        _save_drag_along_velocity_plot(args)
        _save_aero_sideways_components_plot(args)
    end
    _save_trajectory_comparison_plot(args, planet, initial_time, spice_path)
    _save_planet_fixed_trajectory_comparison_plot(args, planet, initial_time, spice_path)
    _save_planet_fixed_axis_trajectory_plot(args, planet, initial_time, spice_path)
    _save_orbital_elements_comparison_plot(args, planet, initial_time, spice_path)
    _save_inertial_position_error_plot(args, initial_time, spice_path)
    _save_inertial_velocity_error_plot(args, initial_time, spice_path)
    _save_planet_fixed_position_error_plot(args, planet, initial_time, spice_path)
    _save_planet_fixed_velocity_error_plot(args, planet, initial_time, spice_path)
    _save_rtn_position_error_plot(args, initial_time, spice_path)
    _save_rtn_velocity_error_plot(args, initial_time, spice_path)
    return (spice_periapsis_event_path=spice_periapsis_event_path, spice_apoapsis_event_path=spice_apoapsis_event_path)
end

planet = Mars("", SPICE_PATH)
smoke_mode = get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"
num_legs = _parse_positive_env_int("SPACEAGORA_ODYSSEY_CHECKPOINT_ORBITS", smoke_mode ? 2 : 50)
orbits_between_corrections = _parse_positive_env_int("SPACEAGORA_ODYSSEY_ORBITS_BETWEEN_CORRECTIONS", 5)
reset_delay_s = _parse_positive_env_float("SPACEAGORA_ODYSSEY_RESET_DELAY_S", DEFAULT_RESET_DELAY_S)
initial_time = InitialTime(year=2001, month=11, day=6, hour=10, minute=5, second=12.7)

ic = _mars_odyssey_initial_condition_from_spice(initial_time, SPICE_PATH)

spacecraft = make_three_body_spacecraft(
    bus_dims=(2.2, 2.6, 1.7),
    panel_dims=(0.01, 5.5 / 1.35, 2.6),
    bus_mass=391.0,
    panel_mass_each=10.0,
    panel_offset_y=2.6 / 2.0 + 5.5 / 4.0,
    ic=ic,
    reflection_coefficient=0.9,
    prop_mass=50.0,
    id=100
)

mars_harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "GMM3.csv")
sun_gravity = NBodyGravityModel(body_names=("Sun",), primary_body_name="Mars", planet=planet)
solar_radiation_pressure = SolarRadiationPressureModel(spacecraft.root.reflection_coefficient, spacecraft.root.ref_area)
dynamic_effectors = smoke_mode ? (
    InverseSquaredGravityModel(),
    sun_gravity
) : (
    InverseSquaredGravityModel(),
    sun_gravity,
    GravitationalHarmonicsModel(50, 50, mars_harmonics_file, planet),
    solar_radiation_pressure,
    AerodynamicCoefficientfM()
)
density_model = GRAMAtmosphereModel(planet_name="mars")
root_output_dir = joinpath(REPO_ROOT, "output", "odyssey_orbit_checkpoint")
checkpoint_dir = joinpath(root_output_dir, "checkpoints")
final_results_dir = joinpath(root_output_dir, "final")

base_args = make_example_config(
    planet=planet,
    spacecraft=spacecraft,
    mission_time=ORIGINAL_MISSION_TIME_S,
    initial_time=initial_time,
    dynamic_effectors=dynamic_effectors,
    density_model=density_model,
    orientation_sim=false,
    keplerian=smoke_mode,
    EI_km=160.0,
    results_directory=final_results_dir
)

odyssey_schedule_base = odyssey_campaign_maneuvers(1:311; planet=planet)
odyssey_schedule = (
    maneuver_orbit_number=odyssey_schedule_base.maneuver_orbit_number .+ 1,
    maneuver_Δv=odyssey_schedule_base.maneuver_Δv
)

thruster = BaseThrusterModel(
    thrust=[4.0],
    direction=[0.0],
    Δv=[0.0],
    start_burn_time=[-1.0],
    stop_burn_time=[-1.0],
    Isp=[300.0]
)
guidance_effector = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
    maneuver_orbit_number=odyssey_schedule.maneuver_orbit_number,
    maneuver_Δv=odyssey_schedule.maneuver_Δv
)
args = SimulationConfiguration(
    file_paths=base_args.file_paths,
    simulation_settings=base_args.simulation_settings,
    mission_configuration=base_args.mission_configuration,
    environment_model=base_args.environment_model,
    dynamics_model=base_args.dynamics_model,
    guidance_model=GuidanceModel(guidance_effectors=(guidance_effector,), guidance_rates=[30.0]),
    navigation_model=base_args.navigation_model,
    control_model=ControlModel(control_effectors=(thruster,), control_rates=[10.0]),
    initial_time=base_args.initial_time,
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-8,
        abstol_orbit=1e-8,
        dt_max_orbit=30.0,
        reltol_atmosphere=1e-8,
        abstol_atmosphere=1e-8,
        dt_max_atmosphere=1.0
    )
)

mkpath(root_output_dir)
mkpath(checkpoint_dir)
clear_settings = _runtime_settings(args.simulation_settings, final_results_dir, checkpoint_dir, 1.0, false)
SimulationEngine._clear_checkpoint!(_with_runtime_settings(args, clear_settings))

result_parts = String[]
leg_rows = NamedTuple[]
current_t = 0.0
current_u = SimulationEngine.build_initial_conditions(args)

sim_elapsed = @elapsed begin
    for leg_idx in 1:num_legs
        global current_t, current_u
        time_to_apo = _time_to_correction_apoapsis_s(current_u, planet, orbits_between_corrections)
        apoapsis_t = current_t + time_to_apo
        planned_correction_t = apoapsis_t + reset_delay_s
        target_t = min(planned_correction_t, ORIGINAL_MISSION_TIME_S)
        leg_start_t = current_t
        leg_start_mass = Float64(current_u.sc[1].mass)
        part_idx = 1
        resume = leg_idx > 1
        final_part_args = nothing
        burn_stop_t = NaN

        while true
            part_dir = joinpath(root_output_dir, "leg_$(lpad(leg_idx, 3, '0'))_part_$(lpad(part_idx, 3, '0'))")
            part_settings = _runtime_settings(
                args.simulation_settings,
                part_dir,
                checkpoint_dir,
                max(1.0, target_t - current_t),
                resume
            )
            part_args = _mission_time_config(_with_runtime_settings(args, part_settings), target_t)
            final_part_args = part_args
            part_elapsed = @elapsed run_simulation(part_args; return_solution=false)
            csv_path = joinpath(part_dir, "simulation_results.csv")
            _append_part_csv!(result_parts, csv_path)

            ckpt = SimulationEngine._load_checkpoint(part_args)
            ckpt === nothing && throw(ArgumentError("No checkpoint found after leg $leg_idx part $part_idx."))
            burn_stop_t = _checkpoint_active_burn_stop(ckpt)
            if isfinite(burn_stop_t) && burn_stop_t + BURN_STOP_EPS_S > target_t && burn_stop_t < ORIGINAL_MISSION_TIME_S
                target_t = min(burn_stop_t + BURN_STOP_EPS_S, ORIGINAL_MISSION_TIME_S)
                current_t = ckpt.t
                resume = true
                part_idx += 1
                println("Leg $leg_idx part $(part_idx - 1) reached t=$(round(ckpt.t, digits=3)) s in $(round(part_elapsed, digits=3)) s; extending to finish active burn at t=$(round(target_t, digits=3)) s.")
                continue
            end

            current_t = ckpt.t
            correction_applied = current_t >= planned_correction_t - 1e-6
            if correction_applied
                current_u = _reset_checkpoint_position_velocity_from_spice!(part_args, initial_time, SPICE_PATH)
                _replace_last_result_sample_state!(csv_path, part_args, current_u, current_t)
            else
                current_u = deepcopy(ckpt.u)
            end
            leg_end_mass = Float64(current_u.sc[1].mass)
            push!(
                leg_rows,
                (
                    leg=leg_idx,
                    orbits_between_corrections=orbits_between_corrections,
                    correction_applied=correction_applied,
                    start_time_s=leg_start_t,
                    apoapsis_time_s=apoapsis_t,
                    reset_time_s=current_t,
                    reset_delay_s=current_t - apoapsis_t,
                    start_mass_kg=leg_start_mass,
                    end_mass_kg=leg_end_mass,
                    burn_stop_time_s=burn_stop_t,
                    parts=part_idx,
                    checkpoint_data=SimulationEngine._checkpoint_paths(part_args).data
                )
            )
            println(
                "Leg $leg_idx complete: start=$(round(leg_start_t, digits=3)) s, " *
                "correction_apoapsis=$(round(apoapsis_t, digits=3)) s, reset=$(round(current_t, digits=3)) s, " *
                "orbits_between_corrections=$orbits_between_corrections, " *
                "correction_applied=$correction_applied, " *
                "mass=$(round(leg_start_mass, digits=6)) -> $(round(leg_end_mass, digits=6)) kg, parts=$part_idx"
            )
            break
        end

        if current_t >= ORIGINAL_MISSION_TIME_S
            println("Reached original Odyssey mission time; stopping after leg $leg_idx.")
            break
        end
    end
end

final_csv_path = joinpath(final_results_dir, "simulation_results.csv")
stitched = _stitch_result_parts(result_parts, final_csv_path)
summary_path = joinpath(root_output_dir, "leg_summary.csv")
CSV.write(summary_path, DataFrame(leg_rows))
println("Saved stitched checkpoint-restart results ($(nrow(stitched)) samples) to $(abspath(final_csv_path))")
println("Saved leg summary to $(abspath(summary_path))")
println("COMPUTATIONAL TIME = $(sim_elapsed) s")

final_args = _with_runtime_settings(
    _mission_time_config(args, current_t),
    SimulationSettings(
        results=true,
        verbose=false,
        results_directory=final_results_dir,
        generate_plots=false,
        generate_filenames=false,
        normalize=false,
        save_csv=true
    )
)
_save_final_plots(final_args, final_args.environment_model.planet, initial_time, SPICE_PATH, odyssey_schedule, smoke_mode)
