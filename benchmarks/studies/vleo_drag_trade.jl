const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_PATH = joinpath(REPO_ROOT, "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
const DEFAULT_OUT_DIR = joinpath(REPO_ROOT, "output", "vleo_drag_trade")

using CSV
using DataFrames
using Dates
using Distributed
using LinearAlgebra
using Printf
using Statistics
using StaticArrays
import GRAMSuite

if !isdefined(@__MODULE__, :SpaceAGORA)
    @eval using SpaceAGORA
end
if myid() == 1 && !isdefined(@__MODULE__, :Plots)
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
    @eval using Plots
end

const SM = SpaceAGORA.SimulationModel
const TV = SpaceAGORA.TelemetryVerification

const VLEO_CHANNELS = (:low, :nominal, :high)
const VLEO_OUTER_ROUTES = (:auto, :none, :threads, :process)
const FULL_PERIGEES_KM = collect(150.0:10.0:300.0)
const FULL_APOGEES_KM = [350.0, 550.0, 800.0]
const REFERENCE_INITIAL_TIME = SM.InitialTime(year=2026, month=3, day=20, hour=12, minute=0, second=0.0)
const REFERENCE_INCLINATION_DEG = 53.0
const OUTER_ROUTE_STATE_FILENAME = "outer_route_state.toml"

const REFERENCE_VEHICLE = (
    name="Starlink-like mid-size reference",
    mass_total_kg=800.0,
    bus_mass_kg=760.0,
    panel_mass_each_kg=20.0,
    bus_dims_m=(2.6, 1.8, 0.35),
    panel_dims_m=(0.02, 4.0, 1.0),
    panel_offset_y_m=2.5,
)

const REFERENCE_SOLAR_CASE = (
    daily_f10=120.0,
    mean_f10=120.0,
    ap=10.0,
)

const CHANNEL_COLORS = Dict(
    "low" => "#1d4ed8",
    "nominal" => "#111827",
    "high" => "#b91c1c",
)

const _STUDY_PLANET_CACHE = Ref{Any}(nothing)
const _OUTER_ROUTE_WARNING_EMITTED = Ref(false)

function parse_cli_vleo(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || error("Unsupported argument '$arg'. Use --key=value or --flag.")
        body = arg[3:end]
        if occursin("=", body)
            k, v = split(body, "=", limit=2)
            opts[k] = v
        else
            if i < length(args) && !startswith(args[i + 1], "--")
                i += 1
                opts[body] = args[i]
            else
                opts[body] = "true"
            end
        end
        i += 1
    end
    return opts
end

@inline parse_bool_vleo(x::AbstractString) = lowercase(strip(String(x))) in ("1", "true", "yes", "on")

@inline function _parse_outer_route_vleo(value)::Symbol
    token = lowercase(strip(String(value)))
    route = Symbol(token)
    route in VLEO_OUTER_ROUTES || throw(ArgumentError(
        "Unsupported --outer-route='$value'. Use one of: auto, none, threads, process."
    ))
    return route
end

@inline function _parse_positive_int_vleo(name::AbstractString, value)::Int
    parsed = try
        parse(Int, String(value))
    catch
        throw(ArgumentError("$name must be an integer, got '$value'."))
    end
    parsed >= 1 || throw(ArgumentError("$name must be >= 1, got $parsed."))
    return parsed
end

@inline function _machine_parallel_class_vleo()::Symbol
    cpu_threads = Sys.CPU_THREADS
    if cpu_threads >= 24
        return :large
    elseif cpu_threads >= 12
        return :medium
    end
    return :small
end

@inline function _study_planet()
    cached = _STUDY_PLANET_CACHE[]
    if cached === nothing
        cached = SM.Earth("", SPICE_PATH)
        _STUDY_PLANET_CACHE[] = cached
    end
    return cached
end

@inline function _outer_route_state_path(out_dir_abs::String)::String
    return joinpath(out_dir_abs, OUTER_ROUTE_STATE_FILENAME)
end

@inline function _vleo_process_workers_target(total_cases::Int, override::Union{Nothing, Int})::Int
    target = override === nothing ? max(1, Sys.CPU_THREADS - 1) : override
    return max(1, min(total_cases, target))
end

function _trapz(times::AbstractVector{<:Real}, values::AbstractVector{<:Real})::Float64
    n = min(length(times), length(values))
    n <= 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        dt = Float64(times[i + 1]) - Float64(times[i])
        dt <= 0.0 && continue
        acc += 0.5 * dt * (Float64(values[i]) + Float64(values[i + 1]))
    end
    return acc
end

@inline function _speed_norm(vx::AbstractVector{<:Real}, vy::AbstractVector{<:Real}, vz::AbstractVector{<:Real})::Vector{Float64}
    n = min(length(vx), length(vy), length(vz))
    speeds = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        speeds[i] = sqrt(Float64(vx[i])^2 + Float64(vy[i])^2 + Float64(vz[i])^2)
    end
    return speeds
end

@inline function _initial_period_s(planet, perigee_km::Float64, apogee_km::Float64)::Float64
    rp = planet.Rp_e + perigee_km * 1e3
    ra = planet.Rp_e + apogee_km * 1e3
    a = 0.5 * (rp + ra)
    return 2pi * sqrt(a^3 / planet.μ)
end

function _case_outer_features(perigee_km::Float64, apogee_km::Float64)::SpaceAGORA.OuterRouteFeatures
    planet = _study_planet()
    period_s = _initial_period_s(planet, perigee_km, apogee_km)
    return SpaceAGORA.OuterRouteFeatures(
        category="deterministic",
        n_sats=1,
        n_links=3,
        max_links_per_sat=3,
        mission_time_s=1.1 * period_s,
        has_nbody=false,
        has_srp=false,
        harmonics_degree=0,
        has_control=false,
        orientation_on=false,
        density_family="gram_point",
        solver_mode="auto",
        dt_max_orbit_s=20.0,
        control_rate_s=0.0,
        guidance_rate_s=0.0,
        navigation_rate_s=0.0,
        gram_surrogate_enabled=false,
        gram_static_grid_enabled=false,
        control_effector_count=0,
        thermal_enabled=true,
        dynamic_effector_count=2,
        effector_cost_class="heavy",
    )
end

function _select_case_outer_route(
    requested_route::Symbol,
    route_state::SpaceAGORA.OuterRouteState,
    features::SpaceAGORA.OuterRouteFeatures,
)::Symbol
    machine_class = _machine_parallel_class_vleo()
    threads_available = Threads.nthreads() > 1
    if requested_route == :auto
        return SpaceAGORA.select_outer_route!(
            route_state,
            features;
            machine_class=machine_class,
            threads_available=threads_available,
            parallel_enabled=true,
        )
    end

    candidates = SpaceAGORA.outer_route_candidates(
        features;
        machine_class=machine_class,
        threads_available=threads_available,
        parallel_enabled=true,
    )
    if !(requested_route in candidates) && !_OUTER_ROUTE_WARNING_EMITTED[]
        _OUTER_ROUTE_WARNING_EMITTED[] = true
        @warn "Requested outer route '$requested_route' is outside the recommended candidate set $(candidates) for this native GRAM study. Proceeding anyway."
    end
    return requested_route
end

function _base_case_env_pairs(route::Symbol)::Vector{Pair{String, String}}
    outer_active = route == :none ? "0" : "1"
    return Pair{String, String}[
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => outer_active,
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_MULTIBODY_PARALLEL" => "off",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_RHS_BATCH_PARALLEL" => "off",
        "SPACEAGORA_GRAM_STATIC_GRID" => "0",
        "SPACEAGORA_GRAM_OFFLINE_SURROGATE" => "off",
        "SPACEAGORA_SAVE_BUNDLE" => "0",
    ]
end

@inline function _initial_et_s(args)::Float64
    return SM.ephemerides_time_seconds(args.initial_time, args.environment_model.ephemerides_model)
end

function _recompute_drag_force_ii(
    args,
    t::Float64,
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    altitude_m::Float64,
    latitude_deg::Float64,
    longitude_deg::Float64,
)::SVector{3, Float64}
    density_model = args.environment_model.density_model
    planet = args.environment_model.planet
    ephemerides_model = args.environment_model.ephemerides_model
    spacecraft = args.dynamics_model.spacecraft[1]

    lat_rad = deg2rad(latitude_deg)
    lon_rad = deg2rad(longitude_deg)
    rho, T, wind = GRAMSuite.point_density_state(density_model.core, altitude_m, lat_rad, lon_rad, t, false)
    if !isfinite(rho) || rho <= eps(Float64) || !isfinite(T) || T <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    et = _initial_et_s(args) + t
    l_pi = SM.planet_frame_lpi(planet, et, ephemerides_model)
    pos_pp, vel_pp = SpaceAGORA.SimulationEngine.r_intor_p!(pos_ii, vel_ii, planet, et, ephemerides_model)
    h_pp = cross(pos_pp, vel_pp)
    h_pp_mag = norm(h_pp)
    if !isfinite(h_pp_mag) || h_pp_mag <= eps(Float64)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    sound_velocity = sqrt(planet.γ * planet.R * T)
    uD, uN, uE = SpaceAGORA.SimulationEngine.latlongtoNED((altitude_m, lat_rad, lon_rad))
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    vel_pp_rw = vel_pp + wind_pp
    vel_pp_rw_mag = norm(vel_pp_rw)
    if !isfinite(vel_pp_rw_mag) || vel_pp_rw_mag <= eps(Float64)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    mach = vel_pp_rw_mag / sound_velocity
    S = sqrt(planet.γ * 0.5) * mach
    drag_pp_hat = -(vel_pp_rw / vel_pp_rw_mag)
    q = 0.5 * rho * vel_pp_rw_mag^2
    l_pi_t = l_pi'
    θ_body = acos(clamp(vel_pp_rw[1] / vel_pp_rw_mag, -1.0, 1.0))

    drag_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for body in spacecraft.links
        α_body = if body.root
            pi / 2
        else
            body_frame_velocity = SM.rot(body.q) * SVector{3, Float64}(1.0, 0.0, 0.0)
            atan(body_frame_velocity[1], body_frame_velocity[3])
        end
        CL_body, CD_body, CS_body, _, _, _ = SM.DynamicEffectors.AerodynamicEffectors.aerodynamic_coefficient_fM(body, T, S, α_body, 0.0, θ_body)
        drag_pp_body = q * CD_body * body.ref_area * drag_pp_hat
        drag_ii .+= l_pi_t * drag_pp_body
    end
    return SVector{3, Float64}(drag_ii)
end

function _make_reference_spacecraft(planet, perigee_km::Float64, apogee_km::Float64)
    ic = SM.InitialCondition(
        ra=planet.Rp_e + apogee_km * 1e3,
        rp=planet.Rp_e + perigee_km * 1e3,
        i=REFERENCE_INCLINATION_DEG,
        ω=0.0,
        Ω=0.0,
        ν=180.0,
    )
    return TV.make_three_body_spacecraft(
        bus_dims=REFERENCE_VEHICLE.bus_dims_m,
        panel_dims=REFERENCE_VEHICLE.panel_dims_m,
        bus_mass=REFERENCE_VEHICLE.bus_mass_kg,
        panel_mass_each=REFERENCE_VEHICLE.panel_mass_each_kg,
        panel_offset_y=REFERENCE_VEHICLE.panel_offset_y_m,
        ic=ic,
        prop_mass=0.0,
        id=1,
    )
end

function _build_case_args(; perigee_km::Float64, apogee_km::Float64, channel::Symbol, results_directory::String)
    planet = SM.Earth("", SPICE_PATH)
    spacecraft = _make_reference_spacecraft(planet, perigee_km, apogee_km)
    period_s = _initial_period_s(planet, perigee_km, apogee_km)
    density_model = SM.GRAMAtmosphereModel(
        planet_name="earth",
        initial_time=REFERENCE_INITIAL_TIME,
        gram_density_channel=channel,
        earth_daily_f10=REFERENCE_SOLAR_CASE.daily_f10,
        earth_mean_f10=REFERENCE_SOLAR_CASE.mean_f10,
        earth_ap=REFERENCE_SOLAR_CASE.ap,
    )
    return TV.make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=1.1 * period_s,
        initial_time=REFERENCE_INITIAL_TIME,
        dynamic_effectors=(SM.InverseSquaredGravityModel(), SM.AerodynamicCoefficientfM()),
        density_model=density_model,
        ephemerides_model=SM.SimpleEphemeridesModel(),
        orientation_sim=false,
        keplerian=false,
        EI_km=900.0,
        verbose=false,
        results=true,
        results_directory=results_directory,
    )
end

@inline _case_id(perigee_km::Float64, apogee_km::Float64, channel::Symbol) =
    @sprintf("earth_vleo_p%03.0f_a%03.0f_%s", perigee_km, apogee_km, String(channel))

function _case_metrics(df::DataFrame, args, perigee_km::Float64, apogee_km::Float64, channel::Symbol)
    nrow(df) >= 3 || throw(ArgumentError("Expected at least 3 saved samples, got $(nrow(df)) for $(_case_id(perigee_km, apogee_km, channel))."))

    times = Float64.(df.time)
    altitude_m = Float64.(df.sc1_altitude)
    latitude_deg = Float64.(df.sc1_latitude_deg)
    longitude_deg = Float64.(df.sc1_longitude_deg)
    pos_x = Float64.(df.sc1_pos_1)
    pos_y = Float64.(df.sc1_pos_2)
    pos_z = Float64.(df.sc1_pos_3)
    vel_x = Float64.(df.sc1_vel_1)
    vel_y = Float64.(df.sc1_vel_2)
    vel_z = Float64.(df.sc1_vel_3)
    mass_kg = Float64.(df.sc1_mass)

    peri_idx = argmin(altitude_m)
    peri_idx < nrow(df) || throw(ArgumentError("Could not find a post-periapsis apoapsis sample for $(_case_id(perigee_km, apogee_km, channel))."))
    apo_rel_idx = argmax(@view altitude_m[peri_idx:end])
    apo_idx = peri_idx + apo_rel_idx - 1
    apo_idx > peri_idx || throw(ArgumentError("Could not find a distinct next apoapsis sample for $(_case_id(perigee_km, apogee_km, channel))."))

    drag_vectors = Vector{SVector{3, Float64}}(undef, apo_idx)
    @inbounds for i in 1:apo_idx
        pos_ii = SVector{3, Float64}(pos_x[i], pos_y[i], pos_z[i])
        vel_ii = SVector{3, Float64}(vel_x[i], vel_y[i], vel_z[i])
        drag_vectors[i] = _recompute_drag_force_ii(
            args,
            times[i],
            pos_ii,
            vel_ii,
            altitude_m[i],
            latitude_deg[i],
            longitude_deg[i],
        )
    end
    drag_x = [drag_vectors[i][1] for i in 1:apo_idx]
    drag_y = [drag_vectors[i][2] for i in 1:apo_idx]
    drag_z = [drag_vectors[i][3] for i in 1:apo_idx]
    speeds = _speed_norm(vel_x, vel_y, vel_z)
    drag_norm = _speed_norm(drag_x, drag_y, drag_z)
    mass0 = mass_kg[1]
    along_drag_accel = Vector{Float64}(undef, apo_idx)
    @inbounds for i in 1:apo_idx
        speed = speeds[i]
        if speed <= eps(Float64)
            along_drag_accel[i] = 0.0
        else
            along_drag_accel[i] = -(
                drag_x[i] * vel_x[i] +
                drag_y[i] * vel_y[i] +
                drag_z[i] * vel_z[i]
            ) / speed / mass0
        end
    end

    drag_impulse_ns = _trapz(@view(times[1:apo_idx]), @view(drag_norm[1:apo_idx]))
    required_reboost_dv_mps = _trapz(@view(times[1:apo_idx]), along_drag_accel)
    period_s = _initial_period_s(args.environment_model.planet, perigee_km, apogee_km)
    apogee_speed_next_mps = speeds[apo_idx]
    apogee_speed_delta_mps = apogee_speed_next_mps - speeds[1]

    return (
        case_id=_case_id(perigee_km, apogee_km, channel),
        perigee_km=perigee_km,
        apogee_km=apogee_km,
        channel=String(channel),
        period_s=period_s,
        periapsis_time_s=times[peri_idx],
        next_apoapsis_time_s=times[apo_idx],
        observed_periapsis_km=altitude_m[peri_idx] * 1e-3,
        initial_apogee_km=altitude_m[1] * 1e-3,
        next_apogee_km=altitude_m[apo_idx] * 1e-3,
        drag_impulse_ns=drag_impulse_ns,
        required_reboost_dv_mps=required_reboost_dv_mps,
        required_reboost_dv_per_day_mps=required_reboost_dv_mps * 86400.0 / period_s,
        apogee_loss_km=(altitude_m[apo_idx] - altitude_m[1]) * 1e-3,
        apogee_speed_next_mps=apogee_speed_next_mps,
        apogee_speed_delta_mps=apogee_speed_delta_mps,
    )
end

function _run_case(;
    perigee_km::Float64,
    apogee_km::Float64,
    channel::Symbol,
    route::Symbol=:none,
    apply_env::Bool=true,
)
    env_pairs = _base_case_env_pairs(route)
    return mktempdir() do tmp
        results_directory = joinpath(tmp, "results")
        args = _build_case_args(
            perigee_km=perigee_km,
            apogee_km=apogee_km,
            channel=channel,
            results_directory=results_directory,
        )
        run_case = () -> SpaceAGORA.run_simulation(args)
        if apply_env
            withenv(env_pairs...) do
                run_case()
            end
        else
            run_case()
        end
        csv_path = joinpath(results_directory, "simulation_results.csv")
        isfile(csv_path) || throw(ArgumentError("Expected study run to write $csv_path"))
        df = CSV.read(csv_path, DataFrame)
        return _case_metrics(df, args, perigee_km, apogee_km, channel)
    end
end

function _prime_vleo_worker_bindings!()::Nothing
    Base.invokelatest(
        SM.GRAMAtmosphereModel;
        planet_name="earth",
        initial_time=REFERENCE_INITIAL_TIME,
        gram_density_channel=:nominal,
        earth_daily_f10=REFERENCE_SOLAR_CASE.daily_f10,
        earth_mean_f10=REFERENCE_SOLAR_CASE.mean_f10,
        earth_ap=REFERENCE_SOLAR_CASE.ap,
    )
    return nothing
end

function _ensure_vleo_workers!(target_workers::Int)::Nothing
    missing_workers = target_workers - nworkers()
    if missing_workers > 0
        worker_project = joinpath(REPO_ROOT, ".AGORA")
        addprocs(
            missing_workers;
            exeflags=`--startup-file=no --project=$(worker_project)`
        )
    end

    script_path = abspath(@__FILE__)
    for w in workers()
        remotecall_wait(w) do
            if !isdefined(Main, :run_vleo_drag_trade)
                include(script_path)
            end
            return nothing
        end
        remotecall_wait(_prime_vleo_worker_bindings!, w)
    end
    return nothing
end

function _execute_case_task(task; apply_env::Bool=true)
    started_ns = time_ns()
    try
        row = _run_case(
            perigee_km=Float64(task.perigee_km),
            apogee_km=Float64(task.apogee_km),
            channel=task.channel,
            route=task.route,
            apply_env=apply_env,
        )
        elapsed_s = (time_ns() - started_ns) * 1e-9
        return (
            success=true,
            row=merge(row, (outer_route=String(task.route), solve_wall_s=elapsed_s)),
            elapsed_s=elapsed_s,
            error_message="",
        )
    catch err
        elapsed_s = (time_ns() - started_ns) * 1e-9
        err_text = sprint(showerror, err, catch_backtrace())
        return (
            success=false,
            row=nothing,
            elapsed_s=elapsed_s,
            error_message=err_text,
        )
    end
end

@inline function _record_case_feedback!(
    route_state::SpaceAGORA.OuterRouteState,
    task,
    result,
)::Nothing
    SpaceAGORA.record_outer_route_feedback!(
        route_state,
        task.features;
        route=task.route,
        successes=result.success ? 1 : 0,
        failures=result.success ? 0 : 1,
        elapsed_success_s=result.success ? Float64(result.elapsed_s) : 0.0,
    )
    return nothing
end

function _run_process_tasks!(results, tasks::Vector, process_workers::Union{Nothing, Int}, route_state)::Nothing
    isempty(tasks) && return nothing
    target_workers = _vleo_process_workers_target(length(tasks), process_workers)
    _ensure_vleo_workers!(target_workers)
    payload = pmap(task -> Main._execute_case_task(task; apply_env=true), tasks)
    @inbounds for i in eachindex(tasks)
        task = tasks[i]
        result = payload[i]
        results[task.index] = result
        _record_case_feedback!(route_state, task, result)
    end
    return nothing
end

function _run_serial_tasks!(results, tasks::Vector, route_state)::Nothing
    @inbounds for task in tasks
        result = _execute_case_task(task; apply_env=true)
        results[task.index] = result
        _record_case_feedback!(route_state, task, result)
    end
    return nothing
end

function _run_thread_tasks!(results, tasks::Vector, route_state)::Nothing
    isempty(tasks) && return nothing
    env_pairs = _base_case_env_pairs(:threads)
    local_results = Vector{Any}(undef, length(tasks))
    withenv(env_pairs...) do
        Threads.@threads for idx in eachindex(tasks)
            local_results[idx] = _execute_case_task(tasks[idx]; apply_env=false)
        end
    end
    @inbounds for i in eachindex(tasks)
        task = tasks[i]
        result = local_results[i]
        results[task.index] = result
        _record_case_feedback!(route_state, task, result)
    end
    return nothing
end

@inline function _value_tol(a::Float64, b::Float64)::Float64
    return max(1e-9, 1e-6 * max(abs(a), abs(b), 1.0))
end

function _assert_ordered(low::Float64, nominal::Float64, high::Float64, label::AbstractString)
    tol_low = _value_tol(low, nominal)
    low <= nominal + tol_low || throw(ArgumentError("$label violated low <= nominal: $low vs $nominal"))
    tol_high = _value_tol(nominal, high)
    nominal <= high + tol_high || throw(ArgumentError("$label violated nominal <= high: $nominal vs $high"))
    return nothing
end

function validate_vleo_drag_trade!(cases_df::DataFrame)
    for apogee_km in sort(unique(Float64.(cases_df.apogee_km)))
        for perigee_km in sort(unique(Float64.(cases_df.perigee_km)))
            mask = (cases_df.apogee_km .== apogee_km) .& (cases_df.perigee_km .== perigee_km)
            local_df = sort(cases_df[mask, :], :channel)
            if nrow(local_df) == length(VLEO_CHANNELS)
                by_channel = Dict(String(row.channel) => row for row in eachrow(local_df))
                haskey(by_channel, "low") || continue
                haskey(by_channel, "nominal") || continue
                haskey(by_channel, "high") || continue
                _assert_ordered(
                    Float64(by_channel["low"].drag_impulse_ns),
                    Float64(by_channel["nominal"].drag_impulse_ns),
                    Float64(by_channel["high"].drag_impulse_ns),
                    "drag_impulse_ns at apogee=$(apogee_km) km, perigee=$(perigee_km) km",
                )
                _assert_ordered(
                    Float64(by_channel["low"].required_reboost_dv_mps),
                    Float64(by_channel["nominal"].required_reboost_dv_mps),
                    Float64(by_channel["high"].required_reboost_dv_mps),
                    "required_reboost_dv_mps at apogee=$(apogee_km) km, perigee=$(perigee_km) km",
                )
            end
        end

        for channel in string.(VLEO_CHANNELS)
            mask = (cases_df.apogee_km .== apogee_km) .& (cases_df.channel .== channel)
            local_df = sort(cases_df[mask, :], :perigee_km, rev=true)
            if nrow(local_df) > 1
                for i in 1:(nrow(local_df) - 1)
                    dv_a = Float64(local_df.required_reboost_dv_mps[i])
                    dv_b = Float64(local_df.required_reboost_dv_mps[i + 1])
                    dv_a <= dv_b + _value_tol(dv_a, dv_b) || throw(ArgumentError(
                        "required_reboost_dv_mps is not nondecreasing as perigee decreases for apogee=$(apogee_km) km channel=$channel"
                    ))

                    loss_a = abs(Float64(local_df.apogee_loss_km[i]))
                    loss_b = abs(Float64(local_df.apogee_loss_km[i + 1]))
                    loss_a <= loss_b + _value_tol(loss_a, loss_b) || throw(ArgumentError(
                        "|apogee_loss_km| is not nondecreasing as perigee decreases for apogee=$(apogee_km) km channel=$channel"
                    ))
                end
            end
        end
    end
    return nothing
end

function _write_markdown_table(io, df::DataFrame, cols::Vector{Symbol})
    headers = String.(cols)
    println(io, "| ", join(headers, " | "), " |")
    println(io, "| ", join(fill("---", length(cols)), " | "), " |")
    for row in eachrow(df)
        cells = String[]
        for col in cols
            value = row[col]
            if value isa AbstractFloat
                push!(cells, @sprintf("%.4f", Float64(value)))
            else
                push!(cells, string(value))
            end
        end
        println(io, "| ", join(cells, " | "), " |")
    end
end

function _write_report(report_path::String, cases_df::DataFrame, summary_df::DataFrame, plot_paths::Vector{String})
    top_df = first(sort(cases_df, :required_reboost_dv_per_day_mps, rev=true), min(5, nrow(cases_df)))
    route_counts = combine(groupby(cases_df, :outer_route), nrow => :count)
    open(report_path, "w") do io
        println(io, "# VLEO Drag Trade Study")
        println(io)
        println(io, "## Summary")
        println(io)
        println(io, "- Scope: one Starlink-like reference vehicle, Earth GRAM deterministic `low/nominal/high`, perigee `150-300 km`, apogee `350/550/800 km`.")
        println(io, "- Atmosphere epoch: `2026-03-20 12:00:00 UTC`.")
        println(io, "- Solar inputs: `DailyF10=$(REFERENCE_SOLAR_CASE.daily_f10)`, `MeanF10=$(REFERENCE_SOLAR_CASE.mean_f10)`, `Ap=$(REFERENCE_SOLAR_CASE.ap)`.")
        println(io, "- Aerodynamics: `AerodynamicCoefficientfM()` geometry-based free-molecular model, not a fixed scalar `Cd`.")
        println(io, "- Validation: channel ordering and perigee-trend checks passed for the completed matrix.")
        println(io, "- Outer routes used: `$(join(["$(row.outer_route)=$(row.count)" for row in eachrow(route_counts)], ", "))`.")
        println(io)
        println(io, "## Reference Vehicle")
        println(io)
        println(io, "- Vehicle: $(REFERENCE_VEHICLE.name)")
        println(io, "- Total mass: `$(REFERENCE_VEHICLE.mass_total_kg) kg`")
        println(io, "- Bus mass: `$(REFERENCE_VEHICLE.bus_mass_kg) kg`")
        println(io, "- Panel mass each: `$(REFERENCE_VEHICLE.panel_mass_each_kg) kg`")
        println(io, "- Bus dimensions: `$(REFERENCE_VEHICLE.bus_dims_m) m`")
        println(io, "- Panel dimensions: `$(REFERENCE_VEHICLE.panel_dims_m) m`")
        println(io, "- Panel offset Y: `$(REFERENCE_VEHICLE.panel_offset_y_m) m`")
        println(io, "- Inclination: `$(REFERENCE_INCLINATION_DEG) deg`")
        println(io, "- Fixed attitude: `orientation_sim=false`")
        println(io)
        println(io, "## Highest-Drag Cases")
        println(io)
        _write_markdown_table(io, top_df[:, [:case_id, :perigee_km, :apogee_km, :channel, :required_reboost_dv_per_day_mps, :drag_impulse_ns]], [:case_id, :perigee_km, :apogee_km, :channel, :required_reboost_dv_per_day_mps, :drag_impulse_ns])
        println(io)
        println(io, "## Artifacts")
        println(io)
        println(io, "- Cases CSV: `$(joinpath(dirname(report_path), "cases.csv"))`")
        println(io, "- Summary CSV: `$(joinpath(dirname(report_path), "summary.csv"))`")
        println(io, "- Outer-route state: `$(joinpath(dirname(report_path), OUTER_ROUTE_STATE_FILENAME))`")
        println(io, "- Completed cases: `$(nrow(summary_df))`")
        if !isempty(plot_paths)
            println(io, "- Plot artifacts:")
            for path in plot_paths
                println(io, "  - `$(path)`")
            end
        end
    end
    return report_path
end

function _metric_plot_path(out_dir::String, basename::String)::String
    plots_dir = joinpath(out_dir, "plots")
    mkpath(plots_dir)
    return joinpath(plots_dir, basename * ".png")
end

function _plot_metric(cases_df::DataFrame, metric::Symbol, ylabel::String, out_path::String)
    apogees = sort(unique(Float64.(cases_df.apogee_km)))
    plot_count = max(1, length(apogees))
    plt = Plots.plot(layout=(plot_count, 1), size=(1200, 320 * plot_count), legend=:topright)
    for (idx, apogee_km) in enumerate(apogees)
        local_df = cases_df[cases_df.apogee_km .== apogee_km, :]
        for channel in string.(VLEO_CHANNELS)
            channel_df = sort(local_df[local_df.channel .== channel, :], :perigee_km)
            nrow(channel_df) == 0 && continue
            Plots.plot!(
                plt,
                channel_df.perigee_km,
                channel_df[!, metric];
                subplot=idx,
                label=channel,
                color=get(CHANNEL_COLORS, channel, "#374151"),
                marker=:circle,
                linewidth=2,
                xlabel=idx == plot_count ? "Perigee [km]" : "",
                ylabel=ylabel,
                title=@sprintf("Apogee %.0f km", apogee_km),
            )
        end
    end
    Plots.savefig(plt, out_path)
    return out_path
end

function _load_outer_route_state_vleo(path::String)::SpaceAGORA.OuterRouteState
    state = SpaceAGORA.OuterRouteState()
    if isfile(path)
        SpaceAGORA.ParallelProfiles.load_outer_route_state!(state, path; replace=true)
    end
    return state
end

function _save_outer_route_state_vleo(state::SpaceAGORA.OuterRouteState, path::String)
    return SpaceAGORA.ParallelProfiles.save_outer_route_state(
        state,
        path;
        metadata=Dict(
            "study" => "vleo_drag_trade",
            "updated_by" => "benchmarks/studies/vleo_drag_trade.jl",
            "updated_utc" => string(now(UTC)),
        ),
    )
end

function _build_case_tasks(
    perigee_grid::Vector{Float64},
    apogee_grid::Vector{Float64},
    requested_outer_route::Symbol,
    route_state::SpaceAGORA.OuterRouteState,
)::Vector{NamedTuple}
    tasks = NamedTuple[]
    task_index = 0
    for apogee_km in apogee_grid
        for perigee_km in perigee_grid
            features = _case_outer_features(perigee_km, apogee_km)
            route = _select_case_outer_route(requested_outer_route, route_state, features)
            for channel in VLEO_CHANNELS
                task_index += 1
                push!(tasks, (
                    index=task_index,
                    case_id=_case_id(perigee_km, apogee_km, channel),
                    perigee_km=perigee_km,
                    apogee_km=apogee_km,
                    channel=channel,
                    route=route,
                    features=features,
                ))
            end
        end
    end
    return tasks
end

function _route_counts_string(tasks::Vector)::String
    counts = Dict(route => 0 for route in (:none, :threads, :process))
    for task in tasks
        counts[task.route] = get(counts, task.route, 0) + 1
    end
    active = String[]
    for route in (:none, :threads, :process)
        count = get(counts, route, 0)
        count > 0 && push!(active, "$(route)=$(count)")
    end
    return isempty(active) ? "none=0" : join(active, ", ")
end

function run_vleo_drag_trade(;
    smoke::Bool=false,
    out_dir::String=DEFAULT_OUT_DIR,
    generate_plots::Bool=true,
    perigees_km::Union{Nothing, AbstractVector{<:Real}}=nothing,
    apogees_km::Union{Nothing, AbstractVector{<:Real}}=nothing,
    outer_route::Symbol=:auto,
    process_workers::Union{Nothing, Int}=nothing,
)
    out_dir_abs = isabspath(out_dir) ? normpath(out_dir) : normpath(joinpath(pwd(), out_dir))
    mkpath(out_dir_abs)

    perigee_grid = perigees_km === nothing ? (smoke ? [150.0] : FULL_PERIGEES_KM) : Float64.(collect(perigees_km))
    apogee_grid = apogees_km === nothing ? (smoke ? [350.0] : FULL_APOGEES_KM) : Float64.(collect(apogees_km))
    route_state_path = _outer_route_state_path(out_dir_abs)
    route_state = _load_outer_route_state_vleo(route_state_path)
    tasks = _build_case_tasks(perigee_grid, apogee_grid, outer_route, route_state)
    results = Vector{Any}(undef, length(tasks))

    println("VLEO drag trade study")
    println("perigees_km = $(perigee_grid)")
    println("apogees_km = $(apogee_grid)")
    println("channels = $(collect(VLEO_CHANNELS))")
    println("requested_outer_route = $(outer_route)")
    println("selected_outer_routes = $(_route_counts_string(tasks))")

    serial_tasks = [task for task in tasks if task.route == :none]
    thread_tasks = [task for task in tasks if task.route == :threads]
    process_tasks = [task for task in tasks if task.route == :process]

    for task in serial_tasks
        println("Running $(task.case_id) [route=$(task.route)]")
    end
    for task in thread_tasks
        println("Queued $(task.case_id) [route=$(task.route)]")
    end
    for task in process_tasks
        println("Queued $(task.case_id) [route=$(task.route)]")
    end

    try
        _run_serial_tasks!(results, serial_tasks, route_state)
        _run_thread_tasks!(results, thread_tasks, route_state)
        _run_process_tasks!(results, process_tasks, process_workers, route_state)
    finally
        _save_outer_route_state_vleo(route_state, route_state_path)
    end

    rows = NamedTuple[]
    for task in tasks
        result = results[task.index]
        result === nothing && throw(ArgumentError("Missing result for $(task.case_id)."))
        result.success || throw(ErrorException("Study case $(task.case_id) failed: $(result.error_message)"))
        push!(rows, result.row)
    end

    cases_df = sort(DataFrame(rows), [:apogee_km, :perigee_km, :channel])
    validate_vleo_drag_trade!(cases_df)

    summary_df = select(
        cases_df,
        :case_id,
        :perigee_km,
        :apogee_km,
        :channel,
        :outer_route,
        :drag_impulse_ns,
        :required_reboost_dv_mps,
        :required_reboost_dv_per_day_mps,
        :apogee_loss_km,
        :apogee_speed_next_mps,
        :apogee_speed_delta_mps,
        :solve_wall_s,
    )

    cases_csv_path = joinpath(out_dir_abs, "cases.csv")
    summary_csv_path = joinpath(out_dir_abs, "summary.csv")
    CSV.write(cases_csv_path, cases_df)
    CSV.write(summary_csv_path, summary_df)

    plot_paths = String[]
    if generate_plots
        push!(plot_paths, _plot_metric(cases_df, :required_reboost_dv_mps, "Required reboost dv per orbit [m/s]", _metric_plot_path(out_dir_abs, "required_reboost_dv_per_orbit")))
        push!(plot_paths, _plot_metric(cases_df, :required_reboost_dv_per_day_mps, "Required reboost dv per day [m/s]", _metric_plot_path(out_dir_abs, "required_reboost_dv_per_day")))
        push!(plot_paths, _plot_metric(cases_df, :apogee_loss_km, "Apogee change [km]", _metric_plot_path(out_dir_abs, "apogee_loss")))
        push!(plot_paths, _plot_metric(cases_df, :apogee_speed_delta_mps, "Apogee speed delta [m/s]", _metric_plot_path(out_dir_abs, "apogee_speed_delta")))
    end

    report_path = joinpath(out_dir_abs, "report.md")
    _write_report(report_path, cases_df, summary_df, plot_paths)

    println("Saved cases CSV: $cases_csv_path")
    println("Saved summary CSV: $summary_csv_path")
    println("Saved report: $report_path")
    println("Saved outer-route state: $route_state_path")
    println("Saved plots: $(length(plot_paths))")

    return (
        cases_df=cases_df,
        summary_df=summary_df,
        cases_csv_path=cases_csv_path,
        summary_csv_path=summary_csv_path,
        report_path=report_path,
        plot_paths=plot_paths,
        outer_route_state_path=route_state_path,
    )
end

function run_vleo_drag_trade_cli(args::Vector{String}=copy(ARGS))
    opts = parse_cli_vleo(args)
    smoke = haskey(opts, "smoke") ? parse_bool_vleo(opts["smoke"]) : false
    out_dir = get(opts, "out-dir", DEFAULT_OUT_DIR)
    generate_plots = haskey(opts, "plots") ? parse_bool_vleo(opts["plots"]) : true
    outer_route = haskey(opts, "outer-route") ? _parse_outer_route_vleo(opts["outer-route"]) : :auto
    process_workers = haskey(opts, "process-workers") ? _parse_positive_int_vleo("--process-workers", opts["process-workers"]) : nothing
    return run_vleo_drag_trade(
        smoke=smoke,
        out_dir=out_dir,
        generate_plots=generate_plots,
        outer_route=outer_route,
        process_workers=process_workers,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_vleo_drag_trade_cli()
end
