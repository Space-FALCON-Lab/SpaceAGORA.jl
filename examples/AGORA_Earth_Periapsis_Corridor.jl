include(joinpath(@__DIR__, "common.jl"))
using CSV
using DataFrames
using LinearAlgebra
using StaticArrays
using TOML
using Dates
using SHA
using DiffEqBase: ContinuousCallback, DiscreteCallback, terminate!

include(joinpath(@__DIR__, "periapsis_corridor_support.jl"))
include(joinpath(@__DIR__, "periapsis_corridor_mean_elements.jl"))
include(joinpath(@__DIR__, "periapsis_corridor_controller.jl"))
include(joinpath(@__DIR__, "periapsis_corridor_plots.jl"))
using .PeriapsisCorridorSupport
using .PeriapsisCorridorMeanElements
using .PeriapsisCorridorController

struct CorridorBallisticDrag <: AbstractForceTorqueModel
    ballistic_coefficient_kg_m2::Float64
    reference_mass_kg::Float64
    atmosphere::CorridorAtmosphereFit
    # :fit evaluates `atmosphere` (US76 piecewise fit); :sampled reads the
    # engine's sampled density (`environment.atmosphere`, e.g. Earth-GRAM).
    density_source::Symbol
end

CorridorBallisticDrag(beta, mass, atmosphere) =
    CorridorBallisticDrag(beta, mass, atmosphere, :fit)

mutable struct CorridorGuidancePolicy <: SimulationModel.AbstractGuidanceModel
    config::CorridorConfig
    phase::Symbol
    policy::Symbol
    max_dv_mps::Float64
    max_burn_duration_s::Float64
    remaining_budget_mps::Float64
    deadband_m::Float64
    target_altitude_m::Float64
    lead_angle_rad::Float64
    periapsis_count::Int
    armed::Bool
    last_selection::Any
    decisions::Vector{NamedTuple}
    burn_plans::Vector{NamedTuple}
    failures::Vector{String}
end

mutable struct CorridorFixedBurnGuidance <: SimulationModel.AbstractGuidanceModel
    delta_v_mps::Float64
    source_orbit::Int64
    armed::Bool
end

function SimulationModel.GuidanceHooks.calcGuidanceEffect!(guidance::CorridorFixedBurnGuidance,
        u, p, t::Float64, sat_idx::Int64)
    commands = p.shared_buffers.maneuver_commands
    (sat_idx < 1 || sat_idx > length(commands)) && return nothing
    commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(
        valid=true, delta_v_mps=guidance.armed ? abs(guidance.delta_v_mps) : 0.0,
        direction_rad=guidance.delta_v_mps >= 0.0 ? 0.0 : π,
        source_orbit=guidance.source_orbit)
    return nothing
end

function SimulationModel.GuidanceHooks.calcGuidanceEffect!(policy::CorridorGuidancePolicy,
        u, p, t::Float64, sat_idx::Int64)
    commands = p.shared_buffers.maneuver_commands
    (sat_idx < 1 || sat_idx > length(commands)) && return nothing
    if !policy.armed
        commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(
            valid=true, source_orbit=policy.periapsis_count + 1)
        return nothing
    end

    r = SVector{3, Float64}(u.sc[sat_idx].pos)
    v = SVector{3, Float64}(u.sc[sat_idx].vel)
    mass = Float64(u.sc[sat_idx].mass)
    vehicle = policy.config.values["vehicle"]
    exhaust_velocity = Float64(vehicle["specific_impulse_s"]) * 9.80665
    consumed_at_limit = Float64(vehicle["thrust_n"]) * policy.max_burn_duration_s / exhaust_velocity
    duration_limited_dv = consumed_at_limit < mass ?
        exhaust_velocity * log(mass / (mass - consumed_at_limit)) : Inf
    # Budget is reduced by delivered delta-v (rocket equation on the mass
    # history; earlier burns have completed by the next pre-apoapsis request).
    wet_mass = Float64(vehicle["dry_mass_kg"]) + Float64(vehicle["propellant_mass_kg"])
    delivered_so_far = mass < wet_mass ? exhaust_velocity * log(wet_mass / mass) : 0.0
    policy.remaining_budget_mps = max(0.0,
        Float64(policy.config.values["comparison"]["total_delta_v_budget_mps"]) - delivered_so_far)
    max_allowed_dv = min(policy.max_dv_mps, policy.remaining_budget_mps, duration_limited_dv)
    predictor = dv -> _predict_next_physical_periapsis(policy.config, policy.phase,
        p.args, r, v, Float64(u.sc[sat_idx].mass), Float64(dv); elapsed_s=Float64(t))
    mission = policy.config.values["mission"]
    choice = select_corridor_maneuver(predictor;
        lower_m=Float64(mission["corridor_min_altitude_m"]),
        upper_m=Float64(mission["corridor_max_altitude_m"]),
        margin_m=Float64(mission["safety_margin_m"]),
        target_m=policy.target_altitude_m, policy=policy.policy,
        max_dv_mps=max_allowed_dv,
        deadband_m=policy.deadband_m,
        monotone_tolerance_m=Float64(policy.config.values["comparison"]["predictor_monotone_tolerance_m"]))
    policy.last_selection = choice
    target_pass = policy.periapsis_count + 1
    push!(policy.decisions, (request_time_s=Float64(t), source_orbit=target_pass,
        target_pass=target_pass, policy=String(policy.policy),
        delta_v_mps=choice.delta_v_mps,
        predicted_periapsis_altitude_m=choice.predicted_periapsis_altitude_m,
        target_periapsis_altitude_m=choice.target_periapsis_altitude_m,
        trigger_reason=String(choice.reason), feasible=choice.feasible,
        predictor_evaluations=choice.evaluations))
    commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(
        valid=true, delta_v_mps=choice.feasible ? abs(choice.delta_v_mps) : 0.0,
        direction_rad=choice.delta_v_mps >= 0.0 ? 0.0 : π,
        source_orbit=target_pass)
    policy.armed = false
    return nothing
end

function SimulationModel.environment_requirements(::CorridorBallisticDrag)
    return SimulationModel.EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
end

function SimulationModel.wrench(model::CorridorBallisticDrag, state::SimulationModel.StateSample,
                                environment::SimulationModel.EnvironmentSample, t::Float64)
    planet_frame = environment.planet_frame
    if planet_frame === nothing
        zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
        return zero3, zero3
    end
    # Geodetic altitude above the planet ellipsoid: the engine's own sample, the
    # same altitude it hands to GRAM, so the fit and GRAM share one definition.
    altitude_m = planet_frame.alt_m
    density = if model.density_source === :sampled
        environment.atmosphere === nothing && throw(ErrorException(
            "CorridorBallisticDrag(:sampled) needs a sampled environment atmosphere, got none"))
        environment.atmosphere.rho_kg_m3
    else
        corridor_density(model.atmosphere, altitude_m)
    end
    density > eps(Float64) || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    # `_planet_relative_state` removes planet rotation before this hook. The fit
    # is evaluated at geodetic altitude, matching the corridor, and assumes zero
    # winds. The sampled (GRAM) density is evaluated by the engine at the same
    # geodetic altitude and its sampled wind is deliberately NOT applied (same zero-wind assumption).
    # Outside the fit nodes the configured below/above behavior applies.
    velocity_pp = planet_frame.vel_pp
    speed = norm(velocity_pp)
    speed > 0.0 || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    beta_current = model.ballistic_coefficient_kg_m2 * state.mass_kg / model.reference_mass_kg
    force_pp = -0.5 * density * state.mass_kg * speed^2 /
        beta_current * (velocity_pp / speed)
    force_ii = SVector{3, Float64}(transpose(planet_frame.l_pi) * force_pp)
    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

const _CORRIDOR_EPOCH = InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)

# Loader hook for the GRAM backend; tests replace it to simulate a missing GRAM.
const _GRAM_LOADER = Ref{Function}(ensure_gramsuite_loaded!)
# GRAM density evaluations spent inside controller predictor runs (needs
# SPACEAGORA_GRAM_PROFILE, which the GRAM backend sets for its runs).
const _GRAM_PREDICTOR_EVALS = Ref(0)

_gram_density_evaluations() =
    Int(SimulationModel.SimulationCallbacks._gram_runtime_stats_snapshot().density_calls)

"""`t0` advanced by `elapsed_s` (millisecond resolution), as an `InitialTime`."""
function _shift_initial_time(t0, elapsed_s::Float64)
    elapsed_s == 0.0 && return t0
    whole = floor(Int, Float64(t0.second))
    ms0 = round(Int, (Float64(t0.second) - whole) * 1000)
    base = DateTime(Int(t0.year), Int(t0.month), Int(t0.day), Int(t0.hour), Int(t0.minute), whole) +
        Millisecond(ms0 + round(Int, elapsed_s * 1000))
    return InitialTime(year=Dates.year(base), month=Dates.month(base), day=Dates.day(base),
        hour=Dates.hour(base), minute=Dates.minute(base),
        second=Dates.second(base) + Dates.millisecond(base) / 1000)
end

_utc_string(t) = string(Int(t.year), "-", Int(t.month), "-", Int(t.day), " ", Int(t.hour), ":",
    Int(t.minute), ":", Float64(t.second), " UTC")

const _GRAM_OFF_VALUES = ("off", "none", "0", "false", "no")

"""The corridor study needs direct per-state GRAM sampling, not an interpolating cache."""
function _corridor_gram_require_direct_sampling()
    for (name, default) in (("SPACEAGORA_GRAM_TRACK_CACHE", "off"), ("SPACEAGORA_GRAM_SEGMENT_CACHE", "off"),
                            ("SPACEAGORA_VACUUM_GRAM_CACHE", "off"), ("SPACEAGORA_DENSITY_FREEZE_PER_STEP", "off"))
        value = lowercase(strip(get(ENV, name, default)))
        value in _GRAM_OFF_VALUES || throw(ErrorException(
            "atmosphere backend \"gram\" requires direct GRAM sampling, but $name=$value interpolates or freezes density; unset it"))
    end
    return nothing
end

const _GRAM_HINT = "Build native GRAM with `julia --project=. scripts/ensure_gram_native.jl` and see docs/src/user/gramsuite_setup.md."

"""Load GRAMSuite and its SpaceAGORA extension, or fail with a clear message (never a fallback)."""
function _corridor_gram_load!()
    try
        _GRAM_LOADER[]()
    catch err
        throw(ErrorException("atmosphere backend \"gram\" is unavailable: GRAMSuite could not be loaded " *
            "($(sprint(showerror, err))). $_GRAM_HINT The US76 fit was NOT substituted."))
    end
    Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) === nothing && throw(ErrorException(
        "atmosphere backend \"gram\" is unavailable: the SpaceAGORAGRAMSuiteExt extension is not loaded. $_GRAM_HINT The US76 fit was NOT substituted."))
    return nothing
end

"""
Earth-GRAM model at `epoch`, or an `ErrorException` saying what is missing.
There is no fallback to the US76 fit: a GRAM request that cannot be honored fails.
"""
function _corridor_gram_density_model(config::CorridorConfig, epoch)
    _corridor_gram_load!()
    _corridor_gram_require_direct_sampling()
    try
        # Furnishes the SPICE kernels GRAM's ephemeris hook needs (solar time,
        # season); dynamics keep SimpleEphemeridesModel and the local planet.
        Earth("", SPICE_PATH)
        # invokelatest: the extension may have been loaded just now (newer world age).
        return Base.invokelatest(SimulationModel.GRAMAtmosphereModel;
            planet_name=String(config.values["atmosphere_gram"]["planet_name"]), initial_time=epoch)
    catch err
        throw(ErrorException("atmosphere backend \"gram\" is unavailable: GRAM or SPICE data could not be " *
            "initialized ($(sprint(showerror, err))). $_GRAM_HINT The US76 fit was NOT substituted."))
    end
end

"""Manifest record of what GRAM was used (strings/numbers/bools only)."""
function _corridor_gram_provenance(config::CorridorConfig, model, epoch)
    core = model.core
    gs = parentmodule(typeof(core))
    library = try String(gs._find_gram_library(core.gram_root)) catch; "unavailable" end
    sha = isfile(library) ? bytes2hex(open(sha256, library)) : "unavailable"
    project_toml = joinpath(REPO_ROOT, "data", "GRAMSuite.jl", "Project.toml")
    gram_pkg_version = try String(TOML.parsefile(project_toml)["version"]) catch; "unavailable" end
    git_rev = try strip(read(`git -C $(joinpath(REPO_ROOT, "data", "GRAMSuite.jl")) rev-parse HEAD`, String)) catch; "unavailable" end
    return Dict{String, Any}(
        "backend" => "gram",
        "gram_suite_directory" => basename(core.gram_root),
        "gram_root" => String(core.gram_root),
        "gram_data_root" => String(core.gram_data_root),
        "spice_root" => String(core.spice_root),
        "gram_library_path" => library,
        "gram_library_sha256" => sha,
        "gramsuite_jl_version" => gram_pkg_version,
        "gramsuite_jl_git_revision" => git_rev,
        "gram_version_note" => "GRAM exposes no version string through the repo; the directory name and library hash identify it",
        "planet_name" => String(core.planet_name),
        "configured_run_epoch_utc" => _utc_string(_CORRIDOR_EPOCH),
        "gram_core_epoch_utc" => _utc_string(core.initial_time),
        "gram_epoch_matches_run_epoch" => _utc_string(core.initial_time) == _utc_string(epoch),
        "predictor_epoch_rule" => "run epoch + elapsed time at the maneuver request, so predictor GRAM sampling matches the truth atmosphere at that time",
        "entry_interface_km" => Float64(config.values["atmosphere_gram"]["entry_interface_km"]),
        "gram_density_ceiling_m" => 2000.0e3,
        "density_quantity" => "GRAM mean density only; gram_seed and perturbation scales do not enter the trajectory (deterministic run)",
        "altitude_passed_to_gram" => "geodetic height (Bowring) above the planet ellipsoid, repo sample path; the corridor and the US76 fit use the same geodetic altitude",
        "wind_applied_to_drag" => false,
        "sampling" => "direct point-to-point (track cache, vacuum cache and per-step freeze off)",
        "ephemerides_model" => "SimpleEphemeridesModel (dynamics); SPICE only for GRAM's solar-time/season hook",
        "parallelism" => "serial; GRAM calls take the global GRAM lock",
    )
end

"""Direct GRAM mean-density evaluation at an event, through the engine's own `getDensity`."""
_gram_density_at(args, alt_m, lat_rad, lon_rad, t_s) = Float64(
    SimulationModel.EnvironmentModels.getDensity(args.environment_model.density_model,
        Float64(alt_m), Float64(lat_rad), Float64(lon_rad), Float64(t_s), true, (args=args,))[1])

function _assert_case_dir_backend(case_dir::AbstractString, backend::AbstractString)
    path = joinpath(case_dir, "run_manifest.toml")
    isfile(path) || return nothing
    existing = String(get(TOML.parsefile(path), "atmosphere_backend", "us76_piecewise"))
    existing == backend || throw(ErrorException(
        "$case_dir already holds a \"$existing\" run; refusing to mix it with backend \"$backend\". Use a different --output."))
    return nothing
end
function _corridor_cli(args)
    smoke = "--smoke" in args || lowercase(get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "false")) in ("1", "true", "yes", "on")
    case_name = "favorable"
    output_dir = joinpath(REPO_ROOT, "output", "periapsis_corridor")
    horizon_override = nothing
    config_path = joinpath(@__DIR__, "periapsis_corridor_config.toml")
    atmosphere = nothing
    output_given = false
    for arg in args
        startswith(arg, "--case=") && (case_name = split(arg, "=", limit=2)[2])
        startswith(arg, "--output=") && (output_dir = abspath(split(arg, "=", limit=2)[2]); output_given = true)
        startswith(arg, "--atmosphere=") && (atmosphere = split(arg, "=", limit=2)[2])
        startswith(arg, "--horizon-s=") && (horizon_override = parse(Float64, split(arg, "=", limit=2)[2]))
        startswith(arg, "--config=") && (config_path = abspath(split(arg, "=", limit=2)[2]))
    end
    valid_cases = ("vacuum", "unfavorable", "favorable", "near_passive", "frequent", "all")
    case_name in valid_cases || throw(ArgumentError("--case must be one of $(join(valid_cases, ", "))"))
    horizon_override === nothing || horizon_override > 0.0 || throw(ArgumentError("--horizon-s must be positive"))
    atmosphere === nothing || atmosphere in CORRIDOR_ATMOSPHERE_BACKENDS ||
        throw(ArgumentError("--atmosphere must be one of $(join(CORRIDOR_ATMOSPHERE_BACKENDS, ", "))"))
    # A GRAM run defaults to its own tree so it can never land next to US76 results.
    atmosphere == "gram" && !output_given && (output_dir = joinpath(REPO_ROOT, "output", "periapsis_corridor_gram"))
    return (; smoke, case_name, output_dir, horizon_override, config_path, atmosphere)
end

"""
`(t, r_ii) -> (geodetic altitude m, geodetic latitude rad)` for a run's `args`:
the planet-fixed position at `t` goes through the engine's `rtolatlong`.
"""
function _geodetic_of(args)
    env = args.environment_model
    planet = env.planet
    et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model)
    return function (t::Float64, r_ii)
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + t, env.ephemerides_model)
        return geodetic_altitude_latitude(lpi * r_ii, planet)
    end
end

function _apsis_observation_callback(observations::Vector{NamedTuple}, args;
                                    controller::Union{Nothing, CorridorGuidancePolicy}=nothing)
    geodetic_of = _geodetic_of(args)
    condition(u, t, integrator) = dot(u.sc[1].pos, u.sc[1].vel)
    function record!(integrator, kind::String)
        integrator.t > 1e-3 || return nothing
        sc = integrator.u.sc[1]
        push!(observations, (time_s=Float64(integrator.t), kind=kind,
            radius_m=norm(sc.pos), physical_altitude_m=geodetic_of(Float64(integrator.t), sc.pos)[1]))
        kind == "periapsis" && controller !== nothing && (controller.periapsis_count += 1)
        return nothing
    end
    affect_pos!(integrator) = record!(integrator, "periapsis")
    affect_neg!(integrator) = record!(integrator, "apoapsis")
    return ContinuousCallback(condition, affect_pos!, affect_neg!; save_positions=(false, false))
end

function _preapoapsis_guidance_callback(controller::CorridorGuidancePolicy)
    target_nu = π - controller.lead_angle_rad
    function true_anomaly(u, planet)
        elements = _corridor_element_summary(SVector{3, Float64}(u.sc[1].pos),
            SVector{3, Float64}(u.sc[1].vel), planet.μ)
        return elements.nu
    end
    condition(u, t, integrator) = sin(true_anomaly(u, integrator.p.args.environment_model.planet) - target_nu)
    function arm_and_schedule!(integrator)
        dot(integrator.u.sc[1].pos, integrator.u.sc[1].vel) > 0.0 || return nothing
        controller.armed = true
        SimulationModel.SimulationCallbacks.schedule_event_driven_thruster_controls!(integrator, 1)
        controller.armed = false
        if controller.last_selection !== nothing && !controller.last_selection.feasible
            push!(controller.failures, "unsafe predictor result at t=$(integrator.t): $(controller.last_selection.reason) $(controller.last_selection.error)")
            terminate!(integrator)
        elseif controller.last_selection !== nothing && abs(controller.last_selection.delta_v_mps) > 0.0
            selection = controller.last_selection
            plans = integrator.p.shared_buffers.maneuver_burn_plans
            plan = isempty(plans) ? nothing : plans[1]
            valid_plan = plan !== nothing && plan.valid &&
                plan.source_orbit == controller.periapsis_count + 1 &&
                plan.start_burn_s >= integrator.t && plan.stop_burn_s > plan.start_burn_s &&
                plan.stop_burn_s - plan.start_burn_s <= controller.max_burn_duration_s
            if !valid_plan
                push!(controller.failures, "finite-thrust plan rejected or stale at t=$(integrator.t) for source orbit $(controller.periapsis_count + 1)")
                integrator.p.shared_buffers.maneuver_commands[1] =
                    SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true,
                        source_orbit=controller.periapsis_count + 1)
                terminate!(integrator)
                return nothing
            end
            push!(controller.burn_plans, (request_time_s=Float64(integrator.t),
                target_pass=controller.periapsis_count + 1,
                execution_start_time_s=plan.start_burn_s,
                execution_stop_time_s=plan.stop_burn_s,
                target_periapsis_altitude_m=selection.target_periapsis_altitude_m,
                predicted_periapsis_altitude_m=selection.predicted_periapsis_altitude_m,
                requested_delta_v_mps=abs(selection.delta_v_mps),
                direction_rad=selection.delta_v_mps >= 0.0 ? 0.0 : π,
                trigger_reason=String(selection.reason)))
        end
        return nothing
    end
    return ContinuousCallback(condition, arm_and_schedule!, nothing;
        save_positions=(false, false))
end

function _fixed_burn_schedule_callback(controller::CorridorFixedBurnGuidance,
                                       scheduled::Base.RefValue{Bool})
    function schedule!(integrator)
        scheduled[] && return nothing
        controller.armed = true
        SimulationModel.SimulationCallbacks.schedule_event_driven_thruster_controls!(integrator, 1)
        controller.armed = false
        if abs(controller.delta_v_mps) > 0.0
            plans = integrator.p.shared_buffers.maneuver_burn_plans
            plan = isempty(plans) ? nothing : plans[1]
            plan !== nothing && plan.valid && plan.source_orbit == controller.source_orbit ||
                throw(ErrorException("finite-thrust predictor could not schedule the requested burn"))
        end
        scheduled[] = true
        return nothing
    end
    condition(u, t, integrator) = false
    initialize!(callback, u, t, integrator) = schedule!(integrator)
    return DiscreteCallback(condition, schedule!; initialize=initialize!,
        save_positions=(false, false))
end

function _corridor_mean_elements(config::CorridorConfig, phase::Symbol, planet, harmonics)
    j2 = -sqrt(5.0) * harmonics.C[3, 1]
    if phase === :vacuum
        vacuum = config.values["vacuum_orbit"]
        mean = (a=Float64(vacuum["mean_semi_major_axis_m"]),
            e=Float64(vacuum["mean_eccentricity"]),
            i=deg2rad(Float64(vacuum["inclination_deg"])),
            raan=deg2rad(Float64(vacuum["raan_deg"])),
            omega=deg2rad(Float64(vacuum["argument_of_periapsis_deg"])),
            nu=deg2rad(Float64(vacuum["mean_true_anomaly_deg"])))
        mapped = mean_to_osculating(mean; mu=harmonics.gm_m3s2,
            radius=harmonics.reference_radius_m, j2=j2)
        mapped.converged || throw(ErrorException("vacuum mean-to-osculating map did not converge"))
        return mapped.elements
    end
    orbit = config.values["orbit"]
    rp = planet.Rp_e + Float64(orbit["initial_periapsis_altitude_m"])
    ra = planet.Rp_e + Float64(orbit["apoapsis_altitude_m"])
    a, e = 0.5 * (ra + rp), (ra - rp) / (ra + rp)
    omega_deg = phase === :unfavorable ? 0.0 : Float64(orbit["argument_of_periapsis_deg"])
    # This high-e Cartesian seed is deliberately osculating. The J2 map is
    # validated only for the separate low-e vacuum diagnostic.
    return (a=a, e=e, i=deg2rad(Float64(orbit["inclination_deg"])),
        raan=deg2rad(Float64(orbit["raan_deg"])), omega=deg2rad(omega_deg),
        nu=deg2rad(Float64(orbit["true_anomaly_deg"])))
end

function _corridor_args(config::CorridorConfig, phase::Symbol, run_horizon_s::Float64,
                        results_dir::String; initial_state=nothing, initial_mass_kg=nothing,
                        controller::Union{Nothing, SimulationModel.AbstractGuidanceModel}=nothing,
                        epoch=_CORRIDOR_EPOCH)
    gravity_path = joinpath(REPO_ROOT, String(config.values["gravity"]["coefficient_file"]))
    planet_seed = make_no_gram_planet(:earth)
    gravity_seed = GravitationalHarmonicsModel(3, 0, gravity_path, planet_seed)
    # Use the GM carried by the sourced coefficient file for element conversion,
    # initial-condition construction, and the truth gravity model.
    planet = Earth(μ=gravity_seed.gm_m3s2)
    harmonics = GravitationalHarmonicsModel(3, 0, gravity_path, planet)
    osculating = _corridor_mean_elements(config, phase, planet, harmonics)
    initial_r, initial_v = initial_state === nothing ?
        elements_to_cartesian(osculating; mu=harmonics.gm_m3s2) : initial_state
    initial = CartesianInitialCondition(initial_r, initial_v)

    vehicle = config.values["vehicle"]
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.05, 2.05, 2.8), panel_dims=(0.01, 0.2, 0.2),
        bus_mass=Float64(vehicle["dry_mass_kg"]), panel_mass_each=0.0,
        panel_offset_y=0.2, ic=initial,
        prop_mass=initial_mass_kg === nothing ? Float64(vehicle["propellant_mass_kg"]) :
            Float64(initial_mass_kg) - Float64(vehicle["dry_mass_kg"]), id=1)
    atmosphere = fitted_atmosphere_parameters(config)
    use_gram = corridor_atmosphere_backend(config) == "gram"
    # The engine density model always extends its end layers; the drag effector
    # (the force actually integrated) honors the configured below/above behavior.
    nodes = atmosphere.node_altitudes_m
    density = use_gram ? _corridor_gram_density_model(config, epoch) :
        PiecewiseExponentialAtmosphereModel(
            nodes, atmosphere.node_density_kg_m3[1:end-1], atmosphere.scale_heights_m;
            temperature_k=atmosphere.temperature_k)
    drag = CorridorBallisticDrag(Float64(vehicle["ballistic_coefficient_kg_m2"]),
        Float64(vehicle["dry_mass_kg"]) + Float64(vehicle["propellant_mass_kg"]),
        atmosphere, use_gram ? :sampled : :fit)
    base = make_example_config(
        planet=planet, spacecraft=spacecraft,
        mission_time=run_horizon_s, initial_time=epoch,
        dynamic_effectors=(harmonics, drag), density_model=density,
        ephemerides_model=SimpleEphemeridesModel(), orientation_sim=false,
        keplerian=false,
        EI_km=use_gram ? Float64(config.values["atmosphere_gram"]["entry_interface_km"]) : 160.0,
        verbose=false, results=false,
        results_directory=results_dir)
    integ = config.values["integration"]
    tolerances = IntegrationTolerances(
        reltol=Float64(integ["reltol"]), abstol=Float64(integ["abstol"]),
        reltol_orbit=Float64(integ["reltol"]), abstol_orbit=Float64(integ["abstol"]),
        reltol_atmosphere=Float64(integ["reltol"]), abstol_atmosphere=Float64(integ["abstol"]),
        dt_max=Float64(integ["dt_max_orbit_s"]), dt_max_orbit=Float64(integ["dt_max_orbit_s"]),
        dt_max_atmosphere=Float64(integ["dt_max_atmosphere_s"]))
    guidance_model, control_model = if controller === nothing
        base.guidance_model, base.control_model
    else
        thruster = BaseThrusterModel([Float64(vehicle["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0],
            [Float64(vehicle["specific_impulse_s"])])
        GuidanceModel((controller,), [1.0e12]), ControlModel((thruster,), [1.0e12])
    end
    args = SimulationConfiguration(
        file_paths=base.file_paths, simulation_settings=base.simulation_settings,
        mission_configuration=base.mission_configuration, environment_model=base.environment_model,
        dynamics_model=base.dynamics_model, guidance_model=guidance_model,
        navigation_model=base.navigation_model, control_model=control_model,
        initial_time=base.initial_time, integration_tolerances=tolerances,
        solver_config=SolverConfig(solver_mode=Symbol(integ["solver_mode"])))
    return controller === nothing ? args : (args, only(control_model.control_effectors))
end

function _predict_next_physical_periapsis(config::CorridorConfig, phase::Symbol, source_args,
                                          r0::SVector{3, Float64}, v0::SVector{3, Float64},
                                          mass_kg::Float64, delta_v_mps::Float64;
                                          elapsed_s::Float64=0.0)
    planet = source_args.environment_model.planet
    oe = _corridor_element_summary(r0, v0, planet.μ)
    a, e = Float64(oe.a), Float64(oe.e)
    a > 0.0 && 0.0 <= e < 1.0 || throw(ArgumentError("predictor state is not elliptic"))
    period_s = 2π * sqrt(a^3 / planet.μ)
    fixed_guidance = CorridorFixedBurnGuidance(delta_v_mps, 1, false)
    # The predictor restarts at t = 0 from the current state. A time-dependent
    # atmosphere (GRAM) must still be sampled at the true date, so its epoch is
    # the run epoch plus the elapsed time; the US76 fit has no time dependence
    # and keeps the original epoch (behavior unchanged).
    epoch = corridor_atmosphere_backend(config) == "gram" ?
        _shift_initial_time(source_args.initial_time, elapsed_s) : _CORRIDOR_EPOCH
    args, _ = _corridor_args(config, phase, 1.2 * period_s,
        source_args.simulation_settings.results_directory;
        initial_state=(r0, v0), initial_mass_kg=mass_kg,
        controller=fixed_guidance, epoch=epoch)
    scheduled = Ref(false)
    altitude = Ref(NaN)
    geodetic_of = _geodetic_of(args)
    radial(u) = dot(u.sc[1].pos, u.sc[1].vel)
    radial_condition(u, t, integrator) = radial(u)
    function record_next_periapsis!(integrator)
        if scheduled[]
            altitude[] = geodetic_of(Float64(integrator.t), integrator.u.sc[1].pos)[1]
            terminate!(integrator)
        end
        return nothing
    end
    schedule_callback = _fixed_burn_schedule_callback(fixed_guidance, scheduled)
    peri_callback = ContinuousCallback(radial_condition, record_next_periapsis!, nothing;
        save_positions=(false, false))
    evals_before = _gram_density_evaluations()
    solution = run_simulation(args; return_solution=true, isolate_state=false,
        extra_callbacks=(schedule_callback, peri_callback))
    _GRAM_PREDICTOR_EVALS[] += _gram_density_evaluations() - evals_before
    isfinite(altitude[]) || throw(ErrorException("finite-thrust truth predictor did not reach the next physical periapsis (retcode=$(solution.retcode))"))
    return altitude[]
end

"""Earth radius (planet `Rp_e`) and GM (gravity coefficient file) used by the corridor runs."""
function _earth_radius_and_gm(config::CorridorConfig)
    planet = make_no_gram_planet(:earth)
    path = joinpath(REPO_ROOT, String(config.values["gravity"]["coefficient_file"]))
    return planet.Rp_e, GravitationalHarmonicsModel(3, 0, path, planet).gm_m3s2
end

function _calibrate_physical_initial_periapsis(config::CorridorConfig, phase::Symbol,
                                                results_dir::String; target_altitude_m::Float64=97_500.0)
    calibration = config.values["calibration"]
    tolerance_m = Float64(calibration["tolerance_m"])
    max_iterations = Int(calibration["max_iterations"])
    working = CorridorConfig(config.path, deepcopy(config.values))
    records = NamedTuple[]
    for iteration in 1:max_iterations
        orbit = working.values["orbit"]
        radius_m, gm_m3s2 = _earth_radius_and_gm(working)
        rp = Float64(orbit["initial_periapsis_altitude_m"]) + radius_m
        ra = Float64(orbit["apoapsis_altitude_m"]) + radius_m
        a0 = 0.5 * (rp + ra)
        period_s = 2π * sqrt(a0^3 / gm_m3s2)
        calibration_horizon_s = Float64(calibration["horizon_fraction_of_period"]) * period_s
        args = _corridor_args(working, phase, calibration_horizon_s, results_dir)
        solution = run_simulation(args; return_solution=true, isolate_state=false)
        events = physical_apsis_events(solution, _geodetic_of(args))
        periapsis = findfirst(event -> event.kind == "periapsis", events)
        periapsis === nothing && throw(ErrorException("physical initial-condition calibration did not reach a periapsis"))
        measured_altitude_m = events[periapsis].physical_altitude_m
        seed_altitude_m = Float64(orbit["initial_periapsis_altitude_m"])
        error_m = target_altitude_m - measured_altitude_m
        push!(records, (iteration=iteration, seed_osculating_periapsis_altitude_m=seed_altitude_m,
            measured_physical_periapsis_altitude_m=measured_altitude_m,
            target_physical_periapsis_altitude_m=target_altitude_m,
            residual_m=measured_altitude_m - target_altitude_m,
            calibration_elapsed_s=solution.t[end]))
        abs(error_m) <= tolerance_m && return working, records
        orbit["initial_periapsis_altitude_m"] = seed_altitude_m + error_m
    end
    last_residual = records[end].residual_m
    throw(ErrorException("physical initial-periapsis calibration did not converge within $(tolerance_m) m; residual=$(last_residual) m"))
end

function _corridor_element_summary(r::SVector{3, Float64}, v::SVector{3, Float64}, mu::Float64)
    rmag, v2 = norm(r), dot(v, v)
    hvec = cross(r, v)
    hmag = norm(hvec)
    energy = 0.5 * v2 - mu / rmag
    a = -mu / (2.0 * energy)
    evec = cross(v, hvec) / mu - r / rmag
    e = norm(evec)
    i = acos(clamp(hvec[3] / hmag, -1.0, 1.0))
    node = cross(SVector(0.0, 0.0, 1.0), hvec)
    nmag = norm(node)
    raan = nmag > 1e-12 ? mod(atan(node[2], node[1]), 2π) : 0.0
    omega = if nmag > 1e-12 && e > 1e-12
        mod(atan(dot(cross(node, evec), hvec) / (nmag * e * hmag), dot(node, evec) / (nmag * e)), 2π)
    else
        0.0
    end
    nu = e > 1e-12 ? mod(atan(dot(cross(evec, r), hvec) / (e * rmag * hmag), dot(evec, r) / (e * rmag)), 2π) : 0.0
    return (; a, e, i, raan, omega, nu, periapsis_radius_m=a * (1.0 - e))
end

"""
Cross-check of the geodetic helpers against the engine's own routines at the
equator, 50 degrees and the pole (geodetic altitude 97.5 km): the helper must
equal the engine's `rtolatlong`, invert the engine's ellipsoid constructor, and
the radius round trip must close. Returns the largest altitude discrepancy (m);
throws above 1 mm.
"""
function _geodetic_helper_crosscheck(planet)
    ctor = PeriapsisCorridorSupport._geodetic_position
    worst = 0.0
    for lat_deg in (0.0, 50.0, 90.0)
        lat = deg2rad(lat_deg)
        r_pp = ctor(lat, 0.3, 97_500.0, planet)
        h, _ = geodetic_altitude_latitude(r_pp, planet)
        engine = SimulationModel.SimulationCallbacks.rtolatlong(r_pp, planet)[1]
        worst = max(worst, abs(h - engine), abs(h - 97_500.0),
            abs(geodetic_altitude_at_radius(r_pp, geodetic_radius_m(97_500.0, lat, planet), planet) - 97_500.0))
    end
    worst <= 1e-3 || throw(ErrorException("geodetic helper disagrees with the engine routine by $worst m"))
    return worst
end

"""
`(min, max)` altitude band the passes are checked against. Every corridor case
uses the configured geodetic corridor. The vacuum case has no corridor; its
diagnostic band stays RADIAL (spherical radius above `Rp_e`, centered on the mean
periapsis radius) because it validates the J2 mean-to-osculating map in radius,
and a geodetic band would add the periapsis-latitude drift (up to ~21 km of
ellipsoid height) that has nothing to do with the map.
"""
function _band_limits(config::CorridorConfig, case_name, planet)
    if case_name == "vacuum"
        vacuum = config.values["vacuum_orbit"]
        center_m = Float64(vacuum["mean_semi_major_axis_m"]) *
            (1.0 - Float64(vacuum["mean_eccentricity"])) - planet.Rp_e
        half_width_m = Float64(vacuum["diagnostic_band_half_width_m"])
        return center_m - half_width_m, center_m + half_width_m
    end
    mission = config.values["mission"]
    return Float64(mission["corridor_min_altitude_m"]), Float64(mission["corridor_max_altitude_m"])
end

function _pass_rows(events, args, config::CorridorConfig, gravity_model, case_name;
                    controller=nothing, burn_rows=NamedTuple[])
    planet = args.environment_model.planet
    env = args.environment_model
    atmosphere = fitted_atmosphere_parameters(config)
    backend = corridor_atmosphere_backend(config)
    mission = config.values["mission"]
    corridor_min_m, corridor_max_m = _band_limits(config, case_name, planet)
    rows = NamedTuple[]
    previous_context = nothing
    for event in events
        event.kind == "periapsis" || continue
        r = event.position_ii_m
        v = event.velocity_ii_mps
        et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model)
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + event.time_s, env.ephemerides_model)
        velocity_pp = lpi * v - cross(planet.ω, lpi * r)
        fit_density = corridor_density(atmosphere, event.physical_altitude_m)
        # Altitude and latitude are the event's own (geodetic, engine routine);
        # longitude is only needed to sample GRAM at the same point.
        geodetic_alt_m, latitude_rad = event.physical_altitude_m, event.latitude_rad
        longitude_rad = SimulationModel.SimulationCallbacks.rtolatlong(lpi * r, planet)[3]
        density = backend == "gram" ?
            _gram_density_at(args, geodetic_alt_m, latitude_rad, longitude_rad, event.time_s) : fit_density
        wet_mass = Float64(config.values["vehicle"]["dry_mass_kg"]) +
            Float64(config.values["vehicle"]["propellant_mass_kg"])
        beta = Float64(config.values["vehicle"]["ballistic_coefficient_kg_m2"]) * event.mass_kg / wet_mass
        drag = drag_context_ratio(density, norm(velocity_pp), beta,
                                  gravity_model.gm_m3s2, event.radius_m)
        elements = _corridor_element_summary(r, v, gravity_model.gm_m3s2)
        j2 = -sqrt(5.0) * gravity_model.C[3, 1]
        j3 = -sqrt(7.0) * gravity_model.C[4, 1]
        # Radius-space corridor: the geodetic band at this event's latitude (a
        # latitude-dependent radius band). The vacuum diagnostic band is radial.
        rp_min_radius_m, rp_max_radius_m = case_name == "vacuum" ?
            (planet.Rp_e + corridor_min_m, planet.Rp_e + corridor_max_m) :
            corridor_radius_band(latitude_rad, corridor_min_m, corridor_max_m, planet)
        context = dimensionless_context(elements.a, elements.e, elements.i, j2, j3,
            gravity_model.reference_radius_m, planet.Rp_e,
            rp_min_radius_m, rp_max_radius_m)
        apo_event_idx = findfirst(candidate -> candidate.kind == "apoapsis" && candidate.time_s > event.time_s, events)
        if apo_event_idx === nothing
            apo_event_idx = findlast(candidate -> candidate.kind == "apoapsis" && candidate.time_s < event.time_s, events)
        end
        physical_apoapsis_altitude = apo_event_idx === nothing ? NaN : events[apo_event_idx].physical_altitude_m
        mean_diag = try
            result = osculating_to_mean((a=elements.a, e=elements.e, i=elements.i,
                raan=elements.raan, omega=elements.omega, nu=elements.nu);
                mu=gravity_model.gm_m3s2, radius=gravity_model.reference_radius_m, j2=j2)
            result.converged ? result : nothing
        catch
            nothing
        end
        mean_elements = mean_diag === nothing ? nothing : mean_diag.elements
        # Vacuum diagnostic: radius at mean true anomaly 0 mapped to osculating
        # by the J2 short-period map, versus the physical radial-event radius.
        map_periapsis_radius_m = NaN
        if case_name == "vacuum" && mean_elements !== nothing
            try
                osc = mean_to_osculating(merge(mean_elements, (nu=0.0,));
                    mu=gravity_model.gm_m3s2, radius=gravity_model.reference_radius_m, j2=j2).elements
                map_periapsis_radius_m = osc.a * (1.0 - osc.e^2) / (1.0 + osc.e * cos(osc.nu))
            catch
            end
        end
        decision = controller === nothing ? nothing : findfirst(
            item -> item.target_pass == length(rows) + 1, controller.decisions)
        decision_row = decision === nothing ? nothing : controller.decisions[decision]
        pass_burns = filter(item -> item.target_pass <= length(rows) + 1, burn_rows)
        cumulative_dv = sum((item.delivered_delta_v_mps for item in pass_burns); init=0.0)
        delta_pi_j2 = previous_context === nothing ? NaN : context.Pi_J2 - previous_context.Pi_J2
        delta_pi_j3 = previous_context === nothing ? NaN : context.Pi_J3 - previous_context.Pi_J3
        previous_context = context
        push!(rows, (pass=length(rows) + 1, time_s=event.time_s,
            physical_periapsis_altitude_m=event.physical_altitude_m,
            osculating_periapsis_altitude_m=geodetic_altitude_at_radius(
                lpi * r, elements.periapsis_radius_m, planet),
            spherical_periapsis_altitude_m=event.radius_m - planet.Rp_e,
            corridor_radius_min_m=rp_min_radius_m, corridor_radius_max_m=rp_max_radius_m,
            physical_periapsis_radius_m=event.radius_m,
            map_periapsis_radius_m=map_periapsis_radius_m,
            map_minus_physical_periapsis_radius_m=map_periapsis_radius_m - event.radius_m,
            physical_apoapsis_altitude_m=physical_apoapsis_altitude,
            corridor_coordinate=corridor_coordinate(
                case_name == "vacuum" ? event.radius_m - planet.Rp_e : event.physical_altitude_m,
                corridor_min_m, corridor_max_m),
            a_m=elements.a, e=elements.e, i_rad=elements.i, raan_rad=elements.raan,
            omega_rad=elements.omega, nu_rad=elements.nu, p_m=context.p_m,
            mean_a_m=mean_elements === nothing ? NaN : mean_elements.a,
            mean_e=mean_elements === nothing ? NaN : mean_elements.e,
            mean_i_rad=mean_elements === nothing ? NaN : mean_elements.i,
            mean_omega_rad=mean_elements === nothing ? NaN : mean_elements.omega,
            mean_k=mean_elements === nothing ? NaN : mean_elements.e * cos(mean_elements.omega),
            mean_h=mean_elements === nothing ? NaN : mean_elements.e * sin(mean_elements.omega),
            mean_elements_valid=mean_diag !== nothing,
            mean_coordinate_converged=mean_diag !== nothing,
            mean_model_certifies_corridor=mean_diag !== nothing && mean_diag.corridor_safety_certified,
            reference_radius_over_a=context.reference_radius_over_a,
            scale_height_over_periapsis=backend == "gram" ? NaN :
                corridor_scale_height(atmosphere, event.physical_altitude_m) / event.radius_m,
            atmosphere_backend=backend, density_kg_m3=density,
            us76_fit_density_at_state_kg_m3=fit_density,
            density_over_fit_density=density / fit_density,
            geodetic_altitude_m=geodetic_alt_m, latitude_rad=latitude_rad, longitude_rad=longitude_rad,
            Pi_J2=context.Pi_J2, Pi_J3=context.Pi_J3,
            delta_Pi_J2=delta_pi_j2, delta_Pi_J3=delta_pi_j3,
            signed_Pi_J3=context.signed_Pi_J3, J3_sign=context.J3_sign,
            Gamma_3B_J2=context.Gamma_3B_J2, Pi_drag=drag.Pi_drag,
            dynamic_pressure_pa=drag.dynamic_pressure_pa,
            predicted_periapsis_altitude_m=decision_row === nothing ? NaN : decision_row.predicted_periapsis_altitude_m,
            decision=decision_row === nothing ? "passive_coast" : decision_row.trigger_reason,
            mass_kg=event.mass_kg, beta_current_kg_m2=beta,
            cumulative_delivered_delta_v_mps=cumulative_dv))
    end
    return rows
end

function _actual_burn_rows(sol, policy, config, events)
    policy === nothing && return NamedTuple[]
    isp = Float64(config.values["vehicle"]["specific_impulse_s"])
    dry_mass = Float64(config.values["vehicle"]["dry_mass_kg"])
    rows = NamedTuple[]
    for plan in policy.burn_plans
        start, stop = plan.execution_start_time_s, plan.execution_stop_time_s
        scheduled = isfinite(start) && isfinite(stop) && stop > start
        started = scheduled && start < sol.t[end] - 1e-8
        complete = scheduled && stop <= sol.t[end] + 1e-8
        # Planned before the horizon end but starting at or after it: not executed within the horizon.
        beyond_horizon = scheduled && !started
        mass_start = started ? try
            Float64(sol(start).sc[1].mass)
        catch
            NaN
        end : NaN
        actual_stop = started ? min(stop, sol.t[end]) : NaN
        mass_stop = started ? try
            Float64(sol(actual_stop).sc[1].mass)
        catch
            NaN
        end : NaN
        used = started && isfinite(mass_start) && isfinite(mass_stop) ? mass_start - mass_stop : 0.0
        delivered = started && isfinite(mass_start) && isfinite(mass_stop) && mass_start > mass_stop > 0.0 ?
            isp * 9.80665 * log(mass_start / mass_stop) : 0.0
        achieved = findfirst(event -> event.kind == "periapsis" &&
            count(candidate -> candidate.kind == "periapsis" && candidate.time_s <= event.time_s, events) == plan.target_pass,
            events)
        achieved_altitude = achieved === nothing ? NaN : events[achieved].physical_altitude_m
        push!(rows, merge(plan, (delivered_delta_v_mps=delivered,
            propellant_used_kg=used, mass_at_start_kg=mass_start,
            mass_at_stop_kg=mass_stop, burn_started=started, burn_beyond_horizon=beyond_horizon,
            execution_complete=complete,
            execution_stop_actual_time_s=actual_stop,
            achieved_periapsis_altitude_m=achieved_altitude,
            propellant_remaining_after_burn_kg=isfinite(mass_stop) ? mass_stop - dry_mass : NaN)))
    end
    return rows
end

function _write_case_outputs(output_dir, case_name, args, sol, config, gravity_model,
        requested_horizon_s, callback_events; smoke=false, calibration_records=NamedTuple[],
        controller=nothing, run_metrics=NamedTuple(), atmosphere_provenance=Dict{String, Any}())
    mkpath(output_dir)
    events = physical_apsis_events(sol, _geodetic_of(args))
    burn_rows = _actual_burn_rows(sol, controller, config, events)
    rows = _pass_rows(events, args, config, gravity_model, case_name;
        controller=controller, burn_rows=burn_rows)
    passes_path = joinpath(output_dir, "passes.csv")
    pass_columns = [:pass, :time_s, :physical_periapsis_altitude_m, :osculating_periapsis_altitude_m,
        :physical_periapsis_radius_m, :map_periapsis_radius_m, :map_minus_physical_periapsis_radius_m,
        :physical_apoapsis_altitude_m, :corridor_coordinate, :a_m, :e, :i_rad, :raan_rad,
        :omega_rad, :nu_rad, :p_m, :mean_a_m, :mean_e, :mean_i_rad, :mean_omega_rad,
        :mean_k, :mean_h, :mean_elements_valid, :mean_coordinate_converged,
        :mean_model_certifies_corridor,
        :reference_radius_over_a, :scale_height_over_periapsis,
        :Pi_J2, :Pi_J3, :delta_Pi_J2, :delta_Pi_J3, :signed_Pi_J3,
        :J3_sign, :Gamma_3B_J2, :Pi_drag, :dynamic_pressure_pa, :predicted_periapsis_altitude_m,
        :decision, :mass_kg, :cumulative_delivered_delta_v_mps,
        :atmosphere_backend, :density_kg_m3, :us76_fit_density_at_state_kg_m3, :density_over_fit_density,
        :spherical_periapsis_altitude_m, :corridor_radius_min_m, :corridor_radius_max_m, :geodetic_altitude_m, :latitude_rad, :longitude_rad]
    pass_df = isempty(rows) ? DataFrame([name => Float64[] for name in pass_columns]) : DataFrame(rows)
    CSV.write(passes_path, pass_df)
    burn_columns = [:request_time_s, :target_pass, :execution_start_time_s,
        :execution_stop_time_s, :target_periapsis_altitude_m,
        :predicted_periapsis_altitude_m, :requested_delta_v_mps,
        :delivered_delta_v_mps, :propellant_used_kg, :mass_at_start_kg,
        :mass_at_stop_kg, :burn_started, :burn_beyond_horizon, :execution_complete,
        :execution_stop_actual_time_s, :achieved_periapsis_altitude_m,
        :propellant_remaining_after_burn_kg,
        :trigger_reason]
    burn_df = isempty(burn_rows) ? DataFrame([name => Any[] for name in burn_columns]) : DataFrame(burn_rows)
    CSV.write(joinpath(output_dir, "burns.csv"), burn_df)
    # Altitude compared with the band: geodetic for every corridor case, spherical
    # radius above Rp_e for the vacuum diagnostic band (see `_band_limits`).
    band_altitude(row) = case_name == "vacuum" ? row.spherical_periapsis_altitude_m :
        row.physical_periapsis_altitude_m
    min_peri = isempty(rows) ? NaN : minimum(band_altitude(row) for row in rows)
    max_peri = isempty(rows) ? NaN : maximum(band_altitude(row) for row in rows)
    mission = config.values["mission"]
    corridor_min_m, corridor_max_m = _band_limits(config, case_name, args.environment_model.planet)
    completed = sol.t[end] + 1e-6 >= requested_horizon_s
    retcode_success = string(sol.retcode) == "Success"
    breach_count = count(row -> band_altitude(row) < corridor_min_m ||
                                 band_altitude(row) > corridor_max_m, rows)
    a0 = if case_name == "vacuum"
        Float64(config.values["vacuum_orbit"]["mean_semi_major_axis_m"])
    else
        ra = Float64(config.values["orbit"]["apoapsis_altitude_m"]) + args.environment_model.planet.Rp_e
        rp = Float64(config.values["orbit"]["initial_periapsis_altitude_m"]) + args.environment_model.planet.Rp_e
        0.5 * (ra + rp)
    end
    period0 = 2π * sqrt(a0^3 / gravity_model.gm_m3s2)
    expected_events = max(1, floor(Int, requested_horizon_s / period0) - 1)
    callback_periapsides = filter(event -> event.kind == "periapsis", callback_events)
    callback_periapsis_count = length(callback_periapsides)
    event_time_errors = Float64[]
    for event in events
        event.kind == "periapsis" || continue
        isempty(callback_periapsides) && continue
        push!(event_time_errors, minimum(abs(event.time_s - observed.time_s) for observed in callback_periapsides))
    end
    max_event_time_error_s = isempty(event_time_errors) ? Inf : maximum(event_time_errors)
    event_coverage_complete = length(rows) == callback_periapsis_count && length(rows) >= expected_events &&
        max_event_time_error_s <= 1.0
    controller_success = controller === nothing ||
        isempty(controller.failures) && all(decision -> decision.feasible, controller.decisions) &&
        all(row -> row.burn_started || row.burn_beyond_horizon, burn_rows)
    full_mission_verified = case_name != "vacuum" && sol.t[end] + 1e-6 >= 5_184_000.0 && retcode_success &&
        event_coverage_complete && breach_count == 0 && controller_success
    first_breach_time_s = findfirst(row -> band_altitude(row) < corridor_min_m ||
        band_altitude(row) > corridor_max_m, rows)
    first_breach_time_s = first_breach_time_s === nothing ? NaN : rows[first_breach_time_s].time_s
    minimum_margin_m = isempty(rows) ? NaN : minimum(min(
        band_altitude(row) - corridor_min_m,
        corridor_max_m - band_altitude(row)) for row in rows)
    summary = (case=case_name, requested_horizon_s=requested_horizon_s,
        elapsed_horizon_s=sol.t[end], full_60day_objective=5_184_000.0,
        horizon_completed=completed, solver_retcode=string(sol.retcode), retcode_success=retcode_success,
        expected_periapsis_events=expected_events, callback_periapsis_events=callback_periapsis_count,
        max_event_time_error_s=max_event_time_error_s,
        event_coverage_complete=event_coverage_complete, full_mission_verified=full_mission_verified,
        controller_success=controller_success,
        control_policy=controller === nothing ? "passive_coast" : String(controller.policy),
        controller_failure_count=controller === nothing ? 0 : length(controller.failures),
        controller_failures=controller === nothing ? "" : join(controller.failures, " | "),
        controller_decision_count=controller === nothing ? 0 : length(controller.decisions),
        first_physical_breach_time_s=first_breach_time_s,
        minimum_physical_corridor_margin_m=minimum_margin_m,
        physical_pass_fraction_inside_corridor=isempty(rows) ? 0.0 :
            (length(rows) - breach_count) / length(rows),
        smoke=smoke, physical_periapsis_events=length(rows),
        geodetic_helper_crosscheck_max_diff_m=_geodetic_helper_crosscheck(args.environment_model.planet),
        corridor_altitude_definition=case_name == "vacuum" ?
            "vacuum diagnostic band: spherical radius above Earth.Rp_e" :
            "geodetic height above the planet ellipsoid",
        diagnostic_corridor_min_altitude_m=corridor_min_m,
        diagnostic_corridor_max_altitude_m=corridor_max_m,
        physical_breach_count=breach_count, minimum_physical_periapsis_altitude_m=min_peri,
        maximum_physical_periapsis_altitude_m=max_peri,
        maintenance_burn_count=count(row -> row.burn_started, burn_rows),
        burns_beyond_horizon_count=count(row -> row.burn_beyond_horizon, burn_rows),
        completed_burn_count=count(row -> row.execution_complete, burn_rows),
        initial_periapsis_calibration_iterations=length(calibration_records),
        initial_periapsis_calibration_final_residual_m=isempty(calibration_records) ? NaN : calibration_records[end].residual_m,
        delivered_delta_v_mps=sum((row.delivered_delta_v_mps for row in burn_rows
                                   if isfinite(row.delivered_delta_v_mps)); init=0.0),
        propellant_remaining_kg=Float64(sol.u[end].sc[1].mass) -
            Float64(config.values["vehicle"]["dry_mass_kg"]),
        mission_result=(smoke ? "smoke_only" : case_name == "vacuum" ? "vacuum_diagnostic_only" :
            !retcode_success || !completed ? "truncated_failure" :
            !controller_success ? "controller_execution_failure" :
            full_mission_verified ? "physical_corridor_pass" :
            sol.t[end] < 5_184_000.0 ? "reduced_horizon_diagnostic" : "physical_corridor_failure"))
    summary = merge(summary, (atmosphere_backend=corridor_atmosphere_backend(config),), run_metrics)
    summary_df = DataFrame([summary])
    CSV.write(joinpath(output_dir, "summary.csv"), summary_df)
    CSV.write(joinpath(output_dir, "comparison.csv"), summary_df)
    write_manifest(joinpath(output_dir, "run_manifest.toml"), config;
        case_name=case_name, gravity_path=joinpath(REPO_ROOT, String(config.values["gravity"]["coefficient_file"])),
        solver_settings=Dict("reltol" => args.integration_tolerances.reltol_orbit,
            "abstol" => args.integration_tolerances.abstol_orbit,
            "dt_max_orbit_s" => args.integration_tolerances.dt_max_orbit,
            "dt_max_atmosphere_s" => args.integration_tolerances.dt_max_atmosphere,
            "solver_mode" => string(args.solver_config.solver_mode)),
        initial_calibration=calibration_records,
        actual_initial_state=(position_ii_m=collect(sol.u[1].sc[1].pos),
            velocity_ii_mps=collect(sol.u[1].sc[1].vel), mass_kg=Float64(sol.u[1].sc[1].mass)),
        source_files=(joinpath(@__DIR__, "AGORA_Earth_Periapsis_Corridor.jl"),
            joinpath(@__DIR__, "periapsis_corridor_support.jl"),
            joinpath(@__DIR__, "periapsis_corridor_mean_elements.jl"),
            joinpath(@__DIR__, "periapsis_corridor_controller.jl"),
            joinpath(@__DIR__, "periapsis_corridor_plots.jl"),
            joinpath(@__DIR__, "periapsis_corridor_config.toml"),
            joinpath(REPO_ROOT, "Project.toml"), joinpath(REPO_ROOT, "Manifest.toml")),
        atmosphere_provenance=atmosphere_provenance,
        run_status=full_mission_verified ? "full_mission_verified" : completed ? "requested_horizon_completed" : "truncated")
    return summary
end

function run_corridor_case(config::CorridorConfig, case_name::AbstractString,
                           output_root::AbstractString; kwargs...)
    corridor_atmosphere_backend(config) == "gram" || return _run_corridor_case(config, case_name, output_root; kwargs...)
    # The extension may load now, so run in the latest world age (otherwise the
    # engine cannot see the extension's methods). The profile counters give the
    # GRAM density-evaluation counts in the summary.
    case_name == "vacuum" && throw(ArgumentError(
        "the vacuum case has no atmosphere; atmosphere backend \"gram\" does not apply to it"))
    _corridor_gram_load!()
    return withenv("SPACEAGORA_GRAM_PROFILE" => "1") do
        Base.invokelatest(_run_corridor_case, config, case_name, output_root; kwargs...)
    end
end

function _run_corridor_case(config::CorridorConfig, case_name::AbstractString,
                            output_root::AbstractString; horizon_s::Union{Nothing, Float64}=nothing,
                            smoke::Bool=false, phase_override::Union{Nothing, Symbol}=nothing)
    case_name in ("vacuum", "unfavorable", "favorable", "near_passive", "frequent") ||
        throw(ArgumentError("unknown implemented case $case_name"))
    phase = case_name == "vacuum" ? :vacuum :
        case_name in ("unfavorable", "near_passive", "frequent") ? :unfavorable : :favorable
    # The sweep driver passes :favorable to run a controller case at the configured
    # argument of periapsis instead of the unfavorable case's forced 0 degrees.
    phase_override === nothing || (phase = phase_override)
    backend = corridor_atmosphere_backend(config)
    backend == "gram" && case_name == "vacuum" && throw(ArgumentError(
        "the vacuum case has no atmosphere; atmosphere backend \"gram\" does not apply to it"))
    case_config = CorridorConfig(config.path, deepcopy(config.values))
    phase === :unfavorable && (case_config.values["orbit"]["argument_of_periapsis_deg"] = 0.0)
    controller_policy = case_name == "near_passive" ? :near_passive :
        case_name == "frequent" ? :frequent : nothing
    atmosphere = case_name == "vacuum" ? NoAtmosphereModel() : nothing
    # The first implementation slice keeps the same sourced gravity/orbit setup
    # while vacuum and aero configurations use separate runs and output trees.
    full_horizon = Float64(case_config.values["mission"]["horizon_s"])
    requested = something(horizon_s, full_horizon)
    run_horizon = if smoke
        radius_m, gm_m3s2 = _earth_radius_and_gm(case_config)
        if phase === :vacuum
            a = Float64(case_config.values["vacuum_orbit"]["mean_semi_major_axis_m"])
        else
            rp = Float64(case_config.values["orbit"]["initial_periapsis_altitude_m"]) + radius_m
            ra = Float64(case_config.values["orbit"]["apoapsis_altitude_m"]) + radius_m
            a = 0.5 * (ra + rp)
        end
        Float64(case_config.values["comparison"]["smoke_orbits"]) * 2π * sqrt(a^3 / gm_m3s2)
    else
        requested
    end
    case_dir = joinpath(output_root, case_name)
    _assert_case_dir_backend(case_dir, backend)
    calibration_records = NamedTuple[]
    calibration_start_ns = time_ns()
    if phase !== :vacuum
        case_config, calibration_records = _calibrate_physical_initial_periapsis(
            case_config, phase, case_dir;
            target_altitude_m=Float64(case_config.values["mission"]["reference_altitude_m"]))
    end
    calibration_wall_s = (time_ns() - calibration_start_ns) / 1e9
    comparison = case_config.values["comparison"]
    controller = controller_policy === nothing ? nothing : CorridorGuidancePolicy(
        case_config, phase, controller_policy,
        Float64(comparison["max_burn_delta_v_mps"]),
        Float64(comparison["max_burn_duration_s"]),
        Float64(comparison["total_delta_v_budget_mps"]),
        Float64(controller_policy === :frequent ? comparison["frequent_deadband_m"] : comparison["near_passive_deadband_m"]),
        Float64(comparison["target_periapsis_altitude_m"]),
        deg2rad(Float64(comparison["planning_lead_deg"])), 0, false, nothing,
        NamedTuple[], NamedTuple[], String[])
    built = _corridor_args(case_config, phase, run_horizon, case_dir; controller=controller)
    args, thruster = controller === nothing ? (built, nothing) : built
    if case_name == "vacuum"
        vacuum_effectors = Tuple(effector for effector in args.dynamics_model.dynamic_effectors
                                 if !(effector isa CorridorBallisticDrag))
        args = SimulationConfiguration(
            file_paths=args.file_paths, simulation_settings=args.simulation_settings,
            mission_configuration=args.mission_configuration, environment_model=EnvironmentModel(
                planet=args.environment_model.planet, EI=args.environment_model.EI,
                density_model=NoAtmosphereModel(), thermal_model=args.environment_model.thermal_model,
                ephemerides_model=args.environment_model.ephemerides_model,
                topography=false, wind=false),
            dynamics_model=DynamicsModel(args.dynamics_model.spacecraft, vacuum_effectors),
            guidance_model=args.guidance_model,
            navigation_model=args.navigation_model, control_model=args.control_model,
            initial_time=args.initial_time, integration_tolerances=args.integration_tolerances,
            solver_config=args.solver_config)
    end
    callback_events = NamedTuple[]
    apsis_callback = _apsis_observation_callback(callback_events, args; controller=controller)
    extra_callbacks = controller === nothing ? (apsis_callback,) :
        (apsis_callback, _preapoapsis_guidance_callback(controller))
    provenance = backend == "gram" ?
        _corridor_gram_provenance(case_config, args.environment_model.density_model, _CORRIDOR_EPOCH) :
        Dict{String, Any}("backend" => backend)
    backend == "gram" && SimulationModel.SimulationCallbacks._gram_runtime_stats_reset!()
    _GRAM_PREDICTOR_EVALS[] = 0
    run_start_ns = time_ns()
    sol = run_simulation(args; return_solution=true, isolate_state=false,
                         extra_callbacks=extra_callbacks)
    run_wall_s = (time_ns() - run_start_ns) / 1e9
    gram_evals = backend == "gram" ? _gram_density_evaluations() : 0
    # Epoch the engine actually used for GRAM (it may rebuild the model at the run epoch).
    if backend == "gram"
        provenance["gram_epoch_used_by_engine_utc"] = try
            _utc_string(sol.prob.p.args.environment_model.density_model.core.initial_time)
        catch
            "unavailable"
        end
    end
    run_metrics = (calibration_wall_s=calibration_wall_s, run_wall_s=run_wall_s,
        gram_density_evaluations=gram_evals,
        gram_predictor_density_evaluations=backend == "gram" ? _GRAM_PREDICTOR_EVALS[] : 0,
        predictor_evaluations=controller === nothing ? 0 :
            sum((d.predictor_evaluations for d in controller.decisions); init=0))
    gravity = only(filter(eff -> eff isa GravitationalHarmonicsModel,
                          args.dynamics_model.dynamic_effectors))
    summary = _write_case_outputs(case_dir, case_name, args, sol, case_config, gravity,
        run_horizon, callback_events; smoke=smoke, calibration_records=calibration_records,
        controller=controller, run_metrics=run_metrics, atmosphere_provenance=provenance)
    println("case=$(case_name) atmosphere=$(backend) wall_s=$(round(run_wall_s; digits=2)) elapsed_s=$(summary.elapsed_horizon_s) periapsis_events=$(summary.physical_periapsis_events) breaches=$(summary.physical_breach_count) result=$(summary.mission_result)")
    smoke && println("Smoke mode completed; it does not evaluate 60-day passive performance.")
    return summary
end

function main(args=ARGS)
    cli = _corridor_cli(args)
    config = load_corridor_config(cli.config_path)
    cli.atmosphere === nothing || (config.values["atmosphere"]["backend"] = cli.atmosphere)
    corridor_atmosphere_backend(config) == "gram" && !haskey(config.values, "atmosphere_gram") &&
        throw(ArgumentError("backend \"gram\" requires [atmosphere_gram] in the configuration"))
    horizon = cli.horizon_override
    cases = cli.case_name == "all" ?
        ("vacuum", "favorable", "unfavorable", "near_passive", "frequent") : (cli.case_name,)
    # The vacuum case has no atmosphere, so a GRAM comparison sweep skips it.
    corridor_atmosphere_backend(config) == "gram" && (cases = filter(!=("vacuum"), cases))
    summaries = map(cases) do case_name
        run_corridor_case(config, case_name, cli.output_dir; horizon_s=horizon, smoke=cli.smoke)
    end
    PeriapsisCorridorPlots.write_corridor_plots(cli.output_dir, config)
    CSV.write(joinpath(cli.output_dir, "comparison.csv"), DataFrame(summaries))
    required_cases = cli.case_name == "all" ?
        ("favorable", "near_passive", "frequent") :
        cli.case_name in ("favorable", "near_passive", "frequent") ? (cli.case_name,) : ()
    if !cli.smoke && something(horizon, Float64(config.values["mission"]["horizon_s"])) >= 5_184_000.0
        failed = filter(summary -> summary.case in required_cases && !summary.full_mission_verified, summaries)
        isempty(failed) || throw(ErrorException("required full-horizon corridor case(s) failed; artifacts were written: " *
            join(("$(summary.case)=$(summary.mission_result)" for summary in failed), ", ")))
    end
    return summaries
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
