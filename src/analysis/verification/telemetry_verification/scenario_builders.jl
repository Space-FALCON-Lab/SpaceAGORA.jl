@inline function _period_seconds(planet, ra::Float64, rp::Float64)::Float64
    a = 0.5 * (ra + rp)
    return 2π * sqrt(a^3 / planet.μ)
end

@inline function _planet_from_name(planet_name::String)
    key = lowercase(strip(planet_name))
    if key == "mars"
        return Mars("", SPICE_PATH)
    elseif key == "venus"
        return Venus("", SPICE_PATH)
    elseif key == "earth"
        return Earth("", SPICE_PATH)
    elseif key == "moon"
        return Moon("", SPICE_PATH)
    end
    throw(ArgumentError("Unsupported planet '$planet_name' in telemetry benchmark."))
end

@inline function _base_gravity_effector(gravity_model::Symbol)
    if gravity_model == :inverse_squared
        return InverseSquaredGravityModel()
    elseif gravity_model == :inverse_squared_j2
        return InverseSquaredJ2GravityModel()
    end
    throw(ArgumentError("Unsupported gravity model '$gravity_model'"))
end

@inline function _harmonics_order(degree::Int, order::Int)::Int
    if degree <= 0
        return 0
    end
    if order == 0
        return 0
    end
    if order < 0
        return degree
    end
    return min(order, degree)
end

@inline function _nbody_primary_name(planet_name::String)::String
    key = lowercase(strip(planet_name))
    if key == "earth"
        return "Earth"
    elseif key == "mars"
        return "Mars"
    elseif key == "venus"
        return "Venus"
    elseif key == "moon"
        return "Moon"
    elseif key == "titan"
        return "Titan"
    end
    throw(ArgumentError("Unsupported N-body primary planet '$planet_name'"))
end

@inline function _telemetry_j2_source_for_scenario(scenario_name::String)::Symbol
    default_raw = lowercase(strip(get(ENV, "SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT", "file_c20")))
    default_source = default_raw in ("planet", "planet_j2") ? :planet_j2 : :file_c20

    scenario_set_raw = strip(get(ENV, "SPACEAGORA_TELEMETRY_J2_SOURCE_PLANET_SCENARIOS", ""))
    isempty(scenario_set_raw) && return default_source

    names = Set{String}()
    for tok in split(scenario_set_raw, ',')
        t = lowercase(strip(tok))
        isempty(t) || push!(names, t)
    end
    return lowercase(strip(scenario_name)) in names ? :planet_j2 : default_source
end

@inline function _telemetry_coefficients_normalized_for_scenario(scenario_name::String)::Bool
    default_raw = lowercase(strip(get(ENV, "SPACEAGORA_TELEMETRY_HARMONICS_NORMALIZED_DEFAULT", "true")))
    default_normalized = !(default_raw in ("0", "false", "no", "off"))

    scenario_set_raw = strip(get(ENV, "SPACEAGORA_TELEMETRY_HARMONICS_UNNORMALIZED_SCENARIOS", ""))
    isempty(scenario_set_raw) && return default_normalized

    names = Set{String}()
    for tok in split(scenario_set_raw, ',')
        t = lowercase(strip(tok))
        isempty(t) || push!(names, t)
    end
    return lowercase(strip(scenario_name)) in names ? false : default_normalized
end

struct ScaledAerodynamicCoefficientfM <: AbstractForceTorqueModel
    model::AerodynamicCoefficientfM
    cd_scale::Float64
    function ScaledAerodynamicCoefficientfM(model::AerodynamicCoefficientfM, cd_scale::Float64)
        cd_scale > 0.0 || throw(ArgumentError("ScaledAerodynamicCoefficientfM.cd_scale must be > 0.0, got $cd_scale"))
        return new(model, cd_scale)
    end
end

@inline _dynamic_effector_threadsafe(::ScaledAerodynamicCoefficientfM)::Bool = true

function SimulationModel.calcForceTorque(
    model::ScaledAerodynamicCoefficientfM,
    x::AbstractVector,
    param::SimulationModel.ODEParams,
    i::Int64
)
    f, τ = SimulationModel.calcForceTorque(model.model, x, param, i)
    return model.cd_scale .* f, model.cd_scale .* τ
end

function _scenario_dynamic_effectors(
    cfg::AbstractScenarioConfig,
    planet,
    spacecraft;
    cd_scale::Float64=1.0,
    cr_override::Union{Nothing, Float64}=nothing
)
    effectors = Any[]

    if cfg.gravity_harmonics_degree > 0
        harmonics_file = cfg.gravity_harmonics_file
        isempty(harmonics_file) && throw(ArgumentError(
            "Scenario $(cfg.name) sets gravity_harmonics_degree=$(cfg.gravity_harmonics_degree) but does not provide gravity_harmonics_file."
        ))
        isfile(harmonics_file) || throw(ArgumentError("Missing gravity_harmonics_file for $(cfg.name): $harmonics_file"))
        push!(
            effectors,
            GravitationalHarmonicsModel(
                cfg.gravity_harmonics_degree,
                _harmonics_order(cfg.gravity_harmonics_degree, cfg.gravity_harmonics_order),
                harmonics_file,
                planet;
                coefficients_normalized=_telemetry_coefficients_normalized_for_scenario(cfg.name),
                j2_source=_telemetry_j2_source_for_scenario(cfg.name)
            )
        )
    else
        push!(effectors, _base_gravity_effector(cfg.gravity_model))
    end

    if !isempty(cfg.nbody_bodies)
        push!(
            effectors,
            NBodyGravityModel(
                body_names=Tuple(cfg.nbody_bodies),
                primary_body_name=_nbody_primary_name(cfg.planet_name),
                planet=planet
            )
        )
    end

    if cfg.srp_enabled
        area_m2 = cfg.srp_area_m2 > 0.0 ? cfg.srp_area_m2 : Float64(spacecraft.root.ref_area)
        area_m2 > 0.0 || throw(ArgumentError("Scenario $(cfg.name) SRP area must be > 0.0 m^2"))
        cr_value = cr_override === nothing ? cfg.srp_cr : Float64(cr_override)
        push!(effectors, SolarRadiationPressureModel(cr_value, area_m2))
    end

    if cfg.drag_enabled
        aero = AerodynamicCoefficientfM()
        if isapprox(cd_scale, 1.0; rtol=0.0, atol=1e-12)
            push!(effectors, aero)
        else
            push!(effectors, ScaledAerodynamicCoefficientfM(aero, cd_scale))
        end
    end
    return Tuple(effectors)
end

@inline function _scenario_density_model(cfg::AbstractScenarioConfig)
    cfg.drag_enabled || return SimulationModel.NoAtmosphereModel()
    if cfg.atmosphere_truth.atmosphere_model == "tabulated_flight"
        return _make_tabulated_flight_density_model(cfg.initial_time, cfg.atmosphere_truth)
    end
    return _make_required_gram_density_model(cfg.planet_name, cfg.initial_time, cfg.atmosphere_truth)
end

# Flight-measured density table (see build script in the lab notes): one row
# per (pass, leg, 1 km altitude bin) with mean density, its standard error,
# and the pass periapsis UTC. Loaded into per-pass in/out profiles keyed by
# elapsed time from the scenario epoch; nearest pass answers each query, so
# archive gaps fall back to the neighboring pass (counted and logged here).
function _make_tabulated_flight_density_model(
    initial_time::InitialTime,
    truth::AtmosphereTruthConfig
)
    path = truth.tabulated_flight_file
    isfile(path) || throw(ArgumentError("tabulated_flight_file not found: $path"))
    tbl = DataFrame(Arrow.Table(path))
    for col in ("P", "leg", "alt_km", "rho_kgm3", "sigma_kgm3", "t_peri_utc")
        hasproperty(tbl, Symbol(col)) || throw(ArgumentError("tabulated_flight_file missing column '$col'"))
    end
    epoch = DateTime(
        Int(initial_time.year), Int(initial_time.month), Int(initial_time.day),
        Int(initial_time.hour), Int(initial_time.minute)
    ) + Millisecond(round(Int, 1000 * Float64(initial_time.second)))
    pass_ids = sort(unique(Int.(tbl.P)))
    peri_el = Float64[]
    alt_pairs = NTuple{2, Vector{Float64}}[]
    log_pairs = NTuple{2, Vector{Float64}}[]
    sig_pairs = NTuple{2, Vector{Float64}}[]
    for pid in pass_ids
        sub = tbl[Int.(tbl.P) .== pid, :]
        t_peri = DateTime(first(sub.t_peri_utc)[1:19])
        push!(peri_el, Float64(Dates.value(t_peri - epoch)) / 1000.0)
        alts = (Float64[], Float64[]); logs = (Float64[], Float64[]); sigs = (Float64[], Float64[])
        for r in eachrow(sub)
            li = r.leg == "in" ? 1 : 2
            rho = Float64(r.rho_kgm3)
            rho > 0.0 || continue
            push!(alts[li], Float64(r.alt_km) * 1000.0)
            push!(logs[li], log(rho))
            push!(sigs[li], Float64(r.sigma_kgm3) / rho)   # sigma of log-density
        end
        for li in (1, 2)
            ord = sortperm(alts[li])
            permute!(alts[li], ord); permute!(logs[li], ord); permute!(sigs[li], ord)
        end
        push!(alt_pairs, alts); push!(log_pairs, logs); push!(sig_pairs, sigs)
    end
    ord = sortperm(peri_el)
    n_gaps = isempty(pass_ids) ? 0 : (maximum(pass_ids) - minimum(pass_ids) + 1 - length(pass_ids))
    println("tabulated_flight: $(length(pass_ids)) passes, $(nrow(tbl)) bins, " *
            "$(n_gaps) archive gaps (nearest-pass fallback), sigma_scale=$(truth.tabulated_flight_sigma)")
    return SimulationModel.TabulatedFlightAtmosphereModel(
        peri_el[ord], alt_pairs[ord], log_pairs[ord], sig_pairs[ord],
        truth.tabulated_flight_sigma, 3.4, 188.92
    )
end

Base.@kwdef struct _GRAMOfflineSurrogateFallbackBase
    planet_name::String
end

function SimulationModel.EnvironmentModels._gram_point_density(
    model::_GRAMOfflineSurrogateFallbackBase,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    # Library-less fallback cannot query native GRAM point model; return vacuum-like fallback.
    return 0.0, 200.0, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function _is_gram_library_missing_error(err)::Bool
    msg = lowercase(sprint(showerror, err))
    return occursin("gram shared library", msg) || occursin("gram_lib", msg)
end

@inline function _libraryless_gram_surrogate_enabled()::Bool
    return _safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_ALLOW_GRAM_OFFLINE_NO_LIB", "1"), true)
end

function _try_libraryless_gram_surrogate(planet_name::String, truth::AtmosphereTruthConfig)
    # gram_offline_surrogate="off" opts the scenario out of every surrogate
    # path, including this library-missing fallback: a benchmark pinned to the
    # native model must fail rather than silently fly the frozen-epoch grid.
    truth.gram_offline_surrogate == "off" && return nothing
    _libraryless_gram_surrogate_enabled() || return nothing
    planet_key = lowercase(strip(planet_name))
    surrogate_file = try
        Base.invokelatest(SimulationModel.EnvironmentModels._gram_default_surrogate_file, planet_key)
    catch
        ""
    end
    isempty(surrogate_file) && return nothing
    isfile(surrogate_file) || return nothing
    base_model = _GRAMOfflineSurrogateFallbackBase(planet_name=planet_key)
    return GRAMAtmosphereModelSurrogate(base_model, surrogate_file, nothing)
end

function _make_required_gram_density_model(
    planet_name::String,
    initial_time::InitialTime,
    truth::AtmosphereTruthConfig
)
    try
        return Base.invokelatest(
            GRAMAtmosphereModel;
            planet_name=planet_name,
            initial_time=initial_time,
            seed=truth.gram_seed,
            gram_min_relative_step_size=truth.gram_min_relative_step_size,
            gram_perturbation_scales=truth.gram_perturbation_scales,
            mars_map_year=truth.mars_map_year,
            mars_mgcm_dust_levels=truth.mars_mgcm_dust_levels,
            mars_dust_storm=truth.mars_dust_storm,
            mars_f107=truth.mars_f107,
            mars_wind_scales=truth.mars_wind_scales,
            mars_mola_heights=truth.mars_mola_heights,
            mars_min_max=truth.mars_min_max
        )
    catch err
        msg = sprint(showerror, err)
        if _is_gram_library_missing_error(err)
            offline_model = _try_libraryless_gram_surrogate(planet_name, truth)
            if offline_model !== nothing
                @warn "GRAM shared library unavailable; using library-less GRAM offline surrogate fallback for telemetry." planet=planet_name surrogate_file=offline_model.surrogate_file
                return offline_model
            end
        end
        throw(ErrorException("GRAM atmosphere initialization failed for planet=$planet_name initial_time=$(initial_time): $msg"))
    end
end

@inline function _make_spacecraft(cfg::SpacecraftConfig, ic::AbstractInitialCondition)
    return make_three_body_spacecraft(
        bus_dims=cfg.bus_dims,
        panel_dims=cfg.panel_dims,
        bus_mass=cfg.bus_mass_kg,
        panel_mass_each=cfg.panel_mass_each_kg,
        panel_offset_y=cfg.panel_offset_y_m,
        ic=ic,
        prop_mass=cfg.prop_mass_kg,
        id=cfg.id,
        bus_ram_face=cfg.bus_ram_face
    )
end

function _with_environment_wind(args::SimulationConfiguration, include_wind::Bool)::SimulationConfiguration
    env = args.environment_model
    env_updated = EnvironmentModel(
        planet=env.planet,
        EI=env.EI,
        density_model=env.density_model,
        ephemerides_model=env.ephemerides_model,
        thermal_model=env.thermal_model,
        topography=env.topography,
        topo_degree=env.topo_degree,
        topo_order=env.topo_order,
        wind=include_wind
    )
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=args.simulation_settings,
        mission_configuration=args.mission_configuration,
        environment_model=env_updated,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

@inline function _has_campaign_maneuvers(cfg::OrbitEventsScenarioConfig)::Bool
    return !isempty(cfg.maneuver_orbit_numbers)
end

function _with_campaign_maneuvers(args::SimulationConfiguration, cfg::OrbitEventsScenarioConfig)::SimulationConfiguration
    !_has_campaign_maneuvers(cfg) && return args

    n_sats = length(args.dynamics_model.spacecraft)
    n_sats == 1 || throw(ArgumentError(
        "Telemetry campaign maneuvers currently support exactly one spacecraft; got $n_sats for scenario $(cfg.name)."
    ))
    cfg.maneuver_thrust_n > 0.0 || throw(ArgumentError("maneuvers.thrust_n must be > 0.0 for scenario $(cfg.name)"))
    cfg.maneuver_isp_s > 0.0 || throw(ArgumentError("maneuvers.isp_s must be > 0.0 for scenario $(cfg.name)"))
    cfg.maneuver_guidance_rate_s > 0.0 || throw(ArgumentError("maneuvers.guidance_rate_s must be > 0.0 for scenario $(cfg.name)"))
    cfg.maneuver_control_rate_s > 0.0 || throw(ArgumentError("maneuvers.control_rate_s must be > 0.0 for scenario $(cfg.name)"))

    thruster = BaseThrusterModel(
        thrust=fill(cfg.maneuver_thrust_n, n_sats),
        direction=fill(0.0, n_sats),
        Δv=fill(0.0, n_sats),
        start_burn_time=fill(-1.0, n_sats),
        stop_burn_time=fill(-1.0, n_sats),
        Isp=fill(cfg.maneuver_isp_s, n_sats)
    )
    # Diagnostic replay scaling: convert flight apoapsis altitudes to radii
    # with the equatorial radius; the flight/sim RATIO is insensitive to the
    # (identical) altitude-to-radius convention at the <0.2% level.
    flight_apo_radius_m = if cfg.maneuver_replay_scale_mode == "flight_apoapsis_ratio"
        planet = args.environment_model.planet
        println("maneuver_replay_scale context=$(cfg.name) mode=flight_apoapsis_ratio burns=$(length(cfg.maneuver_flight_apoapsis_alt_m))")
        Float64[alt + planet.Rp_e for alt in cfg.maneuver_flight_apoapsis_alt_m]
    else
        Float64[]
    end
    guidance_effector = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=cfg.maneuver_orbit_numbers,
        maneuver_Δv=cfg.maneuver_delta_v_mps,
        maneuver_flight_apoapsis_radius_m=flight_apo_radius_m
    )
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=args.simulation_settings,
        mission_configuration=args.mission_configuration,
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=GuidanceModel(
            guidance_effectors=(guidance_effector,),
            guidance_rates=[cfg.maneuver_guidance_rate_s]
        ),
        navigation_model=args.navigation_model,
        control_model=ControlModel(
            control_effectors=(thruster,),
            control_rates=[cfg.maneuver_control_rate_s]
        ),
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

function _with_orbit_mission(
    args::SimulationConfiguration,
    target_orbits::Int,
    mission_time_s::Float64
)::SimulationConfiguration
    mc = args.mission_configuration
    mission_cfg = MissionConfiguration(
        mission_type=MissionOrbits,
        keplerian=mc.keplerian,
        number_of_orbits=target_orbits,
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

# Body-mean-equator inertial axes expressed in J2000: z along the body pole,
# x along the ascending node of the body equator on the J2000 equator.
function _body_equator_frame_rotation(pole_j2000::SVector{3, Float64})::SMatrix{3, 3, Float64, 9}
    ẑ = pole_j2000 / norm(pole_j2000)
    node = SVector{3, Float64}(-ẑ[2], ẑ[1], 0.0)
    node_mag = norm(node)
    node_mag <= 1e-12 && return SMatrix{3, 3, Float64}(1.0I)
    x̂ = node / node_mag
    ŷ = cross(ẑ, x̂)
    return hcat(x̂, ŷ, ẑ)
end

# Flight-dynamics products for Venus/Mars scenarios publish osculating elements
# referenced to the body mean equator, while the engine propagates in J2000
# (dynamics_rhs build_initial_conditions). Elements flagged body_equator_inertial
# are therefore converted to a J2000 Cartesian state here, using the same pole
# model (planet_frame_lpi) the propagation itself uses.
function _initial_condition_in_j2000(
    ic::InitialCondition,
    planet,
    initial_time,
    element_frame::Symbol
)::SimulationModel.AbstractInitialCondition
    element_frame === :j2000 && return ic
    element_frame === :body_equator_inertial || throw(ArgumentError(
        "Unsupported element_frame=$element_frame for orbit-element initial conditions."
    ))
    model = SimulationModel.SpiceEphemeridesModel()
    et = SimulationModel.ephemerides_time_seconds(initial_time, model)
    l_pi = SimulationModel.planet_frame_lpi(planet, et, model)
    pole_j2000 = SVector{3, Float64}(l_pi[3, 1], l_pi[3, 2], l_pi[3, 3])
    rot = _body_equator_frame_rotation(pole_j2000)
    r_body, v_body = SimulationEngine.orbitalelemtorv(ic, planet)
    return CartesianInitialCondition(
        rot * SVector{3, Float64}(r_body),
        rot * SVector{3, Float64}(v_body)
    )
end

# Initial condition for an orbit-events scenario. The exact NAV-kernel Cartesian
# state (initial_state_j2000_m) takes precedence when the manifest provides it;
# the published osculating elements are then documentation only. Otherwise the
# elements are used, converted to J2000 axes if flagged body_equator_inertial.
function _scenario_initial_condition(
    cfg::OrbitEventsScenarioConfig,
    planet
)::SimulationModel.AbstractInitialCondition
    state = cfg.initial_state_j2000_m
    if state !== nothing
        println("initial_state_j2000: kernel Cartesian override active context=$(cfg.name)")
        return CartesianInitialCondition(
            SVector{3, Float64}(state[1], state[2], state[3]),
            SVector{3, Float64}(state[4], state[5], state[6])
        )
    end
    rp_m = planet.Rp_e + cfg.rp_altitude_m
    ic_elements = InitialCondition(
        ra=cfg.ra_m,
        rp=rp_m,
        i=cfg.i_deg,
        ω=cfg.aop_deg,
        Ω=cfg.raan_deg,
        ν=cfg.ta_deg
    )
    return _initial_condition_in_j2000(ic_elements, planet, cfg.initial_time, cfg.element_frame)
end

function _make_orbit_args(
    cfg::OrbitEventsScenarioConfig,
    target_orbits::Int;
    cd_scale::Float64=1.0,
    cr_override::Union{Nothing, Float64}=nothing
)::SimulationConfiguration
    planet = _planet_from_name(cfg.planet_name)
    rp_m = planet.Rp_e + cfg.rp_altitude_m
    ic = _scenario_initial_condition(cfg, planet)

    spacecraft = _make_spacecraft(cfg.spacecraft, ic)
    dynamic_effectors = _scenario_dynamic_effectors(
        cfg,
        planet,
        spacecraft;
        cd_scale=cd_scale,
        cr_override=cr_override
    )
    density_model = _scenario_density_model(cfg)

    period_s = _period_seconds(planet, cfg.ra_m, rp_m)
    mission_time_s = target_orbits * period_s

    base = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time_s,
        initial_time=cfg.initial_time,
        dynamic_effectors=dynamic_effectors,
        density_model=density_model,
        orientation_sim=false,
        keplerian=false,
        EI_km=cfg.EI_km,
        verbose=false
    )
    with_orbits = _with_orbit_mission(base, target_orbits, mission_time_s)
    with_wind = _with_environment_wind(with_orbits, cfg.include_wind)
    return _with_campaign_maneuvers(with_wind, cfg)
end

function _make_time_aligned_args(
    cfg::TimeAlignedScenarioConfig,
    mission_time_s::Float64,
    ic::AbstractInitialCondition;
    cd_scale::Float64=1.0,
    cr_override::Union{Nothing, Float64}=nothing
)::SimulationConfiguration
    planet = _planet_from_name(cfg.planet_name)
    spacecraft = _make_spacecraft(cfg.spacecraft, ic)
    dynamic_effectors = _scenario_dynamic_effectors(
        cfg,
        planet,
        spacecraft;
        cd_scale=cd_scale,
        cr_override=cr_override
    )
    density_model = _scenario_density_model(cfg)

    base = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time_s,
        initial_time=cfg.initial_time,
        dynamic_effectors=dynamic_effectors,
        density_model=density_model,
        orientation_sim=false,
        keplerian=false,
        EI_km=cfg.EI_km,
        verbose=false
    )
    return _with_environment_wind(base, cfg.include_wind)
end

@inline function _has_high_fidelity_effectors(args::SimulationConfiguration)::Bool
    return any(
        eff -> (
            eff isa GravitationalHarmonicsModel ||
            eff isa NBodyGravityModel ||
            eff isa SolarRadiationPressureModel
        ),
        args.dynamics_model.dynamic_effectors
    )
end

function _with_study_settings(args::SimulationConfiguration; quick::Bool=false)::SimulationConfiguration
    hf = _has_high_fidelity_effectors(args)
    rel_orbit = min(quick ? 5e-7 : 1e-7, STRICT_REL_ORBIT)
    abs_orbit = min(quick ? 5e-9 : 1e-9, STRICT_ABS_ORBIT)
    rel_atm = min(quick ? 1e-6 : 1e-7, STRICT_REL_ATM)
    abs_atm = min(quick ? 1e-8 : 1e-9, STRICT_ABS_ATM)
    dt_orbit_base = min(quick ? (hf ? 180.0 : 240.0) : (hf ? 60.0 : 120.0), STRICT_DT_ORBIT)
    dt_atm_base = min(quick ? (hf ? 2.0 : 5.0) : (hf ? 0.2 : 0.5), STRICT_DT_ATM)
    dt_orbit_env = _parse_positive_float_env("SPACEAGORA_TELEMETRY_DT_MAX_ORBIT")
    dt_atm_env = _parse_positive_float_env("SPACEAGORA_TELEMETRY_DT_MAX_ATM")
    dt_orbit = dt_orbit_env === nothing ? dt_orbit_base : min(dt_orbit_env, STRICT_DT_ORBIT)
    dt_atm = dt_atm_env === nothing ? dt_atm_base : min(dt_atm_env, STRICT_DT_ATM)
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            results_directory=args.simulation_settings.results_directory,
            generate_plots=false,
            generate_filenames=false,
            normalize=false,
            save_csv=true
        ),
        mission_configuration=MissionConfiguration(
            mission_type=args.mission_configuration.mission_type,
            keplerian=args.mission_configuration.keplerian,
            number_of_orbits=args.mission_configuration.number_of_orbits,
            mission_time=args.mission_configuration.mission_time,
            orientation_sim=args.mission_configuration.orientation_sim,
            num_steps_to_save=2000,
            data_rate=args.mission_configuration.data_rate
        ),
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=rel_orbit,
            abstol_orbit=abs_orbit,
            dt_max_orbit=dt_orbit,
            reltol_atmosphere=rel_atm,
            abstol_atmosphere=abs_atm,
            dt_max_atmosphere=dt_atm
        )
    )
end

function _save_fields_for_study()
    pos_getter = (u, t, integrator) -> begin
        out = Vector{SVector{3, Float64}}(undef, length(u.sc))
        @inbounds for i in eachindex(u.sc)
            out[i] = SVector{3, Float64}(u.sc[i].pos)
        end
        return out
    end
    vel_getter = (u, t, integrator) -> begin
        out = Vector{SVector{3, Float64}}(undef, length(u.sc))
        @inbounds for i in eachindex(u.sc)
            out[i] = SVector{3, Float64}(u.sc[i].vel)
        end
        return out
    end

    return [
        SaveField(:position, pos_getter; per_satellite=true, column_prefix="pos"),
        SaveField(:velocity, vel_getter; per_satellite=true, column_prefix="vel")
    ]
end
