const _SPACEAGORA_RL_LIVE_BACKEND_MODES = (:spaceagora_physics, :spaceagora_full_physics, :spaceagora_marsgram)
const _SPACEAGORA_RL_LOAD_LOCK = ReentrantLock()
const _SPACEAGORA_RL_SPICE_SETUP_LOCK = ReentrantLock()
const _SPACEAGORA_RL_MARS_CACHE = Dict{String,Any}()
const _SPACEAGORA_RL_SPACEAGORA_MODULE = Ref{Union{Nothing,Module}}(nothing)
const _SPACEAGORA_RL_GRAMSUITE_LOADED = Ref(false)
const _SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS_FLOOR = 5_000_000
const _SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS_PER_PASS = 250_000

_spaceagora_base_root() = dirname(package_root())
_spaceagora_spice_path() = joinpath(_spaceagora_base_root(), "data", "GRAMSuite.jl", "GRAM Suite 2.0", "SPICE")
_spaceagora_gram_root() = joinpath(_spaceagora_base_root(), "data", "GRAMSuite.jl", "GRAM Suite 2.0")
_spaceagora_mars_harmonics_file() = joinpath(_spaceagora_base_root(), "data", "Gravity_harmonics_data", "Mars50c.csv")
function _spaceagora_repo_path(path::AbstractString)
    isempty(path) && return ""
    return isabspath(path) ? normpath(path) : normpath(joinpath(_spaceagora_base_root(), path))
end

function _spaceagora_harmonics_order(degree::Int, order::Int)
    degree <= 0 && return 0
    order == 0 && return 0
    order < 0 && return degree
    return min(order, degree)
end

_is_spaceagora_live_backend(mode::Symbol) = mode in _SPACEAGORA_RL_LIVE_BACKEND_MODES

function _load_spaceagora!(; load_gramsuite::Bool=true)
    lock(_SPACEAGORA_RL_LOAD_LOCK) do
        base_root = _spaceagora_base_root()
        gram_project = joinpath(base_root, "data", "GRAMSuite.jl")
        paths = load_gramsuite ? (gram_project, base_root) : (base_root,)
        for path in paths
            if !(path in LOAD_PATH)
                pushfirst!(LOAD_PATH, path)
            end
        end
        if load_gramsuite && !_SPACEAGORA_RL_GRAMSUITE_LOADED[]
            Base.require(Base.PkgId(Base.UUID("b50455af-6a46-4eae-bf92-8039261dd674"), "GRAMSuite"))
            _SPACEAGORA_RL_GRAMSUITE_LOADED[] = true
        end
        if _SPACEAGORA_RL_SPACEAGORA_MODULE[] === nothing
            _SPACEAGORA_RL_SPACEAGORA_MODULE[] =
                Base.require(Base.PkgId(Base.UUID("afbfb69f-5c0b-4832-b760-43725dff8540"), "SpaceAGORA"))
        end
        return _SPACEAGORA_RL_SPACEAGORA_MODULE[]::Module
    end
end

function _initial_time_from_datetime(dt::DateTime; spaceagora=nothing)
    SM = getproperty(spaceagora === nothing ? _load_spaceagora!() : spaceagora, :SimulationModel)
    return Base.invokelatest(getproperty(SM, :InitialTime);
        year=Int32(Dates.year(dt)),
        month=Int16(Dates.month(dt)),
        day=Int16(Dates.day(dt)),
        hour=Int16(Dates.hour(dt)),
        minute=Int16(Dates.minute(dt)),
        second=Float32(Dates.second(dt) + Dates.millisecond(dt) / 1.0e3),
    )
end

function _spaceagora_initial_condition(spaceagora, config, state, action::AerobrakingAction, planet)
    SM = getproperty(spaceagora, :SimulationModel)
    hp = clamp(periapsis_after_action_m(config, state, action), 50e3, 180e3)
    return Base.invokelatest(
        getproperty(SM, :InitialCondition);
        ra=state.apoapsis_radius_m,
        rp=planet.Rp_e + hp,
        i=rad2deg(state.inclination_rad),
        ω=rad2deg(state.argument_of_periapsis_rad),
        Ω=rad2deg(state.raan_rad),
        ν=rad2deg(state.true_anomaly_rad),
    )
end

function _spaceagora_campaign_mission_time_s(config, state, initial_periapsis_altitude_m::Real,
                                             campaign_max_passes::Integer)
    period = orbital_period_s(config, state.apoapsis_radius_m, Float64(initial_periapsis_altitude_m))
    pass_count = max(1, Int(campaign_max_passes))
    return max(2period, 1.25 * (pass_count + 1) * period)
end

function _spaceagora_density_model(spaceagora, config, state)
    SM = getproperty(spaceagora, :SimulationModel)
    model = config.spaceagora_atmosphere_model
    if model == :tabulated_flight
        TV = getproperty(spaceagora, :TelemetryVerification)
        path = _spaceagora_repo_path(config.spaceagora_tabulated_flight_file)
        truth = Base.invokelatest(
            getproperty(TV, :AtmosphereTruthConfig);
            assumption_id="odyssey_tabulated_flight_rl",
            atmosphere_model="tabulated_flight",
            atmosphere_dataset="Odyssey accelerometer per-pass density profiles (PDS odya_0001, Tolson)",
            space_weather_model="measured (flight)",
            solar_flux_model="measured (flight)",
            tabulated_flight_file=path,
            tabulated_flight_sigma=config.spaceagora_tabulated_flight_sigma,
        )
        return Base.invokelatest(
            getproperty(TV, :_make_tabulated_flight_density_model),
            _initial_time_from_datetime(state.epoch; spaceagora=spaceagora),
            truth,
        )
    end
    model in (:gram, :marsgram) ||
        throw(ArgumentError("unsupported SpaceAGORA atmosphere model: $(config.spaceagora_atmosphere_model)"))
    rand_cfg = config.randomization_config
    perturbation_scale = rand_cfg.marsgram_perturbation_scale
    return Base.invokelatest(
        getproperty(SM, :GRAMAtmosphereModel);
        planet_name="mars",
        gram_root_directory=_spaceagora_gram_root(),
        spice_directory=_spaceagora_spice_path(),
        seed=state.gram_seed,
        initial_time=_initial_time_from_datetime(state.epoch; spaceagora=spaceagora),
        gram_perturbation_scales=perturbation_scale,
    )
end

_spaceagora_live_needs_gramsuite(config) = config.spaceagora_atmosphere_model in (:gram, :marsgram)

function _spaceagora_physics_solver_maxiters(campaign_max_passes::Integer)
    pass_cap = max(1, Int(campaign_max_passes))
    default_maxiters = max(
        _SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS_FLOOR,
        _SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS_PER_PASS * pass_cap,
    )
    raw = strip(get(ENV, "SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS",
                    get(ENV, "SPACEAGORA_SOLVER_MAXITERS", "")))
    isempty(raw) && return default_maxiters
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS must be a positive integer, got '$raw'."))
    end
    parsed > 0 || throw(ArgumentError("SPACEAGORA_RL_PHYSICS_SOLVER_MAXITERS must be a positive integer, got $parsed."))
    return parsed
end

function _spaceagora_mars(spaceagora, spice_path::AbstractString)
    key = normpath(String(spice_path))
    return lock(_SPACEAGORA_RL_SPICE_SETUP_LOCK) do
        get!(_SPACEAGORA_RL_MARS_CACHE, key) do
            SM = getproperty(spaceagora, :SimulationModel)
            Base.invokelatest(getproperty(SM, :Mars), "", key)
        end
    end
end

mutable struct SpaceAGORAPhysicsSimulationTemplate
    planet::Any
    spacecraft::Any
    sun_gravity::Any
    harmonics::Any
    srp::Any
    aero_base::Any
    solver_cfg::Any
end

function _spaceagora_physics_simulation_template(
    config,
    state,
    action::AerobrakingAction;
    campaign_max_passes::Integer=config.termination_config.max_passes,
)
    spaceagora = _load_spaceagora!(; load_gramsuite=_spaceagora_live_needs_gramsuite(config))
    SM = getproperty(spaceagora, :SimulationModel)
    TV = getproperty(spaceagora, :TelemetryVerification)
    planet = deepcopy(_spaceagora_mars(spaceagora, _spaceagora_spice_path()))
    ic = _spaceagora_initial_condition(spaceagora, config, state, action, planet)
    spacecraft = Base.invokelatest(
        getproperty(TV, :make_three_body_spacecraft);
        bus_dims=(2.2, 2.6, 1.7),
        panel_dims=(0.01, 5.5 / 1.35, 2.6),
        bus_mass=391.0,
        panel_mass_each=10.0,
        panel_offset_y=2.6 / 2.0 + 5.5 / 4.0,
        ic=ic,
        reflection_coefficient=0.9,
        prop_mass=50.0,
        id=100,
    )
    sun_gravity = Base.invokelatest(
        getproperty(SM, :NBodyGravityModel);
        body_names=("Sun",),
        primary_body_name="Mars",
        planet=planet,
    )
    harmonics = Base.invokelatest(
        getproperty(SM, :GravitationalHarmonicsModel),
        config.spaceagora_gravity_harmonics_degree,
        _spaceagora_harmonics_order(
            config.spaceagora_gravity_harmonics_degree,
            config.spaceagora_gravity_harmonics_order,
        ),
        _spaceagora_repo_path(config.spaceagora_gravity_harmonics_file),
        planet;
        coefficients_normalized=true,
        j2_source=:file_c20,
    )
    srp = Base.invokelatest(
        getproperty(SM, :SolarRadiationPressureModel),
        spacecraft.root.reflection_coefficient,
        spacecraft.root.ref_area,
    )
    aero_base = Base.invokelatest(getproperty(SM, :AerodynamicCoefficientfM))
    solver_cfg = Base.invokelatest(
        getproperty(SM, :SolverConfig);
        solver_mode=:split_imex,
        maxiters=_spaceagora_physics_solver_maxiters(campaign_max_passes),
        split_imex_solver=:kencarp4,
    )
    return SpaceAGORAPhysicsSimulationTemplate(
        planet,
        spacecraft,
        sun_gravity,
        harmonics,
        srp,
        aero_base,
        solver_cfg,
    )
end

function _spaceagora_physics_simulation_configuration(config,
                                                      state,
                                                      action::AerobrakingAction;
                                                      prediction::Bool=false,
                                                      campaign_max_passes::Integer=config.termination_config.max_passes,
                                                      simulation_template::Union{Nothing,SpaceAGORAPhysicsSimulationTemplate}=nothing)
    spaceagora = _load_spaceagora!(; load_gramsuite=_spaceagora_live_needs_gramsuite(config))
    SM = getproperty(spaceagora, :SimulationModel)
    TV = getproperty(spaceagora, :TelemetryVerification)

    template = simulation_template === nothing ?
               _spaceagora_physics_simulation_template(
                   config,
                   state,
                   action;
                   campaign_max_passes=campaign_max_passes,
               ) :
               simulation_template
    planet = template.planet
    initial_time = _initial_time_from_datetime(state.epoch; spaceagora=spaceagora)
    periapsis_after_action = clamp(periapsis_after_action_m(config, state, action), 50e3, 180e3)
    ic = _spaceagora_initial_condition(spaceagora, config, state, action, planet)
    spacecraft = template.spacecraft
    spacecraft.initial_condition = ic
    sun_gravity = template.sun_gravity
    harmonics = template.harmonics
    srp = template.srp
    aero_base = template.aero_base
    aero = (isapprox(state.aerodynamic_cd_scale, 1.0; rtol=0.0, atol=1e-12) &&
            isapprox(state.aerodynamic_cl_scale, 1.0; rtol=0.0, atol=1e-12)) ?
           aero_base :
           Base.invokelatest(
               getproperty(TV, :ScaledAerodynamicCoefficientfM),
               aero_base,
               state.aerodynamic_cd_scale,
               state.aerodynamic_cl_scale,
           )
    density_model = _spaceagora_density_model(spaceagora, config, state)
    mission_time = _spaceagora_campaign_mission_time_s(
        config,
        state,
        periapsis_after_action,
        campaign_max_passes,
    )

    solver_cfg = template.solver_cfg
    base_args = Base.invokelatest(
        getproperty(TV, :make_example_config);
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time,
        initial_time=initial_time,
        dynamic_effectors=(sun_gravity, harmonics, srp, aero),
        density_model=density_model,
        orientation_sim=false,
        keplerian=false,
        EI_km=160.0,
        verbose=false,
        results=false,
        results_directory=joinpath(package_root(), "outputs", "spaceagora_physics_tmp"),
        solver_config=solver_cfg,
    )

    args = Base.invokelatest(
        getproperty(SM, :SimulationConfiguration);
        file_paths=base_args.file_paths,
        simulation_settings=base_args.simulation_settings,
        mission_configuration=Base.invokelatest(
            getproperty(SM, :MissionConfiguration);
            mission_type=base_args.mission_configuration.mission_type,
            keplerian=base_args.mission_configuration.keplerian,
            number_of_orbits=base_args.mission_configuration.number_of_orbits,
            mission_time=base_args.mission_configuration.mission_time,
            orientation_sim=base_args.mission_configuration.orientation_sim,
            num_steps_to_save=base_args.mission_configuration.num_steps_to_save,
            data_rate=10.0,
        ),
        environment_model=base_args.environment_model,
        dynamics_model=base_args.dynamics_model,
        guidance_model=base_args.guidance_model,
        navigation_model=base_args.navigation_model,
        control_model=base_args.control_model,
        initial_time=base_args.initial_time,
        integration_tolerances=Base.invokelatest(
            getproperty(SM, :IntegrationTolerances);
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=30.0,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_atmosphere=5.0,
        ),
        solver_config=solver_cfg,
    )

    metadata = (
        backend_mode=config.backend_mode,
        prediction=prediction,
        first_action_delta_v_mps=action.delta_v_mps,
        initial_periapsis_altitude_m=periapsis_after_action,
        campaign_max_passes=Int(campaign_max_passes),
        mission_time_s=mission_time,
        gravity_harmonics_degree=config.spaceagora_gravity_harmonics_degree,
        gravity_harmonics_order=_spaceagora_harmonics_order(
            config.spaceagora_gravity_harmonics_degree,
            config.spaceagora_gravity_harmonics_order,
        ),
        gravity_harmonics_file=_spaceagora_repo_path(config.spaceagora_gravity_harmonics_file),
        gravity_harmonics_coefficients_normalized=true,
        gravity_harmonics_j2_source=:file_c20,
        density_model=config.spaceagora_atmosphere_model == :tabulated_flight ?
                      :TabulatedFlightAtmosphereModel :
                      :GRAMAtmosphereModel,
        tabulated_flight_file=config.spaceagora_atmosphere_model == :tabulated_flight ?
                              _spaceagora_repo_path(config.spaceagora_tabulated_flight_file) :
                              "",
        tabulated_flight_sigma=config.spaceagora_tabulated_flight_sigma,
        aerodynamic_cd_scale=state.aerodynamic_cd_scale,
        aerodynamic_cl_scale=state.aerodynamic_cl_scale,
        initial_true_anomaly_deg=rad2deg(state.true_anomaly_rad),
        solver_mode=:split_imex,
        split_imex_solver=:kencarp4,
    )
    return args, metadata
end

function build_simulation_configuration(config, state, action::AerobrakingAction)
    if _is_spaceagora_live_backend(config.backend_mode)
        _, metadata = _spaceagora_physics_simulation_configuration(
            config,
            state,
            action;
            prediction=true,
            campaign_max_passes=1,
        )
        return metadata
    end
    return (
        backend_mode = config.backend_mode,
        phase = config.phase,
        pass_index = state.pass_index + 1,
        apoapsis_radius_m = state.apoapsis_radius_m,
        periapsis_altitude_m = state.periapsis_altitude_m,
        action_delta_v_mps = action.delta_v_mps,
        action_phi_deg = action.phi_deg,
    )
end
