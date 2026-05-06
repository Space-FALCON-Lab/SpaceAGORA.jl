using StaticArrays

# Generic spacecraft identical across all bodies so perturbation comparisons are
# not confounded by vehicle differences.
const _SC_BUS_DIMS      = (2.2, 2.6, 1.7)
const _SC_PANEL_DIMS    = (0.01, 4.0, 2.0)
const _SC_BUS_MASS      = 400.0   # kg dry
const _SC_PANEL_MASS    = 10.0    # kg each
const _SC_PANEL_OFFSET  = _SC_BUS_DIMS[2] / 2.0 + _SC_PANEL_DIMS[2] / 4.0
const _SC_PROP_MASS     = 50.0    # kg propellant
const _SC_CR            = 0.9

_make_study_spacecraft(ic) = make_three_body_spacecraft(
    bus_dims=_SC_BUS_DIMS,
    panel_dims=_SC_PANEL_DIMS,
    bus_mass=_SC_BUS_MASS,
    panel_mass_each=_SC_PANEL_MASS,
    panel_offset_y=_SC_PANEL_OFFSET,
    ic=ic,
    reflection_coefficient=_SC_CR,
    prop_mass=_SC_PROP_MASS,
)

function _make_planet(body_name::Symbol)
    return if body_name == :earth
        Earth("", SPICE_PATH)
    elseif body_name == :mars
        Mars("", SPICE_PATH)
    elseif body_name == :venus
        Venus("", SPICE_PATH)
    elseif body_name == :titan
        Titan("", SPICE_PATH)
    else
        error("Unsupported body: $body_name")
    end
end

function _make_dynamic_effectors(pert_level::Symbol, cfg::BodyStudyConfig, planet, harmonics_model, spacecraft)
    base_drag = AerodynamicCoefficientfM()
    gravity   = InverseSquaredGravityModel()
    srp = SolarRadiationPressureModel(spacecraft.root.reflection_coefficient, spacecraft.root.ref_area)
    nbody = NBodyGravityModel(
        body_names=cfg.nbody_names,
        primary_body_name=string(cfg.planet_name),
        planet=planet,
    )

    return if pert_level == :point_mass
        (gravity, base_drag)
    elseif pert_level == :harmonics
        (gravity, harmonics_model, base_drag)
    elseif pert_level == :srp
        (gravity, srp, base_drag)
    elseif pert_level == :nbody
        (gravity, nbody, base_drag)
    else  # :full
        (gravity, harmonics_model, srp, nbody, base_drag)
    end
end

function make_perturbation_study_config(
    cfg::BodyStudyConfig,
    ic_params::NamedTuple,
    pert_level::Symbol,
    planet,
    harmonics_model,
    gram_model,
    run_id::Int;
    results_directory::String=OUTPUT_DIR,
    smoke_mode::Bool=false
)
    Rp_e = planet.Rp_e
    r_peri = Rp_e + ic_params.h_peri_m
    r_apo  = Rp_e + ic_params.h_apo_m
    ic = InitialCondition(
        ra=r_apo,
        rp=r_peri,
        i=ic_params.inc_deg,
        ω=ic_params.omega_deg,
        Ω=ic_params.RAAN_deg,
    )

    spacecraft = _make_study_spacecraft(ic)
    effectors  = _make_dynamic_effectors(pert_level, cfg, planet, harmonics_model, spacecraft)

    density_model  = gram_model
    corridor_model = CorridorManeuverGuidanceModel(cfg.h_corridor_min_m, cfg.h_corridor_max_m)
    thruster = BaseThrusterModel(
        thrust=[4.0],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0],
    )

    run_dir = joinpath(results_directory, "run_$(run_id)")

    initial_time = InitialTime(year=2000, month=1, day=1, hour=12, minute=0, second=0.0)

    mission_time = smoke_mode ? min(cfg.mission_time_s, 7200.0) : cfg.mission_time_s

    base = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time,
        initial_time=initial_time,
        dynamic_effectors=effectors,
        density_model=density_model,
        orientation_sim=false,
        keplerian=false,
        EI_km=cfg.EI_km,
        verbose=false,
        results=false,
        results_directory=run_dir,
    )

    return SimulationConfiguration(
        file_paths=base.file_paths,
        simulation_settings=base.simulation_settings,
        mission_configuration=base.mission_configuration,
        environment_model=base.environment_model,
        dynamics_model=base.dynamics_model,
        guidance_model=GuidanceModel(
            guidance_effectors=(corridor_model,),
            guidance_rates=[30.0],
        ),
        navigation_model=base.navigation_model,
        control_model=ControlModel(
            control_effectors=(thruster,),
            control_rates=[10.0],
        ),
        initial_time=base.initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=30.0,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_atmosphere=1.0,
        ),
    )
end
