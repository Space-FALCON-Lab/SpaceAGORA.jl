# Core-owned constructors only: this fixture must not import the provider.
using SpaceAGORA.SimulationModel: AerodynamicCoefficientfM

function missing_provider_configuration(directory)
    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0), panel_dims=(0.02, 2.0, 1.0),
        bus_mass=500.0, panel_mass_each=10.0, panel_offset_y=2.0,
        ic=InitialCondition(ra=planet.Rp_e + 120e3, rp=planet.Rp_e + 120e3,
            i=0.0, ω=0.0, Ω=0.0, ν=0.0), prop_mass=0.0, id=1)
    policy = AerobrakingEnergyDepletionConfig(
        guidance_modes=(:max_energy_depletion,),
        max_energy_submodes=(:heat_rate, :structural_load),
        controlled_panel_links=(2, 3), heat_rate_limit_w_cm2=0.15,
        structural_load_limit_pa=0.5, min_alpha_rad=1e-4, max_alpha_rad=pi/2)
    state = AerobrakingEnergyDepletionState(num_sats=1)
    guidance = AerobrakingEnergyDepletionGuidanceModel(policy, state)
    control = AerobrakingEnergyDepletionControlModel(policy, state)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=true, verbose=false, generate_plots=false, results_directory=directory),
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=0.1, keplerian=false, orientation_sim=false,
            num_steps_to_save=25),
        environment_model=EnvironmentModel(planet=planet, EI=160.0,
            density_model=ExponentialAtmosphereModel(1e-8, 120e3, 20e3; temperature_k=200.0),
            ephemerides_model=SimpleEphemeridesModel(), wind=false, topography=false,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet)),
        dynamics_model=DynamicsModel([spacecraft],
            (InverseSquaredGravityModel(), AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel((guidance,), [1/3]),
        control_model=ControlModel((control,), [1/3]),
        initial_time=InitialTime(year=2024, month=1, day=1),
        solver_config=SolverConfig(solver_mode=:tsit5, parallel=false))
end
