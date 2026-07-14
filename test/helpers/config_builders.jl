function build_config(;
    spacecraft::SpacecraftModel,
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    keplerian::Bool=true,
    simulation_settings::SimulationSettings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::SimulationModel.InitialTime=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    ephemerides_model=SpiceEphemeridesModel(),
    planet=EARTH
)
    environment_model = EnvironmentModel(
        planet=planet,
        EI=EI_km,
        density_model=density_model,
        ephemerides_model=ephemerides_model,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    )

    return SimulationConfiguration(
        simulation_settings=simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel([spacecraft], dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end

function build_config_multi(;
    spacecraft::Vector{SpacecraftModel},
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    keplerian::Bool=true,
    simulation_settings::SimulationSettings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::SimulationModel.InitialTime=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    ephemerides_model=SpiceEphemeridesModel(),
    planet=EARTH
)
    environment_model = EnvironmentModel(
        planet=planet,
        EI=EI_km,
        density_model=density_model,
        ephemerides_model=ephemerides_model,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    )

    return SimulationConfiguration(
        simulation_settings=simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end
