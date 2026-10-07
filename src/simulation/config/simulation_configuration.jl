# Final scenario container and shallow copy-with-overrides helper.
# Included inside SimConfig by configuration.jl.

# 2.8. Simulation Configuration
@kwdef struct SimulationConfiguration{P <: AbstractPlanet, D <: AbstractDensityModel, E <: AbstractEphemeridesModel, T <: AbstractThermalModel, DM <: Tuple}
    file_paths::FilePaths = FilePaths() # File paths for data and results
    simulation_settings::SimulationSettings = SimulationSettings() # General simulation settings
    mission_configuration::MissionConfiguration = MissionConfiguration() # Mission-specific configuration
    environment_model::EnvironmentModel{P, D, E, T} # Physical environment models
    dynamics_model::DynamicsModel{DM} # Dynamics models to use for the simulation, e.g., drag, n-body gravity, gravity harmonics, etc. that calculate forces/torques on the spacecraft
    guidance_model::GuidanceModel = GuidanceModel() # Guidance models to use for the simulation, e.g., for calculating control inputs based on the state vector
    navigation_model::NavigationModel = NavigationModel() # Navigation models to use for the simulation, e.g., for calculating state estimates based on sensor data
    control_model::ControlModel = ControlModel() # Control models to use for the simulation, e.g., for calculating control inputs based on the state vector
    initial_time::InitialTime # Initial time for the simulation
    integration_tolerances::IntegrationTolerances = IntegrationTolerances() # Tolerances for the numerical integrator
    solver_config::Union{Nothing, SolverConfig} = nothing # nothing = read from env at run time
    external_propagators::Tuple = () # SimulationModel.ExternalPropagation owners of externally propagated spacecraft (shadow entries); () = none
end # struct SimulationConfiguration

"""
    _with_configuration(args::SimulationConfiguration; overrides...)

Rebuild a configuration with named field overrides, preserving all other
fields and their references. This is a shallow update; run_simulation owns
mutable-state isolation. Infer model types again when models are replaced.
"""
function _with_configuration(args::SimulationConfiguration; overrides...)
    names = fieldnames(typeof(args))
    fields = NamedTuple{names}(map(name -> getfield(args, name), names))
    return SimulationConfiguration(; merge(fields, (; overrides...))...)
end
