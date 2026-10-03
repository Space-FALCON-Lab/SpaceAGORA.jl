# Owns one-run scenario configuration. Keep the existing SimConfig module identity.
module SimConfig
    export SimulationConfiguration, InitialTime, IntegrationTolerances, FilePaths, SimulationSettings, MissionConfiguration, EnvironmentModel, MissionType, MissionTime, MissionOrbits
    export SolverConfig
    using ..AbstractTypes: AbstractPlanet, AbstractDensityModel, AbstractThermalModel, AbstractEphemeridesModel
    using ..SpacecraftModels: DynamicsModel, GuidanceModel, ControlModel, NavigationModel
    using ..EphemeridesModels: SpiceEphemeridesModel
    using ..Planets: Earth

    include("run_settings.jl")
    include("solver_settings.jl")
    include("environment_settings.jl")
    include("simulation_configuration.jl")
end # module SimConfig
