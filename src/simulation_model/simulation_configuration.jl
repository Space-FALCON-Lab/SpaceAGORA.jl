module SimConfig
    export SimulationConfiguration, InitialTime, IntegrationTolerances, FilePaths, SimulationSettings, MissionConfiguration, EnvironmentModel
    using ..AbstractTypes: AbstractPlanet, AbstractDensityModel, AbstractThermalModel
    using ..PhysicalModel: DynamicsModel
    using ..Planets: Earth

    @kwdef struct InitialTime
        year::Int32 = 2000
        month::Int16 = 1
        day::Int16 = 1
        hour::Int16 = 0
        minute::Int16 = 0
        second::Float32 = 0.0
    end # struct InitialTime
    
    @kwdef struct IntegrationTolerances
        reltol::Float64 = 1e-9
        abstol::Float64 = 1e-11
        reltol_orbit::Float64 = 1e-9
        abstol_orbit::Float64 = 1e-11
        reltol_atmosphere::Float64 = 1e-9
        abstol_atmosphere::Float64 = 1e-11
        reltol_quaternion::Float64 = 1e-9
        abstol_quaternion::Float64 = 1e-11
        dt_max::Float64 = 1.0
        dt_max_orbit = 10.0
        dt_max_atmosphere = 1.0
    end # struct IntegrationTolerances

    @kwdef struct FilePaths
        results::String = "Results" # Directory to save results
        GRAM::String = "GRAM_Data" # Directory for GRAM atmospheric model data
        SPICE::String = "GRAM_Data/SPICE" # Directory for SPICE kernels
        topography_harmonics::String = "Topography_harmonics_data" # Directory for topography harmonics data (move to planet?)
        gravity_harmonics::String = "Gravity_harmonics_data" # Directory for gravity harmonics data (move to planet?)
    end # struct FilePaths

    @kwdef struct SimulationSettings
        # Misc simulation parameters
        results::Bool = true # Whether to save simulation results to a file
        verbose::Bool = false # Whether to print detailed simulation logs
        results_directory::String = "output" # Directory to save results
        generate_plots::Bool = true # Whether to generate plots after simulation
        generate_filenames::Bool = false # Whether to generate filenames with specifics of simulation parameters
        normalize::Bool = true # Whether to normalize the state vector for propagation
        save_csv::Bool = true # Whether to save results in CSV format in addition to feather
    end # struct SimulationSettings

    # TODO: Convert mission_type to an abstract type and have terminal conditions be functions that take in the state vector and return a boolean for whether the simulation should terminate. This will allow for more flexible mission configurations without needing to add new fields to the struct every time we want to add a new type of mission.
    @kwdef struct MissionConfiguration
        # Mission setup
        mission_type::String = "Time" # Indicator of the termination condition type (Time, number of orbits, etc.)
        keplerian::Bool = true # Whether to include step 2 (drag passage) as separate step or keep same integration parameters the whole time
        number_of_orbits::Int = 1 # Number of orbits to propagate for (if mission_type is "Orbits")
        mission_time::Float64 = 90.0*60.0*20.0 # Total mission time in seconds (if mission_type is "Time")
        orientation_sim::Bool = true # Whether to simulate orientation dynamics (if false, only position and velocity are simulated)
        num_steps_to_save::Int = 1000 # Number of time steps to store in memory during the simulation before writing to a file
    end # struct MissionConfiguration

    # TODO: Convert all the strings to abstract types to avoid needing if-else statements in complete passage and other functions. This will also make it easier to add new models in the future without needing to change the main code.
    @kwdef struct EnvironmentModel{P <: AbstractPlanet, D <: AbstractDensityModel, T <: AbstractThermalModel}
        # Physical environment model
        planet::P # Planet for which to run the simulation (used for gravity model, atmospheric model, etc.)
        density_model::D # Atmospheric model to use (Constant, Exponential, GRAM, NRLMSISE-00)
        topography::Bool = false # Whether to include topography in the simulation for altitude calculation
        topo_degree::Int = 90 # Maximum degree of spherical harmonics for topography
        topo_order::Int = 90 # Maximum order of spherical harmonics for topography
        wind::Bool = true # Whether to include wind in the simulation for atmospheric effects
        thermal_model::T # Thermal model to use (Maxwellian heat transfer, Convective and Radiative)
    end # struct EnvironmentModel

    @kwdef struct SimulationConfiguration{P <: AbstractPlanet, D <: AbstractDensityModel, T <: AbstractThermalModel, DM <: Tuple}
        file_paths::FilePaths = FilePaths() # File paths for data and results
        simulation_settings::SimulationSettings = SimulationSettings() # General simulation settings
        mission_configuration::MissionConfiguration = MissionConfiguration() # Mission-specific configuration
        environment_model::EnvironmentModel{P, D, T} # Physical environment models
        dynamics_model::DynamicsModel{DM} # Dynamics models to use for the simulation, e.g., drag, n-body gravity, gravity harmonics, etc. that calculate forces/torques on the spacecraft
        initial_time::InitialTime # Initial time for the simulation
        integration_tolerances::IntegrationTolerances = IntegrationTolerances() # Tolerances for the numerical integrator
    end # struct SimulationConfiguration
    
end # module SimConfig