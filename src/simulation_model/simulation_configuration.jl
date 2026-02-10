module SimConfig
    export SimulationConfiguration, InitialCondition, InitialTime, IntegrationTolerances, FilePaths, SimulationSettings, MissionConfiguration, PhysicalModels
    using ..AbstractTypes: Planet
    using ..PhysicalModel: DynamicsModel
    using ..Planets: Earth

    @kwdef struct InitialCondition
        a::Float64 = 0.0 # Semimajor axis (m)
        e::Float64 = 0.0 # Eccentricity (nd)
        i::Float64 = 0.0 # Inclination (deg)
        ω::Float64 = 0.0 # Argument of periapsis (deg)
        Ω::Float64 = 0.0 # RAAN (deg)
        ν::Float64 = 0.0 # True anomaly (deg)
    end # struct InitialCondition

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
        GRAM::String = "GRAM_data" # Directory for GRAM atmospheric model data
        SPICE::String = "GRAM_data/SPICE" # Directory for SPICE kernels
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

    @kwdef struct MissionConfiguration
        # Mission setup
        mission_type::String = "Time" # Indicator of the termination condition type (Time, number of orbits, etc.)
        keplerian::Bool = true # Whether to include step 2 (drag passage) as separate step or keep same integration parameters the whole time
        number_of_orbits::Int = 1 # Number of orbits to propagate for (if mission_type is "Orbits")
        mission_time::Float64 = 90*60.0 # Total mission time in seconds (if mission_type is "Time")
        orientation_sim::Bool = true # Whether to simulate orientation dynamics (if false, only position and velocity are simulated)
        num_steps_to_save::Int = 1000 # Number of time steps to store in memory during the simulation before writing to a file
    end # struct MissionConfiguration

    @kwdef struct PhysicalModels
        # Physical environment model
        planet::Planet = Earth() # Planet for which to run the simulation (used for gravity model, atmospheric model, etc.)
        gravity_model::String = "Inverse squared and J2 effect" # Gravity model to use in the simulation (Constant, Inverse squared, Inverse squared with J2 effect)
        density_model::String = "GRAM" # Atmospheric model to use (Constant, Exponential, GRAM, NRLMSISE-00)
        topography::Bool = false # Whether to include topography in the simulation for altitude calculation
        topo_degree::Int = 90 # Maximum degree of spherical harmonics for topography
        topo_order::Int = 90 # Maximum order of spherical harmonics for topography
        wind::Bool = false # Whether to include wind in the simulation for atmospheric effects
        aerodynamic_model::String = "Mach-dependent" # Aerodynamic model to use (Constant, Mach-dependent, No-Ballistic flight with axial coefficient)
        thermal_model::String = "Maxwellian heat transfer" # Thermal model to use (Maxwellian heat transfer, Convective and Radiative)
    end # struct PhysicalModels

    @kwdef struct SimulationConfiguration
        file_paths::FilePaths = FilePaths() # File paths for data and results
        simulation_settings::SimulationSettings = SimulationSettings() # General simulation settings
        mission_configuration::MissionConfiguration = MissionConfiguration() # Mission-specific configuration
        physical_models::PhysicalModels = PhysicalModels() # Physical environment models
        dynamics_model::DynamicsModel = nothing # Dynamics models to use for the simulation
        initial_conditions::InitialCondition = InitialCondition() # Initial orbital elements for the simulation
        initial_time::InitialTime = InitialTime() # Initial time for the simulation
        integration_tolerances::IntegrationTolerances = IntegrationTolerances() # Tolerances for the numerical integrator
    end # struct SimulationConfiguration
    
    # Convenience constructors
    function InitialCondition(;ra::Float64, rp::Float64, i::Float64, ω::Float64, Ω::Float64, ν::Float64=180.0)
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        return InitialCondition(a, e, i, ω, Ω, ν) # Set true anomaly to 180 degrees by default
    end
end # module SimConfig