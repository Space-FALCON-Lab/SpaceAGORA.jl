module SimConfig
    export SimulationConfiguration, InitialTime, IntegrationTolerances, FilePaths, SimulationSettings, MissionConfiguration, EnvironmentModel, MissionType, MissionTime, MissionOrbits
    using ..AbstractTypes: AbstractPlanet, AbstractDensityModel, AbstractThermalModel
    using ..PhysicalModel: DynamicsModel, GuidanceModel, ControlModel, NavigationModel
    using ..Planets: Earth

    @enum MissionType::UInt8 begin
        MissionTime = 0x01
        MissionOrbits = 0x02
    end

    const _deprecated_mission_type_input_warned = Ref(false)
    @inline _warn_deprecated_config_enabled() = get(ENV, "SPACEAGORA_WARN_DEPRECATED_CONFIG", "1") == "1"
    @inline function _warn_deprecated_mission_type_input!(mission_type)
        if !_warn_deprecated_config_enabled() || _deprecated_mission_type_input_warned[]
            return nothing
        end
        _deprecated_mission_type_input_warned[] = true
        @warn "Passing mission_type=$(repr(mission_type)) as String/Symbol is deprecated; pass MissionType (MissionTime/MissionOrbits) instead."
        return nothing
    end

    @inline function _parse_mission_type(mission_type::MissionType)::MissionType
        return mission_type
    end

    @inline function _parse_mission_type(mission_type::Symbol)::MissionType
        return _parse_mission_type(String(mission_type))
    end

    @inline function _parse_mission_type(mission_type::AbstractString)::MissionType
        key = lowercase(strip(mission_type))
        if key == "time"
            _warn_deprecated_mission_type_input!(mission_type)
            return MissionTime
        elseif key == "orbits" || key == "orbit"
            _warn_deprecated_mission_type_input!(mission_type)
            return MissionOrbits
        end
        throw(ArgumentError("Invalid mission_type=$(repr(mission_type)). Valid mission types: \"Time\", \"Orbits\"."))
    end

    # Backward-compatible comparisons for downstream code still using string/symbol checks.
    @inline function Base.:(==)(lhs::MissionType, rhs::AbstractString)
        try
            return lhs == _parse_mission_type(rhs)
        catch
            return false
        end
    end
    @inline Base.:(==)(lhs::AbstractString, rhs::MissionType) = rhs == lhs
    @inline Base.:(==)(lhs::MissionType, rhs::Symbol) = lhs == String(rhs)
    @inline Base.:(==)(lhs::Symbol, rhs::MissionType) = rhs == lhs

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
        reltol_orbit::Float64 = 1e-6
        abstol_orbit::Float64 = 1e-8
        reltol_atmosphere::Float64 = 1e-7
        abstol_atmosphere::Float64 = 1e-9
        reltol_quaternion::Float64 = 1e-9
        abstol_quaternion::Float64 = 1e-11
        reltol_mass::Float64 = 1e-8
        abstol_mass::Float64 = 1e-10
        reltol_heat_load::Float64 = 1e-7
        abstol_heat_load::Float64 = 1e-9
        reltol_angular_rate::Float64 = 1e-8
        abstol_angular_rate::Float64 = 1e-10
        dt_max::Float64 = 1.0
        dt_max_orbit::Float64 = 30.0
        dt_max_atmosphere::Float64 = 1.0
    end # struct IntegrationTolerances

    @kwdef struct FilePaths
        results::String = "Results" # Directory to save results
        GRAM::String = "data/GRAMSuite.jl/GRAM Suite 2.0" # Directory for GRAM atmospheric model data
        SPICE::String = "data/GRAMSuite.jl/GRAM Suite 2.0/SPICE" # Directory for SPICE kernels
        topography_harmonics::String = "data/Topography_harmonics_data" # Directory for topography harmonics data (move to planet?)
        gravity_harmonics::String = "data/Gravity_harmonics_data" # Directory for gravity harmonics data (move to planet?)
    end # struct FilePaths

    @kwdef struct SimulationSettings
        # Misc simulation parameters
        results::Bool = true # Whether to save simulation results to a file
        verbose::Bool = false # Whether to print detailed simulation logs
        results_directory::String = "output" # Directory to save results
        generate_plots::Bool = true # Whether to generate plots after simulation
        generate_filenames::Bool = false # Whether to generate filenames with specifics of simulation parameters
        normalize::Bool = false # Legacy compatibility flag; typed run_simulation propagates SI-state directly
        save_csv::Bool = true # Whether to save results in CSV format in addition to feather
        checkpoint_enabled::Bool = false # Periodically checkpoint state for restart safety
        checkpoint_interval_s::Float64 = 300.0 # Checkpoint cadence in seconds of simulated time
        checkpoint_directory::String = "" # Empty => use results_directory/checkpoints
        resume_from_checkpoint::Bool = false # Resume run from latest checkpoint if present
    end # struct SimulationSettings

    @kwdef struct MissionConfiguration
        # Mission setup
        mission_type::MissionType = MissionTime # Indicator of the termination condition type (Time, number of orbits, etc.)
        keplerian::Bool = true # Whether to include step 2 (drag passage) as separate step or keep same integration parameters the whole time
        number_of_orbits::Int = 1 # Number of orbits to propagate for (if mission_type is "Orbits")
        mission_time::Float64 = 90.0*60.0*20.0*10.0 # Total mission time in seconds (if mission_type is "Time")
        orientation_sim::Bool = false # Whether to simulate orientation dynamics (if false, only position and velocity are simulated)
        num_steps_to_save::Int = 1000 # Number of time steps to store in memory during the simulation before writing to a file
        data_rate::Float64 = 10.0 # Fixed data output rate, used in saveat argument in solve

        function MissionConfiguration(
            mission_type::MissionType,
            keplerian::Bool,
            number_of_orbits::Integer,
            mission_time::Real,
            orientation_sim::Bool,
            num_steps_to_save::Integer,
            data_rate::Float64=10.0
        )
            number_of_orbits > 0 || throw(ArgumentError("MissionConfiguration.number_of_orbits must be > 0; got $number_of_orbits."))
            mission_time > 0 || throw(ArgumentError("MissionConfiguration.mission_time must be > 0; got $mission_time."))
            num_steps_to_save > 0 || throw(ArgumentError("MissionConfiguration.num_steps_to_save must be > 0; got $num_steps_to_save."))
            return new(
                mission_type,
                keplerian,
                Int(number_of_orbits),
                Float64(mission_time),
                orientation_sim,
                Int(num_steps_to_save),
                data_rate
            )
        end
    end # struct MissionConfiguration

    # Backward-compatible constructor for existing string/symbol call sites.
    function MissionConfiguration(;
        mission_type::Union{MissionType, AbstractString, Symbol}=MissionTime,
        keplerian::Bool=true,
        number_of_orbits::Integer=1,
        mission_time::Real=90.0*60.0*20.0*10.0,
        orientation_sim::Bool=false,
        num_steps_to_save::Integer=1000,
        data_rate::Float64=10.0
    )
        return MissionConfiguration(
            _parse_mission_type(mission_type),
            keplerian,
            number_of_orbits,
            mission_time,
            orientation_sim,
            num_steps_to_save,
            data_rate
        )
    end

    # TODO: Convert all the strings to abstract types to avoid needing if-else statements in complete passage and other functions. This will also make it easier to add new models in the future without needing to change the main code.
    @kwdef struct EnvironmentModel{P <: AbstractPlanet, D <: AbstractDensityModel, T <: AbstractThermalModel}
        # Physical environment model
        planet::P # Planet for which to run the simulation (used for gravity model, atmospheric model, etc.)
        EI::Float64 # Entry Interface altitude in km (used for determining when to start applying atmospheric effects)
        density_model::D # Atmospheric model to use (Constant, Exponential, GRAM, NRLMSISE-00)
        topography::Bool = false # Whether to include topography in the simulation for altitude calculation
        topo_degree::Int = 90 # Maximum degree of spherical harmonics for topography
        topo_order::Int = 90 # Maximum order of spherical harmonics for topography
        wind::Bool = true # Whether to include wind in the simulation for atmospheric effects
        thermal_model::T # Thermal model to use (Maxwellian heat transfer, Convective and Radiative)

        function EnvironmentModel(
            planet::P,
            EI::Real,
            density_model::D,
            topography::Bool,
            topo_degree::Integer,
            topo_order::Integer,
            wind::Bool,
            thermal_model::T
        ) where {P <: AbstractPlanet, D <: AbstractDensityModel, T <: AbstractThermalModel}
            EI >= 0 || throw(ArgumentError("EnvironmentModel.EI must be >= 0 km; got $EI."))
            topo_degree >= 0 || throw(ArgumentError("EnvironmentModel.topo_degree must be >= 0; got $topo_degree."))
            topo_order >= 0 || throw(ArgumentError("EnvironmentModel.topo_order must be >= 0; got $topo_order."))
            return new{P, D, T}(
                planet,
                Float64(EI),
                density_model,
                topography,
                Int(topo_degree),
                Int(topo_order),
                wind,
                thermal_model
            )
        end
    end # struct EnvironmentModel

    @kwdef struct SimulationConfiguration{P <: AbstractPlanet, D <: AbstractDensityModel, T <: AbstractThermalModel, DM <: Tuple}
        file_paths::FilePaths = FilePaths() # File paths for data and results
        simulation_settings::SimulationSettings = SimulationSettings() # General simulation settings
        mission_configuration::MissionConfiguration = MissionConfiguration() # Mission-specific configuration
        environment_model::EnvironmentModel{P, D, T} # Physical environment models
        dynamics_model::DynamicsModel{DM} # Dynamics models to use for the simulation, e.g., drag, n-body gravity, gravity harmonics, etc. that calculate forces/torques on the spacecraft
        guidance_model::GuidanceModel # Guidance models to use for the simulation, e.g., for calculating control inputs based on the state vector
        navigation_model::NavigationModel # Navigation models to use for the simulation, e.g., for calculating state estimates based on sensor data
        control_model::ControlModel # Control models to use for the simulation, e.g., for calculating control inputs based on the state vector
        initial_time::InitialTime # Initial time for the simulation
        integration_tolerances::IntegrationTolerances = IntegrationTolerances() # Tolerances for the numerical integrator
    end # struct SimulationConfiguration
    
end # module SimConfig
