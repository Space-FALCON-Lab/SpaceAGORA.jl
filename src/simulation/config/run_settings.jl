# Run epoch, duration, sampling, paths and output/checkpoint settings.
# Included inside SimConfig by configuration.jl.

## 1. Mission Type (Duration, Orbits)
# 1.1. Definition
@enum MissionType::UInt8 begin
    MissionTime = 0x01
    MissionOrbits = 0x02
end

# 1.2. Warn when mission_type is passed as a String or Symbol instead of the MissionType enum.
const _deprecated_mission_type_input_warned = Ref(false) # Tracks whether the deprecation warning for mission_type input has been issued.
@inline _warn_deprecated_config_enabled() = get(ENV, "SPACEAGORA_WARN_DEPRECATED_CONFIG", "1") == "1"
@inline function _warn_deprecated_mission_type_input!(mission_type)
    if !_warn_deprecated_config_enabled() || _deprecated_mission_type_input_warned[]
        return nothing
    end
    _deprecated_mission_type_input_warned[] = true
    @warn "Passing mission_type=$(repr(mission_type)) as String/Symbol is deprecated; pass MissionType (MissionTime/MissionOrbits) instead."
    return nothing
end

# 1.3. Convert different input types (String, Symbol) to the MissionType enum.
@inline function _parse_mission_type(mission_type::MissionType)::MissionType # Directly returns the input if it's already a MissionType.
    return mission_type
end

@inline function _parse_mission_type(mission_type::Symbol)::MissionType # Converts a Symbol to a MissionType by first converting it to a String.
    return _parse_mission_type(String(mission_type))
end

@inline function _parse_mission_type(mission_type::AbstractString)::MissionType # Parses a string to determine the corresponding MissionType.
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

# 1.4. Backward-compatible comparisons for downstream code still using string/symbol checks.
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

# 2.2. Initial Time
@kwdef struct InitialTime
    year::Int32 = 2000
    month::Int16 = 1
    day::Int16 = 1
    hour::Int16 = 0
    minute::Int16 = 0
    second::Float32 = 0.0
end # struct InitialTime

# 2.4. File Paths
@kwdef struct FilePaths
    results::String = "Results" # Directory to save results
    GRAM::String = "data/GRAMSuite.jl/GRAM Suite 2.0" # Directory for GRAM atmospheric model data
    SPICE::String = "data/GRAMSuite.jl/GRAM Suite 2.0/SPICE" # Directory for SPICE kernels
    topography_harmonics::String = "data/Topography_harmonics_data" # Directory for topography harmonics data (move to planet?)
    gravity_harmonics::String = "data/Gravity_harmonics_data" # Directory for gravity harmonics data (move to planet?)
end # struct FilePaths

# 2.5. Simulation Settings
@kwdef struct SimulationSettings
    # Misc simulation parameters
    results::Bool = true # Whether to save simulation results to a file
    verbose::Bool = false # Whether to print detailed simulation logs
    results_directory::String = "output" # Directory to save results
    generate_plots::Bool = true # Whether to generate plots after simulation
    generate_filenames::Bool = false # Whether to generate filenames with specifics of simulation parameters
    normalize::Bool = false # Legacy compatibility flag; typed run_simulation propagates SI-state directly
    save_csv::Bool = true # Whether to save results in CSV format in addition to feather
    save_visualization_scene::Bool = false # Write the viewer scene sidecar and the link_pose save field (off by default; see SceneVisualization)
    checkpoint_enabled::Bool = false # Periodically checkpoint state for restart safety
    checkpoint_interval_s::Float64 = 300.0 # Checkpoint cadence in seconds of simulated time
    checkpoint_directory::String = "" # Empty => use results_directory/checkpoints
    resume_from_checkpoint::Bool = false # Resume run from latest checkpoint if present
    articulated_live_pose_loads::Bool = false # Articulated spacecraft: aero and facet SRP per link at the live link poses, per-body gravity gradient (default off: loads act on the root body from the configured geometry)
end # struct SimulationSettings

# 2.6. Mission Configuration
# i) Struct
struct MissionConfiguration
    # Mission setup
    mission_type::MissionType # Indicator of the termination condition type (Time, number of orbits, etc.)
    keplerian::Bool # Whether to include step 2 (drag passage) as separate step or keep same integration parameters the whole time
    number_of_orbits::Int # Number of orbits to propagate for (if mission_type is "Orbits")
    mission_time::Float64 # Total mission time in seconds (if mission_type is "Time")
    orientation_sim::Bool # Whether to simulate orientation dynamics (if false, only position and velocity are simulated)
    num_steps_to_save::Int # Number of time steps to store in memory during the simulation before writing to a file
    data_rate::Float64 # Fixed data output rate in seconds, used for saveat in solve

    # ii) Inner constructor with data type conversion and validation
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
        data_rate > 0.0 || throw(ArgumentError("MissionConfiguration.data_rate must be > 0.0; got $data_rate."))
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

# ii) Outer constructor with default values and type parsing
function MissionConfiguration(;
    mission_type::Union{MissionType, AbstractString, Symbol}=MissionTime,
    keplerian::Bool=true,
    number_of_orbits::Integer=1,
    mission_time::Real=90.0*60.0*20.0*10.0,
    orientation_sim::Bool=false,
    num_steps_to_save::Integer=1000,
    data_rate::Float64=10.0
) # Constructor
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
