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
    interlink_model::Union{Nothing, InterlinkModel} = nothing
    scheduling_policy_model::SchedulingPolicyModel = SchedulingPolicyModel(SemimajorAxisRatePolicy()) # or a scheduling policy alone
end # struct SimulationConfiguration

# Preserve the positional constructor used before interlinks were added.
SimulationConfiguration(file_paths, simulation_settings, mission_configuration,
    environment_model, dynamics_model, guidance_model, navigation_model, control_model,
    initial_time, integration_tolerances, solver_config, external_propagators) =
    SimulationConfiguration(file_paths, simulation_settings, mission_configuration,
        environment_model, dynamics_model, guidance_model, navigation_model, control_model,
        initial_time, integration_tolerances, solver_config, external_propagators,
        nothing, SchedulingPolicyModel(SemimajorAxisRatePolicy()))

# A scheduling policy alone is taken as the run's scheduler.
SimulationConfiguration(file_paths, simulation_settings, mission_configuration,
    environment_model, dynamics_model, guidance_model, navigation_model, control_model,
    initial_time, integration_tolerances, solver_config, external_propagators,
    interlink_model, scheduling_policy_model::AbstractSchedulingPolicy) =
    SimulationConfiguration(file_paths, simulation_settings, mission_configuration,
        environment_model, dynamics_model, guidance_model, navigation_model, control_model,
        initial_time, integration_tolerances, solver_config, external_propagators,
        interlink_model, SchedulingPolicyModel(scheduling_policy_model))

"""
    _with_fields(x; changes...)

A copy of `x` with the named fields changed. Every other field is passed on as
it is (the same object, not a copy) through the type's own constructor, so the
constructor's checks run and a field added to the type later is carried without
editing the call. A name that is not a field of `x` throws an `ArgumentError`.

This is how configuration structs are copied: name only the fields that change,
never pass another instance's fields one by one. The constructor takes the fields
by keyword, or in field order for a type whose keyword constructor computes fields
rather than taking them (`SpacecraftModel`). That field-order constructor checks
nothing: a changed `SpacecraftModel` field is not validated, so change only fields
whose values are already valid (the copy without GNC changes only the GNC
containers). It is a shallow update; run_simulation owns mutable-state isolation.
"""
function _with_fields(x; changes...)
    T = typeof(x)
    names = fieldnames(T)
    for name in keys(changes)
        name in names || throw(ArgumentError(
            "$(nameof(T)) has no field `$(name)` to change; its fields are $(join(names, ", "))."))
    end
    values = map(name -> haskey(changes, name) ? changes[name] : getfield(x, name), names)
    constructor = Base.typename(T).wrapper
    _rebuilt_in_field_order(T) && return constructor(values...)
    return constructor(; NamedTuple{names}(values)...)
end

_rebuilt_in_field_order(::Type) = false
_rebuilt_in_field_order(::Type{SpacecraftModel}) = true

"""
    _with_configuration(args::SimulationConfiguration; overrides...)

Rebuild a configuration with named field overrides, preserving all other
fields and their references (see `_with_fields`). This is a shallow update;
run_simulation owns mutable-state isolation. Infer model types again when
models are replaced.
"""
_with_configuration(args::SimulationConfiguration; overrides...) = _with_fields(args; overrides...)
