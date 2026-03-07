# Canonical aggregator: no behavior ownership.
module SimulationModel

# --- Top-Level Dependencies ---
using LinearAlgebra
using StaticArrays
using CSV
using DataFrames
using Reexport

const SPICE_LOCK = ReentrantLock()
const GRAM_LOCK = ReentrantLock()

# --- Utils ---
include(joinpath(@__DIR__, "..", "core", "numerics", "quaternion_utils.jl"))
include(joinpath(@__DIR__, "..", "core", "state", "reference_system_config.jl"))
include(joinpath(@__DIR__, "..", "environment", "ephemerides", "planet_shapes.jl"))

# --- Submodules ---
# We include the files, which define their own modules.
# We then @reexport their public APIs.

# 1. Core abstract types
include(joinpath(@__DIR__, "..", "core", "types", "abstract_types.jl"))
@reexport using .AbstractTypes

	include(joinpath(@__DIR__, "..", "environment", "ephemerides", "planets.jl"))
	@reexport using .Planets
	include(joinpath(@__DIR__, "..", "environment", "ephemerides", "ephemerides_models.jl"))
	@reexport using .EphemeridesModels

# 2. Simple hardware data structs
include(joinpath(@__DIR__, "..", "vehicle", "spacecraft", "components.jl"))
@reexport using .Components

# 3. Main container structs (Link, Joint, Model)
include(joinpath(@__DIR__, "..", "vehicle", "spacecraft", "model.jl"))
@reexport using .PhysicalModel

# 4. Functions for building the model (add_...!)
include(joinpath(@__DIR__, "..", "vehicle", "spacecraft", "assembly.jl"))
@reexport using .Assembly

# 5. Functions for rotations and frames
include(joinpath(@__DIR__, "..", "vehicle", "kinematics", "kinematics.jl"))
@reexport using .Kinematics

# 6. Generic and actuator-specific hook surfaces
include(joinpath(@__DIR__, "..", "vehicle", "actuators", "actuator_hooks.jl"))
@reexport using .ActuatorHooks
include(joinpath(@__DIR__, "..", "vehicle", "actuators", "thruster", "thruster_hooks.jl"))
@reexport using .ThrusterHooks

# 7. Spacecraft structural ownership helpers (mass, inertia, geometry)
include(joinpath(@__DIR__, "..", "vehicle", "structure", "structure_models.jl"))
@reexport using .Structure

# --- Config types ---
include(joinpath(@__DIR__, "..", "core", "state", "simulation_configuration.jl"))
@reexport using .SimConfig

# --- Physical Models ---
include(joinpath(@__DIR__, "..", "environment", "physical_models.jl"))
@reexport using .EnvironmentModels

include(joinpath(@__DIR__, "..", "core", "types", "runtime_types.jl"))
@reexport using .ConfigTypes

# Shared parallel policy used by callbacks and dynamic effectors.
include(joinpath(@__DIR__, "..", "parallel", "policy", "parallel_policy.jl"))

# --- Rotational Dynamics ---
include(joinpath(@__DIR__, "..", "dynamics", "rotational", "rotational_models.jl"))
@reexport using .DynamicsRotational

# --- Translational Dynamics ---
include(joinpath(@__DIR__, "..", "dynamics", "translational", "translational_models.jl"))
@reexport using .DynamicsTranslational

# --- Coupled Dynamic Force/Torque Effectors ---
include(joinpath(@__DIR__, "..", "dynamics", "coupled", "force_torque_models.jl"))
@reexport using .DynamicEffectors

# --- IO Owners ---
include(joinpath(@__DIR__, "..", "io", "config", "io_config.jl"))
@reexport using .IOConfig
include(joinpath(@__DIR__, "..", "io", "serialization", "io_serialization.jl"))
@reexport using .IOSerialization
include(joinpath(@__DIR__, "..", "io", "logging", "io_logging.jl"))
@reexport using .IOLogging
include(joinpath(@__DIR__, "..", "io", "outputs", "io_outputs.jl"))
@reexport using .IOOutputs

# --- Mission Policy Owners ---
include(joinpath(@__DIR__, "..", "mission", "operations", "aerobraking_policy", "policy_types.jl"))
@reexport using .AerobrakingPolicy

# --- Guidance Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "guidance", "guidance_hooks.jl"))
@reexport using .GuidanceHooks
# --- Navigation Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "navigation", "navigation_hooks.jl"))
@reexport using .NavigationHooks
# --- Control Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "control", "control_hooks.jl"))
@reexport using .ControlHooks

	# --- Vehicle Thermal Models ---
	include(joinpath(@__DIR__, "..", "vehicle", "thermal", "thermal_models_module.jl"))
	@reexport using .VehicleThermalModels

	# --- No-GRAM onboarding presets ---
	include(joinpath(@__DIR__, "..", "core", "state", "no_gram_presets.jl"))
	@reexport using .NoGramPresets

# --- Integrator Callbacks ---
include(joinpath(@__DIR__, "..", "simulation", "callbacks", "callbacks.jl"))
@reexport using .SimulationCallbacks
end # module SimulationModel
