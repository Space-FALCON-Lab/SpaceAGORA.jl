# Canonical aggregator: no behavior ownership.
module SimulationModel

# --- Top-Level Dependencies ---
using LinearAlgebra
using StaticArrays
using CSV
using DataFrames
using Reexport

isdefined(parentmodule(@__MODULE__), :RuntimeServices) ||
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "..", "simulation", "runtime_services.jl"))

# --- Utils ---
include(joinpath(@__DIR__, "..", "core", "numerics", "quaternion_utils.jl"))
include(joinpath(@__DIR__, "..", "core", "state", "reference_system_config.jl"))
include(joinpath(@__DIR__, "..", "environment", "ephemerides", "planet_shapes.jl"))
include(joinpath(@__DIR__, "..", "core", "utils", "orbital_elements.jl"))

# --- Submodules ---
# We include the files, which define their own modules.
# We then @reexport their public APIs.

# 1. Core abstract types
include(joinpath(@__DIR__, "..", "core", "types", "abstract_types.jl"))
@reexport using .AbstractTypes
include(joinpath(@__DIR__, "..", "core", "types", "effector_sampling.jl"))
@reexport using .EffectorSampling

include(joinpath(@__DIR__, "..", "gnc", "command_types.jl"))
@reexport using .CommandTypes

include(joinpath(@__DIR__, "..", "gnc", "hypr", "hypr_utils.jl"))

include(joinpath(@__DIR__, "..", "vehicle", "robotics", "robotics.jl"))
@reexport using .Robotics
include(joinpath(@__DIR__, "..", "gnc", "robotics", "robot_arm_planning.jl"))
@reexport using .RobotArmPlanning
include(joinpath(@__DIR__, "..", "dynamics", "multibody_cloth", "cloth_multibody.jl"))
@reexport using .ClothMultibody
include(joinpath(@__DIR__, "..", "dynamics", "multibody_cloth", "cloth_robot_arm_dynamics.jl"))
@reexport using .ClothRobotArmDynamics

include(joinpath(@__DIR__, "..", "vehicle", "actuators", "thruster", "thruster_models_module.jl"))
@reexport using .ThrusterModels
include(joinpath(@__DIR__, "..", "gnc", "guidance", "guidance_models.jl"))
@reexport using .GuidanceModels

	include(joinpath(@__DIR__, "..", "environment", "ephemerides", "planets.jl"))
	@reexport using .Planets
	include(joinpath(@__DIR__, "..", "environment", "ephemerides", "ephemerides_models.jl"))
	@reexport using .EphemeridesModels

# 2. Simple hardware data structs
include(joinpath(@__DIR__, "..", "vehicle", "spacecraft", "components.jl"))
@reexport using .Components

# 3. Main container structs (Link, Joint, Model)
include(joinpath(@__DIR__, "..", "vehicle", "spacecraft", "model.jl"))
@reexport using .SpacecraftModels

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

include(joinpath(@__DIR__, "..", "core", "types", "compat_model_codes.jl"))
@reexport using .LegacyModelCodes

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
include(joinpath(@__DIR__, "..", "environment", "gravity", "gravity_effectors.jl"))
include(joinpath(@__DIR__, "..", "environment", "aerodynamics", "aerodynamic_effectors.jl"))

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

# --- Navigation Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "navigation", "navigation_hooks.jl"))
@reexport using .NavigationHooks
# --- Guidance Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "guidance", "guidance_hooks.jl"))
@reexport using .GuidanceHooks
# --- Control Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "control", "control_hooks.jl"))
@reexport using .ControlHooks

	# --- Vehicle Thermal Models ---
	include(joinpath(@__DIR__, "..", "vehicle", "thermal", "thermal_models_module.jl"))
	@reexport using .VehicleThermalModels

	# --- No-GRAM onboarding presets ---
	include(joinpath(@__DIR__, "..", "core", "state", "no_gram_presets.jl"))
	@reexport using .NoGramPresets

	# --- ORACLE scenario option types ---
	include(joinpath(@__DIR__, "..", "core", "types", "oracle_types.jl"))

export OracleOptions, _validate_options, _with
export ORACLE_PAPER_TARGET_ALTITUDES_KM, ORACLE_PAPER_TARGET_INCLINATIONS_DEG
export ORACLE_PAPER_HELPER_COUNTS, ORACLE_PAPER_FIXED_HELPER_ALTITUDE_KM, ORACLE_PAPER_FIXED_HELPER_INCLINATION_DEG

# --- Integrator Callbacks ---
include(joinpath(@__DIR__, "..", "simulation", "callbacks", "callbacks.jl"))
@reexport using .SimulationCallbacks

# rtn_dcm_from_inertial and rvtoorbitalelement are defined in reference_system.jl
# (included inside SimulationCallbacks) but not exported there; expose them here.
import .SimulationCallbacks: rtn_dcm_from_inertial, rvtoorbitalelement
export rtn_dcm_from_inertial, rvtoorbitalelement
end # module SimulationModel
