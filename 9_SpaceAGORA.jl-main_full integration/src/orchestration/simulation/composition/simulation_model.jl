# Canonical aggregator: no behavior ownership.
module SimulationModel

# --- Top-Level Dependencies ---
using LinearAlgebra
using StaticArrays
using CSV
using DataFrames
using Reexport

isdefined(parentmodule(@__MODULE__), :RuntimeServices) ||
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "..", "runtime_services.jl"))

isdefined(parentmodule(@__MODULE__), :SimulationLifecycle) ||
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "interfaces", "simulation_lifecycle.jl"))

isdefined(parentmodule(@__MODULE__), :CallbackServices) ||
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "..", "callback_services.jl"))

isdefined(parentmodule(@__MODULE__), :VisualizerServices) ||
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "..", "..", "..", "interfaces", "analysis", "visualization", "visualizer_services.jl"))

# --- Utils ---
include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "numerics", "quaternion_utils.jl"))
using .QuaternionMath
include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "numerics", "orientation_quaternion_utils.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "state", "reference_system_config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "models", "environment", "ephemerides", "planet_shapes.jl"))

# --- Submodules ---
# We include the files, which define their own modules.
# We then @reexport their public APIs.

# 1. Core abstract types
include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "types", "abstract_types.jl"))
@reexport using .AbstractTypes
include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "types", "effector_sampling.jl"))
@reexport using .EffectorSampling

include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "command_types.jl"))
@reexport using .CommandTypes

include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "planning", "adapters", "hypr", "support.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "planning", "adapters", "hypr", "hypr_utils.jl"))

include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "robotics", "robotics.jl"))
@reexport using .Robotics
include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "planning", "robot_arm", "planning.jl"))
@reexport using .RobotArmPlanning
include(joinpath(@__DIR__, "..", "..", "..", "models", "dynamics", "compliant_multibody", "compliant_multibody.jl"))
# Keep the defining module identity for existing serialized objects and HYPR.
const CompliantMultibody = ClothMultibody
@reexport using .CompliantMultibody
include(joinpath(@__DIR__, "..", "..", "..", "models", "dynamics", "compliant_multibody", "compliant_robot_arm_dynamics.jl"))
const CompliantRobotArmDynamics = ClothRobotArmDynamics
@reexport using .CompliantRobotArmDynamics

include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "actuators", "thruster", "thruster_models_module.jl"))
@reexport using .ThrusterModels
include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "guidance", "guidance_models.jl"))
@reexport using .GuidanceModels

	include(joinpath(@__DIR__, "..", "..", "..", "models", "environment", "ephemerides", "planets.jl"))
	@reexport using .Planets
	include(joinpath(@__DIR__, "..", "..", "..", "models", "environment", "ephemerides", "ephemerides_models.jl"))
	@reexport using .EphemeridesModels

include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "numerics", "geodesy.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "numerics", "orbital_geometry.jl"))

# 2. Simple hardware data structs
include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "spacecraft", "components.jl"))
@reexport using .Components

# 3. Main container structs (Link, Joint, Model)
include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "spacecraft", "model.jl"))
@reexport using .SpacecraftModels

include(joinpath(@__DIR__, "..", "..", "..", "models", "satellite_interlink", "1_constellation_graph_models.jl"))
@reexport using .InterLinkModels
include(joinpath(@__DIR__, "..", "..", "..", "models", "satellite_interlink", "2_dynamic_programming_models.jl"))
@reexport using .SchedulingPolicies

# 4. Functions for building the model (add_...!)
include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "spacecraft", "assembly.jl"))
@reexport using .Assembly

# Articulated-body dynamics library (needs only the spacecraft model types; the runtime
# types below hold its per-spacecraft run data).
include(joinpath(@__DIR__, "..", "..", "..", "models", "dynamics", "coupled", "articulated_body.jl"))
@reexport using .ArticulatedBody

# Compliant multibody attachments mounted on spacecraft links (the type itself lives in
# SpacecraftModels; this is its runtime and dynamics).
include(joinpath(@__DIR__, "..", "..", "..", "models", "dynamics", "compliant_multibody", "compliant_attachment_dynamics.jl"))
@reexport using .CompliantAttachmentDynamics

# Generic hooks for spacecraft an external integrator owns for a run (internal name space, not re-exported).
include(joinpath(@__DIR__, "..", "..", "..", "models", "dynamics", "coupled", "external_propagation.jl"))

# Private owner for shared coordinate transforms; preserve legacy qualified bindings.
module FrameTransforms
using ..Geodesy
using ..OrbitalGeometry
using ..EphemeridesModels: ephemerides_requires_spice, planet_frame_lpi
include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "interfaces", "reference_system.jl"))
export _EARTH_HIGH_PREC_BODY_FIXED_FRAME, _EARTH_FALLBACK_BODY_FIXED_FRAME
export _spice_lock, _spice_frame_lock, r_intor_p!, r_pintor_i
export _spice_body_fixed_frame, _body_fixed_state_xform
export _j2000_to_body_fixed_state, _body_fixed_to_j2000_state, _planet_flattening
export orbitalelemtorv, _wrap_2pi, _safe_acos, _rvtoorbitalelement_core, rvtoorbitalelement
export rtoalfadeltar, alfadeltartor, latlongtor, latlongtoOE, rtolatlong, rtolatlongrad
export latlongtoNED, orbital_elements_to_lvlh_quaternion, rotate_vector_by_quaternion
export rtn_dcm_from_inertial, _rtn_rate_rad_s, inertial_to_rtn_relative_state
export rtn_to_inertial_relative_state, rtn_accel_to_inertial
end

# 5. Functions for rotations and frames
include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "kinematics", "kinematics.jl"))
@reexport using .Kinematics

# 6. Generic and actuator-specific hook surfaces
include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "actuators", "thruster", "thruster_hooks.jl"))
@reexport using .ThrusterHooks

# 7. Spacecraft structural ownership helpers (mass, inertia, geometry)
include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "structure", "structure_models.jl"))
@reexport using .Structure

# --- Config types ---
include(joinpath(@__DIR__, "..", "config", "configuration.jl"))
@reexport using .SimConfig

# --- Physical Models ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "environment", "physical_models.jl"))
@reexport using .EnvironmentModels

include(joinpath(@__DIR__, "..", "..", "..", "models", "environment", "terrain", "terrain_models.jl"))
@reexport using .TerrainModels

include(joinpath(@__DIR__, "..", "..", "..", "foundation", "core", "types", "compat_model_codes.jl"))
@reexport using .LegacyModelCodes

include(joinpath(@__DIR__, "..", "state", "runtime_types.jl"))
@reexport using .ConfigTypes

# Shared parallel policy used by callbacks and dynamic effectors.
include(joinpath(@__DIR__, "..", "..", "parallel", "policy", "parallel_policy.jl"))

# --- Rotational Dynamics ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "dynamics", "rotational", "rotational_models.jl"))
@reexport using .DynamicsRotational

# --- Translational Dynamics ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "dynamics", "translational", "translational_models.jl"))
@reexport using .DynamicsTranslational

# --- Coupled Dynamic Force/Torque Effectors ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "dynamics", "coupled", "force_torque_models.jl"))
@reexport using .DynamicEffectors
include(joinpath(@__DIR__, "..", "..", "parallel", "routing", "effector_routing.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "models", "environment", "gravity", "gravity_effectors.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "models", "environment", "aerodynamics", "aerodynamic_effectors.jl"))

# Analytic cost model for routing decisions. Included after the effector owners
# because its work-count methods dispatch on their concrete types.
include(joinpath(@__DIR__, "..", "..", "parallel", "cost", "parallel_cost.jl"))

# --- IO Owners ---
include(joinpath(@__DIR__, "..", "..", "..", "interfaces", "io", "paths", "io_paths.jl"))
@reexport using .IOConfig
include(joinpath(@__DIR__, "..", "..", "..", "interfaces", "io", "checkpoints", "checkpoint_storage.jl"))
@reexport using .IOSerialization
include(joinpath(@__DIR__, "..", "..", "..", "interfaces", "io", "outputs", "io_outputs.jl"))
@reexport using .IOOutputs

# --- Mission Policy Owners ---
include(joinpath(@__DIR__, "..", "..", "mission", "operations", "aerobraking_policy", "policy_types.jl"))
@reexport using .AerobrakingPolicy

# --- Navigation Effectors ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "navigation", "navigation_hooks.jl"))
@reexport using .NavigationHooks
# --- Typed EDG platform services ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "guidance", "aerobraking", "typed_edg", "services.jl"))
# --- Guidance Effectors ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "guidance", "guidance_hooks.jl"))
@reexport using .GuidanceHooks
# --- Typed EDG prediction and decisions ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "guidance", "aerobraking", "typed_edg", "algorithms.jl"))
# --- Control Effectors ---
include(joinpath(@__DIR__, "..", "..", "..", "models", "gnc", "control", "control_hooks.jl"))
@reexport using .ControlHooks

	# --- Vehicle Thermal Models ---
	include(joinpath(@__DIR__, "..", "..", "..", "models", "vehicle", "thermal", "thermal_models_module.jl"))
	@reexport using .VehicleThermalModels

	# --- No-GRAM onboarding presets ---
	include(joinpath(@__DIR__, "..", "config", "presets", "no_gram_presets.jl"))
	@reexport using .NoGramPresets

# --- Integrator Callbacks ---
include(joinpath(@__DIR__, "..", "..", "..", "interfaces", "analysis", "visualization", "scene", "scene_visualization.jl"))
@reexport using .SceneVisualization
include(joinpath(@__DIR__, "..", "callbacks", "callbacks.jl"))
@reexport using .SimulationCallbacks
include(joinpath(@__DIR__, "..", "config", "example_configuration.jl"))
@reexport using .ExampleConfiguration
end # module SimulationModel
