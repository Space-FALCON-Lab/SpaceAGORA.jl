# Canonical aggregator: no behavior ownership.
module SimulationModel

## 1. Top-Level Dependencies
using LinearAlgebra
using StaticArrays
using CSV
using DataFrames
using Reexport

## 2. Core: Utilities and Abstract Types
include(joinpath(@__DIR__, "..", "1_core", "numerics", "quaternion_utils.jl"))
include(joinpath(@__DIR__, "..", "1_core", "state", "reference_system_config.jl"))
include(joinpath(@__DIR__, "..", "1_core", "types", "abstract_types.jl"))
@reexport using .AbstractTypes
include(joinpath(@__DIR__, "..", "1_core", "types", "effector_sampling.jl"))
@reexport using .EffectorSampling

## 3. Environment: Shapes and Ephemerides
include(joinpath(@__DIR__, "..", "1.1_environment", "ephemerides", "planet_shapes.jl"))
include(joinpath(@__DIR__, "..", "1.1_environment", "ephemerides", "planets.jl"))
@reexport using .Planets
include(joinpath(@__DIR__, "..", "1.1_environment", "ephemerides", "ephemerides_models.jl"))
@reexport using .EphemeridesModels

## 4. Vehicle: Robotics, Spacecraft, Kinematics, Actuators, Structure, and Thermal
include(joinpath(@__DIR__, "..", "1.2_vehicle", "robotics", "robotics.jl"))
@reexport using .Robotics
include(joinpath(@__DIR__, "..", "1.2_vehicle", "actuators", "thruster", "thruster_models_module.jl"))
@reexport using .ThrusterModels
include(joinpath(@__DIR__, "..", "1.2_vehicle", "spacecraft", "components.jl"))
@reexport using .Components
include(joinpath(@__DIR__, "..", "1.2_vehicle", "spacecraft", "model.jl"))
@reexport using .SpacecraftModels
include(joinpath(@__DIR__, "..", "1.2_vehicle", "spacecraft", "assembly.jl"))
@reexport using .Assembly
include(joinpath(@__DIR__, "..", "1.2_vehicle", "kinematics", "kinematics.jl"))
@reexport using .Kinematics
include(joinpath(@__DIR__, "..", "1.2_vehicle", "actuators", "actuator_hooks.jl"))
@reexport using .ActuatorHooks
include(joinpath(@__DIR__, "..", "1.2_vehicle", "actuators", "thruster", "thruster_hooks.jl"))
@reexport using .ThrusterHooks
include(joinpath(@__DIR__, "..", "1.2_vehicle", "structure", "structure_models.jl"))
@reexport using .Structure
include(joinpath(@__DIR__, "..", "1.2_vehicle", "thermal", "thermal_models_module.jl"))
@reexport using .VehicleThermalModels

## 5. GNC: Commands, Planning, and Guidance Models
include(joinpath(@__DIR__, "..", "1.5_gnc", "command_types.jl"))
@reexport using .CommandTypes
include(joinpath(@__DIR__, "..", "1.5_gnc", "hypr", "hypr_utils.jl"))
include(joinpath(@__DIR__, "..", "1.5_gnc", "robotics", "robot_arm_planning.jl"))
@reexport using .RobotArmPlanning
include(joinpath(@__DIR__, "..", "1.5_gnc", "guidance", "guidance_models.jl"))
@reexport using .GuidanceModels

## 6. Core: Simulation Configuration
include(joinpath(@__DIR__, "..", "1_core", "state", "simulation_configuration.jl"))
@reexport using .SimConfig

## 7. Environment: Physical Models
include(joinpath(@__DIR__, "..", "1.1_environment", "physical_models.jl"))
@reexport using .EnvironmentModels

## 8. Core: Runtime Types and No-GRAM Onboarding Presets
include(joinpath(@__DIR__, "..", "1_core", "types", "compat_model_codes.jl"))
@reexport using .LegacyModelCodes
include(joinpath(@__DIR__, "..", "1_core", "types", "runtime_types.jl"))
@reexport using .ConfigTypes
include(joinpath(@__DIR__, "..", "1_core", "state", "no_gram_presets.jl"))
@reexport using .NoGramPresets

## 9. Parallel: Shared Policy
include(joinpath(@__DIR__, "..", "4_parallel", "policy", "parallel_policy.jl"))

## 10. Simulation: Constellation Types
# Included here (before force_torque_models.jl) so LaserLinkEffectors can reference constellation_struct.
include(joinpath(@__DIR__, "..", "5_simulation", "constellation.jl"))
@reexport using .Constellations

## 11. Dynamics: Multibody, Rotational, Translational, and Coupled
include(joinpath(@__DIR__, "..", "1.3_dynamics", "multibody_cloth", "cloth_multibody.jl"))
@reexport using .ClothMultibody
include(joinpath(@__DIR__, "..", "1.3_dynamics", "multibody_cloth", "cloth_robot_arm_dynamics.jl"))
@reexport using .ClothRobotArmDynamics
include(joinpath(@__DIR__, "..", "1.3_dynamics", "rotational", "rotational_models.jl"))
@reexport using .DynamicsRotational
include(joinpath(@__DIR__, "..", "1.3_dynamics", "translational", "translational_models.jl"))
@reexport using .DynamicsTranslational
include(joinpath(@__DIR__, "..", "1.3_dynamics", "coupled", "force_torque_models.jl"))
@reexport using .DynamicEffectors

## 12. Environment: Effectors
include(joinpath(@__DIR__, "..", "1.1_environment", "gravity", "gravity_effectors.jl"))
include(joinpath(@__DIR__, "..", "1.1_environment", "aerodynamics", "aerodynamic_effectors.jl"))

## 13. IO: Configuration, Serialization, Logging, and Outputs
include(joinpath(@__DIR__, "..", "3_io", "config", "io_config.jl"))
@reexport using .IOConfig
include(joinpath(@__DIR__, "..", "3_io", "serialization", "io_serialization.jl"))
@reexport using .IOSerialization
include(joinpath(@__DIR__, "..", "3_io", "logging", "io_logging.jl"))
@reexport using .IOLogging
include(joinpath(@__DIR__, "..", "3_io", "outputs", "io_outputs.jl"))
@reexport using .IOOutputs

## 14. Mission: Operations Policy
include(joinpath(@__DIR__, "..", "1.6_mission", "operations", "aerobraking_policy", "policy_types.jl"))
@reexport using .AerobrakingPolicy

## 15. GNC: Navigation, Guidance, and Control Hooks
include(joinpath(@__DIR__, "..", "1.5_gnc", "navigation", "navigation_hooks.jl"))
@reexport using .NavigationHooks
include(joinpath(@__DIR__, "..", "1.5_gnc", "guidance", "guidance_hooks.jl"))
@reexport using .GuidanceHooks
include(joinpath(@__DIR__, "..", "1.5_gnc", "control", "control_hooks.jl"))
@reexport using .ControlHooks

## 16. Simulation: Integrator Callbacks
include(joinpath(@__DIR__, "..", "5_simulation", "callbacks", "callbacks.jl"))
@reexport using .SimulationCallbacks
end # module SimulationModel
