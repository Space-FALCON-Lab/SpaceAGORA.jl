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
include(joinpath(@__DIR__, "..", "environment", "ephemerides", "planet_shapes.jl"))

# --- Submodules ---
# We include the files, which define their own modules.
# We then @reexport their public APIs.

# 1. Core abstract types
include(joinpath(@__DIR__, "..", "core", "types", "abstract_types.jl"))
@reexport using .AbstractTypes

include(joinpath(@__DIR__, "..", "environment", "ephemerides", "planets.jl"))
@reexport using .Planets

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

# 6. Functions for calculating effectors (thrusters, etc.)
include(joinpath(@__DIR__, "..", "vehicle", "actuators", "thruster_effectors.jl"))
@reexport using .Effectors

# 7. Spacecraft structural analysis helpers (get_COM, get_inertia, etc.)
include(joinpath(@__DIR__, "..", "vehicle", "spacecraft", "spacecraft_analysis.jl"))
@reexport using .Analysis

# --- Config types ---
include(joinpath(@__DIR__, "..", "core", "state", "simulation_configuration.jl"))
@reexport using .SimConfig

include(joinpath(@__DIR__, "..", "core", "types", "runtime_types.jl"))
@reexport using .ConfigTypes

# Shared parallel policy used by callbacks and dynamic effectors.
include(joinpath(@__DIR__, "..", "parallel", "policy", "parallel_policy.jl"))

# --- Rotational Dynamics ---
include(joinpath(@__DIR__, "..", "dynamics", "rotational", "rotational_models.jl"))
@reexport using .DynamicsRotational

# --- Dynamic Effectors ---
include(joinpath(@__DIR__, "..", "dynamics", "models", "dynamic_effectors.jl"))
@reexport using .DynamicEffectors

# --- Guidance Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "guidance", "effectors.jl"))
@reexport using .GuidanceEffectors
# --- Navigation Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "navigation", "effectors.jl"))
@reexport using .NavigationEffectors
# --- Control Effectors ---
include(joinpath(@__DIR__, "..", "gnc", "control", "effectors.jl"))
@reexport using .ControlEffectors

# --- Physical Models ---
include(joinpath(@__DIR__, "..", "environment", "physical_models.jl"))
@reexport using .EnvironmentModels

# --- Vehicle Thermal Models ---
include(joinpath(@__DIR__, "..", "vehicle", "thermal", "thermal_models_module.jl"))
@reexport using .VehicleThermalModels

# --- Integrator Callbacks ---
include(joinpath(@__DIR__, "..", "simulation", "callbacks", "callbacks.jl"))
@reexport using .SimulationCallbacks
end # module SimulationModel
